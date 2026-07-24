/*
 * PyMOL glTF 2.0 / GLB export
 *
 * Exports the current scene as a GLB (binary glTF 2.0) file.
 * Uses the ray-tracing primitive pipeline to generate triangle meshes,
 * then serializes them into the GLB container format.
 *
 * (c) 2026 Schrodinger, Inc.
 */

#include "os_predef.h"
#include "os_std.h"

#include "CGO.h"
#include "Feedback.h"
#include "GLTF.h"
#include "MemoryDebug.h"
#include "Ray.h"
#include "Scene.h"
#include "Setting.h"
#include "Sphere.h"
#include "Util.h"
#include "Version.h"

#include <nlohmann/json.hpp>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <limits>
#include <memory>
#include <vector>

using json = nlohmann::json;

/* GLB constants */
static constexpr std::uint32_t GLB_MAGIC = 0x46546C67; /* "glTF" */
static constexpr std::uint32_t GLB_VERSION = 2;
static constexpr std::uint32_t GLB_CHUNK_JSON = 0x4E4F534A; /* "JSON" */
static constexpr std::uint32_t GLB_CHUNK_BIN = 0x004E4942;  /* "BIN\0" */
static constexpr std::uint32_t GLB_HEADER_SIZE = 12;
static constexpr std::uint32_t GLB_CHUNK_HEADER_SIZE = 8;

static constexpr float TRANS_PRECISION = 0.01f;

/* A GLB container addresses everything with 32-bit offsets */
static constexpr std::uint64_t GLB_MAX_SIZE = 0xFFFFFFFFull;

/**
 * Convert a display-space (sRGB) color channel to linear space.
 * @param c channel value, expected in the 0-1 range
 * @return the linear-space channel value
 * @note PyMOL color components are display-space values handed straight to
 *   OpenGL, whereas glTF 2.0 defines COLOR_0 and baseColorFactor as linear.
 */
static float srgbToLinear(float c)
{
  if (c <= 0.0f)
    return 0.0f;
  if (c >= 1.0f)
    return 1.0f;
  if (c <= 0.04045f)
    return c / 12.92f;
  return std::pow((c + 0.055f) / 1.055f, 2.4f);
}

/**
 * Collects triangles that share the same material, grouped by transparency
 * level. Each MeshGroup becomes a separate glTF mesh with its own material.
 */
struct MeshGroup {
  float transparency = 0.0f;
  std::vector<float> positions; /* x,y,z per vertex */
  std::vector<float> normals;   /* x,y,z per vertex */
  std::vector<float> colors;    /* r,g,b per vertex */
  std::vector<std::uint32_t> indices;
  std::uint32_t vertex_count = 0;

  /* bounding box for POSITION accessor */
  float min_pos[3] = {std::numeric_limits<float>::max(),
      std::numeric_limits<float>::max(), std::numeric_limits<float>::max()};
  float max_pos[3] = {std::numeric_limits<float>::lowest(),
      std::numeric_limits<float>::lowest(),
      std::numeric_limits<float>::lowest()};

  /**
   * Append a vertex with position, normal, and color, updating the bounding
   * box.
   * @param pos xyz position (3 floats)
   * @param norm xyz normal (3 floats)
   * @param col rgb color (3 floats, 0-1 range, display space)
   * @note The color is converted from display (sRGB) to linear space, which
   *   is what glTF 2.0 expects for COLOR_0.
   */
  void addVertex(const float* pos, const float* norm, const float* col)
  {
    positions.push_back(pos[0]);
    positions.push_back(pos[1]);
    positions.push_back(pos[2]);
    normals.push_back(norm[0]);
    normals.push_back(norm[1]);
    normals.push_back(norm[2]);
    colors.push_back(srgbToLinear(col[0]));
    colors.push_back(srgbToLinear(col[1]));
    colors.push_back(srgbToLinear(col[2]));

    for (int i = 0; i < 3; i++) {
      if (pos[i] < min_pos[i])
        min_pos[i] = pos[i];
      if (pos[i] > max_pos[i])
        max_pos[i] = pos[i];
    }
    vertex_count++;
  }

  /**
   * Append a triangle defined by three vertex indices.
   * @param i0 first vertex index
   * @param i1 second vertex index
   * @param i2 third vertex index
   */
  void addTriangle(std::uint32_t i0, std::uint32_t i1, std::uint32_t i2)
  {
    indices.push_back(i0);
    indices.push_back(i1);
    indices.push_back(i2);
  }

  bool empty() const { return indices.empty(); }
};

/**
 * Find or create a MeshGroup for the given transparency level. Groups are
 * matched with TRANS_PRECISION tolerance so nearly-identical transparencies
 * share a single glTF material.
 * @param groups vector of existing mesh groups
 * @param trans transparency value (0 = opaque, 1 = fully transparent)
 * @return reference to the matching or newly created group
 */
static MeshGroup& findOrCreateGroup(std::vector<MeshGroup>& groups, float trans)
{
  for (auto& g : groups) {
    if (fabsf(g.transparency - trans) < TRANS_PRECISION)
      return g;
  }
  groups.emplace_back();
  groups.back().transparency = trans;
  return groups.back();
}

/**
 * Tessellate a sphere primitive into triangles using the pre-computed
 * SphereRec mesh at the current sphere_quality setting.
 * @param group mesh group to append geometry to
 * @param prim sphere primitive (uses v1 for center, r1 for radius, c1 for
 *   color)
 * @param G global context (for sphere_quality setting and SphereRec lookup)
 */
static void addSphere(MeshGroup& group, const CPrimitive* prim, PyMOLGlobals* G)
{
  int sq = SettingGetGlobal_i(G, cSetting_sphere_quality);
  SphereRec* sp = G->Sphere->Sphere[sq];

  /* Use the pre-computed triangle mesh (Tri array) for spheres.
   * sp->Tri contains triangle indices into sp->dot[] */
  std::uint32_t base = group.vertex_count;

  /* Add all sphere vertices */
  for (int i = 0; i < sp->nDot; i++) {
    float pos[3] = {prim->v1[0] + prim->r1 * sp->dot[i][0],
        prim->v1[1] + prim->r1 * sp->dot[i][1],
        prim->v1[2] + prim->r1 * sp->dot[i][2]};
    float norm[3] = {sp->dot[i][0], sp->dot[i][1], sp->dot[i][2]};
    group.addVertex(pos, norm, prim->c1);
  }

  /* Add triangles */
  for (int i = 0; i < sp->NTri; i++) {
    int i0 = sp->Tri[i * 3 + 0];
    int i1 = sp->Tri[i * 3 + 1];
    int i2 = sp->Tri[i * 3 + 2];
    group.addTriangle(base + i0, base + i1, base + i2);
  }
}

#define GLB_MAX_EDGE 50
#define CONE_MIN_RADIUS 10e-6f // Adopted from COLLADA.cpp

/**
 * Convert a CGO round nub cap (emitted as a triangle strip) into indexed
 * triangles and append them to a mesh group.
 * @param group mesh group to append geometry to
 * @param cgo CGO containing the cap geometry (CGO_BEGIN/NORMAL/VERTEX ops)
 * @param cap_color rgb color for all cap vertices (3 floats)
 * @note The CGO is expected to contain a single triangle strip. Winding
 *   order alternates per the GL_TRIANGLE_STRIP convention.
 */
static void addCylinderCGOCap(
    MeshGroup& group, CGO* cgo, const float* cap_color)
{
  if (!cgo)
    return;

  int strip_vert_count = 0;
  float cur_normal[3] = {0, 0, 0};

  for (auto it = cgo->begin(); !it.is_stop(); ++it) {
    auto pc = it.data();
    int op = it.op_code();
    switch (op) {
    case CGO_BEGIN:
      strip_vert_count = 0;
      break;
    case CGO_NORMAL:
      cur_normal[0] = pc[0];
      cur_normal[1] = pc[1];
      cur_normal[2] = pc[2];
      break;
    case CGO_VERTEX: {
      float pos[3] = {pc[0], pc[1], pc[2]};
      group.addVertex(pos, cur_normal, cap_color);
      strip_vert_count++;

      if (strip_vert_count >= 3) {
        std::uint32_t cur = group.vertex_count - 1;
        if (strip_vert_count % 2 == 1) {
          group.addTriangle(cur - 2, cur - 1, cur);
        } else {
          group.addTriangle(cur - 1, cur - 2, cur);
        }
      }
      break;
    }
    }
  }
}

/**
 * Tessellate a cylinder, sausage, or cone primitive into triangles. Generates
 * the shaft as a ring of quads, plus flat or round caps at each end depending
 * on the primitive's cap types and the stick_round_nub setting.
 * @param group mesh group to append geometry to
 * @param prim cylinder/sausage/cone primitive (uses v1/v2 for endpoints,
 *   r1/r2 for radii, c1/c2 for endpoint colors, cap1/cap2 for cap types)
 * @param G global context (for stick_quality, stick_overlap, stick_nub,
 *   stick_round_nub settings)
 * @note Sausage primitives force round caps on both ends. Cone primitives
 *   allow different radii (r1, r2) at each endpoint.
 */
static void addCylinderWithCaps(
    MeshGroup& group, const CPrimitive* prim, PyMOLGlobals* G)
{
  float d[3], t[3], p0[3], p1[3], p2[3];
  float vv1[3], vv2[3], vvv1[3], vvv2[3];
  float x[GLB_MAX_EDGE + 1], y[GLB_MAX_EDGE + 1];
  float overlap, overlap2, nub, nub2, r2;

  int nEdge = SettingGetGlobal_i(G, cSetting_stick_quality);
  bool stick_round_nub = SettingGetGlobal_i(G, cSetting_stick_round_nub);

  cCylCap captype[2] = {prim->cap1, prim->cap2};
  std::unique_ptr<CGO> cgocap[2];

  if (prim->type == cPrimSausage) {
    captype[0] = captype[1] = cCylCapRound;
  }

  overlap = prim->r1 * SettingGetGlobal_f(G, cSetting_stick_overlap);
  nub = prim->r1 * SettingGetGlobal_f(G, cSetting_stick_nub);
  if (prim->type == cPrimCone) {
    r2 = prim->r2;
    overlap2 = prim->r2 * SettingGetGlobal_f(G, cSetting_stick_overlap);
    nub2 = prim->r2 * SettingGetGlobal_f(G, cSetting_stick_nub);
  } else {
    r2 = prim->r1;
    overlap2 = overlap;
    nub2 = nub;
  }

  if (nEdge > GLB_MAX_EDGE)
    nEdge = GLB_MAX_EDGE;
  subdivide(nEdge, x, y);

  p0[0] = prim->v2[0] - prim->v1[0];
  p0[1] = prim->v2[1] - prim->v1[1];
  p0[2] = prim->v2[2] - prim->v1[2];
  normalize3f(p0);

  copy3f(prim->v1, vv1);
  copy3f(vv1, vvv1);
  copy3f(prim->v2, vv2);
  copy3f(vv2, vvv2);

  d[0] = vv2[0] - vv1[0];
  d[1] = vv2[1] - vv1[1];
  d[2] = vv2[2] - vv1[2];
  get_divergent3f(d, t);
  cross_product3f(d, t, p1);
  normalize3f(p1);
  cross_product3f(d, p1, p2);
  normalize3f(p2);

  if (captype[0] == cCylCapRound) {
    if (stick_round_nub) {
      cgocap[0].reset(CGONew(G));
      CGORoundNub(cgocap[0].get(), prim->v1, p0, p1, p2, -1, nEdge, prim->r1);
    } else {
      for (int i = 0; i < 3; i++) {
        vv1[i] -= p0[i] * overlap;
        vvv1[i] = vv1[i] - p0[i] * nub;
      }
    }
  }
  if (captype[1] == cCylCapRound) {
    if (stick_round_nub) {
      cgocap[1].reset(CGONew(G));
      CGORoundNub(cgocap[1].get(), prim->v2, p0, p1, p2, 1, nEdge, prim->r1);
    } else {
      for (int i = 0; i < 3; i++) {
        vv2[i] += p0[i] * overlap2;
        vvv2[i] = vv2[i] + p0[i] * nub2;
      }
    }
  }

  /* Shaft */
  std::uint32_t base = group.vertex_count;

  for (int c = nEdge; c >= 0; c--) {
    float v[3];
    v[0] = p1[0] * x[c] + p2[0] * y[c];
    v[1] = p1[1] * x[c] + p2[1] * y[c];
    v[2] = p1[2] * x[c] + p2[2] * y[c];

    float bot[3] = {vv1[0] + v[0] * prim->r1, vv1[1] + v[1] * prim->r1,
        vv1[2] + v[2] * prim->r1};
    float top[3] = {vv2[0] + v[0] * r2, vv2[1] + v[1] * r2, vv2[2] + v[2] * r2};
    float norm[3] = {v[0], v[1], v[2]};

    group.addVertex(bot, norm, prim->c1);
    group.addVertex(top, norm, prim->c2);
  }

  for (int c = 0; c < nEdge; c++) {
    std::uint32_t bl = base + c * 2;
    std::uint32_t tl = base + c * 2 + 1;
    std::uint32_t br = base + (c + 1) * 2;
    std::uint32_t tr = base + (c + 1) * 2 + 1;
    group.addTriangle(bl, br, tl);
    group.addTriangle(tl, br, tr);
  }

  /* Flat caps */
  if (prim->r1 > CONE_MIN_RADIUS && !cgocap[0] && captype[0] != cCylCapNone) {
    float neg_p0[3] = {-p0[0], -p0[1], -p0[2]};
    std::uint32_t center_idx = group.vertex_count;
    group.addVertex(vvv1, neg_p0, prim->c1);
    for (int i = 0; i < nEdge; i++) {
      group.addTriangle(center_idx, base + i * 2, base + (i + 1) * 2);
    }
  }

  if (r2 > CONE_MIN_RADIUS && !cgocap[1] && captype[1] != cCylCapNone) {
    std::uint32_t center_idx = group.vertex_count;
    group.addVertex(vvv2, p0, prim->c2);
    for (int i = 0; i < nEdge; i++) {
      group.addTriangle(center_idx, base + (i + 1) * 2 + 1, base + i * 2 + 1);
    }
  }

  /* CGO round caps */
  for (int capIdx = 0; capIdx < 2; capIdx++) {
    if (cgocap[capIdx]) {
      const float* cap_color = (capIdx == 0) ? prim->c1 : prim->c2;
      addCylinderCGOCap(group, cgocap[capIdx].get(), cap_color);
    }
  }
}

/**
 * Add a triangle primitive's three vertices and one triangle index.
 * @param group mesh group to append geometry to
 * @param prim triangle primitive (uses v1-v3, n1-n3, c1-c3)
 * @note Winding order is flipped when TriangleReverse() returns true.
 */
static void addTriangle(MeshGroup& group, const CPrimitive* prim)
{
  std::uint32_t base = group.vertex_count;

  group.addVertex(prim->v1, prim->n1, prim->c1);
  group.addVertex(prim->v2, prim->n2, prim->c2);
  group.addVertex(prim->v3, prim->n3, prim->c3);

  if (TriangleReverse(const_cast<CPrimitive*>(prim))) {
    group.addTriangle(base, base + 2, base + 1);
  } else {
    group.addTriangle(base, base + 1, base + 2);
  }
}

/**
 * Round a size up to the next multiple of 4.
 * @param size byte count to align
 * @return size rounded up to 4-byte boundary
 * @note Required by the GLB spec for chunk alignment.
 */
static std::uint64_t alignTo4(std::uint64_t size)
{
  return (size + 3) & ~UINT64_C(3);
}

/**
 * @return true if the host stores multi-byte integers least significant byte
 *   first, in which case 32-bit arrays can be copied verbatim into the GLB
 *   binary chunk.
 */
static bool isLittleEndianHost()
{
  const std::uint32_t one = 1;
  std::uint8_t first_byte;
  std::memcpy(&first_byte, &one, 1);
  return first_byte == 1;
}

/**
 * Write a uint32 in little-endian byte order.
 * @param[out] dest destination, must have room for 4 bytes
 * @param val 32-bit value to write
 * @return pointer just past the written bytes
 */
static char* writeU32(char* dest, std::uint32_t val)
{
  dest[0] = static_cast<char>(val & 0xFF);
  dest[1] = static_cast<char>((val >> 8) & 0xFF);
  dest[2] = static_cast<char>((val >> 16) & 0xFF);
  dest[3] = static_cast<char>((val >> 24) & 0xFF);
  return dest + 4;
}

/**
 * Copy an array of 32-bit elements (float or uint32) in little-endian byte
 * order.
 * @param[out] dest destination, must have room for 4 * data.size() bytes
 * @param data source data
 */
template <typename T>
static void writeU32Array(char* dest, const std::vector<T>& data)
{
  static_assert(sizeof(T) == 4, "expected 32-bit elements");

  if (isLittleEndianHost()) {
    std::memcpy(dest, data.data(), data.size() * sizeof(T));
    return;
  }

  for (const T& value : data) {
    std::uint32_t bits;
    std::memcpy(&bits, &value, sizeof(T));
    dest = writeU32(dest, bits);
  }
}

/**
 * Byte offset and length of a single attribute's data within the binary
 * buffer. Corresponds to one glTF bufferView.
 */
struct BufferViewInfo {
  std::uint32_t offset;
  std::uint32_t length;
};

/**
 * Buffer layout info for one MeshGroup (positions, normals, colors, indices).
 */
struct GroupBufferInfo {
  BufferViewInfo positions;
  BufferViewInfo normals;
  BufferViewInfo colors;
  BufferViewInfo indices;
};

/**
 * Compute where each mesh group's vertex attributes and indices live inside
 * the single contiguous binary buffer of the GLB BIN chunk.
 *
 * Layout per group: [positions][normals][colors][indices]
 *
 * @param groups mesh groups containing vertex/index data
 * @param[out] infos byte offset and length for each group's attributes
 * @return the total binary buffer size in bytes
 * @note All sections are naturally 4-byte aligned (float and uint32 data).
 *   The result may exceed what a GLB container can address; the caller is
 *   responsible for rejecting it.
 */
static std::uint64_t computeBufferLayout(
    const std::vector<MeshGroup>& groups, std::vector<GroupBufferInfo>& infos)
{
  std::uint64_t offset = 0;
  infos.resize(groups.size());

  auto place = [&offset](BufferViewInfo& view, std::uint64_t length) {
    view.offset = static_cast<std::uint32_t>(offset);
    view.length = static_cast<std::uint32_t>(length);
    offset += length;
  };

  for (size_t g = 0; g < groups.size(); g++) {
    const auto& mesh = groups[g];
    auto& info = infos[g];

    place(info.positions, mesh.positions.size() * sizeof(float));
    place(info.normals, mesh.normals.size() * sizeof(float));
    place(info.colors, mesh.colors.size() * sizeof(float));
    place(info.indices, mesh.indices.size() * sizeof(std::uint32_t));
  }

  return offset;
}

/**
 * Serialize all mesh groups' vertex attributes and indices into the GLB BIN
 * chunk payload, following the layout from computeBufferLayout().
 * @param groups mesh groups containing vertex/index data
 * @param infos buffer layout info from computeBufferLayout()
 * @param[out] dest destination, must have room for the full buffer
 */
static void writeBinaryBuffer(const std::vector<MeshGroup>& groups,
    const std::vector<GroupBufferInfo>& infos, char* dest)
{
  for (size_t g = 0; g < groups.size(); g++) {
    const auto& mesh = groups[g];
    const auto& info = infos[g];

    writeU32Array(dest + info.positions.offset, mesh.positions);
    writeU32Array(dest + info.normals.offset, mesh.normals);
    writeU32Array(dest + info.colors.offset, mesh.colors);
    writeU32Array(dest + info.indices.offset, mesh.indices);
  }
}

/**
 * Build the glTF 2.0 JSON descriptor for the scene. Creates one mesh, node,
 * and PBR material per MeshGroup, with POSITION/NORMAL/COLOR_0 vertex
 * attributes and UNSIGNED_INT indices.
 * @param groups mesh groups (one per transparency level)
 * @param infos buffer layout info from buildBinaryBuffer()
 * @param buffer_size total size of the binary buffer in bytes
 * @return serialized JSON string
 */
static std::string buildGLTFJson(const std::vector<MeshGroup>& groups,
    const std::vector<GroupBufferInfo>& infos, std::uint32_t buffer_size)
{
  json root;

  /* Asset */
  root["asset"] = {
      {"version", "2.0"},
      {"generator", "PyMOL " _PyMOL_VERSION},
  };

  /* Scene */
  root["scene"] = 0;

  /* Build arrays */
  json nodes = json::array();
  json meshes = json::array();
  json materials = json::array();
  json accessors = json::array();
  json buffer_views = json::array();

  int accessor_idx = 0;
  int buffer_view_idx = 0;

  for (size_t g = 0; g < groups.size(); g++) {
    const auto& mesh = groups[g];
    const auto& info = infos[g];

    if (mesh.empty())
      continue;

    int mat_idx = static_cast<int>(materials.size());

    /* Material */
    json mat;
    json pbr;
    pbr["baseColorFactor"] = {1.0, 1.0, 1.0, 1.0 - mesh.transparency};
    pbr["metallicFactor"] = 0.0;
    pbr["roughnessFactor"] = 0.7;
    mat["pbrMetallicRoughness"] = pbr;

    if (mesh.transparency > TRANS_PRECISION) {
      mat["alphaMode"] = "BLEND";
    } else {
      mat["alphaMode"] = "OPAQUE";
    }
    /* PyMOL draws this geometry without back-face culling */
    mat["doubleSided"] = true;
    materials.push_back(mat);

    /* Buffer views: positions */
    int bv_pos = buffer_view_idx++;
    buffer_views.push_back({
        {"buffer", 0},
        {"byteOffset", info.positions.offset},
        {"byteLength", info.positions.length},
        {"target", 34962}, /* ARRAY_BUFFER */
    });

    /* Buffer views: normals */
    int bv_norm = buffer_view_idx++;
    buffer_views.push_back({
        {"buffer", 0},
        {"byteOffset", info.normals.offset},
        {"byteLength", info.normals.length},
        {"target", 34962},
    });

    /* Buffer views: colors */
    int bv_col = buffer_view_idx++;
    buffer_views.push_back({
        {"buffer", 0},
        {"byteOffset", info.colors.offset},
        {"byteLength", info.colors.length},
        {"target", 34962},
    });

    /* Buffer views: indices */
    int bv_idx = buffer_view_idx++;
    buffer_views.push_back({
        {"buffer", 0},
        {"byteOffset", info.indices.offset},
        {"byteLength", info.indices.length},
        {"target", 34963}, /* ELEMENT_ARRAY_BUFFER */
    });

    /* Accessor: positions */
    int acc_pos = accessor_idx++;
    accessors.push_back({
        {"bufferView", bv_pos},
        {"componentType", 5126}, /* FLOAT */
        {"count", mesh.vertex_count},
        {"type", "VEC3"},
        {"min", {mesh.min_pos[0], mesh.min_pos[1], mesh.min_pos[2]}},
        {"max", {mesh.max_pos[0], mesh.max_pos[1], mesh.max_pos[2]}},
    });

    /* Accessor: normals */
    int acc_norm = accessor_idx++;
    accessors.push_back({
        {"bufferView", bv_norm},
        {"componentType", 5126},
        {"count", mesh.vertex_count},
        {"type", "VEC3"},
    });

    /* Accessor: colors */
    int acc_col = accessor_idx++;
    accessors.push_back({
        {"bufferView", bv_col},
        {"componentType", 5126},
        {"count", mesh.vertex_count},
        {"type", "VEC3"},
    });

    /* Accessor: indices */
    int acc_idx = accessor_idx++;
    std::uint32_t max_index =
        *std::max_element(mesh.indices.begin(), mesh.indices.end());
    accessors.push_back({
        {"bufferView", bv_idx},
        {"componentType", 5125}, /* UNSIGNED_INT */
        {"count", mesh.indices.size()},
        {"type", "SCALAR"},
        {"min", {0}},
        {"max", {max_index}},
    });

    /* Mesh */
    json primitive;
    primitive["attributes"] = {
        {"POSITION", acc_pos},
        {"NORMAL", acc_norm},
        {"COLOR_0", acc_col},
    };
    primitive["indices"] = acc_idx;
    primitive["material"] = mat_idx;
    primitive["mode"] = 4; /* TRIANGLES */

    json mesh_obj;
    mesh_obj["primitives"] = json::array({primitive});
    int mesh_idx = static_cast<int>(meshes.size());
    meshes.push_back(mesh_obj);

    /* Node */
    nodes.push_back({{"mesh", mesh_idx}});
  }

  /* Scene: reference all nodes */
  json scene;
  json children = json::array();
  for (size_t i = 0; i < nodes.size(); i++) {
    children.push_back(static_cast<int>(i));
  }
  scene["nodes"] = children;
  root["scenes"] = json::array({scene});

  root["nodes"] = nodes;
  root["meshes"] = meshes;
  root["materials"] = materials;
  root["accessors"] = accessors;
  root["bufferViews"] = buffer_views;

  /* Buffer */
  root["buffers"] = json::array({{{"byteLength", buffer_size}}});

  return root.dump();
}

void RayRenderGLB(CRay* I, int width, int height, char** vla_ptr, float front,
    float back, float fov)
{
  PyMOLGlobals* G = I->G;

  /* Ray trace - expand all objects into primitives */
  RayExpandPrimitives(I);
  RayTransformFirst(I, 0, false);

  /* Collect primitives into mesh groups by transparency */
  std::vector<MeshGroup> groups;

  for (int a = 0; a < I->NPrimitive; a++) {
    CPrimitive* prim = I->Primitive + a;

    MeshGroup& group = findOrCreateGroup(groups, prim->trans);

    switch (prim->type) {
    case cPrimTriangle:
      addTriangle(group, prim);
      break;

    case cPrimSphere:
      addSphere(group, prim, G);
      break;

    case cPrimCylinder:
    case cPrimSausage:
    case cPrimCone:
      addCylinderWithCaps(group, prim, G);
      break;

    case cPrimCharacter:
    case cPrimEllipsoid:
      /* Not supported - skip silently (same as COLLADA) */
      break;
    }
  }

  /* Remove empty groups */
  groups.erase(std::remove_if(groups.begin(), groups.end(),
                   [](const MeshGroup& g) { return g.empty(); }),
      groups.end());

  if (groups.empty()) {
    PRINTFB(G, FB_Ray, FB_Warnings)
    " GLB-Warning: No geometry to export.\n" ENDFB(G);
    VLASize(*vla_ptr, char, 0);
    return;
  }

  /* Lay out the binary buffer without materializing it yet */
  std::vector<GroupBufferInfo> buffer_infos;
  std::uint64_t bin_len = computeBufferLayout(groups, buffer_infos);
  std::uint64_t bin_padded = alignTo4(bin_len);

  if (bin_padded > GLB_MAX_SIZE) {
    PRINTFB(G, FB_Ray, FB_Errors)
    " GLB-Error: Geometry is too large for the GLB container "
    "(%llu bytes, limit is %llu).\n",
        (unsigned long long) bin_padded,
        (unsigned long long) GLB_MAX_SIZE ENDFB(G);
    VLASize(*vla_ptr, char, 0);
    return;
  }

  /* Build JSON */
  std::string json_str =
      buildGLTFJson(groups, buffer_infos, static_cast<std::uint32_t>(bin_len));

  /* Pad JSON to 4-byte alignment with spaces */
  std::uint64_t json_len = json_str.size();
  std::uint64_t json_padded = alignTo4(json_len);

  std::uint64_t total = GLB_HEADER_SIZE + GLB_CHUNK_HEADER_SIZE + json_padded +
                        GLB_CHUNK_HEADER_SIZE + bin_padded;

  if (total > GLB_MAX_SIZE) {
    PRINTFB(G, FB_Ray, FB_Errors)
    " GLB-Error: Scene is too large for the GLB container "
    "(%llu bytes, limit is %llu).\n",
        (unsigned long long) total, (unsigned long long) GLB_MAX_SIZE ENDFB(G);
    VLASize(*vla_ptr, char, 0);
    return;
  }

  std::uint32_t total_size = static_cast<std::uint32_t>(total);

  /* Assemble the GLB directly in the output VLA */
  char* vla = *vla_ptr;

  VLASize(vla, char, total_size);

  char* out = vla;

  /* Header */
  out = writeU32(out, GLB_MAGIC);
  out = writeU32(out, GLB_VERSION);
  out = writeU32(out, total_size);

  /* JSON chunk */
  out = writeU32(out, static_cast<std::uint32_t>(json_padded));
  out = writeU32(out, GLB_CHUNK_JSON);
  std::memcpy(out, json_str.data(), json_len);
  std::memset(out + json_len, ' ', json_padded - json_len);
  out += json_padded;

  /* BIN chunk */
  out = writeU32(out, static_cast<std::uint32_t>(bin_padded));
  out = writeU32(out, GLB_CHUNK_BIN);
  writeBinaryBuffer(groups, buffer_infos, out);
  std::memset(out + bin_len, 0, bin_padded - bin_len);

  *vla_ptr = vla;

  PRINTFB(G, FB_Ray, FB_Actions)
  " GLB: exported %d mesh group(s), %u bytes total.\n", (int) groups.size(),
      total_size ENDFB(G);
}
