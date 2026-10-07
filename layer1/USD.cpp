/*
 * PyMOL OpenUSD export
 *
 * The renderer writes an ASCII USD layer. USDZ packaging is handled by
 * modules/pymol/exporting.py.
 */

#include "Ray.h"

#include "Basis.h"
#include "Color.h"
#include "MemoryDebug.h"
#include "Setting.h"
#include "Vector.h"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <initializer_list>
#include <iomanip>
#include <map>
#include <ostream>
#include <streambuf>
#include <string>
#include <unordered_map>
#include <vector>

namespace
{

constexpr float USD_EPSILON = 1.e-6F;

/// Extent of a collapsed ellipsoid semi-axis, relative to the largest one
constexpr float USD_DEGENERATE_AXIS_SCALE = 1.e-3F;

/**
 * Radial segments of a tessellated cylinder or cone.
 *
 * The ray tracer draws these solids analytically, so PyMOL's quality settings
 * can only raise the exported resolution above the floor.
 */
constexpr int USD_SEGMENTS_MIN = 24;
constexpr int USD_SEGMENTS_MAX = 100;

/// Distance below which two solid ends count as one joint
constexpr float USD_JOINT_TOLERANCE = 1.e-3F;

/**
 * Stream buffer which appends to a char VLA, so that the scene is only ever
 * held once in memory.
 */
class UsdVLAStreamBuf : public std::streambuf
{
public:
  UsdVLAStreamBuf(char** vla_ptr, ov_size* size_ptr)
      : m_vla_ptr(vla_ptr)
      , m_size_ptr(size_ptr)
  {
    setp(m_buffer, m_buffer + sizeof(m_buffer));
  }

  ~UsdVLAStreamBuf() override { Flush(); }

protected:
  int_type overflow(int_type ch) override
  {
    Flush();
    if (!traits_type::eq_int_type(ch, traits_type::eof())) {
      *pptr() = traits_type::to_char_type(ch);
      pbump(1);
    }
    return traits_type::not_eof(ch);
  }

  int sync() override
  {
    Flush();
    return 0;
  }

private:
  void Flush()
  {
    const auto pending = static_cast<ov_size>(pptr() - pbase());
    if (pending) {
      const ov_size size = *m_size_ptr;
      VLACheck(*m_vla_ptr, char, size + pending);
      std::memcpy(*m_vla_ptr + size, pbase(), pending);
      *m_size_ptr = size + pending;
      (*m_vla_ptr)[*m_size_ptr] = '\0';
    }
    setp(m_buffer, m_buffer + sizeof(m_buffer));
  }

  char** m_vla_ptr;
  ov_size* m_size_ptr;
  char m_buffer[1 << 16];
};

float UsdClamp(float value)
{
  return std::clamp(value, 0.F, 1.F);
}

void UsdWriteVec3(std::ostream& out, const float* value)
{
  out << '(' << value[0] << ", " << value[1] << ", " << value[2] << ')';
}

void UsdWriteColor(std::ostream& out, const float* color)
{
  out << '(' << UsdClamp(color[0]) << ", " << UsdClamp(color[1]) << ", "
      << UsdClamp(color[2]) << ')';
}

bool UsdColorsEqual(const float* lhs, const float* rhs)
{
  return std::fabs(lhs[0] - rhs[0]) < USD_EPSILON &&
         std::fabs(lhs[1] - rhs[1]) < USD_EPSILON &&
         std::fabs(lhs[2] - rhs[2]) < USD_EPSILON;
}

/**
 * One UsdPreviewSurface material per distinct color and opacity.
 *
 * AR Quick Look ignores primvar readers, so displayColor cannot drive the
 * material. Colors are quantized to keep smooth ramps and spectra to a few
 * hundred materials.
 */
class UsdMaterials
{
  static constexpr int COLOR_LEVELS = 63;
  static constexpr int OPACITY_LEVELS = 100;

public:
  /// Index of the material for the given color and opacity
  int Get(const float* color, float opacity)
  {
    std::uint32_t key = std::lround(UsdClamp(opacity) * OPACITY_LEVELS);
    for (int i = 0; i < 3; ++i) {
      key = key << 6 | std::lround(UsdClamp(color[i]) * COLOR_LEVELS);
    }

    const auto inserted =
        m_indices.emplace(key, static_cast<int>(m_keys.size()));
    if (inserted.second) {
      m_keys.push_back(key);
    }
    return inserted.first->second;
  }

  void GetColor(int material, float* color) const
  {
    const auto key = m_keys[material];
    for (int i = 0; i < 3; ++i) {
      color[i] = float(key >> (12 - 6 * i) & 63) / COLOR_LEVELS;
    }
  }

  float GetOpacity(int material) const
  {
    return float(m_keys[material] >> 18) / OPACITY_LEVELS;
  }

  static void WriteBinding(std::ostream& out, int material)
  {
    out << "rel material:binding = </PyMOLScene/Materials/Material_" << material
        << ">\n";
  }

  void Write(std::ostream& out) const
  {
    if (m_keys.empty()) {
      return;
    }

    out << "\n"
        << "    def Scope \"Materials\"\n"
        << "    {\n";

    for (int i = 0; i < static_cast<int>(m_keys.size()); ++i) {
      const float opacity = GetOpacity(i);
      float color[3];
      GetColor(i, color);
      const std::string path =
          "/PyMOLScene/Materials/Material_" + std::to_string(i);

      out << (i ? "\n" : "") << "        def Material \"Material_" << i
          << "\"\n"
          << "        {\n"
          << "            token outputs:surface.connect = <" << path
          << "/PreviewSurface.outputs:surface>\n"
          << "\n"
          << "            def Shader \"PreviewSurface\"\n"
          << "            {\n"
          << "                uniform token info:id = \"UsdPreviewSurface\"\n"
          << "                color3f inputs:diffuseColor = ";
      UsdWriteVec3(out, color);
      out << "\n";

      // Renderers treat any material with an opacity as translucent
      if (opacity < 1.F) {
        out << "                float inputs:opacity = " << opacity << "\n";
      }

      out << "                float inputs:roughness = 0.35\n"
          << "                token outputs:surface\n"
          << "            }\n"
          << "        }\n";
    }

    out << "    }\n";
  }

private:
  std::unordered_map<std::uint32_t, int> m_indices;
  std::vector<std::uint32_t> m_keys;
};

void UsdWriteMaterialBinding(std::ostream& out, UsdMaterials& materials,
    const float* color, float transparency)
{
  out << "        ";
  UsdMaterials::WriteBinding(out, materials.Get(color, 1.F - transparency));
  out << "        color3f[] primvars:displayColor = [";
  UsdWriteColor(out, color);
  out << "] (\n"
      << "            interpolation = \"constant\"\n"
      << "        )\n"
      << "        float[] primvars:displayOpacity = ["
      << UsdClamp(1.F - transparency) << "] (\n"
      << "            interpolation = \"constant\"\n"
      << "        )\n";
}

struct UsdQuaternion {
  float real;
  float imaginary[3];
};

UsdQuaternion UsdRotationFromZAxis(const float* direction)
{
  // Quaternion which rotates the USD analytic primitives' +Z axis onto
  // direction. For the anti-parallel case, rotate 180 degrees around X.
  const float dot = direction[2];
  if (dot < -1.F + USD_EPSILON) {
    return {0.F, {1.F, 0.F, 0.F}};
  }

  UsdQuaternion result{1.F + dot, {-direction[1], direction[0], 0.F}};
  const float length = std::sqrt(result.real * result.real +
                                 result.imaginary[0] * result.imaginary[0] +
                                 result.imaginary[1] * result.imaginary[1]);

  if (length < USD_EPSILON) {
    return {1.F, {0.F, 0.F, 0.F}};
  }

  result.real /= length;
  result.imaginary[0] /= length;
  result.imaginary[1] /= length;
  return result;
}

void UsdWriteAnalyticTransform(
    std::ostream& out, const float* start, const float* end, float& height)
{
  float direction[3] = {
      end[0] - start[0], end[1] - start[1], end[2] - start[2]};
  height = std::sqrt(direction[0] * direction[0] + direction[1] * direction[1] +
                     direction[2] * direction[2]);

  if (height > USD_EPSILON) {
    direction[0] /= height;
    direction[1] /= height;
    direction[2] /= height;
  } else {
    direction[0] = 0.F;
    direction[1] = 0.F;
    direction[2] = 1.F;
  }

  const float center[3] = {(start[0] + end[0]) * 0.5F,
      (start[1] + end[1]) * 0.5F, (start[2] + end[2]) * 0.5F};
  const auto orientation = UsdRotationFromZAxis(direction);

  out << "        float3 xformOp:translate = ";
  UsdWriteVec3(out, center);
  out << "\n"
      << "        quatf xformOp:orient = (" << orientation.real << ", "
      << orientation.imaginary[0] << ", " << orientation.imaginary[1] << ", "
      << orientation.imaginary[2] << ")\n"
      << "        uniform token[] xformOpOrder = [\"xformOp:translate\", "
         "\"xformOp:orient\"]\n";
}

void UsdWriteSphere(std::ostream& out, int& index, UsdMaterials& materials,
    const float* center, float radius, const float* color, float transparency)
{
  out << "\n"
      << "    def Sphere \"Sphere_" << index++ << "\" (\n"
      << "        prepend apiSchemas = [\"MaterialBindingAPI\"]\n"
      << "    )\n"
      << "    {\n"
      << "        double radius = " << radius << "\n"
      << "        float3[] extent = [(" << -radius << ", " << -radius << ", "
      << -radius << "), (" << radius << ", " << radius << ", " << radius
      << ")]\n"
      << "        float3 xformOp:translate = ";
  UsdWriteVec3(out, center);
  out << "\n"
      << "        uniform token[] xformOpOrder = [\"xformOp:translate\"]\n";
  UsdWriteMaterialBinding(out, materials, color, transparency);
  out << "    }\n";
}

/**
 * Write a Z-aligned UsdGeomCylinder, UsdGeomCone or UsdGeomCapsule spanning
 * start to end.
 *
 * A capsule's height covers its cylindrical spine only, and its two
 * hemispherical caps reach one radius further along the axis, which is
 * exactly how PyMOL draws a round capped solid.
 *
 * @param schema "Cylinder", "Cone" or "Capsule", also the prim name prefix
 * @param round_caps whether the schema extends one radius beyond each end
 */
void UsdWriteAnalyticSolid(std::ostream& out, int& index,
    UsdMaterials& materials, const char* schema, bool round_caps,
    const float* start, const float* end, float radius, const float* color,
    float transparency)
{
  float height;

  out << "\n"
      << "    def " << schema << " \"" << schema << '_' << index++ << "\" (\n"
      << "        prepend apiSchemas = [\"MaterialBindingAPI\"]\n"
      << "    )\n"
      << "    {\n";
  UsdWriteAnalyticTransform(out, start, end, height);

  const float extent_z = height * 0.5F + (round_caps ? radius : 0.F);

  out << "        uniform token axis = \"Z\"\n"
      << "        double height = " << height << "\n"
      << "        double radius = " << radius << "\n"
      << "        float3[] extent = [(" << -radius << ", " << -radius << ", "
      << -extent_z << "), (" << radius << ", " << radius << ", " << extent_z
      << ")]\n";
  UsdWriteMaterialBinding(out, materials, color, transparency);
  out << "    }\n";
}

/**
 * Vertex and face data of one exported UsdGeomMesh.
 *
 * Tessellated solids buffer their vertices, while the scene mesh reads them
 * straight from the ray tracer's arrays, so the writer pulls one vertex at a
 * time instead of taking flat arrays.
 */
struct UsdMeshSource {
  virtual ~UsdMeshSource() = default;
  virtual std::size_t VertexCount() const = 0;
  virtual std::size_t FaceCount() const = 0;
  /// Number of corners of the given face
  virtual int FaceSize(std::size_t face) const = 0;
  /// Vertex of the given corner, counted over all faces
  virtual int Index(std::size_t corner) const = 0;
  virtual const float* Point(std::size_t vertex) const = 0;
  virtual const float* Normal(std::size_t vertex) const = 0;
  virtual const float* Color(std::size_t vertex) const = 0;
  virtual float Opacity(std::size_t vertex) const = 0;
};

/// Mesh whose vertices and faces are buffered in vectors
struct UsdVectorMesh : UsdMeshSource {
  std::vector<float> points;
  std::vector<float> normals;
  std::vector<float> colors;
  std::vector<float> opacities;
  std::vector<int> counts;
  std::vector<int> indices;

  /// Opacity of the vertices added next
  float opacity = 1.F;

  std::size_t VertexCount() const override { return points.size() / 3; }
  std::size_t FaceCount() const override { return counts.size(); }
  int FaceSize(std::size_t face) const override { return counts[face]; }
  int Index(std::size_t corner) const override { return indices[corner]; }

  float Opacity(std::size_t vertex) const override { return opacities[vertex]; }

  const float* Point(std::size_t vertex) const override
  {
    return points.data() + 3 * vertex;
  }

  const float* Normal(std::size_t vertex) const override
  {
    return normals.data() + 3 * vertex;
  }

  const float* Color(std::size_t vertex) const override
  {
    return colors.data() + 3 * vertex;
  }

  void AddVertex(const float* point, const float* normal, const float* color)
  {
    points.insert(points.end(), point, point + 3);
    normals.insert(normals.end(), normal, normal + 3);
    colors.insert(colors.end(), color, color + 3);
    opacities.push_back(opacity);
  }

  void AddFace(std::initializer_list<int> face)
  {
    counts.push_back(static_cast<int>(face.size()));
    indices.insert(indices.end(), face.begin(), face.end());
  }
};

/// Write one UsdGeomMesh, shared by all exported mesh kinds
void UsdWriteMesh(std::ostream& out, int& index, UsdMaterials& materials,
    const char* name, const UsdMeshSource& mesh)
{
  const auto vertex_count = mesh.VertexCount();
  const auto face_count = mesh.FaceCount();

  if (!vertex_count || !face_count) {
    return;
  }

  float extent_min[3];
  float extent_max[3];

  copy3f(mesh.Point(0), extent_min);
  copy3f(mesh.Point(0), extent_max);

  for (std::size_t vertex = 1; vertex < vertex_count; ++vertex) {
    const auto* point = mesh.Point(vertex);

    for (int axis = 0; axis < 3; ++axis) {
      extent_min[axis] = std::min(extent_min[axis], point[axis]);
      extent_max[axis] = std::max(extent_max[axis], point[axis]);
    }
  }

  out << "\n"
      << "    def Mesh \"" << name << '_' << index++ << "\" (\n"
      << "        prepend apiSchemas = [\"MaterialBindingAPI\"]\n"
      << "    )\n"
      << "    {\n"
      << "        float3[] extent = [";
  UsdWriteVec3(out, extent_min);
  out << ", ";
  UsdWriteVec3(out, extent_max);
  out << "]\n"
      << "        int[] faceVertexCounts = [";

  // Material of each vertex, and the materials used by the mesh, whose
  // colors and opacities make up the indexed displayColor and displayOpacity
  std::vector<int> vertex_materials(vertex_count);
  std::map<int, int> palette;
  for (std::size_t vertex = 0; vertex < vertex_count; ++vertex) {
    vertex_materials[vertex] =
        materials.Get(mesh.Color(vertex), mesh.Opacity(vertex));
    palette.emplace(vertex_materials[vertex], 0);
  }

  int palette_size = 0;
  for (auto& entry : palette) {
    entry.second = palette_size++;
  }

  // Faces of each material, which takes the color of the first corner. An
  // average would leave the path of a color ramp and multiply materials.
  std::map<int, std::vector<int>> material_faces;
  std::size_t corner_count = 0;
  for (std::size_t face = 0; face < face_count; ++face) {
    const int size = mesh.FaceSize(face);
    material_faces[vertex_materials[mesh.Index(corner_count)]].push_back(face);

    out << (face ? ", " : "") << size;
    corner_count += size;
  }

  out << "]\n"
      << "        int[] faceVertexIndices = [";
  for (std::size_t corner = 0; corner < corner_count; ++corner) {
    out << (corner ? ", " : "") << mesh.Index(corner);
  }

  out << "]\n"
      << "        point3f[] points = [\n";
  for (std::size_t vertex = 0; vertex < vertex_count; ++vertex) {
    out << "            ";
    UsdWriteVec3(out, mesh.Point(vertex));
    out << ",\n";
  }

  out << "        ]\n"
      << "        normal3f[] normals = [\n";
  for (std::size_t vertex = 0; vertex < vertex_count; ++vertex) {
    out << "            ";
    UsdWriteVec3(out, mesh.Normal(vertex));
    out << ",\n";
  }

  out << "        ] (\n"
      << "            interpolation = \"vertex\"\n"
      << "        )\n"
      << "        color3f[] primvars:displayColor = [";
  for (const auto& [material, entry] : palette) {
    float color[3];
    materials.GetColor(material, color);
    out << (entry ? ", " : "");
    UsdWriteVec3(out, color);
  }

  out << "] (\n"
      << "            interpolation = \"vertex\"\n"
      << "        )\n"
      << "        int[] primvars:displayColor:indices = [";
  for (std::size_t vertex = 0; vertex < vertex_count; ++vertex) {
    out << (vertex ? ", " : "") << palette[vertex_materials[vertex]];
  }

  out << "]\n"
      << "        float[] primvars:displayOpacity = [";
  for (const auto& [material, entry] : palette) {
    out << (entry ? ", " : "") << materials.GetOpacity(material);
  }

  out << "] (\n"
      << "            interpolation = \"vertex\"\n"
      << "        )\n"
      << "        int[] primvars:displayOpacity:indices = [";
  for (std::size_t vertex = 0; vertex < vertex_count; ++vertex) {
    out << (vertex ? ", " : "") << palette[vertex_materials[vertex]];
  }

  out << "]\n"
      << "        uniform token subdivisionScheme = \"none\"\n"
      // The ray tracer flips the normal of a back face hit, so it draws an
      // open mesh from both sides, while USD culls back faces by default
      << "        uniform bool doubleSided = 1\n";

  if (material_faces.size() == 1) {
    out << "        ";
    UsdMaterials::WriteBinding(out, material_faces.begin()->first);
  } else {
    out << "        uniform token subsetFamily:materialBind:familyType = "
           "\"partition\"\n";

    for (const auto& [material, faces] : material_faces) {
      out << "\n"
          << "        def GeomSubset \"Material_" << material << "\" (\n"
          << "            prepend apiSchemas = [\"MaterialBindingAPI\"]\n"
          << "        )\n"
          << "        {\n"
          << "            uniform token elementType = \"face\"\n"
          << "            uniform token familyName = \"materialBind\"\n"
          << "            int[] indices = [";
      for (std::size_t i = 0; i < faces.size(); ++i) {
        out << (i ? ", " : "") << faces[i];
      }
      out << "]\n"
          << "            ";
      UsdMaterials::WriteBinding(out, material);
      out << "        }\n";
    }
  }

  out << "    }\n";
}

/**
 * Add a tessellated cylinder, cone or truncated cone (frustum).
 *
 * @param segments radial segments, see UsdSolidSegments
 * @param start center of the r1 end
 * @param end center of the r2 end
 * @param cap1 cap of the r1 end, only a flat cap closes the mesh
 * @param cap2 cap of the r2 end
 */
void UsdAddCone(UsdVectorMesh& mesh, int segments, const float* start,
    const float* end, float r1, float r2, const float* c1, const float* c2,
    cCylCap cap1, cCylCap cap2, float opacity)
{
  float axis[3];
  float side[3];
  float up[3];

  subtract3f(end, start, axis);
  const float height = length3f(axis);

  if (!(height > USD_EPSILON) || !(r1 > USD_EPSILON)) {
    return;
  }

  // The normals below assume a unit axis
  normalize3f(axis);
  get_system1f3f(axis, side, up);

  assert(segments >= 3 && segments <= USD_SEGMENTS_MAX);

  const bool pointed = !(r2 > USD_EPSILON);
  float radial[3 * USD_SEGMENTS_MAX];
  float slope[3 * USD_SEGMENTS_MAX];

  for (int k = 0; k < segments; ++k) {
    const auto angle = static_cast<float>(2.0 * PI * k / segments);
    const float cosine = std::cos(angle);
    const float sine = std::sin(angle);

    for (int i = 0; i < 3; ++i) {
      radial[3 * k + i] = side[i] * cosine + up[i] * sine;
      slope[3 * k + i] = radial[3 * k + i] * height + axis[i] * (r1 - r2);
    }

    normalize3f(slope + 3 * k);
  }

  mesh.opacity = opacity;

  auto add_ring = [&](const float* center, float radius, const float* normal,
                      const float* color) {
    for (int k = 0; k < segments; ++k) {
      float point[3];
      scale3f(radial + 3 * k, radius, point);
      add3f(center, point, point);
      mesh.AddVertex(point, normal ? normal : slope + 3 * k, color);
    }
  };

  // Lateral faces between two rings, winding counter-clockwise as seen from
  // outside, so that the right handed face normals agree with the authored
  // ones
  auto add_band = [&](const float* p1, float ra, const float* p2, float rb,
                      const float* color) {
    const auto base = static_cast<int>(mesh.VertexCount());
    add_ring(p1, ra, nullptr, color);
    add_ring(p2, rb, nullptr, color);

    for (int k = 0; k < segments; ++k) {
      const int a = base + k;
      const int b = base + (k + 1) % segments;
      if (rb > USD_EPSILON) {
        mesh.AddFace({a, b, b + segments, a + segments});
      } else {
        mesh.AddFace({a, b, a + segments});
      }
    }
  };

  // Each face takes one material, so two colors meet halfway like PyMOL's
  // half bonds, rather than blending along the axis like the ray tracer
  if (UsdColorsEqual(c1, c2)) {
    add_band(start, r1, end, r2, c1);
  } else {
    float middle[3];
    average3f(start, end, middle);
    const float rm = (r1 + r2) * 0.5F;
    add_band(start, r1, middle, rm, c1);
    add_band(middle, rm, end, r2, c2);
  }

  if (cap1 == cCylCapFlat) {
    const auto base = static_cast<int>(mesh.VertexCount());
    float normal[3];

    invert3f3f(axis, normal);
    mesh.AddVertex(start, normal, c1);
    add_ring(start, r1, normal, c1);

    for (int k = 0; k < segments; ++k) {
      const int next = (k + 1) % segments;
      mesh.AddFace({base, base + 1 + next, base + 1 + k});
    }
  }

  if (cap2 == cCylCapFlat && !pointed) {
    const auto base = static_cast<int>(mesh.VertexCount());

    mesh.AddVertex(end, axis, c2);
    add_ring(end, r2, axis, c2);

    for (int k = 0; k < segments; ++k) {
      const int next = (k + 1) % segments;
      mesh.AddFace({base, base + 1 + k, base + 1 + next});
    }
  }
}

/**
 * Add a tessellated hemisphere, e.g. the round cap of a solid.
 *
 * As a cap it shares the barrel's equator ring, which UsdAddCone leaves open
 * for a round cap, so no buried surface darkens a transparent solid.
 *
 * @param segments radial segments, see UsdSolidSegments
 * @param center center of the capped end
 * @param axis unit axis of the solid, pointing from start to end
 * @param sign +1 for the dome beyond end, -1 for the dome beyond start
 */
void UsdAddHemisphere(UsdVectorMesh& mesh, int segments, const float* center,
    const float* axis, float sign, float radius, const float* color,
    float opacity)
{
  if (!(radius > USD_EPSILON)) {
    return;
  }

  float unit[3];
  float side[3];
  float up[3];
  float pole[3];

  // get_system1f3f normalizes its first argument in place
  copy3f(axis, unit);
  get_system1f3f(unit, side, up);

  // (side, up, pole) stays right handed for either end, so the winding below
  // holds for both, and the equator lands on the barrel's own ring vertices
  scale3f(up, sign, up);
  scale3f(unit, sign, pole);

  assert(segments >= 3 && segments <= USD_SEGMENTS_MAX);

  const int stacks = std::max(2, segments / 4);
  const auto base = static_cast<int>(mesh.VertexCount());

  mesh.opacity = opacity;

  for (int j = 0; j < stacks; ++j) {
    const auto latitude = static_cast<float>(0.5 * PI * j / stacks);
    const float ring = std::cos(latitude);
    const float rise = std::sin(latitude);

    for (int k = 0; k < segments; ++k) {
      const auto angle = static_cast<float>(2.0 * PI * k / segments);
      const float cosine = std::cos(angle);
      const float sine = std::sin(angle);
      float normal[3];
      float point[3];

      for (int i = 0; i < 3; ++i) {
        normal[i] = (side[i] * cosine + up[i] * sine) * ring + pole[i] * rise;
      }

      scale3f(normal, radius, point);
      add3f(center, point, point);
      mesh.AddVertex(point, normal, color);
    }
  }

  const int apex = base + stacks * segments;
  float tip[3];

  scale3f(pole, radius, tip);
  add3f(center, tip, tip);
  mesh.AddVertex(tip, pole, color);

  // Faces wind counter-clockwise as seen from outside, as in UsdAddCone
  for (int j = 0; j + 1 < stacks; ++j) {
    const int ring = base + j * segments;
    for (int k = 0; k < segments; ++k) {
      const int next = (k + 1) % segments;
      mesh.AddFace(
          {ring + k, ring + next, ring + segments + next, ring + segments + k});
    }
  }

  const int top = base + (stacks - 1) * segments;
  for (int k = 0; k < segments; ++k) {
    mesh.AddFace({top + k, top + (k + 1) % segments, apex});
  }
}

/// Add a tessellated sphere, as two hemispheres
void UsdAddSphere(UsdVectorMesh& mesh, int segments, const float* center,
    float radius, const float* color, float opacity)
{
  const float axis[3] = {0.F, 0.F, 1.F};
  UsdAddHemisphere(mesh, segments, center, axis, 1.F, radius, color, opacity);
  UsdAddHemisphere(mesh, segments, center, axis, -1.F, radius, color, opacity);
}

/**
 * Build the upper 3x3 rows of an ellipsoid's transform.
 *
 * Semi-axes can legitimately collapse: a non-positive-definite ANISOU tensor
 * makes RepEllipsoid drop an eigenvector, and ellipsoid_scale = 0 zeroes the
 * radius. A zero row would make the emitted matrix singular, which renderers
 * reject, so collapsed axes get a small deterministic extent instead.
 *
 * @param axes three consecutive transformed unit axes
 * @param scales relative axis scales, one per axis, largest of which is 1
 * @param radius largest semi-axis length
 * @param[out] rows three consecutive scaled, right-handed axes
 * @return false if the ellipsoid has no visible extent and must be skipped
 */
bool UsdEllipsoidRows(
    const float* axes, const float* scales, float radius, float* rows)
{
  float directions[9];
  float lengths[3];
  float max_length = 0.F;
  int valid[3];
  int valid_count = 0;

  for (int axis = 0; axis < 3; ++axis) {
    auto* direction = directions + 3 * axis;
    const float length = radius * scales[axis];

    copy3f(axes + 3 * axis, direction);
    lengths[axis] = 0.F;

    if (!(length > 0.F) || !std::isfinite(length) ||
        !(length3f(direction) > USD_EPSILON)) {
      continue;
    }

    normalize3f(direction);
    lengths[axis] = length;
    max_length = std::max(max_length, length);
    valid[valid_count++] = axis;
  }

  if (max_length <= USD_EPSILON) {
    return false;
  }

  // Complete the collapsed axes into an orthonormal frame, so that the
  // fallback extents point somewhere meaningful.
  if (valid_count == 1) {
    const int first = valid[0];
    get_system1f3f(directions + 3 * first, directions + 3 * ((first + 1) % 3),
        directions + 3 * ((first + 2) % 3));
  } else if (valid_count == 2) {
    const int missing = 3 - valid[0] - valid[1];
    cross_product3f(directions + 3 * valid[0], directions + 3 * valid[1],
        directions + 3 * missing);
    normalize3f(directions + 3 * missing);
  }

  const float fallback = max_length * USD_DEGENERATE_AXIS_SCALE;

  for (int axis = 0; axis < 3; ++axis) {
    if (!(lengths[axis] >= fallback)) {
      lengths[axis] = fallback;
    }
    scale3f(directions + 3 * axis, lengths[axis], rows + 3 * axis);
  }

  const double determinant = determinant33f(rows);

  // The rows are unit axes scaled by lengths, so their product is the
  // determinant of a perfectly orthogonal frame and makes the test relative
  const double orthogonal = static_cast<double>(lengths[0]) *
                            static_cast<double>(lengths[1]) *
                            static_cast<double>(lengths[2]);

  if (!(std::fabs(determinant) > USD_EPSILON * orthogonal)) {
    // Axes which are not linearly independent, e.g. after a shearing object
    // matrix. An axis aligned frame keeps the transform invertible.
    for (int axis = 0; axis < 3; ++axis) {
      zero3f(rows + 3 * axis);
      rows[3 * axis + axis] = lengths[axis];
    }
  } else if (determinant < 0.0) {
    // The axes are eigenvectors with arbitrary signs and may form a
    // left-handed basis. Negating one axis leaves the ellipsoid unchanged.
    for (int i = 6; i < 9; ++i) {
      rows[i] = -rows[i];
    }
  }

  return true;
}

/**
 * @param center transformed ellipsoid center
 * @copydetails UsdEllipsoidRows
 */
void UsdWriteEllipsoid(std::ostream& out, int& index, UsdMaterials& materials,
    const float* center, const float* axes, const float* scales, float radius,
    const float* color, float transparency)
{
  float rows[9];

  if (!UsdEllipsoidRows(axes, scales, radius, rows)) {
    return;
  }

  out << "\n"
      << "    def Sphere \"Ellipsoid_" << index++ << "\" (\n"
      << "        prepend apiSchemas = [\"MaterialBindingAPI\"]\n"
      << "    )\n"
      << "    {\n"
      << "        double radius = 1\n"
      << "        float3[] extent = [(-1, -1, -1), (1, 1, 1)]\n"
      << "        matrix4d xformOp:transform = (";
  for (int axis = 0; axis < 3; ++axis) {
    const auto* row = rows + 3 * axis;
    out << "\n            (" << row[0] << ", " << row[1] << ", " << row[2]
        << ", 0),";
  }
  out << "\n            (" << center[0] << ", " << center[1] << ", "
      << center[2] << ", 1)"
      << "\n        )\n"
      << "        uniform token[] xformOpOrder = [\"xformOp:transform\"]\n";
  UsdWriteMaterialBinding(out, materials, color, transparency);
  out << "    }\n";
}

/**
 * Add a tessellated ellipsoid, a unit sphere transformed by the rows of
 * UsdEllipsoidRows.
 *
 * @copydetails UsdWriteEllipsoid
 */
void UsdAddEllipsoid(UsdVectorMesh& mesh, int segments, const float* center,
    const float* axes, const float* scales, float radius, const float* color,
    float opacity)
{
  float rows[9];

  if (!UsdEllipsoidRows(axes, scales, radius, rows)) {
    return;
  }

  // Normals transform with the inverse transpose, whose rows are the cross
  // products of the other two rows, up to the positive determinant
  float normal_rows[9];
  cross_product3f(rows + 3, rows + 6, normal_rows);
  cross_product3f(rows + 6, rows, normal_rows + 3);
  cross_product3f(rows, rows + 3, normal_rows + 6);

  const float origin[3] = {0.F, 0.F, 0.F};
  const auto base = mesh.VertexCount();
  UsdAddSphere(mesh, segments, origin, 1.F, color, opacity);

  for (auto vertex = base; vertex < mesh.VertexCount(); ++vertex) {
    float* point = mesh.points.data() + 3 * vertex;
    float* normal = mesh.normals.data() + 3 * vertex;
    float unit[3];
    copy3f(point, unit);

    for (int i = 0; i < 3; ++i) {
      point[i] = center[i] + unit[0] * rows[i] + unit[1] * rows[3 + i] +
                 unit[2] * rows[6 + i];
      normal[i] = unit[0] * normal_rows[i] + unit[1] * normal_rows[3 + i] +
                  unit[2] * normal_rows[6 + i];
    }

    normalize3f(normal);
  }
}

/**
 * The scene's triangles, read straight from the ray tracer's arrays.
 *
 * Only the primitive index and the resolved winding are kept per triangle, so
 * that a large surface does not need a second copy of its geometry.
 */
class UsdSceneTriangles : public UsdMeshSource
{
public:
  UsdSceneTriangles(const CRay* ray, const CBasis* basis)
      : m_ray(ray)
      , m_basis(basis)
  {
    for (int i = 0; i < ray->NPrimitive; ++i) {
      auto& primitive = ray->Primitive[i];
      if (primitive.type != cPrimTriangle) {
        continue;
      }

      m_primitives.push_back(i);
      m_reverse.push_back(TriangleReverse(&primitive) != 0);
    }
  }

  std::size_t VertexCount() const override { return 3 * m_primitives.size(); }
  std::size_t FaceCount() const override { return m_primitives.size(); }
  int FaceSize(std::size_t) const override { return 3; }

  int Index(std::size_t corner) const override
  {
    return static_cast<int>(corner);
  }

  const float* Point(std::size_t vertex) const override
  {
    return m_basis->Vertex + 3 * (Primitive(vertex).vert + Corner(vertex));
  }

  const float* Normal(std::size_t vertex) const override
  {
    return m_basis->Normal + 3 * (m_basis->Vert2Normal[Primitive(vertex).vert] +
                                     1 + Corner(vertex));
  }

  const float* Color(std::size_t vertex) const override
  {
    const auto& primitive = Primitive(vertex);
    const float* colors[3] = {primitive.c1, primitive.c2, primitive.c3};
    return colors[Corner(vertex)];
  }

  float Opacity(std::size_t vertex) const override
  {
    return 1.F - Primitive(vertex).tr[Corner(vertex)];
  }

private:
  const CPrimitive& Primitive(std::size_t vertex) const
  {
    return m_ray->Primitive[m_primitives[vertex / 3]];
  }

  /// Corner of the primitive, in counter-clockwise winding order
  int Corner(std::size_t vertex) const
  {
    static constexpr int order[2][3] = {{0, 1, 2}, {0, 2, 1}};
    return order[m_reverse[vertex / 3]][vertex % 3];
  }

  const CRay* m_ray;
  const CBasis* m_basis;
  std::vector<int> m_primitives;
  std::vector<bool> m_reverse;
};

/// Far end of a cylinder, sausage or cone, the near end being basis->Vertex
void UsdSolidEnd(const CBasis* basis, const CPrimitive& primitive, float* end)
{
  const auto* vertex = basis->Vertex + 3 * primitive.vert;
  const auto* direction =
      basis->Normal + 3 * basis->Vert2Normal[primitive.vert];

  for (int i = 0; i < 3; ++i) {
    end[i] = vertex[i] + direction[i] * primitive.l1;
  }
}

/**
 * Spatial index of the ends of all exported solids.
 *
 * PyMOL leaves the joint between abutting solids uncapped, e.g. between the
 * two halves of a stick. Such a joint is sealed by its neighbour, so in the
 * AR layer it may keep the closed analytic prim, while an exposed uncapped
 * end needs an open mesh. The buried end disc only goes unnoticed while the
 * solid is opaque, so transparent solids always take the mesh.
 */
class UsdJointIndex
{
  struct Endpoint {
    float point[3];
    float radius;
    int primitive;
  };

public:
  UsdJointIndex(const CRay* ray, const CBasis* basis)
  {
    if (!AnyUncapped(ray)) {
      return;
    }

    for (int i = 0; i < ray->NPrimitive; ++i) {
      const auto& primitive = ray->Primitive[i];
      const auto* vertex = basis->Vertex + 3 * primitive.vert;
      float end[3];

      switch (primitive.type) {
      case cPrimSphere:
        Add(vertex, primitive.r1, i);
        break;
      case cPrimCylinder:
      case cPrimSausage:
      case cPrimCone:
        UsdSolidEnd(basis, primitive, end);
        Add(vertex, primitive.r1, i);
        Add(end, primitive.type == cPrimCone ? primitive.r2 : primitive.r1, i);
        break;
      }
    }
  }

  /// True if another solid seals a disc of the given radius at point
  bool IsCovered(const float* point, float radius, int primitive) const
  {
    for (int dx = -1; dx <= 1; ++dx) {
      for (int dy = -1; dy <= 1; ++dy) {
        for (int dz = -1; dz <= 1; ++dz) {
          const auto range = m_cells.equal_range(Key(point, dx, dy, dz));

          for (auto it = range.first; it != range.second; ++it) {
            const auto& endpoint = m_endpoints[it->second];

            if (endpoint.primitive != primitive &&
                endpoint.radius >= radius - USD_JOINT_TOLERANCE &&
                diffsq3f(endpoint.point, point) <=
                    USD_JOINT_TOLERANCE * USD_JOINT_TOLERANCE) {
              return true;
            }
          }
        }
      }
    }

    return false;
  }

private:
  static bool AnyUncapped(const CRay* ray)
  {
    for (int i = 0; i < ray->NPrimitive; ++i) {
      const auto& primitive = ray->Primitive[i];

      // a sausage is round capped, its cap fields are never assigned
      if (primitive.type == cPrimCylinder &&
          (primitive.cap1 == cCylCapNone || primitive.cap2 == cCylCapNone)) {
        return true;
      }

      // the ray tracer only draws a flat cone cap
      if (primitive.type == cPrimCone && primitive.cap1 != cCylCapFlat) {
        return true;
      }
    }

    return false;
  }

  static std::int64_t Key(const float* point, int dx, int dy, int dz)
  {
    const std::int64_t cell[3] = {
        static_cast<std::int64_t>(std::floor(point[0] / USD_JOINT_TOLERANCE)) +
            dx,
        static_cast<std::int64_t>(std::floor(point[1] / USD_JOINT_TOLERANCE)) +
            dy,
        static_cast<std::int64_t>(std::floor(point[2] / USD_JOINT_TOLERANCE)) +
            dz};

    return (cell[0] * 73856093) ^ (cell[1] * 19349663) ^ (cell[2] * 83492791);
  }

  void Add(const float* point, float radius, int primitive)
  {
    if (!(radius > USD_EPSILON)) {
      return;
    }

    m_cells.emplace(Key(point, 0, 0, 0), static_cast<int>(m_endpoints.size()));
    m_endpoints.push_back({{point[0], point[1], point[2]}, radius, primitive});
  }

  std::vector<Endpoint> m_endpoints;
  std::unordered_multimap<std::int64_t, int> m_cells;
};

/**
 * Radial segments to tessellate a solid with, from the quality setting which
 * PyMOL's own renderers use for that kind of solid.
 */
int UsdSolidSegments(PyMOLGlobals* G, int setting)
{
  return std::clamp(
      SettingGetGlobal_i(G, setting), USD_SEGMENTS_MIN, USD_SEGMENTS_MAX);
}

/**
 * Replace ramp colors with the ramp's color at each primitive vertex.
 *
 * A ramp color is a negative color index which the ray tracer resolves per
 * hit point. The exporter has no hit points, so it resolves the ramp at the
 * model space vertices instead.
 */
void UsdResolveRampedColors(CRay* ray)
{
  for (int i = 0; i < ray->NPrimitive; ++i) {
    auto& primitive = ray->Primitive[i];
    if (!primitive.ramped) {
      continue;
    }

    float* colors[3] = {primitive.c1, primitive.c2, primitive.c3};
    const float* points[3] = {primitive.v1, primitive.v2, primitive.v3};
    int count = 1;

    switch (primitive.type) {
    case cPrimTriangle:
      count = 3;
      break;
    case cPrimCylinder:
    case cPrimSausage:
    case cPrimCone:
      count = 2;
      break;
    }

    for (int j = 0; j < count; ++j) {
      if (colors[j][0] <= cColorExtCutoff) {
        ColorGetRamped(ray->G, static_cast<int>(colors[j][0] - 0.1F), points[j],
            colors[j], -1);
      }
    }

    primitive.ramped = 0;
  }
}

/**
 * Scale the scene to fit within one meter and stand it centered on the
 * ground (y = 0), which is what AR viewers expect.
 */
void UsdWriteFitTransform(std::ostream& out, CRay* ray)
{
  RayComputeBox(ray);

  const float* lo = ray->min_box;
  const float* hi = ray->max_box;
  const float size = std::max({hi[0] - lo[0], hi[1] - lo[1], hi[2] - lo[2]});

  if (!(size > USD_EPSILON)) {
    return;
  }

  const float scale = 1.F / size;
  const float offset[3] = {
      -(lo[0] + hi[0]) * 0.5F, -lo[1], -(lo[2] + hi[2]) * 0.5F};

  out << "    float3 xformOp:scale = (" << scale << ", " << scale << ", "
      << scale << ")\n"
      << "    float3 xformOp:translate = ";
  UsdWriteVec3(out, offset);
  out << "\n"
      << "    uniform token[] xformOpOrder = [\"xformOp:scale\", "
         "\"xformOp:translate\"]\n";
}

} // namespace

/**
 * Generate an ASCII USD layer of the displayed geometry and append it to
 * `vla_ptr`.
 *
 * Coordinates are in Angstrom and in camera space, or in the original model
 * space with geometry_export_mode = 1, and the layer declares
 * metersPerUnit = 1e-10.
 *
 * Solids are tessellated, since importers like Blender's drop the material
 * of an analytic prim. With `ar` the layer is meant for AR viewers instead:
 * solids keep their compact analytic prims where possible, and the scene is
 * scaled to fit within one meter.
 */
void RayRenderUSDA(CRay* ray, char** vla_ptr, bool ar)
{
  const bool identity =
      SettingGetGlobal_i(ray->G, cSetting_geometry_export_mode) == 1;

  if (!RayExpandPrimitives(ray) || !RayTransformFirst(ray, 0, identity)) {
    return;
  }

  UsdResolveRampedColors(ray);

  ov_size count = 0;
  UsdVLAStreamBuf buffer(vla_ptr, &count);
  std::ostream out(&buffer);

  out << std::setprecision(6) << "#usda 1.0\n"
      << "(\n"
      << "    defaultPrim = \"PyMOLScene\"\n"
      << "    documentation = \"Exported from PyMOL\"\n"
      << "    metersPerUnit = " << (ar ? "1" : "1e-10") << "\n"
      << "    upAxis = \"Y\"\n"
      << ")\n"
      << "\n"
      << "def Xform \"PyMOLScene\"\n"
      << "{\n";

  if (ar) {
    UsdWriteFitTransform(out, ray);
  }

  UsdMaterials materials;
  UsdVectorMesh solids;
  int index = 0;
  const auto* basis = ray->Basis + 1;
  const UsdJointIndex joints(ray, basis);
  const int cylinder_segments =
      UsdSolidSegments(ray->G, cSetting_stick_quality);
  const int cone_segments = UsdSolidSegments(ray->G, cSetting_cone_quality);

  // Spheres are the bulk of a tessellated scene. The default sphere_quality
  // of 1 gives 16 segments.
  const int sphere_segments =
      std::clamp(8 * (SettingGetGlobal_i(ray->G, cSetting_sphere_quality) + 1),
          12, USD_SEGMENTS_MAX);

  for (int i = 0; i < ray->NPrimitive; ++i) {
    const auto& primitive = ray->Primitive[i];
    const auto* vertex = basis->Vertex + 3 * primitive.vert;
    const float opacity = 1.F - primitive.trans;

    const bool transparent = primitive.trans > USD_EPSILON;

    // An analytic UsdGeomCylinder or UsdGeomCone is always closed, so an end
    // which the ray tracer leaves open needs the mesh instead. Only while the
    // solid is opaque does a neighbour, or the whole sphere of a round cap,
    // hide the extra disc.
    const auto sealed = [&](const float* point, float radius, cCylCap cap) {
      return cap == cCylCapFlat ||
             (!transparent &&
                 (cap == cCylCapRound || joints.IsCovered(point, radius, i)));
    };

    switch (primitive.type) {
    case cPrimSphere:
      if (ar) {
        UsdWriteSphere(out, index, materials, vertex, primitive.r1,
            primitive.c1, primitive.trans);
      } else {
        UsdAddSphere(solids, sphere_segments, vertex, primitive.r1,
            primitive.c1, opacity);
      }
      break;
    case cPrimEllipsoid: {
      const auto* axes = basis->Normal + 3 * basis->Vert2Normal[primitive.vert];
      if (ar) {
        UsdWriteEllipsoid(out, index, materials, vertex, axes, primitive.n0,
            primitive.r1, primitive.c1, primitive.trans);
      } else {
        UsdAddEllipsoid(solids, sphere_segments, vertex, axes, primitive.n0,
            primitive.r1, primitive.c1, opacity);
      }
      break;
    }
    case cPrimCylinder:
    case cPrimSausage: {
      float end[3];
      UsdSolidEnd(basis, primitive, end);

      // A sausage is round capped no matter what the primitive says, and its
      // cap fields are never assigned
      const bool sausage = primitive.type == cPrimSausage;
      const cCylCap cap1 = sausage ? cCylCapRound : primitive.cap1;
      const cCylCap cap2 = sausage ? cCylCapRound : primitive.cap2;
      const bool round1 = cap1 == cCylCapRound;
      const bool round2 = cap2 == cCylCapRound;

      // An analytic prim takes a single color
      const bool one_color = UsdColorsEqual(primitive.c1, primitive.c2);

      // A capsule is the whole round capped solid as one closed surface, so
      // it needs no separate cap domes to bury inside it
      const bool capsule = ar && one_color && round1 && round2;

      if (capsule) {
        UsdWriteAnalyticSolid(out, index, materials, "Capsule", true, vertex,
            end, primitive.r1, primitive.c1, primitive.trans);
        break;
      }

      if (ar && one_color && sealed(vertex, primitive.r1, cap1) &&
          sealed(end, primitive.r1, cap2)) {
        UsdWriteAnalyticSolid(out, index, materials, "Cylinder", false, vertex,
            end, primitive.r1, primitive.c1, primitive.trans);
      } else {
        UsdAddCone(solids, cylinder_segments, vertex, end, primitive.r1,
            primitive.r1, primitive.c1, primitive.c2, cap1, cap2, opacity);
      }

      float axis[3];
      subtract3f(end, vertex, axis);
      const float height = length3f(axis);

      // Domes close the open barrel. A whole analytic sphere buries a
      // hemisphere in the barrel instead, which a transparent solid would
      // show as a darker cap the ray tracer does not draw.
      const bool dome = (!ar || transparent) && height > USD_EPSILON;

      if (dome) {
        scale3f(axis, 1.F / height, axis);
      }

      const auto add_cap = [&](const float* center, float sign,
                               const float* color) {
        if (dome) {
          UsdAddHemisphere(solids, cylinder_segments, center, axis, sign,
              primitive.r1, color, opacity);
        } else if (ar) {
          UsdWriteSphere(out, index, materials, center, primitive.r1, color,
              primitive.trans);
        } else {
          UsdAddSphere(
              solids, cylinder_segments, center, primitive.r1, color, opacity);
        }
      };

      if (round1) {
        add_cap(vertex, -1.F, primitive.c1);
      }
      if (round2) {
        add_cap(end, 1.F, primitive.c2);
      }
      break;
    }
    case cPrimCone: {
      float end[3];
      UsdSolidEnd(basis, primitive, end);

      // CRay::cone3fv is the only producer of cPrimCone and orders the radii
      assert(primitive.r1 >= primitive.r2);

      // The ray tracer only draws a flat cone cap, a round one is ignored
      const auto cap1 =
          primitive.cap1 == cCylCapFlat ? cCylCapFlat : cCylCapNone;
      const auto cap2 =
          primitive.cap2 == cCylCapFlat ? cCylCapFlat : cCylCapNone;

      // UsdGeomCone always tapers to a point and takes a single color
      if (ar && primitive.r2 <= USD_EPSILON &&
          sealed(vertex, primitive.r1, cap1) &&
          UsdColorsEqual(primitive.c1, primitive.c2)) {
        UsdWriteAnalyticSolid(out, index, materials, "Cone", false, vertex, end,
            primitive.r1, primitive.c1, primitive.trans);
      } else {
        UsdAddCone(solids, cone_segments, vertex, end, primitive.r1,
            primitive.r2, primitive.c1, primitive.c2, cap1, cap2, opacity);
      }
      break;
    }
    }
  }

  UsdWriteMesh(out, index, materials, "Solids", solids);

  const UsdSceneTriangles triangles(ray, basis);
  UsdWriteMesh(out, index, materials, "Mesh", triangles);

  materials.Write(out);
  out << "}\n";
  out.flush();
}
