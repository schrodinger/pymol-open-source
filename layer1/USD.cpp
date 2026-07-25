/*
 * PyMOL OpenUSD export
 *
 * The renderer writes an ASCII USD layer. Packaging and conversion to a
 * binary Crate layer are handled by modules/pymol/exporting.py.
 */

#include "Ray.h"

#include "Basis.h"
#include "Feedback.h"
#include "MemoryDebug.h"
#include "Setting.h"
#include "Vector.h"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <initializer_list>
#include <iomanip>
#include <ostream>
#include <streambuf>
#include <unordered_map>
#include <vector>

namespace {

constexpr float USD_EPSILON = 1.e-6F;

/// Extent of a collapsed ellipsoid semi-axis, relative to the largest one
constexpr float USD_DEGENERATE_AXIS_SCALE = 1.e-3F;

/**
 * Radial segments of a solid which cannot use an analytic UsdGeom prim.
 *
 * The ray tracer draws these solids analytically and has no tessellation of
 * its own, so PyMOL's quality settings can only raise the exported resolution
 * above the floor, never coarsen it below what a smooth render looks like.
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

void UsdWriteMaterialBinding(
    std::ostream& out, const float* color, float transparency)
{
  out << "        rel material:binding = </PyMOLScene/Material>\n"
      << "        color3f[] primvars:displayColor = [";
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

  UsdQuaternion result{
      1.F + dot, {-direction[1], direction[0], 0.F}};
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

void UsdWriteAnalyticTransform(std::ostream& out, const float* start,
    const float* end, float& height)
{
  float direction[3] = {
      end[0] - start[0], end[1] - start[1], end[2] - start[2]};
  height = std::sqrt(direction[0] * direction[0] +
      direction[1] * direction[1] + direction[2] * direction[2]);

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

void UsdWriteSphere(std::ostream& out, int& index, const float* center,
    float radius, const float* color, float transparency)
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
  UsdWriteMaterialBinding(out, color, transparency);
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
void UsdWriteAnalyticSolid(std::ostream& out, int& index, const char* schema,
    bool round_caps, const float* start, const float* end, float radius,
    const float* color, float transparency)
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
  UsdWriteMaterialBinding(out, color, transparency);
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
  /// Whether opacity varies per vertex rather than over the whole prim
  virtual bool VertexOpacity() const = 0;
  /// Opacity of the given vertex, or of the prim if not VertexOpacity()
  virtual float Opacity(std::size_t vertex) const = 0;
};

/// Mesh whose vertices and faces are buffered in vectors
struct UsdVectorMesh : UsdMeshSource {
  std::vector<float> points;
  std::vector<float> normals;
  std::vector<float> colors;
  std::vector<int> counts;
  std::vector<int> indices;
  float opacity = 1.F;

  std::size_t VertexCount() const override { return points.size() / 3; }
  std::size_t FaceCount() const override { return counts.size(); }
  int FaceSize(std::size_t face) const override { return counts[face]; }
  int Index(std::size_t corner) const override { return indices[corner]; }
  bool VertexOpacity() const override { return false; }
  float Opacity(std::size_t) const override { return opacity; }

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
  }

  void AddFace(std::initializer_list<int> face)
  {
    counts.push_back(static_cast<int>(face.size()));
    indices.insert(indices.end(), face.begin(), face.end());
  }
};

/// Write one UsdGeomMesh, shared by all exported mesh kinds
void UsdWriteMesh(std::ostream& out, int& index, const char* name,
    const UsdMeshSource& mesh)
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

  std::size_t corner_count = 0;
  for (std::size_t face = 0; face < face_count; ++face) {
    const int size = mesh.FaceSize(face);
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
      << "        color3f[] primvars:displayColor = [\n";
  for (std::size_t vertex = 0; vertex < vertex_count; ++vertex) {
    out << "            ";
    UsdWriteColor(out, mesh.Color(vertex));
    out << ",\n";
  }

  out << "        ] (\n"
      << "            interpolation = \"vertex\"\n"
      << "        )\n"
      << "        float[] primvars:displayOpacity = [";
  if (mesh.VertexOpacity()) {
    for (std::size_t vertex = 0; vertex < vertex_count; ++vertex) {
      out << (vertex ? ", " : "") << UsdClamp(mesh.Opacity(vertex));
    }
  } else {
    out << UsdClamp(mesh.Opacity(0));
  }

  out << "] (\n"
      << "            interpolation = \""
      << (mesh.VertexOpacity() ? "vertex" : "constant") << "\"\n"
      << "        )\n"
      << "        rel material:binding = </PyMOLScene/Material>\n"
      << "        uniform token subdivisionScheme = \"none\"\n"
      << "    }\n";
}

/**
 * Write a tessellated cone or truncated cone (frustum).
 *
 * UsdGeomCone always tapers to a point and takes a single color, so a cone
 * with a second radius or with two endpoint colors needs an explicit mesh.
 * Colors are interpolated along the axis, matching the ray tracer.
 *
 * @param segments radial segments, see UsdSolidSegments
 * @param start center of the r1 end
 * @param end center of the r2 end
 * @param cap1 cap of the r1 end, the ray tracer only draws flat caps
 * @param cap2 cap of the r2 end
 */
void UsdWriteConeMesh(std::ostream& out, int& index, int segments,
    const float* start, const float* end, float r1, float r2, const float* c1,
    const float* c2, cCylCap cap1, cCylCap cap2, float transparency)
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

  UsdVectorMesh mesh;
  mesh.opacity = 1.F - transparency;

  auto add_ring = [&](const float* center, float radius, const float* normal,
                      const float* color) {
    for (int k = 0; k < segments; ++k) {
      float point[3];
      scale3f(radial + 3 * k, radius, point);
      add3f(center, point, point);
      mesh.AddVertex(point, normal ? normal : slope + 3 * k, color);
    }
  };

  add_ring(start, r1, nullptr, c1);
  add_ring(end, r2, nullptr, c2);

  // Faces wind counter-clockwise as seen from outside, so that the right
  // handed face normals agree with the authored ones.
  for (int k = 0; k < segments; ++k) {
    const int next = (k + 1) % segments;
    if (pointed) {
      mesh.AddFace({k, next, segments + k});
    } else {
      mesh.AddFace({k, next, segments + next, segments + k});
    }
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

  UsdWriteMesh(out, index, "Cone", mesh);
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
    get_system1f3f(directions + 3 * first,
        directions + 3 * ((first + 1) % 3), directions + 3 * ((first + 2) % 3));
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
      static_cast<double>(lengths[1]) * static_cast<double>(lengths[2]);

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
void UsdWriteEllipsoid(std::ostream& out, int& index, const float* center,
    const float* axes, const float* scales, float radius, const float* color,
    float transparency)
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
  UsdWriteMaterialBinding(out, color, transparency);
  out << "    }\n";
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
  bool VertexOpacity() const override { return true; }

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
    return m_basis->Normal +
        3 * (m_basis->Vert2Normal[Primitive(vertex).vert] + 1 + Corner(vertex));
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

void UsdWriteMaterial(std::ostream& out)
{
  out << "    def Material \"Material\"\n"
      << "    {\n"
      << "        token outputs:surface.connect = "
         "</PyMOLScene/Material/PreviewSurface.outputs:surface>\n"
      << "\n"
      << "        def Shader \"DisplayColorReader\"\n"
      << "        {\n"
      << "            uniform token info:id = \"UsdPrimvarReader_float3\"\n"
      << "            string inputs:varname = \"displayColor\"\n"
      << "            float3 outputs:result\n"
      << "        }\n"
      << "\n"
      << "        def Shader \"DisplayOpacityReader\"\n"
      << "        {\n"
      << "            uniform token info:id = \"UsdPrimvarReader_float\"\n"
      << "            string inputs:varname = \"displayOpacity\"\n"
      << "            float outputs:result\n"
      << "        }\n"
      << "\n"
      << "        def Shader \"PreviewSurface\"\n"
      << "        {\n"
      << "            uniform token info:id = \"UsdPreviewSurface\"\n"
      << "            color3f inputs:diffuseColor.connect = "
         "</PyMOLScene/Material/DisplayColorReader.outputs:result>\n"
      << "            float inputs:opacity.connect = "
         "</PyMOLScene/Material/DisplayOpacityReader.outputs:result>\n"
      << "            float inputs:roughness = 0.35\n"
      << "            token outputs:surface\n"
      << "        }\n"
      << "    }\n";
}

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
 * PyMOL leaves the joint between abutting solids uncapped: the two halves of
 * a two colored stick meet with cCylCapNone, and a bond whose atom already
 * carries a round cap from a thicker bond is uncapped as well. Such a joint
 * is sealed by its neighbour, so it may keep the compact analytic prim, while
 * a genuinely exposed uncapped end needs an open mesh.
 *
 * The closed analytic prim buries an end disc inside the union, which only
 * goes unnoticed while the geometry is opaque. Transparent solids therefore
 * ignore this index and take the open mesh path, so that a viewing ray
 * crosses the same surfaces the ray tracer shows.
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
        static_cast<std::int64_t>(
            std::floor(point[0] / USD_JOINT_TOLERANCE)) + dx,
        static_cast<std::int64_t>(
            std::floor(point[1] / USD_JOINT_TOLERANCE)) + dy,
        static_cast<std::int64_t>(
            std::floor(point[2] / USD_JOINT_TOLERANCE)) + dz};

    return (cell[0] * 73856093) ^ (cell[1] * 19349663) ^ (cell[2] * 83492791);
  }

  void Add(const float* point, float radius, int primitive)
  {
    if (!(radius > USD_EPSILON)) {
      return;
    }

    m_cells.emplace(
        Key(point, 0, 0, 0), static_cast<int>(m_endpoints.size()));
    m_endpoints.push_back(
        {{point[0], point[1], point[2]}, radius, primitive});
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
 * Warn once about geometry colored by a color ramp.
 *
 * A ramp color is encoded as a large negative color value which the ray
 * tracer resolves per hit from the impact point. The exporter has no such
 * point, so UsdWriteColor clamps the encoding to black.
 */
void UsdWarnRamped(const CRay* ray)
{
  for (int i = 0; i < ray->NPrimitive; ++i) {
    if (ray->Primitive[i].ramped) {
      PRINTFB(ray->G, FB_Ray, FB_Warnings)
        " USD-Warning: ramp colors depend on the viewing ray and are not "
        "resolved, affected geometry is exported black.\n" ENDFB(ray->G);
      return;
    }
  }
}

} // namespace

void RayRenderUSDA(CRay* ray, char** vla_ptr)
{
  const bool identity =
      SettingGetGlobal_i(ray->G, cSetting_geometry_export_mode) == 1;

  if (!RayExpandPrimitives(ray) ||
      !RayTransformFirst(ray, 0, identity)) {
    return;
  }

  UsdWarnRamped(ray);

  ov_size count = 0;
  UsdVLAStreamBuf buffer(vla_ptr, &count);
  std::ostream out(&buffer);

  out << std::setprecision(9)
      << "#usda 1.0\n"
      << "(\n"
      << "    defaultPrim = \"PyMOLScene\"\n"
      << "    documentation = \"Exported from PyMOL\"\n"
      << "    metersPerUnit = 1e-10\n"
      << "    upAxis = \"Y\"\n"
      << ")\n"
      << "\n"
      << "def Xform \"PyMOLScene\"\n"
      << "{\n";

  UsdWriteMaterial(out);

  int index = 0;
  const auto* basis = ray->Basis + 1;
  const UsdJointIndex joints(ray, basis);
  const int cylinder_segments =
      UsdSolidSegments(ray->G, cSetting_stick_quality);
  const int cone_segments = UsdSolidSegments(ray->G, cSetting_cone_quality);

  for (int i = 0; i < ray->NPrimitive; ++i) {
    const auto& primitive = ray->Primitive[i];
    const auto* vertex = basis->Vertex + 3 * primitive.vert;

    // An analytic UsdGeomCylinder or UsdGeomCone is always closed, so an end
    // which the ray tracer leaves open needs the mesh writer instead. A
    // neighbour only hides the extra disc while the solid is opaque.
    const auto sealed = [&](const float* point, float radius, cCylCap cap) {
      return cap != cCylCapNone ||
          (!(primitive.trans > USD_EPSILON) &&
              joints.IsCovered(point, radius, i));
    };

    switch (primitive.type) {
    case cPrimSphere:
      UsdWriteSphere(out, index, vertex, primitive.r1, primitive.c1,
          primitive.trans);
      break;
    case cPrimEllipsoid:
      UsdWriteEllipsoid(out, index, vertex,
          basis->Normal + 3 * basis->Vert2Normal[primitive.vert],
          primitive.n0, primitive.r1, primitive.c1, primitive.trans);
      break;
    case cPrimCylinder:
    case cPrimSausage: {
      float end[3];
      UsdSolidEnd(basis, primitive, end);

      // A sausage is round capped no matter what the primitive says, and its
      // cap fields are never assigned
      const bool sausage = primitive.type == cPrimSausage;
      const bool round1 = sausage || primitive.cap1 == cCylCapRound;
      const bool round2 = sausage || primitive.cap2 == cCylCapRound;
      const cCylCap cap1 = sausage ? cCylCapNone : primitive.cap1;
      const cCylCap cap2 = sausage ? cCylCapNone : primitive.cap2;

      // The ray tracer blends the two endpoint colors along the axis, which
      // no single colored analytic prim can express
      const bool one_color = UsdColorsEqual(primitive.c1, primitive.c2);

      // A capsule is the whole round capped solid as one closed surface, so
      // it needs no separate cap spheres to bury inside it
      const bool capsule = one_color && round1 && round2;

      if (capsule) {
        UsdWriteAnalyticSolid(out, index, "Capsule", true, vertex, end,
            primitive.r1, primitive.c1, primitive.trans);
      } else if (one_color && (round1 || sealed(vertex, primitive.r1, cap1)) &&
          (round2 || sealed(end, primitive.r1, cap2))) {
        UsdWriteAnalyticSolid(out, index, "Cylinder", false, vertex, end,
            primitive.r1, primitive.c1, primitive.trans);
      } else {
        UsdWriteConeMesh(out, index, cylinder_segments, vertex, end,
            primitive.r1, primitive.r1, primitive.c1, primitive.c2, cap1, cap2,
            primitive.trans);
      }

      if (round1 && !capsule) {
        UsdWriteSphere(out, index, vertex, primitive.r1, primitive.c1,
            primitive.trans);
      }
      if (round2 && !capsule) {
        UsdWriteSphere(out, index, end, primitive.r1, primitive.c2,
            primitive.trans);
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

      if (primitive.r2 <= USD_EPSILON &&
          sealed(vertex, primitive.r1, cap1) &&
          UsdColorsEqual(primitive.c1, primitive.c2)) {
        UsdWriteAnalyticSolid(out, index, "Cone", false, vertex, end,
            primitive.r1, primitive.c1, primitive.trans);
      } else {
        UsdWriteConeMesh(out, index, cone_segments, vertex, end, primitive.r1,
            primitive.r2, primitive.c1, primitive.c2, cap1, cap2,
            primitive.trans);
      }
      break;
    }
    }
  }

  const UsdSceneTriangles triangles(ray, basis);
  UsdWriteMesh(out, index, "Mesh", triangles);

  out << "}\n";
  out.flush();
}
