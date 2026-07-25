/*
 * PyMOL OpenUSD export
 *
 * The renderer writes an ASCII USD layer. Packaging and conversion to a
 * binary Crate layer are handled by modules/pymol/exporting.py.
 */

#include "Ray.h"

#include "Basis.h"
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
#include <ostream>
#include <streambuf>
#include <unordered_map>
#include <vector>

namespace {

constexpr float USD_EPSILON = 1.e-6F;

/// Extent of a collapsed ellipsoid semi-axis, relative to the largest one
constexpr float USD_DEGENERATE_AXIS_SCALE = 1.e-3F;

/// Radial segments of a cone which cannot use the analytic UsdGeomCone
constexpr int USD_CONE_SEGMENTS = 24;

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
 * Write a Z-aligned UsdGeomCylinder or UsdGeomCone spanning start to end.
 * @param schema "Cylinder" or "Cone", also used as the prim name prefix
 */
void UsdWriteAnalyticSolid(std::ostream& out, int& index, const char* schema,
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
  out << "        uniform token axis = \"Z\"\n"
      << "        double height = " << height << "\n"
      << "        double radius = " << radius << "\n"
      << "        float3[] extent = [(" << -radius << ", " << -radius << ", "
      << -height * 0.5F << "), (" << radius << ", " << radius << ", "
      << height * 0.5F << ")]\n";
  UsdWriteMaterialBinding(out, color, transparency);
  out << "    }\n";
}

/**
 * Write a tessellated cone or truncated cone (frustum).
 *
 * UsdGeomCone always tapers to a point and takes a single color, so a cone
 * with a second radius or with two endpoint colors needs an explicit mesh.
 * Colors are interpolated along the axis, matching the ray tracer.
 *
 * @param start center of the r1 end
 * @param end center of the r2 end
 * @param cap1 cap of the r1 end, the ray tracer only draws flat caps
 * @param cap2 cap of the r2 end
 */
void UsdWriteConeMesh(std::ostream& out, int& index, const float* start,
    const float* end, float r1, float r2, const float* c1, const float* c2,
    cCylCap cap1, cCylCap cap2, float transparency)
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

  const bool pointed = !(r2 > USD_EPSILON);
  float radial[3 * USD_CONE_SEGMENTS];
  float slope[3 * USD_CONE_SEGMENTS];

  for (int k = 0; k < USD_CONE_SEGMENTS; ++k) {
    const auto angle =
        static_cast<float>(2.0 * PI * k / USD_CONE_SEGMENTS);
    const float cosine = std::cos(angle);
    const float sine = std::sin(angle);

    for (int i = 0; i < 3; ++i) {
      radial[3 * k + i] = side[i] * cosine + up[i] * sine;
      slope[3 * k + i] = radial[3 * k + i] * height + axis[i] * (r1 - r2);
    }

    normalize3f(slope + 3 * k);
  }

  std::vector<float> points;
  std::vector<float> normals;
  std::vector<float> colors;
  std::vector<int> counts;
  std::vector<int> indices;

  auto add_vertex = [&](const float* point, const float* normal,
                        const float* color) {
    points.insert(points.end(), point, point + 3);
    normals.insert(normals.end(), normal, normal + 3);
    colors.insert(colors.end(), color, color + 3);
  };

  auto add_face = [&](std::initializer_list<int> face) {
    counts.push_back(static_cast<int>(face.size()));
    indices.insert(indices.end(), face.begin(), face.end());
  };

  auto add_ring = [&](const float* center, float radius, const float* normal,
                      const float* color) {
    for (int k = 0; k < USD_CONE_SEGMENTS; ++k) {
      float point[3];
      scale3f(radial + 3 * k, radius, point);
      add3f(center, point, point);
      add_vertex(point, normal ? normal : slope + 3 * k, color);
    }
  };

  add_ring(start, r1, nullptr, c1);
  add_ring(end, r2, nullptr, c2);

  // Faces wind counter-clockwise as seen from outside, so that the right
  // handed face normals agree with the authored ones.
  for (int k = 0; k < USD_CONE_SEGMENTS; ++k) {
    const int next = (k + 1) % USD_CONE_SEGMENTS;
    if (pointed) {
      add_face({k, next, USD_CONE_SEGMENTS + k});
    } else {
      add_face({k, next, USD_CONE_SEGMENTS + next, USD_CONE_SEGMENTS + k});
    }
  }

  if (cap1 == cCylCapFlat) {
    const auto base = static_cast<int>(points.size() / 3);
    float normal[3];

    invert3f3f(axis, normal);
    add_vertex(start, normal, c1);
    add_ring(start, r1, normal, c1);

    for (int k = 0; k < USD_CONE_SEGMENTS; ++k) {
      const int next = (k + 1) % USD_CONE_SEGMENTS;
      add_face({base, base + 1 + next, base + 1 + k});
    }
  }

  if (cap2 == cCylCapFlat && !pointed) {
    const auto base = static_cast<int>(points.size() / 3);

    add_vertex(end, axis, c2);
    add_ring(end, r2, axis, c2);

    for (int k = 0; k < USD_CONE_SEGMENTS; ++k) {
      const int next = (k + 1) % USD_CONE_SEGMENTS;
      add_face({base, base + 1 + k, base + 1 + next});
    }
  }

  float extent_min[3];
  float extent_max[3];

  copy3f(points.data(), extent_min);
  copy3f(points.data(), extent_max);

  for (std::size_t i = 3; i < points.size(); i += 3) {
    for (int axis_index = 0; axis_index < 3; ++axis_index) {
      extent_min[axis_index] =
          std::min(extent_min[axis_index], points[i + axis_index]);
      extent_max[axis_index] =
          std::max(extent_max[axis_index], points[i + axis_index]);
    }
  }

  out << "\n"
      << "    def Mesh \"Cone_" << index++ << "\" (\n"
      << "        prepend apiSchemas = [\"MaterialBindingAPI\"]\n"
      << "    )\n"
      << "    {\n"
      << "        float3[] extent = [";
  UsdWriteVec3(out, extent_min);
  out << ", ";
  UsdWriteVec3(out, extent_max);
  out << "]\n"
      << "        int[] faceVertexCounts = [";
  for (std::size_t i = 0; i < counts.size(); ++i) {
    out << (i ? ", " : "") << counts[i];
  }
  out << "]\n"
      << "        int[] faceVertexIndices = [";
  for (std::size_t i = 0; i < indices.size(); ++i) {
    out << (i ? ", " : "") << indices[i];
  }
  out << "]\n"
      << "        point3f[] points = [\n";
  for (std::size_t i = 0; i < points.size(); i += 3) {
    out << "            ";
    UsdWriteVec3(out, points.data() + i);
    out << ",\n";
  }
  out << "        ]\n"
      << "        normal3f[] normals = [\n";
  for (std::size_t i = 0; i < normals.size(); i += 3) {
    out << "            ";
    UsdWriteVec3(out, normals.data() + i);
    out << ",\n";
  }
  out << "        ] (\n"
      << "            interpolation = \"vertex\"\n"
      << "        )\n"
      << "        color3f[] primvars:displayColor = [\n";
  for (std::size_t i = 0; i < colors.size(); i += 3) {
    out << "            ";
    UsdWriteColor(out, colors.data() + i);
    out << ",\n";
  }
  out << "        ] (\n"
      << "            interpolation = \"vertex\"\n"
      << "        )\n"
      << "        float[] primvars:displayOpacity = ["
      << UsdClamp(1.F - transparency) << "] (\n"
      << "            interpolation = \"constant\"\n"
      << "        )\n"
      << "        rel material:binding = </PyMOLScene/Material>\n"
      << "        uniform token subdivisionScheme = \"none\"\n"
      << "    }\n";
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

/// One exported triangle, with its winding already resolved
struct UsdTriangle {
  const CPrimitive* primitive;
  const float* vertex;
  const float* normal;
  int order[3];
};

std::vector<UsdTriangle> UsdCollectTriangles(
    const CRay* ray, const CBasis* basis)
{
  std::vector<UsdTriangle> triangles;

  for (int i = 0; i < ray->NPrimitive; ++i) {
    auto& primitive = ray->Primitive[i];
    if (primitive.type != cPrimTriangle) {
      continue;
    }

    const bool reverse =
        TriangleReverse(const_cast<CPrimitive*>(&primitive));
    triangles.push_back({&primitive, basis->Vertex + 3 * primitive.vert,
        basis->Normal + 3 * basis->Vert2Normal[primitive.vert] + 3,
        {0, reverse ? 2 : 1, reverse ? 1 : 2}});
  }

  return triangles;
}

void UsdWriteTriangleMesh(
    std::ostream& out, int& index, const CRay* ray, const CBasis* basis)
{
  const auto triangles = UsdCollectTriangles(ray, basis);
  const auto triangle_count = triangles.size();
  if (!triangle_count) {
    return;
  }

  float extent_min[3] = {0.F, 0.F, 0.F};
  float extent_max[3] = {0.F, 0.F, 0.F};
  bool have_extent = false;

  for (const auto& triangle : triangles) {
    for (int j = 0; j < 9; j += 3) {
      for (int axis = 0; axis < 3; ++axis) {
        if (!have_extent) {
          extent_min[axis] = extent_max[axis] = triangle.vertex[j + axis];
        } else {
          extent_min[axis] =
              std::min(extent_min[axis], triangle.vertex[j + axis]);
          extent_max[axis] =
              std::max(extent_max[axis], triangle.vertex[j + axis]);
        }
      }
      have_extent = true;
    }
  }

  out << "\n"
      << "    def Mesh \"Mesh_" << index++ << "\" (\n"
      << "        prepend apiSchemas = [\"MaterialBindingAPI\"]\n"
      << "    )\n"
      << "    {\n"
      << "        float3[] extent = [";
  UsdWriteVec3(out, extent_min);
  out << ", ";
  UsdWriteVec3(out, extent_max);
  out << "]\n"
      << "        int[] faceVertexCounts = [";
  for (std::size_t i = 0; i < triangle_count; ++i) {
    out << (i ? ", 3" : "3");
  }
  out << "]\n"
      << "        int[] faceVertexIndices = [";
  for (std::size_t i = 0; i < triangle_count * 3; ++i) {
    out << (i ? ", " : "") << i;
  }
  out << "]\n"
      << "        point3f[] points = [\n";

  for (const auto& triangle : triangles) {
    for (const int vert : triangle.order) {
      out << "            ";
      UsdWriteVec3(out, triangle.vertex + 3 * vert);
      out << ",\n";
    }
  }

  out << "        ]\n"
      << "        normal3f[] normals = [\n";
  for (const auto& triangle : triangles) {
    for (const int vert : triangle.order) {
      out << "            ";
      UsdWriteVec3(out, triangle.normal + 3 * vert);
      out << ",\n";
    }
  }
  out << "        ] (\n"
      << "            interpolation = \"vertex\"\n"
      << "        )\n"
      << "        color3f[] primvars:displayColor = [\n";

  for (const auto& triangle : triangles) {
    const float* colors[3] = {triangle.primitive->c1, triangle.primitive->c2,
        triangle.primitive->c3};
    for (const int vert : triangle.order) {
      out << "            ";
      UsdWriteColor(out, colors[vert]);
      out << ",\n";
    }
  }
  out << "        ] (\n"
      << "            interpolation = \"vertex\"\n"
      << "        )\n"
      << "        float[] primvars:displayOpacity = [";

  bool first_opacity = true;
  for (const auto& triangle : triangles) {
    for (const int vert : triangle.order) {
      out << (first_opacity ? "" : ", ")
          << UsdClamp(1.F - triangle.primitive->tr[vert]);
      first_opacity = false;
    }
  }
  out << "] (\n"
      << "            interpolation = \"vertex\"\n"
      << "        )\n"
      << "        rel material:binding = </PyMOLScene/Material>\n"
      << "        uniform token subdivisionScheme = \"none\"\n"
      << "    }\n";
}

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

} // namespace

void RayRenderUSDA(CRay* ray, char** vla_ptr)
{
  const bool identity =
      SettingGetGlobal_i(ray->G, cSetting_geometry_export_mode) == 1;

  if (!RayExpandPrimitives(ray) ||
      !RayTransformFirst(ray, 0, identity)) {
    return;
  }

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

  for (int i = 0; i < ray->NPrimitive; ++i) {
    const auto& primitive = ray->Primitive[i];
    const auto* vertex = basis->Vertex + 3 * primitive.vert;

    // An analytic UsdGeomCylinder or UsdGeomCone is always closed, so an end
    // which the ray tracer leaves open needs the mesh writer instead
    const auto sealed = [&](const float* point, float radius, cCylCap cap) {
      return cap != cCylCapNone || joints.IsCovered(point, radius, i);
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

      // A sausage is round capped no matter what the primitive says
      const bool round1 =
          primitive.type == cPrimSausage || primitive.cap1 == cCylCapRound;
      const bool round2 =
          primitive.type == cPrimSausage || primitive.cap2 == cCylCapRound;

      if ((round1 || sealed(vertex, primitive.r1, primitive.cap1)) &&
          (round2 || sealed(end, primitive.r1, primitive.cap2))) {
        if (UsdColorsEqual(primitive.c1, primitive.c2)) {
          UsdWriteAnalyticSolid(out, index, "Cylinder", vertex, end,
              primitive.r1, primitive.c1, primitive.trans);
        } else {
          const float midpoint[3] = {(vertex[0] + end[0]) * 0.5F,
              (vertex[1] + end[1]) * 0.5F,
              (vertex[2] + end[2]) * 0.5F};
          UsdWriteAnalyticSolid(out, index, "Cylinder", vertex, midpoint,
              primitive.r1, primitive.c1, primitive.trans);
          UsdWriteAnalyticSolid(out, index, "Cylinder", midpoint, end,
              primitive.r1, primitive.c2, primitive.trans);
        }
      } else {
        UsdWriteConeMesh(out, index, vertex, end, primitive.r1, primitive.r1,
            primitive.c1, primitive.c2, primitive.cap1, primitive.cap2,
            primitive.trans);
      }

      if (round1) {
        UsdWriteSphere(out, index, vertex, primitive.r1, primitive.c1,
            primitive.trans);
      }
      if (round2) {
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

      if (primitive.r2 <= USD_EPSILON &&
          sealed(vertex, primitive.r1, cap1) &&
          UsdColorsEqual(primitive.c1, primitive.c2)) {
        UsdWriteAnalyticSolid(out, index, "Cone", vertex, end, primitive.r1,
            primitive.c1, primitive.trans);
      } else {
        UsdWriteConeMesh(out, index, vertex, end, primitive.r1, primitive.r2,
            primitive.c1, primitive.c2, primitive.cap1, primitive.cap2,
            primitive.trans);
      }
      break;
    }
    }
  }

  UsdWriteTriangleMesh(out, index, ray, basis);
  out << "}\n";
  out.flush();
}
