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
#include <cstring>
#include <iomanip>
#include <ostream>
#include <streambuf>
#include <vector>

namespace {

constexpr float USD_EPSILON = 1.e-6F;

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

void UsdWriteCylinder(std::ostream& out, int& index, const float* start,
    const float* end, float radius, const float* color, float transparency)
{
  float height;

  out << "\n"
      << "    def Cylinder \"Cylinder_" << index++ << "\" (\n"
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

void UsdWriteCone(std::ostream& out, int& index, const float* start,
    const float* end, float radius, const float* color, float transparency)
{
  float height;

  out << "\n"
      << "    def Cone \"Cone_" << index++ << "\" (\n"
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
 * @param center transformed ellipsoid center
 * @param axes three consecutive transformed unit axes
 * @param scales relative axis scales, one per axis, largest of which is 1
 * @param radius largest semi-axis length
 */
void UsdWriteEllipsoid(std::ostream& out, int& index, const float* center,
    const float* axes, const float* scales, float radius, const float* color,
    float transparency)
{
  float rows[9];

  for (int axis = 0; axis < 3; ++axis) {
    const auto* source = axes + 3 * axis;
    const float length = radius * scales[axis];
    for (int i = 0; i < 3; ++i) {
      rows[3 * axis + i] = source[i] * length;
    }
  }

  // The axes are eigenvectors with arbitrary signs and may form a left-handed
  // basis. Negating one axis leaves the ellipsoid unchanged.
  if (determinant33f(rows) < 0.0) {
    for (int i = 6; i < 9; ++i) {
      rows[i] = -rows[i];
    }
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

  for (int i = 0; i < ray->NPrimitive; ++i) {
    const auto& primitive = ray->Primitive[i];
    const auto* vertex = basis->Vertex + 3 * primitive.vert;

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
      const auto* direction =
          basis->Normal + 3 * basis->Vert2Normal[primitive.vert];
      const float end[3] = {vertex[0] + direction[0] * primitive.l1,
          vertex[1] + direction[1] * primitive.l1,
          vertex[2] + direction[2] * primitive.l1};

      if (UsdColorsEqual(primitive.c1, primitive.c2)) {
        UsdWriteCylinder(out, index, vertex, end, primitive.r1,
            primitive.c1, primitive.trans);
      } else {
        const float midpoint[3] = {(vertex[0] + end[0]) * 0.5F,
            (vertex[1] + end[1]) * 0.5F,
            (vertex[2] + end[2]) * 0.5F};
        UsdWriteCylinder(out, index, vertex, midpoint, primitive.r1,
            primitive.c1, primitive.trans);
        UsdWriteCylinder(out, index, midpoint, end, primitive.r1,
            primitive.c2, primitive.trans);
      }

      if (primitive.type == cPrimSausage ||
          primitive.cap1 == cCylCapRound) {
        UsdWriteSphere(out, index, vertex, primitive.r1, primitive.c1,
            primitive.trans);
      }
      if (primitive.type == cPrimSausage ||
          primitive.cap2 == cCylCapRound) {
        UsdWriteSphere(out, index, end, primitive.r1, primitive.c2,
            primitive.trans);
      }
      break;
    }
    case cPrimCone: {
      const auto* direction =
          basis->Normal + 3 * basis->Vert2Normal[primitive.vert];
      const float end[3] = {vertex[0] + direction[0] * primitive.l1,
          vertex[1] + direction[1] * primitive.l1,
          vertex[2] + direction[2] * primitive.l1};

      // CRay::cone3fv is the only producer of cPrimCone and orders the radii
      assert(primitive.r1 >= primitive.r2);

      if (primitive.r2 <= USD_EPSILON) {
        UsdWriteCone(out, index, vertex, end, primitive.r1, primitive.c1,
            primitive.trans);
      } else {
        // UsdGeomCone has a point at one end. Preserve a non-degenerate
        // frustum as a cylinder until a mesh representation is added.
        UsdWriteCylinder(out, index, vertex, end,
            (primitive.r1 + primitive.r2) * 0.5F, primitive.c1,
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
