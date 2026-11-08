/**
 * @file
 * Isosurface with VTKm
 *
 * (c) 2020 Schrodinger, Inc.
 */

#include "ContourSurf.h"
#include "CarveHelper.h"
#include "Feedback.h"
#include "Isosurf.h"
#include "Tetsurf.h"
#include "Util.h"
#include "marching_cubes.h"

static constexpr size_t vertices_per_tri = 3;
static constexpr size_t floats_per_trivertex = 3 + 3; // xyz + normal
static constexpr size_t floats_per_tri =
    floats_per_trivertex * vertices_per_tri;

/**
 * Value to fill into ObjectSurfaceState::N
 */
inline unsigned get_num_per_strip(cIsosurfaceMode mode)
{
  switch (mode) {
  case cIsosurfaceMode::dots:
    return 1;
  case cIsosurfaceMode::lines:
    return 2;
  default:
    // triangles have 3 normals and 3 vertices
    return 6;
  }
}

/**
 * Resize and fill the "num" array, accoring to mode and number of floats.
 * @return number of primitives (e.g. triangles)
 */
static size_t fill_num_array(
    pymol::vla<int>& num, size_t float_count, cIsosurfaceMode mode)
{
  const auto num_per_strip = get_num_per_strip(mode);
  const auto n_tri = float_count / (3 * num_per_strip);
  num.resize(n_tri + 1);
  num[n_tri] = 0;
  std::fill_n(num.data(), n_tri, num_per_strip);
  return n_tri;
}

/**
 * Get triangle winding indices for the requested side
 */
static const int* get_winding_indices(cIsosurfaceSide side)
{
  static const int indices_winding_front[] = {0, 2, 1};
  static const int indices_winding_back[] = {0, 1, 2};
  return side == cIsosurfaceSide::front ? indices_winding_front
                                        : indices_winding_back;
}

/**
 * Flip the side for negative iso-levels
 */
static cIsosurfaceSide get_adjusted_side(float level, cIsosurfaceSide side)
{
  return level >= 0 ? side
                    : (side == cIsosurfaceSide::front) ? cIsosurfaceSide::back
                                                       : cIsosurfaceSide::front;
}


class PyMOLMcField : public mc::Field
{
  const Isofield* m_field = nullptr;
  int m_offset[3]{};
  int m_dim[3]{};

public:
  PyMOLMcField(const Isofield* field, const int* range)
      : m_field(field)
  {
    if (!range) {
      copy3(field->dimensions, m_dim);
    } else {
      copy3(range, m_offset);
      m_dim[0] = range[3] - range[0];
      m_dim[1] = range[4] - range[1];
      m_dim[2] = range[5] - range[2];
    }
  }

  size_t xDim() const override { return m_dim[0]; }
  size_t yDim() const override { return m_dim[1]; }
  size_t zDim() const override { return m_dim[2]; }

  float get(size_t x, size_t y, size_t z) const override
  {
    return m_field->data->get<float>(
        x + m_offset[0], y + m_offset[1], z + m_offset[2]);
  }

  mc::Point get_point(size_t x, size_t y, size_t z) const override
  {
    x += m_offset[0];
    y += m_offset[1];
    z += m_offset[2];
    return {
        m_field->points->get<float>(x, y, z, 0),
        m_field->points->get<float>(x, y, z, 1),
        m_field->points->get<float>(x, y, z, 2),
    };
  }
};

/**
 * Isosurface using code based on
 * https://github.com/ilastik/marching_cubes
 */
static int ContourSurfVolumeMcBasic(PyMOLGlobals* G, Isofield* field,
    float level,
    pymol::vla<int>& num,           //
    pymol::vla<float>& vert,        //
    const int* range,               //
    cIsosurfaceMode mode,           //
    const CarveHelper* carvehelper, //
    cIsosurfaceSide side)
{
  switch (mode) {
  case cIsosurfaceMode::triangles_grad_normals:
  case cIsosurfaceMode::triangles_tri_normals:
    break;
  default:
    PRINTFB(G, FB_Isosurface, FB_Warnings)
    " %s: Mode not implemented: %d\n", __func__, int(mode) ENDFB(G);
    return -1;
  }

  PyMOLMcField pmcfield(field, range);
  auto mesh = mc::march(
      pmcfield, level, mode == cIsosurfaceMode::triangles_grad_normals);

  if (mode == cIsosurfaceMode::triangles_tri_normals) {
    calculateNormals(mesh);
  }

  assert(mesh.normals);

  side = get_adjusted_side(level, side);
  const auto* indices_winding = get_winding_indices(side);
  int const normal_dir = int(side);
  assert(normal_dir == 1 || normal_dir == -1);

  size_t vert_size = 0;
  for (size_t i = 0; i < mesh.faceCount; ++i) {
    vert.check(vert_size + 6 * 3 - 1);

    for (size_t j = 0; j < 3; ++j) {
      const int jj = indices_winding[j];

      auto const& normal = mesh.normals[mesh.faces[i * 3 + jj]];
      vert[vert_size++] = normal.x * normal_dir;
      vert[vert_size++] = normal.y * normal_dir;
      vert[vert_size++] = normal.z * normal_dir;

      auto const& point = mesh.vertices[mesh.faces[i * 3 + jj]];
      vert[vert_size++] = point.x;
      vert[vert_size++] = point.y;
      vert[vert_size++] = point.z;
    }

    if (carvehelper && carvehelper->is_excluded(
                           &vert[vert_size - 3 - floats_per_trivertex * 0],
                           &vert[vert_size - 3 - floats_per_trivertex * 1],
                           &vert[vert_size - 3 - floats_per_trivertex * 2])) {
      vert_size -= floats_per_tri;
    }
  }

  vert.resize(vert_size);

  return fill_num_array(num, vert.size(), mode);
}

/**
 * Generate an isosurface with VTK-m. If VTK-m is not available, fall back to
 * the legacy Marching Tetrahedra implementation.
 *
 * @see TetsurfVolume
 */
int ContourSurfVolume(PyMOLGlobals* G, Isofield* field, float level,
    pymol::vla<int>& num,           //
    pymol::vla<float>& vert,        //
    const int* range,               //
    cIsosurfaceMode mode,           //
    const CarveHelper* carvehelper, //
    cIsosurfaceSide side)
{
  auto type = static_cast<cIsosurfaceAlgorithm>(
      SettingGet<int>(G, cSetting_isosurface_algorithm));
  int n_tri = 0;

  switch (type) {
  case cIsosurfaceAlgorithm::MARCHING_CUBES_VTKM: // value 0 kept for old sessions
  case cIsosurfaceAlgorithm::MARCHING_CUBES_BASIC:
    n_tri = ContourSurfVolumeMcBasic(
        G, field, level, num, vert, range, mode, carvehelper, side);
    if (n_tri >= 0)
      break;
  case cIsosurfaceAlgorithm::MARCHING_TETRAHEDRA:
    n_tri = TetsurfVolume(
        G, field, level, num, vert, range, mode, carvehelper, side);
    break;
  default:
    PRINTFB(G, FB_Isosurface, FB_Errors)
    " %s: Invalid surface_type: %d\n", __func__, int(type) ENDFB(G);
  }

  return n_tri;
}
