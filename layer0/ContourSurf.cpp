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

#ifdef _PYMOL_VTKM
#ifdef _WIN32
// error C2039: '_copysign': is not a member of 'std'
// https://docs.microsoft.com/en-us/troubleshoot/cpp/c2653-c2039-error-reference-function
namespace std
{
using ::copysign;
}
#endif

#include <vtkm/Version.h>
#include <vtkm/cont/DataSetBuilderUniform.h>
#include <vtkm/filter/Contour.h>

#include <cassert>
#include <set>
#include <vector>

#if VTKM_VERSION_MINOR < 6 && VTKM_VERSION_PATCH > 0
// 1.5.1
#define ReadPortal GetPortalConstControl
#define make_ArrayHandleMove(v) make_ArrayHandle(v, vtkm::CopyFlag::Off)
#endif

/**
 * Make a copy of a field (with operator()) and expose an array-like data access
 * API (operator[], empty(), size()).
 */
template <typename T, size_t Dim = 3> class CopyFieldFunctor
{
  std::vector<T> m_data;

public:
  const bool empty() const { return m_data.empty(); }
  const size_t size() const { return m_data.size() / Dim; }
  const T* operator[](size_t pos) const { return m_data.data() + Dim * pos; }

  template <typename VecU, typename Storage>
  void operator()(const vtkm::cont::ArrayHandle<VecU, Storage>&)
  {
    assert(false);
  }

  template <typename U, typename Storage>
  void operator()(
      const vtkm::cont::ArrayHandle<vtkm::Vec<U, Dim>, Storage>& array)
  {
    using VecType = vtkm::VecTraits<typename vtkm::Vec<U, Dim>>;

    auto const& portal = array.ReadPortal();
    auto const numIndices = portal.GetNumberOfValues();

    m_data.reserve(numIndices * Dim);

    for (vtkm::Id index = 0; index < numIndices; ++index) {
      const auto& vec = portal.Get(index);
      assert(VecType::GetNumberOfComponents(vec) == Dim);
      for (vtkm::IdComponent c = 0; c < Dim; ++c) {
        m_data.emplace_back(VecType::GetComponent(vec, c));
      }
    }
  }
};

/**
 * Copies a surface to a flat array which corresponds to
 * ObjectSurfaceState::V (with constant ObjectSurfaceState::N)
 */
class CopySurfaceFunctor
{
  using CellSetType = vtkm::cont::CellSetSingleType<>;
  using NormalData = CopyFieldFunctor<float>;

  const CellSetType& m_cellset;
  const NormalData& m_normals;
  const cIsosurfaceSide m_side;
  const int* m_indices_winding;
  const CarveHelper* m_carvehelper;
  const cIsosurfaceMode m_mode;

  pymol::vla<float>& m_out; //!< corresponds to ObjectSurfaceState::V
  size_t m_out_size = 0;

  //! Set of already added lines or dots indices
  std::set<std::pair<size_t, size_t>> m_indices_done;

  //! Add vertex XYZ to the out array
  template <typename PointsPortal>
  void add_vertex(PointsPortal const& pointsPortal, size_t pointIndex)
  {
    using VecType = vtkm::VecTraits<typename PointsPortal::ValueType>;
    const auto& point = pointsPortal.Get(pointIndex);
    assert(VecType::GetNumberOfComponents(point) == 3);
    for (vtkm::IdComponent c = 0; c < 3; ++c) {
      *m_out.check(m_out_size++) = VecType::GetComponent(point, c);
    }
  }

  //! Add triangle corners to the out array if they have not beed added yet and
  //! if not carved.
  template <typename PointsPortal, typename ConnectivityPortal>
  void add_unique_dots_for_triangle(PointsPortal const& pointsPortal,
      ConnectivityPortal const& connPortal, size_t cellIndex)
  {
    for (size_t j = 0; j < vertices_per_tri; ++j) {
      size_t index = connPortal.Get(cellIndex * vertices_per_tri + j);
      if (m_indices_done.emplace(index, 0).second) {
        add_vertex(pointsPortal, index);

        if (m_carvehelper &&
            m_carvehelper->is_excluded(&m_out[m_out_size - 3])) {
          m_out_size -= 3;
        }
      }
    }
  }

  //! Add triangle edges to the out array if they have not beed added yet and if
  //! not carved.
  template <typename PointsPortal, typename ConnectivityPortal>
  void add_unique_lines_for_triangle(PointsPortal const& pointsPortal,
      ConnectivityPortal const& connPortal, size_t cellIndex)
  {
    for (int j = 0; j < vertices_per_tri; ++j) {
      size_t index1 = connPortal.Get(cellIndex * vertices_per_tri + j);
      size_t index2 = connPortal.Get(cellIndex * vertices_per_tri + //
                                     (j + 1) % vertices_per_tri);

      if (index2 > index1) {
        std::swap(index1, index2);
      }

      if (m_indices_done.emplace(index1, index2).second) {
        add_vertex(pointsPortal, index1);
        add_vertex(pointsPortal, index2);

        if (m_carvehelper && m_carvehelper->is_excluded( //
                                 &m_out[m_out_size - 3], //
                                 &m_out[m_out_size - 6])) {
          m_out_size -= 6;
        }
      }
    }
  }

  //! Add a triangle to the out array if not carved.
  template <typename PointsPortal, typename ConnectivityPortal>
  void add_triangle(PointsPortal const& pointsPortal,
      ConnectivityPortal const& connPortal, size_t cellIndex)
  {
    for (int j = 0; j < vertices_per_tri; ++j) {
      const int jj = m_indices_winding[j];
      const auto index = connPortal.Get(cellIndex * vertices_per_tri + jj);

      // add normal
      assert(index < m_normals.size());
      for (unsigned c = 0; c != 3; ++c) {
        assert(int(m_side) == 1 || int(m_side) == -1);
        *m_out.check(m_out_size++) = m_normals[index][c] * -int(m_side);
      }

      add_vertex(pointsPortal, index);
    }

    if (m_carvehelper &&
        m_carvehelper->is_excluded(
            &m_out[m_out_size - 3 - floats_per_trivertex * 0],
            &m_out[m_out_size - 3 - floats_per_trivertex * 1],
            &m_out[m_out_size - 3 - floats_per_trivertex * 2])) {
      m_out_size -= floats_per_tri;
    }
  }

public:
  CopySurfaceFunctor(const CellSetType& cellset, const NormalData& normals,
      pymol::vla<float>& out, cIsosurfaceSide side,
      const CarveHelper* carvehelper, cIsosurfaceMode mode)
      : m_cellset(cellset)
      , m_normals(normals)
      , m_out(out)
      , m_side(side)
      , m_carvehelper(carvehelper)
      , m_mode(mode)
  {
    m_indices_winding = get_winding_indices(side);
  }

  template <typename VecT, typename Storage>
  void operator()(const vtkm::cont::ArrayHandle<VecT, Storage>& array)
  {
    static_assert(VecT::NUM_COMPONENTS == 3, "");
    static_assert(
        std::is_floating_point<typename VecT::ComponentType>::value, "");

    auto const& pointsPortal = array.ReadPortal();
    auto const& connPortal =
        m_cellset
            .GetConnectivityArray(
                vtkm::TopologyElementTagCell(), vtkm::TopologyElementTagPoint())
            .ReadPortal();

    auto const nTri = m_cellset.GetNumberOfCells();
    auto const out_size_expected = nTri * floats_per_tri;
    m_out.reserve(out_size_expected);

    for (vtkm::Id i = 0; i < nTri; ++i) {
      assert(m_cellset.GetNumberOfPointsInCell(i) == vertices_per_tri);

      switch (m_mode) {
      case cIsosurfaceMode::triangles_grad_normals:
      case cIsosurfaceMode::triangles_tri_normals:
        add_triangle(pointsPortal, connPortal, i);
        break;
      case cIsosurfaceMode::lines:
        add_unique_lines_for_triangle(pointsPortal, connPortal, i);
        break;
      case cIsosurfaceMode::dots:
        add_unique_dots_for_triangle(pointsPortal, connPortal, i);
        break;
      }
    }

    m_out.resize(m_out_size);

    assert(m_out_size <= out_size_expected);
  }
};

/**
 * VTK-m implementation of ContourSurfVolume
 */
static int ContourSurfVolumeVtkm(PyMOLGlobals* G, Isofield* field, float level,
    pymol::vla<int>& num,           //
    pymol::vla<float>& vert,        //
    const int* range,               //
    cIsosurfaceMode const mode,     //
    const CarveHelper* carvehelper, //
    cIsosurfaceSide side)
{
  int range_store[6] = {0, 0, 0};
  if (!range) {
    copy3(field->dimensions, range_store + 3);
    range = range_store;
  }

  const CField& data = *field->data;     // CFieldTyped<float>(n_dims=3)
  const CField& points = *field->points; // CFieldTyped<float>(n_dim=4)

  std::string const pointFieldName = "pointdata";
  auto const pointDimensions =
      vtkm::Id3(range[3] - range[0], range[4] - range[1], range[5] - range[2]);
  auto const n_points =
      pointDimensions[0] * pointDimensions[1] * pointDimensions[2];

  // field as flat array in [z][y][x] order with shape pointDimensions
  std::vector<vtkm::Float32> pointdata;
  std::vector<vtkm::Vec3f> coorddata;
  pointdata.reserve(n_points);
  coorddata.reserve(n_points);

  for (size_t z = range[2]; z != range[5]; ++z) {
    for (size_t y = range[1]; y != range[4]; ++y) {
      for (size_t x = range[0]; x != range[3]; ++x) {
        pointdata.emplace_back(data.get<float>(x, y, z));
        coorddata.emplace_back(            //
            points.get<float>(x, y, z, 0), //
            points.get<float>(x, y, z, 1), //
            points.get<float>(x, y, z, 2));
      }
    }
  }

  auto dataSet = vtkm::cont::DataSetBuilderUniform().Create(pointDimensions);
  dataSet.GetCoordinateSystem().SetData(
      vtkm::cont::make_ArrayHandleMove(std::move(coorddata)));
  dataSet.AddField(vtkm::cont::make_Field(pointFieldName,
      vtkm::cont::Field::Association::POINTS, pointdata, vtkm::CopyFlag::Off));

  vtkm::filter::Contour filter;
  filter.SetIsoValue(level);
  filter.SetActiveField(pointFieldName);
  filter.SetGenerateNormals(mode == cIsosurfaceMode::triangles_grad_normals ||
                            mode == cIsosurfaceMode::triangles_tri_normals);
  filter.SetComputeFastNormalsForStructured(
      mode == cIsosurfaceMode::triangles_tri_normals);
  auto outputData = filter.Execute(dataSet);

  auto normalscopy = CopyFieldFunctor<float>();

  if (filter.GetGenerateNormals()) {
    outputData.GetPointField(filter.GetNormalArrayName())
        .GetData()
        .CastAndCall(normalscopy);
  }

  auto const& vertices = outputData.GetCoordinateSystem();
  auto const& cellsetsingletype =
      outputData.GetCellSet().Cast<vtkm::cont::CellSetSingleType<>>();

  side = get_adjusted_side(level, side);

  vtkm::cont::CastAndCall(
      vertices.GetData(), CopySurfaceFunctor{cellsetsingletype, normalscopy,
                              vert, side, carvehelper, mode});

  return fill_num_array(num, vert.size(), mode);
}

#endif

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
  case cIsosurfaceAlgorithm::MARCHING_CUBES_VTKM:
#ifdef _PYMOL_VTKM
    n_tri = ContourSurfVolumeVtkm(
        G, field, level, num, vert, range, mode, carvehelper, side);
    break;
#elif !defined(_WEBGL)
    PRINTFB(G, FB_Isosurface, FB_Warnings)
    " %s: VTKm not available, falling back to internal implementation\n",
        __func__ ENDFB(G);
#endif
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
