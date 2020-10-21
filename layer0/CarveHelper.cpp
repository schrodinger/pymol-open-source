/**
 * @file
 * Isosurface carving
 *
 * (c) 2020 Schrodinger, Inc.
 */

#include "CarveHelper.h"
#include "Map.h"

/**
 * Positive cutoff value keeps triangles within the cutoff radius (a single
 * vertex outside the radius will discard the triangle). A negative cutoff value
 * inverts the criteria (discards triangle if all vertices are within -cutoff).
 *
 * @param cutoff Cutoff for inclusion (if positive) or exclusion (if negative)
 * @param vertices Flat vertex array
 * @param n_vertices Vertex count
 */
CarveHelper::CarveHelper(PyMOLGlobals* G, float cutoff, const float* vertices,
    std::size_t n_vertices)
    : m_cutoff(cutoff)
    , m_vertices(vertices)
{
  if (cutoff < 0) {
    m_cutoff = -cutoff;
    m_avoid_flag = true;
  }

  m_voxelmap.reset(MapNew(G, -m_cutoff, vertices, n_vertices, nullptr));
}

/**
 * True if the given vertex is within the (positive) cutoff radius of any atom
 */
bool CarveHelper::is_within(const float* v) const
{
  for (const auto j : MapEIter(*m_voxelmap, v)) {
    if (within3f(m_vertices + 3 * j, v, m_cutoff)) {
      return true;
    }
  }

  return false;
}

/**
 * True if the given triangle should be excluded
 */
bool CarveHelper::is_excluded(
    const float* v0, const float* v1, const float* v2) const
{
  bool within = is_within(v0) && is_within(v1) && is_within(v2);
  return m_avoid_flag ? within : !within;
}

/**
 * True if the given line should be excluded
 */
bool CarveHelper::is_excluded(const float* v0, const float* v1) const
{
  bool within = is_within(v0) && is_within(v1);
  return m_avoid_flag ? within : !within;
}

/**
 * True if the given vertex should be excluded
 */
bool CarveHelper::is_excluded(const float* v0) const
{
  bool within = is_within(v0);
  return m_avoid_flag ? within : !within;
}
