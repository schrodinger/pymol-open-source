/**
 * @file
 * Isosurface carving
 *
 * (c) 2020 Schrodinger, Inc.
 */

#pragma once

#include <cstddef>
#include <memory>

struct PyMOLGlobals;
struct MapType;

/**
 * Data structure for "carving" isosurfaces around a set of vertices.
 */
class CarveHelper
{
  std::unique_ptr<MapType> m_voxelmap;
  const float* m_vertices;
  float m_cutoff;
  bool m_avoid_flag = false;

  bool is_within(const float* v) const;

public:
  CarveHelper(PyMOLGlobals* G, float cutoff, const float* vertices,
      std::size_t n_vertices);

  bool is_excluded(const float* v0) const;
  bool is_excluded(const float* v0, const float* v1) const;
  bool is_excluded(const float* v0, const float* v1, const float* v2) const;
};
