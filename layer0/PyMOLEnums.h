/**
 * Enums which are used across the API and implementations
 *
 * (c) 2020 Schrodinger, Inc.
 */

#pragma once

enum class cIsomeshMode {
  isomesh = 0,
  isodot = 1,
  gradient = 3,
};

enum class cIsosurfaceMode {
  dots = 0,
  lines = 1,
  triangles_tri_normals = 2,
  triangles_grad_normals = 3,
};

enum class cIsosurfaceSide {
  front = 1,
  back = -1,
};
