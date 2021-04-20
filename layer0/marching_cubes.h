/**
 * @file Marching cubes
 *
 * Based on https://github.com/ilastik/marching_cubes
 *
 * License: BSD-3-Clause
 *
 * Copyright 2018, The ilastik development team
 * Copyright 2020, Schrodinger, Inc.
 */

#ifndef MARCHING_CUBES_H
#define MARCHING_CUBES_H

#include <memory>

namespace mc
{

struct Point {
  float x, y, z;
  float operator[](size_t pos) const { return (&x)[pos]; }
  float& operator[](size_t pos) { return (&x)[pos]; }
};

class Field
{
public:
  virtual ~Field() = default;
  virtual size_t xDim() const = 0;
  virtual size_t yDim() const = 0;
  virtual size_t zDim() const = 0;
  virtual float get(size_t x, size_t y, size_t z) const = 0;
  virtual Point get_point(size_t x, size_t y, size_t z) const = 0;
  Point get_gradient(size_t x, size_t y, size_t z) const;
};

/**
 * Isosurface mesh with faces, vertices and normals
 */
struct Mesh {
  size_t vertexCount = 0;            //!< the number of vertices
  std::unique_ptr<Point[]> vertices; //!< vertex positions
  std::unique_ptr<Point[]> normals;  //!< normal direction of each vertex

  size_t faceCount = 0; //!< the number of faces
  std::unique_ptr<size_t[]>
      faces; //!< the faces given by 3 vertex indices (length = faceCount * 3)
};

Mesh march(const Field& volume, float isoLevel, bool gradient_normals = true);

void calculateNormals(Mesh& mesh);

} // namespace mc
#endif
