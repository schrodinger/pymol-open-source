
#pragma once

/*
A* -------------------------------------------------------------------
B* This file contains source code for the PyMOL computer program
C* Copyright (c) Schrodinger, LLC.
D* -------------------------------------------------------------------
E* It is unlawful to modify or remove this copyright notice.
F* -------------------------------------------------------------------
G* Please see the accompanying LICENSE file for further information.
H* -------------------------------------------------------------------
I* Additional authors of this source file include:
-*
-*
-*
Z* -------------------------------------------------------------------
*/

/* Map - a 3-D hash object for optimizing neighbor searches */

#include "Vector.h"

#include <vector>

struct PyMOLGlobals;

using MapFlag_t = int;

struct MapType {
  PyMOLGlobals* G;
  float Div;
  float recipDiv;
  Vector3i Dim;
  int D1D2;
  Vector3i iMin, iMax;
  std::vector<int> Head;
  std::vector<int> Link;
  std::vector<int> EHead;
  std::vector<int> EList;
  std::vector<int> EMask;
  int NVert;
  Vector3f Max, Min;
  MapType(PyMOLGlobals* G, float range, const float* vert, int nVert,
      const float* extent = nullptr, const MapFlag_t* flag = nullptr);
  int size() const;
};

struct MapCacheType {
  std::vector<int> Cache; // cached indices (0 = not cached, 1 = cached)
  std::vector<int> CacheLink; // linked list of cached indices
  int CacheStart = -1; // linked list head of cached indices

  // TODO: all caches should have valid map.
  // Current usages should be in a std::optional.
  MapCacheType() = default;
  MapCacheType(const MapType& I);

  /**
   * @param idx map index
   * @return true if the given map index is cached
   */
  bool cached(std::size_t idx) const;

  /**
   * Mark the given map index as cached
   * @param idx map index
   */
  void cache(std::size_t idx);
};

#define MapBorder 2

int MapSetupExpress(MapType* I);
int MapSetupExpressPerp(MapType* I, const float* vert, float front,
    int negative_start, const int* spanner);

#define MapFree(I) delete (I)

#define MapFirst(m, a, b, c)                                                   \
  (m->Head.data() + ((a) * m->D1D2) + ((b) * m->Dim[2]) + (c))

#define MapEStart(m, a, b, c)                                                  \
  (m->EHead.data() + ((a) * m->D1D2) + ((b) * m->Dim[2]) + (c))

#define MapNext(m, a) (m->Link[a])
void MapLocus(const MapType* map, const float* v, int* a, int* b, int* c);
int* MapLocusEStart(MapType* map, const float* v);

void MapCacheReset(MapCacheType& M);

float MapGetSeparation(PyMOLGlobals* G, float range, const float* mx,
    const float* mn, float* diagonal);
float MapGetDiv(MapType* I);

/* special routines for raytracing */

int MapInsideXY(MapType* I, const float* v, int* a, int* b, int* c);
int MapSetupExpressXY(MapType* I, int n_vert, int negative_start);

int MapSetupExpressXYVert(
    MapType* I, float* vert, int n_vert, int negative_start);

/**
 * Range iteration over points in proximity of a 3D query point
 */
class MapEIter
{
  const int* m_elist = nullptr;
  int m_i = 0;

public:
  MapEIter() = default;

  /**
   * @param map Map to query
   * @param v 3D query point
   * @param excl If true, exclude `v` if it's outside the grid
   */
  MapEIter(MapType& map, const float* v, bool excl = true);

  bool operator!=(MapEIter const& other) const { return m_i != other.m_i; }

  int operator*() const { return m_elist[m_i]; }

  MapEIter& operator++()
  {
    if (m_elist[++m_i] < 0) {
      m_i = 0;
    }
    return *this;
  }

  MapEIter begin() const { return *this; }
  MapEIter end() const { return {}; }
};

bool MapAnyWithin(
    MapType& map, const float* v_map, const float* v_query, float cutoff);
