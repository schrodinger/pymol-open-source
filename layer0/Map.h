

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
#ifndef _H_Map
#define _H_Map


/* Map - a 3-D hash object for optimizing neighbor searches */

#include"Vector.h"
#include"PyMOLGlobals.h"

struct MapType {
  PyMOLGlobals *G;
  float Div;
  float recipDiv;
  Vector3i Dim;
  int D1D2;
  Vector3i iMin, iMax;
  int* Head = nullptr;
  int* Link = nullptr;
  int* EHead = nullptr;
  int* EList = nullptr;
  int* EMask = nullptr;
  int NVert;
  int NEElem = 0;
  Vector3f Max, Min;
  int group_id;
  int block_base;

  ~MapType();
};

typedef struct {
  PyMOLGlobals *G;
  int *Cache, *CacheLink, CacheStart;
  int block_base;

} MapCache;

#define MapBorder 2

MapType *MapNew(PyMOLGlobals * G, float range, const float *vert, int nVert,
    const float* extent = nullptr);
MapType *MapNewCached(PyMOLGlobals * G, float range, const float *vert, int nVert,
                      const float *extent, int group_id, int block_id);

using MapFlag_t = int;
MapType *MapNewFlagged(PyMOLGlobals * G, float range, const float *vert, int nVert,
                       const float *extent, const MapFlag_t *flag);
int MapSetupExpress(MapType * I);
int MapSetupExpressPerp(MapType * I, const float *vert, float front, int nVertHint,
			int negative_start, const int *spanner);

#define MapFree(I) delete(I)

#define MapFirst(m,a,b,c) (m->Head + ((a) * m->D1D2) + ((b)*m->Dim[2]) + (c))

#define MapEStart(m,a,b,c) (m->EHead + ((a) * m->D1D2) + ((b)*m->Dim[2]) + (c))

#define MapNext(m,a) (*(m->Link+(a)))
void MapLocus(const MapType * map, const float *v, int *a, int *b, int *c);
int *MapLocusEStart(MapType * map, const float *v);

#define MapCache(m,a) {(m)->Cache[a]=1;(m)->CacheLink[a]=(m)->CacheStart;(m)->CacheStart=a;}
#define MapCached(m,a) ((m)->Cache[a])

int MapCacheInit(MapCache * M, MapType * I, int group_id, int block_base);
void MapCacheReset(MapCache * M);
void MapCacheFree(MapCache * M, int group_id, int block_base);

float MapGetSeparation(PyMOLGlobals * G, float range, const float *mx, const float *mn,
                       float *diagonal);
float MapGetDiv(MapType * I);


/* special routines for raytracing */

int MapInsideXY(MapType * I, const float *v, int *a, int *b, int *c);
int MapSetupExpressXY(MapType * I, int n_vert, int negative_start);

int MapSetupExpressXYVert(MapType * I, float *vert, int n_vert, int negative_start);

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

#endif
