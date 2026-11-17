

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
-*   Larry Coopet (various optimizations)
-*
-*
Z* ------------------------------------------------------------------- */
#include "os_predef.h"
#include "os_python.h"
#include "os_std.h"

#include "Base.h"
#include "Err.h"
#include "Feedback.h"
#include "Map.h"
#include "MemoryDebug.h"
#include "Setting.h"

#include "pymol/algorithm.h"

float MapGetDiv(MapType* I)
{
  return I->Div;
}

MapCacheType::MapCacheType(const MapType& map)
{
  Cache.assign(map.NVert, 0);
  CacheLink.resize(map.NVert);
}

bool MapCacheType::cached(std::size_t idx) const
{
  return Cache[idx] != 0;
}

void MapCacheType::cache(std::size_t idx)
{
  Cache[idx] = 1;
  CacheLink[idx] = CacheStart;
  CacheStart = idx;
}

void MapCacheReset(MapCacheType& M)
{
  int i = M.CacheStart;
  auto* cachep = M.Cache.data();
  auto* clinkp = M.CacheLink.data();
  int i1 = 0, i2 = 0, i3 = 0, i4 = 0, ii;
  while (i >= 0) { /* believe it or not, unrolling gives us almost 10%!!! */
    ii = clinkp[i];
    i1 = i;
    i = ii;
    if (ii >= 0) {
      ii = clinkp[ii];
      i2 = i;
      i = ii;
    }
    cachep[i1] = 0; /* this doesn't look safe, but it is. i1-i4 are always valid
                       indices */
    if (ii >= 0) {
      ii = clinkp[ii];
      i3 = i;
      i = ii;
    }
    cachep[i2] = 0; /* this doesn't look safe, but it is. i1-i4 are always valid
                       indices */
    if (ii >= 0) {
      ii = clinkp[ii];
      i4 = i;
      i = ii;
    }
    cachep[i3] = 0; /* this doesn't look safe, but it is. i1-i4 are always valid
                       indices */
    cachep[i4] = 0; /* this doesn't look safe, but it is. i1-i4 are always valid
                       indices */
  }
  M.CacheStart = -1;
}

#define MapSafety 0.01F

int MapInsideXY(MapType* I, const float* v, int* a, int* b, int* c)
{ /* special version for ray-tracing */
  int atmp, btmp, ctmp;
  const float iDiv = I->recipDiv;

  atmp = (int) ((v[0] - I->Min[0]) * iDiv) + MapBorder;
  btmp = (int) ((v[1] - I->Min[1]) * iDiv) + MapBorder;
  ctmp = (int) ((v[2] - I->Min[2]) * iDiv) + MapBorder + 1;

  if (atmp < I->iMin[0]) {
    if ((I->iMin[0] - atmp) > 1)
      return (false);
    else
      atmp = I->iMin[0];
  } else if (atmp > I->iMax[0]) {
    if ((atmp - I->iMax[0]) > 1)
      return (false);
    else
      atmp = I->iMax[0];
  }

  if (btmp < I->iMin[1]) {
    if ((I->iMin[1] - btmp) > 1)
      return (false);
    else
      btmp = I->iMin[1];
  } else if (btmp > I->iMax[1]) {
    if ((btmp - I->iMax[1]) > 1)
      return (false);
    else
      btmp = I->iMax[1];
  }

  if (!*(I->EMask.data() + I->Dim[1] * atmp + btmp))
    return (false);

  if (ctmp < I->iMin[2])
    ctmp = I->iMin[2];
  else if (ctmp > I->iMax[2])
    ctmp = I->iMax[2];

  *a = atmp;
  *b = btmp;
  *c = ctmp;

  return (true);
}

#define ELIST_GROW_FACTOR 3

int MapSetupExpressXY(MapType* I, int n_vert, int negative_start)
{ /* setup a list of XY neighbors for each square */
  PyMOLGlobals* G = I->G;
  int a, b, c, flag;
  int d, e, i;
  unsigned int mapSize;
  int dim2;
  int n_alloc = n_vert * 15; /* emprical est. */
  int ok = true;

  PRINTFD(G, FB_Map)
  " MapSetupExpressXY-Debug: entered.\n" ENDFD;

  mapSize = I->Dim[0] * I->Dim[1] * I->Dim[2];
  I->EHead.assign(mapSize, 0);
  if (ok) {
    I->EList.reserve(n_alloc);
  }
  if (ok)
    I->EMask.assign(I->Dim[0] * I->Dim[1], 0);

  I->EList.push_back(0);
  int n = 1;
  dim2 = I->Dim[2];

  for (a = I->iMin[0]; ok && a <= I->iMax[0]; a++) {
    for (b = I->iMin[1]; ok && b <= I->iMax[1]; b++) {
      for (c = I->iMin[2]; ok && c <= I->iMax[2];
           c++) { /* a better alternative exists... */
        int* iPtr1 =
            (I->Head.data() + ((a - 1) * I->D1D2) + ((b - 1) * dim2) + c);

        int st = n;
        flag = false;

        for (d = a - 1; d <= a + 1; d++) {
          /*int    *iPtr2  = (I->Head + (d * I->D1D2) + ((b-1)*dim2) + c); */
          int* iPtr2 = iPtr1;

          for (e = b - 1; e <= b + 1; e++) {
            /*i      = *MapFirst(I,d,e,c); */
            i = *iPtr2;
            if (i >= 0) {
              flag = true;
              while (i >= 0) {
                I->EList.push_back(i);
                n++;
                i = MapNext(I, i);
              }
            }

            iPtr2 += dim2;
          }

          iPtr1 += I->D1D2;
        }

        if (ok && flag) {
          *(I->EMask.data() + I->Dim[1] * a + b) = true;
          *(MapEStart(I, a, b, c)) = negative_start ? -st : st;
          I->EList.push_back(-1);
          n++;
        }
      }
    }
  }

  PRINTFB(G, FB_Map, FB_Blather)
  " MapSetupExpressXY: %d rows in express table\n", n ENDFB(G);

  if (ok) {
    I->EList.shrink_to_fit();
  }
  PRINTFD(G, FB_Map)
  " MapSetupExpressXY-Debug: leaving...\n" ENDFD;
  return ok;
}

int MapSetupExpressXYVert(
    MapType* I, float* vert, int n_vert, int negative_start)
{ /* setup a list of XY neighbors for each square */
  PyMOLGlobals* G = I->G;
  int h, n, a, b, c;
  int j, k, dim2;
  float* v;
  int *eBase, *hBase;
  int n_alloc = n_vert * 15; /* emprical est. */
  int ok = true;

  PRINTFD(G, FB_Map)
  " MapSetupExpressXYVert-Debug: entered n_vert = %d negative_start = %d\n",
      n_vert, negative_start ENDFD;

  /*mapSize        = I->Dim[0]*I->Dim[1]*I->Dim[2]; */
  I->EHead.assign(I->Dim[0] * I->Dim[1] * I->Dim[2], 0);
  if (ok)
    I->EMask.assign(I->Dim[0] * I->Dim[1], 0);
  if (ok) {
    I->EList.clear();
    I->EList.reserve(n_alloc);
  }

  I->EList.push_back(0);
  n = 1;
  v = vert;
  dim2 = I->Dim[2];

  for (h = 0; h < n_vert; h++) {
    MapLocus(I, v, &j, &k, &c);

    eBase = I->EHead.data() + ((j - 1) * I->D1D2) + ((k - 1) * dim2) + c;
    hBase = I->Head.data() + (((j - 1) - 1) * I->D1D2);

    for (a = j - 1; ok && a <= j + 1; a++) {
      int* ePtr1 = eBase;

      for (b = k - 1; ok && b <= k + 1; b++) {
        if (*ePtr1 == 0) { /* if we haven't yet expanded this voxel... */
          int st = n;
          int flag = false;
          int* hPtr1 = hBase + dim2 * (b - 1) + (c - 1);

          int d, e, f;

          for (d = a - 1; ok && d <= a + 1; d++) {
            int* hPtr2 = hPtr1;
            for (e = b - 1; ok && e <= b + 1; e++) {
              int* hPtr3 = hPtr2;
              for (f = c - 1; ok && f <= c + 1; f++) {
                /*                register int i = *MapFirst(I,d,e,f); */
                int i = *hPtr3;

                if (i > -1) {
                  flag = true;
                  while (ok && i > -1) {
                    I->EList.push_back(i);
                    n++;
                    i = MapNext(I, i);
                  }
                }
                hPtr3++;
              }
              hPtr2 += dim2;
            }
            hPtr1 += I->D1D2;
          }

          if (flag) {
            *(I->EMask.data() + I->Dim[1] * a + b) = true;
            *(MapEStart(I, a, b, c)) = negative_start ? -st : st;
            I->EList.push_back(-1);
            n++;
          }
        }

        ePtr1 += dim2;
      }

      eBase += I->D1D2;
      hBase += I->D1D2;
    }

    v += 3;
  }

  PRINTFB(G, FB_Map, FB_Blather)
  " MapSetupExpressXYVert: %d rows in express table\n", n ENDFB(G);

  if (ok) {
    I->EList.shrink_to_fit();
  }
  PRINTFD(G, FB_Map)
  " MapSetupExpressXYVert-Debug: leaving...\n" ENDFD;
  return ok;
}

int MapSetupExpressPerp(MapType* I, const float* vert, float front,
    int nVertHint, int negative_start, const int* spanner)
{
  PyMOLGlobals* G = I->G;
  int a, b, c, i;

  unsigned int mapSize;
  int st;
  int n_alloc = nVertHint * 15; /* emprical est. */
  int ok = true;

  int iMin0 = I->iMin[0];
  int iMin1 = I->iMin[1];
  int iMax0 = I->iMax[0];
  int iMax1 = I->iMax[1];
  float iDiv = I->recipDiv;
  float min0 = I->Min[0] * iDiv;
  float min1 = I->Min[1] * iDiv;
  float base0, base1;
  float perp_factor, premult;
  int *emask, dim1, *link, *ptr1, *ptr2;

  PRINTFD(G, FB_Map)
  " MapSetupExpress-Debug: entered.\n" ENDFD;

  mapSize = I->Dim[0] * I->Dim[1] * I->Dim[2];
  I->EHead.assign(mapSize, 0);
  if (ok) {
    I->EList.reserve(n_alloc);
  }
  if (ok)
    I->EMask.assign(I->Dim[0] * I->Dim[1], 0);

  emask = I->EMask.data();
  dim1 = I->Dim[1];
  link = I->Link.data();
  premult = -front * iDiv;

  I->EList.push_back(0);
  int n = 1;

  for (a = (iMin0 - 1); ok && a <= (iMax0 + 1); a++)
    for (b = (iMin1 - 1); ok && b <= (iMax1 + 1); b++)
      for (c = (I->iMin[2] - 1); ok && c <= (I->iMax[2] + 1); c++) {

        int d, e, f;

        /* compute a "shadow" mask for all vertices */

        i = *MapFirst(I, a, b, c);
        while (i >= 0) {
          const float* v0 = vert + 3 * i;
          perp_factor = premult / v0[2];
          base0 = v0[0] * perp_factor;
          base1 = v0[1] * perp_factor;

          d = (int) (base0 - min0);
          e = (int) (base1 - min1);

          d += MapBorder;
          e += MapBorder;

          if (d < iMin0) {
            d = iMin0;
          } else if (d > iMax0) {
            d = iMax0;
          }

          if (e < iMin1) {
            e = iMin1;
          } else if (e > iMax1) {
            e = iMax1;
          }
          i = link[i];
          ptr2 = (ptr1 = emask + dim1 * (d - 1) + (e - 1));
          *(ptr2++) = true;
          *(ptr2++) = true;
          *(ptr2++) = true;
          ptr2 = (ptr1 += dim1);
          *(ptr2++) = true;
          *(ptr2++) = true;
          *(ptr2++) = true;
          ptr2 = (ptr1 += dim1);
          *(ptr2++) = true;
          *(ptr2++) = true;
          *(ptr2++) = true;
        }

        {
          const int am1 = a - 1, ap1 = a + 1, bm1 = b - 1, bp1 = b + 1,
                    cm1 = c - 1, cp1 = c + 1;
          const int dim2 = I->Dim[2];
          int flag = false;
          int* hPtr1 = I->Head.data() + ((am1) *I->D1D2) + ((bm1) *dim2) + cm1;
          st = n;
          for (d = am1; ok && d <= ap1; d++) {
            int* hPtr2 = hPtr1;
            for (e = bm1; ok && e <= bp1; e++) {
              int* hPtr3 = hPtr2;
              for (f = cm1; ok && f <= cp1; f++) {
                i = *(hPtr3++);
                /*                i=*MapFirst(I,d,e,f); */
                if (i >= 0) {
                  flag = true;
                  while (ok && i >= 0) {
                    if ((!spanner) || (f == c) || spanner[i]) {
                      /* for non-voxel-spanners, only spread in the XY plane
                       * (memory use ~ 9X instead of 27X -- a big difference!)
                       */
                      I->EList.push_back(i);
                      n++;
                    }
                    i = link[i];
                  }
                }
              }
              hPtr2 += dim2;
            }
            hPtr1 += I->D1D2;
          }
          if (ok && flag) {
            *(MapEStart(I, a, b, c)) = negative_start ? -st : st;
            I->EList.push_back(-1);
            n++;
          }
        }
      }
  PRINTFB(G, FB_Map, FB_Blather)
  " MapSetupExpressPerp: %d rows in express table \n", n ENDFB(G);
  if (ok) {
    I->EList.shrink_to_fit();
  }
  PRINTFD(G, FB_Map)
  " MapSetupExpress-Debug: leaving...n=%d\n", n ENDFD;
  return ok;
}

int MapSetupExpress(MapType* I, bool interruptible)
{ /* setup a list of neighbors for each square */
  PyMOLGlobals* G = I->G;
  int n = 0;
  int c, d, e, f, i, cm1, cp2, D1D2 = I->D1D2, D2 = I->Dim[2];
  int mx2 = I->iMax[2];
  int* link = I->Link.data();
  int st, flag;
  int *i_ptr3, *i_ptr4, *i_ptr5;
  int mx0 = I->iMax[0], mx1 = I->iMax[1], a, am1, ap2, *i_ptr1, b, bm1, bp2,
      *i_ptr2;
  int ok = true;

  PRINTFD(G, FB_Map)
  " MapSetupExpress-Debug: entered.\n" ENDFD;

  auto mapSize = I->Dim[0] * I->Dim[1] * I->Dim[2];
  I->EHead.assign(mapSize, 0);
  std::vector<int> e_list;
  e_list.reserve(1000);

  e_list.push_back(0);
  n = 1;
  for (a = (I->iMin[0] - 1); ok && a <= mx0; a++) {
    am1 = a - 1;
    ap2 = a + 2;
    i_ptr1 = I->Head.data() + am1 * D1D2;
    for (b = (I->iMin[1] - 1); ok && b <= mx1; b++) {
      bm1 = b - 1;
      bp2 = b + 2;
      i_ptr2 = i_ptr1 + bm1 * D2;
      for (c = (I->iMin[2] - 1); ok && c <= mx2; c++) {
        st = n;
        cm1 = c - 1;
        cp2 = c + 2;
        flag = false;
        i_ptr5 = (i_ptr4 = (i_ptr3 = i_ptr2 + cm1));

        for (d = am1; ok && d < ap2; d++) {
          for (e = bm1; ok && e < bp2; e++) {
            for (f = cm1; ok && f < cp2; f++) {
              /*i=*MapFirst(I,d,e,f); */
              if ((i = *(i_ptr5++)) >= 0) {
                flag = true;
                do {
                  e_list.push_back(i);
                  n++;
                  /*i=MapNext(I,i); */
                } while (ok && (i = link[i]) >= 0);
              }
              ok &= !(interruptible && G->Interrupt);
            }
            if (ok)
              i_ptr5 = (i_ptr4 += D2);
          }
          if (ok)
            i_ptr5 = (i_ptr4 = (i_ptr3 += D1D2));
        }
        if (ok) {
          if (flag) {
            *(MapEStart(I, a, b, c)) = st;
            e_list.push_back(-1);
            n++;
          } else {
            *(MapEStart(I, a, b, c)) = 0;
          }
        }
      }
    }
  }
  if (ok) {
    I->EList = std::move(e_list);
    I->EList.shrink_to_fit();
  } else {
    // Interrupted: EHead may already have starts into e_list, which is
    // discarded. Leave a consistent "no neighbors" table instead.
    I->EHead.assign(mapSize, 0);
    I->EList.clear();
  }
  PRINTFD(G, FB_Map)
  " MapSetupExpress-Debug: leaving...n=%d\n", n ENDFD;
  return ok;
}

/**
 * Get the grid indices for a 3D query point. If the point is outside the grid,
 * get the indices of the closest grid cell (clamp `v` to grid boundary).
 * @param v 3D query point
 * @param[out] a,b,c
 */
void MapLocus(const MapType* I, const float* v, int* a, int* b, int* c)
{
  int at, bt, ct;
  float invDiv = I->recipDiv;

  at = (int) ((v[0] - I->Min[0]) * invDiv) + MapBorder;
  bt = (int) ((v[1] - I->Min[1]) * invDiv) + MapBorder;
  ct = (int) ((v[2] - I->Min[2]) * invDiv) + MapBorder;

  *a = std::clamp(at, I->iMin[0], I->iMax[0]);
  *b = std::clamp(bt, I->iMin[1], I->iMax[1]);
  *c = std::clamp(ct, I->iMin[2], I->iMax[2]);
}

/**
 * Get EList start index for points in proximity to `v`.
 * Clamps `v` to grid boundaries.
 * @param v 3D query point
 * @return pointer to EList start index
 */
int* MapLocusEStart(MapType* I, const float* v)
{
  int a, b, c;
  MapLocus(I, v, &a, &b, &c);
  return MapEStart(I, a, b, c);
}

/**
 * Get the grid indices for a 3D query point. If the point is outside the grid,
 * return false and leave (a,b,c) in an unspecified state.
 * @param v 3D query point
 * @param[out] a,b,c Grid indices, but only if function returned true.
 * @return True if `v` is within grid boundaries
 */
static bool MapExclLocus(MapType* I, const float* v, int* a, int* b, int* c)
{
  float invDiv = I->recipDiv;

  *a = (int) (((v[0] - I->Min[0]) * invDiv) + MapBorder);
  if (*a < I->iMin[0])
    return (0);
  else if (*a > I->iMax[0])
    return (0);
  *b = (int) (((v[1] - I->Min[1]) * invDiv) + MapBorder);
  if (*b < I->iMin[1])
    return (0);
  else if (*b > I->iMax[1])
    return (0);
  *c = (int) (((v[2] - I->Min[2]) * invDiv) + MapBorder);
  if (*c < I->iMin[2])
    return (0);
  else if (*c > I->iMax[2])
    return (0);
  return (1);
}

/**
 * Return EList start index for points in proximity to `v`.
 * Return 0 if `v` is outside the grid or there are no points.
 * @param v 3D query point
 */
static int MapExclLocusEStart(MapType* map, const float* v)
{
  int h, k, l;
  if (!MapExclLocus(map, v, &h, &k, &l))
    return 0;
  return *(MapEStart(map, h, k, l));
}

float MapGetSeparation(PyMOLGlobals* G, float range, const float* mx,
    const float* mn, float* diagonal)
{
  float maxSize;
  float size, maxSubDiv, divSize, subDiv[3];
  float maxCubed, subDivCubed;
  int a;
  maxSize = SettingGetGlobal_i(G, cSetting_hash_max);

  maxCubed = maxSize * maxSize * maxSize;

  /* find longest axis: diagonal = max-min, for
   * each axis and find the largest */
  subtract3f(mx, mn, diagonal);
  diagonal[0] = (float) fabs(diagonal[0]);
  diagonal[1] = (float) fabs(diagonal[1]);
  diagonal[2] = (float) fabs(diagonal[2]);
  /* find largest */
  size = diagonal[0];
  if (diagonal[1] > size)
    size = diagonal[1];
  if (diagonal[2] > size)
    size = diagonal[2];
  /* err check size and diagonal */
  if (size == 0.0) {
    diagonal[0] = 1.0;
    diagonal[1] = 1.0;
    diagonal[2] = 1.0;
    size = 1.0;
  }

  /* compute maximum number of subdivisions */
  maxSubDiv = (float) (size / (range + MapSafety));
  if (maxSubDiv < 1.0F)
    maxSubDiv = 1.0F;

  /* find resulting divSize */
  divSize = size / maxSubDiv;
  if (divSize < MapSafety)
    divSize = MapSafety;
  for (a = 0; a < 3; a++) {
    subDiv[a] = (float) ((int) ((diagonal[a] / divSize) + 0.5F));
    subDiv[a] = (subDiv[a] < 1.0F) ? 1.0F : subDiv[a];
  }
  subDivCubed = subDiv[0] * subDiv[1] * subDiv[2];

  if (subDivCubed > maxCubed) {
    divSize = (float) (divSize / pow(maxCubed / subDivCubed, 0.33333F));
  } else if (subDivCubed < maxCubed) {
    divSize = (float) (divSize * pow(subDivCubed / maxCubed, 0.33333F));
  }

  if (divSize < (range + MapSafety))
    divSize = range + MapSafety;

  if (Feedback(G, FB_Map, FB_Debugging)) {
    PRINTF
    " MapGetSeparation: range %8.3f divSize %8.3f size %8.3f\n", range, divSize,
        size ENDF(G);
    /*    dump3f(mx,"mx");
       dump3f(mn,"mn");
       dump3f(diagonal,"diagonal");
       printf("%8.3f\n",divSize);
       printf("divSize %8.3f\n",divSize);
     */
  }
  return (divSize);
}

MapType::MapType(PyMOLGlobals* G, float range, const float* vert, int nVert,
    const float* extent, const int* flag)
{
  auto I = this;
  int mapSize;
  int h, k, l;
  int* list;
  const float* v;
  int firstFlag;
  Vector3f diagonal;
  int ok = true;

  PRINTFD(G, FB_Map)
  " MapNew-Debug: entered.\n" ENDFD;
  /* Initialize */
  I->G = G;

  // empirical limit to avoid crash in PYMOL-3002
  const float SANITY_LIMIT = 1e10;

  // Insane vertices are ignored in Maps now.
  auto is_sane = [SANITY_LIMIT](const float* vert) {
    return vert[0] > -SANITY_LIMIT && vert[0] < SANITY_LIMIT &&
           vert[1] > -SANITY_LIMIT && vert[1] < SANITY_LIMIT &&
           vert[2] > -SANITY_LIMIT && vert[2] < SANITY_LIMIT;
  };

  /* initialize an empty cache for the map */
  I->Link.assign(nVert, -1);

  /* map extents; set if valid, otherwise determine based on the flagged
   * vertices */
  if (extent) {
    /* valid, so copy */
    I->Min[0] = extent[0];
    I->Max[0] = extent[1];
    I->Min[1] = extent[2];
    I->Max[1] = extent[3];
    I->Min[2] = extent[4];
    I->Max[2] = extent[5];
  } else {
    /* blank, so determine */
    I->Min[0] = 0.0F;
    I->Max[0] = 0.0F;
    I->Min[1] = 0.0F;
    I->Max[1] = 0.0F;
    I->Min[2] = 0.0F;
    I->Max[2] = 0.0F;

    /* flag is an array of ints, one per vertex to signify inclusion for
     * consideration in this map */
    if (flag) {
      firstFlag = true;
      v = vert;
      for (int a = 0; a < nVert; a++) {
        /* if we consider this vertex */
        if (flag[a]) {
          /* first-time setup*/
          if (firstFlag) {
            for (int vertIdx = 0; vertIdx < nVert; vertIdx++) {
              if (is_sane(v)) {
                std::copy_n(v, 3, I->Min);
                std::copy_n(v, 3, I->Max);
                break;
              }
            }
            firstFlag = false;
          } else {
            /* min/max extents, over all vertices */
            for (int c = 0; c < 3; c++) {
              if (is_sane(v)) {
                if (I->Min[c] > v[c])
                  I->Min[c] = v[c];
                if (I->Max[c] < v[c])
                  I->Max[c] = v[c];
              }
            }
          }
        }
        v += 3;
      }
    } else {
      /* no flag: do all vertices in the list */
      if (nVert) {
        v = vert;
        int a = 0;

        while (a++ < nVert) {
          if (is_sane(v)) {
            std::copy_n(v, 3, I->Min);
            std::copy_n(v, 3, I->Max);
            break;
          }
        }

        v += a * 3;
        for (a = 1; a < nVert; a++) {
          for (int c = 0; c < 3; c++) {
            if (is_sane(v)) {
              if (I->Min[c] > v[c])
                I->Min[c] = v[c];
              if (I->Max[c] < v[c])
                I->Max[c] = v[c];
            }
          }
          v += 3;
        }
      }
    }
  }

  /* sanity check */
  for (int a = 0; a < 3; a++) {
    if (I->Min[a] > I->Max[a]) {
      std::swap(I->Max[a], I->Min[a]);
    }

    if (I->Min[a] < -SANITY_LIMIT) {
      PRINTFB(G, FB_Map, FB_Warnings)
      " %s-Warning: clamping Min %e -> %e\n", __FUNCTION__, I->Min[a],
          -SANITY_LIMIT ENDFB(G);
      I->Min[a] = -SANITY_LIMIT;
    }
    if (I->Max[a] > SANITY_LIMIT) {
      PRINTFB(G, FB_Map, FB_Warnings)
      " %s-Warning: clamping Max %e -> %e\n", __FUNCTION__, I->Max[a],
          SANITY_LIMIT ENDFB(G);
      I->Max[a] = SANITY_LIMIT;
    }
  }

  if (Feedback(G, FB_Map, FB_Debugging)) {
    printf(" MapSetup: %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f\n", I->Min[0],
        I->Min[1], I->Min[2], I->Max[0], I->Max[1], I->Max[2]);
  }
  /* interesting */
  for (int c = 0; c < 3; c++) {
    I->Min[c] -= MapSafety;
    I->Max[c] += MapSafety;
  }
  /* pad the boundaries by "range" */
  if (range <
      0.0) { /* negative range is a flag to expand edges using "range". */
    range = -range;
    for (int c = 0; c < 3; c++) {
      I->Min[c] -= range;
      I->Max[c] += range;
    }
  }

  /* compute final box size ..................... */
  I->Div = MapGetSeparation(G, range, I->Max, I->Min, diagonal);
  I->recipDiv = 1.0F / (I->Div); /* cache this */

  /* add borders to avoid special edge cases */
  I->Dim[0] = (int) ((diagonal[0] / I->Div) + 1 + (2 * MapBorder));
  I->Dim[1] = (int) ((diagonal[1] / I->Div) + 1 + (2 * MapBorder));
  I->Dim[2] = (int) ((diagonal[2] / I->Div) + 1 + (2 * MapBorder));

  if (Feedback(G, FB_Map, FB_Debugging)) {
    printf(" MapSetup: nVert: %d\n", nVert);
    printf(" MapSetup: I->Div: %8.3f\n", I->Div);
    printf(" MapSetup: %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f\n", I->Min[0],
        I->Min[1], I->Min[2], I->Max[0], I->Max[1], I->Max[2]);
    printf(" MapSetup: %8d %8d %8d\n", I->Dim[0], I->Dim[1], I->Dim[2]);
  }

  I->D1D2 = I->Dim[1] * I->Dim[2];

  I->iMin[0] = MapBorder;
  I->iMin[1] = MapBorder;
  I->iMin[2] = MapBorder;

  I->iMax[0] = I->Dim[0] - (1 + MapBorder);
  I->iMax[1] = I->Dim[1] - (1 + MapBorder);
  I->iMax[2] = I->Dim[2] - (1 + MapBorder);

  /* compute size and allocate */
  mapSize = I->Dim[0] * I->Dim[1] * I->Dim[2];
  I->Head.assign(mapSize, -1);

  I->NVert = nVert;

  PRINTFD(G, FB_Map)
  " MapNew-Debug: creating 3D hash...\n" ENDFD;

  /* create 3-D hash of the vertices */
  if (flag) {
    v = vert;
    for (int a = 0; a < nVert; a++) {
      if (flag[a])
        if (MapExclLocus(I, v, &h, &k, &l)) {
          list = MapFirst(I, h, k, l);
          I->Link[a] = *list;
          *list = a; /*add to top of list */
        }
      v += 3;
    }
  } else {
    v = vert;
    for (int a = 0; a < nVert; a++) {
      if (MapExclLocus(I, v, &h, &k, &l)) {
        list = MapFirst(I, h, k, l);
        /*        printf("LINK %d %d %d %d %5.2f %5.2f
         * %5.2f\n",a,h,k,l,v[0],v[1],v[2]); */
        I->Link[a] = *list;
        *list = a; /*add to top of list */
      }
      v += 3;
    }
  }

  PRINTFD(G, FB_Map)
  " MapNew-Debug: leaving...\n" ENDFD;
}

MapEIter::MapEIter(MapType& map, const float* v, bool excl)
{
  if (map.EList.empty()) {
    // Not interruptible: callers can't tell an interrupted setup from
    // "no neighbors" and would silently produce wrong results.
    MapSetupExpress(&map, false);
  }

  m_elist = map.EList.data();

  if (excl) {
    m_i = MapExclLocusEStart(&map, v);
  } else {
    m_i = *MapLocusEStart(&map, v);
  }
}

/**
 * True if `v_query` is within `cutoff` of any point in the map.
 *
 * @param map A hash map
 * @param v_map The points used to build the map
 * @param v_query A query point
 * @param cutoff The distance cutoff
 */
bool MapAnyWithin(
    MapType& map, const float* v_map, const float* v_query, float cutoff)
{
  for (const auto j : MapEIter(map, v_query)) {
    if (within3f(v_map + 3 * j, v_query, cutoff)) {
      return true;
    }
  }
  return false;
}

int MapType::size() const
{
  return EList.size();
}
