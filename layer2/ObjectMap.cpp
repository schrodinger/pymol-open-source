
/* 
A* -------------------------------------------------------------------
B* This file contains source code for the PyMOL computer program
C* copyright 1998-2000 by Warren Lyford Delano of DeLano Scientific. 
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

#include <stdint.h>
#include <algorithm>

#include"os_python.h"
#include"os_numpy.h"
#include"os_std.h"
#include"os_gl.h"

#include"OOMac.h"
#include"ObjectMap.h"
#include"Base.h"
#include"MemoryDebug.h"
#include"Map.h"
#include"Parse.h"
#include"Isosurf.h"
#include"Vector.h"
#include"Color.h"
#include"main.h"
#include"Scene.h"
#include"PConv.h"
#include"Word.h"
#include"Vector.h"
#include"PyMOLGlobals.h"
#include"Matrix.h"
#include"Util.h"
#include"ShaderMgr.h"
#include"CGO.h"
#include"File.h"
#include"Executive.h"
#include"Field.h"

#define n_space_group_numbers 231
static const char * space_group_numbers[] = {
"",
"P 1",         "P -1",        "P 2",         "P 21",        "C 2",
"P M",         "P C",         "C M",         "C C",         "P 2/M",
"P 21/M",      "C 2/M",       "P 2/C",       "P 21/C",      "C 2/C",
"P 2 2 2",     "P 2 2 21",    "P 21 21 2",   "P 21 21 21",  "C 2 2 21",
"C 2 2 2",     "F 2 2 2",     "I 2 2 2",     "I 21 21 21",  "P M M 2",
"P M C 21",    "P C C 2",     "P M A 2",     "P C A 21",    "P N C 2",
"P M N 21",    "P B A 2",     "P N A 21",    "P N N 2",     "C M M 2",
"C M C 21",    "C C C 2",     "A M M 2",     "A B M 2",     "A M A 2",
"A B A 2",     "F M M 2",     "F D D 2",     "I M M 2",     "I B A 2",
"I M A 2",     "P M M M",     "P N N N",     "P C C M",     "P B A N",
"P M M A",     "P N N A",     "P M N A",     "P C C A",     "P B A M",
"P C C N",     "P B C M",     "P N N M",     "P M M N",     "P B C N",
"P B C A",     "P N M A",     "C M C M",     "C M C A",     "C M M M",
"C C C M",     "C M M A",     "C C C A",     "F M M M",     "F D D D",
"I M M M",     "I B A M",     "I B C A",     "I M M A",     "P 4",
"P 41",        "P 42",        "P 43",        "I 4",         "I 41",
"P -4",        "I -4",        "P 4/M",       "P 42/M",      "P 4/N",
"P 42/N",      "I 4/M",       "I 41/A",      "P 4 2 2",     "P 4 21 2",
"P 41 2 2",    "P 41 21 2",   "P 42 2 2",    "P 42 21 2",   "P 43 2 2",
"P 43 21 2",   "I 4 2 2",     "I 41 2 2",    "P 4 M M",     "P 4 B M",
"P 42 C M",    "P 42 N M",    "P 4 C C",     "P 4 N C",     "P 42 M C",
"P 42 B C",    "I 4 M M",     "I 4 C M",     "I 41 M D",    "I 41 C D",
"P -4 2 M",    "P -4 2 C",    "P -4 21 M",   "P -4 21 C",   "P -4 M 2",
"P -4 C 2",    "P -4 B 2",    "P -4 N 2",    "I -4 M 2",    "I -4 C 2",
"I -4 2 M",    "I -4 2 D",    "P 4/M M M",   "P 4/M C C",   "P 4/N B M",
"P 4/N N C",   "P 4/M B M",   "P 4/M N C",   "P 4/N M M",   "P 4/N C C",
"P 42/M M C",  "P 42/M C M",  "P 42/N B C",  "P 42/N N M",  "P 42/M B C",
"P 42/M N M",  "P 42/N M C",  "P 42/N C M",  "I 4/M M M",   "I 4/M C M",
"I 41/A M D",  "I 41/A C D",  "P 3",         "P 31",        "P 32",
"R 3",         "P -3",        "R -3",        "P 3 1 2",     "P 3 2 1",
"P 31 1 2",    "P 31 2 1",    "P 32 1 2",    "P 32 2 1",    "R 3 2",
"P 3 M 1",     "P 3 1 M",     "P 3 C 1",     "P 3 1 C",     "R 3 M",
"R 3 C",       "P -3 1 M",    "P -3 1 C",    "P -3 M 1",    "P -3 C 1",
"R -3 M",      "R -3 C",      "P 6",         "P 61",        "P 65",
"P 62",        "P 64",        "P 63",        "P -6",        "P 6/M",
"P 63/M",      "P 6 2 2",     "P 61 2 2",    "P 65 2 2",    "P 62 2 2",
"P 64 2 2",    "P 63 2 2",    "P 6 M M",     "P 6 C C",     "P 63 C M",
"P 63 M C",    "P -6 M 2",    "P -6 C 2",    "P -6 2 M",    "P -6 2 C",
"P 6/M M M",   "P 6/M C C",   "P 63/M C M",  "P 63/M M C",  "P 2 3",
"F 2 3",       "I 2 3",       "P 21 3",      "I 21 3",      "P M -3",
"P N -3",      "F M -3",      "F D -3",      "I M -3",      "P A -3",
"I A -3",      "P 4 3 2",     "P 42 3 2",    "F 4 3 2",     "F 41 3 2",
"I 4 3 2",     "P 43 3 2",    "P 41 3 2",    "I 41 3 2",    "P -4 3 M",
"F -4 3 M",    "I -4 3 M",    "P -4 3 N",    "F -4 3 C",    "I -4 3 D",
"P M -3 M",    "P N -3 N",    "P M -3 N",    "P N -3 M",    "F M -3 M",
"F M -3 C",    "F D -3 M",    "F D -3 C",    "I M -3 M",    "I A -3 D",
};

/* MapState::ValidXtal -- determines whether the MapState's Xtal type passed in is valid
 * PARAMS
 *  ms, MapState
 * RETURNS
 *   true/false
 */
int ObjectMapStateValidXtal(ObjectMapState * ms)
{
  if(ms && ms->Active) {
    switch (ms->MapSource) {
      /* these are the only formats which are known to contain xtal
         information */
    case cMapSourceCrystallographic:
    case cMapSourceCCP4:
    case cMapSourceBRIX:
    case cMapSourceGRD:
      return true;
      break;

    }
  }
  return false;
}

/* Map::ValidXtal -- detemines whether the Map is valid
 * PARAMS
 *  I, Map
 *  state, map's state
 * RETURNS
 *  true/false
 */
int ObjectMapValidXtal(ObjectMap * I, int state)
{
  if((state >= 0) && (state < I->State.size())) {
    ObjectMapState *ms = &I->State[state];
    return ObjectMapStateValidXtal(ms);
  }
  return false;
}

/* MapState::GetExcludedStats -- 
 * PARARMS
 *   G, usual PyMOL globals
 *   ms, MapState
 *   vert_vla, variable length array of vertices
 *   beyond, radius of exclusion
 *   within, radius of inlcusion
 *   level, map level
 *
 * RETURNS
 *  number of exluded pts?
 */
int ObjectMapStateGetExcludedStats(PyMOLGlobals * G, ObjectMapState * ms, float *vert_vla,
                                   float beyond, float within, float *level)
{
  double sum = 0.0, sumsq = 0.0;
  float mean, stdev;
  int cnt = 0;
  int list_size;
  float cutoff = beyond;
  MapType *voxelmap = NULL;

  /* size of the VLA */
  if(vert_vla) {
    list_size = VLAGetSize(vert_vla) / 3;
  } else {
    list_size = 0;
  }
  if(cutoff < within)
    cutoff = within;

  /* make a new map from the VLA .............. */
  if(list_size)
    voxelmap = MapNew(G, -cutoff, vert_vla, list_size, NULL);

  if(voxelmap || (!list_size)) {
    int a, b, c;
    int h, k, l, i, j;
    const int *fdim = ms->FDim;
    int within_flag, within_default = false;
    int beyond_flag;

    const Isofield *field = ms->Field.get();
    if(list_size)
      MapSetupExpress(voxelmap);

    within_flag = true;
    beyond_flag = true;

    if(within < R_SMALL4)
      within_default = true;
    for(c = 0; c < fdim[2]; c++) {
      for(b = 0; b < fdim[1]; b++) {
        for(a = 0; a < fdim[0]; a++) {
          if(list_size) {
            within_flag = within_default;
            beyond_flag = true;

            const float* v = F4Ptr(field->points, a, b, c, 0);

            MapLocus(voxelmap, v, &h, &k, &l);
            i = *(MapEStart(voxelmap, h, k, l));
            if(i) {
              j = voxelmap->EList[i++];
              while(j >= 0) {
                if(!within_flag) {
                  if(within3f(vert_vla + 3 * j, v, within)) {
                    within_flag = true;
                  }
                }
                if(within3f(vert_vla + 3 * j, v, beyond)) {
                  beyond_flag = false;
                  break;
                }
                j = voxelmap->EList[i++];
              }
            }
          }

          if(within_flag && beyond_flag) {      /* point isn't too close to any vertex */
            const float f_val = F3(field->data, a, b, c);
            sum += f_val;
            sumsq += (f_val * f_val);
            cnt++;
          }
        }
      }
    }
    if(voxelmap)
      MapFree(voxelmap);
  }
  if(cnt) {
    mean = (float) (sum / cnt);
    stdev = (float) sqrt1d((sumsq - (sum * sum / cnt)) / (cnt));
    level[1] = mean;
    level[0] = mean - stdev;
    level[2] = mean + stdev;
  }
  return cnt;
}

/* MapState::ObjectMapStateGetDataRange -- get min/max from the scalar field
 * PARAMS
 * G
 *    usual PyMOLGlobals
 * ms
 *    I
 * RETURNS
 * npts, but stores range in min/max
 */
int ObjectMapStateGetDataRange(PyMOLGlobals * G, ObjectMapState * ms, float *min,
                               float *max)
{
float max_val = 0.0F, min_val = 0.0F;
  CField *data = ms->Field->data.get();
  int cnt = data->dim[0] * data->dim[1] * data->dim[2];
  float *raw_data = (float *) data->data.data();
  if(cnt) {
    int a;
    min_val = (max_val = *(raw_data++));
    for(a = 1; a < cnt; a++) {
      double f_val = *(raw_data++);
      if(min_val > f_val)
        min_val = f_val;
      if(max_val < f_val)
        max_val = f_val;
    }
  }
  *min = min_val;
  *max = max_val;
  return cnt;
}

/* MapState::ObjectMapStateGetHistogram -- compute a map histogram
 *
 * INPUT PARAMS
 *
 * h_points - number of histogram points
 * limit - limit the data to (mean - limit * stdev, mean + limit * stdev)
 *         if limit <= 0, don't trim the histogram
 * min_arg, max_arg - limit the data, has priority over "limit" param;
 *         Unused if min_arg==max_arg
 *
 * OUTPUT PARAMS
 *
 * histogram - output histogram buffer. first four values are
 *             minimum, maximum, mean and stdev, followed by n_points
 *             of non-normalized histogram counts.
 */
int ObjectMapStateGetHistogram(PyMOLGlobals * G, ObjectMapState * ms, 
                               int n_points, float limit, float *histogram,
                               float min_arg, float max_arg)
{
  float max_val = 0.0f, min_val = 0.0f;
  float sum = 0.0f, sumsq = 0.0f;
  float min_his, max_his, irange, mean, stdev;
  int pos;
  CField *data = ms->Field->data.get();
  int cnt = data->dim[0] * data->dim[1] * data->dim[2];
  float *raw_data = (float *) data->data.data();
  if(cnt) {
    int a;

    // compute min/max/mean/stdev
    sum = min_val = (max_val = *(raw_data++));
    sumsq = sum*sum;
    for(a = 1; a < cnt; a++) {
      double f_val = *(raw_data++);
      if(min_val > f_val)
        min_val = f_val;
      if(max_val < f_val)
        max_val = f_val;
      sum += f_val;
      sumsq += f_val*f_val;
    }
    mean = (float) (sum / cnt);
    stdev = (float) sqrt1d((sumsq - (sum * sum / cnt)) / (cnt));

    // adjust min/max to limit
    if (min_arg != max_arg) {
      min_his = min_arg;
      max_his = max_arg;
    } else if (limit > 0.0F) {
      min_his = mean - limit*stdev;
      if (min_his < min_val)
        min_his = min_val;
      max_his = mean + limit*stdev;
      if (max_his > max_val)
        max_his = max_val;
    } else {
      min_his = min_val;
      max_his = max_val;
    }

    // Compute the histogram
    if(n_points > 0) {
      irange = (float)(n_points-1) / (max_his - min_his);
      for (a = 0; a < n_points; a++)
        histogram[a+4] = 0.0f;
      raw_data = (float *) data->data.data();
      for (a = 0; a < cnt; a++) {
        double f_val = *(raw_data++);
        pos = (int)(irange * (f_val-min_his));
        if (pos >= 0 && pos < n_points) {
          histogram[pos+4] += 1.0;
        }
      }
    }
    histogram[0] = min_his;
    histogram[1] = max_his;
    histogram[2] = mean;
    histogram[3] = stdev;
  } else {
    histogram[0] = 0.0;
    histogram[1] = 1.0;
    histogram[2] = 1.0;
    histogram[3] = 1.0;
  }
  return cnt;
}

/**
 * @see ObjectMapStateInterpolate for parameter description
 * @param state Object state (can be -2 for current state)
 * @return True if the map state exists and `result` has been written to
 */
int ObjectMapInterpolate(ObjectMap * I, int state, const float *array, float *result, int *flag,
                         int n)
{
  int ok = false;
  float txf_buffer[3];
  float *txf = txf_buffer;

  ObjectMapState *ms = ObjectMapGetState(I, state);

  if(ms) {
    const double *matrix = ObjectStateGetInvMatrix(ms);

    if(matrix) {
      /* we have to back-transform points */
      if(n > 1) {
        txf = pymol::malloc<float>(3 * n);
      }

      const float *src = array;
      float *dst = txf;
      array = txf;

      for (int nn = n; nn--; src += 3, dst += 3) {
        transform44d3f(matrix, src, dst);
      }
    }

    ok = ObjectMapStateInterpolate(ms, array, result, flag, n);
    ok = true;
  }

  if(txf != txf_buffer)
    FreeP(txf);

  return (ok);
}

static void ObjectMapStateTrim(PyMOLGlobals * G, ObjectMapState * ms,
                              float *mn, float *mx, int quiet)
{
  int div[3];
  int min[3];
  int fdim[4];
  int new_min[3], new_max[3], new_fdim[3];
  int a, b, c, d, e, f;
  float *vt, *vr;
  float v[3];
  float grid[3];
  Isofield *field;
  float orig_size = 1.0F;
  float new_size = 1.0F;

  if(ObjectMapStateValidXtal(ms)) {
    float tst[3], frac_tst[3];
    float frac_mn[3];
    float frac_mx[3];
    int hit_flag = false;

    /* compute the limiting box extents in fractional space */

    for(a = 0; a < 8; a++) {
      tst[0] = (a & 0x1) ? mn[0] : mx[0];
      tst[1] = (a & 0x2) ? mn[1] : mx[1];
      tst[2] = (a & 0x4) ? mn[2] : mx[2];
      transform33f3f(ms->Symmetry->Crystal.realToFrac(), tst, frac_tst);
      if(!a) {
        copy3f(frac_tst, frac_mn);
        copy3f(frac_tst, frac_mx);
      } else {
        for(b = 0; b < 3; b++) {
          frac_mn[b] = (frac_mn[b] > frac_tst[b]) ? frac_tst[b] : frac_mn[b];
          frac_mx[b] = (frac_mx[b] < frac_tst[b]) ? frac_tst[b] : frac_mx[b];
        }
      }
    }
    for(a = 0; a < 3; a++) {
      div[a] = ms->Div[a];
      min[a] = ms->Min[a];
      fdim[a] = ms->FDim[a];
    }
    fdim[3] = 3;

    {
      int first_flag[3] = { false, false, false };
      for(d = 0; d < 3; d++) {
        int tst_min, tst_max;
        float v_min, v_max;
        for(c = 0; c < (fdim[d] - 1); c++) {
          tst_min = c + min[d];
          tst_max = tst_min + 1;
          v_min = tst_min / ((float) div[d]);
          v_max = tst_max / ((float) div[d]);
          if(((v_min >= frac_mn[d]) && (v_min <= frac_mx[d])) ||
             ((v_max >= frac_mn[d]) && (v_max <= frac_mx[d]))) {
            if(!first_flag[d]) {
              first_flag[d] = true;
              new_min[d] = tst_min;
              new_max[d] = tst_max;
            } else {
              new_min[d] = (new_min[d] > tst_min) ? tst_min : new_min[d];
              new_max[d] = (new_max[d] < tst_max) ? tst_max : new_max[d];
            }
          }
        }
        new_fdim[d] = (new_max[d] - new_min[d]) + 1;
      }
      hit_flag = first_flag[0] && first_flag[1] & first_flag[2];
    }

    if(hit_flag)
      hit_flag = (new_fdim[0] != fdim[0]) || (new_fdim[1] != fdim[1]) ||
        (new_fdim[2] != fdim[2]);

    if(hit_flag) {
      orig_size = fdim[0] * fdim[1] * fdim[2];
      new_size = new_fdim[0] * new_fdim[1] * new_fdim[2];

      field = new Isofield(G, new_fdim);
      field->save_points = ms->Field->save_points;

      for(c = 0; c < new_fdim[2]; c++) {
        f = c + (new_min[2] - min[2]);
        for(b = 0; b < new_fdim[1]; b++) {
          e = b + (new_min[1] - min[1]);
          for(a = 0; a < new_fdim[0]; a++) {
            d = a + (new_min[0] - min[0]);
            vt = F4Ptr(field->points, a, b, c, 0);
            vr = F4Ptr(ms->Field->points, d, e, f, 0);
            copy3f(vr, vt);
            F3(field->data, a, b, c) = F3(ms->Field->data, d, e, f);
          }
        }
      }
      for(a = 0; a < 3; a++) {
        ms->Min[a] = new_min[a];
        ms->Max[a] = new_max[a];
        ms->FDim[a] = new_fdim[a];
      }
      ms->Field.reset(field);

      /* compute new extents */
      v[2] = (ms->Min[2]) / ((float) ms->Div[2]);
      v[1] = (ms->Min[1]) / ((float) ms->Div[1]);
      v[0] = (ms->Min[0]) / ((float) ms->Div[0]);

      auto const& cryst = ms->Symmetry->Crystal;

      transform33f3f(cryst.fracToReal(), v, ms->ExtentMin);

      v[2] = ((ms->FDim[2] - 1) + ms->Min[2]) / ((float) ms->Div[2]);
      v[1] = ((ms->FDim[1] - 1) + ms->Min[1]) / ((float) ms->Div[1]);
      v[0] = ((ms->FDim[0] - 1) + ms->Min[0]) / ((float) ms->Div[0]);

      transform33f3f(cryst.fracToReal(), v, ms->ExtentMax);

      /* new corner */
      {
        float vv[3];
        d = 0;
        for(c = 0; c < ms->FDim[2]; c += (ms->FDim[2] - 1)) {
          v[2] = (c + ms->Min[2]) / ((float) ms->Div[2]);
          for(b = 0; b < ms->FDim[1]; b += (ms->FDim[1] - 1)) {
            v[1] = (b + ms->Min[1]) / ((float) ms->Div[1]);
            for(a = 0; a < ms->FDim[0]; a += (ms->FDim[0] - 1)) {
              v[0] = (a + ms->Min[0]) / ((float) ms->Div[0]);
              transform33f3f(cryst.fracToReal(), v, vv);
              copy3f(vv, ms->Corner + 3 * d);
              d++;
            }
          }
        }
      }
    }
  } else {                      /* not a crystal map */
    int hit_flag = false;
    float *origin = ms->Origin.data();

    for(a = 0; a < 3; a++) {
      min[a] = ms->Min[a];
      grid[a] = ms->Grid[a];
      fdim[a] = ms->FDim[a];
    }
    fdim[3] = 3;

    {
      int first_flag[3] = { false, false, false };
      for(d = 0; d < 3; d++) {
        int tst_min, tst_max;
        float v_min, v_max;
        for(c = 0; c < (fdim[d] - 1); c++) {
          tst_min = c + min[d];
          tst_max = tst_min + 1;
          v_min = origin[d] + grid[d] * (tst_min);
          v_max = origin[d] + grid[d] * (tst_max);
          if(((v_min >= mn[d]) && (v_min <= mx[d])) ||
             ((v_max >= mn[d]) && (v_max <= mx[d]))) {
            if(!first_flag[d]) {
              first_flag[d] = true;
              hit_flag = true;
              new_min[d] = tst_min;
              new_max[d] = tst_max;
            } else {
              new_min[d] = (new_min[d] > tst_min) ? tst_min : new_min[d];
              new_max[d] = (new_max[d] < tst_max) ? tst_max : new_max[d];
            }
          }
        }
        new_fdim[d] = (new_max[d] - new_min[d]) + 1;
      }
      hit_flag = first_flag[0] && first_flag[1] & first_flag[2];
    }

    if(hit_flag)
      hit_flag = (new_fdim[0] != fdim[0]) || (new_fdim[1] != fdim[1]) ||
        (new_fdim[2] != fdim[2]);

    if(hit_flag) {

      orig_size = fdim[0] * fdim[1] * fdim[2];
      new_size = new_fdim[0] * new_fdim[1] * new_fdim[2];

      field = new Isofield(G, new_fdim);
      field->save_points = ms->Field->save_points;

      for(c = 0; c < new_fdim[2]; c++) {
        f = c + (new_min[2] - min[2]);
        for(b = 0; b < new_fdim[1]; b++) {
          e = b + (new_min[1] - min[1]);
          for(a = 0; a < new_fdim[0]; a++) {
            d = a + (new_min[0] - min[0]);
            vt = F4Ptr(field->points, a, b, c, 0);
            vr = F4Ptr(ms->Field->points, d, e, f, 0);
            copy3f(vr, vt);
            F3(field->data, a, b, c) = F3(ms->Field->data, d, e, f);
          }
        }
      }
      for(a = 0; a < 3; a++) {
        ms->Min[a] = new_min[a];
        ms->Max[a] = new_max[a];
        ms->FDim[a] = new_fdim[a];
        if(!ms->Dim.empty())
          ms->Dim[a] = new_fdim[a];
      }

      ms->Field.reset(field);

      for(e = 0; e < 3; e++) {
        ms->ExtentMin[e] = ms->Origin[e] + ms->Grid[e] * ms->Min[e];
        ms->ExtentMax[e] = ms->Origin[e] + ms->Grid[e] * ms->Max[e];
      }

      d = 0;
      for(c = 0; c < ms->FDim[2]; c += ms->FDim[2] - 1) {
        v[2] = ms->Origin[2] + ms->Grid[2] * (c + ms->Min[2]);

        for(b = 0; b < ms->FDim[1]; b += ms->FDim[1] - 1) {
          v[1] = ms->Origin[1] + ms->Grid[1] * (b + ms->Min[1]);

          for(a = 0; a < ms->FDim[0]; a += ms->FDim[0] - 1) {
            v[0] = ms->Origin[0] + ms->Grid[0] * (a + ms->Min[0]);
            copy3f(v, ms->Corner + 3 * d);
            d++;
          }
        }
      }
    }
  }
  if(!quiet) {
    PRINTFB(G, FB_ObjectMap, FB_Actions)
      " ObjectMap: Map volume reduced by %2.0f%%.\n",
      (100 * (orig_size - new_size)) / orig_size ENDFB(G);
  }
}

static void ObjectMapStateDouble(PyMOLGlobals * G, ObjectMapState * ms)
{
  int div[3];
  int min[3];
  int max[3];
  int fdim[4];
  int a, b, c;
  float v[3], vr[3];
  float *vt;
  float x, y, z;
  float grid[3];

  Isofield *field;

  if(ObjectMapStateValidXtal(ms)) {
    for(a = 0; a < 3; a++) {
      div[a] = ms->Div[a] * 2;
      min[a] = ms->Min[a] * 2;
      max[a] = ms->Max[a] * 2;
      fdim[a] = ms->FDim[a] * 2 - 1;
    }
    fdim[3] = 3;
    field = new Isofield(G, fdim);
    field->save_points = ms->Field->save_points;
    for(c = 0; c < fdim[2]; c++) {
      v[2] = (c + min[2]) / ((float) div[2]);
      z = (c & 0x1) ? 0.5F : 0.0F;
      for(b = 0; b < fdim[1]; b++) {
        v[1] = (b + min[1]) / ((float) div[1]);
        y = (b & 0x1) ? 0.5F : 0.0F;
        for(a = 0; a < fdim[0]; a++) {
          v[0] = (a + min[0]) / ((float) div[0]);
          transform33f3f(ms->Symmetry->Crystal.fracToReal(), v, vr);
          x = (a & 0x1) ? 0.5F : 0.0F;
          vt = F4Ptr(field->points, a, b, c, 0);
          copy3f(vr, vt);
          if((a & 0x1) || (b & 0x1) || (c & 0x1)) {
            F3(field->data, a, b, c) = FieldInterpolatef(ms->Field->data.get(),
                                                         a / 2, b / 2, c / 2, x, y, z);
          } else {
            F3(field->data, a, b, c) = F3(ms->Field->data, a / 2, b / 2, c / 2);
          }
        }
      }
    }
    for(a = 0; a < 3; a++) {
      ms->Min[a] = min[a];
      ms->Max[a] = max[a];
      ms->FDim[a] = fdim[a];
      ms->Div[a] = div[a];
    }

    ms->Field.reset(field);
  } else {
    for(a = 0; a < 3; a++) {
      grid[a] = ms->Grid[a] / 2.0F;
      min[a] = ms->Min[a] * 2;
      max[a] = ms->Max[a] * 2;
      fdim[a] = ms->FDim[a] * 2 - 1;
    }
    fdim[3] = 3;

    field = new Isofield(G, fdim);
    field->save_points = ms->Field->save_points;

    for(c = 0; c < fdim[2]; c++) {
      v[2] = ms->Origin[2] + grid[2] * (c + min[2]);
      z = (c & 0x1) ? 0.5F : 0.0F;
      for(b = 0; b < fdim[1]; b++) {
        v[1] = ms->Origin[1] + grid[1] * (b + min[1]);
        y = (b & 0x1) ? 0.5F : 0.0F;
        for(a = 0; a < fdim[0]; a++) {
          v[0] = ms->Origin[0] + grid[0] * (a + min[0]);
          x = (a & 0x1) ? 0.5F : 0.0F;
          vt = F4Ptr(field->points, a, b, c, 0);
          copy3f(v, vt);
          if((a & 0x1) || (b & 0x1) || (c & 0x1)) {
            F3(field->data, a, b, c) = FieldInterpolatef(ms->Field->data.get(),
                                                         a / 2, b / 2, c / 2, x, y, z);
          } else {
            F3(field->data, a, b, c) = F3(ms->Field->data, a / 2, b / 2, c / 2);
          }
        }
      }
    }
    for(a = 0; a < 3; a++) {
      ms->Min[a] = min[a];
      ms->Max[a] = max[a];
      ms->FDim[a] = fdim[a];
      if(!ms->Dim.empty())
        ms->Dim[a] = fdim[a];
      if(!ms->Grid.empty())
        ms->Grid[a] = grid[a];
    }
    ms->Field.reset(field);
  }
}

static void ObjectMapStateHalve(PyMOLGlobals * G, ObjectMapState * ms, int smooth)
{
  int div[3];
  int min[3];
  int max[3];
  int fdim[4];
  int a, b, c;
  float v[3], vr[3];
  float *vt;
  float x, y, z;
  float grid[3];

  Isofield *field;

  if(ObjectMapStateValidXtal(ms)) {
    int *old_div, *old_min, *old_max;
    int a_2, b_2, c_2;

    for(a = 0; a < 3; a++) {
      div[a] = ms->Div[a] / 2;
      /* note that when Div is odd, a one-cell extrapolation will
         occur in order to preserve map size and get us onto a even
         number of sampling intervals for the cell */

      min[a] = ms->Min[a] / 2;
      max[a] = ms->Max[a] / 2;
      while((min[a] * 2) < ms->Min[a])
        min[a]++;
      while((max[a] * 2) > ms->Max[a])
        max[a]--;

      fdim[a] = (max[a] - min[a]) + 1;
    }
    fdim[3] = 3;
    old_div = ms->Div;
    old_min = ms->Min;
    old_max = ms->Max;

    if(smooth)
      FieldSmooth3f(ms->Field->data.get());

    field = new Isofield(G, fdim);
    field->save_points = ms->Field->save_points;

    /*
       printf("old_div %d %d %d\n",old_div[0],old_div[1],old_div[2]);
       printf("old_min %d %d %d\n",old_min[0],old_min[1],old_min[2]);
       printf("min %d %d %d\n",min[0],min[1],min[2]);
       printf("%d %d %d\n",ms->FDim[0],ms->FDim[1],ms->FDim[2]);
     */

    for(c = 0; c < fdim[2]; c++) {
      v[2] = (c + min[2]) / ((float) div[2]);
      c_2 = (c + min[2]) * 2 - old_min[2];
      z = (v[2] - ((c_2 + old_min[2]) / (float) old_div[2])) * old_div[2];
      if(c_2 >= old_max[2]) {
        c_2 = old_max[2] - 1;
        z = (v[2] - ((c_2 + old_min[2]) / (float) old_div[2])) * old_div[2];
      }
      for(b = 0; b < fdim[1]; b++) {
        v[1] = (b + min[1]) / ((float) div[1]);
        b_2 = (b + min[1]) * 2 - old_min[1];
        y = (v[1] - ((b_2 + old_min[1]) / (float) old_div[1])) * old_div[1];
        if(b_2 >= old_max[1]) {
          b_2 = old_max[1] - 1;
          y = (v[1] - ((b_2 + old_min[1]) / (float) old_div[1])) * old_div[1];
        }
        for(a = 0; a < fdim[0]; a++) {
          v[0] = (a + min[0]) / ((float) div[0]);
          a_2 = (a + min[0]) * 2 - old_min[0];
          x = (v[0] - ((a_2 + old_min[0]) / (float) old_div[0])) * old_div[0];
          if(a_2 >= old_max[0]) {
            a_2 = old_max[0] - 1;
            x = (v[0] - ((a_2 + old_min[0]) / (float) old_div[0])) * old_div[0];
          }
          transform33f3f(ms->Symmetry->Crystal.fracToReal(), v, vr);
          vt = F4Ptr(field->points, a, b, c, 0);
          copy3f(vr, vt);
          F3(field->data, a, b, c) = FieldInterpolatef(ms->Field->data.get(),
                                                       a_2, b_2, c_2, x, y, z);
        }
      }
    }
    for(a = 0; a < 3; a++) {
      ms->Min[a] = min[a];
      ms->Max[a] = max[a];
      ms->FDim[a] = fdim[a];
      ms->Div[a] = div[a];
    }

    ms->Field.reset(field);

    /* compute new extents */
    v[2] = (ms->Min[2]) / ((float) ms->Div[2]);
    v[1] = (ms->Min[1]) / ((float) ms->Div[1]);
    v[0] = (ms->Min[0]) / ((float) ms->Div[0]);

    transform33f3f(ms->Symmetry->Crystal.fracToReal(), v, ms->ExtentMin);

    v[2] = ((ms->FDim[2] - 1) + ms->Min[2]) / ((float) ms->Div[2]);
    v[1] = ((ms->FDim[1] - 1) + ms->Min[1]) / ((float) ms->Div[1]);
    v[0] = ((ms->FDim[0] - 1) + ms->Min[0]) / ((float) ms->Div[0]);

    transform33f3f(ms->Symmetry->Crystal.fracToReal(), v, ms->ExtentMax);

    /* new corner */
    {
      float vv[3];
      int d = 0;

      for(c = 0; c < ms->FDim[2]; c += (ms->FDim[2] - 1)) {
        v[2] = (c + ms->Min[2]) / ((float) ms->Div[2]);
        for(b = 0; b < ms->FDim[1]; b += (ms->FDim[1] - 1)) {
          v[1] = (b + ms->Min[1]) / ((float) ms->Div[1]);
          for(a = 0; a < ms->FDim[0]; a += (ms->FDim[0] - 1)) {
            v[0] = (a + ms->Min[0]) / ((float) ms->Div[0]);
            transform33f3f(ms->Symmetry->Crystal.fracToReal(), v, vv);
            copy3f(vv, ms->Corner + 3 * d);
            d++;
          }
        }
      }
    }
  } else {
    for(a = 0; a < 3; a++) {
      grid[a] = ms->Grid[a] * 2.0F;
      min[a] = ms->Min[a] / 2;
      max[a] = ms->Max[a] / 2;
      fdim[a] = (ms->FDim[a] + 1) / 2;
    }
    fdim[3] = 3;

    field = new Isofield(G, fdim);
    field->save_points = ms->Field->save_points;

    for(c = 0; c < fdim[2]; c++) {
      v[2] = ms->Origin[2] + grid[2] * (c + min[2]);
      for(b = 0; b < fdim[1]; b++) {
        v[1] = ms->Origin[1] + grid[1] * (b + min[1]);
        for(a = 0; a < fdim[0]; a++) {
          v[0] = ms->Origin[0] + grid[0] * (a + min[0]);
          vt = F4Ptr(field->points, a, b, c, 0);
          copy3f(v, vt);
          F3(field->data, a, b, c) = F3(ms->Field->data, a * 2, b * 2, c * 2);
        }
      }
    }
    for(a = 0; a < 3; a++) {
      ms->Min[a] = min[a];
      ms->Max[a] = max[a];
      ms->FDim[a] = fdim[a];
      if(!ms->Dim.empty())
        ms->Dim[a] = fdim[a];
      if(!ms->Grid.empty())
        ms->Grid[a] = grid[a];
    }
    ms->Field.reset(field);

  }
}

pymol::Result<> ObjectMapTrim(
    ObjectMap* I, int state, float* mn, float* mx, int quiet)
{
  int update = false;

  /* TO DO: convert mn and mx into map local coordinates if map itself is transformed...  */

  if(state < 0) {
    for(auto& ms : I->State) {
      if(ms.Active) {
        ObjectMapStateTrim(I->G, &ms, mn, mx, quiet);
          update = true;
      }
    }
  } else if((state >= 0) && (state < I->State.size()) && (I->State[state].Active)) {
    ObjectMapStateTrim(I->G, &I->State[state], mn, mx, quiet);
  } else {
    return pymol::make_error("Invalid state.");
  }
  if(update)
    ObjectMapUpdateExtents(I);
  return {};
}

pymol::Result<> ObjectMapDouble(ObjectMap* I, int state)
{
  if(state < 0) {
    for(auto& state : I->State) {
      if(state.Active)
        ObjectMapStateDouble(I->G, &state);
    }
  } else if((state >= 0) && (state < I->State.size()) && (I->State[state].Active)) {
    ObjectMapStateDouble(I->G, &I->State[state]);
  } else {
    return pymol::make_error("Invalidate state.");
  }
  return {};
}

pymol::Result<> ObjectMapHalve(ObjectMap * I, int state, int smooth)
{
  if(state < 0) {
    for(auto& state : I->State) {
      if(state.Active)
        ObjectMapStateHalve(I->G, &state, smooth);
    }
  } else if((state >= 0) && (state < I->State.size()) && (I->State[state].Active)) {
    ObjectMapStateHalve(I->G, &I->State[state], smooth);
  } else {
    return pymol::make_error("Invalidate state.");
  }
  ObjectMapUpdateExtents(I);

  return {};
}

int ObjectMapStateContainsPoint(ObjectMapState * ms, float *point)
{
  int result = false;
  float x, y, z;
  int x_floor, y_floor, z_floor;
  int x_ceil, y_ceil, z_ceil;

  if(ObjectMapStateValidXtal(ms)) {
    float frac[3];

    transform33f3f(ms->Symmetry->Crystal.realToFrac(), point, frac);

    x = (ms->Div[0] * frac[0]);
    y = (ms->Div[1] * frac[1]);
    z = (ms->Div[2] * frac[2]);
    x_floor = floor(x);
    x_ceil = ceil(x);
    y_floor = floor(y);
    y_ceil = ceil(y);
    z_floor = floor(z);
    z_ceil = ceil(z);

    if((x_floor >= ms->Min[0]) && (x_ceil <= ms->Max[0]) &&
       (y_floor >= ms->Min[1]) && (y_ceil <= ms->Max[1]) &&
       (z_floor >= ms->Min[2]) && (z_ceil <= ms->Max[2]))
      result = true;
  } else {
    x = (point[0] - ms->Origin[0]) / ms->Grid[0];
    y = (point[1] - ms->Origin[1]) / ms->Grid[1];
    z = (point[2] - ms->Origin[2]) / ms->Grid[2];

    x_floor = floor(x);
    x_ceil = ceil(x);
    y_floor = floor(y);
    y_ceil = ceil(y);
    z_floor = floor(z);
    z_ceil = ceil(z);

    if((x_floor >= ms->Min[0]) && (x_ceil <= ms->Max[0]) &&
       (y_floor >= ms->Min[1]) && (y_ceil <= ms->Max[1]) &&
       (z_floor >= ms->Min[2]) && (z_ceil <= ms->Max[2]))
      result = true;

    if((x >= ms->Min[0]) && (x <= ms->Max[0]) &&
       (y >= ms->Min[1]) && (y <= ms->Max[1]) && (z >= ms->Min[2]) && (z <= ms->Max[2]))
      result = true;
  }
  return (result);
}

/**
 * @param array Coordinate array of length `3*n`
 * @param[out] result Array of length `n` which will be populated with map
 * values at coorinate positions (linear interpolated between grid points)
 * @param[out] flag Array of length `n` which will be populated with booleans
 * indicating if points were within map bounds (optional, can be NULL)
 * @return False if any coordinate was out of bounds
 */
int ObjectMapStateInterpolate(ObjectMapState * ms, const float *array, float *result, int *flag,
                              int n)
{
  int ok = true;
  const float *inp;
  int a, b, c;
  float x, y, z;
  inp = array;

  if(ObjectMapStateValidXtal(ms)) {
    float frac[3];

    while(n--) {
      /* get the fractional coordinate */
      transform33f3f(ms->Symmetry->Crystal.realToFrac(), inp, frac);

      inp += 3;

      /* compute the effective lattice offset as a function of cell spacing */

      x = (ms->Div[0] * frac[0]);
      y = (ms->Div[1] * frac[1]);
      z = (ms->Div[2] * frac[2]);

      /* now separate the integral and fractional parts for interpolation */

      a = (int) floor(x + R_SMALL8);
      b = (int) floor(y + R_SMALL8);
      c = (int) floor(z + R_SMALL8);
      x -= a;
      y -= b;
      z -= c;

      if(flag)
        *flag = 1;
      /* wrap into the map */

      if(a < ms->Min[0]) {
        if(x < 0.99F) {
          ok = false;
          if(flag)
            *flag = 0;
        }
        x = 0.0F;
        a = ms->Min[0];
      } else if(a >= ms->FDim[0] + ms->Min[0] - 1) {
        if(x > 0.01F) {
          ok = false;
          if(flag)
            *flag = 0;
        }
        x = 0.0F;
        a = ms->FDim[0] + ms->Min[0] - 1;
      }

      if(b < ms->Min[1]) {
        if(y < 0.99F) {
          ok = false;
          if(flag)
            *flag = 0;
        }
        y = 0.0F;
        b = ms->Min[1];
      } else if(b >= ms->FDim[1] + ms->Min[1] - 1) {
        if(y > 0.01F) {
          ok = false;
          if(flag)
            *flag = 0;
        }
        y = 0.0F;
        b = ms->FDim[1] + ms->Min[1] - 1;
      }

      if(c < ms->Min[2]) {
        if(z < 0.99F) {
          ok = false;
          if(flag)
            *flag = 0;
        }
        z = 0.0F;
        c = ms->Min[2];
      } else if(c >= ms->FDim[2] + ms->Min[2] - 1) {
        if(z > 0.01) {
          ok = false;
          if(flag)
            *flag = 0;
        }
        z = 0.0F;
        c = ms->FDim[2] + ms->Min[2] - 1;
      }
      /*      printf("%d %d %d %8.3f %8.3f %8.3f\n",a,b,c,x,y,z); */
      *(result++) = FieldInterpolatef(ms->Field->data.get(),
                                      a - ms->Min[0],
                                      b - ms->Min[1], c - ms->Min[2], x, y, z);
      if(flag)
        flag++;

    }
  } else {

    while(n--) {

      x = (inp[0] - ms->Origin[0]) / ms->Grid[0];
      y = (inp[1] - ms->Origin[1]) / ms->Grid[1];
      z = (inp[2] - ms->Origin[2]) / ms->Grid[2];
      inp += 3;

      a = (int) floor(x + R_SMALL8);
      b = (int) floor(y + R_SMALL8);
      c = (int) floor(z + R_SMALL8);
      x -= a;
      y -= b;
      z -= c;

      if(flag)
        *flag = 1;
      if(a < ms->Min[0]) {
        x = 0.0F;
        a = ms->Min[0];
        ok = false;
        if(flag)
          *flag = 0;
      } else if(a >= ms->Max[0]) {
        x = 1.0F;
        a = ms->Max[0] - 1;
        ok = false;
        if(flag)
          *flag = 0;
      }

      if(b < ms->Min[1]) {
        y = 0.0F;
        b = ms->Min[1];
        ok = false;
        if(flag)
          *flag = 0;
      } else if(b >= ms->Max[1]) {
        y = 1.0F;
        b = ms->Max[1] - 1;
        ok = false;
        if(flag)
          *flag = 0;
      }

      if(c < ms->Min[2]) {
        z = 0.0F;
        c = ms->Min[2];
        ok = false;
        if(flag)
          *flag = 0;
      } else if(c >= ms->Max[2]) {
        z = 1.0F;
        c = ms->Max[2] - 1;
        ok = false;
        if(flag)
          *flag = 0;
      }
      /*      printf("%d %d %d %8.3f %8.3f %8.3f\n",a,b,c,x,y,z); */
      *(result++) = FieldInterpolatef(ms->Field->data.get(),
                                      a - ms->Min[0],
                                      b - ms->Min[1], c - ms->Min[2], x, y, z);
      if(flag)
        flag++;
    }
  }
  return (ok);
}

#ifndef _PYMOL_NOPY
static int ObjectMapNumPyArrayToMapState(PyMOLGlobals * G, ObjectMapState * I,
                                         PyObject * ary, int quiet);
#endif

void ObjectMapRegeneratePoints(ObjectMap * om){
  int i;
  for (i=0; i<om->State.size();i++){
    ObjectMapStateRegeneratePoints(&om->State[i]);
  }
}

void ObjectMapStateRegeneratePoints(ObjectMapState * ms)
{
  int a, b, c, e;
  float v[3], vr[3];

  if(ObjectMapStateValidXtal(ms)) {
    for(c = 0; c < ms->FDim[2]; c++) {
      v[2] = (c + ms->Min[2]) / ((float) ms->Div[2]);
      for(b = 0; b < ms->FDim[1]; b++) {
        v[1] = (b + ms->Min[1]) / ((float) ms->Div[1]);
        for(a = 0; a < ms->FDim[0]; a++) {
          v[0] = (a + ms->Min[0]) / ((float) ms->Div[0]);
          transform33f3f(ms->Symmetry->Crystal.fracToReal(), v, vr);
          for(e = 0; e < 3; e++)
            F4(ms->Field->points, a, b, c, e) = vr[e];
        }
      }
    }
  } else {
    for(c = 0; c < ms->FDim[2]; c++) {
      v[2] = ms->Origin[2] + ms->Grid[2] * (c + ms->Min[2]);
      for(b = 0; b < ms->FDim[1]; b++) {
        v[1] = ms->Origin[1] + ms->Grid[1] * (b + ms->Min[1]);
        for(a = 0; a < ms->FDim[0]; a++) {
          v[0] = ms->Origin[0] + ms->Grid[0] * (a + ms->Min[0]);
          for(e = 0; e < 3; e++) {
            F4(ms->Field->points, a, b, c, e) = v[e];
          }
        }
      }
    }
  }
}

static PyObject *ObjectMapStateAsPyList(ObjectMapState * I)
{
  PyObject *result = NULL;

  result = PyList_New(16);
  PyList_SetItem(result, 0, PyInt_FromLong(I->Active));
  if(I->Symmetry) {
    PyList_SetItem(result, 1, SymmetryAsPyList(I->Symmetry.get()));
  } else {
    PyList_SetItem(result, 1, PConvAutoNone(Py_None));
  }
  if(!I->Origin.empty()) {
    PyList_SetItem(result, 2, PConvFloatArrayToPyList(I->Origin.data(), 3));
  } else {
    PyList_SetItem(result, 2, PConvAutoNone(Py_None));
  }
  if(!I->Range.empty()) {
    PyList_SetItem(result, 3, PConvFloatArrayToPyList(I->Range.data(), 3));
  } else {
    PyList_SetItem(result, 3, PConvAutoNone(Py_None));
  }
  if(!I->Dim.empty()) {
    PyList_SetItem(result, 4, PConvIntArrayToPyList(I->Dim.data(), 3));
  } else {
    PyList_SetItem(result, 4, PConvAutoNone(Py_None));
  }
  if(!I->Grid.empty()) {
    PyList_SetItem(result, 5, PConvFloatArrayToPyList(I->Grid.data(), 3));
  } else {
    PyList_SetItem(result, 5, PConvAutoNone(Py_None));
  }
  PyList_SetItem(result, 6, PConvFloatArrayToPyList(I->Corner, 24));
  PyList_SetItem(result, 7, PConvFloatArrayToPyList(I->ExtentMin, 3));
  PyList_SetItem(result, 8, PConvFloatArrayToPyList(I->ExtentMax, 3));
  PyList_SetItem(result, 9, PyInt_FromLong(I->MapSource));

  PyList_SetItem(result, 10, PConvIntArrayToPyList(I->Div, 3));
  PyList_SetItem(result, 11, PConvIntArrayToPyList(I->Min, 3));
  PyList_SetItem(result, 12, PConvIntArrayToPyList(I->Max, 3));
  PyList_SetItem(result, 13, PConvIntArrayToPyList(I->FDim, 4));

  PyList_SetItem(result, 14, IsosurfAsPyList(I->G, I->Field.get()));
  PyList_SetItem(result, 15, ObjectStateAsPyList(I));
  return (PConvAutoNone(result));
}

static PyObject *ObjectMapAllStatesAsPyList(ObjectMap * I)
{
  PyObject *result = NULL;
  int a;
  result = PyList_New(I->State.size());
  for(a = 0; a < I->State.size(); a++) {
    if(I->State[a].Active) {
      PyList_SetItem(result, a, ObjectMapStateAsPyList(&I->State[a]));
    } else {
      PyList_SetItem(result, a, PConvAutoNone(NULL));
    }
  }
  return (PConvAutoNone(result));

}

static int ObjectMapStateCopy(const ObjectMapState * src, ObjectMapState * I)
{
  int ok = true;
  if(ok) {
    I->Active = src->Active;
    if(I->Active) {
      I->Symmetry = src->Symmetry;
      I->Origin = src->Origin;
      I->Range = src->Range;
      I->Grid = src->Grid;
      I->Dim = src->Dim;
      std::copy_n(src->Corner, 24, I->Corner);

      copy3f(src->ExtentMin, I->ExtentMin);
      copy3f(src->ExtentMax, I->ExtentMax);

      I->MapSource = src->MapSource;

      copy3f(src->Div, I->Div);
      copy3f(src->Min, I->Min);
      copy3f(src->Max, I->Max);
      copy3f(src->FDim, I->FDim);

      I->Field = src->Field;
      if(ok)
        ObjectMapStateRegeneratePoints(I);
    }
  }
  return (ok);
}

ObjectMapState::ObjectMapState(const ObjectMapState& src) : CObjectState(src)
{
  ObjectMapStateCopy(&src, this);
}

ObjectMapState& ObjectMapState::operator=(const ObjectMapState& src)
{
  CObjectState::operator=(src);
  ObjectMapStateCopy(&src, this);
  return *this;
}

static int ObjectMapStateFromPyList(PyMOLGlobals * G, ObjectMapState * I, PyObject * list)
{
  int ok = true;
  int ll = 0;
  PyObject *tmp;
  if(ok)
    ok = (list != NULL);
  if(ok) {
    if(!PyList_Check(list))
      I->Active = false;
    else {
      if(ok)
        ok = PyList_Check(list);
      if(ok)
        ll = PyList_Size(list);
      /* TO SUPPORT BACKWARDS COMPATIBILITY...
         Always check ll when adding new PyList_GetItem's */

      if(ok)
        ok = PConvPyIntToInt(PyList_GetItem(list, 0), &I->Active);
      if(ok) {
        tmp = PyList_GetItem(list, 1);
        if(tmp == Py_None)
          I->Symmetry = NULL;
        else {
          I->Symmetry.reset(SymmetryNewFromPyList(G, tmp));
          ok = I->Symmetry != nullptr;
        }
	CPythonVal_Free(tmp);
      }
      if(ok) {
        tmp = PyList_GetItem(list, 2);
        if(CPythonVal_IsNone(tmp))
          I->Origin.clear();
        else
          ok = PConvFromPyObject(G, tmp, I->Origin);
	CPythonVal_Free(tmp);
      }
      if(ok) {
        tmp = PyList_GetItem(list, 3);
        if(CPythonVal_IsNone(tmp))
          I->Range.clear();
        else
          ok = PConvFromPyObject(G, tmp, I->Range);
	CPythonVal_Free(tmp);
      }
      if(ok) {
        tmp = PyList_GetItem(list, 4);
        if(CPythonVal_IsNone(tmp))
          I->Dim.clear();
        else
          ok = PConvFromPyObject(G, tmp, I->Dim);
	CPythonVal_Free(tmp);
      }
      if(ok) {
        tmp = PyList_GetItem(list, 5);
        if(CPythonVal_IsNone(tmp))
          I->Grid.clear();
        else
          ok = PConvFromPyObject(G, tmp, I->Grid);
	CPythonVal_Free(tmp);
      }
      if(ok)
        ok = PConvPyListToFloatArrayInPlace(PyList_GetItem(list, 6), I->Corner, 24);
      if(ok)
        ok = PConvPyListToFloatArrayInPlace(PyList_GetItem(list, 7), I->ExtentMin, 3);
      if(ok)
        ok = PConvPyListToFloatArrayInPlace(PyList_GetItem(list, 8), I->ExtentMax, 3);
      if(ok)
        ok = PConvPyIntToInt(PyList_GetItem(list, 9), &I->MapSource);
      if(ok)
        ok = PConvPyListToIntArrayInPlace(PyList_GetItem(list, 10), I->Div, 3);
      if(ok)
        ok = PConvPyListToIntArrayInPlace(PyList_GetItem(list, 11), I->Min, 3);
      if(ok)
        ok = PConvPyListToIntArrayInPlace(PyList_GetItem(list, 12), I->Max, 3);
      if(ok)
        ok = CPythonVal_PConvPyListToIntArrayInPlace_From_List(G, list, 13, I->FDim, 4);
      if(ok){
	tmp = CPythonVal_PyList_GetItem(G, list, 14);
        I->Field.reset(IsosurfNewFromPyList(G, tmp));
        ok = I->Field != nullptr;
	CPythonVal_Free(tmp);
      }
      if(ok && (ll > 15)){
	tmp = CPythonVal_PyList_GetItem(G, list, 15);
        ok = ObjectStateFromPyList(G, tmp, I);
	CPythonVal_Free(tmp);
      }
      if(ok)
        ObjectMapStateRegeneratePoints(I);
    }
  }
  return (ok);
}

static int ObjectMapAllStatesFromPyList(ObjectMap * I, PyObject * list)
{
  int ok = true;
  int a;
  if(ok)
    ok = PyList_Check(list);
  if(ok) {
    I->State.resize(PyList_Size(list), ObjectMapState(I->G));
    for(a = 0; a < I->State.size(); a++) {
      CPythonVal *val = CPythonVal_PyList_GetItem(I->G, list, a);
      ok = ObjectMapStateFromPyList(I->G, &I->State[a], val);
      CPythonVal_Free(val);
      if(!ok)
        break;
    }
  }
  return (ok);
}

PyObject *ObjectMapAsPyList(ObjectMap * I)
{
  PyObject *result = NULL;

  result = PyList_New(3);
  PyList_SetItem(result, 0, ObjectAsPyList(I));
  PyList_SetItem(result, 1, PyInt_FromLong(I->State.size()));
  PyList_SetItem(result, 2, ObjectMapAllStatesAsPyList(I));

  return (PConvAutoNone(result));
}

int ObjectMapNewFromPyList(PyMOLGlobals * G, PyObject * list, ObjectMap ** result)
{
  int ok = true;
  ObjectMap *I = NULL;
  (*result) = NULL;

  if(ok)
    ok = (list != NULL);
  if(ok)
    ok = PyList_Check(list);
  /* TO SUPPORT BACKWARDS COMPATIBILITY...
     Always check ll when adding new PyList_GetItem's */
  I = new ObjectMap(G);
  if(ok)
    ok = (I != NULL);

  if(ok){
    auto *val = PyList_GetItem(list, 0);
    ok = ObjectFromPyList(G, val, I);
  }  
  if(ok){
    CPythonVal *val = CPythonVal_PyList_GetItem(G, list, 2);    
    ok = ObjectMapAllStatesFromPyList(I, val);
    CPythonVal_Free(val);
  }
  if(ok) {
    (*result) = I;
    ObjectMapUpdateExtents(I);
  } else {
    /* cleanup? */
  }

  return (ok);
}

int ObjectMapNewCopy(PyMOLGlobals * G, const ObjectMap * src, ObjectMap ** result,
                     int source_state, int target_state)
{
  int ok = true;
  ObjectMap *I = NULL;
  I = new ObjectMap(G);
  if(ok)
    ok = (I != NULL);
  if(ok)
    ok = ObjectCopyHeader(I, src);
  if(ok) {
    if(source_state == -1) {    /* all states */
      int state;
      VecCheckEmplace(I->State, I->State.size(), I->G);
      for(state = 0; state < src->State.size(); state++) {
        I->State[state] = src->State[state];
      }
    } else {
      if(target_state < 0)
        target_state = 0;
      if(source_state < 0)
        source_state = 0;
      VecCheckEmplace(I->State, target_state, G);
      if(source_state < src->State.size()) {
        I->State[target_state] = src->State[source_state];
      } else {
        ok = false;
        /* to do */
      }
    }
  }
  if(ok)
    *result = I;
  return ok;
}

ObjectMapState *ObjectMapStatePrime(ObjectMap * I, int state)
{
  if(state < 0)
    state = I->State.size();
  if(I->State.size() <= state) {
    VecCheckEmplace(I->State, state, I->G);
  }
  return &I->State[state];
}

void ObjectMapUpdateExtents(ObjectMap * I)
{
  int a;
  float *min_ext, *max_ext;
  float tr_min[3], tr_max[3];
  I->ExtentFlag = false;

  for(a = 0; a < I->State.size(); a++) {
    ObjectMapState *ms = &I->State[a];
    if(ms->Active) {
      if(!ms->Matrix.empty()) {
        transform44d3f(ms->Matrix.data(), ms->ExtentMin, tr_min);
        transform44d3f(ms->Matrix.data(), ms->ExtentMax, tr_max);
        {
          float tmp;
          int a;
          for(a = 0; a < 3; a++)
            if(tr_min[a] > tr_max[a]) {
              tmp = tr_min[a];
              tr_min[a] = tr_max[a];
              tr_max[a] = tmp;
            }
        }
        min_ext = tr_min;
        max_ext = tr_max;
      } else {
        min_ext = ms->ExtentMin;
        max_ext = ms->ExtentMax;
      }

      if(!I->ExtentFlag) {
        copy3f(min_ext, I->ExtentMin);
        copy3f(max_ext, I->ExtentMax);
        I->ExtentFlag = true;
      } else {
        min3f(min_ext, I->ExtentMin, I->ExtentMin);
        max3f(max_ext, I->ExtentMax, I->ExtentMax);
      }
    }
  }

  if(I->TTTFlag && I->ExtentFlag) {
    const float *ttt;
    double tttd[16];
    if(ObjectGetTTT(I, &ttt, -1)) {
      convertTTTfR44d(ttt, tttd);
      MatrixTransformExtentsR44d3f(tttd,
                                   I->ExtentMin, I->ExtentMax,
                                   I->ExtentMin, I->ExtentMax);
    }
  }

  PRINTFD(I->G, FB_ObjectMap)
    " ObjectMapUpdateExtents-DEBUG: ExtentFlag %d\n", I->ExtentFlag ENDFD;
}

void ObjectMapStateClamp(ObjectMapState * I, float clamp_floor, float clamp_ceiling)
{
  int a, b, c;
  float *fp;

  for(a = 0; a < I->FDim[0]; a++)
    for(b = 0; b < I->FDim[1]; b++)
      for(c = 0; c < I->FDim[2]; c++) {
        fp = F3Ptr(I->Field->data, a, b, c);
        if(*fp < clamp_floor)
          *fp = clamp_floor;
        else if(*fp > clamp_ceiling)
          *fp = clamp_ceiling;
      }
}

int ObjectMapStateSetBorder(ObjectMapState * I, float level)
{
  int result = true;
  int a, b, c;

  c = I->FDim[2] - 1;
  for(a = 0; a < I->FDim[0]; a++)
    for(b = 0; b < I->FDim[1]; b++) {
      F3(I->Field->data, a, b, 0) = level;
      F3(I->Field->data, a, b, c) = level;
    }

  a = I->FDim[0] - 1;
  for(b = 0; b < I->FDim[1]; b++)
    for(c = 0; c < I->FDim[2]; c++) {
      F3(I->Field->data, 0, b, c) = level;
      F3(I->Field->data, a, b, c) = level;
    }

  b = I->FDim[1] - 1;
  for(a = 0; a < I->FDim[0]; a++)
    for(c = 0; c < I->FDim[2]; c++) {
      F3(I->Field->data, a, 0, c) = level;
      F3(I->Field->data, a, b, c) = level;
    }
  return (result);
}

void ObjectMapStatePurge(PyMOLGlobals * G, ObjectMapState * I)
{
  ObjectStatePurge(I);
  I->Field = nullptr;
  I->Origin.clear();
  I->Dim.clear();
  I->Range.clear();
  I->Grid.clear();
  I->shaderCGO = nullptr;

  I->Symmetry = nullptr;

  I->Active = false;
}

void ObjectMap::update()
{
  auto I = this;
  if(!I->ExtentFlag) {
    ObjectMapUpdateExtents(I);
    if(I->ExtentFlag)
      SceneInvalidate(I->G);
  }
}

void ObjectMap::invalidate(cRep_t rep, cRepInv_t level, int state)
{
  auto I = this;
  if(level >= cRepInvExtents) {
    I->ExtentFlag = false;
  }
  if((rep < 0) || (rep == cRepDot)) {
    int a;
    for(a = 0; a < I->State.size(); a++) {
      if(I->State[a].Active)
        I->State[a].have_range = false;
      I->State[a].shaderCGO = nullptr;
    }
  }
  SceneInvalidate(I->G);
}

/* Has no prototype */
static CGO* ObjectMapCGOGenerate(PyMOLGlobals *G, float* corner)
{
  CGO cgo_value(G);
  CGO* const shaderCGO = &cgo_value;
  CGOBegin(shaderCGO, GL_LINES);

  CGOVertexv(shaderCGO, corner + 3 * 0);
  CGOVertexv(shaderCGO, corner + 3 * 1);

  CGOVertexv(shaderCGO, corner + 3 * 0);
  CGOVertexv(shaderCGO, corner + 3 * 2);

  CGOVertexv(shaderCGO, corner + 3 * 2);
  CGOVertexv(shaderCGO, corner + 3 * 3);

  CGOVertexv(shaderCGO, corner + 3 * 1);
  CGOVertexv(shaderCGO, corner + 3 * 3);

  CGOVertexv(shaderCGO, corner + 3 * 0);
  CGOVertexv(shaderCGO, corner + 3 * 4);

  CGOVertexv(shaderCGO, corner + 3 * 1);
  CGOVertexv(shaderCGO, corner + 3 * 5);

  CGOVertexv(shaderCGO, corner + 3 * 2);
  CGOVertexv(shaderCGO, corner + 3 * 6);

  CGOVertexv(shaderCGO, corner + 3 * 3);
  CGOVertexv(shaderCGO, corner + 3 * 7);

  CGOVertexv(shaderCGO, corner + 3 * 4);
  CGOVertexv(shaderCGO, corner + 3 * 5);

  CGOVertexv(shaderCGO, corner + 3 * 4);
  CGOVertexv(shaderCGO, corner + 3 * 6);

  CGOVertexv(shaderCGO, corner + 3 * 6);
  CGOVertexv(shaderCGO, corner + 3 * 7);

  CGOVertexv(shaderCGO, corner + 3 * 5);
  CGOVertexv(shaderCGO, corner + 3 * 7);

  CGOEnd(shaderCGO);

  CGOStop(shaderCGO);

  return CGOOptimizeToVBONotIndexedWithReturnedData(shaderCGO, 0, false);
}

void ObjectMap::render(RenderInfo * info)
{
  auto I = this;
  int state = info->state;
  CRay *ray = info->ray;
  auto pick = info->pick;
  const RenderPass pass = info->pass;
  ObjectMapState *ms = NULL;

  if(pass != RenderPass::Antialias)
    return;

  for(StateIterator iter(G, I->Setting.get(), state, I->State.size());
      iter.next();) {
    state = iter.state;
    if(I->State[state].Active)
      ms = &I->State[state];

    if(ms) {
      float *corner = ms->Corner;
      float tr_corner[24];
      ObjectPrepareContext(I, info);

      if(!ms->Matrix.empty()) {    /* transform the corners before drawing */
        int a;
        for(a = 0; a < 8; a++) {
          transform44d3f(ms->Matrix.data(), corner + 3 * a, tr_corner + 3 * a);
        }
        corner = tr_corner;
      }

      if((I->visRep & cRepExtentBit)) {
        if(ray) {
          const float *vc;
          float radius = ray->PixelRadius / 1.4142F;
          vc = ColorGet(G, I->Color);
          ray->color3fv(vc);
          ray->sausage3fv(corner + 3 * 0, corner + 3 * 1, radius, vc, vc);
          ray->sausage3fv(corner + 3 * 0, corner + 3 * 2, radius, vc, vc);
          ray->sausage3fv(corner + 3 * 2, corner + 3 * 3, radius, vc, vc);
          ray->sausage3fv(corner + 3 * 1, corner + 3 * 3, radius, vc, vc);
          ray->sausage3fv(corner + 3 * 0, corner + 3 * 4, radius, vc, vc);
          ray->sausage3fv(corner + 3 * 1, corner + 3 * 5, radius, vc, vc);
          ray->sausage3fv(corner + 3 * 2, corner + 3 * 6, radius, vc, vc);
          ray->sausage3fv(corner + 3 * 3, corner + 3 * 7, radius, vc, vc);
          ray->sausage3fv(corner + 3 * 4, corner + 3 * 5, radius, vc, vc);
          ray->sausage3fv(corner + 3 * 4, corner + 3 * 6, radius, vc, vc);
          ray->sausage3fv(corner + 3 * 6, corner + 3 * 7, radius, vc, vc);
          ray->sausage3fv(corner + 3 * 5, corner + 3 * 7, radius, vc, vc);
        } else if(G->HaveGUI && G->ValidContext) {
          if(pick) {
#ifndef PURE_OPENGL_ES_2
          } else if (!info->use_shaders) {
            // immediate
            ObjectUseColor(I);
            glDisable(GL_LIGHTING);
            glBegin(GL_LINES);
            glVertex3fv(corner + 3 * 0);
            glVertex3fv(corner + 3 * 1);

            glVertex3fv(corner + 3 * 0);
            glVertex3fv(corner + 3 * 2);

            glVertex3fv(corner + 3 * 2);
            glVertex3fv(corner + 3 * 3);

            glVertex3fv(corner + 3 * 1);
            glVertex3fv(corner + 3 * 3);

            glVertex3fv(corner + 3 * 0);
            glVertex3fv(corner + 3 * 4);

            glVertex3fv(corner + 3 * 1);
            glVertex3fv(corner + 3 * 5);

            glVertex3fv(corner + 3 * 2);
            glVertex3fv(corner + 3 * 6);

            glVertex3fv(corner + 3 * 3);
            glVertex3fv(corner + 3 * 7);

            glVertex3fv(corner + 3 * 4);
            glVertex3fv(corner + 3 * 5);

            glVertex3fv(corner + 3 * 4);
            glVertex3fv(corner + 3 * 6);

            glVertex3fv(corner + 3 * 6);
            glVertex3fv(corner + 3 * 7);

            glVertex3fv(corner + 3 * 5);
            glVertex3fv(corner + 3 * 7);

            glEnd();
            glEnable(GL_LIGHTING);
          } else {
#endif
            // shader
            if (!ms->shaderCGO) {
              ms->shaderCGO.reset(ObjectMapCGOGenerate(G, corner));
            }

            if (ms->shaderCGO) {
              CShaderPrg* shaderPrg = G->ShaderMgr->Enable_DefaultShader(info->pass);
              if (shaderPrg) {
                shaderPrg->SetLightingEnabled(0);

                CGORenderGL(ms->shaderCGO.get(), ColorGet(G, I->Color),
                    NULL, NULL, info, NULL);
                shaderPrg->Disable();
              }
            }
          }
        }
      }

      if((I->visRep & cRepDotBit)) {
        if(!ms->have_range) {
          double sum = 0.0, sumsq = 0.0;
          CField *data = ms->Field->data.get();
          int cnt = data->dim[0] * data->dim[1] * data->dim[2];
          float *raw_data = (float *) data->data.data();
          int a;
          for(a = 0; a < cnt; a++) {
            double f_val = *(raw_data++);
            sum += f_val;
            sumsq += (f_val * f_val);
          }
          if(cnt) {
            float mean, stdev;
            mean = (float) (sum / cnt);
            stdev = (float) sqrt1d((sumsq - (sum * sum / cnt)) / (cnt));
            ms->high_cutoff = mean + stdev;
            ms->low_cutoff = mean - stdev;
            ms->have_range = true;
          }
        }
        if(ms->have_range && SettingGet_b(G, NULL, I->Setting.get(), cSetting_dot_normals)) {
          IsofieldComputeGradients(G, ms->Field.get());
        }
        if(ms->have_range) {
          int a;
          CField *data = ms->Field->data.get();
          int cnt = data->dim[0] * data->dim[1] * data->dim[2];
          CField *points = ms->Field->points.get();
          CField *gradients = NULL;

          if(SettingGet_b(G, NULL, I->Setting.get(), cSetting_dot_normals)) {
            gradients = ms->Field->gradients.get();
          }
          if(data && points) {
            float *raw_data = (float *) data->data.data();
            float raw_point[3], *raw_point_ptr = (float *) points->data.data();

#define RAW_POINT_TRANSFORM(ptr, v3f) { \
  if(!ms->Matrix.empty()) \
    transform44d3f(ms->Matrix.data(), ptr, v3f); \
  else \
    copy3f(ptr, v3f); \
  ptr += 3; \
}

            float *raw_gradient = NULL;
            float high_cut = ms->high_cutoff, low_cut = ms->low_cutoff;
            float width =
              SettingGet_f(G, NULL, I->Setting.get(), cSetting_dot_width);

            if(ray) {
              float radius = ray->PixelRadius * width / 1.4142F;
              float vc[3];
              int color = I->Color;
              int ramped = ColorCheckRamped(G, I->Color);

              {
                const float *tmp = ColorGet(G, I->Color);
                copy3f(tmp, vc);
              }

              for(a = 0; a < cnt; a++) {
                float f_val = *(raw_data++);
                RAW_POINT_TRANSFORM(raw_point_ptr, raw_point);
                if((f_val >= high_cut) || (f_val <= low_cut)) {
                  if(ramped) {
                    ColorGetRamped(G, color, raw_point, vc, state);
                    ray->color3fv(vc);
                  }
                  ray->sphere3fv(raw_point, radius);
                }
              }
            } else if(G->HaveGUI && G->ValidContext) {
              if(pick) {
              } else if (ALWAYS_IMMEDIATE_OR(!info->use_shaders)) {
                if(gradients) {
                  raw_gradient = (float *) gradients->data.data();
                } else {
                  glDisable(GL_LIGHTING);
                }
                {
                  int ramped = ColorCheckRamped(G, I->Color);
                  float vc[3];
                  int color = I->Color;
                  float gt[3];

                  glPointSize(width);
                  glDisable(GL_POINT_SMOOTH);
                  glBegin(GL_POINTS);
                  ObjectUseColor(I);
                  for(a = 0; a < cnt; a++) {
                    float f_val = *(raw_data++);
                    RAW_POINT_TRANSFORM(raw_point_ptr, raw_point);
                    if(f_val >= high_cut) {
                      if(raw_gradient) {
                        normalize23f(raw_gradient, gt);
                        invert3f(gt);
                        glNormal3fv(gt);
                      }
                      if(ramped) {
                        ColorGetRamped(G, color, raw_point, vc, state);
                        glColor3fv(vc);
                      }
                      glVertex3fv(raw_point);
                    } else if(f_val <= low_cut) {
                      if(raw_gradient) {
                        normalize23f(raw_gradient, gt);
                        glNormal3fv(gt);
                      }
                      if(ramped) {
                        ColorGetRamped(G, color, raw_point, vc, state);
                        glColor3fv(vc);
                      }
                      glVertex3fv(raw_point);
                    }
                    if(raw_gradient)
                      raw_gradient += 3;
                  }
                  glEnd();
                glEnable(GL_POINT_SMOOTH);
                }
              }
            }
          }
        }
      }
    }
  }
}

ObjectMapState::ObjectMapState(PyMOLGlobals* G)
    : CObjectState(G)
{
  auto I = this;
  ObjectMapStatePurge(G, I);
  ObjectStateInit(G, I);
  I->Symmetry.reset(new CSymmetry(G));
  I->Field = NULL;
  I->Origin.clear();
  I->Dim.clear();
  I->Range.clear();
  I->Grid.clear();
  I->MapSource = cMapSourceUndefined;
  I->have_range = false;
}

int ObjectMap::getNFrame() const
{
  return State.size();
}

/*========================================================================*/
CObjectState* ObjectMap::_getObjectState(int state)
{
  if (!State[state].Active)
    return nullptr;
  return &State[state];
}

/*========================================================================*/
CSymmetry const* ObjectMap::getSymmetry(int state) const
{
  auto* ms = ObjectMapGetState(const_cast<ObjectMap*>(this), state);

  if (ms) {
    return ms->Symmetry.get();
  }

  return nullptr;
}

/*========================================================================*/
bool ObjectMap::setSymmetry(CSymmetry const& symmetry, int state)
{
  bool success = false;

  for (StateIterator iter(G, Setting.get(), state, State.size()); iter.next();) {
    auto& oms = State[iter.state];
    if (oms.Active) {
      oms.Symmetry.reset(new CSymmetry(symmetry));
      success = true;
    }
  }

  if (success) {
    ObjectMapRegeneratePoints(this);
  }

  return success;
}

/*========================================================================*/
ObjectMap::ObjectMap(PyMOLGlobals * G) : pymol::CObject(G)
{
  auto I = this;
  I->type = cObjectMap;

  I->visRep = cRepExtentBit;
}


/*========================================================================*/
ObjectMapState *ObjectMapNewStateFromDesc(PyMOLGlobals * G, ObjectMap * I,
                                          ObjectMapDesc * inp_md, int state, int quiet)
{
  int ok = true;
  float v[3];
  int a, b, c, d;
  float *fp;
  ObjectMapState *ms = NULL;
  ObjectMapDesc _md, *md;
  ms = ObjectMapStatePrime(I, state);

  md = &_md;
  *(md) = *(inp_md);

  if(I) {
    ms->Origin = std::vector<float>(3);
    ms->Range = std::vector<float>(3);
    ms->Grid = std::vector<float>(3);
    ms->MapSource = cMapSourceDesc;
  }
  switch (md->mode) {
  case cObjectMap_OrthoMinMaxGrid:     /* Orthorhombic: min, max, spacing, centered over range  */

    subtract3f(md->MaxCorner, md->MinCorner, v);
    for(a = 0; a < 3; a++) {
      if(v[a] < 0.0)
        std::swap(md->MaxCorner[a], md->MinCorner[a]);
    };
    subtract3f(md->MaxCorner, md->MinCorner, v);
    for(a = 0; a < 3; a++) {
      md->Dim[a] = (int) (v[a] / md->Grid[a]);
      if(md->Dim[a] < 1)
        md->Dim[a] = 1;
      if((md->Dim[a] * md->Grid[a]) < v[a])
        md->Dim[a]++;
    }

    PRINTFB(I->G, FB_ObjectMap, FB_Blather)
      " ObjectMap: Dim %d %d %d\n", md->Dim[0], md->Dim[1], md->Dim[2]
      ENDFB(I->G);

    average3f(md->MaxCorner, md->MinCorner, v);
    for(a = 0; a < 3; a++) {
      md->MinCorner[a] = v[a] - 0.5F * md->Dim[a] * md->Grid[a];
    }

    if(Feedback(I->G, FB_ObjectMap, FB_Blather)) {
      dump3f(md->MinCorner, " ObjectMap: MinCorner:");
      dump3f(md->MaxCorner, " ObjectMap: MaxCorner:");
      dump3f(md->Grid, " ObjectMap: Grid:");
    }

    /* now populate the map data structure */

    copy3f(md->MinCorner, ms->Origin.data());
    copy3f(md->Grid, ms->Grid.data());
    for(a = 0; a < 3; a++)
      ms->Range[a] = md->Grid[a] * (md->Dim[a] - 1);

    /* these maps start at zero */
    for(a = 0; a < 3; a++) {
      ms->Min[a] = 0;
      ms->Max[a] = md->Dim[a] - 1;
      ms->Div[a] = md->Dim[a] - 1;
    }

    /* define corners */

    for(a = 0; a < 8; a++)
      copy3f(ms->Origin.data(), ms->Corner + 3 * a);

    d = 0;
    for(c = 0; c < 2; c++) {
      {
        v[2] = (c ? ms->Range[2] : 0.0F);
        for(b = 0; b < 2; b++) {
          v[1] = (b ? ms->Range[1] : 0.0F);
          for(a = 0; a < 2; a++) {
            v[0] = (a ? ms->Range[0] : 0.0F);
            add3f(v, ms->Corner + 3 * d, ms->Corner + 3 * d);
            d++;
          }
        }
      }
    }
    for(a = 0; a < 3; a++)
      ms->FDim[a] = ms->Max[a] - ms->Min[a] + 1;
    ms->FDim[3] = 3;

    ms->Field.reset(new Isofield(I->G, ms->FDim));
    if(!ms->Field)
      ok = false;
    else {
      for(a = 0; a < md->Dim[0]; a++) {
        v[0] = md->MinCorner[0] + a * md->Grid[0];
        for(b = 0; b < md->Dim[1]; b++) {
          v[1] = md->MinCorner[1] + b * md->Grid[1];
          for(c = 0; c < md->Dim[2]; c++) {
            v[2] = md->MinCorner[2] + c * md->Grid[2];
            fp = F4Ptr(ms->Field->points, a, b, c, 0);
            copy3f(v, fp);
          }
        }
      }
    }
    break;
  default:
    ok = false;
  }
  if(ok) {
    switch (md->init_mode) {
    case 0:
      for(a = 0; a < md->Dim[0]; a++) {
        for(b = 0; b < md->Dim[1]; b++) {
          for(c = 0; c < md->Dim[2]; c++) {
            F3(ms->Field->data, a, b, c) = 0.0F;
          }
        }
      }
      break;
    case 1:
      for(a = 0; a < md->Dim[0]; a++) {
        for(b = 0; b < md->Dim[1]; b++) {
          for(c = 0; c < md->Dim[2]; c++) {
            F3(ms->Field->data, a, b, c) = 1.0F;
          }
        }
      }
      break;
    case -2:                   /* for testing... */
      for(a = 0; a < md->Dim[0]; a++) {
        for(b = 0; b < md->Dim[1]; b++) {
          for(c = 0; c < md->Dim[2]; c++) {
            F3(ms->Field->data, a, b, c) = (float) sqrt1d(a * a + b * b + c * c);
          }
        }
      }
      break;
    }
  }

  if(ok) {
    copy3f(ms->Origin.data(), ms->ExtentMin);
    copy3f(ms->Origin.data(), ms->ExtentMax);
    add3f(ms->Range.data(), ms->ExtentMax, ms->ExtentMax);
    ObjectMapUpdateExtents(I);
  }
  if(!ok) {
    ErrMessage(I->G, "ObjectMap", "Unable to create map");
    DeleteP(I);
  } else {
    if(!quiet) {
      PRINTFB(I->G, FB_ObjectMap, FB_Actions)
        " ObjectMap: Map created.\n" ENDFB(I->G);
    }
  }

  return (ms);
}


/*========================================================================*/
static float ccp4_next_value(char ** pp, int mode) {
  char * p = *pp;
  switch(mode) {
    case 0:
      *pp += 1;
      return (float) *((int8_t *) p);
    case 1:
      *pp += 2;
      return (float) *((int16_t *) p);
    case 2:
      *pp += 4;
      return *((float *) p);
  }
  printf("ERROR unsupported mode\n");
  return 0.f;
}

/* swaps n*width bytes in memory starting at p */
static void swap_endian(char * p, int n, int width) {
  char tmp, *q, *pstop = p + (n - 1) * width + 1;
  int w2 = width/2, wm1 = width - 1;
  for(; p < pstop; p += w2) {
    for(q = p + wm1; p < q; p++, q--) {
      tmp = *p; *p = *q; *q = tmp;
    }
  }
}

static bool validateCCP4LoadType(int& format)
{
  switch (format) {
  case cLoadTypeCCP4Map:
    format = cLoadTypeCCP4Str;
  case cLoadTypeCCP4Str:
  case cLoadTypeCCP4Unspecified:
  case cLoadTypeMRC:
    return true;
  case cLoadTypeCCP4UnspecifiedStr:
    format = cLoadTypeCCP4Unspecified;
    return true;
  case cLoadTypeMRCStr:
    format = cLoadTypeMRC;
    return true;
  default:
    return false;
  }
}

static int ObjectMapCCP4StrToMap(ObjectMap * I, char *CCP4Str, int bytes, int state,
                                 int quiet, int format)
{
  auto G = I->G;
  char *p;
  int *i;
  size_t bytes_per_pt;
  char *q;
  float dens;
  int a, b, c, d, e;
  float v[3], vr[3], maxd, mind;
  int ok = true;
  int little_endian = 1, map_endian;
  /* CCP4 named from their docs */
  int nc, nr, ns;
  int map_mode;
  int ncstart, nrstart, nsstart;
  int nx, ny, nz;
  float xlen, ylen, zlen, alpha, beta, gamma;
  int ispg; // space group number
  int sym_skip;
  int mapc, mapr, maps;
  int cc[3];
  int n_pts;
  double sum, sumsq;
  float mean, stdev;
  int normalize;
  ObjectMapState *ms;
  int expectation;

  if (!validateCCP4LoadType(format)) {
    ErrMessage(G, __func__, "wrong format");
    return false;
  }

  if(bytes < 256 * sizeof(int)) {
    PRINTFB(I->G, FB_ObjectMap, FB_Errors)
      " ObjectMapCCP4: Map appears to be truncated -- aborting." ENDFB(I->G);
    return (0);
  }

  /* state check */
  if(state < 0)
    state = I->State.size();
  /* alloc/init a new MapState */
  if(I->State.size() <= state) {
    VecCheckEmplace(I->State, state, I->G);
  }
  ms = &I->State[state];

  normalize = SettingGetGlobal_b(I->G, cSetting_normalize_ccp4_maps);

  p = CCP4Str;
  little_endian = *((char *) &little_endian);
  map_endian = (*p || *(p + 1)); // NOTE: this assumes 0x0 < NC < 0x10000

  if(little_endian != map_endian) {
    if(!quiet) {
      PRINTFB(I->G, FB_ObjectMap, FB_Blather)
        " ObjectMapCCP4: Map appears to be reverse endian, swapping...\n" ENDFB(I->G);
    }
    swap_endian(p, 256, sizeof(int));
  }

  i = (int *) p;
  nc = *(i++);                  /* columns */
  nr = *(i++);                  /* rows */
  ns = *(i++);                  /* sections */
  if(!quiet) {
    PRINTFB(I->G, FB_ObjectMap, FB_Blather)
      " ObjectMapCCP4: NC %d   NR %d   NS %d\n", nc, nr, ns ENDFB(I->G);
  }
  map_mode = *(i++);            /* mode */

  switch(map_mode) {
  case 0:
    bytes_per_pt = 1; // int8_t
    break;
  case 1:
    bytes_per_pt = 2; // int16_t
    break;
  case 2:
    bytes_per_pt = 4; // float
    break;
  default:
    PRINTFB(I->G, FB_ObjectMap, FB_Errors)
      "ObjectMapCCP4-ERR: Only map mode 0-2 currently supported (this map is mode %d)",
      map_mode ENDFB(I->G);
    return (0);
  }

  if(!quiet) {
    PRINTFB(I->G, FB_ObjectMap, FB_Blather)
      " ObjectMapCCP4: Map is mode %d.\n", map_mode ENDFB(I->G);
  }
  ncstart = *(i++);
  nrstart = *(i++);
  nsstart = *(i++);

  if(!quiet) {
    PRINTFB(I->G, FB_ObjectMap, FB_Blather)
      " ObjectMapCCP4: NCSTART %d   NRSTART %d   NSSTART  %d\n",
      ncstart, nrstart, nsstart ENDFB(I->G);
  }

  nx = *(i++);
  ny = *(i++);
  nz = *(i++);

  if(!quiet) {
    PRINTFB(I->G, FB_ObjectMap, FB_Blather)
      " ObjectMapCCP4: NX %d   NY %d   NZ  %d \n", nx, ny, nz ENDFB(I->G);
  }

  xlen = *(float *) (i++);
  ylen = *(float *) (i++);
  zlen = *(float *) (i++);

  if(!quiet) {
    PRINTFB(I->G, FB_ObjectMap, FB_Blather)
      " ObjectMapCCP4: X %8.3f   Y %8.3f  Z  %8.3f \n", xlen, ylen, zlen ENDFB(I->G);
  }

  alpha = *(float *) (i++);
  beta = *(float *) (i++);
  gamma = *(float *) (i++);

  if(!quiet) {
    PRINTFB(I->G, FB_ObjectMap, FB_Blather)
      " ObjectMapCCP4: alpha %8.3f   beta %8.3f  gamma %8.3f \n",
      alpha, beta, gamma ENDFB(I->G);
  }

  mapc = *(i++);
  mapr = *(i++);
  maps = *(i++);

  if(!quiet) {
    PRINTFB(I->G, FB_ObjectMap, FB_Blather)
      " ObjectMapCCP4: MAPC %d   MAPR %d  MAPS  %d \n", mapc, mapr, maps ENDFB(I->G);
  }

  // AMIN, AMAX, AMEAN
  mind = *(float *) (i++);
  maxd = *(float *) (i++);
  mean = *(float *) (i++);

  ispg = *(i++);
  sym_skip = *(i++);

  // CCP4 Format:
  // 25-37 skew transformation
  // 38-52 future use "Storage space sometimes used by specific programs (default = 0)"

  // MRC 2000 Format:
  // 25-49 "Extra, user defined storage space (default = 0)"
  // 50-52 ORIGIN

  {
    int skew = *(i++);

    if (format != cLoadTypeMRC && skew) {
      double matrix[16];
      auto i_float = (const float *) i;

      // SKWMAT          Skew matrix S (in order S11, S12, S13, S21 etc)
      copy33f44d(i_float, matrix);
      xx_matrix_invert(matrix, matrix, 4);

      // SKWTRN          Skew translation t
      matrix[3] = i_float[9];
      matrix[7] = i_float[10];
      matrix[11] = i_float[11];

      // Xo(map) = S * (Xo(atoms) - t)
      ObjectStateSetMatrix(ms, matrix);

      PRINTFB(I->G, FB_ObjectMap, FB_Details)
        " ObjectMapCCP4: Applied skew transformation\n"
        ENDFB(I->G);
    }
  }

  // XORIGIN, YORIGIN, ZORIGIN (50 - 52)
  // TODO See "Origin Conventions" in http://situs.biomachina.org/fmap.pdf
  float * mrc2000origin = (float *)(i + 49 - 25);
  if (format != cLoadTypeCCP4Str && lengthsq3f(mrc2000origin) > R_SMALL4) {
    double matrix[] = {
      1., 0., 0., mrc2000origin[0],
      0., 1., 0., mrc2000origin[1],
      0., 0., 1., mrc2000origin[2],
      0., 0., 0., 1.};

    ObjectStateSetMatrix(ms, matrix);

    if (!quiet) {
      PRINTFB(I->G, FB_ObjectMap, FB_Details)
        " ObjectMapCCP4: Applied MRC 2000 ORIGIN %.2f %.2f %.2f\n",
        mrc2000origin[0], mrc2000origin[1], mrc2000origin[2]
          ENDFB(I->G);

      if (ncstart || nrstart || nsstart) {
        PRINTFB(I->G, FB_ObjectMap, FB_Details)
          " ObjectMapCCP4: Ignoring N*START\n" ENDFB(I->G);
      }
    }

    ncstart = 0;
    nrstart = 0;
    nsstart = 0;
  }

  i += 54 - 25;
  stdev = *(float *) (i++);

  if(!quiet) {
    PRINTFB(I->G, FB_ObjectMap, FB_Blather)
      " ObjectMapCCP4: AMIN %f AMAX %f AMEAN %f ARMS %f\n", mind, maxd, mean, stdev ENDFB(I->G);
  }

  n_pts = nc * ns * nr;

  /* at least one EM map encountered lacks NZ, so we'll try to guess it */

  if((nx == ny) && (!nz)) {
    nz = nx;

    if(!quiet) {
      PRINTFB(I->G, FB_ObjectMap, FB_Warnings)
        " ObjectMapCCP4: NZ value is zero, but NX = NY, so we'll guess NZ = NX = NY.\n"
        ENDFB(I->G);
    }
  }

  expectation = sym_skip + sizeof(int) * 256 + bytes_per_pt * n_pts;

  if(!quiet) {
    PRINTFB(I->G, FB_ObjectMap, FB_Blather)
      " ObjectMapCCP4: sym_skip %d bytes %d expectation %d\n",
      sym_skip, bytes, expectation ENDFB(I->G);
  }

  if(bytes < expectation) {
    if(bytes == (expectation - sym_skip)) {
      /* accomodate bogus CCP4 map files with bad symmetry length information */
      if(!quiet) {
        PRINTFB(I->G, FB_ObjectMap, FB_Blather)
          " ObjectMapCCP4: Map has invalid symmetry length -- working around.\n"
          ENDFB(I->G);
      }

      expectation -= sym_skip;
      sym_skip = 0;

    } else {
      PRINTFB(I->G, FB_ObjectMap, FB_Errors)
        " ObjectMapCCP4: Map appears to be truncated -- aborting.\n" ENDFB(I->G);
      return (0);
    }
  }

  q = p + (sizeof(int) * 256) + sym_skip;

  if(little_endian != map_endian && bytes_per_pt > 1) {
    swap_endian(q, n_pts, bytes_per_pt);
  }

  // with normalize == 2, use mean and stdev from file header
  if(normalize == 1 && n_pts > 1) {
    c = n_pts;
    sum = 0.0;
    sumsq = 0.0;
    while(c--) {
      dens = ccp4_next_value(&q, map_mode);
      sumsq += dens * dens;
      sum += dens;
    }
    mean = (float) (sum / n_pts);
    stdev = (float) sqrt1d((sumsq - (sum * sum / n_pts)) / (n_pts - 1));
    if(stdev < 0.000001)
      stdev = 1.0;
  }

  q = p + (sizeof(int) * 256) + sym_skip;
  mapc--;                       /* convert to C indexing... */
  mapr--;
  maps--;

  ms->Div[0] = nx;
  ms->Div[1] = ny;
  ms->Div[2] = nz;

  ms->FDim[mapc] = nc;
  ms->Min[mapc] = ncstart;
  ms->Max[mapc] = nc + ncstart - 1;

  ms->FDim[mapr] = nr;
  ms->Min[mapr] = nrstart;
  ms->Max[mapr] = nr + nrstart - 1;

  ms->FDim[maps] = ns;
  ms->Min[maps] = nsstart;
  ms->Max[maps] = ns + nsstart - 1;

  if(!quiet) {
    if(Feedback(I->G, FB_ObjectMap, FB_Blather)) {
      dump3i(ms->Div, " ObjectMapCCP4: ms->Div");
      dump3i(ms->Min, " ObjectMapCCP4: ms->Min");
      dump3i(ms->Max, " ObjectMapCCP4: ms->Max");
      dump3i(ms->FDim, " ObjectMapCCP4: ms->FDim");
    }
  }

  ms->Symmetry->Crystal.setDims(xlen, ylen, zlen);
  ms->Symmetry->Crystal.setAngles(alpha, beta, gamma);

  if(ispg < n_space_group_numbers) {
    ms->Symmetry->setSpaceGroup(space_group_numbers[ispg]);
  }

  /* -- JV; for vol 
  ms->Mean = mean;
  ms->SD = stdev;
  ms->MaxValue = maxd;
  ms->MinValue = mind;
  */

  ms->FDim[3] = 3;
  if(!(ms->FDim[0] && ms->FDim[1] && ms->FDim[2]))
    ok = false;
  else {
    ms->Field.reset(new Isofield(I->G, ms->FDim));
    ms->MapSource = cMapSourceCCP4;
    ms->Field->save_points = false;

    maxd = -FLT_MAX;
    mind = FLT_MAX;

    for(cc[maps] = 0; cc[maps] < ms->FDim[maps]; cc[maps]++) {
      v[maps] = (cc[maps] + ms->Min[maps]) / ((float) ms->Div[maps]);

      for(cc[mapr] = 0; cc[mapr] < ms->FDim[mapr]; cc[mapr]++) {
        v[mapr] = (cc[mapr] + ms->Min[mapr]) / ((float) ms->Div[mapr]);

        for(cc[mapc] = 0; cc[mapc] < ms->FDim[mapc]; cc[mapc]++) {
          v[mapc] = (cc[mapc] + ms->Min[mapc]) / ((float) ms->Div[mapc]);

          dens = ccp4_next_value(&q, map_mode);

          if(normalize)
            dens = (dens - mean) / stdev;
          F3(ms->Field->data, cc[0], cc[1], cc[2]) = dens;
          if(maxd < dens)
            maxd = dens;
          if(mind > dens)
            mind = dens;
          transform33f3f(ms->Symmetry->Crystal.fracToReal(), v, vr);
          for(e = 0; e < 3; e++)
            F4(ms->Field->points, cc[0], cc[1], cc[2], e) = vr[e];
        }
      }
    }
  }
  if(ok) {
    d = 0;
    for(c = 0; c < 2; c++) {
      v[2] = (c * (ms->FDim[2] - 1) + ms->Min[2]) / ((float) ms->Div[2]);
      for(b = 0; b < 2; b++) {
        v[1] = (b * (ms->FDim[1] - 1) + ms->Min[1]) / ((float) ms->Div[1]);
        for(a = 0; a < 2; a++) {
          v[0] = (a * (ms->FDim[0] - 1) + ms->Min[0]) / ((float) ms->Div[0]);
          transform33f3f(ms->Symmetry->Crystal.fracToReal(), v, vr);
          copy3f(vr, ms->Corner + 3 * d);
          d++;
        }
      }
    }
  }

  if(!quiet) {
    PRINTFB(I->G, FB_ObjectMap, FB_Details)
      " ObjectMapCCP4: Map Size %d x %d x %d\n", ms->FDim[0], ms->FDim[1], ms->FDim[2]
      ENDFB(I->G);
  }

  if(!quiet) {
    if(n_pts > 1) {
      if(normalize) {
        PRINTFB(I->G, FB_ObjectMap, FB_Details)
          " ObjectMapCCP4: Normalizing with mean = %8.6f and stdev = %8.6f.\n",
          mean, stdev ENDFB(I->G);
      } else {

        PRINTFB(I->G, FB_ObjectMap, FB_Details)
          " ObjectMapCCP4: Map will not be normalized.\n ObjectMapCCP4: Current mean = %8.6f and stdev = %8.6f.\n",
          mean, stdev ENDFB(I->G);
      }
    }
  }

  if(ok) {
    v[2] = (ms->Min[2]) / ((float) ms->Div[2]);
    v[1] = (ms->Min[1]) / ((float) ms->Div[1]);
    v[0] = (ms->Min[0]) / ((float) ms->Div[0]);

    transform33f3f(ms->Symmetry->Crystal.fracToReal(), v, ms->ExtentMin);

    v[2] = ((ms->FDim[2] - 1) + ms->Min[2]) / ((float) ms->Div[2]);
    v[1] = ((ms->FDim[1] - 1) + ms->Min[1]) / ((float) ms->Div[1]);
    v[0] = ((ms->FDim[0] - 1) + ms->Min[0]) / ((float) ms->Div[0]);

    transform33f3f(ms->Symmetry->Crystal.fracToReal(), v, ms->ExtentMax);
  }
#ifdef _UNDEFINED
  printf("%d %d %d %d %d %d %d %d %d\n",
         ms->Div[0],
         ms->Min[0],
         ms->Max[0],
         ms->Div[1], ms->Min[1], ms->Max[1], ms->Div[2], ms->Min[2], ms->Max[2]);
  printf("Okay? %d\n", ok);
  fflush(stdout);
#endif
  if(!ok) {
    ErrMessage(I->G, "ObjectMap", "Error reading map");
  } else {
    ms->Active = true;
    ObjectMapUpdateExtents(I);
    if(!quiet) {
      PRINTFB(I->G, FB_ObjectMap, FB_Results)
        " ObjectMap: Map read.  Range: %5.3f to %5.3f\n", mind, maxd ENDFB(I->G);
    }
  }

  return (ok);
}


/*========================================================================*/
ObjectMapState * getObjectMapState(PyMOLGlobals * G, const char * name, int state) {
  auto* om = ExecutiveFindObject<ObjectMap>(G, name);
  return om ? om->getObjectMapState(state) : nullptr;
}


/*========================================================================*/
/**
 * Export to CCP4 or MRC map format.
 *
 * CCP4:
 * - can have SKWMAT/SKWTRN, but almost no other progam supports them (Maestro does, to some extent)
 * - can't have ORIGIN
 *
 * MRC:
 * - can have ORIGIN, most programs support it
 * - can't have SKWMAT/SKWTRN
 * - can't have N*START (or shouldn't, because of weird ORIGIN conventions)
 * - can't have non-orthogonal axes (or shouldn't, because others might ignore it)
 */
std::vector<char> ObjectMapStateToCCP4Str(const ObjectMapState * ms, int quiet, int format)
{
  std::vector<char> buffer;

  if (!ms || !ms->Active)
    return buffer; // empty

  auto G = ms->G;
  auto field = ms->Field->data;

  if (field->type != cFieldFloat ||
      field->base_size != 4) {
    PRINTFB(G, FB_ObjectMap, FB_Errors)
      " MapStateToCCP4-Error: Unsupported field type\n" ENDFB(G);
    return buffer; // empty
  }

  if (!validateCCP4LoadType(format)) {
    ErrMessage(G, __func__, "wrong format");
    return {};
  }

  buffer.resize(1024 + field->size(), 0);
  auto buffer_s = reinterpret_cast<char*>(&buffer.front());
  auto buffer_i = reinterpret_cast<int32_t*>(&buffer.front());
  auto buffer_f = reinterpret_cast<float*>(&buffer.front());

  buffer_i[0] = ms->FDim[2]; // NC / NX
  buffer_i[1] = ms->FDim[1]; // NR / NY
  buffer_i[2] = ms->FDim[0]; // NS / NZ
  buffer_i[3] = 2;           // MODE (2 = 32-bit reals)
  buffer_i[4] = ms->Min[2];  // NCSTART / NXSTART
  buffer_i[5] = ms->Min[1];  // NRSTART / NYSTART
  buffer_i[6] = ms->Min[0];  // NSSTART / NZSTART
  buffer_i[7] = ms->Div[0];  // NX / MX
  buffer_i[8] = ms->Div[1];  // NY / MY
  buffer_i[9] = ms->Div[2];  // NZ / MZ

  if (ms->Div[0] == 0) {
    // fallback
    buffer_i[7] = ms->FDim[0] - 1;  // NX / MX
    buffer_i[8] = ms->FDim[1] - 1;  // NY / MY
    buffer_i[9] = ms->FDim[2] - 1;  // NZ / MZ
  }

  bool cell_fallback = true;

  // Cell
  if (ms->Symmetry) {
    auto const* dim = ms->Symmetry->Crystal.dims();
    auto const* angle = ms->Symmetry->Crystal.angles();
    buffer_f[10] = dim[0];     // X length / CELL A
    buffer_f[11] = dim[1];     // Y length
    buffer_f[12] = dim[2];     // Z length
    buffer_f[13] = angle[0];   // Alpha    / CELL B
    buffer_f[14] = angle[1];   // Beta
    buffer_f[15] = angle[2];   // Gamma

    // check for 1x1x1 dummy cell
    cell_fallback = fabs(lengthsq3f(dim) - 3.f) < R_SMALL4;
  }

  if (cell_fallback) {
    // fallback: orthogonal extents box
    subtract3f(ms->ExtentMax, ms->ExtentMin, buffer_f + 10);    // CELL A
    set3f(buffer_f + 13, 90.f, 90.f, 90.f);                     // CELL B
  }

  buffer_i[16] = 3; // MAPC
  buffer_i[17] = 2; // MAPR
  buffer_i[18] = 1; // MAPS

  // dummy values
  buffer_f[19] = -5.f; // AMIN
  buffer_f[20] =  5.f; // AMAX
  buffer_f[21] =  0.f; // AMEAN

  // Space group number
  if (ms->Symmetry) {
    for (int ispg = 0; ispg < n_space_group_numbers; ++ispg) {
      if (strcmp(ms->Symmetry->spaceGroup(), space_group_numbers[ispg]) == 0) {
        buffer_i[22] = ispg; // ISPG
        break;
      }
    }
  }

  buffer_i[23] = 0; // NSYMBT

  float* b_skwmat = buffer_f + 25;
  float* b_skwtrn = buffer_f + 34;
  float* b_origin = buffer_f + 49;

  // MRC typically only supports orthogonal axes
  const float _f3_90[3]{90.f, 90.f, 90.f};
  bool mrc_possible =
      format != cLoadTypeCCP4Str && equal3f(buffer_f + 13 /* CELL B */, _f3_90);

  if (format == cLoadTypeMRC && !mrc_possible) {
    PRINTFB(G, FB_ObjectMap, FB_Warnings)
    " %s-Warning: MRC expects orthogonal axes\n", __func__ ENDFB(G);
  }

  // skew transformation
  if (!ms->Matrix.empty()) {
    double m[16];
    copy44d(ms->Matrix.data(), m);

    // Skew translation t
    set3f(buffer_f + 34, m[3], m[7], m[11]);    // SKWTRN

    // Skew matrix S (in order S11, S12, S13, S21 etc)
    m[3] = m[7] = m[11] = 0.;
    xx_matrix_invert(m, m, 4);
    copy44d33f(m, buffer_f + 25);               // SKWMAT

    if (mrc_possible) {
      mrc_possible = is_identityf(3, b_skwmat);

      if (format == cLoadTypeMRC && !mrc_possible) {
        PRINTFB(G, FB_ObjectMap, FB_Warnings)
        " %s-Warning: MRC expects orthonormal map\n", __func__ ENDFB(G);
      }
    }

    if (mrc_possible) {
      // translation only -> use origin instead of skew matrix
      copy3(b_skwtrn, b_origin);                // MRC ORIGIN
      zero3(b_skwtrn);
    } else if (format != cLoadTypeMRC) {
      buffer_i[24] = 1;                         // LSKFLG
    }
  }

  // origin (stored with skew transformation)
  if (!ms->Origin.empty() && lengthsq3f(ms->Origin.data()) > R_SMALL4) {
    if (!mrc_possible && !buffer_i[24] /* LSKFLG */) {
      identity33f(b_skwmat);
      buffer_i[24] = 1;                         // LSKFLG
    }

    if (buffer_i[24] /* LSKFLG */) {
      add3f(ms->Origin.data(), b_skwtrn, b_skwtrn);    // add to SKWTRN
    } else {
      add3f(ms->Origin.data(), b_origin, b_origin);    // add to MRC ORIGIN
    }
  }

  // Don't write ORIGIN and N*START. The convention is to ignore N*START if
  // ORIGIN != (0,0,0), likely because IMOD uses ORIGIN but ignores N*START.
  if (lengthsq3f(b_origin) > R_SMALL4 || format == cLoadTypeMRC) {
    // data origin in fractional coordinates
    float n_start_frac[3] = {
        float(buffer_i[6] /* NSSTART */) / (buffer_i[7] /* NX */),
        float(buffer_i[5] /* NRSTART */) / (buffer_i[8] /* NY */),
        float(buffer_i[4] /* NCSTART */) / (buffer_i[9] /* NZ */)};

    // data origin in Angstrom
    float n_start_real[3];
    transform33f3f(ms->Symmetry->Crystal.fracToReal(), n_start_frac,
                   n_start_real);

    // add to MRC origin
    add3f(n_start_real, b_origin, b_origin);

    // get rid of N*START
    buffer_i[4] = 0;
    buffer_i[5] = 0;
    buffer_i[6] = 0;
  }

  if (format == cLoadTypeCCP4Unspecified) {
    format = mrc_possible ? cLoadTypeMRC : cLoadTypeCCP4Str;
  }

  // Character string 'MAP ' to identify file type
  memcpy(buffer_s + 52 * 4, "MAP ", 4); // MAP

  // Machine stamp indicating the machine type which wrote file
  if (1 == *(int32_t*)"\1\0\0\0") {
    buffer_i[53] = 0x00004144;  // MACHST (little endian)
  } else {
    buffer_i[53] = 0x11110000;  // MACHST (big endian)
  }

  buffer_f[54] = 1.f;   // ARMS
  buffer_i[55] = 1;     // NLABL
  sprintf(buffer_s + 56 * 4, "PyMOL %s format=%s", _PyMOL_VERSION,
      (format == cLoadTypeMRC) ? "MRC" : "CCP4"); // 57-256 LABEL(20,10)

  // Map data array follows
  memcpy(buffer_s + 1024, field->data.data(), field->size());

  return buffer;
}


/*========================================================================*/
static int ObjectMapPHIStrToMap(ObjectMap * I, char *PHIStr, int bytes, int state,
                                int quiet)
{

  char *p;
  float dens, dens_rev;
  int a, b, c, d, e;
  float v[3], maxd, mind;
  int ok = true;
  int little_endian = 1;
  /* PHI named from their docs */
  int map_endian = 0;
  int map_dim;
  int map_bytes;

  ObjectMapState *ms;

  char cc[MAXLINELEN];
  char *rev;

  little_endian = *((char *) &little_endian);

  if(state < 0)
    state = I->State.size();
  if(I->State.size() <= state) {
    VecCheckEmplace(I->State, state, I->G);
  }
  ms = &I->State[state];

  maxd = -FLT_MAX;
  mind = FLT_MAX;
  p = PHIStr;

  if(*p)                        /* use FORMATTED IO record to determine map endiness */
    map_endian = 1;
  else
    map_endian = 0;

  p += 4;

  ParseNCopy(cc, p, 20);
  if(!quiet) {
    PRINTFB(I->G, FB_ObjectMap, FB_Details)
      " PHIStrToMap: %s\n", cc ENDFB(I->G);
  }
  p += 20;
  p += 4;

  p += 4;
  ParseNCopy(cc, p, 10);
  if(!quiet) {
    PRINTFB(I->G, FB_ObjectMap, FB_Details)
      " PHIStrToMap: %s\n", cc ENDFB(I->G);
  }
  p += 10;
  ParseNCopy(cc, p, 60);
  if(!quiet) {
    PRINTFB(I->G, FB_ObjectMap, FB_Details)
      " PHIStrToMap: %s\n", cc ENDFB(I->G);
  }
  p += 60;
  p += 4;

  rev = (char *) &dens_rev;

  if(little_endian != map_endian) {
    rev[0] = p[3];
    rev[1] = p[2];
    rev[2] = p[1];
    rev[3] = p[0];
  } else {
    rev[0] = p[0];              /* gotta go char by char because of memory alignment issues ... */
    rev[1] = p[1];
    rev[2] = p[2];
    rev[3] = p[3];
  }

  map_bytes = *((int *) rev);

  map_dim = (int) (pow((map_bytes / 4.0), 1 / 3.0) + 0.5);

  if((4 * map_dim * map_dim * map_dim) != map_bytes)    /* consistency check */
    map_dim = 65;

  if(!quiet) {
    PRINTFB(I->G, FB_ObjectMap, FB_Details)
      " PHIStrToMap: Map Size %d x %d x %d\n", map_dim, map_dim, map_dim ENDFB(I->G);
  }
  p += 4;

  ms->FDim[0] = map_dim;
  ms->FDim[1] = map_dim;
  ms->FDim[2] = map_dim;
  ms->FDim[3] = 3;

  ms->Div[0] = (map_dim - 1) / 2;
  ms->Div[1] = (map_dim - 1) / 2;
  ms->Div[2] = (map_dim - 1) / 2;
  ms->Min[0] = -ms->Div[0];
  ms->Min[1] = -ms->Div[1];
  ms->Min[2] = -ms->Div[2];
  ms->Max[0] = ms->Div[0];
  ms->Max[1] = ms->Div[1];
  ms->Max[2] = ms->Div[2];

  ms->Field.reset(new Isofield(I->G, ms->FDim));
  ms->MapSource = cMapSourceGeneralPurpose;
  ms->Field->save_points = false;

  for(c = 0; c < ms->FDim[2]; c++) {    /* z y x ordering into c b a  so that x = a, etc. */
    for(b = 0; b < ms->FDim[1]; b++) {
      for(a = 0; a < ms->FDim[0]; a++) {

        if(little_endian != map_endian) {
          rev[0] = p[3];
          rev[1] = p[2];
          rev[2] = p[1];
          rev[3] = p[0];
        } else {
          rev[0] = p[0];        /* gotta go char by char because of memory alignment issues ... */
          rev[1] = p[1];
          rev[2] = p[2];
          rev[3] = p[3];
        }
        dens = *((float *) rev);
        F3(ms->Field->data, a, b, c) = dens;
        if(maxd < dens)
          maxd = dens;
        if(mind > dens)
          mind = dens;
        p += 4;
      }
    }
  }
  p += 4;

  p += 4;
  ParseNCopy(cc, p, 16);
  if(!quiet) {
    PRINTFB(I->G, FB_ObjectMap, FB_Details)
      " PHIStrToMap: %s\n", cc ENDFB(I->G);
  }
  p += 16;
  p += 4;

  ms->Grid = std::vector<float>(3);
  p += 4;
  if(little_endian != map_endian) {
    rev[0] = p[3];
    rev[1] = p[2];
    rev[2] = p[1];
    rev[3] = p[0];
  } else {
    rev[0] = p[0];              /* gotta go char by char because of memory alignment issues ... */
    rev[1] = p[1];
    rev[2] = p[2];
    rev[3] = p[3];
  }
  ms->Grid[0] = 1.0F / (*((float *) rev));
  ms->Grid[1] = ms->Grid[0];
  ms->Grid[2] = ms->Grid[0];
  p += 4;

  ms->Origin = std::vector<float>(3);
  if(little_endian != map_endian) {
    rev[0] = p[3];
    rev[1] = p[2];
    rev[2] = p[1];
    rev[3] = p[0];
  } else {
    rev[0] = p[0];              /* gotta go char by char because of memory alignment issues ... */
    rev[1] = p[1];
    rev[2] = p[2];
    rev[3] = p[3];;
  }
  ms->Origin[0] = *((float *) rev);
  p += 4;

  if(little_endian != map_endian) {
    rev[0] = p[3];
    rev[1] = p[2];
    rev[2] = p[1];
    rev[3] = p[0];
  } else {
    rev[0] = p[0];              /* gotta go char by char because of memory alignment issues ... */
    rev[1] = p[1];
    rev[2] = p[2];
    rev[3] = p[3];;
  }
  ms->Origin[1] = *((float *) rev);
  p += 4;
  if(little_endian != map_endian) {
    rev[0] = p[3];
    rev[1] = p[2];
    rev[2] = p[1];
    rev[3] = p[0];
  } else {
    rev[0] = p[0];              /* gotta go char by char because of memory alignment issues ... */
    rev[1] = p[1];
    rev[2] = p[2];
    rev[3] = p[3];;
  }
  ms->Origin[2] = *((float *) rev);
  p += 4;

  p += 4;

  if(ok) {
    for(e = 0; e < 3; e++) {
      ms->ExtentMin[e] = ms->Origin[e] + ms->Grid[e] * ms->Min[e];
      ms->ExtentMax[e] = ms->Origin[e] + ms->Grid[e] * ms->Max[e];
    }
  }

  for(c = 0; c < ms->FDim[2]; c++) {
    v[2] = ms->Origin[2] + ms->Grid[2] * (c + ms->Min[2]);
    for(b = 0; b < ms->FDim[1]; b++) {
      v[1] = ms->Origin[1] + ms->Grid[1] * (b + ms->Min[1]);
      for(a = 0; a < ms->FDim[0]; a++) {
        v[0] = ms->Origin[0] + ms->Grid[0] * (a + ms->Min[0]);
        for(e = 0; e < 3; e++) {
          F4(ms->Field->points, a, b, c, e) = v[e];
        }
      }
    }
  }

  d = 0;
  for(c = 0; c < ms->FDim[2]; c += ms->FDim[2] - 1) {
    v[2] = ms->Origin[2] + ms->Grid[2] * (c + ms->Min[2]);

    for(b = 0; b < ms->FDim[1]; b += ms->FDim[1] - 1) {
      v[1] = ms->Origin[1] + ms->Grid[1] * (b + ms->Min[1]);

      for(a = 0; a < ms->FDim[0]; a += ms->FDim[0] - 1) {
        v[0] = ms->Origin[0] + ms->Grid[0] * (a + ms->Min[0]);
        copy3f(v, ms->Corner + 3 * d);
        d++;
      }
    }
  }

  /* interpolation test code 
     { 
     float test[3];
     float result;
     float cmp;

     for(c=0;c<ms->FDim[0]-1;c++) {
     for(b=0;b<ms->FDim[1]-1;b++) {
     for(a=0;a<ms->FDim[2]-1;a++) {
     for(e=0;e<3;e++) {
     test[e] = (F4(ms->Field->points,a,b,c,e)+
     F4(ms->Field->points,a,b,c,e))/2.0;
     }
     ObjectMapStateInterpolate(ms,test,&result,1);
     cmp = (F3(ms->Field->data,a,b,c)+
     F3(ms->Field->data,a,b,c))/2.0;
     if(fabs(cmp-result)>0.001) {
     printf("%d %d %d\n",a,b,c);
     printf("%8.5f %8.5f\n",
     cmp,
     result);
     }
     }
     }
     }
     }
   */

  if(!ok) {
    ErrMessage(I->G, "ObjectMap", "Error reading map");
  } else {
    ms->Active = true;
    ObjectMapUpdateExtents(I);
    if(!quiet) {
      PRINTFB(I->G, FB_ObjectMap, FB_Results)
        " ObjectMap: Map read.  Range: %5.3f to %5.3f\n", mind, maxd ENDFB(I->G);
    }
  }
  return (ok);
}


/*========================================================================*/
static int ObjectMapXPLORStrToMap(ObjectMap * I, char *XPLORStr, int state, int quiet)
{

  char *p;
  int a, b, c, d, e;
  float v[3], vr[3], dens, maxd, mind;
  char cc[MAXLINELEN];
  int n;
  int ok = true;
  ObjectMapState *ms;

  if(state < 0)
    state = I->State.size();
  if(I->State.size() <= state) {
    VecCheckEmplace(I->State, state, I->G);
  }

  ms = &I->State[state];

  maxd = -FLT_MAX;
  mind = FLT_MAX;
  p = XPLORStr;

  while(*p) {
    p = ParseNCopy(cc, p, 8);
    if(!*cc)
      p = ParseNextLine(p);
    else if(sscanf(cc, "%i", &n) == 1) {
      p = ParseWordCopy(cc, p, MAXLINELEN);
      if(strstr(cc, "!NTITLE") || (!*cc)) {
        p = ParseNextLine(p);
        while(n--) {
          p = ParseNextLine(p);
        }
      } else if(strstr(cc, "REMARKS")) {
        p = ParseNextLine(p);
      } else {
        break;
      }
    }
  }
  if(*p) {                      /* n contains first dimension */
    ms->Div[0] = n;
    if(sscanf(cc, "%i", &ms->Min[0]) != 1)
      ok = false;
    p = ParseNCopy(cc, p, 8);
    if(sscanf(cc, "%i", &ms->Max[0]) != 1)
      ok = false;
    p = ParseNCopy(cc, p, 8);
    if(sscanf(cc, "%i", &ms->Div[1]) != 1)
      ok = false;
    p = ParseNCopy(cc, p, 8);
    if(sscanf(cc, "%i", &ms->Min[1]) != 1)
      ok = false;
    p = ParseNCopy(cc, p, 8);
    if(sscanf(cc, "%i", &ms->Max[1]) != 1)
      ok = false;
    p = ParseNCopy(cc, p, 8);
    if(sscanf(cc, "%i", &ms->Div[2]) != 1)
      ok = false;
    p = ParseNCopy(cc, p, 8);
    if(sscanf(cc, "%i", &ms->Min[2]) != 1)
      ok = false;
    p = ParseNCopy(cc, p, 8);
    if(sscanf(cc, "%i", &ms->Max[2]) != 1)
      ok = false;
    p = ParseNextLine(p);
    p = ParseNCopy(cc, p, 12);

    float cellparams[3];

    if(sscanf(cc, "%f", &cellparams[0]) != 1)
      ok = false;
    p = ParseNCopy(cc, p, 12);
    if(sscanf(cc, "%f", &cellparams[1]) != 1)
      ok = false;
    p = ParseNCopy(cc, p, 12);
    if(sscanf(cc, "%f", &cellparams[2]) != 1)
      ok = false;

    ms->Symmetry->Crystal.setDims(cellparams);

    p = ParseNCopy(cc, p, 12);
    if(sscanf(cc, "%f", &cellparams[0]) != 1)
      ok = false;
    p = ParseNCopy(cc, p, 12);
    if(sscanf(cc, "%f", &cellparams[1]) != 1)
      ok = false;
    p = ParseNCopy(cc, p, 12);
    if(sscanf(cc, "%f", &cellparams[2]) != 1)
      ok = false;

    ms->Symmetry->Crystal.setAngles(cellparams);

    p = ParseNextLine(p);
    p = ParseNCopy(cc, p, 3);
    if(strcmp(cc, "ZYX"))
      ok = false;
    p = ParseNextLine(p);

  } else {
    ok = false;
  }
  if(ok) {

    ms->FDim[0] = ms->Max[0] - ms->Min[0] + 1;
    ms->FDim[1] = ms->Max[1] - ms->Min[1] + 1;
    ms->FDim[2] = ms->Max[2] - ms->Min[2] + 1;
    ms->FDim[3] = 3;
    if(!(ms->FDim[0] && ms->FDim[1] && ms->FDim[2]))
      ok = false;
    else {
      ms->Field.reset(new Isofield(I->G, ms->FDim));
      ms->MapSource = cMapSourceCrystallographic;
      ms->Field->save_points = false;
      for(c = 0; c < ms->FDim[2]; c++) {
        v[2] = (c + ms->Min[2]) / ((float) ms->Div[2]);
        p = ParseNextLine(p);
        for(b = 0; b < ms->FDim[1]; b++) {
          v[1] = (b + ms->Min[1]) / ((float) ms->Div[1]);
          for(a = 0; a < ms->FDim[0]; a++) {
            v[0] = (a + ms->Min[0]) / ((float) ms->Div[0]);
            p = ParseNCopy(cc, p, 12);
            if(!cc[0]) {
              p = ParseNextLine(p);
              p = ParseNCopy(cc, p, 12);
            }
            if(sscanf(cc, "%f", &dens) != 1) {
              ok = false;
            } else {
              F3(ms->Field->data, a, b, c) = dens;
              if(maxd < dens)
                maxd = dens;
              if(mind > dens)
                mind = dens;
            }
            transform33f3f(ms->Symmetry->Crystal.fracToReal(), v, vr);
            for(e = 0; e < 3; e++) {
              F4(ms->Field->points, a, b, c, e) = vr[e];
            }
          }
        }
        p = ParseNextLine(p);
      }
      if(ok) {
        d = 0;
        for(c = 0; c < ms->FDim[2]; c += (ms->FDim[2] - 1)) {
          v[2] = (c + ms->Min[2]) / ((float) ms->Div[2]);
          for(b = 0; b < ms->FDim[1]; b += (ms->FDim[1] - 1)) {
            v[1] = (b + ms->Min[1]) / ((float) ms->Div[1]);
            for(a = 0; a < ms->FDim[0]; a += (ms->FDim[0] - 1)) {
              v[0] = (a + ms->Min[0]) / ((float) ms->Div[0]);
              transform33f3f(ms->Symmetry->Crystal.fracToReal(), v, vr);
              copy3f(vr, ms->Corner + 3 * d);
              d++;
            }
          }
        }
      }
    }
  }

  if(ok) {
    v[2] = (ms->Min[2]) / ((float) ms->Div[2]);
    v[1] = (ms->Min[1]) / ((float) ms->Div[1]);
    v[0] = (ms->Min[0]) / ((float) ms->Div[0]);

    transform33f3f(ms->Symmetry->Crystal.fracToReal(), v, ms->ExtentMin);

    v[2] = ((ms->FDim[2] - 1) + ms->Min[2]) / ((float) ms->Div[2]);
    v[1] = ((ms->FDim[1] - 1) + ms->Min[1]) / ((float) ms->Div[1]);
    v[0] = ((ms->FDim[0] - 1) + ms->Min[0]) / ((float) ms->Div[0]);

    transform33f3f(ms->Symmetry->Crystal.fracToReal(), v, ms->ExtentMax);

  }
#ifdef _UNDEFINED
  printf("%d %d %d %d %d %d %d %d %d\n",
         ms->Div[0],
         ms->Min[0],
         ms->Max[0],
         ms->Div[1], ms->Min[1], ms->Max[1], ms->Div[2], ms->Min[2], ms->Max[2]);
  printf("Okay? %d\n", ok);
  fflush(stdout);
#endif
  if(!ok) {
    ErrMessage(I->G, "ObjectMap", "Error reading map");
  } else {
    ms->Active = true;
    ObjectMapUpdateExtents(I);
    if(!quiet) {
      PRINTFB(I->G, FB_ObjectMap, FB_Results)
        " ObjectMap: Map read.  Range = %5.3f to %5.3f\n", mind, maxd ENDFB(I->G);
    }
  }

  return (ok);
}


/*========================================================================*/
static int ObjectMapFLDStrToMap(ObjectMap * I, char *PHIStr, int bytes, int state,
                                int quiet)
{
  char *p;
  float dens, dens_rev;
  int a, b, c, d, e;
  float v[3], maxd, mind;
  int ok = true;
  int little_endian = 1;
  /* PHI named from their docs */
  int map_endian = 1;
  char cc[MAXLINELEN];
  char *rev;
  int ndim = 0;
  int veclen = 0;
  int nspace = 0;
  ObjectMapState *ms;
  int got_ndim = false;
  int got_dim1 = false;
  int got_dim2 = false;
  int got_dim3 = false;
  int got_data = false;
  int got_veclen = false;
  int got_min_ext = false;
  int got_max_ext = false;
  int got_field = false;
  int got_nspace = false;

  little_endian = *((char *) &little_endian);
  map_endian = little_endian;

  if(state < 0)
    state = I->State.size();
  if(I->State.size() <= state) {
    VecCheckEmplace(I->State, state, I->G);
  }
  ms = &I->State[state];

  maxd = -FLT_MAX;
  mind = FLT_MAX;

  p = PHIStr;

  while((*p) && (!(got_ndim &&
                   got_dim1 &&
                   got_dim2 &&
                   got_dim3 &&
                   got_veclen &&
                   got_min_ext && got_max_ext && got_field && got_nspace && got_data))) {

    if(!got_ndim) {
      ParseWordCopy(cc, p, 4);
      if(!strcmp(cc, "ndim")) {
        p = ParseSkipEquals(p);
        p = ParseWordCopy(cc, p, 50);
        if(sscanf(cc, "%d", &ndim) == 1)
          got_ndim = true;
      }
    }

    if(!got_dim1) {
      ParseWordCopy(cc, p, 4);
      if(!strcmp(cc, "dim1")) {
        p = ParseSkipEquals(p);
        p = ParseWordCopy(cc, p, 50);
        if(sscanf(cc, "%d", &ms->FDim[0]) == 1)
          got_dim1 = true;
      }
    }

    if(!got_dim2) {
      ParseWordCopy(cc, p, 4);
      if(!strcmp(cc, "dim2")) {
        p = ParseSkipEquals(p);
        p = ParseWordCopy(cc, p, 50);
        if(sscanf(cc, "%d", &ms->FDim[1]) == 1)
          got_dim2 = true;
      }
    }

    if(!got_dim3) {
      ParseWordCopy(cc, p, 4);
      if(!strcmp(cc, "dim3")) {
        p = ParseSkipEquals(p);
        p = ParseWordCopy(cc, p, 50);
        if(sscanf(cc, "%d", &ms->FDim[2]) == 1)
          got_dim3 = true;
      }
    }

    if(!got_nspace) {
      ParseWordCopy(cc, p, 6);
      if(!strcmp(cc, "nspace")) {
        p = ParseSkipEquals(p);
        p = ParseWordCopy(cc, p, 50);
        if(sscanf(cc, "%d", &nspace) == 1)
          got_nspace = true;
      }
    }

    if(!got_veclen) {
      ParseWordCopy(cc, p, 6);
      if(!strcmp(cc, "veclen")) {
        p = ParseSkipEquals(p);
        p = ParseWordCopy(cc, p, 50);
        if(sscanf(cc, "%d", &veclen) == 1)
          got_veclen = true;
      }
    }

    if(!got_data) {
      ParseWordCopy(cc, p, 4);
      if(!strcmp(cc, "data")) {
        p = ParseSkipEquals(p);
        p = ParseWordCopy(cc, p, 50);
        if(!strcmp(cc, "float"))
          got_data = true;
      }
    }

    if(!got_min_ext) {
      ParseWordCopy(cc, p, 7);
      if(!strcmp(cc, "min_ext")) {
        p = ParseSkipEquals(p);
        p = ParseWordCopy(cc, p, 50);
        if(sscanf(cc, "%f", &ms->ExtentMin[0]) == 1) {
          p = ParseWordCopy(cc, p, 50);
          if(sscanf(cc, "%f", &ms->ExtentMin[1]) == 1) {
            p = ParseWordCopy(cc, p, 50);
            if(sscanf(cc, "%f", &ms->ExtentMin[2]) == 1) {
              got_min_ext = true;
            }
          }
        }
      }
    }

    if(!got_max_ext) {
      ParseWordCopy(cc, p, 7);
      if(!strcmp(cc, "max_ext")) {
        p = ParseSkipEquals(p);
        p = ParseWordCopy(cc, p, 50);
        if(sscanf(cc, "%f", &ms->ExtentMax[0]) == 1) {
          p = ParseWordCopy(cc, p, 50);
          if(sscanf(cc, "%f", &ms->ExtentMax[1]) == 1) {
            p = ParseWordCopy(cc, p, 50);
            if(sscanf(cc, "%f", &ms->ExtentMax[2]) == 1) {
              got_max_ext = true;
            }
          }
        }
      }
    }

    if(!got_field) {
      ParseWordCopy(cc, p, 5);
      if(!strcmp(cc, "field")) {
        p = ParseSkipEquals(p);
        p = ParseWordCopy(cc, p, 50);
        if(!strcmp(cc, "uniform"))
          got_field = true;
      }
    }
    p = ParseNextLine(p);
  }

  if(got_ndim &&
     got_dim1 &&
     got_dim2 &&
     got_dim3 &&
     got_veclen && got_min_ext && got_max_ext && got_field && got_nspace && got_data) {

    int pass = 0;

    ms->Origin = std::vector<float>(3);
    ms->Range = std::vector<float>(3);
    ms->Grid = std::vector<float>(3);

    copy3f(ms->ExtentMin, ms->Origin.data());
    subtract3f(ms->ExtentMax, ms->ExtentMin, ms->Range.data());
    ms->FDim[3] = 3;

    PRINTFB(I->G, FB_ObjectMap, FB_Details)
      " FLDStrToMap: Map Size %d x %d x %d\n", ms->FDim[0], ms->FDim[1], ms->FDim[2]
      ENDFB(I->G);

    for(a = 0; a < 3; a++) {
      ms->Min[a] = 0;
      ms->Max[a] = ms->FDim[a] - 1;
      ms->Div[a] = ms->FDim[a] - 1;

      if(ms->FDim[a])
        ms->Grid[a] = ms->Range[a] / (ms->Max[a]);
      else
        ms->Grid[a] = 0.0F;
    }

    ms->Field.reset(new Isofield(I->G, ms->FDim));
    ms->MapSource = cMapSourceFLD;
    ms->Field->save_points = false;

    while(1) {
      maxd = -FLT_MAX;
      mind = FLT_MAX;
      p = PHIStr;

      while(*p) {               /* ^L^L sentinel */
        if((p[0] == 12) && (p[1] == 12)) {
          p += 2;
          break;
        }
        p++;
      }

      rev = (char *) &dens_rev;
      for(c = 0; c < ms->FDim[2]; c++) {        /* z y x ordering into c b a  so that x = a, etc. */
        for(b = 0; b < ms->FDim[1]; b++) {
          for(a = 0; a < ms->FDim[0]; a++) {

            if(little_endian != map_endian) {
              rev[0] = p[3];
              rev[1] = p[2];
              rev[2] = p[1];
              rev[3] = p[0];
            } else {
              rev[0] = p[0];    /* gotta go char by char because of memory alignment issues ... */
              rev[1] = p[1];
              rev[2] = p[2];
              rev[3] = p[3];
            }
            dens = *((float *) rev);
            F3(ms->Field->data, a, b, c) = dens;
            if(maxd < dens)
              maxd = dens;
            if(mind > dens)
              mind = dens;
            p += 4;
          }
        }
      }

      /* There's no way to determine the original handedness of input
         field files.  So instead, we simplymake an educated guess about
         whether we're byte-swapped based on the range of the density
         values obtained. */

      if(((maxd / FLT_MAX) > 0.1F) && ((mind / (-FLT_MAX)) > 0.1F)) {
        if(pass == 0) {
          map_endian = (!map_endian);   /* okay, try again swapped */
        } else if(pass == 1) {
          /* didn't help, so resort to original order */
          map_endian = (!map_endian);
        } else {
          break;
        }
      } else {
        break;
      }
      pass++;
    }

    for(c = 0; c < ms->FDim[2]; c++) {
      v[2] = ms->Origin[2] + ms->Grid[2] * (c + ms->Min[2]);
      for(b = 0; b < ms->FDim[1]; b++) {
        v[1] = ms->Origin[1] + ms->Grid[1] * (b + ms->Min[1]);
        for(a = 0; a < ms->FDim[0]; a++) {
          v[0] = ms->Origin[0] + ms->Grid[0] * (a + ms->Min[0]);
          for(e = 0; e < 3; e++) {
            F4(ms->Field->points, a, b, c, e) = v[e];
          }
        }
      }
    }

    d = 0;
    for(c = 0; c < ms->FDim[2]; c += ms->FDim[2] - 1) {
      v[2] = ms->Origin[2] + ms->Grid[2] * (c + ms->Min[2]);
      if(!c)
        v[2] = ms->ExtentMin[2];
      else
        v[2] = ms->ExtentMax[2];

      for(b = 0; b < ms->FDim[1]; b += ms->FDim[1] - 1) {
        v[1] = ms->Origin[1] + ms->Grid[1] * (b + ms->Min[1]);

        if(!b)
          v[1] = ms->ExtentMin[1];
        else
          v[1] = ms->ExtentMax[1];

        for(a = 0; a < ms->FDim[0]; a += ms->FDim[0] - 1) {
          v[0] = ms->Origin[0] + ms->Grid[0] * (a + ms->Min[0]);

          if(!a)
            v[0] = ms->ExtentMin[0];
          else
            v[0] = ms->ExtentMax[0];

          copy3f(v, ms->Corner + 3 * d);
          d++;
        }
      }
    }

    ms->Active = true;
    ObjectMapUpdateExtents(I);
    if(!quiet) {
      PRINTFB(I->G, FB_ObjectMap, FB_Results)
        " ObjectMap: Map read.  Range: %5.3f to %5.3f\n", mind, maxd ENDFB(I->G);
    }
  } else {
    PRINTFB(I->G, FB_ObjectMap, FB_Errors)
      " Error: unable to read FLD file.\n" ENDFB(I->G);
    /*  printf(" got_ndim %d\n",got_ndim);
       printf(" got_dim1 %d\n",got_dim1);
       printf(" got_dim2 %d\n",got_dim2);
       printf(" got_dim3 %d\n",got_dim3);
       printf(" got_veclen %d\n",got_veclen);
       printf(" got_min_ext %d\n",got_min_ext);
       printf(" got_max_ext %d\n",got_max_ext);
       printf(" got_field %d\n",got_field);
       printf(" got_nspace %d\n",got_nspace);
       printf(" got_data %d\n",got_data);
     */
  }

  return (ok);

}


/*========================================================================*/
static int ObjectMapBRIXStrToMap(ObjectMap * I, char *BRIXStr, int bytes, int state,
                                 int quiet)
{
  char *p, *pp;
  float dens;
  int a, b, c, d, e;
  float v[3], vr[3], maxd, mind;
  int ok = true;
  char cc[MAXLINELEN];
  ObjectMapState *ms;
  int got_origin = false;
  int got_extent = false;
  int got_grid = false;
  int got_cell = false;
  int got_prod = false;
  int got_plus = false;
  int got_sigma = false;
  int got_smiley = false;
  float prod = 1.0F;
  float plus = 0.0F;
  float sigma = 0.0F;
  float calc_sigma = 0.0F;
  float calc_mean = 0.0F;
  int normalize;
  int swap_bytes;

  normalize = SettingGetGlobal_b(I->G, cSetting_normalize_o_maps);
  swap_bytes = SettingGetGlobal_b(I->G, cSetting_swap_dsn6_bytes);
  if(state < 0)
    state = I->State.size();
  if(I->State.size() <= state) {
    VecCheckEmplace(I->State, state, I->G);
  }
  ms = &I->State[state];

  maxd = -FLT_MAX;
  mind = FLT_MAX;

  p = BRIXStr;

  ParseNCopy(cc, p, 3);
  got_smiley = (strcmp(cc, ":-)") == 0);

  if(got_smiley) {

    /* ASCII BRIX format header */

    while((*p) && (!(got_origin &&
                     got_extent &&
                     got_grid && got_cell && got_prod && got_plus && got_sigma))) {

      if(!got_origin) {
        pp = ParseWordCopy(cc, p, 6);
        if(WordMatch(I->G, "origin", cc, true) < 0) {
          p = ParseWordCopy(cc, pp, 50);
          if(sscanf(cc, "%d", &ms->Min[0]) == 1) {
            p = ParseWordCopy(cc, p, 50);
            if(sscanf(cc, "%d", &ms->Min[1]) == 1) {
              p = ParseWordCopy(cc, p, 50);
              if(sscanf(cc, "%d", &ms->Min[2]) == 1) {
                got_origin = true;
              }
            }
          }
        }
      }

      if(!got_extent) {
        pp = ParseWordCopy(cc, p, 6);
        if(WordMatch(I->G, "extent", cc, true) < 0) {
          p = ParseWordCopy(cc, pp, 50);
          if(sscanf(cc, "%d", &ms->Max[0]) == 1) {
            p = ParseWordCopy(cc, p, 50);
            if(sscanf(cc, "%d", &ms->Max[1]) == 1) {
              p = ParseWordCopy(cc, p, 50);
              if(sscanf(cc, "%d", &ms->Max[2]) == 1) {
                got_extent = true;
                ms->Max[0] += ms->Min[0] - 1;
                ms->Max[1] += ms->Min[1] - 1;
                ms->Max[2] += ms->Min[2] - 1;
              }
            }
          }
        }
      }

      if(!got_grid) {
        pp = ParseWordCopy(cc, p, 4);
        if(WordMatch(I->G, "grid", cc, true) < 0) {
          p = ParseWordCopy(cc, pp, 50);
          if(sscanf(cc, "%d", &ms->Div[0]) == 1) {
            p = ParseWordCopy(cc, p, 50);
            if(sscanf(cc, "%d", &ms->Div[1]) == 1) {
              p = ParseWordCopy(cc, p, 50);
              if(sscanf(cc, "%d", &ms->Div[2]) == 1) {
                got_grid = true;
              }
            }
          }
        }
      }

      if(!got_cell) {
        pp = ParseWordCopy(cc, p, 4);
        if(WordMatch(I->G, "cell", cc, true) < 0) {
          p = pp;
          got_cell = true;
          float cellparams[6];

          for (int i = 0; i != 6; ++i) {
            p = ParseWordCopy(cc, p, 50);
            if (sscanf(cc, "%f", cellparams + i) != 1) {
              got_cell = false;
              break;
            }
          }

          if (got_cell) {
            ms->Symmetry->Crystal.setDims(cellparams);
            ms->Symmetry->Crystal.setAngles(cellparams + 3);
          }
        }
      }

      if(!got_plus) {
        pp = ParseWordCopy(cc, p, 4);
        if(WordMatch(I->G, "plus", cc, true) < 0) {
          p = ParseWordCopy(cc, pp, 50);
          if(sscanf(cc, "%f", &plus) == 1) {
            got_plus = true;
          }
        }
      }

      if(!got_prod) {
        pp = ParseWordCopy(cc, p, 4);
        if(WordMatch(I->G, "prod", cc, true) < 0) {
          p = ParseWordCopy(cc, pp, 50);
          if(sscanf(cc, "%f", &prod) == 1) {
            got_prod = true;
          }
        }
      }

      if(!got_sigma) {
        pp = ParseWordCopy(cc, p, 5);
        if(WordMatch(I->G, "sigma", cc, true) < 0) {
          p = ParseWordCopy(cc, pp, 50);
          if(sscanf(cc, "%f", &sigma) == 1) {
            got_sigma = true;
          }
        }
      }

      p++;
    }
  } else {
    /* Binary DSN6 format */

    /* first, determine whether or not we need to byte-swap the header... */

    int passed_endian_check = false;
    short int *shint_ptr;
    float scale_factor;
    char tmp_char;
    int a;
    int start_swap_at = 0;
    int stop_swap_at = bytes;

    shint_ptr = (short int *) (p + 18 * 2);
    if(*shint_ptr == 100) {
      passed_endian_check = true;
      start_swap_at = 512;      /* don't byte-swap header */
    }

    if(!swap_bytes) {
      stop_swap_at = 512;
    }
    /* for DSN6, always swap the bricks */

    p = BRIXStr + start_swap_at;
    for(a = start_swap_at; a < stop_swap_at; a += 2) {
      tmp_char = *p;
      *p = *(p + 1);
      *(p + 1) = tmp_char;
      p += 2;
    }

    p = BRIXStr;

    if(*shint_ptr == 100) {
      passed_endian_check = true;
    }

    if(!passed_endian_check) {
      PRINTFB(I->G, FB_ObjectMap, FB_Errors)
        " Error: This looks like a DSN6 map file, but I can't match endianness.\n"
        ENDFB(I->G);
    } else {
      shint_ptr = (short int *) p;

      ms->Min[0] = *(shint_ptr++);
      ms->Min[1] = *(shint_ptr++);
      ms->Min[2] = *(shint_ptr++);
      got_origin = true;

      ms->Max[0] = *(shint_ptr++);
      ms->Max[1] = *(shint_ptr++);
      ms->Max[2] = *(shint_ptr++);
      got_extent = true;
      ms->Max[0] += ms->Min[0] - 1;
      ms->Max[1] += ms->Min[1] - 1;
      ms->Max[2] += ms->Min[2] - 1;

      ms->Div[0] = *(shint_ptr++);
      ms->Div[1] = *(shint_ptr++);
      ms->Div[2] = *(shint_ptr++);
      got_grid = true;

      float cellparams[6] = {
          float(*(shint_ptr++)),
          float(*(shint_ptr++)),
          float(*(shint_ptr++)),
          float(*(shint_ptr++)),
          float(*(shint_ptr++)),
          float(*(shint_ptr++)),
      };
      got_cell = true;

      prod = (float) (*(shint_ptr++)) / 100.0F;
      got_prod = true;

      plus = (float) (*(shint_ptr++));
      got_plus = true;

      scale_factor = (float) (*(shint_ptr++));

      for(a = 0; a < 6; ++a) {
        cellparams[a] /= scale_factor;
      }
      ms->Symmetry->Crystal.setDims(cellparams);
      ms->Symmetry->Crystal.setAngles(cellparams + 3);
    }
  }

  if(got_origin && got_extent && got_grid && got_cell && got_plus && got_prod) {

    PRINTFB(I->G, FB_ObjectMap, FB_Blather)
      " BRIXStrToMap: Prod = %8.3f, Plus = %8.3f\n", prod, plus ENDFB(I->G);

    ms->FDim[0] = ms->Max[0] - ms->Min[0] + 1;
    ms->FDim[1] = ms->Max[1] - ms->Min[1] + 1;
    ms->FDim[2] = ms->Max[2] - ms->Min[2] + 1;
    ms->FDim[3] = 3;
    if(!(ms->FDim[0] && ms->FDim[1] && ms->FDim[2]))
      ok = false;
    else {
      ms->Field.reset(new Isofield(I->G, ms->FDim));
      ms->MapSource = cMapSourceBRIX;
      ms->Field->save_points = false;

      {
        int block_size = 8;
        int xblocks = ((ms->FDim[0] - 1) / block_size) + 1;
        int yblocks = ((ms->FDim[1] - 1) / block_size) + 1;
        int zblocks = ((ms->FDim[2] - 1) / block_size) + 1;
        int cc, bb, aa, ia, ib, ic, xa, xb, xc;
        double sum = 0.0;
        double sumsq = 0.0;
        int n_pts = 0;

        p = BRIXStr + 512;
        for(cc = 0; cc < zblocks; cc++) {
          for(bb = 0; bb < yblocks; bb++) {
            for(aa = 0; aa < xblocks; aa++) {

              for(c = 0; c < block_size; c++) {
                xc = c + cc * block_size;
                ic = xc + ms->Min[2];
                v[2] = ic / ((float) ms->Div[2]);

                for(b = 0; b < block_size; b++) {
                  xb = b + bb * block_size;
                  ib = xb + ms->Min[1];
                  v[1] = ib / ((float) ms->Div[1]);

                  for(a = 0; a < block_size; a++) {
                    xa = a + aa * block_size;
                    ia = xa + ms->Min[0];
                    v[0] = ia / ((float) ms->Div[0]);

                    dens = (((float) (*((unsigned char *) (p++)))) - plus) / prod;
                    if((ia <= ms->Max[0]) && (ib <= ms->Max[1]) && (ic <= ms->Max[2])) {
                      F3(ms->Field->data, xa, xb, xc) = dens;
                      sumsq += dens * dens;
                      sum += dens;
                      n_pts++;
                      if(maxd < dens)
                        maxd = dens;
                      if(mind > dens)
                        mind = dens;
                      transform33f3f(ms->Symmetry->Crystal.fracToReal(), v, vr);
                      for(e = 0; e < 3; e++) {
                        F4(ms->Field->points, xa, xb, xc, e) = vr[e];
                      }

                    }
                  }
                }
              }
            }
          }
        }

        if(n_pts > 1) {
          calc_mean = (float) (sum / n_pts);
          calc_sigma = (float) sqrt1d((sumsq - (sum * sum / n_pts)) / (n_pts - 1));
        }

      }

      if(ok) {
        d = 0;
        for(c = 0; c < ms->FDim[2]; c += (ms->FDim[2] - 1)) {
          v[2] = (c + ms->Min[2]) / ((float) ms->Div[2]);
          for(b = 0; b < ms->FDim[1]; b += (ms->FDim[1] - 1)) {
            v[1] = (b + ms->Min[1]) / ((float) ms->Div[1]);
            for(a = 0; a < ms->FDim[0]; a += (ms->FDim[0] - 1)) {
              v[0] = (a + ms->Min[0]) / ((float) ms->Div[0]);
              transform33f3f(ms->Symmetry->Crystal.fracToReal(), v, vr);
              copy3f(vr, ms->Corner + 3 * d);
              d++;
            }
          }
        }
      }

      if(ok && normalize) {

        if(calc_sigma > R_SMALL8) {
          for(c = 0; c < ms->FDim[2]; c++) {
            for(b = 0; b < ms->FDim[1]; b++) {
              for(a = 0; a < ms->FDim[0]; a++) {
                F3(ms->Field->data, a, b, c) =
                  (F3(ms->Field->data, a, b, c) - calc_mean) / calc_sigma;
              }
            }
          }
        }
      }

      v[2] = (ms->Min[2]) / ((float) ms->Div[2]);
      v[1] = (ms->Min[1]) / ((float) ms->Div[1]);
      v[0] = (ms->Min[0]) / ((float) ms->Div[0]);

      transform33f3f(ms->Symmetry->Crystal.fracToReal(), v, ms->ExtentMin);

      v[2] = ((ms->FDim[2] - 1) + ms->Min[2]) / ((float) ms->Div[2]);
      v[1] = ((ms->FDim[1] - 1) + ms->Min[1]) / ((float) ms->Div[1]);
      v[0] = ((ms->FDim[0] - 1) + ms->Min[0]) / ((float) ms->Div[0]);

      transform33f3f(ms->Symmetry->Crystal.fracToReal(), v, ms->ExtentMax);

      PRINTFB(I->G, FB_ObjectMap, FB_Details)
        " BRIXStrToMap: Map Size %d x %d x %d\n", ms->FDim[0], ms->FDim[1], ms->FDim[2]
        ENDFB(I->G);

      if(got_sigma) {
        PRINTFB(I->G, FB_ObjectMap, FB_Details)
          " BRIXStrToMap: Reported Sigma = %8.3f\n", sigma ENDFB(I->G);
      }

      PRINTFB(I->G, FB_ObjectMap, FB_Details)
        " BRIXStrToMap: Range = %5.6f to %5.6f\n", mind, maxd ENDFB(I->G);

      PRINTFB(I->G, FB_ObjectMap, FB_Details)
        " BRIXStrToMap: Calculated Mean = %8.3f, Sigma = %8.3f\n", calc_mean, calc_sigma
        ENDFB(I->G);

      if(normalize) {
        PRINTFB(I->G, FB_ObjectMap, FB_Details)
          " BRIXStrToMap: Normalizing...\n" ENDFB(I->G);
      }

      ms->Active = true;
      ObjectMapUpdateExtents(I);

      if (!quiet) {
        CrystalDump(&ms->Symmetry->Crystal);
      }
    }
  } else {
    PRINTFB(I->G, FB_ObjectMap, FB_Errors)
      " Error: unable to read BRIX/DSN6 file.\n" ENDFB(I->G);
  }

  return (ok);

}


/*========================================================================*/
static int ObjectMapGRDStrToMap(ObjectMap * I, char *GRDStr, int bytes, int state,
                                int quiet)
{
  /* NOTE: binary GRD reader not yet validated */

  char *p;
  float dens;
  float *f = NULL;
  int a, b, c, d, e;
  float v[3], vr[3], maxd, mind;
  int ok = true;
  char cc[MAXLINELEN];
  ObjectMapState *ms;
  int got_cell = false;
  int got_brick = false;
  int fast_axis = 1;            /* assumes fast X */

  int got_ranges = false;
  int normalize = false;
  float mean, stdev;
  double sum = 0.0;
  double sumsq = 0.0;
  int n_pts = 0;
  int ascii = true;
  int little_endian = 1, map_endian = 0;
  union int_or_char_array {
    int block_len;
    char rev[4];
  };
  union int_or_char_array rev_union;
  rev_union.block_len = 0;

  if(state < 0)
    state = I->State.size();
  if(I->State.size() <= state) {
    VecCheckEmplace(I->State, state, I->G);
  }
  ms = &I->State[state];
  normalize = SettingGetGlobal_b(I->G, cSetting_normalize_grd_maps);
  maxd = -FLT_MAX;
  mind = FLT_MAX;

  p = GRDStr;

  if(bytes > 4) {
    if((!p[0]) || (!p[1]) || (!p[2]) || (!p[3])) {
      /* if zero appears in the first four bytes, then this is a binary map */
      ascii = false;

      little_endian = *((char *) &little_endian);
      map_endian = (*p || *(p + 1));

    }
  }
  /* print map title */

  if(ascii) {
    p = ParseNCopy(cc, p, 100);
  } else {


/* 

according to one site on the internet...GRD binary formats look like this:

        write(14) title
        write(14) scale,oldmid
        write(14) ivary,nbyte,intdat,extent,extent
        write(14) extent,xang,yang,zang,xstart
        write(14) xend,ystart,yend,zstart,zend
        write(14) intx,inty,intz

        do 100 i = 1,intz+1
            do 100 j = 1,inty+1
100            write(14)(phimap(k,j,i),k=1,intx+1)

    where: scale is the inverse grid spacing,
    oldmid are x,y,z coordinates of the grid center,
    intx,inty,intz = grid points per side - 1
    ivary = 0
    nbyte = 4
    intdat = 0
    xang,yang,zang = 90
    extent is the absolute value of the x,y, or z coordinate of the grid
    corners with the largest absolute value
    xstart,xend,ystart, yend,zstart,zend are the fractional limits (-1 to 1)
    of the molecule in each direction.
    title is a descriptive header for the molecule 

However, that doesn't make sense...in the files I see, scale and oldmid don't seem to exist!

So here's another claim which I find more credible...

character*132 toplbl !ascii header 
integer*4 ivary !0 => x index varys most rapidly
integer*4 nbyte !=4, # of bytes in data
integer*4 inddat !=0, floating point data
real*4 xang,yang,zang !=90,90,90 unit cell angles
integer*4 intx,inty,intz !=igrid-1, # of intervals/grid side
real*4 extent !maximum extent of grid
real*4 xstart,xend !beginning, end of grid sides
real*4 ystart,yend !in fractional
real*4 zstart,zend !units of extent
write(14)toplbl write(14)ivary, nbyte, intdat, extent, extent, extent, xang, yang, zang, xstart, xend, ystart, yend, zstart, zend, intx, inty, intz 
do k = 1,igrid 
   do j = 1,igrid 
      write(14)(phimap(i,j,k),i=1,igrid) 
   end do 
end d

*/

    if(!quiet) {
      PRINTFB(I->G, FB_ObjectMap, FB_Warnings)
        " ObjectMapGRD-Warning: Binary GRD reader not yet validated.\n" ENDFB(I->G);
    }

    if(little_endian != map_endian) {
      rev_union.rev[0] = p[3];
      rev_union.rev[1] = p[2];
      rev_union.rev[2] = p[1];
      rev_union.rev[3] = p[0];
    } else {
      rev_union.rev[0] = p[0];  /* gotta go char by char because of memory alignment issues ... */
      rev_union.rev[1] = p[1];
      rev_union.rev[2] = p[2];
      rev_union.rev[3] = p[3];
    }
    ParseNCopy(cc, p + 4, rev_union.block_len);
    /* now flip file to correct endianess */
    if(little_endian != map_endian) {
      char *c = p;
      int cnt = bytes >> 2;
      unsigned char c0, c1, c2, c3;
      while(cnt) {
        c0 = c[0];
        c1 = c[1];
        c2 = c[2];
        c3 = c[3];
        c[0] = c3;
        c[1] = c2;
        c[2] = c1;
        c[3] = c0;
        c += 4;
        cnt--;
      }
    }
    p += 4 + rev_union.block_len;
    f = (float *) p;
  }

  PRINTFB(I->G, FB_ObjectMap, FB_Details)
    " ObjectMap: %s\n", cc ENDFB(I->G);

  if(ascii)
    p = ParseNextLine(p);

  /* skip format -- we're reading float regardless... */
  if(ascii)
    p = ParseNextLine(p);
  else {
    {
      int block_len_check;
      int ivary;
      int nbyte;
      int intdat;

      block_len_check = *((int *) (f++));

      if(rev_union.block_len != block_len_check) {
        PRINTFB(I->G, FB_ObjectMap, FB_Warnings)
          " ObjectMapGRD-Warning: block length not matched -- not a true GRD binary?\n"
          ENDFB(I->G);
      }

      rev_union.block_len = *((int *) (f++));
      ivary = *((int *) (f++));
      nbyte = *((int *) (f++));
      intdat = *((int *) (f++));

      if((ivary) || (nbyte != 4) || (intdat)) {
        if(!quiet) {
          PRINTFB(I->G, FB_ObjectMap, FB_Warnings)
            " ObjectMapGRD-Warning: funky ivary, nbyte, intdat -- not a true GRD binary?\n"
            ENDFB(I->G);
        }
      }
    }

  }

  /* read unit cell */

  if(ok) {
    float cellparams[6];

    if(ascii) {
      p = ParseWordCopy(cc, p, 50);
      if(sscanf(cc, "%f", &cellparams[0]) == 1) {
        p = ParseWordCopy(cc, p, 50);
        if(sscanf(cc, "%f", &cellparams[1]) == 1) {
          p = ParseWordCopy(cc, p, 50);
          if(sscanf(cc, "%f", &cellparams[2]) == 1) {
            p = ParseWordCopy(cc, p, 50);
            if(sscanf(cc, "%f", &cellparams[3]) == 1) {
              p = ParseWordCopy(cc, p, 50);
              if(sscanf(cc, "%f", &cellparams[4]) == 1) {
                p = ParseWordCopy(cc, p, 50);
                if(sscanf(cc, "%f", &cellparams[5]) == 1) {
                  got_cell = true;
                }
              }
            }
          }
        }
      }
      ok = got_cell;
    } else {
      cellparams[0] = *(f++);   /* x-extent */
      cellparams[1] = *(f++);   /* y-extent */
      cellparams[2] = *(f++);   /* z-extent */
      cellparams[3] = *(f++);   /* xang */
      cellparams[4] = *(f++);   /* yang */
      cellparams[5] = *(f++);   /* zang */
      got_cell = 1;
    }

    ms->Symmetry->Crystal.setDims(cellparams);
    ms->Symmetry->Crystal.setAngles(cellparams + 3);
  }

  if(ascii)
    p = ParseNextLine(p);

  if(ascii) {
    /* read brick dimensions */

    if(ok) {
      p = ParseWordCopy(cc, p, 50);
      if(sscanf(cc, "%d", &ms->FDim[0]) == 1) {
        p = ParseWordCopy(cc, p, 50);
        if(sscanf(cc, "%d", &ms->FDim[1]) == 1) {
          p = ParseWordCopy(cc, p, 50);
          if(sscanf(cc, "%d", &ms->FDim[2]) == 1) {
            got_brick = true;
            ms->FDim[0]++;
            ms->FDim[1]++;
            ms->FDim[2]++;      /* stupid 0-based counters */
          }
        }
      }
      ok = got_brick;
    }

    p = ParseNextLine(p);

    /* read ranges */

    if(ok) {
      p = ParseWordCopy(cc, p, 50);
      if(sscanf(cc, "%d", &fast_axis) == 1) {
        p = ParseWordCopy(cc, p, 50);
        if(sscanf(cc, "%d", &ms->Min[0]) == 1) {
          p = ParseWordCopy(cc, p, 50);
          if(sscanf(cc, "%d", &ms->Div[0]) == 1) {
            p = ParseWordCopy(cc, p, 50);
            if(sscanf(cc, "%d", &ms->Min[1]) == 1) {
              p = ParseWordCopy(cc, p, 50);
              if(sscanf(cc, "%d", &ms->Div[1]) == 1) {
                p = ParseWordCopy(cc, p, 50);
                if(sscanf(cc, "%d", &ms->Min[2]) == 1) {
                  p = ParseWordCopy(cc, p, 50);
                  if(sscanf(cc, "%d", &ms->Div[2]) == 1) {
                    got_ranges = true;
                  }
                }
              }
            }
          }
        }
      }
      ok = got_ranges;
    }
  } else {
    float fmin[3];
    float fmax[3];

    fmin[0] = *(f++);
    fmax[0] = *(f++);
    fmin[1] = *(f++);
    fmax[1] = *(f++);
    fmin[2] = *(f++);
    fmax[2] = *(f++);

    dump3f(fmin, "fmin");
    dump3f(fmax, "fmax");
    ms->FDim[0] = *((int *) f++) + 1;
    ms->FDim[1] = *((int *) f++) + 1;
    ms->FDim[2] = *((int *) f++) + 1;

    ms->Div[0] = pymol_roundf((ms->FDim[0] - 1) / (fmax[0] - fmin[0]));
    ms->Div[1] = pymol_roundf((ms->FDim[1] - 1) / (fmax[1] - fmin[1]));
    ms->Div[2] = pymol_roundf((ms->FDim[2] - 1) / (fmax[2] - fmin[2]));

    ms->Min[0] = pymol_roundf(ms->Div[0] * fmin[0]);
    ms->Min[1] = pymol_roundf(ms->Div[1] * fmin[1]);
    ms->Min[2] = pymol_roundf(ms->Div[2] * fmin[2]);

    ms->Max[0] = ms->Min[0] + ms->FDim[0] - 1;
    ms->Max[1] = ms->Min[1] + ms->FDim[1] - 1;
    ms->Max[2] = ms->Min[2] + ms->FDim[2] - 1;
    got_ranges = true;

    {
      int block_len_check;

      block_len_check = *((int *) (f++));
      if(rev_union.block_len != block_len_check) {
        PRINTFB(I->G, FB_ObjectMap, FB_Warnings)
          " ObjectMapGRD-Warning: block length not matched -- not a true GRD binary?\n"
          ENDFB(I->G);
      }
    }
  }

  if(ok) {

    if(ascii) {
      ms->Div[0] = ms->Div[0] - ms->Min[0];
      ms->Div[1] = ms->Div[1] - ms->Min[1];
      ms->Div[2] = ms->Div[2] - ms->Min[2];

      ms->Max[0] = ms->Min[0] + ms->FDim[0] - 1;
      ms->Max[1] = ms->Min[1] + ms->FDim[1] - 1;
      ms->Max[2] = ms->Min[2] + ms->FDim[2] - 1;
    }

    ms->FDim[3] = 3;

    if(Feedback(I->G, FB_ObjectMap, FB_Blather)) {
      dump3i(ms->Div, "ms->Div");
      dump3i(ms->Min, "ms->Min");
      dump3i(ms->Max, "ms->Max");
      dump3i(ms->FDim, "ms->FDim");
    }

    ms->Field.reset(new Isofield(I->G, ms->FDim));
    ms->MapSource = cMapSourceGRD;
    ms->Field->save_points = false;

    switch (fast_axis) {
    case 3:                    /* Fast Y - BROKEN! */
      PRINTFB(I->G, FB_ObjectMap, FB_Warnings)
        " ObjectMapGRD-Warning: fast_axis %d unsupported!\n", fast_axis ENDFB(I->G);
      /* intentional fall though... */
    case 1:                    /* Fast X */
    default:
      for(c = 0; c < ms->FDim[2]; c++) {
        v[2] = (c + ms->Min[2]) / ((float) ms->Div[2]);
        for(b = 0; b < ms->FDim[1]; b++) {
          if(!ascii)
            f++;                /* skip block delimiter */
          v[1] = (b + ms->Min[1]) / ((float) ms->Div[1]);
          for(a = 0; a < ms->FDim[0]; a++) {
            v[0] = (a + ms->Min[0]) / ((float) ms->Div[0]);
            if(ascii) {
              p = ParseNextLine(p);
              p = ParseNCopy(cc, p, 24);
              if(sscanf(cc, "%f", &dens) != 1) {
                ok = false;
              }
            } else {
              dens = *(f++);
            }
            if(ok) {
              sumsq += dens * dens;
              sum += dens;
              n_pts++;
              F3(ms->Field->data, a, b, c) = dens;
              if(maxd < dens)
                maxd = dens;
              if(mind > dens)
                mind = dens;
            }
            transform33f3f(ms->Symmetry->Crystal.fracToReal(), v, vr);
            for(e = 0; e < 3; e++) {
              F4(ms->Field->points, a, b, c, e) = vr[e];
            }
          }
          if(!ascii)
            f++;                /* skip fortran block delimiter */
        }
      }
      break;
    }
  }

  if(n_pts > 1) {
    mean = (float) (sum / n_pts);
    stdev = (float) sqrt1d((sumsq - (sum * sum / n_pts)) / (n_pts - 1));

    if(normalize) {
      PRINTFB(I->G, FB_ObjectMap, FB_Details)
        " ObjectMapGRDStrToMap: Normalizing: mean = %8.6f and stdev = %8.6f.\n",
        mean, stdev ENDFB(I->G);
      if(stdev < R_SMALL8)
        stdev = 1.0F;
      for(c = 0; c < ms->FDim[2]; c++)
        for(b = 0; b < ms->FDim[1]; b++) {
          for(a = 0; a < ms->FDim[0]; a++) {
            dens = F3(ms->Field->data, a, b, c);
            F3(ms->Field->data, a, b, c) = (dens - mean) / stdev;
          }
        }
    } else {
      PRINTFB(I->G, FB_ObjectMap, FB_Details)
        " ObjectMapGRDStrToMap: Mean = %8.6f and stdev = %8.6f.\n",
        mean, stdev ENDFB(I->G);
    }
  }

  if(ok) {
    d = 0;
    for(c = 0; c < ms->FDim[2]; c += (ms->FDim[2] - 1)) {
      v[2] = (c + ms->Min[2]) / ((float) ms->Div[2]);
      for(b = 0; b < ms->FDim[1]; b += (ms->FDim[1] - 1)) {
        v[1] = (b + ms->Min[1]) / ((float) ms->Div[1]);
        for(a = 0; a < ms->FDim[0]; a += (ms->FDim[0] - 1)) {
          v[0] = (a + ms->Min[0]) / ((float) ms->Div[0]);
          transform33f3f(ms->Symmetry->Crystal.fracToReal(), v, vr);
          copy3f(vr, ms->Corner + 3 * d);
          d++;
        }
      }
    }

    v[2] = (ms->Min[2]) / ((float) ms->Div[2]);
    v[1] = (ms->Min[1]) / ((float) ms->Div[1]);
    v[0] = (ms->Min[0]) / ((float) ms->Div[0]);

    transform33f3f(ms->Symmetry->Crystal.fracToReal(), v, ms->ExtentMin);

    v[2] = ((ms->FDim[2] - 1) + ms->Min[2]) / ((float) ms->Div[2]);
    v[1] = ((ms->FDim[1] - 1) + ms->Min[1]) / ((float) ms->Div[1]);
    v[0] = ((ms->FDim[0] - 1) + ms->Min[0]) / ((float) ms->Div[0]);

    transform33f3f(ms->Symmetry->Crystal.fracToReal(), v, ms->ExtentMax);

    PRINTFB(I->G, FB_ObjectMap, FB_Details)
      " GRDXStrToMap: Map Size %d x %d x %d\n", ms->FDim[0], ms->FDim[1], ms->FDim[2]
      ENDFB(I->G);

    ms->Active = true;
    ObjectMapUpdateExtents(I);
    if(!quiet) {
      PRINTFB(I->G, FB_ObjectMap, FB_Results)
        " ObjectMap: Map read.  Range: %5.3f to %5.3f\n", mind, maxd ENDFB(I->G);

      CrystalDump(&ms->Symmetry->Crystal);
    }
  } else {
    PRINTFB(I->G, FB_ObjectMap, FB_Errors)
      " Error: unable to read GRD file.\n" ENDFB(I->G);
  }

  return (ok);

}


/*========================================================================*/
static ObjectMap *ObjectMapReadXPLORStr(PyMOLGlobals * G, ObjectMap * I, char *XPLORStr,
                                        int state, int quiet)
{
  int ok = true;
  int isNew = true;

  if(!I)
    isNew = true;
  else
    isNew = false;
  if(ok) {
    if(isNew) {
      I = (ObjectMap *) new ObjectMap(G);
      isNew = true;
    } else {
      isNew = false;
    }
    ObjectMapXPLORStrToMap(I, XPLORStr, state, quiet);

    SceneChanged(I->G);
    SceneCountFrames(I->G);
  }
  return (I);
}


/*========================================================================*/
static ObjectMap *ObjectMapReadCCP4Str(PyMOLGlobals * G, ObjectMap * I, char *XPLORStr,
                                       int bytes, int state, int quiet,
                                       int format)
{
  int ok = true;
  int isNew = true;

  if(!I)
    isNew = true;
  else
    isNew = false;
  if(ok) {
    if(isNew) {
      I = (ObjectMap *) new ObjectMap(G);
      isNew = true;
    } else {
      isNew = false;
    }
    ObjectMapCCP4StrToMap(I, XPLORStr, bytes, state, quiet, format);
    SceneChanged(G);
    SceneCountFrames(G);
  }
  return (I);
}


/*========================================================================*/
ObjectMap *ObjectMapLoadCCP4(PyMOLGlobals * G, ObjectMap * obj, const char *fname, int state,
                             int is_string, int bytes, int quiet,
                             int format)
{
  ObjectMap *I = NULL;
  char *buffer;
  long size;

  if(!is_string) {
    if (!quiet)
      PRINTFB(G, FB_ObjectMap, FB_Actions)
        " ObjectMapLoadCCP4File: Loading from '%s'.\n", fname ENDFB(G);

    buffer = FileGetContents(fname, &size);

    if(!buffer)
      ErrMessage(G, "ObjectMapLoadCCP4File", "Unable to open file!");
  } else {
    buffer = (char*) fname;
    size = (long) bytes;
  }

  if (buffer) {
    I = ObjectMapReadCCP4Str(G, obj, buffer, size, state, quiet, format);

    if(!is_string)
      mfree(buffer);

    if(!quiet) {
      if(state < 0)
        state = I->State.size() - 1;
      if(state < I->State.size()) {
        ObjectMapState *ms;
        ms = &I->State[state];
        if(ms->Active) {
          CrystalDump(&ms->Symmetry->Crystal);
        }
      }
    }
  }
  return (I);
}


/*========================================================================*/
static ObjectMap *ObjectMapReadFLDStr(PyMOLGlobals * G, ObjectMap * I, char *MapStr,
                                      int bytes, int state, int quiet)
{
  int ok = true;
  int isNew = true;

  if(!I)
    isNew = true;
  else
    isNew = false;
  if(ok) {
    if(isNew) {
      I = (ObjectMap *) new ObjectMap(G);
      isNew = true;
    } else {
      isNew = false;
    }
    ObjectMapFLDStrToMap(I, MapStr, bytes, state, quiet);
    SceneChanged(G);
    SceneCountFrames(G);
  }
  return (I);
}


/*========================================================================*/
static ObjectMap *ObjectMapReadBRIXStr(PyMOLGlobals * G, ObjectMap * I, char *MapStr,
                                       int bytes, int state, int quiet)
{
  int ok = true;
  int isNew = true;

  if(!I)
    isNew = true;
  else
    isNew = false;
  if(ok) {
    if(isNew) {
      I = (ObjectMap *) new ObjectMap(G);
      isNew = true;
    } else {
      isNew = false;
    }
    ObjectMapBRIXStrToMap(I, MapStr, bytes, state, quiet);
    SceneChanged(G);
    SceneCountFrames(G);
  }
  return (I);
}


/*========================================================================*/
static ObjectMap *ObjectMapReadGRDStr(PyMOLGlobals * G, ObjectMap * I, char *MapStr,
                                      int bytes, int state, int quiet)
{
  int ok = true;
  int isNew = true;

  if(!I)
    isNew = true;
  else
    isNew = false;
  if(ok) {
    if(isNew) {
      I = (ObjectMap *) new ObjectMap(G);
      isNew = true;
    } else {
      isNew = false;
    }
    ObjectMapGRDStrToMap(I, MapStr, bytes, state, quiet);
    SceneChanged(G);
    SceneCountFrames(G);
  }
  return (I);
}


/*========================================================================*/
static ObjectMap *ObjectMapReadPHIStr(PyMOLGlobals * G, ObjectMap * I, char *MapStr,
                                      int bytes, int state, int quiet)
{
  int ok = true;
  int isNew = true;

  if(!I)
    isNew = true;
  else
    isNew = false;
  if(ok) {
    if(isNew) {
      I = (ObjectMap *) new ObjectMap(G);
      isNew = true;
    } else {
      isNew = false;
    }
    ObjectMapPHIStrToMap(I, MapStr, bytes, state, quiet);
    SceneChanged(G);
    SceneCountFrames(G);
  }
  return (I);
}


/*========================================================================*/
ObjectMap *ObjectMapLoadPHI(PyMOLGlobals * G, ObjectMap * obj, const char *fname, int state,
                            int is_string, int bytes, int quiet)
{

  ObjectMap *I = NULL;
  long size;
  char *buffer;

  if(!is_string) {
    if (!quiet)
      PRINTFB(G, FB_ObjectMap, FB_Actions)
        " ObjectMapLoadPHIFile: Loading from '%s'.\n", fname ENDFB(G);

    buffer = FileGetContents(fname, &size);

    if(!buffer)
      ErrMessage(G, "ObjectMapLoadPHIFile", "Unable to open file!");
  } else {
    buffer = (char*) fname;
    size = (long) bytes;
  }

  if(buffer) {
    I = ObjectMapReadPHIStr(G, obj, buffer, size, state, quiet);

    if(!is_string)
      mfree(buffer);

  }
  return (I);

}


/*========================================================================*/

static int is_number(char *p)
{
  int result = (*p != 0);
  char c;
  if(result)
    while((c = *(p++))) {
      if(!((c == '.') ||
           (c == '-') ||
           (c == '+') || (c == 'e') || (c == 'E') || ((c >= '0') && (c <= '9'))))
        result = false;
      break;
    }
  return result;
}

static pymol::Result<std::unique_ptr<ObjectMapState>> ObjectMapDXStrToMap(
    PyMOLGlobals* G, const char* DXStr, int bytes, bool quiet)
{
  if (DXStr[0] == '\x1f' && DXStr[1] == '\x8b') {
    return pymol::Error("gzipped data not supported");
  }

  enum {
    ASCII_NUMERIC,
    BINARY_DOUBLE,
    BINARY_FLOAT,
  } data_type = ASCII_NUMERIC;

  int n_items = 0;

  const char *p, *pp;
  float dens;
  int a, b, c, d, e;
  float v[3], maxd, mind;
  /* DX named from their docs */

  int stage = 0;

  char cc[MAXLINELEN];

  auto ms = pymol::make_unique<ObjectMapState>(G);

  ms->Origin = std::vector<float>(3);
  ms->Grid = std::vector<float>(3);

  maxd = -FLT_MAX;
  mind = FLT_MAX;
  p = DXStr;

  /* get the dimensions */

  ms->FDim[3] = 3;

  while ((*p) && stage == 0) {
    pp = p;
    p = ParseNCopy(cc, p, 35);
    if((strcmp(cc, "object 1 class gridpositions counts") == 0) || is_number(cc)) {
      if(is_number(cc))
        p = pp;
      p = ParseWordCopy(cc, p, 10);
      if(sscanf(cc, "%d", &ms->FDim[0]) == 1) {
        p = ParseWordCopy(cc, p, 10);
        if(sscanf(cc, "%d", &ms->FDim[1]) == 1) {
          p = ParseWordCopy(cc, p, 10);
          if(sscanf(cc, "%d", &ms->FDim[2]) == 1) {
            stage = 1;
          }
        }
      }
    }
    p = ParseNextLine(p);
  }

  if (stage != 1) {
    return pymol::Error("missing 'object 1 class gridpositions count' line");
  }

  {
    PRINTFB(G, FB_ObjectMap, FB_Details)
      " DXStrToMap: Dimensions: %d %d %d\n", ms->FDim[0], ms->FDim[1], ms->FDim[2]
      ENDFB(G);
  }

  /* get the origin */

  while ((*p) && stage == 1) {
    pp = p;
    p = ParseNCopy(cc, p, 6);
    if((strcmp(cc, "origin") == 0) || is_number(cc)) {
      if(is_number(cc))
        p = pp;
      p = ParseWordCopy(cc, p, 20);
      if(sscanf(cc, "%f", &ms->Origin[0]) == 1) {
        p = ParseWordCopy(cc, p, 20);
        if(sscanf(cc, "%f", &ms->Origin[1]) == 1) {
          p = ParseWordCopy(cc, p, 20);
          if(sscanf(cc, "%f", &ms->Origin[2]) == 1) {
            stage = 2;
          }
        }
      }
    }
    p = ParseNextLine(p);
  }

  if (stage != 2) {
    return pymol::Error("missing 'origin'");
  }

  {
    PRINTFB(G, FB_ObjectMap, FB_Details)
      " DXStrToMap: Origin %8.3f %8.3f %8.3f\n", ms->Origin[0], ms->Origin[1],
      ms->Origin[2]
      ENDFB(G);
  }

  float delta[9];
  int delta_i = 0;

  while ((*p) && stage == 2) {
    pp = p;
    p = ParseNCopy(cc, p, 5);

    if(strcmp(cc, "delta") != 0) {
      if(is_number(cc)) {
        p = pp;
      } else {
        p = ParseNextLine(p);
        continue;
      }
    }

    if(3 != sscanf(p, " %f %f %f",
          delta + delta_i,
          delta + delta_i + 1,
          delta + delta_i + 2)) {
      return pymol::Error("expected 3 floats");
    }

    p = ParseNextLine(p);
    delta_i += 3;

    if (delta_i == 9) {
      stage = 3;

      if (is_diagonalf(3, delta)) {
        ms->Grid[0] = delta[0];
        ms->Grid[1] = delta[4];
        ms->Grid[2] = delta[8];
      } else {
        if(ms->Matrix.empty())
          ms->Matrix = std::vector<double>(16);

        copy33f44d(delta, ms->Matrix.data());
        ms->Matrix[3] = ms->Origin[0];
        ms->Matrix[7] = ms->Origin[1];
        ms->Matrix[11] = ms->Origin[2];

        ones3f(ms->Grid.data());
        zero3f(ms->Origin.data());
      }
    }
  }

  if (stage != 3) {
    return pymol::Error("missing 'delta'");
  }

  {
    PRINTFB(G, FB_ObjectMap, FB_Details)
      " DXStrToMap: Grid %8.3f %8.3f %8.3f\n", ms->Grid[0], ms->Grid[1], ms->Grid[2]
      ENDFB(G);
  }

  while((*p) && stage == 3) {
    p = ParseNCopy(cc, p, 6);
    if(strcmp(cc, "object") == 0) {
      p = ParseWordCopy(cc, p, 20);
      char cc_type[8] = "";
      char cc_binary[8] = "";
      if (2 <= sscanf(p, " class array type %7s rank %*s items %s %7s", cc_type,
                   cc, cc_binary)) {
        if(sscanf(cc, "%d", &n_items) == 1) {
          if(n_items == ms->FDim[0] * ms->FDim[1] * ms->FDim[2])
            stage = 4;
        }

        if (strcmp(cc_binary, "binary") == 0) {
          if (strcmp(cc_type, "double") == 0) {
            data_type = BINARY_DOUBLE;
          } else if (strcmp(cc_type, "float") == 0) {
            data_type = BINARY_FLOAT;
          } else {
            return pymol::Error(
                pymol::string_format("type '%s' not supported", cc_type));
          }
        }
      }
    } else if(is_number(cc)) {
      n_items = ms->FDim[0] * ms->FDim[1] * ms->FDim[2];
      stage = 4;
      break;
    }
    p = ParseNextLine(p);
  }

  if (stage != 4) {
    return pymol::Error("missing 'object . class array'");
  }

  PRINTFB(G, FB_ObjectMap, FB_Details)
    " DXStrToMap: %d data points.\n", n_items ENDFB(G);

  unsigned malformed_ascii_floats_count = 0;

  {

    ms->Field.reset(new Isofield(G, ms->FDim));
    ms->MapSource = cMapSourceGeneralPurpose;
    ms->Field->save_points = false;

    for(a = 0; a < 3; a++) {
      ms->Div[a] = ms->FDim[a] - 1;
      ms->Min[a] = 0;
      ms->Max[a] = ms->FDim[a] - 1;
    }

    for(a = 0; a < ms->FDim[0]; a++) {
      for(b = 0; b < ms->FDim[1]; b++) {
        for(c = 0; c < ms->FDim[2]; c++) {

          switch (data_type) {
          case BINARY_DOUBLE:
            dens = *reinterpret_cast<double const*>(p);
            p += sizeof(double);
            break;
          case BINARY_FLOAT:
            dens = *reinterpret_cast<float const*>(p);
            p += sizeof(float);
            break;
          case ASCII_NUMERIC:
            p = ParseWordCopy(cc, p, 20);
            if(!cc[0]) {
              p = ParseNextLine(p);
              p = ParseWordCopy(cc, p, 20);
            }
            if (sscanf(cc, "%f", &dens) != 1) {
              ++malformed_ascii_floats_count;
              continue;
            }
            break;
          default:
            assert(false);
          }

          if(maxd < dens)
            maxd = dens;
          if(mind > dens)
            mind = dens;
          F3(ms->Field->data, a, b, c) = dens;
        }
      }

    }

    for(e = 0; e < 3; e++) {
      ms->ExtentMin[e] = ms->Origin[e] + ms->Grid[e] * ms->Min[e];
      ms->ExtentMax[e] = ms->Origin[e] + ms->Grid[e] * ms->Max[e];
    }

    for(c = 0; c < ms->FDim[2]; c++) {
      v[2] = ms->Origin[2] + ms->Grid[2] * (c + ms->Min[2]);
      for(b = 0; b < ms->FDim[1]; b++) {
        v[1] = ms->Origin[1] + ms->Grid[1] * (b + ms->Min[1]);
        for(a = 0; a < ms->FDim[0]; a++) {
          v[0] = ms->Origin[0] + ms->Grid[0] * (a + ms->Min[0]);
          for(e = 0; e < 3; e++) {
            F4(ms->Field->points, a, b, c, e) = v[e];
          }
        }
      }
    }

    d = 0;
    for(c = 0; c < ms->FDim[2]; c += ms->FDim[2] - 1) {
      v[2] = ms->Origin[2] + ms->Grid[2] * (c + ms->Min[2]);

      for(b = 0; b < ms->FDim[1]; b += ms->FDim[1] - 1) {
        v[1] = ms->Origin[1] + ms->Grid[1] * (b + ms->Min[1]);

        for(a = 0; a < ms->FDim[0]; a += ms->FDim[0] - 1) {
          v[0] = ms->Origin[0] + ms->Grid[0] * (a + ms->Min[0]);
          copy3f(v, ms->Corner + 3 * d);
          d++;
        }
      }
    }

    if(!quiet) {
      PRINTFB(G, FB_ObjectMap, FB_Results)
        " ObjectMap: Map read.  Range: %5.3f to %5.3f\n", mind, maxd ENDFB(G);
    }
  }

  if (malformed_ascii_floats_count > 0) {
    PRINTFB(G, FB_ObjectMap, FB_Warnings)
    " %s-Warning: Failed to parse %u values\n", __func__,
        malformed_ascii_floats_count ENDFB(G);
  }

  ms->Active = true;
  return std::move(ms);
}


/*========================================================================*/
ObjectMap* ObjectMapReadDXStr(PyMOLGlobals* G, ObjectMap* I, const char* MapStr,
    int bytes, int state, bool quiet)
{
  auto mapstateresult = ObjectMapDXStrToMap(G, MapStr, bytes, quiet);
  if (!mapstateresult) {
    ErrMessage(G, __func__, mapstateresult.error().what().c_str());
    return nullptr;
  }

  assert(mapstateresult.result());
  assert(mapstateresult.result()->Active);

  if (!I) {
    I = new ObjectMap(G);
  }

  if (state < 0) {
    state = I->State.size();
  }

  VecCheckEmplace(I->State, state, G);

  I->State[state] = std::move(*mapstateresult.result());

  ObjectMapUpdateExtents(I);
  SceneChanged(G);
  SceneCountFrames(G);

  return (I);
}


/*========================================================================*/
ObjectMap *ObjectMapLoadDXFile(PyMOLGlobals * G, ObjectMap * obj, const char *fname, int state,
                               int quiet)
{
  ObjectMap *I = NULL;
  long size;
  char *buffer;

  buffer = FileGetContents(fname, &size);

  if(!buffer) {
    ErrMessage(G, "ObjectMapLoadDXFile", "Unable to open file!");
    PRINTFB(G, FB_ObjectMap, FB_Errors)
      "ObjectMapLoadDXFile: Does '%s' exist?\n", fname ENDFB(G);
  } else {
    if(Feedback(G, FB_ObjectMap, FB_Actions)) {
      printf(" ObjectMapLoadDXFile: Loading from '%s'.\n", fname);
    }

    I = ObjectMapReadDXStr(G, obj, buffer, size, state, quiet);

    mfree(buffer);
  }
  return (I);

}


/*========================================================================*/
static int ObjectMapACNTStrToMap(ObjectMap * I, char *ACNTStr, int bytes, int state,
                                 int quiet)
{

  int n_items = 0;

  char *p;
  float dens;
  int a, b, c, d, e;
  float v[3], maxd, mind;
  int ok = true;

  int stage = 0;

  ObjectMapState *ms;

  char cc[MAXLINELEN];

  if(state < 0)
    state = I->State.size();
  if(I->State.size() <= state) {
    VecCheckEmplace(I->State, state, I->G);
  }

  ms = &I->State[state];

  ms->Origin = std::vector<float>(3);
  ms->Grid = std::vector<float>(3);

  maxd = -FLT_MAX;
  mind = FLT_MAX;
  p = ACNTStr;

  /* skip header */

  p = ParseNextLine(p);

  /* get the origin, spacing, and dimensions */

  ms->FDim[3] = 3;

  /* read the map info */

  {
    int i;
    for(i = 0; i < 3; i++) {
      int ii = (4 - i) % 3;     /* y, x , z */
      p = ParseWordCopy(cc, p, MAXLINELEN);
      if(sscanf(cc, "%f", &ms->Origin[ii]) == 1) {
        p = ParseWordCopy(cc, p, MAXLINELEN);
        if(sscanf(cc, "%f", &ms->Grid[ii]) == 1) {
          p = ParseWordCopy(cc, p, MAXLINELEN);
          if(sscanf(cc, "%d", &ms->FDim[ii]) == 1) {
            p = ParseNextLine(p);
            stage++;
          }
        }
      }
    }
  }

  /* skip the angles (ignored -- should probably check instead */

  p = ParseNextLine(p);

  /* read the data block */

  if(ok && (stage == 3)) {

    PRINTFB(I->G, FB_ObjectMap, FB_Details)
      " ACNTStrToMap: Dimensions: %d %d %d\n", ms->FDim[0], ms->FDim[1], ms->FDim[2]
      ENDFB(I->G);
    PRINTFB(I->G, FB_ObjectMap, FB_Details)
      " ACNTStrToMap: Origin %8.3f %8.3f %8.3f\n", ms->Origin[0], ms->Origin[1],
      ms->Origin[2]
      ENDFB(I->G);
    PRINTFB(I->G, FB_ObjectMap, FB_Details)
      " ACNTStrToMap: Grid %8.3f %8.3f %8.3f\n", ms->Grid[0], ms->Grid[1], ms->Grid[2]
      ENDFB(I->G);

    n_items = ms->FDim[0] * ms->FDim[1] * ms->FDim[2];

    if(ok && (stage == 1)) {
      PRINTFB(I->G, FB_ObjectMap, FB_Details)
        " ACNTStrToMap: %d data points.\n", n_items ENDFB(I->G);
    }

    ms->Field.reset(new Isofield(I->G, ms->FDim));
    ms->MapSource = cMapSourceGeneralPurpose;
    ms->Field->save_points = false;

    for(a = 0; a < 3; a++) {
      ms->Div[a] = ms->FDim[a] - 1;
      ms->Min[a] = 0;
      ms->Max[a] = ms->FDim[a] - 1;
    }

    for(c = 0; c < ms->FDim[2]; c++) {
      for(a = 0; a < ms->FDim[0]; a++) {
        for(b = 0; b < ms->FDim[1]; b++) {
          p = ParseWordCopy(cc, p, MAXLINELEN);
          p = ParseNextLine(p);
          if(sscanf(cc, "%f", &dens) == 1) {
            if(maxd < dens)
              maxd = dens;
            if(mind > dens)
              mind = dens;
            F3(ms->Field->data, a, b, c) = dens;
          } else {
            ok = false;
          }
        }
      }
    }

    for(e = 0; e < 3; e++) {
      ms->ExtentMin[e] = ms->Origin[e] + ms->Grid[e] * ms->Min[e];
      ms->ExtentMax[e] = ms->Origin[e] + ms->Grid[e] * ms->Max[e];
    }

    for(c = 0; c < ms->FDim[2]; c++) {
      v[2] = ms->Origin[2] + ms->Grid[2] * (c + ms->Min[2]);
      for(b = 0; b < ms->FDim[1]; b++) {
        v[1] = ms->Origin[1] + ms->Grid[1] * (b + ms->Min[1]);
        for(a = 0; a < ms->FDim[0]; a++) {
          v[0] = ms->Origin[0] + ms->Grid[0] * (a + ms->Min[0]);
          for(e = 0; e < 3; e++) {
            F4(ms->Field->points, a, b, c, e) = v[e];
          }
        }
      }
    }

    d = 0;
    for(c = 0; c < ms->FDim[2]; c += ms->FDim[2] - 1) {
      v[2] = ms->Origin[2] + ms->Grid[2] * (c + ms->Min[2]);

      for(b = 0; b < ms->FDim[1]; b += ms->FDim[1] - 1) {
        v[1] = ms->Origin[1] + ms->Grid[1] * (b + ms->Min[1]);

        for(a = 0; a < ms->FDim[0]; a += ms->FDim[0] - 1) {
          v[0] = ms->Origin[0] + ms->Grid[0] * (a + ms->Min[0]);
          copy3f(v, ms->Corner + 3 * d);
          d++;
        }
      }
    }

    if(ok)
      stage = 4;
  }

  if(stage != 4)
    ok = false;

  if(!ok) {
    ErrMessage(I->G, "ObjectMap", "Error reading map");
  } else {
    ms->Active = true;
    ObjectMapUpdateExtents(I);
    if(!quiet) {
      PRINTFB(I->G, FB_ObjectMap, FB_Results)
        " ObjectMap: Map read.  Range: %5.3f to %5.3f\n", mind, maxd ENDFB(I->G);
    }
  }
  return (ok);
}

static ObjectMap *ObjectMapReadACNTStr(PyMOLGlobals * G, ObjectMap * I,
                                       char *MapStr, int bytes, int state, int quiet)
{
  int ok = true;
  int isNew = true;

  if(!I)
    isNew = true;
  else
    isNew = false;
  if(ok) {
    if(isNew) {
      I = (ObjectMap *) new ObjectMap(G);
      isNew = true;
    } else {
      isNew = false;
    }
    ObjectMapACNTStrToMap(I, MapStr, bytes, state, quiet);
    SceneChanged(G);
    SceneCountFrames(G);
  }
  return (I);
}

ObjectMap *ObjectMapLoadACNTFile(PyMOLGlobals * G, ObjectMap * obj, const char *fname,
                                 int state, int quiet)
{
  ObjectMap *I = NULL;
  long size;
  char *buffer;
  
  buffer = FileGetContents(fname, &size);

  if(!buffer) {
    ErrMessage(G, "ObjectMapLoadACNTFile", "Unable to open file!");
    PRINTFB(G, FB_ObjectMap, FB_Errors)
      "ObjectMapLoadACNTFile: Does '%s' exist?\n", fname ENDFB(G);
  } else {
    if(Feedback(G, FB_ObjectMap, FB_Actions)) {
      printf(" ObjectMapLoadACNTFile: Loading from '%s'.\n", fname);
    }

    I = ObjectMapReadACNTStr(G, obj, buffer, size, state, quiet);

    mfree(buffer);
  }
  return (I);

}


/*========================================================================*/
ObjectMap *ObjectMapLoadFLDFile(PyMOLGlobals * G, ObjectMap * obj, const char *fname, int state,
                                int quiet)
{
  ObjectMap *I = NULL;
  long size;
  char *buffer;

  buffer = FileGetContents(fname, &size);

  if(!buffer) {
    ErrMessage(G, "ObjectMapLoadFLDFile", "Unable to open file!");
  } else {
    if(Feedback(G, FB_ObjectMap, FB_Actions)) {
      printf(" ObjectMapLoadFLDFile: Loading from '%s'.\n", fname);
    }

    I = ObjectMapReadFLDStr(G, obj, buffer, size, state, quiet);

    mfree(buffer);
  }
  return (I);

}


/*========================================================================*/
ObjectMap *ObjectMapLoadBRIXFile(PyMOLGlobals * G, ObjectMap * obj, const char *fname,
                                 int state, int quiet)
{
  ObjectMap *I = NULL;
  long size;
  char *buffer;

  buffer = FileGetContents(fname, &size);

  if(!buffer) {
    ErrMessage(G, "ObjectMapLoadBRIXFile", "Unable to open file!");
  } else {
    if(Feedback(G, FB_ObjectMap, FB_Actions)) {
      printf(" ObjectMapLoadBRIXFile: Loading from '%s'.\n", fname);
    }

    I = ObjectMapReadBRIXStr(G, obj, buffer, size, state, quiet);

    mfree(buffer);
  }
  return (I);

}


/*========================================================================*/
ObjectMap *ObjectMapLoadGRDFile(PyMOLGlobals * G, ObjectMap * obj, const char *fname, int state,
                                int quiet)
{
  ObjectMap *I = NULL;
  long size;
  char *buffer;

  buffer = FileGetContents(fname, &size);

  if(!buffer) {
    ErrMessage(G, "ObjectMapLoadGRDFile", "Unable to open file!");
  } else {
    if(Feedback(G, FB_ObjectMap, FB_Actions)) {
      printf(" ObjectMapLoadGRDFile: Loading from '%s'.\n", fname);
    }

    I = ObjectMapReadGRDStr(G, obj, buffer, size, state, quiet);

    mfree(buffer);
  }
  return (I);

}


/*========================================================================*/
ObjectMap *ObjectMapLoadXPLOR(PyMOLGlobals * G, ObjectMap * obj, const char *fname,
                              int state, int is_file, int quiet)
{
  ObjectMap *I = NULL;
  long size;
  char *buffer;

  if(is_file) {
    buffer = FileGetContents(fname, &size);

    if(!buffer)
      ErrMessage(G, "ObjectMapLoadXPLOR", "Unable to open file!");
  } else {
    buffer = (char*) fname;
  }

  if (buffer) {
    if((!quiet) && (Feedback(G, FB_ObjectMap, FB_Actions))) {
      if(is_file) {
        printf(" ObjectMapLoadXPLOR: Loading from '%s'.\n", fname);
      } else {
        printf(" ObjectMapLoadXPLOR: Loading...\n");
      }
    }

    I = ObjectMapReadXPLORStr(G, obj, buffer, state, quiet);

    if(is_file)
      mfree(buffer);
    if(!quiet) {
      if(Feedback(G, FB_ObjectMap, FB_Details)) {
        if(state < 0)
          state = I->State.size() - 1;

        if(state < I->State.size()) {
          ObjectMapState *ms;
          ms = &I->State[state];
          if(ms->Active) {
            CrystalDump(&ms->Symmetry->Crystal);
          }
        }
      }
    }
  }
  return (I);

}


/*========================================================================*/
int ObjectMapSetBorder(ObjectMap * I, float level, int state)
{
  for (StateIterator iter(I, state); iter.next();) {
    auto* oms = &I->State[iter.state];
    if (oms->Active) {
      if (!ObjectMapStateSetBorder(oms, level))
        return false;
    }
  }
  return true;
}


/*========================================================================*/
#ifndef _PYMOL_NOPY
static int ObjectMapNumPyArrayToMapState(PyMOLGlobals * G, ObjectMapState * ms,
                                         PyObject * ary, int quiet)
{

  int a, b, c, d, e;
  float v[3], dens, maxd, mind;
  int ok = true;
  void * ptr;

#ifdef _PYMOL_NUMPY
  PyArrayObject * pao = (PyArrayObject *) ary;
  const int itemsize = PyArray_ITEMSIZE(pao);
#endif
  maxd = -FLT_MAX;
  mind = FLT_MAX;
  if(ok) {
    ms->FDim[0] = ms->Dim[0];
    ms->FDim[1] = ms->Dim[1];
    ms->FDim[2] = ms->Dim[2];
    ms->FDim[3] = 3;

    if(!(ms->FDim[0] && ms->FDim[1] && ms->FDim[2]))
      ok = false;
    else {
      ms->Field.reset(new Isofield(G, ms->FDim));
      for(c = 0; c < ms->FDim[2]; c++) {
        v[2] = ms->Origin[2] + ms->Grid[2] * c;
        for(b = 0; b < ms->FDim[1]; b++) {
          v[1] = ms->Origin[1] + ms->Grid[1] * b;
          for(a = 0; a < ms->FDim[0]; a++) {
            v[0] = ms->Origin[0] + ms->Grid[0] * a;
#ifdef _PYMOL_NUMPY
            ptr = PyArray_GETPTR3(pao, a, b, c);
            switch(itemsize) {
              case sizeof(double):
                dens = (float) *((double*)ptr);
                break;
              case sizeof(float):
                dens = *((float*)ptr);
                break;
              default:
                dens = 0.0;
                printf("no itemsize match\n");
            }
#else
            dens = 0.0;
#endif
            F3(ms->Field->data, a, b, c) = dens;
            if(maxd < dens)
              maxd = dens;
            if(mind > dens)
              mind = dens;
            for(e = 0; e < 3; e++)
              F4(ms->Field->points, a, b, c, e) = v[e];
          }
        }
      }
      d = 0;
      for(c = 0; c < ms->FDim[2]; c += (ms->FDim[2] - 1)) {
        v[2] = ms->Origin[2] + ms->Grid[2] * c;
        for(b = 0; b < ms->FDim[1]; b += (ms->FDim[1] - 1)) {
          v[1] = ms->Origin[1] + ms->Grid[1] * b;
          for(a = 0; a < ms->FDim[0]; a += (ms->FDim[0] - 1)) {
            v[0] = ms->Origin[0] + ms->Grid[0] * a;
            copy3f(v, ms->Corner + 3 * d);
            d++;
          }
        }
      }
    }
  }
  if(ok) {
    copy3f(ms->Origin.data(), ms->ExtentMin);
    copy3f(ms->Origin.data(), ms->ExtentMax);
    add3f(ms->Range.data(), ms->ExtentMax, ms->ExtentMax);
  }
  if(!ok) {
    ErrMessage(G, "ObjectMap", "Error reading map");
  } else {
    ms->Active = true;
    if(!quiet) {
      PRINTFB(G, FB_ObjectMap, FB_Results)
        " ObjectMap: Map read.  Range: %5.3f to %5.3f\n", mind, maxd ENDFB(G);
    }
  }
  return (ok);
}
#endif

/*========================================================================*/
ObjectMap *ObjectMapLoadChemPyBrick(PyMOLGlobals * G, ObjectMap * I, PyObject * Map,
                                    int state, int discrete, int quiet)
{
#ifndef _PYMOL_NUMPY
  ErrMessage(G, "ObjectMap", "missing numpy support, required for chempy.brick");
  return NULL;
#endif

#ifdef _PYMOL_NOPY
  return NULL;
#else

  int ok = true;
  int isNew = true;
  PyObject *tmp;
  ObjectMapState *ms;

  if(!I)
    isNew = true;
  else
    isNew = false;

  if(ok) {

    if(isNew) {
      I = (ObjectMap *) new ObjectMap(G);
      isNew = true;
    } else {
      isNew = false;
    }

    if(state < 0)
      state = I->State.size();
    if(I->State.size() <= state) {
      VecCheckEmplace(I->State, state, I->G);
    }
    ms = &I->State[state];

    if(PyObject_HasAttrString(Map, "origin") &&
       PyObject_HasAttrString(Map, "dim") &&
       PyObject_HasAttrString(Map, "range") &&
       PyObject_HasAttrString(Map, "grid") && PyObject_HasAttrString(Map, "lvl")) {
      tmp = PyObject_GetAttrString(Map, "origin");
      if(tmp) {
        PConvFromPyObject(G, tmp, ms->Origin);
        Py_DECREF(tmp);
      } else
        ok = ErrMessage(G, "ObjectMap", "missing brick origin.");
      tmp = PyObject_GetAttrString(Map, "dim");
      if(tmp) {
        PConvFromPyObject(G, tmp, ms->Dim);
        Py_DECREF(tmp);
      } else
        ok = ErrMessage(G, "ObjectMap", "missing brick dimension.");
      tmp = PyObject_GetAttrString(Map, "range");
      if(tmp) {
        PConvFromPyObject(G, tmp, ms->Range);
        Py_DECREF(tmp);
      } else
        ok = ErrMessage(G, "ObjectMap", "missing brick range.");
      tmp = PyObject_GetAttrString(Map, "grid");
      if(tmp) {
        PConvFromPyObject(G, tmp, ms->Grid);
        Py_DECREF(tmp);
      } else
        ok = ErrMessage(G, "ObjectMap", "missing brick grid.");
      tmp = PyObject_GetAttrString(Map, "lvl");
      if(tmp) {

        ObjectMapNumPyArrayToMapState(G, ms, tmp, quiet);

        Py_DECREF(tmp);
      } else
        ok = ErrMessage(G, "ObjectMap", "missing brick density.");

    } else
      ok = ErrMessage(G, "ObjectMap", "missing any brick attribute.");

    SceneChanged(G);
    SceneCountFrames(G);
    if(ok) {
      int a;
      for(a = 0; a < 3; a++) {
        ms->Min[a] = 0;
        ms->Max[a] = ms->Dim[a] - 1;
      }
      ms->Active = true;
      ms->MapSource = cMapSourceChempyBrick;
      ObjectMapUpdateExtents(I);

    }
  }
  return (I);
#endif
}


/*========================================================================*/
ObjectMap *ObjectMapLoadChemPyMap(PyMOLGlobals * G, ObjectMap * I, PyObject * Map,
                                  int state, int discrete, int quiet)
{
#ifdef _PYMOL_NOPY
  return NULL;
#else
  int ok = true;
  int isNew = true;
  float *cobj;
  WordType format;
  float v[3], vr[3], dens, maxd, mind;
  int a, b, c, d, e;
  ObjectMapState *ms;

  maxd = -FLT_MAX;
  mind = FLT_MAX;

  if(!I)
    isNew = true;
  else
    isNew = false;

  if(ok) {

    if(isNew) {
      I = (ObjectMap *) new ObjectMap(G);
      isNew = true;
    } else {
      isNew = false;
    }

    if(state < 0)
      state = I->State.size();
    if(I->State.size() <= state) {
      VecCheckEmplace(I->State, state, I->G);
    }
    ms = &I->State[state];

    float dim[3], angle[3];

    if(!PConvAttrToStrMaxLen(Map, "format", format, sizeof(WordType) - 1))
      ok = ErrMessage(G, "LoadChemPyMap", "bad 'format' parameter.");
    else if(!PConvAttrToFloatArrayInPlace(Map, "cell_dim", dim, 3))
      ok = ErrMessage(G, "LoadChemPyMap", "bad 'cell_dim' parameter.");
    else if(!PConvAttrToFloatArrayInPlace(Map, "cell_ang", angle, 3))
      ok = ErrMessage(G, "LoadChemPyMap", "bad 'cell_ang' parameter.");
    else if(!PConvAttrToIntArrayInPlace(Map, "cell_div", ms->Div, 3))
      ok = ErrMessage(G, "LoadChemPyMap", "bad 'cell_div' parameter.");
    else if(!PConvAttrToIntArrayInPlace(Map, "first", ms->Min, 3))
      ok = ErrMessage(G, "LoadChemPyMap", "bad 'first' parameter.");
    else if(!PConvAttrToIntArrayInPlace(Map, "last", ms->Max, 3))
      ok = ErrMessage(G, "LoadChemPyMap", "bad 'last' parameter.");

    ms->Symmetry->Crystal.setDims(dim);
    ms->Symmetry->Crystal.setAngles(angle);

    if(ok) {
      if(strcmp(format, "CObjectZYXfloat") == 0) {
        ok = PConvAttrToPtr(Map, "c_object", (void **) (void *) &cobj);
        if(!ok)
          ErrMessage(G, "LoadChemPyMap", "pymol::CObject unreadable.");
      } else {
        ok = ErrMessage(G, "LoadChemPyMap", "unsupported format.");
      }
    }
    /* good to go */

    if(ok) {
      if(strcmp(format, "CObjectZYXfloat") == 0) {

        ms->FDim[0] = ms->Max[0] - ms->Min[0] + 1;
        ms->FDim[1] = ms->Max[1] - ms->Min[1] + 1;
        ms->FDim[2] = ms->Max[2] - ms->Min[2] + 1;
        if(Feedback(G, FB_ObjectMap, FB_Actions)) {
          printf(" LoadChemPyMap: CObjectZYXdouble %dx%dx%d\n", ms->FDim[0], ms->FDim[1],
                 ms->FDim[2]);
        }
        ms->FDim[3] = 3;
        if(!(ms->FDim[0] && ms->FDim[1] && ms->FDim[2]))
          ok = false;
        else {
          ms->Field.reset(new Isofield(G, ms->FDim));
          for(c = 0; c < ms->FDim[2]; c++) {
            v[2] = (c + ms->Min[2]) / ((float) ms->Div[2]);
            for(b = 0; b < ms->FDim[1]; b++) {
              v[1] = (b + ms->Min[1]) / ((float) ms->Div[1]);
              for(a = 0; a < ms->FDim[0]; a++) {
                v[0] = (a + ms->Min[0]) / ((float) ms->Div[0]);

                dens = *(cobj++);

                F3(ms->Field->data, a, b, c) = dens;
                if(maxd < dens)
                  maxd = dens;
                if(mind > dens)
                  mind = dens;
                transform33f3f(ms->Symmetry->Crystal.fracToReal(), v, vr);
                for(e = 0; e < 3; e++)
                  F4(ms->Field->points, a, b, c, e) = vr[e];
              }
            }
          }

          if(ok) {
            d = 0;
            for(c = 0; c < ms->FDim[2]; c += (ms->FDim[2] - 1)) {
              v[2] = (c + ms->Min[2]) / ((float) ms->Div[2]);
              for(b = 0; b < ms->FDim[1]; b += (ms->FDim[1] - 1)) {
                v[1] = (b + ms->Min[1]) / ((float) ms->Div[1]);
                for(a = 0; a < ms->FDim[0]; a += (ms->FDim[0] - 1)) {
                  v[0] = (a + ms->Min[0]) / ((float) ms->Div[0]);
                  transform33f3f(ms->Symmetry->Crystal.fracToReal(), v, vr);
                  copy3f(vr, ms->Corner + 3 * d);
                  d++;
                }
              }
            }
          }
        }
      }
    }

    if(ok) {
      CrystalDump(&ms->Symmetry->Crystal);

      v[2] = (ms->Min[2]) / ((float) ms->Div[2]);
      v[1] = (ms->Min[1]) / ((float) ms->Div[1]);
      v[0] = (ms->Min[0]) / ((float) ms->Div[0]);

      transform33f3f(ms->Symmetry->Crystal.fracToReal(), v, ms->ExtentMin);

      v[2] = ((ms->FDim[2] - 1) + ms->Min[2]) / ((float) ms->Div[2]);
      v[1] = ((ms->FDim[1] - 1) + ms->Min[1]) / ((float) ms->Div[1]);
      v[0] = ((ms->FDim[0] - 1) + ms->Min[0]) / ((float) ms->Div[0]);

      transform33f3f(ms->Symmetry->Crystal.fracToReal(), v, ms->ExtentMax);

    }

    if(!ok) {
      ErrMessage(G, "ObjectMap", "Error reading map");
    } else {
      ms->Active = true;
      ObjectMapUpdateExtents(I);
      if(!quiet) {
        PRINTFB(I->G, FB_ObjectMap, FB_Results)
          " ObjectMap: Map read.  Range: %5.3f to %5.3f\n", mind, maxd ENDFB(I->G);
      }
    }

    if(ok) {
      SceneChanged(G);
      SceneCountFrames(G);
    }
  }
  return (I);
#endif
}

void ObjectMapDump(const ObjectMap* om, const char* fname, int state, int quiet)
{
  auto* oms = om->getObjectMapState(state);

  if (!oms) {
    ErrMessage(om->G, __func__, "state out of range");
    return;
  }

  auto file = fopen(fname, "wb");
  if (!file) {
    ErrMessage(om->G, __func__, "can't open file for writing");
    return;
  }

  auto* field = oms->Field.get();

  for (int xi = 0; xi < field->dimensions[0]; xi++) {
    for (int yi = 0; yi < field->dimensions[1]; yi++) {
      for (int zi = 0; zi < field->dimensions[2]; zi++) {

        float x = field->points->get<float>(xi, yi, zi, 0);
        float y = field->points->get<float>(xi, yi, zi, 1);
        float z = field->points->get<float>(xi, yi, zi, 2);

        switch (field->data->type) {
          case cFieldFloat: {
            float value = field->data->get<float>(xi, yi, zi);
            fprintf(file, "%10.4f%10.4f%10.4f%10.4f\n", x, y, z, value);
            break;
          }
          case cFieldInt: {
            int value = field->data->get<int>(xi, yi, zi);
            fprintf(file, "%10.4f%10.4f%10.4f%10d\n", x, y, z, value);
            break;
          }
          default:
            ErrMessage(om->G, __func__, "unknown field type");
            fclose(file);
            return;
        }
      }
    }
  }
  fclose(file);
  if (!quiet) {
    PRINTFB(om->G, FB_ObjectMap, FB_Actions)
      " ObjectMapDump: %s written to %s\n", om->Name, fname ENDFB(om->G);
  }
}

pymol::CObject* ObjectMap::clone() const
{
  return new ObjectMap(*this);
}

