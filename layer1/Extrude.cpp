
/* 
A* -------------------------------------------------------------------
B* This file contains source code for the PyMOL computer program
C* copyright 1998-2000 by Warren Lyford Delano of DeLano Scientific. 
D* -------------------------------------------------------------------
E* It is unlawful to modify or remove this copyright notice.
F* -------------------------------------------------------------------
G* Please see the accompanying LICENSE file for further information. 
H* --------------------------------------------------\-----------------
I* Additional authors of this source file include:
-* Cameron Mura
-* 
-*
Z* -------------------------------------------------------------------
*/

#include"os_predef.h"
#include"os_std.h"
#include"os_gl.h"

#include"Extrude.h"
#include"Base.h"
#include"OOMac.h"
#include"Setting.h"
#include"Feedback.h"

#include "pymol/algorithm.h"

static
void ExtrudeInit(PyMOLGlobals * G, CExtrude * I);

#define CopyArray(dst,src,type,count) memcpy(dst,src,sizeof(type)*(count))

CExtrude *ExtrudeCopyPointsNormalsColors(CExtrude * orig)
{
  int ok = true;
  OOAlloc(orig->G, CExtrude);
  CHECKOK(ok, I);
  if (ok)
    ExtrudeInit(orig->G, I);

  if (ok)
    ok &= ExtrudeAllocPointsNormalsColors(I, orig->N);

  if (ok){
    CopyArray(I->p, orig->p, float, 3 * I->N);
    CopyArray(I->n, orig->n, float, 9 * I->N);
    CopyArray(I->c, orig->c, float, 3 * I->N);
    CopyArray(I->alpha, orig->alpha, float, I->N);
    CopyArray(I->i, orig->i, unsigned int, I->N);
    CopyArray(I->sf, orig->sf, float, I->N);      /* PUTTY: scale factors */
  } else {
    ExtrudeFree(I);
    I = NULL;
  }
  return (I);
}

void ExtrudeInit(PyMOLGlobals * G, CExtrude * I)
{
  I->G = G;

  I->N = 0;
  I->p = NULL;
  I->n = NULL;
  I->c = NULL;
  I->alpha = nullptr;
  I->i = NULL;

  I->sv = NULL;                 /* shape vertices */
  I->sn = NULL;                 /* shape normals */
  I->tv = NULL;                 /* transformed vertices */
  I->tn = NULL;                 /* transformed normals */
  I->Ns = 0;                    /* number of shape points */

  I->sf = NULL;
}

int ExtrudeCircle(CExtrude * I, int n, float size)
{
  int a;
  float *v, *vn;
  int ok = true;
  PRINTFD(I->G, FB_Extrude)
    " ExtrudeCircle-DEBUG: entered.\n" ENDFD;
  /*
  if(n > 50)
  n = 50;*/

  FreeP(I->sv);
  FreeP(I->sn);
  FreeP(I->tv);
  FreeP(I->tn);

  I->sv = pymol::malloc<float>(3 * (n + 1));
  CHECKOK(ok, I->sv);
  if (ok)
    I->sn = pymol::malloc<float>(3 * (n + 1));
  CHECKOK(ok, I->sn);
  if (ok)
    I->tv = pymol::malloc<float>(3 * (n + 1));
  CHECKOK(ok, I->tv);
  if (ok)
    I->tn = pymol::malloc<float>(3 * (n + 1));
  CHECKOK(ok, I->tn);

  if (ok){
    I->Ns = n;
    I->r = size;
    
    v = I->sv;
    vn = I->sn;
    
    for(a = 0; a <= n; a++) {
      *(vn++) = 0.0;
      *(vn++) = (float) cos(a * 2 * PI / n);
      *(vn++) = (float) sin(a * 2 * PI / n);
      *(v++) = 0.0;
      *(v++) = (float) cos(a * 2 * PI / n) * size;
      *(v++) = (float) sin(a * 2 * PI / n) * size;
    }
  }

  if (!ok){
    FreeP(I->sv);
    FreeP(I->sn);
    FreeP(I->tv);
    FreeP(I->tn);
    I->sv = NULL;
    I->sn = NULL;
    I->tv = NULL;
    I->tn = NULL;
  }

  PRINTFD(I->G, FB_Extrude)
    " ExtrudeCircle-DEBUG: exiting...\n" ENDFD;
  return ok;
}

int ExtrudeOval(CExtrude * I, int n, float width, float length)
{
  int a;
  float *v, *vn;
  int ok = true;

  PRINTFD(I->G, FB_Extrude)
    " ExtrudeOval-DEBUG: entered.\n" ENDFD;

  /*  if(n > 50)
      n = 50;*/

  FreeP(I->sv);
  FreeP(I->sn);
  FreeP(I->tv);
  FreeP(I->tn);

  I->sv = pymol::malloc<float>(3 * (n + 1));
  CHECKOK(ok, I->sv);
  if (ok)
    I->sn = pymol::malloc<float>(3 * (n + 1));
  CHECKOK(ok, I->sn);
  if (ok)
    I->tv = pymol::malloc<float>(3 * (n + 1));
  CHECKOK(ok, I->tv);
  if (ok)
    I->tn = pymol::malloc<float>(3 * (n + 1));
  CHECKOK(ok, I->tn);
  I->Ns = n;

  v = I->sv;
  vn = I->sn;

  for(a = 0; a <= n; a++) {
    *(vn++) = 0.0;
    *(vn++) = (float) cos(a * 2 * PI / n) * length;
    *(vn++) = (float) sin(a * 2 * PI / n) * width;
    *(v++) = 0.0;
    *(v++) = (float) cos(a * 2 * PI / n) * width;
    *(v++) = (float) sin(a * 2 * PI / n) * length;
  }

  PRINTFD(I->G, FB_Extrude)
    " ExtrudeOval-DEBUG: exiting...\n" ENDFD;
  if (!ok){
    FreeP(I->sv);
    FreeP(I->sn);
    FreeP(I->tv);
    FreeP(I->tn);
  }
  return ok;
}

int ExtrudeRectangle(CExtrude * I, float width, float length, int mode)
{
  float *v, *vn;
  int ok = true;

  PRINTFD(I->G, FB_Extrude)
    " ExtrudeRectangle-DEBUG: entered...\n" ENDFD;

  switch (mode) {
  case 0:
    I->Ns = 8;
    break;
  default:
    I->Ns = 4;
    break;
  }

  FreeP(I->sv);
  FreeP(I->sn);
  FreeP(I->tv);
  FreeP(I->tn);

  I->sv = pymol::malloc<float>(3 * (I->Ns + 1));
  CHECKOK(ok, I->sv);
  if (ok)
    I->sn = pymol::malloc<float>(3 * (I->Ns + 1));
  CHECKOK(ok, I->sn);
  if (ok)
    I->tv = pymol::malloc<float>(3 * (I->Ns + 1));
  CHECKOK(ok, I->tv);
  if (ok)
    I->tn = pymol::malloc<float>(3 * (I->Ns + 1));
  CHECKOK(ok, I->tn);

  if (!ok){
    FreeP(I->sv);
    FreeP(I->sn);
    FreeP(I->tv);
    FreeP(I->tn);
    I->sv = NULL;
    I->sn = NULL;
    I->tv = NULL;
    I->tn = NULL;
    return ok;
  }

  v = I->sv;
  vn = I->sn;

  if((!mode) || (mode == 1)) {
    *(vn++) = 0.0;
    *(vn++) = 1.0;
    *(vn++) = 0.0;
    *(vn++) = 0.0;
    *(vn++) = 1.0;
    *(vn++) = 0.0;
    *(v++) = 0.0;
    *(v++) = (float) cos(PI / 4) * width;
    *(v++) = (float) -sin(PI / 4) * length;
    *(v++) = 0.0;
    *(v++) = (float) cos(PI / 4) * width;
    *(v++) = (float) sin(PI / 4) * length;
  }

  if((!mode) || (mode == 2)) {
    *(vn++) = 0.0;
    *(vn++) = 0.0;
    *(vn++) = 1.0;
    *(vn++) = 0.0;
    *(vn++) = 0.0;
    *(vn++) = 1.0;
    *(v++) = 0.0;
    *(v++) = (float) cos(PI / 4) * width;
    *(v++) = (float) sin(PI / 4) * length;
    *(v++) = 0.0;
    *(v++) = (float) -cos(PI / 4) * width;
    *(v++) = (float) sin(PI / 4) * length;
  }

  if((!mode) || (mode == 1)) {
    *(vn++) = 0.0;
    *(vn++) = -1.0;
    *(vn++) = 0.0;
    *(vn++) = 0.0;
    *(vn++) = -1.0;
    *(vn++) = 0.0;
    *(v++) = 0.0;
    *(v++) = (float) -cos(PI / 4) * width;
    *(v++) = (float) sin(PI / 4) * length;
    *(v++) = 0.0;
    *(v++) = (float) -cos(PI / 4) * width;
    *(v++) = (float) -sin(PI / 4) * length;
  }

  if((!mode) || (mode == 2)) {

    *(vn++) = 0.0;
    *(vn++) = 0.0;
    *(vn++) = -1.0;
    *(vn++) = 0.0;
    *(vn++) = 0.0;
    *(vn++) = -1.0;
    *(v++) = 0.0;
    *(v++) = (float) -cos(PI / 4) * width;
    *(v++) = (float) -sin(PI / 4) * length;
    *(v++) = 0.0;
    *(v++) = (float) cos(PI / 4) * width;
    *(v++) = (float) -sin(PI / 4) * length;
  }

  PRINTFD(I->G, FB_Extrude)
    " ExtrudeRectangle-DEBUG: exiting...\n" ENDFD;
  return ok;
}

int ExtrudeDumbbell1(CExtrude * I, float width, float length, int mode)
{
  float *v, *vn;
  int ok = true;

  PRINTFD(I->G, FB_Extrude)
    " ExtrudeDumbbell1-DEBUG: entered...\n" ENDFD;

  switch (mode) {
  case 0:
    I->Ns = 4;
    break;
  default:
    I->Ns = 2;
    break;

  }

  FreeP(I->sv);
  FreeP(I->sn);
  FreeP(I->tv);
  FreeP(I->tn);

  I->sv = pymol::malloc<float>(3 * (I->Ns + 1));
  CHECKOK(ok, I->sv);
  if (ok)
    I->sn = pymol::malloc<float>(3 * (I->Ns + 1));
  CHECKOK(ok, I->sn);
  if (ok)
    I->tv = pymol::malloc<float>(3 * (I->Ns + 1));
  CHECKOK(ok, I->tv);
  if (ok)
    I->tn = pymol::malloc<float>(3 * (I->Ns + 1));
  CHECKOK(ok, I->tn);

  if (!ok){
    FreeP(I->sv);
    FreeP(I->sn);
    FreeP(I->tv);
    FreeP(I->tn);
    I->sv = NULL;
    I->sn = NULL;
    I->tv = NULL;
    I->tn = NULL;
  }

  v = I->sv;
  vn = I->sn;

  if((!mode) || (mode == 1)) {  /* top */
    *(vn++) = 0.0;
    *(vn++) = 1.0;
    *(vn++) = 0.0;
    *(vn++) = 0.0;
    *(vn++) = 1.0;
    *(vn++) = 0.0;
    *(v++) = 0.0;
    *(v++) = (float) cos(PI / 4) * width;
    *(v++) = (float) -sin(PI / 4) * length;
    *(v++) = 0.0;
    *(v++) = (float) cos(PI / 4) * width;
    *(v++) = (float) sin(PI / 4) * length;
  }

  if((!mode) || (mode == 2)) {  /* bottom */
    *(vn++) = 0.0;
    *(vn++) = -1.0;
    *(vn++) = 0.0;
    *(vn++) = 0.0;
    *(vn++) = -1.0;
    *(vn++) = 0.0;
    *(v++) = 0.0;
    *(v++) = (float) -cos(PI / 4) * width;
    *(v++) = (float) sin(PI / 4) * length;
    *(v++) = 0.0;
    *(v++) = (float) -cos(PI / 4) * width;
    *(v++) = (float) -sin(PI / 4) * length;
  }

  PRINTFD(I->G, FB_Extrude)
    " ExtrudeDumbbell1-DEBUG: exiting...\n" ENDFD;
  return ok;
}

void ExtrudeDumbbellEdge(CExtrude * I, int samp, int sign, float length)
{
  int a;
  float *n, *p, f, disp;

  PRINTFD(I->G, FB_Extrude)
    " ExtrudeDumbbellEdge-DEBUG: entered.\n" ENDFD;
  disp = (float) (sign * sin(PI / 4) * length);
  p = I->p;
  n = I->n;
  for(a = 0; a < I->N; a++) {
    if(a <= samp)
      f = disp * smooth((a / ((float) samp)), 2);
    else if(a >= (I->N - samp))
      f = disp * smooth(((I->N - a - 1) / ((float) samp)), 2);
    else
      f = disp;
    n += 6;
    (*p++) += *(n++) * f;
    (*p++) += *(n++) * f;
    (*p++) += *(n++) * f;
  }
  PRINTFD(I->G, FB_Extrude)
    " ExtrudeDumbbellEdge-DEBUG: exiting...\n" ENDFD;

}

#if 0
int ExtrudeDumbbell2(CExtrude * I, int n, int sign, float length, float size)
{
  int a;
  float *v, *vn;
  int ok = true;

  PRINTFD(I->G, FB_Extrude)
    " ExtrudeDumbbell2-DEBUG: entered.\n" ENDFD;
  /*  if(n > 50)
      n = 50;*/

  FreeP(I->sv);
  FreeP(I->sn);
  FreeP(I->tv);
  FreeP(I->tn);
  
  I->sv = pymol::malloc<float>(3 * (n + 1));
  CHECKOK(ok, I->sv);
  if (ok)
    I->sn = pymol::malloc<float>(3 * (n + 1));
  CHECKOK(ok, I->sn);
  if (ok)
    I->tv = pymol::malloc<float>(3 * (n + 1));
  CHECKOK(ok, I->tv);
  if (ok)
    I->tn = pymol::malloc<float>(3 * (n + 1));
  CHECKOK(ok, I->tn);

  if (!ok){
    FreeP(I->sv);
    FreeP(I->sn);
    FreeP(I->tv);
    FreeP(I->tn);
    I->sv = NULL;
    I->sn = NULL;
    I->tv = NULL;
    I->tn = NULL;
  }

  I->Ns = n;

  v = I->sv;
  vn = I->sn;

  for(a = 0; a <= n; a++) {
    *(vn++) = 0.0;
    *(vn++) = (float) cos(a * 2 * PI / n);
    *(vn++) = (float) sin(a * 2 * PI / n);
    *(v++) = 0.0;
    *(v++) = (float) cos(a * 2 * PI / n) * size;
    *(v++) = (float) ((sin(a * 2 * PI / n) * size) + (sign * sin(PI / 4) * length));
  }

  PRINTFD(I->G, FB_Extrude)
    " ExtrudeDumbbell2-DEBUG: exiting...\n" ENDFD;
  return ok;
}
#endif

CExtrude *ExtrudeNew(PyMOLGlobals * G)
{
  int ok = true;
  OOAlloc(G, CExtrude);
  CHECKOK(ok, I);
  if (ok)
    ExtrudeInit(G, I);
  return (I);
}

void ExtrudeBuildNormals1f(CExtrude * I)
{
  int a;
  float *v;

  PRINTFD(I->G, FB_Extrude)
    " ExtrudeBuildNormals1f-DEBUG: entered.\n" ENDFD;

  if(I->N) {
    get_system1f3f(I->n, I->n + 3, I->n + 6);   /* first is arbitrary */
    v = I->n + 9;
    for(a = 1; a < I->N; a++) {
      copy3f(v - 6, v + 3);
      get_system2f3f(v, v + 3, v + 6);  /* the rest are relative to first */
      v += 9;
    }
  }

  PRINTFD(I->G, FB_Extrude)
    " ExtrudeBuildNormals1f-DEBUG: exiting...\n" ENDFD;

}

void ExtrudeBuildNormals2f(CExtrude * I)
{
  int a;
  float *v;
  PRINTFD(I->G, FB_Extrude)
    " ExtrudeBuildNormals2f-DEBUG: entered.\n" ENDFD;

  if(I->N) {
    v = I->n;
    for(a = 0; a < I->N; a++) {
      get_system2f3f(v, v + 3, v + 6);
      v += 9;
    }
  }

  PRINTFD(I->G, FB_Extrude)
    " ExtrudeBuildNormals2f-DEBUG: entering...\n" ENDFD;

}

/**
 * Approximate a helix center trace.
 *
 * - shifts points to the helix center
 * - smoothes the points curve (low-pass filter)
 *
 * @pre I has helix geometry
 * @post I has linear geometry along helix axis
 * @post I has normals
 *
 * @param[in,out] I data structure to modify in-place.
 * @param radius Cylindrical helix radius (for optimizing end points)
 * @param sampling Samples per residue
 */
void ExtrudeShiftToAxis(CExtrude* I, float radius, int sampling)
{
  assert(I->N > 1);

  constexpr float loop_radius = 0.2f;
  const int smooth_cycles =
      SettingGet<int>(I->G, cSetting_cartoon_smooth_cylinder_cycles);
  const int smooth_window =
      SettingGet<int>(I->G, cSetting_cartoon_smooth_cylinder_window);

  float p_start[3];
  float p_end[3];

  // original start and end points
  copy3(I->p, p_start);
  copy3(I->p + (I->N - 1) * 3, p_end);

  ExtrudeBuildNormals2f(I);

  // Because segments have open ends, we don't have ideal helix normals for the
  // first and last position. We can reconstruct the desired normals from the
  // next position by an ideal rotation.
  if (I->N > 2) {
#if 0
    // dump residue rotations to figure out ideal rotation
    for (int a = 1; a + 2 < I->N; ++a) {
      const float* base0 = I->n + (a + 0) * 9;
      const float* base1 = I->n + (a + 1) * 9;

      float base0_inv[9];
      float residue_rotation[9];

      assert(fabs(determinant33f(base0) - 1.0) < 1e-3);
      transpose33f33f(base0, base0_inv);
      multiply33f33f(base1, base0_inv, residue_rotation);

      printf("========= a %d\n", a);
      dump33f(residue_rotation, "rotation");
    }
#endif

    // ideal rotation from one alpha helix residue to the next
    static float const residue_rotation[9] = {
        0.224, -0.809, -0.544, 0.809, -0.157, 0.567, -0.544, -0.567, 0.619};
    static float const residue_rotation_inv[9] = {
        0.224, 0.809, -0.544, -0.809, -0.157, -0.567, -0.544, 0.567, 0.619};

    multiply33f33f(residue_rotation_inv, //
        I->n + 9 * sampling, I->n);
    multiply33f33f(residue_rotation, //
        I->n + 9 * (I->N - 1 - sampling), I->n + 9 * (I->N - 1));
  }

  // move points to helix axes
  for (int a = 0; a < I->N; ++a) {
    float* point = I->p + a * 3;
    const float* normal = I->n + a * 9 + 3;

    // distance to move point (distance between C-alpha atom and helix axis)
    float factor = 2.3;

    // keep end points close enough to adjacent segments so they overlap
    if (a == 0 || a + 1 == I->N) {
      factor = std::min(factor, radius - loop_radius);
    }

    float tmp[3];
    scale3f(normal, -factor, tmp);
    add3f(point, tmp, point);
  }

  // window averaging of point positions
  if (I->N > 2 && smooth_window > 0) {
    int const w2 = sampling * smooth_window;

    for (int i = 0; i < smooth_cycles; ++i) {
      std::vector<float> smoothed((I->N - 2) * 3);

      for (int a = 1; a + 1 < I->N; ++a) {
        float* avg = smoothed.data() + (a - 1) * 3;

        for (int j = -w2; j <= w2; ++j) {
          int const k = pymol::clamp(a + j, 0, I->N - 1);
          add3f(I->p + k * 3, avg, avg);
        }

        scale3f(avg, 1. / (2 * w2 + 1), avg);
      }

      std::copy(smoothed.begin(), smoothed.end(), I->p + 3);
    }
  }

  ExtrudeComputeTangents(I);
  ExtrudeBuildNormals1f(I);

  // extend tips for better geometry overlap with adjacent segments
  auto const push_out_tip = [&](int index, float const* pos_orig,
                                int direction) {
    float const protrusion = loop_radius * 2;
    float const* normal = I->n + index * 9;
    float* pos = I->p + index * 3;
    float offset[3];
    subtract3f(pos_orig, pos, offset);
    float const len = project3f(offset, normal, offset) * direction;
    if (len > -protrusion) {
      scale3f(normal, (len + protrusion) * direction, offset);
      add3f(pos, offset, pos);
    }
  };

  push_out_tip(0, p_start, -1);
  push_out_tip(I->N - 1, p_end, 1);
}

#if 0
/* Is this ever used? */
void ExtrudeCGOTraceAxes(CExtrude * I, CGO * cgo)
{
  int a;
  float *v, *n;
  float v0[3];

  if(I->N) {
    CGOColor(cgo, 0.5, 0.5, 0.5);
    CGOBegin(cgo, GL_LINES);
    v = I->p;
    n = I->n;
    for(a = 0; a < I->N; a++) {
      add3f(v, n, v0);
      CGOVertexv(cgo, v0);
      CGOVertexv(cgo, v);
      n += 3;
      add3f(v, n, v0);
      CGOVertexv(cgo, v0);
      CGOVertexv(cgo, v);
      n += 3;
      add3f(v, n, v0);
      CGOVertexv(cgo, v0);
      CGOVertexv(cgo, v);
      n += 3;
      v += 3;
    }
    CGOEnd(cgo);
  }
}

/* Is this ever used? */
void ExtrudeCGOTrace(CExtrude * I, CGO * cgo)
{
  int a;
  float *v;
  if(I->N) {
    CGOColor(cgo, 0.5, 0.5, 0.5);
    CGOBegin(cgo, GL_LINE_STRIP);
    v = I->p;
    for(a = 0; a < I->N; a++) {
      CGOVertexv(cgo, v);
      v += 3;
    }
    CGOEnd(cgo);
  }
}
#endif

int ExtrudeComputeTangents(CExtrude * I)
{
  float *nv, *v1, *v;
  int a;
  int ok = true;

  PRINTFD(I->G, FB_Extrude)
    " ExtrudeComputeTangents-DEBUG: entered.\n" ENDFD;

  nv = pymol::malloc<float>(I->N * 3);
  CHECKOK(ok, nv);
  if (!ok)
    return ok;

  v = nv;
  v1 = I->p + 3;

  for(a = 1; a < I->N; a++) {
    subtract3f(v1, v1 - 3, v);
    normalize3f(v);
    v += 3;
    v1 += 3;
  }

  /* compute tangents */

  v = nv;
  v1 = I->n;

  *(v1++) = *(v++);             /* first segment */
  *(v1++) = *(v++);
  *(v1++) = *(v++);
  v1 += 6;

  for(a = 1; a < (I->N - 1); a++) {

    add3f(v, (v - 3), v1);
    normalize3f(v1);
    v1 += 9;
    v += 3;
  }

  *(v1++) = *(v - 3);           /* last segment */
  *(v1++) = *(v - 2);
  *(v1++) = *(v - 1);

  FreeP(nv);

  PRINTFD(I->G, FB_Extrude)
    " ExtrudeComputeTangents-DEBUG: exiting...\n" ENDFD;
  return ok;
}

#if 0
/* Is this ever used? */
void ExtrudeCGOTraceFrame(CExtrude * I, CGO * cgo)
{
  int a, b;
  float *v;
  float *n;
  float *sv, *tv;
  float v0[3], v1[3];

  if(I->N && I->Ns) {
    CGOColor(cgo, 0.5, 0.5, 0.5);
    {
      CGOBegin(cgo, GL_LINES);
      v = I->p;
      n = I->n;
      for(a = 0; a < I->N; a++) {
	sv = I->sv;
	tv = I->tv;
	for(b = 0; b < I->Ns; b++) {
	  transform33Tf3f(n, sv, tv);
	  sv += 3;
	  tv += 3;
	}
	/* trace shape */
	tv = I->tv;
	add3f(v, tv, v0);
	for(b = 1; b < I->Ns; b++) {
	  tv += 3;
	  add3f(v, tv, v1);
	  CGOVertexv(cgo, v0);
	  CGOVertexv(cgo, v1);
	  copy3f(v1, v0);
	}
	tv = I->tv;
	add3f(v, tv, v1);
	CGOVertexv(cgo, v0);
	CGOVertexv(cgo, v1);
	v += 3;
	n += 9;
      }
      CGOEnd(cgo);
    }
  }
}
#endif

/**
 * Draw flat cap on a tube cartoon (loop, oval, etc.)
 *
 * I: tube instance
 * cgo: CGO to add to
 * index: sampling index in `I`
 * inv_dir: inverse direction of normal if true
 * color: RGB color or NULL to use color from `I`
 */
static
void TubeCapFlat(const CExtrude * I, CGO * cgo, int index, bool inv_dir, const float * color) {
  const float * vertex = I->p + 3 * index;
  const float * base33 = I->n + 9 * index;
  const float * normal = base33;
  float tmp3f[3];
  int b_end = -1, b_incr = -1;

  if (inv_dir) {
    copy3f(normal, tmp3f);
    invert3f(tmp3f);
    normal = tmp3f;
  } else {
    b_end = I->Ns * 2 + 1;
    b_incr = 1;
  }

  CGOBegin(cgo, GL_TRIANGLE_FAN);
  CGOColorv(cgo, color ? color : (I->c + 3 * index));
  CGOAlpha(cgo, I->alpha[index]);
  CGOPickColor(cgo, I->i[index], cPickableAtom);
  CGONormalv(cgo, normal); // (tmp3f again free to use)
  CGOVertexv(cgo, vertex); // center

  // Indexing trickery: going in a loop, visiting first index twice to
  // close the loop. Iterating backwards in case of `inv_dir`.
  for (int b = I->Ns; b != b_end; b += b_incr) {
    transform33Tf3f(base33, I->sv + (b % I->Ns) * 3, tmp3f);
    add3f(vertex, tmp3f, tmp3f);
    CGOVertexv(cgo, tmp3f);
  }

  CGOEnd(cgo);
  CGOPickColor(cgo, -1, cPickableNoPick);
}

/**
 * I: tube instance
 * cgo: CGO to add to
 * cap: 0: no caps, 1: flat caps, 2: round caps
 * color_override: RGB color or NULL to use color from `I`
 * use_spheres: do round caps with spheres instead of triangles
 * dash: if > 0, skip every segment which is a multiple of `dash`
 */
int ExtrudeCGOSurfaceTube(const CExtrude* I, CGO* cgo, cCylCap cap,
    const float* color_override, bool use_spheres, int dash)
{
  int a, b;
  unsigned int *i;
  float *v;
  float *n;
  float *c;
  const float *alpha;
  float *sv, *sn, *tv, *tn, *tv1, *tn1, *TV = NULL, *TN = NULL;
  int start, stop;
  int ok = true;
  PRINTFD(I->G, FB_Extrude)
    " ExtrudeCGOSurfaceTube-DEBUG: entered.\n" ENDFD;

  if(I->N && I->Ns) {
    TV = pymol::malloc<float>(3 * (I->Ns + 1) * I->N);
    CHECKOK(ok, TV);
    if (ok)
      TN = pymol::malloc<float>(3 * (I->Ns + 1) * I->N);
    CHECKOK(ok, TN);
    /* compute transformed shape vertices */

    if (ok){
      tn = TN;
      tv = TV;
      
      sv = I->sv;
      sn = I->sn;
      for(b = 0; b <= I->Ns; b++) {
	if(b == I->Ns) {
	  sv = I->sv;
	  sn = I->sn;
	}
	v = I->p;
	n = I->n;
	
	for(a = 0; a < I->N; a++) {
	  transform33Tf3f(n, sv, tv);
	  add3f(v, tv, tv);
	  tv += 3;
	  transform33Tf3f(n, sn, tn);
	  tn += 3;
	  n += 9;
	  v += 3;
	}
	sv += 3;
	sn += 3;
      }
      
      start = I->Ns / 4;
      stop = 3 * I->Ns / 4;
    }

    // first loop: axial, for dashes (skipping segments)
    for (int a_start = 0, a_incr = (dash ? dash : I->N);
        a_start < I->N - 1;
        a_start += a_incr) {

      int a_end = a_start + a_incr;
      if (a_end > I->N)
        a_end = I->N;

      // second loop: circumferential, setting up axial triangle strips
      for(b = 0; ok && b < I->Ns; b++) {
	if(SettingGetGlobal_i(I->G, cSetting_cartoon_debug) < 1.5)
	  ok &= CGOBegin(cgo, GL_TRIANGLE_STRIP);
	else {
	  ok &= CGOBegin(cgo, GL_LINE_STRIP);
	}
	if (ok){
          c = I->c + a_start * 3;
          alpha = I->alpha + a_start;
          i = I->i + a_start;

          tv = TV + 3 * (a_start + b * I->N);
          tn = TN + 3 * (a_start + b * I->N);

          tv1 = tv + 3 * I->N;
          tn1 = tn + 3 * I->N;

          // third loop: axial, segments within one "dash"
          for(a = a_start; ok && a < a_end; ++a) {
	    if(color_override && (b > start) && (b < stop))
	      ok &= CGOColorv(cgo, color_override);
	    else
	      ok &= CGOColorv(cgo, c);
        if (ok){
          ok &= CGOAlpha(cgo, *alpha);
        }
	    if (ok)
	      ok &= CGOPickColor(cgo, *i, cPickableAtom);
	    if (ok)
	      ok &= CGONormalv(cgo, tn);
	    if (ok)
	      ok &= CGOVertexv(cgo, tv);
	    tn += 3;
	    tv += 3;
	    if (ok)
	      ok &= CGONormalv(cgo, tn1);
	    if (ok)
	      ok &= CGOVertexv(cgo, tv1);
	    tn1 += 3;
	    tv1 += 3;
	    c += 3;
        alpha++;
	    i++;
	  }
	}
	if (ok)
	  ok &= CGOEnd(cgo);
	if (ok)
	  ok &= CGOPickColor(cgo, -1, cPickableNoPick);
      }

      if (cap == cCylCap::Flat) {
        TubeCapFlat(I, cgo, a_start, true, color_override);
        TubeCapFlat(I, cgo, a_end - 1, false, color_override);
      }
    }

    if (ok){
    switch (cap) {
    case cCylCap::Round:
      {
	float p0[3], p1[3], p2[3], z1, z2, normal[3], vertex1[3];
	float c, d, prev, x, y, *v1, nEdge = I->Ns, nEdgeH = 2.f * floor(I->Ns/2.f);

	if (ok){
	  n = I->n;
	  v = I->p;
	  sv = I->sv;
	  tv = I->tv;
	  for(b = 0; b < I->Ns; b++) {
	    transform33Tf3f(n, sv, tv);
	    add3f(v, tv, tv);
	    sv += 3;
	    tv += 3;
	  }
	  
	  copy3f(I->n, p0);
	  invert3f(p0);
	  transform33Tf3f(I->n, I->sv, p1);
	  cross_product3f(p0, p1, p2);
	  normalize3f(p1);
	  normalize3f(p2);
	}
	if (ok){
	  if(color_override)
	    ok &= CGOColorv(cgo, color_override);
	  else
	    ok &= CGOColorv(cgo, I->c);
	}
    if (ok){
      ok &= CGOAlpha(cgo, I->alpha[0]);
    }
	if (ok)
	  ok &= CGOPickColor(cgo, I->i[0], cPickableAtom);
	if (ok){
        v = I->p;
	if (use_spheres){
	  ok &= CGOSphere(cgo, v, I->r); // this matches the Cylinder
	} else {
	  /* If we don't use spheres, then we need to have the rounded cap
	   * line up with the geometry perfectly.  We generate the cap using 
	   * spheracle coordinates, then for the last line (i.e., last=true)
	   * we use the coordinates from the geometry (i.e., exactly 
	   * how they were genereated from I->n and I->sv).  This is so that 
	   * the vertices line up perfectly. */
	  nEdge = I->Ns;
	  ok &= CGOBegin(cgo, GL_TRIANGLE_STRIP);
	  z1 = z2 = 1.f;
	  prev = 1.f;
	  v1 = v;
	  tv = I->tv;
	  for (c = 1; ok && c <= (nEdgeH/2); c++){
	    short last = (c + 1) > (nEdgeH/2);
	    z1 = z2;
	    z2 = (float) cos((c) * PI / ((float)nEdgeH));
	    for (d = 0; d <= nEdge; d++){
	      x = (float) cos((d) * 2 * PI / (float)nEdge) * sin((c-prev) * PI / (float)nEdgeH);
	      y = (float) sin((d) * 2 * PI / (float)nEdge) * sin((c-prev) * PI / (float)nEdgeH);
	      normal[0] = p1[0] * x + p2[0] * y + p0[0] * z1;
	      normal[1] = p1[1] * x + p2[1] * y + p0[1] * z1;
	      normal[2] = p1[2] * x + p2[2] * y + p0[2] * z1;
	      vertex1[0] = v1[0] + normal[0] * I->r;
	      vertex1[1] = v1[1] + normal[1] * I->r;
	      vertex1[2] = v1[2] + normal[2] * I->r;
	      normalize3f(normal);
	      ok &= CGONormalv(cgo, normal);	  
	      if (ok)
		ok &= CGOVertexv(cgo, vertex1);
	      
	      if (ok){
	      if (last){
		float *vert = tv + 3*(((int)(nEdge-d))%((int)nEdge));
		subtract3f(vert, v, normal);
		ok &= CGONormalv(cgo, normal);
		if (ok)
		  ok &= CGOVertexv(cgo, vert);
	      } else {
		x = (float) cos((d) * 2 * PI / (float)nEdge) * sin((c) * PI / (float)nEdgeH);
		y = (float) sin((d) * 2 * PI / (float)nEdge) * sin((c) * PI / (float)nEdgeH);
		normal[0] = p1[0] * x + p2[0] * y + p0[0] * z2;
		normal[1] = p1[1] * x + p2[1] * y + p0[1] * z2;
		normal[2] = p1[2] * x + p2[2] * y + p0[2] * z2;
		vertex1[0] = v1[0] + normal[0] * I->r;
		vertex1[1] = v1[1] + normal[1] * I->r;
		vertex1[2] = v1[2] + normal[2] * I->r;
		normalize3f(normal);
		ok &= CGONormalv(cgo, normal);	  
		if (ok)
		  ok &= CGOVertexv(cgo, vertex1);
	      }
	      }
	    }
	  }
	  if (ok)
	    ok &= CGOEnd(cgo);
	  if (ok)
	    ok &= CGOPickColor(cgo, -1, cPickableNoPick);
	}

	if (ok){
	  n = I->n + 9 * (I->N - 1);
	  v = I->p + 3 * (I->N - 1);
	  sv = I->sv;
	  tv = I->tv;
	  for(b = 0; b < I->Ns; b++) {
	    transform33Tf3f(n, sv, tv);
	    add3f(v, tv, tv);
	    sv += 3;
	    tv += 3;
	  }
	  
	  copy3f(n, p0);
	  transform33Tf3f(n, I->sv, p1);
	  cross_product3f(p0, p1, p2);
	  normalize3f(p1);
	  normalize3f(p2);
	  
	  if(color_override)
	    ok &= CGOColorv(cgo, color_override);
	  else
	    ok &= CGOColorv(cgo, I->c + 3 * (I->N - 1));
	}
    if (ok){
      ok &= CGOAlpha(cgo, I->alpha[I->N - 1]);
    }
	if (ok)
	  ok &= CGOPickColor(cgo, I->i[I->N - 1], cPickableAtom);

	if (ok){
	if (use_spheres){
	  ok &= CGOSphere(cgo, v, I->r); // this matches the Cylinder
	} else {
	  /* If we don't use spheres, then we need to have the rounded cap
	   * line up with the geometry perfectly.  We generate the cap using 
	   * spheracle coordinates, then for the last line (i.e., last=true)
	   * we use the coordinates from the geometry (i.e., exactly 
	   * how they were genereated from I->n and I->sv).  This is so that 
	   * the vertices line up perfectly. */
	  nEdge = I->Ns;
	  ok &= CGOBegin(cgo, GL_TRIANGLE_STRIP);
	  z1 = z2 = 1.f;
	  prev = 1.f;
	  v1 = v;
	  tv = I->tv;
	  for (c = 1; ok && c <= (nEdgeH/2); c++){
	    short last = (c + 1) > (nEdgeH/2);
	    z1 = z2;
	    z2 = (float) cos((c) * PI / ((float)nEdgeH));
	    for (d = 0; d <= nEdge; d++){
	      x = (float) cos((d) * 2 * PI / (float)nEdge) * sin((c-prev) * PI / (float)nEdgeH);
	      y = (float) sin((d) * 2 * PI / (float)nEdge) * sin((c-prev) * PI / (float)nEdgeH);
	      normal[0] = p1[0] * x + p2[0] * y + p0[0] * z1;
	      normal[1] = p1[1] * x + p2[1] * y + p0[1] * z1;
	      normal[2] = p1[2] * x + p2[2] * y + p0[2] * z1;
	      vertex1[0] = v1[0] + normal[0] * I->r;
	      vertex1[1] = v1[1] + normal[1] * I->r;
	      vertex1[2] = v1[2] + normal[2] * I->r;
	      normalize3f(normal);
	      ok &= CGONormalv(cgo, normal);	  
	      if (ok)
		ok &= CGOVertexv(cgo, vertex1);
	      
	      if (ok){
	      if (last){
		float *vert = tv + 3*(((int)(d))%((int)nEdge));
		subtract3f(vert, v, normal);
		ok &= CGONormalv(cgo, normal);
		if (ok)
		  ok &= CGOVertexv(cgo, vert);
	      } else {
		x = (float) cos((d) * 2 * PI / (float)nEdge) * sin((c) * PI / (float)nEdgeH);
		y = (float) sin((d) * 2 * PI / (float)nEdge) * sin((c) * PI / (float)nEdgeH);
		normal[0] = p1[0] * x + p2[0] * y + p0[0] * z2;
		normal[1] = p1[1] * x + p2[1] * y + p0[1] * z2;
		normal[2] = p1[2] * x + p2[2] * y + p0[2] * z2;
		vertex1[0] = v1[0] + normal[0] * I->r;
		vertex1[1] = v1[1] + normal[1] * I->r;
		vertex1[2] = v1[2] + normal[2] * I->r;
		normalize3f(normal);
		ok &= CGONormalv(cgo, normal);	  
		if (ok)
		  ok &= CGOVertexv(cgo, vertex1);
	      }
	      }
	    }
	  }
	  if (ok)
	    ok &= CGOEnd(cgo);
	  if (ok)
	    ok &= CGOPickColor(cgo, -1, cPickableNoPick);
	}
	}
	}
      }
      break;
    }
    }
    FreeP(TV);
    FreeP(TN);
  }

  PRINTFD(I->G, FB_Extrude)
    " ExtrudeCGOSurfaceTube-DEBUG: exiting...\n" ENDFD;
  return ok;
}

int ExtrudeCylindersToCGO(CExtrude * I, CGO * cgo, float tube_radius){
  float *v1, *c1, midc[3], axis[3];
  const float *alpha;
  int a;
  unsigned int *i;
  int ok = true;

  PRINTFD(I->G, FB_Extrude)
    " ExtrudeCylindersToCGO-DEBUG: entered.\n" ENDFD;

  v1 = I->p + 3;
  c1 = I->c + 3;
  alpha = I->alpha + 1;
  i = I->i + 1;
  
  int cap = (cCylShaderBothCapsRound | cCylShaderInterpColor);
  for(a = 1; a < I->N; a++) {
    average3f(c1-3, c1, midc);
    ok &= CGOPickColor(cgo, *(i-1), cPickableAtom);
    subtract3f(v1, v1-3, axis);
    CGOColorv(cgo, c1-3);
    CGOAlpha(cgo, *(alpha-1));
    Pickable pickcolor2 = { *i, cPickableAtom };
    cgo->add<cgo::draw::shadercylinder2ndcolor>(cgo, v1-3, axis, tube_radius, cap, c1, &pickcolor2);
    v1 += 3;
    c1 += 3;
    alpha++;
    i++;
    cap = cCylShaderCap2Round | cCylShaderInterpColor;
  }
  if (ok)
    ok &= CGOPickColor(cgo, 0, cPickableNoPick);

  PRINTFD(I->G, FB_Extrude)
    " ExtrudeCylindersToCGO-DEBUG: exiting...\n" ENDFD;
  return ok;
}

int ExtrudeCGOSurfaceVariableTube(const CExtrude* I, CGO* cgo, cCylCap cap)
{
  int a, b;
  unsigned int *i;
  float *v;
  float *n;
  float *c;
  const float *alpha;
  float *sv, *sn, *tv, *tn, *tv1, *tn1, *TV = NULL, *TN = NULL, *AN = NULL, *an;
  float v0[3];
  float *sf;                    /* PUTTY: scale factor from ExtrudeMakeSausLUT() */
  int ok = true;
  PRINTFD(I->G, FB_Extrude)
    " ExtrudeCGOSurfaceTube-DEBUG: entered.\n" ENDFD;

  if(I->N && I->Ns) {

    TV = pymol::malloc<float>(3 * (I->Ns + 1) * I->N);
    TN = pymol::malloc<float>(3 * (I->Ns + 1) * I->N);
    AN = pymol::malloc<float>(3 * I->N);        /* normals adjusted for changing widths */

    /* compute transformed shape vertices */

    tv = TV;

    sv = I->sv;
    for(b = 0; b <= I->Ns; b++) {
      if(b == I->Ns) {
        sv = I->sv;
      }

      n = I->n;                 /* NOTE: n is not a counter -- it's a 3x3 coordinate system! */
      v = I->p;
      sf = I->sf;               /* PUTTY: scale factors */

      for(a = 0; a < I->N; a++) {
        transform33Tf3f(n, sv, tv);

        *(tv) *= *sf;
        *(tv + 1) *= *sf;
        *(tv + 2) *= *sf;

        add3f(v, tv, tv);
        tv += 3;
        v += 3;
        sf++;
        n += 9;
      }
      sv += 3;
    }

    /* compute transformed normals, taking into account changing radii */

    tn = TN;
    tv = TV;

    sn = I->sn;
    for(b = 0; b <= I->Ns; b++) {

      float d1, d2, r0, r1, r2, x1, x2;

      if(b == I->Ns) {
        sn = I->sn;
      }

      an = AN;
      v = I->p;

      for(a = 0; a < I->N; a++) {
        if((a > 0) && (a < (I->N - 1))) {
          /* compute rises */

          r0 = (float) diff3f(v, tv);
          r1 = (float) (diff3f(v - 3, tv - 3) - r0);
          r2 = (float) (diff3f(v + 3, tv + 3) - r0);

          /* compute runs */

          d1 = (float) diff3f(v - 3, v);
          d2 = (float) diff3f(v + 3, v);

          /* compute x-to-yz weights */

          x1 = r1 / d1;
          x2 = -r2 / d2;

          if(a == 1) {
            an[-3] = x1;
            an[-2] = sn[1];
            an[-1] = sn[2];
            normalize3f(an - 3);
          } else if(a == I->N - 2) {
            an[3] = x2;
            an[4] = sn[1];
            an[5] = sn[2];
            normalize3f(an + 3);
          }
          an[0] = (x1 + x2) / 2.0F;
          an[1] = sn[1];
          an[2] = sn[2];
          normalize3f(an);
        }
        tv += 3;
        v += 3;
        an += 3;
      }

      n = I->n;                 /* NOTE: n is not a counter -- it's a 3x3 coordinate system! */
      an = AN;

      for(a = 0; a < I->N; a++) {
        transform33Tf3f(n, an, tn);
        tn += 3;
        an += 3;
        n += 9;
      }
      sn += 3;
    }

    /* fill in each strip separately */

    tv = TV;
    tn = TN;

    tv1 = TV + 3 * I->N;
    tn1 = TN + 3 * I->N;

    for(b = 0; b < I->Ns; b++) {
      if(SettingGetGlobal_i(I->G, cSetting_cartoon_debug) < 1.5)
        CGOBegin(cgo, GL_TRIANGLE_STRIP);
      else {
        CGOBegin(cgo, GL_LINE_STRIP);
      }
      c = I->c;
      alpha = I->alpha;
      i = I->i;
      for(a = 0; a < I->N; a++) {
        CGOColorv(cgo, c);
        CGOAlpha(cgo, *alpha);
        CGOPickColor(cgo, *i, cPickableAtom);
        CGONormalv(cgo, tn);
        CGOVertexv(cgo, tv);
        tn += 3;
        tv += 3;
        CGONormalv(cgo, tn1);
        CGOVertexv(cgo, tv1);
        tn1 += 3;
        tv1 += 3;
        c += 3;
        alpha++;
        i++;
      }
      CGOEnd(cgo);
      CGOPickColor(cgo, -1, cPickableNoPick);
    }

    if(ok && SettingGetGlobal_i(I->G, cSetting_cartoon_debug) > 3.5) {

      tv = TV;
      tn = TN;

      tv1 = TV + 3 * I->N;
      tn1 = TN + 3 * I->N;

      for(b = 0; b < I->Ns; b++) {
        float vv[3];
	  CGOBegin(cgo, GL_LINES);
	  c = I->c;
      alpha = I->alpha;
	  i = I->i;
	  for(a = 0; a < I->N; a++) {
	    CGOColorv(cgo, c);
        CGOAlpha(cgo, *alpha);
	    copy3f(tn, vv);
	    scale3f(vv, 0.3F, vv);
	    add3f(vv, tv, vv);
	    CGONormalv(cgo, tn);
	    CGOVertexv(cgo, tv);
	    CGOVertexv(cgo, vv);
	    tn += 3;
	    tv += 3;
	    copy3f(tn1, vv);
	    scale3f(vv, 0.3F, vv);
	    add3f(vv, tv1, vv);
	    CGONormalv(cgo, tn1);
	    CGOVertexv(cgo, tv1);
	    CGOVertexv(cgo, vv);
	    tn1 += 3;
	    tv1 += 3;
	    c += 3;
        alpha++;
	    i++;
	  }
	  CGOEnd(cgo);
      }
    }

    if(ok && cap != cCylCap::None) {

      n = I->n;
      v = I->p;
      sf = I->sf;

      sv = I->sv;
      tv = I->tv;
      for(b = 0; b < I->Ns; b++) {
        transform33Tf3f(n, sv, tv);

        *(tv) *= *sf;
        *(tv + 1) *= *sf;
        *(tv + 2) *= *sf;

        add3f(v, tv, tv);
        sv += 3;
        tv += 3;
      }
      CGOBegin(cgo, GL_TRIANGLE_FAN);
      copy3f(I->n, v0);
      invert3f(v0);
      CGOColorv(cgo, I->c);
      CGOAlpha(cgo, I->alpha[0]);
      CGOPickColor(cgo, I->i[0], cPickableAtom);
      CGONormalv(cgo, v0);
      if (ok) {
	CGOVertexv(cgo, v);
	/* trace shape */
	CGOVertexv(cgo, I->tv);
	for(b = I->Ns - 1; b >= 0; b--) {
	  CGOVertexv(cgo, I->tv + b * 3);
	}
	CGOEnd(cgo);
	n = I->n + 9 * (I->N - 1);
	v = I->p + 3 * (I->N - 1);
	sf = I->sf + (I->N - 1);  /* PUTTY */
	
	sv = I->sv;
	tv = I->tv;
	for(b = 0; b < I->Ns; b++) {
	  transform33Tf3f(n, sv, tv);
	  
	  *(tv) *= *(sf);
	  *(tv + 1) *= *(sf);
	  *(tv + 2) *= *(sf);
	  
	  add3f(v, tv, tv);
	  sv += 3;
	  tv += 3;
	}
	
	CGOBegin(cgo, GL_TRIANGLE_FAN);
	CGOColorv(cgo, I->c + 3 * (I->N - 1));
    CGOAlpha(cgo, I->alpha[I->N - 1]);
	CGOPickColor(cgo, I->i[I->N - 1], cPickableAtom);
	CGONormalv(cgo, n);
	CGOVertexv(cgo, v);
	/* trace shape */
	for(b = 0; b < I->Ns; b++) {
	  CGOVertexv(cgo, I->tv + b * 3);
	}
	CGOVertexv(cgo, I->tv);
	CGOEnd(cgo);
      }
      CGOPickColor(cgo, -1, cPickableNoPick);
      FreeP(TV);
      FreeP(TN);
      FreeP(AN);
    }
    
    PRINTFD(I->G, FB_Extrude)
      " ExtrudeCGOSurfaceTube-DEBUG: exiting...\n" ENDFD;
  }
  return ok;
}

int ExtrudeCGOSurfacePolygon(const CExtrude * I, CGO * cgo, cCylCap cap, const float *color_override)
{
  int a, b;
  unsigned int *i;
  float *v;
  float *n;
  float *c;
  const float *alpha;
  float *sv, *sn, *tv, *tn, *tv1, *tn1, *TV = NULL, *TN = NULL;
  float v0[3];
  int ok = true;

  PRINTFD(I->G, FB_Extrude)
    " ExtrudeCGOSurfacePolygon-DEBUG: entered.\n" ENDFD;

  if(I->N && I->Ns) {

    TV = pymol::malloc<float>(3 * (I->Ns + 1) * I->N);
    CHECKOK(ok, TV);
    if (ok)
      TN = pymol::malloc<float>(3 * (I->Ns + 1) * I->N);
    CHECKOK(ok, TN);
    /* compute transformed shape vertices */

    if (ok){
      tn = TN;
      tv = TV;
      
      sv = I->sv;
      sn = I->sn;
      for(b = 0; b <= I->Ns; b++) {
	if(b == I->Ns) {
	  sv = I->sv;
	  sn = I->sn;
	}
	v = I->p;
	n = I->n;
	
	for(a = 0; a < I->N; a++) {
	  transform33Tf3f(n, sv, tv);
	  add3f(v, tv, tv);
	  tv += 3;
	  transform33Tf3f(n, sn, tn);
	  tn += 3;
	  n += 9;
	  v += 3;
	}
	sv += 3;
	sn += 3;
      }
      
      /* fill in each strip separately */
      
      tv = TV;
      tn = TN;
      
      tv1 = TV + 3 * I->N;
      tn1 = TN + 3 * I->N;
    }
    for(b = 0; ok && b < I->Ns; b += 2) {
      if(SettingGetGlobal_i(I->G, cSetting_cartoon_debug) < 1.5)
        ok &= CGOBegin(cgo, GL_TRIANGLE_STRIP);
      else {
        ok &= CGOBegin(cgo, GL_LINE_STRIP);
      }
      if(ok && color_override)
        ok &= CGOColorv(cgo, color_override);
      c = I->c;
      alpha = I->alpha;
      i = I->i;
      for(a = 0; ok && a < I->N; a++) {
        if(!color_override)
          ok &= CGOColorv(cgo, c);
        if (ok){
          ok &= CGOAlpha(cgo, *alpha);
        }
        if (ok)
	  ok &= CGOPickColor(cgo, *i, cPickableAtom);
	if (ok)
	  ok &= CGONormalv(cgo, tn);
	if (ok)
	  ok &= CGOVertexv(cgo, tv);
        tn += 3;
        tv += 3;
	if (ok)
	  ok &= CGONormalv(cgo, tn1);
	if (ok)
	  ok &= CGOVertexv(cgo, tv1);
        tn1 += 3;
        tv1 += 3;
        c += 3;
        alpha++;
        i++;
      }
      tv += 3 * I->N;
      tn += 3 * I->N;
      tv1 += 3 * I->N;
      tn1 += 3 * I->N;
      if (ok)
	ok &= CGOEnd(cgo);
      if (ok)
	ok &= CGOPickColor(cgo, -1, cPickableNoPick);
    }

    if(ok && cap != cCylCap::None) {

      if(color_override)
        ok &= CGOColorv(cgo, color_override);

      if (ok){
	n = I->n;
	v = I->p;
	
	sv = I->sv;
	tv = I->tv;
	for(b = 0; b < I->Ns; b++) {
	  transform33Tf3f(n, sv, tv);
	  add3f(v, tv, tv);
	  sv += 3;
	  tv += 3;
	}
      }
      if (ok)
	ok &= CGOBegin(cgo, GL_TRIANGLE_FAN);
      if (ok){
	copy3f(I->n, v0);
	invert3f(v0);
	if(!color_override)
	  ok &= CGOColorv(cgo, I->c);
    if (ok){
      ok &= CGOAlpha(cgo, I->alpha[0]);
    }
	if (ok)
	  ok &= CGOPickColor(cgo, I->i[0], cPickableAtom);
	if (ok)
	  ok &= CGONormalv(cgo, v0);
      }
      if (ok)
	ok &= CGOVertexv(cgo, v);
      /* trace shape */
      if (ok)
	ok &= CGOVertexv(cgo, I->tv);
      for(b = I->Ns - 1; ok && b >= 0; b--) {
        ok &= CGOVertexv(cgo, I->tv + b * 3);
      }
      if (ok)
	ok &= CGOEnd(cgo);
      if (ok)
	ok &= CGOPickColor(cgo, -1, cPickableNoPick);
      if (ok){
	n = I->n + 9 * (I->N - 1);
	v = I->p + 3 * (I->N - 1);
	
	sv = I->sv;
	tv = I->tv;
	for(b = 0; b < I->Ns; b++) {
	  transform33Tf3f(n, sv, tv);
	  add3f(v, tv, tv);
	  sv += 3;
	  tv += 3;
	}
      }
      if (ok)
	ok &= CGOBegin(cgo, GL_TRIANGLE_FAN);
      if(ok && !color_override)
        ok &= CGOColorv(cgo, I->c + 3 * (I->N - 1));
      if (ok){
        ok &= CGOAlpha(cgo, I->alpha[I->N - 1]);
      }
      if (ok)
	ok &= CGOPickColor(cgo, I->i[I->N - 1], cPickableAtom);
      if (ok)
	ok &= CGONormalv(cgo, n);
      if (ok)
	ok &= CGOVertexv(cgo, v);
      /* trace shape */
      for(b = 0; ok && b < I->Ns; b++) {
        ok &= CGOVertexv(cgo, I->tv + b * 3);
      }
      if (ok)
	ok &= CGOVertexv(cgo, I->tv);
      if (ok)
	ok &= CGOEnd(cgo);
      if (ok)
	ok &= CGOPickColor(cgo, -1, cPickableNoPick);
    }
    FreeP(TV);
    FreeP(TN);
  }

  PRINTFD(I->G, FB_Extrude)
    " ExtrudeCGOSurfacePolygon-DEBUG: exiting...\n" ENDFD;
  return ok;
}

int ExtrudeCGOSurfacePolygonTaper(const CExtrude * I, CGO * cgo, int sampling,
				  const float *color_override)
{
  int a, b;
  unsigned int *i;
  float *v;
  float *n;
  float *c;
  const float *alpha;
  float *sv, *sn, *tv, *tn, *tv1, *tn1, *TV = NULL, *TN = NULL;
  float s0[3];
  float f;
  int subN;
  int ok = true;

  subN = I->N - sampling;

  PRINTFD(I->G, FB_Extrude)
    " ExtrudeCGOSurfacePolygonTaper-DEBUG: entered.\n" ENDFD;

  if(I->N && I->Ns) {

    TV = pymol::malloc<float>(3 * (I->Ns + 1) * I->N);
    CHECKOK(ok, TV);
    if (ok)
      TN = pymol::malloc<float>(3 * (I->Ns + 1) * I->N);
    CHECKOK(ok, TN);
    /* compute transformed shape vertices */

    if (ok){
    tn = TN;
    tv = TV;

    sv = I->sv;
    sn = I->sn;
    for(b = 0; b <= I->Ns; b++) {
      if(b == I->Ns) {
        sv = I->sv;
        sn = I->sn;
      }
      v = I->p;
      n = I->n;

      for(a = 0; a < I->N; a++) {
        if((a >= sampling) && (a < subN)) {

          transform33Tf3f(n, sv, tv);
          add3f(v, tv, tv);
          tv += 3;
          transform33Tf3f(n, sn, tn);
          tn += 3;
          n += 9;
          v += 3;
        } else {
          copy3f(sv, s0);

          if(a >= subN) {
            f = ((I->N - a - 1) / ((float) sampling));
          } else if(a < sampling) {
            f = (a / ((float) sampling));
          } else
            f = 1.0;
          f = smooth(f, 2);
          s0[2] *= f;

          transform33Tf3f(n, s0, tv);
          add3f(v, tv, tv);
          tv += 3;
          transform33Tf3f(n, sn, tn);
          tn += 3;
          n += 9;
          v += 3;

        }
      }
      sv += 3;
      sn += 3;
    }

    /* fill in each strip separately */

    tv = TV;
    tn = TN;

    tv1 = TV + 3 * I->N;
    tn1 = TN + 3 * I->N;
    }
    for(b = 0; ok && b < I->Ns; b += 2) {
      if (ok){
	if(SettingGetGlobal_i(I->G, cSetting_cartoon_debug) < 1.5)
	  ok &= CGOBegin(cgo, GL_TRIANGLE_STRIP);
	else {
	  ok &= CGOBegin(cgo, GL_LINE_STRIP);
	}
      }
      if(ok && color_override)
        ok &= CGOColorv(cgo, color_override);
      c = I->c;
      alpha = I->alpha;
      i = I->i;
      for(a = 0; ok && a < I->N; a++) {
        if(!color_override)
          ok &= CGOColorv(cgo, c);
        if (ok){
          ok &= CGOAlpha(cgo, *alpha);
        }
        if (ok)
	  ok &= CGOPickColor(cgo, *i, cPickableAtom);
        if (ok)
	  ok &= CGONormalv(cgo, tn);
	if (ok)
	  ok &= CGOVertexv(cgo, tv);
        tn += 3;
        tv += 3;
	if (ok)
	  ok &= CGONormalv(cgo, tn1);
	if (ok)
	  ok &= CGOVertexv(cgo, tv1);
        tn1 += 3;
        tv1 += 3;
        c += 3;
        alpha++;
        i++;
      }
      if (ok){
	tv += 3 * I->N;
	tn += 3 * I->N;
	tv1 += 3 * I->N;
	tn1 += 3 * I->N;
	CGOEnd(cgo);
	CGOPickColor(cgo, -1, cPickableNoPick);
      }
    }

    FreeP(TV);
    FreeP(TN);
  }

  PRINTFD(I->G, FB_Extrude)
    " ExtrudeCGOSurfacePolygonTaper-DEBUG: exiting...\n" ENDFD;
  return ok;
}

int ExtrudeCGOSurfaceStrand(const CExtrude * I, CGO * cgo, int sampling, const float *color_override)
{
  int a, b;
  unsigned int *i;
  float *v;
  float *n;
  float *c;
  const float *alpha;
  float *sv, *sn, *tv, *tn, *tv1, *tn1, *TV = NULL, *TN = NULL;
  float v0[3], n0[3], s0[3], z[3] = { 1.0, 0.0, 1.0 };
  int subN;
  int ok = true;
  subN = I->N - sampling;

  PRINTFD(I->G, FB_Extrude)
    " ExtrudeCGOSurfaceStrand-DEBUG: entered.\n" ENDFD;

  if(I->N && I->Ns) {

    TV = pymol::malloc<float>(3 * (I->Ns + 1) * I->N);
    CHECKOK(ok, TV);
    if (ok)
      TN = pymol::malloc<float>(3 * (I->Ns + 1) * I->N);
    CHECKOK(ok, TN);
    /* compute transformed shape vertices */

    if (ok){
      tn = TN;
      tv = TV;
      sv = I->sv;
      sn = I->sn;
      for(b = 0; b <= I->Ns; b++) {
	if(b == I->Ns) {
	  sv = I->sv;
	  sn = I->sn;
	}
	v = I->p;
	n = I->n;
	for(a = 0; a < I->N; a++) {
	  copy3f(sv, s0);
	  if(a == subN) {
	    scale3f(s0, 0.50F, s0);
	  }
	  transform33Tf3f(n, s0, tv);
	  add3f(v, tv, tv);
	  tv += 3;
	  transform33Tf3f(n, sn, tn);
	  tn += 3;
	  n += 9;
	  v += 3;
	}
	sv += 3;
	sn += 3;
      }

      /* fill in each strip of arrow separately */
      
      tv = TV;
      tn = TN;
      
      tv1 = TV + 3 * I->N;
      tn1 = TN + 3 * I->N;
    }

    for(b = 0; ok && b < I->Ns; b += 2) {
      if(SettingGetGlobal_i(I->G, cSetting_cartoon_debug) < 1.5)
        ok &= CGOBegin(cgo, GL_TRIANGLE_STRIP);
      else {
        ok &= CGOBegin(cgo, GL_LINE_STRIP);
      }
      c = I->c;
      alpha = I->alpha;
      i = I->i;
      for(a = 0; ok && a < I->N; a++) {
        if(a < subN) {
	  if(color_override && ((b == 2) || (b == 3) || (b == 6) || (b == 7)))
	    ok &= CGOColorv(cgo, color_override);
	  else
	    ok &= CGOColorv(cgo, c);
      if (ok){
        ok &= CGOAlpha(cgo, *alpha);
      }
	  if (ok)
	    ok &= CGOPickColor(cgo, *i, cPickableAtom);
	  if (ok)
	    ok &= CGONormalv(cgo, tn);
	  if (ok)
	    ok &= CGOVertexv(cgo, tv);
        }
        tn += 3;
        tv += 3;
        if(ok && a < subN) {
          ok &= CGONormalv(cgo, tn1);
	  if (ok)
	    ok &= CGOVertexv(cgo, tv1);
        }
        tn1 += 3;
        tv1 += 3;
        c += 3;
        alpha++;
        i++;
      }
      tv += 3 * I->N;
      tn += 3 * I->N;
      tv1 += 3 * I->N;
      tn1 += 3 * I->N;
      if (ok)
	ok &= CGOEnd(cgo);
      if (ok)
	ok &= CGOPickColor(cgo, -1, cPickableNoPick);
    }

    if(ok) {

      n = I->n;
      v = I->p;

      sv = I->sv;
      tv = I->tv;
      for(b = 0; b < I->Ns; b++) {
        transform33Tf3f(n, sv, tv);
        add3f(v, tv, tv);
        sv += 3;
        tv += 3;
      }

      copy3f(I->n, v0);
      invert3f(v0);
      if (ok){
	if(color_override)
	  ok &= CGOColorv(cgo, color_override);
	else
	  ok &= CGOColorv(cgo, I->c);
       }
      if(ok){
        ok &= CGOAlpha(cgo, I->alpha[0]);
      }
      if (ok)
	ok &= CGOPickColor(cgo, I->i[0], cPickableAtom);

      if (ok)
	ok &= CGOBegin(cgo, GL_TRIANGLE_STRIP);
      if (ok)
	ok &= CGONormalv(cgo, v0);

      ok &= CGOVertexv(cgo, I->tv + 6 * 3);
      ok &= CGOVertexv(cgo, I->tv + 4 * 3);
      ok &= CGOVertexv(cgo, I->tv + 0 * 3);
      ok &= CGOVertexv(cgo, I->tv + 2 * 3);
      if (ok)
	ok &= CGOEnd(cgo);
      if (ok)
	ok &= CGOPickColor(cgo, -1, cPickableNoPick);
    }

    /* now do the arrow part */

    tn = TN;
    tv = TV;

    sv = I->sv;
    sn = I->sn;
    if (ok){
      for(b = 0; b <= I->Ns; b++) {
	if(b == I->Ns) {
	  sv = I->sv;
	  sn = I->sn;
	}
	v = I->p;
	n = I->n;
	
	for(a = 0; a < I->N; a++) {
	  copy3f(sv, s0);
	  s0[2] = s0[2] * ((1.5F * ((I->N - 1) - a)) / sampling);
	  transform33Tf3f(n, s0, tv);
	  add3f(v, tv, tv);
	  tv += 3;
	  copy3f(sn, n0);
	  if(fabs(dot_product3f(sn, z)) > R_SMALL4) {
	    n0[0] += 0.4F;
	    normalize3f(n0);
	  }
	  transform33Tf3f(n, n0, tn);
	  tn += 3;
	  n += 9;
	  v += 3;
	}
	sv += 3;
	sn += 3;
      }
    }

    tv = TV;
    tn = TN;

    tv1 = TV + 3 * I->N;
    tn1 = TN + 3 * I->N;

    for(b = 0; ok && b < I->Ns; b += 2) {
      if (ok){
	if(SettingGetGlobal_i(I->G, cSetting_cartoon_debug) < 1.5)
	  ok &= CGOBegin(cgo, GL_TRIANGLE_STRIP);
	else {
	  ok &= CGOBegin(cgo, GL_LINE_STRIP);
	}
      }
      c = I->c;
      alpha = I->alpha;
      i = I->i;
      for(a = 0; ok && a < I->N; a++) {
        if(a >= (subN - 1)) {
          if(color_override && ((b == 2) || (b == 3) || (b == 6) || (b == 7)))
            ok &= CGOColorv(cgo, color_override);
          else
            ok &= CGOColorv(cgo, c);
        if (ok){
          ok &= CGOAlpha(cgo, *alpha);
        }
	  if (ok)
	    ok &= CGOPickColor(cgo, *i, cPickableAtom);
          if (ok)
	    ok &= CGONormalv(cgo, tn);
	  if (ok)
	    ok &= CGOVertexv(cgo, tv);
        }
        tn += 3;
        tv += 3;
        if(ok && a >= (subN - 1)) {
          ok &= CGONormalv(cgo, tn1);
	  if (ok)
	    ok &= CGOVertexv(cgo, tv1);
        }
        tn1 += 3;
        tv1 += 3;
        c += 3;
        alpha++;
        i++;
      }
      tv += 3 * I->N;
      tn += 3 * I->N;
      tv1 += 3 * I->N;
      tn1 += 3 * I->N;
      if (ok)
	ok &= CGOEnd(cgo);
      if (ok)
	ok &= CGOPickColor(cgo, -1, cPickableNoPick);
    }

    n = I->n + 9 * (subN - 1);
    v = I->p + 3 * (subN - 1);
    sv = I->sv;
    tv = I->tv;
    if (ok){
      for(b = 0; b < I->Ns; b++) {
	copy3f(sv, s0);
	s0[2] = s0[2] * ((s0[2] < 0.f) ? -1.f : 1.5F) ;
	transform33Tf3f(n, s0, tv);
	add3f(v, tv, tv);
	sv += 3;
	tv += 3;
      }
    }

    // Back end/flat surface of the arrow, 
    // now split into two pieces, one for each side
    copy3f(n, v0);
    invert3f(v0);
    if (ok){
      if(color_override)
	ok &= CGOColorv(cgo, color_override);
      else
	ok &= CGOColorv(cgo, I->c + 3 * (subN - 1));
    }
    if (ok){
      ok &= CGOAlpha(cgo, I->alpha[(subN - 1)]);
    }
    if (ok)
      ok &= CGOPickColor(cgo, I->i[(subN - 1)], cPickableAtom);

    // draw first side
    tv = I->tv;
    if (ok)
      ok &= CGOBegin(cgo, GL_TRIANGLE_STRIP);
    if (ok)
      ok &= CGONormalv(cgo, v0);
    ok &= CGOVertexv(cgo, I->tv + 6 * 3);
    ok &= CGOVertexv(cgo, I->tv + 4 * 3);
    ok &= CGOVertexv(cgo, I->tv + 0 * 3);
    ok &= CGOVertexv(cgo, I->tv + 2 * 3);
    if (ok)
      ok &= CGOEnd(cgo);

    // switch vertices to other side of flat surface of arrow
    sv = I->sv;
    tv = I->tv;
    if (ok){
      for(b = 0; b < I->Ns; b++) {
	copy3f(sv, s0);
	s0[2] = s0[2] * ((s0[2] < 0.f) ? 1.f : -1.5F) ;
	transform33Tf3f(n, s0, tv);
	add3f(v, tv, tv);
	sv += 3;
	tv += 3;
      }
    }
    // draw other side
    if (ok)
      ok &= CGOBegin(cgo, GL_TRIANGLE_STRIP);
    if (ok)
      ok &= CGONormalv(cgo, v0);
    ok &= CGOVertexv(cgo, I->tv + 0 * 3);
    ok &= CGOVertexv(cgo, I->tv + 2 * 3);
    ok &= CGOVertexv(cgo, I->tv + 6 * 3);
    ok &= CGOVertexv(cgo, I->tv + 4 * 3);
    if (ok)
      ok &= CGOEnd(cgo);

    if (ok)
      ok &= CGOPickColor(cgo, -1, cPickableNoPick);

    FreeP(TV);
    FreeP(TN);
  }

  PRINTFD(I->G, FB_Extrude)
    " ExtrudeCGOSurfaceStrand-DEBUG: exiting...\n" ENDFD;
  return ok;
}

int ExtrudeComputePuttyScaleFactors(CExtrude * I, ObjectMolecule * obj, int transform,
				    float mean, float stdev, float min, float max,
				    float power, float range,
				    float min_scale, float max_scale, int window)
{
  float *sf;
  int a;
  unsigned int *i;
  AtomInfoType *at;
  float scale = 1.0F;
  float data_range = max - min;
  int ok = true;

  if(I->N && I->Ns) {
    int invalid = false;
    i = I->i;
    sf = I->sf;

    /* guard against invalid inputs that would imply division by zero */

    switch (transform) {
    case cPuttyTransformNormalizedNonlinear:
    case cPuttyTransformNormalizedLinear:
      /* depend on stdev */
      if(stdev < R_SMALL8)
        invalid = true;
      break;
    }
    switch (transform) {
    case cPuttyTransformNormalizedNonlinear:
    case cPuttyTransformRelativeNonlinear:
    case cPuttyTransformScaledNonlinear:
    case cPuttyTransformNormalizedLinear:
    case cPuttyTransformRelativeLinear:
    case cPuttyTransformScaledLinear:
      /* depend on range */

      if(fabs(range) < R_SMALL8)
        invalid = true;
      break;
    }
    switch (transform) {
    case cPuttyTransformRelativeNonlinear:
    case cPuttyTransformRelativeLinear:
      /* depend on data_range */
      if(fabs(data_range) < R_SMALL8)
        invalid = true;
      break;
    }

    if(!invalid) {
      for(a = 0; a < I->N; a++) {
        at = obj->AtomInfo + (*i);
        switch (transform) {
        case cPuttyTransformNormalizedNonlinear:
          /* normalized by Z-score, with the range affecting the distribution width */
          scale = (range + (at->b - mean) / stdev) / range;
          if(scale < 0.0F)
            scale = 0.0F;
          scale = (float) pow(scale, power);
          break;
        case cPuttyTransformRelativeNonlinear:
          scale = (at->b - min) / (data_range * range);
          if(scale < 0.0F)
            scale = 0.0F;
          scale = (float) pow(scale, power);
          *sf = scale;
          break;
        case cPuttyTransformScaledNonlinear:
          scale = at->b / range;
          if(scale < 0.0F)
            scale = 0.0F;
          scale = (float) pow(scale, power);
          *sf = scale;
          break;
        case cPuttyTransformAbsoluteNonlinear:
          scale = at->b;
          if(scale < 0.0F)
            scale = 0.0F;
          scale = (float) pow(scale, power);
          *sf = scale;
          break;
        case cPuttyTransformNormalizedLinear:
          /* normalized by Z-score, with the range affecting the distribution width */
          scale = (range + (at->b - mean) / stdev) / range;
          if(scale < 0.0F)
            scale = 0.0F;
          break;
        case cPuttyTransformRelativeLinear:
          scale = (at->b - min) / (data_range * range);
          if(scale < 0.0F)
            scale = 0.0F;
          *sf = scale;
          break;
        case cPuttyTransformScaledLinear:
          scale = at->b / range;
          if(scale < 0.0F)
            scale = 0.0F;
          *sf = scale;
          break;
        case cPuttyTransformAbsoluteLinear:
          scale = at->b;
          if(scale < 0.0F)
            scale = 0.0F;
          *sf = scale;
          break;
        case cPuttyTransformImpliedRMS:
          if(scale < 0.0F)
            scale = 0.0F;
          scale = (float) (sqrt1d(at->b / 8.0) / PI);
          break;
        }
        if((scale < min_scale) && (min_scale >= 0.0))
          scale = min_scale;
        if((scale > max_scale) && (max_scale >= 0.0))
          scale = max_scale;
        *(sf++) = scale;
        i++;
      }
    } else {
      PRINTFB(I->G, FB_RepCartoon, FB_Warnings)
        " Extrude-Warning: invalid putty settings (division by zero)\n" ENDFB(I->G);
      for(a = 0; a < I->N; a++) {
        *sf = 0.5F;
        sf++;
      }
    }

    PRINTFB(I->G, FB_RepCartoon, FB_Blather)
      " Putty: mean %8.3f stdev %8.3f min %8.3f max %8.3f\n",
      mean, stdev,
      mean + (pow(min_scale, 1.0F / power) * range - range) * stdev,
      mean + (pow(max_scale, 1.0F / power) * range - range) * stdev ENDFB(I->G);
    /* now compute window average */

    {
      float *SF = pymol::malloc<float>(I->N);
      int w, ww;
      float accum;
      int cnt;

      CHECKOK(ok, SF);

      sf = I->sf;

      if (ok){
	for(a = 1; a < (I->N - 1); a++) {
	  accum = 0.0F;
	  cnt = 0;
	  for(w = -window; w <= window; w++) {
	    ww = w + a;
	    if(ww < 0)
	      ww = 0;
	    else if(ww > (I->N - 1))
	      ww = I->N - 1;
	    accum += sf[ww];
	    cnt++;
	  }
	  SF[a] = accum / cnt;
	}
	for(a = 1; a < I->N - 1; a++)
	  sf[a] = SF[a];
	FreeP(SF);
      }
    }
  }
  return (ok);
}

#if 0
#endif

void ExtrudeTruncate(CExtrude * I, int n)
{

  I->N = n;
  /* should free RAM here... */
}

int ExtrudeAllocPointsNormalsColors(CExtrude * I, int n)
{
  int ok = true;
  if(I->N < n) {
    /* reset */
    FreeP(I->p);
    FreeP(I->n);
    FreeP(I->c);
    FreeP(I->alpha);
    FreeP(I->i);
    FreeP(I->sf);               /* PUTTY */
    I->p = pymol::malloc<float>(3 * (n + 1));
    CHECKOK(ok, I->p);
    if (ok)
      I->n = pymol::malloc<float>(9 * (n + 1));
    CHECKOK(ok, I->n);
    if (ok)
      I->c = pymol::malloc<float>(3 * (n + 1));
    CHECKOK(ok, I->c);
    if (ok)
      I->alpha = pymol::malloc<float>(n + 1);
    CHECKOK(ok, I->alpha);
    if (ok)
      I->i = pymol::malloc<unsigned int>(3 * (n + 1));
    CHECKOK(ok, I->i);
    if (ok)
      I->sf = pymol::malloc<float>(n + 1);        /* PUTTY: scale factors */
    CHECKOK(ok, I->sf);
    if (!ok){
      FreeP(I->p);
      FreeP(I->n);
      FreeP(I->c);
      FreeP(I->alpha);
      FreeP(I->i);
      FreeP(I->sf);
    }
  }
  I->N = n;
  return ok;
}

void ExtrudeFree(CExtrude * I)
{
  FreeP(I->p);
  FreeP(I->n);
  FreeP(I->c);
  FreeP(I->alpha);
  FreeP(I->tn);
  FreeP(I->tv);
  FreeP(I->sn);
  FreeP(I->sv);
  FreeP(I->i);
  FreeP(I->sf);
  OOFreeP(I);
}
