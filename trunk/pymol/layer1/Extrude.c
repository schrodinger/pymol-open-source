
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
#include"os_python.h"

#include"os_predef.h"
#include"os_std.h"
#include"os_gl.h"

#include"Extrude.h"
#include"Base.h"
#include"OOMac.h"
#include"Setting.h"
#include"Feedback.h"

void ExtrudeInit(PyMOLGlobals * G, CExtrude * I);

static float smooth(float x, float power)
{

  if(x <= 0.5) {
    if(x <= 0.0)
      x = 0.0;
    return ((float) (0.5F * pow(2.0 * x, power)));
  } else {
    if(x >= 1.0)
      x = 1.0;
    return ((float) (1.0F - (0.5 * pow(2 * (1.0 - x), power))));
  }
}

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
    CopyArray(I->i, orig->i, int, I->N);
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

  I->sv = Alloc(float, 3 * (n + 1));
  CHECKOK(ok, I->sv);
  if (ok)
    I->sn = Alloc(float, 3 * (n + 1));
  CHECKOK(ok, I->sn);
  if (ok)
    I->tv = Alloc(float, 3 * (n + 1));
  CHECKOK(ok, I->tv);
  if (ok)
    I->tn = Alloc(float, 3 * (n + 1));
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

  I->sv = Alloc(float, 3 * (n + 1));
  CHECKOK(ok, I->sv);
  if (ok)
    I->sn = Alloc(float, 3 * (n + 1));
  CHECKOK(ok, I->sn);
  if (ok)
    I->tv = Alloc(float, 3 * (n + 1));
  CHECKOK(ok, I->tv);
  if (ok)
    I->tn = Alloc(float, 3 * (n + 1));
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

  I->sv = Alloc(float, 3 * (I->Ns + 1));
  CHECKOK(ok, I->sv);
  if (ok)
    I->sn = Alloc(float, 3 * (I->Ns + 1));
  CHECKOK(ok, I->sn);
  if (ok)
    I->tv = Alloc(float, 3 * (I->Ns + 1));
  CHECKOK(ok, I->tv);
  if (ok)
    I->tn = Alloc(float, 3 * (I->Ns + 1));
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

  I->sv = Alloc(float, 3 * (I->Ns + 1));
  CHECKOK(ok, I->sv);
  if (ok)
    I->sn = Alloc(float, 3 * (I->Ns + 1));
  CHECKOK(ok, I->sn);
  if (ok)
    I->tv = Alloc(float, 3 * (I->Ns + 1));
  CHECKOK(ok, I->tv);
  if (ok)
    I->tn = Alloc(float, 3 * (I->Ns + 1));
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
  
  I->sv = Alloc(float, 3 * (n + 1));
  CHECKOK(ok, I->sv);
  if (ok)
    I->sn = Alloc(float, 3 * (n + 1));
  CHECKOK(ok, I->sn);
  if (ok)
    I->tv = Alloc(float, 3 * (n + 1));
  CHECKOK(ok, I->tv);
  if (ok)
    I->tn = Alloc(float, 3 * (n + 1));
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

void ExtrudeCGOTraceAxes(CExtrude * I, CGO * cgo)
{
  int a;
  float *v, *n;
  float v0[3];

  if(I->N) {
    CGOColor(cgo, 0.5, 0.5, 0.5);
#ifdef _PYMOL_CGO_DRAWARRAYS
    {
      int nverts = 6*I->N, pl = 0;
      float *vertexVals, *tmp_ptr;
      vertexVals = CGODrawArrays(cgo, GL_LINES, CGO_VERTEX_ARRAY, nverts);
      v = I->p;
      n = I->n;
      for(a = 0; a < I->N; a++) {
	add3f(v, n, v0);
	tmp_ptr = v0;
	vertexVals[pl++] = tmp_ptr[0]; vertexVals[pl++] = tmp_ptr[1]; vertexVals[pl++] = tmp_ptr[2];
	tmp_ptr = v;
	vertexVals[pl++] = tmp_ptr[0]; vertexVals[pl++] = tmp_ptr[1]; vertexVals[pl++] = tmp_ptr[2];
	n += 3;
	add3f(v, n, v0);
	tmp_ptr = v0;
	vertexVals[pl++] = tmp_ptr[0]; vertexVals[pl++] = tmp_ptr[1]; vertexVals[pl++] = tmp_ptr[2];
	tmp_ptr = v;
	vertexVals[pl++] = tmp_ptr[0]; vertexVals[pl++] = tmp_ptr[1]; vertexVals[pl++] = tmp_ptr[2];
	n += 3;
	add3f(v, n, v0);
	tmp_ptr = v0;
	vertexVals[pl++] = tmp_ptr[0]; vertexVals[pl++] = tmp_ptr[1]; vertexVals[pl++] = tmp_ptr[2];
	tmp_ptr = v;
	vertexVals[pl++] = tmp_ptr[0]; vertexVals[pl++] = tmp_ptr[1]; vertexVals[pl++] = tmp_ptr[2];
	n += 3;
	v += 3;
      }
    }
#else
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
#endif
  }
}

void ExtrudeCGOTrace(CExtrude * I, CGO * cgo)
{
  int a;
  float *v;
  if(I->N) {
    CGOColor(cgo, 0.5, 0.5, 0.5);
#ifdef _PYMOL_CGO_DRAWARRAYS
    {
      int nverts = I->N, pl = 0;
      float *vertexVals;
      vertexVals = CGODrawArrays(cgo, GL_LINE_STRIP, CGO_VERTEX_ARRAY, nverts);
      v = I->p;
      for(a = 0; a < I->N; a++) {
	vertexVals[pl] = v[0]; vertexVals[pl+1] = v[1]; vertexVals[pl+2] = v[2];
	v += 3;
      }
    }
#else
    CGOBegin(cgo, GL_LINE_STRIP);
    v = I->p;
    for(a = 0; a < I->N; a++) {
      CGOVertexv(cgo, v);
      v += 3;
    }
    CGOEnd(cgo);
#endif
  }
}

int ExtrudeComputeTangents(CExtrude * I)
{
  float *nv, *v1, *v;
  int a;
  int ok = true;

  PRINTFD(I->G, FB_Extrude)
    " ExtrudeComputeTangents-DEBUG: entered.\n" ENDFD;

  nv = Alloc(float, I->N * 3);
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
#ifdef _PYMOL_CGO_DRAWARRAYS
      int nverts = 2*I->N + 2, pl = 0;
      float *vertexVals, *tmp_ptr;
      vertexVals = CGODrawArrays(cgo, GL_LINES, CGO_VERTEX_ARRAY, nverts);
#else
      CGOBegin(cgo, GL_LINES);
#endif
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
#ifdef _PYMOL_CGO_DRAWARRAYS
	  tmp_ptr = v0;
	  vertexVals[pl++] = tmp_ptr[0]; vertexVals[pl++] = tmp_ptr[1]; vertexVals[pl++] = tmp_ptr[2];
	  tmp_ptr = v1;
	  vertexVals[pl++] = tmp_ptr[0]; vertexVals[pl++] = tmp_ptr[1]; vertexVals[pl++] = tmp_ptr[2];
#else
	  CGOVertexv(cgo, v0);
	  CGOVertexv(cgo, v1);
#endif
	  copy3f(v1, v0);
	}
	tv = I->tv;
	add3f(v, tv, v1);
#ifdef _PYMOL_CGO_DRAWARRAYS
	tmp_ptr = v0;
	vertexVals[pl++] = tmp_ptr[0]; vertexVals[pl++] = tmp_ptr[1]; vertexVals[pl++] = tmp_ptr[2];
	tmp_ptr = v1;
	vertexVals[pl++] = tmp_ptr[0]; vertexVals[pl++] = tmp_ptr[1]; vertexVals[pl++] = tmp_ptr[2];
#else
	CGOVertexv(cgo, v0);
	CGOVertexv(cgo, v1);
#endif
	v += 3;
	n += 9;
      }
#ifndef _PYMOL_CGO_DRAWARRAYS
      CGOEnd(cgo);
#endif
    }
  }
}

int ExtrudeCGOSurfaceTube(CExtrude * I, CGO * cgo, int cap, float *color_override, short use_spheres)
{
  int a, b, *i;
  float *v;
  float *n;
  float *c;
  float *sv, *sn, *tv, *tn, *tv1, *tn1, *TV, *TN;
  float v0[3];
  int start, stop;
  int ok = true;
  PRINTFD(I->G, FB_Extrude)
    " ExtrudeCGOSurfaceTube-DEBUG: entered.\n" ENDFD;

  if(I->N && I->Ns) {
    TV = Alloc(float, 3 * (I->Ns + 1) * I->N);
    CHECKOK(ok, TV);
    if (ok)
      TN = Alloc(float, 3 * (I->Ns + 1) * I->N);
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
      
      start = I->Ns / 4;
      stop = 3 * I->Ns / 4;
    }
    for(b = 0; ok && b < I->Ns; b++) {
#ifdef _PYMOL_CGO_DRAWARRAYS
      GLenum mode;
      if(SettingGetGlobal_i(I->G, cSetting_cartoon_debug) < 1.5)
	mode = GL_TRIANGLE_STRIP;
      else {
	  mode = GL_LINE_STRIP;
      }
      c = I->c;
      i = I->i;
      if (ok) {
	int nverts = 2*I->N, pl = 0, plc = 0, damode = CGO_VERTEX_ARRAY | CGO_NORMAL_ARRAY | CGO_PICK_COLOR_ARRAY, nxtn = 3;
	float *vertexVals = NULL, *normalVals, *colorVals = 0, *nxtVals = 0, *tmp_ptr;
	float *pickColorVals;
	if (color_override && (b > start) && (b < stop))
	  CGOColorv(cgo, color_override);
	else {
	  damode |= CGO_COLOR_ARRAY;
	}
	vertexVals = CGODrawArrays(cgo, mode, damode, nverts);
	CHECKOK(ok, vertexVals);
	if (ok){
	  normalVals = nxtVals = vertexVals + (nverts*3);
	  if (damode & CGO_COLOR_ARRAY){
	    colorVals = nxtVals = normalVals + (nverts*3);
	    nxtn = 4;
	  }
	  pickColorVals = nxtVals + (nverts*nxtn);
	  
	  for(a = 0; a < I->N; a++) {
	    if(colorVals){
	      tmp_ptr = c;
	      colorVals[plc++] = tmp_ptr[0]; colorVals[plc++] = tmp_ptr[1]; colorVals[plc++] = tmp_ptr[2]; colorVals[plc++] = cgo->alpha;
	    }
	    SetCGOPickColor(pickColorVals, nverts, pl, *i, cPickableAtom);	  
	    tmp_ptr = tn;
	    normalVals[pl] = tmp_ptr[0]; normalVals[pl+1] = tmp_ptr[1]; normalVals[pl+2] = tmp_ptr[2];
	    tmp_ptr = tv;
	    vertexVals[pl] = tmp_ptr[0]; vertexVals[pl+1] = tmp_ptr[1]; vertexVals[pl+2] = tmp_ptr[2];
	    pl += 3;
	    tn += 3;
	    tv += 3;
	    if(colorVals){
	      tmp_ptr = c;
	      colorVals[plc++] = tmp_ptr[0]; colorVals[plc++] = tmp_ptr[1]; colorVals[plc++] = tmp_ptr[2]; colorVals[plc++] = cgo->alpha;
	    }
	    SetCGOPickColor(pickColorVals, nverts, pl, *i, cPickableAtom);	  
	    tmp_ptr = tn1;
	    normalVals[pl] = tmp_ptr[0]; normalVals[pl+1] = tmp_ptr[1]; normalVals[pl+2] = tmp_ptr[2];
	    tmp_ptr = tv1;
	    vertexVals[pl] = tmp_ptr[0]; vertexVals[pl+1] = tmp_ptr[1]; vertexVals[pl+2] = tmp_ptr[2];
	    pl += 3;
	    tn1 += 3;
	    tv1 += 3;
	    c += 3;
	    i++;
	  }
	}
      }
#else
      if (ok){
	if(SettingGetGlobal_i(I->G, cSetting_cartoon_debug) < 1.5)
	  ok &= CGOBegin(cgo, GL_TRIANGLE_STRIP);
	else {
	  ok &= CGOBegin(cgo, GL_LINE_STRIP);
	}
	if (ok){
	  c = I->c;
	  i = I->i;
	  for(a = 0; ok && a < I->N; a++) {
	    if(color_override && (b > start) && (b < stop))
	      ok &= CGOColorv(cgo, color_override);
	    else
	      ok &= CGOColorv(cgo, c);
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
	    i++;
	  }
	}
	if (ok)
	  ok &= CGOEnd(cgo);
      }
#endif
    }

    if (ok){
    switch (cap) {
    case 1:
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

#ifdef _PYMOL_CGO_DRAWARRAYS
      if (ok) {
	int nverts = 2 + I->Ns, pl = 0;
	float *vertexVals = NULL, *tmp_ptr;
	copy3f(I->n, v0);
	invert3f(v0);
	if(color_override)
	  ok &= CGOColorv(cgo, color_override);
	else
	  ok &= CGOColorv(cgo, I->c);
	if (ok)
	  ok &= CGOPickColor(cgo, I->i[0], cPickableAtom);
	if (ok)
	  ok &= CGONormalv(cgo, v0);
	if (ok)
	  vertexVals = CGODrawArrays(cgo, GL_TRIANGLE_FAN, CGO_VERTEX_ARRAY, nverts);	      
	CHECKOK(ok, vertexVals);
	if (ok){
	  tmp_ptr = v;
	  vertexVals[pl] = tmp_ptr[0]; vertexVals[pl+1] = tmp_ptr[1]; vertexVals[pl+2] = tmp_ptr[2];
	  pl += 3;
	  tmp_ptr = I->tv;
	  vertexVals[pl] = tmp_ptr[0]; vertexVals[pl+1] = tmp_ptr[1]; vertexVals[pl+2] = tmp_ptr[2];
	  pl += 3;
	  for(b = I->Ns - 1; b >= 0; b--) {
	    tmp_ptr = I->tv + b*3;
	    vertexVals[pl] = tmp_ptr[0]; vertexVals[pl+1] = tmp_ptr[1]; vertexVals[pl+2] = tmp_ptr[2];
	    pl += 3;
	  }
	}
      }
#else
      if (ok){
	ok &= CGOBegin(cgo, GL_TRIANGLE_FAN);
	if (ok){
	  copy3f(I->n, v0);
	  invert3f(v0);
	  if(color_override)
	    ok &= CGOColorv(cgo, color_override);
	  else
	    ok &= CGOColorv(cgo, I->c);
	}
	if (ok)
	  ok &= CGOPickColor(cgo, I->i[0], cPickableAtom);
	if (ok)
	  ok &= CGONormalv(cgo, v0);
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
      }
#endif
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
#ifndef _PYMOL_CGO_DRAWARRAYS
      if (ok)
	ok &= CGOBegin(cgo, GL_TRIANGLE_FAN);
#endif
      if (ok){
	if(color_override)
	  ok &= CGOColorv(cgo, color_override);
	else
	  ok &= CGOColorv(cgo, I->c + 3 * (I->N - 1));
      }
      if (ok)
	ok &= CGOPickColor(cgo, I->i[I->N - 1], cPickableAtom);
      if (ok)
	ok &= CGONormalv(cgo, n);
#ifdef _PYMOL_CGO_DRAWARRAYS
      if (ok) {
	int nverts = I->Ns + 2, pl = 0;
	float *vertexVals, *tmp_ptr;
	vertexVals = CGODrawArrays(cgo, GL_TRIANGLE_FAN, CGO_VERTEX_ARRAY, nverts);      
	CHECKOK(ok, vertexVals);
	if (ok){
	  vertexVals[pl++] = v[0]; vertexVals[pl++] = v[1]; vertexVals[pl++] = v[2];
	  /* trace shape */
	  for(b = 0; b < I->Ns; b++) {
	    tmp_ptr = I->tv + b * 3;
	    vertexVals[pl++] = tmp_ptr[0]; vertexVals[pl++] = tmp_ptr[1]; vertexVals[pl++] = tmp_ptr[2];
	  }
	  tmp_ptr = I->tv;
	  vertexVals[pl++] = tmp_ptr[0]; vertexVals[pl++] = tmp_ptr[1]; vertexVals[pl++] = tmp_ptr[2];
	}
      }
#else
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
#endif
      break;
    case 2:
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

int ExtrudeCylindersToCGO(CExtrude * I, CGO * cgo, float tube_radius, short is_picking){
  float *v1, *c1, midv[3], midc[3];
  int a, *i;
  int ok = true;

  PRINTFD(I->G, FB_Extrude)
    " ExtrudeCylindersToCGO-DEBUG: entered.\n" ENDFD;

  v1 = I->p + 3;
  c1 = I->c + 3;
  i = I->i + 1;

  if (is_picking){
    float first = 2.f;
    for(a = 1; a < I->N; a++) {
      average3f(v1-3, v1, midv);
      average3f(c1-3, c1, midc);
      ok &= CGOPickColor(cgo, *(i-1), cPickableAtom);
      if (ok)
	ok &= CGOCustomCylinderv(cgo, v1-3, midv, tube_radius, c1-3, midc, first, 0.f);
      if (ok)
	ok &= CGOPickColor(cgo, *i, cPickableAtom);
      if (ok)
	ok &= CGOCustomCylinderv(cgo, midv, v1, tube_radius, midc, c1, 0.f, 2.f);
      v1 += 3;
      c1 += 3;
      i++;
      first = 0.f;
    }
  } else {
    if (I->N > 1){
      //      CGOSausage(cgo, v1-3, v1, tube_radius, c1-3, c1);
      ok &= CGOCustomCylinderv(cgo, v1-3, v1, tube_radius, c1-3, c1, 2.f, 2.f);
      v1 += 3;
      c1 += 3;
    }
    for(a = 2; ok && a < I->N; a++) {
      //      CGOSausage(cgo, v1-3, v1, tube_radius, c1-3, c1);
      ok &= CGOCustomCylinderv(cgo, v1-3, v1, tube_radius, c1-3, c1, 0.f, 2.f);
      v1 += 3;
      c1 += 3;
    }
  }

  PRINTFD(I->G, FB_Extrude)
    " ExtrudeCylindersToCGO-DEBUG: exiting...\n" ENDFD;
  return ok;
}

int ExtrudeCGOSurfaceVariableTube(CExtrude * I, CGO * cgo, int cap)
{
  int a, b, *i;
  float *v;
  float *n;
  float *c;
  float *sv, *sn, *tv, *tn, *tv1, *tn1, *TV = NULL, *TN = NULL, *AN = NULL, *an;
  float v0[3];
  float *sf;                    /* PUTTY: scale factor from ExtrudeMakeSausLUT() */
  int start, stop;
  int ok = true;
  PRINTFD(I->G, FB_Extrude)
    " ExtrudeCGOSurfaceTube-DEBUG: entered.\n" ENDFD;

  if(I->N && I->Ns) {

    TV = Alloc(float, 3 * (I->Ns + 1) * I->N);
    TN = Alloc(float, 3 * (I->Ns + 1) * I->N);
    AN = Alloc(float, 3 * I->N);        /* normals adjusted for changing widths */

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

    start = I->Ns / 4;
    stop = 3 * I->Ns / 4;
    for(b = 0; b < I->Ns; b++) {
#ifdef _PYMOL_CGO_DRAWARRAYS
      GLenum mode;
      if(SettingGetGlobal_i(I->G, cSetting_cartoon_debug) < 1.5)
	mode = GL_TRIANGLE_STRIP;
      else {
	mode = GL_LINE_STRIP;
      }
      c = I->c;
      i = I->i;
      if (ok ){
	int nverts = 2*I->N, pl = 0, plc = 0;
	float *vertexVals = NULL, *normalVals, *colorVals = 0, *nxtVals = 0, *tmp_ptr;
	float *pickColorVals;
	vertexVals = CGODrawArrays(cgo, mode, CGO_VERTEX_ARRAY | CGO_NORMAL_ARRAY | CGO_COLOR_ARRAY | CGO_PICK_COLOR_ARRAY, nverts);
	ok &= vertexVals ? 1 : 0;
	if (!ok)
	  break;

	normalVals = nxtVals = vertexVals + (nverts*3);
	colorVals = nxtVals = normalVals + (nverts*3);
	pickColorVals = (nxtVals + (nverts*4));
	
	for(a = 0; a < I->N; a++) {
	  tmp_ptr = c;
	  colorVals[plc++] = tmp_ptr[0]; colorVals[plc++] = tmp_ptr[1]; colorVals[plc++] = tmp_ptr[2];  colorVals[plc++] = cgo->alpha;
	  SetCGOPickColor(pickColorVals, nverts, pl, *i, cPickableAtom);	  
	  tmp_ptr = tn;
	  normalVals[pl] = tmp_ptr[0]; normalVals[pl+1] = tmp_ptr[1]; normalVals[pl+2] = tmp_ptr[2];
	  tmp_ptr = tv;
	  vertexVals[pl] = tmp_ptr[0]; vertexVals[pl+1] = tmp_ptr[1]; vertexVals[pl+2] = tmp_ptr[2];
	  pl += 3;
	  tn += 3;
	  tv += 3;
	  tmp_ptr = c;
	  colorVals[plc++] = tmp_ptr[0]; colorVals[plc++] = tmp_ptr[1]; colorVals[plc++] = tmp_ptr[2];  colorVals[plc++] = cgo->alpha;
	  SetCGOPickColor(pickColorVals, nverts, pl, *i, cPickableAtom);	  
	  tmp_ptr = tn1;
	  normalVals[pl] = tmp_ptr[0]; normalVals[pl+1] = tmp_ptr[1]; normalVals[pl+2] = tmp_ptr[2];
	  tmp_ptr = tv1;
	  vertexVals[pl] = tmp_ptr[0]; vertexVals[pl+1] = tmp_ptr[1]; vertexVals[pl+2] = tmp_ptr[2];
	  pl += 3;
	  tn1 += 3;
	  tv1 += 3;
	  c += 3;
	  i++;
	}
      }
#else
      if(SettingGetGlobal_i(I->G, cSetting_cartoon_debug) < 1.5)
        CGOBegin(cgo, GL_TRIANGLE_STRIP);
      else {
        CGOBegin(cgo, GL_LINE_STRIP);
      }
      c = I->c;
      i = I->i;
      for(a = 0; a < I->N; a++) {
        CGOColorv(cgo, c);
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
        i++;
      }
      CGOEnd(cgo);
#endif
    }

    if(ok && SettingGetGlobal_i(I->G, cSetting_cartoon_debug) > 3.5) {

      tv = TV;
      tn = TN;

      tv1 = TV + 3 * I->N;
      tn1 = TN + 3 * I->N;

      start = I->Ns / 4;
      stop = 3 * I->Ns / 4;
      for(b = 0; b < I->Ns; b++) {
        float vv[3];
#ifdef _PYMOL_CGO_DRAWARRAYS
	int nverts = 4*I->N, pl = 0;
	float *vertexVals, *normalVals, *tmp_ptr;
	vertexVals = CGODrawArrays(cgo, GL_LINES, CGO_VERTEX_ARRAY | CGO_NORMAL_ARRAY, nverts);
	ok &= vertexVals ? 1 : 0;
	if (ok){
	  normalVals = vertexVals + (nverts*3);
#else
	  CGOBegin(cgo, GL_LINES);
#endif
	  c = I->c;
	  i = I->i;
	  for(a = 0; a < I->N; a++) {
	    CGOColorv(cgo, c);
	    copy3f(tn, vv);
	    scale3f(vv, 0.3F, vv);
	    add3f(vv, tv, vv);
#ifdef _PYMOL_CGO_DRAWARRAYS
	    tmp_ptr = tn;
	    normalVals[pl] = tmp_ptr[0]; normalVals[pl+1] = tmp_ptr[1]; normalVals[pl+2] = tmp_ptr[2];
	    tmp_ptr = tv;
	    vertexVals[pl] = tmp_ptr[0]; vertexVals[pl+1] = tmp_ptr[1]; vertexVals[pl+2] = tmp_ptr[2];
	    pl += 3;
	    tmp_ptr = tn;
	    normalVals[pl] = tmp_ptr[0]; normalVals[pl+1] = tmp_ptr[1]; normalVals[pl+2] = tmp_ptr[2];
	    tmp_ptr = vv;
	    vertexVals[pl] = tmp_ptr[0]; vertexVals[pl+1] = tmp_ptr[1]; vertexVals[pl+2] = tmp_ptr[2];
	    pl += 3;
#else
	    CGONormalv(cgo, tn);
	    CGOVertexv(cgo, tv);
	    CGOVertexv(cgo, vv);
#endif
	    tn += 3;
	    tv += 3;
	    copy3f(tn1, vv);
	    scale3f(vv, 0.3F, vv);
	    add3f(vv, tv1, vv);
#ifdef _PYMOL_CGO_DRAWARRAYS
	    tmp_ptr = tn1;
	    normalVals[pl] = tmp_ptr[0]; normalVals[pl+1] = tmp_ptr[1]; normalVals[pl+2] = tmp_ptr[2];
	    tmp_ptr = tv1;
	    vertexVals[pl] = tmp_ptr[0]; vertexVals[pl+1] = tmp_ptr[1]; vertexVals[pl+2] = tmp_ptr[2];
	    pl += 3;
	    tmp_ptr = tn1;
	    normalVals[pl] = tmp_ptr[0]; normalVals[pl+1] = tmp_ptr[1]; normalVals[pl+2] = tmp_ptr[2];
	    tmp_ptr = vv;
	    vertexVals[pl] = tmp_ptr[0]; vertexVals[pl+1] = tmp_ptr[1]; vertexVals[pl+2] = tmp_ptr[2];
	    pl += 3;
#else
	    CGONormalv(cgo, tn1);
	    CGOVertexv(cgo, tv1);
	    CGOVertexv(cgo, vv);
#endif
	    tn1 += 3;
	    tv1 += 3;
	    c += 3;
	    i++;
	  }
#ifndef _PYMOL_CGO_DRAWARRAYS
	  CGOEnd(cgo);
#else
	}
#endif
      }
    }

    if(ok && cap) {

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
#ifndef _PYMOL_CGO_DRAWARRAYS
      CGOBegin(cgo, GL_TRIANGLE_FAN);
#endif
      copy3f(I->n, v0);
      invert3f(v0);
      CGOColorv(cgo, I->c);
      CGOPickColor(cgo, I->i[0], cPickableAtom);
      CGONormalv(cgo, v0);
#ifdef _PYMOL_CGO_DRAWARRAYS
      if (ok) {
	int nverts =  I->Ns + 2, pl = 0;
	float *vertexVals, *tmp_ptr;
	vertexVals = CGODrawArrays(cgo, GL_TRIANGLE_FAN, CGO_VERTEX_ARRAY, nverts);
	ok &= vertexVals ? 1 : 0;
	if (ok){
	  tmp_ptr = v;
	  vertexVals[pl++] = tmp_ptr[0]; vertexVals[pl++] = tmp_ptr[1]; vertexVals[pl++] = tmp_ptr[2];
	  tmp_ptr = I->tv;
	  vertexVals[pl++] = tmp_ptr[0]; vertexVals[pl++] = tmp_ptr[1]; vertexVals[pl++] = tmp_ptr[2];
	  for(b = I->Ns - 1; b >= 0; b--) {
	    tmp_ptr = I->tv + b * 3;
	    vertexVals[pl++] = tmp_ptr[0]; vertexVals[pl++] = tmp_ptr[1]; vertexVals[pl++] = tmp_ptr[2];
	  }
	}
#else
	CGOVertexv(cgo, v);
	/* trace shape */
	CGOVertexv(cgo, I->tv);
	for(b = I->Ns - 1; b >= 0; b--) {
	  CGOVertexv(cgo, I->tv + b * 3);
	}
#endif
#ifndef _PYMOL_CGO_DRAWARRAYS
	CGOEnd(cgo);
#endif
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
	
#ifndef _PYMOL_CGO_DRAWARRAYS
	CGOBegin(cgo, GL_TRIANGLE_FAN);
#endif
	CGOColorv(cgo, I->c + 3 * (I->N - 1));
	CGOPickColor(cgo, I->i[I->N - 1], cPickableAtom);
	CGONormalv(cgo, n);
#ifdef _PYMOL_CGO_DRAWARRAYS
	if (ok) {
	  int nverts = I->Ns + 2, pl = 0;
	  float *vertexVals, *tmp_ptr;
	  vertexVals = CGODrawArrays(cgo, GL_TRIANGLE_FAN, CGO_VERTEX_ARRAY, nverts);
	  ok &= vertexVals ? 1 : 0;
	  if (ok){
	    tmp_ptr = v;
	    vertexVals[pl++] = tmp_ptr[0]; vertexVals[pl++] = tmp_ptr[1]; vertexVals[pl++] = tmp_ptr[2];
	    for(b = 0; b < I->Ns; b++) {
	      tmp_ptr = I->tv + b *3;
	      vertexVals[pl++] = tmp_ptr[0]; vertexVals[pl++] = tmp_ptr[1]; vertexVals[pl++] = tmp_ptr[2];
	    }
	    tmp_ptr = I->tv;
	    vertexVals[pl++] = tmp_ptr[0]; vertexVals[pl++] = tmp_ptr[1]; vertexVals[pl++] = tmp_ptr[2];
	  }
	}
#else
	CGOVertexv(cgo, v);
	/* trace shape */
	for(b = 0; b < I->Ns; b++) {
	  CGOVertexv(cgo, I->tv + b * 3);
	}
	CGOVertexv(cgo, I->tv);
	CGOEnd(cgo);
#endif
      }
      FreeP(TV);
      FreeP(TN);
      FreeP(AN);
    }
    
    PRINTFD(I->G, FB_Extrude)
      " ExtrudeCGOSurfaceTube-DEBUG: exiting...\n" ENDFD;
  }
  return ok;
}

int ExtrudeCGOSurfacePolygon(CExtrude * I, CGO * cgo, int cap, float *color_override)
{
  int a, b, *i;
  float *v;
  float *n;
  float *c;
  float *sv, *sn, *tv, *tn, *tv1, *tn1, *TV, *TN;
  float v0[3];
  int ok = true;

  PRINTFD(I->G, FB_Extrude)
    " ExtrudeCGOSurfacePolygon-DEBUG: entered.\n" ENDFD;

  if(I->N && I->Ns) {

    TV = Alloc(float, 3 * (I->Ns + 1) * I->N);
    CHECKOK(ok, TV);
    if (ok)
      TN = Alloc(float, 3 * (I->Ns + 1) * I->N);
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
#ifdef _PYMOL_CGO_DRAWARRAYS
      int nverts = 2*I->N, pl = 0, plc = 0, damode = CGO_VERTEX_ARRAY | CGO_NORMAL_ARRAY | CGO_PICK_COLOR_ARRAY, nxtn = 3;
      float *vertexVals, *normalVals, *colorVals = 0, *nxtVals = 0, *tmp_ptr;
      float *pickColorVals;
      GLenum mode = GL_LINE_STRIP;
      if(SettingGetGlobal_i(I->G, cSetting_cartoon_debug) < 1.5)
	mode = GL_TRIANGLE_STRIP;
      if (color_override)
	ok &= CGOColorv(cgo, color_override);
      else 
	damode |= CGO_COLOR_ARRAY;
      if (ok)
	vertexVals = CGODrawArrays(cgo, mode, damode, nverts);
      CHECKOK(ok, vertexVals);
      if (!ok){
	break;
      }
      normalVals = nxtVals = vertexVals + (nverts*3);
      if (!color_override){
	colorVals = nxtVals = normalVals + (nverts*3);
	nxtn = 4;
      }
      pickColorVals = (nxtVals + (nverts*nxtn));
#else
      if(SettingGetGlobal_i(I->G, cSetting_cartoon_debug) < 1.5)
        ok &= CGOBegin(cgo, GL_TRIANGLE_STRIP);
      else {
        ok &= CGOBegin(cgo, GL_LINE_STRIP);
      }
#endif
      if(ok && color_override)
        ok &= CGOColorv(cgo, color_override);
      c = I->c;
      i = I->i;
      for(a = 0; ok && a < I->N; a++) {
#ifdef _PYMOL_CGO_DRAWARRAYS
	if(colorVals){
	  tmp_ptr = c;
	  colorVals[plc++] = tmp_ptr[0]; colorVals[plc++] = tmp_ptr[1]; colorVals[plc++] = tmp_ptr[2]; colorVals[plc++] = 1.f;
	}
	SetCGOPickColor(pickColorVals, nverts, pl, *i, cPickableAtom);
	tmp_ptr = tn;
	normalVals[pl] = tmp_ptr[0]; normalVals[pl+1] = tmp_ptr[1]; normalVals[pl+2] = tmp_ptr[2];
	tmp_ptr = tv;
	vertexVals[pl] = tmp_ptr[0]; vertexVals[pl+1] = tmp_ptr[1]; vertexVals[pl+2] = tmp_ptr[2];
	pl += 3;
#else
        if(!color_override)
          ok &= CGOColorv(cgo, c);
        if (ok)
	  ok &= CGOPickColor(cgo, *i, cPickableAtom);
	if (ok)
	  ok &= CGONormalv(cgo, tn);
	if (ok)
	  ok &= CGOVertexv(cgo, tv);
#endif
        tn += 3;
        tv += 3;
#ifdef _PYMOL_CGO_DRAWARRAYS
	if(colorVals){
	  tmp_ptr = c;
	  colorVals[plc++] = tmp_ptr[0]; colorVals[plc++] = tmp_ptr[1]; colorVals[plc++] = tmp_ptr[2]; colorVals[plc++] = 1.f;
	}
	SetCGOPickColor(pickColorVals, nverts, pl, *i, cPickableAtom);
	tmp_ptr = tn1;
	normalVals[pl] = tmp_ptr[0]; normalVals[pl+1] = tmp_ptr[1]; normalVals[pl+2] = tmp_ptr[2];
	tmp_ptr = tv1;
	vertexVals[pl] = tmp_ptr[0]; vertexVals[pl+1] = tmp_ptr[1]; vertexVals[pl+2] = tmp_ptr[2];
	pl += 3;
#else
	if (ok)
	  ok &= CGONormalv(cgo, tn1);
	if (ok)
	  ok &= CGOVertexv(cgo, tv1);
#endif
        tn1 += 3;
        tv1 += 3;
        c += 3;
        i++;
      }
      tv += 3 * I->N;
      tn += 3 * I->N;
      tv1 += 3 * I->N;
      tn1 += 3 * I->N;
#ifndef _PYMOL_CGO_DRAWARRAYS
      if (ok)
	ok &= CGOEnd(cgo);
#endif
    }

    if(ok && cap) {

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
#ifndef _PYMOL_CGO_DRAWARRAYS
      if (ok)
	ok &= CGOBegin(cgo, GL_TRIANGLE_FAN);
#endif
      if (ok){
	copy3f(I->n, v0);
	invert3f(v0);
	if(!color_override)
	  ok &= CGOColorv(cgo, I->c);
	if (ok)
	  ok &= CGOPickColor(cgo, I->i[0], cPickableAtom);
	if (ok)
	  ok &= CGONormalv(cgo, v0);
      }
#ifdef _PYMOL_CGO_DRAWARRAYS
      if (ok){
	int nverts = 2 + I->Ns, pl = 0;
	float *vertexVals, *tmp_ptr;
	vertexVals = CGODrawArrays(cgo, GL_TRIANGLE_FAN, CGO_VERTEX_ARRAY, nverts);	
	CHECKOK(ok, vertexVals);
	if (ok){
	  tmp_ptr = v;
	  vertexVals[pl] = tmp_ptr[0]; vertexVals[pl+1] = tmp_ptr[1]; vertexVals[pl+2] = tmp_ptr[2];
	  pl += 3;
	  tmp_ptr = I->tv;
	  vertexVals[pl] = tmp_ptr[0]; vertexVals[pl+1] = tmp_ptr[1]; vertexVals[pl+2] = tmp_ptr[2];
	  pl += 3;
	  for(b = I->Ns - 1; b >= 0; b--) {
	    tmp_ptr = I->tv + b * 3;
	    vertexVals[pl] = tmp_ptr[0]; vertexVals[pl+1] = tmp_ptr[1]; vertexVals[pl+2] = tmp_ptr[2];
	    pl += 3;
	  }
	}
      }
#else
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
#endif
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
#ifndef _PYMOL_CGO_DRAWARRAYS
      if (ok)
	ok &= CGOBegin(cgo, GL_TRIANGLE_FAN);
#endif
      if(ok && !color_override)
        ok &= CGOColorv(cgo, I->c + 3 * (I->N - 1));
      if (ok)
	ok &= CGOPickColor(cgo, I->i[I->N - 1], cPickableAtom);
      if (ok)
	ok &= CGONormalv(cgo, n);
#ifdef _PYMOL_CGO_DRAWARRAYS
      if (ok) {
	int nverts = I->Ns + 2, pl = 0;
	float *vertexVals, *tmp_ptr;
	vertexVals = CGODrawArrays(cgo, GL_TRIANGLE_FAN, CGO_VERTEX_ARRAY, nverts);      
	CHECKOK(ok, vertexVals);
	if (ok){
	  vertexVals[pl++] = v[0]; vertexVals[pl++] = v[1]; vertexVals[pl++] = v[2];
	  /* trace shape */
	  for(b = 0; b < I->Ns; b++) {
	    tmp_ptr = I->tv + b * 3;
	    vertexVals[pl++] = tmp_ptr[0]; vertexVals[pl++] = tmp_ptr[1]; vertexVals[pl++] = tmp_ptr[2];
	  }
	  tmp_ptr = I->tv;
	  vertexVals[pl++] = tmp_ptr[0]; vertexVals[pl++] = tmp_ptr[1]; vertexVals[pl++] = tmp_ptr[2];
	}
      }
#else
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
#endif
    }
    FreeP(TV);
    FreeP(TN);
  }

  PRINTFD(I->G, FB_Extrude)
    " ExtrudeCGOSurfacePolygon-DEBUG: exiting...\n" ENDFD;
  return ok;
}

int ExtrudeCGOSurfacePolygonTaper(CExtrude * I, CGO * cgo, int sampling,
				  float *color_override)
{
  int a, b, *i;
  float *v;
  float *n;
  float *c;
  float *sv, *sn, *tv, *tn, *tv1, *tn1, *TV = NULL, *TN = NULL;
  float s0[3];
  float f;
  int subN;
  int ok = true;

  subN = I->N - sampling;

  PRINTFD(I->G, FB_Extrude)
    " ExtrudeCGOSurfacePolygonTaper-DEBUG: entered.\n" ENDFD;

  if(I->N && I->Ns) {

    TV = Alloc(float, 3 * (I->Ns + 1) * I->N);
    CHECKOK(ok, TV);
    if (ok)
      TN = Alloc(float, 3 * (I->Ns + 1) * I->N);
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
#ifdef _PYMOL_CGO_DRAWARRAYS
      GLenum mode;
      if(SettingGetGlobal_i(I->G, cSetting_cartoon_debug) < 1.5)
	mode = GL_TRIANGLE_STRIP;
      else {
	mode = GL_LINE_STRIP;
      }
      c = I->c;
      i = I->i;
      if (ok){
	int nverts = 2*I->N, pl = 0, plc = 0, damode = CGO_VERTEX_ARRAY | CGO_NORMAL_ARRAY | CGO_PICK_COLOR_ARRAY, nxtn = 3;
	float *vertexVals, *normalVals, *colorVals = 0, *nxtVals = 0, *tmp_ptr;
	float *pickColorVals;
	if (color_override)
	  ok &= CGOColorv(cgo, color_override);
	else 
	  damode |= CGO_COLOR_ARRAY;
	if (ok)
	  vertexVals = CGODrawArrays(cgo, mode, damode, nverts);
	CHECKOK(ok, vertexVals);
	if (ok){
	  normalVals = nxtVals = vertexVals + (nverts*3);
	  if (damode & CGO_COLOR_ARRAY){
	    colorVals = nxtVals = normalVals + (nverts*3);
	    nxtn = 4;
	  }
	  pickColorVals = (nxtVals + (nverts*nxtn));
	  
	  for(a = 0; a < I->N; a++) {
	    if(colorVals){
	      tmp_ptr = c;
	      colorVals[plc++] = tmp_ptr[0]; colorVals[plc++] = tmp_ptr[1]; colorVals[plc++] = tmp_ptr[2]; colorVals[plc++] = cgo->alpha;
	    }
	    SetCGOPickColor(pickColorVals, nverts, pl, *i, cPickableAtom);	  
	    tmp_ptr = tn;
	    normalVals[pl] = tmp_ptr[0]; normalVals[pl+1] = tmp_ptr[1]; normalVals[pl+2] = tmp_ptr[2];
	    tmp_ptr = tv;
	    vertexVals[pl] = tmp_ptr[0]; vertexVals[pl+1] = tmp_ptr[1]; vertexVals[pl+2] = tmp_ptr[2];
	    pl += 3;
	    tn += 3;
	    tv += 3;
	    if(colorVals){
	      tmp_ptr = c;
	      colorVals[plc++] = tmp_ptr[0]; colorVals[plc++] = tmp_ptr[1]; colorVals[plc++] = tmp_ptr[2]; colorVals[plc++] = cgo->alpha;
	    }
	    SetCGOPickColor(pickColorVals, nverts, pl, *i, cPickableAtom);	  
	    tmp_ptr = tn1;
	    normalVals[pl] = tmp_ptr[0]; normalVals[pl+1] = tmp_ptr[1]; normalVals[pl+2] = tmp_ptr[2];
	    tmp_ptr = tv1;
	    vertexVals[pl] = tmp_ptr[0]; vertexVals[pl+1] = tmp_ptr[1]; vertexVals[pl+2] = tmp_ptr[2];
	    pl += 3;
	    tn1 += 3;
	    tv1 += 3;
	    c += 3;
	    i++;
	  }
	  tv += 3 * I->N;
	  tn += 3 * I->N;
	  tv1 += 3 * I->N;
	  tn1 += 3 * I->N;
	}
      }
#else
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
      i = I->i;
      for(a = 0; ok && a < I->N; a++) {
        if(!color_override)
          ok &= CGOColorv(cgo, c);
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
        i++;
      }
      if (ok){
	tv += 3 * I->N;
	tn += 3 * I->N;
	tv1 += 3 * I->N;
	tn1 += 3 * I->N;
	CGOEnd(cgo);
      }
#endif
    }

    FreeP(TV);
    FreeP(TN);
  }

  PRINTFD(I->G, FB_Extrude)
    " ExtrudeCGOSurfacePolygonTaper-DEBUG: exiting...\n" ENDFD;
  return ok;
}

int ExtrudeCGOSurfaceStrand(CExtrude * I, CGO * cgo, int sampling, float *color_override)
{
  int a, b, *i;
  float *v;
  float *n;
  float *c;
  float *sv, *sn, *tv, *tn, *tv1, *tn1, *TV, *TN;
  float v0[3], n0[3], s0[3], z[3] = { 1.0, 0.0, 1.0 };
  int subN;
  int ok = true;
  subN = I->N - sampling;

  PRINTFD(I->G, FB_Extrude)
    " ExtrudeCGOSurfaceStrand-DEBUG: entered.\n" ENDFD;

  if(I->N && I->Ns) {

    TV = Alloc(float, 3 * (I->Ns + 1) * I->N);
    CHECKOK(ok, TV);
    if (ok)
      TN = Alloc(float, 3 * (I->Ns + 1) * I->N);
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
#ifdef _PYMOL_CGO_DRAWARRAYS
      GLenum mode;
      if(SettingGetGlobal_i(I->G, cSetting_cartoon_debug) < 1.5)
	mode = GL_TRIANGLE_STRIP;
      else {
	mode = GL_LINE_STRIP;
      }
      c = I->c;
      i = I->i;
      {
	int nverts = 2*subN + 2, pl = 0, plc = 0;
	float *vertexVals, *normalVals, *colorVals = 0, *tmp_ptr;
	float *pickColorVals = 0;
	vertexVals = CGODrawArrays(cgo, mode, CGO_VERTEX_ARRAY | CGO_NORMAL_ARRAY | CGO_COLOR_ARRAY | CGO_PICK_COLOR_ARRAY, nverts);
	CHECKOK(ok, vertexVals);
	if (ok){
	  normalVals = vertexVals + (nverts*3);
	  colorVals = normalVals + (nverts*3);
	  pickColorVals = (colorVals + (nverts*4));
	  for(a = 0; a < I->N; a++) {
	    if(a <= subN) {
	      if(color_override && ((b == 2) || (b == 3) || (b == 6) || (b == 7)))
		tmp_ptr = color_override;
	      else
		tmp_ptr = c;
	      colorVals[plc++] = tmp_ptr[0]; colorVals[plc++] = tmp_ptr[1]; colorVals[plc++] = tmp_ptr[2]; colorVals[plc++] = cgo->alpha;
	      SetCGOPickColor(pickColorVals, nverts, pl, *i, cPickableAtom);
	      tmp_ptr = tn;
	      normalVals[pl] = tmp_ptr[0]; normalVals[pl+1] = tmp_ptr[1]; normalVals[pl+2] = tmp_ptr[2];
	      tmp_ptr = tv;
	      vertexVals[pl] = tmp_ptr[0]; vertexVals[pl+1] = tmp_ptr[1]; vertexVals[pl+2] = tmp_ptr[2];
	      pl += 3;
	    }
	    tn += 3;
	    tv += 3;
	    if(a <= subN) {
	      tmp_ptr = &colorVals[plc-4];
	      colorVals[plc++] = tmp_ptr[0]; colorVals[plc++] = tmp_ptr[1]; colorVals[plc++] = tmp_ptr[2]; colorVals[plc++] = cgo->alpha;
	      SetCGOPickColor(pickColorVals, nverts, pl, *i, cPickableAtom);	  
	      tmp_ptr = tn1;
	      normalVals[pl] = tmp_ptr[0]; normalVals[pl+1] = tmp_ptr[1]; normalVals[pl+2] = tmp_ptr[2];
	      tmp_ptr = tv1;
	      vertexVals[pl] = tmp_ptr[0]; vertexVals[pl+1] = tmp_ptr[1]; vertexVals[pl+2] = tmp_ptr[2];
	      pl += 3;
	    }
	    tn1 += 3;
	    tv1 += 3;
	    c += 3;
	    i++;
	  }
	}
	tv += 3 * I->N;
	tn += 3 * I->N;
	tv1 += 3 * I->N;
	tn1 += 3 * I->N;
      }
#else
      if(SettingGetGlobal_i(I->G, cSetting_cartoon_debug) < 1.5)
        ok &= CGOBegin(cgo, GL_TRIANGLE_STRIP);
      else {
        ok &= CGOBegin(cgo, GL_LINE_STRIP);
      }
      c = I->c;
      i = I->i;
      for(a = 0; ok &= a < I->N; a++) {
        if(a <= subN) {
	  if(color_override && ((b == 2) || (b == 3) || (b == 6) || (b == 7)))
	    ok &= CGOColorv(cgo, color_override);
	  else
	    ok &= CGOColorv(cgo, c);
	  if (ok)
	    ok &= CGOPickColor(cgo, *i, cPickableAtom);
	  if (ok)
	    ok &= CGONormalv(cgo, tn);
	  if (ok)
	    ok &= CGOVertexv(cgo, tv);
        }
        tn += 3;
        tv += 3;
        if(ok && a <= subN) {
          ok &= CGONormalv(cgo, tn1);
	  if (ok)
	    ok &= CGOVertexv(cgo, tv1);
        }
        tn1 += 3;
        tv1 += 3;
        c += 3;
        i++;
      }
      tv += 3 * I->N;
      tn += 3 * I->N;
      tv1 += 3 * I->N;
      tn1 += 3 * I->N;
      if (ok)
	ok &= CGOEnd(cgo);
#endif
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

#ifndef _PYMOL_CGO_DRAWARRAYS
      ok &= CGOBegin(cgo, GL_TRIANGLE_FAN);
#endif
      copy3f(I->n, v0);
      invert3f(v0);
      if (ok){
	if(color_override)
	  ok &= CGOColorv(cgo, color_override);
	else
	  ok &= CGOColorv(cgo, I->c);
      }
      if (ok)
	ok &= CGOPickColor(cgo, I->i[0], cPickableAtom);
      if (ok)
	ok &= CGONormalv(cgo, v0);
#ifdef _PYMOL_CGO_DRAWARRAYS
      if (ok) {
	int nverts = (I->Ns/2) + 2, pl = 0;
	float *vertexVals, *tmp_ptr;
	vertexVals = CGODrawArrays(cgo, GL_TRIANGLE_FAN, CGO_VERTEX_ARRAY, nverts);
	CHECKOK(ok, vertexVals);
	if (ok){
	  tmp_ptr = v;
	  vertexVals[pl++] = tmp_ptr[0]; vertexVals[pl++] = tmp_ptr[1]; vertexVals[pl++] = tmp_ptr[2];
	  tmp_ptr = I->tv;
	  vertexVals[pl++] = tmp_ptr[0]; vertexVals[pl++] = tmp_ptr[1]; vertexVals[pl++] = tmp_ptr[2];
	  /* trace shape */
	  for(b = I->Ns - 2; b >= 0; b -= 2) {
	    tmp_ptr = I->tv + b * 3;
	    vertexVals[pl++] = tmp_ptr[0]; vertexVals[pl++] = tmp_ptr[1]; vertexVals[pl++] = tmp_ptr[2];
	  }
	}
      }
#else
      if (ok)
	ok &= CGOVertexv(cgo, v);
      /* trace shape */
      if (ok)
	ok &= CGOVertexv(cgo, I->tv);
      for(b = I->Ns - 2; ok && b >= 0; b -= 2) {
        ok &= CGOVertexv(cgo, I->tv + b * 3);
      }
      if (ok)
	ok &= CGOEnd(cgo);
#endif
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
#ifdef _PYMOL_CGO_DRAWARRAYS
      GLenum mode;
      if(SettingGetGlobal_i(I->G, cSetting_cartoon_debug) < 1.5)
	mode = GL_TRIANGLE_STRIP;
      else {
	  mode = GL_LINE_STRIP;
      }
      c = I->c;
      i = I->i;
      {
	int nverts = 2*(I->N - subN + 1), pl = 0, plc = 0;
	float *vertexVals, *normalVals, *colorVals = 0, *tmp_ptr;
	float *pickColorVals;
	vertexVals = CGODrawArrays(cgo, mode, CGO_VERTEX_ARRAY | CGO_NORMAL_ARRAY | CGO_COLOR_ARRAY | CGO_PICK_COLOR_ARRAY, nverts);
	CHECKOK(ok, vertexVals);
	if (ok){
	  normalVals = vertexVals + (nverts*3);
	  colorVals = normalVals + (nverts*3);
	  pickColorVals = (colorVals + (nverts*4));
	  for(a = 0; a < I->N; a++) {
	    if(a >= (subN - 1)) {
	      if(color_override && ((b == 2) || (b == 3) || (b == 6) || (b == 7)))
		tmp_ptr = color_override;
	      else
		tmp_ptr = c;
	      colorVals[plc++] = tmp_ptr[0]; colorVals[plc++] = tmp_ptr[1]; colorVals[plc++] = tmp_ptr[2]; colorVals[plc++] = cgo->alpha;
	      SetCGOPickColor(pickColorVals, nverts, pl, *i, cPickableAtom);	  
	      tmp_ptr = tn;
	      normalVals[pl] = tmp_ptr[0]; normalVals[pl+1] = tmp_ptr[1]; normalVals[pl+2] = tmp_ptr[2];
	      tmp_ptr = tv;
	      vertexVals[pl] = tmp_ptr[0]; vertexVals[pl+1] = tmp_ptr[1]; vertexVals[pl+2] = tmp_ptr[2];
	      pl += 3;
	    }
	    tn += 3;
	    tv += 3;
	    if (a >= (subN - 1)) {
	      tmp_ptr = &colorVals[plc-4];
	      colorVals[plc++] = tmp_ptr[0]; colorVals[plc++] = tmp_ptr[1]; colorVals[plc++] = tmp_ptr[2]; colorVals[plc++] = cgo->alpha;
	      SetCGOPickColor(pickColorVals, nverts, pl, *i, cPickableAtom);	  
	      tmp_ptr = tn1;
	      normalVals[pl] = tmp_ptr[0]; normalVals[pl+1] = tmp_ptr[1]; normalVals[pl+2] = tmp_ptr[2];
	      tmp_ptr = tv1;
	      vertexVals[pl] = tmp_ptr[0]; vertexVals[pl+1] = tmp_ptr[1]; vertexVals[pl+2] = tmp_ptr[2];
	      pl += 3;
	    }
	    tn1 += 3;
	    tv1 += 3;
	    c += 3;
	    i++;
	  }
	  tv += 3 * I->N;
	  tn += 3 * I->N;
	  tv1 += 3 * I->N;
	  tn1 += 3 * I->N;
	}
      }
#else
      if (ok){
	if(SettingGetGlobal_i(I->G, cSetting_cartoon_debug) < 1.5)
	  ok &= CGOBegin(cgo, GL_TRIANGLE_STRIP);
	else {
	  ok &= CGOBegin(cgo, GL_LINE_STRIP);
	}
      }
      c = I->c;
      i = I->i;
      for(a = 0; ok && a < I->N; a++) {
        if(a >= (subN - 1)) {
          if(color_override && ((b == 2) || (b == 3) || (b == 6) || (b == 7)))
            ok &= CGOColorv(cgo, color_override);
          else
            ok &= CGOColorv(cgo, c);
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
        i++;
      }
      tv += 3 * I->N;
      tn += 3 * I->N;
      tv1 += 3 * I->N;
      tn1 += 3 * I->N;
      if (ok)
	ok &= CGOEnd(cgo);
#endif
    }

    n = I->n + 9 * (subN - 1);
    v = I->p + 3 * (subN - 1);

    sv = I->sv;
    tv = I->tv;

    if (ok){
      for(b = 0; b < I->Ns; b++) {
	copy3f(sv, s0);
	s0[2] = s0[2] * 1.5F;
	transform33Tf3f(n, s0, tv);
	add3f(v, tv, tv);
	sv += 3;
	tv += 3;
      }
    }

#ifndef _PYMOL_CGO_DRAWARRAYS
    if (ok)
      ok &= CGOBegin(cgo, GL_TRIANGLE_FAN);
#endif
    copy3f(n, v0);
    invert3f(v0);
    if (ok){
      if(color_override)
	ok &= CGOColorv(cgo, color_override);
      else
	ok &= CGOColorv(cgo, I->c + 3 * (subN - 1));
    }
    if (ok)
      ok &= CGOPickColor(cgo, I->i[(subN - 1)], cPickableAtom);
    if (ok)
      ok &= CGONormalv(cgo, v0);
#ifdef _PYMOL_CGO_DRAWARRAYS
    if (ok) {
      int nverts = (I->Ns/2) + 2, pl = 0;
      float *vertexVals, *tmp_ptr;
      vertexVals = CGODrawArrays(cgo, GL_TRIANGLE_FAN, CGO_VERTEX_ARRAY, nverts);
      CHECKOK(ok, vertexVals);
      if (ok){
	tmp_ptr = v;
	vertexVals[pl++] = tmp_ptr[0]; vertexVals[pl++] = tmp_ptr[1]; vertexVals[pl++] = tmp_ptr[2];
	
	tv = I->tv;
	for(b = 0; b < I->Ns; b += 2) {
	  tmp_ptr = I->tv + b * 3;
	  vertexVals[pl++] = tmp_ptr[0]; vertexVals[pl++] = tmp_ptr[1]; vertexVals[pl++] = tmp_ptr[2];
	}
	tmp_ptr = I->tv;
	vertexVals[pl++] = tmp_ptr[0]; vertexVals[pl++] = tmp_ptr[1]; vertexVals[pl++] = tmp_ptr[2];
      }
    }
#else
    if (ok)
      ok &= CGOVertexv(cgo, v);
    /* trace shape */
    tv = I->tv;
    for(b = 0; ok && b < I->Ns; b += 2) {
      ok &= CGOVertexv(cgo, I->tv + b * 3);
    }
    if (ok)
      ok &= CGOVertexv(cgo, I->tv);
    if (ok)
      ok &= CGOEnd(cgo);
#endif
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
  int *i;
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
        *sf = 0.0F;
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
      float *SF = Alloc(float, I->N);
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
    FreeP(I->i);
    FreeP(I->sf);               /* PUTTY */
    I->p = Alloc(float, 3 * (n + 1));
    CHECKOK(ok, I->p);
    if (ok)
      I->n = Alloc(float, 9 * (n + 1));
    CHECKOK(ok, I->n);
    if (ok)
      I->c = Alloc(float, 3 * (n + 1));
    CHECKOK(ok, I->c);
    if (ok)
      I->i = Alloc(int, 3 * (n + 1));
    CHECKOK(ok, I->i);
    if (ok)
      I->sf = Alloc(float, n + 1);        /* PUTTY: scale factors */
    CHECKOK(ok, I->sf);
    if (!ok){
      FreeP(I->p);
      FreeP(I->n);
      FreeP(I->c);
      FreeP(I->i);
      FreeP(I->sf);
      I->p = NULL;
      I->n = NULL;
      I->c = NULL;
      I->i = NULL;
      I->sf = NULL;
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
  FreeP(I->tn);
  FreeP(I->tv);
  FreeP(I->sn);
  FreeP(I->sv);
  FreeP(I->i);
  FreeP(I->sf);
  OOFreeP(I);
}
