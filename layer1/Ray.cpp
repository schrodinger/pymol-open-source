
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
-*   Larry Coopet (various optimizations)
-*   Chris Want (RayRenderVRML2, via the public domain )
-*
Z* -------------------------------------------------------------------
*/
#include"os_python.h"

#include"os_predef.h"
#include"os_std.h"

#include"Base.h"
#include"MemoryDebug.h"
#include"Err.h"
#include"Vector.h"
#include"Setting.h"
#include"Ortho.h"
#include"Util.h"
#include"Ray.h"
#include"Triangle.h"
#include"Color.h"
#include"Matrix.h"
#include"P.h"
#include"MemoryCache.h"
#include"Character.h"
#include"Text.h"
#include"PyMOL.h"
#include"Scene.h"
#include"PConv.h"
#include"MyPNG.h"
#include"CGO.h"
#include "Feedback.h"

#define SettingGetfv SettingGetGlobal_3fv

#include"Basis.h"

#ifndef RAY_SMALL
#define RAY_SMALL 0.00001
#endif

/* BASES 
   0 contains untransformed vertices (vector size = 3)
	1 contains transformed vertices (vector size = 3)
	2 contains transformed vertices for shadowing 
*/


/* note: the following value must be at least one greater than the max
   number of lights */
#define MAX_BASIS 12

typedef float float3[3];
typedef float float4[4];

struct _CRayThreadInfo {
  CRay *ray;
  int width, height;
  unsigned int *image;
  float front, back;
  unsigned int fore_mask;
  float *bkrd_top, *bkrd_bottom;
  short bkrd_is_gradient; /* if not gradient, use bkrd_top as bkrd */
  float ambient;
  unsigned int background;
  int border;
  int phase, n_thread;
  int x_start, x_stop;
  int y_start, y_stop;
  unsigned int *edging;
  unsigned int edging_cutoff;
  int perspective;
  float fov, pos[3];
  float *depth;
  
  int bgWidth, bgHeight;
  void *bkrd_data; /* used for image-based background */
};

struct _CRayHashThreadInfo {
  CBasis *basis;
  int *vert2prim;
  CPrimitive *prim;
  int n_prim;
  float *clipBox;
  unsigned int *image;
  unsigned int background;
  size_t bytes;
  int perspective;
  float front;
  int phase;
  float size_hint;
  CRay *ray;
  float *bkrd_top, *bkrd_bottom;
  short bkrd_is_gradient; /* if not gradient, use bkrd_top as bkrd */
  int width, height;
  int opaque_back;
};

struct _CRayAntiThreadInfo {
  unsigned int *image;
  unsigned int *image_copy;
  unsigned int width, height;
  int mag;
  int phase, n_thread;
  CRay *ray;
};

static
void RayRelease(CRay * I);

static void RayApplyMatrix33(unsigned int n, float3 * q, const float m[16], float3 * p);
static void RayApplyMatrixInverse33(unsigned int n, float3 * q, const float m[16], float3 * p);
static void RayTransformNormals33(unsigned int n, float3 * q, const float m[16], float3 * p);
static void RayTransformInverseNormals33(unsigned int n, float3 * q, const float m[16],
                                  float3 * p);
static
void RayProjectTriangle(CRay * I, RayInfo * r, float *light, float *v0, float *n0,
                        float scale);
void RaySetContext(CRay * I, pymol::RenderContext context)
{
  I->context = context;
}

float RayGetScreenVertexScale(CRay * I, float *v1)
{
  /* what size should a screen pixel be at the coordinate provided? */

  float vt[3];
  float ratio;
  RayApplyMatrix33(1, (float3 *) vt, I->ModelView, (float3 *) v1);

  if(I->Ortho) {
    ratio =
      2 * (float) (fabs(I->Pos.z) * tan((I->Fov / 2.0) * cPI / 180.0)) / (I->Height);
  } else {
    float front_size =
      2 * I->Volume[4] * ((float) tan((I->Fov / 2.0F) * PI / 180.0F)) / (I->Height);
    ratio = fabs(front_size * (-vt[2] / I->Volume[4]));
  }
  return ratio;
}

static void RayApplyContextToVertex(CRay * I, float *v)
{
  switch (I->context) {
  case pymol::RenderContext::UnitWindow:
    {
      float tw;
      float th;

      if(I->AspRatio > 1.0F) {
        tw = I->AspRatio;
        th = 1.0F;
      } else {
        th = 1.0F / I->AspRatio;
        tw = 1.0F;
      }

      if(!SettingGetGlobal_b(I->G, cSetting_ortho)) {
        float scale = v[2] + 0.5F;
        scale = I->FrontBackRatio * scale + 1.0F - scale;

        /* z-coodinate is easy... */

        v[2] = v[2] * I->Range[2] - (I->Volume[4] + I->Volume[5]) / 2.0F;
        v[0] -= 0.5F;
        v[1] -= 0.5F;
        v[0] = scale * v[0] * I->Range[0] / tw + (I->Volume[0] + I->Volume[1]) / 2.0F;
        v[1] = scale * v[1] * I->Range[1] / th + (I->Volume[2] + I->Volume[3]) / 2.0F;

        RayApplyMatrixInverse33(1, (float3 *) v, I->ModelView, (float3 *) v);
      } else {
        v[0] += (tw - 1.0F) / 2;
        v[1] += (th - 1.0F) / 2;
        v[0] = v[0] * (I->Range[0] / tw) + I->Volume[0];
        v[1] = v[1] * (I->Range[1] / th) + I->Volume[2];
        v[2] = v[2] * I->Range[2] - (I->Volume[4] + I->Volume[5]) / 2.0F;
        RayApplyMatrixInverse33(1, (float3 *) v, I->ModelView, (float3 *) v);
      }

    }
    break;
  }
}

static void RayApplyContextToNormal(CRay * I, float *v)
{
  switch (I->context) {
  case pymol::RenderContext::UnitWindow:
    RayTransformInverseNormals33(1, (float3 *) v, I->ModelView, (float3 *) v);
    break;
  }
}

int RayGetNPrimitives(CRay * I)
{
  return (I->NPrimitive);
}


/*========================================================================*/

static void RayGetSphereNormal(CRay * I, RayInfo * r)
{

  r->impact[0] = r->base[0];
  r->impact[1] = r->base[1];
  r->impact[2] = r->base[2] - r->dist;

  r->surfnormal[0] = r->impact[0] - r->sphere[0];
  r->surfnormal[1] = r->impact[1] - r->sphere[1];
  r->surfnormal[2] = r->impact[2] - r->sphere[2];

  normalize3f(r->surfnormal);
}

static void RayGetSphereNormalPerspective(CRay * I, RayInfo * r)
{

  r->impact[0] = r->base[0] + r->dist * r->dir[0];
  r->impact[1] = r->base[1] + r->dist * r->dir[1];
  r->impact[2] = r->base[2] + r->dist * r->dir[2];

  r->surfnormal[0] = r->impact[0] - r->sphere[0];
  r->surfnormal[1] = r->impact[1] - r->sphere[1];
  r->surfnormal[2] = r->impact[2] - r->sphere[2];

  normalize3f(r->surfnormal);
}

static void fill(unsigned int *buffer, unsigned int value, size_t cnt)
{
  // std::fill_n(buffer, cnt, value);
  for (auto it_end = buffer + cnt; buffer != it_end;) {
    *(buffer++) = value;
  }
}
#define MIN(x,y) ((x < y) ? x : y)
#define MAX(x,y) ((x > y) ? x : y)

static void add4ucf(unsigned char *ucval, float *fval, float mulv){
  fval[0] += ucval[0] * mulv;
  fval[1] += ucval[1] * mulv;
  fval[2] += ucval[2] * mulv;
  fval[3] += ucval[3] * mulv;
}

static void copy4fuc(float *fval, unsigned char *ucval){
  ucval[0] = (0xFF & (int)pymol_roundf(fval[0]));
  ucval[1] = (0xFF & (int)pymol_roundf(fval[1]));
  ucval[2] = (0xFF & (int)pymol_roundf(fval[2]));
  ucval[3] = (0xFF & (int)pymol_roundf(fval[3]));
}

static float fmodpos(float a, float b){
  float ret = fmod(a, b);
  if (ret < 0.f){
    ret = fmod(-ret, b);
    ret = fmod( b - ret, b);
  }
  return ret;
}

static void compute_background_for_pixel(unsigned char *bkrd, short isOutsideInY,
    int bg_image_mode, const float *bg_image_tilesize, const float *bg_rgb, int bg_image_linear,
    unsigned char *bkgrd_data, int bkgrd_width, int bkgrd_height, float w, float wr,
    float hl, float hpixelx, short blend_bg)
{
  short isBackground = isOutsideInY;
  float wl;
  switch (bg_image_mode){
  case 1: // isCentered
    {
      float tmpx = floor(w - hpixelx);
      isBackground |= (tmpx < 0.f || tmpx > (float)bkgrd_width);
      wl = fmodpos(tmpx, (float)bkgrd_width);
    }
    break;
  case 2: // isTiled
    wl = bkgrd_width * (fmodpos((float)w, bg_image_tilesize[0])/bg_image_tilesize[0]);
    break;
  case 3: // isCenteredRepeated
    wl = fmodpos(floor(w - hpixelx), (float)bkgrd_width) ;
    break;
  default:
    wl = w * wr;
  }
  if (isBackground){
    copy3f(bg_rgb, bkrd);
    bkrd[3] = 1.f;
  } else if (bg_image_linear){
    int wlf = floor(wl), hlf = floor(hl);
    float xp = (wl - wlf) - .5f, yp = (hl - hlf) - .5f;
    float xpa = fabs(xp), ypa = fabs(yp);
    float xpam = 1.f - xpa, ypam = 1.f - ypa;
    float valx[4] = { 0.f, 0.f, 0.f, 0.f }, valy[4] = { 0.f, 0.f, 0.f, 0.f }, val[4] = { 0.f, 0.f, 0.f, 0.f };
    int xd = (xp > 0.f) ? 1 : -1, yd = (yp > 0.f) ? 1 : -1;
    int wlfn = (wlf + xd + bkgrd_width) % bkgrd_width, hlfn = (hlf + yd + bkgrd_height) % bkgrd_height;
    
    add4ucf(&bkgrd_data[(bkgrd_width*hlf + wlf)*4], valx, xpam);
    add4ucf(&bkgrd_data[(bkgrd_width*hlf + wlfn)*4], valx, xpa);
    
    add4ucf(&bkgrd_data[(bkgrd_width*hlfn + wlf)*4], valy, xpam);
    add4ucf(&bkgrd_data[(bkgrd_width*hlfn + wlfn)*4], valy, xpa);
    
    mult4f(valx, ypam, val);
    mult4f(valy, ypa, valy);
    add4f(valy, val, val);
    copy4fuc(val, bkrd);
  } else {
    copy4f(&bkgrd_data[(bkgrd_width*(int)floor(hl) + (int)floor(wl))*4], bkrd);
  }
  // alpha blending with bg_rgb if needed
  if (blend_bg && bkrd[3] < 255){
    float alpha = (255.f - bkrd[3]) / 255.f;
    float tmpadd1[3], tmpadd2[3];
    tmpadd1[0] = bkrd[0]; tmpadd1[1] = bkrd[1]; tmpadd1[2] = bkrd[2];
    mult3f(tmpadd1, 1.f - alpha, tmpadd1);	
    mult3f(bg_rgb, alpha, tmpadd2);
    add3f(tmpadd1, tmpadd2, tmpadd2);
    bkrd[0] = 0xFF & (int)pymol_roundf(tmpadd2[0]);
    bkrd[1] = 0xFF & (int)pymol_roundf(tmpadd2[1]);
    bkrd[2] = 0xFF & (int)pymol_roundf(tmpadd2[2]);
    bkrd[3] = 255;
  }
}

static void fill_background_image(CRay * I, unsigned int *buffer, int width, int height, unsigned int cnt){
  int bkgrd_width = I->bkgrd_data->getWidth(), bkgrd_height = I->bkgrd_data->getHeight();
  unsigned char *bkgrd_data = I->bkgrd_data->bits();
  int bg_image_mode = SettingGetGlobal_i(I->G, cSetting_bg_image_mode);
  int bg_image_linear = SettingGetGlobal_b(I->G, cSetting_bg_image_linear);
  float bg_rgb[3];
  const float *tmpf;
  int w, h;
  unsigned int value;
  float wr = bkgrd_width/(float)width, hr = bkgrd_height/(float)height;
  float hl;
  unsigned int back_mask;
  unsigned char bkrd[4];
  float hpixelx = floor(width/2.f)  - floor(bkgrd_width / 2.f),
    hpixely = floor(height/2.f) - floor(bkgrd_height / 2.f);
  auto bg_image_tilesize = SettingGet<const float*>(I->G, cSetting_bg_image_tilesize);
  short isOutsideInY = 0;
  int opaque_back ;
  opaque_back = SettingGetGlobal_i(I->G, cSetting_ray_opaque_background);
  if(opaque_back < 0)
    opaque_back = SettingGetGlobal_i(I->G, cSetting_opaque_background);

  tmpf = ColorGet(I->G, SettingGet_color(I->G, nullptr, nullptr, cSetting_bg_rgb));
  mult3f(tmpf, 255.f, bg_rgb);

  if(opaque_back) {
    if(I->BigEndian)
      back_mask = 0x000000FF;
    else
      back_mask = 0xFF000000;
  } else {
    back_mask = 0x00000000;
  }
  // tiled or stretched
  for (h=0; h<height; h++){
    switch (bg_image_mode){
    case 1: // isCentered
      {
	float tmpy = floor(h - hpixely);
	isOutsideInY = (tmpy < 0.f || tmpy > (float)bkgrd_height);
	hl = fmodpos(tmpy, (float)bkgrd_height);
      }
      break;
    case 2: // isTiled
      hl = bkgrd_height * (fmodpos((float)h, bg_image_tilesize[1])/bg_image_tilesize[1]);
      break;
    case 3: // isCenteredRepeated
      hl = fmodpos(floor(h - hpixely), (float)bkgrd_height);
      break;
    default:
      hl = h * hr;
      break;
    }
    for (w=0; w<width; w++){
      compute_background_for_pixel(bkrd, isOutsideInY, bg_image_mode, bg_image_tilesize, bg_rgb, bg_image_linear, bkgrd_data, bkgrd_width, bkgrd_height, w, wr, hl, hpixelx, opaque_back );
      if(I->BigEndian){
	value = 
	  ((0xFF & bkrd[0]) << 24) |
	  ((0xFF & bkrd[1]) << 16) |
	  ((0xFF & bkrd[2]) << 8) |
	  (0xFF & bkrd[3]);
      } else {
	value = 
	  ((0xFF & bkrd[3]) << 24) |
	  ((0xFF & bkrd[2]) << 16) |
	  ((0xFF & bkrd[1]) << 8) |
	  (0xFF & bkrd[0]);
      }
      if(opaque_back) {
	value = back_mask | value;
      }
      *(buffer++) = value;
    }
  }
}

static void fill_gradient(CRay * I, int opaque_back, unsigned int *buffer, float *bkrd_bottom, float *bkrd_top, int width, int height, unsigned int cnt)
{
  const float _p499 = 0.499F;
  int w, h;
  unsigned int value;
  float bkrd[3], perc;
  unsigned int back_mask;

  if(opaque_back) {
    if(I->BigEndian)
      back_mask = 0x000000FF;
    else
      back_mask = 0xFF000000;
  } else {
    back_mask = 0x00000000;
  }
  for (h=0; h<height; h++){
    /* for fill_gradient, y is from top to bottom */
    perc = h/(float)height;
    bkrd[0] = bkrd_top[0] + perc * (bkrd_bottom[0] - bkrd_top[0]);
    bkrd[1] = bkrd_top[1] + perc * (bkrd_bottom[1] - bkrd_top[1]);
    bkrd[2] = bkrd_top[2] + perc * (bkrd_bottom[2] - bkrd_top[2]);
    if(I->BigEndian){
      value = back_mask | ((0xFF & ((unsigned int) (bkrd[0] * 255 + _p499))) << 24) |
	((0xFF & ((unsigned int) (bkrd[1] * 255 + _p499))) << 16) |
	((0xFF & ((unsigned int) (bkrd[2] * 255 + _p499))) << 8);
    } else {
      value = back_mask | ((0xFF & ((unsigned int) (bkrd[2] * 255 + _p499))) << 16) |
	((0xFF & ((unsigned int) (bkrd[1] * 255 + _p499))) << 8) |
	((0xFF & ((unsigned int) (bkrd[0] * 255 + _p499))));
    }
    for (w=0; w<width; w++){
      *(buffer++) = value;
    }
  }
}

/*========================================================================*/
static void RayReflectAndTexture(CRay * I, RayInfo * r, int perspective)
{
  if(r->prim->wobble)
    switch (r->prim->wobble) {
    case 1:
      scatter3f(r->surfnormal, I->WobbleParam[0]);
      break;
    case 2:
      wiggle3f(r->surfnormal, r->impact, I->WobbleParam);
      break;
    case 3:
      {
        float3 v;
        float3 n;
        copy3f(r->impact, v);
        RayApplyMatrixInverse33(1, &v, I->ModelView, &v);
        n[0] = (float) cos((v[0] + v[1] + v[2]) * I->WobbleParam[1]);
        n[1] = (float) cos((v[0] - v[1] + v[2]) * I->WobbleParam[1]);
        n[2] = (float) cos((v[0] + v[1] - v[2]) * I->WobbleParam[1]);
        RayTransformNormals33(1, &n, I->ModelView, &n);
        scale3f(n, I->WobbleParam[0], n);
        add3f(n, r->surfnormal, r->surfnormal);
        normalize3f(r->surfnormal);
      }
    case 4:
      {
        float3 v;
        float3 n;
        float *tp = I->WobbleParam;
        copy3f(r->impact, v);
        RayApplyMatrixInverse33(1, &v, I->ModelView, &v);
        n[0] = I->Random[0xFF & (int) ((cos((v[0]) * tp[1]) * 256 * tp[2]))];
        n[1] = I->Random[0xFF & (int) ((cos((v[1]) * tp[1]) * 256 * tp[2] + 96))];
        n[2] = I->Random[0xFF & (int) ((cos((v[2]) * tp[1]) * 256 * tp[2] + 148))];
        RayTransformNormals33(1, &n, I->ModelView, &n);
        scale3f(n, tp[0], n);
        add3f(n, r->surfnormal, r->surfnormal);
        normalize3f(r->surfnormal);
      }
      break;
    case 5:
      {
        float3 v;
        float3 n;
        float *tp = I->WobbleParam;
        copy3f(r->impact, v);
        RayApplyMatrixInverse33(1, &v, I->ModelView, &v);
        n[0] = I->Random[0xFF & (int) ((v[0] * tp[1]) + 0)] +
          I->Random[0xFF & (int) ((v[1] * tp[1]) + 20)] +
          I->Random[0xFF & (int) ((v[2] * tp[1]) + 40)];
        n[1] = I->Random[0xFF & (int) ((-v[0] * tp[1]) + 90)] +
          I->Random[0xFF & (int) ((v[1] * tp[1]) + 100)] +
          I->Random[0xFF & (int) ((-v[2] * tp[1]) + 120)];
        n[2] = I->Random[0xFF & (int) ((v[0] * tp[1]) + 200)] +
          I->Random[0xFF & (int) ((-v[1] * tp[1]) + 70)] +
          I->Random[0xFF & (int) ((v[2] * tp[1]) + 30)];

        n[0] +=
          I->Random[0xFF & ((int) ((v[0] - v[1]) * tp[1]) + 0)] +
          I->Random[0xFF & ((int) ((v[1] - v[2]) * tp[1]) + 20)] +
          I->Random[0xFF & ((int) ((v[2] - v[0]) * tp[1]) + 40)];
        n[1] +=
          I->Random[0xFF & ((int) ((v[0] + v[1]) * tp[1]) + 10)] +
          I->Random[0xFF & ((int) ((v[1] + v[2]) * tp[1]) + 90)] +
          I->Random[0xFF & ((int) ((v[2] + v[0]) * tp[1]) + 30)];
        n[2] +=
          I->Random[0xFF & ((int) ((-v[0] + v[1]) * tp[1]) + 220)] +
          I->Random[0xFF & ((int) ((-v[1] + v[2]) * tp[1]) + 20)] +
          I->Random[0xFF & ((int) ((-v[2] + v[0]) * tp[1]) + 50)];

        n[0] +=
          I->Random[0xFF & ((int) ((v[0] + v[1] + v[2]) * tp[1]) + 5)] +
          I->Random[0xFF & ((int) ((v[0] + v[1] + v[2]) * tp[1]) + 25)] +
          I->Random[0xFF & ((int) ((v[0] + v[1] + v[2]) * tp[1]) + 46)];
        n[1] +=
          I->Random[0xFF & ((int) ((-v[0] - v[1] + v[2]) * tp[1]) + 90)] +
          I->Random[0xFF & ((int) ((-v[0] - v[1] + v[2]) * tp[1]) + 45)] +
          I->Random[0xFF & ((int) ((-v[0] - v[1] + v[2]) * tp[1]) + 176)];
        n[2] +=
          I->Random[0xFF & ((int) ((v[0] + v[1] - v[2]) * tp[1]) + 192)] +
          I->Random[0xFF & ((int) ((v[0] + v[1] - v[2]) * tp[1]) + 223)] +
          I->Random[0xFF & ((int) ((v[0] + v[1] - v[2]) * tp[1]) + 250)];

        RayTransformNormals33(1, &n, I->ModelView, &n);
        scale3f(n, tp[0], n);
        add3f(n, r->surfnormal, r->surfnormal);
        normalize3f(r->surfnormal);
      }
      break;
    }
  if(perspective) {
    r->dotgle = dot_product3f(r->dir, r->surfnormal);
    r->flat_dotgle = -r->dotgle;

    r->reflect[0] = r->dir[0] - (2 * r->dotgle * r->surfnormal[0]);
    r->reflect[1] = r->dir[1] - (2 * r->dotgle * r->surfnormal[1]);
    r->reflect[2] = r->dir[2] - (2 * r->dotgle * r->surfnormal[2]);
  } else {
    r->dotgle = -r->surfnormal[2];
    r->flat_dotgle = r->surfnormal[2];

    r->reflect[0] = -(2 * r->dotgle * r->surfnormal[0]);
    r->reflect[1] = -(2 * r->dotgle * r->surfnormal[1]);
    r->reflect[2] = -1.0F - (2 * r->dotgle * r->surfnormal[2]);
  }
}


/*========================================================================*/
int RayExpandPrimitives(CRay * I)
{
  int a;
  float *v0, *v1, *n0, *n1;
  CBasis *basis;
  int nVert, nNorm;
  float voxel_floor;
  int ok = true;

  nVert = 0;
  nNorm = 0;
  for(a = 0; a < I->NPrimitive; a++) {
    switch (I->Primitive[a].type) {
    case cPrimSphere:
      nVert++;
      break;
    case cPrimEllipsoid:
      nVert++;
      nNorm += 3;
      break;
    case cPrimCone:
    case cPrimCylinder:
    case cPrimSausage:
      nVert++;
      nNorm++;
      break;
    case cPrimTriangle:
    case cPrimCharacter:
      nVert += 3;
      nNorm += 4;
      break;
    }
  }

  basis = I->Basis;

  VLACacheSize(I->G, basis->Vertex, float, 3 * nVert, 0, cCache_basis_vertex);
  VLACacheSize(I->G, basis->Radius, float, nVert, 0, cCache_basis_radius);
  VLACacheSize(I->G, basis->Radius2, float, nVert, 0, cCache_basis_radius2);
  VLACacheSize(I->G, basis->Vert2Normal, int, nVert, 0, cCache_basis_vert2normal);
  VLACacheSize(I->G, basis->Normal, float, 3 * nNorm, 0, cCache_basis_normal);

  I->Vert2Prim.resize(nVert);

  voxel_floor = I->PixelRadius / 2.0F;

  basis->MaxRadius = 0.0F;
  basis->MinVoxel = 0.0F;
  basis->NVertex = nVert;
  basis->NNormal = nNorm;

  nVert = 0;
  nNorm = 0;
  v0 = basis->Vertex;
  n0 = basis->Normal;
  ok &= !I->G->Interrupt;
  for(a = 0; ok && a < I->NPrimitive; a++) {
    switch (I->Primitive[a].type) {
    case cPrimTriangle:
    case cPrimCharacter:

      I->Primitive[a].vert = nVert;
      I->Vert2Prim[nVert] = a;
      I->Vert2Prim[nVert + 1] = a;
      I->Vert2Prim[nVert + 2] = a;
      basis->Radius[nVert] = I->Primitive[a].r1;
      basis->Radius2[nVert] = I->Primitive[a].r1 * I->Primitive[a].r1;  /*necessary?? */
      /*              if(basis->Radius[nVert]>basis->MinVoxel)
         basis->MinVoxel=basis->Radius[nVert]; */
      if(basis->MinVoxel < voxel_floor)
        basis->MinVoxel = voxel_floor;
      basis->Vert2Normal[nVert] = nNorm;
      basis->Vert2Normal[nVert + 1] = nNorm;
      basis->Vert2Normal[nVert + 2] = nNorm;
      n1 = I->Primitive[a].n0;
      (*n0++) = (*n1++);
      (*n0++) = (*n1++);
      (*n0++) = (*n1++);
      n1 = I->Primitive[a].n1;
      (*n0++) = (*n1++);
      (*n0++) = (*n1++);
      (*n0++) = (*n1++);
      n1 = I->Primitive[a].n2;
      (*n0++) = (*n1++);
      (*n0++) = (*n1++);
      (*n0++) = (*n1++);
      n1 = I->Primitive[a].n3;
      (*n0++) = (*n1++);
      (*n0++) = (*n1++);
      (*n0++) = (*n1++);
      nNorm += 4;
      v1 = I->Primitive[a].v1;
      (*v0++) = (*v1++);
      (*v0++) = (*v1++);
      (*v0++) = (*v1++);
      v1 = I->Primitive[a].v2;
      (*v0++) = (*v1++);
      (*v0++) = (*v1++);
      (*v0++) = (*v1++);
      v1 = I->Primitive[a].v3;
      (*v0++) = (*v1++);
      (*v0++) = (*v1++);
      (*v0++) = (*v1++);
      nVert += 3;
      break;
    case cPrimSphere:
      I->Primitive[a].vert = nVert;
      I->Vert2Prim[nVert] = a;
      v1 = I->Primitive[a].v1;
      basis->Radius[nVert] = I->Primitive[a].r1;
      basis->Radius2[nVert] = I->Primitive[a].r1 * I->Primitive[a].r1;  /*precompute */
      if(basis->Radius[nVert] > basis->MaxRadius)
        basis->MaxRadius = basis->Radius[nVert];
      (*v0++) = (*v1++);
      (*v0++) = (*v1++);
      (*v0++) = (*v1++);
      nVert++;
      break;
    case cPrimEllipsoid:
      I->Primitive[a].vert = nVert;
      I->Vert2Prim[nVert] = a;
      v1 = I->Primitive[a].v1;
      basis->Radius[nVert] = I->Primitive[a].r1;
      basis->Radius2[nVert] = I->Primitive[a].r1 * I->Primitive[a].r1;  /*precompute */
      if(basis->Radius[nVert] > basis->MaxRadius)
        basis->MaxRadius = basis->Radius[nVert];
      basis->Vert2Normal[nVert] = nNorm;
      (*v0++) = (*v1++);
      (*v0++) = (*v1++);
      (*v0++) = (*v1++);
      nVert++;
      n1 = I->Primitive[a].n1;
      (*n0++) = (*n1++);
      (*n0++) = (*n1++);
      (*n0++) = (*n1++);
      n1 = I->Primitive[a].n2;
      (*n0++) = (*n1++);
      (*n0++) = (*n1++);
      (*n0++) = (*n1++);
      n1 = I->Primitive[a].n3;
      (*n0++) = (*n1++);
      (*n0++) = (*n1++);
      (*n0++) = (*n1++);
      nNorm += 3;
      break;
    case cPrimCone:
    case cPrimCylinder:
    case cPrimSausage:
      I->Primitive[a].vert = nVert;
      I->Vert2Prim[nVert] = a;
      basis->Radius[nVert] = I->Primitive[a].r1;
      basis->Radius2[nVert] = I->Primitive[a].r1 * I->Primitive[a].r1;  /*precompute */
      if(basis->MinVoxel < voxel_floor)
        basis->MinVoxel = voxel_floor;
      subtract3f(I->Primitive[a].v2, I->Primitive[a].v1, n0);
      I->Primitive[a].l1 = (float) length3f(n0);
      normalize3f(n0);
      n0 += 3;
      basis->Vert2Normal[nVert] = nNorm;
      nNorm++;
      v1 = I->Primitive[a].v1;
      (*v0++) = (*v1++);
      (*v0++) = (*v1++);
      (*v0++) = (*v1++);
      nVert++;
      break;
    }
    ok &= !I->G->Interrupt;
  }
  if(nVert > basis->NVertex) {
    fprintf(stderr, "Error: basis->NVertex exceeded\n");
  }
  PRINTFB(I->G, FB_Ray, FB_Blather)
    " Ray: minvoxel  %8.3f\n Ray: NPrimit  %d nvert %d\n", basis->MinVoxel, I->NPrimitive,
    nVert ENDFB(I->G);
  return ok;
}


/*========================================================================*/
void RayComputeBox(CRay * I)
{

#define minmax(v,r) { \
  xp = v[0] + r;\
  xm = v[0] - r;\
  yp = v[1] + r;\
  ym = v[1] - r;\
  zp = v[2] + r;\
  zm = v[2] - r;\
  if(xmin>xm) xmin = xm;\
  if(xmax<xp) xmax = xp;\
  if(ymin>ym) ymin = ym;\
  if(ymax<yp) ymax = yp;\
  if(zmin>zm) zmin = zm;\
  if(zmax<zp) zmax = zp;\
}

  CPrimitive *prm;
  CBasis *basis1;

  float xmin = 0.0F, ymin = 0.0F, xmax = 0.0F, ymax = 0.0F, zmin = 0.0F, zmax = 0.0F;
  float xp, xm, yp, ym, zp, zm;

  float *v, r;
  float vt[3];
  const float _0 = 0.0F;
  int a;

  basis1 = I->Basis + 1;
  if(basis1->NVertex) {
    xmin = xmax = basis1->Vertex[0];
    ymin = ymax = basis1->Vertex[1];
    zmin = zmax = basis1->Vertex[2];

    for(a = 0; a < I->NPrimitive; a++) {
      prm = I->Primitive + a;

      switch (prm->type) {
      case cPrimTriangle:
      case cPrimCharacter:

        r = _0;
        v = basis1->Vertex + prm->vert * 3;
        minmax(v, r);
        v = basis1->Vertex + prm->vert * 3 + 3;
        minmax(v, r);
        v = basis1->Vertex + prm->vert * 3 + 6;
        minmax(v, r);
        break;
      case cPrimSphere:
      case cPrimEllipsoid:
        r = prm->r1;
        v = basis1->Vertex + prm->vert * 3;
        minmax(v, r);
        break;
      case cPrimCone:
      case cPrimCylinder:
      case cPrimSausage:
        r = prm->r1;
        v = basis1->Vertex + prm->vert * 3;
        minmax(v, r);
        v = basis1->Normal + basis1->Vert2Normal[prm->vert] * 3;
        scale3f(v, prm->l1, vt);
        v = basis1->Vertex + prm->vert * 3;
        add3f(v, vt, vt);
        minmax(vt, r);
        break;
      }                         /* end of switch */
    }
  }
  // without the R_SMALL4 padding, this caused ray tracing failure
  // on OS X (X11) with a flat ObjectCGO (zmin == zmax).
  I->min_box[0] = xmin - R_SMALL4;
  I->min_box[1] = ymin - R_SMALL4;
  I->min_box[2] = zmin - R_SMALL4;
  I->max_box[0] = xmax + R_SMALL4;
  I->max_box[1] = ymax + R_SMALL4;
  I->max_box[2] = zmax + R_SMALL4;
}

int RayTransformFirst(CRay * I, int perspective, int identity)
{
  CBasis *basis0, *basis1;
  CPrimitive *prm;
  int a;
  float *v0;
  int backface_cull;
  int ok = true;
  int two_sided_lighting = SettingGetGlobal_b(I->G, cSetting_two_sided_lighting);

  if(two_sided_lighting<0) {
    if(SettingGetGlobal_i(I->G, cSetting_surface_cavity_mode))
      two_sided_lighting = true;
    else
      two_sided_lighting = false;
  }

  backface_cull = SettingGetGlobal_b(I->G, cSetting_backface_cull);

  if(two_sided_lighting ||
     (SettingGetGlobal_i(I->G, cSetting_transparency_mode) == 1) ||
     (SettingGetGlobal_i(I->G, cSetting_ray_interior_color) != -1) || I->CheckInterior)
    backface_cull = 0;


  basis0 = I->Basis;
  basis1 = I->Basis + 1;

  if (ok){
    VLACacheSize(I->G, basis1->Vertex, float, 3 * basis0->NVertex, 1, cCache_basis_vertex);
    CHECKOK(ok, basis1->Vertex);
  }
  if (ok){
    VLACacheSize(I->G, basis1->Normal, float, 3 * basis0->NNormal, 1, cCache_basis_normal);
    CHECKOK(ok, basis1->Normal);
  }
  if (ok){
    VLACacheSize(I->G, basis1->Precomp, float, 3 * basis0->NNormal, 1,
		 cCache_basis_precomp);
    CHECKOK(ok, basis1->Precomp);
  }
  if (ok){
    VLACacheSize(I->G, basis1->Vert2Normal, int, basis0->NVertex, 1,
		 cCache_basis_vert2normal);
    CHECKOK(ok, basis1->Vert2Normal);
  }
  if (ok){
    VLACacheSize(I->G, basis1->Radius, float, basis0->NVertex, 1, cCache_basis_radius);
    CHECKOK(ok, basis1->Radius);
  }
  if (ok){
    VLACacheSize(I->G, basis1->Radius2, float, basis0->NVertex, 1, cCache_basis_radius2);
    CHECKOK(ok, basis1->Radius2);
  }
  ok &= !I->G->Interrupt;
  if (ok){
    if(identity) {
      UtilCopyMem(basis1->Vertex, basis0->Vertex, basis0->NVertex * sizeof(float) * 3);
    } else {
      RayApplyMatrix33(basis0->NVertex, (float3 *) basis1->Vertex,
		       I->ModelView, (float3 *) basis0->Vertex);
    }
  }
  ok &= !I->G->Interrupt;

  if (ok){
    memcpy(basis1->Radius, basis0->Radius, basis0->NVertex * sizeof(float));
    memcpy(basis1->Radius2, basis0->Radius2, basis0->NVertex * sizeof(float));
    memcpy(basis1->Vert2Normal, basis0->Vert2Normal, basis0->NVertex * sizeof(int));
  }
  ok &= !I->G->Interrupt;
  if (ok){
    basis1->MaxRadius = basis0->MaxRadius;
    basis1->MinVoxel = basis0->MinVoxel;
    basis1->NVertex = basis0->NVertex;
  }
  ok &= !I->G->Interrupt;
  if (ok){
    if(identity) {
      UtilCopyMem(basis1->Normal, basis0->Normal, basis0->NNormal * sizeof(float) * 3);
    } else {
      RayTransformNormals33(basis0->NNormal, (float3 *) basis1->Normal,
			    I->ModelView, (float3 *) basis0->Normal);
    }
    basis1->NNormal = basis0->NNormal;
  }
  ok &= !I->G->Interrupt;
  if(perspective) {
    for(a = 0; ok && a < I->NPrimitive; a++) {
      prm = I->Primitive + a;

      prm = I->Primitive + a;
      switch (prm->type) {
      case cPrimTriangle:
      case cPrimCharacter:
        BasisTrianglePrecomputePerspective(basis1->Vertex + prm->vert * 3,
                                           basis1->Vertex + prm->vert * 3 + 3,
                                           basis1->Vertex + prm->vert * 3 + 6,
                                           basis1->Precomp +
                                           basis1->Vert2Normal[prm->vert] * 3);
        break;
      }
      ok &= !I->G->Interrupt;
    }
  } else {
    for(a = 0; ok && a < I->NPrimitive; a++) {
      prm = I->Primitive + a;
      switch (prm->type) {
      case cPrimTriangle:
      case cPrimCharacter:
        BasisTrianglePrecompute(basis1->Vertex + prm->vert * 3,
                                basis1->Vertex + prm->vert * 3 + 3,
                                basis1->Vertex + prm->vert * 3 + 6,
                                basis1->Precomp + basis1->Vert2Normal[prm->vert] * 3);
        v0 = basis1->Normal + (basis1->Vert2Normal[prm->vert] * 3 + 3);
        prm->cull = (!identity) && backface_cull && ((v0[2] < 0.0F) && (v0[5] < 0.0F)
                                                     && (v0[8] < 0.0F));
        break;
      case cPrimCone:
      case cPrimSausage:
      case cPrimCylinder:
        BasisCylinderSausagePrecompute(basis1->Normal +
                                       basis1->Vert2Normal[prm->vert] * 3,
                                       basis1->Precomp +
                                       basis1->Vert2Normal[prm->vert] * 3);
        break;

      }
      ok &= !I->G->Interrupt;
    }
  }
  return ok;
}


/*========================================================================*/
static int RayTransformBasis(CRay * I, CBasis * basis1, int group_id)
{
  CBasis *basis0;
  int a;
  float *v0, *v1;
  CPrimitive *prm;
  int ok = true;

  basis0 = I->Basis + 1;

  VLACacheSize(I->G, basis1->Vertex, float, 3 * basis0->NVertex, group_id,
               cCache_basis_vertex);
  CHECKOK(ok, basis1->Vertex);
  if (ok)
    VLACacheSize(I->G, basis1->Normal, float, 3 * basis0->NNormal, group_id,
		 cCache_basis_normal);
  CHECKOK(ok, basis1->Normal);
  if (ok)
    VLACacheSize(I->G, basis1->Precomp, float, 3 * basis0->NNormal, group_id,
		 cCache_basis_precomp);
  CHECKOK(ok, basis1->Precomp);
  if (ok)
    VLACacheSize(I->G, basis1->Vert2Normal, int, basis0->NVertex, group_id,
		 cCache_basis_vert2normal);
  CHECKOK(ok, basis1->Vert2Normal);
  if (ok)
    VLACacheSize(I->G, basis1->Radius, float, basis0->NVertex, group_id,
		 cCache_basis_radius);
  CHECKOK(ok, basis1->Radius);
  if (ok)
    VLACacheSize(I->G, basis1->Radius2, float, basis0->NVertex, group_id,
		 cCache_basis_radius2);
  CHECKOK(ok, basis1->Radius2);
  v0 = basis0->Vertex;
  v1 = basis1->Vertex;
  for(a = 0; ok && a < basis0->NVertex; a++) {
    matrix_transform33f3f(basis1->Matrix, v0, v1);
    v0 += 3;
    v1 += 3;
    basis1->Radius[a] = basis0->Radius[a];
    basis1->Radius2[a] = basis0->Radius2[a];
    basis1->Vert2Normal[a] = basis0->Vert2Normal[a];
  }
  if (ok){
    v0 = basis0->Normal;
    v1 = basis1->Normal;
  }
  for(a = 0; ok && a < basis0->NNormal; a++) {
    matrix_transform33f3f(basis1->Matrix, v0, v1);
    v0 += 3;
    v1 += 3;
  }
  if (ok){
    basis1->MaxRadius = basis0->MaxRadius;
    basis1->MinVoxel = basis0->MinVoxel;
    basis1->NVertex = basis0->NVertex;
    basis1->NNormal = basis0->NNormal;
  }
  for(a = 0; ok && a < I->NPrimitive; a++) {
    prm = I->Primitive + a;
    switch (prm->type) {
    case cPrimTriangle:
    case cPrimCharacter:

      BasisTrianglePrecompute(basis1->Vertex + prm->vert * 3,
                              basis1->Vertex + prm->vert * 3 + 3,
                              basis1->Vertex + prm->vert * 3 + 6,
                              basis1->Precomp + basis1->Vert2Normal[prm->vert] * 3);
      break;
    case cPrimCone:
    case cPrimSausage:
    case cPrimCylinder:
      BasisCylinderSausagePrecompute(basis1->Normal + basis1->Vert2Normal[prm->vert] * 3,
                                     basis1->Precomp +
                                     basis1->Vert2Normal[prm->vert] * 3);
      break;

    }
  }
  return ok;
}


/*========================================================================*/
void RayRenderTest(CRay * I, int width, int height, float front, float back, float fov)
{

  PRINTFB(I->G, FB_Ray, FB_Details)
    " RayRenderTest: obtained %i graphics primitives.\n", I->NPrimitive ENDFB(I->G);
}


/*========================================================================*/

G3dPrimitive *RayRenderG3d(CRay * I, int width, int height,
                           float front, float back, float fov, int quiet)
{
  /* generate a rendering stream for Miguel's G3d java rendering engine */

  float scale_x, scale_y;
  int shift_x, shift_y;
  float *d;
  CBasis *base;
  CPrimitive *prim;
  float *vert;
  float vert2[3];
  int a;
  G3dPrimitive *jprim = VLAlloc(G3dPrimitive, 10000), *jp;
  int n_jp = 0;

#define convert_r(r) 2*(int)(r*scale_x);
#define convert_x(x) shift_x + (int)(x*scale_x);
#define convert_y(y) height - (shift_y + (int)(y*scale_y));
#define convert_z(z) -(int)((z+front)*scale_x);
#define convert_col(c) (0xFF000000 | (((int)(c[0]*255.0))<<16) | (((int)(c[1]*255.0))<<8) | (((int)(c[2]*255.0))))

  RayExpandPrimitives(I);
  RayTransformFirst(I, 0, false);

  if(!quiet) {
    PRINTFB(I->G, FB_Ray, FB_Details)
      " RayRenderG3d: processed %i graphics primitives.\n", I->NPrimitive ENDFB(I->G);
  }
  base = I->Basis + 1;

  /* always orthoscopic */

  /* front should give a zero Z, 
     -I->Range[0] should be off the right hand size
     I->Range[1] should be off the top */
  scale_x = width / I->Range[0];
  scale_y = height / I->Range[1];
  shift_x = width / 2;
  shift_y = height / 2;

  for(a = 0; a < I->NPrimitive; a++) {
    prim = I->Primitive + a;
    vert = base->Vertex + 3 * (prim->vert);
    switch (prim->type) {
    case cPrimSphere:
      VLACheck(jprim, G3dPrimitive, n_jp);
      jp = jprim + n_jp;
      jp->op = 1;
      jp->r = convert_r(prim->r1);
      jp->x1 = convert_x(vert[0]);
      jp->y1 = convert_y(vert[1]);
      jp->z1 = convert_z(vert[2]);
      jp->c = convert_col(prim->c1);
      n_jp++;
      break;
    case cPrimSausage:
      VLACheck(jprim, G3dPrimitive, n_jp);
      d = base->Normal + 3 * base->Vert2Normal[prim->vert];
      scale3f(d, prim->l1, vert2);
      add3f(vert, vert2, vert2);

      jp = jprim + n_jp;
      jp->op = 3;
      jp->r = convert_r(prim->r1);
      jp->x1 = convert_x(vert[0]);
      jp->y1 = convert_y(vert[1]);
      jp->z1 = convert_z(vert[2]);
      jp->x2 = convert_x(vert2[0]);
      jp->y2 = convert_y(vert2[1]);
      jp->z2 = convert_z(vert2[2]);
      jp->c = convert_col(prim->c1);
      n_jp++;
      break;
    case cPrimTriangle:
      VLACheck(jprim, G3dPrimitive, n_jp);
      jp = jprim + n_jp;
      jp->op = 2;
      jp->x1 = convert_x(vert[0]);
      jp->y1 = convert_y(vert[1]);
      jp->z1 = convert_z(vert[2]);
      jp->x2 = convert_x(vert[3]);
      jp->y2 = convert_y(vert[4]);
      jp->z2 = convert_z(vert[5]);
      jp->x3 = convert_x(vert[6]);
      jp->y3 = convert_y(vert[7]);
      jp->z3 = convert_z(vert[8]);
      jp->c = convert_col(prim->c1);
      n_jp++;
      break;
    }
  }
  VLASize(jprim, G3dPrimitive, n_jp);
  return jprim;
}

void RayRenderVRML1(CRay * I, int width, int height,
                    char **vla_ptr, float front, float back,
                    float fov, float angle, float z_corr)
{
  char *vla = *vla_ptr;
  ov_size cc = 0;               /* character count */
  OrthoLineType buffer;

  RayExpandPrimitives(I);
  RayTransformFirst(I, 0, false);

  strcpy(buffer, "#VRML V1.0 ascii\n\n");
  UtilConcatVLA(&vla, &cc, buffer);

  UtilConcatVLA(&vla, &cc, "MaterialBinding { value OVERALL }\n");

  sprintf(buffer,
          "Material {\n ambientColor 0 0 0\n diffuseColor 1 1 1\n specularColor 1 1 1\nshininess 0.2\n}\n");
  UtilConcatVLA(&vla, &cc, buffer);

  {
    int a;
    CPrimitive *prim;
    float *vert;
    CBasis *base = I->Basis + 1;

    UtilConcatVLA(&vla, &cc, "Separator {\n");

    UtilConcatVLA(&vla, &cc, "MatrixTransform {\n");
    UtilConcatVLA(&vla, &cc, "matrix 1.0 0.0 0.0 0.0\n");
    UtilConcatVLA(&vla, &cc, "       0.0 1.0 0.0 0.0\n");
    UtilConcatVLA(&vla, &cc, "       0.0 0.0 1.0 0.0\n");
    sprintf(buffer, "    %8.6f %8.6f %8.6f 1.0\n",
            (I->Volume[0] + I->Volume[1]) / 2, (I->Volume[2] + I->Volume[3]) / 2, 0.0F);
    UtilConcatVLA(&vla, &cc, buffer);
    UtilConcatVLA(&vla, &cc, "}\n");

    for(a = 0; a < I->NPrimitive; a++) {
      prim = I->Primitive + a;
      vert = base->Vertex + 3 * (prim->vert);
      switch (prim->type) {
      case cPrimSphere:
        sprintf(buffer,
                "Material {\ndiffuseColor %6.4f %6.4f %6.4f\n}\n\n",
                prim->c1[0], prim->c1[1], prim->c1[2]);
        UtilConcatVLA(&vla, &cc, buffer);
        UtilConcatVLA(&vla, &cc, "Separator {\n");
        sprintf(buffer,
                "Transform {\ntranslation %8.6f %8.6f %8.6f\nscaleFactor %8.6f %8.6f %8.6f\n}\n",
                vert[0], vert[1], vert[2] - z_corr, prim->r1, prim->r1, prim->r1);
        UtilConcatVLA(&vla, &cc, buffer);
        sprintf(buffer, "Sphere {}\n");
        UtilConcatVLA(&vla, &cc, buffer);
        UtilConcatVLA(&vla, &cc, "}\n\n");
        break;
      case cPrimCylinder:
      case cPrimSausage:
      case cPrimCone:
        break;
      case cPrimTriangle:
        break;
      }
    }

    UtilConcatVLA(&vla, &cc, "}\n");
  }

  *vla_ptr = vla;
}

int TriangleReverse(CPrimitive * p)
{
  float s1[3], s2[3], n0[3];

  subtract3f(p->v1, p->v2, s1);
  subtract3f(p->v3, p->v2, s2);
  cross_product3f(s1, s2, n0);

  if(dot_product3f(p->n0, n0) < 0.0F)
    return 0;
  else
    return 1;
}

void RayRenderVRML2(CRay * I, int width, int height,
                    char **vla_ptr, float front, float back,
                    float fov, float angle, float z_corr)
{

  /* 

     From: pymol-users-admin@lists.sourceforge.net on behalf of Chris Want
     Sent: Tuesday, February 07, 2006 1:47 PM
     To: pymol-users@lists.sourceforge.net
     Subject: [PyMOL] VRML patch

     Hi Warren,

     I took your advice and modified the RayRenderVRML2() function to
     support triangles. I also threw out the sphere code that was already
     there and rewrote it (the code there was for VRML1, not
     VRML2). While I was at it, I also implemented export for cylinders
     and sausages.

     The code in the attached patch (diff-ed against cvs, and tested with
     two VRML2 readers) can be regarded as being in the public domain.

     Regards,
     Chris

     cwant_at_ualberta.ca

   */

  char *vla = *vla_ptr;
  ov_size cc = 0;               /* character count */
  OrthoLineType buffer;
  float mid[3];                 /*, wid[3]; */
  float h_fov = cPI * (fov * width) / (180 * height);
  int identity = (SettingGetGlobal_i(I->G, cSetting_geometry_export_mode) == 1);

  RayExpandPrimitives(I);
  RayTransformFirst(I, 0, identity);
  RayComputeBox(I);
  /*
     mid[0] = (I->max_box[0] + I->min_box[0]) / 2.0;
     mid[1] = (I->max_box[1] + I->min_box[1]) / 2.0;
     mid[2] = (I->max_box[2] + I->min_box[2]) / 2.0;
     wid[0]  = (I->max_box[0] - I->min_box[0]);
     wid[1]  = (I->max_box[1] - I->min_box[1]);
     wid[2]  = (I->max_box[2] - I->min_box[2]);
   */
  auto Pos = glm::value_ptr(I->Pos);
  copy3f(Pos, mid);
  UtilConcatVLA(&vla, &cc, "#VRML V2.0 utf8\n"  /* WLD: most VRML2 readers req. utf8 */
                "\n");
  if(!identity) {
    sprintf(buffer, "Viewpoint {\n" " position 0 0 %6.8f\n" " orientation 1 0 0 0\n" " description \"Z view\"\n" " fieldOfView %8.6f\n" /* WLD: use correct FOV */
            "}\n"
            /* WLD: only write the viewpoint which matches PyMOL
               "Viewpoint {\n"
               " position %6.8f 0 0\n"
               " orientation 0 1 0 1.570796\n"
               " description \"X view\"\n"
               "}\n"
               "Viewpoint {\n"
               " position 0 %6.8f 0\n"
               " orientation 0 -0.707106 -0.7071061 3.141592\n"
               " description \"Y view\"\n"
               "}\n" */ ,
            -z_corr,            /* *0.96646  for some reason, PyMOL and C4D cameras differ by about 3.5% ... */
            h_fov
            /*(wid[2] + wid[1]),
               (wid[0] + wid[1]),
               (wid[1] + wid[2]) */
      );
    UtilConcatVLA(&vla, &cc, buffer);
  }
  if(!identity) {
    float light[3];
    auto lightv = SettingGet<const float*>(I->G, cSetting_light);
    copy3f(lightv, light);
    normalize3f(light);
    sprintf(buffer,
            "DirectionalLight {\n"
            " direction %8.6f %8.6f %8.3f\n" "}\n", light[0], light[1], light[2]);
    UtilConcatVLA(&vla, &cc, buffer);
  }
  UtilConcatVLA(&vla, &cc,
                "NavigationInfo {\n" " headlight TRUE\n" " type \"EXAMINE\"\n" "}\n");
  {
    int a, b;
    CPrimitive *prim;
    float *vert;
    int mesh_obj = false, mesh_start = 0;

    CBasis *base = I->Basis + 1;

    for(a = 0; a < I->NPrimitive; a++) {
      prim = I->Primitive + a;
      vert = base->Vertex + 3 * (prim->vert);

      if(prim->type == cPrimTriangle) {
        if(!mesh_obj) {
          /* start mesh */
          mesh_start = a;
          UtilConcatVLA(&vla, &cc,
                        "Shape {\n"
                        " appearance Appearance {\n"
                        "  material Material { diffuseColor 1.0 1.0 1.0 }\n"
                        " }\n"
                        " geometry IndexedFaceSet {\n"
                        "  coord Coordinate {\n" "   point [\n");
          mesh_obj = true;
        }
      } else if(mesh_obj) {
        CPrimitive *cprim;
        int tri = 0;
        /* output connectivity */
        UtilConcatVLA(&vla, &cc, "   ]\n" "  }\n" "  coordIndex [\n");
        for(b = mesh_start; b < a; b++) {
          cprim = I->Primitive + b;
          if(TriangleReverse(cprim))
            sprintf(buffer, "%d %d %d -1,\n", tri, tri + 2, tri + 1);
          else
            sprintf(buffer, "%d %d %d -1,\n", tri, tri + 1, tri + 2);
          UtilConcatVLA(&vla, &cc, buffer);
          tri += 3;
        }

        /* output vertex colors */
        UtilConcatVLA(&vla, &cc,
                      "  ]\n"
                      "  colorPerVertex TRUE\n" "  color Color {\n" "   color [\n");
        for(b = mesh_start; b < a; b++) {
          cprim = I->Primitive + b;
          sprintf(buffer,
                  "%6.4f %6.4f %6.4f,\n"
                  "%6.4f %6.4f %6.4f,\n"
                  "%6.4f %6.4f %6.4f,\n",
                  cprim->c1[0], cprim->c1[1], cprim->c1[2],
                  cprim->c2[0], cprim->c2[1], cprim->c2[2],
                  cprim->c3[0], cprim->c3[1], cprim->c3[2]);
          UtilConcatVLA(&vla, &cc, buffer);
        }

        /* output vertex normals */
        UtilConcatVLA(&vla, &cc,
                      "  ] } \n"
                      "  normalPerVertex TRUE\n" "  normal Normal {\n" "   vector [\n");
        for(b = mesh_start; b < a; b++) {
          cprim = I->Primitive + b;
          {
            float *norm = base->Normal + 3 * base->Vert2Normal[cprim->vert];
            sprintf(buffer, "%6.4f %6.4f %6.4f,\n" "%6.4f %6.4f %6.4f,\n" "%6.4f %6.4f %6.4f,\n", norm[3], norm[4], norm[5],    /* transformed cprim->n1 */
                    norm[6], norm[7], norm[8],  /* transformed cprim->n2 */
                    norm[9], norm[10], norm[11]);       /* transformed cprim->n3 */
            UtilConcatVLA(&vla, &cc, buffer);
          }
        }
        UtilConcatVLA(&vla, &cc, "  ] }\n" "  normalIndex [ \n");
        tri = 0;
        for(b = mesh_start; b < a; b++) {
          cprim = I->Primitive + b;
          if(TriangleReverse(cprim))
            sprintf(buffer, "%d %d %d -1,\n", tri, tri + 2, tri + 1);
          else
            sprintf(buffer, "%d %d %d -1,\n", tri, tri + 1, tri + 2);
          UtilConcatVLA(&vla, &cc, buffer);
          tri += 3;
        }

        /* close mesh */
        UtilConcatVLA(&vla, &cc, " ] \n" " }\n" "}\n");
        mesh_obj = false;
      }

      switch (prim->type) {
      case cPrimSphere:
        sprintf(buffer,
                "Transform {\n"
                " translation %8.6f %8.6f %8.6f\n"
                " children Shape {\n"
                "  geometry Sphere { radius %8.6f }\n"
                "  appearance Appearance {\n"
                "   material Material { diffuseColor %6.4f %6.4f %6.4f \n"
                "                       specularColor 0.8 0.8 0.8 \n"
                "                       shininess 0.8 }\n"
                "  }\n"
                " }\n"
                "}\n",
                vert[0] - mid[0],
                vert[1] - mid[1],
                vert[2] - mid[2], prim->r1, prim->c1[0], prim->c1[1], prim->c1[2]);
        UtilConcatVLA(&vla, &cc, buffer);
        break;
      case cPrimCone:
        /* TO DO */
        break;
      case cPrimCylinder:
      case cPrimSausage:
        {
          float *d, vert2[3], axis[3], angle;
          OrthoLineType geometry;
          /* find the axis and angle that will rotate the y axis onto
           * the direction of the length of the cylinder
           */
          d = base->Normal + 3 * base->Vert2Normal[prim->vert];
          if((d[0] * d[0] + d[2] * d[2]) < 0.000001) {
            /* parallel with y */
            axis[0] = 1.0;
            axis[1] = 0.0;
            axis[2] = 0.0;
            if(d[1] > 0) {
              angle = 0.0;
            } else {
              angle = cPI;
            }
          } else {
            axis[0] = d[2];
            axis[1] = 0.0;
            axis[2] = -d[0];
            normalize3f(axis);
            angle = d[1];
            if(angle > 1.0)
              angle = 1.0;
            else if(angle < -1.0)
              angle = -1.0;
            angle = acos(angle);
          }
          /* vrml cylinders have origin in middle, not tip, that is why we
           * use prim->l1/2
           */
          scale3f(d, prim->l1 / 2, vert2);
          add3f(vert, vert2, vert2);
          if(prim->type == cPrimSausage) {
            OrthoLineType geom_add;
            sprintf(geometry,
                    "  Shape {\n"
                    "   geometry Cylinder {\n"
                    "    radius %8.6f\n"
                    "    height %8.6f\n"
                    "    bottom FALSE\n"
                    "    top    FALSE\n"
                    "   }\n"
                    "   appearance Appearance {\n"
                    "   material Material { diffuseColor %6.4f %6.4f %6.4f \n"
                    "                       specularColor 0.8 0.8 0.8 \n"
                    "                       shininess 0.8 }\n"
                    "   }\n",
                    prim->r1, prim->l1,
                    (prim->c1[0] + prim->c2[0]) / 2,
                    (prim->c1[1] + prim->c2[1]) / 2, (prim->c1[2] + prim->c2[2]) / 2);
            /* WLD: format string split to comply with ISO C89 standards */
            sprintf(geom_add,
                    "  }\n"
                    "  Transform {\n"
                    "   translation 0.0 %8.6f 0.0\n"
                    "   children Shape {\n"
                    "    geometry Sphere { radius %8.6f }\n"
                    "    appearance Appearance {\n"
                    "   material Material { diffuseColor %6.4f %6.4f %6.4f \n"
                    "                       specularColor 0.8 0.8 0.8 \n"
                    "                       shininess 0.8 }\n"
                    "    }\n"
                    "   }\n"
                    "  }\n", prim->l1 / 2, prim->r1, prim->c1[0], prim->c1[1], prim->c1[2]
              );
            strcat(geometry, geom_add);
            /* WLD: format string split to comply with ISO C89 standards */
            sprintf(geom_add,
                    "  Transform {\n"
                    "   translation 0.0 %8.6f 0.0\n"
                    "   children Shape {\n"
                    "    geometry Sphere { radius %8.6f }\n"
                    "    appearance Appearance {\n"
                    "   material Material { diffuseColor %6.4f %6.4f %6.4f \n"
                    "                       specularColor 0.8 0.8 0.8 \n"
                    "                       shininess 0.8 }\n"
                    "    }\n"
                    "   }\n"
                    "  }\n",
                    -prim->l1 / 2, prim->r1, prim->c2[0], prim->c2[1], prim->c2[2]);
            strcat(geometry, geom_add);
          } else {
            sprintf(geometry,
                    "  Shape {\n"
                    "   geometry Cylinder {\n"
                    "    radius %8.6f\n"
                    "    height %8.6f\n"
                    "   }\n"
                    "   appearance Appearance {\n"
                    "   material Material { diffuseColor %6.4f %6.4f %6.4f \n"
                    "                       specularColor 0.8 0.8 0.8 \n"
                    "                       shininess 0.8 }\n"
                    "   }\n"
                    "  }\n",
                    prim->r1, prim->l1,
                    (prim->c1[0] + prim->c2[0]) / 2,
                    (prim->c1[1] + prim->c2[1]) / 2, (prim->c1[2] + prim->c2[2]) / 2);
          }
          sprintf(buffer,
                  "Transform {\n"
                  " translation %8.6f %8.6f %8.6f\n"
                  " rotation %8.6f %8.6f %8.6f %8.6f\n"
                  " children [\n"
                  "%s"
                  " ]\n"
                  "}\n",
                  vert2[0] - mid[0],
                  vert2[1] - mid[1],
                  vert2[2] - mid[2], axis[0], axis[1], axis[2], angle, geometry);
          UtilConcatVLA(&vla, &cc, buffer);
        }
        break;
      case cPrimTriangle:
        /* output coords. connectivity and vertex colors handled above/below */
        sprintf(buffer,
                "%8.6f %8.6f %8.6f,\n"
                "%8.6f %8.6f %8.6f,\n"
                "%8.6f %8.6f %8.6f,\n",
                vert[0] - mid[0], vert[1] - mid[1], vert[2] - mid[2],
                vert[3] - mid[0], vert[4] - mid[1], vert[5] - mid[2],
                vert[6] - mid[0], vert[7] - mid[1], vert[8] - mid[2]);
        UtilConcatVLA(&vla, &cc, buffer);
        break;
      }
    }

    if(mesh_obj) {
      CBasis *base = I->Basis + 1;
      CPrimitive *cprim;
      int tri = 0;
      /* output connectivity */
      UtilConcatVLA(&vla, &cc, "   ]\n" "  }\n" "  coordIndex [\n");
      for(b = mesh_start; b < a; b++) {
        cprim = I->Primitive + b;
        if(TriangleReverse(cprim))
          sprintf(buffer, "%d %d %d -1,\n", tri, tri + 2, tri + 1);
        else
          sprintf(buffer, "%d %d %d -1,\n", tri, tri + 1, tri + 2);
        UtilConcatVLA(&vla, &cc, buffer);
        tri += 3;
      }

      /* output vertex colors */
      UtilConcatVLA(&vla, &cc,
                    "  ]\n" "  colorPerVertex TRUE\n" "  color Color {\n" "   color [\n");
      for(b = mesh_start; b < a; b++) {
        cprim = I->Primitive + b;
        sprintf(buffer,
                "%6.4f %6.4f %6.4f,\n"
                "%6.4f %6.4f %6.4f,\n"
                "%6.4f %6.4f %6.4f,\n",
                cprim->c1[0], cprim->c1[1], cprim->c1[2],
                cprim->c2[0], cprim->c2[1], cprim->c2[2],
                cprim->c3[0], cprim->c3[1], cprim->c3[2]);
        UtilConcatVLA(&vla, &cc, buffer);
      }

      /* output vertex normals */
      UtilConcatVLA(&vla, &cc,
                    "  ] } \n"
                    "  normalPerVertex TRUE\n" "  normal Normal {\n" "   vector [\n");
      for(b = mesh_start; b < a; b++) {
        cprim = I->Primitive + b;
        {
          float *norm = base->Normal + 3 * base->Vert2Normal[cprim->vert];
          sprintf(buffer, "%6.4f %6.4f %6.4f,\n" "%6.4f %6.4f %6.4f,\n" "%6.4f %6.4f %6.4f,\n", norm[3], norm[4], norm[5],      /* transformed cprim->n1 */
                  norm[6], norm[7], norm[8],    /* transformed cprim->n2 */
                  norm[9], norm[10], norm[11]); /* transformed cprim->n3 */
          UtilConcatVLA(&vla, &cc, buffer);
        }
      }
      UtilConcatVLA(&vla, &cc, "  ] }\n" "  normalIndex [ \n");
      tri = 0;
      for(b = mesh_start; b < a; b++) {
        cprim = I->Primitive + b;
        if(TriangleReverse(cprim))
          sprintf(buffer, "%d %d %d -1,\n", tri, tri + 2, tri + 1);
        else
          sprintf(buffer, "%d %d %d -1,\n", tri, tri + 1, tri + 2);
        UtilConcatVLA(&vla, &cc, buffer);
        tri += 3;
      }

      /* close mesh */
      UtilConcatVLA(&vla, &cc, " ] \n" " }\n" "}\n");
      mesh_obj = false;
    }
  }

  *vla_ptr = vla;
}


/* simple write-once/read-many hash for float-3/4 vectors  */

#define VECTOR_HASH_MASK 0xFFFF

typedef struct {
  float key[4];
  int value;
  int next;                     /* 1-based offsets, 0 = terminal */
} VectorHashElem;

typedef struct {
  int first[VECTOR_HASH_MASK + 1];
  VectorHashElem *elem;
  int size;
} VectorHash;

static void VectorHash_Free(VectorHash * I)
{
  if(I) {
    VLAFreeP(I->elem);
  }
  FreeP(I);
}

static VectorHash *VectorHash_New(void)
{
  VectorHash *I = pymol::calloc<VectorHash>(1);
  if(I) {
    I->elem = VLACalloc(VectorHashElem, 100);
    if(!I->elem) {
      VectorHash_Free(I);
      I = nullptr;
    }
  }
  return I;
}

static int VectorHash_GetOrSetKeyValue(VectorHash * I, float *key, float *alpha,
                                       int *value)
{
  unsigned int hash;
  /* returns non-zero if the entry is new */
  {
    unsigned int a, b, c;

    a = ((unsigned int *) key)[0];
    b = ((unsigned int *) key)[1];
    c = ((unsigned int *) key)[2];

    /* Robert Jenkin's 96 to 32 bit hash (public domain) 
       this is probably way overkill */

    a = a - b;
    a = a - c;
    a = a ^ (c >> 13);
    b = b - c;
    b = b - a;
    b = b ^ (a << 8);
    c = c - a;
    c = c - b;
    c = c ^ (b >> 13);
    a = a - b;
    a = a - c;
    a = a ^ (c >> 12);
    b = b - c;
    b = b - a;
    b = b ^ (a << 16);
    c = c - a;
    c = c - b;
    c = c ^ (b >> 5);
    a = a - b;
    a = a - c;
    a = a ^ (c >> 3);
    b = b - c;
    b = b - a;
    b = b ^ (a << 10);
    c = c - a;
    c = c - b;
    c = c ^ (b >> 15);

    /* mix in the fourth key (if present) */

    if(alpha) {
      unsigned int d = *((unsigned int *) alpha);
      c = c + d;
    }

    /* fold those 32 bits to 16 */

    c ^= (c >> 16);

    /* apply mask */

    hash = c & VECTOR_HASH_MASK;
  }
  {
    int offset = I->first[hash];
    while(offset) {
      VectorHashElem *elem = I->elem + offset;
      float *v = elem->key;
      if((v[0] == key[0]) && (v[1] == key[1]) && (v[2] == key[2])) {
        if((!alpha) || (*alpha == v[3])) {
          *value = elem->value; /* matched, so return value */
          return 0;             /* key/value exists */
        }
      }
      offset = elem->next;
    }
    /* not matched -- add new key/value  */
    if(VLACheck(I->elem, VectorHashElem, ++(I->size))) {
      VectorHashElem *elem = I->elem + I->size;
      elem->next = I->first[hash];
      I->first[hash] = I->size;
      copy3f(key, elem->key);
      if(alpha)
        elem->key[3] = *alpha;
      elem->value = *value;
      return 1;                 /* inform caller */
    } else {
      I->size--;
      return -1;
    }
  }
}

static void unique_vector_add(VectorHash * vh, float *vector,
                              float *vector_array, int *vector_count,
                              int *index_array, int *index_count)
{
  int index = *vector_count;
  switch (VectorHash_GetOrSetKeyValue(vh, vector, nullptr, &index)) {
  case 1:
    {
      float *vector_slot = vector_array + 3 * (*vector_count);
      copy3f(vector, vector_slot);
      (*vector_count)++;
    }
    /* INTENTIONAL omission of break */
  case 0:
    index_array[(*index_count)++] = index;
    break;
  }
}

static void unique_color_add(VectorHash * vh, float *vector,
                             float *vector_array, int *vector_count,
                             int *index_array, int *index_count, float alpha)
{
  int index = *vector_count;

  switch (VectorHash_GetOrSetKeyValue(vh, vector, &alpha, &index)) {
  case 1:
    {
      float *vector_slot = vector_array + 4 * (*vector_count);  /* NOTE 4x spacing */
      copy3f(vector, vector_slot);
      vector_slot[3] = alpha;
      (*vector_count)++;
    }
    /* INTENTIONAL omission of break */
  case 0:
    index_array[(*index_count)++] = index;
    break;
  }
}

#define noIDTF_COLOR

typedef struct {
  int face_count;
  int position_count;
  int normal_count;
  int *face_position_list;
  int *face_normal_list;
  int *face_shading_list;
  float *model_position_list;
  float *model_normal_list;
  VectorHash *position_hash;
  VectorHash *normal_hash;
  int color_count;
  int *face_color_list;
  float *model_diffuse_color_list;
  VectorHash *color_hash;
} IdtfResourceMesh;

typedef struct {
  float *color_list;
  int color_count;
  VectorHash *color_hash;
} IdtfMaterial;

static ov_size idtf_dump_file_header(char **vla, ov_size cnt)
{
  UtilConcatVLA(vla, &cnt, "FILE_FORMAT \"IDTF\"\nFORMAT_VERSION 100\n\n");

  UtilConcatVLA(vla, &cnt, "NODE \"VIEW\" {\n");
  UtilConcatVLA(vla, &cnt, "\tNODE_NAME \"DefaultView\"\n");
  UtilConcatVLA(vla, &cnt, "\tPARENT_LIST {\n");
  UtilConcatVLA(vla, &cnt, "\t\tPARENT_COUNT 1\n");
  UtilConcatVLA(vla, &cnt, "\t\tPARENT 0 {\n\t\t\tPARENT_NAME \"<NULL>\"\n");
  UtilConcatVLA(vla, &cnt, "\t\t\tPARENT_TM {\n");
  UtilConcatVLA(vla, &cnt, "\t\t\t\t1.000000 0.000000 0.000000 0.0\n");
  UtilConcatVLA(vla, &cnt, "\t\t\t\t0.000000 1.000000 0.000000 0.0\n");
  UtilConcatVLA(vla, &cnt, "\t\t\t\t0.000000 0.000000 1.000000 0.0\n");
  UtilConcatVLA(vla, &cnt, "\t\t\t\t0.000000 0.000000 0.000000 1.0\n");
  UtilConcatVLA(vla, &cnt, "\t\t\t}\n");
  UtilConcatVLA(vla, &cnt, "\t\t}\n");
  UtilConcatVLA(vla, &cnt, "\t}\n");
  UtilConcatVLA(vla, &cnt, "\tRESOURCE_NAME \"SceneViewResource\"\n");
  UtilConcatVLA(vla, &cnt, "\tVIEW_DATA {\n");
  UtilConcatVLA(vla, &cnt, "\t\tVIEW_TYPE \"PERSPECTIVE\"\n");
  UtilConcatVLA(vla, &cnt, "\t\tVIEW_PROJECTION 34.515877\n");
  UtilConcatVLA(vla, &cnt, "\t}\n");
  UtilConcatVLA(vla, &cnt, "}\n\n");

  UtilConcatVLA(vla, &cnt, "NODE \"LIGHT\"\n");
  UtilConcatVLA(vla, &cnt, "{\n");
  UtilConcatVLA(vla, &cnt, "\tNODE_NAME \"Omni01\"\n");
  UtilConcatVLA(vla, &cnt, "\tPARENT_LIST {\n");
  UtilConcatVLA(vla, &cnt, "\t\tPARENT_COUNT 1\n");
  UtilConcatVLA(vla, &cnt, "\t\tPARENT 0 {\n");
  UtilConcatVLA(vla, &cnt, "\t\t\tPARENT_NAME \"<NULL>\"\n");
  UtilConcatVLA(vla, &cnt, "\t\t\tPARENT_TM {\n");
  UtilConcatVLA(vla, &cnt, "\t\t\t\t1.000000 0.000000 0.000000 0.000000\n");
  UtilConcatVLA(vla, &cnt, "\t\t\t\t0.000000 1.000000 0.000000 0.000000\n");
  UtilConcatVLA(vla, &cnt, "\t\t\t\t0.000000 0.000000 1.000000 0.000000\n");
  UtilConcatVLA(vla, &cnt, "\t\t\t\t50.000000 -50.00000 50.000000 1.000000\n");
  UtilConcatVLA(vla, &cnt, "\t\t\t}\n");
  UtilConcatVLA(vla, &cnt, "\t\t}\n");
  UtilConcatVLA(vla, &cnt, "\t}\n");
  UtilConcatVLA(vla, &cnt, "\tRESOURCE_NAME \"DefaultPointLight\"\n");
  UtilConcatVLA(vla, &cnt, "}\n\n");
  return cnt;
}

static ov_size idtf_dump_model_nodes(char **vla, ov_size cnt,
                                     IdtfResourceMesh * mesh_vla, int n_mesh)
{
  int a;
  IdtfResourceMesh *mesh = mesh_vla;
  for(a = 0; a < n_mesh; a++) {
    OrthoLineType buffer;

    UtilConcatVLA(vla, &cnt, "NODE \"MODEL\" {\n");

    sprintf(buffer, "\tNODE_NAME \"Mesh%d\"\n", a);
    UtilConcatVLA(vla, &cnt, buffer);

    UtilConcatVLA(vla, &cnt, "\tPARENT_LIST {\n");
    UtilConcatVLA(vla, &cnt, "\t\tPARENT_COUNT 1\n");
    UtilConcatVLA(vla, &cnt, "\t\tPARENT 0 {\n");
    UtilConcatVLA(vla, &cnt, "\t\t\tPARENT_NAME \"<NULL>\"\n");
    UtilConcatVLA(vla, &cnt, "\t\t\tPARENT_TM {\n");
    UtilConcatVLA(vla, &cnt, "\t\t\t1.000000 0.000000 0.000000 0.0\n");
    UtilConcatVLA(vla, &cnt, "\t\t\t0.000000 1.000000 0.000000 0.0\n");
    UtilConcatVLA(vla, &cnt, "\t\t\t0.000000 0.000000 1.000000 0.0\n");
    UtilConcatVLA(vla, &cnt, "\t\t\t0.000000 0.000000 0.000000 1.0\n");
    UtilConcatVLA(vla, &cnt, "\t\t\t}\n");
    UtilConcatVLA(vla, &cnt, "\t\t}\n");
    UtilConcatVLA(vla, &cnt, "\t}\n");

    sprintf(buffer, "\tRESOURCE_NAME \"Mesh%d\"\n", a);
    UtilConcatVLA(vla, &cnt, buffer);

    UtilConcatVLA(vla, &cnt, "}\n\n");

    mesh++;
  }
  return cnt;
}

static ov_size idtf_dump_resource_header(char **vla, ov_size cnt)
{

  UtilConcatVLA(vla, &cnt, "RESOURCE_LIST \"VIEW\" {\n");
  UtilConcatVLA(vla, &cnt, "\tRESOURCE_COUNT 1\n");
  UtilConcatVLA(vla, &cnt, "\tRESOURCE 0 {\n");
  UtilConcatVLA(vla, &cnt, "\t\tRESOURCE_NAME \"SceneViewResource\"\n");
  UtilConcatVLA(vla, &cnt, "\t\tVIEW_PASS_COUNT 1\n");
  UtilConcatVLA(vla, &cnt, "\t\tVIEW_ROOT_NODE_LIST {\n");
  UtilConcatVLA(vla, &cnt, "\t\t\tROOT_NODE 0 {\n");
  UtilConcatVLA(vla, &cnt, "\t\t\t\tROOT_NODE_NAME \"<NULL>\"\n");
  UtilConcatVLA(vla, &cnt, "\t\t\t}\n");
  UtilConcatVLA(vla, &cnt, "\t\t}\n");
  UtilConcatVLA(vla, &cnt, "\t}\n");
  UtilConcatVLA(vla, &cnt, "}\n\n");

  UtilConcatVLA(vla, &cnt, "RESOURCE_LIST \"LIGHT\" {\n");
  UtilConcatVLA(vla, &cnt, "\tRESOURCE_COUNT 1\n");
  UtilConcatVLA(vla, &cnt, "\tRESOURCE 0 {\n");
  UtilConcatVLA(vla, &cnt, "\t\tRESOURCE_NAME \"DefaultPointLight\"\n");
  UtilConcatVLA(vla, &cnt, "\t\tLIGHT_TYPE \"POINT\"\n");
  UtilConcatVLA(vla, &cnt, "\t\tLIGHT_COLOR 1.000000 1.000000 1.000000\n");
  UtilConcatVLA(vla, &cnt, "\t\tLIGHT_ATTENUATION 1.000000 0.000000 0.000000\n");
  UtilConcatVLA(vla, &cnt, "\t\tLIGHT_INTENSITY 1.000000\n");
  UtilConcatVLA(vla, &cnt, "\t}\n");
  UtilConcatVLA(vla, &cnt, "}\n\n");

  return cnt;
}

static ov_size idtf_dump_resources(char **vla, ov_size cnt,
                                   IdtfResourceMesh * mesh_vla, int n_mesh,
                                   IdtfMaterial * material)
{
  {
    OrthoLineType buffer;
    int n_color = material->color_count;

    UtilConcatVLA(vla, &cnt, "RESOURCE_LIST \"SHADER\" {\n");

    sprintf(buffer, "\tRESOURCE_COUNT %d\n", n_color);
    UtilConcatVLA(vla, &cnt, buffer);

    {
      int c;
      for(c = 0; c < n_color; c++) {

        sprintf(buffer, "\tRESOURCE %d {\n", c);
        UtilConcatVLA(vla, &cnt, buffer);

        sprintf(buffer, "\t\tRESOURCE_NAME \"Shader%06d\"\n", c);
        UtilConcatVLA(vla, &cnt, buffer);

        sprintf(buffer, "\t\tSHADER_MATERIAL_NAME \"Material%06d\"\n", c);
        UtilConcatVLA(vla, &cnt, buffer);

        UtilConcatVLA(vla, &cnt, "\t\tSHADER_ACTIVE_TEXTURE_COUNT 0\n");
        UtilConcatVLA(vla, &cnt, "\t}\n");
      }
    }
    UtilConcatVLA(vla, &cnt, "}\n\n");
  }

  {
    OrthoLineType buffer;
    int n_color = material->color_count;

    UtilConcatVLA(vla, &cnt, "RESOURCE_LIST \"MATERIAL\" {\n");

    sprintf(buffer, "\tRESOURCE_COUNT %d\n", n_color);
    UtilConcatVLA(vla, &cnt, buffer);

    {
      int c;
      float *fp = material->color_list;

      for(c = 0; c < n_color; c++) {
        sprintf(buffer, "\tRESOURCE %d {\n", c);
        UtilConcatVLA(vla, &cnt, buffer);
        sprintf(buffer, "\t\tRESOURCE_NAME \"Material%06d\"\n", c);
        UtilConcatVLA(vla, &cnt, buffer);

        sprintf(buffer, "\t\tMATERIAL_AMBIENT %0.6f %0.6f %0.6f\n",
                fp[0] * 0, fp[1] * 0, fp[2] * 0);
        UtilConcatVLA(vla, &cnt, buffer);

        sprintf(buffer, "\t\tMATERIAL_DIFFUSE %0.6f %0.6f %0.6f\n", fp[0], fp[1], fp[2]);
        UtilConcatVLA(vla, &cnt, buffer);

        UtilConcatVLA(vla, &cnt, "\t\tMATERIAL_SPECULAR 0.750000 0.750000 0.750000\n");

        sprintf(buffer, "\t\tMATERIAL_EMISSIVE %0.6f %0.6f %0.6f\n",
                fp[0] * 0.13, fp[1] * 0.13, fp[2] * 0.13);
        UtilConcatVLA(vla, &cnt, buffer);
        UtilConcatVLA(vla, &cnt, "\t\tMATERIAL_REFLECTIVITY 0.40000\n");

        sprintf(buffer, "\t\tMATERIAL_OPACITY %0.6f\n", fp[3]);
        UtilConcatVLA(vla, &cnt, buffer);

        UtilConcatVLA(vla, &cnt, "\t}\n");

        fp += 4;
      }
    }
    UtilConcatVLA(vla, &cnt, "}\n\n");
  }

  {

    OrthoLineType buffer;
    UtilConcatVLA(vla, &cnt, "RESOURCE_LIST \"MODEL\" {\n");

    sprintf(buffer, "\tRESOURCE_COUNT %d\n", n_mesh);
    UtilConcatVLA(vla, &cnt, buffer);

    {
      int a;
      IdtfResourceMesh *mesh = mesh_vla;
      for(a = 0; a < n_mesh; a++) {

        sprintf(buffer, "\tRESOURCE %d {\n", a);
        UtilConcatVLA(vla, &cnt, buffer);

        sprintf(buffer, "\t\tRESOURCE_NAME \"Mesh%d\"\n", a);
        UtilConcatVLA(vla, &cnt, buffer);

        UtilConcatVLA(vla, &cnt, "\t\tMODEL_TYPE \"MESH\"\n");
        UtilConcatVLA(vla, &cnt, "\t\tMESH {\n");

        sprintf(buffer, "\t\t\tFACE_COUNT %d\n", mesh->face_count);
        UtilConcatVLA(vla, &cnt, buffer);

        sprintf(buffer, "\t\t\tMODEL_POSITION_COUNT %d\n", mesh->position_count);
        UtilConcatVLA(vla, &cnt, buffer);

        sprintf(buffer, "\t\t\tMODEL_NORMAL_COUNT %d\n", mesh->normal_count);
        UtilConcatVLA(vla, &cnt, buffer);
#ifdef IDTF_COLOR
        sprintf(buffer, "\t\t\tMODEL_DIFFUSE_COLOR_COUNT %d\n", mesh->color_count);
        UtilConcatVLA(vla, &cnt, buffer);

        sprintf(buffer, "\t\t\tMODEL_SPECULAR_COLOR_COUNT %d\n", mesh->color_count);
        UtilConcatVLA(vla, &cnt, buffer);
#else
        UtilConcatVLA(vla, &cnt, "\t\t\tMODEL_DIFFUSE_COLOR_COUNT 0\n");
        UtilConcatVLA(vla, &cnt, "\t\t\tMODEL_SPECULAR_COLOR_COUNT 0\n");
#endif
        UtilConcatVLA(vla, &cnt, "\t\t\tMODEL_TEXTURE_COORD_COUNT 0\n");
        UtilConcatVLA(vla, &cnt, "\t\t\tMODEL_BONE_COUNT 0\n");

        {
          int n_color = material->color_count;

          sprintf(buffer, "\t\t\tMODEL_SHADING_COUNT %d\n", n_color);
          UtilConcatVLA(vla, &cnt, buffer);

          UtilConcatVLA(vla, &cnt, "\t\t\tMODEL_SHADING_DESCRIPTION_LIST {\n");

          {
            int c;
            for(c = 0; c < n_color; c++) {

              sprintf(buffer, "\t\t\t\tSHADING_DESCRIPTION %d {\n", c);
              UtilConcatVLA(vla, &cnt, buffer);

              UtilConcatVLA(vla, &cnt, "\t\t\t\tTEXTURE_LAYER_COUNT 0\n");

              sprintf(buffer, "\t\t\t\tSHADER_ID %d\n", c + 1);
              UtilConcatVLA(vla, &cnt, buffer);

              UtilConcatVLA(vla, &cnt, "\t\t\t\t}\n");
            }
          }
          UtilConcatVLA(vla, &cnt, "\t\t\t}\n");
        }

        {
          int b;
          int *ip = mesh->face_position_list;
          UtilConcatVLA(vla, &cnt, "\t\t\tMESH_FACE_POSITION_LIST {\n");

          for(b = 0; b < mesh->face_count; b++) {
            sprintf(buffer, "\t\t\t%d %d %d\n", ip[0], ip[1], ip[2]);
            UtilConcatVLA(vla, &cnt, buffer);
            ip += 3;
          }
          UtilConcatVLA(vla, &cnt, "\t\t\t}\n");
        }

        {
          int b;
          int *ip = mesh->face_normal_list;
          UtilConcatVLA(vla, &cnt, "\t\t\tMESH_FACE_NORMAL_LIST {\n");

          for(b = 0; b < mesh->face_count; b++) {
            sprintf(buffer, "\t\t\t%d %d %d\n", ip[0], ip[1], ip[2]);
            UtilConcatVLA(vla, &cnt, buffer);
            ip += 3;
          }
          UtilConcatVLA(vla, &cnt, "\t\t\t}\n");
        }

        {
          int b;
          int *ip = mesh->face_shading_list;
          UtilConcatVLA(vla, &cnt, "\t\t\tMESH_FACE_SHADING_LIST {\n");

          for(b = 0; b < mesh->face_count; b++) {
            sprintf(buffer, "\t\t\t%d\n", ip[0]);
            UtilConcatVLA(vla, &cnt, buffer);
            ip++;
          }
          UtilConcatVLA(vla, &cnt, "\t\t\t}\n");
        }

#ifdef IDTF_COLOR
        {
          int b;
          int *ip = mesh->face_color_list;
          UtilConcatVLA(vla, &cnt, "\t\t\tMESH_FACE_DIFFUSE_COLOR_LIST {\n");

          for(b = 0; b < mesh->face_count; b++) {
            sprintf(buffer, "\t\t\t%d %d %d\n", ip[0], ip[1], ip[2]);
            UtilConcatVLA(vla, &cnt, buffer);
            ip += 3;
          }
          UtilConcatVLA(vla, &cnt, "\t\t\t}\n");
        }
        {
          int b;
          int *ip = mesh->face_color_list;
          UtilConcatVLA(vla, &cnt, "\t\t\tMESH_FACE_SPECULAR_COLOR_LIST {\n");

          for(b = 0; b < mesh->face_count; b++) {
            sprintf(buffer, "\t\t\t%d %d %d\n", ip[0], ip[1], ip[2]);
            UtilConcatVLA(vla, &cnt, buffer);
            ip += 3;
          }
          UtilConcatVLA(vla, &cnt, "\t\t\t}\n");
        }
#endif

        {
          int b;
          float *fp = mesh->model_position_list;
          UtilConcatVLA(vla, &cnt, "\t\t\tMODEL_POSITION_LIST {\n");

          for(b = 0; b < mesh->position_count; b++) {
            sprintf(buffer, "\t\t\t\t%1.6f %1.6f %1.6f\n", fp[0], fp[1], fp[2]);
            UtilConcatVLA(vla, &cnt, buffer);
            fp += 3;
          }

          UtilConcatVLA(vla, &cnt, "\t\t\t}\n");
        }

        {
          int b;
          float *fp = mesh->model_normal_list;
          UtilConcatVLA(vla, &cnt, "\t\t\tMODEL_NORMAL_LIST {\n");

          for(b = 0; b < mesh->normal_count; b++) {
            sprintf(buffer, "\t\t\t\t%1.6f %1.6f %1.6f\n", fp[0], fp[1], fp[2]);
            UtilConcatVLA(vla, &cnt, buffer);
            fp += 3;
          }

          UtilConcatVLA(vla, &cnt, "\t\t\t}\n");
        }
#ifdef IDTF_COLOR
        {
          int b;
          float *fp = mesh->model_diffuse_color_list;
          UtilConcatVLA(vla, &cnt, "\t\t\tMODEL_DIFFUSE_COLOR_LIST {\n");

          for(b = 0; b < mesh->color_count; b++) {
            sprintf(buffer, "\t\t\t\t%1.6f %1.6f %1.6f %1.6f\n", fp[0], fp[1], fp[2],
                    fp[3]);
            UtilConcatVLA(vla, &cnt, buffer);
            fp += 4;
          }

          UtilConcatVLA(vla, &cnt, "\t\t\t}\n");
        }
        {
          int b;
          float *fp = mesh->model_diffuse_color_list;
          UtilConcatVLA(vla, &cnt, "\t\t\tMODEL_SPECULAR_COLOR_LIST {\n");

          for(b = 0; b < mesh->color_count; b++) {
            sprintf(buffer, "\t\t\t\t%1.6f %1.6f %1.6f %1.6f\n", fp[0], fp[1], fp[2],
                    fp[3]);
            UtilConcatVLA(vla, &cnt, buffer);
            fp += 4;
          }

          UtilConcatVLA(vla, &cnt, "\t\t\t}\n");
        }
#endif

        UtilConcatVLA(vla, &cnt, "\t\t}\n");
        UtilConcatVLA(vla, &cnt, "\t}\n");

        mesh++;
      }
    }
    UtilConcatVLA(vla, &cnt, "}\n\n");
  }
  return cnt;
}


/*========================================================================*/
void RayRenderIDTF(CRay * I, char **node_vla, char **rsrc_vla)
{
  int identity = (SettingGetGlobal_i(I->G, cSetting_geometry_export_mode) == 1);

  RayExpandPrimitives(I);
  RayTransformFirst(I, 0, identity);

  {
    CBasis *base = I->Basis + 1;
    CPrimitive *prim = I->Primitive;
    int mesh_cnt = 0;
    IdtfResourceMesh *mesh_vla = VLACalloc(IdtfResourceMesh, 1);
    if(mesh_vla) {
      IdtfResourceMesh *mesh = nullptr;
      int a;

      for(a = 0; a < I->NPrimitive; a++) {

        switch (prim->type) {
        case cPrimTriangle:
        case cPrimSphere:
          if(!mesh) {
            /* create a new triangle mesh */
            if(VLACheck(mesh_vla, IdtfResourceMesh, mesh_cnt)) {
              mesh = mesh_vla + mesh_cnt;
              if((mesh->face_position_list = VLACalloc(int, 3)) && (mesh->face_normal_list = VLACalloc(int, 3)) && (mesh->face_shading_list = VLACalloc(int, 1)) &&     /* defaults to zero */
                  
                 (mesh->model_position_list = VLAlloc(float, 3)) &&
                 (mesh->face_color_list = VLACalloc(int, 3)) &&
                 (mesh->model_diffuse_color_list = VLAlloc(float, 4)) &&
                 ((mesh->color_hash = VectorHash_New())) &&
                 (mesh->model_normal_list = VLAlloc(float, 3)) &&
                 ((mesh->position_hash = VectorHash_New())) &&
                 ((mesh->normal_hash = VectorHash_New()))
                ) {
                mesh_cnt++;
              } else {
                mesh = nullptr;
              }
            }
          }
          break;
        default:               /* close/terminate mesh */
          if(mesh) {
            mesh = nullptr;
          }
          break;
        }

        switch (prim->type) {
        case cPrimTriangle:
          if(mesh) {
            if(VLACheck(mesh->face_position_list, int, mesh->face_count * 3 + 2) &&
               VLACheck(mesh->face_normal_list, int, mesh->face_count * 3 + 2) &&
               VLACheck(mesh->face_shading_list, int, mesh->face_count) &&
               VLACheck(mesh->model_position_list, float, (mesh->position_count + 3) * 3)
               && VLACheck(mesh->face_color_list, int, mesh->face_count * 3 + 2)
               && VLACheck(mesh->model_diffuse_color_list, float,
                           (mesh->color_count + 3) * 4)
               && VLACheck(mesh->model_normal_list, float, (mesh->normal_count + 3) * 3)
              ) {

              float *vert = base->Vertex + 3 * (prim->vert);
              float *norm = base->Normal + 3 * base->Vert2Normal[prim->vert] + 3;
              int reverse = TriangleReverse(prim);
              int face_position_count = mesh->face_count * 3;
              int face_normal_count = face_position_count;
              int face_color_count = face_position_count;

              unique_vector_add(mesh->position_hash, vert,
                                mesh->model_position_list, &mesh->position_count,
                                mesh->face_position_list, &face_position_count);
              unique_vector_add(mesh->normal_hash, norm,
                                mesh->model_normal_list, &mesh->normal_count,
                                mesh->face_normal_list, &face_normal_count);
              unique_color_add(mesh->normal_hash, prim->c1,
                               mesh->model_diffuse_color_list, &mesh->color_count,
                               mesh->face_color_list, &face_color_count,
                               1.0F - prim->trans);
              if(reverse) {
                vert += 6;
                norm += 6;
                unique_vector_add(mesh->position_hash, vert,
                                  mesh->model_position_list, &mesh->position_count,
                                  mesh->face_position_list, &face_position_count);
                unique_vector_add(mesh->normal_hash, norm,
                                  mesh->model_normal_list, &mesh->normal_count,
                                  mesh->face_normal_list, &face_normal_count);
                unique_color_add(mesh->normal_hash, prim->c3,
                                 mesh->model_diffuse_color_list, &mesh->color_count,
                                 mesh->face_color_list, &face_color_count,
                                 1.0F - prim->trans);
                vert -= 3;
                norm -= 3;
                unique_vector_add(mesh->position_hash, vert,
                                  mesh->model_position_list, &mesh->position_count,
                                  mesh->face_position_list, &face_position_count);
                unique_vector_add(mesh->normal_hash, norm,
                                  mesh->model_normal_list, &mesh->normal_count,
                                  mesh->face_normal_list, &face_normal_count);
                unique_color_add(mesh->normal_hash, prim->c2,
                                 mesh->model_diffuse_color_list, &mesh->color_count,
                                 mesh->face_color_list, &face_color_count,
                                 1.0F - prim->trans);
              } else {
                vert += 3;
                norm += 3;
                unique_vector_add(mesh->position_hash, vert,
                                  mesh->model_position_list, &mesh->position_count,
                                  mesh->face_position_list, &face_position_count);
                unique_vector_add(mesh->normal_hash, norm,
                                  mesh->model_normal_list, &mesh->normal_count,
                                  mesh->face_normal_list, &face_normal_count);
                unique_color_add(mesh->normal_hash, prim->c2,
                                 mesh->model_diffuse_color_list, &mesh->color_count,
                                 mesh->face_color_list, &face_color_count,
                                 1.0F - prim->trans);
                vert += 3;
                norm += 3;
                unique_vector_add(mesh->position_hash, vert,
                                  mesh->model_position_list, &mesh->position_count,
                                  mesh->face_position_list, &face_position_count);
                unique_vector_add(mesh->normal_hash, norm,
                                  mesh->model_normal_list, &mesh->normal_count,
                                  mesh->face_normal_list, &face_normal_count);
                unique_color_add(mesh->normal_hash, prim->c3,
                                 mesh->model_diffuse_color_list, &mesh->color_count,
                                 mesh->face_color_list, &face_color_count,
                                 1.0F - prim->trans);
              }
              mesh->face_count++;
            }
          }
          break;
        case cPrimSphere:

          break;
        case cPrimCone:
          break;
        case cPrimCylinder:
        case cPrimSausage:
          break;
        }

        /* looping through each primitive */
        prim++;
      }

      /* now we need to consolidate materials for each color
         combination (creating averages as necessary) and update
         mesh->face_shading_list appropriately for each face */

      {
        IdtfMaterial *material = pymol::calloc<IdtfMaterial>(1);
        if(material &&
           (material->color_list = VLAlloc(float, 4)) &&
           (material->color_hash = VectorHash_New())) {

          const float one_third = (1 / 3.0F);
          int a;
          IdtfResourceMesh *mesh = mesh_vla;

          for(a = 0; a < mesh_cnt; a++) {
            int *ip = mesh->face_color_list;
            int shading_count = 0;
            int b;
            for(b = 0; b < mesh->face_count; b++) {
              float *fp0 = mesh->model_diffuse_color_list + 4 * ip[0];
              float *fp1 = mesh->model_diffuse_color_list + 4 * ip[1];
              float *fp2 = mesh->model_diffuse_color_list + 4 * ip[2];
              float avg_rgba[4];

              avg_rgba[0] = (fp0[0] + fp1[0] + fp2[0]) * one_third;
              avg_rgba[1] = (fp0[1] + fp1[1] + fp2[1]) * one_third;
              avg_rgba[2] = (fp0[2] + fp1[2] + fp2[2]) * one_third;
              avg_rgba[3] = (fp0[3] + fp1[3] + fp2[3]) * one_third;

              if(VLACheck(material->color_list, float, material->color_count * 4 + 3)) {
                unique_color_add(material->color_hash, avg_rgba,
                                 material->color_list, &material->color_count,
                                 mesh->face_shading_list, &shading_count, avg_rgba[3]);
              }
              ip += 3;
            }
            mesh++;
          }

          {
            int cnt = 0;
            cnt = idtf_dump_file_header(node_vla, cnt);
            cnt = idtf_dump_model_nodes(node_vla, cnt, mesh_vla, mesh_cnt);
            VLASize((*node_vla), char, cnt + 1);
          }
          {
            int cnt = 0;
            cnt = idtf_dump_resource_header(rsrc_vla, cnt);
            cnt = idtf_dump_resources(rsrc_vla, cnt, mesh_vla, mesh_cnt, material);
            VLASize((*rsrc_vla), char, cnt + 1);
          }

          VLAFreeP(material->color_list);
          VectorHash_Free(material->color_hash);
        }
        FreeP(material);
      }

      {
        /* refactor into a delete method */
        IdtfResourceMesh *mesh = mesh_vla;
        int i;
        for(i = 0; i < mesh_cnt; i++) {
          VLAFreeP(mesh->face_position_list);
          VLAFreeP(mesh->face_normal_list);
          VLAFreeP(mesh->face_shading_list);
          VLAFreeP(mesh->face_color_list);
          VLAFreeP(mesh->model_position_list);
          VLAFreeP(mesh->model_normal_list);
          VLAFreeP(mesh->model_diffuse_color_list);
          VectorHash_Free(mesh->color_hash);
          VectorHash_Free(mesh->position_hash);
          VectorHash_Free(mesh->normal_hash);
          mesh++;
        }
        VLAFreeP(mesh_vla);
      }
    }
  }
}


/*========================================================================*/
void RayRenderObjMtl(CRay * I, int width, int height, char **objVLA_ptr,
                     char **mtlVLA_ptr, float front, float back, float fov,
                     float angle, float z_corr)
{
  char *objVLA = *objVLA_ptr;
  char *mtlVLA = *mtlVLA_ptr;
  int identity = (SettingGetGlobal_i(I->G, cSetting_geometry_export_mode) == 1);

  ov_size oc = 0;               /* obj character count */
  /*  int mc = 0; *//* mtl character count */

  OrthoLineType buffer;

  RayExpandPrimitives(I);
  RayTransformFirst(I, 0, identity);

  {
    int a;
    CPrimitive *prim;
    float *vert, *norm;
    int vc = 0;
    int nc = 0;
    CBasis *base = I->Basis + 1;

    for(a = 0; a < I->NPrimitive; a++) {
      prim = I->Primitive + a;
      vert = base->Vertex + 3 * (prim->vert);
      norm = base->Normal + 3 * base->Vert2Normal[prim->vert] + 3;
      switch (prim->type) {
      case cPrimTriangle:
        sprintf(buffer, "v %8.6f %8.6f %8.6f\n", vert[0], vert[1], vert[2] - z_corr);
        UtilConcatVLA(&objVLA, &oc, buffer);
        sprintf(buffer, "v %8.6f %8.6f %8.6f\n", vert[3], vert[4], vert[5] - z_corr);
        UtilConcatVLA(&objVLA, &oc, buffer);
        sprintf(buffer, "v %8.6f %8.6f %8.6f\n", vert[6], vert[7], vert[8] - z_corr);
        UtilConcatVLA(&objVLA, &oc, buffer);
        sprintf(buffer, "vn %8.6f %8.6f %8.6f\n", norm[0], norm[1], norm[2]);
        UtilConcatVLA(&objVLA, &oc, buffer);
        sprintf(buffer, "vn %8.6f %8.6f %8.6f\n", norm[3], norm[4], norm[5]);
        UtilConcatVLA(&objVLA, &oc, buffer);
        sprintf(buffer, "vn %8.6f %8.6f %8.6f\n", norm[6], norm[7], norm[8]);
        UtilConcatVLA(&objVLA, &oc, buffer);
        if(TriangleReverse(prim)) {
          sprintf(buffer, "f %d//%d %d//%d %d//%d\n",
                  vc + 1, nc + 1, vc + 3, nc + 3, vc + 2, nc + 2);
        } else {
          sprintf(buffer, "f %d//%d %d//%d %d//%d\n",
                  vc + 1, nc + 1, vc + 2, nc + 2, vc + 3, nc + 3);
        }
        UtilConcatVLA(&objVLA, &oc, buffer);
        nc += 3;
        vc += 3;

        /*
           prim->c1[0],prim->c1[1],prim->c1[2])
           prim->c2[0],prim->c2[1],prim->c2[2],
           prim->c3[0],prim->c3[1],prim->c3[2]
           UtilConcatVLA(&vla,&oc,buffer);
           UtilConcatVLA(&vla,&oc,buffer);
         */

        break;
      case cPrimSphere:
        /*      sprintf(buffer,"v %8.6f %8.6f %8.6f\np %d\n",
           vert[0],vert[1],vert[2]-z_corr,vc+1);
           UtilConcatVLA(&objVLA,&oc,buffer); */

        sprintf(buffer, "v %8.6f %8.6f %8.6f\n", vert[0], vert[1], vert[2] - z_corr);
        UtilConcatVLA(&objVLA, &oc, buffer);
        sprintf(buffer, "v %8.6f %8.6f %8.6f\n", vert[0], vert[1], vert[2] - z_corr);
        UtilConcatVLA(&objVLA, &oc, buffer);
        sprintf(buffer, "v %8.6f %8.6f %8.6f\n", vert[0], vert[1], vert[2] - z_corr);
        UtilConcatVLA(&objVLA, &oc, buffer);
        sprintf(buffer, "f %d %d %d\n", vc + 1, vc + 2, vc + 3);
        UtilConcatVLA(&objVLA, &oc, buffer);
        vc += 3;
        break;
      }
    }
  }

  *objVLA_ptr = objVLA;
  *mtlVLA_ptr = mtlVLA;
}


/*========================================================================*/
void RayRenderPOV(CRay * I, int width, int height, char **headerVLA_ptr,
                  char **charVLA_ptr, float front, float back, float fov,
                  float angle, int antialias)
{
  float fog;
  const float *bkrd;
  float fog_start = 0.0F;
  float gamma;
  float *d;
  CBasis *base;
  CPrimitive *prim;
  OrthoLineType buffer;
  float *vert, *norm;
  float vert2[3];
  ov_size cc, hc;
  int a;
  int smooth_color_triangle;
  int mesh_obj = false;
  char *charVLA, *headerVLA;
  char transmit[64];
  float light[3];
  const float *lightv;
  float spec_power = SettingGetGlobal_f(I->G, cSetting_spec_power);
  int identity = (SettingGetGlobal_i(I->G, cSetting_geometry_export_mode) == 1);

  if(spec_power < 0.0F) {
    spec_power = SettingGetGlobal_f(I->G, cSetting_shininess);
  }
  spec_power /= 4.0F;

  charVLA = *charVLA_ptr;
  headerVLA = *headerVLA_ptr;
  smooth_color_triangle = SettingGetGlobal_b(I->G, cSetting_smooth_color_triangle);
  PRINTFB(I->G, FB_Ray, FB_Blather)
    " RayRenderPOV: w %d h %d f %8.3f b %8.3f\n", width, height, front, back ENDFB(I->G);
  if(Feedback(I->G, FB_Ray, FB_Blather)) {
    dump3f(I->Volume, " RayRenderPOV: vol");
    dump3f(I->Volume + 3, " RayRenderPOV: vol");
  }
  cc = 0;
  hc = 0;
  gamma = SettingGetGlobal_f(I->G, cSetting_gamma);
  if(gamma > R_SMALL4)
    gamma = 1.0F / gamma;
  else
    gamma = 1.0F;

  lightv = SettingGetfv(I->G, cSetting_light);
  copy3f(lightv, light);

  fog = SettingGetGlobal_f(I->G, cSetting_ray_trace_fog);
  if(fog < 0.0F)
    fog = SettingGetGlobal_b(I->G, cSetting_depth_cue);
  if(fog != 0.0F) {
    if(fog > 1.0F)
      fog = 1.0F;
    fog_start = SettingGetGlobal_f(I->G, cSetting_ray_trace_fog_start);
    if(fog_start < 0.0F)
      fog_start = SettingGetGlobal_f(I->G, cSetting_fog_start);
    if(fog_start > 1.0F)
      fog_start = 1.0F;
  }

  /* SETUP */

  if(antialias < 0)
    antialias = SettingGetGlobal_i(I->G, cSetting_antialias);

  bkrd = ColorGet(I->G, SettingGet_color(I->G, nullptr, nullptr, cSetting_bg_rgb));
  RayExpandPrimitives(I);
  RayTransformFirst(I, 0, identity);

  PRINTFB(I->G, FB_Ray, FB_Details)
    " RayRenderPovRay: processed %i graphics primitives.\n", I->NPrimitive ENDFB(I->G);
  base = I->Basis + 1;

  {
    int ortho = SettingGetGlobal_b(I->G, cSetting_ortho);
    if(!identity) {
      if(!ortho) {
        sprintf(buffer, "camera {direction<0.0,0.0,%8.3f>\n location <0.0 , 0.0 , 0.0>\n right %12.10f*x up y \n }\n", -57.3F * cos(fov * cPI / (180 * 2.4)) / fov,     /* by trial and error */
                I->Range[0] / I->Range[1]);
      } else {
        sprintf(buffer,
                "camera {orthographic location <0.0 , 0.0 , %12.10f>\nlook_at  <0.0 , 0.0 , -1.0> right %12.10f*x up %12.10f*y}\n",
                front, -I->Range[0], I->Range[1]);
      }
    } else {
      float look[3];
      float loc[3] = { 0.0F, 0.0F, 0.0F };
      SceneViewType view;

      SceneGetCenter(I->G, look);
      SceneGetView(I->G, view);
      loc[2] = -view[18];
      MatrixInvTransformC44fAs33f3f(view, loc, loc);
      add3f(look, loc, loc);

      if(!ortho) {
        sprintf(buffer,
                "camera {angle %12.10f sky<%12.10f,%12.10f,%12.10f>\nlocation<%12.10f,%12.10f,%12.10f>\nlook_at<%12.10f,%12.10f,%12.10f> right %12.10f*x up y }\n",
                fov * I->Range[0] / I->Range[1],
                view[1], view[5], view[9],
                loc[0], loc[1], loc[2],
                look[0], look[1], look[2], -I->Range[0] / I->Range[1]);
      } else {
        sprintf(buffer,
                "camera {orthographic sky<%12.10f,%12.10f,%12.10f>\nlocation<%12.10f,%12.10f,%12.10f>\nlook_at<%12.10f,%12.10f,%12.10f> right %12.10f*x up %12.10f*y}\n",
                view[1], view[5], view[9],
                loc[0], loc[1], loc[2],
                look[0], look[1], look[2], -I->Range[0], I->Range[1]);
      }
    }
    UtilConcatVLA(&headerVLA, &hc, buffer);
  }

  {
    float ambient =
      SettingGetGlobal_f(I->G, cSetting_ambient) + SettingGetGlobal_f(I->G, cSetting_direct);
    float reflect = SettingGetGlobal_f(I->G, cSetting_reflect);

    if(ambient > 0.5)
      ambient = 0.5;

    reflect = 1.2F - 1.5F * ambient;

    sprintf(buffer,
            "#default { finish{phong %8.3f ambient %8.3f diffuse %8.3f phong_size %8.6f}}\n",
            SettingGetGlobal_f(I->G, cSetting_spec_reflect), ambient, reflect, spec_power);
    UtilConcatVLA(&headerVLA, &hc, buffer);
  }

  if(!identity) {
    if(angle) {
      float temp[16];
      identity44f(temp);
      MatrixRotateC44f(temp, (float) -PI * angle / 180, 0.0F, 1.0F, 0.0F);
      MatrixTransformC44fAs33f3f(temp, light, light);
    }
  }

  {
    float lite[3];
    if(!identity) {
      lite[0] = -light[0] * 10000.0F;
      lite[1] = -light[1] * 10000.0F;
      lite[2] = -light[2] * 10000.0F - front;
    } else {
      /* this is not correct unless the camera is in the default orientation */
      float look[3];
      SceneViewType view;

      SceneGetCenter(I->G, look);
      SceneGetView(I->G, view);

      lite[0] = -light[0] * 10000.0F;
      lite[1] = -light[1] * 10000.0F;
      lite[2] = -light[2] * 10000.0F;
      MatrixInvTransformC44fAs33f3f(view, lite, lite);
      add3f(look, lite, lite);
    }
    sprintf(buffer, "light_source{<%6.4f,%6.4f,%6.4f>  rgb<1.0,1.0,1.0>}\n",
            lite[0], lite[1], lite[2]);
    UtilConcatVLA(&headerVLA, &hc, buffer);
  }

  if(!identity) {
    int opaque_back = SettingGetGlobal_i(I->G, cSetting_ray_opaque_background);
    if(opaque_back < 0)
      opaque_back = SettingGetGlobal_i(I->G, cSetting_opaque_background);

    if(opaque_back) {           /* drop a plane into the background for the background color */
      sprintf(buffer,
              "plane{z , %6.4f \n pigment{color rgb<%6.4f,%6.4f,%6.4f>}\n finish{phong 0 specular 0 diffuse 0 ambient 1.0}}\n",
              -back, bkrd[0], bkrd[1], bkrd[2]);
      UtilConcatVLA(&headerVLA, &hc, buffer);
    }
  }

  for(a = 0; a < I->NPrimitive; a++) {
    auto cap1 = cCylCap::Round;
    auto cap2 = cCylCap::Round;
    prim = I->Primitive + a;
    vert = base->Vertex + 3 * (prim->vert);
    if(prim->type == cPrimTriangle) {
      if(smooth_color_triangle)
        if(!mesh_obj) {
          sprintf(buffer, "mesh {\n");
          UtilConcatVLA(&charVLA, &cc, buffer);
          mesh_obj = true;
        }
    } else if(mesh_obj) {
      sprintf(buffer, " pigment{color rgb <1,1,1>}}");
      UtilConcatVLA(&charVLA, &cc, buffer);
      mesh_obj = false;
    }
    switch (prim->type) {
    case cPrimSphere:
      sprintf(buffer, "sphere{<%12.10f,%12.10f,%12.10f>, %12.10f\n",
              vert[0], vert[1], vert[2], prim->r1);
      UtilConcatVLA(&charVLA, &cc, buffer);
      sprintf(buffer, "pigment{color rgb<%6.4f,%6.4f,%6.4f>}}\n",
              prim->c1[0], prim->c1[1], prim->c1[2]);
      UtilConcatVLA(&charVLA, &cc, buffer);
      break;
    case cPrimCylinder:
      cap1 = prim->cap1;
      cap2 = prim->cap2;
    case cPrimSausage:
      d = base->Normal + 3 * base->Vert2Normal[prim->vert];
      scale3f(d, prim->l1, vert2);
      add3f(vert, vert2, vert2);
      sprintf(buffer,
              "cylinder{<%12.10f,%12.10f,%12.10f>,\n<%12.10f,%12.10f,%12.10f>,\n %12.10f\n",
              vert[0], vert[1], vert[2], vert2[0], vert2[1], vert2[2], prim->r1);
      UtilConcatVLA(&charVLA, &cc, buffer);

      if (cap1 != cCylCapFlat && cap2 != cCylCapFlat) {
        UtilConcatVLA(&charVLA, &cc, "open\n");
      }

      sprintf(buffer, "pigment{color rgb<%6.4f1,%6.4f,%6.4f>}}\n",
              (prim->c1[0] + prim->c2[0]) / 2,
              (prim->c1[1] + prim->c2[1]) / 2, (prim->c1[2] + prim->c2[2]) / 2);
      UtilConcatVLA(&charVLA, &cc, buffer);

      if (cap1 == cCylCapRound) {
      sprintf(buffer, "sphere{<%12.10f,%12.10f,%12.10f>, %12.10f\n",
              vert[0], vert[1], vert[2], prim->r1);
      UtilConcatVLA(&charVLA, &cc, buffer);
      sprintf(buffer, "pigment{color rgb<%6.4f1,%6.4f,%6.4f>}}\n",
              prim->c1[0], prim->c1[1], prim->c1[2]);
      UtilConcatVLA(&charVLA, &cc, buffer);
      }

      if (cap2 == cCylCapRound) {
      sprintf(buffer, "sphere{<%12.10f,%12.10f,%12.10f>, %12.10f\n",
              vert2[0], vert2[1], vert2[2], prim->r1);
      UtilConcatVLA(&charVLA, &cc, buffer);
      sprintf(buffer, "pigment{color rgb<%6.4f1,%6.4f,%6.4f>}}\n",
              prim->c2[0], prim->c2[1], prim->c2[2]);
      UtilConcatVLA(&charVLA, &cc, buffer);
      }

      break;
    case cPrimTriangle:
      norm = base->Normal + 3 * base->Vert2Normal[prim->vert] + 3;      /* first normal is the average */

      if(!TriangleDegenerate(vert, norm, vert + 3, norm + 3, vert + 6, norm + 6)) {

        if(smooth_color_triangle) {
          sprintf(buffer,
                  "smooth_color_triangle{<%12.10f,%12.10f,%12.10f>,\n<%12.10f,%12.10f,%12.10f>,\n<%6.4f1,%6.4f,%6.4f>,\n<%12.10f,%12.10f,%12.10f>,\n<%12.10f,%12.10f,%12.10f>,\n<%6.4f1,%6.4f,%6.4f>,\n<%12.10f,%12.10f,%12.10f>,\n<%12.10f,%12.10f,%12.10f>,\n<%6.4f1,%6.4f,%6.4f> }\n",
                  vert[0], vert[1], vert[2], norm[0], norm[1], norm[2], prim->c1[0],
                  prim->c1[1], prim->c1[2], vert[3], vert[4], vert[5], norm[3], norm[4],
                  norm[5], prim->c2[0], prim->c2[1], prim->c2[2], vert[6], vert[7],
                  vert[8], norm[6], norm[7], norm[8], prim->c3[0], prim->c3[1],
                  prim->c3[2]
            );
          UtilConcatVLA(&charVLA, &cc, buffer);
        } else {
          /* nowadays we use mesh2 to generate smooth_color_triangles */

          UtilConcatVLA(&charVLA, &cc, "mesh2 { ");
          sprintf(buffer,
                  "vertex_vectors { 3, <%12.10f,%12.10f,%12.10f>,\n<%12.10f,%12.10f,%12.10f>,\n<%12.10f,%12.10f,%12.10f>}\n normal_vectors { 3,\n<%12.10f,%12.10f,%12.10f>,\n<%12.10f,%12.10f,%12.10f>,\n<%12.10f,%12.10f,%12.10f>}\n",
                  vert[0], vert[1], vert[2], vert[3], vert[4], vert[5], vert[6], vert[7],
                  vert[8], norm[0], norm[1], norm[2], norm[3], norm[4], norm[5], norm[6],
                  norm[7], norm[8]
            );
          UtilConcatVLA(&charVLA, &cc, buffer);

          if(prim->trans > R_SMALL4)
            sprintf(transmit, "transmit %4.6f", prim->trans);
          else
            transmit[0] = 0;

          sprintf(buffer, "texture_list { 3, ");
          UtilConcatVLA(&charVLA, &cc, buffer);

          sprintf(buffer, "texture { pigment{color rgb<%6.4f1,%6.4f,%6.4f> %s}}\n",
                  prim->c1[0], prim->c1[1], prim->c1[2], transmit);
          UtilConcatVLA(&charVLA, &cc, buffer);

          sprintf(buffer, ",texture { pigment{color rgb<%6.4f1,%6.4f,%6.4f> %s}}\n",
                  prim->c2[0], prim->c2[1], prim->c2[2], transmit);
          UtilConcatVLA(&charVLA, &cc, buffer);

          sprintf(buffer, ",texture { pigment{color rgb<%6.4f1,%6.4f,%6.4f> %s}} }\n",
                  prim->c3[0], prim->c3[1], prim->c3[2], transmit);
          UtilConcatVLA(&charVLA, &cc, buffer);

          sprintf(buffer, "face_indices { 1, <0,1,2>, 0, 1, 2 } }\n");
          UtilConcatVLA(&charVLA, &cc, buffer);
        }
      }
      break;
    }
  }

  if(mesh_obj) {
    sprintf(buffer, " pigment{color rgb <1,1,1>}}");
    UtilConcatVLA(&charVLA, &cc, buffer);
    mesh_obj = false;
  }
  *charVLA_ptr = charVLA;
  *headerVLA_ptr = headerVLA;
}


/*========================================================================*/
void RayProjectTriangle(CRay * I, RayInfo * r, float *light, float *v0, float *n0,
                        float scale)
{
  float w2;
  float d1[3], d2[3], d3[3];
  float p1[3], p2[3], p3[3];
  int c = 0;
  const float _0 = 0.0F;
  float *impact = r->impact;

  if(dot_product3f(light, n0 - 3) >= _0)
    c++;
  else if(dot_product3f(light, n0) >= _0)
    c++;
  else if(dot_product3f(light, n0 + 3) >= _0)
    c++;
  else if(dot_product3f(light, n0 + 6) >= _0)
    c++;

  if(c) {

    w2 = 1.0F - (r->tri1 + r->tri2);

    subtract3f(v0, impact, d1);
    subtract3f(v0 + 3, impact, d2);
    subtract3f(v0 + 6, impact, d3);
    project3f(d1, n0, p1);
    project3f(d2, n0 + 3, p2);
    project3f(d3, n0 + 6, p3);
    scale3f(p1, w2, d1);
    scale3f(p2, r->tri1, d2);
    scale3f(p3, r->tri2, d3);
    add3f(d1, d2, d2);
    add3f(d2, d3, d3);
    scale3f(d3, scale, d3);
    if(dot_product3f(r->surfnormal, d3) >= _0)
      add3f(d3, impact, impact);
  }
}

#ifndef _PYMOL_NOPY
static void RayHashSpawn(CRayHashThreadInfo * Thread, int n_thread, int n_total)
{
  int blocked;
  PyObject *info_list;
  int a, c, n = 0;
  CRay *I = Thread->ray;
  PyMOLGlobals *G = I->G;

  blocked = PAutoBlock(G);

  PRINTFB(I->G, FB_Ray, FB_Blather)
    " Ray: filling voxels with %d threads...\n", n_thread ENDFB(I->G);
  while(n < n_total) {
    c = n;
    info_list = PyList_New(n_thread);
    for(a = 0; a < n_thread; a++) {
      if((c + a) < n_total) {
        PyList_SetItem(
            info_list, a, PyCapsule_New(Thread + c + a, nullptr, nullptr));
      } else {
        PyList_SetItem(info_list, a, PConvAutoNone(NULL));
      }
      n++;
    }
    PXDecRef(PYOBJECT_CALLMETHOD
             (G->P_inst->cmd, "_ray_hash_spawn", "OO", info_list, G->P_inst->cmd));
    Py_DECREF(info_list);
  }
  PAutoUnblock(G, blocked);
}
#endif

#ifndef _PYMOL_NOPY
static void RayAntiSpawn(CRayAntiThreadInfo * Thread, int n_thread)
{
  int blocked;
  PyObject *info_list;
  int a;
  CRay *I = Thread->ray;
  PyMOLGlobals *G = I->G;

  blocked = PAutoBlock(G);

  PRINTFB(I->G, FB_Ray, FB_Blather)
    " Ray: antialiasing with %d threads...\n", n_thread ENDFB(I->G);
  info_list = PyList_New(n_thread);
  for(a = 0; a < n_thread; a++) {
    PyList_SetItem(info_list, a, PyCapsule_New(Thread + a, nullptr, nullptr));
  }
  PXDecRef(PYOBJECT_CALLMETHOD
           (G->P_inst->cmd, "_ray_anti_spawn", "OO", info_list, G->P_inst->cmd));
  Py_DECREF(info_list);
  PAutoUnblock(G, blocked);
}
#endif

int RayHashThread(CRayHashThreadInfo * T)
{
  BasisMakeMap(T->basis, T->vert2prim, T->prim, T->n_prim, T->clipBox, T->phase,
               cCache_ray_map, T->perspective, T->front, T->size_hint);

  /* utilize a little extra wasted CPU time in thread 0 which computes the smaller map... */
  if(!T->phase) {
    if (T->ray->bkgrd_data){
      fill_background_image(T->ray, T->image,  T->width, T->height, T->width * (unsigned int) T->height);
    } else if (T->bkrd_is_gradient){
      fill_gradient(T->ray, T->opaque_back, T->image, T->bkrd_top, T->bkrd_bottom, T->width, T->height, T->width * (unsigned int) T->height);      
    } else {
      fill(T->image, T->background, T->bytes);
    }
    RayComputeBox(T->ray);
  }
  return 1;
}

#ifndef _PYMOL_NOPY
static void RayTraceSpawn(CRayThreadInfo * Thread, int n_thread)
{
  int blocked;
  PyObject *info_list;
  int a;
  CRay *I = Thread->ray;
  PyMOLGlobals *G = I->G;

  blocked = PAutoBlock(G);

  PRINTFB(I->G, FB_Ray, FB_Blather)
    " Ray: rendering with %d threads...\n", n_thread ENDFB(I->G);
  info_list = PyList_New(n_thread);
  for(a = 0; a < n_thread; a++) {
    PyList_SetItem(info_list, a, PyCapsule_New(Thread + a, nullptr, nullptr));
  }
  PXDecRef(PYOBJECT_CALLMETHOD
           (G->P_inst->cmd, "_ray_spawn", "OO", info_list, G->P_inst->cmd));
  Py_DECREF(info_list);
  PAutoUnblock(G, blocked);

}
#endif

static int find_edge(unsigned int *ptr, float *depth, unsigned int width,
                     int threshold, int back)
{                               /* can only be called for a pixel NOT on the edge */
  {                             /* color testing */
    int compare0, compare1, compare2, compare3, compare4, compare5, compare6,
      compare7, compare8;
    {
      int back_test, back_two = false;
      compare0 = (signed int) *(ptr);
      compare1 = (signed int) *(ptr - 1);
      back_test = (compare0 == back);
      compare2 = (signed int) *(ptr + 1);
      back_two = back_two || ((compare1 == back) == back_test);
      compare3 = (signed int) *(ptr - width);
      back_two = back_two || ((compare2 == back) == back_test);
      compare4 = (signed int) *(ptr + width);
      back_two = back_two || ((compare3 == back) == back_test);
      compare5 = (signed int) *(ptr - width - 1);
      back_two = back_two || ((compare4 == back) == back_test);
      compare6 = (signed int) *(ptr + width - 1);
      back_two = back_two || ((compare5 == back) == back_test);
      compare7 = (signed int) *(ptr - width + 1);
      back_two = back_two || ((compare6 == back) == back_test);
      compare8 = (signed int) *(ptr + width + 1);
      back_two = back_two || ((compare7 == back) == back_test);
      if(back_two)
	threshold = (threshold >> 1);   // halve threshold for pixels that hit background 
    }
    
    {
      int current;
      unsigned int shift = 0;
      int sum1 = 0, sum2 = 3, sum3 = 0, sum4 = 0, sum5 = 0, sum6 = 0, sum7 = 0, sum8 = 0;
      int a;
      for(a = 0; a < 4; a++) {
	current = ((compare0 >> shift) & 0xFF);
	sum1 += abs(current - ((compare1 >> shift) & 0xFF));
	sum2 += abs(current - ((compare2 >> shift) & 0xFF));
	if(sum1 >= threshold)
	  return 1;
	sum3 += abs(current - ((compare3 >> shift) & 0xFF));
	if(sum2 >= threshold)
	  return 1;
	sum4 += abs(current - ((compare4 >> shift) & 0xFF));
	if(sum3 >= threshold)
	  return 1;
	sum5 += abs(current - ((compare5 >> shift) & 0xFF));
	if(sum4 >= threshold)
	  return 1;
	sum6 += abs(current - ((compare6 >> shift) & 0xFF));
	if(sum5 >= threshold)
	  return 1;
	sum7 += abs(current - ((compare7 >> shift) & 0xFF));
	if(sum6 >= threshold)
	  return 1;
	sum8 += abs(current - ((compare8 >> shift) & 0xFF));
	if(sum7 >= threshold)
	  return 1;
	if(sum8 >= threshold)
	  return 1;
	shift += 8;
      }
    }
  }

  if(depth) {                   /* depth testing */
    float compare0, compare1, compare2, compare3, compare4, compare5, compare6,
      compare7, compare8;
    float dcutoff = threshold / 128.0F;

    compare1 = *(depth - 1);
    compare0 = *(depth);
    compare2 = *(depth + 1);
    if(fabs(compare0 - compare1) > dcutoff)
      return 1;
    compare5 = *(depth - width - 1);
    if(fabs(compare0 - compare2) > dcutoff)
      return 1;
    compare3 = *(depth - width);
    if(fabs(compare0 - compare5) > dcutoff)
      return 1;
    compare7 = *(depth - width + 1);
    if(fabs(compare0 - compare3) > dcutoff)
      return 1;
    compare6 = *(depth + width - 1);
    if(fabs(compare0 - compare7) > dcutoff)
      return 1;
    compare4 = *(depth + width);
    if(fabs(compare0 - compare6) > dcutoff)
      return 1;
    compare8 = *(depth + width + 1);
    if(fabs(compare0 - compare4) > dcutoff)
      return 1;
    if(fabs(compare0 - compare8) > dcutoff)
      return 1;
    /*    if(fabs(compare0-compare1)>0.001F)
       printf("%8.3f \n",compare0-compare1);
       if(fabs(compare0-compare2)>0.001F)
       printf("%8.3f \n",compare0-compare2);
       if(fabs(compare0-compare3)>0.001F)
       printf("%8.3f \n",compare0-compare3);
       if(fabs(compare0-compare4)>0.001F)
       printf("%8.3f \n",compare0-compare4);
       if(fabs(compare0-compare5)>0.001F)
       printf("%8.3f \n",compare0-compare5);
       if(fabs(compare0-compare6)>0.001F)
       printf("%8.3f \n",compare0-compare6);
       if(fabs(compare0-compare7)>0.001F)
       printf("%8.3f \n",compare0-compare7);
       if(fabs(compare0-compare8)>0.001F)
       printf("%8.3f \n",compare0-compare8); */

  }
  return 0;
}

static void RayPrimGetColorRamped(PyMOLGlobals * G, float *matrix, RayInfo * r, float *fc)
{
  float fc1[3], fc2[3], fc3[3];
  float *c1, *c2, *c3, w2;
  float back_pact[3];
  const float _0 = 0.0F, _1 = 1.0F, _01 = 0.1F;
  CPrimitive *lprim = r->prim;
  inverse_transformC44f3f(matrix, r->impact, back_pact);

  switch (lprim->type) {
  case cPrimTriangle:
    w2 = 1.0F - (r->tri1 + r->tri2);
    c1 = lprim->c1;
    if(c1[0] <= _0) {
      ColorGetRamped(G, (int) (c1[0] - _01), back_pact, fc1, -1);
      c1 = fc1;
    }
    c2 = lprim->c2;
    if(c2[0] <= _0) {
      ColorGetRamped(G, (int) (c2[0] - _01), back_pact, fc2, -1);
      c2 = fc2;
    }
    c3 = lprim->c3;
    if(c3[0] <= _0) {
      ColorGetRamped(G, (int) (c3[0] - _01), back_pact, fc3, -1);
      c3 = fc3;
    }
    fc[0] = (c2[0] * r->tri1) + (c3[0] * r->tri2) + (c1[0] * w2);
    fc[1] = (c2[1] * r->tri1) + (c3[1] * r->tri2) + (c1[1] * w2);
    fc[2] = (c2[2] * r->tri1) + (c3[2] * r->tri2) + (c1[2] * w2);
    break;
  case cPrimSphere:
    c1 = lprim->c1;
    if(c1[0] <= _0) {
      ColorGetRamped(G, (int) (c1[0] - _01), back_pact, fc1, -1);
      c1 = fc1;
    }
    copy3f(c1, fc);
    break;
  case cPrimEllipsoid:
    /* TO DO */
    break;
  case cPrimCone:
  case cPrimCylinder:
  case cPrimSausage:
    w2 = r->tri1;
    c1 = lprim->c1;
    if(c1[0] <= _0) {
      ColorGetRamped(G, (int) (c1[0] - _01), back_pact, fc1, -1);
      c1 = fc1;
    }
    c2 = lprim->c2;
    if(c2[0] <= _0) {
      ColorGetRamped(G, (int) (c2[0] - _01), back_pact, fc2, -1);
      c2 = fc2;
    }
    /* TODO : here is where we would set the correct color if we implemented
       non-linear (smoothstep) half bonds. We probably need to add another
       parameter/variable to the primitive to tell this function what interpolation
       to do */
    fc[0] = (c1[0] * (_1 - w2)) + (c2[0] * w2);
    fc[1] = (c1[1] * (_1 - w2)) + (c2[1] * w2);
    fc[2] = (c1[2] * (_1 - w2)) + (c2[2] * w2);
    break;
  default:
    fc[0] = _1;
    fc[1] = _1;
    fc[2] = _1;
    break;
  }
}

int RayTraceThread(CRayThreadInfo * T)
{
  CRay *I = T->ray;
  int x, y, yy;
  float excess = 0.0F;
  float dotgle;
  float bright, direct_cmp, reflect_cmp, fc[4];
  float ambient, direct, lreflect, ft, ffact = 0.0F, ffact1m;
  unsigned int cc0, cc1, cc2, cc3;
  int i;
  RayInfo r1, r2;
  int fogFlag = false;
  int fogRangeFlag = false;
  int opaque_back = 0, orig_opaque_back = 0;
  int n_hit = 0;
  int two_sided_lighting;
  float fog;
  float inter[3] = { 0.0F, 0.0F, 0.0F };
  float fog_start = 0.0F;
  /*    float gamma,inp,sig=1.0F; */
  float persist, persist_inv;
  float new_front;
  int pass;
  unsigned int last_pixel = 0, *pixel;
  int exclude1, exclude2;
  float lit;
  int backface_cull;
  float project_triangle;
  float excl_trans;
  int shadows;
  int trans_shadows;
  int trans_mode;
  float first_excess;
  int pixel_flag;
  float ray_trans_spec, ray_lab_spec;
  float shadow_fudge;
  int label_shadow_mode;
  int interior_color;
  int interior_flag;
  int interior_shadows;
  int interior_wobble;
  int interior_mode;
  float interior_reflect;
  int wobble_save;
  float settingPower, settingReflectPower, settingSpecPower, settingSpecReflect,
    settingSpecDirect;
  float settingSpecDirectPower;
  float invHgt, invFrontMinusBack, inv1minusFogStart, invWdth, invHgtRange;
  float invWdthRange, vol0;
  float vol2;
  CBasis *bp1, *bp2;
  int render_height;
  int offset = 0;
  BasisCallRec BasisCall[MAX_BASIS];
  float border_offset;
  int edge_sampling = false;
  unsigned int edge_avg[4] = { 0, 0, 0, 0 };
  unsigned int edge_alpha_avg[4] = { 0, 0, 0, 0 };
  int edge_cnt = 0;
  float edge_base[2] = { 0.0F, 0.0F };
  float interior_normal[3] = {0.0F, 0.0F, 0.0F};
  float edge_width = 0.35356F;
  float edge_height = 0.35356F;
  float trans_spec_cut, trans_spec_scale, trans_oblique, oblique_power;
  float direct_shade;
  float red_blend = 0.0F;
  float blue_blend = 0.0F;
  float green_blend = 0.0F;
  float trans_cont;
  float pixel_base[3];
  float inv_trans_cont = 1.0F;
  float trans_cutoff, persist_cutoff;
  int trans_cont_flag = false;
  int blend_colors;
  int max_pass;
  float BasisFudge0, BasisFudge1;
  int perspective = T->perspective;
  float eye[3] = { 0.0F, 0.0F, 0.0F };
  float start[3] = { 0.0F, 0.0F, 0.0F };
  float nudge[3] = { 0.0F, 0.0F, 0.0F };
  float back_pact[3];
  float *depth = T->depth;
  float ray_scatter = SettingGetGlobal_f(I->G, cSetting_ray_scatter);
  const float shadow_decay = SettingGetGlobal_f(I->G, cSetting_ray_shadow_decay_factor);
  const float shadow_range = SettingGetGlobal_f(I->G, cSetting_ray_shadow_decay_range);
  const int clip_shadows = SettingGetGlobal_b(I->G, cSetting_ray_clip_shadows);
  const int spec_local = SettingGetGlobal_i(I->G, cSetting_ray_spec_local);
  float legacy = SettingGetGlobal_f(I->G, cSetting_ray_legacy_lighting);
  int spec_count = SettingGetGlobal_i(I->G, cSetting_spec_count);
  const float _0 = 0.0F;
  const float _1 = 1.0F;
  const float _p5 = 0.5F;
  const float _2 = 2.0F;
  const float _255 = 255.0F;
  const float _p499 = 0.499F;
  const float _persistLimit = 0.0001F;
  float legacy_1m = _1 - legacy;
  int n_basis = I->NBasis;
  int bg_image_mode = SettingGetGlobal_i(I->G, cSetting_bg_image_mode);
  int bg_image_linear = SettingGetGlobal_b(I->G, cSetting_bg_image_linear);
  auto bg_image_tilesize = SettingGet<const float*>(I->G, cSetting_bg_image_tilesize);
  float hpixelx = floor(T->width/2.f)  - floor(T->bgWidth / 2.f),
    hpixely = floor(T->height/2.f) - floor(T->bgHeight / 2.f);
  float hl;
  float wr = T->bgWidth/(float)T->width, hr = T->bgHeight/(float)T->height;
  float bg_rgb[3];
  const float *tmpf;
  int chromadepth;

  tmpf = ColorGet(I->G, SettingGet_color(I->G, nullptr, nullptr, cSetting_bg_rgb));
  mult3f(tmpf, 255.f, bg_rgb);

  /*   MemoryDebugDump();
     printf("%d\n",sizeof(CPrimitive));
   */

  {
    float fudge = SettingGetGlobal_f(I->G, cSetting_ray_triangle_fudge);

    BasisFudge0 = 0.0F - fudge;
    BasisFudge1 = 1.0F + fudge;
  }
  if(spec_count < 0) {
    spec_count = SettingGetGlobal_i(I->G, cSetting_light_count);
  }
  /* SETUP */

  /*  if(T->n_thread>1)
     printf(" Ray: Thread %d: Spawned.\n",T->phase+1);
   */

  interior_shadows = SettingGetGlobal_i(I->G, cSetting_ray_interior_shadows);
  interior_wobble = SettingGetGlobal_i(I->G, cSetting_ray_interior_texture);
  interior_color = SettingGetGlobal_i(I->G, cSetting_ray_interior_color);
  interior_reflect = 1.0F - SettingGetGlobal_f(I->G, cSetting_ray_interior_reflect);
  interior_mode = SettingGetGlobal_i(I->G, cSetting_ray_interior_mode);
  label_shadow_mode = SettingGetGlobal_i(I->G, cSetting_label_shadow_mode);
  project_triangle = SettingGetGlobal_f(I->G, cSetting_ray_improve_shadows);
  shadows = SettingGetGlobal_i(I->G, cSetting_ray_shadows);
  trans_shadows = SettingGetGlobal_i(I->G, cSetting_ray_transparency_shadows);

  backface_cull = SettingGetGlobal_i(I->G, cSetting_backface_cull);
  opaque_back = SettingGetGlobal_i(I->G, cSetting_ray_opaque_background);
  if(opaque_back < 0)
    opaque_back = SettingGetGlobal_i(I->G, cSetting_opaque_background);
  orig_opaque_back = opaque_back;

  if (T->bkrd_data){
    opaque_back = 1;
  }

  two_sided_lighting = SettingGetGlobal_i(I->G, cSetting_two_sided_lighting);
  if(two_sided_lighting<0) {
    if(SettingGetGlobal_i(I->G, cSetting_surface_cavity_mode))
      two_sided_lighting = true;
    else
      two_sided_lighting = false;
  }
  
  ray_trans_spec = SettingGetGlobal_f(I->G, cSetting_ray_transparency_specular);
  ray_lab_spec = SettingGetGlobal_f(I->G, cSetting_ray_label_specular);
  trans_cont = SettingGetGlobal_f(I->G, cSetting_ray_transparency_contrast);
  trans_mode = SettingGetGlobal_i(I->G, cSetting_transparency_mode);
  trans_oblique = SettingGetGlobal_f(I->G, cSetting_ray_transparency_oblique);
  oblique_power = SettingGetGlobal_f(I->G, cSetting_ray_transparency_oblique_power);
  trans_cutoff = SettingGetGlobal_f(I->G, cSetting_ray_trace_trans_cutoff);
  persist_cutoff = SettingGetGlobal_f(I->G, cSetting_ray_trace_persist_cutoff);
  chromadepth = SettingGetGlobal_i(I->G, cSetting_chromadepth);

  if(trans_mode == 1)
    two_sided_lighting = true;
  if(trans_cont > 1.0F) {
    trans_cont_flag = true;
    inv_trans_cont = 1.0F / trans_cont;
  }
  ambient = T->ambient;
  /* divide up the reflected light component over all lights */
  {
    float reflect_scale = SceneGetReflectScaleValue(I->G, 10);
    lreflect = reflect_scale * (SettingGetGlobal_f(I->G, cSetting_reflect) - ray_scatter);
    if(lreflect < _0)
      lreflect = _0;
    ray_scatter = ray_scatter * reflect_scale;
  }
  direct = SettingGetGlobal_f(I->G, cSetting_direct);

  /* apply legacy adjustments */

  ambient *= (1.0F - legacy) + (legacy * (0.22F / 0.12F));
  lreflect *= (1.0F - legacy) + (legacy * (0.72F / 0.45F));
  direct *= (1.0F - legacy) + (legacy * (0.24F / 0.45F));

  direct_shade = SettingGetGlobal_f(I->G, cSetting_ray_direct_shade);
  trans_spec_cut = SettingGetGlobal_f(I->G, cSetting_ray_transparency_spec_cut);
  blend_colors = SettingGetGlobal_i(I->G, cSetting_ray_blend_colors);
  max_pass = SettingGetGlobal_i(I->G, cSetting_ray_max_passes);
  if(blend_colors) {
    red_blend = SettingGetGlobal_f(I->G, cSetting_ray_blend_red);
    green_blend = SettingGetGlobal_f(I->G, cSetting_ray_blend_green);
    blue_blend = SettingGetGlobal_f(I->G, cSetting_ray_blend_blue);
  }

  if(trans_spec_cut < _1)
    trans_spec_scale = _1 / (_1 - trans_spec_cut);
  else
    trans_spec_scale = _0;

  /* COOP */
  settingPower = SettingGetGlobal_f(I->G, cSetting_power);
  settingReflectPower = SettingGetGlobal_f(I->G, cSetting_reflect_power);

  SceneGetAdjustedLightValues(I->G,
      &settingSpecReflect,
      &settingSpecPower,
      &settingSpecDirect,
      &settingSpecDirectPower, 10);

  if((interior_color != -1) || (two_sided_lighting) || (trans_mode == 1)
     || I->CheckInterior)
    backface_cull = 0;

  shadow_fudge = SettingGetGlobal_f(I->G, cSetting_ray_shadow_fudge);

  inv1minusFogStart = _1;

  fog = SettingGetGlobal_f(I->G, cSetting_ray_trace_fog);
  if(fog < 0.0F) {
    if(SettingGetGlobal_b(I->G, cSetting_depth_cue)) {
      fog = SettingGetGlobal_f(I->G, cSetting_fog);
    } else
      fog = _0;
  }

  if(fog != _0) {
    fogFlag = true;
    fog_start = SettingGetGlobal_f(I->G, cSetting_ray_trace_fog_start);
    if(fog_start < 0.0F)
      fog_start = SettingGetGlobal_f(I->G, cSetting_fog_start);
    if(fog_start > 1.0F)
      fog_start = 1.0F;
    if(abs(fog_start) > R_SMALL4) {
      fogRangeFlag = true;
      if(fabs(fog_start - 1.0F) < R_SMALL4)     /* prevent div/0 */
        fogFlag = false;
    }
    inv1minusFogStart = _1 / (_1 - fog_start);
  }

  /* ray-trace */

  if(T->border) {
    invHgt = _1 / (float) (T->height - (3.0F + T->border));
    invWdth = _1 / (float) (T->width - (3.0F + T->border));
  } else {

    invHgt = _1 / (float) (T->height);
    invWdth = _1 / (float) (T->width);
  }

  if(perspective) {
    float height_range, width_range;

    /* subpixel-offsets for antialiasing naturally correspond to
       effective pixel sizes at the front of the visible slab... */

    height_range = (T->front) * 2 * ((float) tan((T->fov / 2.0F) * PI / 180.0F));
    width_range = height_range * (I->Range[0] / I->Range[1]);
    invWdthRange = invWdth * width_range;
    invHgtRange = invHgt * height_range;
    vol0 = eye[0] - width_range / 2.0F;
    vol2 = eye[1] - height_range / 2.0F;
  } else {
    invWdthRange = invWdth * I->Range[0];
    invHgtRange = invHgt * I->Range[1];
    vol0 = I->Volume[0];
    vol2 = I->Volume[2];
  }
  invFrontMinusBack = _1 / (T->front - T->back);

  edge_width *= invWdthRange;
  edge_height *= invHgtRange;

  bp1 = I->Basis + 1;
  if(I->NBasis > 2)
    bp2 = I->Basis + 2;
  else
    bp2 = nullptr;

  render_height = T->y_stop - T->y_start;

  if(render_height) {
    offset = (T->phase * render_height / T->n_thread);
    offset = offset - (offset % T->n_thread) + T->phase;
  }
  if((interior_color != -1) || I->CheckInterior) {

    if(interior_color != -1)
      ColorGetEncoded(I->G, interior_color, inter);
    if(bp2) {
      interior_normal[0] = interior_reflect * bp2->LightNormal[0];
      interior_normal[1] = interior_reflect * bp2->LightNormal[1];
      interior_normal[2] = 1.0F + interior_reflect * bp2->LightNormal[2];
    } else {
      interior_normal[0] = 0.0;
      interior_normal[1] = 0.0;
      interior_normal[2] = 1.0F;
    }
    normalize3f(interior_normal);
  }

  r1.base[2] = _0;

  int* vert2prim_ptr = I->Vert2Prim.empty() ? nullptr : I->Vert2Prim.data();

  BasisCall[0].Basis = I->Basis + 1;
  BasisCall[0].rr = &r1;
  BasisCall[0].vert2prim = vert2prim_ptr;
  BasisCall[0].prim = I->Primitive;
  BasisCall[0].shadow = false;
  BasisCall[0].back = T->back;
  BasisCall[0].trans_shadows = trans_shadows;
  BasisCall[0].nearest_shadow = (shadow_decay != _0);
  BasisCall[0].check_interior = ((interior_color != -1) || I->CheckInterior);
  BasisCall[0].fudge0 = BasisFudge0;
  BasisCall[0].fudge1 = BasisFudge1;

  MapCacheInit(&BasisCall[0].cache, I->Basis[1].Map, T->phase, cCache_map_scene_cache);

  if(shadows && (n_basis > 2)) {
    int bc;
    for(bc = 2; bc < n_basis; bc++) {
      BasisCall[bc].Basis = I->Basis + bc;
      BasisCall[bc].rr = &r2;
      BasisCall[bc].vert2prim = vert2prim_ptr;
      BasisCall[bc].prim = I->Primitive;
      BasisCall[bc].shadow = true;
      BasisCall[bc].front = _0;
      BasisCall[bc].back = _0;
      BasisCall[bc].excl_trans = _0;
      BasisCall[bc].trans_shadows = trans_shadows;
      BasisCall[bc].nearest_shadow = (shadow_decay != _0) || (clip_shadows);
      BasisCall[bc].check_interior = false;
      BasisCall[bc].fudge0 = BasisFudge0;
      BasisCall[bc].fudge1 = BasisFudge1;
      BasisCall[bc].label_shadow_mode = label_shadow_mode;
      MapCacheInit(&BasisCall[bc].cache, I->Basis[bc].Map, T->phase,
                   cCache_map_shadow_cache);
    }
  }

  if(T->border) {
    border_offset = -1.50F + T->border / 2.0F;
  } else {
    border_offset = 0.0F;
  }

  unsigned int back_mask = 0x00000000;
  if (T->bkrd_is_gradient){
    if(opaque_back) {
      if(I->BigEndian)
	back_mask = 0x000000FF;
      else
	back_mask = 0xFF000000;
    }
  }
  for(yy = T->y_start; (yy < T->y_stop); yy++) {
    float perc, bkrd[4] = {0.f, 0.f, 0.f, 1.f};
    unsigned int bkrd_value = 0;
    short isOutsideInY = 0;

    if(I->G->Interrupt)
      break;

    y = T->y_start + ((yy - T->y_start) + offset) % (render_height);    /* make sure threads write to different pages */
    if (T->bkrd_data){
      switch (bg_image_mode){
      case 1: // isCentered
	{
	  float tmpy = floor(y - hpixely);
	  isOutsideInY = (tmpy < 0.f || tmpy > (float)T->bgHeight);
	  hl = fmodpos(tmpy, (float)T->bgHeight);
	}
	break;
      case 2: // isTiled
	hl = (float)T->bgHeight * (fmodpos((float)y, bg_image_tilesize[1])/bg_image_tilesize[1]);
	break;
      case 3: // isCenteredRepeated
	hl = fmodpos(floor(y - hpixely), (float)T->bgHeight);
	break;
      default:
	hl = y * hr;
	break;
      }
    } else if (T->bkrd_is_gradient){
      /* for RayTraceThread, y is from bottom to top */
      perc = y/(float)T->height;
      bkrd[0] = T->bkrd_bottom[0] + perc * (T->bkrd_top[0] - T->bkrd_bottom[0]);
      bkrd[1] = T->bkrd_bottom[1] + perc * (T->bkrd_top[1] - T->bkrd_bottom[1]);
      bkrd[2] = T->bkrd_bottom[2] + perc * (T->bkrd_top[2] - T->bkrd_bottom[2]);
      bkrd[3] = 1.f;
      if(T->ray->BigEndian){
	bkrd_value = back_mask | 
	  ((0xFF & ((unsigned int) (bkrd[0] * 255 + _p499))) << 24) |
	  ((0xFF & ((unsigned int) (bkrd[1] * 255 + _p499))) << 16) |
	  ((0xFF & ((unsigned int) (bkrd[2] * 255 + _p499))) << 8);
      } else {
	bkrd_value = back_mask | 
	  ((0xFF & ((unsigned int) (bkrd[2] * 255 + _p499))) << 16) |
	  ((0xFF & ((unsigned int) (bkrd[1] * 255 + _p499))) << 8) |
	  ((0xFF & ((unsigned int) (bkrd[0] * 255 + _p499))));
      }
    } else {
      bkrd_value = T->background;
      bkrd[0] = T->bkrd_top[0];
      bkrd[1] = T->bkrd_top[1];
      bkrd[2] = T->bkrd_top[2];
      if (orig_opaque_back){
	bkrd[3] = 1.f;
      } else {
	bkrd[3] = 0.f;
      }
    }
    if((!T->phase) && !(yy & 0xF)) {    /* don't slow down rendering too much */
      if(T->edging_cutoff) {
        if(T->edging) {
          OrthoBusyFast(I->G, (int) (2.5F * T->height / 3 + 0.5F * y), 4 * T->height / 3);
        } else {
          OrthoBusyFast(I->G, (int) (T->height / 3 + 0.5F * y), 4 * T->height / 3);
        }
      } else {
        OrthoBusyFast(I->G, T->height / 3 + y, 4 * T->height / 3);
      }
    }
    pixel = T->image + (T->width * y) + T->x_start;

    if((y % T->n_thread) == T->phase) { /* this is my scan line */
      pixel_base[1] = ((y + 0.5F + border_offset) * invHgtRange) + vol2;

      for(x = T->x_start; (x < T->x_stop); x++) {
	if (T->bkrd_data){
	  // Need to compute background for every pixel if image-based
	  unsigned char bkrd_uc[4];
	  compute_background_for_pixel(bkrd_uc, isOutsideInY,
              bg_image_mode, bg_image_tilesize, bg_rgb, bg_image_linear,
              (unsigned char*) T->bkrd_data, T->bgWidth, T->bgHeight,
              x, wr, hl, hpixelx, orig_opaque_back );
	  bkrd[0] = bkrd_uc[0] / 255.f;
	  bkrd[1] = bkrd_uc[1] / 255.f;
	  bkrd[2] = bkrd_uc[2] / 255.f;
	  bkrd[3] = bkrd_uc[3] / 255.f;
	  if(T->ray->BigEndian){
	    bkrd_value = back_mask | 
	      ((0xFF & ((unsigned int) (bkrd_uc[0] * 255 + _p499))) << 24) |
	      ((0xFF & ((unsigned int) (bkrd_uc[1] * 255 + _p499))) << 16) |
	      ((0xFF & ((unsigned int) (bkrd_uc[2] * 255 + _p499))) << 8);
	  } else {
	    bkrd_value = back_mask | 
	      ((0xFF & ((unsigned int) (bkrd_uc[2] * 255 + _p499))) << 16) |
	      ((0xFF & ((unsigned int) (bkrd_uc[1] * 255 + _p499))) << 8) |
	      ((0xFF & ((unsigned int) (bkrd_uc[0] * 255 + _p499))));
	  }
	}

        pixel_base[0] = (((x + 0.5F + border_offset)) * invWdthRange) + vol0;

        while(1) {
          if(T->edging) {
            if(!edge_sampling) {
              if(x && y && (x < (T->width - 1)) && (y < (T->height - 1))) {     /* not on the edge... */
                if(find_edge(T->edging + (pixel - T->image),
                             depth + (pixel - T->image),
                             T->width, T->edging_cutoff, bkrd_value)) {
                  unsigned char *pixel_c = (unsigned char *) pixel;
                  unsigned int c1, c2, c3, c4;
                  edge_cnt = 1;
                  edge_sampling = true;

                  edge_avg[0] = (c1 = pixel_c[0]);
                  edge_avg[1] = (c2 = pixel_c[1]);
                  edge_avg[2] = (c3 = pixel_c[2]);
                  edge_avg[3] = (c4 = pixel_c[3]);

                  edge_alpha_avg[0] = c1 * c4;
                  edge_alpha_avg[1] = c2 * c4;
                  edge_alpha_avg[2] = c3 * c4;
                  edge_alpha_avg[3] = c4;

                  edge_base[0] = pixel_base[0];
                  edge_base[1] = pixel_base[1];
                }
              }
            }
            if(edge_sampling) {
              if(edge_cnt == 5) {
                /* done with edging, so store averaged value */

                unsigned char *pixel_c = (unsigned char *) pixel;
                unsigned int c1, c2, c3, c4;

                edge_sampling = false;
                /* done with edging, so store averaged value */

                if(edge_alpha_avg[3]) {
                  c4 = edge_alpha_avg[3];
                  c1 = edge_alpha_avg[0] / c4;
                  c2 = edge_alpha_avg[1] / c4;
                  c3 = edge_alpha_avg[2] / c4;
                  c4 /= edge_cnt;
                } else {
                  c1 = edge_avg[0] / edge_cnt;
                  c2 = edge_avg[1] / edge_cnt;
                  c3 = edge_avg[2] / edge_cnt;
                  c4 = edge_avg[3] / edge_cnt;
                }
                pixel_c[0] = c1;
                pixel_c[1] = c2;
                pixel_c[2] = c3;
                pixel_c[3] = c4;

		/* Need to account for alpha, if exists */
		if (orig_opaque_back && c4 < 255){
		  float bkpart = (255 - c4) / 255.f;
		  float ppart = 1.f - bkpart;
		  pixel_c[0] = pymol_roundf(bkpart * bkrd[0] * 255 + ppart * pixel_c[0]);
		  pixel_c[1] = pymol_roundf(bkpart * bkrd[1] * 255 + ppart * pixel_c[1]);
		  pixel_c[2] = pymol_roundf(bkpart * bkrd[2] * 255 + ppart * pixel_c[2]);
		  pixel_c[3] = 255;
		}

                /* restore X,Y coordinates */
                r1.base[0] = pixel_base[0];
                r1.base[1] = pixel_base[1];

              } else {
                *pixel = 0;
		//                *pixel = bkrd_value;
                switch (edge_cnt) {
                case 1:
                  r1.base[0] = edge_base[0] + edge_width;
                  r1.base[1] = edge_base[1] + edge_height;
                  break;
                case 2:
                  r1.base[0] = edge_base[0] + edge_width;
                  r1.base[1] = edge_base[1] - edge_height;
                  break;
                case 3:
                  r1.base[0] = edge_base[0] - edge_width;
                  r1.base[1] = edge_base[1] + edge_height;
                  break;
                case 4:
                  r1.base[0] = edge_base[0] - edge_width;
                  r1.base[1] = edge_base[1] - edge_height;
                  break;
                }

              }
            }
            if(!edge_sampling)  /* not oversampling this edge or already done... */
              break;
          } else {
            r1.base[0] = pixel_base[0];
            r1.base[1] = pixel_base[1];
          }

          exclude1 = -1;
          exclude2 = -1;
          persist = _1;
          first_excess = _0;
          excl_trans = _0;
          pass = 0;
          new_front = T->front;

          if(perspective) {
            r1.base[2] = -T->front;
            r1.dir[0] = (r1.base[0] - eye[0]);
            r1.dir[1] = (r1.base[1] - eye[1]);
            r1.dir[2] = (r1.base[2] - eye[2]);
            if(BasisCall[0].check_interior) {
              start[0] = r1.base[0];
              start[1] = r1.base[1];
              start[2] = r1.base[2];
            }
            normalize3f(r1.dir);
            {
              float scale = I->max_box[2] / r1.base[2];

              r1.skip[0] = r1.base[0] * scale;
              r1.skip[1] = r1.base[1] * scale;
              r1.skip[2] = I->max_box[2];
            }

          }
          while((persist > _persistLimit) && (pass <= max_pass)) {
            char fogFlagTmp = fogFlag;
            pixel_flag = false;
            BasisCall[0].except1 = exclude1;
            BasisCall[0].except2 = exclude2;
            BasisCall[0].front = new_front;
            BasisCall[0].excl_trans = excl_trans;
            BasisCall[0].interior_flag = false;
            BasisCall[0].pass = pass;

            if(perspective) {
              if(pass) {
                add3f(nudge, r1.base, r1.base);
                copy3f(r1.base, r1.skip);
              }
              BasisCall[0].back_dist = -(T->back + r1.base[2]) / r1.dir[2];
              i = BasisHitPerspective(&BasisCall[0]);
            } else {
              i = BasisHitOrthoscopic(&BasisCall[0]);
            }
            interior_flag = BasisCall[0].interior_flag && (!pass);

            if(((i >= 0) || interior_flag) && (pass < max_pass)) {
	      int n_basis_tmp = r1.prim->no_lighting ? 0 : n_basis;
              pixel_flag = true;
              n_hit++;
              if(((r1.trans = r1.prim->trans) != _0) && trans_cont_flag) {
                r1.trans = (float) pow(r1.trans, inv_trans_cont);
              }
              if(interior_flag) {
                copy3f(interior_normal, r1.surfnormal);
                if(perspective) {
                  copy3f(start, r1.impact);
                  r1.dist = _0;
                } else {
                  copy3f(r1.base, r1.impact);
                  r1.dist = T->front;
                  r1.impact[2] -= T->front;
                }

                if(interior_wobble >= 0) {
                  wobble_save = r1.prim->wobble;        /* This is a no-no for multithreading! */
                  r1.prim->wobble = interior_wobble;

                  RayReflectAndTexture(I, &r1, perspective);

                  r1.prim->wobble = wobble_save;
                } else
                  RayReflectAndTexture(I, &r1, perspective);

                dotgle = -r1.dotgle;
                if (interior_color >= 0 || interior_color <= cColorExtCutoff) {
                  copy3f(inter, fc);
                } else if (interior_color == cColorAtomic) {
                  // Affects spheres+cylinders, but not surfaces+cartoon.
                  copy3f(r1.prim->c1, fc);
                } else {
                  copy3f(r1.prim->ic, fc);
                }
              } else {
                if(!perspective)
                  new_front = r1.dist;

                switch (r1.prim->type) {
                case cPrimTriangle:

                  BasisGetTriangleNormal(bp1, &r1, i, fc, perspective);
                  r1.trans = (float) pow(r1.trans, inv_trans_cont);

                  if(r1.prim->ramped) {
                    RayPrimGetColorRamped(I->G, I->ModelView, &r1, fc);
                  }
                  if(bp2) {
                    RayProjectTriangle(I, &r1, bp2->LightNormal,
                                       bp1->Vertex + i * 3,
                                       bp1->Normal + bp1->Vert2Normal[i] * 3 + 3,
                                       project_triangle);
                  }

                  RayReflectAndTexture(I, &r1, perspective);
                  if(perspective) {
                    BasisGetTriangleFlatDotglePerspective(bp1, &r1, i);
                  } else {
                    BasisGetTriangleFlatDotgle(bp1, &r1, i);
                  }
                  break;
                case cPrimCharacter:
                  BasisGetTriangleNormal(bp1, &r1, i, fc, perspective);

                  r1.trans = CharacterInterpolate(I->G, r1.prim->char_id, fc);
		  fogFlagTmp = false;
                  RayReflectAndTexture(I, &r1, perspective);
                  BasisGetTriangleFlatDotgle(bp1, &r1, i);
                  break;

                case cPrimEllipsoid:

                  BasisGetEllipsoidNormal(bp1, &r1, i, perspective);
                  RayReflectAndTexture(I, &r1, perspective);

                  fc[0] = r1.prim->c1[0];
                  fc[1] = r1.prim->c1[1];
                  fc[2] = r1.prim->c1[2];
                  break;

                default:       /* sphere, cylinder, sausage, etc. */

                  /* must be a sphere (effectively speaking) */

                  if(perspective) {
                    RayGetSphereNormalPerspective(I, &r1);
                  } else {
                    RayGetSphereNormal(I, &r1);
                  }

                  RayReflectAndTexture(I, &r1, perspective);

                  if(r1.prim->ramped) {
                    RayPrimGetColorRamped(I->G, I->ModelView, &r1, fc);
                  } else {
                    switch (r1.prim->type) {
                    case cPrimCylinder:
                    case cPrimSausage:
                    case cPrimCone:
                      ft = r1.tri1;
                      fc[0] = (r1.prim->c1[0] * (_1 - ft)) + (r1.prim->c2[0] * ft);
                      fc[1] = (r1.prim->c1[1] * (_1 - ft)) + (r1.prim->c2[1] * ft);
                      fc[2] = (r1.prim->c1[2] * (_1 - ft)) + (r1.prim->c2[2] * ft);
                      break;
                    default:
                      fc[0] = r1.prim->c1[0];
                      fc[1] = r1.prim->c1[1];
                      fc[2] = r1.prim->c1[2];
                      break;
                    }
                  }
                  break;
                }

                if((trans_oblique != _0) && (r1.trans != _0)) {
                  if((r1.surfnormal[2] > _0) || two_sided_lighting) {
                    float oblique_factor = r1.surfnormal[2];
                    if(oblique_factor < _0)
                      oblique_factor = -oblique_factor;
                    if(oblique_factor != _1) {
                      if(oblique_factor > _p5) {
                        oblique_factor =
                          (float) (_p5 +
                                   _p5 * (_1 -
                                          pow((_1 - oblique_factor) * _2,
                                              oblique_power)));
                      } else {
                        oblique_factor =
                          (float) (_p5 * pow(oblique_factor * _2, oblique_power));
                      }
                    }
                    r1.trans *= (trans_oblique * oblique_factor) + (1.0F - trans_oblique);
                    if(r1.trans < 0.06F)
                      r1.trans = 0.06F; /* don't allow transparent to become opaque */
                  }
                }

                dotgle = -r1.dotgle;

                if(r1.flat_dotgle < _0) {
                  if((!two_sided_lighting) && (BasisCall[0].check_interior)
                     && (interior_mode != 2)) {
                    interior_flag = true;
                    copy3f(interior_normal, r1.surfnormal);
                    if(perspective) {
                      copy3f(start, r1.impact);
                      r1.dist = _0;
                    } else {
                      copy3f(r1.base, r1.impact);
                      r1.impact[2] -= T->front;
                      r1.dist = T->front;
                    }

                    if(interior_wobble >= 0) {
                      wobble_save = r1.prim->wobble;
                      r1.prim->wobble = interior_wobble;
                      RayReflectAndTexture(I, &r1, perspective);
                      r1.prim->wobble = wobble_save;
                    } else {
                      RayReflectAndTexture(I, &r1, perspective);
                    }

                    dotgle = -r1.dotgle;
                    if (interior_color >= 0 ||
                        interior_color <= cColorExtCutoff) {
                      copy3f(inter, fc);
                    } else if (interior_color == cColorAtomic) {
                      // Affects surfaces+cartoon, but not spheres+cylinders.
                      copy3f(r1.prim->c1, fc);
                    } else {
                      copy3f(r1.prim->ic, fc);
                    }
                  }
                }

                if((r1.flat_dotgle < _0) && (!interior_flag)) {
                  if(two_sided_lighting) {
                    dotgle = -dotgle;
                    invert3f(r1.surfnormal);
                  }
                }
              }

              {
                double pow_dotgle;
                float pow_surfnormal2;

                if(settingPower != _1) {
                  pow_dotgle = pow(dotgle, settingPower);
                  pow_surfnormal2 = (float) pow(r1.surfnormal[2], settingPower);
                } else {
                  pow_dotgle = dotgle;
                  pow_surfnormal2 = r1.surfnormal[2];
                }
                direct_cmp = legacy_1m * pow_surfnormal2 +      /* new model */
                  legacy * ((float) (dotgle + pow_dotgle) * _p5);       /* legacy model */
              }

              reflect_cmp = _0;
              if(settingSpecDirect != _0) {
                if(r1.surfnormal[2] > _0) {
                  excess =
                    (float) (pow(r1.surfnormal[2], settingSpecDirectPower) *
                             settingSpecDirect);
                } else {
                  excess = _0;
                }
              } else {
                excess = _0;
              }

              lit = _1;
              if(n_basis_tmp < 3) {
                reflect_cmp = direct_cmp;
              } else {
                int bc;
                CBasis *bp;
                for(bc = 2; bc < n_basis; bc++) {
                  lit = _1;
                  bp = I->Basis + bc;

                  if(shadows && ((!interior_flag) || (interior_shadows)) &&
                     ((r1.prim->type != cPrimCharacter) || (label_shadow_mode & 0x1))) {
                    matrix_transform33f3f(bp->Matrix, r1.impact, r2.base);
                    r2.base[2] -= shadow_fudge;
                    BasisCall[bc].except2 = -1;
                    BasisCall[bc].except1 = i;  /* exclude current prim from shadow comp */
                    if(BasisHitShadow(&BasisCall[bc]) > -1) {
                      if((!clip_shadows) || (bp->LightNormal[2] >= _0) ||
                         ((T->front + r1.impact[2] - (r2.dist * bp->LightNormal[2])) <
                          _0)) {
                        lit = (float) pow(r2.trans, _p5);
                        if((shadow_decay != _0) && (r2.dist > shadow_range)) {
                          if(shadow_decay > 0) {
                            lit +=
                              ((_1 - lit) * (_1 -
                                             _1 / exp((r2.dist - shadow_range) *
                                                      shadow_decay)));
                          } else {
                            lit +=
                              ((_1 - lit) * (_1 -
                                             _1 / pow(r2.dist / shadow_range,
                                                      -shadow_decay)));
                          }
                        }
                      }
                    }
                  }

                  if(lit > _0) {
                    {
                      double pow_dotgle;

                      dotgle = -dot_product3f(r1.surfnormal, bp->LightNormal);
                      if(dotgle < _0)
                        dotgle = _0;

                      if(settingReflectPower != _1)
                        pow_dotgle = pow(dotgle, settingReflectPower);
                      else
                        pow_dotgle = dotgle;

                      reflect_cmp += legacy_1m * ((float) (lit * pow_dotgle)) + /* new model */
                        legacy * ((float) (lit * (dotgle + pow_dotgle) * _p5)); /* legacy model */
                    }

                    if(bc < (spec_count + 2)) {

                      if(ray_scatter != _0)     /* scattered specular light */
                        excess += ray_scatter * dotgle;

                      if(spec_local && perspective) {
                        /* slower, C4D-like local specular */
                        float tmp[3];

                        add3f(r1.surfnormal, r1.surfnormal, tmp);
                        add3f(tmp, bp->LightNormal, tmp);
                        normalize3f(tmp);
                        dotgle = -(dot_product3f(r1.dir, tmp)) * 1.004F;
                        if(dotgle > _1)
                          dotgle = _1;
                        else if(dotgle < _0)
                          dotgle = _0;
                        dotgle = powf(dotgle, 0.29F);
                      } else {
                        dotgle = -dot_product3f(r1.surfnormal, bp->SpecNormal); /* fast OpenGL-like global specular */
                      }
                      if(dotgle < _0)
                        dotgle = _0;
                      if(r1.prim->type != cPrimCharacter) {
                        excess +=
                          (float) (pow(dotgle, settingSpecPower) * settingSpecReflect *
                                   lit);
                      } else {
                        excess +=
                          (float) (pow(dotgle, settingSpecPower) * settingSpecReflect *
                                   lit * ray_lab_spec);
                      }
                    }
                  }
                }
                if (chromadepth > 0) {
                  float d = -(r1.impact[2] + T->front) / (T->back - T->front);
                  const float margin = 1. / 4.;
                  float h = d * (1. + 2. * margin) - margin;
                  if (h < 0.0) h = 0.0;
                  if (h > 1.0) h = 1.0;
                  float hue6 = 6. * 240. / 360. * h;
                  float rgb[3];
                  rgb[0] = fabs(hue6 - 3.) - 1.;
                  rgb[1] = 2. - fabs(hue6 - 2.);
                  rgb[2] = 2. - fabs(hue6 - 4.);
                  clamp3f(rgb);
                  if (chromadepth == 2) {
                      float lum = 0.30 * fc[0] + 0.59 * fc[1] + 0.11 * fc[2];
                      rgb[0] *= lum;
                      rgb[1] *= lum;
                      rgb[2] *= lum;
                  }
                  copy3(rgb, fc);
                }
              }

              if(fc[0] <= ((float) cColorExtCutoff)) {  /* ramped color */
                inverse_transformC44f3f(I->ModelView, r1.impact, back_pact);
                ColorGetRamped(I->G, (int) (fc[0] - 0.1F), back_pact, fc, -1);
              }

              bright = ambient + (((_1 - direct_shade) + direct_shade * lit) * direct * direct_cmp + lreflect * reflect_cmp * (legacy_1m + legacy * direct_cmp));       /* blend legacy */
              if(excess > _1)
                excess = _1;
              if(bright > _1)
                bright = _1;
              else if(bright < _0)
                bright = _0;

              /*                      bright *= (_1-excess); */

              fc[0] = (bright * fc[0] + excess);
              fc[1] = (bright * fc[1] + excess);
              fc[2] = (bright * fc[2] + excess);

              if(n_basis_tmp && fogFlagTmp) {
                if(perspective) {
                  ffact = (T->front + r1.impact[2]) * invFrontMinusBack;
                } else {
                  ffact = (T->front - r1.dist) * invFrontMinusBack;
                }
                if(fogRangeFlag)
                  ffact = (ffact - fog_start) * inv1minusFogStart;

                ffact *= fog;

                if(ffact < _0)
                  ffact = _0;
                if(ffact > _1)
                  ffact = _1;

                ffact1m = _1 - ffact;

		if(orig_opaque_back) {
		  fc[0] = ffact * bkrd[0] + fc[0] * ffact1m;
		  fc[1] = ffact * bkrd[1] + fc[1] * ffact1m;
		  fc[2] = ffact * bkrd[2] + fc[2] * ffact1m;
		  fc[3] = _1;
		} else {
                  fc[3] = ffact1m * (_1 - r1.trans);
		}

                if(!pass) {
                  if(r1.trans < trans_spec_cut) {
                    first_excess = excess * ffact1m * ray_trans_spec;
                  } else {
                    first_excess = excess * ffact1m * ray_trans_spec *
                      trans_spec_scale * (_1 - r1.trans);
                  }
                } else {
                  fc[0] += first_excess;        /* dubious? */
                  fc[1] += first_excess;
                  fc[2] += first_excess;
                }
              } else {
                if(!pass) {
                  if(r1.trans < trans_spec_cut) {
                    first_excess = excess * ray_trans_spec;
                  } else {
                    first_excess = excess * ray_trans_spec *
                      trans_spec_scale * (_1 - r1.trans);
                  }
                } else {
                  fc[0] += first_excess;
                  fc[1] += first_excess;
                  fc[2] += first_excess;
                }
                if(orig_opaque_back) {
                  fc[3] = _1;
		} else {
		  fc[3] = _1 - r1.trans;
                }
              }
            } else if(pass) {
              /* hit nothing, and we're on on second or greater pass,
                 or we're on the last pass of a dead-end loop */
              i = -1;

	      fc[0] = first_excess + bkrd[0];
	      fc[1] = first_excess + bkrd[1];
	      fc[2] = first_excess + bkrd[2];

	      if(orig_opaque_back) {
                fc[3] = _1;
	      } else {
		fc[3] = _0;
		//		fc[3] = bkrd[3];
	      }

              ffact = 1.0F;
              ffact1m = 0.0F;

              pixel_flag = true;
              if(trans_cont_flag)
                persist = (float) pow(persist, trans_cont);

            }

            if(pixel_flag) {
              /*
                 inp    = (fc[0]+fc[1]+fc[2]) * _inv3;
                 if(inp < R_SMALL4) 
                 sig = _1;
                 else
                 sig = (float)(pow(inp,gamma) / inp);

                 cc0 = (uint)(sig * fc[0] * _255);
                 cc1 = (uint)(sig * fc[1] * _255);
                 cc2 = (uint)(sig * fc[2] * _255);
               */

              cc0 = (uint) (fc[0] * _255);
              cc1 = (uint) (fc[1] * _255);
              cc2 = (uint) (fc[2] * _255);

              if(cc0 > 255)
                cc0 = 255;
              if(cc1 > 255)
                cc1 = 255;
              if(cc2 > 255)
                cc2 = 255;

              if(orig_opaque_back) {
                if(I->BigEndian)
                  *pixel = T->fore_mask | (cc0 << 24) | (cc1 << 16) | (cc2 << 8);
                else
                  *pixel = T->fore_mask | (cc2 << 16) | (cc1 << 8) | cc0;
              } else {
                /* use alpha channel for fog with transparent backgrounds */
                cc3 = (uint) (fc[3] * _255);
                if(cc3 > 255)
                  cc3 = 255;

                if(I->BigEndian)
                  *pixel = (cc0 << 24) | (cc1 << 16) | (cc2 << 8) | cc3;
                else
                  *pixel = (cc3 << 24) | (cc2 << 16) | (cc1 << 8) | cc0;
              }
            }
            if(pass && persist < 1.f) {          /* average all four channels */
              float mix_in;
              if(i >= 0) {
                if(fogFlagTmp) {
                  if(trans_cont_flag && (ffact > _p5)) {
                    mix_in =
                      2 * (persist * (_1 - ffact) +
                           ((float) pow(persist, trans_cont) * (ffact - _p5)))
                      * (_1 - r1.trans * ffact);
                  } else {
                    mix_in = persist * (_1 - r1.trans * ffact);
                  }
                } else {
                  mix_in = persist * (_1 - r1.trans);
                }
              } else {
                mix_in = persist;
              }

              persist_inv = _1 - mix_in;

              if(!orig_opaque_back) {
                if(i < 0) {     /* hit nothing -- so don't blend */
                  fc[0] = (float) (0xFF & (last_pixel >> 24));
                  fc[1] = (float) (0xFF & (last_pixel >> 16));
                  fc[2] = (float) (0xFF & (last_pixel >> 8));
                  fc[3] = (float) (0xFF & (last_pixel));
                  if(trans_cont_flag) { /* unless we are increasing contrast */
                    float m;
                    if(I->BigEndian) {
                      m = _1 - (float) (0xFF & (last_pixel)) / _255;
                    } else {
                      m = _1 - (float) (0xFF & (last_pixel >> 24)) / _255;
                    }
                    m = _1 - (float) pow(m, trans_cont);
                    if(I->BigEndian) {
                      fc[3] = m * _255 + _p499;
                    } else {
                      fc[0] = m * _255 + _p499;
                    }
                  }
                } else {        /* hit something -- so keep blend and compute cumulative alpha */

                  fc[0] =
                    (0xFF & ((*pixel) >> 24)) * mix_in +
                    (0xFF & (last_pixel >> 24)) * persist_inv;
                  fc[1] =
                    (0xFF & ((*pixel) >> 16)) * mix_in +
                    (0xFF & (last_pixel >> 16)) * persist_inv;
                  fc[2] =
                    (0xFF & ((*pixel) >> 8)) * mix_in +
                    (0xFF & (last_pixel >> 8)) * persist_inv;
                  fc[3] =
                    (0xFF & ((*pixel))) * mix_in + (0xFF & (last_pixel)) * persist_inv;

                  if(i >= 0) {  /* make sure opaque objects get opaque alpha */
                    float o1, o2;
                    float m;

                    if(I->BigEndian) {
                      o1 = (float) (0xFF & (last_pixel)) / _255;
                      o2 = (float) (0xFF & (*pixel)) / _255;
                    } else {
                      o1 = (float) (0xFF & (last_pixel >> 24)) / _255;
                      o2 = (float) (0xFF & ((*pixel) >> 24)) / _255;
                    }

                    if(o1 < o2) {       /* make sure o1 is largest opacity */
                      m = o1;
                      o1 = o2;
                      o2 = m;
                    }
                    m = o1 + (1.0F - o1) * o2;
                    if(I->BigEndian) {
                      fc[3] = m * _255 + _p499;
                    } else {
                      fc[0] = m * _255 + _p499;
                    }
                  }
                }
              } else {          /* opaque background, so just blend */
                fc[0] =
                  (0xFF & ((*pixel) >> 24)) * mix_in +
                  (0xFF & (last_pixel >> 24)) * persist_inv;
                fc[1] =
                  (0xFF & ((*pixel) >> 16)) * mix_in +
                  (0xFF & (last_pixel >> 16)) * persist_inv;
                fc[2] =
                  (0xFF & ((*pixel) >> 8)) * mix_in +
                  (0xFF & (last_pixel >> 8)) * persist_inv;
                fc[3] =
                  (0xFF & ((*pixel))) * mix_in + (0xFF & (last_pixel)) * persist_inv;
              }

              cc0 = (uint) (fc[0]);
              cc1 = (uint) (fc[1]);
              cc2 = (uint) (fc[2]);
              cc3 = (uint) (fc[3]);

              if(cc0 > 255)
                cc0 = 255;
              if(cc1 > 255)
                cc1 = 255;
              if(cc2 > 255)
                cc2 = 255;
              if(cc3 > 255)
                cc3 = 255;

              *pixel = (cc0 << 24) | (cc1 << 16) | (cc2 << 8) | cc3;

            }

            if(depth && (i >= 0) &&
               (r1.trans < trans_cutoff) && (persist > persist_cutoff)) {
              depth[pixel - T->image] = (T->front + r1.impact[2]);
            }

            if(i >= 0) {
              if(r1.prim->type == cPrimSausage) {       /* carry ray through the stick */
                if(perspective)
                  excl_trans = (2 * r1.surfnormal[2] * r1.prim->r1 / r1.dir[2]);
                else
                  excl_trans = new_front + (2 * r1.surfnormal[2] * r1.prim->r1);
              }

              if((!backface_cull) && (trans_mode != 2))
                persist = persist * r1.trans;
              else {
                if(!r1.prim->no_lighting && (persist < 0.9999F) && (r1.trans > 0.05F)) {
                  /* don't combine transparent surfaces */
                  *pixel = last_pixel;
                } else {
                  persist = persist * r1.trans;
                }
              }
            }

            if(i < 0) {         /* nothing hit */
              break;
            } else {
              if(perspective) {
                if(r1.prim->type != cPrimCharacter) {
                  float extend = r1.dist + 0.00001F;
                  scale3f(r1.dir, extend, nudge);
                } else {
                  float extend = r1.dist;
                  scale3f(r1.dir, extend, nudge);
                }
              }
              last_pixel = *pixel;
              exclude2 = exclude1;
              exclude1 = i;
              pass++;
            }

          }                     /* end of ray while */

          if(blend_colors) {

            float red_min = _0;
            float green_min = _0;
            float blue_min = _0;
            float red_part;
            float green_part;
            float blue_part;

            if(I->BigEndian) {
              fc[0] = (float) (0xFF & (*pixel >> 24));
              fc[1] = (float) (0xFF & (*pixel >> 16));
              fc[2] = (float) (0xFF & (*pixel >> 8));
              cc3 = (0xFF & (*pixel));
            } else {
              cc3 = (0xFF & (*pixel >> 24));
              fc[2] = (float) (0xFF & (*pixel >> 16));
              fc[1] = (float) (0xFF & (*pixel >> 8));
              fc[0] = (float) (0xFF & (*pixel));
            }

            red_part = red_blend * fc[0];
            green_part = green_blend * fc[1];
            blue_part = blue_blend * fc[2];

            red_min = (green_part > blue_part) ? green_part : blue_part;
            green_min = (red_part > blue_part) ? red_part : blue_part;
            blue_min = (green_part > red_part) ? green_part : red_part;

            if(fc[0] < red_min)
              fc[0] = red_min;
            if(fc[1] < green_min)
              fc[1] = green_min;
            if(fc[2] < blue_min)
              fc[2] = blue_min;

            cc0 = (uint) (fc[0]);
            cc1 = (uint) (fc[1]);
            cc2 = (uint) (fc[2]);

            if(cc0 > 255)
              cc0 = 255;
            if(cc1 > 255)
              cc1 = 255;
            if(cc2 > 255)
              cc2 = 255;

            if(I->BigEndian)
              *pixel = (cc0 << 24) | (cc1 << 16) | (cc2 << 8) | cc3;
            else
              *pixel = (cc3 << 24) | (cc2 << 16) | (cc1 << 8) | cc0;
          }

	  {
	    /* If final pixel has alpha, then the background should be 
	       blended into it */
	    unsigned char *pixel_c = (unsigned char *) pixel;
	    if (pixel_c[3] < 255){
	      float pixa = (pixel_c[3]/255.f);
	      float pixam1 = _1 - pixa;
	      float pixam1255 = pixam1 * 255.f;
	      if (T->bkrd_data || T->bkrd_is_gradient){
		pixel_c[0] = (unsigned char)pymol_roundf(pixel_c[0] * pixa + bkrd[0] * pixam1255);
		pixel_c[1] = (unsigned char)pymol_roundf(pixel_c[1] * pixa + bkrd[1] * pixam1255);
		pixel_c[2] = (unsigned char)pymol_roundf(pixel_c[2] * pixa + bkrd[2] * pixam1255);
		pixel_c[3] = (unsigned char)pymol_roundf(pixel_c[3] * pixa + bkrd[3] * pixam1255);
	      }
	    }
	  }

          if(!T->edging)
            break;
          /* if here, then we're edging...
             so accumulate averages */
          {

            unsigned char *pixel_c = (unsigned char *) pixel;
            unsigned int c1, c2, c3, c4;

            edge_avg[0] += (c1 = pixel_c[0]);
            edge_avg[1] += (c2 = pixel_c[1]);
            edge_avg[2] += (c3 = pixel_c[2]);
            edge_avg[3] += (c4 = pixel_c[3]);

            edge_alpha_avg[0] += c1 * c4;
            edge_alpha_avg[1] += c2 * c4;
            edge_alpha_avg[2] += c3 * c4;
            edge_alpha_avg[3] += c4;

            edge_cnt++;
          }

        }                       /* end of edging while */
        pixel++;
      }                         /* end of for */

    }
    /* end of if */
  }                             /* end of for */
  /*  if(T->n_thread>1) 
     printf(" Ray: Thread %d: Complete.\n",T->phase+1); */
  MapCacheFree(&BasisCall[0].cache, T->phase, cCache_map_scene_cache);

  if(shadows && (I->NBasis > 2)) {
    int bc;
    for(bc = 2; bc < I->NBasis; bc++) {
      MapCacheFree(&BasisCall[bc].cache, T->phase, cCache_map_shadow_cache);
    }
  }
  return (n_hit);
}


/* This is both an antialias and a slight blur */


/* for whatever reason, greatly GCC perfers a linear sequence of
   accumulates over a single large expression -- the difference is
   huge: over 10% !!! */

#define combine4by4(var,src,mask) { \
  var =  ((src)[0 ] & mask)   ; \
  var += ((src)[1 ] & mask)   ; \
  var += ((src)[2 ] & mask)   ; \
  var += ((src)[3 ] & mask)   ; \
  var += ((src)[4 ] & mask)   ; \
  var +=(((src)[5 ] & mask)*13) ; \
  var +=(((src)[6 ] & mask)*13) ; \
  var += ((src)[7 ] & mask)   ; \
  var += ((src)[8 ] & mask)    ; \
  var +=(((src)[9 ] & mask)*13) ; \
  var +=(((src)[10] & mask)*13) ; \
  var += ((src)[11] & mask)   ; \
  var += ((src)[12] & mask)   ; \
  var += ((src)[13] & mask)   ; \
  var += ((src)[14] & mask)   ; \
  var += ((src)[15] & mask)   ; \
  var = (var >> 6) & mask; \
}

#define combine5by5(var,src,mask) { \
  var =  ((src)[0 ] & mask)   ; \
  var += ((src)[1 ] & mask)   ; \
  var += ((src)[2 ] & mask)   ; \
  var += ((src)[3 ] & mask)   ; \
  var += ((src)[4 ] & mask)   ; \
  var += ((src)[5 ] & mask)   ; \
  var +=(((src)[6 ] & mask)*5); \
  var +=(((src)[7 ] & mask)*5); \
  var +=(((src)[8 ] & mask)*5); \
  var += ((src)[9 ] & mask)   ; \
  var += ((src)[10] & mask)   ; \
  var +=(((src)[11] & mask)*5); \
  var +=(((src)[12] & mask)*8); \
  var +=(((src)[13] & mask)*5); \
  var += ((src)[14] & mask)   ; \
  var += ((src)[15] & mask)   ; \
  var +=(((src)[16] & mask)*5); \
  var +=(((src)[17] & mask)*5); \
  var +=(((src)[18] & mask)*5); \
  var += ((src)[19] & mask)   ; \
  var += ((src)[20] & mask)   ; \
  var += ((src)[21] & mask)   ; \
  var += ((src)[22] & mask)   ; \
  var += ((src)[23] & mask)   ; \
  var += ((src)[24] & mask)   ; \
  var = (var >> 6) & mask; \
 }

#define combine6by6(var,src,mask) { \
  var =  ((src)[0 ] & mask)   ; \
  var += ((src)[1 ] & mask)   ; \
  var += ((src)[2 ] & mask)   ; \
  var += ((src)[3 ] & mask)   ; \
  var += ((src)[4 ] & mask)   ; \
  var += ((src)[5 ] & mask)   ; \
  var += ((src)[6 ] & mask)   ; \
  var +=(((src)[7 ] & mask)*5); \
  var +=(((src)[8 ] & mask)*7); \
  var +=(((src)[9 ] & mask)*7); \
  var +=(((src)[10] & mask)*5); \
  var += ((src)[11] & mask)   ; \
  var += ((src)[12] & mask)   ; \
  var +=(((src)[13] & mask)*7); \
  var +=(((src)[14] & mask)*8); \
  var +=(((src)[15] & mask)*8); \
  var +=(((src)[16] & mask)*7); \
  var += ((src)[17] & mask)   ; \
  var += ((src)[18] & mask)   ; \
  var +=(((src)[19] & mask)*7); \
  var +=(((src)[20] & mask)*8); \
  var +=(((src)[21] & mask)*8); \
  var +=(((src)[22] & mask)*7); \
  var += ((src)[23] & mask)   ; \
  var += ((src)[24] & mask)   ; \
  var +=(((src)[25] & mask)*5); \
  var +=(((src)[26] & mask)*7); \
  var +=(((src)[27] & mask)*7); \
  var +=(((src)[28] & mask)*5); \
  var += ((src)[29] & mask)   ; \
  var += ((src)[30] & mask)   ; \
  var += ((src)[31] & mask)   ; \
  var += ((src)[32] & mask)   ; \
  var += ((src)[33] & mask)   ; \
  var += ((src)[34] & mask)   ; \
  var += ((src)[35] & mask)   ; \
  var = (var >> 7) & mask; \
 }

#define m00FF 0x00FF
#define mFF00 0xFF00
#define mFFFF 0xFFFF

int RayAntiThread(CRayAntiThreadInfo * T)
{
  int src_row_pixels;

  unsigned int *pSrc;
  unsigned int *pDst;
  /*   unsigned int m00FF=0x00FF,mFF00=0xFF00,mFFFF=0xFFFF; */
  int width;
  int height;
  int x, y, yy;
  unsigned int *p;
  int offset = 0;
  CRay *I = T->ray;

  OrthoBusyFast(I->G, 9, 10);
  width = (T->width / T->mag) - 2;
  height = (T->height / T->mag) - 2;

  src_row_pixels = T->width;

  offset = (T->phase * height) / T->n_thread;
  offset = offset - (offset % T->n_thread) + T->phase;

  for(yy = 0; yy < height; yy++) {
    y = (yy + offset) % height; /* make sure threads write to different pages */

    if((y % T->n_thread) == T->phase) { /* this is my scan line */
      unsigned long c1, c2, c3, c4, a;
      unsigned char *c;

      pSrc = T->image + src_row_pixels * (y * T->mag);
      pDst = T->image_copy + width * y;
      switch (T->mag) {
      case 2:
        {
          for(x = 0; x < width; x++) {

            c = (unsigned char *) (p = pSrc + (x * T->mag));
            c1 = c2 = c3 = c4 = a = 0;

            c4 += (a = c[3]);
            c1 += c[0] * a;
            c2 += c[1] * a;
            c3 += c[2] * a;
            c += 4;
            c4 += (a = c[3]);
            c1 += c[0] * a;
            c2 += c[1] * a;
            c3 += c[2] * a;
            c += 4;
            c4 += (a = c[3]);
            c1 += c[0] * a;
            c2 += c[1] * a;
            c3 += c[2] * a;
            c += 4;
            c4 += (a = c[3]);
            c1 += c[0] * a;
            c2 += c[1] * a;
            c3 += c[2] * a;
            c += 4;

            c = (unsigned char *) (p += src_row_pixels);

            c4 += (a = c[3]);
            c1 += c[0] * a;
            c2 += c[1] * a;
            c3 += c[2] * a;
            c += 4;
            c4 += (a = c[3] * 13);
            c1 += c[0] * a;
            c2 += c[1] * a;
            c3 += c[2] * a;
            c += 4;
            c4 += (a = c[3] * 13);
            c1 += c[0] * a;
            c2 += c[1] * a;
            c3 += c[2] * a;
            c += 4;
            c4 += (a = c[3]);
            c1 += c[0] * a;
            c2 += c[1] * a;
            c3 += c[2] * a;
            c += 4;

            c = (unsigned char *) (p += src_row_pixels);

            c4 += (a = c[3]);
            c1 += c[0] * a;
            c2 += c[1] * a;
            c3 += c[2] * a;
            c += 4;
            c4 += (a = c[3] * 13);
            c1 += c[0] * a;
            c2 += c[1] * a;
            c3 += c[2] * a;
            c += 4;
            c4 += (a = c[3] * 13);
            c1 += c[0] * a;
            c2 += c[1] * a;
            c3 += c[2] * a;
            c += 4;
            c4 += (a = c[3]);
            c1 += c[0] * a;
            c2 += c[1] * a;
            c3 += c[2] * a;
            c += 4;

            c = (unsigned char *) (p += src_row_pixels);

            c4 += (a = c[3]);
            c1 += c[0] * a;
            c2 += c[1] * a;
            c3 += c[2] * a;
            c += 4;
            c4 += (a = c[3]);
            c1 += c[0] * a;
            c2 += c[1] * a;
            c3 += c[2] * a;
            c += 4;
            c4 += (a = c[3]);
            c1 += c[0] * a;
            c2 += c[1] * a;
            c3 += c[2] * a;
            c += 4;
            c4 += (a = c[3]);
            c1 += c[0] * a;
            c2 += c[1] * a;
            c3 += c[2] * a;
            c += 4;

            if(c4) {
              c1 /= c4;
              c2 /= c4;
              c3 /= c4;
            } else {            /* compute straight RGB average */

              c = (unsigned char *) (p = pSrc + (x * T->mag));
              c1 = c2 = c3 = 0;

              c1 += c[0];
              c2 += c[1];
              c3 += c[2];
              c += 4;
              c1 += c[0];
              c2 += c[1];
              c3 += c[2];
              c += 4;
              c1 += c[0];
              c2 += c[1];
              c3 += c[2];
              c += 4;
              c1 += c[0];
              c2 += c[1];
              c3 += c[2];
              c += 4;

              c = (unsigned char *) (p += src_row_pixels);

              c1 += c[0];
              c2 += c[1];
              c3 += c[2];
              c += 4;
              c1 += 13 * c[0];
              c2 += 13 * c[1];
              c3 += 13 * c[2];
              c += 4;
              c1 += 13 * c[0];
              c2 += 13 * c[1];
              c3 += 13 * c[2];
              c += 4;
              c1 += c[0];
              c2 += c[1];
              c3 += c[2];
              c += 4;

              c = (unsigned char *) (p += src_row_pixels);

              c1 += c[0];
              c2 += c[1];
              c3 += c[2];
              c += 4;
              c1 += 13 * c[0];
              c2 += 13 * c[1];
              c3 += 13 * c[2];
              c += 4;
              c1 += 13 * c[0];
              c2 += 13 * c[1];
              c3 += 13 * c[2];
              c += 4;
              c1 += c[0];
              c2 += c[1];
              c3 += c[2];
              c += 4;

              c = (unsigned char *) (p += src_row_pixels);

              c1 += c[0];
              c2 += c[1];
              c3 += c[2];
              c += 4;
              c1 += c[0];
              c2 += c[1];
              c3 += c[2];
              c += 4;
              c1 += c[0];
              c2 += c[1];
              c3 += c[2];
              c += 4;
              c1 += c[0];
              c2 += c[1];
              c3 += c[2];
              c += 4;

              c1 = c1 >> 6;
              c2 = c2 >> 6;
              c3 = c3 >> 6;
            }

            c = (unsigned char *) (pDst++);

            *(c++) = (unsigned char) c1;
            *(c++) = (unsigned char) c2;
            *(c++) = (unsigned char) c3;
            *(c++) = (unsigned char) (c4 >> 6);
          }
        }
        break;
      case 3:
        {
          for(x = 0; x < width; x++) {

            c = (unsigned char *) (p = pSrc + (x * T->mag));
            c1 = c2 = c3 = c4 = a = 0;

            c4 += (a = c[3]);
            c1 += c[0] * a;
            c2 += c[1] * a;
            c3 += c[2] * a;
            c += 4;
            c4 += (a = c[3]);
            c1 += c[0] * a;
            c2 += c[1] * a;
            c3 += c[2] * a;
            c += 4;
            c4 += (a = c[3]);
            c1 += c[0] * a;
            c2 += c[1] * a;
            c3 += c[2] * a;
            c += 4;
            c4 += (a = c[3]);
            c1 += c[0] * a;
            c2 += c[1] * a;
            c3 += c[2] * a;
            c += 4;
            c4 += (a = c[3]);
            c1 += c[0] * a;
            c2 += c[1] * a;
            c3 += c[2] * a;
            c += 4;

            c = (unsigned char *) (p += src_row_pixels);

            c4 += (a = c[3]);
            c1 += c[0] * a;
            c2 += c[1] * a;
            c3 += c[2] * a;
            c += 4;
            c4 += (a = c[3] * 5);
            c1 += c[0] * a;
            c2 += c[1] * a;
            c3 += c[2] * a;
            c += 4;
            c4 += (a = c[3] * 5);
            c1 += c[0] * a;
            c2 += c[1] * a;
            c3 += c[2] * a;
            c += 4;
            c4 += (a = c[3] * 5);
            c1 += c[0] * a;
            c2 += c[1] * a;
            c3 += c[2] * a;
            c += 4;
            c4 += (a = c[3]);
            c1 += c[0] * a;
            c2 += c[1] * a;
            c3 += c[2] * a;
            c += 4;

            c = (unsigned char *) (p += src_row_pixels);

            c4 += (a = c[3]);
            c1 += c[0] * a;
            c2 += c[1] * a;
            c3 += c[2] * a;
            c += 4;
            c4 += (a = c[3] * 5);
            c1 += c[0] * a;
            c2 += c[1] * a;
            c3 += c[2] * a;
            c += 4;
            c4 += (a = c[3] * 8);
            c1 += c[0] * a;
            c2 += c[1] * a;
            c3 += c[2] * a;
            c += 4;
            c4 += (a = c[3] * 5);
            c1 += c[0] * a;
            c2 += c[1] * a;
            c3 += c[2] * a;
            c += 4;
            c4 += (a = c[3]);
            c1 += c[0] * a;
            c2 += c[1] * a;
            c3 += c[2] * a;
            c += 4;

            c = (unsigned char *) (p += src_row_pixels);

            c4 += (a = c[3]);
            c1 += c[0] * a;
            c2 += c[1] * a;
            c3 += c[2] * a;
            c += 4;
            c4 += (a = c[3] * 5);
            c1 += c[0] * a;
            c2 += c[1] * a;
            c3 += c[2] * a;
            c += 4;
            c4 += (a = c[3] * 5);
            c1 += c[0] * a;
            c2 += c[1] * a;
            c3 += c[2] * a;
            c += 4;
            c4 += (a = c[3] * 5);
            c1 += c[0] * a;
            c2 += c[1] * a;
            c3 += c[2] * a;
            c += 4;
            c4 += (a = c[3]);
            c1 += c[0] * a;
            c2 += c[1] * a;
            c3 += c[2] * a;
            c += 4;

            c = (unsigned char *) (p += src_row_pixels);

            c4 += (a = c[3]);
            c1 += c[0] * a;
            c2 += c[1] * a;
            c3 += c[2] * a;
            c += 4;
            c4 += (a = c[3]);
            c1 += c[0] * a;
            c2 += c[1] * a;
            c3 += c[2] * a;
            c += 4;
            c4 += (a = c[3]);
            c1 += c[0] * a;
            c2 += c[1] * a;
            c3 += c[2] * a;
            c += 4;
            c4 += (a = c[3]);
            c1 += c[0] * a;
            c2 += c[1] * a;
            c3 += c[2] * a;
            c += 4;
            c4 += (a = c[3]);
            c1 += c[0] * a;
            c2 += c[1] * a;
            c3 += c[2] * a;
            c += 4;

            if(c4) {
              c1 /= c4;
              c2 /= c4;
              c3 /= c4;
            } else {            /* compute straight RGB average */

              c = (unsigned char *) (p = pSrc + (x * T->mag));
              c1 = c2 = c3 = 0;

              c1 += c[0];
              c2 += c[1];
              c3 += c[2];
              c += 4;
              c1 += c[0];
              c2 += c[1];
              c3 += c[2];
              c += 4;
              c1 += c[0];
              c2 += c[1];
              c3 += c[2];
              c += 4;
              c1 += c[0];
              c2 += c[1];
              c3 += c[2];
              c += 4;
              c1 += c[0];
              c2 += c[1];
              c3 += c[2];
              c += 4;

              c = (unsigned char *) (p += src_row_pixels);

              c1 += c[0];
              c2 += c[1];
              c3 += c[2];
              c += 4;
              c1 += 5 * c[0];
              c2 += 5 * c[1];
              c3 += 5 * c[2];
              c += 4;
              c1 += 5 * c[0];
              c2 += 5 * c[1];
              c3 += 5 * c[2];
              c += 4;
              c1 += 5 * c[0];
              c2 += 5 * c[1];
              c3 += 5 * c[2];
              c += 4;
              c1 += c[0];
              c2 += c[1];
              c3 += c[2];
              c += 4;

              c = (unsigned char *) (p += src_row_pixels);

              c1 += c[0];
              c2 += c[1];
              c3 += c[2];
              c += 4;
              c1 += 5 * c[0];
              c2 += 5 * c[1];
              c3 += 5 * c[2];
              c += 4;
              c1 += 8 * c[0];
              c2 += 8 * c[1];
              c3 += 8 * c[2];
              c += 4;
              c1 += 5 * c[0];
              c2 += 5 * c[1];
              c3 += 5 * c[2];
              c += 4;
              c1 += c[0];
              c2 += c[1];
              c3 += c[2];
              c += 4;

              c = (unsigned char *) (p += src_row_pixels);

              c1 += c[0];
              c2 += c[1];
              c3 += c[2];
              c += 4;
              c1 += 5 * c[0];
              c2 += 5 * c[1];
              c3 += 5 * c[2];
              c += 4;
              c1 += 5 * c[0];
              c2 += 5 * c[1];
              c3 += 5 * c[2];
              c += 4;
              c1 += 5 * c[0];
              c2 += 5 * c[1];
              c3 += 5 * c[2];
              c += 4;
              c1 += c[0];
              c2 += c[1];
              c3 += c[2];
              c += 4;

              c = (unsigned char *) (p += src_row_pixels);

              c1 += c[0];
              c2 += c[1];
              c3 += c[2];
              c += 4;
              c1 += c[0];
              c2 += c[1];
              c3 += c[2];
              c += 4;
              c1 += c[0];
              c2 += c[1];
              c3 += c[2];
              c += 4;
              c1 += c[0];
              c2 += c[1];
              c3 += c[2];
              c += 4;

              c1 = c1 >> 6;
              c2 = c2 >> 6;
              c3 = c3 >> 6;
            }

            c = (unsigned char *) (pDst++);

            *(c++) = (unsigned char) c1;
            *(c++) = (unsigned char) c2;
            *(c++) = (unsigned char) c3;
            *(c++) = (unsigned char) (c4 >> 6);
          }
        }
        break;
      case 4:
        {
          for(x = 0; x < width; x++) {

            c = (unsigned char *) (p = pSrc + (x * T->mag));
            c1 = c2 = c3 = c4 = a = 0;

            c4 += (a = c[3]);
            c1 += c[0] * a;
            c2 += c[1] * a;
            c3 += c[2] * a;
            c += 4;
            c4 += (a = c[3]);
            c1 += c[0] * a;
            c2 += c[1] * a;
            c3 += c[2] * a;
            c += 4;
            c4 += (a = c[3]);
            c1 += c[0] * a;
            c2 += c[1] * a;
            c3 += c[2] * a;
            c += 4;
            c4 += (a = c[3]);
            c1 += c[0] * a;
            c2 += c[1] * a;
            c3 += c[2] * a;
            c += 4;
            c4 += (a = c[3]);
            c1 += c[0] * a;
            c2 += c[1] * a;
            c3 += c[2] * a;
            c += 4;
            c4 += (a = c[3]);
            c1 += c[0] * a;
            c2 += c[1] * a;
            c3 += c[2] * a;
            c += 4;

            c = (unsigned char *) (p += src_row_pixels);

            c4 += (a = c[3]);
            c1 += c[0] * a;
            c2 += c[1] * a;
            c3 += c[2] * a;
            c += 4;
            c4 += (a = c[3] * 5);
            c1 += c[0] * a;
            c2 += c[1] * a;
            c3 += c[2] * a;
            c += 4;
            c4 += (a = c[3] * 7);
            c1 += c[0] * a;
            c2 += c[1] * a;
            c3 += c[2] * a;
            c += 4;
            c4 += (a = c[3] * 7);
            c1 += c[0] * a;
            c2 += c[1] * a;
            c3 += c[2] * a;
            c += 4;
            c4 += (a = c[3] * 5);
            c1 += c[0] * a;
            c2 += c[1] * a;
            c3 += c[2] * a;
            c += 4;
            c4 += (a = c[3]);
            c1 += c[0] * a;
            c2 += c[1] * a;
            c3 += c[2] * a;
            c += 4;

            c = (unsigned char *) (p += src_row_pixels);

            c4 += (a = c[3]);
            c1 += c[0] * a;
            c2 += c[1] * a;
            c3 += c[2] * a;
            c += 4;
            c4 += (a = c[3] * 7);
            c1 += c[0] * a;
            c2 += c[1] * a;
            c3 += c[2] * a;
            c += 4;
            c4 += (a = c[3] * 8);
            c1 += c[0] * a;
            c2 += c[1] * a;
            c3 += c[2] * a;
            c += 4;
            c4 += (a = c[3] * 8);
            c1 += c[0] * a;
            c2 += c[1] * a;
            c3 += c[2] * a;
            c += 4;
            c4 += (a = c[3] * 7);
            c1 += c[0] * a;
            c2 += c[1] * a;
            c3 += c[2] * a;
            c += 4;
            c4 += (a = c[3]);
            c1 += c[0] * a;
            c2 += c[1] * a;
            c3 += c[2] * a;
            c += 4;

            c = (unsigned char *) (p += src_row_pixels);

            c4 += (a = c[3]);
            c1 += c[0] * a;
            c2 += c[1] * a;
            c3 += c[2] * a;
            c += 4;
            c4 += (a = c[3] * 7);
            c1 += c[0] * a;
            c2 += c[1] * a;
            c3 += c[2] * a;
            c += 4;
            c4 += (a = c[3] * 8);
            c1 += c[0] * a;
            c2 += c[1] * a;
            c3 += c[2] * a;
            c += 4;
            c4 += (a = c[3] * 8);
            c1 += c[0] * a;
            c2 += c[1] * a;
            c3 += c[2] * a;
            c += 4;
            c4 += (a = c[3] * 7);
            c1 += c[0] * a;
            c2 += c[1] * a;
            c3 += c[2] * a;
            c += 4;
            c4 += (a = c[3]);
            c1 += c[0] * a;
            c2 += c[1] * a;
            c3 += c[2] * a;
            c += 4;

            c = (unsigned char *) (p += src_row_pixels);

            c4 += (a = c[3]);
            c1 += c[0] * a;
            c2 += c[1] * a;
            c3 += c[2] * a;
            c += 4;
            c4 += (a = c[3] * 5);
            c1 += c[0] * a;
            c2 += c[1] * a;
            c3 += c[2] * a;
            c += 4;
            c4 += (a = c[3] * 7);
            c1 += c[0] * a;
            c2 += c[1] * a;
            c3 += c[2] * a;
            c += 4;
            c4 += (a = c[3] * 7);
            c1 += c[0] * a;
            c2 += c[1] * a;
            c3 += c[2] * a;
            c += 4;
            c4 += (a = c[3] * 5);
            c1 += c[0] * a;
            c2 += c[1] * a;
            c3 += c[2] * a;
            c += 4;
            c4 += (a = c[3]);
            c1 += c[0] * a;
            c2 += c[1] * a;
            c3 += c[2] * a;
            c += 4;

            c = (unsigned char *) (p += src_row_pixels);

            c4 += (a = c[3]);
            c1 += c[0] * a;
            c2 += c[1] * a;
            c3 += c[2] * a;
            c += 4;
            c4 += (a = c[3]);
            c1 += c[0] * a;
            c2 += c[1] * a;
            c3 += c[2] * a;
            c += 4;
            c4 += (a = c[3]);
            c1 += c[0] * a;
            c2 += c[1] * a;
            c3 += c[2] * a;
            c += 4;
            c4 += (a = c[3]);
            c1 += c[0] * a;
            c2 += c[1] * a;
            c3 += c[2] * a;
            c += 4;
            c4 += (a = c[3]);
            c1 += c[0] * a;
            c2 += c[1] * a;
            c3 += c[2] * a;
            c += 4;
            c4 += (a = c[3]);
            c1 += c[0] * a;
            c2 += c[1] * a;
            c3 += c[2] * a;
            c += 4;

            if(c4) {
              c1 /= c4;
              c2 /= c4;
              c3 /= c4;
            } else {            /* compute straight RGB average */

              c = (unsigned char *) (p = pSrc + (x * T->mag));
              c1 = c2 = c3 = 0;

              c1 += c[0];
              c2 += c[1];
              c3 += c[2];
              c += 4;
              c1 += c[0];
              c2 += c[1];
              c3 += c[2];
              c += 4;
              c1 += c[0];
              c2 += c[1];
              c3 += c[2];
              c += 4;
              c1 += c[0];
              c2 += c[1];
              c3 += c[2];
              c += 4;
              c1 += c[0];
              c2 += c[1];
              c3 += c[2];
              c += 4;
              c1 += c[0];
              c2 += c[1];
              c3 += c[2];
              c += 4;

              c = (unsigned char *) (p += src_row_pixels);

              c1 += c[0];
              c2 += c[1];
              c3 += c[2];
              c += 4;
              c1 += 5 * c[0];
              c2 += 5 * c[1];
              c3 += 5 * c[2];
              c += 4;
              c1 += 7 * c[0];
              c2 += 7 * c[1];
              c3 += 7 * c[2];
              c += 4;
              c1 += 7 * c[0];
              c2 += 7 * c[1];
              c3 += 7 * c[2];
              c += 4;
              c1 += 5 * c[0];
              c2 += 5 * c[1];
              c3 += 5 * c[2];
              c += 4;
              c1 += c[0];
              c2 += c[1];
              c3 += c[2];
              c += 4;

              c = (unsigned char *) (p += src_row_pixels);

              c1 += c[0];
              c2 += c[1];
              c3 += c[2];
              c += 4;
              c1 += 7 * c[0];
              c2 += 7 * c[1];
              c3 += 7 * c[2];
              c += 4;
              c1 += 8 * c[0];
              c2 += 8 * c[1];
              c3 += 8 * c[2];
              c += 4;
              c1 += 8 * c[0];
              c2 += 8 * c[1];
              c3 += 8 * c[2];
              c += 4;
              c1 += 7 * c[0];
              c2 += 7 * c[1];
              c3 += 7 * c[2];
              c += 4;
              c1 += c[0];
              c2 += c[1];
              c3 += c[2];
              c += 4;

              c = (unsigned char *) (p += src_row_pixels);

              c1 += c[0];
              c2 += c[1];
              c3 += c[2];
              c += 4;
              c1 += 7 * c[0];
              c2 += 7 * c[1];
              c3 += 7 * c[2];
              c += 4;
              c1 += 8 * c[0];
              c2 += 8 * c[1];
              c3 += 8 * c[2];
              c += 4;
              c1 += 8 * c[0];
              c2 += 8 * c[1];
              c3 += 8 * c[2];
              c += 4;
              c1 += 7 * c[0];
              c2 += 7 * c[1];
              c3 += 7 * c[2];
              c += 4;
              c1 += c[0];
              c2 += c[1];
              c3 += c[2];
              c += 4;

              c = (unsigned char *) (p += src_row_pixels);

              c1 += c[0];
              c2 += c[1];
              c3 += c[2];
              c += 4;
              c1 += 5 * c[0];
              c2 += 5 * c[1];
              c3 += 5 * c[2];
              c += 4;
              c1 += 7 * c[0];
              c2 += 7 * c[1];
              c3 += 7 * c[2];
              c += 4;
              c1 += 7 * c[0];
              c2 += 7 * c[1];
              c3 += 7 * c[2];
              c += 4;
              c1 += 5 * c[0];
              c2 += 5 * c[1];
              c3 += 5 * c[2];
              c += 4;
              c1 += c[0];
              c2 += c[1];
              c3 += c[2];
              c += 4;

              c = (unsigned char *) (p += src_row_pixels);

              c1 += c[0];
              c2 += c[1];
              c3 += c[2];
              c += 4;
              c1 += c[0];
              c2 += c[1];
              c3 += c[2];
              c += 4;
              c1 += c[0];
              c2 += c[1];
              c3 += c[2];
              c += 4;
              c1 += c[0];
              c2 += c[1];
              c3 += c[2];
              c += 4;
              c1 += c[0];
              c2 += c[1];
              c3 += c[2];
              c += 4;

              c1 = c1 >> 7;
              c2 = c2 >> 7;
              c3 = c3 >> 7;
            }

            c = (unsigned char *) (pDst++);

            *(c++) = (unsigned char) c1;
            *(c++) = (unsigned char) c2;
            *(c++) = (unsigned char) c3;
            *(c++) = (unsigned char) (c4 >> 7);
          }
        }
        break;

      }
    }
  }
  return 1;
}

#ifdef PROFILE_BASIS
extern int n_cells;
extern int n_prims;
extern int n_triangles;
extern int n_spheres;
extern int n_cylinders;
extern int n_sausages;
extern int n_skipped;
#endif

float *rayDepthPixels = nullptr;
int rayVolume = 0, rayWidth = 0, rayHeight = 0;

/*========================================================================*/
void RayRender(CRay * I, unsigned int *image, double timing,
               float angle, int antialias, unsigned int *return_bg)
{
  int a, x, y;
  unsigned int *image_copy = nullptr;
  unsigned int back_mask, fore_mask = 0, trace_word = 0;
  unsigned int background;
  size_t buffer_size;
  int orig_opaque_back = 0, opaque_back = 0;
  int n_hit = 0;
  const float *bkrd_ptr;
  float bkrd_top[3], bkrd_bottom[3];
  short bkrd_is_gradient; /* if not gradient, use bkrd_top as bkrd */
  double now;
  int shadows;
  int n_thread;
  int mag = 1;
  int oversample_cutoff;
  int perspective = SettingGetGlobal_i(I->G, cSetting_ray_orthoscopic);
  int n_light = SettingGetGlobal_i(I->G, cSetting_light_count);
  float ambient;
  float *depth = nullptr;
  float front = I->Volume[4];
  float back = I->Volume[5];
  float fov = I->Fov;
  float *pos = glm::value_ptr(I->Pos);
  size_t width = I->Width;
  size_t height = I->Height;
  int ray_trace_mode;
  const float _0 = 0.0F, _p499 = 0.499F;
  int volume;
  const char * bg_image_filename;
  int ok = true;

  if(n_light > 10)
    n_light = 10;

  if(perspective < 0)
    perspective = SettingGetGlobal_b(I->G, cSetting_ortho);
  perspective = !perspective;

  VLACacheSize(I->G, I->Primitive, CPrimitive, I->NPrimitive, 0, cCache_ray_primitive);
#ifdef PROFILE_BASIS
  n_cells = 0;
  n_prims = 0;
  n_triangles = 0;
  n_spheres = 0;
  n_cylinders = 0;
  n_sausages = 0;
  n_skipped = 0;
#endif

  n_thread = SettingGetGlobal_i(I->G, cSetting_max_threads);
  if(n_thread < 1)
    n_thread = 1;
  if(n_thread > PYMOL_MAX_THREADS)
    n_thread = PYMOL_MAX_THREADS;
  opaque_back = SettingGetGlobal_i(I->G, cSetting_ray_opaque_background);
  if(opaque_back < 0)
    opaque_back = SettingGetGlobal_i(I->G, cSetting_opaque_background);

  orig_opaque_back = opaque_back;
  ray_trace_mode = SettingGetGlobal_i(I->G, cSetting_ray_trace_mode);

  shadows = SettingGetGlobal_i(I->G, cSetting_ray_shadows);

  oversample_cutoff = SettingGetGlobal_i(I->G, cSetting_ray_oversample_cutoff);

  if(antialias < 0) {
    antialias = SettingGetGlobal_i(I->G, cSetting_antialias);
  }

  if(ray_trace_mode && (antialias == 1))
    antialias = 2;
  else if(ray_trace_mode && antialias)
    antialias++;

  if(antialias < 0)
    antialias = 0;
  if(antialias > 4)
    antialias = 4;

  if((!antialias) || ray_trace_mode)
    oversample_cutoff = 0;

  mag = antialias;
  if(mag < 1)
    mag = 1;

  if(antialias > 1) {
    width = (width + 2) * mag;
    height = (height + 2) * mag;
    image_copy = image;
    buffer_size = width * height;
    image = CacheAlloc(I->G, unsigned int, buffer_size, 0, cCache_ray_antialias_buffer);
    ErrChkPtr(I->G, image);
  } else {
    buffer_size = width * height;
  }
  if(ray_trace_mode) {
    depth = pymol::calloc<float>(width * height);
  } else if(oversample_cutoff) {
    depth = pymol::calloc<float>(width * height);
  }
  ambient = SettingGetGlobal_f(I->G, cSetting_ambient);

  bkrd_is_gradient = SettingGetGlobal_b(I->G, cSetting_bg_gradient);
  bg_image_filename = SettingGet_s(I->G, nullptr, nullptr, cSetting_bg_image_filename);
  I->bkgrd_data = OrthoBackgroundDataGet(*I->G->Ortho);
  if (!I->bkgrd_data && bg_image_filename && bg_image_filename[0]){
    I->bkgrd_data = MyPNGRead(bg_image_filename);
    bkrd_is_gradient = 0;
  }
  if (I->bkgrd_data){
    opaque_back = 1;
  }

  if (!opaque_back){
    bkrd_is_gradient = 0;
  }
  if (bkrd_is_gradient){
    bkrd_ptr = ColorGet(I->G, SettingGet_color(I->G, nullptr, nullptr, cSetting_bg_rgb_top));
    copy3f(bkrd_ptr, bkrd_top);
    bkrd_ptr = ColorGet(I->G, SettingGet_color(I->G, nullptr, nullptr, cSetting_bg_rgb_bottom));
    copy3f(bkrd_ptr, bkrd_bottom);
  } else {
    bkrd_ptr = ColorGet(I->G, SettingGet_color(I->G, nullptr, nullptr, cSetting_bg_rgb));
    copy3f(bkrd_ptr, bkrd_top);
    copy3f(bkrd_ptr, bkrd_bottom);
  }
  {                             /* adjust bkrd and trace to offset the effect of gamma correction */
    float gamma = SettingGetGlobal_f(I->G, cSetting_gamma);
    float inp;
    float sig;
    inp = (bkrd_top[0] + bkrd_top[1] + bkrd_top[2]) / 3.0F;
    if(inp < R_SMALL4)
      sig = 1.0F;
    else
      sig = (float) (pow(inp, gamma)) / inp;
    bkrd_top[0] *= sig;
    bkrd_top[1] *= sig;
    bkrd_top[2] *= sig;
    if(bkrd_top[0] > 1.0F)
      bkrd_top[0] = 1.0F;
    if(bkrd_top[1] > 1.0F)
      bkrd_top[1] = 1.0F;
    if(bkrd_top[2] > 1.0F)
      bkrd_top[2] = 1.0F;

    if (bkrd_is_gradient) {
      float inp;
      float sig;
      inp = (bkrd_bottom[0] + bkrd_bottom[1] + bkrd_bottom[2]) / 3.0F;
      if(inp < R_SMALL4)
        sig = 1.0F;
      else
        sig = (float) (pow(inp, gamma)) / inp;
      bkrd_bottom[0] *= sig;
      bkrd_bottom[1] *= sig;
      bkrd_bottom[2] *= sig;
      if(bkrd_bottom[0] > 1.0F)
        bkrd_bottom[0] = 1.0F;
      if(bkrd_bottom[1] > 1.0F)
        bkrd_bottom[1] = 1.0F;
      if(bkrd_bottom[2] > 1.0F)
        bkrd_bottom[2] = 1.0F;
    } else {
      copy3f(bkrd_top, bkrd_bottom);      
    }
    if(ray_trace_mode) {
      float inp;
      float sig;
      int trace_color = SettingGetGlobal_color(I->G, cSetting_ray_trace_color);
      float trgb[3];
      const float *trgb_v = ColorGet(I->G, trace_color);

      copy3f(trgb_v, trgb);
      inp = (trgb[0] + trgb[1] + trgb[2]) / 3.0F;
      if(inp < R_SMALL4)
        sig = 1.0F;
      else
        sig = (float) (pow(inp, gamma)) / inp;
      trgb[0] *= sig;
      trgb[1] *= sig;
      trgb[2] *= sig;
      if(trgb[0] > 1.0F)
        trgb[0] = 1.0F;
      if(trgb[1] > 1.0F)
        trgb[1] = 1.0F;
      if(trgb[2] > 1.0F)
        trgb[2] = 1.0F;

      if(I->BigEndian) {
        trace_word =
          ((0xFF & ((unsigned int) (trgb[0] * 255 + _p499))) << 24) |
          ((0xFF & ((unsigned int) (trgb[1] * 255 + _p499))) << 16) |
          ((0xFF & ((unsigned int) (trgb[2] * 255 + _p499))) << 8) | 0xFF;
      } else {
        trace_word =
          0xFF000000 |
          ((0xFF & ((unsigned int) (trgb[2] * 255 + _p499))) << 16) |
          ((0xFF & ((unsigned int) (trgb[1] * 255 + _p499))) << 8) |
          ((0xFF & ((unsigned int) (trgb[0] * 255 + _p499))));
      }
    }
  }
  if(opaque_back) {
    if(I->BigEndian)
      back_mask = 0x000000FF;
    else
      back_mask = 0xFF000000;
    fore_mask = back_mask;
  } else {
    back_mask = 0x00000000;
  }
  if (!bkrd_is_gradient) {
    if(I->BigEndian) {
      background = back_mask |
	((0xFF & ((unsigned int) (bkrd_top[0] * 255 + _p499))) << 24) |
	((0xFF & ((unsigned int) (bkrd_top[1] * 255 + _p499))) << 16) |
	((0xFF & ((unsigned int) (bkrd_top[2] * 255 + _p499))) << 8);
    } else {
      background = back_mask |
	((0xFF & ((unsigned int) (bkrd_top[2] * 255 + _p499))) << 16) |
	((0xFF & ((unsigned int) (bkrd_top[1] * 255 + _p499))) << 8) |
	((0xFF & ((unsigned int) (bkrd_top[0] * 255 + _p499))));
    }
  } else {
    background = back_mask;
  }
  OrthoBusyFast(I->G, 2, 20);

  if (!bkrd_is_gradient) {
    PRINTFB(I->G, FB_Ray, FB_Blather)
      " RayNew: Background = %x %d %d %d\n", background, (int) (bkrd_top[0] * 255),
      (int) (bkrd_top[1] * 255), (int) (bkrd_top[2] * 255)
      ENDFB(I->G);
  }
  if(return_bg)
    *return_bg = background;

  if(!I->NPrimitive) {          /* nothing to render! */
    if (I->bkgrd_data){
      fill_background_image(I, image, width, height, width * (unsigned int) height);
    } else if (bkrd_is_gradient) {
      fill_gradient(I, opaque_back, image, bkrd_top, bkrd_bottom, width, height, width * (unsigned int) height);
    } else {
      fill(image, background, width * (unsigned int) height);
    }
  } else {

    if(I->PrimSizeCnt) {
      float factor = SettingGetGlobal_f(I->G, cSetting_ray_hint_camera);
      I->PrimSize = I->PrimSize / (I->PrimSizeCnt * factor);
      /*      printf("avg dist %8.7f\n",I->PrimSize); */
    } else {
      I->PrimSize = 0.0F;
    }
    ok &= !I->G->Interrupt;
    if (ok)
      ok &= RayExpandPrimitives(I);
    if (ok)
      ok &= RayTransformFirst(I, perspective, false);

    OrthoBusyFast(I->G, 3, 20);

    now = UtilGetSeconds(I->G) - timing;

    PRINTFB(I->G, FB_Ray, FB_Blather)
      " Ray: processed %i graphics primitives in %4.2f sec.\n", I->NPrimitive, now
      ENDFB(I->G);

    if (ok) {                           /* light sources */
      int bc;
      I->NBasis = n_light + 1;
      if(I->NBasis > MAX_BASIS)
	I->NBasis = MAX_BASIS;
      if(I->NBasis < 2)
	I->NBasis = 2;
      for(bc = 2; ok && bc < I->NBasis; bc++) {
        ok &= BasisInit(I->G, I->Basis + bc, bc);
      }
      for(bc = 2; ok && bc < I->NBasis; bc++) {
        {                       /* setup light & rotate if necessary  */
          float light[3];
          const float *lightv;
          switch (bc) {
          default:
          case 2:
            lightv = SettingGetfv(I->G, cSetting_light);
            break;
          case 3:
            lightv = SettingGetfv(I->G, cSetting_light2);
            break;
          case 4:
            lightv = SettingGetfv(I->G, cSetting_light3);
            break;
          case 5:
            lightv = SettingGetfv(I->G, cSetting_light4);
            break;
          case 6:
            lightv = SettingGetfv(I->G, cSetting_light5);
            break;
          case 7:
            lightv = SettingGetfv(I->G, cSetting_light6);
            break;
          case 8:
            lightv = SettingGetfv(I->G, cSetting_light7);
            break;
          case 9:
            lightv = SettingGetfv(I->G, cSetting_light8);
            break;
          case 10:
            lightv = SettingGetfv(I->G, cSetting_light9);
            break;
          }
          copy3f(lightv, light);
          normalize3f(light);

          if(angle) {
            float temp[16];
            identity44f(temp);
            MatrixRotateC44f(temp, (float) -PI * angle / 180, 0.0F, 1.0F, 0.0F);
            MatrixTransformC44fAs33f3f(temp, light, light);
          }

          I->Basis[bc].LightNormal[0] = light[0];
          I->Basis[bc].LightNormal[1] = light[1];
          I->Basis[bc].LightNormal[2] = light[2];
          normalize3f(I->Basis[bc].LightNormal);

          {
            float spec_vector[3];
            copy3f(I->Basis[bc].LightNormal, spec_vector);
            spec_vector[2]--;
            normalize3f(spec_vector);
            copy3f(spec_vector, I->Basis[bc].SpecNormal);
          }
        }

        if(ok && shadows) {           /* don't waste time on shadows unless needed */
          BasisSetupMatrix(I->Basis + bc);
          ok &= RayTransformBasis(I, I->Basis + bc, bc);
        }
	ok &= !I->G->Interrupt;
      }
    }

    OrthoBusyFast(I->G, 4, 20);
#ifndef _PYMOL_NOPY
    if(shadows && (n_thread > 1)) {     /* parallel execution */

      CRayHashThreadInfo *thread_info = pymol::calloc<CRayHashThreadInfo>(I->NBasis);

      /* rendering map */
      int* vert2prim_ptr = I->Vert2Prim.empty() ? nullptr : I->Vert2Prim.data();

      thread_info[0].basis = I->Basis + 1;
      thread_info[0].vert2prim = vert2prim_ptr;
      thread_info[0].prim = I->Primitive;
      thread_info[0].n_prim = I->NPrimitive;
      thread_info[0].clipBox = I->Volume;
      thread_info[0].phase = 0;
      thread_info[0].perspective = perspective;
      thread_info[0].front = front;

      thread_info[0].image = image;
      thread_info[0].bkrd_is_gradient = bkrd_is_gradient;
      thread_info[0].width = width;
      thread_info[0].height = height;
      if (I->bkgrd_data){
	thread_info[0].background = background;
	thread_info[0].opaque_back = opaque_back;
      } else if (bkrd_is_gradient){
	thread_info[0].bkrd_top = bkrd_top;
	thread_info[0].bkrd_bottom = bkrd_bottom;
	thread_info[0].opaque_back = opaque_back;
      } else {
	thread_info[0].background = background;
      }
      thread_info[0].bytes = width * (unsigned int) height;
      thread_info[0].ray = I;   /* for compute box */
      thread_info[0].size_hint = I->PrimSize;
      /* shadow map */

      {
        int bc;
        float factor = SettingGetGlobal_f(I->G, cSetting_ray_hint_shadow);
        for(bc = 2; bc < I->NBasis; bc++) {
          int* vert2prim_ptr = I->Vert2Prim.empty() ? nullptr : I->Vert2Prim.data();
          thread_info[bc - 1].basis = I->Basis + bc;
          thread_info[bc - 1].vert2prim = vert2prim_ptr;
          thread_info[bc - 1].prim = I->Primitive;
          thread_info[bc - 1].n_prim = I->NPrimitive;
          thread_info[bc - 1].clipBox = nullptr;
          thread_info[bc - 1].phase = bc - 1;
          thread_info[bc - 1].perspective = false;
          thread_info[bc - 1].front = _0;
          /* allowing these maps to be more fine helps performance */
          thread_info[bc - 1].size_hint = I->PrimSize * factor;
        }
      }

      /* NOTE that we're not limiting the number of threads in this phase
         under the assumption that it will usually just be a few threads */
      RayHashSpawn(thread_info, n_thread, I->NBasis - 1);

      FreeP(thread_info);
    } else
#endif
    if (ok){ 
#ifdef _PYMOL_NOPY
      n_thread = 1;          /* serial execution */
#endif
      int* vert2prim_ptr = I->Vert2Prim.empty() ? nullptr : I->Vert2Prim.data();
      ok &= BasisMakeMap(I->Basis + 1, vert2prim_ptr, I->Primitive, I->NPrimitive,
			 I->Volume, 0, cCache_ray_map, perspective, front, I->PrimSize);
      if(ok && shadows) {
        int bc;
        float factor = SettingGetGlobal_f(I->G, cSetting_ray_hint_shadow);
        for(bc = 2; ok && bc < I->NBasis; bc++) {
          ok &= BasisMakeMap(I->Basis + bc, vert2prim_ptr, I->Primitive, I->NPrimitive,
			     nullptr, bc - 1, cCache_ray_map, false, _0, I->PrimSize * factor);
        }
      }

      /* serial tasks which RayHashThread does in parallel mode using the first thread */

      if (ok){
	if (I->bkgrd_data){
	  fill_background_image(I, image,  width, height, width * (unsigned int) height);
	} else if (bkrd_is_gradient) {
	  fill_gradient(I, opaque_back, image, bkrd_top, bkrd_bottom, width, height, width * (unsigned int) height);
	} else {
	  fill(image, background, width * (unsigned int) height);
	}
	RayComputeBox(I);
	
      }
    }

    OrthoBusyFast(I->G, 5, 20);
    now = UtilGetSeconds(I->G) - timing;

    if (ok){
      if(shadows) {
	PRINTFB(I->G, FB_Ray, FB_Blather)
	  " Ray: voxels: [%4.2f:%dx%dx%d], [%4.2f:%dx%dx%d], %4.2f sec.\n",
	  I->Basis[1].Map->Div, I->Basis[1].Map->Dim[0],
	  I->Basis[1].Map->Dim[1], I->Basis[1].Map->Dim[2],
	  I->Basis[2].Map->Div, I->Basis[2].Map->Dim[0],
	  I->Basis[2].Map->Dim[2], I->Basis[2].Map->Dim[2], now ENDFB(I->G);
      } else {
	PRINTFB(I->G, FB_Ray, FB_Blather)
	  " Ray: voxels: [%4.2f:%dx%dx%d], %4.2f sec.\n",
	  I->Basis[1].Map->Div, I->Basis[1].Map->Dim[0],
	  I->Basis[1].Map->Dim[1], I->Basis[1].Map->Dim[2], now ENDFB(I->G);
      }
    }
    /* IMAGING */

    if (ok){
      /* now spawn threads as needed */
      CRayThreadInfo *rt = pymol::calloc<CRayThreadInfo>(n_thread);

      int x_start = 0, y_start = 0;
      int x_stop = 0, y_stop = 0;
      float x_test = _0, y_test = _0;
      int x_pixel, y_pixel;

      if(perspective) {
        int c;

        if(I->min_box[2] > -front)
          I->min_box[2] = -front;
        if(I->max_box[2] > -front)
          I->max_box[2] = -front;

        for(c = 0; c < 4; c++) {
          switch (c) {
          case 0:
            x_test = -I->min_box[0] / I->min_box[2];
            y_test = -I->min_box[1] / I->min_box[2];
            break;
          case 1:
            x_test = -I->min_box[0] / I->max_box[2];
            y_test = -I->min_box[1] / I->max_box[2];
            break;
          case 2:
            x_test = -I->max_box[0] / I->min_box[2];
            y_test = -I->max_box[1] / I->min_box[2];
            break;
          case 3:
            x_test = -I->max_box[0] / I->max_box[2];
            y_test = -I->max_box[1] / I->max_box[2];
            break;
          }

          /* project onto back to get the effective range */

          x_pixel =
            (int) (width * (((x_test * I->Volume[5]) - I->Volume[0]) / I->Range[0]));
          y_pixel =
            (int) (height * (((y_test * I->Volume[5]) - I->Volume[2]) / I->Range[1]));

          if(!c) {
            x_start = x_pixel;
            x_stop = x_pixel;
            y_start = y_pixel;
            y_stop = y_pixel;
          } else {
            if(x_start > x_pixel)
              x_start = x_pixel;
            if(x_stop < x_pixel)
              x_stop = x_pixel;
            if(y_start > y_pixel)
              y_start = y_pixel;
            if(y_stop < y_pixel)
              y_stop = y_pixel;
          }
        }
        x_start -= 2;
        x_stop += 2;
        y_start -= 2;
        y_stop += 2;

        /*
           x_start = 0; 
           y_start = 0;
           x_stop = width;
           y_stop = height; */

      } else {
        x_start = (int) ((width * (I->min_box[0] - I->Volume[0])) / I->Range[0]) - 2;
        x_stop = (int) ((width * (I->max_box[0] - I->Volume[0])) / I->Range[0]) + 2;

        y_stop = (int) ((height * (I->max_box[1] - I->Volume[2])) / I->Range[1]) + 2;
        y_start = (int) ((height * (I->min_box[1] - I->Volume[2])) / I->Range[1]) - 2;

      }
      if(x_start < 0)
        x_start = 0;
      if(y_start < 0)
        y_start = 0;
      if(x_stop > width)
        x_stop = width;
      if(y_stop > height)
        y_stop = height;

      for(a = 0; a < n_thread; a++) {
        rt[a].ray = I;
        rt[a].width = width;
        rt[a].height = height;
        rt[a].x_start = x_start;
        rt[a].x_stop = x_stop;
        rt[a].y_start = y_start;
        rt[a].y_stop = y_stop;
        rt[a].image = image;
        rt[a].border = mag - 1;
        rt[a].front = front;
        rt[a].back = back;
        rt[a].fore_mask = fore_mask;
	rt[a].bkrd_is_gradient = bkrd_is_gradient;
	if (bkrd_is_gradient){
	  rt[a].bkrd_top = bkrd_top;
	  rt[a].bkrd_bottom = bkrd_bottom;
	} else {
	  rt[a].bkrd_top = bkrd_top;
	  rt[a].bkrd_bottom = bkrd_top;
	}
        rt[a].ambient = ambient;
        rt[a].background = background;
        rt[a].phase = a;
        rt[a].n_thread = n_thread;
        rt[a].edging = nullptr;
        rt[a].edging_cutoff = oversample_cutoff;        /* info needed for busy indicator */
        rt[a].perspective = perspective;
        rt[a].fov = fov;
        rt[a].pos[2] = pos[2];
        rt[a].depth = depth;
        if (I->bkgrd_data) {
          rt[a].bgWidth = I->bkgrd_data->getWidth();
          rt[a].bgHeight = I->bkgrd_data->getHeight();
        }
        rt[a].bkrd_data = I->bkgrd_data ? I->bkgrd_data->bits() : nullptr;
      }

#ifndef _PYMOL_NOPY
      if(n_thread > 1)
        RayTraceSpawn(rt, n_thread);
      else
#endif
        RayTraceThread(rt);

      if(oversample_cutoff) {   /* perform edge oversampling, if requested */
        unsigned int *edging;

        edging = CacheAlloc(I->G, unsigned int, buffer_size, 0, cCache_ray_edging_buffer);

        memcpy(edging, image, buffer_size * sizeof(unsigned int));

        for(a = 0; a < n_thread; a++) {
          rt[a].edging = edging;
        }

#ifndef _PYMOL_NOPY
        if(n_thread > 1)
          RayTraceSpawn(rt, n_thread);
        else
#endif
          RayTraceThread(rt);

        CacheFreeP(I->G, edging, 0, cCache_ray_edging_buffer, false);
      }
      FreeP(rt);
    }
  }

  if(ok && depth && ray_trace_mode) {
    float *delta = pymol::malloc<float>(3 * width * height);
    int x, y;
    ErrChkPtr(I->G, delta);
    if (ok) {
      int xc, yc;
      float d, dzdx, dzdy, *p, *q, dd;
      p = depth;
      q = delta;
      for(y = 0; y < height; y++)
        for(x = 0; x < width; x++) {
          dzdx = 0.0F;
          dzdy = 0.0F;
          xc = 0;
          yc = 0;
          d = *p;
          if(x) {
            dd = d - p[-1];
            dzdx += dd;
            xc++;
          }
          if(x < (width - 1)) {
            dd = p[1] - d;
            if((!xc) || (fabs(dd) > fabs(dzdx)))
              dzdx = dd;
            xc = 1;
          }
          if(y) {
            dd = d - p[-width];
            dzdy += dd;
            yc++;
          }
          if(y < (height - 1)) {
            dd = p[width] - d;
            if((!yc) || (fabs(dd) > fabs(dzdy)))
              dzdy = dd;
            yc = 1;
          }
          p++;
          *(q++) = dzdx;
          *(q++) = dzdy;
          /*
             if(((x == y )||(x==width/2)||(y==height/2))&&((dzdx!=0.0F) || (dzdy!=0.0F))) {
             printf("%5d %5d : %8.3f %8.3f\n",y,x,dzdx,dzdy);
             } */

          *(q++) = sqrt1f(dzdx * dzdx + dzdy * dzdy);
        }
    }

    if (ok){
      int i;
      {
        const float _1 = 1.0F;

        float invFrontMinusBack = _1 / (front - back);
        float inv1minusFogStart = _1;
        int fogFlag = false;
        float fog_start = 0.0F;
        int fogRangeFlag = false;
        float fog = SettingGetGlobal_f(I->G, cSetting_ray_trace_fog);
        if(fog < 0.0F) {
          if(SettingGetGlobal_b(I->G, cSetting_depth_cue)) {
            fog = SettingGetGlobal_f(I->G, cSetting_fog);
          } else
            fog = _0;
        }

        if(fog != _0) {
          if(fog > 1.0F)
            fog = 1.0F;
          fogFlag = true;
          fog_start = SettingGetGlobal_f(I->G, cSetting_ray_trace_fog_start);
          if(fog_start < 0.0F)
            fog_start = SettingGetGlobal_f(I->G, cSetting_fog_start);
          if(fog_start > 1.0F)
            fog_start = 1.0F;
          if(fabs(fog_start) > R_SMALL4) {
            fogRangeFlag = true;
            if(fabs(fog_start - 1.0F) < R_SMALL4)       /* prevent div/0 */
              fogFlag = false;
          }
          inv1minusFogStart = _1 / (_1 - fog_start);
        }

        if(fogFlag) {           /* make sure we have depth values at every potentially drawn pixel */
          float *tmp = pymol::malloc<float>(width * height);
          float dep;
          float *p, *q;
          int cnt;
          for(i = 0; i < 3; i++) {      /* three passes required */
            p = depth;
            q = tmp;
            for(y = 0; y < height; y++)
              for(x = 0; x < width; x++) {
                if(fabs(*p) < R_SMALL4) {
                  dep = 0.0F;
                  cnt = 0;
                  if(x) {
                    if(fabs(p[-1]) > R_SMALL4) {
                      dep += p[-1];
                      cnt++;
                    }
                  }
                  if(x < (width - 1)) {
                    if(fabs(p[1]) > R_SMALL4) {
                      dep += p[1];
                      cnt++;
                    }
                  }
                  if(y) {
                    if(fabs(p[-width]) > R_SMALL4) {
                      dep += p[-width];
                      cnt++;
                    }
                  }
                  if(y < (height - 1)) {
                    if(fabs(p[width]) > R_SMALL4) {
                      dep += p[width];
                      cnt++;
                    }
                  }
                  if(cnt) {
                    dep /= cnt;
                    *q = dep;
                  }
                } else {
                  *q = *p;
                }
                p++;
                q++;
              }
            p = tmp;
            tmp = depth;
            depth = p;
          }
          FreeP(tmp);
        }
        {
          unsigned int *q = image;
          float *p = delta;
          int width3 = width * 3;

          float slope_f = SettingGetGlobal_f(I->G, cSetting_ray_trace_slope_factor);
          float depth_f = SettingGetGlobal_f(I->G, cSetting_ray_trace_depth_factor);
          float disco_f = SettingGetGlobal_f(I->G, cSetting_ray_trace_disco_factor);
          float diff, max_depth;
          float dot, min_dot, max_slope, max_dz, max_pz;
          float dx, dy, dz, px = 0.0F, py = 0.0F, pz = 0.0F, ddx, ddy;
          const float _8 = 0.08F;
          const float _4 = 0.4F;
          const float _25 = 0.25F;
          const float _m25 = -0.25F;
          float disco_f_625 = disco_f * 0.625F;
          float disco_f_5 = disco_f * 0.5F;
          float disco_f_45 = disco_f * 0.45F;

          {
            float gain =
              I->PixelRadius / (SettingGetGlobal_f(I->G, cSetting_ray_trace_gain) *
                                I->Magnified);
            if(antialias)
              gain /= antialias;
            slope_f *= gain;
            depth_f *= gain;
            disco_f *= gain;
          }
          for(y = 0; y < height; y++){
	    float bkrd[3], perc = 0.;
	    if (bkrd_is_gradient) {
	      /* for RayRender, y is from bottom to top */
	      perc = y/(float)height;
	      bkrd[0] = bkrd_bottom[0] + perc * (bkrd_top[0] - bkrd_bottom[0]);
	      bkrd[1] = bkrd_bottom[1] + perc * (bkrd_top[1] - bkrd_bottom[1]);
	      bkrd[2] = bkrd_bottom[2] + perc * (bkrd_top[2] - bkrd_bottom[2]);
	      if(I->BigEndian) {
		background = back_mask |
		  ((0xFF & ((unsigned int) (bkrd[0] * 255 + _p499))) << 24) |
		  ((0xFF & ((unsigned int) (bkrd[1] * 255 + _p499))) << 16) |
		  ((0xFF & ((unsigned int) (bkrd[2] * 255 + _p499))) << 8);
	      } else {
		background = back_mask |
		  ((0xFF & ((unsigned int) (bkrd[2] * 255 + _p499))) << 16) |
		  ((0xFF & ((unsigned int) (bkrd[1] * 255 + _p499))) << 8) |
		  ((0xFF & ((unsigned int) (bkrd[0] * 255 + _p499))));
	      }
	    }
            for(x = 0; x < width; x++) {
              max_slope = 0.0F;
              max_depth = 0.0F;
              min_dot = 1.0F;
              max_dz = 0.0F;
              max_pz = 0.0F;
              dx = p[0];
              dy = p[1];
              dz = p[2];
              for(i = 0; i < 8; i++) {
                switch (i) {
                case 0:
                  if(x) {
                    px = p[-3];
                    py = p[-2];
                    pz = p[-1];
                  }
                  break;
                case 1:
                  if(x < (width - 1)) {
                    px = p[3];
                    py = p[4];
                    pz = p[5];
                  }
                  break;
                case 2:
                  if(y) {
                    px = p[-width3];
                    py = p[-width3 + 1];
                    pz = p[-width3 + 2];
                  }
                  break;
                case 3:
                  if(y < (height - 1)) {
                    px = p[width3];
                    py = p[width3 + 1];
                    pz = p[width3 + 2];
                  }
                  break;
                case 4:
                  if(x && y) {
                    px = p[-width3 - 3];
                    py = p[-width3 - 2];
                    pz = p[-width3 - 1];
                  }
                  break;
                case 5:
                  if(x && (y < (height - 1))) {
                    px = p[width3 - 3];
                    py = p[width3 - 2];
                    pz = p[width3 - 1];
                  }
                  break;
                case 6:
                  if(y && (x < (width - 1))) {
                    px = p[-width3 + 3];
                    py = p[-width3 + 4];
                    pz = p[-width3 + 5];
                  }
                  break;
                case 7:
                  if((y < (height - 1)) && (x < (width - 1))) {
                    px = p[width3 + 3];
                    py = p[width3 + 4];
                    pz = p[width3 + 5];
                  }
                  break;
                }
                ddx = dx - px;
                ddy = dy - py;
                diff = ddx * ddx + ddy * ddy;
                if(max_depth < diff)
                  max_depth = diff;
                if((dz > R_SMALL4) && (pz > R_SMALL4)) {
                  dot = (dx / dz) * (px / pz) + (dy / dz) * (py / pz);
                  if(dot < min_dot) {
                    min_dot = dot;
                    max_dz = dz;
                    max_pz = pz;
                  }
                }
                /*                if(dz>max_dz) max_dz = dz;
                   if(pz>max_pz) max_pz = pz; */
                diff = fabs(dz - pz);
                if(diff > max_slope)
                  max_slope = diff;
              }
              if((max_slope > (slope_f))        /* depth */
                 ||(max_depth > (depth_f))
                 /* slope */
                 /* gradient discontinuities -- could probably use more tuning... */
                 || ((min_dot < _8) && ((max_dz > disco_f) || (max_pz > disco_f)))
                 || ((min_dot < _4) && ((max_dz > disco_f_625) || (max_pz > disco_f_625)))
                 || ((min_dot < _25) && ((max_dz > disco_f_5) || (max_pz > disco_f_5)))
                 || ((min_dot < _m25) && ((max_dz > disco_f_45) && (max_pz > disco_f_45)))
                ) {
                if(fogFlag) {

                  float ffact = depth[q - image] * invFrontMinusBack;
                  float ffact1m;
                  float fc[4];
                  unsigned int cc0, cc1, cc2, cc3;

                  if(fogRangeFlag)
                    ffact = (ffact - fog_start) * inv1minusFogStart;

                  ffact *= fog;

                  if(ffact < _0)
                    ffact = _0;
                  if(ffact > _1)
                    ffact = _1;

                  ffact1m = _1 - ffact;

                  if(orig_opaque_back) {
                    fc[0] =
                      (0xFF & (background >> 24)) * ffact +
                      (0xFF & (trace_word >> 24)) * ffact1m;
                    fc[1] =
                      (0xFF & (background >> 16)) * ffact +
                      (0xFF & (trace_word >> 16)) * ffact1m;
                    fc[2] =
                      (0xFF & (background >> 8)) * ffact +
                      (0xFF & (trace_word >> 8)) * ffact1m;
                    fc[3] =
                      (0xFF & (background)) * ffact + (0xFF & (trace_word)) * ffact1m;
                  } else {      /* if non-opaque background, then use alpha to blend */
                    fc[1] = (0xFF & (trace_word >> 16));
                    fc[2] = (0xFF & (trace_word >> 8));
                    if(I->BigEndian) {
                      fc[0] = (0xFF & (trace_word >> 24));
                      fc[3] =
                        (0xFF & (background)) * ffact + (0xFF & (trace_word)) * ffact1m;
                    } else {
                      fc[0] =
                        (0xFF & (background >> 24)) * ffact +
                        (0xFF & (trace_word >> 24)) * ffact1m;
                      fc[3] = (0xFF & (trace_word));
                    }
                  }
                  cc0 = (uint) (fc[0]);
                  cc1 = (uint) (fc[1]);
                  cc2 = (uint) (fc[2]);
                  cc3 = (uint) (fc[3]);

                  if(cc0 > 255)
                    cc0 = 255;
                  if(cc1 > 255)
                    cc1 = 255;
                  if(cc2 > 255)
                    cc2 = 255;
                  if(cc3 > 255)
                    cc3 = 255;

                  *q = (cc0 << 24) | (cc1 << 16) | (cc2 << 8) | cc3;
                } else {
                  *q = trace_word;
                }
              } else if(ray_trace_mode == 2) {      /* only draw edge */
                *q = background;
              } else if(ray_trace_mode == 3) {      /* quantize */
                *q = (*q & 0xC0C0C0C0);
                *q = *q | ((*q) >> 2) | ((*q) >> 4) | ((*q) >> 6);
              }
              p += 3;
              q++;
            }
	  }
        }
      }
    }
    FreeP(delta);
  }

  if(ok && antialias > 1) {
    /* now spawn threads as needed */
    CRayAntiThreadInfo *rt = pymol::calloc<CRayAntiThreadInfo>(n_thread);

    for(a = 0; a < n_thread; a++) {
      rt[a].width = width;
      rt[a].height = height;
      rt[a].image = image;
      rt[a].image_copy = image_copy;
      rt[a].phase = a;
      rt[a].mag = mag;          /* fold magnification */
      rt[a].n_thread = n_thread;
      rt[a].ray = I;
    }

#ifndef _PYMOL_NOPY
    if(n_thread > 1)
      RayAntiSpawn(rt, n_thread);
    else
#endif
      RayAntiThread(rt);
    FreeP(rt);
    CacheFreeP(I->G, image, 0, cCache_ray_antialias_buffer, false);
    image = image_copy;
  }

  PRINTFD(I->G, FB_Ray)
    " RayRender: n_hit %d\n", n_hit ENDFD;
#ifdef PROFILE_BASIS

  printf
    ("int n_cells = %d;\nint n_prims = %d;\nint n_triangles = %8.3f;\nint n_spheres = %8.3f;\nint n_cylinders = %8.3f;\nint n_sausages = %8.3f;\nint n_skipped = %8.3f;\n",
     n_cells, n_prims, n_triangles / ((float) n_cells), n_spheres / ((float) n_cells),
     n_cylinders / ((float) n_cells), n_sausages / ((float) n_cells),
     n_skipped / ((float) n_cells));
#endif

  if (ok){
    /* EXPERIMENTAL RAY-VOLUME CODE */
    volume = SettingGetGlobal_b(I->G, cSetting_ray_volume);
    
    if (volume) {
      for(y = 0; y < height; y++) {
	for(x = 0; x < width; x++) {
	  float dd = depth[x+width*y];
	  if (dd == 0.0) dd = -back;
	  depth[x+width*y] = -dd/(back-front) + 0.1;
	}
      }
      if (rayDepthPixels)
	FreeP(rayDepthPixels);
      rayDepthPixels = depth;
      rayWidth = width;
      rayHeight = height;
      rayVolume = 3;
    } else 
      FreeP(depth);
  }
  I->bkgrd_data = nullptr;
}


void RayRenderColorTable(CRay * I, int width, int height, int *image)
{
  int x, y;
  unsigned int r = 0, g = 0, b = 0;
  unsigned int *pixel, mask, *p;

  if(I->BigEndian)
    mask = 0x000000FF;
  else
    mask = 0xFF000000;

  p = (unsigned int *) image;
  for(x = 0; x < width; x++)
    for(y = 0; y < height; y++)
      *(p++) = mask;

  if((width >= 512) && (height >= 512)) {

    for(y = 0; y < 512; y++)
      for(x = 0; x < 512; x++) {
        pixel = (unsigned int *) (image + ((width) * y) + x);
        if(I->BigEndian) {
          *(pixel) = mask | (r << 24) | (g << 16) | (b << 8);
        } else {
          *(pixel) = mask | (b << 16) | (g << 8) | r;
        }
        b = b + 4;
        if(!(0xFF & b)) {
          b = 0;
          g = g + 4;
          if(!(0xFF & g)) {
            g = 0;
            r = r + 4;
          }
        }
      }
  }
}


/*========================================================================*/
void CRay::wobble(int mode, const float *v)
{
  CRay * I = this;
  I->Wobble = mode;
  if(v)
    copy3f(v, I->WobbleParam);
}


/*========================================================================*/
void CRay::transparentf(float v)
{
  CRay * I = this;
  if(v > 1.0F)
    v = 1.0F;
  if(v < 0.0F)
    v = 0.0F;
  I->Trans = v;
}

void CRay::interiorColor3fv(const float *v, int passive)
{
  CRay * I = this;
  I->IntColor[0] = (*v++);
  I->IntColor[1] = (*v++);
  I->IntColor[2] = (*v++);
  if(!passive)
    I->CheckInterior = true;
}


/*========================================================================*/
void CRay::color3fv(const float *v)
{
  CRay * I = this;
  I->CurColor[0] = (*v++);
  I->CurColor[1] = (*v++);
  I->CurColor[2] = (*v++);
}


/*========================================================================*/
int CRay::sphere3fv(const float *v, float r)
{
  CRay * I = this;
  CPrimitive *p;
  int ok = true;
  float *vv;

  VLACacheCheck(I->G, I->Primitive, CPrimitive, I->NPrimitive, 0, cCache_ray_primitive);
  CHECKOK(ok, I->Primitive);
  if (!ok)
    return false;
  p = I->Primitive + I->NPrimitive;

  p->type = cPrimSphere;
  p->r1 = r;
  p->trans = I->Trans;
  p->wobble = I->Wobble;
  p->ramped = (I->CurColor[0] < 0.0F);
  p->no_lighting = 0;

  I->PrimSize += 2 * r;
  I->PrimSizeCnt++;

  /* 
     copy3f(I->WobbleParam,p->wobble_param); */
  vv = p->v1;
  (*vv++) = (*v++);
  (*vv++) = (*v++);
  (*vv++) = (*v++);

  vv = p->c1;
  v = I->CurColor;
  (*vv++) = (*v++);
  (*vv++) = (*v++);
  (*vv++) = (*v++);

  vv = p->ic;
  v = I->IntColor;
  (*vv++) = (*v++);
  (*vv++) = (*v++);
  (*vv++) = (*v++);

  if(I->TTTFlag) {
    p->r1 *= length3f(glm::value_ptr(I->TTT));
    transformTTT44f3f(glm::value_ptr(I->TTT), p->v1, p->v1);
  }

  RayApplyContextToVertex(I, p->v1);

  I->NPrimitive++;
  return true;
}

void RayGetScaledAxes(CRay * I, float *xn, float *yn)
{
  float *v;
  float vt[3];
  float xn0[3] = { 1.0F, 0.0F, 0.0F };
  float yn0[3] = { 0.0F, 1.0F, 0.0F };
  float v_scale;

  v = TextGetPos(I->G);

  if(I->TTTFlag) {
    transformTTT44f3f(glm::value_ptr(I->TTT), v, vt);
  } else {
    copy3f(v, vt);
  }

  v_scale = RayGetScreenVertexScale(I, vt) / I->Sampling;

  RayApplyMatrixInverse33(1, (float3 *) xn0, glm::value_ptr(I->Rotation), (float3 *) xn0);
  RayApplyMatrixInverse33(1, (float3 *) yn0, glm::value_ptr(I->Rotation), (float3 *) yn0);

  scale3f(xn0, v_scale, xn);
  scale3f(yn0, v_scale, yn);
}


/*========================================================================*/
int CRay::character(int char_id)
{
  CRay * I = this;
  CPrimitive *p;
  float *v;
  float vt[3];
  float *vv;
  float width, height;
  float v_scale;
  int ok = true;

  v = TextGetPos(I->G);
  VLACacheCheck(I->G, I->Primitive, CPrimitive, I->NPrimitive + 1, 0,
                cCache_ray_primitive);
  CHECKOK(ok, I->Primitive);
  if (!ok)
    return false;
  p = I->Primitive + I->NPrimitive;

  p->type = cPrimCharacter;
  p->trans = I->Trans;
  p->char_id = char_id;
  p->wobble = I->Wobble;
  p->ramped = 0;
  p->no_lighting = 0;
  /*
     copy3f(I->WobbleParam,p->wobble_param); */

  vv = p->v1;
  (*vv++) = v[0];
  (*vv++) = v[1];
  (*vv++) = v[2];

  if(I->TTTFlag) {
    transformTTT44f3f(glm::value_ptr(I->TTT), p->v1, p->v1);
  }
  /* what's the width of 1 screen window pixel at this point in space? */

  v_scale = RayGetScreenVertexScale(I, p->v1) / I->Sampling;

  RayApplyContextToVertex(I, p->v1);

  {
    float xn[3] = { 1.0F, 0.0F, 0.0F };
    float yn[3] = { 0.0F, 1.0F, 0.0F };
    float zn[3] = { 0.0F, 0.0F, 1.0F };
    float sc[3];
    float scale;
    float xorig, yorig, advance;
    int width_i, height_i;
    CPrimitive *pp = p + 1;

    RayApplyMatrixInverse33(1, (float3 *) xn, glm::value_ptr(I->Rotation), (float3 *) xn);
    RayApplyMatrixInverse33(1, (float3 *) yn, glm::value_ptr(I->Rotation), (float3 *) yn);
    RayApplyMatrixInverse33(1, (float3 *) zn, glm::value_ptr(I->Rotation), (float3 *) zn);

    CharacterGetGeometry(I->G, char_id, &width_i, &height_i, &xorig, &yorig, &advance);
    width = (float) width_i;
    height = (float) height_i;

    scale = v_scale * advance;
    scale3f(xn, scale, vt);     /* advance raster position in 3-space */
    add3f(v, vt, vt);
    TextSetPos(I->G, vt);

    /* position the pixmap relative to raster position */

    /*    scale = ((-xorig)-0.5F)*I->PixelRadius; */
    scale = ((-xorig) - 0.0F) * v_scale;
    scale3f(xn, scale, sc);
    add3f(sc, p->v1, p->v1);

    scale = ((-yorig) - 0.0F) * v_scale;
    scale3f(yn, scale, sc);
    add3f(sc, p->v1, p->v1);

    scale = v_scale * width;
    scale3f(xn, scale, xn);
    scale = v_scale * height;
    scale3f(yn, scale, yn);

    copy3f(zn, p->n0);
    copy3f(zn, p->n1);
    copy3f(zn, p->n2);
    copy3f(zn, p->n3);

    *(pp) = (*p);

    /* define coordinates of first triangle */

    add3f(p->v1, xn, p->v2);
    add3f(p->v1, yn, p->v3);

    I->PrimSize +=
      2 * (diff3f(p->v1, p->v2) + diff3f(p->v1, p->v3) + diff3f(p->v2, p->v3));
    I->PrimSizeCnt += 6;

    /* encode characters coordinates in the colors  */

    zero3f(p->c1);
    set3f(p->c2, width, 0.0F, 0.0F);
    set3f(p->c3, 0.0F, height, 0.0F);

    /* define coordinates of second triangle */

    add3f(yn, xn, pp->v1);
    add3f(p->v1, pp->v1, pp->v1);
    add3f(p->v1, yn, pp->v2);
    add3f(p->v1, xn, pp->v3);

    {
      float *v, *vv;
      vv = p->ic;
      v = I->IntColor;
      (*vv++) = (*v++);
      (*vv++) = (*v++);
      (*vv++) = (*v++);
      vv = pp->ic;
      v = I->IntColor;
      (*vv++) = (*v++);
      (*vv++) = (*v++);
      (*vv++) = (*v++);
    }

    /* encode integral character coordinates into the vertex colors  */

    set3f(pp->c1, width, height, 0.0F);
    set3f(pp->c2, 0.0F, height, 0.0F);
    set3f(pp->c3, width, 0.0F, 0.0F);

  }

  I->NPrimitive += 2;
  return true;
}


/*========================================================================*/
int CRay::cylinder3fv(const cgo::draw::cylinder &cyl){
  return cylinder3fv(cyl.vertex1, cyl.vertex2, cyl.radius, cyl.color1, cyl.color2, 1.0f - Trans, 1.0f - Trans);
}
int CRay::cylinder3fv(const cgo::draw::cylinder &cyl, const float alpha1, const float alpha2){
  return cylinder3fv(cyl.vertex1, cyl.vertex2, cyl.radius, cyl.color1, cyl.color2, alpha1, alpha2);
}

int CRay::cylinder3fv(const float *v1, const float *v2, float r, const float *c1, const float *c2,
                      const float alpha1, const float alpha2)
{
  CRay * I = this;
  CPrimitive *p;
  int ok = true;
  float *vv;

  VLACacheCheck(I->G, I->Primitive, CPrimitive, I->NPrimitive, 0, cCache_ray_primitive);
  CHECKOK(ok, I->Primitive);
  if (!ok)
    return false;
  p = I->Primitive + I->NPrimitive;

  p->type = cPrimCylinder;
  p->r1 = r;
  p->cap1 = cCylCapFlat;
  p->cap2 = cCylCapFlat;
  p->wobble = I->Wobble;
  p->ramped = ((c1[0] < 0.0F) || (c2[0] < 0.0F));
  p->no_lighting = 0;
  /* 
     copy3f(I->WobbleParam,p->wobble_param); */

  vv = p->v1;
  (*vv++) = (*v1++);
  (*vv++) = (*v1++);
  (*vv++) = (*v1++);
  vv = p->v2;
  (*vv++) = (*v2++);
  (*vv++) = (*v2++);
  (*vv++) = (*v2++);

  I->PrimSize += diff3f(p->v1, p->v2) + 2 * r;
  I->PrimSizeCnt++;

  if(I->TTTFlag) {
    p->r1 *= length3f(glm::value_ptr(I->TTT));
    transformTTT44f3f(glm::value_ptr(I->TTT), p->v1, p->v1);
    transformTTT44f3f(glm::value_ptr(I->TTT), p->v2, p->v2);
  }

  RayApplyContextToVertex(I, p->v1);
  RayApplyContextToVertex(I, p->v2);

  vv = p->c1;
  (*vv++) = (*c1++);
  (*vv++) = (*c1++);
  (*vv++) = (*c1++);
  vv = p->c2;
  (*vv++) = (*c2++);
  (*vv++) = (*c2++);
  (*vv++) = (*c2++);

  // FIXME: alpha1 is not used
  p->trans = 1.0 - alpha2;
  {
    float *v;
    vv = p->ic;
    v = I->IntColor;
    (*vv++) = (*v++);
    (*vv++) = (*v++);
    (*vv++) = (*v++);
  }

  I->NPrimitive++;
  return true;
}


/*========================================================================*/
int CRay::customCylinder3fv(const cgo::draw::custom_cylinder& cyl, const float alpha1, const float alpha2){
  return customCylinder3fv(cyl.vertex1, cyl.vertex2, cyl.radius, cyl.color1,
                        cyl.color2, cyl.get_cap1(), cyl.get_cap2(), alpha1, alpha2);
}

int CRay::customCylinder3fv(const cgo::draw::custom_cylinder& cyl){
  return customCylinder3fv(cyl.vertex1, cyl.vertex2, cyl.radius, cyl.color1,
                        cyl.color2, cyl.get_cap1(), cyl.get_cap2(), 1.f - Trans, 1.f - Trans);
}

int CRay::customCylinderAlpha3fv(const cgo::draw::custom_cylinder_alpha& cyl){
  return customCylinder3fv(cyl.vertex1, cyl.vertex2, cyl.radius, cyl.color1,
                        cyl.color2, cyl.get_cap1(), cyl.get_cap2(), cyl.color1[3], cyl.color2[3]);
}

int CRay::customCylinder3fv(const float* v1, const float* v2, float r,
    const float* c1, const float* c2, const cCylCap cap1, const cCylCap cap2)
{
  return customCylinder3fv(
      v1, v2, r, c1, c2, cap1, cap2, 1.f - Trans, 1.f - Trans);
}

int CRay::customCylinder3fv(const float *v1, const float *v2, float r,
                            const float *c1, const float *c2,
                            const cCylCap cap1,
                            const cCylCap cap2,
                            const float alpha1, const float alpha2)
{
  CRay * I = this;
  CPrimitive *p;
  int ok = true;
  float *vv;

  VLACacheCheck(I->G, I->Primitive, CPrimitive, I->NPrimitive, 0, cCache_ray_primitive);
  CHECKOK(ok, I->Primitive);
  if (!ok)
    return false;
  p = I->Primitive + I->NPrimitive;

  p->type = cPrimCylinder;
  p->r1 = r;
  p->cap1 = cap1;
  p->cap2 = cap2;
  p->wobble = I->Wobble;
  p->ramped = ((c1[0] < 0.0F) || (c2[0] < 0.0F));
  p->no_lighting = 0;
  /*
     copy3f(I->WobbleParam,p->wobble_param); */

  vv = p->v1;
  (*vv++) = (*v1++);
  (*vv++) = (*v1++);
  (*vv++) = (*v1++);
  vv = p->v2;
  (*vv++) = (*v2++);
  (*vv++) = (*v2++);
  (*vv++) = (*v2++);

  I->PrimSize += diff3f(p->v1, p->v2) + 2 * r;
  I->PrimSizeCnt++;

  if(I->TTTFlag) {
    p->r1 *= length3f(glm::value_ptr(I->TTT));
    transformTTT44f3f(glm::value_ptr(I->TTT), p->v1, p->v1);
    transformTTT44f3f(glm::value_ptr(I->TTT), p->v2, p->v2);
  }

  RayApplyContextToVertex(I, p->v1);
  RayApplyContextToVertex(I, p->v2);

  vv = p->c1;
  (*vv++) = (*c1++);
  (*vv++) = (*c1++);
  (*vv++) = (*c1++);
  vv = p->c2;
  (*vv++) = (*c2++);
  (*vv++) = (*c2++);
  (*vv++) = (*c2++);
  vv = p->ic;

  // FIXME: alpha1 is not used
  p->trans = 1.0f - alpha2;

  {
    float *v;
    vv = p->ic;
    v = I->IntColor;
    (*vv++) = (*v++);
    (*vv++) = (*v++);
    (*vv++) = (*v++);
  }

  I->NPrimitive++;
  return true;
}

int CRay::cone3fv(const float *v1, const float *v2, float r1, float r2,
	       const float *c1, const float *c2, cCylCap cap1, cCylCap cap2)
{
  CRay * I = this;
  CPrimitive *p;
  float r_max = (r1 > r2) ? r1 : r2;
  float *vv;
  int ok = true;

  if(r2 > r1) {                 /* make sure r1 is always larger */
    std::swap(r1, r2);
    std::swap(c1, c2);
    std::swap(v1, v2);
    std::swap(cap1, cap2);
  }

  VLACacheCheck(I->G, I->Primitive, CPrimitive, I->NPrimitive, 0, cCache_ray_primitive);
  CHECKOK(ok, I->Primitive);
  if (!ok)
    return false;
  p = I->Primitive + I->NPrimitive;

  p->type = cPrimCone;
  p->r1 = r1;
  p->r2 = r2;
  p->trans = I->Trans;
  p->cap1 = cap1;

  if(cap2 >= cCylCapFlat)
    cap2 = cCylCapFlat;
  if(cap1 >= cCylCapFlat)
    cap1 = cCylCapFlat;

  p->cap2 = cap2;
  p->wobble = I->Wobble;
  p->ramped = ((c1[0] < 0.0F) || (c2[0] < 0.0F));
  p->no_lighting = 0;
  /*
     copy3f(I->WobbleParam,p->wobble_param); */

  vv = p->v1;
  (*vv++) = (*v1++);
  (*vv++) = (*v1++);
  (*vv++) = (*v1++);
  vv = p->v2;
  (*vv++) = (*v2++);
  (*vv++) = (*v2++);
  (*vv++) = (*v2++);

  I->PrimSize += diff3f(p->v1, p->v2) + 2 * r_max;
  I->PrimSizeCnt++;

  if(I->TTTFlag) {
    transformTTT44f3f(glm::value_ptr(I->TTT), p->v1, p->v1);
    transformTTT44f3f(glm::value_ptr(I->TTT), p->v2, p->v2);
  }

  RayApplyContextToVertex(I, p->v1);
  RayApplyContextToVertex(I, p->v2);

  vv = p->c1;
  (*vv++) = (*c1++);
  (*vv++) = (*c1++);
  (*vv++) = (*c1++);
  vv = p->c2;
  (*vv++) = (*c2++);
  (*vv++) = (*c2++);
  (*vv++) = (*c2++);
  vv = p->ic;
  {
    float *v;
    vv = p->ic;
    v = I->IntColor;
    (*vv++) = (*v++);
    (*vv++) = (*v++);
    (*vv++) = (*v++);
  }

  I->NPrimitive++;
  return true;
}


/*========================================================================*/
int CRay::sausage3fv(const float *v1, const float *v2, float r, const float *c1, const float *c2)
{
  CRay * I = this;
  CPrimitive *p;
  int ok = true;
  float *vv;

  VLACacheCheck(I->G, I->Primitive, CPrimitive, I->NPrimitive, 0, cCache_ray_primitive);
  CHECKOK(ok, I->Primitive);
  if (!ok)
    return false;
  p = I->Primitive + I->NPrimitive;

  p->type = cPrimSausage;
  p->r1 = r;
  p->trans = I->Trans;
  p->wobble = I->Wobble;
  p->ramped = ((c1[0] < 0.0F) || (c2[0] < 0.0F));
  p->no_lighting = 0;
  /*  
     copy3f(I->WobbleParam,p->wobble_param); */

  vv = p->v1;
  (*vv++) = (*v1++);
  (*vv++) = (*v1++);
  (*vv++) = (*v1++);
  vv = p->v2;
  (*vv++) = (*v2++);
  (*vv++) = (*v2++);
  (*vv++) = (*v2++);

  I->PrimSize += diff3f(p->v1, p->v2) + 2 * r;
  I->PrimSizeCnt++;

  if(I->TTTFlag) {
    p->r1 *= length3f(glm::value_ptr(I->TTT));
    transformTTT44f3f(glm::value_ptr(I->TTT), p->v1, p->v1);
    transformTTT44f3f(glm::value_ptr(I->TTT), p->v2, p->v2);
  }

  RayApplyContextToVertex(I, p->v1);
  RayApplyContextToVertex(I, p->v2);

  vv = p->c1;
  (*vv++) = (*c1++);
  (*vv++) = (*c1++);
  (*vv++) = (*c1++);
  vv = p->c2;
  (*vv++) = (*c2++);
  (*vv++) = (*c2++);
  (*vv++) = (*c2++);
  vv = p->ic;

  {
    float *v = I->IntColor;
    vv = p->ic;
    (*vv++) = (*v++);
    (*vv++) = (*v++);
    (*vv++) = (*v++);
  }

  I->NPrimitive++;
  return true;
}

int CRay::ellipsoid3fv(const float *v, float r, const float *n1, const float *n2, const float *n3)
{
  CRay * I = this;
  CPrimitive *p;
  int ok = true;
  float *vv;

  VLACacheCheck(I->G, I->Primitive, CPrimitive, I->NPrimitive, 0, cCache_ray_primitive);
  CHECKOK(ok, I->Primitive);
  if (!ok)
    return false;
  p = I->Primitive + I->NPrimitive;

  p->type = cPrimEllipsoid;
  p->r1 = r;                    /* maximum extent */
  p->trans = I->Trans;
  p->wobble = I->Wobble;
  p->ramped = (I->CurColor[0] < 0.0F);
  p->no_lighting = 0;

  I->PrimSize += 2 * r;
  I->PrimSizeCnt++;

  vv = p->n0;                   /* storing lengths of the direction vectors in n0 */

  (*vv++) = length3f(n1);
  (*vv++) = length3f(n2);
  (*vv++) = length3f(n3);

  /* normalize the ellipsoid axes */

  vv = p->n1;
  if(p->n0[0] > R_SMALL8) {
    float factor;
    factor = 1.0F / p->n0[0];
    (*vv++) = (*n1++) * factor;
    (*vv++) = (*n1++) * factor;
    (*vv++) = (*n1++) * factor;
  } else {
    (*vv++) = 0.0F;
    (*vv++) = 0.0F;
    (*vv++) = 0.0F;
  }

  vv = p->n2;
  if(p->n0[1] > R_SMALL8) {
    float factor;
    factor = 1.0F / p->n0[1];
    (*vv++) = (*n2++) * factor;
    (*vv++) = (*n2++) * factor;
    (*vv++) = (*n2++) * factor;
  } else {
    (*vv++) = 0.0F;
    (*vv++) = 0.0F;
    (*vv++) = 0.0F;
  }

  vv = p->n3;
  if(p->n0[2] > R_SMALL8) {
    float factor;
    factor = 1.0F / p->n0[2];
    (*vv++) = (*n3++) * factor;
    (*vv++) = (*n3++) * factor;
    (*vv++) = (*n3++) * factor;
  } else {
    (*vv++) = 0.0F;
    (*vv++) = 0.0F;
    (*vv++) = 0.0F;
  }

  vv = p->v1;
  (*vv++) = (*v++);
  (*vv++) = (*v++);
  (*vv++) = (*v++);

  vv = p->c1;
  v = I->CurColor;
  (*vv++) = (*v++);
  (*vv++) = (*v++);
  (*vv++) = (*v++);

  vv = p->ic;
  v = I->IntColor;
  (*vv++) = (*v++);
  (*vv++) = (*v++);
  (*vv++) = (*v++);

  if(I->TTTFlag) {
    p->r1 *= length3f(glm::value_ptr(I->TTT));
    transformTTT44f3f(glm::value_ptr(I->TTT), p->v1, p->v1);
    transform_normalTTT44f3f(glm::value_ptr(I->TTT), p->n1, p->n1);
    transform_normalTTT44f3f(glm::value_ptr(I->TTT), p->n2, p->n2);
    transform_normalTTT44f3f(glm::value_ptr(I->TTT), p->n3, p->n3);
  }

  RayApplyContextToVertex(I, p->v1);
  RayApplyContextToNormal(I, p->n1);
  RayApplyContextToNormal(I, p->n2);
  RayApplyContextToNormal(I, p->n3);

  I->NPrimitive++;
  return true;
}

int CRay::setLastToNoLighting(char no_lighting)
{
  CRay * I = this;
  CPrimitive *p;
  if (!I->NPrimitive)
    return false;
  p = I->Primitive + I->NPrimitive - 1;
  p->no_lighting = no_lighting;
  return true;
}


/*========================================================================*/
int CRay::triangle3fv(
		   const float *v1, const float *v2, const float *v3,
		   const float *n1, const float *n2, const float *n3, const float *c1, const float *c2, const float *c3)
{
  CRay * I = this;
  CPrimitive *p;
  int ok = true;
  float *vv;
  float n0[3] = { 0.f, 0.f, 1.f }, nx[3], s1[3], s2[3], s3[3];
  float l1, l2, l3;
  short normals_exist = n1 && n2 && n3;

  /*  dump3f(v1," v1");
     dump3f(v2," v2");
     dump3f(v3," v3");
     dump3f(n1," n1");
     dump3f(n2," n2");
     dump3f(n3," n3");
     dump3f(c1," c1");
     dump3f(c2," c2");
     dump3f(c3," c3"); */
  VLACacheCheck(I->G, I->Primitive, CPrimitive, I->NPrimitive, 0, cCache_ray_primitive);
  CHECKOK(ok, I->Primitive);
  if (!ok)
    return false;
  p = I->Primitive + I->NPrimitive;

  p->type = cPrimTriangle;
  p->trans = I->Trans;
  p->tr[0] = I->Trans;
  p->tr[1] = I->Trans;
  p->tr[2] = I->Trans;
  p->wobble = I->Wobble;
  p->ramped = ((c1[0] < 0.0F) || (c2[0] < 0.0F) || (c3[0] < 0.0F));
  p->no_lighting = 0;
  /*
     copy3f(I->WobbleParam,p->wobble_param); */

  /* determine exact triangle normal */
  if (normals_exist){
    add3f(n1, n2, nx);
    add3f(n3, nx, nx);
  }
  subtract3f(v1, v2, s1);
  subtract3f(v3, v2, s2);
  subtract3f(v1, v3, s3);
  cross_product3f(s1, s2, n0);
  if (normals_exist){
    if((fabs(n0[0]) < RAY_SMALL) && (fabs(n0[1]) < RAY_SMALL) && (fabs(n0[2]) < RAY_SMALL)) {
      copy3f(nx, n0);
    } /* fall-back */
    else if(dot_product3f(n0, nx) < 0)
      invert3f(n0);
  }
  normalize3f(n0);

  vv = p->n0;
  (*vv++) = n0[0];
  (*vv++) = n0[1];
  (*vv++) = n0[2];

  /* determine maximum distance from vertex to point */
  l1 = (float) length3f(s1);
  l2 = (float) length3f(s2);
  l3 = (float) length3f(s3);
  if(l2 > l1) {
    if(l3 > l2)
      l1 = l3;
    else
      l1 = l2;
  }
  /* store cutoff distance */

  p->r1 = l1 * 0.6F;

  /*  if(l1>20) {
     printf("%8.3f\n",l1);
     printf("%8.3f %8.3f %8.3f\n",s1[0],s1[1],s1[2]);
     printf("%8.3f %8.3f %8.3f\n",s2[0],s2[1],s2[2]);
     printf("%8.3f %8.3f %8.3f\n",s3[0],s3[1],s3[2]);
     } */

  vv = p->v1;
  (*vv++) = (*v1++);
  (*vv++) = (*v1++);
  (*vv++) = (*v1++);
  vv = p->v2;
  (*vv++) = (*v2++);
  (*vv++) = (*v2++);
  (*vv++) = (*v2++);
  vv = p->v3;
  (*vv++) = (*v3++);
  (*vv++) = (*v3++);
  (*vv++) = (*v3++);

  I->PrimSize += diff3f(p->v1, p->v2) + diff3f(p->v1, p->v3) + diff3f(p->v2, p->v3);
  I->PrimSizeCnt += 3;

  vv = p->c1;
  (*vv++) = (*c1++);
  (*vv++) = (*c1++);
  (*vv++) = (*c1++);
  vv = p->c2;
  (*vv++) = (*c2++);
  (*vv++) = (*c2++);
  (*vv++) = (*c2++);
  vv = p->c3;
  (*vv++) = (*c3++);
  (*vv++) = (*c3++);
  (*vv++) = (*c3++);

  {
    float *v = I->IntColor;
    vv = p->ic;
    (*vv++) = (*v++);
    (*vv++) = (*v++);
    (*vv++) = (*v++);
  }

  if (normals_exist){
    vv = p->n1;
    (*vv++) = (*n1++);
    (*vv++) = (*n1++);
    (*vv++) = (*n1++);
    vv = p->n2;
    (*vv++) = (*n2++);
    (*vv++) = (*n2++);
    (*vv++) = (*n2++);
    vv = p->n3;
    (*vv++) = (*n3++);
    (*vv++) = (*n3++);
    (*vv++) = (*n3++);
  } else {
    vv = p->n1;
    (*vv++) = n0[0];
    (*vv++) = n0[1];
    (*vv++) = n0[2];
    vv = p->n2;
    (*vv++) = n0[0];
    (*vv++) = n0[1];
    (*vv++) = n0[2];
    vv = p->n3;
    (*vv++) = n0[0];
    (*vv++) = n0[1];
    (*vv++) = n0[2];
  }

  if(I->TTTFlag) {
    transformTTT44f3f(glm::value_ptr(I->TTT), p->v1, p->v1);
    transformTTT44f3f(glm::value_ptr(I->TTT), p->v2, p->v2);
    transformTTT44f3f(glm::value_ptr(I->TTT), p->v3, p->v3);
    transform_normalTTT44f3f(glm::value_ptr(I->TTT), p->n0, p->n0);
    transform_normalTTT44f3f(glm::value_ptr(I->TTT), p->n1, p->n1);
    transform_normalTTT44f3f(glm::value_ptr(I->TTT), p->n2, p->n2);
    transform_normalTTT44f3f(glm::value_ptr(I->TTT), p->n3, p->n3);
  }

  RayApplyContextToVertex(I, p->v1);
  RayApplyContextToVertex(I, p->v2);
  RayApplyContextToVertex(I, p->v3);
  RayApplyContextToNormal(I, p->n0);
  RayApplyContextToNormal(I, p->n1);
  RayApplyContextToNormal(I, p->n2);
  RayApplyContextToNormal(I, p->n3);

  I->NPrimitive++;
  return true;
}

int CRay::triangleTrans3fv(
			const float *v1, const float *v2, const float *v3,
			const float *n1, const float *n2, const float *n3,
			const float *c1, const float *c2, const float *c3, float t1, float t2, float t3)
{
  CRay * I = this;
  CPrimitive *p;
  int ok = true;
  ok = I->triangle3fv(v1, v2, v3, n1, n2, n3, c1, c2, c3);
  if (!ok)
    return false;
  p = I->Primitive + I->NPrimitive - 1;

  p->tr[0] = t1;
  p->tr[1] = t2;
  p->tr[2] = t3;
  p->trans = (t1 + t2 + t3) / 3.0F;
  return true;
}


/*========================================================================*/
CRay *RayNew(PyMOLGlobals * G, int antialias)
{
  unsigned int test;
  unsigned char *testPtr;
  int a;

  auto I = new CRay();
  I->G = G;
  test = 0xFF000000;
  testPtr = (unsigned char *) &test;
  I->BigEndian = (*testPtr) & 0x01;
  I->Trans = 0.0F;
  I->Wobble = 0;
  I->TTTFlag = false;
  zero3f(I->WobbleParam);
  PRINTFB(I->G, FB_Ray, FB_Blather)
    " RayNew: BigEndian = %d\n", I->BigEndian ENDFB(I->G);

  I->Basis = CacheAlloc(I->G, CBasis, 12, 0, cCache_ray_basis);
  BasisInit(I->G, I->Basis, 0);
  BasisInit(I->G, I->Basis + 1, 1);
  I->NBasis = 2;
  I->Primitive = nullptr;
  I->NPrimitive = 0;
  I->CheckInterior = false;
  if(antialias < 0)
    antialias = SettingGetGlobal_i(I->G, cSetting_antialias);
  I->Sampling = antialias;
  if(I->Sampling < 2)           /* always supersample fonts by at least 2X */
    I->Sampling = 2;
  for(a = 0; a < 256; a++) {
    I->Random[a] = (float) ((rand() / (1.0 + RAND_MAX)) - 0.5);
  }

  I->Wobble = SettingGet_i(I->G, nullptr, nullptr, cSetting_ray_texture);
  {
    auto v = SettingGet<const float*>(I->G, cSetting_ray_texture_settings);
    int color = SettingGetGlobal_color(I->G, cSetting_ray_interior_color);
    copy3f(v, I->WobbleParam);
    v = ColorGet(I->G, color);
    copy3f(v, I->IntColor);
  }

  return (I);
}


/*========================================================================*/

/* BEGIN PROPRIETARY CODE SEGMENT (see disclaimer in "os_proprietary.h") */
#ifdef PYMOL_EVAL
#include "RayEvalMessage.h"
#endif

/* END PROPRIETARY CODE SEGMENT */

void RayPrepare(CRay * I, float v0, float v1, float v2,
                float v3, float v4, float v5,
                float fov, glm::vec3 pos,
                float *mat, const glm::mat4& rotMat, float aspRat,
                int width, int height, float pixel_scale, int ortho,
                float pixel_ratio, float front_back_ratio, float magnified)

/*prepare for vertex calls */
{
  int a;
  if(!I->Primitive)
    I->Primitive = VLACacheAlloc(I->G, CPrimitive, 10000, 3, cCache_ray_primitive);
  I->Volume[0] = v0;
  I->Volume[1] = v1;
  I->Volume[2] = v2;
  I->Volume[3] = v3;
  I->Volume[4] = v4;
  I->Volume[5] = v5;
  I->Range[0] = I->Volume[1] - I->Volume[0];
  I->Range[1] = I->Volume[3] - I->Volume[2];
  I->Range[2] = I->Volume[5] - I->Volume[4];
  I->AspRatio = aspRat;
  I->Width = width;
  I->Height = height;
  CharacterSetRetention(I->G, true);

  if(mat)
    for(a = 0; a < 16; a++)
      I->ModelView[a] = mat[a];
  else {
    identity44f(I->ModelView);
  }
  identity44f(I->ProMatrix);
  if (ortho){
    I->ProMatrix[0] = 2.f/I->Range[0];
    I->ProMatrix[5] = 2.f/I->Range[1];
    I->ProMatrix[10] = -2.f/I->Range[2];
    I->ProMatrix[12] = -(I->Volume[0] + I->Volume[1])/I->Range[0];
    I->ProMatrix[13] = -(I->Volume[2] + I->Volume[3])/I->Range[1];
    I->ProMatrix[14] = -(I->Volume[4] + I->Volume[5])/I->Range[2];
  } else {
    I->ProMatrix[0] = I->Volume[4]/(front_back_ratio*I->Volume[1]);
    I->ProMatrix[5] = I->Volume[4]/(front_back_ratio*I->Volume[3]);
    I->ProMatrix[10] = -(I->Volume[5] + I->Volume[4])/I->Range[2];
    I->ProMatrix[11] = -1.f;
    I->ProMatrix[14] = (-2.f * I->Volume[5] * I->Volume[4])/I->Range[2];
    I->ProMatrix[15] = 0.f;
  }
  I->Rotation = rotMat;
  I->Ortho = ortho;
  if(ortho) {
    I->PixelRadius = (((float) I->Range[0]) / width) * pixel_scale;
  } else {
    I->PixelRadius = (((float) I->Range[0]) / width) * pixel_scale * pixel_ratio;
  }
  I->PixelRatio = pixel_ratio;
  I->Magnified = magnified;
  I->FrontBackRatio = front_back_ratio;
  I->PrimSizeCnt = 0;
  I->PrimSize = 0.0;
  I->Fov = fov;
  I->Pos = pos;

/* BEGIN PROPRIETARY CODE SEGMENT (see disclaimer in "os_proprietary.h") */
#ifdef PYMOL_EVAL
  RayDrawEvalMessage(I);
#endif

/* END PROPRIETARY CODE SEGMENT */

}


/*========================================================================*/

void RaySetTTT(CRay * I, int flag, float *ttt)
{
  I->TTTFlag = flag;
  if(I->TTTFlag) {
    I->TTT = glm::make_mat4(ttt);
  }
}

void RayGetTTT(CRay * I, float *ttt)
{
  if(!I->TTTFlag) {
    identity44f(ttt);
  } else {
    ttt = glm::value_ptr(I->TTT);
  }
}


/*========================================================================*/
void RayRelease(CRay * I)
{
  int a;

  for(a = 0; a < I->NBasis; a++) {
    BasisFinish(&I->Basis[a], a);
  }
  I->NBasis = 0;
  VLACacheFreeP(I->G, I->Primitive, 0, cCache_ray_primitive, false);
}


/*========================================================================*/
void RayFree(CRay * I)
{
  RayRelease(I);
  CharacterSetRetention(I->G, false);
  CacheFreeP(I->G, I->Basis, 0, cCache_ray_basis, false);
  DeleteP(I);
}


/*========================================================================*/
void RayPushTTT(CRay * I)
{
  if(I->TTTFlag) {
    I->TTTStack.push_back(I->TTT);
  }
}

void RayPopTTT(CRay * I)
{
  if(I->TTTStack.empty()){
    I->TTTFlag = false;
  } else {
    I->TTT = I->TTTStack.back();
    I->TTTStack.pop_back();
    I->TTTFlag = true;
  }
}

static void RayApplyMatrix33(unsigned int n, float3 * q, const float m[16], float3 * p)
{
  {
    unsigned int i;
    float m0 = m[0], m4 = m[4], m8 = m[8], m12 = m[12];
    float m1 = m[1], m5 = m[5], m9 = m[9], m13 = m[13];
    float m2 = m[2], m6 = m[6], m10 = m[10], m14 = m[14];
    for(i = 0; i < n; i++) {
      float p0 = p[i][0], p1 = p[i][1], p2 = p[i][2];
      q[i][0] = m0 * p0 + m4 * p1 + m8 * p2 + m12;
      q[i][1] = m1 * p0 + m5 * p1 + m9 * p2 + m13;
      q[i][2] = m2 * p0 + m6 * p1 + m10 * p2 + m14;
    }
  }
}

static void RayApplyMatrixInverse33(unsigned int n, float3 * q, const float m[16], float3 * p)
{
  {
    unsigned int i;
    float m0 = m[0], m4 = m[4], m8 = m[8], m12 = m[12];
    float m1 = m[1], m5 = m[5], m9 = m[9], m13 = m[13];
    float m2 = m[2], m6 = m[6], m10 = m[10], m14 = m[14];
    for(i = 0; i < n; i++) {
      float p0 = p[i][0] - m12, p1 = p[i][1] - m13, p2 = p[i][2] - m14;
      q[i][0] = m0 * p0 + m1 * p1 + m2 * p2;
      q[i][1] = m4 * p0 + m5 * p1 + m6 * p2;
      q[i][2] = m8 * p0 + m9 * p1 + m10 * p2;
    }
  }
}

static void RayTransformNormals33(unsigned int n, float3 * q, const float m[16], float3 * p)
{
  unsigned int i;
  float m0 = m[0], m4 = m[4], m8 = m[8];
  float m1 = m[1], m5 = m[5], m9 = m[9];
  float m2 = m[2], m6 = m[6], m10 = m[10];
  for(i = 0; i < n; i++) {
    float p0 = p[i][0], p1 = p[i][1], p2 = p[i][2];
    q[i][0] = m0 * p0 + m4 * p1 + m8 * p2;
    q[i][1] = m1 * p0 + m5 * p1 + m9 * p2;
    q[i][2] = m2 * p0 + m6 * p1 + m10 * p2;
  }
  for(i = 0; i < n; i++) {      /* renormalize - can we do this to the matrix instead? */
    normalize3f(q[i]);
  }
}

static void RayTransformInverseNormals33(unsigned int n, float3 * q, const float m[16],
                                  float3 * p)
{
  unsigned int i;
  float m0 = m[0], m4 = m[4], m8 = m[8];
  float m1 = m[1], m5 = m[5], m9 = m[9];
  float m2 = m[2], m6 = m[6], m10 = m[10];
  for(i = 0; i < n; i++) {
    float p0 = p[i][0], p1 = p[i][1], p2 = p[i][2];
    q[i][0] = m0 * p0 + m1 * p1 + m2 * p2;
    q[i][1] = m4 * p0 + m5 * p1 + m6 * p2;
    q[i][2] = m8 * p0 + m9 * p1 + m10 * p2;
  }
  for(i = 0; i < n; i++) {      /* renormalize - can we do this to the matrix instead? */
    normalize3f(q[i]);
  }
}
void RayGetScreenVertex(CRay * I, float *v, float *res){
  MatrixTransformC44f4f(I->ModelView, v, res);
  normalize4f(res);
}

void RayAdjustZtoScreenZ(CRay * ray, float *pos, float zarg){
  PyMOLGlobals *G = ray->G;
  float BackSafe = ray->Volume[5], FrontSafe = ray->Volume[4];
  float clipRange = (BackSafe-FrontSafe);
  float z = (zarg + 1.f) / 2.f;
  float zInPreProj = -(z * clipRange + FrontSafe);
  float pos4[4], tpos[4], npos[4];
  float InvModMatrix[16];
  copy3f(pos, pos4);
  pos4[3] = 1.f;
  MatrixTransformC44f4f(ray->ModelView, pos4, tpos);
  normalize4f(tpos);
  /* NEED TO ACCOUNT FOR ORTHO */
  if (SettingGetGlobal_b(G, cSetting_ortho)){
    npos[0] = tpos[0];
    npos[1] = tpos[1];
  } else {
    npos[0] = zInPreProj * tpos[0] / tpos[2];
    npos[1] = zInPreProj * tpos[1] / tpos[2];
  }
  npos[2] = zInPreProj;
  npos[3] = 1.f;
  MatrixInvertC44f(ray->ModelView, InvModMatrix);
  MatrixTransformC44f4f(InvModMatrix, npos, npos);
  normalize4f(npos);
  copy3f(npos, pos);
}

static float RayGetRawDepth(CRay * ray, float *pos)
{
  float pos4[4], tpos[4];
  copy3f(pos, pos4);
  pos4[3] = 1.f;
  MatrixTransformC44f4f(ray->ModelView, pos4, tpos);
  normalize4f(tpos);
  return -tpos[2];
}

void RayAdjustZtoScreenZofPoint(CRay * ray, float *pos, float *zpoint){
  float BackSafe = ray->Volume[5], FrontSafe = ray->Volume[4];
  float clipRange = (BackSafe-FrontSafe);
  float zInClipping = RayGetRawDepth(ray, zpoint);
  float z = ((2.f * (zInClipping - FrontSafe)/ clipRange) - 1.f);
  RayAdjustZtoScreenZ(ray, pos, z);
}

void RaySetPointToWorldScreenRelative(CRay * ray, float *pos, float *screenPt)
{
  float npos[4];
  float PmvMatrix[16], InvPmvMatrix[16];
  int width = ray->Width, height = ray->Height;
  multiply44f44f44f(ray->ModelView, RayGetProMatrix(ray), PmvMatrix);

  npos[0] = (floor(screenPt[0]*width)) /width ;
  npos[1] = (floor(screenPt[1]*height)) /height ;
  npos[2] = 0.f;
  npos[3] = 1.f;
  MatrixInvertC44f(PmvMatrix, InvPmvMatrix);
  MatrixTransformC44f4f(InvPmvMatrix, npos, npos);
  normalize4f(npos);
  RayAdjustZtoScreenZ(ray, npos, screenPt[2]);
  copy3f(npos, pos);
}

float RayGetScaledAllAxesAtPoint(CRay * I, float *pt, float *xn, float *yn, float *zn)
{
  float xn0[3] = { 1.0F, 0.0F, 0.0F };
  float yn0[3] = { 0.0F, 1.0F, 0.0F };
  float zn0[3] = { 0.0F, 0.0F, 1.0F };
  float v_scale;

  v_scale = RayGetScreenVertexScale(I, pt) / I->Sampling;

  RayApplyMatrixInverse33(1, (float3 *) xn0, glm::value_ptr(I->Rotation), (float3 *) xn0);
  RayApplyMatrixInverse33(1, (float3 *) yn0, glm::value_ptr(I->Rotation), (float3 *) yn0);
  RayApplyMatrixInverse33(1, (float3 *) zn0, glm::value_ptr(I->Rotation), (float3 *) zn0);

  scale3f(xn0, v_scale, xn);
  scale3f(yn0, v_scale, yn);
  scale3f(zn0, v_scale, zn);
  return v_scale;
}
float* RayGetProMatrix(CRay * I){
  return I->ProMatrix;
}
