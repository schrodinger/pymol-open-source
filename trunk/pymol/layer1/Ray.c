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
#include"os_predef.h"
#include"os_std.h"

#include"Base.h"
#include"MemoryDebug.h"
#include"Err.h"
#include"Vector.h"
#include"OOMac.h"
#include"Setting.h"
#include"Ortho.h"
#include"Util.h"
#include"Ray.h"
#include"Triangle.h" 
#include"Color.h"
#include"Matrix.h"
#include"P.h"

#ifdef _PYMOL_INLINE
#undef _PYMOL_INLINE
#include"Basis.c"
#define _PYMOL_INLINE
#endif

#ifndef RAY_SMALL
#define RAY_SMALL 0.00001 
#endif

/* BASES 
   0 contains untransformed vertices (vector size = 3)
	1 contains transformed vertices (vector size = 3)
	2 contains transformed vertices for shadowing 
*/

#define MAX_RAY_THREADS 8

typedef float float3[3];
typedef float float4[4];


static int RandomFlag=0;
static float Random[256];

void RayRelease(CRay *I);
void RayTexture(CRay *I,int mode,float *v);
void RayTransparentf(CRay *I,float v);

void RaySetup(CRay *I);
void RayColor3fv(CRay *I,float *v);
void RaySphere3fv(CRay *I,float *v,float r);
void RayCylinder3fv(CRay *I,float *v1,float *v2,float r,float *c1,float *c2);
void RaySausage3fv(CRay *I,float *v1,float *v2,float r,float *c1,float *c2);

void RayTriangle3fv(CRay *I,
						  float *v1,float *v2,float *v3,
						  float *n1,float *n2,float *n3,
						  float *c1,float *c2,float *c3);

void RayApplyMatrix33( unsigned int n, float3 *q, const float m[16],
							float3 *p );
void RayApplyMatrixInverse33( unsigned int n, float3 *q, const float m[16],
                              float3 *p );

void RayExpandPrimitives(CRay *I);
void RayTransformFirst(CRay *I);
void RayTransformBasis(CRay *I,CBasis *B);

int PrimitiveSphereHit(CRay *I,float *v,float *n,float *minDist,int except);

void RayGetSphereNormal(CRay *I,RayInfo *r);

void RayTransformNormals33( unsigned int n, float3 *q, const float m[16],float3 *p );
void RayTransformInverseNormals33( unsigned int n, float3 *q, const float m[16],float3 *p );
void RayReflectAndTexture(CRay *I,RayInfo *r);
void RayProjectTriangle(CRay *I,RayInfo *r,float *light,float *v0,float *n0,float scale);
void RayCustomCylinder3fv(CRay *I,float *v1,float *v2,float r,
                          float *c1,float *c2,int cap1,int cap2);
void RaySetContext(CRay *I,int context)
{
  I->Context=context;
}
void RayApplyContextToNormal(CRay *I,float *v);
void RayApplyContextToVertex(CRay *I,float *v);

void RayApplyContextToVertex(CRay *I,float *v)
{
  switch(I->Context) {
  case 1:
    {
      float tw;
      float th;
      
      if(I->AspRatio>1.0F) {
        tw = I->AspRatio;
        th = 1.0F;
      } else {
        th = 1.0F/I->AspRatio;
        tw = 1.0F;
      }
      v[0]+=(tw-1.0F)/2;
      v[1]+=(th-1.0F)/2;
      v[0]=v[0]*(I->Range[0]/tw)+I->Volume[0];
      v[1]=v[1]*(I->Range[1]/th)+I->Volume[2];
      v[2]=v[2]*I->Range[2]-(I->Volume[4]+I->Volume[5])/2.0F;
      RayApplyMatrixInverse33(1,(float3*)v,I->ModelView,(float3*)v);    
    }
    break;
  }
}
void RayApplyContextToNormal(CRay *I,float *v)
{
  switch(I->Context) {
  case 1:
    RayTransformInverseNormals33(1,(float3*)v,I->ModelView,(float3*)v);    
    break;
  }
}
int RayGetNPrimitives(CRay *I)    
{
  return(I->NPrimitive);
}
/*========================================================================*/
#ifdef _PYMOL_INLINE
__inline__
#endif
void RayGetSphereNormal(CRay *I,RayInfo *r)
{
  
  r->impact[0]=r->base[0]; 
  r->impact[1]=r->base[1]; 
  r->impact[2]=r->base[2]-r->dist;
  
  r->surfnormal[0]=r->impact[0]-r->sphere[0];
  r->surfnormal[1]=r->impact[1]-r->sphere[1];
  r->surfnormal[2]=r->impact[2]-r->sphere[2];
  
  normalize3f(r->surfnormal);
  
}
/*========================================================================*/
#ifdef _PYMOL_INLINE
__inline__
#endif
void RayReflectAndTexture(CRay *I,RayInfo *r)
{
  r->flat_dotgle = r->surfnormal[2];

  if(r->prim->texture)
    switch(r->prim->texture) {
    case 1:
      scatter3f(r->surfnormal,r->prim->texture_param[0]);
      break;
    case 2:
      wiggle3f(r->surfnormal,r->impact,r->prim->texture_param);
      break;
    case 3: 
      {
        float3 v;
        float3 n;
        copy3f(r->impact,v);
        RayApplyMatrixInverse33(1,&v,I->ModelView,&v);
        n[0]=(float)cos((v[0]+v[1]+v[2])*r->prim->texture_param[1]);
        n[1]=(float)cos((v[0]-v[1]+v[2])*r->prim->texture_param[1]);
        n[2]=(float)cos((v[0]+v[1]-v[2])*r->prim->texture_param[1]);
        RayTransformNormals33(1,&n,I->ModelView,&n);
        scale3f(n,r->prim->texture_param[0],n);
        add3f(n,r->surfnormal,r->surfnormal);
        normalize3f(r->surfnormal);
      }
    case 4: 
      {
        float3 v;
        float3 n;
        float *tp = r->prim->texture_param;
        copy3f(r->impact,v);
        RayApplyMatrixInverse33(1,&v,I->ModelView,&v);
        n[0]=Random[0xFF&(int)((cos((v[0])*tp[1])*256*tp[2]))];
        n[1]=Random[0xFF&(int)((cos((v[1])*tp[1])*256*tp[2]+96))];
        n[2]=Random[0xFF&(int)((cos((v[2])*tp[1])*256*tp[2]+148))];
        RayTransformNormals33(1,&n,I->ModelView,&n);
        scale3f(n,tp[0],n);
        add3f(n,r->surfnormal,r->surfnormal);
        normalize3f(r->surfnormal);
      }
      break;
    case 5: 
      {
        float3 v;
        float3 n;
        float *tp = r->prim->texture_param;
        copy3f(r->impact,v);
        RayApplyMatrixInverse33(1,&v,I->ModelView,&v);
        n[0]=Random[0xFF&(int)((v[0]*tp[1])+0)]+
          Random[0xFF&(int)((v[1]*tp[1])+20)]+
          Random[0xFF&(int)((v[2]*tp[1])+40)];
        n[1]=Random[0xFF&(int)((-v[0]*tp[1])+90)]+
          Random[0xFF&(int)((v[1]*tp[1])+100)]+
          Random[0xFF&(int)((-v[2]*tp[1])+120)];
        n[2]=Random[0xFF&(int)((v[0]*tp[1])+200)]+
          Random[0xFF&(int)((-v[1]*tp[1])+70)]+
          Random[0xFF&(int)((v[2]*tp[1])+30)];
        
        n[0]+=
          Random[0xFF&((int)((v[0]-v[1])*tp[1])+0)] +
          Random[0xFF&((int)((v[1]-v[2])*tp[1])+20)] +
          Random[0xFF&((int)((v[2]-v[0])*tp[1])+40)];
        n[1]+=
          Random[0xFF&((int)((v[0]+v[1])*tp[1])+10)]+
          Random[0xFF&((int)((v[1]+v[2])*tp[1])+90)]+
          Random[0xFF&((int)((v[2]+v[0])*tp[1])+30)];
        n[2]+=
          Random[0xFF&((int)((-v[0]+v[1])*tp[1])+220)]+
          Random[0xFF&((int)((-v[1]+v[2])*tp[1])+20)]+
          Random[0xFF&((int)((-v[2]+v[0])*tp[1])+50)];
        
        n[0]+=
          Random[0xFF&((int)((v[0]+v[1]+v[2])*tp[1])+5)]+
          Random[0xFF&((int)((v[0]+v[1]+v[2])*tp[1])+25)]+
          Random[0xFF&((int)((v[0]+v[1]+v[2])*tp[1])+46)];
        n[1]+=
          Random[0xFF&((int)((-v[0]-v[1]+v[2])*tp[1])+90)]+
          Random[0xFF&((int)((-v[0]-v[1]+v[2])*tp[1])+45)]+
          Random[0xFF&((int)((-v[0]-v[1]+v[2])*tp[1])+176)];
        n[2]+=
          Random[0xFF&((int)((v[0]+v[1]-v[2])*tp[1])+192)]+
          Random[0xFF&((int)((v[0]+v[1]-v[2])*tp[1])+223)]+
          Random[0xFF&((int)((v[0]+v[1]-v[2])*tp[1])+250)];

        RayTransformNormals33(1,&n,I->ModelView,&n);
        scale3f(n,tp[0],n);
        add3f(n,r->surfnormal,r->surfnormal);
        normalize3f(r->surfnormal);
      }
      break;
    }
  
  r->dotgle = -r->surfnormal[2]; 
  
  r->reflect[0]= - ( 2 * r->dotgle * r->surfnormal[0] );
  r->reflect[1]= - ( 2 * r->dotgle * r->surfnormal[1] );
  r->reflect[2]= -1.0F - ( 2 * r->dotgle * r->surfnormal[2] );
}
/*========================================================================*/
void RayExpandPrimitives(CRay *I)
{
  int a;
  float *v0,*v1,*n0,*n1;
  CBasis *basis;
  int nVert, nNorm;
 
  nVert=0;
  nNorm=0;
  for(a=0;a<I->NPrimitive;a++) {
	 switch(I->Primitive[a].type) {
	 case cPrimSphere:
		nVert++;
		break;
	 case cPrimCylinder:
    case cPrimSausage:
		nVert++;
		nNorm++;
		break;
	 case cPrimTriangle:
		nVert+=3;
		nNorm+=4;
		break;
	 }
  }

  basis = I->Basis;
  
  VLACheck(basis->Vertex,float,3*nVert);
  VLACheck(basis->Radius,float,nVert);
  VLACheck(basis->Radius2,float,nVert);
  VLACheck(basis->Vert2Normal,int,nVert);
  VLACheck(basis->Normal,float,3*nNorm);
  VLACheck(I->Vert2Prim,int,nVert);
  basis->MaxRadius = 0.0F;
  basis->MinVoxel = 0.0F;
  basis->NVertex=nVert;
  basis->NNormal=nNorm;

  nVert=0;
  nNorm=0;
  v0=basis->Vertex;
  n0=basis->Normal;
  for(a=0;a<I->NPrimitive;a++) {
	 switch(I->Primitive[a].type) {
	 case cPrimTriangle:
		I->Primitive[a].vert=nVert;
		I->Vert2Prim[nVert]=a;
		I->Vert2Prim[nVert+1]=a;
		I->Vert2Prim[nVert+2]=a;
		basis->Radius[nVert]=I->Primitive[a].r1;
		basis->Radius2[nVert]=I->Primitive[a].r1*I->Primitive[a].r1; /*necessary??*/
		/*		if(basis->Radius[nVert]>basis->MinVoxel)
				basis->MinVoxel=basis->Radius[nVert];*/
		if(basis->MinVoxel<0.001F)
		  basis->MinVoxel=0.001F;
		basis->Vert2Normal[nVert]=nNorm;
		basis->Vert2Normal[nVert+1]=nNorm;
		basis->Vert2Normal[nVert+2]=nNorm;
		n1=I->Primitive[a].n0;
		(*n0++)=(*n1++);
		(*n0++)=(*n1++);
		(*n0++)=(*n1++);
		n1=I->Primitive[a].n1;
		(*n0++)=(*n1++);
		(*n0++)=(*n1++);
		(*n0++)=(*n1++);
		n1=I->Primitive[a].n2;
		(*n0++)=(*n1++);
		(*n0++)=(*n1++);
		(*n0++)=(*n1++);
		n1=I->Primitive[a].n3;
		(*n0++)=(*n1++);
		(*n0++)=(*n1++);
		(*n0++)=(*n1++);
		nNorm+=4;
		v1=I->Primitive[a].v1;
		(*v0++)=(*v1++);
		(*v0++)=(*v1++);
		(*v0++)=(*v1++);
		v1=I->Primitive[a].v2;
		(*v0++)=(*v1++);
		(*v0++)=(*v1++);
		(*v0++)=(*v1++);
		v1=I->Primitive[a].v3;
		(*v0++)=(*v1++);
		(*v0++)=(*v1++);
		(*v0++)=(*v1++);
		nVert+=3;
		break;
	 case cPrimSphere:
		I->Primitive[a].vert=nVert;
		I->Vert2Prim[nVert]=a;
		v1=I->Primitive[a].v1;
		basis->Radius[nVert]=I->Primitive[a].r1;
		basis->Radius2[nVert]=I->Primitive[a].r1*I->Primitive[a].r1; /*precompute*/
		if(basis->Radius[nVert]>basis->MaxRadius)
		  basis->MaxRadius=basis->Radius[nVert];
		(*v0++)=(*v1++);
		(*v0++)=(*v1++);
		(*v0++)=(*v1++);
		nVert++;
		break;
	 case cPrimCylinder:
    case cPrimSausage:
		I->Primitive[a].vert=nVert;
		I->Vert2Prim[nVert]=a;
		basis->Radius[nVert]=I->Primitive[a].r1;
		basis->Radius2[nVert]=I->Primitive[a].r1*I->Primitive[a].r1; /*precompute*/
		if(basis->MinVoxel<0.001F)
        basis->MinVoxel=0.001F;
		subtract3f(I->Primitive[a].v2,I->Primitive[a].v1,n0);
		I->Primitive[a].l1=(float)length3f(n0);
		normalize3f(n0);
		n0+=3;
		basis->Vert2Normal[nVert]=nNorm;
		nNorm++;
		v1=I->Primitive[a].v1;
		(*v0++)=(*v1++);
		(*v0++)=(*v1++);
		(*v0++)=(*v1++);
		nVert++;
		break;
	 }
  }
  /*  printf("minvoxel  %8.3f\n",basis->MinVoxel);
		printf("NPrimit  %d nvert %d\n",I->NPrimitive,nVert);*/
}
/*========================================================================*/
static void RayComputeBox(CRay *I,float *mn,float *mx)
{

#define minmax(v,r) { \
  xp = v[0] + r;\
  xm = v[0] - r;\
  yp = v[1] + r;\
  ym = v[1] - r;\
  if(xmin>xm) xmin = xm;\
  if(xmax<xp) xmax = xp;\
  if(ymin>ym) ymin = ym;\
  if(ymax<yp) ymax = yp;\
}

  CPrimitive *prm;
  CBasis *basis1;
  basis1 = I->Basis+1;

  float xmin=0.0F,ymin=0.0F,xmax=0.0F,ymax=0.0F;
  float xp,xm,yp,ym;

  float *v,r;
  float vt[3];
  const float _0 = 0.0F;
  int a;

  if(basis1->NVertex) {
    xmin = xmax = basis1->Vertex[0];
    ymin = ymax = basis1->Vertex[1];
    
    for(a=0;a<I->NPrimitive;a++) {
      prm=I->Primitive+a;
      
      switch(prm->type) 
      {
      case cPrimTriangle:
        r = _0;
        v = basis1->Vertex + prm->vert*3;
        minmax(v,r);
        v = basis1->Vertex + prm->vert*3+3;
        minmax(v,r);
        v = basis1->Vertex + prm->vert*3+6;
        minmax(v,r);
        break;
      case cPrimSphere:
        r = prm->r1;
        v = basis1->Vertex + prm->vert*3;
        minmax(v,r);
        break;
      case cPrimCylinder:
      case cPrimSausage:
        r = prm->r1;
        v = basis1->Vertex + prm->vert*3;
        minmax(v,r);
        v = basis1->Normal+basis1->Vert2Normal[prm->vert]*3;
        scale3f(v,prm->l1,vt);
        v = basis1->Vertex + prm->vert*3;
        add3f(v,vt,vt);
        minmax(vt,r);
        break;
      }	/* end of switch */
	 }
  }
  mn[0] = xmin;
  mn[1] = ymin;
  mx[0] = xmax;
  mx[1] = ymax;
}

void RayTransformFirst(CRay *I)
{
  CBasis *basis0,*basis1;
  CPrimitive *prm;
  int a;
  float *v0;
  int backface_cull;

  backface_cull = (int)SettingGet(cSetting_backface_cull);
  
  if(SettingGet(cSetting_two_sided_lighting)||
     SettingGet(cSetting_ray_interior_color)>=0)
    backface_cull=0;

  basis0 = I->Basis;
  basis1 = I->Basis+1;
  
  VLACheck(basis1->Vertex,float,3*basis0->NVertex);
  VLACheck(basis1->Normal,float,3*basis0->NNormal);
  VLACheck(basis1->Precomp,float,3*basis0->NNormal);
  VLACheck(basis1->Vert2Normal,int,basis0->NVertex);
  VLACheck(basis1->Radius,float,basis0->NVertex);
  VLACheck(basis1->Radius2,float,basis0->NVertex);
  
  RayApplyMatrix33(basis0->NVertex,(float3*)basis1->Vertex,
					  I->ModelView,(float3*)basis0->Vertex);

  for(a=0;a<basis0->NVertex;a++)
	 {
		basis1->Radius[a]=basis0->Radius[a];
		basis1->Radius2[a]=basis0->Radius2[a];
		basis1->Vert2Normal[a]=basis0->Vert2Normal[a];
	 }
  basis1->MaxRadius=basis0->MaxRadius;
  basis1->MinVoxel=basis0->MinVoxel;
  basis1->NVertex=basis0->NVertex;

  RayTransformNormals33(basis0->NNormal,(float3*)basis1->Normal,
					  I->ModelView,(float3*)basis0->Normal);
  
  basis1->NNormal=basis0->NNormal;

  for(a=0;a<I->NPrimitive;a++) {
	 prm=I->Primitive+a;
	 if(prm->type==cPrimTriangle) {
		BasisTrianglePrecompute(basis1->Vertex+prm->vert*3,
										basis1->Vertex+prm->vert*3+3,
										basis1->Vertex+prm->vert*3+6,
										basis1->Precomp+basis1->Vert2Normal[prm->vert]*3);
      v0 = basis1->Normal + (basis1->Vert2Normal[prm->vert]*3 + 3);
      prm->cull = backface_cull&&((v0[2]<0.0F)&&(v0[5]<0.0F)&&(v0[8]<0.0F));
	 }
  }

}
/*========================================================================*/
void RayTransformBasis(CRay *I,CBasis *basis1)
{
  CBasis *basis0;
  int a;
  float *v0,*v1;
  CPrimitive *prm;

  basis0 = I->Basis+1;

  VLACheck(basis1->Vertex,float,3*basis0->NVertex);
  VLACheck(basis1->Normal,float,3*basis0->NNormal);
  VLACheck(basis1->Precomp,float,3*basis0->NNormal);
  VLACheck(basis1->Vert2Normal,int,basis0->NVertex);
  VLACheck(basis1->Radius,float,basis0->NVertex);
  VLACheck(basis1->Radius2,float,basis0->NVertex);
  v0=basis0->Vertex;
  v1=basis1->Vertex;
  for(a=0;a<basis0->NVertex;a++)
	 {
		matrix_transform33f3f(basis1->Matrix,v0,v1);
		v0+=3;
		v1+=3;
		basis1->Radius[a]=basis0->Radius[a];
		basis1->Radius2[a]=basis0->Radius2[a];
		basis1->Vert2Normal[a]=basis0->Vert2Normal[a];
	 }
  v0=basis0->Normal;
  v1=basis1->Normal;
  for(a=0;a<basis0->NNormal;a++)
	 {
		matrix_transform33f3f(basis1->Matrix,v0,v1);
		v0+=3;
		v1+=3;
	 }
  basis1->MaxRadius=basis0->MaxRadius;
  basis1->MinVoxel=basis0->MinVoxel;
  basis1->NVertex=basis0->NVertex;
  basis1->NNormal=basis0->NNormal;


  for(a=0;a<I->NPrimitive;a++) {
	 prm=I->Primitive+a;
	 if(prm->type==cPrimTriangle) {
		BasisTrianglePrecompute(basis1->Vertex+prm->vert*3,
										basis1->Vertex+prm->vert*3+3,
										basis1->Vertex+prm->vert*3+6,
										basis1->Precomp+basis1->Vert2Normal[prm->vert]*3);

	 }
  }
}

/*========================================================================*/
void RayRenderTest(CRay *I,int width,int height,float front,float back,float fov)
{

  PRINTFB(FB_Ray,FB_Details)
    " RayRenderTest: obtained %i graphics primitives.\n",I->NPrimitive 
    ENDFB;
}

/*========================================================================*/
void RayRenderPOV(CRay *I,int width,int height,char **headerVLA_ptr,
                  char **charVLA_ptr,float front,float back,float fov,
                  float angle)
{
  int antialias;
  int fogFlag=false;
  int fogRangeFlag=false;
  float fog;
  float *bkrd;
  float fog_start=0.0F;
  float gamma;
  float *d;
  CBasis *base;
  CPrimitive *prim;
  OrthoLineType buffer;
  float *vert,*norm;
  float vert2[3];
  float light[3],*lightv;
  int cc,hc;
  int a;
  int smooth_color_triangle;
  int mesh_obj = false;
  char *charVLA,*headerVLA;
  char transmit[64];

  charVLA=*charVLA_ptr;
  headerVLA=*headerVLA_ptr;
  smooth_color_triangle=(int)SettingGet(cSetting_smooth_color_triangle);
  PRINTFB(FB_Ray,FB_Details)
    " RayRenderPOV: w %d h %d f %8.3f b %8.3f\n",width,height,front,back
    ENDFB;
  if(Feedback(FB_Ray,FB_Details)) {
    dump3f(I->Volume," RayRenderPOV: vol");
    dump3f(I->Volume+3," RayRenderPOV: vol");
  }
  cc=0;
  hc=0;
  gamma = SettingGet(cSetting_gamma);
  if(gamma>R_SMALL4)
    gamma=1.0F/gamma;
  else
    gamma=1.0F;

  fog = SettingGet(cSetting_ray_trace_fog);
  if(fog!=0.0F) {
    fogFlag=true;
    fog_start = SettingGet(cSetting_ray_trace_fog_start);
    if(fog_start>R_SMALL4) {
      fogRangeFlag=true;
      if(fabs(fog_start-1.0F)<R_SMALL4) /* prevent div/0 */
        fogFlag=false;
    }
  }
  /* SETUP */
  
  antialias = (int)SettingGet(cSetting_antialias);
  bkrd=SettingGetfv(cSetting_bg_rgb);

  RayExpandPrimitives(I);
  RayTransformFirst(I);

  PRINTFB(FB_Ray,FB_Details)
    " RayRenderPovRay: processed %i graphics primitives.\n",I->NPrimitive 
    ENDFB;
  base = I->Basis+1;

  if(!SettingGet(cSetting_ortho)) {
    sprintf(buffer,"camera {direction<0.0,0.0,%8.3f>\n location <0.0 , 0.0 , 0.0>\n right %12.10f*x up y \n }\n",
            -56.6/fov,/* by trial and error */
            I->Range[0]/I->Range[1]);
  } else {
    sprintf(buffer,"camera {orthographic location <0.0 , 0.0 , %12.10f>\nlook_at  <0.0 , 0.0 , -1.0> right %12.10f*x up %12.10f*y}\n",
            front,I->Range[0],I->Range[1]);
  }
  UtilConcatVLA(&headerVLA,&hc,buffer);

  sprintf(buffer,"#default { finish{phong %8.3f ambient %8.3f diffuse %8.3f phong_size %8.6f}}\n",
          SettingGet(cSetting_spec_reflect),
          SettingGet(cSetting_ambient),
          SettingGet(cSetting_reflect)*1.2,
          SettingGet(cSetting_spec_power)/4.0F);
  UtilConcatVLA(&headerVLA,&hc,buffer);

  lightv = SettingGet_3fv(NULL,NULL,cSetting_light);
  copy3f(lightv,light);
  if(angle) {
    float temp[16];
    MatrixLoadIdentity44f(temp);
    MatrixRotate44f3f(temp,(float)-PI*angle/180,0.0F,1.0F,0.0F);
    MatrixTransform44fAs33f3f(temp,light,light);
  }
  sprintf(buffer,"light_source{<%6.4f,%6.4f,%6.4f>  rgb<1.0,1.0,1.0>}\n",
          -light[0]*10000.0F,
          -light[1]*10000.0F,
          -light[2]*10000.0F-front
          );
  UtilConcatVLA(&headerVLA,&hc,buffer);


  sprintf(buffer,"plane{z , %6.4f \n pigment{color rgb<%6.4f,%6.4f,%6.4f>}\n finish{phong 0 specular 0 diffuse 0 ambient 1.0}}\n",-back,bkrd[0],bkrd[1],bkrd[2]);
  UtilConcatVLA(&headerVLA,&hc,buffer);

  for(a=0;a<I->NPrimitive;a++) {
    prim = I->Primitive+a;
    vert = base->Vertex+3*(prim->vert);
    if(prim->type==cPrimTriangle) {
      if(smooth_color_triangle)
        if(!mesh_obj) {
          sprintf(buffer,"mesh {\n");
          UtilConcatVLA(&charVLA,&cc,buffer);        
          mesh_obj=true;
        }
    } else if(mesh_obj) {
      sprintf(buffer," pigment{color rgb <1,1,1>}}");
      UtilConcatVLA(&charVLA,&cc,buffer);     
      mesh_obj=false;
    }
    switch(prim->type) {
	 case cPrimSphere:
      sprintf(buffer,"sphere{<%12.10f,%12.10f,%12.10f>, %12.10f\n",
             vert[0],vert[1],vert[2],prim->r1);
      UtilConcatVLA(&charVLA,&cc,buffer);      
      sprintf(buffer,"pigment{color rgb<%6.4f,%6.4f,%6.4f>}}\n",
              prim->c1[0],prim->c1[1],prim->c1[2]);
      UtilConcatVLA(&charVLA,&cc,buffer);
		break;
	 case cPrimCylinder:
      d=base->Normal+3*base->Vert2Normal[prim->vert];
      scale3f(d,prim->l1,vert2);
      add3f(vert,vert2,vert2);
      sprintf(buffer,"cylinder{<%12.10f,%12.10f,%12.10f>,\n<%12.10f,%12.10f,%12.10f>,\n %12.10f\n",
              vert[0],vert[1],vert[2],
              vert2[0],vert2[1],vert2[2],
              prim->r1);
      UtilConcatVLA(&charVLA,&cc,buffer);
      sprintf(buffer,"pigment{color rgb<%6.4f1,%6.4f,%6.4f>}}\n",
              (prim->c1[0]+prim->c2[0])/2,
              (prim->c1[1]+prim->c2[1])/2,
              (prim->c1[2]+prim->c2[2])/2);
      UtilConcatVLA(&charVLA,&cc,buffer);
		break;
    case cPrimSausage:
      d=base->Normal+3*base->Vert2Normal[prim->vert];
      scale3f(d,prim->l1,vert2);
      add3f(vert,vert2,vert2);
      sprintf(buffer,"cylinder{<%12.10f,%12.10f,%12.10f>,\n<%12.10f,%12.10f,%12.10f>,\n %12.10f\nopen\n",
              vert[0],vert[1],vert[2],
              vert2[0],vert2[1],vert2[2],
              prim->r1);
      UtilConcatVLA(&charVLA,&cc,buffer);
      sprintf(buffer,"pigment{color rgb<%6.4f1,%6.4f,%6.4f>}}\n",
              (prim->c1[0]+prim->c2[0])/2,
              (prim->c1[1]+prim->c2[1])/2,
              (prim->c1[2]+prim->c2[2])/2);
      UtilConcatVLA(&charVLA,&cc,buffer);

      sprintf(buffer,"sphere{<%12.10f,%12.10f,%12.10f>, %12.10f\n",
             vert[0],vert[1],vert[2],prim->r1);
      UtilConcatVLA(&charVLA,&cc,buffer);      
      sprintf(buffer,"pigment{color rgb<%6.4f1,%6.4f,%6.4f>}}\n",
              prim->c1[0],prim->c1[1],prim->c1[2]);
      UtilConcatVLA(&charVLA,&cc,buffer);

      sprintf(buffer,"sphere{<%12.10f,%12.10f,%12.10f>, %12.10f\n",
             vert2[0],vert2[1],vert2[2],prim->r1);
      UtilConcatVLA(&charVLA,&cc,buffer);      
      sprintf(buffer,"pigment{color rgb<%6.4f1,%6.4f,%6.4f>}}\n",
              prim->c2[0],prim->c2[1],prim->c2[2]);
      UtilConcatVLA(&charVLA,&cc,buffer);

		break;
	 case cPrimTriangle:
      norm=base->Normal+3*base->Vert2Normal[prim->vert]+3;/* first normal is the average */      


      if(!TriangleDegenerate(vert,norm,vert+3,norm+3,vert+6,norm+6)) {

        if(smooth_color_triangle) {
          sprintf(buffer,"smooth_color_triangle{<%12.10f,%12.10f,%12.10f>,\n<%12.10f,%12.10f,%12.10f>,\n<%6.4f1,%6.4f,%6.4f>,\n<%12.10f,%12.10f,%12.10f>,\n<%12.10f,%12.10f,%12.10f>,\n<%6.4f1,%6.4f,%6.4f>,\n<%12.10f,%12.10f,%12.10f>,\n<%12.10f,%12.10f,%12.10f>,\n<%6.4f1,%6.4f,%6.4f> }\n",
                  vert[0],vert[1],vert[2],
                  norm[0],norm[1],norm[2],
                  prim->c1[0],prim->c1[1],prim->c1[2],
                  vert[3],vert[4],vert[5],
                  norm[3],norm[4],norm[5],
                  prim->c2[0],prim->c2[1],prim->c2[2],
                  vert[6],vert[7],vert[8],
                  norm[6],norm[7],norm[8],
                  prim->c3[0],prim->c3[1],prim->c3[2]
                  );
          UtilConcatVLA(&charVLA,&cc,buffer);
        } else {
          sprintf(buffer,"smooth_triangle{<%12.10f,%12.10f,%12.10f>,\n<%12.10f,%12.10f,%12.10f>,\n<%12.10f,%12.10f,%12.10f>,\n<%12.10f,%12.10f,%12.10f>,\n<%12.10f,%12.10f,%12.10f>,\n<%12.10f,%12.10f,%12.10f>\n",
                  vert[0],vert[1],vert[2],
                  norm[0],norm[1],norm[2],
                  vert[3],vert[4],vert[5],
                  norm[3],norm[4],norm[5],
                  vert[6],vert[7],vert[8],
                  norm[6],norm[7],norm[8]
                  );
          UtilConcatVLA(&charVLA,&cc,buffer);
          if(prim->trans>R_SMALL4) 
            sprintf(transmit,"transmit %4.6f",prim->trans);
          else
            transmit[0]=0;
          if(equal3f(prim->c1,prim->c2)||equal3f(prim->c1,prim->c3)) {
            sprintf(buffer,"pigment{color rgb<%6.4f1,%6.4f,%6.4f> %s}}\n",
                    prim->c1[0],prim->c1[1],prim->c1[2],transmit);
          } else if(equal3f(prim->c2,prim->c3)) {
            sprintf(buffer,"pigment{color rgb<%6.4f1,%6.4f,%6.4f> %s}}\n",
                    prim->c2[0],prim->c2[1],prim->c2[2],transmit);
          } else {
            sprintf(buffer,"pigment{color rgb<%6.4f1,%6.4f,%6.4f> %s}}\n",
                    (prim->c1[0]+prim->c2[0]+prim->c3[0])/3,
                  (prim->c1[1]+prim->c2[1]+prim->c3[1])/3,
                    (prim->c1[2]+prim->c2[2]+prim->c3[2])/3,transmit);
          }
        UtilConcatVLA(&charVLA,&cc,buffer);
        }
      }
		break;
    }
  }
  
  if(mesh_obj) {
    sprintf(buffer," pigment{color rgb <1,1,1>}}");
    UtilConcatVLA(&charVLA,&cc,buffer);     
    mesh_obj=false;
  }
  *charVLA_ptr=charVLA;
  *headerVLA_ptr=headerVLA;
}

/*========================================================================*/
void RayProjectTriangle(CRay *I,RayInfo *r,float *light,float *v0,float *n0,float scale)
{
  float w2;
  float d1[3],d2[3],d3[3];
  float p1[3],p2[3],p3[3];
  int c=0;

  if(dot_product3f(light,n0-3)>=0.0F) c++;  
  if(dot_product3f(light,n0)>=0.0F) c++;
  if(dot_product3f(light,n0+3)>=0.0F) c++;
  if(dot_product3f(light,n0+6)>=0.0F) c++;
  
  if(c) {

    w2 = 1.0F-(r->tri1+r->tri2);
    
    subtract3f(v0,r->impact,d1);
    project3f(d1,n0,p1);
    scale3f(p1,w2,d1);
    
    subtract3f(v0+3,r->impact,d2);
    project3f(d2,n0+3,p2);
    scale3f(p2,r->tri1,d2);
    
    subtract3f(v0+6,r->impact,d3);
    project3f(d3,n0+6,p3);
    scale3f(p3,r->tri2,d3);
    
    add3f(d1,d2,d2);
    add3f(d2,d3,d3);
    scale3f(d3,scale,d3);
    if(dot_product3f(r->surfnormal,d3)>=0.0F)
      add3f(d3,r->impact,r->impact);
  }
}

static void RayHashSpawn(CRayHashThreadInfo *Thread,int n_thread)
{
  int blocked;
  PyObject *info_list;
  int a;
  blocked = PAutoBlock();

  PRINTFB(FB_Ray,FB_Details)
    " Ray: filling voxels with %d threads...\n",n_thread
  ENDFB;
  info_list = PyList_New(n_thread);
  for(a=0;a<n_thread;a++) {
    PyList_SetItem(info_list,a,PyCObject_FromVoidPtr(Thread+a,NULL));
  }
  PyObject_CallMethod(P_cmd,"_ray_hash_spawn","O",info_list);
  Py_DECREF(info_list);
  PAutoUnblock(blocked);
}

int RayHashThread(CRayHashThreadInfo *T)
{
  BasisMakeMap(T->basis,T->vert2prim,T->prim,T->clipBox);
  return 1;
}

static void RayTraceSpawn(CRayThreadInfo *Thread,int n_thread)
{
  int blocked;
  PyObject *info_list;
  int a;
  blocked = PAutoBlock();

  PRINTFB(FB_Ray,FB_Details)
    " Ray: rendering with %d threads...\n",n_thread
  ENDFB;
  info_list = PyList_New(n_thread);
  for(a=0;a<n_thread;a++) {
    PyList_SetItem(info_list,a,PyCObject_FromVoidPtr(Thread+a,NULL));
  }
  PyObject_CallMethod(P_cmd,"_ray_spawn","O",info_list);
  Py_DECREF(info_list);
  PAutoUnblock(blocked);
  
}

#if 0 
static void RayTraceStitch(CRayThreadInfo *Thread,int n_thread,unsigned int *image)
{
	int cnt,start_at,done_at;
	CRayThreadInfo *T;
	unsigned int *p;
	int phase;
	int x,y;
	unsigned int n_pixel;
	int	wid, hgt, row, col, found, cPhase, threadCnt;
	
	n_pixel = Thread->width * Thread->height;

	for(phase = 0; phase < n_thread; phase++)
	{
		T			= Thread + phase;
		p			= T->image;
		cPhase		= T->phase;
		threadCnt	= T->n_thread;
		wid			= T->width;
		hgt			= T->height;
		
		for(y = 0; (y < hgt); y++)
		{				
			if((y % threadCnt) == cPhase) 
			{
				unsigned int	*pDst	= image + (wid*y);
				for(x = 0; (x < wid); x++)
					*(pDst++) = *(p++);
			}
		}
	}
}
#endif

int RayTraceThread(CRayThreadInfo *T)
{
	CRay *I=NULL;
	int x,y,yy;
	float excess=0.0F;
	float dotgle;
	float bright,direct_cmp,reflect_cmp,fc[4];
	float ambient,direct,lreflect,ft,ffact,ffact1m;
	unsigned int cc0,cc1,cc2,cc3;
	int i;
	RayInfo r1,r2;
	int fogFlag=false;
	int fogRangeFlag=false;
	int opaque_back=0;
	int n_hit=0;
	int two_sided_lighting;
	float fog;
	float *inter=NULL;
	float fog_start=0.0F;
	float gamma,inp,sig=1.0F;
	float persist,persist_inv;
	float new_front;
	int pass;
	unsigned int last_pixel=0,*pixel;
	int exclude;
	float lit;
	int backface_cull;
	float project_triangle;
	float excl_trans;
	int shadows;
	int trans_shadows;
	float first_excess;
	int pixel_flag;
	float ray_trans_spec;
	float shadow_fudge;
	int interior_color;
	int interior_flag;
	int interior_shadows;
	int interior_texture;
	int dummy;
	int texture_save;
	MapCache cache,shadow_cache;
	float		settingPower, settingReflectPower,settingSpecPower,settingSpecReflect, _0, _1, _p5, _255, _persistLimit, _inv3;
	float		invHgt, invFrontMinusBack, inv1minusFogStart,invWdth,invHgtRange;
	register float       invWdthRange,vol0;
	float       vol2;
	CBasis      *bp1,*bp2;
	int render_height;
	int offset=0;

	_0		= 0.0F;
	_1		= 1.0F;
	_p5		= 0.5F;
	_255	= 255.0F;
	_inv3	= _1/3.0F;
	_persistLimit	= 0.0001F;
  	
	/* SETUP */
	
	/*  if(T->n_thread>1)
	printf(" Ray: Thread %d: Spawned.\n",T->phase+1);
	*/
	I = T->ray;
	
	interior_shadows	= (int)SettingGet(cSetting_ray_interior_shadows);
	interior_texture	= (int)SettingGet(cSetting_ray_interior_texture);
	interior_color		= (int)SettingGet(cSetting_ray_interior_color);
	project_triangle	= SettingGet(cSetting_ray_improve_shadows);
	shadows				= (int)SettingGet(cSetting_ray_shadows);
	trans_shadows		= (int)SettingGet(cSetting_ray_transparency_shadows);
	backface_cull		= (int)SettingGet(cSetting_backface_cull);
	opaque_back			= (int)SettingGet(cSetting_ray_opaque_background);
	two_sided_lighting	= (int)SettingGet(cSetting_two_sided_lighting);
	ray_trans_spec		= SettingGet(cSetting_ray_transparency_specular);
	ambient				= SettingGet(cSetting_ambient);
	lreflect			= SettingGet(cSetting_reflect);
	direct				= SettingGet(cSetting_direct);

	/* COOP */
	settingPower		= SettingGet(cSetting_power);
	settingReflectPower	= SettingGet(cSetting_reflect_power);
	settingSpecPower	= SettingGet(cSetting_spec_power);
	settingSpecReflect	= SettingGet(cSetting_spec_reflect);	
	
	MapCacheInit(&cache,I->Basis[1].Map);
	if(shadows&&(I->NBasis>1))
		MapCacheInit(&shadow_cache,I->Basis[2].Map);
    
	if((interior_color>=0)||(two_sided_lighting))
		backface_cull	= 0;

	shadow_fudge = SettingGet(cSetting_ray_shadow_fudge);
	
	gamma = SettingGet(cSetting_gamma);
	if(gamma > R_SMALL4)
		gamma	= _1/gamma;
	else
		gamma	= _1;
	
	inv1minusFogStart	= _1;
	
	fog = SettingGet(cSetting_ray_trace_fog);

	if(fog != _0) 
	{
		fogFlag	= true;
		fog_start = SettingGet(cSetting_ray_trace_fog_start);
		if(fog_start>R_SMALL4) 
		{
			fogRangeFlag=true;
			if(fabs(fog_start - _1) < R_SMALL4) /* prevent div/0 */
				fogFlag = false;
		}
		inv1minusFogStart	= _1 / (_1 - fog_start);
	}

	if(interior_color>=0)
		inter = ColorGet(interior_color);

	/* ray-trace */
	
	invHgt				= _1 / (float) T->height;
	invFrontMinusBack	= _1 / (T->front - T->back);
	invWdth             = _1 / (float) T->width;
	invWdthRange        = invWdth * I->Range[0];
	invHgtRange         = invHgt * I->Range[1];
	vol0 = I->Volume[0];
	vol2 = I->Volume[2];
	bp1  = I->Basis + 1;
	bp2  = I->Basis + 2;
	
   render_height = T->y_stop - T->y_start;

   if(render_height) {
     offset = (T->phase * render_height/T->n_thread);
     offset = offset - (offset % T->n_thread) + T->phase;
   }
     
	r1.base[2]	= _0;
		
	for(yy = T->y_start; (yy < T->y_stop); yy++)
	{
		
      y = T->y_start + ((yy-T->y_start) + offset) % ( render_height); /* make sure threads write to different pages */

		if((!T->phase)&&!(yy & 0xF))
			OrthoBusyFast(y,T->height); /* don't slow down rendering too much */
				
		pixel = T->image + (T->width * y) + T->x_start;
	
		if((y % T->n_thread) == T->phase)	/* this is my scan line */
		{	
			r1.base[1]	= (y * invHgtRange) + vol2;
			
			for(x = T->x_start; (x < T->x_stop); x++)
			{
				exclude		= -1;
				
				r1.base[0]	= (x * invWdthRange)  + vol0;
				
				persist			= _1;
				first_excess	= _0;
				excl_trans		= _0;
				pass			= 0;
				new_front		= T->front;
				
				while((persist > _persistLimit) && (pass < 25))
				{
					pixel_flag		= false;
					
					interior_flag	= (interior_color >= 0);
					
					i	= BasisHit( bp1, &r1, exclude, I->Vert2Prim, I->Primitive,
									false, new_front, T->back, excl_trans,
									trans_shadows, &interior_flag, &cache);
							   
				    if((i >= 0) || interior_flag) 
					{
						pixel_flag		= true;
						n_hit++;
						
						if(interior_flag)
						{
							copy3f(r1.base,r1.impact);
							r1.surfnormal[0]	= _0;
							r1.surfnormal[1]	= _0;
							r1.surfnormal[2]	= _1;
							r1.impact[2]	-= T->front;
							
							if(interior_texture >= 0) 
							{
								texture_save		= r1.prim->texture; /* This is a no-no for multithreading! */
								r1.prim->texture	= interior_texture;
								
								RayReflectAndTexture(I,&r1);
								
								r1.prim->texture	= texture_save;
							}
							else
								RayReflectAndTexture(I,&r1);
	
							dotgle = -r1.dotgle;
							copy3f(inter,fc);
						}
						else
						{
							new_front	= r1.dist;
							if(r1.prim->type==cPrimTriangle) 
							{
								BasisGetTriangleNormal(bp1,&r1,i,fc);
								
								RayProjectTriangle(I, &r1, bp2->LightNormal,
													bp1->Vertex+i*3,
													bp1->Normal+bp1->Vert2Normal[i]*3+3,
													project_triangle);
													
								RayReflectAndTexture(I,&r1);
								BasisGetTriangleFlatDotgle(bp1,&r1,i);
							}
							else 
							{
								RayGetSphereNormal(I,&r1);
								RayReflectAndTexture(I,&r1);
								
								if((r1.prim->type==cPrimCylinder) || (r1.prim->type==cPrimSausage)) 
								{
									ft = r1.tri1;
									fc[0]=(r1.prim->c1[0]*(_1-ft))+(r1.prim->c2[0]*ft);
									fc[1]=(r1.prim->c1[1]*(_1-ft))+(r1.prim->c2[1]*ft);
									fc[2]=(r1.prim->c1[2]*(_1-ft))+(r1.prim->c2[2]*ft);
								}
								else 
								{
									fc[0]=r1.prim->c1[0];
									fc[1]=r1.prim->c1[1];
									fc[2]=r1.prim->c1[2];
								}
							}
							dotgle=-r1.dotgle;
							
							if(r1.flat_dotgle < _0)
							{
								if((!two_sided_lighting) && (interior_color>=0)) 
								{
									interior_flag		= true;
									r1.surfnormal[0]	= _0;
									r1.surfnormal[1]	= _0;
									r1.surfnormal[2]	= _1;
									copy3f(r1.base,r1.impact);
									r1.impact[2]		-= T->front;
									r1.dist				= T->front;
							
									if(interior_texture >= 0)
									{
										texture_save		= r1.prim->texture;
										r1.prim->texture	= interior_texture;
										RayReflectAndTexture(I,&r1);
										r1.prim->texture	= texture_save;
									}
									else
										RayReflectAndTexture(I,&r1);
							
									dotgle	= -r1.dotgle;
									copy3f(inter,fc);
								}
							}
							
							if((dotgle < _0) && (!interior_flag))
							{
								if(two_sided_lighting) 
								{
									dotgle	= -dotgle;
									invert3f(r1.surfnormal);
								}
								else 
									dotgle	= _0;
							}
						}
						
						direct_cmp = (float) ( (dotgle + (pow(dotgle, settingPower))) * _p5 );
					  
						lit = _1;
						
						if(shadows&&((!interior_flag)||(interior_shadows))) 
						{
							matrix_transform33f3f(bp2->Matrix,r1.impact,r2.base);
							r2.base[2]-=shadow_fudge;
							if(BasisHit(bp2,&r2,i,I->Vert2Prim,I->Primitive, true,_0,_0,_0,trans_shadows,&dummy, &shadow_cache) >= 0) 
							{
								lit	= (float) pow(r2.prim->trans, _p5);
							} 
						}
					  
						if(lit>_0)
						{
							dotgle	= -dot_product3f(r1.surfnormal,bp2->LightNormal);
							if(dotgle < _0) dotgle = _0;
							
							reflect_cmp	=(float)(lit * (dotgle + (pow(dotgle, settingReflectPower))) * _p5 );
							dotgle	= -dot_product3f(r1.surfnormal,T->spec_vector);
							if(dotgle < _0) dotgle=_0;
							excess	= (float)( pow(dotgle, settingSpecPower) * settingSpecReflect * lit);
						}
						else 
						{
							excess		= _0;
							reflect_cmp	= _0;
						}
					  
						bright	= ambient + (_1-ambient) * (direct*direct_cmp + (_1-direct)*direct_cmp*lreflect*reflect_cmp);
					  
						if(bright > _1)			bright = _1;
						else if(bright < _0)	bright = _0;
						
						fc[0] = (bright*fc[0]+excess);
						fc[1] = (bright*fc[1]+excess);
						fc[2] = (bright*fc[2]+excess);
				  
						if (fogFlag) 
						{
							ffact = fog*(T->front - r1.dist) * invFrontMinusBack;
							if(fogRangeFlag)
								ffact = (ffact - fog_start) * inv1minusFogStart;
	
							if(ffact<_0)	ffact = _0;
							if(ffact>_1)	ffact = _0;
							
							ffact1m	= _1-ffact;
							
							if(opaque_back) 
							{
								fc[0]	= ffact*T->bkrd[0]+fc[0]*ffact1m;
								fc[1]	= ffact*T->bkrd[1]+fc[1]*ffact1m;
								fc[2]	= ffact*T->bkrd[2]+fc[2]*ffact1m;
							}
							else 
							{
								fc[3]=ffact1m;
							}
							
							if(!pass)
								first_excess = excess*ffact1m*ray_trans_spec;
							else
							{
								fc[0]+=first_excess;
								fc[1]+=first_excess;
								fc[2]+=first_excess;
							}
						}
						else 
						{
							if(!pass)
								first_excess = excess*ray_trans_spec;
							else 
							{
								fc[0]	+= first_excess;
								fc[1]	+= first_excess;
								fc[2]	+= first_excess;
							}
							fc[3]	= _1;
						}
					}
					else if(pass) 
					{
						/* hit nothing, and we're on on second or greater pass */
						fc[0] = first_excess+T->bkrd[0];
						fc[1] = first_excess+T->bkrd[1];
						fc[2] = first_excess+T->bkrd[2];
						fc[3] = _1;

						pixel_flag	= true;
					}
				
					if(pixel_flag)
					{
						inp	= (fc[0]+fc[1]+fc[2]) * _inv3;
						
						if(inp < R_SMALL4) 
							sig = _1;
						else
							sig = (float)(pow(inp,gamma) / inp);
						
						cc0 = (uint)(sig * fc[0] * _255);
						cc1 = (uint)(sig * fc[1] * _255);
						cc2 = (uint)(sig * fc[2] * _255);
						
						if(cc0 > 255) cc0 = 255;
						if(cc1 > 255) cc1 = 255;
						if(cc2 > 255) cc2 = 255;
						
						if(opaque_back) 
						{ 
							if(I->BigEndian) 
								*pixel = T->fore_mask|(cc0<<24)|(cc1<<16)|(cc2<<8);
							else
								*pixel = T->fore_mask|(cc2<<16)|(cc1<<8)|cc0;
						}
						else	/* use alpha channel for fog with transparent backgrounds */
						{
							cc3	= (uint)(fc[3] * _255);
							if(cc3 > 255) cc3 = 255;
							
							if(I->BigEndian)
								*pixel = (cc0<<24)|(cc1<<16)|(cc2<<8)|cc3;
							else
								*pixel = (cc3<<24)|(cc2<<16)|(cc1<<8)|cc0;
						}
					}
	
					if(pass)	/* average all four channels */
					{	
						persist_inv = _1-persist;
						fc[0]	= (0xFF&((*pixel)>>24)) * persist + (0xFF&(last_pixel>>24))*persist_inv;
						fc[1]	= (0xFF&((*pixel)>>16)) * persist + (0xFF&(last_pixel>>16))*persist_inv;
						fc[2]	= (0xFF&((*pixel)>>8))  * persist + (0xFF&(last_pixel>>8))*persist_inv;
						fc[3]	= (0xFF&((*pixel)))     * persist + (0xFF&(last_pixel))*persist_inv;
						
						cc0		= (uint)(fc[0]);
						cc1		= (uint)(fc[1]);
						cc2		= (uint)(fc[2]);
						cc3		= (uint)(fc[3]);
						
						if(cc0 > 255) cc0	= 255;
						if(cc1 > 255) cc1	= 255;
						if(cc2 > 255) cc2	= 255;
						if(cc3 > 255) cc3	= 255;
						
						*pixel = (cc0<<24)|(cc1<<16)|(cc2<<8)|cc3;
						
					}
					
					if(i >= 0)
					{
						if(r1.prim->type == cPrimSausage)	/* carry ray through the stick */
							excl_trans	= new_front+(2*r1.surfnormal[2]*r1.prim->r1);
						
						if(!backface_cull) 
							persist	= persist * r1.prim->trans;
						else 
						{
							if((persist < 0.9999) && (r1.prim->trans))	/* don't combine transparent surfaces */
								*pixel	= last_pixel;
							else
								persist	= persist * r1.prim->trans;
						}

					}
					
					if( i < 0 )	/* nothing hit */
					{
						break;
					}
					else 
					{
						last_pixel	= *pixel;
						exclude		= i;
						pass++;
					}

				} /* end of while */
				
				pixel++;

			}	/* end of for */
			
		}	/* end of if */
		
	}	/* end of for */
	
	/*  if(T->n_thread>1) 
	  printf(" Ray: Thread %d: Complete.\n",T->phase+1);*/
	MapCacheFree(&cache);
	
	if(shadows&&(I->NBasis>1))
		MapCacheFree(&shadow_cache);
	
	return (n_hit);
}





typedef struct
{
	unsigned int				*buff;
	unsigned int				buffSize;
} RenderBuff;

#define kRenderCacheSize		(MAX_RAY_THREADS + 1)
#define	kMaxBuffIndex			(kRenderCacheSize - 1)
static RenderBuff				RenderCache[kRenderCacheSize];

static unsigned int *GetRenderBuffer(unsigned int inBuffSize, int inThreadIndex)
{
	static int RenderCacheInited = 0;

	/* WARREN -- TO DO -- STRIP OUT THREADING STUFF, and ADD A CLEANUP ROUTINE */

   /*	if (inThreadIndex > kMaxBuffIndex)
		printf(" DRSWAT: OoR thread index: %d.\n", inThreadIndex);*/

	if (RenderCacheInited == 0)
	{
		int i;
		for (i = 0; i < kRenderCacheSize; ++i)
		{
			RenderCache[i].buff = NULL;
			RenderCache[i].buffSize = 0;
		}

		RenderCacheInited = 1;
	}

	if (RenderCache[inThreadIndex].buff == NULL)
	{
     /*		printf("DRSWAT: Allocated render buffer first time. Thread: %d. Size: %d.\n", inThreadIndex, inBuffSize); */

		RenderCache[inThreadIndex].buff = (void *)Alloc(char, inBuffSize);
		RenderCache[inThreadIndex].buffSize = inBuffSize;

		return RenderCache[inThreadIndex].buff;
	}
	else
	{
		if (RenderCache[inThreadIndex].buffSize == inBuffSize)
			return RenderCache[inThreadIndex].buff;
		else
		{
			printf("DRSWAT: Render buffer miss. Thread: %d. Size: %d.\n", inThreadIndex, inBuffSize);

			FreeP(RenderCache[inThreadIndex].buff);
			RenderCache[inThreadIndex].buff = (void *)Alloc(char, inBuffSize);
			RenderCache[inThreadIndex].buffSize = inBuffSize;

			return RenderCache[inThreadIndex].buff;
		}
	}
   /*
     printf("DRSWAT: Oh no! No render buffer found.\n");*/
	return NULL;
}

/* this is both an antialias and a slight blur */

/* for whatever reason, greatly GCC perfers a linear sequence of
   accumulates over a single large expression -- the difference is
   huge: over 10% !!! */

#define combine(var,src,mask) { \
  var = ((src[0 ] & mask)*5);\
  var += ((src[1 ] & mask)*5);\
  var += ((src[2 ] & mask)*5);\
  var += ((src[3 ] & mask)*5);\
  var += (src[4 ] & mask);\
  var += (src[5 ] & mask);\
  var += (src[6 ] & mask);\
  var += (src[7 ] & mask);\
  var += (src[8 ] & mask);\
  var += (src[9 ] & mask);\
  var += (src[10] & mask);\
  var += (src[11] & mask);\
  var += (src[12] & mask);\
  var += (src[13] & mask);\
  var += (src[14] & mask);\
  var += (src[15] & mask);\
  var >>= 5; \
  var &= mask; \
 }

/* this is just a simple antialias */

#define edge_combine(var,src,mask) \
{\
  var =  (src[0 ] & mask);\
  var += (src[1 ] & mask);\
  var += (src[2 ] & mask);\
  var += (src[3 ] & mask);\
  var >>= 4; \
  var &= mask;\
 }



static void fill(unsigned int *buffer, unsigned int value,unsigned int cnt)
{
  while(cnt&0xFFFFFFF80) {
    *(buffer++) = value;
    *(buffer++) = value;
    *(buffer++) = value;
    *(buffer++) = value;
    *(buffer++) = value;
    *(buffer++) = value;
    *(buffer++) = value;
    *(buffer++) = value;
    cnt-=0x20;
    *(buffer++) = value;
    *(buffer++) = value;
    *(buffer++) = value;
    *(buffer++) = value;
    *(buffer++) = value;
    *(buffer++) = value;
    *(buffer++) = value;
    *(buffer++) = value;
    *(buffer++) = value;
    *(buffer++) = value;
    *(buffer++) = value;
    *(buffer++) = value;
    *(buffer++) = value;
    *(buffer++) = value;
    *(buffer++) = value;
    *(buffer++) = value;
    *(buffer++) = value;
    *(buffer++) = value;
    *(buffer++) = value;
    *(buffer++) = value;
    *(buffer++) = value;
    *(buffer++) = value;
    *(buffer++) = value;
    *(buffer++) = value;
  }
  while(cnt--) {
    *(buffer++) = value;
  }
}
/*========================================================================*/
void RayRender(CRay *I,int width,int height,unsigned int *image,
               float front,float back,
               double timing,float angle)
{
  int x,y;
  int a;
  unsigned int *p;
  float *v,light[3];
  unsigned int aa;
 unsigned int *image_copy = NULL;
  unsigned int back_mask,fore_mask=0;
  unsigned int background,buffer_size;
  int antialias;
  int opaque_back=0;
  int n_hit=0;
  float *bkrd;
  double now;
  int shadows;
  float spec_vector[3];
  int n_thread;
  float min_box[2],max_box[2];

  n_thread  = (int)SettingGet(cSetting_max_threads);
  if(n_thread<1)
    n_thread=1;
  if(n_thread>MAX_RAY_THREADS)
    n_thread = MAX_RAY_THREADS;
  opaque_back = (int)SettingGet(cSetting_ray_opaque_background);
  BasisSetFudge(SettingGet(cSetting_ray_triangle_fudge));
  shadows = (int)SettingGet(cSetting_ray_shadows);
  antialias = (int)SettingGet(cSetting_antialias);
  if(antialias>0) {
	 width=width*2;
	 height=height*2;
	 image_copy = image;
	 buffer_size = 4*width*height;
	image = GetRenderBuffer(buffer_size * sizeof(char), MAX_RAY_THREADS);
	/* image = (void *)Alloc(char, buffer_size); */
	 ErrChkPtr(image);
  }

  bkrd=SettingGetfv(cSetting_bg_rgb);
  if(opaque_back) {
    if(I->BigEndian)
      back_mask = 0x000000FF;
    else
      back_mask = 0xFF000000;
    fore_mask = back_mask;
  } else {
    if(I->BigEndian) {
      back_mask = 0x00000000;
    } else {
      back_mask = 0x00000000;
    }
  }
  if(I->BigEndian) {
    background = back_mask|
      ((0xFF& ((unsigned int)(bkrd[0]*255))) <<24)|
      ((0xFF& ((unsigned int)(bkrd[1]*255))) <<16)|
      ((0xFF& ((unsigned int)(bkrd[2]*255))) <<8 );
  } else {
    background = back_mask|
      ((0xFF& ((unsigned int)(bkrd[2]*255))) <<16)|
      ((0xFF& ((unsigned int)(bkrd[1]*255))) <<8)|
      ((0xFF& ((unsigned int)(bkrd[0]*255))) );
  }

  PRINTFB(FB_Ray,FB_Blather) 
    " RayNew: Background = %x %d %d %d\n",background,(int)(bkrd[0]*255),
    (int)(bkrd[1]*255),(int)(bkrd[2]*255)
    ENDFB;

  if(!I->NPrimitive) { /* nothing to render! */
    p=(unsigned int*)image; 
    for(x=0;x<width;x++)
      for(y=0;y<height;y++)
        *(p++)=background;
  } else {
    
    RayExpandPrimitives(I);
    RayTransformFirst(I);
    RayComputeBox(I,min_box,max_box);
    
    now = UtilGetSeconds()-timing;

	 PRINTFB(FB_Ray,FB_Blather)
      " Ray: processed %i graphics primitives in %4.2f sec.\n",I->NPrimitive,now
      ENDFB;

    I->NBasis=3; /* light source */
    BasisInit(I->Basis+2);
    
    { /* setup light & rotate if necessary  */
      v=SettingGetfv(cSetting_light);
      copy3f(v,light);
      
      if(angle) {
        float temp[16];
        MatrixLoadIdentity44f(temp);
        MatrixRotate44f3f(temp,(float)-PI*angle/180,0.0F,1.0F,0.0F);
        MatrixTransform44fAs33f3f(temp,light,light);
      }
      
      I->Basis[2].LightNormal[0]=light[0];
      I->Basis[2].LightNormal[1]=light[1];
      I->Basis[2].LightNormal[2]=light[2];
      normalize3f(I->Basis[2].LightNormal);
      
      copy3f(I->Basis[2].LightNormal,spec_vector);
      spec_vector[2]--; /* HUH? */
      normalize3f(spec_vector);
      
    }

    if(shadows) { /* don't waste time on shadows unless needed */
      BasisSetupMatrix(I->Basis+2);
      RayTransformBasis(I,I->Basis+2);
    }

    if(shadows&&(n_thread>1)) { /* parallel execution */

      CRayHashThreadInfo thread_info[2];
      
      thread_info[0].basis = I->Basis+1;
      thread_info[0].vert2prim = I->Vert2Prim;
      thread_info[0].prim = I->Primitive;
      thread_info[0].clipBox = I->Volume;

      thread_info[1].basis = I->Basis+2;
      thread_info[1].vert2prim = I->Vert2Prim;
      thread_info[1].prim = I->Primitive;
      thread_info[1].clipBox = NULL;

      RayHashSpawn(thread_info,2);
      
    } else { /* serial execution */
      BasisMakeMap(I->Basis+1,I->Vert2Prim,I->Primitive,I->Volume);
      if(shadows) {
        BasisMakeMap(I->Basis+2,I->Vert2Prim,I->Primitive,NULL);
      }
    }

    now = UtilGetSeconds()-timing;

#ifdef _MemoryDebug_ON
    if(shadows) {
      PRINTFB(FB_Ray,FB_Details)
        " Ray: voxels: [%4.2f:%dx%dx%d], [%4.2f:%dx%dx%d], %d MB, %4.2f sec.\n",
        I->Basis[1].Map->Div,   I->Basis[1].Map->Dim[0],
        I->Basis[1].Map->Dim[1],I->Basis[1].Map->Dim[2],
        I->Basis[2].Map->Div,   I->Basis[2].Map->Dim[0],
        I->Basis[2].Map->Dim[2],I->Basis[2].Map->Dim[2],
        (int)(MemoryDebugUsage()/(1024.0*1024)),
        now
        ENDFB;
    } else {
      PRINTFB(FB_Ray,FB_Details)
        " Ray: voxels: [%4.2f:%dx%dx%d], %d MB, %4.2f sec.\n",
        I->Basis[1].Map->Div,   I->Basis[1].Map->Dim[0],
        I->Basis[1].Map->Dim[1],I->Basis[1].Map->Dim[2],
        (int)(MemoryDebugUsage()/(1024.0*1024)),
        now
        ENDFB;
    }
#else
    if(shadows) {
      PRINTFB(FB_Ray,FB_Details)
        " Ray: voxels: [%4.2f:%dx%dx%d], [%4.2f:%dx%dx%d], %4.2f sec.\n",
        I->Basis[1].Map->Div,   I->Basis[1].Map->Dim[0],
        I->Basis[1].Map->Dim[1],I->Basis[1].Map->Dim[2],
        I->Basis[2].Map->Div,   I->Basis[2].Map->Dim[0],
        I->Basis[2].Map->Dim[2],I->Basis[2].Map->Dim[2],
        now
        ENDFB;
    } else {
      PRINTFB(FB_Ray,FB_Details)
        " Ray: voxels: [%4.2f:%dx%dx%d], %4.2f sec.\n",
        I->Basis[1].Map->Div,   I->Basis[1].Map->Dim[0],
        I->Basis[1].Map->Dim[1],I->Basis[1].Map->Dim[2],
        now
        ENDFB;
    }

#endif
    /* IMAGING */
 
 
    fill(image,background,width * (unsigned int)height);
        
    {
		/* now spawn threads as needed */
		CRayThreadInfo rt[MAX_RAY_THREADS];
      int x_start,y_start;
      int x_stop,y_stop;

      x_start = (int)((width * (min_box[0] - I->Volume[0]))/I->Range[0]) - 1;
      x_stop  = (int)((width * (max_box[0] - I->Volume[0]))/I->Range[0]) + 2;
      
      y_stop = (int)((height * (max_box[1] - I->Volume[2]))/I->Range[1]) + 2;
      y_start  = (int)((height * (min_box[1] - I->Volume[2]))/I->Range[1]) - 1;
                      
      if(x_start<0) x_start = 0;
      if(y_start<0) y_start = 0;
      if(x_stop>width) x_stop = width;
      if(y_stop>height) y_stop = height;

		for(a=0;a<n_thread;a++) 
		{
			rt[a].ray = I;
			rt[a].width = width;
			rt[a].height = height;
         rt[a].x_start = x_start;
         rt[a].x_stop = x_stop;
         rt[a].y_start = y_start;
         rt[a].y_stop = y_stop;
			rt[a].image = image;
			  
			rt[a].front = front;
			rt[a].back = back;
			rt[a].fore_mask = fore_mask;
			rt[a].bkrd = bkrd;
			rt[a].background = background;
			rt[a].phase = a;
			rt[a].n_thread = n_thread;
			copy3f(spec_vector,rt[a].spec_vector);
			}
		
		if(n_thread<=1)
        RayTraceThread(rt);
		else 
        RayTraceSpawn(rt,n_thread);
    }
  }
  PRINTFD(FB_Ray)
    " RayRender: n_hit %d\n",n_hit
    ENDFD;

#if 1
  
  if(antialias)
    {	
		int		w2;
		int		wid_1, hgt_1;
      register unsigned int part;
      unsigned int acc;
		unsigned int z[16],zm[16];
		
		unsigned int *pSrc;
		unsigned int *pDst;
		const unsigned int m00FF=0x00FF,mFF00=0xFF00,mFFFF=0xFFFF;

		OrthoBusyFast(9,10);
		width	= width/2;
		height	= height/2;
		
		w2	= width * 2;
		hgt_1	= height - 1;
		wid_1	= width - 1;
		
		for(y = 1; y < hgt_1; y++)
		{
			pSrc  	= image + w2 * (y*2 - 1);
         pDst	= image_copy + width * y + 1;	/* Add 1 since x-wise inner loop (below) starts at 1 */
			
			for(x = 1; x < wid_1; x++)
			{
				aa	= 0;
		
				p	= pSrc + (x * 2);
			
            z[12]	= (*(p-1));
				z[4]	= (*(p));
				z[5]	= (*(p+1));
				z[13]	= (*(p+2));
				
				p	+= w2;
				
				z[6]	= (*(p-1));
				z[0]	= (*(p));
				z[1]	= (*(p+1));
				z[7]	= (*(p+2));
				
				p	+= w2;
				
				z[8]	= (*(p-1));
				z[2]	= (*(p));
				z[3]	= (*(p+1));
				z[9]	= (*(p+2));
				
				p	+= w2;
				
				z[14]	= (*(p-1));
				z[10]	= (*(p));
				z[11]	= (*(p+1));
				z[15]	= (*(p+2));

            for( a = 0; a < 16; a += 4 ) 
              {
                zm[a+0 ] = z[a+0] & mFFFF; /* move half to zm */
                zm[a+1 ] = z[a+1] & mFFFF;
                zm[a+2 ] = z[a+2] & mFFFF;
                zm[a+3 ] = z[a+3] & mFFFF;
                
                z[a   ]	>>= 16; /* keep rest in z */
                z[a+1 ]	>>= 16;
                z[a+2 ]	>>= 16;
                z[a+3 ]	>>= 16;
              }

            combine(part,z,m00FF);
            acc=(part<<16);
            combine(part,z,mFF00);
            acc+=(part<<16);
            combine(part,zm,m00FF);
            acc+=part;
            combine(part,zm,mFF00);
            acc+=part;
            *(pDst++) = acc;

			}
		}
		
		/* top and bottom edges */		
		for(y = 0; y < height; y = y + hgt_1)
		{
			pSrc  	= image + w2 * (y*2);

			pDst	= image_copy + width * y;
			
			for(x = 0; x < width; x++)
			{
				aa		= 0;
				p		= pSrc + (x*2);
				
				z[0]	= *p;
				z[1]	= *(p+1);
				
				p		+= w2;
				
				z[2]	= *(p);
				z[3]	= *(p+1);

            zm[0] = z[0] & mFFFF; /* move half to zm */
            zm[1] = z[1] & mFFFF;
            zm[2] = z[2] & mFFFF;
            zm[3] = z[3] & mFFFF;
            
            z[0]	>>= 16; /* keep rest in z */
            z[1]	>>= 16;
            z[2]	>>= 16;
            z[3]	>>= 16;

            edge_combine(part,z,m00FF);
            acc=(part<<16);
            edge_combine(part,z,mFF00);
            acc+=(part<<16);
            edge_combine(part,zm,m00FF);
            acc+=part;
            edge_combine(part,zm,mFF00);
            acc+=part;
            *(pDst++) = acc;
				
			}
		}

		/* left and right edges */
		for(y = 0; y < height; y++)
		{
			pSrc	= image + w2 * (y*2);
			pDst	= image_copy + width * y;
			
			for(x = 0; x < width; x += wid_1, pDst += wid_1 )
			{
				aa		= 0;
				p		= pSrc + (x*2);
				
				z[0]	= *p;
				z[1]	= *(p+1);
				
				p	+= w2;
				
				z[2]	= *p;
				z[3]	= *(p+1);
				
            zm[0] = z[0] & mFFFF; /* move half to zm */
            zm[1] = z[1] & mFFFF;
            zm[2] = z[2] & mFFFF;
            zm[3] = z[3] & mFFFF;
            
            z[0]	>>= 16; /* keep rest in z */
            z[1]	>>= 16;
            z[2]	>>= 16;
            z[3]	>>= 16;

            edge_combine(part,z,m00FF);
            acc=(part<<16);
            edge_combine(part,z,mFF00);
            acc+=(part<<16);
            edge_combine(part,zm,m00FF);
            acc+=part;
            edge_combine(part,zm,mFF00);
            acc+=part;

            *pDst = acc;

			}
		}
		
		/* FreeP(image); */
		image = image_copy;
	}
#else

	if(antialias)
	{	
		int		w2;
		int		isBigEndian	= I->BigEndian;
		int		wid_1, hgt_1;
      unsigned int z[16],zm[16],tot,za;

		unsigned int *pSrc;
		unsigned int *pDst;
		
		OrthoBusyFast(9,10);
		width	= width/2;
		height	= height/2;
		
		w2	= width * 2;
		hgt_1	= height - 1;
		wid_1	= width - 1;
		
		for(y = 1; y < hgt_1; y++)
		{
			pSrc  	= image + w2 * (y*2 - 1);
			pDst	= image_copy + width * y + 1;	/* Add 1 since x-wise inner loop (below) starts at 1 */
			
			for(x = 1; x < wid_1; x++)
			{
				aa	= 0;
		
				p	= pSrc + (x * 2);
				
			    z[12]	= (*(p-1));
				z[4]	= (*(p));
				z[5]	= (*(p+1));
				z[13]	= (*(p+2));
				
				p	+= w2;
				
				z[6]	= (*(p-1));
				z[0]	= (*(p));
				z[1]	= (*(p+1));
				z[7]	= (*(p+2));
				
				p	+= w2;
				
				z[8]	= (*(p-1));
				z[2]	= (*(p));
				z[3]	= (*(p+1));
				z[9]	= (*(p+2));
				
				p	+= w2;
				
				z[14]	= (*(p-1));
				z[10]	= (*(p));
				z[11]	= (*(p+1));
				z[15]	= (*(p+2));

				if(isBigEndian) 
				{ 
					for( a = 0; a < 16; a += 4 ) 
					{
                 /* unroll a bit */
						zm[a+0]	&= 0xFF;
						zm[a+1]	&= 0xFF;
						zm[a+2] &= 0xFF;
						zm[a+3] &= 0xFF;

						z[a+0]	>>= 8;
						z[a+1]	>>= 8;
						z[a+2]	>>= 8;
						z[a+3]	>>= 8;
					}
				}
				else
				{
					for(a=0;a<16;a++) 
						zm[a] = z[a] >> 24; /* copy alpha channel */
				}
		
				tot	= 0;
				for(a = 0; a < 16; a += 4) /* average alpha channel */
				{
					tot += zm[a+0] & 0xFF;
					tot += zm[a+1] & 0xFF;
					tot += zm[a+2] & 0xFF;
					tot += zm[a+3] & 0xFF;
				}
					
				tot += (zm[0] & 0xFF) << 2;
				tot += (zm[1] & 0xFF) << 2;
				tot += (zm[2] & 0xFF) << 2;
				tot += (zm[3] & 0xFF) << 2;

				za = (0xFF & (tot>>5));
		
				tot = 0;
				for(a = 0; a < 16; a += 4)
				{
					tot += z[a+0] & 0xFF;
					tot += z[a+1] & 0xFF;
					tot += z[a+2] & 0xFF;
					tot += z[a+3] & 0xFF;
				}
					
				tot += (z[0] & 0xFF)<<2;
				tot += (z[1] & 0xFF)<<2;
				tot += (z[2] & 0xFF)<<2;
				tot += (z[3] & 0xFF)<<2;
					
				aa	= aa | (0xFF & (tot>>5));
		
				tot = 0;
				for(a = 0; a < 16; a += 4)
				{
					tot += z[a+0] & 0xFF00;
					tot += z[a+1] & 0xFF00;
					tot += z[a+2] & 0xFF00;
					tot += z[a+3] & 0xFF00;
				}
					
				tot += (z[0] & 0xFF00)<<2;
				tot += (z[1] & 0xFF00)<<2;
				tot += (z[2] & 0xFF00)<<2;
				tot += (z[3] & 0xFF00)<<2;
					
				aa	= aa | (0xFF00 & (tot>>5));
		
				tot	= 0;
				for(a = 0; a < 16; a += 4)
				{
					tot	+= z[a+0] & 0xFF0000;
					tot	+= z[a+1] & 0xFF0000;
					tot	+= z[a+2] & 0xFF0000;
					tot	+= z[a+3] & 0xFF0000;
				}
				tot	+= (z[0] & 0xFF0000)<<2;
				tot	+= (z[1] & 0xFF0000)<<2;
				tot	+= (z[2] & 0xFF0000)<<2;
				tot	+= (z[3] & 0xFF0000)<<2;
					
				aa = aa | (0xFF0000&(tot>>5));			 
		
				if(isBigEndian) 
					aa = (aa<<8) | za;
				else 
					aa = aa | (za<<24);
		
				*(pDst++) = aa;		

			}
		}
		
		/* top and bottom edges */		
		for(y = 0; y < height; y = y + hgt_1)
		{
			pSrc  	= image + w2 * (y*2);
			pDst	= image_copy + width * y;
			
			for(x = 0; x < width; x++)
			{
				aa		= 0;
				p		= pSrc + (x*2);
				
				z[0]	= *p;
				z[1]	= *(p+1);
				
				p		+= w2;
				
				z[2]	= *(p);
				z[3]	= *(p+1);
				
				if(isBigEndian) 
				{ 
					for(a=0;a<4;a++) 
					{
						zm[a]=z[a]&0xFF;
						z[a]=z[a]>>8;
					}
				}
				else 
				{
					for(a=0;a<4;a++)
						zm[a]=z[a]>>24; /* copy alpha channel */
				}
				
				tot=0;
				for(a=0;a<4;a++)
					tot+=(zm[a]&0xFF);
				za=(0xFF&(tot>>2));
				
				tot=0;
				for(a=0;a<4;a++)
					tot+=(z[a]&0xFF);
				aa=aa|(0xFF&(tot>>2));
				
				tot=0;
				for(a=0;a<4;a++)
					tot+=(z[a]&0xFF00);
				aa=aa|(0xFF00&(tot>>2));
				
				tot=0;
				for(a=0;a<4;a++) 
					tot+=(z[a]&0xFF0000);
				aa=aa|(0xFF0000&(tot>>2));			 
				
				if(isBigEndian)
					aa=(aa<<8)|za;
				else
					aa=aa|(za<<24);

				*(pDst++) = aa;			 
			}
		}

		/* left and right edges */
		for(y = 0; y < height; y++)
		{
			pSrc	= image + w2 * (y*2);
			pDst	= image_copy + width * y;
			
			for(x = 0; x < width; x += wid_1, pDst += wid_1 )
			{
				aa		= 0;
				p		= pSrc + (x*2);
				
				z[0]	= *p;
				z[1]	= *(p+1);
				
				p	+= w2;
				
				z[2]	= *p;
				z[3]	= *(p+1);
				
				if(isBigEndian) 
				{ 
					for(a=0;a<4;a++) 
					{
						zm[a]	= z[a] & 0xFF;
						z[a]	= z[a] >> 8;
					}
				}
				else
				{
					for(a=0;a<4;a++)
						zm[a]	= z[a]>>24; /* copy alpha channel */
				}
				
				tot=0;
				for(a=0;a<4;a++)
					tot+=(zm[a]&0xFF);
				za=(0xFF&(tot>>2));
				
				tot=0;
				for(a=0;a<4;a++)
					tot+=(z[a]&0xFF);
				aa=aa|(0xFF&(tot>>2));
				
				tot=0;
				for(a=0;a<4;a++)
					tot+=(z[a]&0xFF00);
				aa=aa|(0xFF00&(tot>>2));
				
				tot=0;
				for(a=0;a<4;a++)
					tot+=(z[a]&0xFF0000);
				aa=aa|(0xFF0000&(tot>>2));			 
				
				if(isBigEndian)
					aa=(aa<<8)|za;
				else
					aa=aa|(za<<24);

				*pDst = aa;
			}
		}
		/* FreeP(image); */
		image = image_copy;
	}
#endif

}

void RayRenderColorTable(CRay *I,int width,int height,int *image)
{
  int x,y;
  unsigned int r=0,g=0,b=0;
  unsigned int *pixel,mask,*p;

  if(I->BigEndian)
    mask = 0x000000FF;
  else
    mask = 0xFF000000;

  p=(unsigned int*)image; 
  for(x=0;x<width;x++)
    for(y=0;y<height;y++)
      *(p++)=mask;
  
  if((width>=512)&&(height>=512)) {
    
    
    for(y=0;y<512;y++) 
      for(x=0;x<512;x++)        
        {
          pixel = (unsigned int*) (image+((width)*y)+x);
          if(I->BigEndian) {
            *(pixel)=
              mask|(r<<24)|(g<<16)|(b<<8);
          } else {
            *(pixel)=
              mask|(b<<16)|(g<<8)|r;
          }
          b = b + 4;
          if(!(0xFF&b)) { 
            b=0;
            g=g+4;
            if(!(0xFF&g)) {           
              g=0;
              r=r+4;
            }
          }
        }
  }
}
/*========================================================================*/
void RayTexture(CRay *I,int mode,float *v)
{
  I->Texture=mode;
  if(v) 
    copy3f(v,I->TextureParam);
}
/*========================================================================*/
void RayTransparentf(CRay *I,float v)
{
  I->Trans=v;
}
/*========================================================================*/
void RayColor3fv(CRay *I,float *v)
{
  I->CurColor[0]=(*v++);
  I->CurColor[1]=(*v++);
  I->CurColor[2]=(*v++);
}
/*========================================================================*/
void RaySphere3fv(CRay *I,float *v,float r)
{
  CPrimitive *p;

  float *vv;

  VLACheck(I->Primitive,CPrimitive,I->NPrimitive);
  p = I->Primitive+I->NPrimitive;

  p->type = cPrimSphere;
  p->r1=r;
  p->trans=I->Trans;
  p->texture=I->Texture;
  copy3f(I->TextureParam,p->texture_param);
  vv=p->v1;
  (*vv++)=(*v++);
  (*vv++)=(*v++);
  (*vv++)=(*v++);


  vv=p->c1;
  v=I->CurColor;
  (*vv++)=(*v++);
  (*vv++)=(*v++);
  (*vv++)=(*v++);

  if(I->TTTFlag) {
    transformTTT44f3f(I->TTT,p->v1,p->v1);
  }

  if(I->Context) {
    RayApplyContextToVertex(I,p->v1);
  }

  I->NPrimitive++;
}
/*========================================================================*/
void RayCylinder3fv(CRay *I,float *v1,float *v2,float r,float *c1,float *c2)
{
  CPrimitive *p;

  float *vv;

  VLACheck(I->Primitive,CPrimitive,I->NPrimitive);
  p = I->Primitive+I->NPrimitive;

  p->type = cPrimCylinder;
  p->r1=r;
  p->trans=I->Trans;
  p->texture=I->Texture;
  p->cap1=cCylCapFlat;
  p->cap2=cCylCapFlat;
  copy3f(I->TextureParam,p->texture_param);

  vv=p->v1;
  (*vv++)=(*v1++);
  (*vv++)=(*v1++);
  (*vv++)=(*v1++);
  vv=p->v2;
  (*vv++)=(*v2++);
  (*vv++)=(*v2++);
  (*vv++)=(*v2++);

  if(I->TTTFlag) {
    transformTTT44f3f(I->TTT,p->v1,p->v1);
    transformTTT44f3f(I->TTT,p->v2,p->v2);
  }

  if(I->Context) {
    RayApplyContextToVertex(I,p->v1);
    RayApplyContextToVertex(I,p->v2);
  }

  vv=p->c1;
  (*vv++)=(*c1++);
  (*vv++)=(*c1++);
  (*vv++)=(*c1++);
  vv=p->c2;
  (*vv++)=(*c2++);
  (*vv++)=(*c2++);
  (*vv++)=(*c2++);

  I->NPrimitive++;
}
/*========================================================================*/
void RayCustomCylinder3fv(CRay *I,float *v1,float *v2,float r,
                          float *c1,float *c2,int cap1,int cap2)
{
  CPrimitive *p;

  float *vv;

  VLACheck(I->Primitive,CPrimitive,I->NPrimitive);
  p = I->Primitive+I->NPrimitive;

  p->type = cPrimCylinder;
  p->r1=r;
  p->trans=I->Trans;
  p->texture=I->Texture;
  p->cap1=cap1;
  p->cap2=cap2;
  copy3f(I->TextureParam,p->texture_param);

  vv=p->v1;
  (*vv++)=(*v1++);
  (*vv++)=(*v1++);
  (*vv++)=(*v1++);
  vv=p->v2;
  (*vv++)=(*v2++);
  (*vv++)=(*v2++);
  (*vv++)=(*v2++);

  if(I->TTTFlag) {
    transformTTT44f3f(I->TTT,p->v1,p->v1);
    transformTTT44f3f(I->TTT,p->v2,p->v2);
  }

  if(I->Context) {
    RayApplyContextToVertex(I,p->v1);
    RayApplyContextToVertex(I,p->v2);
  }

  vv=p->c1;
  (*vv++)=(*c1++);
  (*vv++)=(*c1++);
  (*vv++)=(*c1++);
  vv=p->c2;
  (*vv++)=(*c2++);
  (*vv++)=(*c2++);
  (*vv++)=(*c2++);

  I->NPrimitive++;
}
/*========================================================================*/
void RaySausage3fv(CRay *I,float *v1,float *v2,float r,float *c1,float *c2)
{
  CPrimitive *p;

  float *vv;

  VLACheck(I->Primitive,CPrimitive,I->NPrimitive);
  p = I->Primitive+I->NPrimitive;

  p->type = cPrimSausage;
  p->r1=r;
  p->trans=I->Trans;
  p->texture=I->Texture;
  copy3f(I->TextureParam,p->texture_param);

  vv=p->v1;
  (*vv++)=(*v1++);
  (*vv++)=(*v1++);
  (*vv++)=(*v1++);
  vv=p->v2;
  (*vv++)=(*v2++);
  (*vv++)=(*v2++);
  (*vv++)=(*v2++);

  if(I->TTTFlag) {
    transformTTT44f3f(I->TTT,p->v1,p->v1);
    transformTTT44f3f(I->TTT,p->v2,p->v2);
  }

  if(I->Context) {
    RayApplyContextToVertex(I,p->v1);
    RayApplyContextToVertex(I,p->v2);
  }

  vv=p->c1;
  (*vv++)=(*c1++);
  (*vv++)=(*c1++);
  (*vv++)=(*c1++);
  vv=p->c2;
  (*vv++)=(*c2++);
  (*vv++)=(*c2++);
  (*vv++)=(*c2++);

  I->NPrimitive++;
}
/*========================================================================*/
void RayTriangle3fv(CRay *I,
						  float *v1,float *v2,float *v3,
						  float *n1,float *n2,float *n3,
						  float *c1,float *c2,float *c3)
{
  CPrimitive *p;

  float *vv;
  float n0[3],nx[3],s1[3],s2[3],s3[3];
  float l1,l2,l3;

  /*  dump3f(v1," v1");
  dump3f(v2," v2");
  dump3f(v3," v3");
  dump3f(n1," n1");
  dump3f(n2," n2");
  dump3f(n3," n3");
  dump3f(c1," c1");
  dump3f(c2," c2");
  dump3f(c3," c3");*/
  VLACheck(I->Primitive,CPrimitive,I->NPrimitive);
  p = I->Primitive+I->NPrimitive;

  p->type = cPrimTriangle;
  p->trans=I->Trans;
  p->texture=I->Texture;
  copy3f(I->TextureParam,p->texture_param);

  /* determine exact triangle normal */
  add3f(n1,n2,nx);
  add3f(n3,nx,nx);
  subtract3f(v1,v2,s1);
  subtract3f(v3,v2,s2);
  subtract3f(v1,v3,s3);
  cross_product3f(s1,s2,n0);
  if((fabs(n0[0])<RAY_SMALL)&&
	  (fabs(n0[1])<RAY_SMALL)&&
	  (fabs(n0[2])<RAY_SMALL))
	 {copy3f(nx,n0);} /* fall-back */
  else if(dot_product3f(n0,nx)<0)
	 invert3f(n0);
  normalize3f(n0);

  vv=p->n0;
  (*vv++)=n0[0];
  (*vv++)=n0[1];
  (*vv++)=n0[2];

  /* determine maximum distance from vertex to point */
  l1=(float)length3f(s1);
  l2=(float)length3f(s2);
  l3=(float)length3f(s3);
  if(l2>l1) { if(l3>l2)	l1=l3; else	l1=l2;  }
  /* store cutoff distance */

  p->r1=l1*0.6F;

  /*  if(l1>20) {
		printf("%8.3f\n",l1);
		printf("%8.3f %8.3f %8.3f\n",s1[0],s1[1],s1[2]);
		printf("%8.3f %8.3f %8.3f\n",s2[0],s2[1],s2[2]);
		printf("%8.3f %8.3f %8.3f\n",s3[0],s3[1],s3[2]);
		}*/

  vv=p->v1;
  (*vv++)=(*v1++);
  (*vv++)=(*v1++);
  (*vv++)=(*v1++);
  vv=p->v2;
  (*vv++)=(*v2++);
  (*vv++)=(*v2++);
  (*vv++)=(*v2++);
  vv=p->v3;
  (*vv++)=(*v3++);
  (*vv++)=(*v3++);
  (*vv++)=(*v3++);

  vv=p->c1;
  (*vv++)=(*c1++);
  (*vv++)=(*c1++);
  (*vv++)=(*c1++);
  vv=p->c2;
  (*vv++)=(*c2++);
  (*vv++)=(*c2++);
  (*vv++)=(*c2++);
  vv=p->c3;
  (*vv++)=(*c3++);
  (*vv++)=(*c3++);
  (*vv++)=(*c3++);

  vv=p->n1;
  (*vv++)=(*n1++);
  (*vv++)=(*n1++);
  (*vv++)=(*n1++);
  vv=p->n2;
  (*vv++)=(*n2++);
  (*vv++)=(*n2++);
  (*vv++)=(*n2++);
  vv=p->n3;
  (*vv++)=(*n3++);
  (*vv++)=(*n3++);
  (*vv++)=(*n3++);


  if(I->TTTFlag) {
    transformTTT44f3f(I->TTT,p->v1,p->v1);
    transformTTT44f3f(I->TTT,p->v2,p->v2);
    transformTTT44f3f(I->TTT,p->v3,p->v3);
    transform_normalTTT44f3f(I->TTT,p->n0,p->n0);
    transform_normalTTT44f3f(I->TTT,p->n1,p->n1);
    transform_normalTTT44f3f(I->TTT,p->n2,p->n2);
    transform_normalTTT44f3f(I->TTT,p->n3,p->n3);
  }

  if(I->Context) {
    RayApplyContextToVertex(I,p->v1);
    RayApplyContextToVertex(I,p->v2);
    RayApplyContextToVertex(I,p->v3);
    RayApplyContextToNormal(I,p->n0);
    RayApplyContextToNormal(I,p->n1);
    RayApplyContextToNormal(I,p->n2);
    RayApplyContextToNormal(I,p->n3);
  }

  I->NPrimitive++;

}
/*========================================================================*/
CRay *RayNew(void)
{
  unsigned int test;
  unsigned char *testPtr;
  int a;

  OOAlloc(CRay);
  
  test = 0xFF000000;
  testPtr = (unsigned char*)&test;
  I->BigEndian = (*testPtr)&&1;
  I->Trans=0.0F;
  I->Texture=0;
  I->TTTFlag=false;
  zero3f(I->TextureParam);

  PRINTFB(FB_Ray,FB_Blather) 
    " RayNew: BigEndian = %d\n",I->BigEndian
    ENDFB;

  I->Basis=VLAlloc(CBasis,10);
  BasisInit(I->Basis);
  BasisInit(I->Basis+1);
  I->Vert2Prim=VLAlloc(int,1);
  I->NBasis=2;
  I->Primitive=NULL;
  I->NPrimitive=0;
  I->fColor3fv=RayColor3fv;
  I->fSphere3fv=RaySphere3fv;
  I->fCylinder3fv=RayCylinder3fv;
  I->fCustomCylinder3fv=RayCustomCylinder3fv;
  I->fSausage3fv=RaySausage3fv;
  I->fTriangle3fv=RayTriangle3fv;
  I->fTexture=RayTexture;
  I->fTransparentf=RayTransparentf;
  if(!RandomFlag) {
    for(a=0;a<256;a++) {
      Random[a]=(float)((rand()/(1.0+RAND_MAX))-0.5);
    }
    RandomFlag=1;
  }

  return(I);
}
/*========================================================================*/
void RayPrepare(CRay *I,float v0,float v1,float v2,
                float v3,float v4,float v5,
                float *mat,float aspRat,int ray_width)
	  /*prepare for vertex calls */
{
  int a;
  if(!I->Primitive) 
	 I->Primitive=VLAlloc(CPrimitive,100);  
  if(!I->Vert2Prim) 
	 I->Vert2Prim=VLAlloc(int,100);
  I->Volume[0]=v0;
  I->Volume[1]=v1;
  I->Volume[2]=v2;
  I->Volume[3]=v3;
  I->Volume[4]=v4;
  I->Volume[5]=v5;
  I->Range[0]=I->Volume[1]-I->Volume[0];
  I->Range[1]=I->Volume[3]-I->Volume[2];
  I->Range[2]=I->Volume[5]-I->Volume[4];
  I->AspRatio=aspRat;

  if(mat)  
    for(a=0;a<16;a++)
      I->ModelView[a]=mat[a];
  else {
    for(a=0;a<16;a++)
      I->ModelView[a]=0.0F;
    for(a=0;a<3;a++)
      I->ModelView[a*5]=1.0F;
  }
  if(ray_width)
    I->PixelRadius = ((float)I->Range[0])/ray_width;
  else
    I->PixelRadius = 0.15F;
}
/*========================================================================*/

void RaySetTTT(CRay *I,int flag,float *ttt)
{
  I->TTTFlag=flag;
  if(flag) {
    UtilCopyMem(I->TTT,ttt,sizeof(float)*16);
  }
}

/*========================================================================*/
void RayRelease(CRay *I)
{
  int a;

  for(a=0;a<I->NBasis;a++) {
	 BasisFinish(&I->Basis[a]);
  }
  I->NBasis=0;
  VLAFreeP(I->Primitive);
  VLAFreeP(I->Vert2Prim);
}
/*========================================================================*/
void RayFree(CRay *I)
{
  RayRelease(I);
  VLAFreeP(I->Basis);
  VLAFreeP(I->Vert2Prim);
  OOFreeP(I);
}
/*========================================================================*/

void RayApplyMatrix33( unsigned int n, float3 *q, const float m[16],
                          float3 *p )
{
   {
      unsigned int i;
      float m0 = m[0],  m4 = m[4],  m8 = m[8],  m12 = m[12];
      float m1 = m[1],  m5 = m[5],  m9 = m[9],  m13 = m[13];
      float m2 = m[2],  m6 = m[6],  m10 = m[10],  m14 = m[14];
      for (i=0;i<n;i++) {
         float p0 = p[i][0], p1 = p[i][1], p2 = p[i][2];
         q[i][0] = m0 * p0 + m4  * p1 + m8 * p2 + m12;
         q[i][1] = m1 * p0 + m5  * p1 + m9 * p2 + m13;
         q[i][2] = m2 * p0 + m6 * p1 + m10 * p2 + m14;
      }
   }
}

void RayApplyMatrixInverse33( unsigned int n, float3 *q, const float m[16],
                          float3 *p )
{
   {
      unsigned int i;
      float m0 = m[0],  m4 = m[4],  m8 = m[8],  m12 = m[12];
      float m1 = m[1],  m5 = m[5],  m9 = m[9],  m13 = m[13];
      float m2 = m[2],  m6 = m[6],  m10 = m[10],  m14 = m[14];
      for (i=0;i<n;i++) {
         float p0 = p[i][0]-m12, p1 = p[i][1]-m13, p2 = p[i][2]-m14;
         q[i][0] = m0 * p0 + m1  * p1 + m2 * p2;
         q[i][1] = m4 * p0 + m5  * p1 + m6 * p2;
         q[i][2] = m8 * p0 + m9 * p1 + m10 * p2;
      }
   }
}

void RayTransformNormals33( unsigned int n, float3 *q, const float m[16],float3 *p )
{
  unsigned int i;
  float m0 = m[0],  m4 = m[4],  m8 = m[8];
  float m1 = m[1],  m5 = m[5],  m9 = m[9];
  float m2 = m[2],  m6 = m[6],  m10 = m[10];
  for (i=0;i<n;i++) {
    float p0 = p[i][0], p1 = p[i][1], p2 = p[i][2];
    q[i][0] = m0 * p0 + m4  * p1 + m8 * p2;
    q[i][1] = m1 * p0 + m5  * p1 + m9 * p2;
    q[i][2] = m2 * p0 + m6 * p1 + m10 * p2;
  }
  for (i=0;i<n;i++) { /* renormalize - can we do this to the matrix instead? */
    normalize3f(q[i]);
  }
}

void RayTransformInverseNormals33( unsigned int n, float3 *q, const float m[16],float3 *p )
{
  unsigned int i;
  float m0 = m[0],  m4 = m[4],  m8 = m[8];
  float m1 = m[1],  m5 = m[5],  m9 = m[9];
  float m2 = m[2],  m6 = m[6],  m10 = m[10];
  for (i=0;i<n;i++) {
    float p0 = p[i][0], p1 = p[i][1], p2 = p[i][2];
    q[i][0] = m0 * p0 + m1  * p1 + m2 * p2;
    q[i][1] = m4 * p0 + m5  * p1 + m6 * p2;
    q[i][2] = m8 * p0 + m9 * p1 + m10 * p2;
  }
  for (i=0;i<n;i++) { /* renormalize - can we do this to the matrix instead? */
    normalize3f(q[i]);
  }
}



