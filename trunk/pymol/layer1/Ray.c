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
#include"MemoryCache.h"

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

#define MAX_RAY_THREADS 12

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
void RayTransformBasis(CRay *I,CBasis *B,int group_id);

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

static void fill(unsigned int *buffer, unsigned int value,unsigned int cnt)
{
  while(cnt&0xFFFFFF80) {
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
  
  VLACacheCheck(basis->Vertex,float,3*nVert,0,cCache_basis_vertex);
  VLACacheCheck(basis->Radius,float,nVert,0,cCache_basis_radius);
  VLACacheCheck(basis->Radius2,float,nVert,0,cCache_basis_radius2);
  VLACacheCheck(basis->Vert2Normal,int,nVert,0,cCache_basis_vert2normal);
  VLACacheCheck(basis->Normal,float,3*nNorm,0,cCache_basis_normal);

  VLACacheCheck(I->Vert2Prim,int,nVert,0,cCache_ray_vert2prim);

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
static void RayComputeBox(CRay *I)
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

  float xmin=0.0F,ymin=0.0F,xmax=0.0F,ymax=0.0F;
  float xp,xm,yp,ym;

  float *v,r;
  float vt[3];
  const float _0 = 0.0F;
  int a;

  basis1 = I->Basis+1;
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
  I->min_box[0] = xmin;
  I->min_box[1] = ymin;
  I->max_box[0] = xmax;
  I->max_box[1] = ymax;
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
  
  VLACacheCheck(basis1->Vertex,float,3*basis0->NVertex,1,cCache_basis_vertex);
  VLACacheCheck(basis1->Normal,float,3*basis0->NNormal,1,cCache_basis_normal);
  VLACacheCheck(basis1->Precomp,float,3*basis0->NNormal,1,cCache_basis_precomp);
  VLACacheCheck(basis1->Vert2Normal,int,basis0->NVertex,1,cCache_basis_vert2normal);
  VLACacheCheck(basis1->Radius,float,basis0->NVertex,1,cCache_basis_radius);
  VLACacheCheck(basis1->Radius2,float,basis0->NVertex,1,cCache_basis_radius2);
  
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
    switch(prm->type) {
    case cPrimTriangle:
        BasisTrianglePrecompute(basis1->Vertex+prm->vert*3,
                                basis1->Vertex+prm->vert*3+3,
                                basis1->Vertex+prm->vert*3+6,
                                basis1->Precomp+basis1->Vert2Normal[prm->vert]*3);
        v0 = basis1->Normal + (basis1->Vert2Normal[prm->vert]*3 + 3);
        prm->cull = backface_cull&&((v0[2]<0.0F)&&(v0[5]<0.0F)&&(v0[8]<0.0F));
        break;
    case cPrimSausage:
    case cPrimCylinder:
      BasisCylinderSausagePrecompute(basis1->Normal + basis1->Vert2Normal[prm->vert]*3,
                                     basis1->Precomp + basis1->Vert2Normal[prm->vert]*3);
      break;
      
	 }
  }
}
/*========================================================================*/
void RayTransformBasis(CRay *I,CBasis *basis1,int group_id)
{
  CBasis *basis0;
  int a;
  float *v0,*v1;
  CPrimitive *prm;

  basis0 = I->Basis+1;

  VLACacheCheck(basis1->Vertex,float,3*basis0->NVertex,group_id,cCache_basis_vertex);
  VLACacheCheck(basis1->Normal,float,3*basis0->NNormal,group_id,cCache_basis_normal);
  VLACacheCheck(basis1->Precomp,float,3*basis0->NNormal,group_id,cCache_basis_precomp);
  VLACacheCheck(basis1->Vert2Normal,int,basis0->NVertex,group_id,cCache_basis_vert2normal);
  VLACacheCheck(basis1->Radius,float,basis0->NVertex,group_id,cCache_basis_radius);
  VLACacheCheck(basis1->Radius2,float,basis0->NVertex,group_id,cCache_basis_radius2);
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
    switch(prm->type) {
    case cPrimTriangle:
        BasisTrianglePrecompute(basis1->Vertex+prm->vert*3,
                                basis1->Vertex+prm->vert*3+3,
                                basis1->Vertex+prm->vert*3+6,
                                basis1->Precomp+basis1->Vert2Normal[prm->vert]*3);
        break;
    case cPrimSausage:
    case cPrimCylinder:
      BasisCylinderSausagePrecompute(basis1->Normal + basis1->Vert2Normal[prm->vert]*3,
                                     basis1->Precomp + basis1->Vert2Normal[prm->vert]*3);
      break;
      
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
  if(fog<0.0F)
    fog = SettingGet(cSetting_depth_cue);
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

  PRINTFB(FB_Ray,FB_Blather)
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

static void RayAntiSpawn(CRayAntiThreadInfo *Thread,int n_thread)
{
  int blocked;
  PyObject *info_list;
  int a;
  blocked = PAutoBlock();

  PRINTFB(FB_Ray,FB_Blather)
    " Ray: antialiasing with %d threads...\n",n_thread
  ENDFB;
  info_list = PyList_New(n_thread);
  for(a=0;a<n_thread;a++) {
    PyList_SetItem(info_list,a,PyCObject_FromVoidPtr(Thread+a,NULL));
  }
  PyObject_CallMethod(P_cmd,"_ray_anti_spawn","O",info_list);
  Py_DECREF(info_list);
  PAutoUnblock(blocked);
}

int RayHashThread(CRayHashThreadInfo *T)
{
  BasisMakeMap(T->basis,T->vert2prim,T->prim,T->clipBox,T->phase,cCache_ray_map);

  /* utilize a little extra wasted CPU time in thread 0 which computes the smaller map... */

  if(!T->phase) { 
    fill(T->image,T->background,T->bytes);
    RayComputeBox(T->ray);
  }
  return 1;
}

static void RayTraceSpawn(CRayThreadInfo *Thread,int n_thread)
{
  int blocked;
  PyObject *info_list;
  int a;
  blocked = PAutoBlock();

  PRINTFB(FB_Ray,FB_Blather)
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

static int find_edge(unsigned int *ptr,unsigned int width,int threshold)
{
  unsigned int shift = 0;
  int compare[9];
  int sum[9] = {0,0,0,0,0,0,0,0};
  int current;
  int a;

  compare[0] = (signed int)*(ptr);
  compare[1] = (signed int)*(ptr-1);
  compare[2] = (signed int)*(ptr+1);
  compare[3] = (signed int)*(ptr-width);
  compare[4] = (signed int)*(ptr+width);
  compare[5] = (signed int)*(ptr-width-1);
  compare[6] = (signed int)*(ptr+width-1);
  compare[7] = (signed int)*(ptr-width+1);
  compare[8] = (signed int)*(ptr+width+1);
  
  for(a=0;a<4;a++) {
    current = ((compare[0]>>shift)&0xFF);
    sum[1] += abs(current - ((compare[1]>>shift)&0xFF));
    sum[2] += abs(current - ((compare[2]>>shift)&0xFF));
    if(sum[1]>=threshold) return 1;
    sum[3] += abs(current - ((compare[3]>>shift)&0xFF));
    if(sum[2]>=threshold) return 1;
    sum[4] += abs(current - ((compare[4]>>shift)&0xFF));
    if(sum[3]>=threshold) return 1;
    sum[5] += abs(current - ((compare[5]>>shift)&0xFF));
    if(sum[4]>=threshold) return 1;
    sum[6] += abs(current - ((compare[6]>>shift)&0xFF));
    if(sum[5]>=threshold) return 1;
    sum[7] += abs(current - ((compare[7]>>shift)&0xFF));
    if(sum[6]>=threshold) return 1;
    sum[8] += abs(current - ((compare[8]>>shift)&0xFF));
    if(sum[7]>=threshold) return 1;
    if(sum[8]>=threshold) return 1;
    shift+=8;
  }
  return 0;
}

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
	int texture_save;
	float		settingPower, settingReflectPower,settingSpecPower,settingSpecReflect, _0, _1, _p5, _255, _persistLimit, _inv3;
	float		invHgt, invFrontMinusBack, inv1minusFogStart,invWdth,invHgtRange;
	register float       invWdthRange,vol0;
	float       vol2;
	CBasis      *bp1,*bp2;
	int render_height;
	int offset=0;
   BasisCallRec SceneCall,ShadeCall;
   float border_offset;
   int edge_sampling = false;
   unsigned int edge_avg[4];
   int edge_cnt=0;
   float base[2];
   float edge_width = 0.35356F;
   float edge_height = 0.35356F;
   float trans_spec_cut,trans_spec_scale;

   float red_blend=0.0F;
   float blue_blend=0.0F;
   float green_blend=0.0F;
   int blend_colors;

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
   trans_spec_cut = SettingGet(cSetting_ray_transparency_spec_cut);
   blend_colors    = (int)SettingGet(cSetting_ray_blend_colors);
   if(blend_colors) {
     red_blend = SettingGet(cSetting_ray_blend_red);
     green_blend = SettingGet(cSetting_ray_blend_green);
     blue_blend = SettingGet(cSetting_ray_blend_blue);
   }

   if(trans_spec_cut<_1)
     trans_spec_scale = _1/(_1-trans_spec_cut);
   else
     trans_spec_scale = _0;

	/* COOP */
	settingPower		= SettingGet(cSetting_power);
	settingReflectPower	= SettingGet(cSetting_reflect_power);
	settingSpecPower	= SettingGet(cSetting_spec_power);

	settingSpecReflect	= SettingGet(cSetting_spec_reflect);
   if(settingSpecReflect>1.0F) settingSpecReflect = 1.0F;
	if(SettingGet(cSetting_specular)<R_SMALL4) {
     settingSpecReflect = 0.0F;
   }
    
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
   if(fog<0.0F)
     fog = SettingGet(cSetting_depth_cue);
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
	
   if(T->border) {
     invHgt				= _1 / (float) (T->height-(3.0F+T->border));
     invWdth             = _1 / (float) (T->width-(3.0F+T->border));
   } else {

     invHgt				= _1 / (float) (T->height);
     invWdth             = _1 / (float) (T->width);
   }

	invFrontMinusBack	= _1 / (T->front - T->back);
	invWdthRange        = invWdth * I->Range[0];
   invHgtRange         = invHgt * I->Range[1];
   edge_width *= invWdthRange;
   edge_height *= invHgtRange;
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

   SceneCall.Basis = I->Basis + 1;
   SceneCall.rr = &r1;
   SceneCall.vert2prim = I->Vert2Prim;
   SceneCall.prim = I->Primitive;
   SceneCall.shadow = false;
   SceneCall.back = T->back;
   SceneCall.trans_shadows = trans_shadows;
   SceneCall.check_interior = (interior_color >= 0);
	MapCacheInit(&SceneCall.cache,I->Basis[1].Map,T->phase,cCache_map_scene_cache);

   if(shadows&&(I->NBasis>1)) {
     ShadeCall.Basis = I->Basis + 2;
     ShadeCall.rr = &r2;
     ShadeCall.vert2prim = I->Vert2Prim;
     ShadeCall.prim = I->Primitive;
     ShadeCall.shadow = true;
     ShadeCall.front = _0;
     ShadeCall.back = _0;
     ShadeCall.excl_trans = _0;
     ShadeCall.trans_shadows = trans_shadows;
     ShadeCall.check_interior = false;
     MapCacheInit(&ShadeCall.cache,I->Basis[2].Map,T->phase,cCache_map_shadow_cache);     
   }

   if(T->border) {
     border_offset = -1.50F+T->border/2.0F;
   } else {
     border_offset = 0.0F;
   }
	for(yy = T->y_start; (yy < T->y_stop); yy++)
	{
		
      y = T->y_start + ((yy-T->y_start) + offset) % ( render_height); /* make sure threads write to different pages */

		if((!T->phase)&&!(yy & 0xF)) { /* don't slow down rendering too much */
        if(T->edging_cutoff) {
          if(T->edging) {
            OrthoBusyFast((int)(2.5F*T->height/3 + 0.5F*y),4*T->height/3); 
          } else {
            OrthoBusyFast((int)(T->height/3 + 0.5F*y),4*T->height/3); 
          }
        } else {
			OrthoBusyFast(T->height/3 + y,4*T->height/3); 
        }
      }
		pixel = T->image + (T->width * y) + T->x_start;
	
		if((y % T->n_thread) == T->phase)	/* this is my scan line */
		{	
        r1.base[1]	= ((y+0.5F+border_offset) * invHgtRange) + vol2;
			
			for(x = T->x_start; (x < T->x_stop); x++)
			{
				
           r1.base[0]	= (((x+0.5F+border_offset)) * invWdthRange)  + vol0;

            while(1) {
              if(T->edging) {
                if(!edge_sampling) {
                  if(x&&y&&(x<(T->width-1))&&(y<(T->height-1))) { /* not on the edge... */
                    if(find_edge(T->edging + (pixel - T->image),
                                 T->width, T->edging_cutoff)) {
                      unsigned int value;
                      edge_cnt = 1;
                      edge_sampling = true;
                      value = *pixel;
                      edge_avg[0] = value&0xFF;
                      edge_avg[1] = (value>>8)&0xFF;
                      edge_avg[2] = (value>>16)&0xFF;
                      edge_avg[3] = (value>>24)&0xFF;
                      base[0]=r1.base[0];
                      base[1]=r1.base[1];
                    }
                  }
                }
                if(edge_sampling) {
                  if(edge_cnt==5) {
                    edge_sampling=false;
                    /* done with edging, so store averaged value */

                    edge_avg[0]/=edge_cnt;
                    edge_avg[1]/=edge_cnt;
                    edge_avg[2]/=edge_cnt;
                    edge_avg[3]/=edge_cnt;

                    *pixel = (((edge_avg[0]&0xFF)    )|
                              ((edge_avg[1]&0xFF)<<8 )|
                              ((edge_avg[2]&0xFF)<<16)|
                              ((edge_avg[3]&0xFF)<<24));

                    /**pixel = 0xFFFFFFFF-*pixel;*/
                    /* restore X,Y coordinates */
                    r1.base[0]=base[0];
                    r1.base[1]=base[1];

                    /**pixel = 0xFF00FFFF;*/
                  } else {
                    *pixel = T->background;
                    switch(edge_cnt) {
                    case 1:
                      r1.base[0] = base[0]+edge_width;
                      r1.base[1] = base[1]+edge_height;
                      break;
                    case 2:
                      r1.base[0] = base[0]+edge_width;
                      r1.base[1] = base[1]-edge_height;
                      break;
                    case 3:
                      r1.base[0] = base[0]-edge_width;
                      r1.base[1] = base[1]+edge_height;
                      break;
                    case 4:
                      r1.base[0] = base[0]-edge_width;
                      r1.base[1] = base[1]-edge_height;
                      break;
                    }
                  }
                }
                if(!edge_sampling) /* not oversampling this edge or already done... */
                  break;
              }
              
              exclude		= -1;
              persist			= _1;
              first_excess	= _0;
              excl_trans		= _0;
              pass			= 0;
              new_front		= T->front;
              while((persist > _persistLimit) && (pass < 25))
                {
                  pixel_flag		= false;
                  
                  SceneCall.except = exclude;
                  SceneCall.front = new_front;
                  SceneCall.excl_trans = excl_trans;
                  
#if SPLIT_BASIS
                  i	= BasisHitNoShadow( &SceneCall );
#else
                  i	= BasisHit( &SceneCall );
#endif               
                  interior_flag = SceneCall.interior_flag;
                  
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
                      
                      if(shadows && ((!interior_flag)||(interior_shadows))) 
                        {
                          matrix_transform33f3f(bp2->Matrix,r1.impact,r2.base);
                          r2.base[2]-=shadow_fudge;
                          ShadeCall.except = i;
#if SPLIT_BASIS
                          if(BasisHitShadow(&ShadeCall) > -1)
                            lit	= (float) pow(r2.prim->trans, _p5);
#else
                          if(BasisHit(&ShadeCall) > -1)
                            lit	= (float) pow(r2.prim->trans, _p5);
#endif
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
                              fc[3] = ffact1m*(_1 - r1.prim->trans);
                            }
                          
                          if(!pass) {
                            if(r1.prim->trans<trans_spec_cut) {
                              first_excess = excess*ffact1m*ray_trans_spec;
                            } else {
                              first_excess = excess*ffact1m*ray_trans_spec*
                                trans_spec_scale*(_1 - r1.prim->trans);
                            }
                          } else {
                              fc[0]+=first_excess;
                              fc[1]+=first_excess;
                              fc[2]+=first_excess;
                            }
                        }
                      else 
                        {
                          if(!pass) {
                            if(r1.prim->trans<trans_spec_cut) {
                              first_excess = excess*ray_trans_spec;
                            } else {
                              first_excess = excess*ray_trans_spec*
                                trans_spec_scale*(_1 - r1.prim->trans);
                            }
                          } else {
                              fc[0]	+= first_excess;
                              fc[1]	+= first_excess;
                              fc[2]	+= first_excess;
                            }
                          if(opaque_back) {
                            fc[3]	= _1;
                          } else {
                            fc[3] = _1 - r1.prim->trans;
                          }
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
                      
                      if(!opaque_back) {
                        if(i<0) { /* hit nothing -- so don't blend alpha*/
                          fc[0] = (float)(0xFF&(last_pixel>>24));
                          fc[1] = (float)(0xFF&(last_pixel>>16));
                          fc[2] = (float)(0xFF&(last_pixel>>8));
                          fc[3] = (float)(0xFF&(last_pixel));
                        } else { /* hit something -- so keep blend and compute cumulative alpha*/
                          if(i>=0) { /* make sure opaque objects get opaque alpha*/
                            float o1,o2;
                            float m;
                            
                            if(I->BigEndian) {
                              o1 = (float)(0xFF&(last_pixel))/255.0F;
                              o2 = (float)(0xFF&(*pixel))/255.0F;
                            } else {
                              o1 = (float)(0xFF&(last_pixel>>24))/255.0F;
                              o2 = (float)(0xFF&((*pixel)>>24))/255.0F;
                            }
                            
                            if(o1<o2) { /* make sure o1 is largest opacity*/
                              m = o1;
                              o1 = o2;
                              o2 = m;
                            }
                            m = o1 + (1.0F - o1) * o2;
                            if(I->BigEndian) {
                              fc[3]	= m*255.0F + 0.49F;
                            } else {
                              fc[0]	= m*255.0F + 0.49F;
                            }
                          }
                        }
                      }

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
                  
                } /* end of ray while */
              
              if(blend_colors) {
                
                float red_min = _0;
                float green_min = _0;
                float blue_min = _0;
                float red_part;
                float green_part;
                float blue_part;

                if(I->BigEndian) {
                  fc[0] = (float)(0xFF&(*pixel>>24));
                  fc[1] = (float)(0xFF&(*pixel>>16));
                  fc[2] = (float)(0xFF&(*pixel>>8));
                  cc3   =        (0xFF&(*pixel));
                } else {
                  cc3   =        (0xFF&(*pixel>>24));
                  fc[2] = (float)(0xFF&(*pixel>>16));
                  fc[1] = (float)(0xFF&(*pixel>>8));
                  fc[0] = (float)(0xFF&(*pixel));
                }

                red_part = red_blend * fc[0];
                green_part = green_blend * fc[1];
                blue_part = blue_blend * fc[2];
                
                red_min = (green_part>blue_part) ? green_part : blue_part;
                green_min = (red_part>blue_part) ? red_part : blue_part;
                blue_min = (green_part>red_part) ? green_part : red_part;
                
                if(fc[0]<red_min) fc[0] = red_min;
                if(fc[1]<green_min) fc[1] = green_min;
                if(fc[2]<blue_min) fc[2] = blue_min;

                cc0 = (uint)(fc[0]);
                cc1 = (uint)(fc[1]);
                cc2 = (uint)(fc[2]);
                
                if(cc0 > 255) cc0 = 255;
                if(cc1 > 255) cc1 = 255;
                if(cc2 > 255) cc2 = 255;
                
                if(I->BigEndian) 
                  *pixel = (cc0<<24)|(cc1<<16)|(cc2<<8)|cc3;
                else
                  *pixel = (cc3<<24)|(cc2<<16)|(cc1<<8)|cc0;
              }
            
              if(!T->edging) break;
              /* if here, then we're edging...
                 so accumulate averages */
              { 
                unsigned int value;
                value = *pixel;
                edge_avg[0] += value&0xFF;
                edge_avg[1] += (value>>8)&0xFF;
                edge_avg[2] += (value>>16)&0xFF;
                edge_avg[3] += (value>>24)&0xFF;
                
                edge_cnt++;
              }
              
            } /* end of edging while */
            pixel++;
            
         }	/* end of for */
         
		}	/* end of if */
		
	}	/* end of for */
	
	/*  if(T->n_thread>1) 
	  printf(" Ray: Thread %d: Complete.\n",T->phase+1);*/
	MapCacheFree(&SceneCall.cache,T->phase,cCache_map_scene_cache);
	
	if(shadows&&(I->NBasis>1))
		MapCacheFree(&ShadeCall.cache,T->phase,cCache_map_shadow_cache);
	
	return (n_hit);
}

/* this is both an antialias and a slight blur */

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

int RayAntiThread(CRayAntiThreadInfo *T)
{
	int a;
	int		src_row_pixels;
	unsigned int part;
	unsigned int acc;
   unsigned int z[36],zm[36];
	
	unsigned int *pSrc;
	unsigned int *pDst;
   /*   unsigned int m00FF=0x00FF,mFF00=0xFF00,mFFFF=0xFFFF;*/
	int width;
	int height;
	int x,y,yy;
	unsigned int *p;
	int offset = 0;
	
	OrthoBusyFast(9,10);
	width	= (T->width/T->mag) - 2;
	height = (T->height/T->mag) - 2;
	
	src_row_pixels	= T->width;

	offset = (T->phase * height)/T->n_thread;
	offset = offset - (offset % T->n_thread) + T->phase;

	for(yy = 0; yy< height; yy++ )
	{
		y = (yy + offset) % height; /* make sure threads write to different pages */
				
		if((y % T->n_thread) == T->phase)	/* this is my scan line */
        {
          pSrc	= T->image + src_row_pixels * (y*T->mag);
          pDst	= T->image_copy + width * y ;	
          switch(T->mag) {
          case 2:
            for(x = 0; x < width; x++)
              {
                p	= pSrc + (x * T->mag);
                
                z[0 ]	= p[0];
                z[1 ]	= p[1];
                z[2 ]	= p[2];
                z[3 ]	= p[3];
                
                p	+= src_row_pixels;
                
                z[4 ]	= p[0];
                z[5 ]	= p[1];
                z[6 ]	= p[2];
                z[7 ]	= p[3];
                
                p	+= src_row_pixels;
                
                z[8 ]	= p[0];
                z[9 ]	= p[1];
                z[10]	= p[2];
                z[11]	= p[3];
                
                p	+= src_row_pixels;
                
                z[12]	= p[0];
                z[13]	= p[1];
                z[14]	= p[2];
                z[15]	= p[3];
                
                for( a = 0; a < 16; a += 4 ) 
                  {
                    zm[a+0 ] = z[a+0] & mFFFF; /* move half to zm */
                    zm[a+1 ] = z[a+1] & mFFFF;
                    zm[a+2 ] = z[a+2] & mFFFF;
                    zm[a+3 ] = z[a+3] & mFFFF;
                    
                    z[a+0 ] = (z[a+0 ] >> 16) & mFFFF; /* keep rest in z */
                    z[a+1 ] = (z[a+1 ] >> 16) & mFFFF;
                    z[a+2 ] = (z[a+2 ] >> 16) & mFFFF;
                    z[a+3 ] = (z[a+3 ] >> 16) & mFFFF;
                  }
                
                combine4by4(part,z,m00FF);
                acc	= (part<<16);
                combine4by4(part,z,mFF00);
                acc |= (part<<16);
                combine4by4(part,zm,m00FF);
                acc |= part;
                combine4by4(part,zm,mFF00);
                acc |= (acc | part);

                #ifdef _PYMOL_OSX
                acc = optimizer_workaround1u(acc);
                #endif

                *(pDst++) = acc;
              }
            break;
          case 3:
            for(x = 0; x < width; x++)
              {
                p	= pSrc + (x * T->mag);
                
                z[0 ]	= (*(p  ));
                z[1 ]	= (*(p+1));
                z[2 ]	= (*(p+2));
                z[3 ]	= (*(p+3));
                z[4 ]	= (*(p+4));
                
                p	+= src_row_pixels;
                
                z[5 ]	= (*(p  ));
                z[6 ]	= (*(p+1));
                z[7 ]	= (*(p+2));
                z[8 ]	= (*(p+3));
                z[9 ]	= (*(p+4));
                
                p	+= src_row_pixels;
                
                z[10]	= (*(p  ));
                z[11]	= (*(p+1));
                z[12]	= (*(p+2));
                z[13]	= (*(p+3));
                z[14]	= (*(p+4));
                
                p	+= src_row_pixels;
                
                z[15]	= (*(p  ));
                z[16]	= (*(p+1));
                z[17]	= (*(p+2));
                z[18]	= (*(p+3));
                z[19]	= (*(p+4));      
          
                p	+= src_row_pixels;
                
                z[20]	= (*(p  ));
                z[21]	= (*(p+1));
                z[22]	= (*(p+2));
                z[23]	= (*(p+3));
                z[24]	= (*(p+4));                

                for( a = 0; a < 25; a += 5 ) 
                  {
                    zm[a+0 ] = z[a+0] & mFFFF; /* move half to zm */
                    zm[a+1 ] = z[a+1] & mFFFF;
                    zm[a+2 ] = z[a+2] & mFFFF;
                    zm[a+3 ] = z[a+3] & mFFFF;
                    zm[a+4 ] = z[a+4] & mFFFF;
                    
                    z[a   ]	>>= 16; /* keep rest in z */
                    z[a+1 ]	>>= 16;
                    z[a+2 ]	>>= 16;
                    z[a+3 ]	>>= 16;
                    z[a+4 ]	>>= 16;
                  }
                
                combine5by5(part,z,m00FF);
                acc	= (part<<16);
                combine5by5(part,z,mFF00);
                acc	|= (part<<16);
                combine5by5(part,zm,m00FF);
                acc |= part;
                combine5by5(part,zm,mFF00);
                acc |= part;

                #ifdef _PYMOL_OSX
                acc=optimizer_workaround1u(acc);
                #endif
                *(pDst++) = acc;
              }
            break;
          case 4:
            for(x = 0; x < width; x++)
              {
                p	= pSrc + (x * T->mag);
                
                z[0 ]	= (*(p  ));
                z[1 ]	= (*(p+1));
                z[2 ]	= (*(p+2));
                z[3 ]	= (*(p+3));
                z[4 ]	= (*(p+4));
                z[5 ]	= (*(p+5));
                
                p	+= src_row_pixels;
                
                z[6 ]	= (*(p  ));
                z[7 ]	= (*(p+1));
                z[8 ]	= (*(p+2));
                z[9 ]	= (*(p+3));
                z[10]	= (*(p+4));
                z[11]	= (*(p+5));
                
                p	+= src_row_pixels;
                
                z[12]	= (*(p  ));
                z[13]	= (*(p+1));
                z[14]	= (*(p+2));
                z[15]	= (*(p+3));
                z[16]	= (*(p+4));
                z[17]	= (*(p+5));
                
                p	+= src_row_pixels;
                
                z[18]	= (*(p  ));
                z[19]	= (*(p+1));
                z[20]	= (*(p+2));
                z[21]	= (*(p+3));
                z[22]	= (*(p+4));      
                z[23]	= (*(p+5));      
          
                p	+= src_row_pixels;
                
                z[24]	= (*(p  ));
                z[25]	= (*(p+1));
                z[26]	= (*(p+2));
                z[27]	= (*(p+3));
                z[28]	= (*(p+4));                
                z[29]	= (*(p+5));                

                p	+= src_row_pixels;
                
                z[30]	= (*(p  ));
                z[31]	= (*(p+1));
                z[32]	= (*(p+2));
                z[33]	= (*(p+3));
                z[34]	= (*(p+4));                
                z[35]	= (*(p+5));                

                for( a = 0; a < 36; a += 6 ) 
                  {
                    zm[a+0 ] = z[a+0] & mFFFF; /* move half to zm */
                    zm[a+1 ] = z[a+1] & mFFFF;
                    zm[a+2 ] = z[a+2] & mFFFF;
                    zm[a+3 ] = z[a+3] & mFFFF;
                    zm[a+4 ] = z[a+4] & mFFFF;
                    zm[a+5 ] = z[a+5] & mFFFF;
                    
                    z[a   ]	>>= 16; /* keep rest in z */
                    z[a+1 ]	>>= 16;
                    z[a+2 ]	>>= 16;
                    z[a+3 ]	>>= 16;
                    z[a+4 ]	>>= 16;
                    z[a+5 ]	>>= 16;
                  }
                
                combine6by6(part,z,m00FF);
                acc	= (part<<16);
                combine6by6(part,z,mFF00);
                acc	+= (part<<16);
                combine6by6(part,zm,m00FF);
                acc	+= part;
                combine6by6(part,zm,mFF00);
                acc |= part;

                #ifdef _PYMOL_OSX
                acc = optimizer_workaround1u(acc);
                #endif

                *(pDst++) = acc;
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

/*========================================================================*/
void RayRender(CRay *I,int width,int height,unsigned int *image,
               float front,float back,
               double timing,float angle)
{
  int a;
  float *v,light[3];
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
  int mag=1;
  int oversample_cutoff;

#ifdef PROFILE_BASIS
  n_cells = 0;
  n_prims = 0;
  n_triangles = 0;
  n_spheres = 0;
  n_cylinders = 0;
  n_sausages = 0;
  n_skipped = 0;
#endif

  n_thread  = (int)SettingGet(cSetting_max_threads);
  if(n_thread<1)
    n_thread=1;
  if(n_thread>MAX_RAY_THREADS)
    n_thread = MAX_RAY_THREADS;
  opaque_back = (int)SettingGet(cSetting_ray_opaque_background);
  BasisSetFudge(SettingGet(cSetting_ray_triangle_fudge));
  shadows = (int)SettingGet(cSetting_ray_shadows);
  antialias = (int)SettingGet(cSetting_antialias);
  if(antialias<0) antialias=0;
  if(antialias>4) antialias=4;
  mag = antialias;
  if(mag<1) mag=1;

  if(antialias>1) {
    width=(width+2)*mag;
    height=(height+2)*mag;
    image_copy = image;
    buffer_size = mag*mag*width*height;
    image = CacheAlloc(unsigned int,buffer_size,0,cCache_ray_antialias_buffer);
    ErrChkPtr(image);
  } else {
    buffer_size = width*height;
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

  OrthoBusyFast(3,20);

  PRINTFB(FB_Ray,FB_Blather) 
    " RayNew: Background = %x %d %d %d\n",background,(int)(bkrd[0]*255),
    (int)(bkrd[1]*255),(int)(bkrd[2]*255)
    ENDFB;

  if(!I->NPrimitive) { /* nothing to render! */
    fill(image,background,width * (unsigned int)height);
  } else {
    
    RayExpandPrimitives(I);
    RayTransformFirst(I);
    
    now = UtilGetSeconds()-timing;

	 PRINTFB(FB_Ray,FB_Blather)
      " Ray: processed %i graphics primitives in %4.2f sec.\n",I->NPrimitive,now
      ENDFB;

    I->NBasis=3; /* light source */
    BasisInit(I->Basis+2,2);
    
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
      RayTransformBasis(I,I->Basis+2,2);
    }

    if(shadows&&(n_thread>1)) { /* parallel execution */

      CRayHashThreadInfo thread_info[2];
      
      thread_info[0].basis = I->Basis+1;
      thread_info[0].vert2prim = I->Vert2Prim;
      thread_info[0].prim = I->Primitive;
      thread_info[0].clipBox = I->Volume;
      thread_info[0].image = image;
      thread_info[0].background = background;
      thread_info[0].bytes = width * (unsigned int)height;
      thread_info[0].phase = 0;
      thread_info[0].ray = I; /* for compute box */

      thread_info[1].basis = I->Basis+2;
      thread_info[1].vert2prim = I->Vert2Prim;
      thread_info[1].prim = I->Primitive;
      thread_info[1].clipBox = NULL;
      thread_info[1].phase = 1;
      RayHashSpawn(thread_info,2);
      
    } else { /* serial execution */
      BasisMakeMap(I->Basis+1,I->Vert2Prim,I->Primitive,I->Volume,0,cCache_ray_map);
      if(shadows) {
        BasisMakeMap(I->Basis+2,I->Vert2Prim,I->Primitive,NULL,1,cCache_ray_map);
      }

      /* serial tasks which RayHashThread does in parallel mode */

      fill(image,background,width * (unsigned int)height);
      RayComputeBox(I);
       
    }

    OrthoBusyFast(5,20);
    now = UtilGetSeconds()-timing;

#ifdef _MemoryDebug_ON
    if(shadows) {
      PRINTFB(FB_Ray,FB_Blather)
        " Ray: voxels: [%4.2f:%dx%dx%d], [%4.2f:%dx%dx%d], %d MB, %4.2f sec.\n",
        I->Basis[1].Map->Div,   I->Basis[1].Map->Dim[0],
        I->Basis[1].Map->Dim[1],I->Basis[1].Map->Dim[2],
        I->Basis[2].Map->Div,   I->Basis[2].Map->Dim[0],
        I->Basis[2].Map->Dim[2],I->Basis[2].Map->Dim[2],
        (int)(MemoryDebugUsage()/(1024.0*1024)),
        now
        ENDFB;
    } else {
      PRINTFB(FB_Ray,FB_Blather)
        " Ray: voxels: [%4.2f:%dx%dx%d], %d MB, %4.2f sec.\n",
        I->Basis[1].Map->Div,   I->Basis[1].Map->Dim[0],
        I->Basis[1].Map->Dim[1],I->Basis[1].Map->Dim[2],
        (int)(MemoryDebugUsage()/(1024.0*1024)),
        now
        ENDFB;
    }
#else
    if(shadows) {
      PRINTFB(FB_Ray,FB_Blather)
        " Ray: voxels: [%4.2f:%dx%dx%d], [%4.2f:%dx%dx%d], %4.2f sec.\n",
        I->Basis[1].Map->Div,   I->Basis[1].Map->Dim[0],
        I->Basis[1].Map->Dim[1],I->Basis[1].Map->Dim[2],
        I->Basis[2].Map->Div,   I->Basis[2].Map->Dim[0],
        I->Basis[2].Map->Dim[2],I->Basis[2].Map->Dim[2],
        now
        ENDFB;
    } else {
      PRINTFB(FB_Ray,FB_Blather)
        " Ray: voxels: [%4.2f:%dx%dx%d], %4.2f sec.\n",
        I->Basis[1].Map->Div,   I->Basis[1].Map->Dim[0],
        I->Basis[1].Map->Dim[1],I->Basis[1].Map->Dim[2],
        now
        ENDFB;
    }

#endif
    /* IMAGING */
        
    {
		/* now spawn threads as needed */
		CRayThreadInfo rt[MAX_RAY_THREADS];
      int x_start,y_start;
      int x_stop,y_stop;

      x_start = (int)((width * (I->min_box[0] - I->Volume[0]))/I->Range[0]) - 2;
      x_stop  = (int)((width * (I->max_box[0] - I->Volume[0]))/I->Range[0]) + 2;
      
      y_stop = (int)((height * (I->max_box[1] - I->Volume[2]))/I->Range[1]) + 2;
      y_start  = (int)((height * (I->min_box[1] - I->Volume[2]))/I->Range[1]) - 2;
                      
      if(x_start<0) x_start = 0;
      if(y_start<0) y_start = 0;
      if(x_stop>width) x_stop = width;
      if(y_stop>height) y_stop = height;

      oversample_cutoff = (int)SettingGet(cSetting_ray_oversample_cutoff);

      if(!antialias)
        oversample_cutoff = 0;

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
			rt[a].border = mag-1;
			rt[a].front = front;
			rt[a].back = back;
			rt[a].fore_mask = fore_mask;
			rt[a].bkrd = bkrd;
			rt[a].background = background;
			rt[a].phase = a;
			rt[a].n_thread = n_thread;
         rt[a].edging = NULL;
         rt[a].edging_cutoff = oversample_cutoff; /* info needed for busy indicator */

			copy3f(spec_vector,rt[a].spec_vector);
			}
		
		if(n_thread<=1)
        RayTraceThread(rt);
		else 
        RayTraceSpawn(rt,n_thread);

      if(oversample_cutoff) { /* perform edge oversampling, if requested */
        unsigned int *edging;

        edging = CacheAlloc(unsigned int,buffer_size,0,cCache_ray_edging_buffer);
        
        memcpy(edging,image,buffer_size*sizeof(unsigned int));

        for(a=0;a<n_thread;a++) {
          rt[a].edging = edging;
        }

        if(n_thread<=1)
          RayTraceThread(rt);
        else 
          RayTraceSpawn(rt,n_thread);

        CacheFreeP(edging,0,cCache_ray_edging_buffer,false);
      }
    }
  }
  
  if(antialias>1) {
    {
		/* now spawn threads as needed */
		CRayAntiThreadInfo rt[MAX_RAY_THREADS];

		for(a=0;a<n_thread;a++) 
		{
			rt[a].width = width;
			rt[a].height = height;
         rt[a].image = image;
         rt[a].image_copy = image_copy;
         rt[a].phase = a;
         rt[a].mag = mag; /* fold magnification */
         rt[a].n_thread = n_thread;
      }
		
		if(n_thread<=1)
        RayAntiThread(rt);
		else 
        RayAntiSpawn(rt,n_thread);
    }
    CacheFreeP(image,0,cCache_ray_antialias_buffer,false);
    image = image_copy;
  }

  PRINTFD(FB_Ray)
    " RayRender: n_hit %d\n",n_hit
    ENDFD;
#ifdef PROFILE_BASIS

  printf("int n_cells = %d;\nint n_prims = %d;\nint n_triangles = %8.3f;\nint n_spheres = %8.3f;\nint n_cylinders = %8.3f;\nint n_sausages = %8.3f;\nint n_skipped = %8.3f;\n",
         n_cells,
         n_prims,
         n_triangles/((float)n_cells),
         n_spheres/((float)n_cells),
         n_cylinders/((float)n_cells),
         n_sausages/((float)n_cells),
         n_skipped/((float)n_cells));
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

  VLACacheCheck(I->Primitive,CPrimitive,I->NPrimitive,0,cCache_ray_primitive);
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

  VLACacheCheck(I->Primitive,CPrimitive,I->NPrimitive,0,cCache_ray_primitive);
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

  VLACacheCheck(I->Primitive,CPrimitive,I->NPrimitive,0,cCache_ray_primitive);
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

  VLACacheCheck(I->Primitive,CPrimitive,I->NPrimitive,0,cCache_ray_primitive);
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
  VLACacheCheck(I->Primitive,CPrimitive,I->NPrimitive,0,cCache_ray_primitive);
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

  I->Basis=CacheAlloc(CBasis,3,0,cCache_ray_basis);
  BasisInit(I->Basis,0);
  BasisInit(I->Basis+1,1);
  I->Vert2Prim=VLACacheAlloc(int,1,0,cCache_ray_vert2prim);
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
	 I->Primitive=VLACacheAlloc(CPrimitive,10000,0,cCache_ray_primitive);  
  if(!I->Vert2Prim) 
	 I->Vert2Prim=VLACacheAlloc(int,10000,0,cCache_ray_vert2prim);
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
	 BasisFinish(&I->Basis[a],a);
  }
  I->NBasis=0;
  VLACacheFreeP(I->Primitive,0,cCache_ray_primitive,false);
  VLACacheFreeP(I->Vert2Prim,0,cCache_ray_vert2prim,false);
}
/*========================================================================*/
void RayFree(CRay *I)
{
  RayRelease(I);
  CacheFreeP(I->Basis,0,cCache_ray_basis,false);
  VLACacheFreeP(I->Vert2Prim,0,cCache_ray_vert2prim,false);
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



