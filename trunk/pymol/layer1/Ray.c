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
#include<math.h>

#include"Base.h"
#include"MemoryDebug.h"
#include"Err.h"
#include"Vector.h"
#include"OOMac.h"
#include"Setting.h"
#include"Ortho.h"
#include"Util.h"
#include"Ray.h"

#ifndef RAY_SMALL
#define RAY_SMALL 0.00001
#endif

/* BASES 
   0 contains untransformed vertices (vector size = 3)
	1 contains transformed vertices (vector size = 3)
	2 contains transformed vertices for shadowing 
*/

typedef float float3[3];
typedef float float4[4];

void RayRelease(CRay *I);

void RaySetup(CRay *I);
void RayColor3fv(CRay *I,float *v);
void RaySphere3fv(CRay *I,float *v,float r);
void RayCylinder3fv(CRay *I,float *v1,float *v2,float r,float *c1,float *c2);

void RayTriangle3fv(CRay *I,
						  float *v1,float *v2,float *v3,
						  float *n1,float *n2,float *n3,
						  float *c1,float *c2,float *c3);

void RayApplyMatrix33( unsigned int n, float3 *q, const float m[16],
							float3 *p );

void RayExpandPrimitives(CRay *I);
void RayTransformFirst(CRay *I);
void RayTransformBasis(CRay *I,CBasis *B);

int PrimitiveSphereHit(CRay *I,float *v,float *n,float *minDist,int except);

void RayReflectSphere(CRay *I,RayInfo *r);

void RayTransformNormals33( unsigned int n, float3 *q, const float m[16],float3 *p );

/*========================================================================*/
void RayReflectSphere(CRay *I,RayInfo *r)
{
  
  r->impact[0]=r->base[0]; 
  r->impact[1]=r->base[1]; 
  r->impact[2]=r->base[2]-r->dist;
  
  r->surfnormal[0]=r->impact[0]-r->sphere[0];
  r->surfnormal[1]=r->impact[1]-r->sphere[1];
  r->surfnormal[2]=r->impact[2]-r->sphere[2];
  
  normalize3f(r->surfnormal);
  
  r->dotgle = -r->surfnormal[2]; 
  
  r->reflect[0]= - ( 2 * r->dotgle * r->surfnormal[0] );
  r->reflect[1]= - ( 2 * r->dotgle * r->surfnormal[1] );
  r->reflect[2]= -1.0 - ( 2 * r->dotgle * r->surfnormal[2] );
  
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
  basis->MaxRadius = 0.0;
  basis->MinVoxel = 0.0;
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
		if(basis->MinVoxel<0.001)
		  basis->MinVoxel=0.001;
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
		I->Vert2Prim[nVert]=a;
		basis->Radius[nVert]=I->Primitive[a].r1;
		basis->Radius2[nVert]=I->Primitive[a].r1*I->Primitive[a].r1; /*precompute*/
		if(basis->Radius[nVert]>basis->MinVoxel)
		  basis->MinVoxel=basis->Radius[nVert];
		subtract3f(I->Primitive[a].v2,I->Primitive[a].v1,n0);
		I->Primitive[a].l1=length3f(n0);
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
void RayTransformFirst(CRay *I)
{
  CBasis *basis0,*basis1;
  CPrimitive *prm;
  int a;

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
void RayRender(CRay *I,int width,int height,unsigned int *image,float front,float back)
{
  int x,y;
  int a,b;
  unsigned int *p;
  float excess=0.0;
  float dotgle;
  float bright,direct_cmp,reflect_cmp,*v,fc[3];
  float ambient,direct,lreflect,ft;
  unsigned int c[3],aa;
  unsigned int *image_copy = NULL;
  int i;
  unsigned int background,buffer_size,z[12],tot;
  int antialias;
  RayInfo r1,r2;
  double timing;
  
  /* SETUP */
  
  timing = UtilGetSeconds();

  antialias = (int)SettingGet(cSetting_antialias);
  if(antialias>0) {
	 width=width*2;
	 height=height*2;
	 image_copy = image;
	 buffer_size = 4*width*height;
	 image = (void*)Alloc(char,buffer_size);
	 ErrChkPtr(image);
  }

  v=SettingGetfv(cSetting_bg_rgb);
  if(I->BigEndian) {
	 background = 0x000000FF|
		((0xFF& ((unsigned int)(v[0]*255))) <<24)|
		((0xFF& ((unsigned int)(v[1]*255))) <<16)|
		((0xFF& ((unsigned int)(v[2]*255)) <<8 ));
  } else {
	 background = 0xFF000000|
		((0xFF& ((unsigned int)(v[2]*255))) <<16)|
		((0xFF& ((unsigned int)(v[1]*255))) <<8)|
		((0xFF& ((unsigned int)(v[0]*255)) ));
  }

  if(!I->NPrimitive) {
    p=(unsigned int*)image; 
    for(x=0;x<width;x++)
      for(y=0;y<height;y++)
        *(p++)=background;
  } else {
    
    ambient=SettingGet(cSetting_ambient);
    lreflect=SettingGet(cSetting_reflect);
    direct=SettingGet(cSetting_direct);
    
    RayExpandPrimitives(I);
    RayTransformFirst(I);
    
	 printf(" Ray: %i primitives\n",I->NPrimitive);
    BasisMakeMap(I->Basis+1,I->Vert2Prim,I->Primitive,I->Volume);

    I->NBasis=3; /* light source */
    BasisInit(I->Basis+2);
    
    v=SettingGetfv(cSetting_light);
    
    I->Basis[2].LightNormal[0]=v[0];
    I->Basis[2].LightNormal[1]=v[1];
    I->Basis[2].LightNormal[2]=v[2];
    normalize3f(I->Basis[2].LightNormal);
    
    BasisSetupMatrix(I->Basis+2);
    RayTransformBasis(I,I->Basis+2);
    BasisMakeMap(I->Basis+2,I->Vert2Prim,I->Primitive,NULL);

    printf(" Ray: hash spacing: %4.2f/%4.2f\n",
			  I->Basis[1].Map->Div,I->Basis[2].Map->Div);
    
    /* IMAGING */
    
    /* erase buffer */
    
    p=(unsigned int*)image; 
    if(I->BigEndian)
      for(a=0;a<height;a++)
        for(b=0;b<width;b++)
          *p++=0x000000FF;  
    else 
      for(a=0;a<height;a++)
        for(b=0;b<width;b++)
          *p++=0xFF000000;  
    
    /* ray-trace */
	 r1.base[2]=0.0;
    for(x=0;x<width;x++)
      {
        if(!(x&0xF)) OrthoBusyFast(x,width); /* don't slow down rendering too much */
        r1.base[0]=(((float)x)/width)*I->Range[0]+I->Volume[0];
        for(y=0;y<height;y++)
          {
            r1.base[1]=(((float)y)/height)*I->Range[1]+I->Volume[2];
            
            i=BasisHit(I->Basis+1,&r1,-1,I->Vert2Prim,I->Primitive,false,front,back);
            
            if(i>=0) {

              if(r1.prim->type==cPrimTriangle) {
					 BasisReflectTriangle(I->Basis+1,&r1,i,fc);
				  } else {
					 RayReflectSphere(I,&r1);/*,zRay,sphere,dist,surf,impact,reflect,&dotgle);*/
					 if(r1.prim->type==cPrimCylinder) {
						ft = r1.tri1;
						fc[0]=(r1.prim->c1[0]*(1-ft))+(r1.prim->c2[0]*ft);
						fc[1]=(r1.prim->c1[1]*(1-ft))+(r1.prim->c2[1]*ft);
						fc[2]=(r1.prim->c1[2]*(1-ft))+(r1.prim->c2[2]*ft);
					 } else {
						fc[0]=r1.prim->c1[0];
						fc[1]=r1.prim->c1[1];
						fc[2]=r1.prim->c1[2];
					 }
				  }

              dotgle=-r1.dotgle;
              direct_cmp=(dotgle+(pow(dotgle,SettingGet(cSetting_power))))/2.0;
              
              matrix_transform33f3f(I->Basis[2].Matrix,r1.impact,r2.base);
              
              if(BasisHit(I->Basis+2,&r2,i,I->Vert2Prim,I->Primitive,true,0.0,0.0)<0) {
					 
					 dotgle=-dot_product3f(r1.surfnormal,I->Basis[2].LightNormal);
                if(dotgle<0.0) dotgle=0.0;
                reflect_cmp=(dotgle+(pow(dotgle,SettingGet(cSetting_power))))/2.0;
                excess=pow(dotgle,SettingGet(cSetting_spec_power))*
                  SettingGet(cSetting_spec_reflect);
              } else {
                excess=0.0;
                reflect_cmp=0.0;
              }
				  
              
              bright=ambient+(1.0-ambient)*(direct*direct_cmp+
                                            (1.0-direct)*direct_cmp*lreflect*reflect_cmp);
              
              if(bright>1.0) bright=1.0;
              if(bright<0.0) bright=0.0;
				  
				  
				  c[0]=(bright*fc[0]+excess)*255.0;
				  c[1]=(bright*fc[1]+excess)*255.0;
				  c[2]=(bright*fc[2]+excess)*255.0;

              if(c[0]>255.0) c[0]=255.0;
              if(c[1]>255.0) c[1]=255.0;
              if(c[2]>255.0) c[2]=255.0;
              /*			{ int a1,b1,c1;
                        if(MapInsideXY(I->Basis[1].Map,base,&a1,&b1,&c1))				
                        {
                        c[0]=255.0*((float)a1)/I->Basis[1].Map->iMax[0];
                        c[1]=255.0*((float)b1)/I->Basis[1].Map->iMax[1];
                        c[2]=0xFF&(c[0]*c[1]);
                        }}
              */
              
              if(I->BigEndian) {
                *(image+((width)*y)+x)=
                  0x000000FF|(c[0]<<24)|(c[1]<<16)|(c[2]<<8);
              } else {
                *(image+((width)*y)+x)=
                  0xFF000000|(c[2]<<16)|(c[1]<<8)|c[0];
              }
            } else {
              
              *(image+((width)*y)+x)=
                background;
            }
          }
      }
  }

  if(antialias) {
	 OrthoBusyFast(9,10);
	 width=width/2;
	 height=height/2;
	 for(x=1;x<(width-1);x++)
		for(y=1;y<(height-1);y++)
		  {
			 aa=0;

			 p = image+((width*2)*(y*2-1))+(x*2);
			 
			 z[4] = (*(p));
			 z[5] = (*(p+1));
			 p+=(width*2);
			 z[6] = (*(p-1));
			 z[0] = (*(p));
			 z[1] = (*(p+1));
			 z[7] = (*(p+2));
			 p+=(width*2);
			 z[8] = (*(p-1));
			 z[2] = (*(p));
			 z[3] = (*(p+1));
			 z[9] = (*(p+2));
			 p+=(width*2);
			 z[10] = (*(p));
			 z[11] = (*(p+1));

			 if(I->BigEndian) { 
				for(a=0;a<12;a++)
				  z[a]=z[a]>>8;
			 }

			 tot=0;
			 for(a=0;a<12;a++)
				tot+=(z[a]&0xFF);
			 for(a=0;a<4;a++)
				tot+=(z[a]&0xFF);
			 aa=aa|(0xFF&(tot>>4));
			 
			 tot=0;
			 for(a=0;a<12;a++)
				tot+=(z[a]&0xFF00);
			 for(a=0;a<4;a++)
				tot+=(z[a]&0xFF00);
			 aa=aa|(0xFF00&(tot>>4));

			 tot=0;
			 for(a=0;a<12;a++)
				tot+=(z[a]&0xFF0000);
			 for(a=0;a<4;a++)
				tot+=(z[a]&0xFF0000);
			 aa=aa|(0xFF0000&(tot>>4));			 
			 
			 if(I->BigEndian) {
				aa=aa<<8;
				aa=aa|0x000000FF;
			 } else {
				aa=aa|0xFF000000;
			 }
			 
			 *(image_copy+((width)*y)+x) = aa;			 
		  }

	 for(x=0;x<width;x++)
		for(y=0;y<height;y=y+(height-1))
		  {
			 aa=0;

			 p = image+((width*2)*(y*2))+(x*2);
			 
			 z[0] = (*(p));
			 z[1] = (*(p+1));
			 p+=(width*2);
			 z[2] = (*(p));
			 z[3] = (*(p+1));

			 if(I->BigEndian) { 
				for(a=0;a<12;a++)
				  z[a]=z[a]>>8;
			 }

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
			 
			 if(I->BigEndian) {
				aa=aa<<8;
				aa=aa|0x000000FF;
			 } else {
				aa=aa|0xFF000000;
			 }
			 *(image_copy+((width)*y)+x) = aa;			 
		  }

	 for(x=0;x<width;x=x+(width-1))
		for(y=0;y<height;y++)
		  {
			 aa=0;

			 p = image+((width*2)*(y*2))+(x*2);
			 
			 z[0] = (*(p));
			 z[1] = (*(p+1));
			 p+=(width*2);
			 z[2] = (*(p));
			 z[3] = (*(p+1));

			 if(I->BigEndian) { 
				for(a=0;a<12;a++)
				  z[a]=z[a]>>8;
			 }

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
			 
			 if(I->BigEndian) {
				aa=aa<<8;
				aa=aa|0x000000FF;
			 } else {
				aa=aa|0xFF000000;
			 }
			 *(image_copy+((width)*y)+x) = aa;			 
		  }

	 FreeP(image);
	 image=image_copy;
  }

  timing = UtilGetSeconds()-timing;
  printf(" Ray: rendering time %4.2f sec. (%3.1f fph).\n",
			timing,3600/timing);

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

  vv=p->v1;
  (*vv++)=(*v++);
  (*vv++)=(*v++);
  (*vv++)=(*v++);
  vv=p->c1;
  v=I->CurColor;
  (*vv++)=(*v++);
  (*vv++)=(*v++);
  (*vv++)=(*v++);
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

  vv=p->v1;
  (*vv++)=(*v1++);
  (*vv++)=(*v1++);
  (*vv++)=(*v1++);
  vv=p->v2;
  (*vv++)=(*v2++);
  (*vv++)=(*v2++);
  (*vv++)=(*v2++);

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

  VLACheck(I->Primitive,CPrimitive,I->NPrimitive);
  p = I->Primitive+I->NPrimitive;

  p->type = cPrimTriangle;

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
  l1=length3f(s1);
  l2=length3f(s2);
  l3=length3f(s3);
  if(l2>l1) { if(l3>l2)	l1=l3; else	l1=l2;  }
  /* store cutoff distance */

  p->r1=l1*0.6;

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

  I->NPrimitive++;

}
/*========================================================================*/
CRay *RayNew(void)
{
  unsigned int test;
  unsigned char *testPtr;

  OOAlloc(CRay);
  
  test = 0xFF000000;
  testPtr = (unsigned char*)&test;
  I->BigEndian = *testPtr;

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
  I->fTriangle3fv=RayTriangle3fv;
  
  return(I);
}
/*========================================================================*/
void RayPrepare(CRay *I,float v0,float v1,float v2,float v3,float v4,float v5,float *mat)
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

  for(a=0;a<16;a++)
    I->ModelView[a]=mat[a];
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
      for (i=0;i<n;i++) {
         float p0 = p[i][0], p1 = p[i][1], p2 = p[i][2];
         q[i][0] = m0 * p0 + m4  * p1 + m8 * p2 + m12;
         q[i][1] = m1 * p0 + m5  * p1 + m9 * p2 + m13;
      }
   }
   {
      unsigned int i;
      float m2 = m[2],  m6 = m[6],  m10 = m[10],  m14 = m[14];
      float m3 = m[3],  m7 = m[7],  m11 = m[11],  m15 = m[15];
      if (m3==0.0F && m7==0.0F && m11==0.0F && m15==1.0F) {
         /* common case */
         for (i=0;i<n;i++) {
            float p0 = p[i][0], p1 = p[i][1], p2 = p[i][2];
            q[i][2] = m2 * p0 + m6 * p1 + m10 * p2 + m14;
				/*				q[i][3] = 1.0F;*/
         }
      }
      else {
         /* general case */
         for (i=0;i<n;i++) {
            float p0 = p[i][0], p1 = p[i][1], p2 = p[i][2];
            q[i][2] = m2 * p0 + m6 * p1 + m10 * p2 + m14;
				/*				q[i][3] = m3 * p0 + m7 * p1 + m11 * p2 + m15; */
         }
      }
   }
}

void RayTransformNormals33( unsigned int n, float3 *q, const float m[16],float3 *p )
{
   {
      unsigned int i;
      float m0 = m[0],  m4 = m[4],  m8 = m[8];
      float m1 = m[1],  m5 = m[5],  m9 = m[9];
      for (i=0;i<n;i++) {
         float p0 = p[i][0], p1 = p[i][1], p2 = p[i][2];
         q[i][0] = m0 * p0 + m4  * p1 + m8 * p2;
         q[i][1] = m1 * p0 + m5  * p1 + m9 * p2;
      }
   }
   {
      unsigned int i;
      float m2 = m[2],  m6 = m[6],  m10 = m[10];
      float m3 = m[3],  m7 = m[7],  m11 = m[11],  m15 = m[15];
      if (m3==0.0F && m7==0.0F && m11==0.0F && m15==1.0F) {
         /* common case */
         for (i=0;i<n;i++) {
            float p0 = p[i][0], p1 = p[i][1], p2 = p[i][2];
            q[i][2] = m2 * p0 + m6 * p1 + m10 * p2;
         }
      }
      else {
         /* general case */
         for (i=0;i<n;i++) {
            float p0 = p[i][0], p1 = p[i][1], p2 = p[i][2];
            q[i][2] = m2 * p0 + m6 * p1 + m10 * p2;
         }
      }
   }
	{      
	  unsigned int i;
	  for (i=0;i<n;i++) { /* renormalize - can we do this to the matrix instead? */
		 normalize3f(q[i]);
	  }
	}
}



