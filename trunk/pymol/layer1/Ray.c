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
#include<GL/gl.h>
#include<math.h>
#include<values.h>

#include"Base.h"
#include"MemoryDebug.h"
#include"Err.h"
#include"Vector.h"
#include"OOMac.h"
#include"Setting.h"
#include"Ortho.h"

#include"Ray.h"

#ifndef RAY_SMALL
#define RAY_SMALL 0.00001
#endif

/* BASES 

   0 contains untransformed vertices (vector size = 3)
	1 contains transformed vertices (vector size = 3)

 */

typedef GLfloat GLfloat3[3];
typedef GLfloat GLfloat4[4];

void RayRelease(CRay *I);

void RaySetup(CRay *I);
void RayColor3fv(CRay *I,float *v);
void RaySphere3fv(CRay *I,float *v,float r);
void RayCylinder3fv(CRay *I,float *v1,float *v2,float r);
void RayApplyMatrix33( GLuint n, GLfloat3 *q, const GLfloat m[16],
							GLfloat3 *p );

void RayExpandPrimitives(CRay *I);
void RayTransformFirst(CRay *I);
void RayTransformBasis(CRay *I,CBasis *B);

int PrimitiveSphereHit(CRay *I,float *v,float *n,float *minDist,int except);

void RayReflect(CRay *I,float *v,float *n,float *p,float d,
					 float *sn,float *ip,float *nn,float *dotgle);
void RayTransformNormals33( GLuint n, GLfloat3 *q, const GLfloat m[16],GLfloat3 *p );

/*========================================================================*/
void RayReflect(CRay *I,float *v,float *n,float *p,float d,
									 float *sn,float *ip,float *nn,float *dotgle)
{
  
  ip[0]=v[0]+d*n[0];
  ip[1]=v[1]+d*n[1];
  ip[2]=v[2]+d*n[2];
  
  sn[0]=ip[0]-p[0];
  sn[1]=ip[1]-p[1];
  sn[2]=ip[2]-p[2];
  
  normalize3f(sn);
 
  (*dotgle) = n[0]*sn[0]+n[1]*sn[1]+n[2]*sn[2];

  nn[0]=n[0]-2*(*dotgle)*sn[0];
  nn[1]=n[1]-2*(*dotgle)*sn[1];
  nn[2]=n[2]-2*(*dotgle)*sn[2];

}
/*========================================================================*/
void RayExpandPrimitives(CRay *I)
{
  int a;
  float *v0,*v1,*n0;
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
	 }
  }

  basis = I->Basis;
  
  VLACheck(basis->Vertex,float,3*nVert);
  VLACheck(basis->Radius,float,nVert);
  VLACheck(basis->Radius2,float,nVert);
  VLACheck(basis->Vert2Normal,int,nVert);
  VLACheck(basis->Normal,float,3*nNorm);
  VLACheck(I->Vert2Prim,int,nVert);
  basis->MaxRadius=0.0;
  basis->NVertex=nVert;
  basis->NNormal=nNorm;

  nVert=0;
  nNorm=0;
  v0=basis->Vertex;
  n0=basis->Normal;
  for(a=0;a<I->NPrimitive;a++) {
	 switch(I->Primitive[a].type) {
	 case cPrimSphere:
		I->Vert2Prim[nVert]=a;
		v1=I->Primitive[a].v1;
		basis->Radius[nVert]=I->Primitive[nVert].r1;
		basis->Radius2[nVert]=I->Primitive[nVert].r1*I->Primitive[nVert].r1; /*precompute*/
		if(basis->Radius[nVert]>basis->MaxRadius)
		  basis->MaxRadius=basis->Radius[nVert];
		(*v0++)=(*v1++);
		(*v0++)=(*v1++);
		(*v0++)=(*v1++);
		nVert++;
		break;
	 case cPrimCylinder:
		I->Vert2Prim[nVert]=a;
		basis->Radius[nVert]=I->Primitive[nVert].r1;
		basis->Radius2[nVert]=I->Primitive[nVert].r1*I->Primitive[nVert].r1; /*precompute*/
		if(basis->Radius[nVert]>basis->MaxRadius)
		  basis->MaxRadius=basis->Radius[nVert];
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
}
/*========================================================================*/
void RayTransformFirst(CRay *I)
{
  CBasis *basis0,*basis1;
  int a;

  basis0 = I->Basis;
  basis1 = I->Basis+1;
  
  VLACheck(basis1->Vertex,float,3*basis0->NVertex);
  VLACheck(basis1->Normal,float,3*basis0->NNormal);
  VLACheck(basis1->Vert2Normal,int,basis0->NVertex);
  VLACheck(basis1->Radius,float,basis0->NVertex);
  VLACheck(basis1->Radius2,float,basis0->NVertex);
  
  RayApplyMatrix33(basis0->NVertex,(GLfloat3*)basis1->Vertex,
					  I->ModelView,(GLfloat3*)basis0->Vertex);

  for(a=0;a<basis0->NVertex;a++)
	 {
		basis1->Radius[a]=basis0->Radius[a];
		basis1->Radius2[a]=basis0->Radius2[a];
		basis1->Vert2Normal[a]=basis0->Vert2Normal[a];
	 }

  basis1->MaxRadius=basis0->MaxRadius;
  basis1->NVertex=basis0->NVertex;

  RayTransformNormals33(basis0->NNormal,(GLfloat3*)basis1->Normal,
					  I->ModelView,(GLfloat3*)basis0->Normal);
  
  basis1->NNormal=basis0->NNormal;
}
/*========================================================================*/
void RayTransformBasis(CRay *I,CBasis *basis1)
{
  CBasis *basis0;
  int a;
  float *v0,*v1,*n0,*n1;

  basis0 = I->Basis+1;

  VLACheck(basis1->Vertex,float,3*basis0->NVertex);
  VLACheck(basis1->Normal,float,3*basis0->NNormal);
  VLACheck(basis1->Vert2Normal,int,basis0->NVertex);
  VLACheck(basis1->Radius,float,basis0->NVertex);
  VLACheck(basis1->Radius2,float,basis0->NVertex);

  v0=basis0->Vertex;
  v1=basis1->Vertex;
  n0=basis0->Normal;
  n1=basis1->Normal;
  for(a=0;a<basis0->NVertex;a++)
	 {
		transform33f3f(basis1->Matrix,v0,v1);
		v0+=3;
		v1+=3;
		basis1->Radius[a]=basis0->Radius[a];
		basis1->Radius2[a]=basis0->Radius2[a];
		basis1->Vert2Normal[a]=basis0->Vert2Normal[a];
	 }
  for(a=0;a<basis0->NNormal;a++)
	 {
		transform33f3f(basis1->Matrix,n0,n1);
		n0+=3;
		n1+=3;
	 }
  basis1->MaxRadius=basis0->MaxRadius;
  basis1->NVertex=basis0->NVertex;
  basis1->NNormal=basis0->NNormal;
}
/*========================================================================*/
void RayRender(CRay *I,int width,int height,unsigned int *image,float front,float back)
{
  int x,y;
  int a,b;
  unsigned int *p;
  float excess=0.0;
  float base[3],impact[3],reflect[3],surf[3],dist,dotgle,dist2;
  float impact2[3],sphere[3];
  float bright,direct_cmp,reflect_cmp,*v;
  float ambient,direct,lreflect;
  unsigned int c[3],aa;
  unsigned int *image_copy = NULL;
  int i;
  unsigned int background,buffer_size,z[12],tot;
  int antialias;

  float zRay[3] = { 0.0,0.0,-1.0};
  /* SETUP */

  antialias = (int)SettingGet(cSetting_antialias);
  if(antialias>0) {
	 width=width*2;
	 height=height*2;
	 image_copy = image;
	 buffer_size = 4*width*height;
	 image = (GLvoid*)Alloc(char,buffer_size);
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
  ambient=SettingGet(cSetting_ambient);
  lreflect=SettingGet(cSetting_reflect);
  direct=SettingGet(cSetting_direct);

  RayExpandPrimitives(I);
  RayTransformFirst(I);

  BasisMakeMap(I->Basis+1,I->Basis[1].MaxRadius,I->Vert2Prim,I->Primitive,I->Volume);

  I->NBasis=3; /* light source */
  BasisInit(I->Basis+2);

  v=SettingGetfv(cSetting_light);

  I->Basis[2].LightNormal[0]=v[0];
  I->Basis[2].LightNormal[1]=v[1];
  I->Basis[2].LightNormal[2]=v[2];
  normalize3f(I->Basis[2].LightNormal);
  
  BasisSetupMatrix(I->Basis+2);
  RayTransformBasis(I,I->Basis+2);
  BasisMakeMap(I->Basis+2,I->Basis[2].MaxRadius,I->Vert2Prim,I->Primitive,NULL);

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
  base[2]=0.0;
  for(x=0;x<width;x++)
	 {
		if(!(x&0xF)) OrthoBusyFast(x,width); /* don't slow down rendering too much */
		base[0]=(((float)x)/width)*I->Range[0]+I->Volume[0];
		for(y=0;y<height;y++)
		  {
			 base[1]=(((float)y)/height)*I->Range[1]+I->Volume[2];
			 
			 i=BasisHit(I->Basis+1,base,&dist,-1,I->Vert2Prim,I->Primitive,sphere,false,
							front,back); 
			 
			 if(i>=0) {

				RayReflect(I,base,zRay,sphere,dist,surf,impact,reflect,&dotgle);
				dotgle=-dotgle;
				direct_cmp=(dotgle+(pow(dotgle,SettingGet(cSetting_power))))/2.0;
				
				transform33f3f(I->Basis[2].Matrix,impact,impact2);

				if(BasisHit(I->Basis+2,impact2,&dist2,i,
								I->Vert2Prim,I->Primitive,sphere,true,0,0)<0) {
				  
				  dotgle=-dot_product3f(surf,I->Basis[2].LightNormal);
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
				c[0]=(bright*I->Primitive[i].c1[0]+excess)*255.0;
				c[1]=(bright*I->Primitive[i].c1[1]+excess)*255.0;
				c[2]=(bright*I->Primitive[i].c1[2]+excess)*255.0;
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
void RayCylinder3fv(CRay *I,float *v1,float *v2,float r)
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
  v1=I->CurColor;
  (*vv++)=(*v1++);
  (*vv++)=(*v1++);
  (*vv++)=(*v1++);
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
  
  return(I);
}
/*========================================================================*/
void RayPrepare(CRay *I,float v0,float v1,float v2,float v3,float v4,float v5)
	  /*prepare for vertex calls */
{

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

  glGetFloatv(GL_MODELVIEW_MATRIX,I->ModelView);
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

void RayApplyMatrix33( GLuint n, GLfloat3 *q, const GLfloat m[16],
                          GLfloat3 *p )
{
   {
      GLuint i;
      GLfloat m0 = m[0],  m4 = m[4],  m8 = m[8],  m12 = m[12];
      GLfloat m1 = m[1],  m5 = m[5],  m9 = m[9],  m13 = m[13];
      for (i=0;i<n;i++) {
         GLfloat p0 = p[i][0], p1 = p[i][1], p2 = p[i][2];
         q[i][0] = m0 * p0 + m4  * p1 + m8 * p2 + m12;
         q[i][1] = m1 * p0 + m5  * p1 + m9 * p2 + m13;
      }
   }
   {
      GLuint i;
      GLfloat m2 = m[2],  m6 = m[6],  m10 = m[10],  m14 = m[14];
      GLfloat m3 = m[3],  m7 = m[7],  m11 = m[11],  m15 = m[15];
      if (m3==0.0F && m7==0.0F && m11==0.0F && m15==1.0F) {
         /* common case */
         for (i=0;i<n;i++) {
            GLfloat p0 = p[i][0], p1 = p[i][1], p2 = p[i][2];
            q[i][2] = m2 * p0 + m6 * p1 + m10 * p2 + m14;
				/*				q[i][3] = 1.0F;*/
         }
      }
      else {
         /* general case */
         for (i=0;i<n;i++) {
            GLfloat p0 = p[i][0], p1 = p[i][1], p2 = p[i][2];
            q[i][2] = m2 * p0 + m6 * p1 + m10 * p2 + m14;
				/*				q[i][3] = m3 * p0 + m7 * p1 + m11 * p2 + m15; */
         }
      }
   }
}

void RayTransformNormals33( GLuint n, GLfloat3 *q, const GLfloat m[16],GLfloat3 *p )
{
   {
      GLuint i;
      GLfloat m0 = m[0],  m4 = m[4],  m8 = m[8];
      GLfloat m1 = m[1],  m5 = m[5],  m9 = m[9];
      for (i=0;i<n;i++) {
         GLfloat p0 = p[i][0], p1 = p[i][1], p2 = p[i][2];
         q[i][0] = m0 * p0 + m4  * p1 + m8 * p2;
         q[i][1] = m1 * p0 + m5  * p1 + m9 * p2;
      }
   }
   {
      GLuint i;
      GLfloat m2 = m[2],  m6 = m[6],  m10 = m[10];
      GLfloat m3 = m[3],  m7 = m[7],  m11 = m[11],  m15 = m[15];
      if (m3==0.0F && m7==0.0F && m11==0.0F && m15==1.0F) {
         /* common case */
         for (i=0;i<n;i++) {
            GLfloat p0 = p[i][0], p1 = p[i][1], p2 = p[i][2];
            q[i][2] = m2 * p0 + m6 * p1 + m10 * p2;
         }
      }
      else {
         /* general case */
         for (i=0;i<n;i++) {
            GLfloat p0 = p[i][0], p1 = p[i][1], p2 = p[i][2];
            q[i][2] = m2 * p0 + m6 * p1 + m10 * p2;
         }
      }
   }
	{      
	  GLuint i;
	  for (i=0;i<n;i++) { /* renormalize - can we do this to the matrix instead? */
		 normalize3f(q[i]);
	  }
	}
}



