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

#include<stdlib.h>
#include<stdio.h>
#include<math.h>
#include<values.h>

#include"Base.h"
#include"Triangle.h"
#include"Map.h"
#include"MemoryDebug.h"
#include"Err.h"

typedef struct {
  int index;
  int value;
  int next;
} LinkType;

typedef struct {
  int *active; /* active edges */
  int nActive;
  int *edgeStatus;
  int *tri;
  int nTri;
  float *vBlock; /* direction vector for blocked side of an active edge */
  int nVBlock;
  MapType *map;
  LinkType *link;
  int nLink;
} TriangleSurfaceRec;

TriangleSurfaceRec TriangleSurface;

static int TriangleEdgeStatus(int i1,int i2) 
{
  TriangleSurfaceRec *I=&TriangleSurface;
  int l,low,high;
  low = ( i1>i2 ? i2 : i1 );
  high = ( i1>i2 ? i1 : i2 );

  l=I->edgeStatus[low]; 
  while(l) {
	 if(I->link[l].index == high) 
		return(I->link[l].value); /* -1 closed, 0 open, >0 then has index for blocking vector */
	 l=I->link[l].next;
  }
  return(0);
}

static void TriangleEdgeSetStatus(int i1,int i2,int value) 
{
  TriangleSurfaceRec *I=&TriangleSurface;
  int l,low,high;
  low = ( i1>i2 ? i2 : i1 );
  high = ( i1>i2 ? i1 : i2 );

  /*  printf("set: %i %i %i\n",i1,i2,value);*/
  l=I->edgeStatus[low]; 
  while(l) {
	 if(I->link[l].index == high) {
		I->link[l].value = value; 
		return;
	 }
	 l=I->link[l].next;
  }
  if(!l) {
	 VLACheck(I->link,LinkType,I->nLink);
	 I->link[I->nLink].next = I->edgeStatus[low];
	 I->edgeStatus[low] = I->nLink;
	 /*	 printf("offset %i value %i index %i\n",I->nLink,value,high);*/
	 I->link[I->nLink].index = high;
	 I->link[I->nLink].value = value; 
	 I->nLink++;
  }
}

static void TriangleAdd(int i0,int i1,int i2,float *v,float *vn)
{
  TriangleSurfaceRec *I=&TriangleSurface;
  int s01,s12,s02;
  float *v0,*v1,*v2,*n0,*n1,*n2;
  float vt[3],vt2[3],nt[3];

  v0 = v+3*i0;
  v1 = v+3*i1;
  v2 = v+3*i2;
  s01 = TriangleEdgeStatus(i0,i1);
  s02 = TriangleEdgeStatus(i0,i2);
  s12 = TriangleEdgeStatus(i1,i2);
  n0 = vn+3*i0;
  n1 = vn+3*i1;
  n2 = vn+3*i2;
  /* create a new triangle */
  VLACheck(I->tri,int,(I->nTri*3)+2);
  I->tri[I->nTri*3] = i0;
  I->tri[I->nTri*3+1] = i1;
  I->tri[I->nTri*3+2] = i2;
  I->nTri++;
  /*				printf("creating %i %i %i\n",i0,i1,i2);*/
  if(s01) {
	 TriangleEdgeSetStatus(i0,i1,-1);
  } else {
	 VLACheck(I->vBlock,float,(I->nVBlock*3)+2);
	 subtract3f(v0,v1,nt);
	 normalize3f(nt);
	 subtract3f(v2,v1,vt);
	 remove_component3f(vt,nt,I->vBlock+I->nVBlock*3);
	 TriangleEdgeSetStatus(i0,i1,I->nVBlock);
	 I->nVBlock++;
  }
  if(s02) {
	 TriangleEdgeSetStatus(i0,i2,-1);
  } else {
	 VLACheck(I->vBlock,float,(I->nVBlock*3)+2);
	 subtract3f(v0,v2,nt);
	 normalize3f(nt);
	 subtract3f(v1,v0,vt);
	 remove_component3f(vt,nt,I->vBlock+I->nVBlock*3);
	 TriangleEdgeSetStatus(i0,i2,I->nVBlock);
	 I->nVBlock++;
  }
  if(s12) {
	 TriangleEdgeSetStatus(i1,i2,-1);
  } else {
	 VLACheck(I->vBlock,float,(I->nVBlock*3)+2);
	 subtract3f(v1,v2,nt);
	 normalize3f(nt);
	 subtract3f(v0,v1,vt);
	 remove_component3f(vt,nt,I->vBlock+I->nVBlock*3);
	 TriangleEdgeSetStatus(i1,i2,I->nVBlock);
	 I->nVBlock++;
  }
}

int *TrianglePointsToSurface(float *v,float *vn,int n,float cutoff,int *nTriPtr)
{
  TriangleSurfaceRec *I=&TriangleSurface;
  int a;
  MapType *map;
  float *v1,*v2,*v0,vt[3],vt2[3],vt3[3],vt4[3],nt[3],*n0,*n1,*n2;
  int i0,i1,i2,s01,s02,s12,i,j,h,k,l;
  int doneFlag = false;
  float dif,minDist,d1,d2,dp;
  int flag;

  I->nActive = 0;
  I->active=VLAlloc(int,1000);

  I->link=VLAlloc(LinkType,n*2);
  I->nLink = 1;

  I->vBlock=VLAlloc(float,n);
  I->nVBlock = 1;

  I->tri=VLAlloc(int,n);
  I->nTri = 0;

  I->map=MapNew(cutoff,v,n,NULL);
  MapSetupExpress(I->map);
  map=I->map;
  
  I->edgeStatus = Alloc(int,n);
  for(a=0;a<n;a++) {
	 I->edgeStatus[a]=0;
  }
  /* first, let's hit all the easy triangles */

#ifdef _ASDFSD

  for(a=0;a<n;a++) 
	 {
		v0=v+a*3;
		MapLocus(map,v0,&h,&k,&l);
		i=*(MapEStart(map,h,k,l));
		if(i) {
		  j=map->EList[i++];
		  while(j>=0) {
			 if(j!=a) 
				{
				  dif=diff3f(v+3*j,v0);
				  if(dif<minDist)
					 if(TriangleEdgeStatus(a,j)>=0) /* can we put a triangle here? */
						{
						  minDist = dif;
						  i1=a;
						  i2=j;
						  doneFlag=false;
						}
				}
			 j=map->EList[i++];
		  }
		}
	 }
  
  while(!doneFlag) {
	 minDist = MAXFLOAT;
	 doneFlag=true;
	 for(a=0;a<n;a++) 
		if(!I->edgeStatus[a])
		  {
			 v0=v+a*3;
			 MapLocus(map,v0,&h,&k,&l);
			 i=*(MapEStart(map,h,k,l));
			 if(i) {
				j=map->EList[i++];
				while(j>=0) {
				  if(j!=a) 
					 {
						dif=diff3f(v+3*j,v0);
						if(dif<minDist)
						  if(TriangleEdgeStatus(a,j)>=0) /* can we put a triangle here? */
							 {
								minDist = dif;
								i1=a;
								i2=j;
								doneFlag=false;
							 }
					 }
				  j=map->EList[i++];
				}
			 }
		  }
	 if(!doneFlag) {
		VLACheck(I->active,int,I->nActive*2+1);
		I->active[I->nActive*2] = i1;
		I->active[I->nActive*2+1] = i2;
		I->nActive++;
		while(I->nActive&&(I->nTri<80)) {
		  I->nActive--;
		  i1 = I->active[I->nActive*2];
		  i2 = I->active[I->nActive*2+1];
		  s12 = TriangleEdgeStatus(i1,i2);
		  /*		  printf("%i i1:%i i2:%i s12:%i\n",I->nActive,i1,i2,s12);*/
		  if(s12>=0) {

			 /* find neighboring points */
			 
			 minDist = MAXFLOAT;
			 i0=-1;
			 v1=v+i1*3;
			 v2=v+i2*3;
			 MapLocus(map,v1,&h,&k,&l);
			 i=*(MapEStart(map,h,k,l));
			 if(i) {
				j=map->EList[i++];
				while(j>=0) {
				  if((j!=i1)&&(j!=i2))
					 {
						v0 = v+3*j;
						d1 = diff3f(v0,v1);
						d2 = diff3f(v0,v2);
						dif= ( d2 > d1 ? d2 : d1 );
						if(dif<minDist)
						  {
							 s01 = TriangleEdgeStatus(j,i1);
							 s02 = TriangleEdgeStatus(j,i2);
							 /*							 printf("%i %8.3f %i(%i-%i) %i(%i-%i)\n",j,dif,s01,j,i1,s02,j,i2);*/

							 if((s01>=0)&&(s02>=0))
								{

								  flag=true;

								  if(s12>0) {
									 subtract3f(v0,v1,vt);
									 if(dot_product3f(I->vBlock+s12*3,vt)>0.0)
										flag=false;
								  }
								  if(s01>0) {
									 subtract3f(v2,v0,vt);
									 if(dot_product3f(I->vBlock+s01*3,vt)>0.0)
										flag=false;
								  }
								  if(s02>0) {
									 subtract3f(v1,v0,vt);
									 if(dot_product3f(I->vBlock+s02*3,vt)>0.0)
										flag=false;
								  }
								  if(flag) {
									 minDist = dif;
									 i0=j; 
								  }
								}
						  }
					 }
				  j=map->EList[i++];
				}
			 }
			 if(i0<0) {
				TriangleEdgeSetStatus(i1,i2,-2);				
			 } else {
				TriangleAdd(i0,i1,i2,v,vn);
			 }
		  }
		}
	 }
	 break;
  }
#endif

  printf("NTri: %i\n",I->nTri);
  (*nTriPtr)=I->nTri;
  VLAFreeP(I->active);
  VLAFreeP(I->link);
  VLAFreeP(I->vBlock);
  FreeP(I->edgeStatus);
  MapFree(map);
  return(I->tri);
  }



