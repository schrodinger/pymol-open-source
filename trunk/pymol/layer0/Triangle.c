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
#include"Setting.h"

float TestLine[10000];
int NTestLine = 0;

typedef struct {
  int index;
  int value;
  int next;
} LinkType;

typedef struct {
  int *active; /* active edges */
  int nActive;
  int *edgeStatus;
  int *vertActive;
  int *vertWeight;
  int *tri;
  int nTri;
  float *vNormal; /* normal vector for first triangle of an active edge */
  int *vThird;
  int nEdge;
  MapType *map;
  LinkType *link;
  int nLink;
  int N;
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


static int TriangleActivateEdges(int low)
{
  TriangleSurfaceRec *I=&TriangleSurface;
  int l;
  l=I->edgeStatus[low]; 
  while(l) {
    if(I->link[l].value>0) {
      VLACheck(I->active,int,I->nActive*2+1);
      I->active[I->nActive*2] = low;
      I->active[I->nActive*2+1] = I->link[l].index;
      I->nActive++;
    }
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

#define max_edge_len 2.5

static void AddActive(int i1,int i2) {
  int t;

  TriangleSurfaceRec *I=&TriangleSurface;
  if(i1>i2) {
    t=i1;
    i1=i2;
    i2=t;
  }
  VLACheck(I->active,int,I->nActive*2+1);
  I->active[I->nActive*2] = i1;
  I->active[I->nActive*2+1] = i2;
  I->nActive++;
  if(I->vertActive[i1]<0) I->vertActive[i1]=0;
  I->vertActive[i1]++;
  if(I->vertActive[i2]<0) I->vertActive[i2]=0;
  I->vertActive[i2]++;
}
static void TriangleAdd(int i0,int i1,int i2,float *tNorm,float *v,float *vn)
{
  TriangleSurfaceRec *I=&TriangleSurface;
  int s01,s12,s02;
  float *v0,*v1,*v2,*n0,*n1,*n2;
  float vt[3],vt1[3],vt2[3],vt3[3],nt[3];

  v0 = v+3*i0;
  v1 = v+3*i1;
  v2 = v+3*i2;
  s01 = TriangleEdgeStatus(i0,i1);
  s02 = TriangleEdgeStatus(i0,i2);
  s12 = TriangleEdgeStatus(i1,i2);
  n0 = vn+3*i0;
  n1 = vn+3*i1;
  n2 = vn+3*i2;

  /* now, bend the normals a bit so that they line up better with the
     actual triangles drawn */
  add3f(n0,n1,vt3);
  add3f(n2,vt3,vt3);
  I->vertWeight[i0]++;
  scale3f(n0,I->vertWeight[i0],n0);
  add3f(tNorm,n0,n0);
  normalize3f(n0);
  I->vertWeight[i1]++;
  scale3f(n1,I->vertWeight[i1],n1);
  add3f(tNorm,n1,n1);
  normalize3f(n1);
  I->vertWeight[i2]++;
  scale3f(n2,I->vertWeight[i2],n2);
  add3f(tNorm,n2,n2);
  normalize3f(n2);
  
  /* create a new triangle */
  VLACheck(I->tri,int,(I->nTri*3)+2);
  I->tri[I->nTri*3] = i0;
  I->tri[I->nTri*3+1] = i1;
  I->tri[I->nTri*3+2] = i2;
  I->nTri++;
  /*  	printf("creating %i %i %i\n",i0,i1,i2);*/
  if(s01) {
	 TriangleEdgeSetStatus(i0,i1,-1);
    I->vertActive[i0]--; /* deactivate when all active edges are closed */
    I->vertActive[i1]--;

  } else {
	 VLACheck(I->vThird,int,I->nEdge);
	 I->vThird[I->nEdge] = i2;
	 VLACheck(I->vNormal,float,(I->nEdge*3)+2);
	 copy3f(tNorm,I->vNormal+I->nEdge*3);
	 TriangleEdgeSetStatus(i0,i1,I->nEdge);
	 I->nEdge++;
    AddActive(i0,i1);
  }
  if(s02) {
	 TriangleEdgeSetStatus(i0,i2,-1);
    I->vertActive[i0]--; /* deactivate when all active edges are closed */
    I->vertActive[i2]--;

  } else {
	 VLACheck(I->vThird,int,I->nEdge);
	 I->vThird[I->nEdge] = i1;
	 VLACheck(I->vNormal,float,(I->nEdge*3)+2);
	 copy3f(tNorm,I->vNormal+I->nEdge*3);
	 TriangleEdgeSetStatus(i0,i2,I->nEdge);
	 I->nEdge++;
    AddActive(i0,i2);
  }
  if(s12) {
	 TriangleEdgeSetStatus(i1,i2,-1); 
    I->vertActive[i1]--; /* deactivate when all active edges are closed */
    I->vertActive[i2]--;

  } else {
	 VLACheck(I->vThird,int,I->nEdge);
	 I->vThird[I->nEdge] = i0;
	 VLACheck(I->vNormal,float,(I->nEdge*3)+2);
	 copy3f(tNorm,I->vNormal+I->nEdge*3);
	 TriangleEdgeSetStatus(i1,i2,I->nEdge);
	 I->nEdge++;
    AddActive(i1,i2);
  }
}

static void TriangleBuildMinimum(int i1,int i2,float *v,float *vn,int n)
{
  /* this routine builds small, easy triangles... */

  TriangleSurfaceRec *I=&TriangleSurface;
  MapType *map;
  float *v1,*v2,*v0,*v4,vt[3],vt1[3],vt2[3],vt3[3],vt4[3],*n0,*n1,*n2,tNorm[3];
  int i0,s01,s02,s12,i,j,h,k,l,i4;
  float dif,minDist,d1,d2,dp;
  int flag;
  int used = -1;

  map=I->map;
  s12 = TriangleEdgeStatus(i1,i2);
  if(s12>0) used = I->vThird[s12];
  /*  printf("BM: %i %i %i\n",i1,i2,used);*/

  if(s12>=0) {
    minDist = MAXFLOAT;
    i0=-1;
    v1=v+i1*3; v2=v+i2*3;
    MapLocus(map,v1,&h,&k,&l);
	 i=*(MapEStart(map,h,k,l));
    if(i) {
      j=map->EList[i++];
      while(j>=0) {
        if((j!=i1)&&(j!=i2)&&(j!=used))
          {
            v0 = v+3*j;
            d1 = diff3f(v0,v1);
            d2 = diff3f(v0,v2);
            dif= ( d2 > d1 ? d2 : d1 );
				if(dif<minDist)
				  {
					 minDist = dif;
					 i0=j; 
				  }
			 }
		  j=map->EList[i++];
		}
		if(i0>=0) {
		  s01 = TriangleEdgeStatus(i0,i1); s02 = TriangleEdgeStatus(i0,i2);
		  if(I->vertActive[i0]>0) {
			 if(!((s01>0)||(s02>0)))
				i0=-1; /* don't allow non-adjacent joins to active vertices */
		  }
		}
		if(i0>=0) {
		  v0 = v+3*i0;
		  flag=false;
		  if(I->vertActive[i0]) {
			 if((s01>=0)&&(s02>=0)) flag=true;
			 if(flag) { /* are all normals pointing in generally the same direction? */
				n0 = vn+3*i0; n1 = vn+3*i1; n2 = vn+3*i2;							 
				add3f(n0,n1,vt1);
				add3f(n2,vt1,vt2);
				normalize3f(vt2);
				if(dot_product3f(n0,vt2)<0.1) flag=false;
				if(dot_product3f(n1,vt2)<0.1) flag=false;
				if(dot_product3f(n2,vt2)<0.1) flag=false;
			 } /*printf("pass normal sums %i\n",flag);*/
			 if(flag) { /* does the sum of the normals point in the same direction as the triangle? */
				subtract3f(v1,v0,vt3);
			 subtract3f(v2,v0,vt4);
			 cross_product3f(vt3,vt4,tNorm); 
			 normalize3f(tNorm); 							 
			 dp = dot_product3f(vt2,tNorm);
			 if(dp<0) scale3f(tNorm,-1.0,tNorm);
			 if(fabs(dp)<0.1) flag = false;
			 } /*printf("pass tNorm  %i\n",flag);*/
			 if(flag) {
				if(s12>0) if(dot_product3f(I->vNormal+s12*3,tNorm)<0.1) flag=false; 
				if(s01>0) if(dot_product3f(I->vNormal+s01*3,tNorm)<0.1) flag=false; 
				if(s02>0) if(dot_product3f(I->vNormal+s02*3,tNorm)<0.1) flag=false; 
			 } /*printf("pass compare tNorm %i\n",flag);*/
			 if(flag) { /* are all the Blocking vectors pointing outward, and are the triangle normals consistent? */
				if(s12>0) {
				  i4 = I->vThird[s12];
				  v4=v+i4*3;
				  subtract3f(v0,v1,vt1);
				  subtract3f(v4,v1,vt2);
				  subtract3f(v1,v2,vt);
				  normalize3f(vt);
				  remove_component3f(vt1,vt,vt3);
				  remove_component3f(vt2,vt,vt4);
				  normalize3f(vt3);
				  normalize3f(vt4);
				  if(dot_product3f(vt3,vt4)>0.0) flag=false;
				}			 
				if(s01>0) {
				  i4 = I->vThird[s01];
				  v4=v+i4*3;
				  subtract3f(v2,v0,vt1);
				  subtract3f(v4,v0,vt2);
				  subtract3f(v0,v1,vt);
				  normalize3f(vt);
				  remove_component3f(vt1,vt,vt3);
				  remove_component3f(vt2,vt,vt4);
				  normalize3f(vt3);
				  normalize3f(vt4);
				  if(dot_product3f(vt3,vt4)>0.0) flag=false;
				}
				if(s02>0) {
				  i4 = I->vThird[s02];
				  v4=v+i4*3;
				  subtract3f(v1,v0,vt1);
				  subtract3f(v4,v0,vt2);
				  subtract3f(v0,v2,vt);
				  normalize3f(vt);
				  remove_component3f(vt1,vt,vt3);
				  remove_component3f(vt2,vt,vt4);
				  normalize3f(vt3);
				  normalize3f(vt4);
				  if(dot_product3f(vt3,vt4)>0.0) flag=false;
				}
			 } /*printf("pass blocking %i\n",flag);*/
		  }
		  if(flag) TriangleAdd(i0,i1,i2,tNorm,v,vn);
		}
	 }
  }
}

static void TriangleBuildSecond(int i1,int i2,float *v,float *vn,int n)
{
  /* this routine builds small, easy triangles... */

  TriangleSurfaceRec *I=&TriangleSurface;
  MapType *map;
  float *v1,*v2,*v0,*v4,vt[3],vt1[3],vt2[3],vt3[3],vt4[3],*n0,*n1,*n2,tNorm[3];
  int i0,s01,s02,s12,i,j,h,k,l,i4;
  float dif,minDist,d1,d2,dp;
  int flag;
  int used = -1;

  map=I->map;
  s12 = TriangleEdgeStatus(i1,i2);
  if(s12>0) used = I->vThird[s12];
  if(s12>=0) {
    minDist = MAXFLOAT;
    i0=-1;
    v1=v+i1*3; v2=v+i2*3;
    MapLocus(map,v1,&h,&k,&l);
	 i=*(MapEStart(map,h,k,l));
    if(i) {
      j=map->EList[i++];
      while(j>=0) {
        if((j!=i1)&&(j!=i2)&&(j!=used)&&(I->vertActive[j]))
          {
            v0 = v+3*j;
            d1 = diff3f(v0,v1);
            d2 = diff3f(v0,v2);
            dif= ( d2 > d1 ? d2 : d1 );
				if(dif<minDist)
				  {
					 minDist = dif;
					 i0=j; 
				  }
			 }
		  j=map->EList[i++];
		}
		if(i0>=0) {
		  s01 = TriangleEdgeStatus(i0,i1); s02 = TriangleEdgeStatus(i0,i2);
		  if(I->vertActive[i0]>0) {
			 if(!((s01>0)||(s02>0)))
				i0=-1; /* don't allow non-adjacent joins to active vertices */
		  }
		}
		if(i0>=0) {
		  v0 = v+3*i0;
		  flag=false;
		  if(I->vertActive[i0]) {
			 if((s01>=0)&&(s02>=0)) flag=true;
			 if(flag) { /* are all normals pointing in generally the same direction? */
				n0 = vn+3*i0; n1 = vn+3*i1; n2 = vn+3*i2;							 
				add3f(n0,n1,vt1);
				add3f(n2,vt1,vt2);
				normalize3f(vt2);
				if(dot_product3f(n0,vt2)<0.1) flag=false;
				if(dot_product3f(n1,vt2)<0.1) flag=false;
				if(dot_product3f(n2,vt2)<0.1) flag=false;
			 } /*printf("pass normal sums %i\n",flag);*/
			 if(flag) { /* does the sum of the normals point in the same direction as the triangle? */
				subtract3f(v1,v0,vt3);
				subtract3f(v2,v0,vt4);
				cross_product3f(vt3,vt4,tNorm); 
				normalize3f(tNorm); 							 
				dp = dot_product3f(vt2,tNorm);
				if(dp<0) scale3f(tNorm,-1.0,tNorm);
				if(fabs(dp)<0.1) flag = false;
			 } /*printf("pass tNorm  %i\n",flag);*/
			 if(flag) {
				if(s12>0) if(dot_product3f(I->vNormal+s12*3,tNorm)<0.1) flag=false; 
				if(s01>0) if(dot_product3f(I->vNormal+s01*3,tNorm)<0.1) flag=false; 
				if(s02>0) if(dot_product3f(I->vNormal+s02*3,tNorm)<0.1) flag=false; 
			 } /*printf("pass compare tNorm %i\n",flag);*/
			 if(flag) { /* are all the Blocking vectors pointing outward, and are the triangle normals consistent? */
				if(s12>0) {
				  i4 = I->vThird[s12];
				  v4=v+i4*3;
				  subtract3f(v0,v1,vt1);
				  subtract3f(v4,v1,vt2);
				  subtract3f(v1,v2,vt);
				  normalize3f(vt);
				  remove_component3f(vt1,vt,vt3);
				  remove_component3f(vt2,vt,vt4);
				  normalize3f(vt3);
				  normalize3f(vt4);
				  if(dot_product3f(vt3,vt4)>0.0) flag=false;
				}			 
				if(s01>0) {
				  i4 = I->vThird[s01];
				  v4=v+i4*3;
				  subtract3f(v2,v0,vt1);
				  subtract3f(v4,v0,vt2);
				  subtract3f(v0,v1,vt);
				  normalize3f(vt);
				  remove_component3f(vt1,vt,vt3);
				  remove_component3f(vt2,vt,vt4);
				  normalize3f(vt3);
				  normalize3f(vt4);
				  if(dot_product3f(vt3,vt4)>0.0) flag=false;
				}
				if(s02>0) {
				  i4 = I->vThird[s02];
				  v4=v+i4*3;
				  subtract3f(v1,v0,vt1);
				  subtract3f(v4,v0,vt2);
				  subtract3f(v0,v2,vt);
				  normalize3f(vt);
				  remove_component3f(vt1,vt,vt3);
				  remove_component3f(vt2,vt,vt4);
				  normalize3f(vt3);
				  normalize3f(vt4);
				  if(dot_product3f(vt3,vt4)>0.0) flag=false;
				}
			 } /*printf("pass blocking %i\n",flag);*/
		  }
		  if(flag) TriangleAdd(i0,i1,i2,tNorm,v,vn);
		}
	 }
  }
}

static void TriangleBuildSingle(int i1,int i2,float *v,float *vn,int n)
{
  /* this routine builds small, easy triangles... */

  TriangleSurfaceRec *I=&TriangleSurface;
  MapType *map;
  float *v1,*v2,*v0,*v4,vt[3],vt1[3],vt2[3],vt3[3],vt4[3],*n0,*n1,*n2,tNorm[3];
  int i0,s01,s02,s12,i,j,h,k,l,i4;
  float dif,minDist,d1,d2,dp;
  int flag;
  int used = -1;

  map=I->map;
  s12 = TriangleEdgeStatus(i1,i2);
  if(s12>0) used = I->vThird[s12];
  if(s12>=0) {
    minDist = MAXFLOAT;
    i0=-1;
    v1=v+i1*3; v2=v+i2*3;
    MapLocus(map,v1,&h,&k,&l);
	 i=*(MapEStart(map,h,k,l));
    if(i) {
      j=map->EList[i++];
      while(j>=0) {
        if((j!=i1)&&(j!=i2)&&(j!=used)&&(I->vertActive[j]))
          {
            v0 = v+3*j;
            d1 = diff3f(v0,v1);
            d2 = diff3f(v0,v2);
            dif= ( d2 > d1 ? d2 : d1 );
				if(dif<minDist)
				  {
					 minDist = dif;
					 i0=j; 
				  }
			 }
		  j=map->EList[i++];
		}
		if(i0>=0) {
		  v0 = v+3*i0;
		  flag=false;
		  s01 = TriangleEdgeStatus(i0,i1); s02 = TriangleEdgeStatus(i0,i2);
		  if(I->vertActive[i0]) {
			 if((s01>=0)&&(s02>=0)) flag=true;
			 if(flag) { /* are all normals pointing in generally the same direction? */
				n0 = vn+3*i0; n1 = vn+3*i1; n2 = vn+3*i2;							 
				add3f(n0,n1,vt1);
				add3f(n2,vt1,vt2);
				normalize3f(vt2);
				if(dot_product3f(n0,vt2)<0.1) flag=false;
				if(dot_product3f(n1,vt2)<0.1) flag=false;
				if(dot_product3f(n2,vt2)<0.1) flag=false;
			 } /*printf("pass normal sums %i\n",flag);*/
			 if(flag) { /* does the sum of the normals point in the same direction as the triangle? */
				subtract3f(v1,v0,vt3);
				subtract3f(v2,v0,vt4);
				cross_product3f(vt3,vt4,tNorm); 
				normalize3f(tNorm); 							 
				dp = dot_product3f(vt2,tNorm);
				if(dp<0) scale3f(tNorm,-1.0,tNorm);
				if(fabs(dp)<0.1) flag = false;
			 } /*printf("pass tNorm  %i\n",flag);*/
			 if(flag) {
				if(s12>0) if(dot_product3f(I->vNormal+s12*3,tNorm)<0.1) flag=false; 
				if(s01>0) if(dot_product3f(I->vNormal+s01*3,tNorm)<0.1) flag=false; 
				if(s02>0) if(dot_product3f(I->vNormal+s02*3,tNorm)<0.1) flag=false; 
			 } /*printf("pass compare tNorm %i\n",flag);*/
			 if(flag) { /* are all the Blocking vectors pointing outward, and are the triangle normals consistent? */
				if(s12>0) {
				  i4 = I->vThird[s12];
				  v4=v+i4*3;
				  subtract3f(v0,v1,vt1);
				  subtract3f(v4,v1,vt2);
				  subtract3f(v1,v2,vt);
				  normalize3f(vt);
				  remove_component3f(vt1,vt,vt3);
				  remove_component3f(vt2,vt,vt4);
				  normalize3f(vt3);
				  normalize3f(vt4);
				  if(dot_product3f(vt3,vt4)>0.0) flag=false;
				}			 
				if(s01>0) {
				  i4 = I->vThird[s01];
				  v4=v+i4*3;
				  subtract3f(v2,v0,vt1);
				  subtract3f(v4,v0,vt2);
				  subtract3f(v0,v1,vt);
				  normalize3f(vt);
				  remove_component3f(vt1,vt,vt3);
				  remove_component3f(vt2,vt,vt4);
				  normalize3f(vt3);
				  normalize3f(vt4);
				  if(dot_product3f(vt3,vt4)>0.0) flag=false;
				}
				if(s02>0) {
				  i4 = I->vThird[s02];
				  v4=v+i4*3;
				  subtract3f(v1,v0,vt1);
				  subtract3f(v4,v0,vt2);
				  subtract3f(v0,v2,vt);
				  normalize3f(vt);
				  remove_component3f(vt1,vt,vt3);
				  remove_component3f(vt2,vt,vt4);
				  normalize3f(vt3);
				  normalize3f(vt4);
				  if(dot_product3f(vt3,vt4)>0.0) flag=false;
				}
			 } /*printf("pass blocking %i\n",flag);*/
		  }
		  if(flag) TriangleAdd(i0,i1,i2,tNorm,v,vn);
		}
	 }
  }
}


static void TriangleBuildThird(int i1,int i2,float *v,float *vn,int n)
{
  /* this routine builds small, easy triangles... */

  TriangleSurfaceRec *I=&TriangleSurface;
  MapType *map;
  float *v1,*v2,*v0,*v4,vt[3],vt1[3],vt2[3],vt3[3],vt4[3],*n0,*n1,*n2,tNorm[3];
  int i0,s01,s02,s12,i,j,h,k,l,i4;
  float dif,minDist,d1,d2,dp;
  int used = -1;

  map=I->map;
  s12 = TriangleEdgeStatus(i1,i2);
  if(s12>0) used = I->vThird[s12];
  if(s12>=0) {
    minDist = MAXFLOAT;
    i0=-1;
    v1=v+i1*3; v2=v+i2*3;
    MapLocus(map,v1,&h,&k,&l);
	 i=*(MapEStart(map,h,k,l));
    if(i) {
      j=map->EList[i++];
      while(j>=0) {
        if((j!=i1)&&(j!=i2)&&(j!=used)&&(I->vertActive[j]))
          {
            v0 = v+3*j;
            d1 = diff3f(v0,v1);
            d2 = diff3f(v0,v2);
            dif= ( d2 > d1 ? d2 : d1 );
				if(dif<minDist)
				  {
					 minDist = dif;
					 i0=j; 
				  }
			 }
		  j=map->EList[i++];
		}
		if(i0>=0) {
		  v0 = v+3*i0;
		  s01 = TriangleEdgeStatus(i0,i1); s02 = TriangleEdgeStatus(i0,i2);
		  /* if all three edges are active */
		  if((s12>0)&&(s01>0)&&(s02>0)) { 
			 n0 = vn+3*i0; n1 = vn+3*i1; n2 = vn+3*i2;							 
			 add3f(n0,n1,vt1);
			 add3f(n2,vt1,vt2);
			 subtract3f(v1,v0,vt3);
			 subtract3f(v2,v0,vt4);
			 cross_product3f(vt3,vt4,tNorm); 
			 normalize3f(tNorm); 							 
			 dp = dot_product3f(vt2,tNorm);
			 if(dp<0) scale3f(tNorm,-1.0,tNorm);
			 TriangleAdd(i0,i1,i2,tNorm,v,vn);
		  }
		}
	 }
  }
}


static void TriangleBuildLast(int i1,int i2,float *v,float *vn,int n)
{
  /* this routine builds small, easy triangles... */

  TriangleSurfaceRec *I=&TriangleSurface;
  MapType *map;
  float *v1,*v2,*v0,*v4,vt[3],vt1[3],vt2[3],vt3[3],vt4[3],*n0,*n1,*n2,tNorm[3];
  int i0,s01,s02,s12,i,j,h,k,l,i4;
  float dif,minDist,d1,d2,dp;
  int used = -1;

  map=I->map;
  s12 = TriangleEdgeStatus(i1,i2);
  if(s12>0) used = I->vThird[s12];
  if(s12>=0) {
    minDist = MAXFLOAT;
    i0=-1;
    v1=v+i1*3; v2=v+i2*3;
    MapLocus(map,v1,&h,&k,&l);
	 i=*(MapEStart(map,h,k,l));
    if(i) {
      j=map->EList[i++];
      while(j>=0) {
        if((j!=i1)&&(j!=i2)&&(j!=used)&&(I->vertActive[j]))
          {
            v0 = v+3*j;
            d1 = diff3f(v0,v1);
            d2 = diff3f(v0,v2);
            dif= ( d2 > d1 ? d2 : d1 );
				if(dif<minDist)
				  {
					 minDist = dif;
					 i0=j; 
				  }
			 }
		  j=map->EList[i++];
		}
		if(i0>=0) {
		  v0 = v+3*i0;
		  s01 = TriangleEdgeStatus(i0,i1); s02 = TriangleEdgeStatus(i0,i2);
		  /* if all three edges are active */
		  if((s12>0)&&((s01>0)||(s02>0))) { 
			 n0 = vn+3*i0; n1 = vn+3*i1; n2 = vn+3*i2;							 
			 add3f(n0,n1,vt1);
			 add3f(n2,vt1,vt2);
			 subtract3f(v1,v0,vt3);
			 subtract3f(v2,v0,vt4);
			 cross_product3f(vt3,vt4,tNorm); 
			 normalize3f(tNorm); 							 
			 dp = dot_product3f(vt2,tNorm);
			 if(dp<0) scale3f(tNorm,-1.0,tNorm);
			 TriangleAdd(i0,i1,i2,tNorm,v,vn);
		  }
		}
	 }
  }
}



static void FollowActives(float *v,float *vn,int n,int mode)
{
  TriangleSurfaceRec *I=&TriangleSurface;
  int i1,i2;
  int cnt;

  cnt = SettingGet(cSetting_test1);
  
  while(I->nActive&&(I->nTri<cnt)) {
    I->nActive--;
    i1 = I->active[I->nActive*2];
    i2 = I->active[I->nActive*2+1];
	 switch(mode) {
	 case 0:
		TriangleBuildMinimum(i1,i2,v,vn,n);
		break;
	 case 1:
		TriangleBuildSecond(i1,i2,v,vn,n);
		break;
	 case 2:
		TriangleBuildThird(i1,i2,v,vn,n);
		break;
	 case 3:
		TriangleBuildLast(i1,i2,v,vn,n);
		break;
	 }
  }
}

int *TrianglePointsToSurface(float *v,float *vn,int n,float cutoff,int *nTriPtr)
{
  TriangleSurfaceRec *I=&TriangleSurface;
  MapType *map;
  int a,i,j,h,k,l;
  int i1,i2;
  int doneFlag = false;
  int lastTri,lastTri2,lastTri3;
  float dif,minDist,*v0;
  int cnt;

  I->N=n;
  I->nActive = 0;
  I->active=VLAlloc(int,1000);

  I->link=VLAlloc(LinkType,n);
  I->nLink = 1;

  I->nEdge = 1;

  I->vNormal=VLAlloc(float,n);
  I->vThird=VLAlloc(int,n);

  I->tri=VLAlloc(int,n);
  I->nTri = 0;

  I->map=MapNew(cutoff,v,n,NULL);
  MapSetupExpress(I->map);
  map=I->map;
  
  I->edgeStatus = Alloc(int,n);
  for(a=0;a<n;a++) {
	 I->edgeStatus[a]=0;
  }

  I->vertActive = Alloc(int,n);
  for(a=0;a<n;a++) {
	 I->vertActive[a]=-1;
  }

  I->vertWeight = Alloc(int,n);
  for(a=0;a<n;a++) {
	 I->vertWeight[a]=2;
  }

  lastTri3=-1;
  while(lastTri3!=I->nTri) {
	 lastTri3=I->nTri;

	 minDist = MAXFLOAT;
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
							 }
					 }
				  j=map->EList[i++];
				}
			 }
		  }
	 
	 I->nActive=0;
	 VLACheck(I->active,int,I->nActive*2+1);
	 I->active[I->nActive*2] = i1;
	 I->active[I->nActive*2+1] = i2;
	 I->nActive=1;  
	 lastTri=I->nTri;
	 FollowActives(v,vn,n,0);
	 while(lastTri!=I->nTri) {
		lastTri=I->nTri;
		for(a=0;a<n;a++) 
		  if(I->vertActive[a])
			 TriangleActivateEdges(a);
		FollowActives(v,vn,n,0);
	 }
	 
	 lastTri=I->nTri-1;
	 while(lastTri!=I->nTri) {
		lastTri=I->nTri;
		for(a=0;a<n;a++) 
		  if(I->vertActive[a])
			 TriangleActivateEdges(a);
		FollowActives(v,vn,n,1);
	 }
	 
	 lastTri2=I->nTri-1;
	 while(lastTri2!=I->nTri) {
		lastTri2=I->nTri;
		for(a=0;a<n;a++) 
		  if(I->vertActive[a])
			 TriangleActivateEdges(a);
		if(!I->nActive) break;
		lastTri=I->nTri;
		while(I->nTri==lastTri)
		  {
			 if(!I->nActive) break;
			 I->nActive--;
			 i1 = I->active[I->nActive*2];
			 i2 = I->active[I->nActive*2+1];
			 TriangleBuildSingle(i1,i2,v,vn,n);
		  }
		lastTri=I->nTri;
		FollowActives(v,vn,n,1);
		while(lastTri!=I->nTri) {
		  lastTri=I->nTri;
		  for(a=0;a<n;a++) 
			 if(I->vertActive[a])
				TriangleActivateEdges(a);
		  FollowActives(v,vn,n,1);
		}
	 }
	 
	 for(a=0;a<n;a++) 
		if(I->vertActive[a])
		  TriangleActivateEdges(a);
	 FollowActives(v,vn,n,2);

	 lastTri=I->nTri-1;
	 while(lastTri!=I->nTri) {
		lastTri=I->nTri;
		for(a=0;a<n;a++) 
		  if(I->vertActive[a])
			 TriangleActivateEdges(a);
		FollowActives(v,vn,n,3); /* this is a sloppy, forcing tesselation */
	 }

  }
  /*	   for(a=0;a<n;a++) 
		  if(I->vertActive[a])
		  TriangleActivateEdges(a);
		  FollowActives(v,vn,n,2);*/
  /*
      for(a=0;a<n;a++) 
      if(I->vertActive[a])
        printf("active %i\n",a);
    for(a=0;a<n;a++) 
      if(I->vertActive[a])
        TriangleActivateEdges(a);
    FollowActives(v,vn,n,true);
 */
  /*  for(a=0;a<n;a++) 
	 if(I->vertActive[a])
		{
		  TestLine[NTestLine*6]=0;
		  TestLine[NTestLine*6+1]=0;
		  TestLine[NTestLine*6+2]=0;
		  copy3f(v+3*a,TestLine+NTestLine*6+3);
		  NTestLine++;
		  printf("never %i\n",a);
		  }*/
  printf("NTri: %i\n",I->nTri);

  cnt = SettingGet(cSetting_test1);
  if(I->nTri>cnt)
	 I->nTri=cnt;

  (*nTriPtr)=I->nTri;
  VLAFreeP(I->active);
  VLAFreeP(I->link);
  VLAFreeP(I->vNormal);
  VLAFreeP(I->vThird);
  FreeP(I->edgeStatus);
  FreeP(I->vertActive);
  FreeP(I->vertWeight);
  MapFree(map);
  return(I->tri);
}

