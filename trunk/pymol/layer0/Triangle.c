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
static void TriangleAdd(int i0,int i1,int i2,float *v,float *vn)
{
  TriangleSurfaceRec *I=&TriangleSurface;
  int s01,s12,s02;
  float *v0,*v1,*v2,*n0,*n1,*n2;
  float vt[3],vt1[3],vt2[3],tNorm[3],vt3[3],nt[3];

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
  subtract3f(v1,v0,vt1);
  subtract3f(v2,v0,vt2);
  cross_product3f(vt1,vt2,tNorm);
  normalize3f(tNorm);
  add3f(n0,n1,vt3);
  add3f(n2,vt3,vt3);
  if(dot_product3f(tNorm,vt3)<0.0)
    scale3f(tNorm,-1.0,tNorm);
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
	 VLACheck(I->vBlock,float,(I->nVBlock*3)+2);
	 subtract3f(v0,v1,nt);
	 normalize3f(nt);
	 subtract3f(v2,v1,vt);
	 remove_component3f(vt,nt,I->vBlock+I->nVBlock*3);
	 TriangleEdgeSetStatus(i0,i1,I->nVBlock);
	 I->nVBlock++;
    AddActive(i0,i1);
  }
  if(s02) {
	 TriangleEdgeSetStatus(i0,i2,-1);
    I->vertActive[i0]--; /* deactivate when all active edges are closed */
    I->vertActive[i2]--;

  } else {
	 VLACheck(I->vBlock,float,(I->nVBlock*3)+2);
	 subtract3f(v0,v2,nt);
	 normalize3f(nt);
	 subtract3f(v1,v0,vt);
	 remove_component3f(vt,nt,I->vBlock+I->nVBlock*3);
	 TriangleEdgeSetStatus(i0,i2,I->nVBlock);
	 I->nVBlock++;
    AddActive(i0,i2);
  }
  if(s12) {
	 TriangleEdgeSetStatus(i1,i2,-1); 
    I->vertActive[i1]--; /* deactivate when all active edges are closed */
    I->vertActive[i2]--;

  } else {
	 VLACheck(I->vBlock,float,(I->nVBlock*3)+2);
	 subtract3f(v1,v2,nt);
	 normalize3f(nt);
	 subtract3f(v0,v1,vt);
	 remove_component3f(vt,nt,I->vBlock+I->nVBlock*3);
	 TriangleEdgeSetStatus(i1,i2,I->nVBlock);
	 I->nVBlock++;
    AddActive(i1,i2);
  }
}


static void TriangleBuild(int i1,int i2,float *v,float *vn,int n,int force)
{
  TriangleSurfaceRec *I=&TriangleSurface;
  MapType *map;
  float *v1,*v2,*v0,vt[3],vt2[3],vt3[3],vt4[3],*n0,*n1,*n2;
  int i0,s01,s02,s12,i,j,h,k,l,mMinI0,mMinSeal;
  float dif,minDist,mMinDist,d1,d2;
  int flag;
  int c;

  map=I->map;
  s12 = TriangleEdgeStatus(i1,i2);
  /*          printf("%i i1:%i i2:%i s12:%i\n",I->nActive,i1,i2,s12);*/
  if(s12>=0) {
    /*            printf("n1: %8.3f %8.3f %8.3f n2: %8.3f %8.3f %8.3f\n",
                  vn[3*i1],vn[3*i1+1],vn[3*i1+2],
                  vn[3*i2],vn[3*i1+2],vn[3*i2+2]);*/
    
    /* find neighboring points */
    
    minDist = MAXFLOAT;
    mMinDist = MAXFLOAT;
    mMinSeal = 0;
    i0=-1;
    v1=v+i1*3;
    v2=v+i2*3;
    MapLocus(map,v1,&h,&k,&l);
    i=*(MapEStart(map,h,k,l));
    if(i) {
      j=map->EList[i++];
      while(j>=0) {
        if((j!=i1)&&(j!=i2)&&(I->vertActive[j]))
          {
            v0 = v+3*j;
            d1 = diff3f(v0,v1);
            d2 = diff3f(v0,v2);
            dif= ( d2 > d1 ? d2 : d1 );
            if(dif<max_edge_len) {
              if(dif<mMinDist) {
                s01 = TriangleEdgeStatus(j,i1);
                s02 = TriangleEdgeStatus(j,i2);
                if((s12>0)&&((s01>0)||(s02>0))) {
                  subtract3f(v0,v1,vt);
                  if(dot_product3f(I->vBlock+s12*3,vt)<0.0) {
                    if(s01>0) {
                      subtract3f(v2,v0,vt);
                      if(dot_product3f(I->vBlock+s01*3,vt)<0.0) {
                        mMinSeal=1;
                        mMinDist=dif;
                        mMinI0=j;
                      } 
                    }
                    if(s02>0) {
                      subtract3f(v1,v0,vt);
                      if(dot_product3f(I->vBlock+s02*3,vt)<0.0) {
                        mMinSeal=1;
                        mMinDist=dif;
                        mMinI0=j;
                      }
                    }
                  }
                }
              }
              if(dif<minDist)
                {
                  s01 = TriangleEdgeStatus(j,i1);
                  s02 = TriangleEdgeStatus(j,i2);
                  /*                    printf("%i %8.3f %i(%i-%i) %i(%i-%i)\n",j,dif,s01,j,i1,s02,j,i2);*/
                  if((s12>0)&&(s01>0)&&(s02>0)) { /* if all three edges are active */
                    
                    /*                        subtract3f(v0,v1,vt);                        
                                              printf("%8.3f ",dot_product3f(I->vBlock+s12*3,vt));
                                              subtract3f(v2,v0,vt);
                                              printf("%8.3f ",dot_product3f(I->vBlock+s01*3,vt));
                                              subtract3f(v1,v0,vt);
                                              printf("%8.3f\n",dot_product3f(I->vBlock+s02*3,vt));*/
                    
                    if(force) {
                      c=0;
                      subtract3f(v0,v1,vt);
                      if(dot_product3f(I->vBlock+s12*3,vt)<0.0) c++;
                      subtract3f(v2,v0,vt);
                      if(dot_product3f(I->vBlock+s01*3,vt)<0.0) c++;
                      subtract3f(v1,v0,vt);
                      if(dot_product3f(I->vBlock+s02*3,vt)<0.0) c++;
                      if(c>=2) {
                        i0=j;
                        mMinSeal=0;
                        break;
                      }
                    } else {
                      subtract3f(v0,v1,vt);
                      if(dot_product3f(I->vBlock+s12*3,vt)<0.0) {
                        subtract3f(v2,v0,vt);
                        if(dot_product3f(I->vBlock+s01*3,vt)<0.0) {
                          subtract3f(v1,v0,vt);
                          if(dot_product3f(I->vBlock+s02*3,vt)<0.0) {
                            i0=j;
                            mMinSeal=0;
                            break;
                          }
                        }
                      }
                    }
                  }
                  if((s01>=0)&&(s02>=0))
                    {
                      
                      flag=true;
                      
                      n0 = vn+3*j;
                      n1 = vn+3*i1;
                      n2 = vn+3*i2;
                      
                      if(!force) {
                        /* are all normals pointing in generally the same direction? */
                        
                        add3f(n0,n1,vt);
                        add3f(n2,vt,vt2);
                        normalize3f(vt2);
                        
                        if(dot_product3f(n0,vt2)<0.2)
                          flag=false;
                        if(dot_product3f(n1,vt2)<0.2)
                          flag=false;
                        if(dot_product3f(n2,vt2)<0.2)
                          flag=false;
                        /*                         printf("samedir %i\n",flag);*/
                        
                        /* ...and is that direction consistent with plane of the triangle? */
                        
                        
                        subtract3f(v1,v0,vt);
                        subtract3f(v2,v0,vt3);
                        
                        
                        /*
                          printf("%i %i %i\n",j,i1,i2);
                          printf("vt: %8.3f %8.3f %8.3f vt3: %8.3f %8.3f %8.3f\n",
                          vt[0],vt[1],vt[2],
                          vt3[0],vt3[1],vt3[2]);*/
                        
                        cross_product3f(vt,vt3,vt4); 
                        normalize3f(vt4); 
                        if(fabs(dot_product3f(vt2,vt4))<0.2)
                          flag=false;
                        
                        
                        /*
                          printf("n1: %8.3f %8.3f %8.3f n2: %8.3f %8.3f %8.3f\n",
                          vt2[0],vt2[1],vt2[2],
                          vt4[0],vt4[1],vt4[2]);
                          
                          printf("dot %8.3f\n",fabs(dot_product3f(vt2,vt4)));
                          printf("plane check %i\n",flag);*/
                      }
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
                      
                      /*                            printf("overlap check %i\n",flag); */
                      
                      if(flag) {
                        minDist = dif;
                        i0=j; 
                      } else if(i0>=0)
                        i0=-1;
                    }
                }
            }
          }
        j=map->EList[i++];
      }
    }
    if(i0>=0) {
      if(mMinSeal&&(!force)) 
        if((mMinDist*1.1)<minDist) {
          i0=-1;
        }
    }
    if(i0>=0) {
      TriangleAdd(i0,i1,i2,v,vn);
    }
  }
}

static void FollowActives(float *v,float *vn,int n,int force)
{
  TriangleSurfaceRec *I=&TriangleSurface;
  int cnt;
  int i1,i2;

  cnt = SettingGet(cSetting_test1);
  
  while(I->nActive&&(I->nTri<cnt)) {
    I->nActive--;
    i1 = I->active[I->nActive*2];
    i2 = I->active[I->nActive*2+1];
    TriangleBuild(i1,i2,v,vn,n,force);
  }
}

int *TrianglePointsToSurface(float *v,float *vn,int n,float cutoff,int *nTriPtr)
{
  TriangleSurfaceRec *I=&TriangleSurface;
  MapType *map;
  int a,i,j,h,k,l;
  int i1,i2;
  int doneFlag = false;
  float dif,minDist,*v0;

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

  I->vertActive = Alloc(int,n);
  for(a=0;a<n;a++) {
	 I->vertActive[a]=-1;
  }

  I->vertWeight = Alloc(int,n);
  for(a=0;a<n;a++) {
	 I->vertWeight[a]=3;
  }

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
  
  FollowActives(v,vn,n,false);
  for(a=0;a<n;a++) 
    if(I->vertActive[a])
      TriangleActivateEdges(a);
  FollowActives(v,vn,n,false);
    for(a=0;a<n;a++) 
      if(I->vertActive[a])
        printf("active %i\n",a);
    for(a=0;a<n;a++) 
      if(I->vertActive[a])
        TriangleActivateEdges(a);
    FollowActives(v,vn,n,true);
    for(a=0;a<n;a++) 
      if(I->vertActive[a])
        printf("never %i\n",a);
  
  printf("NTri: %i\n",I->nTri);
  (*nTriPtr)=I->nTri;
  VLAFreeP(I->active);
  VLAFreeP(I->link);
  VLAFreeP(I->vBlock);
  FreeP(I->edgeStatus);
  FreeP(I->vertActive);
  FreeP(I->vertWeight);
  MapFree(map);
  return(I->tri);
}




                            /* make sure surface points in the right direction */
                            
                            /* do we have general agreement about the direction? */                          
/*                            
add3f(n1,n0,vt);
add3f(n1,n2,vt2);
if(dot_product3f(vt,vt2)<0.0)
flag=false;

add3f(n0,n2,vt);
add3f(n1,n2,vt2);
if(dot_product3f(vt,vt2)<0.0)
flag=false;

add3f(n1,n0,vt);
add3f(n0,n2,vt2);
if(dot_product3f(vt,vt2)<0.0)
flag=false;

*/                            
                            /* make sure surface isn't orthogonal to any normal */                          
                            /*
                              subtract3f(v1,v0,vt);
                              subtract3f(v2,v0,vt3);
                              cross_product3f(vt,vt3,vt2);
                              normalize3f(vt2);
                              
                              if(fabs(dot_product3f(n0,vt2))<0.1) 
                              flag=false;
                              if(fabs(dot_product3f(n1,vt2))<0.1)
                              flag=false;
                              if(fabs(dot_product3f(n2,vt2))<0.1)
                              flag=false;
                            */
                            /*                          printf("n: %8.3f\n",length3f(vt2));*/
                            
