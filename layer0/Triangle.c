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

/* glossary:

	active edge = an edge with one triangle
	closed edge = an edge with two triangles

	closed vertex = one which is surrounded by a cycle of triangles
	active vertex = (opposite of closed)

*/

#include"os_predef.h"
#include"os_std.h"

#include"Base.h"
#include"Triangle.h"
#include"Map.h"
#include"MemoryDebug.h"
#include"Err.h"
#include"Setting.h"
#include"Ortho.h"
#include"Feedback.h"

typedef struct {
  int index;
  int value;
  int next;
} LinkType;

typedef struct {
  int vert3,tri1;
  int vert4,tri2;
} EdgeRec;

typedef struct {
  PyMOLGlobals *G;
  int *activeEdge; /* active edges */
  int nActive;
  int *edgeStatus;
  int *vertActive;
  int *vertWeight;
  int *tri;
  int nTri;
  float *vNormal; /* normal vector for first triangle of an active edge */
  EdgeRec *edge;
  int nEdge;
  MapType *map;
  MapCache map_cache;
  LinkType *link;
  int nLink;
  int N;
} TriangleSurfaceRec;

int TriangleDegenerate(float *v1,float *n1,float *v2,float *n2,float *v3,float *n3)
{
  float axis1[3];
  float axis2[3];
  float axis3[3];
  float norm[3];
  float dot[3];

  add3f(n1,n2,norm);
  add3f(n3,norm,norm);
  subtract3f(v1,v2,axis1);
  subtract3f(v3,v2,axis2);
  cross_product3f(axis1,axis2,axis3);
  dot[0]=dot_product3f(axis3,n1);
  dot[1]=dot_product3f(axis3,n2);
  dot[2]=dot_product3f(axis3,n3);
  return(!(
           ((dot[0]>0.0F)&&(dot[1]>0.0F)&&(dot[2]>0.0F))||
           ((dot[0]<0.0F)&&(dot[1]<0.0F)&&(dot[2]<0.0F))
           ));
}

static int TriangleEdgeStatus(TriangleSurfaceRec *II,int i1,int i2) 
{
  TriangleSurfaceRec *I=II;
  int l,low,high;
  low = ( i1>i2 ? i2 : i1 );
  high = ( i1>i2 ? i1 : i2 );

  l=I->edgeStatus[low]; 
  while(l) {
	 if(I->link[l].index == high) 
		return(I->link[l].value); /* <0 closed, 0 open, >0 then has index for blocking vector */
	 l=I->link[l].next;
  }
  return(0);
}

static int *TriangleMakeStripVLA(TriangleSurfaceRec *II,float *v,float *vn,int n) 
{
  register TriangleSurfaceRec *I=II;
  int *tFlag,tmp_int;
  int c,a,cc=0;
  int *s,*sc,*strip;
  int state,s01,i0,i1,i2,tc,t0=0;
  int *t;
  int done;
  int dir,dcnt;
  int flag;
  float *v0,*v1,*v2,vt1[3],vt2[3],*tn0,*tn1,*tn2,tn[3],xtn[3];

  strip=VLAlloc(int,I->nTri*4); /* strip VLA is count,vert,vert,...count,vert,vert...zero */
  tFlag=Alloc(int,I->nTri);  
  for(a=0;a<I->nTri;a++) 
	 tFlag[a]=0;
  s=strip;
  done = false;
  while(!done) {
    done = true;
    t=I->tri;
    dir=0;
    for(a=0;a<I->nTri;a++) {
      if(!tFlag[a]) {
        tc=a;
        flag=false;
        dcnt=0;
        while(dcnt<3) {
          i0 = *(t+3*tc+(dir%3));
          i1 = *(t+3*tc+((dir+1)%3));
          state = TriangleEdgeStatus(I,i0,i1);
          if(state) {
            s01 = abs(state);
            t0 = I->edge[s01].tri1;
            if(!tFlag[t0])
              flag=true;
            else if(state<0) {
              t0 = I->edge[s01].tri2;
              if(!tFlag[t0])
                flag=true;
            }
          }
          if(!flag) {
            dir++;
            dcnt++;
          } else {
            c=0;
            sc = s++;
            *(s++)=i0;
            *(s++)=i1;
            while(1) {
              state = TriangleEdgeStatus(I,s[-2],s[-1]);
              if(!state) break;
              s01 = abs(state);
              /*			 printf("a: %i %i %i\n",a,I->edge[s01].tri1,I->edge[s01].tri2);*/
				  
              t0 = I->edge[s01].tri1;
              if(!tFlag[t0])
                i2 = I->edge[s01].vert3;
              else {
                if(state>=0) break;
                t0 = I->edge[s01].tri2;
                i2 = I->edge[s01].vert4;
                /*				printf("second to %i i2 %i  [t0] %i \n",t0,i2,tFlag[t0]);*/
              } 
              if(tFlag[t0]) break;
              *(s++) = i2;
              tFlag[t0]=true;
              c++;
              done=false;
              if((c==1)||(c==2)) { /* make sure vertices follow standard convention */
                
                /* sum normal */
                tn0 = vn+(*(s-3))*3;
                tn1 = vn+(*(s-2))*3;
                tn2 = vn+(*(s-1))*3;
                add3f(tn0,tn1,tn);
                add3f(tn2,tn,tn);

                /* compute right-hand vector */

                v0 = v+(*(s-3))*3;
                v1 = v+(*(s-2))*3;
                v2 = v+(*(s-1))*3;
                subtract3f(v0,v1,vt1);
                subtract3f(v0,v2,vt2);
                cross_product3f(vt1,vt2,xtn);

                if(c&0x1) {
                  /* reorder triangle if necessary */

                  if(dot_product3f(xtn,tn)<0.0) {
                    tmp_int=*(s-3);
                    *(s-3)=*(s-2);
                    *(s-2)=tmp_int;
                  }
                } else {
                  if(dot_product3f(xtn,tn)>0.0) { /* if we're blowing the right hand rule
                                                     on the second triangle, then 
                                                     terminate this strip and try again */
                    tFlag[t0]=false;
                    c--;
                    s--;
                    break;
                  }
                }
              } else { /* continue proofreading handedness... */
                float dp;
                /* sum normal */
                tn0 = vn+(*(s-3))*3;
                tn1 = vn+(*(s-2))*3;
                tn2 = vn+(*(s-1))*3;
                add3f(tn0,tn1,tn);
                add3f(tn2,tn,tn);

                /* compute right-hand vector */

                v0 = v+(*(s-3))*3;
                v1 = v+(*(s-2))*3;
                v2 = v+(*(s-1))*3;
                subtract3f(v0,v1,vt1);
                subtract3f(v0,v2,vt2);
                cross_product3f(vt1,vt2,xtn);
                
                dp = dot_product3f(xtn,tn);
                if(((c&0x1)&&(dp<0.0)) || ((!(c&0x1))&&(dp>0.0))) {
                  /* truncate if right hand rule has been lost */
                  tFlag[t0]=false;
                  c--;
                  s--;
                  break;
                }
              }
            }
            if(!c)
              s=sc;
            else {
              *sc = c; 
              cc+=c;
            }
            /*				if(c>1) printf("strip %i %i\n",c,cc);*/
            dcnt=0;
            tc=t0;			 
            flag=false;
          }
        }
      }
    }

    /* fail-safe check in case of bad connectivity...*/
	 
    for(a=0;a<I->nTri;a++) {
      if(!tFlag[a]) {
        /*		  printf("missed %i %i %i\n",*(I->tri+3*a),*(I->tri+3*a+1), *(I->tri+3*a+2));*/
        *(s++) = 1;
        *(s++) = *(I->tri+3*a);
        *(s++) = *(I->tri+3*a+1);
        *(s++) = *(I->tri+3*a+2);

        /* make sure vertices follow standard convention */
        
        /* sum normal */
        tn0 = vn+(*(s-3))*3;
        tn1 = vn+(*(s-2))*3;
        tn2 = vn+(*(s-1))*3;
        add3f(tn0,tn1,tn);
        add3f(tn2,tn,tn);
        
        /* compute right-hand vector */

        v0 = v+(*(s-3))*3;
        v1 = v+(*(s-2))*3;
        v2 = v+(*(s-1))*3;
        subtract3f(v0,v1,vt1);
        subtract3f(v0,v2,vt2);
        cross_product3f(vt1,vt2,xtn);
        
        /* reorder triangle if necessary */
        
        if(dot_product3f(xtn,tn)<0.0) {
          tmp_int=*(s-3);
          *(s-3)=*(s-2);
          *(s-2)=tmp_int;
        }
      }
    }
    
    *s=0; /* terminate strip list */
  }
  FreeP(tFlag);
  /* shrink strip */
  return(strip);
}

static void TriangleAdjustNormals(TriangleSurfaceRec *II,float *v,float *vn,int n,int final_pass)
{
  register TriangleSurfaceRec *I=II;
  /* points all normals to the average of the intersecting triangles in order to maximum surface smoothness */
  float *tNorm = NULL,*tWght;
  int *vFlag=NULL;
  float *v0,*v1,*v2,*tn,vt1[3],vt2[3],*vn0,*tn0,*tn1,*tn2,*tw;
  int a,*t,i0,i1,i2;
  float tmp[3];
  tNorm=Alloc(float,3*I->nTri);
  tWght=Alloc(float,I->nTri);
  vFlag=Alloc(int,n);
  for(a=0;a<n;a++) {
	 vFlag[a]=0;
  }
  /* first, calculate and store all triangle normals & weights */
  t=I->tri;
  tn=tNorm;
  tw=tWght;
  for(a=0;a<I->nTri;a++) {
	 vFlag[t[0]]=1;
	 vFlag[t[1]]=1;
	 vFlag[t[2]]=1;
	 v0 = v+(*(t++))*3;
	 v1 = v+(*(t++))*3;
	 v2 = v+(*(t++))*3;
	 subtract3f(v1,v0,vt1);
	 subtract3f(v2,v0,vt2);
	 cross_product3f(vt1,vt2,tn);
    *(tw++) = (float)length3f(tn); /* store weight */
	 normalize3f(tn);
	 tn+=3;
  }
  /* clear normals */
  vn0=vn;
  for(a=0;a<n;a++) 
	 if(vFlag[a]) {
		*(vn0++)=0.0;
		*(vn0++)=0.0;
		*(vn0++)=0.0;
	 } else
		vn0+=3;
  
  /* sum */
  tn=tNorm;
  tw=tWght;
  t=I->tri;
  for(a=0;a<I->nTri;a++) {
	 i0 = *(t++);
	 i1 = *(t++);
	 i2 = *(t++);
    scale3f(tn,*tw,tmp);
	 tn0 = vn+i0*3;
	 tn1 = vn+i1*3;
	 tn2 = vn+i2*3;
	 add3f(tmp,tn0,tn0);
	 add3f(tmp,tn1,tn1);
	 add3f(tmp,tn2,tn2);
	 tn+=3;
    tw++;
  }
  /* normalize */
  vn0=vn;
  for(a=0;a<n;a++) {
	 if(vFlag[a])
		normalize3f(vn0);
	 vn0+=3;
  }
  /* now ensure that no normal is allowed to point behind a triangle...*/
  if(final_pass) {
    int repeat = true;
    int max_cyc =5;
    float *va = Alloc(float,3*I->nTri),*va0,*va1,*va2;
    float vt[3];
    while(repeat && max_cyc) {
      repeat = false;
      va0 = va;
      for(a=0;a<n;a++) {
        vFlag[a]=0;
        *(va0++)=0.0F;
        *(va0++)=0.0F;
        *(va0++)=0.0F;
      }
      max_cyc--;
      t=I->tri;
      tn=tNorm;
      for(a=0;a<I->nTri;a++) {
        i0 = *(t++);
        i1 = *(t++);
        i2 = *(t++);
        tn0 = vn+i0*3; /* triangle normal */
        tn1 = vn+i1*3;
        tn2 = vn+i2*3;
        va0 = va+i0*3; /* adjustment */
        va1 = va+i1*3;
        va2 = va+i2*3;
        if(dot_product3f(tn0,tn)<0.0F) {
          remove_component3f(tn0,tn,vt);
          normalize3f(vt);
          add3f(vt,va0,va0);
          vFlag[i0]=true;
          repeat = true;
        }
        if(dot_product3f(tn1,tn)<0.0F) {
          remove_component3f(tn1,tn,vt);
          normalize3f(vt);
          add3f(vt,va1,va1);
          vFlag[i1]=true;
          repeat = true;
        }
        if(dot_product3f(tn2,tn)<0.0F) {
          remove_component3f(tn2,tn,vt);
          normalize3f(vt);
          add3f(vt,va2,va2);
          vFlag[i2]=true;
          repeat = true;
        }
        tn+=3;
      }
      vn0=vn;
      va0=va;
      for(a=0;a<n;a++) {
        if(vFlag[a]) {
          normalize23f(va0,vn0);
        }
        vn0+=3;
        va0+=3;
      }
    }
    FreeP(va);
  }

  FreeP(vFlag);
  FreeP(tWght);
  FreeP(tNorm);
}



static int TriangleActivateEdges(TriangleSurfaceRec *II,int low)
{
  register TriangleSurfaceRec *I=II;
  int l;
  l=I->edgeStatus[low]; 
  while(l) {
    if(I->link[l].value>0) {
      VLACheck(I->activeEdge,int,I->nActive*2+1);
      I->activeEdge[I->nActive*2] = low;
      I->activeEdge[I->nActive*2+1] = I->link[l].index;
      I->nActive++;
    }
	 l=I->link[l].next;
  }
  return(0);
}

static void TriangleDeleteEdge(TriangleSurfaceRec *II,int i1,int i2)
{
  register TriangleSurfaceRec *I=II;
  int l,low,high;
  int prev = 0;
  low = ( i1>i2 ? i2 : i1 );
  high = ( i1>i2 ? i1 : i2 );

  /*  printf("set: %i %i %i\n",i1,i2,value);*/
  l=I->edgeStatus[low]; 
  while(l) {
	 if(I->link[l].index == high) {
       if(prev) {
         I->link[prev].next = I->link[l].next; /* leaks, but that's no big deal */
         return;
       } else {
         I->edgeStatus[low] = I->link[l].next; /* leaks, but that's no big deal */
       }
	 }
     prev = l;
     l=I->link[l].next;
  }
  return;
}

static void TriangleEdgeSetStatus(TriangleSurfaceRec *II,int i1,int i2,int value) 
{
  register TriangleSurfaceRec *I=II;
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

static void TriangleMove(TriangleSurfaceRec *II,int from,int to)
{
  register TriangleSurfaceRec *I=II;
  int i0,i1,i2,s01,s02,s12;

  i0=I->tri[from*3];
  i1=I->tri[from*3+1];
  i2=I->tri[from*3+2];

  s01 = TriangleEdgeStatus(I,i0,i1);
  s02 = TriangleEdgeStatus(I,i0,i2);
  s12 = TriangleEdgeStatus(I,i1,i2);

#define TMoveMacro(x) {\
  if(x>0) { \
	 if(I->edge[x].tri1==from) \
		I->edge[x].tri1=to; \
	 else if(I->edge[x].tri2==from) \
		I->edge[x].tri2=to; \
  } else if(x<0) { \
	 x=-x; \
	 if(I->edge[x].tri1==from) \
		I->edge[x].tri1=to; \
	 else if(I->edge[x].tri2==from) \
		I->edge[x].tri2=to; \
  }}\

  TMoveMacro(s01)
  TMoveMacro(s02)
  TMoveMacro(s12)
  
  I->tri[to*3]=i0;
  I->tri[to*3+1]=i1;
  I->tri[to*3+2]=i2;
}

#define max_edge_len 2.5
static void TriangleRectify(TriangleSurfaceRec *I, int t, float *v,float *vn)
{
  int i0 = I->tri[t*3+0];
  int i1 = I->tri[t*3+1];
  int i2 = I->tri[t*3+2],it;
  float *n0 = vn+3*i0, *n1 = vn+3*i1, *n2 = vn+3*i2;							 
  float *v0 = v +3*i0, *v1 = v +3*i1, *v2 = v +3*i2;
  float tNorm[3], vt1[3],vt2[3],vt[3];

  add3f(n0,n1,tNorm);
  add3f(n2,tNorm,tNorm);

  subtract3f(v1,v0,vt1);
  subtract3f(v2,v0,vt2);
  cross_product3f(vt1,vt2,vt);
  if(dot_product3f(vt,tNorm)<0.0) { /* if wrong, then interchange */
    it=i1; i1=i2; i2=it;
    I->tri[t*3+1] = i1;
    I->tri[t*3+2] = i2;
  }
}

static void AddActive(TriangleSurfaceRec *II,int i1,int i2) {
  int t;

  register TriangleSurfaceRec *I=II;
  if(i1>i2) {
    t=i1;
    i1=i2;
    i2=t;
  }
  VLACheck(I->activeEdge,int,I->nActive*2+1);
  I->activeEdge[I->nActive*2] = i1;
  I->activeEdge[I->nActive*2+1] = i2;
  I->nActive++;
  if(I->vertActive[i1]<0) I->vertActive[i1]=0;
  I->vertActive[i1]++;
  if(I->vertActive[i2]<0) I->vertActive[i2]=0;
  I->vertActive[i2]++;
}

static void TriangleAdd(TriangleSurfaceRec *II,int i0,int i1,int i2,float *tNorm,float *v,float *vn)
{
  register TriangleSurfaceRec *I=II;
  int s01,s12,s02,it;
  float *v0,*v1,*v2,*n0,*n1,*n2,*ft;
  float vt[3],vt1[3],vt2[3],vt3[3];
  int e1,e2,e3,h,j,k,l;
  MapType *map = I->map;
  MapCache *cache = &I->map_cache;

  v0 = v+3*i0;
  v1 = v+3*i1;
  v2 = v+3*i2;
  n0 = vn+3*i0;
  n1 = vn+3*i1;
  n2 = vn+3*i2;

  /* mark this quadrant as visited so that further activity in this
     quadrant won't prolong the rendering process indefinitely */
  
  
  MapLocus(map,v0,&h,&k,&l);
  e1=*(MapEStart(map,h,k,l));
  
  if(e1) {
    j=map->EList[e1];
    MapCache(cache,j);
  }
  
  MapLocus(map,v1,&h,&k,&l);
  e2=*(MapEStart(map,h,k,l));
  if(e2&&(e1!=e2)) {
    j=map->EList[e2];
    MapCache(cache,j);
  }
  
  MapLocus(map,v2,&h,&k,&l);
  e3=*(MapEStart(map,h,k,l));
  if(e3&&(e3!=e2)&&(e3!=e1)) {
    j=map->EList[e3];
    MapCache(cache,j);
  }

  /* make sure the triangle obeys the right hand rule */

  subtract3f(v1,v0,vt1);
  subtract3f(v2,v0,vt2);
  cross_product3f(vt1,vt2,vt);
  if(dot_product3f(vt,tNorm)<0.0) { /* if wrong, then interchange */
	 it=i1; i1=i2; i2=it;
	 ft=v1; v1=v2; v2=ft;
	 ft=n1; n1=n2; n2=ft;
	 }
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

  s01 = TriangleEdgeStatus(I,i0,i1);
  s02 = TriangleEdgeStatus(I,i0,i2);
  s12 = TriangleEdgeStatus(I,i1,i2);
  
  /* create a new triangle */
  VLACheck(I->tri,int,(I->nTri*3)+2);
  I->tri[I->nTri*3] = i0;
  I->tri[I->nTri*3+1] = i1;
  I->tri[I->nTri*3+2] = i2;
  /*  	printf("creating %i %i %i\n",i0,i1,i2);*/
  if(s01) {
	 if(s01>0) {
		I->edge[s01].vert4 = i2;
		I->edge[s01].tri2 = I->nTri;
		TriangleEdgeSetStatus(I,i0,i1,-s01);
		I->vertActive[i0]--; /* deactivate when all active edges are closed */
		I->vertActive[i1]--;
	 } /*else {
      ErrFatal(I->G,"TriangleAdd","Invalid triangle - s01 negative");
      }*/
  } else {
	 VLACheck(I->edge,EdgeRec,I->nEdge);
	 I->edge[I->nEdge].vert3 = i2;
	 I->edge[I->nEdge].tri1 = I->nTri;
	 VLACheck(I->vNormal,float,(I->nEdge*3)+2);
	 copy3f(tNorm,I->vNormal+I->nEdge*3);
	 TriangleEdgeSetStatus(I,i0,i1,I->nEdge);
	 I->nEdge++;
    AddActive(I,i0,i1);
  }
  if(s02) {
	 if(s02>0) {
		I->edge[s02].vert4 = i1;
		I->edge[s02].tri2 = I->nTri;
		TriangleEdgeSetStatus(I,i0,i2,-s02);
		I->vertActive[i0]--; /* deactivate when all active edges are closed */
		I->vertActive[i2]--;
	 } /*else {
      ErrFatal(I->G,"TriangleAdd","Invalid triangle - s02 negative");
      }*/
  } else {
	 VLACheck(I->edge,EdgeRec,I->nEdge);
	 I->edge[I->nEdge].vert3 = i1;
	 I->edge[I->nEdge].tri1 = I->nTri;
	 VLACheck(I->vNormal,float,(I->nEdge*3)+2);
	 copy3f(tNorm,I->vNormal+I->nEdge*3);
	 TriangleEdgeSetStatus(I,i0,i2,I->nEdge);
	 I->nEdge++;
    AddActive(I,i0,i2);
  }
  if(s12) {
	 if(s12>0) {
		I->edge[s12].vert4 = i0;
		I->edge[s12].tri2 = I->nTri;
		TriangleEdgeSetStatus(I,i1,i2,-s12); 
		I->vertActive[i1]--; /* deactivate when all active edges are closed */
		I->vertActive[i2]--;
	 } /*else {
      ErrFatal(I->G,"TriangleAdd","Invalid triangle - s12 negative");
      }*/
  } else {
	 VLACheck(I->edge,EdgeRec,I->nEdge);
	 I->edge[I->nEdge].vert3 = i0;
	 I->edge[I->nEdge].tri1 = I->nTri;
	 VLACheck(I->vNormal,float,(I->nEdge*3)+2);
	 copy3f(tNorm,I->vNormal+I->nEdge*3);
	 TriangleEdgeSetStatus(I,i1,i2,I->nEdge);
	 I->nEdge++;
    AddActive(I,i1,i2);
  }
  I->nTri++;

}

static void TriangleBuildObvious(TriangleSurfaceRec *II,int i1,int i2,float *v,float *vn,int n)
{
  /* this routine builds obvious, easy triagles where the closest point 
  * to the edge is always tried */

  register TriangleSurfaceRec *I=II;
  MapType *map;
  float *v1,*v2,*v0,*v4,vt[3],vt1[3],vt2[3],vt3[3],vt4[3],*n0,*n1,*n2,tNorm[3];
  int i0,s01=0,s02=0,s12,i,j,h,k,l,i4;
  float dif,minDist,d1,d2,dp;
  int flag;
  int used = -1;
  float maxDot,dot,dot1,dot2;
  const float _plus = R_SMALL4, _0=0.0F;
  const float _5 = 0.5F;
  /*  PRINTFD(I->G,FB_Triangle)
      " TriangleBuildObvious-Debug: entered: i1=%d i2=%d n=%d\n",i1,i2,n
      ENDFD;*/

  map=I->map;
  s12 = TriangleEdgeStatus(I,i1,i2);

  /*  PRINTFD(I->G,FB_Triangle)
    " TriangleBuildObvious-Debug: edge status=%d\n",s12
    ENDFD;
  */

  if(s12>0) used = I->edge[s12].vert3;
  if(s12>=0) {
    minDist = MAXFLOAT;
    maxDot = _plus;
    i0=-1;
    v1=v+i1*3; v2=v+i2*3;
    n1 = vn+3*i1; n2 = vn+3*i2;							 
    MapLocus(map,v1,&h,&k,&l);
	 i=*(MapEStart(map,h,k,l));
    if(i) {
      j=map->EList[i++];
      while(j>=0) {
        if((j!=i1)&&(j!=i2)&&(j!=used)) 
          {
            v0 = v+3*j;
            d1 = (float)diff3f(v0,v1);
            d2 = (float)diff3f(v0,v2);
            dif= ( d2 > d1 ? d2 : d1 );
				if(dif<minDist)
				  {
                n0 = vn + 3*j;
                dot1 = dot_product3f(n0,n1);
                dot2 = dot_product3f(n0,n2);
                dot = dot1 + dot2;
                if((dif/minDist)<_5) {
                  minDist = dif;
                  maxDot = dot;
                  i0=j; 
                } else if((dot>_0)&&(dot1>_0)&&(dot2>_0)) {
                  if((i0<0)||(dot>maxDot)) {
                    minDist = dif;
                    maxDot = dot;
                    i0=j; 
                  } else if( _5*(dif/minDist) < (dot/maxDot)) {
                    minDist = dif;
                    maxDot = dot;
                    i0=j; 
                  }
                }
				  }
			 }
		  j=map->EList[i++];
		}
      /*      PRINTFD(I->G,FB_Triangle)
        " TriangleBuildObvious-Debug: i0=%d\n",i0
        ENDFD;*/

		if(i0>=0) {
		  s01 = TriangleEdgeStatus(I,i0,i1); s02 = TriangleEdgeStatus(I,i0,i2);
		  if(I->vertActive[i0]>0) {
			 if(!((s01>0)||(s02>0)))
				i0=-1; /* don't allow non-adjacent joins to active vertices */
		  }
		}

      /* PRINTFD(I->G,FB_Triangle)
        " TriangleBuildObvious-Debug: i0=%d s01=%d s02=%d\n",i0,s01,s02
        ENDFD;*/

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
            /*            if(Feedback(I->G,FB_Triangle,FB_Debugging)) {
                          dump3f(n0,"n0");
                          dump3f(n1,"n1");
                          dump3f(n2,"n2");
                          dump3f(vt2,"vt2");
                          PRINTF " n0.vt2 = %8.3f\n",dot_product3f(n0,vt2) ENDF(I->G);
                          PRINTF " n1.vt2 = %8.3f\n",dot_product3f(n1,vt2) ENDF(I->G);
                          PRINTF " n2.vt2 = %8.3f\n",dot_product3f(n2,vt2) ENDF(I->G);
                          PRINTF " n0.vt2<0.1 %d\n",dot_product3f(n0,vt2)<0.1 ENDF(I->G);
                          PRINTF " n1.vt2<0.1 %d\n",dot_product3f(n1,vt2)<0.1 ENDF(I->G);
                          PRINTF " n2.vt2<0.1 %d\n",dot_product3f(n2,vt2)<0.1 ENDF(I->G);
                          fflush(stdout);
                          }*/
				if(((dot_product3f(n0,vt2))<0.1)||
               ((dot_product3f(n1,vt2))<0.1)||
               ((dot_product3f(n2,vt2))<0.1)) flag = false; 
            /* modified 010916 to effect workaround of apparent bug in GCC's optimizer */
			 } 
          /*          PRINTFD(I->G,FB_Triangle)
                      " TriangleBuildObvious-Debug: past normal sums, flag= %d\n",flag
            ENDFD;*/
			 if(flag) { /* does the sum of the normals point in the same direction as the triangle? */
				subtract3f(v1,v0,vt3);
				subtract3f(v2,v0,vt4);
				cross_product3f(vt3,vt4,tNorm); 
				normalize3f(tNorm); 							 
				dp = dot_product3f(vt2,tNorm);
				if(dp<0) scale3f(tNorm,-1.0F,tNorm);
				if(fabs(dp)<0.1) flag = false;
			 } 
          /*          PRINTFD(I->G,FB_Triangle)
                      " TriangleBuildObvious-Debug: past tNorm, flag= %d\n",flag
            ENDFD;*/
			 if(flag) {
				if(s12>0) if(dot_product3f(I->vNormal+s12*3,tNorm)<0.1) flag=false; 
				if(s01>0) if(dot_product3f(I->vNormal+s01*3,tNorm)<0.1) flag=false; 
				if(s02>0) if(dot_product3f(I->vNormal+s02*3,tNorm)<0.1) flag=false; 
			 } 
          /*          PRINTFD(I->G,FB_Triangle)
                      " TriangleBuildObvious-Debug: past compare tNorm, flag= %d\n",flag
            ENDFD;*/
			 if(flag) { /* are all the Blocking vectors pointing outward, and are the triangle normals consistent? */
				if(s12>0) {
				  i4 = I->edge[s12].vert3;
				  v4=v+i4*3;
				  subtract3f(v0,v1,vt1);
				  subtract3f(v4,v1,vt2);
				  subtract3f(v1,v2,vt);
				  normalize3f(vt);
				  remove_component3f(vt1,vt,vt3);
				  remove_component3f(vt2,vt,vt4);
				  normalize3f(vt3);
				  normalize3f(vt4);
				  if((dot_product3f(vt3,vt4))>0.0) flag=false;
				}			 
				if(s01>0) {
				  i4 = I->edge[s01].vert3;
				  v4=v+i4*3;
				  subtract3f(v2,v0,vt1);
				  subtract3f(v4,v0,vt2);
				  subtract3f(v0,v1,vt);
				  normalize3f(vt);
				  remove_component3f(vt1,vt,vt3);
				  remove_component3f(vt2,vt,vt4);
				  normalize3f(vt3);
				  normalize3f(vt4);
				  if((dot_product3f(vt3,vt4))>0.0) flag=false;
				}
				if(s02>0) {
				  i4 = I->edge[s02].vert3;
				  v4=v+i4*3;
				  subtract3f(v1,v0,vt1);
				  subtract3f(v4,v0,vt2);
				  subtract3f(v0,v2,vt);
				  normalize3f(vt);
				  remove_component3f(vt1,vt,vt3);
				  remove_component3f(vt2,vt,vt4);
				  normalize3f(vt3);
				  normalize3f(vt4);
				  if((dot_product3f(vt3,vt4))>0.0) flag=false;
				}
			 } 
          /*          PRINTFD(I->G,FB_Triangle)
            " TriangleBuildObvious-Debug: past blocking, flag= %d\n",flag
            ENDFD;*/
		  }
		  if(flag) TriangleAdd(I,i0,i1,i2,tNorm,v,vn);
		}
	 }
  }
}

static void TriangleBuildSecondPass(TriangleSurfaceRec *II,int i1,int i2,float *v,float *vn,int n)
{

  /* in this version, the closest active point is tried.  Closed points
	  are skipped. */

  register TriangleSurfaceRec *I=II;
  MapType *map;
  float *v1,*v2,*v0,*v4,vt[3],vt1[3],vt2[3],vt3[3],vt4[3],*n0,*n1,*n2,tNorm[3];
  int i0,s01=0,s02=0,s12,i,j,h,k,l,i4;
  float dif,minDist,d1,d2,dp;
  int flag;
  int used = -1;
  float dot,dot1,dot2,maxDot;
  const float _plus = R_SMALL4, _0=0.0F;
  const float _5 = 0.5F;

  map=I->map;
  s12 = TriangleEdgeStatus(I,i1,i2);
  if(s12>0) used = I->edge[s12].vert3;
  if(s12>=0) {
    minDist = MAXFLOAT;
    maxDot = _plus;
    i0=-1;
    v1=v+i1*3; v2=v+i2*3;
    n1 = vn+3*i1; n2 = vn+3*i2;							 
    MapLocus(map,v1,&h,&k,&l);
	 i=*(MapEStart(map,h,k,l));
    if(i) {
      j=map->EList[i++];
      while(j>=0) {
        if((j!=i1)&&(j!=i2)&&(j!=used)&&(I->vertActive[j])) 
			 /* eliminate closed vertices from consideration - where vertactive is 0 */
          {
            v0 = v+3*j;
            n0 = vn+3*j;
            d1 = (float)diff3f(v0,v1);
            d2 = (float)diff3f(v0,v2);
            dif= ( d2 > d1 ? d2 : d1 );
				if(dif<minDist)
				  {
                dot1 = dot_product3f(n0,n1);
                dot2 = dot_product3f(n0,n2);
                dot = dot1 + dot2;
                if((dif/minDist)<_5) {
                  minDist = dif;
                  maxDot = dot;
                  i0=j; 
                } else if((dot>_0)&&(dot1>_0)&&(dot2>_0)) {
                  if((i0<0)||(dot>maxDot)) {
                    minDist = dif;
                    maxDot = dot;
                    i0=j; 
                  } else if( _5*(dif/minDist) < (dot/maxDot)) {
                    maxDot = dot;
                    minDist = dif;
                    i0=j; 
                  }
                }
				  }
			 }
		  j=map->EList[i++];
		}
		if(i0>=0) {
		  s01 = TriangleEdgeStatus(I,i0,i1); s02 = TriangleEdgeStatus(I,i0,i2);
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
				if(((dot_product3f(n0,vt2))<0.1)||
               ((dot_product3f(n1,vt2))<0.1)||
               ((dot_product3f(n2,vt2))<0.1)) flag=false;
            /* modified 010916 to effect workaround of apparent bug in GCC's optimizer */
			 } /*printf("pass normal sums %i\n",flag);*/
			 if(flag) { /* does the sum of the normals point in the same direction as the triangle? */
				subtract3f(v1,v0,vt3);
				subtract3f(v2,v0,vt4);
				cross_product3f(vt3,vt4,tNorm); 
				normalize3f(tNorm); 							 
				dp = dot_product3f(vt2,tNorm);
				if(dp<0) scale3f(tNorm,-1.0F,tNorm);
				if(fabs(dp)<0.1) flag = false;
			 } /*printf("pass tNorm  %i\n",flag);*/

			 if(flag) {
				if(s12>0) if(dot_product3f(I->vNormal+s12*3,tNorm)<0.1) flag=false; 
				if(s01>0) if(dot_product3f(I->vNormal+s01*3,tNorm)<0.1) flag=false; 
				if(s02>0) if(dot_product3f(I->vNormal+s02*3,tNorm)<0.1) flag=false; 
			 } /*printf("pass compare tNorm %i\n",flag);*/
			 if(flag) { /* are all the Blocking vectors pointing outward, and are the triangle normals consistent? */
				if(s12>0) {
				  i4 = I->edge[s12].vert3;
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
				  i4 = I->edge[s01].vert3;
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
				  i4 = I->edge[s02].vert3;
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
		  if(flag) TriangleAdd(I,i0,i1,i2,tNorm,v,vn);
		}
	 }
  }
}

static void TriangleBuildSecondSecondPass(TriangleSurfaceRec *II,int i1,int i2,float *v,float *vn,int n,float cutoff)
{

  /* in this version, the closest active point is tried.  Closed points
	  are skipped. */

  register TriangleSurfaceRec *I=II;
  MapType *map;
  float *v1,*v2,*v0,*v4,vt[3],vt1[3],vt2[3],vt3[3],vt4[3],*n0,*n1,*n2,tNorm[3];
  int i0,s01=0,s02=0,s12,i,j,h,k,l,i4;
  float dif,minDist,d1,d2,dp;
  int flag;
  int used = -1;
  float dot;
  const float _5 = 0.5;

  map=I->map;
  s12 = TriangleEdgeStatus(I,i1,i2);
  if(s12>0) used = I->edge[s12].vert3;
  if(s12>=0) {
    minDist = MAXFLOAT;
    i0=-1;
    v1=v+i1*3; v2=v+i2*3;
    n1 = vn+3*i1; n2 = vn+3*i2;							 
    MapLocus(map,v1,&h,&k,&l);
	 i=*(MapEStart(map,h,k,l));
    if(i) {
      j=map->EList[i++];
      while(j>=0) {
        if((j!=i1)&&(j!=i2)&&(j!=used)&&(I->vertActive[j])) 
			 /* eliminate closed vertices from consideration - where vertactive is 0 */
          {
            v0 = v+3*j;
            n0 = vn+3*j;
            d1 = (float)diff3f(v0,v1);
            d2 = (float)diff3f(v0,v2);
            dif= ( d2 > d1 ? d2 : d1 );
				if(dif<minDist)
				  {
                dot = dot_product3f(n0,n1) + dot_product3f(n0,n2);
                if((dot>cutoff)||((dif/minDist)<_5)) {
                  minDist = dif;
                  i0=j; 
                }
				  }
			 }
		  j=map->EList[i++];
		}
		if(i0>=0) {
		  s01 = TriangleEdgeStatus(I,i0,i1); s02 = TriangleEdgeStatus(I,i0,i2);
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
				if(((dot_product3f(n0,vt2))<0.1)||
               ((dot_product3f(n1,vt2))<0.1)||
               ((dot_product3f(n2,vt2))<0.1)) flag=false;
            /* modified 010916 to effect workaround of apparent bug in GCC's optimizer */
			 } /*printf("pass normal sums %i\n",flag);*/
			 if(flag) { /* does the sum of the normals point in the same direction as the triangle? */
				subtract3f(v1,v0,vt3);
				subtract3f(v2,v0,vt4);
				cross_product3f(vt3,vt4,tNorm); 
				normalize3f(tNorm); 							 
				dp = dot_product3f(vt2,tNorm);
				if(dp<0) scale3f(tNorm,-1.0F,tNorm);
				if(fabs(dp)<0.1) flag = false;
			 } /*printf("pass tNorm  %i\n",flag);*/
			 if(flag) {
				if(s12>0) if(dot_product3f(I->vNormal+s12*3,tNorm)<0.1) flag=false; 
				if(s01>0) if(dot_product3f(I->vNormal+s01*3,tNorm)<0.1) flag=false; 
				if(s02>0) if(dot_product3f(I->vNormal+s02*3,tNorm)<0.1) flag=false; 
			 } /*printf("pass compare tNorm %i\n",flag);*/
			 if(flag) { /* are all the Blocking vectors pointing outward, and are the triangle normals consistent? */
				if(s12>0) {
				  i4 = I->edge[s12].vert3;
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
				  i4 = I->edge[s01].vert3;
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
				  i4 = I->edge[s02].vert3;
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
		  if(flag) TriangleAdd(I,i0,i1,i2,tNorm,v,vn);
		}
	 }
  }
}



static void TriangleBuildSingle(TriangleSurfaceRec *II,int i1,int i2,float *v,float *vn,int n)
{

  register TriangleSurfaceRec *I=II;
  MapType *map;
  float *v1,*v2,*v0,*v4,vt[3],vt1[3],vt2[3],vt3[3],vt4[3],*n0,*n1,*n2,tNorm[3];
  int i0,s01=0,s02=0,s12,i,j,h,k,l,i4;
  float dif,minDist,d1,d2,dp;
  int flag;
  int used = -1;

  map=I->map;
  s12 = TriangleEdgeStatus(I,i1,i2);
  if(s12>0) used = I->edge[s12].vert3;
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
            d1 = (float)diff3f(v0,v1);
            d2 = (float)diff3f(v0,v2);
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
		  s01 = TriangleEdgeStatus(I,i0,i1); s02 = TriangleEdgeStatus(I,i0,i2);
		  if(I->vertActive[i0]) {
			 if((s01>=0)&&(s02>=0)) flag=true;
			 if(flag) { /* are all normals pointing in generally the same direction? */
				n0 = vn+3*i0; n1 = vn+3*i1; n2 = vn+3*i2;							 
				add3f(n0,n1,vt1);
				add3f(n2,vt1,vt2);
				normalize3f(vt2);
				if(((dot_product3f(n0,vt2))<0.1)|| 
               ((dot_product3f(n1,vt2))<0.1)|| 
               ((dot_product3f(n2,vt2))<0.1)) flag=false;
            /* modified 010916 to effect workaround of apparent bug in GCC's optimizer */
			 } /*printf("pass normal sums %i\n",flag);*/
			 if(flag) { /* does the sum of the normals point in the same direction as the triangle? */
				subtract3f(v1,v0,vt3);
				subtract3f(v2,v0,vt4);
				cross_product3f(vt3,vt4,tNorm); 
				normalize3f(tNorm); 							 
				dp = dot_product3f(vt2,tNorm);
				if(dp<0) scale3f(tNorm,-1.0F,tNorm);
				if(fabs(dp)<0.1) flag = false;
			 } /*printf("pass tNorm  %i\n",flag);*/
			 if(flag) {
				if(s12>0) if(dot_product3f(I->vNormal+s12*3,tNorm)<0.1) flag=false; 
				if(s01>0) if(dot_product3f(I->vNormal+s01*3,tNorm)<0.1) flag=false; 
				if(s02>0) if(dot_product3f(I->vNormal+s02*3,tNorm)<0.1) flag=false; 
			 } /*printf("pass compare tNorm %i\n",flag);*/
			 if(flag) { /* are all the Blocking vectors pointing outward, and are the triangle normals consistent? */
				if(s12>0) {
				  i4 = I->edge[s12].vert3;
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
				  i4 = I->edge[s01].vert3;
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
				  i4 = I->edge[s02].vert3;
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
		  if(flag) TriangleAdd(I,i0,i1,i2,tNorm,v,vn);
		}
	 }
  }
}


static void TriangleBuildThirdPass(TriangleSurfaceRec *II,int i1,int i2,float *v,float *vn,int n)
{
  /* This routine fills in triangles surrounded by three active edges */

  register TriangleSurfaceRec *I=II;
  MapType *map;
  float *v1,*v2,*v0,vt1[3],vt2[3],vt3[3],vt4[3],*n0,*n1,*n2,tNorm[3];
  int i0,s01,s02,s12,i,j,h,k,l;
  float dif,minDist,d1,d2,dp;
  int used = -1;

  map=I->map;
  s12 = TriangleEdgeStatus(I,i1,i2);
  if(s12>0) used = I->edge[s12].vert3;
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
            d1 = (float)diff3f(v0,v1);
            d2 = (float)diff3f(v0,v2);
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
		  s01 = TriangleEdgeStatus(I,i0,i1); s02 = TriangleEdgeStatus(I,i0,i2);
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
			 if(dp<0) scale3f(tNorm,-1.0F,tNorm);
			 TriangleAdd(I,i0,i1,i2,tNorm,v,vn);
		  }
		}
	 }
  }
}


static void TriangleBuildLast(TriangleSurfaceRec *II,int i1,int i2,float *v,float *vn,int n)
{
  /* this routine is a hack to fill in the odd-ball situations */

  register TriangleSurfaceRec *I=II;
  MapType *map;
  float *v1,*v2,*v0,vt1[3],vt2[3],vt3[3],vt4[3],*n0,*n1,*n2,tNorm[3];
  int i0,s01,s02,s12,i,j,h,k,l;
  float dif,minDist,d1,d2,dp;
  int used = -1;
  int both_active;
  map=I->map;
  s12 = TriangleEdgeStatus(I,i1,i2);
  if(s12>0) used = I->edge[s12].vert3;
  if(s12>=0) {
    minDist = MAXFLOAT;
    i0=-1;
	 both_active = (I->vertActive[i1]&&I->vertActive[i2]);
    v1=v+i1*3; v2=v+i2*3;
    MapLocus(map,v1,&h,&k,&l);
	 i=*(MapEStart(map,h,k,l));
    if(i) {
      j=map->EList[i++];
      while(j>=0) {
        if((j!=i1)&&(j!=i2)&&(j!=used)&&(I->vertActive[j]>0))
          {
            v0 = v+3*j;
            d1 = (float)diff3f(v0,v1);
            d2 = (float)diff3f(v0,v2);
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
		  s01 = TriangleEdgeStatus(I,i0,i1); s02 = TriangleEdgeStatus(I,i0,i2);
		  /* if all three edges are active */
		  if(((s12>0)&&(((s01>0)&&(s02>=0))||((s01>=0)&&(s02>0))))||
			  ((s01>0)&&(s02>0))) { 
			 n0 = vn+3*i0; n1 = vn+3*i1; n2 = vn+3*i2;							 
			 add3f(n0,n1,vt1);
			 add3f(n2,vt1,vt2);
			 subtract3f(v1,v0,vt3);
			 subtract3f(v2,v0,vt4);
			 cross_product3f(vt3,vt4,tNorm); 
			 normalize3f(tNorm); 							 
			 dp = dot_product3f(vt2,tNorm);
			 if(dp<0) scale3f(tNorm,-1.0F,tNorm);
			 TriangleAdd(I,i0,i1,i2,tNorm,v,vn);
		  }
		}
	 }
  }
}



static void FollowActives(TriangleSurfaceRec *II,float *v,float *vn,int n,int mode)
{
  register TriangleSurfaceRec *I=II;
  int i1,i2;
  
  PRINTFD(I->G,FB_Triangle)
    " TriangleFollowActives-Debug: entered: n=%6d     mode=%d\n TriangleFollowActives-Debug:       nTri=%6d nActive=%6d\n",
    n,mode,I->nTri,I->nActive
    ENDFD;

  OrthoBusyFast(I->G,(I->N*3)+I->nTri,I->N*5); /* 3/5 to 4/5 */

  while(I->nActive) {
    I->nActive--;
    i1 = I->activeEdge[I->nActive*2];
    i2 = I->activeEdge[I->nActive*2+1];
	 switch(mode) {
	 case 0:
		TriangleBuildObvious(I,i1,i2,v,vn,n);
		break;
	 case 1:
		TriangleBuildSecondPass(I,i1,i2,v,vn,n);
		break;
	 case 2:
		TriangleBuildSecondSecondPass(I,i1,i2,v,vn,n,0.0F);
		break;
	 case 4:
		TriangleBuildThirdPass(I,i1,i2,v,vn,n);
		break;
	 case 5:
		TriangleBuildLast(I,i1,i2,v,vn,n);
		break;
	 }
  }

  PRINTFD(I->G,FB_Triangle)
    " TriangleFollowActives-Debug: exiting: nTri=%6d nActive=%6d\n",
    I->nTri,I->nActive
    ENDFD;

}

static void TriangleFill(TriangleSurfaceRec *II,float *v,float *vn,int n,int first_time)
{
  register TriangleSurfaceRec *I=II;
  int lastTri,lastTri2,lastTri3;
  int a,i,j,h,k,l;
  float dif,minDist,*v0,*n0,*n1;
  int i1,i2=0;
  int n_pass =0;
  int first_vert=0,first_vert_used=0;

  MapType *map;
  MapCache *cache;

  PRINTFD(I->G,FB_Triangle)
    " TriangleFill-Debug: entered: n=%d\n",n
    ENDFD;
  
  map=I->map;
  cache = &I->map_cache;

  lastTri3=-1;
  while(lastTri3!=I->nTri) {
	 lastTri3=I->nTri;
    n_pass++;
    if(n_pass>(int)SettingGet(I->G,cSetting_triangle_max_passes))
      break;

	 I->nActive=0;
    while((!I->nActive)&&(I->nTri==lastTri3))
      {
        i1=-1;
        minDist = MAXFLOAT;
        
        for(a=0;a<n;a++) 
          if(!I->edgeStatus[a])
            {
              v0=v+a*3;
              n0 = vn+3*a;

              MapLocus(map,v0,&h,&k,&l);
              i=*(MapEStart(map,h,k,l));
              if(i) {
                j=map->EList[i++];
                first_vert = j;
                while(j>=0) {
                  if(j!=a) 
                    {
                      dif=(float)diff3f(v+3*j,v0);
                      if(dif<minDist)
                        if(I->vertActive[a]==-1)
                          if(TriangleEdgeStatus(I,a,j)>=0) /* can we put a triangle here? */
                            {
                              n1 = vn + 3*j;
                              if(dot_product3f(n0,n1)>0.5) { /* start with vertices pointing the same way */
                                minDist = dif;
                                i1=a;
                                i2=j;
                                first_vert_used = first_vert;
                              }
                            }
                    }
                  j=map->EList[i++];
                }
              }
            }
        
        if(i1>=0) {

          if(!MapCached(cache,first_vert_used)) {
            MapCache(cache,first_vert_used);
            if(first_time) {
              n_pass=n_pass/2; 
              /* if we've entered a new map quadrant then half the effective number of passes */

              /* this is a very non-obvious way of making sure that we
                 don't prematurely terminate when surfacing
                 discontinuous surfaces that will always require many passes */
            }
          }

          if(I->vertActive[i1]<0) I->vertActive[i1]--;
          VLACheck(I->activeEdge,int,I->nActive*2+1);
          I->activeEdge[I->nActive*2] = i1;
          I->activeEdge[I->nActive*2+1] = i2;
          I->nActive=1;
          lastTri=I->nTri;
          FollowActives(I,v,vn,n,0);
          while(lastTri!=I->nTri) {
            lastTri=I->nTri;
            for(a=0;a<n;a++) 
              if(I->vertActive[a])
                TriangleActivateEdges(I,a);
            FollowActives(I,v,vn,n,0);
          }
        } else break;
      }


    PRINTFD(I->G,FB_Triangle)
      " TriangleFill-Debug: Follow actives 1 nTri=%d\n",I->nTri
      ENDFD;
	 lastTri=I->nTri-1;
	 while(lastTri!=I->nTri) {
		lastTri=I->nTri;
		for(a=0;a<n;a++) 
		  if(I->vertActive[a])
			 TriangleActivateEdges(I,a);
		FollowActives(I,v,vn,n,1);
	 }	 

    lastTri2=I->nTri-1;
    while(lastTri2!=I->nTri) {
      lastTri2=I->nTri;
      for(a=0;a<n;a++) 
        if(I->vertActive[a])
          {
            TriangleActivateEdges(I,a);
            if(I->nActive) {
              PRINTFD(I->G,FB_Triangle)
                " TriangleFill-Debug: build single:     nTri=%d nActive=%d\n",I->nTri,I->nActive
                ENDFD;
              I->nActive--;
              i1 = I->activeEdge[I->nActive*2];
              i2 = I->activeEdge[I->nActive*2+1];
              TriangleBuildSingle(I,i1,i2,v,vn,n);
              PRINTFD(I->G,FB_Triangle)
                " TriangleFill-Debug: follow actives 1: nTri=%d nActive=%d\n",I->nTri,I->nActive
                ENDFD;
              FollowActives(I,v,vn,n,1);
            }
          }
    }

    PRINTFD(I->G,FB_Triangle)
      " TriangleFill-Debug: Follow actives 1 nTri=%d\n",I->nTri
      ENDFD;
	 lastTri=I->nTri-1;
	 while(lastTri!=I->nTri) {
		lastTri=I->nTri;
		for(a=0;a<n;a++) 
		  if(I->vertActive[a])
			 TriangleActivateEdges(I,a);
		FollowActives(I,v,vn,n,2);
	 }	 

    lastTri2=I->nTri-1;
    while(lastTri2!=I->nTri) {
      lastTri2=I->nTri;
      for(a=0;a<n;a++) 
        if(I->vertActive[a])
          {
            TriangleActivateEdges(I,a);
            if(I->nActive) {
              PRINTFD(I->G,FB_Triangle)
                " TriangleFill-Debug: build single:     nTri=%d nActive=%d\n",I->nTri,I->nActive
                ENDFD;
              I->nActive--;
              i1 = I->activeEdge[I->nActive*2];
              i2 = I->activeEdge[I->nActive*2+1];
              TriangleBuildSingle(I,i1,i2,v,vn,n);
              PRINTFD(I->G,FB_Triangle)
                " TriangleFill-Debug: follow actives 2: nTri=%d nActive=%d\n",I->nTri,I->nActive
                ENDFD;
              FollowActives(I,v,vn,n,2);
            }
          }
    }
    
    PRINTFD(I->G,FB_Triangle)
                " TriangleFill-Debug: follow actives 4: nTri=%d nActive=%d\n",I->nTri,I->nActive
      ENDFD;

	 for(a=0;a<n;a++) 
		if(I->vertActive[a])
		  TriangleActivateEdges(I,a);
	 FollowActives(I,v,vn,n,4);

    PRINTFD(I->G,FB_Triangle)
      " TriangleFill-Debug: follow actives 5: nTri=%d nActive=%d\n",I->nTri,I->nActive
      ENDFD;

	 lastTri=I->nTri-1;
	 while(lastTri!=I->nTri) {
		lastTri=I->nTri;
		for(a=0;a<n;a++) 
		  if(I->vertActive[a])
			 TriangleActivateEdges(I,a);
		FollowActives(I,v,vn,n,5); /* this is a sloppy, forcing tesselation */
	 }
  }
    PRINTFD(I->G,FB_Triangle)
      " TriangleFill: leaving... nTri=%d nActive=%d\n",I->nTri,I->nActive
      ENDFD;

}

static void TriangleTxfFolds(TriangleSurfaceRec *II,float *v,float *vn,int n) 
{
  register TriangleSurfaceRec *I=II;
  int a,b,c,d,l,s01,s02,t1,t2;
  float *v0,*v1,*v2,*v3,d10[3],n10[3],d20[3],d30[3],d21[3],d31[3],d32[3];
  float x1020[3],x1030[3],x2132[3],x2032[3],s2[3],s3[3],nt[3];
  float old_dp, new_dp, old_conv,new_conv;
  for(a=0;a<n;a++) { /* first vertex */
    l=I->edgeStatus[a]; 
    while(l) {
      if( (s01 = I->link[l].value)<0 ) { /* closed edge */
        s01 = -s01;
        b=I->link[l].index; /* second vertex */
        v0 = v+a*3;
        v1 = v+b*3;
        c = I->edge[s01].vert3;
        d = I->edge[s01].vert4;
        v2 = v+c*3;
        v3 = v+d*3;
        subtract3f(v1,v0,d10);
        subtract3f(v2,v0,d20);
        subtract3f(v3,v0,d30);
        cross_product3f(d10,d20,x1020);
        cross_product3f(d10,d30,x1030);
        normalize3f(x1020);
        normalize3f(x1030);
        if( (old_dp = dot_product3f(x1020,x1030)) >0.5F) { /* triangles are nearly opposing one another */

          /*
            CGOLinewidth(I->G->DebugCGO,5.0);
            CGOBegin(I->G->DebugCGO,GL_LINES);
            CGOVertexv(I->G->DebugCGO,v0);
            CGOVertexv(I->G->DebugCGO,v1);
            CGOEnd(I->G->DebugCGO);
            CGOLinewidth(I->G->DebugCGO,3.0);
            CGOBegin(I->G->DebugCGO,GL_LINES);            
            CGOVertexv(I->G->DebugCGO,v2);
            CGOVertexv(I->G->DebugCGO,v3);
            CGOEnd(I->G->DebugCGO);
          */

          normalize23f(d10,n10);
          subtract3f(v2,v1,d21);
          subtract3f(v3,v1,d31);
          add3f(d21,d20,s2);
          add3f(d31,d30,s3);
          remove_component3f(s2,n10,s2);
          remove_component3f(s3,n10,s3);
          normalize3f(s2);
          normalize3f(s3);
          if(dot_product3f(s2,s3)>0.5F) {
            /* 2 & 3 on same side of 01 */
            subtract3f(v3,v2,d32);
            cross_product3f(d21,d32,x2132);
            cross_product3f(d20,d32,x2032);
            normalize3f(x2132);
            normalize3f(x2032);
            if( (new_dp = dot_product3f(x2132,x2032)) < old_dp) {
              int legal = true;
              s02 = TriangleEdgeStatus(I,a,d);
              if(s02<0) {
                s02 = -s02;
                if((I->edge[s02].vert3==c) || (I->edge[s02].vert4==c))
                  legal = false;
              }
              s02 = TriangleEdgeStatus(I,b,d);
              if(s02<0) {
                s02 = -s02;
                if((I->edge[s02].vert3==c) || (I->edge[s02].vert4==c))
                  legal = false;
              }
              s02 = TriangleEdgeStatus(I,a,c);
              if(s02<0) {
                s02 = -s02;
                if((I->edge[s02].vert3==d) || (I->edge[s02].vert4==d))
                  legal = false;
              }
              s02 = TriangleEdgeStatus(I,b,c);
              if(s02<0) {
                s02 = -s02;
                if((I->edge[s02].vert3==d) || (I->edge[s02].vert4==d))
                  legal = false;
              }
              if(legal) {

                /* how consistent are the normals (old versus new) ? */

                copy3f(vn+a*3,nt);
                add3f(vn+b*3,nt,nt);
                add3f(vn+c*3,nt,nt);
                old_conv = dot_product3f(x1020,nt);
                copy3f(vn+a*3,nt);
                add3f(vn+b*3,nt,nt);
                add3f(vn+d*3,nt,nt);
                old_conv += -dot_product3f(x1030,nt);
                old_conv = fabs(old_conv);

                copy3f(vn+a*3,nt);
                add3f(vn+c*3,nt,nt);
                add3f(vn+d*3,nt,nt);
                new_conv = dot_product3f(x2032,nt);
                copy3f(vn+b*3,nt);
                add3f(vn+c*3,nt,nt);
                add3f(vn+d*3,nt,nt);
                new_conv += -dot_product3f(x2132,nt);
                new_conv = fabs(new_conv);

                if((old_conv<new_conv)) {
                  /* switch the edges and triangles around */
                  
                  TriangleDeleteEdge(I,a,b);
                  TriangleEdgeSetStatus(I,c,d,-s01);
                  I->edge[s01].vert3 = a;
                  I->edge[s01].vert4 = b;
                  t1 = I->edge[s01].tri1;
                  t2 = I->edge[s01].tri2;
                  {
                    int i;
                    for(i=0;i<3;i++) {
                      if(I->tri[3*t1+i]==b) { /* a b c -> a c d */
                        I->tri[3*t1+i]=d;
                      }
                      if(I->tri[3*t2+i]==a) { /* a b d -> c b d */
                        I->tri[3*t2+i]=c;
                      }
                    }
                    TriangleRectify(I,t1,v,vn);
                    TriangleRectify(I,t2,v,vn);
                  }

                  s01 = TriangleEdgeStatus(I,a,d);
                  if(s01<0) {
                    s01 = -s01;
                    if(I->edge[s01].vert3==b) {
                      I->edge[s01].vert3=c;
                      I->edge[s01].tri1=t1;
                    } else if(I->edge[s01].vert4==b) {
                      I->edge[s01].vert4=c;
                      I->edge[s01].tri2=t1;
                    }
                  }
                  
                  s01 = TriangleEdgeStatus(I,a,c);
                  if(s01<0) {
                    s01 = -s01;
                    if(I->edge[s01].vert3==b) {
                      I->edge[s01].vert3=d;
                      I->edge[s01].tri1=t1;
                    } else if(I->edge[s01].vert4==b) {
                      I->edge[s01].vert4=d;
                      I->edge[s01].tri2=t1;
                    }
                  }
                  
                  s01 = TriangleEdgeStatus(I,b,c);
                  if(s01<0) {
                    s01 = -s01;
                    if(I->edge[s01].vert3==a) {
                      I->edge[s01].vert3=d;
                      I->edge[s01].tri1=t2;
                    } else if(I->edge[s01].vert4==a) {
                      I->edge[s01].vert4=d;
                      I->edge[s01].tri2=t2;
                    }
                  }
                  
                  s01 = TriangleEdgeStatus(I,b,d);
                  if(s01<0) {
                    s01 = -s01;
                    if(I->edge[s01].vert3==a) {
                      I->edge[s01].vert3=c;
                      I->edge[s01].tri1=t2;
                    } else if(I->edge[s01].vert4==a) {
                      I->edge[s01].vert4=c;
                      I->edge[s01].tri2=t2;
                    }
                  }

                  l=I->edgeStatus[a]; /* start vertex over since we've messed with its edges */
                }
              }
            }
          }
        }
      }
      l=I->link[l].next;
    }
  }
  /*  CGOStop(I->G->DebugCGO);*/
}

static void TriangleFixProblems(TriangleSurfaceRec *II,float *v,float *vn,int n) 
{
  register TriangleSurfaceRec *I=II;
  int problemFlag;
  int a,l,e;
  int i0,i1,i2,s01=0,s02=0,s12;
  int *pFlag = NULL;
  int *vFlag = NULL;
  problemFlag=false;

  pFlag=Alloc(int,n);  
  vFlag=Alloc(int,n);  
  for(a=0;a<n;a++) {
	 vFlag[a]=0;
	 if(I->vertActive[a]) {
		pFlag[a]=1;
		problemFlag=true;
	 } else {
		pFlag[a]=0;
	 }
  }
  if(problemFlag) 
	 {
		a=0;
		while(a<I->nTri) {
		  if(((pFlag[I->tri[a*3]]&&(pFlag[I->tri[a*3+1]]))||
				(pFlag[I->tri[a*3+1]]&&(pFlag[I->tri[a*3+2]]))||
				(pFlag[I->tri[a*3]]&&(pFlag[I->tri[a*3+2]]))))
			 {
				i0=I->tri[a*3];
				i1=I->tri[a*3+1];
				i2=I->tri[a*3+2];
				
				s01=TriangleEdgeStatus(I,i0,i1);
				if(s01<0) 
				  {
					 s01=-s01; 
					 if(I->edge[s01].tri2!=a) {
						I->edge[s01].tri1=I->edge[s01].tri2;
						I->edge[s01].vert3=I->edge[s01].vert4;
					 }
				  } else s01=0;
				TriangleEdgeSetStatus(I,i0,i1,s01);

				s02=TriangleEdgeStatus(I,i0,i2);
				if(s02<0) {
				  s02=-s02; 
					 if(I->edge[s02].tri2!=a) {
						I->edge[s02].tri1=I->edge[s02].tri2;
						I->edge[s02].vert3=I->edge[s02].vert4;
					 }
				} else s02=0;
				TriangleEdgeSetStatus(I,i0,i2,s02);

				s12=TriangleEdgeStatus(I,i1,i2);
				if(s12<0) {
				  s12=-s12; 
					 if(I->edge[s12].tri2!=a) {
						I->edge[s12].tri1=I->edge[s12].tri2;
						I->edge[s12].vert3=I->edge[s12].vert4;
					 }
				} else s12=0;
				TriangleEdgeSetStatus(I,i1,i2,s12);

				I->nTri--;
				TriangleMove(I,I->nTri,a);

				vFlag[i0]=true;
				vFlag[i1]=true;
				vFlag[i2]=true;
			 }
		  a++;
		}
		
		/* now go through the complicated step of resetting vertex activities */

		for(a=0;a<n;a++) 
		  if(vFlag[a]) 
			 I->vertActive[a]=-1;
		
		for(a=0;a<n;a++) {
		  l=I->edgeStatus[a]; 
		  while(l) {
			 if(I->link[l].value>0) {
				if(vFlag[a]) {
				  e=a;
				  if(I->vertActive[e]<0) I->vertActive[e]=0;
				  I->vertActive[e]++;
				} 
				if(vFlag[I->link[l].index]) {
				  e=I->link[l].index;
				  if(I->vertActive[e]<0) I->vertActive[e]=0;
				  I->vertActive[e]++;
				}
			 } if (I->link[l].value<0) {
				if(vFlag[a]) {
				  e=a;
				  if(I->vertActive[e]<0) I->vertActive[e]=0;
				} 
				if(vFlag[I->link[l].index]) {
				  e=I->link[l].index;
				  if(I->vertActive[e]<0) I->vertActive[e]=0;
				}
			 }
			 l=I->link[l].next;
		  }
		}
		
		TriangleAdjustNormals(I,v,vn,n,false);
        TriangleFill(I,v,vn,n,false);

	 }
  FreeP(vFlag);
  FreeP(pFlag);
}

static void TriangleBruteForceClosure(TriangleSurfaceRec *II,float *v,float *vn,int n,float cutoff) 
{
  register TriangleSurfaceRec *I=II;
  int a,b,c,d;
  int i0,i1,i2;
  float *v1,*v2,*v0,vt1[3],vt2[3],vt3[3],vt4[3],*n0,*n1,*n2,tNorm[3];
  int *pFlag = NULL;
  int *pair = NULL;
  int pc;
  int *active = NULL;
  int ac;
  int hits;
  int p1,p2;
  float dp;

  active = Alloc(int,n);
  ac = 0;
  pair = Alloc(int,n*2); 
  pc = 0;
  pFlag=Alloc(int,n);  
  for(a=0;a<n;a++) {
	 if(I->vertActive[a]) {
		pFlag[a]=1;
      active[ac]=a;
      ac++;
	 } else {
		pFlag[a]=0;
	 }
  }
  if(ac<80)  /* there is a limit to how much we can brute force... */
	 {
		a=0;
		while(a<I->nTri&&(pc<n)) {
        i0=I->tri[a*3];
        i1=I->tri[a*3+1];
        i2=I->tri[a*3+2];
        if(pFlag[i0]&&pFlag[i1]) {
          if(i0<i1) {
            pair[pc*2] = i0;
            pair[pc*2+1] = i1;
          } else {
            pair[pc*2] = i1;
            pair[pc*2+1] = i0;
          }
          pc++;
        }
        if(pFlag[i1]&&pFlag[i2]) {
          if(i1<i2) {
            pair[pc*2] = i1;
            pair[pc*2+1] = i2;
          } else {
            pair[pc*2] = i2;
            pair[pc*2+1] = i1;
          }
          pc++;
        }
        if(pFlag[i2]&&pFlag[i0]) {
          if(i2<i0) {
            pair[pc*2] = i2;
            pair[pc*2+1] = i0;
          } else {
            pair[pc*2] = i0;
            pair[pc*2+1] = i2;
          }
          pc++;
        }
		  a++;
      }
      PRINTFD(I->G,FB_Triangle)
        " Triangle-BFS: ac %d pc %d\n",ac,pc
        ENDFD;
        
      for(a=0;a<ac;a++) {
        i0 = active[a];
        for(b=a+1;b<ac;b++) {
          i1 = active[b];
          for(c=b+1;c<ac;c++) { /* consider all three-way possibilities */
            i2 = active[c];
            hits = 0;
            for(d=0;d<pc;d++) {
              p1 = *(pair + d*2 );
              p2 = *(pair + d*2 +1);
              if((p1==i0)&&(p2==i1))
                hits++;
              else if((p1==i1)&&(p2==i2))
                hits++;
              else if((p1==i0)&&(p2==i2))
                hits++;
            }
            if(hits>=3) {
              v0=v+i0*3; v1=v+i1*3; v2=v+i2*3;
              if(within3f(v0,v1,cutoff) &&
                 within3f(v1,v2,cutoff) &&
                 within3f(v0,v2,cutoff)) {
                
                n0 = vn+3*i0; n1 = vn+3*i1; n2 = vn+3*i2;							 
                add3f(n0,n1,vt1);
                add3f(n2,vt1,vt2);
                subtract3f(v1,v0,vt3);
                subtract3f(v2,v0,vt4);
                cross_product3f(vt3,vt4,tNorm); 
                normalize3f(tNorm); 							 
                dp = dot_product3f(vt2,tNorm);
                if(dp<0) scale3f(tNorm,-1.0F,tNorm);
                TriangleAdd(I,i0,i1,i2,tNorm,v,vn);
              }
            }
          }
        }
      }
	 }
  FreeP(active);
  FreeP(pair);
  FreeP(pFlag);
}



int *TrianglePointsToSurface(PyMOLGlobals *G,float *v,float *vn,int n,float cutoff,int *nTriPtr,int **stripPtr,float *extent)
{
  register TriangleSurfaceRec *I= NULL;
  int *result = NULL;
  MapType *map;
  int a;
  
  I = Alloc(TriangleSurfaceRec,1);
  if(I) {

    I->G=G;
    I->N=n;
    I->nActive = 0;
    I->activeEdge=VLAlloc(int,1000);

    I->link=VLAlloc(LinkType,n*2);
    I->nLink = 1;

    I->nEdge = 1;

    I->vNormal=VLAlloc(float,n*2);
    I->edge=VLAlloc(EdgeRec,n*2);

    I->tri=VLAlloc(int,n);
    I->nTri = 0;

    I->map=MapNew(I->G,cutoff,v,n,extent);
    MapSetupExpress(I->map);
    map=I->map;
    MapCacheInit(&I->map_cache,map,0,0);

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

    TriangleFill(I,v,vn,n,true);

    if(Feedback(G,FB_Triangle,FB_Debugging)) {
      for(a=0;a<n;a++) 
        if(I->vertActive[a])
          printf(" TrianglePTS-DEBUG: before fix %i %i\n",a,I->vertActive[a]);
    }

    TriangleTxfFolds(I,v,vn,n);

    TriangleFixProblems(I,v,vn,n);  

    if(Feedback(G,FB_Triangle,FB_Debugging)) {
      for(a=0;a<n;a++) 
        if(I->vertActive[a])
          printf(" TrianglePTS-DEBUG: after fix %i %i\n",a,I->vertActive[a]);
    }

    TriangleBruteForceClosure(I,v,vn,n,cutoff*3); /* abandon algorithm, just CLOSE THOSE GAPS! */

    TriangleAdjustNormals(I,v,vn,n,true);

    *(stripPtr) = TriangleMakeStripVLA(I,v,vn,n);

    (*nTriPtr)=I->nTri;
    VLAFreeP(I->activeEdge);
    VLAFreeP(I->link);
    VLAFreeP(I->vNormal);
    VLAFreeP(I->edge);
    FreeP(I->edgeStatus);
    FreeP(I->vertActive);
    FreeP(I->vertWeight);
    MapCacheFree(&I->map_cache,0,0);
    MapFree(map);

    result = I->tri;
  }
  FreeP(I);
  
  return(result);
}

