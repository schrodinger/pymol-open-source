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
-* 
-* 
-*
Z* -------------------------------------------------------------------
*/

#include"os_std.h"
#include"os_gl.h"

#include"Extrude.h"
#include"Base.h"
#include"OOMac.h"
#include"Setting.h"

void ExtrudeInit(CExtrude *I);

void ExtrudeInit(CExtrude *I)
{
  I->N = 0;
  I->p = NULL;
  I->n = NULL;
  I->c = NULL;
  
  I->sv = NULL; /* shape vertices */
  I->sn = NULL; /* shape normals */
  I->tv = NULL; /* transformed vertices */
  I->tn = NULL; /* transformed normals */
  I->Ns = 0; /* number of shape points */

}

void ExtrudeCircle(CExtrude *I, int n,float size)
{
  int a;
  float *v,*vn;

  if(n>20) n=20;
  
  I->sv = Alloc(float,3*(n+1));
  I->sn = Alloc(float,3*(n+1));
  I->tv = Alloc(float,3*(n+1));
  I->tn = Alloc(float,3*(n+1));
  I->Ns = n;
  
  v = I->sv;
  vn = I->sn;

  for(a=0;a<=n;a++)
	 {
      *(vn++) = 0.0;
      *(vn++) = cos(a*2*PI/n);
      *(vn++) = sin(a*2*PI/n);
      *(v++) = 0.0;
      *(v++) = cos(a*2*PI/n)*size;
      *(v++) = sin(a*2*PI/n)*size;
	 }
}

void ExtrudeRectangle(CExtrude *I,float thick,float width)
{
  float *v,*vn;

  I->Ns = 8;
  
  I->sv = Alloc(float,3*(I->Ns+1));
  I->sn = Alloc(float,3*(I->Ns+1));
  I->tv = Alloc(float,3*(I->Ns+1));
  I->tn = Alloc(float,3*(I->Ns+1));
  
  v = I->sv;
  vn = I->sn;

  *(vn++) = 0.0;
  *(vn++) = 1.0;
  *(vn++) = 0.0;
  *(vn++) = 0.0;
  *(vn++) = 1.0;
  *(vn++) = 0.0;
  *(v++) = 0.0;
  *(v++) = cos(PI/4)*thick;
  *(v++) = -sin(PI/4)*width;
  *(v++) = 0.0;
  *(v++) = cos(PI/4)*thick;
  *(v++) = sin(PI/4)*width;

  *(vn++) = 0.0;
  *(vn++) = 0.0;
  *(vn++) = 1.0;
  *(vn++) = 0.0;
  *(vn++) = 0.0;
  *(vn++) = 1.0;
  *(v++) = 0.0;
  *(v++) = cos(PI/4)*thick;
  *(v++) = sin(PI/4)*width;
  *(v++) = 0.0;
  *(v++) = -cos(PI/4)*thick;
  *(v++) = sin(PI/4)*width;

  *(vn++) = 0.0;
  *(vn++) = -1.0;
  *(vn++) = 0.0;
  *(vn++) = 0.0;
  *(vn++) = -1.0;
  *(vn++) = 0.0;
  *(v++) = 0.0;
  *(v++) = -cos(PI/4)*thick;
  *(v++) = sin(PI/4)*width;
  *(v++) = 0.0;
  *(v++) = -cos(PI/4)*thick;
  *(v++) = -sin(PI/4)*width;

  *(vn++) = 0.0;
  *(vn++) = 0.0;
  *(vn++) = -1.0;
  *(vn++) = 0.0;
  *(vn++) = 0.0;
  *(vn++) = -1.0;
  *(v++) = 0.0;
  *(v++) = -cos(PI/4)*thick;
  *(v++) = -sin(PI/4)*width;
  *(v++) = 0.0;
  *(v++) = cos(PI/4)*thick;
  *(v++) = -sin(PI/4)*width;


}

CExtrude *ExtrudeNew(void)
{
  OOAlloc(CExtrude);
  ExtrudeInit(I);
  return(I);
}

void ExtrudeBuildNormals1f(CExtrude *I)
{
  int a;
  float *v;
  if(I->N) {
    get_system1f3f(I->n,I->n+3,I->n+6); /* first is arbitrary */
    v = I->n+9;
    for(a=1;a<I->N;a++)
      {
        copy3f(v-6,v+3);
        get_system2f3f(v,v+3,v+6); /* the rest are relative to first */
        v+=9;
      }
  }
}

void ExtrudeBuildNormals2f(CExtrude *I)
{
  int a;
  float *v;
  if(I->N) {
    v = I->n;
    for(a=0;a<I->N;a++)
      {
        get_system2f3f(v,v+3,v+6); 
        v+=9;
      }
  }
}

void ExtrudeCGOTraceAxes(CExtrude *I,CGO *cgo)
{
  int a;
  float *v,*n;
  float v0[3];

  if(I->N) {
    CGOColor(cgo,0.5,0.5,0.5);
    CGOBegin(cgo,GL_LINES);
    v=I->p;
    n=I->n;
    for(a=0;a<I->N;a++) {
      add3f(v,n,v0);
      CGOVertexv(cgo,v0);      
      CGOVertexv(cgo,v);
      n+=3;
      add3f(v,n,v0);
      CGOVertexv(cgo,v0);      
      CGOVertexv(cgo,v);
      n+=3;
      add3f(v,n,v0);
      CGOVertexv(cgo,v0);      
      CGOVertexv(cgo,v);
      n+=3;
      v+=3;
    }
    CGOEnd(cgo);
  }
}

void ExtrudeCGOTrace(CExtrude *I,CGO *cgo)
{
  int a;
  float *v;

  if(I->N) {
    CGOColor(cgo,0.5,0.5,0.5);
    CGOBegin(cgo,GL_LINE_STRIP);
    v=I->p;
    for(a=0;a<I->N;a++) {
      CGOVertexv(cgo,v);
      v+=3;
    }
    CGOEnd(cgo);
  }
}

void ExtrudeComputeTangents(CExtrude *I)
{
  float *nv,*v1,*v;

  int a;
  nv = Alloc(float,I->N*3);

  v=nv;
  v1=I->p+3;
  
  for(a=1;a<I->N;a++)
    {
      subtract3f(v1,v1-3,v);
      normalize3f(v);
      v+=3;
      v1+=3;
    }
  
  /* compute tangents */
  
  v=nv;
  v1=I->n;

  *(v1++)=*(v++); /* first segment */
  *(v1++)=*(v++);
  *(v1++)=*(v++);
  v1+=6;

  for(a=1;a<(I->N-1);a++)
    {
      
      add3f(v,(v-3),v1);
      normalize3f(v1);		
      v1+=9;
      v+=3;
    }
  
  *(v1++)=*(v-3); /* last segment */
  *(v1++)=*(v-2);
  *(v1++)=*(v-1);

  FreeP(nv);

}



void ExtrudeCGOTraceFrame(CExtrude *I,CGO *cgo)
{
  int a,b;
  float *v;
  float *n;
  float *sv,*tv;
  float v0[3],v1[3];

  if(I->N&&I->Ns) {
    CGOColor(cgo,0.5,0.5,0.5);
    CGOBegin(cgo,GL_LINES);
    v=I->p;
    n=I->n;
    for(a=0;a<I->N;a++) {
      sv = I->sv;
      tv = I->tv;
      for(b=0;b<I->Ns;b++) {
        transform33Tf3f(n,sv,tv);
        sv+=3;
        tv+=3;
      }
      /* trace shape */
      tv = I->tv;
      add3f(v,tv,v0);
      for(b=1;b<I->Ns;b++) {
        tv+=3;
        add3f(v,tv,v1);
        CGOVertexv(cgo,v0);
        CGOVertexv(cgo,v1);
        copy3f(v1,v0);
      }
      tv = I->tv;
      add3f(v,tv,v1)
      CGOVertexv(cgo,v0);
      CGOVertexv(cgo,v1);
      v+=3;
      n+=9;
    }
    CGOEnd(cgo);
  }
}

void ExtrudeCGOSurfaceTube(CExtrude *I,CGO *cgo,int cap)
{
  int a,b;
  float *v;
  float *n;
  float *c;
  float *sv,*sn,*tv,*tn,*tv1,*tn1,*TV,*TN;
  float v0[3];
  
  if(I->N&&I->Ns) {

    TV=Alloc(float,3*(I->Ns+1)*I->N);
    TN=Alloc(float,3*(I->Ns+1)*I->N);
    
    /* compute transformed shape vertices */
    
    tn=TN;
    tv=TV;

    sv = I->sv;
    sn = I->sn;
    for(b=0;b<=I->Ns;b++) {
      if(b==I->Ns) {
        sv = I->sv;
        sn = I->sn;
      }
      v=I->p;
      n=I->n;
      
      for(a=0;a<I->N;a++) {
        transform33Tf3f(n,sv,tv);
        add3f(v,tv,tv);
        tv+=3;
        transform33Tf3f(n,sn,tn);
        tn+=3;
        n+=9;
        v+=3;
      }
      sv+=3;
      sn+=3;
    }
  
    /* fill in each strip separately */

    tv = TV;
    tn = TN;
    
    tv1 = TV+3*I->N;
    tn1 = TN+3*I->N;

    for(b=0;b<I->Ns;b++) {
      CGOBegin(cgo,GL_TRIANGLE_STRIP);
      c = I->c;
      for(a=0;a<I->N;a++) {
        CGOColorv(cgo,c);
        CGONormalv(cgo,tn);
        CGOVertexv(cgo,tv);
        tn+=3;
        tv+=3;
        CGONormalv(cgo,tn1);
        CGOVertexv(cgo,tv1);
        tn1+=3;
        tv1+=3;
        c+=3;
      }
      CGOEnd(cgo);
    }
    
    if(cap) {

      n = I->n;
      v = I->p;

      sv = I->sv;
      tv = I->tv;
      for(b=0;b<I->Ns;b++) {
        transform33Tf3f(n,sv,tv);
        add3f(v,tv,tv)
        sv+=3;
        tv+=3;
      }

      CGOBegin(cgo,GL_TRIANGLE_FAN);
      copy3f(I->n,v0);
      invert3f(v0);
      CGOColorv(cgo,I->c);
      CGONormalv(cgo,v0);
      CGOVertexv(cgo,v);
      /* trace shape */
      tv = I->tv;
      for(b=0;b<I->Ns;b++) {
        CGOVertexv(cgo,tv);
        tv+=3;
      }
      CGOVertexv(cgo,I->tv);
      CGOEnd(cgo);

      n = I->n+9*(I->N-1);
      v = I->p+3*(I->N-1);

      sv = I->sv;
      tv = I->tv;
      for(b=0;b<I->Ns;b++) {
        transform33Tf3f(n,sv,tv);
        add3f(v,tv,tv)
        sv+=3;
        tv+=3;
      }

      CGOBegin(cgo,GL_TRIANGLE_FAN);
      CGOColorv(cgo,I->c+3*(I->N-1));
      CGONormalv(cgo,n);
      CGOVertexv(cgo,v);
      /* trace shape */
      tv = I->tv;
      for(b=0;b<I->Ns;b++) {
        CGOVertexv(cgo,tv);
        tv+=3;
      }
      CGOVertexv(cgo,I->tv);
      CGOEnd(cgo);
      
    }
    FreeP(TV);
    FreeP(TN);
  }
  
}


void ExtrudeCGOSurfacePolygon(CExtrude *I,CGO *cgo,int cap)
{
  int a,b;
  float *v;
  float *n;
  float *c;
  float *sv,*sn,*tv,*tn,*tv1,*tn1,*TV,*TN;
  float v0[3];
  
  if(I->N&&I->Ns) {

    TV=Alloc(float,3*(I->Ns+1)*I->N);
    TN=Alloc(float,3*(I->Ns+1)*I->N);
    
    /* compute transformed shape vertices */
    
    tn=TN;
    tv=TV;

    sv = I->sv;
    sn = I->sn;
    for(b=0;b<=I->Ns;b++) {
      if(b==I->Ns) {
        sv = I->sv;
        sn = I->sn;
      }
      v=I->p;
      n=I->n;
      
      for(a=0;a<I->N;a++) {
        transform33Tf3f(n,sv,tv);
        add3f(v,tv,tv);
        tv+=3;
        transform33Tf3f(n,sn,tn);
        tn+=3;
        n+=9;
        v+=3;
      }
      sv+=3;
      sn+=3;
    }
  
    /* fill in each strip separately */

    tv = TV;
    tn = TN;
    
    tv1 = TV+3*I->N;
    tn1 = TN+3*I->N;

    for(b=0;b<I->Ns;b+=2) {
      if(SettingGet(cSetting_test1)<0.5)
        CGOBegin(cgo,GL_TRIANGLE_STRIP);
      else {
        CGOBegin(cgo,GL_LINE_STRIP);        
        CGODisable(cgo,GL_LIGHTING);
      }
      c = I->c;
      for(a=0;a<I->N;a++) {
        CGOColorv(cgo,c);
        CGONormalv(cgo,tn);
        CGOVertexv(cgo,tv);
        tn+=3;
        tv+=3;
        CGONormalv(cgo,tn1);
        CGOVertexv(cgo,tv1);
        tn1+=3;
        tv1+=3;
        c+=3;
      }
      tv+=3*I->N;
      tn+=3*I->N;
      tv1+=3*I->N;
      tn1+=3*I->N;
      CGOEnd(cgo);
    }
    
    if(cap) {

      n = I->n;
      v = I->p;

      sv = I->sv;
      tv = I->tv;
      for(b=0;b<I->Ns;b++) {
        transform33Tf3f(n,sv,tv);
        add3f(v,tv,tv)
        sv+=3;
        tv+=3;
      }
      CGOBegin(cgo,GL_TRIANGLE_FAN);
      copy3f(I->n,v0);
      invert3f(v0);
      CGOColorv(cgo,I->c);
      CGONormalv(cgo,v0);
      CGOVertexv(cgo,v);
      /* trace shape */
      tv = I->tv;
      for(b=0;b<I->Ns;b+=2) {
        CGOVertexv(cgo,tv);
        tv+=6;
      }
      CGOVertexv(cgo,I->tv);
      CGOEnd(cgo);

      n = I->n+9*(I->N-1);
      v = I->p+3*(I->N-1);

      sv = I->sv;
      tv = I->tv;
      for(b=0;b<I->Ns;b++) {
        transform33Tf3f(n,sv,tv);
        add3f(v,tv,tv)
        sv+=3;
        tv+=3;
      }
      CGOBegin(cgo,GL_TRIANGLE_FAN);
      CGOColorv(cgo,I->c+3*(I->N-1));
      CGONormalv(cgo,n);
      CGOVertexv(cgo,v);
      /* trace shape */
      tv = I->tv;
      for(b=0;b<I->Ns;b+=2) {
        CGOVertexv(cgo,tv);
        tv+=6;
      }
      CGOVertexv(cgo,I->tv);
      CGOEnd(cgo);
      

    }
    FreeP(TV);
    FreeP(TN);
  }
  
}



void ExtrudeTruncate(CExtrude *I,int n)
{
  I->N = n;
  /* should free RAM here... */
}

void ExtrudeAllocPointsNormalsColors(CExtrude *I,int n)
{
  if(I->N<n) {
    /* reset */
    FreeP(I->p);
    FreeP(I->n);
    FreeP(I->c);

    I->p = Alloc(float,3*n);
    I->n = Alloc(float,9*n);
    I->c = Alloc(float,3*n);
  }
  I->N = n;
}

void ExtrudeFree(CExtrude *I)
{
  FreeP(I->p);
  FreeP(I->n);
  FreeP(I->c);
  FreeP(I->tn);
  FreeP(I->tv);
  FreeP(I->sn);
  FreeP(I->sv);

  OOFreeP(I);
}


