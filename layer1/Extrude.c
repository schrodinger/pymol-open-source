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

#include"os_predef.h"
#include"os_std.h"
#include"os_gl.h"

#include"Extrude.h"
#include"Base.h"
#include"OOMac.h"
#include"Setting.h"
#include"Feedback.h"

void ExtrudeInit(CExtrude *I);

static float smooth(float x,float power)
{

  if(x<=0.5) {
    if(x<=0.0) x=0.0;
    return (0.5*pow(2.0*x,power));    
  } else {
    if(x>=1.0) x=1.0;
    return (1.0-(0.5*pow(2*(1.0-x),power)));
  }
}

CExtrude *ExtrudeCopyPointsNormalsColors(CExtrude *orig)
{
  OOAlloc(CExtrude);
  
  ExtrudeInit(I);

  ExtrudeAllocPointsNormalsColors(I,orig->N);

  CopyArray(I->p,orig->p,float,3*I->N);
  CopyArray(I->n,orig->n,float,9*I->N);
  CopyArray(I->c,orig->c,float,3*I->N);
  return(I);
}

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

  PRINTFD(FB_Extrude)
    " ExtrudeCircle-DEBUG: entered.\n"
    ENDFD;
  if(n>20) n=20;
  
  FreeP(I->sv);
  FreeP(I->sn);
  FreeP(I->tv);
  FreeP(I->tn);
  
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

  PRINTFD(FB_Extrude)
    " ExtrudeCircle-DEBUG: exiting...\n"
    ENDFD;

}

void ExtrudeOval(CExtrude *I, int n,float width,float length)
{
  int a;
  float *v,*vn;

  PRINTFD(FB_Extrude)
    " ExtrudeOval-DEBUG: entered.\n"
    ENDFD;

  if(n>20) n=20;
  
  FreeP(I->sv);
  FreeP(I->sn);
  FreeP(I->tv);
  FreeP(I->tn);
  
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
      *(vn++) = cos(a*2*PI/n)*length;
      *(vn++) = sin(a*2*PI/n)*width;
      *(v++) = 0.0;
      *(v++) = cos(a*2*PI/n)*width;
      *(v++) = sin(a*2*PI/n)*length;
	 }

  PRINTFD(FB_Extrude)
    " ExtrudeOval-DEBUG: exiting...\n"
    ENDFD;

}

void ExtrudeRectangle(CExtrude *I,float width,float length,int mode)
{
  float *v,*vn;

  PRINTFD(FB_Extrude)
    " ExtrudeRectangle-DEBUG: entered...\n"
    ENDFD;

  switch(mode) {
  case 0:
    I->Ns = 8;
    break;
  default:
    I->Ns = 4;
    break;
  }

  FreeP(I->sv);
  FreeP(I->sn);
  FreeP(I->tv);
  FreeP(I->tn);
  
  I->sv = Alloc(float,3*(I->Ns+1));
  I->sn = Alloc(float,3*(I->Ns+1));
  I->tv = Alloc(float,3*(I->Ns+1));
  I->tn = Alloc(float,3*(I->Ns+1));
  
  v = I->sv;
  vn = I->sn;

  if((!mode)||(mode==1)) {
    *(vn++) = 0.0;
    *(vn++) = 1.0;
    *(vn++) = 0.0;
    *(vn++) = 0.0;
    *(vn++) = 1.0;
    *(vn++) = 0.0;
    *(v++) = 0.0;
    *(v++) = cos(PI/4)*width;
    *(v++) = -sin(PI/4)*length;
    *(v++) = 0.0;
    *(v++) = cos(PI/4)*width;
    *(v++) = sin(PI/4)*length;
  }

  if((!mode)||(mode==2)) {  
    *(vn++) = 0.0;
    *(vn++) = 0.0;
    *(vn++) = 1.0;
    *(vn++) = 0.0;
    *(vn++) = 0.0;
    *(vn++) = 1.0;
    *(v++) = 0.0;
    *(v++) = cos(PI/4)*width;
    *(v++) = sin(PI/4)*length;
    *(v++) = 0.0;
    *(v++) = -cos(PI/4)*width;
    *(v++) = sin(PI/4)*length;
  }

  if((!mode)||(mode==1)) {
    *(vn++) = 0.0;
    *(vn++) = -1.0;
    *(vn++) = 0.0;
    *(vn++) = 0.0;
    *(vn++) = -1.0;
    *(vn++) = 0.0;
    *(v++) = 0.0;
    *(v++) = -cos(PI/4)*width;
    *(v++) = sin(PI/4)*length;
    *(v++) = 0.0;
    *(v++) = -cos(PI/4)*width;
    *(v++) = -sin(PI/4)*length;
  }

  if((!mode)||(mode==2)) {  
    
    *(vn++) = 0.0;
    *(vn++) = 0.0;
    *(vn++) = -1.0;
    *(vn++) = 0.0;
    *(vn++) = 0.0;
    *(vn++) = -1.0;
    *(v++) = 0.0;
    *(v++) = -cos(PI/4)*width;
    *(v++) = -sin(PI/4)*length;
    *(v++) = 0.0;
    *(v++) = cos(PI/4)*width;
    *(v++) = -sin(PI/4)*length;
  }

  PRINTFD(FB_Extrude)
    " ExtrudeRectangle-DEBUG: exiting...\n"
    ENDFD;

}




void ExtrudeDumbbell1(CExtrude *I,float width,float length,int mode)
{
  float *v,*vn;

  PRINTFD(FB_Extrude)
    " ExtrudeDumbbell1-DEBUG: entered...\n"
    ENDFD;

  switch(mode) {
  case 0:
    I->Ns = 4;
    break;
  default:
    I->Ns = 2;
    break;
    
  }

  FreeP(I->sv);
  FreeP(I->sn);
  FreeP(I->tv);
  FreeP(I->tn);
  
  I->sv = Alloc(float,3*(I->Ns+1));
  I->sn = Alloc(float,3*(I->Ns+1));
  I->tv = Alloc(float,3*(I->Ns+1));
  I->tn = Alloc(float,3*(I->Ns+1));
  
  v = I->sv;
  vn = I->sn;

  if((!mode)||(mode==1)) { /* top */
    *(vn++) = 0.0;
    *(vn++) = 1.0;
    *(vn++) = 0.0;
    *(vn++) = 0.0;
    *(vn++) = 1.0;
    *(vn++) = 0.0;
    *(v++) = 0.0;
    *(v++) = cos(PI/4)*width;
    *(v++) = -sin(PI/4)*length;
    *(v++) = 0.0;
    *(v++) = cos(PI/4)*width;
    *(v++) = sin(PI/4)*length;
  }
  
  if((!mode)||(mode==2)) { /* bottom */
    *(vn++) = 0.0;
    *(vn++) = -1.0;
    *(vn++) = 0.0;
    *(vn++) = 0.0;
    *(vn++) = -1.0;
    *(vn++) = 0.0;
    *(v++) = 0.0;
    *(v++) = -cos(PI/4)*width;
    *(v++) = sin(PI/4)*length;
    *(v++) = 0.0;
    *(v++) = -cos(PI/4)*width;
    *(v++) = -sin(PI/4)*length;
  }

  PRINTFD(FB_Extrude)
    " ExtrudeDumbbell1-DEBUG: exiting...\n"
    ENDFD;

}



void ExtrudeDumbbellEdge(CExtrude *I,int samp,int sign,float length)
{
  int a;
  float *n,*p,f,disp;

  PRINTFD(FB_Extrude)
    " ExtrudeDumbbellEdge-DEBUG: entered.\n"
    ENDFD;
  disp = (sign*sin(PI/4)*length);
  p=I->p;
  n=I->n;
  for(a=0;a<I->N;a++)
	 {
      if(a<=samp)
        f=disp*smooth((a/((float)samp)),2);
      else if(a>=(I->N-samp))
        f=disp*smooth(((I->N-a-1)/((float)samp)),2);
      else 
        f = disp;
      n+=6;
      (*p++) += *(n++)*f;
      (*p++) += *(n++)*f;
      (*p++) += *(n++)*f;
	 }
  PRINTFD(FB_Extrude)
    " ExtrudeDumbbellEdge-DEBUG: exiting...\n"
    ENDFD;

}



void ExtrudeDumbbell2(CExtrude *I, int n,int sign,float length,float size)
{
  int a;
  float *v,*vn;

  PRINTFD(FB_Extrude)
    " ExtrudeDumbbell2-DEBUG: entered.\n"
    ENDFD;
  if(n>20) n=20;
  
  FreeP(I->sv);
  FreeP(I->sn);
  FreeP(I->tv);
  FreeP(I->tn);
  
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
      *(v++) = (sin(a*2*PI/n)*size)+(sign*sin(PI/4)*length);
	 }

  PRINTFD(FB_Extrude)
    " ExtrudeDumbbell2-DEBUG: exiting...\n"
    ENDFD;

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

  PRINTFD(FB_Extrude)
    " ExtrudeBuildNormals1f-DEBUG: entered.\n"
    ENDFD;

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

  PRINTFD(FB_Extrude)
    " ExtrudeBuildNormals1f-DEBUG: exiting...\n"
    ENDFD;

}

void ExtrudeBuildNormals2f(CExtrude *I)
{
  int a;
  float *v;
  PRINTFD(FB_Extrude)
    " ExtrudeBuildNormals2f-DEBUG: entered.\n"
    ENDFD;

  if(I->N) {
    v = I->n;
    for(a=0;a<I->N;a++)
      {
        get_system2f3f(v,v+3,v+6); 
        v+=9;
      }
  }

  PRINTFD(FB_Extrude)
    " ExtrudeBuildNormals2f-DEBUG: entering...\n"
    ENDFD;

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

  PRINTFD(FB_Extrude)
    " ExtrudeComputeTangents-DEBUG: entered.\n"
    ENDFD;

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

  PRINTFD(FB_Extrude)
    " ExtrudeComputeTangents-DEBUG: exiting...\n"
    ENDFD;

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

void ExtrudeCGOSurfaceTube(CExtrude *I,CGO *cgo,int cap,float *color_override)
{
  int a,b;
  float *v;
  float *n;
  float *c;
  float *sv,*sn,*tv,*tn,*tv1,*tn1,*TV,*TN;
  float v0[3];
  int start,stop;
  PRINTFD(FB_Extrude)
    " ExtrudeCGOSurfaceTube-DEBUG: entered.\n"
    ENDFD;

  
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

    start=I->Ns/4;
    stop=3*I->Ns/4;
    for(b=0;b<I->Ns;b++) {
      if(SettingGet(cSetting_cartoon_debug)<1.5)
        CGOBegin(cgo,GL_TRIANGLE_STRIP);
      else {
        CGOBegin(cgo,GL_LINE_STRIP);        
        CGODisable(cgo,GL_LIGHTING);
      }
      c = I->c;
      for(a=0;a<I->N;a++) {
        if(color_override&&(b>start)&&(b<stop))
          CGOColorv(cgo,color_override);
        else
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

    if(SettingGet(cSetting_cartoon_debug)>=1.5) {
      CGOEnable(cgo,GL_LIGHTING);
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
      if(color_override)
        CGOColorv(cgo,color_override);
      else
        CGOColorv(cgo,I->c);
      CGONormalv(cgo,v0);
      CGOVertexv(cgo,v);
      /* trace shape */
      CGOVertexv(cgo,I->tv);
      for(b=I->Ns-1;b>=0;b--) {
        CGOVertexv(cgo,I->tv+b*3);
      }
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
      if(color_override)
        CGOColorv(cgo,color_override);
      else
        CGOColorv(cgo,I->c+3*(I->N-1));
      CGONormalv(cgo,n);
      CGOVertexv(cgo,v);
      /* trace shape */
      for(b=0;b<I->Ns;b++) {
        CGOVertexv(cgo,I->tv+b*3);
      }
      CGOVertexv(cgo,I->tv);
      CGOEnd(cgo);
      
    }
    FreeP(TV);
    FreeP(TN);
  }
  
  PRINTFD(FB_Extrude)
    " ExtrudeCGOSurfaceTube-DEBUG: exiting...\n"
    ENDFD;

}


void ExtrudeCGOSurfacePolygon(CExtrude *I,CGO *cgo,int cap,float *color_override)
{
  int a,b;
  float *v;
  float *n;
  float *c;
  float *sv,*sn,*tv,*tn,*tv1,*tn1,*TV,*TN;
  float v0[3];

  PRINTFD(FB_Extrude)
    " ExtrudeCGOSurfacePolygon-DEBUG: entered.\n"
    ENDFD;
  

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
      if(SettingGet(cSetting_cartoon_debug)<1.5)
        CGOBegin(cgo,GL_TRIANGLE_STRIP);
      else {
        CGOBegin(cgo,GL_LINE_STRIP);        
        CGODisable(cgo,GL_LIGHTING);
      }
      if(color_override)
        CGOColorv(cgo,color_override);
      c = I->c;
      for(a=0;a<I->N;a++) {
        if(!color_override)
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

    if(SettingGet(cSetting_cartoon_debug)>1.5) {
      CGOEnable(cgo,GL_LIGHTING);
    }

    if(cap) {

      if(color_override)
        CGOColorv(cgo,color_override);

      n = I->n;
      v = I->p;

      sv = I->sv;
      tv = I->tv;
      for(b=0;b<I->Ns;b++) {
        transform33Tf3f(n,sv,tv);
        add3f(v,tv,tv);
        sv+=3;
        tv+=3;
      }
      CGOBegin(cgo,GL_TRIANGLE_FAN);
      copy3f(I->n,v0);
      invert3f(v0);
      if(!color_override)
        CGOColorv(cgo,I->c);
      CGONormalv(cgo,v0);
      CGOVertexv(cgo,v);
      /* trace shape */
      CGOVertexv(cgo,I->tv);
      for(b=I->Ns-1;b>=0;b--) {
        CGOVertexv(cgo,I->tv+b*3);
      }
      CGOEnd(cgo);

      n = I->n+9*(I->N-1);
      v = I->p+3*(I->N-1);

      sv = I->sv;
      tv = I->tv;
      for(b=0;b<I->Ns;b++) {
        transform33Tf3f(n,sv,tv);
        add3f(v,tv,tv);
        sv+=3;
        tv+=3;
      }
      CGOBegin(cgo,GL_TRIANGLE_FAN);
      if(!color_override)
        CGOColorv(cgo,I->c+3*(I->N-1));
      CGONormalv(cgo,n);
      CGOVertexv(cgo,v);
      /* trace shape */
      for(b=0;b<I->Ns;b++) {
        CGOVertexv(cgo,I->tv+b*3);
      }
      CGOVertexv(cgo,I->tv);
      CGOEnd(cgo);

    }
    FreeP(TV);
    FreeP(TN);
  }
  
  PRINTFD(FB_Extrude)
    " ExtrudeCGOSurfacePolygon-DEBUG: exiting...\n"
    ENDFD;

}

void ExtrudeCGOSurfacePolygonTaper(CExtrude *I,CGO *cgo,int sampling,float *color_override)
{
  int a,b;
  float *v;
  float *n;
  float *c;
  float *sv,*sn,*tv,*tn,*tv1,*tn1,*TV,*TN;
  float s0[3];
  float f;
  int subN;

  subN=I->N-sampling;

  PRINTFD(FB_Extrude)
    " ExtrudeCGOSurfacePolygon-DEBUG: entered.\n"
    ENDFD;
  

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
        if((a>=sampling)&&(a<subN)) {

          transform33Tf3f(n,sv,tv);
          add3f(v,tv,tv);
          tv+=3;
          transform33Tf3f(n,sn,tn);
          tn+=3;
          n+=9;
          v+=3;
        } else {
          copy3f(sv,s0);

          if(a>=subN) {
            f = ((I->N-a-1)/((float)sampling));
          } else if(a<sampling) {
            f = (a/((float)sampling));
          } else 
            f = 1.0; 
          f=smooth(f,2);
          s0[2]*=f;

          transform33Tf3f(n,s0,tv);
          add3f(v,tv,tv);
          tv+=3;
          transform33Tf3f(n,sn,tn);
          tn+=3;
          n+=9;
          v+=3;
          
        }
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
      if(SettingGet(cSetting_cartoon_debug)<1.5)
        CGOBegin(cgo,GL_TRIANGLE_STRIP);
      else {
        CGOBegin(cgo,GL_LINE_STRIP);        
        CGODisable(cgo,GL_LIGHTING);
      }
      if(color_override)
        CGOColorv(cgo,color_override);
      c = I->c;
      for(a=0;a<I->N;a++) {
        if(!color_override)
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

    if(SettingGet(cSetting_cartoon_debug)>1.5) {
      CGOEnable(cgo,GL_LIGHTING);
    }

    FreeP(TV);
    FreeP(TN);
  }
  
  PRINTFD(FB_Extrude)
    " ExtrudeCGOSurfacePolygon-DEBUG: exiting...\n"
    ENDFD;
  
}





void ExtrudeCGOSurfaceStrand(CExtrude *I,CGO *cgo,int sampling,float *color_override)
{
  int a,b;
  float *v;
  float *n;
  float *c;
  float *sv,*sn,*tv,*tn,*tv1,*tn1,*TV,*TN;
  float v0[3],n0[3],s0[3],z[3]={1.0,0.0,1.0};
  int subN;
  
  subN=I->N-sampling;

  PRINTFD(FB_Extrude)
    " ExtrudeCGOSurfaceStrand-DEBUG: entered.\n"
    ENDFD;
  

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
        copy3f(sv,s0);
        if(a==subN) {
          scale3f(s0,0.50,s0);
        }
        transform33Tf3f(n,s0,tv);
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
      if(SettingGet(cSetting_cartoon_debug)<1.5)
        CGOBegin(cgo,GL_TRIANGLE_STRIP);
      else {
        CGOBegin(cgo,GL_LINE_STRIP);        
        CGODisable(cgo,GL_LIGHTING);
      }
      c = I->c;
      for(a=0;a<I->N;a++) {
        if(a<=subN) {
          if(color_override&&((b==2)||(b==3)||(b==6)||(b==7)))
            CGOColorv(cgo,color_override);
          else
            CGOColorv(cgo,c);
          CGONormalv(cgo,tn);
          CGOVertexv(cgo,tv);
        }
        tn+=3;
        tv+=3;
        if(a<=subN) {
          CGONormalv(cgo,tn1);
          CGOVertexv(cgo,tv1);
        }
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

    if(SettingGet(cSetting_cartoon_debug)>1.5) {
      CGOEnable(cgo,GL_LIGHTING);
    }

    if(1) {

      n = I->n;
      v = I->p;

      sv = I->sv;
      tv = I->tv;
      for(b=0;b<I->Ns;b++) {
        transform33Tf3f(n,sv,tv);
        add3f(v,tv,tv);
        sv+=3;
        tv+=3;
      }
      CGOBegin(cgo,GL_TRIANGLE_FAN);
      copy3f(I->n,v0);
      invert3f(v0);
      if(color_override)
        CGOColorv(cgo,color_override);
      else
        CGOColorv(cgo,I->c);
      CGONormalv(cgo,v0);
      CGOVertexv(cgo,v);
      /* trace shape */
      CGOVertexv(cgo,I->tv);
      for(b=I->Ns-2;b>=0;b-=2) {
        CGOVertexv(cgo,I->tv+b*3);
      }
      CGOEnd(cgo);
    }

    /* now do the arrow part */

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
        copy3f(sv,s0);
        s0[2]=s0[2]*((1.5*((I->N-1)-a))/sampling);
        transform33Tf3f(n,s0,tv);
        add3f(v,tv,tv);
        tv+=3;
        copy3f(sn,n0);
        if(fabs(dot_product3f(sn,z))>R_SMALL4) {
          n0[0]+=0.4;
          normalize3f(n0);
        }
        transform33Tf3f(n,n0,tn);
        tn+=3;
        n+=9;
        v+=3;
      }
      sv+=3;
      sn+=3;
    }

    tv = TV;
    tn = TN;
    
    tv1 = TV+3*I->N;
    tn1 = TN+3*I->N;

    for(b=0;b<I->Ns;b+=2) {
      if(SettingGet(cSetting_cartoon_debug)<1.5)
        CGOBegin(cgo,GL_TRIANGLE_STRIP);
      else {
        CGOBegin(cgo,GL_LINE_STRIP);        
        CGODisable(cgo,GL_LIGHTING);
      }
      c = I->c;
      for(a=0;a<I->N;a++) {
        if(a>=(subN-1)) {
          if(color_override&&((b==2)||(b==3)||(b==6)||(b==7)))
            CGOColorv(cgo,color_override);
          else
            CGOColorv(cgo,c);
          CGONormalv(cgo,tn);
          CGOVertexv(cgo,tv);
        }
        tn+=3;
        tv+=3;
        if(a>=(subN-1)) {
          CGONormalv(cgo,tn1);
          CGOVertexv(cgo,tv1);
        }
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
    
    n = I->n+9*(subN-1);
    v = I->p+3*(subN-1);

    sv = I->sv;
    tv = I->tv;

    for(b=0;b<I->Ns;b++) {
      copy3f(sv,s0);
      s0[2]=s0[2]*1.5;
      transform33Tf3f(n,s0,tv);
      add3f(v,tv,tv);
      sv+=3;
      tv+=3;
    }

    CGOBegin(cgo,GL_TRIANGLE_FAN);
    copy3f(n,v0);
    invert3f(v0);
    if(color_override)
      CGOColorv(cgo,color_override);
    else
      CGOColorv(cgo,I->c+3*(subN-1));
    CGONormalv(cgo,v0);
    CGOVertexv(cgo,v);

    /* trace shape */
    tv = I->tv;
    for(b=0;b<I->Ns;b+=2) {
      CGOVertexv(cgo,I->tv+b*3);
    }
    CGOVertexv(cgo,I->tv);
    CGOEnd(cgo);

    FreeP(TV);
    FreeP(TN);
  }
  
  PRINTFD(FB_Extrude)
    " ExtrudeCGOSurfaceStrand-DEBUG: exiting...\n"
    ENDFD;

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

    I->p = Alloc(float,3*(n+1));
    I->n = Alloc(float,9*(n+1));
    I->c = Alloc(float,3*(n+1));
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


