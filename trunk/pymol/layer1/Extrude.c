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
-* Cameron Mura
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

void ExtrudeInit(PyMOLGlobals *G,CExtrude *I);

static float smooth(float x,float power)
{

  if(x<=0.5) {
    if(x<=0.0) x=0.0;
    return ((float)(0.5F*pow(2.0*x,power)));    
  } else {
    if(x>=1.0) x=1.0;
    return ((float)(1.0F-(0.5*pow(2*(1.0-x),power))));
  }
}

CExtrude *ExtrudeCopyPointsNormalsColors(CExtrude *orig)
{
  OOAlloc(orig->G,CExtrude);
  
  ExtrudeInit(orig->G,I);

  ExtrudeAllocPointsNormalsColors(I,orig->N);

  CopyArray(I->p,orig->p,float,3*I->N);
  CopyArray(I->n,orig->n,float,9*I->N);
  CopyArray(I->c,orig->c,float,3*I->N);
  CopyArray(I->i,orig->i,int,    I->N);
  CopyArray(I->sf,orig->sf,float,I->N); /* PUTTY: scale factors */
     
  return(I);
}

void ExtrudeInit(PyMOLGlobals *G,CExtrude *I)
{
  I->G=G;

  I->N = 0;
  I->p = NULL;
  I->n = NULL;
  I->c = NULL;
  I->i = NULL;

  I->sv = NULL; /* shape vertices */
  I->sn = NULL; /* shape normals */
  I->tv = NULL; /* transformed vertices */
  I->tn = NULL; /* transformed normals */
  I->Ns = 0; /* number of shape points */

  I->sf = NULL;
}

void ExtrudeCircle(CExtrude *I, int n,float size)
{
  int a;
  float *v,*vn;

  PRINTFD(I->G,FB_Extrude)
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
  I->r = size;

  v = I->sv;
  vn = I->sn;

  for(a=0;a<=n;a++)
	 {
      *(vn++) = 0.0;
      *(vn++) = (float)cos(a*2*PI/n);
      *(vn++) = (float)sin(a*2*PI/n);
      *(v++) = 0.0;
      *(v++) = (float)cos(a*2*PI/n)*size;
      *(v++) = (float)sin(a*2*PI/n)*size;
	 }

  PRINTFD(I->G,FB_Extrude)
    " ExtrudeCircle-DEBUG: exiting...\n"
    ENDFD;

}

void ExtrudeOval(CExtrude *I, int n,float width,float length)
{
  int a;
  float *v,*vn;

  PRINTFD(I->G,FB_Extrude)
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
      *(vn++) = (float)cos(a*2*PI/n)*length;
      *(vn++) = (float)sin(a*2*PI/n)*width;
      *(v++) = 0.0;
      *(v++) = (float)cos(a*2*PI/n)*width;
      *(v++) = (float)sin(a*2*PI/n)*length;
	 }

  PRINTFD(I->G,FB_Extrude)
    " ExtrudeOval-DEBUG: exiting...\n"
    ENDFD;

}

void ExtrudeRectangle(CExtrude *I,float width,float length,int mode)
{
  float *v,*vn;

  PRINTFD(I->G,FB_Extrude)
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
    *(v++) = (float)cos(PI/4)*width;
    *(v++) = (float)-sin(PI/4)*length;
    *(v++) = 0.0;
    *(v++) = (float)cos(PI/4)*width;
    *(v++) = (float)sin(PI/4)*length;
  }

  if((!mode)||(mode==2)) {  
    *(vn++) = 0.0;
    *(vn++) = 0.0;
    *(vn++) = 1.0;
    *(vn++) = 0.0;
    *(vn++) = 0.0;
    *(vn++) = 1.0;
    *(v++) = 0.0;
    *(v++) = (float)cos(PI/4)*width;
    *(v++) = (float)sin(PI/4)*length;
    *(v++) = 0.0;
    *(v++) = (float)-cos(PI/4)*width;
    *(v++) = (float)sin(PI/4)*length;
  }

  if((!mode)||(mode==1)) {
    *(vn++) = 0.0;
    *(vn++) = -1.0;
    *(vn++) = 0.0;
    *(vn++) = 0.0;
    *(vn++) = -1.0;
    *(vn++) = 0.0;
    *(v++) = 0.0;
    *(v++) = (float)-cos(PI/4)*width;
    *(v++) = (float)sin(PI/4)*length;
    *(v++) = 0.0;
    *(v++) = (float)-cos(PI/4)*width;
    *(v++) = (float)-sin(PI/4)*length;
  }

  if((!mode)||(mode==2)) {  
    
    *(vn++) = 0.0;
    *(vn++) = 0.0;
    *(vn++) = -1.0;
    *(vn++) = 0.0;
    *(vn++) = 0.0;
    *(vn++) = -1.0;
    *(v++) = 0.0;
    *(v++) = (float)-cos(PI/4)*width;
    *(v++) = (float)-sin(PI/4)*length;
    *(v++) = 0.0;
    *(v++) = (float)cos(PI/4)*width;
    *(v++) = (float)-sin(PI/4)*length;
  }

  PRINTFD(I->G,FB_Extrude)
    " ExtrudeRectangle-DEBUG: exiting...\n"
    ENDFD;

}




void ExtrudeDumbbell1(CExtrude *I,float width,float length,int mode)
{
  float *v,*vn;

  PRINTFD(I->G,FB_Extrude)
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
    *(v++) = (float)cos(PI/4)*width;
    *(v++) = (float)-sin(PI/4)*length;
    *(v++) = 0.0;
    *(v++) = (float)cos(PI/4)*width;
    *(v++) = (float)sin(PI/4)*length;
  }
  
  if((!mode)||(mode==2)) { /* bottom */
    *(vn++) = 0.0;
    *(vn++) = -1.0;
    *(vn++) = 0.0;
    *(vn++) = 0.0;
    *(vn++) = -1.0;
    *(vn++) = 0.0;
    *(v++) = 0.0;
    *(v++) = (float)-cos(PI/4)*width;
    *(v++) = (float)sin(PI/4)*length;
    *(v++) = 0.0;
    *(v++) = (float)-cos(PI/4)*width;
    *(v++) = (float)-sin(PI/4)*length;
  }

  PRINTFD(I->G,FB_Extrude)
    " ExtrudeDumbbell1-DEBUG: exiting...\n"
    ENDFD;

}



void ExtrudeDumbbellEdge(CExtrude *I,int samp,int sign,float length)
{
  int a;
  float *n,*p,f,disp;

  PRINTFD(I->G,FB_Extrude)
    " ExtrudeDumbbellEdge-DEBUG: entered.\n"
    ENDFD;
  disp = (float)(sign*sin(PI/4)*length);
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
  PRINTFD(I->G,FB_Extrude)
    " ExtrudeDumbbellEdge-DEBUG: exiting...\n"
    ENDFD;

}



void ExtrudeDumbbell2(CExtrude *I, int n,int sign,float length,float size)
{
  int a;
  float *v,*vn;

  PRINTFD(I->G,FB_Extrude)
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
      *(vn++) = (float)cos(a*2*PI/n);
      *(vn++) = (float)sin(a*2*PI/n);
      *(v++) = 0.0;
      *(v++) = (float)cos(a*2*PI/n)*size;
      *(v++) = (float)((sin(a*2*PI/n)*size)+(sign*sin(PI/4)*length));
	 }

  PRINTFD(I->G,FB_Extrude)
    " ExtrudeDumbbell2-DEBUG: exiting...\n"
    ENDFD;

}


CExtrude *ExtrudeNew(PyMOLGlobals *G)
{
  OOAlloc(G,CExtrude);
  ExtrudeInit(G,I);
  return(I);
}

void ExtrudeBuildNormals1f(CExtrude *I)
{
  int a;
  float *v;

  PRINTFD(I->G,FB_Extrude)
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

  PRINTFD(I->G,FB_Extrude)
    " ExtrudeBuildNormals1f-DEBUG: exiting...\n"
    ENDFD;

}

void ExtrudeBuildNormals2f(CExtrude *I)
{
  int a;
  float *v;
  PRINTFD(I->G,FB_Extrude)
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

  PRINTFD(I->G,FB_Extrude)
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

  PRINTFD(I->G,FB_Extrude)
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

  PRINTFD(I->G,FB_Extrude)
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
  int a,b,*i;
  float *v;
  float *n;
  float *c;
  float *sv,*sn,*tv,*tn,*tv1,*tn1,*TV,*TN;
  float v0[3];
  int start,stop;
  PRINTFD(I->G,FB_Extrude)
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
      if(SettingGet(I->G,cSetting_cartoon_debug)<1.5)
        CGOBegin(cgo,GL_TRIANGLE_STRIP);
      else {
        CGOBegin(cgo,GL_LINE_STRIP);        
        CGODisable(cgo,GL_LIGHTING);
      }
      c = I->c;
      i = I->i;
      for(a=0;a<I->N;a++) {
        if(color_override&&(b>start)&&(b<stop))
          CGOColorv(cgo,color_override);
        else
          CGOColorv(cgo,c);
        CGOPickColor(cgo,*i,cPickableAtom);
        CGONormalv(cgo,tn);
        CGOVertexv(cgo,tv);
        tn+=3;
        tv+=3;
        CGONormalv(cgo,tn1);
        CGOVertexv(cgo,tv1);
        tn1+=3;
        tv1+=3;
        c+=3;
        i++;
      }
      CGOEnd(cgo);
    }

    if(SettingGet(I->G,cSetting_cartoon_debug)>=1.5) {
      CGOEnable(cgo,GL_LIGHTING);
    }
    
    switch(cap) {
    case 1:
      
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
      CGOPickColor(cgo,I->i[0],cPickableAtom);
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
      CGOPickColor(cgo,I->i[I->N-1],cPickableAtom);
      CGONormalv(cgo,n);
      CGOVertexv(cgo,v);
      /* trace shape */
      for(b=0;b<I->Ns;b++) {
        CGOVertexv(cgo,I->tv+b*3);
      }
      CGOVertexv(cgo,I->tv);
      CGOEnd(cgo);
      break;
    case 2:
      {
        float factor = (float) cos(0.75*PI/I->Ns);
        if(color_override)
          CGOColorv(cgo,color_override);
        else
          CGOColorv(cgo,I->c);
        CGOPickColor(cgo,I->i[0],cPickableAtom);
        v = I->p;
        CGOSphere(cgo,v,I->r*factor);
        v = I->p+3*(I->N-1);
        if(color_override)
          CGOColorv(cgo,color_override);
        else
          CGOColorv(cgo,I->c+3*(I->N-1));
        CGOPickColor(cgo,I->i[I->N-1],cPickableAtom);
        CGOSphere(cgo,v,I->r*factor); 
      }
      break;
    }
    FreeP(TV);
    FreeP(TN);
  }
  
  PRINTFD(I->G,FB_Extrude)
    " ExtrudeCGOSurfaceTube-DEBUG: exiting...\n"
    ENDFD;

}


void ExtrudeCGOSurfaceVariableTube(CExtrude *I,CGO *cgo,int cap)
{
  int a,b,*i;
  float *v;
  float *n;
  float *c;
  float *sv,*sn,*tv,*tn,*tv1,*tn1,*TV=NULL,*TN=NULL,*AN=NULL,*an;
  float v0[3];
  float *sf;                            /* PUTTY: scale factor from ExtrudeMakeSausLUT() */
  int start,stop;
  PRINTFD(I->G,FB_Extrude)
    " ExtrudeCGOSurfaceTube-DEBUG: entered.\n"
    ENDFD;

  
  if(I->N&&I->Ns) {

    TV=Alloc(float,3*(I->Ns+1)*I->N);
    TN=Alloc(float,3*(I->Ns+1)*I->N);
    AN = Alloc(float,3*I->N); /* normals adjusted for changing widths */

    /* compute transformed shape vertices */
    
    tv=TV;

    sv = I->sv;
    for(b=0;b<=I->Ns;b++) {
      if(b==I->Ns) {
        sv = I->sv;
      }

      n=I->n; /* NOTE: n is not a counter -- it's a 3x3 coordinate system! */
      v=I->p;
      sf=I->sf;                         /* PUTTY: scale factors */

      for(a=0;a<I->N;a++) {
        transform33Tf3f(n,sv,tv);
        
        *(tv) *= *sf;
        *(tv+1) *= *sf;
        *(tv+2) *= *sf;
        
        add3f(v,tv,tv);
        tv+=3;
        v+=3;
        sf++;
        n+=9;
      }
      sv+=3;
    }

    /* compute transformed normals, taking into account changing radii */

    tn=TN;
    tv=TV;
    
    sn = I->sn;
    for(b=0;b<=I->Ns;b++) {

      float d1,d2,r0,r1,r2,x1,x2;        
      
      if(b==I->Ns) {
        sn = I->sn;
      }

      an = AN;
      v=I->p;

      for(a=0;a<I->N;a++) { 
        if((a>0) && (a<(I->N-1))) {
          /* compute rises */

          r0 = (float)diff3f(v,tv);
          r1 = (float)(diff3f(v-3,tv-3) - r0);
          r2 = (float)(diff3f(v+3,tv+3) - r0);
          
          /* compute runs */
          
          d1 = (float)diff3f(v-3,v);
          d2 = (float)diff3f(v+3,v);
          
          /* compute x-to-yz weights */
          
          x1 = r1/d1;
          x2 = -r2/d2;
     
          if(a==1) {
            an[-3] = x1;
            an[-2] = sn[1];
            an[-1] = sn[2];
            normalize3f(an-3);
          } else if(a==I->N-2) {
            an[3] = x2;
            an[4] = sn[1];
            an[5] = sn[2];
            normalize3f(an+3);
          }
          an[0] = (x1 + x2)/2.0F;
          an[1] = sn[1];
          an[2] = sn[2];
          normalize3f(an);
        }
        tv+=3;
        v+=3;
        an+=3;
      }
   
      n=I->n; /* NOTE: n is not a counter -- it's a 3x3 coordinate system! */
      an=AN;

      for(a=0;a<I->N;a++) {
        transform33Tf3f(n,an,tn);
        tn+=3;
        an+=3;
        n+=9;
      }
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
      if(SettingGet(I->G,cSetting_cartoon_debug)<1.5)
        CGOBegin(cgo,GL_TRIANGLE_STRIP);
      else {
        CGOBegin(cgo,GL_LINE_STRIP);        
        CGODisable(cgo,GL_LIGHTING);
      }
      c = I->c;
      i = I->i;
      for(a=0;a<I->N;a++) {
        CGOColorv(cgo,c);
        CGOPickColor(cgo,*i,cPickableAtom);
        CGONormalv(cgo,tn);
        CGOVertexv(cgo,tv);
        tn+=3;
        tv+=3;
        CGONormalv(cgo,tn1);
        CGOVertexv(cgo,tv1);
        tn1+=3;
        tv1+=3;
        c+=3;
        i++;
      }
      CGOEnd(cgo);
    }

    if(SettingGet(I->G,cSetting_cartoon_debug)>3.5) {

      tv = TV;
      tn = TN;
      
      tv1 = TV+3*I->N;
      tn1 = TN+3*I->N;
      
      start=I->Ns/4;
      stop=3*I->Ns/4;
      for(b=0;b<I->Ns;b++) {
        float vv[3];
        CGOBegin(cgo,GL_LINES);
        c = I->c;
        i = I->i;
        for(a=0;a<I->N;a++) {
          CGOColorv(cgo,c);
          copy3f(tn,vv);
          scale3f(vv,0.3F,vv)
          add3f(vv,tv,vv);
          CGONormalv(cgo,tn);
          CGOVertexv(cgo,tv);
          CGOVertexv(cgo,vv);
          tn+=3;
          tv+=3;
          copy3f(tn1,vv);
          scale3f(vv,0.3F,vv)
          add3f(vv,tv1,vv);
          CGONormalv(cgo,tn1);
          CGOVertexv(cgo,tv1);
          CGOVertexv(cgo,vv);
          tn1+=3;
          tv1+=3;
          c+=3;
          i++;
        }
        CGOEnd(cgo);
      }
    }
      
    if(SettingGet(I->G,cSetting_cartoon_debug)>=1.5) {
      CGOEnable(cgo,GL_LIGHTING);
    }
    
    if(cap) {

      n = I->n;
      v = I->p;
      sf = I->sf;

      sv = I->sv;
      tv = I->tv;
      for(b=0;b<I->Ns;b++) {
        transform33Tf3f(n,sv,tv);

        *(tv) *= *sf;         
        *(tv+1) *= *sf;
        *(tv+2) *= *sf;

        add3f(v,tv,tv);
        sv+=3;
        tv+=3;
      }
      
      CGOBegin(cgo,GL_TRIANGLE_FAN);
      copy3f(I->n,v0);
      invert3f(v0);
      CGOColorv(cgo,I->c);
      CGOPickColor(cgo,I->i[0],cPickableAtom);
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
      sf = I->sf+(I->N-1);		/* PUTTY */

      sv = I->sv;
      tv = I->tv;
      for(b=0;b<I->Ns;b++) {
        transform33Tf3f(n,sv,tv);

        *(tv) *= *(sf);         
        *(tv+1) *= *(sf);
        *(tv+2) *= *(sf);

        add3f(v,tv,tv)
        sv+=3;
        tv+=3;
      }

      CGOBegin(cgo,GL_TRIANGLE_FAN);
      CGOColorv(cgo,I->c+3*(I->N-1));
      CGOPickColor(cgo,I->i[I->N-1],cPickableAtom);
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
    FreeP(AN);
  }
  
  PRINTFD(I->G,FB_Extrude)
    " ExtrudeCGOSurfaceTube-DEBUG: exiting...\n"
    ENDFD;

}



void ExtrudeCGOSurfacePolygon(CExtrude *I,CGO *cgo,int cap,float *color_override)
{
  int a,b,*i;
  float *v;
  float *n;
  float *c;
  float *sv,*sn,*tv,*tn,*tv1,*tn1,*TV,*TN;
  float v0[3];

  PRINTFD(I->G,FB_Extrude)
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
      if(SettingGet(I->G,cSetting_cartoon_debug)<1.5)
        CGOBegin(cgo,GL_TRIANGLE_STRIP);
      else {
        CGOBegin(cgo,GL_LINE_STRIP);        
        CGODisable(cgo,GL_LIGHTING);
      }
      if(color_override)
        CGOColorv(cgo,color_override);
      c = I->c;
      i = I->i;
      for(a=0;a<I->N;a++) {
        if(!color_override)
          CGOColorv(cgo,c);
        CGOPickColor(cgo,*i,cPickableAtom);
        CGONormalv(cgo,tn);
        CGOVertexv(cgo,tv);
        tn+=3;
        tv+=3;
        CGONormalv(cgo,tn1);
        CGOVertexv(cgo,tv1);
        tn1+=3;
        tv1+=3;
        c+=3;
        i++;
      }
      tv+=3*I->N;
      tn+=3*I->N;
      tv1+=3*I->N;
      tn1+=3*I->N;
      CGOEnd(cgo);
    }

    if(SettingGet(I->G,cSetting_cartoon_debug)>1.5) {
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
      CGOPickColor(cgo,I->i[0],cPickableAtom);
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
      CGOPickColor(cgo,I->i[I->N-1],cPickableAtom);
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
  
  PRINTFD(I->G,FB_Extrude)
    " ExtrudeCGOSurfacePolygon-DEBUG: exiting...\n"
    ENDFD;

}

void ExtrudeCGOSurfacePolygonTaper(CExtrude *I,CGO *cgo,int sampling,float *color_override)
{
  int a,b,*i;
  float *v;
  float *n;
  float *c;
  float *sv,*sn,*tv,*tn,*tv1,*tn1,*TV,*TN;
  float s0[3];
  float f;
  int subN;

  subN=I->N-sampling;

  PRINTFD(I->G,FB_Extrude)
    " ExtrudeCGOSurfacePolygonTaper-DEBUG: entered.\n"
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
      if(SettingGet(I->G,cSetting_cartoon_debug)<1.5)
        CGOBegin(cgo,GL_TRIANGLE_STRIP);
      else {
        CGOBegin(cgo,GL_LINE_STRIP);        
        CGODisable(cgo,GL_LIGHTING);
      }
      if(color_override)
        CGOColorv(cgo,color_override);
      c = I->c;
      i = I->i;
      for(a=0;a<I->N;a++) {
        if(!color_override)
          CGOColorv(cgo,c);
        CGOPickColor(cgo,*i,cPickableAtom);
        CGONormalv(cgo,tn);
        CGOVertexv(cgo,tv);
        tn+=3;
        tv+=3;
        CGONormalv(cgo,tn1);
        CGOVertexv(cgo,tv1);
        tn1+=3;
        tv1+=3;
        c+=3;
        i++;
      }
      tv+=3*I->N;
      tn+=3*I->N;
      tv1+=3*I->N;
      tn1+=3*I->N;
      CGOEnd(cgo);
    }

    if(SettingGet(I->G,cSetting_cartoon_debug)>1.5) {
      CGOEnable(cgo,GL_LIGHTING);
    }

    FreeP(TV);
    FreeP(TN);
  }
  
  PRINTFD(I->G,FB_Extrude)
    " ExtrudeCGOSurfacePolygonTaper-DEBUG: exiting...\n"
    ENDFD;
  
}





void ExtrudeCGOSurfaceStrand(CExtrude *I,CGO *cgo,int sampling,float *color_override)
{
  int a,b,*i;
  float *v;
  float *n;
  float *c;
  float *sv,*sn,*tv,*tn,*tv1,*tn1,*TV,*TN;
  float v0[3],n0[3],s0[3],z[3]={1.0,0.0,1.0};
  int subN;
  
  subN=I->N-sampling;

  PRINTFD(I->G,FB_Extrude)
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
          scale3f(s0,0.50F,s0);
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
      if(SettingGet(I->G,cSetting_cartoon_debug)<1.5)
        CGOBegin(cgo,GL_TRIANGLE_STRIP);
      else {
        CGOBegin(cgo,GL_LINE_STRIP);        
        CGODisable(cgo,GL_LIGHTING);
      }
      c = I->c;
      i = I->i;
      for(a=0;a<I->N;a++) {
        if(a<=subN) {
          if(color_override&&((b==2)||(b==3)||(b==6)||(b==7)))
            CGOColorv(cgo,color_override);
          else
            CGOColorv(cgo,c);
          CGOPickColor(cgo,*i,cPickableAtom);
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
        i++;
      }
      tv+=3*I->N;
      tn+=3*I->N;
      tv1+=3*I->N;
      tn1+=3*I->N;
      CGOEnd(cgo);
    }

    if(SettingGet(I->G,cSetting_cartoon_debug)>1.5) {
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
      CGOPickColor(cgo,I->i[0],cPickableAtom);
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
        s0[2]=s0[2]*((1.5F*((I->N-1)-a))/sampling);
        transform33Tf3f(n,s0,tv);
        add3f(v,tv,tv);
        tv+=3;
        copy3f(sn,n0);
        if(fabs(dot_product3f(sn,z))>R_SMALL4) {
          n0[0]+=0.4F;
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
      if(SettingGet(I->G,cSetting_cartoon_debug)<1.5)
        CGOBegin(cgo,GL_TRIANGLE_STRIP);
      else {
        CGOBegin(cgo,GL_LINE_STRIP);        
        CGODisable(cgo,GL_LIGHTING);
      }
      c = I->c;
      i = I->i;
      for(a=0;a<I->N;a++) {
        if(a>=(subN-1)) {
          if(color_override&&((b==2)||(b==3)||(b==6)||(b==7)))
            CGOColorv(cgo,color_override);
          else
            CGOColorv(cgo,c);
          CGOPickColor(cgo,*i,cPickableAtom);
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
        i++;
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
      s0[2]=s0[2]*1.5F;
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
    CGOPickColor(cgo,I->i[(subN-1)],cPickableAtom);
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
  
  PRINTFD(I->G,FB_Extrude)
    " ExtrudeCGOSurfaceStrand-DEBUG: exiting...\n"
    ENDFD;

}

void ExtrudeComputeScaleFactors(CExtrude *I,ObjectMolecule *obj,int source_field,
                                float mean, float stdev, float power, float range,
                                float min_scale, float max_scale,
                                int window)
{
  float *sf;
  int a;
  int *i;
  AtomInfoType *at;
  float scale;

  if(I->N&&I->Ns) {

    i=I->i;
    sf=I->sf;

    if(stdev>R_SMALL8) {
      for(a=0;a<I->N;a++) {
        at = obj->AtomInfo + (*i);
        
        switch(source_field) {
        default: /* b*/
          scale = (range+(at->b - mean)/stdev)/range;
          if(scale<0.0F)
            scale = 0.0F;
          scale=(float)pow(scale,power);
          if(scale<min_scale)
            scale=min_scale;
          if(scale>max_scale)
            scale=max_scale;
          *sf = scale;
          break;
        }
        sf++;
        i++;
      }
    } else {
      for(a=0;a<I->N;a++) {
        *sf = 1.0F;
        sf++;
      }
    }
   

    PRINTFB(I->G,FB_RepCartoon,FB_Blather)
            " Putty: mean %8.3f stdev %8.3f min %8.3f max %8.3f\n",
      mean,stdev,
      mean+(pow(min_scale,1.0F/power)*range-range)*stdev,
      mean+(pow(max_scale,1.0F/power)*range-range)*stdev
      ENDFB(I->G);
    /* now compute window average */
    
    
    {
      float *SF = Alloc(float,I->N);
      int w,ww;
      float accum;
      int cnt;

      sf=I->sf;
      
      for(a=1;a<(I->N-1);a++) {
        accum = 0.0F;
        cnt=0;
        for(w=-window;w<=window;w++) {
          ww=w + a;
          if(ww<0) ww=0;
          else if(ww>(I->N-1))
            ww=I->N-1;
          accum += sf[ww];
          cnt++;
        }
        SF[a] = accum / cnt;
      }
      for(a=1;a<I->N-1;a++)
        sf[a] = SF[a];
      FreeP(SF);
    }
  }
}

#if 0
/* PUTTY: look-up table of nearest neighbors and scaling factors (derived from nearby B-facts) */
void ExtrudeMakePuttyLUT(CExtrude *I,ObjectMolecule *Sauce_obj)
{
  int a,b;		  /* loop indices */
  float *v;    		  /* extrusion points */
  float *sf;		  /* PUTTY: scale factors */
  float min_dist = 0.0;   /* PUTTY: distance between given extrusion point 'p' and nearest atom */
  float minB = 0.0;	  /* PUTTY: overall minimum b-fact */
  float maxB = 0.0;   	  /* PUTTY: overall maximum b-fact */
  int my_i = 0;		  /* PUTTY: loop index */
  int startAt = 0;	  /* PUTTY: starting atom index for extrusion point nearest neighbor search */
  int endAt = 0;	  /* PUTTY: ending atom index for extrusion point nearest neighbor search */
  int atm2use = 0;	  /* PUTTY: temporary index of atom nearest to given extrusion point 'p' */
  
  PRINTFD(I->G,FB_Extrude)
    " ExtrudeMakePuttyLUT-DEBUG: entered.\n"
    ENDFD;

  if(I->N&&I->Ns) {
    
    float curAtmVrtx[3];
    float ex2atmDiff[3];

    /* Arbitrarily initialize to values for first atom and then extract max and min b-facts */

    minB = ( maxB = Sauce_obj->AtomInfo[0].b );

    for(my_i=0; my_i<((*(Sauce_obj->CSet))->NIndex); my_i++) { /* WLD: uses all atoms in coordinate set... */
      if (Sauce_obj->AtomInfo[my_i].b > maxB)	
        maxB = Sauce_obj->AtomInfo[my_i].b;
      else if (Sauce_obj->AtomInfo[my_i].b < minB)	
        minB = Sauce_obj->AtomInfo[my_i].b;
    }

    PRINTFB(I->G,FB_RepCartoon,FB_Blather)
      "Puttyage: will scale between minB = %.2f and maxB = %.2f for this chain...\n",minB,maxB
      ENDFB(I->G)
    
    for(b=0;b<=I->Ns;b++) {

      v=I->p;
      sf=I->sf;
      
      for(a=0;a<I->N;a++) {
        
        /* NOTE: To fully scan each extrusion point against all atoms, set startAt = 0 and 
         * 	 endAt = ((*(Sauce_obj->CSet))->NIndex). This is VERY slow if have large num 
         * 	 atoms (>500), so I've assumed that atoms of interest for given extrusion pnt
         * 	 will be within +/- 30 of nearest atom to previous extrusion point. Also, could
         * 	 consider just the CA's, but that makes more jagged cartoons... */

        startAt = ( *(z-3) && *(z-3) > 30 ) ? (*(z-3) - 30) : 0 ;
        endAt =   ( *(z-3) && *(z-3) < ((*(Sauce_obj->CSet))->NIndex - 30) ) ? (*(z-3) + 30) : ((*(Sauce_obj->CSet))->NIndex) ;

        /* Arbitrarily initialize atm2use to startAt and min_dist to distance between extrusion 
         * point 'v' and this atom: */

        CoordSetGetAtomVertex((*(Sauce_obj->CSet)),startAt,curAtmVrtx);
        subtract3f(curAtmVrtx,v,ex2atmDiff);
        min_dist = length3f(ex2atmDiff);
        atm2use = startAt;
        
        for(my_i=startAt; my_i<endAt; my_i++) {
          /*			if(strcmp(Sauce_obj->AtomInfo[my_i].name,"CA") == 0) {		*/
          CoordSetGetAtomVertex((*(Sauce_obj->CSet)),my_i,curAtmVrtx);

          /* basic atom info: name, resid #, b-fact */
          /* printf("name: '%s', resi: %d, b = %f)\n",	Sauce_obj->AtomInfo[my_i].name,
           *		   	   			Sauce_obj->AtomInfo[my_i].resi,
           *		   	   			Sauce_obj->AtomInfo[my_i].b);  		*/

          subtract3f(curAtmVrtx,v,ex2atmDiff);	  

          /*		   	   printf("(%.3f, %.3f, %.3f)\n",*ex2atmDiff,*(ex2atmDiff+1),*(ex2atmDiff+2));		*/

          if ( length3f(ex2atmDiff) <= min_dist )  {
            min_dist = length3f(ex2atmDiff);
            atm2use = my_i; 
          }

          /*			}		*/
        }
        
        /* Apply an arbitrarily chosen nonlinear scaling; could try other recipes, exponents, etc.
         * Note: doing this as triplets b/c could imagine utilizing anisotropic scaling factors or
         * something like that in future... */

        *(sf) =   (float)pow((2.0*Sauce_obj->AtomInfo[atm2use].b / (maxB - minB)), 1.5 );
        *(sf+1) = (float)pow((2.0*Sauce_obj->AtomInfo[atm2use].b / (maxB - minB)), 1.5 );
        *(sf+2) = (float)pow((2.0*Sauce_obj->AtomInfo[atm2use].b / (maxB - minB)), 1.5 );
        

        v+=3;
        sf+=3;
      }
    }
    
  }
  
  PRINTFD(I->G,FB_Extrude)
    " ExtrudeMakePuttyLUT-DEBUG: exiting...\n"
    ENDFD;
}
#endif

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
    FreeP(I->i);
    FreeP(I->sf);		/* PUTTY */
    I->p = Alloc(float,3*(n+1));
    I->n = Alloc(float,9*(n+1));
    I->c = Alloc(float,3*(n+1));
    I->i = Alloc(int,3*(n+1));
    I->sf= Alloc(float,n+1); /* PUTTY: scale factors */ 
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
  FreeP(I->i);
  FreeP(I->sf);
  OOFreeP(I);
}


