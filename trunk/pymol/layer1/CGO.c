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

#include"CGO.h"
#include"Base.h"
#include"GL/gl.h"
#include"OOMac.h"
#include"Setting.h"
#include"Sphere.h"

#define CGO_read_int(p) (*((int*)(p++)))
#define CGO_get_int(p) (*((int*)(p)))
#define CGO_write_int(p,i) ((*((int*)(p++)))=(i))


int CGO_sz[] = {
  CGO_NULL_SZ,
  CGO_NULL_SZ,
  CGO_BEGIN_SZ,
  CGO_END_SZ,
  CGO_VERTEX_SZ,
  
  CGO_NORMAL_SZ,
  CGO_COLOR_SZ,
  CGO_SPHERE_SZ,
  CGO_TRIANGLE_SZ,
  CGO_CYLINDER_SZ,
  
  CGO_LINEWIDTH_SZ,
  CGO_WIDTHSCALE_SZ,
  CGO_NULL_SZ,
  CGO_NULL_SZ,
  CGO_NULL_SZ,

  CGO_NULL_SZ
};

typedef void CGO_op(float *);
typedef CGO_op *CGO_op_fn;

static float *CGO_add(CGO *I,int c);
static float *CGO_size(CGO *I,int sz);
static void subdivide( int n, float *x, float *y);
void CGOSimpleCylinder(CGO *I,float *v1,float *v2,float tube_size,float *c1,float *c2);
void CGOSimpleSphere(CGO *I,float *v,float vdw);

CGO *CGONew(void)
{
  OOAlloc(CGO);
  I->op=VLAlloc(float,33);
  I->c=0;
  return(I);
}

CGO *CGONewSized(int size)
{
  OOAlloc(CGO);
  I->op=VLAlloc(float,size+32);
  I->c=0;
  return(I);
}

void CGOFree(CGO *I)
{
  if(I) {
    VLAFreeP(I->op);
  }
  OOFreeP(I);
}

static float *CGO_add(CGO *I,int c)
{
  float *at;
  VLACheck(I->op,float,I->c+c);
  at = I->op+I->c;
  I->c+=c;
  return(at);
}

static float *CGO_size(CGO *I,int sz)
{
  float *at;
  VLASize(I->op,float,sz);
  at=I->op+I->c;
  I->c=sz;
  return(at);
}

/*===== Object Creation Routines =======*/

int CGOFromFloatArray(CGO *I,float *src,int len)
{
  int op;
  int c;
  int ok;
  int all_ok=true;
  int bad_entry=0;
  int sz;
  int a;
  int cc;
  float val;
  float *pc,*save_pc,*tf;
  VLACheck(I->op,float,I->c+len);
  save_pc=I->op+I->c;
  while(len-->0) {
    cc++;
    c=1;
    op = CGO_MASK&((int)(*(src++)));
    sz = CGO_sz[op];
    if(len<sz) 
      break; /* discard short instruction */
    len-=sz;
    pc=save_pc;
    CGO_write_int(pc,op);
    ok=true;
    for(a=0;a<sz;a++) {
      cc++;
      val=*(src++);
      if((FLT_MAX-val)>0.0F) { /* make sure we have a real float */
        *(pc++)=val;
      } else {
        *(pc++)=0.0;
        ok=false;
      }
    }
    if(ok) {
      switch(op) { /* now convert any instructions with int arguments */
      case CGO_BEGIN:
        tf=save_pc+1;
        CGO_write_int(tf,*(tf));
        break;
      }
      save_pc=pc;
      I->c+=sz+1;
    } else {  /* discard illegal instructions */
      if(all_ok) 
        bad_entry = cc;
      all_ok=false;
    }
  }
  return(bad_entry);
}


void CGOBegin(CGO *I,int mode)
{
  float *pc = CGO_add(I,2);
  CGO_write_int(pc,CGO_BEGIN);
  CGO_write_int(pc,mode);
}

void CGOVertex(CGO *I,float v1,float v2,float v3)
{
  float *pc = CGO_add(I,4);
  CGO_write_int(pc,CGO_VERTEX);
  *(pc++)=v1;
  *(pc++)=v2;
  *(pc++)=v3;
}

void CGOVertexv(CGO *I,float *v)
{
  float *pc = CGO_add(I,4);
  CGO_write_int(pc,CGO_VERTEX);
  *(pc++)=*(v++);
  *(pc++)=*(v++);
  *(pc++)=*(v++);
}

void CGOColor(CGO *I,float v1,float v2,float v3)
{
  float *pc = CGO_add(I,4);
  CGO_write_int(pc,CGO_COLOR);
  *(pc++)=v1;
  *(pc++)=v2;
  *(pc++)=v3;
}

void CGOColorv(CGO *I,float *v)
{
  float *pc = CGO_add(I,4);
  CGO_write_int(pc,CGO_COLOR);
  *(pc++)=*(v++);
  *(pc++)=*(v++);
  *(pc++)=*(v++);
}

void CGONormal(CGO *I,float v1,float v2,float v3)
{
  float *pc = CGO_add(I,4);
  CGO_write_int(pc,CGO_NORMAL);
  *(pc++)=v1;
  *(pc++)=v2;
  *(pc++)=v3;
}

void CGONormalv(CGO *I,float *v)
{
  float *pc = CGO_add(I,4);
  CGO_write_int(pc,CGO_NORMAL);
  *(pc++)=*(v++);
  *(pc++)=*(v++);
  *(pc++)=*(v++);
}

void CGOEnd(CGO *I)
{
  float *pc = CGO_add(I,1);
  CGO_write_int(pc,CGO_END);
}

void CGOStop(CGO *I)
{
  /* add enough zeros to prevent overrun in the event of corruption
   * (include more zeros than the longest instruction in the compiler */

  float *pc = CGO_size(I,I->c+32); 

  CGO_write_int(pc,CGO_STOP);
  CGO_write_int(pc,CGO_STOP);
  CGO_write_int(pc,CGO_STOP);
  CGO_write_int(pc,CGO_STOP);
  CGO_write_int(pc,CGO_STOP);
  CGO_write_int(pc,CGO_STOP);
  CGO_write_int(pc,CGO_STOP);
  CGO_write_int(pc,CGO_STOP);

  CGO_write_int(pc,CGO_STOP);
  CGO_write_int(pc,CGO_STOP);
  CGO_write_int(pc,CGO_STOP);
  CGO_write_int(pc,CGO_STOP);
  CGO_write_int(pc,CGO_STOP);
  CGO_write_int(pc,CGO_STOP);
  CGO_write_int(pc,CGO_STOP);
  CGO_write_int(pc,CGO_STOP);

  CGO_write_int(pc,CGO_STOP);
  CGO_write_int(pc,CGO_STOP);
  CGO_write_int(pc,CGO_STOP);
  CGO_write_int(pc,CGO_STOP);
  CGO_write_int(pc,CGO_STOP);
  CGO_write_int(pc,CGO_STOP);
  CGO_write_int(pc,CGO_STOP);
  CGO_write_int(pc,CGO_STOP);
 
  CGO_write_int(pc,CGO_STOP);
  CGO_write_int(pc,CGO_STOP);
  CGO_write_int(pc,CGO_STOP);
  CGO_write_int(pc,CGO_STOP);
  CGO_write_int(pc,CGO_STOP);
  CGO_write_int(pc,CGO_STOP);

  CGO_write_int(pc,CGO_STOP);
  CGO_write_int(pc,CGO_STOP);

}

int CGOCheckComplex(CGO *I)
{
  register float *pc = I->op;
  int fc=0;
  int nEdge;
  int op;
  SphereRec *sp;
  
  sp = Sphere1;

  nEdge= SettingGet(cSetting_stick_quality);

  while((op=(CGO_MASK&CGO_read_int(pc)))) {
    switch(op) {
    case CGO_CYLINDER:
      fc+=3*(3+(nEdge+1)*9)+9;
      break;
    case CGO_SPHERE:
      fc+=(sp->NVertTot*6)+(sp->NStrip*3)+3;
      break;
    }
    pc+=CGO_sz[op];
  }
  return(fc);
}

CGO *CGOSimplify(CGO *I,int est)
{
  CGO *cgo;

  register float *pc = I->op;
  register float *nc;
  register int op;
  float *save_pc;
  int sz;
  
  cgo=CGONewSized(I->c+est);

  while((op=(CGO_MASK&CGO_read_int(pc)))) {
    save_pc=pc;
    switch(op) {
    case CGO_CYLINDER:
      CGOSimpleCylinder(cgo,pc,pc+3,*(pc+6),pc+7,pc+10);
      break;
    case CGO_SPHERE:
      CGOSimpleSphere(cgo,pc,*(pc+3));
      break;
    default:
      sz=CGO_sz[op];
      nc=CGO_add(cgo,sz+1);
      *(nc++)=*(pc-1);
      while(sz--)
        *(nc++)=*(pc++);
    }
    pc=save_pc;
    pc+=CGO_sz[op];
  }
  CGOStop(cgo);
  return(cgo);
}

/* ======== Raytrace Renderer ======== */

void CGORenderRay(CGO *I,CRay *ray) 
{
  register float *pc = I->op;
  register int op;
  int vc;
  float linewidth=1.0;
  float widthscale=0.15;
  float primwidth=0.15;
  float white[] = {1.0,1.0,1.0};
  float zee[] = {0.0,0.0,1.0};

  float *n0,*n1,*n2,*v0,*v1,*v2,*c0,*c1,*c2;
  int mode = -1;
  c0=white;

  while((op=(CGO_MASK&CGO_read_int(pc)))) {
    switch(op) {
    case CGO_BEGIN:
      mode=CGO_get_int(pc);
      vc=0;
      n0=zee;
      break;
    case CGO_END:
      switch(mode) {
      case GL_LINE_LOOP:
        if(vc>1) ray->fCylinder3fv(ray,v0,v2,primwidth,c0,c2);
        break;
      }
      mode=-1;
      break;
    case CGO_WIDTHSCALE:
      widthscale=*pc;
      primwidth=widthscale*linewidth;
      break;
    case CGO_LINEWIDTH:
      linewidth=*pc;
      primwidth=widthscale*linewidth;
      break;
    case CGO_NORMAL:
      n0=pc;
      break;
    case CGO_COLOR:
      c0=pc;
      ray->fColor3fv(ray,c0);
      break;
    case CGO_VERTEX:
      v0=pc;
      switch(mode) {
      case GL_POINTS:
        ray->fSphere3fv(ray,v0,primwidth);
        break;
      case GL_LINES:
        if(vc&0x1) ray->fCylinder3fv(ray,v0,v1,primwidth,c0,c1);
        v1=v0;
        c1=c0;
        break;
      case GL_LINE_STRIP:
        if(vc) ray->fCylinder3fv(ray,v0,v1,primwidth,c0,c1);
        v1=v0;
        c1=c0;
        break;
      case GL_LINE_LOOP:
        if(vc) 
          ray->fCylinder3fv(ray,v0,v1,primwidth,c0,c1);
        else {
          v2=v0;
          c2=c0;
        }
        v1=v0;
        c1=c0;
        break;
      case GL_TRIANGLES:
        if(3*((vc+1)/3)==vc+1)
          ray->fTriangle3fv(ray,v0,v1,v2,n0,n1,n2,c0,c1,c2);
        v2=v1;
        c2=c1;
        n2=n1;
        v1=v0;
        c1=c0;
        n1=n0;
        break;
      case GL_TRIANGLE_STRIP:
        if(vc>1)
          ray->fTriangle3fv(ray,v0,v1,v2,n0,n1,n2,c0,c1,c2);
        v2=v1;
        c2=c1;
        n2=n1;
        v1=v0;
        c1=c0;
        n1=n0;
        break;
      case GL_TRIANGLE_FAN:
        if(vc>1)
          ray->fTriangle3fv(ray,v0,v1,v2,n0,n1,n2,c0,c1,c2);
        else if(!vc) {
          v2=v0;
          c2=c0;
        }
        v1=v0;
        c1=c0;
        n1=n0;
        break;
      }
      vc++;
      break;
    case CGO_SPHERE:
      ray->fColor3fv(ray,c0);
      ray->fSphere3fv(ray,pc,*(pc+3));
      break;
    case CGO_CYLINDER:
      ray->fCylinder3fv(ray,pc,pc+3,*(pc+6),pc+7,pc+10);
      break;
    case CGO_TRIANGLE:
      ray->fTriangle3fv(ray,pc,pc+3,pc+6,pc+9,pc+12,pc+15,pc+18,pc+21,pc+24);
      break;
    }
    pc+=CGO_sz[op];
  }
}

/* ======== GL Rendering ======== */

static void CGO_gl_begin(float *pc)
{
  glBegin(CGO_read_int(pc));
}

static void CGO_gl_end(float *pc)
{
  glEnd();
}

static void CGO_gl_linewidth(float *pc)
{
  glLineWidth(*pc);
}

static void CGO_gl_null(float *pc) {
}
    
/* dispatch table for OpenGL */

CGO_op_fn CGO_gl[] = {
  CGO_gl_null,
  CGO_gl_null,
  CGO_gl_begin,
  CGO_gl_end,
  (CGO_op_fn)glVertex3fv,
  
  (CGO_op_fn)glNormal3fv,
  (CGO_op_fn)glColor3fv,
  CGO_gl_null,
  CGO_gl_null,
  CGO_gl_null,
  
  CGO_gl_linewidth,
  CGO_gl_null,
  CGO_gl_null,
  CGO_gl_null,
  CGO_gl_null,
  
  CGO_gl_null,
};

void CGORenderGL(CGO *I) /* this should be as fast as you can make it...
                          * the ASM loop is about 2X long as raw looped GL calls,
                          * but hopefully superscaler processors won't care */
{
  register float *pc = I->op;
  register int op;

  while((op=(CGO_MASK&CGO_read_int(pc)))) {
    CGO_gl[op](pc);
    pc+=CGO_sz[op];
  }
}

/* translation function which turns cylinders and spheres into triangles */

void CGOSimpleSphere(CGO *I,float *v,float vdw)
{
  SphereRec *sp;
  int *q,*s;
  int b,c;

  sp = Sphere1;
  
  q=sp->Sequence;

  s=sp->StripLen;

  for(b=0;b<sp->NStrip;b++)
    {
      CGOBegin(I,GL_TRIANGLE_STRIP);
      for(c=0;c<(*s);c++)
        {
          CGONormalv(I,sp->dot[*q].v);
          CGOVertex(I,v[0]+vdw*sp->dot[*q].v[0],
                    v[1]+vdw*sp->dot[*q].v[1],
                    v[2]+vdw*sp->dot[*q].v[2]);
          q++;
        }
      CGOEnd(I);
      s++;
    }
}

static void subdivide( int n, float *x, float *y)
{
  int a;
  if(n<3) {n=3;}
  for(a=0;a<=n;a++)
	 {
		x[a]=cos(a*2*PI/n);
		y[a]=sin(a*2*PI/n);
	 }
}

void CGOSimpleCylinder(CGO *I,float *v1,float *v2,float tube_size,float *c1,float *c2)
{

#define MAX_EDGE 50

  float d[3],t[3],p0[3],p1[3],p2[3],v_buf[9],*v;
  float x[50],y[50];
  float overlap;
  float nub;
  int colorFlag;
  int nEdge;
  int c;

  v=v_buf;
  nEdge= SettingGet(cSetting_stick_quality);
  overlap = tube_size*SettingGet(cSetting_stick_overlap);
  nub = tube_size*SettingGet(cSetting_stick_nub);

  if(nEdge>MAX_EDGE)
    nEdge=MAX_EDGE;
  subdivide(nEdge,x,y);

  colorFlag=(c1!=c2)&&c2;

  CGOColorv(I,c1);

  /* direction vector */
  
  p0[0] = (v2[0] - v1[0]);
  p0[1] = (v2[1] - v1[1]);
  p0[2] = (v2[2] - v1[2]);
  
  normalize3f(p0);
  
  v1[0]-=p0[0]*overlap;
  v1[1]-=p0[1]*overlap;
  v1[2]-=p0[2]*overlap;
  
  v2[0]+=p0[0]*overlap;
  v2[1]+=p0[1]*overlap;
  v2[2]+=p0[2]*overlap;
  
  d[0] = (v2[0] - v1[0]);
  d[1] = (v2[1] - v1[1]);
  d[2] = (v2[2] - v1[2]);
  
  t[0] = d[1];
  t[1] = d[2];
  t[2] = -d[0];
  
  cross_product3f(d,t,p1);
  
  normalize3f(p1);
  
  cross_product3f(d,p1,p2);
  
  normalize3f(p2);
  
  /* now we have a coordinate system*/
  
  CGOBegin(I,GL_TRIANGLE_STRIP);
  for(c=0;c<=nEdge;c++)
	 {
		v[0] = p1[0]*x[c] + p2[0]*y[c];
		v[1] = p1[1]*x[c] + p2[1]*y[c];
		v[2] = p1[2]*x[c] + p2[2]*y[c];
		
		v[3] = v1[0] + v[0]*tube_size;
		v[4] = v1[1] + v[1]*tube_size;
		v[5] = v1[2] + v[2]*tube_size;
		
		v[6] = v[3] + d[0];
		v[7] = v[4] + d[1];
		v[8] = v[5] + d[2];
		
      CGONormalv(I,v);
      if(colorFlag) CGOColorv(I,c1);
      CGOVertexv(I,v+3);
      if(colorFlag) CGOColorv(I,c2);
      CGOVertexv(I,v+6);
	 }
  CGOEnd(I);

  v[0] = -p0[0];
  v[1] = -p0[1];
  v[2] = -p0[2];
  
  v[3] = v1[0] - p0[0]*nub;
  v[4] = v1[1] - p0[1]*nub;
  v[5] = v1[2] - p0[2]*nub;
  
  if(colorFlag) CGOColorv(I,c1);
  CGOBegin(I,GL_TRIANGLE_FAN);
  CGONormalv(I,v);
  CGOVertexv(I,v+3);

  for(c=0;c<=nEdge;c++)
	 {
		
		v[0] = p1[0]*x[c] + p2[0]*y[c];
		v[1] = p1[1]*x[c] + p2[1]*y[c];
		v[2] = p1[2]*x[c] + p2[2]*y[c];
		
		v[3] = v1[0] + v[0]*tube_size;
		v[4] = v1[1] + v[1]*tube_size;
		v[5] = v1[2] + v[2]*tube_size;
		
      CGONormalv(I,v);
      CGOVertexv(I,v+3);
	 }
  CGOEnd(I);

  v[0] = p0[0];
  v[1] = p0[1];
  v[2] = p0[2];
  
  v[3] = v2[0] + p0[0]*nub;
  v[4] = v2[1] + p0[1]*nub;
  v[5] = v2[2] + p0[2]*nub;
  
  if(colorFlag) CGOColorv(I,c2);
  CGOBegin(I,GL_TRIANGLE_FAN);
  CGONormalv(I,v);
  CGOVertexv(I,v+3);

  for(c=0;c<=nEdge;c++)
    {
      
      v[0] = p1[0]*x[c] + p2[0]*y[c];
      v[1] = p1[1]*x[c] + p2[1]*y[c];
      v[2] = p1[2]*x[c] + p2[2]*y[c];
      
      v[3] = v2[0] + v[0]*tube_size;
      v[4] = v2[1] + v[1]*tube_size;
      v[5] = v2[2] + v[2]*tube_size;
      
      CGONormalv(I,v);
      CGOVertexv(I,v+3);
    }
  CGOEnd(I);
}


