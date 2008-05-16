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
#include"os_predef.h"
#include"os_std.h"

#include"Base.h"
#include"Vector.h"
#include"Matrix.h"

#ifndef R_SMALL
#define R_SMALL 0.000000001
#endif

#ifndef R_MED
#define R_MED 0.00001
#endif

#ifndef true
#define true 1
#endif

#ifndef false
#define false 0
#endif


#define cPI            3.14159265358979323846  /* pi */

static const float _0		= 0.0F;
static const float _1		= 1.0F;
static const double _d0		= 0.0;
static const double _d1    = 1.0;

void mix3f(float *v1,float *v2,float fxn,float *v3)
{
  float fxn_1 = 1.0F - fxn;
  v3[0] = v1[0] * fxn_1 + v2[0] * fxn;
  v3[1] = v1[1] * fxn_1 + v2[1] * fxn;
  v3[2] = v1[2] * fxn_1 + v2[2] * fxn;
}

void mix3d(double *v1,double *v2,double fxn,double *v3)
{
  double fxn_1 = 1.0F - fxn;
  v3[0] = v1[0] * fxn_1 + v2[0] * fxn;
  v3[1] = v1[1] * fxn_1 + v2[1] * fxn;
  v3[2] = v1[2] * fxn_1 + v2[2] * fxn;
}

unsigned int optimizer_workaround1u(unsigned int value)
{
  return value;
}

float get_random0to1f()
{
  return(rand()/(_1+RAND_MAX));
}

int pymol_roundf(float f)
{
  if(f>0.0F)
    return (int)(f+0.49999F);
  else
    return (int)(f-0.49999F);
}

double slow_sqrt1f(float f) { /* no good as a macro because f is used twice */
  if(f>_0)
	 return(sqrt(f));
  else
	 return _d0;
}

double slow_sqrt1d(double f) { /* no good as a macro because f is used twice */
  if(f>_0)
	 return(sqrt(f));
  else
	 return _d0;
}
void dump3i( int *v, char *prefix ) /* for debugging */
{
  printf("%s %8i %8i %8i\n",prefix,v[0],v[1],v[2]);
}

void dump3f( float *v, char *prefix ) /* for debugging */
{
  printf("%s %8.3f %8.3f %8.3f\n",prefix,v[0],v[1],v[2]);
}

void dump3d( double *v, char *prefix ) /* for debugging */
{
  printf("%s %8.3f %8.3f %8.3f\n",prefix,v[0],v[1],v[2]);
}

void dump4f( float *v, char *prefix ) /* for debugging */
{
  printf("%s %8.3f %8.3f %8.3f %8.3f\n",prefix,v[0],v[1],v[2],v[3]);
}

void dump33f( float *m, char *prefix ) /* for debugging */
{
  if(m) {
    printf("%s:0 %8.3f %8.3f %8.3f\n",prefix,m[0],m[1],m[2]);
    printf("%s:1 %8.3f %8.3f %8.3f\n",prefix,m[3],m[4],m[5]);
    printf("%s:2 %8.3f %8.3f %8.3f\n",prefix,m[6],m[7],m[8]);
  } else {
    printf("%s: (null matrix pointer)\n",prefix);
  }
}

void dump44f( float *m, char *prefix ) /* for debugging */
{
  if(m) {
    if(prefix) {
      printf("%s:0 %8.3f %8.3f %8.3f %8.3f\n",prefix,m[0],m[1],m[2],m[3]);
      printf("%s:1 %8.3f %8.3f %8.3f %8.3f\n",prefix,m[4],m[5],m[6],m[7]);
      printf("%s:2 %8.3f %8.3f %8.3f %8.3f\n",prefix,m[8],m[9],m[10],m[11]);
      printf("%s:3 %8.3f %8.3f %8.3f %8.3f\n",prefix,m[12],m[13],m[14],m[15]);
    } else {
    }
  } else {
    printf("%s: (null matrix pointer)\n",prefix);
  }

}

void dump44d( double *m, char *prefix ) /* for debugging */
{
  if(m) {
    printf("%s:0 %8.3f %8.3f %8.3f %8.3f\n",prefix,m[0],m[1],m[2],m[3]);
    printf("%s:1 %8.3f %8.3f %8.3f %8.3f\n",prefix,m[4],m[5],m[6],m[7]);
    printf("%s:2 %8.3f %8.3f %8.3f %8.3f\n",prefix,m[8],m[9],m[10],m[11]);
    printf("%s:3 %8.3f %8.3f %8.3f %8.3f\n",prefix,m[12],m[13],m[14],m[15]);
  } else {
    printf("%s: (null matrix pointer)\n",prefix);
  }
}

void dump33d( double *m, char *prefix ) /* for debugging */
{
  printf("%s:0 %8.3f %8.3f %8.3f\n",prefix,m[0],m[1],m[2]);
  printf("%s:1 %8.3f %8.3f %8.3f\n",prefix,m[3],m[4],m[5]);
  printf("%s:2 %8.3f %8.3f %8.3f\n",prefix,m[6],m[7],m[8]);
}

void get_divergent3f(float *src,float *dst)
{
  if(src[0]!=_0) {
    *(dst++)=-*(src++);
    *(dst++)= *(src++)+0.1F;
    *(dst++)= *(src++);
  } else if(src[1]!=_0) {
    *(dst++)= *(src++)+0.1F;
    *(dst++)=-*(src++);
    *(dst++)= *(src++);
  } else {
    *(dst++)= *(src++)+0.1F;
    *(dst++)= *(src++);
    *(dst++)=-*(src++);
  }
}

int equal3f(float *v1,float *v2)
{
  return((fabs(v1[0]-v2[0])<R_SMALL)&&
         (fabs(v1[1]-v2[1])<R_SMALL)&&
         (fabs(v1[2]-v2[2])<R_SMALL));
}

void get_random3f(float *x) /* this needs to be fixed as in Tinker */
{
  x[0]=0.5F-(rand()/(_1+RAND_MAX));
  x[1]=0.5F-(rand()/(_1+RAND_MAX));
  x[2]=0.5F-(rand()/(_1+RAND_MAX));
  normalize3f(x);
}

void get_system3f(float *x,float *y,float *z) /* make random system */
{
  get_random3f(x);
  get_divergent3f(x,y);
  cross_product3f(x,y,z);
  normalize3f(z);
  cross_product3f(z,x,y);
  normalize3f(y);
  normalize3f(x);
}

void get_system1f3f(float *x,float *y,float *z) /* make system in direction of x */
{
  get_divergent3f(x,y);
  cross_product3f(x,y,z);
  normalize3f(z);
  cross_product3f(z,x,y);
  normalize3f(y);
  normalize3f(x);
}

void get_system2f3f(float *x,float *y,float *z) /* make system in direction of x */
{
  cross_product3f(x,y,z);
  normalize3f(z);
  cross_product3f(z,x,y);
  normalize3f(y);
  normalize3f(x);
}

void extrapolate3f(float *v1, float *unit, float *result)
{
  float lsq = lengthsq3f(v1);
  float dp = dot_product3f(v1,unit);
  if(dp!=0.0F) {
    float l2 = lsq/dp;
    scale3f(unit,l2,result);
  }
}


void scatter3f(float *v,float weight)
{
  float r[3];
  get_random3f(r);
  scale3f(r,weight,r);
  add3f(r,v,v);
  normalize3f(v);
}

void wiggle3f(float *v,float *p,float *s)
{
  float q[3];
  q[0]=(float)cos((p[0]+p[1]+p[2])*s[1]);
  q[1]=(float)cos((p[0]-p[1]+p[2])*s[1]);
  q[2]=(float)cos((p[0]+p[1]-p[2])*s[1]);
  scale3f(q,s[0],q);
  add3f(q,v,v);
  normalize3f(v);

}
void copy3d3f ( double *v1,float *v2)
{
  v2[0]=(float)v1[0];
  v2[1]=(float)v1[1];
  v2[2]=(float)v1[2];
}
void copy3f3d ( float *v1,double *v2)
{
  v2[0]=(double)v1[0];
  v2[1]=(double)v1[1];
  v2[2]=(double)v1[2];
}


#ifndef USE_VECTOR_MACROS

void  set3f ( float *v1,float x,float y,float z )
{
  v1[0]=x;
  v1[1]=y;
  v1[2]=z;
}

void swap1f(float *f,float *g)
{
  float h;
  h=*f;
  *f=*g;
  *g=h;
}

void zero3f(float *v1)
{
	v1[0]=_0;
	v1[1]=_0;
	v1[2]=_0;
}

void ones3f(float *v1)
{
  v1[0]=_1;
  v1[1]=_1;
  v1[2]=_1;
}

float dot_product3f ( float *v1, float *v2 )
{
  return( v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2]);
}

void copy4f(v1,v2) 
{(v2)[0]=(v1)[0]; (v2)[1]=(v1)[1]; (v2)[2]=(v1)[2]; (v2)[3]=(v1)[3];}

void invert3f (float *v)
{
  v[0]=-v[0];
  v[1]=-v[1];
  v[2]=-v[2];
}

void invert3f3f (float *v1,float *v2)
{
  v2[0]=-v1[0];
  v2[1]=-v1[1];
  v2[2]=-v1[2];
}

void scale3f ( float *v1,float v0,float *v2)
{
  v2[0]=v1[0]*v0;
  v2[1]=v1[1]*v0;
  v2[2]=v1[2]*v0;
}

void copy3f ( float *v1,float *v2)
{
  v2[0]=v1[0];
  v2[1]=v1[1];
  v2[2]=v1[2];
}

void copy3d ( double *v1,double *v2)
{
  v2[0]=v1[0];
  v2[1]=v1[1];
  v2[2]=v1[2];
}



void add3f ( float *v1, float *v2, float *v3 )
{
  v3[0]=v1[0]+v2[0];
  v3[1]=v1[1]+v2[1];
  v3[2]=v1[2]+v2[2];
}

void subtract3f ( float *v1, float *v2, float *v3 )
{
  v3[0]=v1[0]-v2[0];
  v3[1]=v1[1]-v2[1];
  v3[2]=v1[2]-v2[2];
}

void cross_product3f ( float *v1, float *v2, float *cross )
{
  cross[0] = (v1[1]*v2[2]) - (v1[2]*v2[1]);
  cross[1] = (v1[2]*v2[0]) - (v1[0]*v2[2]);
  cross[2] = (v1[0]*v2[1]) - (v1[1]*v2[0]);
}

double lengthsq3f ( float *v1 )
{
  return((v1[0]*v1[0]) + 
         (v1[1]*v1[1]) + 
         (v1[2]*v1[2]));
} 

double length3f ( float *v1 )
{
  return(sqrt1d((v1[0]*v1[0]) + 
					 (v1[1]*v1[1]) + 
					 (v1[2]*v1[2])));
} 

void average3f ( float *v1, float *v2, float *avg )
{
  avg[0] = (v1[0]+v2[0])/2.0F;
  avg[1] = (v1[1]+v2[1])/2.0F;
  avg[2] = (v1[2]+v2[2])/2.0F;
}

#endif


void min3f ( float *v1, float *v2, float *v3 )
{
  (v3)[0] = ( (v1)[0]<(v2)[0] ? (v1)[0] : (v2)[0] ); 
  (v3)[1] = ( (v1)[1]<(v2)[1] ? (v1)[1] : (v2)[1] ); 
  (v3)[2] = ( (v1)[2]<(v2)[2] ? (v1)[2] : (v2)[2] ); 
}

void max3f ( float *v1, float *v2, float *v3 )
{
  (v3)[0] = ( (v1)[0]>(v2)[0] ? (v1)[0] : (v2)[0] ); 
  (v3)[1] = ( (v1)[1]>(v2)[1] ? (v1)[1] : (v2)[1] ); 
  (v3)[2] = ( (v1)[2]>(v2)[2] ? (v1)[2] : (v2)[2] ); 
}


float slow_project3f ( float *v1, float *v2, float *proj )
{
   float dot;

	dot = v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2];
	proj[0] = v2[0] * dot;
	proj[1] = v2[1] * dot;
	proj[2] = v2[2] * dot;
	
	return(dot);
}

double dot_product3d ( double *v1, double *v2 )
{
  return( v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2]);
}

void identity33f( float *m)
{
  m[0]=_1;
  m[1]=_0;
  m[2]=_0;
  m[3]=_0;
  m[4]=_1;
  m[5]=_0;
  m[6]=_0;
  m[7]=_0;
  m[8]=_1;
}

void identity33d( double *m)
{
  m[0]=_d1;
  m[1]=_d0;
  m[2]=_d0;
  m[3]=_d0;
  m[4]=_d1;
  m[5]=_d0;
  m[6]=_d0;
  m[7]=_d0;
  m[8]=_d1;
}

void identity44f ( float *m1 )
{
  int a;
  for(a=0;a<16;a++) m1[a]=_0;
  for(a=0;a<16;a=a+5) m1[a]=_1;
}
void identity44d ( double *m1 )
{
  int a;
  for(a=0;a<16;a++) m1[a]=_d0;
  for(a=0;a<16;a=a+5) m1[a]=_d1;
}


void copy44f ( float *src, float *dst )
{
  *(dst++)=*(src++);
  *(dst++)=*(src++);
  *(dst++)=*(src++);
  *(dst++)=*(src++);

  *(dst++)=*(src++);
  *(dst++)=*(src++);
  *(dst++)=*(src++);
  *(dst++)=*(src++);

  *(dst++)=*(src++);
  *(dst++)=*(src++);
  *(dst++)=*(src++);
  *(dst++)=*(src++);

  *(dst++)=*(src++);
  *(dst++)=*(src++);
  *(dst++)=*(src++);
  *(dst++)=*(src++);
}

void copy44d ( double *src, double *dst )
{
  *(dst++)=*(src++);
  *(dst++)=*(src++);
  *(dst++)=*(src++);
  *(dst++)=*(src++);

  *(dst++)=*(src++);
  *(dst++)=*(src++);
  *(dst++)=*(src++);
  *(dst++)=*(src++);

  *(dst++)=*(src++);
  *(dst++)=*(src++);
  *(dst++)=*(src++);
  *(dst++)=*(src++);

  *(dst++)=*(src++);
  *(dst++)=*(src++);
  *(dst++)=*(src++);
  *(dst++)=*(src++);
}

void copy44d33f ( double *src, float *dst )
{
  *(dst++)=(float)*(src++);
  *(dst++)=(float)*(src++);
  *(dst++)=(float)*(src++);
  src++;

  *(dst++)=(float)*(src++);
  *(dst++)=(float)*(src++);
  *(dst++)=(float)*(src++);
  src++;

  *(dst++)=(float)*(src++);
  *(dst++)=(float)*(src++);
  *(dst++)=(float)*(src++);
  src++;
}

void copy44f33f ( float *src, float *dst )
{
  *(dst++)=*(src++);
  *(dst++)=*(src++);
  *(dst++)=*(src++);
  src++;

  *(dst++)=*(src++);
  *(dst++)=*(src++);
  *(dst++)=*(src++);
  src++;

  *(dst++)=*(src++);
  *(dst++)=*(src++);
  *(dst++)=*(src++);
  src++;
}

void copy44f44d ( float *src, double *dst)
{
  *(dst++)=(double)*(src++);
  *(dst++)=(double)*(src++);
  *(dst++)=(double)*(src++);
  *(dst++)=(double)*(src++);

  *(dst++)=(double)*(src++);
  *(dst++)=(double)*(src++);
  *(dst++)=(double)*(src++);
  *(dst++)=(double)*(src++);

  *(dst++)=(double)*(src++);
  *(dst++)=(double)*(src++);
  *(dst++)=(double)*(src++);
  *(dst++)=(double)*(src++);

  *(dst++)=(double)*(src++);
  *(dst++)=(double)*(src++);
  *(dst++)=(double)*(src++);
  *(dst++)=(double)*(src++);
}

void copy44d44f ( double *src, float *dst)
{
  *(dst++)=(float)*(src++);
  *(dst++)=(float)*(src++);
  *(dst++)=(float)*(src++);
  *(dst++)=(float)*(src++);

  *(dst++)=(float)*(src++);
  *(dst++)=(float)*(src++);
  *(dst++)=(float)*(src++);
  *(dst++)=(float)*(src++);

  *(dst++)=(float)*(src++);
  *(dst++)=(float)*(src++);
  *(dst++)=(float)*(src++);
  *(dst++)=(float)*(src++);

  *(dst++)=(float)*(src++);
  *(dst++)=(float)*(src++);
  *(dst++)=(float)*(src++);
  *(dst++)=(float)*(src++);
}


void copy33f44d ( float *src, double *dst )
{
  const float _0 = 0.0;
  *(dst++)=(double)*(src++);
  *(dst++)=(double)*(src++);
  *(dst++)=(double)*(src++);
  *(dst++)=_0;

  *(dst++)=(double)*(src++);
  *(dst++)=(double)*(src++);
  *(dst++)=(double)*(src++);
  *(dst++)=_0;

  *(dst++)=(double)*(src++);
  *(dst++)=(double)*(src++);
  *(dst++)=(double)*(src++);
  *(dst++)=_0;

  *(dst++)=_0;
  *(dst++)=_0;
  *(dst++)=_0;
  *(dst++)=_1;
}

void copy33f44f ( float *src, float *dst )
{
  const float _0 = 0.0;
  *(dst++)=*(src++);
  *(dst++)=*(src++);
  *(dst++)=*(src++);
  *(dst++)=_0;

  *(dst++)=*(src++);
  *(dst++)=*(src++);
  *(dst++)=*(src++);
  *(dst++)=_0;

  *(dst++)=*(src++);
  *(dst++)=*(src++);
  *(dst++)=*(src++);
  *(dst++)=_0;

  *(dst++)=_0;
  *(dst++)=_0;
  *(dst++)=_0;
  *(dst++)=_1;
}

void transform33f3f (float *m1, float *m2, float *m3) 
{
  float m2r0=m2[0];
  float m2r1=m2[1];
  float m2r2=m2[2];
  m3[0] = m1[ 0] * m2r0 + m1[ 1] * m2r1 + m1[ 2] * m2r2;
  m3[1] = m1[ 3] * m2r0 + m1[ 4] * m2r1 + m1[ 5] * m2r2;
  m3[2] = m1[ 6] * m2r0 + m1[ 7] * m2r1 + m1[ 8] * m2r2;
}

void transpose33f33f ( float  *m1, float  *m2)
{
  m2[0] = m1[0];
  m2[1] = m1[3];
  m2[2] = m1[6];
  m2[3] = m1[1];
  m2[4] = m1[4];
  m2[5] = m1[7];
  m2[6] = m1[2];
  m2[7] = m1[5];
  m2[8] = m1[8];
}

void transpose44f44f ( float  *m1, float  *m2)
{
  m2[0 ] = m1[0 ];
  m2[1 ] = m1[4 ];
  m2[2 ] = m1[8 ];
  m2[3 ] = m1[12];

  m2[4 ] = m1[1 ];
  m2[5 ] = m1[5 ];
  m2[6 ] = m1[9 ];
  m2[7 ] = m1[13];

  m2[8 ] = m1[2 ];
  m2[9 ] = m1[6 ];
  m2[10] = m1[10];
  m2[11] = m1[14];

  m2[12] = m1[3 ];
  m2[13] = m1[7 ];
  m2[14] = m1[11];
  m2[15] = m1[15];
}
void transpose44d44d ( double *m1, double  *m2)
{
  m2[0 ] = m1[0 ];
  m2[1 ] = m1[4 ];
  m2[2 ] = m1[8 ];
  m2[3 ] = m1[12];

  m2[4 ] = m1[1 ];
  m2[5 ] = m1[5 ];
  m2[6 ] = m1[9 ];
  m2[7 ] = m1[13];

  m2[8 ] = m1[2 ];
  m2[9 ] = m1[6 ];
  m2[10] = m1[10];
  m2[11] = m1[14];

  m2[12] = m1[3 ];
  m2[13] = m1[7 ];
  m2[14] = m1[11];
  m2[15] = m1[15];
}

void transform33Tf3f (float *m1, float *m2, float *m3) 
{
  float m2r0=m2[0];
  float m2r1=m2[1];
  float m2r2=m2[2];
  m3[0] = m1[ 0] * m2r0 + m1[ 3] * m2r1 + m1[ 6] * m2r2;
  m3[1] = m1[ 1] * m2r0 + m1[ 4] * m2r1 + m1[ 7] * m2r2;
  m3[2] = m1[ 2] * m2r0 + m1[ 5] * m2r1 + m1[ 8] * m2r2;
}

void transform44f3f (float *m1, float *m2, float *m3)
{
  register float m2r0 = m2[0];
  register float m2r1 = m2[1];
  register float m2r2 = m2[2];
  m3[0] = m1[ 0] * m2r0 + m1[ 1] * m2r1 + m1[ 2] * m2r2 + m1[ 3];
  m3[1] = m1[ 4] * m2r0 + m1[ 5] * m2r1 + m1[ 6] * m2r2 + m1[ 7];
  m3[2] = m1[ 8] * m2r0 + m1[ 9] * m2r1 + m1[10] * m2r2 + m1[11];
}

void transform44d3f (double *m1, float *m2, float *m3)
{
  register double m2r0 = m2[0];
  register double m2r1 = m2[1];
  register double m2r2 = m2[2];
  m3[0] = (float) (m1[ 0] * m2r0 + m1[ 1] * m2r1 + m1[ 2] * m2r2 + m1[ 3]);
  m3[1] = (float) (m1[ 4] * m2r0 + m1[ 5] * m2r1 + m1[ 6] * m2r2 + m1[ 7]);
  m3[2] = (float) (m1[ 8] * m2r0 + m1[ 9] * m2r1 + m1[10] * m2r2 + m1[11]);
}
void transform44d3d (double *m1, double *m2, double *m3)
{
  register double m2r0 = m2[0];
  register double m2r1 = m2[1];
  register double m2r2 = m2[2];
  m3[0] = (float) (m1[ 0] * m2r0 + m1[ 1] * m2r1 + m1[ 2] * m2r2 + m1[ 3]);
  m3[1] = (float) (m1[ 4] * m2r0 + m1[ 5] * m2r1 + m1[ 6] * m2r2 + m1[ 7]);
  m3[2] = (float) (m1[ 8] * m2r0 + m1[ 9] * m2r1 + m1[10] * m2r2 + m1[11]);
}

void transform44d3fas33d3f (double *m1, float *m2, float *m3)
{
  register double m2r0 = m2[0];
  register double m2r1 = m2[1];
  register double m2r2 = m2[2];
  m3[0] = (float) (m1[ 0] * m2r0 + m1[ 1] * m2r1 + m1[ 2] * m2r2);
  m3[1] = (float) (m1[ 4] * m2r0 + m1[ 5] * m2r1 + m1[ 6] * m2r2);
  m3[2] = (float) (m1[ 8] * m2r0 + m1[ 9] * m2r1 + m1[10] * m2r2);
}

void transform44f3fas33f3f (float *m1, float *m2, float *m3)
{
  register float m2r0 = m2[0];
  register float m2r1 = m2[1];
  register float m2r2 = m2[2];
  m3[0] =  (m1[ 0] * m2r0 + m1[ 1] * m2r1 + m1[ 2] * m2r2);
  m3[1] =  (m1[ 4] * m2r0 + m1[ 5] * m2r1 + m1[ 6] * m2r2);
  m3[2] =  (m1[ 8] * m2r0 + m1[ 9] * m2r1 + m1[10] * m2r2);
}

void inverse_transformC44f3f (float *m1, float *m2, float *m3)
{
  register float m2r0 = m2[0] - m1[12];
  register float m2r1 = m2[1] - m1[13];
  register float m2r2 = m2[2] - m1[14];
  m3[0] = (float) (m1[ 0] * m2r0 + m1[ 1] * m2r1 + m1[ 2] * m2r2);
  m3[1] = (float) (m1[ 4] * m2r0 + m1[ 5] * m2r1 + m1[ 6] * m2r2);
  m3[2] = (float) (m1[ 8] * m2r0 + m1[ 9] * m2r1 + m1[10] * m2r2);
}

void inverse_transform44f3f (float *m1, float *m2, float *m3)
{
  register float m2r0 = m2[0] - m1[3];
  register float m2r1 = m2[1] - m1[7];
  register float m2r2 = m2[2] - m1[11];
  m3[0] = (float) (m1[ 0] * m2r0 + m1[ 4] * m2r1 + m1[ 8] * m2r2);
  m3[1] = (float) (m1[ 1] * m2r0 + m1[ 5] * m2r1 + m1[ 9] * m2r2);
  m3[2] = (float) (m1[ 2] * m2r0 + m1[ 6] * m2r1 + m1[10] * m2r2);
}

void inverse_transform44d3f (double *m1, float *m2, float *m3)
{
  register double m2r0 = m2[0] - m1[3];
  register double m2r1 = m2[1] - m1[7];
  register double m2r2 = m2[2] - m1[11];
  m3[0] = (float) (m1[ 0] * m2r0 + m1[ 4] * m2r1 + m1[ 8] * m2r2);
  m3[1] = (float) (m1[ 1] * m2r0 + m1[ 5] * m2r1 + m1[ 9] * m2r2);
  m3[2] = (float) (m1[ 2] * m2r0 + m1[ 6] * m2r1 + m1[10] * m2r2);
}
void inverse_transform44d3d (double *m1, double *m2, double *m3)
{
  register double m2r0 = m2[0] - m1[3];
  register double m2r1 = m2[1] - m1[7];
  register double m2r2 = m2[2] - m1[11];
  m3[0] = (float) (m1[ 0] * m2r0 + m1[ 4] * m2r1 + m1[ 8] * m2r2);
  m3[1] = (float) (m1[ 1] * m2r0 + m1[ 5] * m2r1 + m1[ 9] * m2r2);
  m3[2] = (float) (m1[ 2] * m2r0 + m1[ 6] * m2r1 + m1[10] * m2r2);
}

void transform44f4f (float *m1, float *m2, float *m3)
{
  float m2r0 = m2[0];
  float m2r1 = m2[1];
  float m2r2 = m2[2];
  float m2r3 = m2[3];
  m3[0] = m1[ 0] * m2r0 + m1[ 1] * m2r1 + m1[ 2] * m2r2 + m1[ 3] * m2r3;
  m3[1] = m1[ 4] * m2r0 + m1[ 5] * m2r1 + m1[ 6] * m2r2 + m1[ 7] * m2r3;
  m3[2] = m1[ 8] * m2r0 + m1[ 9] * m2r1 + m1[10] * m2r2 + m1[11] * m2r3;
  m3[3] = m1[12] * m2r0 + m1[13] * m2r1 + m1[14] * m2r2 + m1[15] * m2r3; 
}


void initializeTTT44f ( float *m )
{
  int a;
  for(a=0;a<16;a++)
    m[a]=_0;
  for(a=0;a<4;a++)
    m[4*a+a]=_1;
}


void combineTTT44f44f( float *m1, float *m2, float *m3) 
/* WARNING: this routine is ill-conceived and essentially broken */

/* NOTE: this is NOT equivalent to 4x4 matrix multiplication.
   TTTs are designed for easily creating movies of rotating 
   bodies! */
{
  float m1_homo[16];
  float m2_homo[16];
  float *src,*dst;
  float pre[3],post[3];

  /* convert the existing TTT into a homogenous transformation matrix */
  
  convertTTTfR44f(m1, m1_homo);
  convertTTTfR44f(m2, m2_homo);
  
  /* combine the matrices */

  left_multiply44f44f(m1_homo, m2_homo);

  /* now use the origin from the most recent TTT */

  src = m1+12;
  invert3f3f(src, pre);

  transform44f3fas33f3f(m2_homo, pre, post);

  m2_homo[ 3] += post[0];
  m2_homo[ 7] += post[1];
  m2_homo[11] += post[2];

  dst = m2_homo+12;

  copy3f(src,dst);
  copy44f(m2_homo,m3);

}

void convertTTTfR44d( float *ttt, double *homo)
{
  /* takes the PyMOL-specific TTT matrix and 
     makes a homogenous 4x4 txf matrix homo of it */

  register double ttt_3  = (double)ttt[3];
  register double ttt_7  = (double)ttt[7];
  register double ttt_11 = (double)ttt[11];
  register double ttt_12 = (double)ttt[12];
  register double ttt_13 = (double)ttt[13];
  register double ttt_14 = (double)ttt[14];

  /*  dump44f(ttt,"ttt");*/

  homo[ 0] = (double)ttt[ 0];
  homo[ 1] = (double)ttt[ 1];
  homo[ 2] = (double)ttt[ 2];
  homo[ 4] = (double)ttt[ 4];
  homo[ 5] = (double)ttt[ 5];
  homo[ 6] = (double)ttt[ 6];
  homo[ 8] = (double)ttt[ 8];
  homo[ 9] = (double)ttt[ 9];
  homo[10] = (double)ttt[10];

  homo[ 3] = (homo[ 0] * ttt_12) + (homo[ 1] * ttt_13) + (homo[ 2] * ttt_14) + ttt_3;
  homo[ 7] = (homo[ 4] * ttt_12) + (homo[ 5] * ttt_13) + (homo[ 6] * ttt_14) + ttt_7;
  homo[11] = (homo[ 8] * ttt_12) + (homo[ 9] * ttt_13) + (homo[10] * ttt_14) + ttt_11;

  homo[12] = 0.0;
  homo[13] = 0.0;
  homo[14] = 0.0;
  homo[15] = 1.0;

  /*  dump44d(homo, "homo");*/

}

void convertTTTfR44f( float *ttt, float *homo)
{
  /* takes the PyMOL-specific TTT matrix and 
     makes a homogenous 4x4 txf matrix homo of it */

  register float ttt_3  = ttt[3];
  register float ttt_7  = ttt[7];
  register float ttt_11 = ttt[11];
  register float ttt_12 = ttt[12];
  register float ttt_13 = ttt[13];
  register float ttt_14 = ttt[14];

  homo[ 0] = ttt[ 0];
  homo[ 1] = ttt[ 1];
  homo[ 2] = ttt[ 2];
  homo[ 4] = ttt[ 4];
  homo[ 5] = ttt[ 5];
  homo[ 6] = ttt[ 6];
  homo[ 8] = ttt[ 8];
  homo[ 9] = ttt[ 9];
  homo[10] = ttt[10];

  homo[ 3] = (homo[ 0] * ttt_12) + (homo[ 1] * ttt_13) + (homo[ 2] * ttt_14) + ttt_3;
  homo[ 7] = (homo[ 4] * ttt_12) + (homo[ 5] * ttt_13) + (homo[ 6] * ttt_14) + ttt_7;
  homo[11] = (homo[ 8] * ttt_12) + (homo[ 9] * ttt_13) + (homo[10] * ttt_14) + ttt_11;

  homo[12] = 0.0;
  homo[13] = 0.0;
  homo[14] = 0.0;
  homo[15] = 1.0;

}

void convert44d44f(double *dbl, float *flt)
{
  flt[ 0] = (float)dbl[ 0];
  flt[ 1] = (float)dbl[ 1];
  flt[ 2] = (float)dbl[ 2];
  flt[ 3] = (float)dbl[ 3];
  flt[ 4] = (float)dbl[ 4];
  flt[ 5] = (float)dbl[ 5];
  flt[ 6] = (float)dbl[ 6];
  flt[ 7] = (float)dbl[ 7];
  flt[ 8] = (float)dbl[ 8];
  flt[ 9] = (float)dbl[ 9];
  flt[10] = (float)dbl[10];
  flt[11] = (float)dbl[11];
  flt[12] = (float)dbl[12];
  flt[13] = (float)dbl[13];
  flt[14] = (float)dbl[14];
  flt[15] = (float)dbl[15];
}

void convert44f44d(float *flt, double *dbl)
{
  dbl[ 0] = (double)flt[ 0];
  dbl[ 1] = (double)flt[ 1];
  dbl[ 2] = (double)flt[ 2];
  dbl[ 3] = (double)flt[ 3];
  dbl[ 4] = (double)flt[ 4];
  dbl[ 5] = (double)flt[ 5];
  dbl[ 6] = (double)flt[ 6];
  dbl[ 7] = (double)flt[ 7];
  dbl[ 8] = (double)flt[ 8];
  dbl[ 9] = (double)flt[ 9];
  dbl[10] = (double)flt[10];
  dbl[11] = (double)flt[11];
  dbl[12] = (double)flt[12];
  dbl[13] = (double)flt[13];
  dbl[14] = (double)flt[14];
  dbl[15] = (double)flt[15];
}

void convertR44dTTTf( double *homo, float *ttt )
{
  /* nowadays, homogeneous matrices with (0,0,0,1) in 4th row are TTT
     compatible */
  convert44d44f(homo,ttt);
}

void multiply44d44d44d( double *left, double *right, double *product)
{
  register double rA = right[ 0];
  register double rB = right[ 4];
  register double rC = right[ 8];
  register double rD = right[12];
  
  product[ 0] = left[ 0] * rA + left[ 1] * rB + left[ 2] * rC + left[ 3] * rD;
  product[ 4] = left[ 4] * rA + left[ 5] * rB + left[ 6] * rC + left[ 7] * rD;
  product[ 8] = left[ 8] * rA + left[ 9] * rB + left[10] * rC + left[11] * rD;
  product[12] = left[12] * rA + left[13] * rB + left[14] * rC + left[15] * rD;

  rA = right[ 1];
  rB = right[ 5];
  rC = right[ 9];
  rD = right[13];
  
  product[ 1] = left[ 0] * rA + left[ 1] * rB + left[ 2] * rC + left[ 3] * rD;
  product[ 5] = left[ 4] * rA + left[ 5] * rB + left[ 6] * rC + left[ 7] * rD;
  product[ 9] = left[ 8] * rA + left[ 9] * rB + left[10] * rC + left[11] * rD;
  product[13] = left[12] * rA + left[13] * rB + left[14] * rC + left[15] * rD;

  rA = right[ 2];
  rB = right[ 6];
  rC = right[10];
  rD = right[14];
  
  product[ 2] = left[ 0] * rA + left[ 1] * rB + left[ 2] * rC + left[ 3] * rD;
  product[ 6] = left[ 4] * rA + left[ 5] * rB + left[ 6] * rC + left[ 7] * rD;
  product[10] = left[ 8] * rA + left[ 9] * rB + left[10] * rC + left[11] * rD;
  product[14] = left[12] * rA + left[13] * rB + left[14] * rC + left[15] * rD;

  rA = right[ 3];
  rB = right[ 7];
  rC = right[11];
  rD = right[15];
  
  product[ 3] = left[ 0] * rA + left[ 1] * rB + left[ 2] * rC + left[ 3] * rD;
  product[ 7] = left[ 4] * rA + left[ 5] * rB + left[ 6] * rC + left[ 7] * rD;
  product[11] = left[ 8] * rA + left[ 9] * rB + left[10] * rC + left[11] * rD;
  product[15] = left[12] * rA + left[13] * rB + left[14] * rC + left[15] * rD;
}

void left_multiply44d44d( double *left, double *right)
{
  register double rA = right[ 0];
  register double rB = right[ 4];
  register double rC = right[ 8];
  register double rD = right[12];
  
  right[ 0] = left[ 0] * rA + left[ 1] * rB + left[ 2] * rC + left[ 3] * rD;
  right[ 4] = left[ 4] * rA + left[ 5] * rB + left[ 6] * rC + left[ 7] * rD;
  right[ 8] = left[ 8] * rA + left[ 9] * rB + left[10] * rC + left[11] * rD;
  right[12] = left[12] * rA + left[13] * rB + left[14] * rC + left[15] * rD;

  rA = right[ 1];
  rB = right[ 5];
  rC = right[ 9];
  rD = right[13];
  
  right[ 1] = left[ 0] * rA + left[ 1] * rB + left[ 2] * rC + left[ 3] * rD;
  right[ 5] = left[ 4] * rA + left[ 5] * rB + left[ 6] * rC + left[ 7] * rD;
  right[ 9] = left[ 8] * rA + left[ 9] * rB + left[10] * rC + left[11] * rD;
  right[13] = left[12] * rA + left[13] * rB + left[14] * rC + left[15] * rD;

  rA = right[ 2];
  rB = right[ 6];
  rC = right[10];
  rD = right[14];
  
  right[ 2] = left[ 0] * rA + left[ 1] * rB + left[ 2] * rC + left[ 3] * rD;
  right[ 6] = left[ 4] * rA + left[ 5] * rB + left[ 6] * rC + left[ 7] * rD;
  right[10] = left[ 8] * rA + left[ 9] * rB + left[10] * rC + left[11] * rD;
  right[14] = left[12] * rA + left[13] * rB + left[14] * rC + left[15] * rD;

  rA = right[ 3];
  rB = right[ 7];
  rC = right[11];
  rD = right[15];
  
  right[ 3] = left[ 0] * rA + left[ 1] * rB + left[ 2] * rC + left[ 3] * rD;
  right[ 7] = left[ 4] * rA + left[ 5] * rB + left[ 6] * rC + left[ 7] * rD;
  right[11] = left[ 8] * rA + left[ 9] * rB + left[10] * rC + left[11] * rD;
  right[15] = left[12] * rA + left[13] * rB + left[14] * rC + left[15] * rD;
}

void right_multiply44d44d( double *left, double *right)
{
  register double cA = left[ 0];
  register double cB = left[ 1];
  register double cC = left[ 2];
  register double cD = left[ 3];
  
  left[ 0] = cA * right[ 0] + cB * right[ 4] + cC * right[ 8] + cD * right[12];
  left[ 1] = cA * right[ 1] + cB * right[ 5] + cC * right[ 9] + cD * right[13];
  left[ 2] = cA * right[ 2] + cB * right[ 6] + cC * right[10] + cD * right[14];
  left[ 3] = cA * right[ 3] + cB * right[ 7] + cC * right[11] + cD * right[15];

  cA = left[ 4];
  cB = left[ 5];
  cC = left[ 6];
  cD = left[ 7];
  
  left[ 4] = cA * right[ 0] + cB * right[ 4] + cC * right[ 8] + cD * right[12];
  left[ 5] = cA * right[ 1] + cB * right[ 5] + cC * right[ 9] + cD * right[13];
  left[ 6] = cA * right[ 2] + cB * right[ 6] + cC * right[10] + cD * right[14];
  left[ 7] = cA * right[ 3] + cB * right[ 7] + cC * right[11] + cD * right[15];

  cA = left[ 8];
  cB = left[ 9];
  cC = left[10];
  cD = left[11];
  
  left[ 8] = cA * right[ 0] + cB * right[ 4] + cC * right[ 8] + cD * right[12];
  left[ 9] = cA * right[ 1] + cB * right[ 5] + cC * right[ 9] + cD * right[13];
  left[10] = cA * right[ 2] + cB * right[ 6] + cC * right[10] + cD * right[14];
  left[11] = cA * right[ 3] + cB * right[ 7] + cC * right[11] + cD * right[15];

  cA = left[12];
  cB = left[13];
  cC = left[14];
  cD = left[15];
  
  left[12] = cA * right[ 0] + cB * right[ 4] + cC * right[ 8] + cD * right[12];
  left[13] = cA * right[ 1] + cB * right[ 5] + cC * right[ 9] + cD * right[13];
  left[14] = cA * right[ 2] + cB * right[ 6] + cC * right[10] + cD * right[14];
  left[15] = cA * right[ 3] + cB * right[ 7] + cC * right[11] + cD * right[15];

}



void multiply44f44f44f( float *left, float *right, float *product)
{
  register float rA = right[ 0];
  register float rB = right[ 4];
  register float rC = right[ 8];
  register float rD = right[12];
  
  product[ 0] = left[ 0] * rA + left[ 1] * rB + left[ 2] * rC + left[ 3] * rD;
  product[ 4] = left[ 4] * rA + left[ 5] * rB + left[ 6] * rC + left[ 7] * rD;
  product[ 8] = left[ 8] * rA + left[ 9] * rB + left[10] * rC + left[11] * rD;
  product[12] = left[12] * rA + left[13] * rB + left[14] * rC + left[15] * rD;

  rA = right[ 1];
  rB = right[ 5];
  rC = right[ 9];
  rD = right[13];
  
  product[ 1] = left[ 0] * rA + left[ 1] * rB + left[ 2] * rC + left[ 3] * rD;
  product[ 5] = left[ 4] * rA + left[ 5] * rB + left[ 6] * rC + left[ 7] * rD;
  product[ 9] = left[ 8] * rA + left[ 9] * rB + left[10] * rC + left[11] * rD;
  product[13] = left[12] * rA + left[13] * rB + left[14] * rC + left[15] * rD;

  rA = right[ 2];
  rB = right[ 6];
  rC = right[10];
  rD = right[14];
  
  product[ 2] = left[ 0] * rA + left[ 1] * rB + left[ 2] * rC + left[ 3] * rD;
  product[ 6] = left[ 4] * rA + left[ 5] * rB + left[ 6] * rC + left[ 7] * rD;
  product[10] = left[ 8] * rA + left[ 9] * rB + left[10] * rC + left[11] * rD;
  product[14] = left[12] * rA + left[13] * rB + left[14] * rC + left[15] * rD;

  rA = right[ 3];
  rB = right[ 7];
  rC = right[11];
  rD = right[15];
  
  product[ 3] = left[ 0] * rA + left[ 1] * rB + left[ 2] * rC + left[ 3] * rD;
  product[ 7] = left[ 4] * rA + left[ 5] * rB + left[ 6] * rC + left[ 7] * rD;
  product[11] = left[ 8] * rA + left[ 9] * rB + left[10] * rC + left[11] * rD;
  product[15] = left[12] * rA + left[13] * rB + left[14] * rC + left[15] * rD;
}

void left_multiply44f44f( float *left, float *right)
{
  register float rA = right[ 0];
  register float rB = right[ 4];
  register float rC = right[ 8];
  register float rD = right[12];
  
  right[ 0] = left[ 0] * rA + left[ 1] * rB + left[ 2] * rC + left[ 3] * rD;
  right[ 4] = left[ 4] * rA + left[ 5] * rB + left[ 6] * rC + left[ 7] * rD;
  right[ 8] = left[ 8] * rA + left[ 9] * rB + left[10] * rC + left[11] * rD;
  right[12] = left[12] * rA + left[13] * rB + left[14] * rC + left[15] * rD;

  rA = right[ 1];
  rB = right[ 5];
  rC = right[ 9];
  rD = right[13];
  
  right[ 1] = left[ 0] * rA + left[ 1] * rB + left[ 2] * rC + left[ 3] * rD;
  right[ 5] = left[ 4] * rA + left[ 5] * rB + left[ 6] * rC + left[ 7] * rD;
  right[ 9] = left[ 8] * rA + left[ 9] * rB + left[10] * rC + left[11] * rD;
  right[13] = left[12] * rA + left[13] * rB + left[14] * rC + left[15] * rD;

  rA = right[ 2];
  rB = right[ 6];
  rC = right[10];
  rD = right[14];
  
  right[ 2] = left[ 0] * rA + left[ 1] * rB + left[ 2] * rC + left[ 3] * rD;
  right[ 6] = left[ 4] * rA + left[ 5] * rB + left[ 6] * rC + left[ 7] * rD;
  right[10] = left[ 8] * rA + left[ 9] * rB + left[10] * rC + left[11] * rD;
  right[14] = left[12] * rA + left[13] * rB + left[14] * rC + left[15] * rD;

  rA = right[ 3];
  rB = right[ 7];
  rC = right[11];
  rD = right[15];
  
  right[ 3] = left[ 0] * rA + left[ 1] * rB + left[ 2] * rC + left[ 3] * rD;
  right[ 7] = left[ 4] * rA + left[ 5] * rB + left[ 6] * rC + left[ 7] * rD;
  right[11] = left[ 8] * rA + left[ 9] * rB + left[10] * rC + left[11] * rD;
  right[15] = left[12] * rA + left[13] * rB + left[14] * rC + left[15] * rD;
}

void right_multiply44f44f( float *left, float *right)
{
  register float cA = left[ 0];
  register float cB = left[ 1];
  register float cC = left[ 2];
  register float cD = left[ 3];
  
  left[ 0] = cA * right[ 0] + cB * right[ 4] + cC * right[ 8] + cD * right[12];
  left[ 1] = cA * right[ 1] + cB * right[ 5] + cC * right[ 9] + cD * right[13];
  left[ 2] = cA * right[ 2] + cB * right[ 6] + cC * right[10] + cD * right[14];
  left[ 3] = cA * right[ 3] + cB * right[ 7] + cC * right[11] + cD * right[15];

  cA = left[ 4];
  cB = left[ 5];
  cC = left[ 6];
  cD = left[ 7];
  
  left[ 4] = cA * right[ 0] + cB * right[ 4] + cC * right[ 8] + cD * right[12];
  left[ 5] = cA * right[ 1] + cB * right[ 5] + cC * right[ 9] + cD * right[13];
  left[ 6] = cA * right[ 2] + cB * right[ 6] + cC * right[10] + cD * right[14];
  left[ 7] = cA * right[ 3] + cB * right[ 7] + cC * right[11] + cD * right[15];

  cA = left[ 8];
  cB = left[ 9];
  cC = left[10];
  cD = left[11];
  
  left[ 8] = cA * right[ 0] + cB * right[ 4] + cC * right[ 8] + cD * right[12];
  left[ 9] = cA * right[ 1] + cB * right[ 5] + cC * right[ 9] + cD * right[13];
  left[10] = cA * right[ 2] + cB * right[ 6] + cC * right[10] + cD * right[14];
  left[11] = cA * right[ 3] + cB * right[ 7] + cC * right[11] + cD * right[15];

  cA = left[12];
  cB = left[13];
  cC = left[14];
  cD = left[15];
  
  left[12] = cA * right[ 0] + cB * right[ 4] + cC * right[ 8] + cD * right[12];
  left[13] = cA * right[ 1] + cB * right[ 5] + cC * right[ 9] + cD * right[13];
  left[14] = cA * right[ 2] + cB * right[ 6] + cC * right[10] + cD * right[14];
  left[15] = cA * right[ 3] + cB * right[ 7] + cC * right[11] + cD * right[15];

}


void invert_special44d44d(double *orig, double *inv)
{
  /* inverse of the rotation matrix */

  inv[ 0] = orig[ 0];
  inv[ 1] = orig[ 4];
  inv[ 2] = orig[ 8];
  inv[ 4] = orig[ 1];
  inv[ 5] = orig[ 5];
  inv[ 6] = orig[ 9];
  inv[ 8] = orig[ 2];
  inv[ 9] = orig[ 6];
  inv[10] = orig[10];

  /* invert the translation portion */

  inv[ 3] = -(orig[ 3] * orig[ 0] + orig[ 7] * orig[ 4] + orig[11] * orig[ 8]);
  inv[ 7] = -(orig[ 3] * orig[ 1] + orig[ 7] * orig[ 5] + orig[11] * orig[ 9]);
  inv[11] = -(orig[ 3] * orig[ 2] + orig[ 7] * orig[ 6] + orig[11] * orig[10]);

  inv[12] = 0.0;
  inv[13] = 0.0;
  inv[14] = 0.0;
  inv[15] = 1.0;

}

void invert_special44f44f(float *orig, float *inv)
{
  /* inverse of the rotation matrix */

  inv[ 0] = orig[ 0];
  inv[ 1] = orig[ 4];
  inv[ 2] = orig[ 8];
  inv[ 4] = orig[ 1];
  inv[ 5] = orig[ 5];
  inv[ 6] = orig[ 9];
  inv[ 8] = orig[ 2];
  inv[ 9] = orig[ 6];
  inv[10] = orig[10];

  /* invert the translation portion */

  inv[ 3] = -(orig[ 3] * orig[ 0] + orig[ 7] * orig[ 4] + orig[11] * orig[ 8]);
  inv[ 7] = -(orig[ 3] * orig[ 1] + orig[ 7] * orig[ 5] + orig[11] * orig[ 9]);
  inv[11] = -(orig[ 3] * orig[ 2] + orig[ 7] * orig[ 6] + orig[11] * orig[10]);

  inv[12] = 0.0F;
  inv[13] = 0.0F;
  inv[14] = 0.0F;
  inv[15] = 1.0F;

}

static void normalize3dp( double *v1, double *v2, double *v3 )
{
  double vlen = sqrt1d((v1[0]*v1[0]) + 
                       (v2[0]*v2[0]) + 
                       (v3[0]*v3[0]));
  if(vlen>R_SMALL)
	 {
		v1[0]/=vlen;
		v2[0]/=vlen;
		v3[0]/=vlen;
	 }
  else
	 {
		v1[0]=_0;
		v2[1]=_0;
		v3[2]=_0;
	 }
} 
/* unused at present
static void normalize3df( float *v1, float *v2, float *v3 )
{
  float vlen = (float)sqrt1f((v1[0]*v1[0]) + 
                             (v2[0]*v2[0]) + 
                             (v3[0]*v3[0]));
  if(vlen>R_SMALL)
    {
      v1[0]/=vlen;
      v2[0]/=vlen;
      v3[0]/=vlen;
    }
  else
    {
      v1[0]=_0;
      v2[1]=_0;
      v3[2]=_0;
    }
} 
*/
void scale3d ( double *v1,double v0,double *v2)
{
  v2[0]=v1[0]*v0;
  v2[1]=v1[1]*v0;
  v2[2]=v1[2]*v0;
}

void add3d ( double *v1, double *v2, double *v3 )
{
  v3[0]=v1[0]+v2[0];
  v3[1]=v1[1]+v2[1];
  v3[2]=v1[2]+v2[2];
}

void cross_product3d ( double *v1, double *v2, double *cross )
{
  cross[0] = (v1[1]*v2[2]) - (v1[2]*v2[1]);
  cross[1] = (v1[2]*v2[0]) - (v1[0]*v2[2]);
  cross[2] = (v1[0]*v2[1]) - (v1[1]*v2[0]);
}

void remove_component3d ( double *v1, double *unit, double *result)
{
  double dot;

  dot = v1[0]*unit[0] + v1[1]*unit[1] + v1[2]*unit[2];
  result[0]=v1[0]-unit[0]*dot;
  result[1]=v1[1]-unit[1]*dot;
  result[2]=v1[2]-unit[2]*dot;  
}

void reorient44d(double *matrix)
{
  double tmp[16];
  int a;

  /* restore orthogonality and recondition */

  for(a=0;a<3;a++) {
    normalize3d(matrix);
    normalize3d(matrix+4);
    normalize3d(matrix+8);
    cross_product3d(matrix+4, matrix+8, tmp);
    cross_product3d(matrix+8, matrix, tmp+4);
    cross_product3d(matrix, matrix+4, tmp+8);
    normalize3d(tmp);
    normalize3d(tmp+4);
    normalize3d(tmp+8);
    scale3d(tmp,2.0,tmp);
    scale3d(tmp+4,2.0,tmp+4);
    scale3d(tmp+8,2.0,tmp+8);
    add3d(matrix,tmp,tmp);
    add3d(matrix+4,tmp+4,tmp+4);
    add3d(matrix+8,tmp+8,tmp+8);
    copy3d(tmp,matrix);
    copy3d(tmp+4,matrix+4);
    copy3d(tmp+8,matrix+8);
  }
   
  normalize3d(matrix);
  normalize3d(matrix+4);
  normalize3d(matrix+8);
  
  copy3d(matrix,tmp);
  remove_component3d(matrix+4,tmp,tmp+4);
  cross_product3d(tmp, tmp+4, tmp+8);    
  normalize3d(tmp+4);
  normalize3d(tmp+8);
  
  recondition44d(tmp);

  copy3d(tmp,matrix);
  copy3d(tmp+4,matrix+4);
  copy3d(tmp+8,matrix+8);
  
}

void recondition33d(double *matrix)
{
  normalize3d(matrix);
  normalize3d(matrix+3);
  normalize3d(matrix+6);
  normalize3dp(matrix + 0, matrix + 3, matrix + 6);
  normalize3dp(matrix + 1, matrix + 4, matrix + 7);
  normalize3dp(matrix + 2, matrix + 5, matrix + 8);
  normalize3d(matrix);
  normalize3d(matrix+3);
  normalize3d(matrix+6);
  normalize3dp(matrix + 0, matrix + 3, matrix + 6 );
  normalize3dp(matrix + 1, matrix + 4, matrix + 7 );
  normalize3dp(matrix + 2, matrix + 5, matrix + 8);
  normalize3d(matrix);
  normalize3d(matrix+3);
  normalize3d(matrix+6);
}

void recondition44d(double *matrix)
{
  normalize3d(matrix);
  normalize3d(matrix+4);
  normalize3d(matrix+8);
  normalize3dp(matrix + 0, matrix + 4, matrix + 8 );
  normalize3dp(matrix + 1, matrix + 5, matrix + 9 );
  normalize3dp(matrix + 2, matrix + 6, matrix + 10);
  normalize3d(matrix);
  normalize3d(matrix+4);
  normalize3d(matrix+8);
  normalize3dp(matrix + 0, matrix + 4, matrix + 8 );
  normalize3dp(matrix + 1, matrix + 5, matrix + 9 );
  normalize3dp(matrix + 2, matrix + 6, matrix + 10);
  normalize3d(matrix);
  normalize3d(matrix+4);
  normalize3d(matrix+8);
}

void invert_rotation_only44d44d(double *orig, double *inv)
{
  /* inverse of the rotation matrix */

  inv[ 0] = orig[ 0];
  inv[ 1] = orig[ 4];
  inv[ 2] = orig[ 8];
  inv[ 4] = orig[ 1];
  inv[ 5] = orig[ 5];
  inv[ 6] = orig[ 9];
  inv[ 8] = orig[ 2];
  inv[ 9] = orig[ 6];
  inv[10] = orig[10];

  inv[ 3] = 0.0;
  inv[ 7] = 0.0;
  inv[11] = 0.0;

  inv[12] = 0.0;
  inv[13] = 0.0;
  inv[14] = 0.0;
  inv[15] = 1.0;

}


void transformTTT44f3f (float *m1, float *m2, float *m3)
{
  float m2r0 = m2[0] + m1[12];
  float m2r1 = m2[1] + m1[13];
  float m2r2 = m2[2] + m1[14];
  m3[0] = m1[ 0] * m2r0 + m1[ 1] * m2r1 + m1[ 2] * m2r2 + m1[3];
  m3[1] = m1[ 4] * m2r0 + m1[ 5] * m2r1 + m1[ 6] * m2r2 + m1[7];
  m3[2] = m1[ 8] * m2r0 + m1[ 9] * m2r1 + m1[10] * m2r2 + m1[11];
}

void transform_normalTTT44f3f (float *m1, float *m2, float *m3)
{
  float m2r0 = m2[0];
  float m2r1 = m2[1];
  float m2r2 = m2[2];
  m3[0] = m1[ 0] * m2r0 + m1[ 1] * m2r1 + m1[ 2] * m2r2 ;
  m3[1] = m1[ 4] * m2r0 + m1[ 5] * m2r1 + m1[ 6] * m2r2 ;
  m3[2] = m1[ 8] * m2r0 + m1[ 9] * m2r1 + m1[10] * m2r2 ;
}

void multiply33f33f ( float *m1,float *m2,float *m3) /* m2 and m3 can be the same matrix */
{
  int a;
  float m2r0,m2r1,m2r2;
  for(a=0;a<3;a++) {
    m2r0=m2[a];
    m2r1=m2[3+a];
    m2r2=m2[6+a];
    m3[a  ] = m1[0]*m2r0 + m1[1]*m2r1 + m1[2]*m2r2;
    m3[3+a] = m1[3]*m2r0 + m1[4]*m2r1 + m1[5]*m2r2;
    m3[6+a] = m1[6]*m2r0 + m1[7]*m2r1 + m1[8]*m2r2;
  }
}

void multiply33d33d ( double *m1,double *m2,double *m3) /* m2 and m3 can be the same matrix */
{
  int a;
  double m2r0,m2r1,m2r2;
  for(a=0;a<3;a++) {
    m2r0=m2[a];
    m2r1=m2[3+a];
    m2r2=m2[6+a];
    m3[a  ] = m1[0]*m2r0 + m1[1]*m2r1 + m1[2]*m2r2;
    m3[3+a] = m1[3]*m2r0 + m1[4]*m2r1 + m1[5]*m2r2;
    m3[6+a] = m1[6]*m2r0 + m1[7]*m2r1 + m1[8]*m2r2;
  }
}

void matrix_multiply33f33f ( Matrix33f m1,Matrix33f m2,Matrix33f m3)
{
  multiply33f33f((float*)m1,(float*)m2,(float*)m3);
}

void matrix_multiply33d33d ( Matrix33d m1,Matrix33d m2,Matrix33d m3)
{
  multiply33d33d((double*)m1[0],(double*)m2,(double*)m3);
}

float deg_to_rad(float angle)
{
  return((float)((angle* cPI)/180.0));
}
float rad_to_deg(float angle)
{
  return((float)(180.0 * (angle/ cPI)));
}

void get_rotation_about3f3fTTTf(float angle, float *dir, float *origin, float *ttt)
{
  float rot[9];
  rotation_matrix3f(angle,dir[0],dir[1],dir[2],rot);
  ttt[ 0] = rot[0];
  ttt[ 1] = rot[1];
  ttt[ 2] = rot[2];
  ttt[ 4] = rot[3];
  ttt[ 5] = rot[4];
  ttt[ 6] = rot[5];
  ttt[ 8] = rot[6];
  ttt[ 9] = rot[7];
  ttt[10] = rot[8];
  ttt[12] = -origin[0];
  ttt[13] = -origin[1];
  ttt[14] = -origin[2];
  ttt[ 3] = origin[0];
  ttt[ 7] = origin[1];
  ttt[11] = origin[2];
  ttt[15] = 1.0F;
}

void rotation_to_matrix33f(float *axis, float angle, Matrix33f mat)
{
  rotation_matrix3f(angle,axis[0],axis[1],axis[2],&mat[0][0]);
}

void rotation_matrix3f( float angle, float x, float y, float z,float *m )
{
  /* returns a row-major rotation matrix */

  int a,b;

   /* This function contributed by Erich Boleyn (erich@uruk.org) */
   float mag, s, c;
   float xx, yy, zz, xy, yz, zx, xs, ys, zs, one_c;

   s = (float)sin(angle);
   c = (float)cos(angle);

   mag = (float)sqrt1f(x*x + y*y + z*z);
	if(mag>=R_SMALL) {
	  x /= mag;
	  y /= mag;
	  z /= mag;
	  
#define M(row,col)  m[row*3+col]
	  
	  xx = x * x;
	  yy = y * y;
	  zz = z * z;
	  xy = x * y;
	  yz = y * z;
	  zx = z * x;
	  xs = x * s;
	  ys = y * s;
	  zs = z * s;
	  one_c = _1 - c;
	  
	  M(0,0) = (one_c * xx) + c;
	  M(0,1) = (one_c * xy) - zs;
	  M(0,2) = (one_c * zx) + ys;
	  
	  M(1,0) = (one_c * xy) + zs;
	  M(1,1) = (one_c * yy) + c;
	  M(1,2) = (one_c * yz) - xs;
	  
	  M(2,0) = (one_c * zx) - ys;
	  M(2,1) = (one_c * yz) + xs;
	  M(2,2) = (one_c * zz) + c;
	}	else {
	  for(a=0;a<3;a++)
		 for(b=0;b<3;b++)
			M(a,b)=0;
	  M(0,0)=_1;
	  M(1,1)=_1;
	  M(2,2)=_1;
	}

}


#define get_angle USED_TO_RETURN_DEGREES

float get_dihedral3f( float *v0, float *v1, float *v2, float *v3 )
{
  Vector3f d01,d21,d32,dd1,dd3,pos_d;
  float result=_0;

  subtract3f(v2,v1,d21);
  subtract3f(v0,v1,d01);
  subtract3f(v3,v2,d32);
  if (length3f(d21)<R_SMALL) {
    result = get_angle3f(d01,d32);
  } else { 
    cross_product3f(d21,d01,dd1);
    cross_product3f(d21,d32,dd3);
    if ((length3f(dd1)<R_SMALL)||(length3f(dd3)<R_SMALL)) {/* degenerate cases */
      result = get_angle3f(d01,d32); /* fall back to angle between vectors */
    }    else {
      result = get_angle3f(dd1,dd3);
      cross_product3f(d21,dd1,pos_d);
      if(dot_product3f(dd3,pos_d)<_0)
        result = -result;
    }
  }
  return(result);
}

float get_angle3f( float *v1, float *v2 )
{
  double denom;
  double result;

  denom = sqrt1f(((v1[0]*v1[0]) + 
                  (v1[1]*v1[1]) + 
                  (v1[2]*v1[2]))) *
    sqrt1f(((v2[0]*v2[0]) + 
            (v2[1]*v2[1]) + 
            (v2[2]*v2[2])));
  
  if(denom>R_SMALL) 
    result = ( v1[0]*v2[0] + 
					v1[1]*v2[1] + 
					v1[2]*v2[2] ) / denom;
  else
    result = _0;
  if(result<-_1) 
    result=-_1;
  else if(result>_1) 
    result=_1;
  result = acos(result);
  return((float)result);
} 

void normalize23f( float *v1 , float *v2)
{
  double vlen;
  vlen = length3f(v1);
  if(vlen>R_SMALL)
	 {
		v2[0]=(float)(v1[0]/vlen);
		v2[1]=(float)(v1[1]/vlen);
		v2[2]=(float)(v1[2]/vlen);
	 }
  else
	 {
		v2[0]=_0;
		v2[1]=_0;
		v2[2]=_0;
	 }
}

void clamp3f(float *v1)
{
  if(v1[0]<_0) v1[0]= _0;
  if(v1[0]>_1) v1[0]= _1;
  if(v1[1]<_0) v1[1]= _0;
  if(v1[1]>_1) v1[1]= _1;
  if(v1[2]<_0) v1[2]= _0;
  if(v1[2]>_1) v1[2]= _1;
}

void slow_normalize3f( float *v1 )
{
	double vlen = length3f(v1);
	if(vlen > R_SMALL)
	{
		float	inV	= (float)(_1 / vlen);
		v1[0] *= inV;
		v1[1] *= inV;
		v1[2] *= inV;
	}
	else
	{
		v1[0]=v1[1]=v1[2]=_0;
	}
} 

void normalize3d( double *v1 )
{
  double vlen;
  vlen = length3d(v1);
  if(vlen>R_SMALL)
	 {
		v1[0]/=vlen;
		v1[1]/=vlen;
		v1[2]/=vlen;
	 }
  else
	 {
		v1[0]=_0;
		v1[1]=_0;
		v1[2]=_0;
	 }
} 

double length3d ( double *v1 )
{
  return(sqrt1d((v1[0]*v1[0]) + 
					 (v1[1]*v1[1]) + 
					 (v1[2]*v1[2])));
} 

double slow_diff3f ( float *v1, float *v2 )
{
  register float dx,dy,dz;
  dx = (v1[0]-v2[0]);
  dy = (v1[1]-v2[1]);
  dz = (v1[2]-v2[2]);
  return(sqrt1d(dx*dx + dy*dy + dz*dz));
}

float slow_diffsq3f ( float *v1, float *v2 )
{
  register float dx,dy,dz;
  dx = (v1[0]-v2[0]);
  dy = (v1[1]-v2[1]);
  dz = (v1[2]-v2[2]);

  return ( dx*dx + dy*dy + dz*dz );

}

int slow_within3f(float *v1,float *v2,float dist)
{
  register float dx,dy,dz;
  dx = (v1[0]-v2[0]);
  if(fabs(dx)>dist) return(false);
  dy = (v1[1]-v2[1]);
  if(fabs(dy)>dist) return(false);
  dz = (v1[2]-v2[2]);
  if(fabs(dz)>dist) return(false);
  return((dx*dx + dy*dy + dz*dz)<=(dist*dist));
}

int slow_within3fsq(float *v1,float *v2,float dist,float dist2)
{
  register float dx,dy,dz;
  dx = (v1[0]-v2[0]);
  if(fabs(dx)>dist) return(false);
  dy = (v1[1]-v2[1]);
  if(fabs(dy)>dist) return(false);
  dz = (v1[2]-v2[2]);
  if(fabs(dz)>dist) return(false);
  return((dx*dx + dy*dy + dz*dz)<=(dist2));
}

int slow_within3fret(float *v1,float *v2,float cutoff, float cutoff2, float *diff, float *dist)
{
  register float dx,dy,dz,dist2;
  dx = (float)fabs( (diff[0] = v1[0]-v2[0]) );
  dy = (float)fabs( (diff[1] = v1[1]-v2[1]) );
  if(dx>cutoff) return 0;
  dz = (float)fabs( (diff[2] = v1[2]-v2[2]) );
  dx = dx * dx;
  if(dy>cutoff) return 0;
  dy = dy * dy;
  if(dz>cutoff) return 0;
  if((dist2 = ((dx + dy) + dz*dz))>cutoff2) 
    return 0;
  *dist = (float)sqrt1f(dist2);
  return 1;
}

void slow_remove_component3f ( float *v1, float *unit, float *result)
{
  float dot;

  dot = v1[0]*unit[0] + v1[1]*unit[1] + v1[2]*unit[2];
  result[0]=v1[0]-unit[0]*dot;
  result[1]=v1[1]-unit[1]*dot;
  result[2]=v1[2]-unit[2]*dot;  
}


double distance_line2point3f(float *base,float *normal,float *point,float *alongNormalSq)
{
  float hyp[3],adj[3];
  double result;

  hyp[0] = point[0] - base[0];
  hyp[1] = point[1] - base[1];
  hyp[2] = point[2] - base[2];

  project3f(hyp,normal,adj);

  (*alongNormalSq) = ((adj[0]*adj[0])+(adj[1]*adj[1])+(adj[2]*adj[2])); 
  result=((hyp[0]*hyp[0])+(hyp[1]*hyp[1])+(hyp[2]*hyp[2]))
	 - (*alongNormalSq);
  if(result<=_0) 
	 return(_0);
  else 
	 return(sqrt1d(result));

}

double distance_halfline2point3f(float *base,float *normal,float *point,float *alongNormalSq)
{
  float hyp[3],adj[3];
  double result;

  hyp[0] = point[0] - base[0];
  hyp[1] = point[1] - base[1];
  hyp[2] = point[2] - base[2];


  if(project3f(hyp,normal,adj)>_0)
	 {
		(*alongNormalSq) = ((adj[0]*adj[0])+(adj[1]*adj[1])+(adj[2]*adj[2])); 
		result=((hyp[0]*hyp[0])+(hyp[1]*hyp[1])+(hyp[2]*hyp[2]))
		  - (*alongNormalSq);
		if(result<=_0) 
		  return(_0);
		else 
		  return(sqrt1d(result));
	 } else {
		return(MAXFLOAT);
	 }
}


void matrix_transform33f3f(Matrix33f m1,float *v1,float *v2)
{
  v2[0] = m1[0][0]*v1[0] + m1[0][1]*v1[1] + m1[0][2]*v1[2];
  v2[1] = m1[1][0]*v1[0] + m1[1][1]*v1[1] + m1[1][2]*v1[2];
  v2[2] = m1[2][0]*v1[0] + m1[2][1]*v1[1] + m1[2][2]*v1[2];
}

void matrix_inverse_transform33f3f(Matrix33f m1,float *v1,float *v2)
{
	v2[0] = m1[0][0]*v1[0] + m1[1][0]*v1[1] + m1[2][0]*v1[2];
	v2[1] = m1[0][1]*v1[0] + m1[1][1]*v1[1] + m1[2][1]*v1[2];
	v2[2] = m1[0][2]*v1[0] + m1[1][2]*v1[1] + m1[2][2]*v1[2];
}


#if 0
double matdiffsq ( float *v1, oMatrix5f m, float *v2 )
{
  register double dx,dy,dz;
  float vx,vy,vz;
  
  dx = v2[0] - m[3][0];
  dy = v2[1] - m[3][1];
  dz = v2[2] - m[3][2];

  vx = m[0][0] * dx + m[0][1] * dy + m[0][2] * dz;
  vy = m[1][0] * dx + m[1][1] * dy + m[1][2] * dz;
  vz = m[2][0] * dx + m[2][1] * dy + m[2][2] * dz;

  dx = (v1[0]-(vx+m[4][0]));
  dy = (v1[1]-(vy+m[4][1]));
  dz = (v1[2]-(vz+m[4][2]));

  return( dx*dx + dy*dy + dz*dz );

}
#endif

void transform5f3f (oMatrix5f m, float *v1, float *v2 )
{
  register double dx,dy,dz;
  double vx,vy,vz;
  
  dx = v1[0] - m[3][0];
  dy = v1[1] - m[3][1];
  dz = v1[2] - m[3][2];

  vx = m[0][0] * dx + m[0][1] * dy + m[0][2] * dz;
  vy = m[1][0] * dx + m[1][1] * dy + m[1][2] * dz;
  vz = m[2][0] * dx + m[2][1] * dy + m[2][2] * dz;

  v2[0] = (((float)vx)+m[4][0]);
  v2[1] = (((float)vy)+m[4][1]);
  v2[2] = (((float)vz)+m[4][2]);

}

void transform3d3f ( oMatrix3d m1,float *v1,float *v2)
{
  int b;
  for(b=0;b<3;b++)
	 v2[b] = m1[b][0]*v1[0] +
		m1[b][1]*v1[1] + m1[b][2]*v1[2];
}

void transform33d3f ( Matrix33d m1,float *v1,float *v2)
{
  int b;
  for(b=0;b<3;b++)
	 v2[b] = (float)(m1[b][0]*v1[0] +
		m1[b][1]*v1[1] + m1[b][2]*v1[2]);
}

/*

void matcopy  ( oMatrix5f to,oMatrix5f from )
{
  int a,b;
  for(a=0;a<5;a++)
	 for(b=0;b<3;b++)
		to[a][b] = from[a][b];
}

void mattran ( oMatrix5f nm, oMatrix5f om, int axis, float dist )
{
  matcopy(nm,om);
  nm[4][axis] += dist;
}

void matrot ( oMatrix5f nm, oMatrix5f om, int axis, float angle )
{
  oMatrix5f rm;

  float ca,sa;
  int a,b;

  ca = cos(angle);
  sa = sin(angle);

  switch(axis)
	 {
	 case 0:
		rm[0][0] = _1;		rm[0][1] = _0;		rm[0][2] = _0;
		rm[1][0] = _0;		rm[1][1] =  ca;		rm[1][2] =  sa;
		rm[2][0] = _0;		rm[2][1] = -sa;		rm[2][2] =  ca;
		break;
	 case 1:
		rm[0][0] =  ca;		rm[0][1] = _0;		rm[0][2] = -sa;
		rm[1][0] = _0;		rm[1][1] = _1;		rm[1][2] = _0;
		rm[2][0] =  sa;		rm[2][1] = _0;		rm[2][2] =  ca;
		break;
	 case 2:
		rm[0][0] =  ca;		rm[0][1] =  sa;		rm[0][2] = _0;
		rm[1][0] = -sa;		rm[1][1] =  ca;		rm[1][2] = _0;
		rm[2][0] = _0;		rm[2][1] = _0;		rm[2][2] = _1;
		break;
	 }
  for(a=0;a<3;a++)
	 {
		nm[3][a] = om[3][a];
		nm[4][a] = om[4][a];
		for(b=0;b<3;b++)
		  nm[a][b] = 
			 rm[a][0]*om[0][b] + 
			 rm[a][1]*om[1][b] +
			 rm[a][2]*om[2][b];
	 }
  
  normalize3f(nm[0]);
  normalize3f(nm[1]);
  normalize3f(nm[2]);

}
*/

void rotation_to_matrix(Matrix53f rot,float *axis, float angle)
{
  rotation_matrix3f(angle,axis[0],axis[1],axis[2],&rot[0][0]);
}



static void find_axis( Matrix33d a, float *axis)
{
  doublereal at[3][3],v[3][3],vt[3][3],fv1[3][3];
  integer iv1[3];
  integer ierr;
  integer nm,n,matz;
  doublereal wr[3],wi[3];
  /*p[3][3];*/
  int x,y;

  nm = 3;
  n = 3;
  matz = 1;

  recondition33d(&a[0][0]); /* IMPORTANT! */

  for(x=0;x<3;x++)
	 {
		for(y=0;y<3;y++)
		  {
			 at[y][x] = a[x][y];
		  }
	 } 

  pymol_rg_(&nm,&n,&at[0][0],wr,wi,&matz,&vt[0][0],iv1,&fv1[0][0],&ierr);

  for(x=0;x<3;x++)
	 {
		for(y=0;y<3;y++)
		  {
			 v[y][x] = vt[x][y];
		  }
	 } 

  axis[0]=0.0F;
  axis[1]=0.0F;
  axis[2]=0.0F;

  {
    doublereal max_real = 0.0F, test_real;
    doublereal min_imag = 1.0F, test_imag;
    float test_inp[3],test_out[3];

    for(x=0;x<3;x++) { /* looking for an eigvalue of (1,0) */
      /*      printf("wr %8.3f wi %8.3f\n",wr[x],wi[x]);
      printf("%8.3f %8.3f %8.3f\n",
      v[0][x],v[1][x],v[2][x]);*/
      test_real = fabs(wr[x]);
      test_imag = fabs(wi[x]);
      
      if((test_real>=max_real)&&(test_imag<=min_imag)) {
        for(y=0;y<3;y++)
          test_inp[y] = (float)v[y][x];
        transform33d3f(a,test_inp,test_out); /* confirm that axis is invariant to rotation */
        test_out[0] -= test_inp[0];
        test_out[1] -= test_inp[1];
        test_out[2] -= test_inp[2];
        if((test_out[0]*test_out[0]+
            test_out[1]*test_out[1]+
            test_out[2]*test_out[2])<0.1) {
          for(y=0;y<3;y++)
            axis[y] = test_inp[y];
          max_real = test_real;
          min_imag = test_imag;
        }
      } else {
        /*for(y=0;y<3;y++)
          v[y][x]=_0;*/
      }
    }
  }
  /*
    printf("eigenvectors\n%8.3f %8.3f %8.3f\n",v[0][0],v[0][1],v[0][2]);
    printf("%8.3f %8.3f %8.3f\n",v[1][0],v[1][1],v[1][2]);
    printf("%8.3f %8.3f %8.3f\n",v[2][0],v[2][1],v[2][2]);
    
    
    printf("eigenvalues\n%8.3f %8.3f %8.3f\n",wr[0],wr[1],wr[2]);
    printf("%8.3f %8.3f %8.3f\n",wi[0],wi[1],wi[2]);
  */

    /*    matrix_multiply33d33d(a,v,p);
    
    printf("invariance\n");
    printf("%8.3f %8.3f %8.3f\n",p[0][0],p[0][1],p[0][2]);
    printf("%8.3f %8.3f %8.3f\n",p[1][0],p[1][1],p[1][2]);
    printf("%8.3f %8.3f %8.3f\n",p[2][0],p[2][1],p[2][2]);
    */

}

void matrix_to_rotation(Matrix53f rot,float *axis, float *angle)
{
  float perp[3],tmp[3],rperp[3],dirck[3];
  Matrix33d rot3d;
  Matrix53f rotck;
  int a,b;


#ifdef MATCHK
  printf("starting matrix\n");
  for(a=0;a<3;a++)
	 printf("%8.3f %8.3f %8.3f\n",rot[a][0],rot[a][1],rot[a][2]);
#endif

  for(a=0;a<3;a++)
	 for(b=0;b<3;b++)
		rot3d[a][b] = (double)rot[a][b];

  find_axis(rot3d,axis);


  /* find a perpendicular vector */

  perp[0]=axis[1]*axis[0]-axis[2]*axis[2];
  perp[1]=axis[2]*axis[1]-axis[0]*axis[0];
  perp[2]=axis[0]*axis[2]-axis[1]*axis[1];

  if(length3f(perp)<R_SMALL) {
		tmp[0]=axis[0];
		tmp[1]=-2*axis[1];
		tmp[2]=axis[2];
		cross_product3f(axis,tmp,perp);
  }

  normalize3f(perp);

  transform33d3f(rot3d,perp,rperp);

  *angle = get_angle3f(perp,rperp);

  cross_product3f(perp,rperp,dirck);
  if(((dirck[0]*axis[0])+(dirck[1]*axis[1])+(dirck[2]*axis[2]))<_0)
	 *angle = -*angle;
	 
  /*  printf("angle %8.3f \n",*angle);*/

  rotation_to_matrix(rotck,axis,*angle);


#ifdef MATCHK  
  printf("reconstructed matrix: \n");
  for(a=0;a<3;a++)
	 printf("%8.3f %8.3f %8.3f\n",rotck[a][0],rotck[a][1],rotck[a][2]);
  printf("\n");
#endif

}



