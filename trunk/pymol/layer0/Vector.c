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

#include"os_std.h"

#include"Base.h"
#include"Vector.h"

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

float sqrt1f(float f) { /* no good as a macro because f is used twice */
  if(f>0.0)
	 return(sqrt(f));
  else
	 return(0.0);
}

void dump3f( float *v, char *prefix ) /* for debugging */
{
  printf("%s %8.3f %8.3f %8.3f\n",prefix,v[0],v[1],v[2]);
}

void dump33f( float *m, char *prefix ) /* for debugging */
{
  printf("%s:0 %8.3f %8.3f %8.3f\n",prefix,m[0],m[1],m[2]);
  printf("%s:1 %8.3f %8.3f %8.3f\n",prefix,m[3],m[4],m[5]);
  printf("%s:2 %8.3f %8.3f %8.3f\n",prefix,m[6],m[7],m[8]);
}

void get_divergent3f(float *src,float *dst)
{
  if(src[0]!=0.0) {
    *(dst++)=-*(src++);
    *(dst++)= *(src++)+0.1;
    *(dst++)= *(src++);
  } else if(src[1]!=0.0) {
    *(dst++)= *(src++)+0.1;
    *(dst++)=-*(src++);
    *(dst++)= *(src++);
  } else {
    *(dst++)= *(src++)+0.1;
    *(dst++)= *(src++);
    *(dst++)=-*(src++);
  }
}

void get_random3f(float *x)
{
  x[0]=0.5-(rand()/(1.0+RAND_MAX));
  x[1]=0.5-(rand()/(1.0+RAND_MAX));
  x[2]=0.5-(rand()/(1.0+RAND_MAX));
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

#ifndef USE_VECTOR_MACROS
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

float lengthsq3f ( float *v1 )
{
  return((v1[0]*v1[0]) + 
					 (v1[1]*v1[1]) + 
					 (v1[2]*v1[2]));
} 

float length3f ( float *v1 )
{
  return(sqrt1f((v1[0]*v1[0]) + 
					 (v1[1]*v1[1]) + 
					 (v1[2]*v1[2])));
} 

void average3f ( float *v1, float *v2, float *avg )
{
  avg[0] = (v1[0]+v2[0])/2;
  avg[1] = (v1[1]+v2[1])/2;
  avg[2] = (v1[2]+v2[2])/2;
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


float project3f ( float *v1, float *v2, float *proj )
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

void identity44f ( float *m1 )
{
  int a;
  for(a=0;a<16;a++) m1[a]=0;
  for(a=0;a<16;a=a+5) m1[a]=1.0;
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
void transform33f3f (float *m1, float *m2, float *m3) 
{
  float m2r0=m2[0];
  float m2r1=m2[1];
  float m2r2=m2[2];
  m3[0] = m1[ 0] * m2r0 + m1[ 1] * m2r1 + m1[ 2] * m2r2;
  m3[1] = m1[ 3] * m2r0 + m1[ 4] * m2r1 + m1[ 5] * m2r2;
  m3[2] = m1[ 6] * m2r0 + m1[ 7] * m2r1 + m1[ 8] * m2r2;
}

void transform44f3f (float *m1, float *m2, float *m3)
{
  float m2r0 = m2[0];
  float m2r1 = m2[1];
  float m2r2 = m2[2];
  m3[0] = m1[ 0] * m2r0 + m1[ 1] * m2r1 + m1[ 2] * m2r2 + m1[ 3];
  m3[1] = m1[ 4] * m2r0 + m1[ 5] * m2r1 + m1[ 6] * m2r2 + m1[ 7];
  m3[2] = m1[ 8] * m2r0 + m1[ 9] * m2r1 + m1[10] * m2r2 + m1[11];
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

void transformTTT44f3f (float *m1, float *m2, float *m3)
{
  float m2r0 = m2[0] + m1[12];
  float m2r1 = m2[1] + m1[13];
  float m2r2 = m2[2] + m1[14];
  m3[0] = m1[ 0] * m2r0 + m1[ 1] * m2r1 + m1[ 2] * m2r2 + m1[3];
  m3[1] = m1[ 4] * m2r0 + m1[ 5] * m2r1 + m1[ 6] * m2r2 + m1[7];
  m3[2] = m1[ 8] * m2r0 + m1[ 9] * m2r1 + m1[10] * m2r2 + m1[11];
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
  return((angle* cPI)/180.0);
}
float rad_to_deg(float angle)
{
  return(180.0 * (angle/ cPI));
}

void rotation_to_matrix33f(float *axis, float angle, Matrix33f mat)
{
  rotation_matrix3f(angle,axis[0],axis[1],axis[2],&mat[0][0]);
}

void rotation_matrix3f( float angle, float x, float y, float z,float *m )
{
  int a,b;

   /* This function contributed by Erich Boleyn (erich@uruk.org) */
   float mag, s, c;
   float xx, yy, zz, xy, yz, zx, xs, ys, zs, one_c;

   s = sin(angle);
   c = cos(angle);

   mag = sqrt1f(x*x + y*y + z*z);
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
	  one_c = 1.0F - c;
	  
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
	  M(0,0)=1.0;
	  M(1,1)=1.0;
	  M(2,2)=1.0;
	}

}


#define get_angle USED_TO_RETURN_DEGREES

float get_dihedral3f( float *v0, float *v1, float *v2, float *v3 )
{
  Vector3f d01,d21,d32,dd1,dd3,pos_d;
  float result=0.0;

  subtract3f(v2,v1,d21);
  subtract3f(v0,v1,d01);
  subtract3f(v3,v2,d32);
  if (length3f(d21)<R_SMALL) {
    result = get_angle3f(d01,d32);
  } else { 
    cross_product3f(d21,d01,dd1);
    cross_product3f(d21,d32,dd3);
    if ((length3f(dd1)<R_SMALL)||(length3f(dd3)<R_SMALL)) /* degenerate cases */
      result = get_angle3f(d01,d32); /* fall back to angle between vectors */
    else {
      result = get_angle3f(dd1,dd3);
      cross_product3f(d21,dd1,pos_d);
      if(dot_product3f(dd3,pos_d)<0.0)
        result = -result;
    }
  }
  return(result);
}

float get_angle3f( float *v1, float *v2 )
{
  double denom;
  float result;

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
    result = 0.0;
  
  result = acos(result);
  
  return(result);
} 

void normalize23f( float *v1 , float *v2)
{
  float vlen;
  vlen = length3f(v1);
  if(vlen>R_SMALL)
	 {
		v2[0]=v1[0]/vlen;
		v2[1]=v1[1]/vlen;
		v2[2]=v1[2]/vlen;
	 }
  else
	 {
		v2[0]=0.0;
		v2[1]=0.0;
		v2[2]=0.0;
	 }
}

void normalize3f( float *v1 )
{
  float vlen;
  vlen = length3f(v1);
  if(vlen>R_SMALL)
	 {
		v1[0]/=vlen;
		v1[1]/=vlen;
		v1[2]/=vlen;
	 }
  else
	 {
		v1[0]=0.0;
		v1[1]=0.0;
		v1[2]=0.0;
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
		v1[0]=0.0;
		v1[1]=0.0;
		v1[2]=0.0;
	 }
} 


double length3d ( double *v1 )
{
  return(sqrt1f((v1[0]*v1[0]) + 
					 (v1[1]*v1[1]) + 
					 (v1[2]*v1[2])));
} 

double diff3f ( float *v1, float *v2 )
{
  register float dx,dy,dz;
  dx = (v1[0]-v2[0]);
  dy = (v1[1]-v2[1]);
  dz = (v1[2]-v2[2]);
  return(sqrt1f(dx*dx + dy*dy + dz*dz));
}

double diffsq3f ( float *v1, float *v2 )
{
  register double dx,dy,dz;
  dx = (v1[0]-v2[0]);
  dy = (v1[1]-v2[1]);
  dz = (v1[2]-v2[2]);

  return( dx*dx + dy*dy + dz*dz );

}


int within3f(float *v1,float *v2,float dist)
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


void remove_component3f ( float *v1, float *unit, float *result)
{
  float dot;

  dot = v1[0]*unit[0] + v1[1]*unit[1] + v1[2]*unit[2];
  result[0]=v1[0]-unit[0]*dot;
  result[1]=v1[1]-unit[1]*dot;
  result[2]=v1[2]-unit[2]*dot;  
}


float distance_line2point3f(float *base,float *normal,float *point,float *alongNormalSq)
{
  float hyp[3],adj[3];
  float result;

  hyp[0] = point[0] - base[0];
  hyp[1] = point[1] - base[1];
  hyp[2] = point[2] - base[2];

  project3f(hyp,normal,adj);

  (*alongNormalSq) = ((adj[0]*adj[0])+(adj[1]*adj[1])+(adj[2]*adj[2])); 
  result=((hyp[0]*hyp[0])+(hyp[1]*hyp[1])+(hyp[2]*hyp[2]))
	 - (*alongNormalSq);
  if(result<=0.0) 
	 return(0.0);
  else 
	 return(sqrt1f(result));

}

float distance_halfline2point3f(float *base,float *normal,float *point,float *alongNormalSq)
{
  float hyp[3],adj[3];
  float result;

  hyp[0] = point[0] - base[0];
  hyp[1] = point[1] - base[1];
  hyp[2] = point[2] - base[2];


  if(project3f(hyp,normal,adj)>0.0)
	 {
		(*alongNormalSq) = ((adj[0]*adj[0])+(adj[1]*adj[1])+(adj[2]*adj[2])); 
		result=((hyp[0]*hyp[0])+(hyp[1]*hyp[1])+(hyp[2]*hyp[2]))
		  - (*alongNormalSq);
		if(result<=0.0) 
		  return(0.0);
		else 
		  return(sqrt1f(result));
	 } else {
		return(MAXFLOAT);
	 }
}


void matrix_transform33f3f(Matrix33f m1,float *v1,float *v2)
{
  int b;
  for(b=0;b<3;b++)
	 v2[b] = m1[b][0]*v1[0] +
		m1[b][1]*v1[1] + m1[b][2]*v1[2];
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
void transform5f3f (oMatrix5f m, float *v1, float *v2 )
{
  register double dx,dy,dz;
  float vx,vy,vz;
  
  dx = v1[0] - m[3][0];
  dy = v1[1] - m[3][1];
  dz = v1[2] - m[3][2];

  vx = m[0][0] * dx + m[0][1] * dy + m[0][2] * dz;
  vy = m[1][0] * dx + m[1][1] * dy + m[1][2] * dz;
  vz = m[2][0] * dx + m[2][1] * dy + m[2][2] * dz;

  v2[0] = (vx+m[4][0]);
  v2[1] = (vy+m[4][1]);
  v2[2] = (vz+m[4][2]);

}

void transform3d3f ( oMatrix3d m1,float *v1,float *v2)
{
  int b;
  for(b=0;b<3;b++)
	 v2[b] = m1[b][0]*v1[0] +
		m1[b][1]*v1[1] + m1[b][2]*v1[2];
}


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
		rm[0][0] = 1.0;		rm[0][1] = 0.0;		rm[0][2] = 0.0;
		rm[1][0] = 0.0;		rm[1][1] =  ca;		rm[1][2] =  sa;
		rm[2][0] = 0.0;		rm[2][1] = -sa;		rm[2][2] =  ca;
		break;
	 case 1:
		rm[0][0] =  ca;		rm[0][1] = 0.0;		rm[0][2] = -sa;
		rm[1][0] = 0.0;		rm[1][1] = 1.0;		rm[1][2] = 0.0;
		rm[2][0] =  sa;		rm[2][1] = 0.0;		rm[2][2] =  ca;
		break;
	 case 2:
		rm[0][0] =  ca;		rm[0][1] =  sa;		rm[0][2] = 0.0;
		rm[1][0] = -sa;		rm[1][1] =  ca;		rm[1][2] = 0.0;
		rm[2][0] = 0.0;		rm[2][1] = 0.0;		rm[2][2] = 1.0;
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

void rotation_to_matrix(oMatrix5f rot,float *axis, float angle)
{
  rotation_matrix3f(angle,axis[0],axis[1],axis[2],&rot[0][0]);
}

void matrix_interpolate(oMatrix5f imat,oMatrix5f mat,float *pivot,
								float *axis,float angle,float tAngle,
								int linear,int tLinear,float fxn)
{
  int a;
  float pos[3],rotaxis[3],adj[3],adjdir[3],opp[3],oppdir[3];
  float p0[3],p1[3],center[3];
  float hyplen,adjlen,opplen;
  float tAlpha;

			 /*           ______--------______
           *        /____________          \
           *     /   \   opp     |adj         \
           *   |      \          |        trans   | 
           * (CM)---------------------------------->(CM)
           *    \       \        |               /
           *       \      \hyp   |             /
			  *          \p0   \    |adjdir   p1/
			  *             \    \  |        /
			  *                \  \ |     /     
			  *                   \\|  /
           *            <--------O pivot
			  *                      F-raxis
			  */


  for(a=0;a<3;a++)
	 {
		imat[a][0]=0.0;
		imat[a][1]=0.0;
		imat[a][2]=0.0;
		imat[a][a]=1.0;
	 }

  if(!linear) 
	 rotation_to_matrix(imat,axis,fxn*angle);
  else
	 tLinear = true;
  if(!tLinear)
	 {
		subtract3f(mat[3],pivot,p0);
		subtract3f(mat[4],pivot,p1);

		/*		printf("length match? %8.3f %8.3f\n",length3f(p0),length3f(p1));*/
		hyplen = length3f(p0);
		
		average3f(mat[3],mat[4],center);
		
		cross_product3f(p1,p0,rotaxis);
		normalize3f(rotaxis);
		subtract3f(center,pivot,adjdir);
		normalize3f(adjdir);
		cross_product3f(rotaxis,adjdir,oppdir);
		normalize3f(oppdir);
		
		tAlpha = fabs(0.5-fxn)*tAngle;
		
		opplen = fabs(hyplen * sin(deg_to_rad(tAlpha)));
		adjlen = fabs(hyplen * cos(deg_to_rad(tAlpha)));
		
		scale3f(oppdir,opplen,opp);
		scale3f(adjdir,adjlen,adj);
		
		add3f(pivot,adj,pos);
		
		if(fxn<=0.5)
		  add3f(pos,opp,pos);
		else
		  subtract3f(pos,opp,pos);
	 }
  else
	 {
		for(a=0;a<3;a++)
		  pos[a] = (mat[3][a]*(1.0-fxn))+(fxn*mat[4][a]);
	 }

  for(a=0;a<3;a++)
	 {
		imat[3][a] = mat[3][a];
		imat[4][a] = pos[a];
	 }

}

#ifdef FIND_AXIS
void matrix_to_rotation(oMatrix5f rot,float *axis, float *angle)
{
  float perp[3],tmp[3],rperp[3],dirck[3];
  oMatrix3d rot3d;
  oMatrix5f rotck;
  int a,b;
  
#ifdef MATCHK
  printf("starting matrix\n");
  for(a=0;a<5;a++)
	 printf("%8.3f %8.3f %8.3f\n",rot[a][0],rot[a][1],rot[a][2]);
#endif

  for(a=0;a<3;a++)
	 for(b=0;b<3;b++)
		rot3d[a][b] = rot[a][b];

  find_axis(rot3d,axis);

  /* find a perpendicular vector */

  perp[0]=axis[1]*axis[0]-axis[2]*axis[2];
  perp[1]=axis[2]*axis[1]-axis[0]*axis[0];
  perp[2]=axis[0]*axis[2]-axis[1]*axis[1];

  if(length3f(perp)<R_SMALL)
	 {
		tmp[0]=axis[0];
		tmp[1]=-2*axis[1];
		tmp[2]=axis[2];
		cross_product3f(axis,tmp,perp);
	 }

  normalize3f(perp);

  transform3d3f(rot3d,perp,rperp);

  *angle = get_angle3f(perp,rperp);

  cross_product3f(perp,rperp,dirck);
  if(((dirck[0]*axis[0])+(dirck[1]*axis[1])+(dirck[2]*axis[2]))<0.0)
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
#endif



#ifdef NOT_USED
void find_axis( oMatrix3d a, float *axis)
{
  double at[3][3],v[3][3],vt[3][3],fv1[3][3];
  int iv1[3];
  int ierr;
  int nm,n,matz;
  double wr[3],wi[3];
  /*p[3][3];*/
  int x,y;

  nm = 3;
  n = 3;
  matz = 1;

  for(x=0;x<3;x++)
	 {
		for(y=0;y<3;y++)
		  {
			 at[y][x] = a[x][y];
		  }
	 } 
 
  /*  rg_(&nm,&n,&at[0][0],wr,wi,&matz,&vt[0][0],iv1,&fv1[0][0],&ierr);*/

  for(x=0;x<3;x++)
	 {
		for(y=0;y<3;y++)
		  {
			 v[y][x] = vt[x][y];
		  }
	 } 

  for(x=0;x<3;x++)
	 {
		if((fabs(wr[x]-1.0)<R_MED)&&
			(fabs(wi[x])<R_SMALL))
		  for(y=0;y<3;y++)
			 axis[y] = v[y][x];
		else
		  for(y=0;y<3;y++)
			 v[y][x]=0.0;
	 }

  /*    printf("eigenvectors\n%8.3f %8.3f %8.3f\n",v[0][0],v[0][1],v[0][2]);
  printf("%8.3f %8.3f %8.3f\n",v[1][0],v[1][1],v[1][2]);
  printf("%8.3f %8.3f %8.3f\n",v[2][0],v[2][1],v[2][2]);


  printf("eigenvalues\n%8.3f %8.3f %8.3f\n",wr[0],wr[1],wr[2]);
  printf("%8.3f %8.3f %8.3f\n",wi[0],wi[1],wi[2]);
  
  matrix_multiply33d33d(a,v,p);

  printf("invariance\n");
  printf("%8.3f %8.3f %8.3f\n",p[0][0],p[0][1],p[0][2]);
  printf("%8.3f %8.3f %8.3f\n",p[1][0],p[1][1],p[1][2]);
  printf("%8.3f %8.3f %8.3f\n",p[2][0],p[2][1],p[2][2]);
  */
}
#endif
#endif


