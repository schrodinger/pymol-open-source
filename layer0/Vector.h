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
#ifndef _H_Vector
#define _H_Vector

/* NOTE THIS VERSION USES RADIANS BY DEFAULT! */

#define cPI            3.14159265358979323846  /* pi */

typedef float Vector3f[3]; /* for local vars only - use float* for parameters */
typedef int Vector3i[3];
		  
typedef float Matrix33f[3][3]; 
typedef double Matrix33d[3][3]; 

float deg_to_rad(float angle);
float rad_to_deg(float angle);

void normalize3f( float *v1 );
void normalize3d( double *v1 );
void copy3f( float *v1 , float *v2);
void add3f ( float *v1, float *v2, float *v3 );
void subtract3f ( float *v1, float *v2, float *v3 );

#define average3f(a,b,c) average(a,b,c)

void average3f ( float *v1, float *v2, float *avg );
void scale3f ( float *v1, float v0, float *v2);
float length3f ( float *v1 );


double length3d ( double *v1 );
float lengthsq3f ( float *v1 );
float get_angle3f( float *v1, float *v2 );
double diffsq3f ( float *v1, float *v2 );
double diff3f ( float *v1, float *v2 );


void cross_product3f ( float *v1, float *v2, float *cross );
float dot_product3f ( float *v1, float *v2 );
double dot_product3d ( double *v1, double *v2 );
float project3f ( float *v1, float *v2, float *proj );
void remove_component3f ( float *v1, float *unit, float *result);

float distance_line2point3f(float *base,float *normal,float *point,float *alongNormalSq);
float distance_halfline2point3f(float *base,float *normal,float *point,float *alongNormalSq);

int within3f(float *v1,float *v2,float dist);

/* REVISED Matrix Routines */

void transform33f3f ( Matrix33f m1,float *v1,float *v2);
void rotation_to_matrix33f(float *axis, float angle, Matrix33f mat);
void multiply3d3d ( Matrix33d m1,Matrix33d m2,Matrix33d m3);

/* OLD MATRIX STUFF NEEDS REWORKING */

typedef float *oMatrix5f[5]; /* PHASE THESE OUT! - THEY CAUSE PROBLEMS! */

typedef float *oMatrix3f[3];

typedef float *oMatrix3d[3];

double matdiffsq ( float *v1, oMatrix5f m, float *v2 );
void find_axis( oMatrix3d a, float *axis);
void matcopy ( oMatrix5f to, oMatrix5f from );
void mattran ( oMatrix5f nm, oMatrix5f om, int axis, float dist );
void matrot ( oMatrix5f nm, oMatrix5f om, int axis, float angle );
void matrix_to_rotation(oMatrix5f rot,float *axis, float *angle);
void rotation_to_matrix(oMatrix5f rot,float *axis, float angle);

void matrix_interpolate(oMatrix5f imat,oMatrix5f mat,float *pivot,
								float *axis,float angle,float tAngle,
								int linear,int tLinear,float fxn);


void transform3d3f ( oMatrix3d m1,float *v1,float *v2);
void transform5f3f ( oMatrix5f m, float *v1, float *v2 );



#endif










