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

/* NOTE: All matrices are assumed to be column-major in this module */

#define cPI            3.14159265358979323846  /* pi */

typedef float Vector3f[3]; /* for local vars only - use float* for parameters */
typedef int Vector3i[3];
		  
typedef float Matrix33f[3][3]; 
typedef double Matrix33d[3][3]; 

unsigned int optimizer_workaround1u(unsigned int value);

float get_random0to1f(void);

float deg_to_rad(float angle);
float rad_to_deg(float angle);

double slow_sqrt1f(float f);
double slow_sqrt1d(double d);

void slow_normalize3f( float *v1 );
void normalize23f( float *v1 , float *v2);
void normalize3d( double *v1 );

void clamp3f(float *v1);
void get_divergent3f(float *src,float *dst);
void get_random3f(float *x);
void scatter3f(float *v,float weight);
void wiggle3f(float *v,float *p,float *s);

void get_system3f(float *x,float *y,float *z); /* make random system */
void get_system1f3f(float *x,float *y,float *z); /* make system in direction of x */
void get_system2f3f(float *x,float *y,float *z); /* make system in direction of x, perp to x,y */

double dot_product3d ( double *v1, double *v2 );
float slow_project3f ( float *v1, float *v2, float *proj );
void slow_remove_component3f ( float *v1, float *unit, float *result);

double distance_line2point3f(float *base,float *normal,float *point,float *alongNormalSq);
double distance_halfline2point3f(float *base,float *normal,float *point,float *alongNormalSq);

double slow_diffsq3f ( float *v1, float *v2 );
double slow_diff3f ( float *v1, float *v2 );
int slow_within3f(float *v1,float *v2,float dist);
int slow_within3fsq(float *v1,float *v2,float dist,float dist2);

int equal3f(float *v1,float *v2);

float get_angle3f( float *v1, float *v2 );
float get_dihedral3f( float *v0, float *v1, float *v2, float *v3 );
double length3d ( double *v1 );

void min3f ( float *v1, float *v2, float *v3 );
void max3f ( float *v1, float *v2, float *v3 );

void dump3i( int *v, char *prefix );
void dump3f( float *v, char *prefix );
void dump33f( float *m, char *prefix );
void dump44f( float *m, char *prefix );

void copy44f( float *src,float *dst);

/* REVISED Matrix Routines using  float pointers */

void identity44f ( float *m1 );

void copy44f44f ( float *src, float *dst );

/* in the following matrix multiplies and transformations:
   the last two matrices can be the same matrix! */

void transform33f3f ( float  *m1, float  *m2,  float  *m3 );
void transform33Tf3f ( float  *m1, float  *m2,  float  *m3 ); /* uses transpose */

void transform44f3f ( float  *m1, float  *m2,  float  *m3 );
void transform44f4f ( float  *m1, float  *m2,  float  *m3 );

void multiply33f33f ( float  *m1, float  *m2,  float  *m3 );
void multiply33d33d ( double *m1, double *m2, double  *m3 );

/* as matrix types */

void matrix_transform33f3f ( Matrix33f m1,float *v1,float *v2);
void rotation_to_matrix33f (float *axis, float angle, Matrix33f mat);
void matrix_multiply33f33f ( Matrix33f m1,Matrix33f m2,Matrix33f m3);
void matrix_multiply33d33d ( Matrix33d m1,Matrix33d m2,Matrix33d m3);

/* A 4x4 TTT matrix is really a 3x3 rotation matrix with two translation vectors:
   (1) a pre-translation stored in forth row, first three columns.
   (2) and a post-translation stored in forth column, first three rows.
   There are certain cases where this representation is more convenient.
 */
void combineTTT44f44f( float *m1, float *m2, float *m3);
void transformTTT44f3f ( float *m1, float *m2, float *m3 );
void transform_normalTTT44f3f ( float *m1, float *m2, float *m3 );
void initializeTTT44f ( float *m );
/* end revised matrix routines */

/*------------------------------------------------------------------------*/
/* OLD MATRIX STUFF below NEEDS REWORKING */

void rotation_matrix3f( float angle, float x, float y, float z,float *m );

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

/* macros */

#define USE_VECTOR_MACROS

#ifndef USE_VECTOR_MACROS

float dot_product3f ( float *v1, float *v2 );
void  invert3f ( float *v );
void  scale3f ( float *v1, float v0, float *v2);
void  copy3f( float *src , float *dst);
void  copy4f( float *src , float *dst);
void  add3f ( float *v1, float *v2, float *sum );
void  subtract3f ( float *v1, float *v2, float *v3 );
double lengthsq3f ( float *v1 );
double length3f ( float *v1 );
void  cross_product3f ( float *v1, float *v2, float *cross );
void  average3f ( float *v1, float *v2, float *avg );
void  zero3f ( float *v1 )
void  set3f ( float *v1,float x,float y,float z );
void  swap1f (float *f, float *g);

#else

#define set3f(v1,x,y,z) { (v1)[0]=(x);(v1)[1]=(y);(v1)[2]=(z); }
#define zero3f(v1) { (v1)[0]=0.0;(v1)[1]=0.0;(v1)[2]=0.0; }
#define dot_product3f(v1,v2) ((v1)[0]*(v2)[0] + (v1)[1]*(v2)[1] + (v1)[2]*(v2)[2])
#define invert3f(v) {(v)[0]=-(v)[0]; (v)[1]=-(v)[1]; (v)[2]=-(v)[2];}
#define scale3f(v1,v0,v2) {(v2)[0]=(v1)[0]*(v0); (v2)[1]=(v1)[1]*(v0); (v2)[2]=(v1)[2]*(v0);}
#define copy3f(v1,v2) {(v2)[0]=(v1)[0]; (v2)[1]=(v1)[1]; (v2)[2]=(v1)[2];}
#define copy4f(v1,v2) {(v2)[0]=(v1)[0]; (v2)[1]=(v1)[1]; (v2)[2]=(v1)[2]; (v2)[3]=(v1)[3];}
#define add3f(v1,v2,v3) {(v3)[0]=(v1)[0]+(v2)[0]; (v3)[1]=(v1)[1]+(v2)[1]; (v3)[2]=(v1)[2]+(v2)[2];}
#define subtract3f(v1,v2,v3) {(v3)[0]=(v1)[0]-(v2)[0]; (v3)[1]=(v1)[1]-(v2)[1]; (v3)[2]=(v1)[2]-(v2)[2];}
#define lengthsq3f(v1) (((v1)[0]*(v1)[0]) + ((v1)[1]*(v1)[1]) + ((v1)[2]*(v1)[2]))
#define length3f(v1) (sqrt1f(((v1)[0]*(v1)[0]) + ((v1)[1]*(v1)[1]) + ((v1)[2]*(v1)[2])))
#define average3f(v1,v2,avg) { \
  (avg)[0] = ((v1)[0]+(v2)[0])/2; \
  (avg)[1] = ((v1)[1]+(v2)[1])/2; \
  (avg)[2] = ((v1)[2]+(v2)[2])/2; \
}
#define cross_product3f(v1,v2,cross) { \
  (cross)[0] = ((v1)[1]*(v2)[2]) - ((v1)[2]*(v2)[1]); \
  (cross)[1] = ((v1)[2]*(v2)[0]) - ((v1)[0]*(v2)[2]); \
  (cross)[2] = ((v1)[0]*(v2)[1]) - ((v1)[1]*(v2)[0]); \
}
#define swap1f(f,g) { float h;h=*(f);*(f)=*(g);*(g)=h; }
#endif


#ifdef _PYMOL_INLINE

#ifdef _PYMOL_WIN32
#define __inline__ __inline
#endif

static const float _0f_inline = 0.0F;
static const double _0d_inline = 0.0;
static const float _1f_inline = 1.0F;
static const double _1d_inline = 1.0;
static const float R_SMALL_inline = 0.000000001F;
static const double R_SMALLd_inline = 0.000000001;

#define normalize3f inline_normalize3f
#define sqrt1f inline_sqrt1f
#define sqrt1d inline_sqrt1d
#define diff3f inline_diff3f
#define diffsq3f inline_diffsq3f
#define within3f inline_within3f
#define within3fsq inline_within3fsq
#define remove_component3f inline_remove_component3f
#define project3f inline_project3f

__inline__ static double inline_sqrt1f(float f) { /* no good as a macro because f is used twice */
  if(f>_0f_inline)
	 return(sqrt(f));
  else
	 return(_0d_inline);
}

__inline__ static double inline_sqrt1d(double f) { /* no good as a macro because f is used twice */
  if(f>_0d_inline)
	 return(sqrt(f));
  else
	 return(_0d_inline);
}

__inline__ static void inline_normalize3f( float *v1 )
{
	double vlen = length3f(v1);
	if(vlen > R_SMALLd_inline)
	{
		float	inV	= (float)(_1d_inline / vlen);
		v1[0] *= inV;
		v1[1] *= inV;
		v1[2] *= inV;
	}
	else
	{
		v1[0]=v1[1]=v1[2]=_0f_inline;
	}
}

__inline__ static double inline_diff3f ( float *v1, float *v2 )
{
  register float dx,dy,dz;
  dx = (v1[0]-v2[0]);
  dy = (v1[1]-v2[1]);
  dz = (v1[2]-v2[2]);
  return(sqrt1d(dx*dx + dy*dy + dz*dz));
}

__inline__ static double inline_diffsq3f ( float *v1, float *v2 )
{
  register double dx,dy,dz;
  dx = (v1[0]-v2[0]);
  dy = (v1[1]-v2[1]);
  dz = (v1[2]-v2[2]);

  return( dx*dx + dy*dy + dz*dz );

}

__inline__ static int inline_within3f(float *v1,float *v2,float dist)
{
  register float dx,dy,dz;
  dx = (v1[0]-v2[0]);
  if(fabs(dx)>dist) return(0);
  dy = (v1[1]-v2[1]);
  if(fabs(dy)>dist) return(0);
  dz = (v1[2]-v2[2]);
  if(fabs(dz)>dist) return(0);
  return((dx*dx + dy*dy + dz*dz)<=(dist*dist));
}

__inline__ static int inline_within3fsq(float *v1,float *v2,float dist,float dist2)
{
  register float dx,dy,dz;
  dx = (v1[0]-v2[0]);
  if(fabs(dx)>dist) return(0);
  dy = (v1[1]-v2[1]);
  if(fabs(dy)>dist) return(0);
  dz = (v1[2]-v2[2]);
  if(fabs(dz)>dist) return(0);
  return((dx*dx + dy*dy + dz*dz)<=(dist2));
}


__inline__ static void inline_remove_component3f ( float *v1, float *unit, float *result)
{
  float dot;

  dot = v1[0]*unit[0] + v1[1]*unit[1] + v1[2]*unit[2];
  result[0]=v1[0]-unit[0]*dot;
  result[1]=v1[1]-unit[1]*dot;
  result[2]=v1[2]-unit[2]*dot;  
}

__inline__ static float inline_project3f ( float *v1, float *v2, float *proj )
{
   float dot;

	dot = v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2];
	proj[0] = v2[0] * dot;
	proj[1] = v2[1] * dot;
	proj[2] = v2[2] * dot;
	
	return(dot);
}

#else

#define normalize3f slow_normalize3f
#define sqrt1f slow_sqrt1f
#define sqrt1d slow_sqrt1d
#define diff3f slow_diff3f
#define diffsq3f slow_diffsq3f
#define within3f slow_within3f
#define within3fsq slow_within3fsq
#define project3f slow_project3f
#define remove_component3f slow_remove_component3f

#endif

#endif










