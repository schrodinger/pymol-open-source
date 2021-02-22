

/* 
A* -------------------------------------------------------------------
B* This file contains source code for the PyMOL computer program
C* Copyright (c) Schrodinger, LLC. 
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

#include"os_predef.h"
#include"os_gl.h"
#include<math.h>

#include <cassert>

/* NOTE THIS VERSION USES RADIANS BY DEFAULT! */


/* NOTE: Matrices are assumed to be row-major (C-like not
 * OpenGL-like) unless explictly labeled as per the following
 * conventions:
 * 
 * row-major:    33f, 33d, 44f, 44d, R33f, R33d, R44f, R44d
 * column-major: C33f, C33d, C44f, C44d
 */

#define cPI            3.14159265358979323846   /* pi */

#define GET_BIT(val, bit)       (((val) >> (bit)) & 1)
#define SET_BIT(val, bit)       (val) |= (1 << (bit))
#define SET_BIT_OFF(val, bit)   (val) &= ~(1 << (bit))
#define SET_BIT_TO(val, bit, v) {if(v) SET_BIT(val, bit); else SET_BIT_OFF(val, bit);}

short countBits(unsigned long bits);
short countBitsInt(int bits);

typedef float Vector3f[3];      /* for local vars only - use float* for parameters */
typedef float Vector4f[4];
typedef int Vector3i[3];

typedef float Matrix33f[3][3];
typedef double Matrix33d[3][3];
typedef float Matrix53f[5][3];
typedef double Matrix53d[5][3];

unsigned int optimizer_workaround1u(unsigned int value);

float get_random0to1f(void);

float deg_to_rad(float angle);
float rad_to_deg(float angle);

void normalize23f(const float *v1, float *v2);
void normalize2f(float *v1);
void normalize4f(float *v1);

void clamp3f(float *v1);
void get_divergent3f(const float *src, float *dst);
void get_random3f(float *x);
void scatter3f(float *v, float weight);
void wiggle3f(float *v, const float *p, const float *s);
void extrapolate3f(const float *v1, const float *unit, float *result);

void mix3f(const float *v1, const float *v2, float fxn, float *v3);
void mix3d(const double *v1, const double *v2, double fxn, double *v3);

void get_system3f(float *x, float *y, float *z);        /* make random system */
void get_system1f3f(float *x, float *y, float *z);      /* make system in direction of x */
void get_system2f3f(float *x, float *y, float *z);      /* make system in direction of x, perp to x,y */

double distance_line2point3f(const float *base, const float *normal, const float *point,
                             float *alongNormalSq);
double distance_halfline2point3f(const float *base, const float *normal, const float *point,
                                 float *alongNormalSq);

int equal3f(const float *v1, const float *v2);

int pymol_roundf(float f);

float get_angle3f(const float *v1, const float *v2);
float get_dihedral3f(const float *v0, const float *v1, const float *v2, const float *v3);

void min3f(const float *v1, const float *v2, float *v3);
void max3f(const float *v1, const float *v2, float *v3);

void dump3i(const int *v, const char *prefix);
void dump2f(const float *v, const char *prefix);
void dump3f(const float *v, const char *prefix);
void dump3d(const double *v, const char *prefix);
void dump4f(const float *v, const char *prefix);
void dump33f(const float *m, const char *prefix);
void dump33d(const double *m, const char *prefix);
void dump44f(const float *m, const char *prefix);
void dump44d(const double *m, const char *prefix);

void copy44f(const float *src, float *dst);
void copy44d(const double *src, double *dst);

void identity33f(float *m1);
void identity33d(double *m);
void identity44f(float *m1);
void identity44d(double *m1);

bool is_identityf(int n, const float *m, float threshold=1.0E-6F);
bool is_allclosef(int nrow,
    const float *m1, int ncol1,
    const float *m2, int ncol2, float threshold=1.0E-6F);
bool is_diagonalf(int nrow,
    const float *m, int ncol=0, float threshold=1.0E-6F);
double determinant33f(const float *m, int ncol=3);

void glOrtho44f(float *m1, 
		GLfloat left, GLfloat right, 
		GLfloat bottom, GLfloat top,
		GLfloat nearVal, GLfloat farVal);
void glFrustum44f(float *m1,
		  GLfloat left, GLfloat right, 
		  GLfloat bottom, GLfloat top,
		  GLfloat nearVal, GLfloat farVal);

void copy44f44f(const float *src, float *dst);
void copy44d44f(const double *src, float *dst);
void copy44f44d(const float *src, double *dst);

void copy44d33f(const double *src, float *dst);
void copy44f33f(const float *src, float *dst);
void copy33f44d(const float *src, double *dst);
void copy33f44f(const float *src, float *dst);
void copy3d3f(const double *v1, float *v2);
void copy3f3d(const float *v1, double *v2);


/* in the following matrix multiplies and transformations:
   the last two matrices can be the same matrix! */

void transpose33f33f(const float *m1, float *m2);
void transpose33d33d(const double *m1, double *m2);
void transpose44f44f(const float *m1, float *m2);
void transpose44d44d(const double *m1, double *m2);

void transform33f3f(const float *m1, const float *m2, float *m3);
void transform33Tf3f(const float *m1, const float *m2, float *m3);  /* uses transpose */

void transform44f3f(const float *m1, const float *m2, float *m3);
void transform44f4f(const float *m1, const float *m2, float *m3);

void transform44d3f(const double *m1, const float *m2, float *m3);
void transform44d3d(const double *m1, const double *m2, double *m3);
void inverse_transformC44f3f(const float *m1, const float *m2, float *m3);
void inverse_transform44f3f(const float *m1, const float *m2, float *m3);
void inverse_transform44d3f(const double *m1, const float *m2, float *m3);
void inverse_transform44d3d(const double *m1, const double *m2, double *m3);
void transform44f3fas33f3f(const float *m1, const float *m2, float *m3);
void transform44d3fas33d3f(const double *m1, const float *m2, float *m3);

void multiply33f33f(const float *m1, const float *m2, float *m3);
void multiply33d33d(const double *m1, const double *m2, double *m3);


/* as matrix types */

void matrix_transform33f3f(const Matrix33f m1, const float *v1, float *v2);
void matrix_inverse_transform33f3f(const Matrix33f m1, const float *v1, float *v2);

void rotation_to_matrix33f(const float *axis, float angle, Matrix33f mat);
void matrix_multiply33f33f(Matrix33f m1, Matrix33f m2, Matrix33f m3);
void matrix_multiply33d33d(Matrix33d m1, Matrix33d m2, Matrix33d m3);


/* A 4x4 TTT matrix is really a 3x3 rotation matrix with two translation vectors:
   (1) a pre-translation stored in forth row, first three columns.
   (2) and a post-translation stored in forth column, first three rows.
   There are certain cases where this representation is more convenient.
 */
void combineTTT44f44f(const float *m1, const float *m2, float *m3);
void transformTTT44f3f(const float *m1, const float *m2, float *m3);
void transform_normalTTT44f3f(const float *m1, const float *m2, float *m3);
void initializeTTT44f(float *m);

void multiply44d44d44d(const double *left, const double *right, double *product);
void left_multiply44d44d(const double *left, double *right);
void right_multiply44d44d(double *left, const double *right);

void multiply44f44f44f(const float *left, const float *right, float *product);
void left_multiply44f44f(const float *left, float *right);
void right_multiply44f44f(float *left, const float *right);

void reorient44d(double *matrix);

void recondition33d(double *matrix);
void recondition44d(double *matrix);


/* invert a 4x4 homogenous that contains just rotation & tranlation
  (e.g. no scaling & fourth row is 0,0,0,1) */
void invert_special44d44d(const double *original, double *inv);
void invert_special44f44f(const float *original, float *inv);

void invert_rotation_only44d44d(const double *original, double *inv);

void convertTTTfR44d(const float *ttt, double *homo);
void convertTTTfR44f(const float *ttt, float *homo);
void convertR44dTTTf(const double *homo, float *ttt);
void convert44d44f(const double *dbl, float *flt);
void convert44f44d(const float *flt, double *dbl);

void get_rotation_about3f3fTTTf(float angle, const float *dir, const float *origin, float *ttt);


/* end revised matrix routines */


/*------------------------------------------------------------------------*/


/* OLD MATRIX STUFF below NEEDS REWORKING */

void rotation_matrix3f(float angle, float x, float y, float z, float *m);

typedef float *oMatrix5f[5];    /* PHASE THESE OUT! - THEY CAUSE PROBLEMS! */

typedef float *oMatrix3f[3];

typedef float *oMatrix3d[3];


/*void matcopy ( oMatrix5f to, oMatrix5f from );
  void mattran ( oMatrix5f nm, oMatrix5f om, int axis, float dist );
  void matrot ( oMatrix5f nm, oMatrix5f om, int axis, float angle );*/

void matrix_to_rotation(Matrix53f rot, float *axis, float *angle);
void rotation_to_matrix(Matrix53f rot, const float *axis, float angle);

void transform3d3f(const oMatrix3d m1, const float *v1, float *v2);
void transform33d3f(const Matrix33d m1, const float *v1, float *v2);
void transform5f3f(const oMatrix5f m, const float *v1, float *v2);

void mult4f(const float *vsrc, float val, float *vdest);
void mult3f(const float *vsrc, float val, float *vdest);
float max3(float val1, float val2, float val3);
float ave3(float val1, float val2, float val3);
float ave2(float val1, float val2);
void white4f(float *rgba, float value);
void add4f(const float *v1, const float *v2, float *sum);

int countchrs(const char *str, char ch);

float smooth(float x, float power);

void subdivide(int n, float *x, float *y);

//-------------------------------------------------------------------------
// Small inline functions
// (many of these were macros up to PyMOL 1.7.6)
//-------------------------------------------------------------------------

inline void set3f(float * v1, float x, float y, float z) {
  v1[0] = x;
  v1[1] = y;
  v1[2] = z;
}

template <typename T>
inline void zero3(T * v1) {
  v1[0] = 0;
  v1[1] = 0;
  v1[2] = 0;
}

#define zero3f zero3
#define zero3i zero3

inline void ones3f(float * v1) {
  v1[0] = 1.0F;
  v1[1] = 1.0F;
  v1[2] = 1.0F;
}

/**
 * Sets the values of the second vector to the additive inverses of each
 * respective value of the first vector, leaving the first vector unchanged.
 * @param v1 an array of three constant floats
 * @param v2 an array of three constant floats
 */
inline void invert3f3f(const float * v1, float * v2) {
  v2[0] = -v1[0];
  v2[1] = -v1[1];
  v2[2] = -v1[2];
}

inline void invert3f(float * v) {
  invert3f3f(v, v);
}

template <typename S, typename D>
void copyN(const S * src, D * dst, int N) {
  for (int i = 0; i < N; ++i)
    dst[i] = src[i];
}

template <typename S, typename D>
void copy2(const S * src, D * dst) {
  dst[0] = src[0];
  dst[1] = src[1];
}

template <typename S, typename D>
void copy3(const S * src, D * dst) {
  dst[0] = src[0];
  dst[1] = src[1];
  dst[2] = src[2];
}

template <typename S, typename D>
void copy4(const S * src, D * dst) {
  dst[0] = src[0];
  dst[1] = src[1];
  dst[2] = src[2];
  dst[3] = src[3];
}

#define copy2f copy2<float, float>
#define copy3f copy3
#define copy3d copy3<double, double>
#define copy4f copy4

namespace pymol
{
//! Dot product of two 3-dimensional vectors
template <typename T> T dot_product3(const T* v1, const T* v2)
{
  return v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2];
}

//! Square root of v, or 0 for negative v.
template <typename T> T sqrt1(T v)
{
  return (v > 0) ? sqrt(v) : 0;
}

//! v3 = v1 + v2
template <typename S1, typename S2, typename D>
void add3(const S1* v1, const S2* v2, D* v3)
{
  v3[0] = v1[0] + v2[0];
  v3[1] = v1[1] + v2[1];
  v3[2] = v1[2] + v2[2];
}

//! v3 = v1 - v2
template <typename S1, typename S2, typename D>
void subtract3(const S1* v1, const S2* v2, D* v3)
{
  v3[0] = v1[0] - v2[0];
  v3[1] = v1[1] - v2[1];
  v3[2] = v1[2] - v2[2];
}

//! Squared length of a 3-dimensional vector
template <typename T> T lengthsq3(const T* v1)
{
  return (v1[0] * v1[0]) + (v1[1] * v1[1]) + (v1[2] * v1[2]);
}

//! Length of a 3-dimensional vector
template <typename T> T length3(const T* v1)
{
  return sqrt1(lengthsq3(v1));
}

template <typename S, typename F, typename D>
void scale3(const S* v1, F v0, D* v2)
{
  v2[0] = v1[0] * v0;
  v2[1] = v1[1] * v0;
  v2[2] = v1[2] * v0;
}

/**
 * Scale v1 to length 1, unless it's a null pointer, then make it length 0.
 */
template <typename T> void normalize3(T* v1)
{
  auto const vlen = length3(v1);
  if (vlen > 1e-8) {
    scale3(v1, 1 / vlen, v1);
  } else {
    v1[0] = v1[1] = v1[2] = 0;
  }
}

/**
 * Cross product of two 3-dimensional vectors.
 * @param[out] cross Result buffer, must not overlap with v1 or v2
 */
template <typename T> void cross_product3(const T* v1, const T* v2, T* cross)
{
  assert(v1 != cross);
  assert(v2 != cross);
  cross[0] = (v1[1] * v2[2]) - (v1[2] * v2[1]);
  cross[1] = (v1[2] * v2[0]) - (v1[0] * v2[2]);
  cross[2] = (v1[0] * v2[1]) - (v1[1] * v2[0]);
}

/**
 * Same as:
 *
 *     project3f(v1, unit, result);
 *     subtract3f(v1, result, result);
 */
template <typename T>
void remove_component3(const T *v1, const T *unit, T *result)
{
  auto dot = v1[0] * unit[0] + v1[1] * unit[1] + v1[2] * unit[2];
  result[0] = v1[0] - unit[0] * dot;
  result[1] = v1[1] - unit[1] * dot;
  result[2] = v1[2] - unit[2] * dot;
}

template <typename T = double, typename U>
T diffsq3(const U* const v1, const U* const v2)
{
  auto const dx = T(v1[0]) - v2[0];
  auto const dy = T(v1[1]) - v2[1];
  auto const dz = T(v1[2]) - v2[2];
  return dz * dz + dy * dy + dx * dx;
}

template <typename T = double, typename U>
T diff3(const U* const v1, const U* const v2)
{
  return sqrt1(diffsq3<T>(v1, v2));
}

void meanNx3(float const* data, size_t N, float* out);
} // namespace pymol

#define dot_product3f pymol::dot_product3<float>
#define dot_product3d pymol::dot_product3<double>
#define add3f pymol::add3<float, float, float>
#define add3d pymol::add3<double, double, double>
#define subtract3f pymol::subtract3<float, float, float>
#define lengthsq3f pymol::lengthsq3<float>
#define sqrt1f pymol::sqrt1<float>
#define sqrt1d pymol::sqrt1<double>
#define length3f pymol::length3<float>
#define length3d pymol::length3<double>
#define scale3f pymol::scale3<float, float, float>
#define scale3d pymol::scale3<double, double, double>
#define normalize3f pymol::normalize3<float>
#define normalize3d pymol::normalize3<double>
#define cross_product3f pymol::cross_product3<float>
#define cross_product3d pymol::cross_product3<double>
#define remove_component3f pymol::remove_component3<float>
#define diffsq3f pymol::diffsq3<float, float>
#define diff3f pymol::diff3<double, float>

inline void average3f(const float * v1, const float * v2, float * avg) {
  (avg)[0] = ((v1)[0]+(v2)[0])/2;
  (avg)[1] = ((v1)[1]+(v2)[1])/2;
  (avg)[2] = ((v1)[2]+(v2)[2])/2;
}

inline float length2f(const float * v1) {
  return sqrt1f((v1[0] * v1[0]) + (v1[1] * v1[1]));
}

inline bool within3f(const float *v1, const float *v2, float dist)
{
  float dx, dy, dz, dist2;
  dx = (float) fabs(v1[0] - v2[0]);
  dy = (float) fabs(v1[1] - v2[1]);
  if(dx > dist)
    return (0);
  dz = (float) fabs(v1[2] - v2[2]);
  dx = dx * dx;
  if(dy > dist)
    return (0);
  dy = dy * dy;
  dist2 = dist * dist;
  if(dz > dist)
    return (0);
  return (((dx + dy) + dz * dz) <= dist2);
}

inline bool within3fret(const float *v1, const float *v2, float cutoff,
                                         const float cutoff2, float *diff, float *dist)
{
  float dx, dy, dz, dist2;
  dx = (float) fabs((diff[0] = v1[0] - v2[0]));
  dy = (float) fabs((diff[1] = v1[1] - v2[1]));
  if(dx > cutoff)
    return 0;
  dz = (float) fabs((diff[2] = v1[2] - v2[2]));
  dx = dx * dx;
  if(dy > cutoff)
    return 0;
  dy = dy * dy;
  if(dz > cutoff)
    return 0;
  if((dist2 = ((dx + dy) + dz * dz)) > cutoff2)
    return 0;
  *dist = (float) sqrt1f(dist2);
  return 1;
}

/**
 * Get the shortest position vector for a plane with position v1 and normal vector v2.
 * @param v1 Point on plane
 * @param v2 Normal vector of plane
 * @param[out] proj Point on plane, collinear with v2
 * @return Length of proj
 */
inline float project3f(const float *v1, const float *v2, float *proj)
{
  float dot;

  dot = v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2];
  proj[0] = v2[0] * dot;
  proj[1] = v2[1] * dot;
  proj[2] = v2[2] * dot;

  return (dot);
}

#endif
