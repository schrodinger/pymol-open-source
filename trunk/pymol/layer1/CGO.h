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
#ifndef _H_CGO
#define _H_CGO

#include"Ray.h"

/* Compiled Graphics Library for simple graphics objects
   in floating point three-space, with the goal of achieving
   quick and easy rendering in multiple environments without the
   headaches of OpenGL arrays.

*/

/* Supported functions:
 * stop
 * null
 * begin
     GL_POINTS, 
     GL_LINES, GL_LINE_LOOP, GL_LINE_STRIP,
     GL_TRIANGLE, GL_TRIANGLE_STRIP, GL_TRIANGLE_FAN
 * end
 * vertex 
 * normal 
 * color 
 * sphere   * currently for ray-tracing only
 * triangle * currently for ray-tracing only
 * cylinder * currently for ray-tracing only
 * linewidth
 * primwidth * ray-tracing 
 */

typedef struct {
  float *op;
  int c;
} CGO;

/* instructions and data segment sizes */

#define CGO_STOP                 0x00
#define CGO_STOP_SZ              0
#define CGO_NULL                 0x01
#define CGO_NULL_SZ              0
#define CGO_BEGIN                0x02
#define CGO_BEGIN_SZ             1
#define CGO_END                  0x03
#define CGO_END_SZ               0
#define CGO_VERTEX               0x04
#define CGO_VERTEX_SZ            3
#define CGO_NORMAL               0x05
#define CGO_NORMAL_SZ            3
#define CGO_COLOR                0x06
#define CGO_COLOR_SZ             3
#define CGO_SPHERE               0x07
#define CGO_SPHERE_SZ            4
#define CGO_TRIANGLE             0x08
#define CGO_TRIANGLE_SZ          27
#define CGO_CYLINDER             0x09
#define CGO_CYLINDER_SZ          13
#define CGO_LINEWIDTH            0x10
#define CGO_LINEWIDTH_SZ         1
#define CGO_WIDTHSCALE           0x11
#define CGO_WIDTHSCALE_SZ        1

#define CGO_MASK                 0x0F
               
CGO *CGONew(void);
CGO *CGONewSized(int size);

void CGOFree(CGO *I);
CGO *CGOSimplify(CGO *I,int est);

void CGOReserve(CGO *ptr,int est);

int CGOCheckComplex(CGO *I);

int CGOFromFloatArray(CGO *I,float *src,int len);

void CGOBegin(CGO *I,int mode);
void CGOEnd(CGO *I);

void CGOVertex(CGO *I,float v1,float v2,float v3);
void CGOVertexv(CGO *I,float *v);
void CGOColor(CGO *I,float v1,float v2,float v3);
void CGOColorv(CGO *I,float *v);
void CGONormal(CGO *I,float v1,float v2,float v3);
void CGONormalv(CGO *I,float *v);

void CGOStop(CGO *I);

void CGORenderGL(CGO *I);
void CGORenderRay(CGO *I,CRay *ray);

#endif
