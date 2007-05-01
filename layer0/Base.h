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
#ifndef _H_Base
#define _H_Base

#include"os_limits.h"
#include"os_types.h"

#ifndef true
#define true 1
#endif

#ifndef false
#define false 0
#endif

#ifndef NULL
#define NULL ((void*)0)
#endif

#ifndef PI
#define PI 3.14159265358979323846
#endif


#ifndef uchar
#define uchar unsigned char
#endif

#ifndef uint
#define uint unsigned int
#endif

#define MAX_VDW 2.5F  /* this has to go */


#ifndef MAXFLOAT
#define MAXFLOAT FLT_MAX
#endif

#ifndef R_SMALL4
#define R_SMALL4 0.0001F
#endif

#ifndef R_SMALL8
#define R_SMALL8 0.00000001F
#endif

typedef struct { 
  int index; /* atom index.
                NOTE: that first record contains the list count...not pick info */
  int bond; /* bond index, 
               >=0 for bond
               -1 for atom
               -2 for label
               -3 for gadget */
} Pickable;

#define cPickableAtom -1
#define cPickableLabel -2
#define cPickableGadget -3

typedef struct {
  void *object;
  int state;
  /*  int instance; */  /* to come... */
} PickContext;

typedef struct {
  Pickable src;
  PickContext context;
} Picking;

typedef struct {
  int mode;
  int x,y,w,h;
  Picking *picked;
} Multipick;

typedef struct {
  int mode;
  float pos[3];
  float offset[3];
} LabPosType;

/* not a global, but CRay widely used and Ray.h definitely isn't a
 * lightweight include... */

typedef struct _CRay               CRay;

typedef struct {
  int state;
  CRay *ray;
  Picking **pick;
  int pass,slot;
  int width_scale_flag;
  float front, back, stereo_front;
  float fog_start, fog_end;
  float view_normal[3];
  float width_scale;
  float vertex_scale; /* how large is a screen pixel in model space at the origin */
  float *pmv_matrix;
  int sampling; /* are we supersampling? */
} RenderInfo;

#define MAXLINELEN 1024

#define PYMOL_MAX_THREADS 32

#ifndef _PYMOL_NO_XRAY
#define _PYMOL_XRAY
#endif

#endif
