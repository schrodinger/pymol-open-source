

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
#ifndef _H_Base
#define _H_Base

#include "os_limits.h"
#include "os_types.h"
#include "Picking.h"
#include <vector>

#ifndef PI
#define PI 3.14159265358979323846
#endif

typedef unsigned char uchar;

typedef unsigned int uint;

#define MAX_VDW 2.5F            /* this has to go */

#ifndef R_SMALL4
#define R_SMALL4 0.0001F
#endif

#ifndef R_SMALL8
#define R_SMALL8 0.00000001F
#endif

#define cPuttyTransformNormalizedNonlinear 0
#define cPuttyTransformRelativeNonlinear   1
#define cPuttyTransformScaledNonlinear     2
#define cPuttyTransformAbsoluteNonlinear   3

#define cPuttyTransformNormalizedLinear    4
#define cPuttyTransformRelativeLinear      5
#define cPuttyTransformScaledLinear        6
#define cPuttyTransformAbsoluteLinear      7

#define cPuttyTransformImpliedRMS          8

typedef struct LabPosType {
  int mode;
  float pos[3];
  float offset[3];
} LabPosType;

typedef struct RefPosType {
  float coord[3];
  int specified;
} RefPosType;


/* not a global, but CRay widely used and Ray.h definitely isn't a
 * lightweight include... */

typedef struct _CRay CRay;


/* likewise */

#ifndef CGO_DEFINED
class CGO;
#define CGO_DEFINED
#endif

struct RenderInfo {
  int state;
  CRay *ray;
  CGO *alpha_cgo;
  PickColorManager* pick = nullptr;
  int pass;
  int width_scale_flag;
  float front, back, stereo_front;
  float fog_start, fog_end;
  float view_normal[3];
  float width_scale;
  float vertex_scale;           /* how large is a screen pixel in model space at the origin */
  int sampling;                 /* are we supersampling? */
  int ortho;                    /* orthoscopic projection? */
  int line_lighting;            /* line lighting */
  int dynamic_width;
  float dynamic_width_factor, dynamic_width_min, dynamic_width_max;
  int texture_font_size;
  int use_shaders;
};

#define MAXLINELEN 1024

#define PYMOL_MAX_THREADS 125

#ifndef _PYMOL_NO_XRAY
#define _PYMOL_XRAY
#endif

// no MMLIBS in open-source
#ifndef NO_MMLIBS
#define NO_MMLIBS
#endif

#endif
