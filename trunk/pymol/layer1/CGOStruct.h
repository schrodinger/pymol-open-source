
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
#ifndef _H_CGOStruct
#define _H_CGOStruct

#include"Base.h"

struct _CGO {
  PyMOLGlobals *G;
  float *op;
  int c;
  int z_flag;
  float z_min, z_max;
  float z_vector[3];
  float alpha;
  int *i_start, i_size;
  short has_begin_end;
  int current_pick_color_index, current_pick_color_bond;
  float current_accessibility;
  short has_draw_buffers, has_draw_cylinder_buffers, has_draw_sphere_buffers;
  float normal[3], color[3], texture[2];
  uchar pickColor[4];
  short use_shader, cgo_shader_ub_color, cgo_shader_ub_normal;
  short debug;
  short enable_shaders;
  short no_pick;
};

#endif
