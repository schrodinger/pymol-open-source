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
#ifndef _H_ObjectSlice
#define _H_ObjectSlice

#include"ObjectMap.h"

typedef struct {
  int Active;
  char MapName[ObjNameMax];
  int MapState;
  int displayList;
  CCrystal Crystal;
  int RefreshFlag;

  float opacity;  /* this should be handled by the object transparency setting instead */
  
  float origin[3]; /* the origin of the plane */
  float system[9]; /* x, y, and z of the system */
  float grid;   /* sampling interval for the map */

  int min[2],max[2]; /* extents of the arrays */

  /* the data is normalized for easier ploting */
  int n_points;

  float *values;
  float *points;
  int   *flags;
  float *colors;
  float *normals;

  int   n_strips;
  int   *strips;

  /* orig_data = values*norm_factor+norm_constant */

  float norm_factor;
  float norm_constant;

  float ExtentMin[3];
  float ExtentMax[3];

  int ExtentFlag;

  CGO *UnitCellCGO;
} ObjectSliceState;

typedef struct ObjectSlice {
  CObject Obj;
  ObjectSliceState *State;
  int NState;
} ObjectSlice;

ObjectSlice *ObjectSliceFromMap(ObjectSlice * obj, ObjectMap* map,float grid,int state,int map_state);
/*void ObjectSliceDump(ObjectSlice *I,char *fname,int state);*/

PyObject *ObjectSliceAsPyList(ObjectSlice *I);
int ObjectSliceNewFromPyList(PyObject *list,ObjectSlice **result);

int ObjectSliceSetOpacity(ObjectSlice *I,float opacity,int state);
ObjectSliceState *ObjectSliceStateGetActive(ObjectSlice *I,int state);
void ObjectSliceStateValue2RGB(ObjectSliceState * s,float normalized_value,float * result);
int ObjectSliceGetOrigin(ObjectSlice *I,int state,float *origin);
void ObjectSliceDrag(ObjectSlice *I, int state, int mode, float *pt, float *mov, float *z_dir);
int ObjectSliceGetVertex(ObjectSlice *I,int index,int base, float *v);

#endif

