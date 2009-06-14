
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
  PyMOLGlobals *G;
  /* stored in a session */

  int Active;
  ObjectNameType MapName;
  int MapState;
  float MapMean;
  float MapStdev;

  float ExtentMin[3];
  float ExtentMax[3];
  int ExtentFlag;

  float origin[3];              /* the origin of the plane */
  float system[9];              /* x, y, and z of the system */

  /* not stored in session */

  int RefreshFlag;
  int min[2], max[2];           /* extents of the arrays */
  int displayList;
  int displayListInvalid;
  float last_scale;

  /* the data is normalized for easier ploting */
  int n_points;

  float *values;
  float *points;
  int *flags;
  float *colors;
  float *normals;

  int n_strips;
  int *strips;

} ObjectSliceState;

typedef struct ObjectSlice {
  CObject Obj;
  ObjectSliceState *State;
  int NState;
} ObjectSlice;

ObjectSlice *ObjectSliceFromMap(PyMOLGlobals * G, ObjectSlice * obj, ObjectMap * map,
                                int state, int map_state);

/*void ObjectSliceDump(ObjectSlice *I,char *fname,int state);*/

PyObject *ObjectSliceAsPyList(ObjectSlice * I);
int ObjectSliceNewFromPyList(PyMOLGlobals * G, PyObject * list, ObjectSlice ** result);

ObjectSliceState *ObjectSliceStateGetActive(ObjectSlice * I, int state);
void ObjectSliceStateValue2RGB(ObjectSliceState * s, float normalized_value,
                               float *result);
int ObjectSliceGetOrigin(ObjectSlice * I, int state, float *origin);
void ObjectSliceDrag(ObjectSlice * I, int state, int mode, float *pt, float *mov,
                     float *z_dir);
int ObjectSliceGetVertex(ObjectSlice * I, int index, int base, float *v);

#endif
