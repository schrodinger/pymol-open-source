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
#ifndef _H_ObjectGadgetRamp
#define _H_ObjectGadgetRamp

#include"os_python.h"

#include"PyMOLObject.h"
#include"CGO.h"

#include"GadgetSet.h"
#include"ObjectMap.h"

#define dRampNone 0
#define cRampMap 1
#define cRampProp 2

typedef struct ObjectGadgetRamp {
  ObjectGadget Gadget;

  int RampType; 

  int NColor;
  float *Level;
  float *Color;
  int var_index;

  /* cRampMap */

  char SrcName[ObjNameMax];
  int SrcState;

  /* fields below are not saved in session */
  ObjectMap *Map;

  float border;
  float width;
  float height;
  float bar_width;
  float bar_height;
  float text_height;
  float text_raise;
  float text_border;
  float text_scale_h;
  float text_scale_v;
  float x,y;

} ObjectGadgetRamp;

ObjectGadgetRamp *ObjectGadgetRampNew(void);

ObjectGadgetRamp *ObjectGadgetRampMapNewAsDefined(ObjectMap *map,PyObject *level,
                                               PyObject *color,int map_state);

int ObjectGadgetRampInterpolate(ObjectGadgetRamp *I,float level,float *color);
int ObjectGadgetRampInterVertex(ObjectGadgetRamp *I,float *pos,float *color);

PyObject *ObjectGadgetRampAsPyList(ObjectGadgetRamp *I);
int ObjectGadgetRampNewFromPyList(PyObject *list,ObjectGadgetRamp **result);

void ObjectGadgetRampUpdate(ObjectGadgetRamp *I);
void ObjectGadgetRampFree(ObjectGadgetRamp *I);

#endif












