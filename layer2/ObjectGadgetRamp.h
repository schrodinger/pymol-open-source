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

  int NLevel;
  float *Level;
  float *Color;
  int var_index;

  /* cRampMap */

  char SrcName[ObjNameMax];
  int SrcState;

  int CalcMode;

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


#define cRAMP_TRADITIONAL 1
#define cRAMP_SLUDGE 2
#define cRAMP_OCEAN 3
#define cRAMP_HOT 4
#define cRAMP_GRAYABLE 5
#define cRAMP_RAINBOW 6
#define cRAMP_AFMHOT 7
#define cRAMP_GRAYSCALE 8

ObjectGadgetRamp *ObjectGadgetRampNew(PyMOLGlobals *G);

ObjectGadgetRamp *ObjectGadgetRampMapNewAsDefined(PyMOLGlobals *G,ObjectMap *map,PyObject *level,
                                                  PyObject *color,int map_state,float *vert_vla,
                                                  float beyond,float within,float sigma,int zero);

int ObjectGadgetRampInterpolate(ObjectGadgetRamp *I,float level,float *color);
int ObjectGadgetRampInterVertex(ObjectGadgetRamp *I,float *pos,float *color);

PyObject *ObjectGadgetRampAsPyList(ObjectGadgetRamp *I);
int ObjectGadgetRampNewFromPyList(PyMOLGlobals *G,PyObject *list,ObjectGadgetRamp **result,int version);

void ObjectGadgetRampUpdate(ObjectGadgetRamp *I);
void ObjectGadgetRampFree(ObjectGadgetRamp *I);

#endif












