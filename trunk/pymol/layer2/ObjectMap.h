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
#ifndef _H_ObjectMap
#define _H_ObjectMap

#include<Python.h>

#include"Object.h"
#include"Crystal.h"
#include"Isosurf.h"

typedef struct ObjectMap {
  Object Obj;
  CCrystal *Crystal;
  int Div[3],Min[3],Max[3],FDim[4];
  Isofield *Field;
  float Corner[8][3];
  int *Dim;
  float *Origin;
  float *Range;
  float *Grid;
} ObjectMap;

ObjectMap *ObjectMapNew(void);
ObjectMap *ObjectMapLoadXPLORFile(ObjectMap *obj,char *fname,int frame);
ObjectMap *ObjectMapReadXPLORStr(ObjectMap *I,char *XPLORStr,int frame);
int ObjectMapXPLORStrToMap(ObjectMap *I,char *XPLORStr,int frame);
ObjectMap *ObjectMapLoad(ObjectMap *obj,char *fname,int frame);
ObjectMap *ObjectMapLoadChemPyBrick(ObjectMap *I,PyObject *Map,
                                    int frame,int discrete);
ObjectMap *ObjectMapLoadCObject(ObjectMap *obj,int frame);
ObjectMap *ObjectMapLoadChemPyMap(ObjectMap *I,PyObject *Map,
                                  int frame,int discrete);

#endif











