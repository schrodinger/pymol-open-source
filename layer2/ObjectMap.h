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

#include"Object.h"
#include"Crystal.h"
#include"Isosurf.h"

typedef struct ObjectMap {
  Object Obj;
  CCrystal *Crystal;
  int Div[3],Min[3],Max[3],FDim[4];
  Isofield *Field;
  int *N;
  float *V;
  int ResurfaceFlag;
} ObjectMap;

ObjectMap *ObjectMapNew(void);
ObjectMap *ObjectMapLoadXPLORFile(ObjectMap *obj,char *fname,int frame);
ObjectMap *ObjectMapReadXPLORStr(ObjectMap *I,char *XPLORStr,int frame);
int ObjectMapXPLORStrToMap(ObjectMap *I,char *XPLORStr,int frame);

#endif











