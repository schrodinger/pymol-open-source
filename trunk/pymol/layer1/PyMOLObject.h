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
#ifndef _H_PyMOLObject
#define _H_PyMOLObject

/* literally a 3-D object...also an object object */

#include"Base.h"
#include"Ray.h"
#include"Rep.h"
#include"Setting.h"

#define ObjNameMax        255
#define cObjectMolecule     1
#define cObjectMap          2
#define cObjectMesh         3
#define cObjectDist         4
#define cObjectCallback     5
#define cObjectCGO          6
#define cObjectSurface      7
#define cObjectGadget       8
#define cObjectCalculator   9

/* 
   the object base class is in the process of being converted to support
   states explicitly (an unfortunate early omission), which will allow
   for simplified implementation of future multi-state objects.
 */

typedef struct CObject {
  void (*fUpdate)(struct CObject *I); /* update representations */
  void (*fRender)(struct CObject *I,int frame,CRay *ray,Pickable **pick,int pass);
  void (*fFree)(struct CObject *I);
  int  (*fGetNFrame)(struct CObject *I);
  void (*fDescribeElement)(struct CObject *I,int index,char *buffer);
  void (*fInvalidate)(struct CObject *I,int rep,int level,int state);
  CSetting **(*fGetSettingHandle)(struct CObject *I,int state);
  int type;
  char Name[ObjNameMax];
  int Color;
  int RepVis[cRepCnt]; /* currently used only by non atomic objects */
  float ExtentMin[3],ExtentMax[3];
  int ExtentFlag,TTTFlag;
  float TTT[16]; /* translate, transform, translate matrix */
  CSetting *Setting;
  int Enabled; /* read-only... maintained by Scene */
  int Context; /* 0 = Camera, 1 = Unit Window, 2 = Scaled Window */
} CObject;

void ObjectInit(CObject *I);
void ObjectPurge(CObject *I);
void ObjectSetName(CObject *I,char *name);
void ObjectFree(CObject *I);
void ObjectUseColor(CObject *I);
void ObjectSetRepVis(CObject *I,int rep,int state);
void ObjectPrepareContext(CObject *I,CRay *ray);
void ObjectCombineTTT(CObject *I,float *ttt);
void ObjectSetTTTOrigin(CObject *I,float *origin);
void ObjectResetTTT(CObject *I);
PyObject *ObjectAsPyList(CObject *I);
int ObjectFromPyList(PyObject *list,CObject *I);
int ObjectGetCurrentState(CObject *I,int ignore_all_states);

#endif



