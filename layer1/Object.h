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
#ifndef _H_Object
#define _H_Object

/* literally a 3-D object...also an object object */

#include"Base.h"
#include"Ray.h"
#include"Rep.h"

#define ObjNameMax 255
#define cObjectMolecule 1
#define cObjectMap 2
#define cObjectMesh 3
#define cObjectDist 4

typedef struct Object {
  void (*fUpdate)(struct Object *I); /* update representations */
  void (*fRender)(struct Object *I,int frame,CRay *ray,Pickable **pick);
  void (*fFree)(struct Object *I);
  int  (*fGetNFrame)(struct Object *I);
  void (*fDescribeElement)(struct Object *I,int index);
  int type;
  char Name[ObjNameMax];
  int Color;
  int RepVis[cRepCnt]; /* currently used only by non atomic objects */
} Object;


void ObjectInit(Object *I);
void ObjectSetName(Object *I,char *name);
void ObjectFree(Object *I);
void ObjectUseColor(Object *I);
void ObjectSetRepVis(Object *I,int rep,int state);

#endif



