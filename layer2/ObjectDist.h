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
#ifndef _H_ObjectDist
#define _H_ObjectDist

#include"Object.h"
#include"AtomInfo.h"
#include"Vector.h"
#include"Color.h"

typedef struct ObjectDist {
  Object Obj;
  struct DistSet **DSet;
  int NDSet;
  AtomInfoType *AtomInfo;
  int NAtom;
  int CurDSet;
  char Color[3];
} ObjectDist;

ObjectDist *ObjectDistNew(int sele1,int sele2,int mode,float cutoff);
void ObjectDistInvalidateRep(ObjectDist *I,int rep);

#endif











