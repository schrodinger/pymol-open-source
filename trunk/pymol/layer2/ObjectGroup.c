/* 
A* -------------------------------------------------------------------
B* This file contains source code for the PyMOL computer program
C* copyright 1998-2006 by Warren Lyford Delano of DeLano Scientific. 
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

#include"os_predef.h"
#include"os_std.h"
#include"os_gl.h"

#include"OOMac.h"
#include"ObjectGroup.h"
#include"Base.h"
#include"MemoryDebug.h"
#include"PConv.h"

int ObjectGroupNewFromPyList(PyMOLGlobals *G,PyObject *list,ObjectGroup **result,int version)
{
#ifdef _PYMOL_NOPY
  return 0;
#else
  int ok = true;
  ObjectGroup *I=NULL;
  (*result) = NULL;
  if(ok) ok=(list!=Py_None);
  if(ok) ok=PyList_Check(list);

  I=ObjectGroupNew(G);
  if(ok) ok = (I!=NULL);
  if(ok) ok = ObjectFromPyList(G,PyList_GetItem(list,0),&I->Obj);
 
  return(ok);
#endif
}

PyObject *ObjectGroupAsPyList(ObjectGroup *I)
{
#ifdef _PYMOL_NOPY
  return NULL;
#else

  PyObject *result=NULL;

  result = PyList_New(1);
  PyList_SetItem(result,0,ObjectAsPyList(&I->Obj));

  return(PConvAutoNone(result));  
#endif
}


/*========================================================================*/

static void ObjectGroupFree(ObjectGroup *I) {
   ObjectPurge(&I->Obj);
  OOFreeP(I);
}


/*========================================================================*/
ObjectGroup *ObjectGroupNew(PyMOLGlobals *G)
{
  OOAlloc(G,ObjectGroup);

  ObjectInit(G,(CObject*)I);

  I->Obj.type = cObjectGroup;
  I->Obj.fFree = (void (*)(struct CObject *))ObjectGroupFree;
  I->Obj.fRender = NULL;
  I->OpenOrClosed = false;
  return(I);
}

