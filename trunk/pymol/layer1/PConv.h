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
#ifndef _H_PConv
#define _H_PConv

#include<Python.h>


/* Convenient conversion routines for C<->Python data interchange
   
   Note that all of these routines assume that we have the global
   interpreter lock - blocking all other threads.
   
   There are three ways to get it:
   
   - call PBlock() [followe by PUnblock() when done]
   
   - call PBlockAndUnlockAPI - [followed by PLockAPIAndUnblock() when
   done]
   
   - or in response to a call to the PM API, you will have the main
   Python thread by default.  [Note that within an
   APIEntry(),APIExit() block the lock is released, so these
   functions should be called outside of that block].

*/

PyObject *PConvFloatVLAToPyList(float *f);

void PConvFloatToPyDictItem(PyObject *dict,char *key,float f);
void PConvStringToPyDictItem(PyObject *dict,char *key,char *f);
void PConvIntToPyDictItem(PyObject *dict,char *key,int i);

void PConvFloat3ToPyObjAttr(PyObject *obj,char *attr,float *v);
void PConvFloatToPyObjAttr(PyObject *obj,char *attr,float f);
void PConvIntToPyObjAttr(PyObject *obj,char *attr,int i);
void PConvInt2ToPyObjAttr(PyObject *obj,char *attr,int *v);
void PConvStringToPyObjAttr(PyObject *obj,char *attr,char *f);

int PConvPyObjectToFloat(PyObject *object,float *value);
int PConvPyObjectToInt(PyObject *object,int *value);
int PConvPyObjectToStrMaxLen(PyObject *object,char *value,int ln);
int PConvPyObjectToStrMaxClean(PyObject *object,char *value,int ln);

PyObject *PConvStringListToPyList(int l,char **str);

void PConv44PyListTo44f(PyObject *src,float *dest); /* note loss of precision */


#endif







