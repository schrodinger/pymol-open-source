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

#include"os_std.h"
#include"os_memory.h"
#include"os_python.h"

#include"champ.h"

static PyObject *AutoNone(PyObject *result) /* automatically own Py_None */
{
  if(result==Py_None)
    Py_INCREF(result);
  else if(result==NULL) {
    result=Py_None;
    Py_INCREF(result);
  } 
  return(result);
}

static PyObject *_new(PyObject *self,      PyObject *args)
{
  return(PyCObject_FromVoidPtr((void*)ChampNew(),(void*)(void*)ChampFree));
}

static PyObject *_memory_dump(PyObject *self,      PyObject *args)
{
  os_memory_dump();
  return(AutoNone(Py_None));
}

static PyMethodDef champ_methods[] = {
  {"new",	              _new,                   METH_VARARGS },
  {"memory_dump",         _memory_dump,           METH_VARARGS },
  {NULL,		              NULL}     /* sentinel */        
};

void init_champ(void)
{
  Py_InitModule("_champ", champ_methods);
}

