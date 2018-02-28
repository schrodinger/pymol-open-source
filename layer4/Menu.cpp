

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
#include"os_python.h"

#include"os_predef.h"

#include"os_std.h"

#include "Menu.h"
#include "PopUp.h"
#include "P.h"
#include "Ortho.h"

void MenuActivate(PyMOLGlobals * G, int x, int y, int last_x, int last_y, int passive,
                  const char *name, const char *sele)
{
#ifndef _PYMOL_NOPY
  PyObject *list;

  PBlock(G);

  list = PYOBJECT_CALLMETHOD(P_menu, name, "Os", G->P_inst->cmd, sele);
  if(PyErr_Occurred())
    PyErr_Print();
  if(list) {
    PopUpNew(G, x, y, last_x, last_y, passive, list, NULL);
    Py_DECREF(list);
  }
  PUnblock(G);
#endif

}

void MenuActivate3fv(PyMOLGlobals * G, int x, int y, int last_x, int last_y, int passive,
                     const char *name, const float *xyz)
{
#ifndef _PYMOL_NOPY
  PyObject *list;

  PBlock(G);

  list =
    PYOBJECT_CALLMETHOD(P_menu, name, "O(fff)(ii)", G->P_inst->cmd,
        xyz[0], xyz[1], xyz[2], x, y);
  if(PyErr_Occurred())
    PyErr_Print();
  if(list) {
    PopUpNew(G, x, y, last_x, last_y, passive, list, NULL);
    Py_DECREF(list);
  }
  PUnblock(G);
#endif
}

void MenuActivate2Arg(PyMOLGlobals * G, int x, int y, int last_x, int last_y, int passive,
                      const char *name, const char *sele1, const char *sele2)
{
#ifndef _PYMOL_NOPY
  PyObject *list;

  PBlock(G);

  list = PYOBJECT_CALLMETHOD(P_menu, name, "Oss", G->P_inst->cmd, sele1, sele2);
  if(PyErr_Occurred())
    PyErr_Print();
  if(list) {
    PopUpNew(G, x, y, last_x, last_y, passive, list, NULL);
    Py_DECREF(list);
  }
  PUnblock(G);
#endif
}

Block *MenuActivate1Arg(PyMOLGlobals * G, int x, int y, int last_x, int last_y, int passive,
			const char *name, const char *arg1)
{
  Block *block = NULL;
#ifndef _PYMOL_NOPY
  PyObject *list;
  PBlock(G);

  list = PYOBJECT_CALLMETHOD(P_menu, name, "Os", G->P_inst->cmd, arg1);
  if(PyErr_Occurred())
    PyErr_Print();
  if(list) {
    block = PopUpNew(G, x, y, last_x, last_y, passive, list, NULL);
    Py_DECREF(list);
  }
  PUnblock(G);
#endif
  return block;
}

void MenuActivate0Arg(PyMOLGlobals * G, int x, int y, int last_x, int last_y, int passive,
                      const char *name)
{
#ifndef _PYMOL_NOPY
  PyObject *list;

  PBlock(G);

  list = PYOBJECT_CALLMETHOD(P_menu, (char*)name, "O", G->P_inst->cmd);
  if(PyErr_Occurred())
    PyErr_Print();
  if(list) {
    PopUpNew(G, x, y, last_x, last_y, passive, list, NULL);
    Py_DECREF(list);
  }
  PUnblock(G);
#endif
}
