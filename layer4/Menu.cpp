

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

/**
 * Shows a menu list. A GUI that has its own menus sets pymol.menu.popup
 * (the Qt window does) and gets the list; otherwise, or when that returns
 * false, the menu is drawn in the viewport.
 * @return the viewport popup, nullptr when the GUI took the menu
 */
static Block* MenuShow(PyMOLGlobals* G, int x, int y, int last_x, int last_y,
    int passive, PyObject* list)
{
  bool taken = false;
  PyObject* popup = PyObject_GetAttrString(P_menu, "popup");
  if (popup && popup != Py_None) {
    PyObject* result = PyObject_CallFunctionObjArgs(popup, list, nullptr);
    taken = result && PyObject_IsTrue(result);
    Py_XDECREF(result);
  }
  Py_XDECREF(popup);
  if (PyErr_Occurred())
    PyErr_Print();
  return taken ? nullptr
               : PopUpNew(G, x, y, last_x, last_y, passive, list, nullptr);
}

void MenuActivate(PyMOLGlobals * G, int x, int y, int last_x, int last_y, int passive,
                  const char *name, const char *sele)
{
  PyObject *list;

  PBlock(G);

  list = PYOBJECT_CALLMETHOD(P_menu, name, "Os", G->P_inst->cmd, sele);
  if(PyErr_Occurred())
    PyErr_Print();
  if(list) {
    MenuShow(G, x, y, last_x, last_y, passive, list);
    Py_DECREF(list);
  }
  PUnblock(G);

}

void MenuActivate3fv(PyMOLGlobals * G, int x, int y, int last_x, int last_y, int passive,
                     const char *name, const float *xyz)
{
  PyObject *list;

  PBlock(G);

  list =
    PYOBJECT_CALLMETHOD(P_menu, name, "O(fff)(ii)", G->P_inst->cmd,
        xyz[0], xyz[1], xyz[2], x, y);
  if(PyErr_Occurred())
    PyErr_Print();
  if(list) {
    MenuShow(G, x, y, last_x, last_y, passive, list);
    Py_DECREF(list);
  }
  PUnblock(G);
}

void MenuActivate2Arg(PyMOLGlobals * G, int x, int y, int last_x, int last_y, int passive,
                      const char *name, const char *sele1, const char *sele2)
{
  PyObject *list;

  PBlock(G);

  list = PYOBJECT_CALLMETHOD(P_menu, name, "Oss", G->P_inst->cmd, sele1, sele2);
  if(PyErr_Occurred())
    PyErr_Print();
  if(list) {
    MenuShow(G, x, y, last_x, last_y, passive, list);
    Py_DECREF(list);
  }
  PUnblock(G);
}

Block *MenuActivate1Arg(PyMOLGlobals * G, int x, int y, int last_x, int last_y, int passive,
			const char *name, const char *arg1)
{
  Block *block = nullptr;
  PyObject *list;
  PBlock(G);

  list = PYOBJECT_CALLMETHOD(P_menu, name, "Os", G->P_inst->cmd, arg1);
  if(PyErr_Occurred())
    PyErr_Print();
  if(list) {
    block = MenuShow(G, x, y, last_x, last_y, passive, list);
    Py_DECREF(list);
  }
  PUnblock(G);
  return block;
}

void MenuActivate0Arg(PyMOLGlobals * G, int x, int y, int last_x, int last_y, int passive,
                      const char *name)
{
  PyObject *list;

  PBlock(G);

  list = PYOBJECT_CALLMETHOD(P_menu, (char*)name, "O", G->P_inst->cmd);
  if(PyErr_Occurred())
    PyErr_Print();
  if(list) {
    MenuShow(G, x, y, last_x, last_y, passive, list);
    Py_DECREF(list);
  }
  PUnblock(G);
}
