
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

#include<string.h>
#include<Python.h>

#include "Menu.h"
#include "PopUp.h"
#include "P.h"
#include "Ortho.h"


void MenuActivate(int x,int y,char *name,char *sele)
{

  PyObject *list;

  PBlock(); /* menu doesn't currently call API, so leave it locked */

  list = PyObject_CallMethod(P_menu,name,"s",sele); 
  if(PyErr_Occurred()) PyErr_Print();
  if(list) {
    PopUpNew(x,y,list);
    Py_DECREF(list);
  }
  PUnblock();
}

