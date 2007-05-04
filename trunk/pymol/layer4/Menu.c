
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

#include"os_predef.h"
#include"os_python.h"

#include"os_std.h"

#include "Menu.h"
#include "PopUp.h"
#include "P.h"
#include "Ortho.h"


void MenuActivate(PyMOLGlobals *G,int x,int y,int last_x,int last_y,int passive,
                  char *name,char *sele)
{
#ifndef _PYMOL_NOPY
  PyObject *list;

  PBlock(G); 
 
  list = PyObject_CallMethod(P_menu,name,"Os",G->P_inst->cmd,sele); 
  if(PyErr_Occurred()) PyErr_Print();
  if(list) {
    PopUpNew(G,x,y,last_x,last_y,passive,list,NULL);
    Py_DECREF(list);
  }
  PUnblock(G);
#endif

}

void MenuActivate2Arg(PyMOLGlobals *G,int x,int y,int last_x,int last_y,int passive,
                      char *name,char *sele1,char *sele2)
{
#ifndef _PYMOL_NOPY
  PyObject *list;

  PBlock(G); 

  list = PyObject_CallMethod(P_menu,name,"Oss",G->P_inst->cmd,sele1,sele2); 
  if(PyErr_Occurred()) PyErr_Print();
  if(list) {
    PopUpNew(G,x,y,last_x,last_y,passive,list,NULL);
    Py_DECREF(list);
  }
  PUnblock(G);
#endif
}

void MenuActivate0Arg(PyMOLGlobals *G,int x,int y,int last_x,int last_y,int passive,
                      char *name)
{
#ifndef _PYMOL_NOPY
  PyObject *list;

  PBlock(G); 

  list = PyObject_CallMethod(P_menu,name,"O",G->P_inst->cmd);
  if(PyErr_Occurred()) PyErr_Print();
  if(list) {
    PopUpNew(G,x,y,last_x,last_y,passive,list,NULL);
    Py_DECREF(list);
  }
  PUnblock(G);
#endif
}

