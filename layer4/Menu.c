
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
#include "PM.h"
#include "PopUp.h"
#include "PUtils.h"
#include "Ortho.h"

extern PyThreadState *_save;

void MenuActivate(int x,int y,char *name,char *sele)
{
  /* assumes a locked API and unblocked python threads (GLUT event) */

  OrthoLineType buffer;
  PyObject *list;

  PBlock(&_save);

  sprintf(buffer,"pymol_menu = pmm.%s('%s')",name,sele);
  PyRun_SimpleString(buffer);
  list = PyDict_GetItemString(PM_Globals,"pymol_menu");
  PopUpNew(x,y,list);

  PUnblock(&_save);
}

