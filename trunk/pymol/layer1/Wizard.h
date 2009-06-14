
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
#ifndef _H_Wizard
#define _H_Wizard

#include"Ortho.h"

int WizardInit(PyMOLGlobals * G);
int WizardActive(PyMOLGlobals * G);
void WizardSet(PyMOLGlobals * G, PyObject * wiz, int replace);
PyObject *WizardGet(PyMOLGlobals * G);
void WizardFree(PyMOLGlobals * G);
Block *WizardGetBlock(PyMOLGlobals * G);
void WizardRefresh(PyMOLGlobals * G);
int WizardDoPick(PyMOLGlobals * G, int bondFlag);
int WizardDoSelect(PyMOLGlobals * G, char *name);
void WizardPurgeStack(PyMOLGlobals * G);
PyObject *WizardGetStack(PyMOLGlobals * G);
int WizardSetStack(PyMOLGlobals * G, PyObject * wiz);
int WizardDoKey(PyMOLGlobals * G, unsigned char k, int x, int y, int mod);
int WizardDoSpecial(PyMOLGlobals * G, int k, int x, int y, int mod);
int WizardDoScene(PyMOLGlobals * G);
int WizardDoState(PyMOLGlobals * G);
int WizardDoFrame(PyMOLGlobals * G);
int WizardDoDirty(PyMOLGlobals * G);
int WizardDoView(PyMOLGlobals * G, int force);
int WizardDoPosition(PyMOLGlobals * G, int force);

void WizardDirty(PyMOLGlobals * G);
int WizardUpdate(PyMOLGlobals * G);

#endif
