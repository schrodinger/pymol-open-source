
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

#ifdef _PYMOL_NOPY
#define WizardDoKey(...) false
#define WizardDoPick(...) false
#define WizardDoPosition(...) false
#define WizardDoScene(...) false
#define WizardDoSelect(...) false
#define WizardDoSpecial(...) false
#define WizardDoView(...) false
#define WizardFree(...) false
#define WizardGetBlock(...) nullptr
#define WizardInit(...) false
#define WizardRefresh(...)
#define WizardUpdate(...) false
#else

#include"Ortho.h"
#include "P.h"
#include "Result.h"
#include <vector>

int WizardInit(PyMOLGlobals * G);
int WizardActive(PyMOLGlobals * G);
pymol::Result<> WizardSet(PyMOLGlobals * G, PyObject * wiz, bool replace);
PyObject* WizardGet(PyMOLGlobals * G);
void WizardFree(PyMOLGlobals * G);
Block *WizardGetBlock(PyMOLGlobals * G);
void WizardRefresh(PyMOLGlobals * G);
int WizardDoPick(PyMOLGlobals * G, int bondFlag, int state=0);
int WizardDoSelect(PyMOLGlobals * G, const char* name, int state=0);
void WizardPurgeStack(PyMOLGlobals * G);
PyObject* WizardGetStack(PyMOLGlobals * G);
pymol::Result<> WizardSetStack(PyMOLGlobals * G, PyObject * wiz);
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

std::vector<unique_PyObject_ptr_auto_gil> WizardGetWizardCopies(PyMOLGlobals* G);
void WizardSetWizards(PyMOLGlobals* G, const std::vector<unique_PyObject_ptr_auto_gil>& wizs);

#endif
#endif
