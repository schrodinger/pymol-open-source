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

void WizardInit(void);
int WizardActive(void);
void WizardSet(PyObject *wiz,int replace);
PyObject *WizardGet(void);
void WizardFree(void);
Block *WizardGetBlock(void);
void WizardRefresh(void);
int WizardDoPick(int bondFlag);
int WizardDoSelect(char *name);
void WizardPurgeStack(void);
PyObject *WizardGetStack(void);
int WizardSetStack(PyObject *wiz);
int WizardDoKey(unsigned char k, int x, int y, int mod);
int WizardDoSpecial(int k, int x, int y, int mod);
int WizardDoScene(void);
int WizardDoState(void);
int WizardDoFrame(void);

void WizardDirty(void);
int WizardUpdate(void);

#endif
