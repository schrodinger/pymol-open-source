/* 
A* -------------------------------------------------------------------
B* This file contains source code for the PyMOL computer program
C* copyright 1998-2000 by Warrn Lyford Delano of DeLano Scientific. 
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

#ifndef _H_ClassPyMOL
#define _H_ClassPyMOL

#ifndef COption_DEFINED
typedef struct _COption COption;
#define COption_DEFINED
#endif

/* creation and destruction */

typedef struct _ClassPyMOL ClassPyMOL;

ClassPyMOL *ClassPyMOLNew(void);
ClassPyMOL *ClassPyMOLNewWithOptions(COption *option);

void ClassPyMOLFree(ClassPyMOL *I);

COption *ClassPyMOLOptionsNew(void);
void ClassPyMOLOptionsFree(COption *option);

/* developer/transient privates */

struct _PyMOLGlobals *ClassPyMOLGetGlobals(ClassPyMOL *I);

#endif
