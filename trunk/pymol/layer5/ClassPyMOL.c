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

#include "os_std.h"
#include "MemoryDebug.h"

#include "ClassPyMOL.h"
#include "PyMOLGlobals.h"

PyMOLGlobals *TempPyMOLGlobals;

typedef struct _ClassPyMOL {
  PyMOLGlobals *G;
} _ClassPyMOL;

ClassPyMOL *ClassPyMOLNew(void)
{
  ClassPyMOL *result = NULL;

  /* allocate global container */

  if( (result = Calloc(ClassPyMOL,1)) ) {
    if( (result->G = Calloc(PyMOLGlobals,1)) ) {

      /* temporary global pointer for the transition period */

      TempPyMOLGlobals=result->G;

      /* continue initialization */

    } else {
      FreeP(result);
    }
  }
  return result;
}

void ClassPyMOLFree(ClassPyMOL *I)
{
  /* take PyMOL down gracefully */

  FreeP(I->G);
  FreeP(I);

}

struct _PyMOLGlobals *ClassPyMOLGetGlobals(ClassPyMOL *I)
{
  return I->G;
}
