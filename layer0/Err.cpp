

/* 
A* -------------------------------------------------------------------
B* This file contains source code for the PyMOL computer program
C* Copyright (c) Schrodinger, LLC. 
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
#include"os_std.h"

#include "PyMOLGlobals.h"
#include"Err.h"
#include"Ortho.h"
#include"Feedback.h"

void ErrFatal(const PyMOLGlobals *, const char *where, const char *what)
{
  fprintf(stderr, "%s-Error: %s\n", where, what);
  fflush(stderr);
  exit(1);
}

int ErrMessage(PyMOLGlobals * G, const char *where, const char *what)
{
    /* unclassified errors are assigned to the Executive catch-all */

  PRINTFB(G, FB_Executive, FB_Errors)
    "%s-Error: %s\n", where, what
    ENDFB(G);
  return (0);
}

void ErrPointer(const PyMOLGlobals *, const char *file, int line)
{
  fprintf(stderr, "NULL-POINTER-ERROR: in %s line %i\n", file, line);
  printf
    ("****************************************************************************\n");
  printf
    ("*** EEK!  PyMOL just ran out of memory and crashed.  To get around this, ***\n");
  printf
    ("*** you may need to reduce the quality, size, or complexity of the scene ***\n");
  printf
    ("*** that you are viewing or rendering.    Sorry for the inconvenience... ***\n");
  printf
    ("****************************************************************************\n");
  exit(EXIT_FAILURE);
}
