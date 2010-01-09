

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

#include"Err.h"
#include"Ortho.h"
#include"Feedback.h"

void ErrFatal(PyMOLGlobals * G, const char *where, const char *what)
{
  fprintf(stderr, "%s-Error: %s\n", where, what);
  fflush(stderr);
  exit(1);
}

int ErrMessage(PyMOLGlobals * G, const char *where, const char *what)
{
  char buffer[1024];
  if(Feedback(G, FB_Executive, FB_Errors)) {

    /* unclassified errors are assigned to the Executive catch-all */

    sprintf(buffer, "%s-Error: %s\n", where, what);
    OrthoAddOutput(G, buffer);
    OrthoRestorePrompt(G);
  }
  return (0);
}

void ErrPointer(PyMOLGlobals * G, const char *file, int line)
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
