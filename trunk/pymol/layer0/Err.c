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
#include"os_std.h"

#include"Err.h"
#include"Ortho.h"

void ErrFatal(const char *where,const char *what)
{
  fprintf(stderr,"%s-ERR: %s\n",where,what);
  fflush(stderr);
  exit(1);
}

int ErrMessage(const char *where,const char *what)
{
  char buffer[1024];
  sprintf(buffer,"%s-ERR: %s\n",where,what);
  OrthoAddOutput(buffer);
  OrthoRestorePrompt();
  return(0);
}

int ErrOk(const char *where,const char *what)
{
  char buffer[1024];
  sprintf(buffer,"%s: %s\n",where,what);
  OrthoAddOutput(buffer);
  OrthoRestorePrompt();
  return(1);
}

void ErrPointer(const char *file,int line)
{
  fprintf(stderr,"NULL-POINTER-ERR: in %s line %i\n",file,line);
  fflush(stderr);
    while(1);
  exit(1);
}





