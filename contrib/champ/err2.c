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

#include"err2.h"

void err_fatal(const char *where,const char *what)
{
  fprintf(stderr,"%s-ERR: %s\n",where,what);
  fflush(stderr);
  exit(1);
}

int err_message(const char *where,const char *what)
{
  printf("%s-ERR: %s\n",where,what);
  fflush(stderr);
  return(0);
}

