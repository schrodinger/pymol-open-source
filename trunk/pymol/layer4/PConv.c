
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

#include<stdlib.h>
#include<Python.h>
#include<signal.h>
#include<string.h>
#include<sys/types.h>
#include<sys/time.h>
#include<unistd.h>

#include"MemoryDebug.h"
#include"Base.h"
#include"PConv.h"
#include"PUtils.h"

void PConv44To44f(PyObject *src,float *dest)
{
  PBlock();

  PUnblock();
}
