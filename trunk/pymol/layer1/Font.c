
/* 
A* -------------------------------------------------------------------
B* This file contains source code for the PyMOL computer program
C* copyright 1998-2003 by Warren Lyford Delano of DeLano Scientific. 
D* -------------------------------------------------------------------
E* It is unlawful to modify or remove this copyright notice.
F* -------------------------------------------------------------------
G* Please see the accompanying LICENSE file for further information. 
H* --------------------------------------------------\-----------------
I* Additional authors of this source file include:
-* 
-* 
-*
Z* -------------------------------------------------------------------
*/

#include"Util.h"
#include"Font.h"

int FontInit(PyMOLGlobals * G, CFont * I)
{
  UtilZeroMem(I, sizeof(CFont));
  I->G = G;
  return 1;
}

void FontPurge(CFont * I)
{

}
