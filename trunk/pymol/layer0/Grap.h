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

#ifndef _H_Grap
#define _H_Grap

/* Graphics wrapper macros and routines.
 * Not widely used at present, but that will change. */

void GrapMoveTo(int x,int y);
void GrapDrawStr(char *c,int x,int y);
void GrapDrawSubStrSafe(char *c,int x,int y,int start,int n);
void GrapDrawSubStrFast(char *c,int x,int y,int start,int n);
void GrapContStr(char *c);
int GrapMeasureStr(char *c);


#endif

