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

#ifndef _H_Match
#define _H_Match

typedef struct {
  float **mat;
  int *pair;
  int na,nb;
} CMatch;

CMatch *MatchNew(unsigned int na,unsigned int nb);
void MatchFree(CMatch *I);
float MatchAlign(CMatch *I,float gap_penalty,float ext_penalty,int max_gap);

#endif
