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

#include "PyMOLGlobals.h"

typedef struct {
  PyMOLGlobals *G;
  float **smat;
  float **mat;
  float **da, **db;
  int na,nb;
  int *pair;
  float score; /* result */
  int n_pair;
} CMatch;

CMatch *MatchNew(PyMOLGlobals *G,unsigned int na,unsigned int nb,int dist_mats);
int MatchResidueToCode(CMatch *I,int *vla,int n);
int MatchMatrixFromFile(CMatch *I,char *fname,int quiet);
int MatchPreScore(CMatch *I,int *vla1,int n1,int *vla2,int n2,int quiet);
void MatchFree(CMatch *I);
int MatchAlign(CMatch *I,float gap_penalty,float ext_penalty,
               int max_gap,int max_skip,int quiet,
               int window,float ante);

#endif
