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
#ifndef _H_AtomInfo
#define _H_AtomInfo


#include"Rep.h"

#define cResnLen 20
#define cResiLen 5
#define cAtomNameLen 4
#define cSegiLen 4


typedef char Chain[2];

typedef char SegIdent[cSegiLen+1];
typedef char ResIdent[cResiLen+1];
typedef char ResName[cResnLen+1];
typedef char AtomName[cAtomNameLen+1];

typedef struct AtomInfoType {
  int resv;
  Chain chain;
  ResIdent resi;
  SegIdent segi;
  ResName resn;
  AtomName name;
  int ludiType;
  int customType;
  int customFlag;
  int priority;
  float b,q,vdw;
  signed char hetatm;
  short int model; /* remaining items only used during selection */
  int atom; 
  int selEntry;
  short int visRep[cRepCnt];
  int color;
  int tmpID; /* used for reading conect records */
} AtomInfoType;

int *AtomInfoGetSortedIndex(AtomInfoType *rec,int n,int **outdex);
void AtomInfoAssignParameters(AtomInfoType *I);
void AtomInfoFreeSortedIndexes(int *index,int *outdex);
int AtomInfoMatch(AtomInfoType *at1,AtomInfoType *at2);
int AtomInfoCompare(AtomInfoType *at1,AtomInfoType *at2);

#endif
