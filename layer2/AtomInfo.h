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

#define cResnLen 5
#define cResiLen 5
#define cAtomNameLen 4
#define cSegiLen 4
#define cTextTypeLen 20
#define cLabelTypeLen 20

#define cAtomInfoTetrahedral 4
#define cAtomInfoPlaner 3
#define cAtomInfoLinear 2
#define cAtomInfoSingle 1
#define cAtomInfoNone 5

typedef char Chain[2];

typedef char SegIdent[cSegiLen+1];
typedef char ResIdent[cResiLen+1];
typedef char ResName[cResnLen+1];
typedef char AtomName[cAtomNameLen+1];
typedef char TextType[cTextTypeLen+1];
typedef char LabelType[cLabelTypeLen+1];

#define cAtomInfoNoType -9999

typedef struct AtomInfoType {
  int resv;
  Chain chain;
  Chain alt;
  ResIdent resi;
  SegIdent segi;
  ResName resn;
  AtomName name;
  AtomName elem;
  TextType textType;
  LabelType label;
  int hydrogen;
  int ludiType; /* this will go away soon */
  int customType;
  int priority;
  float b,q,vdw,partialCharge;
  int formalCharge;
  signed char hetatm;
  short int model; /* remaining items only used during selection */
  int atom;
  int selEntry;
  short int visRep[cRepCnt];
  int color;
  int id; /* used for reading conect records */
  unsigned int flags;
  int bonded;
  int chemFlag;
  int geom;
  int valence;
} AtomInfoType;

int *AtomInfoGetSortedIndex(AtomInfoType *rec,int n,int **outdex);
void AtomInfoAssignParameters(AtomInfoType *I);
void AtomInfoFreeSortedIndexes(int *index,int *outdex);
void AtomInfoPrimeColors(void);
int AtomInfoGetColor(AtomInfoType *at1);

int AtomInfoMatch(AtomInfoType *at1,AtomInfoType *at2);
int AtomInfoAltMatch(AtomInfoType *at1,AtomInfoType *at2);
int AtomInfoCompare(AtomInfoType *at1,AtomInfoType *at2);

#endif
