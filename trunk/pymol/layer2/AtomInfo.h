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


#define cAN_H   1
#define cAN_C   6
#define cAN_N   7
#define cAN_O   8
#define cAN_F   9
#define cAN_Na 11
#define cAN_Mg 12
#define cAN_P  15
#define cAN_S  16
#define cAN_Cl 17
#define cAN_K  19
#define cAN_Ca 20
#define cAN_Zn 30
#define cAN_Br 35
#define cAN_I  53

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
  int customType;
  int priority;
  float b,q,vdw,partialCharge;
  int formalCharge;
  signed char hetatm;
  short int model; 
  int atom;
  int selEntry;
  short int visRep[cRepCnt];
  int color;
  int id; 
  unsigned int flags;
  signed char bonded; /* be careful not to write at these as (int*) */
  signed char chemFlag;
  signed char geom;
  signed char valence;
  signed char deleteFlag;
  signed char masked;
  signed char protected;
  signed char protons;
} AtomInfoType;

int *AtomInfoGetSortedIndex(AtomInfoType *rec,int n,int **outdex);
void AtomInfoAssignParameters(AtomInfoType *I);
void AtomInfoFreeSortedIndexes(int *index,int *outdex);
void AtomInfoPrimeColors(void);
int AtomInfoGetColor(AtomInfoType *at1);
int AtomInfoGetExpectedValence(AtomInfoType *I);

int AtomInfoMatch(AtomInfoType *at1,AtomInfoType *at2);
int AtomInfoAltMatch(AtomInfoType *at1,AtomInfoType *at2);
int AtomInfoCompare(AtomInfoType *at1,AtomInfoType *at2);
float AtomInfoGetBondLength(AtomInfoType *ai1,AtomInfoType *ai2);
int AtomInfoSameResidue(AtomInfoType *at1,AtomInfoType *at2);

void AtomInfoBracketResidue(AtomInfoType *ai0,int n0,AtomInfoType *ai,int *st,int *nd);
void AtomInfoBracketResidueFast(AtomInfoType *ai0,int n0,int cur,int *st,int *nd);

void AtomInfoUniquefyNames(AtomInfoType *atInfo0,int n0,AtomInfoType *atInfo1,int n1);
void AtomInfoCombine(AtomInfoType *dst,AtomInfoType *src);
int AtomInfoNameOrder(AtomInfoType *at1,AtomInfoType *at2);

#endif
