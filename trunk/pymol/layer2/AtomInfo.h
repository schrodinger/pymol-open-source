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
#define cResiLen 4
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
  float b,vdw;
  int model; /* remaining items only used during selection */
  int atom; 
  int selEntry;
  int visRep[cRepCnt];
  int color;
} AtomInfoType;

#endif
