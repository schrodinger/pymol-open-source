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
#ifndef _H_Character
#define _H_Character

#include "Pixmap.h"
#include "PyMOLGlobals.h"

typedef unsigned char CharColor[4];


typedef struct {
  int text_id;
  short int ch;
  short int height;
  CharColor color;
} CharInfo; /* currently 5 short ints */

typedef struct {
  unsigned short int data[6];
} CharData;

typedef union {
  CharInfo i;
  CharData d;
} CharUnion;

typedef struct {
  unsigned short hash_code;
  CharUnion u;
} CharFngrprnt;

typedef struct {
  int Active;
  CPixmap Pixmap;
  int Width;
  int Height;
  int Next,Prev,HashNext,HashPrev;
  CharFngrprnt Fngrprnt;
} CharRec;

struct _CCharacter {
  int MaxAlloc;
  int LastFree;
  int NewestUsed;
  int OldestUsed;
  int NUsed;
  int TargetMaxUsage; /* don't store more than this many pixmaps in RAM */
  int *Hash;
  int RetainAll;
  CharRec *Char;
};

int CharacterInit(PyMOLGlobals *G);
void CharacterFree(PyMOLGlobals *G);

int CharacterGetNew(PyMOLGlobals *G);
int CharacterGetWidth(PyMOLGlobals *G,int id);
int CharacterGetHeight(PyMOLGlobals *G,int id);

int CharacterNewFromBitmap(PyMOLGlobals *G,int width, int height,
                           unsigned char *bitmap,
                           CharFngrprnt *fprnt);

int CharacterFind(PyMOLGlobals *G,CharFngrprnt *fprnt);

float CharacterInterpolate(PyMOLGlobals *G,int id,float *v);
void CharacterSetRetention(PyMOLGlobals *G,int retail_all);

#endif

