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

typedef struct {
  int MaxAlloc;
  int LastFree;
  int NewestUsed;
  int OldestUsed;
  int NUsed;
  int TargetMaxUsage; /* don't store more than this many pixmaps in RAM */
  int *Hash;
  int RetainAll;
  CharRec *Char;
} CCharacter;

void CharacterInit(void);
void CharacterFree(void);

int CharacterGetNew(void);
int CharacterGetWidth(int id);
int CharacterGetHeight(int id);

int CharacterNewFromBitmap(int width, int height,
                           unsigned char *bitmap,
                           CharFngrprnt *fprnt);

int CharacterFind(CharFngrprnt *fprnt);

float CharacterInterpolate(int id,float *v);
void CharacterSetRetention(int retail_all);

#endif

