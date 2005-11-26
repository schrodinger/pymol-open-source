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
  int text_id; /* 32 bits */
  short int ch; /* 16 bits */
  short int size; /* 16 bits */
  CharColor color; /* 32 bits */
  short int flat; /* 16 bits */
} CharInfo; /* 7 short ints */

typedef struct {
  unsigned short int data[7];
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
  float Advance;
  float XOrig,YOrig;
  int Next,Prev,HashNext,HashPrev;
  CharFngrprnt Fngrprnt;
  float extent[2]; /* texture extent */
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
int CharacterGetGeometry(PyMOLGlobals *G,int id,
                         int *width, int *height, 
                         float *xorig, float *yorig, float *advance);

int CharacterNewFromBitmap(PyMOLGlobals *G,int width, int height,
                           unsigned char *bitmap,
                           float x_orig, float y_orig, float advance,
                           CharFngrprnt *fprnt,int sampling);

int CharacterNewFromBytemap(PyMOLGlobals *G, int width, int height,
                            int pitch, unsigned char *bytemap, 
                            float x_orig, float y_orig, float advance,
                            CharFngrprnt *fprnt);

int CharacterFind(PyMOLGlobals *G,CharFngrprnt *fprnt);

float CharacterInterpolate(PyMOLGlobals *G,int id,float *v);
void CharacterSetRetention(PyMOLGlobals *G,int retail_all);
unsigned char *CharacterGetPixmapBuffer(PyMOLGlobals *G,int id);
void CharacterRenderOpenGL(PyMOLGlobals *G,RenderInfo *info,int id);

#endif

