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
#ifndef _H_Color
#define _H_Color

#include"Rep.h"
#include"Vector.h"

typedef char ColorName[64];

typedef struct {
  ColorName Name;
  Vector3f Color;
} ColorRec;

typedef struct  {
  ColorRec *Color;
  int NColor;
} CColor;

void ColorInit(void);
void ColorFree(void);

int ColorGetIndex(char *name);
float *ColorGet(int index); /* pointer maybe invalid after creating a new color */
float *ColorGetNamed(char *name);
void ColorDef(char *name,float *v);
int ColorGetNColor(void);
char *ColorGetName(int index);
int ColorGetStatus(int index);

#endif

