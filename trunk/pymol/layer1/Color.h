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

#include"os_python.h"

#include"Rep.h"
#include"Vector.h"

typedef char ColorName[64];

#define cColorGadgetRamp 1

typedef struct {
  ColorName Name;
  Vector3f Color,Clamped;
  int ClampedFlag;
  int Custom;
} ColorRec;

typedef struct {
  ColorName Name;
  void *Ptr;
  int Type;
} ExtRec;

typedef struct  {
  ColorRec *Color;
  int NColor;
  ExtRec *Ext;
  int NExt;
  unsigned int *ColorTable;
  int BigEndian;
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
void ColorReset(void);

int ColorGetRamped(int index,float *vertex,float *color);
int ColorCheckRamped(int index);

void ColorRegisterExt(char *name,void *extPtr,int type);
void ColorForgetExt(char *name);

PyObject *ColorAsPyList(void);
int ColorFromPyList(PyObject *list);

int ColorExtFromPyList(PyObject *list);
PyObject *ColorExtAsPyList(void);
int ColorTableLoad(char *fname,int quiet);
void ColorUpdateClamp(int index);

#endif

