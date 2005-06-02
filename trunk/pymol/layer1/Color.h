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
#include"PyMOLGlobals.h"

#include"OVLexicon.h"
#include"OVOneToOne.h"

typedef char ColorName[24];

#define cColorGadgetRamp  1

#define cColorNewAuto     -2
#define cColorCurAuto     -3
#define cColorAtomic      -4

typedef struct {
  ColorName Name;
  Vector3f Color,Clamped;
  char ClampedFlag;
  char Custom;
} ColorRec;

typedef struct {
  ColorName Name;
  void *Ptr;
  int Type;
} ExtRec;

struct _CColor {
  ColorRec *Color;
  int NColor;
  ExtRec *Ext;
  int NExt;
  unsigned int *ColorTable;
  int BigEndian;
  OVLexicon *Lex;
  OVOneToOne *Idx;
};

int ColorInit(PyMOLGlobals *G);
void ColorFree(PyMOLGlobals *G);

int ColorGetNext(PyMOLGlobals *G);
int ColorGetCurrent(PyMOLGlobals *G);
int ColorGetIndex(PyMOLGlobals *G,char *name);
float *ColorGet(PyMOLGlobals *G,int index); /* pointer maybe invalid after creating a new color */
float *ColorGetNamed(PyMOLGlobals *G,char *name);
void ColorDef(PyMOLGlobals *G,char *name,float *v);
int ColorGetNColor(PyMOLGlobals *G);
char *ColorGetName(PyMOLGlobals *G,int index);
int ColorGetStatus(PyMOLGlobals *G,int index);
void ColorReset(PyMOLGlobals *G);

int ColorGetRamped(PyMOLGlobals *G,int index,float *vertex,float *color);
int ColorCheckRamped(PyMOLGlobals *G,int index);

struct ObjectGadgetRamp* ColorGetRamp(PyMOLGlobals *G,int index);


void ColorRegisterExt(PyMOLGlobals *G,char *name,void *extPtr,int type);
void ColorForgetExt(PyMOLGlobals *G,char *name);

PyObject *ColorAsPyList(PyMOLGlobals *G);
int ColorFromPyList(PyMOLGlobals *G,PyObject *list);

int ColorExtFromPyList(PyMOLGlobals *G,PyObject *list);
PyObject *ColorExtAsPyList(PyMOLGlobals *G);
int ColorTableLoad(PyMOLGlobals *G,char *fname,int quiet);
void ColorUpdateClamp(PyMOLGlobals *G,int index);
void ColorGetBkrdContColor(PyMOLGlobals *G,float *rgb, int invert_flag);
unsigned int ColorGet32BitWord(PyMOLGlobals *G,float *rgba);

#endif

