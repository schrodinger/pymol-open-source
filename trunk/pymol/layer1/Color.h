
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
#include"Word.h"
#include"OVLexicon.h"
#include"OVOneToOne.h"

#define cColorGadgetRamp  1

#define cColorDefault     -1
#define cColorNewAuto     -2
#define cColorCurAuto     -3
#define cColorAtomic      -4
#define cColorObject      -5
#define cColorFront       -6
#define cColorBack        -7

#define cColorExtCutoff (-10)

#define cColor_TRGB_Bits  0x40000000
#define cColor_TRGB_Mask  0xC0000000

typedef struct {
  int Name;
  Vector3f Color, LutColor;
  char LutColorFlag;
  char Custom, Fixed;
  /* not saved */
  int old_session_index;
} ColorRec;

typedef struct {
  int Name;
  void *Ptr;
  int Type;
  /* not saved */
  int old_session_index;
} ExtRec;

struct _CColor {
  ColorRec *Color;
  int NColor;
  ExtRec *Ext;
  int NExt;
  int LUTActive;
  unsigned int *ColorTable;
  float Gamma;
  int BigEndian;
  OVLexicon *Lex;
  OVOneToOne *Idx;
  float RGBColor[3];            /* save global float for returning (float*) */
  char RGBName[10];
  /* not stored */
  int HaveOldSessionColors;
  int HaveOldSessionExtColors;
  float Front[3], Back[3];
};

int ColorInit(PyMOLGlobals * G);
void ColorFree(PyMOLGlobals * G);

int ColorGetNext(PyMOLGlobals * G);
int ColorGetCurrent(PyMOLGlobals * G);
int ColorGetIndex(PyMOLGlobals * G, char *name);
int ColorConvertOldSessionIndex(PyMOLGlobals * G, int index);
void ColorUpdateFront(PyMOLGlobals * G, float *back);

float *ColorGet(PyMOLGlobals * G, int index);   /* pointer maybe invalid after creating a new color */
float *ColorGetRaw(PyMOLGlobals * G, int index);        /* pointer maybe invalid after creating a new color */

float *ColorGetSpecial(PyMOLGlobals * G, int index);
float *ColorGetNamed(PyMOLGlobals * G, char *name);
void ColorDef(PyMOLGlobals * G, char *name, float *v, int mode, int quiet);
int ColorGetNColor(PyMOLGlobals * G);
char *ColorGetName(PyMOLGlobals * G, int index);
int ColorGetStatus(PyMOLGlobals * G, int index);
void ColorReset(PyMOLGlobals * G);

int ColorGetRamped(PyMOLGlobals * G, int index, float *vertex, float *color, int state);
int ColorCheckRamped(PyMOLGlobals * G, int index);

struct ObjectGadgetRamp *ColorGetRamp(PyMOLGlobals * G, int index);

void ColorRegisterExt(PyMOLGlobals * G, char *name, void *extPtr, int type);
void ColorForgetExt(PyMOLGlobals * G, char *name);

PyObject *ColorAsPyList(PyMOLGlobals * G);
int ColorFromPyList(PyMOLGlobals * G, PyObject * list, int partial_restore);

int ColorExtFromPyList(PyMOLGlobals * G, PyObject * list, int partial_restore);
PyObject *ColorExtAsPyList(PyMOLGlobals * G);
int ColorTableLoad(PyMOLGlobals * G, char *fname, float gamma, int quiet);
void ColorUpdateFromLut(PyMOLGlobals * G, int index);
int ColorLookupColor(PyMOLGlobals * G, float *color);
void ColorGetBkrdContColor(PyMOLGlobals * G, float *rgb, int invert_flag);
unsigned int ColorGet32BitWord(PyMOLGlobals * G, float *rgba);
int ColorGetEncoded(PyMOLGlobals * G, int index, float *color);

#endif
