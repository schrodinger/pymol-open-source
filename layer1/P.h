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
#ifndef _H_P
#define _H_P

#include"os_python.h"
#include"AtomInfo.h"

#define cLockAPI 1
#define cLockInbox 2
#define cLockOutbox 3

#define cPLog_pml_lf    0
#define cPLog_pml       1
#define cPLog_pym       2
#define cPLog_no_flush  3

#ifdef _PYMOL_NOPY

#define PRunString(x)

#define PAutoBlock() 1
#define PAutoUnblock(a)

#define PBlock()
#define PUnblock()

#define PLockAPIAsGlut()
#define PUnlockAPIAsGlut()

#define PBlockAndUnlockAPI()
#define PLockAPIAndUnblock()

#define PFlush()
#define PFlushFast()
#define PParse(s)
#define PDo(s)

#define PLog(a,b)
#define PLogFlush()

#define PIsGlutThread() 1
#define PComplete(a,b) 0

#define PSGIStereo(a)
#define PPovrayRender(a,b,c,d,e,f) 0

#define PTruthCallStr(a,b,c)

#define PSleep(a)

#define PFree()
#define PInit(G)
#define PInitEmbedded(a,b)
#define PGetOptions(a)

#define PAlterAtom(a,b,c,d,e) 0
#define PLabelAtom(a,b,c) 0
#define PAlterAtomState(a,b,c,d,e,f) 0

#else

void PInit(PyMOLGlobals *G);
void PInitEmbedded(int argc,char **argv);

  struct PyMOLOptionRec;

void PGetOptions(CPyMOLOptions *rec);

void PFree(void);
void PExit(int code);
void PParse(char *str); /* only accepts one command */
void PDo(char *str); /* accepts multple commands seperated by newlines */


int PAlterAtom(AtomInfoType *at,char *expr,int read_only,char *model,int index);
int PLabelAtom(AtomInfoType *at,char *expr,int index);
int PAlterAtomState(float *v,char *expr,int read_only,AtomInfoType *at,char *model, int index);


void PLog(char *str,int lf);
void PLogFlush(void);

void PSleep(int usec);

void PLockAPIAsGlut(void);
void PUnlockAPIAsGlut(void);

void PBlock(void);
void PUnblock(void);

int PAutoBlock(void);
void PAutoUnblock(int flag);

void PBlockAndUnlockAPI(void);
void PLockAPIAndUnblock(void);

void PFlush(void);
void PFlushFast(void);
void PXDecRef(PyObject *obj);

void PSGIStereo(int flag);
void PDefineFloat(char *name,float value);

void PRunString(char *str);
void PDumpTraceback(PyObject *err);
void PDumpException(void);

int PComplete(char *str,int buf_size);

int PTruthCallStr(PyObject *object,char *method,char *argument);
int PTruthCallStr0(PyObject *object,char *method);
int PTruthCallStr1i(PyObject *object,char *method,int argument);
int PTruthCallStr4i(PyObject *object,char *method,int a1,int a2,int a3,int a4);
int PPovrayRender(char *header,char *inp,char *file,int width,int height,int antialias);
int PIsGlutThread(void);

PyObject *PGetFontDict(float size,int face,int style);
PyObject *GetBondsDict(void);

extern PyObject *P_globals; /* used by main */

extern PyObject *P_cmd; /* used by Ray and main */
extern PyObject *P_menu; /* used by Menu */
extern PyObject *P_xray; /* used by Symmetry */
extern PyObject *P_chempy; /* used by CoordSet and Selector for construction of models */
extern PyObject *P_models; /* used by Selector for construction of models */
extern PyObject *P_setting; /* used by Setting.c */
extern PyObject *P_embed; /* not set by PyMOL -- must be set by host context */

extern int P_glut_thread_keep_out;
extern unsigned int P_glut_thread_id;

#endif
#endif






