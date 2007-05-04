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
#include"PyMOLGlobals.h"

#define cLockAPI 1
#define cLockInbox 2
#define cLockOutbox 3

#define cPLog_pml_lf    0
#define cPLog_pml       1
#define cPLog_pym       2
#define cPLog_no_flush  3

#ifdef _PYMOL_NOPY

#define PRunStringInstance(x)
#define PRunStringModule(x)

#define PAutoBlock(G) 1
#define PAutoUnblock(G,a)

#define PBlock(G)
#define PUnblock(G)

#define PBlockLegacy()
#define PUnblockLegacy()

#define PLockAPIAsGlut(G,block_if_busy)
#define PUnlockAPIAsGlut(G)
#define PUnlockAPIAsGlutNoFlush(G)

#define PLockStatus()
#define PLockStatusAttempt() 1
#define PUnlockStatus()

#define PBlockAndUnlockAPI(G)
#define PLockAPIAndUnblock(G)

#define PFlush()
#define PFlushFast()
#define PParse(G,s)
#define PDo(G,s)

#define PLog(G,a,b)
#define PLogFlush(G)

#define PIsGlutThread() 1
#define PComplete(G,a,b) 0

#define PSGIStereo(a)
#define PPovrayRender(a,b,c,d,e,f) 0

#define PTruthCallStr(a,b,c)

#define PSleep(G,a)
#define PSleepWhileBusy(G,a)
#define PSleepUnlocked(G,a)

#define PFree()
#define PInit(G,a)
#define PSetupEmbedded(G,a,b)
#define PGetOptions(a)

#define PAlterAtom(G,a,b,c,d,e,f) 0
#define PLabelAtom(G,a,b,c,d) 0
#define PAlterAtomState(G,a,b,c,d,e,f,g) 0

#else

void PInit(PyMOLGlobals *G,int global_instance);
void PSetupEmbedded(PyMOLGlobals *G,int argc,char **argv);

  struct PyMOLOptionRec;

void PGetOptions(CPyMOLOptions *rec);

void PFree(void);
void PExit(PyMOLGlobals *G,int code);
void PParse(PyMOLGlobals *G,char *str); /* only accepts one command */
void PDo(PyMOLGlobals *G,char *str); /* accepts multple commands seperated by newlines */


int PAlterAtom(PyMOLGlobals *G,AtomInfoType *at,char *expr,
               int read_only,char *model,int index,PyObject *space);
int PLabelAtom(PyMOLGlobals *G,AtomInfoType *at,char *model,char *expr,int index);
int PAlterAtomState(PyMOLGlobals *G,float *v,char *expr,int read_only,
                    AtomInfoType *at,char *model, int index,PyObject *space);


void PLog(PyMOLGlobals *G,char *str,int lf);
void PLogFlush(PyMOLGlobals *G);

void PSleep(PyMOLGlobals *G,int usec);
void PSleepWhileBusy(PyMOLGlobals *G,int usec);
void PSleepUnlocked(PyMOLGlobals *G,int usec);

int PLockAPIAsGlut(PyMOLGlobals *G,int block_if_busy);
void PUnlockAPIAsGlut(PyMOLGlobals *G);
void PUnlockAPIAsGlutNoFlush(PyMOLGlobals *G);

void PLockStatus(void);
int  PLockStatusAttempt(void);
void PUnlockStatus(void);

void PBlock(PyMOLGlobals *G);
void PUnblock(PyMOLGlobals *G);

void PBlockLegacy(void);
void PUnblockLegacy(void);

int PAutoBlock(PyMOLGlobals *G);
void PAutoUnblock(PyMOLGlobals *G,int flag);

void PBlockAndUnlockAPI(PyMOLGlobals *G);
void PLockAPIAndUnblock(PyMOLGlobals *G);

void PFlush(PyMOLGlobals *G);
void PFlushFast(PyMOLGlobals *G);
void PXDecRef(PyObject *obj);

void PSGIStereo(PyMOLGlobals *G,int flag);
void PDefineFloat(PyMOLGlobals *G,char *name,float value);

void PRunStringModule(PyMOLGlobals *G,char *str);
void PRunStringInstance(PyMOLGlobals *G,char *str);

void PDumpTraceback(PyObject *err);
void PDumpException(void);

int PComplete(PyMOLGlobals *G,char *str,int buf_size);

int PTruthCallStr(PyObject *object,char *method,char *argument);
int PTruthCallStr0(PyObject *object,char *method);
int PTruthCallStr1i(PyObject *object,char *method,int argument);
int PTruthCallStr4i(PyObject *object,char *method,int a1,int a2,int a3,int a4);
int PPovrayRender(PyMOLGlobals *G,char *header,char *inp,char *file,int width,int height,int antialias);
int PIsGlutThread(void);

PyObject *PGetFontDict(PyMOLGlobals *G,float size,int face,int style);

/* instance-specific Python object, containers, and closures */
struct _CP_inst {
  PyObject *obj;
  PyObject *dict; 
  PyObject *exec;
  PyObject *cmd; 
  PyObject *parse; /* parse closure */
  PyObject *complete; /* complete partial command / TAB action */
  PyObject *cmd_do;
};

/* PyObject *GetBondsDict(void); */

/* all of the following Python objects must be invariant global
   modules & module dictionaries for the application */

extern PyObject *P_menu; /* used by Menu */
extern PyObject *P_xray; /* used by Symmetry */
extern PyObject *P_chempy; /* used by CoordSet and Selector for construction of models */
extern PyObject *P_models; /* used by Selector for construction of models */
extern PyObject *P_setting; /* used by Setting.c */

/* unused? */
extern PyObject *P_embed; /* not set by PyMOL -- must be set by host context */

extern int P_glut_thread_keep_out;
extern unsigned int P_glut_thread_id;

#endif
#endif






