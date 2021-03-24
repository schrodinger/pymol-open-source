
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
#include"ObjectMolecule.h"
#include"CoordSet.h"
#include"PyMOLGlobals.h"

#include "pymol/zstring_view.h"

#define cLockAPI 1
#define cLockInbox 2
#define cLockOutbox 3

#define cPLog_pml_lf    0
#define cPLog_pml       1
#define cPLog_pym       2
#define cPLog_no_flush  3

#define cPType_string          1
#define cPType_int             2
#define cPType_int_as_string   3
#define cPType_float           4
#define cPType_stereo          5
#define cPType_char_as_type    6
#define cPType_model           7
#define cPType_index           8
#define cPType_int_custom_type 9
#define cPType_xyz_float      10
#define cPType_settings       11
#define cPType_properties     12
#define cPType_state          13
#define cPType_schar          14
#define cPType_uint32         15

#define NUM_ATOM_PROPERTIES    41

#define cPRunType_alter          1
#define cPRunType_alter_state    2
#define cPRunType_label          3


int PLabelExprUsesVariable(PyMOLGlobals * G, const char *expr, const char *var);

int PLabelAtomAlt(PyMOLGlobals * G, AtomInfoType * at, const char *model, const char *expr,
                  int index);

#ifdef _PYMOL_NOPY

#define PRunStringInstance(G,x)
#define PRunStringModule(G,x)

#define PAutoBlock(G) 1
#define PAutoUnblock(G,a)

#define PBlock(G)
#define PUnblock(G)

#define PLockAPIAsGlut(G,block_if_busy)
#define PUnlockAPIAsGlut(G)
#define PUnlockAPIAsGlutNoFlush(G)

#define PLockAPI(G)
#define PUnlockAPI(G)

#define PLockStatus(G)
#define PLockStatusAttempt()G 1
#define PUnlockStatusG()

#define PBlockAndUnlockAPI(G)
#define PBlockAndUnlockAPI(G)
#define PLockAPIAndUnblock(G)
#define PTryLockAPIAndUnblock(G)
#define PFlush(G)
#define PFlushFast(G)
#define PParse(G,s)
#define PDo(G,s)

#define PLog(G,a,b)
#define PLogFlush(G)

#define PIsGlutThread() 1
#define PComplete(G,a,b) 0

#define PPovrayRender(a,b,c,d,e,f,g) 0

#define PTruthCallStr(a,b,c)

#define PSleep(G,a)
#define PSleepWhileBusy(G,a)
#define PSleepUnlocked(G,a)

#define PFree(G)
#define PInit(G,a)
#define PSetupEmbedded(G,a,b)
#define PConvertOptions(a,b)
#define PGetOptions(a)

#define PAlterAtom(G,a,b,c,d,e,f,g,h) 0
#define PLabelAtom(G,a,b,c,d,e,f) 0

#define PAlterAtomState(G,a,b,c,d,e,f,g,h,i,j,k) 0

#else

ov_status PCacheSet(PyMOLGlobals * G, PyObject * entry, PyObject * output);
ov_status PCacheGet(PyMOLGlobals * G,
                    PyObject ** result_output, PyObject ** result_entry,
                    PyObject * input);

void PInit(PyMOLGlobals * G, int global_instance);
void PSetupEmbedded(PyMOLGlobals * G, int argc, char **argv);

struct PyMOLOptionRec;

void PConvertOptions(CPyMOLOptions * rec, PyObject * options);
void PGetOptions(CPyMOLOptions * rec);

void PFree(PyMOLGlobals * G);
void PExit(PyMOLGlobals * G, int code);
void PParse(PyMOLGlobals * G, pymol::zstring_view str_view);       /* only accepts one command */
void PDo(PyMOLGlobals * G, const char *str);  /* accepts multple commands seperated by newlines */

int PAlterAtom(PyMOLGlobals * G, ObjectMolecule *obj, CoordSet *cs, PyObject *expr_co,
               int read_only, int atm, PyObject * space);
int PLabelAtom(PyMOLGlobals * G, ObjectMolecule *obj, CoordSet *cs, PyObject *expr_co, int atm);
int PAlterAtomState(PyMOLGlobals * G, PyObject *expr_co, int read_only,
                    ObjectMolecule *obj, CoordSet *cs, int atm, int idx,
                    int state, PyObject * space);

void PLog(PyMOLGlobals * G, pymol::zstring_view str, int lf);
void PLogFlush(PyMOLGlobals * G);

void PSleep(PyMOLGlobals * G, int usec);
void PSleepWhileBusy(PyMOLGlobals * G, int usec);
void PSleepUnlocked(PyMOLGlobals * G, int usec);

int PLockAPI(PyMOLGlobals * G, int block_if_busy);
void PUnlockAPI(PyMOLGlobals * G);

int PLockAPIAsGlut(PyMOLGlobals * G, int block_if_busy);
void PUnlockAPIAsGlut(PyMOLGlobals * G);
void PUnlockAPIAsGlutNoFlush(PyMOLGlobals * G);

void PLockStatus(PyMOLGlobals * G);
int PLockStatusAttempt(PyMOLGlobals * G);
void PUnlockStatus(PyMOLGlobals * G);

void PBlock(PyMOLGlobals * G);
void PUnblock(PyMOLGlobals * G);

int PAutoBlock(PyMOLGlobals * G);
void PAutoUnblock(PyMOLGlobals * G, int flag);

void PBlockAndUnlockAPI(PyMOLGlobals * G);
void PLockAPIAndUnblock(PyMOLGlobals * G);
int PTryLockAPIAndUnblock(PyMOLGlobals * G);

int PFlush(PyMOLGlobals * G);
int PFlushFast(PyMOLGlobals * G);
void PXDecRef(PyObject * obj);
PyObject *PXIncRef(PyObject * obj);

/**
 * Like `Py_INCREF` but returns the input argument for convenience.
 * Does not accept NULL, unlike `PXIncRef`.
 */
inline PyObject* PIncRef(PyObject* obj)
{
  Py_INCREF(obj);
  return obj;
}

void PDefineFloat(PyMOLGlobals * G, const char *name, float value);

void PErrPrintIfOccurred(PyMOLGlobals*);

void PRunStringModule(PyMOLGlobals * G, const char *str);
void PRunStringInstance(PyMOLGlobals * G, const char *str);

void PDumpTraceback(PyObject * err);
void PDumpException(void);

int PComplete(PyMOLGlobals * G, char *str, int buf_size);

int PTruthCallStr(PyObject * object, const char *method, const char *argument);
int PTruthCallStr0(PyObject * object, const char *method);
int PTruthCallStr1i(PyObject * object, const char *method, int argument);
int PTruthCallStr1s(PyObject * object, const char *method, const char *argument);
int PTruthCallStr4i(PyObject * object, const char *method, int a1, int a2, int a3, int a4);
int PPovrayRender(PyMOLGlobals * G, const char *header, const char *inp, const char *file, int width,
                  int height, int antialias);
int PIsGlutThread(void);

PyObject *PGetFontDict(PyMOLGlobals * G, float size, int face, int style);

typedef struct {
  long id;
  PyThreadState *state;
} SavedThreadRec;

struct SettingPropertyWrapperObject;

struct WrapperObject : PyObject {
  //  PyObject* dict;
  ObjectMolecule *obj;
  CoordSet *cs;
  AtomInfoType *atomInfo;
  int atm;
  int idx;
  int state;
  short read_only; // set for PLabelAtom
  PyMOLGlobals * G;
  PyObject *dict;
  SettingPropertyWrapperObject* settingWrapperObject;
#ifdef _PYMOL_IP_PROPERTIES
  SettingPropertyWrapperObject* propertyWrapperObject;
#endif
};

struct SettingPropertyWrapperObject : PyObject {
  WrapperObject* wobj;
};

/* instance-specific Python object, containers, closures, and threads */

#define MAX_SAVED_THREAD ((PYMOL_MAX_THREADS)+3)

struct _CP_inst {
  /* instance-specific storage */

  PyObject *obj;
  PyObject *dict;
  PyObject *exec;
  PyObject *cmd;
  PyObject *parse;              /* parse closure */
  PyObject *complete;           /* complete partial command / TAB action */
  PyObject *cmd_do;
  PyObject *colortype;          /* backwards compatible iterate/alter color type */

  PyObject *cache;

  /* locks and threads */

  PyObject *lock;               /* API locks */
  PyObject *lock_attempt;
  PyObject *unlock;

  PyObject *lock_api_status;        /* status locks */
  PyObject *lock_api_glut;          /* GLUT locks */

  int glut_thread_keep_out;
  SavedThreadRec savedThread[MAX_SAVED_THREAD];
};

using unique_PyObject_ptr_auto_gil = std::unique_ptr<PyObject, pymol::pyobject_delete_auto_gil>;

namespace pymol
{

/**
 * Handy RAII applications for acquiring the GIL
 */

class pblock
{
  PyMOLGlobals* m_G{nullptr};

public:
  pblock(PyMOLGlobals* G) : m_G(G) { PBlock(m_G); }
  ~pblock() { PUnblock(m_G); }
};

class pautoblock
{
  PyMOLGlobals* m_G;
  int m_blocked{};

public:
  pautoblock(PyMOLGlobals* G) : m_G(G), m_blocked(PAutoBlock(m_G)) {}
  ~pautoblock() { PAutoUnblock(m_G, m_blocked); }
};

/**
 * Shares ownership of a managed Python object.
 * @param obj Borrowed reference whose ownership will be transferred
 * @return owning pointer to Python object
 */

unique_PyObject_ptr_auto_gil make_auto_gil(PyObject* obj);

} // namespace pymol

/**
 * Makes Deep Copy of PyObject
 * @param input source PyObject to be copied
 * @return new PyObject
 */

unique_PyObject_ptr_auto_gil PDeepCopy(PyMOLGlobals* G, PyObject* input);

/* PyObject *GetBondsDict(void); */


/* all of the following Python objects must be invariant global
   modules & module dictionaries for the application */

extern PyObject *P_menu;        /* used by Menu */
extern PyObject *P_xray;        /* used by Symmetry */
extern PyObject *P_chempy;      /* used by CoordSet and Selector for construction of models */
extern PyObject *P_models;      /* used by Selector for construction of models */
extern PyObject *P_setting;     /* used by Setting.c */
extern PyTypeObject *P_wrapper;     /* used by P.c for lazy-loading settings/properties/attributes */
extern PyObject *P_CmdException;    /* pymol.CmdException */
extern PyObject *P_QuietException;  /* pymol.parsing.CmdException */
extern PyObject *P_IncentiveOnlyException; /* pymol.IncentiveOnlyException */

#endif
#endif
