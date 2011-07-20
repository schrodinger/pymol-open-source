
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


/* meaning of defines 

_PYMOL_MONOLITHIC: means that we're building PyMOL and its Python C
dependencies as one C library.  That means we need to explicitly call
the initialization functions for these libraries on startup.

*/
#include"os_python.h"

#include"os_predef.h"
#include"os_std.h"
#include"Base.h"


/* BEGIN PROPRIETARY CODE SEGMENT (see disclaimer in "os_proprietary.h") */
#ifdef WIN32
#include<windows.h>
#include<process.h>
#ifndef WIN64
#if 0
#include<winappc.h>
#endif
#endif
#endif

/* END PROPRIETARY CODE SEGMENT */

#ifndef _PYMOL_NOPY
#ifdef _PYMOL_MINGW
#define putenv _putenv
#endif
#include"os_python.h"
#endif

#include"os_std.h"
#include"os_time.h"
#include"os_unix.h"

#include"MemoryDebug.h"
#include"Base.h"
#include"Err.h"
#include"P.h"
#include"PConv.h"
#include"Ortho.h"
#include"Cmd.h"
#include"main.h"
#include"AtomInfo.h"
#include"CoordSet.h"
#include"Util.h"
#include"Executive.h"
#include"PyMOLOptions.h"
#include"PyMOL.h"

static int label_copy_text(char *dst, char *src, int len, int max)
{
  dst += len;
  while(len < max) {
    if(!*src)
      break;
    *(dst++) = *(src++);
    len++;
  }
  *dst = 0;
  return len;
}

static int label_next_token(WordType dst, char **expr)
{
  char *p = *expr;
  char *q = dst;
  char ch;
  int tok_len = 0;
  int tok_max = sizeof(WordType) - 1;

  /* skip leading whitespace (if any) */

  while((ch = *p)) {
    if(ch > 33)
      break;
    p++;
  }

  /* copy the token */

  while((ch = *p)) {
    if(((ch >= 'a') && (ch <= 'z')) ||
       ((ch >= 'A') && (ch <= 'Z')) || ((ch >= '0') && (ch <= '9')) || ((ch == '_'))) {
      if(tok_len < tok_max) {
        *(q++) = ch;
        tok_len++;
      }
    } else {
      break;
    }
    p++;
  }
  *q = 0;
  if(p != *expr)
    *expr = p;
  else if(*p)
    *expr = p + 1;              /* always advance input by at least one character */

  /* let caller know whether we read anything */

  return (q != dst);
}

/* MMSTEREO_STALE is how this mmstereo field gets intialized, and it means
   that mmstereo information needs to be computed */
#define MMSTEREO_STALE 0                  // ' '
#define MMSTEREO_NO_CHIRALITY 127         // ' '
#define PYMOL_MMSTEREO_CHIRALITY_R 1      // 'R'
#define PYMOL_MMSTEREO_CHIRALITY_S 2      // 'S'
#define PYMOL_MMSTEREO_E 11               // 'E'
#define PYMOL_MMSTEREO_Z 12               // 'Z'
#define PYMOL_MMSTEREO_P 13               // 'P'
#define PYMOL_MMSTEREO_M 14               // 'M'
#define PYMOL_MMSTEREO_CHIRALITY_ANR 50   // 'r'
#define PYMOL_MMSTEREO_CHIRALITY_ANS 51   // 's'
#define PYMOL_MMSTEREO_UNDEF 99
#define PYMOL_MMSTEREO_ST_INDEFINITE 100
#define PYMOL_MMSTEREO_GEOM_INDEFINITE 101
#define PYMOL_MMSTEREO_AN_GEOM_INDEFINITE 102

char convertStereoToChar(int stereo){
  switch (stereo){
  case PYMOL_MMSTEREO_CHIRALITY_R:
    return 'R';
  case PYMOL_MMSTEREO_CHIRALITY_S:
    return 'S';
  case  PYMOL_MMSTEREO_E:
    return 'E';
  case  PYMOL_MMSTEREO_Z:
    return 'Z';
  case  PYMOL_MMSTEREO_P:
    return 'P';
  case  PYMOL_MMSTEREO_M:
    return 'M';
  case PYMOL_MMSTEREO_CHIRALITY_ANR:
    return 'r';
  case PYMOL_MMSTEREO_CHIRALITY_ANS:
    return 's';
  case PYMOL_MMSTEREO_UNDEF:
  case PYMOL_MMSTEREO_ST_INDEFINITE:
  case PYMOL_MMSTEREO_GEOM_INDEFINITE:
  case PYMOL_MMSTEREO_AN_GEOM_INDEFINITE:
    return '?';
  }
  return ' ';
}
int convertCharToStereo(char stereo){
  switch (stereo){
  case 'R':
    return PYMOL_MMSTEREO_CHIRALITY_R;
  case 'S':
    return PYMOL_MMSTEREO_CHIRALITY_S;
  case 'r':
    return PYMOL_MMSTEREO_CHIRALITY_ANR;
  case 's':
    return PYMOL_MMSTEREO_CHIRALITY_ANS;
  case 'E':
    return PYMOL_MMSTEREO_E;
  case 'Z':
    return PYMOL_MMSTEREO_Z;
  case 'P':
    return PYMOL_MMSTEREO_P;
  case 'M':
    return PYMOL_MMSTEREO_M;
  case '?':
    return PYMOL_MMSTEREO_GEOM_INDEFINITE;
  }
  return MMSTEREO_NO_CHIRALITY;
}


int PLabelExprUsesVariable(PyMOLGlobals * G, char *expr, char *var)
{
  OrthoLineType buffer;
  char ch, quote = 0;
  int escaped = false;
  while((ch = *(expr++))) {
    if(!quote) {
      if(ch == '\'') {
        quote = ch;
      } else if(ch == '"') {
        quote = ch;
      } else if((ch < 33) || (ch == '+') || (ch == '(') || (ch == ')')) {
        /* nop */
      } else if(ch > 32) {
        WordType tok;
        expr--;
        if(label_next_token(tok, &expr)) {
          buffer[0] = 0;
          if(!strcmp(tok, var)) {
	    return 1;
	  }
	}
      }
    } else {
      if(ch == quote) {
        quote = 0;
      } else if(ch == '\\') {
        if(!escaped) {
          escaped = true;
        } else {
          escaped = false;
        }
      }
    }
  }
  return 0;
}

int PLabelAtomAlt(PyMOLGlobals * G, AtomInfoType * at, char *model, char *expr, int index)
{
  /* alternate C implementation which bypasses Python expressions -- works
     only for simple label formats "..."+property+... */

  int result = true;
  OrthoLineType label;
  int label_len = 0;
  int label_max = sizeof(OrthoLineType);
  OrthoLineType buffer;
  char ch, quote = 0;
  int escaped = false;

  label[0] = 0;
  while((ch = *(expr++))) {
    if(!quote) {
      if(ch == '\'') {
        quote = ch;
      } else if(ch == '"') {
        quote = ch;
      } else if((ch < 33) || (ch == '+') || (ch == '(') || (ch == ')')) {
        /* nop */
      } else if(ch > 32) {
        WordType tok;
        expr--;
        if(label_next_token(tok, &expr)) {
          /* brain-dead linear string matching */
          buffer[0] = 0;
          if(!strcmp(tok, "model")) {
            label_len = label_copy_text(label, model, label_len, label_max);
          } else if(!strcmp(tok, "index")) {
            sprintf(buffer, "%d", index + 1);
          } else if(!strcmp(tok, "type")) {
            if(at->hetatm)
              label_len = label_copy_text(label, "HETATM", label_len, label_max);
            else
              label_len = label_copy_text(label, "ATOM", label_len, label_max);
          } else if(!strcmp(tok, "name")) {
            label_len = label_copy_text(label, at->name, label_len, label_max);
          } else if(!strcmp(tok, "resn")) {
            label_len = label_copy_text(label, at->resn, label_len, label_max);
          } else if(!strcmp(tok, "resi")) {
            label_len = label_copy_text(label, at->resi, label_len, label_max);
          } else if(!strcmp(tok, "resv")) {
            sprintf(buffer, "%d", at->resv);
          } else if(!strcmp(tok, "chain")) {
            label_len = label_copy_text(label, at->chain, label_len, label_max);
          } else if(!strcmp(tok, "alt")) {
            label_len = label_copy_text(label, at->alt, label_len, label_max);
          } else if(!strcmp(tok, "segi")) {
            label_len = label_copy_text(label, at->segi, label_len, label_max);
          } else if(!strcmp(tok, "ss")) {
            label_len = label_copy_text(label, at->ssType, label_len, label_max);
          } else if(!strcmp(tok, "vdw")) {
            sprintf(buffer, "%1.2f", at->vdw);
          } else if(!strcmp(tok, "elec_radius")) {
            sprintf(buffer, "%1.2f", at->elec_radius);
          } else if(!strcmp(tok, "text_type")) {
            char null_st[1] = "";
            char *st = null_st;
            if(at->textType)
              st = OVLexicon_FetchCString(G->Lexicon, at->textType);
            label_len = label_copy_text(label, st, label_len, label_max);
          } else if(!strcmp(tok, "elem")) {
            label_len = label_copy_text(label, at->elem, label_len, label_max);
          } else if(!strcmp(tok, "geom")) {
            sprintf(buffer, "%d", at->geom);
          } else if(!strcmp(tok, "valence")) {
            sprintf(buffer, "%d", at->valence);
          } else if(!strcmp(tok, "rank")) {
            sprintf(buffer, "%d", at->rank);
          } else if(!strcmp(tok, "flags")) {
            if(at->flags) {
              sprintf(buffer, "%X", at->flags);
            } else {
              strcpy(buffer, "0");
            }
          } else if(!strcmp(tok, "q")) {
            sprintf(buffer, "%1.2f", at->q);
          } else if(!strcmp(tok, "b")) {
            sprintf(buffer, "%1.2f", at->b);
          } else if(!strcmp(tok, "numeric_type")) {
            if(at->customType != cAtomInfoNoType)
              sprintf(buffer, "%d", at->customType);
            else {
              strcpy(buffer, "?");
            }
          } else if(!strcmp(tok, "partial_charge")) {
            sprintf(buffer, "%1.3f", at->partialCharge);
          } else if(!strcmp(tok, "formal_charge")) {
            sprintf(buffer, "%d", at->formalCharge);
          } else if(!strcmp(tok, "stereo")) {
	    sprintf(buffer, "%c", convertStereoToChar(at->mmstereo));
          } else if(!strcmp(tok, "color")) {
            sprintf(buffer, "%d", at->color);
          } else if(!strcmp(tok, "cartoon")) {
            sprintf(buffer, "%d", at->cartoon);
          } else if(!strcmp(tok, "ID")) {
            sprintf(buffer, "%d", at->id);
          } else if(!strcmp(tok, "str")) {
            /* nop */
          }
          if(buffer[0]) {
            label_len = label_copy_text(label, buffer, label_len, label_max);
          }
        } else {
          label_len = label_copy_text(label, "?", label_len, label_max);
          label_len = label_copy_text(label, tok, label_len, label_max);
        }
      } else {
        if(label_len < label_max) {
          label[label_len] = '?';
          label_len++;
        }
      }
    } else {
      if(ch == quote) {
        quote = 0;
      } else if(ch == '\\') {
        if(!escaped) {
          escaped = true;
        } else {
          if(label_len < label_max) {
            label[label_len] = ch;
            label_len++;
          }
          escaped = false;
        }
      } else {
        if(label_len < label_max) {
          label[label_len] = ch;
          label_len++;
          label[label_len] = 0;
        }
      }
    }
  }

  if(result) {
    if(at->label) {
      OVLexicon_DecRef(G->Lexicon, at->label);
      at->label = 0;
    }

    if(label[0]) {
      OVreturn_word ret = OVLexicon_GetFromCString(G->Lexicon, label);
      if(OVreturn_IS_OK(ret)) {
        at->label = ret.word;
      }
    }
  } else {
    ErrMessage(G, "Label", "Aborting on error. Labels may be incomplete.");
  }
  return (result);
}

#ifndef _PYMOL_NOPY


/* all of the following Python objects must be invariant & global for the application */


/* these are module / module properties -- global and static for a given interpreter */


/* local to this C code module */

static PyObject *P_pymol = NULL;
static PyObject *P_pymol_dict = NULL;   /* must be refomed into globals and instance properties */
static PyObject *P_cmd = NULL;

static PyObject *P_povray = NULL;
static PyObject *P_traceback = NULL;
static PyObject *P_parser = NULL;

static PyObject *P_main = NULL;
static PyObject *P_vfont = NULL;


/* used elsewhere */

PyObject *P_menu = NULL;        /* menu definitions are currently global */
PyObject *P_xray = NULL;        /* okay as global */
PyObject *P_chempy = NULL;      /* okay as global */
PyObject *P_models = NULL;      /* okay as global */
PyObject *P_setting = NULL;     /* okay as global -- just used for names */


/* BEGIN PROPRIETARY CODE SEGMENT (see disclaimer in "os_proprietary.h") */
#ifdef WIN32
static PyObject *P_time = NULL;
static PyObject *P_sleep = NULL;
#endif

/* END PROPRIETARY CODE SEGMENT */

static void PUnlockAPIWhileBlocked(PyMOLGlobals * G);
static void PLockAPIWhileBlocked(PyMOLGlobals * G);

#define P_log_file_str "_log_file"

#define xxxPYMOL_NEW_THREADS 1

unsigned int PyThread_get_thread_ident(void);   /* critical functionality */

void PRunStringModule(PyMOLGlobals * G, char *str);

void PLockStatus(PyMOLGlobals * G)
{                               /* assumes we have the GIL */
  PXDecRef(PyObject_CallFunction(G->P_inst->lock_status, "O", G->P_inst->cmd));
}

int PLockStatusAttempt(PyMOLGlobals * G)
{                               /* assumes we have the GIL */
  int result = true;
  PyObject *got_lock =
    PyObject_CallFunction(G->P_inst->lock_status_attempt, "O", G->P_inst->cmd);
  if(got_lock) {
    if(!PyInt_AsLong(got_lock)) {
      result = false;
    }
    Py_DECREF(got_lock);
  }
  return result;
}

void PUnlockStatus(PyMOLGlobals * G)
{                               /* assumes we have the GIL */
  PXDecRef(PyObject_CallFunction(G->P_inst->unlock_status, "O", G->P_inst->cmd));
}

static void PLockGLUT(PyMOLGlobals * G)
{                               /* assumes we have the GIL */
  PXDecRef(PyObject_CallFunction(G->P_inst->lock_glut, "O", G->P_inst->cmd));
}

static void PUnlockGLUT(PyMOLGlobals * G)
{                               /* assumes we have the GIL */
  PXDecRef(PyObject_CallFunction(G->P_inst->unlock_glut, "O", G->P_inst->cmd));
}

unsigned int P_glut_thread_id = -1;


/* enables us to keep glut out if by chance it grabs the API
 * in the middle of a nested API based operation */

void PCatchInit(void);
void my_interrupt(int a);
char *getprogramname(void);


/*
  PyObject *GetBondsDict(PyMOLGlobals *G)
  {
  PyObject *result = NULL;
  result = PyObject_GetAttrString(P_chempy,"bonds");
  if(!result) ErrMessage(G,"PyMOL","can't find 'chempy.bonds.bonds'");
  return(result);
  }
*/

PyObject *PGetFontDict(PyMOLGlobals * G, float size, int face, int style)
{                               /* assumes we have a valid interpreter lock */
  PyObject *result = NULL;

  if(!P_vfont) {
    PRunStringModule(G, "import vfont\n");
    P_vfont = PyDict_GetItemString(P_pymol_dict, "vfont");
    Py_XINCREF(P_vfont);
  }
  if(!P_vfont) {
    PRINTFB(G, FB_Python, FB_Errors)
      " PyMOL-Error: can't find module 'vfont'" ENDFB(G);
  } else {
    result = PyObject_CallMethod(P_vfont, "get_font", "fii", size, face, style);
  }
  return (PConvAutoNone(result));
}

int PComplete(PyMOLGlobals * G, char *str, int buf_size)
{
  int ret = false;
  PyObject *result;
  char *st2;
  PBlockAndUnlockAPI(G);
  if(G->P_inst->complete) {
    result = PyObject_CallFunction(G->P_inst->complete, "s", str);
    if(result) {
      if(PyString_Check(result)) {
        st2 = PyString_AsString(result);
        UtilNCopy(str, st2, buf_size);
        ret = true;
      }
      Py_DECREF(result);
    }
  }
  PLockAPIAndUnblock(G);
  return (ret);
}

int PTruthCallStr0(PyObject * object, char *method)
{
  int result = false;
  PyObject *tmp;
  tmp = PyObject_CallMethod(object, method, "");
  if(tmp) {
    if(PyObject_IsTrue(tmp))
      result = 1;
    Py_DECREF(tmp);
  }
  return (result);
}

int PTruthCallStr(PyObject * object, char *method, char *argument)
{
  int result = false;
  PyObject *tmp;
  tmp = PyObject_CallMethod(object, method, "s", argument);
  if(tmp) {
    if(PyObject_IsTrue(tmp))
      result = 1;
    Py_DECREF(tmp);
  }
  return (result);
}

int PTruthCallStr1i(PyObject * object, char *method, int argument)
{
  int result = false;
  PyObject *tmp;
  tmp = PyObject_CallMethod(object, method, "i", argument);
  if(tmp) {
    if(PyObject_IsTrue(tmp))
      result = 1;
    Py_DECREF(tmp);
  }
  return (result);
}

int PTruthCallStr1s(PyObject * object, char *method, char *argument)
{
  int result = false;
  PyObject *tmp;
  tmp = PyObject_CallMethod(object, method, "s", argument);
  if(tmp) {
    if(PyObject_IsTrue(tmp))
      result = 1;
    Py_DECREF(tmp);
  }
  return (result);
}

int PTruthCallStr4i(PyObject * object, char *method, int a1, int a2, int a3, int a4)
{
  int result = false;
  PyObject *tmp;
  tmp = PyObject_CallMethod(object, method, "iiii", a1, a2, a3, a4);
  if(tmp) {
    if(PyObject_IsTrue(tmp))
      result = 1;
    Py_DECREF(tmp);
  }
  return (result);
}

PyObject *PXIncRef(PyObject * obj)
{
  if(!obj)
    obj = Py_None;
  Py_XINCREF(obj);
  return obj;
}

void PXDecRef(PyObject * obj)
{
  Py_XDECREF(obj);
}

OV_STATIC ov_status CacheCreateEntry(PyObject ** result, PyObject * input)
{
  ov_status status = OV_STATUS_FAILURE;
  if(input && PyTuple_Check(input)) {
    ov_size tuple_size = PyTuple_Size(input);
    ov_size tot_size = tuple_size;
    PyObject *hash_code = PyTuple_New(tuple_size);
    PyObject *entry = PyList_New(6);
    if(hash_code && entry) {
      /* compute hash codes & total input size */
      ov_size i;
      status = OV_STATUS_SUCCESS;
      for(i = 0; i < tuple_size; i++) {
        PyObject *item = PyTuple_GetItem(input, i);
        long hash_long;
        if(item != Py_None) {
          /* here we are assuming that different Python versions will
           * hash tuples of ints & floats in exactly the same way (at least to 31 bits of significance) */
          hash_long = 0x7FFFFFFF & PyObject_Hash(item); /* pos 32 bit # to preserve 32-bit/64-bit compat */
        } else {
          hash_long = 0;        /* None doesn't hash consistently from Python version to version */
        }
        PyTuple_SetItem(hash_code, i, PyInt_FromLong(hash_long));
        if(PyTuple_Check(item)) {
          tot_size += PyTuple_Size(item);
        }
      }
      PyList_SetItem(entry, 0, PyInt_FromLong(tot_size));
      PyList_SetItem(entry, 1, hash_code);
      PyList_SetItem(entry, 2, PXIncRef(input));
      PyList_SetItem(entry, 3, PXIncRef(NULL));
      PyList_SetItem(entry, 4, PyInt_FromLong(0));      /* access count */
      PyList_SetItem(entry, 5, PyFloat_FromDouble(0.0));        /* timestamp */
    }
    if(!OV_OK(status)) {
      PXDecRef(hash_code);
      PXDecRef(entry);
    } else {
      *result = entry;
    }
  }
  if(PyErr_Occurred())
    PyErr_Print();
  return status;
}

ov_status PCacheSet(PyMOLGlobals * G, PyObject * entry, PyObject * output)
{
  ov_status status = OV_STATUS_FAILURE;
  if(G->P_inst->cache && output) {
    ov_size tuple_size = PyTuple_Size(output);
    ov_size tot_size = tuple_size + PyInt_AsLong(PyList_GetItem(entry, 0));
    status = OV_STATUS_SUCCESS;
    {
      ov_size i;
      for(i = 0; i < tuple_size; i++) {
        PyObject *item = PyTuple_GetItem(output, i);
        if(PyTuple_Check(item)) {
          tot_size += PyTuple_Size(item);
        }
      }
    }
    PyList_SetItem(entry, 0, PyInt_FromLong(tot_size)); /* update total size */
    PyList_SetItem(entry, 3, PXIncRef(output));
    PXDecRef(PyObject_CallMethod(G->P_inst->cmd, "_cache_set",
                                 "OiO", entry, SettingGetGlobal_i(G, cSetting_cache_max),
                                 G->P_inst->cmd));
    /* compute the hash codes */
  }
  if(PyErr_Occurred())
    PyErr_Print();
  return status;
}

ov_status PCacheGet(PyMOLGlobals * G,
                    PyObject ** result_output, PyObject ** result_entry, PyObject * input)
{
  ov_status status = OV_STATUS_NO;
  if(G->P_inst->cache) {
    PyObject *entry = NULL;
    PyObject *output = NULL;

    if(OV_OK(CacheCreateEntry(&entry, input))) {
      output = PyObject_CallMethod(G->P_inst->cmd, "_cache_get",
                                   "OOO", entry, Py_None, G->P_inst->cmd);
      if(output == Py_None) {
        Py_DECREF(output);
        output = NULL;
      } else {
        status = OV_STATUS_YES;
      }
    }
    /* compute the hash codes */
    if(OV_OK(status)) {
      *result_entry = entry;
      *result_output = output;
    } else {
      PXDecRef(entry);
      PXDecRef(output);
    }
  }
  if(PyErr_Occurred())
    PyErr_Print();
  return status;
}

void PSleepWhileBusy(PyMOLGlobals * G, int usec)
{
#ifndef WIN32
  struct timeval tv;
  PRINTFD(G, FB_Threads)
    " PSleep-DEBUG: napping.\n" ENDFD;
  tv.tv_sec = 0;
  tv.tv_usec = usec;
  select(0, NULL, NULL, NULL, &tv);
  PRINTFD(G, FB_Threads)
    " PSleep-DEBUG: nap over.\n" ENDFD;
#else
  /* BEGIN PROPRIETARY CODE SEGMENT (see disclaimer in "os_proprietary.h") */
  PBlock(G);
  PXDecRef(PyObject_CallFunction(P_sleep, "f", usec / 1000000.0));
  PUnblock(G);
  /* END PROPRIETARY CODE SEGMENT */
#endif
}

void PSleepUnlocked(PyMOLGlobals * G, int usec)
{                               /* can only be called by the glut process */
#ifndef WIN32
  struct timeval tv;
  PRINTFD(G, FB_Threads)
    " PSleep-DEBUG: napping.\n" ENDFD;
  tv.tv_sec = 0;
  tv.tv_usec = usec;
  select(0, NULL, NULL, NULL, &tv);
  PRINTFD(G, FB_Threads)
    " PSleep-DEBUG: nap over.\n" ENDFD;
#else
  /* BEGIN PROPRIETARY CODE SEGMENT (see disclaimer in "os_proprietary.h") */
  PBlock(G);
  PXDecRef(PyObject_CallFunction(P_sleep, "f", usec / 1000000.0));
  PUnblock(G);
  /* END PROPRIETARY CODE SEGMENT */
#endif
}

void PSleep(PyMOLGlobals * G, int usec)
{                               /* can only be called by the glut process */
#ifndef WIN32
  struct timeval tv;
  PUnlockAPIAsGlut(G);
  PRINTFD(G, FB_Threads)
    " PSleep-DEBUG: napping.\n" ENDFD;
  tv.tv_sec = 0;
  tv.tv_usec = usec;
  select(0, NULL, NULL, NULL, &tv);
  PRINTFD(G, FB_Threads)
    " PSleep-DEBUG: nap over.\n" ENDFD;
  PLockAPIAsGlut(G, true);
#else
  /* BEGIN PROPRIETARY CODE SEGMENT (see disclaimer in "os_proprietary.h") */
  PBlockAndUnlockAPI(G);
  PXDecRef(PyObject_CallFunction(P_sleep, "f", usec / 1000000.0));
  PLockAPIAndUnblock(G);
  /* END PROPRIETARY CODE SEGMENT */
#endif

}

static PyObject *PCatchWrite(PyObject * self, PyObject * args);

void my_interrupt(int a)
{
  exit(EXIT_FAILURE);
}

void PDumpTraceback(PyObject * err)
{
  PyObject_CallMethod(P_traceback, "print_tb", "O", err);
}

void PDumpException()
{
  PyObject_CallMethod(P_traceback, "print_exc", "");
}

int PAlterAtomState(PyMOLGlobals * G, float *v, char *expr, int read_only,
                    AtomInfoType * at, char *model, int index, PyObject * space)

/* assumes Blocked python interpreter */
{
  PyObject *dict;
  int result = true;
  float f[3];
  PyObject *x_id1, *x_id2 = NULL, *y_id1, *y_id2 = NULL, *z_id1, *z_id2 = NULL;
  char atype[7], mmstereotype[2];
  PyObject *flags_id1 = NULL, *flags_id2 = NULL;
  int flags;
  dict = PyDict_New();

  if(at) {
    if(at->hetatm)
      strcpy(atype, "HETATM");
    else
      strcpy(atype, "ATOM");

    /* immutables */
    PConvStringToPyDictItem(dict, "model", model);
    PConvIntToPyDictItem(dict, "index", index + 1);

    PConvStringToPyDictItem(dict, "type", atype);
    PConvStringToPyDictItem(dict, "name", at->name);
    PConvStringToPyDictItem(dict, "resn", at->resn);
    PConvStringToPyDictItem(dict, "resi", at->resi);
    PConvIntToPyDictItem(dict, "resv", at->resv);       /* subordinate to resi */
    PConvStringToPyDictItem(dict, "chain", at->chain);
    PConvStringToPyDictItem(dict, "alt", at->alt);
    PConvStringToPyDictItem(dict, "segi", at->segi);
    PConvStringToPyDictItem(dict, "elem", at->elem);
    PConvStringToPyDictItem(dict, "ss", at->ssType);

    {
      char null_st[1] = "";
      char *st = null_st;

      if(at->textType)
        st = OVLexicon_FetchCString(G->Lexicon, at->textType);
      PConvStringToPyDictItem(dict, "text_type", st);

      st = null_st;
      if(at->label)
        st = OVLexicon_FetchCString(G->Lexicon, at->label);
      PConvStringToPyDictItem(dict, "label", st);
    }
    PConvIntToPyDictItem(dict, "numeric_type", at->customType);
    PConvFloatToPyDictItem(dict, "q", at->q);
    PConvFloatToPyDictItem(dict, "b", at->b);
    PConvFloatToPyDictItem(dict, "vdw", at->vdw);
    PConvFloatToPyDictItem(dict, "elec_radius", at->elec_radius);
    PConvFloatToPyDictItem(dict, "partial_charge", at->partialCharge);
    PConvIntToPyDictItem(dict, "formal_charge", at->formalCharge);
    mmstereotype[0] = convertStereoToChar(at->mmstereo);
    mmstereotype[1] = 0;
    PConvStringToPyDictItem(dict, "stereo", mmstereotype);
    PConvIntToPyDictItem(dict, "cartoon", at->cartoon);
    PConvIntToPyDictItem(dict, "color", at->color);
    PConvIntToPyDictItem(dict, "ID", at->id);
    PConvIntToPyDictItem(dict, "rank", at->rank);

    /* mutables */
    flags_id1 = PConvIntToPyDictItem(dict, "flags", at->flags);
  }
  x_id1 = PConvFloatToPyDictItem(dict, "x", v[0]);
  y_id1 = PConvFloatToPyDictItem(dict, "y", v[1]);
  z_id1 = PConvFloatToPyDictItem(dict, "z", v[2]);
  PXDecRef(PyRun_String(expr, Py_single_input, space, dict));
  if(PyErr_Occurred()) {
    PyErr_Print();
    result = false;
  } else if(!read_only) {
    if(result) {
      if(!(x_id2 = PyDict_GetItemString(dict, "x")))
        result = false;
      else if(!(y_id2 = PyDict_GetItemString(dict, "y")))
        result = false;
      else if(!(z_id2 = PyDict_GetItemString(dict, "z")))
        result = false;
      else if(at) {
        if(!(flags_id2 = PyDict_GetItemString(dict, "flags")))
          result = false;
      }
      if(PyErr_Occurred()) {
        PyErr_Print();
        result = false;
        ErrMessage(G, "AlterState", "Aborting on error. Assignment may be incomplete.");
      }
    }
    if(result) {
      f[0] = (float) PyFloat_AsDouble(x_id2);
      f[1] = (float) PyFloat_AsDouble(y_id2);
      f[2] = (float) PyFloat_AsDouble(z_id2);
      if(at)
        if(flags_id1 != flags_id2) {
          if(!PConvPyObjectToInt(flags_id2, &flags))
            result = false;
          else
            at->flags = flags;
        }
      if(PyErr_Occurred()) {
        PyErr_Print();
        result = false;
        ErrMessage(G, "AlterState", "Aborting on error. Assignment may be incomplete.");
      } else {
        v[0] = f[0];
        v[1] = f[1];
        v[2] = f[2];
      }
    }

  }
  Py_DECREF(dict);
  return result;
}

int PAlterAtom(PyMOLGlobals * G,
               AtomInfoType * at, char *expr, int read_only,
               char *model, int index, PyObject * space)
{
  /* assumes Blocked python interpreter */
  WordType buf;
  AtomName name;
  PyObject *name_id1, *name_id2 = NULL;
  ElemName elem;
  PyObject *elem_id1, *elem_id2 = NULL;
  ResName resn;
  PyObject *resn_id1, *resn_id2 = NULL;
  ResIdent resi;
  PyObject *resi_id1, *resi_id2 = NULL;
  int resv;
  PyObject *resv_id1, *resv_id2 = NULL;
  Chain chain;
  PyObject *chain_id1, *chain_id2 = NULL;
  Chain alt;
  PyObject *alt_id1, *alt_id2 = NULL;
  SegIdent segi;
  PyObject *flags_id1, *flags_id2 = NULL;
  int flags;
  PyObject *segi_id1, *segi_id2 = NULL;
  PyObject *text_type_id1, *text_type_id2 = NULL;
  SSType ssType;
  PyObject *ss_id1, *ss_id2 = NULL;
  char atype[7], mmstereotype[2];
  PyObject *type_id1, *type_id2 = NULL;
  float b, q, partialCharge, vdw, elec_radius;
  PyObject *b_id1, *b_id2 = NULL;
  PyObject *q_id1, *q_id2 = NULL;
  PyObject *partial_charge_id1, *partial_charge_id2 = NULL;
  PyObject *vdw_id1, *vdw_id2 = NULL;
  PyObject *elec_radius_id1, *elec_radius_id2 = NULL;
  int formalCharge, numericType;
  PyObject *formal_charge_id1, *formal_charge_id2 = NULL;
  PyObject *numeric_type_id1, *numeric_type_id2 = NULL;
  char stereo[2];
  PyObject *stereo_id1, *stereo_id2 = NULL;
  int cartoon;
  PyObject *cartoon_id1, *cartoon_id2 = NULL;
  int color;
  PyObject *color_id1, *color_id2 = NULL;
  PyObject *label_id1, *label_id2 = NULL;
  int id;
  PyObject *ID_id1, *ID_id2 = NULL;
  int rank;
  PyObject *rank_id1, *rank_id2 = NULL;
  PyObject *state_id1, *state_id2 = NULL;
  int state;
  PyObject *dict;
  int result = true;

  if(at->hetatm)
    strcpy(atype, "HETATM");
  else
    strcpy(atype, "ATOM");

  /* PBlockAndUnlockAPI() is not safe, thus these
   * expressions must not call the PyMOL API...
   * what if "at" is destroyed by another thread? */

  dict = PyDict_New();

  /* immutables */
  PConvStringToPyDictItem(dict, "model", model);
  PConvIntToPyDictItem(dict, "index", index + 1);

  /* mutables */
  type_id1 = PConvStringToPyDictItem(dict, "type", atype);
  name_id1 = PConvStringToPyDictItem(dict, "name", at->name);
  resn_id1 = PConvStringToPyDictItem(dict, "resn", at->resn);
  flags_id1 = PConvIntToPyDictItem(dict, "flags", at->flags);
  resi_id1 = PConvStringToPyDictItem(dict, "resi", at->resi);
  resv_id1 = PConvIntToPyDictItem(dict, "resv", at->resv);      /* subordinate to resi */
  chain_id1 = PConvStringToPyDictItem(dict, "chain", at->chain);
  alt_id1 = PConvStringToPyDictItem(dict, "alt", at->alt);
  segi_id1 = PConvStringToPyDictItem(dict, "segi", at->segi);
  elem_id1 = PConvStringToPyDictItem(dict, "elem", at->elem);
  ss_id1 = PConvStringToPyDictItem(dict, "ss", at->ssType);
  numeric_type_id1 = PConvIntToPyDictItem(dict, "numeric_type", at->customType);
  q_id1 = PConvFloatToPyDictItem(dict, "q", at->q);
  b_id1 = PConvFloatToPyDictItem(dict, "b", at->b);
  vdw_id1 = PConvFloatToPyDictItem(dict, "vdw", at->vdw);
  elec_radius_id1 = PConvFloatToPyDictItem(dict, "elec_radius", at->elec_radius);
  partial_charge_id1 = PConvFloatToPyDictItem(dict, "partial_charge", at->partialCharge);
  formal_charge_id1 = PConvIntToPyDictItem(dict, "formal_charge", at->formalCharge);
  mmstereotype[0] = convertStereoToChar(at->mmstereo);
  mmstereotype[1] = 0;
  stereo_id1 = PConvStringToPyDictItem(dict, "stereo", mmstereotype);
  cartoon_id1 = PConvIntToPyDictItem(dict, "cartoon", at->cartoon);

  {
    char null_st[1] = "";
    char *st = null_st;

    if(at->textType)
      st = OVLexicon_FetchCString(G->Lexicon, at->textType);
    text_type_id1 = PConvStringToPyDictItem(dict, "text_type", st);

    st = null_st;
    if(at->label)
      st = OVLexicon_FetchCString(G->Lexicon, at->label);
    label_id1 = PConvStringToPyDictItem(dict, "label", st);
  }

  color_id1 = PConvIntToPyDictItem(dict, "color", at->color);
  ID_id1 = PConvIntToPyDictItem(dict, "ID", at->id);
  state_id1 = PConvIntToPyDictItem(dict, "state", at->discrete_state);
  rank_id1 = PConvIntToPyDictItem(dict, "rank", at->rank);

  PXDecRef(PyRun_String(expr, Py_single_input, space, dict));
  if(PyErr_Occurred()) {
    ErrMessage(G, "Alter", "Aborting on error. Assignment may be incomplete.");
    PyErr_Print();
    result = false;
  } else if(read_only) {
    result = true;
  }
  if(PyErr_Occurred()) {
    PyErr_Print();
    result = false;
  } else if(!read_only) {

    if(result) {
      /* get new object IDs */

      if(!(type_id2 = PyDict_GetItemString(dict, "type")))
        result = false;
      else if(!(name_id2 = PyDict_GetItemString(dict, "name")))
        result = false;
      else if(!(elem_id2 = PyDict_GetItemString(dict, "elem")))
        result = false;
      else if(!(resn_id2 = PyDict_GetItemString(dict, "resn")))
        result = false;
      else if(!(flags_id2 = PyDict_GetItemString(dict, "flags")))
        result = false;
      else if(!(resi_id2 = PyDict_GetItemString(dict, "resi")))
        result = false;
      else if(!(resv_id2 = PyDict_GetItemString(dict, "resv")))
        result = false;
      else if(!(segi_id2 = PyDict_GetItemString(dict, "segi")))
        result = false;
      else if(!(alt_id2 = PyDict_GetItemString(dict, "alt")))
        result = false;
      else if(!(chain_id2 = PyDict_GetItemString(dict, "chain")))
        result = false;
      else if(!(text_type_id2 = PyDict_GetItemString(dict, "text_type")))
        result = false;
      else if(!(ss_id2 = PyDict_GetItemString(dict, "ss")))
        result = false;
      else if(!(b_id2 = PyDict_GetItemString(dict, "b")))
        result = false;
      else if(!(q_id2 = PyDict_GetItemString(dict, "q")))
        result = false;
      else if(!(vdw_id2 = PyDict_GetItemString(dict, "vdw")))
        result = false;
      else if(!(elec_radius_id2 = PyDict_GetItemString(dict, "elec_radius")))
        result = false;
      else if(!(partial_charge_id2 = PyDict_GetItemString(dict, "partial_charge")))
        result = false;
      else if(!(formal_charge_id2 = PyDict_GetItemString(dict, "formal_charge")))
        result = false;
      else if(!(stereo_id2 = PyDict_GetItemString(dict, "stereo")))
        result = false;
      else if(!(cartoon_id2 = PyDict_GetItemString(dict, "cartoon")))
        result = false;
      else if(!(color_id2 = PyDict_GetItemString(dict, "color")))
        result = false;
      else if(!(label_id2 = PyDict_GetItemString(dict, "label")))
        result = false;
      if(!(numeric_type_id2 = PyDict_GetItemString(dict, "numeric_type")))
        result = false;
      if(!(ID_id2 = PyDict_GetItemString(dict, "ID")))
        result = false;
      if(!(state_id2 = PyDict_GetItemString(dict, "state")))
        result = false;
      if(!(rank_id2 = PyDict_GetItemString(dict, "rank")))
        result = false;

      if(PyErr_Occurred()) {
        PyErr_Print();
        result = false;
      }
    }
    if(result) {
      if(type_id1 != type_id2) {
        if(!PConvPyObjectToStrMaxLen(type_id2, atype, 6))
          result = false;
        else
          at->hetatm = ((atype[0] == 'h') || (atype[0] == 'H'));
      }
      if(name_id1 != name_id2) {
        if(!PConvPyObjectToStrMaxLen(name_id2, name, sizeof(AtomName) - 1))
          result = false;
        else
          strcpy(at->name, name);
      }
      if(elem_id1 != elem_id2) {
        if(!PConvPyObjectToStrMaxLen(elem_id2, elem, sizeof(ElemName) - 1))
          result = false;
        else {
          strcpy(at->elem, elem);
          AtomInfoAssignParameters(G, at);
        }
      }
      if(resn_id1 != resn_id2) {
        if(!PConvPyObjectToStrMaxLen(resn_id2, resn, sizeof(ResName) - 1))
          result = false;
        else
          strcpy(at->resn, resn);
      }
      if(resi_id1 != resi_id2) {
        if(!PConvPyObjectToStrMaxLen(resi_id2, resi, sizeof(ResIdent) - 1))
          result = false;
        else {
          if(strcmp(at->resi, resi) != 0)
            at->resv = AtomResvFromResi(resi);
          strcpy(at->resi, resi);
        }
      } else if(resv_id1 != resv_id2) {
        if(!PConvPyObjectToInt(resv_id2, &resv))
          result = false;
        else {
          sprintf(buf, "%d", resv);
          buf[sizeof(ResIdent) - 1] = 0;
          strcpy(at->resi, buf);
        }

      }
      if(segi_id1 != segi_id2) {
        if(!PConvPyObjectToStrMaxLen(segi_id2, segi, sizeof(SegIdent) - 1))
          result = false;
        else
          strcpy(at->segi, segi);

      }
      if(chain_id1 != chain_id2) {
        if(!PConvPyObjectToStrMaxLen(chain_id2, chain, sizeof(Chain) - 1))
          result = false;
        else
          strcpy(at->chain, chain);
      }
      if(alt_id1 != alt_id2) {
        if(!PConvPyObjectToStrMaxLen(alt_id2, alt, sizeof(Chain) - 1))
          result = false;
        else
          strcpy(at->alt, alt);
      }
      if(text_type_id1 != text_type_id2) {

        OrthoLineType temp;
        if(at->textType) {
          OVLexicon_DecRef(G->Lexicon, at->textType);
        }
        at->textType = 0;

        if(PConvPyObjectToStrMaxLen(text_type_id2, temp, sizeof(OrthoLineType) - 1)) {
          if(temp[0]) {
            OVreturn_word result = OVLexicon_GetFromCString(G->Lexicon, temp);
            if(OVreturn_IS_OK(result)) {
              at->textType = result.word;
            }
          }
        }
      }
      if(ss_id1 != ss_id2) {
        if(!PConvPyObjectToStrMaxLen(ss_id2, ssType, sizeof(SSType) - 1))
          result = false;
        else {
          strcpy(at->ssType, ssType);
          at->ssType[0] = toupper(at->ssType[0]);
        }
      }
      if(b_id1 != b_id2) {
        if(!PConvPyObjectToFloat(b_id2, &b))
          result = false;
        else
          at->b = b;
      }
      if(q_id1 != q_id2) {
        if(!PConvPyObjectToFloat(q_id2, &q))
          result = false;
        else
          at->q = q;
      }
      if(vdw_id1 != vdw_id2) {
        if(!PConvPyObjectToFloat(vdw_id2, &vdw))
          result = false;
        else
          at->vdw = vdw;

      }
      if(elec_radius_id1 != elec_radius_id2) {
        if(!PConvPyObjectToFloat(elec_radius_id2, &elec_radius))
          result = false;
        else
          at->elec_radius = elec_radius;
      }
      if(partial_charge_id1 != partial_charge_id2) {
        if(!PConvPyObjectToFloat(partial_charge_id2, &partialCharge))
          result = false;
        else
          at->partialCharge = partialCharge;

      }
      if(formal_charge_id1 != formal_charge_id2) {
        if(!PConvPyObjectToInt(formal_charge_id2, &formalCharge))
          result = false;
        else {
          at->formalCharge = formalCharge;
          at->chemFlag = false; /* invalidate chemistry info for this atom */
        }

      }
      if(stereo_id1 != stereo_id2) {
        if(!PConvPyObjectToStrMaxLen(stereo_id2, stereo, 2))
          result = false;
        else
          at->mmstereo = convertStereoToChar(stereo[0]);
      }
      if(cartoon_id1 != cartoon_id2) {
        if(!PConvPyObjectToInt(cartoon_id2, &cartoon))
          result = false;
        else
          at->cartoon = cartoon;
      }
      if(color_id1 != color_id2) {
        if(!PConvPyObjectToInt(color_id2, &color))
          result = false;
        else
          at->color = color;
      }
      if(state_id1 != state_id2) {
        if(!PConvPyObjectToInt(state_id2, &state))
          result = false;
        else
          at->discrete_state = state;
      }
      if(label_id1 != label_id2) {
        OrthoLineType temp;
        if(at->label) {
          OVLexicon_DecRef(G->Lexicon, at->label);
        }
        at->label = 0;

        if(PConvPyObjectToStrMaxLen(label_id2, temp, sizeof(OrthoLineType) - 1)) {
          if(temp[0]) {
            OVreturn_word result = OVLexicon_GetFromCString(G->Lexicon, temp);
            if(OVreturn_IS_OK(result)) {
              at->label = result.word;
            }
          }
        }
      }
      if(flags_id1 != flags_id2) {
        if(!PConvPyObjectToInt(flags_id2, &flags))
          result = false;
        else
          at->flags = flags;
      }

      if(numeric_type_id1 != numeric_type_id2) {
        if(!PConvPyObjectToInt(numeric_type_id2, &numericType))
          result = false;
        else
          at->customType = numericType;
      }
      if(ID_id1 != ID_id2) {
        if(!PConvPyObjectToInt(ID_id2, &id))
          result = false;
        else
          at->id = id;
      }
      if(rank_id1 != rank_id2) {
        if(!PConvPyObjectToInt(rank_id2, &rank))
          result = false;
        else
          at->rank = rank;
      }

      if(PyErr_Occurred()) {
        PyErr_Print();
        result = false;
      }
    }
    if(!result) {
      ErrMessage(G, "Alter", "Aborting on error. Assignment may be incomplete.");
    }
  }
  Py_DECREF(dict);
  return (result);
}

int PLabelAtom(PyMOLGlobals * G, AtomInfoType * at, char *model, char *expr, int index)
{
  PyObject *dict;
  PyObject *P_inst_dict = G->P_inst->dict;
  int result;
  OrthoLineType label;
  char atype[7], mmstereotype[2];
  OrthoLineType buffer;
  if(at->hetatm)
    strcpy(atype, "HETATM");
  else
    strcpy(atype, "ATOM");
  PBlock(G);
  dict = PyDict_New();

  PConvStringToPyDictItem(dict, "model", model);
  PConvIntToPyDictItem(dict, "index", index + 1);
  PConvStringToPyDictItem(dict, "type", atype);
  PConvStringToPyDictItem(dict, "name", at->name);
  PConvStringToPyDictItem(dict, "resn", at->resn);
  PConvStringToPyDictItem(dict, "resi", at->resi);
  PConvIntToPyDictItem(dict, "resv", at->resv);
  PConvStringToPyDictItem(dict, "chain", at->chain);
  PConvStringToPyDictItem(dict, "alt", at->alt);
  PConvStringToPyDictItem(dict, "segi", at->segi);
  PConvStringToPyDictItem(dict, "ss", at->ssType);
  PConvFloatToPyDictItem(dict, "vdw", at->vdw);
  PConvFloatToPyDictItem(dict, "elec_radius", at->elec_radius);
  {
    char null_st[1] = "";
    char *st = null_st;

    if(at->textType)
      st = OVLexicon_FetchCString(G->Lexicon, at->textType);
    PConvStringToPyDictItem(dict, "text_type", st);

    st = null_st;
    if(at->label)
      st = OVLexicon_FetchCString(G->Lexicon, at->label);
    PConvStringToPyDictItem(dict, "label", st);
  }
  PConvStringToPyDictItem(dict, "elem", at->elem);
  PConvIntToPyDictItem(dict, "geom", at->geom);
  PConvIntToPyDictItem(dict, "valence", at->valence);
  PConvIntToPyDictItem(dict, "rank", at->rank);
  if(at->flags) {
    sprintf(buffer, "%X", at->flags);
    PConvStringToPyDictItem(dict, "flags", buffer);
  } else {
    PConvStringToPyDictItem(dict, "flags", "0");
  }
  PConvFloatToPyDictItem(dict, "q", at->q);
  PConvFloatToPyDictItem(dict, "b", at->b);
  if(at->customType != cAtomInfoNoType)
    PConvIntToPyDictItem(dict, "numeric_type", at->customType);
  else
    PConvStringToPyDictItem(dict, "numeric_type", "?");
  PConvFloatToPyDictItem(dict, "partial_charge", at->partialCharge);
  PConvIntToPyDictItem(dict, "formal_charge", at->formalCharge);
  mmstereotype[0] = convertStereoToChar(at->mmstereo);
  mmstereotype[1] = 0;
  PConvStringToPyDictItem(dict, "stereo", mmstereotype);
  PConvIntToPyDictItem(dict, "color", at->color);
  PConvIntToPyDictItem(dict, "cartoon", at->cartoon);
  PConvIntToPyDictItem(dict, "ID", at->id);
  PXDecRef(PyRun_String(expr, Py_single_input, P_inst_dict, dict));
  if(PyErr_Occurred()) {
    PyErr_Print();
    result = false;
  } else {
    result = true;
    if(!PConvPyObjectToStrMaxLen(PyDict_GetItemString(dict, "label"),
                                 label, sizeof(OrthoLineType) - 1))
      result = false;
    if(PyErr_Occurred()) {
      PyErr_Print();
      result = false;
    }
    if(result) {
      if(at->label) {
        OVLexicon_DecRef(G->Lexicon, at->label);
      }
      at->label = 0;

      if(label[0]) {
        OVreturn_word ret = OVLexicon_GetFromCString(G->Lexicon, label);
        if(OVreturn_IS_OK(ret)) {
          /*printf("alloc'd %d [%s]\n",OVLexicon_GetNActive(G->Lexicon),label); */
          at->label = ret.word;
        }
      }
    } else {
      ErrMessage(G, "Label", "Aborting on error. Labels may be incomplete.");
    }
  }
  Py_DECREF(dict);
  PUnblock(G);
  return (result);
}

void PUnlockAPIAsGlut(PyMOLGlobals * G)
{                               /* must call with unblocked interpreter */
  PRINTFD(G, FB_Threads)
    " PUnlockAPIAsGlut-DEBUG: entered as thread 0x%x\n", PyThread_get_thread_ident()
    ENDFD;
  PBlock(G);
  PXDecRef(PyObject_CallFunction(G->P_inst->unlock, "iO", 0, G->P_inst->cmd));  /* NOTE this may flush the command buffer! */
  PLockStatus(G);
  PyMOL_PopValidContext(G->PyMOL);
  PUnlockStatus(G);
  PUnlockGLUT(G);
  PUnblock(G);
}

void PUnlockAPIAsGlutNoFlush(PyMOLGlobals * G)
{                               /* must call with unblocked interpreter */
  PRINTFD(G, FB_Threads)
    " PUnlockAPIAsGlut-DEBUG: entered as thread 0x%x\n", PyThread_get_thread_ident()
    ENDFD;
  PBlock(G);
  PXDecRef(PyObject_CallFunction(G->P_inst->unlock, "iO", -1, G->P_inst->cmd)); /* prevents flushing of the buffer */
  PLockStatus(G);
  PyMOL_PopValidContext(G->PyMOL);
  PUnlockStatus(G);
  PUnlockGLUT(G);
  PUnblock(G);
}

static int get_api_lock(PyMOLGlobals * G, int block_if_busy)
{
  int result = true;

  if(block_if_busy) {

    PXDecRef(PyObject_CallFunction(G->P_inst->lock, "O", G->P_inst->cmd));

  } else {                      /* not blocking if PyMOL is busy */

    PyObject *got_lock =
      PyObject_CallFunction(G->P_inst->lock_attempt, "O", G->P_inst->cmd);

    if(got_lock) {
      if(!PyInt_AsLong(got_lock)) {
        if(!G) {                /* impossible (unless stack trashed?) */
          result = false;
        } else {
          PLockStatus(G);
          if(PyMOL_GetBusy(G->PyMOL, false))
            result = false;
          PUnlockStatus(G);
          if(!G) {              /* impossible (unless stack trashed?) */
            result = false;
          } else {
            if(result) {        /* didn't get lock, but not busy, so block and wait for lock */
              PXDecRef(PyObject_CallFunction(G->P_inst->lock, "O", G->P_inst->cmd));
            }
          }
        }
      }
      Py_DECREF(got_lock);
    }
  }
  return result;
}

int PLockAPIAsGlut(PyMOLGlobals * G, int block_if_busy)
{
  PRINTFD(G, FB_Threads)
    "*PLockAPIAsGlut-DEBUG: entered as thread 0x%x\n", PyThread_get_thread_ident()
    ENDFD;

  PBlock(G);

  PLockGLUT(G);

  PLockStatus(G);
  PyMOL_PushValidContext(G->PyMOL);
  PUnlockStatus(G);

  PRINTFD(G, FB_Threads)
    "#PLockAPIAsGlut-DEBUG: acquiring lock as thread 0x%x\n", PyThread_get_thread_ident()
    ENDFD;

  if(!get_api_lock(G, block_if_busy)) {
    PLockStatus(G);
    PyMOL_PopValidContext(G->PyMOL);
    PUnlockStatus(G);
    PUnlockGLUT(G);
    PUnblock(G);
    return false;               /* busy -- so allow main to update busy status display (if any) */
  }

  while(G->P_inst->glut_thread_keep_out) {
    /* IMPORTANT: keeps the glut thread out of an API operation... */
    /* NOTE: the keep_out variable can only be changed or read by the thread
       holding the API lock, therefore it is safe even through increment
       isn't atomic. */
    PRINTFD(G, FB_Threads)
      "-PLockAPIAsGlut-DEBUG: glut_thread_keep_out 0x%x\n", PyThread_get_thread_ident()
      ENDFD;

    PXDecRef(PyObject_CallFunction(G->P_inst->unlock, "iO", -1, G->P_inst->cmd));       /* prevent buffer flushing */
#ifndef WIN32
    {
      struct timeval tv;

      PUnblock(G);
      tv.tv_sec = 0;
      tv.tv_usec = 50000;
      select(0, NULL, NULL, NULL, &tv);
      PBlock(G);
    }
#else
    /* BEGIN PROPRIETARY CODE SEGMENT (see disclaimer in "os_proprietary.h") */
    PXDecRef(PyObject_CallFunction(P_sleep, "f", 0.050));
    /* END PROPRIETARY CODE SEGMENT */
#endif

    if(!get_api_lock(G, block_if_busy)) {
      /* return false-- allow main to update busy status display (if any) */
      PLockStatus(G);
      PyMOL_PopValidContext(G->PyMOL);
      PUnlockStatus(G);
      PUnlockGLUT(G);
      PUnblock(G);
      return false;
    }
  }

  PUnblock(G);                  /* API is now locked, so we can free up Python... */

  PRINTFD(G, FB_Threads)
    "=PLockAPIAsGlut-DEBUG: acquired\n" ENDFD;
  return true;
}


/* THESE CALLS ARE REQUIRED FOR MONOLITHIC COMPILATION TO SUCCEED UNDER WINDOWS. */
#ifndef _PYMOL_EMBEDDED

/* 
 *  void        initExtensionClass(void);
 *   void        initsglite(void);
 */
void init_champ(void);
void init_opengl(void);
void init_opengl_num(void);
void init_glu(void);
void init_glu_num(void);
void init_glut(void);
void initopenglutil(void);
void initopenglutil_num(void);
#endif

#ifdef _PYMOL_MONOLITHIC
#ifndef _PYMOL_EMBEDDED

/* BEGIN PROPRIETARY CODE SEGMENT (see disclaimer in "os_proprietary.h") */
#ifdef WIN32
#ifdef _PYMOL_NUMPY_INIT
void init_numpy();
void initmultiarray();
void initarrayfns();
void initlapack_lite();
void initumath();
void initranlib();
#endif
#endif

/* END PROPRIETARY CODE SEGMENT */
#endif
#endif

#ifdef _PYMOL_MONOLITHIC
#ifndef _PYMOL_EMBEDDED

/*
 * void        initExtensionClass(void);
 * void        initsglite(void);
 */
void init_champ(void);
#ifdef _PYMOL_PYOMM
void init_pyomm(void);
#endif
void init_opengl(void);
void init_opengl_num(void);
void init_glu(void);
void init_glu_num(void);
void init_glut(void);
void initopenglutil(void);
void initopenglutil_num(void);
#endif
#endif


/* BEGIN PROPRIETARY CODE SEGMENT (see disclaimer in "os_proprietary.h") */
#ifdef WIN32
static int IsSecurityRequired()
{
  DWORD WindowsVersion = GetVersion();
  DWORD WindowsMajorVersion = (DWORD) (LOBYTE(LOWORD(WindowsVersion)));
  DWORD WindowsMinorVersion = (DWORD) (HIBYTE(LOWORD(WindowsVersion)));

  if(WindowsVersion >= 0x80000000)
    return FALSE;

  return TRUE;
}
#endif

/* END PROPRIETARY CODE SEGMENT */

void PSetupEmbedded(PyMOLGlobals * G, int argc, char **argv)
{
  /* This routine is called if we are running with an embedded Python interpreter */
  PyObject *args, *pymol;

  /* BEGIN PROPRIETARY CODE SEGMENT (see disclaimer in "os_proprietary.h") */
#ifdef WIN32

  /* Windows PyMOL now ships with Python 2.5 for both
     32 and 64 bit */

#ifndef EMBEDDED_PYTHONHOME
#define EMBEDDED_PYTHONHOME "\\py27"
#endif

  {                             /* Automatically hide the window if this process was started as a
                                   vanilla console application (by double-clicking).
                                   Conveniently, this doesn't hide the window when launched from a
                                   command window. */
    HWND hwndFound;
    if(hwndFound = FindWindow(NULL, argv[0])) {
      ShowWindow(hwndFound, SW_HIDE);
    }
  }

  {                             /* if PYMOL_PATH and/or PYTHONHOME isn't in the environment coming
                                   in, then the user may simply have clicked PyMOL.exe, in which
                                   case we need to consult the registry regarding the location of
                                   the install */

    static char line1[8092];
    static char line2[8092];

    {                           /* If no PYMOL_PATH specific, but we were launched with an
                                 * absolute path, then try using that path first.  With embedded
                                 * builds, the .EXE should always be located at the root of
                                 * $PYMOL_PATH */

      char *pymol_path = getenv("PYMOL_PATH");
      if((!pymol_path) && (argc > 0) && argv[0][0] && (argv[0][1] == ':')
         && (argv[0][2] == '\\')) {

        char *p;
        strcpy(line1, "PYMOL_PATH=");
        strcat(line1, argv[0]);
        p = line1 + strlen(line1);
        while(*p != '\\') {
          *p = 0;
          p--;
        }
        *p = 0;
        putenv(line1);
      }
    }

    {
      OrthoLineType path_buffer;
      HKEY phkResult;
      int lpcbData;
      int lpType = REG_SZ;
      int r1, r2;
      char *pymol_path;
      char *pythonhome;
      int pythonhome_set = false;
      int restart_flag = false;

      pymol_path = getenv("PYMOL_PATH");
      pythonhome = getenv("PYTHONHOME");
      if((!pymol_path) || (!pythonhome)) {
        lpcbData = sizeof(OrthoLineType) - 1;
        r1 = RegOpenKeyEx(HKEY_CLASSES_ROOT,
#ifdef PYMOL_EVAL
			"Software\\PyMOL\\PyMOL Eval\\PYMOL_PATH",
#else
			"Software\\PyMOL\\PyMOL\\PYMOL_PATH",
#endif
                          0, KEY_EXECUTE, &phkResult);
        if(r1 != ERROR_SUCCESS) {
          r1 = RegOpenKeyEx(HKEY_CURRENT_USER,
#ifdef PYMOL_EVAL
			"Software\\PyMOL\\PyMOL Eval\\PYMOL_PATH",
#else
			"Software\\PyMOL\\PyMOL\\PYMOL_PATH",
#endif
                            0, KEY_EXECUTE, &phkResult);
        }
        if(r1 == ERROR_SUCCESS) {
          r2 = RegQueryValueEx(phkResult, "", NULL, &lpType, path_buffer, &lpcbData);
          if(r2 == ERROR_SUCCESS) {
            /* use environment variable PYMOL_PATH first, registry entry
               second */
            if(!pymol_path) {
              strcpy(line1, "PYMOL_PATH=");
              strcat(line1, path_buffer);
              _putenv(line1);
              if(!pythonhome) { /* only set PYTHONHOME if already
                                   setting new PYMOL_PATH */
                pythonhome_set = true;
                strcpy(line2, "PYTHONHOME=");
                strcat(line2, path_buffer);
                strcat(line2, EMBEDDED_PYTHONHOME);
                restart_flag = true;
                _putenv(line2);
              }
            }
          }
          RegCloseKey(phkResult);
        }
        /* this allows us to just specify PYMOL_PATH with no registry entries */
        if((!pythonhome_set) && (pymol_path) && (!pythonhome)) {
          strcpy(line2, "PYTHONHOME=");
          strcat(line2, pymol_path);
          strcat(line2, EMBEDDED_PYTHONHOME);
          _putenv(line2);
          restart_flag = true;
        }
      }
      if(restart_flag && getenv("PYMOL_PATH") && getenv("PYTHONHOME")) {

        /* now that we have the environment defined, restart the process
         * so that Python can use the new environment.  If we don't do
         * this, then Python won't see the new environment vars. Why not? */

        /* note that we use CreateProcesss to launch the console
         * application instead of exec or spawn in order to hide the
         * console window. Otherwise a console window might appear, and
         * that would suck. */

        char command[8092];
        static char cmd_line[8092];
        char *p, *q;
        int a;

        /* copy arguments, installing quotes around them */

        sprintf(command, "%s\\pymol.exe", getenv("PYMOL_PATH"));
        p = cmd_line;

        sprintf(p, "\"%s\"", command);
        p += strlen(p);
        *(p++) = ' ';
        *p = 0;

        for(a = 1; a <= argc; a++) {
          q = argv[a];
          if(q) {
            if(*q != '"') {     /* add quotes if not present */
              *(p++) = '"';
              while(*q) {
                *(p++) = *(q++);
              }
              *(p++) = '"';
            } else {
              while(*q) {
                *(p++) = *(q++);
              }
            }
            *(p++) = 32;
            *p = 0;
          }
        }

        {
          LPSECURITY_ATTRIBUTES lpSA = NULL;
          PSECURITY_DESCRIPTOR lpSD = NULL;
          STARTUPINFO si;
          PROCESS_INFORMATION pi;
          HANDLE hProcess = GetCurrentProcess();

          ZeroMemory(&si, sizeof(STARTUPINFO));
          si.cb = sizeof(STARTUPINFO);
          si.dwFlags = STARTF_USESHOWWINDOW;
          si.wShowWindow = SW_HIDE;

          if(IsSecurityRequired()) {
            lpSD = GlobalAlloc(GPTR, SECURITY_DESCRIPTOR_MIN_LENGTH);
            InitializeSecurityDescriptor(lpSD, SECURITY_DESCRIPTOR_REVISION);
            SetSecurityDescriptorDacl(lpSD, -1, 0, 0);

            lpSA = GlobalAlloc(GPTR, sizeof(SECURITY_ATTRIBUTES));
            lpSA->nLength = sizeof(SECURITY_ATTRIBUTES);
            lpSA->lpSecurityDescriptor = lpSD;
            lpSA->bInheritHandle = TRUE;
          }

          if(CreateProcessA(NULL, (LPTSTR) cmd_line, lpSA, NULL, TRUE,
                            0, NULL, NULL, &si, &pi)) {

            WaitForSingleObject(pi.hProcess, INFINITE);
          } else {
            printf("ERROR: Unable to restart PyMOL process with new environment:\n");
            system("set");      /* dump the environment. */
            printf("CreateProcess failed, code %d: %s\n", GetLastError(), cmd_line);
            printf("PyMOL will now terminate.\n");
          }

          if(lpSA != NULL)
            GlobalFree(lpSA);
          if(lpSD != NULL)
            GlobalFree(lpSD);
          _exit(0);
        }
      }
    }
  }
#endif
  /* END PROPRIETARY CODE SEGMENT */

  /* compatibility for old compile-time defines */

#ifdef _PYMOL_SETUP_PY21
#ifndef _PYMOL_SETUP_PY_EXT
#define _PYMOL_SETUP_PY_EXT
#endif
#endif
#ifdef _PYMOL_SETUP_PY22
#ifndef _PYMOL_SETUP_PY_EXT
#define _PYMOL_SETUP_PY_EXT
#endif
#endif
#ifdef _PYMOL_SETUP_PY23
#ifndef _PYMOL_SETUP_PY_EXT
#define _PYMOL_SETUP_PY_EXT
#endif
#endif
#ifdef _PYMOL_SETUP_PY24
#ifndef _PYMOL_SETUP_PY_EXT
#define _PYMOL_SETUP_PY_EXT
#endif
#endif
#ifdef _PYMOL_SETUP_PY25
#ifndef _PYMOL_SETUP_PY_EXT
#define _PYMOL_SETUP_PY_EXT
#endif
#endif
#ifdef _PYMOL_SETUP_PY26
#ifndef _PYMOL_SETUP_PY_EXT
#define _PYMOL_SETUP_PY_EXT
#endif
#endif

  /* should we set up PYTHONHOME in the ext directory? */

#ifdef _PYMOL_SETUP_PY_EXT
  {
    static char line1[8092];
    static char line2[8092];
    if(!getenv("PYMOL_PATH")) { /* if PYMOL_PATH isn't defined... */

      /* was our startup path absolute? */

      if((argc > 0) && (argv[0][0] == '/')) {
        /* PYMOL was started with an absolute path, so try using that... */
        char *p;
        strcpy(line1, "PYMOL_PATH=");
        strcat(line1, argv[0]);
        p = line1 + strlen(line1);
        while(*p != '/') {
          *p = 0;
          p--;
        }
        *p = 0;
        putenv(line1);
      } else if((argc > 0) && getenv("PWD")
                && ((argv[0][0] == '.') || (strstr(argv[0], "/")))) {
        /* was the path relative? */
        char *p;
        strcpy(line1, "PYMOL_PATH=");
        strcat(line1, getenv("PWD"));
        strcat(line1, "/");
        strcat(line1, argv[0]);
        p = line1 + strlen(line1);
        while(*p != '/') {
          *p = 0;
          p--;
        }
        *p = 0;
        putenv(line1);
      } else {                  /* otherwise, just try using the current working directory */
        if(getenv("PWD")) {
          strcpy(line1, "PYMOL_PATH=");
          strcat(line1, getenv("PWD"));
          putenv(line1);
        }
      }
    }

    /* now set PYTHONHOME so that we use the right binary libraries for
       this executable */

    if(getenv("PYMOL_PATH")) {
      strcpy(line2, "PYTHONHOME=");
      strcat(line2, getenv("PYMOL_PATH"));
      strcat(line2, "/ext");
      putenv(line2);
    }
  }
#endif

#ifndef _PYMOL_EMBEDDED
  Py_Initialize();
  PyEval_InitThreads();
  PyUnicode_SetDefaultEncoding("utf-8");        /* is this safe & legal? */
#endif

  init_cmd();

#ifdef _PYMOL_MONOLITHIC
#ifndef _PYMOL_EMBEDDED
  /*
   * initExtensionClass();
   * initsglite();
   */
  /* initialize champ */
  init_champ();

#ifdef _PYMOL_PYOMM
  init_pyomm();
#endif

  /* BEGIN PROPRIETARY CODE SEGMENT (see disclaimer in "os_proprietary.h") */
#ifdef WIN32
#ifdef _PYMOL_NUMPY_INIT
  /* initialize numeric python */
  init_numpy();
  initmultiarray();
  initarrayfns();
  initlapack_lite();
  initumath();
  initranlib();
#endif
#endif

  /* END PROPRIETARY CODE SEGMENT */
#endif
#endif

  PyRun_SimpleString("import os\n");
  PyRun_SimpleString("import sys\n");
  /* BEGIN PROPRIETARY CODE SEGMENT (see disclaimer in "os_proprietary.h") */
#ifdef WIN32
  {
    /* getenv('PYMOL_PATH') and os.environ['PYMOL_PATH'] aren't
       automatically synchronized on Windows, so here we do the job
       manually... */

    char *pymol_path = getenv("PYMOL_PATH");
    if(pymol_path) {
      PyObject *os = PyImport_AddModule("os");  /* borrowed ref */
      char *buffer = Alloc(char, strlen(pymol_path) + 100);
      if(os && buffer) {
        PyObject *envir = PyObject_GetAttrString(os, "environ");
        if(envir) {
          if(!PTruthCallStr1s(envir, "has_key", "PYMOL_PATH")) {
            sprintf(buffer, "os.environ['PYMOL_PATH']=r'''%s'''\n", pymol_path);
            PyRun_SimpleString(buffer);
          }
        }
        PXDecRef(envir);
      }
      FreeP(buffer);
    }
  }
  /* ultimate fallback -- try using the current working directory */
  PyRun_SimpleString
    ("if not os.environ.has_key('PYMOL_PATH'): os.environ['PYMOL_PATH']=os.getcwd()\n");
#endif
  /* END PROPRIETARY CODE SEGMENT */

#ifdef _PYMOL_SETUP_TCLTK83
  /* used by semistatic pymol */
  PyRun_SimpleString
    ("if os.path.exists(os.environ['PYMOL_PATH']+'/ext/lib/tcl8.3'): os.environ['TCL_LIBRARY']=os.environ['PYMOL_PATH']+'/ext/lib/tcl8.3'\n");
  PyRun_SimpleString
    ("if os.path.exists(os.environ['PYMOL_PATH']+'/ext/lib/tk8.3'): os.environ['TK_LIBRARY']=os.environ['PYMOL_PATH']+'/ext/lib/tk8.3'\n");
#endif

#ifdef _PYMOL_SETUP_TCLTK84
  /* used by semistatic pymol */
  PyRun_SimpleString
    ("if os.path.exists(os.environ['PYMOL_PATH']+'/ext/lib/tcl8.4'): os.environ['TCL_LIBRARY']=os.environ['PYMOL_PATH']+'/ext/lib/tcl8.4'\n");
  PyRun_SimpleString
    ("if os.path.exists(os.environ['PYMOL_PATH']+'/ext/lib/tk8.4'): os.environ['TK_LIBRARY']=os.environ['PYMOL_PATH']+'/ext/lib/tk8.4'\n");
#endif

#ifdef _PYMOL_SETUP_TCLTK85
  /* used by semistatic pymol */
  PyRun_SimpleString
    ("if os.path.exists(os.environ['PYMOL_PATH']+'/ext/lib/tcl8.5'): os.environ['TCL_LIBRARY']=os.environ['PYMOL_PATH']+'/ext/lib/tcl8.5'\n");
  PyRun_SimpleString
    ("if os.path.exists(os.environ['PYMOL_PATH']+'/ext/lib/tk8.5'): os.environ['TK_LIBRARY']=os.environ['PYMOL_PATH']+'/ext/lib/tk8.5'\n");
#endif

  /* BEGIN PROPRIETARY CODE SEGMENT (see disclaimer in "os_proprietary.h") */
#ifdef WIN32
  PyRun_SimpleString
    ("if (os.environ['PYMOL_PATH']+'/modules') not in sys.path: sys.path.insert(0,os.environ['PYMOL_PATH']+'/modules')\n");
#endif
  /* END PROPRIETARY CODE SEGMENT */

  P_main = PyImport_AddModule("__main__");
  if(!P_main)
    ErrFatal(G, "PyMOL", "can't find '__main__'");

  /* inform PyMOL's other half that we're launching embedded-style */
  PyObject_SetAttrString(P_main, "pymol_launch", PyInt_FromLong(4));

  args = PConvStringListToPyList(argc, argv);   /* prepare our argument list */
  if(!args)
    ErrFatal(G, "PyMOL", "can't process arguments.");

  /* copy arguments to __main__.pymol_argv */
  PyObject_SetAttrString(P_main, "pymol_argv", args);
  PyRun_SimpleString
    ("import __main__\nif not hasattr(sys,'argv'): sys.argv=__main__.pymol_argv");

  PyRun_SimpleString("if (os.environ['PYMOL_PATH']+'/modules') not in sys.path: sys.path.insert(0,os.environ['PYMOL_PATH']+'/modules')\n");     /* needed for semistatic pymol */

  PyRun_SimpleString("import pymol");   /* create the global PyMOL namespace */

  pymol = PyImport_AddModule("pymol");  /* get it */
  if(!pymol)
    ErrFatal(G, "PyMOL", "can't find module 'pymol'");

}

void PConvertOptions(CPyMOLOptions * rec, PyObject * options)
{
  char *load_str;

  rec->pmgui = !PyInt_AsLong(PyObject_GetAttrString(options, "no_gui"));
  rec->internal_gui = PyInt_AsLong(PyObject_GetAttrString(options, "internal_gui"));
  rec->internal_feedback =
    PyInt_AsLong(PyObject_GetAttrString(options, "internal_feedback"));
  rec->show_splash = PyInt_AsLong(PyObject_GetAttrString(options, "show_splash"));
  rec->security = PyInt_AsLong(PyObject_GetAttrString(options, "security"));
  rec->game_mode = PyInt_AsLong(PyObject_GetAttrString(options, "game_mode"));
  rec->force_stereo = PyInt_AsLong(PyObject_GetAttrString(options, "force_stereo"));
  rec->winX = PyInt_AsLong(PyObject_GetAttrString(options, "win_x"));
  rec->winY = PyInt_AsLong(PyObject_GetAttrString(options, "win_y"));
  rec->winPX = PyInt_AsLong(PyObject_GetAttrString(options, "win_px"));
  rec->winPY = PyInt_AsLong(PyObject_GetAttrString(options, "win_py"));
  rec->blue_line = PyInt_AsLong(PyObject_GetAttrString(options, "blue_line"));
  rec->external_gui = PyInt_AsLong(PyObject_GetAttrString(options, "external_gui"));
  rec->siginthand = PyInt_AsLong(PyObject_GetAttrString(options, "sigint_handler"));
  rec->reuse_helper = PyInt_AsLong(PyObject_GetAttrString(options, "reuse_helper"));
  rec->auto_reinitialize =
    PyInt_AsLong(PyObject_GetAttrString(options, "auto_reinitialize"));
  rec->keep_thread_alive =
    PyInt_AsLong(PyObject_GetAttrString(options, "keep_thread_alive"));
  rec->quiet = PyInt_AsLong(PyObject_GetAttrString(options, "quiet"));
#ifdef _PYMOL_IP_EXTRAS
  rec->incentive_product = true;
  PyObject_SetAttrString(options, "incentive_product", PyInt_FromLong(1));
#else
  rec->incentive_product =
    PyInt_AsLong(PyObject_GetAttrString(options, "incentive_product"));
#endif
  rec->multisample = PyInt_AsLong(PyObject_GetAttrString(options, "multisample"));
  rec->window_visible = PyInt_AsLong(PyObject_GetAttrString(options, "window_visible"));
  rec->read_stdin = PyInt_AsLong(PyObject_GetAttrString(options, "read_stdin"));
  rec->presentation = PyInt_AsLong(PyObject_GetAttrString(options, "presentation"));
  rec->defer_builds_mode =
    PyInt_AsLong(PyObject_GetAttrString(options, "defer_builds_mode"));
  rec->full_screen = PyInt_AsLong(PyObject_GetAttrString(options, "full_screen"));
  load_str = PyString_AsString(PyObject_GetAttrString(options, "after_load_script"));
  rec->sphere_mode = PyInt_AsLong(PyObject_GetAttrString(options, "sphere_mode"));
  rec->stereo_capable = PyInt_AsLong(PyObject_GetAttrString(options, "stereo_capable"));
  rec->stereo_mode = PyInt_AsLong(PyObject_GetAttrString(options, "stereo_mode"));
  rec->zoom_mode = PyInt_AsLong(PyObject_GetAttrString(options, "zoom_mode"));
  rec->no_quit = PyInt_AsLong(PyObject_GetAttrString(options, "no_quit"));

  if(load_str) {
    if(load_str[0]) {
      UtilNCopy(rec->after_load_script, load_str, PYMOL_MAX_OPT_STR);
    }
  }
  if(PyErr_Occurred()) {
    PyErr_Print();
  }
}

void PGetOptions(CPyMOLOptions * rec)
{
  PyObject *pymol, *invocation, *options;

  pymol = PyImport_AddModule("pymol");  /* borrowed reference!!! */
  if(!pymol) {
    fprintf(stderr, "PyMOL-ERROR: can't find module 'pymol'");
    exit(EXIT_FAILURE);
  }

  invocation = PyObject_GetAttrString(pymol, "invocation");     /* get a handle to the invocation module */
  if(!invocation) {
    fprintf(stderr, "PyMOL-ERROR: can't find module 'invocation'");
    exit(EXIT_FAILURE);
  }

  options = PyObject_GetAttrString(invocation, "options");
  if(!options) {
    fprintf(stderr, "PyMOL-ERROR: can't get 'invocation.options'.");
    exit(EXIT_FAILURE);
  }

  PConvertOptions(rec, options);
  Py_XDECREF(invocation);
  Py_XDECREF(options);
}

void PRunStringModule(PyMOLGlobals * G, char *str)
{                               /* runs a string in the namespace of the pymol global module */
  PXDecRef(PyObject_CallFunction(G->P_inst->exec, "Os", P_pymol, str));
}

void PRunStringInstance(PyMOLGlobals * G, char *str)
{                               /* runs a string in the namespace of the pymol instance */
  PXDecRef(PyObject_CallFunction(G->P_inst->exec, "Os", G->P_inst->obj, str));
}

void PInit(PyMOLGlobals * G, int global_instance)
{
  PyObject *sys, *pcatch;

#ifdef PYMOL_NEW_THREADS
  PyEval_InitThreads();
#endif

  /* BEGIN PROPRIETARY CODE SEGMENT (see disclaimer in "os_proprietary.h") */
#ifdef WIN32
#ifdef _PYMOL_MONOLITHIC
#ifndef _PYMOL_EMBEDDED
#define _PYMOL_INIT_MODULES
#endif
#endif
#endif
  /* END PROPRIETARY CODE SEGMENT */

#ifdef _PYMOL_INIT_MODULES
  /* Win32 module build: includes pyopengl, numpy, and sglite */
  /* sglite 
   * initExtensionClass();
   * initsglite();
   */
  init_champ();
#ifdef _PYMOL_PYOMM
  init_pyomm();
#endif

  /* initialize numeric python */
#ifdef _PYMOL_NUMPY_INIT
  init_numpy();
  initmultiarray();
  initarrayfns();
  initlapack_lite();
  initumath();
  initranlib();
#endif
  /* initialize PyOpenGL */
#endif

  if(global_instance) {
    PCatchInit();               /* setup standard-output catch routine */
  }

  /* assumes that pymol module has been loaded via Python or PyRun_SimpleString */

  P_pymol = PyImport_AddModule("pymol");        /* get it */
  if(!P_pymol)
    ErrFatal(G, "PyMOL", "can't find module 'pymol'");
  P_pymol_dict = PyModule_GetDict(P_pymol);
  Py_XINCREF(P_pymol_dict);
  if(!P_pymol_dict)
    ErrFatal(G, "PyMOL", "can't find globals for 'pymol'");

  if(global_instance) {         /* if global singleton PyMOL... */
    G->P_inst = Calloc(CP_inst, 1);
    G->P_inst->obj = P_pymol;
    G->P_inst->dict = P_pymol_dict;
    {
      int a;
      SavedThreadRec *str = G->P_inst->savedThread;
      for(a = 0; a < MAX_SAVED_THREAD; a++) {
        (str++)->id = -1;
      }
    }
  }

  {
    G->P_inst->exec = PyDict_GetItemString(P_pymol_dict, "exec_str");
    Py_XINCREF(G->P_inst->exec);
    if(!G->P_inst->exec)
      ErrFatal(G, "PyMOL", "can't find 'pymol.exec_str()'");

    sys = PyDict_GetItemString(P_pymol_dict, "sys");
    Py_XINCREF(sys);
    if(!sys)
      ErrFatal(G, "PyMOL", "can't find 'pymol.sys'");

    if(global_instance) {

      /* implies global singleton pymol, so set up the global handle */
      PyDict_SetItemString(P_pymol_dict, "_COb",
                           PyCObject_FromVoidPtr((void *) &SingletonPyMOLGlobals, NULL));

      pcatch = PyImport_AddModule("pcatch");
      if(!pcatch)
        ErrFatal(G, "PyMOL", "can't find module 'pcatch'");
      PyObject_SetAttrString(sys, "stdout", pcatch);
      PyObject_SetAttrString(sys, "stderr", pcatch);
    }

    PRunStringModule(G, "import traceback\n");
    P_traceback = PyDict_GetItemString(P_pymol_dict, "traceback");
    Py_XINCREF(P_traceback);
    if(!P_traceback)
      ErrFatal(G, "PyMOL", "can't find 'traceback'");

    PRunStringModule(G, "import cmd\n");
    P_cmd = PyDict_GetItemString(P_pymol_dict, "cmd");
    Py_XINCREF(P_cmd);
    if(!P_cmd)
      ErrFatal(G, "PyMOL", "can't find 'cmd'");

    if(global_instance) {
      /* implies global singleton pymol, so set up the global handle */
      PyObject_SetAttrString(P_cmd, "_COb",
                             PyCObject_FromVoidPtr((void *) &SingletonPyMOLGlobals,
                                                   NULL));

      /* cmd module is itself the api for the global PyMOL instance */
      G->P_inst->cmd = P_cmd;
    }

    PyObject_SetAttrString(G->P_inst->cmd, "_pymol", G->P_inst->obj);

    /* right now, all locks are global -- eventually some of these may
       become instance-specific in order to improve concurrency */

    G->P_inst->lock = PyObject_GetAttrString(G->P_inst->cmd, "lock");
    if(!G->P_inst->lock)
      ErrFatal(G, "PyMOL", "can't find 'cmd.lock()'");

    G->P_inst->lock_attempt = PyObject_GetAttrString(G->P_inst->cmd, "lock_attempt");
    if(!G->P_inst->lock_attempt)
      ErrFatal(G, "PyMOL", "can't find 'cmd.lock_attempt()'");

    G->P_inst->unlock = PyObject_GetAttrString(G->P_inst->cmd, "unlock");
    if(!G->P_inst->unlock)
      ErrFatal(G, "PyMOL", "can't find 'cmd.unlock()'");

    G->P_inst->lock_c = PyObject_GetAttrString(G->P_inst->cmd, "lock_c");
    if(!G->P_inst->lock_c)
      ErrFatal(G, "PyMOL", "can't find 'cmd.lock_c()'");

    G->P_inst->unlock_c = PyObject_GetAttrString(G->P_inst->cmd, "unlock_c");
    if(!G->P_inst->unlock_c)
      ErrFatal(G, "PyMOL", "can't find 'cmd.unlock_c()'");

    G->P_inst->lock_status = PyObject_GetAttrString(G->P_inst->cmd, "lock_status");
    if(!G->P_inst->lock_status)
      ErrFatal(G, "PyMOL", "can't find 'cmd.lock_status()'");

    G->P_inst->lock_status_attempt =
      PyObject_GetAttrString(G->P_inst->cmd, "lock_status_attempt");
    if(!G->P_inst->lock_status_attempt)
      ErrFatal(G, "PyMOL", "can't find 'cmd.lock_status_attempt()'");

    G->P_inst->unlock_status = PyObject_GetAttrString(G->P_inst->cmd, "unlock_status");
    if(!G->P_inst->unlock_status)
      ErrFatal(G, "PyMOL", "can't find 'cmd.unlock_status()'");

    G->P_inst->lock_glut = PyObject_GetAttrString(G->P_inst->cmd, "lock_glut");
    if(!G->P_inst->lock_glut)
      ErrFatal(G, "PyMOL", "can't find 'cmd.lock_glut()'");

    G->P_inst->unlock_glut = PyObject_GetAttrString(G->P_inst->cmd, "unlock_glut");
    if(!G->P_inst->unlock_glut)
      ErrFatal(G, "PyMOL", "can't find 'cmd.unlock_glut()'");

    /* 'do' command */

    G->P_inst->cmd_do = PyObject_GetAttrString(G->P_inst->cmd, "do");
    if(!G->P_inst->cmd_do)
      ErrFatal(G, "PyMOL", "can't find 'cmd.do()'");

    /* cache */
    G->P_inst->cache = PyObject_GetAttrString(G->P_inst->obj, "_cache");

    /* invariant stuff */

    PRunStringModule(G, "import menu\n");
    P_menu = PyDict_GetItemString(P_pymol_dict, "menu");
    Py_XINCREF(P_menu);
    if(!P_menu)
      ErrFatal(G, "PyMOL", "can't find module 'menu'");

    PRunStringModule(G, "import setting\n");
    P_setting = PyDict_GetItemString(P_pymol_dict, "setting");
    Py_XINCREF(P_setting);
    if(!P_setting)
      ErrFatal(G, "PyMOL", "can't find module 'setting'");

    PRunStringModule(G, "import povray\n");
    P_povray = PyDict_GetItemString(P_pymol_dict, "povray");
    Py_XINCREF(P_povray);
    if(!P_povray)
      ErrFatal(G, "PyMOL", "can't find module 'povray'");

#ifdef _PYMOL_XRAY
    PRunStringModule(G, "import xray\n");
    P_xray = PyDict_GetItemString(P_pymol_dict, "xray");
    Py_XINCREF(P_xray);
    if(!P_xray)
      ErrFatal(G, "PyMOL", "can't find module 'xray'");
#endif

    /* BEGIN PROPRIETARY CODE SEGMENT (see disclaimer in "os_proprietary.h") */
#ifdef WIN32
    PRunStringModule(G, "import time\n");
    P_time = PyDict_GetItemString(P_pymol_dict, "time");
    Py_XINCREF(P_time);
    if(!P_time)
      ErrFatal(G, "PyMOL", "can't find module 'time'");

    P_sleep = PyObject_GetAttrString(P_time, "sleep");
    Py_XINCREF(P_sleep);
    if(!P_sleep)
      ErrFatal(G, "PyMOL", "can't find 'time.sleep()'");
#endif
    /* END PROPRIETARY CODE SEGMENT */

    PRunStringModule(G, "import parser\n");
    P_parser = PyDict_GetItemString(P_pymol_dict, "parser");
    Py_XINCREF(P_parser);
    if(!P_parser)
      ErrFatal(G, "PyMOL", "can't find module 'parser'");

    {
      PyObject *fn_closure = PyObject_GetAttrString(P_parser, "new_parse_closure");
      G->P_inst->parse = PyObject_CallFunction(fn_closure, "O", G->P_inst->cmd);
      PXDecRef(fn_closure);
      if(!G->P_inst->parse)
        ErrFatal(G, "PyMOL", "can't create 'parse' function closure");
    }

    {
      PyObject *fn_closure = PyObject_GetAttrString(P_parser, "new_complete_closure");
      G->P_inst->complete = PyObject_CallFunction(fn_closure, "O", G->P_inst->cmd);
      PXDecRef(fn_closure);
      if(!G->P_inst->complete)
        ErrFatal(G, "PyMOL", "can't create 'complete' function closure");
    }

    PRunStringModule(G, "import chempy");
    P_chempy = PyDict_GetItemString(P_pymol_dict, "chempy");
    Py_XINCREF(P_chempy);
    if(!P_chempy)
      ErrFatal(G, "PyMOL", "can't find 'chempy'");

    PRunStringModule(G, "from chempy.bonds import bonds");      /* load bond dictionary */

    PRunStringModule(G, "from chempy import models");
    P_models = PyDict_GetItemString(P_pymol_dict, "models");
    Py_XINCREF(P_models);
    if(!P_models)
      ErrFatal(G, "PyMOL", "can't find 'chempy.models'");

    PRunStringModule(G, "import util\n");
    PRunStringModule(G, "import preset\n");
    PRunStringModule(G, "import contrib\n");

    PRunStringModule(G, "import string\n");

    /* backwards compatibility */

    PRunStringModule(G, "pm = cmd\n");
    PRunStringModule(G, "pmu = util\n");

    PRunStringModule(G, "glutThread = thread.get_ident()");

    P_glut_thread_id = PyThread_get_thread_ident();

#ifndef WIN32
    if(G->Option->siginthand) {
      signal(SIGINT, my_interrupt);
    }
#endif

    /* required environment variables */

    PyRun_SimpleString("import os");

    PyRun_SimpleString
      ("if not os.environ.has_key('PYMOL_DATA'): os.environ['PYMOL_DATA']=os.environ['PYMOL_PATH']+'/data'");
    PyRun_SimpleString("os.environ['TUT']=os.environ['PYMOL_DATA']+'/tut'");

    PyRun_SimpleString
      ("if not os.environ.has_key('PYMOL_SCRIPTS'): os.environ['PYMOL_SCRIPTS']=os.environ['PYMOL_PATH']+'/scripts'");

  }
}

int PPovrayRender(PyMOLGlobals * G, char *header, char *inp, char *file, int width,
                  int height, int antialias)
{
  PyObject *result;
  int ok;
  PBlock(G);
  result =
    PyObject_CallMethod(P_povray, "render_from_string", "sssiii", header, inp, file,
                        width, height, antialias);
  ok = PyObject_IsTrue(result);
  Py_DECREF(result);
  PUnblock(G);
  return (ok);
}

void PSGIStereo(PyMOLGlobals * G, int flag)
{
  int blocked;
  blocked = PAutoBlock(G);
  if(flag)
    PRunStringModule(G, "cmd._sgi_stereo(1)");
  else
    PRunStringModule(G, "cmd._sgi_stereo(0)");
  if(blocked)
    PUnblock(G);
}

void PFree(void)
{
}

void PExit(PyMOLGlobals * G, int code)
{
  ExecutiveDelete(G, "all");
  PBlock(G);
#ifndef _PYMOL_NO_MAIN
  if(G->Main) {
    MainFree();
  }
#endif

  /* we're having trouble with threading errors after calling Py_Exit,
     so for the time being, let's just take the process down at this
     point, instead of allowing PyExit to be called. */

  exit(EXIT_SUCCESS);

  Py_Exit(code);
}

void PParse(PyMOLGlobals * G, char *str)
{
  OrthoCommandIn(G, str);
}

void PDo(PyMOLGlobals * G, char *str)
{                               /* assumes we already hold the re-entrant API lock */
  int blocked;
  blocked = PAutoBlock(G);
  Py_XDECREF(PyObject_CallFunction(G->P_inst->cmd_do, "s", str));
  PAutoUnblock(G, blocked);
}

void PLog(PyMOLGlobals * G, char *str, int format)

/* general log routine can write PML 

   or PYM commands to appropriate log file */
{
  int mode;
  int a;
  int blocked;
  PyObject *log;
  OrthoLineType buffer = "";
  mode = (int) SettingGet(G, cSetting_logging);
  if(mode) {
    blocked = PAutoBlock(G);
    log = PyDict_GetItemString(P_pymol_dict, P_log_file_str);
    if(log && (log != Py_None)) {
      if(format == cPLog_no_flush) {
        PyObject_CallMethod(log, "write", "s", str);    /* maximize responsiveness (for real-time) */
      } else {
        switch (mode) {
        case cPLog_pml:        /* .pml file */
          switch (format) {
          case cPLog_pml_lf:
            strcpy(buffer, str);
            break;
          case cPLog_pml:
          case cPLog_pym:
            strcpy(buffer, str);
            strcat(buffer, "\n");
            break;
          }
          break;
        case cPLog_pym:        /* .pym file */
          if((str[0] == '_') && (str[1]) == ' ')
            str += 2;
          switch (format) {
          case cPLog_pml_lf:
            a = strlen(str);
            while(a) {          /* trim CR/LF etc. */
              if(*(str + a) >= 32)
                break;
              *(str + a) = 0;
              a--;
            }
          case cPLog_pml:
            strcpy(buffer, "cmd.do('''");
            strcat(buffer, str);
            strcat(buffer, "''')\n");
            break;
          case cPLog_pym:
            strcpy(buffer, str);
            strcat(buffer, "\n");
            break;
          }
        }
        PyObject_CallMethod(log, "write", "s", buffer);
        PyObject_CallMethod(log, "flush", "");
      }
    }
    PAutoUnblock(G, blocked);
  }
}

void PLogFlush(PyMOLGlobals * G)
{
  int mode;
  PyObject *log;
  int blocked;
  mode = (int) SettingGet(G, cSetting_logging);
  if(mode) {
    blocked = PAutoBlock(G);
    log = PyDict_GetItemString(P_pymol_dict, P_log_file_str);
    if(log && (log != Py_None)) {
      PyObject_CallMethod(log, "flush", "");
    }
    PAutoUnblock(G, blocked);
  }
}

int PFlush(PyMOLGlobals * G)
{
  /* NOTE: ASSUMES unblocked Python threads and a locked API */
  PyObject *err;
  char buffer[OrthoLineLength + 1];
  int did_work = false;
  if(OrthoCommandWaiting(G)) {
    did_work = true;
    PBlock(G);
    if(!(PIsGlutThread() && G->P_inst->glut_thread_keep_out)) {
      /* don't run if we're currently banned */
      while(OrthoCommandOut(G, buffer)) {
        OrthoCommandNest(G, 1);
        PUnlockAPIWhileBlocked(G);
        if(PyErr_Occurred()) {
          PyErr_Print();
          PRINTFB(G, FB_Python, FB_Errors)
            " PFlush: Uncaught exception.  PyMOL may have a bug.\n" ENDFB(G);
        }
        PXDecRef(PyObject_CallFunction(G->P_inst->parse, "si", buffer, 0));
        err = PyErr_Occurred();
        if(err) {
          PyErr_Print();
          PRINTFB(G, FB_Python, FB_Errors)
            " PFlush: Uncaught exception.  PyMOL may have a bug.\n" ENDFB(G);
        }
        PLockAPIWhileBlocked(G);
        /* make sure no commands left at this level */
        while(OrthoCommandWaiting(G))
          PFlushFast(G);
        OrthoCommandNest(G, -1);
      }
    }
    PUnblock(G);
  }
  return did_work;
}

int PFlushFast(PyMOLGlobals * G)
{
  /* NOTE: ASSUMES we currently have blocked Python threads and an unlocked API */
  PyObject *err;
  char buffer[OrthoLineLength + 1];
  int did_work = false;
  while(OrthoCommandOut(G, buffer)) {
    OrthoCommandNest(G, 1);
    did_work = true;
    PRINTFD(G, FB_Threads)
      " PFlushFast-DEBUG: executing '%s' as thread 0x%x\n", buffer,
      PyThread_get_thread_ident()
      ENDFD;
    if(PyErr_Occurred()) {
      PyErr_Print();
      PRINTFB(G, FB_Python, FB_Errors)
        " PFlushFast: Uncaught exception.  PyMOL may have a bug.\n" ENDFB(G);
    }
    PXDecRef(PyObject_CallFunction(G->P_inst->parse, "si", buffer, 0));
    err = PyErr_Occurred();
    if(err) {
      PyErr_Print();
      PRINTFB(G, FB_Python, FB_Errors)
        " PFlushFast: Uncaught exception.  PyMOL may have a bug.\n" ENDFB(G);
    }
    /* make sure no commands left at this level */
    while(OrthoCommandWaiting(G))
      PFlushFast(G);
    OrthoCommandNest(G, -1);
  }
  return did_work;
}

void PBlockLegacy()
{
  PBlock(SingletonPyMOLGlobals);
}

void PUnblockLegacy()
{
  PUnblock(SingletonPyMOLGlobals);
}

void PBlock(PyMOLGlobals * G)
{

  if(!PAutoBlock(G)) {
    // int *p = 0;
    //  *p = 0;
    ErrFatal(G, "PBlock", "Threading error detected.  Terminating...");
  }
}

int PAutoBlock(PyMOLGlobals * G)
{
#ifndef _PYMOL_EMBEDDED
  int a, id;
  SavedThreadRec *SavedThread = G->P_inst->savedThread;
  /* synchronize python */

  id = PyThread_get_thread_ident();
  PRINTFD(G, FB_Threads)
    " PAutoBlock-DEBUG: search 0x%x (0x%x, 0x%x, 0x%x)\n", id,
    SavedThread[MAX_SAVED_THREAD - 1].id,
    SavedThread[MAX_SAVED_THREAD - 2].id, SavedThread[MAX_SAVED_THREAD - 3].id ENDFD;
  a = MAX_SAVED_THREAD - 1;
  while(a) {
    if(!((SavedThread + a)->id - id)) {
      /* astoundingly, equality test fails on ALPHA even 
       * though the ints are equal. Must be some kind of optimizer bug
       * or mis-assumption */

      PRINTFD(G, FB_Threads)
        " PAutoBlock-DEBUG: seeking global lock 0x%x\n", id ENDFD;

#ifdef PYMOL_NEW_THREADS

      PyEval_AcquireLock();

      PRINTFD(G, FB_Threads)
        " PAutoBlock-DEBUG (NewThreads): restoring 0x%x\n", id ENDFD;

      PyThreadState_Swap((SavedThread + a)->state);

#else
      PRINTFD(G, FB_Threads)
        " PAutoBlock-DEBUG: restoring 0x%x\n", id ENDFD;

      PyEval_RestoreThread((SavedThread + a)->state);
#endif

      PRINTFD(G, FB_Threads)
        " PAutoBlock-DEBUG: restored 0x%x\n", id ENDFD;

      PRINTFD(G, FB_Threads)
        " PAutoBlock-DEBUG: clearing 0x%x\n", id ENDFD;

      PXDecRef(PyObject_CallFunction(G->P_inst->lock_c, "O", G->P_inst->cmd));
      SavedThread[a].id = -1;
      /* this is the only safe time we can change things */
      PXDecRef(PyObject_CallFunction(G->P_inst->unlock_c, "O", G->P_inst->cmd));

      PRINTFD(G, FB_Threads)
        " PAutoBlock-DEBUG: blocked 0x%x (0x%x, 0x%x, 0x%x)\n",
        PyThread_get_thread_ident(), SavedThread[MAX_SAVED_THREAD - 1].id,
        SavedThread[MAX_SAVED_THREAD - 2].id, SavedThread[MAX_SAVED_THREAD - 3].id ENDFD;

      return 1;
    }
    a--;
  }
  PRINTFD(G, FB_Threads)
    " PAutoBlock-DEBUG: 0x%x not found, thus already blocked.\n",
    PyThread_get_thread_ident()
    ENDFD;
  return 0;
#else
  return 1;
#endif
}

int PIsGlutThread(void)
{
  return (PyThread_get_thread_ident() == P_glut_thread_id);
}

void PUnblock(PyMOLGlobals * G)
{
#ifndef _PYMOL_EMBEDDED
  int a;
  SavedThreadRec *SavedThread = G->P_inst->savedThread;
  /* NOTE: ASSUMES a locked API */
  PRINTFD(G, FB_Threads)
    " PUnblock-DEBUG: entered as thread 0x%x\n", PyThread_get_thread_ident()
    ENDFD;

  /* reserve a space while we have a lock */
  PXDecRef(PyObject_CallFunction(G->P_inst->lock_c, "O", G->P_inst->cmd));
  a = MAX_SAVED_THREAD - 1;
  while(a) {
    if((SavedThread + a)->id == -1) {
      (SavedThread + a)->id = PyThread_get_thread_ident();
#ifdef PYMOL_NEW_THREADS
      (SavedThread + a)->state = PyThreadState_Get();
#endif
      break;
    }
    a--;
  }
  PRINTFD(G, FB_Threads)
    " PUnblock-DEBUG: 0x%x stored in slot %d\n", (SavedThread + a)->id, a ENDFD;
  PXDecRef(PyObject_CallFunction(G->P_inst->unlock_c, "O", G->P_inst->cmd));
#ifdef PYMOL_NEW_THREADS
  PyThreadState_Swap(NULL);
  PyEval_ReleaseLock();
#else
  (SavedThread + a)->state = PyEval_SaveThread();
#endif
#endif
}

void PAutoUnblock(PyMOLGlobals * G, int flag)
{
  if(flag)
    PUnblock(G);
}

void PBlockAndUnlockAPI(PyMOLGlobals * G)
{
  PBlock(G);
  PXDecRef(PyObject_CallFunction(G->P_inst->unlock, "iO", 0, G->P_inst->cmd));
}

int PLockAPI(PyMOLGlobals * G, int block_if_busy)
{
  int result = true;
  PBlock(G);
  if(block_if_busy) {
    PXDecRef(PyObject_CallFunction(G->P_inst->lock, "O", G->P_inst->cmd));
  } else {                      /* not blocking if PyMOL is busy */

    PyObject *got_lock =
      PyObject_CallFunction(G->P_inst->lock_attempt, "O", G->P_inst->cmd);

    if(got_lock) {
      result = PyInt_AsLong(got_lock);
      Py_DECREF(got_lock);
    }
  }
  PUnblock(G);
  return result;
}

void PUnlockAPI(PyMOLGlobals * G)
{
  PBlock(G);
  PXDecRef(PyObject_CallFunction(G->P_inst->unlock, "iO", 0, G->P_inst->cmd));
  PUnblock(G);
}

static void PUnlockAPIWhileBlocked(PyMOLGlobals * G)
{
  PXDecRef(PyObject_CallFunction(G->P_inst->unlock, "iO", -1, G->P_inst->cmd));
}

static void PLockAPIWhileBlocked(PyMOLGlobals * G)
{
  PXDecRef(PyObject_CallFunction(G->P_inst->lock, "O", G->P_inst->cmd));
}

int PTryLockAPIAndUnblock(PyMOLGlobals * G)
{
  int result = get_api_lock(G, false);
  if(result) {
    PUnblock(G);
  }
  return result;
}

void PLockAPIAndUnblock(PyMOLGlobals * G)
{
  PXDecRef(PyObject_CallFunction(G->P_inst->lock, "O", G->P_inst->cmd));
  PUnblock(G);
}

void PDefineFloat(PyMOLGlobals * G, char *name, float value)
{
  char buffer[OrthoLineLength];
  sprintf(buffer, "%s = %f\n", name, value);
  PBlock(G);
  PRunStringModule(G, buffer);
  PUnblock(G);
}


/* This function is called by the interpreter to get its own name */
char *getprogramname(void)
{
  return ("PyMOL");
}


/* A static module */

static PyObject *PCatchWrite(PyObject * self, PyObject * args)
{
  char *str;
  PyArg_ParseTuple(args, "s", &str);
  if(str[0]) {
    if(SingletonPyMOLGlobals) {
      if(Feedback(SingletonPyMOLGlobals, FB_Python, FB_Output)) {
        OrthoAddOutput(SingletonPyMOLGlobals, str);
      }
    }
  }
  return PConvAutoNone(Py_None);
}

static PyObject *PCatchWritelines(PyObject * self, PyObject * args)
{
  PyObject *seq;
  int len;
  PyArg_ParseTuple(args, "O", &seq);
  if(seq && PySequence_Check(seq)) {
    if((len = PySequence_Size(seq)) > 0) {
      int i;
      for(i = 0; i < len; i++) {
        PyObject *obj = PySequence_GetItem(seq, i);
        if(obj && PyString_Check(obj)) {
          char *str = PyString_AsString(obj);
          if(SingletonPyMOLGlobals) {
            if(Feedback(SingletonPyMOLGlobals, FB_Python, FB_Output)) {
              OrthoAddOutput(SingletonPyMOLGlobals, str);
            }
          }
        }
        Py_XDECREF(obj);
      }
    }
  }
  return PConvAutoNone(Py_None);
}

static PyObject *PCatchFlush(PyObject * self, PyObject * args)
{
  fflush(stdout);
  fflush(stderr);
  return PConvAutoNone(Py_None);
}

static PyMethodDef PCatch_methods[] = {
  {"writelines", PCatchWritelines, METH_VARARGS},
  {"write", PCatchWrite, METH_VARARGS},
  {"flush", PCatchFlush, METH_VARARGS},
  {NULL, NULL}                  /* sentinel */
};

void PCatchInit(void)
{
  PyImport_AddModule("pcatch");
  Py_InitModule("pcatch", PCatch_methods);
}
#endif
