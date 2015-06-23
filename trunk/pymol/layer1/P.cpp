
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
#include"os_proprietary.h"
#include<process.h>
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

static int label_copy_text(char *dst, const char *src, int len, int max)
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

static int label_next_token(WordType dst, const char **expr)
{
  const char *p = *expr;
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

#if 0
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
#endif


int PLabelExprUsesVariable(PyMOLGlobals * G, const char *expr, const char *var)
{
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

int PLabelAtomAlt(PyMOLGlobals * G, AtomInfoType * at, const char *model, const char *expr, int index)
{
  /* alternate C implementation which bypasses Python expressions -- works
     only for simple label formats "..."+property+... */

  int result = false;
  OrthoLineType label;
  int label_len = 0;
  int label_max = sizeof(OrthoLineType);
  OrthoLineType buffer;
  char ch, quote = 0;
  int escaped = false;
  const char *origexpr = expr;

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
	int tokresult = true;
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
            label_len = label_copy_text(label, LexStr(G, at->chain), label_len, label_max);
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
            const char *st = LexStr(G, at->textType);
            label_len = label_copy_text(label, st, label_len, label_max);
          } else if(!strcmp(tok, "custom")) {
            const char *st = LexStr(G, at->custom);
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
          } else {
	    tokresult = false;
	  }
          if(buffer[0]) {
            label_len = label_copy_text(label, buffer, label_len, label_max);
          }
        } else {
	  if (tok[0]){
	    label_len = label_copy_text(label, "?", label_len, label_max);
	    label_len = label_copy_text(label, tok, label_len, label_max);
	  } else {
	    tokresult = false;
	  }
	}
	result |= tokresult;
      } else {
        if(label_len < label_max) {
          label[label_len] = '?';
          label_len++;
	  result = true;
        }
      }
    } else {
      if(ch == quote) {
        quote = 0;
	result = true;
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

  if (!result && !label[0]){
    // if label is not set, just use expression as a string for label
    strncpy(label, origexpr, OrthoLineLength);
    result = true;
  }

  LexDec(G, at->label);
  at->label = result ? LexIdx(G, label) : 0;

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

static PyMappingMethods wrapperMappingMethods;
static PyTypeObject Wrapper_Type = {
  PyObject_HEAD_INIT(NULL)
};

static PyObject* PyObject_GenericGetAttrOrItem(PyObject *o, PyObject *key) {
  PyObject *ret = PyObject_GenericGetAttr(o, key);
  if (!PyErr_Occurred())
    return ret;
  PyErr_Clear();
  return PyObject_GetItem(o, key);
}

int PyObject_GenericSetAttrAsItem(PyObject *o, PyObject *key, PyObject *value) {
  return PyObject_SetItem(o, key, value);
}

PyObject * WrapperObjectSubScript(PyObject *obj, PyObject *key){
  WrapperObject *wobj = (WrapperObject*)obj;
  char *aprop;
  AtomPropertyInfo *ap;
  PyObject *ret = NULL;
  bool borrowed = false;

  if (!wobj || !wobj->obj){
    /* ERROR */
    PRINTFB(wobj->G, FB_Python, FB_Errors)
      "Error: wrappers cannot be used outside of the iterate/alter/alter_state commands\n" ENDFB(wobj->G);
    return NULL;
  }
  PyObject *keyobj = PyObject_Str(key);
  aprop = PyString_AS_STRING(keyobj);
  ap = PyMOL_GetAtomPropertyInfo(wobj->G->PyMOL, aprop);
  Py_DECREF(keyobj);
  if (ap){
    switch (ap->Ptype){
    case cPType_string:
      {
	char *val = (char*)(((char*)wobj->atomInfo) + ap->offset);
	ret = PyString_FromString(val);
      }
      break;
    case cPType_schar:
      {
	signed char val = *(signed char*)(((char*)wobj->atomInfo) + ap->offset);
	ret = PyInt_FromLong((long)val);
      }
      break;
    case cPType_int:
      {
	int val = *(int*)(((char*)wobj->atomInfo) + ap->offset);
	ret = PyInt_FromLong((long)val);
      }
      break;
    case cPType_int_as_string:
      {
	int val = *(int*)(((char*)wobj->atomInfo) + ap->offset);
	const char *st = LexStr(wobj->G, val);
	ret = PyString_FromString(st);
      }
      break;
    case cPType_float:
      {
	float val = *(float*)(((char*)wobj->atomInfo) + ap->offset);
	ret = PyFloat_FromDouble(val);
      }
      break;
    case cPType_stereo:
      {
	char val = *(char*)(((char*)wobj->atomInfo) + ap->offset);	
	char mmstereotype[2];
	mmstereotype[0] = convertStereoToChar(val);
	mmstereotype[1] = 0;
	ret = PyString_FromString(mmstereotype);
      }
      break;
    case cPType_char_as_type:
      {
	char val = *(char*)(((char*)wobj->atomInfo) + ap->offset);
	char atype[7];
	if(val)
	  strcpy(atype, "HETATM");
	else
	  strcpy(atype, "ATOM");
	ret = PyString_FromString(atype);
      }
      break;
    case cPType_model:
      {
	if (wobj->model){
	  ret = PyString_FromString(wobj->model);
	}
      }
      break;
    case cPType_index:
      {
	ret = PyInt_FromLong((long)wobj->index + 1);
      }
      break;
    case cPType_int_custom_type:
      {
	int val = *(int*)(((char*)wobj->atomInfo) + ap->offset);
	if(val != cAtomInfoNoType){
	  ret = PyInt_FromLong((long)val);
	} else {
	  ret = PyString_FromString("?");
	}
      }
      break;
    case cPType_xyz_float:
      {
	if (wobj->v){
	  ret = PyFloat_FromDouble(wobj->v[ap->offset]);
	} else {
	  PRINTFB(wobj->G, FB_Python, FB_Errors)
	    " PLabelAtom/PAlterAtom: Warning: x/y/z cannot be set in alter/iterate/label function\n" ENDFB(wobj->G);
	}
      }
      break;
    case cPType_settings:
    case cPType_properties:
      PRINTFB(wobj->G, FB_Python, FB_Errors)
        " PLabelAtom/PAlterAtom: Warning: Accessing settings and properties not supported in Open-Source PyMOL\n" ENDFB(wobj->G);
      break;
    case cPType_state:
      {
	if (wobj->v){
	  ret = PyInt_FromLong((long)wobj->state);
	} else {
	  PRINTFB(wobj->G, FB_Python, FB_Errors)
	    " PLabelAtom/PAlterAtom: Warning: state cannot be set in alter/iterate/label function\n" ENDFB(wobj->G);
	}
      }
    }
  } else {
    /* if not an atom property, check if local variable in dict */
    ret = PyDict_GetItem(wobj->dict, key);
    borrowed = true;
  }

  if (borrowed)
    PXIncRef(ret);
  return ret;
}

float PyObject_as_Number(PyObject *obj){
  if (PyInt_Check(obj)){
    return PyInt_AS_LONG(obj);
  } else if (PyFloat_Check(obj)){
    return PyFloat_AS_DOUBLE(obj);
  } 
  return 0.f;
}

int WrapperObjectAssignSubScript(PyObject *obj, PyObject *key, PyObject *val){
  WrapperObject *wobj = (WrapperObject*)obj;
  if (!wobj || !wobj->obj){
    /* ERROR */
    PRINTFB(wobj->G, FB_Python, FB_Errors)
      "Error: wrappers cannot be used outside of the iterate/alter/alter_state commands\n" ENDFB(wobj->G);
    return -1;
  }
  {
    PyObject *keyobj = PyObject_Str(key);
    char *aprop = PyString_AS_STRING(keyobj);
    AtomPropertyInfo *ap;
    ap = PyMOL_GetAtomPropertyInfo(wobj->G->PyMOL, aprop);
    Py_DECREF(keyobj);
    if (ap){
      short changed = false;
      if (wobj->read_only){
	PRINTFB(wobj->G, FB_Python, FB_Errors)
	  " PIterateAtom/PLabelAtom: Warning: properties cannot be set in iterate/label function\n" ENDFB(wobj->G);
	return -1;
      }
      if (wobj->v){ // PAlterAtomState, must be setting x/y/z or flags
	if (ap->Ptype != cPType_xyz_float && ap->id != ATOM_PROP_FLAGS){
	  PRINTFB(wobj->G, FB_Python, FB_Errors)
	    " PAlterAtomState: Warning: properties cannot be set in alter_state function\n" ENDFB(wobj->G);
	  return -1;
	}
      }

      switch (ap->Ptype){
      case cPType_string:
	{
          PyObject *valobj = PyObject_Str(val);
	  char *valstr = PyString_AS_STRING(valobj);
	  char *dest = (char*)(((char*)wobj->atomInfo) + ap->offset);
	  if (strlen(valstr) > ap->maxlen){
	    strncpy(dest, valstr, ap->maxlen);
	  } else {
	    strcpy(dest, valstr);
	  }
          Py_DECREF(valobj);
	  changed = true;
	}
	break;
      case cPType_schar:
	{
	  int valint = PyInt_AsLong(val);
	  signed char *dest;
	  if (valint == -1 && PyErr_Occurred())
	    break;
	  dest = (signed char*)(((char*)wobj->atomInfo) + ap->offset);
	  *dest = valint;
	  changed = true;
	}
        break;
      case cPType_int:
	{
	  int valint = PyInt_AsLong(val);
	  int *dest;
	  if (valint == -1 && PyErr_Occurred())
	    break;
	  dest = (int*)(((char*)wobj->atomInfo) + ap->offset);
	  *dest = valint;
	  changed = true;
	}
	break;
      case cPType_int_as_string:
	{
	  int *dest = (int*)(((char*)wobj->atomInfo) + ap->offset);
          PyObject *valobj = PyObject_Str(val);
	  char *valstr = PyString_AS_STRING(valobj);
	  LexDec(wobj->G, *dest);
	  *dest = LexIdx(wobj->G, valstr);
          Py_DECREF(valobj);
	  changed = true;
	}
	break;
      case cPType_float:
	{
	  float valfloat = PyObject_as_Number(val);
	  float *dest = (float*)(((char*)wobj->atomInfo) + ap->offset);
	  *dest = valfloat;
	  changed = true;
	}
	break;	
      case cPType_stereo:
	{
          PyObject *valobj = PyObject_Str(val);
	  char *stereo = PyString_AS_STRING(valobj);
	  if (stereo){
	    wobj->atomInfo->mmstereo = convertStereoToChar(stereo[0]);
	  }
          Py_DECREF(valobj);
	}
	break;
      case cPType_char_as_type:
	{
	  signed char *dest = (signed char *)(((char*)wobj->atomInfo) + ap->offset);
          PyObject *valobj = PyObject_Str(val);
          char *valstr = PyString_AS_STRING(valobj);
	  *dest = ((valstr[0] == 'h') || (valstr[0] == 'H'));
          Py_DECREF(valobj);
	  changed = true;
	}
	break;
      case cPType_model:
	PRINTFB(wobj->G, FB_Python, FB_Errors)
	  " PAlterAtom: Warning: cannot change model value for individual atoms, please rename object\n" ENDFB(wobj->G);
	break;
      case cPType_index:
	PRINTFB(wobj->G, FB_Python, FB_Errors)
	  " PAlterAtom: Warning: cannot change index value for atoms\n" ENDFB(wobj->G);
	break;
      case cPType_int_custom_type:
	{
          PyObject *valobj = PyObject_Str(val);
	  char *valstr = PyString_AS_STRING(valobj);
	  int *dest = (int*)(((char*)wobj->atomInfo) + ap->offset);
	  if (valstr[0] == '?'){
	    *dest = cAtomInfoNoType;
	  } else {
	    int valint = PyInt_AS_LONG(val);
	    *dest = valint;
	  }
          Py_DECREF(valobj);
	  changed = true;
	}
	break;
      case cPType_xyz_float:
	if (wobj->v){
	  float valfloat = PyObject_as_Number(val);
	  wobj->v[ap->offset] = valfloat;
	} else {
	  PRINTFB(wobj->G, FB_Python, FB_Errors)
	    " PAlterAtom: Warning: can only set x/y/z in alter_state\n" ENDFB(wobj->G);
	}
	break;
      case cPType_settings:
      case cPType_properties:
	PRINTFB(wobj->G, FB_Python, FB_Errors)
	  " PLabelAtom/PAlterAtom/PIterateAtom: Warning: settings or properties object cannot be set\n" ENDFB(wobj->G);
      }
      if (changed){
	switch (ap->id){
	case ATOM_PROP_ELEM:
	  AtomInfoAssignParameters(wobj->G, wobj->atomInfo);
	  break;
	case ATOM_PROP_RESI:
	  wobj->atomInfo->resv = AtomResvFromResi(wobj->atomInfo->resi);
	  break;
	case ATOM_PROP_RESV:
	  {
	    WordType buf;
	    sprintf(buf, "%d", wobj->atomInfo->resv);
	    buf[sizeof(ResIdent) - 1] = 0;
	    strcpy(wobj->atomInfo->resi, buf);
	  }
	  break;
	case ATOM_PROP_SS:
	  wobj->atomInfo->ssType[0] = toupper(wobj->atomInfo->ssType[0]);
	  break;
	case ATOM_PROP_FORMAL_CHARGE:
	  wobj->atomInfo->chemFlag = false;
	  break;
	}
      }
    } else {
      /* if not an atom property, then its a local variable, store it */
      PyDict_SetItem(wobj->dict, key, val);
    }
  }
  return 0; /* 0 success, -1 failure */
}

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

void PLockStatus(PyMOLGlobals * G)
{                               /* assumes we have the GIL */
  PXDecRef(PYOBJECT_CALLFUNCTION(G->P_inst->lock_status, "O", G->P_inst->cmd));
}

int PLockStatusAttempt(PyMOLGlobals * G)
{                               /* assumes we have the GIL */
  int result = true;
  PyObject *got_lock =
    PYOBJECT_CALLFUNCTION(G->P_inst->lock_status_attempt, "O", G->P_inst->cmd);
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
  PXDecRef(PYOBJECT_CALLFUNCTION(G->P_inst->unlock_status, "O", G->P_inst->cmd));
}

static void PLockGLUT(PyMOLGlobals * G)
{                               /* assumes we have the GIL */
  PXDecRef(PYOBJECT_CALLFUNCTION(G->P_inst->lock_glut, "O", G->P_inst->cmd));
}

static void PUnlockGLUT(PyMOLGlobals * G)
{                               /* assumes we have the GIL */
  PXDecRef(PYOBJECT_CALLFUNCTION(G->P_inst->unlock_glut, "O", G->P_inst->cmd));
}

unsigned int P_glut_thread_id = -1;


/* enables us to keep glut out if by chance it grabs the API
 * in the middle of a nested API based operation */

void PCatchInit(void);
void my_interrupt(int a);


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
    PRunStringModule(G, "import pymol.vfont\n");
    P_vfont = PyDict_GetItemString(P_pymol_dict, "vfont");
    Py_XINCREF(P_vfont);
  }
  if(!P_vfont) {
    PRINTFB(G, FB_Python, FB_Errors)
      " PyMOL-Error: can't find module 'vfont'" ENDFB(G);
  } else {
    result = PYOBJECT_CALLMETHOD(P_vfont, "get_font", "fii", size, face, style);
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
    result = PYOBJECT_CALLFUNCTION(G->P_inst->complete, "s", str);
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

int PTruthCallStr0(PyObject * object, const char *method)
{
  int result = false;
  PyObject *tmp;
  tmp = PYOBJECT_CALLMETHOD(object, method, "");
  if(tmp) {
    if(PyObject_IsTrue(tmp))
      result = 1;
    Py_DECREF(tmp);
  }
  return (result);
}

int PTruthCallStr(PyObject * object, const char *method, const char *argument)
{
  int result = false;
  PyObject *tmp;
  tmp = PYOBJECT_CALLMETHOD(object, method, "s", argument);
  if(tmp) {
    if(PyObject_IsTrue(tmp))
      result = 1;
    Py_DECREF(tmp);
  }
  return (result);
}

int PTruthCallStr1i(PyObject * object, const char *method, int argument)
{
  int result = false;
  PyObject *tmp;
  tmp = PYOBJECT_CALLMETHOD(object, method, "i", argument);
  if(tmp) {
    if(PyObject_IsTrue(tmp))
      result = 1;
    Py_DECREF(tmp);
  }
  return (result);
}

int PTruthCallStr1s(PyObject * object, const char *method, const char *argument)
{
  int result = false;
  PyObject *tmp;
  tmp = PYOBJECT_CALLMETHOD(object, method, "s", argument);
  if(tmp) {
    if(PyObject_IsTrue(tmp))
      result = 1;
    Py_DECREF(tmp);
  }
  return (result);
}

int PTruthCallStr4i(PyObject * object, const char *method, int a1, int a2, int a3, int a4)
{
  int result = false;
  PyObject *tmp;
  tmp = PYOBJECT_CALLMETHOD(object, method, "iiii", a1, a2, a3, a4);
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
    PXDecRef(PYOBJECT_CALLMETHOD(G->P_inst->cmd, "_cache_set",
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
      output = PYOBJECT_CALLMETHOD(G->P_inst->cmd, "_cache_get",
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
  PXDecRef(PYOBJECT_CALLFUNCTION(P_sleep, "f", usec / 1000000.0));
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
  PXDecRef(PYOBJECT_CALLFUNCTION(P_sleep, "f", usec / 1000000.0));
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
  PXDecRef(PYOBJECT_CALLFUNCTION(P_sleep, "f", usec / 1000000.0));
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
  PYOBJECT_CALLMETHOD(P_traceback, "print_tb", "O", err);
}

void PDumpException()
{
  PYOBJECT_CALLMETHOD(P_traceback, "print_exc", "");
}

int PAlterAtomState(PyMOLGlobals * G, float *v, PyCodeObject *expr_co, int read_only,
                    ObjectMolecule *obj, CoordSet *cs, AtomInfoType * at, const char *model, int index, int csindex, int state, PyObject * space)

/* assumes Blocked python interpreter */
{
  int result = true;

  G->P_inst->wrapperObject->obj = obj;
  G->P_inst->wrapperObject->cs = cs;
  G->P_inst->wrapperObject->atomInfo = at;
  G->P_inst->wrapperObject->model = model;
  G->P_inst->wrapperObject->index = index;
  G->P_inst->wrapperObject->csindex = csindex;
  G->P_inst->wrapperObject->read_only = read_only;
  G->P_inst->wrapperObject->v = v;
  G->P_inst->wrapperObject->state = state + 1;

  PXDecRef(PyEval_EvalCode(expr_co, space, (PyObject*)G->P_inst->wrapperObject));
  WrapperObjectReset(G->P_inst->wrapperObject);

  if(PyErr_Occurred()) {
    PyErr_Print();
    result = false;
  }
  return result;
}

int PAlterAtom(PyMOLGlobals * G,
               ObjectMolecule *obj, CoordSet *cs, AtomInfoType * at, PyCodeObject *expr_co, int read_only,
               const char *model, int index, PyObject * space)
{
  int result = true;

  G->P_inst->wrapperObject->obj = obj;
  G->P_inst->wrapperObject->cs = cs;
  G->P_inst->wrapperObject->atomInfo = at;
  G->P_inst->wrapperObject->model = model;
  G->P_inst->wrapperObject->index = index;
  G->P_inst->wrapperObject->read_only = read_only;
  G->P_inst->wrapperObject->v = NULL;
  G->P_inst->wrapperObject->state = -1;

  PXDecRef(PyEval_EvalCode(expr_co, space, (PyObject*)G->P_inst->wrapperObject));
  WrapperObjectReset(G->P_inst->wrapperObject);

  if(PyErr_Occurred()) {
    PyErr_Print();
    result = false;
  }
  return result;
}

int PLabelAtom(PyMOLGlobals * G, ObjectMolecule *obj, CoordSet *cs, AtomInfoType * at, PyCodeObject *expr_co, const char *model, int index)
{
  int result = true;
  PyObject *P_inst_dict = G->P_inst->dict;
  PyObject *resultPyObject;
  OrthoLineType label;

  G->P_inst->wrapperObject->obj = obj;
  G->P_inst->wrapperObject->cs = cs;
  G->P_inst->wrapperObject->atomInfo = at;
  G->P_inst->wrapperObject->model = model;
  G->P_inst->wrapperObject->index = index;
  G->P_inst->wrapperObject->read_only = true;
  G->P_inst->wrapperObject->v = NULL;
  G->P_inst->wrapperObject->state = -1;

  if (!expr_co){
    // unsetting label
    if(at->label) {
      OVLexicon_DecRef(G->Lexicon, at->label);
    }
    at->label = 0;
    return true;
  }
  resultPyObject = PyEval_EvalCode(expr_co, P_inst_dict, (PyObject*)G->P_inst->wrapperObject);
  WrapperObjectReset(G->P_inst->wrapperObject);

  if(PyErr_Occurred()) {
    PyErr_Print();
    result = false;
  } else {
    result = true;
    if(!PConvPyObjectToStrMaxLen(resultPyObject,
                                 label, sizeof(OrthoLineType) - 1))
      result = false;
    if(PyErr_Occurred()) {
      PyErr_Print();
      result = false;
    }
    if(result) {
      LexDec(G, at->label);
      at->label = LexIdx(G, label);
    } else {
      ErrMessage(G, "Label", "Aborting on error. Labels may be incomplete.");
    }
  }
  PXDecRef(resultPyObject);
  return (result);
}

void PUnlockAPIAsGlut(PyMOLGlobals * G)
{                               /* must call with unblocked interpreter */
  PRINTFD(G, FB_Threads)
    " PUnlockAPIAsGlut-DEBUG: entered as thread %ld\n", PyThread_get_thread_ident()
    ENDFD;
  PBlock(G);
  PXDecRef(PYOBJECT_CALLFUNCTION(G->P_inst->unlock, "iO", 0, G->P_inst->cmd));  /* NOTE this may flush the command buffer! */
  PLockStatus(G);
  PyMOL_PopValidContext(G->PyMOL);
  PUnlockStatus(G);
  PUnlockGLUT(G);
  PUnblock(G);
}

void PUnlockAPIAsGlutNoFlush(PyMOLGlobals * G)
{                               /* must call with unblocked interpreter */
  PRINTFD(G, FB_Threads)
    " PUnlockAPIAsGlut-DEBUG: entered as thread %ld\n", PyThread_get_thread_ident()
    ENDFD;
  PBlock(G);
  PXDecRef(PYOBJECT_CALLFUNCTION(G->P_inst->unlock, "iO", -1, G->P_inst->cmd)); /* prevents flushing of the buffer */
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

    PXDecRef(PYOBJECT_CALLFUNCTION(G->P_inst->lock, "O", G->P_inst->cmd));

  } else {                      /* not blocking if PyMOL is busy */

    PyObject *got_lock =
      PYOBJECT_CALLFUNCTION(G->P_inst->lock_attempt, "O", G->P_inst->cmd);

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
              PXDecRef(PYOBJECT_CALLFUNCTION(G->P_inst->lock, "O", G->P_inst->cmd));
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
    "*PLockAPIAsGlut-DEBUG: entered as thread %ld\n", PyThread_get_thread_ident()
    ENDFD;

  PBlock(G);

  PLockGLUT(G);

  PLockStatus(G);
  PyMOL_PushValidContext(G->PyMOL);
  PUnlockStatus(G);

  PRINTFD(G, FB_Threads)
    "#PLockAPIAsGlut-DEBUG: acquiring lock as thread %ld\n", PyThread_get_thread_ident()
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
      "-PLockAPIAsGlut-DEBUG: glut_thread_keep_out %ld\n", PyThread_get_thread_ident()
      ENDFD;

    PXDecRef(PYOBJECT_CALLFUNCTION(G->P_inst->unlock, "iO", -1, G->P_inst->cmd));       /* prevent buffer flushing */
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
    PXDecRef(PYOBJECT_CALLFUNCTION(P_sleep, "f", 0.050));
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

#ifdef __cplusplus
extern "C" {
#endif

/* 
 *  void        initExtensionClass(void);
 *   void        initsglite(void);
 */
void init_champ(void);

#ifdef __cplusplus
}
#endif
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
#ifdef __cplusplus
extern "C" {
#endif

/*
 * void        initExtensionClass(void);
 * void        initsglite(void);
 */
void init_champ(void);
#ifdef _PYMOL_PYOMM
void init_pyomm(void);
#endif

#ifdef __cplusplus
}
#endif
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
    if(hwndFound = FindWindowA(NULL, argv[0])) {
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
      DWORD lpcbData;
      DWORD lpType = REG_SZ;
      int r1, r2;
      char *pymol_path;
      char *pythonhome;
      int pythonhome_set = false;
      int restart_flag = false;

      pymol_path = getenv("PYMOL_PATH");
      pythonhome = getenv("PYTHONHOME");
      if((!pymol_path) || (!pythonhome)) {
        lpcbData = sizeof(OrthoLineType) - 1;
        r1 = RegOpenKeyExA(HKEY_CLASSES_ROOT,
#ifdef PYMOL_EVAL
			"Software\\PyMOL\\PyMOL Eval\\PYMOL_PATH",
#else
			"Software\\PyMOL\\PyMOL\\PYMOL_PATH",
#endif
                          0, KEY_EXECUTE, &phkResult);
        if(r1 != ERROR_SUCCESS) {
          r1 = RegOpenKeyExA(HKEY_CURRENT_USER,
#ifdef PYMOL_EVAL
			"Software\\PyMOL\\PyMOL Eval\\PYMOL_PATH",
#else
			"Software\\PyMOL\\PyMOL\\PYMOL_PATH",
#endif
                            0, KEY_EXECUTE, &phkResult);
        }
        if(r1 == ERROR_SUCCESS) {
          r2 = RegQueryValueExA(phkResult, "", NULL, &lpType, (LPBYTE) path_buffer, &lpcbData);
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
          STARTUPINFOA si;
          PROCESS_INFORMATION pi;
          HANDLE hProcess = GetCurrentProcess();

          ZeroMemory(&si, sizeof(STARTUPINFOA));
          si.cb = sizeof(STARTUPINFOA);
          si.dwFlags = STARTF_USESHOWWINDOW;
          si.wShowWindow = SW_HIDE;

          if(IsSecurityRequired()) {
            lpSD = GlobalAlloc(GPTR, SECURITY_DESCRIPTOR_MIN_LENGTH);
            InitializeSecurityDescriptor(lpSD, SECURITY_DESCRIPTOR_REVISION);
            SetSecurityDescriptorDacl(lpSD, -1, 0, 0);

            lpSA = (LPSECURITY_ATTRIBUTES) GlobalAlloc(GPTR, sizeof(SECURITY_ATTRIBUTES));
            lpSA->nLength = sizeof(SECURITY_ATTRIBUTES);
            lpSA->lpSecurityDescriptor = lpSD;
            lpSA->bInheritHandle = TRUE;
          }

          if(CreateProcessA(NULL, (LPSTR) cmd_line, lpSA, NULL, TRUE,
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

void PRunStringModule(PyMOLGlobals * G, const char *str)
{                               /* runs a string in the namespace of the pymol global module */
  PXDecRef(PYOBJECT_CALLFUNCTION(G->P_inst->exec, "Os", P_pymol, str));
}

void PRunStringInstance(PyMOLGlobals * G, const char *str)
{                               /* runs a string in the namespace of the pymol instance */
  PXDecRef(PYOBJECT_CALLFUNCTION(G->P_inst->exec, "Os", G->P_inst->obj, str));
}

void WrapperObjectReset(WrapperObject *wo){
  wo->obj = NULL;
  wo->cs = NULL;
  wo->atomInfo = NULL;
  wo->model = NULL;
  wo->v = NULL;
  PyDict_Clear(wo->dict);
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

    PRunStringModule(G, "import pymol.cmd\n");
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

    PRunStringModule(G, "import pymol.menu\n");
    P_menu = PyDict_GetItemString(P_pymol_dict, "menu");
    Py_XINCREF(P_menu);
    if(!P_menu)
      ErrFatal(G, "PyMOL", "can't find module 'menu'");

    PRunStringModule(G, "import pymol.setting\n");
    P_setting = PyDict_GetItemString(P_pymol_dict, "setting");
    Py_XINCREF(P_setting);
    if(!P_setting)
      ErrFatal(G, "PyMOL", "can't find module 'setting'");

    PRunStringModule(G, "import pymol.povray\n");
    P_povray = PyDict_GetItemString(P_pymol_dict, "povray");
    Py_XINCREF(P_povray);
    if(!P_povray)
      ErrFatal(G, "PyMOL", "can't find module 'povray'");

#ifdef _PYMOL_XRAY
    PRunStringModule(G, "import pymol.xray\n");
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

    PRunStringModule(G, "import pymol.parser\n");
    P_parser = PyDict_GetItemString(P_pymol_dict, "parser");
    Py_XINCREF(P_parser);
    if(!P_parser)
      ErrFatal(G, "PyMOL", "can't find module 'parser'");

    {
      PyObject *fn_closure = PyObject_GetAttrString(P_parser, "new_parse_closure");
      G->P_inst->parse = PYOBJECT_CALLFUNCTION(fn_closure, "O", G->P_inst->cmd);
      PXDecRef(fn_closure);
      if(!G->P_inst->parse)
        ErrFatal(G, "PyMOL", "can't create 'parse' function closure");
    }

    {
      PyObject *fn_closure = PyObject_GetAttrString(P_parser, "new_complete_closure");
      G->P_inst->complete = PYOBJECT_CALLFUNCTION(fn_closure, "O", G->P_inst->cmd);
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

    PRunStringModule(G, "import pymol.util\n");
    PRunStringModule(G, "import pymol.preset\n");
    PRunStringModule(G, "import pymol.contrib\n");

    PRunStringModule(G, "import string\n");

    /* backwards compatibility */

    PRunStringModule(G, "pm = cmd\n");
    PRunStringModule(G, "pmu = util\n");

    PRunStringModule(G, "glutThread = thread.get_ident()");

    P_glut_thread_id = (unsigned int)PyThread_get_thread_ident();

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

    Wrapper_Type.tp_name = "wrapper.Wrapper";
    Wrapper_Type.tp_basicsize = sizeof(WrapperObject);
    Wrapper_Type.tp_flags = Py_TPFLAGS_DEFAULT;
    wrapperMappingMethods.mp_length = NULL;
    wrapperMappingMethods.mp_subscript = &WrapperObjectSubScript;
    wrapperMappingMethods.mp_ass_subscript = &WrapperObjectAssignSubScript;
    Wrapper_Type.tp_as_mapping = &wrapperMappingMethods;
    
    if (PyType_Ready(&Wrapper_Type) < 0){
      PRINTFB(G, FB_Python, FB_Errors)
	" PInit: Wrapper_Type, settingWrapper_Type, propertyWrapper_Type not ready\n" ENDFB(G);
      return;
    }
    Py_INCREF(&Wrapper_Type);
    G->P_inst->wrapperObject = (WrapperObject *)PyType_GenericNew(&Wrapper_Type, Py_None, Py_None);
    G->P_inst->wrapperObject->G = G;
    G->P_inst->wrapperObject->dict = PyDict_New();
    Py_INCREF(G->P_inst->wrapperObject);
  }
}

int PPovrayRender(PyMOLGlobals * G, const char *header, const char *inp, const char *file, int width,
                  int height, int antialias)
{
  PyObject *result;
  int ok;
  PBlock(G);
  result =
    PYOBJECT_CALLMETHOD(P_povray, "render_from_string", "sssiii", header, inp, file,
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

  exit(code);

  Py_Exit(code);
}

void PParse(PyMOLGlobals * G, const char *str)
{
  OrthoCommandIn(G, str);
}

void PDo(PyMOLGlobals * G, const char *str)
{                               /* assumes we already hold the re-entrant API lock */
  int blocked;
  PyObject *ret ;
  blocked = PAutoBlock(G);
  ret = PYOBJECT_CALLFUNCTION(G->P_inst->cmd_do, "s", str);
  Py_XDECREF(ret);
  PAutoUnblock(G, blocked);
}

void PLog(PyMOLGlobals * G, const char *str, int format)

/* general log routine can write PML 

   or PYM commands to appropriate log file */
{
  int mode;
  int a = sizeof(OrthoLineType) - 15;
  int blocked;
  PyObject *log;
  OrthoLineType buffer = "";
  mode = SettingGetGlobal_i(G, cSetting_logging);
  if(mode) {
    blocked = PAutoBlock(G);
    log = PyDict_GetItemString(P_pymol_dict, P_log_file_str);
    if(log && (log != Py_None)) {
      if(format == cPLog_no_flush) {
        PYOBJECT_CALLMETHOD(log, "write", "s", str);    /* maximize responsiveness (for real-time) */
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
            while(a && str[a - 1] < 32) a--; /* trim CR/LF etc. */
          case cPLog_pml:
            strcpy(buffer, "cmd.do('''");
            strncat(buffer, str, a);
            strcat(buffer, "''')\n");
            break;
          case cPLog_pym:
            strcpy(buffer, str);
            strcat(buffer, "\n");
            break;
          }
        }
        PYOBJECT_CALLMETHOD(log, "write", "s", buffer);
        PYOBJECT_CALLMETHOD(log, "flush", "");
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
  mode = SettingGetGlobal_i(G, cSetting_logging);
  if(mode) {
    blocked = PAutoBlock(G);
    log = PyDict_GetItemString(P_pymol_dict, P_log_file_str);
    if(log && (log != Py_None)) {
      PYOBJECT_CALLMETHOD(log, "flush", "");
    }
    PAutoUnblock(G, blocked);
  }
}

int PFlush(PyMOLGlobals * G)
{
  /* NOTE: ASSUMES unblocked Python threads and a locked API */
  PyObject *err;
  int did_work = false;
  if(OrthoCommandWaiting(G)) {
    did_work = true;
    PBlock(G);
    if(!(PIsGlutThread() && G->P_inst->glut_thread_keep_out)) {
      /* don't run if we're currently banned */
      char *buffer = 0;
      int size;
      while((size = OrthoCommandOutSize(G))){
	if (!buffer){
	  buffer = VLACalloc(char, size);
	} else {
	  VLACheck(buffer, char, size);
	}
	OrthoCommandSetBusy(G, true);
	OrthoCommandOut(G, buffer);
        OrthoCommandNest(G, 1);
        PUnlockAPIWhileBlocked(G);
        if(PyErr_Occurred()) {
          PyErr_Print();
          PRINTFB(G, FB_Python, FB_Errors)
            " PFlush: Uncaught exception.  PyMOL may have a bug.\n" ENDFB(G);
        }
        PXDecRef(PYOBJECT_CALLFUNCTION(G->P_inst->parse, "si", buffer, 0));
        err = PyErr_Occurred();
        if(err) {
          PyErr_Print();
          PRINTFB(G, FB_Python, FB_Errors)
            " PFlush: Uncaught exception.  PyMOL may have a bug.\n" ENDFB(G);
        }
        PLockAPIWhileBlocked(G);
	OrthoCommandSetBusy(G, false);
        /* make sure no commands left at this level */
        while(OrthoCommandWaiting(G))
          PFlushFast(G);
        OrthoCommandNest(G, -1);
      }
      if (buffer)
	VLAFreeP(buffer);
    }
    PUnblock(G);
  }
  return did_work;
}

int PFlushFast(PyMOLGlobals * G)
{
  /* NOTE: ASSUMES we currently have blocked Python threads and an unlocked API */
  PyObject *err;
  int did_work = false;
  char *buffer = 0;
  int size;
  while((size = OrthoCommandOutSize(G))){
    if (!buffer){
      buffer = VLACalloc(char, size);
    } else {
      VLACheck(buffer, char, size);
    }
    OrthoCommandSetBusy(G, true);
    OrthoCommandOut(G, buffer);
    OrthoCommandNest(G, 1);
    did_work = true;
    PRINTFD(G, FB_Threads)
      " PFlushFast-DEBUG: executing '%s' as thread %ld\n", buffer,
      PyThread_get_thread_ident()
      ENDFD;
    if(PyErr_Occurred()) {
      PyErr_Print();
      PRINTFB(G, FB_Python, FB_Errors)
        " PFlushFast: Uncaught exception.  PyMOL may have a bug.\n" ENDFB(G);
    }
    PXDecRef(PYOBJECT_CALLFUNCTION(G->P_inst->parse, "si", buffer, 0));
    err = PyErr_Occurred();
    if(err) {
      PyErr_Print();
      PRINTFB(G, FB_Python, FB_Errors)
        " PFlushFast: Uncaught exception.  PyMOL may have a bug.\n" ENDFB(G);
    }
    OrthoCommandSetBusy(G, false);
    /* make sure no commands left at this level */
    while(OrthoCommandWaiting(G))
      PFlushFast(G);
    OrthoCommandNest(G, -1);
  }
  if (buffer)
    VLAFreeP(buffer);

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

      PXDecRef(PYOBJECT_CALLFUNCTION(G->P_inst->lock_c, "O", G->P_inst->cmd));
      SavedThread[a].id = -1;
      /* this is the only safe time we can change things */
      PXDecRef(PYOBJECT_CALLFUNCTION(G->P_inst->unlock_c, "O", G->P_inst->cmd));

      PRINTFD(G, FB_Threads)
        " PAutoBlock-DEBUG: blocked %ld (%d, %d, %d)\n",
        PyThread_get_thread_ident(), SavedThread[MAX_SAVED_THREAD - 1].id,
        SavedThread[MAX_SAVED_THREAD - 2].id, SavedThread[MAX_SAVED_THREAD - 3].id ENDFD;

      return 1;
    }
    a--;
  }
  PRINTFD(G, FB_Threads)
    " PAutoBlock-DEBUG: %ld not found, thus already blocked.\n",
    PyThread_get_thread_ident()
    ENDFD;
  return 0;
#else
  return 1;
#endif
}

int PIsGlutThread(void)
{
  return (((unsigned int)PyThread_get_thread_ident()) == P_glut_thread_id);
}

void PUnblock(PyMOLGlobals * G)
{
#ifndef _PYMOL_EMBEDDED
  int a;
  SavedThreadRec *SavedThread = G->P_inst->savedThread;
  /* NOTE: ASSUMES a locked API */
  PRINTFD(G, FB_Threads)
    " PUnblock-DEBUG: entered as thread %ld\n", PyThread_get_thread_ident()
    ENDFD;

  /* reserve a space while we have a lock */
  PXDecRef(PYOBJECT_CALLFUNCTION(G->P_inst->lock_c, "O", G->P_inst->cmd));
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
  PXDecRef(PYOBJECT_CALLFUNCTION(G->P_inst->unlock_c, "O", G->P_inst->cmd));
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
  PXDecRef(PYOBJECT_CALLFUNCTION(G->P_inst->unlock, "iO", 0, G->P_inst->cmd));
}

int PLockAPI(PyMOLGlobals * G, int block_if_busy)
{
  int result = true;
  PBlock(G);
  if(block_if_busy) {
    PXDecRef(PYOBJECT_CALLFUNCTION(G->P_inst->lock, "O", G->P_inst->cmd));
  } else {                      /* not blocking if PyMOL is busy */

    PyObject *got_lock =
      PYOBJECT_CALLFUNCTION(G->P_inst->lock_attempt, "O", G->P_inst->cmd);

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
  PXDecRef(PYOBJECT_CALLFUNCTION(G->P_inst->unlock, "iO", 0, G->P_inst->cmd));
  PUnblock(G);
}

static void PUnlockAPIWhileBlocked(PyMOLGlobals * G)
{
  PXDecRef(PYOBJECT_CALLFUNCTION(G->P_inst->unlock, "iO", -1, G->P_inst->cmd));
}

static void PLockAPIWhileBlocked(PyMOLGlobals * G)
{
  PXDecRef(PYOBJECT_CALLFUNCTION(G->P_inst->lock, "O", G->P_inst->cmd));
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
  PXDecRef(PYOBJECT_CALLFUNCTION(G->P_inst->lock, "O", G->P_inst->cmd));
  PUnblock(G);
}

void PDefineFloat(PyMOLGlobals * G, const char *name, float value)
{
  char buffer[OrthoLineLength];
  sprintf(buffer, "%s = %f\n", name, value);
  PBlock(G);
  PRunStringModule(G, buffer);
  PUnblock(G);
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
