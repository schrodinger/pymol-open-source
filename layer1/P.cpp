
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

#ifdef _PYMOL_MINGW
#define putenv _putenv
#endif

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
#include "Lex.h"
#include "Seeker.h"

#ifdef _PYMOL_IP_PROPERTIES
#include"Property.h"
#endif

#include <memory>

/**
 * Use with functions which return a PyObject - if they return NULL, then an
 * unexpected exception has occured.
 */
#define ASSERT_PYOBJECT_NOT_NULL(obj, return_value)                            \
  if (!obj) {                                                                  \
    PyErr_Print();                                                             \
    return return_value;                                                       \
  }

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
            label_len = label_copy_text(label, LexStr(G, at->name), label_len, label_max);
          } else if(!strcmp(tok, "resn")) {
            label_len = label_copy_text(label, LexStr(G, at->resn), label_len, label_max);
          } else if(!strcmp(tok, "resi")) {
            sprintf(buffer, "%d%c", at->resv, at->inscode);
          } else if(!strcmp(tok, "resv")) {
            sprintf(buffer, "%d", at->resv);
          } else if(!strcmp(tok, "chain")) {
            label_len = label_copy_text(label, LexStr(G, at->chain), label_len, label_max);
          } else if(!strcmp(tok, "alt")) {
            label_len = label_copy_text(label, at->alt, label_len, label_max);
          } else if(!strcmp(tok, "segi")) {
            label_len = label_copy_text(label, LexStr(G, at->segi), label_len, label_max);
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
            strcpy(buffer, AtomInfoGetStereoAsStr(at));
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

static PyObject *P_vfont = NULL;

/* module import helper */

static
PyObject * PImportModuleOrFatal(const char * name) {
  PyObject * mod = PyImport_ImportModule(name);
  if(!mod) {
    fprintf(stderr, "PyMOL-Error: can't find '%s'\n", name);
    exit(EXIT_FAILURE);
  }
  return mod;
}

static
PyObject * PGetAttrOrFatal(PyObject * o, const char * name) {
  PyObject * attr = PyObject_GetAttrString(o, name);
  if(!attr) {
    fprintf(stderr, "PyMOL-Error: can't find '%s'\n", name);
    exit(EXIT_FAILURE);
  }
  return attr;
}

/* used elsewhere */

PyObject *P_menu = NULL;        /* menu definitions are currently global */
PyObject *P_xray = NULL;        /* okay as global */
PyObject *P_chempy = NULL;      /* okay as global */
PyObject *P_models = NULL;      /* okay as global */
PyObject *P_setting = NULL;     /* okay as global -- just used for names */
PyObject *P_CmdException = nullptr;
PyObject *P_QuietException = nullptr;
PyObject *P_IncentiveOnlyException = nullptr;

static PyMappingMethods wrapperMappingMethods, settingMappingMethods;
static PyTypeObject Wrapper_Type = {
  PyVarObject_HEAD_INIT(NULL, 0)
  "wrapper.Wrapper",            /* tp_name */
  0,                            /* tp_basicsize */
};
static PyTypeObject settingWrapper_Type = {
  PyVarObject_HEAD_INIT(NULL, 0)
  "wrapper.SettingWrapper",     /* tp_name */
  0,                            /* tp_basicsize */
};

#ifdef _PYMOL_IP_PROPERTIES
#endif

/**
 * If `wob` is not in a valid state (outside iterate-family context), raise
 * an error and return false.
 */
static bool check_wrapper_scope(WrapperObject * wobj) {
  if (wobj && wobj->obj)
    return true;

  PyErr_SetString(PyExc_RuntimeError,
      "wrappers cannot be used outside the iterate-family commands");

  return false;
}

/**
 * key: Python int (setting index) or str (setting name)
 *
 * Return the setting index or -1 for unknown `key`
 *
 * Raise LookupError if `key` doesn't name a known setting
 */
static int get_and_check_setting_index(PyMOLGlobals * G, PyObject * key) {
  int setting_id;

  if(PyInt_Check(key)) {
    setting_id = PyInt_AS_LONG(key);
  } else {
    key = PyObject_Str(key);
    setting_id = SettingGetIndex(G, PyString_AS_STRING(key));
    Py_DECREF(key);
  }

  if (setting_id < 0 || setting_id >= cSetting_INIT) {
    PyErr_SetString(PyExc_LookupError, "unknown setting");
    return -1;
  }

  return setting_id;
}

/**
 * Access a setting with iterate et. al.
 *
 * s[key]
 *
 * obj: `s` object in iterate-family namespace
 *
 * Raise LookupError if `key` doesn't name a known setting
 */
static
PyObject *SettingWrapperObjectSubScript(PyObject *obj, PyObject *key){
  auto& wobj = reinterpret_cast<SettingPropertyWrapperObject*>(obj)->wobj;
  int setting_id;
  PyObject *ret = NULL;

  if (!check_wrapper_scope(wobj)) {
    return NULL;
  }

  auto G = wobj->G;

  if ((setting_id = get_and_check_setting_index(G, key)) == -1) {
    return NULL;
  }

  if (wobj->idx >= 0){
    // atom-state level
    ret = SettingGetIfDefinedPyObject(G, wobj->cs, wobj->idx, setting_id);
  }

  if (!ret){
    // atom level
    ret = SettingGetIfDefinedPyObject(G, wobj->atomInfo, setting_id);

    if (!ret) {
      // object-state, object, or global
      ret = SettingGetPyObject(G,
          wobj->cs ? wobj->cs->Setting.get() : NULL,
          wobj->obj->Setting.get(), setting_id);
    }
  }
  return PConvAutoNone(ret);
}

/**
 * Set an atom or atom-state level setting with alter or alter_state.
 *
 * s[key] = val
 *
 * obj: `s` object in cmd.alter/cmd.alter_state namespace
 *
 * Return 0 on success or -1 on failure.
 *
 * Raise TypeError if setting not modifiable in the current context, and
 * LookupError if `key` doesn't name a known setting
 */
static
int SettingWrapperObjectAssignSubScript(PyObject *obj, PyObject *key, PyObject *val){
  auto& wobj = reinterpret_cast<SettingPropertyWrapperObject*>(obj)->wobj;

  if (!check_wrapper_scope(wobj)) {
    return -1;
  }

  int setting_id;
  auto G = wobj->G;

  if (wobj->read_only){
    PyErr_SetString(PyExc_TypeError, "Use alter/alter_state to modify settings");
    return -1;
  }

  if ((setting_id = get_and_check_setting_index(G, key)) == -1) {
    return -1;
  }

  if (wobj->idx >= 0) {
    // atom-state level
    if(!SettingLevelCheck(G, setting_id, cSettingLevel_astate)) {
      PyErr_SetString(PyExc_TypeError,
          "only atom-state level settings can be set in alter_state function");
      return -1; // failure
    } else if (CoordSetSetSettingFromPyObject(G, wobj->cs, wobj->idx, setting_id, val)) {
    }
  } else {
    // atom level
    if(!SettingLevelCheck(G, setting_id, cSettingLevel_atom)) {
      PyErr_SetString(PyExc_TypeError,
          "only atom-level settings can be set in alter function");
      return -1; // failure
    } else if (AtomInfoSetSettingFromPyObject(G, wobj->atomInfo, setting_id, val)) {
      AtomInfoSettingGenerateSideEffects(G, wobj->obj, setting_id, wobj->atm);
    }
  }

  return 0; // success
}

#ifdef _PYMOL_IP_PROPERTIES
#endif

/**
 * Python iterator over atom or atom-state setting indices
 */
static PyObject* SettingWrapperObjectIter(PyObject *self)
{
  auto& wobj = reinterpret_cast<SettingPropertyWrapperObject*>(self)->wobj;

  if (!check_wrapper_scope(wobj)) {
    return NULL;
  }

  int unique_id = wobj->atomInfo->unique_id;

  if (wobj->idx >= 0) {
    unique_id =
      wobj->cs->atom_state_setting_id ?
      wobj->cs->atom_state_setting_id[wobj->idx] : 0;
  }

  PyObject * items = SettingUniqueGetIndicesAsPyList(wobj->G, unique_id);
  PyObject * iter = PyObject_GetIter(items);
  Py_XDECREF(items);

  return iter;
}

#ifdef _PYMOL_IP_PROPERTIES
#endif

/**
 * Allows attribute-like syntax for item lookups
 *
 * o.key -> o[key] if `key` is not an attribute of `o`
 */
static PyObject* PyObject_GenericGetAttrOrItem(PyObject *o, PyObject *key) {
  PyObject *ret = PyObject_GenericGetAttr(o, key);
  if (!PyErr_Occurred())
    return ret;
  PyErr_Clear();
  return PyObject_GetItem(o, key);
}

/**
 * Allows attribute-like syntax for item assignment
 *
 * `o.key = value` -> `o[key] = value`
 */
static
int PyObject_GenericSetAttrAsItem(PyObject *o, PyObject *key, PyObject *value) {
  return PyObject_SetItem(o, key, value);
}

/**
 * Generic getter for member variable pointer at struct byte offset
 */
template <typename T, typename S>
static T* get_member_pointer(S* instance, size_t offset)
{
  return reinterpret_cast<T*>(reinterpret_cast<char*>(instance) + offset);
}
template <typename T, typename S>
static T const* get_member_pointer(S const* instance, size_t offset)
{
  return reinterpret_cast<T const*>(
      reinterpret_cast<char const*>(instance) + offset);
}

/**
 * iterate-family namespace implementation: lookup
 *
 * Raise NameError if state attributes are accessed outside of iterate_state
 */
static
PyObject * WrapperObjectSubScript(PyObject *obj, PyObject *key){

  static PyObject * pystr_HETATM        = PyString_InternFromString("HETATM");
  static PyObject * pystr_ATOM          = PyString_InternFromString("ATOM");
  static PyObject * pystr_QuestionMark  = PyString_InternFromString("?");

  auto wobj = static_cast<WrapperObject*>(obj);

  if (!check_wrapper_scope(wobj))
    return NULL;

  PyMOLGlobals* G = wobj->G;
  PyObject* ret = nullptr;

  auto const keyobj = unique_PyObject_ptr(PyObject_Str(key));
  auto const aprop = PyString_AS_STRING(keyobj.get());
  auto const ap = PyMOL_GetAtomPropertyInfo(G->PyMOL, aprop);

  if (ap) {
#ifdef _PYMOL_IP_EXTRAS
    switch (ap->id) {
    case ATOM_PROP_STEREO:
      if (ObjectMoleculeUpdateMMStereoInfoForState(G, wobj->obj, wobj->state - 1) < 0) {
        PyErr_SetString(P_CmdException,
            "please install rdkit or set SCHRODINGER variable");
        return NULL;
      }
      break;
    case ATOM_PROP_TEXT_TYPE:
#ifndef NO_MMLIBS
      ObjectMoleculeUpdateAtomTypeInfoForState(G, wobj->obj, wobj->state - 1, 1, 0);
#endif
      break;
    }
#endif

    switch (ap->Ptype){
    case cPType_string:
      ret = PyUnicode_FromString(
          get_member_pointer<char>(wobj->atomInfo, ap->offset));
      break;
    case cPType_schar:
      ret = PyLong_FromLong(
          *get_member_pointer<signed char>(wobj->atomInfo, ap->offset));
      break;
    case cPType_int:
      ret =
          PyLong_FromLong(*get_member_pointer<int>(wobj->atomInfo, ap->offset));
      break;
    case cPType_uint32:
      ret = PyLong_FromUnsignedLong(
          *get_member_pointer<uint32_t>(wobj->atomInfo, ap->offset));
      break;
    case cPType_int_as_string:
      ret = PyUnicode_FromString(LexStr(wobj->G,
          *get_member_pointer<lexborrow_t>(wobj->atomInfo, ap->offset)));
      break;
    case cPType_float:
      ret = PyFloat_FromDouble(
          *get_member_pointer<float>(wobj->atomInfo, ap->offset));
      break;
    case cPType_char_as_type:
      ret = PIncRef(wobj->atomInfo->hetatm ? pystr_HETATM : pystr_ATOM);
      break;
    case cPType_model:
      ret = PyUnicode_FromString(wobj->obj->Name);
      break;
    case cPType_index:
      ret = PyLong_FromLong(wobj->atm + 1);
      break;
    case cPType_int_custom_type: {
      auto val = *get_member_pointer<int>(wobj->atomInfo, ap->offset);
      if (val != cAtomInfoNoType) {
        ret = PyLong_FromLong(val);
      } else {
        ret = PIncRef(pystr_QuestionMark);
      }
    } break;
    case cPType_xyz_float:
      if (wobj->idx < 0) {
        PyErr_SetString(PyExc_NameError,
            "x/y/z only available in iterate_state and alter_state");
      } else {
        ret = PyFloat_FromDouble(wobj->cs->coordPtr(wobj->idx)[ap->offset]);
      }
      break;
    case cPType_settings:
      if (!wobj->settingWrapperObject) {
        wobj->settingWrapperObject = static_cast<SettingPropertyWrapperObject*>(
            PyType_GenericNew(&settingWrapper_Type, Py_None, Py_None));
        wobj->settingWrapperObject->wobj = wobj;
      }
      ret = PIncRef(wobj->settingWrapperObject);
      break;
    case cPType_properties:
#ifndef _PYMOL_IP_PROPERTIES
      PyErr_SetString(P_IncentiveOnlyException,
          "'properties/p' not supported in Open-Source PyMOL");
#else
      static_assert(false, "");
#endif
      break;
    case cPType_state:
      ret = PyLong_FromLong(wobj->state);
      break;
    default:
      switch (ap->id) {
      case ATOM_PROP_RESI: {
        char resi[8];
        AtomResiFromResv(resi, sizeof(resi), wobj->atomInfo);
        ret = PyUnicode_FromString(resi);
      } break;
      case ATOM_PROP_STEREO: {
        auto mmstereotype = AtomInfoGetStereoAsStr(wobj->atomInfo);
        ret = PyUnicode_FromString(mmstereotype);
      } break;
      case ATOM_PROP_ONELETTER: {
        const char* st = LexStr(G, wobj->atomInfo->resn);
        char abbr[2] = {SeekerGetAbbr(G, st, 'O', 'X'), 0};
        ret = PyUnicode_FromString(abbr);
      } break;
      default:
        PyErr_SetString(PyExc_SystemError, "unhandled atom property type");
      }
    }
  } else {
    /* if not an atom property, check if local variable in dict */
    if (wobj->dict) {
      ret = PyDict_GetItem(wobj->dict, key); // Borrowed reference
    }
    if (ret) {
      Py_INCREF(ret);
    } else {
      PyErr_SetNone(PyExc_KeyError);
    }
  }

  return ret;
}

/**
 * Make IPython happy, which may call f_locals.get("__tracebackhide__", 0) and
 * crash if there is no get() method (the wrapper object is f_locals).
 */
static PyObject* WrapperObject_get(PyObject* self, PyObject* args)
{
  auto nargs = PyTuple_Size(args);
  assert(0 < nargs && nargs < 3);

  // Could call WrapperObjectSubScript here, but we don't really need that. To
  // fix the IPython issue, it's sufficient to return default or None.

  if (nargs == 2) {
    return PIncRef(PyTuple_GET_ITEM(args, 1));
  }

  Py_RETURN_NONE;
}

static PyMethodDef wrapperMethods[] = {
    {"get", WrapperObject_get, METH_VARARGS, nullptr},
    {nullptr},
};

/**
 * iterate-family namespace implementation: assignment
 *
 * Raise TypeError for read-only variables
 */
static
int WrapperObjectAssignSubScript(PyObject *obj, PyObject *key, PyObject *val){
  auto wobj = static_cast<WrapperObject*>(obj);

  if (!check_wrapper_scope(wobj)) {
    return -1;
  }

  auto G = wobj->G;

  const auto keyobj = unique_PyObject_ptr(PyObject_Str(key));
  const char* const aprop = PyString_AS_STRING(keyobj.get());
  const AtomPropertyInfo* ap = PyMOL_GetAtomPropertyInfo(G->PyMOL, aprop);

  if (ap) {
    if (wobj->read_only) {
      PyErr_SetString(
          PyExc_TypeError, "Use alter/alter_state to modify values");
      return -1;
    }

#ifdef _PYMOL_IP_EXTRAS
    if (wobj->cs) {
      switch (ap->id) {
      case ATOM_PROP_STEREO:
        wobj->cs->validMMStereo = MMPYMOLX_PROP_STATE_USER;
        break;
      case ATOM_PROP_TEXT_TYPE:
        wobj->cs->validTextType = MMPYMOLX_PROP_STATE_USER;
        break;
      }
    }
#endif

    bool changed = false;

    switch (ap->Ptype) {
    case cPType_string: {
      PyObject* valobj = PyObject_Str(val);
      const char* valstr = PyString_AS_STRING(valobj);
      char* dest = get_member_pointer<char>(wobj->atomInfo, ap->offset);
      if (strlen(valstr) > ap->maxlen) {
        strncpy(dest, valstr, ap->maxlen);
      } else {
        strcpy(dest, valstr);
      }
      Py_DECREF(valobj);
      changed = true;
    } break;
    case cPType_schar: {
      int valint = PyInt_AsLong(val);
      if (valint == -1 && PyErr_Occurred())
        return -1;
      *get_member_pointer<signed char>(wobj->atomInfo, ap->offset) = valint;
      changed = true;
    } break;
    case cPType_int: {
      int valint = PyInt_AsLong(val);
      if (valint == -1 && PyErr_Occurred())
        return -1;
      *get_member_pointer<int>(wobj->atomInfo, ap->offset) = valint;
      changed = true;
    } break;
    case cPType_uint32: {
      auto valint = PyLong_AsUnsignedLong(val);
      if (valint == -1 && PyErr_Occurred())
        return -1;
      *get_member_pointer<uint32_t>(wobj->atomInfo, ap->offset) = valint;
      changed = true;
    } break;
    case cPType_int_as_string: {
      auto dest = get_member_pointer<lexidx_t>(wobj->atomInfo, ap->offset);
      const auto valobj = unique_PyObject_ptr(PyObject_Str(val));
      const char* valstr = PyString_AS_STRING(valobj.get());
      LexAssign(G, *dest, valstr);
      changed = true;
    } break;
    case cPType_float:
      if (!PConvPyObjectToFloat(
              val, get_member_pointer<float>(wobj->atomInfo, ap->offset))) {
        return -1;
      }
      changed = true;
      break;
    case cPType_char_as_type: {
      const auto valobj = unique_PyObject_ptr(PyObject_Str(val));
      const char* valstr = PyString_AS_STRING(valobj.get());
      wobj->atomInfo->hetatm = ((valstr[0] == 'h') || (valstr[0] == 'H'));
      changed = true;
    } break;
    case cPType_int_custom_type: {
      const auto valobj = unique_PyObject_ptr(PyObject_Str(val));
      const char* valstr = PyString_AS_STRING(valobj.get());
      auto* dest = get_member_pointer<int>(wobj->atomInfo, ap->offset);
      if (valstr[0] == '?') {
        *dest = cAtomInfoNoType;
      } else {
        int valint = PyInt_AS_LONG(val);
        *dest = valint;
      }
      changed = true;
    } break;
    case cPType_xyz_float:
      if (wobj->idx < 0) {
        PyErr_SetString(PyExc_NameError, "x/y/z only available in alter_state");
        return -1;
      } else {
        float* v = wobj->cs->coordPtr(wobj->idx) + ap->offset;
        if (!PConvPyObjectToFloat(val, v)) {
          return -1;
        }
      }
      break;
    default:
      switch (ap->id) {
      case ATOM_PROP_RESI:
        if (PConvPyIntToInt(val, &wobj->atomInfo->resv)) {
          wobj->atomInfo->inscode = '\0';
        } else {
          const auto valobj = unique_PyObject_ptr(PyObject_Str(val));
          wobj->atomInfo->setResi(PyString_AS_STRING(valobj.get()));
        }
        break;
      case ATOM_PROP_STEREO: {
        const auto valobj = unique_PyObject_ptr(PyObject_Str(val));
        const char* valstr = PyString_AS_STRING(valobj.get());
        AtomInfoSetStereo(wobj->atomInfo, valstr);
      } break;
      default:
        PyErr_Format(PyExc_TypeError, "'%s' is read-only", aprop);
        return -1;
      }
    }

    if (changed) {
      switch (ap->id) {
      case ATOM_PROP_ELEM:
        wobj->atomInfo->protons = 0;
        wobj->atomInfo->vdw = 0;
        AtomInfoAssignParameters(wobj->G, wobj->atomInfo);
        break;
      case ATOM_PROP_RESV:
        wobj->atomInfo->inscode = '\0';
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
    if (!wobj->dict) {
      wobj->dict = PyDict_New();
    }
    PyDict_SetItem(wobj->dict, key, val);
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

/**
 * Acquire the status lock (blocking)
 * @pre GIL
 */
void PLockStatus(PyMOLGlobals * G)
{                               /* assumes we have the GIL */
  assert(PyGILState_Check());

  PXDecRef(PYOBJECT_CALLMETHOD(G->P_inst->lock_api_status, "acquire", nullptr));
}

/**
 * Acquire the status lock (non-blocking)
 * @return True if lock could be acquired
 * @pre GIL
 */
int PLockStatusAttempt(PyMOLGlobals * G)
{                               /* assumes we have the GIL */
  assert(PyGILState_Check());

  unique_PyObject_ptr got_lock(
      PYOBJECT_CALLMETHOD(G->P_inst->lock_api_status, "acquire", "i", 0));

  ASSERT_PYOBJECT_NOT_NULL(got_lock, true);

  return PyObject_IsTrue(got_lock.get());
}

/**
 * Release the status lock
 * @pre GIL
 */
void PUnlockStatus(PyMOLGlobals * G)
{                               /* assumes we have the GIL */
  assert(PyGILState_Check());

  PXDecRef(PYOBJECT_CALLMETHOD(G->P_inst->lock_api_status, "release", nullptr));
}

/**
 * Acquire GLUT lock
 * @pre GIL
 */
static void PLockGLUT(PyMOLGlobals * G)
{                               /* assumes we have the GIL */
  PXDecRef(PYOBJECT_CALLMETHOD(G->P_inst->lock_api_glut, "acquire", nullptr));
}

/**
 * Release GLUT lock
 * @pre GIL
 */
static void PUnlockGLUT(PyMOLGlobals * G)
{                               /* assumes we have the GIL */
  PXDecRef(PYOBJECT_CALLMETHOD(G->P_inst->lock_api_glut, "release", nullptr));
}

static long P_glut_thread_id = -1;


/* enables us to keep glut out if by chance it grabs the API
 * in the middle of a nested API based operation */

static
void PCatchInit(void);


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
  assert(PyGILState_Check());

  PyObject *result = NULL;

  if(!P_vfont) {
    P_vfont = PyImport_ImportModule("pymol.vfont");
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
  assert(!PyGILState_Check());

  int ret = false;
  PyObject *result;
  const char *st2;
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
  assert(PyGILState_Check());

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
  assert(PyGILState_Check());

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
  assert(PyGILState_Check());

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
  assert(PyGILState_Check());

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
  assert(PyGILState_Check());

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
  assert(PyGILState_Check());

  if(!obj)
    obj = Py_None;
  Py_XINCREF(obj);
  return obj;
}

void PXDecRef(PyObject * obj)
{
  assert(PyGILState_Check());

  Py_XDECREF(obj);
}

static ov_status CacheCreateEntry(PyObject ** result, PyObject * input)
{
  assert(PyGILState_Check());

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
  assert(PyGILState_Check());

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
  assert(PyGILState_Check());

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
  assert(!PyGILState_Check());

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
  assert(!PyGILState_Check());

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
  assert(!PyGILState_Check());

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
static PyObject *PCatch_install(PyObject * self, PyObject * args);

static
void my_interrupt(int a)
{
  exit(EXIT_FAILURE);
}

void PDumpTraceback(PyObject * err)
{
  assert(PyGILState_Check());

  PYOBJECT_CALLMETHOD(P_traceback, "print_tb", "O", err);
}

void PDumpException()
{
  assert(PyGILState_Check());

  PYOBJECT_CALLMETHOD(P_traceback, "print_exc", "");
}

static
WrapperObject * WrapperObjectNew() {
  auto wobj = static_cast<WrapperObject*>(
      PyType_GenericNew(&Wrapper_Type, Py_None, Py_None));
  wobj->dict = NULL;
  wobj->settingWrapperObject = NULL;
#ifdef _PYMOL_IP_PROPERTIES
  wobj->propertyWrapperObject = NULL;
#endif
  return wobj;
}

/**
 * @pre GIL
 */
int PAlterAtomState(PyMOLGlobals * G, PyObject *expr_co, int read_only,
                    ObjectMolecule *obj, CoordSet *cs, int atm, int idx,
                    int state, PyObject * space)
{
  assert(PyGILState_Check());

  int result = true;

  auto wobj = WrapperObjectNew();
  wobj->G = G;
  wobj->obj = obj;
  wobj->cs = cs;
  wobj->atomInfo = obj->AtomInfo + atm;
  wobj->atm = atm;
  wobj->idx = idx;
  wobj->read_only = read_only;
  wobj->state = state + 1;

  PXDecRef(PyEval_EvalCode((PyObject*) expr_co, space, (PyObject*) wobj));
  Py_DECREF(wobj);

  if(PyErr_Occurred()) {
    result = false;
  }
  return result;
}

int PAlterAtom(PyMOLGlobals* G, ObjectMolecule* obj, CoordSet* cs,
    PyObject* expr_co, int read_only, int atm, PyObject* space)
{
  int state = (obj->DiscreteFlag ? obj->AtomInfo[atm].discrete_state : 0) - 1;
  return PAlterAtomState(G, expr_co, read_only, obj, cs, atm, /* idx */ -1, state, space);
}

/**
 * String conversion which takes "label_digits" setting into account.
 */
static
int PLabelPyObjectToStrMaxLen(PyMOLGlobals * G, PyObject * obj, char *buffer, int maxlen)
{
  assert(PyGILState_Check());

  if (obj && PyFloat_Check(obj)) {
    snprintf(buffer, maxlen + 1, "%.*f",
        SettingGetGlobal_i(G, cSetting_label_digits),
        PyFloat_AsDouble(obj));
    return true;
  }
  return PConvPyObjectToStrMaxLen(obj, buffer, maxlen);
}

int PLabelAtom(PyMOLGlobals* G, ObjectMolecule* obj, CoordSet* cs,
    PyObject* expr_co, int atm)
{
  assert(PyGILState_Check());

  PyObject *P_inst_dict = G->P_inst->dict;
  OrthoLineType label;
  AtomInfoType * ai = obj->AtomInfo + atm;

  if (!expr_co){
    // unsetting label
    LexAssign(G, ai->label, 0);
    return true;
  }

  auto wobj = WrapperObjectNew();
  wobj->G = G;
  wobj->obj = obj;
  wobj->cs = cs;
  wobj->atomInfo = ai;
  wobj->atm = atm;
  wobj->idx = -1;
  wobj->read_only = true;

  if (obj->DiscreteFlag) {
    wobj->state = obj->AtomInfo[atm].discrete_state;
  } else {
    wobj->state = 0;
  }

  auto resultPyObject =
      unique_PyObject_ptr(PyEval_EvalCode(expr_co, P_inst_dict, wobj));

  if (PyErr_Occurred()) {
    return false;
  }

  if (!PLabelPyObjectToStrMaxLen(
          G, resultPyObject.get(), label, sizeof(OrthoLineType) - 1)) {
    if (!PyErr_Occurred()) {
      ErrMessage(G, "Label", "Aborting on error. Labels may be incomplete.");
    }

    return false;
  }

  LexAssign(G, ai->label, label);
  return true;
}

/**
 * - Release API lock with flushing
 * - Pop valid context
 * - Release GLUT lock
 * @pre no GIL
 */
void PUnlockAPIAsGlut(PyMOLGlobals * G)
{                               /* must call with unblocked interpreter */
  assert(!PyGILState_Check());

  PBlock(G);
  PXDecRef(PYOBJECT_CALLFUNCTION(G->P_inst->unlock, "iO", 0, G->P_inst->cmd));  /* NOTE this may flush the command buffer! */
  PLockStatus(G);
  PyMOL_PopValidContext(G->PyMOL);
  PUnlockStatus(G);
  PUnlockGLUT(G);
  PUnblock(G);
}

/**
 * - Release API lock without flushing
 * - Pop valid context
 * - Release GLUT lock
 * @pre no GIL
 */
void PUnlockAPIAsGlutNoFlush(PyMOLGlobals * G)
{                               /* must call with unblocked interpreter */
  assert(!PyGILState_Check());

  PBlock(G);
  PXDecRef(PYOBJECT_CALLFUNCTION(G->P_inst->unlock, "iO", -1, G->P_inst->cmd)); /* prevents flushing of the buffer */
  PLockStatus(G);
  PyMOL_PopValidContext(G->PyMOL);
  PUnlockStatus(G);
  PUnlockGLUT(G);
  PUnblock(G);
}

/**
 * Acquire API lock
 * @param block_if_busy If false and PyMOL is busy, do not block
 * @return False if API lock could not be acquired
 * @pre GIL
 */
static int get_api_lock(PyMOLGlobals * G, int block_if_busy)
{
  assert(PyGILState_Check());

  if (!block_if_busy) {
    unique_PyObject_ptr got_lock(
        PYOBJECT_CALLFUNCTION(G->P_inst->lock_attempt, "O", G->P_inst->cmd));

    ASSERT_PYOBJECT_NOT_NULL(got_lock, false);

    if (PyObject_IsTrue(got_lock.get())) {
      return true;
    }

    PLockStatus(G);
    bool busy = PyMOL_GetBusy(G->PyMOL, false);
    PUnlockStatus(G);

    if (busy) {
      return false;
    }
  }

  PXDecRef(PYOBJECT_CALLFUNCTION(G->P_inst->lock, "O", G->P_inst->cmd));
  return true;
}

/**
 * - Acquire GLUT lock
 * - Push valid context
 * - Acquire API lock
 * - Wait for glut_thread_keep_out == 0
 *
 * @param block_if_busy If false and PyMOL is busy, API lock is non-blocking
 * @return False if locks could not be acquired
 * @pre no GIL
 */
int PLockAPIAsGlut(PyMOLGlobals * G, int block_if_busy)
{
  assert(!PyGILState_Check());

  PBlock(G);

  PLockGLUT(G);

  PLockStatus(G);
  PyMOL_PushValidContext(G->PyMOL);
  PUnlockStatus(G);

  while (true) {
    if(!get_api_lock(G, block_if_busy)) {
      PLockStatus(G);
      PyMOL_PopValidContext(G->PyMOL);
      PUnlockStatus(G);
      PUnlockGLUT(G);
      PUnblock(G);
      return false;               /* busy -- so allow main to update busy status display (if any) */
    }

    if (G->P_inst->glut_thread_keep_out == 0) {
      break;
    }

    /* IMPORTANT: keeps the glut thread out of an API operation... */
    /* NOTE: the keep_out variable can only be changed or read by the thread
       holding the API lock, therefore it is safe even through increment
       isn't atomic. */

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
  }

  PUnblock(G);                  /* API is now locked, so we can free up Python... */

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


void PConvertOptions(CPyMOLOptions * rec, PyObject * options)
{
  assert(PyGILState_Check());

  const char *load_str;

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
  rec->launch_status = PyInt_AsLong(PyObject_GetAttrString(options, "launch_status"));
  rec->gldebug = PyInt_AsLong(PyObject_GetAttrString(options, "gldebug"));
  rec->openvr_stub = PyInt_AsLong(PyObject_GetAttrString(options, "openvr_stub"));

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
  assert(PyGILState_Check());

  PyObject *pymol, *invocation, *options;

  pymol = PImportModuleOrFatal("pymol");
  invocation = PGetAttrOrFatal(pymol, "invocation");     /* get a handle to the invocation module */
  options = PGetAttrOrFatal(invocation, "options");

  PConvertOptions(rec, options);
  Py_XDECREF(invocation);
  Py_XDECREF(options);
  Py_XDECREF(pymol);
}

/**
 * Runs a string in the namespace of the pymol global module
 * @pre GIL
 */
void PRunStringModule(PyMOLGlobals * G, const char *str)
{
  assert(PyGILState_Check());

  PXDecRef(PYOBJECT_CALLFUNCTION(G->P_inst->exec, "Os", P_pymol, str));
}

/**
 * Runs a string in the namespace of the pymol instance
 * @pre GIL
 */
void PRunStringInstance(PyMOLGlobals * G, const char *str)
{
  assert(PyGILState_Check());

  PXDecRef(PYOBJECT_CALLFUNCTION(G->P_inst->exec, "Os", G->P_inst->obj, str));
}

static void WrapperObjectDealloc(PyObject* self)
{
  auto wo = static_cast<WrapperObject*>(self);
  Py_XDECREF(wo->settingWrapperObject);
#ifdef _PYMOL_IP_PROPERTIES
  Py_XDECREF(wo->propertyWrapperObject);
#endif
  Py_XDECREF(wo->dict);
  self->ob_type->tp_free(self);
}

void PInit(PyMOLGlobals * G, int global_instance)
{
  assert(PyGILState_Check());

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

  if(true /* global_instance */) {
    PCatchInit();               /* setup standard-output catch routine */
  }

  P_pymol = PImportModuleOrFatal("pymol");
  P_pymol_dict = PyModule_GetDict(P_pymol);
  Py_XINCREF(P_pymol_dict);
  if(!P_pymol_dict)
    ErrFatal(G, "PyMOL", "can't find globals for 'pymol'");

  if(global_instance) {         /* if global singleton PyMOL... */
    G->P_inst = pymol::calloc<CP_inst>(1);
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
    G->P_inst->exec = PGetAttrOrFatal(P_pymol, "exec_str");

    if(global_instance) {
      PCatch_install(NULL, NULL);
    }

    P_traceback = PImportModuleOrFatal("traceback");
    P_cmd = PImportModuleOrFatal("pymol.cmd");

    if(global_instance) {
      assert(SingletonPyMOLGlobals);
      /* implies global singleton pymol, so set up the global handle */
      PyObject_SetAttrString(P_cmd, "_COb",
          PyCapsule_New(&SingletonPyMOLGlobals, nullptr, nullptr));

      /* cmd module is itself the api for the global PyMOL instance */
      G->P_inst->cmd = P_cmd;
    }

    /* right now, all locks are global -- eventually some of these may
       become instance-specific in order to improve concurrency */

    G->P_inst->lock = PGetAttrOrFatal(G->P_inst->cmd, "lock");
    G->P_inst->lock_attempt = PGetAttrOrFatal(G->P_inst->cmd, "lock_attempt");
    G->P_inst->unlock = PGetAttrOrFatal(G->P_inst->cmd, "unlock");
    G->P_inst->lock_api_status = PGetAttrOrFatal(G->P_inst->cmd, "lock_api_status");
    G->P_inst->lock_api_glut = PGetAttrOrFatal(G->P_inst->cmd, "lock_api_glut");

    /* 'do' command */

    G->P_inst->cmd_do = PGetAttrOrFatal(G->P_inst->cmd, "do");

    /* cache */
    G->P_inst->cache = PyObject_GetAttrString(G->P_inst->obj, "_cache");

    /* invariant stuff */

    P_menu = PImportModuleOrFatal("pymol.menu");
    P_setting = PImportModuleOrFatal("pymol.setting");
    P_povray = PImportModuleOrFatal("pymol.povray");

#ifdef _PYMOL_XRAY
    P_xray = PImportModuleOrFatal("pymol.xray");
#endif

    /* BEGIN PROPRIETARY CODE SEGMENT (see disclaimer in "os_proprietary.h") */
#ifdef WIN32
    P_time = PImportModuleOrFatal("time");
    P_sleep = PGetAttrOrFatal(P_time, "sleep");
#endif
    /* END PROPRIETARY CODE SEGMENT */

    P_parser = PImportModuleOrFatal("pymol.parser");

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

    {
      PyObject *fn_closure = PGetAttrOrFatal(P_pymol, "_colortype");
      G->P_inst->colortype = PYOBJECT_CALLFUNCTION(fn_closure, "O", G->P_inst->cmd);
      PXDecRef(fn_closure);
    }

    P_chempy = PImportModuleOrFatal("chempy");
    P_models = PImportModuleOrFatal("chempy.models");
    P_CmdException = PGetAttrOrFatal(P_pymol, "CmdException");
    P_QuietException = PGetAttrOrFatal(P_cmd, "QuietException");
    P_IncentiveOnlyException = PGetAttrOrFatal(P_pymol, "IncentiveOnlyException");

    /* backwards compatibility */

    PRunStringModule(G, "glutThread = thread.get_ident()");

    P_glut_thread_id = PyThread_get_thread_ident();

#ifndef WIN32
    if(G->Option->siginthand) {
      signal(SIGINT, my_interrupt);
    }
#endif

  if (!Wrapper_Type.tp_basicsize) {
    Wrapper_Type.tp_basicsize = sizeof(WrapperObject);
    Wrapper_Type.tp_dealloc = &WrapperObjectDealloc;
    Wrapper_Type.tp_flags = Py_TPFLAGS_DEFAULT;
    wrapperMappingMethods.mp_length = NULL;
    wrapperMappingMethods.mp_subscript = &WrapperObjectSubScript;
    wrapperMappingMethods.mp_ass_subscript = &WrapperObjectAssignSubScript;
    Wrapper_Type.tp_as_mapping = &wrapperMappingMethods;
    Wrapper_Type.tp_methods = wrapperMethods;
    
    settingWrapper_Type.tp_basicsize = sizeof(SettingPropertyWrapperObject);
    settingWrapper_Type.tp_flags = Py_TPFLAGS_DEFAULT;
    settingWrapper_Type.tp_iter = &SettingWrapperObjectIter;
    settingMappingMethods.mp_length = NULL;
    settingMappingMethods.mp_subscript = &SettingWrapperObjectSubScript;
    settingMappingMethods.mp_ass_subscript = &SettingWrapperObjectAssignSubScript;
    settingWrapper_Type.tp_as_mapping = &settingMappingMethods;
    settingWrapper_Type.tp_getattro = PyObject_GenericGetAttrOrItem;
    settingWrapper_Type.tp_setattro = PyObject_GenericSetAttrAsItem;
    
#ifdef _PYMOL_IP_PROPERTIES
#endif

    if (PyType_Ready(&Wrapper_Type) < 0
        || PyType_Ready(&settingWrapper_Type) < 0
        ){
      PRINTFB(G, FB_Python, FB_Errors)
	" PInit: Wrapper_Type, settingWrapper_Type, propertyWrapper_Type not ready\n" ENDFB(G);
      return;
    }
    Py_INCREF(&Wrapper_Type);
    Py_INCREF(&settingWrapper_Type);
  }
  }

#ifdef _PYMOL_NO_MSGPACKC
  // fallback MMTF support
  PyRun_SimpleString(
      "import pymol.importing;"
      "pymol.importing.loadfunctions.setdefault('mmtf',"
      "pymol.importing.load_mmtf)");
#endif
}

int PPovrayRender(PyMOLGlobals * G, const char *header, const char *inp, const char *file, int width,
                  int height, int antialias)
{
  assert(!PyGILState_Check());

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

void PFree(PyMOLGlobals * G)
{
  assert(PyGILState_Check());

  PXDecRef(G->P_inst->parse);
  PXDecRef(G->P_inst->complete);
  PXDecRef(G->P_inst->colortype);
}

/**
 * Terminates the application with a call to `exit(code)`.
 */
void PExit(PyMOLGlobals * G, int code)
{
  assert(!PyGILState_Check());

  ExecutiveDelete(G, "all");
  PBlock(G);

  PyMOL_PushValidContext(G->PyMOL);
  PyMOL_Stop(G->PyMOL);
  PyMOL_PopValidContext(G->PyMOL);

#ifndef _PYMOL_NO_MAIN
  if(G->Main) {
    MainFree();
  }
#endif

  PyMOL_Free(G->PyMOL);

#if 1
  /* we're having trouble with threading errors after calling Py_Exit,
     so for the time being, let's just take the process down at this
     point, instead of allowing PyExit to be called. */

  exit(code);
#else

  Py_Exit(code);
#endif
}

/**
 * Add `str` to the command queue
 */
void PParse(PyMOLGlobals * G, pymol::zstring_view str_view)
{
  OrthoCommandIn(G, str_view.c_str());
}

/**
 * Call `cmd.do(str)`
 *
 * Effectively calls `PParse(str.splitlines())` with logging and echo.
 */
void PDo(PyMOLGlobals * G, const char *str)
{                               /* assumes we already hold the re-entrant API lock */
  int blocked;
  PyObject *ret ;
  blocked = PAutoBlock(G);
  ret = PYOBJECT_CALLFUNCTION(G->P_inst->cmd_do, "s", str);
  Py_XDECREF(ret);
  PAutoUnblock(G, blocked);
}

/**
 * Write `str` to the log file (if one is open).
 *
 * str: command or expression to log
 * format: cPLog_pml (`str` is PyMOL command)
 *         cPLog_pym (`str` is Python expression)
 *         cPLog_no_flush (write `str` as is)
 *         cPLog_pml_lf (unused TODO remove?)
 *
 * See also equivalent Python impelemtation: cmd.log()
 */
void PLog(PyMOLGlobals * G, pymol::zstring_view str_view, int format)
{
  int mode;
  int a = sizeof(OrthoLineType) - 15;
  int blocked;
  PyObject *log;
  auto str = str_view.c_str();
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
            if (str[0] == '/') {
              strncat(buffer, str + 1, a - 1);
              strcat(buffer, "\n");
            } else {
              strcpy(buffer, "cmd.do('''");
              char * b = buffer + strlen(buffer);
              for (; a && *str; --a) {
                if (*str == '\\' || *str == '\'') {
                  *(b++) = '\\';
                }
                *(b++) = *(str++);
              }
              strcpy(b, "''')\n");
            }
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

/**
 * Call flush() on the open log file handle
 */
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

#define PFLUSH_REPORT_UNCAUGHT_EXCEPTIONS(G)                                   \
  if (PyErr_Occurred()) {                                                      \
    PyErr_Print();                                                             \
    PRINTFB(G, FB_Python, FB_Errors)                                           \
    " %s: Uncaught exception.  PyMOL may have a bug.\n", __func__ ENDFB(G);    \
  }

/**
 * Execute commands from the command queue.
 * If the current thread holds the GIL, do nothing. Otherwise the queue will
 * be empty after this call.
 * @return False if the queue was empty
 */
int PFlush(PyMOLGlobals * G)
{
  /* NOTE: ASSUMES unblocked Python threads and a locked API */
  if(!OrthoCommandWaiting(G))
    return false;

  if (PAutoBlock(G)) {
    if(!(PIsGlutThread() && G->P_inst->glut_thread_keep_out)) {
      /* don't run if we're currently banned */
      auto ortho = G->Ortho;
      while(!OrthoCommandIsEmpty(*ortho)){
        auto buffer = OrthoCommandOut(*ortho);
        OrthoCommandSetBusy(G, true);
        OrthoCommandNest(G, 1);

        // TODO if we release here, another thread might acquire the API lock
        // and block us when we try to lock again. (see PYMOL-3249)
        // Example: set_key F1, ray async=1
        //          => hitting F1 will block GUI updates, even though ray
        //             is running async
        // PUnlockAPIWhileBlocked(G);

        PFLUSH_REPORT_UNCAUGHT_EXCEPTIONS(G);
        PXDecRef(PYOBJECT_CALLFUNCTION(G->P_inst->parse, "si", buffer.c_str(), 0));
        PFLUSH_REPORT_UNCAUGHT_EXCEPTIONS(G);

        // TODO see above
        // PLockAPIWhileBlocked(G);

        OrthoCommandSetBusy(G, false);
        /* make sure no commands left at this level */
        while(OrthoCommandWaiting(G))
          PFlushFast(G);
        OrthoCommandNest(G, -1);
      }
    }
    PUnblock(G);
  }
  return true;
}

/**
 * Execute commands from the command queue.
 * @return False if the queue was empty
 * @pre GIL
 * @post queue is empty
 */
int PFlushFast(PyMOLGlobals * G)
{
  assert(PyGILState_Check());

  /* NOTE: ASSUMES we currently have blocked Python threads and an unlocked API */
  int did_work = false;
  auto ortho = G->Ortho;
  while(!OrthoCommandIsEmpty(*ortho)){
    auto buffer = OrthoCommandOut(*ortho);
    OrthoCommandSetBusy(G, true);
    OrthoCommandNest(G, 1);
    did_work = true;

    PFLUSH_REPORT_UNCAUGHT_EXCEPTIONS(G);
    PXDecRef(PYOBJECT_CALLFUNCTION(G->P_inst->parse, "si", buffer.c_str(), 0));
    PFLUSH_REPORT_UNCAUGHT_EXCEPTIONS(G);

    OrthoCommandSetBusy(G, false);
    /* make sure no commands left at this level */
    while(OrthoCommandWaiting(G))
      PFlushFast(G);
    OrthoCommandNest(G, -1);
  }
  return did_work;
}

/**
 * Acquire the GIL.
 * @pre no GIL
 * @post GIL
 */
void PBlock(PyMOLGlobals * G)
{
  assert(!PyGILState_Check());

  if(!PAutoBlock(G)) {
    ErrFatal(G, "PBlock", "Threading error detected.  Terminating...");
  }

  assert(PyGILState_Check());
}

/**
 * Acquire the GIL.
 * Return false if the current thread already holds the GIL.
 * @post GIL
 */
int PAutoBlock(PyMOLGlobals * G)
{
#ifndef _PYMOL_EMBEDDED
  SavedThreadRec *SavedThread = G->P_inst->savedThread;

  auto id = PyThread_get_thread_ident();

  for (auto a = MAX_SAVED_THREAD - 1; a; --a) {
    if (SavedThread[a].id == id) {
      assert(!PyGILState_Check());

      // Aquire GIL
      PyEval_RestoreThread(SavedThread[a].state);

      SavedThread[a].id = -1;

      assert(PyGILState_Check());
      return 1;
    }
  }

  assert(PyGILState_Check());
  return 0;
#else
  return 1;
#endif
}

/**
 * Can be called anytime, no lock conditions
 * @return True if the current thread is the GUI thread
 */
int PIsGlutThread(void)
{
  return (PyThread_get_thread_ident() == P_glut_thread_id);
}

/**
 * Release GIL
 * @pre GIL
 * @post no GIL
 */
void PUnblock(PyMOLGlobals * G)
{
#ifndef _PYMOL_EMBEDDED
  assert(PyGILState_Check());

  SavedThreadRec *SavedThread = G->P_inst->savedThread;
  auto a = MAX_SAVED_THREAD - 1;
  /* NOTE: ASSUMES a locked API */

  /* reserve a space while we have a lock */
  for (; a; --a) {
    if (SavedThread[a].id == -1) {
      SavedThread[a].id = PyThread_get_thread_ident();
      break;
    }
  }

  // Release GIL
  SavedThread[a].state = PyEval_SaveThread();

  assert(!PyGILState_Check());
#endif
}

/**
 * Release GIL if flag=true
 */
void PAutoUnblock(PyMOLGlobals * G, int flag)
{
  if(flag)
    PUnblock(G);
}

/**
 * Acquire GIL, release API lock and flush
 * @pre no GIL
 * @post GIL
 * @post no API lock
 */
void PBlockAndUnlockAPI(PyMOLGlobals * G)
{
  PBlock(G);
  PXDecRef(PYOBJECT_CALLFUNCTION(G->P_inst->unlock, "iO", 0, G->P_inst->cmd));
}

/**
 * Acquire API lock
 * @param block_if_busy If false and PyMOL is busy, do not block
 * @return False if API lock could not be acquired
 * @pre no GIL
 */
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

/**
 * Release API lock with flushing
 * @pre no GIL
 * @pre API lock
 */
void PUnlockAPI(PyMOLGlobals * G)
{
  PBlock(G);
  PXDecRef(PYOBJECT_CALLFUNCTION(G->P_inst->unlock, "iO", 0, G->P_inst->cmd));
  PUnblock(G);
}

/**
 * Release API lock without flushing
 * @pre GIL
 * @pre API lock
 * @post no API lock
 */
static void PUnlockAPIWhileBlocked(PyMOLGlobals * G)
{
  PXDecRef(PYOBJECT_CALLFUNCTION(G->P_inst->unlock, "iO", -1, G->P_inst->cmd));
}

/**
 * Acquire API lock
 * @pre GIL
 * @post GIL
 * @post API lock
 */
static void PLockAPIWhileBlocked(PyMOLGlobals * G)
{
  PXDecRef(PYOBJECT_CALLFUNCTION(G->P_inst->lock, "O", G->P_inst->cmd));
}

/**
 * Try to acquire API lock (non-blocking if busy) and release GIL
 * @return False if API lock could not be acquired
 * @pre GIL
 * @post no GIL if API lock could be acquired
 */
int PTryLockAPIAndUnblock(PyMOLGlobals * G)
{
  int result = get_api_lock(G, false);
  if(result) {
    PUnblock(G);
  }
  return result;
}

/**
 * Acquire API lock and release the GIL
 * @pre GIL
 * @post no GIL
 * @post API Lock
 */
void PLockAPIAndUnblock(PyMOLGlobals * G)
{
  assert(PyGILState_Check());

  PXDecRef(PYOBJECT_CALLFUNCTION(G->P_inst->lock, "O", G->P_inst->cmd));
  PUnblock(G);
}

/**
 * Assign a variable in the pymol namespace, like
 * @verbatim
   import pymol
   setattr(pymol, name, value)
   @endverbatim
 * @pre no GIL
 */
void PDefineFloat(PyMOLGlobals * G, const char *name, float value)
{
  assert(!PyGILState_Check());

  char buffer[OrthoLineLength];
  sprintf(buffer, "%s = %f\n", name, value);
  PBlock(G);
  PRunStringModule(G, buffer);
  PUnblock(G);
}

/**
 * If the error indicator is set (means: If a Python exception was raised),
 * then clear it and print the traceback.
 *
 * Special case is CmdException, which will be printed without a traceback.
 */
void PErrPrintIfOccurred(PyMOLGlobals* G)
{
  assert(PyGILState_Check());

  PyObject *type = nullptr, *value = nullptr, *traceback = nullptr;
  PyErr_Fetch(&type, &value, &traceback);

  if (!type) {
    return;
  }

  if (!value || !PyErr_GivenExceptionMatches(type, P_CmdException)) {
    PyErr_Restore(type, value, traceback);
    PyErr_Print();
    return;
  }

  Py_XDECREF(traceback);

  if (PyObject* strobj = PyObject_Str(value)) {
    const char* str = PyUnicode_AsUTF8(strobj);
    assert(str);

    G->Feedback->addColored(str, FB_Errors);
    G->Feedback->add("\n");

    Py_DECREF(strobj);
  } else {
    assert(PyErr_Occurred());
    PyErr_Print();
  }

  Py_DECREF(type);
  Py_DECREF(value);
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
          const char *str = PyString_AsString(obj);
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

static PyObject *PCatchIsAtty(PyObject * self, PyObject * args)
{
  Py_RETURN_FALSE;
}

static PyObject *PCatch_install(PyObject * self, PyObject * args)
{
  PyRun_SimpleString(
      "import sys, pcatch\n"
      "if sys.stdout is not pcatch:"
      "pcatch.closed = False;"
      "pcatch.encoding = 'UTF-8';"
      "sys.stderr = sys.stdout = pcatch");
  return PConvAutoNone(Py_None);
}

static PyMethodDef PCatch_methods[] = {
  {"writelines", PCatchWritelines, METH_VARARGS},
  {"write", PCatchWrite, METH_VARARGS},
  {"flush", PCatchFlush, METH_VARARGS},
  {"isatty", PCatchIsAtty, METH_VARARGS}, // called by pip.main(["install", "..."])
  {"_install", PCatch_install, METH_VARARGS},
  {NULL, NULL}                  /* sentinel */
};

void PCatchInit(void)
{
  assert(PyGILState_Check());

  static struct PyModuleDef moduledef = { PyModuleDef_HEAD_INIT,
    "pcatch", NULL, -1, PCatch_methods };
  PyObject * pcatch = PyModule_Create(&moduledef);
  if (pcatch) {
    PyDict_SetItemString(PyImport_GetModuleDict(), "pcatch", pcatch);
    Py_DECREF(pcatch);
  }
}

namespace pymol
{
GIL_Ensure::GIL_Ensure() {
  state = PyGILState_Ensure();
}
GIL_Ensure::~GIL_Ensure()
{
  PyGILState_Release(state);
}
} // namespace pymol

#endif
