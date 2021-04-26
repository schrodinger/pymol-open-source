
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
#ifndef _H_Setting
#define _H_Setting

#include <vector>
#include <string>

#include"os_python.h"
#include"PyMOLGlobals.h"
#include"OVOneToOne.h"
#include"pymol/memory.h"

typedef char SettingName[255];

/**
 * Setting record for atom/astate/bond/bstate level settings
 */
struct SettingUniqueEntry {
  int setting_id;
  union {
    int int_;
    float float_;
    float float3_[3];
  } value;
  int next;                     /* for per-atom setting lists & memory management */
};

struct CSettingUnique {
  OVOneToOne *id2offset;
  OVOneToOne *old2new;
  SettingUniqueEntry *entry;
  int n_alloc, next_free;
};

/**
 * Setting record for global/object/ostate level settings
 */
struct SettingRec {
  union {
    int int_;
    float float_;
    float float3_[3];
    std::string * str_;
  };

  bool defined;
  bool changed;

private:
  void setChanged() {
    defined = true;
    changed = true;
  }

public:
  void set_i(int value) {
    int_ = value;
    setChanged();
  }

  void set_f(float value) {
    float_ = value;
    setChanged();
  }

  void set_3f(float v0, float v1, float v2) {
    float3_[0] = v0;
    float3_[1] = v1;
    float3_[2] = v2;
    setChanged();
  }

  void set_3f(const float * value) {
    set_3f(value[0], value[1], value[2]);
  }

  void set_s(const char * value) {
    if (!value) {
      delete_s();
    } else if (!str_) {
      str_ = new std::string(value);
    } else {
      str_->assign(value);
    }
    setChanged();
  }

  void delete_s() {
    if (str_) {
      delete str_;
      str_ = NULL;
    }
  }
};

#define cSetting_tuple      -1 // for get_setting_tuple
#define cSetting_blank       0
#define cSetting_boolean     1
#define cSetting_int         2
#define cSetting_float       3
#define cSetting_float3      4
#define cSetting_color       5
#define cSetting_string      6


/* Atomic Settings */

void SettingUniqueDetachChain(PyMOLGlobals * G, int index);

/* New API 
 * NOTE: get commands are not range-checked, so be careful
 * in contrast, set commands expand the current list 
 */

int SettingUniqueSetTypedValue(PyMOLGlobals * G, int unique_id, int setting_id,
			       int setting_type, const void *value);

/*
 * bool overload
 */
inline
int SettingUniqueSetTypedValue(PyMOLGlobals * G, int unique_id, int setting_id,
			       int setting_type, const bool *value) {
  int i = *value;
  return SettingUniqueSetTypedValue(G, unique_id, setting_id, setting_type, &i);
}

bool SettingUniqueUnset(PyMOLGlobals * G, int unique_id, int setting_id);

#ifndef _PYMOL_NOPY
bool SettingUniqueSetPyObject(PyMOLGlobals * G, int unique_id, int setting_id, PyObject *value);
#endif

int SettingUniqueCheck(PyMOLGlobals * G, int unique_id, int setting_id);
PyObject *SettingUniqueGetPyObject(PyMOLGlobals * G, int unique_id, int index);

void SettingUniqueResetAll(PyMOLGlobals * G);
PyObject *SettingUniqueAsPyList(PyMOLGlobals * G);
int SettingUniqueFromPyList(PyMOLGlobals * G, PyObject * list, int partial_restore);
int SettingUniqueConvertOldSessionID(PyMOLGlobals * G, int old_unique_id);

int SettingUniqueCopyAll(PyMOLGlobals * G, int src_unique_id, int dst_unique_id);
void SettingInitGlobal(PyMOLGlobals * G, int alloc, int reset_gui, int use_default);
void SettingStoreDefault(PyMOLGlobals * G);
void SettingPurgeDefault(PyMOLGlobals * G);

void SettingFreeGlobal(PyMOLGlobals * G);

CSetting *SettingNew(PyMOLGlobals * G);
CSetting* SettingCopyAll(PyMOLGlobals* G, const CSetting* src, CSetting* dst);
void SettingFreeP(CSetting * I);
void SettingCheckHandle(PyMOLGlobals * G, pymol::copyable_ptr<CSetting>& handle);

#define SettingSet_b SettingSet_i
int SettingSet_i(CSetting * I, int index, int value);
int SettingSet_f(CSetting * I, int index, float value);
int SettingSet_s(CSetting * I, int index, const char *value);
int SettingSet_3fv(CSetting * I, int index, const float *value);

int SettingGetTextValue(PyMOLGlobals * G, const CSetting * set1, const CSetting * set2, int index,
                        char *buffer);
const char * SettingGetTextPtr(PyMOLGlobals * G, const CSetting * set1, const CSetting * set2,
                               int index, char *buffer);

int SettingUnset(CSetting * I, int index);

void SettingRestoreDefault(CSetting * I, int index, const CSetting * src=NULL);

int SettingGetType(int index);
inline int SettingGetType(PyMOLGlobals *, int index) {
  return SettingGetType(index);
}

template <typename V> inline int SettingGetType();
template <> inline int SettingGetType<bool>()           { return cSetting_boolean; }
template <> inline int SettingGetType<int>()            { return cSetting_int; }
template <> inline int SettingGetType<float>()          { return cSetting_float; }
template <> inline int SettingGetType<const float*>()   { return cSetting_float3; }
template <> inline int SettingGetType<float*>()         { return cSetting_float3; }
template <> inline int SettingGetType<const char*>()    { return cSetting_string; }
template <> inline int SettingGetType<char*>()          { return cSetting_string; }

#define SettingSetGlobal_b SettingSet<bool>
#define SettingSetGlobal_i SettingSet<int>
#define SettingSetGlobal_s SettingSet<const char *>
#define SettingSetGlobal_f SettingSet<float>

int SettingSetSmart_i(PyMOLGlobals * G, CSetting * set1, CSetting * set2, int index,
                      int value);

int SettingSet_color(CSetting * I, int index, const char *value);

int SettingSetFromString(PyMOLGlobals * G, CSetting * I, int index, const char *st);
int SettingStringToTypedValue(PyMOLGlobals * G, int index, const char *st, int *type,
                              int *value);

#ifndef _PYMOL_NOPY
int SettingSetFromTuple(PyMOLGlobals * G, CSetting * I, int index, PyObject * tuple);
PyObject *SettingGetPyObject(PyMOLGlobals * G, const CSetting * set1, const CSetting * set2, int index);
PyObject *SettingGetTuple(PyMOLGlobals * G, const CSetting * set1, const CSetting * set2, int index);       /* (type,(value,)) */
PyObject *SettingGetSettingIndices();
PyObject *SettingUniqueGetIndicesAsPyList(PyMOLGlobals * G, int unique_id);
#endif

std::vector<int> SettingGetUpdateList(PyMOLGlobals * G, const char * name="", int state=0);

void SettingGenerateSideEffects(PyMOLGlobals * G, int index, const char *sele, int state, int quiet);

int SettingGetIndex(PyMOLGlobals * G, const char *name);
int SettingGetName(PyMOLGlobals * G, int index, SettingName name);
const char * SettingGetName(int index);

PyObject *SettingAsPyList(CSetting * I, bool incl_blacklisted=false);
int SettingFromPyList(CSetting * I, PyObject * list);

int SettingSetGlobalsFromPyList(PyMOLGlobals * G, PyObject * list);
PyObject *SettingGetGlobalsAsPyList(PyMOLGlobals * G);


CSetting *SettingNewFromPyList(PyMOLGlobals * G, PyObject * list);

int SettingCheckFontID(PyMOLGlobals * G, CSetting * set1, CSetting * set2, int font_id);

// The following defines the enum with all cSetting_<settingname> indices
#include "SettingInfo.h"

#define cStereo_default              0  /* stereo_mode=0 only used for startup */
#define cStereo_quadbuffer           1
#define cStereo_crosseye             2
#define cStereo_walleye              3
#define cStereo_geowall              4
#define cStereo_sidebyside           5
#define cStereo_stencil_by_row       6
#define cStereo_stencil_by_column    7
#define cStereo_stencil_checkerboard 8
#define cStereo_stencil_custom       9  /* for hardware developers to use */
#define cStereo_anaglyph            10
#define cStereo_dynamic             11  /* dynamic polarization */
#define cStereo_clone_dynamic       12
#define cStereo_openvr              13

struct CSetting {
  PyMOLGlobals* G;
  SettingRec info[cSetting_INIT]{};

  CSetting(PyMOLGlobals*);
  CSetting(const CSetting&);
  CSetting& operator=(const CSetting&);
  ~CSetting();
};

namespace pymol
{
struct CObject;
}

/**
 * State index iterator which iterates either over a single state (state >= 0),
 * the current state (state == -2), or all states (state == -1). Takes
 * static singletons into account. Zero iterations if state >= nstate.
 *
 *     StateIterator iter(G, I->Setting, state, I->NState);
 *     while (iter.next()) {
 *         printf("in state %d\n", iter.state);
 *     }
 */
class StateIterator {
  StateIndex_t end;

public:
  StateIndex_t state;

  StateIterator(PyMOLGlobals*, CSetting*, StateIndex_t, int nstate);
  StateIterator(pymol::CObject*, StateIndex_t);

  bool next() {
    return (++state < end);
  };
};

/*
 * Setting levels (global, object, ...)
 */

enum {
  cSettingLevel_unused = 0,
  cSettingLevel_global,
  cSettingLevel_object,
  cSettingLevel_ostate,
  cSettingLevel_atom,
  cSettingLevel_astate,
  cSettingLevel_bond,
  cSettingLevel_bstate
};

// Setting level info table
extern const struct SettingLevelInfoType {
  // name of this level, for feedback
  const char * name;
  // bitmask of valid (sub-)levels
  unsigned char mask;
} SettingLevelInfo[];

const char * SettingLevelGetName(unsigned index);
bool SettingLevelCheckMask(PyMOLGlobals * G, int index, unsigned char mask);
bool SettingLevelCheck(PyMOLGlobals * G, int index, unsigned char level);

bool CPyMOLInitSetting(OVLexicon * Lex, OVOneToOne * Setting);
extern "C" OVreturn_word get_setting_id(CPyMOL * I, const char *setting);

/*
 * Overloaded setters for templatted programming
 */

inline void SettingSet(CSetting * s, int i, bool v)         { SettingSet_b(s, i, v); }
inline void SettingSet(CSetting * s, int i, int v)          { SettingSet_i(s, i, v); }
inline void SettingSet(CSetting * s, int i, long int v)     { SettingSet_i(s, i, v); }
inline void SettingSet(CSetting * s, int i, float v)        { SettingSet_f(s, i, v); }
inline void SettingSet(CSetting * s, int i, const char *v)  { SettingSet_s(s, i, v); }
inline void SettingSet(CSetting * s, int i, const float *v) { SettingSet_3fv(s, i, v); }

template <typename V>
void SettingUniqueSet(PyMOLGlobals * G, int uid, int index, V value) {
  SettingUniqueSetTypedValue(G, uid, index, SettingGetType<V>(), &value);
}

template <typename V> void SettingSet(PyMOLGlobals * G, pymol::copyable_ptr<CSetting>& handle, int index, V value) {
  SettingCheckHandle(G, handle);
  SettingSet(handle.get(), index, value);
}

// global setting
template <typename V> bool SettingSet(PyMOLGlobals * G, int index, V value) {
  SettingSet(G->Setting, index, value);
  return true;
}

// light setting array
extern int light_setting_indices[];

/*
 * Getters for templatted programming
 */

const CSetting * _SettingGetFirstDefined(int index,
    PyMOLGlobals * G,
    const CSetting * set1,
    const CSetting * set2);

template <typename V> V _SettingGet(int index, const CSetting *);
template <typename V> V SettingGet(PyMOLGlobals * G,
    const CSetting * set1,
    const CSetting * set2, int index) {
  return SettingGet<V>(index, _SettingGetFirstDefined(index, G, set1, set2));
}
template <typename V> V SettingGet(PyMOLGlobals * G, int index) {
  return SettingGet<V>(index, G->Setting);
}

template <typename V> V SettingGet(int index, const CSetting* set)
{
  using T = typename std::conditional<std::is_enum<V>::value, int, V>::type;
  return static_cast<V>(_SettingGet<T>(index, set));
}

/**
 * Get setting value if it's defined in the given set, and assign value to
 * `out` variable and return true. Otherwise return false and leave `out`
 * untouched.
 */
template <typename V>
bool SettingGetIfDefined(const CSetting * s, int index, V * out) {
  if (s && s->info[index].defined) {
    *out = SettingGet<V>(index, s);
    return true;
  }
  return false;
}

/**
 * Return setting value if it's defined in the given set, or `default_` otherwise.
 */
template <typename V>
V SettingGetWD(const CSetting* s, int index, V default_) {
  V out;
  if (SettingGetIfDefined<V>(s, index, &out))
    return out;
  return default_;
}

#define SettingGet_b            SettingGet<bool>
#define SettingGet_i            SettingGet<int>
#define SettingGet_color        SettingGet<int>
#define SettingGet_f            SettingGet<float>
#define SettingGet_s            SettingGet<const char *>
#define SettingGet_3fv          SettingGet<const float *>

#define SettingGetGlobal_b      SettingGet<bool>
#define SettingGetGlobal_i      SettingGet<int>
#define SettingGetGlobal_color  SettingGet<int>
#define SettingGetGlobal_f      SettingGet<float>
#define SettingGetGlobal_s      SettingGet<const char *>
#define SettingGetGlobal_3fv    SettingGet<const float *>

#define SettingGetIfDefined_b(G, s, i, o) SettingGetIfDefined(s, i, o)
#define SettingGetIfDefined_i(G, s, i, o) SettingGetIfDefined(s, i, o)

/*
 * templatted unique settings
 */

bool SettingUniqueGetTypedValuePtr(PyMOLGlobals * G, int unique_id, int index,
    int setting_type, void * out);

/*
 * bool overload
 */
inline
bool SettingUniqueGetTypedValuePtr(PyMOLGlobals * G, int unique_id, int index,
    int setting_type, bool * out) {
  int i = *out;
  bool r = SettingUniqueGetTypedValuePtr(G, unique_id, index, setting_type, &i);
  *out = i;
  return r;
}

/**
 * SettingGetIfDefined() equivalent for unique settings.
 */
template <typename V>
bool SettingUniqueGetIfDefined(PyMOLGlobals * G, int unique_id, int index, V * out) {
  return SettingUniqueGetTypedValuePtr(G, unique_id, index,
      SettingGetType<V>(), out);
}

#endif
