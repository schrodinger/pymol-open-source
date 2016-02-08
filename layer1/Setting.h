
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

typedef char SettingName[255];

/*
 * Setting record for atom/astate/bond/bstate level settings
 */
typedef struct {
  int setting_id;
  union {
    int int_;
    float float_;
    float float3_[3];
  } value;
  int next;                     /* for per-atom setting lists & memory management */
} SettingUniqueEntry;

struct _CSettingUnique {
  OVOneToOne *id2offset;
  OVOneToOne *old2new;
  SettingUniqueEntry *entry;
  int n_alloc, next_free;
};

/*
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
    if (!str_) {
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

struct _CSetting {
  PyMOLGlobals *G;
  ov_size size;
  SettingRec *info;
};

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

void SettingUniqueSet_b(PyMOLGlobals * G, int unique_id, int setting_id, int value);
void SettingUniqueSet_i(PyMOLGlobals * G, int unique_id, int setting_id, int value);
void SettingUniqueSet_f(PyMOLGlobals * G, int unique_id, int setting_id, float value);
void SettingUniqueSet_color(PyMOLGlobals * G, int unique_id, int setting_id, int value);
int SettingUniqueSetTypedValue(PyMOLGlobals * G, int unique_id, int setting_id,
			       int setting_type, const void *value);
#ifndef _PYMOL_NOPY
bool SettingUniqueSetPyObject(PyMOLGlobals * G, int unique_id, int setting_id, PyObject *value);
#endif

int SettingUniqueCheck(PyMOLGlobals * G, int unique_id, int setting_id);
int SettingUniqueGet_b(PyMOLGlobals * G, int unique_id, int setting_id, int *value);
int SettingUniqueGet_i(PyMOLGlobals * G, int unique_id, int setting_id, int *value);
int SettingUniqueGet_f(PyMOLGlobals * G, int unique_id, int setting_id, float *value);
int SettingUniqueGet_color(PyMOLGlobals * G, int unique_id, int setting_id, int *value);
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
void SettingFreeP(CSetting * I);
void SettingInit(PyMOLGlobals * G, CSetting * I);
void SettingPurge(CSetting * I);
void SettingCheckHandle(PyMOLGlobals * G, CSetting ** handle);

#define SettingSet_b SettingSet_i
int SettingSet_i(CSetting * I, int index, int value);
int SettingSet_color_from_3f(CSetting * I, int index, float *vals);
int SettingSet_f(CSetting * I, int index, float value);
int SettingSet_s(CSetting * I, int index, const char *value);
int SettingSet_3f(CSetting * I, int index, float value1, float value2, float value3);
int SettingSet_3fv(CSetting * I, int index, const float *value);

int SettingGetTextValue(PyMOLGlobals * G, CSetting * set1, CSetting * set2, int index,
                        char *buffer);
const char * SettingGetTextPtr(PyMOLGlobals * G, CSetting * set1, CSetting * set2,
                               int index, char *buffer);

int SettingUnset(CSetting * I, int index);

void SettingRestoreDefault(CSetting * I, int index, const CSetting * src=NULL);

bool SettingIsDefaultZero(int index);

int SettingGetType(PyMOLGlobals * G, int index);        /* based on global types, always succeeds */

int SettingGetGlobal_b(PyMOLGlobals * G, int index);    /* always succeed */
int SettingGetGlobal_i(PyMOLGlobals * G, int index);    /* always succeed */
float SettingGetGlobal_f(PyMOLGlobals * G, int index);  /* always succeed */
char *SettingGetGlobal_s(PyMOLGlobals * G, int index);  /* always succeeds */
int SettingGetGlobal_color(PyMOLGlobals * G, int index);        /* always succeed */

void SettingGetGlobal_3f(PyMOLGlobals * G, int index, float *value);    /* always succeeds */
float *SettingGetGlobal_3fv(PyMOLGlobals * G, int index);       /* always succeed */

int SettingSetGlobal_b(PyMOLGlobals * G, int index, int value);
int SettingSetGlobal_i(PyMOLGlobals * G, int index, int value);
int SettingSetGlobal_s(PyMOLGlobals * G, int index, const char *value);
int SettingSetGlobal_f(PyMOLGlobals * G, int index, float value);
int SettingSetGlobal_3f(PyMOLGlobals * G, int index, float value1, float value2,
                        float value3);
int SettingSetSmart_i(PyMOLGlobals * G, CSetting * set1, CSetting * set2, int index,
                      int value);

/* more to come */

int SettingGetIfDefined_i(PyMOLGlobals * G, CSetting * set1, int index, int *value);
int SettingGetIfDefined_b(PyMOLGlobals * G, CSetting * set1, int index, int *value);
int SettingGetIfDefined_f(PyMOLGlobals * G, CSetting * set1, int index, float *value);
int SettingGetIfDefined_s(PyMOLGlobals * G, CSetting * set1, int index, char **value);
int SettingGetIfDefined_3fv(PyMOLGlobals * G, CSetting * set1, int index, float **value);
int SettingGetIfDefined_color(PyMOLGlobals * G, CSetting * set1, int index, int *value);


/* more to come */

int SettingSet_color(CSetting * I, int index, const char *value);

int SettingGet_b(PyMOLGlobals * G, CSetting * set1, CSetting * set2, int index);
int SettingGet_i(PyMOLGlobals * G, CSetting * set1, CSetting * set2, int index);
float SettingGet_f(PyMOLGlobals * G, CSetting * set1, CSetting * set2, int index);
char *SettingGet_s(PyMOLGlobals * G, CSetting * set1, CSetting * set2, int index);
void SettingGet_3f(PyMOLGlobals * G, CSetting * set1, CSetting * set2, int index,
                   float *value);
float *SettingGet_3fv(PyMOLGlobals * G, CSetting * set1, CSetting * set2, int index);
int SettingGet_color(PyMOLGlobals * G, CSetting * set1, CSetting * set2, int index);

int SettingSetFromString(PyMOLGlobals * G, CSetting * I, int index, const char *st);
int SettingStringToTypedValue(PyMOLGlobals * G, int index, const char *st, int *type,
                              int *value);

#ifndef _PYMOL_NOPY
int SettingSetFromTuple(PyMOLGlobals * G, CSetting * I, int index, PyObject * tuple);
PyObject *SettingGetPyObject(PyMOLGlobals * G, CSetting * set1, CSetting * set2, int index);
PyObject *SettingGetTuple(PyMOLGlobals * G, CSetting * set1, CSetting * set2, int index);       /* (type,(value,)) */
PyObject *SettingGetSettingIndices();
PyObject *SettingUniqueGetIndicesAsPyList(PyMOLGlobals * G, int unique_id);
#endif

std::vector<int> SettingGetUpdateList(PyMOLGlobals * G, const char *, int);

void SettingGenerateSideEffects(PyMOLGlobals * G, int index, const char *sele, int state, int quiet);


/* Legacy API below */

int SettingGetIndex(PyMOLGlobals * G, const char *name);
int SettingGetName(PyMOLGlobals * G, int index, SettingName name);
const char * SettingGetName(int index);

PyObject *SettingAsPyList(CSetting * I);
int SettingFromPyList(CSetting * I, PyObject * list);

int SettingSetGlobalsFromPyList(PyMOLGlobals * G, PyObject * list);
PyObject *SettingGetGlobalsAsPyList(PyMOLGlobals * G);


CSetting *SettingNewFromPyList(PyMOLGlobals * G, PyObject * list);


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

/*
 * State index iterator which iterates either over a single state (state >= 0),
 * the current state (state == -2), or all states (state == -1). Takes
 * static singletons into account. Zero iterations if state >= nstate.
 *
 * StateIterator iter(G, I->Obj.Setting, state, I->NState);
 * while(iter.next()) {
 *   printf("in state %d\n", iter.state);
 * }
 */
class StateIterator {
  int end;

public:
  int state;

  StateIterator(PyMOLGlobals * G, CSetting * set, int state_, int nstate);

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

const char * SettingLevelGetName(PyMOLGlobals * G, int index);
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

inline void SettingUniqueSet(PyMOLGlobals * G, int uid, int i, bool v)          { SettingUniqueSet_b(G, uid, i, v); }
inline void SettingUniqueSet(PyMOLGlobals * G, int uid, int i, int v)           { SettingUniqueSet_i(G, uid, i, v); }
inline void SettingUniqueSet(PyMOLGlobals * G, int uid, int i, float v)         { SettingUniqueSet_f(G, uid, i, v); }

template <typename V> void SettingSet(PyMOLGlobals * G, CSetting ** handle, int index, V value) {
  SettingCheckHandle(G, handle);
  SettingSet(*handle, index, value);
}

// global setting
template <typename V> void SettingSet(PyMOLGlobals * G, int index, V value) {
  SettingSet(G->Setting, index, value);
}

#endif
