
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
#ifndef _H_ObjectMap
#define _H_ObjectMap

#include"os_python.h"

#include"CGO.h"
#include"PyMOLObject.h"
#include"Symmetry.h"
#include"Isosurf.h"
#include"vla.h"
#include"Result.h"

#define cMapSourceUndefined 0


/* in order to achieve better backwards compatibility
 * with session files, future map source should be either:
 * cMapSourceCrystallographic or cMapSourceGeneralPurpose */

#define cMapSourceCrystallographic 1
#define cMapSourceCCP4 2
#define cMapSourceGeneralPurpose 3
#define cMapSourceDesc 4
#define cMapSourceFLD 5
#define cMapSourceBRIX 6
#define cMapSourceGRD 7
#define cMapSourceChempyBrick 8
#define cMapSourceVMDPlugin 9
#define cMapSourceObsolete   10

struct ObjectMapState : public CObjectState {
  int Active = false;
  pymol::copyable_ptr<CSymmetry> Symmetry;
  int Div[3];                   /* NOTE: Div is only meaningful for maps defined relative to a unit cell */
  int Min[3], Max[3];           /* valid min and max indices, required. */
  int FDim[4];                  /* Array dimensions with 3 in fourth slot, required */
  int MapSource;
  pymol::copyable_ptr<Isofield> Field;
  float Corner[24];
  std::vector<int> Dim;                     /* this field is redundant and should be eliminated -- if exists, must match FDim */
  std::vector<float> Origin;                /* Origin for non-xtal maps */
  std::vector<float> Range;                 /* Range for non-xtal maps */
  std::vector<float> Grid;                  /* Spacing for non-xtal maps */
  float ExtentMin[3], ExtentMax[3];
  float Mean, SD; /* -- JV for vol */
  pymol::cache_ptr<CGO> shaderCGO;
  /* below not stored */

  int have_range = false;
  float high_cutoff, low_cutoff;
  ObjectMapState(PyMOLGlobals* G);
  ObjectMapState(const ObjectMapState&);
  ObjectMapState& operator=(const ObjectMapState&);
};

struct ObjectMap : public pymol::CObject {
  using StateT = ObjectMapState;

  std::vector<ObjectMapState> State;
  ObjectMap(PyMOLGlobals* G);

  /// Typed version of getObjectState
  StateT* getObjectMapState(int state)
  {
    return static_cast<StateT*>(getObjectState(state));
  }
  const StateT* getObjectMapState(int state) const
  {
    return static_cast<const StateT*>(getObjectState(state));
  }

  // virtual methods
  void update() override;
  void render(RenderInfo* info) override;
  void invalidate(cRep_t rep, cRepInv_t level, int state) override;
  int getNFrame() const override;
  pymol::CObject* clone() const override;
  CSymmetry const* getSymmetry(int state = 0) const override;
  bool setSymmetry(CSymmetry const& symmetry, int state = 0) override;

protected:
  CObjectState* _getObjectState(int state) override;
};

#define cObjectMap_OrthoMinMaxGrid 0

typedef struct ObjectMapDesc {  /* information for creating a new map */
  int mode;
  float Grid[3];
  int Dim[3];
  float MinCorner[3], MaxCorner[3];
  int init_mode;                /* -1 = nothing
                                   0 = zeros
                                   1 = ones */
} ObjectMapDesc;


int ObjectMapNewCopy(PyMOLGlobals * G, const ObjectMap * src, ObjectMap ** result,
                     int source_state, int target_state);
ObjectMapState *ObjectMapNewStateFromDesc(PyMOLGlobals * G, ObjectMap * I,
                                          ObjectMapDesc * md, int state, int quiet);
int ObjectMapStateGetExcludedStats(PyMOLGlobals * G, ObjectMapState * ms, float *vert_vla,
                                   float beyond, float within, float *level);

int ObjectMapValidXtal(ObjectMap * I, int state);
int ObjectMapStateValidXtal(ObjectMapState * ms);

void ObjectMapStateClamp(ObjectMapState * I, float clamp_floor, float clamp_ceiling);

ObjectMap *ObjectMapLoadXPLOR(PyMOLGlobals * G, ObjectMap * obj, const char *fname,
                              int state, int is_file, int quiet);

ObjectMap *ObjectMapLoadCCP4(PyMOLGlobals * G, ObjectMap * obj, const char *fname,
                             int state, int is_string, int bytes, int quiet, int);

ObjectMap *ObjectMapLoadPHI(PyMOLGlobals * G, ObjectMap * obj, const char *fname, int state,
                            int is_string, int bytes, int quiet);

ObjectMap* ObjectMapReadDXStr(PyMOLGlobals*, ObjectMap*, const char* MapStr,
    int bytes, int state, bool quiet);

ObjectMap *ObjectMapLoadDXFile(PyMOLGlobals * G, ObjectMap * obj, const char *fname, int state,
                               int quiet);
ObjectMap *ObjectMapLoadFLDFile(PyMOLGlobals * G, ObjectMap * obj, const char *fname, int state,
                                int quiet);
ObjectMap *ObjectMapLoadBRIXFile(PyMOLGlobals * G, ObjectMap * obj, const char *fname,
                                 int state, int quiet);
ObjectMap *ObjectMapLoadGRDFile(PyMOLGlobals * G, ObjectMap * obj, const char *fname, int state,
                                int quiet);
ObjectMap *ObjectMapLoadACNTFile(PyMOLGlobals * G, ObjectMap * obj, const char *fname,
                                 int state, int quiet);

ObjectMap *ObjectMapLoadChemPyBrick(PyMOLGlobals * G, ObjectMap * I, PyObject * Map,
                                    int state, int discrete, int quiet);
ObjectMap *ObjectMapLoadChemPyMap(PyMOLGlobals * G, ObjectMap * I, PyObject * Map,
                                  int state, int discrete, int quiet);
pymol::Result<> ObjectMapDouble(ObjectMap * I, int state);
pymol::Result<> ObjectMapHalve(ObjectMap * I, int state, int smooth);
pymol::Result<> ObjectMapTrim(ObjectMap * I, int state, float *mn, float *mx, int quiet);
int ObjectMapSetBorder(ObjectMap * I, float level, int state);
int ObjectMapStateSetBorder(ObjectMapState * I, float level);
void ObjectMapStatePurge(PyMOLGlobals * G, ObjectMapState * I);
int ObjectMapStateInterpolate(ObjectMapState * ms, const float *array, float *result, int *flag,
                              int n);
int ObjectMapStateContainsPoint(ObjectMapState * ms, float *point);
ObjectMapState *ObjectMapStatePrime(ObjectMap * I, int state);
void ObjectMapUpdateExtents(ObjectMap * I);

#define ObjectMapStateGetActive(I, state) (I)->getObjectMapState(state)
#define ObjectMapGetState(I, state) (I)->getObjectMapState(state)

PyObject *ObjectMapAsPyList(ObjectMap * I);
int ObjectMapNewFromPyList(PyMOLGlobals * G, PyObject * list, ObjectMap ** result);

int ObjectMapInterpolate(ObjectMap * I, int state, const float *array, float *result, int *flag,
                         int n);

void ObjectMapRegeneratePoints(ObjectMap * om);
void ObjectMapStateRegeneratePoints(ObjectMapState * ms);
int ObjectMapStateGetDataRange(PyMOLGlobals * G, ObjectMapState * ms, float *min,
                               float *max);
int ObjectMapStateGetHistogram(PyMOLGlobals * G, ObjectMapState * ms,
                               int n_points, float limit, float *histogram,
                               float min_arg, float max_arg);

void ObjectMapDump(const ObjectMap* I, const char* fname, int state, int quiet);

ObjectMapState * getObjectMapState(PyMOLGlobals * G, const char * name, int state);

/*========================================================================*/
std::vector<char> ObjectMapStateToCCP4Str(const ObjectMapState * ms, int quiet, int);

template <typename T>
std::vector<char> ObjectMapGetCCP4Str(PyMOLGlobals * G, const T * ptr, int state, int quiet, int format) {
  return ObjectMapStateToCCP4Str(getObjectMapState(G, ptr, state), quiet, format);
}

#endif
