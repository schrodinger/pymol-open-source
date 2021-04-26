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
#include"os_python.h"
#include"os_predef.h"
#include"os_std.h"
#include"os_gl.h"

#include <set>
#include <algorithm>
#include <map>
#include <cassert>
#include <clocale>
#include <vector>

#include "pymol/type_traits.h"

#ifdef PYMOL_OPENMP
#include <omp.h>
#endif

#include"Version.h"
#include"main.h"
#include"Base.h"
#include"OOMac.h"
#include"Executive.h"
#include"SpecRec.h"
#include"ObjectMap.h"
#include"ObjectMesh.h"
#include"ObjectDist.h"
#include"ObjectSurface.h"
#include"ObjectSlice.h"
#include"ObjectAlignment.h"
#include"ObjectGroup.h"
#include"ObjectVolume.h"
#include"ObjectCallback.h"
#include"ObjectMap.h"
#include"ObjectMolecule3.h"
#include"ListMacros.h"
#include"List.h"
#include"MyPNG.h"
#include"Ortho.h"
#include"Scene.h"
#include"ScenePicking.h"
#include"SceneRay.h"
#include"Selector.h"
#include"Vector.h"
#include"Color.h"
#include"Setting.h"
#include"Matrix.h"
#include"P.h"
#include"PConv.h"
#include"Match.h"
#include"ObjectCGO.h"
#include"Util.h"
#include"Util2.h"
#include"Wizard.h"
#include"ScrollBar.h"
#include"Movie.h"
#include"ObjectGadgetRamp.h"
#include"SculptCache.h"
#include"Control.h"
#include"Menu.h"
#include"Map.h"
#include"Editor.h"
#include"RepDot.h"
#include"Seq.h"
#include"Text.h"
#include"PyMOL.h"
#include"PyMOLOptions.h"
#include"Tracker.h"
#include"TrackerList.h"
#include"Word.h"
#include"main.h"
#include"Parse.h"
#include"PlugIOManager.h"
#include "Lex.h"
#include "List.h"
#include "AtomIterators.h"

#include"OVContext.h"
#include"OVLexicon.h"
#include"OVOneToOne.h"
#include"OVOneToAny.h"

#include"ShaderMgr.h"
#include"File.h"
#include"FileStream.h"
#include"ExecutiveLoad.h"

#include "MovieScene.h"
#include "Texture.h"

#ifdef _PYMOL_OPENVR
#include "OpenVRMode.h"
#endif

#ifndef _PYMOL_NOPY
#include "ce_types.h"
#endif

#include <glm/glm.hpp>
#include <glm/gtc/type_ptr.hpp>

#ifdef WIN32
#define mkstemp _mktemp_s
#endif

#define cTempRectSele "_rect"
#define cLeftButSele "lb"
#define cIndicateSele "indicate"

#ifndef NO_MMLIBS
#include "mmpymolx.h"
#endif

/**
 * Macro with error handling to get a selection ID from a non-empty expression.
 * @param expression Selection expression
 * @param[out] tmpsele Name of the pymol::Result<SelectorTmp> instance
 * @param[out] seleID Name of the SelectorID_t instance
 */
#define SETUP_SELE(expression, tmpsele, seleID)                                \
  auto tmpsele = SelectorTmp::make(G, expression);                             \
  p_return_if_error(tmpsele);                                                  \
  SelectorID_t seleID = tmpsele->getIndex();                                 \
  if (seleID < 0) {                                                            \
    return pymol::Error("This should not happen - PyMOL may have a bug");      \
  }                                                                            \
  assert(seleID != cSelectionInvalid)

/**
 * Macro with error handling to get a selection ID from a non-empty expression
 * `sN`. Assigns `tmpseleN` and `seleN`.
 */
#define SETUP_SELE_DEFAULT(N) SETUP_SELE(s##N, tmpsele##N, sele##N)

/**
 * Macro with error handling to get a selection ID from a non-empty expression
 * `sN`. Assigns `tmpseleN` and `seleN`.
 * Adds "Selection N: " prefix to error messages.
 * @param same_value If `sN` is "same" then use this value for `seleN`
 */
#define SETUP_SELE_DEFAULT_PREFIXED(N, same_value)                             \
  pymol::Result<SelectorTmp> tmpsele##N;                                       \
  SelectorID_t sele##N = same_value;                                           \
  if (!WordMatchExact(G, s##N, cKeywordSame, true)) {                          \
    tmpsele##N = SelectorTmp::make(G, s##N);                                   \
    p_return_if_error_prefixed(tmpsele##N, "Selection " #N ": ");              \
    sele##N = (tmpsele##N)->getIndex();                                        \
  }                                                                            \
  if (sele##N == cSelectionInvalid) {                                          \
    return pymol::Error("Invalid selection " #N);                              \
  }

/* routines that still need to be updated for Tracker list iteration

Low priority:

ExecutiveSculptIterate
ExecutiveSculptActivate
ExecutiveSculptDeactivate
ExecutiveRenameObjectAtoms
ExecutiveSpheroid

*/

int ExecutiveGetNamesListFromPattern(PyMOLGlobals * G, const char *name,
                                     int allow_partial, int expand_groups);

pymol::TrackerAdapter<SpecRec> ExecutiveGetSpecRecsFromPattern(PyMOLGlobals* G,
    pymol::zstring_view str, bool allow_partial, bool expand_groups)
{
  return pymol::TrackerAdapter<SpecRec>(
      G->Executive->Tracker, ExecutiveGetNamesListFromPattern(
                                 G, str.c_str(), allow_partial, expand_groups));
}

/**
 * Get the list of objects which match the pattern, or all objects
 * if pattern equals 'all'.
 * @param str pattern string provided by user
 * @return List of borrowed object pointers
 */
static std::vector<pymol::CObject*> ExecutiveGetObjectsFromPattern(
    PyMOLGlobals* G, pymol::zstring_view pattern)
{
  std::vector<pymol::CObject*> objects;

  for (auto& rec : ExecutiveGetSpecRecsFromPattern(G, pattern)) {
    switch (rec.type) {
    case cExecObject:
      objects.push_back(rec.obj);
      break;
    case cExecAll:
      for (auto& rec : pymol::make_list_adapter(G->Executive->Spec)) {
        if (rec.type == cExecObject) {
          objects.push_back(rec.obj);
        }
      }
    }
  }

  return objects;
}

static void ExecutiveSpecEnable(PyMOLGlobals * G, SpecRec * rec, int parents, int log);
static void ExecutiveSetAllRepVisMask(PyMOLGlobals * G, int repmask, int state);
static int ExecutiveSetObjectMatrix2(PyMOLGlobals * G, pymol::CObject * obj, int state,
                                     double *matrix);
static int ExecutiveGetObjectMatrix2(PyMOLGlobals * G, pymol::CObject * obj, int state,
                                     double **matrix, int incl_ttt);
static pymol::Result<> ExecutiveTransformObjectSelection2(
    PyMOLGlobals* G, pymol::CObject* obj, int state, const char* s1, int log,
    const float* matrix, int homogenous, int global);

/**
 * ObjectIterator methods
 */
void ObjectIterator::reset() {
  rec = G->Executive->Spec;

  // DEBUG assume first element is always "all"
  if (rec->type != cExecAll)
    printf("Error: first SpecRec is not cExecAll\n");
}

bool ObjectIterator::next() {
  if (!rec || !(rec = rec->next))
    return false;

  if (rec->type != cExecObject)
    return next();

  return true;
}

pymol::CObject * ObjectIterator::getObject() {
  return rec->obj;
}

/**
 * Find object of given type, or delete object if it exists but has the wrong
 * type.
 * @param name Object name
 * @return NULL if object can't be found or has the wrong type
 */
template <typename ObjectT>
ObjectT* ExecutiveFindOrDeleteObject(PyMOLGlobals* G, pymol::zstring_view name)
{
  auto anyObj = ExecutiveFindObjectByName(G, name.c_str());
  auto obj = dynamic_cast<ObjectT*>(anyObj);
  if (anyObj && !obj) {
    // incompatible object with the same name
    ExecutiveDelete(G, name.c_str());
  }
  return obj;
}

/**
 * True if `rec` and all its parent groups are enabled
 */
static bool SpecRecIsEnabled(const SpecRec * rec) {
  while (rec->visible && (rec = rec->group)) {}
  return !rec;
}

bool ExecutiveIsObjectType(const SpecRec& rec, int cObjectType)
{
  return rec.type == cExecObject && rec.obj->type == cObjectType;
}

/* ======================================================= */

static void ReportEnabledChange(PyMOLGlobals * G, SpecRec *rec){
#ifdef _PYMOL_LIB
  if (G->CallbackObject && G->enabledCallback){
    G->enabledCallback(G->CallbackObject, rec->name, rec->visible);
  }
#endif
  OrthoInvalidateDoDraw(G);
  ExecutiveInvalidateSelectionIndicatorsCGO(G);
}

int ExecutiveGroupMotionModify(PyMOLGlobals *G, pymol::CObject *group, int action, 
                                int index, int count, int target, int freeze)
{
  CExecutive *I = G->Executive;
  int result = true;
  CTracker *I_Tracker = I->Tracker;
  int list_id = ExecutiveGetExpandedGroupList(G,group->Name);
  int iter_id = TrackerNewIter(I_Tracker, 0, list_id);
  SpecRec *rec;
  while(TrackerIterNextCandInList(I_Tracker, iter_id, (TrackerRef **) (void *) &rec)) {
    if(rec) {
      switch (rec->type) {
      case cExecObject:
        if(rec->obj->type != cObjectGroup) {
          ObjectMotionModify(rec->obj, action, index, count, target, freeze, true);
        }
        break;
      }
    }
  }
  TrackerDelList(I_Tracker, list_id);
  TrackerDelIter(I_Tracker, iter_id);
  return result;
}

int ExecutiveGroupMotion(PyMOLGlobals *G, pymol::CObject *group,int action, int first,
                         int last, float power, float bias,
                         int simple, float linear, int wrap,
                         int hand, int window, int cycles, int state, int quiet)
{
  CExecutive *I = G->Executive;
  int result = true;
  CTracker *I_Tracker = I->Tracker;
  int list_id = ExecutiveGetExpandedGroupList(G,group->Name);
  int iter_id = TrackerNewIter(I_Tracker, 0, list_id);
  SpecRec *rec;
  while(TrackerIterNextCandInList(I_Tracker, iter_id, (TrackerRef **) (void *) &rec)) {
    if(rec) {
      switch (rec->type) {
      case cExecObject:
        if(rec->obj->type != cObjectGroup) {
          ObjectMotion(rec->obj,action,first,last,power,bias,simple,linear,wrap,hand,window,cycles,state,quiet);
        }
        break;
      }
    }
  }
  TrackerDelList(I_Tracker, list_id);
  TrackerDelIter(I_Tracker, iter_id);
  return result;
}

int ExecutiveGroupCombineTTT(PyMOLGlobals *G, pymol::CObject *group, const float *ttt, int reverse_order, int store)
{
  CExecutive *I = G->Executive;
  int result = true;
  CTracker *I_Tracker = I->Tracker;
  int list_id = ExecutiveGetExpandedGroupList(G,group->Name);
  int iter_id = TrackerNewIter(I_Tracker, 0, list_id);
  SpecRec *rec;
  while(TrackerIterNextCandInList(I_Tracker, iter_id, (TrackerRef **) (void *) &rec)) {
    if(rec) {
      switch (rec->type) {
      case cExecObject:
        if(rec->obj->type != cObjectGroup) {
          ObjectCombineTTT(rec->obj, ttt, reverse_order, store);
        }
        break;
      }
    }
  }
  TrackerDelList(I_Tracker, list_id);
  TrackerDelIter(I_Tracker, iter_id);
  return result;
}

int ExecutiveGroupTranslateTTT(PyMOLGlobals *G, pymol::CObject *group, const float *v, int store)
{
  CExecutive *I = G->Executive;
  int result = true;
  CTracker *I_Tracker = I->Tracker;
  int list_id = ExecutiveGetExpandedGroupList(G,group->Name);
  int iter_id = TrackerNewIter(I_Tracker, 0, list_id);
  SpecRec *rec;
  while(TrackerIterNextCandInList(I_Tracker, iter_id, (TrackerRef **) (void *) &rec)) {
    if(rec) {
      switch (rec->type) {
      case cExecObject:
        if(rec->obj->type != cObjectGroup) {
          ObjectTranslateTTT(rec->obj,v, store);
        }
        break;
      }
    }
  }
  TrackerDelList(I_Tracker, list_id);
  TrackerDelIter(I_Tracker, iter_id);
  return result;
}


int ExecutiveMotionView(PyMOLGlobals *G, int action, int first,
              int last, float power, float bias,
              int simple, float linear, const char *name, int wrap,
              int hand, int window, int cycles,
              const char *scene_name, float scene_cut, int state, int quiet, int autogen)
{
  int ok = true;

  CExecutive *I = G->Executive;

  if(wrap < 0) {
    wrap = SettingGetGlobal_b(G, cSetting_movie_loop);
  }

  if((!name)||(!name[0])||(!strcmp(name,cKeywordNone))||
     (!strcmp(name,cKeywordAll))||(!strcmp(name,cKeywordSame))) { 
    if(autogen) { 
      /* if autogenerated, then use current settings */
      power  = SettingGetGlobal_f(G, cSetting_motion_power);
      bias   = SettingGetGlobal_f(G, cSetting_motion_bias);
      linear = SettingGetGlobal_f(G, cSetting_motion_linear);
      hand   = SettingGetGlobal_i(G, cSetting_motion_hand);
    }
    /* camera */

    ok = MovieView(G, action, first, last, power,
                   bias, true, linear, wrap, hand, window, cycles,
                   scene_name, scene_cut, state, quiet);

    if(name && name[0] && strcmp(name, cKeywordNone)) {
      /* also do all other objects */
      SpecRec *rec = NULL;
      while(ListIterate(I->Spec, rec, next)) {
        switch(rec->type) {
        case cExecObject:
          if(autogen) {
            power  = SettingGet_f(G, NULL, rec->obj->Setting.get(), cSetting_motion_power);
            bias   = SettingGet_f(G, NULL, rec->obj->Setting.get(), cSetting_motion_bias);
            simple = SettingGet_i(G, NULL, rec->obj->Setting.get(), cSetting_motion_simple);
            linear = SettingGet_f(G, NULL, rec->obj->Setting.get(), cSetting_motion_linear);
            hand   = SettingGet_i(G, NULL, rec->obj->Setting.get(), cSetting_motion_hand);
          }
          if((ObjectGetSpecLevel(rec->obj,0)>=0)||(!strcmp(name,cKeywordAll))) {
            ok = ObjectMotion(rec->obj, action, first, last, power, bias,
                              simple < 0 ? 0 : 1, 
                              linear, wrap, hand, window, cycles, state, quiet);
          }
          break;
        }
      }
    }
  } else { /* pattern */
    CTracker *I_Tracker = I->Tracker;
    SpecRec *rec = NULL;
    int list_id = ExecutiveGetNamesListFromPattern(G, name, true, true);
    int iter_id = TrackerNewIter(I_Tracker, 0, list_id);
    while(TrackerIterNextCandInList(I_Tracker, iter_id, (TrackerRef **) (void *) &rec)) {
      if(rec) {

        switch (rec->type) {
        case cExecObject: 
          if(autogen) {
            power  = SettingGet_f(G, NULL, rec->obj->Setting.get(), cSetting_motion_power);
            bias   = SettingGet_f(G, NULL, rec->obj->Setting.get(), cSetting_motion_bias);
            simple = SettingGet_i(G, NULL, rec->obj->Setting.get(), cSetting_motion_simple);
            linear = SettingGet_f(G, NULL, rec->obj->Setting.get(), cSetting_motion_linear);
            hand   = SettingGet_i(G, NULL, rec->obj->Setting.get(), cSetting_motion_hand);
          }
          
          ok = ObjectMotion(rec->obj, action, first, last, power, bias,
                            simple < 0 ? 0 : simple, 
                              linear, wrap, hand, window, cycles, state, quiet);
          break;
        }
      }
    }
    TrackerDelList(I_Tracker, list_id);
    TrackerDelIter(I_Tracker, iter_id);

    // fix for PYMOL-1465
    OrthoReshape(G, -1, -1, false);
  }
  ExecutiveCountMotions(G);
  return ok;
}

pymol::Result<> ExecutiveMotionViewModify(PyMOLGlobals* G, int action,
    int index, int count, int target, const char* name, int freeze, int quiet)
{
  CExecutive *I = G->Executive;
  if((!name)||(!name[0])||
     (!strcmp(name,cKeywordNone))||
     (!strcmp(name,cKeywordSame))||
     (!strcmp(name,cKeywordAll))) {
    /* camera */
    if(MovieGetSpecLevel(G,0)>=0) {
      MovieViewModify(G, action, index, count, target, true, true);
    }
    if((!name) || strcmp(name, cKeywordNone)) {
      /* also do all other objects */
      SpecRec *rec = NULL;
      while(ListIterate(I->Spec, rec, next)) {
        switch(rec->type) {
        case cExecObject:
          if(ObjectGetSpecLevel(rec->obj,0)>=0) {
            ObjectMotionModify(rec->obj, action, index, count, target, true, true);
          }        
          break;
        }
      }
      ExecutiveMotionTrim(G);
    } else {
      ExecutiveMotionExtend(G,true);
    }
    
    if((!freeze) && SettingGetGlobal_i(G,cSetting_movie_auto_interpolate)) {
      ExecutiveMotionReinterpolate(G);
    }
  } else { /* pattern */
    CTracker *I_Tracker = I->Tracker;
    SpecRec *rec = NULL;
    int list_id = ExecutiveGetNamesListFromPattern(G, name, true, true);
    int iter_id = TrackerNewIter(I_Tracker, 0, list_id);
    while(TrackerIterNextCandInList(I_Tracker, iter_id, (TrackerRef **) (void *) &rec)) {
      if(rec) {
        switch (rec->type) {
        case cExecObject: 
          if(ObjectGetSpecLevel(rec->obj,0)>=0) {
            /* only modify objects with motion matrices */
            ObjectMotionModify(rec->obj, action, index, count, target, freeze, false);
          }
        }
      }
    }
    TrackerDelList(I_Tracker, list_id);
    TrackerDelIter(I_Tracker, iter_id);
  }
  ExecutiveCountMotions(G);
  SceneCountFrames(G);
  return {};
}

void ExecutiveMotionReinterpolate(PyMOLGlobals * G)
{
  CExecutive *I = G->Executive;
  SpecRec *rec = NULL;
  while(ListIterate(I->Spec, rec, next)) {
    switch(rec->type) {
    case cExecAll:
      if(MovieGetSpecLevel(G,0)>=0) {
        MovieViewReinterpolate(G);
      }
      break;
    case cExecObject:
      if(ObjectGetSpecLevel(rec->obj,0)>=0) {
        ObjectMotionReinterpolate(rec->obj);
      }        
      break;
    }
  }
}

void ExecutiveMotionTrim(PyMOLGlobals * G)
{
  int n_frame = MovieGetLength(G);
  CExecutive *I = G->Executive;
  SpecRec *rec = NULL;
  while(ListIterate(I->Spec, rec, next)) {
    switch(rec->type) {
    case cExecObject:
      if(ObjectGetSpecLevel(rec->obj,0)>=0) {
        ObjectMotionTrim(rec->obj,n_frame);
      }        
      break;
    }
  }
}

void ExecutiveMotionExtend(PyMOLGlobals * G, int freeze)
{
  int n_frame = 0;
  int max_length = 0;
  CExecutive *I = G->Executive;
  SpecRec *rec = NULL;
  if(MovieGetSpecLevel(G,-1)>0)
    n_frame = MovieGetLength(G);
  while(ListIterate(I->Spec, rec, next)) {
    switch(rec->type) {
    case cExecObject:
      if(ObjectGetSpecLevel(rec->obj,-1)>0) {
        int length = ObjectMotionGetLength(rec->obj);
        if(max_length < length)
          max_length = length;
      }
      break;
    }
  }
  if(max_length) {
    if(n_frame < max_length)
      MovieViewTrim(G,max_length);
    while(ListIterate(I->Spec, rec, next)) {
      switch(rec->type) {
      case cExecObject:
        if(ObjectGetSpecLevel(rec->obj,-1)>0) {
          ObjectMotionTrim(rec->obj,max_length);
        }        
        break;
      }
    }
  }
  if((!freeze) && SettingGetGlobal_i(G,cSetting_movie_auto_interpolate)) {
    ExecutiveMotionReinterpolate(G);
  }

}

int ExecutiveCountMotions(PyMOLGlobals * G)
{
  int count = 0;
  CExecutive *I = G->Executive;
  if(MovieGetLength(G)) {
    SpecRec *rec = NULL;
    while(ListIterate(I->Spec, rec, next)) {
      switch(rec->type) {
      case cExecAll:
        if(MovieGetSpecLevel(G,0)>=0)
          count++;
        break;
      case cExecObject:
        if(ObjectGetSpecLevel(rec->obj,0)>=0)
          count++;
        break;
      }
    }
  }

  if (count < 1 && SceneGetNFrame(G) > 1)
    count = 1;
  
  if(count != I->LastMotionCount) {
    if(SettingGetGlobal_i(G,cSetting_movie_panel)) {
      OrthoDoViewportWhenReleased(G);
    }
    I->LastMotionCount = count;
  }
  
  return (count);
}

void ExecutiveMotionDraw(PyMOLGlobals * G, BlockRect *rect, int expected ORTHOCGOARG)
{
  CExecutive *I = G->Executive;
  SpecRec *rec = NULL;
  int frames = MovieGetLength(G);
  BlockRect draw_rect = *rect;
  int count = 0;
  int height = rect->top - rect->bottom;
  while(ListIterate(I->Spec, rec, next)) {
    switch(rec->type) {
    case cExecAll:
      if(MovieGetSpecLevel(G,0)>=0) {
        int presentation = SettingGetGlobal_b(G, cSetting_presentation);
        if(presentation) { 
          expected = 1;
        }
        draw_rect.top = rect->top - (height * count) / expected;
        draw_rect.bottom = rect->top - (height * (count + 1)) / expected;
        MovieDrawViewElem(G,&draw_rect,frames ORTHOCGOARGVAR);
        count++;
        if(presentation) { 
          goto done;
        }
          
      }
      break;
    case cExecObject:
      if(ObjectGetSpecLevel(rec->obj,0)>=0) {
        draw_rect.top = rect->top - (height * count) / expected;
        draw_rect.bottom = rect->top - (height * (count + 1)) / expected;
        ObjectDrawViewElem(rec->obj,&draw_rect,frames ORTHOCGOARGVAR);
        count++;
      }
      break;
    }
  }
 done:
  return;
}

void ExecutiveMotionMenuActivate(PyMOLGlobals * G, BlockRect *rect, int expected, int passive, 
                                 int x, int y, int same)
{
  CExecutive *I = G->Executive;
  SpecRec *rec = NULL;
  BlockRect draw_rect = *rect;
  int count = 0;
  int height = rect->top - rect->bottom;
  if(same) {
    if(MovieGetSpecLevel(G,0)>=0) {
      int n_frame = MovieGetLength(G);
      int frame = MovieXtoFrame(G, &draw_rect, n_frame, x, false);
      WordType frame_str = "0";
      if((frame>=0) && (frame<n_frame)) {
        sprintf(frame_str,"%d",frame+1);
      }
      MenuActivate2Arg(
          G, x, y, x, y, passive, "obj_motion", cKeywordSame, frame_str);
    }
  } else {
  while(ListIterate(I->Spec, rec, next)) {
    switch(rec->type) {
    case cExecAll:
      if(MovieGetSpecLevel(G,0)>=0) {
        draw_rect.top = rect->top - (height * count) / expected;
        draw_rect.bottom = rect->top - (height * (count + 1)) / expected;
        if((y>draw_rect.bottom) && (y<draw_rect.top)) {
          int n_frame = MovieGetLength(G);
          int frame = MovieXtoFrame(G, &draw_rect, n_frame, x, false);
          WordType frame_str = "0";
          if((frame>=0) && (frame<n_frame)) {
            sprintf(frame_str,"%d",frame+1);
          }
          MenuActivate1Arg(G, x, y, x, y, passive, "camera_motion",frame_str);
          goto done;
        }
        count++;
      }
      break;
    case cExecObject:
      if(ObjectGetSpecLevel(rec->obj,0)>=0) {
        draw_rect.top = rect->top - (height * count) / expected;
        draw_rect.bottom = rect->top - (height * (count + 1)) / expected;
        if((y>draw_rect.bottom) && (y<draw_rect.top)) {
          int n_frame = MovieGetLength(G);
          int frame = MovieXtoFrame(G, &draw_rect, n_frame, x, false);
          WordType frame_str = "0";
          if((frame>=0) && (frame<n_frame)) {
            sprintf(frame_str,"%d",frame+1);
          }
          MenuActivate2Arg(G, x, y, x, y, passive, "obj_motion", rec->obj->Name, frame_str);
          goto done;
        }
        count++;
      }
      break;
    }
  }
  }
 done:
  return;
}

void ExecutiveMotionClick(PyMOLGlobals * G, BlockRect *rect,int mode, int expected, int x, int y, int nearest)
{
  CExecutive *I = G->Executive;
  SpecRec *rec = NULL;
  BlockRect draw_rect = *rect;
  int count = 0;
  int height = rect->top - rect->bottom;
  while(ListIterate(I->Spec, rec, next)) {
    switch(rec->type) {
    case cExecAll:
      if(MovieGetSpecLevel(G,0)>=0) {
        draw_rect.top = rect->top - (height * count) / expected;
        draw_rect.bottom = rect->top - (height * (count + 1)) / expected;
        if((y>=draw_rect.bottom) && (y<=draw_rect.top)) {
          MoviePrepareDrag(G,&draw_rect,NULL,mode,x,y,nearest);
          goto done;
        }
        count++;
      }
      break;
    case cExecObject:
      if(ObjectGetSpecLevel(rec->obj,0)>=0) {
        MoviePrepareDrag(G,rect,NULL,mode,x,y,nearest);
        draw_rect.top = rect->top - (height * count) / expected;
        draw_rect.bottom = rect->top - (height * (count + 1)) / expected;
        if((y>=draw_rect.bottom) && (y<=draw_rect.top)) {
          MoviePrepareDrag(G,&draw_rect,rec->obj,mode,x,y,nearest);
          goto done;
        }
        count++;
      }
      break;
    }
  }
 done:
  return;
}

int ExecutiveReference(PyMOLGlobals * G, int action, const char *sele, int state, int quiet)
{
  int result = -1;
  int s1 = SelectorIndexByName(G, sele);
  if(s1 >= 0) {
    ObjectMoleculeOpRec op;
    ObjectMoleculeOpRecInit(&op);

    switch (action) {
    case 1:
      op.code = OMOP_ReferenceStore;
      break;
    case 2:
      op.code = OMOP_ReferenceRecall;
      break;
    case 3:
      op.code = OMOP_ReferenceValidate;
      break;
    case 4:
      op.code = OMOP_ReferenceSwap;
      break;
    }
    op.i1 = state;
    op.i2 = 0;
    ExecutiveObjMolSeleOp(G, s1, &op);
    result = op.i2;
  }
  return result;
}

/**
 * Does the object and state validation which is common to isosurface, isomesh,
 * and volume creation.
 */
template <typename ObjectT>
static pymol::Result<> EtcHelper(PyMOLGlobals* G,             //
    const char* obj_name, int& obj_state, ObjectT*& origObj,  //
    const char* map_name, int& map_state, ObjectMap*& mapObj, //
    int& multi)
{
  if (obj_state < -3) {
    return pymol::make_error("Invalid state ", obj_state + 1);
  }

  if (map_state < -3) {
    return pymol::make_error("Invalid source_state ", map_state + 1);
  }

  mapObj = ExecutiveFindObject<ObjectMap>(G, map_name);
  if (!mapObj) {
    return pymol::make_error("Map object \"", map_name, "\" not found");
  }

  // if object with same name but different type exists, overwrite it
  auto anyOrigObj = ExecutiveFindObjectByName(G, obj_name);
  origObj = dynamic_cast<ObjectT*>(anyOrigObj);
  if (anyOrigObj && !origObj) {
    ExecutiveDelete(G, obj_name);
  }

  switch (obj_state) {
  case -1:
    // all states
    obj_state = 0;
    map_state = -1;
    break;
  case -2:
    // current state
    obj_state = SceneGetState(G);
    if (map_state < 0)
      map_state = obj_state;
    break;
  case -3:
    // append mode
    obj_state = origObj ? origObj->getNFrame() : 0;
    if (map_state < 0)
      map_state = obj_state;
    break;
  }

  switch (map_state) {
  case -1:
    // all states
    map_state = 0;
    multi = true;
    break;
  case -2:
    // current state
    map_state = SceneGetState(G);
    break;
  case -3:
    // append mode
    map_state = mapObj->getNFrame() - 1;
    break;
  }

  return {};
}

pymol::Result<> ExecutiveIsosurfaceEtc(PyMOLGlobals * G,
                           const char *surf_name, const char *map_name, float lvl,
                           const char *sele, float fbuf, int state,
                           float carve, int map_state, int side,
                           int quiet, int surf_mode)
{
  int c;
  ObjectSurface *obj = nullptr, *origObj = nullptr;
  ObjectMap *mapObj;
  float mn[3] = { 0, 0, 0 };
  float mx[3] = { 15, 15, 15 };
  pymol::vla<float> vert_vla;
  ObjectMapState *ms;
  int multi = false;

  auto res0 = EtcHelper(
      G, surf_name, state, origObj, map_name, map_state, mapObj, multi);
  if (!res0) {
    return res0.error();
  }

  {
    while(1) {
      ms = ObjectMapStateGetActive(mapObj, map_state);
      if(ms) {
        int box_mode = (sele && sele[0]) ? 1 : 0;
        switch (box_mode) {
        case 0:                /* using map to get extents */
          for(c = 0; c < 3; c++) {
            mn[c] = ms->Corner[c];
            mx[c] = ms->Corner[3 * 7 + c];
          }
          if(!ms->Matrix.empty()) {
            transform44d3f(ms->Matrix.data(), mn, mn);
            transform44d3f(ms->Matrix.data(), mx, mx);
            {
              float tmp;
              int a;
              for(a = 0; a < 3; a++)
                if(mn[a] > mx[a]) {
                  tmp = mn[a];
                  mn[a] = mx[a];
                  mx[a] = tmp;
                }
            }
          }
          carve = 0.0F;
          break;
        case 1:                /* using selection to get extents */
          {
            auto tmpsele = SelectorTmp2::make(G, sele);
            p_return_if_error(tmpsele);
            auto s1 = tmpsele->getName();
            ExecutiveGetExtent(G, s1, mn, mx, false, -1, false);  /* TODO state */
            if(carve != 0.0F) {
              vert_vla = pymol::vla_take_ownership(
                  ExecutiveGetVertexVLA(G, s1, state));
              if(fbuf <= R_SMALL4)
                fbuf = fabs(carve);
            }
          }
          for(c = 0; c < 3; c++) {
            mn[c] -= fbuf;
            mx[c] += fbuf;
          }
          break;
        }
        PRINTFB(G, FB_CCmd, FB_Blather)
          " Isosurface: buffer %8.3f carve %8.3f\n", fbuf, carve ENDFB(G);
        obj =
          ObjectSurfaceFromBox(G, origObj, mapObj,
                                           map_state, state, mn, mx, lvl,
                                           static_cast<cIsosurfaceMode>(surf_mode),
                                           carve, std::move(vert_vla),
                                           static_cast<cIsosurfaceSide>(side),
                                           quiet);
        /* copy the map's TTT */
        ExecutiveMatrixCopy2(G, mapObj, obj, 1, 1, -1, -1, false, 0, quiet);

        if(!origObj) {
          ObjectSetName(obj, surf_name);
          ExecutiveManageObject(G, obj, -1, quiet);
        }
        if(SettingGetGlobal_b(G, cSetting_isomesh_auto_state))
          if(obj)
            ObjectGotoState(obj, state);
        if(!quiet) {
          PRINTFB(G, FB_ObjectSurface, FB_Actions)
            " Isosurface: created \"%s\", setting level to %5.3f\n", surf_name, lvl
            ENDFB(G);
        }
      } else if(!multi) {
        return pymol::make_error(
            "state ", map_state + 1, " not present in map \"", map_name, "\"");
      }
      if(multi) {
        origObj = obj;
        map_state++;
        state++;
        if(map_state >= mapObj->State.size())
          break;
      } else {
        break;
      }
    }
  }
  return {};
}

pymol::Result<> ExecutiveIsomeshEtc(PyMOLGlobals * G,
                        const char *mesh_name, const char *map_name, float lvl,
                        const char *sele, float fbuf, int state,
                        float carve, int map_state, int quiet,
                        int mesh_mode_, float alt_lvl)
{
  auto const mesh_mode = static_cast<cIsomeshMode>(mesh_mode_);
  ObjectMesh *obj = nullptr, *origObj = nullptr;
  ObjectMap *mapObj;
  float mn[3] = { 0, 0, 0 };
  float mx[3] = { 15, 15, 15 };
  float *vert_vla = NULL;
  int multi = false;
  ObjectMapState *ms;
  ObjectMolecule *sele_obj = NULL;
  CSymmetry *symm;

  auto res0 = EtcHelper(
      G, mesh_name, state, origObj, map_name, map_state, mapObj, multi);
  if (!res0) {
    return res0.error();
  }

  {
    while(1) {
      ms = ObjectMapStateGetActive(mapObj, map_state);
      if(ms) {
        int box_mode = (sele && sele[0]) ? 1 : 0;
        switch (box_mode) {
        case 0:                /* do the whole map */
          {
            int c;
            for(c = 0; c < 3; c++) {
              mn[c] = ms->Corner[c];
              mx[c] = ms->Corner[3 * 7 + c];
            }
          }
          if(!ms->Matrix.empty()) {
            transform44d3f(ms->Matrix.data(), mn, mn);
            transform44d3f(ms->Matrix.data(), mx, mx);
            {
              float tmp;
              int a;
              for(a = 0; a < 3; a++)
                if(mn[a] > mx[a]) {
                  tmp = mn[a];
                  mn[a] = mx[a];
                  mx[a] = tmp;
                }
            }
          }
          carve = -0.0;         /* impossible */
          break;
        case 1:                /* just do area around selection */
	  /* determine the selected object */
          {
            auto tmpsele = SelectorTmp2::make(G, sele);
            p_return_if_error(tmpsele);
            int sele1 = tmpsele->getIndex();
            if(sele1 >= 0)
              sele_obj = SelectorGetSingleObjectMolecule(G, sele1);
            auto s1 = tmpsele->getName();
            ExecutiveGetExtent(G, s1, mn, mx, false, -1, false);  /* TODO state */
            if(carve != 0.0) {
              vert_vla = ExecutiveGetVertexVLA(G, s1, state);
              if(fbuf <= R_SMALL4)
                fbuf = fabs(carve);
            }
          }
          {
            int c;
            for(c = 0; c < 3; c++) {
              mn[c] -= fbuf;
              mx[c] += fbuf;
            }
          }
          break;
        }
        PRINTFB(G, FB_CCmd, FB_Blather)
          " Isomesh: buffer %8.3f carve %8.3f \n", fbuf, carve ENDFB(G);

        symm = NULL;
        if(sele_obj &&  ObjectMapValidXtal(mapObj, state)) {
          if(SettingGet_b(G, NULL, sele_obj->Setting.get(), cSetting_map_auto_expand_sym)
              && (sele_obj->Symmetry)) {
            // legacy default: take symmetry from molecular object
            symm = sele_obj->Symmetry.get();
          } else if(SettingGet_b(G, NULL, mapObj->Setting.get(), cSetting_map_auto_expand_sym)) {
            // fallback: take symmetry from map state
            symm = ms->Symmetry.get();
          }
        }

        if(symm) {
          obj = ObjectMeshFromXtalSym(G, origObj, mapObj,
                                                  symm,
                                                  map_state, state, mn, mx, lvl,
                                                  mesh_mode, carve, vert_vla, alt_lvl,
                                                  quiet);
        } else {
          obj = NULL;
        }
        if(!obj) {
          obj = ObjectMeshFromBox(G, origObj, mapObj,
                                              map_state, state, mn, mx, lvl, mesh_mode,
                                              carve, vert_vla, alt_lvl, quiet);
        }
        /* copy the map's TTT */
        ExecutiveMatrixCopy2(G, mapObj, obj, 1, 1, -1, -1, false, 0, quiet);

        if(!origObj) {
          ObjectSetName(obj, mesh_name);
          ExecutiveManageObject(G, obj, false, quiet);
        }

        if(SettingGetGlobal_b(G, cSetting_isomesh_auto_state))
          if(obj)
            ObjectGotoState(obj, state);
        if(!quiet) {
          if (mesh_mode != cIsomeshMode::gradient) {
            PRINTFB(G, FB_ObjectMesh, FB_Actions)
              " Isomesh: created \"%s\", setting level to %5.3f\n", mesh_name, lvl
              ENDFB(G);
          } else {
            PRINTFB(G, FB_ObjectMesh, FB_Actions)
              " Gradient: created \"%s\"\n", mesh_name ENDFB(G);
          }
        }
      } else if(!multi) {
        return pymol::make_error(
            "state ", map_state + 1, " not present in map \"", map_name, "\"");
      }
      if(multi) {
        origObj = obj;
        map_state++;
        state++;
        if(map_state >= mapObj->State.size())
          break;
      } else {
        break;
      }
    }
  }
  return {};
}

pymol::Result<>
ExecutiveVolume(PyMOLGlobals * G, const char *volume_name, const char *map_name,
                        float lvl,
                        const char *sele, float fbuf, int state,
                        float carve, int map_state, int quiet)
{
  ObjectVolume *obj = nullptr, *origObj = nullptr;
  ObjectMap *mapObj;
  float mn[3] = { 0, 0, 0 };
  float mx[3] = { 15, 15, 15 };
  float *vert_vla = NULL;
  int multi = false;
  ObjectMapState *ms;
  ObjectMolecule *sele_obj = NULL;
  CSymmetry *symm;

  /* Goal: make a new volume map from the object or selection
   *
   * If the user specifies an VOLUME OBJECT NAME that already exists, then
   * check to make sure it's a volume.  If it is, keep it, else delete it
   * 
   * If the user specifies a MAP OBJECT NAME that already exists, then
   * make sure it's really a map.  Keep it if it is (to append the state)
   * otherwise, delete it.
   *
   * Given that the MAP is valid, find the current states for the MAP object
   * and the scene/protein.
   *
   * If the user is loading a map with N-states, into an N-state object, we need
   * to loop over each state.  So,
   *
   * FOR EACH state in the MAP:
   *   - determine the extents of the MAP for THIS STATE
   *
   *   - create the VOLUME object from the MAP/XTAL SYMM
   *   - if VOLUME CREATION didn't work, then try creating from a BOX
   *
   *   - transform (translate/rotate) the VOLUME vis so that it
   *     matches the map's ObjectMatrix.
   *
   *   - if this is a new VOLUME, (eg., origObj==NULL) then set the internal
   *     PyMOL state to manage the object and add the name to the object list
   *
   *   - update the map state for the next loop and re-iterate
   */

  auto res0 = EtcHelper(
      G, volume_name, state, origObj, map_name, map_state, mapObj, multi);
  if (!res0) {
    return res0.error();
  }

  {
    /* do for each state */
    while(1) {
      ms = ObjectMapStateGetActive(mapObj, map_state);
      if(ms) {
	/* determine extents */
        int box_mode = (sele && sele[0]) ? 1 : 0;
        switch (box_mode) {
        case 0:                /* do the whole map */
          {
            int c;
            for(c = 0; c < 3; c++) {
              mn[c] = ms->Corner[c];
              mx[c] = ms->Corner[3 * 7 + c];
            }
          }
          if(!ms->Matrix.empty()) {
            transform44d3f(ms->Matrix.data(), mn, mn);
            transform44d3f(ms->Matrix.data(), mx, mx);
            {
              float tmp;
              int a;
              for(a = 0; a < 3; a++)
                if(mn[a] > mx[a]) {
                  tmp = mn[a];
                  mn[a] = mx[a];
                  mx[a] = tmp;
                }
            }
          }
          carve = -0.0;         /* impossible */
          break;
        case 1:                /* just do area around selection */
          {
            auto tmpsele = SelectorTmp2::make(G, sele);
            p_return_if_error(tmpsele);
            int sele1 = tmpsele->getIndex();
            if(sele1 >= 0)
              sele_obj = SelectorGetSingleObjectMolecule(G, sele1);
            auto const s1 = tmpsele->getName();
            ExecutiveGetExtent(G, s1, mn, mx, false, -1, false);  /* TODO state */
            if(carve != 0.0) {
              vert_vla = ExecutiveGetVertexVLA(G, s1, state);
              if(fbuf <= R_SMALL4)
                fbuf = fabs(carve);
            }
          }
          {
            int c;
            for(c = 0; c < 3; c++) {
              mn[c] -= fbuf;
              mx[c] += fbuf;
            }
          }
          break;
        }
        PRINTFB(G, FB_CCmd, FB_Blather)
          " Volume: buffer %8.3f carve %8.3f \n", fbuf, carve ENDFB(G);

        symm = NULL;
        if(sele_obj && ObjectMapValidXtal(mapObj, state)) {
          if(SettingGet_b(G, NULL, sele_obj->Setting.get(), cSetting_map_auto_expand_sym)
              && (sele_obj->Symmetry)) {
            // legacy default: take symmetry from molecular object
            symm = sele_obj->Symmetry.get();
          } else if(SettingGet_b(G, NULL, mapObj->Setting.get(), cSetting_map_auto_expand_sym)) {
            // fallback: take symmetry from map state
            symm = ms->Symmetry.get();
          }
        }

        if(symm) {
          obj = ObjectVolumeFromXtalSym(G, origObj, mapObj,
                                                  symm,
                                                  map_state, state, mn, mx, lvl,
                                                  box_mode, carve, vert_vla,
                                                  quiet);
        } else {
          obj = NULL;
        }

        if(!obj) {
          obj = ObjectVolumeFromBox(G, origObj, mapObj,
                                              map_state, state, mn, mx, lvl, box_mode,
                                              carve, vert_vla, quiet);
        }
        /* copy the map's TTT */
        ExecutiveMatrixCopy2(G, mapObj, obj, 1, 1, -1, -1, false, 0, quiet);
	/* set the object name
	 * manage the object in the UI */
        if(!origObj) {
          ObjectSetName(obj, volume_name);
          ExecutiveManageObject(G, obj, false, quiet);
        }

        if(SettingGetGlobal_b(G, cSetting_isomesh_auto_state))
          if(obj)
            ObjectGotoState(obj, state);
        if(!quiet) {
	  PRINTFB(G, FB_ObjectVolume, FB_Actions)
	    " Volume: created \"%s\"\n", volume_name
	    ENDFB(G);
        }
      } else if(!multi) {
        return pymol::make_error(
            "state ", map_state + 1, " not present in map \"", map_name, "\"");
      }
      if(multi) {
        origObj = obj;
        map_state++;
        state++;
        if(map_state >= mapObj->State.size())
          break;
      } else {
        break;
      }
    }
  }
  return {};
}

std::string ExecutivePreparePseudoatomName(
    PyMOLGlobals* G, pymol::zstring_view object_name)
{
  std::string new_object_name;
  if (object_name.empty()) {
    new_object_name = ExecutiveGetUnusedName(G, "pseudo");
  } else {
    ObjectNameType valid_name{};
    assert(object_name.size() < sizeof(ObjectNameType));
    std::copy_n(object_name.c_str(), object_name.size(), valid_name);
    ObjectMakeValidName(G, valid_name);
    new_object_name = valid_name;
  }
  return new_object_name;
}

pymol::Result<> ExecutivePseudoatom(PyMOLGlobals* G, pymol::zstring_view object_name_view,
    const char* sele, const char* name, const char* resn, const char* resi,
    const char* chain, const char* segi, const char* elem, float vdw,
    int hetatm, float b, float q, const char* label, const float* pos, int color,
    int state, int mode, int quiet)
{
  pymol::Result<SelectorTmp> s1;

  auto object_name = object_name_view.c_str();
  auto obj = ExecutiveFindObject<ObjectMolecule>(G, object_name);

  int is_new = false;
  int sele_index = -1;
  float local_pos[3];

  if(sele && sele[0]) {
    if(WordMatchExact(G, cKeywordCenter, sele, true)) {
      SceneGetCenter(G, local_pos);
      pos = local_pos;
    } else if(WordMatchExact(G, cKeywordOrigin, sele, true)) {
      SceneOriginGet(G, local_pos);
      pos = local_pos;
    } else {
      s1 = SelectorTmp::make(G, sele);
      p_return_if_error(s1);
      sele_index = s1->getIndex();
      assert(sele_index >= 0);
    }
  }
  if(!obj) {
    /* new object */
    is_new = true;
    obj = new ObjectMolecule(G, false);
    ObjectSetName(obj, object_name);
  }

    if(ObjectMoleculeAddPseudoatom(obj, sele_index, name, resn, resi, chain,
                                   segi, elem, vdw, hetatm, b, q, label, pos, color,
                                   state, mode, quiet)) {
      if(is_new) {
        ExecutiveDelete(G, object_name);        /* just in case */
        ExecutiveManageObject(G, obj, false, true);
      } else {
        ExecutiveUpdateObjectSelection(G, obj);
      }
    }
  return {};
}

static void ExecutiveInvalidateGridSlots(PyMOLGlobals * G)
{
  CExecutive *I = G->Executive;
  I->ValidGridSlots = false;
}

static void ExecutiveSetGridSlot(SpecRec *rec, int new_grid_slot){
  if (rec->grid_slot != new_grid_slot){
    CGOFree(rec->gridSlotSelIndicatorsCGO);
    rec->gridSlotSelIndicatorsCGO = NULL;
    rec->grid_slot = new_grid_slot;
  }
}
static void ExecutiveUpdateGridSlots(PyMOLGlobals * G, int force)
{
  CExecutive *I = G->Executive;
  int grid_slot_count = 0;
  int grid_by_group = 1;        /* grid slots are inherited this many levels */

  ExecutiveUpdateGroups(G, false);
  if(force || (!I->ValidGridSlots)) {
    CTracker *I_Tracker = I->Tracker;
    I->ValidGridSlots = true;
    {
      SpecRec *rec = NULL;
      while(ListIterate(I->Spec, rec, next)) {
	ExecutiveSetGridSlot(rec, 0);
        if(rec->type == cExecObject) {
          /* make sure every object (potentially) needing a grid slot gets one */
          switch (rec->obj->type) {
          case cObjectMolecule:
          case cObjectMap:
          case cObjectMesh:
          case cObjectMeasurement:
          case cObjectCallback:
          case cObjectCGO:
          case cObjectSurface:
          case cObjectSlice:
          case cObjectGadget:
          case cObjectGroup:
          case cObjectVolume:
	    ExecutiveSetGridSlot(rec, ++grid_slot_count);
            break;
          }
        }
      }
    }

    if(grid_by_group) {
      SpecRec *rec = NULL, *group_rec = NULL;
      while(ListIterate(I->Spec, rec, next)) {
        OVreturn_word result;
        if(OVreturn_IS_OK
           ((result = OVLexicon_BorrowFromCString(I->Lex, rec->group_name)))) {
          if(OVreturn_IS_OK((result = OVOneToOne_GetForward(I->Key, result.word)))) {
            if(TrackerGetCandRef(I_Tracker, result.word,
                                 (TrackerRef **) (void *) &group_rec)) {
              int grid_slot_group_depth = grid_by_group;
              {
                SpecRec *check_rec = group_rec;
                while(check_rec && grid_slot_group_depth) {
                  if(grid_slot_group_depth == 1)
		    ExecutiveSetGridSlot(rec, check_rec->grid_slot);
                  if(check_rec == rec) {        /* cycle */
                    break;
                  } else {
                    check_rec = check_rec->group;
                    grid_slot_group_depth--;
                  }
                }
              }
            }
          }
        }
      }
    }

    {
      SpecRec *rec = NULL;
      while(ListIterate(I->Spec, rec, next)) {
        if(rec->type == cExecObject) {
          int obj_slot = SettingGet_i(G, rec->obj->Setting.get(), NULL, cSetting_grid_slot);
          if(obj_slot == -1) {
            rec->obj->grid_slot = rec->grid_slot;
          } else
            rec->obj->grid_slot = obj_slot;
        }
      }
    }
  }
}

void ExecutiveInvalidatePanelList(PyMOLGlobals * G)
{
  CExecutive *I = G->Executive;
  I->Panel.clear();
  ExecutiveInvalidateGridSlots(G);
}


/**
 * Add all members which belong to group
 */
static void PanelListGroup(
    CExecutive* I, SpecRec const* group, unsigned level, bool hide_underscore)
{
  for (auto& rec : pymol::make_list_adapter(I->Spec)) {
    if (rec.group != group) {
      continue;
    }

    assert(!rec.in_panel);

    if (rec.isHiddenNotRecursive(hide_underscore)) {
      continue;
    }

    if (!level) {
      assert(!rec.group_name[0]);
      rec.group_name[0] = 0;     /* force open any cycles which have been created... */
    }

    I->Panel.emplace_back(&rec, level);
    rec.in_panel = true;

    if (auto obj_group = dynamic_cast<ObjectGroup*>(rec.obj)) {
      auto& panelitem = I->Panel.back();
      panelitem.is_group = true;
      if (obj_group->OpenOrClosed) {
        panelitem.is_open = true;
        PanelListGroup(I, &rec, level + 1, hide_underscore);
      }
    }
  }
}

static void ExecutiveUpdatePanelList(PyMOLGlobals * G)
{
  CExecutive *I = G->Executive;
  int hide_underscore = SettingGetGlobal_b(G, cSetting_hide_underscore_names);
  if (I->Panel.empty()) {
    for (auto& rec : pymol::make_list_adapter(I->Spec)) {
      rec.in_panel = false;
    }

    /* brute-force & inefficient -- need to optimize algorithm */
    PanelListGroup(I, nullptr, 0, hide_underscore);
  }
}

void ExecutiveInvalidateSceneMembers(PyMOLGlobals * G)
{
  CExecutive *I = G->Executive;
  I->ValidSceneMembers = false;
}

void ExecutiveUpdateSceneMembers(PyMOLGlobals * G)
{
  CExecutive *I = G->Executive;
  ExecutiveUpdateGroups(G, false);
  ExecutiveUpdateGridSlots(G, false);
  if(!I->ValidSceneMembers) {
    SpecRec *rec = NULL;
    while(ListIterate(I->Spec, rec, next)) {
      if(rec->type == cExecObject) {
        int visible = rec->visible;
        SpecRec *group_rec = rec->group;
        while(visible && group_rec) {   /* visibility is a group issue... */
          if(!group_rec->visible)
            visible = false;
          else
            group_rec = group_rec->group;
        }
        if(rec->in_scene && !visible) {
          rec->in_scene = SceneObjectDel(G, rec->obj, true);
        } else if(visible && !rec->in_scene) {
          rec->in_scene = SceneObjectAdd(G, rec->obj);
        }
      }
    }
    I->ValidSceneMembers = true;
  }
}

void ExecutiveInvalidateGroups(PyMOLGlobals * G, bool force)
{
  auto I = G->Executive;
  if(!(force || I->ValidGroups)) {
    return;
  }
  for (auto& rec : pymol::make_list_adapter(G->Executive->Spec)) {
    rec.group = nullptr;
    if (ExecutiveIsObjectType(rec, cObjectGroup)) {
      if (rec.group_member_list_id) {
        TrackerDelList(I->Tracker, rec.group_member_list_id);
      }
      rec.group_member_list_id = 0;        /* not a list */
    }
  }
  I->ValidGroups = false;
  ExecutiveInvalidateSceneMembers(G);
  ExecutiveInvalidatePanelList(G);
  /* any changes to group structure means that we need to check scene
     members */
}

void ExecutiveUpdateGroups(PyMOLGlobals * G, bool force)
{
  CExecutive *I = G->Executive;

  if(force || (!I->ValidGroups)) {
    CTracker *I_Tracker = I->Tracker;

    /* first, get rid of existing group lists */

    if(force || I->ValidGroups)
      ExecutiveInvalidateGroups(G, true);

    /* create empty lists for each group (also init grid_slot) */

    for (auto& rec : pymol::make_list_adapter(G->Executive->Spec)) {
      rec.group = nullptr;
      if (ExecutiveIsObjectType(rec, cObjectGroup)) {
          rec.group_member_list_id = TrackerNewList(I_Tracker, nullptr);
      }
    }

    /* iterate through and populate groups lists with their members */

    for (auto& rec : pymol::make_list_adapter(G->Executive->Spec)) {
      OVreturn_word result;
      if(OVreturn_IS_OK
         ((result = OVLexicon_BorrowFromCString(I->Lex, rec.group_name)))) {
        if(OVreturn_IS_OK((result = OVOneToOne_GetForward(I->Key, result.word)))) {
          SpecRec* group_rec = nullptr;
          if(TrackerGetCandRef(I_Tracker, result.word,
                                 (TrackerRef **) (void *) &group_rec)) {
            int cycle = false;
            {                 /* don't close infinite loops */
              SpecRec *check_rec = group_rec;
              while(check_rec) {
                if(check_rec == &rec) {
                  cycle = true;
                  break;
                } else {
                  check_rec = check_rec->group;
                }
              }
            }
            if(!cycle) {
              rec.group = group_rec;
              TrackerLink(I_Tracker, rec.cand_id, group_rec->group_member_list_id, 1);
            }
          }
        }
      }
    }

    /* note that it is possible to have infinite loops -- these must be
       allowed for later in the group expansion routine(s) */
    I->ValidGroups = true;
    ExecutiveInvalidatePanelList(G);
  }
}

static int ExecutiveGetObjectParentList(PyMOLGlobals * G, SpecRec * child)
{
  int list_id = 0;
  ExecutiveUpdateGroups(G, false);
  {
    CExecutive *I = G->Executive;
    CTracker *I_Tracker = I->Tracker;
    int priority = 1;           /* generations removed from child */
    int repeat_flag = true;
    SpecRec *group_rec = NULL;

    list_id = TrackerNewList(I_Tracker, NULL);
    while(child && child->group && repeat_flag) {
      OVreturn_word result;
      repeat_flag = false;
      if(OVreturn_IS_OK
         ((result = OVLexicon_BorrowFromCString(I->Lex, child->group_name)))) {
        if(OVreturn_IS_OK((result = OVOneToOne_GetForward(I->Key, result.word)))) {
          if(TrackerGetCandRef(I_Tracker, result.word,
                               (TrackerRef **) (void *) &group_rec)) {
            if(TrackerLink(I_Tracker, result.word, list_id, priority++)) {
              /* checking this prevents infinite loops */
              if(group_rec->group) {
                repeat_flag = true;
                child = group_rec;
              }
            }
          }
        }
      }
    }
  }
  return list_id;
}

pymol::TrackerAdapter<SpecRec> ExecutiveGetSpecRecParents(
    PyMOLGlobals* G, SpecRec& rec)
{
  return pymol::TrackerAdapter<SpecRec>(
      G->Executive->Tracker, ExecutiveGetObjectParentList(G, &rec));
}

int ExecutiveVdwFit(PyMOLGlobals * G, const char *s1, int state1, const char *s2, int state2,
                    float buffer, int quiet)
{
  SelectorTmp tmpsele1(G, s1);
  SelectorTmp tmpsele2(G, s2);
  int sele1 = tmpsele1.getIndex();
  int sele2 = tmpsele2.getIndex();

  int ok = true;

  if((sele1 >= 0) && (sele2 >= 0)) {
    ok = SelectorVdwFit(G, sele1, state1, sele2, state2, buffer, quiet);
  } else {
    ok = false;
  }
  return ok;
}

static int get_op_cnt(PyMOLGlobals * G)
{
  int result = 5;
  if(!strcmp(SettingGetGlobal_s(G, cSetting_button_mode_name), "3-Button Motions"))
    result = 6;
  return result;
}

static int ExecutiveAddKey(CExecutive * I, SpecRec * rec)
{
  int ok = false;
  OVreturn_word result;
  if(OVreturn_IS_OK((result = OVLexicon_GetFromCString(I->Lex, rec->name)))) {
    if(OVreturn_IS_OK(OVOneToOne_Set(I->Key, result.word, rec->cand_id))) {
      ok = true;
    }
  }
  return ok;
}

static int ExecutiveDelKey(CExecutive * I, SpecRec * rec)
{
  int ok = false;
  OVreturn_word result;
  if(OVreturn_IS_OK((result = OVLexicon_BorrowFromCString(I->Lex, rec->name)))) {
    if(OVreturn_IS_OK(OVLexicon_DecRef(I->Lex, result.word)) &&
       OVreturn_IS_OK(OVOneToOne_DelForward(I->Key, result.word))) {
      ok = true;
    }
  }
  return ok;
}

static SpecRec *ExecutiveUnambiguousNameMatch(PyMOLGlobals * G, const char *name)
{
  CExecutive *I = G->Executive;
  SpecRec *result = NULL;
  SpecRec *rec = NULL;
  int best = 0;
  int wm;
  int ignore_case = SettingGetGlobal_b(G, cSetting_ignore_case);

  while(ListIterate(I->Spec, rec, next)) {
    wm = WordMatch(G, name, rec->name, ignore_case);
    if(wm < 0) {                /* exact match, so this is valid */
      result = rec;
      best = wm;
      break;
    } else if((wm > 0) && (best < wm)) {
      result = rec;
      best = wm;
    } else if((wm > 0) && (best == wm)) {       /* ambiguous match... no good */
      result = NULL;
    }
  }
  return (result);
}

static SpecRec *ExecutiveAnyCaseNameMatch(PyMOLGlobals * G, const char *name)
{
  CExecutive *I = G->Executive;
  SpecRec *result = NULL;
  SpecRec *rec = NULL;

  int ignore_case = SettingGetGlobal_b(G, cSetting_ignore_case);
  while(ListIterate(I->Spec, rec, next)) {
    if(WordMatchExact(G, name, rec->name, ignore_case)) {
      result = rec;
      break;
    }
  }
  return (result);
}

/**
 * Scroll the i'th match in the object menu panel to the top.
 * Scroll to last match if i < 0 and to first match if i > #-1.
 * Open groups if hit is inside.
 * Highlight the hit (same as mouse click highlight).
 *
 * Returns the number of hits
 */
int ExecutiveScrollTo(PyMOLGlobals * G, const char * name, int i) {
  CExecutive *I = G->Executive;
  int pos = 0, numhits = 0;
  ObjectGroup *group;
  SpecRec *tmp, *spec = NULL, *first = NULL;
  int ignore_case = SettingGetGlobal_b(G, cSetting_ignore_case);
  int j, lendiff, plen = strlen(name);

  ok_assert(1, I->Spec);

  // i'th substring match, skip the "all" item
  for(tmp = I->Spec->next; tmp; tmp = tmp->next) {
    lendiff = strlen(tmp->name) - plen;
    for(j = 0; j <= lendiff; j++)
      if(WordMatchNoWild(G, name, tmp->name + j, ignore_case)) {
        if(numhits++ == i || i < 0)
          spec = tmp;
        if(!first)
          first = tmp;
        break;
      }
    tmp->hilight = 0;
  }

  // if i was out of range
  if(!spec)
    spec = first;

  ok_assert(1, spec);

  // flash button until panel is clicked for the next time
  spec->hilight = 1;

  // open parent groups
  for(tmp = spec->group; tmp; tmp = tmp->group) {
    if(!(tmp->type == cExecObject &&
         tmp->obj->type == cObjectGroup))
      break;
    group = (ObjectGroup *) tmp->obj;
    if(!group->OpenOrClosed) {
      group->OpenOrClosed = 1;
      ExecutiveInvalidatePanelList(G);
    }
  }

  // in case any parent got opened
  ExecutiveUpdatePanelList(G);

  // scroll that record to the top
  for (auto const& panelitem : I->Panel) {
    if (panelitem.spec == spec) {
      I->m_ScrollBar.setValueNoCheck(pos);
      return numhits;
    }
    pos++;
  }

ok_except1:
  return numhits;
}

void ExecutiveUpdateColorDepends(PyMOLGlobals * G, ObjectMolecule * mol)
{
  CExecutive *I = G->Executive;
  SpecRec *rec = NULL;

  while(ListIterate(I->Spec, rec, next)) {
    if(rec->type == cExecObject) {
      if(rec->obj->type == cObjectGadget) {
        ObjectGadget *gadget = (ObjectGadget *) rec->obj;
        if(gadget->GadgetType == cGadgetRamp) {
          ObjectGadgetRamp *ramp = (ObjectGadgetRamp *) gadget;
          if(ramp->RampType == cRampMol) {
            if(ramp->Mol == mol) {
              ExecutiveInvalidateRep(G, cKeywordAll, cRepAll, cRepInvColor);
              break;
            }
          }
        }
      }
    }
  }
}

void ExecutiveUpdateCoordDepends(PyMOLGlobals * G, ObjectMolecule * mol)
{                               /* nasty, ugly, inefficient hack */

  CExecutive *I = G->Executive;
  SpecRec *rec = NULL;
  ObjectGadget *gadget;
  int done_inv_all = false;
  int dynamic_measures = SettingGet_b(G, mol ? mol->Setting.get() : NULL, NULL,
      cSetting_dynamic_measures);

  while(ListIterate(I->Spec, rec, next)) {
    if(rec->type == cExecObject) {
      switch(rec->obj->type) {
      case cObjectGadget:
        if(done_inv_all)
          break;
        gadget = (ObjectGadget *) rec->obj;
        if(gadget->GadgetType == cGadgetRamp) {
          ObjectGadgetRamp *ramp = (ObjectGadgetRamp *) gadget;
          if(ramp->RampType == cRampMol) {
            if(ramp->Mol == mol) {
              ExecutiveInvalidateRep(G, cKeywordAll, cRepAll, cRepInvColor);
              done_inv_all = true;
            }
          }
        }
        break;
      case cObjectMeasurement:
        if(dynamic_measures)
          ObjectDistMoveWithObject((ObjectDist*) rec->obj, mol);
        break;
      case cObjectAlignment:
        rec->obj->invalidate(
            cRepAll, cRepInvRep, cSelectorUpdateTableAllStates);
        break;
      }
    }
  }
}

int ExecutiveValidNamePattern(PyMOLGlobals * G, const char *name)
{
  int result = false;
  CWordMatcher *matcher;
  CWordMatchOptions options;
  const char *wildcard = SettingGetGlobal_s(G, cSetting_wildcard);

  WordMatchOptionsConfigNameList(&options,
                                 *wildcard, SettingGetGlobal_b(G, cSetting_ignore_case));
  matcher = WordMatcherNew(G, name, &options, false);
  if(matcher) {                 /* this appears to be a pattern */
    result = true;
    WordMatcherFree(matcher);
  } else if(ExecutiveUnambiguousNameMatch(G, name)) {
    /* this does not appear to be a pattern, so it is an unambiguous partial name? */
    result = true;
  }
  return result;

}

#define cExecNoExpand false
#define cExecExpandGroups true
#define cExecExpandKeepGroups 2

static void ExecutiveExpandGroupsInList(PyMOLGlobals * G, int list_id, int expand_groups)
{
  CExecutive *I = G->Executive;
  CTracker *I_Tracker = I->Tracker;
  int new_member_added = true;
  SpecRec *rec;
  ExecutiveUpdateGroups(G, false);
  while(new_member_added) {     /* keep adding til we can't add no more */
    int iter_id = TrackerNewIter(I_Tracker, 0, list_id);
    int cand_id;
    new_member_added = false;
    if(iter_id) {
      while((cand_id = TrackerIterNextCandInList(I_Tracker, iter_id,
                                                 (TrackerRef **) (void *) &rec))) {
        if(rec && (rec->type == cExecObject) &&
           rec->group_member_list_id && (rec->obj->type == cObjectGroup)) {
          int group_iter_id = TrackerNewIter(I_Tracker, 0, rec->group_member_list_id);
          int group_cand_id;
          SpecRec *group_rec;
          if(group_iter_id) {
            while((group_cand_id = TrackerIterNextCandInList(I_Tracker, group_iter_id,
                                                             (TrackerRef **) (void *)
                                                             &group_rec))) {
              if(group_rec && group_cand_id) {
                if(TrackerLink(I_Tracker, group_cand_id, list_id, 1))
                  new_member_added = true;
              }
            }
            TrackerDelIter(I_Tracker, group_iter_id);
          }
        }
      }
      TrackerDelIter(I_Tracker, iter_id);
    }
  }
  /* now purge all groups from the expanded list */
  if(expand_groups != cExecExpandKeepGroups) {
    int iter_id = TrackerNewIter(I_Tracker, 0, list_id);
    int cand_id;
    while((cand_id = TrackerIterNextCandInList(I_Tracker, iter_id,
                                               (TrackerRef **) (void *) &rec))) {
      if(rec && (rec->type == cExecObject) && (rec->obj->type == cObjectGroup)) {
        TrackerUnlink(I_Tracker, cand_id, list_id);
      }
    }
  }
}


/* DON'T FORGET TO RELEASE LIST WHEN DONE!!! */
int ExecutiveGetNamesListFromPattern(PyMOLGlobals * G, const char *name,
                                            int allow_partial, int expand_groups)
{
  CExecutive *I = G->Executive;
  int result = 0;
  CWordMatcher *matcher;
  CWordMatchOptions options;
  CTracker *I_Tracker = I->Tracker;
  const char *wildcard = SettingGetGlobal_s(G, cSetting_wildcard);
  int iter_id = TrackerNewIter(I_Tracker, 0, I->all_names_list_id);
  int cand_id;
  int group_found = false;
  SpecRec *rec;

  if (!name)
    return -1;

  // sanity check: name patterns are not object selections, bail if
  // parenthesis or operators in pattern
  if (strchr(name, '(') || strchr(name, ')') || strchr(name, '|')) {
    PRINTFB(G, FB_Executive, FB_Errors)
      " Names-Pattern-Error: Pattern looks like an atom selection"
      " (has parenthesis or operators), this is not supported for"
      " object name patterns.\n" ENDFB(G);
    return -1;
  }

  // special case: allow "not ..."
  bool match_not = false;
  if (WordMatchNoWild(G, "not ", name, false)) {
    match_not = true;
    name += 4;
  } else if (name[0] == '!') {
    match_not = true;
    name += 1;
  }

  // skip whitespace
  while (name[0] == ' ') {
    ++name;
  }

  bool match_enabled = WordMatchExact(G, "enabled", name, false);

  // ignore % and ? prefixes
  while(name[0] && (name[0] == '%' || name[0] == '?'))
    name++;

  WordMatchOptionsConfigNameList(&options,
                                 *wildcard, SettingGetGlobal_b(G, cSetting_ignore_case));
  matcher = WordMatcherNew(G, name, &options, /* force= */ match_not);
  if(matcher || match_enabled) {
    if(iter_id) {
      while((cand_id = TrackerIterNextCandInList(I_Tracker, iter_id,
                                                 (TrackerRef **) (void *) &rec))) {
        if(rec && !(rec->type == cExecAll)) {
          bool test = match_enabled ? SpecRecIsEnabled(rec) :
            WordMatcherMatchAlpha(matcher, rec->name);
          if(test ^ match_not) {
            if((rec->type == cExecObject) && (rec->obj->type == cObjectGroup))
              group_found = true;
            if(!result)
              result = TrackerNewList(I_Tracker, NULL);
            if(result) {
              TrackerLink(I_Tracker, cand_id, result, 1);
            }
          }
        }
      }
    }
  } else if((rec = ExecutiveFindSpec(G, name))) {       /* only one name in list */
    if((rec->type == cExecObject) && (rec->obj->type == cObjectGroup))
      group_found = true;
    result = TrackerNewList(I_Tracker, NULL);
    TrackerLink(I_Tracker, rec->cand_id, result, 1);
  } else if(allow_partial && (rec = ExecutiveUnambiguousNameMatch(G, name))) {
    if((rec->type == cExecObject) && (rec->obj->type == cObjectGroup))
      group_found = true;
    result = TrackerNewList(I_Tracker, NULL);
    TrackerLink(I_Tracker, rec->cand_id, result, 1);
  }
  if(matcher)
    WordMatcherFree(matcher);
  if(iter_id)
    TrackerDelIter(I->Tracker, iter_id);
  if(group_found && expand_groups) {
    ExecutiveExpandGroupsInList(G, result, expand_groups);
  }
  return result;
}

static int ExecutiveUngroup(PyMOLGlobals* G, const char* members, int quiet)
{
  for (auto& rec : ExecutiveGetSpecRecsFromPattern(G, members)) {
    rec.group_name[0] = 0;
    rec.group = nullptr;
  }
  ExecutiveInvalidateGroups(G, true);
  return true;
}

/**
 * Removes recs from a group rec
 *
 * @param[in] rec group rec to remove recs from
 * @param[out] discarded discarded recs
 *
 * Note: Caller is responsible for calling ExecutivePurgeSpec on discarded recs
 */

static void ExecutiveGroupPurge(PyMOLGlobals* G, SpecRec* rec, std::vector<DiscardedRec>* const discarded = nullptr)
{
  std::vector<DiscardedRec> recs;
  auto ignore_case = SettingGet<bool>(G, cSetting_ignore_case);
  SpecRec *rec2 = NULL;
  while (ListIterate(G->Executive->Spec, rec2, next)) {
    if ((rec2->group == rec)
       || WordMatchExact(G, rec2->group_name, rec->name, ignore_case)) {
      bool save = discarded != nullptr;
      auto result = ExecutiveDelete(G, rec2->name, save);
      if (discarded && result) {
        discarded->insert(discarded->end(), result->begin(), result->end());
      }
      rec2 = nullptr;
    }
  }
}

int ExecutiveGroup(PyMOLGlobals * G, pymol::zstring_view nameView,
    pymol::zstring_view membersView, int action, int quiet)
{
  auto name = nameView.c_str();
  auto members = membersView.c_str();
  if (action == cExecutiveGroupUngroup) {
    // Up to PyMOL 2.3 the member argument was ignored and 'name' used for ungrouping
    if (name[0]) {
      members = name;
    }
    return ExecutiveUngroup(G, members, quiet);
  }

  int ok = true;
  CExecutive *I = G->Executive;
  auto ignore_case = SettingGet<bool>(G, cSetting_ignore_case);

  ObjectNameType valid_name;
  UtilNCopy(valid_name, name, sizeof(ObjectNameType));
  ObjectMakeValidName(G, valid_name);

  pymol::CObject *obj = ExecutiveFindObjectByName(G, valid_name);

  if(obj && (obj->type != cObjectGroup)) {
      PRINTFB(G, FB_Executive, FB_Errors)
        " Group-Error: object '%s' is not a group object.", name ENDFB(G);
      ok = false;
  } else {
    if((!obj) && (action == cExecutiveGroupAdd)) {
      obj = new ObjectGroup(G);
      if(obj) {
        ObjectSetName(obj, valid_name);
        ExecutiveManageObject(G, obj, false, true);
      }
    }
  }
  if((!members[0]) && ((action == cExecutiveGroupOpen) ||
                       (action == cExecutiveGroupClose) ||
                       (action == cExecutiveGroupToggle) ||
                       (action == cExecutiveGroupEmpty) ||
                       (action == cExecutiveGroupPurge) ||
                       (action == cExecutiveGroupExcise))) {
    ExecutiveUpdateGroups(G, false);
    {
      CTracker *I_Tracker = I->Tracker;
      int list_id = ExecutiveGetNamesListFromPattern(G, name, true, false);
      int iter_id = TrackerNewIter(I_Tracker, 0, list_id);
      SpecRec *rec;

      while(TrackerIterNextCandInList(I_Tracker, iter_id, (TrackerRef **) (void *) &rec)) {
        if(rec) {
          ObjectGroup *objGroup = NULL;
          if((rec->type == cExecObject) && (rec->obj->type == cObjectGroup)) {
            objGroup = (ObjectGroup *) rec->obj;
          }

          switch (action) {
          case cExecutiveGroupOpen:
            if(objGroup)
              objGroup->OpenOrClosed = 1;
            break;
          case cExecutiveGroupClose:
            if(objGroup)
              objGroup->OpenOrClosed = 0;
            break;
          case cExecutiveGroupToggle:
            if(objGroup)
              objGroup->OpenOrClosed = !objGroup->OpenOrClosed;
            break;
          case cExecutiveGroupEmpty:
            if(objGroup) {
              SpecRec *rec2 = NULL;
              while(ListIterate(I->Spec, rec2, next)) {
                if((rec2->group == rec)
                   || WordMatchExact(G, rec2->group_name, rec->name, ignore_case)) {
                  rec2->group = NULL;
                  rec2->group_name[0] = 0;
                }
              }
            }
            break;
          case cExecutiveGroupPurge:
            if(objGroup) {
              ExecutiveGroupPurge(G, rec);
            }
            break;
          case cExecutiveGroupExcise:
            if(objGroup) {

              if(rec->group_name[0]) {
                /* cascade group members up to the surrounding group */
                SpecRec *rec2 = NULL;
                while(ListIterate(I->Spec, rec2, next)) {
                  if((rec2->group == rec) ||
                     WordMatch(G, rec->name, rec2->group_name, ignore_case)) {
                    strcpy(rec2->group_name, rec->group_name);
                    rec2->group = rec->group;
                  }
                }
              } else if((rec->type == cExecObject) && (rec->obj->type == cObjectGroup)) {
                /* and/or delete their group membership */
                SpecRec *rec2 = NULL;
                while(ListIterate(I->Spec, rec2, next)) {
                  if((rec2->group == rec) ||
                     WordMatch(G, rec->name, rec2->group_name, ignore_case)) {
                    rec2->group_name[0] = 0;
                    rec2->group = nullptr;
                  }
                }
              }
              ExecutiveDelete(G, rec->name);
            }
            break;
          }
        }
      }
      TrackerDelList(I_Tracker, list_id);
      TrackerDelIter(I_Tracker, iter_id);
      ExecutiveInvalidateGroups(G, true);
    }
  } else {
    if(obj && (obj->type == cObjectGroup)) {
      ObjectGroup *objGroup = (ObjectGroup *) obj;
      switch (action) {
      case cExecutiveGroupOpen:
        objGroup->OpenOrClosed = 1;
        break;
      case cExecutiveGroupClose:
        objGroup->OpenOrClosed = 0;
        break;
      case cExecutiveGroupToggle:
        objGroup->OpenOrClosed = !objGroup->OpenOrClosed;
        break;
      }
      if(members[0] && (action != cExecutiveGroupRemove))
        action = cExecutiveGroupAdd;

      switch (action) {
      case cExecutiveGroupAdd:
        {
          for (auto& rec : ExecutiveGetSpecRecsFromPattern(G, members, true, false)) {
            if (rec.type != cExecObject || (rec.type == cExecObject && rec.obj != obj)) {
              UtilNCopy(rec.group_name, valid_name, sizeof(WordType));
              if(!quiet) {
                PRINTFB(G, FB_Executive, FB_Actions)
                  " Executive: adding '%s' to group '%s'.\n", rec.name, rec.group_name
                  ENDFB(G);
              }
            }
          }
        }
        break;
      }

      ExecutiveInvalidateGroups(G, true);
    }
  }
  return ok;
}

static int ExecutiveGetUniqueIDAtomVLADict(PyMOLGlobals * G,
                                            ExecutiveObjectOffset ** return_vla,
                                            OVOneToOne ** return_dict)
{
  CExecutive *I = G->Executive;
  OVOneToOne *o2o = OVOneToOne_New(G->Context->heap);
  ExecutiveObjectOffset *vla = VLAlloc(ExecutiveObjectOffset, 1000);
  int n_oi = 0;
  {
    SpecRec *rec = NULL;
    while(ListIterate(I->Spec, rec, next)) {
      if(rec->type == cExecObject) {
        if(rec->obj->type == cObjectMolecule) {
          ObjectMolecule *obj = (ObjectMolecule *) rec->obj;
          int a, id, n_atom = obj->NAtom;
          const AtomInfoType *ai = obj->AtomInfo.data();
          for(a = 0; a < n_atom; a++) {
            if((id = ai->unique_id)) {
              if(OVOneToOne_GetForward(o2o, id).status == OVstatus_NOT_FOUND) {
                if(OVreturn_IS_OK(OVOneToOne_Set(o2o, id, n_oi))) {
                  VLACheck(vla, ExecutiveObjectOffset, n_oi);
                  vla[n_oi].obj = obj;
                  vla[n_oi].atm = a;
                  n_oi++;
                }
              }
            }
            ai++;
          }
        }
      }
    }
  }
  *return_dict = o2o;
  VLASize(vla, ExecutiveObjectOffset, n_oi);
  *return_vla = vla;
  return 1;
}

int ExecutiveDrawCmd(PyMOLGlobals * G, int width, int height, int antialias,
                     int entire_window, int quiet)
{
  CExecutive *I = G->Executive;
  if((width <= 0) && (height <= 0)) {
    SceneGetWidthHeight(G, &width, &height);
  }
  if(antialias < 0)
    antialias = SettingGetGlobal_i(G, cSetting_antialias);
  if(entire_window) {
    SceneInvalidateCopy(G, false);
    OrthoDirty(G);
    I->CaptureFlag = true;
  } else {
    if(SettingGetGlobal_i(G, cSetting_draw_mode) == -1) {
      ExecutiveSetSettingFromString(G, cSetting_draw_mode, "-2", "", -1, true, true);
      SceneUpdate(G, false);
    }
    SceneDeferImage(G, width, height, NULL, antialias, -1.0, cMyPNG_FormatPNG, quiet);
  }
  return 1;
}

int ExecutiveMatrixCopy2(PyMOLGlobals * G,
                         pymol::CObject * source_obj, pymol::CObject * target_obj,
                         int source_mode, int target_mode,
                         int source_state, int target_state,
                         int target_undo, int log, int quiet)
{
  /*  mode 0: raw coordinates, as per the txf history
     mode 1: object TTT matrix
     mode 2: state matrix */

  int ok = true;
  int copy_ttt_too = false;
  int matrix_mode = SettingGetGlobal_i(G, cSetting_matrix_mode);
  if(matrix_mode < 0)
    matrix_mode = 0; /* for now */

  if((source_mode < 0) && (target_mode < 0)) {
    copy_ttt_too = true;
  }
  if(source_mode < 0)
    source_mode = matrix_mode;
  if(target_mode < 0)
    target_mode = matrix_mode;

  switch (source_mode) {
  case 0:                      /* txf history is the source matrix */
    {
      double *history = NULL;
      int found = ExecutiveGetObjectMatrix2(G, source_obj, source_state, &history, false);
      if(found) {
        switch (target_mode) {
        case 0:                /* apply changes to coordinates in the target object */
          {
            double temp_inverse[16];
            if(target_undo) {
              double *target_history = NULL;
              int target_found = ExecutiveGetObjectMatrix2(G, source_obj,
                                                           target_state,
                                                           &target_history,
                                                           false);
              if(target_found && target_history) {
                invert_special44d44d(target_history, temp_inverse);
                if(history) {
                  right_multiply44d44d(temp_inverse, history);
                  history = temp_inverse;
                } else {
                  history = temp_inverse;
                }
              }
              {
                float historyf[16];
                if(history) {
                  convert44d44f(history, historyf);
                } else {
                  identity44f(historyf);
                }
                ExecutiveTransformObjectSelection2(G, target_obj, target_state,
                                                   "", log, historyf, true, false);
              }
            }
            if(copy_ttt_too) {
              const float *tttf;
              int found = ObjectGetTTT(source_obj, &tttf, -1);
              if(found) {
                ObjectSetTTT(target_obj, tttf, -1, -1);
                target_obj->invalidate(cRepNone, cRepInvExtents, -1);
              }
            }
          }
          break;
        case 1:                /* applying changes to the object's TTT matrix */
          if(history) {
            float tttf[16];
            convertR44dTTTf(history, tttf);
            ObjectSetTTT(target_obj, tttf, -1, -1);
          } else {
            ObjectSetTTT(target_obj, NULL, -1, -1);
          }
          target_obj->invalidate(cRepNone, cRepInvExtents, -1);
          break;
        case 2:                /* applying changes to the state matrix */
          ok = ExecutiveSetObjectMatrix2(G, target_obj, target_state, history);
          break;
        }
        break;
      }
    }
    break;
  case 1:                      /* from the TTT matrix */
    {
      /* note that for now we're forcing states to be -1 */
      /* in the future, we may have per-state TTTs -- though right now the
         view matrices serve that purpose */

      const float *tttf;
      int found = ObjectGetTTT(source_obj, &tttf, -1);
      if(found) {
        switch (target_mode) {
        case 0:                /* coordinates & history unsupported.. */
          /* should complain */
          break;
        case 1:                /* TTT */
          ObjectSetTTT(target_obj, tttf, -1, -1);
          target_obj->invalidate(cRepNone, cRepInvExtents, -1);
          break;
        case 2:                /* State */
          if(tttf) {
            double homo[16];
            convertTTTfR44d(tttf, homo);
            ok = ExecutiveSetObjectMatrix2(G, target_obj, -1, homo);
          } else {
            ok = ExecutiveSetObjectMatrix2(G, target_obj, -1, NULL);
          }
          break;
        }
      }
    }
    break;
  case 2:                      /* from the state matrix */
    {
      double *homo;
      int found = ExecutiveGetObjectMatrix2(G, source_obj, source_state, &homo, false);
      if(found) {
        switch (target_mode) {
        case 0:                /* coordinates & history */
          /* TODO */
          break;
        case 1:                /* TTT */
          if(homo) {
            float tttf[16];
            convertR44dTTTf(homo, tttf);
            ObjectSetTTT(target_obj, tttf, -1, -1);
            target_obj->invalidate(cRepNone, cRepInvExtents, -1);
          } else {
            ObjectSetTTT(target_obj, NULL, -1, -1);
            target_obj->invalidate(cRepNone, cRepInvExtents, -1);
          }
          break;
        case 2:                /* State */
          ok = ExecutiveSetObjectMatrix2(G, target_obj, target_state, homo);
          if(copy_ttt_too) {
            const float *tttf;
            int found = ObjectGetTTT(source_obj, &tttf, -1);
            if(found) {
              ObjectSetTTT(target_obj, tttf, -1, -1);
              target_obj->invalidate(cRepNone, cRepInvExtents, -1);
            }
          }
          break;
        }
      }
    }
    break;
  }
  SceneInvalidate(G);
  return ok;
}

int ExecutiveMatrixCopy(PyMOLGlobals * G,
                        const char *source_name, const char *target_name,
                        int source_mode, int target_mode,
                        int source_state, int target_state,
                        int target_undo, int log, int quiet)
{
  /*  mode 0: raw coordinates, as per the txf history
     mode 1: object TTT matrix
     mode 2: state matrix
     mode 3 (source only): camera matrix transformation */

  CExecutive *I = G->Executive;
  CTracker *I_Tracker = I->Tracker;
  SpecRec *src_rec = NULL;
  int ok = true;
  int copy_ttt_too = false;
  int matrix_mode = SettingGetGlobal_i(G, cSetting_matrix_mode);
  if(matrix_mode < 0)
    matrix_mode = 0; /* for now */

  if((source_mode < 0) && (target_mode < 0)) {
    copy_ttt_too = true;
  }

  if(source_mode < 0)
    source_mode = matrix_mode;
  if(target_mode < 0)
    target_mode = matrix_mode;
  if(source_name[0] == 0) {
    source_mode = 3;
    target_undo = 0;
  } else
    src_rec = ExecutiveFindSpec(G, source_name);

  if (source_mode != 3 && !src_rec) {
    PRINTFB(G, FB_Executive, FB_Warnings)
      " %s-Warning: Can't find source object '%s'.\n", __FUNCTION__, source_name
      ENDFB(G);
  }

  int list_id = ExecutiveGetNamesListFromPattern(G, target_name, true, cExecExpandKeepGroups);

  if (!list_id) {
    PRINTFB(G, FB_Executive, FB_Warnings)
      " %s-Warning: No match for target '%s'.\n", __FUNCTION__, target_name
      ENDFB(G);
  }

  switch (source_mode) {
  case 0:                      /* txf history is the source matrix */
    {
      double *history = NULL;
      int found = ExecutiveGetObjectMatrix(G, source_name, source_state, &history, false);
      if(found) {

        int iter_id = TrackerNewIter(I_Tracker, 0, list_id);
        SpecRec *rec;

        while(TrackerIterNextCandInList(I_Tracker, iter_id,
                                        (TrackerRef **) (void *) &rec)) {
          if(rec && (rec != src_rec)) {
            switch (rec->type) {
            case cExecObject:

              switch (target_mode) {
              case 0:          /* apply changes to coordinates in the target object */
                {
                  double temp_inverse[16];
                  if(target_undo) {
                    double *target_history = NULL;
                    int target_found = ExecutiveGetObjectMatrix(G, rec->name,
                                                                target_state,
                                                                &target_history,
                                                                false);
                    if(target_found && target_history) {
                      invert_special44d44d(target_history, temp_inverse);
                      if(history) {
                        right_multiply44d44d(temp_inverse, history);
                        history = temp_inverse;
                      } else {
                        history = temp_inverse;
                      }
                    }
                  }
                  {
                    float historyf[16];
                    if(history) {
                      convert44d44f(history, historyf);
                    } else {
                      identity44f(historyf);
                    }
                    ExecutiveTransformObjectSelection(G, rec->name, target_state,
                                                      "", log, historyf, true, false);
                  }
                  if(copy_ttt_too) {
                    const float *tttf;
                    int found = ExecutiveGetObjectTTT(G, source_name, &tttf, -1, quiet);
                    if(found) {
                      ExecutiveSetObjectTTT(G, rec->name, tttf, -1, quiet, -1);
                    }
                  }
                }
                break;
              case 1:          /* applying changes to the object's TTT matrix */
                if(history) {
                  float tttf[16];
                  convertR44dTTTf(history, tttf);
                  ExecutiveSetObjectTTT(G, rec->name, tttf, -1, quiet, -1);
                } else {
                  ExecutiveSetObjectTTT(G, rec->name, NULL, -1, quiet, -1);
                }
                /* to do: logging, return values, etc. */
                break;
              case 2:          /* applying changes to the state matrix */
                ok = ExecutiveSetObjectMatrix(G, rec->name, target_state, history);
                break;
              }
              break;
            }
          }
        }
        TrackerDelIter(I_Tracker, iter_id);
      }
    }
    break;
  case 1:                      /* from the TTT matrix */
    {
      /* note that for now we're forcing states to be -1 */
      /* in the future, we may have per-state TTTs -- though right now the
         view matrices serve that purpose */

      const float *tttf;
      int found = ExecutiveGetObjectTTT(G, source_name, &tttf, -1, quiet);
      if(found) {

        int iter_id = TrackerNewIter(I_Tracker, 0, list_id);
        SpecRec *rec;

        while(TrackerIterNextCandInList(I_Tracker, iter_id,
                                        (TrackerRef **) (void *) &rec)) {
          if(rec && (rec != src_rec)) {

            switch (rec->type) {
            case cExecObject:

              switch (target_mode) {
              case 0:          /* coordinates & history unsupported.. */
                /* should complain */
                break;
              case 1:          /* TTT */
                ExecutiveSetObjectTTT(G, rec->name, tttf, -1, quiet, -1);
                break;
              case 2:          /* State */
                if(tttf) {
                  double homo[16];
                  convertTTTfR44d(tttf, homo);
                  ok = ExecutiveSetObjectMatrix(G, rec->name, -1, homo);
                } else {
                  ok = ExecutiveSetObjectMatrix(G, rec->name, -1, NULL);
                }
                break;
              }
            }
          }
        }

        TrackerDelIter(I_Tracker, iter_id);
      }
    }
    break;
  case 2:                      /* from the state matrix */
    {
      double *homo;
      int found = ExecutiveGetObjectMatrix(G, source_name, source_state, &homo, false);
      if(found) {

        int iter_id = TrackerNewIter(I_Tracker, 0, list_id);
        SpecRec *rec;

        while(TrackerIterNextCandInList(I_Tracker, iter_id,
                                        (TrackerRef **) (void *) &rec)) {
          if(rec && (rec != src_rec)) {
            switch (rec->type) {
            case cExecObject:

              switch (target_mode) {
              case 0:          /* coordinates & history */
                /* TODO */
                break;
              case 1:          /* TTT */
                if(homo) {
                  float tttf[16];
                  convertR44dTTTf(homo, tttf);
                  ExecutiveSetObjectTTT(G, rec->name, tttf, -1, quiet, -1);
                } else {
                  ExecutiveSetObjectTTT(G, rec->name, NULL, -1, quiet, -1);
                }
                break;
              case 2:          /* State */
                ok = ExecutiveSetObjectMatrix(G, rec->name, target_state, homo);
                if(copy_ttt_too) {
                  const float *tttf;
                  int found = ExecutiveGetObjectTTT(G, source_name, &tttf, -1, quiet);
                  if(found) {
                    ExecutiveSetObjectTTT(G, rec->name, tttf, -1, quiet, -1);
                  }
                }
                break;
              }
            }
          }
        }
        TrackerDelIter(I_Tracker, iter_id);
      }
    }
    break;
  case 3:                      /* camera */
    {
      SceneViewType view;
      double homo[16], *history;
      int iter_id = TrackerNewIter(I_Tracker, 0, list_id);
      SpecRec *rec;
      SceneGetView(G, view);
      homo[0] = view[0];
      homo[1] = view[4];
      homo[2] = view[8];
      homo[3] = -(view[0] * view[19] + view[4] * view[20] + view[8] * view[21]);
      homo[4] = view[1];
      homo[5] = view[5];
      homo[6] = view[9];
      homo[7] = -(view[1] * view[19] + view[5] * view[20] + view[9] * view[21]);
      homo[8] = view[2];
      homo[9] = view[6];
      homo[10] = view[10];
      homo[11] = -(view[2] * view[19] + view[6] * view[20] + view[10] * view[21]);
      homo[12] = 0.0;
      homo[13] = 0.0;
      homo[14] = 0.0;
      homo[15] = 1.0;
      history = homo;
      while(TrackerIterNextCandInList(I_Tracker, iter_id, (TrackerRef **) (void *) &rec)) {
        if(rec && (rec != src_rec)) {
          switch (rec->type) {
          case cExecObject:

            switch (target_mode) {
            case 0:            /* apply changes to coordinates in the target object */
              {
                double temp_inverse[16];
                if(target_undo) {
                  double *target_history = NULL;
                  int target_found = ExecutiveGetObjectMatrix(G, rec->name,
                                                              target_state,
                                                              &target_history,
                                                              false);
                  if(target_found && target_history) {
                    invert_special44d44d(target_history, temp_inverse);
                    if(history) {
                      right_multiply44d44d(temp_inverse, history);
                      history = temp_inverse;
                    } else {
                      history = temp_inverse;
                    }
                  }
                }
                {
                  float historyf[16];
                  if(history) {
                    convert44d44f(history, historyf);
                  } else {
                    identity44f(historyf);
                  }
                  ExecutiveTransformObjectSelection(G, rec->name, target_state,
                                                    "", log, historyf, true, false);
                }
              }
              break;
            case 1:            /* applying changes to the object's TTT matrix */
              if(history) {
                float tttf[16];
                convertR44dTTTf(history, tttf);
                ExecutiveSetObjectTTT(G, rec->name, tttf, -1, quiet, -1);
              } else {
                ExecutiveSetObjectTTT(G, rec->name, NULL, -1, quiet, -1);
              }
              /* to do: logging, return values, etc. */
              break;
            case 2:            /* applying changes to the state matrix */
              ok = ExecutiveSetObjectMatrix(G, rec->name, target_state, history);
              break;
            }
            break;
          }
        }
        TrackerDelIter(I_Tracker, iter_id);
      }
    }
    break;
  }

  TrackerDelList(I_Tracker, list_id);

  SceneInvalidate(G);
  return ok;
}

void ExecutiveInvalidateMapDependents(
    PyMOLGlobals* G, const char* map_name, const char* new_name)
{
  CExecutive *I = G->Executive;
  SpecRec *rec = NULL;
  while(ListIterate(I->Spec, rec, next)) {
    if(rec->type == cExecObject) {
      switch (rec->obj->type) {
      case cObjectMesh:
        ObjectMeshInvalidateMapName((ObjectMesh *) rec->obj, map_name, new_name);
        break;
      case cObjectSurface:
        ObjectSurfaceInvalidateMapName((ObjectSurface *) rec->obj, map_name, new_name);
        break;
      case cObjectVolume:
        ObjectVolumeInvalidateMapName((ObjectVolume *) rec->obj, map_name, new_name);
        break;
      }
    }
  }
  SceneInvalidate(G);
}

pymol::Result<> ExecutiveResetMatrix(
    PyMOLGlobals* G, const char* name, int mode, int state, int log, int quiet)
{
  CExecutive *I = G->Executive;
  CTracker *I_Tracker = I->Tracker;
  int list_id = ExecutiveGetNamesListFromPattern(G, name, true, true);
  int iter_id = TrackerNewIter(I_Tracker, 0, list_id);
  SpecRec *rec;
  int matrix_mode = SettingGetGlobal_i(G, cSetting_matrix_mode);
  if(matrix_mode < 0)
    matrix_mode = 0; /* for now */

  if(mode < 0)
    mode = matrix_mode;

  bool obj_found = false;
  while(TrackerIterNextCandInList(I_Tracker, iter_id, (TrackerRef **) (void *) &rec)) {
    if(rec && (rec->type == cExecObject)) {
      /*  pymol::CObject *obj = ExecutiveFindObjectByName(G,name); */
      pymol::CObject *obj = rec->obj;
      if(obj) {
        obj_found = true;
        switch (obj->type) {
        case cObjectMolecule:
          switch (mode) {
          case 0:              /* transformations already applied to the coordinates */
            for (StateIterator iter(rec->obj, state); iter.next();) {
              double *history = NULL;
              bool found = ExecutiveGetObjectMatrix2(G, rec->obj, iter.state, &history, false);
              if(found && history) {
                double temp_inverse[16];
                float historyf[16];
                invert_special44d44d(history, temp_inverse);
                convert44d44f(temp_inverse, historyf);
                ExecutiveTransformObjectSelection2(G, rec->obj, iter.state, "",
                                                  log, historyf, true, false);
              }
            }
            break;
          case 1:              /* operate on the TTT display matrix */
            ObjectResetTTT(obj,SettingGetGlobal_b(G,cSetting_movie_auto_store));
            obj->invalidate(cRepNone, cRepInvExtents, -1);

            break;
          case 2:              /* applying changes to the state matrix */
            {
              double ident[16];
              identity44d(ident);
              ExecutiveSetObjectMatrix(G, rec->name, state, ident);
            }
            break;
          }
          break;
        default:
          if (auto* objstate = obj->getObjectState(state)) {
            ObjectStateResetMatrix(objstate);
            obj->invalidate(cRepNone, cRepInvExtents, state);
          }
          break;
        }
      }
    }
  }
  if (!obj_found) {
    return pymol::make_error("No object found");
  }
  return {};
}

static double ret_mat[16];      /* UGH ..not thread-safe */

static int ExecutiveGetObjectMatrix2(PyMOLGlobals * G, pymol::CObject * obj, int state,
                                     double **matrix, int incl_ttt)
{
  /* right now, this only makes sense for molecule objects -- but in
     time all objects should have per-state matrices */

  int ok = false;
  if(state < 0) {
    /* to do -- TTT only */
  } else {
    if (auto* objstate = obj->getObjectState(state)) {
      *matrix = ObjectStateGetMatrix(objstate);
      ok = true;
    }

    if(ok && incl_ttt) {
      const float *ttt;
      double tttd[16];
      if(ObjectGetTTT(obj, &ttt, -1)) {
        convertTTTfR44d(ttt, tttd);
        if(*matrix) {
          copy44d(*matrix, ret_mat);
        } else {
          identity44d(ret_mat);
        }
        left_multiply44d44d(tttd, ret_mat);
        *matrix = ret_mat;
      }
    }
  }
  return ok;
}

int ExecutiveGetObjectMatrix(PyMOLGlobals * G, const char *name, int state, double **matrix,
                             int incl_ttt)
{
  int ok = false;
  pymol::CObject *obj = ExecutiveFindObjectByName(G, name);
  if(obj) {
    return ExecutiveGetObjectMatrix2(G, obj, state, matrix, incl_ttt);
  }
  return ok;
}

static int ExecutiveSetObjectMatrix2(PyMOLGlobals * G, pymol::CObject * obj, int state,
                                     double *matrix)
{
  /* -1 for the TTT matrix, 0 or greater for the state matrix */

  /* right now, this only makes sense for molecule objects -- but in
     time all objects should have per-state matrices */
  int ok = false;
  if(state < 0) {

  } else {
    if (auto* objstate = obj->getObjectState(state)) {
      ObjectStateSetMatrix(objstate, matrix);
      ok = true;
    }
  }
  return ok;
}

int ExecutiveSetObjectMatrix(PyMOLGlobals * G, const char *name, int state, double *matrix)
{
  int ok = false;
  pymol::CObject *obj = ExecutiveFindObjectByName(G, name);
  if(obj) {
    return ExecutiveSetObjectMatrix2(G, obj, state, matrix);
  }
  return ok;
}

static int ExecutiveCountNames(PyMOLGlobals * G)
{
  int count = 0;
  CExecutive *I = G->Executive;
  SpecRec *rec = NULL;
  while(ListIterate(I->Spec, rec, next))
    count++;

  return (count);
}

static int ReorderOrderFn(PyMOLGlobals* G, const SpecRec* const* rec, int l, int r)
{
  return (WordCompare(G, rec[l]->name, rec[r]->name, true) <= 0);
}

static int SpecRecListPopulate(SpecRec ** list, SpecRec * first, const char * group_name)
{
  /* Add items to list such that group items are ordered behind their group records.
   * author: Thomas Holder, 2013-08
   * runtime: O(N * G) where N is the number of SpecRecs and G the number of groups
   * args: list: pointer to the end of the list
   *       first: pointer to G->Executive->Spec (root of linked list)
   *       group_name: if empty string, append ungrouped items, otherwise append group items
   * returns: number of appended items
   */
  SpecRec *rec;
  int a = 0;
  for(rec = first; rec; rec = rec->next) {
    if(!strcmp(group_name, rec->group_name)) {
      list[a++] = rec;
      if(rec->type == cExecObject && rec->obj->type == cObjectGroup)
        a += SpecRecListPopulate(list + a, first, rec->name);
    }
  }
  return a;
}

pymol::Result<> ExecutiveOrder(PyMOLGlobals * G, pymol::zstring_view s1View, int sort, int location)
{
  CExecutive *I = G->Executive;
  CTracker *I_Tracker = I->Tracker;
  auto s1 = s1View.c_str();
  CWordList *word_list = WordListNew(G, s1);
  int n_names = ExecutiveCountNames(G);

  if(n_names) {
    SpecRec **list, **subset, **sorted;
    int *index = NULL;
    int n_sel;
    int source_row = -1;
    int min_row = -1;
    list = pymol::malloc<SpecRec *>(n_names);
    subset = pymol::calloc<SpecRec *>(n_names);
    sorted = pymol::calloc<SpecRec *>(n_names);
    index = pymol::malloc<int>(n_names);
    if(list && subset) {
      /* create an array of current names */
      {
        int a = 0;
        /* copy all names into array */
        /* update Thomas Holder 2013-08: order group members behind group.
         * fixes PYMOL-1382 (eventually).
         * I don't know if this is the best place for the fix, but it's my best guess.
         */
#ifndef NDEBUG
        auto const list_size =
#endif
            SpecRecListPopulate(list, I->Spec, "");
        assert(list_size == n_names);
        /* unlink them */
        for(a = 0; a < n_names; a++) {
          list[a]->next = NULL;
        }
      }
      /* transfer matching names to the subset array */
      {
        int a;
        int entry;
        int min_entry = word_list->n_word;
        const char *word = NULL;
        int word_iter = 0;
        while(WordListIterate(G, word_list, &word, &word_iter)) {
          int list_id = ExecutiveGetNamesListFromPattern(G, word, true, false);
          SpecRec *rec = NULL;
          entry = word_iter - 1;
          for(a = n_names - 1; a > 0; a--) {    /* skipping zeroth */
            int iter_id = TrackerNewIter(I_Tracker, 0, list_id);
            while(TrackerIterNextCandInList
                  (I_Tracker, iter_id, (TrackerRef **) (void *) &rec)) {
              if(rec == list[a]) {
                if((a < min_row) || (min_row < 0))
                  min_row = a;
                if(entry <= min_entry) {
                  source_row = a;       /* where will new list be inserted... */
                  min_entry = entry;
                }
                /* ensure that each record appears only once */
                rec->next = subset[entry];
                subset[entry] = rec;
                list[a] = NULL;
              }
            }
            TrackerDelIter(I_Tracker, iter_id);
          }
          TrackerDelList(I_Tracker, list_id);
        }
        if(word_list->n_word && WordMatchExact(G, word_list->start[0], cKeywordAll, true))
          location = -1;        /* set to top if "all" is first in list */
      }
      /* expand the selected entries */
      {
        SpecRec *rec, *last;
        int b;
        n_sel = 0;
        for(b = 0; b < word_list->n_word; b++) {
          rec = subset[b];
          while(rec) {
            sorted[n_sel++] = rec;
            last = rec;
            rec = rec->next;
            last->next = NULL;
          }
        }
      }
      /* sort the selected entries, if requested */
      if(sort) {
        UtilCopyMem(subset, sorted, sizeof(SpecRec *) * n_sel);
        {
          int a;
          UtilSortIndexGlobals(G, n_sel, subset, index,
                               (UtilOrderFnGlobals *) ReorderOrderFn);
          for(a = 0; a < n_sel; a++) {
            sorted[a] = subset[index[a]];
          }
        }
      }
      /* reassemble the list using the new order */
      {
        SpecRec *spec = NULL;
        SpecRec *last = NULL;
        int a, b;
        int flag;
        for(a = 0; a < n_names; a++) {
          flag = false;
          if(sorted) {          /* not yet added */
            switch (location) {
            case -1:           /* top */
              if(a == 1)
                flag = true;
              break;
            case -2:           /* upper */
              if(min_row >= 0) {
                if(a == min_row)
                  flag = true;
              } else if(!list[a])
                flag = true;
              break;
            case 0:            /* current */
              if(source_row >= 0) {
                if(a == source_row)
                  flag = true;
              } else if(!list[a])
                flag = true;
              break;
            }
          }
          if(flag) {
            for(b = 0; b < n_sel; b++) {
              if(sorted[b]) {
                if(last)
                  last->next = sorted[b];
                last = sorted[b];
                if(!spec)
                  spec = last;
              }
            }
            FreeP(sorted);
          }
          if(list[a]) {
            if(last)
              last->next = list[a];
            last = list[a];
            if(!spec)
              spec = last;
          }
        }
        if(sorted) {            /* still not yet readded? */
          for(b = 0; b < n_sel; b++) {
            if(sorted[b]) {
              if(last)
                last->next = sorted[b];
              last = sorted[b];
              if(!spec)
                spec = last;
            }
          }
        }
        I->Spec = spec;
        OrthoDirty(G);
        SeqChanged(G);
      }
      FreeP(index);
      FreeP(sorted);
      FreeP(list);
      FreeP(subset);
    }
    ExecutiveInvalidatePanelList(G);
  }
  WordListFreeP(word_list);
  return {};
}

pymol::vla<ObjectMolecule*> ExecutiveGetObjectMoleculeVLA(PyMOLGlobals * G, const char *sele)
{
  ObjectMolecule **result = NULL;
  int s1 = SelectorIndexByName(G, sele);
  if(s1 >= 0) {
    ObjectMoleculeOpRec op;
    ObjectMoleculeOpRecInit(&op);
    op.code = OMOP_GetObjects;
    op.obj1VLA = (ObjectMolecule **) VLAlloc(ObjectMolecule *, 10);
    op.i1 = 0;
    ExecutiveObjMolSeleOp(G, s1, &op);
    result = (ObjectMolecule **) op.obj1VLA;
    VLASize(result, ObjectMolecule *, op.i1);
  }
  return pymol::vla_take_ownership(result);
}


/* #define ExecLineHeight 18 */
#define ExecClickMargin DIP2PIXEL(2)
#define ExecTopMargin 0
#define ExecToggleMargin DIP2PIXEL(2)
#define ExecLeftMargin DIP2PIXEL(1)
#define ExecRightMargin 0
#define ExecToggleWidth DIP2PIXEL(17)
#define ExecToggleSize DIP2PIXEL(16)
#define ExecToggleTextShift DIP2PIXEL(4)

int ExecutiveSetDrag(PyMOLGlobals * G, const char *name, int quiet,int mode)
{
  char drag_name[] = cEditorDrag;
  int set_flag = false;
  int need_sele = true;
  int result = true;
  if(name[0]) {
    pymol::CObject *obj = ExecutiveFindObjectByName(G, name);
    if(obj) {
      EditorSetDrag(G, obj, -1, quiet, SceneGetState(G));
      set_flag = true;
    } else {
      SpecRec *rec = ExecutiveFindSpec(G, name);
      if(rec) {
        if(rec->type == cExecSelection) {
          SelectorCreate(G, drag_name, name, NULL, true, NULL);
          need_sele = false;
          {
            int sele = SelectorIndexByName(G, drag_name);
            ObjectMolecule *objMol = SelectorGetSingleObjectMolecule(G, sele);
            if(objMol) {
              if(mode>0) 
                sele = -1; /* force drag by matrix */
              EditorSetDrag(G, objMol, sele, quiet, SceneGetState(G));
              set_flag = true;
            } else {
              PRINTFB(G, FB_Executive, FB_Errors)
                " Drag-Error: selection spans more than one object.\n" ENDFB(G);
            }
          }
        } else if(rec->type == cExecObject) {
          switch (rec->obj->type) {
          case cObjectGroup:
            PRINTFB(G, FB_Executive, FB_Errors)
              " Drag-Error: cannot drag group objects yet.\n" ENDFB(G);
            break;
            
          }
          result = false;
        }
      }
    }
    result = set_flag;
    if(!result) {
      EditorInactivate(G);
      PRINTFB(G, FB_Executive, FB_Errors)
        " Drag-Error: invalid or empty selection." ENDFB(G);
    } else if(EditorDraggingObjectMatrix(G)) {
      SelectorCreate(G, drag_name, "none", NULL, true, NULL);    
    } else if(need_sele && (obj->type == cObjectMolecule) && (!EditorDraggingObjectMatrix(G))) {
      SelectorCreate(G, drag_name, obj->Name, (ObjectMolecule*)(void*)obj, true, NULL);     /* for indication only */
    }
  } else {
    EditorInactivate(G);
  }
  return result;
}

int ExecutivePop(PyMOLGlobals * G, const char *target, const char *source, int quiet)
{
  int ok = true;
  int src;
  int result = 0;

  ExecutiveDelete(G, target);
  if(ExecutiveFindObjectMoleculeByName(G, source)) {
    ok = false;
    PRINTFB(G, FB_Executive, FB_Errors)
      " Pop-Error: source selection '%s' can't be an object.\n", source ENDFB(G);

  } else {
    src = SelectorIndexByName(G, source);
    if(src < 0)
      ok = false;
    if(!ok) {
      PRINTFB(G, FB_Executive, FB_Errors)
        " Pop-Error: invalid source selection name '%s'\n", source ENDFB(G);
    } else {
      ObjectMoleculeOpRec op;

      ObjectMoleculeOpRecInit(&op);
      op.code = OMOP_Pop;
      SelectorCreateEmpty(G, target, true);
      op.i1 = SelectorIndexByName(G, target);
      op.i2 = 1;
      op.i3 = 0;
      ExecutiveObjMolSeleOp(G, src, &op);
      result = op.i3;
    }
  }
  if(!result)
    ExecutiveDelete(G, target);
  if(!ok)
    return -1;
  else
    return result;
}

/**
 * Return the selector index of the "active" alignment.
 */
int ExecutiveGetActiveAlignmentSele(PyMOLGlobals * G)
{
  const char* alignment = ExecutiveGetActiveAlignment(G);
  if (alignment && alignment[0]) {
    return SelectorIndexByName(G, alignment);
  }
  return -1;
}

/**
 * Return the name of the "active" alignment. That is:
 *
 * 1) The "seq_view_alignment" setting if it's set
 * 2) or the name of the first enabled alignment object
 * 3) or NULL
 */
const char* ExecutiveGetActiveAlignment(PyMOLGlobals* G)
{
  const char *alignment = SettingGetGlobal_s(G, cSetting_seq_view_alignment);
  if(alignment && alignment[0]) {       /* explicit alignment setting name */
    return alignment;
  } else {                      /* otherwise, use the first active alignment */
    SpecRec *rec = NULL;
    CExecutive *I = G->Executive;
    while(ListIterate(I->Spec, rec, next)) {
      if(rec->visible) {
        if(rec->type == cExecObject)
          if(rec->obj->type == cObjectAlignment) {
            return rec->obj->Name;
          }
      }
    }
  }
  return nullptr;
}

int ExecutiveGetActiveSele(PyMOLGlobals * G)
{
  ObjectNameType name;
  if(ExecutiveGetActiveSeleName(G, name, false, false))
    return SelectorIndexByName(G, name);
  else
    return -1;

}

int ExecutiveGetActiveSeleName(PyMOLGlobals * G, char *name, int create_new, int log)
{
  /* TODO: cache/optimize to avoid table scan */

  int result = false;
  SpecRec *rec = NULL;
  CExecutive *I = G->Executive;
  while(ListIterate(I->Spec, rec, next)){
    if(rec->type == cExecSelection){
      if(rec->visible) {
        strcpy(name, rec->name);
        result = true;
      }
    }
  }
  if((!result) && create_new) {
    if(SettingGetGlobal_b(G, cSetting_auto_number_selections)) {
      int sel_num = SettingGetGlobal_i(G, cSetting_sel_counter) + 1;

      SettingSetGlobal_i(G, cSetting_sel_counter, sel_num);
      sprintf(name, "sel%02d", sel_num);
      SelectorCreateEmpty(G, name, -1);
      if(log) {
        if(SettingGetGlobal_i(G, cSetting_logging)) {
          OrthoLineType buf2;
          sprintf(buf2, "cmd.select('%s','none')\n", name);
          PLog(G, buf2, cPLog_no_flush);
        }
      }
    } else {
      sprintf(name, "sele");
      SelectorCreateEmpty(G, name, -1);
      if(log) {
        OrthoLineType buf2;
        sprintf(buf2, "cmd.select('%s','none')\n", name);
        PLog(G, buf2, cPLog_no_flush);
      }
    }
  }
  return result;
}

int ExecutiveGetActiveSeleName(PyMOLGlobals* G, std::string& name, int create_new, int log)
{
  ObjectNameType name_arg;
  auto result = ExecutiveGetActiveSeleName(G, name_arg, create_new, log);
  name = name_arg;
  return result;
}

pymol::Result<> ExecutiveFixChemistry(PyMOLGlobals * G, const char *s1, const char *s2, int invalidate, int quiet)
{
  SETUP_SELE_DEFAULT_PREFIXED(1, cSelectionInvalid);
  SETUP_SELE_DEFAULT_PREFIXED(2, sele1);

  SpecRec *rec = NULL;
  CExecutive *I = G->Executive;

  assert((sele1 >= 0) && (sele2 >= 0));
  {
    while(ListIterate(I->Spec, rec, next)) {
      if(rec->type == cExecObject)
        if(rec->obj->type == cObjectMolecule) {
          ObjectMoleculeFixChemistry((ObjectMolecule *) rec->obj, sele1, sele2,
                                     invalidate);
        }
    }
  }
  return {};
}

pymol::Result<> ExecutiveSetObjectColor(PyMOLGlobals * G, const char *name, const char *color, int quiet)
{
  int col_ind = ColorGetIndex(G, color);
  auto obj = ExecutiveFindObjectByName(G, name);
  if(obj) {
    obj->Color = col_ind;
  } else {
    return pymol::make_error("Object ", name, " not found.");
  }
  return {};
}

int ExecutiveGetObjectColorIndex(PyMOLGlobals * G, const char *name)
{
  int result = -1;
  pymol::CObject *obj = NULL;
  obj = ExecutiveFindObjectByName(G, name);
  if(obj) {
    result = obj->Color;
  }
  return (result);
}

pymol::Result<std::array<float, 3>> ExecutiveGetAtomVertex(PyMOLGlobals * G, const char *s1, int state, int index)
{
  auto tmpsele1 = SelectorTmp::make(G, s1);
  p_return_if_error(tmpsele1);

  switch (tmpsele1->getAtomCount()) {
  case 0:
    return pymol::Error("Empty Selection");
  case 1:
    return SelectorGetSingleAtomVertex(G, tmpsele1->getIndex(), state);
  }

  assert(tmpsele1->getAtomCount() > 0);
  return pymol::Error("More than one atom found");
}

void ExecutiveMakeUnusedName(PyMOLGlobals * G, char * prefix, int length,
                             bool alwaysnumber, int start,
                             const char * pattern) {

  if (!prefix[0])
    strcpy(prefix, "obj");

  int prefixlen = strlen(prefix);
  int suffixlen = length - prefixlen;
  char * end = prefix + prefixlen;

  for (int cnt = start; alwaysnumber || ExecutiveValidName(G, prefix); ++cnt) {
    snprintf(end, suffixlen, pattern, cnt);
    alwaysnumber = false;
  }
}

std::string ExecutiveGetUnusedName(PyMOLGlobals * G, const char * prefix,
                                   bool alwaysnumber) {

  OrthoLineType unused_name;
  strcpy(unused_name, prefix);

  ObjectMakeValidName(G, unused_name);

  ExecutiveMakeUnusedName(G, unused_name, OrthoLineLength, alwaysnumber);

  return std::string(unused_name);
}

int ExecutiveProcessObjectName(PyMOLGlobals * G, const char *proposed, char *actual)
{
  int result = true;
  UtilNCopy(actual, proposed, sizeof(ObjectNameType));
  if(SettingGetGlobal_b(G, cSetting_validate_object_names))
    ObjectMakeValidName(G, actual);
  if(SettingGetGlobal_b(G, cSetting_auto_rename_duplicate_objects) || !proposed[0]) {
    ExecutiveMakeUnusedName(G, actual, sizeof(ObjectNameType), false, 2, "_%d");
  }
  return result;
}

pymol::Result<std::vector<DiscardedRec>> ExecutiveSetName(
        PyMOLGlobals * G, pymol::zstring_view old_name_view, pymol::zstring_view new_name_view, bool save)
{
  std::vector<DiscardedRec> discarded;
  SpecRec *rec = NULL;
  CExecutive *I = G->Executive;
  int found = false;
  auto ignore_case = SettingGet<bool>(G, cSetting_ignore_case);
  auto old_name = old_name_view.c_str();
  auto new_name = new_name_view.c_str();

  ObjectNameType name;
  UtilNCopy(name, new_name, sizeof(ObjectNameType));

  if (ObjectMakeValidName(name)) {
    PRINTFB(G, FB_Executive, FB_Warnings)
      " Warning: Invalid characters in '%s' have been replaced or stripped\n",
      name ENDFB(G);
  }

  if(!name[0]) {
    return pymol::make_error("Blank names not allowed.");
  } else if(WordMatchExact(G, name, cKeywordSame, ignore_case) || SelectorNameIsKeyword(G, name)) {
    return pymol::make_error("Name ", name, " is a selection keyword.");
  }
  if(!WordMatchExact(G, name, old_name, ignore_case)) {
    while(ListIterate(I->Spec, rec, next)) {
      if(found)
        break;
      switch (rec->type) {
      case cExecObject:
        if(WordMatchExact(G, rec->obj->Name, old_name, ignore_case)) {
          ExecutiveDelKey(I, rec);
          auto r = ExecutiveDelete(G, name, save);
          discarded = std::move(*r);
          assert(r);
          ObjectSetName(rec->obj, name);
          UtilNCopy(rec->name, rec->obj->Name, WordLength);
          ExecutiveAddKey(I, rec);
          if(rec->obj->type == cObjectMolecule) {
            /*
               SelectorDelete(G,old_name);
               ExecutiveUpdateObjectSelection(G,rec->obj);
             */
            SelectorSetName(G, name, old_name);
            SceneChanged(G);
            SeqChanged(G);
          }
          if (rec->obj->type == cObjectMap)
            ExecutiveInvalidateMapDependents(G, old_name, name);

          found = true;
        }
        break;
      case cExecSelection:
        if(WordMatchExact(G, rec->name, old_name, ignore_case)) {
          if(SelectorSetName(G, name, old_name)) {
            auto r = ExecutiveDelete(G, name, save); /* just in case */
            discarded = std::move(*r);
            assert(r);
            ExecutiveDelKey(I, rec);
            UtilNCopy(rec->name, name, WordLength);
            ExecutiveAddKey(I, rec);
            found = true;
            OrthoDirty(G);
          }
        }
        break;
      }
    }
    if(!found)
      return pymol::make_error("Could not find object named ", name);
    else {
      rec = NULL;
      int old_name_len = strlen(old_name);
      int new_name_len = strlen(name);
      ObjectNameType childname;
      UtilNCopy(childname, name, sizeof(ObjectNameType));
      while(ListIterate(I->Spec, rec, next)) {
        if(WordMatchExact(G, rec->group_name, old_name, ignore_case)) {
          UtilNCopy(rec->group_name, name, WordLength);
          // rename group members for group_auto_mode
          if (strncmp(rec->name, old_name, old_name_len) == 0 && rec->name[old_name_len] == '.') {
            UtilNCopy(childname + new_name_len, rec->name + old_name_len, sizeof(ObjectNameType) - new_name_len);
            ExecutiveSetName(G, rec->name, childname);
          }
        }
      }
      ExecutiveInvalidateGroups(G, false);
    }
  }
  return discarded;
}


/**
 * Load any file type which is implemented in C.
 *
 * fname:       File name, can be empty if content is provided
 * content:     File contents, if not NULL
 * content_length:      Length of "content", if not NULL
 * content_format:      File type code
 * object_name:         New object name
 * state:       Object state to start loading new coordsets in
 * zoom:        Zoom on new loaded atoms
 * discrete:    Make discrete states
 * finish:      update object selection & zoom
 * multiplex:   Split new states into objects
 * quiet:       Suppress feedback
 * object_props:        names of object properties to load
 * atom_props:          names of atom properties to load
 */
pymol::Result<>
ExecutiveLoad(PyMOLGlobals * G,
              const char *fname,
              const char *content, int content_length,
              cLoadType_t content_format,
              const char *object_name_proposed,
              int state, int zoom,
              int discrete, int finish, int multiplex, int quiet,
              const char * plugin_arg,
              const char * object_props,
              const char * atom_props,
              bool mimic)
{
  auto res = ExecutiveLoadPrepareArgs(G, fname, content, content_length,
      content_format, object_name_proposed, state, zoom, discrete, finish,
      multiplex, quiet, plugin_arg, object_props, atom_props, mimic);
  p_return_if_error(res);
  return ExecutiveLoad(G, res.result());
}

/**
 * Prepare arguments for ExecutiveLoad
 */
pymol::Result<ExecutiveLoadArgs>
ExecutiveLoadPrepareArgs(PyMOLGlobals * G,
                  pymol::null_safe_zstring_view fname,
                  const char *content, int content_length,
                  cLoadType_t content_format,
                  const char *object_name_proposed,
                  int state, int zoom,
                  int discrete, int finish, int multiplex, int quiet,
                  const char * plugin_arg,
                  const char * object_props,
                  const char * atom_props,
                  bool mimic)
{
  ExecutiveLoadArgs args;
  bool fname_null_ok = false;

  // validate proposed object name
  ObjectNameType object_name = "";
  ExecutiveProcessObjectName(G, object_name_proposed, object_name);
  args.object_name = object_name;

  // Ensure correct float parsing with scanf. It's possible to change this from
  // Python, so don't rely on a persistent global value.
  std::setlocale(LC_NUMERIC, "C");

  if (!object_props) object_props = SettingGetGlobal_s(G, cSetting_load_object_props_default);
  if (!atom_props)   atom_props   = SettingGetGlobal_s(G, cSetting_load_atom_props_default);

  // multiplex -2 -> "multiplex" setting
  // multiplex -1 -> file type dependant default
  // multiplex 1 -> split entries (don't try to load into existing object)
  if(multiplex == -2) {
    multiplex = SettingGetGlobal_i(G, cSetting_multiplex);
  }

  switch (content_format) {
  // string loading functions
  case cLoadTypePDBStr:
  case cLoadTypeVDBStr:
  case cLoadTypeCIFStr:
  case cLoadTypeMMTFStr:
  case cLoadTypeMAEStr:
  case cLoadTypeXPLORStr:
  case cLoadTypeCCP4Str:
  case cLoadTypeCCP4UnspecifiedStr:
  case cLoadTypeMRCStr:
  case cLoadTypePHIStr:
  case cLoadTypeMMDStr:
  case cLoadTypeMOLStr:
  case cLoadTypeMOL2Str:
  case cLoadTypeSDF2Str:
  case cLoadTypeXYZStr:
  case cLoadTypeDXStr:
    if (!content) {
      return pymol::Error("content is NULL");
    }
    fname_null_ok = true;
    break;
  case cLoadTypePQR:
  case cLoadTypePDBQT:
  case cLoadTypePDB:
  case cLoadTypeCIF:
  case cLoadTypeMMTF:
  case cLoadTypeMAE:
  case cLoadTypeXPLORMap:
  case cLoadTypeCCP4Map:
  case cLoadTypeCCP4Unspecified:
  case cLoadTypeMRC:
  case cLoadTypePHIMap:
  case cLoadTypeMMD:
  case cLoadTypeMOL:
  case cLoadTypeMOL2:
  case cLoadTypeSDF2:
  case cLoadTypeXYZ:
  case cLoadTypeDXMap:
    if (content) {
      fname_null_ok = true;
      break;
    }

    if (fname.empty()) {
      break;
    }

    try {
      args.content = pymol::file_get_contents(fname);
      PRINTFB(G, FB_Executive, FB_Blather)
        " %s: Loading from %s.\n", __func__, fname.c_str() ENDFB(G);
    } catch (...) {
      return pymol::Error(
          pymol::string_format("Unable to open file '%s'", fname.c_str()));
    }

    break;

  // molfile_plugin based formats
  case cLoadTypeCUBEMap:
    args.plugin = "cube";
    break;
  case cLoadTypeSpider:
    args.plugin = "spider";
    break;
  case cLoadTypeXTC:
    args.plugin = "xtc";
    break;
  case cLoadTypeTRR:
    args.plugin = "trr";
    break;
  case cLoadTypeGRO:
    args.plugin = "gro";
    break;
  case cLoadTypeG96:
    args.plugin = "g96";
    break;
  case cLoadTypeTRJ2:
    args.plugin = "trj";
    break;
  case cLoadTypeDCD:
    args.plugin = "dcd";
    break;
  case cLoadTypeDTR:
    args.plugin = "dtr";
    break;
  case cLoadTypeCMS:
    args.plugin = "mae";
    break;
  default:
    if (plugin_arg) {
      args.plugin = plugin_arg;
      // Hackish way to request a "sub-plugin". For example the "cube" plugin
      // can load maps or molecules:
      // load example.cube, molobj, format=plugin:cube:1
      // load example.cube, mapobj, format=plugin:cube:4
      auto pos = args.plugin.find(':');
      if (pos != std::string::npos) {
        args.plugin_mask = atoi(args.plugin.c_str() + pos + 1);
        args.plugin.resize(pos);
      }
    }
    break;
  }

  if (fname.empty() && !fname_null_ok) {
    return pymol::Error("This format requires a filename to load");
  }

  if (!args.plugin.empty()) {
    content_format = cLoadTypePlugin;
  }

  if (content) {
    assert(args.content.empty());
    args.content = std::string(content, content_length);
  }

  args.content_format = content_format;
  args.fname = fname.c_str();
  args.state = state;
  args.zoom = zoom;
  args.discrete = discrete;
  args.finish = finish;
  args.multiplex = multiplex;
  args.quiet = quiet;
  args.object_props = object_props;
  args.atom_props = atom_props;
  args.mimic = mimic;
  return args;
}

pymol::Result<> ExecutiveLoad(PyMOLGlobals* G, ExecutiveLoadArgs const& args)
{
  pymol::CObject* origObj = nullptr;
  const char* fname = args.fname.c_str();
  const char* content = args.content.data();
  int size = args.content.size();
  const char* object_name = args.object_name.c_str();
#ifdef _PYMOL_IP_PROPERTIES
  const char* object_props = args.object_props.c_str();
  const char* atom_props = args.atom_props.c_str();
#endif
  const char* plugin = args.plugin.c_str();
  auto content_format = args.content_format;
  auto state = args.state;
  auto zoom = args.zoom;
  auto discrete = args.discrete;
  auto finish = args.finish;
  auto multiplex = args.multiplex;
  auto quiet = args.quiet;

  if (multiplex != 1) {
    origObj = ExecutiveGetExistingCompatible(G, object_name, content_format);
  }

  // file type dependent multiplex and discrete default
  if(discrete < 0) {
    if(multiplex == 1) {
      discrete = 0;
    } else {
      switch (content_format) {
        case cLoadTypeMOL2:
        case cLoadTypeMOL2Str:
          discrete = -1;        /* content-dependent behavior... */
        case cLoadTypeSDF2:     /* SDF files currently default to discrete */
        case cLoadTypeSDF2Str:
        case cLoadTypeMAE:      /* Let MAE handle itself */
        case cLoadTypeMAEStr:
          break;
        default:
          discrete = 0;
          break;
      }
    }
  }

  int ok = true;
  OrthoLineType buf = "";
  int pdb_variant = PDB_VARIANT_DEFAULT;
  pymol::CObject *obj = NULL;

  // downstream file type reading functions
  switch (content_format) {
  case cLoadTypePQR:
    pdb_variant = PDB_VARIANT_PQR;
  case cLoadTypePDBQT:
    if (content_format == cLoadTypePDBQT)
      pdb_variant = PDB_VARIANT_PDBQT;
  case cLoadTypeVDBStr:
    if (ok && content_format == cLoadTypeVDBStr)
      pdb_variant = PDB_VARIANT_VDB;
  case cLoadTypePDB:
  case cLoadTypePDBStr:
    ok = ExecutiveProcessPDBFile(G, origObj, fname, content, object_name,
        state, discrete, finish, buf, pdb_variant,
        quiet, multiplex, zoom);
    break;
  case cLoadTypeCIF:
  case cLoadTypeCIFStr: {
    auto res =
        ObjectMoleculeReadCifStr(G, static_cast<ObjectMolecule*>(origObj),
            content, state, discrete, quiet, multiplex, zoom);
    p_return_if_error(res);
    obj = res.result();
  } break;
  case cLoadTypeMMTF:
  case cLoadTypeMMTFStr:
    obj = ObjectMoleculeReadMmtfStr(G, (ObjectMolecule *) origObj,
        content, size, state, discrete, quiet, multiplex, zoom);
    break;
  case cLoadTypeMAE:
  case cLoadTypeMAEStr:
#ifndef _PYMOL_IP_EXTRAS
    return pymol::Error::make<pymol::Error::INCENTIVE_ONLY>(
        "'mae' format not supported by this PyMOL build");
#else
#error ""
#endif
    break;
  case cLoadTypeTOP:
    if(origObj) {
      /* always reinitialize topology objects from scratch */
      ExecutiveDelete(G, origObj->Name);
      origObj = NULL;
    }
    obj = ObjectMoleculeLoadTOPFile(G, NULL, fname, state, discrete);
    break;
  case cLoadTypeTRJ:
    if(origObj) {
      ObjectMoleculeLoadTRJFile(G, (ObjectMolecule *) origObj, fname, state,
          1, 1, 1, -1, -1, NULL, 1, NULL, quiet);
    } else {
      return pymol::Error("must load object topology before loading trajectory!");
    }
    break;
  case cLoadTypeCRD:
    if(origObj) {
      ObjectMoleculeLoadRSTFile(G, (ObjectMolecule *) origObj, fname, state, quiet, 1);
    } else {
      return pymol::Error("must load object topology before loading coordinate file!");
    }
    break;
  case cLoadTypeRST:
    if(origObj) {
      ObjectMoleculeLoadRSTFile(G, (ObjectMolecule *) origObj, fname, state, quiet, 0);
    } else {
      return pymol::Error("must load object topology before loading restart file!");
    }
    break;
  case cLoadTypePMO:
    return pymol::Error("PMO format no longer supported.");
  case cLoadTypeDXStr:
    obj = ObjectMapReadDXStr(
        G, dynamic_cast<ObjectMap*>(origObj), content, size, state, quiet);
    break;
  case cLoadTypeDXMap:
    obj = ObjectMapLoadDXFile(G, (ObjectMap *) origObj, fname,
        state, quiet);
    break;
  case cLoadTypeFLDMap:
    obj = ObjectMapLoadFLDFile(G, (ObjectMap *) origObj, fname,
        state, quiet);
    break;
  case cLoadTypeBRIXMap:
    obj = ObjectMapLoadBRIXFile(G, (ObjectMap *) origObj, fname,
        state, quiet);
    break;
  case cLoadTypeGRDMap:
    obj = ObjectMapLoadGRDFile(G, (ObjectMap *) origObj, fname,
        state, quiet);
    break;
  case cLoadTypeACNTMap:
    obj = ObjectMapLoadACNTFile(G, (ObjectMap *) origObj, fname,
        state, quiet);
    break;
  case cLoadTypeXPLORMap:
  case cLoadTypeXPLORStr:
    obj = ObjectMapLoadXPLOR(G, (ObjectMap *) origObj, content,
        state, false, quiet);
    break;
  case cLoadTypePHIMap:
  case cLoadTypePHIStr:
    obj = ObjectMapLoadPHI(G, (ObjectMap *) origObj, content,
        state, true, size, quiet);
    break;
  case cLoadTypeCCP4Map:
  case cLoadTypeCCP4Str:
  case cLoadTypeCCP4Unspecified:
  case cLoadTypeCCP4UnspecifiedStr:
  case cLoadTypeMRC:
  case cLoadTypeMRCStr:
    obj = ObjectMapLoadCCP4(G, (ObjectMap *) origObj, content,
        state, true, size, quiet, content_format);
    break;
  case cLoadTypeCGO:
    obj = ObjectCGOFromFloatArray(G, (ObjectCGO *) origObj,
        (float *) content, size, state,
        quiet);
    break;
  case cLoadTypeMOL:
  case cLoadTypeMOLStr:
  case cLoadTypeMOL2:
  case cLoadTypeMOL2Str:
  case cLoadTypeSDF2:
  case cLoadTypeSDF2Str:
  case cLoadTypeXYZ:
  case cLoadTypeXYZStr:
  case cLoadTypeMMD:
  case cLoadTypeMMDStr:
    {
      const char * next_entry = content;
      char new_name[WordLength] = "";

      OVLexicon *loadproplex = NULL;
      bool loadpropertiesall = false;

      // (some of) these file types support multiple molecules per file,
      // and we support to load them into separate objects (multiplex).
      do {
        obj = ObjectMoleculeReadStr(G, (ObjectMolecule *) origObj,
            &next_entry, content_format,
            state, discrete,
            quiet, multiplex, new_name,
            loadpropertiesall, loadproplex);

        if(new_name[0]) {
          // multiplexing
          ObjectSetName(obj, new_name);
          ExecutiveDelete(G, obj->Name);       // just in case there is a collision
          ExecutiveManageObject(G, obj, zoom, true);
          new_name[0] = 0;
          obj = NULL;
        }
      } while(next_entry);

      OVLexicon_Del(loadproplex);
    }
    break;
  default:
    if(plugin[0]) {
      obj = PlugIOManagerLoad(G, origObj ? &origObj : NULL, fname, state, quiet,
          plugin, args.plugin_mask);
    } else {
        return pymol::make_error("Unable to read that file type from C (",
            content_format, ", ", plugin, ")");
    }
  }

  if(origObj && obj) {
    if(finish)
      ExecutiveUpdateObjectSelection(G, origObj);

    if(fname)
      sprintf(buf, " CmdLoad: \"%s\" appended into object \"%s\", state %d.\n",
          fname, object_name, state + 1);
  } else if(obj) {
    ObjectSetName(obj, object_name);
    ExecutiveManageObject(G, obj, zoom, true);

    if(fname)
      sprintf(buf, " CmdLoad: \"%s\" loaded as \"%s\".\n", fname, obj->Name);
    else
      sprintf(buf, " CmdLoad: loaded as \"%s\".\n", obj->Name);
  }

  if(!quiet && buf[0]) {
    PRINTFB(G, FB_Executive, FB_Actions)
      "%s", buf ENDFB(G);
  }

  return {};
}

/* ExecutiveGetExistingCompatible
 *
 * PARAMS
 *  oname -- object name
 *  type  -- new object type
 *
 * RETURNS
 *  Base-class object ptr
 *
 * SIDE EFFECTS
 *  If an object with the same name but different type already exists,
 *  then it is deleted.
 */
pymol::CObject* ExecutiveGetExistingCompatible(PyMOLGlobals * G, const char* oname, cLoadType_t type)
{
  pymol::CObject *origObj = NULL;
  origObj = ExecutiveFindObjectByName(G, oname);
  /* check for existing object of right type, delete if not */
  if(origObj) {
    int new_type = -1;
    switch (type) {
    case cLoadTypePlugin:
      // let PlugIOManager delete incompatible objects
      return origObj;
    case cLoadTypeChemPyModel:
    case cLoadTypePDB:
    case cLoadTypePDBStr:
    case cLoadTypeVDBStr:
    case cLoadTypeCIF:
    case cLoadTypeCIFStr:
    case cLoadTypeMMTF:
    case cLoadTypeMMTFStr:
    case cLoadTypeXYZ:
    case cLoadTypeXYZStr:
    case cLoadTypeMOL:
    case cLoadTypeMOLStr:
    case cLoadTypeMMD:
    case cLoadTypeMMDSeparate:
    case cLoadTypeMMDStr:
    case cLoadTypeTOP:
    case cLoadTypeTRJ:
    case cLoadTypeCRD:
    case cLoadTypeRST:
    case cLoadTypeMOL2:
    case cLoadTypeMOL2Str:
    case cLoadTypeSDF2:
    case cLoadTypeSDF2Str:
    case cLoadTypePQR:
    case cLoadTypePDBQT:
    case cLoadTypeXTC:
    case cLoadTypeDTR:
    case cLoadTypeTRR:
    case cLoadTypeGRO:
    case cLoadTypeTRJ2:
    case cLoadTypeG96:
    case cLoadTypeDCD:
      new_type = cObjectMolecule;
      break;
    case cLoadTypeChemPyBrick:
    case cLoadTypeChemPyMap:
    case cLoadTypeXPLORMap:
    case cLoadTypeXPLORStr:
    case cLoadTypeCCP4Map:
    case cLoadTypeCCP4Str:
    case cLoadTypeCCP4Unspecified:
    case cLoadTypeCCP4UnspecifiedStr:
    case cLoadTypeMRC:
    case cLoadTypeMRCStr:
    case cLoadTypeFLDMap:
    case cLoadTypeBRIXMap:
    case cLoadTypeGRDMap:
    case cLoadTypeDXMap:
    case cLoadTypeDXStr:
      new_type = cObjectMap;
      break;
    case cLoadTypeCallback:
      new_type = cObjectCallback;
      break;
    case cLoadTypeCGO:
      new_type = cObjectCGO;
      break;
    }
    if(new_type == -1 || new_type != origObj->type) {
      ExecutiveDelete(G, origObj->Name);
      origObj = NULL;
    }
  }
  return origObj;
}

int ExecutiveProcessPDBFile(PyMOLGlobals * G, pymol::CObject * origObj,
                            const char *fname, const char *buffer,
                            const char *oname, int frame, int discrete, int finish,
                            OrthoLineType buf, int variant, int quiet,
                            int multiplex, int zoom)
{
  int ok = true;
  pymol::CObject *obj;
  char pdb_name[WordLength] = "";
  char cur_name[WordLength] = "";
  const char *next_pdb = NULL;
  int repeat_flag = true;
  int n_processed = 0;
  PDBInfoRec pdb_info_rec, *pdb_info = NULL;
  int model_number;
  pymol::CObject *deferred_zoom_obj = NULL;

  UtilZeroMem(&pdb_info_rec, sizeof(PDBInfoRec));
  pdb_info = &pdb_info_rec;
  pdb_info->multiplex = multiplex;
  pdb_info->variant = variant;

  while(repeat_flag && ok) {
    const char *start_at = buffer;
    int is_repeat_pass = false;
    int eff_frame = frame;
    int is_new = false;

    if(next_pdb) {
      start_at = next_pdb;
      is_repeat_pass = true;
    }

    repeat_flag = false;
    next_pdb = NULL;
    if(!origObj) {

      is_new = true;
      pdb_name[0] = 0;
      model_number = 0;
      obj = ObjectMoleculeReadPDBStr(G, (ObjectMolecule *) origObj,
                                                 start_at, eff_frame, discrete,
                                                 pdb_name,
                                                 &next_pdb, pdb_info, quiet,
                                                 &model_number);

    } else {
      model_number = 0;
      ObjectMoleculeReadPDBStr(G, (ObjectMolecule *) origObj,
                               start_at, eff_frame, discrete,
                               pdb_name, &next_pdb, pdb_info, quiet, &model_number);

      if(finish) {
        ExecutiveUpdateObjectSelection(G, origObj);
        ExecutiveDoZoom(G, origObj, false, zoom, quiet);
      }
      if(eff_frame < 0)
        eff_frame = ((ObjectMolecule *) origObj)->NCSet - 1;
      if(buf) {
        if(fname)
          sprintf(buf, " CmdLoad: \"%s\" appended into object \"%s\", state %d.\n",
                  fname, oname, eff_frame + 1);
        else
          sprintf(buf, " CmdLoad: PDB-string appended into object \"%s\", state %d.\n",
                  oname, eff_frame + 1);
      }
      obj = origObj;
    }

    if(obj) {
      if(next_pdb) {
        /* NOTE: if set, assume that multiple PDBs are present in the file */
        repeat_flag = true;
      }
    }

    if(is_new) {
      if(obj) {
        if(next_pdb) {
          if(pdb_name[0] == 0) {
            if(cur_name[0]) {
              sprintf(pdb_name, "%s_%04d", cur_name, n_processed + 1);
            } else {
              sprintf(pdb_name, "%s_%04d", oname, n_processed + 1);
            }
          } else if(multiplex > 0) {
            if(pdb_info->multi_object_status == 1) {    /* this is a multi-object PDB file */
              strcpy(cur_name, pdb_name);
            } else if(cur_name[0] == 0) {
              strcpy(cur_name, oname);
            }
            if(model_number > 0) {
              sprintf(pdb_name, "%s_%04d", cur_name, model_number);
            } else {
              sprintf(pdb_name, "%s_%04d", cur_name, n_processed + 1);
            }
          }
          ObjectSetName(obj, pdb_name);
          ExecutiveDelete(G, obj->Name); /* just in case */
        } else {
          if(is_repeat_pass) {
            if(pdb_name[0] == 0) {
              if(cur_name[0]) {
                sprintf(pdb_name, "%s_%04d", cur_name, n_processed + 1);
              } else {
                sprintf(pdb_name, "%s_%04d", oname, n_processed + 1);
              }
            } else if(multiplex > 0) {
              if(pdb_info->multi_object_status == 1) {  /* this is a multi-object PDB file */
                strcpy(cur_name, pdb_name);
              } else if(cur_name[0] == 0) {
                strcpy(cur_name, oname);
              }
              if(model_number > 0) {
                sprintf(pdb_name, "%s_%04d", cur_name, model_number);
              } else {
                sprintf(pdb_name, "%s_%04d", cur_name, n_processed + 1);
              }
            }
            ObjectSetName(obj, pdb_name);       /* from PDB */
            ExecutiveDelete(G, obj->Name);       /* just in case */
          } else {
            ObjectSetName(obj, oname);  /* from filename/parameter */
          }
        }

        if(obj) {
          int do_zoom = repeat_flag ? 0 : zoom;
          if(do_zoom != zoom)
            deferred_zoom_obj = obj;
          else
            deferred_zoom_obj = NULL;
          ExecutiveManageObject(G, obj, do_zoom, true);
          if(eff_frame < 0)
            eff_frame = ((ObjectMolecule *) obj)->NCSet - 1;
          if(buf) {
            if(n_processed < 1) {
              if(fname)
                sprintf(buf, " CmdLoad: \"%s\" loaded as \"%s\".\n", fname, oname);
              else
                sprintf(buf,
                        " CmdLoad: PDB-string loaded into object \"%s\", state %d.\n",
                        oname, eff_frame + 1);
            } else {
              if(fname) {
                sprintf(buf, " CmdLoad: loaded %d objects from \"%s\".\n",
                        n_processed + 1, fname);
              } else {
                sprintf(buf, " CmdLoad: loaded %d objects from string.\n",
                        n_processed + 1);
              }
            }
          }

        }
      }
    }

    if(obj) {
      n_processed++;
    }
  }

  if(deferred_zoom_obj) {
    ExecutiveDoZoom(G, deferred_zoom_obj, true, zoom, true);
  }

  return ok;
}

pymol::Result<> ExecutiveAssignSS(PyMOLGlobals* G,
    const char* target, int state, const char* context, int preserve,
    ObjectMolecule* single_object, int quiet)
{
  if (!target[0]) {
    return pymol::make_error("selection must not be empty");
  }
  pymol::Result<SelectorTmp> targetTmp;
  pymol::Result<SelectorTmp> contextTmp;
  int sele0 = -1;
  int sele1 = -1;
  sele0 = SelectorIndexByName(G, target);
  if (sele0 < 0) {
    targetTmp = SelectorTmp::make(G, target);
    p_return_if_error(targetTmp);
    sele0 = targetTmp->getIndex();
    assert(sele0 >= 0);
  }
  if (!context || !context[0]) {
    sele1 = sele0;
  } else {
    contextTmp = SelectorTmp::make(G, context);
    p_return_if_error(contextTmp);
    sele1 = contextTmp->getIndex();
    assert(sele1 >= 0);
  }
  SelectorAssignSS(G, sele0, sele1, state, preserve, single_object, quiet);
  return {};
}

static int * getRepArrayFromBitmask(int visRep);

PyObject *ExecutiveGetVisAsPyDict(PyMOLGlobals * G)
{
  assert(PyGILState_Check());

  PyObject *result = NULL, *list;
  CExecutive *I = G->Executive;
  SpecRec *rec = NULL;
  result = PyDict_New();
  while(ListIterate(I->Spec, rec, next)) {
    if(rec->name[0] != '_') {
      list = PyList_New(4);
      PyList_SetItem(list, 0, PyInt_FromLong(rec->visible));

      PyList_SetItem(list, 1, PyList_New(0));

      if(rec->type != cExecObject) {
        PyList_SetItem(list, 2, PConvAutoNone(Py_None));
        PyList_SetItem(list, 3, PConvAutoNone(Py_None));
      } else {
        auto vla = getRepArrayFromBitmask(rec->obj->visRep);
        PyList_SetItem(list, 2, PConvIntVLAToPyList(vla));
        VLAFreeP(vla);

        PyList_SetItem(list, 3, PyInt_FromLong(rec->obj->Color));
      }

      PyDict_SetItemString(result, rec->name, list);
      Py_DECREF(list);
    }
  }
  return (result);
}

static int * getRepArrayFromBitmask(int visRep) {
  int n_vis = 0;
  int *RepVis = VLACalloc(int, cRepCnt);

  for(int a = 0; a < cRepCnt; a++)
    if(GET_BIT(visRep, a))
      RepVis[n_vis++] = a;

  VLASize(RepVis, int, n_vis);
  return RepVis;
}

#ifdef _PYMOL_LIB
/**
 * Returns a list (VLA) of enabled atom representations (AtomInfoType.visRep)
 * in selection (e.g. {cRepLine, cRepNonbonded})
 */
int *ExecutiveGetRepsInSceneForObject(PyMOLGlobals *G, const char *name){
  int visRep = 0;

  for (SeleAtomIterator iter(G, name); iter.next();) {
    AtomInfoType * ai = iter.getAtomInfo();
    visRep |= ai->visRep;
  }

  return getRepArrayFromBitmask(visRep);
}

/**
 * Returns a list (VLA) of enabled object representations (rec->obj->visRep)
 * (e.g. {cRepCell, cRepExtent})
 */
int *ExecutiveGetRepsForObject(PyMOLGlobals *G, const char *name){
  SpecRec *rec = ExecutiveFindSpec(G, (char*)name);

  if(rec)
    return getRepArrayFromBitmask(rec->obj->visRep);

  return NULL;
}
#endif

int ExecutiveSetVisFromPyDict(PyMOLGlobals * G, PyObject * dict)
{
#ifdef _PYMOL_NOPY
  return 0;
#else
  assert(PyGILState_Check());

  int ok = true;
  WordType name;
  PyObject *key, *list, *col;
  PyObject *vis_list = NULL;
  Py_ssize_t pos = 0;
  SpecRec *rec, *grec, **recstack = NULL;
  int n_vis;
  int rep;
  int a;
  int ll = 0;
  if(ok)
    ok = (dict != NULL);
  if(ok)
    ok = PyDict_Check(dict);
  if(ok) {

    SceneObjectDel(G, NULL, true);    /* remove all objects from scene */
    ExecutiveInvalidateSceneMembers(G);

    // stack for putative visible records
    recstack = pymol::calloc<SpecRec*>(PyDict_Size(dict) + 1);

    while(PyDict_Next(dict, &pos, &key, &list)) {
      if(!PConvPyStrToStr(key, name, sizeof(WordType))) {
        ok = false;
      } else {
        rec = ExecutiveFindSpec(G, name);
        if(rec) {
          if(ok)
            ok = (list != NULL);
          if(ok)
            ok = PyList_Check(list);
          if(ok)
            ll = PyList_Size(list);
          if(ok)
            ok = (ll >= 2);
          if(ok)
            ok = PConvPyObjectToInt(PyList_GetItem(list, 0), &rec->visible);

          /* before version 1.8 item 1 was rec reps (repOn) */

          if(ok && (rec->type == cExecObject)) {        /* object properties */

            if(ll > 2) {        /* object visibility */
              vis_list = PyList_GetItem(list, 2);
              if(ok)
                ok = (vis_list != NULL);
              if(ok) {
                if(PyList_Check(vis_list)) {
                  n_vis = PyList_Size(vis_list);
                  rec->obj->visRep = 0;
                  for(a = 0; a < n_vis; a++) {
                    if(PConvPyObjectToInt(PyList_GetItem(vis_list, a), &rep)) {
                      if((rep >= 0) && (rep < cRepCnt))
                        SET_BIT(rec->obj->visRep, rep);
                    }
                  }
                } else if (PyInt_Check(vis_list)) {
                  PConvPyObjectToInt(vis_list, &rec->obj->visRep);
                }
              }
            }
            if(ll > 3) {        /* object color */
              col = PyList_GetItem(list, 3);
              if(ok)
                ok = (col != NULL);
              if(ok)
                if(PyInt_Check(col)) {
                  ok = PConvPyObjectToInt(col, &rec->obj->Color);
                  rec->obj->invalidate(cRepAll, cRepInvColor, -1);
                }
            }
          }
          if(rec->visible && (rec->type == cExecObject)) {
            (*(++recstack)) = rec;
          }
        }
      }
    }

    // add visible objects to scene
    for(; (rec = *recstack); recstack--) {
      // check visibility of all parent groups
      for(grec = rec; grec->visible && (grec = grec->group););
      if(!grec) {
        // ok, no invisible parent found
        rec->in_scene = SceneObjectAdd(G, rec->obj);
        ExecutiveInvalidateSceneMembers(G);
      }
    }
    mfree(recstack);
  }
  return ok;
#endif
}

/**
 * returns a pointer to the data in a volume or map object
 */
CField * ExecutiveGetVolumeField(PyMOLGlobals * G, const char * objName, int state) {
  ObjectMapState *oms;
  pymol::CObject *obj;

  obj = ExecutiveFindObjectByName(G, objName);
  ok_assert(1, obj);

  switch (obj->type) {
  case cObjectVolume:
    return ObjectVolumeGetField((ObjectVolume *) obj);
  case cObjectMap:
    oms = ObjectMapGetState((ObjectMap *) obj, state);
    ok_assert(1, oms && oms->Field);
    return oms->Field->data.get();
  }

ok_except1:
  return NULL;
}

/**
 * returns allocated memory
 */
pymol::Result<std::vector<float>>
ExecutiveGetHistogram(PyMOLGlobals * G, const char * objName, int n_points, float min_val, float max_val) {
  pymol::CObject *obj;
  ObjectMapState *oms = NULL;

  obj = ExecutiveFindObjectByName(G, objName);

  if (!obj) {
    return pymol::make_error("could not find object ", objName);
  }

  switch (obj->type) {
  case cObjectMap:
    oms = ObjectMapGetState((ObjectMap *) obj, 0);
    break;
  case cObjectVolume:
    oms = ObjectVolumeGetMapState((ObjectVolume *) obj);
    break;
  default:
    return pymol::make_error("object type must be map or volume");
  }

  if(oms) {
    auto hist = std::vector<float>(n_points + 4);
    float range = SettingGet_f(G, obj->Setting.get(), NULL, cSetting_volume_data_range);
    ObjectMapStateGetHistogram(
        G, oms, n_points, range, hist.data(), min_val, max_val);
    return hist;
  }

  return pymol::make_error("failed to get map state");
}

PyObject* ExecutiveGetVolumeRamp(PyMOLGlobals * G, const char * objName) {
#ifdef _PYMOL_NOPY
  return NULL;
#else
  pymol::CObject *obj;
  PyObject* result = NULL;

  PRINTFD(G, FB_Executive) "Executive-GetVolumeRamp Entered.\n" ENDFD;

  obj = ExecutiveFindObjectByName(G, objName);
  if(obj && obj->type==cObjectVolume) {
    result = ObjectVolumeGetRamp((ObjectVolume *) obj);
  }

  PRINTFD(G, FB_Executive) "Executive-GetVolumeRamp Exited.\n" ENDFD;

  return result;

#endif
}

pymol::Result<> ExecutiveSetVolumeRamp(PyMOLGlobals * G, const char * objName, std::vector<float> ramp_list) {

  auto obj = ExecutiveFindObject<ObjectVolume>(G, objName);
  if(obj) {
    return ObjectVolumeSetRamp(obj, std::move(ramp_list));
  }

  return pymol::make_error("Object ", objName, " not found");
}

pymol::Result<> ExecutiveIsolevel(
    PyMOLGlobals* G, const char* name, float level, int state, int quiet)
{
  auto obj = ExecutiveFindObjectByName(G, name);
  if(obj) {
    switch (obj->type) {
    case cObjectMesh:
      ObjectMeshSetLevel((ObjectMesh *) obj, level, state, quiet);
      SceneChanged(G);
      return {};
    case cObjectSurface:
      ObjectSurfaceSetLevel((ObjectSurface *) obj, level, state, quiet);
      SceneChanged(G);
      return {};
    default:
      return pymol::make_error("Object ", name, " is of wrong type.");
    }
  }
  return pymol::make_error("Object not found");
}

pymol::Result<float> ExecutiveGetIsolevel(
    PyMOLGlobals* G, const char* name, int state)
{
  auto obj = ExecutiveFindObjectByName(G, name);
  if (obj) {
    switch (obj->type) {
    case cObjectMesh:
      return ObjectMeshGetLevel((ObjectMesh*) obj, state);
    case cObjectSurface:
      return ObjectSurfaceGetLevel((ObjectSurface*) obj, state);
    default:
      return pymol::make_error("Object ", name, " is of wrong type.");
    }
  }
  return pymol::make_error("Object not found");
}

pymol::Result<std::pair<float, float>> ExecutiveSpectrum(PyMOLGlobals* G,
    pymol::zstring_view s_s1, pymol::zstring_view s_expr, float min, float max, int first, int last,
    pymol::zstring_view s_prefix, int digits, int byres, int quiet)
{
  std::pair<float, float> ret;
  int n_color, n_atom;
  ObjectMoleculeOpRec op;
  WordType buffer;
  std::vector<int> color_index;
  std::vector<float> value;
  int a, b;
  char pat[] = "%0Xd";
  int pref_len;
  char *at;
  float range;

  auto s1 = s_s1.c_str();
  auto prefix = s_prefix.c_str();
  auto expr = s_expr.c_str();

  auto tmpsele1 = SelectorTmp::make(G, s1);
  p_return_if_error(tmpsele1);
  int sele1 = tmpsele1->getIndex();

  if(sele1 >= 0) {

    if(digits > 9)
      digits = 9;
    pat[2] = ('0' + digits);
    UtilNCopy(buffer, prefix, sizeof(WordType) - digits);

    pref_len = strlen(prefix);
    at = buffer + pref_len;

    n_color = abs(first - last) + 1;
    if(n_color) {
      color_index.resize(n_color);
      for(a = 0; a < n_color; a++) {
        b = first + ((last - first) * a) / (n_color - 1);
        sprintf(at, pat, b);
        color_index[a] = ColorGetIndex(G, buffer);
      }

      // set up iterator
      SeleAtomIterator iter(G, sele1);
      SelectorUpdateTable(G, cSelectorUpdateTableAllStates, -1);

      // count atoms
      for (n_atom = 0; iter.next();) {
        ++n_atom;
      }

      if(n_atom) {
        value.resize(n_atom);

        if(WordMatchExact(G, "count", expr, true)) {
          for(a = 0; a < n_atom; a++) {
            value[a] = (float) a + 1;
          }
        } else {
          if (WordMatchExact(G, "pc", expr, true)) {
            expr = "partial_charge";
          } else if (WordMatchExact(G, "resi", expr, true)) {
            expr = "resv";
          }

          // look up expression definition
          auto ap = PyMOL_GetAtomPropertyInfo(G->PyMOL, expr);
          if (!ap) {
            return pymol::make_error("Unknown expression: ", expr);
          }

          // for enumerated values
          std::map<size_t, unsigned int> enumerated_values;
          union {
            size_t value_e;
            char   value_s[sizeof(size_t)];
          };

          for (a = 0, iter.reset(); iter.next(); ++a) {
            const auto ai = iter.getAtomInfo();
            const auto raw_ptr = reinterpret_cast<const char*>(ai) + ap->offset;

            // numeric values
            switch (ap->Ptype) {
              case cPType_float:
                value[a] = *reinterpret_cast<const float*>(raw_ptr);
                continue;
              case cPType_int:
              case cPType_int_custom_type:
                value[a] = *reinterpret_cast<const int*>(raw_ptr);
                continue;
              case cPType_uint32:
                value[a] = *reinterpret_cast<const uint32_t*>(raw_ptr);
                continue;
              case cPType_schar:
                value[a] = *reinterpret_cast<const signed char*>(raw_ptr);
                continue;
              case cPType_char_as_type:
                value[a] = ai->hetatm;
                continue;
              case cPType_index:
                value[a] = iter.getAtm() + 1.f;
                continue;
            }

            // enumerated values
            switch (ap->Ptype) {
              case cPType_int_as_string:
                value_e = LexNumeric(*reinterpret_cast<const lexidx_t*>(raw_ptr));
                break;
              case cPType_string:
                // works for small strings
                strncpy(value_s, raw_ptr, sizeof(value_e));
                break;
              case cPType_model:
                value_e = (size_t) iter.obj;
                break;
              default:
                return pymol::make_error("Unsupported Ptype for expr: ", expr);
            }

            // lookup or insert value
            auto& e = enumerated_values[value_e];
            if (e == 0)
              e = enumerated_values.size();
            value[a] = e - 1.f;
          }

          if (!quiet && !enumerated_values.empty()) {
            PRINTFB(G, FB_Executive, FB_Actions)
              " Spectrum: Expression is non-numeric, enumerating values\n" ENDFB(G);
          }
        }

        if(max < min) {
          max = value[0];
          min = value[0];
          for(a = 1; a < n_atom; a++) {
            if(value[a] < min)
              min = value[a];
            if(value[a] > max)
              max = value[a];
          }
        }
        range = max - min;

        if(!quiet) {
          PRINTFB(G, FB_Executive, FB_Actions)
            " Spectrum: range (%8.5f to %8.5f).\n", min, max ENDFB(G);
        }
        if(range == 0.0F)
          range = 1.0F;
        ret.first = min;
        ret.second = max;

        op.code = OMOP_Spectrum;
        op.i1 = n_color - 1;
        op.i2 = n_atom;
        op.i3 = 0;
        op.i4 = byres;
        op.ii1 = color_index.data();
        op.ff1 = value.data();
        op.f1 = min;
        op.f2 = range;

        ExecutiveObjMolSeleOp(G, sele1, &op);

        op.code = OMOP_INVA;
        op.i1 = cRepBitmask;
        op.i2 = cRepInvColor;
        ExecutiveObjMolSeleOp(G, sele1, &op);

      }
    }
  }
  return ret;
}

static int fStrOrderFn(const char * const* array, int l, int r) {
  return strcmp(array[l], array[r]) < 0;
}

/**
 * Returns an VLA with pointers into G->Lexicon
 */
pymol::Result<std::vector<const char*>> ExecutiveGetChains(
    PyMOLGlobals* G, const char* s1, int state)
{
  std::set<lexidx_t> chains;
  int c = 0;
  ObjectMoleculeOpRec op;

  SETUP_SELE_DEFAULT(1);

  {

    ObjectMoleculeOpRecInit(&op);
    op.code = OMOP_GetChains;
    op.ii1 = (int*) (void*) &chains; // pointer pack
    op.i1 = 0;
    ExecutiveObjMolSeleOp(G, sele1, &op);
    std::vector<const char*> result(chains.size());
    for (const auto& chain : chains) {
      result[c++] = LexStr(G, chain);
    }
    // sort the array
    UtilSortInPlace(G, result.data(), chains.size(), sizeof(char *),
        (UtilOrderFn *) fStrOrderFn);
    return result;
  }
}

int ExecutiveValidateObjectPtr(PyMOLGlobals * G, pymol::CObject * ptr, int object_type)
{
  /* this routine needs to be sped up significantly... */

  CExecutive *I = G->Executive;
  int ok = false;
  SpecRec *rec = NULL;

  while(ListIterate(I->Spec, rec, next)) {
    if(rec->obj == ptr) {
      if(rec->type == cExecObject) {
        if((!object_type) || (rec->obj->type == object_type)) {
          ok = true;
          break;
        }
      }
    }
  }
  return (ok);
}

pymol::Result<> ExecutiveRampNew(PyMOLGlobals* G, const char* name,
    const char* src_name, pymol::vla<float> range, pymol::vla<float> color,
    int src_state, const char* sele, float beyond, float within, float sigma,
    int zero, int calc_mode, int quiet)
{
  ObjectGadgetRamp *obj = NULL;
  ObjectGadgetRamp *origRamp = NULL;
  pymol::CObject *src_obj = NULL;
  pymol::CObject *origObj = ExecutiveFindObjectByName(G, name);
  float *vert_vla = NULL;
  int rampType = -1;

  if (origObj &&
      origObj->type == cObjectGadget &&
      ((ObjectGadget*)origObj)->GadgetType == cGadgetRamp) {
    origRamp = (ObjectGadgetRamp*)origObj;
    rampType = origRamp->RampType;
  } else if (!range || !(color || calc_mode)) {
     return pymol::make_error("Missing 'range' or 'color' to create new ramp.");
  }

  if (src_name && src_name[0]) {
    if (WordMatchExact(G, src_name, cKeywordNone, true)) {
      rampType = cRampNone;
    } else {
      src_obj = ExecutiveFindObjectByName(G, src_name);
      if(src_obj) {
        switch (src_obj->type) {
          case cObjectMap:
            rampType = cRampMap;
            break;
          case cObjectMolecule:
            rampType = cRampMol;
            break;
          default:
            pymol::make_error(src_name, " is not a map or molecule.");
        }
      } else {
        return pymol::make_error(src_name, " not found.");
      }
    }
  }

  switch (rampType) {
    case cRampMap:
      /* mapping this ramp from a selection */
      if(sele && sele[0]) {
        auto tmpsele = SelectorTmp::make(G, sele);
        p_return_if_error(tmpsele);
        sele = tmpsele->getName();
        assert(sele[0]);

        vert_vla = ExecutiveGetVertexVLA(G, sele, src_state);
      }
      obj = ObjectGadgetRampMapNewAsDefined(G, origRamp, (ObjectMap *) src_obj,
            std::move(range), std::move(color), src_state,
            vert_vla, beyond, within,
            sigma, zero, calc_mode);
      VLAFreeP(vert_vla);
      break;
    case cRampNone:
    case cRampMol:
      obj = ObjectGadgetRampMolNewAsDefined(G, origRamp, (ObjectMolecule *) src_obj,
            std::move(range), std::move(color), src_state,
            calc_mode);
      break;
    default:
      return pymol::make_error("Missing 'name' to create new ramp.");
  }

  if (!obj)
    return pymol::make_error("Object not found");

  if (obj != origRamp) {
    ExecutiveDelete(G, name);
    ObjectSetName(obj, name);
    ColorRegisterExt(G, obj->Name, obj);
    ExecutiveManageObject(G, obj, false, quiet);
  }

  ExecutiveInvalidateRep(G, cKeywordAll, cRepAll, cRepInvColor);      /* recolor everything */
  return {};
}

static int ExecutiveSetNamedEntries(PyMOLGlobals * G, PyObject * names, int version,
                                    int part_rest, int part_sess)
{
  CExecutive *I = G->Executive;
  int ok = true;
  int skip = false;
  int a = 0, l = 0, ll = 0;
  PyObject *cur, *el;
  SpecRec *rec = NULL;
  int extra_int;
  int incomplete = false;
  ObjectNameType new_name;

  if(ok)
    ok = (names != NULL);
  if(ok)
    ok = PyList_Check(names);
  if(ok)
    l = PyList_Size(names);

  while(ok && (a < l)) {
    cur = PyList_GetItem(names, a);
    if(cur != Py_None) {        /* skip over None w/o aborting */
      skip = false;
      rec = NULL;
      ListElemCalloc(G, rec, SpecRec);
      rec->next = NULL;
      rec->name[0] = 0;
      if(ok)
        ok = PyList_Check(cur);
      if(ok)
        ll = PyList_Size(cur);
      if(ok)
        ok = PConvPyStrToStr(PyList_GetItem(cur, 0), rec->name, sizeof(WordType));
      if(ok)
        ok = PConvPyIntToInt(PyList_GetItem(cur, 1), &rec->type);
      if(ok)
        ok = CPythonVal_PConvPyIntToInt_From_List(G, cur, 2, &rec->visible);

      /* before version 1.8 item 3 was rec reps (repOn) */

      if(ok)
        ok = PConvPyIntToInt(PyList_GetItem(cur, 4), &extra_int);
      switch (rec->type) {
      case cExecObject:
        if(!ok)
          break;

        el = PyList_GetItem(cur, 5);

        switch (extra_int) {
        case cObjectMolecule:
          ok = ObjectMoleculeNewFromPyList(G, el, (ObjectMolecule **) (void *) &rec->obj);
          break;
        case cObjectMeasurement:
          ok = ObjectDistNewFromPyList(G, el, (ObjectDist **) (void *) &rec->obj);
          break;
        case cObjectMap:
          ok = ObjectMapNewFromPyList(G, el, (ObjectMap **) (void *) &rec->obj);
          break;
        case cObjectMesh:
          ok = ObjectMeshNewFromPyList(G, el, (ObjectMesh **) (void *) &rec->obj);
          break;
        case cObjectSlice:
          ok = ObjectSliceNewFromPyList(G, el, (ObjectSlice **) (void *) &rec->obj);
          break;
        case cObjectSurface:
          ok = ObjectSurfaceNewFromPyList(G, el, (ObjectSurface **) (void *) &rec->obj);
          break;
        case cObjectCGO:
          ok = ObjectCGONewFromPyList(G, el, (ObjectCGO **) (void *) &rec->obj, version);
          break;
        case cObjectGadget:
          ok = ObjectGadgetNewFromPyList(G, el, (ObjectGadget **) (void *) &rec->obj, version);
          break;
        case cObjectAlignment:
          ok = ObjectAlignmentNewFromPyList(G, el, (ObjectAlignment **) (void *) &rec->obj,
                                              version);
          break;
        case cObjectGroup:
          if(part_rest) {
            // if group already exists, do not create new one
            pymol::CObject *obj = ExecutiveFindObjectByName(G, rec->name);
            if(obj && obj->type == cObjectGroup) {
              skip = 1;
              break;
            }
          }
          ok = ObjectGroupNewFromPyList(G, el, (ObjectGroup **) (void *) &rec->obj, version);
          break;
        case cObjectVolume:
          ok = ObjectVolumeNewFromPyList(G, el, (ObjectVolume **) (void *) &rec->obj);
          break;
#ifndef _PYMOL_NOPY
        case cObjectCallback:
          // skip dummy entries from old sessions and failed-to-pickle sessions
          skip = !ObjectCallbackNewFromPyList(G, el, (ObjectCallback **) (void *) &rec->obj);
          break;
#endif
        default:
          PRINTFB(G, FB_Executive, FB_Errors)
            " Executive: skipping unrecognized object \"%s\" of type %d.\n",
            rec->name, extra_int ENDFB(G);
          skip = true;
          break;
        }

        CPythonVal_Free(el);

        break;
      case cExecSelection:     // on the first pass, just create an entry in the rec list
        rec->sele_color = extra_int;
        if(part_rest || part_sess) {    // don't attempt to restore selections with partial sessions
          skip = true;
        }
        break;
      }

      if(ll > 6) {
        if(ok){
          ok = PConvPyStrToStr(PyList_GetItem(cur, 6), rec->group_name, sizeof(WordType));
	}
      }

      if(PyErr_Occurred()) {
        PRINTFB(G, FB_Executive, FB_Errors)
          "ExectiveSetNamedEntries-Error: after object \"%s\".\n", rec->name ENDFB(G);
        PyErr_Print();
      }

      if(ok && !skip) {
        if(part_rest && ExecutiveProcessObjectName(G, rec->name, new_name)) {
          // rename duplicates
          strcpy(rec->obj->Name, new_name);
          strcpy(rec->name, new_name);
        }

        // replace existing object (unless auto_rename_duplicate_objects=1)
        if(ExecutiveValidName(G, rec->name)) {
          ExecutiveDelete(G, rec->name);
        }

        switch (rec->type) {
        case cExecObject:
          if(rec->visible) {
            rec->in_scene = SceneObjectAdd(G, rec->obj);
            ExecutiveInvalidateSceneMembers(G);
          }
          ExecutiveUpdateObjectSelection(G, rec->obj);
          break;
        }

        rec->cand_id = TrackerNewCand(I->Tracker, (TrackerRef *) rec);
        TrackerLink(I->Tracker, rec->cand_id, I->all_names_list_id, 1);

        switch (rec->type) {
        case cExecObject:
          TrackerLink(I->Tracker, rec->cand_id, I->all_obj_list_id, 1);
          break;
        case cExecSelection:
          TrackerLink(I->Tracker, rec->cand_id, I->all_sel_list_id, 1);
          break;
        }
        ListAppend(I->Spec, rec, next, SpecRec);
        ExecutiveAddKey(I, rec);
        ExecutiveInvalidateGroups(G, false);
        ExecutiveInvalidatePanelList(G);
      } else {
        ListElemFree(rec);
      }
    }
    a++;
    if(!ok) {
      incomplete = true;
      ok = true;
    }
  }
  return (!incomplete);
}

static int ExecutiveSetSelectionsFromPyList(PyMOLGlobals * G, PyObject * names)
{
  /* must already have objects loaded at this point... */

  int ok = true;
  int a = 0, l = 0;
  PyObject *cur;
  SpecRec *rec = NULL;
  int extra;
  int incomplete = false;

  if(ok)
    ok = (names != NULL);
  if(ok)
    ok = PyList_Check(names);
  if(ok)
    l = PyList_Size(names);
  while(ok && (a < l)) {
    cur = PyList_GetItem(names, a);
    if(cur != Py_None) {        /* skip over None w/o aborting */
      rec = NULL;
      ListElemCalloc(G, rec, SpecRec);
      rec->next = NULL;

      if(ok)
        ok = PyList_Check(cur);
      if(ok)
        ok = PConvPyStrToStr(PyList_GetItem(cur, 0), rec->name, sizeof(WordType));
      if(ok)
        ok = PConvPyIntToInt(PyList_GetItem(cur, 1), &rec->type);
      if(ok)
        ok = PConvPyIntToInt(PyList_GetItem(cur, 2), &rec->visible);
      if(ok)
        ok = PConvPyIntToInt(PyList_GetItem(cur, 4), &extra);
      switch (rec->type) {
      case cExecSelection:
        ok = SelectorFromPyList(G, rec->name, PyList_GetItem(cur, 5));
        break;
      }
      ListElemFree(rec);
    }
    a++;
    if(!ok) {
      incomplete = true;
      ok = true;
    }
  }
  return (!incomplete);
}

static PyObject *ExecutiveGetExecObjectAsPyList(PyMOLGlobals * G, SpecRec * rec)
{

  PyObject *result = NULL;
  int recobjtype = rec->obj->type;
  switch (recobjtype){
  case cObjectMesh:
    { /* If a mesh no longer has its dependent map, then it gets saved as a CGO */
      int allMapsExist = ObjectMeshAllMapsInStatesExist((ObjectMesh *) rec->obj);
      if (!allMapsExist){
	recobjtype = cObjectCGO;
      }
    }
  }
  result = PyList_New(7);
  PyList_SetItem(result, 0, PyString_FromString(rec->obj->Name));
  PyList_SetItem(result, 1, PyInt_FromLong(cExecObject));
  PyList_SetItem(result, 2, PyInt_FromLong(rec->visible));
  /* before version 1.8 item 3 was rec reps (repOn) */
  PyList_SetItem(result, 3, PConvAutoNone(NULL));
  PyList_SetItem(result, 4, PyInt_FromLong(recobjtype));
  switch (rec->obj->type) {
  case cObjectGadget:
    PyList_SetItem(result, 5, ObjectGadgetAsPyList((ObjectGadget *) rec->obj));
    break;
  case cObjectMolecule:
    PyList_SetItem(result, 5, ObjectMoleculeAsPyList((ObjectMolecule *) rec->obj));
    break;
  case cObjectMeasurement:
    PyList_SetItem(result, 5, ObjectDistAsPyList((ObjectDist *) rec->obj));
    break;
  case cObjectMap:
    PyList_SetItem(result, 5, ObjectMapAsPyList((ObjectMap *) rec->obj));
    break;
  case cObjectMesh:
    PyList_SetItem(result, 5, ObjectMeshAsPyList((ObjectMesh *) rec->obj));
    break;
  case cObjectSlice:
    PyList_SetItem(result, 5, ObjectSliceAsPyList((ObjectSlice *) rec->obj));
    break;
  case cObjectSurface:
    PyList_SetItem(result, 5, ObjectSurfaceAsPyList((ObjectSurface *) rec->obj));
    break;
  case cObjectCGO:
    PyList_SetItem(result, 5, ObjectCGOAsPyList((ObjectCGO *) rec->obj));
    break;
  case cObjectAlignment:
    PyList_SetItem(result, 5, ObjectAlignmentAsPyList((ObjectAlignment *) rec->obj));
    break;
  case cObjectGroup:
    PyList_SetItem(result, 5, ObjectGroupAsPyList((ObjectGroup *) rec->obj));
    break;
  case cObjectVolume:
    PyList_SetItem(result, 5, ObjectVolumeAsPyList((ObjectVolume *) rec->obj));
    break;
  case cObjectCallback:
    PyList_SetItem(result, 5, ObjectCallbackAsPyList((ObjectCallback *) rec->obj));
    break;
  default:
    PyList_SetItem(result, 5, PConvAutoNone(NULL));
    break;
  }
  PyList_SetItem(result, 6, PyString_FromString(rec->group_name));

  return (result);
}

static PyObject *ExecutiveGetExecSeleAsPyList(PyMOLGlobals * G, SpecRec * rec)
{
  PyObject *result = NULL;
  int sele;

  sele = SelectorIndexByName(G, rec->name);
  if(sele >= 0) {
    result = PyList_New(7);
    PyList_SetItem(result, 0, PyString_FromString(rec->name));
    PyList_SetItem(result, 1, PyInt_FromLong(cExecSelection));
    PyList_SetItem(result, 2, PyInt_FromLong(rec->visible));
    /* before version 1.8 item 3 was rec reps (repOn) */
    PyList_SetItem(result, 3, PConvAutoNone(NULL));
    PyList_SetItem(result, 4, PyInt_FromLong(-1));
    PyList_SetItem(result, 5, SelectorAsPyList(G, sele));
    PyList_SetItem(result, 6, PyString_FromString(rec->group_name));

  }
  return (PConvAutoNone(result));
}

static PyObject *ExecutiveGetNamedEntries(PyMOLGlobals * G, int list_id, int partial)
{
  CExecutive *I = G->Executive;
  CTracker *I_Tracker = I->Tracker;
  PyObject *result = NULL;
  int count = 0, total_count = 0;
  int iter_id = 0;
  SpecRec *rec = NULL, *list_rec = NULL;

  SelectorUpdateTable(G, cSelectorUpdateTableAllStates, -1);

  if(list_id) {
    total_count = TrackerGetNCandForList(I_Tracker, list_id);
    iter_id = TrackerNewIter(I_Tracker, 0, list_id);
  } else {
    total_count = ExecutiveCountNames(G);
  }
  result = PyList_New(total_count);

  /* critical reliance on short-circuit behavior */

  while((iter_id && TrackerIterNextCandInList(I_Tracker, iter_id,
                                              (TrackerRef **) (void *) &list_rec)) ||
        ((!iter_id) && ListIterate(I->Spec, rec, next))) {

    if(list_id)
      rec = list_rec;
    if(count >= total_count)
      break;
    if(rec) {
      switch (rec->type) {
      case cExecObject:
        PyList_SetItem(result, count, ExecutiveGetExecObjectAsPyList(G, rec));
        break;
      case cExecSelection:
        if(!partial) {
          PyList_SetItem(result, count, ExecutiveGetExecSeleAsPyList(G, rec));
        } else {
          /* cannot currently save selections in partial sessions */
          PyList_SetItem(result, count, PConvAutoNone(NULL));
        }
        break;
      default:
        PyList_SetItem(result, count, PConvAutoNone(NULL));
        break;
      }
    } else {
      PyList_SetItem(result, count, PConvAutoNone(NULL));
    }
    count++;
  }

  while(count < total_count) {  /* insure that all members of outgoing list are defined */
    PyList_SetItem(result, count, PConvAutoNone(NULL));
    count++;
  }

  if(iter_id) {
    TrackerDelIter(I_Tracker, iter_id);
  }
  return (PConvAutoNone(result));
}

#ifdef PYMOL_EVAL
#include "ExecutiveEvalMessage.h"
#endif

int ExecutiveGetSession(PyMOLGlobals * G, PyObject * dict, const char *names, int partial,
                        int quiet)
{
  assert(PyGILState_Check());

  int list_id = 0;
  SceneViewType sv;
  PyObject *tmp;

  if(names && names[0]) {
    list_id = ExecutiveGetNamesListFromPattern(G, names, true, cExecExpandKeepGroups);
  }

  tmp = MovieScenesAsPyList(G);
  PyDict_SetItemString(dict, "moviescenes", tmp);
  Py_XDECREF(tmp);

  tmp = PyInt_FromLong(_PyMOL_VERSION_int);
  PyDict_SetItemString(dict, "version", tmp);
  Py_XDECREF(tmp);

  tmp = ExecutiveGetNamedEntries(G, list_id, partial);
  PyDict_SetItemString(dict, "names", tmp);
  Py_XDECREF(tmp);

  tmp = ColorAsPyList(G);
  PyDict_SetItemString(dict, "colors", tmp);
  Py_XDECREF(tmp);

  tmp = ColorExtAsPyList(G);
  PyDict_SetItemString(dict, "color_ext", tmp);
  Py_XDECREF(tmp);

  tmp = SettingUniqueAsPyList(G);
  PyDict_SetItemString(dict, "unique_settings", tmp);
  Py_XDECREF(tmp);

  if(partial) {                 /* mark this as a partial session */

    PyDict_SetItemString(dict, "partial", PConvAutoNone(Py_None));

  } else {

    /* none of the following information is saved in partial sessions */

    tmp = SelectorSecretsAsPyList(G);
    PyDict_SetItemString(dict, "selector_secrets", tmp);
    Py_XDECREF(tmp);

    tmp = SettingGetGlobalsAsPyList(G);
    PyDict_SetItemString(dict, "settings", tmp);
    Py_XDECREF(tmp);

    SceneGetView(G, sv);
    tmp = PConvFloatArrayToPyList(sv, cSceneViewSize);
    PyDict_SetItemString(dict, "view", tmp);
    Py_XDECREF(tmp);

    tmp = MovieAsPyList(G);
    PyDict_SetItemString(dict, "movie", tmp);
    Py_XDECREF(tmp);

    tmp = EditorAsPyList(G);
    PyDict_SetItemString(dict, "editor", tmp);
    Py_XDECREF(tmp);

    tmp = MainAsPyList(G);
    PyDict_SetItemString(dict, "main", tmp);
    Py_XDECREF(tmp);

#ifdef PYMOL_EVAL
    ExecutiveEvalMessage(G, dict);
#endif

  }

  return true;
}

static void ExecutiveMigrateSession(PyMOLGlobals * G, int session_version)
{
  if (session_version < 1700) {
    if (SettingGetGlobal_i(G, cSetting_seq_view_label_color) == 0 /* white */) {
      SettingSetGlobal_i(G, cSetting_seq_view_label_color, cColorFront);
    }
  }
  if(session_version < 100) {
    /* migrate lighting model */
    SettingSetGlobal_f(G, cSetting_direct, 1.8 * SettingGetGlobal_f(G, cSetting_direct));
    SettingSetGlobal_f(G, cSetting_reflect,
                       0.5 * SettingGetGlobal_f(G, cSetting_reflect));
    SettingSetGlobal_f(G, cSetting_ambient,
                       1.166 * SettingGetGlobal_f(G, cSetting_ambient));
    SettingSetGlobal_f(G, cSetting_gamma, 0.769 * SettingGetGlobal_f(G, cSetting_gamma));

    /* try best to meet existing expectations with existing sessions */
    SettingSetGlobal_f(G, cSetting_ray_legacy_lighting, 1.0F);

    /* force use of movie_delay in preference to movie_fps */

    SettingSetGlobal_f(G, cSetting_movie_fps, 0.0F);

    /* and labels */
    SettingSetGlobal_i(G, cSetting_label_digits, 2);
  }
  if(session_version < 99) {
    SettingSetGlobal_f(G, cSetting_cartoon_ladder_mode, 0);
    SettingSetGlobal_f(G, cSetting_cartoon_tube_cap, 0);
    SettingSetGlobal_f(G, cSetting_cartoon_nucleic_acid_mode, 1);
    {
      float old_sulfur[3] = { 1.0, 0.5, 0.0 };
      ColorDef(G, "sulfur", old_sulfur, 0, true);
    }
  }
  if(session_version < 98) {
    /* produce expected rendering quality & performance with old sessions */
    SettingSetGlobal_b(G, cSetting_ray_orthoscopic, 1);
  }
  if(session_version < 96) {
    SettingSetGlobal_f(G, cSetting_ray_transparency_contrast, 1.0F);
  }
  if(session_version < 95) {
    {
      /* adjust fog to reflect current importance of seeing to the Z-slab center w/o fog */

      float fog_start = SettingGetGlobal_f(G, cSetting_fog_start);
      float ray_trace_fog_start = SettingGetGlobal_f(G, cSetting_ray_trace_fog_start);
      if((fog_start == 0.40F) || (fog_start == 0.35F) || (fog_start == 0.30F)) {
        SettingSetGlobal_f(G, cSetting_fog_start, 0.45F);
      }
      if((ray_trace_fog_start == 0.45F) || (ray_trace_fog_start == 0.40F)
         || (ray_trace_fog_start == 0.35F)) {
        SettingSetGlobal_f(G, cSetting_ray_trace_fog_start, 0.50F);
      }
    }

    {                           /* adjust GUI width */

      int gui_width = SettingGetGlobal_i(G, cSetting_internal_gui_width);

      if(gui_width == 160) {
        SettingSetGlobal_i(G, cSetting_internal_gui_width, 220);
      }
    }

    {                           /* enable antialiasing */

      int antialias = SettingGetGlobal_i(G, cSetting_antialias);

      if(antialias == 0) {
        SettingSetGlobal_i(G, cSetting_antialias, 1);
      }

    }
  }
}

int ExecutiveSetSessionNoMLock(PyMOLGlobals* G, PyObject* session)
{
  auto mLocked = MovieLocked(G);
  auto res = ExecutiveSetSession(G, session, false, true);
  MovieSetLock(G, mLocked);
  return res;
}

int ExecutiveSetSession(PyMOLGlobals * G, PyObject * session,
                        int partial_restore, int quiet)
{
  assert(PyGILState_Check());

  int ok = true;
  int incomplete = false;
  PyObject *tmp;
  SceneViewType sv;
  int version = -1, version_full;
  int migrate_sessions = SettingGetGlobal_b(G, cSetting_session_migration);
  char active[WordLength] = "";
  int have_active = false;
  int partial_session = false;

  if(!partial_restore) {        /* if user has requested partial restore */
    ExecutiveDelete(G, "all");
    ColorReset(G);
  }

  if(!session || !PyDict_Check(session)) {
    PRINTFB(G, FB_Executive, FB_Errors)
      "Error: not a dict\n" ENDFB(G);
    return 0;
  }

  if(ok) {                      /* if session is partial, then don't error about missing stuff */
    tmp = PyDict_GetItemString(session, "partial");
    if(tmp) {
      partial_session = true;
    }
  }

  if (partial_restore) {
    G->SettingUnique->old2new = OVOneToOne_New(G->Context->heap);
  }

  if(ok) {
    tmp = PyDict_GetItemString(session, "version");
    if(tmp) {
      ok = PConvPyIntToInt(tmp, &version);
      if(ok) {
	version_full = version;
	while (version_full < 210)  /* any version less than 2.1 (account for next major version 2) should be 4 digits, otherwise 3 */
	  version_full *= 10;
        float version_float = version_full / (version >= 1000000 ? 1000000.f : 1000.f);
        if(version > _PyMOL_VERSION_int) {
          if(!quiet) {
            PRINTFB(G, FB_Executive, FB_Errors)
              "Warning: This session was created with a newer version of PyMOL (%1.6f).\n",
              version_float ENDFB(G);
            if(SettingGetGlobal_i(G, cSetting_session_version_check)) {
              PRINTFB(G, FB_Executive, FB_Errors)
                "Error: Please update first -- see http://www.pymol.org\n" ENDFB(G);
              ok = false;
            } else {
              PRINTFB(G, FB_Executive, FB_Errors)
                "Warning: Some content may not load completely.\n" ENDFB(G);
            }
          }
        } else {
          if(!quiet) {
            PRINTFB(G, FB_Executive, FB_Details)
              " Executive: Loading version %1.6f session...\n", version_float ENDFB(G);
          }
        }
      }
    }
  }
#ifndef PYMOL_EVAL
  if(ok) {
    tmp = PyDict_GetItemString(session, "eval_nag");
    if(tmp) {
      ok = PyString_Check(tmp);
      if(ok) {
        const char *st = PyString_AsString(tmp);
        if(st) {
          if(Feedback(G, FB_Nag, FB_Warnings)) {
            OrthoAddOutput(G, st);
          }
        }
      }
    }
  }
#endif

  if(ok) {                      /* colors must be restored before settings and objects */
    tmp = PyDict_GetItemString(session, "colors");
    if(tmp) {
      ok = ColorFromPyList(G, tmp, partial_restore);
    }

    if(tmp || (!partial_restore)) {     /* ignore missing if partial restore */

      if(PyErr_Occurred()) {
        PyErr_Print();
        ok = false;
      }
      if(!ok) {
        PRINTFB(G, FB_Executive, FB_Errors)
          "ExectiveSetSession-Error: after colors.\n" ENDFB(G);
      }
      if(!ok) {
        incomplete = true;
        ok = true;              /* keep trying...don't give up */
      }
    }
  }

  if(ok) {
    tmp = PyDict_GetItemString(session, "color_ext");
    if(tmp) {
      ok = ColorExtFromPyList(G, tmp, partial_restore);
    }

    if(tmp || (!partial_session)) {     /* ignore missing if partial restore */
      if(PyErr_Occurred()) {
        PyErr_Print();
        ok = false;
      }
      if(!ok) {
        PRINTFB(G, FB_Executive, FB_Errors)
          "ExectiveSetSession-Error: after color_ext.\n" ENDFB(G);
      }
      if(!ok) {
        incomplete = true;
        ok = true;              /* keep trying...don't give up */
      }
    }
  }
  if(ok) {
    tmp = PyDict_GetItemString(session, "settings");
    if(tmp && !partial_restore) {
      SettingInitGlobal(G, false, false, /* use_defaults */ true);
      SettingSetGlobalsFromPyList(G, tmp);
    }

    if(tmp || (!(partial_restore | partial_session))) { /* ignore missing if partial restore */
      if(PyErr_Occurred()) {
        PyErr_Print();
        ok = false;
      }
      if(!ok) {
        PRINTFB(G, FB_Executive, FB_Errors)
          "ExectiveSetSession-Error: after settings.\n" ENDFB(G);
      }
      if(!ok) {
        incomplete = true;
        ok = true;              /* keep trying...don't give up */
      }
    }
  }
  if(ok) {
    tmp = PyDict_GetItemString(session, "unique_settings");
    if(tmp) {
      ok = SettingUniqueFromPyList(G, tmp, partial_restore);
    }

    if(tmp || (!partial_session)) {     /* ignore missing if partial restore */
      if(PyErr_Occurred()) {
        PyErr_Print();
        ok = false;
      }
      if(!ok) {
        PRINTFB(G, FB_Executive, FB_Errors)
          "ExectiveSetSession-Error: after settings.\n" ENDFB(G);
      }
      if(!ok) {
        incomplete = true;
        ok = true;              /* keep trying...don't give up */
      }
    }
  }
  if(ok) {
    tmp = PyDict_GetItemString(session, "names");
    if(tmp) {
      if(ok)
        ok = ExecutiveSetNamedEntries(G, tmp, version, partial_restore, partial_session);
      if(!(partial_restore || partial_session)) {
        if(ok)
          ok = ExecutiveSetSelectionsFromPyList(G, tmp);
        if(ok)
          have_active = ExecutiveGetActiveSeleName(G, active, false, false);
      }
    }
    if(PyErr_Occurred()) {
      PyErr_Print();
      ok = false;
    }
    if(!ok) {
      PRINTFB(G, FB_Executive, FB_Errors)
        "ExectiveSetSession-Error: after names.\n" ENDFB(G);
    }
    if(!ok) {
      incomplete = true;
      ok = true;                /* keep trying...don't give up */
    }
  }
  if(ok && !(partial_restore)) {
    tmp = PyDict_GetItemString(session, "selector_secrets");
    if(tmp) {
      if(ok)
        ok = SelectorSecretsFromPyList(G, tmp);
    }

    if(tmp || (!(partial_restore | partial_session))) { /* ignore missing if partial restore */
      if(PyErr_Occurred()) {
        PyErr_Print();
        ok = false;
      }
      if(!ok) {
        PRINTFB(G, FB_Executive, FB_Errors)
          "ExectiveSetSession-Error: after selector secrets.\n" ENDFB(G);
      }
      if(!ok) {
        incomplete = true;
        ok = true;              /* keep trying...don't give up */
      }
    }
  }
  if(ok && !(partial_restore)) {
    tmp = PyDict_GetItemString(session, "view");
    if(tmp) {
      ok = PConvPyListToFloatArrayInPlace(tmp, sv, cSceneViewSize);
    }
    if(tmp || (!(partial_restore | partial_session))) { /* ignore missing if partial restore */
      if(ok)
        SceneSetView(G, sv, true, 0, 0);
      if(PyErr_Occurred()) {
        PyErr_Print();
        ok = false;
      }
      if(!ok) {
        PRINTFB(G, FB_Executive, FB_Errors)
          "ExectiveSetSession-Error: after view.\n" ENDFB(G);
      }
      if(!ok) {
        incomplete = true;
        ok = true;              /* keep trying...don't give up */
      }
    }
  }

  if (ok && !partial_restore) {
    tmp = CPythonVal_PyDict_GetItemString(G, session, "moviescenes");
    if (tmp) {
      MovieScenesFromPyList(G, tmp);
      CPythonVal_Free(tmp);
    }
  }

  if(ok && !(partial_restore)) {
    int warning;
    tmp = PyDict_GetItemString(session, "movie");
    if(tmp) {
      ok = MovieFromPyList(G, tmp, &warning);
    }
    if(tmp || (!(partial_restore | partial_session))) { /* ignore missing if partial restore */

      if(PyErr_Occurred()) {
        PyErr_Print();
        ok = false;
      }
      if(!ok) {
        PRINTFB(G, FB_Executive, FB_Errors)
          "ExectiveSetSession-Error: after movie.\n" ENDFB(G);
      }
      if(!ok) {
        incomplete = true;
        ok = true;              /* keep trying...don't give up */
      }
    }
  }

  if(ok && !(partial_restore)) {
    tmp = PyDict_GetItemString(session, "editor");
    if(tmp) {
      ok = EditorFromPyList(G, tmp);
    }
    if(tmp || (!(partial_restore | partial_session))) { /* ignore missing if partial restore */

      if(PyErr_Occurred()) {
        PyErr_Print();
        ok = false;
      }
      if(!ok) {
        PRINTFB(G, FB_Executive, FB_Errors)
          "ExectiveSetSession-Error: after editor.\n" ENDFB(G);
      }
      if(!ok) {
        incomplete = true;
        ok = true;              /* keep trying...don't give up */
      }
    }
  }

  if (ok) {
    PParse(G, "cmd.mouse(quiet=1)");
  }
  // Do not load viewport size when we have a GUI
  if(ok) {
    tmp = PyDict_GetItemString(session, "main");
    if(tmp) {
      if (!G->HaveGUI &&
          /* PYMOL-775 added suspend_updates check, but does it make sense? */
          !SettingGetGlobal_b(G, cSetting_suspend_updates) &&
          !partial_restore) {
        ok = MainFromPyList(G, tmp);
#ifndef _PYMOL_NOPY
      } else if (!quiet) {
        int viewport[2];
        PConvPyListToIntArrayInPlace(tmp, viewport, 2);
        PRINTFB(G, FB_Executive, FB_Actions)
          " Session was saved with: viewport %d, %d\n",
          viewport[0], viewport[1] ENDFB(G);
#endif
      }
    }
    if(tmp || (!(partial_restore | partial_session))) { /* ignore missing if partial restore */
      if(PyErr_Occurred()) {
        PyErr_Print();
        ok = false;
      }
      if(!ok) {
        PRINTFB(G, FB_Executive, FB_Errors)
          "ExectiveSetSession-Error: after main.\n" ENDFB(G);
      }
      if(!ok) {
        incomplete = true;
        ok = true;              /* keep trying...don't give up */
      }
    }
  }

  if(ok && migrate_sessions) {  /* migrate sessions */
    tmp = PyDict_GetItemString(session, "version");
    if(tmp) {
      ok = PConvPyIntToInt(tmp, &version);
      if(ok) {
        ExecutiveMigrateSession(G, version);
      }
    }
  }
  if(ok) {
    if(have_active)
      ExecutiveSetObjVisib(G, active, true, false);
  }

  OVOneToOne_DEL_AUTO_NULL(G->SettingUnique->old2new);

  if(incomplete) {
    PRINTFB(G, FB_Executive, FB_Warnings)
      "ExectiveSetSession-Warning: restore may be incomplete.\n" ENDFB(G);
  }
  G->ShaderMgr->Set_Reload_Bits(RELOAD_ALL_SHADERS);
  OrthoBackgroundTextureNeedsUpdate(G);
  ExecutiveInvalidateSelectionIndicatorsCGO(G);
  OrthoInvalidateDoDraw(G);
  SceneChanged(G);

  G->Color->HaveOldSessionColors = false;
  G->Color->HaveOldSessionExtColors = false;

  return (ok);
}


#define ExecScrollBarMargin DIP2PIXEL(1)
#define ExecScrollBarWidth DIP2PIXEL(13)

/**
 * @param sele object name (or single-object atom selection expression)
 * @param state object state (only for maps, ignored for molecules)
 * @return true if symmetry is defined
 */
pymol::Result<bool>
ExecutiveGetSymmetry(PyMOLGlobals * G, const char *sele, int state, float *a, float *b, float *c,
                        float *alpha, float *beta, float *gamma, char *sgroup)
{
  // try `sele` as an object name
  auto obj = ExecutiveFindObjectByName(G, sele);

  // odd feature: accept atom selections to select single object
  if (!obj) {
    auto tmpsele1 = SelectorTmp::make(G, sele);
    p_return_if_error(tmpsele1);
    obj = SelectorGetSingleObjectMolecule(G, tmpsele1->getIndex());
    if (!obj) {
      return pymol::make_error("selection must refer to exactly one object");
    }
  }

  CSymmetry const* symm = obj->getSymmetry(state);

  if(symm) {
    *a = symm->Crystal.dims()[0];
    *b = symm->Crystal.dims()[1];
    *c = symm->Crystal.dims()[2];
    *alpha = symm->Crystal.angles()[0];
    *beta = symm->Crystal.angles()[1];
    *gamma = symm->Crystal.angles()[2];
    UtilNCopy(sgroup, symm->spaceGroup(), sizeof(WordType));
    return true;
  }

  return false;
}

/**
 * Set symmetry for one or more objects.
 * @param names Object name pattern
 * @param state Object state (supports all=-1 and current=-2)
 * @return True if symmetry could be set for at least one object
 */
static bool ExecutiveSetSymmetry(PyMOLGlobals* G, pymol::zstring_view names,
    int const state, CSymmetry const& symmetry, bool const quiet = false)
{
  auto objects = ExecutiveGetObjectsFromPattern(G, names);
  bool success = false;

  for (auto* obj : objects) {
    if (!obj->setSymmetry(symmetry, state)) {
      PRINTFB(G, FB_Executive, FB_Warnings)
      " %s-Warning: Cannot set symmetry for '%s' (type %s)\n", __func__,
          obj->Name, typeid(obj).name() ENDFB(G);
      continue;
    }

    success = true;

    if (!quiet) {
      PRINTFB(G, FB_Executive, FB_Details)
      " %s-Details: Updated symmetry for '%s'\n", __func__, obj->Name ENDFB(G);
    }
  }

  return success;
}

pymol::Result<> ExecutiveSetSymmetry(PyMOLGlobals* G, const char* sele,
    int state, float a, float b, float c, float alpha, float beta, float gamma,
    const char* sgroup, int quiet)
{
  /* create a new symmetry object for copying */
  CSymmetry symmetry(G);
  symmetry.Crystal.setDims(a, b, c);
  symmetry.Crystal.setAngles(alpha, beta, gamma);
  symmetry.setSpaceGroup(sgroup);

  if (!ExecutiveSetSymmetry(G, sele, state, symmetry, quiet)) {
    return pymol::Error("no object selected");
  }

  return {};
}

pymol::Result<> ExecutiveSymmetryCopy(PyMOLGlobals* G, const char* source_name,
    const char* target_name, int source_state, int target_state, int quiet)
{

  /* Copy the symmetry info from source to target; currently maps can have
   * multiple states for symmetry, but ObjectMolecule cannot */

  auto const* source_obj = ExecutiveFindObjectByName(G, source_name);
  if (!source_obj) {
    return pymol::Error("source object not found");
  }

  auto const* source_symm = source_obj->getSymmetry(source_state);
  if (!source_symm) {
    return pymol::Error(pymol::string_format(
        "no symmetry in object '%s' state %d", source_name, source_state));
  }

  if (!ExecutiveSetSymmetry(
          G, target_name, target_state, *source_symm, quiet)) {
    return pymol::Error("target object not found");
  }

  return {};
}

/**
 * Performs a window average of coordinate states (trajectories).
 *
 * @param pbc Consider periodic boundary conditions
 */
pymol::Result<> ExecutiveSmooth(PyMOLGlobals* G, const char* selection,
    int cycles, int window, int first, int last, int ends, int quiet,
    float dist_cutoff, bool pbc)
{
  SETUP_SELE(selection, tmpsele1, sele);
  const char *name = tmpsele1->getName();

  ObjectMoleculeOpRec op;
  int range, offset;
  int end_skip = 0;
  bool loop = false;

  PRINTFD(G, FB_Executive)
    " %s: entered %s,%d,%d,%d,%d,%d\n", __func__, name, cycles, first, last, window,
    ends ENDFD;

  /* count the number of states over which to smooth */
  int max_state = ExecutiveCountStates(G, name) - 1;
  if(last < 0)
    last = max_state;
  if(first < 0)
    first = 0;
  if(last < first) {
    std::swap(last, first);
  }
  if(last > max_state)
    last = max_state;

  int const n_state = last - first + 1;

  int const backward = window / 2;
  int const forward = window / 2;

  switch (ends) {
  case 0:
    end_skip = 1;
    break;
  case 1:
    end_skip = 0;
    break;
  case 2:
    end_skip = backward;
    break;
  case 3:                    /* cyclic averaging */
    end_skip = 0;
    loop = true;
    break;
  default:
    end_skip = 0;
    break;
  }

  if(ends) {
    range = (last - first) + 1;
    offset = 0;
  } else {
    range = (last - end_skip) - (first + end_skip) + 1;
    offset = end_skip;
  }

  PRINTFD(G, FB_Executive)
      " %s: first %d last %d n_state %d backward %d forward %d range %d\n", __func__,
      first, last, n_state, backward, forward, range ENDFD;

  auto const window_abs = std::abs(window);
  if (window_abs < 2) {
    return pymol::Error("window must be at least size 2");
  }

  if (n_state < window_abs) {
    if (!quiet) {
      PRINTFB(G, FB_Executive, FB_Warnings)
        " %s: Window size (%d) larger than number of frames (%d)\n", __func__,
        window_abs, n_state ENDFB(G);
    }
    return {};
  }

  float const dist_cutoff_sq =
      dist_cutoff > 0 ? (dist_cutoff * dist_cutoff) : -1;

  /* determine storage req */
  ObjectMoleculeOpRecInit(&op);
  op.code = OMOP_CountAtoms;
  op.i1 = 0;
  ExecutiveObjMolSeleOp(G, sele, &op);
  auto const n_atom = size_t(op.i1);

  // TODO Could we call SelectorCountAtoms instead?
  assert(n_atom == SelectorCountAtoms(G, sele, cStateAll));

  if (!n_atom) {
    return {};
  }

  auto pbc_obj = pbc ? SelectorGetFirstObjectMolecule(G, sele) : nullptr;
  if (pbc_obj) {
    ObjectMoleculePBCUnwrap(*pbc_obj);
  }

  auto const n_index = n_atom * n_state;
  auto coord0 = std::vector<float>(3 * n_index);
  auto flag0 = std::vector<int>(n_index);
  auto flag1 = std::vector<int>(n_index);

  /* get the data */
  if(!quiet) {
    PRINTFB(G, FB_Executive, FB_Actions)
      " Smooth: copying coordinates to temporary arrays..\n" ENDFB(G);
  }
  op.code = OMOP_CSetIdxGetAndFlag;
  op.i1 = n_atom;
  op.i2 = 0;
  op.cs1 = first;
  op.cs2 = last;
  op.vv1 = coord0.data();
  op.ii1 = flag0.data();
  op.nvv1 = 0;
  ExecutiveObjMolSeleOp(G, sele, &op);

  PRINTFD(G, FB_Executive)
          " %s: got %d %d\n", __func__, op.i2, op.nvv1 ENDFD;

  // make a copy of the input coordinates
  auto coord1 = coord0;

  for (int a = 0; a < cycles; ++a) {
    if(!quiet) {
      PRINTFB(G, FB_Executive, FB_Actions)
        " Smooth: smoothing (pass %d)...\n", a + 1 ENDFB(G);
    }

    // use coordiantes from previous pass as input
    std::swap(coord0, coord1);

    // views for friendly indexing
    auto const* coord0x3 = pymol::reshape<3>(coord0.data());
    auto* const coord1x3 = pymol::reshape<3>(coord1.data());

#ifdef PYMOL_OPENMP
#pragma omp parallel for
#endif
    for (int b = 0; b < range; ++b) {
      int const st_b = b + offset;

      if (st_b < end_skip || st_b >= (n_state - end_skip)) {
        continue;
      }

      for (int c = 0; c < n_atom; ++c) {
        auto const index_c = n_atom * st_b + c;

        if (!flag0[index_c]) {
          continue;
        }

        float sum[3] = {};
        int cnt = 0;
        int st_prev = 0;
        for (int d = -backward; d <= forward; ++d) {
          int st = st_b + d;
          if(loop) {
            if(st < 0) {
              st = n_state + st;
            } else if(st >= n_state) {
              st = st - n_state;
            }
          } else {
            if(st < 0) {
              st = 0;
            } else if(st >= n_state) {
              st = n_state - 1;
            }
          }

          auto const index = n_atom * st + c;
          if (!flag0[index]) {
            continue;
          }

          float const* v0 = coord0x3[index];

          if (dist_cutoff_sq > 0 && cnt) {
            auto const* v0_prev = coord0x3[n_atom * st_prev + c];
            if (diffsq3f(v0, v0_prev) > dist_cutoff_sq) {
              if (d <= 0) {
                scale3f(v0, cnt, sum);
              } else {
                for (;d <= forward; ++d) {
                  add3f(sum, v0_prev, sum);
                  ++cnt;
                }
                break;
              }
            }
          }

          /* atom's avg position */
          add3f(sum, v0, sum);
          ++cnt;
          st_prev = st;
        }
        if(cnt) {
          flag1[index_c] = true;
          scale3f(sum, 1.0F / cnt, coord1x3[index_c]);
        }
      }
    }
  }

  if(!quiet) {
    PRINTFB(G, FB_Executive, FB_Actions)
      " Smooth: updating coordinates...\n" ENDFB(G);
  }

  /* set the new coordinates */

  op.code = OMOP_CSetIdxSetFlagged;
  op.i1 = n_atom;
  op.i2 = 0;
  if(ends) {
    op.cs1 = first;
    op.cs2 = last;
    op.vv1 = coord1.data();
    op.ii1 = flag1.data();
  } else {
    op.cs1 = first + end_skip;
    op.cs2 = last - end_skip;
    op.vv1 = coord1.data() + (end_skip * n_atom * 3);
    op.ii1 = flag1.data() + (end_skip * n_atom);
  }
  op.nvv1 = 0;

  ExecutiveObjMolSeleOp(G, sele, &op);
  PRINTFD(G, FB_Executive)
          " %s: put %d %d\n", __func__, op.i2, op.nvv1 ENDFD;

  if (pbc_obj) {
    ObjectMoleculePBCWrap(*pbc_obj);
  }

  return {};
}


/*========================================================================*/
int ExecutiveDebug(PyMOLGlobals * G, const char *name)
{
  ObjectMolecule *obj;
  ObjectMoleculeBPRec bp;
  int a;

  obj = (ObjectMolecule *) ExecutiveFindObjectByName(G, name);
  if(obj) {
    ObjectMoleculeInitBondPath(obj, &bp);
    ObjectMoleculeGetBondPaths(obj, 0, 10, &bp);
    for(a = 0; a < bp.n_atom; a++) {
      printf("%d %d %d\n", a, bp.list[a], bp.dist[bp.list[a]]);
    }

    ObjectMoleculePurgeBondPath(obj, &bp);
  }
  return (1);
}


/*========================================================================*/
int ***ExecutiveGetBondPrint(PyMOLGlobals * G, const char *name, int max_bond, int max_type,
                             int *dim)
{
  int ***result = NULL;
  pymol::CObject *obj;
  ObjectMolecule *objMol;

  obj = ExecutiveFindObjectByName(G, name);
  if(obj->type == cObjectMolecule) {
    objMol = (ObjectMolecule *) obj;
    result = ObjectMoleculeGetBondPrint(objMol, max_bond, max_type, dim);
  }
  return (result);
}


/*========================================================================*/
#define cMapOperatorMinimum 0
#define cMapOperatorMaximum 1
#define cMapOperatorSum     2
#define cMapOperatorAverage 3
#define cMapOperatorDifference 4
#define cMapOperatorCopy     5
#define cMapOperatorUnique   6

pymol::Result<> ExecutiveMapSet(PyMOLGlobals* G, const char* name,
    int operator_, const char* operands, int target_state, int source_state,
    int zoom, int quiet)
{
  CExecutive *I = G->Executive;
  CTracker *I_Tracker = I->Tracker;
  int isNew = false;
  ObjectMap *target = ExecutiveFindObjectMapByName(G, name);
  ObjectMap *first_operand = NULL;
  int src_state_start = 0, src_state_stop = 0;
  int list_id = ExecutiveGetNamesListFromPattern(G, operands, true, true);

  if(target_state < 0)          /* if we're targeting all states, 0 is the offset */
    target_state = 0;

  /* first, figure out what the range of input states is */

  if(source_state < 0) {        /* all source states */
    int iter_id = TrackerNewIter(I_Tracker, 0, list_id);
    int max_n_state = 0;
    SpecRec *rec;
    while(TrackerIterNextCandInList(I_Tracker, iter_id, (TrackerRef **) (void *) &rec)) {
      if(rec) {
        switch (rec->type) {
        case cExecObject:
          if(rec->obj->type == cObjectMap) {
            ObjectMap *obj = (ObjectMap *) rec->obj;
            if(obj->State.size() > max_n_state)
              max_n_state = obj->State.size();        /* count states */
          }
        }
      }
    }
    TrackerDelIter(I_Tracker, iter_id);
    src_state_start = 0;
    src_state_stop = max_n_state;
  } else {
    src_state_start = source_state;
    src_state_stop = source_state + 1;
  }

  {
    /* next, find the first operand */

    OrthoLineType first_op_st;
    ParseWordCopy(first_op_st, operands, sizeof(OrthoLineType) - 1);    /* copy the first operand */
    {
      int sub_list_id = ExecutiveGetNamesListFromPattern(G, first_op_st, true, true);
      int sub_iter_id = TrackerNewIter(I_Tracker, 0, sub_list_id);
      SpecRec *rec;

      while(TrackerIterNextCandInList(I_Tracker, sub_iter_id,
                                      (TrackerRef **) (void *) &rec)) {
        if(rec) {
          switch (rec->type) {
          case cExecObject:
            if(rec->obj->type == cObjectMap) {
              ObjectMap *obj = (ObjectMap *) rec->obj;
              first_operand = obj;
            }
            break;
          }
        }
        if(first_operand)
          break;
      }
      TrackerDelList(I_Tracker, sub_list_id);
      TrackerDelIter(I_Tracker, sub_iter_id);
    }
  }

  {

    /* okay, next thing we need to worry about is where we're putting the data.

       Case 1. If the map already exists, then we'll use the existing map points for storing
       the result.

       Case 2. If the operation implies a copy of existing map geometry, then we'll create that
       copy first before performing the calulation.

       Case 3. If the operation implies a new map geometry, then we need to compute that geometry and
       create the map.

     */

    if(!target) {               /* target map doesn't exist... */
      int need_union_geometry = false;
      int need_first_geometry = false;
      switch (operator_) {
      case cMapOperatorSum:
      case cMapOperatorAverage:
      case cMapOperatorMinimum:
      case cMapOperatorMaximum:
      case cMapOperatorDifference:
        need_union_geometry = true;
        break;
      case cMapOperatorUnique:
      case cMapOperatorCopy:
        need_first_geometry = true;
        break;
      }

      if(need_union_geometry) {
        int src_state, trg_state;
        ObjectMapDesc desc;
        target = new ObjectMap(G);

        ObjectSetName(target, name);
        isNew = true;

        for(src_state = src_state_start; src_state < src_state_stop; src_state++) {
          trg_state = src_state + target_state;
          desc.mode = cObjectMap_OrthoMinMaxGrid;       /* Orthorhombic: min, max, 
                                                           spacing, 
                                                           centered over range  */
          desc.init_mode = 0;   /* zeros */
          desc.Grid[0] = 1.0F;
          desc.Grid[1] = 1.0F;
          desc.Grid[2] = 1.0F;

          {
            int iter_id = TrackerNewIter(I_Tracker, 0, list_id);
            float grid_sum[3] = { 0.0F, 0.0F, 0.0F };
            int grid_count = 0;
            int first_extent = true;

            SpecRec *rec;
            while(TrackerIterNextCandInList(I_Tracker, iter_id,
                                            (TrackerRef **) (void *) &rec)) {
              if(rec) {
                /* compute an average grid and get the effective "union" extent */
                switch (rec->type) {
                case cExecObject:
                  if(rec->obj->type == cObjectMap) {

                    ObjectMap *obj = (ObjectMap *) rec->obj;
                    ObjectMapState *ms = &obj->State[src_state];
                    if(src_state < obj->State.size()) {
                      if(ms->Active) {
                        if(first_extent) {
                          copy3f(ms->ExtentMin, desc.MinCorner);
                          copy3f(ms->ExtentMax, desc.MaxCorner);
                          first_extent = false;
                        } else {
                          int b;
                          for(b = 0; b < 3; b++) {
                            if(ms->ExtentMin[b] < desc.MinCorner[b])
                              desc.MinCorner[b] = ms->ExtentMin[b];
                            if(ms->ExtentMax[b] > desc.MaxCorner[b])
                              desc.MaxCorner[b] = ms->ExtentMax[b];
                          }
                        }
                        if(!ObjectMapStateValidXtal(ms)) {
                          /* other general-purpose maps currently handled */
                          int b;
                          for(b = 0; b < 3; b++) {
                            grid_sum[b] += ms->Grid[b];
                          }
                          grid_count++;
                        }
                      }
                    }
                  }
                }
              }
            }
            TrackerDelIter(I_Tracker, iter_id);
            if(grid_count) {
              int b;
              for(b = 0; b < 3; b++) {
                desc.Grid[b] = grid_sum[b] / grid_count;
              }
            }
            if(!first_extent) {
              add3f(desc.Grid, desc.MaxCorner, desc.MaxCorner);
              subtract3f(desc.MinCorner, desc.Grid, desc.MinCorner);
              ObjectMapNewStateFromDesc(G, target, &desc, trg_state, quiet);
              target->State[trg_state].Active = true;
            }
          }
        }
        /* need union geometry */
      } else if(need_first_geometry) {
        if(first_operand) {
          if(ObjectMapNewCopy(G, first_operand, &target, source_state, target_state)) {
            if(target) {
              ObjectSetName(target, name);
              isNew = true;
            }
          }
        }
      }
    }
  }

  if(!target) {
    TrackerDelList(I_Tracker, list_id);
    return pymol::make_error("Cannot find or construct target map.");
  }

  /* now do the actual operation */

  int src_state;
  for(src_state = src_state_start; src_state < src_state_stop; src_state++) {
    int trg_state = src_state + target_state;
    ObjectMapState *ms;
    VecCheckEmplace(target->State, trg_state, G);

    ms = &target->State[target_state];
    if(ms->Active) {
      int iter_id = TrackerNewIter(I_Tracker, 0, list_id);
      int n_pnt = (ms->Field->points->size() / ms->Field->points->base_size) / 3;
      float *pnt = (float *) ms->Field->points->data.data();
      float *r_value = pymol::malloc<float>(n_pnt);
      float *l_value = pymol::calloc<float>(n_pnt);
      int *present = pymol::calloc<int>(n_pnt);
      int *inside = pymol::malloc<int>(n_pnt);
      SpecRec *rec;

      while(TrackerIterNextCandInList(I_Tracker, iter_id,
                                      (TrackerRef **) (void *) &rec)) {
        if(rec) {
          if(rec->type == cExecObject) {
            if(rec->obj->type == cObjectMap) {
              ObjectMap *obj = (ObjectMap *) rec->obj;
              if (ObjectMapInterpolate(obj, src_state, pnt, r_value, inside, n_pnt))
              {
                int a;
                float *rv = r_value;
                float *lv = l_value;
                int *flg = inside;
                int *pre = present;

                switch (operator_) {
                case cMapOperatorCopy:
                  for(a = 0; a < n_pnt; a++) {
                    if(flg) {
                      *lv = *rv;
                    }
                    rv++;
                    lv++;
                    flg++;
                  }
                  break;
                case cMapOperatorMinimum:
                  for(a = 0; a < n_pnt; a++) {
                    if(flg) {
                      if(*pre) {
                        if(*lv > *rv)
                          *lv = *rv;
                      } else {        /* first map */
                        *pre = 1;
                        *lv = *rv;
                      }

                    }
                    rv++;
                    lv++;
                    flg++;
                    pre++;
                  }
                  break;
                case cMapOperatorMaximum:
                  for(a = 0; a < n_pnt; a++) {
                    if(flg) {
                      if(*pre) {
                        if(*lv < *rv)
                          *lv = *rv;
                      } else {        /* first map */
                        *pre = 1;
                        *lv = *rv;
                      }
                    }
                    rv++;
                    lv++;
                    flg++;
                    pre++;
                  }
                  break;
                case cMapOperatorSum:
                  for(a = 0; a < n_pnt; a++) {
                    if(flg) {
                      *lv += *rv;
                    }
                    rv++;
                    lv++;
                    flg++;
                  }
                  break;
                case cMapOperatorAverage:
                  for(a = 0; a < n_pnt; a++) {
                    if(flg) {
                      *lv += *rv;
                    }
                    (*pre)++;
                    rv++;
                    lv++;
                    flg++;
                    pre++;
                  }
                  break;
                case cMapOperatorDifference:
                  if(obj != first_operand) {
                    for(a = 0; a < n_pnt; a++) {
                      if(flg) {
                        *lv -= *rv;
                      }
                      rv++;
                      lv++;
                      flg++;
                    }
                  } else {
                    for(a = 0; a < n_pnt; a++) {
                      if(flg) {
                        *lv += *rv;
                      }
                      rv++;
                      lv++;
                      flg++;
                    }
                  }
                  break;
                case cMapOperatorUnique:
                  if(obj != first_operand) {
                    for(a = 0; a < n_pnt; a++) {
                      if(flg) {
                        *lv -= *rv;
                      }
                      rv++;
                      lv++;
                      flg++;
                    }
                  } else {
                    for(a = 0; a < n_pnt; a++) {
                      if(flg) {
                        *lv += *rv;
                      }
                      rv++;
                      lv++;
                      flg++;
                    }
                  }

                  break;
                }
              }
            }
          }
        }
      }

      {
        int a;
        float *lv = l_value;
        int *pre = present;

        switch (operator_) {
        case cMapOperatorUnique:
          lv = l_value;
          for(a = 0; a < n_pnt; a++) {
            if(*lv < 0.0F)
              *lv = 0.0F;
            lv++;
          }
          break;
        case cMapOperatorAverage:
          lv = l_value;
          pre = present;
          for(a = 0; a < n_pnt; a++) {
            if(*pre)
              *lv /= *pre;
            lv++;
            pre++;
          }
        }
      }

      /* copy after calculation so that operand can include target */

      memcpy(ms->Field->data->data.data(), l_value, n_pnt * sizeof(float));

      FreeP(present);
      FreeP(l_value);
      FreeP(r_value);
      FreeP(inside);
      TrackerDelIter(I_Tracker, iter_id);
    }
  }
  

  /* and finally, update */

  if(target) {
    ObjectMapUpdateExtents(target);
    if(isNew) {
      ExecutiveManageObject(G, target, -1, quiet);
    } else {
      ExecutiveDoZoom(G, target, false, zoom, true);
    }
    SceneChanged(G);
  }
  TrackerDelList(I_Tracker, list_id);

  return {};
}


const char * ExecutiveMapGenerate(PyMOLGlobals * G, const char * name, const char * reflection_file, const char * tempFile,
				  char * amplitudes, const char * phases, const char * weights, double reso_low,
				  double reso_high, const char * space_group, double cell[6], int quiet, int zoom)
{
  /* FIXME: this should be returned in memory!! */
  /* In the meantime, use mkstemp and load in Python */
  int ok;

  ok = 0;

  if (weights && (!strncmp(weights,"None",4))) 
    weights=NULL;

  /* printf("Passing to primex driver: space_group=%s, cell=[%f %f %f %f %f %f], reso_high=%f, rseo_low=%f, refl_file=%s, ampl=%s, phases=%s, weights=%s, map_file=%s", space_group, cell[0], cell[1], cell[2], cell[3], cell[4], cell[5], reso_high, reso_low, reflection_file, amplitudes, phases, weights, tempFile); */

#ifndef NO_MMLIBS
  ok = !(primex_pymol_driver2(space_group, cell, reso_high, reso_low, reflection_file, amplitudes,
			    phases, weights, tempFile)); 
#else
  PRINTFB(G, FB_Executive, FB_Errors)
    " Error: MTZ map loading not supported in this PyMOL build.\n" ENDFB(G);
#endif

  if (!ok) 
    return NULL;
  else
    return (const char*) tempFile;

}

pymol::Result<>
ExecutiveMapNew(PyMOLGlobals * G, const char *name, int type, float grid_spacing,
                const char *s0, float buffer,
                const float *minCorner,
                const float *maxCorner, int state, int have_corners,
                int quiet, int zoom, int normalize, float clamp_floor,
                float clamp_ceiling, float resolution)
{
  pymol::CObject *origObj = NULL;
  ObjectMap *objMap;
  ObjectMapState *ms = NULL;
  int a;
  float v[3];
  ObjectMapDesc _md, *md;

  SETUP_SELE_DEFAULT(0);
  const char* sele = tmpsele0->getName();

  int isNew = true;
  int n_state;
  int valid_extent = false;
  int st;
  int st_once_flag = true;
  int n_st;
  int extent_state;
  int clamp_flag = (clamp_floor <= clamp_ceiling);

  md = &_md;

  float grid[3] = {grid_spacing, grid_spacing, grid_spacing};

  if((state == -2) || (state == -3))    /* TO DO: support per-object states */
    state = SceneGetState(G);

  /* remove object if it already exists */

  origObj = ExecutiveFindObjectByName(G, name);

  if(origObj) {
    if(origObj->type != cObjectMap) {
      ExecutiveDelete(G, origObj->Name);
    } else {
      isNew = false;
    }
  }

  n_st = ExecutiveCountStates(G, NULL);

  for(st = 0; st < n_st; st++) {
    if(state == -1)
      st_once_flag = false;     /* each state, separate map, separate extent */
    if(!st_once_flag)
      state = st;
    extent_state = state;
    if(state <= -3)
      extent_state = -1;
    if(strlen(sele) && (!have_corners)) {
      valid_extent = ExecutiveGetExtent(G, sele, md->MinCorner,
                                        md->MaxCorner, true, extent_state, false);
      /* TODO restrict to state */
    } else {
      valid_extent = 1;
      copy3f(minCorner, md->MinCorner);
      copy3f(maxCorner, md->MaxCorner);
    }
    copy3f(grid, md->Grid);

    subtract3f(md->MaxCorner, md->MinCorner, v);
    for(a = 0; a < 3; a++) {
      if(v[a] < 0.0)
        std::swap(md->MaxCorner[a], md->MinCorner[a]);
    };
    subtract3f(md->MaxCorner, md->MinCorner, v);

    if(buffer < -0.5F) { // buffer == -1
      buffer = SettingGetGlobal_f(G, cSetting_gaussian_resolution);
    }

    if(buffer > 0.0F) {
      for(a = 0; a < 3; a++) {
        md->MinCorner[a] -= buffer;
        md->MaxCorner[a] += buffer;
      }
    }
    md->mode = cObjectMap_OrthoMinMaxGrid;
    md->init_mode = -1;         /* no initialization */

    /* validate grid */
    for(a = 0; a < 3; a++)
      if(md->Grid[a] <= R_SMALL8)
        md->Grid[a] = R_SMALL8;

    if(isNew)
      objMap = new ObjectMap(G);
    else
      objMap = (ObjectMap *) origObj;
    if(objMap) {
      int once_flag = true;
      n_state = SelectorCountStates(G, sele0);
      if(valid_extent)
        for(a = 0; a < n_state; a++) {
          if(state == -5)
            once_flag = false;        /* api: state=-4 = each state, separate map, shared extent */
          if(state == -4)
            state = -1;       /* api: state=-3 all states, but one map */
          if(!once_flag)
            state = a;
          ms = ObjectMapNewStateFromDesc(G, objMap, md, state, quiet);
          if(!ms) {
            return pymol::make_error("Invalid state ", state);
          }

          switch (type) {
          case 0:          /* vdw */
            SelectorMapMaskVDW(G, sele0, ms, 0.0F, state);
            break;
          case 1:          /* coulomb */
            SelectorMapCoulomb(G, sele0, ms, 0.0F, state, false, false, 1.0F);
            break;
          case 2:          /* gaussian */
            SelectorMapGaussian(G, sele0, ms, 0.0F, state, normalize, false, quiet, resolution);
            break;
          case 3:          /* coulomb_neutral */
            SelectorMapCoulomb(G, sele0, ms, 0.0F, state, true, false, 1.0F);
            break;
          case 4:          /* coulomb_local */
            SelectorMapCoulomb(G, sele0, ms,
                               SettingGetGlobal_f(G, cSetting_coulomb_cutoff), state,
                               false, true, 2.0F);
            break;
          case 5:          /* gaussian_max */
            SelectorMapGaussian(G, sele0, ms, 0.0F, state, normalize, true, quiet, resolution);
            break;
          }
          if(!ms->Active)
            ObjectMapStatePurge(G, ms);
          else if(clamp_flag) {
            ObjectMapStateClamp(ms, clamp_floor, clamp_ceiling);
          }
          
          if(once_flag)
            break;
        }

      ObjectSetName(objMap, name);
      ObjectMapUpdateExtents(objMap);
      if(isNew) {
        ExecutiveManageObject(G, objMap, -1, quiet);
      } else {
        ExecutiveDoZoom(G, objMap, false, zoom, true);
      }
      isNew = false;
      origObj = objMap;
    }
    SceneChanged(G);

    if(st_once_flag)
      break;
  }
  return {};
}


/*========================================================================*/
int ExecutiveSculptIterateAll(PyMOLGlobals * G)
{
  int active = false;
  float center_array[8] = { 0.0F, 0.0F, 0.0F, 0.0F };
  float *center = center_array;
  CExecutive *I = G->Executive;
  SpecRec *rec = NULL;
  ObjectMolecule *objMol;
  CGOReset(G->DebugCGO);

  if(SettingGetGlobal_b(G, cSetting_sculpting)) {
    if(!SettingGetGlobal_b(G, cSetting_sculpt_auto_center))
      center = NULL;

    while(ListIterate(I->Spec, rec, next)) {
      if(rec->type == cExecObject) {
        if(rec->obj->type == cObjectMolecule) {
          objMol = (ObjectMolecule *) rec->obj;
          if(SettingGet_b(G, NULL, objMol->Setting.get(), cSetting_sculpting)) {
            constexpr int state = -2; // current state
            ObjectMoleculeSculptIterate(objMol, state,
                                        SettingGet_i(G, NULL, objMol->Setting.get(),
                                                     cSetting_sculpting_cycles), center);
            active = true;
          }
        }
      }
    }
    if(center && (center[3] > 1.0F)) {
      float pos[3];
      SceneGetCenter(G, pos);
      center[3] = 1.0F / center[3];
      scale3f(center, center[3], center);
      center[7] = 1.0F / center[7];
      scale3f(center + 4, center[7], center + 4);
      subtract3f(center, center + 4, center);
      add3f(pos, center, center);
      ExecutiveCenter(G, NULL, -1, true, false, center, true);
    }
  }
  if (active){
    EditorInvalidateShaderCGO(G);
  }
  return (active);
}


/*========================================================================*/
float ExecutiveSculptIterate(PyMOLGlobals * G, const char *name, int state, int n_cycle)
{
  pymol::CObject *obj = ExecutiveFindObjectByName(G, name);
  CExecutive *I = G->Executive;
  SpecRec *rec = NULL;
  ObjectMolecule *objMol;
  float total_strain = 0.0F;

  if(WordMatchExact(G, name, cKeywordAll, true)) {
    while(ListIterate(I->Spec, rec, next)) {
      if(rec->type == cExecObject) {
        if(rec->obj->type == cObjectMolecule) {
          objMol = (ObjectMolecule *) rec->obj;
          total_strain += ObjectMoleculeSculptIterate(objMol, state, n_cycle, NULL);
        }
      }
    }
  } else if(!obj) {
    PRINTFB(G, FB_Executive, FB_Errors)
      "Executive-Error: object %s not found.\n", name ENDFB(G);
  } else if(obj->type != cObjectMolecule) {
    PRINTFB(G, FB_Executive, FB_Errors)
      "Executive-Error: object %s is not a molecular object.\n", name ENDFB(G);
  } else {
    total_strain =
      ObjectMoleculeSculptIterate((ObjectMolecule *) obj, state, n_cycle, NULL);
  }
  return (total_strain);
}


/*========================================================================*/
int ExecutiveSculptActivate(PyMOLGlobals * G, const char *name, int state, int match_state,
                            int match_by_segment)
{
  pymol::CObject *obj = ExecutiveFindObjectByName(G, name);
  SpecRec *rec = NULL;
  ObjectMolecule *objMol;
  CExecutive *I = G->Executive;
  int ok = true;
  if(state < 0)
    state = SceneGetState(G);

  if(WordMatchExact(G, name, cKeywordAll, true)) {
    while(ListIterate(I->Spec, rec, next)) {
      if(rec->type == cExecObject) {
        if(rec->obj->type == cObjectMolecule) {
          objMol = (ObjectMolecule *) rec->obj;
          ObjectMoleculeSculptImprint(objMol, state, match_state, match_by_segment);
        }
      }
    }
  } else if(!obj) {
    PRINTFB(G, FB_Executive, FB_Errors)
      "Executive-Error: object %s not found.\n", name ENDFB(G);
    ok = false;
  } else if(obj->type != cObjectMolecule) {
    PRINTFB(G, FB_Executive, FB_Errors)
      "Executive-Error: object %s is not a molecular object.\n", name ENDFB(G);
    ok = false;
  } else {
    ObjectMoleculeSculptImprint((ObjectMolecule *) obj, state, match_state,
                                match_by_segment);
  }
  return (ok);
}


/*========================================================================*/
int ExecutiveSculptDeactivate(PyMOLGlobals * G, const char *name)
{
  pymol::CObject *obj = ExecutiveFindObjectByName(G, name);
  SpecRec *rec = NULL;
  ObjectMolecule *objMol;
  CExecutive *I = G->Executive;

  int ok = true;

  if(WordMatchExact(G, name, cKeywordAll, true)) {
    while(ListIterate(I->Spec, rec, next)) {
      if(rec->type == cExecObject) {
        if(rec->obj->type == cObjectMolecule) {
          objMol = (ObjectMolecule *) rec->obj;
          ObjectMoleculeSculptClear(objMol);
        }
      }
    }
  } else if(!obj) {
    PRINTFB(G, FB_Executive, FB_Errors)
      "Executive-Error: object %s not found.\n", name ENDFB(G);
    ok = false;
  } else if(obj->type != cObjectMolecule) {
    PRINTFB(G, FB_Executive, FB_Errors)
      "Executive-Error: object %s is not a molecular object.\n", name ENDFB(G);
    ok = false;
  } else {
    ObjectMoleculeSculptClear((ObjectMolecule *) obj);
  }
  return (ok);
}


/*========================================================================*/
pymol::Result<> ExecutiveSetGeometry(
    PyMOLGlobals* G, const char* s1, int geom, int valence)
{
  SETUP_SELE_DEFAULT(1);

  int count = 0;
  auto recs = pymol::make_list_adapter(G->Executive->Spec);
  for (auto& rec : recs) {
    auto obj = ExecutiveIsObjectType(rec, cObjectMolecule) ? rec.obj : nullptr;
    if(obj) {
      count += ObjectMoleculeSetGeometry(
        G, static_cast<ObjectMolecule*>(obj), sele1, geom, valence);
    }
  }
  if (count == 0) {
    return pymol::make_error("Empty selection.");
  }
  return {};
}


/*========================================================================*/
int ExecutiveMapSetBorder(PyMOLGlobals * G, const char *name, float level, int state)
{
  CExecutive *I = G->Executive;
  int result = true;
  CTracker *I_Tracker = I->Tracker;
  int list_id = ExecutiveGetNamesListFromPattern(G, name, true, true);
  int iter_id = TrackerNewIter(I_Tracker, 0, list_id);
  SpecRec *rec;

  while(TrackerIterNextCandInList(I_Tracker, iter_id, (TrackerRef **) (void *) &rec)) {
    if(rec) {
      switch (rec->type) {
      case cExecObject:
        if(rec->obj->type == cObjectMap) {
          ObjectMap *obj = (ObjectMap *) rec->obj;
          result = ObjectMapSetBorder(obj, level, state);

          if(result) {
            ExecutiveInvalidateMapDependents(G, obj->Name);
          }
        }
        break;
      }
    }
  }

  TrackerDelList(I_Tracker, list_id);
  TrackerDelIter(I_Tracker, iter_id);
  return result;
}

pymol::Result<> ExecutiveMapDouble(PyMOLGlobals* G, const char* name, int state)
{
  CExecutive *I = G->Executive;
  CTracker *I_Tracker = I->Tracker;
  int list_id = ExecutiveGetNamesListFromPattern(G, name, true, true);
  int iter_id = TrackerNewIter(I_Tracker, 0, list_id);
  SpecRec *rec;

  while(TrackerIterNextCandInList(I_Tracker, iter_id, (TrackerRef **) (void *) &rec)) {
    if(rec) {
      switch (rec->type) {
      case cExecObject:
        if(rec->obj->type == cObjectMap) {
          ObjectMap *obj = (ObjectMap *) rec->obj;
          auto result = ObjectMapDouble(obj, state);
          if(result) {
            ExecutiveInvalidateMapDependents(G, obj->Name);
          } else {
            return result;
          }
          if(result && rec->visible)
            SceneChanged(G);
        }
        break;
      }
    }
  }

  TrackerDelList(I_Tracker, list_id);
  TrackerDelIter(I_Tracker, iter_id);
  return {};
}

pymol::Result<> ExecutiveMapHalve(
    PyMOLGlobals* G, const char* name, int state, int smooth)
{
  CExecutive *I = G->Executive;
  CTracker *I_Tracker = I->Tracker;
  int list_id = ExecutiveGetNamesListFromPattern(G, name, true, true);
  int iter_id = TrackerNewIter(I_Tracker, 0, list_id);
  SpecRec *rec;

  while(TrackerIterNextCandInList(I_Tracker, iter_id, (TrackerRef **) (void *) &rec)) {
    if(rec) {
      switch (rec->type) {
      case cExecObject:
        if(rec->obj->type == cObjectMap) {
          ObjectMap *obj = (ObjectMap *) rec->obj;
          auto result = ObjectMapHalve(obj, state, smooth);
          if(result) {
            ExecutiveInvalidateMapDependents(G, obj->Name);
          } else {
            return result;
          }
          if(result && rec->visible)
            SceneChanged(G);
        }
        break;
      }
    }
  }

  TrackerDelList(I_Tracker, list_id);
  TrackerDelIter(I_Tracker, iter_id);
  return {};
}

pymol::Result<> ExecutiveMapTrim(PyMOLGlobals* G, const char* name,
    const char* sele, float buffer, int map_state, int sele_state, int quiet)
{
  auto s1 = SelectorTmp2::make(G, sele);
  CExecutive *I = G->Executive;
  float mn[3], mx[3];
  if(ExecutiveGetExtent(G, s1->getName(), mn, mx, true, sele_state, false)) {
    CTracker *I_Tracker = I->Tracker;
    int list_id = ExecutiveGetNamesListFromPattern(G, name, true, true);
    int iter_id = TrackerNewIter(I_Tracker, 0, list_id);
    SpecRec *rec;

    {
      int a;
      float t;
      for(a = 0; a < 3; a++) {
        mn[a] -= buffer;
        mx[a] += buffer;
        if(mn[a] > mx[a]) {
          t = mn[a];
          mn[a] = mx[a];
          mx[a] = t;
        }
      }
    }
    while(TrackerIterNextCandInList(I_Tracker, iter_id, (TrackerRef **) (void *) &rec)) {
      if(rec) {
        switch (rec->type) {
        case cExecObject:
          if(rec->obj->type == cObjectMap) {
            ObjectMap *obj = (ObjectMap *) rec->obj;
            auto result = ObjectMapTrim(obj, map_state, mn, mx, quiet);
            if(result) {
              ExecutiveInvalidateMapDependents(G, obj->Name);
            } else {
              return result;
            }
            if(result && rec->visible)
              SceneChanged(G);
          }
          break;
        }
      }
    }
    TrackerDelList(I_Tracker, list_id);
    TrackerDelIter(I_Tracker, iter_id);
  }
  return {};
}

pymol::Result<NetSelect> ExecutiveSelectRect(PyMOLGlobals * G, BlockRect * rect, int mode)
{
  NetSelect result;
  pymol::Result<int> selectorResult;
  int beforeSelected = 0;
  Multipick smp;
  char selName[WordLength] = cLeftButSele;
  char prefix[3] = "";
  int log_box = 0;
  const char *sel_mode_kw = "";

  int logging = SettingGet<int>(G, cSetting_logging);
  if(logging)
    log_box = SettingGetGlobal_b(G, cSetting_log_box_selections);
  /*  if(logging==cPLog_pml)
     strcpy(prefix,"_ "); */
  smp.x = rect->left;
  smp.y = rect->bottom;
  smp.w = rect->right - rect->left;
  smp.h = rect->top - rect->bottom;
  SceneMultipick(G, &smp);
  if (!smp.picked.empty()) {
    selectorResult = SelectorCreate(G, cTempRectSele, nullptr, nullptr, 1, &smp);
    if(log_box)
      SelectorLogSele(G, cTempRectSele);
    switch (mode) {
    case cButModeRect:
      if(mode == cButModeRect) {
        selectorResult = SelectorCreate(G, cLeftButSele, cTempRectSele, nullptr, 1, nullptr);
        if(log_box) {
          auto buf2 = pymol::string_format("%scmd.select(\"%s\",\"%s\",enable=1)\n", prefix, cLeftButSele,
                  cTempRectSele);
          PLog(G, buf2, cPLog_no_flush);
        }
      }
      break;
    case cButModeSeleSetBox:
    case cButModeSeleAddBox:
    case cButModeSeleSubBox:
      ExecutiveGetActiveSeleName(G, selName, true, SettingGetGlobal_i(G, cSetting_logging));
      sel_mode_kw = SceneGetSeleModeKeyword(G);
      /* intentional omission of break! */
    case cButModeRectAdd:
    case cButModeRectSub:
      beforeSelected = SelectorCountAtoms(G, SelectorTmp{G, selName}.getIndex(), -1);
      if(SelectorIndexByName(G, selName) >= 0) {
        if((mode == cButModeRectAdd) || (mode == cButModeSeleAddBox)) {
          auto buffer = pymol::string_format("(?%s or %s(%s))", selName, sel_mode_kw, cTempRectSele);
          selectorResult = SelectorCreate(G, selName, buffer.c_str(), nullptr, 0, nullptr);
          if(log_box) {
            auto buf2 = pymol::string_format("%scmd.select(\"%s\",\"(%s)\",enable=1)\n", prefix, selName,
                    buffer);
            PLog(G, buf2, cPLog_no_flush);
          }
        } else if((mode == cButModeRectSub) || (mode == cButModeSeleSubBox)) {
          auto buffer = pymol::string_format("(%s(?%s) and not %s(%s))", sel_mode_kw, selName, sel_mode_kw,
                  cTempRectSele);
          selectorResult = SelectorCreate(G, selName, buffer.c_str(), nullptr, 0, nullptr);
          if(log_box) {
            auto buf2 = pymol::string_format("%scmd.select(\"%s\",\"%s\",enable=1)\n", prefix, selName,
                    buffer);
            PLog(G, buf2, cPLog_no_flush);
          }
        } else {
          auto buffer = pymol::string_format("(%s(?%s))", sel_mode_kw, cTempRectSele);
          selectorResult = SelectorCreate(G, selName, buffer.c_str(), nullptr, 0, nullptr);
          if(log_box) {
            auto buf2 = pymol::string_format("%scmd.select(\"%s\",\"%s\",enable=1)\n", prefix, selName,
                    buffer);
            PLog(G, buf2, cPLog_no_flush);
          }
        }
      } else {
        if((mode == cButModeRectAdd) || (mode == cButModeSeleAddBox)) {
          auto buffer = pymol::string_format("%s(?%s)", sel_mode_kw, cTempRectSele);
          selectorResult = SelectorCreate(G, selName, buffer.c_str(), nullptr, 0, nullptr);
          if(log_box) {
            auto buf2 = pymol::string_format("%scmd.select(\"%s\",\"%s\",enable=1)\n", prefix, selName,
                    buffer);
            PLog(G, buf2, cPLog_no_flush);
          }
        } else if((mode == cButModeRectSub) || (mode == cButModeSeleSubBox)) {
          selectorResult = SelectorCreate(G, selName, "(none)", nullptr, 0, nullptr);
          if(log_box) {
            auto buf2 = pymol::string_format("%scmd.select(\"%s\",\"(none)\",enable=1)\n", prefix, selName);
            PLog(G, buf2, cPLog_no_flush);
          }
        } else {
          auto buffer = pymol::string_format("%s(?%s)", sel_mode_kw, cTempRectSele);
          selectorResult = SelectorCreate(G, selName, buffer.c_str(), nullptr, 0, nullptr);
          if(log_box) {
            auto buf2 = pymol::string_format("%scmd.select(\"%s\",\"%s\",enable=1)\n", prefix, selName,
                    buffer);
            PLog(G, buf2, cPLog_no_flush);
          }
        }
      }
      if(SettingGetGlobal_b(G, cSetting_auto_show_selections)) {
        ExecutiveSetObjVisib(G, selName, true, false);
      }
      break;
    }
    if(log_box) {
      auto buf2 = pymol::string_format("%scmd.delete(\"%s\")\n", prefix, cTempRectSele);
      PLog(G, buf2, cPLog_no_flush);
      PLogFlush(G);
    }
    ExecutiveDelete(G, cTempRectSele);
    WizardDoSelect(G, selName);
  } else {
    switch (mode) {
    case cButModeSeleSetBox:
      {
        ObjectNameType name;
        if(ExecutiveGetActiveSeleName(G, name, false, SettingGetGlobal_i(G, cSetting_logging))) {
          ExecutiveSetObjVisib(G, name, 0, false);
          if(SettingGetGlobal_i(G, cSetting_logging)) {
            auto buf2 = pymol::string_format("cmd.disable('%s')\n", name);
            PLog(G, buf2, cPLog_no_flush);
          }
        }
      }
      break;
    }
  }
  if (selectorResult) {
    result.currSelected = *selectorResult;
    result.netSelected = *selectorResult - beforeSelected;
    return result;
  } else {
    return selectorResult.error();
  }
}

pymol::Result<> ExecutiveTranslateAtom(
    PyMOLGlobals* G, const char* sele, const float* v, int state, int mode, int log)
{
  SETUP_SELE(sele, tmpsele0, sele0);

  {
    auto obj0 = SelectorGetSingleObjectMolecule(G, sele0);
    if(!obj0) {
      return pymol::make_error("Selection isn't a single atom.");
    } else {
      auto i0 = ObjectMoleculeGetAtomIndex(obj0, sele0);
      if(i0 < 0) {
        return pymol::make_error("Selection isn't a single atom.");
      } else {
        ObjectMoleculeMoveAtom(obj0, state, i0, v, mode, log);
      }
    }
  }
  return {};
}

/**
 * Performs a TTT-altering function whose transformation is optionally stored as a keyframe
 * @param name name of object or pattern
 * @param store to store transformation as keyframe
 * @param func function to call (ex. ExecutiveSetObjectTTT)
 * @param args arguments to func
 */

template<typename Func, typename... Args>
void ExecutiveObjectFuncTTT(PyMOLGlobals* G, pymol::zstring_view name, int store, Func func, Args... args)
{
  auto I = G->Executive;
  if (name.empty() || name == cKeywordAll || name == cKeywordSame) {
    for (auto& rec : pymol::make_list_adapter(I->Spec)) {
      switch (rec.type) {
      case cExecObject:
        {
          pymol::CObject *obj = rec.obj;
          if (ObjectGetSpecLevel(rec.obj, 0) >= 0 || name == cKeywordAll) {
            func(obj, args...);
            obj->invalidate(cRepNone, cRepInvExtents, -1);
          }
        }
        break;
      }
    }
    if(store && SettingGet<bool>(G,cSetting_movie_auto_interpolate)) {
      ExecutiveMotionReinterpolate(G);
    }
  } else { /* pattern */
    for (auto& rec : ExecutiveGetSpecRecsFromPattern(G, name)) {
      switch (rec.type) {
      case cExecObject:
        {
          pymol::CObject *obj = rec.obj;
          func(obj, args...);
          obj->invalidate(cRepNone, cRepInvExtents, -1);
        }
        break;
      }
    }
    if (store && SettingGet<bool>(G,cSetting_movie_auto_interpolate)) {
      ExecutiveMotionReinterpolate(G);
    }
  }
  SceneInvalidate(G);
}
pymol::Result<> ExecutiveCombineObjectTTT(PyMOLGlobals* G,
    pymol::zstring_view name, const float* ttt, int reverse_order, int store)
{
  ExecutiveObjectFuncTTT(G, name, store, ObjectCombineTTT, ttt, reverse_order, store);
  return {};
}

pymol::Result<> ExecutiveTranslateObjectTTT(PyMOLGlobals * G, pymol::zstring_view name, const float *trans, int store, int quiet)
{
  ExecutiveObjectFuncTTT(G, name, store, ObjectTranslateTTT, trans, store);
  return {};
}


pymol::Result<> ExecutiveSetObjectTTT(PyMOLGlobals * G, pymol::zstring_view name, const float *ttt, int state, int quiet, int store)
{
  ExecutiveObjectFuncTTT(G, name, store, ObjectSetTTT, ttt, state, store);
  return {};
}

int ExecutiveGetObjectTTT(PyMOLGlobals * G, const char *name, const float **ttt, int state, int quiet)
{
  pymol::CObject *obj = ExecutiveFindObjectByName(G, name);
  int ok = true;

  if(!obj) {
    PRINTFB(G, FB_ObjectMolecule, FB_Errors)
      "Error: object %s not found.\n", name ENDFB(G);
    ok = false;
  } else {
    ObjectGetTTT(obj, ttt, state);
  }
  return (ok);
}

pymol::Result<> ExecutiveTransformSelection(PyMOLGlobals* G,
    int state, const char* s1, int log, const float* ttt, int homogenous)
{
  SETUP_SELE_DEFAULT(1);

  auto vla = pymol::vla_take_ownership(SelectorGetObjectMoleculeVLA(G, sele1));
  if(!vla)
    return pymol::make_error("Could not find selection");
  for(auto obj : vla) {
    ObjectMoleculeTransformSelection(obj, state, sele1, ttt, log,
        tmpsele1->getName(), homogenous, true);
  }
  SceneInvalidate(G);
  return {};
}

pymol::Result<> ExecutiveTransformObjectSelection2(PyMOLGlobals* G,
    pymol::CObject* obj, int state, const char* s1, int log, const float* matrix,
    int homogenous, int global)
{
  int ok = true;

  switch (obj->type) {
  case cObjectMolecule:
    {
      int sele = -1;
      ObjectMolecule *objMol = (ObjectMolecule *) obj;

      if(s1 && s1[0]) {
        sele = SelectorIndexByName(G, s1);
        if(sele < 0)
          ok = false;
      }
      if(!ok) {
        return pymol::make_error("Selection object ", s1, " not found.");
      } else {
        ObjectMoleculeTransformSelection(objMol, state, sele, matrix, log, s1, homogenous,
                                         global);
      }
      EditorDihedralInvalid(G, objMol);
      SceneInvalidate(G);
    }
    break;
  default:
    if (auto* objstate = obj->getObjectState(state)) {
      double matrixd[116];
      if(homogenous) {
        convert44f44d(matrix, matrixd);
      } else {
        convertTTTfR44d(matrix, matrixd);
      }
      ObjectStateTransformMatrix(objstate, matrixd);
      obj->invalidate(cRepNone, cRepInvExtents, state);
    }
    break;
  }
  return {};
}

pymol::Result<> ExecutiveTransformObjectSelection(PyMOLGlobals* G,
    const char* name, int state, const char* s1, int log, const float* matrix,
    int homogenous, int global)
{
  pymol::CObject *obj = ExecutiveFindObjectByName(G, name);
  if(obj) {
    return ExecutiveTransformObjectSelection2(G, obj, state, s1, log, matrix, homogenous,
                                              global);
  }
  return {};
}

int ExecutiveValidName(PyMOLGlobals * G, const char *name)
{
  if(!ExecutiveFindSpec(G, name)) {
    int ignore_case = SettingGetGlobal_b(G, cSetting_ignore_case);

    if(!WordMatchExact(G, name, cKeywordAll, ignore_case))
      if(!WordMatchExact(G, name, cKeywordSame, ignore_case))
        if(!WordMatchExact(G, name, cKeywordCenter, ignore_case))
          if(!WordMatchExact(G, name, cKeywordOrigin, ignore_case))
            return false;
  }
  return true;
}

int ExecutivePhiPsi(PyMOLGlobals * G, const char *s1, ObjectMolecule *** objVLA, int **iVLA,
                    float **phiVLA, float **psiVLA, int state)
{
  SelectorTmp tmpsele1(G, s1);
  int sele1 = tmpsele1.getIndex();

  int result = false;
  ObjectMoleculeOpRec op1;
  if(sele1 >= 0) {
    ObjectMoleculeOpRecInit(&op1);
    op1.i1 = 0;
    op1.i2 = state;
    op1.obj1VLA = VLAlloc(ObjectMolecule *, 1000);
    op1.i1VLA = VLAlloc(int, 1000);
    op1.f1VLA = VLAlloc(float, 1000);
    op1.f2VLA = VLAlloc(float, 1000);
    op1.code = OMOP_PhiPsi;
    ExecutiveObjMolSeleOp(G, sele1, &op1);
    result = op1.i1;
    VLASize(op1.i1VLA, int, op1.i1);
    VLASize(op1.obj1VLA, ObjectMolecule *, op1.i1);
    VLASize(op1.f1VLA, float, op1.i1);
    VLASize(op1.f2VLA, float, op1.i1);
    *iVLA = op1.i1VLA;
    *objVLA = op1.obj1VLA;
    *phiVLA = op1.f1VLA;
    *psiVLA = op1.f2VLA;
  } else {
    *objVLA = NULL;
    *iVLA = NULL;
    *phiVLA = NULL;
    *psiVLA = NULL;
  }
  return (result);
}

int ExecutiveAlign(PyMOLGlobals * G, const char *s1, const char *s2, const char *mat_file, float gap,
                   float extend, int max_gap, int max_skip, float cutoff, int cycles,
                   int quiet, const char *oname, int state1, int state2,
                   ExecutiveRMSInfo * rms_info, int transform, int reset, float seq_wt,
                   float radius, float scale, float base, float coord_wt, float expect,
                   int window, float ante)
{
  int sele1 = SelectorIndexByName(G, s1);
  int sele2 = SelectorIndexByName(G, s2);
  int *vla1 = NULL;
  int *vla2 = NULL;
  int na, nb;
  int c;
  int ok = true;
  int use_sequence = (mat_file && mat_file[0] && (seq_wt != 0.0F));
  int use_structure = (seq_wt >= 0.0F); /* negative seq_wt means sequence only! */
  ObjectMolecule *mobile_obj = NULL;
  CMatch *match = NULL;

  if(!use_structure)
    window = 0;

  if((scale == 0.0F) && (seq_wt == 0.0F) && (ante < 0.0F) && window)
    ante = window;

  if(ante < 0.0F)
    ante = 0.0F;

  if((sele1 >= 0)) {
    mobile_obj = SelectorGetSingleObjectMolecule(G, sele1);
    if(!mobile_obj) {
      ok = false;
      PRINTFB(G, FB_Executive, FB_Errors)
        " ExecutiveAlign: mobile selection must derive from one object only.\n" ENDFB(G);
    }
  }
  if(ok && (sele1 >= 0) && (sele2 >= 0) && rms_info) {
    vla1 = SelectorGetResidueVLA(G, sele1, use_structure, NULL);
    vla2 = SelectorGetResidueVLA(G, sele2, use_structure, mobile_obj);
    if(vla1 && vla2) {
      na = VLAGetSize(vla1) / 3;
      nb = VLAGetSize(vla2) / 3;
      if(na && nb) {
        match = MatchNew(G, na, nb, window);
        if(match) {
          if(use_sequence) {
            if(ok)
              ok = MatchResidueToCode(match, vla1, na);
            if(ok)
              ok = MatchResidueToCode(match, vla2, nb);
            if(ok)
              ok = MatchMatrixFromFile(match, mat_file, quiet);
            if(ok)
              ok = MatchPreScore(match, vla1, na, vla2, nb, quiet);
          }
          if(use_structure) {
	    /* avoid degenerate alignments */
	    ok = ((na>1) && (nb>1) && ok);
            if(ok) {
              ok = SelectorResidueVLAsTo3DMatchScores(G, match,
                                                      vla1, na, state1,
                                                      vla2, nb, state2, seq_wt,
                                                      radius, scale, base,
                                                      coord_wt, expect);
	    } else {
	      PRINTFB(G, FB_Executive, FB_Errors)
		" ExecutiveAlign: No alignment found.\n" ENDFB(G);
	    }
	  }
          if(ok)
            ok = MatchAlign(match, gap, extend, max_gap, max_skip, quiet, window, ante);
          if(ok) {
            rms_info->raw_alignment_score = match->score;
            rms_info->n_residues_aligned = match->n_pair;
            if(match->pair) {

              c = SelectorCreateAlignments(G, match->pair,
                                           sele1, vla1, sele2, vla2,
                                           "_align1", "_align2", false, false);

              if(c) {
                int mode = 2;
                if(!quiet) {
                  PRINTFB(G, FB_Executive, FB_Actions)
                    " %s: %d atoms aligned.\n", __func__, c ENDFB(G);
                }
                if(oname && oname[0] && reset)
                  ExecutiveDelete(G, oname);
                if(!transform)
                  mode = 1;
                ok = ExecutiveRMS(G, "_align1", "_align2", mode, cutoff, cycles,
                                  quiet, oname, state1, state2, false, 0, rms_info);
              } else {
                if(!quiet) {
                  PRINTFB(G, FB_Executive, FB_Actions)
                    " ExecutiveAlign-Error: atomic alignment failed (mismatched identifiers?).\n"
                    ENDFB(G);
                }
                ok = false;
              }
            }
          }
          MatchFree(match);
        }
      } else {
        ok = false;
        PRINTFB(G, FB_Executive, FB_Errors)
          " ExecutiveAlign: invalid selections for alignment.\n" ENDFB(G);
      }
    }
    ExecutiveUpdateCoordDepends(G, mobile_obj); //Updates dynamic_measures - see PYMOL-3090
  }

  VLAFreeP(vla1);
  VLAFreeP(vla2);
  return ok;
}

/**
 * Implementation of `cmd.find_pairs()`
 *
 * @param s1 atom selection expression
 * @param s2 atom selection expression
 * @param state1 object state
 * @param state2 object state
 * @param mode 0: find all pairs, 1: find hydrogen bonds
 * @param cutoff maximum distance in Angstrom
 * @param h_angle hydrogen bond angle cutoff for mode=1
 * @param[out] indexVLA pairs of atom indices (0-based)
 * @param[out] objVLA pairs of objects
 * @return number of pairs
 */
pymol::Result<int>
ExecutivePairIndices(PyMOLGlobals * G, const char *s1, const char *s2, int state1, int state2,
                         int mode, float cutoff, float h_angle,
                         int **indexVLA, ObjectMolecule *** objVLA)
{
  SETUP_SELE_DEFAULT_PREFIXED(1, cSelectionInvalid);
  SETUP_SELE_DEFAULT_PREFIXED(2, sele1);

  return SelectorGetPairIndices(G, sele1, state1, sele2, state2,
                                    mode, cutoff, h_angle, indexVLA, objVLA);
}

pymol::Result<int> ExecutiveCartoon(PyMOLGlobals* G, int type, const char* s1)
{
  SETUP_SELE_DEFAULT(1);

  ObjectMoleculeOpRec op1;

  ObjectMoleculeOpRecInit(&op1);
  op1.i2 = 0;
  {
    op1.code = OMOP_Cartoon;
    op1.i1 = type;
    op1.i2 = 0;
    op1.i3 = 0;
    ExecutiveObjMolSeleOp(G, sele1, &op1);
    if (op1.i3>0){
        op1.code = OMOP_INVA;
        op1.i1 = cRepCartoonBit;
        op1.i2 = cRepInvRep;
        ExecutiveObjMolSeleOp(G, sele1, &op1);
    }
  }
  return (op1.i2);
}


/*========================================================================*/
float *ExecutiveGetVertexVLA(PyMOLGlobals * G, const char *s1, int state)
{
  /* returns NULL if none found */

  float *result = NULL;
  ObjectMoleculeOpRec op1;
  int sele1;
  sele1 = SelectorIndexByName(G, s1);
  if(sele1 >= 0) {
    ObjectMoleculeOpRecInit(&op1);
    op1.nvv1 = 0;
    op1.vv1 = VLAlloc(float, 1000);
    if(state >= 0) {
      op1.cs1 = state;
      op1.code = OMOP_SingleStateVertices;
    } else {
      op1.code = OMOP_VERT;
    }
    ExecutiveObjMolSeleOp(G, sele1, &op1);
    VLASize(op1.vv1, float, op1.nvv1 * 3);
    result = op1.vv1;
  }
  return (result);
}


/*========================================================================*/
#ifndef _PYMOL_NOPY
/**
 * @pre GIL
 */
PyObject *ExecutiveGetSettingOfType(PyMOLGlobals * G, int index,
                                    const char *object, int state, int type)
{
  assert(PyGILState_Check());

  PyObject *result = NULL;
  pymol::CObject *obj = NULL;
  CSetting *set_ptr1 = NULL, *set_ptr2 = NULL;
  pymol::copyable_ptr<CSetting>* handle = nullptr;

  if(object && object[0]) {
      obj = ExecutiveFindObjectByName(G, object);
      if(!obj)
        return PyErr_Format(P_CmdException, "object \"%s\" not found", object);
    handle = obj->getSettingHandle(-1);
    if(handle)
      set_ptr1 = handle->get();
    if(state >= 0) {
      handle = obj->getSettingHandle(state);
      if(handle)
        set_ptr2 = handle->get();
      else {
        return PyErr_Format(
            P_CmdException, "object \"%s\" lacks state %d", object, state + 1);
      }
    }
  }
  {
    switch (type) {
    case cSetting_boolean:
      {
        int value = SettingGet_b(G, set_ptr2, set_ptr1, index);
        result = PyBool_FromLong(value);
      }
      break;
    case cSetting_int:
      {
        int value = SettingGet_i(G, set_ptr2, set_ptr1, index);
        result = Py_BuildValue("i", value);
      }
      break;
    case cSetting_float:
      {
        float value = SettingGet_f(G, set_ptr2, set_ptr1, index);
        result = PyFloat_FromDouble(pymol::pretty_f2d(value));
      }
      break;
    case cSetting_float3:
      {
        const float * value = SettingGet_3fv(G, set_ptr2, set_ptr1, index);
        if (value) {
          result = Py_BuildValue("fff", pymol::pretty_f2d(value[0]),
              pymol::pretty_f2d(value[1]), pymol::pretty_f2d(value[2]));
        } else {
          PyErr_SetNone(PyExc_ValueError);
        }
      }
      break;
    case cSetting_color:
      {
        int value = SettingGet_color(G, set_ptr2, set_ptr1, index);
        result = Py_BuildValue("i", value);
      }
      break;
    case cSetting_string:
      {
        OrthoLineType buffer = "";
        result = Py_BuildValue("s",
            SettingGetTextPtr(G, set_ptr2, set_ptr1, index, buffer));
      }
      break;
    case cSetting_tuple:
      result = SettingGetTuple(G, set_ptr2, set_ptr1, index);
      break;
    default:
      PyErr_Format(PyExc_ValueError, "invalid setting type %d", type);
      break;
    }
  }
  return (result);
}
#endif


/*========================================================================*/
void ExecutiveSetLastObjectEdited(PyMOLGlobals * G, pymol::CObject * o)
{
  CExecutive *I = G->Executive;
  I->LastEdited = o;
}


/*========================================================================*/
pymol::CObject *ExecutiveGetLastObjectEdited(PyMOLGlobals * G)
{
  CExecutive *I = G->Executive;
  return (I->LastEdited);
}


/*========================================================================*/
int ExecutiveSaveUndo(PyMOLGlobals * G, const char *s1, int state)
  {
  int sele1;
  ObjectMoleculeOpRec op1;

  if(state < 0)
    state = SceneGetState(G);
  sele1 = SelectorIndexByName(G, s1);
  ObjectMoleculeOpRecInit(&op1);
  op1.i2 = 0;
  if(sele1 >= 0) {
    op1.code = OMOP_SaveUndo;
    op1.i1 = state;
    ExecutiveObjMolSeleOp(G, sele1, &op1);
  }
  return (op1.i2);
}


/*========================================================================*/
pymol::Result<> ExecutiveSetTitle(PyMOLGlobals * G, const char *name, int state, const char *text)
{
  ObjectMolecule *obj;
  obj = ExecutiveFindObjectMoleculeByName(G, name);
  if(!obj) {
    return pymol::make_error("Object ", name, " not found.");
  } else {
    auto res = ObjectMoleculeSetStateTitle(obj, state, text);
    if(!res){
      return res;
    }
  }
  SceneDirty(G);
  return {};
}


/*========================================================================*/
const char *ExecutiveGetTitle(PyMOLGlobals * G, const char *name, int state)
{
  ObjectMolecule *obj;
  obj = ExecutiveFindObjectMoleculeByName(G, name);
  if(!obj) {
    PRINTFB(G, FB_ObjectMolecule, FB_Errors)
      "Error: object %s not found.\n", name ENDFB(G);
    return NULL;
  }
  return ObjectMoleculeGetStateTitle(obj, state);
}


/*========================================================================*/
void ExecutiveHideSelections(PyMOLGlobals * G)
{
  CExecutive *I = G->Executive;
  SpecRec *rec = NULL;

  while(ListIterate(I->Spec, rec, next)) {
    if(rec->type == cExecSelection) {
      if(rec->visible) {
        rec->visible = false;
        SceneInvalidate(G);
        SeqDirty(G);
	ReportEnabledChange(G, rec);
      }
    }
  }
}

void ExecutiveInvalidateSelectionIndicators(PyMOLGlobals *G){
  CExecutive *I = G->Executive;
  ExecutiveInvalidateSelectionIndicatorsCGO(G);
  I->selectorTextureSize = 0;
}

void ExecutiveInvalidateSelectionIndicatorsCGO(PyMOLGlobals *G){
  CExecutive *I = G->Executive;
  SpecRec *rec = NULL;
  if (I){
    if (I->selIndicatorsCGO){
      CGOFree(I->selIndicatorsCGO);
      I->selIndicatorsCGO = 0;
    }
    while(ListIterate(I->Spec, rec, next)) {
      if(rec->type == cExecObject) {
	  CGOFree(rec->gridSlotSelIndicatorsCGO);	  
      }
    }
  }
}

static void ExecutiveRegenerateTextureForSelector(PyMOLGlobals *G, int round_points, int *widths_arg){
  CExecutive *I = G->Executive;
  unsigned char *temp_buffer = pymol::malloc<unsigned char>(widths_arg[0] * widths_arg[0] * 4);
  int a, b;
  float mid_point, disty, distx, dist, wminusd;
  float widths[] = { widths_arg[0]/2.f, widths_arg[1]/2.f, widths_arg[2]/2.f };
  unsigned char *q;
  mid_point = ((widths_arg[0]-1.f)/2.f);
  TextureInitTextTexture(G);
  if (I->selectorTextureAllocatedSize < widths_arg[0]){
    // need to re-allocate a new part of the texture for the selection indicator
    TextureGetPlacementForNewSubtexture(G, widths_arg[0], widths_arg[0], &I->selectorTexturePosX, &I->selectorTexturePosY);
    I->selectorTextureAllocatedSize = widths_arg[0];
  }
  q = temp_buffer;
  if (round_points){
    for(b = 0; b < widths_arg[0]; b++) {
      disty = fabs(mid_point - b);
      for(a = 0; a < widths_arg[0]; a++) {
	distx = fabs(mid_point - a);
	dist = sqrt(distx*distx + disty*disty);
	wminusd = widths[0] - dist;
	if (dist < widths[2]){
	  // white
	  q[0] = q[1] = q[2] = q[3] = 255;
	} else if (dist < widths[1]){
	  // black
	  q[0] = q[1] = q[2] = 0; q[3] = 255;
	} else if (fabs(wminusd) < .5f) {
	  // color plus transparency
	  q[0] = 255;
	  q[1] = 51;
	  q[2] = 153;
	  q[3] = (unsigned char)((wminusd + .5)*255);
	} else if (dist < widths[0]) {
	  // color (pink by default 1.f, .2f, .6f)
	  q[0] = 255;
	  q[1] = 51;
	  q[2] = 153;
          q[3] = 255;
	} else {
	  // black
	  q[0] = q[1] = q[2] = q[3] = 0;
	}
	//	printf("%3d ", q[1]);
	q += 4;
      }
      //      printf("\n");
    }
  } else {
    for(b = 0; b < widths_arg[0]; b++) {
      dist = disty = fabs(mid_point - b);
      for(a = 0; a < widths_arg[0]; a++) {
	distx = fabs(mid_point - a);
	dist = disty;
	if (distx > dist){
	  dist = distx;
	}
	if (dist < widths[2]){
	  // white
	  q[0] = q[1] = q[2] = q[3] = 255;
	} else if (dist < widths[1]){
	  // black
	  q[0] = q[1] = q[2] = 0; q[3] = 255;
	} else {
	  // color (pink by default 1.f, .2f, .6f)
	  q[0] = 255;
	  q[1] = 51;
	  q[2] = 153;
          q[3] = 255;
	}
	//	printf("%3d ", q[1]);
	q += 4;
      }
      //      printf("\n");
    }
  }
  glTexSubImage2D(GL_TEXTURE_2D, 0, I->selectorTexturePosX, I->selectorTexturePosY,
		  widths_arg[0], widths_arg[0], GL_RGBA, GL_UNSIGNED_BYTE, temp_buffer);  
  FreeP(temp_buffer);
}

static void ExecutiveRenderIndicatorCGO(PyMOLGlobals * G, CGO *selIndicatorsCGO){
  CExecutive *I = G->Executive;
  CShaderPrg *shaderPrg;
  float text_texture_dim = TextureGetTextTextureSize(G);
  float textureScale;
  int no_depth = (int) SettingGetGlobal_f(G, cSetting_selection_overlay);
  shaderPrg = G->ShaderMgr->Enable_IndicatorShader();
  if (!shaderPrg)return;
  glEnable(GL_POINT_SPRITE);
  glEnable(GL_VERTEX_PROGRAM_POINT_SIZE);
  shaderPrg->SetLightingEnabled(0);
  shaderPrg->SetAttrib4fLocation("a_Color", 1.f, 1.f, 1.f, 1.f);
  shaderPrg->Set1f("g_pointSize", DIP2PIXEL(I->selectorTextureSize));
  shaderPrg->Set2f("textureLookup", I->selectorTexturePosX/text_texture_dim, I->selectorTexturePosY/text_texture_dim);
  textureScale = I->selectorTextureSize/text_texture_dim ;
  shaderPrg->Set2f("textureScale", textureScale, textureScale);
  int v[4];
  glGetIntegerv(GL_VIEWPORT, v);
  shaderPrg->Set4f("viewport", v[0], v[1], v[2], v[3]);
  if(no_depth)
    glDisable(GL_DEPTH_TEST);
  CGORenderGL(selIndicatorsCGO, NULL, NULL, NULL, NULL, NULL);
  if(no_depth)
    glEnable(GL_DEPTH_TEST);
  glDisable(GL_VERTEX_PROGRAM_POINT_SIZE);
  glDisable(GL_POINT_SPRITE);
  shaderPrg->Disable();
}

static void ExecutiveSetupIndicatorPassGLImmediate(PyMOLGlobals * G, SpecRec *rec, int pass, float gl_width, int width){
#ifndef PURE_OPENGL_ES_2
  switch (pass){
  case 0:
    if(rec->sele_color < 0)
      glColor3f(1.0F, 0.2F, 0.6F);
    else
      glColor3fv(ColorGet(G, rec->sele_color));
    glPointSize(gl_width);
    break;
  case 1:
    if(width > 2) {
      switch (width) {
      case 1: case 2: case 3:
	glPointSize(1.0F);
	break;
      case 4:
	glPointSize(2.0F);
	break;
      case 5:
	glPointSize(3.0F);
	break;
      case 6: case 7: case 8: case 9:
	glPointSize(4.0F);
	break;
      default:
	glPointSize(6.0F);
	break;
      }
      glColor3f(0.0F, 0.0F, 0.0F);
      break;
    }
  case 2:
    if(width > 4) {
      if(width > 5) {
	glPointSize(2.0F);
      } else {
	glPointSize(1.0F);
      }
      glColor3f(1.0F, 1.0F, 1.0F);
    }
    break;
  }
#endif
}

/*========================================================================*/
void ExecutiveRenderSelections(PyMOLGlobals * G, int curState, int slot, GridInfo *grid)
{
  CExecutive *I = G->Executive;
  SpecRec *rec = NULL;
  int any_active = false;
  int no_depth = (int) SettingGetGlobal_f(G, cSetting_selection_overlay);
  short use_shader = (short) SettingGetGlobal_b(G, cSetting_use_shaders);
  float min_width;
  float gl_width;
  int width;
  int max_width = (int) SettingGetGlobal_f(G, cSetting_selection_width_max);
  float width_scale = SettingGetGlobal_f(G, cSetting_selection_width_scale);
  int round_points = SettingGetGlobal_b(G, cSetting_selection_round_points);
  short inv_indicators = false;
  min_width = SettingGetGlobal_f(G, cSetting_selection_width);
  if(width_scale >= 0.0F) {
    width = (int) ((width_scale *
		    fabs(SettingGetGlobal_f(G, cSetting_stick_radius)) /
		    SceneGetScreenVertexScale(G, NULL)));
    if(width < min_width)
      width = (int) min_width;
    else if(width > max_width)
      width = (int) max_width;
  } else
    width = (int) min_width;
  if(round_points) {
    width = (int) (width * 1.44F);
  }
  gl_width = (float) width;
  if(width > 6) {   /* keep it even above 6 */
    if(width & 0x1) {
      width--;
      gl_width = (float) width;
    }
  }

  if (use_shader){
    if (I->selectorTextureSize!=(int)gl_width || round_points != I->selectorIsRound ){
      inv_indicators = true;
    } else {
      if (slot){
	while(ListIterate(I->Spec, rec, next)) {
	  if(rec->type == cExecObject) {
	    if (SceneGetDrawFlagGrid(G, grid, rec->obj->grid_slot)){
	      if (rec->gridSlotSelIndicatorsCGO){
		ExecutiveRenderIndicatorCGO(G, rec->gridSlotSelIndicatorsCGO);
	      }
	    }
	  }
	}
	rec = NULL;
      } else if (I->selIndicatorsCGO){
	ExecutiveRenderIndicatorCGO(G, I->selIndicatorsCGO);
	return;
      }
    }
  } else {
    inv_indicators = true;
  }
  if (inv_indicators){
      CGOFree(I->selIndicatorsCGO);
    if (slot){
      while(ListIterate(I->Spec, rec, next)) {
	if(rec->type == cExecObject) {
	  if (SceneGetDrawFlagGrid(G, grid, rec->obj->grid_slot)){
	      CGOFree(rec->gridSlotSelIndicatorsCGO);
	  }
	}
      }
    }
  }
  while(ListIterate(I->Spec, rec, next)) {
    if(rec->type == cExecSelection) {
      if (rec->visible) {
        any_active = true;
        break;
      }
    }
  }

  if(any_active) {
    SpecRec *rec1;
    int sele;
    int vis_only = SettingGetGlobal_b(G, cSetting_selection_visible_only);
    int fog = SettingGetGlobal_b(G, cSetting_depth_cue) && SettingGetGlobal_f(G, cSetting_fog);

    (void)fog;

    rec = NULL;

    if(round_points) {
      glEnable(GL_POINT_SMOOTH);
      glAlphaFunc(GL_GREATER, 0.5F);
      glEnable(GL_ALPHA_TEST);
      glHint(GL_POINT_SMOOTH_HINT, GL_NICEST);
    } else {
      glDisable(GL_POINT_SMOOTH);
      glDisable(GL_ALPHA_TEST);
      glHint(GL_POINT_SMOOTH_HINT, GL_FASTEST);
    }

    if (use_shader && (I->selectorTextureSize!=(int)gl_width || round_points != I->selectorIsRound )){
      int level_widths[] = { 0, 0, 0 };
      level_widths[0] = (int)gl_width;
      if(width > 2) {
	switch (width) {
	case 1: case 2: case 3:
	  level_widths[1] = 1;
	  break;
	case 4:
	  level_widths[1] = 2;
	  break;
	case 5:
	  level_widths[1] = 3;
	  break;
	case 6: case 7: case 8: case 9:
	  level_widths[1] = 4;
	  break;
	default:
	  level_widths[1] = 6;
	  break;
	}
      }
      if(width > 4) {
	if(width > 5) {
	  level_widths[2] = 2;
	} else {
	  level_widths[2] = 1;
	}
      }
      ExecutiveInvalidateSelectionIndicatorsCGO(G);
      ExecutiveRegenerateTextureForSelector(G, round_points, level_widths);
      I->selectorTextureSize = gl_width;
      I->selectorIsRound = round_points;
    }

    while(ListIterate(I->Spec, rec, next)) {
      if(rec->type == cExecSelection) {

        if(rec->visible) {
          int enabled = true;
          SpecRec *group_rec = rec->group;
          while(enabled && group_rec) {
            if(!group_rec->visible)
              enabled = false;
            else
              group_rec = group_rec->group;
          }

          if(enabled) {

            sele = SelectorIndexByName(G, rec->name);   /* TODO: speed this up */
            if(sele >= 0) {
	      int pass, totpass = use_shader ? 1 : 3;

	      if (!use_shader){
		if(no_depth)
		  glDisable(GL_DEPTH_TEST);
		glDisable(GL_FOG);
	      }
	      for (pass=0; pass<totpass; pass++){
		if (use_shader){ // there is only one pass anyway
		  if (!slot){
		    I->selIndicatorsCGO = CGONew(G);
		    CGODotwidth(I->selIndicatorsCGO, gl_width);
		    CGOBegin(I->selIndicatorsCGO, GL_POINTS);
		  }
		} else {
		  ExecutiveSetupIndicatorPassGLImmediate(G, rec, pass, gl_width, width);
		  glBegin(GL_POINTS);
		}
		rec1 = NULL;
		while(ListIterate(I->Spec, rec1, next)) {
		  if(rec1->type == cExecObject) {
		    if(rec1->obj->type == cObjectMolecule) {
		      if (SceneGetDrawFlagGrid(G, grid, rec1->obj->grid_slot)){
			if (!use_shader || !slot){
			  ObjectMoleculeRenderSele((ObjectMolecule *) rec1->obj, curState, sele,
						   vis_only SELINDICATORARGVAR);
			} else if (!rec1->gridSlotSelIndicatorsCGO){
			  if (!use_shader){
			    ObjectMoleculeRenderSele((ObjectMolecule *) rec1->obj, curState, sele,
						     vis_only, NULL);
			  } else {
			    rec1->gridSlotSelIndicatorsCGO = CGONew(G);
			    CGODotwidth(rec1->gridSlotSelIndicatorsCGO, gl_width);
			    CGOBegin(rec1->gridSlotSelIndicatorsCGO, GL_POINTS);
			    ObjectMoleculeRenderSele((ObjectMolecule *) rec1->obj, curState, sele,
						     vis_only, rec1->gridSlotSelIndicatorsCGO);
			    CGOEnd(rec1->gridSlotSelIndicatorsCGO);
			    CGOStop(rec1->gridSlotSelIndicatorsCGO);
                            CGO* optimized = CGOOptimizeToVBONotIndexedNoShader(
                                rec1->gridSlotSelIndicatorsCGO);
                            CGOFree(rec1->gridSlotSelIndicatorsCGO);
                            rec1->gridSlotSelIndicatorsCGO = optimized;
                            assert(rec1->gridSlotSelIndicatorsCGO->use_shader);
                            ExecutiveRenderIndicatorCGO(G, rec1->gridSlotSelIndicatorsCGO);
			  }
			}
		      } else {
			CGOFree(rec1->gridSlotSelIndicatorsCGO);
		      }
		    }
		  }
		}
		if (use_shader){
		  if (!slot){
		    CGOEnd(I->selIndicatorsCGO);
		    CGOStop(I->selIndicatorsCGO);
                    CGO* optimized =
                        CGOOptimizeToVBONotIndexedNoShader(I->selIndicatorsCGO);
                    CGOFree(I->selIndicatorsCGO);
                    if ((I->selIndicatorsCGO = optimized)) {
                      assert(I->selIndicatorsCGO->use_shader);
                      ExecutiveRenderSelections(G, curState, slot, grid);
                    }
                    return;
		  }
		} else {
		  glEnd();
		}
	      }
	      if(fog)
		glEnable(GL_FOG);
	      if(no_depth)
		glEnable(GL_DEPTH_TEST);
            }
          }
        }
      }
    }
    if(!use_shader && round_points) {
      glAlphaFunc(GL_GREATER, 0.05F);
    }
  }
}

/*========================================================================*/

/**
 * Macro with error handling to get a single atom vertex and associated temp
 * selection. Assigns `vertexN`, `tmpseleN`, and `resultN`
 */
#define SETUP_SINGLE_ATOM_VERTEX(N)                                            \
  auto tmpsele##N = SelectorTmp::make(G, s##N);                                \
  p_return_if_error_prefixed(tmpsele##N, "Selection " #N ": ");                \
  auto result##N =                                                             \
      SelectorGetSingleAtomVertex(G, (tmpsele##N)->getIndex(), state);         \
  p_return_if_error_prefixed(result##N, "Selection " #N ": ");                 \
  const float* vertex##N = (result##N).result().data()

/**
 * Get the distance between two atoms.
 * Implementation of `cmd.get_distance()`
 *
 * @param s0 single-atom selection expression
 * @param s1 single-atom selection expression
 * @param state object state
 * @return distance in Angstrom
 */
pymol::Result<float> ExecutiveGetDistance(
    PyMOLGlobals* G, const char* s1, const char* s2, int state)
{
  SETUP_SINGLE_ATOM_VERTEX(1);
  SETUP_SINGLE_ATOM_VERTEX(2);
  return static_cast<float>(diff3f(vertex1, vertex2));
}


/*========================================================================*/
/**
 * Get the angle between tree atoms s0-s1-s2.
 * Implementation of `cmd.get_angle()`
 *
 * @param s0 single-atom selection expression
 * @param s1 single-atom selection expression
 * @param s2 single-atom selection expression
 * @param state object state
 * @return angle in degrees
 */
pymol::Result<float> ExecutiveGetAngle(
    PyMOLGlobals* G, const char* s1, const char* s2, const char* s3, int state)
{
  SETUP_SINGLE_ATOM_VERTEX(1);
  SETUP_SINGLE_ATOM_VERTEX(2);
  SETUP_SINGLE_ATOM_VERTEX(3);

  float d1[3], d2[3];
  subtract3f(vertex1, vertex2, d1);
  subtract3f(vertex3, vertex2, d2);
  return rad_to_deg(get_angle3f(d1, d2));
}


/*========================================================================*/
/**
 * Get the dihedral angle between four atoms.
 * Implementation of `cmd.get_dihedral()`
 *
 * @param s0 single-atom selection expression
 * @param s1 single-atom selection expression
 * @param s2 single-atom selection expression
 * @param s3 single-atom selection expression
 * @param state object state
 * @return dihedral angle in degrees
 */
pymol::Result<float> ExecutiveGetDihe(PyMOLGlobals* G, const char* s1,
    const char* s2, const char* s3, const char* s4, int state)
{
  SETUP_SINGLE_ATOM_VERTEX(1);
  SETUP_SINGLE_ATOM_VERTEX(2);
  SETUP_SINGLE_ATOM_VERTEX(3);
  SETUP_SINGLE_ATOM_VERTEX(4);
  return rad_to_deg(get_dihedral3f(vertex1, vertex2, vertex3, vertex4));
}


/*========================================================================*/
pymol::Result<> ExecutiveSetDihe(PyMOLGlobals* G, const char* s1,
    const char* s2, const char* s3, const char* s4, float value, int state,
    int quiet)
{
  SETUP_SINGLE_ATOM_VERTEX(1);
  SETUP_SINGLE_ATOM_VERTEX(2);
  SETUP_SINGLE_ATOM_VERTEX(3);
  SETUP_SINGLE_ATOM_VERTEX(4);
  float current = rad_to_deg(get_dihedral3f(vertex1, vertex2, vertex3, vertex4));
  float change = value - current;
  auto save_state = SceneGetState(G);
  SceneSetFrame(G, -1, state);        /* KLUDGE ALERT!
                                       * necessary because the editor 
                                       * can only work on the current state...this
                                       * needs to be changed.*/
  EditorSelect(G,
      tmpsele3->getName(),
      tmpsele2->getName(),
      "", "", false, true, true);
  EditorTorsion(G, change);
  SceneSetFrame(G, -1, save_state);
  if(!quiet) {
    PRINTFB(G, FB_Editor, FB_Actions)
      " SetDihedral: adjusted to %5.3f\n", value ENDFB(G);
  }
  return {};
}


/*========================================================================*/
pymol::Result<float> ExecutiveGetArea(
    PyMOLGlobals* G, const char* sele, int state, bool load_b)
{
  SETUP_SELE(sele, tmpsele0, sele0);

  auto obj0 = SelectorGetSingleObjectMolecule(G, sele0);
  if (!obj0) {
    if (SelectorCountAtoms(G, sele0, state) > 0)
      return pymol::Error("Selection must be within a single object");
    return 0.f;
  }

  auto cs = obj0->getCoordSet(state);
  if (!cs)
    return pymol::Error("Invalid state");

  auto rep = (RepDot*) RepDotDoNew(cs, cRepDotAreaType, state);
  if (!rep)
    return pymol::Error("Can't get dot representation.");

  if (load_b) {
    /* zero out B-values within selection */
    ObjectMoleculeOpRec op;
    ObjectMoleculeOpRecInit(&op);
    op.code = OMOP_SetB;
    op.f1 = 0.0;
    op.i1 = 0;
    ExecutiveObjMolSeleOp(G, sele0, &op);
  }

  float const* const area = rep->A;
  int const* const ati = rep->Atom;

  AtomInfoType* ai = nullptr;
  int known_member = -1;
  bool is_member = false;
  float result = 0.f;

  for (int a = 0; a < rep->N; ++a) {
    if (known_member != ati[a]) {
      known_member = ati[a];
      ai = obj0->AtomInfo + known_member;
      is_member = SelectorIsMember(G, ai->selEntry, sele0);
    }

    if (is_member) {
      result += area[a];
      if (load_b)
        ai->b += area[a];
    }
  }

  delete rep;
  return result;
}

/*========================================================================*/
/**
 * Implementation of `cmd.get_names()`
 *
 * @param mode see modules/pymol/querying.py
 * @param enabled_only only include enabled names
 * @param s0 atom selection expression, limits results to molecular objects and
 * named selections which have at least one atom in the selection.
 */
pymol::Result<std::vector<const char*>> ExecutiveGetNames(
    PyMOLGlobals* G, int mode, int enabled_only, const char* s0)
{
  enum {
    cGetNames_all = 0,
    cGetNames_objects = 1,
    cGetNames_selections = 2,
    cGetNames_public = 3,
    cGetNames_public_objects = 4,
    cGetNames_public_selections = 5,
    cGetNames_public_nongroup_objects = 6,
    cGetNames_public_group_objects = 7,
    cGetNames_nongroup_objects = 8,
    cGetNames_group_objects = 9,
  };

  bool include_objects = //
      mode == cGetNames_all || mode == cGetNames_objects ||
      mode == cGetNames_public || mode == cGetNames_public_objects;
  bool include_selections =
      mode == cGetNames_all || mode == cGetNames_selections ||
      mode == cGetNames_public || mode == cGetNames_public_selections;
  bool include_group_objects =
      mode == cGetNames_group_objects || mode == cGetNames_public_group_objects;
  bool include_nongroup_objects = //
      mode == cGetNames_nongroup_objects ||
      mode == cGetNames_public_nongroup_objects;
  bool public_only =
      mode >= cGetNames_public && mode <= cGetNames_public_group_objects;

  CExecutive *I = G->Executive;
  SpecRec *rec = NULL;
  std::vector<const char*> result;
  pymol::Result<SelectorTmp> tmpsele0;
  SelectorID_t sele0 = cSelectionInvalid;
  int incl_flag;
  if(s0[0]) {
    tmpsele0 = SelectorTmp::make(G, s0);
    p_return_if_error(tmpsele0);
    sele0 = tmpsele0->getIndex();
    assert(sele0 != cSelectionInvalid);
  }

  while(ListIterate(I->Spec, rec, next)) {
    incl_flag = 0;
    if((rec->type == cExecObject
        && (include_objects
            || (rec->obj->type != cObjectGroup && include_nongroup_objects)
            || (rec->obj->type == cObjectGroup && include_group_objects)))
       || (rec->type == cExecSelection && include_selections)
      ) {
      if(!public_only || rec->name[0] != '_') {
        if((!enabled_only) || (rec->visible)) {
          incl_flag = 0;
          if(sele0 < 0)
            incl_flag = 1;
          else
            switch (rec->type) {
            case cExecObject:
              if(rec->obj->type == cObjectMolecule) {
                int a;
                ObjectMolecule *obj_mol = (ObjectMolecule *) rec->obj;
                const AtomInfoType *ai = obj_mol->AtomInfo.data();
                for(a = 0; a < obj_mol->NAtom; a++) {
                  if(SelectorIsMember(G, ai->selEntry, sele0)) {
                    incl_flag = 1;
                    break;
                  }
                  ai++;
                }
              }
              break;
            case cExecSelection:
              if(SelectorCheckIntersection(G, sele0, SelectorIndexByName(G, rec->name))) {
                incl_flag = 1;
                break;
              }
              break;
            }
          if(incl_flag) {
            result.push_back(rec->name);
          }
        }
      }
    }
  }
  return (result);
}


/*========================================================================*/
/**
 * Return true if `name` is the name of an object molecule or a named selection
 */
bool ExecutiveIsMoleculeOrSelection(PyMOLGlobals * G, const char *name)
{
  if (!strcmp(name, cKeywordAll))
    return true;

  SpecRec *rec = ExecutiveFindSpec(G, name);
  if (rec && (
        (rec->type == cExecObject && rec->obj->type == cObjectMolecule) ||
        (rec->type == cExecSelection)))
    return true;

  return false;
}


/*========================================================================*/
/**
 * @param name Name of the object or selection
 * @return Type label like "object:molecule" or "selection"
 */
pymol::Result<char const*> ExecutiveGetType(PyMOLGlobals* G, const char* name)
{
  auto rec = ExecutiveFindSpec(G, name);

  if (!rec) {
    return pymol::Error("object not found");
  }

  if (rec->type == cExecSelection) {
    return "selection";
  }

  if (rec->type != cExecObject) {
    return "";
  }

  switch (rec->obj->type) {
  case cObjectMolecule:
    return "object:molecule";
  case cObjectMap:
    return "object:map";
  case cObjectMesh:
    return "object:mesh";
  case cObjectSlice:
    return "object:slice";
  case cObjectSurface:
    return "object:surface";
  case cObjectMeasurement:
    return "object:measurement";
  case cObjectCGO:
    return "object:cgo";
  case cObjectGroup:
    return "object:group";
  case cObjectVolume:
    return "object:volume";
  case cObjectAlignment:
    return "object:alignment";
  case cObjectGadget:
    return "object:ramp";
  default:
    return "object:";
  }
}

/*========================================================================*/
pymol::Result<> ExecutiveUpdateCmd(PyMOLGlobals* G, const char* s0,
    const char* s1, int sta0, int sta1, int method, int quiet)
{
  SETUP_SELE_DEFAULT_PREFIXED(0, cSelectionInvalid);
  SETUP_SELE_DEFAULT_PREFIXED(1, cSelectionInvalid);

  return SelectorUpdateCmd(G, sele0, sele1, sta0, sta1, method, quiet);
}


/*========================================================================*/
pymol::Result<> ExecutiveRenameObjectAtoms(
    PyMOLGlobals* G, const char* selection, int force, int quiet)
{
  SETUP_SELE(selection, tmpsele1, sele);

  {
    ObjectMoleculeOpRec op;
    ObjectMoleculeOpRecInit(&op);
    op.code = OMOP_RenameAtoms;
    op.i1 = 0;
    op.i2 = force;
    ExecutiveObjMolSeleOp(G, sele, &op);

    if(!quiet) {
      PRINTFB(G, FB_Executive, FB_Actions)
        " Rename: renamed %d atoms.\n", op.i1 ENDFB(G);
    }
  }
  return {};
}


/*========================================================================*/
pymol::Result<> ExecutiveFuse(PyMOLGlobals* G, const char* s0, const char* s1,
    int mode, int recolor, int move_flag)
{
  int i0 = -1;
  int i1 = -1;
  SelectorID_t sele2 = cSelectionInvalid;
  ObjectMolecule *obj0, *obj1;
  ObjectMoleculeOpRec op;
#define tmp_fuse_sele "tmp_fuse_sele"

  SETUP_SELE_DEFAULT_PREFIXED(0, cSelectionInvalid);
  SETUP_SELE_DEFAULT_PREFIXED(1, cSelectionInvalid);

  {
    {
      EditorInactivate(G);
      obj0 = SelectorGetSingleObjectMolecule(G, sele0);
      obj1 = SelectorGetSingleObjectMolecule(G, sele1);
      if(obj0)
        i0 = ObjectMoleculeGetAtomIndex(obj0, sele0);
      if(obj1)
        i1 = ObjectMoleculeGetAtomIndex(obj1, sele1);
      if(obj0 && obj1 && (i0 >= 0) && (i1 >= 0) && (obj0 != obj1)) {
        ObjectMoleculeVerifyChemistry(obj0, -1);
        ObjectMoleculeVerifyChemistry(obj1, -1);

        SelectorCreate(G, tmp_fuse_sele, NULL, obj0, 1, NULL);
        sele2 = SelectorIndexByName(G, tmp_fuse_sele);
        if(mode) {
          ObjectMoleculeOpRecInit(&op);
          op.code = OMOP_PrepareFromTemplate;
          op.ai = obj1->AtomInfo + i1;
          op.i1 = mode;
          op.i2 = 0;
          op.i3 = recolor;
          if(recolor)
            op.i4 = obj1->Color;
          ExecutiveObjMolSeleOp(G, sele2, &op);
        }
        SelectorDelete(G, tmp_fuse_sele);

        p_return_if_error(
            ObjectMoleculeFuse(obj1, i1, obj0, i0, mode != 3, move_flag));
      }
    }
  }
  return {};
}


/*========================================================================*/
pymol::Result<> ExecutiveSpheroid(PyMOLGlobals * G, const char *name, int average)
{                               /* EXPERIMENTAL */
  CExecutive *I = G->Executive;
  pymol::CObject *os = NULL;

  if(strlen(name)) {
    os = ExecutiveFindObjectByName(G, name);
    if(!os)
      return pymol::make_error("Object not found.");
    else if(os->type != cObjectMolecule) {
      return pymol::make_error("Bad object type.");
    }
  }

  if(os || (!strlen(name))) {   /* sort one or all */
    SpecRec* rec = nullptr;
    while(ListIterate(I->Spec, rec, next)) {
      if(rec->type == cExecObject)
        if(rec->obj->type == cObjectMolecule)
          if((!os) || (rec->obj == os)) {
            auto obj = (ObjectMolecule *) rec->obj;
            ObjectMoleculeCreateSpheroid(obj, average);
            obj->invalidate(cRepAll, cRepInvRep, -1);
          }
    }
    SceneChanged(G);
  }
  return {};
}


/*========================================================================*/
void ExecutiveRebuildAll(PyMOLGlobals * G)
{
  CExecutive *I = G->Executive;
  SpecRec *rec = NULL;
  PRINTFD(G, FB_Executive)
    " ExecutiveRebuildAll: entered.\n" ENDFD;
  auto defer_builds_mode = SettingGet<bool>(G, cSetting_defer_builds_mode);
  while(ListIterate(I->Spec, rec, next)) {
    if(rec->type == cExecObject) {
      auto level = cRepInvAll;
      switch (rec->obj->type) {
      case cObjectMeasurement:
        ObjectDistInvalidateRep((ObjectDist *) rec->obj, cRepAll);
        break;
      case cObjectMolecule:
        level = defer_builds_mode ? cRepInvPurge : cRepInvRep;
      case cObjectSurface:
      case cObjectMesh:
      case cObjectSlice:
      case cObjectAlignment:
      case cObjectCGO:
        rec->obj->invalidate(cRepAll, level, -1);
        break;
      }
    }
  }
  SeqChanged(G);
  SceneChanged(G);
}


/*========================================================================*/
void ExecutiveRebuildAllObjectDist(PyMOLGlobals * G)
{
  CExecutive *I = G->Executive;
  SpecRec *rec = NULL;
  while(ListIterate(I->Spec, rec, next)) {
    if(rec->type == cExecObject) {
      if(rec->obj->type == cObjectMeasurement) {
        ObjectDistInvalidateRep((ObjectDist *) rec->obj, cRepAll);
      }
    }
  }
  SceneInvalidate(G);
}


/*========================================================================*/
void ExecutiveUndo(PyMOLGlobals * G, int dir)
{
  CExecutive *I = G->Executive;
  ObjectMolecule *obj = NULL, *compObj;
  SpecRec *rec = NULL;

  auto o = ExecutiveGetLastObjectEdited(G);
  PRINTFB(G, FB_Executive, FB_Debugging)
  " ExecutiveUndo: last object %p\n", (void *) o ENDFB(G);
  if(o)
    if(o->type == cObjectMolecule)
      obj = (ObjectMolecule *) o;
  /* make sure this is still a real object */
  if(obj) {
    while(ListIterate(I->Spec, rec, next)) {
      if(rec->type == cExecObject)
        if(rec->obj->type == cObjectMolecule) {
          compObj = (ObjectMolecule *) rec->obj;
          if(obj == compObj) {
            ObjectMoleculeUndo(obj, dir);
            break;
          }
        }
    }
  }

}


/*========================================================================*/
pymol::Result<> ExecutiveSort(PyMOLGlobals* G, const char* name)
{
  CExecutive *I = G->Executive;
  ObjectMolecule *obj;
  SpecRec *rec = NULL;
  ObjectMoleculeOpRec op;
  int sele;
  int ok = true;
  if((!name) || (!name[0]))
    name = cKeywordAll;

  {
    CTracker *I_Tracker = I->Tracker;
    int list_id = ExecutiveGetNamesListFromPattern(G, name, true, true);
    int iter_id = TrackerNewIter(I_Tracker, 0, list_id);
    int changed = false;
    while(TrackerIterNextCandInList(I_Tracker, iter_id, (TrackerRef **) (void *) &rec)) {
      if(rec) {
        switch (rec->type) {
        case cExecAll:
          rec = NULL;
          while(ListIterate(I->Spec, rec, next)) {
            if((rec->type == cExecObject) && (rec->obj->type == cObjectMolecule)) {
              obj = (ObjectMolecule *) rec->obj;
              if (ok)
		ok &= ObjectMoleculeSort(obj);
	      if (ok){
		changed = true;
		sele = SelectorIndexByName(G, rec->name);
		if(sele >= 0) {
		  ObjectMoleculeOpRecInit(&op);
		  op.code = OMOP_INVA;
		  op.i1 = cRepCartoonBit | cRepRibbonBit;
		  op.i2 = cRepInvRep;
		  ExecutiveObjMolSeleOp(G, sele, &op);
		}
	      }
            }
          }
          break;
        case cExecSelection:
          sele = SelectorIndexByName(G, rec->name);
          if(sele >= 0) {
            op.code = OMOP_Sort;
            ExecutiveObjMolSeleOp(G, sele, &op);
            ObjectMoleculeOpRecInit(&op);
            op.code = OMOP_INVA;
            op.i1 = cRepCartoonBit | cRepRibbonBit;
            op.i2 = cRepInvRep;
            ExecutiveObjMolSeleOp(G, sele, &op);
            ObjectMoleculeOpRecInit(&op);
          }
          break;
        case cExecObject:
          if(rec->obj->type == cObjectMolecule) {
            obj = (ObjectMolecule *) rec->obj;
            if (ok)
	      ok &= ObjectMoleculeSort(obj);
            changed = true;
            sele = SelectorIndexByName(G, rec->name);
            if(sele >= 0) {
              ObjectMoleculeOpRecInit(&op);
              op.code = OMOP_INVA;
              op.i1 = cRepCartoonBit | cRepRibbonBit;
              op.i2 = cRepInvRep;
              ExecutiveObjMolSeleOp(G, sele, &op);
            }
          }
          break;
        }
      }
    }
    TrackerDelList(I_Tracker, list_id);
    TrackerDelIter(I_Tracker, iter_id);
    if(changed)
      SceneChanged(G);
  }
  return {};
}


/*========================================================================*/
pymol::Result<> ExecutiveRemoveAtoms(PyMOLGlobals * G, const char *s1, int quiet)
{
  SETUP_SELE_DEFAULT(1);

  CExecutive *I = G->Executive;
  SpecRec *rec = NULL;
  ObjectMolecule *obj = NULL;
  ObjectMoleculeOpRec op;

  {
    while(ListIterate(I->Spec, rec, next)) {
      if(rec->type == cExecObject) {
        if(rec->obj->type == cObjectMolecule) {
          ObjectMoleculeOpRecInit(&op);
          op.code = OMOP_Remove;
          op.i1 = 0;
          obj = (ObjectMolecule *) rec->obj;
          ObjectMoleculeVerifyChemistry(obj, -1);       /* remember chemistry for later */
          ObjectMoleculeSeleOp(obj, sele1, &op);
          if(op.i1) {
            if(!quiet) {
              PRINTFD(G, FB_Editor)
                " ExecutiveRemove-Debug: purging %i of %i atoms in %s\n",
                op.i1, obj->NAtom, obj->Name ENDFD;
            }
            ObjectMoleculePurge(obj);
            if(!quiet) {
              PRINTFB(G, FB_Editor, FB_Actions)
                " Remove: eliminated %d atoms in model \"%s\".\n",
                op.i1, obj->Name ENDFB(G);
            }
          }
        }
      }
    }
  }
  EditorRemoveStale(G);
  return {};
}


/*========================================================================*/
pymol::Result<> ExecutiveAddHydrogens(
    PyMOLGlobals* G, const char* s1, int quiet, int state, bool legacy)
{
  ObjectMoleculeOpRec op;

  if (legacy) {
    PRINTFB(G, FB_Executive, FB_Warnings)
      " %s-Warning: legacy mode was removed\n", __FUNCTION__ ENDFB(G);
  }

  SETUP_SELE_DEFAULT(1);

  {
    ObjectMoleculeOpRecInit(&op);
    op.code = OMOP_AddHydrogens;
    op.i1 = state;
    ExecutiveObjMolSeleOp(G, sele1, &op);
  }
  return {};
}


/*========================================================================*/
void ExecutiveFixHydrogens(PyMOLGlobals * G, const char *s1, int quiet)
{
  int sele1;
  ObjectMoleculeOpRec op;

  sele1 = SelectorIndexByName(G, s1);
  if(sele1 >= 0) {
    ObjectMoleculeOpRecInit(&op);
    op.code = OMOP_FixHydrogens;
    ExecutiveObjMolSeleOp(G, sele1, &op);
  }
}


/*========================================================================*/
pymol::Result<> ExecutiveFlag(PyMOLGlobals * G, int flag, const char* sele, int action, int quiet)
{
  enum {
    cFlagActionReset = 0,
    cFlagActionSet = 1,
    cFlagActionClear = 2,
  };

  if (flag < 0 || flag > 31) {
    return pymol::make_error("flag ", flag, " out of range [0, 31]");
  }

  ObjectMoleculeOpRec op;
  ObjectMoleculeOpRecInit(&op);
  switch (action) {
  case cFlagActionReset:
    op.code = OMOP_Flag;
    break;
  case cFlagActionSet:
    op.code = OMOP_FlagSet;
    break;
  case cFlagActionClear:
    op.code = OMOP_FlagClear;
    break;
  default:
    return pymol::make_error("invalid action ", action);
  }

  op.i1 = 1u << flag;
  op.i2 = ~op.i1;
  op.i3 = 0;
  op.i4 = 0;

  if ((op.i1 & cAtomFlag_exfoliate) && action != cFlagActionClear) {
    PRINTFB(G, FB_Executive, FB_Warnings)
    "The 'exfoliate' flag is deprecated. "
    "Use 'hide surface, (%s)' instead.\n",
        sele ENDFB(G);
  }

  SETUP_SELE(sele, s1, sele1);

  ExecutiveObjMolSeleOp(G, sele1, &op);
  if(Feedback(G, FB_Executive, FB_Actions)) {
    if(!quiet) {
      switch (action) {
      case cFlagActionReset:
        if(op.i3) {
          PRINTF " Flag: flag %d is set in %d of %d atoms.\n", flag, op.i3,
            op.i4 ENDF(G);
        } else {
          PRINTF " Flag: flag %d cleared on all atoms.\n", flag ENDF(G);
        }
        break;
      case cFlagActionSet:
        PRINTF " Flag: flag %d set on %d atoms.\n", flag, op.i3 ENDF(G);
        break;
      case cFlagActionClear:
        PRINTF " Flag: flag %d cleared on %d atoms.\n", flag, op.i3 ENDF(G);
        break;
      }
    }
  }
  if(SettingGetGlobal_b(G, cSetting_auto_indicate_flags)) {
    auto buffer = pymol::string_format("(flag %d)", flag);
    SelectorCreate(G, cIndicateSele, buffer.c_str(), NULL, true, NULL);
    ExecutiveSetObjVisib(G, cIndicateSele, true, false);
    SceneInvalidate(G);
  }
  return {};
}


/*========================================================================*/
float ExecutiveOverlap(PyMOLGlobals * G, const char *s1, int state1, const char *s2, int state2,
                       float adjust)
{
  SelectorTmp tmpsele1(G, s1);
  SelectorTmp tmpsele2(G, s2);
  int sele1 = tmpsele1.getIndex();
  int sele2 = tmpsele2.getIndex();

  float result = 0.0;

  if(state1 < 0)
    state1 = 0;
  if(state2 < 0)
    state2 = 0;

  if((sele1 >= 0) && (sele2 >= 0))
    result = SelectorSumVDWOverlap(G, sele1, state1, sele2, state2, adjust);

  return (result);
}


/*========================================================================*/
pymol::Result<> ExecutiveProtect(PyMOLGlobals * G, const char *s1, int mode, int quiet)
{
  ObjectMoleculeOpRec op;

  SETUP_SELE_DEFAULT(1);

  {
    ObjectMoleculeOpRecInit(&op);
    op.code = OMOP_Protect;
    op.i1 = mode;
    op.i2 = 0;
    ExecutiveObjMolSeleOp(G, sele1, &op);
    if(!quiet) {
      if(Feedback(G, FB_Executive, FB_Actions)) {
        if(op.i2) {
          if(mode) {
            PRINTF " Protect: %d atoms protected from movement.\n", op.i2 ENDF(G);
          } else {
            PRINTF " Protect: %d atoms deprotected.\n", op.i2 ENDF(G);
          }
        }
      }
    }
  }
  return {};
}


/*========================================================================*/
pymol::Result<> ExecutiveMask(
    PyMOLGlobals* G, const char* s1, int mode, int quiet)
{
  ObjectMoleculeOpRec op;

  SETUP_SELE_DEFAULT(1);

  {
    ObjectMoleculeOpRecInit(&op);
    op.code = OMOP_Mask;
    op.i1 = mode;
    op.i2 = 0;
    ExecutiveObjMolSeleOp(G, sele1, &op);
    if(!quiet) {
      if(Feedback(G, FB_Executive, FB_Actions)) {
        if(op.i2) {
          if(mode) {
            PRINTF " Mask: %d atoms masked (cannot be picked or selected).\n",
              op.i2 ENDF(G);
          } else {
            PRINTF " Mask: %d atoms unmasked.\n", op.i2 ENDF(G);
          }
        }
      }
    }
    op.code = OMOP_INVA;        /* need to invalidate all pickable representations */
    op.i1 = cRepsAtomMask;
    op.i2 = cRepInvPick;
    ExecutiveObjMolSeleOp(G, sele1, &op);
  }
  return {};
}


/*========================================================================*/
/**
 * flag > 0:  Set stereo_mode and turn stereo on
 * flag = 0:  Turn off stereo
 * flag = -1: Swap eyes (stereo_shift *= -1)
 * flag = -2: Turn on stereo with current stereo_mode
 * flag = -3: Turn on chromadepth and turn off stereo
 */
pymol::Result<> ExecutiveStereo(PyMOLGlobals * G, int flag)
{
  switch (flag) {
  case -3:
    SettingSet(G, cSetting_chromadepth, 1);
    SceneSetStereo(G, 0);
    break;
  case -1:
    SettingSetGlobal_f(G, cSetting_stereo_shift, -SettingGetGlobal_f(G, cSetting_stereo_shift));
    break;
  default:
    SettingSet(G, cSetting_chromadepth, 0);

    if (flag == cStereo_quadbuffer && !G->StereoCapable) {
      return pymol::make_error("no 'quadbuffer' support detected (force with 'pymol -S')");
    }

    if (flag == cStereo_openvr) {
#ifdef _PYMOL_OPENVR
      OpenVRInit(G);
      OpenVRFeedback(G);
#else
      return pymol::make_error("'openvr' stereo mode not available in this build");
#endif
    }

    if (flag > 0) {
      SettingSet(G, cSetting_stereo_mode, flag);
    }

    SceneSetStereo(G, flag != 0);
  }

  // for chromadepth
  G->ShaderMgr->Set_Reload_Bits(RELOAD_VARIABLES);

  SceneDirty(G);
  return {};
}


/*========================================================================*/
pymol::Result<> ExecutiveRevalence(PyMOLGlobals* G, const char* s1,
    const char* s2, const char* src, int target_state, int source_state,
    int reset, int quiet)
{
  SETUP_SELE_DEFAULT_PREFIXED(1, cSelectionInvalid);
  SETUP_SELE_DEFAULT_PREFIXED(2, sele1);

  {
    if(src && src[0]) {
      SETUP_SELE(src, tmpsele3, sele3);

      {
        ObjectMolecule *obj3 = SelectorGetSingleObjectMolecule(G, sele3);
        if(!obj3) {
          return pymol::make_error("Revalence can only source a single object at a time.");
        } else {
          ObjectMoleculeOpRec op;

          ObjectMoleculeOpRecInit(&op);
          op.code = OMOP_RevalenceFromSource;
          op.i1 = sele1;
          op.i2 = sele2;
          op.i3 = target_state;
          op.obj3 = obj3;
          op.i4 = sele3;
          op.i5 = source_state;
          op.i6 = quiet;

          ExecutiveObjMolSeleOp(G, sele1, &op);

          /*
             if(ObjectMoleculeXferValences(obj1,sele1,sele2,target_state,obj3,sele3,source_state,quiet)) {
             ObjectMoleculeVerifyChemistry(obj1,target_state);
             ObjectMoleculeInvalidate(obj1,cRepAll,cRepInvBonds,target_state);
             }
           */

        }
      }
    } else {                    /* guess valences */
      ObjectMoleculeOpRec op;

      ObjectMoleculeOpRecInit(&op);
      op.code = OMOP_RevalenceByGuessing;
      op.i1 = sele1;
      op.i2 = sele2;
      op.i3 = target_state;
      op.i4 = reset;
      op.i6 = quiet;

      ExecutiveObjMolSeleOp(G, sele1, &op);

    }
  }
  return {};
}


/*========================================================================*/
pymol::Result<> ExecutiveBond(PyMOLGlobals* G, const char* s1, const char* s2,
    int order, int mode, int quiet, pymol::zstring_view symop)
{
  int cnt;
  CExecutive *I = G->Executive;
  SpecRec *rec = NULL;
  int flag = false;

  SETUP_SELE_DEFAULT_PREFIXED(1, cSelectionInvalid);
  SETUP_SELE_DEFAULT_PREFIXED(2, cSelectionInvalid);

  {
    ObjectMolecule *obj1 = SelectorGetSingleObjectMolecule(G, sele1);
    ObjectMolecule *obj2 = SelectorGetSingleObjectMolecule(G, sele2);
    if((!obj1) || (!obj2) || (obj1 != obj2)) {
      if((!quiet) && (mode == 1)) {
        PRINTFB(G, FB_Editor, FB_Warnings)
          "Editor-Warning: bonds cannot be created between objects, only within.\n"
          ENDFB(G);
      }
    }
    while(ListIterate(I->Spec, rec, next)) {
      if(rec->type == cExecObject) {
        if(rec->obj->type == cObjectMolecule) {
          switch (mode) {
          case 1:              /* add */
            cnt = ObjectMoleculeAddBond((ObjectMolecule *) rec->obj, sele1, sele2, order, symop);
            if(cnt) {
              if(!quiet) {
                PRINTFB(G, FB_Editor, FB_Actions)
                  " Bond: %d bonds added to model \"%s\".\n", cnt, rec->obj->Name
                  ENDFB(G);
                flag = true;
              }
            }
            break;
          case 2:              /* adjust */
            cnt =
              ObjectMoleculeAdjustBonds((ObjectMolecule *) rec->obj, sele1, sele2, 1,
                                        order, symop);
            if(cnt) {
              if(!quiet) {
                PRINTFB(G, FB_Editor, FB_Actions)
                  " Valence: %d bond valences adjusted in model \"%s\".\n", cnt,
                  rec->obj->Name ENDFB(G);
                flag = true;
              }
            }
            break;
          case 0:              /* remove */
          default:
            cnt = ObjectMoleculeRemoveBonds((ObjectMolecule *) rec->obj, sele1, sele2);
            if(cnt) {
              if(!quiet) {
                PRINTFB(G, FB_Editor, FB_Actions)
                  " Unbond: %d bonds removed from model \"%s\".\n",
                  cnt, rec->obj->Name ENDFB(G);
              }
              flag = true;
            }
          }
        }
      }
    }
    if(!flag) {
      if(!quiet) {
        switch (mode) {
        case 1:
          PRINTFB(G, FB_Editor, FB_Warnings)
            "Bond-Warning: no bonds added." ENDFB(G);
          break;
        case 2:
          PRINTFB(G, FB_Editor, FB_Warnings)
            "Valence-Warning: no bond valences changed." ENDFB(G);
          break;
        case 0:
        default:
          PRINTFB(G, FB_Editor, FB_Warnings)
            "Unbond-Warning: no bonds removed." ENDFB(G);
          break;
        }
      }
    }
  }
  return {};
}

/**
 * brief Returns an error if quiet flag is not set
 * @param quiet quiet flag
 * @param msgs error messages
 * @return pymol::make_error with conditionally filled message
 */

template<typename... Ts>
static pymol::Error ErrorMsgIfQuiet(bool quiet, Ts&&... msgs){
  if(!quiet){
    return pymol::make_error(std::forward<Ts>(msgs)...);
  } else {
    return pymol::make_error("");
  }
}

/*========================================================================*/
pymol::Result<float> ExecutiveAngle(PyMOLGlobals* G,
    const char* nam, const char* s1, const char* s2, const char* s3, int mode,
    int labels, int reset, int zoom, int quiet, int state, int state1,
    int state2, int state3)
{
  SETUP_SELE_DEFAULT_PREFIXED(1, cSelectionInvalid);
  SETUP_SELE_DEFAULT_PREFIXED(2, sele1);
  SETUP_SELE_DEFAULT_PREFIXED(3, sele2);

  auto obj = ExecutiveFindOrDeleteObject<ObjectDist>(G, nam);
  auto need_manage = !obj;
  float result = -1.0F;

  obj = ObjectDistNewFromAngleSele(G, obj,
                                     sele1, sele2, sele3,
                                     mode, labels, &result, reset, state,
                                     state1, state2, state3);

  // ObjectDistNewFrom... always succeeds
  assert(obj);

  if (need_manage) {
        ObjectSetName(obj, nam);
        ExecutiveManageObject(G, obj, zoom, quiet);
        if(!labels)
          ExecutiveSetRepVisib(G, nam, cRepLabel, 0);
  }

  return rad_to_deg(result);
}


/*========================================================================*/
pymol::Result<float> ExecutiveDihedral(PyMOLGlobals* G, const char* nam,
    const char* s1, const char* s2, const char* s3, const char* s4, int mode,
    int labels, int reset, int zoom, int quiet, int state)
{
  SETUP_SELE_DEFAULT_PREFIXED(1, cSelectionInvalid);
  SETUP_SELE_DEFAULT_PREFIXED(2, sele1);
  SETUP_SELE_DEFAULT_PREFIXED(3, sele2);
  SETUP_SELE_DEFAULT_PREFIXED(4, sele3);

  auto obj = ExecutiveFindOrDeleteObject<ObjectDist>(G, nam);
  auto need_manage = !obj;
  float result = -1.0F;

  obj = ObjectDistNewFromDihedralSele(G, obj,
                                        sele1, sele2, sele3, sele4,
                                        mode, labels, &result, reset, state);

  // ObjectDistNewFrom... always succeeds
  assert(obj);

  if (need_manage) {
        ObjectSetName(obj, nam);
        ExecutiveManageObject(G, obj, zoom, quiet);
        if(!labels)
          ExecutiveSetRepVisib(G, nam, cRepLabel, 0);
  }

  return rad_to_deg(result);
}

/**
 * Create a distance measurement object
 *
 * @param nam: name of measurement object to create or add to
 * @param s1: selection expression
 * @param s2: selection expression or "same" keyword (shortcut for s1 = s2)
 * @param mode: 0 (any), 1 (bonds), 2 (hbonds), 3 (distance_exclusion),
 * 4 (centroids), 5 (pi-*), 6 (pi-pi), 7 (pi-cat), 8 (vdw-dist-ratio)
 * @param cutoff: Distance cutoff in Angstrom, or vdw-distance ratio cutoff (mode 8)
 * @param labels: Boolean flag to show distance labels
 * @param reset: Boolean flag to clear a distance object if it already exists
 * @return average of measured distances in Angstrom
 */
pymol::Result<float> ExecutiveDistance(PyMOLGlobals* G, const char* nam,
    const char* s1, const char* s2, int mode, float cutoff, int labels,
    int quiet, int reset, int state, int zoom, int state1, int state2)
{
  if (strcmp(s1, s2) == 0) {
    s2 = cKeywordSame;
  }

  SETUP_SELE_DEFAULT_PREFIXED(1, cSelectionInvalid);
  SETUP_SELE_DEFAULT_PREFIXED(2, sele1);

  auto obj = ExecutiveFindOrDeleteObject<ObjectDist>(G, nam);
  auto need_manage = !obj;
  float result = -1.0F;

  /* create a new distance from the two selections */
  obj = ObjectDistNewFromSele(G, obj, sele1, sele2, mode, cutoff, labels, reset,
      &result, state, state1, state2);

  // ObjectDistNewFrom... always succeeds
  assert(obj);

  if (need_manage) {
    switch (mode) {
    case 6: // pi-pi
      SettingSet(cSetting_dash_color, "0x06c5ff" /* light blue */, obj);
      break;
    case 7: // pi-cat
      SettingSet(cSetting_dash_color, "0x468c00" /* dark green */, obj);
      break;
    case 8: // "clashes"
      SettingSet(cSetting_dash_color, "0xff8800" /* light red */, obj);
      break;
    }

      ObjectSetName(obj, nam);
      ExecutiveManageObject(G, obj, zoom, quiet);
      if(!labels)
        ExecutiveSetRepVisib(G, nam, cRepLabel, 0);
  }

  return result;
}

/*========================================================================*/
char *ExecutiveNameToSeqAlignStrVLA(PyMOLGlobals * G, const char *name, int state, int format,
                                    int quiet)
{
  char *result = NULL;
  if((!name) || (!name[0]) || (strcmp(name, "(all)") == 0)) {
    /* use current alignment as the default */
    name = SettingGetGlobal_s(G, cSetting_seq_view_alignment);
    if(name[0] == 0) {
      SpecRec *rec = NULL;
      CExecutive *I = G->Executive;
      while(ListIterate(I->Spec, rec, next)) {
        if(rec->visible) {
          if(rec->type == cExecObject)
            if(rec->obj->type == cObjectAlignment) {
              name = rec->obj->Name;
              break;
            }
        }
      }
    }
  }
  if(!name) {
    ErrMessage(G, " Executive", "invalid alignment object name.");
  } else {
    pymol::CObject *obj = ExecutiveFindObjectByName(G, name);

    if(!obj) {
      ErrMessage(G, " Executive", "alignment object not found.");
    } else if(obj->type != cObjectAlignment) {
      ErrMessage(G, " Executive", "invalid object type.");
    } else {
      ObjectAlignmentAsStrVLA(G, (ObjectAlignment *) obj, state, format, &result);
    }
  }
  return (result);
}


/*========================================================================*/
pymol::Result<> ExecutiveSeleToObject(PyMOLGlobals* G, const char* name,
    const char* s1, int source, int target, int discrete, int zoom, int quiet,
    int singletons, int copy_properties)
{
  SelectorTmp tmpsele1(G, s1);
  int sele1 = tmpsele1.getIndex();

  int ok = false;
  ObjectNameType valid_name;

  UtilNCopy(valid_name, name, sizeof(valid_name));
  if(SettingGetGlobal_b(G, cSetting_validate_object_names)) {
    ObjectMakeValidName(G, valid_name);
    name = valid_name;
  }
  {
    int exists = (ExecutiveFindObjectMoleculeByName(G, name) != NULL);

    if(sele1 >= 0) {
      ok = SelectorCreateObjectMolecule(G, sele1, name, target,
                                        source, discrete, false, quiet, singletons, copy_properties);
      if(ok) {
        int sele2 = SelectorIndexByName(G, name);
        ObjectMolecule *old_obj, *new_obj;
        old_obj = SelectorGetFirstObjectMolecule(G, sele1);     /* get at least one object */
        new_obj = SelectorGetSingleObjectMolecule(G, sele2);

        /* first we need to make sure that the object being moved
           matches the target with respect to both the TTT and the
           object's state matrix (if any) */

        if(old_obj && new_obj) {
          ExecutiveMatrixCopy(G, old_obj->Name, new_obj->Name, 1, 1,    /* TTT mode */
                              source, target, false, 0, quiet);

          ExecutiveMatrixCopy(G, old_obj->Name, new_obj->Name, 2, 2,    /* Object state mode */
                              source, target, false, 0, quiet);

#ifdef _PYMOL_IP_PROPERTIES
#endif

          ExecutiveDoZoom(G, new_obj, !exists, zoom, true);
        }
      }
    }
  }
  if(ok){
    return {};
  } else {
    return pymol::make_error("Failed to Create Object");
  }
}


/*========================================================================*/
pymol::Result<> ExecutiveCopy(
    PyMOLGlobals* G, const char* src, const char* dst, int zoom)
{
  const pymol::CObject *os = ExecutiveFindObjectByName(G, src);
  if(!os) {
    return pymol::make_error("Object not found.");
  }

  pymol::CObject* oDst = os->clone();

  if(!oDst) {
    return pymol::make_error("Failed to create copy");
  }

  strcpy(oDst->Name, dst);
  ExecutiveManageObject(G, oDst, zoom, false);

  PRINTFB(G, FB_Executive, FB_Actions)
    " Executive: object %s created.\n", oDst->Name ENDFB(G);
  SceneChanged(G);
  return {};
}


/*========================================================================*/
pymol::Result<> ExecutiveOrient(PyMOLGlobals * G, const char *sele,
                     int state, float animate, int complete, float buffer, int quiet)
{
  double egval[3], egvali[3];
  double evect[3][3];
  float m[4][4], mt[4][4];
  float t[3];
  const float _0 = 0.0F;
  int a, b;

  double mi[16];

  SelectorTmp tmpsele(G, sele);
  sele = tmpsele.getName();

  if(!ExecutiveGetMoment(G, sele, mi, state)) {
    return {};
  }

  if(!MatrixEigensolveC33d(G, mi, egval, egvali, (double *) (void *) evect)) {

    normalize3d(evect[0]);
    normalize3d(evect[1]);
    normalize3d(evect[2]);

    for(a = 0; a < 3; a++) {
      for(b = 0; b < 3; b++) {
        m[a][b] = (float) evect[b][a];  /* fill columns */
      }
    }

    for(a = 0; a < 3; a++) {    /* expand to 4x4 */
      m[3][a] = 0;
      m[a][3] = 0;
    }
    m[3][3] = 1.0;

    normalize3f(m[0]);          /* cross normalization (probably unnec.)  */
    normalize3f(m[1]);
    normalize3f(m[2]);

    for(a = 0; a < 3; a++)      /* convert to row-major */
      for(b = 0; b < 3; b++)
        mt[a][b] = m[b][a];

    cross_product3f(mt[0], mt[1], t);   /* insure right-handed matrix */
    if(dot_product3f(t, mt[2]) < 0.0) {
      mt[2][0] = -mt[2][0];
      mt[2][1] = -mt[2][1];
      mt[2][2] = -mt[2][2];
    }

    for(a = 0; a < 3; a++)      /* convert back to column major */
      for(b = 0; b < 3; b++)
        m[a][b] = mt[b][a];

    if(animate < 0.0F) {
      if(SettingGetGlobal_b(G, cSetting_animation))
        animate = SettingGetGlobal_f(G, cSetting_animation_duration);
      else
        animate = 0.0F;
    }
    if(animate != 0.0F)
      ScenePrimeAnimation(G);

    {
      float old_mat[16];
      float new_mat[16];
      float x, y, z;
      copy44f(SceneGetMatrix(G), old_mat);

      SceneSetMatrix(G, m[0]);  /* load matrix */

      /* there must  be a more elegant to get the PC on X and the SC
       * on Y then what is shown below, but I couldn't get it to work.
       * I tried swapping the eigen-columns around but either that is 
       * a bogus approach (?) or my code was buggy.  Hence the following...*/

      if((egval[0] < egval[2]) && (egval[2] < egval[1])) {      /* X < Z < Y */
        SceneRotate(G, 90, 1, 0, 0);    /*1<-->2 */
      } else if((egval[1] < egval[0]) && (egval[0] < egval[2])) {       /* Y < X < Z */
        SceneRotate(G, 90, 0, 0, 1);    /*0<-->1 */
      } else if((egval[1] < egval[2]) && (egval[2] < egval[0])) {       /* Y < Z < X */
        SceneRotate(G, 90, 0, 1, 0);    /*1<-->2 */
        SceneRotate(G, 90, 0, 0, 1);    /*0<-->1 */
      } else if((egval[2] < egval[1]) && (egval[1] < egval[0])) {       /* Z < Y < X */
        SceneRotate(G, 90, 0, 1, 0);    /*0<-->2 */
      } else if((egval[2] < egval[0]) && (egval[0] < egval[1])) {       /* Z < X < Y */
        SceneRotate(G, 90, 0, 1, 0);    /*0<-->2 */
        SceneRotate(G, 90, 1, 0, 0);    /*0<-->1 */
      }

      /* now choose orientation that has the least perturbation from the starting matrix */

      copy44f(SceneGetMatrix(G), new_mat);

      x = old_mat[0] * new_mat[0] + old_mat[4] * new_mat[4] + old_mat[8] * new_mat[8];
      y = old_mat[1] * new_mat[1] + old_mat[5] * new_mat[5] + old_mat[9] * new_mat[9];
      z = old_mat[2] * new_mat[2] + old_mat[6] * new_mat[6] + old_mat[10] * new_mat[10];

      if((x > _0) && (y < _0) && (z < _0)) {
        SceneRotate(G, 180, 1, 0, 0);
      } else if((x < _0) && (y > _0) && (z < _0)) {
        SceneRotate(G, 180, 0, 1, 0);
      } else if((x < _0) && (y < _0) && (z > _0)) {
        SceneRotate(G, 180, 0, 0, 1);
      }
    }

    /* X < Y < Z  - do nothing - that's what we want */

    ExecutiveWindowZoom(G, sele, buffer, state, complete, false, quiet);
    if(animate != 0.0F)
      SceneLoadAnimation(G, animate, 0);

  }
  return {};
}

pymol::Result<> ExecutiveMove(
    PyMOLGlobals* G, pymol::zstring_view axis, float dist)
{
  switch (axis[0]) {
  case 'x':
    SceneTranslate(G, dist, 0.0, 0.0);
    break;
  case 'y':
    SceneTranslate(G, 0.0, dist, 0.0);
    break;
  case 'z':
    SceneTranslate(G, 0.0, 0.0, dist);
    break;
  default:
    return pymol::make_error("Axis must be x, y, or z");
  }
  return {};
}

/*========================================================================*/
pymol::Result<> ExecutiveLabel(PyMOLGlobals * G, const char *str1, const char *expr, int quiet, int eval_mode)
{
  int sele1;
  ObjectMoleculeOpRec op1;
  int cnt;

  SelectorTmp s1(G, str1);
  sele1 = s1.getIndex();
  if(sele1 >= 0) {
    ObjectMoleculeOpRecInit(&op1);
    op1.code = OMOP_LABL;
    op1.s1 = expr;
    op1.i1 = 0;
    op1.i2 = eval_mode;

#ifndef _PYMOL_NOPY
    pymol::pautoblock gil(G);
#endif

    if (!ExecutiveObjMolSeleOp(G, sele1, &op1)) {
      return pymol::Error();
    }

    cnt = op1.i1;
    op1.code = OMOP_VISI;
    op1.i1 = cRepLabelBit;
    op1.i2 = cVis_SHOW;
    ExecutiveObjMolSeleOp(G, sele1, &op1);
    op1.code = OMOP_INVA;
    op1.i2 = cRepInvVisib;
    ExecutiveObjMolSeleOp(G, sele1, &op1);

    if(!quiet) {
      {
	const char *unlabelledstr = "";
	if (cnt<0){ /* if negative, say unlabelled */
	  cnt = -cnt;
	  unlabelledstr = "un";
	}
      PRINTFB(G, FB_Executive, FB_Actions)
        " Label: %slabelled %i atoms.\n", unlabelledstr, cnt ENDFB(G);
      }
    }
  } else {
    return pymol::make_error("No atoms selected");
  }
  return {};
}


/*========================================================================*/
#ifdef _WEBGL
#else
pymol::Result<int> ExecutiveIterate(PyMOLGlobals * G, const char *str1, const char *expr, int read_only, int quiet,
                     PyObject * space)
#endif
{
  ObjectMoleculeOpRec op1;
  ObjectMoleculeOpRecInit(&op1);
#ifdef _WEBGL
#endif
  SelectorTmp tmpsele1(G, str1);
  int sele1 = tmpsele1.getIndex();
  op1.i1 = 0;
  if(sele1 >= 0) {
    op1.code = OMOP_ALTR;
    op1.i1 = 0;
    op1.i2 = read_only;
#ifdef _WEBGL
#else
    op1.s1 = expr;
    op1.py_ob1 = space;
#endif

    if (!ExecutiveObjMolSeleOp(G, sele1, &op1)) {
      return pymol::Error();
    }

    if(!quiet) {
      if(!read_only) {
        PRINTFB(G, FB_Executive, FB_Actions)
          " Alter: modified %i atoms.\n", op1.i1 ENDFB(G);
      } else {
        PRINTFB(G, FB_Executive, FB_Actions)
          " Iterate: iterated over %i atoms.\n", op1.i1 ENDFB(G);
      }
    }
    if (!read_only) {
      SeqChanged(G);
    }
  } else {
    if(!quiet) {
      PRINTFB(G, FB_Executive, FB_Warnings)
        " %s: No atoms selected.\n", __func__ ENDFB(G);
    }
  }
  return (op1.i1);
}

SelectArgs ExecutiveSelectPrepareArgs(PyMOLGlobals* G, pymol::zstring_view sname, pymol::zstring_view sele)
{
  SelectArgs args;
  args.sname = sname.c_str();
  args.sele = sele.c_str();

  // selection in first argument (cmd.select("expr"))
  if (args.sele.empty()) {
    args.sele = sname.c_str();
    args.sname = SettingGet<bool>(G, cSetting_auto_number_selections) ? "" : "sele";
  }

  if (args.sname.empty()) {
    auto sel_num = SettingGet<int>(G, cSetting_sel_counter) + 1;
    SettingSet<int>(G, cSetting_sel_counter, sel_num);
    args.sname = pymol::string_format("sel%02u", static_cast<unsigned>(sel_num));
  }
  return args;
}

/**
 * cmd.select() implementation
 *
 * @param name Name of selection to create, or selection expression if `sele` is
 * empty. If `sele` is empty, `name` will become "sele" or "sel01", depending on
 * `auto_number_selections` and `sel_counter`. If `name` is empty, it will
 * become "sel01", depending on `sel_counter` but independent of
 * `auto_number_selections`.
 * @param sele Selection expression
 * @param enable Enable the named selection's SpecRec if 1, disable if 0, keep
 * unchanged (if exists, otherwise enable) if -1.
 * @param merge Discard existing named selection if 0. Merge with existing named
 * selection if 1, merge only if existing one is enabled if 2.
 * @param state Object state for state (coordinate) dependent expressions
 * @param domain Existing selection name, same as selecting `(sele) and (domain)`
 * @return Number of selected atoms
 */
pymol::Result<int> ExecutiveSelect(PyMOLGlobals* G, const SelectArgs& sargs,
    int enable, int quiet, int merge, int state,
    const char* domain)
{
  const char* name = sargs.sname.c_str();
  const char* sele = sargs.sele.c_str();
  // bail if name not available
  if (ExecutiveFindObjectByName(G, name)) {
    return pymol::make_error("name conflicts with an object");
  }

  // merge with existing selection
  std::string selebuf;
  if (merge) {
    if (merge == 2) {
      // merge if exists and active
      selebuf = pymol::join_to_string("(", sele, ") or ??", name);
    } else {
      // merge if exists
      selebuf = pymol::join_to_string("(", sele, ") or ?", name);
    }
    sele = selebuf.c_str();
  }

  auto res = SelectorCreateWithStateDomain(
      G, name, sele, nullptr, quiet, nullptr, state, domain);

  p_return_if_error(res);

  if (enable == 1) {
    ExecutiveSetObjVisib(G, name, 1, 0);
  } else if (enable == 0) {
    ExecutiveSetObjVisib(G, name, 0, 0);
  }

  SceneInvalidate(G);
  SeqDirty(G);

  return res.result();
}

/*========================================================================*/
pymol::Result<int> ExecutiveSelectList(PyMOLGlobals* G, const char* sele_name,
    const char* oname, const int* list, size_t list_len, int state, int mode,
    int quiet)
{
  enum {
    MODE_INDEX = 0,
    MODE_ID = 1,
    MODE_RANK = 2,
  };

  auto obj = ExecutiveFindObject<ObjectMolecule>(G, oname);
  if (!obj) {
    return pymol::Error("object not found");
  }

  std::vector<int> idx_list;
  idx_list.reserve(list_len);

  if (mode == MODE_INDEX) {
    for (size_t i = 0; i != list_len; ++i) {
      // convert 1-based indices to 0-based array offsets
      idx_list.push_back(list[i] - 1);
    }
  } else if (mode == MODE_ID || mode == MODE_RANK) {
    CoordSet const* cs = obj->getCoordSet(state);
    std::set<int> const list_set(list, list + list_len);

    for (int atm = 0; atm < obj->NAtom; ++atm) {
      auto const& ai = obj->AtomInfo[atm];
      int const index = (mode == MODE_ID) ? ai.id : ai.rank;

      if (list_set.count(index)) {
        if (!cs || cs->atmToIdx(atm) >= 0) {
          idx_list.push_back(atm);
        }
      }
    }
  } else {
    return pymol::Error("invalid mode");
  }

  return SelectorCreateOrderedFromObjectIndices(
      G, sele_name, obj, idx_list.data(), idx_list.size());
}


/*========================================================================*/
pymol::Result<int> ExecutiveIterateList(PyMOLGlobals* G, const char* str1,
    PyObject* list, int read_only, int quiet, PyObject* space)
{
#ifdef _PYMOL_NOPY
  return pymol::make_error("Iterate List not available.");
#else
  assert(PyGILState_Check());

  int ok = true;
  int n_eval = 0;
  SelectorTmp s1(G, str1);
  int sele0 = s1.getIndex();
  PyObject *entry = NULL;
  ObjectMolecule *obj = NULL;
  if(sele0 >= 0)
    obj = SelectorGetSingleObjectMolecule(G, sele0);
  if(obj) {
    int n_atom = obj->NAtom;
    int list_len = 0;
    int a;
    int index = 0;
    const char *expr = NULL;
    if(ok)
      ok = PyList_Check(list);
    if(ok) {
      list_len = PyList_Size(list);
      for(a = 0; a < list_len; a++) {
        if(ok)
          entry = PyList_GetItem(list, a);
        if(ok)
          ok = PyList_Check(entry);
        if(ok)
          ok = (PyList_Size(entry) == 2);
        if(ok)
          ok = PConvPyIntToInt(PyList_GetItem(entry, 0), &index);
        if(ok)
          ok = PConvPyStrToStrPtr(PyList_GetItem(entry, 1), &expr);
        if(ok)
          ok = ((index <= n_atom) && (index > 0));
        if(ok)
        {
	  CoordSet *cs = NULL;
	  if(obj->DiscreteFlag && obj->DiscreteCSet) {
	    cs = obj->DiscreteCSet[index - 1];
	  } else if (obj->NCSet == 1){
	    cs = obj->CSet[0];
	  }
          auto expr_co =
              unique_PyObject_ptr(Py_CompileString(expr, "", Py_single_input));
          ok = expr_co && PAlterAtom(G, obj, cs, expr_co.get(), read_only,
                              index - 1, space);
        }
        if(ok)
          n_eval++;
        else
          break;
      }
    }
  } else {
    return pymol::make_error("Selection cannot span more than one object.");
  }
  if(ok) {
    if(!quiet) {
      if(!read_only) {
        PRINTFB(G, FB_Executive, FB_Actions)
          " AlterList: modified %i atoms.\n", n_eval ENDFB(G);
      } else {
        PRINTFB(G, FB_Executive, FB_Actions)
          " IterateList: iterated over %i atoms.\n", n_eval ENDFB(G);
      }
    }
    if (!read_only) {
      SeqChanged(G);
    }
  }
  if(!ok)
    return pymol::make_error("An error occurred.");
  else
    return n_eval;
#endif
}


/*========================================================================*/
#ifdef _WEBGL
#else
pymol::Result<int> ExecutiveIterateState(PyMOLGlobals* G, int state,
    const char* str1, const char* expr, int read_only,
    int quiet, PyObject* space)
#endif
{
#ifdef _WEBGL
#endif
  SelectorTmp tmpsele1(G, str1);
  int sele1 = tmpsele1.getIndex();
  if(sele1 >= 0) {
    int start_state = 0, stop_state = 0;
    ObjectMoleculeOpRec op1;

    if(state >= 0) {
      start_state = state;
      stop_state = state + 1;
    } else {
      if((state == -2) || (state == -3)) {      /* current state, TO DO: effective object state */
        state = SceneGetState(G);
        start_state = state;
        stop_state = state + 1;
      } else if(state == -1) {  /* all states (for the selection) */
        start_state = 0;
        stop_state = SelectorCountStates(G, sele1);
      }
    }
    ObjectMoleculeOpRecInit(&op1);
    op1.i1 = 0;

    for(state = start_state; state < stop_state; state++) {
      op1.code = OMOP_AlterState;
#ifdef _WEBGL
#else
      op1.s1 = expr;
      op1.py_ob1 = space;
#endif
      op1.i2 = state;
      op1.i3 = read_only;
      if (!ExecutiveObjMolSeleOp(G, sele1, &op1)) {
        return pymol::Error();
      }
    }
    if(!read_only) {
      // for dynamic_measures
      ExecutiveUpdateCoordDepends(G, NULL);
      SeqChanged(G);
    }
    if(!quiet) {
      if(!read_only) {
        PRINTFB(G, FB_Executive, FB_Actions)
          " AlterState: modified %i atom coordinate states.\n", op1.i1 ENDFB(G);
      } else {
        PRINTFB(G, FB_Executive, FB_Actions)
          " IterateState: iterated over %i atom coordinate states.\n", op1.i1 ENDFB(G);
      }
    }
    return op1.i1;
  } else {
    if(!quiet) {
      PRINTFB(G, FB_Executive, FB_Warnings)
        "ExecutiveIterateState: No atoms selected.\n" ENDFB(G);
    }
    return 0;
  }
}

typedef struct {
  int priority;
  float vertex[3];
  AtomInfoType *ai;
} FitVertexRec;

static int fVertexOrdered(const FitVertexRec * array, int l, int r)
{
  return (array[l].priority <= array[r].priority);
}

static int fAtomOrdered(PyMOLGlobals * G, const AtomInfoType* const* array, int l, int r)
{
  return (AtomInfoCompare(G, array[l], array[r]));
}

static int fAtomIDOrdered(const AtomInfoType* const* array, int l, int r)
{
  return (array[l]->id <= array[r]->id);
}

static int fAtomRankOrdered(const AtomInfoType* const* array, int l, int r)
{
  return (array[l]->rank <= array[r]->rank);
}

static int fAtomTemp1Ordered(const AtomInfoType* const* array, int l, int r)
{
  return (array[l]->temp1 <= array[r]->temp1);
}

static void PackSortedIndices(int n, int *x, int rec_size, void *data)
{
  int a;
  for(a = 0; a < n; a++) {
    if(a != x[a]) {
      memcpy(((char *) data) + (a * rec_size),
             ((char *) data) + (x[a] * rec_size), rec_size);
    }
  }
}


/*========================================================================*/
int ExecutiveRMS(PyMOLGlobals * G, const char *s1, const char *s2, int mode, float refine,
                 int max_cyc, int quiet, const char *oname, int state1, int state2,
                 int ordered_selections, int matchmaker, ExecutiveRMSInfo * rms_info)
{
  /* mode 0: measure rms without superposition
     mode 1: measure rms with trial superposition
     mode 2: measure rms with actual superposition */

  int sele1, sele2;
  float rms = -1.0;
  int a, b;
  float inv, *f, *f1, *f2;
  ObjectMoleculeOpRec op1;
  ObjectMoleculeOpRec op2;
  OrthoLineType buffer;
  int *flag;
  int ok = true;
  int repeat;
  float v1[3], *v2;
  ObjectAlignment *align_to_update = NULL;

  bool ignore_case = SettingGetGlobal_b(G, cSetting_ignore_case);
  bool ignore_case_chain = SettingGetGlobal_b(G, cSetting_ignore_case_chain);

  int matrix_mode = SettingGetGlobal_i(G, cSetting_matrix_mode);
  if(matrix_mode < 0)
    matrix_mode = 0; /* for now */

  if(matchmaker == -1) {
    /* matchmaker -1 is the same as matchmaker 0 except that the
       selections are not pre-matched prior to calling of this routine */
    matchmaker = 0;
  }

  sele1 = SelectorIndexByName(G, s1);
  sele2 = SelectorIndexByName(G, s2);

  /* this function operates on stored coordinates -- thus transformation 
     matrices will need to be applied to the resulting atoms */

  // get coordinates
  {
    auto sele = sele1;
    auto state = state1;
    auto op = &op1;

    // for both selections
    do {
      ObjectMoleculeOpRecInit(op);

      if(sele >= 0) {
        if(state < 0) {
          op->code = OMOP_AVRT;
        } else {
          op->code = OMOP_StateVRT;
          op->i1 = state;
        }

        op->nvv1 = 0;                       // length of vc1 (number of atoms with coordinates)
        op->vc1 = VLACalloc(int, 1000);     // number of states per atom
        op->vv1 = VLACalloc(float, 1000);   // coordinates (sum over states)

        if(mode == 0)
          op->i2 = true;            /* if measuring current coordinates, then get global txfd values */

        if(matchmaker || (oname && oname[0]))
          op->ai1VLA = VLACalloc(AtomInfoType*, 1000);

        if(ordered_selections)
          op->vp1 = VLAlloc(int, 1000);     // selection member "priority"? (MemberType::tag)

        ExecutiveObjMolSeleOp(G, sele, op);

        for(a = 0; a < op->nvv1; a++) {
          inv = (float) op->vc1[a]; /* average over coordinate sets */
          if(inv) {
            f = op->vv1 + (a * 3);
            scale3f(f, 1.F / inv, f);
          }
        }
      }

      // second iteration done
      if (sele == sele2)
        break;

      sele = sele2;
      state = state2;
      op = &op2;

    } while (true);
  }

  if(op1.vv1 && op2.vv1) {
    if(op1.nvv1 && op2.nvv1) {
      ObjectMolecule *mobile_obj = NULL;

      int n_pair = 0;

      if(!(mobile_obj = SelectorGetSingleObjectMolecule(G, sele1))) {
        if(mode != 2) {
          PRINTFB(G, FB_Executive, FB_Warnings)
            "Executive-Warning: Mobile selection spans more than one object.\n" ENDFB(G);
        } else {
          PRINTFB(G, FB_Executive, FB_Errors)
            "Executive-Error: Mobile selection spans more than one object. Aborting.\n"
            ENDFB(G);
          ok = false;
        }
      }

      if(ok && op1.nvv1 && op2.nvv1 && (matchmaker > 0)) {
        /* matchmaker 0 is the default... by internal atom ordering only */
        int *idx1 = pymol::malloc<int>(op1.nvv1);
        int *idx2 = pymol::malloc<int>(op2.nvv1);
        int sort_flag = false;
        if(!(idx1 && idx2))
          ok = false;
        else {
          switch (matchmaker) {
          case 1:              /* by atom info-based ordering */
            UtilSortIndexGlobals(G, op1.nvv1, op1.ai1VLA, idx1,
                                 (UtilOrderFnGlobals *) fAtomOrdered);
            UtilSortIndexGlobals(G, op2.nvv1, op2.ai1VLA, idx2,
                                 (UtilOrderFnGlobals *) fAtomOrdered);
            sort_flag = true;
            break;
          case 2:              /* by matching atom identifiers */
            UtilSortIndex(op1.nvv1, op1.ai1VLA, idx1, (UtilOrderFn *) fAtomIDOrdered);
            UtilSortIndex(op2.nvv1, op2.ai1VLA, idx2, (UtilOrderFn *) fAtomIDOrdered);
            sort_flag = true;
            break;
          case 3:              /* by matching atom ranks */
            UtilSortIndex(op1.nvv1, op1.ai1VLA, idx1, (UtilOrderFn *) fAtomRankOrdered);
            UtilSortIndex(op2.nvv1, op2.ai1VLA, idx2, (UtilOrderFn *) fAtomRankOrdered);
            sort_flag = true;
            break;
          case 4:              /* by internal atom indexes (stored in temp1 kludge field) */
            UtilSortIndex(op1.nvv1, op1.ai1VLA, idx1, (UtilOrderFn *) fAtomTemp1Ordered);
            UtilSortIndex(op2.nvv1, op2.ai1VLA, idx2, (UtilOrderFn *) fAtomTemp1Ordered);
            sort_flag = true;
            break;
          }
          if(sort_flag) {
            /* GOD this is SO ugly! */

            if(op1.vv1) {
              float *tmp = VLAlloc(float, op1.nvv1 * 3);
              if(!tmp)
                ok = false;
              else {
                UtilApplySortedIndices(op1.nvv1, idx1, 3 * sizeof(float), op1.vv1, tmp);
                VLAFreeP(op1.vv1);
                op1.vv1 = tmp;
              }
            }
            if(op1.vc1) {
              int *tmp = VLAlloc(int, op1.nvv1);
              if(!tmp)
                ok = false;
              else {
                UtilApplySortedIndices(op1.nvv1, idx1, sizeof(int), op1.vc1, tmp);
                VLAFreeP(op1.vc1);
                op1.vc1 = tmp;
              }
            }
            if(op1.vp1) {
              int *tmp = VLAlloc(int, op1.nvv1);
              if(!tmp)
                ok = false;
              else {
                UtilApplySortedIndices(op1.nvv1, idx1, sizeof(int), op1.vp1, tmp);
                VLAFreeP(op1.vp1);
                op1.vp1 = tmp;
              }
            }
            if(op1.ai1VLA) {
              AtomInfoType **tmp = VLACalloc(AtomInfoType *, op1.nvv1);
              if(!tmp)
                ok = false;
              else {
                UtilApplySortedIndices(op1.nvv1, idx1, sizeof(AtomInfoType *), op1.ai1VLA,
                                       tmp);
                VLAFreeP(op1.ai1VLA);
                op1.ai1VLA = tmp;
              }
            }

            if(op2.vv1) {
              float *tmp = VLAlloc(float, op2.nvv1 * 3);
              if(!tmp)
                ok = false;
              else {
                UtilApplySortedIndices(op2.nvv1, idx2, 3 * sizeof(float), op2.vv1, tmp);
                VLAFreeP(op2.vv1);
                op2.vv1 = tmp;
              }
            }
            if(op2.vc1) {
              int *tmp = VLAlloc(int, op2.nvv1);
              if(!tmp)
                ok = false;
              else {
                UtilApplySortedIndices(op2.nvv1, idx2, sizeof(int), op2.vc1, tmp);
                VLAFreeP(op2.vc1);
                op2.vc1 = tmp;
              }
            }
            if(op2.vp1) {
              int *tmp = VLAlloc(int, op2.nvv1);
              if(!tmp)
                ok = false;
              else {
                UtilApplySortedIndices(op2.nvv1, idx2, sizeof(int), op2.vp1, tmp);
                VLAFreeP(op2.vp1);
                op2.vp1 = tmp;
              }
            }
            if(op2.ai1VLA) {
              AtomInfoType **tmp = VLACalloc(AtomInfoType *, op2.nvv1);
              if(!tmp)
                ok = false;
              else {
                UtilApplySortedIndices(op2.nvv1, idx2, sizeof(AtomInfoType *), op2.ai1VLA,
                                       tmp);
                VLAFreeP(op2.ai1VLA);
                op2.ai1VLA = tmp;
              }
            }

          }
        }

        if(matchmaker != 0) {
          int n1 = 0, n2 = 0, c1 = 0, c2 = 0;
          int cmp;

          while((n1 < op1.nvv1) && (n2 < op2.nvv1)) {
            cmp = 0;
            switch (matchmaker) {
            case 1:            /* insure that AtomInfoType matches */
              if(AtomInfoMatch(G, op1.ai1VLA[n1], op2.ai1VLA[n2], ignore_case, ignore_case_chain))
                cmp = 0;
              else
                cmp = AtomInfoCompare(G, op1.ai1VLA[n1], op2.ai1VLA[n2]);
              printf("%d-%d %d-%d: %d\n", c1, n1, c2, n2, cmp);
              break;
            case 2:            /* ID */
            case 3:            /* rank */
              {
                int val1;
                int val2;

                switch (matchmaker) {
                case 2:        /* ID */
                  val1 = op1.ai1VLA[n1]->id;
                  val2 = op2.ai1VLA[n2]->id;
                  break;
                case 3:        /* rank */
                  val1 = op1.ai1VLA[n1]->rank;
                  val2 = op2.ai1VLA[n2]->rank;
                  break;
                case 4:        /* index (via temp1) */
                  val1 = op1.ai1VLA[n1]->temp1;
                  val2 = op2.ai1VLA[n2]->temp1;
                  break;
                default:
                  val1 = 0;
                  val2 = 0;
                  break;
                }
                if(val1 == val2)
                  cmp = 0;
                else if(val1 < val2)
                  cmp = -1;
                else
                  cmp = 1;
              }
              break;
            }
            if(!cmp) {          /* match found */
              idx1[c1++] = n1++;
              idx2[c2++] = n2++;
              n_pair++;
            } else if(cmp < 0) {        /* op1 below op2 */
              n1++;
            } else {            /* op2 below op1 */
              n2++;
            }
          }

          if(n_pair) {
            if(op1.vv1)
              PackSortedIndices(n_pair, idx1, 3 * sizeof(float), op1.vv1);
            if(op1.vc1)
              PackSortedIndices(n_pair, idx1, sizeof(int), op1.vc1);
            if(op1.vp1)
              PackSortedIndices(n_pair, idx1, sizeof(int), op1.vp1);
            if(op1.ai1VLA)
              PackSortedIndices(n_pair, idx1, sizeof(AtomInfoType *), op1.ai1VLA);

            if(op2.vv1)
              PackSortedIndices(n_pair, idx2, 3 * sizeof(float), op2.vv1);
            if(op2.vc1)
              PackSortedIndices(n_pair, idx2, sizeof(int), op2.vc1);
            if(op2.vp1)
              PackSortedIndices(n_pair, idx2, sizeof(int), op2.vp1);
            if(op2.ai1VLA)
              PackSortedIndices(n_pair, idx2, sizeof(AtomInfoType *), op2.ai1VLA);
          }
        }
        FreeP(idx1);
        FreeP(idx2);
      } else if(op1.nvv1 != op2.nvv1) {
        sprintf(buffer, "Atom counts between selections don't match (%d vs %d)",
                op1.nvv1, op2.nvv1);
        ErrMessage(G, "ExecutiveRMS", buffer);
        n_pair = 0;
        ok = false;
      } else {
        n_pair = 0;
        for(a = 0; a < op1.nvv1; ++a) { // for atoms in selection
          if (op1.vc1[a] && op2.vc1[a]) { // check state counts
            if (n_pair < a) { // copy over if necessary
              copy3(op1.vv1 + 3 * a, op1.vv1 + 3 * n_pair);
              copy3(op2.vv1 + 3 * a, op2.vv1 + 3 * n_pair);
              if(op1.ai1VLA) op1.ai1VLA[n_pair] = op1.ai1VLA[a];
              if(op2.ai1VLA) op2.ai1VLA[n_pair] = op2.ai1VLA[a];
              if(op1.vp1) op1.vp1[n_pair] = op1.vp1[a];
              if(op2.vp1) op2.vp1[n_pair] = op2.vp1[a];
              op1.vc1[n_pair] = op1.vc1[a];
              op2.vc1[n_pair] = op2.vc1[a];
            }
            ++n_pair;
          }
        }
      }

      if(n_pair) {
        /* okay -- we're on track to do an alignment */

        if(ordered_selections && op1.vp1 && op2.vp1) {
          /* if we expected ordered selections and have priorities, 
             then we may need to sort vertices */

          int sort_flag1 = false, sort_flag2 = false;
          int well_defined1 = true, well_defined2 = true;

          for(a = 0; a < (n_pair - 1); a++) {
            /*          printf("op1 vertex %d priority %d\n",a,op1.vp1[a]);
               printf("op2 vertex %d priority %d\n",a,op2.vp1[a]); */

            if(op1.vp1[a] > op1.vp1[a + 1])
              sort_flag1 = true;
            else if(op1.vp1[a] == op1.vp1[a + 1])
              well_defined1 = false;
            if(op2.vp1[a] > op2.vp1[a + 1])
              sort_flag2 = true;
            else if(op2.vp1[a] == op2.vp1[a + 1])
              well_defined2 = false;
          }

          if(sort_flag1 || sort_flag2) {
            if(!(well_defined1 || well_defined2)) {
              PRINTFB(G, FB_Executive, FB_Warnings)
                "Executive-Warning: Ordering requested but not well defined.\n" ENDFB(G);
            } else {
              FitVertexRec *vert = pymol::malloc<FitVertexRec>(n_pair);

              if(sort_flag1) {
                float *src, *dst;
                src = op1.vv1;
                for(a = 0; a < n_pair; a++) {
                  vert[a].priority = op1.vp1[a];
                  dst = vert[a].vertex;
                  copy3f(src, dst);
                  src += 3;
                }
                UtilSortInPlace(G, vert, n_pair, sizeof(FitVertexRec),
                                (UtilOrderFn *) fVertexOrdered);
                dst = op1.vv1;
                for(a = 0; a < n_pair; a++) {
                  src = vert[a].vertex;
                  copy3f(src, dst);
                  dst += 3;
                }
              }

              if(sort_flag2) {
                float *src, *dst;
                src = op2.vv1;
                for(a = 0; a < n_pair; a++) {
                  vert[a].priority = op2.vp1[a];
                  dst = vert[a].vertex;
                  copy3f(src, dst);
                  src += 3;
                }
                UtilSortInPlace(G, vert, n_pair, sizeof(FitVertexRec),
                                (UtilOrderFn *) fVertexOrdered);
                dst = op2.vv1;
                for(a = 0; a < n_pair; a++) {
                  src = vert[a].vertex;
                  copy3f(src, dst);
                  dst += 3;
                }
              }

              FreeP(vert);
            }
          }
        }

        if(rms_info) {
          rms_info->initial_n_atom = n_pair;
          rms_info->n_cycles_run = 0;
          rms_info->final_n_atom = n_pair;      /* in case there is no refinement */
        }

        if(mode != 0) {
          rms = MatrixFitRMSTTTf(G, n_pair, op1.vv1, op2.vv1, NULL, op2.ttt);
          if(rms_info) {
            rms_info->initial_rms = rms;
            rms_info->final_rms = rms;
          }
          repeat = true;
          b = 0;
          while(repeat) {
            repeat = false;
            b++;
            if(b > max_cyc)
              break;
            if((refine > R_SMALL4) && (rms > R_SMALL4)) {
              int n_next = n_pair;
              AtomInfoType **ai1, **ai2;

              flag = pymol::malloc<int>(n_pair);

              if(flag) {
                for(a = 0; a < n_pair; a++) {
                  MatrixTransformTTTfN3f(1, v1, op2.ttt, op1.vv1 + (a * 3));
                  v2 = op2.vv1 + (a * 3);
                  if((diff3f(v1, v2) / rms) > refine) {
                    flag[a] = false;
                    repeat = true;
                  } else
                    flag[a] = true;
                }
                f1 = op1.vv1;
                f2 = op2.vv1;
                ai1 = op1.ai1VLA;
                ai2 = op2.ai1VLA;
                for(a = 0; a < n_pair; a++) {
                  if(!flag[a]) {
                    n_next--;
                  } else {
                    copy3f(op1.vv1 + (3 * a), f1);
                    copy3f(op2.vv1 + (3 * a), f2);
                    f1 += 3;
                    f2 += 3;
                    if(ai1 && ai2) {    /* make sure we keep track of which atoms are aligned */
                      *(ai1++) = op1.ai1VLA[a];
                      *(ai2++) = op2.ai1VLA[a];
                    }
                  }
                }
                if(!quiet && (n_next != n_pair)) {
                  PRINTFB(G, FB_Executive, FB_Actions)
                    " %s: %d atoms rejected during cycle %d (RMSD=%0.2f).\n", __func__,
                    n_pair - n_next, b, rms ENDFB(G);
                }
                n_pair = n_next;
                FreeP(flag);
                if(n_pair) {
                  rms = MatrixFitRMSTTTf(G, n_pair, op1.vv1, op2.vv1, NULL, op2.ttt);
                  if(rms_info) {
                    rms_info->n_cycles_run = b;
                    rms_info->final_n_atom = n_pair;
                    rms_info->final_rms = rms;
                  }
                } else
                  break;
              }
            }
          }
        } else {                /* mode == 0 -- simple RMS, with no coordinate movement */
          rms = MatrixGetRMS(G, n_pair, op1.vv1, op2.vv1, NULL);
          if(rms_info) {
            rms_info->initial_rms = rms;
            rms_info->final_rms = rms;
          }
        }
      }
      if(!n_pair) {
        PRINTFB(G, FB_Executive, FB_Results)
          " Executive: Error -- no atoms left after refinement!\n" ENDFB(G);
        ok = false;
      }

      if(ok) {
        if(!quiet) {
          PRINTFB(G, FB_Executive, FB_Results)
            " Executive: RMSD = %8.3f (%d to %d atoms)\n", rms, n_pair, n_pair ENDFB(G);
        }
        if(oname && oname[0]) {
            int align_state = state2;
            ObjectMolecule *trg_obj = SelectorGetSingleObjectMolecule(G, sele2);

            if(align_state < 0) {
              align_state = SceneGetState(G);
            }

            /* we're going to create/update an alignment object */

            {
              /* Get unique ids and construct the alignment vla */
              pymol::vla<int> align_vla(n_pair * 3);

              {
                int *id_p = align_vla.data();
                int i;
                for(i = 0; i < n_pair; i++) {
                  id_p[0] = AtomInfoCheckUniqueID(G, op2.ai1VLA[i]);    /* target */
                  id_p[1] = AtomInfoCheckUniqueID(G, op1.ai1VLA[i]);
                  id_p[2] = 0;
                  id_p += 3;
                }
                VLASize(align_vla, int, n_pair * 3);
              }
              {
                ObjectAlignment *obj = NULL;

                /* does object already exist? */
                {
                  pymol::CObject *execObj = ExecutiveFindObjectByName(G, oname);
                  if(execObj && (execObj->type != cObjectAlignment))
                    ExecutiveDelete(G, oname);
                  else
                    obj = (ObjectAlignment *) execObj;
                }
                obj =
                  ObjectAlignmentDefine(G, obj, align_vla, align_state, true, trg_obj,
                                        mobile_obj);
                obj->Color = ColorGetIndex(G, "yellow");
                ObjectSetName(obj, oname);
                ExecutiveManageObject(G, obj, 0, quiet);
                align_to_update = obj;
                SceneInvalidate(G);
              }
            }
        }
        if(ok && mode == 2) {
          if(matrix_mode>0) {

            ObjectMolecule *src_obj, *trg_obj;
            src_obj = SelectorGetFirstObjectMolecule(G, sele1); /* get at least one object */
            trg_obj = SelectorGetSingleObjectMolecule(G, sele2);

            /* first we need to make sure that the object being moved
               matches the target with respect to both the TTT and the
               object's state matrix (if any) */

            if(src_obj && trg_obj) {
              ExecutiveMatrixCopy(G, trg_obj->Name, src_obj->Name, 1, 1,        /* TTT mode */
                                  state2, state1, false, 0, quiet);

              ExecutiveMatrixCopy(G, trg_obj->Name, src_obj->Name, 2, 2,        /* Object state mode */
                                  state2, state1, false, 0, quiet);

              switch (matrix_mode) {
              case 1:          /* TTTs */
                ExecutiveCombineObjectTTT(G, src_obj->Name, op2.ttt, true, -1);
                break;
              case 2:
                {
                  double homo[16], *src_homo;
                  convertTTTfR44d(op2.ttt, homo);
                  if(ExecutiveGetObjectMatrix
                     (G, src_obj->Name, state1, &src_homo, false)) {
                    left_multiply44d44d(src_homo, homo);
                    ExecutiveSetObjectMatrix(G, src_obj->Name, state1, homo);
                  }
                }
                break;
              }
              /* next we need to update the object's TTT matrix to reflect
                 the transformation */
            }
          } else {              /* matrix_mode is zero -- legacy behavior */
            /* this will transform the actual coordinates */
            op2.code = OMOP_TTTF;
            ExecutiveObjMolSeleOp(G, sele1, &op2);
          }
        }
      }
    } else {
      ErrMessage(G, __func__, "No atoms selected.");
      ok = false;
    }
  }

  if(align_to_update) {
    align_to_update->update();
  }

  VLAFreeP(op1.vv1);
  VLAFreeP(op2.vv1);
  VLAFreeP(op1.vc1);
  VLAFreeP(op2.vc1);
  VLAFreeP(op1.vp1);
  VLAFreeP(op2.vp1);
  VLAFreeP(op1.ai1VLA);
  VLAFreeP(op2.ai1VLA);
  return (ok);
}


/*========================================================================*/
/**
 * Implementation of `cmd.identify()`
 *
 * @param s1 atom selection expression
 * @param mode 0 for index only, 1 for (obj, index)
 * @param[out] indexVLA
 * @param[out] objVLA
 * @return number of atoms or -1 on selection error
 */
int ExecutiveIdentifyObjects(PyMOLGlobals * G, const char *s1, int mode, int **indexVLA,
                             ObjectMolecule *** objVLA)
{
  SelectorTmp tmpsele1(G, s1);
  int sele1 = tmpsele1.getIndex();
  ObjectMoleculeOpRec op2;
  if(sele1 >= 0) {
    ObjectMoleculeOpRecInit(&op2);
    op2.code = OMOP_IdentifyObjects;
    if (mode != 0) {
      op2.obj1VLA = VLAlloc(ObjectMolecule*, 1000);
    }
    op2.i1VLA = VLAlloc(int, 1000);
    op2.i1 = 0;
    ExecutiveObjMolSeleOp(G, sele1, &op2);
    VLASize(op2.i1VLA, int, op2.i1);
    if (mode != 0) {
      VLASize(op2.obj1VLA, ObjectMolecule*, op2.i1);
    }
    (*indexVLA) = op2.i1VLA;
    (*objVLA) = op2.obj1VLA;
  } else {
    return -1;
  }
  return (op2.i1);
}


/*========================================================================*/
int ExecutiveIndex(PyMOLGlobals * G, const char *s1, int mode, int **indexVLA,
                   ObjectMolecule *** objVLA)
{
  ObjectMoleculeOpRec op2;

  SelectorTmp tmpsele1(G, s1);
  int sele1 = tmpsele1.getIndex();

  if(sele1 >= 0) {
    ObjectMoleculeOpRecInit(&op2);
    op2.code = OMOP_Index;
    op2.obj1VLA = VLAlloc(ObjectMolecule *, 1000);
    op2.i1VLA = VLAlloc(int, 1000);
    op2.i1 = 0;
    ExecutiveObjMolSeleOp(G, sele1, &op2);
    VLASize(op2.i1VLA, int, op2.i1);
    VLASize(op2.obj1VLA, ObjectMolecule *, op2.i1);
    (*indexVLA) = op2.i1VLA;
    (*objVLA) = op2.obj1VLA;
  } else {
    return -1; // invalid selection
  }
  return (op2.i1);
}


/*========================================================================*/
/**
 * Fit states or calculate ensemble RMSD
 *
 * @param s1 atom selection expression
 * @param target reference state
 * @param mode 2=intra_fit, 1=intra_rms, 0=intra_rms_cur
 * @param mix intra_fit only, average the prior target coordinates
 * @param pbc Consider periodic boundary conditions
 */
pymol::Result<pymol::vla<float>> ExecutiveRMSStates(
    PyMOLGlobals * G, const char *s1, int target, int mode, int quiet, int mix,
    bool pbc)
{
  SelectorTmp tmpsele1(G, s1);
  int sele1 = tmpsele1.getIndex();

  ObjectMoleculeOpRec op1;
  ObjectMoleculeOpRec op2;
  float *result = NULL;
  int ok = true;

  ObjectMoleculeOpRecInit(&op1);
  ObjectMoleculeOpRecInit(&op2);
  op1.vv1 = NULL;
  op2.vv1 = NULL;

  ObjectMolecule* obj = SelectorGetSingleObjectMolecule(G, sele1);

  if (!obj) {
    if(mode != 2) {
      PRINTFB(G, FB_Executive, FB_Warnings)
        "Executive-Warning: Mobile selection spans more than one object.\n" ENDFB(G);
    } else {
      return pymol::make_error("Mobile selection spans more than one object.");
    }
  }

  if (target == cStateCurrent) {
    target = obj ? obj->getCurrentState() : SceneGetState(G);
  }

  if (target < 0) {
    target = 0;
  }

  if (mode != 2) {
    pbc = false;
  }

  if(ok && sele1 >= 0) {
    op1.code = OMOP_SVRT;
    op1.nvv1 = 0;
    op1.i1 = target;
    op1.vv1 = (float *) VLAMalloc(1000, sizeof(float), 5, 0);
    op1.i1VLA = VLAlloc(int, 1000);
    ExecutiveObjMolSeleOp(G, sele1, &op1);

    if (pbc) {
      ObjectMoleculePBCUnwrap(*obj);
    }

    op2.vv2 = op1.vv1;
    op2.nvv2 = op1.nvv1;
    op2.i1VLA = op1.i1VLA;
    op2.i2 = target;
    op2.i1 = mode;
    op2.i3 = mix;
    op2.f1VLA = VLAlloc(float, 10);
    VLASize(op2.f1VLA, float, 0);       /* failsafe */
    op2.vv1 = (float *) VLAMalloc(1000, sizeof(float), 5, 0);
    op2.code = OMOP_SFIT;
    op2.nvv1 = 0;
    ExecutiveObjMolSeleOp(G, sele1, &op2);
    result = op2.f1VLA;
    VLAFreeP(op1.i1VLA);
    VLAFreeP(op2.vv1);

    if (pbc) {
      float center[3];
      pymol::meanNx3(op1.vv1, op1.nvv1, center);
      ObjectMoleculePBCWrap(*obj, center);
    }

    VLAFreeP(op1.vv1);

    if (mode == 2) {
      ExecutiveUpdateCoordDepends(G, obj);
    }
  }
  return pymol::vla_take_ownership(result);
}


/*========================================================================*/
float ExecutiveRMSPairs(PyMOLGlobals* G, const std::vector<SelectorTmp>& sele,
    int mode, bool quiet)
{
  int sele1, sele2;
  int a, c;
  float rms = -1.0, inv, *f;
  OrthoLineType buffer;

  ObjectMoleculeOpRec op1;
  ObjectMoleculeOpRec op2;
  OrthoLineType combi, s1;

  ObjectMoleculeOpRecInit(&op1);
  ObjectMoleculeOpRecInit(&op2);
  op1.nvv1 = 0;
  op1.vc1 = (int *) VLAMalloc(1000, sizeof(int), 5, 1);
  op1.vv1 = (float *) VLAMalloc(1000, sizeof(float), 5, 1);     /* auto-zero */
  op1.code = OMOP_AVRT;

  op2.nvv1 = 0;
  op2.vc1 = (int *) VLAMalloc(1000, sizeof(int), 5, 1);
  op2.vv1 = (float *) VLAMalloc(1000, sizeof(float), 5, 1);     /* auto-zero */
  op2.code = OMOP_AVRT;

  strcpy(combi, "(");
  c = 0;
  auto pairs = sele.size() / 2;
  for(a = 0; a < pairs; a++) {
    sele1 = sele[c].getIndex();
    if(sele1 >= 0)
      ExecutiveObjMolSeleOp(G, sele1, &op1);
    strcat(combi, sele[c].getName());
    if(a < (pairs - 1))
      strcat(combi, " or ");
    c++;
    sele2 = sele[c].getIndex();
    if(sele2 >= 0)
      ExecutiveObjMolSeleOp(G, sele2, &op2);
    c++;
  }
  strcat(combi, ")");
  for(a = 0; a < op1.nvv1; a++) {
    inv = (float) op1.vc1[a];
    if(inv) {
      f = op1.vv1 + (a * 3);
      inv = 1.0F / inv;
      *(f++) *= inv;
      *(f++) *= inv;
      *(f++) *= inv;
    }
  }
  for(a = 0; a < op2.nvv1; a++) {
    inv = (float) op2.vc1[a];
    if(inv) {
      f = op2.vv1 + (a * 3);
      inv = 1.0F / inv;
      *(f++) *= inv;
      *(f++) *= inv;
      *(f++) *= inv;
    }
  }
  if(op1.vv1 && op2.vv1) {
    if(op1.nvv1 != op2.nvv1) {
      sprintf(buffer, "Atom counts between selection sets don't match (%d != %d).",
              op1.nvv1, op2.nvv1);
      ErrMessage(G, __func__, buffer);
    } else if(op1.nvv1) {
      if(mode != 0)
        rms = MatrixFitRMSTTTf(G, op1.nvv1, op1.vv1, op2.vv1, NULL, op2.ttt);
      else
        rms = MatrixGetRMS(G, op1.nvv1, op1.vv1, op2.vv1, NULL);

      if (!quiet)
      PRINTFB(G, FB_Executive, FB_Results)
        " %s: RMSD = %8.3f (%d to %d atoms)\n", __func__, rms, op1.nvv1, op2.nvv1 ENDFB(G);

      op2.code = OMOP_TTTF;
      SelectorGetTmp(G, combi, s1);
      sele1 = SelectorIndexByName(G, s1);
      ExecutiveObjMolSeleOp(G, sele1, &op2);
      SelectorFreeTmp(G, s1);
    } else {
      ErrMessage(G, __func__, "No atoms selected.");
    }
  }
  VLAFreeP(op1.vv1);
  VLAFreeP(op2.vv1);
  VLAFreeP(op1.vc1);
  VLAFreeP(op2.vc1);
  return (rms);
}


/*========================================================================*/
void ExecutiveUpdateObjectSelection(PyMOLGlobals * G, pymol::CObject * obj)
{
  if(obj->type == cObjectMolecule) {
    SelectorUpdateObjectSele(G, (ObjectMolecule *) obj);
  }
}


/*========================================================================*/
/**
 * Reset camera view or object TTT matrix. Stores key frames for modified
 * objects if `movie_auto_store=on`.
 *
 * @param name Empty, "all", "same", or object name pattern
 *
 * - empty name: Reset camera view
 * - "all": Reset TTT matrices and store key frames for all objects
 * - "same": Reset TTT matrices and store key frames for objects which currently have any key frames
 * - pattern: Reset TTT matrices and store key frames for objects which match the pattern
 */
pymol::Result<> ExecutiveReset(PyMOLGlobals* G, pymol::zstring_view name)
{
  if (name.empty()) {
    SceneResetMatrix(G);
    ExecutiveWindowZoom(G, cKeywordAll, 0.0, -1, 0, 0, true);   /* reset does all states */
    return {};
  }

  bool do_reset_all = name == cKeywordAll;
  auto store = SettingGet<bool>(G, cSetting_movie_auto_store);

  /**
   * @param any_spec_level If false, then filter for objects with spec level >= 0
   */
  auto reset_rec = [&](SpecRec& rec, bool any_spec_level = true) {
    pymol::CObject* obj = rec.obj;
    if (rec.type == cExecObject &&
        (any_spec_level || ObjectGetSpecLevel(obj, 0) >= 0)) {
      ObjectResetTTT(obj, store);
      obj->invalidate(cRepNone, cRepInvExtents, -1);
    }
  };

  if (do_reset_all || name == cKeywordSame) {
    for (auto& rec : pymol::make_list_adapter(G->Executive->Spec)) {
      reset_rec(rec, do_reset_all);
    }
  } else {
    for (auto& rec : ExecutiveGetSpecRecsFromPattern(G, name)) {
      reset_rec(rec);
    }
  }

  if (store && SettingGet<bool>(G, cSetting_movie_auto_interpolate)) {
    ExecutiveMotionReinterpolate(G);
  }

  SceneInvalidate(G);

  return {};
}


/*========================================================================*/
void ExecutiveDrawNow(PyMOLGlobals * G)
{
  CExecutive *I = G->Executive;

  if(PyMOL_GetIdleAndReady(G->PyMOL) && !SettingGetGlobal_b(G, cSetting_suspend_deferred))
    OrthoExecDeferred(G);
  if(!SettingGetGlobal_b(G, cSetting_suspend_updates)){
    int stereo_mode = SettingGetGlobal_i(G, cSetting_stereo_mode);
    int stereo = SettingGetGlobal_i(G, cSetting_stereo);
    if(G->HaveGUI && G->ValidContext) {
      glMatrixMode(GL_MODELVIEW);       /* why is this necessary?  is it? */
    }

    ExecutiveUpdateSceneMembers(G);
    SceneUpdate(G, false);
    if(WizardUpdate(G))
      SceneUpdate(G, false);
    if (stereo){
      switch (stereo_mode) {
      case cStereo_geowall:
	{
	  int width = G->Option->winX;
	  int height = G->Option->winY;
	  glViewport(0, 0, width / 2, height);
	  OrthoDoDraw(G, OrthoRenderMode::GeoWallLeft);
	  OrthoDoDraw(G, OrthoRenderMode::GeoWallRight);
	  glViewport(0, 0, width, height);
	}
	break;
#ifdef _PYMOL_OPENVR
      case cStereo_openvr:
        {
          Block* scene_block = SceneGetBlock(G);
          int scene_width = scene_block->rect.right - scene_block->rect.left;
          int scene_height = scene_block->rect.top - scene_block->rect.bottom;
          OpenVRFrameStart(G);
          float matrix[16];
          SceneGetModel2WorldMatrix(G, matrix);
          OpenVRHandleInput(G, scene_block->rect.left, scene_block->rect.bottom, scene_width, scene_height, matrix);
          OrthoDoDraw(G, OrthoRenderMode::VR);
          if (SettingGetGlobal_b(G, cSetting_openvr_cut_laser) && OpenVRIsScenePickerActive(G)) {
            int x = scene_block->rect.left + scene_width / 2;
            int y = scene_block->rect.bottom + scene_height / 2;
            float atomWorldPos[3];
            ScenePickAtomInWorld(G, x, y, atomWorldPos);
            OpenVRUpdateScenePickerLength(G, atomWorldPos);
          }
          OpenVRFrameFinish(G);
          PyMOL_NeedRedisplay(G->PyMOL);
        }
        break;
#endif
      default:
	OrthoDoDraw(G, OrthoRenderMode::Main);
	break;
      }
    } else {
      OrthoDoDraw(G, OrthoRenderMode::Main);
    }

    if(G->HaveGUI && G->ValidContext) {
      if(I->CaptureFlag) {
        I->CaptureFlag = false;
        SceneCaptureWindow(G);
      }
    }
    PyMOL_NeedSwap(G->PyMOL);
  }

  //  PRINTFD(G, FB_Executive)
  //    " ExecutiveDrawNow: leaving.\n" ENDFD;
}


/*========================================================================*/
int ExecutiveCountStates(PyMOLGlobals * G, const char *s1)
{
  CExecutive *I = G->Executive;
  int sele1;
  int result = 0;
  int n_state;
  SpecRec* list_rec = nullptr;
  if((!s1) || (!s1[0]))
    s1 = cKeywordAll;
  for (const auto& rec : ExecutiveGetSpecRecsFromPattern(G, s1)) {
    switch (rec.type) {
    case cExecAll:
      while(ListIterate(I->Spec, list_rec, next)) {
        if(list_rec->type == cExecObject) {
            n_state = list_rec->obj->getNFrame();
            if(result < n_state)
              result = n_state;
        }
      }
      break;
    case cExecSelection:
      sele1 = SelectorIndexByName(G, rec.name);
      if(sele1 >= 0) {
        SelectorUpdateTable(G, cSelectorUpdateTableAllStates, -1);
        n_state = SelectorGetSeleNCSet(G, sele1);
        if(result < n_state)
          result = n_state;
      }
      break;
    case cExecObject:
        n_state = rec.obj->getNFrame();
        if(result < n_state)
          result = n_state;
      break;
    }
  }
  return (result);
}


/*========================================================================*/
int ExecutiveRay(PyMOLGlobals * G, int width, int height, int mode,
                 float angle, float shift, int quiet, int defer, int antialias)
{
  if((mode == 0) && G->HaveGUI && SettingGetGlobal_b(G, cSetting_auto_copy_images)) {
    /* force deferred behavior if copying image to clipboard */
    defer = 1;
  }

  ExecutiveUpdateSceneMembers(G);

  if(defer && (mode == 0)) {
    SceneDeferRay(G, width, height, mode, angle, shift, quiet, true, antialias);
  } else {
    SceneRay(G, width, height, mode, NULL, NULL, angle, shift, quiet, NULL, true,
               antialias);
  }
  return 1;
}


/*========================================================================*/
int *ExecutiveGetG3d(PyMOLGlobals * G)
{
  int *result = NULL;
  SceneRay(G, 0, 0, 3, NULL, NULL, 0.0F, 0.0F, true, (G3dPrimitive **) (void *) &result,
           false, -1);
  return result;
}

int ExecutiveSetBondSettingFromString(PyMOLGlobals * G,
                                      int index, const char *value,
                                      const char *s1, const char *s2, int state,
                                      int quiet, int updates)
{
  CExecutive *I = G->Executive;
  SpecRec *rec = NULL;
  ObjectMolecule *obj = NULL;
  int sele1, sele2;
  SettingName name;
  int ok = true;
  int side_effects = false;
  int value_storage[3], *value_ptr;
  float float_storage[3];
  int value_type = 0;
  PRINTFD(G, FB_Executive)
    " %s: entered. '%s' '%s'\n", __func__, s1, s2 ENDFD;
  sele1 = SelectorIndexByName(G, s1);
  sele2 = SelectorIndexByName(G, s2);
  value_ptr = &value_storage[0];
  if((sele1 >= 0) && (sele2 >= 0)) {
    int have_value = false;
    int type = SettingGetType(G, index);
    switch (type) {
    case cSetting_boolean:
      {
        if((!*value) || (*value == '0') || (*value == 'F') || WordMatchExact(G, value, "on", true)
           || WordMatchExact(G, value, "false", true))
          *(value_ptr) = 0;
        else
          *(value_ptr) = 1;
        value_type = cSetting_boolean;
        have_value = true;
      }
      break;
    case cSetting_int:
      {
        if(sscanf(value, "%d", value_ptr) == 1) {
          value_type = cSetting_int;
          have_value = true;
        } else {
          ok = false;
        }
      }
      break;
    case cSetting_float:
      {
        if(sscanf(value, "%f", &float_storage[0]) == 1) {
          value_ptr = (int*) (void*) &float_storage[0];
          value_type = cSetting_float;
          have_value = true;
        } else {
          ok = false;
        }
      }
      break;
    case cSetting_float3:
      if(sscanf(value, "%f%f%f", &float_storage[0],
                &float_storage[1], &float_storage[2]) == 3) {
        value_ptr = (int*) (void*) &float_storage[0];
        value_type = cSetting_float3;
        have_value = true;
      } else {
        ok = false;
      }
      break;
    case cSetting_color:
      {
        int color_index = ColorGetIndex(G, value);
        if((color_index < 0) && (color_index > cColorExtCutoff))
          color_index = 0;
        *(value_ptr) = color_index;
        value_type = cSetting_color;
        have_value = true;
      }
      break;
      /* cSetting_string? */
    default:
      ok = false;
      break;
    }
    
    if(ok && have_value) {
      rec = NULL;
      while((ListIterate(I->Spec, rec, next))) {
        if((rec->type == cExecObject) && (rec->obj->type == cObjectMolecule)) {
          obj = (ObjectMolecule *) rec->obj;
          {
            int a, nBond = obj->NBond;
            int nSet = 0;
            BondType *bi = obj->Bond.data();
            const AtomInfoType *ai1, *ai2, *ai = obj->AtomInfo.data();
            for(a = 0; a < nBond; a++) {
              ai1 = ai + bi->index[0];
              ai2 = ai + bi->index[1];
              if((SelectorIsMember(G, ai1->selEntry, sele1) &&
                  SelectorIsMember(G, ai2->selEntry, sele2)) ||
                 (SelectorIsMember(G, ai2->selEntry, sele1) &&
                  SelectorIsMember(G, ai1->selEntry, sele2))) {
                
                int uid = AtomInfoCheckUniqueBondID(G, bi);
		int isset;
                bi->has_setting = true;
                isset = SettingUniqueSetTypedValue(G, uid, index, value_type, value_ptr);
                if(updates && isset)
                  side_effects = true;
                nSet++;
              }
              bi++;
            }
            if(nSet && !quiet) {
              SettingGetName(G, index, name);
              PRINTF
                " Setting: %s set for %d bonds in object \"%s\".\n",
                name, nSet, obj->Name ENDF(G);
            }
          }
        }
      }
    }
  }
  if(side_effects) {
    SettingGenerateSideEffects(G, index, s1, state, quiet); /* not strickly correct */
    /*    SettingGenerateSideEffects(G,index,s2,state); */
  }
  return (ok);
}
/*========================================================================*/
/**
 * @pre GIL
 */
PyObject *ExecutiveGetBondSetting(PyMOLGlobals * G, int index, 
				  char *s1, const char *s2, int state, int quiet, int updates)
{
#ifdef _PYMOL_NOPY
  return 0;
#else
  assert(PyGILState_Check());

  CExecutive *I = G->Executive;
  SpecRec *rec = NULL;
  ObjectMolecule *obj = NULL;
  int sele1, sele2;
  SettingName name;
  PyObject *result = PyList_New(0);

  sele1 = SelectorIndexByName(G, s1);
  sele2 = SelectorIndexByName(G, s2);

  if((sele1 >= 0) && (sele2 >= 0)) {
    while((ListIterate(I->Spec, rec, next))) {
      if((rec->type == cExecObject) && (rec->obj->type == cObjectMolecule)) {
	obj = (ObjectMolecule *) rec->obj;
	{
	  int a, nBond = obj->NBond ;
	  int nSet = 0;
	  const BondType *bi = obj->Bond.data();
	  const AtomInfoType *ai1, *ai2, *ai = obj->AtomInfo.data();

	  PyObject *pyObjList = NULL;
	  PyObject *pyBondList = NULL;

	  for(a = 0; a < nBond; a++) {
	    ai1 = ai + bi->index[0];
	    ai2 = ai + bi->index[1];
	    if((SelectorIsMember(G, ai1->selEntry, sele1) &&
		SelectorIsMember(G, ai2->selEntry, sele2)) ||
	       (SelectorIsMember(G, ai2->selEntry, sele1) &&
		SelectorIsMember(G, ai1->selEntry, sele2))) {
	      PyObject *pyBondInfo = PyList_New(3);
	      PyObject *bond_setting_value = NULL;
	      if (!pyObjList){
		pyObjList = PyList_New(2);
		pyBondList = PyList_New(0);
		PyList_SetItem(pyObjList, 0, PyString_FromString(obj->Name));
		PyList_SetItem(pyObjList, 1, pyBondList);
		PyList_Append(result, pyObjList);
		Py_DECREF(pyObjList);
	      }
	      PyList_SetItem(pyBondInfo, 0, PyInt_FromLong((long)bi->index[0]+1));
	      PyList_SetItem(pyBondInfo, 1, PyInt_FromLong((long)bi->index[1]+1));
	      if (bi->has_setting){
		bond_setting_value = SettingUniqueGetPyObject(G, bi->unique_id, index);
	      }
	      PyList_SetItem(pyBondInfo, 2, PConvAutoNone(bond_setting_value));
	      PyList_Append(pyBondList, pyBondInfo);
	      Py_DECREF(pyBondInfo);
	      nSet++;
	    }
	    bi++;
	  }
	  if(nSet && !quiet) {
	    SettingGetName(G, index, name);
	    PRINTF
	      " Getting: %s for %d bonds in object \"%s\".\n",
	      name, nSet, obj->Name ENDF(G);
	  }
	}
      }
    }
  }
  return result;
#endif
}
/*========================================================================*/
int ExecutiveSetBondSetting(PyMOLGlobals * G, int index, PyObject * tuple,
                            const char *s1, const char *s2, int state, int quiet, int updates)
{
#ifdef _PYMOL_NOPY
  return 0;
#else

  CExecutive *I = G->Executive;
  SpecRec *rec = NULL;
  ObjectMolecule *obj = NULL;
  int sele1, sele2;
  SettingName name = "";
  int unblock;
  int ok = true;
  int side_effects = false;
  union {
    int value_ptr[1];
    float float_storage[1];
  };
  int value_type = 0;

  PRINTFD(G, FB_Executive)
    " %s: entered. '%s' '%s'\n", __func__, s1, s2 ENDFD;
  unblock = PAutoBlock(G);
  sele1 = SelectorIndexByName(G, s1);
  sele2 = SelectorIndexByName(G, s2);
  if((sele1 >= 0) && (sele2 >= 0)) {
    int have_value = false;
    int type = PyInt_AsLong(PyTuple_GetItem(tuple, 0));
    PyObject *value = PyTuple_GetItem(tuple, 1);
    if(value) {
      switch (type) {
      case cSetting_boolean:
        *(value_ptr) = PyInt_AsLong(value);
        value_type = cSetting_boolean;
        have_value = true;
        break;
      case cSetting_int:
        *(value_ptr) = PyInt_AsLong(value);
        value_type = cSetting_int;
        have_value = true;
        break;
      case cSetting_float:
        float_storage[0] = PyFloat_AsDouble(value);
        value_type = cSetting_float;
        have_value = true;
        break;
      case cSetting_color:
        {
          int color_index =
            ColorGetIndex(G, PyString_AsString(value));
          if((color_index < 0) && (color_index > cColorExtCutoff))
            color_index = 0;
          *(value_ptr) = color_index;
          value_type = cSetting_color;
          have_value = true;
        }
        break;
      }
      if(have_value) {
        rec = NULL;
        while((ListIterate(I->Spec, rec, next))) {
          if((rec->type == cExecObject) && (rec->obj->type == cObjectMolecule)) {
            obj = (ObjectMolecule *) rec->obj;
            {
              int a, nBond = obj->NBond;
              int nSet = 0;
              BondType *bi = obj->Bond.data();
              const AtomInfoType *ai1, *ai2, *ai = obj->AtomInfo.data();
              for(a = 0; a < nBond; a++) {
                ai1 = ai + bi->index[0];
                ai2 = ai + bi->index[1];
                if((SelectorIsMember(G, ai1->selEntry, sele1) &&
                    SelectorIsMember(G, ai2->selEntry, sele2)) ||
                   (SelectorIsMember(G, ai2->selEntry, sele1) &&
                    SelectorIsMember(G, ai1->selEntry, sele2))) {

                  int uid = AtomInfoCheckUniqueBondID(G, bi);
                  bi->has_setting = true;
                  SettingUniqueSetTypedValue(G, uid, index, value_type, value_ptr);
                  if(updates)
                    side_effects = true;
                  nSet++;
                }
                bi++;
              }
              if(nSet && !quiet) {
                SettingGetName(G, index, name);
                PRINTF
                  " Setting: %s set for %d bonds in object \"%s\".\n",
                  name, nSet, obj->Name ENDF(G);
              }
            }
          }
        }
      }
    }
  }
  if(side_effects) {
    SettingGenerateSideEffects(G, index, s1, state, quiet);  /* not strictly correct */
    /*    SettingGenerateSideEffects(G,index,s2,state); */
  }

  if(!SettingLevelCheck(G, index, cSettingLevel_bond)) {
    if (!name[0])
      SettingGetName(G, index, name);

    PRINTFB(G, FB_Setting, FB_Warnings)
      " Setting-Warning: '%s' is not a bond-level setting\n", name
      ENDFB(G);
  }

  PAutoUnblock(G, unblock);
  return (ok);
#endif
}


/*========================================================================*/
/**
 * @return Always true
 */
int ExecutiveUnsetBondSetting(PyMOLGlobals * G, int index, const char *s1, const char *s2,
                              int state, int quiet, int updates)
{
  CExecutive *I = G->Executive;
  SpecRec *rec = NULL;
  ObjectMolecule *obj = NULL;
  SettingName name;
  /*  int unblock; */
  int ok = true;
  int side_effects = false;
  int sele1, sele2;
  PRINTFD(G, FB_Executive)
    " %s: entered. sele '%s' '%s'\n", __func__, s1, s2 ENDFD;
  /* unblock = PAutoBlock(G); */
  sele1 = SelectorIndexByName(G, s1);
  sele2 = SelectorIndexByName(G, s2);
  if((sele1 >= 0) && (sele2 >= 0)) {
    rec = NULL;
    while((ListIterate(I->Spec, rec, next))) {
      if((rec->type == cExecObject) && (rec->obj->type == cObjectMolecule)) {
        obj = (ObjectMolecule *) rec->obj;
        {
          int nSet = 0;
          BondType *bi = obj->Bond.data();
          BondType *bi_end = bi + obj->NBond;
          AtomInfoType *ai1, *ai2, *ai = obj->AtomInfo.data();
          for(; bi != bi_end; ++bi) {
            if(!bi->has_setting)
              continue;
            ai1 = ai + bi->index[0];
            ai2 = ai + bi->index[1];
            if((SelectorIsMember(G, ai1->selEntry, sele1) &&
                SelectorIsMember(G, ai2->selEntry, sele2)) ||
               (SelectorIsMember(G, ai2->selEntry, sele1) &&
                SelectorIsMember(G, ai1->selEntry, sele2))) {
              int uid = AtomInfoCheckUniqueBondID(G, bi);
              if(!SettingUniqueUnset(G, uid, index))
                continue;
              if(updates)
                side_effects = true;
              nSet++;
            }
          }
          if(nSet && !quiet) {
            SettingGetName(G, index, name);
            PRINTF
              " Setting: %s unset for %d bonds in object \"%s\".\n",
              name, nSet, rec->obj->Name ENDF(G);
          }
        }
      }
    }
  }
  if(side_effects) {
    SettingGenerateSideEffects(G, index, s1, state, quiet);
    /*    SettingGenerateSideEffects(G,index,s2,state); */
  }
  /* PAutoUnblock(G, unblock); */
  return (ok);
}


/*========================================================================*/
pymol::Result<> ExecutiveSetSetting(PyMOLGlobals * G, int index, PyObject * tuple,
                        pymol::zstring_view preSele, int state, int quiet, int updates)
{
#ifdef _PYMOL_NOPY
  return pymol::make_error("Python not available.");
#else

  CExecutive *I = G->Executive;
  SpecRec *rec = NULL;
  ObjectMolecule *obj = NULL;
  int sele1;
  ObjectMoleculeOpRec op;
  OrthoLineType value;
  pymol::copyable_ptr<CSetting>* handle = nullptr;
  SettingName name = "";
  int nObj = 0;
  const char* sele = preSele.c_str();
  int ok = true;

  pymol::Result<SelectorTmp2> s1;

  if (!preSele.empty())
  {
    s1 = SelectorTmp2::make(G, preSele.c_str());
    p_return_if_error(s1);
    sele = s1->getName();
  }

  PRINTFD(G, FB_Executive)
    " %s: entered. sele \"%s\" updates=%d index=%d\n", __func__, sele, updates, index ENDFD;

  if(!quiet) {
    SettingGetName(G, index, name);
  }

  pymol::pautoblock unblock(G);

  if((!sele) || (sele[0] == 0)) {       /* global setting */
    ok = SettingSetFromTuple(G, NULL, index, tuple);
    if(ok) {
      if(!quiet) {
        if(Feedback(G, FB_Setting, FB_Actions)) {
          SettingGetTextValue(G, NULL, NULL, index, value);
          PRINTF " Setting: %s set to %s.\n", name, value ENDF(G);
        }
      }
      if(updates) {
        SettingGenerateSideEffects(G, index, NULL, state, quiet);
      }
    }
  }
  else {
    unsigned char levelmask = 0;
    int side_effects = false;

    CTracker *I_Tracker = I->Tracker;
    int list_id = ExecutiveGetNamesListFromPattern(G, sele, true, true);
    int iter_id = TrackerNewIter(I_Tracker, 0, list_id);
    while(TrackerIterNextCandInList(I_Tracker, iter_id, (TrackerRef **) (void *) &rec)) {
      if(rec) {
        switch (rec->type) {
        case cExecAll:
          levelmask |= SettingLevelInfo[state < 0 ? cSettingLevel_object : cSettingLevel_ostate].mask;
          rec = NULL;
          while(ListIterate(I->Spec, rec, next)) {
            if(rec->type == cExecObject) {
              {
                handle = rec->obj->getSettingHandle(state);
                if(handle) {
                  SettingCheckHandle(G, *handle);
                  ok = SettingSetFromTuple(G, handle->get(), index, tuple);
                  if(updates)
                    side_effects = true;
                  nObj++;
                }
              }
            }
          }
          if(Feedback(G, FB_Setting, FB_Actions)) {
            if(nObj && handle) {
              SettingGetTextValue(G, handle->get(), NULL, index, value);
              if(!quiet) {
                if(state < 0) {
                  PRINTF
                    " Setting: %s set to %s in %d objects.\n", name, value, nObj ENDF(G);
                } else {
                  PRINTF
                    " Setting: %s set to %s in %d objects, state %d.\n",
                    name, value, nObj, state + 1 ENDF(G);
                }
              }
            }
          }
          break;
        case cExecSelection:
          if (SettingLevelCheckMask(G, index, SettingLevelInfo[cSettingLevel_bond].mask)) {
            // handle bond-level settings (PYMOL-2726)
            ok = ExecutiveSetBondSetting(G, index, tuple, sele, sele, state, quiet, false);
            if (updates)
              side_effects = true;
            sele1 = -1;
          } else {
            levelmask |= SettingLevelInfo[cSettingLevel_atom].mask;
            sele1 = SelectorIndexByName(G, rec->name);
          }

          if(sele1 >= 0) {
            int have_atomic_value = false;
            int type = PyInt_AsLong(PyTuple_GetItem(tuple, 0));
            PyObject *value = PyTuple_GetItem(tuple, 1);
            if(value) {
              ObjectMoleculeOpRecInit(&op);
              op.code = OMOP_SetAtomicSetting;
              op.i1 = index;
              op.ii1 = &op.i3;
              switch (type) {
              case cSetting_boolean:
                *(op.ii1) = PyInt_AsLong(value);
                op.i2 = cSetting_boolean;
                have_atomic_value = true;
                break;
              case cSetting_int:
                *(op.ii1) = PyInt_AsLong(value);
                op.i2 = cSetting_int;
                have_atomic_value = true;
                break;
              case cSetting_float:
                *(float *) op.ii1 = (float) PyFloat_AsDouble(value);
                op.i2 = cSetting_float;
                have_atomic_value = true;
                break;
              case cSetting_float3:
		{
		  PConvPyListOrTupleToFloatArrayInPlace(value, op.ttt, 3);
		  op.mat1 = op.ttt; // for passing (float**)
		  op.ii1 = (int*) &op.mat1;
		  op.i2 = cSetting_float3;
		  have_atomic_value = true;
		}
                break;
              case cSetting_color:
                {
                  int color_index =
                    ColorGetIndex(G, PyString_AsString(value));
                  if((color_index < 0) && (color_index > cColorExtCutoff)) {
                    switch (color_index) {
                    case cColorAtomic:
                      color_index = -1;
                      break;
                    case cColorFront:
                    case cColorBack:
                    case cColorDefault:
                      break;
                    default:
                      color_index = 0;
                      break;
                    }
                  }
                  *(op.ii1) = color_index;
                  op.i2 = cSetting_color;
                  have_atomic_value = true;
                }
                break;
              }
              if(have_atomic_value) {
                rec = NULL;
                while((ListIterate(I->Spec, rec, next))) {
                  if((rec->type == cExecObject) && (rec->obj->type == cObjectMolecule)) {
                    obj = (ObjectMolecule *) rec->obj;
                    op.i4 = 0;
                    ObjectMoleculeSeleOp(obj, sele1, &op);
                    if(op.i4) {
                      if(updates)
                        side_effects = true;
                      if(!quiet) {
                        PRINTF
                          " Setting: %s set for %d atoms in object \"%s\".\n",
                          name, op.i4, rec->obj->Name ENDF(G);
                      }
                    }
                  }
                }
              }
            }
          }
          break;
        case cExecObject:
          levelmask |= SettingLevelInfo[state < 0 ? cSettingLevel_object : cSettingLevel_ostate].mask;
          {
            handle = rec->obj->getSettingHandle(state);
            if(handle) {
              SettingCheckHandle(G, *handle);
              ok = SettingSetFromTuple(G, handle->get(), index, tuple);
              if(ok) {
                if(updates)
                  side_effects = true;
                if(!quiet) {
                  if(state < 0) {       /* object-specific */
                    if(Feedback(G, FB_Setting, FB_Actions)) {
                      SettingGetTextValue(G, handle->get(), NULL, index, value);
                      PRINTF
                        " Setting: %s set to %s in object \"%s\".\n",
                        name, value, rec->obj->Name ENDF(G);
                    }
                  } else {      /* state-specific */
                    if(Feedback(G, FB_Setting, FB_Actions)) {
                      SettingGetTextValue(G, handle->get(), NULL, index, value);
                      PRINTF
                        " Setting: %s set to %s in object \"%s\", state %d.\n",
                        name, value, rec->obj->Name, state + 1 ENDF(G);
                    }
                  }
                }
              }
            }
          }
          break;
        }
      }
    }
    TrackerDelList(I_Tracker, list_id);
    TrackerDelIter(I_Tracker, iter_id);

    if(side_effects)
      SettingGenerateSideEffects(G, index, sele, state, quiet);

    if(!SettingLevelCheckMask(G, index, levelmask)) {
      if(!name[0])
        SettingGetName(G, index, name);

      PRINTFB(G, FB_Setting, FB_Warnings)
        " Setting-Warning: '%s' is a %s-level setting\n", name, SettingLevelGetName(index)
        ENDFB(G);
    }
  }

  if (!ok) {
    return pymol::make_error("Error");
  }
  return {};
#endif
}

int ExecutiveGetSettingFromString(PyMOLGlobals * G, PyMOLreturn_value *result, 
                                  int index, const char *sele,
                                  int state, int quiet)
{
  pymol::CObject *obj = NULL;
  CSetting *set_ptr1 = NULL, *set_ptr2 = NULL;
  pymol::copyable_ptr<CSetting>* handle = nullptr;
  int ok = true;
  int type;
  type = SettingGetType(G, index);
  if(sele)
    if(sele[0]) {
      obj = ExecutiveFindObjectByName(G, sele);
      if(!obj)
        ok = false;
    }
  if(!ok) {
    PRINTFB(G, FB_Executive, FB_Errors)
      " %s-Error: sele \"%s\" not found.\n", __func__, sele ENDFB(G);
    ok = false;
  } else if(obj) {
    handle = obj->getSettingHandle(-1);
    if(handle)
      set_ptr1 = handle->get();
    if(state >= 0) {
      handle = obj->getSettingHandle(state);
      if(handle)
        set_ptr2 = handle->get();
      else {
        PRINTFB(G, FB_Executive, FB_Errors)
          " %s-Error: sele \"%s\" lacks state %d.\n", __func__, sele, state + 1
          ENDFB(G);
        ok = false;
      }
    }
  }
  if(ok) {
    switch (type) {
    case cSetting_boolean:
      {
        int value = SettingGet_b(G, set_ptr2, set_ptr1, index);
	result->type = PYMOL_RETURN_VALUE_IS_INT;
	result->int_value = value;
      }
      break;
    case cSetting_int:
      {
        int value = SettingGet_i(G, set_ptr2, set_ptr1, index);
	result->type = PYMOL_RETURN_VALUE_IS_INT;
	result->int_value = value;
      }
      break;
    case cSetting_float:
      {
        float value = SettingGet_f(G, set_ptr2, set_ptr1, index);
	result->type = PYMOL_RETURN_VALUE_IS_FLOAT;
	result->float_value = value;
      }
      break;
    case cSetting_float3:
      {
	result->type = PYMOL_RETURN_VALUE_IS_FLOAT_ARRAY;
	result->float_array = VLAlloc(float, 3);
	result->array_length = 3;
	copy3f(SettingGet<const float *>(G, set_ptr2, set_ptr1, index),
	    result->float_array);
      }
      break;
    case cSetting_color:
      {
        int value = SettingGet_color(G, set_ptr2, set_ptr1, index);
	result->type = PYMOL_RETURN_VALUE_IS_INT;
	result->int_value = value;
      }
      break;
    case cSetting_string:
      {
        OrthoLineType buffer = "";
	result->type = PYMOL_RETURN_VALUE_IS_STRING;
	result->string = mstrdup(SettingGetTextPtr(G, set_ptr2, set_ptr1, index, buffer));
      }
      break;
    default:
      break;
    }
  }
  return (ok);
}

int ExecutiveSetSettingFromString(PyMOLGlobals * G,
                                  int index, const char *value, const char *sele,
                                  int state, int quiet, int updates)
{
  CExecutive *I = G->Executive;
  SpecRec *rec = NULL;
  ObjectMolecule *obj = NULL;
  int sele1;
  ObjectMoleculeOpRec op;
  OrthoLineType value2;
  pymol::copyable_ptr<CSetting>* handle = nullptr;
  SettingName name;
  int nObj = 0;
  int ok = true;

  PRINTFD(G, FB_Executive)
    " %s: entered. sele \"%s\"\n", __func__, sele ENDFD;
  if(sele[0] == 0) {            /* global setting */
    ok = SettingSetFromString(G, NULL, index, value);
    if(ok) {
      if(!quiet) {
        if(Feedback(G, FB_Setting, FB_Actions)) {
          SettingGetTextValue(G, NULL, NULL, index, value2);
          SettingGetName(G, index, name);
          PRINTF " Setting: %s set to %s.\n", name, value2 ENDF(G);
        }
      }
      if(updates)
        SettingGenerateSideEffects(G, index, sele, state, quiet);
    }
  }
  else {
    CTracker *I_Tracker = I->Tracker;
    int list_id = ExecutiveGetNamesListFromPattern(G, sele, true, true);
    int iter_id = TrackerNewIter(I_Tracker, 0, list_id);
    while(TrackerIterNextCandInList(I_Tracker, iter_id, (TrackerRef **) (void *) &rec)) {
      if(rec) {
        switch (rec->type) {
        case cExecAll:
          rec = NULL;
          while(ListIterate(I->Spec, rec, next)) {
            if(rec->type == cExecObject) {
              {
                handle = rec->obj->getSettingHandle(state);
                if(handle) {
                  SettingCheckHandle(G, *handle);
                  ok = SettingSetFromString(G, handle->get(), index, value);
                  if(updates)
                    SettingGenerateSideEffects(G, index, rec->name, state, quiet);
                  nObj++;
                }
              }
            }
          }
          if(Feedback(G, FB_Setting, FB_Actions)) {
            if(nObj && handle) {
              SettingGetTextValue(G, handle->get(), NULL, index, value2);
              SettingGetName(G, index, name);
              if(!quiet) {
                if(state < 0) {
                  PRINTF
                    " Setting: %s set to %s in %d objects.\n", name, value2, nObj ENDF(G);
                } else {
                  PRINTF
                    " Setting: %s set to %s in %d objects, state %d.\n",
                    name, value2, nObj, state + 1 ENDF(G);
                }
              }
            }
          }
          break;
        case cExecSelection:
          /* this code has not yet been tested... */

          sele1 = SelectorIndexByName(G, rec->name);
          if(sele1 >= 0) {
            int type;
            int value_store;
            if(SettingStringToTypedValue(G, index, value, &type, &value_store)) {
              ObjectMoleculeOpRecInit(&op);
              op.code = OMOP_SetAtomicSetting;
              op.i1 = index;
              op.i2 = type;
              op.ii1 = &value_store;
              rec = NULL;
              while((ListIterate(I->Spec, rec, next))) {
                if((rec->type == cExecObject) && (rec->obj->type == cObjectMolecule)) {
                  obj = (ObjectMolecule *) rec->obj;
                  op.i4 = 0;
                  ObjectMoleculeSeleOp(obj, sele1, &op);
                  if(op.i4) {
                    if(updates)
                      SettingGenerateSideEffects(G, index, rec->name, state, quiet);
                    if(!quiet) {
                      SettingGetName(G, index, name);
                      PRINTF
                        " Setting: %s set for %d atoms in object \"%s\".\n",
                        name, op.i4, rec->obj->Name ENDF(G);
                    }
                  }
                }
              }
            }
          }
          break;
        case cExecObject:
          {
            handle = rec->obj->getSettingHandle(state);
            if(handle) {
              SettingCheckHandle(G, *handle);
              ok = SettingSetFromString(G, handle->get(), index, value);
              if(ok) {
                if(updates)
                  SettingGenerateSideEffects(G, index, sele, state, quiet);
                if(!quiet) {
                  if(state < 0) {       /* object-specific */
                    if(Feedback(G, FB_Setting, FB_Actions)) {
                      SettingGetTextValue(G, handle->get(), NULL, index, value2);
                      SettingGetName(G, index, name);
                      PRINTF
                        " Setting: %s set to %s in object \"%s\".\n",
                        name, value2, rec->obj->Name ENDF(G);
                    }
                  } else {      /* state-specific */
                    if(Feedback(G, FB_Setting, FB_Actions)) {
                      SettingGetTextValue(G, handle->get(), NULL, index, value2);
                      SettingGetName(G, index, name);
                      PRINTF
                        " Setting: %s set to %s in object \"%s\", state %d.\n",
                        name, value2, rec->obj->Name, state + 1 ENDF(G);
                    }
                  }
                }
              }
            }
          }
          break;
        }
      }
    }
    TrackerDelList(I_Tracker, list_id);
    TrackerDelIter(I_Tracker, iter_id);
  }
  return (ok);
}

int ExecutiveSetObjSettingFromString(PyMOLGlobals * G,
                                     int index, const char *value, pymol::CObject * obj,
                                     int state, int quiet, int updates)
{
  OrthoLineType value2;
  pymol::copyable_ptr<CSetting>* handle = nullptr;
  SettingName name;
  int ok = true;

  PRINTFD(G, FB_Executive)
    " ExecutiveSetObjSettingFromString: entered \n" ENDFD;
  if(!obj) {                    /* global */
    ok = SettingSetFromString(G, NULL, index, value);
    if(ok) {
      if(!quiet) {
        if(Feedback(G, FB_Setting, FB_Actions)) {
          SettingGetTextValue(G, NULL, NULL, index, value2);
          SettingGetName(G, index, name);
          PRINTF " Setting: %s set to %s.\n", name, value2 ENDF(G);
        }
      }
      if(updates)
        SettingGenerateSideEffects(G, index, obj->Name, state, quiet);
    }
  } else {                      /* based on a single object */
    {
      handle = obj->getSettingHandle(state);
      if(handle) {
        SettingCheckHandle(G, *handle);
        ok = SettingSetFromString(G, handle->get(), index, value);
        if(ok) {
          if(updates)
            SettingGenerateSideEffects(G, index, obj->Name, state, quiet);
          if(!quiet) {
            if(state < 0) {     /* object-specific */
              if(Feedback(G, FB_Setting, FB_Actions)) {
                SettingGetTextValue(G, handle->get(), NULL, index, value2);
                SettingGetName(G, index, name);
                PRINTF
                  " Setting: %s set to %s in object \"%s\".\n",
                  name, value2, obj->Name ENDF(G);
              }
            } else {            /* state-specific */
              if(Feedback(G, FB_Setting, FB_Actions)) {
                SettingGetTextValue(G, handle->get(), NULL, index, value2);
                SettingGetName(G, index, name);
                PRINTF
                  " Setting: %s set to %s in object \"%s\", state %d.\n",
                  name, value2, obj->Name, state + 1 ENDF(G);
              }
            }
          }
        }
      }
    }
  }
  return (ok);
}


/*========================================================================*/
/**
 * Restore setting to default value.
 *
 * The effective value will be the value at the next higher level (e.g. unset
 * state-level -> use object-level) or the value stored in `G->Default` for
 * global settings.
 *
 * @return Always true
 */
pymol::Result<> ExecutiveUnsetSetting(PyMOLGlobals * G, int index, pymol::zstring_view preSele,
                          int state, int quiet, int updates)
{
  CExecutive *I = G->Executive;
  SpecRec *rec = NULL;
  ObjectMolecule *obj = NULL;
  int sele1;
  ObjectMoleculeOpRec op;
  pymol::copyable_ptr<CSetting>* handle = nullptr;
  const char * name = SettingGetName(index);
  int nObj = 0;
  int ok = true;
  const char* sele = preSele.c_str();

  pymol::Result<SelectorTmp2> s1;

  if (!preSele.empty()) {
    s1 = SelectorTmp2::make(G, preSele.c_str());
    p_return_if_error(s1);
    sele = s1->getName();
  }

  if(sele[0] == 0) {
    SettingRestoreDefault(G->Setting, index, G->Default);
    if (!quiet && Feedback(G, FB_Executive, FB_Actions)) {
      OrthoLineType value = "";
      SettingGetTextValue(G, nullptr, nullptr, index, value);
      PRINTF " Setting: %s restored to default (%s)\n", name, value ENDF(G);
    }
  }
  else {
    // Undefine per-object, per-state, or per-atom settings.
    CTracker *I_Tracker = I->Tracker;
    int list_id = ExecutiveGetNamesListFromPattern(G, sele, true, true);
    int iter_id = TrackerNewIter(I_Tracker, 0, list_id);
    while(TrackerIterNextCandInList(I_Tracker, iter_id, (TrackerRef **) (void *) &rec)) {
      if(rec) {
        switch (rec->type) {
        case cExecAll:
          rec = NULL;
          while(ListIterate(I->Spec, rec, next)) {
            if(rec->type == cExecObject) {
              {
                handle = rec->obj->getSettingHandle(state);
                if (handle && *handle && SettingUnset(handle->get(), index)) {
                  nObj++;
                }
              }
            }
          }
          if(Feedback(G, FB_Setting, FB_Actions)) {
            if(nObj && handle) {
              if(!quiet) {
                if(state < 0) {
                  PRINTF " Setting: %s unset in %d objects.\n", name, nObj ENDF(G);
                } else {
                  PRINTF
                    " Setting: %s unset in %d objects, state %d.\n",
                    name, nObj, state + 1 ENDF(G);
                }
              }
            }
          }
          break;
        case cExecSelection:
          if (SettingLevelCheckMask(G, index, SettingLevelInfo[cSettingLevel_bond].mask)) {
            // handle bond-level settings (PYMOL-2726)
            ok = ExecutiveUnsetBondSetting(G, index, sele, sele, state, quiet, false);
            sele1 = -1;
          } else {
            sele1 = SelectorIndexByName(G, rec->name);
          }

          if(sele1 >= 0) {
            ObjectMoleculeOpRecInit(&op);
            op.code = OMOP_SetAtomicSetting;
            op.i1 = index;
            op.i2 = cSetting_blank;
            op.ii1 = NULL;

            rec = NULL;
            while((ListIterate(I->Spec, rec, next))) {
              if((rec->type == cExecObject) && (rec->obj->type == cObjectMolecule)) {
                obj = (ObjectMolecule *) rec->obj;
                op.i4 = 0;
                ObjectMoleculeSeleOp(obj, sele1, &op);
                if(op.i4) {
                  if(!quiet) {
                    PRINTF
                      " Setting: %s unset for %d atoms in object \"%s\".\n",
                      name, op.i4, rec->obj->Name ENDF(G);
                  }
                }
              }
            }
          }
          break;
        case cExecObject:
          {
            handle = rec->obj->getSettingHandle(state);
            if (handle && *handle && SettingUnset(handle->get(), index)) {
                if(!quiet) {
                  if(state < 0) {       /* object-specific */
                    if(Feedback(G, FB_Setting, FB_Actions)) {
                      PRINTF
                        " Setting: %s unset in object \"%s\".\n",
                        name, rec->obj->Name ENDF(G);
                    }
                  } else {      /* state-specific */
                    if(Feedback(G, FB_Setting, FB_Actions)) {
                      PRINTF
                        " Setting: %s unset in object \"%s\", state %d.\n",
                        name, rec->obj->Name, state + 1 ENDF(G);
                    }
                  }
                }
            }
          }
          break;
        }
      }
    }
    TrackerDelList(I_Tracker, list_id);
    TrackerDelIter(I_Tracker, iter_id);
  }
  if(updates)
    SettingGenerateSideEffects(G, index, sele, state, quiet);
  if (!ok) {
    return pymol::make_error("Error");
  }
  return {};
}


/*========================================================================*/
pymol::Result<> ExecutiveColorFromSele(
    PyMOLGlobals* G, const char* sele, const char* color, int flags, int quiet)
{
  auto s1 = SelectorTmp2::make(G, sele);
  p_return_if_error(s1);
  return ExecutiveColor(G, s1->getName(), color, flags, quiet);
}

pymol::Result<> ExecutiveColor(
    PyMOLGlobals* G, const char* name, const char* color, int flags, int quiet)
{
  /* flags: 
     0x1 -- ignore or suppress selection name matches
   */

  CExecutive *I = G->Executive;
  int col_ind;
  int ok = false;
  col_ind = ColorGetIndex(G, color);
  if((!name) || (!name[0]))
    name = cKeywordAll;
  if(col_ind == -1) {
    return pymol::Error("Unknown color.");
  } else {
    CTracker *I_Tracker = I->Tracker;
    SpecRec *rec = NULL;
    int n_atm = 0;
    int n_obj = 0;

    int list_id = ExecutiveGetNamesListFromPattern(G, name, true, true);
    int iter_id = TrackerNewIter(I_Tracker, 0, list_id);
    while(TrackerIterNextCandInList(I_Tracker, iter_id, (TrackerRef **) (void *) &rec)) {
      if(rec) {
        switch (rec->type) {
        case cExecSelection:
        case cExecObject:
        case cExecAll:
          if((rec->type == cExecSelection) ||   /* coloring a selection */
             (rec->type == cExecAll) || /* coloring all */
             ((rec->type == cExecObject) &&     /* coloring object and its backing selection */
              (rec->obj->type == cObjectMolecule))) {
            if(!(flags & 0x1)) {
              int sele = SelectorIndexByName(G, rec->name);
              ObjectMoleculeOpRec op;
              if(sele >= 0) {
                ok = true;
                ObjectMoleculeOpRecInit(&op);
                op.code = OMOP_COLR;
                op.i1 = col_ind;
                op.i2 = n_atm;
                ExecutiveObjMolSeleOp(G, sele, &op);
                n_atm = op.i2;
                op.code = OMOP_INVA;
                op.i1 = cRepBitmask;
                op.i2 = cRepInvColor;
                ExecutiveObjMolSeleOp(G, sele, &op);
              }
            }
          }
          break;
        }

        switch (rec->type) {    /* sets object color */
        case cExecObject:
          rec->obj->Color = col_ind;
          rec->obj->invalidate(cRepAll, cRepInvColor, -1);
          n_obj++;
          ok = true;
          SceneInvalidate(G);
          break;
        case cExecAll:
          rec = NULL;
          while(ListIterate(I->Spec, rec, next)) {
            if(rec->type == cExecObject) {
              rec->obj->Color = col_ind;
              rec->obj->invalidate(cRepAll, cRepInvColor, -1);
              n_obj++;
              ok = true;
              SceneInvalidate(G);
            }
          }
          break;
        }
      }
    }
    TrackerDelList(I_Tracker, list_id);
    TrackerDelIter(I_Tracker, iter_id);

    if(n_obj || n_atm) {
      char atms[] = "s";
      char objs[] = "s";
      if(n_obj < 2)
        objs[0] = 0;
      if(n_atm < 2)
        atms[0] = 0;
      if(!quiet) {

        if(n_obj && n_atm) {
          PRINTFB(G, FB_Executive, FB_Actions)
            " Executive: Colored %d atom%s and %d object%s.\n", n_atm, atms, n_obj, objs
            ENDFB(G);
        } else if(n_obj) {
          PRINTFB(G, FB_Executive, FB_Actions)
            " Executive: Colored %d object%s.\n", n_obj, objs ENDFB(G);
        } else {
          PRINTFB(G, FB_Executive, FB_Actions)
            " Executive: Colored %d atom%s.\n", n_atm, atms ENDFB(G);
        }
      }
    }
  }
  return {};
}


/*========================================================================*/
const char *ExecutiveFindBestNameMatch(PyMOLGlobals * G, const char *name)
{
  const char *result;
  CExecutive *I = G->Executive;
  SpecRec *rec = NULL, *best_rec = NULL;
  int best;
  int wm;
  auto ignore_case = SettingGet<bool>(G, cSetting_ignore_case);

  best = 0;
  result = name;

  while(ListIterate(I->Spec, rec, next)) {
    wm = WordMatch(G, name, rec->name, ignore_case);
    if(wm < 0) {
      best_rec = rec;
      best = wm;
      break;
    } else if((best > 0) && (best < wm)) {
      best_rec = rec;
      best = wm;
    }
  }
  if(best_rec)
    result = best_rec->name;
  return (result);
}


/*========================================================================*/
static int count_objects(PyMOLGlobals * G, int public_only)
{
  int count = 0;
  CExecutive *I = G->Executive;
  SpecRec *rec = NULL;
  while(ListIterate(I->Spec, rec, next)) {
    if(rec->type == cExecObject) {
      if(!public_only)
        count++;
      else if(rec->obj->Name[0] != '_')
        count++;
    }
  }
  return count;
}

SpecRec* ExecutiveFindSpec(PyMOLGlobals* G, pymol::zstring_view name_view)
{
  CExecutive *I = G->Executive;
  SpecRec *rec = NULL;
  const char* name = name_view.c_str();
  // ignore % prefix
  if(name[0] && name[0] == '%')
    name++;
  {                             /* first, try for perfect, case-specific match */
    OVreturn_word result;
    if(OVreturn_IS_OK((result = OVLexicon_BorrowFromCString(I->Lex, name)))) {
      if(OVreturn_IS_OK((result = OVOneToOne_GetForward(I->Key, result.word)))) {
        if(!TrackerGetCandRef(I->Tracker, result.word, (TrackerRef **) (void *) &rec)) {
          rec = NULL;
        }
      }
    }
    if(!rec) {                  /* otherwise try partial/case-nonspecific match */
      rec = ExecutiveAnyCaseNameMatch(G, name);
    }
  }
  return (rec);
}


/*========================================================================*/
bool ExecutiveObjMolSeleOp(PyMOLGlobals * G, int sele, ObjectMoleculeOpRec * op)
{
  CExecutive *I = G->Executive;
  SpecRec *rec = NULL;
  ObjectMolecule *obj = NULL;
  int update_table = true;

	/* if we're given a valid selection */
  if(sele >= 0) {
		/* iterate over all the objects in the global list */
    while(ListIterate(I->Spec, rec, next)) {
      if(rec->type == cExecObject) {
        if(rec->obj->type == cObjectMolecule) {
					/* if the objects are valid molecules, then perform the operation in op_code */
          obj = (ObjectMolecule *) rec->obj;

          switch (op->code) {
          case OMOP_RenameAtoms:
            {
              int result = SelectorRenameObjectAtoms(G, obj, sele, op->i2, update_table);
              if(result > 0)
                op->i1 += result;
              update_table = false;
            }
            break;
          default:
						/* all other cases, perform the operation on obj */
            if (!ObjectMoleculeSeleOp(obj, sele, op)) {
              return false;
            }
            break;
          }
        }
      }
    }
  }

  return true;
}


/*========================================================================*/
int ExecutiveGetCameraExtent(PyMOLGlobals * G, const char *name, float *mn, float *mx,
                             int transformed, int state)
{
  int sele;
  ObjectMoleculeOpRec op;
  int flag = false;

  if((state == -2) || (state == -3))    /* TO DO: support per-object states */
    state = SceneGetState(G);

  PRINTFD(G, FB_Executive)
    " %s: name %s state %d\n", __func__, name, state ENDFD;

  sele = SelectorIndexByName(G, name);

  if(sele >= 0) {
    ObjectMoleculeOpRecInit(&op);
    if(state < 0) {
      op.code = OMOP_CameraMinMax;
    } else {
      op.code = OMOP_CSetCameraMinMax;
      op.cs1 = state;
    }
    op.v1[0] = FLT_MAX;
    op.v1[1] = FLT_MAX;
    op.v1[2] = FLT_MAX;
    op.v2[0] = -FLT_MAX;
    op.v2[1] = -FLT_MAX;
    op.v2[2] = -FLT_MAX;
    op.i1 = 0;
    op.i2 = transformed;
    op.mat1 = SceneGetMatrix(G);

    ExecutiveObjMolSeleOp(G, sele, &op);

    PRINTFD(G, FB_Executive)
      " %s: minmax over %d vertices\n", __func__, op.i1 ENDFD;
    if(op.i1)
      flag = true;
  }
  copy3f(op.v1, mn);
  copy3f(op.v2, mx);

  PRINTFD(G, FB_Executive)
    " %s: returning %d\n", __func__, flag ENDFD;

  return (flag);
}


/*========================================================================*/
int ExecutiveGetExtent(PyMOLGlobals * G, const char *name, float *mn, float *mx,
                       int transformed, int state, int weighted)
{
  int sele;
  ObjectMoleculeOpRec op, op2;
  CExecutive *I = G->Executive;
  pymol::CObject *obj;
  int result = false;
  float f1, f2, fmx;
  int a;

  if(WordMatchExact(G, cKeywordCenter, name, true)) {
    SceneGetCenter(G, mn);
    copy3f(mn, mx);
    return 1;
  }
  if(WordMatchExact(G, cKeywordOrigin, name, true)) {
    SceneOriginGet(G, mn);
    copy3f(mn, mx);
    return 1;
  }

  PRINTFD(G, FB_Executive)
    " %s: name %s state %d\n", __func__, name, state ENDFD;

  ObjectMoleculeOpRecInit(&op);
  ObjectMoleculeOpRecInit(&op2);

  if((state == -2) || (state == -3)) {  /* we want the currently displayed state */
    state = SceneGetState(G);
    op.include_static_singletons = true;        /* make sure we get the static singletons too */
    op2.include_static_singletons = true;
  }

  op2.i1 = 0;
  op2.v1[0] = -1.0;
  op2.v1[1] = -1.0;
  op2.v1[2] = -1.0;
  op2.v2[0] = 1.0;
  op2.v2[1] = 1.0;
  op2.v2[2] = 1.0;

  {
    auto matched_recs = ExecutiveGetSpecRecsFromPattern(G, name);
    int have_atoms_flag = false;
    int have_extent_flag = false;

    /* first, compute atomic extents */

    if(weighted) {
      op2.i1 = 0;

      op2.v1[0] = 0.0F;
      op2.v1[1] = 0.0F;
      op2.v1[2] = 0.0F;

      op.i1 = 0;

      op.v1[0] = FLT_MAX;
      op.v1[1] = FLT_MAX;
      op.v1[2] = FLT_MAX;

      op.v2[0] = -FLT_MAX;
      op.v2[1] = -FLT_MAX;
      op.v2[2] = -FLT_MAX;
    }

    /* first, handle molecular objects */

    for (auto& recref : matched_recs) {
      auto* rec = &recref;
          switch (rec->type) {
          case cExecObject:
            if (rec->obj->type != cObjectMolecule &&
                rec->obj->type != cObjectAlignment)
              break;
          case cExecSelection:
          case cExecAll:
            if(rec->type == cExecAll)
              sele = SelectorIndexByName(G, cKeywordAll);
            else
              sele = SelectorIndexByName(G, rec->name);
            if(sele >= 0) {
              if(state < 0) {
                op.code = OMOP_MNMX;
              } else {
                op.code = OMOP_CSetMinMax;
                op.cs1 = state;
              }
              op.i2 = transformed;
              ExecutiveObjMolSeleOp(G, sele, &op);
              if(op.i1) {
                have_atoms_flag = true;
              }
              PRINTFD(G, FB_Executive)
                " %s: minmax over %d vertices\n", __func__, op.i1 ENDFD;
            }

            if(weighted) {
              if(state < 0)
                op2.code = OMOP_SUMC;
              else {
                op2.code = OMOP_CSetSumVertices;
                op2.cs1 = state;
              }
              op2.i2 = transformed;
              ExecutiveObjMolSeleOp(G, sele, &op2);
            }
            break;
          }
    }
    if(have_atoms_flag)
      have_extent_flag = true;

    /* now handle nonmolecular objects */

    for (auto& recref : matched_recs) {
      auto* rec = &recref;
          switch (rec->type) {
          case cExecAll:
            rec = NULL;
            while(ListIterate(I->Spec, rec, next)) {
              if(rec->type == cExecObject) {
                obj = rec->obj;
                if(!obj->ExtentFlag) {
                  switch (obj->type) {
                  case cObjectMap:
                  case cObjectMesh:
                  case cObjectSurface:
                    if(!rec->obj->ExtentFlag) {
                      /* allow object to update extents, if necessary */
                      rec->obj->update();
                    }
                  }
                }
                if(obj->ExtentFlag)
                  switch (obj->type) {
                  case cObjectMolecule:
                    break;
                    /* intentional fall-through */
                  default:
                    if(!have_extent_flag) {
                      copy3f(obj->ExtentMin, op.v1);
                      copy3f(obj->ExtentMax, op.v2);
                      have_extent_flag = true;
                    } else {
                      min3f(obj->ExtentMin, op.v1, op.v1);
                      max3f(obj->ExtentMax, op.v2, op.v2);
                    }
                    break;
                  }
              }
            }
            break;
          case cExecObject:
            obj = rec->obj;
            if(!obj->ExtentFlag) {
              switch (obj->type) {
              case cObjectMap:
              case cObjectMesh:
              case cObjectSurface:
                if(!rec->obj->ExtentFlag) {
                  /* allow object to update extents, if necessary */
                  rec->obj->update();
                }
              }
            }
            if(obj->ExtentFlag)
              switch (obj->type) {
              case cObjectMolecule:    /* will have been handled above... */
                break;
              default:
                if(!have_extent_flag) {
                  copy3f(obj->ExtentMin, op.v1);
                  copy3f(obj->ExtentMax, op.v2);
                  have_extent_flag = true;
                } else {
                  min3f(obj->ExtentMin, op.v1, op.v1);
                  max3f(obj->ExtentMax, op.v2, op.v2);
                }
                break;
              }
            break;
          }
    }

    if(have_atoms_flag && weighted) {
      if(op2.i1) {
        op2.v1[0] /= op2.i1;    /* compute average */
        op2.v1[1] /= op2.i1;
        op2.v1[2] /= op2.i1;

        for(a = 0; a < 3; a++) {        /* this puts origin at the weighted center */
          f1 = op2.v1[a] - op.v1[a];
          f2 = op.v2[a] - op2.v1[a];
          if(f1 > f2)
            fmx = f1;
          else
            fmx = f2;
          op.v1[a] = op2.v1[a] - fmx;
          op.v2[a] = op2.v1[a] + fmx;
        }
      }
    }

    if(have_extent_flag) {
      copy3f(op.v1, mn);
      copy3f(op.v2, mx);
    } else {
      zero3f(mn);
      zero3f(mx);
    }

    result = have_extent_flag;

  }

  PRINTFD(G, FB_Executive)
    " %s: returning %d\n", __func__, result ENDFD;

  return result;
}


/*========================================================================*/
static int ExecutiveGetMaxDistance(PyMOLGlobals * G, const char *name, float *pos, float *dev,
                                   int transformed, int state)
{
  int sele;
  ObjectMoleculeOpRec op, op2;
  CExecutive *I = G->Executive;
  pymol::CObject *obj;
  int flag = false;
  float f1, fmx = 0.0F;

  if((state == -2) || (state == -3))    /* TO DO: support per-object states */
    state = SceneGetState(G);

  PRINTFD(G, FB_Executive)
    " %s: name %s state %d\n", __func__, name, state ENDFD;

  ObjectMoleculeOpRecInit(&op);
  ObjectMoleculeOpRecInit(&op2);

  {
    auto matched_recs = ExecutiveGetSpecRecsFromPattern(G, name);

    op2.i1 = 0;
    op2.v1[0] = -1.0;
    op2.v1[1] = -1.0;
    op2.v1[2] = -1.0;
    op2.v2[0] = 1.0;
    op2.v2[1] = 1.0;
    op2.v2[2] = 1.0;

    {
      /* first handle molecular objects */

      for (auto& recref : matched_recs) {
        auto* rec = &recref;
          switch (rec->type) {
          case cExecObject:
          case cExecSelection:
          case cExecAll:
            if(rec->type == cExecAll)
              sele = SelectorIndexByName(G, cKeywordAll);
            else
              sele = SelectorIndexByName(G, rec->name);
            if(sele >= 0) {
              if(state < 0) {
                op.code = OMOP_MaxDistToPt;
              } else {
                op.code = OMOP_CSetMaxDistToPt;
                op.cs1 = state;
              }
              op.v1[0] = pos[0];
              op.v1[1] = pos[1];
              op.v1[2] = pos[2];
              op.i1 = 0;
              op.f1 = 0.0F;
              op.i2 = transformed;
              ExecutiveObjMolSeleOp(G, sele, &op);
              fmx = op.f1;

              if(op.i1)
                flag = true;
            }
            break;
          }
      }
    }

    {
      /* now handle nonmolecular objects */

      for (auto& recref : matched_recs) {
        auto* rec = &recref;
          switch (rec->type) {
          case cExecAll:
            rec = NULL;
            while(ListIterate(I->Spec, rec, next)) {
              if(rec->type == cExecObject) {
                obj = rec->obj;
                if(obj->ExtentFlag) {
                  switch (obj->type) {
                  case cObjectMolecule:
                    break;
                  default:
                    if(obj->ExtentFlag) {
                      f1 = (float) diff3f(obj->ExtentMin, pos);
                      if(fmx < f1)
                        fmx = f1;
                      f1 = (float) diff3f(obj->ExtentMax, pos);
                      if(fmx < f1)
                        fmx = f1;
                      flag = true;
                      break;
                    }
                  }
                }
              }
            }
            break;
          case cExecObject:
            obj = rec->obj;
            switch (rec->obj->type) {
            case cObjectMolecule:
              break;
            default:
              if(obj->ExtentFlag) {
                f1 = (float) diff3f(obj->ExtentMin, pos);
                if(fmx < f1)
                  fmx = f1;
                f1 = (float) diff3f(obj->ExtentMax, pos);
                if(fmx < f1)
                  fmx = f1;
                flag = true;
              }
              break;
            }
          }
        }
    }
  }
  *dev = fmx;
  return (flag);
}


/*========================================================================*/
pymol::Result<> ExecutiveWindowZoom(PyMOLGlobals* G,
    const char* name, float buffer, int state, int inclusive, float animate,
    int quiet)
{
  float center[3], radius;
  float mn[3], mx[3], df[3];
  int sele0;

  PRINTFD(G, FB_Executive)
    " ExecutiveWindowZoom-DEBUG: entered\n" ENDFD;
  if(ExecutiveGetExtent(G, name, mn, mx, true, state, true)) {
    if(buffer != 0.0F) {
      mx[0] += buffer;
      mx[1] += buffer;
      mx[2] += buffer;
      mn[0] -= buffer;
      mn[1] -= buffer;
      mn[2] -= buffer;
    }
    subtract3f(mx, mn, df);
    average3f(mn, mx, center);
    if(inclusive) {
      if(!ExecutiveGetMaxDistance(G, name, center, &radius, true, state))
        radius = 0.0;
      radius += buffer;
    } else {
      radius = df[0];
      if(radius < df[1])
        radius = df[1];
      if(radius < df[2])
        radius = df[2];
      radius = radius / 2.0F;
    }
    if(radius < MAX_VDW)
      radius = MAX_VDW;
    PRINTFD(G, FB_Executive)
      " %s: zooming with radius %8.3f...state %d\n", __func__, radius, state ENDFD;
    PRINTFD(G, FB_Executive)
      " %s: on center %8.3f %8.3f %8.3f...\n", __func__, center[0],
      center[1], center[2]
      ENDFD;
    if(animate < 0.0F) {
      if(SettingGetGlobal_b(G, cSetting_animation))
        animate = SettingGetGlobal_f(G, cSetting_animation_duration);
      else
        animate = 0.0F;
    }
    if(animate != 0.0F)
      ScenePrimeAnimation(G);
    SceneOriginSet(G, center, false);
    SceneWindowSphere(G, center, radius);
    if(animate != 0.0F)
      SceneLoadAnimation(G, animate, 0);
    else
      SceneAbortAnimation(G);
    SceneInvalidate(G);
  } else {

    sele0 = SelectorIndexByName(G, name);
    if(sele0 > 0) {             /* any valid selection except "all" */
      /* no longer an error to zoom on an empty selection -- just has no effect */
      if(!quiet) {
        PRINTFB(G, FB_Executive, FB_Warnings)
          "ExecutiveWindowZoom-Warning: selection doesn't specify any coordinates.\n"
          ENDFB(G);
      }
    } else if(ExecutiveValidName(G, name)) {
      PRINTFD(G, FB_Executive)
        " ExecutiveWindowZoom-DEBUG: name valid, but no extents -- using default view\n"
        ENDFD;
      SceneSetDefaultView(G);
      SceneInvalidate(G);
    } else {
      return pymol::make_error( __func__, "selection or object unknown.");
    }
  }
  return {};
}


/*========================================================================*/
pymol::Result<> ExecutiveCenter(PyMOLGlobals* G, const char* name,
    int state, int origin, float animate, float* pos, int quiet)
{
  float center[3];
  float mn[3], mx[3];
  int sele0;
  int have_center = false;

  if(name && ExecutiveGetExtent(G, name, mn, mx, true, state, true)) {
    average3f(mn, mx, center);
    have_center = true;
    PRINTFD(G, FB_Executive)
      " %s: centering state %d\n", __func__, state ENDFD;
    PRINTFD(G, FB_Executive)
      " %s: on center %8.3f %8.3f %8.3f...\n", __func__, center[0],
      center[1], center[2]
      ENDFD;
  } else if(pos) {
    have_center = true;
    copy3f(pos, center);
  }
  if(have_center) {
    if(animate < 0.0F) {
      if(SettingGetGlobal_b(G, cSetting_animation))
        animate = SettingGetGlobal_f(G, cSetting_animation_duration);
      else
        animate = 0.0F;
    }

    if(animate != 0.0F)
      ScenePrimeAnimation(G);
    if(origin)
      SceneOriginSet(G, center, false);
    SceneRelocate(G, center);
    SceneInvalidate(G);
    if(animate != 0.0F)
      SceneLoadAnimation(G, animate, 0);
  } else {
    sele0 = SelectorIndexByName(G, name);
    if(sele0 >= 0) {            /* any valid selection except "all" */
      if(!quiet) {
        /* no longer an error to center on an empty selection -- just have no effect */
        PRINTFB(G, FB_Executive, FB_Warnings)
          "ExecutiveCenter-Warning: selection doesn't specify any coordinates.\n"
          ENDFB(G);
      }
    } else if(ExecutiveValidName(G, name)) {
      SceneSetDefaultView(G);
      SceneInvalidate(G);
    } else {
      return pymol::make_error("Selection or object unknown.");
    }
  }
  return {};
}


/*========================================================================*/
pymol::Result<> ExecutiveOrigin(PyMOLGlobals* G, const char* sele, int preserve,
    const char* oname, const float* pos, int state)
{

  float center[3];
  float mn[3], mx[3];
  pymol::CObject *obj = NULL;
  int have_center = false;
  if(oname && oname[0]) {
    obj = ExecutiveFindObjectByName(G, oname);
    if(!obj)
      return pymol::make_error("Object ", oname, " not found.");
  }
  if(sele && sele[0]) {
    auto s1 = SelectorTmp2::make(G, sele);
    auto has_extent = ExecutiveGetExtent(G, s1->getName(), mn, mx, true, state, true);
    if(!has_extent) {
      return pymol::make_error("Could not determine extent of selection.");
    }
    average3f(mn, mx, center);
    have_center = true;
  } else if(pos) {
    copy3f(pos, center);
    have_center = true;
  }

  if (!have_center) {
    return pymol::make_error("Center could not be determined.");
  }

  if(obj) {
    ObjectSetTTTOrigin(obj, center);
    PRINTFB(G, FB_Executive, FB_Blather)
      " %s: origin for %s set to %8.3f %8.3f %8.3f\n", __func__,
      oname, center[0], center[1], center[2]
      ENDFB(G);
  } else {
    PRINTFB(G, FB_Executive, FB_Blather)
      " %s: scene origin set to %8.3f %8.3f %8.3f\n", __func__,
      center[0], center[1], center[2]
      ENDFB(G);
    SceneOriginSet(G, center, preserve);
  }
  SceneInvalidate(G);
  return {};
}


/*========================================================================*/
int ExecutiveGetMoment(PyMOLGlobals* G, const char* name, double* mi, int state)
{
  int sele;
  ObjectMoleculeOpRec op;
  int a, b;
  int c = 0;

  if((state == -2) || (state == -3))    /* TO DO: support per-object states */
    state = SceneGetState(G);

  sele = SelectorIndexByName(G, name);
  if(sele >= 0) {
    ObjectMoleculeOpRecInit(&op);
    if(state < 0) {
      op.code = OMOP_SUMC;
    } else {
      op.code = OMOP_CSetSumVertices;
      op.cs1 = state;
    }

    op.v1[0] = 0.0;
    op.v1[1] = 0.0;
    op.v1[2] = 0.0;
    op.i1 = 0;
    op.i2 = 0;                  /* untransformed...is this right? */

    ExecutiveObjMolSeleOp(G, sele, &op);

    if(op.i1) {                 /* any vertices? */
      c += op.i1;
      scale3f(op.v1, 1.0F / op.i1, op.v1);      /* compute raw average */
      if(state < 0) {
        op.code = OMOP_MOME;
      } else {
        op.code = OMOP_CSetMoment;
        op.cs1 = state;
      }
      for(a = 0; a < 3; a++)
        for(b = 0; b < 3; b++)
          op.d[a][b] = 0.0;
      ExecutiveObjMolSeleOp(G, sele, &op);
      {
        double *p = mi;
        for(a = 0; a < 3; a++)
          for(b = 0; b < 3; b++)
            *(p++) = op.d[a][b];
      }
    }
  } else {
    identity33d(mi);
  }

  return (c);
}


/*========================================================================*/
pymol::Result<bool>
ExecutiveSetObjVisib(PyMOLGlobals * G, pymol::zstring_view name, int onoff, int parents)
{
  CExecutive *I = G->Executive;
  bool changed = false;
  PRINTFD(G, FB_Executive)
    " ExecutiveSetObjVisib: entered.\n" ENDFD;
  {
    CTracker *I_Tracker = I->Tracker;
    SpecRec *rec;
    int list_id = ExecutiveGetNamesListFromPattern(G, name.data(), true, false);
    int iter_id = TrackerNewIter(I_Tracker, 0, list_id);
    while(TrackerIterNextCandInList(I_Tracker, iter_id, (TrackerRef **) (void *) &rec)) {

      if(rec) {
        if (!changed && rec->visible != onoff) {
          changed = true;
        }
        switch (rec->type) {
        case cExecAll:
          {
            bool const suppress_hidden =
              SettingGet<bool>(G, cSetting_suppress_hidden);
            bool const hide_underscore =
              SettingGet<bool>(G, cSetting_hide_underscore_names);
            if (suppress_hidden && hide_underscore)
              ExecutiveUpdateGroups(G, false);
            SpecRec *tRec = NULL;
            while(ListIterate(I->Spec, tRec, next)) {
              if(onoff != tRec->visible) {
                if(tRec->type == cExecObject) {
                  if(tRec->visible) {
                    tRec->in_scene = SceneObjectDel(G, tRec->obj, true);
                    ExecutiveInvalidateSceneMembers(G);
                    tRec->visible = !tRec->visible;
		    ReportEnabledChange(G, rec);
                  } else {
                    if (!(suppress_hidden && tRec->isHidden(hide_underscore))) {
                      tRec->in_scene = SceneObjectAdd(G, tRec->obj);
                      ExecutiveInvalidateSceneMembers(G);
                      tRec->visible = !tRec->visible;
		      ReportEnabledChange(G, rec);
                    }
                  }
                } else if((tRec->type != cExecSelection) || (!onoff)) {  /* hide all selections, but show all */
                  tRec->visible = !tRec->visible;
                  ReportEnabledChange(G, rec);
                }
              }
            }
          }
          break;
        case cExecObject:
          /*
             if(rec->visible!=onoff) {
             if(rec->visible) {
             rec->in_scene = SceneObjectDel(G,rec->obj);                           
             ExecutiveInvalidateSceneMembers(G);
             } else {
             rec->in_scene = SceneObjectAdd(G,rec->obj);
             ExecutiveInvalidateSceneMembers(G);
             }
             rec->visible=!rec->visible;
             }
           */
          if(onoff) {           /* enable */
            ExecutiveSpecEnable(G, rec, parents, false);
          } else {              /* disable */
            if(rec->visible) {
              if(rec->in_scene)
                rec->in_scene = SceneObjectDel(G, rec->obj, true);
              rec->visible = false;
              ExecutiveInvalidateSceneMembers(G);
	      ReportEnabledChange(G, rec);
            }
          }
          break;
        case cExecSelection:
          if(rec->visible != onoff) {
	    int previousVisible = rec->visible;
            rec->visible = !rec->visible;
            if(rec->visible)
              if(SettingGetGlobal_b(G, cSetting_active_selections)) {
                ExecutiveHideSelections(G);
                rec->visible = true;
              }
            SceneInvalidate(G);
            SeqDirty(G);
	    if (previousVisible!=rec->visible){
	      ReportEnabledChange(G, rec);
	    }
          }
          break;
        }
      }
    }
    TrackerDelList(I_Tracker, list_id);
    TrackerDelIter(I_Tracker, iter_id);
  }
  PRINTFD(G, FB_Executive)
    " ExecutiveSetObjVisib: leaving...\n" ENDFD;
  return changed;
}


/*========================================================================*/
/**
 * Full screen state fallback in case we can't get get the state from
 * the window manager.
 */
static bool _is_full_screen = false;

/*========================================================================*/
/**
 * Get the fullscreen state from the window manager or return -1 if
 * not available.
 */
bool ExecutiveIsFullScreen(PyMOLGlobals * G) {
  if(!G->HaveGUI || !G->ValidContext)
    return false;

  int flag = -1;

#if defined(GLUT_FULL_SCREEN)
  flag = glutGet(GLUT_FULL_SCREEN);
#endif

  PRINTFD(G, FB_Executive)
    " %s: flag=%d fallback=%d.\n", __func__,
    flag, _is_full_screen ENDFD;

  if (flag > -1)
    return flag;
  return _is_full_screen;
}

/*========================================================================*/
void ExecutiveFullScreen(PyMOLGlobals * G, int flag)
{
  if(!G->HaveGUI)
    return;

  int wm_flag = ExecutiveIsFullScreen(G);

  if(flag < 0) {
    flag = !wm_flag;
  }

  _is_full_screen = (flag != 0);

#ifndef _PYMOL_NO_GLUT
  if(G->HaveGUI && G->ValidContext) {
    if (flag) {
#ifndef GLUT_FULL_SCREEN
      if(wm_flag < 1) {
        CExecutive *I = G->Executive;
        I->oldPX = p_glutGet(P_GLUT_WINDOW_X);
        I->oldPY = p_glutGet(P_GLUT_WINDOW_Y);
        I->oldWidth = p_glutGet(P_GLUT_WINDOW_WIDTH);
        I->oldHeight = p_glutGet(P_GLUT_WINDOW_HEIGHT);
      }
#endif

      p_glutFullScreen();
    } else {
#ifndef GLUT_FULL_SCREEN
      // freeglut < 2.6
      CExecutive *I = G->Executive;
      p_glutReshapeWindow(I->oldWidth, I->oldHeight);
      p_glutPositionWindow(I->oldPX, I->oldPY);
#elif !defined(GLUT_HAS_MULTI)
      // freeglut < 2.8
      if(wm_flag)
        glutFullScreenToggle();
#else
      glutLeaveFullScreen();
#endif
    }
  }
#endif

  PyMOL_NeedReshape(G->PyMOL, flag, 0, 0, 0, 0);
  SceneChanged(G);
}


/*========================================================================*/
static
void fInvalidateRepMask(pymol::CObject * obj, cRepBitmask_t repmask, int state=-1) {
    for (auto a = cRep_t(0); a < cRepCnt; ++a) {
      if ((1 << a) & repmask)
        obj->invalidate(a, cRepInvVisib, state);
    }
}


/*========================================================================*/
pymol::Result<>
ExecutiveToggleRepVisib(PyMOLGlobals * G, const char *name, int rep)
{
  int sele = -1;
  SpecRec *tRec;
  ObjectMoleculeOpRec op;
  OrthoLineType tmpname;

  PRINTFD(G, FB_Executive)
    " ExecutiveToggleRepVisib: entered.\n" ENDFD;

  tRec = ExecutiveFindSpec(G, name);

  if(rep == -2) {
    // special case: toggle object visibility (should that be in this function?)
    if(tRec) {
      ExecutiveSetObjVisib(G, name, !tRec->visible, 0);
    } else {
      return pymol::make_error(name, " not found.");
    }
  } else if(tRec && tRec->type == cExecObject &&
      tRec->obj->type != cObjectMolecule) {
    // non-atom object
    tRec->obj->visRep ^= rep;
    fInvalidateRepMask(tRec->obj, rep, 0);
    SceneChanged(G);
  } else if(SelectorGetTmp(G, name, tmpname) >= 0) {
    // atom selection
    sele = SelectorIndexByName(G, tmpname);
    if(sele >= 0) {
          ObjectMoleculeOpRecInit(&op);

          op.code = OMOP_CheckVis;
          op.i1 = rep;
          op.i2 = false;
          ExecutiveObjMolSeleOp(G, sele, &op);
          op.i2 = !op.i2;

          op.code = OMOP_VISI;
          op.i1 = rep;
          ExecutiveObjMolSeleOp(G, sele, &op);
          op.code = OMOP_INVA;
          op.i2 = cRepInvVisib;
          ExecutiveObjMolSeleOp(G, sele, &op);

    }
    SelectorFreeTmp(G, tmpname);
  }
  PRINTFD(G, FB_Executive)
    " ExecutiveToggleRepVisib: leaving...\n" ENDFD;
  return {};
}


/*========================================================================*/
pymol::Result<>
ExecutiveSetRepVisib(PyMOLGlobals * G, pymol::zstring_view name, int rep, int state)
{
  int repmask = (rep == cRepAll) ? cRepBitmask : (1 << rep);
  return ExecutiveSetRepVisMask(G, name, repmask, state);
}

pymol::Result<>
ExecutiveSetRepVisMaskFromSele(PyMOLGlobals* G, pymol::zstring_view sele, int repmask, int state)
{
  if(sele[0] == '@') {
    // DEPRECATED
    sele = cKeywordAll;
    repmask = cRepBitmask;
  }
  auto s1 = SelectorTmp2::make(G, sele.c_str());
  p_return_if_error(s1);
  return ExecutiveSetRepVisMask(G, s1->getName(), repmask, state);
}

pymol::Result<>
ExecutiveSetRepVisMask(PyMOLGlobals * G, pymol::zstring_view name, int repmask, int state)
{
  PRINTFD(G, FB_Executive)
    " ExecutiveSetRepVisib: entered.\n" ENDFD;

  {
    CExecutive *I = G->Executive;
    CTracker *I_Tracker = I->Tracker;
    SpecRec *rec = NULL;
    int list_id = ExecutiveGetNamesListFromPattern(G, name.data(), true, true);
    int iter_id = TrackerNewIter(I_Tracker, 0, list_id);
    while(TrackerIterNextCandInList(I_Tracker, iter_id, (TrackerRef **) (void *) &rec)) {
      if(rec) {
        /* per-atom */

        switch (rec->type) {
        case cExecObject:
          if (rec->obj->type != cObjectMolecule &&
              rec->obj->type != cObjectAlignment)
            break;
        case cExecSelection:
          {
            int sele = SelectorIndexByName(G, rec->name);
            if(sele >= 0) {
              ObjectMoleculeOpRec op;
              ObjectMoleculeOpRecInit(&op);
              op.code = OMOP_VISI;
              op.i1 = repmask;
              op.i2 = state;
              ExecutiveObjMolSeleOp(G, sele, &op);
              op.code = OMOP_INVA;
              if (state == cVis_AS)
                op.i1 = cRepBitmask;
              op.i2 = cRepInvVisib;
              ExecutiveObjMolSeleOp(G, sele, &op);
            }
          }
          break;
        }

        /* per-object/name */

        switch (rec->type) {
        case cExecObject:
          ObjectSetRepVisMask(rec->obj, repmask, state);
          fInvalidateRepMask(rec->obj, repmask, 0);
          SceneChanged(G);
          break;
        case cExecAll:
          ExecutiveSetAllRepVisMask(G, repmask, state);
          break;
        }
      }
    }
    TrackerDelList(I_Tracker, list_id);
    TrackerDelIter(I_Tracker, iter_id);
  }
  PRINTFD(G, FB_Executive)
    " ExecutiveSetRepVisib: leaving...\n" ENDFD;
  return {};
}


/*========================================================================*/
pymol::Result<>
ExecutiveSetOnOffBySele(PyMOLGlobals * G, pymol::zstring_view sname, int onoff)
{
  SelectorTmp2 tmp{G, sname.data()};
  const char* name = tmp.getName();

  auto tRec = ExecutiveFindSpec(G, name);
  if(!tRec && sname == cKeywordAll) {
    ExecutiveSetObjVisib(G, name, onoff, false);
  }
  if(tRec) {
    int sele = tmp.getIndex();
    if(sele >= 0) {
      ObjectMoleculeOpRec op;
      ObjectMoleculeOpRecInit(&op);

      op.code = OMOP_OnOff;
      op.i1 = onoff;
      ExecutiveObjMolSeleOp(G, sele, &op);
    }
  }

  return {};
}


/*========================================================================*/
/**
 * @param repmask rep bit mask
 * @param state 0 (hide), 1 (show), 2 (as)
 */
static void ExecutiveSetAllRepVisMask(PyMOLGlobals * G, int repmask, int state)
{
  ObjectMoleculeOpRec op;
  ObjectMolecule *obj;
  int sele;

  CExecutive *I = G->Executive;
  SpecRec *rec = NULL;
  PRINTFD(G, FB_Executive)
    " ExecutiveSetAllRepVisib: entered.\n" ENDFD;
  while(ListIterate(I->Spec, rec, next)) {
    if(rec->type == cExecObject) {
      if(rec->type == cExecObject) {
        switch (rec->obj->type) {
        case cObjectMolecule:
          obj = (ObjectMolecule *) rec->obj;
          sele = SelectorIndexByName(G, obj->Name);
          ObjectMoleculeOpRecInit(&op);

          op.code = OMOP_VISI;
          op.i1 = repmask;
          op.i2 = state;
          ObjectMoleculeSeleOp(obj, sele, &op);
          op.code = OMOP_INVA;
          if (state == cVis_AS)
            op.i1 = cRepBitmask;
          op.i2 = cRepInvVisib;
          ObjectMoleculeSeleOp(obj, sele, &op);
          break;
        default:
          ObjectSetRepVisMask(rec->obj, repmask, state);
          fInvalidateRepMask(rec->obj, repmask, -1);
          SceneInvalidate(G);
          break;
        }
      }
    }
  }
  PRINTFD(G, FB_Executive)
    " ExecutiveSetAllRepVisib: leaving...\n" ENDFD;

}


/*========================================================================*/
pymol::Result<> ExecutiveInvalidateRep(
    PyMOLGlobals* G, const char* str1, cRep_t rep, cRepInv_t level)
{
  CExecutive *I = G->Executive;
  ObjectMoleculeOpRec op;
  SpecRec *rec = NULL;
  const char* name = "";

  SelectorTmp2 s1;
  if(str1 && !WordMatchExact(G, str1, "all", true)) {
    s1 = SelectorTmp2(G, str1);
    name = s1.getName();
  } else {
    name = str1;
  }

  if((!name) || (!name[0]))
    name = cKeywordAll;
  {
    CTracker *I_Tracker = I->Tracker;
    int list_id = ExecutiveGetNamesListFromPattern(G, name, true, true);
    int iter_id = TrackerNewIter(I_Tracker, 0, list_id);

    while(TrackerIterNextCandInList(I_Tracker, iter_id, (TrackerRef **) (void *) &rec)) {
      if(rec) {
        switch (rec->type) {
        case cExecSelection:
        case cExecObject:
          {
            int sele = SelectorIndexByName(G, rec->name);
            if(sele >= 0) {
              ObjectMoleculeOpRecInit(&op);
              op.code = OMOP_INVA;
              op.i1 = (rep == cRepAll) ? cRepBitmask : (1 << rep);
              op.i2 = level;
              ExecutiveObjMolSeleOp(G, sele, &op);
            } else {
              rec->obj->invalidate(rep, level, -1);
            }
          }
          break;
        case cExecAll:
          rec = NULL;
          while(ListIterate(I->Spec, rec, next)) {
            if(rec->type == cExecObject) {
              rec->obj->invalidate(rep, level, -1);
            }
          }
          SceneInvalidate(G);
          break;
        }
      }
    }
    TrackerDelList(I_Tracker, list_id);
    TrackerDelIter(I_Tracker, iter_id);
  }
  return {};
}

int ExecutiveCheckGroupMembership(PyMOLGlobals * G, int list_id, pymol::CObject * obj)
{
  CExecutive *I = G->Executive;
  int result = false;
  CTracker *I_Tracker = I->Tracker;
  int iter_id = TrackerNewIter(I_Tracker, 0, list_id);
  if(iter_id) {
    SpecRec *rec = NULL;
    while(TrackerIterNextCandInList(I_Tracker, iter_id, (TrackerRef **) (void *) &rec)) {
      if(rec && (rec->type == cExecObject) && (rec->obj == obj)) {
        result = true;
        break;
      }
    }
    TrackerDelIter(I_Tracker, iter_id);
  }
  return result;
}

int ExecutiveGetExpandedGroupListFromPattern(PyMOLGlobals * G, const char *name)
{
  CExecutive *I = G->Executive;
  int result = 0;
  CWordMatcher *matcher;
  CWordMatchOptions options;
  CTracker *I_Tracker = I->Tracker;
  const char *wildcard = SettingGetGlobal_s(G, cSetting_wildcard);
  int iter_id = TrackerNewIter(I_Tracker, 0, I->all_names_list_id);
  int cand_id;
  SpecRec *rec;

  WordMatchOptionsConfigNameList(&options,
                                 *wildcard, SettingGetGlobal_b(G, cSetting_ignore_case));
  matcher = WordMatcherNew(G, name, &options, false);
  if(matcher) {
    if(iter_id) {
      while((cand_id = TrackerIterNextCandInList(I_Tracker, iter_id,
                                                 (TrackerRef **) (void *) &rec))) {
        if(rec && !(rec->type == cExecAll)) {
          if(WordMatcherMatchAlpha(matcher, rec->name)) {
            if((rec->type == cExecObject) && (rec->obj->type == cObjectGroup)) {
              if(!result)
                result = TrackerNewList(I_Tracker, NULL);
              if(result) {
                TrackerLink(I_Tracker, cand_id, result, 1);
              }
            }
          }
        }
      }
    }
  } else if((rec = ExecutiveFindSpec(G, name))) {       /* only one name in list */
    if((rec->type == cExecObject) && (rec->obj->type == cObjectGroup)) {
      result = TrackerNewList(I_Tracker, NULL);
      TrackerLink(I_Tracker, rec->cand_id, result, 1);
    }
  } else if((rec = ExecutiveUnambiguousNameMatch(G, name))) {
    if((rec->type == cExecObject) && (rec->obj->type == cObjectGroup)) {
      result = TrackerNewList(I_Tracker, NULL);
      TrackerLink(I_Tracker, rec->cand_id, result, 1);
    }
  }
  if(matcher)
    WordMatcherFree(matcher);
  if(iter_id)
    TrackerDelIter(I->Tracker, iter_id);
  if(result)
    ExecutiveExpandGroupsInList(G, result, cExecExpandGroups);
  return result;
}

int ExecutiveGetExpandedGroupList(PyMOLGlobals * G, const char *name)
{
  CExecutive *I = G->Executive;
  int result = 0;
  int list_id = 0;

  SpecRec *rec = ExecutiveFindSpec(G, name);
  ExecutiveUpdateGroups(G, false);
  if(rec && (rec->type == cExecObject) && (rec->obj->type == cObjectGroup)) {
    list_id = rec->group_member_list_id;
  }
  if(list_id) {
    result = TrackerNewListCopy(I->Tracker, list_id, NULL);
    ExecutiveExpandGroupsInList(G, result, cExecExpandGroups);
  }
  return result;
}

void ExecutiveFreeGroupList(PyMOLGlobals * G, int list_id)
{
  CExecutive *I = G->Executive;
  TrackerDelList(I->Tracker, list_id);
}


/*========================================================================*/
pymol::CObject *ExecutiveFindObjectByName(PyMOLGlobals * G, pymol::zstring_view name)
{
  pymol::CObject *obj = NULL;
  SpecRec *rec = ExecutiveFindSpec(G, name.c_str());
  if(rec && (rec->type == cExecObject)) {
    obj = rec->obj;
  }
  return (obj);
}

/*========================================================================*/
/* returns NULL if none found */
pymol::CObject ** ExecutiveFindObjectsByType(PyMOLGlobals * G, int objType) {
  CExecutive *I = G->Executive;
  SpecRec *rec = NULL;
  int n = 0;
  pymol::CObject** rVal = VLAlloc(pymol::CObject*, 1);

  /* loop over all known objects */
  while(ListIterate(I->Spec,rec,next)) {
    /* make sure it exists and is the right type */
    if (rec->obj && rec->type==cExecObject) {
      if(rec->obj->type==objType) {
	/* this could be optimized */
	VLACheck(rVal, pymol::CObject*, n);
	*(rVal+n) = rec->obj;
	n++;
      }
    }
  }
  VLASize(rVal, pymol::CObject*, n);
  if (n==0) {
	VLAFree(rVal);
    return NULL;
  }
  else
    return rVal;
}

/*========================================================================*/
Block *ExecutiveGetBlock(PyMOLGlobals * G)
{
  CExecutive *I = G->Executive;
  return (I);
}


/**
 * Makes a custom symmetry operation label.
 */
static std::string make_symexp_segi_label(int a, int x, int y, int z)
{
  SegIdent seg = "_000";

  if (a > 61) {
    // beyond what can be encoded with a single alphanumeric
    // character (PYMOL-2475)
  } else if (a > 35) {
    seg[0] = 'a' + (a - 36);
  } else if (a > 25) {
    seg[0] = '0' + (a - 26);
  } else {
    seg[0] = 'A' + a;
  }

  if (x > 0) {
    seg[1] = 'A' + x - 1;
  } else if (x < 0) {
    seg[1] = 'Z' + x + 1;
  }

  if (y > 0) {
    seg[2] = 'A' + y - 1;
  } else if (y < 0) {
    seg[2] = 'Z' + y + 1;
  }

  if (z > 0) {
    seg[3] = 'A' + z - 1;
  } else if (z < 0) {
    seg[3] = 'Z' + z + 1;
  }

  return seg;
}

/*========================================================================*/
/**
 * Generates symmetry mates of object @a name as new objects, with all states.
 *
 * Only includes mates within @a cutoff of selection @a s1, and only in a range
 * of +/-1 unit cells. The cutoff applies across states, so a symmetry atom will
 * be included if it's within any atom of @a s1 in any state.
 *
 * @param name Prefix for new objects
 * @param oname Object to replicate
 * @param s1 Atom selection (for cutoff)
 * @param cutoff Distance cutoff from selection
 * @param segi (bool) If true, write symmetry operation code to segment
 * identifier
 */
void ExecutiveSymExp(PyMOLGlobals * G, const char *name,
                     const char *oname, const char *s1, float cutoff, int segi, int quiet)
{                               /* TODO state */
  SelectorTmp tmpsele1(G, s1);
  auto sele = tmpsele1.getIndex();

	/* object to expand */
  auto const* const obj = ExecutiveFindObject<ObjectMolecule>(G, oname);

  if(!(obj && sele)) {
    ErrMessage(G, __func__, "Invalid object");
    return;
  }

  // Get symmetry from first coordset
  const CSymmetry* sym = nullptr;
  for (int b = 0; b < obj->NCSet && !sym; ++b) {
    if (auto const* cs = obj->CSet[b]) {
      sym = cs->getSymmetry();
    }
  }

  if (!sym) {
    ErrMessage(G, __func__, "No symmetry loaded!");
    return;
  }

  int const nsymmat = sym->getNSymMat();

  if (nsymmat < 1) {
    ErrMessage(G, __func__, "unknown space group");
    return;
  }

  if (!quiet) {
      PRINTFB(G, FB_Executive, FB_Actions)
        " ExecutiveSymExp: Generating symmetry mates...\n" ENDFB(G);
  }

  bool const matrix_mode = SettingGet<int>(G, cSetting_matrix_mode) > 0;

  // Get the center of mass for this object/selection (over all states!)
  ObjectMoleculeOpRec op;
  ObjectMoleculeOpRecInit(&op);
  op.code = OMOP_SUMC;
  ExecutiveObjMolSeleOp(G, sele, &op);
  auto tc_vec = glm::make_vec3(op.v1) / float(std::max(1, op.i1));
  auto tc = glm::value_ptr(tc_vec);

  /* Transformation: RealToFrac (3x3) * tc (3x1) = (3x1) */
  transform33f3f(sym->Crystal.realToFrac(), tc, tc);

  /* 2.  Copy the coordinates for the atoms in this selection into op */
  op.code = OMOP_VERT;
  op.nvv1 = 0;
  op.vv1 = VLAlloc(float, 10000);
  ExecutiveObjMolSeleOp(G, sele, &op);
  auto const vv1 = pymol::vla_take_ownership(op.vv1);

  /* op.nvv1 is the number of atom coordinates we copied in the previous step */
  if (!op.nvv1) {
    ErrMessage(G, __func__, "No atoms indicated!");
    return;
  }

  auto map =
      std::unique_ptr<MapType>(MapNew(G, -cutoff, vv1.data(), op.nvv1, NULL));
  if (!map) {
    ErrMessage(G, __func__, "map is NULL");
    return;
  }

  std::vector<glm::vec3> cs_centers(obj->NCSet);
  for (int b = 0; b < obj->NCSet; ++b) {
    if (auto const* cs = obj->CSet[b]) {
      CoordSetGetAverage(cs, glm::value_ptr(cs_centers[b]));
    }
  }

  /* go out no more than one lattice step in each direction: -1, 0, +1 */
  for (int x = -1; x < 2; ++x) {
    for (int y = -1; y < 2; ++y) {
      for (int z = -1; z < 2; ++z) {
        for (int a = 0; a < nsymmat; a++) {
          /* make a copy of the original */
          auto new_obj = ObjectMoleculeCopy(obj);
          bool keepFlag = false;

          for (int b = 0; b < new_obj->NCSet; ++b) {
            auto* cs = new_obj->CSet[b];
            if (!cs) {
              continue;
            }

            float mat[16];
            copy33f44f(sym->Crystal.realToFrac(), mat);
            left_multiply44f44f(sym->getSymMat(a), mat);

            float ts[3];
            transform44f3f(mat, glm::value_ptr(cs_centers[b]), ts);

            /* compute the effective translation resulting
               from application of the symmetry operator so
               that we can shift it into the cell of the
               target selection */
            for (int c = 0; c < 3; c++) {
              ts[c] = std::round(tc[c] - ts[c]);
            }
            float const shift[3] = {ts[0] + x, ts[1] + y, ts[2] + z};
            float m[16];
            identity44f(m);
            m[3] = shift[0];
            m[7] = shift[1];
            m[11] = shift[2];

            left_multiply44f44f(const_cast<float const*>(m), mat);
            copy33f44f(sym->Crystal.fracToReal(), m);
            left_multiply44f44f(const_cast<float const*>(m), mat);

            if (is_identityf(4, mat)) {
              continue;
            }

            double mat_d[16];
            copy44f44d(mat, mat_d);
            ObjectStateLeftCombineMatrixR44d(cs, mat_d);

            if (!matrix_mode) {
              CoordSetTransform44f(cs, mat);
            }

            if (keepFlag) {
              continue;
            }

            /* for each coordinate in this coordinate set */
            for (unsigned idx = 0; idx < cs->NIndex; ++idx) {
              const auto* v2 = cs->coordPtr(idx);

              if (matrix_mode) {
                transform44f3f(mat, v2, ts);
                v2 = ts;
              }

              if (MapAnyWithin(*map, vv1.data(), v2, cutoff)) {
                keepFlag = true;
                break;
              }
            }

            // CIF-style symop label (e.g. "1_555")
            cs->setTitle(pymol::string_format("%d_%.0f%.0f%.0f", a + 1,
                shift[0] + 5, shift[1] + 5, shift[2] + 5));
          }

          if (!keepFlag) {
            DeleteP(new_obj);
            continue;
          }

          if (segi) {
            auto seg = make_symexp_segi_label(a, x, y, z);
            lexidx_t segi = LexIdx(G, seg.c_str());
            for (unsigned atm = 0; atm < new_obj->NAtom; ++atm) {
              LexAssign(G, new_obj->AtomInfo[atm].segi, segi);
            }
            LexDec(G, segi);
          }

          // manage the new object

          auto const new_name =
              pymol::string_format("%s%02d%02d%02d%02d", name, a, x, y, z);

          PRINTFB(G, FB_Executive, FB_Blather)
          "Making new object: %s from name=%s, a=%d, x=%d, y=%d, z=%d\n",
              new_name.c_str(), name, a, x, y, z ENDFB(G);

          ObjectSetName(new_obj, new_name.c_str());
          ExecutiveDelete(G, new_obj->Name);
          ExecutiveManageObject(G, new_obj, false, quiet);
        }
      }
    }
  }
}

void ExecutivePurgeSpec(PyMOLGlobals * G, SpecRec * rec, bool save)
{
  CExecutive *I = G->Executive;

  if (!save) {
    CGOFree(rec->gridSlotSelIndicatorsCGO);
  }

  ExecutiveInvalidateGroups(G, false);
  ExecutiveInvalidatePanelList(G);
  switch (rec->type) {
  case cExecObject:
    if(I->LastEdited == rec->obj)
      I->LastEdited = NULL;
    if(rec->obj->type == cObjectMolecule)
      if(EditorIsAnActiveObject(G, (ObjectMolecule *) rec->obj))
        EditorInactivate(G);
    SeqChanged(G);
    if(rec->visible) {
      SceneObjectDel(G, rec->obj, false);
      ExecutiveInvalidateSceneMembers(G);
    }
    ExecutiveDelKey(I, rec);
    SelectorDelete(G, rec->name);
    if (!save) {
      DeleteP(rec->obj);
    }
    TrackerDelCand(I->Tracker, rec->cand_id);
    break;
  case cExecSelection:
    if(rec->visible) {
      SceneInvalidate(G);
      SeqDirty(G);
    }
    ExecutiveDelKey(I, rec);
    SelectorDelete(G, rec->name);
    TrackerDelCand(I->Tracker, rec->cand_id);
    break;
  }
}


/*========================================================================*/
pymol::Result<std::vector<DiscardedRec>> ExecutiveDelete(PyMOLGlobals * G, pymol::zstring_view nameView, bool save)
{
  std::vector<DiscardedRec> discardedRecs;
  CExecutive *I = G->Executive;
  auto name = nameView.c_str();
  std::vector<OrderRec> specPositions;
  if (save) {
    specPositions = ExecutiveGetOrderOf(G, nameView);
  }

  auto getRecPos = [&] (SpecRec* rec) {
    auto it = std::find_if(specPositions.begin(), specPositions.end(),
                [&](OrderRec& oRec) {
                  return oRec.name == rec->name;
                });
    return it == specPositions.end() ? -1 : it->pos;
  };
  auto deleteObjRec = [&] (SpecRec* rec) {
    if(save) {
      if(rec->obj->type == cObjectGroup) {
        ExecutiveGroupPurge(G, rec, &discardedRecs);
      }
      ExecutivePurgeSpec(G, rec, save);
      auto rec_pos = getRecPos(rec);
      auto rec_out = ListDetachT(I->Spec, rec);
      SceneObjectRemove(G, rec->obj);
      assert(rec_pos);
      discardedRecs.emplace_back(rec_out, rec_pos);
    } else {
      if(rec->obj->type == cObjectGroup) {
        /* remove members of the group */
        ExecutiveGroup(G, rec->name, "", cExecutiveGroupPurge, true);
      }
      ExecutivePurgeSpec(G, rec, save);
      ListDelete(I->Spec, rec, next, SpecRec);        /* NOTE: order N in list length! TO FIX */
    }
  };

  auto deleteSeleRec = [&] (SpecRec* rec) {
    ExecutivePurgeSpec(G, rec, save);
    if (save) {
      auto rec_pos = getRecPos(rec);
      auto rec_out = ListDetachT(I->Spec, rec);
      assert(rec_pos);
      discardedRecs.emplace_back(rec_out, rec_pos);
    } else {
      ListDelete(I->Spec, rec, next, SpecRec);        /* NOTE: order N in list length! TO FIX */
    }
  };

  SpecRec *rec = NULL;
  CTracker *I_Tracker = I->Tracker;
  int list_id = ExecutiveGetNamesListFromPattern(G, name, false, false);
  int iter_id = TrackerNewIter(I_Tracker, 0, list_id);
  while(TrackerIterNextCandInList(I_Tracker, iter_id, (TrackerRef **) (void *) &rec)) {
    if(rec) {
      switch (rec->type) {
      case cExecSelection:
        deleteSeleRec(rec);
        break;
      case cExecObject:
        deleteObjRec(rec);
        break;
      case cExecAll:
        rec = nullptr;
        while(ListIterate(I->Spec, rec, next)) {
          switch (rec->type) {
          case cExecAll:
            break;
          case cExecObject:
            deleteObjRec(rec);
            rec = nullptr;
            break;
          case cExecSelection:
            deleteSeleRec(rec);
            rec = nullptr;
            break;
          }
        }
        SelectorDefragment(G);
        break;
      }
    }
  }
  TrackerDelList(I_Tracker, list_id);
  TrackerDelIter(I_Tracker, iter_id);

  // fix for PYMOL-757 - Note: This should probably go somewhere else, but we
  // couldn't figure out the correct place so far
  ExecutiveUpdateGroups(G, false);
  return discardedRecs;
}

void ExecutiveReAddSpec(PyMOLGlobals* G, std::vector<DiscardedRec>& specs)
{
  auto I = G->Executive;
  for(auto& spec : specs){
    auto rec = spec.rec;
    auto pos = spec.pos;
    rec->cand_id = TrackerNewCand(I->Tracker, (TrackerRef *) rec);
    TrackerLink(I->Tracker, rec->cand_id, I->all_names_list_id, 1);
    TrackerLink(I->Tracker, rec->cand_id, I->all_obj_list_id, 1);
    ListInsertAt(I->Spec, rec, pos);
    ExecutiveAddKey(I, rec);
    ExecutiveInvalidatePanelList(G);
    if(rec->type == cExecObject) {
      rec->in_scene = SceneObjectAdd(G, rec->obj);
    }
    ExecutiveInvalidateSceneMembers(G);
    ExecutiveUpdateGroups(G, true);
  }
  specs.clear();
}


/*========================================================================*/
void ExecutiveDump(PyMOLGlobals * G, const char *fname, const char *obj, int state, int quiet)
{
  SpecRec *rec = NULL;
  CExecutive *I = G->Executive;

  SceneUpdate(G, false);

  while(ListIterate(I->Spec, rec, next)) {
    if(rec->type == cExecObject) {
      if(strcmp(rec->obj->Name, obj) == 0)
        break;
    }
  }
  if(rec) {
    if(rec->obj->type == cObjectMesh) {
      ObjectMeshDump((ObjectMesh *) rec->obj, fname, state, quiet);
    } else if(rec->obj->type == cObjectSurface) {
      ObjectSurfaceDump((ObjectSurface *) rec->obj, fname, state, quiet);
    } else if (rec->obj->type == cObjectMap) {
      ObjectMapDump((ObjectMap *) rec->obj, fname, state, quiet);
    } else {
      ErrMessage(G, __func__, "Invalid object type for this operation.");
    }
  } else {
    ErrMessage(G, __func__, "Object not found.");
  }
}


/*========================================================================*/
void ExecutiveMemoryDump(PyMOLGlobals * G)
{
  CExecutive *I = G->Executive;
  fprintf(stderr, " Executive: %d candidate(s) %d list(s) %d link(s).\n",
          TrackerGetNCand(I->Tracker),
          TrackerGetNList(I->Tracker), TrackerGetNLink(I->Tracker));
}

/**
 * The `zoom` argument takes the following values:
 *
 *   -1 = use `auto_zoom` setting
 *    0 = do nothing
 *    1 = zoom `obj` if `is_new` is ture, otherwise do nothing
 *    2 = zoom `obj`
 *    3 = zoom current state of `obj`
 *    4 = zoom all
 *    5 = zoom `obj` if it is the only object
 *
 */
void ExecutiveDoZoom(PyMOLGlobals * G, pymol::CObject * obj, int is_new, int zoom, int quiet)
{
  if (zoom) {
    if(zoom < 0) {
      zoom = SettingGetGlobal_i(G, cSetting_auto_zoom);
      if(zoom < 0) {
        zoom = 1;
      }
    }
    switch (zoom) {
    case 1:                    /* zoom when new */
      if(is_new)
        ExecutiveWindowZoom(G, obj->Name, 0.0, -1, 0, 0, quiet);       
      /* (all states) */
      break;
    case 2:                    /* zoom always */
      ExecutiveWindowZoom(G, obj->Name, 0.0, -1, 0, 0, quiet);  
      /* (all states) */
      break;
    case 3:                    /* always zoom current state */
      ExecutiveWindowZoom(G, obj->Name, 0.0, ObjectGetCurrentState(obj, false), 0, 0, quiet);  
      break;
    case 4:                    /* zoom all objects */
      ExecutiveWindowZoom(G, cKeywordAll, 0.0, -1, 0, 0, quiet);
      break;
    case 5:                    /* zoom first object only */
      if(count_objects(G, true) == 1)
        ExecutiveWindowZoom(G, obj->Name, 0.0, -1, 0, 0, quiet);        /* (all states) */
      break;
    }
  }
}

static void ExecutiveDoAutoGroup(PyMOLGlobals * G, SpecRec * rec)
{
  CExecutive *I = G->Executive;
  int auto_mode = SettingGetGlobal_i(G, cSetting_group_auto_mode);
  if(auto_mode && (rec->name[0] != '_')) {
    auto ignore_case = SettingGet<bool>(G, cSetting_ignore_case);
    char *period = rec->name + strlen(rec->name);
    SpecRec *found_group = NULL;
    WordType seek_group_name;
    UtilNCopy(seek_group_name, rec->name, sizeof(WordType));

    while((period > rec->name) && (!found_group)) {
      period--;
      if(*period == '.') {
        seek_group_name[period - rec->name] = 0;
        {
          SpecRec *group_rec = NULL;
          while(ListIterate(I->Spec, group_rec, next)) {
            if((group_rec->type == cExecObject) && (group_rec->obj->type == cObjectGroup)) {
              if(WordMatchExact(G, group_rec->name, seek_group_name, ignore_case)) {
                found_group = group_rec;
                strcpy(rec->group_name, seek_group_name);
                break;
              }
            }
          }
        }

        if((!found_group) && (auto_mode == 2)) {
          auto obj = new ObjectGroup(G);
          if(obj) {
            ObjectSetName(obj, seek_group_name);
            strcpy(rec->group_name, obj->Name);
            ExecutiveManageObject(G, obj, false, true);
            ExecutiveInvalidateGroups(G, false);
            break;
          }
        }
      }
    }
    if(found_group)
      ExecutiveInvalidateGroups(G, false);
  }
}


/*========================================================================*/
/**
 * Manages an object. Adds it to the list of spec records and related trackers,
 * and to the scene for rendering.
 *
 * Also handles:
 * - auto_dss
 * - group_auto_mode
 * - auto_defer_builds
 *
 * If the object is already managed, then only do `auto_dss` and
 * `auto_defer_builds`.
 *
 * If an object with the same name exists, then delete it and re-use the
 * existing spec rec to manage the new object.
 *
 * @param obj Object to manage. Executive takes ownership.
 * @param zoom Zoom the camera, see ExecutiveDoZoom for valid values.
 */
void ExecutiveManageObject(PyMOLGlobals * G, pymol::CObject * obj, int zoom, int quiet)
{
  SpecRec *rec = NULL;
  CExecutive *I = G->Executive;
  int exists = false;
  int previousVisible;
  int previousObjType = 0;

  if(SettingGetGlobal_b(G, cSetting_auto_hide_selections))
    ExecutiveHideSelections(G);
  while(ListIterate(I->Spec, rec, next)) {
    if(rec->obj == obj) {
      exists = true;
    }
  }
  if(!exists) {
    if (WordMatchExact(G, cKeywordAll, obj->Name, true)) {
      PRINTFB(G, FB_Executive, FB_Warnings)
        " Executive: object name \"%s\" is illegal -- renamed to 'all_'.\n", obj->Name
        ENDFB(G);
      strcat(obj->Name, "_");   /* don't allow object named "all" */
    } else if (SelectorNameIsKeyword(G, obj->Name)) {
      PRINTFB(G, FB_Executive, FB_Warnings)
        " Executive-Warning: name \"%s\" collides with a selection language keyword.\n",
        obj->Name ENDFB(G);
    }

    while(ListIterate(I->Spec, rec, next)) {
      if(rec->type == cExecObject) {
        if(strcmp(rec->obj->Name, obj->Name) == 0)
          break;
      }
    }
    if(rec) {                   /* another object of this type already exists */
      /* purge it */
      SceneObjectDel(G, rec->obj, false);
      ExecutiveInvalidateSceneMembers(G);
      previousObjType = rec->obj->type;
      DeleteP(rec->obj);
    } else {
      if(!quiet)
        if(obj->Name[0] != '_') {       /* suppress internal objects */
          PRINTFB(G, FB_Executive, FB_Actions)
            " Executive: object \"%s\" created.\n", obj->Name ENDFB(G);
        }
    }
    if(!rec)
      ListElemCalloc(G, rec, SpecRec);

    strcpy(rec->name, obj->Name);
    rec->type = cExecObject;
    rec->obj = obj;
    previousVisible = rec->visible;
    if (previousObjType == rec->obj->type) {
      // skip
    } else if (rec->obj->type == cObjectMap) {
      rec->visible = 0;
    } else {
      rec->visible = 1;
      /*      SceneObjectAdd(G,obj); */
    }
    if (previousVisible!=rec->visible){
      ReportEnabledChange(G, rec);
    }

    // Is this a new rec?
    if (!rec->cand_id) {
      rec->cand_id = TrackerNewCand(I->Tracker, (TrackerRef *) (void *) rec);

      TrackerLink(I->Tracker, rec->cand_id, I->all_names_list_id, 1);
      TrackerLink(I->Tracker, rec->cand_id, I->all_obj_list_id, 1);
      ListAppend(I->Spec, rec, next, SpecRec);
      ExecutiveAddKey(I, rec);
      ExecutiveInvalidatePanelList(G);

      ExecutiveDoAutoGroup(G, rec);
    }

    if(rec->visible) {
      rec->in_scene = SceneObjectAdd(G, obj);
      ExecutiveInvalidateSceneMembers(G);
    }
  }

  ExecutiveUpdateObjectSelection(G, obj);

  if(SettingGetGlobal_b(G, cSetting_auto_dss)) {
    if(obj->type == cObjectMolecule) {
      ObjectMolecule *objMol = (ObjectMolecule *) obj;
      if(objMol->NCSet == 1) {
        ExecutiveAssignSS(G, obj->Name, 0, NULL, true, objMol, true);
      }
    }
  }

  {
    int n_state = obj->getNFrame();
    int defer_limit = SettingGetGlobal_i(G, cSetting_auto_defer_builds);
    if((defer_limit >= 0) && (n_state >= defer_limit)) {
      int defer_builds = SettingGetGlobal_b(G, cSetting_defer_builds_mode);
      if(!defer_builds)
        SettingSetGlobal_b(G, cSetting_defer_builds_mode, 3);
    }
  }

  ExecutiveDoZoom(G, obj, !exists, zoom, true);

  SeqChanged(G);
  OrthoInvalidateDoDraw(G);
}


/*========================================================================*/
void ExecutiveManageSelection(PyMOLGlobals * G, const char *name)
{

  SpecRec* rec = nullptr;
  CExecutive *I = G->Executive;
  bool hidden_name = name[0] == '_';
  bool hide_all =
      !hidden_name && (SettingGet<bool>(G, cSetting_active_selections) ||
                       SettingGet<bool>(G, cSetting_auto_hide_selections));

  for (auto& rec_ : pymol::make_list_adapter(I->Spec)) {
    if (rec_.type == cExecSelection) {
      if (!rec && pymol::zstring_view(rec_.name) == name) {
        rec = &rec_;
      } else if (hide_all) {
        rec_.setEnabled(G, false);
      }
    }
  }

  if(!rec) {
    ListElemCalloc(G, rec, SpecRec);
    strcpy(rec->name, name);
    rec->type = cExecSelection;
    rec->next = NULL;
    rec->sele_color = -1;
    assert(!rec->visible);
    rec->cand_id = TrackerNewCand(I->Tracker, (TrackerRef *) (void *) rec);
    TrackerLink(I->Tracker, rec->cand_id, I->all_names_list_id, 1);
    TrackerLink(I->Tracker, rec->cand_id, I->all_sel_list_id, 1);
    ListAppend(I->Spec, rec, next, SpecRec);
    ExecutiveAddKey(I, rec);
    ExecutiveInvalidatePanelList(G);
  }
  assert(rec);
  if(!hidden_name) {
    if(SettingGetGlobal_b(G, cSetting_auto_show_selections) && !rec->visible) {
      rec->visible = true;
  ReportEnabledChange(G, rec);
    }
  }
  if(rec->visible)
    SceneInvalidate(G);
  ExecutiveDoAutoGroup(G, rec);
  SeqDirty(G);
}


/*========================================================================*/
int CExecutive::click(int button, int x, int y, int mod)
{
  PyMOLGlobals *G = m_G;
  CExecutive *I = G->Executive;
  int n, a;
  SpecRec *rec = NULL;
  int t, xx;
  int pass = false;
  int skip;
  int ExecLineHeight = DIP2PIXEL(SettingGetGlobal_i(G, cSetting_internal_gui_control_size));
#ifndef NDEBUG
  int hide_underscore = SettingGetGlobal_b(G, cSetting_hide_underscore_names);
#endif
  int op_cnt = get_op_cnt(G);

  if(y < I->HowFarDown) {
    if(SettingGet<InternalGUIMode>(cSetting_internal_gui_mode, G->Setting) != InternalGUIMode::Default)
      return SceneGetBlock(G)->click(button, x, y, mod);
  }

  switch(button) {
    case P_GLUT_BUTTON_SCROLL_FORWARD:
      I->m_ScrollBar.moveBy(-1);
      return 1;
    case P_GLUT_BUTTON_SCROLL_BACKWARD:
      I->m_ScrollBar.moveBy(1);
      return 1;
  }

  n = ((I->rect.top - y) - (ExecTopMargin + ExecClickMargin)) / ExecLineHeight;
  a = n;
  xx = (x - rect.left);
  if(I->ScrollBarActive) {
    if((x - rect.left) <
       (ExecScrollBarWidth + ExecScrollBarMargin + ExecToggleMargin)) {
      pass = 1;
      I->m_ScrollBar.click(button, x, y, mod);
    }
    xx -= (ExecScrollBarWidth + ExecScrollBarMargin);
  }
  skip = I->NSkip;
  if(!pass) {
    I->RecoverPressed = NULL;
    /* while(ListIterate(I->Spec,rec,next)) { */
    for (auto& panelitem : I->Panel) {
      auto* const panel = &panelitem;
      rec = panel->spec;

      assert(rec->name[0] != '_' || !hide_underscore);
      {
        if(skip) {
          skip--;
        } else {
          if(!a) {
            t = ((rect.right - ExecRightMargin) - x - 1) / ExecToggleWidth;
            if(t < op_cnt) {
              int my = rect.top - (ExecTopMargin + n * ExecLineHeight) - 3;
              int mx = rect.right - (ExecRightMargin + t * ExecToggleWidth);

#if 0
              // prefix name with % as guarantee to not match a selection keyword
              WordType namesele = "%";
              UtilNCopy(namesele + 1, (rec->type == cExecObject) ?
                  rec->obj->Name : rec->name, WordLength - 1);
#else
              char *namesele = (rec->type == cExecObject) ? rec->obj->Name : rec->name;
#endif

              t = (op_cnt - t) - 1;
              switch (t) {
              case 0:          /* action */
                switch (rec->type) {
                case cExecAll:
                  MenuActivate(G, mx, my, x, y, false, "all_action", namesele);
                  break;
                case cExecSelection:
                  MenuActivate(G, mx, my, x, y, false, "sele_action", namesele);
                  break;
                case cExecObject:
                  switch (rec->obj->type) {
                  case cObjectGroup:
                    MenuActivate(G, mx, my, x, y, false, "group_action", namesele);
                    break;
                  case cObjectMolecule:
                    MenuActivate(G, mx, my, x, y, false, "mol_action", namesele);
                    break;
                  case cObjectMap:
                    MenuActivate(G, mx, my, x, y, false, "map_action", namesele);
                    break;
                  case cObjectSurface:
                    MenuActivate(G, mx, my, x, y, false, "surface_action",
                                 namesele);
                    break;
                  case cObjectMesh:
                    MenuActivate(G, mx, my, x, y, false, "mesh_action", namesele);
                    break;
                  case cObjectMeasurement:
                  case cObjectCGO:
                  case cObjectCallback:
                  case cObjectAlignment:
		  case cObjectVolume:
                    MenuActivate(G, mx, my, x, y, false, "simple_action", namesele);
                    break;
                  case cObjectSlice:
                    MenuActivate(G, mx, my, x, y, false, "slice_action", namesele);
                    break;
                  case cObjectGadget:
                    MenuActivate(G, mx, my, x, y, false, "ramp_action", namesele);
                    break;
                  }
                  break;
                }
                break;
              case 1:
                switch (rec->type) {
                case cExecAll:
                  MenuActivate(G, mx, my, x, y, false, "mol_show", cKeywordAll);
                  break;
                case cExecSelection:
                  MenuActivate(G, mx, my, x, y, false, "mol_show", namesele);
                  break;
                case cExecObject:
                  switch (rec->obj->type) {
                  case cObjectGroup:
                  case cObjectMolecule:
                    MenuActivate(G, mx, my, x, y, false, "mol_show", namesele);
                    break;
                  case cObjectCGO:
                  case cObjectAlignment:
                    MenuActivate(G, mx, my, x, y, false, "cgo_show", namesele);
                    break;
                  case cObjectMeasurement:
                    MenuActivate(G, mx, my, x, y, false, "measurement_show",
                                 namesele);
                    break;
                  case cObjectMap:
                    MenuActivate(G, mx, my, x, y, false, "map_show", namesele);
                    break;
                  case cObjectMesh:
                    MenuActivate(G, mx, my, x, y, false, "mesh_show", namesele);
                    break;
                  case cObjectSurface:
                    MenuActivate(G, mx, my, x, y, false, "surface_show", namesele);
                    break;
                  case cObjectSlice:
                    MenuActivate(G, mx, my, x, y, false, "slice_show", namesele);
                    break;
		  case cObjectVolume:
                    MenuActivate(G, mx, my, x, y, false, "volume_show", namesele);
                    break;

                  }
                  break;
                }
                break;
              case 2:
                switch (rec->type) {
                case cExecAll:
                  MenuActivate(G, mx, my, x, y, false, "mol_hide", cKeywordAll);
                  break;
                case cExecSelection:
                  MenuActivate(G, mx, my, x, y, false, "mol_hide", namesele);
                  break;
                case cExecObject:
                  switch (rec->obj->type) {
                  case cObjectGroup:
                  case cObjectMolecule:
                    MenuActivate(G, mx, my, x, y, false, "mol_hide", namesele);
                    break;
                  case cObjectCGO:
                  case cObjectAlignment:
                    MenuActivate(G, mx, my, x, y, false, "cgo_hide", namesele);
                    break;
                  case cObjectMeasurement:
                    MenuActivate(G, mx, my, x, y, false, "measurement_hide",
                                 namesele);
                    break;
                  case cObjectMap:
                    MenuActivate(G, mx, my, x, y, false, "map_hide", namesele);
                    break;
                  case cObjectMesh:
                    MenuActivate(G, mx, my, x, y, false, "mesh_hide", namesele);
                    break;
                  case cObjectSurface:
                    MenuActivate(G, mx, my, x, y, false, "surface_hide", namesele);
                    break;
                  case cObjectSlice:
                    MenuActivate(G, mx, my, x, y, false, "slice_hide", namesele);
                    break;
		  case cObjectVolume:
                    MenuActivate(G, mx, my, x, y, false, "volume_hide", namesele);
                    break;
		    

                  }
                  break;
                }
                break;
              case 3:
                switch (rec->type) {
                case cExecAll:
                  MenuActivate(G, mx, my, x, y, false, "mol_labels", "(all)");
                  break;
                case cExecSelection:
                  MenuActivate(G, mx, my, x, y, false, "mol_labels", namesele);
                  break;
                case cExecObject:
                  switch (rec->obj->type) {
                  case cObjectGroup:
                  case cObjectMolecule:
                    MenuActivate(G, mx, my, x, y, false, "mol_labels", namesele);
                    break;
                  case cObjectMeasurement:
                    break;
                  case cObjectMap:
                  case cObjectSurface:
                  case cObjectMesh:
                  case cObjectSlice:
                    break;
                  }
                  break;
                }
                break;
              case 4:
                switch (rec->type) {
                case cExecAll:
                case cExecSelection:
                  MenuActivate(G, mx, my, x, y, false, "mol_color", namesele);
                  break;
                case cExecObject:
                  switch (rec->obj->type) {
                  case cObjectGroup:
                  case cObjectMolecule:
                    MenuActivate(G, mx, my, x, y, false, "mol_color", namesele);
                    break;
                  case cObjectMap:
                  case cObjectCGO:
                  case cObjectAlignment:
                    MenuActivate(G, mx, my, x, y, false, "general_color", namesele);
                    break;
                  case cObjectMesh:
                    MenuActivate(G, mx, my, x, y, false, "mesh_color", namesele);
                    break;
                  case cObjectSurface:
                    MenuActivate2Arg(G, mx, my, x, y, false, "mesh_color", namesele, "surface");
                    break;
                  case cObjectMeasurement:
                    MenuActivate(G, mx, my, x, y, false, "measurement_color", namesele);
                    break;
                  case cObjectSlice:
                    MenuActivate(G, mx, my, x, y, false, "slice_color", namesele);
                    break;
                  case cObjectVolume:
                    MenuActivate(G, mx, my, x, y, false, "vol_color", namesele);
                    break;
                  case cObjectGadget:
                    MenuActivate(G, mx, my, x, y, false, "ramp_color", namesele);
                    break;
                  }
                  break;
                }
                break;
              case 5:
                switch (rec->type) {
                case cExecAll:
                  MenuActivate0Arg(G, mx, my, x, y, false, "camera_motion");
                  break;
                case cExecSelection: /* not relevant at present */
		  break;
                case cExecObject: /* this kinds of objects can be animated */
                  switch (rec->obj->type) {
                  case cObjectGroup: /* however, this just groups motions for a set of objects */
                  case cObjectMolecule:
		  case cObjectMeasurement:
		  case cObjectMap:
		  case cObjectSurface:
		  case cObjectCGO:
		  case cObjectMesh:
                    MenuActivate(G, mx, my, x, y, false, "obj_motion", namesele);
                    break;
                  }
                  break;
                }
                break;
              }
            } else {            /* clicked in variable area */

              if(((panel->is_group) && (((xx) - 1) / DIP2PIXEL(8)) > (panel->nest_level + 1)) ||
                 ((!panel->is_group) && (((xx) - 1) / DIP2PIXEL(8)) > panel->nest_level)) {
                /* clicked on name */

                rec->hilight = 1;
                switch (button) {
                case P_GLUT_LEFT_BUTTON:
                  if (I->DragMode != ExecutiveDragMode::Off) {
                    break;
                  }
                  I->Pressed = n;
                  I->OldVisibility = rec->visible;
                  I->Over = n;
                  I->DragMode = ExecutiveDragMode::Visibility;
                  I->ToggleMode = ExecutiveToggleMode::DeferVisibility;
                  I->LastChanged = NULL;
                  I->LastZoomed = NULL;
                  if(mod == (cOrthoSHIFT | cOrthoCTRL)) {
                    I->ToggleMode = ExecutiveToggleMode::HoverActivate;
                    if(!rec->visible) {
                      ExecutiveSpecSetVisibility(G, rec, !I->OldVisibility, mod, false);
                    }
                    if(rec != I->LastZoomed)
                      ExecutiveWindowZoom(G, rec->name, 0.0F, -1, false, -1.0F, true);
                    I->LastZoomed = rec;
                    I->LastChanged = rec;
                  } else if(mod & cOrthoSHIFT) {
                    ExecutiveSpecSetVisibility(G, rec, !I->OldVisibility, mod, false);
                    I->ToggleMode = ExecutiveToggleMode::ImmediateVisibility;
                  } else if(mod & cOrthoCTRL) {
                    I->ToggleMode = ExecutiveToggleMode::HoverActivate;
                    if(!rec->visible) {
                      ExecutiveSpecSetVisibility(G, rec, !I->OldVisibility, mod, false);
                    }
                    I->LastChanged = rec;
                  }
                  I->PressedWhat = 1;
                  I->OverWhat = 1;
                  break;
                case P_GLUT_MIDDLE_BUTTON:
                  if (I->DragMode != ExecutiveDragMode::Off) {
                    break;
                  }
                  I->Pressed = n;
                  I->OldVisibility = rec->visible;
                  I->Over = n;
                  I->LastOver = I->Over;
                  I->DragMode = ExecutiveDragMode::VisibilityWithCamera;
                  I->ToggleMode = ExecutiveToggleMode::DeferVisibility;
                  I->LastChanged = rec;
                  I->LastZoomed = NULL;
                  if(mod & cOrthoCTRL) {
                    I->ToggleMode = ExecutiveToggleMode::ZoomActivateDeactivatePrevious;
                    ExecutiveWindowZoom(G, rec->name, 0.0F, -1, false, -1.0F, true);
                    I->LastZoomed = rec;
                    if(mod & cOrthoSHIFT) {     /* exclusive */
                      I->ToggleMode = ExecutiveToggleMode::ZoomExclusiveActivate;
                      ExecutiveSetObjVisib(G, cKeywordAll, false, false);       
                      /* need to log this */
                      if(!rec->visible)
                        ExecutiveSpecSetVisibility(G, rec, true, 0, true);
                    }
                  } else {
                    I->ToggleMode = ExecutiveToggleMode::CenterActivateDeactivatePrevious;
                    ExecutiveCenter(G, rec->name, -1, true, -1.0F, NULL, true);
                  }
                  if(!rec->visible) {
                    ExecutiveSpecSetVisibility(G, rec, !rec->visible, mod, false);
                    I->LastChanged = rec;
                  }
                  I->PressedWhat = 1;
                  I->OverWhat = 1;
                  break;
                case P_GLUT_RIGHT_BUTTON:
                  if (I->DragMode != ExecutiveDragMode::Off) {
                    break;
                  }
                  I->DragMode = ExecutiveDragMode::Reorder;      /* reorder */
                  I->Pressed = n;
                  I->Over = n;
                  I->LastOver = I->Over;
                  I->PressedWhat = 1;
                  I->OverWhat = 1;
                  break;
                }
                OrthoGrab(G, this);
                OrthoDirty(G);
              } else if(panel->is_group) {
                /* clicked on group control */
                rec->hilight = 2;
                I->DragMode = ExecutiveDragMode::Visibility;
                I->Pressed = n;
                I->Over = n;
                I->LastOver = I->Over;
                I->PressedWhat = 2;
                I->OverWhat = 2;

                OrthoGrab(G, this);
                OrthoDirty(G);
              }
            }
          }
          a--;
        }
      }
    }
  }
  PyMOL_NeedRedisplay(G->PyMOL);
  return (1);
}


/*========================================================================*/
static void ExecutiveSpecEnable(PyMOLGlobals * G, SpecRec * rec, int parents, int log)
{
  if(log && SettingGetGlobal_b(G, cSetting_logging)) {
    OrthoLineType buffer = "";
    sprintf(buffer, "cmd.enable('%s',%d)", rec->obj->Name, parents);
    PLog(G, buffer, cPLog_pym);
  }

  if(!rec->visible) {
    rec->visible = true;
    ReportEnabledChange(G, rec);
  }
  if(!rec->in_scene) {
    rec->in_scene = SceneObjectAdd(G, rec->obj);
  }

  if(parents) {
    CExecutive *I = G->Executive;
    CTracker *I_Tracker = I->Tracker;
    int list_id = ExecutiveGetObjectParentList(G, rec);
    if(list_id) {
      int iter_id = TrackerNewIter(I_Tracker, 0, list_id);
      SpecRec *parent_rec = NULL;

      while(TrackerIterNextCandInList(I_Tracker, iter_id,
                                      (TrackerRef **) (void *) &parent_rec)) {
        if(rec) {
          switch (parent_rec->type) {
          case cExecObject:
            if(!parent_rec->in_scene) {
              parent_rec->in_scene = SceneObjectAdd(G, parent_rec->obj);
            }
            if(!parent_rec->visible) {
              parent_rec->visible = true;
	      ReportEnabledChange(G, parent_rec);
            }
          }
        }
      }
      TrackerDelIter(I_Tracker, iter_id);
    }
    TrackerDelList(I_Tracker, list_id);
  }
  ExecutiveInvalidateSceneMembers(G);
}

void ExecutiveSpecSetVisibility(PyMOLGlobals * G, SpecRec * rec,
                                int new_vis, int mod, int parents)
{
  std::string buffer;
  int logging = SettingGetGlobal_i(G, cSetting_logging);
  if(rec->type == cExecObject) {
    if(rec->visible && !new_vis) {
      if(logging)
        buffer = pymol::string_format("cmd.disable('%s')", rec->obj->Name);
      SceneObjectDel(G, rec->obj, true);
      ExecutiveInvalidateSceneMembers(G);
      if (rec->visible != new_vis){
	rec->visible = new_vis;
	ReportEnabledChange(G, rec);
      }
    } else if((!rec->visible) && new_vis) {
      ExecutiveSpecEnable(G, rec, parents, logging);
      /*
         if(logging) 
         sprintf(buffer,"cmd.enable('%s')",rec->obj->Name);
         rec->in_scene = SceneObjectAdd(G,rec->obj);
         rec->visible=new_vis;
         ExecutiveInvalidateSceneMembers(G);
       */
    }
    SceneChanged(G);
    if(logging && buffer[0]) {
      PLog(G, buffer, cPLog_pym);
    }
  } else if(rec->type == cExecAll) {
    if(SettingGetGlobal_i(G, cSetting_logging)) {
      if(rec->visible)
        buffer = "cmd.disable('all')";
      else
        buffer = "cmd.enable('all')";
      PLog(G, buffer, cPLog_pym);
    }
    ExecutiveSetObjVisib(G, cKeywordAll, !rec->visible, false);
  } else if(rec->type == cExecSelection) {
    if(mod & cOrthoCTRL) {
      buffer = pymol::string_format("cmd.enable('%s')", rec->name);
      PLog(G, buffer, cPLog_pym);
      if (!rec->visible){
	rec->visible = true;
	ReportEnabledChange(G, rec);
      }
    } else {

      if(rec->visible && !new_vis) {
        if(SettingGetGlobal_i(G, cSetting_logging))
          buffer = pymol::string_format("cmd.disable('%s')", rec->name);
      } else if((!rec->visible) && new_vis) {
        buffer = pymol::string_format("cmd.enable('%s')", rec->name);
      }
      if(new_vis && SettingGetGlobal_b(G, cSetting_active_selections)) {
        ExecutiveHideSelections(G);
      }
      if(SettingGetGlobal_i(G, cSetting_logging)) {
        PLog(G, buffer, cPLog_pym);
      }
      if (rec->visible != new_vis){
	rec->visible = new_vis;
	ReportEnabledChange(G, rec);
      }
    }
    SceneChanged(G);
  }
}

int ExecutiveAssignAtomTypes(PyMOLGlobals * G, const char *s1, int format, int state, int quiet)
{
  int result = 0;
  int sele1;

  sele1 = SelectorIndexByName(G, s1);
  if(state < 0)
    state = 0;
  {
    if(sele1 >= 0) {
      result = SelectorAssignAtomTypes(G, sele1, state, quiet, format);
    }
  }
  return (result);
}

int CExecutive::release(int button, int x, int y, int mod)
{
  PyMOLGlobals *G = m_G;
  CExecutive *I = G->Executive;
  SpecRec *rec = NULL;
  int pass = false;
  int skip;
  int xx;
#ifndef NDEBUG
  int hide_underscore = SettingGetGlobal_b(G, cSetting_hide_underscore_names);
#endif
  if(y < I->HowFarDown) {
    if(SettingGet<InternalGUIMode>(cSetting_internal_gui_mode, G->Setting) != InternalGUIMode::Default)
      return SceneGetBlock(G)->release(button, x, y, mod);
  }

  xx = (x - rect.left);
  if(I->ScrollBarActive) {
    if((x - rect.left) <
       (ExecScrollBarWidth + ExecScrollBarMargin + ExecToggleMargin)) {
      pass = 1;
      I->m_ScrollBar.release(button, x, y, mod);
      OrthoUngrab(G);
    }
    xx -= (ExecScrollBarWidth + ExecScrollBarMargin);
  }

  skip = I->NSkip;

  if(!pass) {
    I->drag(x, y, mod);    /* incorporate final changes in cursor position */
    switch (I->DragMode) {
    case ExecutiveDragMode::Visibility:

      /*while(ListIterate(I->Spec,rec,next)) { */
      for (auto& panelitem : I->Panel) {
        auto* const panel = &panelitem;
        rec = panel->spec;

        assert(rec->name[0] != '_' || !hide_underscore);
        {
          if(skip) {
            skip--;
          } else {
            if((I->PressedWhat == 1) &&
               (((panel->is_group) && ((xx - 1) / DIP2PIXEL(8)) > (panel->nest_level + 1)) ||
                ((!panel->is_group) && ((xx - 1) / DIP2PIXEL(8)) > panel->nest_level))) {
              /* over name */
              if(rec->hilight == 1) {
                if(rec->type == cExecSelection) {
                  ExecutiveSpecSetVisibility(G, rec, !I->OldVisibility, 0, false);
                } else {
                  ExecutiveSpecSetVisibility(G, rec, !I->OldVisibility, mod, true);
                }
              }
            } else if((I->PressedWhat == 2) && (panel->is_group)) {
              if(rec->hilight == 2) {
                ObjectGroup *obj = (ObjectGroup *) rec->obj;
                OrthoLineType buf2;
                sprintf(buf2, "cmd.group(\"%s\",action='%s')\n", rec->obj->Name,
                        obj->OpenOrClosed ? "close" : "open");
                PLog(G, buf2, cPLog_no_flush);
                ExecutiveGroup(G, rec->obj->Name, "", 5, 1);

              }
              /* over group control */
            }
          }
        }
      }
      break;
    case ExecutiveDragMode::Reorder:
      if(I->ReorderFlag) {
        I->ReorderFlag = false;
        PLog(G, I->ReorderLog, cPLog_no_flush);
      }
      break;
    }
  }

  {
    SpecRec *rec = NULL;
    while(ListIterate(I->Spec, rec, next)) {
      rec->hilight = 0;
    }
  }

  I->Over = -1;
  I->Pressed = -1;
  I->DragMode = ExecutiveDragMode::Off;
  I->PressedWhat = 0;
  OrthoUngrab(G);
  PyMOL_NeedRedisplay(G->PyMOL);
  return (1);
}


/*========================================================================*/
int CExecutive::drag(int x, int y, int mod)
{
  PyMOLGlobals *G = m_G;
  CExecutive *I = G->Executive;
  int xx, t;
  int ExecLineHeight = DIP2PIXEL(SettingGetGlobal_i(G, cSetting_internal_gui_control_size));
#ifndef NDEBUG
  int hide_underscore = SettingGetGlobal_b(G, cSetting_hide_underscore_names);
#endif
  int op_cnt = get_op_cnt(G);

  ExecutiveUpdateGroups(G, false);
  ExecutiveUpdatePanelList(G);
  if(y < I->HowFarDown) {
    if(SettingGet<InternalGUIMode>(cSetting_internal_gui_mode, G->Setting) != InternalGUIMode::Default)
      return SceneGetBlock(G)->drag(x, y, mod);
  }

  if(I->DragMode != ExecutiveDragMode::Off) {
    xx = (x - rect.left);
    t = ((rect.right - ExecRightMargin) - x) / ExecToggleWidth;
    if(I->ScrollBarActive) {
      xx -= (ExecScrollBarWidth + ExecScrollBarMargin);
    }

    {
      int row_offset;
      if((xx >= 0) && (t >= op_cnt)) {
        row_offset = ((rect.top - y) -
                      (ExecTopMargin + ExecClickMargin)) / ExecLineHeight;
        I->Over = row_offset;
      } else {
        I->Over = -1;
        row_offset = -1;
        {
          SpecRec *rec = NULL;
          while(ListIterate(I->Spec, rec, next))
            rec->hilight = 0;
        }
      }

      if(I->RecoverPressed) {
        SpecRec *rec = NULL;
        int skip = I->NSkip;
        int row = 0;
        for (auto& panelitem : I->Panel) {
          auto* const panel = &panelitem;
          rec = panel->spec;

          assert(rec->name[0] != '_' || !hide_underscore);
          {
            if(skip) {
              skip--;
            } else {
              if(rec == I->RecoverPressed) {
                I->Pressed = row;
                I->RecoverPressed = NULL;
              }
              row++;
            }
          }
        }
      }

      if(I->PressedWhat == 2) {
        I->OverWhat = 0;
      }
      if(I->Over >= 0) {
        SpecRec *rec = NULL;
        int skip = I->NSkip;
        int row = 0;
        switch (I->DragMode) {
          case ExecutiveDragMode::Visibility:

          for (auto& panelitem : I->Panel) {
            auto* const panel = &panelitem;
            rec = panel->spec;

            assert(rec->name[0] != '_' || !hide_underscore);
            {
              if(skip) {
                skip--;
              } else {

                rec->hilight = 0;
                if((I->PressedWhat == 1) &&     /* name button */
                   (((row >= I->Over) && (row <= I->Pressed)) ||
                    ((row >= I->Pressed) && (row <= I->Over)))) {

                  if(((panel->is_group) && ((xx - 1) / DIP2PIXEL(8)) > (panel->nest_level + 1)) ||
                     ((!panel->is_group) && ((xx - 1) / DIP2PIXEL(8)) > panel->nest_level)) {
                    /* dragged over name */

                    I->OverWhat = 1;

                    switch (I->ToggleMode) {
                    case ExecutiveToggleMode::DeferVisibility:
                      if(row || (row == I->Pressed))
                        rec->hilight = 1;
                      break;
                    case ExecutiveToggleMode::ImmediateVisibility:
                      if(row)
                        ExecutiveSpecSetVisibility(G, rec, !I->OldVisibility, mod, false);
                      break;
                    case ExecutiveToggleMode::HoverActivate:
                      if((row == I->Over) && row) {
                        if(I->LastChanged != rec) {
                          ExecutiveSpecSetVisibility(G, I->LastChanged, false, mod,
                                                     false);
                        }
                        if(!rec->visible) {
                          ExecutiveSpecSetVisibility(G, rec, true, mod, false);
                          I->LastChanged = rec;
                        }
                        if(mod == (cOrthoSHIFT | cOrthoCTRL)) {
                          if(rec != I->LastZoomed)
                            ExecutiveWindowZoom(G, rec->name, 0.0F, -1, false, -1.0F,
                                                true);
                          I->LastZoomed = rec;
                        }
                      }
                      break;
                    }
                  }
                } else if((row == I->Pressed) && (I->PressedWhat == 2)) {
                  if(!((panel->is_group) && ((xx - 1) / DIP2PIXEL(8)) > (panel->nest_level + 1))) {

                    /* on group control */

                    I->OverWhat = 2;
                    if(I->PressedWhat == I->OverWhat) {
                      rec->hilight = 2;
                    }
                  } else {
                    I->OverWhat = 0;
                  }
                }
                row++;
              }
            }
          }
          break;
        case ExecutiveDragMode::Reorder:                /* right buttonBROKEN */
          {
            if((I->Over != I->Pressed) && (I->LastOver != I->Over)) {
              SpecRec *new_rec = NULL;
              SpecRec *mov_rec = NULL;

              for (auto& panelitem : I->Panel) {
                rec = panelitem.spec;
                assert(rec->name[0] != '_' || !hide_underscore);
                {
                  {
                    if(skip) {
                      skip--;
                    } else {
                      if(row == I->Pressed) {
                        mov_rec = rec;
                      }
                      if(row == I->Over) {
                        new_rec = rec;
                      }
                      row++;
                    }
                  }
                }
              }
              {
                int group_flag = false;
                int order_flag = false;
                int location = 0;
                char *first = NULL, *second = NULL;
                int is_child = false;

                if(mov_rec && (!new_rec) && (I->Over > I->Pressed) && mov_rec->group) {
                  first = mov_rec->group->name;
                  second = mov_rec->name;
                  order_flag = true;
                  strcpy(mov_rec->group_name, mov_rec->group->group_name);
                  group_flag = true;
                } else if(mov_rec && new_rec) {
                  if (new_rec->isChildOf(mov_rec)) {
                    /* do nothing when a group is dragged over one of its members */
                  } else {

                    if(I->Pressed < I->Over) {  /* downward */
                      first = new_rec->name;
                      second = mov_rec->name;
                      order_flag = true;
                    } else {    /* upward */
                      first = mov_rec->name;
                      second = new_rec->name;
                      order_flag = true;
                      location = -2;    /* upper */
                    }

                    if(mov_rec->group == new_rec->group) {      /* reordering within a group level */
                      if((new_rec->type == cExecObject)
                         && (new_rec->obj->type == cObjectGroup)) {
                        ObjectGroup *group = (ObjectGroup *) new_rec->obj;
                        if(group->OpenOrClosed && !is_child) {
                          /* put inside an open group */
                          strcpy(mov_rec->group_name, new_rec->name);
                          order_flag = false;
                          group_flag = true;
                        }
                      }
                    } else if((mov_rec->group != new_rec) && (new_rec->group != mov_rec)) {

                      if((new_rec->type == cExecObject)
                         && (new_rec->obj->type == cObjectGroup)) {
                        ObjectGroup *group = (ObjectGroup *) new_rec->obj;
                        if(group->OpenOrClosed && !is_child) {
                          /* put inside group */
                          strcpy(mov_rec->group_name, new_rec->name);
                        } else {
                          strcpy(mov_rec->group_name, new_rec->group_name);
                        }

                        /* WLD TEST */

                        if(I->Pressed < I->Over) {      /* downward */
                          first = new_rec->name;
                          second = mov_rec->name;
                          order_flag = true;
                        } else {        /* upward */
                          first = mov_rec->name;
                          second = new_rec->name;
                          order_flag = true;
                        }

                        /* WLD END */

                      } else
                        strcpy(mov_rec->group_name, new_rec->group_name);
                      group_flag = true;
                    }
                  }
                }

                if(group_flag) {
                  OrthoLineType buf;
                  if(mov_rec->group_name[0]) {
                    sprintf(buf, "group %s, %s\n", mov_rec->group_name, mov_rec->name);
                  } else {
                    sprintf(buf, "ungroup %s\n", mov_rec->name);
                  }
                  PLog(G, buf, cPLog_no_flush);
                  ExecutiveInvalidateGroups(G, false);
                  I->RecoverPressed = mov_rec;
                  I->Pressed = 0;
                }
                if(order_flag && first && second) {
                  OrthoLineType order_input;
                  sprintf(order_input, "%s %s", first, second);
                  ExecutiveOrder(G, order_input, false, location);
                  sprintf(I->ReorderLog, "cmd.order(\"%s\",location=\"%s\")\n",
                          order_input, (location == -2) ? "upper" : "current");
                  PLog(G, I->ReorderLog, cPLog_no_flush);
                  I->RecoverPressed = mov_rec;
                  I->Pressed = 0;
                }
              }
            }
          }
          break;
        case ExecutiveDragMode::VisibilityWithCamera:                /* middle button */
          for (auto& panelitem : I->Panel) {
            rec = panelitem.spec;
            assert(rec->name[0] != '_' || !hide_underscore);
            {
              if(skip) {
                skip--;
              } else {
                rec->hilight = 0;
                if(((row >= I->Over) && (row <= I->Pressed)) ||
                   ((row >= I->Pressed) && (row <= I->Over))) {
                  switch (I->ToggleMode) {
                   case ExecutiveToggleMode::CenterActivateDeactivatePrevious:      /* center and activate, while deactivating previous */
                    if((row == I->Over) && row) {
                      if(I->LastChanged != rec) {
                        ExecutiveSpecSetVisibility(G, I->LastChanged, false, mod, false);
                        ExecutiveCenter(G, rec->name, -1, true, -1.0F, NULL, true);
                        if(!rec->visible)
                          ExecutiveSpecSetVisibility(G, rec, true, mod, false);
                        I->LastChanged = rec;
                      }
                      rec->hilight = 0;
                    }
                    break;
                   case ExecutiveToggleMode::ZoomActivateDeactivatePrevious:      /* zoom and activate, while deactivating previous */
                    if((row == I->Over) && row) {
                      if(I->LastChanged != rec) {
                        ExecutiveSpecSetVisibility(G, I->LastChanged, false, mod, false);
                        ExecutiveWindowZoom(G, rec->name, 0.0F, -1, false, -1.0F, true);
                        if(!rec->visible)
                          ExecutiveSpecSetVisibility(G, rec, true, mod, false);
                        I->LastChanged = rec;
                      }
                      rec->hilight = 1;
                    }
                    break;
                   case ExecutiveToggleMode::ZoomExclusiveActivate:      /* zoom and make only object enabled */
                    if((row == I->Over) && row) {
                      if(rec != I->LastZoomed) {
                        ExecutiveSpecSetVisibility(G, I->LastZoomed, false, mod, false);
                        ExecutiveWindowZoom(G, rec->name, 0.0F, -1, false, -1.0F, true);
                        I->LastZoomed = rec;
                        ExecutiveSpecSetVisibility(G, rec, true, 0, false);
                      }
                      rec->hilight = 1;
                    }
                  }
                }
                row++;
              }
            }
          }
          break;
        }
        I->LastOver = I->Over;

      } else if(I->LastChanged)
        ExecutiveSpecSetVisibility(G, I->LastChanged, false, mod, false);
      OrthoDirty(G);
    }
  }
  return (1);
}

static void draw_button(int x2, int y2, int w, int h, const float *light, const float *dark,
                        const float *inside ORTHOCGOARG)
{
  if (orthoCGO){
    CGOColorv(orthoCGO, light);
    CGOBegin(orthoCGO, GL_TRIANGLE_STRIP);
    CGOVertex(orthoCGO, x2, y2, 0.f);
    CGOVertex(orthoCGO, x2, y2 + h, 0.f);
    CGOVertex(orthoCGO, x2 + w, y2, 0.f);
    CGOVertex(orthoCGO, x2 + w, y2 + h, 0.f);
    CGOEnd(orthoCGO);
  } else {
    glColor3fv(light);
    glBegin(GL_POLYGON);
    glVertex2i(x2, y2);
    glVertex2i(x2, y2 + h);
    glVertex2i(x2 + w, y2 + h);
    glVertex2i(x2 + w, y2);
    glEnd();
  }

  if (orthoCGO){
    CGOColorv(orthoCGO, dark);
    CGOBegin(orthoCGO, GL_TRIANGLE_STRIP);
    CGOVertex(orthoCGO, x2 + 1, y2, 0.f);
    CGOVertex(orthoCGO, x2 + 1, y2 + h - 1, 0.f);
    CGOVertex(orthoCGO, x2 + w, y2, 0.f);
    CGOVertex(orthoCGO, x2 + w, y2 + h - 1, 0.f);
    CGOEnd(orthoCGO);
  } else {
    glColor3fv(dark);
    glBegin(GL_POLYGON);
    glVertex2i(x2 + 1, y2);
    glVertex2i(x2 + 1, y2 + h - 1);
    glVertex2i(x2 + w, y2 + h - 1);
    glVertex2i(x2 + w, y2);
    glEnd();
  }

  if(inside) {
    if (orthoCGO){
      CGOColorv(orthoCGO, inside);
      CGOBegin(orthoCGO, GL_TRIANGLE_STRIP);
      CGOVertex(orthoCGO, x2 + 1, y2 + 1, 0.f);
      CGOVertex(orthoCGO, x2 + 1, y2 + h - 1, 0.f);
      CGOVertex(orthoCGO, x2 + w - 1, y2 + 1, 0.f);
      CGOVertex(orthoCGO, x2 + w - 1, y2 + h - 1, 0.f);
      CGOEnd(orthoCGO);
    } else {
      glColor3fv(inside);
      glBegin(GL_POLYGON);
      glVertex2i(x2 + 1, y2 + 1);
      glVertex2i(x2 + 1, y2 + h - 1);
      glVertex2i(x2 + w - 1, y2 + h - 1);
      glVertex2i(x2 + w - 1, y2 + 1);
      glEnd();
    }
  } else {                      /* rainbow */
    if (orthoCGO){
      CGOBegin(orthoCGO, GL_TRIANGLE_STRIP);
      CGOColor(orthoCGO, 0.1F, 1.0F, 0.1F); // green
      CGOVertex(orthoCGO, x2 + 1, y2 + h - 1, 0.f);
      CGOColor(orthoCGO, 1.0F, 1.0F, 0.1F);  // yellow
      CGOVertex(orthoCGO, x2 + w - 1, y2 + h - 1, 0.f);
      CGOColor(orthoCGO, 1.f, 0.1f, 0.1f); // red
      CGOVertex(orthoCGO, x2 + 1, y2 + 1, 0.f);
      CGOColor(orthoCGO, 0.1F, 0.1F, 1.0F);  // blue
      CGOVertex(orthoCGO, x2 + w - 1, y2 + 1, 0.f);
      CGOEnd(orthoCGO);
    } else {
      glBegin(GL_POLYGON);
      glColor3f(1.0F, 0.1F, 0.1F);
      glVertex2i(x2 + 1, y2 + 1);
      glColor3f(0.1F, 1.0F, 0.1F);
      glVertex2i(x2 + 1, y2 + h - 1);
      glColor3f(1.0F, 1.0F, 0.1F);
      glVertex2i(x2 + w - 1, y2 + h - 1);
      glColor3f(0.1F, 0.1F, 1.0F);
      glVertex2i(x2 + w - 1, y2 + 1);
      glEnd();
    }
  }
}

#ifndef _PYMOL_NOPY
static void draw_button_char(PyMOLGlobals * G, int x2, int y2, char ch ORTHOCGOARG)
{
  TextSetColor3f(G, 0.0F, 0.0F, 0.0F);
  TextSetPos2i(G, x2 + ExecToggleTextShift, y2);
  TextDrawChar(G, ch ORTHOCGOARGVAR);
}
#endif

/**
 * Get the text color for the object name.
 * @param obj Object pointer, may be NULL
 * @param mode 0: default, 1: first carbon atom color, 2: object color
 * @param bg_rgb Background color
 * @param default_rgb Default text color
 * @return pointer to borrowed RGB array
 */
static const float* getNameColor(
    const pymol::CObject* obj, int mode, const float* bg_rgb, const float* default_rgb)
{
  enum {
    cNameColorMode_default = 0,
    cNameColorMode_carbon = 1,
    cNameColorMode_object = 2,
  };

  if (mode == cNameColorMode_default || !obj) {
    return default_rgb;
  }

  int color = cColorDefault;

  switch (mode) {
  case cNameColorMode_carbon:
    // First carbon atom color
    // Assuming that there is typically a carbon atom within the first few
    // atoms, this procedure should be cheap (O(1)).
    if (obj->type == cObjectMolecule) {
      auto objmol = static_cast<const ObjectMolecule*>(obj);
      for (auto ai = objmol->AtomInfo.data(), ai_end = ai + objmol->NAtom;
           ai != ai_end; ++ai) {
        if (ai->protons == cAN_C) {
          color = ai->color;
          break;
        }
      }
    }
    break;
  case cNameColorMode_object:
    // Object color
    color = obj->Color;
    break;
  }

  if (color != cColorDefault) {
    const float* rgb = ColorGet(obj->G, color);

    // only use if it's different from the background color
    if (!within3f(bg_rgb, rgb, 0.1f)) {
      return rgb;
    }
  }

  return default_rgb;
}

/*========================================================================*/
void CExecutive::draw(CGO* orthoCGO)
{
  PyMOLGlobals *G = m_G;
  int x, y, xx, x2, y2;
  WordType ch;
  char *c = NULL;
  float enabledColor[3] = { 0.5F, 0.5F, 0.5F };
  float cloakedColor[3] = { 0.35F, 0.35F, 0.35F };
  float pressedColor[3] = { 0.7F, 0.7F, 0.7F };
  float disabledColor[3] = { 0.25F, 0.25F, 0.25F };
  float lightEdge[3] = { 0.6F, 0.6F, 0.6F };
  float darkEdge[3] = { 0.35F, 0.35F, 0.35F };
  float captionColor[3] = { 0.3F, 0.9F, 0.3F };

#ifndef _PYMOL_NOPY
  float activeColor[3] = { 0.9F, 0.9F, 1.0F };
  float toggleColor3[3] = { 0.6F, 0.6F, 0.8F };
#endif

  SpecRec *rec = NULL;
  CExecutive *I = G->Executive;
  int n_ent;
  int n_disp;
  int skip = 0;
  int row = -1;
  int ExecLineHeight = DIP2PIXEL(SettingGetGlobal_i(G, cSetting_internal_gui_control_size));
  int text_lift = (ExecLineHeight / 2) - DIP2PIXEL(5);
  int op_cnt = get_op_cnt(G);
  int full_names = SettingGetGlobal_b(G, cSetting_group_full_member_names);
  int arrows = SettingGetGlobal_b(G, cSetting_group_arrow_prefix);
  auto name_color_mode = SettingGet<int>(G, cSetting_internal_gui_name_color_mode);

  ExecutiveUpdatePanelList(G);
  
  /* if we're running with a GUI and have a valid panel */
  if(G->HaveGUI && G->ValidContext && ((rect.right - rect.left) > 6)
     && !I->Panel.empty()) {
    int max_char;
    int nChar;

    /* count entries
     * do we have enough structures to warrant a scroll bar? */
    n_ent = 0;
    for (auto& panelitem : I->Panel) {
      rec = panelitem.spec;
      assert(rec && (rec->name[0] != '_' ||
                        !SettingGet<bool>(G, cSetting_hide_underscore_names)));
        n_ent++;
    }

    n_disp =
      ((rect.top - rect.bottom) - (ExecTopMargin)) / ExecLineHeight;
    if(n_disp < 1)
      n_disp = 1;

    /* we need a scrollbar */
    if(n_ent > n_disp) {
      int bar_maxed = I->m_ScrollBar.isMaxed();
      if(!I->ScrollBarActive) {
        I->m_ScrollBar.setLimits(n_ent, n_disp);
        if(bar_maxed) {
          I->m_ScrollBar.maxOut();
          I->NSkip = static_cast<int>(I->m_ScrollBar.getValue());
        } else {
          I->m_ScrollBar.setValue(0);
          I->NSkip = 0;
        }
      } else {
        I->m_ScrollBar.setLimits(n_ent, n_disp);
        if(bar_maxed)
          I->m_ScrollBar.maxOut();
        I->NSkip = static_cast<int>(I->m_ScrollBar.getValue());
      }
      I->ScrollBarActive = 1;
    } else {
      I->ScrollBarActive = 0;
      I->NSkip = 0;
    }

    /* determination of longest string based on internal_gui_size, etc... */
    max_char =
      (((rect.right - rect.left) -
        (ExecLeftMargin + ExecRightMargin + 4)) - (op_cnt * ExecToggleWidth));
    if(I->ScrollBarActive) {
      max_char -= (ExecScrollBarMargin + ExecScrollBarWidth);
    }
    max_char /= DIP2PIXEL(8);

    /* fill and outline the entire block */
    if(SettingGet<InternalGUIMode>(cSetting_internal_gui_mode, G->Setting) == InternalGUIMode::Default) {
      if (orthoCGO)
	CGOColorv(orthoCGO, BackColor);
#ifndef PURE_OPENGL_ES_2
      else
	glColor3fv(BackColor);
#endif
      fill(orthoCGO);
      drawLeftEdge(orthoCGO);
    }

    /* draw the scroll bar */
    if(I->ScrollBarActive) {
      I->m_ScrollBar.setBox(rect.top - ExecScrollBarMargin,
                      rect.left + ExecScrollBarMargin,
                      rect.bottom + 2,
                      rect.left + ExecScrollBarMargin + ExecScrollBarWidth);
      I->m_ScrollBar.draw(orthoCGO);
    }

    x = rect.left + ExecLeftMargin;
    y = (rect.top - ExecLineHeight) - ExecTopMargin;
    /*    xx = rect.right-ExecRightMargin-ExecToggleWidth*(cRepCnt+op_cnt); */
#ifndef _PYMOL_NOPY
    xx = rect.right - ExecRightMargin - ExecToggleWidth * (op_cnt);
#else
    xx = rect.right - ExecRightMargin;
#endif

    if(I->ScrollBarActive) {
      x += ExecScrollBarWidth + ExecScrollBarMargin;
    }
    skip = I->NSkip;

    /* for each object in the Panel... */
    for (auto& panelitem : I->Panel) {
      auto* const panel = &panelitem;
      rec = panel->spec;
      {
        if(skip) {
          skip--;
        } else {
	  /* setup the X,Y offsets for this entry */
          row++;
          x2 = xx;
          y2 = y;
          nChar = max_char;

          if((x - ExecToggleMargin) - (xx - ExecToggleMargin) > -10) {
            x2 = x + 10;
          }
#ifndef _PYMOL_NOPY
	  /*
	   * The ASHLC menus; these access Python so,
	   * protect this block from non-python instances
	   */
          {
            int a;
            float toggleColor[3] = { 0.5F, 0.5F, 1.0F };
            float toggleColor2[3] = { 0.4F, 0.4F, 0.6F };
            float toggleDarkEdge[3] = { 0.3F, 0.3F, 0.5F };
            float toggleLightEdge[3] = { 0.7F, 0.7F, 0.9F };

            glColor3fv(toggleColor);
            for(a = 0; a < op_cnt; a++) {
              switch (a) {
              case 0:
                /*
                   glColor3fv(toggleColor);
                   glBegin(GL_POLYGON);
                   glVertex2i(x2,y2+(ExecToggleSize)/2);
                   glVertex2i(x2+(ExecToggleSize)/2,y2);
                   glVertex2i(x2+ExecToggleSize,y2+(ExecToggleSize)/2);
                   glVertex2i(x2+(ExecToggleSize)/2,y2+ExecToggleSize);
                   glEnd();
                 */

		/* the infamous ASHLC! */

                draw_button(x2, y2, ExecToggleSize, (ExecLineHeight - 1),
                            toggleLightEdge, toggleDarkEdge, toggleColor ORTHOCGOARGVAR);

                draw_button_char(G, x2, y2 + text_lift, 'A' ORTHOCGOARGVAR);
                break;
              case 1:
                draw_button(x2, y2, ExecToggleSize, (ExecLineHeight - 1),
                            toggleLightEdge, toggleDarkEdge, toggleColor3 ORTHOCGOARGVAR);

                draw_button_char(G, x2, y2 + text_lift, 'S' ORTHOCGOARGVAR);
                break;
              case 2:
                draw_button(x2, y2, ExecToggleSize, (ExecLineHeight - 1),
                            toggleLightEdge, toggleDarkEdge, toggleColor2 ORTHOCGOARGVAR);
                draw_button_char(G, x2, y2 + text_lift, 'H' ORTHOCGOARGVAR);
                break;
              case 3:
                draw_button(x2, y2, ExecToggleSize, (ExecLineHeight - 1),
                            toggleLightEdge, toggleDarkEdge, toggleColor ORTHOCGOARGVAR);
                draw_button_char(G, x2, y2 + text_lift, 'L' ORTHOCGOARGVAR);
                break;
              case 4:
                draw_button(x2, y2, ExecToggleSize, (ExecLineHeight - 1),
                            toggleLightEdge, toggleDarkEdge, NULL ORTHOCGOARGVAR);
                draw_button_char(G, x2, y2 + text_lift, 'C' ORTHOCGOARGVAR);
                break;
              case 5:
                {
                  float *button_color = toggleColor2;
                  int spec_level = 0;
                  /* choose color / brightness based on specification level */
                  
                  if(rec->type == cExecAll) {
                    spec_level = MovieGetSpecLevel(G,SceneGetFrame(G));
                  } else if(rec->type == cExecObject) {
                    spec_level = ObjectGetSpecLevel(rec->obj,SceneGetFrame(G));
                  }
                  switch(spec_level) {
                    case 1:
                      button_color = toggleColor3;
                      break;
                    case 2:
                      button_color = activeColor;
                      break;
                  }
                  draw_button(x2, y2, ExecToggleSize, (ExecLineHeight - 1),
                              toggleLightEdge, toggleDarkEdge, button_color ORTHOCGOARGVAR);
                }

                draw_button_char(G, x2, y2 + text_lift, 'M' ORTHOCGOARGVAR);
                break;
              }
              x2 += ExecToggleWidth;
            }
          }
	  /* end ASHLC */
#endif

          {
            int x3 = x;
            int hidden_prefix = false;

            TextSetColor(G, TextColor);
            TextSetPos2i(G, x3 + DIP2PIXEL(2), y2 + text_lift);

            if((rec->type == cExecObject) ||
               (rec->type == cExecAll) || (rec->type == cExecSelection)) {

              y2 = y;
              x2 = xx;
              if((x - ExecToggleMargin) - (xx - ExecToggleMargin) > DIP2PIXEL(-10)) {
                x2 = x + DIP2PIXEL(10);
              }
              x3 += panel->nest_level * DIP2PIXEL(8);
              TextSetPos2i(G, x3 + DIP2PIXEL(2), y2 + text_lift);
              nChar -= panel->nest_level;

              const float* but_color = disabledColor;
              {
                int but_width = (x2 - x3) - 1;

		/* drawing a group +/- NAME */
                if(panel->is_group) {
                  const int button_width = DIP2PIXEL(15);
                  if((rec->hilight == 2) && (I->Over == I->Pressed)) {
                    draw_button(x3, y2, button_width, (ExecLineHeight - 1), lightEdge, darkEdge,
                                pressedColor ORTHOCGOARGVAR);
                  } else if(panel->is_open) {
                    draw_button(x3, y2, button_width, (ExecLineHeight - 1), lightEdge, darkEdge,
                                disabledColor ORTHOCGOARGVAR);
                  } else {
                    draw_button(x3, y2, button_width, (ExecLineHeight - 1), lightEdge, darkEdge,
                                disabledColor ORTHOCGOARGVAR);
                  }
                  TextSetPos2i(G, x3 + DIP2PIXEL(4), y2 + text_lift);
                  if(panel->is_open) {
                    TextDrawChar(G, '-' ORTHOCGOARGVAR);
                  } else {
                    TextDrawChar(G, '+' ORTHOCGOARGVAR);
                  }

                  but_width -= DIP2PIXEL(16);
                  x3 += DIP2PIXEL(16);
                  nChar -= 2;

                  TextSetPos2i(G, x3 + DIP2PIXEL(2), y2 + text_lift);
                }

                if((rec->hilight == 1) || ((row == I->Over) && (I->OverWhat == 1))) {
		  /* button hull */
                  but_color = pressedColor;
                } else if(rec->visible) {
                  int enabled = true;
                  SpecRec *group_rec = rec->group;
                  while(enabled && group_rec) {
                    if(!group_rec->visible)
                      enabled = false;
                    else
                      group_rec = group_rec->group;
                  }

                  if(enabled) {
                    but_color = enabledColor;
                  } else {
                    but_color = cloakedColor;
                  }
                }

                draw_button(x3, y2, but_width, (ExecLineHeight - 1), lightEdge,
                    darkEdge, but_color ORTHOCGOARGVAR);
              }

              TextSetColor(G, getNameColor(rec->obj, name_color_mode, but_color,
                                  TextColor));

	      /* object name */
              c = rec->name;
	      
	      /* parse out the prefix if group.foo */
              if(!full_names) {
                if(rec->group) {        /* if prefix matches group name, then truncate */
                  char *p = c, *q = rec->group->name;
                  while((*p == *q) && (*q)) {
                    p++;
                    q++;
                  }
                  if((*p) && (!*q) && (*p == '.')) {
                    hidden_prefix = true;
                    c = p;
                  }
                }
              }

	      /* wrap selection names with "(" and ")" */
              if(rec->type == cExecSelection)
                if((nChar--) > 0) {
                  TextDrawChar(G, '(' ORTHOCGOARGVAR);
                }
            }

            if(c) {
              if(hidden_prefix) {
		/* ^.name */
                if(arrows && ((nChar--) > 0)) {
                  TextDrawChar(G, '^' ORTHOCGOARGVAR);
                  TextSetPos2i(G, x3 + DIP2PIXEL(2), y2 + text_lift);
                  TextDrawChar(G, '|' ORTHOCGOARGVAR);
                }
              }

	      /* draw the object name, char by char */
              while(*c) {
                if((nChar--) > 0) {
                  TextDrawChar(G, *(c++) ORTHOCGOARGVAR);
		}
                else
                  break;
              }
            }

	    /* SELECTIONS: wrap selection names with "(" and ")" */
            if(rec->type == cExecSelection) {
              if((nChar--) > 0) {
                TextDrawChar(G, ')' ORTHOCGOARGVAR);
              }
            }

	    /* OBJECTS: output any label captions, like state number, state title */
            if(rec->type == cExecObject) {
		/* get this object's "caption" that goes on its title line,
		 * currently, this is "state-title [curState/nState]" */
              c = rec->obj->getCaption(ch, WordLength);
	      /* now print the caption */
              if(c && c[0] && nChar > 1) {
                TextSetColor(G, captionColor);
                TextSetPos2i(G, x + DIP2PIXEL(2) + DIP2PIXEL(8) * (max_char - nChar), y2 + text_lift);
                if((nChar--) > 0)
                  TextDrawChar(G, ' ' ORTHOCGOARGVAR);
                while(*c && nChar > 0) {
		    /* allow color encoding for names */
                  if(TextSetColorFromCode(G, c, captionColor)) {
		      c += 4;
                  } else {
                    TextDrawChar(G, *(c++) ORTHOCGOARGVAR);
                    --nChar;
		  }
                }
              }
            }
          }

          y -= ExecLineHeight;
          if(y < (rect.bottom))
            break;
        }
      }
    }
    I->HowFarDown = y;
  }
}


/*========================================================================*/
int ExecutiveIterateObject(PyMOLGlobals * G, pymol::CObject ** obj, void **hidden)
{
  SpecRec **rec = (SpecRec **) hidden, *I_Spec = G->Executive->Spec;
  while(ListIterate(I_Spec, (*rec), next)) {
    if((*rec)->type == cExecObject)
      break;
  }
  if(*rec)
    (*obj) = (*rec)->obj;
  else
    (*obj) = NULL;
  return ((*rec) != NULL);
}


/*========================================================================*/
int ExecutiveIterateObjectMolecule(PyMOLGlobals * G, ObjectMolecule ** obj, void **hidden)
{
  SpecRec **rec = (SpecRec **) hidden, *I_Spec = G->Executive->Spec;
  while(ListIterate(I_Spec, (*rec), next)) {
    if(((*rec)->type == cExecObject) && ((*rec)->obj->type == cObjectMolecule))
      break;
  }
  if(*rec)
    (*obj) = (ObjectMolecule *) (*rec)->obj;
  else
    (*obj) = NULL;
  return ((*rec) != NULL);
}


/*========================================================================*/
void CExecutive::reshape(int width, int height)
{
  PyMOLGlobals *G = m_G;
  CExecutive *I = G->Executive;

  Block::reshape(width, height);

  I->Width = rect.right - rect.left + 1;
  I->Height = rect.top - rect.bottom + 1;

}


/*========================================================================*/
pymol::Result<> ExecutiveReinitialize(PyMOLGlobals * G, int what, pymol::zstring_view inPattern)
{
  CExecutive *I = G->Executive;
#ifndef _PYMOL_NOPY
  int blocked = false;
#endif
  /* reinitialize PyMOL */
  const char* pattern = inPattern.data();
  if(what == 2)
    pattern = NULL;

  if(pattern && (!pattern[0]))
    pattern = NULL;
  if(!pattern) {

    switch (what) {
    case 0:                    /* everything */
      ExecutiveDelete(G, cKeywordAll);
      ColorReset(G);
      SettingInitGlobal(G, false, false, true);
      ColorUpdateFrontFromSettings(G);
      MovieReset(G);
      EditorInactivate(G);
      ControlRock(G, 0);
      OrthoReshape(G, -1, -1, false);
      MovieScenesInit(G);

#ifndef _PYMOL_NOPY
      blocked = PAutoBlock(G);
      PRunStringInstance(G, "cmd.view('*','clear')");
      PRunStringInstance(G, "cmd.config_mouse(\"three_button\")");
      WizardSet(G, NULL, false);
      PAutoUnblock(G, blocked);
#endif

      SculptCachePurge(G);
      SceneReinitialize(G);
      SelectorReinit(G);
      SeqChanged(G);
      break;
    case 1:                    /* settings */
      SettingInitGlobal(G, false, false, true);
      ExecutiveRebuildAll(G);
      break;
    case 2:                    /* store_defaults */
      SettingStoreDefault(G);
      break;
    case 3:                    /* original_settings */
      SettingInitGlobal(G, false, false, false);
      ExecutiveRebuildAll(G);
      break;
    case 4:                    /* purge_defaults */
      SettingPurgeDefault(G);
      break;
      /* reinitialize is called with what + 5 to reset internal_gui if necessary (PYMOL-1425) */
    case 5:
    case 6:
      if (G->Default){
	SettingSetGlobal_i(G, cSetting_internal_gui, SettingGet_i(G, G->Default, NULL, cSetting_internal_gui));
	SettingGenerateSideEffects(G, cSetting_internal_gui, NULL, -1, 0);
      }
      break; 
    }
    SceneUpdateStereo(G);
  } else {
    {
      CTracker *I_Tracker = I->Tracker;
      int list_id = ExecutiveGetNamesListFromPattern(G, pattern, true, true);
      int iter_id = TrackerNewIter(I_Tracker, 0, list_id);
      SpecRec *rec;

      while(TrackerIterNextCandInList(I_Tracker, iter_id, (TrackerRef **) (void *) &rec)) {
        if(rec) {
          switch (rec->type) {
          case cExecObject:
            switch (what) {
            case 0:
            case 1:
              if(rec->obj->Setting.get()) {
                ObjectPurgeSettings(rec->obj);
                rec->obj->invalidate(cRepAll, cRepInvAll, -1);
                SceneInvalidate(G);
                SeqChanged(G);
              }
              break;
            }
          }
        }
      }
      TrackerDelList(I_Tracker, list_id);
      TrackerDelIter(I_Tracker, iter_id);
    }

    /* to do */
  }

  return {};
}


/*========================================================================*/
int ExecutiveInit(PyMOLGlobals * G)
{
  CExecutive *I = NULL;
  if((I = (G->Executive = new CExecutive(G)))) {

    SpecRec *rec = NULL;

    ListInit(I->Spec);
    I->Tracker = TrackerNew(G);
    I->all_names_list_id = TrackerNewList(I->Tracker, NULL);
    I->all_obj_list_id = TrackerNewList(I->Tracker, NULL);
    I->all_sel_list_id = TrackerNewList(I->Tracker, NULL);
    I->active = true;
    OrthoAttach(G, I, cOrthoTool);
#ifndef GLUT_FULL_SCREEN
    I->oldWidth = 640;
    I->oldHeight = 480;
#endif

    I->Lex = OVLexicon_New(G->Context->heap);
    I->Key = OVOneToOne_New(G->Context->heap);

    /* create "all" entry */

    ListElemCalloc(G, rec, SpecRec);

    strcpy(rec->name, cKeywordAll);
    rec->type = cExecAll;
    rec->visible = true;
    rec->next = NULL;
    rec->cand_id = TrackerNewCand(I->Tracker, (TrackerRef *) (void *) rec);
    TrackerLink(I->Tracker, rec->cand_id, I->all_names_list_id, 1);
    ListAppend(I->Spec, rec, next, SpecRec);
    ExecutiveAddKey(I, rec);

    return 1;
  } else
    return 0;

}


/*========================================================================*/
void ExecutiveFree(PyMOLGlobals * G)
{
  CExecutive *I = G->Executive;
  SpecRec *rec = NULL;
    CGOFree(I->selIndicatorsCGO);
  while(ListIterate(I->Spec, rec, next)) {
    if(rec->type == cExecObject)
      DeleteP(rec->obj);
  }
  ListFree(I->Spec, next, SpecRec);
  if(I->Tracker)
    TrackerFree(I->Tracker);
  OVLexicon_DEL_AUTO_NULL(I->Lex);
  OVOneToOne_DEL_AUTO_NULL(I->Key);

  ExecutiveUniqueIDAtomDictInvalidate(G);

  DeleteP(G->Executive);
}

#ifdef _undefined

matrix checking code ...

double mt[3][3], mt2[3][3], pr[3][3], im[3][3], em[3][3];
printf("normalized matrix \n");
for(a = 0; a < 3; a++) {
  for(b = 0; b < 3; b++) {
    printf("%12.3f ", evect[a][b]);
  }
  printf("\n");
}

printf("\n");

printf("tensor \n");
for(a = 0; a < 3; a++) {
  for(b = 0; b < 3; b++) {
    printf("%12.3f ", mi[a][b]);
  }
  printf("\n");
}

printf("\n");

for(a = 0; a < 3; a++) {
  for(b = 0; b < 3; b++) {
    mt[a][b] = evect[a][b];
  }
}

for(a = 0; a < 3; a++) {
  for(b = 0; b < 3; b++) {
    mt2[a][b] = evect[b][a];
  }
}

matrix_multiply33d33d(mt, mt2, pr);
printf("self product 1 \n");
for(a = 0; a < 3; a++) {
  for(b = 0; b < 3; b++) {
    printf("%8.3f ", pr[a][b]);
  }
  printf("\n");
}

printf("\n");

matrix_multiply33d33d(mt, mi, im);
matrix_multiply33d33d(im, mt2, pr);
printf("diagonal product 1 \n");
for(a = 0; a < 3; a++) {
  for(b = 0; b < 3; b++) {
    printf("%8.3f ", pr[a][b]);
  }
  printf("\n");
}

printf("\n");

for(a = 0; a < 3; a++) {
  for(b = 0; b < 3; b++) {
    em[a][b] = 0.0;
  }
  em[a][a] = egval[a];
}

matrix_multiply33d33d(mt2, em, im);
matrix_multiply33d33d(im, mt, pr);
printf("diagonal product 4 \n");
for(a = 0; a < 3; a++) {
  for(b = 0; b < 3; b++) {
    printf("%8.3f ", pr[a][b]);
  }
  printf("\n");
}

printf("\n");
#endif

PyObject * ExecutiveCEAlign(PyMOLGlobals * G, PyObject * listA, PyObject * listB, int lenA, int lenB, float d0, float d1, int windowSize, int gapMax) {
#ifdef _PYMOL_NOPY
  return NULL;
#else
  int i=0;
  int smaller;
  double **dmA ,**dmB, **S;
  int bufferSize;
  pcePoint coordsA, coordsB;
  pathCache paths = NULL;
  PyObject * result;

  smaller = lenA < lenB ? lenA : lenB;
	
  /* get the coodinates from the Python objects */
  coordsA = (pcePoint) getCoords(listA, lenA);
  coordsB = (pcePoint) getCoords(listB, lenB);
	
  /* calculate the distance matrix for each protein */
  dmA = (double**) calcDM(coordsA, lenA);
  dmB = (double**) calcDM(coordsB, lenB);
	
  /* calculate the CE Similarity matrix */
  S = (double**) calcS(dmA, dmB, lenA, lenB, windowSize);
	
  /* find the best path through the CE Sim. matrix */
  bufferSize = 0;

  /* the following line HANGS PyMOL */
  paths = (pathCache) findPath(S, dmA, dmB, lenA, lenB, d0, d1, windowSize, gapMax, &bufferSize);
	
  /* Get the optimal superposition here... */
  result = (PyObject*) findBest(coordsA, coordsB, paths, bufferSize, smaller, windowSize);
	
  /* release memory */
  free(coordsA);
  free(coordsB);
  for ( i = 0; i < bufferSize; ++i )
    free(paths[i]);
  free(paths);
	
  /* distance matrices	 */
  for ( i = 0; i < lenA; i++ )
    free( dmA[i] );
  free(dmA);
	
  for  ( i = 0; i < lenB; i++ )
    free( dmB[i] );
  free(dmB);
	
  /* similarity matrix */
  for ( i = 0; i < lenA; i++ )
    free( S[i] );
  free(S);
	
  return (PyObject*) result;
#endif
}

char *ExecutiveGetObjectNames(PyMOLGlobals * G, int mode, const char *name, int enabled_only, int *numstrs){
  char *res;
  int size=0, stlen;
  CExecutive *I = G->Executive;
  CTracker *I_Tracker = I->Tracker;
  *numstrs = 0;
  {
    int list_id = ExecutiveGetNamesListFromPattern(G, name, true, true);
    int iter_id = TrackerNewIter(I_Tracker, 0, list_id);
    SpecRec *rec;
    res = VLAlloc(char, 1000);
    while(TrackerIterNextCandInList(I_Tracker, iter_id, (TrackerRef **) (void *) &rec)) {
      if((rec->type == cExecObject
	  && (((!mode) || (mode == cObjectTypeObjects) || (mode == cObjectTypePublic) || (mode == cObjectTypePublicObjects))
	      || ((rec->obj->type != cObjectGroup) && ((mode == cObjectTypePublicNonGroupObjects) || (mode == cObjectTypeNonGroupObjects)))
	      || ((rec->obj->type == cObjectGroup) && ((mode == cObjectTypePublicGroupObjects) || (mode == cObjectTypeGroupObjects)))))
	 || (rec->type == cExecSelection
	     && ((!mode) || (mode == cObjectTypeSelections) || (mode == cObjectTypePublic) || (mode == cObjectTypePublicSelections)))
	 ) {
	if((mode < 3) || (mode > 7) || (mode == cObjectTypeGroupObjects) || (rec->name[0] != '_')) {
	  if((!enabled_only) || (rec->visible)) {
	    stlen = strlen(rec->name);	    
            VLACheck(res, char, size + stlen + 1);
            strcpy(res + size, rec->name);
            size += stlen + 1;
	    *numstrs += 1;
	  }
	}
      }
    }
  }
  if (!size){
    VLAFreeP(res);
    res = NULL;
  } else {
    VLASize(res, char, size);
  }
  return (res);
}

/**
 * Get the coord set for the given object name and state index. If "omp" is
 * not NULL, then also store a pointer to the object molecule.
 */
CoordSet * ExecutiveGetCoordSet(PyMOLGlobals * G, const char * name, int state, ObjectMolecule ** omp) {
  ObjectMolecule * om = NULL;
  CoordSet * cs = NULL;

  ok_assert(1, om = ExecutiveFindObject<ObjectMolecule>(G, name));
  ok_assert(1, cs = om->getCoordSet(state));

ok_except1:
  if (omp != NULL)
    *omp = om;
  return cs;
}

#ifdef _PYMOL_IP_PROPERTIES
#endif // _PYMOL_IP_PROPERTIES


void ExecutiveUniqueIDAtomDictInvalidate(PyMOLGlobals * G) {
  CExecutive *I = G->Executive;
  if (I->m_eoo) {
    OVOneToOne_DEL_AUTO_NULL(I->m_id2eoo);
    VLAFreeP(I->m_eoo);
  }
}

const ExecutiveObjectOffset * ExecutiveUniqueIDAtomDictGet(PyMOLGlobals * G, int i) {
  CExecutive *I = G->Executive;
  OVreturn_word offset;

  if (!I->m_eoo)
    ExecutiveGetUniqueIDAtomVLADict(G, &I->m_eoo, &I->m_id2eoo);

  if(!OVreturn_IS_OK(offset = OVOneToOne_GetForward(I->m_id2eoo, i)))
    return NULL;

  return I->m_eoo + offset.word;
}

pymol::Result<> ExecutiveSetFeedbackMask(
    PyMOLGlobals* G, int action, unsigned int sysmod, unsigned char mask)
{
  switch (action) {
  case 0:
    G->Feedback->setMask(sysmod, mask);
    break;
  case 1:
    G->Feedback->enable(sysmod, mask);
    break;
  case 2:
    G->Feedback->disable(sysmod, mask);
    break;
  case 3:
    G->Feedback->push();
    break;
  case 4:
    G->Feedback->pop();
    break;
  }
  return {};
}

pymol::Result<> ExecutiveSliceNew(PyMOLGlobals* G, const char* slice_name,
    const char* map_name, int state, int map_state)
{
  int multi = false;
  pymol::CObject *obj = NULL, *mObj, *origObj;
  ObjectMap *mapObj;
  ObjectMapState *ms;

  origObj = ExecutiveFindObjectByName(G, slice_name);
  if(origObj) {
    if(origObj->type != cObjectSlice) {
      return pymol::make_error("Object ", slice_name, " is not an ObjectSlice.");
    }
  }

  mObj = ExecutiveFindObjectByName(G, map_name);
  if(mObj) {
    if(mObj->type != cObjectMap)
      mObj = NULL;
  }
  if(mObj) {
    mapObj = (ObjectMap *) mObj;
    if(state == -1) {
      multi = true;
      state = 0;
      map_state = 0;
    } else if(state == -2) {
      state = SceneGetState(G);
      if(map_state < 0)
        map_state = state;
    } else if(state == -3) {    /* append mode */
      state = 0;
      if(origObj)
        state = origObj->getNFrame();
    } else {
      if(map_state == -1) {
        map_state = 0;
        multi = true;
      } else {
        multi = false;
      }
    }
    while(1) {
      if(map_state == -2)
        map_state = SceneGetState(G);
      if(map_state == -3)
        map_state = mapObj->getNFrame() - 1;
      ms = ObjectMapStateGetActive(mapObj, map_state);
      if(ms) {
        obj = ObjectSliceFromMap(G, (ObjectSlice *) origObj, mapObj,
                                             state, map_state);

        if(!origObj) {
          ObjectSetName(obj, slice_name);
          ExecutiveManageObject(G, obj, -1, false);
        }
        PRINTFB(G, FB_ObjectMesh, FB_Actions)
          " SliceMap: created \"%s\".\n", slice_name ENDFB(G);

      } else if(!multi) {
        return pymol::make_error("State ", map_state + 1, " not present in map ", map_name);
      }
      if(multi) {
        origObj = obj;
        map_state++;
        state++;
        if(map_state >= mapObj->State.size())
          break;
      } else {
        break;
      }
    }
  } else {
    return pymol::make_error("Map or brick object ", map_name, " not found.");
  }
  return {};
}

pymol::Result<> ExecutiveLoadCoordset(
    PyMOLGlobals* G, pymol::zstring_view oname, PyObject* model, int frame)
{
  auto origObj = ExecutiveFindObjectByName(G, oname.c_str());
  if(!origObj || origObj->type != cObjectMolecule) {
    return pymol::make_error("Invalid object molecule.");
  }

  PBlock(G);
  auto obj = ObjectMoleculeLoadCoords(G, (ObjectMolecule *) origObj, model, frame);
  PUnblock(G);
  if (!obj) {
    return pymol::make_error("Load Coordset Error");
  }
  if(frame < 0)
    frame = obj->NCSet - 1;

  PRINTFB(G, FB_Executive, FB_Actions)
    " CmdLoad: Coordinates appended into object \"%s\", state %d.\n",
    oname.c_str(), frame + 1 ENDFB(G);

  return {};
}

/**
 * Discard all bonds and do distance based bonding.
 * Implementation of `cmd.rebond()`
 *
 * @param oname object name
 * @param state object state, negative values fall back to current state
 */
pymol::Result<> ExecutiveRebond(
    PyMOLGlobals* G, const char* oname, int state, bool pbc)
{
  auto obj = ExecutiveFindObjectMoleculeByName(G, oname);
  if (!obj) {
    return pymol::make_error("cannot find object");
  }

  auto cs = obj->getCoordSet(state);
  if (!cs) {
    return pymol::make_error("no such state");
  }

  ObjectMoleculeRemoveBonds(obj, 0, 0);
  ObjectMoleculeConnect(obj, cs, true, 3, pbc);
  obj->invalidate(cRepAll, cRepInvAll, -1);

  return {};
}

bool ExecutiveIsSpecRecType(
    PyMOLGlobals* G, pymol::zstring_view name, ExecRec_t execType)
{
  for (const auto& rec : pymol::make_list_adapter(G->Executive->Spec)) {
    if (rec.name == name) {
      return rec.type == execType;
    }
  }
  return false;
}

pymol::Result<> ExecutiveAddBondByIndices(PyMOLGlobals* G,
    pymol::zstring_view oname, unsigned int atm1, unsigned int atm2, int order)
{
  auto obj = ExecutiveFindObject<ObjectMolecule>(G, oname.c_str());
  if (!obj) {
    return pymol::make_error("Cannot find object ", oname);
  }
  return ObjectMoleculeAddBondByIndices(obj, atm1, atm2, order);
}

pymol::Result<> ExecutiveLoadTraj(PyMOLGlobals* G, pymol::zstring_view oname,
    pymol::zstring_view fname, int frame, int type, int interval, int average,
    int start, int stop, int max, pymol::zstring_view str1, int image,
    const float* shift, pymol::zstring_view plugin, int quiet)
{
  auto  s1 = SelectorTmp::make(G, str1.c_str());
  p_return_if_error(s1);
  bool ok = true;
  auto origObj = ExecutiveFindObjectByName(G, oname.c_str());
  if (!origObj) {
    return pymol::make_error("Must load object topology before loading trajectory.");
  }
  if (origObj->type != cObjectMolecule) {
    return pymol::make_error("Object '", oname, "' is not a molecular object.");
  }
  if ((type == cLoadTypeTRJ) && (!plugin.empty()))
    type = cLoadTypeTRJ2;
  switch (type) {
  case cLoadTypeTRJ: /* this is the ascii AMBER trajectory format... */
    PRINTFD(G, FB_CCmd) " ExecutiveLoadTraj-DEBUG: loading TRJ\n" ENDFD;
    ObjectMoleculeLoadTRJFile(G, (ObjectMolecule*) origObj, fname.c_str(), frame,
        interval, average, start, stop, max, s1->getName(), image, shift, quiet);
    PRINTFB(G, FB_Executive, FB_Actions)
        " ExecutiveLoadTraj: \"%s\" appended into object \"%s\".\n"
        " ExecutiveLoadTraj: %d total states in the object.\n",
        fname.c_str(), oname.c_str(), ((ObjectMolecule*) origObj)->NCSet ENDFB(G);
    break;
  default:
    ok = PlugIOManagerLoadTraj(G, (ObjectMolecule*) origObj, fname.c_str(), frame,
        interval, average, start, stop, max, s1->getName(), image, shift, quiet,
        plugin.c_str());
  }
  if(ok) {
    return {};
  } else {
    return pymol::make_error("Could not load trajectory");
  }
}

pymol::Result<ExecutiveRMSInfo> ExecutiveFit(PyMOLGlobals* G,
    pymol::zstring_view str1, pymol::zstring_view str2, int mode, int cutoff,
    int cycles, int quiet, pymol::zstring_view object, int state1, int state2,
    int matchmaker)
{
  SelectorTmp s1(G, str1.c_str());
  SelectorTmp s2(G, str2.c_str());
  ExecutiveRMSInfo rmsinfo;
  ExecutiveRMS(G, s1.getName(), s2.getName(), mode, cutoff, cycles, quiet,
      object.c_str(), state1, state2, false, matchmaker, &rmsinfo);
  return rmsinfo;
}

std::string ExecutiveGetGroupMemberNames(PyMOLGlobals* G, pymol::zstring_view groupName)
{
  std::string str;
  for (auto& rec : pymol::make_list_adapter(G->Executive->Spec)) {
    if (groupName == rec.group_name) {
      str += std::string(rec.name) + " ";
    }
  }
  return str;
}

/**
 * Set the `visible` property and emit `ReportEnabledChange` if the property
 * changed. If it doesn't change, do nothing.
 */
void SpecRec::setEnabled(PyMOLGlobals* G, bool enabled)
{
  if (enabled != visible) {
    visible = enabled;
    ReportEnabledChange(G, this);
  }
}

std::vector<OrderRec> ExecutiveGetOrderOf(PyMOLGlobals* G, pymol::zstring_view nameListView)
{
  auto I = G->Executive;
  std::vector<OrderRec> recs;
  for (auto& rec : ExecutiveGetSpecRecsFromPattern(G, nameListView, true, false)) {
    auto pos = ListGetPosition(I->Spec, &rec);
    recs.emplace_back(rec.name, *pos);
  }
  auto sort_by_pos =
    [](const OrderRec& a, const OrderRec& b)
    { return a.pos < b.pos; };
  std::sort(recs.begin(), recs.end(), sort_by_pos);
  return recs;
}

void ExecutiveSetOrderOf(PyMOLGlobals* G, const std::vector<OrderRec>& recs)
{
  auto I = G->Executive;
  for (const auto& oRec : recs) {
    auto rec = ExecutiveFindSpec(G, oRec.name);
    auto detachedRec = ListDetachT(I->Spec, rec);
    ListInsertAt(I->Spec, detachedRec, oRec.pos);
  }
  if (!recs.empty()) {
    ExecutiveInvalidatePanelList(G);
  }
}

pymol::Result<> ExecutiveBackgroundColor(PyMOLGlobals* G, pymol::zstring_view color)
{
  SettingSet_color(G->Setting, cSetting_bg_rgb, color.c_str());
  SettingGenerateSideEffects(G, cSetting_bg_rgb, nullptr, -1, 0);
  return {};
}
