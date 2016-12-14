/*
 * PyMOL Scene C++ implementation
 *
 * Author: Thomas Holder
 * (c) 2014 Schrodinger, Inc.
 */

#include "os_python.h"
#include "PyMOLGlobals.h"
#include "Util2.h"
#include "Movie.h"
#include "P.h"
#include "PConv.h"
#include "PConvArgs.h"
#include "Setting.h"
#include "Scene.h"
#include "View.h"
#include "Selector.h"
#include "Executive.h"
#include "Err.h"

#include <algorithm>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <vector>

#include "MovieScene.h"

enum {
  STORE_VIEW   = (1 << 0),
  STORE_ACTIVE = (1 << 1),
  STORE_COLOR  = (1 << 2),
  STORE_REP    = (1 << 3),
  STORE_FRAME  = (1 << 4)
};

void SceneSetNames(PyMOLGlobals * G, const std::vector<std::string> &list);

/*
 * Struct to hold scene stored atom properties
 */
class MovieSceneAtom {
public:
  int color;
  int visRep;
};

/*
 * Struct to hold scene stored object properties
 */
class MovieSceneObject {
public:
  int color;
  int visRep;
};

/*
 * Struct to hold all scene data
 */
class MovieScene {
public:
  // bitmask, features stored in this scene
  int storemask;

  // global state or movie frame
  int frame;

  // text to display (with message wizard)
  std::string message;

  // camera view
  SceneViewType view;

  // atom properties (color, rep, etc.)
  std::map<int, MovieSceneAtom> atomdata;

  // objects properties (enabled, color, reps, etc.)
  std::map<std::string, MovieSceneObject> objectdata;
};

/*
 * Replacement for pymol._scene_dict and pymol._scene_order
 */
class CMovieScenes {
  int scene_counter;

public:
  std::map<std::string, MovieScene> dict;
  std::vector<std::string> order;

  CMovieScenes() {
    scene_counter = 1;
  }

  std::string getUniqueKey();
};

/*
 * Get read-only pointer to G->scenes->order
 */
const std::vector<std::string> & MovieSceneGetOrder(PyMOLGlobals * G) {
  return G->scenes->order;
}

/*
 * Get a unique scene key
 */
std::string CMovieScenes::getUniqueKey()
{
  char key[16];

  for (;; ++scene_counter) {
    sprintf(key, "%03d", scene_counter);

    if (dict.find(key) == dict.end())
      return key;
  }
}

/*
 * Change the ordering of scenes
 *
 * Examples:
 * a b c d e
 * move d behind b with:    order("b d", "current")
 * a b d c e
 * move c to the top with:  order("c", "top")
 * c a b d e
 *
 * names: space separated list of names to (optionally) sort and to move to
 *        given location
 * location: current|top|botton
 */
bool MovieSceneOrder(PyMOLGlobals * G, const char * names, bool sort,
    const char * location)
{
  std::vector<std::string> names_list;
  std::vector<std::string> new_order;
  bool is_all = false;

  if (strcmp("*", names) == 0) {
    is_all = true;
    names_list = G->scenes->order;
  } else {
    names_list = strsplit(names);

    // check that all given names are existing scenes
    for (auto it = names_list.begin(); it != names_list.end(); ++it) {
      if (G->scenes->dict.find(*it) == G->scenes->dict.end()) {
        PRINTFB(G, FB_Scene, FB_Errors)
          " Error: scene '%s' is not defined.\n", it->c_str() ENDFB(G);
        return false;
      }
    }
  }

  if (names_list.empty()) {
    return true;
  }

  if (sort) {
    std::sort(names_list.begin(), names_list.end(), strlessnat);
  }

  if (is_all) {
    new_order = names_list;
  } else {
    std::set<std::string> names_set(names_list.begin(), names_list.end());

    // sanity check: unique keys?
    if (names_set.size() != names_list.size()) {
      PRINTFB(G, FB_Scene, FB_Errors)
        " Error: duplicated keys.\n" ENDFB(G);
      return false;
    }

    char loc = location ? location[0] : 'c';

    // sanity check: valid location identifier?
    if (loc != 't' && loc != 'c' && loc != 'b') {
      PRINTFB(G, FB_Scene, FB_Errors)
        " Error: invalid location '%s'.\n", location ENDFB(G);
      return false;
    }

    if (loc == 't' /* top */) {
      new_order.insert(new_order.end(), names_list.begin(), names_list.end());
    }

    for (auto it = G->scenes->order.begin(); it != G->scenes->order.end(); ++it) {
      if (!names_set.count(*it)) {
        new_order.push_back(*it);
      } else if (loc == 'c' /* current */ && *it == names_list[0]) {
        new_order.insert(new_order.end(), names_list.begin(), names_list.end());
      }
    }

    if (loc == 'b' /* bottom */) {
      new_order.insert(new_order.end(), names_list.begin(), names_list.end());
    }
  }

  G->scenes->order = new_order;
  SceneSetNames(G, G->scenes->order);

  return true;
}

/*
 * Store a scene
 *
 * name:            name (key) of the scene to store or update
 *                  ("new"/empty = unique key)
 * message:         wizard message to display with this scene
 * store_view:      store the camera view
 * store_color:     store colors
 * store_active:    store enabled/disabled
 * store_rep:       store reps
 * store_frame:     store movie frame
 */
static bool MovieSceneStore(PyMOLGlobals * G, const char * name,
    const char * message,
    bool store_view,
    bool store_color,
    bool store_active,
    bool store_rep,
    bool store_frame,
    const char * sele)
{
  std::string key(name);

  // new key?
  if (key.empty() || key == "new") {
    key = G->scenes->getUniqueKey();
    G->scenes->order.push_back(key);
  } else if (G->scenes->dict.find(key) == G->scenes->dict.end()) {
    G->scenes->order.push_back(key);
  }

  SceneSetNames(G, G->scenes->order);

  // set scene_current_name
  SettingSetGlobal_s(G, cSetting_scene_current_name, key.c_str());

  MovieScene &scene = G->scenes->dict[key];

  // storemask
  scene.storemask = (
      (store_view ? STORE_VIEW : 0) |
      (store_active ? STORE_ACTIVE : 0) |
      (store_color ? STORE_COLOR : 0) |
      (store_rep ? STORE_REP : 0) |
      (store_frame ? STORE_FRAME : 0));

  // message
  scene.message = message ? message : "";

  // camera view
  SceneGetView(G, scene.view);

  // frame
  scene.frame = SceneGetFrame(G);

  // atomdata
  if (store_color || store_rep) {

    // fill atomdata dict
    for (SeleAtomIterator iter(G, sele); iter.next();) {

      // don't store atom data for disabled objects
      if (!((CObject*)iter.obj)->Enabled)
        continue;

      AtomInfoType * ai = iter.getAtomInfo();
      int unique_id = AtomInfoCheckUniqueID(G, ai);
      MovieSceneAtom &sceneatom = scene.atomdata[unique_id];

      sceneatom.color = ai->color;
      sceneatom.visRep = ai->visRep;
    }
  }

  // objectdata
  for (ObjectIterator iter(G); iter.next();) {
    CObject * obj = iter.getObject();
    MovieSceneObject &sceneobj = scene.objectdata[obj->Name];

    sceneobj.color = obj->Color;
    sceneobj.visRep = obj->visRep;

    // store the "enabled" state in the first bit of visRep
    SET_BIT_TO(sceneobj.visRep, 0, obj->Enabled);
  }

  PRINTFB(G, FB_Scene, FB_Details)
    " scene: scene stored as \"%s\".\n", key.c_str() ENDFB(G);

  return true;
}

/*
 * Display a message (with the message wizard)
 */
static void MovieSceneRecallMessage(PyMOLGlobals * G, const std::string &message)
{
#ifndef _PYMOL_NOPY
  // we can't just call python functions because we might be in the wrong
  // thread. Instead, pass a parsable python string to PParse.
  // To avoid syntax errors, replace all single quotes in *message*.

  std::string pystr = "/cmd.scene_recall_message(r'''" + message + "''')";
  std::replace(pystr.begin() + 30, pystr.end() - 4, '\'', '`');
  PParse(G, pystr.c_str());
#endif
}

/*
 * Set the frame or state, depending on whether a movie is defined and/or
 * playing, and depending on the scene_frame_mode setting.
 */
static void MovieSceneRecallFrame(PyMOLGlobals * G, int frame)
{
  int mode = 4;

  if (MoviePlaying(G)) {
    mode = 10; // seek scene
  } else if (frame == SceneGetFrame(G)) {
    return;
  } else {
    int scene_frame_mode = SettingGetGlobal_i(G, cSetting_scene_frame_mode);
    if(scene_frame_mode == 0 || (scene_frame_mode < 0 && MovieDefined(G)))
      return;
  }

#ifdef _PYMOL_NOPY
  SceneSetFrame(G, mode, frame);
#else
  // PBlock fails with SceneSetFrame. Workaround: call from Python
  PXDecRef(PYOBJECT_CALLMETHOD(G->P_inst->cmd, "set_frame", "ii", frame + 1, mode));
#endif
}

/*
 * Recall a scene
 *
 * name:            name (key) of the scene to recall
 * animate:         animation duration, use scene_animation_duration if -1
 * store_view:      restore the camera view
 * store_color:     restore colors
 * store_active:    restore enabled/disabled
 * store_rep:       restore reps
 */
bool MovieSceneRecall(PyMOLGlobals * G, const char * name, float animate,
    bool recall_view,
    bool recall_color,
    bool recall_active,
    bool recall_rep,
    bool recall_frame,
    const char * sele)
{
  auto it = G->scenes->dict.find(name);

  if (it == G->scenes->dict.end()) {
    PRINTFB(G, FB_Scene, FB_Errors)
      " Error: scene '%s' is not defined.\n", name
      ENDFB(G);
    return false;
  }

  // set scene_current_name
  SettingSetGlobal_s(G, cSetting_scene_current_name, name);

  MovieScene &scene = it->second;

  // recall features if requested and stored
  recall_view &= (scene.storemask & STORE_VIEW) != 0;
  recall_active &= (scene.storemask & STORE_ACTIVE) != 0;
  recall_color &= (scene.storemask & STORE_COLOR) != 0;
  recall_rep &= (scene.storemask & STORE_REP) != 0;
  recall_frame &= (scene.storemask & STORE_FRAME) != 0;

  // keep track of changes
  // (obj -> repbitmask) stores the rep bits for all reps that changed,
  // a value of zero means just to invalidate color.
  std::map<CObject*, int> objectstoinvalidate;

  // atomdata
  if (recall_color || recall_rep) {

    // fill atomdata dict
    for (SeleAtomIterator iter(G, sele); iter.next();) {
      AtomInfoType * ai = iter.getAtomInfo();
      auto it = scene.atomdata.find(ai->unique_id);
      if (it == scene.atomdata.end())
        continue;

      MovieSceneAtom &sceneatom = it->second;

      if (recall_color) {
        if (ai->color != sceneatom.color)
          objectstoinvalidate[(CObject*) iter.obj];

        ai->color = sceneatom.color;
      }

      if (recall_rep) {
        int changed = (ai->visRep ^ sceneatom.visRep) & cRepsAtomMask;
        if (changed)
          objectstoinvalidate[(CObject*) iter.obj] |= changed;

        ai->visRep = sceneatom.visRep;
      }
    }
  }

  // disable all objects
  if (recall_active) {
    // need to control SpecRec
    ExecutiveSetObjVisib(G, "*", false, false);
  }

  // objectdata
  for (ObjectIterator iter(G); iter.next();) {
    CObject * obj = iter.getObject();
    auto it = scene.objectdata.find(obj->Name);
    if (it == scene.objectdata.end())
      continue;

    MovieSceneObject &sceneobj = it->second;

    if (recall_color) {
      if (obj->Color != sceneobj.color)
        objectstoinvalidate[obj];

      obj->Color = sceneobj.color;
    }

    if (recall_rep) {
      int changed = (obj->visRep ^ sceneobj.visRep) & cRepsObjectMask;
      if (changed)
        objectstoinvalidate[obj] |= changed;

      obj->visRep = sceneobj.visRep;
    }

    // "enabled" state is first bit in visRep
    int enabled = GET_BIT(sceneobj.visRep, 0);
    if(recall_active && enabled) {
      // need to control SpecRec
      ExecutiveSetObjVisib(G, obj->Name, enabled, false);
    }
  }

  // invalidate
  for (auto it = objectstoinvalidate.begin();
      it != objectstoinvalidate.end(); ++it) {
    it->first->invalidate(cRepAll, it->second ? cRepInvVisib : cRepInvColor, -1);
  }

  // camera view
  if (recall_view) {
    if (animate < -0.5) // == -1
      animate = SettingGetGlobal_f(G, cSetting_scene_animation_duration);

    SceneSetView(G, scene.view, true, animate, 1);
  }

  // message
  MovieSceneRecallMessage(G, scene.message);

  // frame
  if (recall_frame) {
    MovieSceneRecallFrame(G, scene.frame);
  }

  return true;
}

/*
 * Rename or delete a scene
 *
 * name: name to rename or delete, or "*" to delete all
 * new_name: new scene name to rename, or NULL to delete
 */
static bool MovieSceneRename(PyMOLGlobals * G, const char * name, const char * new_name = NULL) {

  if (strcmp(name, "*") == 0) {
    // delete all scenes
    G->scenes->dict.clear();
    G->scenes->order.clear();
    SceneSetNames(G, G->scenes->order);
    return true;
  }

  if (!new_name) {
    new_name = "";
  } else if (strcmp(name, new_name) == 0) {
    return true;
  }

  auto it = G->scenes->dict.find(name);

  if (it != G->scenes->dict.end()) {
    if (new_name[0])
      std::swap(G->scenes->dict[new_name], it->second);
    G->scenes->dict.erase(it);

    // does a scene named "new_name" already exist?
    auto old_new = std::find(G->scenes->order.begin(), G->scenes->order.end(), new_name);

    // replace in or remove from "G->scenes->order" list
    auto v_it = std::find(G->scenes->order.begin(), G->scenes->order.end(), name);
    if (v_it == G->scenes->order.end()) {
      printf("this is a bug, name must be in G->scenes->order\n");
    } else {
      if (new_name[0]) {
        // rename
        v_it->assign(new_name);

        // overwritten existing key?
        if (old_new != G->scenes->order.end())
          G->scenes->order.erase(old_new);
      } else {
        // remove
        G->scenes->order.erase(v_it);
      }
    }

    SceneSetNames(G, G->scenes->order);

    // update scene_current_name
    if (0 == strcmp(name, SettingGetGlobal_s(G,
            cSetting_scene_current_name))) {
      SettingSetGlobal_s(G, cSetting_scene_current_name, new_name);
    }

    return true;
  }

  return false;
}

/*
 * Print current scene order
 */
static bool MovieScenePrintOrder(PyMOLGlobals * G) {
  PRINTFB(G, FB_Scene, FB_Details)
    " scene: current order:\n" ENDFB(G);

  for (auto it = G->scenes->order.begin(); it != G->scenes->order.end(); ++it) {
    PRINTFB(G, FB_Scene, FB_Details)
      " %s", it->c_str() ENDFB(G);
  }

  PRINTFB(G, FB_Scene, FB_Details)
    "\n" ENDFB(G);

  return true;
}

/*
 * Based on the "scene_current_name" setting, get the next or previous key.
 *
 * If the "scene_loop" setting is false and the key is out of range, return
 * an empty string.
 *
 * next: true = next, false = previous
 */
static const char * MovieSceneGetNextKey(PyMOLGlobals * G, bool next) {
  const char * current_name = SettingGetGlobal_s(G, cSetting_scene_current_name);
  int scene_loop = SettingGetGlobal_b(G, cSetting_scene_loop);

  if (!current_name[0])
    scene_loop = true;

  auto it = std::find(G->scenes->order.begin(), G->scenes->order.end(), current_name);

  if (next) {
    if (it < G->scenes->order.end() - 1) {
      ++it;
    } else if (scene_loop) {
      it = G->scenes->order.begin();
    } else {
      return "";
    }
  } else {
    if (it != G->scenes->order.begin() && it != G->scenes->order.end()) {
      --it;
    } else if (scene_loop) {
      it = G->scenes->order.end() - 1;
    } else {
      return "";
    }
  }

  return it->c_str();
}

/*
 * Move the current scene (scene_current_name) before or after "key"
 */
static bool MovieSceneOrderBeforeAfter(PyMOLGlobals * G, const char * key, bool before)
{
  const char * location = NULL;
  const char * key2 = SettingGetGlobal_s(G, cSetting_scene_current_name);

  if (before) {
    auto it = std::find(G->scenes->order.begin(), G->scenes->order.end(), key);
    if (it == G->scenes->order.begin()) {
      location = "top";
      key = "";
    } else {
      key = (it - 1)->c_str();
    }
  }

  // order = key + ' ' + key2
  std::string order(key);
  order.append(" ").append(key2);

  MovieSceneOrder(G, order.c_str(), false, location);
  return true;
}

/*
 * C implementation of the "scene" command
 */
bool MovieSceneFunc(PyMOLGlobals * G, const char * key,
    const char * action,
    const char * message,
    bool store_view,
    bool store_color,
    bool store_active,
    bool store_rep,
    bool store_frame,
    float animate,
    const char * new_key,
    bool hand,
    const char * sele)
{
  std::string prev_name;
  short beforeafter = 0;
  bool status = false;

  PRINTFB(G, FB_Scene, FB_Blather)
    " MovieScene: key=%s action=%s message=%s store_view=%d store_color=%d"
    " store_active=%d store_rep=%d animate=%f new_key=%s hand=%d\n",
    key, action, message, store_view, store_color, store_active, store_rep,
    animate, new_key, hand
    ENDFB(G);

  // insert_before, insert_after
  if (strncmp(action, "insert_", 7) == 0) {
    prev_name = SettingGetGlobal_s(G, cSetting_scene_current_name);
    if (!prev_name.empty())
      beforeafter = (action[7] == 'b') ? 1 : 2;
    action = "store";
  }

  if (strcmp(action, "next") == 0 ||
      strcmp(action, "previous") == 0) {
    ok_assert(NOSCENES, G->scenes->order.size());

    key = MovieSceneGetNextKey(G, action[0] == 'n');
    action = "recall";
  } else if (strcmp(action, "start") == 0) {
    ok_assert(NOSCENES, G->scenes->order.size());

    key = G->scenes->order[0].c_str();
    action = "recall";
  } else if (strcmp(key, "auto") == 0) {
    key = SettingGetGlobal_s(G, cSetting_scene_current_name);
  }

  if (strcmp(action, "recall") == 0) {
    if (strcmp(key, "*") == 0)
      return MovieScenePrintOrder(G);

    if (!key[0]) {
      // empty key -> put up blank screen
      SettingSetGlobal_s(G, cSetting_scene_current_name, "");
      ExecutiveSetObjVisib(G, "*", false, false);
      MovieSceneRecallMessage(G, "");
    } else {
      status = MovieSceneRecall(G, key, animate, store_view, store_color,
          store_active, store_rep, store_frame, sele);
    }

  } else if (strcmp(action, "store") == 0) {
    status = MovieSceneStore(G, key, message, store_view, store_color,
        store_active, store_rep, store_frame, sele);

    // insert_before, insert_after
    if (status && beforeafter)
      status = MovieSceneOrderBeforeAfter(G, prev_name.c_str(), beforeafter == 1);

  } else if (strcmp(action, "delete") == 0) {
    status = MovieSceneRename(G, key, NULL);
  } else if (strcmp(action, "rename") == 0) {
    status = MovieSceneRename(G, key, new_key);
  } else if (strcmp(action, "order") == 0) {
    status = MovieSceneOrder(G, key);
  } else if (strcmp(action, "sort") == 0) {
    status = MovieSceneOrder(G, key, true);
  } else if (strcmp(action, "first") == 0) {
    status = MovieSceneOrder(G, key, false, "top");
  } else {
    PRINTFB(G, FB_Scene, FB_Errors)
      " Error: invalid action '%s'\n", action ENDFB(G);
  }

  // trigger GUI updates (scene buttons, Tcl/Tk menu)
  SettingSetGlobal_b(G, cSetting_scenes_changed, 1);
  SettingGenerateSideEffects(G, cSetting_scenes_changed, NULL, 0, true);

  return status;

ok_exceptNOSCENES:
  PRINTFB(G, FB_Scene, FB_Errors)
    " Error: no scenes\n" ENDFB(G);
  return false;
}

/*
 * Init/Free to set up PyMOLGlobals in PyMOL_Start
 */

void MovieScenesInit(PyMOLGlobals * G) {
  MovieScenesFree(G);
  G->scenes = new CMovieScenes;
}

void MovieScenesFree(PyMOLGlobals * G) {
  if (G->scenes) {
    delete G->scenes;
    G->scenes = NULL;
  }
}

/*
 * PConvToPyObject/PConvFromPyObject
 *
 * Convertion to/from Python objects for all MovieScene types
 */

static PyObject * PConvToPyObject(const MovieSceneAtom &v) {
  return PConvArgsToPyList(v.color, v.visRep);
}

static PyObject * PConvToPyObject(const MovieSceneObject &v) {
  return PConvArgsToPyList(v.color, v.visRep);
}

static PyObject * PConvToPyObject(const MovieScene &v) {
  PyObject * obj = PyList_New(6);
  PyList_SET_ITEM(obj, 0, PConvToPyObject(v.storemask));
  PyList_SET_ITEM(obj, 1, PConvToPyObject(v.frame));
  PyList_SET_ITEM(obj, 2, PConvToPyObject(v.message.c_str()));
  PyList_SET_ITEM(obj, 3, PConvToPyObject(v.view, cSceneViewSize));
  PyList_SET_ITEM(obj, 4, PConvToPyObject(v.atomdata));
  PyList_SET_ITEM(obj, 5, PConvToPyObject(v.objectdata));
  return obj;
}

static bool PConvFromPyObject(PyMOLGlobals *, PyObject * obj, MovieSceneAtom &out) {
  return PConvArgsFromPyList(NULL, obj, out.color, out.visRep);
}

static bool PConvFromPyObject(PyMOLGlobals *, PyObject * obj, MovieSceneObject &out) {
  return PConvArgsFromPyList(NULL, obj, out.color, out.visRep);
}

static bool PConvFromPyObject(PyMOLGlobals * G, PyObject * obj, MovieScene &out) {
  std::map<int, MovieSceneAtom> atomdata_old_ids;

  if (!G) {
    printf(" Error: G is NULL\n");
    return false;
  }

  if (!PConvArgsFromPyList(NULL, obj,
        out.storemask,
        out.frame,
        out.message,
        out.view,
        atomdata_old_ids,
        out.objectdata))
    /* ignore */;

  // restore atomdata dict but with converted ids
  PConvFromPyObject(G, PyList_GetItem(obj, 4), atomdata_old_ids);
  for (auto it = atomdata_old_ids.begin(); it != atomdata_old_ids.end(); ++it) {
    int unique_id = SettingUniqueConvertOldSessionID(G, it->first);
    std::swap(out.atomdata[unique_id], it->second);
  }

  return true;
}

/*
 * For get_session
 */

PyObject * MovieScenesAsPyList(PyMOLGlobals * G) {
  return PConvArgsToPyList(G->scenes->order, G->scenes->dict);
}

void MovieScenesFromPyList(PyMOLGlobals * G, PyObject * o) {
  // delete existing scenes
  MovieSceneRename(G, "*");

  PConvArgsFromPyList(G, o, G->scenes->order, G->scenes->dict);

  SceneSetNames(G, G->scenes->order);
}

// vi:sw=2:expandtab:cindent
