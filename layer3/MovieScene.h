/*
 * PyMOL Scene C++ implementation
 *
 * (c) 2014 Schrodinger, Inc.
 */

#ifndef _H_MOVIESCENE
#define _H_MOVIESCENE

#include "os_python.h"
#include "Base.h"
#include "Result.h"
#include "SceneView.h"
#include "Image.h"

#include <vector>
#include <string>

struct PyMOLGlobals;

enum {
  cMovieSceneStackDefault = 0,
  cMovieSceneStackUndo = 1,
  cMovieSceneStack_SIZE
};

enum {
  STORE_VIEW   = (1 << 0),
  STORE_ACTIVE = (1 << 1),
  STORE_COLOR  = (1 << 2),
  STORE_REP    = (1 << 3),
  STORE_FRAME  = (1 << 4),
  STORE_THUMBNAIL = (1 << 5)
};

struct MovieSceneFuncArgs
{
  std::string key;
  std::string action;
  std::string message;
  bool store_view = true;
  bool store_color = true;
  bool store_active = true;
  bool store_rep = true;
  bool store_frame = true;
  bool store_thumbnail = true;
  float animate = -1.0f;
  std::string new_key;
  bool hand = true;
  std::string sele = "all";
  std::size_t stack = cMovieSceneStackDefault;
  bool quiet = false;
};

/**
 * Struct to hold scene stored atom properties
 */
class MovieSceneAtom {
public:
  int color;
  int visRep;
};

/**
 * Struct to hold scene stored object properties
 */
class MovieSceneObject {
public:
  int color;
  int visRep;
};

/**
 * Struct to hold all scene data
 */
class MovieScene {
public:

  MovieScene() = default;
  MovieScene(const MovieScene&) = default;
  MovieScene& operator=(const MovieScene&) = default;
  MovieScene(MovieScene&&) = default;
  MovieScene& operator=(MovieScene&&) = default;
  ~MovieScene() = default;

  /// bitmask, features stored in this scene
  int storemask;

  /// global state or movie frame
  int frame;

  /// text to display (with message wizard)
  std::string message;

  /// camera view
  SceneViewType view;

  /// A png buffer storing a thumbnail of the scene
  pymol::Image thumbnail;
  constexpr static int THUMBNAIL_WIDTH = 220;
  constexpr static int THUMBNAIL_HEIGHT = 124;

  /// atom properties (color, rep, etc.)
  std::map<int, MovieSceneAtom> atomdata;

  /// objects properties (enabled, color, reps, etc.)
  std::map<std::string, MovieSceneObject> objectdata;
};

/**
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

pymol::Result<> MovieSceneFunc(PyMOLGlobals* G, const MovieSceneFuncArgs& args);

pymol::Result<> MovieSceneStore(PyMOLGlobals * G, const char * name,
    const char * message,
    bool store_view,
    bool store_color,
    bool store_active,
    bool store_rep,
    bool store_frame,
    bool store_thumbnail,
    const char * sele,
    std::size_t stack,
    bool quiet);

pymol::Result<> MovieSceneRename(PyMOLGlobals * G, const char * name, const char * new_name = nullptr);

pymol::Result<> MovieSceneRecall(PyMOLGlobals * G, const char * name, float animate = -1.0,
    bool recall_view = true,
    bool recall_color = true,
    bool recall_active = true,
    bool recall_rep = true,
    bool recall_frame = true,
    const char * sele = "all",
    size_t stack = cMovieSceneStackDefault);

pymol::Result<> MovieSceneRecallImpl(PyMOLGlobals* G, const MovieScene& scene,
    float animate, bool recall_view, bool recall_color, bool recall_active,
    bool recall_rep, bool recall_frame, const char* sele);

pymol::Result<> MovieSceneRecallAllInstant(PyMOLGlobals * G, pymol::zstring_view name,
    std::size_t stack = cMovieSceneStackDefault);

pymol::Result<> MovieSceneDelete(PyMOLGlobals * G, const char * name,
    size_t stack = cMovieSceneStackDefault);

pymol::Result<> MovieSceneOrder(PyMOLGlobals * G, const char * names,
    bool sort = false,
    const char * location = nullptr /* "current" */);
pymol::Result<> MovieSceneOrder(PyMOLGlobals* G, std::vector<std::string> names,
    bool sort = false, const char* location = nullptr);

const std::vector<std::string> & MovieSceneGetOrder(PyMOLGlobals * G);
std::vector<unsigned char> MovieSceneGetThumbnail(
    PyMOLGlobals* G, pymol::zstring_view name);
std::string MovieSceneGetMessage(PyMOLGlobals* G, pymol::zstring_view name);
pymol::Result<> MovieSceneSetMessage(
    PyMOLGlobals* G, pymol::zstring_view name, pymol::zstring_view message);
void MovieScenesInit(PyMOLGlobals * G);
void MovieScenesFree(PyMOLGlobals * G);

PyObject * MovieScenesAsPyList(PyMOLGlobals * G);
void MovieScenesFromPyList(PyMOLGlobals * G, PyObject * o);

#endif
