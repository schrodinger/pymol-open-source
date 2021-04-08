
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


/* 
   NOTICE:

   Important thread safety tip:

   PM operations which will ultimately call GLUT can only be called by
   the main GLUT thread (with some exceptions, such as the simple
   drawing operations which seem to be thread safe).

   Thus, pm.py needs to guard against direct invocation of certain _cmd
   (API) calls from outside threads [Examples: _cmd.png(..),
   _cmd.mpng(..) ].  Instead, it needs to hand them over to the main
   thread by way of a cmd.mdo(...)  statement.

   Note that current, most glut operations have been pushed into the
   main event and redraw loops to avoid these kinds of problems - so
   I'm not sure how important this really is anymore.

*/


/* TODO: Put in some exception handling and reporting for the
 * python calls, especially, PyArg_ParseTuple()
 */

#ifndef _PYMOL_NOPY
#define PY_SSIZE_T_CLEAN
#include"os_python.h"
#include"PyMOLGlobals.h"
#include"PyMOLOptions.h"
#include"os_predef.h"
#include"os_gl.h"
#include"os_std.h"
#include"Version.h"
#include"MemoryDebug.h"
#include"MemoryUsage.h"
#include"Err.h"
#include"Util.h"
#include"Cmd.h"
#include"ButMode.h"
#include"Ortho.h"
#include"ObjectMolecule.h"
#include"ObjectMolecule3.h"
#include"ObjectMesh.h"
#include"ObjectMap.h"
#include"ObjectCallback.h"
#include"ObjectCGO.h"
#include"ObjectSurface.h"
#include"ObjectSlice.h"
#include"Executive.h"
#include"ExecutivePython.h"
#include"Selector.h"
#include"main.h"
#include"Scene.h"
#include"SceneRay.h"
#include"Setting.h"
#include"Movie.h"
#include"P.h"
#include"PConv.h"
#include"Control.h"
#include"Editor.h"
#include"Wizard.h"
#include"SculptCache.h"
#include"TestPyMOL.h"
#include"Color.h"
#include"Seq.h"
#include"PyMOL.h"
#include"Movie.h"
#include"OVContext.h"
#include"PlugIOManager.h"
#include"Seeker.h"
#include"ListMacros.h"
#include"ObjectAlignment.h"

#include "MovieScene.h"
#include "CifFile.h"

#include "MoleculeExporter.h"

#define tmpSele "_tmp"
#define tmpSele1 "_tmp1"
#define tmpSele2 "_tmp2"

static int flush_count = 0;

#ifndef _PYMOL_NO_MAIN
static int run_only_once = true;
#endif

#define API_SETUP_PYMOL_GLOBALS \
  G = _api_get_pymol_globals(self)

#define API_SETUP_ARGS(G, self, args, ...)                                     \
  if (!PyArg_ParseTuple(args, __VA_ARGS__))                                    \
    return nullptr;                                                            \
  G = _api_get_pymol_globals(self);                                            \
  API_ASSERT(G);

/*
 * C-level tests
 */
#ifdef _PYMOL_CTEST
#include "TestCmdTest2.h"
#else
static PyObject* CmdTest2(PyObject*, PyObject*)
{
  PyErr_SetString(PyExc_NotImplementedError, "compile with --testing");
  return nullptr;
}
#endif

/**
 * Start a headless singleton instance in the current thread.
 *
 * Unlike when calling `pymol.finish_launching()`, there is no event loop,
 * so animations, continuous sculpting and modal draw are not supported.
 *
 * After calling this, SingletonPyMOLGlobals will be available.
 */
static void launch_library_singleton() {
  PyRun_SimpleString(
      "import pymol.invocation, pymol2\n"
      "pymol.invocation.parse_args(['pymol', '-cqk'])\n"
      "pymol2.SingletonPyMOL().start()");
}

/// Creating non-singleton instances disables auto-library mode
static bool auto_library_mode_disabled = false;

/**
 * Get the PyMOLGlobals pointer from the `self` object (_self._COb in Python).
 *
 * If _COb is None, launch a headless singleton ("library mode").
 */
static PyMOLGlobals * _api_get_pymol_globals(PyObject * self) {
  if(self == Py_None) {
    if (auto_library_mode_disabled) {
      PyErr_SetString(PyExc_RuntimeError, "Missing PyMOL instance");
      return nullptr;
    }

    launch_library_singleton();
    return SingletonPyMOLGlobals;
  }

  if (self && PyCapsule_CheckExact(self)) {
    auto G_handle =
        reinterpret_cast<PyMOLGlobals**>(PyCapsule_GetPointer(self, nullptr));
    if(G_handle) { \
      return *G_handle;
    } \
  }

  return NULL;
}

/**
 * Reports an error but keeps going
 */
#define API_HANDLE_ERROR \
   if (PyErr_Occurred()) PyErr_Print(); \
   fprintf(stderr,"API-Error: in %s line %d.\n",__FILE__,__LINE__);

/**
 * If `x` is false, raises CmdException("x")
 */
#define API_ASSERT(x)                                                          \
  if (!(x)) {                                                                  \
    if (!PyErr_Occurred())                                                     \
      PyErr_SetString(P_CmdException ? P_CmdException : PyExc_Exception, #x);  \
    return nullptr;                                                            \
  }

/* NOTE: the glut_thread_keep_out variable can only be changed by the thread
   holding the API lock, therefore this is safe even through increment
   isn't (necessarily) atomic. */

static void APIEnter(PyMOLGlobals * G)
{                               /* assumes API is locked */
  PRINTFD(G, FB_API)
    " APIEnter-DEBUG: as thread %ld.\n", PyThread_get_thread_ident()
    ENDFD;

  if(G->Terminating) {          /* try to bail */


/* BEGIN PROPRIETARY CODE SEGMENT (see disclaimer in "os_proprietary.h") */
#ifdef WIN32
    abort();
#endif


/* END PROPRIETARY CODE SEGMENT */
    exit(0);
  }

  if(!PIsGlutThread())
    G->P_inst->glut_thread_keep_out++;
  PUnblock(G);
}

static int APIEnterNotModal(PyMOLGlobals * G)
{                               /* assumes API is locked */
  if(PyMOL_GetModalDraw(G->PyMOL)) {
    return false;
  } else {
    APIEnter(G);
    return true;
  }
}

static void APIEnterBlocked(PyMOLGlobals * G)
{                               /* assumes API is locked */

  PRINTFD(G, FB_API)
    " APIEnterBlocked-DEBUG: as thread %ld.\n", PyThread_get_thread_ident()
    ENDFD;

  if(G->Terminating) {          /* try to bail */


/* BEGIN PROPRIETARY CODE SEGMENT (see disclaimer in "os_proprietary.h") */
#ifdef WIN32
    abort();
#endif


/* END PROPRIETARY CODE SEGMENT */
    exit(0);
  }

  if(!PIsGlutThread())
    G->P_inst->glut_thread_keep_out++;
}

static int APIEnterBlockedNotModal(PyMOLGlobals * G)
{                               /* assumes API is locked */
  if(PyMOL_GetModalDraw(G->PyMOL)) {
    return false;
  } else {
    APIEnterBlocked(G);
    return true;
  }
}

static void APIExit(PyMOLGlobals * G)
{                               /* assumes API is locked */
  PBlock(G);
  if(!PIsGlutThread())
    G->P_inst->glut_thread_keep_out--;
  PRINTFD(G, FB_API)
    " APIExit-DEBUG: as thread %ld.\n", PyThread_get_thread_ident()
    ENDFD;
}

static void APIExitBlocked(PyMOLGlobals * G)
{                               /* assumes API is locked */

  if(!PIsGlutThread())
    G->P_inst->glut_thread_keep_out--;
  PRINTFD(G, FB_API)
    " APIExitBlocked-DEBUG: as thread %ld.\n", PyThread_get_thread_ident()
    ENDFD;
}

/**
 * Return None (pymol.constants.DEFAULT_SUCCESS)
 */
static PyObject *APISuccess(void)
{                               /* success returns None */
  return PConvAutoNone(Py_None);
}

/**
 * Return -1 (pymol.constants.DEFAULT_ERROR)
 */
static PyObject *APIFailure(void)
{                               /* returns -1: a general unspecified
                                 * error */
  return (Py_BuildValue("i", -1));
}

/**
 * Raise CmdException(msg).
 */
static PyObject* APIFailure(PyMOLGlobals* G, const char* msg = nullptr)
{
  if (msg) {
    PyErr_SetString(P_CmdException, msg);
  } else {
    PyErr_SetNone(P_CmdException);
  }
  return nullptr;
}

/**
 * Raise `error` as an exception.
 */
static PyObject* APIFailure(PyMOLGlobals* G, const pymol::Error& error)
{
  if (PyErr_Occurred()) {
    return nullptr;
  }

  PyObject* exc_type;
  switch (error.code()) {
  case pymol::Error::QUIET:
    exc_type = P_QuietException;
    break;
  case pymol::Error::MEMORY:
    exc_type = PyExc_MemoryError;
    break;
  case pymol::Error::INCENTIVE_ONLY:
    exc_type = P_IncentiveOnlyException;
    break;
  default:
    exc_type = P_CmdException;
  }

  PyErr_SetString(exc_type, error.what().c_str());
  return nullptr;
}

static PyObject *APIResultCode(int code)
{                               /* innteger result code
                                 * (could be a value, a
                                 * count, or a boolean) */
  return (Py_BuildValue("i", code));
}

/**
 * If `ok` is true, return None (DEFAULT_SUCCESS).
 * Else return -1 (DEFAULT_ERROR).
 */
static PyObject *APIResultOk(int ok)
{
  if(ok)
    return APISuccess();
  else
    return APIFailure();
}

/**
 * If `ok` is true, return None (DEFAULT_SUCCESS).
 * Else raise CmdException().
 */
static PyObject* APIResultOk(PyMOLGlobals* G, bool ok)
{
  if (ok)
    return APISuccess();
  return APIFailure(G);
}

/**
 * If `res` is true, return res.result().
 * Else raise CmdException(res.error()).
 */
template <typename T>
PyObject* APIResult(PyMOLGlobals* G, pymol::Result<T>& res)
{
  if (res)
    return PConvToPyObject(res.result());
  return APIFailure(G, res.error());
}

static PyObject *APIIncRef(PyObject * result)
{
  Py_INCREF(result);
  return (result);
}

static PyObject *APIAutoNone(PyObject * result)
{                               /* automatically owned Py_None */
  if(result == Py_None)
    Py_INCREF(result);
  else if(result == NULL) {
    result = Py_None;
    Py_INCREF(result);
  }
  return (result);
}

static PyObject *CmdGetModalDraw(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  int ok = false;
  int status = false;
  ok = PyArg_ParseTuple(args, "O", &self);
  if(ok) {
    API_SETUP_PYMOL_GLOBALS;
    ok = (G != NULL);
  } else {
    API_HANDLE_ERROR;
  }
  if(ok) {
    APIEnterBlocked(G);
    status = PyMOL_GetModalDraw(G->PyMOL);
    APIExitBlocked(G);
  }
  return APIResultCode(status);
}

static PyObject *CmdPseudoatom(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  char *object_name, *sele, *label;
  char *name, *resn, *resi, *chain, *segi, *elem;
  float vdw;
  int hetatm, color;
  float b, q;
  PyObject *pos;
  int state, mode, quiet;

  API_SETUP_ARGS(G, self, args, "OssssssssfiffsOiiii", &self,
                 &object_name, &sele, &name, &resn, &resi, &chain,
                 &segi, &elem, &vdw, &hetatm, &b, &q, &label, &pos, &color,
                 &state, &mode, &quiet);
  float pos_array[3], *pos_ptr = NULL;
  if(pos && PyTuple_Check(pos) && (PyTuple_Size(pos) == 3))
    if(PyArg_ParseTuple(pos, "fff", pos_array, pos_array + 1, pos_array + 2))
      pos_ptr = pos_array;

  API_ASSERT(APIEnterBlockedNotModal(G));
  auto pseudoatom_name = ExecutivePreparePseudoatomName(G, object_name);
  auto result = ExecutivePseudoatom(G, pseudoatom_name, sele, name, resn,
      resi, chain, segi, elem, vdw, hetatm, b, q, label, pos_ptr, color, state,
      mode, quiet);
  APIExitBlocked(G);
  return APIResult(G, result);
}

static PyObject *CmdFixChemistry(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  char *str2, *str3;
  int quiet;
  int invalidate;
  API_SETUP_ARGS(G, self, args, "Ossii", &self, &str2, &str3, &invalidate, &quiet);
  API_ASSERT(APIEnterNotModal(G));
  auto result = ExecutiveFixChemistry(G, str2, str3, invalidate, quiet);
  APIExit(G);
  return APIResult(G, result);
}

static PyObject *CmdRayAntiThread(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  PyObject *py_thread_info;
  API_SETUP_ARGS(G, self, args, "OO", &self, &py_thread_info);

  auto thread_info = reinterpret_cast<CRayAntiThreadInfo*>(
      PyCapsule_GetPointer(py_thread_info, nullptr));
  API_ASSERT(thread_info);

  PUnblock(G);
  RayAntiThread(thread_info);
  PBlock(G);

  return APISuccess();
}

static PyObject *CmdRayHashThread(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  PyObject *py_thread_info;
  API_SETUP_ARGS(G, self, args, "OO", &self, &py_thread_info);

  auto thread_info = reinterpret_cast<CRayHashThreadInfo*>(
      PyCapsule_GetPointer(py_thread_info, nullptr));
  API_ASSERT(thread_info);

  PUnblock(G);
  RayHashThread(thread_info);
  PBlock(G);

  return APISuccess();
}

static PyObject *CmdRayTraceThread(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  PyObject *py_thread_info;
  API_SETUP_ARGS(G, self, args, "OO", &self, &py_thread_info);

  auto thread_info = reinterpret_cast<CRayThreadInfo*>(
      PyCapsule_GetPointer(py_thread_info, nullptr));
  API_ASSERT(thread_info);

  PUnblock(G);
  RayTraceThread(thread_info);
  PBlock(G);

  return APISuccess();
}

static PyObject *CmdCoordSetUpdateThread(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  PyObject *py_thread_info;
  API_SETUP_ARGS(G, self, args, "OO", &self, &py_thread_info);

  auto thread_info = reinterpret_cast<CCoordSetUpdateThreadInfo*>(
      PyCapsule_GetPointer(py_thread_info, nullptr));
  API_ASSERT(thread_info);

  PUnblock(G);
  CoordSetUpdateThread(thread_info);
  PBlock(G);

  return APISuccess();
}

static PyObject *CmdObjectUpdateThread(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  PyObject *py_thread_info;
  API_SETUP_ARGS(G, self, args, "OO", &self, &py_thread_info);

  auto thread_info = reinterpret_cast<CObjectUpdateThreadInfo*>(
      PyCapsule_GetPointer(py_thread_info, nullptr));
  API_ASSERT(thread_info);

  PUnblock(G);
  SceneObjectUpdateThread(thread_info);
  PBlock(G);

  return APISuccess();
}

static PyObject *CmdGetMovieLocked(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  API_SETUP_ARGS(G, self, args, "O", &self);
  return APIResultCode(MovieLocked(G));
}

static PyObject *CmdFakeDrag(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  int ok = false;
  ok = PyArg_ParseTuple(args, "O", &self);
  if(ok) {
    API_SETUP_PYMOL_GLOBALS;
    ok = (G != NULL);
  } else {
    API_HANDLE_ERROR;
  }
  if(ok) {
    PyMOL_NeedFakeDrag(G->PyMOL);
  }
  return APISuccess();
}

static PyObject *CmdDelColorection(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  int ok = false;
  PyObject *list;
  char *prefix;
  API_SETUP_ARGS(G, self, args, "OOs", &self, &list, &prefix);
  API_ASSERT(APIEnterBlockedNotModal(G));
    ok = SelectorColorectionFree(G, list, prefix);
    APIExitBlocked(G);
  return APIResultOk(G, ok);
}

static PyObject *CmdSetColorection(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  int ok = false;
  char *prefix;
  PyObject *list;
  API_SETUP_ARGS(G, self, args, "OOs", &self, &list, &prefix);
  API_ASSERT(APIEnterBlockedNotModal(G));
    ok = SelectorColorectionApply(G, list, prefix);
    APIExitBlocked(G);
  return APIResultOk(G, ok);
}

static PyObject *CmdGetColorection(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  PyObject *result = NULL;
  const char* prefix;
  API_SETUP_ARGS(G, self, args, "Os", &self, &prefix);
  API_ASSERT(APIEnterBlockedNotModal(G));
    result = SelectorColorectionGet(G, prefix);
    APIExitBlocked(G);
  return (APIAutoNone(result));
}

static PyObject *CmdSetRawAlignment(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = nullptr;
  const char *alnname;
  const char *guidename;
  PyObject *raw;
  int state = 0, quiet = 1;

  API_SETUP_ARGS(G, self, args, "sOsii" "O",
        &alnname, &raw, &guidename, &state, &quiet,
        &self);
  API_ASSERT(APIEnterNotModal(G));
  auto result = ExecutiveSetRawAlignment(G, alnname, raw, guidename, state, quiet);
  APIExit(G);
  return APIResult(G, result);
}

static PyObject* GetRawAlignment(PyMOLGlobals* G,
    const ObjectAlignment* alnobj,
    bool active_only,
    int state)
{
  if (state >= alnobj->getNFrame()) {
    PyErr_Format(PyExc_IndexError, "state %d >= NState %d", state, alnobj->getNFrame());
    return nullptr;
  }

  const auto& vla = alnobj->State[state].alignVLA;

  if (!vla) {
    PyErr_Format(PyExc_IndexError, "state %d not valid", state);
    return nullptr;
  }

  auto hide_underscore = SettingGet<bool>(G, cSetting_hide_underscore_names);
  const auto vla_len = VLAGetSize(vla);

  PyObject * raw = PyList_New(0);

  for (size_t i = 0; i < vla_len; ++i) {
    PyObject * col = PyList_New(0);

    for (int id; (id = vla[i]); ++i) {
      auto eoo = ExecutiveUniqueIDAtomDictGet(G, id);
      if (eoo
          && (!active_only || eoo->obj->Enabled)
          && (!hide_underscore || eoo->obj->Name[0] != '_')) {
        PyObject * idx = Py_BuildValue("si", eoo->obj->Name, eoo->atm + 1);
        PyList_Append(col, idx);
        Py_DECREF(idx);
      }
    }

    if (PyList_Size(col) > 0) {
      PyList_Append(raw, col);
    }

    Py_DECREF(col);
  }

  return raw;
}

static PyObject *CmdGetRawAlignment(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  const char *name;
  int active_only;
  int state = 0;
  PyObject *result = NULL;
  API_SETUP_ARGS(G, self, args, "Osi|i", &self, &name, &active_only, &state);
  APIEnterBlocked(G);
  {
    if (!name[0]) {
      name = ExecutiveGetActiveAlignment(G);
    }
    if (name && name[0]) {
      pymol::CObject *obj = ExecutiveFindObjectByName(G, name);
      if (obj && obj->type == cObjectAlignment) {
        result = GetRawAlignment(G, (ObjectAlignment*) obj, active_only, state);
      } else {
        PyErr_Format(PyExc_KeyError, "no such alignment: '%s'", name);
      }
    }
    APIExitBlocked(G);
  }
  if(!result && !PyErr_Occurred()) {
    return APIFailure(G);
  } else
    return result;
}

static PyObject *CmdGetOrigin(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  int ok = false;
  float origin[3];
  char *object;
  ok = PyArg_ParseTuple(args, "Os", &self, &object);
  if(ok) {
    API_SETUP_PYMOL_GLOBALS;
    ok = (G != NULL);
  } else {
    API_HANDLE_ERROR;
  }
  if(ok && (ok = APIEnterBlockedNotModal(G))) {
    if((!object) || (!object[0])) {
      SceneOriginGet(G, origin);
    } else {
      pymol::CObject *obj = ExecutiveFindObjectByName(G, object);
      if(!obj) {
        ok = false;
      } else {
        if(obj->TTTFlag) {
          origin[0] = -obj->TTT[12];
          origin[1] = -obj->TTT[13];
          origin[2] = -obj->TTT[14];
        } else {
          SceneOriginGet(G, origin);    /* otherwise, return scene origin */
        }
      }
    }
    APIExitBlocked(G);
  }
  if(ok) {
    return (Py_BuildValue("fff", origin[0], origin[1], origin[2]));
  } else {
    return APIFailure();
  }
}

static PyObject * CmdFindMolfilePlugin(PyObject * self, PyObject * args)
{
  PyMOLGlobals * G = NULL;
  const char * ext = NULL;
  int mask = 0;
  if (!PyArg_ParseTuple(args, "Os|i", &self, &ext, &mask)) {
    API_HANDLE_ERROR;
  } else {
    API_SETUP_PYMOL_GLOBALS;
    if (G) {
      APIEnter(G);
      const char * plugin = PlugIOManagerFindPluginByExt(G, ext, mask);
      PyObject * result = PyString_FromString(plugin ? plugin : "");
      APIExit(G);
      return APIAutoNone(result);
    }
  }
  return APIAutoNone(NULL);
}

static PyObject * CmdGetCCP4Str(PyObject * self, PyObject * args)
{
  PyMOLGlobals * G = NULL;
  const char * name = NULL;
  int state = 0;
  int quiet = 1;
  int format = cLoadTypeCCP4Unspecified;
  if (!PyArg_ParseTuple(args, "Osii|i", &self, &name, &state, &quiet, &format)) {
    API_HANDLE_ERROR;
  } else {
    API_SETUP_PYMOL_GLOBALS;
    if (G) {
      APIEnter(G);
      auto v = ObjectMapGetCCP4Str(G, name, state, quiet, format);
      PyObject * result = v.empty() ? NULL :
        PyBytes_FromStringAndSize(&v.front(), v.size());

      APIExit(G);
      return APIAutoNone(result);
    }
  }
  return APIAutoNone(NULL);
}

static PyObject * CmdGetVolumeField(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  PyObject *result = NULL;
  int state = 0;
  int ok = false;
  char* objName;
  short copy = 1;
  ok = PyArg_ParseTuple(args, "Os|ih", &self, &objName, &state, &copy);

  if(ok) {
    API_SETUP_PYMOL_GLOBALS;
    ok = (G != NULL);
  } else {
    API_HANDLE_ERROR;
  }
  if(ok && (ok = APIEnterBlockedNotModal(G))) {
    CField * field = ExecutiveGetVolumeField(G, objName, state);
    if (field) {
      result = FieldAsNumPyArray(field, copy);
    }
    APIExitBlocked(G);
  }

  if(!result) {
    return APIFailure();
  } else
    return result;
}

static PyObject * CmdGetVolumeHistogram(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  char* objName;
  float min_val = 0.f, max_val = 0.f;
  int n_points = 64;
  API_SETUP_ARGS(G, self, args, "Os|i(ff)", &self, &objName, &n_points,
      &min_val, &max_val);
  API_ASSERT(APIEnterBlockedNotModal(G));

  auto res = ExecutiveGetHistogram(G, objName, n_points, min_val, max_val);

  APIExitBlocked(G);
  return APIResult(G, res);
}

static PyObject * CmdGetVolumeRamp(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  PyObject *result = NULL;
  int ok = false;
  char* objName;
  ok = PyArg_ParseTuple(args, "Os", &self, &objName);

  if(ok) {
    API_SETUP_PYMOL_GLOBALS;
    ok = (G != NULL);
  } else {
    API_HANDLE_ERROR;
  }
  if(ok && (ok = APIEnterBlockedNotModal(G))) {
    result = ExecutiveGetVolumeRamp(G,objName);
    APIExitBlocked(G);
  }

  if(!result) {
    return APIFailure();
  } else
    return result;
}

static PyObject * CmdSetVolumeRamp(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  char* objName;
  PyObject *ramp_list;
  std::vector<float> float_array;

  API_SETUP_ARGS(G, self, args, "OsO", &self, &objName, &ramp_list);

  if (!PyList_Check(ramp_list) ||
      !PConvFromPyObject(G, ramp_list, float_array)) {
    return APIFailure(G, pymol::make_error("Invalid color array"));
  }

  API_ASSERT(APIEnterBlockedNotModal(G));
  auto result =
    ExecutiveSetVolumeRamp(G, objName, std::move(float_array));
  APIExitBlocked(G);
  return APIResult(G, result);
}

static PyObject *CmdGetVis(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  PyObject *result = NULL;
  API_SETUP_ARGS(G, self, args, "O", &self);
  API_ASSERT(APIEnterBlockedNotModal(G));
    result = ExecutiveGetVisAsPyDict(G);
    APIExitBlocked(G);
  return (APIAutoNone(result));
}

static PyObject *CmdSetVis(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  int ok = false;
  PyObject *visDict;
  API_SETUP_ARGS(G, self, args, "OO", &self, &visDict);
  API_ASSERT(APIEnterBlockedNotModal(G));
    ok = ExecutiveSetVisFromPyDict(G, visDict);
    APIExitBlocked(G);
  return APIResultOk(G, ok);
}

static PyObject *CmdReinitialize(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  int what;
  char *object;
  API_SETUP_ARGS(G, self, args, "Ois", &self, &what, &object);
  API_ASSERT(APIEnterNotModal(G));
  auto res = ExecutiveReinitialize(G, what, object);
  APIExit(G);
  return APIResult(G, res);
}

static PyObject *CmdSpectrum(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  char *str1, *expr, *prefix;
  float min, max;
  int digits, start, stop, byres;
  int quiet;
  API_SETUP_ARGS(G, self, args, "Ossffiisiii", &self, &str1, &expr, &min, &max,
      &start, &stop, &prefix, &digits, &byres, &quiet);
  API_ASSERT(APIEnterNotModal(G));

  auto res = ExecutiveSpectrum(
      G, str1, expr, min, max, start, stop, prefix, digits, byres, quiet);

  APIExit(G);
  return APIResult(G, res);
}

static PyObject *CmdMDump(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  int ok = false;
  ok = PyArg_ParseTuple(args, "O", &self);
  if(ok) {
    API_SETUP_PYMOL_GLOBALS;
    ok = (G != NULL);
  } else {
    API_HANDLE_ERROR;
  }
  if(ok) {
    APIEnter(G);
    MovieDump(G);
    APIExit(G);
  }
  return APIResultOk(ok);
}

static PyObject *CmdAccept(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  int ok = false;
  ok = PyArg_ParseTuple(args, "O", &self);
  if(ok) {
    API_SETUP_PYMOL_GLOBALS;
    ok = (G != NULL);
  } else {
    API_HANDLE_ERROR;
  }
  if(ok && (ok = APIEnterNotModal(G))) {
    MovieSetLock(G, false);
    PRINTFB(G, FB_Movie, FB_Actions)
      " Movie: Risk accepted by user.  Movie commands have been enabled.\n" ENDFB(G);
    APIExit(G);
  }
  return APIResultOk(ok);

}

static PyObject *CmdDecline(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  int ok = false;
  ok = PyArg_ParseTuple(args, "O", &self);
  if(ok) {
    API_SETUP_PYMOL_GLOBALS;
    ok = (G != NULL);
  } else {
    API_HANDLE_ERROR;
  }
  if(ok && (ok = APIEnterNotModal(G))) {
    MovieReset(G);
    PRINTFB(G, FB_Movie, FB_Actions)
      " Movie: Risk declined by user.  Movie commands have been deleted.\n" ENDFB(G);
    APIExit(G);
  }
  return APIResultOk(ok);

}

static PyObject *CmdSetSymmetry(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  char *str1, *str2;
  int state;
  float a, b, c, alpha, beta, gamma;
  int quiet;

  API_SETUP_ARGS(G, self, args, "Osiffffffsi", &self, &str1, &state, &a, &b, &c,
      &alpha, &beta, &gamma, &str2, &quiet);
  API_ASSERT(APIEnterNotModal(G));
  auto res = ExecutiveSetSymmetry(
      G, str1, state, a, b, c, alpha, beta, gamma, str2, quiet);
  APIExit(G);
  return APIResult(G, res);
}

static PyObject *CmdGetSymmetry(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  char *str1;
  float a, b, c, alpha, beta, gamma;
  int state;
  WordType sg;
  API_SETUP_ARGS(G, self, args, "Osi", &self, &str1, &state);

  APIEnter(G);
  auto res = ExecutiveGetSymmetry(
      G, str1, state, &a, &b, &c, &alpha, &beta, &gamma, sg);
  APIExit(G);

  if (!res) {
    return APIFailure(G, res.error());
  }

  if (!res.result()) {
    return APIAutoNone(nullptr);
  }

  return Py_BuildValue("[fff fff s]", a, b, c, alpha, beta, gamma, sg);
}

static PyObject *CmdSmooth(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  char *sele;
  int cycles, window, first, last, ends, quiet;
  float cutoff = -1;
  int pbc = true;
  API_SETUP_ARGS(G, self, args, "Osiiiiii|fi", &self, &sele, &cycles, &window,
      &first, &last, &ends, &quiet, &cutoff, &pbc);
  API_ASSERT(APIEnterNotModal(G));
  auto result = ExecutiveSmooth(
      G, sele, cycles, window, first, last, ends, quiet, cutoff, pbc);
  APIExit(G);
  return APIResult(G, result);
}

static PyObject *CmdGetSession(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  PyObject *dict;
  int partial, quiet;
  const char* names;
  int binary = -1;
  float version = -1.f;

  API_SETUP_ARGS(G, self, args, "OOsii|if", &self, &dict, &names, &partial,
      &quiet, &binary, &version);
  API_ASSERT(-1 <= binary && binary <= 1);

  APIEnterBlocked(G);

  const auto binary_orig = SettingGet<bool>(G, cSetting_pse_binary_dump);
  if (binary != -1)
    SettingSet(G, cSetting_pse_binary_dump, bool(binary));

  const auto version_orig = SettingGet<float>(G, cSetting_pse_export_version);
  if (version >= 0.f)
    SettingSet(G, cSetting_pse_export_version, version);

  ExecutiveGetSession(G, dict, names, partial, quiet);

  SettingSet(G, cSetting_pse_binary_dump, binary_orig);
  SettingSet(G, cSetting_pse_export_version, version_orig);

  APIExitBlocked(G);

  if (PyErr_Occurred()) {
    return nullptr;
  }

  return APISuccess();
}

static PyObject *CmdSetSession(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  int quiet, partial;
  PyObject *obj;
  API_SETUP_ARGS(G, self, args, "OOii", &self, &obj, &partial, &quiet);
  API_ASSERT(APIEnterBlockedNotModal(G));
  bool ok = ExecutiveSetSession(G, obj, partial, quiet);
  APIExitBlocked(G);
  return APIResultOk(G, ok);
}

static PyObject *CmdSetName(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  char *str1, *str2;
  API_SETUP_ARGS(G, self, args, "Oss", &self, &str1, &str2);
  API_ASSERT(APIEnterNotModal(G));
  auto result = ExecutiveSetName(G, str1, str2);
  APIExit(G);
  if (result) {
    return APISuccess();
  }
  return APIFailure(G, result.error());
}

static PyObject *CmdGetBondPrint(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  int ok = false;
  char *str1;
  int ***array = NULL;
  PyObject *result = NULL;
  int int1, int2;
  int dim[3];
  ok = PyArg_ParseTuple(args, "Osii", &self, &str1, &int1, &int2);
  if(ok) {
    API_SETUP_PYMOL_GLOBALS;
    ok = (G != NULL);
  } else {
    API_HANDLE_ERROR;
  }
  if(ok && (ok = APIEnterNotModal(G))) {
    array = ExecutiveGetBondPrint(G, str1, int1, int2, dim);
    APIExit(G);
    if(array) {
      result = PConv3DIntArrayTo3DPyList(array, dim);
      FreeP(array);
    }
  }
  return (APIAutoNone(result));
}

static PyObject *CmdDebug(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  int ok = false;
  char *str1;
  ok = PyArg_ParseTuple(args, "Os", &self, &str1);
  if(ok) {
    API_SETUP_PYMOL_GLOBALS;
    ok = (G != NULL);
  } else {
    API_HANDLE_ERROR;
  }
  if(ok && (ok = APIEnterNotModal(G))) {
    ok = ExecutiveDebug(G, str1);
    APIExit(G);
  }
  return APIResultOk(ok);
}

static PyObject *CmdSculptPurge(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  int ok = false;
  ok = PyArg_ParseTuple(args, "O", &self);
  if(ok) {
    API_SETUP_PYMOL_GLOBALS;
    ok = (G != NULL);
  } else {
    API_HANDLE_ERROR;
  }
  if(ok && (ok = APIEnterNotModal(G))) {
    SculptCachePurge(G);
    APIExit(G);
  }
  return APIResultOk(ok);
}

static PyObject *CmdScene(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = nullptr;

  MovieSceneFuncArgs margs;
  const char *key, *action, *message = NULL, *new_key = NULL;
  const char * sele = "all";

  API_SETUP_ARGS(G, self, args, "Oss|zbbbbbfzbs", &self, &key, &action,
      &message, &margs.store_view, &margs.store_color, &margs.store_active, &margs.store_rep,
      &margs.store_frame, &margs.animate, &new_key, &margs.hand, &sele);
  API_ASSERT(APIEnterBlockedNotModal(G));

  margs.key = key;
  margs.action = action;
  margs.message = message ? message : "";
  margs.new_key = new_key ? new_key : "";
  margs.sele = sele;
  auto res = MovieSceneFunc(G, margs);
  APIExitBlocked(G);
  return APIResult(G, res);
}

static PyObject *CmdSceneOrder(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;

  const char *location;
  unsigned char sort;
  PyObject* pynames = nullptr;

  API_SETUP_ARGS(G, self, args, "OObs", &self, &pynames, &sort, &location);

  std::vector<std::string> names;
  API_ASSERT(PConvFromPyObject(G, pynames, names));

  API_ASSERT(APIEnterBlockedNotModal(G));

  auto result = MovieSceneOrder(G, std::move(names), sort, location);

  APIExitBlocked(G);
  return APIResult(G, result);
}

static PyObject *CmdGetSceneOrder(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  PyObject * result = NULL;

  API_SETUP_ARGS(G, self, args, "O", &self);
  APIEnterBlocked(G);

  result = PConvToPyObject(MovieSceneGetOrder(G));

  APIExitBlocked(G);
  return (APIAutoNone(result));
}

static PyObject *CmdSculptDeactivate(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  int ok = false;
  char *str1;
  ok = PyArg_ParseTuple(args, "Os", &self, &str1);
  if(ok) {
    API_SETUP_PYMOL_GLOBALS;
    ok = (G != NULL);
  } else {
    API_HANDLE_ERROR;
  }
  if(ok && (ok = APIEnterNotModal(G))) {
    ok = ExecutiveSculptDeactivate(G, str1);
    APIExit(G);
  }
  return APIResultOk(ok);
}

static PyObject *CmdSculptActivate(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  int ok = false;
  int int1, int2, int3;
  char *str1;
  ok = PyArg_ParseTuple(args, "Osiii", &self, &str1, &int1, &int2, &int3);
  if(ok) {
    API_SETUP_PYMOL_GLOBALS;
    ok = (G != NULL);
  } else {
    API_HANDLE_ERROR;
  }
  if(ok && (ok = APIEnterNotModal(G))) {
    ok = ExecutiveSculptActivate(G, str1, int1, int2, int3);
    APIExit(G);
  }
  return APIResultOk(ok);
}

static PyObject *CmdSculptIterate(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  int ok = false;
  int int1, int2;
  char *str1;
  float total_strain = 0.0F;
  ok = PyArg_ParseTuple(args, "Osii", &self, &str1, &int1, &int2);
  if(ok) {
    API_SETUP_PYMOL_GLOBALS;
    ok = (G != NULL);
  } else {
    API_HANDLE_ERROR;
  }
  if(ok && (ok = APIEnterNotModal(G))) {
    total_strain = ExecutiveSculptIterate(G, str1, int1, int2);
    APIExit(G);
  }
  return PyFloat_FromDouble((double) total_strain);
}

static PyObject *CmdSetObjectTTT(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  float ttt[16];
  int quiet;
  char *name;
  int state;
  API_SETUP_ARGS(G, self, args, "Os(ffffffffffffffff)ii",
                            &self, &name,
                            &ttt[0], &ttt[1], &ttt[2], &ttt[3],
                            &ttt[4], &ttt[5], &ttt[6], &ttt[7],
                            &ttt[8], &ttt[9], &ttt[10], &ttt[11],
                            &ttt[12], &ttt[13], &ttt[14], &ttt[15],
                            &state, &quiet);
  API_ASSERT(APIEnterNotModal(G));
  auto result = ExecutiveSetObjectTTT(G, name, ttt, state, quiet, SettingGetGlobal_i(G, cSetting_movie_auto_store));
  APIExit(G);
  return APIResult(G, result);
}

static PyObject *CmdTranslateObjectTTT(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  float mov[3];
  char *name;
  API_SETUP_ARGS(G, self, args, "Os(fff)",
                            &self, &name,
                            &mov[0], &mov[1], &mov[2]);
  API_ASSERT(APIEnterNotModal(G));
  auto result = ExecutiveTranslateObjectTTT(G, name, mov, SettingGetGlobal_i(G, cSetting_movie_auto_store), true);
  APIExit(G);
  return APIResult(G, result);
}

static PyObject *CmdCombineObjectTTT(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  char *name;
  PyObject *m;
  float ttt[16];
  API_SETUP_ARGS(G, self, args, "OsO", &self, &name, &m);
  API_ASSERT(APIEnterNotModal(G))
  pymol::Result<> result;
    if(PConvPyListToFloatArrayInPlace(m, ttt, 16) > 0) {
    result = ExecutiveCombineObjectTTT(G, name, ttt, false, -1);
    } else {
    result = pymol::Error{"Bad Matrix"};
    }
  APIExit(G);
  return APIResult(G, result);
}

static PyObject *CmdGetColor(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  char *name;
  int mode;
  int a, nc, nvc;
  const float *rgb;
  int index;
  PyObject *result = NULL;
  API_SETUP_ARGS(G, self, args, "Osi", &self, &name, &mode);
  APIEnterBlocked(G);
  {
    switch (mode) {
    case 0:                    /* by name or index, return floats */
      index = ColorGetIndex(G, name);
      if(index >= 0) {
        rgb = ColorGet(G, index);
        result = Py_BuildValue("fff", rgb[0], rgb[1], rgb[2]);
      }
      break;
    case 1:                    /* get color names with NO NUMBERS in their names */
      nc = ColorGetNColor(G);
      nvc = 0;
      for(a = 0; a < nc; a++) {
        if(ColorGetStatus(G, a) == 1)
          nvc++;
      }
      result = PyList_New(nvc);
      nvc = 0;
      for(a = 0; a < nc; a++) {
        if(ColorGetStatus(G, a) == 1) {
          const char *color_name = ColorGetName(G, a);
          if (!color_name) {
            color_name = "";
          }
          PyObject* tup = Py_BuildValue("si", color_name, a);
          PyList_SetItem(result, nvc++, tup);
        }
      }
      break;
    case 2:                    /* get all colors */
      nc = ColorGetNColor(G);
      nvc = 0;
      for(a = 0; a < nc; a++) {
        if(ColorGetStatus(G, a) != 0)
          nvc++;
      }
      result = PyList_New(nvc);
      nvc = 0;
      for(a = 0; a < nc; a++) {
        if(ColorGetStatus(G, a)) {
          const char* color_name = ColorGetName(G, a);
          if (!color_name) {
            color_name = "";
          }
          PyObject* tup = Py_BuildValue("si", color_name, a);
          PyList_SetItem(result, nvc++, tup);
        }
      }
      break;
    case 3:                    /* get a single color index */
      result = PyInt_FromLong(ColorGetIndex(G, name));
      break;
    case 4:                    /* by name or index, return floats including negative R for special colors */
      index = ColorGetIndex(G, name);
      rgb = ColorGetSpecial(G, index);
      result = Py_BuildValue("fff", rgb[0], rgb[1], rgb[2]);
      break;
    }
    APIExitBlocked(G);
  }
  return (APIAutoNone(result));
}

static PyObject *CmdGetChains(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;

  char *str1;
  int int1;
  API_SETUP_ARGS(G, self, args, "Osi", &self, &str1, &int1);
  APIEnter(G);
  auto res = ExecutiveGetChains(G, str1, int1);
  APIExit(G);
  return APIResult(G, res);
}

static PyObject *CmdGetClickString(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = nullptr;
  int reset = 0;
  API_SETUP_ARGS(G, self, args, "O|i", &self, &reset);
  APIEnter(G);
  char* clickstr = PyMOL_GetClickString(G->PyMOL, reset);
  APIExit(G);

  if (!clickstr) {
    return APIFailure(G, "not click-ready");
  }

  PyObject* result = PyUnicode_FromString(clickstr);
  pymol::free(clickstr);
  return result;
}

static PyObject *CmdRampNew(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  char *name;
  int ok = false;
  char *map;
  int state;
  char *sele;
  float beyond, within;
  float sigma, *range_vla = NULL;
  float *color_vla = NULL;
  int zero, quiet, calc_mode = 0;
  PyObject *range, *color;
  API_SETUP_ARGS(G, self, args, "OssOOisfffii", &self, &name, &map, &range, &color,
                        &state, &sele, &beyond, &within, &sigma, &zero, &quiet);
  API_ASSERT(APIEnterNotModal(G));
  pymol::Result<> result = pymol::Error();
  ok = G != nullptr;
      if (PyList_Size(range) > 0)
        ok = PConvPyListToFloatVLA(range, &range_vla);

    if(ok) {
      if(PyList_Check(color)) {
        if (PyList_Size(color) > 0)
          ok = PConvPyList3ToFloatVLA(color, &color_vla);
      } else if(PyInt_Check(color)) {
        ok = PConvPyIntToInt(color, &calc_mode);
      }
    }
  if(ok) {
    result = ExecutiveRampNew(G, name, map,
        pymol::vla_take_ownership(range_vla),
        pymol::vla_take_ownership(color_vla), state, sele, beyond, within,
        sigma, zero, calc_mode, quiet);
  }
  APIExit(G);
  return APIResult(G, result);
}

/* Synthesize a map from possibly weighted coefficients */
/* See MapCalculate for the actual calculation */
static PyObject * CmdMapGenerate(PyObject * self, PyObject * args)
{
  PyMOLGlobals * G = NULL;
  int ok = false;

  char * name, * reflection_file, * tempFile, * amplitudes, * phases, * weights, *space_group;
  const char * cResult = NULL;
  int  quiet, zoom;
  double reso_high, reso_low, cell[6];


  ok = PyArg_ParseTuple(args, "Ossssszddsddddddii", &self, &name, &reflection_file, &tempFile,
			&amplitudes, &phases, &weights, &reso_low, &reso_high,
			&space_group, &cell[0], &cell[1], &cell[2], &cell[3], &cell[4],
			&cell[5],&quiet, &zoom);

  if(ok) {
    API_SETUP_PYMOL_GLOBALS;
    ok = (G != NULL);
  } else {
    API_HANDLE_ERROR;
  }

  if (ok && (ok = APIEnterNotModal(G))) {

    if (ok) {
      PRINTFB(G, FB_CCmd, FB_Blather)
	" Cmd-Update: Start ExecutiveMapGenerate."
      ENDFB(G);

      cResult = (const char*) ExecutiveMapGenerate(G, name, reflection_file, tempFile, 
						   amplitudes, phases, weights, 
						   reso_low, reso_high, 
						   space_group, cell, quiet, zoom);

      PRINTFB(G, FB_CCmd, FB_Blather)
	" Cmd-Update: Finished ExecutiveMapGenerate."
      ENDFB(G);
    }
    APIExit(G);
  }

  return APIAutoNone(Py_BuildValue("s", cResult));
}

static PyObject *CmdMapNew(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  char *name;
  float minCorner[3], maxCorner[3];
  float grid;
  float buffer, floor, ceiling, resolution;
  int type;
  int state;
  int have_corners;
  int quiet, zoom;
  int normalize;
  char *selection;
  API_SETUP_ARGS(G, self, args, "Osifsf(ffffff)iiiiifff",
                 &self, &name, &type, &grid, &selection, &buffer,
                        &minCorner[0], &minCorner[1], &minCorner[2],
                        &maxCorner[0], &maxCorner[1], &maxCorner[2],
                        &state, &have_corners, &quiet, &zoom, &normalize,
                        &floor, &ceiling, &resolution);
  API_ASSERT(APIEnterNotModal(G));
  auto result = ExecutiveMapNew(G, name, type, grid, selection, buffer,
      minCorner, maxCorner, state, have_corners, quiet, zoom, normalize, floor,
      ceiling, resolution);
  APIExit(G);
  return APIResult(G, result);
}

static PyObject *CmdMapSetBorder(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  char *name;
  float level;
  int state;
  int ok = false;
  ok = PyArg_ParseTuple(args, "Osfi", &self, &name, &level, &state);
  if(ok) {
    API_SETUP_PYMOL_GLOBALS;
    ok = (G != NULL);
  } else {
    API_HANDLE_ERROR;
  }
  if(ok && (ok = APIEnterNotModal(G))) {
    ok = ExecutiveMapSetBorder(G, name, level, state);
    APIExit(G);
  }
  return APIResultOk(ok);
}

static PyObject *CmdMapSet(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  char *name, *operands;
  int target_state, source_state, operator_;
  int zoom, quiet;

  API_SETUP_ARGS(G, self, args, "Osisiiii", &self, &name, &operator_, &operands,
      &target_state, &source_state, &zoom, &quiet);
  API_ASSERT(APIEnterNotModal(G));
  auto result = ExecutiveMapSet(
      G, name, operator_, operands, target_state, source_state, zoom, quiet);
  APIExit(G);
  return APIResult(G, result);
}

static PyObject *CmdMapTrim(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  char *name, *sele;
  int map_state, sele_state;
  float buffer;
  int quiet;
  API_SETUP_ARGS(G, self, args, "Ossfiii", &self, &name, &sele, &buffer,
      &map_state, &sele_state, &quiet);
  API_ASSERT(APIEnterNotModal(G));
  auto result = ExecutiveMapTrim(G, name, sele, buffer, map_state, sele_state, quiet);
  APIExit(G);
  return APIResult(G, result);
}

static PyObject *CmdMapDouble(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  char *name;
  int state;
  API_SETUP_ARGS(G, self, args, "Osi", &self, &name, &state);
  API_ASSERT(APIEnterNotModal(G));
  auto result = ExecutiveMapDouble(G, name, state);
  APIExit(G);
  return APIResult(G, result);
}

static PyObject *CmdMapHalve(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  char *name;
  int state;
  int smooth;
  API_SETUP_ARGS(G, self, args, "Osii", &self, &name, &state, &smooth);
  API_ASSERT(APIEnterNotModal(G));
  auto result = ExecutiveMapHalve(G, name, state, smooth);
  APIExit(G);
  return APIResult(G, result);
}

static PyObject *CmdGetRenderer(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  char *vendor = NULL, *renderer = NULL, *version = NULL;
  API_SETUP_ARGS(G, self, args, "O", &self);
  APIEnter(G);
    SceneGetCardInfo(G, &vendor, &renderer, &version);
    APIExit(G);
  return Py_BuildValue("(sss)", vendor, renderer, version);
}

#include <PyMOLBuildInfo.h>

static PyObject *CmdGetVersion(PyObject * self, PyObject * args)
{
  return Py_BuildValue("(sdiisi)",
      _PyMOL_VERSION,
      _PyMOL_VERSION_double,
      _PyMOL_VERSION_int,
#ifdef _PyMOL_BUILD_DATE
      _PyMOL_BUILD_DATE,
      _PYMOL_BUILD_GIT_SHA,
      0
#else
      0, "", 0
#endif
      );
}

static PyObject* CmdMemoryUsage(PyObject* self, PyObject* args)
{
  return PConvToPyObject(pymol::memory_usage());
}

static PyObject* CmdMemoryAvailable(PyObject* self, PyObject* args)
{
  return PConvToPyObject(pymol::memory_available());
}

static PyObject* CmdGetCapabilities(PyObject*, PyObject*)
{
  static PyObject* caps = nullptr;

  if (!caps) {
    caps = PySet_New(nullptr);
#ifndef _PYMOL_NO_MAIN
    // compiled with --glut
    PySet_Add(caps, PConvToPyObject("glut"));
#endif
#ifndef _PYMOL_NO_MSGPACKC
    // fast MMTF import and export
    PySet_Add(caps, PConvToPyObject("mmtf"));
#endif
#ifdef _HAVE_LIBXML
    // COLLADA export
    PySet_Add(caps, PConvToPyObject("collada"));
#endif
#ifdef _PYMOL_CTEST
    // compiled with --testing
    PySet_Add(caps, PConvToPyObject("testing"));
#endif
#ifdef _PYMOL_OPENVR
    // openvr stereo support (compiled with --openvr)
    PySet_Add(caps, PConvToPyObject("openvr"));
#endif
#ifdef _PYMOL_VMD_PLUGINS
    // VMD molfile plugins
    PySet_Add(caps, PConvToPyObject("vmdplugins"));
#endif
#ifdef _PYMOL_NUMPY
    // numpy support (cmd.get_coords, cmd.get_volume_field)
    PySet_Add(caps, PConvToPyObject("numpy"));
#endif
#ifdef _PYMOL_IP_PROPERTIES
    // object and atom level properties, incentive feature since PyMOL 1.6
    PySet_Add(caps, PConvToPyObject("properties"));
#endif
#ifndef NO_MMLIBS
    // mmlibs atom typing and stereochemistry, incentive feature PyMOL 1.4-1.8
    PySet_Add(caps, PConvToPyObject("mmlibs"));
#endif
#ifdef _PYMOL_INCENTIVE
    // Incentive PyMOL
    PySet_Add(caps, PConvToPyObject("incentive"));
#endif // _PYMOL_INCENTIVE
  }

  Py_INCREF(caps);
  return caps;
}

static PyObject *CmdTranslateAtom(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  char *str1;
  int state, log, mode;
  float v[3];
  API_SETUP_ARGS(G, self, args, "Osfffiii", &self, &str1, v, v + 1, v + 2,
      &state, &mode, &log);
  API_ASSERT(APIEnterNotModal(G));
  auto result = ExecutiveTranslateAtom(G, str1, v, state, mode, log);
  APIExit(G);
  return APIResult(G, result);
}

static PyObject *CmdMatrixCopy(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  char *source_name, *target_name;
  int source_mode, target_mode;
  int source_state, target_state, target_undo;
  int log;
  int quiet;
  API_SETUP_ARGS(G, self, args, "Ossiiiiiii", &self, &source_name, &target_name,
      &source_mode, &target_mode, &source_state, &target_state, &target_undo,
      &log, &quiet);
  API_ASSERT(APIEnterNotModal(G));
  ExecutiveMatrixCopy(G, source_name, target_name, source_mode,
      target_mode, source_state, target_state, target_undo, log, quiet);
  APIExit(G);
  return APISuccess();
}

static PyObject *CmdResetMatrix(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  char *name;
  int mode;
  int state;
  int log;
  int quiet;
  API_SETUP_ARGS(G, self, args, "Osiiii", &self, &name, &mode, &state, &log, &quiet);
  API_ASSERT(APIEnterNotModal(G));
  auto result = ExecutiveResetMatrix(G, name, mode, state, log, quiet);
  APIExit(G);
  return APIResult(G, result);
}

static PyObject *CmdTransformObject(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  char *name, *sele;
  int state, log;
  PyObject *m;
  float matrix[16];
  int homo;
  API_SETUP_ARGS(G, self, args, "OsiOisi", &self, &name, &state, &m, &log, &sele, &homo);
  API_ASSERT(APIEnterNotModal(G));
  pymol::Result<> result;
  if(PConvPyListToFloatArrayInPlace(m, matrix, 16) > 0) {
    result = ExecutiveTransformObjectSelection(G, name,
                                                 state, sele, log, matrix, homo, true);
  } else {
    result = pymol::Error{"Bad Matrix"};
  }
  APIExit(G);
  return APIResult(G, result);
}

static PyObject *CmdTransformSelection(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  char *sele;
  int state, log;
  int homo;
  PyObject *m;
  float ttt[16];
  API_SETUP_ARGS(G, self, args, "OsiOii", &self, &sele, &state, &m, &log, &homo);
  API_ASSERT(APIEnterNotModal(G));
  pymol::Result<> result;
  if(PConvPyListToFloatArrayInPlace(m, ttt, 16) > 0) {
    result = ExecutiveTransformSelection(G, state, sele, log, ttt, homo);
  } else {
    result = pymol::Error{"Bad Matrix"};
  }
  APIExit(G);
  return APIResult(G, result);
}

static PyObject *CmdLoadColorTable(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  const char* str1;
  float gamma;
  int quiet;
  API_SETUP_ARGS(G, self, args, "Osfi", &self, &str1, &gamma, &quiet);
  API_ASSERT(APIEnterNotModal(G));
  bool ok = ColorTableLoad(G, str1, gamma, quiet);
  APIExit(G);
  return APIResultOk(G, ok);
}

static PyObject *CmdLoadPNG(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  char *str1;
  int ok = false;
  int quiet;
  int movie, stereo;
  API_SETUP_ARGS(G, self, args, "Osiii", &self, &str1, &movie, &stereo, &quiet);
  API_ASSERT(APIEnterNotModal(G));
    ok = SceneLoadPNG(G, str1, movie, stereo, quiet);
    APIExit(G);
  return APIResultOk(G, ok);
}

static PyObject *CmdBackgroundColor(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  char *str1;
  API_SETUP_ARGS(G, self, args, "Os", &self, &str1);
  API_ASSERT(APIEnterNotModal(G));
  ExecutiveBackgroundColor(G, str1);
  APIExit(G);
  return APISuccess();
}

static PyObject *CmdGetPosition(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  float v[3] = { 0.0F, 0.0F, 0.0F };
  API_SETUP_ARGS(G, self, args, "O", &self);
  APIEnter(G);
    SceneGetCenter(G, v);
    APIExit(G);
  return PConvFloatArrayToPyList(v, 3);
}

static PyObject *CmdGetMoviePlaying(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  PyObject *result = NULL;
  int ok = false;
  ok = PyArg_ParseTuple(args, "O", &self);;
  if(ok) {
    API_SETUP_PYMOL_GLOBALS;
    ok = (G != NULL);
    result = PyInt_FromLong(MoviePlaying(G));
  } else {
    API_HANDLE_ERROR;
  }
  return (APIAutoNone(result));
}

static PyObject *CmdGetPhiPsi(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  char *str1;
  int state;
  PyObject *result = Py_None;
  pymol::vla<ObjectMolecule*> oVLA;
  pymol::vla<int> iVLA;
  pymol::vla<float> pVLA, sVLA;
  int l = 0;
  int a;
  API_SETUP_ARGS(G, self, args, "Osi", &self, &str1, &state);
  APIEnter(G);
  {
    l = ExecutivePhiPsi(G, str1, &oVLA, &iVLA, &pVLA, &sVLA, state);
    APIExit(G);
    if(iVLA) {
      result = PyDict_New();
      for(a = 0; a < l; a++) {
        auto key = Py_BuildValue("si", oVLA[a]->Name, iVLA[a] + 1 /* 1-based */);
        auto value = Py_BuildValue("ff", pVLA[a], sVLA[a]);
        PyDict_SetItem(result, key, value);
        Py_DECREF(key);
        Py_DECREF(value);
      }
    } else {
      result = PyDict_New();
    }
  }
  return (APIAutoNone(result));
}

static PyObject *CmdAlign(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  char *str2, *str3, *mfile, *oname;
  OrthoLineType s2 = "", s3 = "";
  int ok = false;
  int quiet, cycles, max_skip;
  float cutoff, gap, extend, seq;
  int state1, state2;
  int max_gap, transform, reset, window;
  float radius, scale, base, coord, expect, ante;
  ExecutiveRMSInfo rms_info;

  ok = PyArg_ParseTuple(args, "Ossfiffissiiiiiiffffffif", &self, &str2, &str3,
                        &cutoff, &cycles, &gap, &extend, &max_gap, &oname,
                        &mfile, &state1, &state2, &quiet, &max_skip,
                        &transform, &reset, &seq, &radius, &scale, &base,
                        &coord, &expect, &window, &ante);

  if(ok) {
    API_SETUP_PYMOL_GLOBALS;
    ok = (G != NULL);
  } else {
    API_HANDLE_ERROR;
  }
  if(ok) {
    PRINTFD(G, FB_CCmd)
      "CmdAlign-DEBUG %s %s\n", str2, str3 ENDFD;
    if((ok = APIEnterNotModal(G))) {
      ok = ((SelectorGetTmp(G, str2, s2) >= 0) && (SelectorGetTmp(G, str3, s3) >= 0));
      if(ok) {
        ok = ExecutiveAlign(G, s2, s3,
                       mfile, gap, extend, max_gap,
                       max_skip, cutoff,
                       cycles, quiet, oname, state1, state2,
                       &rms_info, transform, reset, seq,
                       radius, scale, base, coord, expect, window, ante);
      }
      SelectorFreeTmp(G, s2);
      SelectorFreeTmp(G, s3);
      APIExit(G);
    }
  }
  if(ok) {
    return Py_BuildValue("(fiififi)",
                         rms_info.final_rms,
                         rms_info.final_n_atom,
                         rms_info.n_cycles_run,
                         rms_info.initial_rms,
                         rms_info.initial_n_atom,
                         rms_info.raw_alignment_score, rms_info.n_residues_aligned);
  } else {
    return APIFailure();
  }
}

static PyObject *CmdGetCoordsAsNumPy(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  char *str1;
  int state = 0;
  OrthoLineType s1;
  PyObject *result = NULL;

  API_SETUP_ARGS(G, self, args, "Os|i", &self, &str1, &state);
  API_ASSERT(str1[0]);
  APIEnterBlocked(G);

  if(SelectorGetTmp(G, str1, s1) >= 0) {
    int sele1 = SelectorIndexByName(G, s1);
    if(sele1 >= 0) {
      result = SelectorGetCoordsAsNumPy(G, sele1, state);
    }
    SelectorFreeTmp(G, s1);
  }

  APIExitBlocked(G);
  return (APIAutoNone(result));
}

static PyObject *CmdGetCoordSetAsNumPy(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  PyObject *result = NULL;
  CoordSet *cs;
  int state = 0;
  char *name;
  short copy = 1;

  API_SETUP_ARGS(G, self, args, "Os|ih", &self, &name, &state, &copy);

  ok_assert(2, name[0]);
  ok_assert(2, state >= 0);

  APIEnterBlocked(G);

  ok_assert(1, cs = ExecutiveGetCoordSet(G, name, state));
  result = CoordSetAsNumPyArray(cs, copy);

ok_except1:
  APIExitBlocked(G);
ok_except2:
  return (APIAutoNone(result));
}

static PyObject *CmdGetSettingUpdates(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  PyObject *result = NULL;
  int state;
  char *name;

  if(!PyArg_ParseTuple(args, "Osi", &self, &name, &state)) {
    API_HANDLE_ERROR;
    ok_raise(2);
  }

  API_SETUP_PYMOL_GLOBALS;
  ok_assert(2, G && APIEnterBlockedNotModal(G));

  result = PConvToPyObject(SettingGetUpdateList(G, name, state));

  APIExitBlocked(G);
ok_except2:
  return (APIAutoNone(result));
}

static PyObject *CmdGetSettingIndices(PyObject * self, PyObject * args)
{
  return SettingGetSettingIndices();
}

static PyObject *CmdGetView(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  SceneViewType view;
  API_SETUP_ARGS(G, self, args, "O", &self);
  {
    APIEnter(G);
    SceneGetView(G, view);
    APIExit(G);
    return (Py_BuildValue
            ("(fffffffffffffffffffffffff)", view[0], view[1], view[2], view[3],
             /* 4x4 mat */
             view[4], view[5], view[6], view[7], view[8], view[9], view[10], view[11],
             view[12], view[13], view[14], view[15], view[16], view[17], view[18],
             /* pos */
             view[19], view[20], view[21],      /* origin */
             view[22], view[23],        /* clip */
             view[24]           /* orthoscopic */
            ));
  }
}

static PyObject *CmdGetViewPort(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  int width, height;
  API_SETUP_ARGS(G, self, args, "O", &self);
  APIEnter(G);
  SceneGetWidthHeight(G, &width, &height);
  APIExit(G);
  return Py_BuildValue("ii", width, height);
}

static PyObject *CmdSetView(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  SceneViewType view;
  int quiet;
  float animate;
  int hand;
  API_SETUP_ARGS(G, self, args, "O(fffffffffffffffffffffffff)ifi",
                            &self,
                            &view[0], &view[1], &view[2], &view[3],     /* 4x4 mat */
                            &view[4], &view[5], &view[6], &view[7],
                            &view[8], &view[9], &view[10], &view[11],
                            &view[12], &view[13], &view[14], &view[15],
                            &view[16], &view[17], &view[18],    /* pos */
                            &view[19], &view[20], &view[21],    /* origin */
                            &view[22], &view[23],       /* clip */
                            &view[24],  /* orthoscopic */
                            &quiet, &animate, &hand);
  API_ASSERT(APIEnterNotModal(G));
  SceneSetView(G, view, quiet, animate, hand);        /* TODO STATUS */
  APIExit(G);
  return APISuccess();
}

static PyObject *CmdGetState(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  int result = 0;
  int ok = false;
  ok = PyArg_ParseTuple(args, "O", &self);
  if(ok) {
    API_SETUP_PYMOL_GLOBALS;
    ok = (G != NULL);
  } else {
    API_HANDLE_ERROR;
  }
  if(ok) {
    result = SceneGetState(G);  /* shouldn't this be +1? */
  }
  return (APIResultCode(result));
}

static PyObject *CmdGetEditorScheme(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  int result = 0;
  int ok = false;
  ok = PyArg_ParseTuple(args, "O", &self);
  if(ok) {
    API_SETUP_PYMOL_GLOBALS;
    ok = (G != NULL);
  } else {
    API_HANDLE_ERROR;
  }
  if(ok) {
    result = EditorGetScheme(G);
  }
  return (APIResultCode(result));
}

static PyObject *CmdGetFrame(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  int result = 0;
  int ok = false;
  ok = PyArg_ParseTuple(args, "O", &self);
  if(ok) {
    API_SETUP_PYMOL_GLOBALS;
    ok = (G != NULL);
  } else {
    API_HANDLE_ERROR;
  }
  if(ok) {
    result = SceneGetFrame(G) + 1;
  }
  return (APIResultCode(result));
}

static PyObject *CmdSetTitle(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  char *str1, *str2;
  int int1;
  API_SETUP_ARGS(G, self, args, "Osis", &self, &str1, &int1, &str2);
  API_ASSERT(APIEnterNotModal(G));
  ExecutiveSetTitle(G, str1, int1, str2);
  APIExit(G);
  return APISuccess();
}

static PyObject *CmdGetTitle(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  char *str1;
  int int1;
  PyObject *result = Py_None;
  API_SETUP_ARGS(G, self, args, "Osi", &self, &str1, &int1);
  APIEnter(G);
    const char *str2 = ExecutiveGetTitle(G, str1, int1);
    APIExit(G);
    if(str2)
      result = PyString_FromString(str2);
  return (APIAutoNone(result));
}

static PyObject *CmdGetArea(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  char *str1;
  int int1, int2;
  API_SETUP_ARGS(G, self, args, "Osii", &self, &str1, &int1, &int2);
  APIEnter(G);
  auto res = ExecutiveGetArea(G, str1, int1, int2);
  APIExit(G);
  return APIResult(G, res);
}

static PyObject *CmdPushUndo(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  char *str0;
  int state;
  OrthoLineType s0 = "";
  int ok = false;
  ok = PyArg_ParseTuple(args, "Osi", &self, &str0, &state);
  if(ok) {
    API_SETUP_PYMOL_GLOBALS;
    ok = (G != NULL);
  } else {
    API_HANDLE_ERROR;
  }
  if(ok && (ok = APIEnterNotModal(G))) {
    if(str0[0])
      ok = (SelectorGetTmp(G, str0, s0) >= 0);
    if(ok)
      ok = ExecutiveSaveUndo(G, s0, state);
      if(s0[0])
	SelectorFreeTmp(G, s0);
    APIExit(G);
  }
  return APIResultOk(ok);
}

static PyObject *CmdGetType(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  char *str1;
  API_SETUP_ARGS(G, self, args, "Os", &self, &str1);
  APIEnter(G);
  auto res = ExecutiveGetType(G, str1);
  APIExit(G);
  return APIResult(G, res);
}

static PyObject *CmdGetObjectSettings(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  pymol::CObject *obj = NULL;
  PyObject *result = NULL;
  const char *oname;
  int state = -1;

  if (!PyArg_ParseTuple(args, "Os|i", &self, &oname, &state)) {
    API_HANDLE_ERROR;
    ok_raise(1);
  }

  API_SETUP_PYMOL_GLOBALS;
  ok_assert(1, G);
  APIEnterBlocked(G);

  obj = ExecutiveFindObjectByName(G, oname);

  if(!obj) {
    ErrMessage(G, "GetObjectSettings", "named object not found.");
  } else {
    auto handle = obj->getSettingHandle(-1);

    if (state != -1) {
      auto handle_state = obj->getSettingHandle(state);

      // only accept handle if different from object-level settings
      handle = (handle_state == handle) ? NULL : handle_state;
    }

    if (handle) {
      result = SettingAsPyList(handle->get(), true);
    }
  }

  APIExitBlocked(G);
ok_except1:
  return APIAutoNone(result);
}

static PyObject *CmdGetUnusedName(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  char * prefix = NULL;
  int alwaysnumber = false;
  int ok = false;
  ok = PyArg_ParseTuple(args, "Osi", &self, &prefix, &alwaysnumber);
  if(ok) {
    API_SETUP_PYMOL_GLOBALS;
    ok = (G != NULL);
  } else {
    API_HANDLE_ERROR;
  }
  if (ok && (ok = APIEnterNotModal(G))) {
    auto result = PConvToPyObject(ExecutiveGetUnusedName(G, prefix, alwaysnumber));
    APIExit(G);
    return result;
  } else {
    return APIResultOk(ok);
  }
}

static PyObject *CmdGetDragObjectName(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  PyObject *result = Py_None;
  int ok = false;
  ok = PyArg_ParseTuple(args, "O", &self);
  if(ok) {
    API_SETUP_PYMOL_GLOBALS;
    ok = (G != NULL);
  } else {
    API_HANDLE_ERROR;
  }
  if(ok) {
    char *name = NULL;
    APIEnter(G);
    {
      pymol::CObject *obj = EditorDragObject(G);    
      if(obj)
        name = obj->Name;
    }
    APIExit(G);
    if(name) 
      result = PyString_FromString(name);
    else
      result = PyString_FromString("");
  }
  return (APIAutoNone(result));
}

static PyObject *CmdGetLegalName(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  WordType name;
  char *str0;
  API_SETUP_ARGS(G, self, args, "Os", &self, &str0);
    APIEnter(G);
    UtilNCopy(name, str0, sizeof(WordType));
  ObjectMakeValidName(G, name, true /* quiet */);
    APIExit(G);
  return PyString_FromString(name);
}

static PyObject *CmdGetNames(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  int int1, int2;
  char *str0;
  API_SETUP_ARGS(G, self, args, "Oiis", &self, &int1, &int2, &str0);
  APIEnter(G);
  auto res = ExecutiveGetNames(G, int1, int2, str0);
  APIExit(G);
  return APIResult(G, res);
}

static PyObject *CmdInterrupt(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  int int1;
  int ok = false;
  ok = PyArg_ParseTuple(args, "Oi", &self, &int1);
  if(ok) {
    API_SETUP_PYMOL_GLOBALS;
    ok = (G != NULL);
  } else {
    API_HANDLE_ERROR;
  }
  if(ok) {
    PyMOL_SetInterrupt(G->PyMOL, int1);
  }
  return APIResultOk(ok);
}

static PyObject *CmdInvert(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  int int1;
  API_SETUP_ARGS(G, self, args, "Oi", &self, &int1);
  API_ASSERT(APIEnterNotModal(G));
  auto res = EditorInvert(G, int1);
    APIExit(G);
  return APIResult(G, res);
}

static PyObject *CmdTorsion(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  float float1;
  API_SETUP_ARGS(G, self, args, "Of", &self, &float1);
  API_ASSERT(APIEnterNotModal(G));
  auto result = EditorTorsion(G, float1);
  APIExit(G);
  return APIResult(G, result);
}

static PyObject *CmdUndo(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  int int1;
  int ok = false;
  ok = PyArg_ParseTuple(args, "Oi", &self, &int1);
  if(ok) {
    API_SETUP_PYMOL_GLOBALS;
    ok = (G != NULL);
  } else {
    API_HANDLE_ERROR;
  }
  if(ok && (ok = APIEnterNotModal(G))) {
    ExecutiveUndo(G, int1);     /* TODO STATUS */
    APIExit(G);
  }
  return APIResultOk(ok);
}

static PyObject *CmdMask(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  char *str1;
  int int1, quiet;

  API_SETUP_ARGS(G, self, args, "Osii", &self, &str1, &int1, &quiet);
  API_ASSERT(APIEnterNotModal(G));
  auto result = ExecutiveMask(G, str1, int1, quiet);
  APIExit(G);
  return APIResult(G, result);
}

static PyObject *CmdProtect(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  char *str1;
  int int1, int2;

  API_SETUP_ARGS(G, self, args, "Osii", &self, &str1, &int1, &int2);
  API_ASSERT(APIEnterNotModal(G));
  auto result = ExecutiveProtect(G, str1, int1, int2);
  APIExit(G);
  return APIResult(G, result);
}

static PyObject *CmdButton(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  int i1, i2;
  int ok = false;
  ok = PyArg_ParseTuple(args, "Oii", &self, &i1, &i2);
  if(ok) {
    API_SETUP_PYMOL_GLOBALS;
    ok = (G != NULL);
  } else {
    API_HANDLE_ERROR;
  }
  if(ok && (ok = APIEnterNotModal(G))) {
    ButModeSet(G, i1, i2);      /* TODO STATUS */
    APIExit(G);
  }
  return APIResultOk(ok);
}

static PyObject *CmdFeedback(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  int i1, i2, result = 0;
  int ok = false;
  ok = PyArg_ParseTuple(args, "Oii", &self, &i1, &i2);
  if(ok) {
    API_SETUP_PYMOL_GLOBALS;
    ok = (G != NULL);
  } else {
    API_HANDLE_ERROR;
  }
  if(ok) {
    /* NO API Entry for performance,
     *feedback (MACRO) just accesses a safe global */
    result = Feedback(G, i1, i2);
  }
  return Py_BuildValue("i", result);
}

static PyObject *CmdSetFeedbackMask(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  int i1, i2, i3;
  API_SETUP_ARGS(G, self, args, "Oiii", &self, &i1, &i2, &i3);
  API_ASSERT(APIEnterNotModal(G));
  ExecutiveSetFeedbackMask(
      G, i1, static_cast<unsigned char>(i2), static_cast<unsigned char>(i3));
  APIExit(G);
  return APISuccess();
}

static PyObject *CmdPop(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  char *str1, *str2;
  int quiet;
  int result = 0;
  int ok = false;
  ok = PyArg_ParseTuple(args, "Ossi", &self, &str1, &str2, &quiet);
  if(ok) {
    API_SETUP_PYMOL_GLOBALS;
    ok = (G != NULL);
  } else {
    API_HANDLE_ERROR;
  }
  if(ok && (ok = APIEnterNotModal(G))) {
    result = ExecutivePop(G, str1, str2, quiet);
    APIExit(G);
  } else
    result = -1;
  return (APIResultCode(result));

}

static PyObject *CmdFlushNow(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  int ok = false;
  ok = PyArg_ParseTuple(args, "O", &self);;
  if(ok) {
    API_SETUP_PYMOL_GLOBALS;
    ok = (G != NULL);
  } else {
    API_HANDLE_ERROR;
  }
  if(ok && G->Ready) {
    /* only called by the GLUT thread with unlocked API, blocked interpreter */
    if(flush_count < 8) {       /* prevent super-deep recursion */
      flush_count++;
      PFlushFast(G);
      flush_count--;
    } else {
      PRINTFB(G, FB_CCmd, FB_Warnings)
        " Cmd: PyMOL lagging behind API requests...\n" ENDFB(G);
    }
  }
  return APISuccess();
}

static PyObject *CmdWaitQueue(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  PyObject *result = NULL;
  int ok = false;
  ok = PyArg_ParseTuple(args, "O", &self);;

  if(ok) {
    API_SETUP_PYMOL_GLOBALS;
    ok = (G != NULL);
  } else {
    API_HANDLE_ERROR;
  }
  if(ok) {
    int waiting = 0;

    /* called by non-GLUT thread with unlocked API, blocked interpreter */
    if(!G->Terminating) {
      if(APIEnterBlockedNotModal(G)) {
        if(OrthoCommandWaiting(G)
           || (flush_count > 1))
          waiting = 1;          /* commands are waiting or we're in nested execution */
        APIExitBlocked(G);
      } else {
        waiting = 1;            /* we're performing a "modal" task... */
      }
    } else {
      waiting = 1;              /* we're shutting down... */
    }
    result = PyInt_FromLong(waiting);
  }
  return APIAutoNone(result);
}

static PyObject *CmdWaitDeferred(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  PyObject *result = NULL;
  int ok = false;
  ok = PyArg_ParseTuple(args, "O", &self);;
  if(ok) {
    API_SETUP_PYMOL_GLOBALS;
    ok = (G != NULL);
  } else {
    API_HANDLE_ERROR;
  }
  if(ok) {
    if(!G->Terminating) {
      if(APIEnterBlockedNotModal(G)) {
        if(OrthoDeferredWaiting(G))
          result = PyInt_FromLong(1);
        else
          result = PyInt_FromLong(0);
        APIExitBlocked(G);
      }
    }
  }
  return APIAutoNone(result);
}

static PyObject *CmdPaste(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  PyObject *list, *str;
  const char *st;
  int l, a;
  int ok = false;
  ok = PyArg_ParseTuple(args, "OO", &self, &list);
  if(ok) {
    API_SETUP_PYMOL_GLOBALS;
    ok = (G != NULL);
  } else {
    API_HANDLE_ERROR;
  }
  if(ok) {
    if(!list)
      ok = false;
    else if(!PyList_Check(list))
      ok = false;
    else {
      l = PyList_Size(list);
      for(a = 0; a < l; a++) {
        str = PyList_GetItem(list, a);
        if(str) {
          if(PyString_Check(str)) {
            st = PyString_AsString(str);
            if((ok = APIEnterNotModal(G))) {
              OrthoPasteIn(G, st);
              if(a < (l - 1))
                OrthoPasteIn(G, "\n");
              APIExit(G);
            }
          } else {
            ok = false;
          }
        }
      }
    }
  }
  return APIResultOk(ok);
}

static PyObject *CmdGetVRML(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  PyObject *result = NULL;
  int ver;
  API_SETUP_ARGS(G, self, args, "Oi", &self, &ver);
  {
    char *vla = NULL;
    API_ASSERT(APIEnterNotModal(G));
      SceneRay(G, 0, 0, (ver == 1) ? 6 : 4,     /* VRML1 or 2? */
               NULL, &vla, 0.0F, 0.0F, false, NULL, false, -1);
      APIExit(G);
    if(vla) {
      result = Py_BuildValue("s", vla);
    }
    VLAFreeP(vla);
  }
  return (APIAutoNone(result));
}

/**
 * Return a COLLADA string or None on failure
 */
static PyObject *CmdGetCOLLADA(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  PyObject *result = NULL;
  int ver;
  char *vla = NULL;

  API_SETUP_ARGS(G, self, args, "Oi", &self, &ver);
  API_ASSERT(APIEnterNotModal(G));

  SceneRay(G, 0, 0, 8,     /* mode 8 = COLLADA */
      NULL, &vla, 0.0F, 0.0F, false, NULL, false, -1);
  APIExit(G);

  if (vla && vla[0]) {
    result = Py_BuildValue("s", vla);
  }

  VLAFreeP(vla);

  return APIAutoNone(result);
}


static PyObject *CmdGetIdtf(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  PyObject *result = NULL;
  API_SETUP_ARGS(G, self, args, "O", &self);
  {
    char *node = NULL, *rsrc = NULL;
    API_ASSERT(APIEnterNotModal(G));
      SceneRay(G, 0, 0, cSceneRay_MODE_IDTF,
               &node, &rsrc, 0.0F, 0.0F, false, NULL, false, -1);
      APIExit(G);
    if(node && rsrc) {
      result = Py_BuildValue("(ss)", node, rsrc);
    }
    VLAFreeP(node);
    VLAFreeP(rsrc);
  }
  return (APIAutoNone(result));
}

static PyObject *CmdGetPovRay(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  PyObject *result = NULL;
  API_SETUP_ARGS(G, self, args, "O", &self);
  {
    char *header = NULL, *geom = NULL;
    API_ASSERT(APIEnterNotModal(G));
      SceneRay(G, 0, 0, 1, &header, &geom, 0.0F, 0.0F, false, NULL, false, -1);
      APIExit(G);
    if(header && geom) {
      result = Py_BuildValue("(ss)", header, geom);
    }
    VLAFreeP(header);
    VLAFreeP(geom);
  }
  return (APIAutoNone(result));
}

static PyObject *CmdGetMtlObj(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  PyObject *result = NULL;
  API_SETUP_ARGS(G, self, args, "O", &self);
  {
    char *obj = NULL, *mtl = NULL;
    API_ASSERT(APIEnterNotModal(G));
      SceneRay(G, 0, 0, 5, &obj, &mtl, 0.0F, 0.0F, false, NULL, false, -1);
      APIExit(G);
    if(obj && mtl) {
      result = Py_BuildValue("(ss)", mtl, obj);
    }
    VLAFreeP(obj);
    VLAFreeP(mtl);
  }
  return (APIAutoNone(result));
}

static PyObject *CmdGetWizard(PyObject * self, PyObject * args)
{
  PyMOLGlobals* G = nullptr;
  API_SETUP_ARGS(G, self, args, "O", &self);
  API_ASSERT(APIEnterNotModal(G));
  auto result = WizardGet(G);
  APIExit(G);
  PyObject* res = result ? result : Py_None;
  return APIIncRef(res);
}

static PyObject *CmdGetWizardStack(PyObject * self, PyObject * args)
{
  PyMOLGlobals* G = nullptr;
  API_SETUP_ARGS(G, self, args, "O", &self);
  API_ASSERT(APIEnterBlockedNotModal(G));
  auto result = WizardGetStack(G);
  APIExitBlocked(G);
  return result;
}

static PyObject *CmdSetWizard(PyObject * self, PyObject * args)
{
  PyMOLGlobals* G = nullptr;
  PyObject* obj;
  int replace;
  API_SETUP_ARGS(G, self, args, "OOi", &self, &obj, &replace);
  API_ASSERT(APIEnterNotModal(G));
  if (!obj) {
    return APIFailure(G, "Invalid wizard.");
  }
  auto result = WizardSet(G, obj, replace);
  APIExit(G);
  return APIResult(G, result);
}

static PyObject *CmdSetWizardStack(PyObject * self, PyObject * args)
{
  PyMOLGlobals* G = nullptr;
  PyObject* obj;
  API_SETUP_ARGS(G, self, args, "OO", &self, &obj);
  API_ASSERT(APIEnterNotModal(G));
  if(!obj) {
    return APIFailure(G, "Invalid wizard.");
  }
  auto result = WizardSetStack(G, obj);
  APIExit(G);
  return APIResult(G, result);
}

static PyObject *CmdRefreshWizard(PyObject * self, PyObject * args)
{
  PyMOLGlobals* G = nullptr;
  API_SETUP_ARGS(G, self, args, "O", &self);
  API_ASSERT(APIEnterNotModal(G));
  WizardRefresh(G);
  OrthoInvalidateDoDraw(G);
  OrthoDirty(G);
  APIExit(G);
  return APISuccess();
}

static PyObject *CmdDirtyWizard(PyObject * self, PyObject * args)
{
  PyMOLGlobals* G = nullptr;
  API_SETUP_ARGS(G, self, args, "O", &self);
  API_ASSERT(APIEnterNotModal(G));
  WizardDirty(G);
  APIExit(G);
  return APISuccess();
}

static PyObject *CmdSplash(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  int query;
  int result = 1;
  int ok = false;
  ok = PyArg_ParseTuple(args, "Oi", &self, &query);
  if(ok) {
    API_SETUP_PYMOL_GLOBALS;
    ok = (G != NULL);
  } else {
    API_HANDLE_ERROR;
  }
  if(!query) {
    if(ok && (ok = APIEnterNotModal(G))) {
      OrthoSplash(G);
      APIExit(G);
    }
  } else {
    /* BEGIN PROPRIETARY CODE SEGMENT (see disclaimer in "os_proprietary.h") */
#ifdef PYMOL_EVAL
    result = 2;
#else
#ifdef _PYMOL_IP_EXTRAS
    result = 0;
#endif
#endif
    /* END PROPRIETARY CODE SEGMENT */
  }
  return APIResultCode(result);
}

static PyObject *CmdCls(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  int ok = false;
  ok = PyArg_ParseTuple(args, "O", &self);
  if(ok) {
    API_SETUP_PYMOL_GLOBALS;
    ok = (G != NULL);
  } else {
    API_HANDLE_ERROR;
  }
  if(ok && (ok = APIEnterNotModal(G))) {
    OrthoClear(G);
    APIExit(G);
  }
  return APISuccess();
}

static PyObject *CmdDump(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  char *str1, *str2;
  int state, quiet;
  int ok = false;
  ok = PyArg_ParseTuple(args, "Ossii", &self, &str1, &str2, &state, &quiet);
  if(ok) {
    API_SETUP_PYMOL_GLOBALS;
    ok = (G != NULL);
  } else {
    API_HANDLE_ERROR;
  }
  if(ok && (ok = APIEnterNotModal(G))) {
    ExecutiveDump(G, str1, str2, state, quiet);       /* TODO STATUS */
    APIExit(G);
  }
  return APIResultOk(ok);
}

static PyObject *CmdIsomesh(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  char *mesh_name, *map_name, *sele;
  float lvl, fbuf, alt_lvl;
  int mesh_mode;
  int state = -1;
  float carve;
  int map_state;
  int quiet;
  API_SETUP_ARGS(G, self, args, "Osssffiifiif", &self, &mesh_name, &map_name,
      &sele, &fbuf, &lvl, &mesh_mode, &state, &carve, &map_state, &quiet,
      &alt_lvl);
  API_ASSERT(APIEnterNotModal(G));
  auto res = ExecutiveIsomeshEtc(G, mesh_name, map_name, lvl, sele, fbuf, state,
      carve, map_state, quiet, mesh_mode, alt_lvl);
  APIExit(G);
  return APIResult(G, res);
}

static PyObject *CmdSliceNew(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  char *slice;
  char *map;
  int state, map_state;

  API_SETUP_ARGS(G, self, args, "Ossii", &self, &slice, &map, &state, &map_state);
  API_ASSERT(APIEnterNotModal(G));
  auto result = ExecutiveSliceNew(G, slice, map, state, map_state);
  APIExit(G);
  return APIResult(G, result);
}

static PyObject *CmdIsosurface(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  char *surf_name, *map_name, *sele;
  float lvl, fbuf;
  int surf_mode;
  int state = -1;
  float carve;
  int map_state = 0;
  int side;
  int quiet;

  API_SETUP_ARGS(G, self, args, "Osssffiifiii", &self, &surf_name, &map_name,
                        &sele, &fbuf, &lvl, &surf_mode, &state, &carve, &map_state,
                        &side, &quiet);
  API_ASSERT(APIEnterNotModal(G));
  auto res = ExecutiveIsosurfaceEtc(G, surf_name, map_name, lvl, sele, fbuf,
      state, carve, map_state, side, quiet, surf_mode);
  APIExit(G);
  return APIResult(G, res);
}

static PyObject *CmdSymExp(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  char *str1, *str2, *str3;
  float cutoff;
  pymol::CObject *mObj;
  int segi;
  int quiet;
  /* oper 0 = all, 1 = sele + buffer, 2 = vector */

  int ok = false;
  ok =
    PyArg_ParseTuple(args, "Osssfii", &self, &str1, &str2, &str3, &cutoff, &segi, &quiet);
  if(ok) {
    API_SETUP_PYMOL_GLOBALS;
    ok = (G != NULL);
  } else {
    API_HANDLE_ERROR;
  }
  if(ok && (ok = APIEnterNotModal(G))) {
    mObj = ExecutiveFindObjectByName(G, str2);
    if(mObj) {
      if(mObj->type != cObjectMolecule) {
        mObj = NULL;
        ok = false;
      }
    }
    if(mObj) {
      ExecutiveSymExp(G, str1, str2, str3, cutoff, segi, quiet);        /* TODO STATUS */
    }
    APIExit(G);
  }
  return APIResultOk(ok);
}

static PyObject *CmdSymmetryCopy(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  char *source_name, *target_name;
  int source_state, target_state;
  int quiet;
  API_SETUP_ARGS(G, self, args, "Ossiii", &self,
                        &source_name, &target_name,
                        &source_state, &target_state, &quiet);
  API_ASSERT(APIEnterNotModal(G));
  auto res = ExecutiveSymmetryCopy(G,
			  source_name, target_name,
			  source_state, target_state, quiet);
    APIExit(G);
  return APIResult(G, res);
}

static PyObject *CmdOverlap(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  char *str1, *str2;
  int state1, state2;
  float overlap = -1.0;
  float adjust;
  int ok = false;
  ok = PyArg_ParseTuple(args, "Ossiif", &self, &str1, &str2, &state1, &state2, &adjust);
  if(ok) {
    API_SETUP_PYMOL_GLOBALS;
    ok = (G != NULL);
  } else {
    API_HANDLE_ERROR;
  }
  if(ok && (ok = APIEnterNotModal(G))) {
    overlap = ExecutiveOverlap(G, str1, state1, str2, state2, adjust);
    APIExit(G);
  }
  if (overlap < 0.f) {
    return APIFailure(G);
  }
  return (Py_BuildValue("f", overlap));
}

static PyObject *CmdDist(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  char *name, *str1, *str2;
  float cutoff;
  int labels, quiet;
  int mode, reset, state, zoom;
  int state1, state2;
  API_SETUP_ARGS(G, self, args, "Osssifiiiiiii", &self, &name, &str1, &str2,
      &mode, &cutoff, &labels, &quiet, &reset, &state, &zoom, &state1, &state2);
  API_ASSERT(APIEnterNotModal(G));
  auto res = ExecutiveDistance(G, name, str1,
      str2, mode, cutoff, labels, quiet, reset, state, zoom, state1, state2);
  APIExit(G);
  return APIResult(G, res);
}

static PyObject *CmdAngle(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  char *name, *str1, *str2, *str3;
  int labels, quiet;
  int mode;
  int reset, zoom;
  int state;
  int state1, state2, state3;
  API_SETUP_ARGS(G, self, args, "Ossssiiiiiiiii", &self, &name, &str1, &str2,
      &str3, &mode, &labels, &reset, &zoom, &quiet, &state, &state1, &state2,
      &state3);
  API_ASSERT(APIEnterNotModal(G));
  auto res =
      ExecutiveAngle(G, name, str1, str2, str3,
          mode, labels, reset, zoom, quiet, state, state1, state2, state3);
  APIExit(G);
  return APIResult(G, res);
}

static PyObject *CmdDihedral(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  char *name, *str1, *str2, *str3, *str4;
  int labels, quiet;
  int mode;
  int reset, zoom;
  int state;
  API_SETUP_ARGS(G, self, args, "Osssssiiiiii", &self, &name, &str1, &str2,
      &str3, &str4, &mode, &labels, &reset, &zoom, &quiet, &state);
  API_ASSERT(APIEnterNotModal(G));
  auto res = ExecutiveDihedral(G, name,
      str1, str2, str3, str4, mode, labels, reset, zoom, quiet, state);
  APIExit(G);
  return APIResult(G, res);
}

static PyObject *CmdBond(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  char *str1, *str2;
  int order, mode;
  int quiet;
  const char* symop = "";
  API_SETUP_ARGS(G, self, args, "Ossiii|s", &self, &str1, &str2, &order, &mode, &quiet, &symop);
  API_ASSERT(APIEnterNotModal(G));
  auto result = ExecutiveBond(G, str1, str2, order, mode, quiet, symop);
  APIExit(G);
  return APIResult(G, result);
}

static PyObject* CmdAddBond(PyObject* self, PyObject* args)
{
  PyMOLGlobals* G = nullptr;
  const char* oname;
  int atm1, atm2;
  int order;
  API_SETUP_ARGS(G, self, args, "Osiii", &self, &oname, &atm1, &atm2, &order);
  APIEnterBlocked(G);

  auto result = ExecutiveAddBondByIndices(G, oname, atm1, atm2, order);

  APIExitBlocked(G);
  return APIResult(G, result);
}

static PyObject* CmdRebond(PyObject* self, PyObject* args)
{
  PyMOLGlobals* G = nullptr;
  const char* oname;
  int state;
  int pbc = 0;
  API_SETUP_ARGS(G, self, args, "Osi|i", &self, &oname, &state, &pbc);
  API_ASSERT(APIEnterNotModal(G));
  auto res = ExecutiveRebond(G, oname, state, pbc);
  APIExit(G);
  return APIResult(G, res);
}

static PyObject *CmdRevalence(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  char *sele1, *sele2, *source;
  int source_state, target_state, reset;
  int quiet;
  API_SETUP_ARGS(G, self, args, "Osssiiii", &self, &sele1, &sele2, &source,
                        &target_state, &source_state, &reset, &quiet);
  API_ASSERT(APIEnterNotModal(G));
  auto result = ExecutiveRevalence(
      G, sele1, sele2, source, target_state, source_state, reset, quiet);
  APIExit(G);
  return APIResult(G, result);
}

static PyObject* CmdPBCUnwrap(PyObject* self, PyObject* args)
{
  PyMOLGlobals* G = nullptr;
  const char* oname;
  int bymol = true;
  API_SETUP_ARGS(G, self, args, "Osi", &self, &oname, &bymol);
  API_ASSERT(APIEnterNotModal(G));

  auto obj = ExecutiveFindObjectMoleculeByName(G, oname);
  if (!obj) {
    APIExit(G);
    return APIFailure(G, "cannot find object");
  }

  ObjectMoleculePBCUnwrap(*obj, bymol);

  APIExit(G);
  return APISuccess();
}

static PyObject* CmdPBCWrap(PyObject* self, PyObject* args)
{
  PyMOLGlobals* G = nullptr;
  const char* oname;
  PyObject* pycenter = nullptr;
  API_SETUP_ARGS(G, self, args, "OsO", &self, &oname, &pycenter);

  std::vector<float> center;
  if (pycenter != Py_None) {
    API_ASSERT(PConvFromPyObject(G, pycenter, center) && center.size() == 3);
  }

  API_ASSERT(APIEnterNotModal(G));

  auto obj = ExecutiveFindObjectMoleculeByName(G, oname);
  if (!obj) {
    APIExit(G);
    return APIFailure(G, "cannot find object");
  }

  ObjectMoleculePBCWrap(*obj, center.empty() ? nullptr : center.data());

  APIExit(G);
  return APISuccess();
}

static PyObject *CmdVdwFit(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  char *str1, *str2;
  int state1, state2, quiet;
  float buffer;
  int ok = false;
  ok =
    PyArg_ParseTuple(args, "Osisifi", &self, &str1, &state1, &str2, &state2, &buffer,
                     &quiet);
  if(ok) {
    API_SETUP_PYMOL_GLOBALS;
    ok = (G != NULL);
  } else {
    API_HANDLE_ERROR;
  }
  if(ok && (ok = APIEnterNotModal(G))) {
    ok = ExecutiveVdwFit(G, str1, state1, str2, state2, buffer, quiet);
    APIExit(G);
  }
  return APIResultOk(ok);
}

static PyObject *CmdLabel(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  char *str1, *str2;
  int quiet;
  API_SETUP_ARGS(G, self, args, "Ossi", &self, &str1, &str2, &quiet);
  API_ASSERT(APIEnterBlockedNotModal(G));
  ExecutiveLabel(G, str1, str2, quiet, cExecutiveLabelEvalOn);
  APIExitBlocked(G);
  if (PyErr_Occurred()) {
    return nullptr;
  }
  return APISuccess();
}

static PyObject *CmdLabel2(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  char *str1, *str2;
  int quiet;
  API_SETUP_ARGS(G, self, args, "Ossi", &self, &str1, &str2, &quiet);
  API_ASSERT(APIEnterNotModal(G));
  ExecutiveLabel(G, str1, str2, quiet, cExecutiveLabelEvalAlt);
  APIExit(G);
  return APISuccess();
}

static PyObject *CmdAlter(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  char *str1, *str2;
  int read_only, quiet;
  PyObject *space;
  API_SETUP_ARGS(G, self, args, "OssiiO", &self, &str1, &str2, &read_only, &quiet, &space);
  API_ASSERT(APIEnterBlockedNotModal(G));
  pymol::Result<int> result{-1};
  if (read_only) {
    result = ExecutiveIterate(
        G, str1, str2, read_only, quiet, space);
  } else {
    result = ExecutiveIterate(G, str1, str2,
        read_only, quiet, space);
  }
  APIExitBlocked(G);
  return APIResult(G, result);
}

static PyObject *CmdAlterList(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  char *str1;
  int quiet;
  PyObject *space;
  PyObject *list;
  API_SETUP_ARGS(G, self, args, "OsOiO", &self, &str1, &list, &quiet, &space);
  API_ASSERT(APIEnterBlockedNotModal(G));
  auto result = ExecutiveIterateList(G, str1, list, false, quiet, space);
  APIExitBlocked(G);
  return APIResult(G, result);
}

static PyObject *CmdSelectList(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  char *str1, *sele_name;
  int quiet;
  int mode;
  int state;
  PyObject *list;
  API_SETUP_ARGS(G, self, args, "OssO!iii", &self, &sele_name, &str1,
      &PyList_Type, &list, &state, &mode, &quiet);
  std::vector<int> int_array;
  API_ASSERT(PConvFromPyObject(G, list, int_array));
  API_ASSERT(APIEnterNotModal(G));
  auto result = ExecutiveSelectList(G, sele_name, str1, int_array.data(),
      int_array.size(), state, mode, quiet);
  SceneInvalidate(G);
  SeqDirty(G);
  APIExit(G);
  return APIResult(G, result);
}

static PyObject *CmdAlterState(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  char *str1, *str2;
  int state, read_only, quiet;
  PyObject *obj;
  API_SETUP_ARGS(G, self, args, "OissiiO", &self, &state, &str1, &str2, &read_only, &quiet, &obj);
  API_ASSERT(APIEnterBlockedNotModal(G));
  auto result = ExecutiveIterateState(G, state, str1, str2, read_only,  quiet, obj);
  APIExitBlocked(G);
  return APIResult(G, result);
}

static PyObject *CmdCopy(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  char *str1, *str2;
  int zoom;
  API_SETUP_ARGS(G, self, args, "Ossi", &self, &str1, &str2, &zoom);
  API_ASSERT(APIEnterNotModal(G));
  auto result = ExecutiveCopy(G, str1, str2, zoom);
  APIExit(G);
  return APIResult(G, result);
}

static PyObject *CmdRecolor(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  char *str1;
  cRep_t rep;
  API_SETUP_ARGS(G, self, args, "Osi", &self, &str1, &rep);
  API_ASSERT(APIEnterNotModal(G));
  auto result = ExecutiveInvalidateRep(G, str1, rep, cRepInvColor);
  APIExit(G);
  return APIResult(G, result);
}

static PyObject *CmdRebuild(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  const char* str1;
  cRep_t rep;
  API_SETUP_ARGS(G, self, args, "Osi", &self, &str1, &rep);
  API_ASSERT(APIEnterNotModal(G));

  pymol::Result<> res;
  if (WordMatchExact(G, str1, cKeywordAll, true)) {
    ExecutiveRebuildAll(G);
  } else {
    const cRepInv_t level = SettingGet<bool>(G, cSetting_defer_builds_mode)
                                ? cRepInvPurgeAll
                                : cRepInvAll;
    res = ExecutiveInvalidateRep(G, str1, rep, level);
  }

  APIExit(G);
  return APIResult(G, res);
}

static PyObject *CmdResetRate(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  API_SETUP_ARGS(G, self, args, "O", &self);
  API_ASSERT(APIEnterNotModal(G));
  ButModeResetRate(G);
  APIExit(G);
  return APISuccess();
}

static PyObject *CmdReady(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  int ok = false;
  ok = PyArg_ParseTuple(args, "O", &self);
  if(ok) {
    API_SETUP_PYMOL_GLOBALS;
    ok = (G != NULL);
  } else {
    API_HANDLE_ERROR;
  }
  if(ok) {
    return (APIResultCode(G->Ready));
  } else {
    return (APIResultCode(0));
  }
}

static PyObject *CmdMem(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  int ok = false;
  ok = PyArg_ParseTuple(args, "O", &self);
  if(ok) {
    API_SETUP_PYMOL_GLOBALS;
    ok = (G != NULL);
  } else {
    API_HANDLE_ERROR;
  }
  if(ok) {
    OVHeap_Dump(G->Context->heap, 0);
    SelectorMemoryDump(G);
    ExecutiveMemoryDump(G);
  }
  return APISuccess();
}

/**
 * Simple `glViewport` wrapper to call from Python without depending on
 * the heavy PyOpenGL module.
 */
static PyObject *Cmd_glViewport(PyObject * self, PyObject * args)
{
  int x, y, width, height;

  if(!PyArg_ParseTuple(args, "iiii", &x, &y, &width, &height)) {
    API_HANDLE_ERROR;
  } else {
    glViewport(x, y, width, height);
  }

  return APIIncRef(Py_None);
}

static void PyMOLGlobalsCapsuleDestructor(PyObject* self)
{
  assert(self != Py_None);
  auto G = _api_get_pymol_globals(self);
  assert(G);
  PyMOL_Free(G->PyMOL);
}

/**
 * Create a `_COb` instance.
 *
 * @param pymol The `pymol` module or a `pymol2.PyMOL` instance
 * @param options Options like `pymol.invocation.options`
 * @param singleton (bool) Whether this instance should be registered as
 * the SingletonPyMOLGlobals (raises RuntimeError if it already exists)
 */
static PyObject *Cmd_New(PyObject * self, PyObject * args)
{
  PyObject *pymol = NULL;       /* pymol object instance */
  PyObject *pyoptions = Py_None;
  int singleton = false;

  if (!PyArg_ParseTuple(args, "O|Op", &pymol, &pyoptions, &singleton)) {
    return nullptr;
  }

  if (singleton && SingletonPyMOLGlobals) {
    PyErr_SetString(PyExc_RuntimeError, "Singleton already exists");
    return nullptr;
  }

  CPyMOLOptions* options = PyMOLOptions_New();
  assert(options);

  if (pyoptions != Py_None) {
    PConvertOptions(options, pyoptions);
  }

  CPyMOL* I = PyMOL_NewWithOptions(options);
  PyMOLOptions_Free(options);

  if (!I) {
    PyErr_SetString(PyExc_Exception, "PyMOL_NewWithOptions failed");
    return nullptr;
  }

  PyMOLGlobals* G = PyMOL_GetGlobals(I);
  assert(G);

  if (singleton) {
    assert(!SingletonPyMOLGlobals);
    SingletonPyMOLGlobals = G;
  } else {
    // Creating a non-singleton instance disables auto-library mode
    auto_library_mode_disabled = true;
  }

  G->P_inst = pymol::calloc<CP_inst>(1);
  G->P_inst->obj = pymol;
  G->P_inst->dict = PyObject_GetAttrString(pymol, "__dict__");
  Py_DECREF(G->P_inst->dict); // borrow reference

  PyObject* tmp = PyCapsule_New(I, nullptr, nullptr);
  PyObject_SetAttrString(pymol, "__pymol__", tmp);
  Py_DECREF(tmp);

  for (int a = 0; a < MAX_SAVED_THREAD; ++a) {
    G->P_inst->savedThread[a].id = -1;
  }

  return PyCapsule_New(
      PyMOL_GetGlobalsHandle(I), nullptr, PyMOLGlobalsCapsuleDestructor);
}

static PyObject *Cmd_Start(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  PyObject *cmd = NULL;
  int ok = true;
  ok = PyArg_ParseTuple(args, "OO", &self, &cmd);
  if(ok) {
    API_SETUP_PYMOL_GLOBALS;
    ok = (G && G->PyMOL);
  } else {
    API_HANDLE_ERROR;
  }
  if(ok) {
    G->P_inst->cmd = cmd;
    PyMOL_StartWithPython(G->PyMOL);
  }
  return APIResultOk(ok);
}

static PyObject *Cmd_Stop(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  int ok = true;
  ok = PyArg_ParseTuple(args, "O", &self);
  if(ok) {
    API_SETUP_PYMOL_GLOBALS;
    ok = (G && G->PyMOL);
  } else {
    API_HANDLE_ERROR;
  }
  if(ok) {
    PyMOL_Stop(G->PyMOL);
  }
  return APIResultOk(ok);
}

static PyObject *Cmd_Idle(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  int ok = true;
  int result = false;
  ok = PyArg_ParseTuple(args, "O", &self);
  if(ok) {
    API_SETUP_PYMOL_GLOBALS;
    ok = (G && G->PyMOL);
  } else {
    API_HANDLE_ERROR;
  }
  if(ok && PTryLockAPIAndUnblock(G)) {
    result = PyMOL_Idle(G->PyMOL);
    PBlockAndUnlockAPI(G);
  }
  return APIResultCode(result);
}

static PyObject *Cmd_Reshape(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  int ok = true;
  int width, height, force;
  ok = PyArg_ParseTuple(args, "Oiii", &self, &width, &height, &force);
  if(ok) {
    API_SETUP_PYMOL_GLOBALS;
    ok = (G && G->PyMOL);
  } else {
    API_HANDLE_ERROR;
  }
  if(ok && PTryLockAPIAndUnblock(G)) {
    PyMOL_Reshape(G->PyMOL, width, height, force);
    PBlockAndUnlockAPI(G);
  }
  return APIResultOk(ok);
}

static PyObject *Cmd_GetRedisplay(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  int ok = true;
  int result = false;
  int reset;
  ok = PyArg_ParseTuple(args, "Oi", &self, &reset);
  if(ok) {
    API_SETUP_PYMOL_GLOBALS;
    ok = (G && G->PyMOL);
  } else {
    API_HANDLE_ERROR;
  }
  if(ok && PTryLockAPIAndUnblock(G)) {
    result = PyMOL_GetRedisplay(G->PyMOL, reset);
    PBlockAndUnlockAPI(G);
  }
  return APIResultCode(result);
}

static PyObject *Cmd_Draw(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  int ok = true;
  ok = PyArg_ParseTuple(args, "O", &self);
  if(ok) {
    API_SETUP_PYMOL_GLOBALS;
    ok = (G && G->PyMOL);
  } else {
    API_HANDLE_ERROR;
  }
  if(ok && PTryLockAPIAndUnblock(G)) {
    PyMOL_Draw(G->PyMOL);
    PBlockAndUnlockAPI(G);
  }
  return APIResultOk(ok);
}

static PyObject *Cmd_Button(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  int ok = true;
  int button, state;
  int x, y, modifiers;
  ok = PyArg_ParseTuple(args, "Oiiiii", &self, &button, &state, &x, &y, &modifiers);
  if(ok) {
    API_SETUP_PYMOL_GLOBALS;
    ok = (G && G->PyMOL);
  } else {
    API_HANDLE_ERROR;
  }
  if(ok && PTryLockAPIAndUnblock(G)) {
    PyMOL_Button(G->PyMOL, button, state, x, y, modifiers);
    PBlockAndUnlockAPI(G);
  }
  return APIResultOk(ok);
}

static PyObject *Cmd_Drag(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  int ok = true;
  int x, y, modifiers;
  ok = PyArg_ParseTuple(args, "Oiii", &self, &x, &y, &modifiers);
  if(ok) {
    API_SETUP_PYMOL_GLOBALS;
    ok = (G && G->PyMOL);
  } else {
    API_HANDLE_ERROR;
  }
  if(ok && PTryLockAPIAndUnblock(G)) {
    PyMOL_Drag(G->PyMOL, x, y, modifiers);
    PBlockAndUnlockAPI(G);
  }
  return APIResultOk(ok);
}

static PyObject *Cmd_Sdof(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  float tx, ty, tz, rx, ry, rz;
  int ok = false;
  ok = PyArg_ParseTuple(args, "Offffff", &self, &tx, &ty, &tz, &rx, &ry, &rz);
  if(ok) {
    API_SETUP_PYMOL_GLOBALS;
    ok = (G != NULL);
  } else {
    API_HANDLE_ERROR;
  }
  if(ok) {
    PLockStatus(G);
    ControlSdofUpdate(G, tx, ty, tz, rx, ry, rz);
    PUnlockStatus(G);
  }
  return APIResultOk(ok);
}

static PyObject *CmdRunPyMOL(PyObject * self, PyObject * args)
{
#ifdef _PYMOL_NO_MAIN
  PyErr_SetString(PyExc_NotImplementedError, "compile with --glut");
  return nullptr;
#else

  if(run_only_once) {
    run_only_once = false;

    // _PYMOL_MODULE stuff
    {
      int block_input_hook = false;
      if(!PyArg_ParseTuple(args, "Oi", &self, &block_input_hook))
        block_input_hook = false;

      main_shared(block_input_hook);
    }
  }
#endif

  return APISuccess();
}

static PyObject *CmdCountStates(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  char *str1;
  OrthoLineType s1;
  int ok = false;
  int count = 0;
  API_SETUP_ARGS(G, self, args, "Os", &self, &str1);
  APIEnter(G);
    ok = (SelectorGetTmp2(G, str1, s1) >= 0);
    count = ExecutiveCountStates(G, s1);
    if(count < 0)
      ok = false;
    SelectorFreeTmp(G, s1);
    APIExit(G);
  return ok ? APIResultCode(count) : APIFailure(G);
}

static PyObject *CmdCountFrames(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  API_SETUP_ARGS(G, self, args, "O", &self);
  APIEnter(G);
  SceneCountFrames(G);
  int result = SceneGetNFrame(G, NULL);
  APIExit(G);
  return (APIResultCode(result));
}

static PyObject *CmdGetMovieLength(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  int result = 0;
  API_SETUP_ARGS(G, self, args, "O", &self);
  APIEnter(G);
    result = MovieGetLength(G);
    APIExit(G);
  return (APIResultCode(result));
}

static PyObject *CmdIdentify(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  char *str1;
  int mode;
  PyObject *result = nullptr;
  pymol::vla<int> iVLA;
  pymol::vla<ObjectMolecule*> oVLA;

  API_SETUP_ARGS(G, self, args, "Osi", &self, &str1, &mode);
  APIEnter(G);
  int l = ExecutiveIdentifyObjects(G, str1, mode, &iVLA, &oVLA);
  APIExit(G);

  if (l < 0) {
    return APIFailure(G, "invalid selection");
  }

  if(!iVLA) {
    result = PyList_New(0);
  } else if(!mode) {
    result = PConvIntVLAToPyList(iVLA);
  } else {                  /* object mode */
    result = PyList_New(l);
    for (int a = 0; a < l; ++a) {
      auto tuple = Py_BuildValue("si", oVLA[a]->Name, iVLA[a]);
      PyList_SetItem(result, a, tuple);
    }
  }

  return result;
}

static PyObject *CmdIndex(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  char *str1;
  int mode;
  pymol::vla<ObjectMolecule*> oVLA;
  pymol::vla<int> iVLA;

  API_SETUP_ARGS(G, self, args, "Osi", &self, &str1, &mode);
  APIEnter(G);
  int l = ExecutiveIndex(G, str1, mode, &iVLA, &oVLA);
  APIExit(G);

  if (l == -1) {
    return APIFailure(G, "invalid selection");
  }

  PyObject* result = PyList_New(l);
  for (int a = 0; a < l; ++a) {
    auto tuple = Py_BuildValue("si", oVLA[a]->Name, iVLA[a] + 1 /* 1-based */);
    PyList_SetItem(result, a, tuple);
  }
  return result;
}

static PyObject *CmdFindPairs(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  char *str1, *str2;
  int state1, state2;
  float cutoff;
  float angle;
  int mode;
  pymol::vla<int> iVLA;
  pymol::vla<ObjectMolecule*> oVLA;

  API_SETUP_ARGS(G, self, args, "Ossiiiff", &self, &str1, &str2, &state1, &state2, &mode,
                     &cutoff, &angle);
  APIEnter(G);
  auto res = ExecutivePairIndices(
      G, str1, str2, state1, state2, mode, cutoff, angle, &iVLA, &oVLA);
  APIExit(G);

  if (!res) {
    return APIFailure(G, res.error());
  }

  int l = res.result();
  PyObject* result = PyList_New(l);
  for (int a = 0; a < l; ++a) {
    auto tuple = Py_BuildValue("(si)(si)", //
        oVLA[a * 2]->Name, iVLA[a * 2] + 1 /* 1-based */, oVLA[a * 2 + 1]->Name,
        iVLA[a * 2 + 1] + 1 /* 1-based */);
    PyList_SetItem(result, a, tuple);
  }
  return result;
}

static PyObject *CmdSystem(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  char *str1;
  int ok = false;
  int async;
  ok = PyArg_ParseTuple(args, "Osi", &self, &str1, &async);
  if(ok) {
    API_SETUP_PYMOL_GLOBALS;
    ok = (G != NULL);
  } else {
    API_HANDLE_ERROR;
  }
  if(ok) {
    if(async) {
      PUnblock(G);              /* free up PyMOL and the API */
      ok = system(str1);
      PBlock(G);
    } else if((ok = APIEnterNotModal(G))) {     /* keep PyMOL locked */
      ok = system(str1);
      APIExit(G);
    }
  }
  return APIResultOk(ok);
}

static PyObject *CmdGetFeedback(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  int ok = false;
  ok = PyArg_ParseTuple(args, "O", &self);
  if(ok) {
    API_SETUP_PYMOL_GLOBALS;
    ok = (G != NULL);
  } else {
    API_HANDLE_ERROR;
  }
  if(ok) {
    if(G->Ready) {
      PyObject *result = NULL;

      if(G->Terminating) {      /* try to bail */
        /* BEGIN PROPRIETARY CODE SEGMENT (see disclaimer in "os_proprietary.h") */
#ifdef WIN32
        abort();
#endif
        /* END PROPRIETARY CODE SEGMENT */
        exit(0);
      }

      /* ALLOWED DURING MODAL DRAWING */
      APIEnterBlocked(G);
      auto buffer = OrthoFeedbackOut(G, *G->Ortho);
      APIExitBlocked(G);
      if(!buffer.empty())
        result = Py_BuildValue("s", buffer.c_str());
      return (APIAutoNone(result));
    }
  }
  return (APIAutoNone(NULL));
}

static PyObject *CmdGetSeqAlignStr(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  char *str1;
  char *seq = NULL;
  int state;
  int format;
  int quiet;
  PyObject *result = NULL;
  int ok = false;
  ok = PyArg_ParseTuple(args, "Osiii", &self, &str1, &state, &format, &quiet);
  if(ok) {
    API_SETUP_PYMOL_GLOBALS;
    ok = (G != NULL);
  } else {
    API_HANDLE_ERROR;
  }
  if(ok) {
    APIEnter(G);
    seq = ExecutiveNameToSeqAlignStrVLA(G, str1, state, format, quiet);
    APIExit(G);
    if(seq)
      result = Py_BuildValue("s", seq);
    VLAFreeP(seq);
  }
  return (APIAutoNone(result));
}

static PyObject *CmdGetStr(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  PyObject *result = NULL;
  pymol::vla<char> vla;
  char *format;
  char *sele;
  int state;
  char *ref;
  int ref_state;
  int quiet;
  int multi;

  API_SETUP_ARGS(G, self, args, "Ossisiii", &self, &format, &sele, &state, &ref,
      &ref_state, &multi, &quiet);
  APIEnter(G);
  vla = MoleculeExporterGetStr(G, format, sele, state,
      ref, ref_state, multi, quiet);
  APIExit(G);

  if (vla) {
    result = PyBytes_FromStringAndSize(vla, vla.size());
  }

  return APIAutoNone(result);
}

static PyObject *CmdGetModel(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  char *str1;
  int state;
  char *ref_object;
  int ref_state;
  OrthoLineType s1;
  PyObject *result = NULL;
  int ok = false;
  API_SETUP_ARGS(G, self, args, "Osisi", &self, &str1, &state, &ref_object, &ref_state);
    if(!ref_object[0])
      ref_object = NULL;
  APIEnterBlocked(G);
      ok = (SelectorGetTmp(G, str1, s1) >= 0);
      if(ok)
        result = ExecutiveSeleToChemPyModel(G, s1, state, ref_object, ref_state);
      SelectorFreeTmp(G, s1);
      APIExitBlocked(G);
  if (!result)
    return APIFailure(G);
  return result;
}

static PyObject *CmdGetBonds(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  PyObject *result = NULL;
  char *sele;
  int state;

  API_SETUP_ARGS(G, self, args, "Osi", &self, &sele, &state);
  APIEnter(G);

  result = MoleculeExporterGetPyBonds(G, sele, state);
  APIExit(G);

  return APIAutoNone(result);
}

static PyObject *CmdCreate(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  char *str1, *str2;
  int target, source, discrete, quiet;
  int singletons;
  int copy_properties = 0;
  int zoom;
  API_SETUP_ARGS(G, self, args, "Ossiiiiiii", &self, &str1, &str2, &source,
                        &target, &discrete, &zoom, &quiet, &singletons, &copy_properties);
  API_ASSERT(APIEnterNotModal(G));
  auto res = ExecutiveSeleToObject(G, str1,
      str2, source, target, discrete, zoom, quiet, singletons, copy_properties);
  APIExit(G);
  return APIResult(G, res);
}

static PyObject *CmdOrient(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  char *str1;
  int state;
  float animate;
  int quiet = false;            /* TODO */
  API_SETUP_ARGS(G, self, args, "Osif", &self, &str1, &state, &animate);
  API_ASSERT(APIEnterNotModal(G));
  auto result = ExecutiveOrient(G, str1, state, animate, false, 0.0F, quiet);
  APIExit(G);
  return APIResult(G, result);
}

static PyObject *CmdFitPairs(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  PyObject *list;
  int quiet = 0;
  API_SETUP_ARGS(G, self, args, "OOi", &self, &list, &quiet);
  API_ASSERT(APIEnterBlockedNotModal(G));
  auto result = ExecutiveFitPairs(G, list, quiet);
  APIExitBlocked(G);
  return APIResult(G, result);
}

static PyObject *CmdIntraFit(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  char *str1;
  int state;
  int mode;
  int quiet;
  int mix;
  int pbc = true;
  API_SETUP_ARGS(G, self, args, "Osiiii|i", &self, &str1, &state, &mode, &quiet, &mix, &pbc);
  API_ASSERT(APIEnterNotModal(G));
  auto fVLA = ExecutiveRMSStates(G, str1, state, mode, quiet, mix, pbc);
  APIExit(G);
  PyObject* result = nullptr;
  if(fVLA) {
    result = PConvFloatVLAToPyList(fVLA.result().data());
  }
  return APIAutoNone(result);
}

static PyObject *CmdGetAtomCoords(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  char *str1;
  int state;
  int quiet;
  API_SETUP_ARGS(G, self, args, "Osii", &self, &str1, &state, &quiet);
  APIEnter(G);
  auto result = ExecutiveGetAtomVertex(G, str1, state, quiet);
  APIExit(G);
  return APIResult(G, result);
}

static PyObject *CmdFit(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  char *str1, *str2;
  int mode;
  int quiet;
  float cutoff;
  int state1, state2;
  int matchmaker, cycles;
  char *object;
  API_SETUP_ARGS(G, self, args, "Ossiiiiifis", &self, &str1, &str2, &mode,
      &state1, &state2, &quiet, &matchmaker, &cutoff, &cycles, &object);
  API_ASSERT(APIEnterNotModal(G));
  auto result = ExecutiveFit(G, str1, str2, mode, cutoff, cycles, quiet, object, state1, state2, matchmaker);
  APIExit(G);
  if(result) {
    return PConvToPyObject(result.result().final_rms);
  } else {
    return APIFailure(G, result.error());
  }
}

static PyObject *CmdUpdate(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  char *str1, *str2;
  int int1, int2;
  int matchmaker, quiet;
  API_SETUP_ARGS(G, self, args, "Ossiiii", &self, &str1, &str2, &int1, &int2,
      &matchmaker, &quiet);
  API_ASSERT(APIEnterNotModal(G));
  auto result = ExecutiveUpdateCmd(G, str1, str2, int1, int2, matchmaker, quiet); 
  APIExit(G);
  return APIResult(G, result);
}

static PyObject *CmdDirty(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  API_SETUP_ARGS(G, self, args, "O", &self);
  API_ASSERT(APIEnterNotModal(G));
  OrthoDirty(G);
  APIExit(G);
  return APISuccess();
}

static PyObject *CmdGetObjectList(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  char *str1;
  OrthoLineType s1;
  int ok = false;
  PyObject *result = NULL;

  ok = PyArg_ParseTuple(args, "Os", &self, &str1);

  if(ok) {
    API_SETUP_PYMOL_GLOBALS;
    ok = (G != NULL);
  } else {
    API_HANDLE_ERROR;
  }
  if(ok && (ok = APIEnterBlockedNotModal(G))) {
    ok = (SelectorGetTmp(G, str1, s1) >= 0);
    auto list = ExecutiveGetObjectMoleculeVLA(G, s1);
    if(list) {
      unsigned int size = VLAGetSize(list);
      result = PyList_New(size);
      if(result) {
        unsigned int a;
        for(a = 0; a < size; a++) {
          PyList_SetItem(result, a, PyString_FromString(list[a]->Name));
        }
      }
    }
    SelectorFreeTmp(G, s1);
    APIExitBlocked(G);
  }
  return (APIAutoNone(result));
}

static PyObject *CmdGetDistance(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  char *str1, *str2;
  int int1;
  API_SETUP_ARGS(G, self, args, "Ossi", &self, &str1, &str2, &int1);
  APIEnter(G);
  auto res = ExecutiveGetDistance(G, str1, str2, int1);
  APIExit(G);
  return APIResult(G, res);
}

static PyObject *CmdGetAngle(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  char *str1, *str2, *str3;
  int int1;
  API_SETUP_ARGS(G, self, args, "Osssi", &self, &str1, &str2, &str3, &int1);
  APIEnter(G);
  auto res = ExecutiveGetAngle(G, str1, str2, str3, int1);
  APIExit(G);
  return APIResult(G, res);
}

static PyObject *CmdGetDihe(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  char *str1, *str2, *str3, *str4;
  int int1;
  API_SETUP_ARGS(G, self, args, "Ossssi", &self, &str1, &str2, &str3, &str4, &int1);
  APIEnter(G);
  auto res = ExecutiveGetDihe(G, str1, str2, str3, str4, int1);
  APIExit(G);
  return APIResult(G, res);
}

static PyObject *CmdSetDihe(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  char *str1, *str2, *str3, *str4;
  float float1;
  int int1;
  int quiet;
  API_SETUP_ARGS(G, self, args, "Ossssfii", &self, &str1, &str2, &str3, &str4,
      &float1, &int1, &quiet);
  API_ASSERT(APIEnterNotModal(G));
  auto result = ExecutiveSetDihe(G, str1, str2, str3, str4, float1, int1, quiet);
  APIExit(G);
  return APIResult(G, result);
}

static PyObject *CmdDo(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  char *str1;
  int log;
  int echo;
  API_SETUP_ARGS(G, self, args, "Osii", &self, &str1, &log, &echo);
  API_ASSERT(APIEnterNotModal(G));

    if(str1[0] != '_') {        /* suppress internal call-backs */
      if(strncmp(str1, "cmd._", 5) && (strncmp(str1, "_cmd.", 5))) {
        if(echo) {
          OrthoAddOutput(G, "PyMOL>");
          OrthoAddOutput(G, str1);
          OrthoNewLine(G, NULL, true);
        }
        if((str1[0] == 'P') && (str1[1] == 'y') && (str1[2] == 'M') &&
           (str1[3] == 'O') && (str1[4] == 'L') && (str1[5] == '>')) {
          /* ignore pasted-in /PyMOL>\s?/ for sake of end-user convenience */
          str1 += 6;
          if(str1[0] == ' ')
            str1++;
        }
        if(log)
          if(WordMatch(G, str1, "quit", true) == 0)     /* don't log quit */
            PLog(G, str1, cPLog_pml);
      }
      PParse(G, str1);
    } else if(str1[1] == ' ') {
      /* "_ command" suppresses echoing of command, but it is still logged */
      if(log)
        if(WordMatch(G, str1 + 2, "quit", true) == 0)   /* don't log quit */
          PLog(G, str1 + 2, cPLog_pml);
      PParse(G, str1 + 2);
    } else {
      PParse(G, str1);
    }
    APIExit(G);
  return APISuccess();
}

static PyObject *CmdRock(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  int int1;
  API_SETUP_ARGS(G, self, args, "Oi", &self, &int1);
  API_ASSERT(APIEnterNotModal(G));
  auto result = ControlRock(G, int1);
  APIExit(G);
  return APIResult(G, result);
}

static PyObject *CmdBusyDraw(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  int int1;
  int ok = false;
  ok = PyArg_ParseTuple(args, "Oi", &self, &int1);
  if(ok) {
    API_SETUP_PYMOL_GLOBALS;
    ok = (G != NULL);
  } else {
    API_HANDLE_ERROR;
  }
  if(ok && (ok = APIEnterNotModal(G))) {
    if(SettingGetGlobal_b(G, cSetting_show_progress)) {
      OrthoBusyDraw(G, int1);
    }
    APIExit(G);
  }
  return APIResultOk(ok);
}

static PyObject *CmdSetBusy(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  int int1;
  int ok = false;
  ok = PyArg_ParseTuple(args, "Oi", &self, &int1);
  if(ok) {
    API_SETUP_PYMOL_GLOBALS;
    ok = (G != NULL);
  } else {
    API_HANDLE_ERROR;
  }
  if(ok) {
    PLockStatus(G);
    PyMOL_SetBusy(G->PyMOL, int1);
    PUnlockStatus(G);
  }
  return APIResultOk(ok);
}

static PyObject *CmdGetBusy(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  int result = 0;
  int ok = false;
  int int1;
  ok = PyArg_ParseTuple(args, "Oi", &self, &int1);
  if(ok) {
    API_SETUP_PYMOL_GLOBALS;
    ok = (G != NULL);
  } else {
    API_HANDLE_ERROR;
  }
  if(ok) {
    PLockStatus(G);
    result = PyMOL_GetBusy(G->PyMOL, int1);
    PUnlockStatus(G);
  }
  return (APIResultCode(result));
}

static PyObject *CmdGetProgress(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  int ok = false;
  int reset;                    /* TO DO */
  ok = PyArg_ParseTuple(args, "Oi", &self, &reset);
  if(ok) {
    API_SETUP_PYMOL_GLOBALS;
    ok = (G != NULL);
  } else {
    API_HANDLE_ERROR;
  }
  if(ok) {
    if(G->Ready && (!SettingGetGlobal_b(G, cSetting_sculpting))) {

      /* assumes status is already locked */

      float result = -1.0F;
      float value = 0.0F, range = 1.0F;
      int ok = false;
      int int1;
      int offset;
      int progress[PYMOL_PROGRESS_SIZE];

      ok = PyArg_ParseTuple(args, "Oi", &self, &int1);
      if(ok) {
        if(PyMOL_GetBusy(G->PyMOL, false)) {
          PyMOL_GetProgress(G->PyMOL, progress, false);

          for(offset = PYMOL_PROGRESS_FAST; offset >= PYMOL_PROGRESS_SLOW; offset -= 2) {
            if(progress[offset + 1]) {
              float old_value = value;
              float old_range = range;

              range = (float) (progress[offset + 1]);
              value = (float) (progress[offset]);

              value += (1.0F / range) * (old_value / old_range);

              result = value / range;
            }
          }
        }
      }
      return (PyFloat_FromDouble((double) result));
    }
  }
  return (PyFloat_FromDouble(-1.0));
}

static PyObject *CmdGetMoment(PyObject * self, PyObject * args)
{                               /* missing? */
  PyMOLGlobals *G = NULL;
  double moment[16];
  PyObject *result;
  char *str1;
  int ok = false;
  int state;

  ok = PyArg_ParseTuple(args, "Osi", &self, &str1, &state);
  if(ok) {
    API_SETUP_PYMOL_GLOBALS;
    ok = (G != NULL);
  } else {
    API_HANDLE_ERROR;
  }
  if(ok) {
    APIEnter(G);
    ExecutiveGetMoment(G, str1, moment, state);
    APIExit(G);
  }
  result = Py_BuildValue("(ddd)(ddd)(ddd)",
                         moment[0], moment[1], moment[2],
                         moment[3], moment[4], moment[5],
                         moment[6], moment[7], moment[8]);

  return result;
}

static PyObject *CmdGetSettingType(PyObject *, PyObject * args)
{
  int index;
  if (!PyArg_ParseTuple(args, "i", &index)) {
    return nullptr;
  }
  return PyLong_FromLong(SettingGetType(index));
}

static PyObject* CmdGetSettingLevel(PyObject*, PyObject* args)
{
  unsigned index;
  API_ASSERT(PyArg_ParseTuple(args, "I", &index));
  API_ASSERT(index < cSetting_INIT);
  return PyUnicode_FromString(SettingLevelGetName(index));
}

static PyObject *CmdGetSettingOfType(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  int int1, int2, int3;
  char *str1;
  API_SETUP_ARGS(G, self, args, "Oisii", &self, &int1, &str1, &int2, &int3);
  APIEnterBlocked(G);
  // may set the Python error indicator and return NULL
  auto result = ExecutiveGetSettingOfType(G, int1, str1, int2, int3);
  APIExitBlocked(G);
  return result;
}

static PyObject *CmdSetFrame(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  int mode, frm;
  API_SETUP_ARGS(G, self, args, "Oii", &self, &mode, &frm);
  API_ASSERT(APIEnterNotModal(G));
  SceneSetFrame(G, mode, frm);
  APIExit(G);
  return APISuccess();
}

static PyObject *CmdFrame(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  int frm,trigger;
  API_SETUP_ARGS(G, self, args, "Oii", &self, &frm, &trigger);
  API_ASSERT(APIEnterNotModal(G));
  SceneSetFrame(G, trigger ? 4 : 0, frm);
  APIExit(G);
  return APISuccess();
}

static PyObject *CmdStereo(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  int i1;
  API_SETUP_ARGS(G, self, args, "Oi", &self, &i1);
  API_ASSERT(APIEnterNotModal(G));
  auto res = ExecutiveStereo(G, i1);
  APIExit(G);
  return APIResult(G, res);
}

static PyObject *CmdReset(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  char *obj;
  API_SETUP_ARGS(G, self, args, "Os", &self, &obj);
  API_ASSERT(APIEnterNotModal(G));
  auto result = ExecutiveReset(G, obj);
  APIExit(G);
  return APIResult(G, result);
}

#if 0
static PyObject *CmdSetMatrix(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  float m[16];
  int ok = false;
  ok = PyArg_ParseTuple(args, "Offffffffffffffff", &self,
                        &m[0], &m[1], &m[2], &m[3],
                        &m[4], &m[5], &m[6], &m[7],
                        &m[8], &m[9], &m[10], &m[11], &m[12], &m[13], &m[14], &m[15]);
  if(ok) {
    API_SETUP_PYMOL_GLOBALS;
    ok = (G != NULL);
  } else {
    API_HANDLE_ERROR;
  }
  if(ok && (ok = APIEnterNotModal(G))) {
    SceneSetMatrix(G, m);       /* TODO STATUS */
    APIExit(G);
  }
  return APIResultOk(ok);
}
#endif

static PyObject *CmdGetMinMax(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  float mn[3], mx[3];
  char *str1;
  int state;
  OrthoLineType s1;
  int flag;

  int ok = false;
  API_SETUP_ARGS(G, self, args, "Osi", &self, &str1, &state);
  APIEnter(G);
  {
    ok = (SelectorGetTmp2(G, str1, s1) >= 0);
    flag = ExecutiveGetExtent(G, s1, mn, mx, true, state, false);
    SelectorFreeTmp(G, s1);
    APIExit(G);
    if(flag)
      return Py_BuildValue("[[fff],[fff]]", mn[0], mn[1], mn[2], mx[0], mx[1], mx[2]);
  }
  return Py_BuildValue("[[fff],[fff]]", -0.5, -0.5, -0.5, 0.5, 0.5, 0.5);
}

#if 0
static PyObject *CmdGetMatrix(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  float *f;
  PyObject *result = NULL;
  int ok = false;
  ok = PyArg_ParseTuple(args, "O", &self);
  if(ok) {
    API_SETUP_PYMOL_GLOBALS;
    ok = (G != NULL);
  } else {
    API_HANDLE_ERROR;
  }
  if(ok && (ok = APIEnterNotModal(G))) {
    f = SceneGetMatrix(G);
    APIExit(G);
    result = Py_BuildValue("ffffffffffffffff",
                           f[0], f[1], f[2], f[3],
                           f[4], f[5], f[6], f[7],
                           f[8], f[9], f[10], f[11], f[12], f[13], f[14], f[15]
      );
  }
  return APIAutoNone(result);
}
#endif

static PyObject *CmdGetObjectMatrix(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  PyObject *result = NULL;
  char *name;
  double *history = NULL;
  int ok = false;
  int found;
  int state;
  int incl_ttt = true;
  ok = PyArg_ParseTuple(args, "Osi|i", &self, &name, &state, &incl_ttt);
  if(ok) {
    API_SETUP_PYMOL_GLOBALS;
    ok = (G != NULL);
  } else {
    API_HANDLE_ERROR;
  }
  if(ok) {
    APIEnter(G);
    found = ExecutiveGetObjectMatrix(G, name, state, &history, incl_ttt);
    APIExit(G);
    if(found) {
      if(history)
        result = Py_BuildValue("dddddddddddddddd",
                               history[0], history[1], history[2], history[3],
                               history[4], history[5], history[6], history[7],
                               history[8], history[9], history[10], history[11],
                               history[12], history[13], history[14], history[15]
          );
      else
        result = Py_BuildValue("dddddddddddddddd",
                               1.0, 0.0, 0.0, 0.0,
                               0.0, 1.0, 0.0, 0.0,
                               0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0);
    }
  }
  return APIAutoNone(result);
}

static PyObject *CmdGetObjectTTT(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  PyObject *result = NULL;
  const char *name;
  int state, quiet;
  const float *ttt = NULL;

  if (!PyArg_ParseTuple(args, "Osii", &self, &name, &state, &quiet)) {
    API_HANDLE_ERROR;
    ok_raise(1);
  }

  API_SETUP_PYMOL_GLOBALS;
  ok_assert(1, G);
  APIEnter(G);

  ExecutiveGetObjectTTT(G, name, &ttt, state, quiet);
  if (ttt)
    result = PConvFloatArrayToPyList(ttt, 16);

  APIExit(G);
ok_except1:
  return APIAutoNone(result);
}

static PyObject *CmdMDo(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  char *cmd;
  int frame;
  int append;
  int ok = false;
  ok = PyArg_ParseTuple(args, "Oisi", &self, &frame, &cmd, &append);
  if(ok) {
    API_SETUP_PYMOL_GLOBALS;
    ok = (G != NULL);
  } else {
    API_HANDLE_ERROR;
  }
  if(ok && (ok = APIEnterNotModal(G))) {
    if(frame < 0) {
      if(frame == -1) {
        frame = SceneGetFrame(G);
      } else {
        frame = MovieGetLength(G) + 2 + frame;
        if(frame<0)
          frame = 0;
      }
    }
    if(append) {
      MovieAppendCommand(G, frame, cmd);
    } else {
      MovieSetCommand(G, frame, cmd);
    }
    APIExit(G);
  }
  return APIResultOk(ok);
}

static PyObject *CmdMPlay(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  int cmd;
  int ok = false;
  ok = PyArg_ParseTuple(args, "Oi", &self, &cmd);
  if(ok) {
    API_SETUP_PYMOL_GLOBALS;
    ok = (G != NULL);
  } else {
    API_HANDLE_ERROR;
  }
  if(ok && (ok = APIEnterNotModal(G))) {
    MoviePlay(G, cmd);
    APIExit(G);
  }
  return APIResultOk(ok);
}

static PyObject *CmdMMatrix(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  int cmd;
  int ok = false;
  ok = PyArg_ParseTuple(args, "Oi", &self, &cmd);
  if(ok) {
    API_SETUP_PYMOL_GLOBALS;
    ok = (G != NULL);
  } else {
    API_HANDLE_ERROR;
  }
  if(ok && (ok = APIEnterNotModal(G))) {
    ok = MovieMatrix(G, cmd);
    APIExit(G);
  }
  return APIResultOk(ok);
}

static PyObject *CmdMClear(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  int ok = false;
  ok = PyArg_ParseTuple(args, "O", &self);
  if(ok) {
    API_SETUP_PYMOL_GLOBALS;
    ok = (G != NULL);
  } else {
    API_HANDLE_ERROR;
  }
  if(ok && (ok = APIEnterNotModal(G))) {
    MovieClearImages(G);
    APIExit(G);
  }
  return APISuccess();
}

static PyObject *CmdRefreshLater(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  int ok = false;
  ok = PyArg_ParseTuple(args, "O", &self);
  if(ok) {
    API_SETUP_PYMOL_GLOBALS;
    ok = (G != NULL);
  } else {
    API_HANDLE_ERROR;
  }
  if(ok && (ok = APIEnterNotModal(G))) {
    SceneInvalidate(G);
    APIExit(G);
  }
  return APISuccess();
}

static PyObject *CmdRefresh(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  int ok = false;
  ok = PyArg_ParseTuple(args, "O", &self);
  if(ok) {
    API_SETUP_PYMOL_GLOBALS;
    ok = (G != NULL);
  } else {
    API_HANDLE_ERROR;
  }
  if(ok && (ok = APIEnterNotModal(G))) {
    SceneInvalidateCopy(G, false);
    ExecutiveDrawNow(G);        /* TODO STATUS */
    APIExit(G);
  }
  return APISuccess();
}

static PyObject *CmdRefreshNow(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  int ok = false;
  ok = PyArg_ParseTuple(args, "O", &self);
  if(ok) {
    API_SETUP_PYMOL_GLOBALS;
    ok = (G != NULL);
  } else {
    API_HANDLE_ERROR;
  }
  if(ok && (ok = APIEnterNotModal(G))) {
    PyMOL_PushValidContext(G->PyMOL);   /* we're trusting the caller on this... */
    SceneInvalidateCopy(G, false);
    ExecutiveDrawNow(G);        /* TODO STATUS */
#ifndef _PYMOL_NO_MAIN
    if(G->Main) {
      MainRefreshNow();
    }
#endif
    PyMOL_PopValidContext(G->PyMOL);
    APIExit(G);
  }
  return APISuccess();
}

static PyObject* CmdPushValidContext(PyObject* self, PyObject* args)
{
  assert(PIsGlutThread());
  PyMOLGlobals* G = nullptr;
  API_SETUP_ARGS(G, self, args, "O", &self);
  PyMOL_PushValidContext(G->PyMOL);
  return APISuccess();
}

static PyObject* CmdPopValidContext(PyObject* self, PyObject* args)
{
  assert(PIsGlutThread());
  PyMOLGlobals* G = nullptr;
  API_SETUP_ARGS(G, self, args, "O", &self);
  PyMOL_PopValidContext(G->PyMOL);
  return APISuccess();
}

static PyObject *CmdPNG(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  char* filename = nullptr;
  int quiet;
  int result = 0;
  int width, height, ray;
  int prior, format;
  float dpi;

  API_SETUP_ARGS(G, self, args, "Oziifiiii", &self, &filename, &width, &height,
      &dpi, &ray, &quiet, &prior, &format);
  API_ASSERT(APIEnterNotModal(G));

  // if `filename` is None, then return a PNG buffer
  std::vector<unsigned char> pngbuf;

  {
    // with prior=1 other arguments (width, height, ray) are ignored

    if(!prior) {
      if(ray || (!G->HaveGUI && (!SceneGetCopyType(G) || width || height))) {
        prior = SceneRay(G, width, height, SettingGetGlobal_i(G, cSetting_ray_default_renderer),
                 NULL, NULL, 0.0F, 0.0F, false, NULL, true, -1);
      } else if(width || height) {
        prior = !SceneDeferImage(
            G, width, height, filename, -1, dpi, quiet, format);
        result = bool(filename);
      } else if(!SceneGetCopyType(G)) {
        ExecutiveDrawNow(G);      /* TODO STATUS */
      }
    }

    if(!result) {
      if (ScenePNG(G, filename, dpi, quiet, prior, format,
              filename ? nullptr : &pngbuf))
        result = 1;             /* signal success by returning 1 instead of 0, or -1 for error  */
    }
    APIExit(G);
  }

  if (!filename) {
    if (pngbuf.empty()) {
      return APIFailure(G, "getting png buffer failed");
    }

    return PyBytes_FromStringAndSize(
        reinterpret_cast<const char*>(pngbuf.data()), pngbuf.size());
  }

  return APIResultCode(result);
}

static PyObject *CmdMPNG(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  char *str1;
  int int1, int2, int3, int4, format, mode, quiet;
  int ok = false;
  int width = 0, height = 0;
  ok = PyArg_ParseTuple(args, "Osiiiiiiiii", &self, &str1, &int1, &int2,
                        &int3, &int4, &format, &mode, &quiet,
                        &width, &height);
  if(ok) {
    API_SETUP_PYMOL_GLOBALS;
    ok = (G != NULL);
  } else {
    API_HANDLE_ERROR;
  }
  if(ok && (ok = APIEnterNotModal(G))) {
    ok = MoviePNG(G, str1, SettingGetGlobal_b(G, cSetting_cache_frames),
                  int1, int2, int3, int4, format, mode, quiet,
                  width, height);
    /* TODO STATUS */
    APIExit(G);
  }
  return APIResultOk(G, ok);
}

static PyObject *CmdMSet(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  char *str1;
  int start_from,freeze;
  API_SETUP_ARGS(G, self, args, "Osii", &self, &str1, &start_from,&freeze);
  API_ASSERT(APIEnterNotModal(G));
  MovieSet(G, str1, start_from, freeze);
  APIExit(G);
  return APISuccess();
}

static PyObject *CmdMModify(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  char *object;
  int action,index,count,target,freeze,quiet;
  API_SETUP_ARGS(G, self, args, "Oiiiisii", &self, &action, &index, &count,
      &target, &object, &freeze, &quiet);
  API_ASSERT(APIEnterNotModal(G));
  auto result = ExecutiveMotionViewModify(G,action,index,count,target,object,freeze,quiet);
  APIExit(G);
  return APIResult(G, result);
}

static PyObject *CmdMView(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  int ok = false;
  int action, first, last, simple, wrap, window, cycles, quiet, state, autogen;
  float power, bias, linear, hand, scene_cut;
  char *object, *scene_name;
  ok = PyArg_ParseTuple(args, "Oiiiffifsiiiisfiii", &self, &action, &first, &last, &power,
                        &bias, &simple, &linear, &object, &wrap, &hand,
                        &window, &cycles, &scene_name, &scene_cut, &quiet, &state, &autogen);
  if(ok) {
    API_SETUP_PYMOL_GLOBALS;
    ok = (G != NULL);
  } else {
    API_HANDLE_ERROR;
  }
  if(ok && (ok = APIEnterNotModal(G))) {
    ok = ExecutiveMotionView(G, action, first, last, power, bias, simple, 
                             linear, object, wrap, hand, window, cycles, 
                             scene_name, scene_cut, state, quiet, autogen);
    APIExit(G);
  }
  return APIResultOk(ok);
}

static PyObject *CmdViewport(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  int w, h;
  API_SETUP_ARGS(G, self, args, "Oii", &self, &w, &h);
  API_ASSERT(APIEnterNotModal(G));
  {
    {
      if(!(w < 1 && h < 1 && ExecutiveIsFullScreen(G))) {

        if(((w > 0) && (h <= 0)) || ((h > 0) && (w <= 0))) {
          int cw, ch;
          SceneGetWidthHeight(G, &cw, &ch);
          if(h <= 0) {
            h = (w * ch) / cw;
          }
          if(w <= 0) {
            w = (h * cw) / ch;
          }
        }

        if((w > 0) && (h > 0)) {
          if(w < 10)
            w = 10;
          if(h < 10)
            h = 10;

          if(SettingGetGlobal_b(G, cSetting_internal_gui)) {
            w += DIP2PIXEL(SettingGetGlobal_i(G, cSetting_internal_gui_width));
          }

          if(SettingGetGlobal_i(G, cSetting_internal_feedback)) {
            h += (SettingGetGlobal_i(G, cSetting_internal_feedback) - 1) * cOrthoLineHeight +
              cOrthoBottomSceneMargin;
          }

          h += MovieGetPanelHeight(G);
        } else {
          w = -1;
          h = -1;
        }
      } else {
        w = 0;
        h = 0;
      }

#ifndef _PYMOL_NO_MAIN
        if(G->Main) {
          MainDoReshape(w, h);    /* should be moved into Executive */
        }
      else
#endif
      {
        PyMOL_NeedReshape(G->PyMOL, 2, 0, 0, w, h);
      }
      APIExit(G);
    }

  }
  return APISuccess();
}

static PyObject *CmdFlag(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  char *str1;
  int flag;
  int action;
  int quiet;
  API_SETUP_ARGS(G, self, args, "Oisii", &self, &flag, &str1, &action, &quiet);
  API_ASSERT(APIEnterNotModal(G));
  auto result = ExecutiveFlag(G, flag, str1, action, quiet);
  APIExit(G);
  return APIResult(G, result);
}

static PyObject *CmdColor(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  char *str1, *color;
  int flags;
  int quiet;

  API_SETUP_ARGS(G, self, args, "Ossii", &self, &color, &str1, &flags, &quiet);
  API_ASSERT(APIEnterNotModal(G));
  auto result = ExecutiveColorFromSele(G, str1, color, flags, quiet);
  APIExit(G);
  return APIResult(G, result);
}

static PyObject *CmdColorDef(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  const char* color;
  float v[3];
  int mode;
  int quiet;
  API_SETUP_ARGS(G, self, args, "Osfffii", &self, &color, v, v + 1, v + 2, &mode, &quiet);
  API_ASSERT(APIEnterNotModal(G));
  ColorDef(G, color, v, mode, quiet);
  APIExit(G);
  return APISuccess();
}

static PyObject *CmdDraw(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  int int1, int2;
  int quiet, antialias;

  API_SETUP_ARGS(G, self, args, "Oiiii", &self, &int1, &int2, &antialias, &quiet);
  API_ASSERT(APIEnterNotModal(G));

  bool entire_window = false;

  if (antialias == -2) {
    // capture action
    entire_window = true;
    int1 = 0;
    int2 = 0;
    antialias = 0;
  }

  bool ok = ExecutiveDrawCmd(G, int1, int2, antialias, entire_window, quiet);
  APIExit(G);
  return APIResultOk(G, ok);
}

static PyObject *CmdRay(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  int w, h, mode;
  float angle, shift;
  int quiet;
  int antialias;
  API_SETUP_ARGS(G, self, args, "Oiiiffii", &self, &w, &h,
                        &antialias, &angle, &shift, &mode, &quiet);
  API_ASSERT(APIEnterNotModal(G));
  {
    if(mode < 0)
      mode = SettingGetGlobal_i(G, cSetting_ray_default_renderer);
    ExecutiveRay(G, w, h, mode, angle, shift, quiet, false, antialias); /* TODO STATUS */
    APIExit(G);
  }
  return APISuccess();
}

static PyObject *CmdClip(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  char *sname;
  float dist;
  char *str1;
  int state;
  API_SETUP_ARGS(G, self, args, "Osfsi", &self, &sname, &dist, &str1, &state);
  API_ASSERT(APIEnterNotModal(G));
  SelectorTmp2 s1(G, str1);
  auto result = SceneClipFromMode(G, sname, dist, s1.getName(), state);
  APIExit(G);
  return APIResult(G, result);
}

static PyObject *CmdMove(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  char *sname;
  float dist;
  API_SETUP_ARGS(G, self, args, "Osf", &self, &sname, &dist);
  API_ASSERT(APIEnterNotModal(G));
  auto result = ExecutiveMove(G, sname, dist);
  APIExit(G);
  return APIResult(G, result);
}

static PyObject *CmdTurn(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  char *sname;
  float angle;
  API_SETUP_ARGS(G, self, args, "Osf", &self, &sname, &angle);
  API_ASSERT(APIEnterNotModal(G));
  SceneRotateAxis(G, angle, sname[0]);
  APIExit(G);
  return APISuccess();
}

static PyObject *CmdUnset(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  int index;
  char *str;
  int state;
  int quiet;
  int updates;
  API_SETUP_ARGS(G, self, args, "Oisiii", &self, &index, &str, &state, &quiet, &updates);
  API_ASSERT(APIEnterNotModal(G));
  auto result = ExecutiveUnsetSetting(G, index, str, state, quiet, updates);
  APIExit(G);
  return APIResult(G, result);
}

static PyObject *CmdUnsetBond(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  int index;
  char *str3, *str4;
  int state;
  int quiet;
  int updates;
  API_SETUP_ARGS(G, self, args, "Oissiii", &self, &index, &str3, &str4, &state, &quiet,
                     &updates);
  API_ASSERT(APIEnterNotModal(G));
  // TODO move selection handling to Executive
  auto res = [&]() -> pymol::Result<> {
    auto tmpsele1 = SelectorTmp::make(G, str3);
    p_return_if_error(tmpsele1);
    auto tmpsele2 = SelectorTmp::make(G, str4);
    p_return_if_error(tmpsele2);
    ExecutiveUnsetBondSetting(G, index, tmpsele1->getName(),
        tmpsele2->getName(), state, quiet, updates);
    return {};
  }();
  APIExit(G);
  return APIResult(G, res);
}

static PyObject *CmdSet(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  int index;
  PyObject *value;
  char *str;
  int state;
  int quiet;
  int updates;
  API_SETUP_ARGS(G, self, args, "OiOsiii", &self, &index, &value, &str, &state, &quiet,
                 &updates);
  API_ASSERT(APIEnterNotModal(G));
  auto result = ExecutiveSetSetting(G, index, value, str, state, quiet, updates);
  APIExit(G);
  return APIResult(G, result);
}

static PyObject *CmdSetBond(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  int index;
  PyObject *value;
  char *str3, *str4;
  int state;
  int quiet;
  int updates;
  API_SETUP_ARGS(G, self, args, "OiOssiii", &self, &index, &value, &str3, &str4, &state,
                     &quiet, &updates);
  API_ASSERT(APIEnterNotModal(G));
  // TODO move selection handling to Executive
  auto res = [&]() -> pymol::Result<> {
    auto tmpsele1 = SelectorTmp::make(G, str3);
    p_return_if_error(tmpsele1);
    auto tmpsele2 = SelectorTmp::make(G, str4);
    p_return_if_error(tmpsele2);
    if (!ExecutiveSetBondSetting(G, index, value, tmpsele1->getName(),
            tmpsele2->getName(), state, quiet, updates)) {
      return pymol::Error();
    }
    return {};
  }();
  APIExit(G);
  return APIResult(G, res);
}

static PyObject *CmdGetBond(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  int index;
  char *str3, *str4;
  int state;
  int quiet;
  int updates;
  API_SETUP_ARGS(G, self, args, "Oissiii", &self, &index, &str3, &str4, &state,
                     &quiet, &updates);
  APIEnterBlocked(G);
  // TODO move selection handling to Executive
  auto res = [&]() -> pymol::Result<PyObject*> {
    auto tmpsele1 = SelectorTmp::make(G, str3);
    p_return_if_error(tmpsele1);
    auto tmpsele2 = SelectorTmp::make(G, str4);
    p_return_if_error(tmpsele2);
    return ExecutiveGetBondSetting(G, index,
        (/* TODO */ char*) tmpsele1->getName(), tmpsele2->getName(), state,
        quiet, updates);
  }();
  APIExitBlocked(G);
  return APIResult(G, res);
}

static PyObject *CmdDelete(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  char *sname;
  API_SETUP_ARGS(G, self, args, "Os", &self, &sname);
  API_ASSERT(APIEnterNotModal(G));
  ExecutiveDelete(G, sname);
  APIExit(G);
  return APISuccess();
}

static PyObject *CmdCartoon(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  char *sname;
  int type;
  API_SETUP_ARGS(G, self, args, "Osi", &self, &sname, &type);
  API_ASSERT(APIEnterNotModal(G));
  auto result = ExecutiveCartoon(G, type, sname);
  APIExit(G);
  return APIResult(G, result);
}

static PyObject *CmdShowHide(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  const char *sname;
  int rep;
  int state;
  API_SETUP_ARGS(G, self, args, "Osii", &self, &sname, &rep, &state);
  API_ASSERT(APIEnterNotModal(G));
  if(sname[0] == '@') {
    // DEPRECATED
    sname = cKeywordAll;
    rep = cRepBitmask;
  }
  pymol::Result<> res;
  {
    auto tmpsele1 = SelectorTmp2::make(G, sname);
    if (tmpsele1) {
      ExecutiveSetRepVisMask(G, tmpsele1->getName(), rep, state);
    } else {
      res = tmpsele1.error_move();
    }
  }
  APIExit(G);
  return APIResult(G, res);
}

static PyObject *CmdOnOffBySele(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  char *sname;
  int onoff;
  API_SETUP_ARGS(G, self, args, "Osi", &self, &sname, &onoff);
  API_ASSERT(APIEnterNotModal(G));
  ExecutiveSetOnOffBySele(G, sname, onoff);
  APIExit(G);
  return APISuccess();
}

static PyObject *CmdOnOff(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  char *name;
  int state;
  int parents = 0;
  API_SETUP_ARGS(G, self, args, "Osii", &self, &name, &state, &parents);
  API_ASSERT(APIEnterNotModal(G));
  ExecutiveSetObjVisib(G, name, state, parents);
  APIExit(G);
  return APISuccess();
}

static PyObject *CmdToggle(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  char *sname;
  int rep;
  API_SETUP_ARGS(G, self, args, "Osi", &self, &sname, &rep);
  API_ASSERT(APIEnterNotModal(G));
  auto result = ExecutiveToggleRepVisib(G, sname, rep);
  APIExit(G);
  return APIResult(G, result);
}

static PyObject *CmdQuit(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  int ok = false;
  int code = EXIT_SUCCESS;
  ok = PyArg_ParseTuple(args, "O|i", &self, &code);
  if(ok) {
    API_SETUP_PYMOL_GLOBALS;
    ok = (G != NULL);
  } else {
    API_HANDLE_ERROR;
  }
  if(ok) {
    if(!APIEnterNotModal(G)) {  /* override modal state to enable quit */
      PyMOL_SetModalDraw(G->PyMOL, NULL);
      APIEnter(G);
    }
    if(!G->Option->no_quit) {
      G->Terminating = true;
      PExit(G, code);
    } else {
      OrthoAddOutput(G, "Cmd-Error: cannot quit from within this context.\n");
    }
    APIExit(G);
  }
  return APISuccess();
}

static PyObject *CmdFullScreen(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  int flag = 0;
  int ok = false;
  API_SETUP_ARGS(G, self, args, "Oi", &self, &flag);
  API_ASSERT(APIEnterNotModal(G));
    ExecutiveFullScreen(G, flag);       /* TODO STATUS */
  APIExit(G);
  return APIResultOk(G, ok);
}

static PyObject *CmdGroup(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  char *gname, *names;
  int quiet, action;
  int ok = false;
  API_SETUP_ARGS(G, self, args, "Ossii", &self, &gname, &names, &action, &quiet);
  API_ASSERT(APIEnterNotModal(G));
  {
    ok = ExecutiveGroup(G, gname, names, action, quiet);
    APIExit(G);
  }
  return APIResultOk(G, ok);
}

static PyObject *CmdSelect(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  char *sname, *sele;
  int quiet;
  int state = 0;
  char *domain;
  int enable = -1;
  int merge = 0;
  API_SETUP_ARGS(G, self, args, "Ossiis|ii", &self, &sname, &sele, &quiet,
      &state, &domain, &enable, &merge);
  API_ASSERT(APIEnterNotModal(G));

  auto seleargs = ExecutiveSelectPrepareArgs(G, sname, sele);
  auto res =
      ExecutiveSelect(G, seleargs, enable, quiet, merge, state, domain);

  APIExit(G);
  return APIResult(G, res);
}

static PyObject *CmdFinishObject(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  const char* oname;

  API_SETUP_ARGS(G, self, args, "Os", &self, &oname);
  API_ASSERT(APIEnterNotModal(G));

  auto origObj = ExecutiveFindObjectByName(G, oname);
  auto ok = bool(origObj);
  if(ok) {
      if(origObj->type == cObjectMolecule) {
        ObjectMoleculeUpdateIDNumbers((ObjectMolecule *) origObj);
        ObjectMoleculeUpdateNonbonded((ObjectMolecule *) origObj);
        origObj->invalidate(cRepAll, cRepInvAll, -1);
      }
      ExecutiveUpdateObjectSelection(G, origObj);       /* TODO STATUS */
  }
  APIExit(G);
  return APIResultOk(G, ok);
}

static PyObject *CmdLoadObject(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  char *oname;
  PyObject *model;
  int frame, type;
  int finish, discrete;
  int quiet;
  int zoom;
  API_SETUP_ARGS(G, self, args, "OsOiiiiii", &self, &oname, &model, &frame,
      &type, &finish, &discrete, &quiet, &zoom);
  API_ASSERT(APIEnterNotModal(G));
  ExecutiveLoadObject(G, oname, model, frame, type, finish, discrete, quiet, zoom);
  APIExit(G);
  return APISuccess();
}

static PyObject *CmdSetStateOrder(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  char *oname;
  PyObject *order;
  pymol::CObject *obj = NULL;
  int *int_array = NULL;
  int ok = false;

  if(!PyArg_ParseTuple(args, "OsO", &self, &oname, &order)) {
    API_HANDLE_ERROR;
    ok_raise(1);
  }

  ok_assert(1, PyList_Check(order));

  API_SETUP_PYMOL_GLOBALS;
  ok_assert(1, G && APIEnterNotModal(G));

  obj = ExecutiveFindObjectByName(G, oname);
  if(!obj || obj->type != cObjectMolecule) {
    ErrMessage(G, "SetStateOrder", "named object molecule not found.");
    ok_raise(2);
  }

  if(PConvPyListToIntArray(order, &int_array)) {
    int len = PyList_Size(order);

    PBlock(G);
    ok = ObjectMoleculeSetStateOrder((ObjectMolecule *) obj, int_array, len);
    PUnblock(G);

    FreeP(int_array);
  } else {
    ErrMessage(G, "SetStateOrder", "not an integer list.");
    ok_raise(2);
  }

  APIExit(G);
  return APIResultOk(ok);
ok_except2:
  APIExit(G);
ok_except1:
  return APIFailure();
}

static PyObject *CmdLoadCoords(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  char *str1;
  int state = 0;
  OrthoLineType s1;
  PyObject *coords = NULL;

  API_SETUP_ARGS(G, self, args, "OsO|i", &self, &str1, &coords, &state);
  pymol::Result<> result;
  if(!str1[0]) {
    return APIFailure(G, "selection is empty");
  }
  API_ASSERT(APIEnterBlockedNotModal(G));

  if(SelectorGetTmp(G, str1, s1) >= 0) {
    int sele1 = SelectorIndexByName(G, s1);
    if(sele1 >= 0) {
      result = SelectorLoadCoords(G, coords, sele1, state);
    }
    SelectorFreeTmp(G, s1);
  }

  APIExitBlocked(G);

  if (!result && PyErr_Occurred()) {
    return nullptr;
  }

  return APIResult(G, result);
}

static PyObject *CmdLoadCoordSet(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  const char* oname;
  PyObject *model;
  int frame;

  API_SETUP_ARGS(G, self, args, "OsOi", &self, &oname, &model, &frame);
  API_ASSERT(APIEnterNotModal(G));
  auto result = ExecutiveLoadCoordset(G, oname, model, frame);
  APIExit(G);
  return APISuccess();
}

static PyObject *CmdLoad(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  char *fname, *oname;
  char *object_props = NULL;
  char *atom_props = NULL;
  char *plugin = NULL;
  int frame;
  cLoadType_t type;
  int finish, discrete;
  int quiet;
  int multiplex;
  int zoom;
  Py_ssize_t bytes;
  int mimic;
  const char* contents;

  API_SETUP_ARGS(G, self, args, "Oszz#iiiiiii|zzzi", &self,
                        &oname, &fname, &contents, &bytes, &frame, &type,
                        &finish, &discrete, &quiet, &multiplex, &zoom,
                        &plugin, &object_props, &atom_props, &mimic);
  API_ASSERT(APIEnterNotModal(G));

  auto result = ExecutiveLoad(G,
                         fname, contents, bytes, type,
                         oname, frame, zoom,
                         discrete, finish, multiplex, quiet, plugin);

  OrthoRestorePrompt(G);
  APIExit(G);
  
  return APIResult(G, result);
}

static PyObject *CmdLoadTraj(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  char *fname, *oname;
  int frame, type;
  int interval, average, start, stop, max, image;
  char *str1;
  float shift[3];
  char *plugin = NULL;
  int quiet = 0;                /* TODO */
  API_SETUP_ARGS(G, self, args, "Ossiiiiiiisifffs", &self, &oname, &fname,
      &frame, &type, &interval, &average, &start, &stop, &max, &str1, &image,
      &shift[0], &shift[1], &shift[2], &plugin);
  API_ASSERT(APIEnterNotModal(G));
  auto result = ExecutiveLoadTraj(G, oname, fname, frame, type,
      interval, average, start, stop, max, str1, image, shift, plugin, quiet);
  APIExit(G);
  return APIResult(G, result);
}

static PyObject *CmdOrigin(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  char *str1, *obj;
  float v[3];
  int state;
  API_SETUP_ARGS(G, self, args, "Oss(fff)i", &self, &str1, &obj, v, v + 1, v + 2, &state);
  API_ASSERT(APIEnterNotModal(G));
  auto result = ExecutiveOrigin(G, str1, 1, obj, v, state);
  APIExit(G);
  return APIResult(G, result);
}

static PyObject *CmdSort(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  char *name;
  API_SETUP_ARGS(G, self, args, "Os", &self, &name);
  API_ASSERT(APIEnterNotModal(G));
  auto result = ExecutiveSort(G, name);
  APIExit(G);
  return APIResult(G, result);
}

static PyObject *CmdAssignSS(PyObject * self, PyObject * args)


/* EXPERIMENTAL */
{
  PyMOLGlobals *G = NULL;
  int state, quiet;
  char *str1, *str2;
  int preserve;
  API_SETUP_ARGS(G, self, args, "Osisii", &self, &str1, &state, &str2, &preserve, &quiet);
  API_ASSERT(APIEnterNotModal(G));
  auto result =
      ExecutiveAssignSS(G, str1, state, str2, preserve, nullptr, quiet);
  APIExit(G);
  return APIResult(G, result);
}

static PyObject *CmdSpheroid(PyObject * self, PyObject * args)


/* EXPERIMENTAL */
{
  PyMOLGlobals *G = NULL;
  char *name;
  int average;
  API_SETUP_ARGS(G, self, args, "Osi", &self, &name, &average);
  API_ASSERT(APIEnterNotModal(G));
  auto result = ExecutiveSpheroid(G, name, average);
  APIExit(G);
  return APIResult(G, result);
}

static PyObject *CmdTest(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  /* regression tests */

  int ok = false;
  int code;
  int group;

  ok = PyArg_ParseTuple(args, "Oii", &self, &group, &code);
  if(ok) {
    API_SETUP_PYMOL_GLOBALS;
    ok = (G != NULL);
  } else {
    API_HANDLE_ERROR;
  }
  if(ok && (ok = APIEnterNotModal(G))) {
    PRINTFB(G, FB_CCmd, FB_Details)
      " Cmd: initiating test %d-%d.\n", group, code ENDFB(G);
    ok = TestPyMOLRun(G, group, code);
    PRINTFB(G, FB_CCmd, FB_Details)
      " Cmd: concluding test %d-%d.\n", group, code ENDFB(G);
    APIExit(G);
  }
  return APIResultOk(ok);
}

static PyObject *CmdCenter(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  char *str1;
  int state;
  int origin;
  float animate;
  int quiet = false;            /* TODO */
  API_SETUP_ARGS(G, self, args, "Osiif", &self, &str1, &state, &origin, &animate);
  API_ASSERT(APIEnterNotModal(G));
  pymol::Result<> res;
  {
    auto tmpsele1 = SelectorTmp2::make(G, str1);
    if (tmpsele1) {
      ExecutiveCenter(G, tmpsele1->getName(), state, origin, animate, nullptr, quiet);
    } else {
      res = tmpsele1.error_move();
    }
  }
  APIExit(G);
  return APIResult(G, res);
}

static PyObject *CmdReference(PyObject * self, PyObject * args)
{
  int ok = false;
  PyMOLGlobals *G = NULL;
  OrthoLineType s1;
  int action;
  char *sele1;
  int state;
  int quiet;
  ok = PyArg_ParseTuple(args, "Oisii", &self, &action, &sele1, &state, &quiet);
  if(ok) {
    API_SETUP_PYMOL_GLOBALS;
    ok = (G != NULL);
  } else {
    API_HANDLE_ERROR;
  }
  if(ok && (ok = APIEnterNotModal(G))) {
    ok = (SelectorGetTmp(G, sele1, s1) >= 0);
    if(ok)
      ok = ExecutiveReference(G, action, s1, state, quiet);
    SelectorFreeTmp(G, s1);
    APIExit(G);
  }
  return APIResultOk(ok);
}

static PyObject *CmdZoom(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  char *str1;
  float buffer;
  int state;
  int inclusive;
  float animate;
  int quiet = false;            /* TODO */
  API_SETUP_ARGS(G, self, args, "Osfiif", &self, &str1, &buffer, &state, &inclusive, &animate);
  API_ASSERT(APIEnterNotModal(G));
  SelectorTmp2 s1(G, str1);
  ExecutiveWindowZoom(G, s1.getName(), buffer, state, inclusive, animate, quiet);
  APIExit(G);
  return APISuccess();
}

static PyObject *CmdIsolevel(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  float level;
  int state;
  char *name;
  int query, quiet;
  API_SETUP_ARGS(G, self, args, "Osfiii", &self, &name, &level, &state, &query, &quiet);
  if(!query) {
    API_ASSERT(APIEnterNotModal(G));
    auto result = ExecutiveIsolevel(G, name, level, state, quiet);
    APIExit(G);
    return APIResult(G, result);
  } else {
    APIEnter(G);
    auto result = ExecutiveGetIsolevel(G, name, state);
    APIExit(G);
    return APIResult(G, result);
  }
}

static PyObject *CmdHAdd(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  char *str1;
  int quiet;
  int state;
  int legacy;
  API_SETUP_ARGS(G, self, args, "Osiii", &self, &str1, &quiet, &state, &legacy);
  API_ASSERT(APIEnterNotModal(G));
  auto result = ExecutiveAddHydrogens(G, str1, quiet, state, legacy);
  APIExit(G);
  return APIResult(G, result);
}

static PyObject *CmdSetObjectColor(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  char *name, *color;
  int quiet;
  API_SETUP_ARGS(G, self, args, "Ossi", &self, &name, &color, &quiet);
  API_ASSERT(APIEnterNotModal(G));
  auto result = ExecutiveSetObjectColor(G, name, color, quiet);
  APIExit(G);
  return APIResult(G, result);
}

static PyObject *CmdGetObjectColorIndex(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  char *str1;
  API_SETUP_ARGS(G, self, args, "Os", &self, &str1);
  APIEnter(G);
  auto result = ExecutiveGetObjectColorIndex(G, str1);
  APIExit(G);
  return (APIResultCode(result));
}

static PyObject *CmdRemove(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  char *str1;
  int quiet;
  API_SETUP_ARGS(G, self, args, "Osi", &self, &str1, &quiet);
  API_ASSERT(APIEnterNotModal(G));
  auto result = ExecutiveRemoveAtoms(G, str1, quiet);
  APIExit(G);
  return APIResult(G, result);
}

static PyObject *CmdRemovePicked(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  int i1;
  int quiet;
  API_SETUP_ARGS(G, self, args, "Oii", &self, &i1, &quiet);
  API_ASSERT(APIEnterNotModal(G));
  auto result = EditorRemove(G, i1, quiet);
  APIExit(G);
  return APIResult(G, result);
}

static PyObject *CmdHFill(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  int quiet;
  API_SETUP_ARGS(G, self, args, "Oi", &self, &quiet);
  API_ASSERT(APIEnterNotModal(G));
  auto result = EditorHFill(G, quiet);
  APIExit(G);
  return APIResult(G, result);
}

static PyObject *CmdHFix(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  int quiet;
  char *str1;
  API_SETUP_ARGS(G, self, args, "Osi", &self, &str1, &quiet);
  API_ASSERT(APIEnterNotModal(G));
  SelectorTmp2 s1(G, str1);
  auto result = EditorHFix(G, s1.getName(), quiet);
  APIExit(G);
  return APIResult(G, result);
}

static PyObject *CmdCycleValence(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  int quiet;
  API_SETUP_ARGS(G, self, args, "Oi", &self, &quiet);
  API_ASSERT(APIEnterNotModal(G));
  auto result = EditorCycleValence(G, quiet);
  APIExit(G);
  return APIResult(G, result);
}

static PyObject *CmdReplace(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  int i1, i2;
  char *str1, *str2;
  int quiet;
  API_SETUP_ARGS(G, self, args, "Osiisi", &self, &str1, &i1, &i2, &str2, &quiet);
  API_ASSERT(APIEnterNotModal(G));
  auto result = EditorReplace(G, str1, i1, i2, str2, quiet);
  APIExit(G);
  return APIResult(G, result);
}

static PyObject *CmdSetGeometry(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  int geom, valence;
  char *sele;
  API_SETUP_ARGS(G, self, args, "Osii", &self, &sele, &geom, &valence);
  API_ASSERT(APIEnterNotModal(G));
  auto result = ExecutiveSetGeometry(G, sele, geom, valence);
  APIExit(G);
  return APIResult(G, result);
}

static PyObject *CmdAttach(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  int i1, i2;
  char *str1;
  int quiet;
  char *name;
  API_SETUP_ARGS(G, self, args, "Osiis", &self, &str1, &i1, &i2, &name, &quiet);
  API_ASSERT(APIEnterNotModal(G));
  auto result = EditorAttach(G, str1, i1, i2, name, quiet);
  APIExit(G);
  return APIResult(G, result);
}

static PyObject *CmdFuse(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  char *str1, *str2;
  int mode;
  int recolor;
  int move_flag;
  API_SETUP_ARGS(G, self, args, "Ossiii", &self, &str1, &str2, &mode, &recolor,
      &move_flag);
  API_ASSERT(APIEnterNotModal(G));
  auto result = ExecutiveFuse(
      G, str1, str2, mode, recolor, move_flag);
  APIExit(G);
  return APIResult(G, result);
}

static PyObject *CmdUnpick(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  API_SETUP_ARGS(G, self, args, "O", &self);
  API_ASSERT(APIEnterNotModal(G));
  EditorInactivate(G);
  APIExit(G);
  return APISuccess();
}

static PyObject *CmdEdit(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  char *str0, *str1, *str2, *str3;
  int pkresi, pkbond;
  int quiet;
  API_SETUP_ARGS(G, self, args, "Ossssiii", &self, &str0, &str1, &str2, &str3,
      &pkresi, &pkbond, &quiet);
  API_ASSERT(APIEnterNotModal(G));
  pymol::Result<> result;
  if(!str0[0]) {
    EditorInactivate(G);
  } else {
    result = EditorSelect(G, str0, str1, str2, str3, pkresi, pkbond, quiet);
  }
  APIExit(G);
  return APIResult(G, result);
}

static PyObject *CmdDrag(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  char *str0;
  OrthoLineType s0 = "";
  int ok = false;
  int quiet;
  int mode;
  ok = PyArg_ParseTuple(args, "Osii", &self, &str0, &quiet,&mode);
  if(ok) {
    API_SETUP_PYMOL_GLOBALS;
    ok = (G != NULL);
  } else {
    API_HANDLE_ERROR;
  }
  if(ok && (ok = APIEnterNotModal(G))) {
    ok = (SelectorGetTmp2(G, str0, s0) >= 0);
    if(ok) {
      ok = ExecutiveSetDrag(G, s0, quiet,mode);
      SelectorFreeTmp(G, s0);
    }
    APIExit(G);
  }
  return APIResultOk(ok);
}

static PyObject *CmdRename(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  char *str1;
  int int1, int2;

  API_SETUP_ARGS(G, self, args, "Osii", &self, &str1, &int1, &int2);
  API_ASSERT(APIEnterNotModal(G));
  auto result = ExecutiveRenameObjectAtoms(G, str1, int1, int2);
  APIExit(G);
  return APIResult(G, result);
}

static PyObject *CmdOrder(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  char *str1;
  int int1, int2;

  API_SETUP_ARGS(G, self, args, "Osii", &self, &str1, &int1, &int2);
  APIEnterNotModal(G);
  auto result = ExecutiveOrder(G, str1, int1, int2);
  APIExit(G);
  return APIResult(G, result);
}

static PyObject *CmdScrollTo(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  char *name;
  int i = 0, r = -1;

  if(!PyArg_ParseTuple(args, "Os|i", &self, &name, &i)) {
    API_HANDLE_ERROR;
    ok_raise(1);
  }

  ok_assert(1, name && name[0]);

  API_SETUP_PYMOL_GLOBALS;
  ok_assert(2, G && APIEnterBlockedNotModal(G));

  r = ExecutiveScrollTo(G, name, i);

ok_except2:
  APIExitBlocked(G);
ok_except1:
  return Py_BuildValue("i", r);
}

static PyObject *CmdWindow(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  int int1, x, y, width, height;
  API_SETUP_ARGS(G, self, args, "Oiiiii", &self, &int1, &x, &y, &width, &height);
  API_ASSERT(G->HaveGUI);
  API_ASSERT(APIEnterNotModal(G));
  {
#ifndef _PYMOL_NO_MAIN
    switch (int1) {
    case 0:
    case 1:
      if(G->Main)
        MainSetWindowVisibility(int1);
      break;
    case 2:                    /* position */
      if(G->Main)
        MainSetWindowPosition(G, x, y);
      break;
    case 3:                    /* size */
      if((width == 0) && (height == 0) && (x != 0) && (y != 0)) {
        width = x;
        height = y;
      }
      if(G->Main)
        MainSetWindowSize(G, width, height);
      break;
    case 4:                    /* position and size */
      if(G->Main) {
        MainSetWindowPosition(G, x, y);
        MainSetWindowSize(G, width, height);
      }
      break;
    case 5:                    /* maximize -- 
                                   should use the window manager, 
                                   but GLUT doesn't provide for that */
      if(G->Main)
        MainMaximizeWindow(G);
      break;
    case 6:
      if(G->Main)
        MainCheckWindowFit(G);
      break;
    }
#endif

    APIExit(G);
  }
  return APISuccess();
}

static PyObject *CmdCEAlign(PyObject *self, PyObject *args)
{
  PyMOLGlobals * G = NULL;
  int ok = false;
  int windowSize = 8, gap_max=30;
  float d0=3.0, d1=4.0;
  PyObject *listA, *listB, *result;
  Py_ssize_t lenA, lenB;

  /* Unpack the arguments from Python */

  ok = PyArg_ParseTuple(args, "OOO|ffii", &self, &listA, &listB, &d0, &d1, &windowSize, &gap_max);

  /* Handle errors */

  if(ok) {
    API_SETUP_PYMOL_GLOBALS;
    ok = (G != NULL);
  } else {
    API_HANDLE_ERROR;
  }

  /* Get the list lengths */

  lenA = PyList_Size(listA);
  if (lenA < 1) {
    result = NULL;
    ok = false;
  }
	
  if(ok)
    lenB = PyList_Size(listB);
  if (ok && lenB < 1) {
    result = NULL;
    ok = false;
  }

  /* Call CEAlign */

  if(ok) {
    APIEnterBlocked(G);
    result = (PyObject*) ExecutiveCEAlign(G, listA, listB, lenA, lenB, d0, d1, windowSize, gap_max);
    APIExitBlocked(G);
  }
  return result;
}

static PyObject *CmdVolume(PyObject *self, PyObject *args)
{ 
  PyMOLGlobals *G = NULL;
  char *volume_name, *map_name, *sele;
  float lvl, fbuf;
  int state = -1;
  float carve;
  int map_state;
  int quiet;
  API_SETUP_ARGS(G, self, args, "Osssffifii", &self, &volume_name, &map_name,
      &sele, &fbuf, &lvl, &state, &carve, &map_state, &quiet);
  API_ASSERT(APIEnterNotModal(G));
  auto res = ExecutiveVolume(G, volume_name, map_name, lvl, sele, fbuf, state,
      carve, map_state, quiet);
  APIExit(G);
  return APIResult(G, res);
}

static PyObject *CmdAssignAtomTypes(PyObject *self, PyObject *args)
{ 
  PyMOLGlobals *G = NULL;
  char *sele;
  int state = -1;
  int ok = false;
  int format;
  int quiet;
  OrthoLineType s1;

  API_SETUP_ARGS(G, self, args, "Osiii", &self, &sele, &format, &state, &quiet);
  {
    APIEnterBlocked(G);
    ok = (SelectorGetTmp(G, sele, s1) >= 0);
    if(ok){
      /* format : 1: mol/sybyl, 2: macromodel/mmd */
      ok = ExecutiveAssignAtomTypes(G, s1, format, state, quiet);
      SelectorFreeTmp(G, s1);
    }
    APIExitBlocked(G);
  }
  return APIResultOk(G, ok);
}

static PyObject *CmdSetDiscrete(PyObject * self, PyObject * args)
{
  const char *name;
  int discrete;
  bool status = false;

  if (!PyArg_ParseTuple(args, "Osi", &self, &name, &discrete)) {
    API_HANDLE_ERROR;
  } else {
    PyMOLGlobals *G = NULL;
    API_SETUP_PYMOL_GLOBALS;

    if (G && APIEnterBlockedNotModal(G)) {
      ObjectMolecule *obj = ExecutiveFindObjectMoleculeByName(G, name);

      if (!obj) {
        PRINTFB(G, FB_Executive, FB_Errors)
          " Executive-Error: object '%s' not found.\n", name ENDFB(G);
      } else {
        status = ObjectMoleculeSetDiscrete(G, obj, discrete);
      }

      APIExitBlocked(G);
    }
  }

  return APIResultOk(status);
}

static PyObject *CmdCountDiscrete(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  char *str1;
  OrthoLineType s1;
  pymol::vla<ObjectMolecule*> list;
  int discrete = 0;

  API_SETUP_ARGS(G, self, args, "Os", &self, &str1);
  APIEnterBlocked(G);
  ok_assert(2, SelectorGetTmp(G, str1, s1) >= 0);

  if((list = ExecutiveGetObjectMoleculeVLA(G, s1))) {
    unsigned int i, size = VLAGetSize(list);
    for(i = 0; i < size; i++)
      if(list[i]->DiscreteFlag)
        discrete++;
  }

  SelectorFreeTmp(G, s1);
ok_except2:
  APIExitBlocked(G);
  return Py_BuildValue("i", discrete);
}

/**
 * Experimental - SUBJECT TO CHANGE
 */
static PyObject *CmdCifGetArray(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  const char *name, *key, *dtype = "";
  ObjectMolecule *obj;
  PyObject *ret = NULL;

  ok_assert(1, PyArg_ParseTuple(args, "Oss|s", &self, &name, &key, &dtype));
  API_SETUP_PYMOL_GLOBALS;
  ok_assert(1, G);
  APIEnterBlocked(G);

  obj = ExecutiveFindObjectMoleculeByName(G, name);

  if (!obj) {
    PRINTFB(G, FB_Executive, FB_Errors)
      " Executive-Error: object '%s' not found.\n", name ENDFB(G);
  } else if (!obj->m_cifdata) {
    PRINTFB(G, FB_Executive, FB_Warnings)
      " Executive-Warning: no cif data for object '%s'\n"
      " ! The 'cif_keepinmemory' setting needs to be set prior to loading a cif file.\n",
      name ENDFB(G);
  } else {
    auto* arr = obj->m_cifdata->get_arr(key);
    if (!arr) {
      PRINTFB(G, FB_Executive, FB_Details)
        " Executive-Details: key '%s' not in cif data for object '%s'.\n", key, name ENDFB(G);
    } else {
      switch (dtype[0]) {
        case 'i': ret = PConvToPyObject(arr->to_vector<int>()); break;
        case 'f': ret = PConvToPyObject(arr->to_vector<double>()); break;
        default:  ret = PConvToPyObject(arr->to_vector<const char*>()); break;
      }
    }
  }

  APIExitBlocked(G);
  return APIAutoNone(ret);
ok_except1:
  API_HANDLE_ERROR;
  return APIAutoNone(NULL);
}

static PyObject* CmdM2ioFirstBlockProperties(PyObject* self, PyObject* args)
{
  const char* contents = nullptr;
  Py_ssize_t length;
  if (!PyArg_ParseTuple(args, "s#", &contents, &length))
    return nullptr;
#ifdef _PYMOL_IP_EXTRAS
#else
  PyErr_SetNone(P_IncentiveOnlyException);
  return nullptr;
#endif
}

static PyMethodDef Cmd_methods[] = {
  {"glViewport", Cmd_glViewport, METH_VARARGS},
  {"_new", Cmd_New, METH_VARARGS},
  {"_start", Cmd_Start, METH_VARARGS},
  {"_stop", Cmd_Stop, METH_VARARGS},
  {"_idle", Cmd_Idle, METH_VARARGS},
  {"_popValidContext", CmdPopValidContext, METH_VARARGS},
  {"_pushValidContext", CmdPushValidContext, METH_VARARGS},
  {"_reshape", Cmd_Reshape, METH_VARARGS},
  {"_getRedisplay", Cmd_GetRedisplay, METH_VARARGS},
  {"_draw", Cmd_Draw, METH_VARARGS},
  {"_button", Cmd_Button, METH_VARARGS},
  {"_drag", Cmd_Drag, METH_VARARGS},
  {"_sdof", Cmd_Sdof, METH_VARARGS},
  {"accept", CmdAccept, METH_VARARGS},
  {"align", CmdAlign, METH_VARARGS},
  {"alter", CmdAlter, METH_VARARGS},
  {"alter_list", CmdAlterList, METH_VARARGS},
  {"alter_state", CmdAlterState, METH_VARARGS},
  {"angle", CmdAngle, METH_VARARGS},
  {"assign_atom_types", CmdAssignAtomTypes, METH_VARARGS},
  {"attach", CmdAttach, METH_VARARGS},
  {"bg_color", CmdBackgroundColor, METH_VARARGS},
  {"bond", CmdBond, METH_VARARGS},
  {"add_bond", CmdAddBond, METH_VARARGS},
  {"rebond", CmdRebond, METH_VARARGS},
  {"busy_draw", CmdBusyDraw, METH_VARARGS},
  {"button", CmdButton, METH_VARARGS},
  /*  {"cache",                 CmdCache,                METH_VARARGS }, */
  {"cartoon", CmdCartoon, METH_VARARGS},
  {"cealign", CmdCEAlign, METH_VARARGS},
  {"center", CmdCenter, METH_VARARGS},
  {"cif_get_array", CmdCifGetArray, METH_VARARGS},
  {"clip", CmdClip, METH_VARARGS},
  {"cls", CmdCls, METH_VARARGS},
  {"color", CmdColor, METH_VARARGS},
  {"colordef", CmdColorDef, METH_VARARGS},
  {"combine_object_ttt", CmdCombineObjectTTT, METH_VARARGS},
  {"coordset_update_thread", CmdCoordSetUpdateThread, METH_VARARGS},
  {"copy", CmdCopy, METH_VARARGS},
  {"create", CmdCreate, METH_VARARGS},
  {"count_states", CmdCountStates, METH_VARARGS},
  {"count_frames", CmdCountFrames, METH_VARARGS},
  {"count_discrete", CmdCountDiscrete, METH_VARARGS},
  {"cycle_valence", CmdCycleValence, METH_VARARGS},
  {"debug", CmdDebug, METH_VARARGS},
  {"decline", CmdDecline, METH_VARARGS},
  {"del_colorection", CmdDelColorection, METH_VARARGS},
  {"fake_drag", CmdFakeDrag, METH_VARARGS},
  {"delete", CmdDelete, METH_VARARGS},
  {"dirty", CmdDirty, METH_VARARGS},
  {"dirty_wizard", CmdDirtyWizard, METH_VARARGS},
  {"dihedral", CmdDihedral, METH_VARARGS},
  {"dist", CmdDist, METH_VARARGS},
  {"do", CmdDo, METH_VARARGS},
  {"draw", CmdDraw, METH_VARARGS},
  {"drag", CmdDrag, METH_VARARGS},
  {"dump", CmdDump, METH_VARARGS},
  {"edit", CmdEdit, METH_VARARGS},
  {"torsion", CmdTorsion, METH_VARARGS},
  {"feedback", CmdFeedback, METH_VARARGS},
  {"find_pairs", CmdFindPairs, METH_VARARGS},
  {"find_molfile_plugin", CmdFindMolfilePlugin, METH_VARARGS},
  {"finish_object", CmdFinishObject, METH_VARARGS},
  {"fit", CmdFit, METH_VARARGS},
  {"fit_pairs", CmdFitPairs, METH_VARARGS},
  {"fix_chemistry", CmdFixChemistry, METH_VARARGS},
  {"flag", CmdFlag, METH_VARARGS},
  {"frame", CmdFrame, METH_VARARGS},
  {"flush_now", CmdFlushNow, METH_VARARGS},
  {"delete_colorection", CmdDelColorection, METH_VARARGS},
  {"dss", CmdAssignSS, METH_VARARGS},
  {"full_screen", CmdFullScreen, METH_VARARGS},
  {"fuse", CmdFuse, METH_VARARGS},
  {"get_angle", CmdGetAngle, METH_VARARGS},
  {"get_area", CmdGetArea, METH_VARARGS},
  {"get_atom_coords", CmdGetAtomCoords, METH_VARARGS},
  {"get_bond_print", CmdGetBondPrint, METH_VARARGS},
  {"get_busy", CmdGetBusy, METH_VARARGS},
  {"get_chains", CmdGetChains, METH_VARARGS},
  {"get_click_string", CmdGetClickString, METH_VARARGS},
  {"get_collada", CmdGetCOLLADA, METH_VARARGS},
  {"get_color", CmdGetColor, METH_VARARGS},
  {"get_colorection", CmdGetColorection, METH_VARARGS},
  {"get_coords", CmdGetCoordsAsNumPy, METH_VARARGS},
  {"get_coordset", CmdGetCoordSetAsNumPy, METH_VARARGS},
  {"get_distance", CmdGetDistance, METH_VARARGS},
  {"get_dihe", CmdGetDihe, METH_VARARGS},
  {"get_drag_object_name", CmdGetDragObjectName, METH_VARARGS},
  {"get_editor_scheme", CmdGetEditorScheme, METH_VARARGS},
  {"get_frame", CmdGetFrame, METH_VARARGS},
  {"get_feedback", CmdGetFeedback, METH_VARARGS},
  {"get_idtf", CmdGetIdtf, METH_VARARGS},
  {"get_legal_name", CmdGetLegalName, METH_VARARGS},
  {"get_m2io_first_block_properties", CmdM2ioFirstBlockProperties, METH_VARARGS},
//  {"get_matrix", CmdGetMatrix, METH_VARARGS},
  {"get_min_max", CmdGetMinMax, METH_VARARGS},
  {"get_mtl_obj", CmdGetMtlObj, METH_VARARGS},
  {"get_model", CmdGetModel, METH_VARARGS},
  {"get_bonds", CmdGetBonds, METH_VARARGS},
  {"get_modal_draw", CmdGetModalDraw, METH_VARARGS},
  {"get_moment", CmdGetMoment, METH_VARARGS},
  {"get_movie_length", CmdGetMovieLength, METH_VARARGS},
  {"get_movie_locked", CmdGetMovieLocked, METH_VARARGS},
  {"get_movie_playing", CmdGetMoviePlaying, METH_VARARGS},
  {"get_names", CmdGetNames, METH_VARARGS},
  {"get_object_color_index", CmdGetObjectColorIndex, METH_VARARGS},
  {"get_object_matrix", CmdGetObjectMatrix, METH_VARARGS},
  {"get_object_ttt", CmdGetObjectTTT, METH_VARARGS},
  {"get_object_settings", CmdGetObjectSettings, METH_VARARGS},
  {"get_origin", CmdGetOrigin, METH_VARARGS},
  {"get_position", CmdGetPosition, METH_VARARGS},
  {"get_povray", CmdGetPovRay, METH_VARARGS},
  {"get_progress", CmdGetProgress, METH_VARARGS},
  {"get_phipsi", CmdGetPhiPsi, METH_VARARGS},
  {"get_renderer", CmdGetRenderer, METH_VARARGS},
  {"get_raw_alignment", CmdGetRawAlignment, METH_VARARGS},
  {"get_seq_align_str", CmdGetSeqAlignStr, METH_VARARGS},
  {"get_session", CmdGetSession, METH_VARARGS},
  {"get_setting_of_type", CmdGetSettingOfType, METH_VARARGS},
  {"get_setting_type", CmdGetSettingType, METH_VARARGS},
  {"get_setting_level", CmdGetSettingLevel, METH_VARARGS},
  {"get_setting_updates", CmdGetSettingUpdates, METH_VARARGS},
  {"get_setting_indices", CmdGetSettingIndices, METH_VARARGS},
  {"get_object_list", CmdGetObjectList, METH_VARARGS},
  {"get_symmetry", CmdGetSymmetry, METH_VARARGS},
  {"get_state", CmdGetState, METH_VARARGS},
  {"get_str", CmdGetStr, METH_VARARGS},
  {"get_title", CmdGetTitle, METH_VARARGS},
  {"get_type", CmdGetType, METH_VARARGS},
  {"get_unused_name", CmdGetUnusedName, METH_VARARGS},
  {"get_version", CmdGetVersion, METH_VARARGS},
  {"get_view", CmdGetView, METH_VARARGS},
  {"get_viewport", CmdGetViewPort, METH_VARARGS},
  {"get_vis", CmdGetVis, METH_VARARGS},
  {"get_capabilities", CmdGetCapabilities, METH_NOARGS, "Get a set of compiled-in capabilities"},
  {"get_ccp4str", CmdGetCCP4Str, METH_VARARGS},
  {"get_volume_field", CmdGetVolumeField, METH_VARARGS},
  {"get_volume_histogram", CmdGetVolumeHistogram, METH_VARARGS},
  {"get_volume_ramp", CmdGetVolumeRamp, METH_VARARGS},
  {"set_volume_ramp", CmdSetVolumeRamp, METH_VARARGS},
  {"get_vrml", CmdGetVRML, METH_VARARGS},
  {"get_wizard", CmdGetWizard, METH_VARARGS},
  {"get_wizard_stack", CmdGetWizardStack, METH_VARARGS},
  {"group", CmdGroup, METH_VARARGS},
  {"h_add", CmdHAdd, METH_VARARGS},
  {"h_fill", CmdHFill, METH_VARARGS},
  {"h_fix", CmdHFix, METH_VARARGS},
  {"identify", CmdIdentify, METH_VARARGS},
  {"index", CmdIndex, METH_VARARGS},
  {"intrafit", CmdIntraFit, METH_VARARGS},
  {"invert", CmdInvert, METH_VARARGS},
  {"interrupt", CmdInterrupt, METH_VARARGS},
  {"isolevel", CmdIsolevel, METH_VARARGS},
  {"isomesh", CmdIsomesh, METH_VARARGS},
  {"isosurface", CmdIsosurface, METH_VARARGS},
  {"wait_deferred", CmdWaitDeferred, METH_VARARGS},
  {"wait_queue", CmdWaitQueue, METH_VARARGS},
  {"label", CmdLabel, METH_VARARGS},
  {"label2", CmdLabel2, METH_VARARGS},
  {"load", CmdLoad, METH_VARARGS},
  {"load_color_table", CmdLoadColorTable, METH_VARARGS},
  {"load_coords", CmdLoadCoords, METH_VARARGS},
  {"load_coordset", CmdLoadCoordSet, METH_VARARGS},
  {"load_png", CmdLoadPNG, METH_VARARGS},
  {"load_object", CmdLoadObject, METH_VARARGS},
  {"load_traj", CmdLoadTraj, METH_VARARGS},
  {"map_generate", CmdMapGenerate, METH_VARARGS},
  {"map_new", CmdMapNew, METH_VARARGS},
  {"map_double", CmdMapDouble, METH_VARARGS},
  {"map_halve", CmdMapHalve, METH_VARARGS},
  {"map_set", CmdMapSet, METH_VARARGS},
  {"map_set_border", CmdMapSetBorder, METH_VARARGS},
  {"map_trim", CmdMapTrim, METH_VARARGS},
  {"mask", CmdMask, METH_VARARGS},
  {"mclear", CmdMClear, METH_VARARGS},
  {"mdo", CmdMDo, METH_VARARGS},
  {"mdump", CmdMDump, METH_VARARGS},
  {"mem", CmdMem, METH_VARARGS},
  {"memory_available", CmdMemoryAvailable, METH_VARARGS},
  {"memory_usage", CmdMemoryUsage, METH_VARARGS},
  {"mmodify", CmdMModify, METH_VARARGS},
  {"move", CmdMove, METH_VARARGS},
  {"mset", CmdMSet, METH_VARARGS},
  {"mplay", CmdMPlay, METH_VARARGS},
  {"mpng_", CmdMPNG, METH_VARARGS},
  {"mmatrix", CmdMMatrix, METH_VARARGS},
  {"mview", CmdMView, METH_VARARGS},
  {"object_update_thread", CmdObjectUpdateThread, METH_VARARGS},
  {"origin", CmdOrigin, METH_VARARGS},
  {"orient", CmdOrient, METH_VARARGS},
  {"onoff", CmdOnOff, METH_VARARGS},
  {"onoff_by_sele", CmdOnOffBySele, METH_VARARGS},
  {"order", CmdOrder, METH_VARARGS},
  {"scrollto", CmdScrollTo, METH_VARARGS},
  {"overlap", CmdOverlap, METH_VARARGS},
  {"paste", CmdPaste, METH_VARARGS},
  {"png", CmdPNG, METH_VARARGS},
  {"pop", CmdPop, METH_VARARGS},
  {"protect", CmdProtect, METH_VARARGS},
  {"pseudoatom", CmdPseudoatom, METH_VARARGS},
#if 1
  {"push_undo", CmdPushUndo, METH_VARARGS},
#endif
  {"pbc_unwrap", CmdPBCUnwrap, METH_VARARGS},
  {"pbc_wrap", CmdPBCWrap, METH_VARARGS},
  {"quit", CmdQuit, METH_VARARGS},
  {"ray_trace_thread", CmdRayTraceThread, METH_VARARGS},
  {"ray_hash_thread", CmdRayHashThread, METH_VARARGS},
  {"ray_anti_thread", CmdRayAntiThread, METH_VARARGS},
  {"ramp_new", CmdRampNew, METH_VARARGS},
  {"ready", CmdReady, METH_VARARGS},
  {"rebuild", CmdRebuild, METH_VARARGS},
  {"recolor", CmdRecolor, METH_VARARGS},
  {"reference", CmdReference, METH_VARARGS},
  {"refresh", CmdRefresh, METH_VARARGS},
  {"refresh_later", CmdRefreshLater, METH_VARARGS},
  {"refresh_now", CmdRefreshNow, METH_VARARGS},
  {"refresh_wizard", CmdRefreshWizard, METH_VARARGS},
  {"remove", CmdRemove, METH_VARARGS},
  {"remove_picked", CmdRemovePicked, METH_VARARGS},
  {"render", CmdRay, METH_VARARGS},
  {"rename", CmdRename, METH_VARARGS},
  {"replace", CmdReplace, METH_VARARGS},
  {"reinitialize", CmdReinitialize, METH_VARARGS},
  {"reset", CmdReset, METH_VARARGS},
  {"reset_rate", CmdResetRate, METH_VARARGS},
  {"reset_matrix", CmdResetMatrix, METH_VARARGS},
  {"revalence", CmdRevalence, METH_VARARGS},
  {"rock", CmdRock, METH_VARARGS},
  {"runpymol", CmdRunPyMOL, METH_VARARGS},
  {"select", CmdSelect, METH_VARARGS},
  {"select_list", CmdSelectList, METH_VARARGS},
  {"set", CmdSet, METH_VARARGS},
  {"set_bond", CmdSetBond, METH_VARARGS},
  {"get_bond", CmdGetBond, METH_VARARGS},
  {"scene", CmdScene, METH_VARARGS},
  {"scene_order", CmdSceneOrder, METH_VARARGS},
  {"get_scene_order", CmdGetSceneOrder, METH_VARARGS},
  {"sculpt_deactivate", CmdSculptDeactivate, METH_VARARGS},
  {"sculpt_activate", CmdSculptActivate, METH_VARARGS},
  {"sculpt_iterate", CmdSculptIterate, METH_VARARGS},
  {"sculpt_purge", CmdSculptPurge, METH_VARARGS},
  {"set_raw_alignment", CmdSetRawAlignment, METH_VARARGS},
  {"set_busy", CmdSetBusy, METH_VARARGS},
  {"set_colorection", CmdSetColorection, METH_VARARGS},
  {"set_dihe", CmdSetDihe, METH_VARARGS},
  {"set_discrete", CmdSetDiscrete, METH_VARARGS},
  {"set_feedback", CmdSetFeedbackMask, METH_VARARGS},
  {"set_frame", CmdSetFrame, METH_VARARGS},
  {"set_name", CmdSetName, METH_VARARGS},
  {"set_geometry", CmdSetGeometry, METH_VARARGS},
//  {"set_matrix", CmdSetMatrix, METH_VARARGS},
  {"set_object_ttt", CmdSetObjectTTT, METH_VARARGS},
  {"set_object_color", CmdSetObjectColor, METH_VARARGS},
  {"set_session", CmdSetSession, METH_VARARGS},
  {"set_state_order", CmdSetStateOrder, METH_VARARGS},
  {"set_symmetry", CmdSetSymmetry, METH_VARARGS},
  {"set_title", CmdSetTitle, METH_VARARGS},
  {"set_wizard", CmdSetWizard, METH_VARARGS},
  {"set_wizard_stack", CmdSetWizardStack, METH_VARARGS},
  {"set_view", CmdSetView, METH_VARARGS},
  {"set_vis", CmdSetVis, METH_VARARGS},
  {"showhide", CmdShowHide, METH_VARARGS},
  {"slice_new", CmdSliceNew, METH_VARARGS},
  {"smooth", CmdSmooth, METH_VARARGS},
  {"sort", CmdSort, METH_VARARGS},
  {"spectrum", CmdSpectrum, METH_VARARGS},
  {"spheroid", CmdSpheroid, METH_VARARGS},
  {"splash", CmdSplash, METH_VARARGS},
  {"stereo", CmdStereo, METH_VARARGS},
  {"system", CmdSystem, METH_VARARGS},
  {"symexp", CmdSymExp, METH_VARARGS},
  {"symmetry_copy", CmdSymmetryCopy, METH_VARARGS},
  {"test", CmdTest, METH_VARARGS},
  {"test2", CmdTest2, METH_VARARGS},
  {"toggle", CmdToggle, METH_VARARGS},
  {"matrix_copy", CmdMatrixCopy, METH_VARARGS},
  {"transform_object", CmdTransformObject, METH_VARARGS},
  {"transform_selection", CmdTransformSelection, METH_VARARGS},
  {"translate_atom", CmdTranslateAtom, METH_VARARGS},
  {"translate_object_ttt", CmdTranslateObjectTTT, METH_VARARGS},
  {"turn", CmdTurn, METH_VARARGS},
  {"viewport", CmdViewport, METH_VARARGS},
  {"vdw_fit", CmdVdwFit, METH_VARARGS},
  {"volume", CmdVolume, METH_VARARGS},
  {"undo", CmdUndo, METH_VARARGS},
  {"unpick", CmdUnpick, METH_VARARGS},
  {"unset", CmdUnset, METH_VARARGS},
  {"unset_bond", CmdUnsetBond, METH_VARARGS},
  {"update", CmdUpdate, METH_VARARGS},
  {"window", CmdWindow, METH_VARARGS},
  {"zoom", CmdZoom, METH_VARARGS},
  {NULL, NULL}                  /* sentinel */
};

#ifdef __cplusplus
extern "C" {
#endif

PyObject * PyInit__cmd(void)
{
  static struct PyModuleDef moduledef = {
    PyModuleDef_HEAD_INIT,
    "pymol._cmd",
    "DO NOT USE",
    -1,
    Cmd_methods };
  return PyModule_Create(&moduledef);
}

// Required for AxPyMOL
void init_cmd(void)
{
  auto _cmd = PyInit__cmd();
  if (_cmd) {
    PyDict_SetItemString(PyImport_GetModuleDict(), "pymol._cmd", _cmd);
    Py_DECREF(_cmd);
  }
}

#ifdef __cplusplus
}
#endif

#else
typedef int this_file_is_no_longer_empty;
#endif
