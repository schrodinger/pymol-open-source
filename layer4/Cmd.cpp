
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
#include"os_python.h"
#include"PyMOLGlobals.h"
#include"PyMOLOptions.h"
#include"os_predef.h"
#include"os_gl.h"
#include"os_std.h"
#include"Version.h"
#include"MemoryDebug.h"
#include"Err.h"
#include"Util.h"
#include"Cmd.h"
#include"ButMode.h"
#include"Ortho.h"
#include"ObjectMolecule.h"
#include"ObjectMesh.h"
#include"ObjectMap.h"
#include"ObjectCallback.h"
#include"ObjectCGO.h"
#include"ObjectSurface.h"
#include"ObjectSlice.h"
#include"Executive.h"
#include"Selector.h"
#include"main.h"
#include"Scene.h"
#include"SceneRay.h"
#include"Setting.h"
#include"Movie.h"
#include"Export.h"
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
#include"MacPyMOL.h"
#include"ObjectAlignment.h"

#include "MovieScene.h"
#include "CifFile.h"

#include "MoleculeExporter.h"

#define tmpSele "_tmp"
#define tmpSele1 "_tmp1"
#define tmpSele2 "_tmp2"

static int flush_count = 0;

#ifndef _PYMOL_NO_MAIN
#ifndef _PYMOL_WX_GLUT
static int run_only_once = true;
#endif
#endif
#ifdef _PYMOL_WX_GLUT
#ifndef _PYMOL_OLD_ACTIVEX
#ifndef _PYMOL_EMBEDDED
static int run_only_once = true;
#endif
#endif
#endif

#define API_SETUP_PYMOL_GLOBALS \
  G = _api_get_pymol_globals(self)

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

/*
 * Start a headless singleton instance in the current thread.
 *
 * Unlike when calling `pymol.finish_launching()`, there is no event loop,
 * so animations, continuous sculpting and modal draw are not supported.
 *
 * After calling this, SingletonPyMOLGlobals will be available.
 */
static void launch_library_singleton() {
  PyRun_SimpleString(
      "print(' PyMOL not running, entering library mode (experimental)')\n"
      "import pymol.invocation, pymol2\n"
      "pymol.invocation.parse_args(['pymol', '-cqk'])\n"
      "pymol2.SingletonPyMOL().start()");
}

/*
 * Get the PyMOLGlobals pointer from the `self` object (_self._COb in Python).
 *
 * If _COb is None, launch a headless singleton ("library mode").
 */
static PyMOLGlobals * _api_get_pymol_globals(PyObject * self) {
  if(self == Py_None) {
    launch_library_singleton();
    return SingletonPyMOLGlobals;
  }

  if(self && PyCObject_Check(self)) { \
    PyMOLGlobals **G_handle = (PyMOLGlobals**)PyCObject_AsVoidPtr(self); \
    if(G_handle) { \
      return *G_handle;
    } \
  }

  return NULL;
}

#define API_HANDLE_ERROR \
   if (PyErr_Occurred()) PyErr_Print(); \
   fprintf(stderr,"API-Error: in %s line %d.\n",__FILE__,__LINE__);


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

static PyObject *APISuccess(void)
{                               /* success returns None */
  return PConvAutoNone(Py_None);
}

static PyObject *APIFailure(void)
{                               /* returns -1: a general unspecified
                                 * error */
  return (Py_BuildValue("i", -1));
}

static PyObject *APIResultCode(int code)
{                               /* innteger result code
                                 * (could be a value, a
                                 * count, or a boolean) */
  return (Py_BuildValue("i", code));
}

static PyObject *APIResultOk(int ok)
{
  if(ok)
    return APISuccess();
  else
    return APIFailure();
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
  OrthoLineType s1;
  char *name, *resn, *resi, *chain, *segi, *elem;
  float vdw;
  int hetatm, color;
  float b, q;
  PyObject *pos;
  int state, mode, quiet;
  int ok = false;

  ok = PyArg_ParseTuple(args, "OssssssssfiffsOiiii", &self,
                        &object_name, &sele, &name, &resn, &resi, &chain,
                        &segi, &elem, &vdw, &hetatm, &b, &q, &label, &pos, &color,
                        &state, &mode, &quiet);
  if(ok) {
    API_SETUP_PYMOL_GLOBALS;
    ok = (G != NULL);
  } else {
    API_HANDLE_ERROR;
  }
  if(ok) {
    float pos_array[3], *pos_ptr = NULL;
    if(ok) {
      if(pos && PyTuple_Check(pos) && (PyTuple_Size(pos) == 3))
        if(PyArg_ParseTuple(pos, "fff", pos_array, pos_array + 1, pos_array + 2))
          pos_ptr = pos_array;
    }
    if((ok = APIEnterBlockedNotModal(G))) {
      if(sele[0])
        ok = (SelectorGetTmp2(G, sele, s1) >= 0);
      else
        s1[0] = 0;
      if(ok) {
        ok = ExecutivePseudoatom(G, object_name, s1,
                                 name, resn, resi, chain, segi, elem,
                                 vdw, hetatm, b, q, label, pos_ptr,
                                 color, state, mode, quiet);
      }
      if(sele[0])
        SelectorFreeTmp(G, s1);
      APIExitBlocked(G);
    }
  }
  return APIResultOk(ok);
}

static PyObject *CmdFixChemistry(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  char *str2, *str3;
  int ok = false;
  int quiet;
  int invalidate;
  ok = PyArg_ParseTuple(args, "Ossii", &self, &str2, &str3, &invalidate, &quiet);
  if(ok) {
    API_SETUP_PYMOL_GLOBALS;
    ok = (G != NULL);
  } else {
    API_HANDLE_ERROR;
  }
  if(ok && (ok = APIEnterNotModal(G))) {
    if(ok)
      ok = ExecutiveFixChemistry(G, str2, str3, invalidate, quiet);
    APIExit(G);
  }
  return APIResultOk(ok);
}

static PyObject *CmdRayAntiThread(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  int ok = false;
  PyObject *py_thread_info;

  CRayAntiThreadInfo *thread_info = NULL;

  ok = PyArg_ParseTuple(args, "OO", &self, &py_thread_info);
  if(ok)
    ok = PyCObject_Check(py_thread_info);
  if(ok)
    ok = ((thread_info = (CRayAntiThreadInfo *) PyCObject_AsVoidPtr(py_thread_info)) != NULL);
  if(ok) {
    API_SETUP_PYMOL_GLOBALS;
    ok = (G != NULL);
  } else {
    API_HANDLE_ERROR;
  }
  if(ok) {
    PUnblock(G);
    RayAntiThread(thread_info);
    PBlock(G);
  }
  return APIResultOk(ok);
}

static PyObject *CmdRayHashThread(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  int ok = false;
  PyObject *py_thread_info;

  CRayHashThreadInfo *thread_info = NULL;

  ok = PyArg_ParseTuple(args, "OO", &self, &py_thread_info);
  if(ok)
    ok = PyCObject_Check(py_thread_info);
  if(ok)
    ok = ((thread_info = (CRayHashThreadInfo*) PyCObject_AsVoidPtr(py_thread_info)) != NULL);
  if(ok) {
    API_SETUP_PYMOL_GLOBALS;
    ok = (G != NULL);
  }
  if(ok) {
    PUnblock(G);
    RayHashThread(thread_info);
    PBlock(G);
  }
  return APIResultOk(ok);
}

static PyObject *CmdRayTraceThread(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  int ok = false;
  PyObject *py_thread_info;

  CRayThreadInfo *thread_info = NULL;

  ok = PyArg_ParseTuple(args, "OO", &self, &py_thread_info);
  if(ok)
    ok = PyCObject_Check(py_thread_info);
  if(ok)
    ok = ((thread_info = (CRayThreadInfo*) PyCObject_AsVoidPtr(py_thread_info)) != NULL);
  if(ok) {
    API_SETUP_PYMOL_GLOBALS;
    ok = (G != NULL);
  } else {
    API_HANDLE_ERROR;
  }
  if(ok) {
    PUnblock(G);
    RayTraceThread(thread_info);
    PBlock(G);
  }
  return APIResultOk(ok);
}

static PyObject *CmdCoordSetUpdateThread(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  int ok = false;
  PyObject *py_thread_info;

  CCoordSetUpdateThreadInfo *thread_info = NULL;

  ok = PyArg_ParseTuple(args, "OO", &self, &py_thread_info);
  if(ok)
    ok = PyCObject_Check(py_thread_info);
  if(ok)
    ok = ((thread_info = (CCoordSetUpdateThreadInfo*) PyCObject_AsVoidPtr(py_thread_info)) != NULL);
  if(ok) {
    API_SETUP_PYMOL_GLOBALS;
    ok = (G != NULL);
  }
  if(ok) {
    PUnblock(G);
    CoordSetUpdateThread(thread_info);
    PBlock(G);
  }
  return APIResultOk(ok);
}

static PyObject *CmdObjectUpdateThread(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  int ok = false;
  PyObject *py_thread_info;

  CObjectUpdateThreadInfo *thread_info = NULL;

  ok = PyArg_ParseTuple(args, "OO", &self, &py_thread_info);
  if(ok)
    ok = PyCObject_Check(py_thread_info);
  if(ok)
    ok = ((thread_info = (CObjectUpdateThreadInfo*) PyCObject_AsVoidPtr(py_thread_info)) != NULL);
  if(ok) {
    API_SETUP_PYMOL_GLOBALS;
    ok = (G != NULL);
  } else {
    API_HANDLE_ERROR;
  }
  if(ok) {
    PUnblock(G);
    SceneObjectUpdateThread(thread_info);
    PBlock(G);
  }
  return APIResultOk(ok);
}

static PyObject *CmdGetMovieLocked(PyObject * self, PyObject * args)
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
    return APIResultCode(MovieLocked(G));
  } else {
    return APIResultOk(ok);
  }
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
  ok = PyArg_ParseTuple(args, "OOs", &self, &list, &prefix);
  if(ok) {
    API_SETUP_PYMOL_GLOBALS;
    ok = (G != NULL);
  } else {
    API_HANDLE_ERROR;
  }
  if(ok && (ok = APIEnterBlockedNotModal(G))) {
    ok = SelectorColorectionFree(G, list, prefix);
    APIExitBlocked(G);
  }
  return APIResultOk(ok);
}

static PyObject *CmdSetColorectionName(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  int ok = false;
  PyObject *list;
  char *prefix, *new_prefix;
  ok = PyArg_ParseTuple(args, "OOss", &self, &list, &prefix, &new_prefix);
  if(ok) {
    API_SETUP_PYMOL_GLOBALS;
    ok = (G != NULL);
  } else {
    API_HANDLE_ERROR;
  }
  if(ok && (ok = APIEnterBlockedNotModal(G))) {
    ok = SelectorColorectionSetName(G, list, prefix, new_prefix);
    APIExitBlocked(G);
  }
  return APIResultOk(ok);
}

static PyObject *CmdSetColorection(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  int ok = false;
  char *prefix;
  PyObject *list;
  ok = PyArg_ParseTuple(args, "OOs", &self, &list, &prefix);
  if(ok) {
    API_SETUP_PYMOL_GLOBALS;
    ok = (G != NULL);
  } else {
    API_HANDLE_ERROR;
  }
  if(ok && (ok = APIEnterBlockedNotModal(G))) {
    ok = SelectorColorectionApply(G, list, prefix);
    APIExitBlocked(G);
  }
  return APIResultOk(ok);
}

static PyObject *CmdGetColorection(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  PyObject *result = NULL;
  int ok = false;
  char *prefix;
  ok = PyArg_ParseTuple(args, "Os", &self, &prefix);
  if(ok) {
    API_SETUP_PYMOL_GLOBALS;
    ok = (G != NULL);
  } else {
    API_HANDLE_ERROR;
  }
  if(ok && (ok = APIEnterBlockedNotModal(G))) {
    result = SelectorColorectionGet(G, prefix);
    APIExitBlocked(G);
  }
  return (APIAutoNone(result));
}

static PyObject *CmdSetRawAlignment(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = nullptr;
  const char *alnname;
  const char *guidename;
  PyObject *raw;
  int state = 0, quiet = 1;
  ObjectMolecule *guide = nullptr;

  if(!PyArg_ParseTuple(args, "sOsii" "O",
        &alnname, &raw, &guidename, &state, &quiet,
        &self)) {
    return nullptr;
  }

  API_SETUP_PYMOL_GLOBALS;
  if(G == nullptr) {
    PyErr_BadInternalCall();
    return nullptr;
  }

  if (guidename[0]) {
    guide = ExecutiveFindObjectMoleculeByName(G, guidename);
  }

  if(!PyList_Check(raw)) {
    PyErr_SetString(PyExc_TypeError, "alignment must be list");
    return nullptr;
  }

  auto n_cols = PyList_Size(raw);

  pymol::vla<int> align_vla(n_cols * 3);
  size_t vla_offset = 0;

  for(size_t c = 0; c < n_cols; ++c) {
    PyObject * col = PyList_GetItem(raw, c);

    if(!PyList_Check(col)) {
      PyErr_SetString(PyExc_TypeError, "columns must be list");
      return nullptr;
    }

    auto n_idx = PyList_Size(col);

    for(size_t i = 0; i < n_idx; ++i) {
      const char * model;
      int index;

      PyObject * idx = PyList_GetItem(col, i);

      if(!PyArg_ParseTuple(idx, "si", &model, &index)) {
        PyErr_SetString(PyExc_TypeError, "indices must be (str, int)");
        return nullptr;
      }

      ObjectMolecule * mol = ExecutiveFindObjectMoleculeByName(G, model);

      if(!mol) {
        PyErr_Format(PyExc_KeyError, "object '%s' not found", model);
        return nullptr;
      }

      if (!guide) {
        guide = mol;
      }

      if (index < 1 || mol->NAtom < index) {
        PyErr_Format(PyExc_IndexError, "index ('%s', %d) out of range", model, index);
        return nullptr;
      }

      auto uid = AtomInfoCheckUniqueID(G, mol->AtomInfo + index - 1);
      *(align_vla.check(vla_offset++)) = uid;
    }

    *(align_vla.check(vla_offset++)) = 0;
  }

  align_vla.resize(vla_offset);

  // does alignment object already exist?
  auto cobj = ExecutiveFindObjectByName(G, alnname);
  if (cobj && cobj->type != cObjectAlignment) {
    ExecutiveDelete(G, cobj->Name);
    cobj = nullptr;
  }

  // create alignment object
  cobj = (CObject*) ObjectAlignmentDefine(G, (ObjectAlignment*) cobj,
      align_vla.data(), state, true, guide, nullptr);

  // manage alignment object
  ObjectSetName(cobj, alnname);
  ExecutiveManageObject(G, cobj, 0, quiet);
  SceneInvalidate(G);

  // make available as selection FIXME find better solution
  cobj->update();

  return APISuccess();
}

static PyObject* GetRawAlignment(PyMOLGlobals* G,
    const ObjectAlignment* alnobj,
    bool active_only,
    int state)
{
  if (state >= alnobj->NState) {
    PyErr_Format(PyExc_IndexError, "state %d >= NState %d", state, alnobj->NState);
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
          && (!active_only || eoo->obj->Obj.Enabled)
          && (!hide_underscore || eoo->obj->Obj.Name[0] != '_')) {
        PyObject * idx = Py_BuildValue("si", eoo->obj->Obj.Name, eoo->atm + 1);
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
  int ok = false;
  const char *name;
  int active_only;
  int state = 0;
  PyObject *result = NULL;
  ok = PyArg_ParseTuple(args, "Osi|i", &self, &name, &active_only, &state);
  if(ok) {
    API_SETUP_PYMOL_GLOBALS;
    ok = (G != NULL);
  } else {
    API_HANDLE_ERROR;
  }
  if(ok && (ok = APIEnterBlockedNotModal(G))) {
    if (!name[0]) {
      name = ExecutiveGetActiveAlignment(G);
    }
    if (name && name[0]) {
      CObject *obj = ExecutiveFindObjectByName(G, name);
      if (obj && obj->type == cObjectAlignment) {
        result = GetRawAlignment(G, (ObjectAlignment*) obj, active_only, state);
      } else {
        PyErr_Format(PyExc_KeyError, "no such alignment: '%s'", name);
      }
    }
    APIExitBlocked(G);
  }
  if(!result && !PyErr_Occurred()) {
    return APIFailure();
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
      CObject *obj = ExecutiveFindObjectByName(G, object);
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
    if (G && APIEnterNotModal(G)) {
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
  if (!PyArg_ParseTuple(args, "Osii", &self, &name, &state, &quiet)) {
    API_HANDLE_ERROR;
  } else {
    API_SETUP_PYMOL_GLOBALS;
    if (G && APIEnterNotModal(G)) {
      auto v = ObjectMapGetCCP4Str(G, name, state, quiet);
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
  PyObject *result = NULL;
  int ok = false;
  char* objName;
  float min_val = 0.f, max_val = 0.f;
  int n_points = 64;
  ok = PyArg_ParseTuple(args, "Os|i(ff)", &self, &objName, &n_points, &min_val, &max_val);

  if(ok) {
    API_SETUP_PYMOL_GLOBALS;
    ok = (G != NULL);
  } else {
    API_HANDLE_ERROR;
  }
  if(ok && (ok = APIEnterBlockedNotModal(G))) {
    float * hist = ExecutiveGetHistogram(G, objName, n_points, min_val, max_val);
    if (hist) {
      result = PConvFloatArrayToPyList(hist, n_points + 4);
      mfree(hist);
    }
    APIExitBlocked(G);
  }

  if(!result) {
    return APIFailure();
  } else
    return result;
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
  int ok = false;
  char* objName;
  PyObject *ramp_list;
  float *float_array;
  int list_len;

  if(!PyArg_ParseTuple(args, "OsO", &self, &objName, &ramp_list)) {
    API_HANDLE_ERROR;
    ok_raise(1);
  }

  API_SETUP_PYMOL_GLOBALS;
  ok_assert(1, G && APIEnterBlockedNotModal(G));

  ok_assert(2, PyList_Check(ramp_list));
  ok_assert(2, list_len = PyList_Size(ramp_list));
  ok_assert(2, PConvPyListToFloatArray(ramp_list, &float_array));

  ok = ExecutiveSetVolumeRamp(G, objName, float_array, list_len);

  if(!ok)
    mfree(float_array);

ok_except2:
  APIExitBlocked(G);
ok_except1:
  return APIResultOk(ok);
}

static PyObject *CmdGetVis(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  PyObject *result = NULL;
  int ok = false;
  ok = PyArg_ParseTuple(args, "O", &self);
  if(ok) {
    API_SETUP_PYMOL_GLOBALS;
    ok = (G != NULL);
  } else {
    API_HANDLE_ERROR;
  }
  if(ok && (ok = APIEnterBlockedNotModal(G))) {
    result = ExecutiveGetVisAsPyDict(G);
    APIExitBlocked(G);
  }
  return (APIAutoNone(result));
}

static PyObject *CmdSetVis(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  int ok = false;
  PyObject *visDict;
  ok = PyArg_ParseTuple(args, "OO", &self, &visDict);
  if(ok) {
    API_SETUP_PYMOL_GLOBALS;
    ok = (G != NULL);
  } else {
    API_HANDLE_ERROR;
  }
  if(ok && (ok = APIEnterBlockedNotModal(G))) {
    ok = ExecutiveSetVisFromPyDict(G, visDict);
    APIExitBlocked(G);
  }
  return APIResultOk(ok);
}

static PyObject *CmdReinitialize(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  int ok = false;
  int what;
  char *object;
  ok = PyArg_ParseTuple(args, "Ois", &self, &what, &object);
  if(ok) {
    API_SETUP_PYMOL_GLOBALS;
    ok = (G != NULL);
  } else {
    API_HANDLE_ERROR;
  }
  if(ok && (ok = APIEnterNotModal(G))) {
    ok = ExecutiveReinitialize(G, what, object);
    APIExit(G);
  }
  return APIResultOk(ok);

}

static PyObject *CmdSpectrum(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  char *str1, *expr, *prefix;
  float min, max;
  int digits, start, stop, byres;
  int quiet;
  int ok = false;
  float min_ret, max_ret;
  PyObject *result = Py_None;
  ok = PyArg_ParseTuple(args, "Ossffiisiii", &self, &str1, &expr,
                        &min, &max, &start, &stop, &prefix, &digits, &byres, &quiet);
  if(ok) {
    API_SETUP_PYMOL_GLOBALS;
    ok = (G != NULL);
  } else {
    API_HANDLE_ERROR;
  }
  if(ok && (ok = APIEnterNotModal(G))) {
    if(ok) {
      ok =
        ExecutiveSpectrum(G, str1, expr, min, max, start, stop, prefix, digits, byres,
                          quiet, &min_ret, &max_ret);
    }
    APIExit(G);
    if(ok) {
      result = Py_BuildValue("ff", min_ret, max_ret);
    }
  }
  return (APIAutoNone(result));
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
  if(ok && (ok = APIEnterNotModal(G))) {
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
  int ok = false;
  char *str1, *str2;
  int state;
  OrthoLineType s1;
  float a, b, c, alpha, beta, gamma;

  ok = PyArg_ParseTuple(args, "Osiffffffs", &self, &str1, &state, &a, &b, &c,
                        &alpha, &beta, &gamma, &str2);
  if(ok) {
    API_SETUP_PYMOL_GLOBALS;
    ok = (G != NULL);
  } else {
    API_HANDLE_ERROR;
  }
  if(ok && (ok = APIEnterNotModal(G))) {
    ok = (SelectorGetTmp2(G, str1, s1) >= 0);
    if(ok)
      ok = ExecutiveSetSymmetry(G, s1, state, a, b, c, alpha, beta, gamma, str2);
    SelectorFreeTmp(G, s1);
    APIExit(G);
  }
  return APIResultOk(ok);
}

static PyObject *CmdGetSymmetry(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  int ok = false;
  char *str1;
  OrthoLineType s1;
  float a, b, c, alpha, beta, gamma;
  int state;
  WordType sg;
  PyObject *result = NULL;
  int defined;
  ok = PyArg_ParseTuple(args, "Osi", &self, &str1, &state);
  if(ok) {
    API_SETUP_PYMOL_GLOBALS;
    ok = (G != NULL);
  } else {
    API_HANDLE_ERROR;
  }
  if(ok && (ok = APIEnterNotModal(G))) {
    ok = (SelectorGetTmp2(G, str1, s1) >= 0);
    if(ok)
      ok = ExecutiveGetSymmetry(G, s1, state, &a, &b, &c, &alpha, &beta, &gamma, sg, &defined);
    APIExit(G);
    if(ok) {
      if(defined) {
        result = PyList_New(7);
        if(result) {
          PyList_SetItem(result, 0, PyFloat_FromDouble(a));
          PyList_SetItem(result, 1, PyFloat_FromDouble(b));
          PyList_SetItem(result, 2, PyFloat_FromDouble(c));
          PyList_SetItem(result, 3, PyFloat_FromDouble(alpha));
          PyList_SetItem(result, 4, PyFloat_FromDouble(beta));
          PyList_SetItem(result, 5, PyFloat_FromDouble(gamma));
          PyList_SetItem(result, 6, PyString_FromString(sg));
        }
      } else {                  /* no symmetry defined, then return empty list */
        result = PyList_New(0);
      }
    }
    SelectorFreeTmp(G, s1);
  }
  return (APIAutoNone(result));
}

static PyObject *CmdSmooth(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  int ok = false;
  char *str1;
  int int1, int2, int3, int4, int5, int6;
  ok =
    PyArg_ParseTuple(args, "Osiiiiii", &self, &str1, &int1, &int2, &int3, &int4, &int5,
                     &int6);
  if(ok) {
    API_SETUP_PYMOL_GLOBALS;
    ok = (G != NULL);
  } else {
    API_HANDLE_ERROR;
  }
  if(ok && (ok = APIEnterNotModal(G))) {
    if(ok)
      ok = ExecutiveSmooth(G, str1, int1, int2, int3, int4, int5, int6);
    APIExit(G);
  }
  return APIResultOk(ok);
}

static PyObject *CmdGetSession(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  int ok = false;
  PyObject *dict;
  int partial, quiet;
  char *names;
  ok = PyArg_ParseTuple(args, "OOsii", &self, &dict, &names, &partial, &quiet);
  if(ok) {
    API_SETUP_PYMOL_GLOBALS;
    ok = (G != NULL);
  } else {
    API_HANDLE_ERROR;
  }
  if(ok && (ok = APIEnterBlockedNotModal(G))) {
    ok = ExecutiveGetSession(G, dict, names, partial, quiet);
    APIExitBlocked(G);
  }
  return APIResultOk(ok);
}

static PyObject *CmdSetSession(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  int ok = false;
  int quiet, partial;
  PyObject *obj;

  ok = PyArg_ParseTuple(args, "OOii", &self, &obj, &partial, &quiet);
  if(ok) {
    API_SETUP_PYMOL_GLOBALS;
    ok = (G != NULL);
  } else {
    API_HANDLE_ERROR;
  }
  if(ok && (ok = APIEnterBlockedNotModal(G))) {
    ok = ExecutiveSetSession(G, obj, partial, quiet);
    APIExitBlocked(G);
  }
  return APIResultOk(ok);
}

static PyObject *CmdSetName(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  int ok = false;
  char *str1, *str2;
  ok = PyArg_ParseTuple(args, "Oss", &self, &str1, &str2);
  if(ok) {
    API_SETUP_PYMOL_GLOBALS;
    ok = (G != NULL);
  } else {
    API_HANDLE_ERROR;
  }
  if(ok && (ok = APIEnterNotModal(G))) {
    ok = ExecutiveSetName(G, str1, str2);
    APIExit(G);
  }
  return APIResultOk(ok);
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

static PyObject *CmdPGlutGetRedisplay(PyObject * self, PyObject * args)
{
#ifdef _PYMOL_PRETEND_GLUT
#ifndef _PYMOL_NO_GLUT
  return (APIResultCode(p_glutGetRedisplay()));
#else
  return (APIResultCode(0));
#endif
#else
  return (APIResultCode(0));
#endif
}

static PyObject *CmdPGlutEvent(PyObject * self, PyObject * args)
{
  int ok = false;
#ifdef _PYMOL_PRETEND_GLUT
#ifndef _PYMOL_NO_GLUT
  PyMOLGlobals *G = NULL;
  p_glut_event ev;
  ok = PyArg_ParseTuple(args, "Oiiiiii", &self, &ev.event_code,
                        &ev.x, &ev.y, &ev.input, &ev.state, &ev.mod);
  if(ok) {
    API_SETUP_PYMOL_GLOBALS;
    ok = (G != NULL);
  } else {
    API_HANDLE_ERROR;
  }
  if(ok) {
    PUnblock(G);
    p_glutHandleEvent(&ev);
    PBlock(G);
  }
#endif
#endif
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
  PyMOLGlobals *G = NULL;
  int ok = false;

  const char *key, *action, *message = NULL, *new_key = NULL;
  unsigned char store_view = true,
                store_color = true,
                store_active = true,
                store_rep = true,
                store_frame = true,
                hand = true;
  float animate = -1.0;
  const char * sele = "all";

  if(!PyArg_ParseTuple(args, "Oss|zbbbbbfzbs", &self, &key, &action, &message,
        &store_view, &store_color, &store_active, &store_rep, &store_frame,
        &animate, &new_key, &hand, &sele)) {
    API_HANDLE_ERROR;
    ok_raise(2);
  }

  API_SETUP_PYMOL_GLOBALS;
  ok_assert(2, G && APIEnterBlockedNotModal(G));

  ok = MovieSceneFunc(G, key, action, message,
      store_view, store_color, store_active, store_rep, store_frame,
      animate, new_key, hand, sele);

  APIExitBlocked(G);
ok_except2:
  return APIResultOk(ok);
}

static PyObject *CmdSceneOrder(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  int ok = false;

  const char *names, *location;
  unsigned char sort;

  if(!PyArg_ParseTuple(args, "Osbs", &self, &names, &sort, &location)) {
    API_HANDLE_ERROR;
    ok_raise(2);
  }

  API_SETUP_PYMOL_GLOBALS;
  ok_assert(2, G && APIEnterBlockedNotModal(G));

  ok = MovieSceneOrder(G, names, sort, location);

  APIExitBlocked(G);
ok_except2:
  return APIResultOk(ok);
}

static PyObject *CmdGetSceneOrder(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  PyObject * result = NULL;

  if(!PyArg_ParseTuple(args, "O", &self)) {
    API_HANDLE_ERROR;
    ok_raise(2);
  }

  API_SETUP_PYMOL_GLOBALS;
  ok_assert(2, G && APIEnterBlockedNotModal(G));

  result = PConvToPyObject(MovieSceneGetOrder(G));

  APIExitBlocked(G);
ok_except2:
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
  int ok = PyArg_ParseTuple(args, "Os(ffffffffffffffff)ii",
                            &self, &name,
                            &ttt[0], &ttt[1], &ttt[2], &ttt[3],
                            &ttt[4], &ttt[5], &ttt[6], &ttt[7],
                            &ttt[8], &ttt[9], &ttt[10], &ttt[11],
                            &ttt[12], &ttt[13], &ttt[14], &ttt[15],
                            &state, &quiet);
  if(ok) {
    API_SETUP_PYMOL_GLOBALS;
    ok = (G != NULL);
  } else {
    API_HANDLE_ERROR;
  }
  if(ok && (ok = APIEnterNotModal(G))) {
    ExecutiveSetObjectTTT(G, name, ttt, state, quiet, SettingGetGlobal_i(G, cSetting_movie_auto_store));
    APIExit(G);
  }
  return APIResultOk(ok);
}

static PyObject *CmdTranslateObjectTTT(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  float mov[3];
  char *name;
  int ok = PyArg_ParseTuple(args, "Os(fff)",
                            &self, &name,
                            &mov[0], &mov[1], &mov[2]);
  if(ok) {
    API_SETUP_PYMOL_GLOBALS;
    ok = (G != NULL);
  } else {
    API_HANDLE_ERROR;
  }
  if(ok && (ok = APIEnterNotModal(G))) {
    ok = ExecutiveTranslateObjectTTT(G, name, mov, SettingGetGlobal_i(G, cSetting_movie_auto_store), true);
    APIExit(G);
  }
  return APIResultOk(ok);
}

static PyObject *CmdCombineObjectTTT(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  char *name;
  PyObject *m;
  float ttt[16];
  int ok = false;
  ok = PyArg_ParseTuple(args, "OsO", &self, &name, &m);
  if(ok) {
    API_SETUP_PYMOL_GLOBALS;
    ok = (G != NULL);
  } else {
    API_HANDLE_ERROR;
  }
  if(ok) {
    if(PConvPyListToFloatArrayInPlace(m, ttt, 16) > 0) {
      if((ok = APIEnterNotModal(G))) {
        ok = ExecutiveCombineObjectTTT(G, name, ttt, false, -1);
        APIExit(G);
      }
    } else {
      PRINTFB(G, FB_CCmd, FB_Errors)
        "CmdCombineObjectTTT-Error: bad matrix\n" ENDFB(G);
      ok = false;
    }
  }
  return APIResultOk(ok);
}

static PyObject *CmdGetColor(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  char *name;
  int mode;
  int ok = false;
  int a, nc, nvc;
  const float *rgb;
  int index;
  PyObject *result = NULL;
  PyObject *tup;
  ok = PyArg_ParseTuple(args, "Osi", &self, &name, &mode);
  if(ok) {
    API_SETUP_PYMOL_GLOBALS;
    ok = (G != NULL);
  } else {
    API_HANDLE_ERROR;
  }
  if(ok && (ok = APIEnterBlockedNotModal(G))) {
    switch (mode) {
    case 0:                    /* by name or index, return floats */
      index = ColorGetIndex(G, name);
      if(index >= 0) {
        rgb = ColorGet(G, index);
        tup = PyTuple_New(3);
        PyTuple_SetItem(tup, 0, PyFloat_FromDouble(*(rgb++)));
        PyTuple_SetItem(tup, 1, PyFloat_FromDouble(*(rgb++)));
        PyTuple_SetItem(tup, 2, PyFloat_FromDouble(*rgb));
        result = tup;
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
          WordType buffer;
          tup = PyTuple_New(2);
          {
            const char *color_name = ColorGetName(G, a);
            if(color_name) {
              strcpy(buffer, color_name);
            } else {
              buffer[0] = 0;
            }
          }
          PyTuple_SetItem(tup, 0, PyString_FromString(buffer));
          PyTuple_SetItem(tup, 1, PyInt_FromLong(a));
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
          WordType buffer;
          tup = PyTuple_New(2);
          {
            const char *color_name = ColorGetName(G, a);
            if(color_name) {
              strcpy(buffer, color_name);
            } else {
              buffer[0] = 0;
            }
          }
          PyTuple_SetItem(tup, 0, PyString_FromString(buffer));
          PyTuple_SetItem(tup, 1, PyInt_FromLong(a));
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
      tup = PyTuple_New(3);
      PyTuple_SetItem(tup, 0, PyFloat_FromDouble(*(rgb++)));
      PyTuple_SetItem(tup, 1, PyFloat_FromDouble(*(rgb++)));
      PyTuple_SetItem(tup, 2, PyFloat_FromDouble(*rgb));
      result = tup;
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
  PyObject *result = NULL;
  char **vla = NULL;
  int ok = false;
  ok = PyArg_ParseTuple(args, "Osi", &self, &str1, &int1);
  if(ok) {
    API_SETUP_PYMOL_GLOBALS;
    ok = (G != NULL);
  } else {
    API_HANDLE_ERROR;
  }
  if(ok && (ok = APIEnterNotModal(G))) {
    vla = (char**) ExecutiveGetChains(G, str1, int1);
    APIExit(G);

    if(vla) {
      result = PConvStringListToPyList(VLAGetSize(vla), vla);
      VLAFreeP(vla);
    }
  }
  if(result) {
    return (APIAutoNone(result));
  } else {
    return APIFailure();
  }
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
  OrthoLineType s1;
  PyObject *range, *color;
  ok = PyArg_ParseTuple(args, "OssOOisfffii", &self, &name, &map, &range, &color,
                        &state, &sele, &beyond, &within, &sigma, &zero, &quiet);
  if(ok) {
    API_SETUP_PYMOL_GLOBALS;
    ok = (G != NULL);
  } else {
    API_HANDLE_ERROR;
  }
  if(ok && (ok = APIEnterNotModal(G))) {
    ok = (SelectorGetTmp(G, sele, s1) >= 0);
    if(ok)
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
    if(ok)
      ok = ExecutiveRampNew(G, name, map, range_vla,
                            color_vla, state, s1, beyond, within, sigma,
                            zero, calc_mode, quiet);
    SelectorFreeTmp(G, s1);
    APIExit(G);
  }
  return APIResultOk(ok);
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
  float grid[3];
  float buffer, floor, ceiling, resolution;
  int type;
  int state;
  int have_corners;
  int quiet, zoom;
  int normalize;
  char *selection;
  OrthoLineType s1 = "";
  int ok = false;
  ok = PyArg_ParseTuple(args, "Osifsf(ffffff)iiiiifff",
                        &self, &name, &type, &grid[0], &selection, &buffer,
                        &minCorner[0], &minCorner[1], &minCorner[2],
                        &maxCorner[0], &maxCorner[1], &maxCorner[2],
                        &state, &have_corners, &quiet, &zoom, &normalize,
                        &floor, &ceiling, &resolution);
  if(ok) {
    API_SETUP_PYMOL_GLOBALS;
    ok = (G != NULL);
  } else {
    API_HANDLE_ERROR;
  }
  if(ok && (ok = APIEnterNotModal(G))) {
    grid[1] = grid[0];
    grid[2] = grid[0];
    ok = (SelectorGetTmp(G, selection, s1) >= 0);
    if(ok)
      ok = ExecutiveMapNew(G, name, type, grid, s1, buffer,
                           minCorner, maxCorner, state, have_corners, quiet, zoom,
                           normalize, floor, ceiling, resolution);
    SelectorFreeTmp(G, s1);
    APIExit(G);

  }
  return APIResultOk(ok);
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
  int ok = false;

  ok =
    PyArg_ParseTuple(args, "Osisiiii", &self, &name, &operator_, &operands, &target_state,
                     &source_state, &zoom, &quiet);
  if(ok) {
    API_SETUP_PYMOL_GLOBALS;
    ok = (G != NULL);
  } else {
    API_HANDLE_ERROR;
  }
  if(ok && (ok = APIEnterNotModal(G))) {
    ok =
      ExecutiveMapSet(G, name, operator_, operands, target_state, source_state, zoom,
                      quiet);
    APIExit(G);
  }
  return APIResultOk(ok);
}

static PyObject *CmdMapTrim(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  char *name, *sele;
  int map_state, sele_state;
  int ok = false;
  float buffer;
  int quiet;
  OrthoLineType s1;
  ok = PyArg_ParseTuple(args, "Ossfiii", &self, &name, &sele, &buffer,
                        &map_state, &sele_state, &quiet);
  if(ok) {
    API_SETUP_PYMOL_GLOBALS;
    ok = (G != NULL);
  } else {
    API_HANDLE_ERROR;
  }
  if(ok && (ok = APIEnterNotModal(G))) {
    ok = (SelectorGetTmp2(G, sele, s1) >= 0);
    ok = ExecutiveMapTrim(G, name, s1, buffer, map_state, sele_state, quiet);
    SelectorFreeTmp(G, s1);
    APIExit(G);
  }
  return APIResultOk(ok);
}

static PyObject *CmdMapDouble(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  char *name;
  int state;
  int ok = false;
  ok = PyArg_ParseTuple(args, "Osi", &self, &name, &state);
  if(ok) {
    API_SETUP_PYMOL_GLOBALS;
    ok = (G != NULL);
  } else {
    API_HANDLE_ERROR;
  }
  if(ok && (ok = APIEnterNotModal(G))) {
    ok = ExecutiveMapDouble(G, name, state);
    APIExit(G);
  }
  return APIResultOk(ok);
}

static PyObject *CmdMapHalve(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  char *name;
  int state;
  int ok = false;
  int smooth;
  ok = PyArg_ParseTuple(args, "Osii", &self, &name, &state, &smooth);
  if(ok) {
    API_SETUP_PYMOL_GLOBALS;
    ok = (G != NULL);
  } else {
    API_HANDLE_ERROR;
  }
  if(ok && (ok = APIEnterNotModal(G))) {
    ok = ExecutiveMapHalve(G, name, state, smooth);
    APIExit(G);
  }
  return APIResultOk(ok);
}

static PyObject *CmdGetRenderer(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  char *vendor = NULL, *renderer = NULL, *version = NULL;
  int ok = false;
  ok = PyArg_ParseTuple(args, "O", &self);
  if(ok) {
    API_SETUP_PYMOL_GLOBALS;
    ok = (G != NULL);
  } else {
    API_HANDLE_ERROR;
  }
  if(ok && (ok = APIEnterNotModal(G))) {
    SceneGetCardInfo(G, &vendor, &renderer, &version);
    APIExit(G);
  }
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

static PyObject *CmdTranslateAtom(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  char *str1;
  int state, log, mode;
  float v[3];
  int ok = false;
  ok =
    PyArg_ParseTuple(args, "Osfffiii", &self, &str1, v, v + 1, v + 2, &state, &mode,
                     &log);
  if(ok) {
    API_SETUP_PYMOL_GLOBALS;
    ok = (G != NULL);
  } else {
    API_HANDLE_ERROR;
  }
  if(ok && (ok = APIEnterNotModal(G))) {
    ok = ExecutiveTranslateAtom(G, str1, v, state, mode, log);
    APIExit(G);
  }
  return APIResultOk(ok);
}

static PyObject *CmdMatrixCopy(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  char *source_name, *target_name;
  int source_mode, target_mode;
  int source_state, target_state, target_undo;
  int log;
  int quiet;
  int ok = false;
  ok = PyArg_ParseTuple(args, "Ossiiiiiii", &self,
                        &source_name, &target_name,
                        &source_mode, &target_mode,
                        &source_state, &target_state, &target_undo, &log, &quiet);
  if(ok) {
    API_SETUP_PYMOL_GLOBALS;
    ok = (G != NULL);
  } else {
    API_HANDLE_ERROR;
  }
  if(ok && (ok = APIEnterNotModal(G))) {
    ExecutiveMatrixCopy(G,
                        source_name, target_name,
                        source_mode, target_mode,
                        source_state, target_state, target_undo, log, quiet);
    APIExit(G);
  }
  return APIResultOk(ok);
}

static PyObject *CmdResetMatrix(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  char *name;
  int mode;
  int state;
  int log;
  int quiet;
  int ok = false;
  ok = PyArg_ParseTuple(args, "Osiiii", &self, &name, &mode, &state, &log, &quiet);
  if(ok) {
    API_SETUP_PYMOL_GLOBALS;
    ok = (G != NULL);
  } else {
    API_HANDLE_ERROR;
  }
  if(ok && (ok = APIEnterNotModal(G))) {
    ExecutiveResetMatrix(G, name, mode, state, log, quiet);
    APIExit(G);
  }
  return APIResultOk(ok);
}

static PyObject *CmdTransformObject(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  char *name, *sele;
  int state, log;
  PyObject *m;
  float matrix[16];
  int homo;
  int ok = false;
  ok = PyArg_ParseTuple(args, "OsiOisi", &self, &name, &state, &m, &log, &sele, &homo);
  if(ok) {
    API_SETUP_PYMOL_GLOBALS;
    ok = (G != NULL);
  } else {
    API_HANDLE_ERROR;
  }
  if(ok) {
    if(PConvPyListToFloatArrayInPlace(m, matrix, 16) > 0) {
      if((ok = APIEnterNotModal(G))) {
        int matrix_mode = SettingGetGlobal_i(G, cSetting_matrix_mode);
        if(matrix_mode<0) matrix_mode = 0;
        if((matrix_mode == 0) || (sele[0] != 0)) {
          ok = ExecutiveTransformObjectSelection(G, name,
                                                 state, sele, log, matrix, homo, true);
        } else {
          /* state? */
          ok = ExecutiveCombineObjectTTT(G, name, matrix, false, 
                                         SettingGetGlobal_i(G, cSetting_movie_auto_store));
        }
        APIExit(G);
      }
    } else {
      PRINTFB(G, FB_CCmd, FB_Errors)
        "CmdTransformObject-DEBUG: bad matrix\n" ENDFB(G);
      ok = false;
    }
  }
  return APIResultOk(ok);
}

static PyObject *CmdTransformSelection(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  char *sele;
  int state, log;
  int homo;
  PyObject *m;
  float ttt[16];
  int ok = false;
  ok = PyArg_ParseTuple(args, "OsiOii", &self, &sele, &state, &m, &log, &homo);
  if(ok) {
    API_SETUP_PYMOL_GLOBALS;
    ok = (G != NULL);
  } else {
    API_HANDLE_ERROR;
  }
  if(ok && (ok = APIEnterNotModal(G))) {
    if(PConvPyListToFloatArrayInPlace(m, ttt, 16) > 0) {
      ok = ExecutiveTransformSelection(G, state, sele, log, ttt, homo);
    } else {
      PRINTFB(G, FB_CCmd, FB_Errors)
        "CmdTransformSelection-DEBUG: bad matrix\n" ENDFB(G);
      ok = false;
    }
    APIExit(G);
  }
  return APIResultOk(ok);
}

static PyObject *CmdLoadColorTable(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  char *str1;
  float gamma;
  int ok = false;
  int quiet;
  ok = PyArg_ParseTuple(args, "Osfi", &self, &str1, &gamma, &quiet);
  if(ok) {
    API_SETUP_PYMOL_GLOBALS;
    ok = (G != NULL);
  } else {
    API_HANDLE_ERROR;
  }
  if(ok && (ok = APIEnterNotModal(G))) {
    ok = ColorTableLoad(G, str1, gamma, quiet);
    APIExit(G);
  }
  return APIResultOk(ok);
}

static PyObject *CmdLoadPNG(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  char *str1;
  int ok = false;
  int quiet;
  int movie, stereo;
  ok = PyArg_ParseTuple(args, "Osiii", &self, &str1, &movie, &stereo, &quiet);
  if(ok) {
    API_SETUP_PYMOL_GLOBALS;
    ok = (G != NULL);
  } else {
    API_HANDLE_ERROR;
  }
  if(ok && (ok = APIEnterNotModal(G))) {
    ok = SceneLoadPNG(G, str1, movie, stereo, quiet);
    APIExit(G);
  }
  return APIResultOk(ok);
}

static PyObject *CmdBackgroundColor(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  char *str1;
  int ok = false;
  ok = PyArg_ParseTuple(args, "Os", &self, &str1);
  if(ok) {
    API_SETUP_PYMOL_GLOBALS;
    ok = (G != NULL);
  } else {
    API_HANDLE_ERROR;
  }
  if(ok && (ok = APIEnterNotModal(G))) {
    ok = SettingSet_color(G->Setting, cSetting_bg_rgb, str1);
    SettingGenerateSideEffects(G, cSetting_bg_rgb, NULL, -1, 0);
    APIExit(G);
  }
  return APIResultOk(ok);
}

static PyObject *CmdGetPosition(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  PyObject *result;
  float v[3] = { 0.0F, 0.0F, 0.0F };
  int ok = false;
  ok = PyArg_ParseTuple(args, "O", &self);
  if(ok) {
    API_SETUP_PYMOL_GLOBALS;
    ok = (G != NULL);
  } else {
    API_HANDLE_ERROR;
  }
  if(ok && (ok = APIEnterNotModal(G))) {
    SceneGetCenter(G, v);
    APIExit(G);
  }
  result = PConvFloatArrayToPyList(v, 3);
  return (APIAutoNone(result));
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
  PyObject *key = Py_None;
  PyObject *value = Py_None;
  int *iVLA = NULL;
  float *pVLA = NULL, *sVLA = NULL;
  int l = 0;
  int *i;
  ObjectMolecule **o, **oVLA = NULL;
  int a;
  float *s, *p;
  int ok = PyArg_ParseTuple(args, "Osi", &self, &str1, &state);
  if(ok) {
    API_SETUP_PYMOL_GLOBALS;
    ok = (G != NULL);
  } else {
    API_HANDLE_ERROR;
  }
  if(ok && (ok = APIEnterNotModal(G))) {
    l = ExecutivePhiPsi(G, str1, &oVLA, &iVLA, &pVLA, &sVLA, state);
    APIExit(G);
    if(iVLA) {
      result = PyDict_New();
      i = iVLA;
      o = oVLA;
      p = pVLA;
      s = sVLA;
      for(a = 0; a < l; a++) {
        key = PyTuple_New(2);
        PyTuple_SetItem(key, 1, PyInt_FromLong(*(i++) + 1));    /* +1 for index */
        PyTuple_SetItem(key, 0, PyString_FromString((*(o++))->Obj.Name));
        value = PyTuple_New(2);
        PyTuple_SetItem(value, 0, PyFloat_FromDouble(*(p++)));  /* +1 for index */
        PyTuple_SetItem(value, 1, PyFloat_FromDouble(*(s++)));
        PyDict_SetItem(result, key, value);
        Py_DECREF(key);
        Py_DECREF(value);
      }
    } else {
      result = PyDict_New();
    }
    VLAFreeP(iVLA);
    VLAFreeP(oVLA);
    VLAFreeP(sVLA);
    VLAFreeP(pVLA);
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

  if(!PyArg_ParseTuple(args, "Os|i", &self, &str1, &state)) {
    API_HANDLE_ERROR;
    ok_raise(2);
  }

  ok_assert(2, str1[0]);
  API_SETUP_PYMOL_GLOBALS;
  ok_assert(2, G && APIEnterBlockedNotModal(G));

  if(SelectorGetTmp(G, str1, s1) >= 0) {
    int sele1 = SelectorIndexByName(G, s1);
    if(sele1 >= 0) {
      int unblock = PAutoBlock(G);
      result = SelectorGetCoordsAsNumPy(G, sele1, state);
      PAutoUnblock(G, unblock);
    }
    SelectorFreeTmp(G, s1);
  }

  APIExitBlocked(G);
ok_except2:
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

  if(!PyArg_ParseTuple(args, "Os|ih", &self, &name, &state, &copy)) {
    API_HANDLE_ERROR;
    ok_raise(2);
  }

  ok_assert(2, name[0]);
  ok_assert(2, state >= 0);

  API_SETUP_PYMOL_GLOBALS;
  ok_assert(2, G && APIEnterBlockedNotModal(G));

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
  int ok = false;
  ok = PyArg_ParseTuple(args, "O", &self);
  if(ok) {
    API_SETUP_PYMOL_GLOBALS;
    ok = (G != NULL);
  } else {
    API_HANDLE_ERROR;
  }
  if(ok && (ok = APIEnterNotModal(G))) {
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
  } else {
    return (APIAutoNone(NULL));
  }
}

static PyObject *CmdGetViewPort(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  int width, height;
  int ok = false;
  ok = PyArg_ParseTuple(args, "O", &self);
  if(ok) {
    API_SETUP_PYMOL_GLOBALS;
    ok = (G != NULL);
  } else {
    API_HANDLE_ERROR;
  }
  if(ok && (ok = APIEnterNotModal(G))) {
    SceneGetWidthHeight(G, &width, &height);
    APIExit(G);
    return (Py_BuildValue
            ("(ii)", width, height
            ));
  } else {
    return (APIAutoNone(NULL));
  }
}

static PyObject *CmdSetView(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  SceneViewType view;
  int quiet;
  float animate;
  int hand;
  int ok = PyArg_ParseTuple(args, "O(fffffffffffffffffffffffff)ifi",
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
  if(ok) {
    API_SETUP_PYMOL_GLOBALS;
    ok = (G != NULL);
  } else {
    API_HANDLE_ERROR;
  }
  if(ok && (ok = APIEnterNotModal(G))) {
    SceneSetView(G, view, quiet, animate, hand);        /* TODO STATUS */
    APIExit(G);
  }
  return APIResultOk(ok);
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
  int ok = false;
  ok = PyArg_ParseTuple(args, "Osis", &self, &str1, &int1, &str2);
  if(ok) {
    API_SETUP_PYMOL_GLOBALS;
    ok = (G != NULL);
  } else {
    API_HANDLE_ERROR;
  }
  if(ok && (ok = APIEnterNotModal(G))) {
    ok = ExecutiveSetTitle(G, str1, int1, str2);
    APIExit(G);
  }
  return APIResultOk(ok);
}

static PyObject *CmdGetTitle(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  char *str1;
  int int1;
  int ok = false;
  PyObject *result = Py_None;
  ok = PyArg_ParseTuple(args, "Osi", &self, &str1, &int1);
  if(ok) {
    API_SETUP_PYMOL_GLOBALS;
    ok = (G != NULL);
  } else {
    API_HANDLE_ERROR;
  }
  if(ok && (ok = APIEnterNotModal(G))) {
    const char *str2 = ExecutiveGetTitle(G, str1, int1);
    APIExit(G);
    if(str2)
      result = PyString_FromString(str2);
  }
  return (APIAutoNone(result));
}

static PyObject *CmdExportCoords(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  void *result;
  char *str1;
  int int1;
  PyObject *py_result = Py_None;
  int ok = false;
  ok = PyArg_ParseTuple(args, "Osi", &self, &str1, &int1);
  if(ok) {
    API_SETUP_PYMOL_GLOBALS;
    ok = (G != NULL);
  } else {
    API_HANDLE_ERROR;
  }
  if(ok && (ok = APIEnterNotModal(G))) {
    result = ExportCoordsExport(G, str1, int1, 0);
    APIExit(G);
    if(result)
      py_result = PyCObject_FromVoidPtr(result, (void (*)(void *)) ExportCoordsFree);
  }
  return (APIAutoNone(py_result));
}

static PyObject *CmdImportCoords(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  char *str1;
  int int1;
  PyObject *cObj;
  void *mmdat = NULL;
  int ok = false;
  ok = PyArg_ParseTuple(args, "OsiO", &self, &str1, &int1, &cObj);
  if(ok) {
    API_SETUP_PYMOL_GLOBALS;
    ok = (G != NULL);
  } else {
    API_HANDLE_ERROR;
  }
  if(ok) {
    if(PyCObject_Check(cObj))
      mmdat = PyCObject_AsVoidPtr(cObj);
    if((ok = APIEnterNotModal(G))) {
      if(mmdat)
        ok = ExportCoordsImport(G, str1, int1, (ExportCoords*) mmdat, 0);
      APIExit(G);
    }
  }
  return APIResultOk(ok);
}

static PyObject *CmdGetArea(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  char *str1;
  int int1, int2;
  float result = -1.0;
  int ok = false;
  ok = PyArg_ParseTuple(args, "Osii", &self, &str1, &int1, &int2);
  if(ok) {
    API_SETUP_PYMOL_GLOBALS;
    ok = (G != NULL);
  } else {
    API_HANDLE_ERROR;
  }
  if(ok && (ok = APIEnterNotModal(G))) {
    result = ExecutiveGetArea(G, str1, int1, int2);
    APIExit(G);

  }
  return (Py_BuildValue("f", result));

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
  WordType type = "";
  int ok = false;
  ok = PyArg_ParseTuple(args, "Os", &self, &str1);
  if(ok) {
    API_SETUP_PYMOL_GLOBALS;
    ok = (G != NULL);
  } else {
    API_HANDLE_ERROR;
  }
  if(ok && (ok = APIEnterNotModal(G))) {
    ok = ExecutiveGetType(G, str1, type);
    APIExit(G);
  }
  if(ok)
    return (Py_BuildValue("s", type));
  else
    return APIResultOk(ok);
}

static PyObject *CmdGetObjectSettings(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  CObject *obj = NULL;
  PyObject *result = NULL;
  const char *oname;
  int state = -1;

  if (!PyArg_ParseTuple(args, "Os|i", &self, &oname, &state)) {
    API_HANDLE_ERROR;
    ok_raise(1);
  }

  API_SETUP_PYMOL_GLOBALS;
  ok_assert(1, G && APIEnterNotModal(G));

  obj = ExecutiveFindObjectByName(G, oname);

  if(!obj) {
    ErrMessage(G, "GetObjectSettings", "named object not found.");
  } else if (obj->fGetSettingHandle) {
    auto handle = obj->fGetSettingHandle(obj, -1);

    if (state != -1) {
      auto handle_state = obj->fGetSettingHandle(obj, state);

      // only accept handle if different from object-level settings
      handle = (handle_state == handle) ? NULL : handle_state;
    }

    if (handle) {
      result = SettingAsPyList(*handle, true);
    }
  }

  APIExit(G);
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
      CObject *obj = EditorDragObject(G);    
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
  PyObject *result = Py_None;
  int ok = false;
  char *str0;
  ok = PyArg_ParseTuple(args, "Os", &self, &str0);
  if(ok) {
    API_SETUP_PYMOL_GLOBALS;
    ok = (G != NULL);
  } else {
    API_HANDLE_ERROR;
  }
  if(ok) {
    APIEnter(G);
    UtilNCopy(name, str0, sizeof(WordType));
    ObjectMakeValidName(G, name);
    APIExit(G);
    result = PyString_FromString(name);
  }
  return (APIAutoNone(result));
}

static PyObject *CmdGetNames(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  int int1, int2;
  char *vla = NULL;
  OrthoLineType s0 = "";
  PyObject *result = Py_None;
  int ok = false;
  char *str0;
  ok = PyArg_ParseTuple(args, "Oiis", &self, &int1, &int2, &str0);
  if(ok) {
    API_SETUP_PYMOL_GLOBALS;
    ok = (G != NULL);
  } else {
    API_HANDLE_ERROR;
  }
  if(ok && (ok = APIEnterNotModal(G))) {
    if(str0[0])
      ok = (SelectorGetTmp(G, str0, s0) >= 0);
    vla = ExecutiveGetNames(G, int1, int2, s0);
    if(s0[0])
      SelectorFreeTmp(G, s0);
    APIExit(G);
    result = PConvStringVLAToPyList(vla);
    VLAFreeP(vla);
  }
  return (APIAutoNone(result));
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
  int ok = false;
  ok = PyArg_ParseTuple(args, "Oi", &self, &int1);
  if(ok) {
    API_SETUP_PYMOL_GLOBALS;
    ok = (G != NULL);
  } else {
    API_HANDLE_ERROR;
  }
  if(ok && (ok = APIEnterNotModal(G))) {
    ok = ExecutiveInvert(G, int1);
    APIExit(G);
  }
  return APIResultOk(ok);
}

static PyObject *CmdTorsion(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  float float1;
  int ok = false;
  ok = PyArg_ParseTuple(args, "Of", &self, &float1);
  if(ok) {
    API_SETUP_PYMOL_GLOBALS;
    ok = (G != NULL);
  } else {
    API_HANDLE_ERROR;
  }
  if(ok && (ok = APIEnterNotModal(G))) {
    ok = EditorTorsion(G, float1);
    APIExit(G);
  }
  return APIResultOk(ok);
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

  int ok = false;
  ok = PyArg_ParseTuple(args, "Osii", &self, &str1, &int1, &quiet);
  if(ok) {
    API_SETUP_PYMOL_GLOBALS;
    ok = (G != NULL);
  } else {
    API_HANDLE_ERROR;
  }
  if(ok && (ok = APIEnterNotModal(G))) {
    ExecutiveMask(G, str1, int1, quiet);  /* TODO STATUS */
    APIExit(G);
  }
  return APIResultOk(ok);
}

static PyObject *CmdProtect(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  char *str1;
  int int1, int2;

  int ok = false;
  ok = PyArg_ParseTuple(args, "Osii", &self, &str1, &int1, &int2);
  if(ok) {
    API_SETUP_PYMOL_GLOBALS;
    ok = (G != NULL);
  } else {
    API_HANDLE_ERROR;
  }
  if(ok && (ok = APIEnterNotModal(G))) {
    ExecutiveProtect(G, str1, int1, int2);        /* TODO STATUS */
    APIExit(G);

  }
  return APIResultOk(ok);
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
  int ok = false;
  ok = PyArg_ParseTuple(args, "Oiii", &self, &i1, &i2, &i3);
  if(ok) {
    API_SETUP_PYMOL_GLOBALS;
    ok = (G != NULL);
  } else {
    API_HANDLE_ERROR;
  }
  if(ok && (ok = APIEnterNotModal(G))) {
    switch (i1) {               /* TODO STATUS */
    case 0:
      FeedbackSetMask(G, i2, (uchar) i3);
      break;
    case 1:
      FeedbackEnable(G, i2, (uchar) i3);
      break;
    case 2:
      FeedbackDisable(G, i2, (uchar) i3);
      break;
    case 3:
      FeedbackPush(G);
      break;
    case 4:
      FeedbackPop(G);
      break;
    }
    APIExit(G);
  }
  return APIResultOk(ok);
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
  int ok = false;
  int ver;
  ok = PyArg_ParseTuple(args, "Oi", &self, &ver);
  if(ok) {
    API_SETUP_PYMOL_GLOBALS;
    ok = (G != NULL);
  } else {
    API_HANDLE_ERROR;
  }
  if(ok) {
    char *vla = NULL;
    if((ok = APIEnterNotModal(G))) {
      SceneRay(G, 0, 0, (ver == 1) ? 6 : 4,     /* VRML1 or 2? */
               NULL, &vla, 0.0F, 0.0F, false, NULL, false, -1);
      APIExit(G);
    }
    if(vla) {
      result = Py_BuildValue("s", vla);
    }
    VLAFreeP(vla);
  }
  return (APIAutoNone(result));
}

/*
 * Return a COLLADA string or None on failure
 */
static PyObject *CmdGetCOLLADA(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  PyObject *result = NULL;
  int ver;
  char *vla = NULL;

  ok_assert(1, PyArg_ParseTuple(args, "Oi", &self, &ver));
  API_SETUP_PYMOL_GLOBALS;
  ok_assert(1, G && APIEnterNotModal(G));

  SceneRay(G, 0, 0, 8,     /* mode 8 = COLLADA */
      NULL, &vla, 0.0F, 0.0F, false, NULL, false, -1);

  ok_assert(2, vla && vla[0]);
  result = Py_BuildValue("s", vla);

ok_except2:
  VLAFreeP(vla);

  APIExit(G);
  return APIAutoNone(result);
ok_except1:
  API_HANDLE_ERROR;
  return APIAutoNone(NULL);
}


static PyObject *CmdGetIdtf(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  PyObject *result = NULL;
  int ok = false;
  ok = PyArg_ParseTuple(args, "O", &self);
  if(ok) {
    API_SETUP_PYMOL_GLOBALS;
    ok = (G != NULL);
  } else {
    API_HANDLE_ERROR;
  }
  if(ok) {
    char *node = NULL, *rsrc = NULL;
    if((ok = APIEnterNotModal(G))) {
      SceneRay(G, 0, 0, cSceneRay_MODE_IDTF,
               &node, &rsrc, 0.0F, 0.0F, false, NULL, false, -1);
      APIExit(G);
    }
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
  int ok = false;
  ok = PyArg_ParseTuple(args, "O", &self);
  if(ok) {
    API_SETUP_PYMOL_GLOBALS;
    ok = (G != NULL);
  } else {
    API_HANDLE_ERROR;
  }
  if(ok) {
    char *header = NULL, *geom = NULL;
    if((ok = APIEnterNotModal(G))) {
      SceneRay(G, 0, 0, 1, &header, &geom, 0.0F, 0.0F, false, NULL, false, -1);
      APIExit(G);
    }
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
  int ok = false;
  ok = PyArg_ParseTuple(args, "O", &self);
  if(ok) {
    API_SETUP_PYMOL_GLOBALS;
    ok = (G != NULL);
  } else {
    API_HANDLE_ERROR;
  }
  if(ok) {
    char *obj = NULL, *mtl = NULL;
    if((ok = APIEnterNotModal(G))) {
      SceneRay(G, 0, 0, 5, &obj, &mtl, 0.0F, 0.0F, false, NULL, false, -1);
      APIExit(G);
    }
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
  PyMOLGlobals *G = NULL;
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
    result = WizardGet(G);
    APIExit(G);
  }
  if(!result)
    result = Py_None;
  return APIIncRef(result);
}

static PyObject *CmdGetWizardStack(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  PyObject *result = NULL;
  int ok = false;
  ok = PyArg_ParseTuple(args, "O", &self);
  if(ok) {
    API_SETUP_PYMOL_GLOBALS;
    ok = (G != NULL);
  } else {
    API_HANDLE_ERROR;
  }
  if(ok && (ok = APIEnterBlockedNotModal(G))) {
    result = WizardGetStack(G);
    APIExitBlocked(G);
  }
  if(!result)
    result = Py_None;
  return APIIncRef(result);
}

static PyObject *CmdSetWizard(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;

  PyObject *obj;
  int ok = false;
  int replace;
  ok = PyArg_ParseTuple(args, "OOi", &self, &obj, &replace);
  if(ok) {
    API_SETUP_PYMOL_GLOBALS;
    ok = (G != NULL);
  } else {
    API_HANDLE_ERROR;
  }
  if(ok) {
    if(!obj)
      ok = false;
    else if((ok = APIEnterNotModal(G))) {
      WizardSet(G, obj, replace);       /* TODO STATUS */
      APIExit(G);
    }
  }
  return APIResultOk(ok);
}

static PyObject *CmdSetWizardStack(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;

  PyObject *obj;
  int ok = false;
  ok = PyArg_ParseTuple(args, "OO", &self, &obj);
  if(ok) {
    API_SETUP_PYMOL_GLOBALS;
    ok = (G != NULL);
  } else {
    API_HANDLE_ERROR;
  }
  if(ok) {
    if(!obj)
      ok = false;
    else if((ok = APIEnterNotModal(G))) {
      WizardSetStack(G, obj);   /* TODO STATUS */
      APIExit(G);
    }
  }
  return APIResultOk(ok);
}

static PyObject *CmdRefreshWizard(PyObject * self, PyObject * args)
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
    WizardRefresh(G);
    OrthoInvalidateDoDraw(G);
    OrthoDirty(G);
    APIExit(G);
  }
  return APISuccess();
}

static PyObject *CmdDirtyWizard(PyObject * self, PyObject * args)
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
    WizardDirty(G);
    APIExit(G);
  }
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
  int ok = false;
  ok = PyArg_ParseTuple(args, "Oss", &self, &str1, &str2);
  if(ok) {
    API_SETUP_PYMOL_GLOBALS;
    ok = (G != NULL);
  } else {
    API_HANDLE_ERROR;
  }
  if(ok && (ok = APIEnterNotModal(G))) {
    ExecutiveDump(G, str1, str2);       /* TODO STATUS */
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
  int box_mode;
  float carve;
  int ok = false;
  int map_state;
  int quiet;
  /* oper 0 = all, 1 = sele + buffer, 2 = vector */

  ok = PyArg_ParseTuple(args, "Ossisffiifiif", &self, &mesh_name, &map_name, &box_mode,
                        &sele, &fbuf, &lvl, &mesh_mode, &state, &carve, &map_state,
                        &quiet, &alt_lvl);
  if(ok) {
    API_SETUP_PYMOL_GLOBALS;
    ok = (G != NULL);
  } else {
    API_HANDLE_ERROR;
  }
  if(ok && (ok = APIEnterNotModal(G))) {
    ok = ExecutiveIsomeshEtc(G, mesh_name, map_name, lvl, sele, fbuf,
                             state, carve, map_state, quiet, mesh_mode, box_mode,
                             alt_lvl);

    APIExit(G);
  }
  return APIResultOk(ok);
}

static int ExecutiveSliceNew(PyMOLGlobals * G, char *slice_name,
                             char *map_name, int state, int map_state)
{
  int ok = true;
  int multi = false;
  CObject *obj = NULL, *mObj, *origObj;
  ObjectMap *mapObj;
  ObjectMapState *ms;

  origObj = ExecutiveFindObjectByName(G, slice_name);
  if(origObj) {
    if(origObj->type != cObjectSlice) {
      ExecutiveDelete(G, slice_name);
      origObj = NULL;
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
        if(origObj->fGetNFrame)
          state = origObj->fGetNFrame(origObj);
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
        map_state = ObjectMapGetNStates(mapObj) - 1;
      ms = ObjectMapStateGetActive(mapObj, map_state);
      if(ms) {
        obj = (CObject *) ObjectSliceFromMap(G, (ObjectSlice *) origObj, mapObj,
                                             state, map_state);

        if(!origObj) {
          ObjectSetName(obj, slice_name);
          ExecutiveManageObject(G, (CObject *) obj, -1, false);
        }
        PRINTFB(G, FB_ObjectMesh, FB_Actions)
          " SliceMap: created \"%s\".\n", slice_name ENDFB(G);

      } else if(!multi) {
        PRINTFB(G, FB_ObjectSlice, FB_Warnings)
          " SliceMap-Warning: state %d not present in map \"%s\".\n", map_state + 1,
          map_name ENDFB(G);
        ok = false;
      }
      if(multi) {
        origObj = obj;
        map_state++;
        state++;
        if(map_state >= mapObj->NState)
          break;
      } else {
        break;
      }
    }
  } else {
    PRINTFB(G, FB_ObjectSlice, FB_Errors)
      " SliceMap: Map or brick object \"%s\" not found.\n", map_name ENDFB(G);
    ok = false;
  }
  return ok;
}

static PyObject *CmdSliceNew(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  int ok = false;
  char *slice;
  char *map;
  int state, map_state;

  ok = PyArg_ParseTuple(args, "Ossii", &self, &slice, &map, &state, &map_state);
  if(ok) {
    API_SETUP_PYMOL_GLOBALS;
    ok = (G != NULL);
  } else {
    API_HANDLE_ERROR;
  }
  if(ok && (ok = APIEnterNotModal(G))) {
    ok = ExecutiveSliceNew(G, slice, map, state, map_state);
    APIExit(G);
  }
  return APIResultOk(ok);
}

static PyObject *CmdIsosurface(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  char *surf_name, *map_name, *sele;
  float lvl, fbuf;
  int surf_mode;
  int state = -1;
  int box_mode;
  float carve;
  int ok = false;
  int map_state = 0;
  int side;
  int quiet;
  /* box_mode 0 = all, 1 = sele + buffer, 2 = vector */

  ok = PyArg_ParseTuple(args, "Ossisffiifiii", &self, &surf_name, &map_name, &box_mode,
                        &sele, &fbuf, &lvl, &surf_mode, &state, &carve, &map_state,
                        &side, &quiet);
  if(ok) {
    API_SETUP_PYMOL_GLOBALS;
    ok = (G != NULL);
  } else {
    API_HANDLE_ERROR;
  }
  if(ok && (ok = APIEnterNotModal(G))) {

    ok = ExecutiveIsosurfaceEtc(G, surf_name, map_name, lvl, sele, fbuf, state,
                                carve, map_state, side, quiet, surf_mode, box_mode);

    APIExit(G);
  }
  return APIResultOk(ok);
}

static PyObject *CmdSymExp(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  char *str1, *str2, *str3;
  float cutoff;
  CObject *mObj;
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
  int source_mode, target_mode;
  int source_state, target_state, target_undo;
  int log;
  int quiet;
  int ok = false;

  ok = PyArg_ParseTuple(args, "Ossiiiiiii", &self,
                        &source_name, &target_name,
                        &source_mode, &target_mode,
                        &source_state, &target_state, &target_undo, &log, &quiet);
  if(ok) {
    API_SETUP_PYMOL_GLOBALS;
    ok = (G != NULL);
  } else {
    API_HANDLE_ERROR;
  }
  if(ok && (ok = APIEnterNotModal(G))) {
    ExecutiveSymmetryCopy(G,
			  source_name, target_name,
			  source_mode, target_mode,
			  source_state, target_state, target_undo, log, quiet);
    APIExit(G);
  }

  return APIResultOk(ok);
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
  return (Py_BuildValue("f", overlap));
}

static PyObject *CmdDist(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  char *name, *str1, *str2;
  float cutoff, result = -1.0;
  int labels, quiet;
  int mode, reset, state, zoom;
  int state1, state2;
  int ok = false;
  ok = PyArg_ParseTuple(args, "Osssifiiiiiii", &self, &name, &str1,
                        &str2, &mode, &cutoff, &labels, &quiet, &reset, &state, &zoom,
                        &state1, &state2);
  if(ok) {
    API_SETUP_PYMOL_GLOBALS;
    ok = (G != NULL);
  } else {
    API_HANDLE_ERROR;
  }
  if(ok && (ok = APIEnterNotModal(G))) {
    ok = ExecutiveDist(G, &result, name, str1, str2, mode, cutoff,
        labels, quiet, reset, state, zoom, state1, state2);
    APIExit(G);
  }
  if(!ok)
    return APIFailure();
  else
    return (Py_BuildValue("f", result));
}

static PyObject *CmdAngle(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  char *name, *str1, *str2, *str3;
  float result = -999.0;
  int labels, quiet;
  int mode;
  int ok = false;
  int reset, zoom;
  int state;
  int state1, state2, state3;
  ok = PyArg_ParseTuple(args, "Ossssiiiiiiiii", &self,
                        &name, &str1, &str2, &str3,
                        &mode, &labels, &reset, &zoom, &quiet, &state,
                        &state1, &state2, &state3);
  if(ok) {
    API_SETUP_PYMOL_GLOBALS;
    ok = (G != NULL);
  } else {
    API_HANDLE_ERROR;
  }
  if(ok && (ok = APIEnterNotModal(G))) {
    ok = ExecutiveAngle(G, &result, name, str1, str2, str3,
        mode, labels, reset, zoom, quiet, state,
        state1, state2, state3);
    APIExit(G);
  }
  return (Py_BuildValue("f", result));
}

static PyObject *CmdDihedral(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  char *name, *str1, *str2, *str3, *str4;
  float result = -999.0;
  int labels, quiet;
  int mode;
  int ok = false;
  int reset, zoom;
  int state;
  ok = PyArg_ParseTuple(args, "Osssssiiiiii", &self,
                        &name, &str1, &str2, &str3, &str4,
                        &mode, &labels, &reset, &zoom, &quiet, &state);
  if(ok) {
    API_SETUP_PYMOL_GLOBALS;
    ok = (G != NULL);
  } else {
    API_HANDLE_ERROR;
  }
  if(ok && (ok = APIEnterNotModal(G))) {
    ok = ExecutiveDihedral(G, &result, name, str1, str2, str3, str4,
        mode, labels, reset, zoom, quiet, state);
    APIExit(G);
  }
  return (Py_BuildValue("f", result));
}

static PyObject *CmdBond(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  char *str1, *str2;
  int order, mode;
  int quiet;
  int ok = false;
  ok = PyArg_ParseTuple(args, "Ossiii", &self, &str1, &str2, &order, &mode, &quiet);
  if(ok) {
    API_SETUP_PYMOL_GLOBALS;
    ok = (G != NULL);
  } else {
    API_HANDLE_ERROR;
  }
  if(ok && (ok = APIEnterNotModal(G))) {
    ok = ExecutiveBond(G, str1, str2, order, mode, quiet);
    APIExit(G);
  }
  return APIResultOk(ok);
}

static PyObject *CmdRevalence(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  char *sele1, *sele2, *source;
  int source_state, target_state, reset;
  int quiet;
  int ok = false;
  ok = PyArg_ParseTuple(args, "Osssiiii", &self, &sele1, &sele2, &source,
                        &target_state, &source_state, &reset, &quiet);
  if(ok) {
    API_SETUP_PYMOL_GLOBALS;
    ok = (G != NULL);
  } else {
    API_HANDLE_ERROR;
  }
  if(ok && (ok = APIEnterNotModal(G))) {
    ok = ExecutiveRevalence(G, sele1, sele2, source, target_state, source_state, reset, quiet);
    APIExit(G);
  }
  return APIResultOk(ok);
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
  OrthoLineType s1;
  int quiet;
  int ok = false;
  ok = PyArg_ParseTuple(args, "Ossi", &self, &str1, &str2, &quiet);
  if(ok) {
    API_SETUP_PYMOL_GLOBALS;
    ok = (G != NULL);
  } else {
    API_HANDLE_ERROR;
  }
  if(ok && (ok = APIEnterNotModal(G))) {
    ok = (SelectorGetTmp(G, str1, s1) >= 0);
    if(ok)
      ok = ExecutiveLabel(G, s1, str2, quiet, cExecutiveLabelEvalOn);
    SelectorFreeTmp(G, s1);
    APIExit(G);
  }
  return APIResultOk(ok);
}

static PyObject *CmdLabel2(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  char *str1, *str2;
  OrthoLineType s1;
  int quiet;
  int ok = false;
  ok = PyArg_ParseTuple(args, "Ossi", &self, &str1, &str2, &quiet);
  if(ok) {
    API_SETUP_PYMOL_GLOBALS;
    ok = (G != NULL);
  } else {
    API_HANDLE_ERROR;
  }
  if(ok && (ok = APIEnterNotModal(G))) {
    ok = (SelectorGetTmp(G, str1, s1) >= 0);
    if(ok)
      ok = ExecutiveLabel(G, s1, str2, quiet, cExecutiveLabelEvalAlt);
    APIExit(G);
  }
  return APIResultOk(ok);
}

static PyObject *CmdAlter(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  char *str1, *str2;
  int i1, quiet;
  int result = 0;
  int ok = false;
  PyObject *space;
  ok = PyArg_ParseTuple(args, "OssiiO", &self, &str1, &str2, &i1, &quiet, &space);
  if(ok) {
    API_SETUP_PYMOL_GLOBALS;
    ok = (G != NULL);
  } else {
    API_HANDLE_ERROR;
  }
  if(ok && (ok = APIEnterNotModal(G))) {
    result = ExecutiveIterate(G, str1, str2, i1, quiet, space);   /* TODO STATUS */
    APIExit(G);
  }
  return Py_BuildValue("i", result);
}

static PyObject *CmdAlterList(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  char *str1;
  OrthoLineType s1;
  int quiet;
  int result = 0;
  int ok = false;
  PyObject *space;
  PyObject *list;
  ok = PyArg_ParseTuple(args, "OsOiO", &self, &str1, &list, &quiet, &space);
  if(ok) {
    API_SETUP_PYMOL_GLOBALS;
    ok = (G != NULL);
  } else {
    API_HANDLE_ERROR;
  }
  if(ok && (ok = APIEnterBlockedNotModal(G))) {
    ok = (SelectorGetTmp(G, str1, s1) >= 0);
    result = ExecutiveIterateList(G, s1, list, false, quiet, space);    /* TODO STATUS */
    SelectorFreeTmp(G, s1);
    APIExitBlocked(G);
  }
  return Py_BuildValue("i", result);
}

static PyObject *CmdSelectList(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  char *str1, *sele_name;
  OrthoLineType s1;
  int quiet;
  int result = 0;
  int ok = false;
  int mode;
  int state;
  PyObject *list;
  ok =
    PyArg_ParseTuple(args, "OssOiii", &self, &sele_name, &str1, &list, &state, &mode,
                     &quiet);
  if(ok) {
    API_SETUP_PYMOL_GLOBALS;
    ok = (G != NULL);
  } else {
    API_HANDLE_ERROR;
  }
  if(ok && (ok = APIEnterBlockedNotModal(G))) {
    int *int_array = NULL;
    ok = (SelectorGetTmp(G, str1, s1) >= 0);
    if(ok)
      ok = PyList_Check(list);
    if(ok)
      ok = PConvPyListToIntArray(list, &int_array);
    if(ok) {
      int list_len = PyList_Size(list);
      result =
        ExecutiveSelectList(G, sele_name, s1, int_array, list_len, state, mode, quiet);
      SceneInvalidate(G);
      SeqDirty(G);
    }
    FreeP(int_array);
    APIExitBlocked(G);
  }
  return Py_BuildValue("i", result);
}

static PyObject *CmdAlterState(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  char *str1, *str2;
  int i1, i2, i3, quiet;
  int result = -1;
  PyObject *obj;
  int ok = false;
  ok =
    PyArg_ParseTuple(args, "OissiiiO", &self, &i1, &str1, &str2, &i2, &i3, &quiet, &obj);
  if(ok) {
    API_SETUP_PYMOL_GLOBALS;
    ok = (G != NULL);
  } else {
    API_HANDLE_ERROR;
  }
  if(ok && (ok = APIEnterNotModal(G))) {
    result = ExecutiveIterateState(G, i1, str1, str2, i2, i3, quiet, obj);
    APIExit(G);
  }
  return PyInt_FromLong(result);
}

static PyObject *CmdCopy(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  char *str1, *str2;
  int ok = false;
  int zoom;
  ok = PyArg_ParseTuple(args, "Ossi", &self, &str1, &str2, &zoom);
  if(ok) {
    API_SETUP_PYMOL_GLOBALS;
    ok = (G != NULL);
  } else {
    API_HANDLE_ERROR;
  }
  if(ok && (ok = APIEnterNotModal(G))) {
    ExecutiveCopy(G, str1, str2, zoom); /* TODO STATUS */
    APIExit(G);
  }
  return APIResultOk(ok);
}

static PyObject *CmdRecolor(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  char *str1;
  OrthoLineType s1;
  int ok = false;
  int rep = -1;
  ok = PyArg_ParseTuple(args, "Osi", &self, &str1, &rep);
  if(ok) {
    API_SETUP_PYMOL_GLOBALS;
    ok = (G != NULL);
  } else {
    API_HANDLE_ERROR;
  }
  if(ok) {
    PRINTFD(G, FB_CCmd)
      " CmdRecolor: called with %s.\n", str1 ENDFD;

    if((ok = APIEnterNotModal(G))) {
      if(WordMatchExact(G, str1, "all", true))
        ExecutiveInvalidateRep(G, str1, rep, cRepInvColor);
      else {
        ok = (SelectorGetTmp2(G, str1, s1) >= 0);
        ExecutiveInvalidateRep(G, s1, rep, cRepInvColor);
        SelectorFreeTmp(G, s1);
      }
      APIExit(G);
    }
  } else {
    ok = -1;                    /* special error convention */
  }
  return APIResultOk(ok);
}

static PyObject *CmdRebuild(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  char *str1;
  OrthoLineType s1;
  int ok = false;
  int rep = -1;
  ok = PyArg_ParseTuple(args, "Osi", &self, &str1, &rep);
  if(ok) {
    API_SETUP_PYMOL_GLOBALS;
    ok = (G != NULL);
  } else {
    API_HANDLE_ERROR;
  }
  if(ok) {
    PRINTFD(G, FB_CCmd)
      " CmdRebuild: called with %s.\n", str1 ENDFD;

    if((ok = APIEnterNotModal(G))) {
      if(WordMatchExact(G, str1, "all", true))
        ExecutiveRebuildAll(G);
      else {
        ok = (SelectorGetTmp2(G, str1, s1) >= 0);
        if(SettingGetGlobal_b(G, cSetting_defer_builds_mode))
          ExecutiveInvalidateRep(G, s1, rep, cRepInvPurge);
        else
          ExecutiveInvalidateRep(G, s1, rep, cRepInvAll);
        SelectorFreeTmp(G, s1);
      }
      APIExit(G);
    }
  } else {
    ok = -1;                    /* special error convention */
  }
  return APIResultOk(ok);
}

static PyObject *CmdResetRate(PyObject * self, PyObject * args)
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
    ButModeResetRate(G);
    APIExit(G);
  }
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

static PyObject *Cmd_GetGlobalCObject(PyObject * self, PyObject * args)
{
  return PyCObject_FromVoidPtr((void *) &SingletonPyMOLGlobals, NULL);
}

/*
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

static PyObject *Cmd_New(PyObject * self, PyObject * args)
{
  PyObject *result = NULL;
  PyObject *pymol = NULL;       /* pymol object instance */
  CPyMOLOptions *options = PyMOLOptions_New();

  if(options) {
    {
      int ok = true;
      PyObject *pyoptions = NULL;
      ok = PyArg_ParseTuple(args, "OO", &pymol, &pyoptions);
      if(!pyoptions) {
        options->show_splash = false;
      } else {
        PConvertOptions(options, pyoptions);
      }
      {
        CPyMOL *I = PyMOL_NewWithOptions(options);
        PyMOLGlobals *G = PyMOL_GetGlobals(I);
        if(I) {

          G->P_inst = pymol::calloc<CP_inst>(1);
          G->P_inst->obj = pymol;
          G->P_inst->dict = PyObject_GetAttrString(pymol, "__dict__");
          Py_DECREF(G->P_inst->dict); // borrow reference
          {
            /* store the PyMOL struct as a CObject */
            PyObject *tmp = PyCObject_FromVoidPtr(I, NULL);
            PyObject_SetAttrString(pymol, "__pymol__", tmp);
            Py_DECREF(tmp);
          }
          {
            int a;
            SavedThreadRec *str = G->P_inst->savedThread;
            for(a = 0; a < MAX_SAVED_THREAD; a++) {
              (str++)->id = -1;
            }
          }
          result = PyCObject_FromVoidPtr((void *) PyMOL_GetGlobalsHandle(I), NULL);
        }
      }
    }
    PyMOLOptions_Free(options);
  }
  return APIAutoNone(result);
}

static PyObject *Cmd_Del(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  int ok = true;
  ok = PyArg_ParseTuple(args, "O", &self);
  if(ok) {
    API_SETUP_PYMOL_GLOBALS;
    ok = (G != NULL);
  } else {
    API_HANDLE_ERROR;
  }
  if(ok) {
    /* leaking Px */
    PyMOL_Free(G->PyMOL);
  }
  return APIResultOk(ok);
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
  if(ok) {
    PLockAPIAndUnblock(G);
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
  if(ok) {
    PLockAPIAndUnblock(G);
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
  if(ok) {
    PLockAPIAndUnblock(G);
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
  if(ok) {
    PLockAPIAndUnblock(G);
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
  if(ok) {
    PLockAPIAndUnblock(G);
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
  if(ok) {
    PLockAPIAndUnblock(G);
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
  exit(0);
#else
#ifndef _PYMOL_WX_GLUT

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
#endif

  return APISuccess();
}

static PyObject *CmdRunWXPyMOL(PyObject * self, PyObject * args)
{
#ifdef _PYMOL_WX_GLUT
#ifndef _PYMOL_OLD_ACTIVEX
#ifndef _PYMOL_EMBEDDED
  if(run_only_once) {
    run_only_once = false;
    was_main();
  }
#endif
#endif
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
  ok = PyArg_ParseTuple(args, "Os", &self, &str1);
  if(ok) {
    API_SETUP_PYMOL_GLOBALS;
    ok = (G != NULL);
  } else {
    API_HANDLE_ERROR;
  }
  if(ok && (ok = APIEnterNotModal(G))) {
    ok = (SelectorGetTmp2(G, str1, s1) >= 0);
    count = ExecutiveCountStates(G, s1);
    if(count < 0)
      ok = false;
    SelectorFreeTmp(G, s1);
    APIExit(G);
  }
  return ok ? APIResultCode(count) : APIFailure();
}

static PyObject *CmdCountFrames(PyObject * self, PyObject * args)
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
  if(ok && (ok = APIEnterNotModal(G))) {
    SceneCountFrames(G);
    result = SceneGetNFrame(G, NULL);
    APIExit(G);
  }
  return (APIResultCode(result));
}

static PyObject *CmdGetMovieLength(PyObject * self, PyObject * args)
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
  if(ok && (ok = APIEnterNotModal(G))) {
    result = MovieGetLength(G);
    APIExit(G);
  }
  return (APIResultCode(result));
}

static PyObject *CmdIdentify(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  char *str1;
  OrthoLineType s1;
  int mode;
  int a, l = 0;
  PyObject *result = Py_None;
  PyObject *tuple;
  int *iVLA = NULL, *i;
  ObjectMolecule **oVLA = NULL, **o;
  int ok = false;
  ok = PyArg_ParseTuple(args, "Osi", &self, &str1, &mode);
  if(ok) {
    API_SETUP_PYMOL_GLOBALS;
    ok = (G != NULL);
  } else {
    API_HANDLE_ERROR;
  }
  if(ok && (ok = APIEnterNotModal(G))) {
    ok = (SelectorGetTmp(G, str1, s1) >= 0);
    if(ok) {
      if(!mode) {
        iVLA = ExecutiveIdentify(G, s1, mode);
      } else {
        l = ExecutiveIdentifyObjects(G, s1, mode, &iVLA, &oVLA);
      }
    }
    SelectorFreeTmp(G, s1);
    APIExit(G);
    if(iVLA) {
      if(!mode) {
        result = PConvIntVLAToPyList(iVLA);
      } else {                  /* object mode */
        result = PyList_New(l);
        i = iVLA;
        o = oVLA;
        for(a = 0; a < l; a++) {
          tuple = PyTuple_New(2);
          PyTuple_SetItem(tuple, 1, PyInt_FromLong(*(i++)));
          PyTuple_SetItem(tuple, 0, PyString_FromString((*(o++))->Obj.Name));
          PyList_SetItem(result, a, tuple);
        }
      }
    } else {
      result = PyList_New(0);
    }
  }
  VLAFreeP(iVLA);
  VLAFreeP(oVLA);
  if(!ok) {
    if(result && (result != Py_None)) {
      Py_DECREF(result);
    }
    return APIFailure();
  } else {
    return (APIAutoNone(result));
  }
}

static PyObject *CmdIndex(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  char *str1;
  int mode;
  PyObject *result = Py_None;
  PyObject *tuple = Py_None;
  int *iVLA = NULL;
  int l = 0;
  int *i;
  ObjectMolecule **o, **oVLA = NULL;
  int a;

  int ok = false;
  ok = PyArg_ParseTuple(args, "Osi", &self, &str1, &mode);
  if(ok) {
    API_SETUP_PYMOL_GLOBALS;
    ok = (G != NULL);
  } else {
    API_HANDLE_ERROR;
  }
  if(ok && (ok = APIEnterNotModal(G))) {
    l = ExecutiveIndex(G, str1, mode, &iVLA, &oVLA);
    APIExit(G);
    if(iVLA) {
      result = PyList_New(l);
      i = iVLA;
      o = oVLA;
      for(a = 0; a < l; a++) {
        tuple = PyTuple_New(2);
        PyTuple_SetItem(tuple, 1, PyInt_FromLong(*(i++) + 1));  /* +1 for index */
        PyTuple_SetItem(tuple, 0, PyString_FromString((*(o++))->Obj.Name));
        PyList_SetItem(result, a, tuple);
      }
    } else {
      result = PyList_New(0);
    }
    VLAFreeP(iVLA);
    VLAFreeP(oVLA);
  }
  if(!ok) {
    if(result && (result != Py_None)) {
      Py_DECREF(result);
    }
    return APIFailure();
  } else {
    return (APIAutoNone(result));
  }
}

static PyObject *CmdFindPairs(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  char *str1, *str2;
  int state1, state2;
  float cutoff;
  float angle;
  int mode;
  OrthoLineType s1, s2;
  PyObject *result = Py_None;
  PyObject *tuple = Py_None;
  PyObject *tuple1 = Py_None;
  PyObject *tuple2 = Py_None;
  int *iVLA = NULL;
  int l;
  int *i;
  ObjectMolecule **o, **oVLA = NULL;
  int a;

  int ok = false;
  ok =
    PyArg_ParseTuple(args, "Ossiiiff", &self, &str1, &str2, &state1, &state2, &mode,
                     &cutoff, &angle);
  if(ok) {
    API_SETUP_PYMOL_GLOBALS;
    ok = (G != NULL);
  } else {
    API_HANDLE_ERROR;
  }
  if(ok && (ok = APIEnterNotModal(G))) {
    ok = ((SelectorGetTmp(G, str1, s1) >= 0) && (SelectorGetTmp(G, str2, s2) >= 0));
    l =
      ExecutivePairIndices(G, s1, s2, state1, state2, mode, cutoff, angle, &iVLA, &oVLA);
    SelectorFreeTmp(G, s1);
    SelectorFreeTmp(G, s2);
    APIExit(G);

    if(iVLA && oVLA) {
      result = PyList_New(l);
      i = iVLA;
      o = oVLA;
      for(a = 0; a < l; a++) {
        tuple1 = PyTuple_New(2);
        PyTuple_SetItem(tuple1, 0, PyString_FromString((*(o++))->Obj.Name));
        PyTuple_SetItem(tuple1, 1, PyInt_FromLong(*(i++) + 1)); /* +1 for index */
        tuple2 = PyTuple_New(2);
        PyTuple_SetItem(tuple2, 0, PyString_FromString((*(o++))->Obj.Name));
        PyTuple_SetItem(tuple2, 1, PyInt_FromLong(*(i++) + 1)); /* +1 for index */
        tuple = PyTuple_New(2);
        PyTuple_SetItem(tuple, 0, tuple1);
        PyTuple_SetItem(tuple, 1, tuple2);
        PyList_SetItem(result, a, tuple);
      }
    } else {
      result = PyList_New(0);
    }
    VLAFreeP(iVLA);
    VLAFreeP(oVLA);
  }
  return (APIAutoNone(result));
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
  if(ok && (ok = APIEnterNotModal(G))) {
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

  ok_assert(1, PyArg_ParseTuple(args, "Ossisiii", &self,
        &format, &sele, &state, &ref, &ref_state, &multi, &quiet));
  API_SETUP_PYMOL_GLOBALS;
  ok_assert(1, G && APIEnterNotModal(G));

  vla = MoleculeExporterGetStr(G, format, sele, state,
      ref, ref_state, multi, quiet);

  ok_assert(2, vla);
  result = PyBytes_FromStringAndSize(vla, VLAGetSize(vla));

ok_except2:
  APIExit(G);
  return APIAutoNone(result);
ok_except1:
  API_HANDLE_ERROR;
  return APIAutoNone(NULL);
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
  ok = PyArg_ParseTuple(args, "Osisi", &self, &str1, &state, &ref_object, &ref_state);
  if(ok) {
    API_SETUP_PYMOL_GLOBALS;
    ok = (G != NULL);
  } else {
    API_HANDLE_ERROR;
  }
  if(ok) {
    if(!ref_object[0])
      ref_object = NULL;
    if((ok = APIEnterBlockedNotModal(G))) {
      ok = (SelectorGetTmp(G, str1, s1) >= 0);
      if(ok)
        result = ExecutiveSeleToChemPyModel(G, s1, state, ref_object, ref_state);
      SelectorFreeTmp(G, s1);
      APIExitBlocked(G);
    }
  }
  return (APIAutoNone(result));
}

static PyObject *CmdGetBonds(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  PyObject *result = NULL;
  char *sele;
  int state;

  ok_assert(1, PyArg_ParseTuple(args, "Osi", &self, &sele, &state));
  API_SETUP_PYMOL_GLOBALS;

  ok_assert(1, G && APIEnterNotModal(G));
  result = MoleculeExporterGetPyBonds(G, sele, state);
  APIExit(G);

  if (0) {
ok_except1:
    API_HANDLE_ERROR;
  }

  return APIAutoNone(result);
}

static PyObject *CmdCreate(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  char *str1, *str2;
  int target, source, discrete, quiet;
  int singletons;
  int copy_properties;
  int ok = false;
  int zoom;
  ok = PyArg_ParseTuple(args, "Ossiiiiii", &self, &str1, &str2, &source,
                        &target, &discrete, &zoom, &quiet, &singletons);
  if(ok) {
    API_SETUP_PYMOL_GLOBALS;
    ok = (G != NULL);
  } else {
    API_HANDLE_ERROR;
  }
  if(ok && (ok = APIEnterNotModal(G))) {
    ok = ExecutiveSeleToObject(G, str1, str2,
                                 source, target, discrete, zoom, quiet, singletons, copy_properties);
    APIExit(G);
  }
  return APIResultOk(ok);
}

static PyObject *CmdOrient(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  double m[16];
  char *str1;
  OrthoLineType s1;
  int state;
  int ok = false;
  float animate;
  int quiet = false;            /* TODO */
  ok = PyArg_ParseTuple(args, "Osif", &self, &str1, &state, &animate);
  if(ok) {
    API_SETUP_PYMOL_GLOBALS;
    ok = (G != NULL);
  } else {
    API_HANDLE_ERROR;
  }
  if(ok && (ok = APIEnterNotModal(G))) {
    ok = (SelectorGetTmp(G, str1, s1) >= 0);
    if(ExecutiveGetMoment(G, s1, m, state))
      ExecutiveOrient(G, s1, m, state, animate, false, 0.0F, quiet);    /* TODO STATUS */
    SelectorFreeTmp(G, s1);
    APIExit(G);
  }
  return APIResultOk(ok);
}

static PyObject *CmdFitPairs(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  PyObject *list;
  WordType *word = NULL;
  int ln = 0;
  int a;
  PyObject *result = NULL;
  float valu = -1.0F;
  int ok = false;
  int quiet = 0;
  ok = PyArg_ParseTuple(args, "OOi", &self, &list, &quiet);
  if(ok) {
    API_SETUP_PYMOL_GLOBALS;
    ok = (G != NULL);
  } else {
    API_HANDLE_ERROR;
  }
  if(ok && (ok = APIEnterBlockedNotModal(G))) {
    ln = PyObject_Length(list);
    if(ln) {
      if(ln & 0x1)
        ok = ErrMessage(G, "FitPairs", "must supply an even number of selections.");
    } else
      ok = false;

    if(ok) {
      word = pymol::malloc<WordType>(ln);

      a = 0;
      while(a < ln) {
        PyObject * item = PySequence_GetItem(list, a);
        SelectorGetTmp(G, PyString_AsString(item), word[a]);
        Py_DECREF(item);
        a++;
      }
      if((ok = APIEnterNotModal(G))) {
        valu = ExecutiveRMSPairs(G, word, ln / 2, 2, quiet);
        APIExit(G);
      }
      result = Py_BuildValue("f", valu);
      for(a = 0; a < ln; a++)
        SelectorFreeTmp(G, word[a]);
      FreeP(word);
    }
    APIExitBlocked(G);
  }
  return APIAutoNone(result);
}

static PyObject *CmdIntraFit(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  char *str1;
  int state;
  int mode;
  int quiet;
  int mix;
  float *fVLA = NULL;
  PyObject *result = Py_None;
  int ok = false;
  ok = PyArg_ParseTuple(args, "Osiiii", &self, &str1, &state, &mode, &quiet, &mix);
  if(ok) {
    API_SETUP_PYMOL_GLOBALS;
    ok = (G != NULL);
  } else {
    API_HANDLE_ERROR;
  }
  if(ok) {
    if(state < 0)
      state = 0;
    if((ok = APIEnterNotModal(G))) {
      fVLA = ExecutiveRMSStates(G, str1, state, mode, quiet, mix);
      APIExit(G);
    }
    if(fVLA) {
      result = PConvFloatVLAToPyList(fVLA);
      VLAFreeP(fVLA);
    }
  }
  return APIAutoNone(result);
}

static PyObject *CmdGetAtomCoords(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  char *str1;
  int state;
  int quiet;
  OrthoLineType s1;
  float vertex[3];
  PyObject *result = Py_None;
  int ok = false;
  ok = PyArg_ParseTuple(args, "Osii", &self, &str1, &state, &quiet);
  if(ok) {
    API_SETUP_PYMOL_GLOBALS;
    ok = (G != NULL);
  } else {
    API_HANDLE_ERROR;
  }
  if(ok && (ok = APIEnterNotModal(G))) {
    ok = (SelectorGetTmp(G, str1, s1) >= 0);
    if(ok)
      ok = ExecutiveGetAtomVertex(G, s1, state, quiet, vertex);
    SelectorFreeTmp(G, s1);
    APIExit(G);
    if(ok) {
      result = PConvFloatArrayToPyList(vertex, 3);
    }
  }
  return APIAutoNone(result);
}

static PyObject *CmdFit(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  char *str1, *str2;
  int mode;
  int quiet;
  OrthoLineType s1, s2;
  PyObject *result;
  float cutoff;
  int state1, state2;
  int ok = false;
  int matchmaker, cycles;
  char *object;
  ExecutiveRMSInfo rms_info;
  ok = PyArg_ParseTuple(args, "Ossiiiiifis", &self, &str1, &str2, &mode,
                        &state1, &state2, &quiet, &matchmaker, &cutoff, &cycles, &object);
  if(ok) {
    API_SETUP_PYMOL_GLOBALS;
    ok = (G != NULL);
  } else {
    API_HANDLE_ERROR;
  }
  if(ok && (ok = APIEnterNotModal(G))) {
    ok = ((SelectorGetTmp(G, str1, s1) >= 0) && (SelectorGetTmp(G, str2, s2) >= 0));
    if(ok)
      ok = ExecutiveRMS(G, s1, s2, mode,
                        cutoff, cycles, quiet, object, state1, state2,
                        false, matchmaker, &rms_info);
    SelectorFreeTmp(G, s1);
    SelectorFreeTmp(G, s2);
    APIExit(G);
  }
  if(ok) {
    result = Py_BuildValue("f", rms_info.final_rms);
  } else {
    result = Py_BuildValue("f", -1.0F);
  }
  return result;
}

static PyObject *CmdUpdate(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  char *str1, *str2;
  int int1, int2;
  int ok = false;
  int matchmaker, quiet;
  ok =
    PyArg_ParseTuple(args, "Ossiiii", &self, &str1, &str2, &int1, &int2, &matchmaker,
                     &quiet);
  if(ok) {
    API_SETUP_PYMOL_GLOBALS;
    ok = (G != NULL);
  } else {
    API_HANDLE_ERROR;
  }
  if(ok && (ok = APIEnterNotModal(G))) {
    ExecutiveUpdateCmd(G, str1, str2, int1, int2, matchmaker, quiet);       /* TODO STATUS */
    APIExit(G);
  }
  return APIResultOk(ok);
}

static PyObject *CmdDirty(PyObject * self, PyObject * args)
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
    PRINTFD(G, FB_CCmd)
      " CmdDirty: called.\n" ENDFD;
    if((ok = APIEnterNotModal(G))) {
      OrthoDirty(G);
      APIExit(G);
    }
  }
  return APISuccess();
}

static PyObject *CmdGetObjectList(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  char *str1;
  OrthoLineType s1;
  int ok = false;
  ObjectMolecule **list = NULL;
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
    list = ExecutiveGetObjectMoleculeVLA(G, s1);
    if(list) {
      unsigned int size = VLAGetSize(list);
      result = PyList_New(size);
      if(result) {
        unsigned int a;
        for(a = 0; a < size; a++) {
          PyList_SetItem(result, a, PyString_FromString(list[a]->Obj.Name));
        }
      }
      VLAFreeP(list);
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
  float result;
  int int1;
  int ok = false;
  ok = PyArg_ParseTuple(args, "Ossi", &self, &str1, &str2, &int1);

  if(ok) {
    API_SETUP_PYMOL_GLOBALS;
    ok = (G != NULL);
  } else {
    API_HANDLE_ERROR;
  }
  if(ok && (ok = APIEnterNotModal(G))) {
    ok = ExecutiveGetDistance(G, str1, str2, &result, int1);
    APIExit(G);
  }

  if(ok) {
    return (Py_BuildValue("f", result));
  } else {
    return APIFailure();
  }
}

static PyObject *CmdGetAngle(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  char *str1, *str2, *str3;
  float result;
  int int1;
  int ok = false;
  ok = PyArg_ParseTuple(args, "Osssi", &self, &str1, &str2, &str3, &int1);

  if(ok) {
    API_SETUP_PYMOL_GLOBALS;
    ok = (G != NULL);
  } else {
    API_HANDLE_ERROR;
  }
  if(ok && (ok = APIEnterNotModal(G))) {
    ok = ExecutiveGetAngle(G, str1, str2, str3, &result, int1);
    APIExit(G);
  }

  if(ok) {
    return (Py_BuildValue("f", result));
  } else {
    return APIFailure();
  }
}

static PyObject *CmdGetDihe(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  char *str1, *str2, *str3, *str4;
  float result;
  int int1;
  int ok = false;
  ok = PyArg_ParseTuple(args, "Ossssi", &self, &str1, &str2, &str3, &str4, &int1);

  if(ok) {
    API_SETUP_PYMOL_GLOBALS;
    ok = (G != NULL);
  } else {
    API_HANDLE_ERROR;
  }
  if(ok && (ok = APIEnterNotModal(G))) {
    ok = ExecutiveGetDihe(G, str1, str2, str3, str4, &result, int1);
    APIExit(G);
  }

  if(ok) {
    return (Py_BuildValue("f", result));
  } else {
    return APIFailure();
  }
}

static PyObject *CmdSetDihe(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  char *str1, *str2, *str3, *str4;
  float float1;
  int int1;
  int quiet;
  int ok = false;
  ok =
    PyArg_ParseTuple(args, "Ossssfii", &self, &str1, &str2, &str3, &str4, &float1, &int1,
                     &quiet);
  if(ok) {
    API_SETUP_PYMOL_GLOBALS;
    ok = (G != NULL);
  } else {
    API_HANDLE_ERROR;
  }
  if(ok && (ok = APIEnterNotModal(G))) {
    ok = ExecutiveSetDihe(G, str1, str2, str3, str4, float1, int1, quiet);
    APIExit(G);
  }
  return APIResultOk(ok);
}

static PyObject *CmdDo(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  char *str1;
  int log;
  int ok = false;
  int echo;
  ok = PyArg_ParseTuple(args, "Osii", &self, &str1, &log, &echo);
  if(ok) {
    API_SETUP_PYMOL_GLOBALS;
    ok = (G != NULL);
  } else {
    API_HANDLE_ERROR;
  }
  if(ok && (ok = APIEnterNotModal(G))) {
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
  }
  return APIResultOk(ok);
}

static PyObject *CmdRock(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  int int1;
  int result = -1;
  int ok = false;

  ok = PyArg_ParseTuple(args, "Oi", &self, &int1);
  if(ok) {
    API_SETUP_PYMOL_GLOBALS;
    ok = (G != NULL);
  } else {
    API_HANDLE_ERROR;
  }
  if(ok && (ok = APIEnterNotModal(G))) {
    result = ControlRock(G, int1);
    APIExit(G);
  }
  return APIResultCode(result);
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
  if(ok && (ok = APIEnterNotModal(G))) {
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
  int index, type = -1;
  if (PyArg_ParseTuple(args, "i", &index)) {
    type = SettingGetType(index);
  }
  return PyInt_FromLong(type);
}

static PyObject *CmdGetSettingTuple(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  PyObject *result = Py_None;
  int int1, int2;
  char *str1;
  int ok = false;
  ok = PyArg_ParseTuple(args, "Oisi", &self, &int1, &str1, &int2);      /* setting, object, state */
  if(ok) {
    API_SETUP_PYMOL_GLOBALS;
    ok = (G != NULL);
  } else {
    API_HANDLE_ERROR;
  }
  if(ok) {
    APIEnterBlocked(G);
    result = ExecutiveGetSettingTuple(G, int1, str1, int2);
    APIExitBlocked(G);
  }
  return APIAutoNone(result);
}

static PyObject *CmdGetSettingOfType(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  PyObject *result = Py_None;
  int int1, int2, int3;
  char *str1;
  int ok = false;
  ok = PyArg_ParseTuple(args, "Oisii", &self, &int1, &str1, &int2, &int3);      /* setting, object, state */
  if(ok) {
    API_SETUP_PYMOL_GLOBALS;
    ok = (G != NULL);
  } else {
    API_HANDLE_ERROR;
  }
  if(ok) {
    APIEnterBlocked(G);
    result = ExecutiveGetSettingOfType(G, int1, str1, int2, int3);
    APIExitBlocked(G);
  }
  return APIAutoNone(result);
}

static PyObject *CmdGetSettingText(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  PyObject *result = Py_None;
  int int1, int2;
  char *str1;
  int ok = false;
  ok = PyArg_ParseTuple(args, "Oisi", &self, &int1, &str1, &int2);      /* setting, object, state */
  if(ok) {
    API_SETUP_PYMOL_GLOBALS;
    ok = (G != NULL);
  } else {
    API_HANDLE_ERROR;
  }
  if(ok) {
    APIEnterBlocked(G);
    result = ExecutiveGetSettingText(G, int1, str1, int2);
    APIExitBlocked(G);
  }
  return APIAutoNone(result);
}

static PyObject *CmdExportDots(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  PyObject *result = NULL;
  PyObject *cObj;
  ExportDotsObj *obj;
  char *str1;
  int int1;

  int ok = false;
  ok = PyArg_ParseTuple(args, "Osi", &self, &str1, &int1);
  if(ok) {
    API_SETUP_PYMOL_GLOBALS;
    ok = (G != NULL);
  } else {
    API_HANDLE_ERROR;
  }
  if(ok && (ok = APIEnterNotModal(G))) {
    obj = ExportDots(G, str1, int1);
    APIExit(G);
    if(obj) {
      cObj = PyCObject_FromVoidPtr(obj, (void (*)(void *)) ExportDeleteMDebug);
      if(cObj) {
        result = Py_BuildValue("O", cObj);      /* I think this */
        Py_DECREF(cObj);        /* transformation is unnecc. */
      }
    }
  }
  return APIAutoNone(result);
}

static PyObject *CmdSetFrame(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  int mode, frm;
  int ok = false;
  ok = PyArg_ParseTuple(args, "Oii", &self, &mode, &frm);
  if(ok) {
    API_SETUP_PYMOL_GLOBALS;
    ok = (G != NULL);
  } else {
    API_HANDLE_ERROR;
  }
  if(ok && (ok = APIEnterNotModal(G))) {
    SceneSetFrame(G, mode, frm);
    APIExit(G);
  }
  return APIResultOk(ok);
}

static PyObject *CmdFrame(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  int frm,trigger;
  int ok = false;
  ok = PyArg_ParseTuple(args, "Oii", &self, &frm, &trigger);
  if(ok) {
    API_SETUP_PYMOL_GLOBALS;
    ok = (G != NULL);
  } else {
    API_HANDLE_ERROR;
  }
  if(ok) {
    if(frm < 0)
      frm = 0;
    if((ok = APIEnterNotModal(G))) {
      SceneSetFrame(G, trigger ? 4 : 0, frm);
      APIExit(G);
    }
  }
  return APIResultOk(ok);
}

static PyObject *CmdStereo(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  int i1;
  int ok = false;

  ok = PyArg_ParseTuple(args, "Oi", &self, &i1);
  if(ok) {
    API_SETUP_PYMOL_GLOBALS;
    ok = (G != NULL);
  } else {
    API_HANDLE_ERROR;
  }
  if(ok && (ok = APIEnterNotModal(G))) {
    ok = ExecutiveStereo(G, i1);
    APIExit(G);
  }
  return APIResultOk(ok);
}

static PyObject *CmdReset(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  int cmd;
  int ok = false;
  char *obj;
  ok = PyArg_ParseTuple(args, "Ois", &self, &cmd, &obj);
  if(ok) {
    API_SETUP_PYMOL_GLOBALS;
    ok = (G != NULL);
  } else {
    API_HANDLE_ERROR;
  }
  if(ok && (ok = APIEnterNotModal(G))) {
    ok = ExecutiveReset(G, cmd, obj);
    APIExit(G);
  }
  return APIResultOk(ok);
}

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

static PyObject *CmdGetMinMax(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  float mn[3], mx[3];
  char *str1;
  int state;
  OrthoLineType s1;
  PyObject *result = Py_None;
  int flag;

  int ok = false;
  ok = PyArg_ParseTuple(args, "Osi", &self, &str1, &state);
  if(ok) {
    API_SETUP_PYMOL_GLOBALS;
    ok = (G != NULL);
  } else {
    API_HANDLE_ERROR;
  }
  if(ok && (ok = APIEnterNotModal(G))) {
    ok = (SelectorGetTmp2(G, str1, s1) >= 0);
    flag = ExecutiveGetExtent(G, s1, mn, mx, true, state, false);
    SelectorFreeTmp(G, s1);
    APIExit(G);
    if(flag)
      result = Py_BuildValue("[[fff],[fff]]", mn[0], mn[1], mn[2], mx[0], mx[1], mx[2]);
    else
      result = Py_BuildValue("[[fff],[fff]]", -0.5, -0.5, -0.5, 0.5, 0.5, 0.5);

  }
  return APIAutoNone(result);
}

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
  if(ok && (ok = APIEnterNotModal(G))) {
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
  ok_assert(1, G && APIEnterNotModal(G));

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

static PyObject *CmdPNG(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  char *str1;
  int ok = false;
  int quiet;
  int result = 0;
  int width, height, ray;
  int prior, format;
  float dpi;
  ok =
    PyArg_ParseTuple(args, "Osiifiiii", &self, &str1, &width, &height, &dpi, &ray, &quiet,
                     &prior, &format);
  if(ok) {
    API_SETUP_PYMOL_GLOBALS;
    ok = (G != NULL);
  } else {
    API_HANDLE_ERROR;
  }
  if(ok && (ok = APIEnterNotModal(G))) {
    // with prior=1 other arguments (width, height, ray) are ignored

    if(!prior) {
      if(ray || (!G->HaveGUI && (!SceneGetCopyType(G) || width || height))) {
        prior = SceneRay(G, width, height, SettingGetGlobal_i(G, cSetting_ray_default_renderer),
                 NULL, NULL, 0.0F, 0.0F, false, NULL, true, -1);
      } else if(width || height) {
        SceneDeferImage(G, width, height, str1, -1, dpi, quiet, format);
        result = 1;
      } else if(!SceneGetCopyType(G)) {
        ExecutiveDrawNow(G);      /* TODO STATUS */
      }
    }

    if(!result) {
      PyMOL_PushValidContext(G->PyMOL); // PyQt hack?
      if(ScenePNG(G, str1, dpi, quiet, prior, format))
        result = 1;             /* signal success by returning 1 instead of 0, or -1 for error  */
      PyMOL_PopValidContext(G->PyMOL);
    }
    APIExit(G);
  }
  if(!ok)
    result = -1;
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
    PyMOL_PushValidContext(G->PyMOL); // PyQt hack?
    ok = MoviePNG(G, str1, SettingGetGlobal_b(G, cSetting_cache_frames),
                  int1, int2, int3, int4, format, mode, quiet,
                  width, height);
    PyMOL_PopValidContext(G->PyMOL);
    /* TODO STATUS */
    APIExit(G);
  }
  return APIResultOk(ok);
}

static PyObject *CmdMSet(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  char *str1;
  int ok = false;
  int start_from,freeze;
  ok = PyArg_ParseTuple(args, "Osii", &self, &str1, &start_from,&freeze);
  if(ok) {
    API_SETUP_PYMOL_GLOBALS;
    ok = (G != NULL);
  } else {
    API_HANDLE_ERROR;
  }
  if(ok && (ok = APIEnterNotModal(G))) {
    MovieAppendSequence(G, str1, start_from,freeze);
    SceneCountFrames(G);
    APIExit(G);
  }

  // fix for PYMOL-1465
  // force GUI update for movie panel
  if(G->HaveGUI)
  OrthoReshape(G, -1, -1, false);

  return APIResultOk(ok);
}

static PyObject *CmdMModify(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  char *object;
  int ok = false;
  int action,index,count,target,freeze,quiet;
  ok = PyArg_ParseTuple(args, "Oiiiisii",  &self, &action, &index,
                        &count, &target, &object, &freeze, &quiet);
  if(ok) {
    API_SETUP_PYMOL_GLOBALS;
    ok = (G != NULL);
  } else {
    API_HANDLE_ERROR;
  }
  if(ok && (ok = APIEnterNotModal(G))) {
    ExecutiveMotionViewModify(G,action,index,count,target,object,freeze,quiet);
    SceneCountFrames(G);
    APIExit(G);
  }
  return APIResultOk(ok);
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
  int ok = false;
  ok = PyArg_ParseTuple(args, "Oii", &self, &w, &h);
  if(ok) {
    API_SETUP_PYMOL_GLOBALS;
    ok = (G != NULL);
  } else {
    API_HANDLE_ERROR;
  }
  if(ok) {
    if((ok = APIEnterNotModal(G))) {
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
  return APIResultOk(ok);
}

static PyObject *CmdFlag(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  char *str1;
  int flag;
  int action;
  int quiet;
  OrthoLineType s1;
  int ok = false;
  ok = PyArg_ParseTuple(args, "Oisii", &self, &flag, &str1, &action, &quiet);
  if(ok) {
    API_SETUP_PYMOL_GLOBALS;
    ok = (G != NULL);
  } else {
    API_HANDLE_ERROR;
  }
  if(ok && (ok = APIEnterNotModal(G))) {
    ok = (SelectorGetTmp(G, str1, s1) >= 0);
    ExecutiveFlag(G, flag, s1, action, quiet);
    SelectorFreeTmp(G, s1);
    APIExit(G);
  }
  return APIResultOk(ok);
}

static PyObject *CmdColor(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  char *str1, *color;
  int flags;
  OrthoLineType s1;
  int ok = false;
  int quiet;

  ok = PyArg_ParseTuple(args, "Ossii", &self, &color, &str1, &flags, &quiet);
  if(ok) {
    API_SETUP_PYMOL_GLOBALS;
    ok = (G != NULL);
  } else {
    API_HANDLE_ERROR;
  }
  if(ok && (ok = APIEnterNotModal(G))) {
    ok = (SelectorGetTmp2(G, str1, s1) >= 0);
    if(ok) {
      ok = ExecutiveColor(G, s1, color, flags, quiet);
    }
    SelectorFreeTmp(G, s1);
    APIExit(G);
  }
  return APIResultOk(ok);
}

static PyObject *CmdColorDef(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  char *color;
  float v[3];
  int ok = false;
  int mode;
  int quiet;
  ok = PyArg_ParseTuple(args, "Osfffii", &self, &color, v, v + 1, v + 2, &mode, &quiet);
  if(ok) {
    API_SETUP_PYMOL_GLOBALS;
    ok = (G != NULL);
  } else {
    API_HANDLE_ERROR;
  }
  if(ok && (ok = APIEnterNotModal(G))) {
    ColorDef(G, color, v, mode, quiet);
    APIExit(G);
  }
  return APIResultOk(ok);
}

static PyObject *CmdDraw(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  int int1, int2;
  int quiet, antialias;
  int ok = false;

  ok = PyArg_ParseTuple(args, "Oiiii", &self, &int1, &int2, &antialias, &quiet);

  if(ok) {
    API_SETUP_PYMOL_GLOBALS;
    ok = (G != NULL);
  } else {
    API_HANDLE_ERROR;
  }
  if(ok && (ok = APIEnterNotModal(G))) {
    if(antialias == -2) {
      ok = ExecutiveDrawCmd(G, 0, 0, 0, true, quiet);   /* capture action */
    } else {

      ok = ExecutiveDrawCmd(G, int1, int2, antialias, false, quiet);
    }
    APIExit(G);
  }
  return APIResultOk(ok);
}

static PyObject *CmdRay(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  int w, h, mode;
  float angle, shift;
  int ok = false;
  int quiet;
  int antialias;
  ok = PyArg_ParseTuple(args, "Oiiiffii", &self, &w, &h,
                        &antialias, &angle, &shift, &mode, &quiet);
  if(ok) {
    API_SETUP_PYMOL_GLOBALS;
    ok = (G != NULL);
  } else {
    API_HANDLE_ERROR;
  }
  if(ok && (ok = APIEnterNotModal(G))) {
    if(mode < 0)
      mode = SettingGetGlobal_i(G, cSetting_ray_default_renderer);
    ExecutiveRay(G, w, h, mode, angle, shift, quiet, false, antialias); /* TODO STATUS */
    APIExit(G);
  }
  return APIResultOk(ok);
}

static PyObject *CmdClip(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  char *sname;
  float dist;
  char *str1;
  int state;
  OrthoLineType s1;
  int ok = false;
  ok = PyArg_ParseTuple(args, "Osfsi", &self, &sname, &dist, &str1, &state);
  if(ok) {
    API_SETUP_PYMOL_GLOBALS;
    ok = (G != NULL);
  } else {
    API_HANDLE_ERROR;
  }
  if(ok && (ok = APIEnterNotModal(G))) {
    ok = (SelectorGetTmp2(G, str1, s1) >= 0);
    switch (sname[0]) {         /* TODO STATUS */
    case 'N':
    case 'n':
      SceneClip(G, 0, dist, s1, state);
      break;
    case 'f':
    case 'F':
      SceneClip(G, 1, dist, s1, state);
      break;
    case 'm':
    case 'M':
      SceneClip(G, 2, dist, s1, state);
      break;
    case 's':
    case 'S':
      SceneClip(G, 3, dist, s1, state);
      break;
    case 'a':
    case 'A':
      SceneClip(G, 4, dist, s1, state);
      break;
    }
    SelectorFreeTmp(G, s1);
    APIExit(G);
  }
  return APIResultOk(ok);

}

static PyObject *CmdMove(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  char *sname;
  float dist;
  int ok = false;
  ok = PyArg_ParseTuple(args, "Osf", &self, &sname, &dist);
  if(ok) {
    API_SETUP_PYMOL_GLOBALS;
    ok = (G != NULL);
  } else {
    API_HANDLE_ERROR;
  }
  if(ok && (ok = APIEnterNotModal(G))) {
    switch (sname[0]) {         /* TODO STATUS */
    case 'x':
      SceneTranslate(G, dist, 0.0, 0.0);
      break;
    case 'y':
      SceneTranslate(G, 0.0, dist, 0.0);
      break;
    case 'z':
      SceneTranslate(G, 0.0, 0.0, dist);
      break;
    }
    APIExit(G);
  }
  return APIResultOk(ok);
}

static PyObject *CmdTurn(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  char *sname;
  float angle;
  int ok = false;
  ok = PyArg_ParseTuple(args, "Osf", &self, &sname, &angle);
  if(ok) {
    API_SETUP_PYMOL_GLOBALS;
    ok = (G != NULL);
  } else {
    API_HANDLE_ERROR;
  }
  if(ok && (ok = APIEnterNotModal(G))) {
    switch (sname[0]) {         /* TODO STATUS */
    case 'x':
      SceneRotate(G, angle, 1.0, 0.0, 0.0);
      break;
    case 'y':
      SceneRotate(G, angle, 0.0, 1.0, 0.0);
      break;
    case 'z':
      SceneRotate(G, angle, 0.0, 0.0, 1.0);
      break;
    }
    APIExit(G);
  }
  return APIResultOk(ok);
}

static PyObject *CmdUnset(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  int index;
  int tmpFlag = false;
  char *str3;
  int state;
  int quiet;
  int updates;
  OrthoLineType s1;
  int ok = false;
  ok = PyArg_ParseTuple(args, "Oisiii", &self, &index, &str3, &state, &quiet, &updates);
  if(ok) {
    API_SETUP_PYMOL_GLOBALS;
    ok = (G != NULL);
  } else {
    API_HANDLE_ERROR;
  }
  if(ok && (ok = APIEnterNotModal(G))) {
    s1[0] = 0;
    if(!strcmp(str3, "all")) {
      strcpy(s1, str3);
    } else if(str3[0] != 0) {
      tmpFlag = true;
      ok = (SelectorGetTmp2(G, str3, s1) >= 0);
    }
    if(ok)
      ok = ExecutiveUnsetSetting(G, index, s1, state, quiet, updates);
    if(tmpFlag)
      SelectorFreeTmp(G, s1);
    APIExit(G);
  }
  return APIResultOk(ok);
}

static PyObject *CmdUnsetBond(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  int index;
  char *str3, *str4;
  int state;
  int quiet;
  int updates;
  OrthoLineType s1, s2;
  int ok = false;
  ok =
    PyArg_ParseTuple(args, "Oissiii", &self, &index, &str3, &str4, &state, &quiet,
                     &updates);
  if(ok) {
    API_SETUP_PYMOL_GLOBALS;
    ok = (G != NULL);
  } else {
    API_HANDLE_ERROR;
  }
  if(ok && (ok = APIEnterNotModal(G))) {
    s1[0] = 0;
    s2[0] = 0;
    ok = (SelectorGetTmp(G, str3, s1) >= 0);
    ok = (SelectorGetTmp(G, str4, s2) >= 0) && ok;
    if(ok)
      ok = ExecutiveUnsetBondSetting(G, index, s1, s2, state, quiet, updates);
    SelectorFreeTmp(G, s1);
    SelectorFreeTmp(G, s2);
    APIExit(G);
  }
  return APIResultOk(ok);
}

static PyObject *CmdSet(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  int index;
  int tmpFlag = false;
  PyObject *value;
  char *str3;
  int state;
  int quiet;
  int updates;
  OrthoLineType s1;
  int ok = false;
  ok =
    PyArg_ParseTuple(args, "OiOsiii", &self, &index, &value, &str3, &state, &quiet,
                     &updates);
  if(ok) {
    API_SETUP_PYMOL_GLOBALS;
    ok = (G != NULL);
  } else {
    API_HANDLE_ERROR;
  }
  if(ok && (ok = APIEnterNotModal(G))) {
    s1[0] = 0;
    if(!strcmp(str3, "all")) {
      strcpy(s1, str3);
    } else if(str3[0] != 0) {
      tmpFlag = true;
      ok = (SelectorGetTmp2(G, str3, s1) >= 0);
    }
    if(ok)
      ok = ExecutiveSetSetting(G, index, value, s1, state, quiet, updates);
    if(tmpFlag)
      SelectorFreeTmp(G, s1);
    APIExit(G);
  }
  return APIResultOk(ok);
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
  OrthoLineType s1, s2;
  int ok = false;
  ok =
    PyArg_ParseTuple(args, "OiOssiii", &self, &index, &value, &str3, &str4, &state,
                     &quiet, &updates);
  if(ok) {
    API_SETUP_PYMOL_GLOBALS;
    ok = (G != NULL);
  } else {
    API_HANDLE_ERROR;
  }
  if(ok && (ok = APIEnterNotModal(G))) {
    s1[0] = 0;
    s2[0] = 0;
    ok = (SelectorGetTmp(G, str3, s1) >= 0);
    ok = (SelectorGetTmp(G, str4, s2) >= 0) && ok;
    if(ok)
      ok = ExecutiveSetBondSetting(G, index, value, s1, s2, state, quiet, updates);
    SelectorFreeTmp(G, s1);
    SelectorFreeTmp(G, s2);
    APIExit(G);
  }
  return APIResultOk(ok);
}

static PyObject *CmdGetBond(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  PyObject *result = Py_None;
  int index;
  char *str3, *str4;
  int state;
  int quiet;
  int updates;
  OrthoLineType s1, s2;
  int ok = false;
  ok =
    PyArg_ParseTuple(args, "Oissiii", &self, &index, &str3, &str4, &state,
                     &quiet, &updates);
  if(ok) {
    API_SETUP_PYMOL_GLOBALS;
    ok = (G != NULL);
  } else {
    API_HANDLE_ERROR;
  }
  if(ok && (ok = APIEnterNotModal(G))) {
    s1[0] = 0;
    s2[0] = 0;
    ok = (SelectorGetTmp(G, str3, s1) >= 0);
    ok = (SelectorGetTmp(G, str4, s2) >= 0) && ok;
    if(ok)
      result = ExecutiveGetBondSetting(G, index, s1, s2, state, quiet, updates);
    SelectorFreeTmp(G, s1);
    SelectorFreeTmp(G, s2);
    APIExit(G);
  }
  return APIAutoNone(result);
}

static PyObject *CmdDelete(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  char *sname;

  int ok = false;
  ok = PyArg_ParseTuple(args, "Os", &self, &sname);
  if(ok) {
    API_SETUP_PYMOL_GLOBALS;
    ok = (G != NULL);
  } else {
    API_HANDLE_ERROR;
  }
  if(ok && (ok = APIEnterNotModal(G))) {
    ExecutiveDelete(G, sname);  /* TODO STATUS */
    APIExit(G);
  }
  return APIResultOk(ok);
}

static PyObject *CmdCartoon(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  char *sname;
  int type;
  OrthoLineType s1;
  int ok = false;
  ok = PyArg_ParseTuple(args, "Osi", &self, &sname, &type);
  if(ok) {
    API_SETUP_PYMOL_GLOBALS;
    ok = (G != NULL);
  } else {
    API_HANDLE_ERROR;
  }
  if(ok && (ok = APIEnterNotModal(G))) {
    ok = (SelectorGetTmp(G, sname, s1) >= 0);
    if(ok)
      ExecutiveCartoon(G, type, s1);    /* TODO STATUS */
    SelectorFreeTmp(G, s1);
    APIExit(G);
  }
  return APIResultOk(ok);
}

static PyObject *CmdShowHide(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  const char *sname;
  int rep;
  int state;
  OrthoLineType s1;
  int ok = false;
  ok = PyArg_ParseTuple(args, "Osii", &self, &sname, &rep, &state);
  if(ok) {
    API_SETUP_PYMOL_GLOBALS;
    ok = (G != NULL);
  } else {
    API_HANDLE_ERROR;
  }
  if(ok && (ok = APIEnterNotModal(G))) {        /* TODO STATUS */
    if(sname[0] == '@') {
      // DEPRECATED
      sname = cKeywordAll;
      rep = cRepBitmask;
    }
    {
      ok = (SelectorGetTmp2(G, sname, s1) >= 0);
      ExecutiveSetRepVisMask(G, s1, rep, state);
      SelectorFreeTmp(G, s1);
    }
    APIExit(G);
  }
  return APIResultOk(ok);
}

static PyObject *CmdOnOffBySele(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  char *sname;
  int onoff;
  OrthoLineType s1;
  int ok = false;
  ok = PyArg_ParseTuple(args, "Osi", &self, &sname, &onoff);
  if(ok) {
    API_SETUP_PYMOL_GLOBALS;
    ok = (G != NULL);
  } else {
    API_HANDLE_ERROR;
  }
  if(ok && (ok = APIEnterNotModal(G))) {        /* TODO STATUS */
    ok = (SelectorGetTmp2(G, sname, s1) >= 0);
    if(ok)
      ok = ExecutiveSetOnOffBySele(G, s1, onoff);
    SelectorFreeTmp(G, s1);
    APIExit(G);
  }
  return APIResultOk(ok);
}

static PyObject *CmdOnOff(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  char *name;
  int state;
  int parents = 0;
  int ok = false;
  ok = PyArg_ParseTuple(args, "Osii", &self, &name, &state, &parents);
  if(ok) {
    API_SETUP_PYMOL_GLOBALS;
    ok = (G != NULL);
  } else {
    API_HANDLE_ERROR;
  }
  if(ok && (ok = APIEnterNotModal(G))) {        /* TODO STATUS */
    ExecutiveSetObjVisib(G, name, state, parents);
    APIExit(G);
  }
  return APIResultOk(ok);
}

static PyObject *CmdToggle(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  char *sname;
  int rep;
  int ok = false;
  ok = PyArg_ParseTuple(args, "Osi", &self, &sname, &rep);
  if(ok) {
    API_SETUP_PYMOL_GLOBALS;
    ok = (G != NULL);
  } else {
    API_HANDLE_ERROR;
  }
  if(ok && (ok = APIEnterNotModal(G))) {
    ok = ExecutiveToggleRepVisib(G, sname, rep);
    APIExit(G);
  }
  return APIResultOk(ok);
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
  ok = PyArg_ParseTuple(args, "Oi", &self, &flag);
  if(ok) {
    API_SETUP_PYMOL_GLOBALS;
    ok = (G != NULL);
  } else {
    API_HANDLE_ERROR;
  }
  if(ok && (ok = APIEnterNotModal(G))) {
    ExecutiveFullScreen(G, flag);       /* TODO STATUS */
    APIExit(G);
  }
  return APIResultOk(ok);
}

static PyObject *CmdGroup(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  char *gname, *names;
  int quiet, action;
  int ok = false;
  ok = PyArg_ParseTuple(args, "Ossii", &self, &gname, &names, &action, &quiet);
  if(ok) {
    API_SETUP_PYMOL_GLOBALS;
    ok = (G != NULL);
  } else {
    API_HANDLE_ERROR;
  }
  if(ok && (ok = APIEnterNotModal(G))) {
    ok = ExecutiveGroup(G, gname, names, action, quiet);
    APIExit(G);
  }
  return APIResultOk(ok);
}

static PyObject *CmdSelect(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  char *sname, *sele;
  int quiet;
  int ok = false;
  int count = 0;
  int state = 0;
  char *domain;
  ok = PyArg_ParseTuple(args, "Ossiis", &self, &sname, &sele, &quiet, &state, &domain);
  if(ok) {
    API_SETUP_PYMOL_GLOBALS;
    ok = (G != NULL);
  } else {
    API_HANDLE_ERROR;
  }
  if(ok && (ok = APIEnterNotModal(G))) {
    if(!domain[0])
      domain = NULL;
    if(ExecutiveFindObjectByName(G, sname)) {   /* name conflicts with an object */
      count = -1;
    } else {
      count =
        SelectorCreateWithStateDomain(G, sname, sele, NULL, quiet, NULL, state, domain);
    }
    if(count < 0)
      ok = false;
    SceneInvalidate(G);
    SeqDirty(G);
    APIExit(G);
  }
  return ok ? APIResultCode(count) : APIFailure();
}

static PyObject *CmdFinishObject(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  char *oname;
  CObject *origObj = NULL;

  int ok = false;
  ok = PyArg_ParseTuple(args, "Os", &self, &oname);

  if(ok) {
    API_SETUP_PYMOL_GLOBALS;
    ok = (G != NULL);
  } else {
    API_HANDLE_ERROR;
  }
  if(ok && (ok = APIEnterNotModal(G))) {
    origObj = ExecutiveFindObjectByName(G, oname);
    if(origObj) {
      if(origObj->type == cObjectMolecule) {
        ObjectMoleculeUpdateIDNumbers((ObjectMolecule *) origObj);
        ObjectMoleculeUpdateNonbonded((ObjectMolecule *) origObj);
        ObjectMoleculeInvalidate((ObjectMolecule *) origObj, cRepAll, cRepInvAll, -1);
      }
      ExecutiveUpdateObjectSelection(G, origObj);       /* TODO STATUS */
    } else
      ok = false;
    APIExit(G);
  }
  return APIResultOk(ok);
}

static PyObject *CmdLoadObject(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  char *oname;
  PyObject *model;
  CObject *origObj = NULL, *obj;
  OrthoLineType buf;
  int frame, type;
  int finish, discrete;
  int quiet;
  int ok = false;
  int zoom;
  ok = PyArg_ParseTuple(args, "OsOiiiiii", &self, &oname, &model, &frame, &type,
                        &finish, &discrete, &quiet, &zoom);
  if(ok) {
    API_SETUP_PYMOL_GLOBALS;
    ok = (G != NULL);
  } else {
    API_HANDLE_ERROR;
  }
  if(ok && (ok = APIEnterNotModal(G))) {
    ObjectNameType valid_name = "";
  
    buf[0] = 0;
    ExecutiveProcessObjectName(G, oname, valid_name);

    origObj = ExecutiveFindObjectByName(G, valid_name);

    /* TODO check for existing object of wrong type */

    switch (type) {
    case cLoadTypeChemPyModel:
      if(origObj){
        if(origObj->type != cObjectMolecule) {
          ExecutiveDelete(G, valid_name);
          origObj = NULL;
        } else {
	  discrete = 1;
	}
      }
      PBlock(G);                /*PBlockAndUnlockAPI(); */
      obj = (CObject *) ObjectMoleculeLoadChemPyModel(G, (ObjectMolecule *)
                                                      origObj, model, frame, discrete);
      PUnblock(G);              /*PLockAPIAndUnblock(); */
      if(!origObj) {
        if(obj) {
          ObjectSetName(obj, valid_name);
          ExecutiveManageObject(G, obj, zoom, quiet);
          if(frame < 0)
            frame = ((ObjectMolecule *) obj)->NCSet - 1;
          sprintf(buf, " CmdLoad: ChemPy-model loaded into object \"%s\", state %d.\n",
                  valid_name, frame + 1);
        }
      } else if(origObj) {
        if(finish)
          ExecutiveUpdateObjectSelection(G, origObj);
        if(frame < 0)
          frame = ((ObjectMolecule *) origObj)->NCSet - 1;
        sprintf(buf, " CmdLoad: ChemPy-model appended into object \"%s\", state %d.\n",
                valid_name, frame + 1);
      }
      break;
    case cLoadTypeChemPyBrick:
      if(origObj)
        if(origObj->type != cObjectMap) {
          ExecutiveDelete(G, valid_name);
          origObj = NULL;
        }
      PBlock(G);                /*PBlockAndUnlockAPI(); */
      obj =
        (CObject *) ObjectMapLoadChemPyBrick(G, (ObjectMap *) origObj, model, frame,
                                             discrete, quiet);
      PUnblock(G);              /*PLockAPIAndUnblock(); */
      if(!origObj) {
        if(obj) {
          ObjectSetName(obj, valid_name);
          ExecutiveManageObject(G, obj, zoom, quiet);
          sprintf(buf, " CmdLoad: chempy.brick loaded into object \"%s\"\n", valid_name);
        }
      } else if(origObj) {
        sprintf(buf, " CmdLoad: chempy.brick appended into object \"%s\"\n", valid_name);
      }
      break;
    case cLoadTypeChemPyMap:
      if(origObj)
        if(origObj->type != cObjectMap) {
          ExecutiveDelete(G, valid_name);
          origObj = NULL;
        }
      PBlock(G);                /*PBlockAndUnlockAPI(); */
      obj =
        (CObject *) ObjectMapLoadChemPyMap(G, (ObjectMap *) origObj, model, frame,
                                           discrete, quiet);
      PUnblock(G);              /*PLockAPIAndUnblock(); */
      if(!origObj) {
        if(obj) {
          ObjectSetName(obj, valid_name);
          ExecutiveManageObject(G, obj, zoom, quiet);
          sprintf(buf, " CmdLoad: chempy.map loaded into object \"%s\"\n", valid_name);
        }
      } else if(origObj) {
        sprintf(buf, " CmdLoad: chempy.map appended into object \"%s\"\n", valid_name);
      }
      break;
    case cLoadTypeCallback:
      if(origObj)
        if(origObj->type != cObjectCallback) {
          ExecutiveDelete(G, valid_name);
          origObj = NULL;
        }
      PBlock(G);                /*PBlockAndUnlockAPI(); */
      obj = (CObject *) ObjectCallbackDefine(G, (ObjectCallback *) origObj, model, frame);
      PUnblock(G);              /*PLockAPIAndUnblock(); */
      if(!origObj) {
        if(obj) {
          ObjectSetName(obj, valid_name);
          ExecutiveManageObject(G, obj, zoom, quiet);
          sprintf(buf, " CmdLoad: pymol.callback loaded into object \"%s\"\n",
                  valid_name);
        }
      } else if(origObj) {
        sprintf(buf, " CmdLoad: pymol.callback appended into object \"%s\"\n",
                valid_name);
      }
      break;
    case cLoadTypeCGO:
      if(origObj)
        if(origObj->type != cObjectCGO) {
          ExecutiveDelete(G, valid_name);
          origObj = NULL;
        }
      PBlock(G);                /*PBlockAndUnlockAPI(); */
      obj = (CObject *) ObjectCGODefine(G, (ObjectCGO *) origObj, model, frame);
      PUnblock(G);              /*PLockAPIAndUnblock(); */
      if(!origObj) {
        if(obj) {
          ObjectSetName(obj, valid_name);
          ExecutiveManageObject(G, obj, zoom, quiet);
          sprintf(buf, " CmdLoad: CGO loaded into object \"%s\"\n", valid_name);
        }
      } else if(origObj) {
        sprintf(buf, " CmdLoad: CGO appended into object \"%s\"\n", valid_name);
      }
      break;

    }
    if(origObj && !quiet) {
      PRINTFB(G, FB_Executive, FB_Actions)
        "%s", buf ENDFB(G);
      OrthoRestorePrompt(G);
    }
    APIExit(G);
  }
  return APIResultOk(ok);
}

static PyObject *CmdSetStateOrder(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  char *oname;
  PyObject *order;
  CObject *obj = NULL;
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
  int result = false, state = 0;
  OrthoLineType s1;
  PyObject *coords = NULL;

  if(!PyArg_ParseTuple(args, "OsO|i", &self, &str1, &coords, &state)) {
    API_HANDLE_ERROR;
    ok_raise(2);
  }

  ok_assert(2, str1[0]);
  API_SETUP_PYMOL_GLOBALS;
  ok_assert(2, G && APIEnterBlockedNotModal(G));

  if(SelectorGetTmp(G, str1, s1) >= 0) {
    int sele1 = SelectorIndexByName(G, s1);
    if(sele1 >= 0) {
      int unblock = PAutoBlock(G);
      result = SelectorLoadCoords(G, coords, sele1, state);
      PAutoUnblock(G, unblock);
    }
    SelectorFreeTmp(G, s1);
  }

  APIExitBlocked(G);
ok_except2:
  return APIResultOk(result);
}

static PyObject *CmdLoadCoordSet(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  char *oname;
  PyObject *model;
  CObject *origObj = NULL;
  ObjectMolecule *obj;
  int frame;

  if(!PyArg_ParseTuple(args, "OsOi", &self, &oname, &model, &frame)) {
    API_HANDLE_ERROR;
    ok_raise(1);
  }

  API_SETUP_PYMOL_GLOBALS;
  ok_assert(1, G && APIEnterNotModal(G));

  origObj = ExecutiveFindObjectByName(G, oname);
  if(!origObj || origObj->type != cObjectMolecule) {
    ErrMessage(G, "LoadCoords", "named object molecule not found.");
    ok_raise(2);
  }

  PBlock(G);              /*PBlockAndUnlockAPI(); */
  obj = ObjectMoleculeLoadCoords(G, (ObjectMolecule *) origObj, model, frame);
  PUnblock(G);            /*PLockAPIAndUnblock(); */
  ok_assert(2, obj);

  if(frame < 0)
    frame = obj->NCSet - 1;

  PRINTFB(G, FB_Executive, FB_Actions)
    " CmdLoad: Coordinates appended into object \"%s\", state %d.\n",
    oname, frame + 1 ENDFB(G);
  OrthoRestorePrompt(G);

  APIExit(G);
  return APISuccess();
ok_except2:
  APIExit(G);
ok_except1:
  return APIFailure();
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
  int ok = false;
  int multiplex;
  int zoom;
  int bytes;
  int mimic;
  ok = PyArg_ParseTuple(args, "Oss#iiiiiii|zzzi", &self,
                        &oname, &fname, &bytes, &frame, &type,
                        &finish, &discrete, &quiet, &multiplex, &zoom,
                        &plugin, &object_props, &atom_props, &mimic);
  if(ok) {
    API_SETUP_PYMOL_GLOBALS;
    ok = (G != NULL);
  } else {
    API_HANDLE_ERROR;
  }
  if(ok && (ok = APIEnterNotModal(G))) {
    PRINTFD(G, FB_CCmd)
      "CmdLoad-DEBUG %s %s %d %d %d %d\n",
      oname, fname, frame, type, finish, discrete ENDFD;

    ok = ExecutiveLoad(G,
                         fname, bytes, type,
                         oname, frame, zoom,
                         discrete, finish, multiplex, quiet, plugin);

    OrthoRestorePrompt(G);
    APIExit(G);
  }
  return APIResultOk(ok);
}

static PyObject *CmdLoadTraj(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  char *fname, *oname;
  CObject *origObj = NULL;
  OrthoLineType buf;
  int frame, type;
  int interval, average, start, stop, max, image;
  OrthoLineType s1;
  char *str1;
  int ok = false;
  float shift[3];
  int quiet = 0;                /* TODO */
  char *plugin = NULL;
  ok = PyArg_ParseTuple(args, "Ossiiiiiiisifffs", &self, &oname, &fname, &frame, &type,
                        &interval, &average, &start, &stop, &max, &str1,
                        &image, &shift[0], &shift[1], &shift[2], &plugin);

  buf[0] = 0;
  if(ok) {
    API_SETUP_PYMOL_GLOBALS;
    ok = (G != NULL);
  } else {
    API_HANDLE_ERROR;
  }
  if(ok && (ok = APIEnterNotModal(G))) {
    if(str1[0])
      ok = (SelectorGetTmp(G, str1, s1) >= 0);
    else
      s1[0] = 0;                /* no selection */
    origObj = ExecutiveFindObjectByName(G, oname);
    /* check for existing object of right type, delete if not */
    if(origObj) {
      if(origObj->type != cObjectMolecule) {
        ExecutiveDelete(G, origObj->Name);
        origObj = NULL;
      }
    }
    if((type == cLoadTypeTRJ) && (plugin[0]))
      type = cLoadTypeTRJ2;
    /*printf("plugin %s %d\n",plugin,type); */
    if(origObj) {
      switch (type) {
      case cLoadTypeTRJ:       /* this is the ascii AMBER trajectory format... */
        PRINTFD(G, FB_CCmd) " CmdLoadTraj-DEBUG: loading TRJ\n" ENDFD;
        ObjectMoleculeLoadTRJFile(G, (ObjectMolecule *) origObj, fname, frame,
                                  interval, average, start, stop, max, s1, image, shift,
                                  quiet);
        /* if(finish)
           ExecutiveUpdateObjectSelection(G,origObj); unnecc */
        sprintf(buf,
                " CmdLoadTraj: \"%s\" appended into object \"%s\".\n CmdLoadTraj: %d total states in the object.\n",
                fname, oname, ((ObjectMolecule *) origObj)->NCSet);
        break;
      default:
        ok = PlugIOManagerLoadTraj(G, (ObjectMolecule *) origObj, fname, frame,
                              interval, average, start, stop, max, s1, image, shift,
                              quiet, plugin);
      }
    } else {
      PRINTFB(G, FB_CCmd, FB_Errors)
        "CmdLoadTraj-Error: must load object topology before loading trajectory.\n"
        ENDFB(G);
    }
    if(origObj) {
      PRINTFB(G, FB_Executive, FB_Actions)
        "%s", buf ENDFB(G);
      OrthoRestorePrompt(G);
    }
    SelectorFreeTmp(G, s1);
    APIExit(G);
  }
  return APIResultOk(ok);
}

static PyObject *CmdOrigin(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  char *str1, *obj;
  OrthoLineType s1;
  float v[3];
  int ok = false;
  int state;
  ok = PyArg_ParseTuple(args, "Oss(fff)i", &self, &str1, &obj, v, v + 1, v + 2, &state);
  if(ok) {
    API_SETUP_PYMOL_GLOBALS;
    ok = (G != NULL);
  } else {
    API_HANDLE_ERROR;
  }
  if(ok && (ok = APIEnterNotModal(G))) {
    if(str1[0])
      ok = (SelectorGetTmp2(G, str1, s1) >= 0);
    else
      s1[0] = 0;                /* no selection */
    ok = ExecutiveOrigin(G, s1, 1, obj, v, state);      /* TODO STATUS */
    if(str1[0])
      SelectorFreeTmp(G, s1);
    APIExit(G);
  }
  return APIResultOk(ok);
}

static PyObject *CmdSort(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  char *name;
  int ok = false;
  ok = PyArg_ParseTuple(args, "Os", &self, &name);
  if(ok) {
    API_SETUP_PYMOL_GLOBALS;
    ok = (G != NULL);
  } else {
    API_HANDLE_ERROR;
  }
  if(ok && (ok = APIEnterNotModal(G))) {
    ExecutiveSort(G, name);     /* TODO STATUS */
    APIExit(G);
  }
  return APIResultOk(ok);
}

static PyObject *CmdAssignSS(PyObject * self, PyObject * args)


/* EXPERIMENTAL */
{
  PyMOLGlobals *G = NULL;
  int ok = false;
  int state, quiet;
  char *str1, *str2;
  int preserve;
  OrthoLineType s1, s2;
  ok = PyArg_ParseTuple(args, "Osisii", &self, &str1, &state, &str2, &preserve, &quiet);
  if(ok) {
    API_SETUP_PYMOL_GLOBALS;
    ok = (G != NULL);
  } else {
    API_HANDLE_ERROR;
  }
  if(ok && (ok = APIEnterNotModal(G))) {
    ok = ((SelectorGetTmp(G, str1, s1) >= 0) && (SelectorGetTmp(G, str2, s2) >= 0));
    if(ok)
      ok = ExecutiveAssignSS(G, s1, state, s2, preserve, NULL, quiet);
    APIExit(G);
  }
  return APIResultOk(ok);
}

static PyObject *CmdSpheroid(PyObject * self, PyObject * args)


/* EXPERIMENTAL */
{
  PyMOLGlobals *G = NULL;
  char *name;
  int ok = false;
  int average;
  ok = PyArg_ParseTuple(args, "Osi", &self, &name, &average);
  if(ok) {
    API_SETUP_PYMOL_GLOBALS;
    ok = (G != NULL);
  } else {
    API_HANDLE_ERROR;
  }
  if(ok && (ok = APIEnterNotModal(G))) {
    ExecutiveSpheroid(G, name, average);        /* TODO STATUS */
    APIExit(G);
  }
  return APIResultOk(ok);
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
  OrthoLineType s1;
  int state;
  int origin;
  int ok = false;
  float animate;
  int quiet = false;            /* TODO */
  ok = PyArg_ParseTuple(args, "Osiif", &self, &str1, &state, &origin, &animate);
  if(ok) {
    API_SETUP_PYMOL_GLOBALS;
    ok = (G != NULL);
  } else {
    API_HANDLE_ERROR;
  }
  if(ok && (ok = APIEnterNotModal(G))) {
    ok = (SelectorGetTmp2(G, str1, s1) >= 0);
    if(ok)
      ok = ExecutiveCenter(G, s1, state, origin, animate, NULL, quiet);
    SelectorFreeTmp(G, s1);
    APIExit(G);
  }
  return APIResultOk(ok);
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
  OrthoLineType s1;
  float buffer;
  int state;
  int inclusive;
  int ok = false;
  float animate;
  int quiet = false;            /* TODO */
  ok =
    PyArg_ParseTuple(args, "Osfiif", &self, &str1, &buffer, &state, &inclusive, &animate);
  if(ok) {
    API_SETUP_PYMOL_GLOBALS;
    ok = (G != NULL);
  } else {
    API_HANDLE_ERROR;
  }
  if(ok && (ok = APIEnterNotModal(G))) {
    ok = (SelectorGetTmp2(G, str1, s1) >= 0);
    if(ok)
      ok = ExecutiveWindowZoom(G, s1, buffer, state, inclusive, animate, quiet);
    SelectorFreeTmp(G, s1);
    APIExit(G);
  }
  return APIResultOk(ok);
}

static PyObject *CmdIsolevel(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  float level, result = 0.0F;
  int state;
  char *name;
  int query, quiet;
  int ok = false;
  ok = PyArg_ParseTuple(args, "Osfiii", &self, &name, &level, &state, &query, &quiet);
  if(ok) {
    API_SETUP_PYMOL_GLOBALS;
    ok = (G != NULL);
  } else {
    API_HANDLE_ERROR;
  }
  if(ok && (ok = APIEnterNotModal(G))) {
    ok = ExecutiveIsolevel(G, name, level, state, query, &result, quiet);
    APIExit(G);
  }
  if(!query)
    return APIResultOk(ok);
  else
    return PyFloat_FromDouble((double) result);
}

static PyObject *CmdHAdd(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  char *str1;
  int quiet;
  int ok = false;
  int state;
  int legacy;
  ok = PyArg_ParseTuple(args, "Osiii", &self, &str1, &quiet, &state, &legacy);
  if(ok) {
    API_SETUP_PYMOL_GLOBALS;
    ok = (G != NULL);
  } else {
    API_HANDLE_ERROR;
  }
  if(ok && (ok = APIEnterNotModal(G))) {
#ifndef _PYMOL_NO_UNDO
#endif
    ExecutiveAddHydrogens(G, str1, quiet, state, legacy);
#ifndef _PYMOL_NO_UNDO
#endif
    APIExit(G);
  }
  return APIResultOk(ok);
}

static PyObject *CmdSetObjectColor(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  char *name, *color;
  int quiet;
  int ok = false;
  ok = PyArg_ParseTuple(args, "Ossi", &self, &name, &color, &quiet);
  if(ok) {
    API_SETUP_PYMOL_GLOBALS;
    ok = (G != NULL);
  } else {
    API_HANDLE_ERROR;
  }
  if(ok && (ok = APIEnterNotModal(G))) {
    ok = ExecutiveSetObjectColor(G, name, color, quiet);
    APIExit(G);
  }
  return (APIResultOk(ok));
}

static PyObject *CmdGetObjectColorIndex(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  char *str1;
  int result = -1;
  int ok = false;
  ok = PyArg_ParseTuple(args, "Os", &self, &str1);
  if(ok) {
    API_SETUP_PYMOL_GLOBALS;
    ok = (G != NULL);
  } else {
    API_HANDLE_ERROR;
  }
  if(ok && (ok = APIEnterNotModal(G))) {
    result = ExecutiveGetObjectColorIndex(G, str1);
    APIExit(G);
  }
  return (APIResultCode(result));
}

static PyObject *CmdRemove(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  char *str1;
  OrthoLineType s1;
  int quiet;
  int ok = false;
  ok = PyArg_ParseTuple(args, "Osi", &self, &str1, &quiet);
  if(ok) {
    API_SETUP_PYMOL_GLOBALS;
    ok = (G != NULL);
  } else {
    API_HANDLE_ERROR;
  }
  if(ok && (ok = APIEnterNotModal(G))) {
    ok = (SelectorGetTmp(G, str1, s1) >= 0);
    ExecutiveRemoveAtoms(G, s1, quiet); /* TODO STATUS */
    SelectorFreeTmp(G, s1);
    APIExit(G);
  }
  return APIResultOk(ok);
}

static PyObject *CmdRemovePicked(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  int i1;
  int ok = false;
  int quiet;
  ok = PyArg_ParseTuple(args, "Oii", &self, &i1, &quiet);
  if(ok) {
    API_SETUP_PYMOL_GLOBALS;
    ok = (G != NULL);
  } else {
    API_HANDLE_ERROR;
  }
  if(ok && (ok = APIEnterNotModal(G))) {
    EditorRemove(G, i1, quiet); /* TODO STATUS */
    APIExit(G);
  }
  return APIResultOk(ok);
}

static PyObject *CmdHFill(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  int ok = false;
  int quiet;
  ok = PyArg_ParseTuple(args, "Oi", &self, &quiet);
  if(ok) {
    API_SETUP_PYMOL_GLOBALS;
    ok = (G != NULL);
  } else {
    API_HANDLE_ERROR;
  }
  if(ok && (ok = APIEnterNotModal(G))) {
    EditorHFill(G, quiet);      /* TODO STATUS */
    APIExit(G);
  }
  return APIResultOk(ok);
}

static PyObject *CmdHFix(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  int ok = false;
  int quiet;
  OrthoLineType s1;
  char *str1;
  ok = PyArg_ParseTuple(args, "Osi", &self, &str1, &quiet);
  if(ok) {
    API_SETUP_PYMOL_GLOBALS;
    ok = (G != NULL);
  } else {
    API_HANDLE_ERROR;
  }
  if(ok && (ok = APIEnterNotModal(G))) {
    ok = (SelectorGetTmp(G, str1, s1) >= 0);
    EditorHFix(G, s1, quiet);   /* TODO STATUS */
    SelectorFreeTmp(G, s1);
    APIExit(G);
  }
  return APIResultOk(ok);
}

static PyObject *CmdCycleValence(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  int ok = false;
  int quiet;
  ok = PyArg_ParseTuple(args, "Oi", &self, &quiet);
  if(ok) {
    API_SETUP_PYMOL_GLOBALS;
    ok = (G != NULL);
  } else {
    API_HANDLE_ERROR;
  }
  if(ok && (ok = APIEnterNotModal(G))) {
    EditorCycleValence(G, quiet);       /* TODO STATUS */
    APIExit(G);
  }
  return APIResultOk(ok);
}

static PyObject *CmdReplace(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  int i1, i2;
  char *str1, *str2;
  int ok = false;
  int quiet;
  ok = PyArg_ParseTuple(args, "Osiisi", &self, &str1, &i1, &i2, &str2, &quiet);
  if(ok) {
    API_SETUP_PYMOL_GLOBALS;
    ok = (G != NULL);
  } else {
    API_HANDLE_ERROR;
  }
  if(ok && (ok = APIEnterNotModal(G))) {
    EditorReplace(G, str1, i1, i2, str2, quiet);        /* TODO STATUS */
    APIExit(G);
  }
  return APIResultOk(ok);
}

static PyObject *CmdSetGeometry(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  int i1, i2;
  char *str1;
  int ok = false;
  ok = PyArg_ParseTuple(args, "Osii", &self, &str1, &i1, &i2);
  if(ok) {
    API_SETUP_PYMOL_GLOBALS;
    ok = (G != NULL);
  } else {
    API_HANDLE_ERROR;
  }
  if(ok && (ok = APIEnterNotModal(G))) {
    ok = ExecutiveSetGeometry(G, str1, i1, i2);   /* TODO STATUS */
    APIExit(G);
  }
  return APIResultOk(ok);
}

static PyObject *CmdAttach(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  int i1, i2;
  char *str1;
  int ok = false;
  int quiet;
  char *name;
  ok = PyArg_ParseTuple(args, "Osiis", &self, &str1, &i1, &i2, &name, &quiet);
  if(ok) {
    API_SETUP_PYMOL_GLOBALS;
    ok = (G != NULL);
  } else {
    API_HANDLE_ERROR;
  }
  if(ok && (ok = APIEnterNotModal(G))) {
    EditorAttach(G, str1, i1, i2, name, quiet); /* TODO STATUS */
    APIExit(G);
  }
  return APIResultOk(ok);
}

static PyObject *CmdFuse(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  char *str1, *str2;
  int mode;
  int recolor;
  int ok = false;
  int move_flag;
  ok = PyArg_ParseTuple(args, "Ossiii", &self, &str1, &str2, &mode, &recolor, &move_flag);
  if(ok) {
    API_SETUP_PYMOL_GLOBALS;
    ok = (G != NULL);
  } else {
    API_HANDLE_ERROR;
  }
  if(ok && (ok = APIEnterNotModal(G))) {
    ExecutiveFuse(G, str1, str2, mode, recolor, move_flag); /* TODO STATUS */
    APIExit(G);
  }
  return APIResultOk(ok);
}

static PyObject *CmdUnpick(PyObject * self, PyObject * args)
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
    EditorInactivate(G);        /* TODO STATUS */
    APIExit(G);
  }
  return APISuccess();
}

static PyObject *CmdEdit(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  char *str0, *str1, *str2, *str3;
  OrthoLineType s0 = "";
  OrthoLineType s1 = "";
  OrthoLineType s2 = "";
  OrthoLineType s3 = "";
  int pkresi, pkbond;
  int ok = false;
  int quiet;
  ok =
    PyArg_ParseTuple(args, "Ossssiii", &self, &str0, &str1, &str2, &str3, &pkresi,
                     &pkbond, &quiet);
  if(ok) {
    API_SETUP_PYMOL_GLOBALS;
    ok = (G != NULL);
  } else {
    API_HANDLE_ERROR;
  }
  if(ok && (ok = APIEnterNotModal(G))) {
    if(!str0[0]) {
      EditorInactivate(G);
    } else {
      ok = (SelectorGetTmp(G, str0, s0) >= 0);
      if(str1[0])
        ok = (SelectorGetTmp(G, str1, s1) >= 0);
      if(str2[0])
        ok = (SelectorGetTmp(G, str2, s2) >= 0);
      if(str3[0])
        ok = (SelectorGetTmp(G, str3, s3) >= 0);
      ok = EditorSelect(G, s0, s1, s2, s3, pkresi, pkbond, quiet);
      if(s0[0])
        SelectorFreeTmp(G, s0);
      if(s1[0])
        SelectorFreeTmp(G, s1);
      if(s2[0])
        SelectorFreeTmp(G, s2);
      if(s3[0])
        SelectorFreeTmp(G, s3);
    }
    APIExit(G);
  }
  return APIResultOk(ok);
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
  OrthoLineType s1;

  int ok = false;
  ok = PyArg_ParseTuple(args, "Osii", &self, &str1, &int1, &int2);
  if(ok) {
    API_SETUP_PYMOL_GLOBALS;
    ok = (G != NULL);
  } else {
    API_HANDLE_ERROR;
  }
  if(ok && (ok = APIEnterNotModal(G))) {
    ok = (SelectorGetTmp2(G, str1, s1) >= 0);
    ExecutiveRenameObjectAtoms(G, s1, int1, int2);      /* TODO STATUS */
    SelectorFreeTmp(G, s1);
    APIExit(G);
  }
  return APIResultOk(ok);
}

static PyObject *CmdOrder(PyObject * self, PyObject * args)
{
  PyMOLGlobals *G = NULL;
  char *str1;
  int int1, int2;

  int ok = false;
  ok = PyArg_ParseTuple(args, "Osii", &self, &str1, &int1, &int2);
  if(ok) {
    API_SETUP_PYMOL_GLOBALS;
    ok = (G != NULL);
  } else {
    API_HANDLE_ERROR;
  }
  if(ok && (ok = APIEnterNotModal(G))) {
    ok = ExecutiveOrder(G, str1, int1, int2);
    APIExit(G);
  }
  return APIResultOk(ok);
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
  int ok = false;
  ok = PyArg_ParseTuple(args, "Oiiiii", &self, &int1, &x, &y, &width, &height);
  if(ok) {
    API_SETUP_PYMOL_GLOBALS;
    ok = (G != NULL);
  } else {
    API_HANDLE_ERROR;
  }
  if(G->HaveGUI && ok && (ok = APIEnterNotModal(G))) {
    switch (int1) {
#ifndef _PYMOL_NO_MAIN
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
#endif
#ifdef _MACPYMOL_XCODE
    default:
      MacPyMOL_doWindow(int1, x, y, width, height);
      break;
#endif
    }

    APIExit(G);
  }
  return APIResultOk(ok);
}

#ifdef _PYMOL_IP_EXTRAS
#include "IncentiveCopyToClipboard.h"
#endif

static PyObject *CmdCopyImage(PyObject * self, PyObject * args)
{                               /* should come in as GLUT thread just to be safe... */
  PyMOLGlobals *G = NULL;
  int ok = false;
  int quiet = true;
  ok = PyArg_ParseTuple(args, "Oi", &self, &quiet);
  if(ok) {
    API_SETUP_PYMOL_GLOBALS;
    ok = (G != NULL);
  } else {
    API_HANDLE_ERROR;
  }
  if(ok) {
    if(G->HaveGUI) {

#ifdef _PYMOL_IP_EXTRAS
      ok = IncentiveCopyToClipboard(G, quiet);
#else
#ifdef PYMOL_EVAL
      PRINTFB(G, FB_Scene, FB_Warnings)
        " Warning: Clipboard image transfers disabled in Evaluation Builds.\n" ENDFB(G);
#endif
#endif
    } else {
      ok = false;
    }
  }
  return APIResultOk(ok);
}

static PyObject *CmdGetCThreadingAPI(PyObject * self, PyObject * args)
{
  PyObject *result = PyList_New(2);
  PyList_SetItem(result, 0, PyCObject_FromVoidPtr((void *) PBlock, NULL));
  PyList_SetItem(result, 1, PyCObject_FromVoidPtr((void *) PUnblock, NULL));
  return result;
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

static PyObject *CmdVolumeColor(PyObject * self, PyObject * args)
{
  PyMOLGlobals * G = NULL;
  char * volume_name;
  PyObject * oColors;
  float * colors;
  int ncolors = 0;
  int ok = false;

  /* volume_name = (string) name
   * oColors = array of N-RGBA colors as, [ (r,g,b,a), (r,g,b,a), ..., (r,g,b,a) ],
   *   which will be unpacked to an array of N*4 floats
   */
  ok = PyArg_ParseTuple(args, "OsO", &self, &volume_name, &oColors);

  if(ok) {
    API_SETUP_PYMOL_GLOBALS;
    ok = (G != NULL);
  } else {
    API_HANDLE_ERROR;
  }

  ncolors = PyList_Size(oColors);

  PRINTFB(G, FB_ObjectVolume, FB_Blather)
    " CmdVolumeColor-Warning: ncolors=%d were passed in.\n", ncolors
  ENDFB(G);

  ok = (ncolors!=0);

  if((ok) && (ok = APIEnterNotModal(G))) {
    ok = PConvPyListToFloatVLA(oColors, &colors);
    
    if (ok) 
      ok = ExecutiveVolumeColor(G, volume_name, colors, ncolors);

    APIExit(G);
  }
  return APIResultOk(ok);
}

static PyObject *CmdVolume(PyObject *self, PyObject *args)
{ 
  PyMOLGlobals *G = NULL;
  char *volume_name, *map_name, *sele;
  float lvl, fbuf, alt_lvl;
  int mesh_mode;
  int state = -1;
  int box_mode;
  float carve;
  int ok = false;
  int map_state;
  int quiet;

  ok = PyArg_ParseTuple(args, "Ossisffiifiif", &self, &volume_name, &map_name, &box_mode,
                        &sele, &fbuf, &lvl, &mesh_mode, &state, &carve, &map_state,
                        &quiet, &alt_lvl);
  if(ok) {
    API_SETUP_PYMOL_GLOBALS;
    ok = (G != NULL);
  } else {
    API_HANDLE_ERROR;
  }
  if(ok && (ok = APIEnterNotModal(G))) {
    ok = ExecutiveVolume(G, volume_name, map_name, lvl, sele, fbuf,
                         state, carve, map_state, quiet, mesh_mode, box_mode,
                         alt_lvl);
    APIExit(G);
  }
  return APIResultOk(ok);
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
  PyObject *result = NULL;

  ok = PyArg_ParseTuple(args, "Osiii", &self, &sele, &format, &state, &quiet);
  if(ok) {
    API_SETUP_PYMOL_GLOBALS;
    ok = (G != NULL);
  } else {
    API_HANDLE_ERROR;
  }
  if(ok) {
    APIEnterBlocked(G);
    ok = (SelectorGetTmp(G, sele, s1) >= 0);
    if(ok){
      /* format : 1: mol/sybyl, 2: macromodel/mmd */
      result = PyInt_FromLong(
          ExecutiveAssignAtomTypes(G, s1, format, state, quiet));
      SelectorFreeTmp(G, s1);
    }
    APIExitBlocked(G);
  }
  return (APIAutoNone(result));
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
  ObjectMolecule **list;
  int discrete = 0;

  ok_assert(1, PyArg_ParseTuple(args, "Os", &self, &str1));
  API_SETUP_PYMOL_GLOBALS;
  ok_assert(1, G && APIEnterBlockedNotModal(G));
  ok_assert(2, SelectorGetTmp(G, str1, s1) >= 0);

  if((list = ExecutiveGetObjectMoleculeVLA(G, s1))) {
    unsigned int i, size = VLAGetSize(list);
    for(i = 0; i < size; i++)
      if(list[i]->DiscreteFlag)
        discrete++;
    VLAFreeP(list);
  }

  SelectorFreeTmp(G, s1);
ok_except2:
  APIExitBlocked(G);
  return Py_BuildValue("i", discrete);
ok_except1:
  API_HANDLE_ERROR;
  return APIAutoNone(NULL);
}

/*
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
  ok_assert(1, G && APIEnterBlockedNotModal(G));

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
    const cif_array *arr = obj->m_cifdata->get_arr(key);
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

static PyMethodDef Cmd_methods[] = {
  {"_get_c_threading_api", CmdGetCThreadingAPI, METH_VARARGS},
  {"_del", Cmd_Del, METH_VARARGS},
  {"_get_global_C_object", Cmd_GetGlobalCObject, METH_VARARGS},
  {"glViewport", Cmd_glViewport, METH_VARARGS},
  {"_new", Cmd_New, METH_VARARGS},
  {"_start", Cmd_Start, METH_VARARGS},
  {"_stop", Cmd_Stop, METH_VARARGS},
  {"_idle", Cmd_Idle, METH_VARARGS},
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
  {"copy_image", CmdCopyImage, METH_VARARGS},
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
  {"export_dots", CmdExportDots, METH_VARARGS},
  {"export_coords", CmdExportCoords, METH_VARARGS},
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
  {"get_matrix", CmdGetMatrix, METH_VARARGS},
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
  {"get_setting_tuple", CmdGetSettingTuple, METH_VARARGS},
  {"get_setting_text", CmdGetSettingText, METH_VARARGS},
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
  {"import_coords", CmdImportCoords, METH_VARARGS},
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
  {"p_glut_event", CmdPGlutEvent, METH_VARARGS},
  {"p_glut_get_redisplay", CmdPGlutGetRedisplay, METH_VARARGS},
  {"paste", CmdPaste, METH_VARARGS},
  {"png", CmdPNG, METH_VARARGS},
  {"pop", CmdPop, METH_VARARGS},
  {"protect", CmdProtect, METH_VARARGS},
  {"pseudoatom", CmdPseudoatom, METH_VARARGS},
  {"push_undo", CmdPushUndo, METH_VARARGS},
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
  {"runwxpymol", CmdRunWXPyMOL, METH_VARARGS},
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
  {"set_colorection_name", CmdSetColorectionName, METH_VARARGS},
  {"set_dihe", CmdSetDihe, METH_VARARGS},
  {"set_discrete", CmdSetDiscrete, METH_VARARGS},
  {"set_feedback", CmdSetFeedbackMask, METH_VARARGS},
  {"set_frame", CmdSetFrame, METH_VARARGS},
  {"set_name", CmdSetName, METH_VARARGS},
  {"set_geometry", CmdSetGeometry, METH_VARARGS},
  {"set_matrix", CmdSetMatrix, METH_VARARGS},
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
  {"volume_color", CmdVolumeColor, METH_VARARGS},
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

void init_cmd(void)
{
#if PY_MAJOR_VERSION < 3
  PyUnicode_SetDefaultEncoding("utf-8");
  Py_InitModule4("pymol._cmd",
                 Cmd_methods,
                 "PyMOL _cmd internal API -- PRIVATE: DO NOT USE!",
                 PyCObject_FromVoidPtr((void *) &SingletonPyMOLGlobals, NULL),
                 PYTHON_API_VERSION);
#endif
}

#if PY_MAJOR_VERSION >= 3
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
#endif

#ifdef __cplusplus
}
#endif

#else
typedef int this_file_is_no_longer_empty;
#endif
