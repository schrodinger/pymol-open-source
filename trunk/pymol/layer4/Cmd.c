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

#include<Python.h>

#include"os_std.h"

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
#include"Executive.h"
#include"Selector.h"
#include"main.h"
#include"Scene.h"
#include"Setting.h"
#include"Movie.h"
#include"Export.h"
#include"P.h"
#include"PConv.h"
#include"Control.h"
#include"Editor.h"
#include"Wizard.h"

#define cLoadTypePDB 0
#define cLoadTypeMOL 1
#define cLoadTypeSDF 2
#define cLoadTypeMOLStr 3
#define cLoadTypeMMD 4
#define cLoadTypeMMDSeparate 5
#define cLoadTypeMMDStr 6
#define cLoadTypeXPLORMap 7
#define cLoadTypeChemPyModel 8
#define cLoadTypePDBStr 9
#define cLoadTypeChemPyBrick 10
#define cLoadTypeChemPyMap 11
#define cLoadTypeCallback 12
#define cLoadTypeCGO 13
#define cLoadTypeR3D 14
#define cLoadTypeXYZ 15
#define cLoadTypeCCP4Map 18
#define cLoadTypePMO  19

#define tmpSele "_tmp"
#define tmpSele1 "_tmp1"
#define tmpSele2 "_tmp2"

static int flush_count = 0;

int run_only_once = true;

/* NOTE: the glut_thread_keep_out variable can only be changed by the thread
   holding the API lock, therefore this is safe even through increment
   isn't (necessarily) atomic. */

int PyThread_get_thread_ident();

static void APIEntry(void) /* assumes API is locked */
{
  PRINTFD(FB_API)
    " APIEntry-DEBUG: as thread 0x%x.\n",PyThread_get_thread_ident()
    ENDFD;

if(PyMOLTerminating) {/* try to bail */
#ifdef WIN32
	abort();
#endif
    exit(0);
	}
 P_glut_thread_keep_out++;  
  PUnblock();
}

static PyObject *APISuccess(void)
{
  return(Py_BuildValue("i",1));
}

static PyObject *APIFailure(void)
{
  Py_INCREF(Py_None);
  return(Py_None);
}

static PyObject *APIStatus(int status) /* status/integer return */
{
  return(Py_BuildValue("i",status));
}

static PyObject *APIIncRef(PyObject *result) /* automatically own Py_None */
{
  Py_INCREF(result);
  return(result);
}
static PyObject *APIAutoNone(PyObject *result) /* automatically own Py_None */
{
  if(result==Py_None)
    Py_INCREF(result);
  else if(result==NULL) {
    result=Py_None;
    Py_INCREF(result);
  } 
  return(result);
}


static void APIExit(void) /* assumes API is locked */
{
  PBlock();
  P_glut_thread_keep_out--;
  PRINTFD(FB_API)
    " APIExit-DEBUG: as thread 0x%x.\n",PyThread_get_thread_ident()
    ENDFD;
}

static PyObject *CmdAlign(PyObject *self,   PyObject *args);
static PyObject *CmdAlter(PyObject *self,   PyObject *args);
static PyObject *CmdAlterState(PyObject *self,   PyObject *args);
static PyObject *CmdAttach(PyObject *self, 	PyObject *args);
static PyObject *CmdBackgroundColor(PyObject *dummy, PyObject *args);
static PyObject *CmdBond(PyObject *dummy, PyObject *args);
static PyObject *CmdButton(PyObject *dummy, PyObject *args);
static PyObject *CmdCartoon(PyObject *self, 	PyObject *args);
static PyObject *CmdClip(PyObject *self, 	PyObject *args);
static PyObject *CmdCls(PyObject *self, 	PyObject *args);
static PyObject *CmdColor(PyObject *self, PyObject *args);
static PyObject *CmdColorDef(PyObject *self, 	PyObject *args);
static PyObject *CmdCombineObjectTTT(PyObject *self, 	PyObject *args);
static PyObject *CmdCopy(PyObject *self, PyObject *args);
static PyObject *CmdCountStates(PyObject *self, PyObject *args);
static PyObject *CmdCreate(PyObject *self, PyObject *args);
static PyObject *CmdCycleValence(PyObject *self, PyObject *args);
static PyObject *CmdDelete(PyObject *self, PyObject *args);
static PyObject *CmdDirty(PyObject *self, 	PyObject *args);
static PyObject *CmdDist(PyObject *dummy, PyObject *args);
static PyObject *CmdDistance(PyObject *dummy, PyObject *args);
static PyObject *CmdDo(PyObject *self, 	PyObject *args);
static PyObject *CmdDump(PyObject *self, 	PyObject *args);
static PyObject *CmdEdit(PyObject *self, 	PyObject *args);
static PyObject *CmdTorsion(PyObject *self, PyObject *args);
static PyObject *CmdExportDots(PyObject *self, PyObject *args);
static PyObject *CmdFeedback(PyObject *dummy, PyObject *args);
static PyObject *CmdFindPairs(PyObject *dummy, PyObject *args);
static PyObject *CmdFit(PyObject *dummy, PyObject *args);
static PyObject *CmdFitPairs(PyObject *dummy, PyObject *args);
static PyObject *CmdFlag(PyObject *self, 	PyObject *args);
static PyObject *CmdFlushNow(PyObject *self, 	PyObject *args);
static PyObject *CmdFocus(PyObject *self, 	PyObject *args);
static PyObject *CmdFullScreen(PyObject *self,PyObject *args);
static PyObject *CmdFuse(PyObject *self, 	PyObject *args);
static PyObject *CmdHAdd(PyObject *self, PyObject *args);
static PyObject *CmdIdentify(PyObject *dummy, PyObject *args);
static PyObject *CmdIndex(PyObject *dummy, PyObject *args);
static PyObject *CmdIntraFit(PyObject *dummy, PyObject *args);
static PyObject *CmdInvert(PyObject *self, PyObject *args);
static PyObject *CmdIsomesh(PyObject *self, 	PyObject *args);
static PyObject *CmdFinishObject(PyObject *self, PyObject *args);
static PyObject *CmdFrame(PyObject *self, PyObject *args);
static PyObject *CmdGet(PyObject *self, 	PyObject *args);
static PyObject *CmdGetArea(PyObject *self, 	PyObject *args);
static PyObject *CmdGetColor(PyObject *self, 	PyObject *args);
static PyObject *CmdGetDihe(PyObject *self, 	PyObject *args);
static PyObject *CmdGetFeedback(PyObject *dummy, PyObject *args);
static PyObject *CmdGetFrame(PyObject *self, 	PyObject *args);
static PyObject *CmdSetGeometry(PyObject *self, 	PyObject *args);
static PyObject *CmdGetPDB(PyObject *dummy, PyObject *args);
static PyObject *CmdGetMatrix(PyObject *self, 	PyObject *args);
static PyObject *CmdGetMinMax(PyObject *self, 	PyObject *args);
static PyObject *CmdGetModel(PyObject *dummy, PyObject *args);
static PyObject *CmdGetMoment(PyObject *self, 	PyObject *args);
static PyObject *CmdGetNames(PyObject *self, 	PyObject *args);
static PyObject *CmdGetPhiPsi(PyObject *self, 	PyObject *args);
static PyObject *CmdGetPosition(PyObject *self, 	PyObject *args);
static PyObject *CmdGetPovRay(PyObject *dummy, PyObject *args);
static PyObject *CmdGetRenderer(PyObject *self,  PyObject *args);
static PyObject *CmdGetSetting(PyObject *self, 	PyObject *args);
static PyObject *CmdGetSettingText(PyObject *self, 	PyObject *args);
static PyObject *CmdGetSettingTuple(PyObject *self, 	PyObject *args);
static PyObject *CmdGetSettingUpdates(PyObject *self, 	PyObject *args);
static PyObject *CmdGetState(PyObject *self, 	PyObject *args);
static PyObject *CmdGetType(PyObject *self, 	PyObject *args);
static PyObject *CmdGetWizard(PyObject *self, PyObject *args);
static PyObject *CmdGetView(PyObject *self, 	PyObject *args);
static PyObject *CmdMask(PyObject *self, PyObject *args);
static PyObject *CmdMem(PyObject *self, 	PyObject *args);
static PyObject *CmdLabel(PyObject *self,   PyObject *args);
static PyObject *CmdLoad(PyObject *self, 	PyObject *args);
static PyObject *CmdLoadCoords(PyObject *self, PyObject *args);
static PyObject *CmdLoadObject(PyObject *self, PyObject *args);
static PyObject *CmdLoadPNG(PyObject *self, PyObject *args);
static PyObject *CmdMapSetBorder(PyObject *self, 	PyObject *args);
static PyObject *CmdMClear(PyObject *self, 	PyObject *args);
static PyObject *CmdMDo(PyObject *self, 	PyObject *args);
static PyObject *CmdMMatrix(PyObject *self, 	PyObject *args);
static PyObject *CmdMove(PyObject *self, 	PyObject *args);
static PyObject *CmdMPlay(PyObject *self, 	PyObject *args);
static PyObject *CmdMPNG(PyObject *self, 	PyObject *args);
static PyObject *CmdMSet(PyObject *self, 	PyObject *args);
static PyObject *CmdMultiSave(PyObject *self, 	PyObject *args);
static PyObject *CmdExportCoords(PyObject *self, 	PyObject *args);
static PyObject *CmdImportCoords(PyObject *self, 	PyObject *args);
static PyObject *CmdOrigin(PyObject *self, PyObject *args);
static PyObject *CmdOnOff(PyObject *self, 	PyObject *args);
static PyObject *CmdOrient(PyObject *dummy, PyObject *args);
static PyObject *CmdOverlap(PyObject *self, 	PyObject *args);
static PyObject *CmdPaste(PyObject *self, 	PyObject *args);
static PyObject *CmdPNG(PyObject *self, 	PyObject *args);
static PyObject *CmdProtect(PyObject *self, PyObject *args);
static PyObject *CmdQuit(PyObject *self, 	PyObject *args);
static PyObject *CmdRay(PyObject *self, 	PyObject *args);
static PyObject *CmdRebuild(PyObject *self, PyObject *args);
static PyObject *CmdRecolor(PyObject *self, PyObject *args);
static PyObject *CmdHFill(PyObject *self, PyObject *args);
static PyObject *CmdRemove(PyObject *self, PyObject *args);
static PyObject *CmdRemovePicked(PyObject *self, PyObject *args);
static PyObject *CmdRename(PyObject *self, 	PyObject *args);
static PyObject *CmdReplace(PyObject *self, PyObject *args);
static PyObject *CmdReset(PyObject *self, PyObject *args);
static PyObject *CmdResetRate(PyObject *dummy, PyObject *args);
static PyObject *CmdRefresh(PyObject *self, 	PyObject *args);
static PyObject *CmdRefreshNow(PyObject *self, 	PyObject *args);
static PyObject *CmdReady(PyObject *dummy, PyObject *args);
static PyObject *CmdRock(PyObject *self, PyObject *args);
static PyObject *CmdRunPyMOL(PyObject *dummy, PyObject *args);
static PyObject *CmdSelect(PyObject *self, PyObject *args);
static PyObject *CmdSetMatrix(PyObject *self, 	PyObject *args);
static PyObject *CmdSet(PyObject *self, 	PyObject *args);
static PyObject *CmdLegacySet(PyObject *self, 	PyObject *args);
static PyObject *CmdSetDihe(PyObject *self, 	PyObject *args);
static PyObject *CmdSetFeedbackMask(PyObject *dummy, PyObject *args);
static PyObject *CmdSetFrame(PyObject *self, PyObject *args);
static PyObject *CmdSetTitle(PyObject *self, PyObject *args);
static PyObject *CmdSetView(PyObject *self, 	PyObject *args);
static PyObject *CmdSetWizard(PyObject *self, PyObject *args);
static PyObject *CmdRefreshWizard(PyObject *dummy, PyObject *args);
static PyObject *CmdShowHide(PyObject *self, 	PyObject *args);
static PyObject *CmdSort(PyObject *dummy, PyObject *args);
static PyObject *CmdSplash(PyObject *dummy, PyObject *args);
static PyObject *CmdSpheroid(PyObject *dummy, PyObject *args);
static PyObject *CmdStereo(PyObject *self, PyObject *args);
static PyObject *CmdSystem(PyObject *dummy, PyObject *args);
static PyObject *CmdSymExp(PyObject *dummy, PyObject *args);
static PyObject *CmdTest(PyObject *self, 	PyObject *args);
static PyObject *CmdTransformObject(PyObject *self, 	PyObject *args);
static PyObject *CmdTranslateAtom(PyObject *self, PyObject *args);
static PyObject *CmdTurn(PyObject *self, 	PyObject *args);
static PyObject *CmdViewport(PyObject *self, 	PyObject *args);
static PyObject *CmdZoom(PyObject *self, PyObject *args);
static PyObject *CmdUnpick(PyObject *dummy, PyObject *args);
static PyObject *CmdUpdate(PyObject *dummy, PyObject *args);
static PyObject *CmdWaitQueue(PyObject *self, 	PyObject *args);
static PyObject *CmdUndo(PyObject *self, 	PyObject *args);
static PyObject *CmdPushUndo(PyObject *self, 	PyObject *args);


static PyMethodDef Cmd_methods[] = {
	{"align",	              CmdAlign,                METH_VARARGS },
	{"alter",	              CmdAlter,                METH_VARARGS },
	{"alter_state",           CmdAlterState,           METH_VARARGS },
	{"attach",                CmdAttach,               METH_VARARGS },
   {"bg_color",              CmdBackgroundColor,      METH_VARARGS },
	{"bond",                  CmdBond,                 METH_VARARGS },
   {"button",                CmdButton,               METH_VARARGS },
   {"cartoon",               CmdCartoon,              METH_VARARGS },
	{"clip",	                 CmdClip,                 METH_VARARGS },
	{"cls",	                 CmdCls,                  METH_VARARGS },
	{"color",	              CmdColor,                METH_VARARGS },
	{"colordef",	           CmdColorDef,             METH_VARARGS },
   {"combine_object_ttt",    CmdCombineObjectTTT,     METH_VARARGS },
	{"copy",                  CmdCopy,                 METH_VARARGS },
	{"create",                CmdCreate,               METH_VARARGS },
	{"count_states",          CmdCountStates,          METH_VARARGS },
	{"cycle_valence",         CmdCycleValence,         METH_VARARGS },
	{"delete",                CmdDelete,               METH_VARARGS },
	{"dirty",                 CmdDirty,                METH_VARARGS },
	{"distance",	           CmdDistance,             METH_VARARGS },
	{"dist",    	           CmdDist,                 METH_VARARGS },
	{"do",	                 CmdDo,                   METH_VARARGS },
	{"dump",	                 CmdDump,                 METH_VARARGS },
   {"edit",                  CmdEdit,                 METH_VARARGS },
   {"torsion",               CmdTorsion,              METH_VARARGS },
	{"export_dots",           CmdExportDots,           METH_VARARGS },
	{"export_coords",         CmdExportCoords,         METH_VARARGS },
	{"feedback",              CmdFeedback,             METH_VARARGS },
	{"find_pairs",            CmdFindPairs,            METH_VARARGS },
	{"finish_object",         CmdFinishObject,         METH_VARARGS },
	{"fit",                   CmdFit,                  METH_VARARGS },
	{"fit_pairs",             CmdFitPairs,             METH_VARARGS },
	{"flag",                  CmdFlag,                 METH_VARARGS },
	{"frame",	              CmdFrame,                METH_VARARGS },
   {"flush_now",             CmdFlushNow,             METH_VARARGS },
   {"focus",                 CmdFocus,                METH_VARARGS },
   {"full_screen",           CmdFullScreen,           METH_VARARGS },
   {"fuse",                  CmdFuse,                 METH_VARARGS },
	{"get",	                 CmdGet,                  METH_VARARGS },
	{"get_area",              CmdGetArea,              METH_VARARGS },
	{"get_color",             CmdGetColor,             METH_VARARGS },
	{"get_dihe",              CmdGetDihe,              METH_VARARGS },
	{"get_frame",             CmdGetFrame,             METH_VARARGS },
	{"get_feedback",          CmdGetFeedback,          METH_VARARGS },
	{"get_matrix",	           CmdGetMatrix,            METH_VARARGS },
	{"get_min_max",           CmdGetMinMax,            METH_VARARGS },
	{"get_model",	           CmdGetModel,             METH_VARARGS },
	{"get_moment",	           CmdGetMoment,            METH_VARARGS },
   {"get_names",             CmdGetNames,             METH_VARARGS },
	{"get_position",	        CmdGetPosition,          METH_VARARGS },
	{"get_povray",	           CmdGetPovRay,            METH_VARARGS },
	{"get_pdb",	              CmdGetPDB,               METH_VARARGS },
   {"get_phipsi",            CmdGetPhiPsi,            METH_VARARGS },
   {"get_renderer",          CmdGetRenderer,          METH_VARARGS },
	{"get_setting",           CmdGetSetting,           METH_VARARGS },
	{"get_setting_tuple",     CmdGetSettingTuple,      METH_VARARGS },
	{"get_setting_text",      CmdGetSettingText,       METH_VARARGS },
   {"get_setting_updates",   CmdGetSettingUpdates,    METH_VARARGS },
	{"get_state",             CmdGetState,             METH_VARARGS },
	{"get_type",              CmdGetType,              METH_VARARGS },
   {"get_view",              CmdGetView,              METH_VARARGS },
   {"get_wizard",            CmdGetWizard,            METH_VARARGS },
	{"h_add",                 CmdHAdd,                 METH_VARARGS },
	{"h_fill",                CmdHFill,                METH_VARARGS },
   {"identify",              CmdIdentify,             METH_VARARGS },
	{"import_coords",         CmdImportCoords,         METH_VARARGS },
   {"index",                 CmdIndex,                METH_VARARGS },
	{"intrafit",              CmdIntraFit,             METH_VARARGS },
   {"invert",                CmdInvert,               METH_VARARGS },
	{"isomesh",	              CmdIsomesh,              METH_VARARGS },
   {"wait_queue",            CmdWaitQueue,            METH_VARARGS },
   {"label",                 CmdLabel,                METH_VARARGS },
	{"load",	                 CmdLoad,                 METH_VARARGS },
	{"load_coords",           CmdLoadCoords,           METH_VARARGS },
	{"load_png",              CmdLoadPNG,              METH_VARARGS },
	{"load_object",           CmdLoadObject,           METH_VARARGS },
   {"map_set_border",        CmdMapSetBorder,         METH_VARARGS },
	{"mask",	                 CmdMask,                 METH_VARARGS },
	{"mclear",	              CmdMClear,               METH_VARARGS },
	{"mdo",	                 CmdMDo,                  METH_VARARGS },
	{"mem",	                 CmdMem,                  METH_VARARGS },
	{"move",	                 CmdMove,                 METH_VARARGS },
	{"mset",	                 CmdMSet,                 METH_VARARGS },
	{"mplay",	              CmdMPlay,                METH_VARARGS },
	{"mpng_",	              CmdMPNG,                 METH_VARARGS },
	{"mmatrix",	              CmdMMatrix,              METH_VARARGS },
	{"multisave",             CmdMultiSave,            METH_VARARGS },
	{"origin",	              CmdOrigin,               METH_VARARGS },
	{"orient",	              CmdOrient,               METH_VARARGS },
	{"onoff",                 CmdOnOff,                METH_VARARGS },
	{"overlap",               CmdOverlap,              METH_VARARGS },
	{"paste",	              CmdPaste,                METH_VARARGS },
	{"png",	                 CmdPNG,                  METH_VARARGS },
	{"protect",	              CmdProtect,              METH_VARARGS },
	{"push_undo",	           CmdPushUndo,             METH_VARARGS },
	{"quit",	                 CmdQuit,                 METH_VARARGS },
	{"ready",                 CmdReady,                METH_VARARGS },
   {"rebuild",               CmdRebuild,              METH_VARARGS },
   {"recolor",               CmdRecolor,              METH_VARARGS },
	{"refresh",               CmdRefresh,              METH_VARARGS },
	{"refresh_now",           CmdRefreshNow,           METH_VARARGS },
	{"refresh_wizard",        CmdRefreshWizard,        METH_VARARGS },
	{"remove",	              CmdRemove,               METH_VARARGS },
	{"remove_picked",         CmdRemovePicked,         METH_VARARGS },
	{"render",	              CmdRay,                  METH_VARARGS },
   {"rename",                CmdRename,               METH_VARARGS },
   {"replace",               CmdReplace,              METH_VARARGS },
	{"reset",                 CmdReset,                METH_VARARGS },
	{"reset_rate",	           CmdResetRate,            METH_VARARGS },
	{"rock",	                 CmdRock,                 METH_VARARGS },
	{"runpymol",	           CmdRunPyMOL,             METH_VARARGS },
	{"select",                CmdSelect,               METH_VARARGS },
	{"set",	                 CmdSet,                  METH_VARARGS },
	{"legacy_set",            CmdLegacySet,            METH_VARARGS },
	{"set_dihe",              CmdSetDihe,              METH_VARARGS },
	{"set_feedback",          CmdSetFeedbackMask,      METH_VARARGS },
   {"set_geometry",          CmdSetGeometry,          METH_VARARGS },
	{"set_title",             CmdSetTitle,             METH_VARARGS },
	{"set_wizard",            CmdSetWizard,            METH_VARARGS },
   {"set_view",              CmdSetView,              METH_VARARGS },
	{"setframe",	           CmdSetFrame,             METH_VARARGS },
	{"showhide",              CmdShowHide,             METH_VARARGS },
	{"set_matrix",	           CmdSetMatrix,            METH_VARARGS },
	{"sort",                  CmdSort,                 METH_VARARGS },
   {"spheroid",              CmdSpheroid,             METH_VARARGS },
	{"splash",                CmdSplash,               METH_VARARGS },
	{"stereo",	              CmdStereo,               METH_VARARGS },
	{"system",	              CmdSystem,               METH_VARARGS },
	{"symexp",	              CmdSymExp,               METH_VARARGS },
	{"test",	                 CmdTest,                 METH_VARARGS },
	{"transform_object",      CmdTransformObject,      METH_VARARGS },
	{"translate_atom",        CmdTranslateAtom,        METH_VARARGS },
	{"turn",	                 CmdTurn,                 METH_VARARGS },
	{"viewport",              CmdViewport,             METH_VARARGS },
	{"undo",                  CmdUndo,                 METH_VARARGS },
	{"unpick",                CmdUnpick,               METH_VARARGS },
	{"update",                CmdUpdate,               METH_VARARGS },
	{"zoom",	                 CmdZoom,                 METH_VARARGS },
	{NULL,		              NULL}     /* sentinel */        
};

static PyObject *CmdCombineObjectTTT(PyObject *self, 	PyObject *args)
{
  char *name;
  PyObject *m;
  float ttt[16];
  int ok = false;
  ok = PyArg_ParseTuple(args,"sO",&name,&m);
  if(ok) {
    if(PConvPyListToFloatArrayInPlace(m,ttt,16)) {
      APIEntry();
      ok = ExecutiveCombineObjectTTT(name,ttt);
      APIExit();
    } else {
      PRINTFB(FB_CCmd,FB_Errors)
        "CmdCombineObjectTTT-Error: bad matrix\n"
        ENDFB;
      ok=false;
    }
  }
  return(APIStatus(ok));
}

static PyObject *CmdGetColor(PyObject *self, PyObject *args)
{
  char *name;
  int mode;
  int ok = false;
  int a,nc,nvc;
  float *rgb;
  int index;
  PyObject *result = NULL;
  PyObject *tup;
  ok = PyArg_ParseTuple(args,"si",&name,&mode);
  if(ok) {
    APIEntry();
    switch(mode) {
    case 0: /* by name or index, return floats */
      index = ColorGetIndex(name);
      if(index>=0) {
        rgb = ColorGet(index);
        tup = PyTuple_New(3);
        PyTuple_SetItem(tup,0,PyFloat_FromDouble(*(rgb++)));
        PyTuple_SetItem(tup,1,PyFloat_FromDouble(*(rgb++)));
        PyTuple_SetItem(tup,2,PyFloat_FromDouble(*rgb));
        result=tup;
      }
      break;
    case 1: /* get color names with NO NUMBERS in their names */
      PBlock();
      nc=ColorGetNColor();
      nvc=0;
      for(a=0;a<nc;a++) {
        if(ColorGetStatus(a))
          nvc++;
      }
      result = PyList_New(nvc);
      nvc=0;
      for(a=0;a<nc;a++) {
        if(ColorGetStatus(a)) {
          tup = PyTuple_New(2);
          PyTuple_SetItem(tup,0,PyString_FromString(ColorGetName(a)));
          PyTuple_SetItem(tup,1,PyInt_FromLong(a));
          PyList_SetItem(result,nvc++,tup);
        }
      }
      PUnblock();
      break;
    }
    APIExit();
  }
  return(APIAutoNone(result));
}

static PyObject *CmdMultiSave(PyObject *self, PyObject *args)
{
  char *name,*object;
  int append,state;
  int ok = false;
  ok = PyArg_ParseTuple(args,"ssii",&name,&object,&state,&append);
  if(ok) {
    APIEntry();
    ok = ExecutiveMultiSave(name,object,state,append);
    APIExit();
  }
  return(APIStatus(ok));
}

static PyObject *CmdMapSetBorder(PyObject *self, PyObject *args)
{
  char *name;
  float level;
  int ok = false;
  ok = PyArg_ParseTuple(args,"sf",&name,&level);
  if(ok) {
    APIEntry();
    ok = ExecutiveMapSetBorder(name,level);
    APIExit();
  }
  return(APIStatus(ok));
}

static PyObject *CmdGetRenderer(PyObject *self, PyObject *args)
{
  char *vendor,*renderer,*version;
  APIEntry();
  SceneGetCardInfo(&vendor,&renderer,&version);
  APIExit();
  return Py_BuildValue("(sss)",vendor,renderer,version);
}

static PyObject *CmdTranslateAtom(PyObject *self, PyObject *args)
{
  char *str1;
  int state,log,mode;
  float v[3];
  OrthoLineType s1;
  int ok = false;
  ok = PyArg_ParseTuple(args,"sfffiii",&str1,v,v+1,v+2,&state,&mode,&log);
  if(ok) {
    SelectorGetTmp(str1,s1);
    APIEntry();
    ok = ExecutiveTranslateAtom(s1,v,state,mode,log);
    APIExit();
    SelectorFreeTmp(s1);
  }
  return(APIStatus(ok));
}

static PyObject *CmdTransformObject(PyObject *self, PyObject *args)
{
  char *name,*sele;
  int state,log;
  PyObject *m;
  float ttt[16];
  int ok = false;
  ok = PyArg_ParseTuple(args,"siOis",&name,&state,&m,&log,&sele);
  if(ok) {
    if(PConvPyListToFloatArrayInPlace(m,ttt,16)) {
      APIEntry();
      ok = ExecutiveTransformObjectSelection(name,state,sele,log,ttt);
      APIExit();
    } else {
      PRINTFB(FB_CCmd,FB_Errors)
        "CmdTransformObject-DEBUG: bad matrix\n"
        ENDFB;
      ok=false;
    }
  }
  return(APIStatus(ok));
}

static PyObject *CmdLoadPNG(PyObject *self, PyObject *args)
{
  char *str1;
  int ok = false;
  ok = PyArg_ParseTuple(args,"s",&str1);
  if(ok) {
    APIEntry();
    ok = SceneLoadPNG(str1,true);
    APIExit();
  }
  return(APIStatus(ok));
}

static PyObject *CmdBackgroundColor(PyObject *self, PyObject *args)
{
  char *str1;
  int ok = false;
  ok = PyArg_ParseTuple(args,"s",&str1);
  if(ok) {
    APIEntry();
    ok = SettingSetfv(cSetting_bg_rgb,ColorGetNamed(str1));
    APIExit();
  }
  return(APIStatus(ok));
}

static PyObject *CmdGetPosition(PyObject *self, 	PyObject *args)
{
  PyObject *result;
  float v[3];
  APIEntry();
  SceneGetPos(v);
  APIExit();
  result=PConvFloatArrayToPyList(v,3);
  return(result);
}

static PyObject *CmdGetPhiPsi(PyObject *self, 	PyObject *args)
{
  char *str1;
  OrthoLineType s1;
  int state;
  PyObject *result = Py_None;
  PyObject *key = Py_None;
  PyObject *value = Py_None;
  int *iVLA=NULL;
  float *pVLA,*sVLA;
  int l;
  int *i;
  ObjectMolecule **o,**oVLA=NULL;
  int a;
  float *s,*p;
  int ok =  PyArg_ParseTuple(args,"si",&str1,&state);
  if(ok) {
    APIEntry();
    SelectorGetTmp(str1,s1);
    l = ExecutivePhiPsi(s1,&oVLA,&iVLA,&pVLA,&sVLA,state);
    SelectorFreeTmp(s1);
    APIExit();
    if(iVLA) {
      result=PyDict_New();
      i = iVLA;
      o = oVLA;
      p = pVLA;
      s = sVLA;
      for(a=0;a<l;a++) {
        key = PyTuple_New(2);      
        PyTuple_SetItem(key,1,PyInt_FromLong(*(i++)+1)); /* +1 for index */
        PyTuple_SetItem(key,0,PyString_FromString((*(o++))->Obj.Name));
        value = PyTuple_New(2);      
        PyTuple_SetItem(value,0,PyFloat_FromDouble(*(p++))); /* +1 for index */
        PyTuple_SetItem(value,1,PyFloat_FromDouble(*(s++)));
        PyDict_SetItem(result,key,value);
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
  return(APIAutoNone(result));
}

static PyObject *CmdAlign(PyObject *self, 	PyObject *args) {
  char *str2,*str3;
  OrthoLineType s2="",s3="";
  float result = -1.0;
  int ok = false;
  ok = PyArg_ParseTuple(args,"ss",&str2,&str3);
  if(ok) {
    PRINTFD(FB_CCmd)
      "CmdAlign-DEBUG %s %s\n",
      str2,str3
      ENDFD;
    
    APIEntry();
    SelectorGetTmp(str2,s2);
    SelectorGetTmp(str3,s3);
    result = ExecutiveAlign(s2,s3);
    SelectorFreeTmp(s2);
    SelectorFreeTmp(s3);
    APIExit();
  }
  return Py_BuildValue("f",result);
}

static PyObject *CmdGetSettingUpdates(PyObject *self, 	PyObject *args)
{
  PyObject *result = NULL;
  APIEntry();
  result = SettingGetUpdateList(NULL);
  APIExit();
  return(APIAutoNone(result));
}

static PyObject *CmdGetView(PyObject *self, 	PyObject *args)
{
  SceneViewType view;
  APIEntry();
  SceneGetView(view);
  APIExit();
  return(Py_BuildValue("(fffffffffffffffffffffffff)",
                   view[ 0],view[ 1],view[ 2],view[ 3], /* 4x4 mat */
                   view[ 4],view[ 5],view[ 6],view[ 7],
                   view[ 8],view[ 9],view[10],view[11],
                   view[12],view[13],view[14],view[15],
                   view[16],view[17],view[18], /* pos */
                   view[19],view[20],view[21], /* origin */
                   view[22],view[23], /* clip */
                   view[24] /* orthoscopic*/
                       ));
}
static PyObject *CmdSetView(PyObject *self, 	PyObject *args)
{
  SceneViewType view;
  int ok=PyArg_ParseTuple(args,"(fffffffffffffffffffffffff)",
                   &view[ 0],&view[ 1],&view[ 2],&view[ 3], /* 4x4 mat */
                   &view[ 4],&view[ 5],&view[ 6],&view[ 7],
                   &view[ 8],&view[ 9],&view[10],&view[11],
                   &view[12],&view[13],&view[14],&view[15],
                   &view[16],&view[17],&view[18], /* pos */
                   &view[19],&view[20],&view[21], /* origin */
                   &view[22],&view[23], /* clip */
                   &view[24] /* orthoscopic*/
                   );
  if(ok) {
    APIEntry();
    SceneSetView(view); /* TODO STATUS */
    APIExit();
  }
  return(APIStatus(ok));
}
static PyObject *CmdGetState(PyObject *self, 	PyObject *args)
{
  return(APIStatus(SceneGetState()));
}

static PyObject *CmdGetFrame(PyObject *self, 	PyObject *args)
{
  return(APIStatus(SceneGetFrame()));
}

static PyObject *CmdSetTitle(PyObject *self, PyObject *args)
{
  char *str1,*str2;
  int int1;
  int ok = false;
  ok = PyArg_ParseTuple(args,"sis",&str1,&int1,&str2);
  if(ok) {
    APIEntry();
    ok = ExecutiveSetTitle(str1,int1,str2);
    APIExit();
  }
  return(APIStatus(ok));
}

static PyObject *CmdExportCoords(PyObject *self, 	PyObject *args)
{
  void *result;
  char *str1;
  int int1;
  PyObject *py_result = Py_None;
  int ok = false;
  ok = PyArg_ParseTuple(args,"si",&str1,&int1);
  if(ok) {
    APIEntry();
    result = ExportCoordsExport(str1,int1,0);
    APIExit();
    if(result) 
      py_result = PyCObject_FromVoidPtr(result,(void(*)(void*))ExportCoordsFree);
  }
  return(APIAutoNone(py_result));
}

static PyObject *CmdImportCoords(PyObject *self, 	PyObject *args)
{
  char *str1;
  int int1;
  PyObject *cObj;
  void *mmdat=NULL;
  int ok = false;
  ok = PyArg_ParseTuple(args,"siO",&str1,&int1,&cObj);
  if(ok) {
    if(PyCObject_Check(cObj))
      mmdat = PyCObject_AsVoidPtr(cObj);
    APIEntry();
    if(mmdat)
      ok = ExportCoordsImport(str1,int1,mmdat,0);
    APIExit();
  }
  return(APIStatus(ok));
}

static PyObject *CmdGetArea(PyObject *self, 	PyObject *args)
{
  char *str1;
  int int1,int2;
  OrthoLineType s1="";
  float result = -1.0;
  int ok=false;
  ok = PyArg_ParseTuple(args,"sii",&str1,&int1,&int2);
  if(ok) {
    APIEntry();
    if(str1[0]) SelectorGetTmp(str1,s1);
    result = ExecutiveGetArea(s1,int1,int2);
    if(s1[0]) SelectorFreeTmp(s1);
    APIExit();
  }
  return(Py_BuildValue("f",result));

}

static PyObject *CmdPushUndo(PyObject *self, 	PyObject *args)
{
  char *str0;
  int state;
  OrthoLineType s0="";
  int ok=false;
  ok = PyArg_ParseTuple(args,"si",&str0,&state);
  if(ok) {
    APIEntry();
    if(str0[0]) SelectorGetTmp(str0,s0);
    ok = ExecutiveSaveUndo(s0,state);
    if(s0[0]) SelectorFreeTmp(s0);
    APIExit();
  }
  return(APIStatus(ok));
}

static PyObject *CmdGetType(PyObject *self, 	PyObject *args)
{
  char *str1;
  WordType type = "";
  int ok=false;
  ok = PyArg_ParseTuple(args,"s",&str1);
  if (ok) {
    APIEntry();
    ok = ExecutiveGetType(str1,type);
    APIExit();
  } 
  if(ok) 
    return(Py_BuildValue("s",type));
  else
    return(APIStatus(ok));
}
static PyObject *CmdGetNames(PyObject *self, 	PyObject *args)
{
  int int1;
  char *vla = NULL;
  PyObject *result = Py_None;
  int ok=false;
  ok = PyArg_ParseTuple(args,"i",&int1);
  if(ok) {
    APIEntry();
    vla = ExecutiveGetNames(int1);
    APIExit();
    result = PConvStringVLAToPyList(vla);
    VLAFreeP(vla);
  }
  return(APIAutoNone(result));
}

static PyObject *CmdInvert(PyObject *self, PyObject *args)
{
  char *str0,*str1;
  int int1;
  OrthoLineType s0="",s1="";
  int ok=false;
  ok = PyArg_ParseTuple(args,"ssi",&str0,&str1,&int1);
  if(ok) {
    APIEntry();
    if(str0[0]) SelectorGetTmp(str0,s0);
    if(str1[0]) SelectorGetTmp(str1,s1);
    ok = ExecutiveInvert(s0,s1,int1);
    if(s0[0]) SelectorFreeTmp(s0);
    if(s1[0]) SelectorFreeTmp(s1);
    APIExit();
  }
  return(APIStatus(ok));
}

static PyObject *CmdTorsion(PyObject *self, PyObject *args)
{
  float float1;
  int ok=false;
  ok = PyArg_ParseTuple(args,"f",&float1);
  if (ok) {
    APIEntry();
    ok = EditorTorsion(float1);
    APIExit();
  }
  return(APIStatus(ok));
}

static PyObject *CmdUndo(PyObject *self, PyObject *args)
{
  int int1;
  int ok=false;
  ok = PyArg_ParseTuple(args,"i",&int1);
  if (ok) {
    APIEntry();
    ExecutiveUndo(int1); /* TODO STATUS */
    APIExit();
  }
  return(APIStatus(ok));
}

static PyObject *CmdMask(PyObject *self, PyObject *args)
{
  char *str1;
  int int1;
  OrthoLineType s1;

  int ok=false;
  ok = PyArg_ParseTuple(args,"si",&str1,&int1);
  if (ok) {
    APIEntry();
    SelectorGetTmp(str1,s1);
    ExecutiveMask(s1,int1); /* TODO STATUS */
    SelectorFreeTmp(s1);
    APIExit();
  }
  return(APIStatus(ok));
}

static PyObject *CmdProtect(PyObject *self, PyObject *args)
{
  char *str1;
  int int1;
  OrthoLineType s1;

  int ok=false;
  ok = PyArg_ParseTuple(args,"si",&str1,&int1);
  if (ok) {
    APIEntry();
    SelectorGetTmp(str1,s1);
    ExecutiveProtect(s1,int1); /* TODO STATUS */
    SelectorFreeTmp(s1);
    APIExit();
  }
  return(APIStatus(ok));
}


static PyObject *CmdButton(PyObject *self, 	PyObject *args)
{
  int i1,i2;
  int ok=false;
  ok = PyArg_ParseTuple(args,"ii",&i1,&i2);
  if (ok) {
    APIEntry();
    ButModeSet(i1,i2); /* TODO STATUS */
    APIExit();
  }
  return(APIStatus(ok));  
}

static PyObject *CmdFeedback(PyObject *self, 	PyObject *args)
{
  int i1,i2,result = 0;
  int ok=false;
  ok = PyArg_ParseTuple(args,"ii",&i1,&i2);
  if (ok) {
    /* NO API Entry for performance,
     *feedback (MACRO) just accesses a safe global */
    result = Feedback(i1,i2);
  }
  return Py_BuildValue("i",result);
}

static PyObject *CmdSetFeedbackMask(PyObject *self, 	PyObject *args)
{
  int i1,i2,i3;
  int ok=false;
  ok = PyArg_ParseTuple(args,"iii",&i1,&i2,&i3);
  if (ok) {
    APIEntry();
    switch(i1) { /* TODO STATUS */
    case 0: 
      FeedbackSetMask(i2,(uchar)i3);
      break;
    case 1:
      FeedbackEnable(i2,(uchar)i3);
      break;
    case 2:
      FeedbackDisable(i2,(uchar)i3);
      break;
    case 3:
      FeedbackPush();
      break;
    case 4:
      FeedbackPop();
      break;
    }
    APIExit();
  }
  return(APIStatus(ok));  
}

static PyObject *CmdFocus(PyObject *self, 	PyObject *args)
{
  /* BROKEN */
  APIEntry();
  ExecutiveFocus();
  APIExit();
  return(APIFailure());
}

static PyObject *CmdFlushNow(PyObject *self, 	PyObject *args)
{
  /* only called by the GLUT thread with unlocked API, blocked interpreter */
  if(flush_count<8) { /* prevent super-deep recursion */
    flush_count++;
    PFlushFast();
    flush_count--;
  } else {
    PRINTFB(FB_CCmd,FB_Warnings)
      " Cmd: PyMOL lagging behind API requests...\n"
      ENDFB;
  }
  return(APISuccess());  
}

static PyObject *CmdWaitQueue(PyObject *self, 	PyObject *args)
{
  /* called by non-GLUT thread with unlocked API, blocked interpreter */
  PyObject *result;
  if(OrthoCommandWaiting()||(flush_count>1)) 
    result = PyInt_FromLong(1);
  else
    result = PyInt_FromLong(0);
  return result;
}

static PyObject *CmdPaste(PyObject *dummy, PyObject *args)
{
  PyObject *list,*str;
  char *st;
  int l,a;
  int ok=false;
  ok = PyArg_ParseTuple(args,"O",&list);
  if(ok) {
    if(!list) 
      ok=false;
    else if(!PyList_Check(list)) 
      ok=false;
    else
      {
        l=PyList_Size(list);
        for(a=0;a<l;a++) {
          str = PyList_GetItem(list,a);
          if(str) {
            if(PyString_Check(str)) {
              st = PyString_AsString(str);
              APIEntry();
              OrthoPasteIn(st);
              APIExit();
            } else {
              ok = false;
            }
          }
        }
      }
  }
  return(APIStatus(ok));  
}

static PyObject *CmdGetPovRay(PyObject *dummy, PyObject *args)
{
  PyObject *result = NULL;
  char *header=NULL,*geom=NULL;
  APIEntry();
  SceneRay(0,0,1,&header,&geom);
  if(header&&geom) {
    result = Py_BuildValue("(ss)",header,geom);
  }
  VLAFreeP(header);
  VLAFreeP(geom);
  APIExit();
  return(APIAutoNone(result));
}

static PyObject *CmdGetWizard(PyObject *dummy, PyObject *args)
{
  PyObject *result;
  APIEntry();
  result = WizardGet();
  APIExit();
  if(!result)
    result=Py_None;
  return APIIncRef(result);
}

static PyObject *CmdSetWizard(PyObject *dummy, PyObject *args)
{
  
  PyObject *obj;
  int ok=false;
  ok = PyArg_ParseTuple(args,"O",&obj);
  if(ok) {
    if(!obj)
      ok=false;
    else
      {
        APIEntry();
        WizardSet(obj); /* TODO STATUS */
        APIExit();
      }
  }
  return(APIStatus(ok));  
}

static PyObject *CmdRefreshWizard(PyObject *dummy, PyObject *args)
{
  
  APIEntry();
  WizardRefresh();
  APIExit();
  return(APISuccess());  
}

static PyObject *CmdSplash(PyObject *dummy, PyObject *args)
{
  APIEntry();
  OrthoSplash();
  APIExit();
  return(APISuccess());  
}

static PyObject *CmdCls(PyObject *dummy, PyObject *args)
{
  APIEntry();
  OrthoClear();
  APIExit();
  return(APISuccess());  
}

static PyObject *CmdDump(PyObject *dummy, PyObject *args)
{
  char *str1,*str2;
  int ok=false;
  ok = PyArg_ParseTuple(args,"ss",&str1,&str2);
  if (ok) {
    APIEntry();
    ExecutiveDump(str1,str2); /* TODO STATUS */
    APIExit();
  }
  return(APIStatus(ok));  
}

static PyObject *CmdIsomesh(PyObject *self, 	PyObject *args) {
  char *str1,*str2,*str3;
  float lvl,fbuf;
  int dotFlag;
  int c,state=-1;
  OrthoLineType s1;
  int oper,frame;
  float carve;
  Object *obj,*mObj,*origObj;
  ObjectMap *mapObj;
  float mn[3] = { 0,0,0};
  float mx[3] = { 15,15,15};
  float *vert_vla = NULL;
  int ok = false;
  /* oper 0 = all, 1 = sele + buffer, 2 = vector */

  ok = PyArg_ParseTuple(args,"sisisffiif",&str1,&frame,&str2,&oper,
                   &str3,&fbuf,&lvl,&dotFlag,&state,&carve);
  if (ok) {
    APIEntry();

    origObj=ExecutiveFindObjectByName(str1);  
    if(origObj) {
      if(origObj->type!=cObjectMesh) {
        ExecutiveDelete(str1);
        origObj=NULL;
      }
    }
    
    mObj=ExecutiveFindObjectByName(str2);  
    if(mObj) {
      if(mObj->type!=cObjectMap)
        mObj=NULL;
    }
    if(mObj) {
      mapObj = (ObjectMap*)mObj;
      switch(oper) {
      case 0:
        for(c=0;c<3;c++) {
          mn[c] = mapObj->Corner[0][c];
          mx[c] = mapObj->Corner[7][c];
        }
        carve = false; /* impossible */
        break;
      case 1:
        SelectorGetTmp(str3,s1);
        ExecutiveGetExtent(s1,mn,mx,false);
        if(carve>=0.0) {
          vert_vla = ExecutiveGetVertexVLA(s1,state);
          if(fbuf<=R_SMALL4)
            fbuf = carve;
        }
        SelectorFreeTmp(s1);
        for(c=0;c<3;c++) {
          mn[c]-=fbuf;
          mx[c]+=fbuf;
        }
        break;
      }
      PRINTFB(FB_CCmd,FB_Blather)
        " Isomesh: buffer %8.3f carve %8.3f \n",fbuf,carve
        ENDFB;
      obj=(Object*)ObjectMeshFromBox((ObjectMesh*)origObj,mapObj,state,mn,mx,lvl,dotFlag,
                                     carve,vert_vla);
      if(!origObj) {
        ObjectSetName(obj,str1);
        ExecutiveManageObject((Object*)obj);
      }
      if(SettingGet(cSetting_isomesh_auto_state))
        if(obj) ObjectGotoState((ObjectMolecule*)obj,state);
      PRINTFB(FB_ObjectMesh,FB_Actions)
        " Isomesh: created \"%s\", setting level to %5.3f\n",str1,lvl
        ENDFB;
    } else {
      PRINTFB(FB_ObjectMesh,FB_Errors)
        " Isomesh: Map or brick object '%s' not found.\n",str2
        ENDFB;
      ok=false;
    }
    APIExit();
  }
  return(APIStatus(ok));  
}

static PyObject *CmdSymExp(PyObject *self, 	PyObject *args) {
  char *str1,*str2,*str3;
  OrthoLineType s1;
  float cutoff;
  Object *mObj;
  /* oper 0 = all, 1 = sele + buffer, 2 = vector */

  int ok=false;
  ok = PyArg_ParseTuple(args,"sssf",&str1,&str2,&str3,&cutoff);
  if (ok) {
    APIEntry();
    mObj=ExecutiveFindObjectByName(str2);  
    if(mObj) {
      if(mObj->type!=cObjectMolecule) {
        mObj=NULL;
        ok = false;
      }
    }
    if(mObj) {
      SelectorGetTmp(str3,s1);
      ExecutiveSymExp(str1,str2,s1,cutoff); /* TODO STATUS */
      SelectorFreeTmp(s1);
    }
    APIExit();
  }
  return(APIStatus(ok));  
}

static PyObject *CmdOverlap(PyObject *dummy, PyObject *args)
{
  char *str1,*str2;
  int state1,state2;
  float overlap = -1.0;
  float adjust;
  OrthoLineType s1,s2;
  int ok=false;
  ok = PyArg_ParseTuple(args,"ssiif",&str1,&str2,&state1,&state2,&adjust);
  if (ok) {
    APIEntry();
    SelectorGetTmp(str1,s1);
    SelectorGetTmp(str2,s2);
    overlap = ExecutiveOverlap(s1,state1,s2,state2,adjust);
    SelectorFreeTmp(s1);
    SelectorFreeTmp(s2);
    APIExit();
  }
  return(Py_BuildValue("f",overlap));
}

static PyObject *CmdDist(PyObject *dummy, PyObject *args)
{
  char *name,*str1,*str2;
  float cutoff,result=-1.0;
  int mode;
  OrthoLineType s1,s2;
  int ok=false;
  ok = PyArg_ParseTuple(args,"sssif",&name,&str1,&str2,&mode,&cutoff);
  if (ok) {
    APIEntry();
    SelectorGetTmp(str1,s1);
    SelectorGetTmp(str2,s2);
    result = ExecutiveDist(name,s1,s2,mode,cutoff);
    SelectorFreeTmp(s1);
    SelectorFreeTmp(s2);
    APIExit();
  }
  return(Py_BuildValue("f",result));
}

static PyObject *CmdBond(PyObject *dummy, PyObject *args)
{
  char *str1,*str2;
  int order,mode;
  OrthoLineType s1,s2;
  int ok=false;
  ok = PyArg_ParseTuple(args,"ssii",&str1,&str2,&order,&mode);
  if (ok) {
    APIEntry();
    SelectorGetTmp(str1,s1);
    SelectorGetTmp(str2,s2);
    ExecutiveBond(s1,s2,order,mode); /* TODO STATUS */
    SelectorFreeTmp(s1);
    SelectorFreeTmp(s2);
    APIExit();
  }
  return(APIStatus(ok));  
}

static PyObject *CmdDistance(PyObject *dummy, PyObject *args)
{
  char *str1,*str2;
  OrthoLineType s1,s2;
  float dist=-1.0;
  int ok=false;
  ok = PyArg_ParseTuple(args,"ss",&str1,&str2);
  if (ok) {
    APIEntry();
    SelectorGetTmp(str1,s1);
    SelectorGetTmp(str2,s2);
    dist = ExecutiveDistance(s1,s2);
    SelectorFreeTmp(s1);
    SelectorFreeTmp(s2);
    APIExit();
  }
  return(Py_BuildValue("f",dist));
}

static PyObject *CmdLabel(PyObject *self,   PyObject *args)
{
  char *str1,*str2;
  OrthoLineType s1;

  int ok=false;
  ok = PyArg_ParseTuple(args,"ss",&str1,&str2);
  if (ok) {
    APIEntry();
    SelectorGetTmp(str1,s1);
    ExecutiveLabel(s1,str2); /* TODO STATUS */
    SelectorFreeTmp(s1);
    APIExit();
  }
  return(APIStatus(ok));

}

static PyObject *CmdAlter(PyObject *self,   PyObject *args)
{
  char *str1,*str2;
  int i1;
  OrthoLineType s1;
  int result=0;
  int ok=false;
  ok = PyArg_ParseTuple(args,"ssi",&str1,&str2,&i1);
  if (ok) {
    APIEntry();
    SelectorGetTmp(str1,s1);
    result=ExecutiveIterate(s1,str2,i1); /* TODO STATUS */
    SelectorFreeTmp(s1);
    APIExit();
  }
  return Py_BuildValue("i",result);

}

static PyObject *CmdAlterState(PyObject *self,   PyObject *args)
{
  char *str1,*str2;
  int i1,i2;
  OrthoLineType s1;

  int ok=false;
  ok = PyArg_ParseTuple(args,"issi",&i1,&str1,&str2,&i2);
  if (ok) {
    APIEntry();
    SelectorGetTmp(str1,s1);
    ExecutiveIterateState(i1,s1,str2,i2); /* TODO STATUS */
    SelectorFreeTmp(s1);
    APIExit();
  }
  return(APIStatus(ok));

}

static PyObject *CmdCopy(PyObject *self,   PyObject *args)
{
  char *str1,*str2;
  int ok=false;
  ok = PyArg_ParseTuple(args,"ss",&str1,&str2);
  if (ok) {
    APIEntry();
    ExecutiveCopy(str1,str2); /* TODO STATUS */
    APIExit();
  }
  return(APIStatus(ok));
}

static PyObject *CmdRecolor(PyObject *self,   PyObject *args)
{
  char *str1;
  OrthoLineType s1;
  int ok=true;
  int rep=-1;
  ok = PyArg_ParseTuple(args,"si",&str1,&rep);
  PRINTFD(FB_CCmd)
    " CmdRebuild: called with %s.\n",str1
    ENDFD;

  if (ok) {
    APIEntry();
    if(WordMatch(str1,"all",true)<0)
      ExecutiveInvalidateRep(str1,rep,cRepInvColor);
    else {
      SelectorGetTmp(str1,s1);
      ExecutiveInvalidateRep(s1,rep,cRepInvColor);
      SelectorFreeTmp(s1); 
    }
    APIExit();
  } else {
    ok = -1; /* special error convention */
  }
  return(APIStatus(ok));
}

static PyObject *CmdRebuild(PyObject *self,   PyObject *args)
{
  char *str1;
  OrthoLineType s1;
  int ok=true;
  int rep=-1;
  ok = PyArg_ParseTuple(args,"si",&str1,&rep);
  PRINTFD(FB_CCmd)
    " CmdRebuild: called with %s.\n",str1
    ENDFD;

  if (ok) {
    APIEntry();
    if(WordMatch(str1,"all",true)<0)
      ExecutiveRebuildAll();
    else {
      SelectorGetTmp(str1,s1);
      ExecutiveInvalidateRep(s1,rep,cRepInvAll);
      SelectorFreeTmp(s1); 
    }
    APIExit();
  } else {
    ok = -1; /* special error convention */
  }
  return(APIStatus(ok));
}

static PyObject *CmdResetRate(PyObject *dummy, PyObject *args)
{
  APIEntry();
  ButModeResetRate();
  APIExit();
  return(APISuccess());
}

static PyObject *CmdReady(PyObject *dummy, PyObject *args)
{
  return(APIStatus(PyMOLReady));
}

static PyObject *CmdMem(PyObject *dummy, PyObject *args)
{
  MemoryDebugDump();
  return(APISuccess());
}

static PyObject *CmdRunPyMOL(PyObject *dummy, PyObject *args)
{
  if(run_only_once) {
    run_only_once=false;
#ifdef _PYMOL_MODULE
    was_main();
#endif
  }
  return(APISuccess());
}

static PyObject *CmdCountStates(PyObject *dummy, PyObject *args)
{
  char *str1;
  OrthoLineType s1;
  int ok=false;
  ok = PyArg_ParseTuple(args,"s",&str1);
  if (ok) {
    APIEntry();
    SelectorGetTmp(str1,s1);
    ok = ExecutiveCountStates(s1);
    SelectorFreeTmp(s1); 
    APIExit();
  } else {
    ok = -1; /* special error convention */
  }
  return(APIStatus(ok));
}

static PyObject *CmdIdentify(PyObject *dummy, PyObject *args)
{
  char *str1;
  OrthoLineType s1;
  int mode;
  int a,l=0;
  PyObject *result = Py_None;
  PyObject *tuple;
  int *iVLA=NULL,*i;
  ObjectMolecule **oVLA=NULL,**o;
  int ok=false;
  ok = PyArg_ParseTuple(args,"si",&str1,&mode);
  if (ok) {
    APIEntry();
    SelectorGetTmp(str1,s1);
    if(!mode) {
      iVLA=ExecutiveIdentify(s1,mode);
    } else {
      l = ExecutiveIdentifyObjects(s1,mode,&iVLA,&oVLA);
    }
    SelectorFreeTmp(s1);
    APIExit();
    if(iVLA) {
      if(!mode) {
        result=PConvIntVLAToPyList(iVLA);
      } else { /* object mode */
        result=PyList_New(l);
        i = iVLA;
        o = oVLA;
        for(a=0;a<l;a++) {
          tuple = PyTuple_New(2);      
          PyTuple_SetItem(tuple,1,PyInt_FromLong(*(i++))); 
          PyTuple_SetItem(tuple,0,PyString_FromString((*(o++))->Obj.Name));
          PyList_SetItem(result,a,tuple);
        }
      }
    } else {
      result = PyList_New(0);
    }
  }
  VLAFreeP(iVLA);
  VLAFreeP(oVLA);
  return(APIAutoNone(result));
}

static PyObject *CmdIndex(PyObject *dummy, PyObject *args)
{
  char *str1;
  OrthoLineType s1;
  int mode;
  PyObject *result = Py_None;
  PyObject *tuple = Py_None;
  int *iVLA=NULL;
  int l;
  int *i;
  ObjectMolecule **o,**oVLA=NULL;
  int a;

  int ok=false;
  ok = PyArg_ParseTuple(args,"si",&str1,&mode);
  if (ok) {
    APIEntry();
    SelectorGetTmp(str1,s1);
    l = ExecutiveIndex(s1,mode,&iVLA,&oVLA);
    SelectorFreeTmp(s1);
    APIExit();
    if(iVLA) {
      result=PyList_New(l);
      i = iVLA;
      o = oVLA;
      for(a=0;a<l;a++) {
        tuple = PyTuple_New(2);      
        PyTuple_SetItem(tuple,1,PyInt_FromLong(*(i++)+1)); /* +1 for index */
        PyTuple_SetItem(tuple,0,PyString_FromString((*(o++))->Obj.Name));
        PyList_SetItem(result,a,tuple);
      }
    } else {
      result = PyList_New(0);
    }
    VLAFreeP(iVLA);
    VLAFreeP(oVLA);
  }
  return(APIAutoNone(result));
}


static PyObject *CmdFindPairs(PyObject *dummy, PyObject *args)
{
  char *str1,*str2;
  int state1,state2;
  float cutoff;
  float angle;
  int mode;
  OrthoLineType s1,s2;
  PyObject *result = Py_None;
  PyObject *tuple = Py_None;
  PyObject *tuple1 = Py_None;
  PyObject *tuple2 = Py_None;
  int *iVLA=NULL;
  int l;
  int *i;
  ObjectMolecule **o,**oVLA=NULL;
  int a;
  
  int ok=false;
  ok = PyArg_ParseTuple(args,"ssiiiff",&str1,&str2,&state1,&state2,&mode,&cutoff,&angle);
  if (ok) {
    APIEntry();
    SelectorGetTmp(str1,s1);
    SelectorGetTmp(str2,s2);
    l = ExecutivePairIndices(s1,s2,state1,state2,mode,cutoff,angle,&iVLA,&oVLA);
    SelectorFreeTmp(s1);
    SelectorFreeTmp(s2);
    APIExit();
    
    if(iVLA&&oVLA) {
      result=PyList_New(l);
      i = iVLA;
      o = oVLA;
      for(a=0;a<l;a++) {
        tuple1 = PyTuple_New(2);      
        PyTuple_SetItem(tuple1,0,PyString_FromString((*(o++))->Obj.Name));
        PyTuple_SetItem(tuple1,1,PyInt_FromLong(*(i++)+1)); /* +1 for index */
        tuple2 = PyTuple_New(2);
        PyTuple_SetItem(tuple2,0,PyString_FromString((*(o++))->Obj.Name));
        PyTuple_SetItem(tuple2,1,PyInt_FromLong(*(i++)+1)); /* +1 for index */
        tuple = PyTuple_New(2);
        PyTuple_SetItem(tuple,0,tuple1);
        PyTuple_SetItem(tuple,1,tuple2);
        PyList_SetItem(result,a,tuple);
      }
    } else {
      result = PyList_New(0);
    }
    VLAFreeP(iVLA);
    VLAFreeP(oVLA);
  }
  return(APIAutoNone(result));
}

static PyObject *CmdSystem(PyObject *dummy, PyObject *args)
{
  char *str1;
  int ok=false;
  ok = PyArg_ParseTuple(args,"s",&str1);
  if (ok) {
    APIEntry();
    ok = system(str1);
    APIExit();
  }
  return(APIStatus(ok));
}

static PyObject *CmdGetFeedback(PyObject *dummy, PyObject *args)
{
  PyObject *result = NULL;
  OrthoLineType buffer;
  int ok;

  if(PyMOLTerminating) { /* try to bail */
#ifdef WIN32
	abort();
#endif
    exit(0);
  }
  ok = OrthoFeedbackOut(buffer); 
  if(ok) result = Py_BuildValue("s",buffer);
  return(APIAutoNone(result));
}

static PyObject *CmdGetPDB(PyObject *dummy, PyObject *args)
{
  char *str1;
  char *pdb = NULL;
  int state;
  OrthoLineType s1 = "";
  PyObject *result = NULL;
  int ok=false;
  ok = PyArg_ParseTuple(args,"si",&str1,&state);
  if (ok) {
    APIEntry();
    SelectorGetTmp(str1,s1);
    pdb=ExecutiveSeleToPDBStr(s1,state,true);
    SelectorFreeTmp(s1);
    APIExit();
    if(pdb) result = Py_BuildValue("s",pdb);
    FreeP(pdb);
  }
  return(APIAutoNone(result));
}

static PyObject *CmdGetModel(PyObject *dummy, PyObject *args)
{
  char *str1;
  int state;
  OrthoLineType s1;
  PyObject *result = NULL;
  int ok=false;
  ok = PyArg_ParseTuple(args,"si",&str1,&state);
  if (ok) {
    APIEntry();
    SelectorGetTmp(str1,s1);
    result=ExecutiveSeleToChemPyModel(s1,state);
    SelectorFreeTmp(s1);
    APIExit();
  }
  return(APIAutoNone(result));
}

static PyObject *CmdCreate(PyObject *dummy, PyObject *args)
{
  char *str1,*str2;
  int target,source;
  OrthoLineType s1;
  int ok=false;
  ok = PyArg_ParseTuple(args,"ssii",&str1,&str2,&source,&target);
  if (ok) {
    APIEntry();
    SelectorGetTmp(str2,s1);
    ExecutiveSeleToObject(str1,s1,source,target); /* TODO STATUS */
    SelectorFreeTmp(s1);
    APIExit();
  }
  return(APIStatus(ok));
}


static PyObject *CmdOrient(PyObject *dummy, PyObject *args)
{
  Matrix33d m;
  char *str1;
  OrthoLineType s1;

  int ok=false;
  ok = PyArg_ParseTuple(args,"s",&str1);
  if (ok) {
    APIEntry();
    SelectorGetTmp(str1,s1);
    if(ExecutiveGetMoment(s1,m))
      ExecutiveOrient(s1,m); /* TODO STATUS */
    else
      ok=false;
    SelectorFreeTmp(s1); 
    APIExit();
  }
  return(APIStatus(ok));
}

static PyObject *CmdFitPairs(PyObject *dummy, PyObject *args)
{
  PyObject *list;
  WordType *word = NULL;
  int ln=0;
  int a;
  PyObject *result = NULL;
  float valu;
  int ok=false;
  ok = PyArg_ParseTuple(args,"O",&list);
  if(ok) {
    ln = PyObject_Length(list);
    if(ln) {
      if(ln&0x1)
        ok=ErrMessage("FitPairs","must supply an even number of selections.");
    } else ok=false;
    
    if(ok) {
      word = Alloc(WordType,ln);
      
      a=0;
      while(a<ln) {
        SelectorGetTmp(PyString_AsString(PySequence_GetItem(list,a)),word[a]);
        a++;
      }
      APIEntry();
      valu = ExecutiveRMSPairs(word,ln/2,2);
      APIExit();
      result=Py_BuildValue("f",valu);
      for(a=0;a<ln;a++)
        SelectorFreeTmp(word[a]);
      FreeP(word);
    }
  }
  return APIAutoNone(result);
}

static PyObject *CmdIntraFit(PyObject *dummy, PyObject *args)
{
  char *str1;
  int state;
  int mode;
  OrthoLineType s1;
  float *fVLA;
  PyObject *result=Py_None;
  int ok=false;
  ok = PyArg_ParseTuple(args,"sii",&str1,&state,&mode);
  if(state<0) state=0;
  if (ok) {
    APIEntry();
    SelectorGetTmp(str1,s1);
    fVLA=ExecutiveRMSStates(s1,state,mode);
    SelectorFreeTmp(s1);
    APIExit();
    if(fVLA) {
      result=PConvFloatVLAToPyList(fVLA);
      VLAFreeP(fVLA);
    }
  }
  return APIAutoNone(result);
}

static PyObject *CmdFit(PyObject *dummy, PyObject *args)
{
  char *str1,*str2;
  int mode;
  OrthoLineType s1,s2;
  PyObject *result;
  float tmp_result = -1.0;
  
  int ok=false;
  ok = PyArg_ParseTuple(args,"ssi",&str1,&str2,&mode);
  if (ok) {
    APIEntry();
    SelectorGetTmp(str1,s1);
    SelectorGetTmp(str2,s2);
    tmp_result=ExecutiveRMS(s1,s2,mode,0.0);
    SelectorFreeTmp(s1);
    SelectorFreeTmp(s2);
    APIExit();
  }
  result=Py_BuildValue("f",tmp_result);
  return result;
}

static PyObject *CmdUpdate(PyObject *dummy, PyObject *args)
{
  char *str1,*str2;
  int int1,int2;
  OrthoLineType s1,s2;
  int ok=false;
  ok = PyArg_ParseTuple(args,"ssii",&str1,&str2,&int1,&int2);
  if (ok) {
    APIEntry();
    SelectorGetTmp(str1,s1);
    SelectorGetTmp(str2,s2);
    ExecutiveUpdateCmd(s1,s2,int1,int2); /* TODO STATUS */
    SelectorFreeTmp(s1);
    SelectorFreeTmp(s2);
    APIExit();
  }
  return(APIStatus(ok));
}

static PyObject *CmdDirty(PyObject *self, 	PyObject *args)
{
  PRINTFD(FB_CCmd)
    " CmdDirty: called.\n"
    ENDFD;
  APIEntry();
  OrthoDirty();
  APIExit();
  return(APISuccess());
}

static PyObject *CmdGetDihe(PyObject *self, 	PyObject *args)
{
  char *str1,*str2,*str3,*str4;
  float result;
  int int1;
  OrthoLineType s1,s2,s3,s4;
  int ok=false;
  ok = PyArg_ParseTuple(args,"ssssi",&str1,&str2,&str3,&str4,&int1);
  
  if (ok) {
    APIEntry();
    SelectorGetTmp(str1,s1);
    SelectorGetTmp(str2,s2);
    SelectorGetTmp(str3,s3);
    SelectorGetTmp(str4,s4);
    ok = ExecutiveGetDihe(s1,s2,s3,s4,&result,int1);
    SelectorFreeTmp(s1);
    SelectorFreeTmp(s2);
    SelectorFreeTmp(s3);
    SelectorFreeTmp(s4);
    APIExit();
  }
  
  if(ok) {
    return(Py_BuildValue("f",result));
  } else {
    return APIFailure();
  }
}

static PyObject *CmdSetDihe(PyObject *self, 	PyObject *args)
{
  char *str1,*str2,*str3,*str4;
  float float1;
  int int1;
  OrthoLineType s1,s2,s3,s4;
  int ok=false;
  ok = PyArg_ParseTuple(args,"ssssfi",&str1,&str2,&str3,&str4,&float1,&int1);
  if (ok) {
    APIEntry();
    SelectorGetTmp(str1,s1);
    SelectorGetTmp(str2,s2);
    SelectorGetTmp(str3,s3);
    SelectorGetTmp(str4,s4);
    ok = ExecutiveSetDihe(s1,s2,s3,s4,float1,int1);
    SelectorFreeTmp(s1);
    SelectorFreeTmp(s2);
    SelectorFreeTmp(s3);
    SelectorFreeTmp(s4);
    APIExit();
  }
  return(APIStatus(ok));
}

static PyObject *CmdDo(PyObject *self, 	PyObject *args)
{
  char *str1;

  int ok=false;
  ok = PyArg_ParseTuple(args,"s",&str1);
  if (ok) {
    APIEntry();
    if(str1[0]!='_') { /* suppress internal call-backs */
      if(strncmp(str1,"cmd._",5)) {
        OrthoAddOutput("PyMOL>");
        OrthoAddOutput(str1);
        OrthoNewLine(NULL);
        if(WordMatch(str1,"quit",true)==0) /* don't log quit */
          PLog(str1,cPLog_pml);
      }
      PParse(str1);
    } else if(str1[1]==' ') { /* "_ command" suppresses echoing of command, but it is still logged */
      if(WordMatch(str1+2,"quit",true)==0) /* don't log quit */
        PLog(str1+2,cPLog_pml);
      PParse(str1+2);    
    } else {
      PParse(str1);
    }
    APIExit();
  }
  return(APIStatus(ok));
}

static PyObject *CmdRock(PyObject *self, PyObject *args)
{
  APIEntry();
  ControlRock(-1);
  APIExit();
  return APISuccess();
}

static PyObject *CmdGetMoment(PyObject *self, 	PyObject *args)
{
  Matrix33d m;
  PyObject *result;

  char *str1;
  int ok=false;
  ok = PyArg_ParseTuple(args,"s",&str1);
  if (ok) {
    APIEntry();
    ExecutiveGetMoment(str1,m);
    APIExit();
  }
  result = Py_BuildValue("(ddd)(ddd)(ddd)", 
								 m[0][0],m[0][1],m[0][2],
								 m[1][0],m[1][1],m[1][2],
								 m[2][0],m[2][1],m[2][2]);
  return result;
}

static PyObject *CmdGetSetting(PyObject *self, 	PyObject *args)
{
  PyObject *result = Py_None;
  char *str1;
  float value;
  int ok=false;
  ok = PyArg_ParseTuple(args,"s",&str1);
  if (ok) {
    APIEntry();
    value=SettingGetNamed(str1);
    APIExit();
    result = Py_BuildValue("f", SettingGetNamed(str1));
  }
  return APIAutoNone(result);
}

static PyObject *CmdGetSettingTuple(PyObject *self, 	PyObject *args)
{
  PyObject *result = Py_None;
  int int1,int2;
  char *str1;
  int ok = false;
  ok = PyArg_ParseTuple(args,"isi",&int1,&str1,&int2); /* setting, object, state */
  if (ok) {
    APIEntry();
    result =  ExecutiveGetSettingTuple(int1,str1,int2);
    APIExit();
  }
  return APIAutoNone(result);
}

static PyObject *CmdGetSettingText(PyObject *self, 	PyObject *args)
{
  PyObject *result = Py_None;
  int int1,int2;
  char *str1;
  int ok=false;
  ok = PyArg_ParseTuple(args,"isi",&int1,&str1,&int2); /* setting, object, state */
  if (ok) {
    APIEntry();
    result =  ExecutiveGetSettingText(int1,str1,int2);
    APIExit();
  }
  return APIAutoNone(result);
}

static PyObject *CmdExportDots(PyObject *self, 	PyObject *args)
{
  PyObject *result=NULL;
  PyObject *cObj;
  ExportDotsObj *obj;
  char *str1;
  int int1;
  
  int ok=false;
  ok = PyArg_ParseTuple(args,"si",&str1,&int1);
  if (ok) {
    APIEntry();
    obj = ExportDots(str1,int1-1);
    APIExit();
    if(obj) 
      {
        cObj = PyCObject_FromVoidPtr(obj,(void(*)(void*))ExportDeleteMDebug);
        if(cObj) {
          result = Py_BuildValue("O",cObj); /* I think this */
          Py_DECREF(cObj); /* transformation is unnecc. */
        } 
      }
  }
  return APIAutoNone(result);
}

static PyObject *CmdSetFrame(PyObject *self, PyObject *args)
{
  int mode,frm;
  int ok=false;
  ok = PyArg_ParseTuple(args,"ii",&mode,&frm);
  if (ok) {
    APIEntry();
    SceneSetFrame(mode,frm);
    APIExit();
  }
  return(APIStatus(ok));
}

static PyObject *CmdFrame(PyObject *self, PyObject *args)
{
  int frm;
  int ok=false;
  ok = PyArg_ParseTuple(args,"i",&frm);
  frm--;
  if(frm<0) frm=0;
  if (ok) {
    APIEntry();
    SceneSetFrame(0,frm);
    APIExit();
  }
  return(APIStatus(ok));
}

static PyObject *CmdStereo(PyObject *self, PyObject *args)
{
  int i1;
  int ok=false;
  
  if(StereoCapable) {
    ok = PyArg_ParseTuple(args,"i",&i1);
    if (ok) {
      APIEntry();
      ExecutiveStereo(i1); /* TODO STATUS */
      APIExit();
    }
  }
  return APIStatus(ok);
}

static PyObject *CmdReset(PyObject *self, PyObject *args)
{
  int cmd;
  int ok=false;
  char *obj;
  ok = PyArg_ParseTuple(args,"is",&cmd,&obj);
  if (ok) {
    APIEntry();
    ok = ExecutiveReset(cmd,obj); 
    APIExit();
  }
  return(APIStatus(ok));
}

static PyObject *CmdSetMatrix(PyObject *self, 	PyObject *args)
{
  float m[16];
  int ok=false;
  ok = PyArg_ParseTuple(args,"ffffffffffffffff",
						 &m[0],&m[1],&m[2],&m[3],
						 &m[4],&m[5],&m[6],&m[7],
						 &m[8],&m[9],&m[10],&m[11],
						 &m[12],&m[13],&m[14],&m[15]);
  if (ok) {
    APIEntry();
    SceneSetMatrix(m); /* TODO STATUS */
    APIExit();
  }
  return(APIStatus(ok));
}

static PyObject *CmdGetMinMax(PyObject *self, 	PyObject *args)
{
  float mn[3],mx[3];
  char *str1;
  int state;
  OrthoLineType s1;
  PyObject *result = Py_None;
  int flag;

  int ok=false;
  ok = PyArg_ParseTuple(args,"si",&str1,&state); /* state currently ignored */
  if (ok) {
    APIEntry();
    SelectorGetTmp(str1,s1);
    flag = ExecutiveGetExtent(s1,mn,mx,false);
    SelectorFreeTmp(s1);
    if(flag) 
      result = Py_BuildValue("[[fff],[fff]]", 
                             mn[0],mn[1],mn[2],
                             mx[0],mx[1],mx[2]);
    else 
      result = Py_BuildValue("[[fff],[fff]]", 
                             -0.5,-0.5,-0.5,
                             0.5,0.5,0.5);
    
    APIExit();
  }
  return APIAutoNone(result);
}

static PyObject *CmdGetMatrix(PyObject *self, 	PyObject *args)
{
  float *f;
  PyObject *result;
  
  APIEntry();
  f=SceneGetMatrix();
  APIExit();
  result = Py_BuildValue("ffffffffffffffff", 
								 f[0],f[1],f[2],f[3],
								 f[4],f[5],f[6],f[7],
								 f[8],f[9],f[10],f[11],
								 f[12],f[13],f[14],f[15]
								 );
  return result;
}

static PyObject *CmdMDo(PyObject *self, 	PyObject *args)
{
  char *cmd;
  int frame;
  int append;
  int ok=false;
  ok = PyArg_ParseTuple(args,"isi",&frame,&cmd,&append);
  if (ok) {
    APIEntry();
    if(append) {
      MovieAppendCommand(frame,cmd);
    } else {
      MovieSetCommand(frame,cmd);
    }
    APIExit();
  }
  return(APIStatus(ok));
}

static PyObject *CmdMPlay(PyObject *self, 	PyObject *args)
{
  int cmd;
  int ok=false;
  ok = PyArg_ParseTuple(args,"i",&cmd);
  if (ok) {
    APIEntry();
    MoviePlay(cmd);
    APIExit();
  }
  return(APIStatus(ok));
}

static PyObject *CmdMMatrix(PyObject *self, 	PyObject *args)
{
  int cmd;
  int ok=false;
  ok = PyArg_ParseTuple(args,"i",&cmd);
  if (ok) {
    APIEntry();
    MovieMatrix(cmd);
    APIExit();
  }
  return(APIStatus(ok));
}

static PyObject *CmdMClear(PyObject *self, 	PyObject *args)
{
  APIEntry();
  MovieClearImages();
  APIExit();
  return(APISuccess());
}

static PyObject *CmdRefresh(PyObject *self, 	PyObject *args)
{
  APIEntry();
  ExecutiveDrawNow(); /* TODO STATUS */
  APIExit();
  return(APISuccess());
}

static PyObject *CmdRefreshNow(PyObject *self, 	PyObject *args)
{
  APIEntry();
  ExecutiveDrawNow(); /* TODO STATUS */
  MainRefreshNow();
  APIExit();
  return(APISuccess());
}

static PyObject *CmdPNG(PyObject *self, 	PyObject *args)
{
  char *str1;
  int ok=false;
  ok = PyArg_ParseTuple(args,"s",&str1);
  if (ok) {
    APIEntry();
    ExecutiveDrawNow();		 /* TODO STATUS */
    ScenePNG(str1);
    APIExit();
  }
  return(APIStatus(ok));
}

static PyObject *CmdMPNG(PyObject *self, 	PyObject *args)
{
  char *str1;
  int ok=false;
  ok = PyArg_ParseTuple(args,"s",&str1);
  if (ok) {
    APIEntry();
    MoviePNG(str1,(int)SettingGet(cSetting_cache_frames));
    /* TODO STATUS */
    APIExit();
  }
  return(APIStatus(ok));
}

static PyObject *CmdMSet(PyObject *self, 	PyObject *args)
{
  char *str1;
  int ok=false;
  ok = PyArg_ParseTuple(args,"s",&str1);
  if (ok) {
    APIEntry();
    MovieSequence(str1);
    SceneCountFrames();
    APIExit();
  }
  return(APIStatus(ok));
}

static PyObject *CmdViewport(PyObject *self, 	PyObject *args)
{
  int w,h;
  int ok=false;
  ok = PyArg_ParseTuple(args,"ii",&w,&h);
  if(ok) {
    if((w>0)&&(h>0)) {
      if(w<10) w=10;
      if(h<10) h=10;
      if(SettingGet(cSetting_internal_gui)) {
        if(!SettingGet(cSetting_full_screen))
           w+=(int)SettingGet(cSetting_internal_gui_width);
      }
      if(SettingGet(cSetting_internal_feedback)) {
        if(!SettingGet(cSetting_full_screen))
          h+=(int)(SettingGet(cSetting_internal_feedback)-1)*cOrthoLineHeight +
            cOrthoBottomSceneMargin;
      }
    } else {
      w=-1;
      h=-1;
    }
    APIEntry();
    MainDoReshape(w,h); /* should be moved into Executive */
    APIExit();
  }
  return(APIStatus(ok));
}

static PyObject *CmdFlag(PyObject *self, 	PyObject *args)
{
  char *str1;
  int flag;
  int action;
  OrthoLineType s1;
  int ok=false;
  ok = PyArg_ParseTuple(args,"isi",&flag,&str1,&action);
  if (ok) {
    APIEntry();
    SelectorGetTmp(str1,s1);
    ExecutiveFlag(flag,s1,action);
    SelectorFreeTmp(s1);
    APIExit();
  }
  return(APIStatus(ok));
}

static PyObject *CmdColor(PyObject *self, 	PyObject *args)
{
  char *str1,*color;
  int flags;
  OrthoLineType s1;
  int ok=false;
  ok = PyArg_ParseTuple(args,"ssi",&color,&str1,&flags);
  if (ok) {
    APIEntry();
    SelectorGetTmp(str1,s1);
    ExecutiveColor(s1,color,flags);
    SelectorFreeTmp(s1);
    APIExit();
  }
  return(APIStatus(ok));
}

static PyObject *CmdColorDef(PyObject *self, 	PyObject *args)
{
  char *color;
  float v[3];
  int ok=false;
  ok = PyArg_ParseTuple(args,"sfff",&color,v,v+1,v+2);
  if (ok) {
    APIEntry();
    ColorDef(color,v);
    APIExit();
  }
  return(APIStatus(ok));
}

static PyObject *CmdRay(PyObject *self, 	PyObject *args)
{
  int w,h,mode;

  int ok=false;
  ok = PyArg_ParseTuple(args,"iii",&w,&h,&mode);
  if (ok) {
    APIEntry();
    if(mode<0)
      mode=(int)SettingGet(cSetting_ray_default_renderer);
    ExecutiveRay(w,h,mode); /* TODO STATUS */
    APIExit();
  }
  return(APIStatus(ok));
}

static PyObject *CmdClip(PyObject *self, 	PyObject *args)
{
  char *sname;
  float dist;
  int ok=false;
  ok = PyArg_ParseTuple(args,"sf",&sname,&dist);
  if (ok) {
    APIEntry();
    switch(sname[0]) { /* TODO STATUS */
    case 'n':
      SceneClip(0,dist);
      break;
    case 'f':
      SceneClip(1,dist);
      break;
    case 'm':
      SceneClip(2,dist);
      break;
    case 's':
      SceneClip(3,dist);
      break;
    }
    APIExit();
  }
  return(APIStatus(ok));

}

static PyObject *CmdMove(PyObject *self, 	PyObject *args)
{
  char *sname;
  float dist;
  int ok=false;
  ok = PyArg_ParseTuple(args,"sf",&sname,&dist);
  if (ok) {
    APIEntry();
    switch(sname[0]) { /* TODO STATUS */
    case 'x':
      SceneTranslate(dist,0.0,0.0);
      break;
    case 'y':
      SceneTranslate(0.0,dist,0.0);
      break;
    case 'z':
      SceneTranslate(0.0,0.0,dist);
      break;
    }
    APIExit();
  }
  return(APIStatus(ok));
}

static PyObject *CmdTurn(PyObject *self, 	PyObject *args)
{
  char *sname;
  float angle;
  int ok=false;
  ok = PyArg_ParseTuple(args,"sf",&sname,&angle);
  if (ok) {
    APIEntry();
    switch(sname[0]) { /* TODO STATUS */
    case 'x':
      SceneRotate(angle,1.0,0.0,0.0);
      break;
    case 'y':
      SceneRotate(angle,0.0,1.0,0.0);
      break;
    case 'z':
      SceneRotate(angle,0.0,0.0,1.0);
      break;
    }
    APIExit();
  }
  return(APIStatus(ok));
}

static PyObject *CmdLegacySet(PyObject *self, 	PyObject *args)
{
  char *sname, *value;
  int ok=false;
  ok = PyArg_ParseTuple(args,"ss",&sname,&value);
  if (ok) {
    APIEntry();
    ok = SettingSetNamed(sname,value); 
    APIExit();
  }
  return(APIStatus(ok));
}

static PyObject *CmdSet(PyObject *self, 	PyObject *args)
{
  int index;
  int tmpFlag=false;
  PyObject *value;
  char *str3;
  int state;
  int quiet;
  int updates;
  OrthoLineType s1;
  int ok=false;
  ok = PyArg_ParseTuple(args,"iOsiii",&index,&value,&str3,&state,&quiet,&updates);
  s1[0]=0;
  if (ok) {
    APIEntry();
    if(!strcmp(str3,"all")) {
      strcpy(s1,str3);
    } else if(str3[0]!=0) {
      tmpFlag=true;
      SelectorGetTmp(str3,s1);
    }
    ok = ExecutiveSetSetting(index,value,s1,state,quiet,updates);
    if(tmpFlag) 
      SelectorFreeTmp(s1);
    APIExit();
  }
  return(APIStatus(ok));
}

static PyObject *CmdGet(PyObject *self, 	PyObject *args)
{
  float f;
  char *sname;
  PyObject *result = Py_None;

  int ok=false;
  ok = PyArg_ParseTuple(args,"s",&sname);
  if (ok) {
    APIEntry();
    f=SettingGetNamed(sname);
    APIExit();
    result = Py_BuildValue("f", f);
  }
  return APIAutoNone(result);
}

static PyObject *CmdDelete(PyObject *self, 	PyObject *args)
{
  char *sname;

  int ok=false;
  ok = PyArg_ParseTuple(args,"s",&sname);
  if (ok) {
    APIEntry();
    ExecutiveDelete(sname); /* TODO STATUS */
    APIExit();
  }
  return(APIStatus(ok));
}

static PyObject *CmdCartoon(PyObject *self, 	PyObject *args)
{
  char *sname;
  int type;
  OrthoLineType s1;
  int ok=false;
  ok = PyArg_ParseTuple(args,"si",&sname,&type);
  if (ok) {
    APIEntry();
    SelectorGetTmp(sname,s1);
    ExecutiveCartoon(type,s1); /* TODO STATUS */
    SelectorFreeTmp(s1);
    APIExit();
  }
  return(APIStatus(ok));
}

static PyObject *CmdShowHide(PyObject *self, 	PyObject *args)
{
  char *sname;
  int rep;
  int state;
  OrthoLineType s1;
  int ok=false;
  ok = PyArg_ParseTuple(args,"sii",&sname,&rep,&state);
  if (ok) { /* TODO STATUS */
    APIEntry();
    if(sname[0]=='@') {
      ExecutiveSetAllVisib(state);
    } else {
      SelectorGetTmp(sname,s1);
      ExecutiveSetRepVisib(s1,rep,state);
      SelectorFreeTmp(s1);
    }
    APIExit();
  }
  return(APIStatus(ok));
}

static PyObject *CmdOnOff(PyObject *self, 	PyObject *args)
{
  char *name;
  int state;
  int ok=false;
  ok = PyArg_ParseTuple(args,"si",&name,&state);
  if (ok) { /* TODO STATUS */
    APIEntry();
    ExecutiveSetObjVisib(name,state);
    APIExit();
  }
  return(APIStatus(ok));
}

static PyObject *CmdQuit(PyObject *self, 	PyObject *args)
{
  APIEntry();
  PExit(EXIT_SUCCESS);
  APIExit();
  return(APISuccess());
}
static PyObject *CmdFullScreen(PyObject *self,PyObject *args)
{
  int flag = 0;
  int ok=false;
  ok = PyArg_ParseTuple(args,"i",&flag);
  if (ok) {
    APIEntry();
    ExecutiveFullScreen(flag); /* TODO STATUS */
    APIExit();
  }
  return(APIStatus(ok));
}
static PyObject *CmdSelect(PyObject *self, PyObject *args)
{
  char *sname,*sele;
  int quiet;
  int ok=false;
  ok = PyArg_ParseTuple(args,"ssi",&sname,&sele,&quiet);
  if (ok) {
    APIEntry();
    ok = SelectorCreate(sname,sele,NULL,quiet,NULL);
    SceneDirty();
    APIExit();
  } else {
    ok=-1;
  }
  return APIStatus(ok);
}

static PyObject *CmdFinishObject(PyObject *self, PyObject *args)
{
  char *oname;
  Object *origObj = NULL;

  int ok=false;
  ok = PyArg_ParseTuple(args,"s",&oname);

  if (ok) {
    APIEntry();
    origObj=ExecutiveFindObjectByName(oname);
    if(origObj) 
      ExecutiveUpdateObjectSelection(origObj); /* TODO STATUS */
    else
      ok=false;
    APIExit();
  }
  return(APIStatus(ok));
}

static PyObject *CmdLoadObject(PyObject *self, PyObject *args)
{
  char *oname;
  PyObject *model;
  Object *origObj = NULL,*obj;
  OrthoLineType buf;
  int frame,type;
  int finish,discrete;
  int ok=false;
  ok = PyArg_ParseTuple(args,"sOiiii",&oname,&model,&frame,&type,&finish,&discrete);
  buf[0]=0;
  if (ok) {
    APIEntry();
    origObj=ExecutiveFindObjectByName(oname);
    
    /* TODO check for existing object of wrong type */
    
    switch(type) {
    case cLoadTypeChemPyModel:
      if(origObj)
        if(origObj->type!=cObjectMolecule) {
          ExecutiveDelete(oname);
          origObj=NULL;
        }
      PBlock(); /*PBlockAndUnlockAPI();*/
      obj=(Object*)ObjectMoleculeLoadChemPyModel((ObjectMolecule*)
                                                 origObj,model,frame,discrete);
      PUnblock(); /*PLockAPIAndUnblock();*/
      if(!origObj) {
        if(obj) {
          ObjectSetName(obj,oname);
          ExecutiveManageObject(obj);
          if(frame<0)
            frame = ((ObjectMolecule*)obj)->NCSet-1;
          sprintf(buf," CmdLoad: ChemPy-model loaded into object \"%s\", state %d.\n",
                  oname,frame+1);		  
        }
      } else if(origObj) {
        if(finish)
          ExecutiveUpdateObjectSelection(origObj);
        if(frame<0)
          frame = ((ObjectMolecule*)origObj)->NCSet-1;
        sprintf(buf," CmdLoad: ChemPy-model appended into object \"%s\", state %d.\n",
                oname,frame+1);
      }
      break;
    case cLoadTypeChemPyBrick:
      if(origObj)
        if(origObj->type!=cObjectMap) {
          ExecutiveDelete(oname);
          origObj=NULL;
        }
      PBlock(); /*PBlockAndUnlockAPI();*/
      obj=(Object*)ObjectMapLoadChemPyBrick((ObjectMap*)origObj,model,frame,discrete);
      PUnblock(); /*PLockAPIAndUnblock();*/
      if(!origObj) {
        if(obj) {
          ObjectSetName(obj,oname);
          ExecutiveManageObject(obj);
          sprintf(buf," CmdLoad: chempy.brick loaded into object \"%s\"\n",
                  oname);		  
        }
      } else if(origObj) {
        sprintf(buf," CmdLoad: chempy.brick appended into object \"%s\"\n",
                oname);
      }
      break;
    case cLoadTypeChemPyMap:
      if(origObj)
        if(origObj->type!=cObjectMap) {
          ExecutiveDelete(oname);
          origObj=NULL;
        }
      PBlock(); /*PBlockAndUnlockAPI();*/
      obj=(Object*)ObjectMapLoadChemPyMap((ObjectMap*)origObj,model,frame,discrete);
      PUnblock(); /*PLockAPIAndUnblock();*/
      if(!origObj) {
        if(obj) {
          ObjectSetName(obj,oname);
          ExecutiveManageObject(obj);
          sprintf(buf," CmdLoad: chempy.map loaded into object \"%s\"\n",
                  oname);		  
        }
      } else if(origObj) {
        sprintf(buf," CmdLoad: chempy.map appended into object \"%s\"\n",
                oname);
      }
      break;
    case cLoadTypeCallback:
      if(origObj)
        if(origObj->type!=cObjectCallback) {
          ExecutiveDelete(oname);
          origObj=NULL;
        }
      PBlock(); /*PBlockAndUnlockAPI();*/
      obj=(Object*)ObjectCallbackDefine((ObjectCallback*)origObj,model,frame);
      PUnblock(); /*PLockAPIAndUnblock();*/
      if(!origObj) {
        if(obj) {
          ObjectSetName(obj,oname);
          ExecutiveManageObject(obj);
          sprintf(buf," CmdLoad: pymol.callback loaded into object \"%s\"\n",
                  oname);		  
        }
      } else if(origObj) {
        sprintf(buf," CmdLoad: pymol.callback appended into object \"%s\"\n",
                oname);
      }
      break;
    case cLoadTypeCGO:
      if(origObj)
        if(origObj->type!=cObjectCGO) {
          ExecutiveDelete(oname);
          origObj=NULL;
        }
      PBlock(); /*PBlockAndUnlockAPI();*/
      obj=(Object*)ObjectCGODefine((ObjectCGO*)origObj,model,frame);
      PUnblock(); /*PLockAPIAndUnblock();*/
      if(!origObj) {
        if(obj) {
          ObjectSetName(obj,oname);
          ExecutiveManageObject(obj);
          sprintf(buf," CmdLoad: CGO loaded into object \"%s\"\n",
                  oname);		  
        }
      } else if(origObj) {
        sprintf(buf," CmdLoad: CGO appended into object \"%s\"\n",
                oname);
      }
      break;
      
    }
    if(origObj) {
      PRINTFB(FB_Executive,FB_Actions) 
        "%s",buf
        ENDFB;
      OrthoRestorePrompt();
    }
    APIExit();
  }
  return(APIStatus(ok));
}

static PyObject *CmdLoadCoords(PyObject *self, PyObject *args)
{
  char *oname;
  PyObject *model;
  Object *origObj = NULL,*obj;
  OrthoLineType buf;
  int frame,type;
  int ok=false;

  buf[0]=0;

  ok = PyArg_ParseTuple(args,"sOii",&oname,&model,&frame,&type);

  if (ok) {
    APIEntry();
    origObj=ExecutiveFindObjectByName(oname);
    
    /* TODO check for existing object of wrong type */
    if(!origObj) {
      ErrMessage("LoadCoords","named object not found.");
      ok=false;
    } else 
      {
        switch(type) {
        case cLoadTypeChemPyModel:
          PBlock(); /*PBlockAndUnlockAPI();*/
          obj=(Object*)ObjectMoleculeLoadCoords((ObjectMolecule*)origObj,model,frame);
          PUnblock(); /*PLockAPIAndUnblock();*/
          if(frame<0)
            frame=((ObjectMolecule*)obj)->NCSet-1;
          sprintf(buf," CmdLoad: Coordinates appended into object \"%s\", state %d.\n",
                  oname,frame+1);
          break;
        }
      }
    if(origObj) {
      PRINTFB(FB_Executive,FB_Actions) 
        "%s",buf
        ENDFB;
      OrthoRestorePrompt();
    }
    APIExit();
  }
  return(APIStatus(ok));
}

static PyObject *CmdLoad(PyObject *self, PyObject *args)
{
  char *fname,*oname;
  Object *origObj = NULL,*obj;
  OrthoLineType buf;
  int frame,type;
  int finish,discrete;
  int new_type;
  int ok=false;
  ok = PyArg_ParseTuple(args,"ssiiii",&oname,&fname,&frame,&type,&finish,&discrete);

  buf[0]=0;
  PRINTFD(FB_CCmd)
    "CmdLoad-DEBUG %s %s %d %d %d %d\n",
    oname,fname,frame,type,finish,discrete
    ENDFD;
  if (ok) {
    APIEntry();
    origObj=ExecutiveFindObjectByName(oname);
    /* check for existing object of right type, delete if not */
    if(origObj) {
      new_type = -1;
      switch(type) {
      case cLoadTypeChemPyModel:
      case cLoadTypePDB:
      case cLoadTypeXYZ:
      case cLoadTypePDBStr:
      case cLoadTypeMOL:
      case cLoadTypeMOLStr:
      case cLoadTypeMMD:
      case cLoadTypeMMDSeparate:
      case cLoadTypeMMDStr:
      case cLoadTypePMO:
        new_type = cObjectMolecule;
        break;
      case cLoadTypeChemPyBrick:
      case cLoadTypeChemPyMap:
      case cLoadTypeXPLORMap:
      case cLoadTypeCCP4Map:
        new_type = cObjectMap;
        break;
      case cLoadTypeCallback:
        new_type = cObjectCallback;
        break;
      case cLoadTypeCGO:
        new_type = cObjectCGO;
        break;
      }
      if (new_type!=origObj->type) {
        ExecutiveDelete(origObj->Name);
        origObj=NULL;
      }
    }
    
    
    switch(type) {
    case cLoadTypePDB:
      PRINTFD(FB_CCmd) " CmdLoad-DEBUG: loading PDB\n" ENDFD;
      if(!origObj) {
        obj=(Object*)ObjectMoleculeLoadPDBFile(NULL,fname,frame,discrete);
        if(obj) {
          ObjectSetName(obj,oname);
          ExecutiveManageObject(obj);
          if(frame<0)
            frame = ((ObjectMolecule*)obj)->NCSet-1;
          sprintf(buf," CmdLoad: \"%s\" loaded into object \"%s\", state %d.\n",
                  fname,oname,frame+1);
        }
      } else {
        ObjectMoleculeLoadPDBFile((ObjectMolecule*)origObj,fname,frame,discrete);
        if(finish)
          ExecutiveUpdateObjectSelection(origObj);
        if(frame<0)
          frame = ((ObjectMolecule*)origObj)->NCSet-1;
        sprintf(buf," CmdLoad: \"%s\" appended into object \"%s\", state %d.\n",
                fname,oname,frame+1);
      }
      break;
    case cLoadTypePMO:
      PRINTFD(FB_CCmd) " CmdLoad-DEBUG: loading PMO\n" ENDFD;
      if(!origObj) {
        obj=(Object*)ObjectMoleculeLoadPMOFile(NULL,fname,frame,discrete);
        if(obj) {
          ObjectSetName(obj,oname);
          ExecutiveManageObject(obj);
          if(frame<0)
            frame = ((ObjectMolecule*)obj)->NCSet-1;
          sprintf(buf," CmdLoad: \"%s\" loaded into object \"%s\", state %d.\n",
                  fname,oname,frame+1);
        }
      } else {
        ObjectMoleculeLoadPMOFile((ObjectMolecule*)origObj,fname,frame,discrete);
        if(finish)
          ExecutiveUpdateObjectSelection(origObj);
        if(frame<0)
          frame = ((ObjectMolecule*)origObj)->NCSet-1;
        sprintf(buf," CmdLoad: \"%s\" appended into object \"%s\", state %d.\n",
                fname,oname,frame+1);
      }
      break;
    case cLoadTypeXYZ:
      PRINTFD(FB_CCmd) " CmdLoad-DEBUG: loading XYZStr\n" ENDFD;
      if(!origObj) {
        obj=(Object*)ObjectMoleculeLoadXYZFile(NULL,fname,frame,discrete);
        if(obj) {
          ObjectSetName(obj,oname);
          ExecutiveManageObject(obj);
          if(frame<0)
            frame = ((ObjectMolecule*)obj)->NCSet-1;
          sprintf(buf," CmdLoad: \"%s\" loaded into object \"%s\", state %d.\n",
                  fname,oname,frame+1);
        }
      } else {
        ObjectMoleculeLoadXYZFile((ObjectMolecule*)origObj,fname,frame,discrete);
        if(finish)
          ExecutiveUpdateObjectSelection(origObj);
        if(frame<0)
          frame = ((ObjectMolecule*)origObj)->NCSet-1;
        sprintf(buf," CmdLoad: \"%s\" appended into object \"%s\", state %d.\n",
                fname,oname,frame+1);
      }
      break;
    case cLoadTypePDBStr:
      PRINTFD(FB_CCmd) " CmdLoad-DEBUG: loading PDBStr\n" ENDFD;
      obj=(Object*)ObjectMoleculeReadPDBStr((ObjectMolecule*)origObj,fname,frame,discrete);
      if(!origObj) {
        if(obj) {
          ObjectSetName(obj,oname);
          ExecutiveManageObject(obj);
          if(frame<0)
            frame = ((ObjectMolecule*)obj)->NCSet-1;
          sprintf(buf," CmdLoad: PDB-string loaded into object \"%s\", state %d.\n",
                  oname,frame+1);		  
        }
      } else if(origObj) {
        if(finish)
          ExecutiveUpdateObjectSelection(origObj);
        if(frame<0)
          frame = ((ObjectMolecule*)origObj)->NCSet-1;
        sprintf(buf," CmdLoad: PDB-string appended into object \"%s\", state %d.\n",
                oname,frame+1);
      }
      break;
    case cLoadTypeMOL:
      PRINTFD(FB_CCmd) " CmdLoad-DEBUG: loading MOL\n" ENDFD;
      obj=(Object*)ObjectMoleculeLoadMOLFile((ObjectMolecule*)origObj,fname,frame,discrete);
      if(!origObj) {
        if(obj) {
          ObjectSetName(obj,oname);
          ExecutiveManageObject(obj);
          if(frame<0)
            frame = ((ObjectMolecule*)obj)->NCSet-1;
          sprintf(buf," CmdLoad: \"%s\" loaded into object \"%s\", state %d.\n",
                  fname,oname,frame+1);		  
        }
      } else if(origObj) {
        if(finish)
          ExecutiveUpdateObjectSelection(origObj);
        if(frame<0)
          frame = ((ObjectMolecule*)origObj)->NCSet-1;
        sprintf(buf," CmdLoad: \"%s\" appended into object \"%s\", state %d.\n",
                fname,oname,frame+1);
      }
      break;
    case cLoadTypeMOLStr:
      PRINTFD(FB_CCmd) " CmdLoad-DEBUG: reading MOLStr\n" ENDFD;
      obj=(Object*)ObjectMoleculeReadMOLStr((ObjectMolecule*)origObj,fname,frame,discrete);
      if(!origObj) {
        if(obj) {
          ObjectSetName(obj,oname);
          ExecutiveManageObject(obj);
          if(frame<0)
            frame = ((ObjectMolecule*)obj)->NCSet-1;
          sprintf(buf," CmdLoad: MOL-string loaded into object \"%s\", state %d.\n",
                  oname,frame+1);		  
        }
      } else if(origObj) {
        if(finish)
          ExecutiveUpdateObjectSelection(origObj);
        if(frame<0)
          frame = ((ObjectMolecule*)origObj)->NCSet-1;
        sprintf(buf," CmdLoad: MOL-string appended into object \"%s\", state %d.\n",
                oname,frame+1);
      }
      break;
    case cLoadTypeMMD:
      PRINTFD(FB_CCmd) " CmdLoad-DEBUG: loading MMD\n" ENDFD;
      obj=(Object*)ObjectMoleculeLoadMMDFile((ObjectMolecule*)origObj,fname,frame,NULL,discrete);
      if(!origObj) {
        if(obj) {
          ObjectSetName(obj,oname);
          ExecutiveManageObject(obj);
          if(frame<0)
            frame = ((ObjectMolecule*)obj)->NCSet-1;
          sprintf(buf," CmdLoad: \"%s\" loaded into object \"%s\", state %d.\n",
                  fname,oname,frame+1);		  
        }
      } else if(origObj) {
        if(finish)
          ExecutiveUpdateObjectSelection(origObj);
        if(frame<0)
          frame = ((ObjectMolecule*)origObj)->NCSet-1;
        sprintf(buf," CmdLoad: \"%s\" appended into object \"%s\", state %d.\n",
                fname,oname,frame+1);
      }
      break;
    case cLoadTypeMMDSeparate:
      ObjectMoleculeLoadMMDFile((ObjectMolecule*)origObj,fname,frame,oname,discrete);
      break;
    case cLoadTypeMMDStr:
      PRINTFD(FB_CCmd) " CmdLoad-DEBUG: loading MMDStr\n" ENDFD;
      obj=(Object*)ObjectMoleculeReadMMDStr((ObjectMolecule*)origObj,fname,frame,discrete);
      if(!origObj) {
        if(obj) {
          ObjectSetName(obj,oname);
          ExecutiveManageObject(obj);
          if(frame<0)
            frame = ((ObjectMolecule*)obj)->NCSet-1;
          sprintf(buf," CmdLoad: MMD-string loaded into object \"%s\", state %d.\n",
                  oname,frame+1);		  
        }
      } else if(origObj) {
        if(finish)
          ExecutiveUpdateObjectSelection(origObj);
        if(frame<0)
          frame = ((ObjectMolecule*)origObj)->NCSet-1;
        sprintf(buf," CmdLoad: MMD-string appended into object \"%s\", state %d\n",
                oname,frame+1);
      }
      break;
    case cLoadTypeXPLORMap:
      PRINTFD(FB_CCmd) " CmdLoad-DEBUG: loading XPLORMap\n" ENDFD;
      if(!origObj) {
        obj=(Object*)ObjectMapLoadXPLORFile(NULL,fname,frame);
        if(obj) {
          ObjectSetName(obj,oname);
          ExecutiveManageObject((Object*)obj);
          sprintf(buf," CmdLoad: \"%s\" loaded into object \"%s\".\n",fname,oname);
        }
      } else {
        ObjectMapLoadXPLORFile((ObjectMap*)origObj,fname,frame);
        sprintf(buf," CmdLoad: \"%s\" appended into object \"%s\".\n",
                fname,oname);
      }
      break;
    case cLoadTypeCCP4Map:
      PRINTFD(FB_CCmd) " CmdLoad-DEBUG: loading CCP4Map\n" ENDFD;
      if(!origObj) {
        obj=(Object*)ObjectMapLoadCCP4File(NULL,fname,frame);
        if(obj) {
          ObjectSetName(obj,oname);
          ExecutiveManageObject((Object*)obj);
          sprintf(buf," CmdLoad: \"%s\" loaded into object \"%s\".\n",fname,oname);
        }
      } else {
        ObjectMapLoadXPLORFile((ObjectMap*)origObj,fname,frame);
        sprintf(buf," CmdLoad: \"%s\" appended into object \"%s\".\n",
                fname,oname);
      }
      break;
    }
    if(origObj) {
      PRINTFB(FB_Executive,FB_Actions) 
        "%s",buf
        ENDFB;
      OrthoRestorePrompt();
    }
    APIExit();
  }
  return(APIStatus(ok));
}

static PyObject *CmdOrigin(PyObject *self, PyObject *args)
{
  char *str1,*obj;
  OrthoLineType s1;
  
  int ok=false;
  ok = PyArg_ParseTuple(args,"ss",&str1,&obj);
  if (ok) {
    APIEntry();
    SelectorGetTmp(str1,s1);
    ok = ExecutiveCenter(s1,1,obj); /* TODO STATUS */
    SelectorFreeTmp(s1);
    APIExit();
  }
  return(APIStatus(ok));
}

static PyObject *CmdSort(PyObject *self, PyObject *args)
{
  char *name;
  int ok=false;
  ok = PyArg_ParseTuple(args,"s",&name);
  if (ok) {
    APIEntry();
    ExecutiveSort(name); /* TODO STATUS */
    APIExit();
  }
  return(APIStatus(ok));
}

static PyObject *CmdSpheroid(PyObject *self, PyObject *args)
/* EXPERIMENTAL */
{
  char *name;
  int ok=false;
  ok = PyArg_ParseTuple(args,"s",&name);
  if (ok) {
    APIEntry();
    ExecutiveSpheroid(name); /* TODO STATUS */
    APIExit();
  }
  return(APIStatus(ok));
}

static PyObject *CmdTest(PyObject *self, PyObject *args)
{
  int ok=true;
  /*  ok = PyArg_ParseTuple(args,"i",&int1);*/
  return(APIStatus(ok));
}

static PyObject *CmdZoom(PyObject *self, PyObject *args)
{
  char *str1;
  OrthoLineType s1;
  float buffer;

  int ok=false;
  ok = PyArg_ParseTuple(args,"sf",&str1,&buffer);
  if (ok) {
    APIEntry();
    SelectorGetTmp(str1,s1);
    ExecutiveWindowZoom(s1,buffer); /* TODO STATUS */
    SelectorFreeTmp(s1);
    APIExit();
  }
  return(APIStatus(ok));
}

static PyObject *CmdHAdd(PyObject *self, PyObject *args)
{
  char *str1;
  OrthoLineType s1;

  int ok=false;
  ok = PyArg_ParseTuple(args,"s",&str1);
  if (ok) {
    APIEntry();
    SelectorGetTmp(str1,s1);
    ExecutiveAddHydrogens(s1); /* TODO STATUS */
    SelectorFreeTmp(s1);
    SelectorGetTmp(str1,s1);
    ExecutiveAddHydrogens(s1); /* TODO STATUS */
    SelectorFreeTmp(s1);
    SelectorGetTmp(str1,s1);
    ExecutiveAddHydrogens(s1); /* TODO STATUS */
    SelectorFreeTmp(s1);
    SelectorGetTmp(str1,s1);
    ExecutiveAddHydrogens(s1); /* TODO STATUS */
    SelectorFreeTmp(s1);
    APIExit();
  }
  return(APIStatus(ok));
}

static PyObject *CmdRemove(PyObject *self, PyObject *args)
{
  char *str1;
  OrthoLineType s1;

  int ok=false;
  ok = PyArg_ParseTuple(args,"s",&str1);
  if (ok) {
    APIEntry();
    SelectorGetTmp(str1,s1);
    ExecutiveRemoveAtoms(s1); /* TODO STATUS */
    SelectorFreeTmp(s1);
    APIExit();
  }
  return(APIStatus(ok));
}

static PyObject *CmdRemovePicked(PyObject *self, PyObject *args)
{
  int i1;
  int ok=false;
  ok = PyArg_ParseTuple(args,"i",&i1);
  if (ok) {
    APIEntry();
    EditorRemove(i1); /* TODO STATUS */
    APIExit();
  }
  return(APIStatus(ok));
}

static PyObject *CmdHFill(PyObject *self, PyObject *args)
{
  int ok = true;
  APIEntry();
  EditorHFill(); /* TODO STATUS */
  APIExit();
  return(APIStatus(ok));
}

static PyObject *CmdCycleValence(PyObject *self, PyObject *args)
{
  int ok = true;
  APIEntry();
  EditorCycleValence();  /* TODO STATUS */
  APIExit();
  return(APIStatus(ok));
}

static PyObject *CmdReplace(PyObject *self, 	PyObject *args)
{
  int i1,i2;
  char *str1;
  int ok=false;
  ok = PyArg_ParseTuple(args,"sii",&str1,&i1,&i2);
  if (ok) {
    APIEntry();
    EditorReplace(str1,i1,i2);  /* TODO STATUS */
    APIExit();
  }
  return(APIStatus(ok));  
}

static PyObject *CmdSetGeometry(PyObject *self, 	PyObject *args)
{
  int i1,i2;
  char *str1;
  OrthoLineType s1;
  int ok=false;
  ok = PyArg_ParseTuple(args,"sii",&str1,&i1,&i2);
  if (ok) {
    APIEntry();
    SelectorGetTmp(str1,s1);
    ok = ExecutiveSetGeometry(s1,i1,i2);  /* TODO STATUS */
    SelectorFreeTmp(s1);
    APIExit();
  }
  return(APIStatus(ok));  
}

static PyObject *CmdAttach(PyObject *self, 	PyObject *args)
{
  int i1,i2;
  char *str1;
  int ok=false;
  ok = PyArg_ParseTuple(args,"sii",&str1,&i1,&i2);
  if (ok) {
    APIEntry();
    EditorAttach(str1,i1,i2);  /* TODO STATUS */
    APIExit();
  }
  return(APIStatus(ok));  
}

static PyObject *CmdFuse(PyObject *self, 	PyObject *args)
{
  char *str1,*str2;
  int mode;
  OrthoLineType s1,s2;

  int ok=false;
  ok = PyArg_ParseTuple(args,"ssi",&str1,&str2,&mode);
  if (ok) {
    APIEntry();
    SelectorGetTmp(str1,s1);
    SelectorGetTmp(str2,s2);
    ExecutiveFuse(s1,s2,mode);  /* TODO STATUS */
    SelectorFreeTmp(s1);
    SelectorFreeTmp(s2);
    APIExit();
  }
  return(APIStatus(ok));  
}

static PyObject *CmdUnpick(PyObject *self, 	PyObject *args)
{
  APIEntry();
  EditorInactive();  /* TODO STATUS */
  APIExit();
  return(APISuccess());  
}

static PyObject *CmdEdit(PyObject *self, 	PyObject *args)
{
  char *str0,*str1,*str2,*str3;
  OrthoLineType s0 = "";
  OrthoLineType s1 = "";
  OrthoLineType s2 = "";
  OrthoLineType s3 = "";
  int pkresi;
  int ok=false;
  ok = PyArg_ParseTuple(args,"ssssi",&str0,&str1,&str2,&str3,&pkresi);
  if (ok) {
    APIEntry();
    if(str0[0]) SelectorGetTmp(str0,s0);
    if(str1[0]) SelectorGetTmp(str1,s1);
    if(str2[0]) SelectorGetTmp(str2,s2);
    if(str3[0]) SelectorGetTmp(str3,s3);
    ok = EditorSelect(s0,s1,s2,s3,pkresi);
    if(s0[0]) SelectorFreeTmp(s0);
    if(s1[0]) SelectorFreeTmp(s1);
    if(s2[0]) SelectorFreeTmp(s2);
    if(s3[0]) SelectorFreeTmp(s3);
    
    APIExit();
  }
  return APIStatus(ok);
}

static PyObject *CmdRename(PyObject *self, 	PyObject *args)
{
  char *str1;
  int int1;
  OrthoLineType s1;

  int ok=false;
  ok = PyArg_ParseTuple(args,"si",&str1,&int1);
  if (ok) {
    APIEntry();
    SelectorGetTmp(str1,s1);
    ExecutiveRenameObjectAtoms(s1,int1); /* TODO STATUS */
    SelectorFreeTmp(s1);
    APIExit();
  }
  return(APIStatus(ok));  
}


void init_cmd(void)
{
  PyImport_AddModule("_cmd");
  Py_InitModule("_cmd", Cmd_methods);
}

