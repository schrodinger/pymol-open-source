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

#define tmpSele "_tmp"
#define tmpSele1 "_tmp1"
#define tmpSele2 "_tmp2"

int flush_count = 0;

static void APIEntry(void)
{
  P_glut_thread_keep_out++;
  PUnblock();
}

static void APIExit(void)
{
  PBlock();
  P_glut_thread_keep_out--;
}

static PyObject *CmdAlter(PyObject *self,   PyObject *args);
static PyObject *CmdAlterState(PyObject *self,   PyObject *args);
static PyObject *CmdBond(PyObject *dummy, PyObject *args);
static PyObject *CmdButton(PyObject *dummy, PyObject *args);
static PyObject *CmdClip(PyObject *self, 	PyObject *args);
static PyObject *CmdCls(PyObject *self, 	PyObject *args);
static PyObject *CmdColor(PyObject *self, PyObject *args);
static PyObject *CmdColorDef(PyObject *self, 	PyObject *args);
static PyObject *CmdCopy(PyObject *self, PyObject *args);
static PyObject *CmdCountStates(PyObject *self, PyObject *args);
static PyObject *CmdCreate(PyObject *self, PyObject *args);
static PyObject *CmdDelete(PyObject *self, PyObject *args);
static PyObject *CmdDirty(PyObject *self, 	PyObject *args);
static PyObject *CmdDist(PyObject *dummy, PyObject *args);
static PyObject *CmdDistance(PyObject *dummy, PyObject *args);
static PyObject *CmdDo(PyObject *self, 	PyObject *args);
static PyObject *CmdDump(PyObject *self, 	PyObject *args);
static PyObject *CmdExportDots(PyObject *self, PyObject *args);
static PyObject *CmdFit(PyObject *dummy, PyObject *args);
static PyObject *CmdFitPairs(PyObject *dummy, PyObject *args);
static PyObject *CmdFlag(PyObject *self, 	PyObject *args);
static PyObject *CmdFlushNow(PyObject *self, 	PyObject *args);
static PyObject *CmdIdentify(PyObject *dummy, PyObject *args);
static PyObject *CmdIntraFit(PyObject *dummy, PyObject *args);
static PyObject *CmdIsomesh(PyObject *self, 	PyObject *args);
static PyObject *CmdFinishObject(PyObject *self, PyObject *args);
static PyObject *CmdFrame(PyObject *self, PyObject *args);
static PyObject *CmdGet(PyObject *self, 	PyObject *args);
static PyObject *CmdGetPDB(PyObject *dummy, PyObject *args);
static PyObject *CmdGetMatrix(PyObject *self, 	PyObject *args);
static PyObject *CmdGetMinMax(PyObject *self, 	PyObject *args);
static PyObject *CmdGetModel(PyObject *dummy, PyObject *args);
static PyObject *CmdGetFeedback(PyObject *dummy, PyObject *args);
static PyObject *CmdGetMoment(PyObject *self, 	PyObject *args);
static PyObject *CmdGetSetting(PyObject *self, 	PyObject *args);
static PyObject *CmdMask(PyObject *self, PyObject *args);
static PyObject *CmdMem(PyObject *self, 	PyObject *args);
static PyObject *CmdLabel(PyObject *self,   PyObject *args);
static PyObject *CmdLoad(PyObject *self, 	PyObject *args);
static PyObject *CmdLoadCoords(PyObject *self, PyObject *args);
static PyObject *CmdLoadObject(PyObject *self, PyObject *args);
static PyObject *CmdMClear(PyObject *self, 	PyObject *args);
static PyObject *CmdMDo(PyObject *self, 	PyObject *args);
static PyObject *CmdMMatrix(PyObject *self, 	PyObject *args);
static PyObject *CmdMove(PyObject *self, 	PyObject *args);
static PyObject *CmdMPlay(PyObject *self, 	PyObject *args);
static PyObject *CmdMPNG(PyObject *self, 	PyObject *args);
static PyObject *CmdMSet(PyObject *self, 	PyObject *args);
static PyObject *CmdOrigin(PyObject *self, PyObject *args);
static PyObject *CmdOnOff(PyObject *self, 	PyObject *args);
static PyObject *CmdOrient(PyObject *dummy, PyObject *args);
static PyObject *CmdOverlap(PyObject *self, 	PyObject *args);
static PyObject *CmdPaste(PyObject *self, 	PyObject *args);
static PyObject *CmdPNG(PyObject *self, 	PyObject *args);
static PyObject *CmdProtect(PyObject *self, PyObject *args);
static PyObject *CmdQuit(PyObject *self, 	PyObject *args);
static PyObject *CmdReset(PyObject *self, PyObject *args);
static PyObject *CmdRay(PyObject *self, 	PyObject *args);
static PyObject *CmdRemove(PyObject *self, PyObject *args);
static PyObject *CmdResetRate(PyObject *dummy, PyObject *args);
static PyObject *CmdRefresh(PyObject *self, 	PyObject *args);
static PyObject *CmdRefreshNow(PyObject *self, 	PyObject *args);
static PyObject *CmdReady(PyObject *dummy, PyObject *args);
static PyObject *CmdRock(PyObject *self, PyObject *args);
static PyObject *CmdRunPyMOL(PyObject *dummy, PyObject *args);
static PyObject *CmdSelect(PyObject *self, PyObject *args);
static PyObject *CmdSetMatrix(PyObject *self, 	PyObject *args);
static PyObject *CmdSet(PyObject *self, 	PyObject *args);
static PyObject *CmdSetFrame(PyObject *self, PyObject *args);
static PyObject *CmdShowHide(PyObject *self, 	PyObject *args);
static PyObject *CmdSort(PyObject *dummy, PyObject *args);
static PyObject *CmdSplash(PyObject *dummy, PyObject *args);
static PyObject *CmdStereo(PyObject *self, PyObject *args);
static PyObject *CmdSystem(PyObject *dummy, PyObject *args);
static PyObject *CmdSymExp(PyObject *dummy, PyObject *args);
static PyObject *CmdTest(PyObject *self, 	PyObject *args);
static PyObject *CmdTurn(PyObject *self, 	PyObject *args);
static PyObject *CmdViewport(PyObject *self, 	PyObject *args);
static PyObject *CmdZoom(PyObject *self, PyObject *args);
static PyObject *CmdWaitQueue(PyObject *self, 	PyObject *args);

static PyMethodDef Cmd_methods[] = {
	{"alter",	     CmdAlter,        METH_VARARGS },
	{"alter_state",  CmdAlterState,   METH_VARARGS },
	{"bond",         CmdBond,         METH_VARARGS },
   {"button",       CmdButton,       METH_VARARGS },
	{"clip",	        CmdClip,         METH_VARARGS },
	{"cls",	        CmdCls,          METH_VARARGS },
	{"color",	     CmdColor,        METH_VARARGS },
	{"colordef",	  CmdColorDef,     METH_VARARGS },
	{"copy",         CmdCopy,         METH_VARARGS },
	{"create",       CmdCreate,       METH_VARARGS },
	{"count_states", CmdCountStates,  METH_VARARGS },
	{"delete",       CmdDelete,       METH_VARARGS },
	{"dirty",        CmdDirty,        METH_VARARGS },
	{"distance",	  CmdDistance,     METH_VARARGS },
	{"dist",    	  CmdDist,         METH_VARARGS },
	{"do",	        CmdDo,           METH_VARARGS },
	{"dump",	        CmdDump,         METH_VARARGS },
	{"export_dots",  CmdExportDots,   METH_VARARGS },
	{"finish_object",CmdFinishObject, METH_VARARGS },
	{"fit",          CmdFit,          METH_VARARGS },
	{"fit_pairs",    CmdFitPairs,     METH_VARARGS },
	{"flag",         CmdFlag,         METH_VARARGS },
	{"frame",	     CmdFrame,        METH_VARARGS },
   {"flush_now",    CmdFlushNow,     METH_VARARGS },
	{"get",	        CmdGet,          METH_VARARGS },
	{"get_feedback", CmdGetFeedback,  METH_VARARGS },
	{"get_matrix",	  CmdGetMatrix,    METH_VARARGS },
	{"get_min_max",  CmdGetMinMax,    METH_VARARGS },
	{"get_model",	  CmdGetModel,     METH_VARARGS },
	{"get_moment",	  CmdGetMoment,    METH_VARARGS },
	{"get_pdb",	     CmdGetPDB,       METH_VARARGS },
	{"get_setting",  CmdGetSetting,   METH_VARARGS },
	{"identify",     CmdIdentify,     METH_VARARGS },
	{"intrafit",     CmdIntraFit,     METH_VARARGS },
	{"isomesh",	     CmdIsomesh,      METH_VARARGS },
   {"wait_queue",   CmdWaitQueue,    METH_VARARGS },
   {"label",        CmdLabel,        METH_VARARGS },
	{"load",	        CmdLoad,         METH_VARARGS },
	{"load_coords",  CmdLoadCoords,   METH_VARARGS },
	{"load_object",  CmdLoadObject,   METH_VARARGS },
	{"mask",	        CmdMask,         METH_VARARGS },
	{"mclear",	     CmdMClear,       METH_VARARGS },
	{"mdo",	        CmdMDo,          METH_VARARGS },
	{"mem",	        CmdMem,          METH_VARARGS },
	{"move",	        CmdMove,         METH_VARARGS },
	{"mset",	        CmdMSet,         METH_VARARGS },
	{"mplay",	     CmdMPlay,        METH_VARARGS },
	{"mpng_",	     CmdMPNG,         METH_VARARGS },
	{"mmatrix",	     CmdMMatrix,      METH_VARARGS },
	{"origin",	     CmdOrigin,       METH_VARARGS },
	{"orient",	     CmdOrient,       METH_VARARGS },
	{"onoff",        CmdOnOff,        METH_VARARGS },
	{"overlap",      CmdOverlap,      METH_VARARGS },
	{"paste",	     CmdPaste,        METH_VARARGS },
	{"png",	        CmdPNG,          METH_VARARGS },
	{"protect",	     CmdProtect,      METH_VARARGS },
	{"quit",	        CmdQuit,         METH_VARARGS },
	{"ready",        CmdReady,        METH_VARARGS },
	{"refresh",      CmdRefresh,      METH_VARARGS },
	{"refresh_now",  CmdRefreshNow,   METH_VARARGS },
	{"remove",	     CmdRemove,       METH_VARARGS },
	{"render",	     CmdRay,          METH_VARARGS },
	{"reset",        CmdReset,        METH_VARARGS },
	{"reset_rate",	  CmdResetRate,    METH_VARARGS },
	{"rock",	        CmdRock,         METH_VARARGS },
	{"runpymol",	  CmdRunPyMOL,     METH_VARARGS },
	{"select",       CmdSelect,       METH_VARARGS },
	{"set",	        CmdSet,          METH_VARARGS },
	{"setframe",	  CmdSetFrame,     METH_VARARGS },
	{"showhide",     CmdShowHide,     METH_VARARGS },
	{"set_matrix",	  CmdSetMatrix,    METH_VARARGS },
	{"sort",         CmdSort,         METH_VARARGS },
	{"splash",       CmdSplash,       METH_VARARGS },
	{"stereo",	     CmdStereo,       METH_VARARGS },
	{"system",	     CmdSystem,       METH_VARARGS },
	{"symexp",	     CmdSymExp,       METH_VARARGS },
	{"test",	        CmdTest,         METH_VARARGS },
	{"turn",	        CmdTurn,         METH_VARARGS },
	{"viewport",     CmdViewport,     METH_VARARGS },
	{"zoom",	        CmdZoom,         METH_VARARGS },
	{NULL,		     NULL}		/* sentinel */
};

static PyObject *CmdMask(PyObject *self, PyObject *args)
{
  char *str1;
  int int1;
  OrthoLineType s1;

  PyArg_ParseTuple(args,"si",&str1,&int1);
  APIEntry();
  SelectorGetTmp(str1,s1);
  ExecutiveMask(s1,int1);
  SelectorFreeTmp(s1);
  APIExit();
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject *CmdProtect(PyObject *self, PyObject *args)
{
  char *str1;
  int int1;
  OrthoLineType s1;

  PyArg_ParseTuple(args,"si",&str1,&int1);
  APIEntry();
  SelectorGetTmp(str1,s1);
  ExecutiveProtect(s1,int1);
  SelectorFreeTmp(s1);
  APIExit();
  Py_INCREF(Py_None);
  return Py_None;
}


static PyObject *CmdButton(PyObject *self, 	PyObject *args)
{
  int i1,i2;
  PyArg_ParseTuple(args,"ii",&i1,&i2);
  APIEntry();
  ButModeSet(i1,i2);
  APIExit();
  Py_INCREF(Py_None);
  return Py_None;  
}

static PyObject *CmdFlushNow(PyObject *self, 	PyObject *args)
{
  /* only called by the GLUT thread with unlocked API */
  P_glut_thread_keep_out++;
  /*  if(!flush_count) {
      flush_count++;*/
  PFlushFast();
    /*    flush_count--;
          }*/
  P_glut_thread_keep_out--;
  Py_INCREF(Py_None);
  return Py_None;  
}

static PyObject *CmdWaitQueue(PyObject *self, 	PyObject *args)
{
  PyObject *result;
  P_glut_thread_keep_out++;
  if(OrthoCommandWaiting()) 
    result = PyInt_FromLong(1);
  else
    result = PyInt_FromLong(0);
  P_glut_thread_keep_out--;
  return result;
}

static PyObject *CmdPaste(PyObject *dummy, PyObject *args)
{
  PyObject *list,*str;
  char *st;
  int l,a;
  APIEntry();
  PyArg_ParseTuple(args,"O",&list);
  if(list) 
    if(PyList_Check(list)) 
      {
        l=PyList_Size(list);
        for(a=0;a<l;a++) {
          str = PyList_GetItem(list,a);
          if(str)
            if(PyString_Check(str)) {
              st = PyString_AsString(str);
              OrthoPasteIn(st);
            }
        }
      }
  APIExit();
  Py_INCREF(Py_None);
  return Py_None;  
}

static PyObject *CmdSplash(PyObject *dummy, PyObject *args)
{
  APIEntry();
  OrthoSplash();
  APIExit();
  Py_INCREF(Py_None);
  return Py_None;  
}

static PyObject *CmdCls(PyObject *dummy, PyObject *args)
{
  APIEntry();
  OrthoClear();
  APIExit();
  Py_INCREF(Py_None);
  return Py_None;  
}

static PyObject *CmdDump(PyObject *dummy, PyObject *args)
{
  char *str1,*str2;
  PyArg_ParseTuple(args,"ss",&str1,&str2);
  APIEntry();
  ExecutiveDump(str1,str2);
  APIExit();
  Py_INCREF(Py_None);
  return Py_None;  
}

static PyObject *CmdIsomesh(PyObject *self, 	PyObject *args) {
  char *str1,*str2,*str3,*str4;
  float lvl,fbuf;
  int dotFlag;
  int c;
  OrthoLineType s1;
  int oper,frame;
  Object *obj,*mObj;
  ObjectMap *mapObj;
  float mn[3] = { 0,0,0};
  float mx[3] = { 15,15,15};
  OrthoLineType buf;

  /* oper 0 = all, 1 = sele + buffer, 2 = vector */

  PyArg_ParseTuple(args,"sisissfi",&str1,&frame,&str2,&oper,&str3,&str4,&lvl,&dotFlag);
  APIEntry();
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
      break;
    case 1:
      SelectorGetTmp(str3,s1);
      ExecutiveGetExtent(s1,mn,mx);
      SelectorFreeTmp(s1);
      if(sscanf(str4,"%f",&fbuf)==1) {
        for(c=0;c<3;c++) {
          mn[c]-=fbuf;
          mx[c]+=fbuf;
        }
      }
      break;
    }
    obj=ExecutiveFindObjectByName(str1);  
    if(obj) {
      ExecutiveDelete(obj->Name);
      obj=NULL;
    }
    obj=(Object*)ObjectMeshFromBox(mapObj,mn,mx,lvl,dotFlag);
    if(obj) {
      ObjectSetName(obj,str1);
      ExecutiveManageObject((Object*)obj);
      sprintf(buf," Mesh: created \"%s\", setting level to %5.3f\n",str1,lvl);
      OrthoAddOutput(buf);
    }
  } else {
    sprintf(buf,"Map object '%s' not found.",str2);
    ErrMessage("Mesh",buf);
  }
  APIExit();
  Py_INCREF(Py_None);
  return Py_None;  
}

static PyObject *CmdSymExp(PyObject *self, 	PyObject *args) {
  char *str1,*str2,*str3;
  OrthoLineType s1;
  float cutoff;
  Object *mObj;
  /* oper 0 = all, 1 = sele + buffer, 2 = vector */

  PyArg_ParseTuple(args,"sssf",&str1,&str2,&str3,&cutoff);
  APIEntry();
  mObj=ExecutiveFindObjectByName(str2);  
  if(mObj) {
    if(mObj->type!=cObjectMolecule) {
      mObj=NULL;
    }
  }
  if(mObj) {
    SelectorGetTmp(str3,s1);
    ExecutiveSymExp(str1,str2,s1,cutoff);
    SelectorFreeTmp(s1);
  }
  APIExit();
  Py_INCREF(Py_None);
  return Py_None;  
}



static PyObject *CmdOverlap(PyObject *dummy, PyObject *args)
{
  char *str1,*str2;
  int state1,state2;
  float overlap;
  OrthoLineType s1,s2;
  PyObject *result;
  PyArg_ParseTuple(args,"ssii",&str1,&str2,&state1,&state2);
  APIEntry();
  SelectorGetTmp(str1,s1);
  SelectorGetTmp(str2,s2);
  overlap = ExecutiveOverlap(s1,state1,s2,state2);
  SelectorFreeTmp(s1);
  SelectorFreeTmp(s2);
  APIExit();
  result = Py_BuildValue("f",overlap);
  return result;
}

static PyObject *CmdDist(PyObject *dummy, PyObject *args)
{
  char *name,*str1,*str2;
  float cutoff;
  int mode;
  OrthoLineType s1,s2;
  PyArg_ParseTuple(args,"sssif",&name,&str1,&str2,&mode,&cutoff);
  APIEntry();
  SelectorGetTmp(str1,s1);
  SelectorGetTmp(str2,s2);
  ExecutiveDist(name,s1,s2,mode,cutoff);
  SelectorFreeTmp(s1);
  SelectorFreeTmp(s2);
  APIExit();
  Py_INCREF(Py_None);
  return Py_None;  
}

static PyObject *CmdBond(PyObject *dummy, PyObject *args)
{
  char *str1,*str2;
  int order,mode;
  OrthoLineType s1,s2;
  PyArg_ParseTuple(args,"ssii",&str1,&str2,&order,&mode);
  APIEntry();
  SelectorGetTmp(str1,s1);
  SelectorGetTmp(str2,s2);
  ExecutiveBond(s1,s2,order,mode);
  SelectorFreeTmp(s1);
  SelectorFreeTmp(s2);
  APIExit();
  Py_INCREF(Py_None);
  return Py_None;  
}

static PyObject *CmdDistance(PyObject *dummy, PyObject *args)
{
  char *str1,*str2;
  OrthoLineType s1,s2;
  float dist;
  PyObject *result;
  PyArg_ParseTuple(args,"ss",&str1,&str2);
  APIEntry();
  SelectorGetTmp(str1,s1);
  SelectorGetTmp(str2,s2);
  dist = ExecutiveDistance(s1,s2);
  SelectorFreeTmp(s1);
  SelectorFreeTmp(s2);
  APIExit();
  result = Py_BuildValue("f",dist);
  return result;
}

static PyObject *CmdLabel(PyObject *self,   PyObject *args)
{
  char *str1,*str2;
  OrthoLineType s1;

  PyArg_ParseTuple(args,"ss",&str1,&str2);
  APIEntry();
  SelectorGetTmp(str1,s1);
  ExecutiveLabel(s1,str2);
  SelectorFreeTmp(s1);
  APIExit();
  Py_INCREF(Py_None);
  return Py_None;

}

static PyObject *CmdAlter(PyObject *self,   PyObject *args)
{
  char *str1,*str2;
  OrthoLineType s1;

  PyArg_ParseTuple(args,"ss",&str1,&str2);
  APIEntry();
  SelectorGetTmp(str1,s1);
  ExecutiveAlter(s1,str2);
  SelectorFreeTmp(s1);
  APIExit();
  Py_INCREF(Py_None);
  return Py_None;

}

static PyObject *CmdAlterState(PyObject *self,   PyObject *args)
{
  char *str1,*str2;
  int i1;
  OrthoLineType s1;

  PyArg_ParseTuple(args,"iss",&i1,&str1,&str2);
  APIEntry();
  SelectorGetTmp(str1,s1);
  ExecutiveAlterState(i1,s1,str2);
  SelectorFreeTmp(s1);
  APIExit();
  Py_INCREF(Py_None);
  return Py_None;

}

static PyObject *CmdCopy(PyObject *self,   PyObject *args)
{
  char *str1,*str2;
  PyArg_ParseTuple(args,"ss",&str1,&str2);
  APIEntry();
  ExecutiveCopy(str1,str2);
  APIExit();
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject *CmdResetRate(PyObject *dummy, PyObject *args)
{
  PyObject *result = NULL;
  result=Py_None;
  APIEntry();
  ButModeResetRate();
  APIExit();
  Py_INCREF(result);
  return(result);
}

static PyObject *CmdReady(PyObject *dummy, PyObject *args)
{
  PyObject *result = NULL;
  result = Py_BuildValue("i",PyMOLReady);
  return(result);
}

static PyObject *CmdMem(PyObject *dummy, PyObject *args)
{
  PyObject *result = NULL;
  MemoryDebugDump();
  result = Py_None;
  Py_INCREF(result);
  return(result);
}

static PyObject *CmdRunPyMOL(PyObject *dummy, PyObject *args)
{
#ifdef _PYMOL_MODULE
  was_main();
#endif
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject *CmdCountStates(PyObject *dummy, PyObject *args)
{
  char *str1;
  OrthoLineType s1;
  PyObject *result;
  int states;

  PyArg_ParseTuple(args,"s",&str1);
  APIEntry();
  SelectorGetTmp(str1,s1);
  states = ExecutiveCountStates(s1);
  SelectorFreeTmp(s1); 
  APIExit();
  result = Py_BuildValue("i",states);
  return(result);
}

static PyObject *CmdIdentify(PyObject *dummy, PyObject *args)
{
  char *str1;
  OrthoLineType s1;
  int mode;
  PyObject *result = Py_None;
  int *iVLA=NULL;
  PyArg_ParseTuple(args,"si",&str1,&mode);
  APIEntry();
  SelectorGetTmp(str1,s1);
  iVLA=ExecutiveIdentify(s1,mode);
  SelectorFreeTmp(s1);
  APIExit();
  if(iVLA) {
    result=PConvIntVLAToPyList(iVLA);
    VLAFreeP(iVLA);
  } else {
    result = PyList_New(0);
  }
  if(result==Py_None) Py_INCREF(result);
  return(result);
}

static PyObject *CmdSystem(PyObject *dummy, PyObject *args)
{
  char *str1;
  PyObject *result = NULL;
  int code;

  PyArg_ParseTuple(args,"s",&str1);
  APIEntry();
  code = system(str1);
  APIExit();
  result = Py_BuildValue("i",code);
  return(result);
}

static PyObject *CmdGetFeedback(PyObject *dummy, PyObject *args)
{
  OrthoLineType buffer;
  PyObject *result = NULL;
  int code;

  code = OrthoFeedbackOut(buffer);
  if(code)
    result = Py_BuildValue("s",buffer);
  if(!result) {
	result=Py_None;
	Py_INCREF(result);
  }
  return(result);
}


static PyObject *CmdGetPDB(PyObject *dummy, PyObject *args)
{
  char *str1;
  char *pdb = NULL;
  int state;
  OrthoLineType s1;

  PyObject *result = NULL;
  
  PyArg_ParseTuple(args,"si",&str1,&state);
  APIEntry();
  SelectorGetTmp(str1,s1);
  pdb=ExecutiveSeleToPDBStr(s1,state,true);
  SelectorFreeTmp(s1);
  APIExit();
  if(pdb)
    result = Py_BuildValue("s",pdb);
  if(!result) {
    result=Py_None;
    Py_INCREF(result);
  }
  FreeP(pdb);
  return(result);
}

static PyObject *CmdGetModel(PyObject *dummy, PyObject *args)
{
  char *str1;
  int state;
  OrthoLineType s1;
  PyObject *model;
  PyObject *result = NULL;
  
  PyArg_ParseTuple(args,"si",&str1,&state);
  APIEntry();
  SelectorGetTmp(str1,s1);
  model=ExecutiveSeleToChemPyModel(s1,state);
  SelectorFreeTmp(s1);
  APIExit();
  if(model) {
    result = model;
  }
  if(!result) {
    result=Py_None;
    Py_INCREF(result);
  }
  return(model);
}

static PyObject *CmdCreate(PyObject *dummy, PyObject *args)
{
  char *str1,*str2;
  int target,source;
  OrthoLineType s1;

  PyObject *result = NULL;
  PyArg_ParseTuple(args,"ssii",&str1,&str2,&source,&target);
  APIEntry();
  SelectorGetTmp(str2,s1);
  ExecutiveSeleToObject(str1,s1,source,target);
  SelectorFreeTmp(s1);
  APIExit();
  result=Py_None;
  Py_INCREF(result);
  return(result);
}


static PyObject *CmdOrient(PyObject *dummy, PyObject *args)
{
  Matrix33d m;
  char *str1;
  OrthoLineType s1;

  PyArg_ParseTuple(args,"s",&str1);
  APIEntry();
  SelectorGetTmp(str1,s1);
  if(ExecutiveGetMoment(s1,m))
    ExecutiveOrient(s1,m);
  SelectorFreeTmp(s1); 
  APIExit();
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject *CmdFitPairs(PyObject *dummy, PyObject *args)
{
  PyObject *list;
  WordType *word = NULL;
  int ln=0;
  int a;
  int ok=true;
  PyObject *result = NULL;
  float valu;
  PyArg_ParseTuple(args,"O",&list);
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
  if(!result) {
    result=Py_None;
    Py_INCREF(result);
  }
  return result;
}

static PyObject *CmdIntraFit(PyObject *dummy, PyObject *args)
{
  char *str1;
  int state;
  int mode;
  OrthoLineType s1;
  float *fVLA;
  PyObject *result=Py_None;
  PyArg_ParseTuple(args,"sii",&str1,&state,&mode);
  if(state<0) state=0;
  APIEntry();
  SelectorGetTmp(str1,s1);
  fVLA=ExecutiveRMSStates(s1,state,mode);
  SelectorFreeTmp(s1);
  APIExit();
  if(fVLA) {
    result=PConvFloatVLAToPyList(fVLA);
    VLAFreeP(fVLA);
  }
  if(result==Py_None) Py_INCREF(result);
  return result;
}

static PyObject *CmdFit(PyObject *dummy, PyObject *args)
{
  char *str1,*str2;
  int mode;
  OrthoLineType s1,s2;
  PyObject *result;
  PyArg_ParseTuple(args,"ssi",&str1,&str2,&mode);
  APIEntry();
  SelectorGetTmp(str1,s1);
  SelectorGetTmp(str2,s2);
  result=Py_BuildValue("f",ExecutiveRMS(s1,s2,mode));
  SelectorFreeTmp(s1);
  SelectorFreeTmp(s2);
  APIExit();
  return result;
}

static PyObject *CmdDirty(PyObject *self, 	PyObject *args)
{
  APIEntry();
  OrthoDirty();
  APIExit();
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject *CmdDo(PyObject *self, 	PyObject *args)
{
  char *str1;

  PyArg_ParseTuple(args,"s",&str1);
  APIEntry();
  if(str1[0]!='_') {
    OrthoAddOutput("PyMOL>");
    OrthoAddOutput(str1);
  }
  PParse(str1);
  OrthoNewLine(NULL);
  APIExit();
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject *CmdRock(PyObject *self, PyObject *args)
{
  APIEntry();
  ControlRock(-1);
  APIExit();
  return Py_None;
}

static PyObject *CmdGetMoment(PyObject *self, 	PyObject *args)
{
  Matrix33d m;
  PyObject *result;

  char *str1;
  PyArg_ParseTuple(args,"s",&str1);
  APIEntry();
  ExecutiveGetMoment(str1,m);
  APIExit();
  result = Py_BuildValue("(ddd)(ddd)(ddd)", 
								 m[0][0],m[0][1],m[0][2],
								 m[1][0],m[1][1],m[1][2],
								 m[2][0],m[2][1],m[2][2]);
  return result;
}

static PyObject *CmdGetSetting(PyObject *self, 	PyObject *args)
{
  PyObject *result;
  char *str1;
  float value;
  PyArg_ParseTuple(args,"s",&str1);
  APIEntry();
  value=SettingGetNamed(str1);
  APIExit();
  result = Py_BuildValue("f", SettingGetNamed(str1));
  return result;
}

static PyObject *CmdExportDots(PyObject *self, 	PyObject *args)
{
  PyObject *result=NULL;
  PyObject *cObj;
  ExportDotsObj *obj;
  char *str1;
  int int1;
  
  PyArg_ParseTuple(args,"si",&str1,&int1);
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
  if(!result)
	 {
		result = Py_None;
		Py_INCREF(result);
	 }
  return result;
}

static PyObject *CmdSetFrame(PyObject *self, PyObject *args)
{
  int mode,frm;
  PyArg_ParseTuple(args,"ii",&mode,&frm);
  APIEntry();
  SceneSetFrame(mode,frm);
  APIExit();
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject *CmdFrame(PyObject *self, PyObject *args)
{
  int frm;
  PyArg_ParseTuple(args,"i",&frm);
  frm--;
  if(frm<0) frm=0;
  APIEntry();
  SceneSetFrame(0,frm);
  APIExit();
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject *CmdStereo(PyObject *self, PyObject *args)
{
  PyObject *result;
  int i1;
  if(StereoCapable) {
  PyArg_ParseTuple(args,"i",&i1);
  APIEntry();
  ExecutiveStereo(i1);
  APIExit();
  result=Py_BuildValue("i",1);
  } else {
  result=Py_BuildValue("i",0);
  }
  return result;
}

static PyObject *CmdReset(PyObject *self, PyObject *args)
{
  int cmd;
  PyArg_ParseTuple(args,"i",&cmd);
  APIEntry();
  ExecutiveReset(cmd);
  APIExit();
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject *CmdSetMatrix(PyObject *self, 	PyObject *args)
{
  float m[16];
  PyArg_ParseTuple(args,"ffffffffffffffff",
						 &m[0],&m[1],&m[2],&m[3],
						 &m[4],&m[5],&m[6],&m[7],
						 &m[8],&m[9],&m[10],&m[11],
						 &m[12],&m[13],&m[14],&m[15]);
  APIEntry();
  SceneSetMatrix(m);
  APIExit();
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject *CmdGetMinMax(PyObject *self, 	PyObject *args)
{
  float mn[3],mx[3];
  char *str1;
  int state;
  OrthoLineType s1;
  PyObject *result;
  int flag;

  PyArg_ParseTuple(args,"si",&str1,&state); /* state currently ignored */
  APIEntry();
  SelectorGetTmp(str1,s1);
  flag = ExecutiveGetExtent(s1,mn,mx);
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
  return result;
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
  PyArg_ParseTuple(args,"is",&frame,&cmd);
  APIEntry();
  MovieSetCommand(frame,cmd);
  APIExit();
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject *CmdMPlay(PyObject *self, 	PyObject *args)
{
  int cmd;
  PyArg_ParseTuple(args,"i",&cmd);
  APIEntry();
  MoviePlay(cmd);
  APIExit();
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject *CmdMMatrix(PyObject *self, 	PyObject *args)
{
  int cmd;
  PyArg_ParseTuple(args,"i",&cmd);
  APIEntry();
  MovieMatrix(cmd);
  APIExit();
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject *CmdMClear(PyObject *self, 	PyObject *args)
{
  APIEntry();
  MovieClearImages();
  APIExit();
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject *CmdRefresh(PyObject *self, 	PyObject *args)
{
  APIEntry();
  ExecutiveDrawNow();
  APIExit();
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject *CmdRefreshNow(PyObject *self, 	PyObject *args)
{
  APIEntry();
  ExecutiveDrawNow();
  MainRefreshNow();
  APIExit();
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject *CmdPNG(PyObject *self, 	PyObject *args)
{
  char *str1;
  PyArg_ParseTuple(args,"s",&str1);
  APIEntry();
  ExecutiveDrawNow();		
  ScenePNG(str1);
  APIExit();
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject *CmdMPNG(PyObject *self, 	PyObject *args)
{
  char *str1;
  PyArg_ParseTuple(args,"s",&str1);
  APIEntry();
  MoviePNG(str1,SettingGet(cSetting_cache_frames));
  APIExit();
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject *CmdMSet(PyObject *self, 	PyObject *args)
{
  char *str1;
  PyArg_ParseTuple(args,"s",&str1);
  APIEntry();
  MovieSequence(str1);
  SceneCountFrames();
  APIExit();
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject *CmdViewport(PyObject *self, 	PyObject *args)
{
  int w,h;
  PyArg_ParseTuple(args,"ii",&w,&h);
  if(w<10) w=10;
  if(h<10) h=10;

  w+=cOrthoRightSceneMargin;
  h+=cOrthoBottomSceneMargin;
  APIEntry();
  MainDoReshape(w,h); /* should be moved into Executive */
  APIExit();
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject *CmdFlag(PyObject *self, 	PyObject *args)
{
  char *str1;
  int flag;
  OrthoLineType s1;
  PyArg_ParseTuple(args,"is",&flag,&str1);
  APIEntry();
  SelectorGetTmp(str1,s1);
  ExecutiveFlag(flag,s1);
  SelectorFreeTmp(s1);
  APIExit();
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject *CmdColor(PyObject *self, 	PyObject *args)
{
  char *str1,*color;
  int flags;
  OrthoLineType s1;
  PyArg_ParseTuple(args,"ssi",&color,&str1,&flags);
  APIEntry();
  SelectorGetTmp(str1,s1);
  ExecutiveColor(s1,color,flags);
  SelectorFreeTmp(s1);
  APIExit();
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject *CmdColorDef(PyObject *self, 	PyObject *args)
{
  char *color;
  float v[3];
  PyArg_ParseTuple(args,"sfff",&color,v,v+1,v+2);
  APIEntry();
  ColorDef(color,v);
  APIExit();
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject *CmdRay(PyObject *self, 	PyObject *args)
{
  APIEntry();
  ExecutiveRay();
  APIExit();
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject *CmdClip(PyObject *self, 	PyObject *args)
{
  char *sname;
  float dist;
  PyArg_ParseTuple(args,"sf",&sname,&dist);
  APIEntry();
  switch(sname[0]) {
  case 'n':
	 SceneClip(0,dist);
	 break;
  case 'f':
	 SceneClip(1,dist);
	 break;
  }
  APIExit();
  Py_INCREF(Py_None);
  return Py_None;

}

static PyObject *CmdMove(PyObject *self, 	PyObject *args)
{
  char *sname;
  float dist;
  PyArg_ParseTuple(args,"sf",&sname,&dist);
  APIEntry();
  switch(sname[0]) {
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
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject *CmdTurn(PyObject *self, 	PyObject *args)
{
  char *sname;
  float angle;
  PyArg_ParseTuple(args,"sf",&sname,&angle);
  APIEntry();
  switch(sname[0]) {
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
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject *CmdSet(PyObject *self, 	PyObject *args)
{
  char *sname,*value;
  PyArg_ParseTuple(args,"ss",&sname,&value);
  APIEntry();
  ExecutiveSetSetting(sname,value);
  APIExit();
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject *CmdGet(PyObject *self, 	PyObject *args)
{
  float f;
  char *sname;
  PyObject *result;

  PyArg_ParseTuple(args,"s",&sname);
  APIEntry();
  f=SettingGetNamed(sname);
  APIExit();
  result = Py_BuildValue("f", f);
  return result;
}

static PyObject *CmdDelete(PyObject *self, 	PyObject *args)
{
  char *sname;

  PyArg_ParseTuple(args,"s",&sname);
  APIEntry();
  ExecutiveDelete(sname);
  APIExit();
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject *CmdShowHide(PyObject *self, 	PyObject *args)
{
  char *sname;
  int rep;
  int state;
  OrthoLineType s1;
  PyArg_ParseTuple(args,"sii",&sname,&rep,&state);
  APIEntry();
  if(sname[0]=='!') {
	 ExecutiveSetAllVisib(state);
  } else {
    SelectorGetTmp(sname,s1);
	 ExecutiveSetRepVisib(s1,rep,state);
	 SelectorFreeTmp(s1);
  }
  APIExit();
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject *CmdOnOff(PyObject *self, 	PyObject *args)
{
  char *name;
  int state;
  PyArg_ParseTuple(args,"si",&name,&state);
  APIEntry();
  ExecutiveSetObjVisib(name,state);
  APIExit();
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject *CmdQuit(PyObject *self, 	PyObject *args)
{
  APIEntry();
  PExit(EXIT_SUCCESS);
  APIExit();
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject *CmdSelect(PyObject *self, PyObject *args)
{
  char *sname,*sele;

  PyArg_ParseTuple(args,"ss",&sname,&sele);
  APIEntry();
  SelectorCreate(sname,sele,NULL,false);
  OrthoDirty();
  APIExit();
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject *CmdFinishObject(PyObject *self, PyObject *args)
{
  char *oname;
  Object *origObj = NULL;

  PyArg_ParseTuple(args,"s",&oname);

  APIEntry();
  origObj=ExecutiveFindObjectByName(oname);

  if(origObj) 
    ExecutiveUpdateObjectSelection(origObj);
  APIExit();
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject *CmdLoadObject(PyObject *self, PyObject *args)
{
  char *oname;
  PyObject *model;
  Object *origObj = NULL,*obj;
  OrthoLineType buf;
  int frame,type;
  int finish,discrete;

  buf[0]=0;

  PyArg_ParseTuple(args,"sOiiii",&oname,&model,&frame,&type,&finish,&discrete);

  APIEntry();
  origObj=ExecutiveFindObjectByName(oname);
  
      /* TODO check for existing object of wrong type */
  
  switch(type) {
  case cLoadTypeChemPyModel:
    PBlockAndUnlockAPI();
	 obj=(Object*)ObjectMoleculeLoadChemPyModel((ObjectMolecule*)origObj,model,frame,discrete);
    PLockAPIAndUnblock();
	 if(!origObj) {
	   if(obj) {
		 ObjectSetName(obj,oname);
		 ExecutiveManageObject(obj);
       if(frame<0)
         frame = ((ObjectMolecule*)obj)->NCSet-1;
		 sprintf(buf," CmdLoad: ChemPy-model loaded into object \"%s\", frame %d.\n",
               oname,frame+1);		  
	   }
	 } else if(origObj) {
      if(finish)
      ExecutiveUpdateObjectSelection(origObj);
       if(frame<0)
         frame = ((ObjectMolecule*)origObj)->NCSet-1;
		sprintf(buf," CmdLoad: ChemPy-model appended into object \"%s\", frame %d.\n",
              oname,frame+1);
	 }
	 break;
  }
  if(origObj) {
	 OrthoAddOutput(buf);
	 OrthoRestorePrompt();
  }
  APIExit();
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject *CmdLoadCoords(PyObject *self, PyObject *args)
{
  char *oname;
  PyObject *model;
  Object *origObj = NULL,*obj;
  OrthoLineType buf;
  int frame,type;

  buf[0]=0;

  PyArg_ParseTuple(args,"sOii",&oname,&model,&frame,&type);

  APIEntry();
  origObj=ExecutiveFindObjectByName(oname);
  
      /* TODO check for existing object of wrong type */
  if(!origObj)
    ErrMessage("LoadCoords","named object not found.");
  else 
    {
      switch(type) {
      case cLoadTypeChemPyModel:
        PBlockAndUnlockAPI();
        obj=(Object*)ObjectMoleculeLoadCoords((ObjectMolecule*)origObj,model,frame);
        PLockAPIAndUnblock();
        if(frame<0)
          frame=((ObjectMolecule*)obj)->NCSet-1;
        sprintf(buf," CmdLoad: Coordinates appended into object \"%s\", state %d.\n",
                oname,frame+1);
        break;
      }
    }
  if(origObj) {
	 OrthoAddOutput(buf);
	 OrthoRestorePrompt();
  }
  APIExit();
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject *CmdLoad(PyObject *self, PyObject *args)
{
  char *fname,*oname;
  Object *origObj = NULL,*obj;
  OrthoLineType buf;
  int frame,type;
  int finish,discrete;

  buf[0]=0;

  PyArg_ParseTuple(args,"ssiiii",&oname,&fname,&frame,&type,&finish,&discrete);

  APIEntry();
  origObj=ExecutiveFindObjectByName(oname);

      /* TODO check for existing object of wrong type */
  
  switch(type) {
  case cLoadTypePDB:
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
  case cLoadTypePDBStr:
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
  }
  if(origObj) {
	 OrthoAddOutput(buf);
	 OrthoRestorePrompt();
  }
  APIExit();
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject *CmdOrigin(PyObject *self, PyObject *args)
{
  char *str1;
  OrthoLineType s1;

  PyArg_ParseTuple(args,"s",&str1);
  APIEntry();
  SelectorGetTmp(str1,s1);
  ExecutiveCenter(s1,1);
  SelectorFreeTmp(s1);
  APIExit();
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject *CmdSort(PyObject *self, PyObject *args)
{
  char *name;
  PyArg_ParseTuple(args,"s",&name);
  APIEntry();
  ExecutiveSort(name);
  APIExit();
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject *CmdTest(PyObject *self, PyObject *args)
{
  Object *obj;
  APIEntry();
  obj=ExecutiveFindObjectByName("test");
  if(obj) ObjectMoleculeInferChemFromNeighGeom((ObjectMolecule*)obj,0);
  if(obj) ObjectMoleculeInferChemForProtein((ObjectMolecule*)obj,0);
  APIExit();
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject *CmdZoom(PyObject *self, PyObject *args)
{
  char *str1;
  OrthoLineType s1;

  PyArg_ParseTuple(args,"s",&str1);
  APIEntry();
  SelectorGetTmp(str1,s1);
  ExecutiveWindowZoom(s1);
  SelectorFreeTmp(s1);
  APIExit();
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject *CmdRemove(PyObject *self, PyObject *args)
{
  char *str1;
  OrthoLineType s1;

  PyArg_ParseTuple(args,"s",&str1);
  APIEntry();
  SelectorGetTmp(str1,s1);
  ExecutiveRemoveAtoms(s1);
  SelectorFreeTmp(s1);
  APIExit();
  Py_INCREF(Py_None);
  return Py_None;
}

void init_cmd(void)
{
  PyImport_AddModule("_cmd");
  Py_InitModule("_cmd", Cmd_methods);
}

