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
 
#include<stdlib.h>
#include<Python.h>

#include"MemoryDebug.h"
#include"Err.h"
#include"Util.h"
#include"PM.h"
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
#include"PUtils.h"
#include"Control.h"

#define tmpSele "_tmp"
#define tmpSele1 "_tmp1"
#define tmpSele2 "_tmp2"

PyObject *PM_Globals = NULL;

static PyObject *PMAlter(PyObject *self,   PyObject *args);
static PyObject *PMClip(PyObject *self, 	PyObject *args);
static PyObject *PMColor(PyObject *self, PyObject *args);
static PyObject *PMColorDef(PyObject *self, 	PyObject *args);
static PyObject *PMCopy(PyObject *self, PyObject *args);
static PyObject *PMCountStates(PyObject *self, PyObject *args);
static PyObject *PMDelete(PyObject *self, PyObject *args);
static PyObject *PMDirty(PyObject *self, 	PyObject *args);
static PyObject *PMDistance(PyObject *dummy, PyObject *args);
static PyObject *PMDo(PyObject *self, 	PyObject *args);
static PyObject *PMDump(PyObject *self, 	PyObject *args);
static PyObject *PMExportDots(PyObject *self, PyObject *args);
static PyObject *PMFit(PyObject *dummy, PyObject *args);
static PyObject *PMFitPairs(PyObject *dummy, PyObject *args);
static PyObject *PMIntraFit(PyObject *dummy, PyObject *args);
static PyObject *PMIsomesh(PyObject *self, 	PyObject *args);
static PyObject *PMFrame(PyObject *self, PyObject *args);
static PyObject *PMGet(PyObject *self, 	PyObject *args);
static PyObject *PMGetPDB(PyObject *dummy, PyObject *args);
static PyObject *PMGetFeedback(PyObject *dummy, PyObject *args);
static PyObject *PMGetGlobals(PyObject *dummy, PyObject *args);
static PyObject *PMGetMatrix(PyObject *self, 	PyObject *args);
static PyObject *PMGetMoment(PyObject *self, 	PyObject *args);
static PyObject *PMMem(PyObject *self, 	PyObject *args);
static PyObject *PMLoad(PyObject *self, 	PyObject *args);
static PyObject *PMMClear(PyObject *self, 	PyObject *args);
static PyObject *PMMDo(PyObject *self, 	PyObject *args);
static PyObject *PMMMatrix(PyObject *self, 	PyObject *args);
static PyObject *PMMove(PyObject *self, 	PyObject *args);
static PyObject *PMMPlay(PyObject *self, 	PyObject *args);
static PyObject *PMMPNG(PyObject *self, 	PyObject *args);
static PyObject *PMMSet(PyObject *self, 	PyObject *args);
static PyObject *PMOrigin(PyObject *self, PyObject *args);
static PyObject *PMOnOff(PyObject *self, 	PyObject *args);
static PyObject *PMOrient(PyObject *dummy, PyObject *args);
static PyObject *PMOverlap(PyObject *self, 	PyObject *args);
static PyObject *PMPNG(PyObject *self, 	PyObject *args);
static PyObject *PMSelect(PyObject *self, PyObject *args);
static PyObject *PMShowHide(PyObject *self, 	PyObject *args);
static PyObject *PMSetMatrix(PyObject *self, 	PyObject *args);
static PyObject *PMQuit(PyObject *self, 	PyObject *args);
static PyObject *PMReset(PyObject *self, PyObject *args);
static PyObject *PMRay(PyObject *self, 	PyObject *args);
static PyObject *PMResetRate(PyObject *dummy, PyObject *args);
static PyObject *PMRefresh(PyObject *self, 	PyObject *args);
static PyObject *PMRefreshNow(PyObject *self, 	PyObject *args);
static PyObject *PMReady(PyObject *dummy, PyObject *args);
static PyObject *PMRock(PyObject *self, PyObject *args);
static PyObject *PMRunPyMOL(PyObject *dummy, PyObject *args);
static PyObject *PMSystem(PyObject *dummy, PyObject *args);
static PyObject *PMSet(PyObject *self, 	PyObject *args);
static PyObject *PMSetFrame(PyObject *self, PyObject *args);
static PyObject *PMSetGlobals(PyObject *dummy, PyObject *args);
static PyObject *PMSort(PyObject *dummy, PyObject *args);
static PyObject *PMStereo(PyObject *self, PyObject *args);
static PyObject *PMTurn(PyObject *self, 	PyObject *args);
static PyObject *PMViewport(PyObject *self, 	PyObject *args);
static PyObject *PMZoom(PyObject *self, PyObject *args);

static PyMethodDef PM_methods[] = {
	{"alter",	     PMAlter,        METH_VARARGS },
	{"clip",	        PMClip,         METH_VARARGS },
	{"color",	     PMColor,        METH_VARARGS },
	{"colordef",	  PMColorDef,     METH_VARARGS },
	{"copy",         PMCopy,         METH_VARARGS },
	{"count_states", PMCountStates,  METH_VARARGS },
	{"delete",       PMDelete,       METH_VARARGS },
	{"dirty",        PMDirty,        METH_VARARGS },
	{"distance",	  PMDistance,     METH_VARARGS },
	{"do",	        PMDo,           METH_VARARGS },
	{"dump",	        PMDump,         METH_VARARGS },
	{"export_dots",  PMExportDots,   METH_VARARGS },
	{"fit",          PMFit,          METH_VARARGS },
	{"fit_pairs",    PMFitPairs,     METH_VARARGS },
	{"frame",	     PMFrame,        METH_VARARGS },
	{"get",	        PMGet,          METH_VARARGS },
	{"get_feedback", PMGetFeedback,  METH_VARARGS },
	{"get_globals",  PMGetGlobals,   METH_VARARGS },
	{"get_matrix",	  PMGetMatrix,    METH_VARARGS },
	{"get_moment",	  PMGetMoment,    METH_VARARGS },
	{"get_pdb",	     PMGetPDB,       METH_VARARGS },
	{"intrafit",     PMIntraFit,     METH_VARARGS },
	{"isomesh",	     PMIsomesh,      METH_VARARGS },
	{"load",	        PMLoad,         METH_VARARGS },
	{"mclear",	     PMMClear,       METH_VARARGS },
	{"mdo",	        PMMDo,          METH_VARARGS },
	{"mem",	        PMMem,          METH_VARARGS },
	{"move",	        PMMove,         METH_VARARGS },
	{"mset",	        PMMSet,         METH_VARARGS },
	{"mplay",	     PMMPlay,        METH_VARARGS },
	{"mpng_",	     PMMPNG,         METH_VARARGS },
	{"mmatrix",	     PMMMatrix,      METH_VARARGS },
	{"origin",	     PMOrigin,       METH_VARARGS },
	{"orient",	     PMOrient,       METH_VARARGS },
	{"onoff",        PMOnOff,        METH_VARARGS },
	{"overlap",      PMOverlap,      METH_VARARGS },
	{"png",	        PMPNG,          METH_VARARGS },
	{"quit",	        PMQuit,         METH_VARARGS },
	{"ready",        PMReady,        METH_VARARGS },
	{"refresh",      PMRefresh,      METH_VARARGS },
	{"refresh_now",  PMRefreshNow,   METH_VARARGS },
	{"render",	     PMRay,          METH_VARARGS },
	{"reset",        PMReset,        METH_VARARGS },
	{"reset_rate",	  PMResetRate,    METH_VARARGS },
	{"rock",	        PMRock,         METH_VARARGS },
	{"runpymol",	  PMRunPyMOL,     METH_VARARGS },
	{"select",       PMSelect,       METH_VARARGS },
	{"set",	        PMSet,          METH_VARARGS },
	{"setframe",	  PMSetFrame,     METH_VARARGS },
	{"showhide",     PMShowHide,     METH_VARARGS },
	{"set_globals",  PMSetGlobals,   METH_VARARGS },
	{"set_matrix",	  PMSetMatrix,    METH_VARARGS },
	{"sort",         PMSort,         METH_VARARGS },
	{"stereo",	     PMStereo,       METH_VARARGS },
	{"system",	     PMSystem,       METH_VARARGS },
	{"turn",	        PMTurn,         METH_VARARGS },
	{"viewport",     PMViewport,     METH_VARARGS },
	{"zoom",	        PMZoom,         METH_VARARGS },
	{NULL,		     NULL}		/* sentinel */
};

static PyObject *PMDump(PyObject *dummy, PyObject *args)
{
  char *str1,*str2;
  PyArg_ParseTuple(args,"ss",&str1,&str2);
  ExecutiveDump(str1,str2);
  Py_INCREF(Py_None);
  return Py_None;  
}

static PyObject *PMIsomesh(PyObject *self, 	PyObject *args) {
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
      ExecutiveGetBBox(s1,mn,mx);
      SelectorFreeTmp(s1);
      if(sscanf(str4,"%f",&fbuf)==1) {
        for(c=0;c<3;c++) {
          mn[c]-=fbuf;
          mx[c]+=fbuf;
        }
      }
      break;
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
  Py_INCREF(Py_None);
  return Py_None;  
}


static PyObject *PMOverlap(PyObject *dummy, PyObject *args)
{
  char *str1,*str2;
  int state1,state2;
  OrthoLineType s1,s2;
  PyObject *result;
  PyArg_ParseTuple(args,"ssii",&str1,&str2,&state1,&state2);
  SelectorGetTmp(str1,s1);
  SelectorGetTmp(str2,s2);
  result = Py_BuildValue("f",ExecutiveOverlap(s1,state1,s2,state2));
  SelectorFreeTmp(s1);
  SelectorFreeTmp(s2);
  return result;
}

static PyObject *PMDistance(PyObject *dummy, PyObject *args)
{
  char *str1,*str2;
  OrthoLineType s1,s2;
  PyObject *result;
  PyArg_ParseTuple(args,"ss",&str1,&str2);
  SelectorGetTmp(str1,s1);
  SelectorGetTmp(str2,s2);
  result = Py_BuildValue("f",ExecutiveDistance(s1,s2));
  SelectorFreeTmp(s1);
  SelectorFreeTmp(s2);
  return result;
}

static PyObject *PMAlter(PyObject *self,   PyObject *args)
{
  char *str1,*str2;
  OrthoLineType s1;

  PyArg_ParseTuple(args,"ss",&str1,&str2);
  SelectorGetTmp(str1,s1);
  ExecutiveAlter(s1,str2);
  SelectorFreeTmp(s1);
  Py_INCREF(Py_None);
  return Py_None;

}

static PyObject *PMCopy(PyObject *self,   PyObject *args)
{
  char *str1,*str2;
  PyArg_ParseTuple(args,"ss",&str1,&str2);
  ExecutiveCopy(str1,str2);
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject *PMResetRate(PyObject *dummy, PyObject *args)
{
  PyObject *result = NULL;
  result=Py_None;
  ButModeResetRate();
  Py_INCREF(result);
  return(result);
}

static PyObject *PMReady(PyObject *dummy, PyObject *args)
{
  PyObject *result = NULL;
  result = Py_BuildValue("i",PyMOLReady);
  return(result);
}

static PyObject *PMMem(PyObject *dummy, PyObject *args)
{
  PyObject *result = NULL;
  MemoryDebugDump();
  result = Py_None;
  Py_INCREF(result);
  return(result);
}

static PyObject *PMRunPyMOL(PyObject *dummy, PyObject *args)
{
  int gui;
  PyArg_ParseTuple(args,"i",&gui);
#ifdef _PYMOL_MODULE
  was_main(gui);
#endif
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject *PMCountStates(PyObject *dummy, PyObject *args)
{
  char *str1;
  OrthoLineType s1;
  PyObject *result;
  PyArg_ParseTuple(args,"s",&str1);
  SelectorGetTmp(str1,s1);
  result = Py_BuildValue("i",ExecutiveCountStates(s1));
  SelectorFreeTmp(s1); 
  return(result);
}

static PyObject *PMSystem(PyObject *dummy, PyObject *args)
{
  char *str1;
  PyObject *result = NULL;
  PyArg_ParseTuple(args,"s",&str1);
  result = Py_BuildValue("i",system(str1));
  return(result);
}

static PyObject *PMGetFeedback(PyObject *dummy, PyObject *args)
{
  OrthoLineType buffer;
  PyObject *result = NULL;

  if(OrthoFeedbackOut(buffer))
	result = Py_BuildValue("s",buffer);
  if(!result) {
	result=Py_None;
	Py_INCREF(result);
  }
  return(result);
}


static PyObject *PMGetPDB(PyObject *dummy, PyObject *args)
{
  char *str1;
  char *pdb = NULL;
  int state;
  OrthoLineType s1;

  PyObject *result = NULL;
  
  PyArg_ParseTuple(args,"si",&str1,&state);
  SelectorGetTmp(str1,s1);
  pdb=ExecutiveSeleToPDBStr(s1,state,true);
  SelectorFreeTmp(s1);
  if(pdb)
    result = Py_BuildValue("s",pdb);
  if(!result) {
    result=Py_None;
    Py_INCREF(result);
  }
  FreeP(pdb);
  return(result);
}


static PyObject *PMOrient(PyObject *dummy, PyObject *args)
{
  Matrix33d m;
  char *str1;
  OrthoLineType s1;

  PyArg_ParseTuple(args,"s",&str1);
  SelectorGetTmp(str1,s1);
  if(ExecutiveGetMoment(s1,m))
    ExecutiveOrient(s1,m);
  SelectorFreeTmp(s1); 
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject *PMFitPairs(PyObject *dummy, PyObject *args)
{
  PyObject *list;
  WordType *word = NULL;
  int ln=0;
  int a;
  int ok=true;
  PyObject *result = NULL;
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
    result=Py_BuildValue("f",ExecutiveRMSPairs(word,ln/2,2));
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

static PyObject *PMIntraFit(PyObject *dummy, PyObject *args)
{
  char *str1;
  int state;
  int mode;
  OrthoLineType s1;
  float *fVLA;
  PyObject *result=Py_None;
  PyArg_ParseTuple(args,"sii",&str1,&state,&mode);
  SelectorGetTmp(str1,s1);
  if(state<0) state=0;
  fVLA=ExecutiveRMSStates(s1,state,mode);
  if(fVLA) {
    result=PFloatVLAToPyList(fVLA);
    VLAFreeP(fVLA);
  }
  SelectorFreeTmp(s1);
  if(result==Py_None) Py_INCREF(result);
  return result;
}

static PyObject *PMFit(PyObject *dummy, PyObject *args)
{
  char *str1,*str2;
  int mode;
  OrthoLineType s1,s2;
  PyObject *result;
  PyArg_ParseTuple(args,"ssi",&str1,&str2,&mode);
  SelectorGetTmp(str1,s1);
  SelectorGetTmp(str2,s2);
  result=Py_BuildValue("f",ExecutiveRMS(s1,s2,mode));
  SelectorFreeTmp(s1);
  SelectorFreeTmp(s2);
  return result;
}

static PyObject *PMDirty(PyObject *self, 	PyObject *args)
{
  OrthoDirty();
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject *PMDo(PyObject *self, 	PyObject *args)
{
  char *str1;

  PyArg_ParseTuple(args,"s",&str1);
  if(str1[0]!='_') {
    OrthoAddOutput("PyMOL>");
    OrthoAddOutput(str1);
  }
  PParse(str1);
  OrthoNewLine(NULL);
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject *PMRock(PyObject *self, PyObject *args)
{
  ControlRock(-1);
  return Py_None;
}

static PyObject *PMGetMoment(PyObject *self, 	PyObject *args)
{
  Matrix33d m;
  PyObject *result;

  char *str1;
  PyArg_ParseTuple(args,"s",&str1);

  ExecutiveGetMoment(str1,m);
  result = Py_BuildValue("(ddd)(ddd)(ddd)", 
								 m[0][0],m[0][1],m[0][2],
								 m[1][0],m[1][1],m[1][2],
								 m[2][0],m[2][1],m[2][2]);
  return result;
}

static PyObject *PMExportDots(PyObject *self, 	PyObject *args)
{
  PyObject *result=NULL;
  PyObject *cObj;
  ExportDotsObj *obj;
  char *str1;
  int int1;
  
  PyArg_ParseTuple(args,"si",&str1,&int1);
  obj = ExportDots(str1,int1-1);
  if(obj) 
	 {
		cObj = PyCObject_FromVoidPtr(obj,(void(*)(void*))ExportDeleteMDebug);
		if(cObj) {
		  result = Py_BuildValue("O",cObj);
		  Py_DECREF(cObj); /* IMPORTANT */
		} 
	 }
  if(!result)
	 {
		result = Py_None;
		Py_INCREF(result);
	 }
  return result;
}

static PyObject *PMSetFrame(PyObject *self, PyObject *args)
{
  int mode,frm;
  PyArg_ParseTuple(args,"ii",&mode,&frm);
  SceneSetFrame(mode,frm);
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject *PMFrame(PyObject *self, PyObject *args)
{
  int frm;
  PyArg_ParseTuple(args,"i",&frm);
  frm--;
  if(frm<0) frm=0;
  SceneSetFrame(0,frm);
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject *PMStereo(PyObject *self, PyObject *args)
{
  PyObject *result;
  int i1;
  if(StereoCapable) {
  PyArg_ParseTuple(args,"i",&i1);
  ExecutiveStereo(i1);
  result=Py_BuildValue("i",1);
  } else {
  result=Py_BuildValue("i",0);
  }
  return result;
}

static PyObject *PMReset(PyObject *self, PyObject *args)
{
  int cmd;
  PyArg_ParseTuple(args,"i",&cmd);
  ExecutiveReset(cmd);
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject *PMSetMatrix(PyObject *self, 	PyObject *args)
{
  float m[16];
  PyArg_ParseTuple(args,"ffffffffffffffff",
						 &m[0],&m[1],&m[2],&m[3],
						 &m[4],&m[5],&m[6],&m[7],
						 &m[8],&m[9],&m[10],&m[11],
						 &m[12],&m[13],&m[14],&m[15]);
  SceneSetMatrix(m);
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject *PMGetMatrix(PyObject *self, 	PyObject *args)
{
  float *f;
  PyObject *result;

  f=SceneGetMatrix();
  result = Py_BuildValue("ffffffffffffffff", 
								 f[0],f[1],f[2],f[3],
								 f[4],f[5],f[6],f[7],
								 f[8],f[9],f[10],f[11],
								 f[12],f[13],f[14],f[15]
								 );
  return result;
}

static PyObject *PMMDo(PyObject *self, 	PyObject *args)
{
  char *cmd;
  int frame;
  PyArg_ParseTuple(args,"is",&frame,&cmd);
  MovieSetCommand(frame,cmd);
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject *PMMPlay(PyObject *self, 	PyObject *args)
{
  int cmd;
  PyArg_ParseTuple(args,"i",&cmd);
  MoviePlay(cmd);
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject *PMMMatrix(PyObject *self, 	PyObject *args)
{
  int cmd;
  PyArg_ParseTuple(args,"i",&cmd);
  MovieMatrix(cmd);
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject *PMSetGlobals(PyObject *dummy, PyObject *args)
{
  PyObject *result = NULL;
  PyObject *globals;
  
  PyArg_ParseTuple(args, "O", &globals);
  if(globals) {
	 PM_Globals = globals;
  }

  result = Py_None;
  Py_INCREF(result);
  return result;
}

static PyObject *PMGetGlobals(PyObject *dummy, PyObject *args)
{
  PyObject *result = PM_Globals;
  Py_INCREF(result);
  return result;
}

static PyObject *PMMClear(PyObject *self, 	PyObject *args)
{
  MovieClearImages();
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject *PMRefresh(PyObject *self, 	PyObject *args)
{
  ExecutiveDrawNow();
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject *PMRefreshNow(PyObject *self, 	PyObject *args)
{
  ExecutiveDrawNow();
  MainRefreshNow();
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject *PMPNG(PyObject *self, 	PyObject *args)
{
  char *str1;
  PyArg_ParseTuple(args,"s",&str1);
  ExecutiveDrawNow();
  ScenePNG(str1);
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject *PMMPNG(PyObject *self, 	PyObject *args)
{
  char *str1;
  PyArg_ParseTuple(args,"s",&str1);
  MoviePNG(str1,SettingGet(cSetting_cache_frames));
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject *PMMSet(PyObject *self, 	PyObject *args)
{
  char *str1;
  PyArg_ParseTuple(args,"s",&str1);
  MovieSequence(str1);
  SceneCountFrames();
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject *PMViewport(PyObject *self, 	PyObject *args)
{
  int w,h;
  PyArg_ParseTuple(args,"ii",&w,&h);
  if(w<10) w=10;
  if(h<10) h=10;

  w+=cOrthoRightSceneMargin;
  h+=cOrthoBottomSceneMargin;
  MainDoReshape(w,h); /* should be moved into Executive */
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject *PMColor(PyObject *self, 	PyObject *args)
{
  char *str1,*color;
  int flags;
  OrthoLineType s1;
  PyArg_ParseTuple(args,"ssi",&color,&str1,&flags);
  SelectorGetTmp(str1,s1);
  ExecutiveColor(s1,color,flags);
  SelectorFreeTmp(s1);
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject *PMColorDef(PyObject *self, 	PyObject *args)
{
  char *color;
  float v[3];
  PyArg_ParseTuple(args,"sfff",&color,v,v+1,v+2);
  ColorDef(color,v);
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject *PMRay(PyObject *self, 	PyObject *args)
{
  PyThreadState *_save;
  Py_UNBLOCK_THREADS;

  ExecutiveRay();

  Py_BLOCK_THREADS;
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject *PMClip(PyObject *self, 	PyObject *args)
{
  char *sname;
  float dist;
  PyArg_ParseTuple(args,"sf",&sname,&dist);
  switch(sname[0]) {
  case 'n':
	 SceneClip(0,dist);
	 break;
  case 'f':
	 SceneClip(1,dist);
	 break;
  }
  Py_INCREF(Py_None);
  return Py_None;

}

static PyObject *PMMove(PyObject *self, 	PyObject *args)
{
  char *sname;
  float dist;
  PyArg_ParseTuple(args,"sf",&sname,&dist);
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
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject *PMTurn(PyObject *self, 	PyObject *args)
{
  char *sname;
  float angle;
  PyArg_ParseTuple(args,"sf",&sname,&angle);
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
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject *PMSet(PyObject *self, 	PyObject *args)
{
  char *sname,*value;
  PyArg_ParseTuple(args,"ss",&sname,&value);
  ExecutiveSetSetting(sname,value);
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject *PMGet(PyObject *self, 	PyObject *args)
{
  float f;
  char *sname;
  PyObject *result;

  PyArg_ParseTuple(args,"s",&sname);
  f=SettingGetNamed(sname);
  result = Py_BuildValue("f", f);
  return result;
}

static PyObject *PMDelete(PyObject *self, 	PyObject *args)
{
  char *sname;

  PyArg_ParseTuple(args,"s",&sname);
  ExecutiveDelete(sname);
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject *PMShowHide(PyObject *self, 	PyObject *args)
{
  char *sname;
  int rep;
  int state;
  OrthoLineType s1;
  PyArg_ParseTuple(args,"sii",&sname,&rep,&state);
  if(sname[0]=='!') {
	 ExecutiveSetAllVisib(state);
  } else {
    SelectorGetTmp(sname,s1);
	 ExecutiveSetRepVisib(s1,rep,state);
	 SelectorFreeTmp(s1);
  }
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject *PMOnOff(PyObject *self, 	PyObject *args)
{
  char *name;
  int state;
  PyArg_ParseTuple(args,"si",&name,&state);
  ExecutiveSetObjVisib(name,state);
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject *PMQuit(PyObject *self, 	PyObject *args)
{
  PExit(EXIT_SUCCESS);
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject *PMSelect(PyObject *self, PyObject *args)
{
  char *sname,*sele;

  PyArg_ParseTuple(args,"ss",&sname,&sele);
  SelectorCreate(sname,sele,NULL);
  OrthoDirty();
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject *PMLoad(PyObject *self, PyObject *args)
{
  char *fname,*oname;
  Object *origObj = NULL,*obj;
  OrthoLineType buf;
  int frame,type;

  buf[0]=0;
  #define cLoadTypePDB 0
  #define cLoadTypeMOL 1
  #define cLoadTypeSDF 2
  #define cLoadTypeMOLStr 3
  #define cLoadTypeMMD 4
  #define cLoadTypeMMDSeparate 5
  #define cLoadTypeMMDStr 6
  #define cLoadTypeXPLORMap 7

  PyArg_ParseTuple(args,"ssii",&oname,&fname,&frame,&type);

  origObj=ExecutiveFindObjectByName(oname);

      /* TODO check for existing object of wrong type */
  
  switch(type) {
  case cLoadTypePDB:
	 if(!origObj) {
		obj=(Object*)ObjectMoleculeLoadPDBFile(NULL,fname,frame);
		if(obj) {
		  ObjectSetName(obj,oname);
		  ExecutiveManageObject(obj);
		  sprintf(buf," PMLoad: \"%s\" loaded into object \"%s\".\n",fname,oname);
		}
	 } else {
		ObjectMoleculeLoadPDBFile((ObjectMolecule*)origObj,fname,frame);
		ExecutiveUpdateObjectSelection(origObj);
		sprintf(buf," PMLoad: \"%s\" appended into object \"%s\".\n",fname,oname);
	 }
	 break;
  case cLoadTypeMOL:
	 obj=(Object*)ObjectMoleculeLoadMOLFile((ObjectMolecule*)origObj,fname,frame);
	 if(!origObj) {
	   if(obj) {
		 ObjectSetName(obj,oname);
		 ExecutiveManageObject(obj);
		 sprintf(buf," PMLoad: \"%s\" loaded into object \"%s\".\n",fname,oname);		  
	   }
	 } else if(origObj) {
		ExecutiveUpdateObjectSelection(origObj);
		sprintf(buf," PMLoad: \"%s\" appended into object \"%s\".\n",fname,oname);
	 }
	 break;
  case cLoadTypeMOLStr:
	 obj=(Object*)ObjectMoleculeReadMOLStr((ObjectMolecule*)origObj,fname,frame);
	 if(!origObj) {
	   if(obj) {
		 ObjectSetName(obj,oname);
		 ExecutiveManageObject(obj);
		 sprintf(buf," PMLoad: MOL-string loaded into object \"%s\".\n",oname);		  
	   }
	 } else if(origObj) {
		ExecutiveUpdateObjectSelection(origObj);
		sprintf(buf," PMLoad: MOL-string appended into object \"%s\".\n",oname);
	 }
	 break;
  case cLoadTypeMMD:
	 obj=(Object*)ObjectMoleculeLoadMMDFile((ObjectMolecule*)origObj,fname,frame,NULL);
	 if(!origObj) {
	   if(obj) {
        ObjectSetName(obj,oname);
        ExecutiveManageObject(obj);
        sprintf(buf," PMLoad: \"%s\" loaded into object \"%s\".\n",fname,oname);		  
	   }
	 } else if(origObj) {
		ExecutiveUpdateObjectSelection(origObj);
		sprintf(buf," PMLoad: \"%s\" appended into object \"%s\".\n",fname,oname);
	 }
    break;
  case cLoadTypeMMDSeparate:
	 ObjectMoleculeLoadMMDFile((ObjectMolecule*)origObj,fname,frame,oname);
    break;
  case cLoadTypeMMDStr:
	 obj=(Object*)ObjectMoleculeReadMMDStr((ObjectMolecule*)origObj,fname,frame);
	 if(!origObj) {
	   if(obj) {
		 ObjectSetName(obj,oname);
		 ExecutiveManageObject(obj);
		 sprintf(buf," PMLoad: MMD-string loaded into object \"%s\".\n",oname);		  
	   }
	 } else if(origObj) {
		ExecutiveUpdateObjectSelection(origObj);
		sprintf(buf," PMLoad: MMD-string appended into object \"%s\".\n",oname);
	 }
	 break;
  case cLoadTypeXPLORMap:
	 if(!origObj) {
		obj=(Object*)ObjectMapLoadXPLORFile(NULL,fname,frame);
		if(obj) {
		  ObjectSetName(obj,oname);
		  ExecutiveManageObject((Object*)obj);
		  sprintf(buf," PMLoad: \"%s\" loaded into object \"%s\".\n",fname,oname);
		}
	 } else {
		ObjectMapLoadXPLORFile((ObjectMap*)origObj,fname,frame);
		sprintf(buf," PMLoad: \"%s\" appended into object \"%s\".\n",
				  fname,oname);
	 }
	 break;
  }
  if(origObj) {
	 OrthoAddOutput(buf);
	 OrthoRestorePrompt();
  }
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject *PMOrigin(PyObject *self, PyObject *args)
{
  char *str1;
  OrthoLineType s1;

  PyArg_ParseTuple(args,"s",&str1);
  SelectorGetTmp(str1,s1);
  ExecutiveCenter(s1,1);
  SelectorFreeTmp(s1);
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject *PMSort(PyObject *self, PyObject *args)
{
  char *name;
  PyArg_ParseTuple(args,"s",&name);
  ExecutiveSort(name);
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject *PMZoom(PyObject *self, PyObject *args)
{
  char *str1;
  OrthoLineType s1;

  PyArg_ParseTuple(args,"s",&str1);
  SelectorGetTmp(str1,s1);
  ExecutiveCenter(s1,1);
  ExecutiveWindowZoom(s1);
  SelectorFreeTmp(s1);
  Py_INCREF(Py_None);
  return Py_None;
}

#ifdef _PYMOL_MODULE
void PMInit(void) {}
void init_pm(void);
void init_pm(void)
#else
void PMInit(void)
#endif
{
  PyImport_AddModule("_pm");
  Py_InitModule("_pm", PM_methods);
}
