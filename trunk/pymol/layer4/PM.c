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

   Thus, pm.py needs to guard against direct invocation of certain _pm
   (API) calls from outside threads [Examples: _pm.png(..),
   _pmmpng(..) ].  Instead, it needs to hand them over to the main
   thread by way of a pm.mdo(...)  statement.

   Note that current, most glut operations have been pushed into the
   main event and redraw loops to avoid these kinds of problems - so
   I'm not sure how important this really is anymore.

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
#include"PConv.h"
#include"Control.h"

#define tmpSele "_tmp"
#define tmpSele1 "_tmp1"
#define tmpSele2 "_tmp2"

PyObject *PM_Globals = NULL;

static void APIEntry(void)
{
  P_glut_thread_keep_out++;
  PUnblock();
}

static void APIExit(void)
{
  P_glut_thread_keep_out--;
  PBlock();
}

static PyObject *PMAlter(PyObject *self,   PyObject *args);
static PyObject *PMClip(PyObject *self, 	PyObject *args);
static PyObject *PMColor(PyObject *self, PyObject *args);
static PyObject *PMColorDef(PyObject *self, 	PyObject *args);
static PyObject *PMCopy(PyObject *self, PyObject *args);
static PyObject *PMCountStates(PyObject *self, PyObject *args);
static PyObject *PMDelete(PyObject *self, PyObject *args);
static PyObject *PMDirty(PyObject *self, 	PyObject *args);
static PyObject *PMDist(PyObject *dummy, PyObject *args);
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
static PyObject *PMQuit(PyObject *self, 	PyObject *args);
static PyObject *PMReset(PyObject *self, PyObject *args);
static PyObject *PMRay(PyObject *self, 	PyObject *args);
static PyObject *PMResetRate(PyObject *dummy, PyObject *args);
static PyObject *PMRefresh(PyObject *self, 	PyObject *args);
static PyObject *PMRefreshNow(PyObject *self, 	PyObject *args);
static PyObject *PMReady(PyObject *dummy, PyObject *args);
static PyObject *PMRock(PyObject *self, PyObject *args);
static PyObject *PMRunPyMOL(PyObject *dummy, PyObject *args);
static PyObject *PMSelect(PyObject *self, PyObject *args);
static PyObject *PMSetMatrix(PyObject *self, 	PyObject *args);
static PyObject *PMSet(PyObject *self, 	PyObject *args);
static PyObject *PMSetFrame(PyObject *self, PyObject *args);
static PyObject *PMSetGlobals(PyObject *dummy, PyObject *args);
static PyObject *PMShowHide(PyObject *self, 	PyObject *args);
static PyObject *PMSort(PyObject *dummy, PyObject *args);
static PyObject *PMStereo(PyObject *self, PyObject *args);
static PyObject *PMSystem(PyObject *dummy, PyObject *args);
static PyObject *PMSymExp(PyObject *dummy, PyObject *args);
static PyObject *PMTest(PyObject *self, 	PyObject *args);
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
	{"dist",    	  PMDist,         METH_VARARGS },
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
	{"symexp",	     PMSymExp,       METH_VARARGS },
	{"test",	        PMTest,         METH_VARARGS },
	{"turn",	        PMTurn,         METH_VARARGS },
	{"viewport",     PMViewport,     METH_VARARGS },
	{"zoom",	        PMZoom,         METH_VARARGS },
	{NULL,		     NULL}		/* sentinel */
};

static PyObject *PMDump(PyObject *dummy, PyObject *args)
{
  char *str1,*str2;
  PyArg_ParseTuple(args,"ss",&str1,&str2);
  APIEntry();
  ExecutiveDump(str1,str2);
  APIExit();
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
  APIExit();
  Py_INCREF(Py_None);
  return Py_None;  
}

static PyObject *PMSymExp(PyObject *self, 	PyObject *args) {
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



static PyObject *PMOverlap(PyObject *dummy, PyObject *args)
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

static PyObject *PMDist(PyObject *dummy, PyObject *args)
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

static PyObject *PMDistance(PyObject *dummy, PyObject *args)
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

static PyObject *PMAlter(PyObject *self,   PyObject *args)
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

static PyObject *PMCopy(PyObject *self,   PyObject *args)
{
  char *str1,*str2;
  PyArg_ParseTuple(args,"ss",&str1,&str2);
  APIEntry();
  ExecutiveCopy(str1,str2);
  APIExit();
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject *PMResetRate(PyObject *dummy, PyObject *args)
{
  PyObject *result = NULL;
  result=Py_None;
  APIEntry();
  ButModeResetRate();
  APIExit();
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

static PyObject *PMSystem(PyObject *dummy, PyObject *args)
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

static PyObject *PMGetFeedback(PyObject *dummy, PyObject *args)
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


static PyObject *PMGetPDB(PyObject *dummy, PyObject *args)
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


static PyObject *PMOrient(PyObject *dummy, PyObject *args)
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

static PyObject *PMFitPairs(PyObject *dummy, PyObject *args)
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

static PyObject *PMIntraFit(PyObject *dummy, PyObject *args)
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

static PyObject *PMFit(PyObject *dummy, PyObject *args)
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

static PyObject *PMDirty(PyObject *self, 	PyObject *args)
{
  APIEntry();
  OrthoDirty();
  APIExit();
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject *PMDo(PyObject *self, 	PyObject *args)
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

static PyObject *PMRock(PyObject *self, PyObject *args)
{
  APIEntry();
  ControlRock(-1);
  APIExit();
  return Py_None;
}

static PyObject *PMGetMoment(PyObject *self, 	PyObject *args)
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

static PyObject *PMExportDots(PyObject *self, 	PyObject *args)
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
  APIEntry();
  SceneSetFrame(mode,frm);
  APIExit();
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject *PMFrame(PyObject *self, PyObject *args)
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

static PyObject *PMStereo(PyObject *self, PyObject *args)
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

static PyObject *PMReset(PyObject *self, PyObject *args)
{
  int cmd;
  PyArg_ParseTuple(args,"i",&cmd);
  APIEntry();
  ExecutiveReset(cmd);
  APIExit();
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
  APIEntry();
  SceneSetMatrix(m);
  APIExit();
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject *PMGetMatrix(PyObject *self, 	PyObject *args)
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

static PyObject *PMMDo(PyObject *self, 	PyObject *args)
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

static PyObject *PMMPlay(PyObject *self, 	PyObject *args)
{
  int cmd;
  PyArg_ParseTuple(args,"i",&cmd);
  APIEntry();
  MoviePlay(cmd);
  APIExit();
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject *PMMMatrix(PyObject *self, 	PyObject *args)
{
  int cmd;
  PyArg_ParseTuple(args,"i",&cmd);
  APIEntry();
  MovieMatrix(cmd);
  APIExit();
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
  APIEntry();
  MovieClearImages();
  APIExit();
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject *PMRefresh(PyObject *self, 	PyObject *args)
{
  APIEntry();
  ExecutiveDrawNow();
  APIExit();
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject *PMRefreshNow(PyObject *self, 	PyObject *args)
{
  APIEntry();
  ExecutiveDrawNow();
  MainRefreshNow();
  APIExit();
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject *PMPNG(PyObject *self, 	PyObject *args)
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

static PyObject *PMMPNG(PyObject *self, 	PyObject *args)
{
  char *str1;
  PyArg_ParseTuple(args,"s",&str1);
  APIEntry();
  MoviePNG(str1,SettingGet(cSetting_cache_frames));
  APIExit();
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject *PMMSet(PyObject *self, 	PyObject *args)
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

static PyObject *PMViewport(PyObject *self, 	PyObject *args)
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

static PyObject *PMColor(PyObject *self, 	PyObject *args)
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

static PyObject *PMColorDef(PyObject *self, 	PyObject *args)
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

static PyObject *PMRay(PyObject *self, 	PyObject *args)
{
  APIEntry();
  ExecutiveRay();
  APIExit();
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject *PMClip(PyObject *self, 	PyObject *args)
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

static PyObject *PMMove(PyObject *self, 	PyObject *args)
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

static PyObject *PMTurn(PyObject *self, 	PyObject *args)
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

static PyObject *PMSet(PyObject *self, 	PyObject *args)
{
  char *sname,*value;
  PyArg_ParseTuple(args,"ss",&sname,&value);
  APIEntry();
  ExecutiveSetSetting(sname,value);
  APIExit();
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject *PMGet(PyObject *self, 	PyObject *args)
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

static PyObject *PMDelete(PyObject *self, 	PyObject *args)
{
  char *sname;

  PyArg_ParseTuple(args,"s",&sname);
  APIEntry();
  ExecutiveDelete(sname);
  APIExit();
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

static PyObject *PMOnOff(PyObject *self, 	PyObject *args)
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
  APIEntry();
  SelectorCreate(sname,sele,NULL);
  OrthoDirty();
  APIExit();
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

  APIEntry();
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
  APIExit();
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject *PMOrigin(PyObject *self, PyObject *args)
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

static PyObject *PMSort(PyObject *self, PyObject *args)
{
  char *name;
  PyArg_ParseTuple(args,"s",&name);
  APIEntry();
  ExecutiveSort(name);
  APIExit();
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject *PMTest(PyObject *self, PyObject *args)
{
  Object *obj;
  APIEntry();
  obj=ExecutiveFindObjectByName("test");
  if(obj) ObjectMoleculeBlindSymMovie((ObjectMolecule*)obj);
  APIExit();
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject *PMZoom(PyObject *self, PyObject *args)
{
  char *str1;
  OrthoLineType s1;

  PyArg_ParseTuple(args,"s",&str1);
  APIEntry();
  SelectorGetTmp(str1,s1);
  ExecutiveCenter(s1,1);
  ExecutiveWindowZoom(s1);
  SelectorFreeTmp(s1);
  APIExit();
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
