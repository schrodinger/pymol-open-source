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

#include"Err.h"
#include"Util.h"
#include"PM.h"
#include"ObjectMolecule.h"
#include"Executive.h"
#include"Selector.h"
#include"main.h"
#include"Scene.h"
#include"Setting.h"
#include"Movie.h"
#include"Export.h"
#include"PUtils.h"

#define tmpSele "_tmp"

PyObject *PM_Globals = NULL;

static PyObject *PMQuit(PyObject *self, 	PyObject *args);
static PyObject *PMLoad(PyObject *self, 	PyObject *args);
static PyObject *PMOrigin(PyObject *self, PyObject *args);
static PyObject *PMZoom(PyObject *self, PyObject *args);
static PyObject *PMSelect(PyObject *self, PyObject *args);
static PyObject *PMDelete(PyObject *self, PyObject *args);
static PyObject *PMColor(PyObject *self, PyObject *args);
static PyObject *PMSet(PyObject *self, 	PyObject *args);
static PyObject *PMGet(PyObject *self, 	PyObject *args);
static PyObject *PMRay(PyObject *self, 	PyObject *args);
static PyObject *PMTurn(PyObject *self, 	PyObject *args);
static PyObject *PMClip(PyObject *self, 	PyObject *args);
static PyObject *PMMove(PyObject *self, 	PyObject *args);
static PyObject *PMSetGlobals(PyObject *dummy, PyObject *args);
static PyObject *PMMSet(PyObject *self, 	PyObject *args);
static PyObject *PMRefresh(PyObject *self, 	PyObject *args);
static PyObject *PMMClear(PyObject *self, 	PyObject *args);
static PyObject *PMMDo(PyObject *self, 	PyObject *args);
static PyObject *PMMPlay(PyObject *self, 	PyObject *args);
static PyObject *PMMMatrix(PyObject *self, 	PyObject *args);
static PyObject *PMViewport(PyObject *self, 	PyObject *args);
static PyObject *PMMPNG(PyObject *self, 	PyObject *args);
static PyObject *PMPNG(PyObject *self, 	PyObject *args);
static PyObject *PMShowHide(PyObject *self, 	PyObject *args);
static PyObject *PMOnOff(PyObject *self, 	PyObject *args);
static PyObject *PMSetMatrix(PyObject *self, 	PyObject *args);
static PyObject *PMGetMatrix(PyObject *self, 	PyObject *args);
static PyObject *PMReset(PyObject *self, PyObject *args);
static PyObject *PMFrame(PyObject *self, PyObject *args);
static PyObject *PMExportDots(PyObject *self, PyObject *args);

static PyMethodDef PM_methods[] = {
	{"quit",	  PMQuit,   METH_VARARGS },
	{"load",	  PMLoad,   METH_VARARGS },
	{"showhide",    PMShowHide,   METH_VARARGS },
	{"onoff",    PMOnOff,   METH_VARARGS },
	{"select",    PMSelect,   METH_VARARGS },
	{"viewport",    PMViewport,   METH_VARARGS },
	{"refresh",    PMRefresh,   METH_VARARGS },
	{"delete",    PMDelete,   METH_VARARGS },
	{"origin",	 PMOrigin,   METH_VARARGS },
	{"zoom",	 PMZoom,   METH_VARARGS },
	{"color",	 PMColor,   METH_VARARGS },
	{"set",	 PMSet,   METH_VARARGS },
	{"get",	 PMGet,   METH_VARARGS },
	{"render",	 PMRay,   METH_VARARGS },
	{"reset",	 PMReset,   METH_VARARGS },
	{"frame",	 PMFrame,   METH_VARARGS },
	{"turn",	 PMTurn,   METH_VARARGS },
	{"clip",	 PMClip,   METH_VARARGS },
	{"move",	 PMMove,   METH_VARARGS },
	{"set_globals",	 PMSetGlobals,   METH_VARARGS },
	{"mset",	 PMMSet,   METH_VARARGS },
	{"mclear",	 PMMClear,   METH_VARARGS },
	{"mplay",	 PMMPlay,   METH_VARARGS },
	{"mpng",	 PMMPNG,   METH_VARARGS },
	{"mdo",	 PMMDo,   METH_VARARGS },
	{"mmatrix",	 PMMMatrix,   METH_VARARGS },
	{"png",	 PMPNG,   METH_VARARGS },
	{"set_matrix",	 PMSetMatrix,   METH_VARARGS },
	{"get_matrix",	 PMGetMatrix,   METH_VARARGS },
	{"export_dots", PMExportDots,   METH_VARARGS },
	{NULL,		NULL}		/* sentinel */
};

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
  int a;
  double md[16];
  float m[16];
  PyArg_ParseTuple(args,"dddddddddddddddd",
						 &md[0],&md[1],&md[2],&md[3],
						 &md[4],&md[5],&md[6],&md[7],
						 &md[8],&md[9],&md[10],&md[11],
						 &md[12],&md[13],&md[14],&md[15]);
  for(a=0;a<16;a++)
	 m[a]=md[a];
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
  int flag;
  PyArg_ParseTuple(args,"si",&str1,&flag);
  MoviePNG(str1,flag);
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
  MainReshape(w,h); /* should be moved into Executive */
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject *PMColor(PyObject *self, 	PyObject *args)
{
  char *sname,*color;
  int flags;
  PyArg_ParseTuple(args,"ssi",&color,&sname,&flags);
  if(sname[0]=='(') {
	 SelectorCreate(tmpSele,sname,NULL);
	 ExecutiveColor(tmpSele,color,flags);
	 ExecutiveDelete(tmpSele);
  } else {
	 ExecutiveColor(sname,color,flags);
  }
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject *PMRay(PyObject *self, 	PyObject *args)
{
  /*  PyArg_ParseTuple(args,"ss",&sname,&value);
		ExecutiveSetSetting(sname,value);*/
  ExecutiveRay();
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject *PMClip(PyObject *self, 	PyObject *args)
{
  char *sname;
  char *value;
  float dist;
  PyArg_ParseTuple(args,"ss",&sname,&value);
  sscanf(value,"%f",&dist);
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
  char *value;
  float dist;
  PyArg_ParseTuple(args,"ss",&sname,&value);
  sscanf(value,"%f",&dist);
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
  char *value;
  float angle;
  PyArg_ParseTuple(args,"ss",&sname,&value);
  sscanf(value,"%f",&angle);
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
  PyArg_ParseTuple(args,"sii",&sname,&rep,&state);
  if(sname[0]=='!') {
	 ExecutiveSetAllVisib(state);
  } else if(sname[0]=='(') {
	 SelectorCreate(tmpSele,sname,NULL);
	 ExecutiveSetRepVisib(tmpSele,rep,state);
	 ExecutiveDelete(tmpSele);
  } else {
	 ExecutiveSetRepVisib(sname,rep,state);
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
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject *PMLoad(PyObject *self, PyObject *args)
{
  char *fname,*oname;
  ObjectMolecule *objMol = NULL,*obj;
  OrthoLineType buf;
  int frame,type;

  buf[0]=0;
  #define cLoadTypePDB 0
  #define cLoadTypeMOL 1
  #define cLoadTypeSDF 2
  #define cLoadTypeMOLStr 3

  PyArg_ParseTuple(args,"ssii",&oname,&fname,&frame,&type);

  objMol=(ObjectMolecule*)ExecutiveFindObjectByName(oname);
  
  switch(type) {
  case cLoadTypePDB:
	 if(!objMol) {
		objMol=ObjectMoleculeLoadPDBFile(NULL,fname,frame);
		fflush(stdout);
		if(objMol) {
		  ObjectSetName((Object*)objMol,oname);
		  ExecutiveManageObject((Object*)objMol);
		  sprintf(buf,"PMLoad: \"%s\" loaded into object \"%s\".\n",
					 fname,oname);
		}
	 } else {
		ObjectMoleculeLoadPDBFile(objMol,fname,frame);
		ExecutiveUpdateObjectSelection((Object*)objMol);
		sprintf(buf,"PMAppend: \"%s\" appended into object \"%s\".\n",
				  fname,oname);
	 }
	 break;
  case cLoadTypeMOL:
	 obj=ObjectMoleculeLoadMOLFile(objMol,fname,frame);
	 if(!objMol) {
		ObjectSetName((Object*)obj,oname);
		ExecutiveManageObject((Object*)obj);
		sprintf(buf,"PMLoad: \"%s\" loaded into object \"%s\".\n",
				  fname,oname);		  
	 } else if(objMol) {
		ExecutiveUpdateObjectSelection((Object*)objMol);
		sprintf(buf,"PMAppend: \"%s\" appended into object \"%s\".\n",
				  fname,oname);
	 }
	 break;
  case cLoadTypeMOLStr:
	 obj=ObjectMoleculeReadMOLStr(objMol,fname,frame);
	 if(!objMol) {
		ObjectSetName((Object*)obj,oname);
		ExecutiveManageObject((Object*)obj);
		sprintf(buf,"PMLoad: MOL-string loaded into object \"%s\".\n",
				  oname);		  
	 } else if(objMol) {
		ExecutiveUpdateObjectSelection((Object*)objMol);
		sprintf(buf,"PMAppend: MOL-string appended into object \"%s\".\n",
				  oname);
	 }
	 break;
  }

  if(objMol) {
	 OrthoAddOutput(buf);
	 OrthoRestorePrompt();
  }

  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject *PMOrigin(PyObject *self, PyObject *args)
{
  char *name;
  PyArg_ParseTuple(args,"s",&name);
  ExecutiveCenter(name,1);
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject *PMZoom(PyObject *self, PyObject *args)
{
  char *name;
  PyArg_ParseTuple(args,"s",&name);
  ExecutiveCenter(name,1);
  ExecutiveWindowZoom(name);
  Py_INCREF(Py_None);
  return Py_None;
}

void PMInit(void)

{
  PyImport_AddModule("pmx");
  Py_InitModule("pmx", PM_methods);
}
