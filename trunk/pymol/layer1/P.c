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

#ifdef WIN32
#include<windows.h>
#endif

#include"os_std.h"
#include"os_time.h"
#include"os_unix.h"

#include<Python.h>


#include"MemoryDebug.h"
#include"Base.h"
#include"Err.h"
#include"P.h"
#include"PConv.h"
#include"Ortho.h"
#include"Cmd.h"
#include"main.h"
#include"AtomInfo.h"
#include"CoordSet.h"
#include"Util.h"

PyObject *P_globals = NULL;

PyObject *P_cmd = NULL;
PyObject *P_menu = NULL;
PyObject *P_xray = NULL;
PyObject *P_parser = NULL;
PyObject *P_setting = NULL;

PyObject *P_chempy = NULL;
PyObject *P_models = NULL;

PyObject *P_complete = NULL;

PyObject *P_exec = NULL;
PyObject *P_parse = NULL;
PyObject *P_lock = NULL;
PyObject *P_unlock = NULL;

PyObject *P_time = NULL;
PyObject *P_sleep = NULL;

unsigned int PyThread_get_thread_ident(void); /* critical functionality */

typedef struct {
  int id;
  PyThreadState *state;
} SavedThreadRec;

#define MAX_SAVED_THREAD 255

int NSavedThread = 0;
SavedThreadRec SavedThread[MAX_SAVED_THREAD];

int P_glut_thread_keep_out = 0; 
/* enables us to keep glut out if by chance it grabs the API
 * in the middle of a nested API based operation */

void PCatchInit(void);
void my_interrupt(int a);
char *getprogramname(void);

int PComplete(char *str,int buf_size)
{
  int ret = false;
  PyObject *result;
  char *st2;
  PBlockAndUnlockAPI();
  if(P_complete) {
    fflush(stdout);
    result = PyObject_CallFunction(P_complete,"s",str);
    if(result) {
      if(PyString_Check(result)) {
        st2 = PyString_AsString(result);
        UtilNCopy(str,st2,buf_size);
        ret=true;
      }
      Py_DECREF(result);
    }
  }
  PLockAPIAndUnblock();
  return(ret);
}

int PTruthCallStr(PyObject *object,char *method,char *argument)
{
  int result = false;
  PyObject *tmp;
  tmp = PyObject_CallMethod(object,method,"s",argument);
  if(tmp) {
    if(PyObject_IsTrue(tmp))
      result = 1;
    Py_DECREF(tmp);
  }
  return(result);
}
                                       
void PXDecRef(PyObject *obj)
{
  Py_XDECREF(obj);
}

void PSleep(int usec)
{ /* can only be called by the glut process */
#ifndef WIN32
  struct timeval tv;
  PUnlockAPIAsGlut();
  tv.tv_sec=0;
  tv.tv_usec=usec; 
  select(0,NULL,NULL,NULL,&tv);
  PLockAPIAsGlut();
#else
  PBlockAndUnlockAPI();
  PXDecRef(PyObject_CallFunction(P_sleep,"f",usec/1000000.0));
  PLockAPIAndUnblock();
#endif

}

static PyObject *PCatchWrite(PyObject *self, 	PyObject *args);

void my_interrupt(int a)
{
  exit(EXIT_FAILURE);
}

void PDumpTraceback(PyObject *err)
{
  PyObject *traceback;
  traceback = PyObject_GetAttrString(P_globals,"traceback");
  if(!traceback) 
    ErrFatal("PyMOL","can't find module 'traceback'");
  else {
    PyObject_CallMethod(traceback,"print_tb","o",err);
  }
}

int PAlterAtomState(float *v,char *expr,int read_only)
{
  PyObject *dict; /* TODO: this function badly need error checking code */
  int result=false;
  float f[3];
  PBlockAndUnlockAPI();

  dict = PyDict_New();

  PConvFloatToPyDictItem(dict,"x",v[0]);
  PConvFloatToPyDictItem(dict,"y",v[1]);
  PConvFloatToPyDictItem(dict,"z",v[2]);
  PyRun_String(expr,Py_single_input,P_globals,dict);
  if(PyErr_Occurred()) {
    PyErr_Print();
    result=false;
  } else if(!read_only) {
    f[0]=PyFloat_AsDouble(PyDict_GetItemString(dict,"x"));
    f[1]=PyFloat_AsDouble(PyDict_GetItemString(dict,"y"));
    f[2]=PyFloat_AsDouble(PyDict_GetItemString(dict,"z"));
    if(PyErr_Occurred()) {
      PyErr_Print();
      result=false;
      ErrMessage("AlterState","Aborting on error. Assignment may be incomplete.");
    } else {
      v[0]=f[0];
      v[1]=f[1];
      v[2]=f[2];
      result=true;
    }
  } else {
    result=true;
  }
  Py_DECREF(dict);
  PLockAPIAndUnblock();
  return result;
}

int PAlterAtom(AtomInfoType *at,char *expr,int read_only)
{
  AtomName name;
  ResName resn;
  ResIdent resi;
  Chain chain,alt;
  SegIdent segi;
  TextType textType;
  char atype[7];
  float b,q,partialCharge,vdw;
  int formalCharge,numericType;
  
  PyObject *dict;
  int result=false;
  
  if(at->hetatm)
    strcpy(atype,"HETATM");
  else
    strcpy(atype,"ATOM");
  PBlockAndUnlockAPI();

  dict = PyDict_New();

  PConvStringToPyDictItem(dict,"type",atype);
  PConvStringToPyDictItem(dict,"name",at->name);
  PConvStringToPyDictItem(dict,"resn",at->resn);
  PConvStringToPyDictItem(dict,"resi",at->resi);
  PConvStringToPyDictItem(dict,"chain",at->chain);
  PConvStringToPyDictItem(dict,"alt",at->alt);
  PConvStringToPyDictItem(dict,"segi",at->segi);
  PConvStringToPyDictItem(dict,"text_type",at->textType);
  if(at->customType!=cAtomInfoNoType)
    PConvIntToPyDictItem(dict,"numeric_type",at->customType);
  PConvFloatToPyDictItem(dict,"q",at->q);
  PConvFloatToPyDictItem(dict,"b",at->b);
  PConvFloatToPyDictItem(dict,"vdw",at->vdw);
  PConvFloatToPyDictItem(dict,"partial_charge",at->partialCharge);
  PConvIntToPyDictItem(dict,"formal_charge",at->formalCharge);
  
  PyRun_String(expr,Py_single_input,P_globals,dict);
  if(PyErr_Occurred()) {
    PyErr_Print();
    result=false;
  } else if(read_only) {
    result=true;
  } else {
    result=true;
    if(!PConvPyObjectToStrMaxLen(PyDict_GetItemString(dict,"type"),atype,6)) 
      result=false;
    else if(!PConvPyObjectToStrMaxLen(PyDict_GetItemString(dict,"name"),name,sizeof(AtomName)-1))
      result=false;
    else if(!PConvPyObjectToStrMaxLen(PyDict_GetItemString(dict,"resn"),resn,sizeof(ResName)-1))
      result=false;
    else if(!PConvPyObjectToStrMaxLen(PyDict_GetItemString(dict,"resi"),resi,sizeof(ResIdent)-1))
      result=false;
    else if(!PConvPyObjectToStrMaxLen(PyDict_GetItemString(dict,"segi"),segi,sizeof(SegIdent)-1))
      result=false;
    else if(!PConvPyObjectToStrMaxLen(PyDict_GetItemString(dict,"alt"),alt,sizeof(Chain)-1))
      result=false;
    else if(!PConvPyObjectToStrMaxLen(PyDict_GetItemString(dict,"chain"),chain,sizeof(Chain)-1))
      result=false;
    else if(!PConvPyObjectToStrMaxLen(PyDict_GetItemString(dict,"text_type"),textType,sizeof(TextType)-1))
      result=false;
    else if(!PConvPyObjectToFloat(PyDict_GetItemString(dict,"b"),&b))
      result=false;
    else if(!PConvPyObjectToFloat(PyDict_GetItemString(dict,"q"),&q))
      result=false;
    else if(!PConvPyObjectToFloat(PyDict_GetItemString(dict,"vdw"),&vdw))
      result=false;
    else if(!PConvPyObjectToFloat(PyDict_GetItemString(dict,"partial_charge"),&partialCharge))
      result=false;
    else if(!PConvPyObjectToInt(PyDict_GetItemString(dict,"formal_charge"),&formalCharge))
      result=false;
    if(!PConvPyObjectToInt(PyDict_GetItemString(dict,"numeric_type"),&numericType))
      numericType = cAtomInfoNoType;
    if(PyErr_Occurred()) {
      PyErr_Print();
      result=false;
    }
    if(result) { 
      at->hetatm=((atype[0]=='h')||(atype[0]=='H'));
      strcpy(at->name,name);
      strcpy(at->chain,chain);
      strcpy(at->resn,resn);
      if(strcmp(at->resi,resi)!=0)
        if(!sscanf(resi,"%i",&at->resv))
          at->resv=1;
      strcpy(at->resi,resi);
      strcpy(at->segi,segi);
      strcpy(at->chain,chain);
      strcpy(at->textType,textType);
      strcpy(at->alt,alt);
      if(numericType!=cAtomInfoNoType)
        at->customType = numericType;
      at->b=b;
      at->q=q;
      at->vdw=vdw;
      at->partialCharge=partialCharge;
      at->formalCharge=formalCharge;
    } else {
      ErrMessage("Alter","Aborting on error. Assignment may be incomplete.");
    }
  } 
  Py_DECREF(dict);
  PLockAPIAndUnblock();
  return(result);
}

int PLabelAtom(AtomInfoType *at,char *expr)
{
  PyObject *dict;
  int result;
  LabelType label;
  char atype[7];
  OrthoLineType buffer;
  if(at->hetatm)
    strcpy(atype,"HETATM");
  else
    strcpy(atype,"ATOM");
  PBlockAndUnlockAPI();

  dict = PyDict_New();

  PConvStringToPyDictItem(dict,"type",atype);
  PConvStringToPyDictItem(dict,"name",at->name);
  PConvStringToPyDictItem(dict,"resn",at->resn);
  PConvStringToPyDictItem(dict,"resi",at->resi);
  PConvStringToPyDictItem(dict,"chain",at->chain);
  PConvStringToPyDictItem(dict,"alt",at->alt);
  PConvStringToPyDictItem(dict,"segi",at->segi);
  PConvFloatToPyDictItem(dict,"vdw",at->vdw);
  PConvStringToPyDictItem(dict,"text_type",at->textType);
  PConvStringToPyDictItem(dict,"elem",at->elem);
  PConvIntToPyDictItem(dict,"geom",at->geom);
  PConvIntToPyDictItem(dict,"valence",at->valence);
  if(at->flags) {
    sprintf(buffer,"%X",at->flags);
    PConvStringToPyDictItem(dict,"flags",buffer);
  } else {
    PConvStringToPyDictItem(dict,"flags","0");
  }
  PConvFloatToPyDictItem(dict,"q",at->q);
  PConvFloatToPyDictItem(dict,"b",at->b);
  if(at->customType!=cAtomInfoNoType)
    PConvIntToPyDictItem(dict,"numeric_type",at->customType);
  PConvFloatToPyDictItem(dict,"partial_charge",at->partialCharge);
  PConvIntToPyDictItem(dict,"formal_charge",at->formalCharge);
  PConvIntToPyDictItem(dict,"id",at->id);
  PyRun_String(expr,Py_single_input,P_globals,dict);
  if(PyErr_Occurred()) {
    PyErr_Print();
    result=false;
  } else {
    result=true;
    if(!PConvPyObjectToStrMaxLen(PyDict_GetItemString(dict,"label"),label,sizeof(LabelType)-1))
      result=false;
    if(PyErr_Occurred()) {
      PyErr_Print();
      result=false;
    }
    if(result) { 
      strcpy(at->label,label);
    } else {
      ErrMessage("Label","Aborting on error. Labels may be incomplete.");
    }
  }
  Py_DECREF(dict);
  PLockAPIAndUnblock();
  return(result);
}

void PUnlockAPIAsGlut(void) /* must call with unblocked interpreter */
{
  PRINTFD(FB_Threads)
    " PUnlockAPIAsGlut-DEBUG: entered as thread 0x%x\n",PyThread_get_thread_ident()
    ENDFD;
  PBlock();
  PXDecRef(PyObject_CallFunction(P_unlock,NULL));
  PUnblock();
}

void PLockAPIAsGlut(void) /* must call with an unblocked interpreter */
{
  PRINTFD(FB_Threads)
    " PLockAPIAsGlut-DEBUG: entered as thread 0x%x\n",PyThread_get_thread_ident()
    ENDFD;

  PBlock();
  PXDecRef(PyObject_CallFunction(P_lock,NULL));
  while(P_glut_thread_keep_out) {
    /* IMPORTANT: keeps the glut thread out of an API operation... */
    /* NOTE: the keep_out variable can only be changed by the thread
       holding the API lock, therefore it is safe even through increment
       isn't atomic. */
    
    PXDecRef(PyObject_CallFunction(P_unlock,NULL));
#ifndef WIN32
    { 
      struct timeval tv;
      PUnblock();
      tv.tv_sec=0;
      tv.tv_usec=50000; 
      select(0,NULL,NULL,NULL,&tv);
      PBlock();
    } 
#else
    PXDecRef(PyObject_CallFunction(P_sleep,"f",0.050));
#endif
    PXDecRef(PyObject_CallFunction(P_lock,NULL));
  }
  PUnblock();
}

/* THESE CALLS ARE REQUIRED FOR MONOLITHIC COMPILATION TO SUCCEED UNDER WINDOWS. */

#ifdef _PYMOL_MONOLITHIC
void	initExtensionClass();
void	initsglite();
void    init_opengl();
void    init_opengl_num();
void    init_glu();
void    init_glu_num();
void    init_glut();
void    initopenglutil();
void    initopenglutil_num();
#endif

void PInitEmbedded(int argc,char **argv)
{
  /* This routine is called if we are running with an embedded Python interpreter */
  
  PyObject *args,*pymol,*invocation;

#ifdef WIN32
  OrthoLineType path_buffer,command;
  HKEY phkResult;
  int lpcbData;
  int lpType = REG_SZ;
  int r1,r2;
#endif

  Py_Initialize();
  PyEval_InitThreads();

  init_cmd();
#ifdef _PYMOL_MONOLITHIC
	initExtensionClass();
	initsglite();
    init_opengl();
    init_opengl_num();
    init_glu();
    init_glu_num();
    init_glut();
    initopenglutil();
	/* missing initopenglutil_num()???  WLD 3/26/2001*/
#endif

  PyRun_SimpleString("import os\n");
  PyRun_SimpleString("import sys\n");
#ifdef WIN32
  PyRun_SimpleString("if not os.environ.has_key('PYTHONPATH'): os.environ['PYTHONPATH']=''\n");

lpcbData = sizeof(OrthoLineType)-1;
r1=RegOpenKeyEx(HKEY_CLASSES_ROOT,"Software\\DeLano Scientific\\PyMOL\\PYMOL_PATH",0,KEY_EXECUTE,&phkResult);
  if(r1==ERROR_SUCCESS) {
	  r2 = RegQueryValueEx(phkResult,"",NULL,&lpType,path_buffer,&lpcbData);
	  if (r2==ERROR_SUCCESS) {
			/* use environment variable PYMOL_PATH first, registry entry second */
			sprintf(command,"_registry_pymol_path = r'''%s'''\n",path_buffer);
			PyRun_SimpleString(command);
			PyRun_SimpleString("print _registry_pymol_path\n");
			PyRun_SimpleString("if not os.environ.has_key('PYMOL_PATH'): os.environ['PYMOL_PATH']=_registry_pymol_path\n");
	  }
	RegCloseKey(phkResult);
	} 
/*  PyRun_SimpleString("if not os.environ.has_key('PYMOL_PATH'): os.environ['PYTHONPATH']=os.environ['PYTHONPATH']+';'+os.getcwd()+'/modules'\n");*/
  PyRun_SimpleString("if not os.environ.has_key('PYMOL_PATH'): os.environ['PYMOL_PATH']=os.getcwd()\n");
#endif
  PyRun_SimpleString("sys.path.append(os.environ['PYMOL_PATH']+'/modules')\n");
  PyRun_SimpleString("import pymol"); /* create the global PyMOL namespace */

  pymol = PyImport_AddModule("pymol"); /* get it */

  if(!pymol) ErrFatal("PyMOL","can't find module 'pymol'");

  invocation = PyObject_GetAttrString(pymol,"invocation"); /* get a handle to the invocation module */
  if(!pymol) ErrFatal("PyMOL","can't find module 'invocation'");

  args = PConvStringListToPyList(argc,argv); /* prepare our argument list */
  if(!pymol) ErrFatal("PyMOL","can't process arguments.");

  PXDecRef(PyObject_CallMethod(invocation,"parse_args","O",args)); /* parse the arguments */

}

void PGetOptions(int *pmgui,int *internal_gui,int *show_splash)
{
  PyObject *pymol,*invocation,*options;

  pymol = PyImport_AddModule("pymol"); /* get it */
  if(!pymol) ErrFatal("PyMOL","can't find module 'pymol'");

  invocation = PyObject_GetAttrString(pymol,"invocation"); /* get a handle to the invocation module */
  if(!pymol) ErrFatal("PyMOL","can't find module 'invocation'");

  options = PyObject_GetAttrString(invocation,"options");
  if(!pymol) ErrFatal("PyMOL","can't get 'invocation.options'.");

  (*pmgui) = ! PyInt_AsLong(PyObject_GetAttrString(options,"no_gui"));
  (*internal_gui) = PyInt_AsLong(PyObject_GetAttrString(options,"internal_gui"));
  (*show_splash) = PyInt_AsLong(PyObject_GetAttrString(options,"show_splash"));
  
}

void PRunString(char *str) /* runs a string in the global PyMOL module namespace */
{
  PXDecRef(PyObject_CallFunction(P_exec,"s",str));
}

void PInit(void) 
{
  PyObject *pymol,*sys,*pcatch;

  PCatchInit();   /* setup standard-output catch routine */

/* assumes that pymol module has been loaded */

  pymol = PyImport_AddModule("pymol"); /* get it */
  if(!pymol) ErrFatal("PyMOL","can't find module 'pymol'");
  P_globals = PyModule_GetDict(pymol);
  if(!P_globals) ErrFatal("PyMOL","can't find globals for 'pymol'");
  P_exec = PyDict_GetItemString(P_globals,"exec_str");
  if(!P_exec) ErrFatal("PyMOL","can't find 'pymol.exec_str()'");

  sys = PyDict_GetItemString(P_globals,"sys");
  if(!sys) ErrFatal("PyMOL","can't find 'pymol.sys'");
  pcatch = PyImport_AddModule("pcatch"); 
  if(!pcatch) ErrFatal("PyMOL","can't find module 'pcatch'");
  PyObject_SetAttrString(sys,"stdout",pcatch);
  PyObject_SetAttrString(sys,"stderr",pcatch);

  PRunString("import cmd\n");  
  P_cmd = PyDict_GetItemString(P_globals,"cmd");
  if(!P_cmd) ErrFatal("PyMOL","can't find 'cmd'");

  P_lock = PyObject_GetAttrString(P_cmd,"lock");
  if(!P_lock) ErrFatal("PyMOL","can't find 'pm.lock()'");

  P_unlock = PyObject_GetAttrString(P_cmd,"unlock");
  if(!P_unlock) ErrFatal("PyMOL","can't find 'pm.unlock()'");

  PRunString("import menu\n");  
  P_menu = PyDict_GetItemString(P_globals,"menu");
  if(!P_menu) ErrFatal("PyMOL","can't find module 'menu'");

  PRunString("import setting\n");  
  P_setting = PyDict_GetItemString(P_globals,"setting");
  if(!P_setting) ErrFatal("PyMOL","can't find module 'setting'");

#ifdef _PYMOL_XRAY
  PRunString("import xray\n");  
  P_xray = PyDict_GetItemString(P_globals,"xray");
  if(!P_xray) ErrFatal("PyMOL","can't find module 'xray'");
#endif

#ifdef WIN32
  PRunString("import time\n");  
  P_time = PyDict_GetItemString(P_globals,"time");
  if(!P_time) ErrFatal("PyMOL","can't find module 'time'");

  P_sleep = PyObject_GetAttrString(P_time,"sleep");
  if(!P_sleep) ErrFatal("PyMOL","can't find 'time.sleep()'");
#endif

  PRunString("import parser\n");  
  P_parser = PyDict_GetItemString(P_globals,"parser");
  if(!P_parser) ErrFatal("PyMOL","can't find module 'parser'");

  P_parse = PyObject_GetAttrString(P_parser,"parse");
  if(!P_parse) ErrFatal("PyMOL","can't find 'parser.parse()'");

  P_complete = PyObject_GetAttrString(P_parser,"complete");
  if(!P_complete) ErrFatal("PyMOL","can't find 'parser.complete()'");

  PRunString("import chempy"); 
  P_chempy = PyDict_GetItemString(P_globals,"chempy");
  if(!P_chempy) ErrFatal("PyMOL","can't find 'chempy'");

  PRunString("from chempy import models"); 
  P_models = PyDict_GetItemString(P_globals,"models");
  if(!P_models) ErrFatal("PyMOL","can't find 'chempy.models'");

  PRunString("import util\n");  
#ifdef _PYMOL_XRAY
  PRunString("import sglite\n"); 
#endif
  PRunString("import string\n"); 
  PRunString("import traceback\n"); 

  /* backwards compatibility */

  PRunString("pm = cmd\n");  
  PRunString("pmu = util\n");  

  PRunString("glutThread = thread.get_ident()");

#ifndef WIN32
  signal(SIGINT,my_interrupt);
#endif

}

void PStereoOff(void) 
{
  PBlock();
  PRunString("pm._stereo(0)");
  PUnblock();
}

void PFree(void)
{
}

void PExit(int code)
{
  PBlock();
  MainFree();
  Py_Exit(code);
}

void PParse(char *str) 
{
  OrthoCommandIn(str);
}

void PFlush(void) {  
  /* NOTE: ASSUMES we current have unblocked Python threads and a locked API */
  char buffer[OrthoLineLength+1];
  while(OrthoCommandOut(buffer)) {
    PBlockAndUnlockAPI();
    PXDecRef(PyObject_CallFunction(P_parse,"s",buffer));
    PLockAPIAndUnblock();
  }
}

void PFlushFast(void) {
  /* NOTE: ASSUMES we currently have blocked Python threads and an unlocked API */ 
 char buffer[OrthoLineLength+1];
  while(OrthoCommandOut(buffer)) {
   PXDecRef(PyObject_CallFunction(P_parse,"s",buffer));
  }
}

void PBlock(void)
{
  int a,id;
  /* synchronize python */

  PRINTFD(FB_Threads)
    " PBlock-DEBUG: entered as thread 0x%x\n",PyThread_get_thread_ident()
    ENDFD;

  id = PyThread_get_thread_ident();
  a = NSavedThread;
  while(a) {
    a--;
    if(SavedThread[a].id==id) {
      PyEval_RestoreThread(SavedThread[a].state);
      NSavedThread--;
      if(NSavedThread) {
        SavedThread[a] = SavedThread[NSavedThread];
      }
      return;
    }
  }
  ErrFatal("PBlock","oops, no saved thread -- your threads must be tangled!");
}

int PAutoBlock(void)
{
  int a,id;
  /* synchronize python */

  id = PyThread_get_thread_ident();
  a = NSavedThread;
  while(a) {
    a--;
    if(SavedThread[a].id==id) {
      PyEval_RestoreThread(SavedThread[a].state);
      NSavedThread--;
      if(NSavedThread) {
        SavedThread[a] = SavedThread[NSavedThread];
      }
      PRINTFD(FB_Threads)
        " PAutoBlock-DEBUG: saved thread 0x%x\n",PyThread_get_thread_ident()
        ENDFD;
      return 1;
    }
  }
  return 0;
}

void PAutoUnblock(int flag)
{
  if(flag) PUnblock();
}

void PBlockAndUnlockAPI(void)
{
  PBlock();
  PXDecRef(PyObject_CallFunction(P_unlock,NULL));
}

void PLockAPIAndUnblock(void)
{
  PXDecRef(PyObject_CallFunction(P_lock,NULL));
  PUnblock();
}

void PUnblock(void)
{
  /* NOTE: ASSUMES a locked API */

  PRINTFD(FB_Threads)
    " PUnblock-DEBUG: entered as thread 0x%x\n",PyThread_get_thread_ident()
    ENDFD;

  SavedThread[NSavedThread].id = PyThread_get_thread_ident();
  SavedThread[NSavedThread].state = PyEval_SaveThread();  
  NSavedThread++;
}

void PDefineFloat(char *name,float value) {
  char buffer[OrthoLineLength];
  sprintf(buffer,"%s = %f\n",name,value);
  PBlock();
  PRunString(buffer);
  PUnblock();
}

/* This function is called by the interpreter to get its own name */
char *getprogramname(void)
{
	return("PyMOL");
}

/* A static module */

static PyObject *PCatchWrite(PyObject *self, 	PyObject *args)
{
  char *str;
  
  PyArg_ParseTuple(args,"s",&str);
  if(str[0]) {
	 OrthoAddOutput(str);
  }
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject *PCatchFlush(PyObject *self, 	PyObject *args)
{
  fflush(stdout);
  fflush(stderr);
  Py_INCREF(Py_None);
  return Py_None;
}

static PyMethodDef PCatch_methods[] = {
	{"write",	  PCatchWrite,   METH_VARARGS},
	{"flush",	  PCatchFlush,   METH_VARARGS},
	{NULL,		NULL}		/* sentinel */
};

void PCatchInit(void)
{
	PyImport_AddModule("pcatch");
	Py_InitModule("pcatch", PCatch_methods);
}


