
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
#include<signal.h>
#include<string.h>
#include<sys/types.h>
#include<sys/time.h>
#include<unistd.h>

#include"MemoryDebug.h"
#include"Base.h"
#include"PUtils.h"
#include"Ortho.h"
#include"PM.h"
#include"main.h"
#include"AtomInfo.h"
#include"CoordSet.h"


void PSleep(int usec)
{ /* assumes threads have already been unblocked */
  struct timeval tv;

  tv.tv_sec=0;
  tv.tv_usec=usec; 
  select(0,NULL,NULL,NULL,&tv);
}

static PyObject *PCatchWrite(PyObject *self, 	PyObject *args);
void PCatchInit(void);
void PMInit(void);
void PMQuit(void);
void my_interrupt(int a);

char *getprogramname(void);

void my_interrupt(int a)
{
  exit(EXIT_FAILURE);
}

PyObject *PFloatVLAToPyList(float *f)
{
  int a,l;
  PyObject *result = Py_None;
  l=VLAGetSize(f);
  result=PyList_New(l);
  for(a=0;a<l;a++) {
    PyList_SetItem(result,a,PyFloat_FromDouble((double)f[a]));
  }
  return(result);
}

void PLock(int lock,PyThreadState **save)
{
  PyThreadState *_save;

  _save = (*save);
  Py_BLOCK_THREADS;
  switch(lock) {
  case cLockAPI:
	 PyRun_SimpleString("pm.lock()");
	 break;
  }
  Py_UNBLOCK_THREADS;
  (*save)=_save;
}

int PAlterAtom(AtomInfoType *at,char *expr)
{
  char atype[255],name[255],resi[255],chain[255],resn[255],segi[255];
  float b,q;
  PyObject *input,*output;
  int result;
  if(at->hetatm)
    strcpy(atype,"HETATM");
  else
    strcpy(atype,"ATOM");
  input = Py_BuildValue("[ssssssfffffs]",
                        expr,
                        atype,
                        at->name,
                        at->resn,
                        at->chain,
                        at->resi,
                        0.0,0.0,0.0,
                        at->q,
                        at->b,
                        at->segi);
  PyDict_SetItemString(PM_Globals,"alter_inp",input);
  Py_DECREF(input);
  result = PyRun_SimpleString("alter_out = pm._alter_do(alter_inp)");
  output = PyDict_GetItemString(PM_Globals,"alter_out");
  if(output) {
    strcpy(atype,PyString_AsString(PyList_GetItem(output,0)));
    strcpy(name,PyString_AsString(PyList_GetItem(output,1)));
    strcpy(resn,PyString_AsString(PyList_GetItem(output,2)));
    strcpy(chain,PyString_AsString(PyList_GetItem(output,3)));
    strcpy(resi,PyString_AsString(PyList_GetItem(output,4)));
    q=PyFloat_AsDouble(PyList_GetItem(output,8));
    b=PyFloat_AsDouble(PyList_GetItem(output,9));
    strcpy(segi,PyString_AsString(PyList_GetItem(output,10)));
    strcpy(at->name,name);
    strcpy(at->resi,resi);
    strcpy(at->chain,chain);
    strcpy(at->resn,resn);
    at->b = b;
    at->q = q;
    strcpy(at->segi,segi);
    at->hetatm = (strcmp(atype,"HETATM")==0);
    /*    printf("%s %s %s %s %s %8.3f %8.3f %s\n",
          atype,name,resi,chain,resn,q,b,segi);*/
  }
  return(result);
}

void PUnlock(int lock,PyThreadState **save)
{
  PyThreadState *_save;
  _save = (*save);
  Py_BLOCK_THREADS;
  switch(lock) {
  case cLockAPI:
	 PyRun_SimpleString("pm.unlock()\n");
	 break;
  }
  Py_UNBLOCK_THREADS;
  (*save)=_save;
}

void PInit(void)
{
  /* Initialize the Python interpreter.  Required. */

#ifndef _PYMOL_MODULE
  Py_Initialize();
  PyEval_InitThreads();
#endif

  /* Add a static module */
  PCatchInit();
  PMInit();

  PyRun_SimpleString("import os\n");
  PyRun_SimpleString("import sys\n");
  PyRun_SimpleString("sys.path.append(os.environ['PYMOL_PATH']+'/modules')\n");
  /* redirect output to our catch routine*/
  PyRun_SimpleString("import _pm\n"); /* the API */
  PyRun_SimpleString("import pcatch\n");
  
   PyRun_SimpleString("sys.stdout = pcatch\n");
   /*		PyRun_SimpleString("sys.stderr = pcatch\n");*/

  PyRun_SimpleString("_pm.set_globals(globals())");
  
  PyRun_SimpleString("import pm\n"); 
  PyRun_SimpleString("import pmu\n");  
  PyRun_SimpleString("import pmm\n");  
  PyRun_SimpleString("import string\n"); 
  PyRun_SimpleString("import sglite\n"); 

#ifndef _PYMOL_MODULE
  PyRun_SimpleString("import thread\n"); 
  PyRun_SimpleString("import threading\n"); 
  PyRun_SimpleString("lock_oq = threading.RLock()");
  PyRun_SimpleString("lock_iq = threading.RLock()");
  PyRun_SimpleString("lock_api = threading.RLock()");
  PyRun_SimpleString("glutThread = thread.get_ident()");
#endif

  PyRun_SimpleString("from pmp import *\n");

#ifndef _PYMOL_MODULE
  PyRun_SimpleString("pm.setup_global_locks()");
  if(PMGUI) {
    PyRun_SimpleString("sys.argv=['pymol']\n");
    PyRun_SimpleString("_t=threading.Thread(target=execfile,args=(os.environ['PYMOL_PATH']+'/modules/pmg.py',globals(),locals()))\n_t.setDaemon(1)\n_t.start()"); 
  }
#endif

  signal(SIGINT,my_interrupt);
}

void PStereoOff(void) 
{
  PyRun_SimpleString("pm._stereo(0)");
}

void PFree(void)
{
}

void PExit(int code)
{
  Py_Finalize();
  MainFree();
  Py_Exit(code);
}

void PParse(char *str) 
{
  OrthoCommandIn(str);
}

void PFlush(PyThreadState **save) {  
  /* NOTE: ASSUMES we current have unblocked Python threads and a locked API */
  char buffer[OrthoLineLength+1];
  if(OrthoCommandOut(buffer)) {

	PyThreadState *_save;
	_save = (*save);

	PUnlock(cLockAPI,&_save);
	Py_BLOCK_THREADS;

	PyDict_SetItemString(PM_Globals,"pymol_cmd",PyString_FromString(buffer));
	PyRun_SimpleString("pmp_cmd[pmp_nest] = pymol_cmd");
	PyRun_SimpleString("exec(pymol,globals(),globals())");

	Py_UNBLOCK_THREADS;
	PLock(cLockAPI,&_save);

  }
}

void PBlock(PyThreadState **save)
{
  /* NOTE: ASSUMES we current have unblocked Python threads and a locked API */

  PyThreadState *_save;
  _save = (*save);
  Py_BLOCK_THREADS;
  PyRun_SimpleString("pm.unlock()");
  (*save)=_save;
}

void PUnblock(PyThreadState **save)
{
  /* NOTE: ASSUMES we current have blocked Python threads and an unlocked API */
  PyThreadState *_save;
  _save = (*save);
  PyRun_SimpleString("pm.lock()");
  Py_UNBLOCK_THREADS;
  (*save)=_save;
}

void PDefineFloat(char *name,float value) {
  char buffer[OrthoLineLength];
  sprintf(buffer,"%s = %f\n",name,value);
  PyRun_SimpleString(buffer);
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

static PyMethodDef PCatch_methods[] = {
	{"write",	  PCatchWrite,   METH_VARARGS},
	{NULL,		NULL}		/* sentinel */
};


void PCatchInit(void)
{
	PyImport_AddModule("pcatch");
	Py_InitModule("pcatch", PCatch_methods);
}






