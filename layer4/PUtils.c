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
/* Example of embedding Python in another program */

#include<stdlib.h>
#include<Python.h>
#include<signal.h>
#include <sys/types.h>
#include <sys/time.h>
#include <unistd.h>

#include"Base.h"
#include"PUtils.h"
#include"Ortho.h"
#include"PM.h"
#include"main.h"

void PSleep(int usec)
{ /* assumes threads have already been unblocked */
  struct timeval tv;

  PyThreadState *_save;
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
#endif

#ifdef _PYMOL_THREADS
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
  PyRun_SimpleString("import string\n"); 

#ifndef _PYMOL_MODULE
  PyRun_SimpleString("import thread\n"); 
  PyRun_SimpleString("import threading\n"); 
  PyRun_SimpleString("lock_oq = threading.RLock()");
  PyRun_SimpleString("lock_iq = threading.RLock()");
  PyRun_SimpleString("lock_api = threading.RLock()");
  PyRun_SimpleString("import thread\n"); 
  PyRun_SimpleString("glutThread = thread.get_ident()");
#endif


  PyRun_SimpleString("from pmp import *\n");

#ifndef _PYMOL_MODULE
  PyRun_SimpleString("pm.setup_global_locks()");
  /*  PyRun_SimpleString("execfile(os.environ['PYMOL_PATH']+'/modules/pmx.py',globals(),globals())\n"); */
  PyRun_SimpleString("sys.argv=['pymol']\n");
  /*  PyRun_SimpleString("thread.start_new_thread(execfile,(os.environ['PYMOL_PATH']+'/modules/pmg.py',globals(),locals()))\n"); */
  PyRun_SimpleString("_t=threading.Thread(target=execfile,args=(os.environ['PYMOL_PATH']+'/modules/pmg.py',globals(),locals()))\n_t.setDaemon(1)\n_t.start()"); 
#endif

  signal(SIGINT,my_interrupt);
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






