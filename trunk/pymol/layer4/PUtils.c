
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

PyObject *P_pm = NULL;
PyObject *P_pmm = NULL;
PyObject *P_pmx = NULL;

PyThreadState *P_glut_thread_state; /* this is the state for the main GUI thread */
PyThreadState *P_api_thread_state; /* this is the thread state for an alternate thread */
int P_glut_thread_active = 1;
int P_glut_thread_keep_out = 0; /* enables us to keep glut out if by chance it grabs the API
                                        * in the middle of a nested API based operation */

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


int PAlterAtom(AtomInfoType *at,char *expr)
{
  char atype[255],name[255],resi[255],chain[255],resn[255],segi[255];
  float b,q;
  PyObject *output;
  int result;

  if(at->hetatm)
    strcpy(atype,"HETATM");
  else
    strcpy(atype,"ATOM");
  PBlockAndUnlockAPI();
  output = PyObject_CallMethod(P_pm,"_alter_do","[ssssssfffffs]",
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
    if(output) {
      result = 0;
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
      Py_DECREF(output);
    } else {
      result = -1;
    }
    PLockAPIAndUnblock();
    return(result);
}

void PUnlockAPIAsGlut(void)
{
  PyEval_RestoreThread(P_glut_thread_state); /* grab python */
  PyRun_SimpleString("pm.unlock()\n");
  P_glut_thread_state = PyEval_SaveThread(); /* release python */
}

void PLockAPIAsGlut(void)
{
  PyEval_RestoreThread(P_glut_thread_state); /* grab python */
  PyRun_SimpleString("pm.lock()");
  while(P_glut_thread_keep_out) { /* IMPORTANT: keeps the glut thread out of an API operation... */
	 PyRun_SimpleString("pm.unlock()");
    P_glut_thread_state = PyEval_SaveThread(); /* release python */
    PSleep(50000); /* wait 50 msec */
    PyEval_RestoreThread(P_glut_thread_state); /* grab python */
	 PyRun_SimpleString("pm.lock()");    
  }
  P_glut_thread_state = PyEval_SaveThread(); /* release python */
  P_glut_thread_active = 1; /* if we come in on a glut event - then it is the active thread */
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
  P_pm = PyDict_GetItemString(PM_Globals,"pm");

  PyRun_SimpleString("import pmm\n");  
  P_pmm = PyDict_GetItemString(PM_Globals,"pmm");

  PyRun_SimpleString("import pmx\n");  
  P_pmx = PyDict_GetItemString(PM_Globals,"pmx");

  PyRun_SimpleString("import pmu\n");  
  PyRun_SimpleString("import sglite\n"); 
  PyRun_SimpleString("import string\n"); 
  

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

void PFlush(void) {  
  /* NOTE: ASSUMES we current have unblocked Python threads and a locked API */
  char buffer[OrthoLineLength+1];
  if(OrthoCommandOut(buffer)) {

   PBlockAndUnlockAPI();

	PyDict_SetItemString(PM_Globals,"pymol_cmd",PyString_FromString(buffer));
	PyRun_SimpleString("pmp_cmd[pmp_nest] = pymol_cmd");
	PyRun_SimpleString("exec(pymol,globals(),globals())");

   PLockAPIAndUnblock();
  }
}

void PBlock(void)
{
  /* synchronize python */

  if(P_glut_thread_active)
    PyEval_RestoreThread(P_glut_thread_state);
  else 
    PyEval_RestoreThread(P_api_thread_state);

}

void PBlockAndUnlockAPI(void)
{
  PBlock();
  PyRun_SimpleString("pm.unlock()");
}

void PLockAPIAndUnblock(void)
{
  PyRun_SimpleString("pm.lock()");
  PUnblock();
}

void PUnblock(void)
{
  /* NOTE: ASSUMES a locked API */

  PyObject *is_glut;

  /* P_glut_thread_active will not change as long as lock holds */

  is_glut = PyObject_CallMethod(P_pm,"is_glut_thread","");
  P_glut_thread_active = PyInt_AsLong(is_glut);
  Py_DECREF(is_glut);

  /* allow python to run async */

  if(P_glut_thread_active) 
    P_glut_thread_state = PyEval_SaveThread();
  else
    P_api_thread_state = PyEval_SaveThread();
}

void PDefineFloat(char *name,float value) {
  char buffer[OrthoLineLength];
  sprintf(buffer,"%s = %f\n",name,value);
  PBlock();
  PyRun_SimpleString(buffer);
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

static PyMethodDef PCatch_methods[] = {
	{"write",	  PCatchWrite,   METH_VARARGS},
	{NULL,		NULL}		/* sentinel */
};


void PCatchInit(void)
{
	PyImport_AddModule("pcatch");
	Py_InitModule("pcatch", PCatch_methods);
}


