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

#include"PUtils.h"
#include"Ortho.h"
#include"PM.h"
#include"main.h"

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

void PInit(void)
{
  /* Initialize the Python interpreter.  Required. */
  Py_Initialize();

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
  PyRun_SimpleString("import pmx\n"); /* the API */
  PyRun_SimpleString("import pcatch\n");
  
  PyRun_SimpleString("sys.stdout = pcatch\n");
  PyRun_SimpleString("sys.stderr = pcatch\n");
  /*	PyRun_SimpleString("execfile(\"modules/pymol.py\",globals(),globals());");*/

  PyRun_SimpleString("pmx.set_globals(globals())");
  
  PyRun_SimpleString("import pm\n"); 
  PyRun_SimpleString("import pmu\n");
  
  PyRun_SimpleString("import string\n"); 
  PyRun_SimpleString("from pmp import *\n");
  
  PyRun_SimpleString("import pmg\n"); 
  signal(SIGINT,my_interrupt);
}

void PExit(int code)
{
  Py_Finalize();
  MainFree();
  Py_Exit(code);
}

void PParse(char *str) {

  PyDict_SetItemString(PM_Globals,"pymol_cmd",PyString_FromString(str));

  PyRun_SimpleString("pmp_cmd[pmp_nest] = pymol_cmd");
  PyRun_SimpleString("exec(pymol,globals(),globals())");

}

void PRunStr(char *str)
{
  PyRun_SimpleString(str);
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






