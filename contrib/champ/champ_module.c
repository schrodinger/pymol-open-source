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

#include"os_std.h"
#include"os_memory.h"
#include"os_python.h"
#include"const.h"
#include"champ.h"
#include"list.h"
#include"vla.h"

static PyObject *RetNone()
{
  return(Py_BuildValue("(iO)",0,Py_None));
}

static PyObject *RetObj(int ok,PyObject *result)
{
  PyObject *ret;
  if(result==Py_None) {
    Py_INCREF(result);
  } else if(result==NULL) {
    result=Py_None;
    Py_INCREF(result);
  } 
  ret = Py_BuildValue("(iO)",!ok,result);
  Py_DECREF(result);
  return(ret);
}

static PyObject *RetInt(int ok,int result) /* status/integer return */
{
  return(Py_BuildValue("(ii)",!ok,result));
}

static PyObject *_new(PyObject *self,      PyObject *args)
{
  return(PyCObject_FromVoidPtr((void*)ChampNew(),(void*)(void*)ChampFree));
}

static PyObject *_memory_dump(PyObject *self,      PyObject *args)
{
  os_memory_dump();
  return(RetNone());
}

static PyObject *insert_smiles(PyObject *self,      PyObject *args)
{
  int ok=true;
  int result = 0;
  char *str1;
  PyObject *O;
  CChamp *I;
  ok = PyArg_ParseTuple(args,"Os",&O,&str1);
  ok = PyCObject_Check(O);
  if(ok) {
    I = PyCObject_AsVoidPtr(O);
    result = ChampSmiToPat(I,str1);
  }
  return(RetInt(ok,result));
}

static PyObject *get_pattern_list(PyObject *self,      PyObject *args)
{
  int ok=true;
  PyObject *result = NULL;
  PyObject *O;
  int list_index;
  CChamp *I;
  int i,c;
  ok = PyArg_ParseTuple(args,"Oi",&O,&list_index);
  ok = PyCObject_Check(O);
  if(ok) {
    I = PyCObject_AsVoidPtr(O);
    i = list_index;
    c = 0;
    while(i) {
      i = I->Int[i].link;
      c++;
    }
    result = PyList_New(c);
    i = list_index;
    c = 0;
    while(i) {
      PyList_SetItem(result,c,PyInt_FromLong(I->Int[i].value));
      i = I->Int[i].link;
      c++;
    }
  }
  return(RetObj(ok,result));
}


static PyObject *get_smiles(PyObject *self,      PyObject *args)
{
  int ok=true;
  PyObject *result = NULL;
  PyObject *O;
  int pat_index;
  CChamp *I;
  char *smi;

  ok = PyArg_ParseTuple(args,"Oi",&O,&pat_index);
  ok = PyCObject_Check(O);
  if(ok) {
    I = PyCObject_AsVoidPtr(O);
    smi = ChampPatToSmiVLA(I,pat_index);
    result = PyString_FromString(smi);
    vla_free(smi);
  }
  return(RetObj(ok,result));
}


static PyObject *extend_list_smiles(PyObject *self,      PyObject *args)
{
  int ok=true;
  int result = 0;
  int list_index,pat_index;
  int a,l;
  PyObject *O,*list;
  CChamp *I;
  ok = PyArg_ParseTuple(args,"OiO",&O,&list_index,&list);
  ok = PyCObject_Check(O);
  ok = PyList_Check(list);
  if(ok) {
    I = PyCObject_AsVoidPtr(O);
    l = PyList_Size(list);
    for(a=0;a<l;a++) {
      pat_index = ChampSmiToPat(I,PyString_AsString(PyList_GetItem(list,a)));
      if(!pat_index) {
        ok=false;
        break;
      }
      list_index = ListElemPush(&I->Int,list_index);
      I->Int[list_index].value = pat_index;
    }
  }
  if(ok) result=list_index;
  return(RetInt(ok,result));
}


static PyObject *sss_1v1_b(PyObject *self,      PyObject *args)
{
  int ok=true;
  int result = 0;
  int int1,int2;

  PyObject *O;
  CChamp *I;
  ok = PyArg_ParseTuple(args,"Oii",&O,&int1,&int2);
  ok = PyCObject_Check(O);
  if(ok) {
    I = PyCObject_AsVoidPtr(O);
    result = ChampSSS_1V1_B(I,int1,int2);
  }
  return(RetInt(ok,result));
}

static PyObject *sss_1vN_n(PyObject *self,      PyObject *args)
{
  int ok=true;
  int result = 0;
  int int1,int2;

  PyObject *O;
  CChamp *I;
  ok = PyArg_ParseTuple(args,"Oii",&O,&int1,&int2);
  ok = PyCObject_Check(O);
  if(ok) {
    I = PyCObject_AsVoidPtr(O);
    result = ChampSSS_1VN_N(I,int1,int2);
  }
  return(RetInt(ok,result));
}

static PyObject *insert_model(PyObject *self,      PyObject *args)
{
  int ok=true;
  int result = 0;

  PyObject *O,*M;
  CChamp *I;
  ok = PyArg_ParseTuple(args,"OO",&O,&M);
  ok = PyCObject_Check(O);
  if(ok) {
    I = PyCObject_AsVoidPtr(O);
    result = ChampModelToPat(I,M);
  }
  return(RetInt(ok,result));
}

static PyMethodDef champ_methods[] = {
  {"new",	                 _new,                   METH_VARARGS },
  {"memory_dump",            _memory_dump,           METH_VARARGS },
  {"get_smiles",             get_smiles,             METH_VARARGS },
  {"insert_smiles",          insert_smiles,          METH_VARARGS },
  {"insert_model",           insert_model,           METH_VARARGS },
  {"extend_list_smiles",     extend_list_smiles,     METH_VARARGS },
  {"get_pattern_list",       get_pattern_list,       METH_VARARGS },
  {"sss_1v1_b",              sss_1v1_b,              METH_VARARGS },
  {"sss_1vN_n",              sss_1vN_n,              METH_VARARGS },
  {NULL,		                 NULL}     /* sentinel */        
};

void init_champ(void);
void init_champ(void)
{
  Py_InitModule("_champ", champ_methods);
}

