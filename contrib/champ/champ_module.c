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

static PyObject *RetStatus(int ok)
{
  return(Py_BuildValue("(iO)",!ok,Py_None));
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

static PyObject *list_new(PyObject *self,      PyObject *args)
{
  int ok=true;
  int result = 0;
  PyObject *O;
  CChamp *I;
  ok = PyArg_ParseTuple(args,"O",&O);
  ok = PyCObject_Check(O);
  if(ok) {
    I = PyCObject_AsVoidPtr(O);
    result = ListElemNewZero(&I->Int); 
  }
  return(RetInt(ok,result));
}

static PyObject *list_free(PyObject *self,      PyObject *args)
{
  int ok=true;
  int list_index,list_handle,purge;
  int i;
  PyObject *O;
  CChamp *I;
  ok = PyArg_ParseTuple(args,"Oii",&O,&list_handle,&purge);
  ok = PyCObject_Check(O);
  if(ok) {
    I = PyCObject_AsVoidPtr(O);
    list_index = I->Int[list_handle].link;
    i = list_index;
    while(i) {
      if(purge) ChampPatFree(I,I->Int[i].value);
      i = I->Int[i].link;
    }
    ListElemFreeChain(I->Int,list_index);
  }
  return(RetStatus(ok));
}

static PyObject *list_get_pattern_list(PyObject *self,      PyObject *args)
{
  int ok=true;
  PyObject *result = NULL;
  PyObject *O;
  int list_index,list_handle;
  CChamp *I;
  int i,c;
  ok = PyArg_ParseTuple(args,"Oi",&O,&list_handle);
  ok = PyCObject_Check(O);
  if(ok) {
    I = PyCObject_AsVoidPtr(O);
    list_index = I->Int[list_handle].link;
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
    smi = ChampPatToSmiVLA(I,pat_index,NULL);
    result = PyString_FromString(smi);
    vla_free(smi);
  }
  return(RetObj(ok,result));
}


static PyObject *list_get_smiles_list(PyObject *self,      PyObject *args)
{
  int ok=true;
  PyObject *result = NULL,*str;
  PyObject *O;
  int list_index,list_handle;
  CChamp *I;
  int i,c;
  char *smi = NULL;
  ok = PyArg_ParseTuple(args,"Oi",&O,&list_handle);
  ok = PyCObject_Check(O);
  if(ok) {
    I = PyCObject_AsVoidPtr(O);
    list_index = I->Int[list_handle].link;
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
      smi = ChampPatToSmiVLA(I,I->Int[i].value,smi);
      str = PyString_FromString(smi);
      PyList_SetItem(result,c,str);
      i = I->Int[i].link;
      c++;
    }
    vla_free(smi);
  }
  return(RetObj(ok,result));
}


static PyObject *list_prepend_smiles_list(PyObject *self,      PyObject *args)
{
  int ok=true;
  int list_index,list_handle,pat_index;
  int a,l;
  PyObject *O,*smi_list;
  CChamp *I;
  ok = PyArg_ParseTuple(args,"OiO",&O,&list_handle,&smi_list);
  ok = PyCObject_Check(O);
  ok = PyList_Check(smi_list);
  if(ok) {
    I = PyCObject_AsVoidPtr(O);
    l = PyList_Size(smi_list);
    for(a=l-1;a>=0;a--) {
      pat_index = ChampSmiToPat(I,PyString_AsString(PyList_GetItem(smi_list,a)));
      if(!pat_index) {
        ok=false;
        break;
      }
      list_index = ListElemPush(&I->Int,I->Int[list_handle].link);
      I->Int[list_index].value = pat_index;
      I->Int[list_handle].link = list_index;
    }
  }
  return(RetStatus(ok));
}


static PyObject *match_1v1_b(PyObject *self,      PyObject *args)
{
  int ok=true;
  int result = 0;
  int pat1,pat2;

  PyObject *O;
  CChamp *I;
  ok = PyArg_ParseTuple(args,"Oii",&O,&pat1,&pat2);
  ok = PyCObject_Check(O);
  if(ok) {
    I = PyCObject_AsVoidPtr(O);
    result = ChampMatch_1V1_B(I,pat1,pat2);
  }
  return(RetInt(ok,result));
}


static PyObject *match_1v1_map(PyObject *self,      PyObject *args)
{
  int ok=true;
  PyObject *result = NULL;
  PyObject *set,*pairs,*tmpl,*targ;
  int pat1,pat2;
  int limit;
  int mat_start;
  int mat_cnt,pair_cnt;
  int a,b,i,j;
  ListMatch *mat;

  PyObject *O;
  CChamp *I;
  ok = PyArg_ParseTuple(args,"Oiii",&O,&pat1,&pat2,&limit);
  ok = PyCObject_Check(O);
  if(ok) {
    I = PyCObject_AsVoidPtr(O);
    mat_start = ChampMatch_1V1_Map(I,pat1,pat2,limit);
    mat_cnt = 0;
    i = mat_start;
    while(i) {
      mat_cnt++;
      i = I->Match[i].link;
    }
    result = PyList_New(mat_cnt);
    i = mat_start;
    for(a=0;a<mat_cnt;a++) {
      mat = I->Match + i;
      set = PyList_New(2); /* atom matches, bond matches */

      pairs = PyList_New(2); /* paired atoms */
      pair_cnt = 0;
      j = mat->atom;
      while(j) {
        pair_cnt++;
        j = I->Int2[j].link;
      }
      tmpl = PyList_New(pair_cnt);
      targ = PyList_New(pair_cnt);
      j = mat->atom;
      for(b=0;b<pair_cnt;b++) {
        PyList_SetItem(tmpl,b,PyInt_FromLong(I->Atom[I->Int2[j].value[0]].index));
        PyList_SetItem(targ,b,PyInt_FromLong(I->Atom[I->Int2[j].value[1]].index));
        j = I->Int2[j].link;
      }
      PyList_SetItem(pairs,0,tmpl);
      PyList_SetItem(pairs,1,targ);
      PyList_SetItem(set,0,pairs);

      pairs = PyList_New(2); /* paired bonds */
      pair_cnt = 0;
      j = mat->bond;
      while(j) {
        pair_cnt++;
        j = I->Int2[j].link;
      }
      tmpl = PyList_New(pair_cnt);
      targ = PyList_New(pair_cnt);
      j = mat->bond;
      for(b=0;b<pair_cnt;b++) {
        PyList_SetItem(tmpl,b,PyInt_FromLong(I->Bond[I->Int2[j].value[0]].index));
        PyList_SetItem(targ,b,PyInt_FromLong(I->Bond[I->Int2[j].value[1]].index));
        j = I->Int2[j].link;
      }
      PyList_SetItem(pairs,0,tmpl);
      PyList_SetItem(pairs,1,targ);
      PyList_SetItem(set,1,pairs);

      PyList_SetItem(result,a,set);
      i = I->Match[i].link;      
    }
  }
  return(RetObj(ok,result));
}

static PyObject *match_1vN_n(PyObject *self,      PyObject *args)
{
  int ok=true;
  int result = 0;
  int list_handle,list_index;
  int pattern;

  PyObject *O;
  CChamp *I;
  ok = PyArg_ParseTuple(args,"Oii",&O,&pattern,&list_handle);
  ok = PyCObject_Check(O);
  if(ok) {
    I = PyCObject_AsVoidPtr(O);
    list_index = I->Int[list_handle].link;
    result = ChampMatch_1VN_N(I,pattern,list_index);
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
  {"new",	                     _new,                   METH_VARARGS },
  {"memory_dump",               _memory_dump,           METH_VARARGS },
  {"get_smiles",                get_smiles,             METH_VARARGS },
  {"insert_smiles",             insert_smiles,          METH_VARARGS },
  {"insert_model",              insert_model,           METH_VARARGS },
  {"list_prepend_smiles_list",  list_prepend_smiles_list, METH_VARARGS },
  {"list_get_pattern_list",     list_get_pattern_list,       METH_VARARGS },
  {"list_get_smiles_list",      list_get_smiles_list,       METH_VARARGS },
  {"list_free",                 list_new,               METH_VARARGS },
  {"list_new",                  list_new,              METH_VARARGS },
  {"match_1v1_b",                 match_1v1_b,              METH_VARARGS },
  {"match_1v1_map",             match_1v1_map,           METH_VARARGS},
  {"match_1vN_n",                 match_1vN_n,              METH_VARARGS },
  {NULL,		                    NULL}     /* sentinel */        
};

void init_champ(void);
void init_champ(void)
{
  Py_InitModule("_champ", champ_methods);
}

