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

/*static PyObject *RetNone()
{
  return(Py_BuildValue("(iO)",0,Py_None));
}
*/

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

typedef void(*CObjectDestruct)(void*) ;

static PyObject *_new(PyObject *self,      PyObject *args)
{
  return(PyCObject_FromVoidPtr((void*)ChampNew(),(CObjectDestruct)ChampFree));
}

static PyObject *_memory_dump(PyObject *self,      PyObject *args)
{
  int ok=true;
  PyObject *O;
  CChamp *I;
  ok = PyArg_ParseTuple(args,"O",&O);
  ok = PyCObject_Check(O);
  if(ok) {
    I = PyCObject_AsVoidPtr(O);
    ChampMemoryDump(I);    
  }
  return(RetStatus(ok));
}

static PyObject *insert_pattern_string(PyObject *self,      PyObject *args)
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


static PyObject *pattern_free(PyObject *self,      PyObject *args)
{
  int ok=true;
  int int1;
  PyObject *O;
  CChamp *I;
  ok = PyArg_ParseTuple(args,"Oi",&O,&int1);
  ok = PyCObject_Check(O);
  if(ok) {
    I = PyCObject_AsVoidPtr(O);
    ChampPatFree(I,int1);
  }
  return(RetStatus(ok));
}


static PyObject *pattern_orient_bonds(PyObject *self,      PyObject *args)
{
  int ok=true;
  int int1;
  PyObject *O;
  CChamp *I;
  ok = PyArg_ParseTuple(args,"Oi",&O,&int1);
  ok = PyCObject_Check(O);
  if(ok) {
    I = PyCObject_AsVoidPtr(O);
    ChampOrientBonds(I,int1);
  }
  return(RetStatus(ok));
}

static PyObject *pattern_detect_chirality(PyObject *self,      PyObject *args)
{
  int ok=true;
  int int1;
  PyObject *O;
  CChamp *I;
  ok = PyArg_ParseTuple(args,"Oi",&O,&int1);
  ok = PyCObject_Check(O);
  if(ok) {
    I = PyCObject_AsVoidPtr(O);
    ChampDetectChirality(I,int1);
  }
  return(RetStatus(ok));
}


static PyObject *pattern_generalize(PyObject *self,      PyObject *args)
{
  int ok=true;
  int int1;
  PyObject *O;
  CChamp *I;
  ok = PyArg_ParseTuple(args,"Oi",&O,&int1);
  ok = PyCObject_Check(O);
  if(ok) {
    I = PyCObject_AsVoidPtr(O);
    ChampGeneralize(I,int1);
  }
  return(RetStatus(ok));
}

static PyObject *pattern_get_cycle(PyObject *self,      PyObject *args)
{
  int ok=true;
  int int1;
  PyObject *result = NULL;
  PyObject *O,*l1,*l2;
  CChamp *I;
  ListPat *pat;
  ListAtom *at;
  ListBond *bd;
  int a,b;
  int n_atom,n_bond;
  ok = PyArg_ParseTuple(args,"Oi",&O,&int1);
  ok = PyCObject_Check(O);
  if(ok) {
    I = PyCObject_AsVoidPtr(O);
    pat = I->Pat+int1;

    n_atom = ListLen(I->Atom,pat->atom);
    at = I->Atom + pat->atom;
    l1 = PyList_New(n_atom);
    for(a=0;a<n_atom;a++) {
      PyList_SetItem(l1,a,PyInt_FromLong(at->cycle));
      at = I->Atom + at->link;
    }

    n_bond = ListLen(I->Bond,pat->bond);
    l2 = PyList_New(n_bond);
    bd = I->Bond + pat->bond;
    for(b=0;b<n_bond;b++) {
      PyList_SetItem(l2,b,PyInt_FromLong(bd->cycle));
      bd = I->Bond + bd->link;
    }

    result = PyList_New(2);
    PyList_SetItem(result,0,l1);
    PyList_SetItem(result,1,l2);
  }
  return(RetObj(ok,result));
}


static PyObject *pattern_get_class(PyObject *self,      PyObject *args)
{
  int ok=true;
  int int1;
  PyObject *result = NULL;
  PyObject *O,*l1,*l2;
  CChamp *I;
  ListPat *pat;
  ListAtom *at;
  ListBond *bd;
  int a,b;
  int n_atom,n_bond;
  ok = PyArg_ParseTuple(args,"Oi",&O,&int1);
  ok = PyCObject_Check(O);
  if(ok) {
    I = PyCObject_AsVoidPtr(O);
    pat = I->Pat+int1;
    
    n_atom = ListLen(I->Atom,pat->atom);
    at = I->Atom + pat->atom;
    l1 = PyList_New(n_atom);
    for(a=0;a<n_atom;a++) {
      PyList_SetItem(l1,a,PyInt_FromLong(at->class));
      at = I->Atom + at->link;
    }

    n_bond = ListLen(I->Bond,pat->bond);
    l2 = PyList_New(n_bond);
    bd = I->Bond + pat->bond;
    for(b=0;b<n_bond;b++) {
      PyList_SetItem(l2,b,PyInt_FromLong(bd->class));
      bd = I->Bond + bd->link;
    }

    result = PyList_New(2);
    PyList_SetItem(result,0,l1);
    PyList_SetItem(result,1,l2);
  }
  return(RetObj(ok,result));
}

static PyObject *pattern_get_codes(PyObject *self,      PyObject *args)
{
  int ok=true;
  int int1,ai,bi;
  PyObject *result = NULL;
  PyObject *O,*l1,*l2;
  CChamp *I;
  ListPat *pat;
  ListAtom *at;
  ListBond *bd;
  int a,b;
  int n_atom,n_bond;
  char code[255],atom[10];

  ok = PyArg_ParseTuple(args,"Oi",&O,&int1);
  ok = PyCObject_Check(O);
  if(ok) {
    I = PyCObject_AsVoidPtr(O);
    pat = I->Pat+int1;
    n_atom = ListLen(I->Atom,pat->atom);
    ai = pat->atom;
    l1 = PyList_New(n_atom);
    for(a=0;a<n_atom;a++) {
      at = I->Atom + ai;

      if(at->class&cH_Aliphatic)
        code[0]='A';
      else if(at->class&cH_Aromatic)
        code[0]='R';
      else
        code[0]='P';
      if(at->cycle&(cH_Ring3|cH_Ring4|cH_Ring5|cH_Ring6|cH_Ring7))
        code[1]='C';
      else 
        code[1]='X';
      code[2]=0;
      ChampAtomToString(I,ai,atom);
      if(atom[0]>96)
        atom[0]-=32;
      strcat(code,atom);
      PyList_SetItem(l1,a,PyString_FromString(code));
      ai = at->link;
    }

    n_bond = ListLen(I->Bond,pat->bond);
    l2 = PyList_New(n_bond);
    bi = pat->bond;
    for(b=0;b<n_bond;b++) {
      bd = I->Bond + bi;
      if(bd->class&cH_Aliphatic)
        code[0]='A';
      else if(bd->class&cH_Aromatic)
        code[0]='R';
      else
        code[0]='P';
      if(bd->cycle&(cH_Ring3|cH_Ring4|cH_Ring5|cH_Ring6|cH_Ring7))
        code[1]='C';
      else 
        code[1]='X';

      switch(bd->order) {
      case cH_Single: code[2]='-'; break;
      case cH_Double: code[2]='='; break;
      case cH_Triple: code[2]='#'; break;
      }
      code[3]=0;

      PyList_SetItem(l2,b,PyString_FromString(code));
      bi = bd->link;
    }

    result = PyList_New(2);
    PyList_SetItem(result,0,l1);
    PyList_SetItem(result,1,l2);
  }
  return(RetObj(ok,result));
}


static PyObject *pattern_get_tags(PyObject *self,      PyObject *args)
{
  int ok=true;
  int int1,ai,bi;
  PyObject *result = NULL;
  PyObject *O,*l1,*l2,*l3;
  CChamp *I;
  ListPat *pat;
  ListAtom *at;
  ListBond *bd;
  int a,b,c;
  int n_atom,n_bond;
  int bit_mask;
  int bit_no;
  int n_bits;

  ok = PyArg_ParseTuple(args,"Oi",&O,&int1);
  ok = PyCObject_Check(O);
  if(ok) {
    I = PyCObject_AsVoidPtr(O);
    pat = I->Pat+int1;
    n_atom = ListLen(I->Atom,pat->atom);
    ai = pat->atom;
    l1 = PyList_New(n_atom);
    for(a=0;a<n_atom;a++) {
      at = I->Atom + ai;

      n_bits = 0;
      bit_mask = at->tag;
      while(bit_mask) { /* count bits */
        if(bit_mask&0x1) 
          n_bits++;
        bit_mask=((bit_mask>>1)&0x7FFFFFFF);
      }
      l3 = PyList_New(n_bits);
      bit_no = 0;
      bit_mask = at->tag;
      c=0;
      for(bit_no=0;bit_no<32;bit_no++) {
        if(bit_mask&0x1) {
          PyList_SetItem(l3,c,PyInt_FromLong(bit_no));
          c++;
        }
        bit_mask=((bit_mask>>1)&0x7FFFFFFF);
      }
      PyList_SetItem(l1,a,l3);

      ai = at->link;
    }

    n_bond = ListLen(I->Bond,pat->bond);
    l2 = PyList_New(n_bond);
    bi = pat->bond;
    for(b=0;b<n_bond;b++) {
      bd = I->Bond + bi;
      n_bits = 0;
      bit_mask = bd->tag;
      while(bit_mask) { /* count bits */
        if(bit_mask&0x1) 
          n_bits++;
        bit_mask=((bit_mask>>1)&0x7FFFFFFF);
      }
      l3 = PyList_New(n_bits);
      bit_no = 0;
      bit_mask = bd->tag;
      c=0;
      for(bit_no=0;bit_no<32;bit_no++) {
        if(bit_mask&0x1) {
          PyList_SetItem(l3,c,PyInt_FromLong(bit_no));
          c++;
        }
        bit_mask=((bit_mask>>1)&0x7FFFFFFF);
      }
      PyList_SetItem(l2,b,l3);
      bi = bd->link;
    }

    result = PyList_New(2);
    PyList_SetItem(result,0,l1);
    PyList_SetItem(result,1,l2);
  }
  return(RetObj(ok,result));
}

static PyObject *pattern_get_ext_indices_with_tags(PyObject *self,      PyObject *args)
{
  int ok=true;
  int int1,ai,bi;
  PyObject *result = NULL;
  PyObject *O,*l1,*l2,*l3,*l4;
  CChamp *I;
  ListPat *pat;
  ListAtom *at;
  ListBond *bd;
  int a,b,c;
  int n_atom,n_bond;
  int bit_mask;
  int bit_no;
  int n_bits;

  ok = PyArg_ParseTuple(args,"Oi",&O,&int1);
  ok = PyCObject_Check(O);
  if(ok) {
    I = PyCObject_AsVoidPtr(O);
    pat = I->Pat+int1;
    n_atom = ListLen(I->Atom,pat->atom);
    ai = pat->atom;
    l1 = PyList_New(n_atom);
    for(a=0;a<n_atom;a++) {
      at = I->Atom + ai;

      n_bits = 0;
      bit_mask = at->tag;
      while(bit_mask) { /* count bits */
        if(bit_mask&0x1) 
          n_bits++;
        bit_mask=((bit_mask>>1)&0x7FFFFFFF);
      }
      l4 = PyList_New(2);
      l3 = PyList_New(n_bits);
      bit_no = 0;
      bit_mask = at->tag;
      c=0;
      for(bit_no=0;bit_no<32;bit_no++) {
        if(bit_mask&0x1) {
          PyList_SetItem(l3,c,PyInt_FromLong(bit_no));
          c++;
        }
        bit_mask=((bit_mask>>1)&0x7FFFFFFF);
      }
      PyList_SetItem(l1,a,l4);
      PyList_SetItem(l4,0,PyInt_FromLong(at->ext_index));
      PyList_SetItem(l4,1,l3);
      ai = at->link;
    }

    n_bond = ListLen(I->Bond,pat->bond);
    l2 = PyList_New(n_bond);
    bi = pat->bond;
    for(b=0;b<n_bond;b++) {
      bd = I->Bond + bi;
      n_bits = 0;
      bit_mask = bd->tag;
      while(bit_mask) { /* count bits */
        if(bit_mask&0x1) 
          n_bits++;
        bit_mask=((bit_mask>>1)&0x7FFFFFFF);
      }
      l4 = PyList_New(2);
      l3 = PyList_New(n_bits);
      bit_no = 0;
      bit_mask = bd->tag;
      c=0;
      for(bit_no=0;bit_no<32;bit_no++) {
        if(bit_mask&0x1) {
          PyList_SetItem(l3,c,PyInt_FromLong(bit_no));
          c++;
        }
        bit_mask=((bit_mask>>1)&0x7FFFFFFF);
      }
      PyList_SetItem(l2,b,l4);
      PyList_SetItem(l4,0,PyInt_FromLong(bd->ext_index));
      PyList_SetItem(l4,1,l3);

      bi = bd->link;
    }

    result = PyList_New(2);
    PyList_SetItem(result,0,l1);
    PyList_SetItem(result,1,l2);
  }
  return(RetObj(ok,result));
}

static PyObject *pattern_get_tag_masks(PyObject *self,      PyObject *args)
{
  int ok=true;
  int int1,ai,bi;
  PyObject *result = NULL;
  PyObject *O,*l1,*l2;
  CChamp *I;
  ListPat *pat;
  ListAtom *at;
  ListBond *bd;
  int a,b;
  int n_atom,n_bond;

  ok = PyArg_ParseTuple(args,"Oi",&O,&int1);
  ok = PyCObject_Check(O);
  if(ok) {
    I = PyCObject_AsVoidPtr(O);
    pat = I->Pat+int1;
    n_atom = ListLen(I->Atom,pat->atom);
    ai = pat->atom;
    l1 = PyList_New(n_atom);
    for(a=0;a<n_atom;a++) {
      at = I->Atom + ai;

      PyList_SetItem(l1,a,PyInt_FromLong(at->tag));

      ai = at->link;
    }

    n_bond = ListLen(I->Bond,pat->bond);
    l2 = PyList_New(n_bond);
    bi = pat->bond;
    for(b=0;b<n_bond;b++) {
      bd = I->Bond + bi;
      PyList_SetItem(l2,b,PyInt_FromLong(bd->tag));
      bi = bd->link;
    }

    result = PyList_New(2);
    PyList_SetItem(result,0,l1);
    PyList_SetItem(result,1,l2);
  }
  return(RetObj(ok,result));
}


static PyObject *pattern_clear_tags(PyObject *self,      PyObject *args)
{
  int ok=true;
  int int1,ai,bi;
  PyObject *O;
  CChamp *I;
  ListPat *pat;
  ListAtom *at;
  ListBond *bd;

  ok = PyArg_ParseTuple(args,"Oi",&O,&int1);
  ok = PyCObject_Check(O);
  if(ok) {
    I = PyCObject_AsVoidPtr(O);
    pat = I->Pat+int1;
    ai = pat->atom;
    while(ai) {
      at = I->Atom + ai;
      at->tag = 0;
      ai = at->link;
    }

    bi = pat->bond;
    while(bi) {
      bd = I->Bond + bi;
      bd->tag = 0;
      bi = bd->link;
    }
  }
  return(RetStatus(ok));
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

static PyObject *list_get_pattern_indices(PyObject *self,      PyObject *args)
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

static PyObject *pattern_get_string(PyObject *self,      PyObject *args)
{
  int ok=true;
  PyObject *result = NULL;
  PyObject *O;
  int pat_index;
  int mode;
  CChamp *I;
  char *smi;

  ok = PyArg_ParseTuple(args,"Oii",&O,&pat_index,&mode);
  ok = PyCObject_Check(O);
  if(ok) {
    I = PyCObject_AsVoidPtr(O);
    smi = ChampPatToSmiVLA(I,pat_index,NULL,mode);
    result = PyString_FromString(smi);
    vla_free(smi);
  }
  return(RetObj(ok,result));
}

static PyObject *pattern_dump(PyObject *self,      PyObject *args)
{
  int ok=true;
  PyObject *result = NULL;
  PyObject *O;
  int pat_index;
  CChamp *I;

  ok = PyArg_ParseTuple(args,"Oi",&O,&pat_index);
  ok = PyCObject_Check(O);
  if(ok) {
    I = PyCObject_AsVoidPtr(O);
    ChampPatDump(I,pat_index);
  }
  return(RetObj(ok,result));
}

static PyObject *pattern_get_atom_symbols(PyObject *self,      PyObject *args)
{
  int ok=true;
  int int1,ai;
  PyObject *result = NULL;
  PyObject *O;
  CChamp *I;
  ListPat *pat;
  ListAtom *at;
  int a;
  int n_atom;
  char atom[255];

  ok = PyArg_ParseTuple(args,"Oi",&O,&int1);
  ok = PyCObject_Check(O);
  if(ok) {
    I = PyCObject_AsVoidPtr(O);
    pat = I->Pat+int1;
    n_atom = ListLen(I->Atom,pat->atom);
    ai = pat->atom;
    result = PyList_New(n_atom);
    for(a=0;a<n_atom;a++) {
      at = I->Atom + ai;
      ChampAtomToString(I,ai,atom);
      PyList_SetItem(result,a,PyString_FromString(atom));
      ai = at->link;
    }
  }
  return(RetObj(ok,result));
}

static PyObject *pattern_get_atom_names(PyObject *self,      PyObject *args)
{
  int ok=true;
  int int1,ai;
  PyObject *result = NULL;
  PyObject *O;
  CChamp *I;
  ListPat *pat;
  ListAtom *at;
  int a;
  int n_atom;

  ok = PyArg_ParseTuple(args,"Oi",&O,&int1);
  ok = PyCObject_Check(O);
  if(ok) {
    I = PyCObject_AsVoidPtr(O);
    pat = I->Pat+int1;
    n_atom = ListLen(I->Atom,pat->atom);
    ai = pat->atom;
    result = PyList_New(n_atom);
    for(a=0;a<n_atom;a++) {
      at = I->Atom + ai;
      PyList_SetItem(result,a,PyString_FromString(at->name));
      ai = at->link;
    }
  }
  return(RetObj(ok,result));
}


static PyObject *list_get_pattern_strings(PyObject *self,      PyObject *args)
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
      smi = ChampPatToSmiVLA(I,I->Int[i].value,smi,0);
      str = PyString_FromString(smi);
      PyList_SetItem(result,c,str);
      i = I->Int[i].link;
      c++;
    }
    vla_free(smi);
  }
  return(RetObj(ok,result));
}

 
static PyObject *list_prepend_pattern_strings(PyObject *self,      PyObject *args)
{
  int ok=true;
  int list_handle,pat_index;
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
      I->Int[list_handle].link = ListElemPushInt(&I->Int,I->Int[list_handle].link,pat_index);
    }
  }
  return(RetStatus(ok));
}
 

static PyObject *list_prepend_pattern_index(PyObject *self,      PyObject *args)
{
  int ok=true;
  int list_handle,pat_index;
  PyObject *O;
  CChamp *I;
  ok = PyArg_ParseTuple(args,"Oii",&O,&list_handle,&pat_index);
  ok = PyCObject_Check(O);
  if(ok) {
    I = PyCObject_AsVoidPtr(O);
    if(pat_index)
      I->Int[list_handle].link = ListElemPushInt(&I->Int,I->Int[list_handle].link,pat_index);
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

static PyObject *match_1v1_n(PyObject *self,      PyObject *args)
{
  int ok=true;
  int result = 0;
  int pat1,pat2;
  int limit;
  int tag;
  PyObject *O;
  CChamp *I;
  ok = PyArg_ParseTuple(args,"Oiiii",&O,&pat1,&pat2,&limit,&tag);
  ok = PyCObject_Check(O);
  if(ok) {
    I = PyCObject_AsVoidPtr(O);
    result = ChampMatch_1V1_N(I,pat1,pat2,limit,tag);
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
  int tag;
  ListMatch *mat;

  PyObject *O;
  CChamp *I;
  ok = PyArg_ParseTuple(args,"Oiiii",&O,&pat1,&pat2,&limit,&tag);
  ok = PyCObject_Check(O);
  if(ok) {
    I = PyCObject_AsVoidPtr(O);
    mat_start = ChampMatch_1V1_Map(I,pat1,pat2,limit,tag);

    ChampPatReindex(I,pat1);
    ChampPatReindex(I,pat2);

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
        /*
                printf("%2d %2d %4d %4d\n",I->Int2[j].value[0],I->Int2[j].value[1],
                       I->Atom[I->Int2[j].value[0]].index,
                       I->Atom[I->Int2[j].value[1]].index);*/
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

static PyObject *exact_1vN_n(PyObject *self,      PyObject *args)
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
    result = ChampExact_1VN_N(I,pattern,list_index);
  }
  return(RetInt(ok,result));
}


static PyObject *match_Nv1_n(PyObject *self,      PyObject *args)
{
  int ok=true;
  int result = 0;
  int list_handle,list_index;
  int pattern;
  int limit,tag;
  PyObject *O;
  CChamp *I;
  ok = PyArg_ParseTuple(args,"Oiiii",&O,&list_handle,&pattern,&limit,&tag);
  ok = PyCObject_Check(O);
  if(ok) {
    I = PyCObject_AsVoidPtr(O);
    list_index = I->Int[list_handle].link;
    result = ChampMatch_NV1_N(I,list_index,pattern,limit,tag);
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
  {"insert_pattern_string",     insert_pattern_string,          METH_VARARGS },
  {"insert_model",              insert_model,           METH_VARARGS },
  {"pattern_free",              pattern_free,          METH_VARARGS },
  {"pattern_clear_tags",        pattern_clear_tags,      METH_VARARGS },
  {"pattern_get_cycle",   pattern_get_cycle,       METH_VARARGS },
  {"pattern_get_class",   pattern_get_class,       METH_VARARGS },
  {"pattern_get_codes",   pattern_get_codes,       METH_VARARGS },
  {"pattern_get_string",        pattern_get_string,       METH_VARARGS },
  {"pattern_get_tags",    pattern_get_tags,       METH_VARARGS },
  {"pattern_get_tag_masks",    pattern_get_tag_masks,       METH_VARARGS },
  {"pattern_get_atom_symbols",  pattern_get_atom_symbols,     METH_VARARGS },
  {"pattern_get_atom_names",  pattern_get_atom_names,     METH_VARARGS },
  {"pattern_get_ext_indices_with_tags",  pattern_get_ext_indices_with_tags, METH_VARARGS },
  {"pattern_dump",  pattern_dump,     METH_VARARGS },
  {"pattern_orient_bonds",   pattern_orient_bonds,       METH_VARARGS },
  {"pattern_detect_chirality",   pattern_detect_chirality,       METH_VARARGS },
  {"pattern_generalize",   pattern_generalize,     METH_VARARGS },
  {"list_prepend_pattern_strings",  list_prepend_pattern_strings, METH_VARARGS },
  {"list_prepend_pattern_index",  list_prepend_pattern_index, METH_VARARGS },
  {"list_get_pattern_indices",  list_get_pattern_indices,       METH_VARARGS },
  {"list_get_pattern_strings",      list_get_pattern_strings,       METH_VARARGS },
  {"list_free",                 list_free,               METH_VARARGS },
  {"list_new",                  list_new,              METH_VARARGS },
  {"match_1v1_b",                 match_1v1_b,              METH_VARARGS },
  {"match_1v1_map",             match_1v1_map,           METH_VARARGS},
  {"match_1v1_n",               match_1v1_n,             METH_VARARGS},
  {"match_1vN_n",               match_1vN_n,              METH_VARARGS },
  {"match_Nv1_n",               match_Nv1_n,              METH_VARARGS },
  /*  {"map_1v1_to_indexed_strings",  map_1v1_to_indexed_strings, METH_VARARGS },*/
  {"exact_1vN_n",               exact_1vN_n,              METH_VARARGS },

  {NULL,		                    NULL}     /* sentinel */        
};

void init_champ(void);
void init_champ(void)
{
  Py_InitModule("_champ", champ_methods);
}

