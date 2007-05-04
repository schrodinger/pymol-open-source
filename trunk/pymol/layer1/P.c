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

/* meaning of defines 

_PYMOL_MONOLITHIC: means that we're building PyMOL and its Python C
dependencies as one C library.  That means we need to explicitly call
the initialization functions for these libraries on startup.

*/

#ifndef _PYMOL_NOPY

#include"os_predef.h"
#include"Base.h"

#define MAX_SAVED_THREAD ((PYMOL_MAX_THREADS)+3)


/* BEGIN PROPRIETARY CODE SEGMENT (see disclaimer in "os_proprietary.h") */
#ifdef WIN32
#include<windows.h>
#include<process.h>
#include<winappc.h>
#endif
/* END PROPRIETARY CODE SEGMENT */

#include"os_python.h"

#include"os_std.h"
#include"os_time.h"
#include"os_unix.h"

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
#include"Executive.h"
#include"PyMOLOptions.h"
#include"PyMOL.h"



/* all of the following Python objects must be invariant & global for the application */

/* local to this C code module */

static PyObject *P_pymol = NULL;
static PyObject *P_pymol_dict = NULL; /* must be refomed into globals and instance properties */
static PyObject *P_cmd = NULL; 

/* used elsewhere */

PyObject *P_menu = NULL; /* menu definitions are currently global */
PyObject *P_xray = NULL; /* okay as global */
PyObject *P_chempy = NULL; /* okay as global */
PyObject *P_models = NULL; /* okay as global */
PyObject *P_setting = NULL; /* must be reformed somehow */
PyObject *P_embed = NULL;

/* local to this module */

static PyObject *P_povray = NULL;
static PyObject *P_traceback = NULL;
static PyObject *P_parser = NULL; 

static PyObject *P_lock = NULL; /* API locks */
static PyObject *P_lock_attempt = NULL;
static PyObject *P_unlock = NULL;

static PyObject *P_lock_c = NULL; /* C locks */
static PyObject *P_unlock_c = NULL;

static PyObject *P_lock_status = NULL; /* status locks */
static PyObject *P_lock_status_attempt = NULL; /* status locks */
static PyObject *P_unlock_status = NULL;

static PyObject *P_lock_glut = NULL; /* GLUT locks */
static PyObject *P_unlock_glut = NULL;

/* BEGIN PROPRIETARY CODE SEGMENT (see disclaimer in "os_proprietary.h") */
#ifdef WIN32
static PyObject *P_time = NULL;
static PyObject *P_sleep = NULL;
#endif
/* END PROPRIETARY CODE SEGMENT */

static PyObject *P_main = NULL;
static PyObject *P_vfont = NULL;

#define P_log_file_str "_log_file"

#define xxxPYMOL_NEW_THREADS

unsigned int PyThread_get_thread_ident(void); /* critical functionality */

typedef struct {
  int id;
  PyThreadState *state;
} SavedThreadRec;

void PRunStringModule(PyMOLGlobals *G,char *str);

void PLockStatus(void) /* assumes we have the GIL */
{
  PXDecRef(PyObject_CallFunction(P_lock_status,NULL));
}

int PLockStatusAttempt(void) /* assumes we have the GIL */
{
  int result = true;
  PyObject *got_lock = PyObject_CallFunction(P_lock_status_attempt,NULL);
  if(got_lock) {
    if(!PyInt_AsLong(got_lock)) {
      result = false;
    }
    Py_DECREF(got_lock);
  }
  return result;
}
void PUnlockStatus(void) /* assumes we have the GIL */
{
  PXDecRef(PyObject_CallFunction(P_unlock_status,NULL));
}

static void PLockGLUT(void) /* assumes we have the GIL */
{
  PXDecRef(PyObject_CallFunction(P_lock_glut,NULL));
}

static void PUnlockGLUT(void) /* assumes we have the GIL */
{
  PXDecRef(PyObject_CallFunction(P_unlock_glut,NULL));
}

static SavedThreadRec SavedThread[MAX_SAVED_THREAD];

int P_glut_thread_keep_out = 0; 
unsigned int P_glut_thread_id;

/* enables us to keep glut out if by chance it grabs the API
 * in the middle of a nested API based operation */

void PCatchInit(void);
void my_interrupt(int a);
char *getprogramname(void);

/*
PyObject *GetBondsDict(PyMOLGlobals *G)
{
  PyObject *result = NULL;
  result = PyObject_GetAttrString(P_chempy,"bonds");
  if(!result) ErrMessage(G,"PyMOL","can't find 'chempy.bonds.bonds'");
  return(result);
}
*/

PyObject *PGetFontDict(PyMOLGlobals *G, float size,int face,int style)
{ /* assumes we have a valid interpreter lock */
  PyObject *result = NULL; 

  if(!P_vfont) {
    PRunStringModule(G,"import vfont\n");  
    P_vfont = PyDict_GetItemString(P_pymol_dict,"vfont");
  }
  if(!P_vfont) {
    PRINTFB(G,FB_Python,FB_Errors)
      " PyMOL-Error: can't find module 'vfont'"
      ENDFB(G);
  }
  else {
    result = PyObject_CallMethod(P_vfont,"get_font","fii",size,face,style);
  }
  return(PConvAutoNone(result));
}

int PComplete(PyMOLGlobals *G,char *str,int buf_size)
{
  int ret = false;
  PyObject *result;
  char *st2;
  PBlockAndUnlockAPI(G);
  if(G->P_inst->complete) {
    result = PyObject_CallFunction(G->P_inst->complete,"s",str);
    if(result) {
      if(PyString_Check(result)) {
        st2 = PyString_AsString(result);
        UtilNCopy(str,st2,buf_size);
        ret=true;
      }
      Py_DECREF(result);
    }
  }
  PLockAPIAndUnblock(G);
  return(ret);
}

int PTruthCallStr0(PyObject *object,char *method)
{
  int result = false;
  PyObject *tmp;
  tmp = PyObject_CallMethod(object,method,"");
  if(tmp) {
    if(PyObject_IsTrue(tmp))
      result = 1;
    Py_DECREF(tmp);
  }
  return(result);
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

int PTruthCallStr1i(PyObject *object,char *method,int argument)
{
  int result = false;
  PyObject *tmp;
  tmp = PyObject_CallMethod(object,method,"i",argument);
  if(tmp) {
    if(PyObject_IsTrue(tmp))
      result = 1;
    Py_DECREF(tmp);
  }
  return(result);
}
int PTruthCallStr4i(PyObject *object,char *method,int a1,int a2,int a3,int a4)
{
  int result = false;
  PyObject *tmp;
  tmp = PyObject_CallMethod(object,method,"iiii",a1,a2,a3,a4);
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

void PSleepWhileBusy(PyMOLGlobals *G,int usec)
{
#ifndef WIN32
  struct timeval tv;
  PRINTFD(G,FB_Threads)
    " PSleep-DEBUG: napping.\n"
  ENDFD;
  tv.tv_sec=0;
  tv.tv_usec=usec; 
  select(0,NULL,NULL,NULL,&tv);
  PRINTFD(G,FB_Threads)
    " PSleep-DEBUG: nap over.\n"
  ENDFD;
#else
/* BEGIN PROPRIETARY CODE SEGMENT (see disclaimer in "os_proprietary.h") */
  PBlock(G);
  PXDecRef(PyObject_CallFunction(P_sleep,"f",usec/1000000.0));
  PUnblock(G);
/* END PROPRIETARY CODE SEGMENT */
#endif
}

void PSleepUnlocked(PyMOLGlobals *G,int usec)
{ /* can only be called by the glut process */
#ifndef WIN32
  struct timeval tv;
  PRINTFD(G,FB_Threads)
    " PSleep-DEBUG: napping.\n"
  ENDFD;
  tv.tv_sec=0;
  tv.tv_usec=usec; 
  select(0,NULL,NULL,NULL,&tv);
  PRINTFD(G,FB_Threads)
    " PSleep-DEBUG: nap over.\n"
  ENDFD;
#else
/* BEGIN PROPRIETARY CODE SEGMENT (see disclaimer in "os_proprietary.h") */
  PBlock(G);
  PXDecRef(PyObject_CallFunction(P_sleep,"f",usec/1000000.0));
  PUnblock(G);
/* END PROPRIETARY CODE SEGMENT */
#endif
}

void PSleep(PyMOLGlobals *G,int usec)
{ /* can only be called by the glut process */
#ifndef WIN32
  struct timeval tv;
  PUnlockAPIAsGlut(G);
  PRINTFD(G,FB_Threads)
    " PSleep-DEBUG: napping.\n"
  ENDFD;
  tv.tv_sec=0;
  tv.tv_usec=usec; 
  select(0,NULL,NULL,NULL,&tv);
  PRINTFD(G,FB_Threads)
    " PSleep-DEBUG: nap over.\n"
  ENDFD;
  PLockAPIAsGlut(G,true);
#else
/* BEGIN PROPRIETARY CODE SEGMENT (see disclaimer in "os_proprietary.h") */
  PBlockAndUnlockAPI(G);
  PXDecRef(PyObject_CallFunction(P_sleep,"f",usec/1000000.0));
  PLockAPIAndUnblock(G);
/* END PROPRIETARY CODE SEGMENT */
#endif

}

static PyObject *PCatchWrite(PyObject *self, 	PyObject *args);

void my_interrupt(int a)
{
  exit(EXIT_FAILURE);
}

void PDumpTraceback(PyObject *err)
{
  PyObject_CallMethod(P_traceback,"print_tb","O",err);
}

void PDumpException()
{
  PyObject_CallMethod(P_traceback,"print_exc","");
}

int PAlterAtomState(PyMOLGlobals *G,float *v,char *expr,int read_only,
                    AtomInfoType *at,char *model,int index,
                    PyObject *space) 
     /* assumes Blocked python interpreter*/
{
  PyObject *dict; 
  int result=true;
  float f[3];
  PyObject *x_id1,*x_id2=NULL,*y_id1,*y_id2=NULL,*z_id1,*z_id2=NULL;
  char atype[7];
  PyObject *flags_id1=NULL,*flags_id2=NULL;
  int flags;
  dict = PyDict_New();
  
  if(at) {
    if(at->hetatm)
      strcpy(atype,"HETATM");
    else
      strcpy(atype,"ATOM");

    /* immutables */
    PConvStringToPyDictItem(dict,"model",model);
    PConvIntToPyDictItem(dict,"index",index+1);
    
    PConvStringToPyDictItem(dict,"type",atype);
    PConvStringToPyDictItem(dict,"name",at->name);
    PConvStringToPyDictItem(dict,"resn",at->resn);
    PConvStringToPyDictItem(dict,"resi",at->resi);
    PConvIntToPyDictItem(dict,"resv",at->resv); /* subordinate to resi */
    PConvStringToPyDictItem(dict,"chain",at->chain);
    PConvStringToPyDictItem(dict,"alt",at->alt);
    PConvStringToPyDictItem(dict,"segi",at->segi);
    PConvStringToPyDictItem(dict,"elem",at->elem);
    PConvStringToPyDictItem(dict,"ss",at->ssType);

    {
      char null_st[1] = "";
      char *st = null_st;
      
      if(at->textType) st = OVLexicon_FetchCString(G->Lexicon,at->textType);
      PConvStringToPyDictItem(dict,"text_type",st);
      
      st = null_st;
      if(at->label) st = OVLexicon_FetchCString(G->Lexicon,at->label);
      PConvStringToPyDictItem(dict,"label",st);
    }
    PConvIntToPyDictItem(dict,"numeric_type",at->customType);
    PConvFloatToPyDictItem(dict,"q",at->q);
    PConvFloatToPyDictItem(dict,"b",at->b);
    PConvFloatToPyDictItem(dict,"vdw",at->vdw);
    PConvFloatToPyDictItem(dict,"elec_radius",at->elec_radius);
    PConvFloatToPyDictItem(dict,"partial_charge",at->partialCharge);
    PConvIntToPyDictItem(dict,"formal_charge",at->formalCharge);
    PConvIntToPyDictItem(dict,"cartoon",at->cartoon);
    PConvIntToPyDictItem(dict,"color",at->color);
    PConvIntToPyDictItem(dict,"ID",at->id);
    PConvIntToPyDictItem(dict,"rank",at->rank);

    /* mutables */
    flags_id1 = PConvIntToPyDictItem(dict,"flags",at->flags);
  }
  x_id1 = PConvFloatToPyDictItem(dict,"x",v[0]);
  y_id1 = PConvFloatToPyDictItem(dict,"y",v[1]);
  z_id1 = PConvFloatToPyDictItem(dict,"z",v[2]);
  PXDecRef(PyRun_String(expr,Py_single_input,space,dict));
  if(PyErr_Occurred()) {
    PyErr_Print();
    result=false;
  } else if(!read_only) {
    if(result) {
      if(!(x_id2 = PyDict_GetItemString(dict,"x")))
        result=false;
      else if(!(y_id2 = PyDict_GetItemString(dict,"y")))
        result=false;
      else if(!(z_id2 = PyDict_GetItemString(dict,"z")))
        result=false;
      else if(at) {
        if(!(flags_id2 = PyDict_GetItemString(dict,"flags")))
          result=false;
      }
      if(PyErr_Occurred()) {
        PyErr_Print();
        result=false;
        ErrMessage(G,"AlterState","Aborting on error. Assignment may be incomplete.");
      }
    }
    if(result) {
      f[0]=(float)PyFloat_AsDouble(x_id2);
      f[1]=(float)PyFloat_AsDouble(y_id2);
      f[2]=(float)PyFloat_AsDouble(z_id2);
      if(at) 
        if(flags_id1!=flags_id2) {
          if(!PConvPyObjectToInt(flags_id2,&flags))
            result=false;
          else
            at->flags = flags;
        }
      if(PyErr_Occurred()) {
        PyErr_Print();
        result=false;
        ErrMessage(G,"AlterState","Aborting on error. Assignment may be incomplete.");
      } else {
        v[0]=f[0];
        v[1]=f[1];
        v[2]=f[2];
      }
    }
    
  }
  Py_DECREF(dict);
  return result;
}

int PAlterAtom(PyMOLGlobals *G,
               AtomInfoType *at,char *expr,int read_only,
               char *model,int index,PyObject *space)
{
  /* assumes Blocked python interpreter*/
  WordType buf;
  AtomName name;
  PyObject *name_id1,*name_id2=NULL;
  AtomName elem;
  PyObject *elem_id1,*elem_id2=NULL;
  ResName resn;
  PyObject *resn_id1,*resn_id2=NULL;
  ResIdent resi;
  PyObject *resi_id1,*resi_id2=NULL;
  int resv;
  PyObject *resv_id1,*resv_id2=NULL;
  Chain chain;
  PyObject *chain_id1,*chain_id2=NULL;
  Chain alt;
  PyObject *alt_id1,*alt_id2=NULL;
  SegIdent segi;
  PyObject *flags_id1,*flags_id2=NULL;
  int flags;
  PyObject *segi_id1,*segi_id2=NULL;
  PyObject *text_type_id1,*text_type_id2=NULL;
  SSType ssType;
  PyObject *ss_id1,*ss_id2=NULL;
  char atype[7];
  PyObject *type_id1,*type_id2=NULL;
  float b,q,partialCharge,vdw,elec_radius;
  PyObject *b_id1,*b_id2=NULL;
  PyObject *q_id1,*q_id2=NULL;
  PyObject *partial_charge_id1,*partial_charge_id2=NULL;
  PyObject *vdw_id1,*vdw_id2=NULL;
  PyObject *elec_radius_id1,*elec_radius_id2=NULL;
  int formalCharge,numericType;
  PyObject *formal_charge_id1,*formal_charge_id2=NULL;
  PyObject *numeric_type_id1,*numeric_type_id2=NULL;
  int cartoon;
  PyObject *cartoon_id1,*cartoon_id2=NULL;
  int color;
  PyObject *color_id1,*color_id2=NULL;
  PyObject *label_id1,*label_id2=NULL;
  int id;
  PyObject *ID_id1,*ID_id2=NULL;
  int rank;
  PyObject *rank_id1,*rank_id2=NULL;
  PyObject *state_id1,*state_id2=NULL;
  int state;
  PyObject *dict;
  int result=true;
  
  if(at->hetatm)
    strcpy(atype,"HETATM");
  else
    strcpy(atype,"ATOM");

  /* PBlockAndUnlockAPI() is not safe, thus these
   * expressions must not call the PyMOL API...
   * what if "at" is destroyed by another thread? */

  dict = PyDict_New();

  /* immutables */
  PConvStringToPyDictItem(dict,"model",model);
  PConvIntToPyDictItem(dict,"index",index+1);

  /* mutables */
  type_id1 = PConvStringToPyDictItem(dict,"type",atype);
  name_id1 = PConvStringToPyDictItem(dict,"name",at->name);
  resn_id1 = PConvStringToPyDictItem(dict,"resn",at->resn);
  flags_id1 = PConvIntToPyDictItem(dict,"flags",at->flags);
  resi_id1 = PConvStringToPyDictItem(dict,"resi",at->resi);
  resv_id1 = PConvIntToPyDictItem(dict,"resv",at->resv); /* subordinate to resi */
  chain_id1 = PConvStringToPyDictItem(dict,"chain",at->chain);
  alt_id1 = PConvStringToPyDictItem(dict,"alt",at->alt);
  segi_id1 = PConvStringToPyDictItem(dict,"segi",at->segi);
  elem_id1 = PConvStringToPyDictItem(dict,"elem",at->elem);
  ss_id1 = PConvStringToPyDictItem(dict,"ss",at->ssType);
  numeric_type_id1 = PConvIntToPyDictItem(dict,"numeric_type",at->customType);
  q_id1 = PConvFloatToPyDictItem(dict,"q",at->q);
  b_id1 = PConvFloatToPyDictItem(dict,"b",at->b);
  vdw_id1 = PConvFloatToPyDictItem(dict,"vdw",at->vdw);
  elec_radius_id1 = PConvFloatToPyDictItem(dict,"elec_radius",at->elec_radius);
  partial_charge_id1 = PConvFloatToPyDictItem(dict,"partial_charge",at->partialCharge);
  formal_charge_id1 = PConvIntToPyDictItem(dict,"formal_charge",at->formalCharge);
  cartoon_id1 = PConvIntToPyDictItem(dict,"cartoon",at->cartoon);

  {
    char null_st[1] = "";
    char *st = null_st;

    if(at->textType) st = OVLexicon_FetchCString(G->Lexicon,at->textType);
    text_type_id1 = PConvStringToPyDictItem(dict,"text_type",st);

    st = null_st;
    if(at->label) st = OVLexicon_FetchCString(G->Lexicon,at->label);
    label_id1 = PConvStringToPyDictItem(dict,"label",st);
  }

  color_id1 = PConvIntToPyDictItem(dict,"color",at->color);
  ID_id1 = PConvIntToPyDictItem(dict,"ID",at->id);
  state_id1 = PConvIntToPyDictItem(dict,"state",at->discrete_state);
  rank_id1 = PConvIntToPyDictItem(dict,"rank",at->rank);

  PXDecRef(PyRun_String(expr,Py_single_input,space,dict));
  if(PyErr_Occurred()) {
    ErrMessage(G,"Alter","Aborting on error. Assignment may be incomplete.");
    PyErr_Print();
    result=false;
  } else if(read_only) {
    result=true;
  } if(PyErr_Occurred()) {
    PyErr_Print();
    result=false;
  } else if(!read_only) {

    if(result) {
      /* get new object IDs */
      
      if(!(type_id2 = PyDict_GetItemString(dict,"type")))
        result=false;
      else if(!(name_id2 = PyDict_GetItemString(dict,"name")))
        result=false;
      else if(!(elem_id2 = PyDict_GetItemString(dict,"elem")))
        result=false;
      else if(!(resn_id2 = PyDict_GetItemString(dict,"resn")))
        result=false;
      else if(!(flags_id2 = PyDict_GetItemString(dict,"flags")))
        result=false;
      else if(!(resi_id2 = PyDict_GetItemString(dict,"resi")))
        result=false;
      else if(!(resv_id2 = PyDict_GetItemString(dict,"resv")))
        result=false;
      else if(!(segi_id2 = PyDict_GetItemString(dict,"segi")))
        result=false;
      else if(!(alt_id2 = PyDict_GetItemString(dict,"alt")))
        result=false;
      else if(!(chain_id2 = PyDict_GetItemString(dict,"chain")))
        result=false;
      else if(!(text_type_id2 = PyDict_GetItemString(dict,"text_type")))
        result=false;
      else if(!(ss_id2 = PyDict_GetItemString(dict,"ss")))
        result=false;
      else if(!(b_id2 = PyDict_GetItemString(dict,"b")))
        result=false;
      else if(!(q_id2 = PyDict_GetItemString(dict,"q")))
        result=false;
      else if(!(vdw_id2=PyDict_GetItemString(dict,"vdw")))
        result=false;
      else if(!(elec_radius_id2=PyDict_GetItemString(dict,"elec_radius")))
        result=false;
      else if(!(partial_charge_id2 = PyDict_GetItemString(dict,"partial_charge")))
        result=false;
      else if(!(formal_charge_id2 = PyDict_GetItemString(dict,"formal_charge")))
        result=false;
      else if(!(cartoon_id2 = PyDict_GetItemString(dict,"cartoon")))
        result=false;
      else if(!(color_id2=PyDict_GetItemString(dict,"color")))
        result=false;
      else if(!(label_id2=PyDict_GetItemString(dict,"label")))
        result=false;
      if(!(numeric_type_id2 = PyDict_GetItemString(dict,"numeric_type")))
        result=false;
      if(!(ID_id2 = PyDict_GetItemString(dict,"ID")))
        result=false;
      if(!(state_id2 = PyDict_GetItemString(dict,"state")))
        result=false;
      if(!(rank_id2 = PyDict_GetItemString(dict,"rank")))
        result=false;

      if(PyErr_Occurred()) {
        PyErr_Print();
        result=false;
      }
    }
    if(result) {
      if(type_id1!=type_id2) {
        if(!PConvPyObjectToStrMaxLen(type_id2,atype,6))
          result=false;
        else
          at->hetatm=((atype[0]=='h')||(atype[0]=='H'));
      }
      if(name_id1!=name_id2) {
        if(!PConvPyObjectToStrMaxLen(name_id2,name,sizeof(AtomName)-1))
          result=false;
        else
          strcpy(at->name,name);
      }
      if(elem_id1!=elem_id2) {
        if(!PConvPyObjectToStrMaxLen(elem_id2,elem,sizeof(AtomName)-1)) 
          result=false;
        else {
          strcpy(at->elem,elem);
          AtomInfoAssignParameters(G,at);
        }
      }
      if(resn_id1!=resn_id2) {
        if(!PConvPyObjectToStrMaxLen(resn_id2,resn,sizeof(ResName)-1))
          result=false;
        else
          strcpy(at->resn,resn);
      }
      if(resi_id1!=resi_id2) {
        if(!PConvPyObjectToStrMaxLen(resi_id2,resi,sizeof(ResIdent)-1))
          result=false;
        else {
          if(strcmp(at->resi,resi)!=0)
            at->resv=AtomResvFromResi(resi);
          strcpy(at->resi,resi);
        }
      } else if (resv_id1!=resv_id2) {
        if(!PConvPyObjectToInt(resv_id2,&resv))
          result=false;
        else {
          sprintf(buf,"%d",resv);
          buf[sizeof(ResIdent)-1]=0;
          strcpy(at->resi,buf);
        }
        
      }
      if(segi_id1!=segi_id2) {
        if(!PConvPyObjectToStrMaxLen(segi_id2,segi,sizeof(SegIdent)-1))
          result=false;
        else
          strcpy(at->segi,segi);

      }
      if(chain_id1!=chain_id2) {
        if(!PConvPyObjectToStrMaxLen(chain_id2,chain,sizeof(Chain)-1))
          result=false;
        else
          strcpy(at->chain,chain);
      }
      if(alt_id1!=alt_id2) {
        if(!PConvPyObjectToStrMaxLen(alt_id2,alt,sizeof(Chain)-1))
          result=false;
        else
          strcpy(at->alt,alt);
      }
      if(text_type_id1!=text_type_id2) {
        
        OrthoLineType temp;
        if(at->textType) {
          OVLexicon_DecRef(G->Lexicon,at->textType);
        }
        at->textType = 0;

        if(PConvPyObjectToStrMaxLen(text_type_id2,temp,sizeof(OrthoLineType)-1)) {
          if(temp[0]) {
            OVreturn_word result = OVLexicon_GetFromCString(G->Lexicon,temp);
            if(OVreturn_IS_OK(result)) {
              at->textType = result.word;
            }
          }
        }
      }
      if(ss_id1!=ss_id2) {
        if(!PConvPyObjectToStrMaxLen(ss_id2,ssType,sizeof(SSType)-1))
          result=false;
        else {
          strcpy(at->ssType,ssType);
          at->ssType[0] = toupper(at->ssType[0]);
        }
      }
      if(b_id1!=b_id2) {
        if(!PConvPyObjectToFloat(b_id2,&b))
          result=false;
        else
          at->b=b;
      }
      if(q_id1!=q_id2) {
        if(!PConvPyObjectToFloat(q_id2,&q))
          result=false;
        else
          at->q=q;
      }
      if(vdw_id1!=vdw_id2) {
        if(!PConvPyObjectToFloat(vdw_id2,&vdw))
          result=false;
        else
          at->vdw=vdw;

      }
      if(elec_radius_id1!=elec_radius_id2) {
        if(!PConvPyObjectToFloat(elec_radius_id2,&elec_radius))
          result=false;
        else
          at->elec_radius=elec_radius;
      }
      if(partial_charge_id1!=partial_charge_id2) {
        if(!PConvPyObjectToFloat(partial_charge_id2,&partialCharge))
          result=false;
        else
          at->partialCharge=partialCharge;

      }
      if(formal_charge_id1!=formal_charge_id2) {
        if(!PConvPyObjectToInt(formal_charge_id2,&formalCharge))
          result=false;
        else {
          at->formalCharge=formalCharge;
          at->chemFlag = false; /* invalidate chemistry info for this atom */
        }

      }
      if(cartoon_id1!=cartoon_id2) {
        if(!PConvPyObjectToInt(cartoon_id2,&cartoon))
          result=false;
        else
          at->cartoon=cartoon;
      }
      if(color_id1!=color_id2) {
        if(!PConvPyObjectToInt(color_id2,&color))
          result=false;
        else
          at->color=color;
      }
      if(state_id1!=state_id2) {
        if(!PConvPyObjectToInt(state_id2,&state))
          result=false;
        else
          at->discrete_state=state;
      }
      if(label_id1!=label_id2) {
        OrthoLineType temp;
        if(at->label) {
          OVLexicon_DecRef(G->Lexicon,at->label);
        }
        at->label = 0;

        if(PConvPyObjectToStrMaxLen(label_id2,temp,sizeof(OrthoLineType)-1)) {
          if(temp[0]) {
            OVreturn_word result = OVLexicon_GetFromCString(G->Lexicon,temp);
            if(OVreturn_IS_OK(result)) {
              at->label = result.word;
            }
          }
        }
      }
      if(flags_id1!=flags_id2) {
        if(!PConvPyObjectToInt(flags_id2,&flags))
          result=false;
        else
          at->flags = flags;
      }

      if(numeric_type_id1!=numeric_type_id2) {
        if(!PConvPyObjectToInt(numeric_type_id2,&numericType))
          result=false;
        else
          at->customType = numericType;
      }
      if(ID_id1!=ID_id2) {
        if(!PConvPyObjectToInt(ID_id2,&id))
          result=false;
        else
          at->id=id;
      }
      if(rank_id1!=rank_id2) {
        if(!PConvPyObjectToInt(rank_id2,&rank))
          result=false;
        else
          at->rank=rank;
      }

      if(PyErr_Occurred()) {
        PyErr_Print();
        result=false;
      }
    }
    if(!result) { 
      ErrMessage(G,"Alter","Aborting on error. Assignment may be incomplete.");
    }
  } 
  Py_DECREF(dict);
  return(result);
}

int PLabelAtom(PyMOLGlobals *G, AtomInfoType *at,char *model,char *expr,int index)
{
  PyObject *dict;
  PyObject *P_inst_dict = G->P_inst->dict;
  int result;
  OrthoLineType label;
  char atype[7];
  OrthoLineType buffer;
  if(at->hetatm)
    strcpy(atype,"HETATM");
  else
    strcpy(atype,"ATOM");
  PBlock(G);
  /* PBlockAndUnlockAPI() is not safe.
   * what if "at" is destroyed by another thread? */
  dict = PyDict_New();

  PConvStringToPyDictItem(dict,"model",model);
  PConvIntToPyDictItem(dict,"index",index+1);
  PConvStringToPyDictItem(dict,"type",atype);
  PConvStringToPyDictItem(dict,"name",at->name);
  PConvStringToPyDictItem(dict,"resn",at->resn);
  PConvStringToPyDictItem(dict,"resi",at->resi);
  PConvIntToPyDictItem(dict,"resv",at->resv); 
  PConvStringToPyDictItem(dict,"chain",at->chain);
  PConvStringToPyDictItem(dict,"alt",at->alt);
  PConvStringToPyDictItem(dict,"segi",at->segi);
  PConvStringToPyDictItem(dict,"ss",at->ssType);
  PConvFloatToPyDictItem(dict,"vdw",at->vdw);
  PConvFloatToPyDictItem(dict,"elec_radius",at->elec_radius);
  {
    char null_st[1] = "";
    char *st = null_st;

    if(at->textType) st = OVLexicon_FetchCString(G->Lexicon,at->textType);
    PConvStringToPyDictItem(dict,"text_type",st);

    st = null_st;
    if(at->label) st = OVLexicon_FetchCString(G->Lexicon,at->label);
    PConvStringToPyDictItem(dict,"label",st);
  }
  PConvStringToPyDictItem(dict,"elem",at->elem);
  PConvIntToPyDictItem(dict,"geom",at->geom);
  PConvIntToPyDictItem(dict,"valence",at->valence);
  PConvIntToPyDictItem(dict,"rank",at->rank);
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
  else
    PConvStringToPyDictItem(dict,"numeric_type","?");  
  PConvFloatToPyDictItem(dict,"partial_charge",at->partialCharge);
  PConvIntToPyDictItem(dict,"formal_charge",at->formalCharge);
  PConvIntToPyDictItem(dict,"color",at->color);
  PConvIntToPyDictItem(dict,"cartoon",at->cartoon);
  PConvIntToPyDictItem(dict,"ID",at->id);
  PXDecRef(PyRun_String(expr,Py_single_input,P_inst_dict,dict));
  if(PyErr_Occurred()) {
    PyErr_Print();
    result=false;
  } else {
    result=true;
    if(!PConvPyObjectToStrMaxLen(PyDict_GetItemString(dict,"label"),
                                 label,sizeof(OrthoLineType)-1))
      result=false;
    if(PyErr_Occurred()) {
      PyErr_Print();
      result=false;
    }
    if(result) { 
      if(at->label) {
        OVLexicon_DecRef(G->Lexicon,at->label);
      }
      at->label = 0;
      
      if(label[0]) {
        OVreturn_word ret = OVLexicon_GetFromCString(G->Lexicon,label);
        if(OVreturn_IS_OK(ret)) {
          /*printf("alloc'd %d [%s]\n",OVLexicon_GetNActive(G->Lexicon),label);*/
          at->label = ret.word;
        }
      }
    } else {
      ErrMessage(G,"Label","Aborting on error. Labels may be incomplete.");
    }
  }
  Py_DECREF(dict);
  PUnblock(G);
  return(result);
}

void PUnlockAPIAsGlut(PyMOLGlobals *G) /* must call with unblocked interpreter */
{
  PRINTFD(G,FB_Threads)
    " PUnlockAPIAsGlut-DEBUG: entered as thread 0x%x\n",PyThread_get_thread_ident()
    ENDFD;
  PBlock(G);
  PXDecRef(PyObject_CallFunction(P_unlock,NULL)); /* NOTE this may flush the command buffer! */
  PLockStatus();
  PyMOL_PopValidContext(G->PyMOL);
  PUnlockStatus();
  PUnlockGLUT();
  PUnblock(G);
}

void PUnlockAPIAsGlutNoFlush(PyMOLGlobals *G) /* must call with unblocked interpreter */
{
  PRINTFD(G,FB_Threads)
    " PUnlockAPIAsGlut-DEBUG: entered as thread 0x%x\n",PyThread_get_thread_ident()
    ENDFD;
  PBlock(G);
  PXDecRef(PyObject_CallFunction(P_unlock,"i",-1)); /* prevents flushing of the buffer */
  PLockStatus();
  PyMOL_PopValidContext(G->PyMOL);
  PUnlockStatus();
  PUnlockGLUT();
  PUnblock(G);
}

static int get_api_lock(PyMOLGlobals *G,int block_if_busy) 
{
  int result = true;

  if(block_if_busy) {
    
    PXDecRef(PyObject_CallFunction(P_lock,NULL));

  } else { /* not blocking if PyMOL is busy */

    PyObject *got_lock = PyObject_CallFunction(P_lock_attempt,NULL);
    
    if(got_lock) {
      if(!PyInt_AsLong(got_lock)) {
        PLockStatus();
        if(PyMOL_GetBusy(G->PyMOL,false))
          result = false;
        PUnlockStatus();
        
        if(result) { /* didn't get lock, but not busy, so block and wait for lock */
          PXDecRef(PyObject_CallFunction(P_lock,NULL));
        }
      }
      Py_DECREF(got_lock);
    }
  }
  return result;
}

int PLockAPIAsGlut(PyMOLGlobals *G,int block_if_busy)
{
  PRINTFD(G,FB_Threads)
    "*PLockAPIAsGlut-DEBUG: entered as thread 0x%x\n",PyThread_get_thread_ident()
    ENDFD;

  PBlock(G);

  PLockGLUT();

  PLockStatus();
  PyMOL_PushValidContext(G->PyMOL);
  PUnlockStatus();

  PRINTFD(G,FB_Threads)
    "#PLockAPIAsGlut-DEBUG: acquiring lock as thread 0x%x\n",PyThread_get_thread_ident()
    ENDFD;
  
  if(!get_api_lock(G,block_if_busy)) {
    PLockStatus();
    PyMOL_PopValidContext(G->PyMOL);
    PUnlockStatus();
    PUnlockGLUT();
    PUnblock(G);
    return false;/* busy -- so allow main to update busy status display (if any) */
  }
   
  while(P_glut_thread_keep_out) {
    /* IMPORTANT: keeps the glut thread out of an API operation... */
    /* NOTE: the keep_out variable can only be changed or read by the thread
       holding the API lock, therefore it is safe even through increment
       isn't atomic. */
    PRINTFD(G,FB_Threads)
      "-PLockAPIAsGlut-DEBUG: glut_thread_keep_out 0x%x\n",PyThread_get_thread_ident()
      ENDFD;
    
    PXDecRef(PyObject_CallFunction(P_unlock,"i",-1)); /* prevent buffer flushing */
#ifndef WIN32
    { 
      struct timeval tv;

      PUnblock(G);
      tv.tv_sec=0;
      tv.tv_usec=50000; 
      select(0,NULL,NULL,NULL,&tv);
      PBlock(G); 
    } 
#else
/* BEGIN PROPRIETARY CODE SEGMENT (see disclaimer in "os_proprietary.h") */
    PXDecRef(PyObject_CallFunction(P_sleep,"f",0.050));
/* END PROPRIETARY CODE SEGMENT */
#endif

    if(!get_api_lock(G,block_if_busy)) {
      /* return false-- allow main to update busy status display (if any) */
      PLockStatus();
      PyMOL_PopValidContext(G->PyMOL);
      PUnlockStatus();
      PUnlockGLUT();
      PUnblock(G);
      return false;
    }
  }


  PUnblock(G); /* API is now locked, so we can free up Python...*/

  PRINTFD(G,FB_Threads)
    "=PLockAPIAsGlut-DEBUG: acquired\n"
    ENDFD;
  return true;
}

/* THESE CALLS ARE REQUIRED FOR MONOLITHIC COMPILATION TO SUCCEED UNDER WINDOWS. */
#ifndef _PYMOL_ACTIVEX
#ifndef _PYMOL_EMBEDDED
void	initExtensionClass(void);
void	initsglite(void);
void  init_champ(void);
void    init_opengl(void);
void    init_opengl_num(void);
void    init_glu(void);
void    init_glu_num(void);
void    init_glut(void);
void    initopenglutil(void);
void    initopenglutil_num(void);
#endif
#endif

#ifdef _PYMOL_MONOLITHIC
#ifndef _PYMOL_ACTIVEX
#ifndef _PYMOL_EMBEDDED
/* BEGIN PROPRIETARY CODE SEGMENT (see disclaimer in "os_proprietary.h") */
#ifdef WIN32
void	init_numpy();
void	initmultiarray();
void	initarrayfns();
void	initlapack_lite();
void	initumath();
void	initranlib();
void  init_champ();
#endif
/* END PROPRIETARY CODE SEGMENT */
#endif
#endif
#endif

#ifdef _PYMOL_MONOLITHIC
#ifndef _PYMOL_ACTIVEX
#ifndef _PYMOL_EMBEDDED
void	initExtensionClass(void);
void	initsglite(void);
void  init_champ(void);
void    init_opengl(void);
void    init_opengl_num(void);
void    init_glu(void);
void    init_glu_num(void);
void    init_glut(void);
void    initopenglutil(void);
void    initopenglutil_num(void);
#endif
#endif
#endif
/* BEGIN PROPRIETARY CODE SEGMENT (see disclaimer in "os_proprietary.h") */
#ifdef WIN32
static int IsSecurityRequired()
{
  DWORD WindowsVersion = GetVersion();
  DWORD WindowsMajorVersion = (DWORD)(LOBYTE(LOWORD(WindowsVersion)));
  DWORD WindowsMinorVersion = (DWORD)(HIBYTE(LOWORD(WindowsVersion)));

  if (WindowsVersion >= 0x80000000) return FALSE;
  
  return TRUE;
}
#endif
/* END PROPRIETARY CODE SEGMENT */

void PSetupEmbedded(PyMOLGlobals *G,int argc,char **argv)
{
  /* This routine is called if we are running with an embedded Python interpreter */
  PyObject *args, *pymol;

/* BEGIN PROPRIETARY CODE SEGMENT (see disclaimer in "os_proprietary.h") */
#ifdef WIN32
    
  { /* automatically hide the window if this process
       was started as a vanilla console application */
    HWND hwndFound;         
    if(hwndFound=FindWindow(NULL, argv[0])) {
      ShowWindow(hwndFound,SW_HIDE);
    }
  }
    
  {/* if PYMOL_PATH and/or PYTHONHOME isn't in the environment coming
      in, then the user may simply have clicked PyMOL.exe, in which
      case we need to consult the registry regarding the location of
      the install */
    
    /* embedded version of PyMOL currently ships with Python 2.4 */
#define EMBEDDED_PYTHONHOME "\\py24"
        
    OrthoLineType path_buffer;
    static char line1[8092];
    static char line2[8092];
    HKEY phkResult;
    int lpcbData;
    int lpType = REG_SZ;
    int r1,r2;
    char *pymol_path;
    char *pythonhome;
    int pythonhome_set = false;
    int restart_flag = false;
        
        
    pymol_path = getenv("PYMOL_PATH");
    pythonhome = getenv("PYTHONHOME");
    if((!pymol_path)||(!pythonhome)) {
      lpcbData = sizeof(OrthoLineType)-1;
      r1=RegOpenKeyEx(HKEY_CLASSES_ROOT,
                      "Software\\DeLano Scientific\\PyMOL\\PYMOL_PATH",
                      0,KEY_EXECUTE,&phkResult);
      if(r1!=ERROR_SUCCESS) {
        r1=RegOpenKeyEx(HKEY_CURRENT_USER,
                        "Software\\DeLano Scientific\\PyMOL\\PYMOL_PATH",
                        0,KEY_EXECUTE,&phkResult);
      }  
      if(r1==ERROR_SUCCESS) {
        r2 = RegQueryValueEx(phkResult,"",NULL,
                             &lpType,path_buffer,&lpcbData);
        if (r2==ERROR_SUCCESS) {
          /* use environment variable PYMOL_PATH first, registry entry
             second */
          if(!pymol_path) {
            strcpy(line1,"PYMOL_PATH=");
            strcat(line1,path_buffer);
            _putenv(line1);
            if(!pythonhome) { /* only set PYTHONHOME if already
                                 setting new PYMOL_PATH */
              pythonhome_set = true;
              strcpy(line2,"PYTHONHOME=");
              strcat(line2,path_buffer);
              strcat(line2,EMBEDDED_PYTHONHOME); 
              restart_flag = true;
              _putenv(line2);
            }
          }
        }
        RegCloseKey(phkResult);
      }
      /* this allows us to just specify PYMOL_PATH with no registry entries */
      if((!pythonhome_set)&&(pymol_path)&&(!pythonhome)) {
        strcpy(line2,"PYTHONHOME=");
        strcat(line2,pymol_path);
        strcat(line2,EMBEDDED_PYTHONHOME);
        _putenv(line2);
        restart_flag = true;
      }
    }
    if(restart_flag && getenv("PYMOL_PATH") && getenv("PYTHONHOME")) { 
            
      /* now that we have the environment defined, restart the process
       * so that Python can use the new environment.  If we don't do
       * this, then Python won't see the new environment vars. Why not? */

      /* note that we use CreateProcesss to launch the console
       * application instead of exec or spawn in order to hide the
       * console window. Otherwise a console window might appear, and
       * that would suck. */
            
      char command[8092];
      static char cmd_line[8092];
      char *p,*q;
      int a;
            
      /* copy arguments, installing quotes around them */
            
      sprintf(command,"%s\\pymol.exe",getenv("PYMOL_PATH"));
      p = cmd_line;
            
      sprintf(p,"\"%s\"",command);
      p+=strlen(p);
      *(p++)=' ';
      *p=0;
            
      for(a=1;a<=argc;a++) {
        q = argv[a];
        if(q) {
          if(*q!='"') { /* add quotes if not present */
            *(p++)='"';
            while(*q) {
              *(p++)=*(q++);
            }
            *(p++)='"'; 
          } else {
            while(*q) {
              *(p++)=*(q++);
            }
          }
          *(p++)=32;
          *p=0;
        }
      }

      {
        LPSECURITY_ATTRIBUTES lpSA = NULL;
        PSECURITY_DESCRIPTOR lpSD = NULL;
        STARTUPINFO si;
        PROCESS_INFORMATION pi;
        HANDLE hProcess = GetCurrentProcess();
                
        ZeroMemory(&si, sizeof(STARTUPINFO));
        si.cb = sizeof(STARTUPINFO);
        si.dwFlags = STARTF_USESHOWWINDOW;
        si.wShowWindow = SW_HIDE; 
                
        if(IsSecurityRequired())
          {
            lpSD = GlobalAlloc(GPTR, SECURITY_DESCRIPTOR_MIN_LENGTH);
            InitializeSecurityDescriptor(lpSD, SECURITY_DESCRIPTOR_REVISION);
            SetSecurityDescriptorDacl(lpSD, -1, 0, 0);
                    
            lpSA = GlobalAlloc(GPTR, sizeof(SECURITY_ATTRIBUTES));
            lpSA->nLength = sizeof(SECURITY_ATTRIBUTES);
            lpSA->lpSecurityDescriptor = lpSD;
            lpSA->bInheritHandle = TRUE;
          }
                
        if(CreateProcess(NULL, (LPTSTR)cmd_line, lpSA, NULL, TRUE,
                         0, NULL, NULL, &si, &pi)) {
          WaitForSingleObject(pi.hProcess, INFINITE);
        }
        if (lpSA != NULL) GlobalFree(lpSA);
        if (lpSD != NULL) GlobalFree(lpSD);
        _exit(0);
      }
    }
  }

#endif
/* END PROPRIETARY CODE SEGMENT */


  /* compatibility for old compile-time defines */

#ifdef _PYMOL_SETUP_PY21
#ifndef _PYMOL_SETUP_PY_EXT
#define _PYMOL_SETUP_PY_EXT
#endif
#endif
#ifdef _PYMOL_SETUP_PY22
#ifndef _PYMOL_SETUP_PY_EXT
#define _PYMOL_SETUP_PY_EXT
#endif
#endif
#ifdef _PYMOL_SETUP_PY23
#ifndef _PYMOL_SETUP_PY_EXT
#define _PYMOL_SETUP_PY_EXT
#endif
#endif
#ifdef _PYMOL_SETUP_PY24
#ifndef _PYMOL_SETUP_PY_EXT
#define _PYMOL_SETUP_PY_EXT
#endif
#endif
#ifdef _PYMOL_SETUP_PY25
#ifndef _PYMOL_SETUP_PY_EXT
#define _PYMOL_SETUP_PY_EXT
#endif
#endif

  /* should we set up PYTHONHOME in the ext directory? */

#ifdef _PYMOL_SETUP_PY_EXT
  {
    static char line1[8092];
    static char line2[8092];

    if(!getenv("PYMOL_PATH")) { /* if PYMOL_PATH isn't defined...*/
    
      /* was our startup path absolute? */
    
      if((argc>0) && (argv[0][0]=='/')) {
        /* PYMOL was started with an absolute path, so try using that... */
        char *p;
        strcpy(line1,"PYMOL_PATH=");
        strcat(line1,argv[0]);
        p=line1 + strlen(line1);
        while(*p!='/') {
          *p=0;
          p--;
        }
        *p=0;
        putenv(line1);
      } else if((argc>0) && getenv("PWD") && ((argv[0][0]=='.')||(strstr(argv[0],"/")))) { 
        /* was the path relative? */
        char *p;
        strcpy(line1,"PYMOL_PATH=");
        strcat(line1,getenv("PWD"));
        strcat(line1,"/");
        strcat(line1,argv[0]);
        p=line1 + strlen(line1);
        while(*p!='/') {
          *p=0;
          p--;
        }
        *p=0;
        putenv(line1);
      } else { /* otherwise, just try using the current working directory */
        if(getenv("PWD")) {
          strcpy(line1,"PYMOL_PATH=");
          strcat(line1,getenv("PWD"));
          putenv(line1);
        }
      }
    }
 
    /* now set PYTHONHOME so that we use the right binary libraries for
       this executable */

    if(getenv("PYMOL_PATH")) {
      strcpy(line2,"PYTHONHOME=");
      strcat(line2,getenv("PYMOL_PATH"));
      strcat(line2,"/ext");
      putenv(line2);
    }
  }
#endif

#ifndef _PYMOL_ACTIVEX
#ifndef _PYMOL_EMBEDDED
  Py_Initialize();
  PyEval_InitThreads();
  PyUnicode_SetDefaultEncoding("utf-8"); /* is this safe & legal? */

#endif
#endif

  init_cmd();
#ifdef _PYMOL_MONOLITHIC
#ifndef _PYMOL_ACTIVEX
#ifndef _PYMOL_EMBEDDED
  initExtensionClass();
  initsglite();
  /* initialize champ */
  init_champ();
/* BEGIN PROPRIETARY CODE SEGMENT (see disclaimer in "os_proprietary.h") */
#ifdef WIN32
  /* initialize numeric python */
  init_numpy();
  initmultiarray();
  initarrayfns();
  initlapack_lite();
  initumath();
  initranlib();
#endif
/* END PROPRIETARY CODE SEGMENT */
  init_opengl();
  init_opengl_num();
  init_glu();
  init_glu_num();
  init_glut();
  initopenglutil();
  initopenglutil_num();
#endif
#endif
#endif
  PyRun_SimpleString("import os\n");
  PyRun_SimpleString("import sys\n");
/* BEGIN PROPRIETARY CODE SEGMENT (see disclaimer in "os_proprietary.h") */
#ifdef WIN32
#if 0
  {
    OrthoLineType path_buffer,command;
    HKEY phkResult;
    int lpcbData;
    int lpType = REG_SZ;
    int r1,r2;
  
    lpcbData = sizeof(OrthoLineType)-1;
    r1=RegOpenKeyEx(HKEY_CLASSES_ROOT,"Software\\DeLano Scientific\\PyMOL\\PYMOL_PATH",0,KEY_EXECUTE,&phkResult);
    if(r1==ERROR_SUCCESS) {
      r2 = RegQueryValueEx(phkResult,"",NULL,&lpType,path_buffer,&lpcbData);
      if (r2==ERROR_SUCCESS) {
        /* use environment variable PYMOL_PATH first, registry entry second */
        sprintf(command,"_registry_pymol_path = r'''%s'''\n",path_buffer);
        PyRun_SimpleString(command);
        PyRun_SimpleString("if not os.environ.has_key('PYMOL_PATH'): os.environ['PYMOL_PATH']=_registry_pymol_path\n");
      }
      RegCloseKey(phkResult);
    }
  } 
#endif
  PyRun_SimpleString("if not os.environ.has_key('PYMOL_PATH'): os.environ['PYMOL_PATH']=os.getcwd()\n");
#endif
/* END PROPRIETARY CODE SEGMENT */

#ifdef _PYMOL_SETUP_TCLTK83
  /* used by semistatic pymol */
  PyRun_SimpleString("if os.path.exists(os.environ['PYMOL_PATH']+'/ext/lib/tcl8.3'): os.environ['TCL_LIBRARY']=os.environ['PYMOL_PATH']+'/ext/lib/tcl8.3'\n");
  PyRun_SimpleString("if os.path.exists(os.environ['PYMOL_PATH']+'/ext/lib/tk8.3'): os.environ['TK_LIBRARY']=os.environ['PYMOL_PATH']+'/ext/lib/tk8.3'\n");
#endif

#ifdef _PYMOL_SETUP_TCLTK84
  /* used by semistatic pymol */
  PyRun_SimpleString("if os.path.exists(os.environ['PYMOL_PATH']+'/ext/lib/tcl8.4'): os.environ['TCL_LIBRARY']=os.environ['PYMOL_PATH']+'/ext/lib/tcl8.4'\n");
  PyRun_SimpleString("if os.path.exists(os.environ['PYMOL_PATH']+'/ext/lib/tk8.4'): os.environ['TK_LIBRARY']=os.environ['PYMOL_PATH']+'/ext/lib/tk8.4'\n");
#endif

#if 0
  /* no longer necessary since we're setting PYTHONHOME */
#ifdef _PYMOL_SETUP_PY21 
  /* used by semistatic pymol */
  PyRun_SimpleString("import string");
  PyRun_SimpleString("sys.path=filter(lambda x:string.find(x,'static/ext-static')<0,sys.path)"); /* clean bogus entries in sys.path */
#endif
#ifdef _PYMOL_SETUP_PY22
  /* used by semistatic pymol */
  PyRun_SimpleString("import string");
  PyRun_SimpleString("sys.path=filter(lambda x:string.find(x,'static/ext')<0,sys.path)"); /* clean bogus entries in sys.path */
#endif
#ifdef _PYMOL_SETUP_PY23
  /* used by semistatic pymol */
  PyRun_SimpleString("import string");
  PyRun_SimpleString("sys.path=filter(lambda x:string.find(x,'static/ext')<0,sys.path)"); /* clean bogus entries in sys.path */
#endif
#endif
/* BEGIN PROPRIETARY CODE SEGMENT (see disclaimer in "os_proprietary.h") */
#ifdef WIN32
  PyRun_SimpleString("if (os.environ['PYMOL_PATH']+'/modules') not in sys.path: sys.path.insert(0,os.environ['PYMOL_PATH']+'/modules')\n");
#endif
/* END PROPRIETARY CODE SEGMENT */

  P_main = PyImport_AddModule("__main__");
  if(!P_main) ErrFatal(G,"PyMOL","can't find '__main__'");

  /* inform PyMOL's other half that we're launching embedded-style */
  PyObject_SetAttrString(P_main,"pymol_launch",PyInt_FromLong(4));

  args = PConvStringListToPyList(argc,argv); /* prepare our argument list */
  if(!args) ErrFatal(G,"PyMOL","can't process arguments.");

  /* copy arguments to __main__.pymol_argv */
  PyObject_SetAttrString(P_main,"pymol_argv",args);
  
  PyRun_SimpleString("if (os.environ['PYMOL_PATH']+'/modules') not in sys.path: sys.path.insert(0,os.environ['PYMOL_PATH']+'/modules')\n"); /* needed for semistatic pymol */

  PyRun_SimpleString("import pymol"); /* create the global PyMOL namespace */

  pymol = PyImport_AddModule("pymol"); /* get it */
  if(!pymol) ErrFatal(G,"PyMOL","can't find module 'pymol'");

}

void PGetOptions(CPyMOLOptions *rec)
{
  PyObject *pymol,*invocation,*options;
  char *load_str;

  pymol = PyImport_AddModule("pymol"); /* get it */
  if(!pymol) {fprintf(stderr,"PyMOL-ERROR: can't find module 'pymol'"); exit(EXIT_FAILURE);}

  invocation = PyObject_GetAttrString(pymol,"invocation"); /* get a handle to the invocation module */
  if(!invocation) {fprintf(stderr,"PyMOL-ERROR: can't find module 'invocation'"); exit(EXIT_FAILURE);}

  options = PyObject_GetAttrString(invocation,"options");
  if(!options) {fprintf(stderr,"PyMOL-ERROR: can't get 'invocation.options'.");exit(EXIT_FAILURE);}

  rec->pmgui = ! PyInt_AsLong(PyObject_GetAttrString(options,"no_gui"));
  rec->internal_gui = PyInt_AsLong(PyObject_GetAttrString(options,"internal_gui"));
  rec->internal_feedback = PyInt_AsLong(PyObject_GetAttrString(options,"internal_feedback"));
  rec->show_splash = PyInt_AsLong(PyObject_GetAttrString(options,"show_splash"));
  rec->security = PyInt_AsLong(PyObject_GetAttrString(options,"security"));
  rec->game_mode = PyInt_AsLong(PyObject_GetAttrString(options,"game_mode"));
  rec->force_stereo = PyInt_AsLong(PyObject_GetAttrString(options,"force_stereo"));
  rec->winX = PyInt_AsLong(PyObject_GetAttrString(options,"win_x"));
  rec->winY = PyInt_AsLong(PyObject_GetAttrString(options,"win_y"));
  rec->winPX = PyInt_AsLong(PyObject_GetAttrString(options,"win_px"));
  rec->winPY = PyInt_AsLong(PyObject_GetAttrString(options,"win_py"));
  rec->blue_line = PyInt_AsLong(PyObject_GetAttrString(options,"blue_line"));
  rec->external_gui = PyInt_AsLong(PyObject_GetAttrString(options,"external_gui"));
  rec->siginthand = PyInt_AsLong(PyObject_GetAttrString(options,"sigint_handler"));
  rec->reuse_helper = PyInt_AsLong(PyObject_GetAttrString(options,"reuse_helper"));
  rec->auto_reinitialize = PyInt_AsLong(PyObject_GetAttrString(options,"auto_reinitialize"));
  rec->keep_thread_alive = PyInt_AsLong(PyObject_GetAttrString(options,"keep_thread_alive"));
  rec->quiet = PyInt_AsLong(PyObject_GetAttrString(options,"quiet"));
#ifdef _IPYMOL
  rec->incentive_product = true;
#else
  rec->incentive_product = PyInt_AsLong(PyObject_GetAttrString(options,"incentive_product"));
#endif
  rec->multisample = PyInt_AsLong(PyObject_GetAttrString(options,"multisample"));
  rec->window_visible = PyInt_AsLong(PyObject_GetAttrString(options,"window_visible"));
  rec->read_stdin = PyInt_AsLong(PyObject_GetAttrString(options,"read_stdin"));
  rec->presentation = PyInt_AsLong(PyObject_GetAttrString(options,"presentation"));
  rec->defer_builds_mode = PyInt_AsLong(PyObject_GetAttrString(options,"defer_builds_mode"));
  rec->full_screen = PyInt_AsLong(PyObject_GetAttrString(options,"full_screen"));
  load_str = PyString_AsString(PyObject_GetAttrString(options,"after_load_script"));
  rec->sphere_mode = PyInt_AsLong(PyObject_GetAttrString(options,"sphere_mode"));
  rec->stereo_capable = PyInt_AsLong(PyObject_GetAttrString(options,"stereo_capable"));
  rec->passive_stereo= PyInt_AsLong(PyObject_GetAttrString(options,"passive_stereo"));
  rec->zoom_mode= PyInt_AsLong(PyObject_GetAttrString(options,"zoom_mode"));
  if(load_str) {
    if(load_str[0]) {
      UtilNCopy(rec->after_load_script,load_str,PYMOL_MAX_OPT_STR);
    }
  }
  if(PyErr_Occurred()) {
    PyErr_Print();
  }
}

void PRunStringModule(PyMOLGlobals *G,char *str) /* runs a string in the namespace of the pymol instance */
{
  PXDecRef(PyObject_CallFunction(G->P_inst->exec,"Os",P_pymol,str));
}

void PRunStringInstance(PyMOLGlobals *G,char *str) /* runs a string in the namespace of the pymol instance */
{
  PXDecRef(PyObject_CallFunction(G->P_inst->exec,"Os",G->P_inst->obj,str));
}

void PInit(PyMOLGlobals *G,int global_instance) 
{
  PyObject *sys,*pcatch;
  int a;

#ifdef PYMOL_NEW_THREADS
   PyEval_InitThreads();
#endif

/* BEGIN PROPRIETARY CODE SEGMENT (see disclaimer in "os_proprietary.h") */
#ifdef WIN32
#ifdef _PYMOL_MONOLITHIC
#ifndef _PYMOL_ACTIVEX
#ifndef _PYMOL_EMBEDDED
#define _PYMOL_INIT_MODULES
#endif
#endif
#endif
#endif
/* END PROPRIETARY CODE SEGMENT */

#ifdef _PYMOL_INIT_MODULES
	/* Win32 module build: includes pyopengl, numpy, and sglite */
	/* sglite */
	initExtensionClass();
	initsglite();
    init_champ();
	/* initialize numeric python */
	init_numpy();
	initmultiarray();
	initarrayfns();
	initlapack_lite();
	initumath();
	initranlib();
	/* initialize PyOpenGL */
    init_opengl();
    init_opengl_num();
    init_glu();
    init_glu_num();
    init_glut();
    initopenglutil();
	initopenglutil_num();
#endif


  for(a=0;a<MAX_SAVED_THREAD;a++) {
    SavedThread[a].id=-1;
  }

  if(global_instance) {
    PCatchInit();   /* setup standard-output catch routine */
  }

/* assumes that pymol module has been loaded via Python or PyRun_SimpleString */

  P_pymol = PyImport_AddModule("pymol"); /* get it */
  if(!P_pymol) ErrFatal(G,"PyMOL","can't find module 'pymol'");
  P_pymol_dict = PyModule_GetDict(P_pymol);
  if(!P_pymol_dict) ErrFatal(G,"PyMOL","can't find globals for 'pymol'");

  if(global_instance) { /* if global singleton PyMOL... */
    G->P_inst = Calloc(CP_inst,1);
    G->P_inst->obj = P_pymol;
    G->P_inst->dict = P_pymol_dict;
  }
  
  {
    G->P_inst->exec = PyDict_GetItemString(P_pymol_dict,"exec_str");
    if(!G->P_inst->exec) ErrFatal(G,"PyMOL","can't find 'pymol.exec_str()'");

    sys = PyDict_GetItemString(P_pymol_dict,"sys");
    if(!sys) ErrFatal(G,"PyMOL","can't find 'pymol.sys'");

    if(global_instance) { 

      /* implies global singleton pymol, so set up the global handle */
      PyDict_SetItemString(P_pymol_dict,"_COb",PyCObject_FromVoidPtr((void*)&TempPyMOLGlobals,NULL));

      pcatch = PyImport_AddModule("pcatch"); 
      if(!pcatch) ErrFatal(G,"PyMOL","can't find module 'pcatch'");
      PyObject_SetAttrString(sys,"stdout",pcatch);
      PyObject_SetAttrString(sys,"stderr",pcatch);
    }

    PRunStringModule(G,"import traceback\n");  
    P_traceback = PyDict_GetItemString(P_pymol_dict,"traceback");
    if(!P_traceback) ErrFatal(G,"PyMOL","can't find 'traceback'");

    PRunStringModule(G,"import cmd\n");  
    P_cmd = PyDict_GetItemString(P_pymol_dict,"cmd");
    if(!P_cmd) ErrFatal(G,"PyMOL","can't find 'cmd'");

    if(global_instance) { 
      /* implies global singleton pymol, so set up the global handle */
      PyObject_SetAttrString(P_cmd,"_COb",PyCObject_FromVoidPtr((void*)&TempPyMOLGlobals,NULL));

      /* cmd module is itself the api for the global PyMOL instance */
      G->P_inst->cmd = P_cmd;
    }

    PyObject_SetAttrString(G->P_inst->cmd,"_pymol",G->P_inst->obj);

    /* right now, all locks are global -- eventually some of these may
       become instance-specific in order to improve concurrency */

    P_lock = PyObject_GetAttrString(P_cmd,"lock");
    if(!P_lock) ErrFatal(G,"PyMOL","can't find 'cmd.lock()'");

    P_lock_attempt = PyObject_GetAttrString(P_cmd,"lock_attempt");
    if(!P_lock_attempt) ErrFatal(G,"PyMOL","can't find 'cmd.lock_attempt()'");

    P_unlock = PyObject_GetAttrString(P_cmd,"unlock");
    if(!P_unlock) ErrFatal(G,"PyMOL","can't find 'cmd.unlock()'");

    P_lock_c = PyObject_GetAttrString(P_cmd,"lock_c");
    if(!P_lock_c) ErrFatal(G,"PyMOL","can't find 'cmd.lock_c()'");

    P_unlock_c = PyObject_GetAttrString(P_cmd,"unlock_c");
    if(!P_unlock_c) ErrFatal(G,"PyMOL","can't find 'cmd.unlock_c()'");

    P_lock_status = PyObject_GetAttrString(P_cmd,"lock_status");
    if(!P_lock_status) ErrFatal(G,"PyMOL","can't find 'cmd.lock_status()'");

    P_lock_status_attempt = PyObject_GetAttrString(P_cmd,"lock_status_attempt");
    if(!P_lock_status_attempt) ErrFatal(G,"PyMOL","can't find 'cmd.lock_status_attempt()'");

    P_unlock_status = PyObject_GetAttrString(P_cmd,"unlock_status");
    if(!P_unlock_status) ErrFatal(G,"PyMOL","can't find 'cmd.unlock_status()'");

    P_lock_glut = PyObject_GetAttrString(P_cmd,"lock_glut");
    if(!P_lock_glut) ErrFatal(G,"PyMOL","can't find 'cmd.lock_glut()'");

    P_unlock_glut = PyObject_GetAttrString(P_cmd,"unlock_glut");
    if(!P_unlock_glut) ErrFatal(G,"PyMOL","can't find 'cmd.unlock_glut()'");
    
    /* 'do' command */

    G->P_inst->cmd_do = PyObject_GetAttrString(G->P_inst->cmd,"do");
    if(!G->P_inst->cmd_do) ErrFatal(G,"PyMOL","can't find 'cmd.do()'");

    /* invariant stuff */

    PRunStringModule(G,"import menu\n");  
    P_menu = PyDict_GetItemString(P_pymol_dict,"menu");
    if(!P_menu) ErrFatal(G,"PyMOL","can't find module 'menu'");

    PRunStringModule(G,"import setting\n");  
    P_setting = PyDict_GetItemString(P_pymol_dict,"setting");
    if(!P_setting) ErrFatal(G,"PyMOL","can't find module 'setting'");

    PRunStringModule(G,"import povray\n");  
    P_povray = PyDict_GetItemString(P_pymol_dict,"povray");
    if(!P_povray) ErrFatal(G,"PyMOL","can't find module 'povray'");

#ifdef _PYMOL_XRAY
    PRunStringModule(G,"import xray\n");  
    P_xray = PyDict_GetItemString(P_pymol_dict,"xray");
    if(!P_xray) ErrFatal(G,"PyMOL","can't find module 'xray'");
#endif

    /* BEGIN PROPRIETARY CODE SEGMENT (see disclaimer in "os_proprietary.h") */
#ifdef WIN32
    PRunStringModule(G,"import time\n");  
    P_time = PyDict_GetItemString(P_pymol_dict,"time");
    if(!P_time) ErrFatal(G,"PyMOL","can't find module 'time'");

    P_sleep = PyObject_GetAttrString(P_time,"sleep");
    if(!P_sleep) ErrFatal(G,"PyMOL","can't find 'time.sleep()'");
#endif
    /* END PROPRIETARY CODE SEGMENT */

    PRunStringModule(G,"import parser\n");  
    P_parser = PyDict_GetItemString(P_pymol_dict,"parser");
    if(!P_parser) ErrFatal(G,"PyMOL","can't find module 'parser'");

    {
      PyObject *fn_closure = PyObject_GetAttrString(P_parser,"new_parse_closure");
      G->P_inst->parse = PyObject_CallFunction(fn_closure,"O",G->P_inst->cmd);
      PXDecRef(fn_closure);
      if(!G->P_inst->parse) ErrFatal(G,"PyMOL","can't create 'parse' function closure");
    }

    {
      PyObject *fn_closure = PyObject_GetAttrString(P_parser,"new_complete_closure");
      G->P_inst->complete = PyObject_CallFunction(fn_closure,"O",G->P_inst->cmd);
      PXDecRef(fn_closure);
      if(!G->P_inst->complete) ErrFatal(G,"PyMOL","can't create 'complete' function closure");
    }

    PRunStringModule(G,"import chempy"); 
    P_chempy = PyDict_GetItemString(P_pymol_dict,"chempy");
    if(!P_chempy) ErrFatal(G,"PyMOL","can't find 'chempy'");

    PRunStringModule(G,"from chempy.bonds import bonds"); /* load bond dictionary */

    PRunStringModule(G,"from chempy import models"); 
    P_models = PyDict_GetItemString(P_pymol_dict,"models");
    if(!P_models) ErrFatal(G,"PyMOL","can't find 'chempy.models'");

    PRunStringModule(G,"import util\n");  
    PRunStringModule(G,"import preset\n");  
    PRunStringModule(G,"import contrib\n");

    PRunStringModule(G,"import string\n"); 

    /* backwards compatibility */

    PRunStringModule(G,"pm = cmd\n");  
    PRunStringModule(G,"pmu = util\n");  

    PRunStringModule(G,"glutThread = thread.get_ident()");

    P_glut_thread_id = PyThread_get_thread_ident();

#ifndef WIN32
    if(G->Option->siginthand) {
      signal(SIGINT,my_interrupt);
    }
#endif

    /* required environment variables */

    PyRun_SimpleString("import os");
    PyRun_SimpleString(
                       "if not os.environ.has_key('PYMOL_DATA'): os.environ['PYMOL_DATA']=os.environ['PYMOL_PATH']+'/data'");
    PyRun_SimpleString(
                       "os.environ['TUT']=os.environ['PYMOL_DATA']+'/tut'");

    PyRun_SimpleString(
                       "if not os.environ.has_key('PYMOL_SCRIPTS'): os.environ['PYMOL_SCRIPTS']=os.environ['PYMOL_PATH']+'/scripts'");

  }
}

int PPovrayRender(PyMOLGlobals *G,char *header,char *inp,char *file,int width,int height,int antialias) 
{
  PyObject *result;
  int ok;
  PBlock(G);
  result = PyObject_CallMethod(P_povray,"render_from_string","sssiii",header,inp,file,width,height,antialias);
  ok = PyObject_IsTrue(result);
  Py_DECREF(result);
  PUnblock(G);
  return(ok);
}

void PSGIStereo(PyMOLGlobals *G,int flag) 
{
  int blocked;
  blocked = PAutoBlock(G);
  if(flag) 
    PRunStringModule(G,"cmd._sgi_stereo(1)");
  else
    PRunStringModule(G,"cmd._sgi_stereo(0)");
  if(blocked) PUnblock(G);
}

void PFree(void)
{
}

void PExit(PyMOLGlobals *G,int code)
{
  ExecutiveDelete(G,"all");
  PBlock(G);
#ifndef _PYMOL_NO_MAIN
  MainFree();
#endif

  /* we're having trouble with threading errors after calling Py_Exit,
     so for the time being, let's just take the process down at this
     point, instead of allowing PyExit to be called. */

  exit(EXIT_SUCCESS);                        

  Py_Exit(code);
}

void PParse(PyMOLGlobals *G,char *str) 
{
  OrthoCommandIn(G,str);
}

void PDo(PyMOLGlobals *G,char *str) /* assumes we already hold the re-entrant API lock */
{
  int blocked;
  blocked = PAutoBlock(G);
  Py_XDECREF(PyObject_CallFunction(G->P_inst->cmd_do,"s",str));
  PAutoUnblock(G,blocked);
}

void PLog(PyMOLGlobals *G,char *str,int format) 
     /* general log routine can write PML 
        or PYM commands to appropriate log file */
{  
  int mode;
  int a;
  int blocked;
  PyObject *log;
  OrthoLineType buffer="";
  mode = (int)SettingGet(G,cSetting_logging);
  if(mode)
    {
      blocked = PAutoBlock(G);
      log = PyDict_GetItemString(P_pymol_dict,P_log_file_str);
      if(log&&(log!=Py_None)) {
        if(format==cPLog_no_flush) {
          PyObject_CallMethod(log,"write","s",str); /* maximize responsiveness (for real-time) */
        } else {
          switch(mode) {
          case cPLog_pml: /* .pml file */
            switch(format) {
            case cPLog_pml_lf:
              strcpy(buffer,str);
              break;
            case cPLog_pml:
            case cPLog_pym:
              strcpy(buffer,str);
              strcat(buffer,"\n");
              break;
            }
            break;
          case cPLog_pym: /* .pym file */
            if((str[0]=='_')&&(str[1])==' ')
              str+=2;
            switch(format) {
            case cPLog_pml_lf:
              a =strlen(str);
              while(a) { /* trim CR/LF etc. */
                if(*(str+a)>=32) break;
                *(str+a)=0;
                a--;
              }
            case cPLog_pml:
              strcpy(buffer,"cmd.do('''");
              strcat(buffer,str);
              strcat(buffer,"''')\n");
              break;
            case cPLog_pym:
              strcpy(buffer,str);
              strcat(buffer,"\n");
              break;
            }
          }
          PyObject_CallMethod(log,"write","s",buffer);        
          PyObject_CallMethod(log,"flush","");
        }
      }
      PAutoUnblock(G,blocked);
    }
}

void PLogFlush(PyMOLGlobals *G)
{
  int mode;
  PyObject *log;
  int blocked;
  mode = (int)SettingGet(G,cSetting_logging);
  if(mode)
    {
      blocked = PAutoBlock(G);
      log = PyDict_GetItemString(P_pymol_dict,P_log_file_str);
      if(log&&(log!=Py_None)) {
        PyObject_CallMethod(log,"flush","");
      }
      PAutoUnblock(G,blocked);
    }
}

void PFlush(PyMOLGlobals *G) {  
  /* NOTE: ASSUMES unblocked Python threads and a locked API */
  PyObject *err;
  char buffer[OrthoLineLength+1];
  while(OrthoCommandOut(G,buffer)) {
    PBlockAndUnlockAPI(G);
    
    PXDecRef(PyObject_CallFunction(G->P_inst->parse,"si",buffer,0));
    err = PyErr_Occurred();
    if(err) {
      PyErr_Print();
      PRINTFB(G,FB_Python,FB_Errors)
        " PFlush: Uncaught exception.  PyMOL may have a bug.\n"
        ENDFB(G);
    }
    PLockAPIAndUnblock(G);
  }
}

void PFlushFast(PyMOLGlobals *G) {
  /* NOTE: ASSUMES we currently have blocked Python threads and an unlocked API */ 
  PyObject *err;
  char buffer[OrthoLineLength+1];
  while(OrthoCommandOut(G,buffer)) {
    PRINTFD(G,FB_Threads)
      " PFlushFast-DEBUG: executing '%s' as thread 0x%x\n",buffer,
      PyThread_get_thread_ident()
      ENDFD;
    PXDecRef(PyObject_CallFunction(G->P_inst->parse,"si",buffer,0));
    err = PyErr_Occurred();
    if(err) {
      PyErr_Print();
      PRINTFB(G,FB_Python,FB_Errors)
        " PFlushFast: Uncaught exception.  PyMOL may have a bug.\n"
        ENDFB(G);
    }
  }
}


void PBlockLegacy()
{
  PBlock(TempPyMOLGlobals);
}

void PUnblockLegacy()
{
  PUnblock(TempPyMOLGlobals);
}

void PBlock(PyMOLGlobals *G)
{

  if(!PAutoBlock(G)) {
    ErrFatal(G,"PBlock","Threading error detected.  Terminating...");
  }
}


int PAutoBlock(PyMOLGlobals *G)
{
#ifndef _PYMOL_ACTIVEX
#ifndef _PYMOL_EMBEDDED
  int a,id;
  /* synchronize python */

  id = PyThread_get_thread_ident();
  PRINTFD(G,FB_Threads)
	 " PAutoBlock-DEBUG: search 0x%x (0x%x, 0x%x, 0x%x)\n",id,
	 SavedThread[MAX_SAVED_THREAD-1].id,
	 SavedThread[MAX_SAVED_THREAD-2].id,
	 SavedThread[MAX_SAVED_THREAD-3].id
	 ENDFD;
  a = MAX_SAVED_THREAD-1;
  while(a) {
    if(!((SavedThread+a)->id-id)) { 
      /* astoundingly, equality test fails on ALPHA even 
       * though the ints are equal. Must be some kind of optimizer bug
       * or mis-assumption */
      
      PRINTFD(G,FB_Threads)
        " PAutoBlock-DEBUG: seeking global lock 0x%x\n",id
      ENDFD;

#ifdef PYMOL_NEW_THREADS

      PyEval_AcquireLock();

      PRINTFD(G,FB_Threads)
        " PAutoBlock-DEBUG: restoring 0x%x\n",id
      ENDFD;
      
      PyThreadState_Swap((SavedThread+a)->state);

#else
      PRINTFD(G,FB_Threads)
        " PAutoBlock-DEBUG: restoring 0x%x\n",id
      ENDFD;
      
      PyEval_RestoreThread((SavedThread+a)->state);
#endif
      
      PRINTFD(G,FB_Threads)
        " PAutoBlock-DEBUG: restored 0x%x\n",id
      ENDFD;

      PRINTFD(G,FB_Threads)
        " PAutoBlock-DEBUG: clearing 0x%x\n",id
      ENDFD;

      PXDecRef(PyObject_CallFunction(P_lock_c,NULL));
      SavedThread[a].id = -1; 
      /* this is the only safe time we can change things */
      PXDecRef(PyObject_CallFunction(P_unlock_c,NULL));
      
      PRINTFD(G,FB_Threads)
        " PAutoBlock-DEBUG: blocked 0x%x (0x%x, 0x%x, 0x%x)\n",PyThread_get_thread_ident(),
        SavedThread[MAX_SAVED_THREAD-1].id,
        SavedThread[MAX_SAVED_THREAD-2].id,
        SavedThread[MAX_SAVED_THREAD-3].id
        ENDFD;

      return 1;
    }
    a--;
  }
  PRINTFD(G,FB_Threads)
    " PAutoBlock-DEBUG: 0x%x not found, thus already blocked.\n",PyThread_get_thread_ident()
    ENDFD;
  return 0;
#else
  return 1;
#endif
#else
  return 1;
#endif
}

int PIsGlutThread(void)
{
  return(PyThread_get_thread_ident()==P_glut_thread_id);
}

void PUnblock(PyMOLGlobals *G)
{
#ifndef _PYMOL_ACTIVEX
#ifndef _PYMOL_EMBEDDED
  int a;
  /* NOTE: ASSUMES a locked API */

  PRINTFD(G,FB_Threads)
    " PUnblock-DEBUG: entered as thread 0x%x\n",PyThread_get_thread_ident()
    ENDFD;

  /* reserve a space while we have a lock */
  PXDecRef(PyObject_CallFunction(P_lock_c,NULL));
  a = MAX_SAVED_THREAD-1;
  while(a) {
    if((SavedThread+a)->id == -1 ) {
      (SavedThread+a)->id = PyThread_get_thread_ident();
#ifdef PYMOL_NEW_THREADS
      (SavedThread+a)->state = PyThreadState_Get();
#endif
      break;
    }
    a--;
  }
  PRINTFD(G,FB_Threads)
    " PUnblock-DEBUG: 0x%x stored in slot %d\n",(SavedThread+a)->id,a
    ENDFD;
  PXDecRef(PyObject_CallFunction(P_unlock_c,NULL));
#ifdef PYMOL_NEW_THREADS
  PyThreadState_Swap(NULL);
  PyEval_ReleaseLock();
#else
  (SavedThread+a)->state = PyEval_SaveThread();  
#endif
  
#endif
#endif

}


void PAutoUnblock(PyMOLGlobals *G,int flag)
{
  if(flag) PUnblock(G);
}

void PBlockAndUnlockAPI(PyMOLGlobals *G)
{
  PBlock(G);
  PXDecRef(PyObject_CallFunction(P_unlock,NULL));
}

void PLockAPIAndUnblock(PyMOLGlobals *G)
{
  PXDecRef(PyObject_CallFunction(P_lock,NULL));
  PUnblock(G);
}

void PDefineFloat(PyMOLGlobals *G,char *name,float value) {
  char buffer[OrthoLineLength];
  sprintf(buffer,"%s = %f\n",name,value);
  PBlock(G);
  PRunStringModule(G,buffer);
  PUnblock(G);
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
    if(TempPyMOLGlobals) {
      if(Feedback(TempPyMOLGlobals,FB_Python,FB_Output)) {
        OrthoAddOutput(TempPyMOLGlobals,str);
      }
    }
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
#else
typedef int this_source_file_is_not_empty;
#endif


