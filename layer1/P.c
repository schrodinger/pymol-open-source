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

#include"os_predef.h"
#ifdef WIN32
#include<windows.h>
#endif
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

PyObject *P_globals = NULL;

PyObject *P_cmd = NULL;
PyObject *P_menu = NULL;
PyObject *P_xray = NULL;
PyObject *P_parser = NULL;
PyObject *P_setting = NULL;
PyObject *P_povray = NULL;
PyObject *P_traceback = NULL;

PyObject *P_chempy = NULL;
PyObject *P_models = NULL;

PyObject *P_complete = NULL;

PyObject *P_exec = NULL;
PyObject *P_parse = NULL;
PyObject *P_lock = NULL; /* API locks */
PyObject *P_unlock = NULL;

PyObject *P_lock_c = NULL; /* C locks */
PyObject *P_unlock_c = NULL;

PyObject *P_time = NULL;
PyObject *P_sleep = NULL;
PyObject *P_main = NULL;
PyObject *P_vfont = NULL;
PyObject *P_embed = NULL;

#define P_log_file_str "_log_file"

#define xxxPYMOL_NEW_THREADS


unsigned int PyThread_get_thread_ident(void); /* critical functionality */

typedef struct {
  int id;
  PyThreadState *state;
} SavedThreadRec;

#define MAX_SAVED_THREAD 16

static SavedThreadRec SavedThread[MAX_SAVED_THREAD];

int P_glut_thread_keep_out = 0; 
unsigned int P_glut_thread_id;

/* enables us to keep glut out if by chance it grabs the API
 * in the middle of a nested API based operation */

void PCatchInit(void);
void my_interrupt(int a);
char *getprogramname(void);

PyObject *GetBondsDict(void)
{
  PyObject *result = NULL;
  result = PyObject_GetAttrString(P_chempy,"bonds");
  if(!result) ErrMessage("PyMOL","can't find 'chempy.bonds.bonds'");
  return(result);
}

PyObject *PGetFontDict(float size,int face,int style)
{ /* assumes we have a valid interpreter lock */
  PyObject *result = NULL; 

  if(!P_vfont) {
    PRunString("import vfont\n");  
    P_vfont = PyDict_GetItemString(P_globals,"vfont");
  }
  if(!P_vfont) {
    PRINTFB(FB_Python,FB_Errors)
      " PyMOL-Error: can't find module 'vfont'"
      ENDFB;
  }
  else {
    result = PyObject_CallMethod(P_vfont,"get_font","fii",size,face,style);
  }
  return(PConvAutoNone(result));
}

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
  PRINTFD(FB_Threads)
    " PSleep-DEBUG: napping.\n"
  ENDFD;
  tv.tv_sec=0;
  tv.tv_usec=usec; 
  select(0,NULL,NULL,NULL,&tv);
  PRINTFD(FB_Threads)
    " PSleep-DEBUG: nap over.\n"
  ENDFD;
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
  PyObject_CallMethod(P_traceback,"print_tb","O",err);
}

void PDumpException()
{
  PyObject_CallMethod(P_traceback,"print_exc","");
}

int PAlterAtomState(float *v,char *expr,int read_only,AtomInfoType *at) 
     /* assumes Blocked python interpreter*/
{
  PyObject *dict; 
  int result=true;
  float f[3];
  PyObject *x_id1,*x_id2=NULL,*y_id1,*y_id2=NULL,*z_id1,*z_id2=NULL;
  char atype[7];
  dict = PyDict_New();

  if(at) {
    if(at->hetatm)
      strcpy(atype,"HETATM");
    else
      strcpy(atype,"ATOM");
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
    PConvStringToPyDictItem(dict,"text_type",at->textType);
    PConvIntToPyDictItem(dict,"numeric_type",at->customType);
    PConvFloatToPyDictItem(dict,"q",at->q);
    PConvFloatToPyDictItem(dict,"b",at->b);
    PConvFloatToPyDictItem(dict,"vdw",at->vdw);
    PConvFloatToPyDictItem(dict,"partial_charge",at->partialCharge);
    PConvIntToPyDictItem(dict,"formal_charge",at->formalCharge);
    PConvIntToPyDictItem(dict,"cartoon",at->cartoon);
    PConvStringToPyDictItem(dict,"label",at->label);
    PConvIntToPyDictItem(dict,"color",at->color);
    PConvIntToPyDictItem(dict,"ID",at->id);
  }
  x_id1 = PConvFloatToPyDictItem(dict,"x",v[0]);
  y_id1 = PConvFloatToPyDictItem(dict,"y",v[1]);
  z_id1 = PConvFloatToPyDictItem(dict,"z",v[2]);
  PyRun_String(expr,Py_single_input,P_globals,dict);
  if(PyErr_Occurred()) {
    PyErr_Print();
    result=false;
  } else if(!read_only) {
    if(result) {
      if(!(x_id2 = PyDict_GetItemString(dict,"x")))
        result=false;
      if(!(y_id2 = PyDict_GetItemString(dict,"y")))
        result=false;
      if(!(z_id2 = PyDict_GetItemString(dict,"z")))
        result=false;
      if(PyErr_Occurred()) {
        PyErr_Print();
        result=false;
        ErrMessage("AlterState","Aborting on error. Assignment may be incomplete.");
      }
    }
    if(result) {
      f[0]=(float)PyFloat_AsDouble(x_id2);
      f[1]=(float)PyFloat_AsDouble(y_id2);
      f[2]=(float)PyFloat_AsDouble(z_id2);
      if(PyErr_Occurred()) {
        PyErr_Print();
        result=false;
        ErrMessage("AlterState","Aborting on error. Assignment may be incomplete.");
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

int PAlterAtom(AtomInfoType *at,char *expr,int read_only,char *model,int index)
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
  PyObject *segi_id1,*segi_id2=NULL;
  TextType textType;
  PyObject *text_type_id1,*text_type_id2=NULL;
  SSType ssType;
  PyObject *ss_id1,*ss_id2=NULL;
  char atype[7];
  PyObject *type_id1,*type_id2=NULL;
  float b,q,partialCharge,vdw;
  PyObject *b_id1,*b_id2=NULL;
  PyObject *q_id1,*q_id2=NULL;
  PyObject *partial_charge_id1,*partial_charge_id2=NULL;
  PyObject *vdw_id1,*vdw_id2=NULL;
  int formalCharge,numericType;
  PyObject *formal_charge_id1,*formal_charge_id2=NULL;
  PyObject *numeric_type_id1,*numeric_type_id2=NULL;
  int cartoon;
  PyObject *cartoon_id1,*cartoon_id2=NULL;
  int color;
  PyObject *color_id1,*color_id2=NULL;
  PyObject *label_id1,*label_id2=NULL;
  LabelType label;
  int id;
  PyObject *ID_id1,*ID_id2=NULL;
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
  resi_id1 = PConvStringToPyDictItem(dict,"resi",at->resi);
  resv_id1 = PConvIntToPyDictItem(dict,"resv",at->resv); /* subordinate to resi */
  chain_id1 = PConvStringToPyDictItem(dict,"chain",at->chain);
  alt_id1 = PConvStringToPyDictItem(dict,"alt",at->alt);
  segi_id1 = PConvStringToPyDictItem(dict,"segi",at->segi);
  elem_id1 = PConvStringToPyDictItem(dict,"elem",at->elem);
  ss_id1 = PConvStringToPyDictItem(dict,"ss",at->ssType);
  text_type_id1 = PConvStringToPyDictItem(dict,"text_type",at->textType);
  numeric_type_id1 = PConvIntToPyDictItem(dict,"numeric_type",at->customType);
  q_id1 = PConvFloatToPyDictItem(dict,"q",at->q);
  b_id1 = PConvFloatToPyDictItem(dict,"b",at->b);
  vdw_id1 = PConvFloatToPyDictItem(dict,"vdw",at->vdw);
  partial_charge_id1 = PConvFloatToPyDictItem(dict,"partial_charge",at->partialCharge);
  formal_charge_id1 = PConvIntToPyDictItem(dict,"formal_charge",at->formalCharge);
  cartoon_id1 = PConvIntToPyDictItem(dict,"cartoon",at->cartoon);
  label_id1 = PConvStringToPyDictItem(dict,"label",at->label);
  color_id1 = PConvIntToPyDictItem(dict,"color",at->color);
  ID_id1 = PConvIntToPyDictItem(dict,"ID",at->id);

  PyRun_String(expr,Py_single_input,P_globals,dict);
  if(PyErr_Occurred()) {
    ErrMessage("Alter","Aborting on error. Assignment may be incomplete.");
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
          AtomInfoAssignParameters(at);
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
            if(!sscanf(resi,"%i",&at->resv))
              at->resv=1;
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
        if(!PConvPyObjectToStrMaxLen(text_type_id2,textType,sizeof(TextType)-1))
          result=false;
        else
          strcpy(at->textType,textType);
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
      if(partial_charge_id1!=partial_charge_id2) {
        if(!PConvPyObjectToFloat(partial_charge_id2,&partialCharge))
          result=false;
        else
          at->partialCharge=partialCharge;

      }
      if(formal_charge_id1!=formal_charge_id2) {
        if(!PConvPyObjectToInt(formal_charge_id2,&formalCharge))
          result=false;
        else
          at->formalCharge=formalCharge;

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
      if(label_id1!=label_id2) {
        if(!PConvPyObjectToStrMaxLen(label_id2,label,sizeof(LabelType)-1))
          result=false;
        else {
          strcpy(at->label,label);
        }
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

      if(PyErr_Occurred()) {
        PyErr_Print();
        result=false;
      }
    }
    if(!result) { 
      ErrMessage("Alter","Aborting on error. Assignment may be incomplete.");
    }
  } 
  Py_DECREF(dict);
  return(result);
}

int PLabelAtom(AtomInfoType *at,char *expr,int index)
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
  PBlock();
  /* PBlockAndUnlockAPI() is not safe.
   * what if "at" is destroyed by another thread? */
  dict = PyDict_New();

  PConvIntToPyDictItem(dict,"index",index+1);
  PConvStringToPyDictItem(dict,"type",atype);
  PConvStringToPyDictItem(dict,"name",at->name);
  PConvStringToPyDictItem(dict,"resn",at->resn);
  PConvStringToPyDictItem(dict,"resi",at->resi);
  PConvStringToPyDictItem(dict,"chain",at->chain);
  PConvStringToPyDictItem(dict,"alt",at->alt);
  PConvStringToPyDictItem(dict,"segi",at->segi);
  PConvStringToPyDictItem(dict,"ss",at->ssType);
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
  else
    PConvStringToPyDictItem(dict,"numeric_type","?");  
  PConvFloatToPyDictItem(dict,"partial_charge",at->partialCharge);
  PConvIntToPyDictItem(dict,"formal_charge",at->formalCharge);
  PConvIntToPyDictItem(dict,"color",at->color);
  PConvIntToPyDictItem(dict,"cartoon",at->cartoon);
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
  PUnblock();
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
    "*PLockAPIAsGlut-DEBUG: entered as thread 0x%x\n",PyThread_get_thread_ident()
    ENDFD;

  PBlock();
  PRINTFD(FB_Threads)
    "#PLockAPIAsGlut-DEBUG: acquiring lock as thread 0x%x\n",PyThread_get_thread_ident()
    ENDFD;
  PXDecRef(PyObject_CallFunction(P_lock,NULL));
  while(P_glut_thread_keep_out) {
    /* IMPORTANT: keeps the glut thread out of an API operation... */
    /* NOTE: the keep_out variable can only be changed by the thread
       holding the API lock, therefore it is safe even through increment
       isn't atomic. */
    PRINTFD(FB_Threads)
      "-PLockAPIAsGlut-DEBUG: glut_thread_keep_out 0x%x\n",PyThread_get_thread_ident()
      ENDFD;
    
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
  PUnblock(); /* API is now locked, so we can free up Python...*/

  PRINTFD(FB_Threads)
    "=PLockAPIAsGlut-DEBUG: acquired\n"
    ENDFD;

}

/* THESE CALLS ARE REQUIRED FOR MONOLITHIC COMPILATION TO SUCCEED UNDER WINDOWS. */
#ifndef _PYMOL_ACTIVEX
#ifndef _EPYMOL
void	initExtensionClass(void);
void	initsglite(void);
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
#ifndef _EPYMOL
#ifdef WIN32
void	init_numpy();
void	initmultiarray();
void	initarrayfns();
void	initlapack_lite();
void	initumath();
void	initranlib();
#endif
#endif
#endif
#endif

#ifdef _PYMOL_MONOLITHIC
#ifndef _PYMOL_ACTIVEX
#ifndef _EPYMOL
void	initExtensionClass(void);
void	initsglite(void);
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

void PInitEmbedded(int argc,char **argv)
{
  /* This routine is called if we are running with an embedded Python interpreter */
  
  PyObject *args,*pymol;

#ifdef WIN32
  OrthoLineType path_buffer,command;
  HKEY phkResult;
  int lpcbData;
  int lpType = REG_SZ;
  int r1,r2;
#endif


#ifdef _PYMOL_SETUP_PY21 
  /* used by semistatic PyMOL */
{
  char line[5000];
  static char line1[5000];
  static char line2[5000];
  static char line3[5000];
  char *pymol_path;

  if(!getenv("PYMOL_PATH")) {
    if(getenv("PWD")) {
	strcpy(line1,"PYMOL_PATH=");
	strcat(line1,getenv("PWD"));
      putenv(line1);
	/* setenv("PYMOL_PATH",getenv("PWD"),1); */
    }
  }

  if(!getenv("PYTHONPATH")) { /* create PYTHONPATH */
    if(getenv("PYMOL_PATH")) {
	strcpy(line2,"PYTHONPATH=");
      strcat(line2,getenv("PYMOL_PATH"));
      strcat(line2,getenv("PYMOL_PATH"));
      strcat(line2,"/ext/lib/python2.1:");
      strcat(line2,"/ext/lib/python2.1/plat-linux2:");
      strcat(line2,getenv("PYMOL_PATH"));
      strcat(line2,"/ext/lib/python2.1/lib-tk:");
      strcat(line2,getenv("PYMOL_PATH"));
      strcat(line2,"/ext/lib/python2.1/lib-dynload");
putenv(line2);
      /*setenv("PYTHONPATH",line2,1);*/
    }
  } else { /* preempt existing PYTHONPATH */
    strcat(line3,"PYTHONPATH=");
      strcpy(line3,getenv("PYMOL_PATH"));
      strcat(line3,"/ext/lib/python2.1:");
      strcat(line3,getenv("PYMOL_PATH"));
      strcat(line3,"/ext/lib/python2.1/plat-linux2:");
      strcat(line3,getenv("PYMOL_PATH"));
      strcat(line3,"/ext/lib/python2.1/lib-tk:");
      strcat(line3,getenv("PYMOL_PATH"));
      strcat(line3,"/ext/lib/python2.1/lib-dynload:");
      strcat(line3,getenv("PYTHONPATH"));
      /*setenv("PYTHONPATH",line3,1);*/
putenv(line3);
  }
}
#endif

#ifdef _PYMOL_SETUP_PY22
  /* used by semistatic PyMOL */
{
static char line1[5000];
static char line2[5000];
static char line3[5000];
  char *pymol_path;

  if(!getenv("PYMOL_PATH")) {
    if(getenv("PWD")) {
	strcpy(line1,"PYMOL_PATH=");
	strcat(line1,getenv("PWD"));
      putenv(line1);
	/* setenv("PYMOL_PATH",getenv("PWD"),1); */
    }
  }

  if(!getenv("PYTHONPATH")) { /* create PYTHONPATH */
    if(getenv("PYMOL_PATH")) {
	strcpy(line2,"PYTHONPATH=");
      strcat(line2,getenv("PYMOL_PATH"));
      strcat(line2,"/ext/lib/python2.2:");
      strcat(line2,getenv("PYMOL_PATH"));
      strcat(line2,"/ext/lib/python2.2/plat-linux2:");
      strcat(line2,getenv("PYMOL_PATH"));
      strcat(line2,"/ext/lib/python2.2/lib-tk:");
      strcat(line2,getenv("PYMOL_PATH"));
      strcat(line2,"/ext/lib/python2.2/lib-dynload");
putenv(line2);
      /* setenv("PYTHONPATH",line2,1);*/
    }
  } else { /* preempt existing PYTHONPATH */
strcpy(line3,"PYTHONPATH=");
      strcat(line3,getenv("PYMOL_PATH"));
      strcat(line3,"/ext/lib/python2.2:");
      strcat(line3,getenv("PYMOL_PATH"));
      strcat(line3,"/ext/lib/python2.2/plat-linux2:");
      strcat(line3,getenv("PYMOL_PATH"));
      strcat(line3,"/ext/lib/python2.2/lib-tk:");
      strcat(line3,getenv("PYMOL_PATH"));
      strcat(line3,"/ext/lib/python2.2/lib-dynload");
      strcat(line3,getenv("PYTHONPATH"));
putenv(line3);
      /*setenv("PYTHONPATH",line3,1);*/
  }
}
#endif



#ifndef _PYMOL_ACTIVEX
#ifndef _EPYMOL
  Py_Initialize();
  PyEval_InitThreads();
#endif
#endif

  init_cmd();
#ifdef _PYMOL_MONOLITHIC
#ifndef _PYMOL_ACTIVEX
#ifndef _EPYMOL
	initExtensionClass();
	initsglite();
#ifdef WIN32
	/* initialize numeric python */
	init_numpy();
	initmultiarray();
	initarrayfns();
	initlapack_lite();
	initumath();
	initranlib();
#endif
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
			PyRun_SimpleString("if not os.environ.has_key('PYMOL_PATH'): os.environ['PYMOL_PATH']=_registry_pymol_path\n");
	  }
	RegCloseKey(phkResult);
	} 
  PyRun_SimpleString("if not os.environ.has_key('PYMOL_PATH'): os.environ['PYMOL_PATH']=os.getcwd()\n");
#endif

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


#ifdef _PYMOL_SETUP_PY21 
/* used by semistatic pymol */
  PyRun_SimpleString("import string");
  PyRun_SimpleString("sys.path=filter(lambda x:string.find(x,'warren/ext-static')<0,sys.path)"); /* clean bogus entries in sys.path */
#endif

#ifdef _PYMOL_SETUP_PY22
/* used by semistatic pymol */
  PyRun_SimpleString("import string");
  PyRun_SimpleString("sys.path=filter(lambda x:string.find(x,'warren/ext')<0,sys.path)"); /* clean bogus entries in sys.path */
#endif

#ifdef _PYMOL_SETUP_PY23
/* used by semistatic pymol */
  PyRun_SimpleString("import string");
  PyRun_SimpleString("sys.path=filter(lambda x:string.find(x,'warren/ext')<0,sys.path)"); /* clean bogus entries in sys.path */
#endif

#ifdef WIN32
  PyRun_SimpleString("if (os.environ['PYMOL_PATH']+'/modules') not in sys.path: sys.path.append(os.environ['PYMOL_PATH']+'/modules')\n");
#endif

  P_main = PyImport_AddModule("__main__");
  if(!P_main) ErrFatal("PyMOL","can't find '__main__'");

  /* inform PyMOL's other half that we're launching embedded-style */
  PyObject_SetAttrString(P_main,"pymol_launch",PyInt_FromLong(4));

  args = PConvStringListToPyList(argc,argv); /* prepare our argument list */
  if(!args) ErrFatal("PyMOL","can't process arguments.");

  /* copy arguments to __main__.pymol_argv */
  PyObject_SetAttrString(P_main,"pymol_argv",args);
  
  PyRun_SimpleString("if (os.environ['PYMOL_PATH']+'/modules') not in sys.path: sys.path.append(os.environ['PYMOL_PATH']+'/modules')\n"); /* needed for semistatic pymol */

  PyRun_SimpleString("import pymol"); /* create the global PyMOL namespace */

  pymol = PyImport_AddModule("pymol"); /* get it */
  if(!pymol) ErrFatal("PyMOL","can't find module 'pymol'");

}

void PGetOptions(PyMOLOptionRec *rec)
{
  PyObject *pymol,*invocation,*options;
  char *load_str;

  pymol = PyImport_AddModule("pymol"); /* get it */
  if(!pymol) ErrFatal("PyMOL","can't find module 'pymol'");

  invocation = PyObject_GetAttrString(pymol,"invocation"); /* get a handle to the invocation module */
  if(!pymol) ErrFatal("PyMOL","can't find module 'invocation'");

  options = PyObject_GetAttrString(invocation,"options");
  if(!pymol) ErrFatal("PyMOL","can't get 'invocation.options'.");

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
  
  load_str = PyString_AsString(PyObject_GetAttrString(options,"after_load_script"));
  if(load_str) {
    if(load_str[0]) {
      UtilNCopy(rec->after_load_script,load_str,PYMOL_MAX_OPT_STR);
    }
  }
  if(PyErr_Occurred()) {
    PyErr_Print();
  }
}

void PRunString(char *str) /* runs a string in the global PyMOL module namespace */
{
  PXDecRef(PyObject_CallFunction(P_exec,"s",str));
}

void PInit(void) 
{
  PyObject *pymol,*sys,*pcatch;
  int a;

#ifdef PYMOL_NEW_THREADS
   PyEval_InitThreads();
#endif

#ifdef WIN32
#ifdef _PYMOL_MONOLITHIC
#ifndef _PYMOL_ACTIVEX
#ifndef _EPYMOL
	/* Win32 module build: includes pyopengl, numpy, and sglite */
	/* sglite */
	initExtensionClass();
	initsglite();
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
#endif
#endif
#endif


  for(a=0;a<MAX_SAVED_THREAD;a++) {
    SavedThread[a].id=-1;
  }

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

  PRunString("import traceback\n");  
  P_traceback = PyDict_GetItemString(P_globals,"traceback");
  if(!P_traceback) ErrFatal("PyMOL","can't find 'traceback'");

  PRunString("import cmd\n");  
  P_cmd = PyDict_GetItemString(P_globals,"cmd");
  if(!P_cmd) ErrFatal("PyMOL","can't find 'cmd'");

  P_lock = PyObject_GetAttrString(P_cmd,"lock");
  if(!P_lock) ErrFatal("PyMOL","can't find 'pm.lock()'");

  P_unlock = PyObject_GetAttrString(P_cmd,"unlock");
  if(!P_unlock) ErrFatal("PyMOL","can't find 'pm.unlock()'");

  P_lock_c = PyObject_GetAttrString(P_cmd,"lock_c");
  if(!P_lock_c) ErrFatal("PyMOL","can't find 'pm.lock_c()'");

  P_unlock_c = PyObject_GetAttrString(P_cmd,"unlock_c");
  if(!P_unlock_c) ErrFatal("PyMOL","can't find 'pm.unlock_c()'");

  PRunString("import menu\n");  
  P_menu = PyDict_GetItemString(P_globals,"menu");
  if(!P_menu) ErrFatal("PyMOL","can't find module 'menu'");

  PRunString("import setting\n");  
  P_setting = PyDict_GetItemString(P_globals,"setting");
  if(!P_setting) ErrFatal("PyMOL","can't find module 'setting'");

  PRunString("import povray\n");  
  P_povray = PyDict_GetItemString(P_globals,"povray");
  if(!P_povray) ErrFatal("PyMOL","can't find module 'povray'");

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

  PRunString("from chempy.bonds import bonds"); /* load bond dictionary */

  PRunString("from chempy import models"); 
  P_models = PyDict_GetItemString(P_globals,"models");
  if(!P_models) ErrFatal("PyMOL","can't find 'chempy.models'");

  PRunString("import util\n");  
  PRunString("import contrib\n");
  /*#ifdef _PYMOL_XRAY
  PRunString("import sglite\n"); 
  #endif*/

  PRunString("import string\n"); 

  /* backwards compatibility */

  PRunString("pm = cmd\n");  
  PRunString("pmu = util\n");  

  PRunString("glutThread = thread.get_ident()");

  P_glut_thread_id = PyThread_get_thread_ident();

  #ifndef WIN32
  if(PyMOLOption->siginthand) {
    signal(SIGINT,my_interrupt);
  }
  #endif

  /* required environment variables */

  PyRun_SimpleString(
"if not os.environ.has_key('PYMOL_DATA'): os.environ['PYMOL_DATA']=os.environ['PYMOL_PATH']+'/data'");
  PyRun_SimpleString(
"os.environ['TUT']=os.environ['PYMOL_DATA']+'/tut'");

  PyRun_SimpleString(
"if not os.environ.has_key('PYMOL_SCRIPTS'): os.environ['PYMOL_SCRIPTS']=os.environ['PYMOL_PATH']+'/scripts'");

}

int PPovrayRender(char *header,char *inp,char *file,int width,int height,int antialias) 
{
  PyObject *result;
  int ok;
  PBlock();
  result = PyObject_CallMethod(P_povray,"render_from_string","sssiii",header,inp,file,width,height,antialias);
  ok = PyObject_IsTrue(result);
  Py_DECREF(result);
  PUnblock();
  return(ok);
}

void PSGIStereo(int flag) 
{
  int blocked;
  blocked = PAutoBlock();
  if(flag) 
    PRunString("cmd._sgi_stereo(1)");
  else
    PRunString("cmd._sgi_stereo(0)");
  if(blocked) PUnblock();
}

void PFree(void)
{
}

void PExit(int code)
{
  ExecutiveDelete("all");
  PBlock();
  MainFree();
  Py_Exit(code);
}

void PParse(char *str) 
{
  OrthoCommandIn(str);
}

void PLog(char *str,int format) 
     /* general log routine can write PML 
        or PYM commands to appropriate log file */
{  
  int mode;
  int a;
  int blocked;
  PyObject *log;
  OrthoLineType buffer="";
  mode = (int)SettingGet(cSetting_logging);
  if(mode)
    {
      blocked = PAutoBlock();
      log = PyDict_GetItemString(P_globals,P_log_file_str);
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
      PAutoUnblock(blocked);
    }
}

void PLogFlush(void)
{
  int mode;
  PyObject *log;
  int blocked;
  mode = (int)SettingGet(cSetting_logging);
  if(mode)
    {
      blocked = PAutoBlock();
      log = PyDict_GetItemString(P_globals,P_log_file_str);
      if(log&&(log!=Py_None)) {
        PyObject_CallMethod(log,"flush","");
      }
      PAutoUnblock(blocked);
    }
}

static void PDeleteAll(void *p)
{ /* assumes blocked and unlocked API */
  if(!PyMOLTerminating)
    PyObject_CallMethod(P_cmd,"delete","s","all");
}

static void PMaintainObjectAll(void) 
/* legacy exception for "del all", which is broken by version 0.86
   "del all" now means what it should in Python: delete an object
   referenced by "all".  In order to support old scripts, we create an
   "all" object and then catch its destruction event...(yes, I know
   this is weak, but its is the most compatible solution...).  */
{
  PyObject *all;

  all = PyDict_GetItemString(P_globals,"all");
  if(all==NULL) {
    all = PyCObject_FromVoidPtr(NULL,PDeleteAll);
    PyDict_SetItemString(P_globals,"all",all);
    Py_DECREF(all);
  }
}


void PFlush(void) {  
  /* NOTE: ASSUMES unblocked Python threads and a locked API */
  PyObject *err;
  char buffer[OrthoLineLength+1];
  while(OrthoCommandOut(buffer)) {
    PBlockAndUnlockAPI();
    PMaintainObjectAll();
    PXDecRef(PyObject_CallFunction(P_parse,"s",buffer));
    err = PyErr_Occurred();
    if(err) {
      PyErr_Print();
      PRINTFB(FB_Python,FB_Errors)
        " PFlush: Uncaught exception.  PyMOL may have a bug.\n"
        ENDFB;
    }
    PLockAPIAndUnblock();
  }
}

void PFlushFast(void) {
  /* NOTE: ASSUMES we currently have blocked Python threads and an unlocked API */ 
  PyObject *err;
  char buffer[OrthoLineLength+1];
  while(OrthoCommandOut(buffer)) {
    PMaintainObjectAll();
    PRINTFD(FB_Threads)
      " PFlushFast-DEBUG: executing '%s' as thread 0x%x\n",buffer,
      PyThread_get_thread_ident()
      ENDFD;
    PXDecRef(PyObject_CallFunction(P_parse,"s",buffer));
    err = PyErr_Occurred();
    if(err) {
      PyErr_Print();
      PRINTFB(FB_Python,FB_Errors)
        " PFlushFast: Uncaught exception.  PyMOL may have a bug.\n"
        ENDFB;
    }
  }
}


void PBlock(void)
{

  if(!PAutoBlock()) {
    ErrFatal("PBlock","Threading error detected.  Terminating...");
  }
}


int PAutoBlock(void)
{
#ifndef _PYMOL_ACTIVEX
#ifndef _EPYMOL
  int a,id;
  /* synchronize python */

  id = PyThread_get_thread_ident();
  PRINTFD(FB_Threads)
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
      
      PRINTFD(FB_Threads)
        " PAutoBlock-DEBUG: seeking global lock 0x%x\n",id
      ENDFD;

#ifdef PYMOL_NEW_THREADS

      PyEval_AcquireLock();

      PRINTFD(FB_Threads)
        " PAutoBlock-DEBUG: restoring 0x%x\n",id
      ENDFD;
      
      PyThreadState_Swap((SavedThread+a)->state);

#else
      PRINTFD(FB_Threads)
        " PAutoBlock-DEBUG: restoring 0x%x\n",id
      ENDFD;
      
      PyEval_RestoreThread((SavedThread+a)->state);
#endif
      
      PRINTFD(FB_Threads)
        " PAutoBlock-DEBUG: restored 0x%x\n",id
      ENDFD;

      PRINTFD(FB_Threads)
        " PAutoBlock-DEBUG: clearing 0x%x\n",id
      ENDFD;

      PXDecRef(PyObject_CallFunction(P_lock_c,NULL));
      SavedThread[a].id = -1; 
      /* this is the only safe time we can change things */
      PXDecRef(PyObject_CallFunction(P_unlock_c,NULL));
      
      PRINTFD(FB_Threads)
        " PAutoBlock-DEBUG: blocked 0x%x (0x%x, 0x%x, 0x%x)\n",PyThread_get_thread_ident(),
        SavedThread[MAX_SAVED_THREAD-1].id,
        SavedThread[MAX_SAVED_THREAD-2].id,
        SavedThread[MAX_SAVED_THREAD-3].id
        ENDFD;

      return 1;
    }
    a--;
  }
  PRINTFD(FB_Threads)
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

void PUnblock(void)
{
#ifndef _PYMOL_ACTIVEX
#ifndef _EPYMOL
  int a;
  /* NOTE: ASSUMES a locked API */

  PRINTFD(FB_Threads)
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
  PRINTFD(FB_Threads)
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
    if(Feedback(FB_Python,FB_Output)) {
      OrthoAddOutput(str);
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


