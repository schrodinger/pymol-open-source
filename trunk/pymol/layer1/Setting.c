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
#include"os_std.h"

#include"Base.h"
#include"OOMac.h"
#include"MemoryDebug.h"
#include"Ortho.h"
#include"Setting.h"
#include"Scene.h"
#include"ButMode.h"
#include"Executive.h"
#include"Editor.h"
#include"P.h"
#include"Util.h"
#include"main.h"
#include"PConv.h"
#include"Wizard.h"
#include"Seq.h"
#include"PyMOLOptions.h"

static void *SettingPtr(CSetting *I,int index,unsigned int size);

int SettingSetSmart_i(PyMOLGlobals *G,CSetting *set1,CSetting *set2,int index, int value)
{
  int dummy;
  if(set1 && SettingGetIfDefined_i(G,set1,index,&dummy)) {
    return SettingSet_i(set1,index,value);
  }
  if(set2 && SettingGetIfDefined_i(G,set2,index,&dummy)) {
    return SettingSet_i(set2,index,value);
  }
  return SettingSetGlobal_i(G,index,value);
}

int SettingSetGlobalsFromPyList(PyMOLGlobals *G,PyObject *list)
{
#ifdef _PYMOL_NOPY
  return 0;
#else
  int ok=true;
  int stereo=SettingGetGlobal_b(G,cSetting_stereo);
  int session_migration=SettingGetGlobal_b(G,cSetting_session_migration);
  int session_version_check=SettingGetGlobal_b(G,cSetting_session_version_check);

  register CSetting *I=G->Setting;
  if(list)
    if(PyList_Check(list)) 
      ok = SettingFromPyList(I,list);
  
  SettingSet_i(I,cSetting_security,G->Security); /* always override Security setting with global variable */
  SettingSet_b(I,cSetting_stereo,stereo); /* preserve current stereo mode */
  SettingSet_b(I,cSetting_session_migration,session_migration); /* preserve current migration info */
  SettingSet_b(I,cSetting_session_version_check,session_version_check);
  return(ok);
#endif
}

PyObject *SettingGetGlobalsPyList(PyMOLGlobals *G)
{
#ifdef _PYMOL_NOPY
  return NULL;
#else
  PyObject *result = NULL;
  register CSetting *I=G->Setting;
  result = SettingAsPyList(I);
  return(PConvAutoNone(result));
#endif
}


#ifndef _PYMOL_NOPY
static PyObject *get_list(CSetting *I,int index)
{
  PyObject *result = NULL;
  int setting_type = I->info[index].type;
  switch(setting_type) {

  case cSetting_boolean:
  case cSetting_int:
  case cSetting_color:
    result = PyList_New(3);
    PyList_SetItem(result,0,PyInt_FromLong(index));
    PyList_SetItem(result,1,PyInt_FromLong(setting_type));
    PyList_SetItem(result,2,PyInt_FromLong(*((int*)(I->data+I->info[index].offset))));
    break;
  case cSetting_float:
    result = PyList_New(3);
    PyList_SetItem(result,0,PyInt_FromLong(index));
    PyList_SetItem(result,1,PyInt_FromLong(setting_type));
    PyList_SetItem(result,2,PyFloat_FromDouble(*((float*)(I->data+I->info[index].offset))));
    break;
  case cSetting_float3:
    result = PyList_New(3);
    PyList_SetItem(result,0,PyInt_FromLong(index));
    PyList_SetItem(result,1,PyInt_FromLong(setting_type));
    PyList_SetItem(result,2,PConvFloatArrayToPyList(((float*)(I->data+I->info[index].offset)),3));
    break;
  case cSetting_string:
    result = PyList_New(3);
    PyList_SetItem(result,0,PyInt_FromLong(index));
    PyList_SetItem(result,1,PyInt_FromLong(setting_type));
    PyList_SetItem(result,2,PyString_FromString(((char*)(I->data+I->info[index].offset))));
    break;
  }
  return(PConvAutoNone(result));
}
#endif

PyObject *SettingAsPyList(CSetting *I)
{
#ifdef _PYMOL_NOPY
  return NULL;
#else

  PyObject *result = NULL;
  int cnt = 0;
  int a;

  if(I) {
    for(a=0;a<cSetting_INIT;a++) {
      if(I->info[a].defined)
        cnt++;
    }
    result = PyList_New(cnt);
    cnt=0;
    for(a=0;a<cSetting_INIT;a++) {
      if(I->info[a].defined) {
        PyList_SetItem(result,cnt,get_list(I,a));
        cnt++;
      }
    }
  }
  return(PConvAutoNone(result));
#endif
}
/*========================================================================*/
#ifndef _PYMOL_NOPY
static int set_list(CSetting *I,PyObject *list)
{
  int ok=true;
  int index;
  int setting_type;
  char *str;
  if(list!=Py_None) {
    if(ok) ok=(list!=NULL);
    if(ok) ok=PyList_Check(list);
    if(ok) ok=PConvPyIntToInt(PyList_GetItem(list,0),&index);
    if(ok) ok=PConvPyIntToInt(PyList_GetItem(list,1),&setting_type);
    if(ok&&(index<cSetting_INIT)) { /* ignore unknown settings */
      switch(index) { 
        /* don't restore the folllowing settings,
           which are inherently system-dependent */
      case cSetting_stereo_double_pump_mono: 
      case cSetting_max_threads:
      case cSetting_session_migration:
        break;
      default:
        if(ok) switch(setting_type) {
        case cSetting_boolean:
        case cSetting_int:
        case cSetting_color:
          ok = PConvPyIntToInt(PyList_GetItem(list,2),
                               (int*)SettingPtr(I,index,sizeof(int)));
          break;
        case cSetting_float:
          ok = PConvPyFloatToFloat(PyList_GetItem(list,2),
                                   (float*)SettingPtr(I,index,sizeof(float)));
          break;
        case cSetting_float3:
          ok = PConvPyListToFloatArrayInPlaceAutoZero(PyList_GetItem(list,2),
                                                      (float*)SettingPtr(I,index,3*sizeof(float)),3);
          break;
        case cSetting_string:
          ok = PConvPyStrToStrPtr(PyList_GetItem(list,2),&str);
          if(ok) {
            strcpy(((char*)SettingPtr(I,index,strlen(str)+1)),str);
          }
          break;
        }
      }
      if(ok) I->info[index].type=setting_type;
    }
  }
  return(ok);
}
#endif
/*========================================================================*/
CSetting *SettingNewFromPyList(PyMOLGlobals *G,PyObject *list)
{
#ifdef _PYMOL_NOPY
  return NULL;
#else

  int ok=true;
  int size;
  int a;
  CSetting *I=NULL;
  if(ok) ok=(list!=NULL);
  if(ok) ok=PyList_Check(list);
  if(ok) {
    I=SettingNew(G);
    size=PyList_Size(list);
    for(a=0;a<size;a++) {
      if(ok) ok=set_list(I,PyList_GetItem(list,a));
    }
  }
  return(I);
#endif
}
/*========================================================================*/
int SettingFromPyList(CSetting *I,PyObject *list)
{
#ifdef _PYMOL_NOPY
  return 0;
#else

  int ok=true;
  int size;
  int a;
  if(ok) ok=(I!=NULL);
  if(ok) ok=PyList_Check(list);
  if(ok) {
    size=PyList_Size(list);
    for(a=0;a<size;a++) {
      if(ok) ok=set_list(I,PyList_GetItem(list,a));
    }
  }
  return(ok);
#endif
}


/*========================================================================*/
#ifndef _PYMOL_NOPY
PyObject *SettingGetUpdateList(PyMOLGlobals *G,CSetting *I)
{ /* assumes blocked interpreter */

  int a;
  int n;
  PyObject *result;

  if(!I) I=G->Setting; /* fall back on global settings */

  n=VLAGetSize(I->info);
  result=PyList_New(0);
  for(a=0;a<n;a++) {
    if(I->info[a].changed) {
      I->info[a].changed=false;
      PyList_Append(result,PyInt_FromLong(a));
    }
  }
  return(result);

}
#endif

/*========================================================================*/
void SettingCheckHandle(PyMOLGlobals *G,CSetting **handle)
{
  if(!*handle)
    *handle=SettingNew(G);
}
/*========================================================================*/
int SettingGetTextValue(PyMOLGlobals *G,CSetting *set1,CSetting *set2,int index,char *buffer) 
/* not range checked */
{
  int type;
  int ok=true;
  int tmp1;
  float *ptr;
  type = SettingGetType(G,index);
  switch(type) {
  case cSetting_boolean:
    if(SettingGet_b(G,set1,set2,index))
      sprintf(buffer,"on");
    else
      sprintf(buffer,"off");      
    break;
  case cSetting_int:
    sprintf(buffer,"%d",SettingGet_i(G,set1,set2,index));
    break;
  case cSetting_float:
    sprintf(buffer,"%1.5f",SettingGet_f(G,set1,set2,index));
    break;
  case cSetting_float3:
    ptr = SettingGet_3fv(G,set1,set2,index);
    sprintf(buffer,"[ %1.5f, %1.5f, %1.5f ]",
            ptr[0],ptr[1],ptr[2]);
    break;
  case cSetting_color:
    tmp1 = SettingGet_color(G,set1,set2,index);
    if(tmp1<0) 
      strcpy(buffer,"default");
    else
      strcpy(buffer,ColorGetName(G,tmp1));
    break;
  case cSetting_string:
    strcpy(buffer,SettingGet_s(G,set1,set2,index));
    break;
  default:
    ok=false;
    break;
  }
  return(ok);
}

/*========================================================================*/
#ifndef _PYMOL_NOPY
int SettingSetTuple(PyMOLGlobals *G,CSetting *I,int index,PyObject *tuple) 
     /* must have interpret locked to make this call */
{
  PyObject *value;
  int type;
  int ok=true;
  if(!I) I=G->Setting; /* fall back on global settings */

  /* this data structure has been pre-checked at the python level... */

  type  = PyInt_AsLong(PyTuple_GetItem(tuple,0));
  value = PyTuple_GetItem(tuple,1);
  switch(type) {
  case cSetting_boolean:
    SettingSet_b(I,index,PyInt_AsLong(PyTuple_GetItem(value,0)));
    break;
  case cSetting_int:
    SettingSet_i(I,index,PyInt_AsLong(PyTuple_GetItem(value,0)));
    break;
  case cSetting_float:
   SettingSet_f(I,index,(float)PyFloat_AsDouble(PyTuple_GetItem(value,0)));
    break;
  case cSetting_float3:
        SettingSet_3f(I,index,
                  (float)PyFloat_AsDouble(PyTuple_GetItem(value,0)),
                  (float)PyFloat_AsDouble(PyTuple_GetItem(value,1)),
                  (float)PyFloat_AsDouble(PyTuple_GetItem(value,2)));
    break;
  case cSetting_color:
    SettingSet_color(I,index,
                     PyString_AsString(PyTuple_GetItem(value,0)));
    break;
  case cSetting_string:
    SettingSet_s(I,index,
                 PyString_AsString(PyTuple_GetItem(value,0)));
    break;
  default:
    ok=false;
    break;
  }
  return(ok);
}
#endif
/*========================================================================*/
#ifndef _PYMOL_NOPY
PyObject *SettingGetTuple(PyMOLGlobals *G,CSetting *set1,CSetting *set2,int index)
{  /* assumes blocked python interpeter */
  PyObject *result = NULL;
  float *ptr;
  int type = SettingGetType(G,index);

  switch(type) {
  case cSetting_boolean:
    result = Py_BuildValue("(i(i))",type,
                           SettingGet_b(G,set1,set2,index));
    break;
  case cSetting_int:
    result = Py_BuildValue("(i(i))",type,
                           SettingGet_i(G,set1,set2,index));
    break;
  case cSetting_float:
    result = Py_BuildValue("(i(f))",type,
                           SettingGet_f(G,set1,set2,index));
    break;
  case cSetting_float3:
    ptr =  SettingGet_3fv(G,set1,set2,index);
    result = Py_BuildValue("(i(fff))",type,
                           ptr[0],ptr[1],ptr[2]);
    break;
  case cSetting_color:
    result = Py_BuildValue("(i(i))",type,
                           SettingGet_color(G,set1,set2,index));
    break;
  case cSetting_string:
    result = Py_BuildValue("(i(s))",type,
                           SettingGet_s(G,set1,set2,index));
    break;
  default:
    Py_INCREF(Py_None);
    result = Py_None;
    break;
  }
  return result;
}
#endif
/*========================================================================*/
#ifndef _PYMOL_NOPY
PyObject *SettingGetDefinedTuple(PyMOLGlobals *G,CSetting *set1,int index)
{  /* Assumes blocked Python interpreter */
  PyObject *result = NULL;
  int defined = true;
  int type = SettingGetType(G,index);
  int int1;
  float float1,*vect1;
  char *str1;
  switch(type) {
  case cSetting_boolean:
    defined = SettingGetIfDefined_b(G,set1,index,&int1);
    if(defined) 
      result = Py_BuildValue("(i(i))",type,int1);
    break;
  case cSetting_int:
    defined = SettingGetIfDefined_i(G,set1,index,&int1);
    if(defined) 
      result = Py_BuildValue("(i(i))",type,int1);
    break;
  case cSetting_float:
    defined = SettingGetIfDefined_f(G,set1,index,&float1);
    if(defined) 
      result = Py_BuildValue("(i(f))",type,float1);
    break;
  case cSetting_float3:
    defined = SettingGetIfDefined_3fv(G,set1,index,&vect1);
    result = Py_BuildValue("(i(fff))",type,
                           vect1[0],vect1[1],vect1[2]);
    break;
  case cSetting_color:
    defined = SettingGetIfDefined_color(G,set1,index,&int1);
    if(defined) 
      result = Py_BuildValue("(i(i))",type,int1);
    break;
  case cSetting_string:
    defined = SettingGetIfDefined_s(G,set1,index,&str1);
    if(defined) 
      result = Py_BuildValue("(i(s))",type,str1);
    break;
  default:
    break;
  }
  if(!defined)  {
    result = Py_BuildValue("(i)",0);
    }
  if(!result) {
    Py_INCREF(Py_None);
    result = Py_None;
  }
  return result;
}
#endif

/*========================================================================*/
CSetting *SettingNew(PyMOLGlobals *G)
{
  OOAlloc(G,CSetting);
  SettingInit(G,I);
  return(I);
}
/*========================================================================*/
void SettingPurge(CSetting *I)
{
  if(I) {
    VLAFreeP(I->data);
    VLAFreeP(I->info);
    I->size=0;
  }
}
/*========================================================================*/
void SettingFreeP(CSetting *I)
{
  if(I) SettingPurge(I);
  OOFreeP(I);
}
/*========================================================================*/
void SettingInit(PyMOLGlobals *G,CSetting *I)
{
  I->G=G;
  I->size=sizeof(int); /* insures offset is never zero, except when undef */
  I->data=VLAlloc(char,10);
  I->info=VLAMalloc(cSetting_INIT,sizeof(SettingRec),5,1); /* auto-zero */
}
/*========================================================================*/
void SettingClear(CSetting *I,int index)
{
  if(I)
    I->info[index].defined = false; 
}
/*========================================================================*/
static void *SettingPtr(CSetting *I,int index,unsigned int size)
{
  SettingRec *sr = I->info+index;
  if(size<sizeof(int)) size=sizeof(int); /* make sure we're word aligned */
  while(size&(sizeof(int)-1)) size++;
                         
  if((!sr->offset)||(sr->max_size<size)) { 
    sr->offset=I->size;
    I->size+=size;
    sr->max_size=size;
    VLACheck(I->data,char,I->size);
  }
  sr->defined = true;
  sr->changed = true;
  return(I->data+sr->offset);
}
/*========================================================================*/
int SettingUnset(CSetting *I,int index)
{
  if(I) {
    SettingRec *sr = I->info+index;
    sr->defined = false;
    sr->changed = true;
  }
  return true;
}
/*========================================================================*/
int SettingGetType(PyMOLGlobals *G,int index)
{
  register CSetting *I=G->Setting;
  return(I->info[index].type);
}

/*========================================================================*/
static int get_i(CSetting *I,int index)
{
  PyMOLGlobals *G=I->G;
  int result;
  switch(I->info[index].type) {
  case cSetting_boolean:
  case cSetting_int:
  case cSetting_color:
    result = (*((int*)(I->data+I->info[index].offset)));
    break;
  case cSetting_float:
    result = (int)(*((float*)(I->data+I->info[index].offset)));
    break;
  default:
    PRINTFB(G,FB_Setting,FB_Errors)
      "Setting-Error: type read mismatch (int) %d\n",index
      ENDFB(G);
    result = 0;
    break;
  }
  return(result);
}
/*========================================================================*/
static int get_b(CSetting *I,int index)
{
  int result;
  PyMOLGlobals *G=I->G;
  switch(I->info[index].type) {
  case cSetting_boolean:
  case cSetting_int:
  case cSetting_color:
    result = (*((int*)(I->data+I->info[index].offset)));
    break;
  case cSetting_float:
    result = (int)(*((float*)(I->data+I->info[index].offset)));
    break;
  default:
    PRINTFB(G,FB_Setting,FB_Errors)
      "Setting-Error: type read mismatch (boolean) %d\n",index
      ENDFB(G);
    result = 0;
  }
  return(result);
}
/*========================================================================*/
static int get_color(CSetting *I,int index)
{
  int result;
  PyMOLGlobals *G=I->G;
  switch(I->info[index].type) {
  case cSetting_boolean:
  case cSetting_int:
  case cSetting_color:
    result = (*((int*)(I->data+I->info[index].offset)));
    break;
  case cSetting_float:
    result = (int)(*((float*)(I->data+I->info[index].offset)));
    break;
  default:
    PRINTFB(G,FB_Setting,FB_Errors)
      "Setting-Error: type read mismatch (color) %d\n",index
      ENDFB(G);
    result = 0;
  }
  return(result);
}
/*========================================================================*/
static float get_f(CSetting *I,int index)
{
  float result;
  PyMOLGlobals *G=I->G;
  switch(I->info[index].type) {
  case cSetting_boolean:
  case cSetting_int:
  case cSetting_color:
    result = (float)(*((int*)(I->data+I->info[index].offset)));
    break;
  case cSetting_float:
    result = (*((float*)(I->data+I->info[index].offset)));
    break;
  default:
    PRINTFB(G,FB_Setting,FB_Errors)
      "Setting-Error: type read mismatch (float) %d\n",index
      ENDFB(G);
    result = 0.0F;
  }
  return(result);
}
/*========================================================================*/
static char *get_s(CSetting *I,int index)
{
  char *result;
  PyMOLGlobals *G=I->G;
  switch(I->info[index].type) {
  case cSetting_string:
    result = ((char*)(I->data+I->info[index].offset));
    break;
  default:
    PRINTFB(G,FB_Setting,FB_Errors)
      "Setting-Error: type read mismatch (string) %d\n",index
      ENDFB(G);
    result = NULL;
  }
  return(result);
}
/*========================================================================*/
int SettingSet_b(CSetting *I,int index, int value)
{
  int ok=true;
  if(I) {
    int setting_type = I->info[index].type;
    PyMOLGlobals *G=I->G;
    switch(setting_type) {
    case cSetting_blank:
    case cSetting_boolean:
    case cSetting_int:
    case cSetting_color:
      VLACheck(I->info,SettingRec,index);
      *((int*)SettingPtr(I,index,sizeof(int))) = value;
      break;
    case cSetting_float:
      *((float*)SettingPtr(I,index,sizeof(float))) = (float)value;
      break;
    default:
      PRINTFB(G,FB_Setting,FB_Errors)
        "Setting-Error: type set mismatch (boolean) %d\n",index
        ENDFB(G);
        ok=false;
    }
    if(setting_type==cSetting_blank)
      I->info[index].type = cSetting_boolean;
  } else {
    ok=false;
  }
  return(ok);
}
/*========================================================================*/
int SettingSet_i(CSetting *I,int index, int value)
{
  int ok=true;
  if(I) {
    PyMOLGlobals *G=I->G;
    int setting_type = I->info[index].type;
    switch(setting_type) {
    case cSetting_blank:
    case cSetting_boolean:
    case cSetting_int:
    case cSetting_color:
      VLACheck(I->info,SettingRec,index);
      *((int*)SettingPtr(I,index,sizeof(int))) = value;
      break;
    case cSetting_float:
      *((float*)SettingPtr(I,index,sizeof(float))) = (float)value;
      break;
    default:
      PRINTFB(G,FB_Setting,FB_Errors)
        "Setting-Error: type set mismatch (integer)\n"
        ENDFB(G);
        ok=false;
    }
    if(setting_type==cSetting_blank)
      I->info[index].type = cSetting_int;
  } else {
    ok=false;
  }
  return(ok);
}
/*========================================================================*/
int SettingSet_color(CSetting *I,int index, char *value)
{
  int ok=true;
  int color_index;
  if(I) {
    PyMOLGlobals *G=I->G;
    color_index=ColorGetIndex(G,value);
    if((color_index==-1)&&(strcmp(value,"-1"))) {
      PRINTFB(G,FB_Setting,FB_Errors)
        "Setting-Error: unknown color '%s'\n",value
        ENDFB(G);
      ok=false;
    
    } else {
    
      int setting_type = I->info[index].type;
      switch(setting_type) {
      case cSetting_blank:
      case cSetting_boolean:
      case cSetting_int:
      case cSetting_color:
        VLACheck(I->info,SettingRec,index);
        *((int*)SettingPtr(I,index,sizeof(int))) = color_index;
        break;
      case cSetting_float:
        *((float*)SettingPtr(I,index,sizeof(float))) = (float)color_index;
        break;
      default:
        PRINTFB(G,FB_Setting,FB_Errors)
          "Setting-Error: type set mismatch (color)\n"
          ENDFB(G);
          ok=false;
      }
      if(setting_type==cSetting_blank)
        I->info[index].type = cSetting_color;
    
      VLACheck(I->info,SettingRec,index);
      *((int*)SettingPtr(I,index,sizeof(int))) = color_index;
      I->info[index].type = cSetting_color;

    }
  }
  return(ok);
}
/*========================================================================*/
int SettingSet_f(CSetting *I,int index, float value)
{
  int ok=true;
  if(I) {
    PyMOLGlobals *G=I->G;
    int setting_type = I->info[index].type;
    switch(setting_type) {
    case cSetting_boolean:
    case cSetting_int:
    case cSetting_color:
      VLACheck(I->info,SettingRec,index);
      *((int*)SettingPtr(I,index,sizeof(int))) = (int)value;
      break;
    case cSetting_blank:
    case cSetting_float:
      *((float*)SettingPtr(I,index,sizeof(float))) = value;
      break;
    default:
      PRINTFB(G,FB_Setting,FB_Errors)
        "Setting-Error: type set mismatch (float)\n"
        ENDFB(G);
        ok=false;
    }
    if(setting_type==cSetting_blank)
      I->info[index].type = cSetting_float;
  } else {
    ok=false;
  }
  return(ok);
}
/*========================================================================*/
int SettingSet_s(CSetting *I,int index, char *value)
{
  int ok=true;
  if(I) {
    PyMOLGlobals *G=I->G;
    int setting_type = I->info[index].type;
    switch(setting_type) {
    case cSetting_blank:
    case cSetting_string:
      VLACheck(I->info,SettingRec,index);
      strcpy(((char*)SettingPtr(I,index,strlen(value)+1)),value);
      I->info[index].type = cSetting_string;
      break;
    default:
      PRINTFB(G,FB_Setting,FB_Errors)
        "Setting-Error: type set mismatch (string)\n"
        ENDFB(G);
        ok=false;
    }
    if(setting_type==cSetting_blank)
      I->info[index].type = cSetting_string;
  } else {
    ok=false;
  }
  return(ok);
}
/*========================================================================*/
int SettingSet_3f(CSetting *I,int index, float value1,float value2,float value3)
{
  int ok=false;
  float *ptr;
  if(I) {
    PyMOLGlobals *G=I->G;
    int setting_type = I->info[index].type;
    switch(setting_type) {
    case cSetting_blank:
    case cSetting_float3:
      VLACheck(I->info,SettingRec,index);
      ptr = (float*)SettingPtr(I,index,sizeof(float)*3);
      ptr[0]=value1;
      ptr[1]=value2;
      ptr[2]=value3;
      break;
    default:
      PRINTFB(G,FB_Setting,FB_Errors)
        "Setting-Error: type set mismatch (float3)\n"
        ENDFB(G);
        ok=false;
    }
    if(setting_type==cSetting_blank)
      I->info[index].type = cSetting_float3;
  } else {
    ok=false;
  }
  return(ok);
}
/*========================================================================*/
int SettingSet_3fv(CSetting *I,int index, float *vector)
{
  float *ptr;
  VLACheck(I->info,SettingRec,index);
  ptr = (float*)SettingPtr(I,index,sizeof(float)*3);
  copy3f(vector,ptr);
  I->info[index].type = cSetting_float3;
  return(true);
}
/*========================================================================*/
int   SettingGetGlobal_b(PyMOLGlobals *G,int index) 
{
  register CSetting *I=G->Setting;
  return(get_b(I,index));
}
/*========================================================================*/
int   SettingGetGlobal_i(PyMOLGlobals *G,int index) 
{
  register CSetting *I=G->Setting;
  return(get_i(I,index));
}
/*========================================================================*/
float SettingGetGlobal_f(PyMOLGlobals *G,int index)
{
  register CSetting *I=G->Setting;
  return(get_f(I,index));
}
/*========================================================================*/
char *SettingGetGlobal_s(PyMOLGlobals *G,int index)
{
  register CSetting *I=G->Setting;
  return(get_s(I,index));
}
/*========================================================================*/
void  SettingGetGlobal_3f(PyMOLGlobals *G,int index,float *value)
{
  register CSetting *I=G->Setting;
  float *ptr;
  ptr = (float*)(I->data+I->info[index].offset);
  copy3f(ptr,value);
}
/*========================================================================*/
float *SettingGetGlobal_3fv(PyMOLGlobals *G,int index)
{
  register CSetting *I=G->Setting;
  return (float*)(I->data+I->info[index].offset);
}
/*========================================================================*/
int   SettingGet_b(PyMOLGlobals *G,CSetting *set1,CSetting *set2,int index)
{
  if(set1) {
    if(set1->info[index].defined) {
      return(get_b(set1,index));
    }
  }
  if(set2) { 
    if(set2->info[index].defined) {
      return(get_b(set2,index));
    }
  }
  return(SettingGetGlobal_i(G,index));
}
/*========================================================================*/
int SettingGetIfDefined_b(PyMOLGlobals *G,CSetting *set1,int index,int *value)
{
  int result=false;
  if(set1) {
    if(set1->info[index].defined) {
      *value=get_b(set1,index);
      result=true;
    }
  }
  return(result);
}
/*========================================================================*/
int SettingGetIfDefined_i(PyMOLGlobals *G,CSetting *set1,int index,int *value)
{
  int result=false;
  if(set1) {
    if(set1->info[index].defined) {
      *value=get_i(set1,index);
      result=true;
    }
  }
  return(result);
}
/*========================================================================*/
int   SettingGetIfDefined_color(PyMOLGlobals *G,CSetting *set1,int index,int *value)
{
  int result=false;
  if(set1) {
    if(set1->info[index].defined) {
      *value=get_color(set1,index);
      result=true;
    }
  }
  return(result);
}
/*========================================================================*/
int   SettingGetIfDefined_f(PyMOLGlobals *G,CSetting *set1,int index,float *value){
  int result=false;
  if(set1) {
    if(set1->info[index].defined) {
      *value=get_f(set1,index);
      result=true;
    }
  }
  return(result);
}
/*========================================================================*/
int   SettingGetIfDefined_3fv(PyMOLGlobals *G,CSetting *set1,int index,float **value)
{
  int result=false;
  if(set1) {
    if(set1->info[index].defined) {
      (*value) = (float*)(set1->data+set1->info[index].offset);
      result=true;
    }
  }
  return(result);
}
/*========================================================================*/
int   SettingGetIfDefined_s(PyMOLGlobals *G,CSetting *set1,int index,char **value)
{
  int result=false;
  if(set1) {
    if(set1->info[index].defined) {
      *value=get_s(set1,index);
      result=true;
    }
  }
  return(result);
}

/*========================================================================*/
int   SettingGet_i(PyMOLGlobals *G,CSetting *set1,CSetting *set2,int index)
{
  if(set1) {
    if(set1->info[index].defined) {
      return(get_i(set1,index));
    }
  }
  if(set2) { 
    if(set2->info[index].defined) {
      return(get_i(set2,index));
    }
  }
  return(SettingGetGlobal_i(G,index));
}
/*========================================================================*/
int   SettingGet_color(PyMOLGlobals *G,CSetting *set1,CSetting *set2,int index)
{
  if(set1) {
    if(set1->info[index].defined) {
      return(get_color(set1,index));
    }
  }
  if(set2) { 
    if(set2->info[index].defined) {
      return(get_color(set2,index));
    }
  }
  return(SettingGetGlobal_i(G,index));
}
/*========================================================================*/
float SettingGet_f(PyMOLGlobals *G,CSetting *set1,CSetting *set2,int index)
{
  if(set1) {
    if(set1->info[index].defined) {
      return(get_f(set1,index));
    }
  }
  if(set2) {
    if(set2->info[index].defined) {
      return(get_f(set2,index));
    }
  }
  return(SettingGetGlobal_f(G,index));
}
/*========================================================================*/
char *SettingGet_s(PyMOLGlobals *G,CSetting *set1,CSetting *set2,int index)
{
  if(set1) {
    if(set1->info[index].defined) {
      return(get_s(set1,index));
    }
  }
  if(set2) {
    if(set2->info[index].defined) {
      return(get_s(set2,index));
    }
  }
  return(SettingGetGlobal_s(G,index));
}
/*========================================================================*/
void  SettingGet_3f(PyMOLGlobals *G,CSetting *set1,CSetting *set2,int index,float *value)
{
  float *ptr;
  if(set1) {
    if(set1->info[index].defined) {
      ptr = (float*)(set1->data+set1->info[index].offset);
      copy3f(ptr,value);
      return;
    }
  }
  if(set2) {
    if(set2->info[index].defined) {
      ptr = (float*)(set2->data+set2->info[index].offset);
      copy3f(ptr,value);
      return;
    }
  }
  SettingGetGlobal_3f(G,index,value);
}
/*========================================================================*/
float *SettingGet_3fv(PyMOLGlobals *G,CSetting *set1,CSetting *set2,int index)
{
  if(set1) {
    if(set1->info[index].defined) {
      return (float*)(set1->data+set1->info[index].offset);
    }
  }
  if(set2) {
    if(set2->info[index].defined) {
      return (float*)(set2->data+set2->info[index].offset);
    }
  }
  return(SettingGetGlobal_3fv(G,index));
}

/*========================================================================*/
/*========================================================================*/
int SettingGetIndex(PyMOLGlobals *G,char *name) /* can be called from any thread state */
{
#ifdef _PYMOL_NOPY
  /* we're going to need a C-based dictionary and settings name list for this situation */
  return 0;
#else
  PyObject *tmp;
  int unblock;
  int index=-1; 
  
  unblock = PAutoBlock();
  if(P_setting) {
    tmp = PyObject_CallMethod(P_setting,"_get_index","s",name);
    if(tmp) {
      if(PyInt_Check(tmp))
        index = PyInt_AsLong(tmp);
      Py_DECREF(tmp);
    }
  }
  PAutoUnblock(unblock);

  return(index);
#endif
}
/*========================================================================*/
int SettingGetName(PyMOLGlobals *G,int index,SettingName name) /* can be called from any thread state */
{
#ifdef _PYMOL_NOPY
  /* we're going to need a C-based dictionary and settings name list for this situation */
  name[0]=0;
  return 0;
#else
  PyObject *tmp;
  int unblock;
  name[0]=0;
  unblock = PAutoBlock();
  if(P_setting) {
    tmp = PyObject_CallMethod(P_setting,"_get_name","i",index);
    if(tmp) {
      if(PyString_Check(tmp))
        UtilNCopy(name,PyString_AsString(tmp),sizeof(SettingName));
      Py_DECREF(tmp);
    }
  }
  PAutoUnblock(unblock);
  return(name[0]!=0);
#endif
}

/*========================================================================*/
void SettingGenerateSideEffects(PyMOLGlobals *G,int index,char *sele,int state)
{
  char all[] = "all";
  char *inv_sele;
  OrthoLineType command;
  if(!sele) {
    inv_sele = all;
  } else if(sele[0]==0) {
    inv_sele = all;
  } else {
    inv_sele = sele;
  }
  switch(index) {
  case cSetting_stereo:
    SceneUpdateStereo(G);
    break;
  case cSetting_pickable:
    ExecutiveInvalidateRep(G,inv_sele,cRepAll,cRepInvAll);
    SceneChanged(G);
    
    break;
  case cSetting_seq_view:
  case cSetting_seq_view_label_spacing:
  case cSetting_seq_view_label_mode:
  case cSetting_seq_view_label_start:
  case cSetting_seq_view_format:
  case cSetting_seq_view_color:
    SeqChanged(G);
    break;
  case cSetting_static_singletons:
    SeqChanged(G);
    break;
  case cSetting_seq_view_location:
  case cSetting_seq_view_overlay:
    PParse("cmd.viewport(-1,-1)");
    break;
  case cSetting_stereo_mode:
    SceneUpdateStereoMode(G);
    break;
  case cSetting_dot_lighting:
  case cSetting_mesh_lighting:
  case cSetting_light:
  case cSetting_fog:
  case cSetting_field_of_view:
  case cSetting_fog_start:
  case cSetting_two_sided_lighting:
  case cSetting_transparency_mode:
  case cSetting_dot_normals:
  case cSetting_mesh_normals:
	 SceneDirty(G);
	 break;
  case cSetting_stereo_shift:
  case cSetting_stereo_angle:
  case cSetting_texture_fonts:
	 SceneDirty(G);
	 break;
  case cSetting_dash_round_ends:
	 SceneDirty(G);
	 break;
  case cSetting_internal_gui_control_size:
    WizardRefresh(G);
    OrthoDirty(G);
    break;
  case cSetting_min_mesh_spacing:
  case cSetting_mesh_mode:
  case cSetting_mesh_type:
  case cSetting_mesh_solvent:
  case cSetting_mesh_quality:
    ExecutiveInvalidateRep(G,inv_sele,cRepMesh,cRepInvRep);
    SceneChanged(G);
    break;
  case cSetting_valence:
  case cSetting_half_bonds:
    ExecutiveInvalidateRep(G,inv_sele,cRepLine,cRepInvRep);
    ExecutiveInvalidateRep(G,inv_sele,cRepCyl,cRepInvRep);
    SceneChanged(G);
    break;
  case cSetting_slice_height_scale:
  case cSetting_slice_height_map:
  case cSetting_slice_grid:
  case cSetting_slice_dynamic_grid:
  case cSetting_slice_dynamic_grid_resolution:
    ExecutiveInvalidateRep(G,inv_sele,cRepSlice,cRepInvRep);    
    SceneChanged(G);
    break;
  case cSetting_label_font_id:
    SceneChanged(G);
    break;
  case cSetting_retain_order:
  case cSetting_pdb_hetatm_sort:
  case cSetting_pdb_insertions_go_first:
    ExecutiveSort(G,cKeywordAll);
    break;
  case cSetting_roving_lines:
  case cSetting_roving_sticks:
  case cSetting_roving_spheres:
  case cSetting_roving_labels:
  case cSetting_roving_selection:
  case cSetting_roving_ribbon:
  case cSetting_roving_cartoon:
  case cSetting_roving_polar_contacts:
  case cSetting_roving_polar_cutoff:
  case cSetting_roving_nonbonded:
  case cSetting_roving_nb_spheres:
  case cSetting_roving_map1_level:
  case cSetting_roving_map2_level:
  case cSetting_roving_map3_level:
  case cSetting_roving_map1_name:
  case cSetting_roving_map2_name:
  case cSetting_roving_map3_name:
  case cSetting_roving_isosurface:
  case cSetting_roving_isomesh:
    SceneRovingChanged(G);
    break;
  case cSetting_roving_byres:
  case cSetting_roving_detail:
    SceneRovingDirty(G);
    break;
  case cSetting_dash_length:
  case cSetting_dash_gap:
  case cSetting_dash_radius:
  case cSetting_dash_width:
    ExecutiveRebuildAllObjectDist(G);
    SceneChanged(G);
    break;
  case cSetting_button_mode:
    SceneDirty(G);
    break;
  case cSetting_stick_radius:
  case cSetting_stick_ball:
  case cSetting_stick_nub:
  case cSetting_stick_ball_ratio:
  case cSetting_stick_fixed_radius:
  case cSetting_stick_quality:
  case cSetting_stick_overlap:
  case cSetting_stick_color:
    ExecutiveInvalidateRep(G,inv_sele,cRepCyl,cRepInvRep);
    SceneChanged(G);
    break;
  case cSetting_clamp_colors:
    ExecutiveInvalidateRep(G,inv_sele,cRepAll,cRepInvColor);
    SceneChanged(G);
    break;
  case cSetting_label_color:
    ExecutiveInvalidateRep(G,inv_sele,cRepLabel,cRepInvRep);
    SceneChanged(G);
    break;
  case cSetting_cartoon_color:
    ExecutiveInvalidateRep(G,inv_sele,cRepCartoon,cRepInvRep);
    SceneChanged(G);
    break;
  case cSetting_ribbon_color:
    ExecutiveInvalidateRep(G,inv_sele,cRepRibbon,cRepInvRep);
    SceneChanged(G);
    break;
  case cSetting_all_states:
    SceneChanged(G);
    break;
  case cSetting_sel_counter:
	 break;
  case cSetting_line_width: /* auto-disable smooth lines if line width > 1 */
    /*    SettingSet(G,cSetting_line_smooth,0);  NO LONGER */
  case cSetting_line_radius:
    ExecutiveInvalidateRep(G,inv_sele,cRepLine,cRepInvRep);
    ExecutiveInvalidateRep(G,inv_sele,cRepNonbonded,cRepInvRep);
    SceneChanged(G);
    break;
  case cSetting_mesh_width: 
  case cSetting_mesh_color: 
    ExecutiveInvalidateRep(G,inv_sele,cRepMesh,cRepInvColor);
    SceneChanged(G);
    break;
  case cSetting_cull_spheres:
  case cSetting_sphere_scale:
  case cSetting_sphere_solvent:
    ExecutiveInvalidateRep(G,inv_sele,cRepSphere,cRepInvRep);
    SceneChanged(G);
    break;
  case cSetting_sphere_quality:
    ExecutiveInvalidateRep(G,inv_sele,cRepCyl,cRepInvRep);
    ExecutiveInvalidateRep(G,inv_sele,cRepNonbondedSphere,cRepInvRep);
    ExecutiveInvalidateRep(G,inv_sele,cRepSphere,cRepInvRep);
    SceneChanged(G);
    break;
  case cSetting_nonbonded_size:
    ExecutiveInvalidateRep(G,inv_sele,cRepNonbonded,cRepInvRep);
    ExecutiveInvalidateRep(G,inv_sele,cRepNonbondedSphere,cRepInvRep);
    SceneChanged(G);
    break;
  case cSetting_mesh_radius:
    ExecutiveInvalidateRep(G,inv_sele,cRepMesh,cRepInvColor);
    SceneChanged(G);
    break;
  case cSetting_surface_color:
  case cSetting_surface_ramp_above_mode:
    ExecutiveInvalidateRep(G,inv_sele,cRepSurface,cRepInvColor);
    SceneChanged(G);
    break;
  case cSetting_dot_color:
    ExecutiveInvalidateRep(G,inv_sele,cRepDot,cRepInvColor);
    SceneChanged(G);
    break;
  case cSetting_sphere_color:
    ExecutiveInvalidateRep(G,inv_sele,cRepSphere,cRepInvColor);
    SceneChanged(G);
    break;
  case cSetting_surface_quality:
  case cSetting_surface_mode:
  case cSetting_surface_type:
  case cSetting_surface_carve_state:
  case cSetting_surface_carve_cutoff:
  case cSetting_surface_carve_selection:
  case cSetting_surface_clear_state:
  case cSetting_surface_clear_cutoff:
  case cSetting_surface_clear_selection:
  case cSetting_surface_trim_cutoff:
  case cSetting_surface_trim_factor:
  case cSetting_surface_solvent:
  case cSetting_surface_proximity:
  case cSetting_cavity_cull:
    ExecutiveInvalidateRep(G,inv_sele,cRepSurface,cRepInvRep);
    SceneChanged(G);
    break;
  case cSetting_solvent_radius:
    ExecutiveInvalidateRep(G,inv_sele,cRepSurface,cRepInvRep);
    ExecutiveInvalidateRep(G,inv_sele,cRepMesh,cRepInvRep);
    ExecutiveInvalidateRep(G,inv_sele,cRepDot,cRepInvRep);
    SceneChanged(G);
    break;
  case cSetting_ribbon_power:
  case cSetting_ribbon_power_b:
  case cSetting_ribbon_sampling:
  case cSetting_ribbon_radius:
  case cSetting_ribbon_width:
  case cSetting_ribbon_throw:
  case cSetting_ribbon_trace:
    ExecutiveInvalidateRep(G,inv_sele,cRepRibbon,cRepInvRep);
    SceneChanged(G);
    break;
  case cSetting_cartoon_side_chain_helper:
    ExecutiveInvalidateRep(G,inv_sele,cRepCartoon,cRepInvRep);
    ExecutiveInvalidateRep(G,inv_sele,cRepLine,cRepInvRep);
    ExecutiveInvalidateRep(G,inv_sele,cRepCyl,cRepInvRep);
    SceneChanged(G);
    break;
  case cSetting_cartoon_transparency:
  case cSetting_cartoon_trace:
  case cSetting_cartoon_refine:
  case cSetting_cartoon_nucleic_acid_mode:
  case cSetting_cartoon_sampling:
  case cSetting_cartoon_loop_quality:
  case cSetting_cartoon_loop_radius:
  case cSetting_cartoon_tube_quality:
  case cSetting_cartoon_tube_radius:
  case cSetting_cartoon_putty_quality:
  case cSetting_cartoon_putty_radius:
  case cSetting_cartoon_putty_scale_min:
  case cSetting_cartoon_putty_scale_max:
  case cSetting_cartoon_putty_scale_power:
  case cSetting_cartoon_power:
  case cSetting_cartoon_power_b:
  case cSetting_cartoon_rect_length:
  case cSetting_cartoon_rect_width:
  case cSetting_cartoon_oval_length:
  case cSetting_cartoon_oval_width:
  case cSetting_cartoon_oval_quality:
  case cSetting_cartoon_round_helices:
  case cSetting_cartoon_flat_sheets:
  case cSetting_cartoon_refine_normals:
  case cSetting_cartoon_smooth_loops:
  case cSetting_cartoon_dumbbell_width:
  case cSetting_cartoon_dumbbell_length:
  case cSetting_cartoon_dumbbell_radius:
  case cSetting_cartoon_fancy_helices:
  case cSetting_cartoon_fancy_sheets:
  case cSetting_cartoon_cylindrical_helices:
  case cSetting_cartoon_refine_tips:
  case cSetting_cartoon_helix_radius:
  case cSetting_cartoon_throw:
  case cSetting_cartoon_debug:
  case cSetting_cartoon_highlight_color:
  case cSetting_cartoon_discrete_colors:
  case cSetting_cartoon_smooth_first:
  case cSetting_cartoon_smooth_last:
  case cSetting_cartoon_smooth_cycles:
  case cSetting_cartoon_flat_cycles:

    ExecutiveInvalidateRep(G,inv_sele,cRepCartoon,cRepInvRep);
    SceneChanged(G);
    break;
  case cSetting_dot_width:
  case cSetting_dot_radius:
  case cSetting_dot_density:
  case cSetting_dot_mode:
  case cSetting_dot_solvent:
  case cSetting_dot_hydrogens:
  case cSetting_trim_dots:
    ExecutiveInvalidateRep(G,inv_sele,cRepDot,cRepInvRep);
    SceneChanged(G);
    break;
  case cSetting_line_smooth:
  case cSetting_transparency:
  case cSetting_sphere_transparency:
  case cSetting_stick_transparency:
  case cSetting_ortho:
  case cSetting_gl_ambient:
  case cSetting_bg_rgb:
  case cSetting_specular:
  case cSetting_specular_intensity:
  case cSetting_cgo_line_width:
  case cSetting_selection_width:
    SceneDirty(G);
    break;
  case cSetting_depth_cue: 
    SceneDirty(G);
    break;
  case cSetting_sculpting:
    OrthoDirty(G);
    break;
  case cSetting_overlay:
  case cSetting_overlay_lines:
  case cSetting_text:
    OrthoDirty(G);
    break;
  case cSetting_internal_gui_mode:
    break;
  case cSetting_internal_gui_width:
  case cSetting_internal_gui:
  case cSetting_internal_feedback:
    sprintf(command,"viewport");
    OrthoCommandIn(G,command);
    break;
  case cSetting_suspend_updates:
    if(!SettingGet(G,cSetting_suspend_updates))
      SceneChanged(G); /* force big update upon resumption */
    break;
  case cSetting_security:
    G->Security = (int)SettingGet(G,cSetting_security);
    break;
  case cSetting_state:
  case cSetting_frame:
    SceneChanged(G);
    break;
  default:
	 break;
  }
}
/*========================================================================*/
int SettingSetfv(PyMOLGlobals *G,int index,float *v)
{
  /* Warren, are these side effects still relevant? */

  register CSetting *I=G->Setting;
  int ok=true;
  switch(index) {
  case cSetting_dot_mode:
    SettingSet_f(I,index,v[0]);
    /*I->Setting[index].Value[0]=v[0];*/
	 break;
  case cSetting_bg_rgb:
  case cSetting_light:
    SettingSet_3fv(I,index,v); 
    /*
      I->Setting[index].Value[0]=v[0];
      I->Setting[index].Value[1]=v[1];
      I->Setting[index].Value[2]=v[2];
    */
	 SceneDirty(G);
	 break;
  case cSetting_valence:
    ExecutiveInvalidateRep(G,"all",cRepLine,cRepInvRep);
    SettingSet_f(I,index,v[0]);
    SceneChanged(G);
    break;
  case cSetting_dash_length:
  case cSetting_dash_gap:
    ExecutiveInvalidateRep(G,"all",cRepDash,cRepInvRep);
    SettingSet_f(I,index,v[0]);
    SceneChanged(G);
    break;
  case cSetting_button_mode:
    SettingSet_f(I,index,v[0]);
    SceneDirty(G);
    break;
  case cSetting_stick_radius:
  case cSetting_stick_quality:
  case cSetting_stick_overlap:
    ExecutiveInvalidateRep(G,"all",cRepCyl,cRepInvRep);
    SettingSet_f(I,index,v[0]);
    /*I->Setting[index].Value[0]=v[0];   */
    SceneChanged(G);
    break;
  case cSetting_label_color:
    ExecutiveInvalidateRep(G,"all",cRepLabel,cRepInvRep);
    SettingSet_f(I,index,v[0]);
    /* I->Setting[index].Value[0]=v[0]; */
    SceneChanged(G);
    break;
  case cSetting_all_states:
    SettingSet_f(I,index,v[0]);
    /* I->Setting[index].Value[0]=v[0];  */
    SceneChanged(G);
    break;
  case cSetting_dot_density:
    SettingSet_f(I,index,v[0]);
    /*I->Setting[index].Value[0]=v[0];*/
	 break;
  case cSetting_sel_counter:
    SettingSet_f(I,index,v[0]);
    /*I->Setting[index].Value[0]=v[0];*/
	 break;
  case cSetting_ortho:
  case cSetting_gl_ambient:
	 SceneDirty(G);
    break;
  case cSetting_overlay:
  case cSetting_text:
    OrthoDirty(G);
  default:
    ok = SettingSet_f(I,index,v[0]);
    /*I->Setting[index].Value[0]=v[0];*/
	 break;
  }
  return(ok);
}
/*========================================================================*/
int   SettingSetGlobal_b(PyMOLGlobals *G,int index,int value) 
{
  return(SettingSet_b(G->Setting,index,value));
}
/*========================================================================*/
int   SettingSetGlobal_i(PyMOLGlobals *G,int index,int value) 
{
  return(SettingSet_i(G->Setting,index,value));
}
/*========================================================================*/
int   SettingSetGlobal_f(PyMOLGlobals *G,int index,float value) 
{
  return(SettingSet_f(G->Setting,index,value));
}
/*========================================================================*/
int SettingSet(PyMOLGlobals *G,int index,float v)
{
  return(SettingSetfv(G,index,&v));
}
/*========================================================================*/
int SettingSetNamed(PyMOLGlobals *G,char *name,char *value)
{
  int ok=true;
  int index = SettingGetIndex(G,name);
  float v,vv[3];
  SettingName realName;
  char buffer[1024] = "";
  if(index>=0) {
    SettingGetName(G,index,realName);
	 switch(index) {
	 case cSetting_dot_mode:
		if(strcmp(value,"molecular")==0) {
		  v=0.0;
		  SettingSetfv(G,index,&v);
		  sprintf(buffer," Setting: %s set to %s\n",realName,value);
		} else if(strcmp(value,"solvent_accessible")==0) {
		  v=1.0;
		  SettingSetfv(G,index,&v);
		  sprintf(buffer," Setting: %s set to %s\n",realName,value);
		} else if(sscanf(value,"%f",&v)==1) {
		  SettingSetfv(G,index,&v);
		  sprintf(buffer," Setting: %s set to %s\n",realName,value);
      }
		break;
	 case cSetting_bg_rgb:
	 case cSetting_light:
		if(sscanf(value,"%f%f%f",vv,vv+1,vv+2)==3) {
		  SettingSetfv(G,index,vv);
		  sprintf(buffer," Setting: %s set to %5.3f %8.3f %8.3f\n",realName,
					 *vv,*(vv+1),*(vv+2));
		}
		break;
	 case cSetting_dot_density:
		sscanf(value,"%f",&v);
		SettingSetfv(G,index,&v);
		sprintf(buffer," Setting: %s set to %d\n",realName,(int)v);
		break;
	 case cSetting_text:
	 case cSetting_overlay:
	 case cSetting_sel_counter:
    case cSetting_dist_counter:
		sscanf(value,"%f",&v);
		SettingSetfv(G,index,&v);
		break;
    case cSetting_line_width: /* auto-disable smooth lines if line width > 1 */
    case cSetting_mesh_width:
      sscanf(value,"%f",&v);
		SettingSetfv(G,index,&v);
		sprintf(buffer," Setting: %s set to %5.3f\n",realName,v);
      SceneDirty(G);
      break;
	 default:
		sscanf(value,"%f",&v);
		SettingSetfv(G,index,&v);
		sprintf(buffer," Setting: %s set to %5.3f\n",realName,v);
		break;
	 }
  } else {
    PRINTFB(G,FB_Setting,FB_Warnings)
      " Error: Non-Existent Settin\n"
      ENDFB(G);
    ok=false;
  }
  if(buffer[0]) {
    PRINTFB(G,FB_Setting,FB_Actions)
      "%s",buffer
      ENDFB(G);
  }
  return(ok);
}
/*========================================================================*/
float SettingGetNamed(PyMOLGlobals *G,char *name)
{
  return(SettingGet(G,SettingGetIndex(G,name)));
}
/*========================================================================*/
float SettingGet(PyMOLGlobals *G,int index)
{
  return(SettingGetGlobal_f(G,index));
}
/*========================================================================*/
float *SettingGetfv(PyMOLGlobals *G,int index)
{
  return(SettingGetGlobal_3fv(G,index));
}
/*========================================================================*/
void SettingFreeGlobal(PyMOLGlobals *G)
{
  register CSetting *I=G->Setting;
  SettingPurge(I);
  FreeP(G->Setting);
}
/*========================================================================*/
void SettingInitGlobal(PyMOLGlobals *G,int alloc,int reset_gui)
{
  register CSetting *I=G->Setting;

  if(alloc || !I) {
    I=(G->Setting=Calloc(CSetting,1));
    SettingInit(G,I);
  }

  SettingSet_f(I,cSetting_bonding_vdw_cutoff, 0.2F);

  SettingSet_f(I,cSetting_min_mesh_spacing, 0.6F);

  SettingSet_i(I,cSetting_dot_density, 2);

  SettingSet_i(I,cSetting_dot_mode, 0);

  SettingSet_f(I,cSetting_solvent_radius, 1.4F);

  SettingSet_i(I,cSetting_sel_counter, 0);

  SettingSet_3f(I,cSetting_bg_rgb, 0.0F, 0.0F, 0.0F);

  SettingSet_f(I,cSetting_ambient, 0.12F);

  SettingSet_f(I,cSetting_direct, 0.25F);

  SettingSet_f(I,cSetting_reflect, 0.90F);

  SettingSet_3f(I,cSetting_light, -0.4F, -0.4F, -1.0F);

  SettingSet_i(I,cSetting_antialias, 1);

  SettingSet_i(I,cSetting_cavity_cull, 10);

  SettingSet_f(I,cSetting_gl_ambient,  0.12F);

  SettingSet_b(I,cSetting_single_image, 0);

  SettingSet_f(I,cSetting_movie_delay, 30.0F);

  SettingSet_f(I,cSetting_ribbon_power, 2.0F);

  SettingSet_f(I,cSetting_ribbon_power_b, 0.5F);

  SettingSet_i(I,cSetting_ribbon_sampling, 1);

  SettingSet_f(I,cSetting_ribbon_radius, 0.0F);

  SettingSet_f(I,cSetting_stick_radius, 0.25F);

  SettingSet_i(I,cSetting_hash_max, 100);

  SettingSet_b(I,cSetting_ortho, 0);

  SettingSet_f(I,cSetting_power, 1.0F);

  SettingSet_f(I,cSetting_spec_reflect, 0.5F);

  SettingSet_f(I,cSetting_spec_power, 60.0F);

  SettingSet_f(I,cSetting_sweep_angle, 20.0F);

  SettingSet_f(I,cSetting_sweep_speed, 0.75F);

  SettingSet_b(I,cSetting_dot_hydrogens, 1);

  SettingSet_f(I,cSetting_dot_radius, 0.0F);

  SettingSet_b(I,cSetting_ray_trace_frames, 0);

  SettingSet_b(I,cSetting_cache_frames, 0);

  SettingSet_b(I,cSetting_trim_dots, 1);

  SettingSet_b(I,cSetting_cull_spheres, 1);

  SettingSet_f(I,cSetting_test1, -56.6F);

  SettingSet_f(I,cSetting_test2, 1.0F);

  SettingSet_f(I,cSetting_surface_best, 0.25F);

  SettingSet_f(I,cSetting_surface_normal, 0.5F);

  SettingSet_i(I,cSetting_surface_quality, 0);

  SettingSet_b(I,cSetting_surface_proximity, 1);

  SettingSet_f(I,cSetting_stereo_angle, 2.1F);

  SettingSet_f(I,cSetting_stereo_shift, 2.0F);

  SettingSet_b(I,cSetting_line_smooth, 1);

  SettingSet_f(I,cSetting_line_width, 1.0F);

  SettingSet_b(I,cSetting_half_bonds, 0);

  SettingSet_i(I,cSetting_stick_quality, 8);

  SettingSet_f(I,cSetting_stick_overlap, 0.2F);

  SettingSet_f(I,cSetting_stick_nub, 0.7F);

  SettingSet_b(I,cSetting_all_states, 0);

  SettingSet_b(I,cSetting_pickable, 1);

  SettingSet_b(I,cSetting_auto_show_lines, 1);

  SettingSet_f(I,cSetting_fast_idle, 20000.0F); /* 1/50 of a sec. */

  SettingSet_f(I,cSetting_no_idle, 10000.0F); /* 1/100 of a sec. */ 

#ifdef _PYMOL_OSX
  SettingSet_f(I,cSetting_slow_idle, 40000.0F); /* 1/25 of a sec. */
#else
  SettingSet_f(I,cSetting_slow_idle, 200000.0F); /* 1/5 of a sec. */
#endif

  SettingSet_f(I,cSetting_idle_delay, 1.5F);

  SettingSet_f(I,cSetting_rock_delay, 30.0F);

  SettingSet_i(I,cSetting_dist_counter, 0);

  SettingSet_f(I,cSetting_dash_length, 0.15F);

  SettingSet_f(I,cSetting_dash_gap, 0.35F);

  SettingSet_b(I,cSetting_auto_zoom, 1);

  SettingSet_i(I,cSetting_overlay, 0);

  SettingSet_i(I,cSetting_text, 0);

  SettingSet_i(I,cSetting_button_mode, 0);

  SettingSet_b(I,cSetting_valence, 0);

  SettingSet_f(I,cSetting_nonbonded_size, 0.25F);

  SettingSet_color(I,cSetting_label_color, "-1");

  SettingSet_f(I,cSetting_ray_trace_fog, -1.0F);

  SettingSet_f(I,cSetting_spheroid_scale, 1.0F);

  SettingSet_f(I,cSetting_ray_trace_fog_start, 0.50F);

  SettingSet_f(I,cSetting_spheroid_smooth, 1.1F);

  SettingSet_f(I,cSetting_spheroid_fill, 1.30F);

  SettingSet_b(I,cSetting_auto_show_nonbonded, 1);

  SettingSet_f(I,cSetting_mesh_radius, 0.000F);
 
#ifdef WIN32
  SettingSet_b(I,cSetting_cache_display, 0);
#else
  SettingSet_b(I,cSetting_cache_display, 1);
#endif

  SettingSet_b(I,cSetting_normal_workaround, 0);

  SettingSet_b(I,cSetting_backface_cull, 1);

  SettingSet_f(I,cSetting_gamma, 1.3F);

  SettingSet_f(I,cSetting_dot_width, 2.0F);

  SettingSet_b(I,cSetting_auto_show_selections, 1);

  SettingSet_b(I,cSetting_auto_hide_selections, 1);

  SettingSet_f(I,cSetting_selection_width, 4.0F);

  SettingSet_f(I,cSetting_selection_overlay, 1.0F);

  SettingSet_b(I,cSetting_static_singletons, 1);

  SettingSet_i(I,cSetting_max_triangles, 1000000); /* no longer used */

  SettingSet_b(I,cSetting_depth_cue, 1);

  SettingSet_f(I,cSetting_specular, 1.0F);

  SettingSet_f(I,cSetting_shininess, 40.0F);

  SettingSet_i(I,cSetting_sphere_quality, 1);

  SettingSet_f(I,cSetting_fog, 1.0F);

  SettingSet_b(I,cSetting_isomesh_auto_state, 0); /* no longer necessary? */

  SettingSet_f(I,cSetting_mesh_width, 1.0F);

  SettingSet_i(I,cSetting_cartoon_sampling, 7);

  SettingSet_f(I,cSetting_cartoon_loop_radius, 0.2F);

  SettingSet_f(I,cSetting_cartoon_loop_quality, 6.0F);

  SettingSet_f(I,cSetting_cartoon_power, 2.0F);

  SettingSet_f(I,cSetting_cartoon_power_b, 0.52F);

  SettingSet_f(I,cSetting_cartoon_rect_length, 1.40F);

  SettingSet_f(I,cSetting_cartoon_rect_width, 0.4F);

  if(reset_gui) {
    SettingSet_i(I,cSetting_internal_gui_width, cOrthoRightSceneMargin);

    SettingSet_b(I,cSetting_internal_gui, 1);
  }

  SettingSet_f(I,cSetting_cartoon_oval_length, 1.35F);

  SettingSet_f(I,cSetting_cartoon_oval_width, 0.25F);

  SettingSet_f(I,cSetting_cartoon_oval_quality, 10.0F);

  SettingSet_f(I,cSetting_cartoon_tube_radius, 0.5F);

  SettingSet_f(I,cSetting_cartoon_tube_quality, 9.0F);

  SettingSet_i(I,cSetting_cartoon_debug, 0);

  SettingSet_f(I,cSetting_ribbon_width, 3.0F);

  SettingSet_f(I,cSetting_dash_width, 3.0F);

  SettingSet_f(I,cSetting_dash_radius, 0.00F);

  SettingSet_f(I,cSetting_cgo_ray_width_scale, -0.15F);

  SettingSet_f(I,cSetting_line_radius, 0.0F);

  SettingSet_b(I,cSetting_cartoon_round_helices, 1);

  SettingSet_b(I,cSetting_cartoon_refine_normals, 1);
  
  SettingSet_b(I,cSetting_cartoon_flat_sheets, 1);

  SettingSet_b(I,cSetting_cartoon_smooth_loops, 0);

  SettingSet_f(I,cSetting_cartoon_dumbbell_length, 1.60F);

  SettingSet_f(I,cSetting_cartoon_dumbbell_width, 0.17F);

  SettingSet_f(I,cSetting_cartoon_dumbbell_radius, 0.16F);

  SettingSet_b(I,cSetting_cartoon_fancy_helices, 0);  

  SettingSet_b(I,cSetting_cartoon_fancy_sheets, 1);  

  SettingSet_b(I,cSetting_ignore_pdb_segi, 0);  

  SettingSet_f(I,cSetting_ribbon_throw, 1.35F);  

  SettingSet_f(I,cSetting_cartoon_throw, 1.35F);  

  SettingSet_i(I,cSetting_cartoon_refine, 5);  

  SettingSet_i(I,cSetting_cartoon_refine_tips, 10);  

  SettingSet_b(I,cSetting_cartoon_discrete_colors, 0);  

  SettingSet_b(I,cSetting_normalize_ccp4_maps, 1);  

  SettingSet_f(I,cSetting_surface_poor, 0.85F);  

  SettingSet_i(I,cSetting_internal_feedback, G->Option->internal_feedback);

  SettingSet_f(I,cSetting_cgo_line_width, 1.00F);

  SettingSet_f(I,cSetting_cgo_line_radius, -0.05F);

  SettingSet_i(I,cSetting_logging, 0); /* 0 = off, 1 = regular (PML), 2 = python (PYM) */

  SettingSet_b(I,cSetting_robust_logs, 0);

  SettingSet_b(I,cSetting_log_box_selections, 1);

  SettingSet_b(I,cSetting_log_conformations, 1);

  SettingSet_f(I,cSetting_valence_size, 0.050F);

  SettingSet_f(I,cSetting_surface_miserable, 2.0F);

  SettingSet_b(I,cSetting_ray_opaque_background, 1);

  SettingSet_f(I,cSetting_transparency, 0.0F);

  SettingSet_i(I,cSetting_ray_texture, 0);

  SettingSet_3f(I,cSetting_ray_texture_settings, 0.1F, 5.0F, 1.0F);

  SettingSet_b(I,cSetting_suspend_updates, 0);

  SettingSet_b(I,cSetting_full_screen, 0);

  SettingSet_i(I,cSetting_surface_mode, 0); /* by flag is the default */

  SettingSet_color(I,cSetting_surface_color,"-1"); /* use atom colors by default */

  SettingSet_i(I,cSetting_mesh_mode,0); /* by flag is the default */

  SettingSet_color(I,cSetting_mesh_color,"-1"); /* use atom colors by default */

  SettingSet_b(I,cSetting_auto_indicate_flags,0); 

  SettingSet_i(I,cSetting_surface_debug,0);

  SettingSet_f(I,cSetting_ray_improve_shadows,0.1F);

  SettingSet_b(I,cSetting_smooth_color_triangle,0);

  SettingSet_i(I,cSetting_ray_default_renderer,0);

  SettingSet_f(I,cSetting_field_of_view,20.0F);

  SettingSet_f(I,cSetting_reflect_power, 1.0F);

  SettingSet_b(I,cSetting_preserve_chempy_ids, 0);

  SettingSet_f(I,cSetting_sphere_scale, 1.0F);

  SettingSet_b(I,cSetting_two_sided_lighting, 0);

  SettingSet_f(I,cSetting_secondary_structure, 2.0F); /* unused? */

  SettingSet_b(I,cSetting_auto_remove_hydrogens, 0);

  SettingSet_b(I,cSetting_raise_exceptions, 1);

  SettingSet_b(I,cSetting_stop_on_exceptions, 0);  

  SettingSet_b(I,cSetting_sculpting, 0);  

  SettingSet_b(I,cSetting_auto_sculpt, 0);  

  SettingSet_f(I,cSetting_sculpt_vdw_scale, 0.97F);  

  SettingSet_f(I,cSetting_sculpt_vdw_scale14, 0.915F);  

  SettingSet_f(I,cSetting_sculpt_vdw_weight, 1.0F);  

  SettingSet_f(I,cSetting_sculpt_vdw_weight14, 0.33F);  

  SettingSet_f(I,cSetting_sculpt_bond_weight, 2.25F);  

  SettingSet_f(I,cSetting_sculpt_angl_weight, 1.0F);  

  SettingSet_f(I,cSetting_sculpt_pyra_weight, 1.0F);  

  SettingSet_f(I,cSetting_sculpt_plan_weight, 1.0F);  

  SettingSet_i(I,cSetting_sculpting_cycles, 10);  

  SettingSet_f(I,cSetting_sphere_transparency, 0.0F);

  SettingSet_color(I,cSetting_sphere_color,"-1"); /* use atom colors by default */

  SettingSet_f(I,cSetting_sculpt_field_mask, (float)0xFF );  /* all terms */

  SettingSet_f(I,cSetting_sculpt_hb_overlap, 1.0F);

  SettingSet_f(I,cSetting_sculpt_hb_overlap_base, 0.35F);

  SettingSet_b(I,cSetting_legacy_vdw_radii, 0);  

  SettingSet_b(I,cSetting_sculpt_memory, 1);
  
  SettingSet_i(I,cSetting_connect_mode, 0);

  SettingSet_b(I,cSetting_cartoon_cylindrical_helices, 0);

  SettingSet_f(I,cSetting_cartoon_helix_radius, 2.25F);

  SettingSet_f(I,cSetting_connect_cutoff, 0.35F);

  SettingSet_b(I,cSetting_save_pdb_ss, 0);

  SettingSet_f(I,cSetting_sculpt_line_weight, 1.0F);

  SettingSet_i(I,cSetting_fit_iterations, 1000);

  SettingSet_f(I,cSetting_fit_tolerance, 0.00001F);

  SettingSet_s(I,cSetting_batch_prefix,"tmp_pymol");

  if(G->StereoCapable) {
    SettingSet_i(I,cSetting_stereo_mode, 1); 
  } else {
    SettingSet_i(I,cSetting_stereo_mode, 2);
  }

  SettingSet_i(I,cSetting_cgo_sphere_quality, 1);

  SettingSet_b(I,cSetting_pdb_literal_names, 0);

  SettingSet_b(I,cSetting_wrap_output, 0);

  SettingSet_f(I,cSetting_fog_start, 0.45F);

  SettingSet_i(I,cSetting_frame,1);

  SettingSet_i(I,cSetting_state,1);

  SettingSet_i(I,cSetting_ray_shadows,1);

  SettingSet_i(I,cSetting_ribbon_trace,0);

  SettingSet_i(I,cSetting_security,1);

  SettingSet_f(I,cSetting_stick_transparency,0.0F);

  SettingSet_i(I,cSetting_ray_transparency_shadows,1);

  SettingSet_i(I,cSetting_session_version_check,0);

  SettingSet_f(I,cSetting_ray_transparency_specular,0.4F);

  SettingSet_b(I,cSetting_stereo_double_pump_mono,0);

  SettingSet_i(I,cSetting_sphere_solvent,0);

  SettingSet_i(I,cSetting_mesh_quality,2);

  SettingSet_i(I,cSetting_mesh_solvent,0);
  
  SettingSet_b(I,cSetting_dot_solvent,0);

  SettingSet_f(I,cSetting_ray_shadow_fudge,0.001F);

  SettingSet_f(I,cSetting_ray_triangle_fudge,0.0000001F);

  SettingSet_i(I,cSetting_debug_pick,0);

  SettingSet_color(I,cSetting_dot_color,"-1"); /* use atom colors by default */

  SettingSet_f(I,cSetting_mouse_limit,100.0F);

  SettingSet_f(I,cSetting_mouse_scale,1.3F);

  SettingSet_i(I,cSetting_transparency_mode,2);

  SettingSet_b(I,cSetting_clamp_colors,1);

  SettingSet_f(I,cSetting_pymol_space_max_red,0.90F);

  SettingSet_f(I,cSetting_pymol_space_max_green,0.75F);

  SettingSet_f(I,cSetting_pymol_space_max_blue,0.90F);

  SettingSet_f(I,cSetting_pymol_space_min_factor,0.15F);

  SettingSet_b(I,cSetting_roving_origin,1);

  SettingSet_f(I,cSetting_roving_sticks,6.0F);

  SettingSet_f(I,cSetting_roving_lines,10.0F);

  SettingSet_f(I,cSetting_roving_spheres,0.0F);

  SettingSet_f(I,cSetting_roving_labels,0.0F);

  SettingSet_f(I,cSetting_roving_delay,0.2F);

  SettingSet_s(I,cSetting_roving_selection,"all");

  SettingSet_b(I,cSetting_roving_byres,1);

  SettingSet_f(I,cSetting_roving_ribbon,-7.0F);
  
  SettingSet_f(I,cSetting_roving_cartoon,0.0F);

  SettingSet_f(I,cSetting_roving_polar_contacts,7.0F);

  SettingSet_f(I,cSetting_roving_polar_cutoff,3.31F);

  SettingSet_f(I,cSetting_roving_nonbonded,0.0F);

  SettingSet_i(I,cSetting_float_labels,0);

  SettingSet_b(I,cSetting_roving_detail,0);

  SettingSet_f(I,cSetting_roving_nb_spheres,8.0F);

  SettingSet_color(I,cSetting_ribbon_color,"-1"); /* use atom colors by default */

  SettingSet_color(I,cSetting_cartoon_color,"-1"); /* use atom colors by default */

  SettingSet_i(I,cSetting_ribbon_smooth,-1);

  SettingSet_b(I,cSetting_auto_color,1);

  SettingSet_i(I,cSetting_auto_color_next,0);

  SettingSet_color(I,cSetting_ray_interior_color,"-1"); /* no color */

  SettingSet_color(I,cSetting_cartoon_highlight_color,"-1"); /* no color */

  /* In order to get electrostatic potentials in kT from the Coulomb equation... 

  PyMOL charges: Q, in e
  PyMOL distances: r, in Angstrom
  Coulomb Constant: K = 8.987552e9 ((J*m)/(C^2))
  Angstrom Convertor: 1 A = 1e-10 m
  Coulomb Convertor: 1 e = 1.60217733e-19 C
  Angstrom Convertor: 1 A = 10e-10 m
  Dielectric Constant: D (unitless)

  ePot = (KQ)/(Dr) = 

  8.987552e9 J*m     1.6021773e-19 C   1.6021773e-19 C      1 A        Q  
  ---------------- * --------------- * --------------- * ---------- * --- =
         C^2              1 e               1 e            1e-10 m     Dr


    2.3070795237e-18 J*A      Q
  = ---------------------- * ---
            e^2               Dr

  Boltzmann Constant: k = 1.380658e-23 (J/K)
  Temperature: 300 Kelvin

  kT = 1.380658e-23 * 300 = 4.141974e-21 J

         2.3070795237e-18 J*A         1 kT             Q
  ePot = --------------------- * ------------------ * ---
                 e^2              4.141974e-21 J       Dr
  
         557.00000 kT*A    Q
  ePot = -------------- * --- which will give kT/e units when applied
             e^2           Dr
  */

  SettingSet_f(I,cSetting_coulomb_units_factor,557.00000F);

  SettingSet_f(I,cSetting_coulomb_dielectric,2.0F);

  SettingSet_b(I,cSetting_ray_interior_shadows,0);

  SettingSet_i(I,cSetting_ray_interior_texture,-1);

  SettingSet_s(I,cSetting_roving_map1_name,"");

  SettingSet_s(I,cSetting_roving_map2_name,"");

  SettingSet_s(I,cSetting_roving_map3_name,"");

  SettingSet_f(I,cSetting_roving_map1_level,1.0F);

  SettingSet_f(I,cSetting_roving_map2_level,2.0F);

  SettingSet_f(I,cSetting_roving_map3_level,3.0F);

  SettingSet_f(I,cSetting_roving_isomesh, 8.0F);

  SettingSet_f(I,cSetting_roving_isosurface, 0.0F);

  SettingSet_f(I,cSetting_scenes_changed, 1.0F);

  SettingSet_f(I,cSetting_gaussian_b_adjust, 0.0F );

  SettingSet_b(I,cSetting_pdb_standard_order, 1);

  SettingSet_i(I,cSetting_cartoon_smooth_first, 1);
  SettingSet_i(I,cSetting_cartoon_smooth_last, 1);
  SettingSet_i(I,cSetting_cartoon_smooth_cycles,2);
  SettingSet_i(I,cSetting_cartoon_flat_cycles,4);

  SettingSet_i(I,cSetting_max_threads, 1);

  SettingSet_i(I,cSetting_show_progress, 1);

  SettingSet_i(I,cSetting_use_display_lists, 0); /* don't make this default
                                                    until we have a way of
                                                    reusing display list 
                                                    identifiers */

  SettingSet_i(I,cSetting_cache_memory, 0);  /* doesn't seem to do any good :( */
  
  SettingSet_i(I,cSetting_simplify_display_lists, 0);

  SettingSet_i(I,cSetting_retain_order,0);

  SettingSet_i(I,cSetting_pdb_hetatm_sort,0);

  SettingSet_i(I,cSetting_pdb_use_ter_records,1);

  SettingSet_i(I,cSetting_cartoon_trace,0);

  SettingSet_i(I,cSetting_ray_oversample_cutoff,64);

  SettingSet_f(I,cSetting_gaussian_resolution,2.0F);

  SettingSet_f(I,cSetting_gaussian_b_floor,0.0F);

  SettingSet_i(I,cSetting_sculpt_nb_interval,17);

  SettingSet_f(I,cSetting_sculpt_tors_weight,0.05F);

  SettingSet_f(I,cSetting_sculpt_tors_tolerance,0.05F);

  SettingSet_b(I,cSetting_stick_ball, false);

  SettingSet_f(I,cSetting_stick_ball_ratio, 1.0F);

  SettingSet_b(I,cSetting_stick_fixed_radius, false);

  SettingSet_f(I,cSetting_cartoon_transparency, 0.0F);

  SettingSet_b(I,cSetting_dash_round_ends,1);

  SettingSet_f(I,cSetting_h_bond_max_angle,63.0F);

  SettingSet_f(I,cSetting_h_bond_cutoff_center,3.6F);

  SettingSet_f(I,cSetting_h_bond_cutoff_edge,3.2F);

  SettingSet_f(I,cSetting_h_bond_power_a,1.6F);

  SettingSet_f(I,cSetting_h_bond_power_b,5.0F);

  SettingSet_f(I,cSetting_h_bond_cone, 180.0F);

  SettingSet_f(I,cSetting_ss_helix_psi_target, -48.0F);
  SettingSet_f(I,cSetting_ss_helix_psi_include, 55.0F); /* 30 */
  SettingSet_f(I,cSetting_ss_helix_psi_exclude, 85.0F);
  
  SettingSet_f(I,cSetting_ss_helix_phi_target, -57.0F);
  SettingSet_f(I,cSetting_ss_helix_phi_include, 55.0F);
  SettingSet_f(I,cSetting_ss_helix_phi_exclude, 85.0F);
  
  SettingSet_f(I,cSetting_ss_strand_psi_target, 124.0F); 
  SettingSet_f(I,cSetting_ss_strand_psi_include, 40.0F); 
  SettingSet_f(I,cSetting_ss_strand_psi_exclude, 90.0F); /* 80 */
  
  SettingSet_f(I,cSetting_ss_strand_phi_target, -129.0F); 
  SettingSet_f(I,cSetting_ss_strand_phi_include, 40.0F); 
  SettingSet_f(I,cSetting_ss_strand_phi_exclude, 100.0F);

  SettingSet_b(I,cSetting_movie_loop, 1);
     
  SettingSet_b(I,cSetting_pdb_retain_ids, 0);

  SettingSet_b(I,cSetting_pdb_no_end_record, 0);

  SettingSet_f(I,cSetting_cgo_dot_width, 2.0F);
  SettingSet_f(I,cSetting_cgo_dot_radius, -1.0F);
  SettingSet_b(I,cSetting_defer_updates,0);
  SettingSet_b(I,cSetting_normalize_o_maps,1);
  SettingSet_b(I,cSetting_swap_dsn6_bytes,1);
  SettingSet_b(I,cSetting_pdb_insertions_go_first,0);
  SettingSet_b(I,cSetting_roving_origin_z,1);
  SettingSet_f(I,cSetting_roving_origin_z_cushion,3.0F);
  SettingSet_f(I,cSetting_specular_intensity,0.8F);
  SettingSet_i(I,cSetting_overlay_lines,5);
  SettingSet_f(I,cSetting_ray_transparency_spec_cut,0.9F);
  SettingSet_b(I,cSetting_internal_prompt,1);
  SettingSet_b(I,cSetting_normalize_grd_maps,0);

  SettingSet_b(I,cSetting_ray_blend_colors,0);
  SettingSet_f(I,cSetting_ray_blend_red,0.17F);
  SettingSet_f(I,cSetting_ray_blend_green,0.25F);
  SettingSet_f(I,cSetting_ray_blend_blue,0.14F);
  SettingSet_f(I,cSetting_png_screen_gamma,2.4F);
  SettingSet_f(I,cSetting_png_file_gamma,1.0F);
  SettingSet_b(I,cSetting_editor_label_fragments,0);
  SettingSet_i(I,cSetting_internal_gui_control_size, 18);
  SettingSet_b(I,cSetting_auto_dss,1);
  SettingSet_i(I,cSetting_transparency_picking_mode,2); /* auto */
  SettingSet_b(I,cSetting_virtual_trackball,1);
  SettingSet_i(I,cSetting_transparency_picking_mode,2); /* auto */
  SettingSet_i(I,cSetting_pdb_reformat_names_mode,0); /*
                                                        0 = no reformatting, 
                                                        1 = pdb compliant,
                                                        2 = amber compliant,
                                                        3 = pdb I/O, but iupac inside
                                                      */
  SettingSet_b(I,cSetting_ray_pixel_scale_to_window,1);
  SettingSet_i(I,cSetting_label_font_id,0);
  SettingSet_b(I,cSetting_pdb_conect_all,0);
  SettingSet_s(I,cSetting_button_mode_name,"");
  SettingSet_i(I,cSetting_surface_type,0);
  SettingSet_b(I,cSetting_dot_normals,1);
  SettingSet_b(I,cSetting_session_migration,1);
  SettingSet_b(I,cSetting_mesh_normals,1);
  SettingSet_b(I,cSetting_mesh_type,0);

  SettingSet_b(I,cSetting_dot_lighting,1);
  SettingSet_b(I,cSetting_mesh_lighting,0);
  SettingSet_b(I,cSetting_surface_solvent,0);
  SettingSet_i(I,cSetting_triangle_max_passes,5);
  SettingSet_f(I,cSetting_ray_interior_reflect,0.4F);
  SettingSet_i(I,cSetting_internal_gui_mode,0);
  SettingSet_s(I,cSetting_surface_carve_selection,"");
  SettingSet_i(I,cSetting_surface_carve_state,0);
  SettingSet_f(I,cSetting_surface_carve_cutoff,0.0F);
  SettingSet_s(I,cSetting_surface_clear_selection,"");
  SettingSet_i(I,cSetting_surface_clear_state,0);
  SettingSet_f(I,cSetting_surface_clear_cutoff,0.0F);
  SettingSet_f(I,cSetting_surface_trim_cutoff,0.2F);
  SettingSet_f(I,cSetting_surface_trim_factor,2.0F);
  SettingSet_i(I,cSetting_ray_max_passes,25);
  SettingSet_b(I,cSetting_active_selections,true);
  SettingSet_f(I,cSetting_ray_transparency_contrast,1.0F);
  SettingSet_b(I,cSetting_seq_view,0);
  SettingSet_i(I,cSetting_mouse_selection_mode,1);
  SettingSet_i(I,cSetting_seq_view_label_spacing,5);
  SettingSet_i(I,cSetting_seq_view_label_start,1);
  SettingSet_i(I,cSetting_seq_view_format,0);
  SettingSet_i(I,cSetting_seq_view_location,0);
  SettingSet_b(I,cSetting_seq_view_overlay,0);
  SettingSet_b(I,cSetting_auto_classify_atoms,1);
  SettingSet_i(I,cSetting_cartoon_nucleic_acid_mode,0);
  SettingSet_color(I,cSetting_seq_view_color, "-1");
  SettingSet_i(I,cSetting_seq_view_label_mode,2);
  SettingSet_i(I,cSetting_surface_ramp_above_mode,0);
  SettingSet_b(I,cSetting_stereo,0);
  SettingSet_i(I,cSetting_wizard_prompt_mode,1);
  SettingSet_f(I,cSetting_coulomb_cutoff,10.0F);
  SettingSet_b(I,cSetting_slice_track_camera, 0);
  SettingSet_f(I,cSetting_slice_height_scale, 1.0F);
  SettingSet_b(I,cSetting_slice_height_map, 0);
  SettingSet_f(I,cSetting_slice_grid, 0.3F);
  SettingSet_b(I,cSetting_slice_dynamic_grid, 0);
  SettingSet_f(I,cSetting_slice_dynamic_grid_resolution, 3.0F);
  SettingSet_b(I,cSetting_pdb_insure_orthogonal, 1);
  SettingSet_f(I,cSetting_ray_direct_shade,0.0F);
  SettingSet_color(I,cSetting_stick_color,"-1");
  SettingSet_f(I,cSetting_cartoon_putty_radius, 0.40F);
  SettingSet_f(I,cSetting_cartoon_putty_quality, 11.0F);
  SettingSet_f(I,cSetting_cartoon_putty_scale_min, 0.6F);
  SettingSet_f(I,cSetting_cartoon_putty_scale_max, 4.0F);
  SettingSet_f(I,cSetting_cartoon_putty_scale_power, 1.5F);
  SettingSet_f(I,cSetting_cartoon_putty_range, 2.0F);
  SettingSet_b(I,cSetting_cartoon_side_chain_helper, 0);
  SettingSet_b(I,cSetting_surface_optimize_subsets, 1);
  SettingSet_i(I,cSetting_multiplex, -1);
  SettingSet_b(I,cSetting_texture_fonts, 0);
}


