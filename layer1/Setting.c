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

CSetting Setting;

static void *SettingPtr(CSetting *I,int index,unsigned int size);


/*========================================================================*/
PyObject *SettingGetUpdateList(CSetting *I)
{
  int a;
  int n;
  int unblock;
  PyObject *result;

  if(!I) I=&Setting; /* fall back on global settings */

  unblock = PAutoBlock();
  n=VLAGetSize(I->info);
  result=PyList_New(0);
  for(a=0;a<n;a++) {
    if(I->info[a].changed) {
      I->info[a].changed=false;
      PyList_Append(result,PyInt_FromLong(a));
    }
  }
  PAutoUnblock(unblock);

  return(result);
}
/*========================================================================*/
void SettingCheckHandle(CSetting **handle)
{
  if(!*handle)
    *handle=SettingNew();
}
/*========================================================================*/
void SettingGetTextValue(CSetting *set1,CSetting *set2,int index,char *buffer) 
/* not range checked */
{
  int type;
  float *ptr;
  type = SettingGetType(index);
  switch(type) {
  case cSetting_boolean:
    if(SettingGet_b(set1,set2,index))
      sprintf(buffer,"on.");
    else
      sprintf(buffer,"off.");      
    break;
  case cSetting_int:
    sprintf(buffer,"%d",SettingGet_i(set1,set2,index));
    break;
  case cSetting_float:
    sprintf(buffer,"%1.5f",SettingGet_f(set1,set2,index));
    break;
  case cSetting_float3:
    ptr = SettingGet_fv(set1,set2,index);
    sprintf(buffer,"[ %1.5f, %1.5f, %1.5f ]",
            ptr[0],ptr[1],ptr[2]);
    break;
  }
}

/*========================================================================*/
int SettingSetTuple(CSetting *I,int index,PyObject *tuple) 
     /* must have interpret locked to make this call */
{
  PyObject *value;
  int type;
  
  if(!I) I=&Setting; /* fall back on global settings */

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
   SettingSet_f(I,index,PyFloat_AsDouble(PyTuple_GetItem(value,0)));
    break;
  case cSetting_float3:
        SettingSet_3f(I,index,
                  PyFloat_AsDouble(PyTuple_GetItem(value,0)),
                  PyFloat_AsDouble(PyTuple_GetItem(value,1)),
                  PyFloat_AsDouble(PyTuple_GetItem(value,2)));
    break;
  }
  return 1;
}
/*========================================================================*/
PyObject *SettingGetTuple(CSetting *set1,CSetting *set2,int index)
{
  PyObject *result;
  float *ptr;
  int type = SettingGetType(index);
  switch(type) {
  case cSetting_boolean:
    result = Py_BuildValue("(i(i))",type,
                           SettingGet_b(set1,set2,index));
    break;
  case cSetting_int:
    result = Py_BuildValue("(i(i))",type,
                           SettingGet_i(set1,set2,index));
    break;
  case cSetting_float:
    result = Py_BuildValue("(i(f))",type,
                           SettingGet_f(set1,set2,index));
    break;
  case cSetting_float3:
    ptr =  SettingGet_fv(set1,set2,index);
    result = Py_BuildValue("(i(fff))",type,
                           ptr[0],ptr[1],ptr[2]);
    break;
  default:
    Py_INCREF(Py_None);
    result = Py_None;
    break;
  }
  return result;
}
/*========================================================================*/
CSetting *SettingNew(void)
{
  OOAlloc(CSetting);
  SettingInit(I);
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
void SettingInit(CSetting *I)
{
  I->size=sizeof(int); /* insures offset is never zero, except when undef */
  I->data=VLAlloc(char,10);
  I->info=VLAMalloc(cSetting_INIT,sizeof(SettingRec),5,1); /* auto-zero */
}
/*========================================================================*/
void SettingClear(CSetting *I,int index)
{
  I->info[index].defined = false; 
}
/*========================================================================*/
static void *SettingPtr(CSetting *I,int index,unsigned int size)
{
  SettingRec *sr = I->info+index;
  if(!sr->offset) {
    sr->offset=I->size;
    I->size+=size;
    VLACheck(I->data,char,I->size);
  }
  sr->defined = true;
  sr->changed = true;
  return(I->data+sr->offset);
}
/*========================================================================*/
int SettingGetType(int index)
{
  CSetting *I=&Setting;
  return(I->info[index].type);
}
/*========================================================================*/
void SettingSet_b(CSetting *I,int index, int value)
{
  if(Setting.info[index].type&&(Setting.info[index].type!=cSetting_boolean)) {
    PRINTFB(FB_Setting,FB_Errors)
      "Setting-Error: type mismatch (boolean)\n"
      ENDFB
  } else {
    VLACheck(I->info,SettingRec,index);
    *((int*)SettingPtr(I,index,sizeof(int))) = value;
    I->info[index].type = cSetting_boolean;
  }
}
/*========================================================================*/
void SettingSet_i(CSetting *I,int index, int value)
{
  if(Setting.info[index].type&&(Setting.info[index].type!=cSetting_int)) {
    PRINTFB(FB_Setting,FB_Errors)
      "Setting-Error: type mismatch (int)\n"
      ENDFB
  } else {
    VLACheck(I->info,SettingRec,index);
    *((int*)SettingPtr(I,index,sizeof(int))) = value;
    I->info[index].type = cSetting_int;
  }
}
/*========================================================================*/
void SettingSet_f(CSetting *I,int index, float value)
{
  if(Setting.info[index].type&&(Setting.info[index].type!=cSetting_float)) {
    PRINTFB(FB_Setting,FB_Errors)
      "Setting-Error: type mismatch (float)\n"
      ENDFB
  } else {
    VLACheck(I->info,SettingRec,index);
    *((float*)SettingPtr(I,index,sizeof(float)))=value;
    I->info[index].type = cSetting_float;
  }
}
/*========================================================================*/
void SettingSet_3f(CSetting *I,int index, float value1,float value2,float value3)
{
  float *ptr;
  if(Setting.info[index].type&&(Setting.info[index].type!=cSetting_float3)) {
    PRINTFB(FB_Setting,FB_Errors)
      "Setting-Error: type mismatch (float3)\n"
      ENDFB
  } else {
    VLACheck(I->info,SettingRec,index);
    ptr = (float*)SettingPtr(I,index,sizeof(float)*3);
    ptr[0]=value1;
    ptr[1]=value2;
    ptr[2]=value3;
    I->info[index].type = cSetting_float3;
  }
}
/*========================================================================*/
void SettingSet_3fv(CSetting *I,int index, float *vector)
{
  float *ptr;
  VLACheck(I->info,SettingRec,index);
  ptr = (float*)SettingPtr(I,index,sizeof(float)*3);
  copy3f(vector,ptr);
  I->info[index].type = cSetting_float3;
}
/*========================================================================*/
int   SettingGetGlobal_b(int index) 
{
  CSetting *I=&Setting;
  return(*((int*)(I->data+I->info[index].offset)));
}
/*========================================================================*/
int   SettingGetGlobal_i(int index) 
{
  CSetting *I=&Setting;
  return(*((int*)(I->data+I->info[index].offset)));
}
/*========================================================================*/
float SettingGetGlobal_f(int index)
{
  CSetting *I=&Setting;
  return(*((float*)(I->data+I->info[index].offset)));
}
/*========================================================================*/
void  SettingGetGlobal_3f(int index,float *value)
{
  CSetting *I=&Setting;
  float *ptr;
  ptr = (float*)(I->data+I->info[index].offset);
  copy3f(ptr,value);
}
/*========================================================================*/
float *SettingGetGlobal_fv(int index)
{
  CSetting *I=&Setting;
  return (float*)(I->data+I->info[index].offset);
}
/*========================================================================*/
int   SettingGet_b  (CSetting *set1,CSetting *set2,int index)
{
  if(set1) {
    if(set1->info[index].defined) {
      return(*((int*)(set1->data+set1->info[index].offset)));
    }
  }
  if(set2) { 
    if(set2->info[index].defined) {
      return(*((int*)(set2->data+set2->info[index].offset)));      
    }
  }
  return(SettingGetGlobal_i(index));
}
/*========================================================================*/
int   SettingGet_i  (CSetting *set1,CSetting *set2,int index)
{
  if(set1) {
    if(set1->info[index].defined) {
      return(*((int*)(set1->data+set1->info[index].offset)));
    }
  }
  if(set2) { 
    if(set2->info[index].defined) {
      return(*((int*)(set2->data+set2->info[index].offset)));      
    }
  }
  return(SettingGetGlobal_i(index));
}
/*========================================================================*/
float SettingGet_f  (CSetting *set1,CSetting *set2,int index)
{
  if(set1) {
    if(set1->info[index].defined) {
      return(*((float*)(set1->data+set1->info[index].offset)));
    }
  }
  if(set2) {
    if(set2->info[index].defined) {
      return(*((float*)(set2->data+set2->info[index].offset)));      
    }
  }
  return(SettingGetGlobal_f(index));
}
/*========================================================================*/
void  SettingGet_3f(CSetting *set1,CSetting *set2,int index,float *value)
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
  SettingGetGlobal_3f(index,value);
}
/*========================================================================*/
float *SettingGet_fv(CSetting *set1,CSetting *set2,int index)
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
  return(SettingGetGlobal_fv(index));
}

/*========================================================================*/
/*========================================================================*/
int SettingGetIndex(char *name) /* can be called from any thread state */
{
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
}
/*========================================================================*/
int SettingGetName(int index,SettingName name) /* can be called from any thread state */
{
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
}

/*========================================================================*/
void SettingGenerateSideEffects(int index,char *sele,int state)
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
  case cSetting_light:
	 SceneDirty();
	 break;
  case cSetting_valence:
    ExecutiveInvalidateRep(inv_sele,cRepLine,cRepInvRep);
    SceneChanged();
    break;
  case cSetting_dash_length:
  case cSetting_dash_gap:
  case cSetting_dash_radius:
  case cSetting_dash_width:
    ExecutiveRebuildAllObjectDist();
    SceneChanged();
    break;
  case cSetting_button_mode:
    SceneDirty();
    break;
  case cSetting_stick_radius:
  case cSetting_stick_quality:
  case cSetting_stick_overlap:
    ExecutiveInvalidateRep(inv_sele,cRepCyl,cRepInvRep);
    SceneChanged();
    break;
  case cSetting_label_color:
    ExecutiveInvalidateRep(inv_sele,cRepLabel,cRepInvRep);
    SceneChanged();
    break;
  case cSetting_all_states:
    SceneChanged();
    break;
  case cSetting_dot_density:
	 break;
  case cSetting_sel_counter:
	 break;
  case cSetting_line_width: /* auto-disable smooth lines if line width > 1 */
    SettingSet(cSetting_line_smooth,0);
    SceneChanged();
    break;
  case cSetting_mesh_width: /* auto-disable smooth lines if line width > 1 */
    ExecutiveInvalidateRep(inv_sele,cRepMesh,cRepInvColor);
    SettingSet(cSetting_line_smooth,0);
    SceneChanged();
    break;
  case cSetting_nonbonded_size:
    ExecutiveInvalidateRep(inv_sele,cRepNonbonded,cRepInvRep);
    ExecutiveInvalidateRep(inv_sele,cRepNonbondedSphere,cRepInvRep);
    SceneChanged();
    break;
  case cSetting_mesh_radius:
    ExecutiveInvalidateRep(inv_sele,cRepMesh,cRepInvColor);
    SceneChanged();
    break;
  case cSetting_surface_quality:
    ExecutiveInvalidateRep(inv_sele,cRepSurface,cRepInvRep);
    SceneChanged();
    break;
  case cSetting_ribbon_power:
  case cSetting_ribbon_power_b:
  case cSetting_ribbon_sampling:
  case cSetting_ribbon_radius:
  case cSetting_ribbon_width:
    ExecutiveInvalidateRep(inv_sele,cRepRibbon,cRepInvRep);
    SceneChanged();
    break;
  case cSetting_cartoon_sampling:
  case cSetting_cartoon_loop_quality:
  case cSetting_cartoon_loop_radius:
  case cSetting_cartoon_tube_quality:
  case cSetting_cartoon_tube_radius:
  case cSetting_cartoon_power:
  case cSetting_cartoon_power_b:
  case cSetting_cartoon_rect_length:
  case cSetting_cartoon_rect_width:
  case cSetting_cartoon_oval_length:
  case cSetting_cartoon_oval_width:
  case cSetting_cartoon_oval_quality:
  case cSetting_cartoon_debug:
    ExecutiveInvalidateRep(inv_sele,cRepCartoon,cRepInvRep);
    SceneChanged();
    break;
  case cSetting_dot_width:
  case cSetting_line_smooth:
  case cSetting_ortho:
  case cSetting_ambient:
  case cSetting_bg_rgb:
  case cSetting_depth_cue:
  case cSetting_specular:
	 SceneDirty();
  case cSetting_overlay:
  case cSetting_text:
    OrthoDirty();
  case cSetting_internal_gui_width:
  case cSetting_internal_gui:
    sprintf(command,"viewport");
    OrthoCommandIn(command);
  default:
	 break;
  }
}
/*========================================================================*/
void SettingSetfv(int index,float *v)
{
  CSetting *I=&Setting;
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
	 SceneDirty();
	 break;
  case cSetting_valence:
    ExecutiveInvalidateRep("all",cRepLine,cRepInvRep);
    SettingSet_f(I,index,v[0]);
    /*   I->Setting[index].Value[0]=v[0];    */
    SceneChanged();
    break;
  case cSetting_dash_length:
  case cSetting_dash_gap:
    ExecutiveInvalidateRep("all",cRepDash,cRepInvRep);
    SettingSet_f(I,index,v[0]);
    /*    I->Setting[index].Value[0]=v[0];    */
    SceneChanged();
    break;
  case cSetting_button_mode:
    SettingSet_f(I,index,v[0]);
    /* I->Setting[index].Value[0]=v[0]; */
    SceneDirty();
    break;
  case cSetting_stick_radius:
  case cSetting_stick_quality:
  case cSetting_stick_overlap:
    ExecutiveInvalidateRep("all",cRepCyl,cRepInvRep);
    SettingSet_f(I,index,v[0]);
    /*I->Setting[index].Value[0]=v[0];   */
    SceneChanged();
    break;
  case cSetting_label_color:
    ExecutiveInvalidateRep("all",cRepLabel,cRepInvRep);
    SettingSet_f(I,index,v[0]);
    /* I->Setting[index].Value[0]=v[0]; */
    SceneChanged();
    break;
  case cSetting_all_states:
    SettingSet_f(I,index,v[0]);
    /* I->Setting[index].Value[0]=v[0];  */
    SceneChanged();
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
  case cSetting_ambient:
	 SceneDirty();
    break;
  case cSetting_overlay:
  case cSetting_text:
    OrthoDirty();
  default:
    SettingSet_f(I,index,v[0]);
    /*I->Setting[index].Value[0]=v[0];*/
	 break;
  }
}
/*========================================================================*/
void SettingSet(int index,float v)
{
  SettingSetfv(index,&v);
}
/*========================================================================*/
void SettingSetNamed(char *name,char *value)
{
  int index = SettingGetIndex(name);
  float v,vv[3];
  SettingName realName;
  char buffer[1024] = "";
  if(index>=0) {
    SettingGetName(index,realName);
	 switch(index) {
	 case cSetting_dot_mode:
		if(strcmp(value,"molecular")==0) {
		  v=0.0;
		  SettingSetfv(index,&v);
		  sprintf(buffer," Setting: %s set to %s\n",realName,value);
		} else if(strcmp(value,"solvent_accessible")==0) {
		  v=1.0;
		  SettingSetfv(index,&v);
		  sprintf(buffer," Setting: %s set to %s\n",realName,value);
		} else if(sscanf(value,"%f",&v)==1) {
		  SettingSetfv(index,&v);
		  sprintf(buffer," Setting: %s set to %s\n",realName,value);
      }
		break;
	 case cSetting_bg_rgb:
	 case cSetting_light:
		if(sscanf(value,"%f%f%f",vv,vv+1,vv+2)==3) {
		  SettingSetfv(index,vv);
		  sprintf(buffer," Setting: %s set to %5.3f %8.3f %8.3f\n",realName,
					 *vv,*(vv+1),*(vv+2));
		}
		break;
	 case cSetting_dot_density:
		sscanf(value,"%f",&v);
		SettingSetfv(index,&v);
		sprintf(buffer," Setting: %s set to %d\n",realName,(int)v);
		break;
	 case cSetting_text:
	 case cSetting_overlay:
	 case cSetting_sel_counter:
    case cSetting_dist_counter:
		sscanf(value,"%f",&v);
		SettingSetfv(index,&v);
		break;
    case cSetting_line_width: /* auto-disable smooth lines if line width > 1 */
    case cSetting_mesh_width:
      sscanf(value,"%f",&v);
      if(v!=1.0)
        SettingSet(cSetting_line_smooth,0);
		SettingSetfv(index,&v);
		sprintf(buffer," Setting: %s set to %5.3f\n",realName,v);
      SceneDirty();
      break;
	 default:
		sscanf(value,"%f",&v);
		SettingSetfv(index,&v);
		sprintf(buffer," Setting: %s set to %5.3f\n",realName,v);
		break;
	 }
  } else {
    PRINTFB(FB_Setting,FB_Warnings)
      " Error: Non-Existent Settin\n"
      ENDFB;
  }
  if(buffer[0]) {
    PRINTFB(FB_Setting,FB_Actions)
      "%s",buffer
      ENDFB;
  }
}
/*========================================================================*/
float SettingGetNamed(char *name)
{
  return(SettingGet(SettingGetIndex(name)));
}
/*========================================================================*/
float SettingGet(int index)
{
  return(SettingGetGlobal_f(index));
}
/*========================================================================*/
float *SettingGetfv(int index)
{
  return(SettingGetGlobal_fv(index));
}
/*========================================================================*/
void SettingFreeGlobal(void)
{
  CSetting *I=&Setting;
  SettingPurge(I);
}
/*========================================================================*/
void SettingInitGlobal(void)
{
  CSetting *I=&Setting;

  SettingInit(I);

  SettingSet_f(I,cSetting_bonding_vdw_cutoff, 0.2);

  SettingSet_f(I,cSetting_min_mesh_spacing, 0.6);

  SettingSet_f(I,cSetting_dot_density, 2);

  SettingSet_f(I,cSetting_dot_mode, 0);

  SettingSet_f(I,cSetting_solvent_radius, 1.4);

  SettingSet_f(I,cSetting_sel_counter, 0.0);

  SettingSet_3f(I,cSetting_bg_rgb, 0.0, 0.0, 0.0);

  SettingSet_f(I,cSetting_ambient, 0.30);

  SettingSet_f(I,cSetting_direct, 0.35);

  SettingSet_f(I,cSetting_reflect, 1.2);

  SettingSet_3f(I,cSetting_light, -0.4, -0.4, -1.0);

  SettingSet_f(I,cSetting_antialias, 0.0);

  SettingSet_f(I,cSetting_cavity_cull, 10.0);

  SettingSet_f(I,cSetting_ambient_scale, 0.4);

  SettingSet_f(I,cSetting_single_image, 0.0);

  SettingSet_f(I,cSetting_movie_delay, 30);

  SettingSet_f(I,cSetting_ribbon_power, 2);

  SettingSet_f(I,cSetting_ribbon_power_b, 0.5);

  SettingSet_f(I,cSetting_ribbon_sampling, 16);

  SettingSet_f(I,cSetting_ribbon_radius, 0.4);

  SettingSet_f(I,cSetting_stick_radius, 0.3);

  SettingSet_f(I,cSetting_hash_max, 80);

  SettingSet_f(I,cSetting_ortho, 0);

  SettingSet_f(I,cSetting_power, 3.0);

  SettingSet_f(I,cSetting_spec_reflect, 0.4);

  SettingSet_f(I,cSetting_spec_power, 40);

  SettingSet_f(I,cSetting_sweep_angle, 15.0);

  SettingSet_f(I,cSetting_sweep_speed, 0.3);

  SettingSet_f(I,cSetting_dot_hydrogens, 1.0);

  SettingSet_f(I,cSetting_dot_size, 0.06);

  SettingSet_f(I,cSetting_ray_trace_frames, 0.0);

  SettingSet_f(I,cSetting_cache_frames, 0.0);

  SettingSet_f(I,cSetting_trim_dots, 1.0);

  SettingSet_f(I,cSetting_cull_spheres, 1.0);

  SettingSet_f(I,cSetting_test1, 0.0);

  SettingSet_f(I,cSetting_test2, 1.0);

  SettingSet_f(I,cSetting_surface_best, 0.2);

  SettingSet_f(I,cSetting_surface_normal, 0.5);

  SettingSet_f(I,cSetting_surface_quality, 0.0);

  SettingSet_f(I,cSetting_surface_proximity, 1.5);

  SettingSet_f(I,cSetting_stereo_angle, 2.1);

  SettingSet_f(I,cSetting_stereo_shift, 2.0);

  SettingSet_f(I,cSetting_line_smooth, 0.0);

  SettingSet_f(I,cSetting_line_width, 1.0);

  SettingSet_f(I,cSetting_half_bonds, 0.0);

  SettingSet_f(I,cSetting_stick_quality, 8.0);

  SettingSet_f(I,cSetting_stick_overlap, 0.2);

  SettingSet_f(I,cSetting_stick_nub, 0.7);

  SettingSet_f(I,cSetting_all_states, 0.0);

  SettingSet_f(I,cSetting_pickable, 1.0);

  SettingSet_f(I,cSetting_autoshow_lines, 1.0);

  SettingSet_f(I,cSetting_fast_idle, 20000);

  SettingSet_f(I,cSetting_no_idle, 5000);

  SettingSet_f(I,cSetting_slow_idle, 200000);

  SettingSet_f(I,cSetting_idle_delay, 1.5);

  SettingSet_f(I,cSetting_rock_delay, 30);

  SettingSet_f(I,cSetting_dist_counter, 0.0);

  SettingSet_f(I,cSetting_dash_length, 0.15);

  SettingSet_f(I,cSetting_dash_gap, 0.35);

  SettingSet_f(I,cSetting_auto_zoom, 1.0);

  SettingSet_f(I,cSetting_overlay, 0.0);

  SettingSet_f(I,cSetting_text, 0.0);

  SettingSet_f(I,cSetting_button_mode, 0.0);

  SettingSet_f(I,cSetting_valence, 0.0);

  SettingSet_f(I,cSetting_nonbonded_size, 0.25);

  SettingSet_f(I,cSetting_label_color, -1.0);

  SettingSet_f(I,cSetting_ray_trace_fog, 1.0);

  SettingSet_f(I,cSetting_spheroid_scale, 1.0);

  SettingSet_f(I,cSetting_ray_trace_fog_start, 0.30);

  SettingSet_f(I,cSetting_spheroid_smooth, 1.1);

  SettingSet_f(I,cSetting_spheroid_fill, 1.30);

  SettingSet_f(I,cSetting_autoshow_nonbonded, 1.0);

  SettingSet_f(I,cSetting_mesh_radius, 0.05);
 
#ifdef WIN32
  SettingSet_f(I,cSetting_cache_display, 0.0);
#else
  SettingSet_f(I,cSetting_cache_display, 1.0);
#endif

  SettingSet_f(I,cSetting_normal_workaround, 0);

  SettingSet_f(I,cSetting_backface_cull, 1);

  SettingSet_f(I,cSetting_gamma, 1.2);

  SettingSet_f(I,cSetting_dot_width, 1.0);

  SettingSet_f(I,cSetting_autoshow_selections, 1.0);

  SettingSet_f(I,cSetting_autohide_selections, 1.0);

  SettingSet_f(I,cSetting_selection_width, 5.0);

  SettingSet_f(I,cSetting_selection_overlay, 1.0);

  SettingSet_f(I,cSetting_static_singletons, 1.0);

  SettingSet_f(I,cSetting_max_triangles, 1000000.0);

  SettingSet_f(I,cSetting_depth_cue, 1.0);

  SettingSet_f(I,cSetting_specular, 0.8);

  SettingSet_f(I,cSetting_shininess, 40.0);

  SettingSet_f(I,cSetting_sphere_quality, 1.0);

  SettingSet_f(I,cSetting_fog, 0.75);

  SettingSet_f(I,cSetting_isomesh_auto_state, 1.0);

  SettingSet_f(I,cSetting_mesh_width, 1.0);

  SettingSet_f(I,cSetting_cartoon_sampling, 7);

  SettingSet_f(I,cSetting_cartoon_loop_radius, 0.2);

  SettingSet_f(I,cSetting_cartoon_loop_quality, 6);

  SettingSet_f(I,cSetting_cartoon_power, 1.8);

  SettingSet_f(I,cSetting_cartoon_power_b, 0.45);

  SettingSet_f(I,cSetting_cartoon_rect_length, 1.40);

  SettingSet_f(I,cSetting_cartoon_rect_width, 0.4);

  SettingSet_f(I,cSetting_internal_gui_width, 160);

  SettingSet_f(I,cSetting_internal_gui, 1.0);

  SettingSet_f(I,cSetting_cartoon_oval_length, 1.50);

  SettingSet_f(I,cSetting_cartoon_oval_width, 0.25);

  SettingSet_f(I,cSetting_cartoon_oval_quality, 10);

  SettingSet_f(I,cSetting_cartoon_tube_radius, 0.9);

  SettingSet_f(I,cSetting_cartoon_tube_quality, 9);

  SettingSet_f(I,cSetting_cartoon_debug, 0);

  SettingSet_f(I,cSetting_ribbon_width, 1.0);

  SettingSet_f(I,cSetting_dash_width, 1.0);

  SettingSet_f(I,cSetting_dash_radius, 0.14);

  
}

