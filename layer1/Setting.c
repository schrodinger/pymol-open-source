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

CSetting Setting;

/*========================================================================*/
int SettingGetIndex(char *name)
{
  CSetting *I=&Setting;
  int index=-1; 
  int a;
  for(a=0;a<I->NSetting;a++)
	 {
		if(strcmp(name,I->Setting[a].Name)==0) 
		  {
			 index=a;
			 break;
		  }
	 }
  return(index);
}
/*========================================================================*/
void SettingSetfv(int index,float *v)
{
  CSetting *I=&Setting;
  switch(index) {
  case cSetting_dot_surface:
	 I->Setting[index].Value[0]=v[0];
	 break;
  case cSetting_bg_rgb:
  case cSetting_light:
	 I->Setting[index].Value[0]=v[0];
	 I->Setting[index].Value[1]=v[1];
	 I->Setting[index].Value[2]=v[2];
	 SceneDirty();
	 break;
  case cSetting_valence:
    ExecutiveInvalidateRep("all",cRepLine,cRepInvRep);
    I->Setting[index].Value[0]=v[0];
    SceneChanged();
    break;
  case cSetting_dash_length:
  case cSetting_dash_gap:
    ExecutiveInvalidateRep("all",cRepDash,cRepInvRep);
    I->Setting[index].Value[0]=v[0];
    SceneChanged();
    break;
  case cSetting_button_mode:
    EditorSetActiveObject(NULL,0);
	 I->Setting[index].Value[0]=v[0];
    SceneDirty();
    break;
  case cSetting_stick_radius:
  case cSetting_stick_quality:
  case cSetting_stick_overlap:
    ExecutiveInvalidateRep("all",cRepCyl,cRepInvRep);
    I->Setting[index].Value[0]=v[0];
    SceneChanged();
    break;
  case cSetting_label_color:
    ExecutiveInvalidateRep("all",cRepLabel,cRepInvRep);
    I->Setting[index].Value[0]=v[0];
    SceneChanged();
    break;
  case cSetting_all_states:
	 I->Setting[index].Value[0]=v[0];
    SceneChanged();
    break;
  case cSetting_dot_density:
	 I->Setting[index].Value[0]=v[0];
	 break;
  case cSetting_sel_counter:
	 I->Setting[index].Value[0]=v[0];
	 break;
  case cSetting_ortho:
  case cSetting_ambient:
	 SceneDirty();
  case cSetting_overlay:
  case cSetting_text:
    OrthoDirty();
  default:
	 I->Setting[index].Value[0]=v[0];
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
  CSetting *I=&Setting;
  int index = SettingGetIndex(name);
  float v,vv[3];
  char buffer[1024] = "";
  if(index>=0) {
	 switch(index) {
	 case cSetting_dot_surface:
		if(strcmp(value,"molecular")==0) {
		  v=0.0;
		  SettingSetfv(index,&v);
		  sprintf(buffer," Setting: %s set to %s\n",I->Setting[index].Name,value);
		} else if(strcmp(value,"solvent_accessible")==0) {
		  v=1.0;
		  SettingSetfv(index,&v);
		  sprintf(buffer," Setting: %s set to %s\n",I->Setting[index].Name,value);
		}
		break;
	 case cSetting_bg_rgb:
	 case cSetting_light:
		if(sscanf(value,"%f%f%f",vv,vv+1,vv+2)==3) {
		  SettingSetfv(index,vv);
		  sprintf(buffer," Setting: %s set to %5.3f %8.3f %8.3f\n",I->Setting[index].Name,
					 *vv,*(vv+1),*(vv+2));
		}
		break;
	 case cSetting_dot_density:
		sscanf(value,"%f",&v);
		SettingSetfv(index,&v);
		sprintf(buffer," Setting: %s set to %d\n",I->Setting[index].Name,(int)v);
		break;
	 case cSetting_text:
	 case cSetting_overlay:
	 case cSetting_sel_counter:
    case cSetting_dist_counter:
		sscanf(value,"%f",&v);
		SettingSetfv(index,&v);
		break;
	 default:
		sscanf(value,"%f",&v);
		SettingSetfv(index,&v);
		sprintf(buffer," Setting: %s set to %5.3f\n",I->Setting[index].Name,v);
		break;
	 }
  } else {
	 OrthoAddOutput(" Error: Non-Existent Setting");
	 OrthoNewLine(NULL);
	 /*
	 VLACheck(I->Setting,SettingRec,I->NSetting);
	 I->Setting[I->NSetting].Value[0] = value;
	 strcpy(I->Setting[I->NSetting].Name,name);
	 I->NSetting++;*/
  }
  if(buffer[0]) {
	 OrthoAddOutput(buffer);
	 OrthoNewLine(NULL);
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
  CSetting *I=&Setting;
  if((index>=0)&&(index<I->NSetting))
	 return(I->Setting[index].Value[0]);
  else
	 return(0.0);
}

float *SettingGetfv(int index)
{
  CSetting *I=&Setting;
  if((index>=0)&&(index<I->NSetting))
	 return(I->Setting[index].Value);
  else
	 return(NULL);
}
/*========================================================================*/
void SettingFree(void)
{
  CSetting *I=&Setting;
  VLAFreeP(I->Setting);
}
/*========================================================================*/
void SettingInit(void)
{
  CSetting *I=&Setting;

  I->Setting=VLAlloc(SettingRec,100);
  I->NSetting=0;

  I->NSetting++;
  I->Setting[cSetting_bonding_vdw_cutoff].Value[0] = 0.2;
  strcpy(I->Setting[cSetting_bonding_vdw_cutoff].Name,
			"bonding_vdw_cutoff");

  I->NSetting++;
  I->Setting[cSetting_min_mesh_spacing].Value[0] = 0.6;
  strcpy(I->Setting[cSetting_min_mesh_spacing].Name,
			"min_mesh_spacing");

  I->NSetting++;
  I->Setting[cSetting_dot_density].Value[0] = 2;
  strcpy(I->Setting[cSetting_dot_density].Name,
			"dot_density");

  I->NSetting++;
  I->Setting[cSetting_dot_surface].Value[0] = 0;
  strcpy(I->Setting[cSetting_dot_surface].Name,
			"dot_surface");

  I->NSetting++;
  I->Setting[cSetting_solvent_radius].Value[0] = 1.4;
  strcpy(I->Setting[cSetting_solvent_radius].Name,
			"solvent_radius");

  I->NSetting++;
  I->Setting[cSetting_sel_counter].Value[0] = 0.0;
  strcpy(I->Setting[cSetting_sel_counter].Name,
			"sel_counter");

  I->NSetting++;
  I->Setting[cSetting_bg_rgb].Value[0] = 0.0;
  I->Setting[cSetting_bg_rgb].Value[1] = 0.0;
  I->Setting[cSetting_bg_rgb].Value[2] = 0.0;
  strcpy(I->Setting[cSetting_bg_rgb].Name,
			"bg_rgb");

  I->NSetting++;
  I->Setting[cSetting_ambient].Value[0] = 0.30;
  strcpy(I->Setting[cSetting_ambient].Name,
			"ambient");

  I->NSetting++;
  I->Setting[cSetting_direct].Value[0] = 0.35;
  strcpy(I->Setting[cSetting_direct].Name,
			"direct");

  I->NSetting++;
  I->Setting[cSetting_reflect].Value[0] = 1.2;
  strcpy(I->Setting[cSetting_reflect].Name,
			"reflect");

  I->NSetting++;
  I->Setting[cSetting_light].Value[0] = -0.4;
  I->Setting[cSetting_light].Value[1] = -0.4;
  I->Setting[cSetting_light].Value[2] = -1.0;
  strcpy(I->Setting[cSetting_light].Name,
			"light");

  I->NSetting++;
  I->Setting[cSetting_antialias].Value[0] = 0.0;
  strcpy(I->Setting[cSetting_antialias].Name,
			"antialias");

  I->NSetting++;
  I->Setting[cSetting_cavity_cull].Value[0] = 10.0;
  strcpy(I->Setting[cSetting_cavity_cull].Name,
			"cavity_cull");

  I->NSetting++;
  I->Setting[cSetting_ambient_scale].Value[0] = 0.4;
  strcpy(I->Setting[cSetting_ambient_scale].Name,
			"ambient_scale");

  I->NSetting++;
  I->Setting[cSetting_single_image].Value[0] = 0.0;
  strcpy(I->Setting[cSetting_single_image].Name,
			"single_image");

  I->NSetting++;
  I->Setting[cSetting_movie_delay].Value[0] = 30;
  strcpy(I->Setting[cSetting_movie_delay].Name,
			"movie_delay");

  I->NSetting++;
  I->Setting[cSetting_ribbon_power].Value[0] = 5;
  strcpy(I->Setting[cSetting_ribbon_power].Name,
			"ribbon_power");

  I->NSetting++;
  I->Setting[cSetting_ribbon_power_b].Value[0] = 2.0;
  strcpy(I->Setting[cSetting_ribbon_power_b].Name,
			"ribbon_power_b");

  I->NSetting++;
  I->Setting[cSetting_ribbon_sampling].Value[0] = 16;
  strcpy(I->Setting[cSetting_ribbon_sampling].Name,
			"ribbon_sampling");

  I->NSetting++;
  I->Setting[cSetting_ribbon_radius].Value[0] = 0.4;
  strcpy(I->Setting[cSetting_ribbon_radius].Name,
			"ribbon_radius");

  I->NSetting++;
  I->Setting[cSetting_stick_radius].Value[0] = 0.3;
  strcpy(I->Setting[cSetting_stick_radius].Name,
			"stick_radius");

  I->NSetting++;
  I->Setting[cSetting_hash_max].Value[0] = 80;
  strcpy(I->Setting[cSetting_hash_max].Name,
			"hash_max");

  I->NSetting++;
  I->Setting[cSetting_ortho].Value[0] = 0;
  strcpy(I->Setting[cSetting_ortho].Name,
			"ortho");

  I->NSetting++;
  I->Setting[cSetting_power].Value[0] = 3.0;
  strcpy(I->Setting[cSetting_power].Name,
			"power");

  I->NSetting++;
  I->Setting[cSetting_spec_reflect].Value[0] = 0.4;
  strcpy(I->Setting[cSetting_spec_reflect].Name,
			"spec_reflect");

  I->NSetting++;
  I->Setting[cSetting_spec_power].Value[0] = 40;
  strcpy(I->Setting[cSetting_spec_power].Name,
			"spec_power");


  I->NSetting++;
  I->Setting[cSetting_sweep_angle].Value[0] = 15.0;
  strcpy(I->Setting[cSetting_sweep_angle].Name,
			"sweep_angle");

  I->NSetting++;
  I->Setting[cSetting_sweep_speed].Value[0] = 1.0;
  strcpy(I->Setting[cSetting_sweep_speed].Name,
			"sweep_speed");

  I->NSetting++;
  I->Setting[cSetting_dot_hydrogens].Value[0] = 1.0;
  strcpy(I->Setting[cSetting_dot_hydrogens].Name,
			"dot_hydrogens");

  I->NSetting++;
  I->Setting[cSetting_dot_size].Value[0] = 0.06;
  strcpy(I->Setting[cSetting_dot_size].Name,
			"dot_size");

  I->NSetting++;
  I->Setting[cSetting_ray_trace_frames].Value[0] = 0.0;
  strcpy(I->Setting[cSetting_ray_trace_frames].Name,
			"ray_trace_frames");

  I->NSetting++;
  I->Setting[cSetting_cache_frames].Value[0] = 0.0;
  strcpy(I->Setting[cSetting_cache_frames].Name,
			"cache_frames");

  I->NSetting++;
  I->Setting[cSetting_trim_dots].Value[0] = 1.0;
  strcpy(I->Setting[cSetting_trim_dots].Name,
			"trim_dots");

  I->NSetting++;
  I->Setting[cSetting_cull_spheres].Value[0] = 1.0;
  strcpy(I->Setting[cSetting_cull_spheres].Name,
			"cull_spheres");

  I->NSetting++;
  I->Setting[cSetting_test1].Value[0] = 1000000.0;
  strcpy(I->Setting[cSetting_test1].Name,
			"test1");

  I->NSetting++;
  I->Setting[cSetting_test2].Value[0] = 1.0;
  strcpy(I->Setting[cSetting_test2].Name,
			"test2");

  I->NSetting++;
  I->Setting[cSetting_surface_best].Value[0] = 0.2;
  strcpy(I->Setting[cSetting_surface_best].Name,
			"surface_best");

  I->NSetting++;
  I->Setting[cSetting_surface_normal].Value[0] = 0.5;
  strcpy(I->Setting[cSetting_surface_normal].Name,
			"surface_normal");

  I->NSetting++;
  I->Setting[cSetting_surface_quality].Value[0] = 0.0;
  strcpy(I->Setting[cSetting_surface_quality].Name,
			"surface_quality");

  I->NSetting++;
  I->Setting[cSetting_surface_proximity].Value[0] = 1.5;
  strcpy(I->Setting[cSetting_surface_proximity].Name,
			"surface_proximity");

  I->NSetting++;
  I->Setting[cSetting_stereo_angle].Value[0] = 2.1;
  strcpy(I->Setting[cSetting_stereo_angle].Name,
			"stereo_angle");

  I->NSetting++;
  I->Setting[cSetting_stereo_shift].Value[0] = 2.0;
  strcpy(I->Setting[cSetting_stereo_shift].Name,
			"stereo_shift");

  I->NSetting++;
  I->Setting[cSetting_line_smooth].Value[0] = 0.0;
  strcpy(I->Setting[cSetting_line_smooth].Name,
			"line_smooth");

  I->NSetting++;
  I->Setting[cSetting_line_width].Value[0] = 1.0;
  strcpy(I->Setting[cSetting_line_width].Name,
			"line_width");

  I->NSetting++;
  I->Setting[cSetting_half_bonds].Value[0] = 0.0;
  strcpy(I->Setting[cSetting_half_bonds].Name,
			"half_bonds");

  I->NSetting++;
  I->Setting[cSetting_stick_quality].Value[0] = 8.0;
  strcpy(I->Setting[cSetting_stick_quality].Name,
			"stick_quality");

  I->NSetting++;
  I->Setting[cSetting_stick_overlap].Value[0] = 0.2;
  strcpy(I->Setting[cSetting_stick_overlap].Name,
			"stick_overlap");

  I->NSetting++;
  I->Setting[cSetting_stick_nub].Value[0] = 0.7;
  strcpy(I->Setting[cSetting_stick_nub].Name,
			"stick_nub");

  I->NSetting++;
  I->Setting[cSetting_all_states].Value[0] = 0.0;
  strcpy(I->Setting[cSetting_all_states].Name,
			"all_states");

  I->NSetting++;
  I->Setting[cSetting_pickable].Value[0] = 1.0;
  strcpy(I->Setting[cSetting_pickable].Name,
			"pickable");

  I->NSetting++;
  I->Setting[cSetting_autoshow_lines].Value[0] = 1.0;
  strcpy(I->Setting[cSetting_autoshow_lines].Name,
			"autoshow_lines");

  I->NSetting++;
  I->Setting[cSetting_fast_idle].Value[0] = 20000;
  strcpy(I->Setting[cSetting_fast_idle].Name,
			"fast_idle");

  I->NSetting++;
  I->Setting[cSetting_no_idle].Value[0] = 5000;
  strcpy(I->Setting[cSetting_no_idle].Name,
			"no_idle");

  I->NSetting++;
  I->Setting[cSetting_slow_idle].Value[0] = 200000;
  strcpy(I->Setting[cSetting_slow_idle].Name,
			"slow_idle");

  I->NSetting++;
  I->Setting[cSetting_idle_delay].Value[0] = 1.5;
  strcpy(I->Setting[cSetting_idle_delay].Name,
			"idle_delay");

  I->NSetting++;
  I->Setting[cSetting_rock_delay].Value[0] = 30;
  strcpy(I->Setting[cSetting_rock_delay].Name,
			"rock_delay");

  I->NSetting++;
  I->Setting[cSetting_dist_counter].Value[0] = 0.0;
  strcpy(I->Setting[cSetting_dist_counter].Name,
			"dist_counter");

  I->NSetting++;
  I->Setting[cSetting_dash_length].Value[0] = 0.2;
  strcpy(I->Setting[cSetting_dash_length].Name,
			"dash_length");

  I->NSetting++;
  I->Setting[cSetting_dash_gap].Value[0] = 0.4;
  strcpy(I->Setting[cSetting_dash_gap].Name,
			"dash_gap");

  I->NSetting++;
  I->Setting[cSetting_auto_zoom].Value[0] = 1.0;
  strcpy(I->Setting[cSetting_auto_zoom].Name,
			"auto_zoom");

  I->NSetting++;
  I->Setting[cSetting_overlay].Value[0] = 0.0;
  strcpy(I->Setting[cSetting_overlay].Name,
			"overlay");

  I->NSetting++;
  I->Setting[cSetting_text].Value[0] = 0.0;
  strcpy(I->Setting[cSetting_text].Name,
			"text");

  I->NSetting++;
  I->Setting[cSetting_button_mode].Value[0] = 0.0;
  strcpy(I->Setting[cSetting_button_mode].Name,
			"button_mode"); 

  I->NSetting++;
  I->Setting[cSetting_valence].Value[0] = 0.0;
  strcpy(I->Setting[cSetting_valence].Name,
			"valence");

  I->NSetting++;
  I->Setting[cSetting_nonbonded_size].Value[0] = 0.25;
  strcpy(I->Setting[cSetting_nonbonded_size].Name,
			"nonbonded_size");

  I->NSetting++;
  I->Setting[cSetting_label_color].Value[0] = -1.0;
  strcpy(I->Setting[cSetting_label_color].Name,
			"label_color");

  I->NSetting++;
  I->Setting[cSetting_ray_trace_fog].Value[0] = 1.0;
  strcpy(I->Setting[cSetting_ray_trace_fog].Name,
			"ray_trace_fog");

  I->NSetting++;
  I->Setting[cSetting_spheroid_scale].Value[0] = 1.0;
  strcpy(I->Setting[cSetting_spheroid_scale].Name,
			"spheroid_scale");

  I->NSetting++;
  I->Setting[cSetting_ray_trace_fog_start].Value[0] = 0.30;
  strcpy(I->Setting[cSetting_ray_trace_fog_start].Name,
			"ray_trace_fog_start");

  I->NSetting++;
  I->Setting[cSetting_spheroid_smooth].Value[0] = 1.1;
  strcpy(I->Setting[cSetting_spheroid_smooth].Name,
			"spheroid_smooth");

  I->NSetting++;
  I->Setting[cSetting_spheroid_fill].Value[0] = 1.30;
  strcpy(I->Setting[cSetting_spheroid_fill].Name,
			"spheroid_fill");

  I->NSetting++;
  I->Setting[cSetting_autoshow_nonbonded].Value[0] = 1.0;
  strcpy(I->Setting[cSetting_autoshow_nonbonded].Name,
			"autoshow_nonbonded");

  I->NSetting++;
#ifdef WIN32
  I->Setting[cSetting_cache_display].Value[0] = 0.0;
#else
  I->Setting[cSetting_cache_display].Value[0] = 1.0;
#endif
  strcpy(I->Setting[cSetting_cache_display].Name,
			"cache_display");

  I->NSetting++;
  I->Setting[cSetting_normal_workaround].Value[0] = 0;
  strcpy(I->Setting[cSetting_normal_workaround].Name,
			"normal_workaround");

}


