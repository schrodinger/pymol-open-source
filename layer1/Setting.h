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
#ifndef _H_Setting
#define _H_Setting

typedef char SettingName[255];

typedef struct {
  SettingName Name;
  float Value[3];
} SettingRec;

typedef struct  {
  SettingRec *Setting;
  int NSetting;
} CSetting;

void SettingInit(void);
void SettingFree(void);

int SettingGetIndex(char *name);
float SettingGet(int index);
void SettingSet(int index,float v);
void SettingSetfv(int index,float *value);
float *SettingGetfv(int index);
void SettingSetNamed(char *name,char *value);
float SettingGetNamed(char *name);

#define cSetting_bonding_vdw_cutoff    0
#define cSetting_min_mesh_spacing      1
#define cSetting_dot_density           2
#define cSetting_dot_surface           3
#define cSetting_solvent_radius        4
#define cSetting_sel_counter           5
#define cSetting_bg_rgb                6
#define cSetting_ambient               7
#define cSetting_direct                8
#define cSetting_reflect               9
#define cSetting_light                10
#define cSetting_power                11
#define cSetting_antialias            12
#define cSetting_cavity_cull          13
#define cSetting_ambient_scale        14
#define cSetting_single_image         15
#define cSetting_min_delay            16
#define cSetting_ribbon_power         17
#define cSetting_ribbon_power_b       18
#define cSetting_ribbon_sampling      19
#define cSetting_ribbon_radius        20
#define cSetting_stick_radius         21
#define cSetting_hash_max             22
#define cSetting_ortho                23
#define cSetting_spec_reflect         24
#define cSetting_spec_power           25
#define cSetting_sweep_angle          26
#define cSetting_sweep_speed          27
#define cSetting_dot_hydrogens        28
#define cSetting_dot_size             29
#define cSetting_ray_trace_frames     30
#define cSetting_cache_frames         31
#define cSetting_trim_dots            32
#define cSetting_cull_spheres         33
#define cSetting_test1                34
#define cSetting_test2                35
#define cSetting_surface_best         36
#define cSetting_surface_normal       37
#define cSetting_surface_quality      38
#define cSetting_surface_proximity    39
#define cSetting_normal_workaround    40 

#endif

