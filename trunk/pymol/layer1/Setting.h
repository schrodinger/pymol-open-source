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

#include<Python.h>

typedef char SettingName[255];

typedef struct {
  int defined;
  int changed;
  int type;
  unsigned int offset;
} SettingRec;

typedef struct {
  unsigned int size;
  char *data;
  SettingRec *info;
} CSetting;

#define cSetting_boolean     1
#define cSetting_int         2
#define cSetting_float       3
#define cSetting_float3      4
#define cSetting_color       5

/* New API 
 * NOTE: get commands are not range-checked, so be careful
 * in contrast, set commands expand the current list 
 */

void SettingInitGlobal(void);
void SettingFreeGlobal(void);

CSetting *SettingNew(void);
void SettingFreeP(CSetting *I);
void SettingInit(CSetting *I);
void SettingPurge(CSetting *I);
void SettingCheckHandle(CSetting **handle);

int SettingSet_b(CSetting *I,int index, int value);
int SettingSet_i(CSetting *I,int index, int value);
int SettingSet_f(CSetting *I,int index, float value);
int SettingSet_3f(CSetting *I,int index, float value1,float value2,float value3);
int SettingSet_3fv(CSetting *I,int index, float *value);

int SettingGetTextValue(CSetting *set1,CSetting *set2,int index,char *buffer);


int SettingSetTuple(CSetting *I,int index,PyObject *tuple);

void SettingClear(CSetting *I,int index); /* don't call this for the global list! */

int SettingGetType(int index); /* based on global types, always succeeds */

int   SettingGetGlobal_b(int index); /* always succeed */
int   SettingGetGlobal_i(int index); /* always succeed */
float SettingGetGlobal_f(int index); /* always succeed */
void  SettingGetGlobal_3f(int index,float *value); /* always succeeds */
float *SettingGetGlobal_fv(int index); /* always succeed */
int SettingSet_color(CSetting *I,int index, char *value);

int   SettingGet_b  (CSetting *set1,CSetting *set2,int index);
int   SettingGet_i  (CSetting *set1,CSetting *set2,int index);
float SettingGet_f  (CSetting *set1,CSetting *set2,int index);
void  SettingGet_3f (CSetting *set1,CSetting *set2,int index,float *value);
float *SettingGet_fv (CSetting *set1,CSetting *set2,int index);
int   SettingGet_color(CSetting *set1,CSetting *set2,int index);

PyObject *SettingGetTuple(CSetting *set1,CSetting *set2,int index); /* (type,(value,)) */

void SettingGenerateSideEffects(int index,char *sele,int state);
PyObject *SettingGetUpdateList(CSetting *I);

/* Legacy API below */

int SettingGetIndex(char *name);
float SettingGet(int index);
int SettingSet(int index,float v);
int SettingSetfv(int index,float *value);
float *SettingGetfv(int index);
int SettingSetNamed(char *name,char *value);
float SettingGetNamed(char *name);
int SettingGetName(int index,SettingName name);

#define cSetting_bonding_vdw_cutoff            0
#define cSetting_min_mesh_spacing              1
#define cSetting_dot_density                   2
#define cSetting_dot_mode                      3
#define cSetting_solvent_radius                4
#define cSetting_sel_counter                   5
#define cSetting_bg_rgb                        6
#define cSetting_ambient                       7
#define cSetting_direct                        8
#define cSetting_reflect                       9
#define cSetting_light                        10
#define cSetting_power                        11
#define cSetting_antialias                    12
#define cSetting_cavity_cull                  13
#define cSetting_gl_ambient                   14
#define cSetting_single_image                 15
#define cSetting_movie_delay                  16
#define cSetting_ribbon_power                 17
#define cSetting_ribbon_power_b               18
#define cSetting_ribbon_sampling              19
#define cSetting_ribbon_radius                20
#define cSetting_stick_radius                 21
#define cSetting_hash_max                     22
#define cSetting_ortho                        23
#define cSetting_spec_reflect                 24
#define cSetting_spec_power                   25
#define cSetting_sweep_angle                  26
#define cSetting_sweep_speed                  27
#define cSetting_dot_hydrogens                28
#define cSetting_dot_radius                   29
#define cSetting_ray_trace_frames             30
#define cSetting_cache_frames                 31
#define cSetting_trim_dots                    32
#define cSetting_cull_spheres                 33
#define cSetting_test1                        34
#define cSetting_test2                        35
#define cSetting_surface_best                 36
#define cSetting_surface_normal               37
#define cSetting_surface_quality              38
#define cSetting_surface_proximity            39
#define cSetting_normal_workaround            40
#define cSetting_stereo_angle                 41
#define cSetting_stereo_shift                 42
#define cSetting_line_smooth                  43
#define cSetting_line_width                   44
#define cSetting_half_bonds                   45
#define cSetting_stick_quality                46
#define cSetting_stick_overlap                47
#define cSetting_stick_nub                    48
#define cSetting_all_states                   49
#define cSetting_pickable                     50
#define cSetting_auto_show_lines              51
#define cSetting_idle_delay                   52
#define cSetting_no_idle                      53
#define cSetting_fast_idle                    54
#define cSetting_slow_idle                    55
#define cSetting_rock_delay                   56
#define cSetting_dist_counter                 57
#define cSetting_dash_length                  58
#define cSetting_dash_gap                     59
#define cSetting_auto_zoom                    60
#define cSetting_overlay                      61
#define cSetting_text                         62
#define cSetting_button_mode                  63
#define cSetting_valence                      64
#define cSetting_nonbonded_size               65
#define cSetting_label_color                  66
#define cSetting_ray_trace_fog                67
#define cSetting_spheroid_scale               68
#define cSetting_ray_trace_fog_start          69
#define cSetting_spheroid_smooth              70
#define cSetting_spheroid_fill                71
#define cSetting_auto_show_nonbonded          72
#define cSetting_cache_display                73
#define cSetting_mesh_radius                  74
#define cSetting_backface_cull                75
#define cSetting_gamma                        76
#define cSetting_dot_width                    77
#define cSetting_auto_show_selections         78
#define cSetting_auto_hide_selections         79
#define cSetting_selection_width              80
#define cSetting_selection_overlay            81
#define cSetting_static_singletons            82
#define cSetting_max_triangles                83
#define cSetting_depth_cue                    84
#define cSetting_specular                     85
#define cSetting_shininess                    86
#define cSetting_sphere_quality               87
#define cSetting_fog                          88
#define cSetting_isomesh_auto_state           89
#define cSetting_mesh_width                   90
#define cSetting_cartoon_sampling             91
#define cSetting_cartoon_loop_radius          92
#define cSetting_cartoon_loop_quality         93
#define cSetting_cartoon_power                94
#define cSetting_cartoon_power_b              95
#define cSetting_cartoon_rect_length          96
#define cSetting_cartoon_rect_width           97
#define cSetting_internal_gui_width           98
#define cSetting_internal_gui                 99
#define cSetting_cartoon_oval_length         100
#define cSetting_cartoon_oval_width          101
#define cSetting_cartoon_oval_quality        102
#define cSetting_cartoon_tube_radius         103
#define cSetting_cartoon_tube_quality        104
#define cSetting_cartoon_debug               105
#define cSetting_ribbon_width                106
#define cSetting_dash_width                  107
#define cSetting_dash_radius                 108
#define cSetting_cgo_ray_width_scale         109
#define cSetting_line_radius                 110
#define cSetting_cartoon_round_helices       111
#define cSetting_cartoon_refine_normals      112
#define cSetting_cartoon_flat_sheets         113
#define cSetting_cartoon_smooth_loops        114
#define cSetting_cartoon_dumbbell_length     115
#define cSetting_cartoon_dumbbell_width      116
#define cSetting_cartoon_dumbbell_radius     117
#define cSetting_cartoon_fancy_helices       118
#define cSetting_cartoon_fancy_sheets        119
#define cSetting_ignore_pdb_segi             120
#define cSetting_ribbon_throw                121
#define cSetting_cartoon_throw               122
#define cSetting_cartoon_refine              123
#define cSetting_cartoon_refine_tips         124
#define cSetting_cartoon_discrete_colors     125
#define cSetting_normalize_ccp4_maps         126
#define cSetting_surface_poor                127
#define cSetting_internal_feedback           128
#define cSetting_cgo_line_width              129
#define cSetting_cgo_line_radius             130
#define cSetting_logging                     131
#define cSetting_robust_logs                 132
#define cSetting_log_box_selections          133
#define cSetting_log_conformations           134
#define cSetting_valence_default             135
#define cSetting_surface_miserable           136
#define cSetting_ray_opaque_background       137
#define cSetting_transparency                138
#define cSetting_ray_texture                 139
#define cSetting_ray_texture_settings        140
#define cSetting_suspend_updates             141
#define cSetting_full_screen                 142
#define cSetting_surface_mode                143
#define cSetting_surface_color               144
#define cSetting_mesh_mode                   145
#define cSetting_mesh_color                  146
#define cSetting_auto_indicate_flags         147
#define cSetting_surface_debug               148
#define cSetting_ray_improve_shadows         149
#define cSetting_smooth_color_triangle       150
#define cSetting_ray_default_renderer        151
#define cSetting_field_of_view               152
#define cSetting_reflect_power               153
#define cSetting_preserve_chempy_ids         154
#define cSetting_sphere_scale                155
#define cSetting_two_sided_lighting          156
#define cSetting_secondary_structure         157
#define cSetting_auto_remove_hydrogens       158
#define cSetting_raise_exceptions            159
#define cSetting_stop_on_exceptions          160
#define cSetting_sculpting                   161
#define cSetting_auto_sculpt                 162
#define cSetting_sculpt_vdw                  163
#define cSetting_sculpt_vdw14                164
#define cSetting_sculpting_cycles            165

#define cSetting_INIT                        166

#endif

