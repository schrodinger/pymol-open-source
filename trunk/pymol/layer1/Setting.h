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

#include"os_python.h"
#include"PyMOLGlobals.h"
#include"OVOneToOne.h"

typedef char SettingName[255];

/* for atomic settings */


typedef struct {
  int setting_id;
  int type; /* must be cSetting_boolean, cSetting_int, cSetting_float, or cSetting_color */
  union {
  int int_;
  float float_;
  } value;
  int next; /* for per-atom setting lists & memory management */
} SettingUniqueEntry;

struct _CSettingUnique {
  OVOneToOne *id2offset;
  OVOneToOne *old2new;
  SettingUniqueEntry *entry;
  int n_alloc, next_free;
};

typedef struct {
  int defined;
  int changed;
  int type;
  ov_diff offset;
  ov_size max_size;
} SettingRec;

struct _CSetting {
  PyMOLGlobals *G;
  ov_size size;
  char *data;
  SettingRec *info;
};

#define cSetting_blank       0
#define cSetting_boolean     1
#define cSetting_int         2
#define cSetting_float       3
#define cSetting_float3      4
#define cSetting_color       5
#define cSetting_string      6

/* Atomic Settings */

void SettingUniqueDetachChain(PyMOLGlobals *G,int index);
/* New API 
 * NOTE: get commands are not range-checked, so be careful
 * in contrast, set commands expand the current list 
 */

void SettingUniqueSet_b(PyMOLGlobals *G,int unique_id,int setting_id,int value);
void SettingUniqueSet_i(PyMOLGlobals *G,int unique_id,int setting_id,int value);
void SettingUniqueSet_f(PyMOLGlobals *G,int unique_id,int setting_id,float value);
void SettingUniqueSet_color(PyMOLGlobals *G,int unique_id,int setting_id,int value);
void SettingUniqueSetTypedValue(PyMOLGlobals *G,int unique_id,int setting_id,int setting_type, void *value);

int SettingUniqueCheck(PyMOLGlobals *G,int unique_id,int setting_id);
int SettingUniqueGet_b(PyMOLGlobals *G,int unique_id,int setting_id,int *value);
int SettingUniqueGet_i(PyMOLGlobals *G,int unique_id,int setting_id,int *value);
int SettingUniqueGet_f(PyMOLGlobals *G,int unique_id,int setting_id,float *value);
int SettingUniqueGet_color(PyMOLGlobals *G,int unique_id,int setting_id,int *value);

void SettingUniqueResetAll(PyMOLGlobals *G);
PyObject *SettingUniqueAsPyList(PyMOLGlobals *G);
int SettingUniqueFromPyList(PyMOLGlobals *G,PyObject *list,int partial_restore);
int SettingUniqueConvertOldSessionID(PyMOLGlobals *G,int old_unique_id);

int SettingUniqueCopyAll(PyMOLGlobals *G, int src_unique_id, int dst_unique_id);
void SettingInitGlobal(PyMOLGlobals *G,int alloc,int reset_gui,int use_default);
void SettingStoreDefault(PyMOLGlobals *G);
void SettingPurgeDefault(PyMOLGlobals *G);

void SettingFreeGlobal(PyMOLGlobals *G);


CSetting *SettingNew(PyMOLGlobals *G);
void SettingFreeP(CSetting *I);
void SettingInit(PyMOLGlobals *G,CSetting *I);
void SettingPurge(CSetting *I);
void SettingCheckHandle(PyMOLGlobals *G,CSetting **handle);

int SettingSet_b(CSetting *I,int index, int value);
int SettingSet_i(CSetting *I,int index, int value);
int SettingSet_f(CSetting *I,int index, float value);
int SettingSet_s(CSetting *I,int index, char *value);
int SettingSet_3f(CSetting *I,int index, float value1,float value2,float value3);
int SettingSet_3fv(CSetting *I,int index, float *value);

int SettingGetTextValue(PyMOLGlobals *G,CSetting *set1,CSetting *set2,int index,char *buffer);

int SettingUnset(CSetting *I,int index);

void SettingClear(CSetting *I,int index); /* don't call this for the global list! */

int SettingGetType(PyMOLGlobals *G,int index); /* based on global types, always succeeds */

int   SettingGetGlobal_b(PyMOLGlobals *G,int index); /* always succeed */
int   SettingGetGlobal_i(PyMOLGlobals *G,int index); /* always succeed */
float SettingGetGlobal_f(PyMOLGlobals *G,int index); /* always succeed */
char *SettingGetGlobal_s(PyMOLGlobals *G,int index); /* always succeeds */
int   SettingGetGlobal_color(PyMOLGlobals *G,int index); /* always succeed */

void  SettingGetGlobal_3f(PyMOLGlobals *G,int index,float *value); /* always succeeds */
float *SettingGetGlobal_3fv(PyMOLGlobals *G,int index); /* always succeed */

int   SettingSetGlobal_b(PyMOLGlobals *G,int index,int value);
int   SettingSetGlobal_i(PyMOLGlobals *G,int index,int value);
int   SettingSetGlobal_f(PyMOLGlobals *G,int index,float value);
int   SettingSetGlobal_3f(PyMOLGlobals *G,int index, float value1,float value2,float value3);
int   SettingSetSmart_i(PyMOLGlobals *G,CSetting *set1,CSetting *set2,int index, int value);
/* more to come */

int SettingGetIfDefined_i(PyMOLGlobals *G,CSetting *set1,int index,int *value);
int SettingGetIfDefined_b(PyMOLGlobals *G,CSetting *set1,int index,int *value);
int SettingGetIfDefined_f(PyMOLGlobals *G,CSetting *set1,int index,float *value);
int SettingGetIfDefined_s(PyMOLGlobals *G,CSetting *set1,int index,char **value);
int SettingGetIfDefined_3fv(PyMOLGlobals *G,CSetting *set1,int index,float **value);
int SettingGetIfDefined_color(PyMOLGlobals *G,CSetting *set1,int index,int *value);

/* more to come */

int SettingSet_color(CSetting *I,int index, char *value);


int   SettingGet_b(PyMOLGlobals *G,CSetting *set1,CSetting *set2,int index);
int   SettingGet_i(PyMOLGlobals *G,CSetting *set1,CSetting *set2,int index);
float SettingGet_f(PyMOLGlobals *G,CSetting *set1,CSetting *set2,int index);
char  *SettingGet_s(PyMOLGlobals *G,CSetting *set1,CSetting *set2,int index);
void  SettingGet_3f(PyMOLGlobals *G,CSetting *set1,CSetting *set2,int index,float *value);
float *SettingGet_3fv(PyMOLGlobals *G,CSetting *set1,CSetting *set2,int index);
int   SettingGet_color(PyMOLGlobals *G,CSetting *set1,CSetting *set2,int index);

int SettingSetFromString(PyMOLGlobals *G,CSetting *I,int index,char *st);
int SettingStringToTypedValue(PyMOLGlobals *G,int index,char *st, int *type, int *value);

#ifndef _PYMOL_NOPY
int SettingSetFromTuple(PyMOLGlobals *G,CSetting *I,int index,PyObject *tuple);
PyObject *SettingGetTuple(PyMOLGlobals *G,CSetting *set1,CSetting *set2,int index); /* (type,(value,)) */
PyObject *SettingGetDefinedTuple(PyMOLGlobals *G,CSetting *set1,int index);
PyObject *SettingGetUpdateList(PyMOLGlobals *G,CSetting *I);
#endif

void SettingGenerateSideEffects(PyMOLGlobals *G,int index,char *sele,int state);

/* Legacy API below */

int SettingGetIndex(PyMOLGlobals *G,char *name);
float SettingGet(PyMOLGlobals *G,int index);
int SettingSet(PyMOLGlobals *G,int index,float v);
int SettingSetfv(PyMOLGlobals *G,int index,float *value);
float *SettingGetfv(PyMOLGlobals *G,int index);
int SettingSetNamed(PyMOLGlobals *G,char *name,char *value);
float SettingGetNamed(PyMOLGlobals *G,char *name);
int SettingGetName(PyMOLGlobals *G,int index,SettingName name);

PyObject *SettingAsPyList(CSetting *I);
int SettingFromPyList(CSetting *I,PyObject *list);

int SettingSetGlobalsFromPyList(PyMOLGlobals *G,PyObject *list);
PyObject *SettingGetGlobalsAsPyList(PyMOLGlobals *G);

/* proposed...
PyObject *SettingGetDefaultsAsPyList(PyMOLGlobals *G);
int SettingSetDefaultsFromPyList(PyMOLGlobals *G,PyObject *list);
*/

CSetting *SettingNewFromPyList(PyMOLGlobals *G,PyObject *list);

/* WARNING: do not delete or change indices
   since they are used in session objects */

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
#define cSetting_valence_size                135
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
#define cSetting_sculpt_vdw_scale            163
#define cSetting_sculpt_vdw_scale14          164
#define cSetting_sculpt_vdw_weight           165
#define cSetting_sculpt_vdw_weight14         166
#define cSetting_sculpt_bond_weight          167
#define cSetting_sculpt_angl_weight          168
#define cSetting_sculpt_pyra_weight          169
#define cSetting_sculpt_plan_weight          170
#define cSetting_sculpting_cycles            171
#define cSetting_sphere_transparency         172
#define cSetting_sphere_color                173
#define cSetting_sculpt_field_mask           174
#define cSetting_sculpt_hb_overlap           175
#define cSetting_sculpt_hb_overlap_base      176
#define cSetting_legacy_vdw_radii            177
#define cSetting_sculpt_memory               178
#define cSetting_connect_mode                179
#define cSetting_cartoon_cylindrical_helices 180
#define cSetting_cartoon_helix_radius        181
#define cSetting_connect_cutoff              182
#define cSetting_save_pdb_ss                 183
#define cSetting_sculpt_line_weight          184
#define cSetting_fit_iterations              185
#define cSetting_fit_tolerance               186
#define cSetting_batch_prefix                187

#define cSetting_stereo_mode                 188

#define cStereo_default              0 
#define cStereo_quadbuffer           1
#define cStereo_crosseye             2 
#define cStereo_walleye              3 
#define cStereo_geowall              4 
#define cStereo_sidebyside           5
#define cStereo_stencil_by_row       6 
#define cStereo_stencil_by_column    7 
#define cStereo_stencil_checkerboard 8 
#define cStereo_stencil_custom       9 /* for hardware developers to use */
#define cStereo_anaglyph            10 /* not yet implemented */
#define cStereo_dynamic             11 /* dynamic polarization */
#define cStereo_clone_dynamic       12
   
#define cSetting_cgo_sphere_quality          189
#define cSetting_pdb_literal_names           190
#define cSetting_wrap_output                 191
#define cSetting_fog_start                   192
#define cSetting_state                       193
#define cSetting_frame                       194
#define cSetting_ray_shadows                 195
#define cSetting_ribbon_trace_atoms          196
#define cSetting_security                    197
#define cSetting_stick_transparency          198 
#define cSetting_ray_transparency_shadows    199
#define cSetting_session_version_check       200
#define cSetting_ray_transparency_specular   201
#define cSetting_stereo_double_pump_mono     202
#define cSetting_sphere_solvent              203
#define cSetting_mesh_quality                204
#define cSetting_mesh_solvent                205
#define cSetting_dot_solvent                 206
#define cSetting_ray_shadow_fudge            207
#define cSetting_ray_triangle_fudge          208
#define cSetting_debug_pick                  209
#define cSetting_dot_color                   210
#define cSetting_mouse_limit                 211
#define cSetting_mouse_scale                 212
#define cSetting_transparency_mode           213
#define cSetting_clamp_colors                214
#define cSetting_pymol_space_max_red         215
#define cSetting_pymol_space_max_green       216
#define cSetting_pymol_space_max_blue        217
#define cSetting_pymol_space_min_factor      218
#define cSetting_roving_origin               219
#define cSetting_roving_lines                220
#define cSetting_roving_sticks               221
#define cSetting_roving_spheres              222
#define cSetting_roving_labels               223
#define cSetting_roving_delay                224
#define cSetting_roving_selection            225
#define cSetting_roving_byres                226
#define cSetting_roving_ribbon               227
#define cSetting_roving_cartoon              228
#define cSetting_roving_polar_contacts       229
#define cSetting_roving_polar_cutoff         230
#define cSetting_roving_nonbonded            231
#define cSetting_float_labels                232
#define cSetting_roving_detail               233
#define cSetting_roving_nb_spheres           234
#define cSetting_ribbon_color                235
#define cSetting_cartoon_color               236
#define cSetting_ribbon_smooth               237
#define cSetting_auto_color                  238
#define cSetting_auto_color_next             239
#define cSetting_ray_interior_color          240
#define cSetting_cartoon_highlight_color     241
#define cSetting_coulomb_units_factor        242
#define cSetting_coulomb_dielectric          243
#define cSetting_ray_interior_shadows        244
#define cSetting_ray_interior_texture        245

#define cSetting_roving_map1_name            246
#define cSetting_roving_map2_name            247
#define cSetting_roving_map3_name            248

#define cSetting_roving_map1_level           249
#define cSetting_roving_map2_level           250
#define cSetting_roving_map3_level           251

#define cSetting_roving_isomesh              252
#define cSetting_roving_isosurface           253
#define cSetting_scenes_changed              254

#define cSetting_gaussian_b_adjust           255
#define cSetting_pdb_standard_order          256

#define cSetting_cartoon_smooth_first        257
#define cSetting_cartoon_smooth_last         258
#define cSetting_cartoon_smooth_cycles       259
#define cSetting_cartoon_flat_cycles         260

#define cSetting_max_threads                 261
#define cSetting_show_progress               262
#define cSetting_use_display_lists           263
#define cSetting_cache_memory                264
#define cSetting_simplify_display_lists      265
#define cSetting_retain_order                266
#define cSetting_pdb_hetatm_sort             267
#define cSetting_pdb_use_ter_records         268
#define cSetting_cartoon_trace_atoms         269
#define cSetting_ray_oversample_cutoff       270
#define cSetting_gaussian_resolution         271
#define cSetting_gaussian_b_floor            272
#define cSetting_sculpt_nb_interval          273
#define cSetting_sculpt_tors_weight          274
#define cSetting_sculpt_tors_tolerance       275
#define cSetting_stick_ball                  276
#define cSetting_stick_ball_ratio            277
#define cSetting_stick_fixed_radius          278
#define cSetting_cartoon_transparency        279
#define cSetting_dash_round_ends             280
#define cSetting_h_bond_max_angle            281
#define cSetting_h_bond_cutoff_center        282
#define cSetting_h_bond_cutoff_edge          283
#define cSetting_h_bond_power_a              284
#define cSetting_h_bond_power_b              285
#define cSetting_h_bond_cone                 286

#define cSetting_ss_helix_psi_target         287 
#define cSetting_ss_helix_psi_include        288
#define cSetting_ss_helix_psi_exclude        289

#define cSetting_ss_helix_phi_target         290
#define cSetting_ss_helix_phi_include        291
#define cSetting_ss_helix_phi_exclude        292

#define cSetting_ss_strand_psi_target          293
#define cSetting_ss_strand_psi_include         294
#define cSetting_ss_strand_psi_exclude         295

#define cSetting_ss_strand_phi_target          296
#define cSetting_ss_strand_phi_include         297
#define cSetting_ss_strand_phi_exclude         298
#define cSetting_movie_loop                    299

#define cSetting_pdb_retain_ids             300
#define cSetting_pdb_no_end_record          301
#define cSetting_cgo_dot_width              302
#define cSetting_cgo_dot_radius             303
#define cSetting_defer_updates              304
#define cSetting_normalize_o_maps           305
#define cSetting_swap_dsn6_bytes            306
#define cSetting_pdb_insertions_go_first    307
#define cSetting_roving_origin_z            308
#define cSetting_roving_origin_z_cushion    309
#define cSetting_specular_intensity         310
#define cSetting_overlay_lines              311
#define cSetting_ray_transparency_spec_cut  312
#define cSetting_internal_prompt            313
#define cSetting_normalize_grd_maps         314
#define cSetting_ray_blend_colors           315
#define cSetting_ray_blend_red              316
#define cSetting_ray_blend_green            317
#define cSetting_ray_blend_blue             318
#define cSetting_png_screen_gamma           319
#define cSetting_png_file_gamma             320
#define cSetting_editor_label_fragments     321
#define cSetting_internal_gui_control_size  322
#define cSetting_auto_dss                   323
#define cSetting_transparency_picking_mode  324
#define cSetting_virtual_trackball          325
#define cSetting_pdb_reformat_names_mode    326
#define cSetting_ray_pixel_scale            327
#define cSetting_label_font_id              328
#define cSetting_pdb_conect_all             329
#define cSetting_button_mode_name           330
#define cSetting_surface_type               331
#define cSetting_dot_normals                332
#define cSetting_session_migration          333
#define cSetting_mesh_normals               334
#define cSetting_mesh_type                  335
#define cSetting_dot_lighting               336
#define cSetting_mesh_lighting              337
#define cSetting_surface_solvent            338 
#define cSetting_triangle_max_passes        339
#define cSetting_ray_interior_reflect       340
#define cSetting_internal_gui_mode          341
#define cSetting_surface_carve_selection    342
#define cSetting_surface_carve_state        343
#define cSetting_surface_carve_cutoff       344
#define cSetting_surface_clear_selection    345
#define cSetting_surface_clear_state        346
#define cSetting_surface_clear_cutoff       347
#define cSetting_surface_trim_cutoff        348
#define cSetting_surface_trim_factor        349
#define cSetting_ray_max_passes             350
#define cSetting_active_selections          351
#define cSetting_ray_transparency_contrast  352
#define cSetting_seq_view                   353
#define cSetting_mouse_selection_mode       354
#define cSetting_seq_view_label_spacing     355
#define cSetting_seq_view_label_start       356
#define cSetting_seq_view_format            357
#define cSetting_seq_view_location          358
#define cSetting_seq_view_overlay           359
#define cSetting_auto_classify_atoms        360
#define cSetting_cartoon_nucleic_acid_mode  361
#define cSetting_seq_view_color             362
#define cSetting_seq_view_label_mode        363
#define cSetting_surface_ramp_above_mode    364
#define cSetting_stereo                     365
#define cSetting_wizard_prompt_mode         366
#define cSetting_coulomb_cutoff             367
#define cSetting_slice_track_camera         368
#define cSetting_slice_height_scale         369 
#define cSetting_slice_height_map           370
#define cSetting_slice_grid                 371
#define cSetting_slice_dynamic_grid         372
#define cSetting_slice_dynamic_grid_resolution 373
#define cSetting_pdb_insure_orthogonal        374
#define cSetting_ray_direct_shade           375
#define cSetting_stick_color                376
#define cSetting_cartoon_putty_radius       377
#define cSetting_cartoon_putty_quality      378
#define cSetting_cartoon_putty_scale_min    379
#define cSetting_cartoon_putty_scale_max    380
#define cSetting_cartoon_putty_scale_power  381
#define cSetting_cartoon_putty_range        382
#define cSetting_cartoon_side_chain_helper  383
#define cSetting_surface_optimize_subsets   384
#define cSetting_multiplex                  385
#define cSetting_texture_fonts              386
#define cSetting_pqr_workarounds            387
#define cSetting_animation                  388
#define cSetting_animation_duration         389
#define cSetting_scene_animation            390
#define cSetting_line_stick_helper          391
#define cSetting_ray_orthoscopic            392
#define cSetting_ribbon_side_chain_helper   393
#define cSetting_selection_width_max        394
#define cSetting_selection_width_scale      395
#define cSetting_scene_current_name         396
#define cSetting_presentation               397
#define cSetting_presentation_mode          398
#define cSetting_pdb_truncate_residue_name  399
#define cSetting_scene_loop                 400
#define cSetting_sweep_mode                 401
#define cSetting_sweep_phase                402
#define cSetting_scene_restart_movie_delay  403
#define cSetting_mouse_restart_movie_delay  404
#define cSetting_angle_size                 405
#define cSetting_angle_label_position       406
#define cSetting_dihedral_size              407
#define cSetting_dihedral_label_position    408
#define cSetting_defer_builds_mode          409
#define cSetting_seq_view_discrete_by_state 410
#define cSetting_scene_animation_duration   411
#define cSetting_wildcard                   412
#define cSetting_atom_name_wildcard         413
#define cSetting_ignore_case                414
#define cSetting_presentation_auto_quit     415
#define cSetting_editor_auto_dihedral       416
#define cSetting_presentation_auto_start    417
#define cSetting_validate_object_names      418
#define cSetting_unused_boolean_def_true    419
#define cSetting_auto_show_spheres          420
#define cSetting_sphere_mode                421
#define cSetting_sphere_point_max_size      422
#define cSetting_sphere_point_size          423
#define cSetting_pdb_honor_model_number     424
#define cSetting_rank_assisted_sorts        425
#define cSetting_ribbon_nucleic_acid_mode   426
#define cSetting_cartoon_ring_mode          427
#define cSetting_cartoon_ring_width         428
#define cSetting_cartoon_ring_color         429
#define cSetting_cartoon_ring_finder        430
#define cSetting_cartoon_tube_cap           431
#define cSetting_cartoon_loop_cap           432
#define cSetting_nvidia_bugs                433
#define cSetting_image_dots_per_inch        434
#define cSetting_opaque_background          435
#define cSetting_draw_frames                436
#define cSetting_show_alpha_checker         437
#define cSetting_matrix_mode                438
#define cSetting_editor_auto_origin         439
#define cSetting_session_file               440
#define cSetting_cgo_transparency           441
#define cSetting_legacy_mouse_zoom          442
#define cSetting_auto_number_selections     443
#define cSetting_sculpt_vdw_vis_mode        444
#define cSetting_sculpt_vdw_vis_min         445
#define cSetting_sculpt_vdw_vis_mid         446
#define cSetting_sculpt_vdw_vis_max         447
#define cSetting_cartoon_ladder_mode        448
#define cSetting_cartoon_ladder_radius      449
#define cSetting_cartoon_ladder_color       450
#define cSetting_cartoon_nucleic_acid_color 451
#define cSetting_cartoon_ring_transparency  452
#define cSetting_label_size                 453
#define cSetting_spec_direct                454
#define cSetting_light_count                455
#define cSetting_light2                     456
#define cSetting_light3                     457
#define cSetting_hide_underscore_names      458
#define cSetting_selection_round_points     459
#define cSetting_distance_exclusion         460
#define cSetting_h_bond_exclusion           461
#define cSetting_label_shadow_mode          462
#define cSetting_light4                     463
#define cSetting_light5                     464
#define cSetting_light6                     465
#define cSetting_light7                     466
#define cSetting_label_outline_color        467
#define cSetting_ray_trace_mode             468
#define cSetting_ray_trace_gain             469
#define cSetting_selection_visible_only     470
#define cSetting_label_position             471
#define cSetting_ray_trace_depth_factor     472
#define cSetting_ray_trace_slope_factor     473
#define cSetting_ray_trace_disco_factor     474
#define cSetting_ray_shadow_decay_factor    475
#define cSetting_ray_interior_mode          476
#define cSetting_ray_legacy_lighting        477
#define cSetting_sculpt_auto_center         478
#define cSetting_pdb_discrete_chains        479
#define cSetting_pdb_unbond_cations         480
#define cSetting_sculpt_tri_scale           481
#define cSetting_sculpt_tri_weight          482
#define cSetting_sculpt_tri_min             483 
#define cSetting_sculpt_tri_max             484
#define cSetting_sculpt_tri_mode            485
#define cSetting_pdb_echo_tags              486
#define cSetting_connect_bonded             487
#define cSetting_spec_direct_power          488
#define cSetting_light8                     489
#define cSetting_light9                     490
#define cSetting_ray_shadow_decay_range     491
#define cSetting_spec_count                 492
#define cSetting_sculpt_min_scale           493
#define cSetting_sculpt_min_weight          494
#define cSetting_sculpt_min_min             495
#define cSetting_sculpt_min_max             496
#define cSetting_sculpt_max_scale           497
#define cSetting_sculpt_max_weight          498
#define cSetting_sculpt_max_min             499
#define cSetting_sculpt_max_max             500
#define cSetting_surface_circumscribe       501
#define cSetting_sculpt_avd_weight          502
#define cSetting_sculpt_avd_gap             503
#define cSetting_sculpt_avd_range           504
#define cSetting_sculpt_avd_excl            505
#define cSetting_async_builds               506
#define cSetting_fetch_path                 507
#define cSetting_cartoon_ring_radius        508
#define cSetting_ray_color_ramps            509
#define cSetting_ray_hint_camera            510
#define cSetting_ray_hint_shadow            511
#define cSetting_stick_valence_scale        512
#define cSetting_seq_view_alignment         513
#define cSetting_seq_view_unaligned_mode    514
#define cSetting_seq_view_unaligned_color   515
#define cSetting_seq_view_fill_char         516
#define cSetting_seq_view_fill_color        517
#define cSetting_seq_view_label_color       518
#define cSetting_surface_carve_normal_cutoff 519
#define cSetting_trace_atoms_mode           520
#define cSetting_session_changed            521
#define cSetting_ray_clip_shadows           522
#define cSetting_mouse_wheel_scale          523
#define cSetting_nonbonded_transparency     524
#define cSetting_ray_spec_local             525
#define cSetting_line_color                 526
#define cSetting_ray_label_specular         527
#define cSetting_mesh_skip                  528
#define cSetting_label_digits               529
#define cSetting_label_distance_digits      530
#define cSetting_label_angle_digits         531
#define cSetting_label_dihedral_digits      532
#define cSetting_surface_negative_visible   533
#define cSetting_surface_negative_color     534
#define cSetting_mesh_negative_visible      535
#define cSetting_mesh_negative_color        536
#define cSetting_group_auto_mode            537
#define cSetting_group_full_member_names    538
#define cSetting_gradient_max_length        539
#define cSetting_gradient_min_length        540
#define cSetting_gradient_min_slope         541
#define cSetting_gradient_normal_min_dot    542
#define cSetting_gradient_step_size         543
#define cSetting_gradient_spacing           544
#define cSetting_gradient_symmetry          545
#define cSetting_ray_trace_color            546
#define cSetting_group_arrow_prefix         547
#define cSetting_suppress_hidden            548
#define cSetting_session_compression        549
#define cSetting_movie_fps                  550
#define cSetting_ray_transparency_oblique   551
#define cSetting_ray_trace_trans_cutoff     552
#define cSetting_ray_trace_persist_cutoff   553
#define cSetting_ray_transparency_oblique_power 554
#define cSetting_ray_scatter                555
#define cSetting_h_bond_from_proton         556
#define cSetting_auto_copy_images           557
#define cSetting_moe_separate_chains        558
#define cSetting_transparency_global_sort   559
#define cSetting_hide_long_bonds            560
#define cSetting_auto_rename_duplicate_objects 561
#define cSetting_pdb_hetatm_guess_valences  562
#define cSetting_ellipsoid_quality          563
#define cSetting_cgo_ellipsoid_quality      564
#define cSetting_movie_animate_by_frame     565
#define cSetting_ramp_blend_nearby_colors   566
#define cSetting_auto_defer_builds          567
#define cSetting_ellipsoid_probability      568
#define cSetting_ellipsoid_scale            569
#define cSetting_ellipsoid_color            570
#define cSetting_ellipsoid_transparency     571
#define cSetting_movie_rock                 572
#define cSetting_cache_mode                 573
#define cSetting_dash_color                 574
#define cSetting_angle_color                575
#define cSetting_dihedral_color             576
#define cSetting_grid_mode                  577
#define cSetting_cache_max                  578
#define cSetting_grid_slot                  579
#define cSetting_grid_max                   580
#define cSetting_cartoon_putty_transform    581
#define cSetting_rock                       582
#define cSetting_cone_quality               583
#define cSetting_pdb_formal_charges         584
#define cSetting_ati_bugs                   585
#define cSetting_geometry_export_mode       586
#define cSetting_mouse_grid                 587
#define cSetting_mesh_cutoff                588
#define cSetting_mesh_carve_selection       589
#define cSetting_mesh_carve_state           590
#define cSetting_mesh_carve_cutoff          591
#define cSetting_mesh_clear_selection       592
#define cSetting_mesh_clear_state           593
#define cSetting_mesh_clear_cutoff          594
#define cSetting_mesh_grid_max              595
#define cSetting_session_cache_optimize     596
#define cSetting_sdof_drag_scale            597
#define cSetting_scene_buttons_mode         598
#define cSetting_scene_buttons              599
#define cSetting_map_auto_expand_sym        600
#define cSetting_image_copy_always          601
#define cSetting_max_ups                    602
#define cSetting_auto_overlay               603
#define cSetting_stick_ball_color           604
#define cSetting_stick_h_scale              605
#define cSetting_sculpt_pyra_inv_weight     606
#define cSetting_keep_alive                 607
#define cSetting_fit_kabsch                 608
#define cSetting_stereo_dynamic_strength    609
#define cSetting_dynamic_width              610
#define cSetting_dynamic_width_factor       611
#define cSetting_dynamic_width_min          612
#define cSetting_dynamic_width_max          613
#define cSetting_draw_mode                  614
#define cSetting_clean_electro_mode         615
#define cSetting_valence_mode               616
#define cSetting_show_frame_rate            617

/* when you add a new setting also remember:
   layer1/Setting.c
   modules/pymol/setting.py
  layer5/PyMOL.c 
*/

/* cSetting_ss_INIT must always be last setting_index +1 */

#define cSetting_INIT                       618

#endif


