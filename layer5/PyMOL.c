/* 
   A* -------------------------------------------------------------------
   B* This file contains source code for the PyMOL computer program
   C* copyright 1998-2000 by Warrn Lyford Delano of DeLano Scientific. 
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

#include "os_std.h"
#include "os_gl.h"

#include "MemoryDebug.h"

#include "Base.h"

#include "OVContext.h"

#include "MemoryDebug.h"
#include "MemoryCache.h"
#include "Err.h"
#include "Util.h"
#include "Selector.h"
#include "Color.h"
#include "Ortho.h"
#include "Scene.h"
#include "PyMOLObject.h"
#include "Executive.h"
#include "Word.h"
#include "RepMesh.h"
#include "ObjectMolecule.h"
#include "Control.h"
#include "Sphere.h"
#include "Setting.h"
#include "Ray.h"
#include "Util.h"
#include "Movie.h"
#include "P.h"
#include "Editor.h"
#include "SculptCache.h"
#include "Isosurf.h"
#include "Tetsurf.h"
#include "PConv.h"
#include "VFont.h"
#include "Wizard.h"
#include "Text.h"
#include "Character.h"
#include "Seq.h"
#include "Seeker.h"
#include "Texture.h"
#include "TestPyMOL.h"
#include "TypeFace.h"
#include "PlugIOManager.h"

#include "PyMOL.h"
#include "PyMOLGlobals.h"
#include "PyMOLOptions.h"

#ifndef _PYMOL_NOPY
PyMOLGlobals *SingletonPyMOLGlobals = NULL;
#endif

/* BEGIN PROPRIETARY CODE SEGMENT (see disclaimer in "os_proprietary.h") */ 
#ifndef _PYMOL_NOPY
#ifdef _MACPYMOL_XCODE
extern int *MacPyMOLReady;
extern CPyMOLOptions *MacPyMOLOption;
#endif
#endif
/* END PROPRIETARY CODE SEGMENT */

#ifdef _MACPYMOL_XCODE
/* BEGIN PROPRIETARY CODE SEGMENT (see disclaimer in "os_proprietary.h") */ 
#define PYMOL_API_LOCK if((I->PythonInitStage) && (!I->ModalDraw) && PLockAPIAsGlut(I->G,true)) {
#define PYMOL_API_LOCK_MODAL if((I->PythonInitStage) && PLockAPIAsGlut(I->G,true)) {
#define PYMOL_API_TRYLOCK PYMOL_API_LOCK
#define PYMOL_API_UNLOCK PUnlockAPIAsGlut(I->G); }
#define PYMOL_API_UNLOCK_NO_FLUSH PUnlockAPIAsGlutNoFlush(I->G); }
/* END PROPRIETARY CODE SEGMENT */
#else 
#ifdef _PYMOL_LIB_HAS_PYTHON
#define PYMOL_API_LOCK if(I->PythonInitStage && (!I->ModalDraw)) { PLockAPIAndUnblock(I->G); {
#define PYMOL_API_LOCK_MODAL if(I->PythonInitStage) { PLockAPIAndUnblock(I->G); {
#define PYMOL_API_TRYLOCK if(I->PythonInitStage && (!I->ModalDraw)) { if(PTryLockAPIAndUnblock(I->G)) {
#define PYMOL_API_UNLOCK PBlockAndUnlockAPI(I->G); }}
#define PYMOL_API_UNLOCK_NO_FLUSH PBlockAndUnlockAPI(I->G); }}
#else
#define PYMOL_API_LOCK if(!I->ModalDraw) {
#define PYMOL_API_LOCK_MODAL {
#define PYMOL_API_TRYLOCK if(!I->ModalDraw) {
#define PYMOL_API_UNLOCK }
#define PYMOL_API_UNLOCK_NO_FLUSH }
#endif
#endif
#define IDLE_AND_READY 3

typedef struct _CPyMOL {
  PyMOLGlobals *G;
  int FakeDragFlag;
  int RedisplayFlag;
  int PassiveFlag;
  int SwapFlag;
  int BusyFlag;
  int InterruptFlag; 
  int ReshapeFlag;
  int ClickReadyFlag;
  int DrawnFlag;
  ObjectNameType ClickedObject;
  int ClickedIndex, ClickedButton, ClickedModifiers, ClickedX, ClickedY;
  int ImageRequestedFlag,ImageReadyFlag;
  int DraggedFlag;
  int Reshape[PYMOL_RESHAPE_SIZE];
  int Progress[PYMOL_PROGRESS_SIZE];
  int ProgressChanged;
  int IdleAndReady;
  int ExpireCount;

  PyMOLModalDrawFn  *ModalDraw;

  PyMOLSwapBuffersFn *SwapFn;

/* Python stuff */
#ifndef _PYMOL_NOPY
  int PythonInitStage;
#endif
  /* dynamically mapped string constants */

  OVLexicon *Lex;
  ov_word lex_pdb, lex_mol2, lex_mol, lex_sdf, lex_xplor, lex_ccp4, lex_phi;
  ov_word lex_string, lex_filename, lex_raw;

  OVOneToOne *Rep;
  ov_word lex_everything, lex_sticks, lex_spheres, lex_surface;
  ov_word lex_labels, lex_nb_spheres, lex_cartoon, lex_ribbon;
  ov_word lex_lines, lex_mesh, lex_dots, lex_dashes, lex_nonbonded;
  ov_word lex_cell, lex_cgo, lex_callback, lex_extent, lex_slice;

  OVOneToOne *Clip;  
  ov_word lex_near, lex_far, lex_move, lex_slab, lex_atoms;

  OVOneToOne *Reinit;  
  ov_word lex_settings;

  OVOneToOne *SelectList;
  ov_word lex_index, lex_id, lex_rank;

  OVOneToOne *Setting;
  ov_word lex_bonding_vdw_cutoff;
  ov_word lex_min_mesh_spacing;
  ov_word lex_dot_density;
  ov_word lex_dot_mode;
  ov_word lex_solvent_radius;
  ov_word lex_sel_counter;
  ov_word lex_bg_rgb;
  ov_word lex_ambient;
  ov_word lex_direct;
  ov_word lex_reflect;
  ov_word lex_light;
  ov_word lex_power;
  ov_word lex_antialias;
  ov_word lex_cavity_cull;
  ov_word lex_gl_ambient;
  ov_word lex_single_image;
  ov_word lex_movie_delay;
  ov_word lex_ribbon_power;
  ov_word lex_ribbon_power_b;
  ov_word lex_ribbon_sampling;
  ov_word lex_ribbon_radius;
  ov_word lex_stick_radius;
  ov_word lex_hash_max;
  ov_word lex_ortho;
  ov_word lex_spec_reflect;
  ov_word lex_spec_power;
  ov_word lex_sweep_angle;
  ov_word lex_sweep_speed;
  ov_word lex_dot_hydrogens;
  ov_word lex_dot_radius;
  ov_word lex_ray_trace_frames;
  ov_word lex_cache_frames;
  ov_word lex_trim_dots;
  ov_word lex_cull_spheres;
  ov_word lex_test1;
  ov_word lex_test2;
  ov_word lex_surface_best;
  ov_word lex_surface_normal;
  ov_word lex_surface_quality;
  ov_word lex_surface_proximity;
  ov_word lex_normal_workaround;
  ov_word lex_stereo_angle;
  ov_word lex_stereo_shift;
  ov_word lex_line_smooth;
  ov_word lex_line_width;
  ov_word lex_half_bonds;
  ov_word lex_stick_quality;
  ov_word lex_stick_overlap;
  ov_word lex_stick_nub;
  ov_word lex_all_states;
  ov_word lex_pickable;
  ov_word lex_auto_show_lines;
  ov_word lex_idle_delay;
  ov_word lex_no_idle;
  ov_word lex_fast_idle;
  ov_word lex_slow_idle;
  ov_word lex_rock_delay;
  ov_word lex_dist_counter;
  ov_word lex_dash_length;
  ov_word lex_dash_gap;
  ov_word lex_auto_zoom;
  ov_word lex_overlay;
  ov_word lex_text;
  ov_word lex_button_mode;
  ov_word lex_valence;
  ov_word lex_nonbonded_size;
  ov_word lex_label_color;
  ov_word lex_ray_trace_fog;
  ov_word lex_spheroid_scale;
  ov_word lex_ray_trace_fog_start;
  ov_word lex_spheroid_smooth;
  ov_word lex_spheroid_fill;
  ov_word lex_auto_show_nonbonded;
  ov_word lex_cache_display;
  ov_word lex_mesh_radius;
  ov_word lex_backface_cull;
  ov_word lex_gamma;
  ov_word lex_dot_width;
  ov_word lex_auto_show_selections;
  ov_word lex_auto_hide_selections;
  ov_word lex_selection_width;
  ov_word lex_selection_overlay;
  ov_word lex_static_singletons;
  ov_word lex_max_triangles;
  ov_word lex_depth_cue;
  ov_word lex_specular;
  ov_word lex_shininess;
  ov_word lex_sphere_quality;
  ov_word lex_fog;
  ov_word lex_isomesh_auto_state;
  ov_word lex_mesh_width;
  ov_word lex_cartoon_sampling;
  ov_word lex_cartoon_loop_radius;
  ov_word lex_cartoon_loop_quality;
  ov_word lex_cartoon_power;
  ov_word lex_cartoon_power_b;
  ov_word lex_cartoon_rect_length;
  ov_word lex_cartoon_rect_width;
  ov_word lex_internal_gui_width;
  ov_word lex_internal_gui;
  ov_word lex_cartoon_oval_length;
  ov_word lex_cartoon_oval_width;
  ov_word lex_cartoon_oval_quality;
  ov_word lex_cartoon_tube_radius;
  ov_word lex_cartoon_tube_quality;
  ov_word lex_cartoon_debug;
  ov_word lex_ribbon_width;
  ov_word lex_dash_width;
  ov_word lex_dash_radius;
  ov_word lex_cgo_ray_width_scale;
  ov_word lex_line_radius;
  ov_word lex_cartoon_round_helices;
  ov_word lex_cartoon_refine_normals;
  ov_word lex_cartoon_flat_sheets;
  ov_word lex_cartoon_smooth_loops;
  ov_word lex_cartoon_dumbbell_length;
  ov_word lex_cartoon_dumbbell_width;
  ov_word lex_cartoon_dumbbell_radius;
  ov_word lex_cartoon_fancy_helices;
  ov_word lex_cartoon_fancy_sheets;
  ov_word lex_ignore_pdb_segi;
  ov_word lex_ribbon_throw;
  ov_word lex_cartoon_throw;
  ov_word lex_cartoon_refine;
  ov_word lex_cartoon_refine_tips;
  ov_word lex_cartoon_discrete_colors;
  ov_word lex_normalize_ccp4_maps;
  ov_word lex_surface_poor;
  ov_word lex_internal_feedback;
  ov_word lex_cgo_line_width;
  ov_word lex_cgo_line_radius;
  ov_word lex_logging;
  ov_word lex_robust_logs;
  ov_word lex_log_box_selections;
  ov_word lex_log_conformations;
  ov_word lex_valence_size;
  ov_word lex_surface_miserable;
  ov_word lex_ray_opaque_background;
  ov_word lex_transparency;
  ov_word lex_ray_texture;
  ov_word lex_ray_texture_settings;
  ov_word lex_suspend_updates;
  ov_word lex_full_screen;
  ov_word lex_surface_mode;
  ov_word lex_surface_color;
  ov_word lex_mesh_mode;
  ov_word lex_mesh_color;
  ov_word lex_auto_indicate_flags;
  ov_word lex_surface_debug;
  ov_word lex_ray_improve_shadows;
  ov_word lex_smooth_color_triangle;
  ov_word lex_ray_default_renderer;
  ov_word lex_field_of_view;
  ov_word lex_reflect_power;
  ov_word lex_preserve_chempy_ids;
  ov_word lex_sphere_scale;
  ov_word lex_two_sided_lighting;
  ov_word lex_secondary_structure;
  ov_word lex_auto_remove_hydrogens;
  ov_word lex_raise_exceptions;
  ov_word lex_stop_on_exceptions;
  ov_word lex_sculpting;
  ov_word lex_auto_sculpt;
  ov_word lex_sculpt_vdw_scale;
  ov_word lex_sculpt_vdw_scale14;
  ov_word lex_sculpt_vdw_weight;
  ov_word lex_sculpt_vdw_weight14;
  ov_word lex_sculpt_bond_weight;
  ov_word lex_sculpt_angl_weight;
  ov_word lex_sculpt_pyra_weight;
  ov_word lex_sculpt_plan_weight;
  ov_word lex_sculpting_cycles;
  ov_word lex_sphere_transparency;
  ov_word lex_sphere_color;
  ov_word lex_sculpt_field_mask;
  ov_word lex_sculpt_hb_overlap;
  ov_word lex_sculpt_hb_overlap_base;
  ov_word lex_legacy_vdw_radii;
  ov_word lex_sculpt_memory;
  ov_word lex_connect_mode;
  ov_word lex_cartoon_cylindrical_helices;
  ov_word lex_cartoon_helix_radius;
  ov_word lex_connect_cutoff;
  ov_word lex_save_pdb_ss;
  ov_word lex_sculpt_line_weight;
  ov_word lex_fit_iterations;
  ov_word lex_fit_tolerance;
  ov_word lex_batch_prefix;
  ov_word lex_stereo_mode;
  ov_word lex_cgo_sphere_quality;
  ov_word lex_pdb_literal_names;
  ov_word lex_wrap_output;
  ov_word lex_fog_start;
  ov_word lex_state;
  ov_word lex_frame;
  ov_word lex_ray_shadows;
  ov_word lex_ribbon_trace_atoms;
  ov_word lex_security;
  ov_word lex_stick_transparency; 
  ov_word lex_ray_transparency_shadows;
  ov_word lex_session_version_check;
  ov_word lex_ray_transparency_specular;
  ov_word lex_stereo_double_pump_mono;
  ov_word lex_sphere_solvent;
  ov_word lex_mesh_quality;
  ov_word lex_mesh_solvent;
  ov_word lex_dot_solvent;
  ov_word lex_ray_shadow_fudge;
  ov_word lex_ray_triangle_fudge;
  ov_word lex_debug_pick;
  ov_word lex_dot_color;
  ov_word lex_mouse_limit;
  ov_word lex_mouse_scale;
  ov_word lex_transparency_mode;
  ov_word lex_clamp_colors;
  ov_word lex_pymol_space_max_red;
  ov_word lex_pymol_space_max_green;
  ov_word lex_pymol_space_max_blue;
  ov_word lex_pymol_space_min_factor;
  ov_word lex_roving_origin;
  ov_word lex_roving_lines;
  ov_word lex_roving_sticks;
  ov_word lex_roving_spheres;
  ov_word lex_roving_labels;
  ov_word lex_roving_delay;
  ov_word lex_roving_selection;
  ov_word lex_roving_byres;
  ov_word lex_roving_ribbon;
  ov_word lex_roving_cartoon;
  ov_word lex_roving_polar_contacts;
  ov_word lex_roving_polar_cutoff;
  ov_word lex_roving_nonbonded;
  ov_word lex_float_labels;
  ov_word lex_roving_detail;
  ov_word lex_roving_nb_spheres;
  ov_word lex_ribbon_color;
  ov_word lex_cartoon_color;
  ov_word lex_ribbon_smooth;
  ov_word lex_auto_color;
  ov_word lex_auto_color_next;
  ov_word lex_ray_interior_color;
  ov_word lex_cartoon_highlight_color;
  ov_word lex_coulomb_units_factor;
  ov_word lex_coulomb_dielectric;
  ov_word lex_ray_interior_shadows;
  ov_word lex_ray_interior_texture;

  ov_word lex_roving_map1_name;
  ov_word lex_roving_map2_name;
  ov_word lex_roving_map3_name;

  ov_word lex_roving_map1_level;
  ov_word lex_roving_map2_level;
  ov_word lex_roving_map3_level;

  ov_word lex_roving_isomesh;
  ov_word lex_roving_isosurface;
  ov_word lex_scenes_changed;

  ov_word lex_gaussian_b_adjust;
  ov_word lex_pdb_standard_order;

  ov_word lex_cartoon_smooth_first;
  ov_word lex_cartoon_smooth_last;
  ov_word lex_cartoon_smooth_cycles;
  ov_word lex_cartoon_flat_cycles;

  ov_word lex_max_threads;
  ov_word lex_show_progress;
  ov_word lex_use_display_lists;
  ov_word lex_cache_memory;
  ov_word lex_simplify_display_lists;
  ov_word lex_retain_order;
  ov_word lex_pdb_hetatm_sort;
  ov_word lex_pdb_use_ter_records;
  ov_word lex_cartoon_trace_atoms;
  ov_word lex_ray_oversample_cutoff;
  ov_word lex_gaussian_resolution;
  ov_word lex_gaussian_b_floor;
  ov_word lex_sculpt_nb_interval;
  ov_word lex_sculpt_tors_weight;
  ov_word lex_sculpt_tors_tolerance;
  ov_word lex_stick_ball;
  ov_word lex_stick_ball_ratio;
  ov_word lex_stick_fixed_radius;
  ov_word lex_cartoon_transparency;
  ov_word lex_dash_round_ends;
  ov_word lex_h_bond_max_angle;
  ov_word lex_h_bond_cutoff_center;
  ov_word lex_h_bond_cutoff_edge;
  ov_word lex_h_bond_power_a;
  ov_word lex_h_bond_power_b;
  ov_word lex_h_bond_cone;

  ov_word lex_ss_helix_psi_target; 
  ov_word lex_ss_helix_psi_include;
  ov_word lex_ss_helix_psi_exclude;

  ov_word lex_ss_helix_phi_target;
  ov_word lex_ss_helix_phi_include;
  ov_word lex_ss_helix_phi_exclude;

  ov_word lex_ss_strand_psi_target;
  ov_word lex_ss_strand_psi_include;
  ov_word lex_ss_strand_psi_exclude;

  ov_word lex_ss_strand_phi_target;
  ov_word lex_ss_strand_phi_include;
  ov_word lex_ss_strand_phi_exclude;
  ov_word lex_movie_loop;

  ov_word lex_pdb_retain_ids;
  ov_word lex_pdb_no_end_record;
  ov_word lex_cgo_dot_width;
  ov_word lex_cgo_dot_radius;
  ov_word lex_defer_updates;
  ov_word lex_normalize_o_maps;
  ov_word lex_swap_dsn6_bytes;
  ov_word lex_pdb_insertions_go_first;
  ov_word lex_roving_origin_z;
  ov_word lex_roving_origin_z_cushion;
  ov_word lex_specular_intensity;
  ov_word lex_overlay_lines;
  ov_word lex_ray_transparency_spec_cut;
  ov_word lex_internal_prompt;
  ov_word lex_normalize_grd_maps;
  ov_word lex_ray_blend_colors;
  ov_word lex_ray_blend_red;
  ov_word lex_ray_blend_green;
  ov_word lex_ray_blend_blue;
  ov_word lex_png_screen_gamma;
  ov_word lex_png_file_gamma;
  ov_word lex_editor_label_fragments;
  ov_word lex_internal_gui_control_size;
  ov_word lex_auto_dss;
  ov_word lex_transparency_picking_mode;
  ov_word lex_virtual_trackball;
  ov_word lex_pdb_reformat_names_mode;
  ov_word lex_ray_pixel_scale;
  ov_word lex_label_font_id;
  ov_word lex_pdb_conect_all;
  ov_word lex_button_mode_name;
  ov_word lex_surface_type;
  ov_word lex_dot_normals;
  ov_word lex_session_migration;
  ov_word lex_mesh_normals;
  ov_word lex_mesh_type;
  ov_word lex_dot_lighting;
  ov_word lex_mesh_lighting;
  ov_word lex_surface_solvent; 
  ov_word lex_triangle_max_passes;
  ov_word lex_ray_interior_reflect;
  ov_word lex_internal_gui_mode;
  ov_word lex_surface_carve_selection;
  ov_word lex_surface_carve_state;
  ov_word lex_surface_carve_cutoff;
  ov_word lex_surface_clear_selection;
  ov_word lex_surface_clear_state;
  ov_word lex_surface_clear_cutoff;
  ov_word lex_surface_trim_cutoff;
  ov_word lex_surface_trim_factor;
  ov_word lex_ray_max_passes;
  ov_word lex_active_selections;
  ov_word lex_ray_transparency_contrast;
  ov_word lex_seq_view;
  ov_word lex_mouse_selection_mode;
  ov_word lex_seq_view_label_spacing;
  ov_word lex_seq_view_label_start;
  ov_word lex_seq_view_format;
  ov_word lex_seq_view_location;
  ov_word lex_seq_view_overlay;
  ov_word lex_auto_classify_atoms;
  ov_word lex_cartoon_nucleic_acid_mode;
  ov_word lex_seq_view_color;
  ov_word lex_seq_view_label_mode;
  ov_word lex_surface_ramp_above_mode;
  ov_word lex_stereo;
  ov_word lex_wizard_prompt_mode;
  ov_word lex_coulomb_cutoff;
  ov_word lex_slice_track_camera;
  ov_word lex_slice_height_scale; 
  ov_word lex_slice_height_map;
  ov_word lex_slice_grid;
  ov_word lex_slice_dynamic_grid;
  ov_word lex_slice_dynamic_grid_resolution;
  ov_word lex_pdb_insure_orthogonal;
  ov_word lex_ray_direct_shade;
  ov_word lex_stick_color;
  ov_word lex_cartoon_putty_radius;
  ov_word lex_cartoon_putty_quality;
  ov_word lex_cartoon_putty_scale_min;
  ov_word lex_cartoon_putty_scale_max;
  ov_word lex_cartoon_putty_scale_power;
  ov_word lex_cartoon_putty_range;
  ov_word lex_cartoon_side_chain_helper;
  ov_word lex_surface_optimize_subsets;
  ov_word lex_multiplex;
  ov_word lex_texture_fonts;
  ov_word lex_pqr_workarounds;
  ov_word lex_animation;
  ov_word lex_animation_duration;
  ov_word lex_scene_animation;
  ov_word lex_line_stick_helper;
  ov_word lex_ray_orthoscopic;
  ov_word lex_ribbon_side_chain_helper;
  ov_word lex_selection_width_max;
  ov_word lex_selection_width_scale;
  ov_word lex_scene_current_name;
  ov_word lex_presentation;
  ov_word lex_presentation_mode;
  ov_word lex_pdb_truncate_residue_name;
  ov_word lex_scene_loop;
  ov_word lex_sweep_mode;
  ov_word lex_sweep_phase;
  ov_word lex_scene_restart_movie_delay;
  ov_word lex_mouse_restart_movie_delay;
  ov_word lex_angle_size;
  ov_word lex_angle_label_position;
  ov_word lex_dihedral_size;
  ov_word lex_dihedral_label_position;
  ov_word lex_defer_builds_mode;
  ov_word lex_seq_view_discrete_by_state;
  ov_word lex_scene_animation_duration;
  ov_word lex_wildcard;
  ov_word lex_atom_name_wildcard;
  ov_word lex_ignore_case;
  ov_word lex_presentation_auto_quit;
  ov_word lex_editor_auto_dihedral;
  ov_word lex_presentation_auto_start;
  ov_word lex_validate_object_names;
  /*  ov_word lex_ray_pixel_scale_limit;*/
  ov_word lex_auto_show_spheres;
  ov_word lex_sphere_mode;
  ov_word lex_sphere_point_max_size;
  ov_word lex_sphere_point_size;
  ov_word lex_pdb_honor_model_number;
  ov_word lex_rank_assisted_sorts;
  ov_word lex_ribbon_nucleic_acid_mode;
  ov_word lex_cartoon_ring_mode;
  ov_word lex_cartoon_ring_width;
  ov_word lex_cartoon_ring_color;
  ov_word lex_cartoon_ring_finder;
  ov_word lex_cartoon_tube_cap;
  ov_word lex_cartoon_loop_cap;
  ov_word lex_nvidia_bugs;
  ov_word lex_image_dots_per_inch;
  ov_word lex_opaque_background;
  ov_word lex_draw_frames;
  ov_word lex_show_alpha_checker;
  ov_word lex_matrix_mode;
  ov_word lex_editor_auto_origin;
  ov_word lex_session_file;
  ov_word lex_cgo_transparency;
  ov_word lex_legacy_mouse_zoom;
  ov_word lex_auto_number_selections;
  ov_word lex_sculpt_vdw_vis_mode;
  ov_word lex_sculpt_vdw_vis_min;
  ov_word lex_sculpt_vdw_vis_mid;
  ov_word lex_sculpt_vdw_vis_max;
  ov_word lex_cartoon_ladder_mode;
  ov_word lex_cartoon_ladder_radius;
  ov_word lex_cartoon_ladder_color;
  ov_word lex_cartoon_nucleic_acid_color;
  ov_word lex_cartoon_ring_transparency;
  ov_word lex_label_size;
  ov_word lex_spec_direct;
  ov_word lex_light_count;
  ov_word lex_light2;
  ov_word lex_light3;
  ov_word lex_hide_underscore_names;
  ov_word lex_selection_round_points;
  ov_word lex_distance_exclusion;
  ov_word lex_h_bond_exclusion;
  ov_word lex_label_shadow_mode;
  ov_word lex_light4;
  ov_word lex_light5;
  ov_word lex_light6;
  ov_word lex_light7;
  ov_word lex_label_outline_color;
  ov_word lex_ray_trace_mode;
  ov_word lex_ray_trace_gain;
  ov_word lex_selection_visible_only;
  ov_word lex_label_position;
  ov_word lex_ray_trace_depth_factor;
  ov_word lex_ray_trace_slope_factor;
  ov_word lex_ray_trace_disco_factor;
  ov_word lex_ray_shadow_decay_factor;
  ov_word lex_ray_interior_mode;
  ov_word lex_ray_legacy_lighting;
  ov_word lex_sculpt_auto_center;
  ov_word lex_pdb_discrete_chains;
  ov_word lex_pdb_unbond_cations; 
  ov_word lex_sculpt_tri_scale;
  ov_word lex_sculpt_tri_weight;
  ov_word lex_sculpt_tri_min;
  ov_word lex_sculpt_tri_max;
  ov_word lex_sculpt_tri_mode;
  ov_word lex_pdb_echo_tags;
  ov_word lex_connect_bonded;
  ov_word lex_spec_direct_power;
  ov_word lex_light8;
  ov_word lex_light9;
  ov_word lex_ray_shadow_decay_range;
  ov_word lex_spec_count;
  ov_word lex_sculpt_min_scale;
  ov_word lex_sculpt_min_weight;
  ov_word lex_sculpt_min_min;
  ov_word lex_sculpt_min_max;
  ov_word lex_sculpt_max_scale;
  ov_word lex_sculpt_max_weight;
  ov_word lex_sculpt_max_min;
  ov_word lex_sculpt_max_max;
  ov_word lex_surface_circumscribe;
  ov_word lex_sculpt_avd_weight;          
  ov_word lex_sculpt_avd_gap;  
  ov_word lex_sculpt_avd_range;
  ov_word lex_sculpt_avd_excl;
  ov_word lex_async_builds;
  ov_word lex_fetch_path;
  ov_word lex_cartoon_ring_radius;
  ov_word lex_ray_color_ramps;
  ov_word lex_ray_hint_camera;
  ov_word lex_ray_hint_shadow;
  ov_word lex_stick_valence_scale;
  ov_word lex_seq_view_alignment;
  ov_word lex_seq_view_unaligned_mode;
  ov_word lex_seq_view_unaligned_color;
  ov_word lex_seq_view_fill_char;
  ov_word lex_seq_view_fill_color;
  ov_word lex_seq_view_label_color;
  ov_word lex_surface_carve_normal_cutoff;
  ov_word lex_trace_atoms_mode;
  ov_word lex_session_changed;
  ov_word lex_ray_clip_shadows;
  ov_word lex_mouse_wheel_scale;
  ov_word lex_nonbonded_transparency;
  ov_word lex_ray_spec_local;
  ov_word lex_line_color;
  ov_word lex_ray_label_specular;
  ov_word lex_mesh_skip;
  ov_word lex_label_digits;
  ov_word lex_label_distance_digits;
  ov_word lex_label_angle_digits;
  ov_word lex_label_dihedral_digits;
  ov_word lex_surface_negative_visible;
  ov_word lex_surface_negative_color;
  ov_word lex_mesh_negative_visible;
  ov_word lex_mesh_negative_color;
  ov_word lex_auto_group;
  ov_word lex_group_auto_mode;
  ov_word lex_group_full_member_names;
  ov_word lex_gradient_max_length;
  ov_word lex_gradient_min_length;
  ov_word lex_gradient_min_slope;
  ov_word lex_gradient_normal_min_dot;
  ov_word lex_gradient_step_size;
  ov_word lex_gradient_spacing;
  ov_word lex_gradient_symmetry;
  ov_word lex_ray_trace_color;
  ov_word lex_group_arrow_prefix;
  ov_word lex_suppress_hidden;
  ov_word lex_session_compression;
  ov_word lex_movie_fps;
  ov_word lex_ray_transparency_oblique;
  ov_word lex_ray_trace_trans_cutoff;
  ov_word lex_ray_trace_persist_cutoff; 
  ov_word lex_ray_transparency_oblique_power;
  ov_word lex_ray_scatter;
  ov_word lex_h_bond_from_proton;
  ov_word lex_auto_copy_images;
  ov_word lex_moe_separate_chains;
  ov_word lex_transparency_global_sort;
  ov_word lex_hide_long_bonds;
  ov_word lex_auto_rename_duplicate_objects; 
  ov_word lex_pdb_hetatm_guess_valences;
  ov_word lex_ellipsoid_quality;
  ov_word lex_cgo_ellipsoid_quality;
  ov_word lex_movie_animate_by_frame;
  ov_word lex_ramp_blend_nearby_colors;
  ov_word lex_auto_defer_builds;
  ov_word lex_ellipsoid_probability;
  ov_word lex_ellipsoid_scale;  
  ov_word lex_ellipsoid_color;       
  ov_word lex_ellipsoid_transparency;
  ov_word lex_movie_rock;
  ov_word lex_cache_mode;
  ov_word lex_dash_color;
  ov_word lex_angle_color;
  ov_word lex_dihedral_color;
  ov_word lex_grid_mode;
  ov_word lex_cache_max;
  ov_word lex_grid_slot;
  ov_word lex_grid_max;
  ov_word lex_cartoon_putty_transform;
  ov_word lex_rock;
  ov_word lex_cone_quality;
  ov_word lex_pdb_formal_charges;
  ov_word lex_ati_bugs;
  ov_word lex_geometry_export_mode;
  ov_word lex_mouse_grid;
  ov_word lex_mesh_cutoff;
  ov_word lex_mesh_carve_selection;
  ov_word lex_mesh_carve_state;
  ov_word lex_mesh_carve_cutoff;
  ov_word lex_mesh_clear_selection;
  ov_word lex_mesh_clear_state;
  ov_word lex_mesh_clear_cutoff;
  ov_word lex_mesh_grid_max;
  ov_word lex_session_cache_optimize;
  ov_word lex_sdof_drag_scale;
  ov_word lex_scene_buttons_mode;
  ov_word lex_scene_buttons;
  ov_word lex_map_auto_expand_sym;
  ov_word lex_image_copy_always;
  ov_word lex_max_ups;
  ov_word lex_auto_overlay;
  ov_word lex_stick_ball_color;
  ov_word lex_stick_h_scale;
  ov_word lex_sculpt_pyra_inv_weight;
  ov_word lex_keep_alive;
  ov_word lex_fit_kabsch;
  ov_word lex_stereo_dynamic_strength;
  ov_word lex_dynamic_width;
  ov_word lex_dynamic_width_factor;
  ov_word lex_dynamic_width_min;
  ov_word lex_dynamic_width_max;
  ov_word lex_draw_mode;
  ov_word lex_clean_electro_mode;
  ov_word lex_valence_mode;
  ov_word lex_show_frame_rate;
} _CPyMOL;

/* convenience functions -- inline */

#ifdef _PYMOL_INLINE
#define CC_INLINE __inline__
#else
#define CC_INLINE
#endif

CC_INLINE static PyMOLstatus get_status_ok(int ok) 
{
  if(ok) 
    return PyMOLstatus_SUCCESS;
  else
    return PyMOLstatus_FAILURE;
}

CC_INLINE static PyMOLreturn_status return_status_ok(int ok)
{
  PyMOLreturn_status result;
  result.status = get_status_ok(ok);
  return result;
}

CC_INLINE static PyMOLreturn_status return_status(int status)
{
  PyMOLreturn_status result;
  result.status = status;
  return result;
}

static OVstatus PyMOL_InitAPI(CPyMOL *I)
{
  OVContext *C = I->G->Context;
  OVreturn_word result;
  I->Lex = OVLexicon_New(C->heap);
  if(!I->Lex) 
    return_OVstatus_FAILURE;

  /* the following preprocessor macros may require GNU's cpp or VC++
     we'll see... */

#define LEX(ARG)  \
  if(!OVreturn_IS_OK( (result= OVLexicon_GetFromCString(I->Lex,#ARG))))  \
    return_OVstatus_FAILURE \
    else \
      I -> lex_ ## ARG = result.word;
  
  LEX(pdb);
  LEX(sdf);
  LEX(mol);
  LEX(mol2);
  LEX(xplor);
  LEX(ccp4);
  LEX(phi);

  LEX(string);
  LEX(filename);
  LEX(raw);

  /* string constants that are accepted on input */

#define LEX_REP(NAME,CODE) LEX(NAME) \
    if(!OVreturn_IS_OK( OVOneToOne_Set(I->Rep,I->lex_ ## NAME, CODE)))  \
      return_OVstatus_FAILURE;

  I->Rep = OVOneToOne_New(C->heap);
  if(!I->Rep)
    return_OVstatus_FAILURE;
  
  LEX_REP(everything,-1);
  LEX_REP(sticks,0);
  LEX_REP(spheres,1);
  LEX_REP(surface,2);
  LEX_REP(labels,3);
  LEX_REP(nb_spheres,4);
  LEX_REP(cartoon,5);
  LEX_REP(ribbon,6);
  LEX_REP(lines,7);
  LEX_REP(mesh,8);
  LEX_REP(dots,9);
  LEX_REP(dashes,10);
  LEX_REP(nonbonded,11);
  LEX_REP(cell,12);
  LEX_REP(cgo,13);
  LEX_REP(callback,14);
  LEX_REP(extent,15);
  LEX_REP(slice,16);
  
  /* workaround for unexplained bug with nested macro on VC6 */

#define LEX_CLIP(NAME,CODE) {if(!OVreturn_IS_OK( (result= OVLexicon_GetFromCString(I->Lex,#NAME))))  \
    return_OVstatus_FAILURE \
    else \
    I -> lex_ ## NAME = result.word;} \
    if(!OVreturn_IS_OK( OVOneToOne_Set(I->Clip,I->lex_ ## NAME, CODE)))  \
      return_OVstatus_FAILURE;
  
  I->Clip = OVOneToOne_New(C->heap);
  if(!I->Clip)
    return_OVstatus_FAILURE;

  LEX_CLIP(near,0);
  LEX_CLIP(far,1);
  LEX_CLIP(move,2);
  LEX_CLIP(slab,3);
  LEX_CLIP(atoms,4);

#define LEX_REINIT(NAME,CODE) {if(!OVreturn_IS_OK( (result= OVLexicon_GetFromCString(I->Lex,#NAME))))  \
    return_OVstatus_FAILURE \
    else \
    I -> lex_ ## NAME = result.word;} \
    if(!OVreturn_IS_OK( OVOneToOne_Set(I->Reinit,I->lex_ ## NAME, CODE)))  \
      return_OVstatus_FAILURE;
  
  I->Reinit = OVOneToOne_New(C->heap);
  if(!I->Reinit)
    return_OVstatus_FAILURE;

  if(!OVreturn_IS_OK( OVOneToOne_Set(I->Reinit,I->lex_everything, 0))) return_OVstatus_FAILURE;
  LEX_REINIT(settings,1);

#define LEX_SELLIST(NAME,CODE) {if(!OVreturn_IS_OK( (result= OVLexicon_GetFromCString(I->Lex,#NAME))))  \
    return_OVstatus_FAILURE \
    else \
    I -> lex_ ## NAME = result.word;} \
    if(!OVreturn_IS_OK( OVOneToOne_Set(I->SelectList,I->lex_ ## NAME, CODE)))  \
      return_OVstatus_FAILURE;
  
  I->SelectList = OVOneToOne_New(C->heap);
  if(!I->SelectList)
    return_OVstatus_FAILURE;

  LEX_SELLIST(index,0);
  LEX_SELLIST(id,1);
  LEX_SELLIST(rank,2);

  I->Setting = OVOneToOne_New(C->heap);
  if(!I->Setting)
    return_OVstatus_FAILURE;

#define LEX_SETTING(NAME,CODE) LEX(NAME) \
    if(!OVreturn_IS_OK( OVOneToOne_Set(I->Setting,I->lex_ ## NAME, CODE)))  \
      return_OVstatus_FAILURE;

  LEX_SETTING(bonding_vdw_cutoff, 0);
  LEX_SETTING(min_mesh_spacing, 1);
  LEX_SETTING(dot_density, 2);
  LEX_SETTING(dot_mode, 3);
  LEX_SETTING(solvent_radius, 4);
  LEX_SETTING(sel_counter, 5);
  LEX_SETTING(bg_rgb, 6);
  LEX_SETTING(ambient, 7);
  LEX_SETTING(direct, 8);
  LEX_SETTING(reflect, 9);
  LEX_SETTING(light, 10);
  LEX_SETTING(power, 11);
  LEX_SETTING(antialias, 12);
  LEX_SETTING(cavity_cull, 13);
  LEX_SETTING(gl_ambient, 14);
  LEX_SETTING(single_image, 15);
  LEX_SETTING(movie_delay, 16);
  LEX_SETTING(ribbon_power, 17);
  LEX_SETTING(ribbon_power_b, 18);
  LEX_SETTING(ribbon_sampling, 19);
  LEX_SETTING(ribbon_radius, 20);
  LEX_SETTING(stick_radius, 21);
  LEX_SETTING(hash_max, 22);
  LEX_SETTING(ortho, 23);
  LEX_SETTING(spec_reflect, 24);
  LEX_SETTING(spec_power, 25);
  LEX_SETTING(sweep_angle, 26);
  LEX_SETTING(sweep_speed, 27);
  LEX_SETTING(dot_hydrogens, 28);
  LEX_SETTING(dot_radius, 29);
  LEX_SETTING(ray_trace_frames, 30);
  LEX_SETTING(cache_frames, 31);
  LEX_SETTING(trim_dots, 32);
  LEX_SETTING(cull_spheres, 33);
  LEX_SETTING(test1, 34);
  LEX_SETTING(test2, 35);
  LEX_SETTING(surface_best, 36);
  LEX_SETTING(surface_normal, 37);
  LEX_SETTING(surface_quality, 38);
  LEX_SETTING(surface_proximity, 39);
  LEX_SETTING(normal_workaround, 40);
  LEX_SETTING(stereo_angle, 41);
  LEX_SETTING(stereo_shift, 42);
  LEX_SETTING(line_smooth, 43);
  LEX_SETTING(line_width, 44);
  LEX_SETTING(half_bonds, 45);
  LEX_SETTING(stick_quality, 46);
  LEX_SETTING(stick_overlap, 47);
  LEX_SETTING(stick_nub, 48);
  LEX_SETTING(all_states, 49);
  LEX_SETTING(pickable, 50);
  LEX_SETTING(auto_show_lines, 51);
  LEX_SETTING(idle_delay, 52);
  LEX_SETTING(no_idle, 53);
  LEX_SETTING(fast_idle, 54);
  LEX_SETTING(slow_idle, 55);
  LEX_SETTING(rock_delay, 56);
  LEX_SETTING(dist_counter, 57);
  LEX_SETTING(dash_length, 58);
  LEX_SETTING(dash_gap, 59);
  LEX_SETTING(auto_zoom, 60);
  LEX_SETTING(overlay, 61);
  LEX_SETTING(text, 62);
  LEX_SETTING(button_mode, 63);
  LEX_SETTING(valence, 64);
  LEX_SETTING(nonbonded_size, 65);
  LEX_SETTING(label_color, 66);
  LEX_SETTING(ray_trace_fog, 67);
  LEX_SETTING(spheroid_scale, 68);
  LEX_SETTING(ray_trace_fog_start, 69);
  LEX_SETTING(spheroid_smooth, 70);
  LEX_SETTING(spheroid_fill, 71);
  LEX_SETTING(auto_show_nonbonded, 72);
  LEX_SETTING(cache_display, 73);
  LEX_SETTING(mesh_radius, 74);
  LEX_SETTING(backface_cull, 75);
  LEX_SETTING(gamma, 76);
  LEX_SETTING(dot_width, 77);
  LEX_SETTING(auto_show_selections, 78);
  LEX_SETTING(auto_hide_selections, 79);
  LEX_SETTING(selection_width, 80);
  LEX_SETTING(selection_overlay, 81);
  LEX_SETTING(static_singletons, 82);
  LEX_SETTING(max_triangles, 83);
  LEX_SETTING(depth_cue, 84);
  LEX_SETTING(specular, 85);
  LEX_SETTING(shininess, 86);
  LEX_SETTING(sphere_quality, 87);
  LEX_SETTING(fog, 88);
  LEX_SETTING(isomesh_auto_state, 89);
  LEX_SETTING(mesh_width, 90);
  LEX_SETTING(cartoon_sampling, 91);
  LEX_SETTING(cartoon_loop_radius, 92);
  LEX_SETTING(cartoon_loop_quality, 93);
  LEX_SETTING(cartoon_power, 94);
  LEX_SETTING(cartoon_power_b, 95);
  LEX_SETTING(cartoon_rect_length, 96);
  LEX_SETTING(cartoon_rect_width, 97);
  LEX_SETTING(internal_gui_width, 98);
  LEX_SETTING(internal_gui, 99);
  LEX_SETTING(cartoon_oval_length, 100);
  LEX_SETTING(cartoon_oval_width, 101);
  LEX_SETTING(cartoon_oval_quality, 102);
  LEX_SETTING(cartoon_tube_radius, 103);
  LEX_SETTING(cartoon_tube_quality, 104);
  LEX_SETTING(cartoon_debug, 105);
  LEX_SETTING(ribbon_width, 106);
  LEX_SETTING(dash_width, 107);
  LEX_SETTING(dash_radius, 108);
  LEX_SETTING(cgo_ray_width_scale, 109);
  LEX_SETTING(line_radius, 110);
  LEX_SETTING(cartoon_round_helices, 111);
  LEX_SETTING(cartoon_refine_normals, 112);
  LEX_SETTING(cartoon_flat_sheets, 113);
  LEX_SETTING(cartoon_smooth_loops, 114);
  LEX_SETTING(cartoon_dumbbell_length, 115);
  LEX_SETTING(cartoon_dumbbell_width, 116);
  LEX_SETTING(cartoon_dumbbell_radius, 117);
  LEX_SETTING(cartoon_fancy_helices, 118);
  LEX_SETTING(cartoon_fancy_sheets, 119);
  LEX_SETTING(ignore_pdb_segi, 120);
  LEX_SETTING(ribbon_throw, 121);
  LEX_SETTING(cartoon_throw, 122);
  LEX_SETTING(cartoon_refine, 123);
  LEX_SETTING(cartoon_refine_tips, 124);
  LEX_SETTING(cartoon_discrete_colors, 125);
  LEX_SETTING(normalize_ccp4_maps, 126);
  LEX_SETTING(surface_poor, 127);
  LEX_SETTING(internal_feedback, 128);
  LEX_SETTING(cgo_line_width, 129);
  LEX_SETTING(cgo_line_radius, 130);
  LEX_SETTING(logging, 131);
  LEX_SETTING(robust_logs, 132);
  LEX_SETTING(log_box_selections, 133);
  LEX_SETTING(log_conformations, 134);
  LEX_SETTING(valence_size, 135);
  LEX_SETTING(surface_miserable, 136);
  LEX_SETTING(ray_opaque_background, 137);
  LEX_SETTING(transparency, 138);
  LEX_SETTING(ray_texture, 139);
  LEX_SETTING(ray_texture_settings, 140);
  LEX_SETTING(suspend_updates, 141);
  LEX_SETTING(full_screen, 142);
  LEX_SETTING(surface_mode, 143);
  LEX_SETTING(surface_color, 144);
  LEX_SETTING(mesh_mode, 145);
  LEX_SETTING(mesh_color, 146);
  LEX_SETTING(auto_indicate_flags, 147);
  LEX_SETTING(surface_debug, 148);
  LEX_SETTING(ray_improve_shadows, 149);
  LEX_SETTING(smooth_color_triangle, 150);
  LEX_SETTING(ray_default_renderer, 151);
  LEX_SETTING(field_of_view, 152);
  LEX_SETTING(reflect_power, 153);
  LEX_SETTING(preserve_chempy_ids, 154);
  LEX_SETTING(sphere_scale, 155);
  LEX_SETTING(two_sided_lighting, 156);
  LEX_SETTING(secondary_structure, 157);
  LEX_SETTING(auto_remove_hydrogens, 158);
  LEX_SETTING(raise_exceptions, 159);
  LEX_SETTING(stop_on_exceptions, 160);
  LEX_SETTING(sculpting, 161);
  LEX_SETTING(auto_sculpt, 162);
  LEX_SETTING(sculpt_vdw_scale, 163);
  LEX_SETTING(sculpt_vdw_scale14, 164);
  LEX_SETTING(sculpt_vdw_weight, 165);
  LEX_SETTING(sculpt_vdw_weight14, 166);
  LEX_SETTING(sculpt_bond_weight, 167);
  LEX_SETTING(sculpt_angl_weight, 168);
  LEX_SETTING(sculpt_pyra_weight, 169);
  LEX_SETTING(sculpt_plan_weight, 170);
  LEX_SETTING(sculpting_cycles, 171);
  LEX_SETTING(sphere_transparency, 172);
  LEX_SETTING(sphere_color, 173);
  LEX_SETTING(sculpt_field_mask, 174);
  LEX_SETTING(sculpt_hb_overlap, 175);
  LEX_SETTING(sculpt_hb_overlap_base, 176);
  LEX_SETTING(legacy_vdw_radii, 177);
  LEX_SETTING(sculpt_memory, 178);
  LEX_SETTING(connect_mode, 179);
  LEX_SETTING(cartoon_cylindrical_helices, 180);
  LEX_SETTING(cartoon_helix_radius, 181);
  LEX_SETTING(connect_cutoff, 182);
  LEX_SETTING(save_pdb_ss, 183);
  LEX_SETTING(sculpt_line_weight, 184);
  LEX_SETTING(fit_iterations, 185);
  LEX_SETTING(fit_tolerance, 186);
  LEX_SETTING(batch_prefix, 187);
  LEX_SETTING(stereo_mode, 188);
  LEX_SETTING(cgo_sphere_quality, 189);
  LEX_SETTING(pdb_literal_names, 190);
  LEX_SETTING(wrap_output, 191);
  LEX_SETTING(fog_start, 192);
  LEX_SETTING(state, 193);
  LEX_SETTING(frame, 194);
  LEX_SETTING(ray_shadows, 195);
  LEX_SETTING(ribbon_trace_atoms, 196);
  LEX_SETTING(security, 197);
  LEX_SETTING(stick_transparency, 198); 
  LEX_SETTING(ray_transparency_shadows, 199);
  LEX_SETTING(session_version_check, 200);
  LEX_SETTING(ray_transparency_specular, 201);
  LEX_SETTING(stereo_double_pump_mono, 202);
  LEX_SETTING(sphere_solvent, 203);
  LEX_SETTING(mesh_quality, 204);
  LEX_SETTING(mesh_solvent, 205);
  LEX_SETTING(dot_solvent, 206);
  LEX_SETTING(ray_shadow_fudge, 207);
  LEX_SETTING(ray_triangle_fudge, 208);
  LEX_SETTING(debug_pick, 209);
  LEX_SETTING(dot_color, 210);
  LEX_SETTING(mouse_limit, 211);
  LEX_SETTING(mouse_scale, 212);
  LEX_SETTING(transparency_mode, 213);
  LEX_SETTING(clamp_colors, 214);
  LEX_SETTING(pymol_space_max_red, 215);
  LEX_SETTING(pymol_space_max_green, 216);
  LEX_SETTING(pymol_space_max_blue, 217);
  LEX_SETTING(pymol_space_min_factor, 218);
  LEX_SETTING(roving_origin, 219);
  LEX_SETTING(roving_lines, 220);
  LEX_SETTING(roving_sticks, 221);
  LEX_SETTING(roving_spheres, 222);
  LEX_SETTING(roving_labels, 223);
  LEX_SETTING(roving_delay, 224);
  LEX_SETTING(roving_selection, 225);
  LEX_SETTING(roving_byres, 226);
  LEX_SETTING(roving_ribbon, 227);
  LEX_SETTING(roving_cartoon, 228);
  LEX_SETTING(roving_polar_contacts, 229);
  LEX_SETTING(roving_polar_cutoff, 230);
  LEX_SETTING(roving_nonbonded, 231);
  LEX_SETTING(float_labels, 232);
  LEX_SETTING(roving_detail, 233);
  LEX_SETTING(roving_nb_spheres, 234);
  LEX_SETTING(ribbon_color, 235);
  LEX_SETTING(cartoon_color, 236);
  LEX_SETTING(ribbon_smooth, 237);
  LEX_SETTING(auto_color, 238);
  LEX_SETTING(auto_color_next, 239);
  LEX_SETTING(ray_interior_color, 240);
  LEX_SETTING(cartoon_highlight_color, 241);
  LEX_SETTING(coulomb_units_factor, 242);
  LEX_SETTING(coulomb_dielectric, 243);
  LEX_SETTING(ray_interior_shadows, 244);
  LEX_SETTING(ray_interior_texture, 245);

  LEX_SETTING(roving_map1_name, 246);
  LEX_SETTING(roving_map2_name, 247);
  LEX_SETTING(roving_map3_name, 248);

  LEX_SETTING(roving_map1_level, 249);
  LEX_SETTING(roving_map2_level, 250);
  LEX_SETTING(roving_map3_level, 251);

  LEX_SETTING(roving_isomesh, 252);
  LEX_SETTING(roving_isosurface, 253);
  LEX_SETTING(scenes_changed, 254);

  LEX_SETTING(gaussian_b_adjust, 255);
  LEX_SETTING(pdb_standard_order, 256);

  LEX_SETTING(cartoon_smooth_first, 257);
  LEX_SETTING(cartoon_smooth_last, 258);
  LEX_SETTING(cartoon_smooth_cycles, 259);
  LEX_SETTING(cartoon_flat_cycles, 260);

  LEX_SETTING(max_threads, 261);
  LEX_SETTING(show_progress, 262);
  LEX_SETTING(use_display_lists, 263);
  LEX_SETTING(cache_memory, 264);
  LEX_SETTING(simplify_display_lists, 265);
  LEX_SETTING(retain_order, 266);
  LEX_SETTING(pdb_hetatm_sort, 267);
  LEX_SETTING(pdb_use_ter_records, 268);
  LEX_SETTING(cartoon_trace_atoms, 269);
  LEX_SETTING(ray_oversample_cutoff, 270);
  LEX_SETTING(gaussian_resolution, 271);
  LEX_SETTING(gaussian_b_floor, 272);
  LEX_SETTING(sculpt_nb_interval, 273);
  LEX_SETTING(sculpt_tors_weight, 274);
  LEX_SETTING(sculpt_tors_tolerance, 275);
  LEX_SETTING(stick_ball, 276);
  LEX_SETTING(stick_ball_ratio, 277);
  LEX_SETTING(stick_fixed_radius, 278);
  LEX_SETTING(cartoon_transparency, 279);
  LEX_SETTING(dash_round_ends, 280);
  LEX_SETTING(h_bond_max_angle, 281);
  LEX_SETTING(h_bond_cutoff_center, 282);
  LEX_SETTING(h_bond_cutoff_edge, 283);
  LEX_SETTING(h_bond_power_a, 284);
  LEX_SETTING(h_bond_power_b, 285);
  LEX_SETTING(h_bond_cone, 286);

  LEX_SETTING(ss_helix_psi_target, 287); 
  LEX_SETTING(ss_helix_psi_include, 288);
  LEX_SETTING(ss_helix_psi_exclude, 289);

  LEX_SETTING(ss_helix_phi_target, 290);
  LEX_SETTING(ss_helix_phi_include, 291);
  LEX_SETTING(ss_helix_phi_exclude, 292);

  LEX_SETTING(ss_strand_psi_target, 293);
  LEX_SETTING(ss_strand_psi_include, 294);
  LEX_SETTING(ss_strand_psi_exclude, 295);

  LEX_SETTING(ss_strand_phi_target, 296);
  LEX_SETTING(ss_strand_phi_include, 297);
  LEX_SETTING(ss_strand_phi_exclude, 298);
  LEX_SETTING(movie_loop, 299);

  LEX_SETTING(pdb_retain_ids, 300);
  LEX_SETTING(pdb_no_end_record, 301);
  LEX_SETTING(cgo_dot_width, 302);
  LEX_SETTING(cgo_dot_radius, 303);
  LEX_SETTING(defer_updates, 304);
  LEX_SETTING(normalize_o_maps, 305);
  LEX_SETTING(swap_dsn6_bytes, 306);
  LEX_SETTING(pdb_insertions_go_first, 307);
  LEX_SETTING(roving_origin_z, 308);
  LEX_SETTING(roving_origin_z_cushion, 309);
  LEX_SETTING(specular_intensity, 310);
  LEX_SETTING(overlay_lines, 311);
  LEX_SETTING(ray_transparency_spec_cut, 312);
  LEX_SETTING(internal_prompt, 313);
  LEX_SETTING(normalize_grd_maps, 314);
  LEX_SETTING(ray_blend_colors, 315);
  LEX_SETTING(ray_blend_red, 316);
  LEX_SETTING(ray_blend_green, 317);
  LEX_SETTING(ray_blend_blue, 318);
  LEX_SETTING(png_screen_gamma, 319);
  LEX_SETTING(png_file_gamma, 320);
  LEX_SETTING(editor_label_fragments, 321);
  LEX_SETTING(internal_gui_control_size, 322);
  LEX_SETTING(auto_dss, 323);
  LEX_SETTING(transparency_picking_mode, 324);
  LEX_SETTING(virtual_trackball, 325);
  LEX_SETTING(pdb_reformat_names_mode, 326);
  LEX_SETTING(ray_pixel_scale, 327);
  LEX_SETTING(label_font_id, 328);
  LEX_SETTING(pdb_conect_all, 329);
  LEX_SETTING(button_mode_name, 330);
  LEX_SETTING(surface_type, 331);
  LEX_SETTING(dot_normals, 332);
  LEX_SETTING(session_migration, 333);
  LEX_SETTING(mesh_normals, 334);
  LEX_SETTING(mesh_type, 335);
  LEX_SETTING(dot_lighting, 336);
  LEX_SETTING(mesh_lighting, 337);
  LEX_SETTING(surface_solvent, 338); 
  LEX_SETTING(triangle_max_passes, 339);
  LEX_SETTING(ray_interior_reflect, 340);
  LEX_SETTING(internal_gui_mode, 341);
  LEX_SETTING(surface_carve_selection, 342);
  LEX_SETTING(surface_carve_state, 343);
  LEX_SETTING(surface_carve_cutoff, 344);
  LEX_SETTING(surface_clear_selection, 345);
  LEX_SETTING(surface_clear_state, 346);
  LEX_SETTING(surface_clear_cutoff, 347);
  LEX_SETTING(surface_trim_cutoff, 348);
  LEX_SETTING(surface_trim_factor, 349);
  LEX_SETTING(ray_max_passes, 350);
  LEX_SETTING(active_selections, 351);
  LEX_SETTING(ray_transparency_contrast, 352);
  LEX_SETTING(seq_view, 353);
  LEX_SETTING(mouse_selection_mode, 354);
  LEX_SETTING(seq_view_label_spacing, 355);
  LEX_SETTING(seq_view_label_start, 356);
  LEX_SETTING(seq_view_format, 357);
  LEX_SETTING(seq_view_location, 358);
  LEX_SETTING(seq_view_overlay, 359);
  LEX_SETTING(auto_classify_atoms, 360);
  LEX_SETTING(cartoon_nucleic_acid_mode, 361);
  LEX_SETTING(seq_view_color, 362);
  LEX_SETTING(seq_view_label_mode, 363);
  LEX_SETTING(surface_ramp_above_mode, 364);
  LEX_SETTING(stereo, 365);
  LEX_SETTING(wizard_prompt_mode, 366);
  LEX_SETTING(coulomb_cutoff, 367);
  LEX_SETTING(slice_track_camera, 368);
  LEX_SETTING(slice_height_scale, 369); 
  LEX_SETTING(slice_height_map, 370);
  LEX_SETTING(slice_grid, 371);
  LEX_SETTING(slice_dynamic_grid, 372);
  LEX_SETTING(slice_dynamic_grid_resolution, 373);
  LEX_SETTING(pdb_insure_orthogonal, 374);
  LEX_SETTING(ray_direct_shade, 375);
  LEX_SETTING(stick_color, 376);
  LEX_SETTING(cartoon_putty_radius, 377);
  LEX_SETTING(cartoon_putty_quality, 378);
  LEX_SETTING(cartoon_putty_scale_min, 379);
  LEX_SETTING(cartoon_putty_scale_max, 380);
  LEX_SETTING(cartoon_putty_scale_power, 381);
  LEX_SETTING(cartoon_putty_range, 382);
  LEX_SETTING(cartoon_side_chain_helper, 383);
  LEX_SETTING(surface_optimize_subsets, 384);
  LEX_SETTING(multiplex, 385);
  LEX_SETTING(texture_fonts, 386);
  LEX_SETTING(pqr_workarounds, 387);
  LEX_SETTING(animation, 388);
  LEX_SETTING(animation_duration, 389);
  LEX_SETTING(scene_animation, 390);
  LEX_SETTING(line_stick_helper, 391);
  LEX_SETTING(ray_orthoscopic, 392);
  LEX_SETTING(ribbon_side_chain_helper, 393);
  LEX_SETTING(selection_width_max, 394);
  LEX_SETTING(selection_width_scale, 395);
  LEX_SETTING(scene_current_name, 396);
  LEX_SETTING(presentation, 397);
  LEX_SETTING(presentation_mode, 398);
  LEX_SETTING(pdb_truncate_residue_name, 399);
  LEX_SETTING(scene_loop, 400);
  LEX_SETTING(sweep_mode, 401);
  LEX_SETTING(sweep_phase, 402);
  LEX_SETTING(scene_restart_movie_delay, 403);
  LEX_SETTING(mouse_restart_movie_delay, 404);
  LEX_SETTING(angle_size, 405);
  LEX_SETTING(angle_label_position, 406);
  LEX_SETTING(dihedral_size, 407);
  LEX_SETTING(dihedral_label_position, 408);
  LEX_SETTING(defer_builds_mode, 409);
  LEX_SETTING(seq_view_discrete_by_state, 410);
  LEX_SETTING(scene_animation_duration, 411);
  LEX_SETTING(wildcard,412);
  LEX_SETTING(atom_name_wildcard,413);
  LEX_SETTING(ignore_case,414);
  LEX_SETTING(presentation_auto_quit,415);
  LEX_SETTING(editor_auto_dihedral,416);
  LEX_SETTING(presentation_auto_start,417);
  LEX_SETTING(validate_object_names,418);
  /*  LEX_SETTING(ray_pixel_scale_limit, 419);*/
  LEX_SETTING(auto_show_spheres, 420);
  LEX_SETTING(sphere_mode, 421);
  LEX_SETTING(sphere_point_max_size, 422);
  LEX_SETTING(sphere_point_size, 423);
  LEX_SETTING(pdb_honor_model_number, 424);
  LEX_SETTING(rank_assisted_sorts, 425);
  LEX_SETTING(ribbon_nucleic_acid_mode, 426);
  LEX_SETTING(cartoon_ring_mode, 427);
  LEX_SETTING(cartoon_ring_width, 428);
  LEX_SETTING(cartoon_ring_color, 429);
  LEX_SETTING(cartoon_ring_finder, 430);
  LEX_SETTING(cartoon_tube_cap, 431);
  LEX_SETTING(cartoon_loop_cap, 432);
  LEX_SETTING(nvidia_bugs, 433);
  LEX_SETTING(image_dots_per_inch, 434);
  LEX_SETTING(opaque_background, 435);
  LEX_SETTING(draw_frames, 436);
  LEX_SETTING(show_alpha_checker, 437);
  LEX_SETTING(matrix_mode, 438);
  LEX_SETTING(editor_auto_origin, 439);
  LEX_SETTING(session_file, 440);
  LEX_SETTING(cgo_transparency, 441);
  LEX_SETTING(legacy_mouse_zoom, 442);
  LEX_SETTING(auto_number_selections, 443);
  LEX_SETTING(sculpt_vdw_vis_mode, 444);
  LEX_SETTING(sculpt_vdw_vis_min, 445);
  LEX_SETTING(sculpt_vdw_vis_mid, 446);
  LEX_SETTING(sculpt_vdw_vis_max, 447);
  LEX_SETTING(cartoon_ladder_mode, 448);
  LEX_SETTING(cartoon_ladder_radius, 449);
  LEX_SETTING(cartoon_ladder_color, 450);
  LEX_SETTING(cartoon_nucleic_acid_color, 451);
  LEX_SETTING(cartoon_ring_transparency, 452);
  LEX_SETTING(label_size, 453);
  LEX_SETTING(spec_direct, 454);
  LEX_SETTING(light_count, 455);
  LEX_SETTING(light2, 456);
  LEX_SETTING(light3, 457);
  LEX_SETTING(hide_underscore_names, 458);
  LEX_SETTING(selection_round_points, 459);
  LEX_SETTING(distance_exclusion, 460);
  LEX_SETTING(h_bond_exclusion, 461);
  LEX_SETTING(label_shadow_mode, 462);
  LEX_SETTING(light4, 463);
  LEX_SETTING(light5, 464);
  LEX_SETTING(light6, 465);
  LEX_SETTING(light7, 466);
  LEX_SETTING(label_outline_color, 467);
  LEX_SETTING(ray_trace_mode,468);
  LEX_SETTING(ray_trace_gain,469);
  LEX_SETTING(selection_visible_only,470);
  LEX_SETTING(label_position,471);
  LEX_SETTING(ray_trace_depth_factor,472);
  LEX_SETTING(ray_trace_slope_factor,473);
  LEX_SETTING(ray_trace_disco_factor,474);
  LEX_SETTING(ray_shadow_decay_factor, 475);
  LEX_SETTING(ray_interior_mode, 476);
  LEX_SETTING(ray_legacy_lighting, 477);
  LEX_SETTING(sculpt_auto_center, 478);
  LEX_SETTING(pdb_discrete_chains, 479);
  LEX_SETTING(pdb_unbond_cations, 480);
  LEX_SETTING(sculpt_tri_scale, 481);
  LEX_SETTING(sculpt_tri_weight, 482);
  LEX_SETTING(sculpt_tri_min, 483);
  LEX_SETTING(sculpt_tri_max, 484);
  LEX_SETTING(sculpt_tri_mode, 485);
  LEX_SETTING(pdb_echo_tags, 486);
  LEX_SETTING(connect_bonded, 487);
  LEX_SETTING(spec_direct_power, 488);
  LEX_SETTING(light8, 489);
  LEX_SETTING(light9, 490);
  LEX_SETTING(ray_shadow_decay_range, 491);
  LEX_SETTING(spec_count, 492);
  LEX_SETTING(sculpt_min_scale, 493);
  LEX_SETTING(sculpt_min_weight, 494);
  LEX_SETTING(sculpt_min_min, 495);
  LEX_SETTING(sculpt_min_max, 496);
  LEX_SETTING(sculpt_max_scale, 497);
  LEX_SETTING(sculpt_max_weight, 498);
  LEX_SETTING(sculpt_max_min, 499);
  LEX_SETTING(sculpt_max_max, 500);
  LEX_SETTING(surface_circumscribe, 501);
  LEX_SETTING(sculpt_avd_weight, 502);
  LEX_SETTING(sculpt_avd_gap, 503);
  LEX_SETTING(sculpt_avd_range, 504);
  LEX_SETTING(sculpt_avd_excl, 505);
  LEX_SETTING(async_builds, 506);
  LEX_SETTING(fetch_path, 507);
  LEX_SETTING(cartoon_ring_radius, 508);
  LEX_SETTING(ray_color_ramps,509);
  LEX_SETTING(ray_hint_camera,510);
  LEX_SETTING(ray_hint_shadow,511);
  LEX_SETTING(stick_valence_scale,512);
  LEX_SETTING(seq_view_alignment,513);
  LEX_SETTING(seq_view_unaligned_mode,514);
  LEX_SETTING(seq_view_unaligned_color,515);
  LEX_SETTING(seq_view_fill_char,516);
  LEX_SETTING(seq_view_fill_color,517);
  LEX_SETTING(seq_view_label_color,518);
  LEX_SETTING(surface_carve_normal_cutoff,519);
  LEX_SETTING(trace_atoms_mode,520);
  LEX_SETTING(session_changed,521);
  LEX_SETTING(ray_clip_shadows,522);
  LEX_SETTING(mouse_wheel_scale,523);
  LEX_SETTING(nonbonded_transparency,524);
  LEX_SETTING(ray_spec_local,525);
  LEX_SETTING(line_color,526);
  LEX_SETTING(ray_label_specular, 527);
  LEX_SETTING(mesh_skip,528);
  LEX_SETTING(label_digits,529);
  LEX_SETTING(label_distance_digits,530);
  LEX_SETTING(label_angle_digits,531);
  LEX_SETTING(label_dihedral_digits,532);
  LEX_SETTING(surface_negative_visible,533);
  LEX_SETTING(surface_negative_color, 534);
  LEX_SETTING(mesh_negative_visible,535);
  LEX_SETTING(mesh_negative_color, 536);
  LEX_SETTING(group_auto_mode, 537);
  LEX_SETTING(group_full_member_names, 538);
  LEX_SETTING(gradient_max_length,539);
  LEX_SETTING(gradient_min_length,540);
  LEX_SETTING(gradient_min_slope,541);
  LEX_SETTING(gradient_normal_min_dot,542);
  LEX_SETTING(gradient_step_size,543);
  LEX_SETTING(gradient_spacing,544);
  LEX_SETTING(gradient_symmetry,545);
  LEX_SETTING(ray_trace_color,546);
  LEX_SETTING(group_arrow_prefix,547);
  LEX_SETTING(suppress_hidden,548);
  LEX_SETTING(session_compression,549);
  LEX_SETTING(movie_fps,550);
  LEX_SETTING(ray_transparency_oblique,551);
  LEX_SETTING(ray_trace_trans_cutoff,552);
  LEX_SETTING(ray_trace_persist_cutoff,553);
  LEX_SETTING(ray_transparency_oblique_power,554);  
  LEX_SETTING(ray_scatter,555);
  LEX_SETTING(h_bond_from_proton,556);
  LEX_SETTING(auto_copy_images, 557);
  LEX_SETTING(moe_separate_chains, 558);
  LEX_SETTING(transparency_global_sort, 559);
  LEX_SETTING(hide_long_bonds, 560);
  LEX_SETTING(auto_rename_duplicate_objects, 561);
  LEX_SETTING(pdb_hetatm_guess_valences, 562);
  LEX_SETTING(ellipsoid_quality, 563);
  LEX_SETTING(cgo_ellipsoid_quality, 564);
  LEX_SETTING(movie_animate_by_frame, 565);
  LEX_SETTING(ramp_blend_nearby_colors, 566);
  LEX_SETTING(auto_defer_builds, 567);
  LEX_SETTING(ellipsoid_probability,568);
  LEX_SETTING(ellipsoid_scale,569);
  LEX_SETTING(ellipsoid_color,570);
  LEX_SETTING(ellipsoid_transparency,571);
  LEX_SETTING(movie_rock,572);
  LEX_SETTING(cache_mode,573);
  LEX_SETTING(dash_color,574);
  LEX_SETTING(angle_color,575);
  LEX_SETTING(dihedral_color,576);
  LEX_SETTING(grid_mode,577);
  LEX_SETTING(cache_max,578);
  LEX_SETTING(grid_slot,579);
  LEX_SETTING(grid_max,580);
  LEX_SETTING(cartoon_putty_transform,581);
  LEX_SETTING(rock,582);
  LEX_SETTING(cone_quality,583);
  LEX_SETTING(pdb_formal_charges,584);
  LEX_SETTING(ati_bugs,585);
  LEX_SETTING(geometry_export_mode,586);
  LEX_SETTING(mouse_grid,587);
  LEX_SETTING(mesh_cutoff,588);
  LEX_SETTING(mesh_carve_selection,589);
  LEX_SETTING(mesh_carve_state,590);
  LEX_SETTING(mesh_carve_cutoff,591);
  LEX_SETTING(mesh_clear_selection,592);
  LEX_SETTING(mesh_clear_state,593);
  LEX_SETTING(mesh_clear_cutoff,594);
  LEX_SETTING(mesh_grid_max,595);
  LEX_SETTING(session_cache_optimize,596);
  LEX_SETTING(sdof_drag_scale,597);
  LEX_SETTING(scene_buttons_mode,598);
  LEX_SETTING(scene_buttons,599);
  LEX_SETTING(map_auto_expand_sym,600);
  LEX_SETTING(image_copy_always,601);
  LEX_SETTING(max_ups,602);
  LEX_SETTING(auto_overlay,603);
  LEX_SETTING(stick_ball_color,604);
  LEX_SETTING(stick_h_scale,605);
  LEX_SETTING(sculpt_pyra_inv_weight,606);
  LEX_SETTING(keep_alive,607);
  LEX_SETTING(fit_kabsch,608);
  LEX_SETTING(stereo_dynamic_strength,609);
  LEX_SETTING(dynamic_width,610);
  LEX_SETTING(dynamic_width_factor,611);
  LEX_SETTING(dynamic_width_min,612);
  LEX_SETTING(dynamic_width_max,613);
  LEX_SETTING(draw_mode,614);
  LEX_SETTING(clean_electro_mode,615);
  LEX_SETTING(valence_mode,616);
  LEX_SETTING(show_frame_rate,617);
  return_OVstatus_SUCCESS;
}

int PyMOL_NewG3DStream(CPyMOL *I,int **array_ptr)
{
  int *return_vla = ExecutiveGetG3d(I->G);
  int result = OVstatus_FAILURE;
  if(return_vla) {
    result = VLAGetSize(return_vla)*(sizeof(G3dPrimitive)/sizeof(int));
  }
  if(array_ptr)
    *array_ptr = return_vla;
  return result;
}

int PyMOL_DelG3DStream(CPyMOL *I,int *array_ptr)
{
  VLAFreeP(array_ptr);
  return OVstatus_SUCCESS;
}


static OVstatus PyMOL_PurgeAPI(CPyMOL *I)
{
  OVOneToOne_DEL_AUTO_NULL(I->Setting);
  OVOneToOne_DEL_AUTO_NULL(I->Clip);
  OVOneToOne_DEL_AUTO_NULL(I->SelectList);
  OVOneToOne_DEL_AUTO_NULL(I->Reinit);
  OVOneToOne_DEL_AUTO_NULL(I->Rep);
  OVLexicon_DEL_AUTO_NULL(I->Lex);
  return_OVstatus_SUCCESS;
}

int PyMOL_FreeResultArray(CPyMOL *I,void *array)
{
  if(array) {
    VLAFreeP(array);
    return PyMOLstatus_SUCCESS;
  } else {
    return PyMOLstatus_FAILURE; 
  }
}

PyMOLreturn_status PyMOL_CmdDraw(CPyMOL *I,int width, int height,
                                 int antialias, int quiet)
{
  PyMOLreturn_status result = { PyMOLstatus_FAILURE };
  PYMOL_API_LOCK
  result.status = get_status_ok(ExecutiveDrawCmd(I->G,width,height,antialias,false,quiet));
  I->ImageRequestedFlag = true;
  I->ImageReadyFlag = false;
  PYMOL_API_UNLOCK
  return result;
}

PyMOLreturn_status PyMOL_CmdCapture(CPyMOL *I,int quiet)
{
  PyMOLreturn_status result = { PyMOLstatus_FAILURE };
  PYMOL_API_LOCK
  result.status = get_status_ok(ExecutiveDrawCmd(I->G,-1,-1,0,true,quiet));
  I->ImageRequestedFlag = true;
  I->ImageReadyFlag = false;
  PYMOL_API_UNLOCK
  return result;
}

PyMOLreturn_status PyMOL_CmdRay(CPyMOL *I,int width, int height,int antialias,
                                float angle, float shift,int renderer, int defer, int quiet)
{
  PyMOLreturn_status result = { PyMOLstatus_FAILURE };
  PYMOL_API_LOCK
  
  if(renderer<0) renderer=SettingGetGlobal_i(I->G,cSetting_ray_default_renderer);
  SceneInvalidateCopy(I->G,true);
  result.status = get_status_ok(ExecutiveRay(I->G,width,height,renderer,angle,shift,quiet,defer,antialias));
  if(defer) {
    I->ImageRequestedFlag = true;
    I->ImageReadyFlag = false;
  } else {
    I->ImageRequestedFlag = false;
    if(SceneHasImage(I->G)) {
      I->ImageReadyFlag = true;
    } else {
      I->ImageReadyFlag = false;
    }
  }
  PYMOL_API_UNLOCK
  return result;
}

PyMOLreturn_status PyMOL_CmdSetView(CPyMOL *I, float *view, int view_len,
                                    float animate, int quiet)
{
  PyMOLreturn_status result = { PyMOLstatus_FAILURE };
  SceneViewType tmp;
  PYMOL_API_LOCK
  if(view_len>=18) {
    int a;
    UtilZeroMem(tmp,sizeof(tmp));
    tmp[15]=1.0F;
    for(a=0;a<3;a++) {
      tmp[a]=view[a];
      tmp[a+4]=view[a+3];
      tmp[a+8]=view[a+6];
      tmp[a+16]=view[a+9];
      tmp[a+19]=view[a+12];
      tmp[a+22]=view[a+15];
    }
    SceneSetView(I->G,tmp,quiet,animate,0);  /* TO DO -- add hand to the API */
    result.status = get_status_ok(true);
  } else {
    result.status = get_status_ok(false);
  }
  PYMOL_API_UNLOCK
  return result;
}

PyMOLreturn_float_array PyMOL_CmdGetView(CPyMOL *I,int quiet)
{
  PyMOLreturn_float_array result = { PyMOLstatus_FAILURE };
  SceneViewType tmp;
  PYMOL_API_LOCK
  result.size = 18;
  result.array = VLAlloc(float,result.size);
  if(result.array) {
    int a;
    SceneGetView(I->G,tmp);
    for(a=0;a<3;a++) {
      result.array[a]=tmp[a];
      result.array[a+3]=tmp[a+4];
      result.array[a+6]=tmp[a+8];
      result.array[a+9]=tmp[a+16];
      result.array[a+12]=tmp[a+19];
      result.array[a+15]=tmp[a+22];
    }
    result.status = get_status_ok(true);
  } else {
    result.status = get_status_ok(false);
  }
  PYMOL_API_UNLOCK
  return result;
}

PyMOLreturn_float_array PyMOL_CmdAlign(CPyMOL *I, char *source, char *target, float cutoff, 
                                 int cycles, float gap, float extend, int max_gap, 
                                 char *object, char *matrix, int source_state, int target_state, 
                                 int quiet, int max_skip, int transform, int reset) 
{
  PyMOLreturn_float_array result = { PyMOLstatus_FAILURE };
  PYMOL_API_LOCK
  OrthoLineType s2="",s3="";
  int ok = false;
  ExecutiveRMSInfo rms_info;
  result.size = 7;
  result.array = VLAlloc(float,result.size);
  if(!result.array) {
    ok=false;
  } else {
    ok = ((SelectorGetTmp(I->G,source,s2)>=0) &&
          (SelectorGetTmp(I->G,target,s3)>=0));
    if(ok) {
      const float _0 = 0.0F; /* GCC compiler bug workaround */
      const float _m1 = -1.0F;
      ok = ExecutiveAlign(I->G,s2,s3,matrix,gap,extend,max_gap,
                          max_skip,cutoff,cycles,quiet,object,
                          source_state-1, target_state-1,
                          &rms_info,transform,reset,
                          _m1, _0, _0, _0, _0, _0, 0, _0);
      if(ok) {
        result.array[0] = rms_info.final_rms;
        result.array[1] = rms_info.final_n_atom;
        result.array[2] = rms_info.n_cycles_run;
        result.array[3] = rms_info.initial_rms;
        result.array[4] = rms_info.initial_n_atom;
        result.array[5] = rms_info.raw_alignment_score;
        result.array[6] = rms_info.n_residues_aligned;
      }
    }
  }
  SelectorFreeTmp(I->G,s2);
  SelectorFreeTmp(I->G,s3);
  if(!ok) { 
    VLAFreeP(result.array);
  }
  result.status = get_status_ok(ok);

  PYMOL_API_UNLOCK
  return result;
 
}

PyMOLreturn_status PyMOL_CmdDelete(CPyMOL *I,char *name,int quiet)
{
  PYMOL_API_LOCK
  ExecutiveDelete(I->G, name);
  PYMOL_API_UNLOCK
 return return_status_ok(true); /* TO DO: return a real result */
}

PyMOLreturn_status PyMOL_CmdZoom(CPyMOL *I,char *selection, float buffer,
                              int state, int complete, float animate, int quiet)
{
  int ok=false;
  PYMOL_API_LOCK
  ok = ExecutiveWindowZoom(I->G, selection, buffer, state-1, 
                               complete, animate, quiet);
  PYMOL_API_UNLOCK
  return return_status_ok(ok);
}

PyMOLreturn_status PyMOL_CmdOrient(CPyMOL *I,char *selection, float buffer, 
                                int state, int complete, float  animate, int quiet)
{
  int ok=true;
  PYMOL_API_LOCK
  double m[16];
  OrthoLineType s1;
  SelectorGetTmp(I->G,selection,s1);
  if(ExecutiveGetMoment(I->G,s1,m,state))
    ExecutiveOrient(I->G,s1,m,state-1,animate,complete,buffer,quiet); /* TODO STATUS */
  else
    ok=false;
  SelectorFreeTmp(I->G,s1);
  PYMOL_API_UNLOCK
  return return_status_ok(ok);
}

PyMOLreturn_status PyMOL_CmdCenter(CPyMOL *I,char *selection, int state, int origin, float animate, int quiet)
{
  int ok = false;
  PYMOL_API_LOCK
  ok = ExecutiveCenter(I->G,selection,state-1,origin,animate,NULL,quiet);
  PYMOL_API_UNLOCK
  return return_status_ok(ok);
}


PyMOLreturn_status PyMOL_CmdOrigin(CPyMOL *I,char *selection, int state, int quiet)
{
  int ok=true;
  PYMOL_API_LOCK
  OrthoLineType s1;
  float v[3] = { 0.0F, 0.0F, 0.0F };
  SelectorGetTmp(I->G,selection,s1);
  ok = ExecutiveOrigin(I->G,s1,true,"",v,state-1); /* TODO STATUS */
  SelectorFreeTmp(I->G,s1);
  PYMOL_API_UNLOCK
  return return_status_ok(ok);
}

PyMOLreturn_status PyMOL_CmdOriginAt(CPyMOL *I,float x, float y, float z, int quiet)
{
  int ok=true;
  PYMOL_API_LOCK
  float v[3];
  v[0]=x;v[1]=y;v[2]=z;
  ok = ExecutiveOrigin(I->G,"",true,"",v,0); /* TODO STATUS */
  PYMOL_API_UNLOCK
  return return_status_ok(ok);
}

static OVreturn_word get_rep_id(CPyMOL *I,char *representation)
{
  OVreturn_word result;

  if(!OVreturn_IS_OK( (result = OVLexicon_BorrowFromCString(I->Lex,representation))))
    return result;
  return OVOneToOne_GetForward(I->Rep,result.word);
}

static OVreturn_word get_setting_id(CPyMOL *I,char *setting)
{
  OVreturn_word result;
  if(!OVreturn_IS_OK( (result = OVLexicon_BorrowFromCString(I->Lex,setting))))
    return result;
  return OVOneToOne_GetForward(I->Setting,result.word);
}

static OVreturn_word get_clip_id(CPyMOL *I,char *clip)
{
  OVreturn_word result;
  if(!OVreturn_IS_OK( (result = OVLexicon_BorrowFromCString(I->Lex,clip))))
    return result;
  return OVOneToOne_GetForward(I->Clip,result.word);
}

static OVreturn_word get_reinit_id(CPyMOL *I,char *reinit)
{
  OVreturn_word result;
  if(!OVreturn_IS_OK( (result = OVLexicon_BorrowFromCString(I->Lex,reinit))))
    return result;
  return OVOneToOne_GetForward(I->Reinit,result.word);
}

static OVreturn_word get_select_list_mode(CPyMOL *I,char *mode)
{
  OVreturn_word result;
  if(!OVreturn_IS_OK( (result = OVLexicon_BorrowFromCString(I->Lex,mode))))
    return result;
  return OVOneToOne_GetForward(I->SelectList,result.word);
}

PyMOLreturn_status PyMOL_CmdClip(CPyMOL *I,char *mode, float amount, char *selection, int state, int quiet)
{
  int ok=true;
  PYMOL_API_LOCK
  OrthoLineType s1;
  OVreturn_word clip_id;
  if(OVreturn_IS_OK( (clip_id= get_clip_id(I,mode)))) {
    SelectorGetTmp(I->G,selection,s1);
    SceneClip(I->G,clip_id.word,amount,s1,state-1);
    SelectorFreeTmp(I->G,s1);
  }
  PYMOL_API_UNLOCK
  return return_status_ok(ok);
}

PyMOLreturn_status PyMOL_CmdLabel(CPyMOL *I,char *selection, char *text, int quiet)
{
  int ok=true;
  PYMOL_API_LOCK
  OrthoLineType s1;
  SelectorGetTmp(I->G,selection,s1);
  ok = ExecutiveLabel(I->G,s1,text,quiet,0);
  SelectorFreeTmp(I->G,s1);
  PYMOL_API_UNLOCK
  return return_status_ok(ok);
}

PyMOLreturn_status PyMOL_CmdSelect(CPyMOL *I,char *name, char *selection, int quiet)
{
  int ok=true;
  PYMOL_API_LOCK
  ok = SelectorCreate(I->G,name,selection,NULL,quiet,NULL);
  PYMOL_API_UNLOCK
  return return_status_ok(ok);
}

PyMOLreturn_status PyMOL_CmdSelectList(CPyMOL *I,char *name, char *object, int *list, int list_len, int state, char *mode, int quiet)
{
  PyMOLreturn_status result = { PyMOLstatus_FAILURE };
  PYMOL_API_LOCK
  OVreturn_word mode_id;
  if(OVreturn_IS_OK( (mode_id= get_select_list_mode(I,mode)))) {
    result.status = ExecutiveSelectList(I->G,name, object,  list, list_len, state-1, mode_id.word, quiet);
  }
  PYMOL_API_UNLOCK
  return result;
}

PyMOLreturn_status PyMOL_CmdShow(CPyMOL *I,char *representation, char *selection, int quiet)
{
  int ok=true;
  PYMOL_API_LOCK
  OrthoLineType s1;
  OVreturn_word rep_id;
  if(OVreturn_IS_OK( (rep_id= get_rep_id(I,representation)))) {
    SelectorGetTmp(I->G,selection,s1);
    ExecutiveSetRepVisib(I->G,s1,rep_id.word,true);
    SelectorFreeTmp(I->G,s1);
  } else {
    ok=false;
  }
  PYMOL_API_UNLOCK
  return return_status_ok(ok);
}

PyMOLreturn_status PyMOL_CmdHide(CPyMOL *I,char *representation, char *selection, int quiet)
{
  int ok=true;
  PYMOL_API_LOCK
  OrthoLineType s1;
  OVreturn_word rep_id;
  if(OVreturn_IS_OK( (rep_id = get_rep_id(I,representation)))) {
    SelectorGetTmp(I->G,selection,s1);
    ExecutiveSetRepVisib(I->G,s1,rep_id.word,false);
    SelectorFreeTmp(I->G,s1);
  } else {
    ok=false;
  }
  PYMOL_API_UNLOCK
  return return_status_ok(ok);
}

PyMOLreturn_status PyMOL_CmdEnable(CPyMOL *I,char *name,int quiet)
{
  int ok = false; 
  PYMOL_API_LOCK
  if(name[0]=='(') {
    OrthoLineType s1;
    ok = (SelectorGetTmp(I->G,name,s1)>=0);
    if(ok) ok = ExecutiveSetOnOffBySele(I->G,s1,true);
    SelectorFreeTmp(I->G,s1);
  }
  ok = ExecutiveSetObjVisib(I->G,name,true,false); /* TO DO: parents */
  PYMOL_API_UNLOCK
  return return_status_ok(ok);
}

PyMOLreturn_status PyMOL_CmdDisable(CPyMOL *I,char *name,int quiet)
{
  int ok = false;
  PYMOL_API_LOCK
  if(name[0]=='(') {
    OrthoLineType s1;
    ok = (SelectorGetTmp(I->G,name,s1)>=0);
    if(ok) ok = ExecutiveSetOnOffBySele(I->G,s1,false);
    SelectorFreeTmp(I->G,s1);
  } else {
    ok = ExecutiveSetObjVisib(I->G,name,false,false);
  }
  PYMOL_API_UNLOCK
  return return_status_ok(ok);
}

PyMOLreturn_status PyMOL_CmdSet(CPyMOL *I,char *setting, char *value, char *selection, int state, int quiet, int side_effects)
{
  int ok=true;
  PYMOL_API_LOCK
  OVreturn_word setting_id;
  if(OVreturn_IS_OK( (setting_id = get_setting_id(I,setting)))) {
    ExecutiveSetSettingFromString(I->G, setting_id.word, value, selection,
                                  state-1, quiet, side_effects);
  }
  PYMOL_API_UNLOCK
  return return_status_ok(ok);
}

PyMOLreturn_status PyMOL_CmdColor(CPyMOL *I,char *color, char *selection, int flags, int quiet)
{
  int ok=true;
  PYMOL_API_LOCK
  OrthoLineType s1 = "";
  
  SelectorGetTmp(I->G,selection,s1);
  ok = ExecutiveColor(I->G,s1,color,flags,quiet);
  SelectorFreeTmp(I->G,s1);
  PYMOL_API_UNLOCK
  return return_status_ok(ok);
}

PyMOLreturn_status PyMOL_CmdReinitialize(CPyMOL *I,char *what, char *object_name)
{
  int ok = true;
  OVreturn_word what_id;
  PYMOL_API_LOCK
  if(OVreturn_IS_OK( (what_id= get_reinit_id(I,what)))) {
    ok = ExecutiveReinitialize(I->G,what_id.word,object_name);
  }
  PYMOL_API_UNLOCK
  return return_status_ok(ok);
}

PyMOLreturn_float PyMOL_CmdGetDistance(CPyMOL *I,
                                        char *selection1,
                                        char *selection2, 
                                        int state,
                                        int quiet)
{
  int ok=true;
  OrthoLineType s1="",s2="";
  PyMOLreturn_float result;
  PYMOL_API_LOCK

  if(ok) ok = (SelectorGetTmp(I->G,selection1,s1)>=0);
  if(ok) ok = (SelectorGetTmp(I->G,selection2,s2)>=0);

  if(ok) {
    ok = ExecutiveGetDistance(I->G,s1,s2,&result.value,state);
    result.status = get_status_ok(ok);    
  } else {
    result.status = PyMOLstatus_FAILURE;
    result.value = -1.0F;
  }
  SelectorFreeTmp(I->G,s1);
  SelectorFreeTmp(I->G,s2);
  PYMOL_API_UNLOCK
  return result;
}

PyMOLreturn_float PyMOL_CmdDistance(CPyMOL *I,
                                    char *name,
                                    char *selection1,
                                    char *selection2, 
                                    int mode,
                                    float cutoff,
                                    int label, 
                                    int reset,
                                    int zoom,
                                    int state,
                                    int quiet)
{
  int ok=true;
  OrthoLineType s1="",s2="";
  PyMOLreturn_float result;
  PYMOL_API_LOCK

  if(ok) ok = (SelectorGetTmp(I->G,selection1,s1)>=0);
  if(ok) ok = (SelectorGetTmp(I->G,selection2,s2)>=0);

  if(ok) {
    ok = ExecutiveDist(I->G,&result.value,name,s1,s2,
                       mode,cutoff,label,quiet,reset,state,zoom);
    result.status = get_status_ok(ok);
  } else {
    result.status = PyMOLstatus_FAILURE;
    result.value = -1.0F;
  }
  SelectorFreeTmp(I->G,s1);
  SelectorFreeTmp(I->G,s2);
  PYMOL_API_UNLOCK
  return result;
}

PyMOLreturn_float PyMOL_CmdGetAngle(CPyMOL *I,
                                     char *selection1,
                                     char *selection2,
                                     char *selection3,
                                     int state,
                                     int quiet)
{
  int ok=true;
  OrthoLineType s1="",s2="",s3="";
  PyMOLreturn_float result;
  PYMOL_API_LOCK

  if(ok) ok = (SelectorGetTmp(I->G,selection1,s1)>=0);
  if(ok) ok = (SelectorGetTmp(I->G,selection2,s2)>=0);
  if(ok) ok = (SelectorGetTmp(I->G,selection3,s3)>=0);

  if(ok) {
    ok = ExecutiveGetAngle(I->G,s1,s2,s3,&result.value,state);
    result.status = get_status_ok(ok);
  } else {
    result.status = PyMOLstatus_FAILURE;
    result.value = 0.0F;
  }
  SelectorFreeTmp(I->G,s1);
  SelectorFreeTmp(I->G,s2);
  SelectorFreeTmp(I->G,s3);
  PYMOL_API_UNLOCK
  return result;
}

PyMOLreturn_float PyMOL_CmdAngle(CPyMOL *I,
                                 char *name,
                                 char *selection1,
                                 char *selection2, 
                                 char *selection3, 
                                 int mode,
                                 int label, 
                                 int reset,
                                 int zoom,
                                 int state,
                                 int quiet)
{
  int ok=true;
  OrthoLineType s1="",s2="",s3="";
  PyMOLreturn_float result;
  PYMOL_API_LOCK

  if(ok) ok = (SelectorGetTmp(I->G,selection1,s1)>=0);
  if(ok) ok = (SelectorGetTmp(I->G,selection2,s2)>=0);
  if(ok) ok = (SelectorGetTmp(I->G,selection3,s3)>=0);

  if(ok) {
    ok = ExecutiveAngle(I->G,&result.value,name,s1,s2,s3,
                       mode,label,reset,zoom,quiet,state);
    result.status = get_status_ok(ok);
  } else {
    result.status = PyMOLstatus_FAILURE;
    result.value = -1.0F;
  }
  SelectorFreeTmp(I->G,s1);
  SelectorFreeTmp(I->G,s2);
  SelectorFreeTmp(I->G,s3);
  PYMOL_API_UNLOCK
  return result;
}


PyMOLreturn_float PyMOL_CmdGetDihedral(CPyMOL *I,
                                        char *selection1,
                                        char *selection2,
                                        char *selection3,
                                        char *selection4,
                                        int state,
                                        int quiet)
{
  int ok=true;
  OrthoLineType s1="",s2="",s3="",s4="";
  PyMOLreturn_float result;
  PYMOL_API_LOCK

  if(ok) ok = (SelectorGetTmp(I->G,selection1,s1)>=0);
  if(ok) ok = (SelectorGetTmp(I->G,selection2,s2)>=0);
  if(ok) ok = (SelectorGetTmp(I->G,selection3,s3)>=0);
  if(ok) ok = (SelectorGetTmp(I->G,selection4,s4)>=0);

   if(ok) {
    ok = ExecutiveGetDihe(I->G,s1,s2,s3,s4,&result.value,state);
    result.status = get_status_ok(ok);
 } else {
    result.status = PyMOLstatus_FAILURE;
    result.value = 0.0F;
  }
  SelectorFreeTmp(I->G,s1);
  SelectorFreeTmp(I->G,s2);
  SelectorFreeTmp(I->G,s3);
  SelectorFreeTmp(I->G,s4);
  PYMOL_API_UNLOCK
  return result;
}

PyMOLreturn_float PyMOL_CmdDihedral(CPyMOL *I,
                                    char *name,
                                    char *selection1,
                                    char *selection2, 
                                    char *selection3, 
                                    char *selection4, 
                                    int mode,
                                    int label, 
                                    int reset,
                                    int zoom,
                                    int state,
                                    int quiet)
{
  int ok=true;
  OrthoLineType s1="",s2="",s3="",s4="";
  PyMOLreturn_float result;
  PYMOL_API_LOCK
    
  if(ok) ok = (SelectorGetTmp(I->G,selection1,s1)>=0);
  if(ok) ok = (SelectorGetTmp(I->G,selection2,s2)>=0);
  if(ok) ok = (SelectorGetTmp(I->G,selection3,s3)>=0);
  if(ok) ok = (SelectorGetTmp(I->G,selection4,s4)>=0);
  
  if(ok) {
    ok = ExecutiveDihedral(I->G,&result.value,name,s1,s2,s3,s4,
                       mode,label,reset,zoom,quiet,state);
    result.status = get_status_ok(ok);
  } else {
    result.status = PyMOLstatus_FAILURE;
    result.value = -1.0F;
  }
  SelectorFreeTmp(I->G,s1);
  SelectorFreeTmp(I->G,s2);
  SelectorFreeTmp(I->G,s3);
  SelectorFreeTmp(I->G,s4);
  PYMOL_API_UNLOCK
  return result;
}


PyMOLreturn_status PyMOL_CmdIsodot(CPyMOL *I, char *name, char *map_name, float level, char *selection,
                                   float buffer, int state, float carve, int source_state, int quiet)
{
  int ok=true;
  OrthoLineType s1="";
  PyMOLreturn_status result = { PyMOLstatus_FAILURE };
  int box_mode = 0;
  PYMOL_API_LOCK
  
  if(selection && selection[0]) {  
    if(ok) ok = (SelectorGetTmp(I->G,selection,s1)>=0);
    if(ok) box_mode = 1;
  }

  if(ok) {
      ok = ExecutiveIsomeshEtc(I->G, name, map_name, level, s1, buffer,
                             state-1, carve, source_state-1, quiet, 1, box_mode, level);
    result.status = get_status_ok(ok);
  } else {
    result.status = PyMOLstatus_FAILURE;
  }
  SelectorFreeTmp(I->G,s1);
  PYMOL_API_UNLOCK
  return result;
  
}

PyMOLreturn_status PyMOL_CmdIsomesh(CPyMOL *I, char *name, char *map_name, float level, char *selection,
                                    float buffer, int state, float carve, int source_state, int quiet)
{
  int ok=true;
  OrthoLineType s1="";
  PyMOLreturn_status result = { PyMOLstatus_FAILURE };
  int box_mode = 0;
  PYMOL_API_LOCK
  
  if(selection && selection[0]) {  
    if(ok) ok = (SelectorGetTmp(I->G,selection,s1)>=0);
    if(ok) box_mode = 1;
  }
  if(ok) {
      ok = ExecutiveIsomeshEtc(I->G, name, map_name, level, s1, buffer,
                             state-1, carve, source_state-1, quiet, 0, box_mode, level);
    result.status = get_status_ok(ok);
  } else {
    result.status = PyMOLstatus_FAILURE;
  }
  SelectorFreeTmp(I->G,s1);
  PYMOL_API_UNLOCK
  return result;
}

PyMOLreturn_status PyMOL_CmdIsosurface(CPyMOL *I, char *name, char *map_name, float level, char *selection,
                                       float buffer, int state, float carve, int source_state, int side, 
                                       int mode, int quiet)
{
  int ok=true;
  OrthoLineType s1="";
  PyMOLreturn_status result = { PyMOLstatus_FAILURE };
  int box_mode = 0;
  PYMOL_API_LOCK
  
  if(selection && selection[0]) {  
    if(ok) ok = (SelectorGetTmp(I->G,selection,s1)>=0);
    if(ok) box_mode = 1;
  }

  if(ok) {
    ok = ExecutiveIsosurfaceEtc(I->G, name, map_name, level, s1, buffer,
                                state-1, carve, source_state-1, side, quiet, mode, box_mode);
    result.status = get_status_ok(ok);
  } else {
    result.status = PyMOLstatus_FAILURE;
  }
  SelectorFreeTmp(I->G,s1);
  PYMOL_API_UNLOCK
  return result;
}

PyMOLreturn_status PyMOL_CmdGradient(CPyMOL *I, char *name, char *map_name, float minimum, float maximum,
                                     char *selection, float buffer, int state, float carve, int source_state, 
                                     int quiet)
{
  int ok=true;
  OrthoLineType s1="";
  PyMOLreturn_status result = { PyMOLstatus_FAILURE };
  int box_mode = 0;
  PYMOL_API_LOCK
  
  if(selection && selection[0]) {  
    if(ok) ok = (SelectorGetTmp(I->G,selection,s1)>=0);
    if(ok) box_mode = 1;
  }

  if(ok) {
      ok = ExecutiveIsomeshEtc(I->G, name, map_name, minimum, s1, buffer,
                             state-1, carve, source_state-1, quiet, 3, box_mode, maximum);
    result.status = get_status_ok(ok);
  } else {
    result.status = PyMOLstatus_FAILURE;
  }
  SelectorFreeTmp(I->G,s1);
  PYMOL_API_UNLOCK
  return result;
}

PyMOLreturn_float PyMOL_CmdIsolevel(CPyMOL *I,char *name, float level, int state, int query, int quiet)
{
  int ok=true;
  OrthoLineType s1="";
  PyMOLreturn_float result;
  PYMOL_API_LOCK

  if(ok) {
    ok = ExecutiveIsolevel(I->G,name,level,state-1,query,&result.value,quiet);

    result.status = get_status_ok(ok);
  } else {
    result.status = PyMOLstatus_FAILURE;
    result.value = 0.0F;
  }
  SelectorFreeTmp(I->G,s1);
  PYMOL_API_UNLOCK
  return result;
}

static int word_count(char *src) /* only works for ascii */
{
  int cnt = 0;
  while((*src)&&((*src)<33)) /* skip leading whitespace */
    src++;
  while(*src) {
    if((*src)>32) {
      cnt++;
      while((*src)&&((*src)>32))
        src++;
    }
    while((*src)&&((*src)<33))
      src++;
  }
  return cnt;
}

static char *next_word(char *src,char *dst,int buf_size) /* only works for ascii */
{
  while((*src)&&((*src)<33)) /* skip leading whitespace */
    src++;
  while(*src) {
    if((*src)>32) {
      while((*src)&&((*src)>32)&&(buf_size>1)) {
        *(dst++) = *(src++);
        buf_size--;
      }
      break;
    }
  }
  dst[0]=0;
  return src;
}

PyMOLreturn_status PyMOL_CmdRampNew(CPyMOL *I,char *name, char *map, float *range, 
                                    int n_level, char *color, int state, char *selection,
                                    float beyond, float within, float sigma,
                                    int zero, int calc_mode, int quiet)

{
  int ok = true;
  PyMOLreturn_status result = { PyMOLstatus_FAILURE };
  OrthoLineType s1="";
  float *color_vla = NULL;
  float *range_vla = NULL;
  PYMOL_API_LOCK
  if(selection && selection[0]) {  
    if(ok) ok = (SelectorGetTmp(I->G,selection,s1)>=0);
  }
  if(ok) {
    if(range&&n_level) {
      range_vla = VLAlloc(float,n_level);
      UtilCopyMem(range_vla,range,sizeof(float)*n_level);
    }
  }

  if(ok&&color) {
    int n_color = word_count(color);
    /* to do */
    if(color&&n_color) {
      color_vla = VLAlloc(float,n_color*3);
      if(color_vla) {
        WordType colorName;
        int a;
        for(a=0;a<n_color;a++) {
          color = next_word(color,colorName,sizeof(colorName));
          {
            float *src = ColorGetNamed(I->G,colorName);
            float *dst = color_vla + 3*a;
            copy3f(src,dst);
          }
        }
      }
    }
  }
  if(ok) {
    ok = ExecutiveRampNew(I->G,name,map,range_vla,
                          color_vla,state,s1,beyond,within,sigma,
                          zero,calc_mode,quiet);
    result.status = get_status_ok(ok);
  } else {
    result.status = PyMOLstatus_FAILURE;
  }
  SelectorFreeTmp(I->G,s1);
  PYMOL_API_UNLOCK

  return result;
}

static PyMOLreturn_status Loader(CPyMOL *I,char *content,  char *content_type, 
                                 int content_length, char *content_format, 
                                 char *object_name, int state, 
                                 int discrete, int finish, 
                                 int quiet, int multiplex, int zoom)
{
  OVreturn_word result;
  int type_code = 0;
  int format_code = 0;
  int ok = true;
  WordType obj_name;
  
  if(!OVreturn_IS_OK( (result= OVLexicon_BorrowFromCString(I->Lex,content_type))))
    ok = false;
  else
    type_code = result.word;

  if(ok) {
    if(!OVreturn_IS_OK( (result= OVLexicon_BorrowFromCString(I->Lex,content_format))))
      ok = false;
    else
      format_code = result.word;
  }

  if(ok) {
    if((type_code != I->lex_filename) &&
       (type_code != I->lex_string) &&
       (type_code != I->lex_raw) &&
       (type_code != I->lex_cgo)) {
      
      ok = false;
    }
  }
  if(ok) {
    
    /* handling of multiplex option */
    
    if(multiplex==-2) /* use setting default value */
      multiplex = SettingGetGlobal_i(I->G,cSetting_multiplex);
    if(multiplex<0) /* default behavior is not to multiplex */
      multiplex = 0;
    
    /* handing of discete option */
    
    if(discrete<0) {/* use default discrete behavior for the file format 
                     * this will be the case for MOL2 and SDF */ 
      if(multiplex==1) /* if also multiplexing, then default discrete
                        * behavior is not load as discrete objects */
        discrete=0;
      else
        discrete=1; /* otherwise, allow discrete to be the default */
    }
    
    { /* if object_name is blank and content is a filename, then 
         compute the object_name from the file prefix */
      if((!object_name[0])&&(type_code == I->lex_filename)) {
        char *start, *stop;
        stop = start = content + strlen(content)-1;
        while(start>content) { /* known path separators */
          if((start[-1]==':')||
             (start[-1]=='\'')||
             (start[-1]=='/'))
            break;
          start--;
        }
        while(stop>start) {
          if(*stop=='.')
            break;
          stop--;
        }
        if(stop==start)
          stop = content + strlen(content);
        if((stop-start) >=sizeof(WordType))
          stop = start+sizeof(WordType)-1;
        { 
          char *p,*q;
          p=start;
          q=obj_name;
          while(p<stop) {
            *(q++)=*(p++);
          }
          *q=0;
          object_name = obj_name;
        }
      }
    }
    {
      int pymol_content_type = cLoadTypeUnknown;
      CObject *existing_object = NULL;

      /* convert text format strings into integral load types */

      if(format_code == I->lex_pdb) {
        if((type_code == I->lex_raw) ||(type_code == I->lex_string))
          pymol_content_type = cLoadTypePDBStr;
        else if( type_code == I->lex_filename)
          pymol_content_type = cLoadTypePDB;

      } else if(format_code == I->lex_mol2) {
        if((type_code == I->lex_raw) || (type_code == I->lex_string))
          pymol_content_type = cLoadTypeMOL2Str;
        else if( type_code == I->lex_filename)
          pymol_content_type = cLoadTypeMOL2;

      } else if(format_code == I->lex_mol) {
        if((type_code == I->lex_raw) || (type_code == I->lex_string))
          pymol_content_type = cLoadTypeMOLStr;
        else if( type_code == I->lex_filename)
          pymol_content_type = cLoadTypeMOL;

      } else if(format_code == I->lex_sdf) {
        if((type_code == I->lex_raw) || (type_code == I->lex_string))
          pymol_content_type = cLoadTypeSDF2Str;
        else if( type_code == I->lex_filename)
          pymol_content_type = cLoadTypeSDF2;

      } else if(format_code == I->lex_ccp4) {
        if((type_code == I->lex_raw) || (type_code == I->lex_string))
          pymol_content_type = cLoadTypeCCP4Str;

      } else if(format_code == I->lex_xplor) {
        if((type_code == I->lex_raw) || (type_code == I->lex_string))
          pymol_content_type = cLoadTypeXPLORStr;
        else if( type_code == I->lex_filename)
          pymol_content_type = cLoadTypeXPLORMap;

      } else if(format_code == I->lex_phi) {
        if((type_code == I->lex_raw) || (type_code == I->lex_string))
          pymol_content_type = cLoadTypePHIStr;
        else if(type_code == I->lex_filename)
          pymol_content_type = cLoadTypePHIMap;

      } else if(format_code == I->lex_cgo) {
        if(type_code == I->lex_cgo) {
          pymol_content_type = cLoadTypeCGO;
        }
      }

      if(pymol_content_type != cLoadTypeUnknown) {
        existing_object = ExecutiveGetExistingCompatible(I->G,
                                                         object_name,
                                                         pymol_content_type);
      }

      /* measure the length if it wasn't provided */

      if(content_length<0) {
        if(type_code == I->lex_string)
          content_length = strlen(content);
      }

      switch(pymol_content_type) {
      case cLoadTypePDB:
      case cLoadTypePDBStr:
      case cLoadTypeMOL:
      case cLoadTypeMOLStr:
      case cLoadTypeMOL2:
      case cLoadTypeMOL2Str:
      case cLoadTypeSDF2:
      case cLoadTypeSDF2Str:
      case cLoadTypeXPLORMap:
      case cLoadTypeXPLORStr:
      case cLoadTypePHIMap:
      case cLoadTypePHIStr:
      case cLoadTypeCCP4Map:
      case cLoadTypeCCP4Str:
      case cLoadTypeCGO:
        ok = ExecutiveLoad(I->G, existing_object, 
                           content, content_length, 
                           pymol_content_type,
                           object_name, 
                           state-1,  zoom,
                           discrete, finish,
                           multiplex, quiet, NULL);
        break;
      default:
        ok=false;
        break;
      }
    }
  }
  return return_status_ok(ok);
}

PyMOLreturn_status PyMOL_CmdLoad(CPyMOL *I,char *content,  
                                 char *content_type,
                                 char *content_format, 
                                 char *object_name, int state, 
                                 int discrete, int finish, 
                                 int quiet, int multiplex, int zoom)
{
  PyMOLreturn_status status = { PyMOLstatus_FAILURE };
  PYMOL_API_LOCK
  status = Loader(I,content, content_type, -1, content_format, object_name, 
                state, discrete, finish, quiet, multiplex, zoom);
  PYMOL_API_UNLOCK
  return status;
}
     
PyMOLreturn_status PyMOL_CmdLoadRaw(CPyMOL *I,char *content, 
                                    int content_length,
                                    char *content_format, 
                                    char *object_name, int state, 
                                    int discrete, int finish, 
                                    int quiet, int multiplex, int zoom)
{
  PyMOLreturn_status status = { PyMOLstatus_FAILURE };
  PYMOL_API_LOCK
  status = Loader(I,content, "raw", content_length, content_format,
                object_name, state, discrete, finish, quiet, multiplex, zoom);
  PYMOL_API_UNLOCK
  return status;
}

PyMOLreturn_status PyMOL_CmdLoadCGO(CPyMOL *I,float *content, 
                                    int content_length,
                                    char *object_name, int state, 
                                    int quiet, int zoom)
{
  PyMOLreturn_status status = { PyMOLstatus_FAILURE };
  PYMOL_API_LOCK
  status = Loader(I, (char*)content, "cgo", content_length, "cgo",
                  object_name, state, 0, 1, quiet, 0, zoom);
  PYMOL_API_UNLOCK
  return status;
}

PyMOLreturn_status PyMOL_CmdCreate(CPyMOL *I, char *name, 
                                   char *selection, int source_state,
                                   int target_state, int discrete, 
                                   int zoom, int quiet, int singletons,
                                   char *extract)
{
  int ok = true;
  PYMOL_API_LOCK
  ok = ExecutiveSeleToObject(I->G, name, selection, source_state, target_state, 
                                 discrete, zoom, quiet, singletons);
  PYMOL_API_UNLOCK
  return return_status_ok(ok);
}

static const CPyMOLOptions Defaults = {
 true, /* pmgui */
#ifndef _PYMOL_NOPY
  true, /* internal_gui*/
#else
  false, 
#endif
#ifndef _PYMOL_NOPY
  true, /* show_splash */
#else
  false,
#endif
#ifndef _PYMOL_NOPY
  1,   /* internal_feedback */
#else
  0, 
#endif
  true, /* security */
  false, /* game mode */
  0, /* force_stereo */
  640, /* winX */
  480, /* winY */
  false, /* blue_line */
  0, /* winPX */
  175, /* winPY */
  true, /* external_gui */
  true, /* siginthand */
  false, /* reuse helper */
  false, /* auto reinitialize */
  false, /* keep thread alive */
  false, /* quiet */
  false, /* incentive product */
  "", /* after_load_script */
  0, /* multisample */
  1, /* window_visible */
  0, /* read_stdin */
  0, /* presentation */
  0, /* defer builds mode */
  0, /* full screen mode */
  -1, /* sphere mode */
  0, /* stereo capable */
  0, /* stereo mode */
  -1, /* zoom mode */
  0, /* no quit */
};

CPyMOLOptions *PyMOLOptions_New(void)
{
  CPyMOLOptions *result = NULL;
  result = Calloc(CPyMOLOptions,1);
  if(result)
    *result = Defaults;
  return result;
}

#ifndef _PYMOL_NOPY
static void init_python(int argc, char *argv[])
{
    Py_Initialize();
    if(argv) {
	  PySys_SetArgv(argc,argv);
    }

/* BEGIN PROPRIETARY CODE SEGMENT (see disclaimer in "os_proprietary.h") */ 
#ifdef _MACPYMOL_XCODE
    /* there appears to be a bug or a race collection in the garbage
       collector of the Python version that ships with Mac OS --
       better to potentially leak a little RAM than crash unexpectedly
       in a _PyObject_GC_Del call.
       BTW: PyMOL doesn't itself need the GC, but end-user code
       might. */
    PyRun_SimpleString("import gc");
    PyRun_SimpleString("gc.disable()");
#endif
/* END PROPRIETARY CODE SEGMENT */

	PyEval_InitThreads();
	
#ifdef _PYMOL_OWN_INTERP
	{ /* NOTE this doesn't work 'cause we can't unpack code in restricted environments! */
	
	   /* get us a brand new interpreter, independent of any other */
	   
	   PyThreadState *tstate = Py_NewInterpreter();
	   
	   /* now release the first interpreter and use the second */
       PyThreadState_Swap(tstate);
	}
#endif

    PyUnicode_SetDefaultEncoding("utf-8"); /* is this safe & legal? */
	PyRun_SimpleString("import sys");
	PyRun_SimpleString("import os");
	PyRun_SimpleString("sys.path.insert(0,os.environ['PYMOL_PATH']+'/modules')");

/* BEGIN PROPRIETARY CODE SEGMENT (see disclaimer in "os_proprietary.h") */ 
#if 0
/* Not currently needed, since we've switched MacPyMOL away from the system Python */
#ifdef _MACPYMOL_XCODE
    { /* add an architecture-dependent search path for platform-specific binary modules
		here we're cheating by using endianness instead of figuring out how to get the
		true architecture from OS X */
		unsigned int val = 0x01020304;
		unsigned int *i_ptr=&val;
		char *c_ptr = (char*)i_ptr;
		if(*c_ptr==0x01) {
			PyRun_SimpleString("sys.path.insert(0,os.environ['PYMOL_PATH']+'/modules/ppc')");
		} else {
			PyRun_SimpleString("sys.path.insert(0,os.environ['PYMOL_PATH']+'/modules/i386')");
		}
	}
#endif
#endif

/* END PROPRIETARY CODE SEGMENT */

	PyRun_SimpleString("import __main__");
    {
		PyObject *P_main = PyImport_AddModule("__main__");
		if(!P_main) printf("PyMOL can't find '__main__'\n");

		/* set up a dry run through 'import pymol' */
 
		PyObject_SetAttrString(P_main,"pymol_launch",PyInt_FromLong(3));
    }
 
    /* initiate PyMOL dry run to create __main__.pymol */
	
 	PyRun_SimpleString("import pymol");
	
    /* parse arguments */

	PyRun_SimpleString("pymol.invocation.parse_args(sys.argv)");

}

/* WARNING: the routine below only works with the global singleton
   model -- not Python-enabled PyMOL instances */

CPyMOLOptions *PyMOLOptions_NewWithPython(int argc, char *argv[])
{
  CPyMOLOptions *result = PyMOLOptions_New();
  
  /* use Python to parse options based on the command line */
  
  init_python(argc,argv);
  PGetOptions(result);
  return result;
}
#endif

void PyMOLOptions_Free(CPyMOLOptions *options)
{
  FreeP(options);
}

void PyMOL_ResetProgress(CPyMOL *I)
{
  I->ProgressChanged = true;
  UtilZeroMem(I->Progress, sizeof(int)*6);
}

void PyMOL_SetProgress(CPyMOL *I,int offset, int current, int range)
{
  switch(offset) {
  case PYMOL_PROGRESS_SLOW:
  case PYMOL_PROGRESS_MED:
  case PYMOL_PROGRESS_FAST:
    if(current!=I->Progress[offset]) {
      I->Progress[offset] = current;
      I->ProgressChanged = true;
    }
    if(range!=I->Progress[offset+1]) {
      I->Progress[offset+1] = range;
      I->ProgressChanged = true;
    }
  }
}

int PyMOL_GetProgress(CPyMOL *I,int *progress,int reset)
{
  int a;
  int result = I->ProgressChanged;
  for(a=0;a<PYMOL_PROGRESS_SIZE;a++) {
    progress[a] = I->Progress[a];
  }
  if(reset) 
    I->ProgressChanged=false;
  return result;
}

int PyMOL_GetProgressChanged(CPyMOL *I,int reset)
{
  int result = I->ProgressChanged;
  if(reset)
    I->ProgressChanged=false;
  return result;
}

static CPyMOL *_PyMOL_New(void)
{
  CPyMOL *result = NULL;

  /* allocate global container */

  if( (result = Calloc(CPyMOL,1)) ) { /* all values initialized to zero */

    if( (result->G = Calloc(PyMOLGlobals,1)) ) {
      
      result->G->PyMOL = result; /* store the instance pointer */

      result->BusyFlag = false;
      result->InterruptFlag = false;
      PyMOL_ResetProgress(result);

      #ifndef _PYMOL_NOPY

      /* for the time being, the first PyMOL object created becomes
         the singleton object -- this is failsafe behavior designed to
         carry us through the transition to fully objectified PyMOL
         (PS note race in assignment covered by pymol2.pymol2_lock) */

      if(!SingletonPyMOLGlobals) {
        SingletonPyMOLGlobals=result->G;
      }
      #endif

      /* continue initialization */

    } else {
      FreeP(result);
    }
  }
  return result;
}
 
static void _PyMOL_Config(CPyMOL *I)
{
    I->G->HaveGUI = I->G->Option->pmgui;
    I->G->Security = I->G->Option->security;
}

CPyMOL *PyMOL_New(void)
{
  CPyMOL *result = _PyMOL_New();
  if(result && result->G) {
    result->G->Option = Calloc(CPyMOLOptions,1);
    if(result->G->Option)
      (*result->G->Option) = Defaults;
    _PyMOL_Config(result);
  }
  return result;
}

CPyMOL *PyMOL_NewWithOptions(CPyMOLOptions *option)
{
  CPyMOL *result = _PyMOL_New();
  if(result && result->G) {
    result->G->Option = Calloc(CPyMOLOptions,1);
    if(result->G->Option)
      *(result->G->Option) = *option;
    _PyMOL_Config(result);
  }
  result->G->StereoCapable = option->stereo_capable;
  return result;
}



void PyMOL_Start(CPyMOL *I)
{
  PyMOLGlobals *G=I->G;

  G->Context = OVContext_New();
  G->Lexicon = OVLexicon_New(G->Context->heap);

  if(OVreturn_IS_ERROR(PyMOL_InitAPI(I))) {
    printf("ERROR: PyMOL internal C API initialization failed.\n");
  }

  MemoryCacheInit(G);
  FeedbackInit(G,G->Option->quiet);
  WordInit(G);
  UtilInit(G);
  ColorInit(G);
  CGORendererInit(G);
  SettingInitGlobal(G,true,true,false);  
  SettingSetGlobal_i(G,cSetting_internal_gui,G->Option->internal_gui);
  SettingSetGlobal_i(G,cSetting_internal_feedback,G->Option->internal_feedback);
  TextureInit(G);
  TypeInit(G);
  TextInit(G);
  CharacterInit(G);
  PlugIOManagerInit(G);
  SphereInit(G);
  OrthoInit(G,G->Option->show_splash);
  WizardInit(G); /* must come after ortho */
  SceneInit(G);
  MovieInit(G);
  SelectorInit(G);
  SeqInit(G);
  SeekerInit(G);
  ButModeInit(G);
  ControlInit(G);
  AtomInfoInit(G);
  SculptCacheInit(G);
  VFontInit(G);
  ExecutiveInit(G);
  IsosurfInit(G);
  TetsurfInit(G);
  EditorInit(G);

#ifdef TRACKER_UNIT_TEST
  TrackerUnitTest(G);
#endif

/* BEGIN PROPRIETARY CODE SEGMENT (see disclaimer in "os_proprietary.h") */ 
#ifdef _MACPYMOL_XCODE
	SettingSetGlobal_b(G,cSetting_stereo_double_pump_mono, true);
    if(G->Option->stereo_capable) {
		SettingSetGlobal_i(G,cSetting_stereo_mode, 1); 
       }
    /*    	SettingSetGlobal_i(G,cSetting_show_progress, 0);  */
#endif
/* END PROPRIETARY CODE SEGMENT */

  I->DrawnFlag = false;
  I->RedisplayFlag = true;
  G->Ready = true; 
}

#ifndef _PYMOL_NOPY

/* BEGIN PROPRIETARY CODE SEGMENT (see disclaimer in "os_proprietary.h") */ 
#ifndef _PYMOL_API_HAS_PYTHON

#ifdef _MACPYMOL_XCODE
#define _PYMOL_API_HAS_PYTHON
#else

#ifdef _PYMOL_LIB_HAS_PYTHON
#define _PYMOL_API_HAS_PYTHON
#endif

#endif
#endif

#ifdef _PYMOL_API_HAS_PYTHON
void init_cmd(void);
/* 
 * void initExtensionClass(void);
 * void initsglite(void);
*/
void init_champ(void);
#endif
/* END PROPRIETARY CODE SEGMENT */

void PyMOL_StartWithPython(CPyMOL *I)
{
	PyMOL_Start(I);
	
	{
		PyObject *P_main = PyImport_AddModule("__main__");
		if(!P_main) printf("PyMOL can't find '__main__'\n");

		/* set up for embedded-style launch */
		PyObject_SetAttrString(P_main,"pymol_launch",PyInt_FromLong(5));
    }
	
	/* initialize our embedded C modules */

#ifdef _PYMOL_API_HAS_PYTHON
    init_cmd();
    /* 
     * initExtensionClass();
     * initsglite();
     */

	init_champ();
#endif
    
	/* launch pymol's Python subsystems */
	
	PyRun_SimpleString("import sys;reload(sys.modules['pymol'])");
	
	/* now locate all the C to Python function hooks and objects we need */
	
#ifdef _MACPYMOL_XCODE	
	PInit(I->G,true);
#else
	PInit(I->G,false);
#endif

	/* and begin the initialization sequence */
	
	I->PythonInitStage = 1;
	
/* BEGIN PROPRIETARY CODE SEGMENT (see disclaimer in "os_proprietary.h") */ 
#ifdef _MACPYMOL_XCODE
  MacPyMOLOption = I->G->Option;
  MacPyMOLReady = &I->G->Ready;
#endif  
/* END PROPRIETARY CODE SEGMENT */
}

#endif


void PyMOL_Stop(CPyMOL *I)
{
  PyMOLGlobals *G=I->G;
  G->Terminating=true;
  TetsurfFree(G);
  IsosurfFree(G);
  WizardFree(G);
  SceneCleanupStereo(G);
  EditorFree(G);
  ExecutiveFree(G);
  VFontFree(G);
  SculptCacheFree(G);
  AtomInfoFree(G);
  ButModeFree(G);
  ControlFree(G);
  SeekerFree(G);
  SeqFree(G);
  SelectorFree(G);
  MovieFree(G);
  SceneFree(G);
  OrthoFree(G);
  SettingFreeGlobal(G);
  CharacterFree(G);
  TextFree(G);
  TypeFree(G);
  TextureFree(G);
  SphereFree(G);
  PlugIOManagerFree(G);
  PFree();
  CGORendererFree(G);
  ColorFree(G);
  UtilFree(G);
  WordFree(G);
  FeedbackFree(G);
  MemoryCacheDone(G);
  
  PyMOL_PurgeAPI(I);
  /*    printf("%d \n", OVLexicon_GetNActive(G->Lexicon));*/
  OVLexicon_Del(G->Lexicon);
  OVContext_Del(G->Context);   
#ifndef _PYMOL_NOPY
  FreeP(G->P_inst);
#endif

}



void PyMOL_Free(CPyMOL *I)
{
#ifndef _PYMOL_ACTIVEX
  PYMOL_API_LOCK
#endif
    /* take PyMOL down gracefully */
  PyMOLOptions_Free(I->G->Option);
  FreeP(I->G);
#ifndef _PYMOL_NOPY
#ifdef _PYMOL_OWN_INTERP
  if(I->PythonInitStage) { /* shut down this interpreter gracefully, then free the GIL for others to use */
    PBlock(G);
    { /* should this be moved into a PDestroy?() to clear out the thread record too? */
      PyThreadState *tstate = PyEval_SaveThread();
      PyEval_AcquireThread(tstate);
      Py_EndInterpreter(tstate);
      PyEval_ReleaseLock();
	}
  }
#endif
#endif
  FreeP(I);
  return;
#ifndef _PYMOL_ACTIVEX
  PYMOL_API_UNLOCK;
#endif
}

struct _PyMOLGlobals *PyMOL_GetGlobals(CPyMOL *I)
{
  return I->G;
}

struct _PyMOLGlobals **PyMOL_GetGlobalsHandle(CPyMOL *I)
{
  return &(I->G);
}

void PyMOL_LockAPIAndUnblock(CPyMOL *I)
{
  PyMOLGlobals *G = I->G;
  PLockAPIAndUnblock(G);
}

void PyMOL_BlockAndUnlockAPI(CPyMOL *I)
{
  PyMOLGlobals *G = I->G;
  PBlockAndUnlockAPI(G);
}

void PyMOL_AdaptToHardware(CPyMOL *I)
{
  PYMOL_API_LOCK
  
  PyMOLGlobals *G = I->G;
  if(G->HaveGUI) {
    PyMOL_PushValidContext(I);
    {
      char *vendor = (char*)glGetString(GL_VENDOR);
      char *renderer =(char*)glGetString(GL_RENDERER);
      char *version = (char*)glGetString(GL_VERSION);    
      if(vendor && version) {
        /* work around broken lighting under Windows GDI Generic */
        if((strcmp(vendor,"Microsoft Corporation")==0) &&
           (strcmp(renderer,"GDI Generic")==0)) {
          ExecutiveSetSettingFromString(I->G, cSetting_light_count, "1","",0,1,0);
          ExecutiveSetSettingFromString(I->G, cSetting_spec_direct, "0.7","",0,1,0);
        }
      }
    }
    PyMOL_PopValidContext(I);
  }
  PYMOL_API_UNLOCK
}

static void setup_gl_state(void) 
{
  /* get us into a well defined GL state */
  
  /*glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();*/
  
  glDisable(GL_ALPHA_TEST);
  glDisable(GL_AUTO_NORMAL);
  glDisable(GL_BLEND);
  glDisable(GL_COLOR_LOGIC_OP);
  glDisable(GL_COLOR_MATERIAL);
  glDisable(GL_CULL_FACE);
  
  glDisable(GL_DEPTH_TEST);
  glDisable(GL_DITHER);
  glDisable(GL_FOG);
  glDisable(GL_LIGHTING);
  glDisable(GL_LIGHT0);
  glDisable(GL_LIGHT1);
  glDisable(GL_LINE_SMOOTH);
  glDisable(GL_NORMALIZE);
  glDisable(GL_POLYGON_SMOOTH);
}

void PyMOL_Draw(CPyMOL *I)
{
  PYMOL_API_LOCK_MODAL
  
  PyMOLGlobals *G = I->G;

  if(I->ModalDraw) {
    if(G->HaveGUI) {
      PyMOL_PushValidContext(I);
      setup_gl_state();
    }

    {
      PyMOLModalDrawFn *fn = I->ModalDraw;
      I->ModalDraw = NULL; /* always resets to NULL! */
      fn(G);
    }
    
    if(G->HaveGUI) {
      PyMOL_PopValidContext(I);
    }
  } else {

    if(I->DraggedFlag) {
      if(ControlIdling(I->G)) {
        ExecutiveSculptIterateAll(I->G);
      }
      I->DraggedFlag = false;
    }

    
    if(G->HaveGUI) {
      
      PyMOL_PushValidContext(I);
      
      setup_gl_state();

      /* BEGIN PROPRIETARY CODE SEGMENT (see disclaimer in "os_proprietary.h") */ 
#ifdef _MACPYMOL_XCODE
      { /* on a mac, this can change if we've switched contexts...*/
        GLboolean state;
        glGetBooleanv(GL_STEREO, &state);
        G->StereoCapable = (int) state;
      }
#endif
      /* END PROPRIETARY CODE SEGMENT */
      
      if(!I->DrawnFlag) {
        SceneSetCardInfo(G,(char*)glGetString(GL_VENDOR),
                         (char*)glGetString(GL_RENDERER),
                         (char*)glGetString(GL_VERSION));
        if(G->Option->show_splash) {
          
          printf(" OpenGL graphics engine:\n");
          printf("  GL_VENDOR: %s\n",(char*)glGetString(GL_VENDOR));
          printf("  GL_RENDERER: %s\n",(char*)glGetString(GL_RENDERER));
          printf("  GL_VERSION: %s\n",(char*)glGetString(GL_VERSION));
          if(Feedback(G,FB_OpenGL,FB_Blather)) {
            printf("  GL_EXTENSIONS: %s\n",(char*)glGetString(GL_EXTENSIONS));
          }

          if(G->StereoCapable) {
            printf("  Quad-buffer stereo 3D capability detected.\n");
          } else if(!G->Option->quiet) {
            if(G->LaunchStatus & cPyMOLGlobals_LaunchStatus_StereoFailed) {
              OrthoAddOutput(G,"Error: The requested stereo 3D visualization mode is not available.");
            }
          }
          
          if(G->LaunchStatus & cPyMOLGlobals_LaunchStatus_MultisampleFailed) {
            OrthoAddOutput(G,"Error: The requested multisampling mode is not available.");
          }
        } 
        I->DrawnFlag = true;
      }
    } else {
      I->DrawnFlag = true;
    }
    
    I->RedisplayFlag = false;
    
    OrthoBusyPrime(G);
    ExecutiveDrawNow(G);
    
    if(I->ImageRequestedFlag) {
      if(SceneHasImage(G)) {
        I->ImageReadyFlag = true;
        I->ImageRequestedFlag = false;
        {
          int w, h;
          SceneGetImageSize(I->G,&w,&h);
        }
        
      } else {
        I->ImageReadyFlag = false;
      }
    } else if(I->ImageReadyFlag) {
      if(!SceneHasImage(G)) 
        I->ImageReadyFlag = false;
    }

    if(G->HaveGUI) PyMOL_PopValidContext(I);
  }
  PYMOL_API_UNLOCK
}

void PyMOL_Key(CPyMOL *I,unsigned char k, int x, int y, int modifiers)
{
  PYMOL_API_LOCK
  
  PyMOLGlobals *G = I->G;

  if(!WizardDoKey(G,k,x,y,modifiers))
    OrthoKey(G,k,x,y,modifiers);
  PYMOL_API_UNLOCK
}


void PyMOL_Special(CPyMOL *I,int k, int x, int y, int modifiers)
{
  PYMOL_API_LOCK
  
  PyMOLGlobals *G = I->G;

  int grabbed = false;
  char buffer[255];
  
  if(!grabbed)
    grabbed = WizardDoKey(G,(unsigned char)k,x,y,modifiers);
  
  switch(k) {
  case P_GLUT_KEY_UP:
  case P_GLUT_KEY_DOWN:
    grabbed=1;
    OrthoSpecial(G,k,x,y,modifiers);
    break;
  case P_GLUT_KEY_LEFT:
  case P_GLUT_KEY_RIGHT:      
    if(OrthoArrowsGrabbed(G)) {
      grabbed=1;
      OrthoSpecial(G,k,x,y,modifiers);
    }
    break;
  }
  
#ifndef _PYMOL_NOPY
  if(!grabbed) {
    sprintf(buffer,"_special %d,%d,%d,%d",k,x,y,modifiers);
    PLog(G,buffer,cPLog_pml);
    PParse(G,buffer);
    PFlush(G);
  }
#endif

  PYMOL_API_UNLOCK
}

void PyMOL_Reshape(CPyMOL *I,int width, int height, int force)
{
  PYMOL_API_LOCK
  PyMOLGlobals *G = I->G;

  G->Option->winX = width;
  G->Option->winY = height;

  OrthoReshape(G,width,height,force);
  PYMOL_API_UNLOCK
}

int PyMOL_Idle(CPyMOL *I)
{
  int did_work = false;
  PYMOL_API_TRYLOCK

  PyMOLGlobals *G = I->G;

  I->DraggedFlag = false;
  if(I->IdleAndReady<IDLE_AND_READY) {
    if(I->DrawnFlag)
      I->IdleAndReady++;
  }
  if(I->FakeDragFlag==1) {
    I->FakeDragFlag = false;
    OrthoFakeDrag(G);
    did_work = true;
  }
  
  if(ControlIdling(G)) {
    ExecutiveSculptIterateAll(G);
    ControlSdofIterate(G);
    did_work = true;
  }
  
  SceneIdle(G); 
  
  if(SceneRovingCheckDirty(G)) {
    SceneRovingUpdate(G);
    did_work = true;
  }
  
#ifndef _PYMOL_NOPY
  if(PFlush(G)) {
    did_work = true;
  }

  if(I->PythonInitStage>0) {
    if(I->PythonInitStage<2) {
      I->PythonInitStage++;
    } else {
      I->PythonInitStage=-1;
      PBlock(G);
      /* BEGIN PROPRIETARY CODE SEGMENT (see disclaimer in "os_proprietary.h") */ 
#ifdef _MACPYMOL_XCODE
      /* restore working directory if asked to */
      PRunStringModule(G,"if os.environ.has_key('PYMOL_WD'): os.chdir(os.environ['PYMOL_WD'])");
      PXDecRef(PyObject_CallMethod(G->P_inst->obj,"launch_gui","O",G->P_inst->obj));
#endif
      /* END PROPRIETARY CODE SEGMENT */
      PXDecRef(PyObject_CallMethod(G->P_inst->obj,"adapt_to_hardware","O",G->P_inst->obj));
      PXDecRef(PyObject_CallMethod(G->P_inst->obj,"exec_deferred","O",G->P_inst->obj));
      PUnblock(G);
      PFlush(G);
    }
  }
#endif

  /* reset the interrupt flag if we're not doing anything */

  if(!(did_work||I->ModalDraw))
    if(PyMOL_GetInterrupt(I,false))
      PyMOL_SetInterrupt(I,false);

  PYMOL_API_UNLOCK_NO_FLUSH
  return (did_work || I->ModalDraw);
}

void PyMOL_ExpireIfIdle(CPyMOL *I)
{
  PYMOL_API_LOCK

  PyMOLGlobals *G = I->G;
  int final_init_done = true;
#ifndef _PYMOL_NOPY
  final_init_done = (I->PythonInitStage==-1);
#endif

  if(!G->HaveGUI) {
    if(final_init_done) {
      if(!OrthoCommandWaiting(G)) {
        if((!G->Option->keep_thread_alive)&&
           (!G->Option->read_stdin)) {
          I->ExpireCount++;
          if(I->ExpireCount==10) {
            PParse(G,"_quit");
          }
        }
      }
    }
  }
  PYMOL_API_UNLOCK;
}

void PyMOL_NeedFakeDrag(CPyMOL *I)
{
  I->FakeDragFlag = true;
}

void PyMOL_NeedRedisplay(CPyMOL *I)
{
  I->RedisplayFlag = true;
}

void PyMOL_NeedSwap(CPyMOL *I)
{
  I->SwapFlag = true;
}

void PyMOL_NeedReshape(CPyMOL *I,int mode, int x, int y, int width, int height)
{
  PyMOLGlobals *G = I->G;
  if(width<0) {
    int h;
    BlockGetSize(SceneGetBlock(G),&width,&h);
    if(SettingGetGlobal_b(G,cSetting_internal_gui))
      width+=SettingGetGlobal_i(G,cSetting_internal_gui_width);
  }
  
  /* if height is negative, force a reshape based on the current height */
  
  if(height<0) { 
    int w;
    int internal_feedback;
    BlockGetSize(SceneGetBlock(G),&w,&height);
    internal_feedback = (int)SettingGet(G,cSetting_internal_feedback);
    if(internal_feedback)
      height+=(internal_feedback-1)*cOrthoLineHeight+cOrthoBottomSceneMargin;
    if(SettingGetGlobal_b(G,cSetting_seq_view)&&!SettingGetGlobal_b(G,cSetting_seq_view_overlay))
      height+=SeqGetHeight(G);
  }

  if(G->HaveGUI) {
    I->ReshapeFlag = true;
    I->Reshape[0] = mode;
    I->Reshape[1] = x;
    I->Reshape[2] = y;
    I->Reshape[3] = width;
    I->Reshape[4] = height;
    PyMOL_NeedRedisplay(I);
  } else {
    /* if no gui, then force immediate reshape */
    PyMOLGlobals *G = I->G;
    
    G->Option->winX = width;
    G->Option->winY = height;
    
    OrthoReshape(G,width,height,true);
  }
}

int PyMOL_GetIdleAndReady(CPyMOL *I)
{
  return (I->IdleAndReady==IDLE_AND_READY);
}

int PyMOL_GetReshape(CPyMOL *I)
{
  return I->ReshapeFlag;
}

PyMOLreturn_int_array PyMOL_GetReshapeInfo(CPyMOL *I,int reset)
{
	PyMOLreturn_int_array result = { PyMOLstatus_SUCCESS, PYMOL_RESHAPE_SIZE, NULL };
	if(reset)
		I->ReshapeFlag = false;
	result.array = VLAlloc(int,PYMOL_RESHAPE_SIZE);
	if(!result.array) {
		result.status = PyMOLstatus_FAILURE;
	} else {
		int a;
		for(a=0;a<PYMOL_RESHAPE_SIZE;a++)
			result.array[a] = I->Reshape[a];
	}
	
  return result;
}

void PyMOL_SetPassive(CPyMOL *I,int onOff)
{
  I->PassiveFlag = onOff;
}

void PyMOL_SetClickReady(CPyMOL *I, char *name, int index, int button, int mod,int x,int y)
{

  if(name && name[0] && (index>=0)) {
    I->ClickReadyFlag = true;
    strcpy(I->ClickedObject,name);
    I->ClickedIndex = index;
    I->ClickedButton = button;
    I->ClickedModifiers = mod;
    I->ClickedX = x;
    I->ClickedY = y;
  } else {
    I->ClickedObject[0] = 0;
    I->ClickReadyFlag = true;
    I->ClickedX = x;
    I->ClickedY = y;
    I->ClickedIndex = index;
    I->ClickedButton = button;
    I->ClickedModifiers = mod;
  }
}

int PyMOL_GetClickReady(CPyMOL *I, int reset)
{
  int result = I->ClickReadyFlag;
  if(reset) {
    I->ClickReadyFlag = false;
  }
  return result;
}

char *PyMOL_GetClickString(CPyMOL *I,int reset)
{
  char *result = NULL;
  PYMOL_API_LOCK
  int ready = I->ClickReadyFlag;
  if(reset)
    I->ClickReadyFlag = false;
  if(ready) {
    result = Alloc(char, OrthoLineLength+1);
    if(result) {
      WordType butstr="left", modstr="";
      result[0]=0;
      switch(I->ClickedButton) {
      case P_GLUT_SINGLE_LEFT:
        strcpy(butstr,"single_left");
        break;
      case P_GLUT_SINGLE_MIDDLE:
        strcpy(butstr,"single_middle");
        break;
      case P_GLUT_SINGLE_RIGHT:
        strcpy(butstr,"single_right");
        break;
      case P_GLUT_DOUBLE_LEFT:
        strcpy(butstr,"double_left");
        break;
      case P_GLUT_DOUBLE_MIDDLE:
      strcpy(butstr,"double_middle");
      break;
      case P_GLUT_DOUBLE_RIGHT:
        strcpy(butstr,"double_right");
        break;
      }
      if(cOrthoCTRL & I->ClickedModifiers) {
        if(modstr[0]) strcat(modstr," ");
        strcat(modstr,"ctrl");
      }
      if(cOrthoALT & I->ClickedModifiers) {
        if(modstr[0]) strcat(modstr," ");
        strcat(modstr,"alt");
      }
      if(cOrthoSHIFT & I->ClickedModifiers) {
        if(modstr[0]) strcat(modstr," ");
        strcat(modstr,"shift");
      }
      if(!I->ClickedObject[0]) {
          sprintf(result,
                  "type=none\nclick=%s\nmod_keys=%s\nx=%d\ny=%d\n",
                  butstr,
                  modstr,
                  I->ClickedX,
                  I->ClickedY);
      } else { 
        ObjectMolecule *obj = ExecutiveFindObjectMoleculeByName(I->G,I->ClickedObject);
        if(obj && (I->ClickedIndex < obj->NAtom)) {
          AtomInfoType *ai = obj->AtomInfo + I->ClickedIndex;
          sprintf(result,
                  "type=object:molecule\nobject=%s\nindex=%d\nrank=%d\nid=%d\nsegi=%s\nchain=%s\nresn=%s\nresi=%s\nname=%s\nalt=%s\nclick=%s\nmod_keys=%s\nx=%d\ny=%d\n",
                  I->ClickedObject,
                  I->ClickedIndex+1,
                  ai->rank,
                  ai->id,
                  ai->segi,
                  ai->chain,
                  ai->resn,
                  ai->resi,
                  ai->name,
                  ai->alt,
                  butstr,
                  modstr,
                  I->ClickedX,
                  I->ClickedY);
        }
      }
    }
  }
  PYMOL_API_UNLOCK
  return(result);
}

int PyMOL_GetImageReady(CPyMOL *I, int reset)
{
  int result = I->ImageReadyFlag;
  if(reset) {
    I->ImageReadyFlag = false;
  }
  return result;
}

PyMOLreturn_int_array PyMOL_GetImageInfo(CPyMOL *I)
{
  PyMOLreturn_int_array result = { PyMOLstatus_SUCCESS, 2, NULL };
  PYMOL_API_LOCK
  result.array = VLAlloc(int,2);
  if(!result.array) {
    result.status = PyMOLstatus_FAILURE;
  } else {
    SceneGetImageSize(I->G,result.array,result.array+1);
  }
  PYMOL_API_UNLOCK
  return result;
}

int PyMOL_GetImageData(CPyMOL *I, 
                       int width, int height,
                       int row_bytes, void *buffer, 
		       int mode, int reset)
{
  int ok=true;
  PYMOL_API_LOCK
  if(reset)
    I->ImageReadyFlag=false;
  ok = SceneCopyExternal(I->G,width, height, row_bytes,(unsigned char *)buffer,mode);
  PYMOL_API_UNLOCK
  return get_status_ok(ok);
}

int PyMOL_FreeResultString(CPyMOL *I,char *st)
{
  PYMOL_API_LOCK
  FreeP(st);
  PYMOL_API_UNLOCK
  return get_status_ok((st!=NULL));
}

int PyMOL_GetRedisplay(CPyMOL *I, int reset)
{
  int result = false;
  
  PYMOL_API_TRYLOCK

  PyMOLGlobals *G = I->G;
  result = I->RedisplayFlag;
  
  if(result) {
    if(SettingGet_b(G,NULL,NULL,cSetting_defer_updates)) {
      result = false;
    } else {
      if(reset)
        I->RedisplayFlag = false;
    }
  }
  PYMOL_API_UNLOCK_NO_FLUSH

  return (result || I->ModalDraw); /* always true when ModalDraw is set */
}

int PyMOL_GetPassive(CPyMOL *I, int reset)
{/* lock intentionally omitted */
  int result = I->PassiveFlag;
  if(reset)
    I->PassiveFlag = false;
  return result;
}

int PyMOL_GetModalDraw(CPyMOL *I)
{
  if(I) 
    return (I->ModalDraw != NULL);
  return false;
}

void PyMOL_SetModalDraw(CPyMOL *I, PyMOLModalDrawFn *fn)
{
  I->ModalDraw = fn;
}

int PyMOL_GetSwap(CPyMOL *I, int reset)
{/* lock intentionally omitted */
  int result = I->SwapFlag;
  if(reset)
    I->SwapFlag = false;
  return result;
}

int PyMOL_GetBusy(CPyMOL *I, int reset)
{/* lock intentionally omitted */
  int result = I->BusyFlag;
  if(reset)
    PyMOL_SetBusy(I,false);
  return result;
}


void PyMOL_SetBusy(CPyMOL *I, int value)
{/* lock intentionally omitted */
  if(!I->BusyFlag) 
    /* if we weren't busy before, then reset the progress indicators */
    PyMOL_ResetProgress(I);

  I->BusyFlag = value;

#if 0
  if(!I->BusyFlag) /* reset the interrupt flag once we're done being busy */
    PyMOL_SetInterrupt(I,false);
#endif
}

int PyMOL_GetInterrupt(CPyMOL *I, int reset)
{ /* lock intentionally omitted */
  if(I) {
    int result = I->InterruptFlag;
    
    if(reset)
      PyMOL_SetInterrupt(I,false);
    return result;
  } else 
    return false;
}

void PyMOL_SetInterrupt(CPyMOL *I, int value)
{/* lock intentionally omitted */
  if(I) {
    I->InterruptFlag = value;
    if(I->G) I->G->Interrupt = value;
  }
}

void PyMOL_Drag(CPyMOL *I,int x, int y, int modifiers)
{
  PYMOL_API_LOCK
  OrthoDrag(I->G,x,y,modifiers);
  I->DraggedFlag = true;
  PYMOL_API_UNLOCK
}

void PyMOL_Button(CPyMOL *I,int button, int state,int x, int y, int modifiers)
{
  PYMOL_API_LOCK
  OrthoButton(I->G,button,state,x,y,modifiers);
  PYMOL_API_UNLOCK
}

void PyMOL_SetSwapBuffersFn(CPyMOL *I, PyMOLSwapBuffersFn *fn)
{
  I->SwapFn = fn;
}

void PyMOL_SwapBuffers(CPyMOL *I)
{
  if(I->SwapFn && I->G->ValidContext) {
    I->SwapFn();
    I->SwapFlag = false;
  } else {
    I->SwapFlag = true;
  }
}

void PyMOL_RunTest(CPyMOL *I, int group, int test)
{
  PYMOL_API_LOCK
  TestPyMOLRun(I->G, group, test);
  PYMOL_API_UNLOCK
}

void PyMOL_PushValidContext(CPyMOL *I)
{
  if(I && I->G)
    I->G->ValidContext++;
}

void PyMOL_PopValidContext(CPyMOL *I)
{
  if(I && I->G && (I->G->ValidContext>0))
    I->G->ValidContext--;
}

void PyMOL_SetDefaultMouse(CPyMOL *I)
{
  PYMOL_API_LOCK
  PyMOLGlobals *G = I->G;

  ButModeSet(G,cButModeLeftNone,cButModeRotXYZ);
  ButModeSet(G,cButModeMiddleNone,cButModeTransXY);
  ButModeSet(G,cButModeRightNone,cButModeTransZ);

  ButModeSet(G,cButModeLeftShft,cButModePotentialClick);
  ButModeSet(G,cButModeMiddleShft,cButModePotentialClick);
  ButModeSet(G,cButModeRightShft,cButModeClipNF);

  ButModeSet(G,cButModeLeftCtrl,cButModePotentialClick);
  ButModeSet(G,cButModeMiddleCtrl,cButModePotentialClick);
  ButModeSet(G,cButModeRightCtrl,cButModePotentialClick);

  ButModeSet(G,cButModeLeftCtSh,cButModePotentialClick);
  ButModeSet(G,cButModeMiddleCtSh,cButModePotentialClick);
  ButModeSet(G,cButModeRightCtSh,cButModePotentialClick);

  ButModeSet(G,cButModeWheelNone,cButModeScaleSlab);
  ButModeSet(G,cButModeWheelShft,cButModeMoveSlab);
  ButModeSet(G,cButModeWheelCtrl,cButModeMoveSlabAndZoom);
  ButModeSet(G,cButModeWheelCtSh,cButModeTransZ);

  ButModeSet(G,cButModeMiddleCtSh,cButModeOrigAt);

  ButModeSet(G,cButModeLeftSingle,cButModeSimpleClick);
  ButModeSet(G,cButModeMiddleSingle,cButModeCent);
  ButModeSet(G,cButModeRightSingle,cButModeSimpleClick);

  ButModeSet(G,cButModeLeftDouble,cButModeSimpleClick);
  ButModeSet(G,cButModeRightDouble,cButModeSimpleClick);
  
  {
    int a;
    for(a=cButModeLeftShftDouble;a<=cButModeRightCtrlAltShftSingle;a++) {
      ButModeSet(G,a,cButModeSimpleClick);
    }
    for(a=cButModeLeftAlt;a<=cButModeRightCtrlAltShft;a++) {
      ButModeSet(G,a,cButModePotentialClick);
    }

  }
  G->Feedback->Mask[FB_Scene] &= ~(FB_Results); /* suppress click messages */
  PYMOL_API_UNLOCK
}

