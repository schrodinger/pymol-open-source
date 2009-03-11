#A* -------------------------------------------------------------------
#B* This file contains source code for the PyMOL computer program
#C* copyright 1998-2000 by Warren Lyford Delano of DeLano Scientific. 
#D* -------------------------------------------------------------------
#E* It is unlawful to modify or remove this copyright notice.
#F* -------------------------------------------------------------------
#G* Please see the accompanying LICENSE file for further information. 
#H* -------------------------------------------------------------------
#I* Additional authors of this source file include:
#-* 
#-* 
#-*
#Z* -------------------------------------------------------------------

if __name__=='pymol.setting':
    
    import traceback
    import string
    import types
    import selector
    from shortcut import Shortcut
    import cmd
    from cmd import _cmd,lock,lock_attempt,unlock,QuietException, \
          _feedback,fb_module,fb_mask, \
          DEFAULT_ERROR, DEFAULT_SUCCESS, _raising, is_ok, is_error        
    from cmd import is_string
    import re
    
    boolean_type = 1
    int_type     = 2
    float_type   = 3
    float3_type  = 4

    self = cmd
    
    class SettingIndex: # number, description, range, notes
        bonding_vdw_cutoff                 = 0
        min_mesh_spacing                   = 1
        dot_density                        = 2
        dot_mode                           = 3
        solvent_radius                     = 4
        sel_counter                        = 5
        bg_rgb                             = 6
        ambient                            = 7
        direct                             = 8
        reflect                            = 9
        light                              = 10 
        power                              = 11 
        antialias                          = 12 
        cavity_cull                        = 13 
        gl_ambient                         = 14 
        single_image                       = 15 
        movie_delay                        = 16 
        ribbon_power                       = 17 
        ribbon_power_b                     = 18 
        ribbon_sampling                    = 19 
        ribbon_radius                      = 20 
        stick_radius                       = 21 
        hash_max                           = 22 
        orthoscopic                        = 23 
        spec_reflect                       = 24 
        spec_power                         = 25 
        sweep_angle                        = 26 
        sweep_speed                        = 27 
        dot_hydrogens                      = 28 
        dot_radius                         = 29 
        ray_trace_frames                   = 30 
        cache_frames                       = 31 
        trim_dots                          = 32 
        cull_spheres                       = 33 
        test1                              = 34 
        test2                              = 35 
        surface_best                       = 36 
        surface_normal                     = 37 
        surface_quality                    = 38 
        surface_proximity                  = 39 
        normal_workaround                  = 40 
        stereo_angle                       = 41 
        stereo_shift                       = 42 
        line_smooth                        = 43 
        line_width                         = 44 
        half_bonds                         = 45 
        stick_quality                      = 46 
        stick_overlap                      = 47 
        stick_nub                          = 48 
        all_states                         = 49 
        pickable                           = 50 
        auto_show_lines                    = 51 
        idle_delay                         = 52 
        no_idle                            = 53 
        fast_idle                          = 54 
        slow_idle                          = 55 
        rock_delay                         = 56 
        dist_counter                       = 57 
        dash_length                        = 58 
        dash_gap                           = 59 
        auto_zoom                          = 60 
        overlay                            = 61 
        text                               = 62 
        button_mode                        = 63 
        valence                            = 64 
        nonbonded_size                     = 65 
        label_color                        = 66 
        ray_trace_fog                      = 67 
        spheroid_scale                     = 68 
        ray_trace_fog_start                = 69 
        spheroid_smooth                    = 70 
        spheroid_fill                      = 71 
        auto_show_nonbonded                = 72 
        cache_display                      = 73 
        mesh_radius                        = 74 
        backface_cull                      = 75 
        gamma                              = 76 
        dot_width                          = 77 
        auto_show_selections               = 78 
        auto_hide_selections               = 79 
        selection_width                    = 80 
        selection_overlay                  = 81 
        static_singletons                  = 82 
        max_triangles                      = 83 
        depth_cue                          = 84 
        specular                           = 85 
        shininess                          = 86 
        sphere_quality                     = 87 
        fog                                = 88 
        isomesh_auto_state                 = 89 
        mesh_width                         = 90 
        cartoon_sampling                   = 91 
        cartoon_loop_radius                = 92 
        cartoon_loop_quality               = 93 
        cartoon_power                      = 94 
        cartoon_power_b                    = 95 
        cartoon_rect_length                = 96 
        cartoon_rect_width                 = 97 
        internal_gui_width                 = 98 
        internal_gui                       = 99 
        cartoon_oval_length                = 100
        cartoon_oval_width                 = 101
        cartoon_oval_quality               = 102
        cartoon_tube_radius                = 103
        cartoon_tube_quality               = 104
        cartoon_debug                      = 105
        ribbon_width                       = 106
        dash_width                         = 107
        dash_radius                        = 108
        cgo_ray_width_scale                = 109
        line_radius                        = 110
        cartoon_round_helices              = 111
        cartoon_refine_normals             = 112
        cartoon_flat_sheets                = 113
        cartoon_smooth_loops               = 114
        cartoon_dumbbell_length            = 115
        cartoon_dumbbell_width             = 116
        cartoon_dumbbell_radius            = 117
        cartoon_fancy_helices              = 118
        cartoon_fancy_sheets               = 119
        ignore_pdb_segi                    = 120
        ribbon_throw                       = 121
        cartoon_throw                      = 122
        cartoon_refine                     = 123
        cartoon_refine_tips                = 124
        cartoon_discrete_colors            = 125
        normalize_ccp4_maps                = 126
        surface_poor                       = 127
        internal_feedback                  = 128
        cgo_line_width                     = 129
        cgo_line_radius                    = 130
        logging                            = 131
        robust_logs                        = 132
        log_box_selections                 = 133
        log_conformations                  = 134
        valence_size                       = 135
        valence_default                    = 135 # deprecated
        surface_miserable                  = 136
        ray_opaque_background              = 137
        transparency                       = 138
        ray_texture                        = 139
        ray_texture_settings               = 140
        suspend_updates                    = 141
        full_screen                        = 142
        surface_mode                       = 143
        surface_color                      = 144
        mesh_mode                          = 145
        mesh_color                         = 146
        auto_indicate_flags                = 147
        surface_debug                      = 148
        ray_improve_shadows                = 149
        smooth_color_triangle              = 150
        ray_default_renderer               = 151
        field_of_view                      = 152
        reflect_power                      = 153
        preserve_chempy_ids                = 154
        sphere_scale                       = 155
        two_sided_lighting                 = 156
        secondary_structure                = 157
        auto_remove_hydrogens              = 158
        raise_exceptions                   = 159
        stop_on_exceptions                 = 160
        sculpting                          = 161
        auto_sculpt                        = 162
        sculpt_vdw_scale                   = 163
        sculpt_vdw_scale14                 = 164
        sculpt_vdw_weight                  = 165
        sculpt_vdw_weight14                = 166
        sculpt_bond_weight                 = 167
        sculpt_angl_weight                 = 168
        sculpt_pyra_weight                 = 169
        sculpt_plan_weight                 = 170
        sculpting_cycles                   = 171
        sphere_transparency                = 172
        sphere_color                       = 173
        sculpt_field_mask                  = 174
        sculpt_hb_overlap                  = 175
        sculpt_hb_overlap_base             = 176
        legacy_vdw_radii                   = 177
        sculpt_memory                      = 178
        connect_mode                       = 179
        cartoon_cylindrical_helices        = 180
        cartoon_helix_radius               = 181
        connect_cutoff                     = 182
    #   save_pdb_ss                        = 183
        sculpt_line_weight                 = 184
        fit_iterations                     = 185
        fit_tolerance                      = 186
        batch_prefix                       = 187
        stereo_mode                        = 188
        cgo_sphere_quality                 = 189
        pdb_literal_names                  = 190
        wrap_output                        = 191
        fog_start                          = 192
        state                              = 193
        frame                              = 194
        ray_shadows                        = 195 # legacy
        ray_shadow                         = 195 # improved ease of use..(for ray_sh <Tab> )
        ribbon_trace_atoms                 = 196
        security                           = 197
        stick_transparency                 = 198
        ray_transparency_shadows           = 199
        session_version_check              = 200
        ray_transparency_specular          = 201
        stereo_double_pump_mono            = 202
        sphere_solvent                     = 203
        mesh_quality                       = 204
        mesh_solvent                       = 205
        dot_solvent                        = 206
        ray_shadow_fudge                   = 207
        ray_triangle_fudge                 = 208
        debug_pick                         = 209
        dot_color                          = 210
        mouse_limit                        = 211
        mouse_scale                        = 212
        transparency_mode                  = 213
        clamp_colors                       = 214
        pymol_space_max_red                = 215
        pymol_space_max_green              = 216
        pymol_space_max_blue               = 217
        pymol_space_min_factor             = 218
        roving_origin                      = 219
        roving_lines                       = 220
        roving_sticks                      = 221
        roving_spheres                     = 222
        roving_labels                      = 223
        roving_delay                       = 224
        roving_selection                   = 225
        roving_byres                       = 226
        roving_ribbon                      = 227
        roving_cartoon                     = 228
        roving_polar_contacts              = 229
        roving_polar_cutoff                = 230
        roving_nonbonded                   = 231
        float_labels                       = 232
        roving_detail                      = 233
        roving_nb_spheres                  = 234
        ribbon_color                       = 235
        cartoon_color                      = 236
        ribbon_smooth                      = 237
        auto_color                         = 238
        ray_interior_color                 = 240
        cartoon_highlight_color            = 241
        coulomb_units_factor               = 242
        coulomb_dielectric                 = 243
        ray_interior_shadows               = 244
        ray_interior_texture               = 245
        roving_map1_name                   = 246
        roving_map2_name                   = 247
        roving_map3_name                   = 248
        roving_map1_level                  = 249
        roving_map2_level                  = 250
        roving_map3_level                  = 251
        roving_isomesh                     = 252
        roving_isosurface                  = 253
        scenes_changed                     = 254
        gaussian_b_adjust                  = 255
        pdb_standard_order                 = 256
        cartoon_smooth_first               = 257
        cartoon_smooth_last                = 258
        cartoon_smooth_cycles              = 259
        cartoon_flat_cycles                = 260
        max_threads                        = 261
        show_progress                      = 262
        use_display_lists                  = 263
        cache_memory                       = 264
        simplify_display_lists             = 265
        retain_order                       = 266
        pdb_hetatm_sort                    = 267
        pdb_use_ter_records                = 268
        cartoon_trace_atoms                = 269
        ray_oversample_cutoff              = 270
        gaussian_resolution                = 271
        gaussian_b_floor                   = 272
        sculpt_nb_interval                 = 273
        sculpt_tors_weight                 = 274
        sculpt_tors_tolerance              = 275
        stick_ball                         = 276
        stick_ball_ratio                   = 277
        stick_fixed_radius                 = 278
        cartoon_transparency               = 279
        dash_round_ends                    = 280
        h_bond_max_angle                   = 281
        h_bond_cutoff_center               = 282
        h_bond_cutoff_edge                 = 283
        h_bond_power_a                     = 284
        h_bond_power_b                     = 285
        h_bond_cone                        = 286
        ss_helix_psi_target                = 287 
        ss_helix_psi_include               = 288
        ss_helix_psi_exclude               = 289
        ss_helix_phi_target                = 290
        ss_helix_phi_include               = 291
        ss_helix_phi_exclude               = 292
        ss_strand_psi_target               = 293
        ss_strand_psi_include              = 294
        ss_strand_psi_exclude              = 295
        ss_strand_phi_target               = 296
        ss_strand_phi_include              = 297
        ss_strand_phi_exclude              = 298
        movie_loop                         = 299
        pdb_retain_ids                     = 300
        pdb_no_end_record                  = 301
        cgo_dot_width                      = 302
        cgo_dot_radius                     = 303
        defer_updates                      = 304
        normalize_o_maps                   = 305
        swap_dsn6_bytes                    = 306
        pdb_insertions_go_first            = 307
        roving_origin_z                    = 308
        roving_origin_z_cushion            = 309
        specular_intensity                 = 310
        overlay_lines                      = 311
        ray_transparency_spec_cut          = 312
        internal_prompt                    = 313
        normalize_grd_maps                 = 314
        ray_blend_colors                   = 315
        ray_blend_red                      = 316
        ray_blend_green                    = 317
        ray_blend_blue                     = 318
        png_screen_gamma                   = 319
        png_file_gamma                     = 320
        editor_label_fragments             = 321
        internal_gui_control_size          = 322
        auto_dss                           = 323
        transparency_picking_mode          = 324
        virtual_trackball                  = 325
        pdb_reformat_names_mode            = 326
        ray_pixel_scale                    = 327
        label_font_id                      = 328
        pdb_conect_all                     = 329
        button_mode_name                   = 330
        surface_type                       = 331
        dot_normals                        = 332
        session_migration                  = 333
        mesh_normals                       = 334
        mesh_type                          = 335
        dot_lighting                       = 336
        mesh_lighting                      = 337
        surface_solvent                    = 338
        triangle_max_passes                = 339
        ray_interior_reflect               = 340
        internal_gui_mode                  = 341
        surface_carve_selection            = 342
        surface_carve_state                = 343
        surface_carve_cutoff               = 344
        surface_clear_selection            = 345
        surface_clear_state                = 346
        surface_clear_cutoff               = 347
        surface_trim_cutoff                = 348
        surface_trim_factor                = 349
        ray_max_passes                     = 350
        active_selections                  = 351
        ray_transparency_contrast          = 352
        seq_view                           = 353
        mouse_selection_mode               = 354
        seq_view_label_spacing             = 355
        seq_view_label_start               = 356
        seq_view_format                    = 357
        seq_view_location                  = 358
        seq_view_overlay                   = 359
        auto_classify_atoms                = 360
        cartoon_nucleic_acid_mode          = 361
        seq_view_color                     = 362
        seq_view_label_mode                = 363
        surface_ramp_above_mode            = 364
        stereo                             = 365
        wizard_prompt_mode                 = 366
        coulomb_cutoff                     = 367
        slice_track_camera                 = 368
        slice_height_scale                 = 369
        slice_height_map                   = 370
        slice_grid                         = 371
        slice_dynamic_grid                 = 372
        slice_dynamic_grid_resolution      = 373
        pdb_insure_orthogonal              = 374
        ray_direct_shade                   = 375
        stick_color                        = 376 
        cartoon_putty_radius               = 377 
        cartoon_putty_quality              = 378 
        cartoon_putty_scale_min            = 379 
        cartoon_putty_scale_max            = 380 
        cartoon_putty_scale_power          = 381 
        cartoon_putty_range                = 382 
        cartoon_side_chain_helper          = 383             
        surface_optimize_subsets           = 384 
        multiplex                          = 385 
        texture_fonts                      = 386 
        pqr_no_chain_id                    = 387 
        animation                          = 388 
        animation_duration                 = 389 
        scene_animation                    = 390 
        line_stick_helper                  = 391 
        ray_orthoscopic                    = 392 
        ribbon_side_chain_helper           = 393 
        selection_width_max                = 394 
        selection_width_scale              = 395 
        scene_current_name                 = 396 
        presentation                       = 397 
        presentation_mode                  = 398 
        pdb_truncate_residue_name          = 399 
        scene_loop                         = 400 
        sweep_mode                         = 401 
        sweep_phase                        = 402 
        scene_restart_movie_delay          = 403 
        mouse_restart_movie_delay          = 404 
        angle_size                         = 405 
        angle_label_position               = 406 
        dihedral_size                      = 407 
        dihedral_label_position            = 408 
        defer_builds_mode                  = 409 
        seq_view_discrete_by_state         = 410 
        scene_animation_duration           = 411 
        wildcard                           = 412 
        atom_name_wildcard                 = 413 
        ignore_case                        = 414 
        presentation_auto_quit             = 415 
        editor_auto_dihedral               = 416 
        presentation_auto_start            = 417 
        validate_object_names              = 418 
#        ray_pixel_scale_limit             = 419 
        auto_show_spheres                  = 420 
        sphere_mode                        = 421 
        sphere_point_max_size              = 422 
        sphere_point_size                  = 423 
        pdb_honor_model_number             = 424 
        rank_assisted_sorts                = 425 
        ribbon_nucleic_acid_mode           = 426 
        cartoon_ring_mode                  = 427 
        cartoon_ring_width                 = 428 
        cartoon_ring_color                 = 429 
        cartoon_ring_finder                = 430 
        cartoon_tube_cap                   = 431 
        cartoon_loop_cap                   = 432 
        nvidia_bugs                        = 433 
        image_dots_per_inch                = 434 
        opaque_background                  = 435 
        draw_frames                        = 436 
        show_alpha_checker                 = 437 
        matrix_mode                        = 438 
        editor_auto_origin                 = 439 
        session_file                       = 440 
        cgo_transparency                   = 441 
        legacy_mouse_zoom                  = 442 
        auto_number_selections             = 443 
        sculpt_vdw_vis_mode                = 444 
        sculpt_vdw_vis_min                 = 445 
        sculpt_vdw_vis_mid                 = 446 
        sculpt_vdw_vis_max                 = 447 
        cartoon_ladder_mode                = 448 
        cartoon_ladder_radius              = 449 
        cartoon_ladder_color               = 450 
        cartoon_nucleic_acid_color         = 451 
        cartoon_ring_transparency          = 452 
        label_size                         = 453 
        spec_direct                        = 454 
        light_count                        = 455 
        light2                             = 456 
        light3                             = 457 
        hide_underscore_names              = 458 
        selection_round_points             = 459 
        distance_exclusion                 = 460 
        h_bond_exclusion                   = 461 
        label_shadow_mode                  = 462 
        light4                             = 463 
        light5                             = 464 
        light6                             = 465 
        light7                             = 466 
        label_outline_color                = 467 
        ray_trace_mode                     = 468 
        ray_trace_gain                     = 469 
        selection_visible_only             = 470 
        label_position                     = 471 
        ray_trace_depth_factor             = 472 
        ray_trace_slope_factor             = 473 
        ray_trace_disco_factor             = 474 
        ray_shadow_decay_factor            = 475 
        ray_interior_mode                  = 476 
        ray_legacy_lighting                = 477 
        sculpt_auto_center                 = 478 
        pdb_discrete_chains                = 479 
        pdb_unbond_cations                 = 480 
        sculpt_tri_scale                   = 481 
        sculpt_tri_weight                  = 482 
        sculpt_tri_min                     = 483 
        sculpt_tri_max                     = 484 
        sculpt_tri_mode                    = 485 
        pdb_echo_tags                      = 486 
        connect_bonded                     = 487 
        spec_direct_power                  = 488 
        light8                             = 489 
        light9                             = 490 
        ray_shadow_decay_range             = 491 
        spec_count                         = 492 
        sculpt_min_scale                   = 493 
        sculpt_min_weight                  = 494 
        sculpt_min_min                     = 495 
        sculpt_min_max                     = 496 
        sculpt_max_scale                   = 497 
        sculpt_max_weight                  = 498 
        sculpt_max_min                     = 499 
        sculpt_max_max                     = 500 
        surface_circumscribe               = 501 
        sculpt_avd_weight                  = 502 
        sculpt_avd_gap                     = 503 
        sculpt_avd_range                   = 504 
        sculpt_avd_excl                    = 505 
        async_builds                       = 506 
        fetch_path                         = 507 
        cartoon_ring_radius                = 508 
        ray_color_ramps                    = 509 
        ray_hint_camera                    = 510 
        ray_hint_shadow                    = 511 
        stick_valence_scale                = 512 
        seq_view_alignment                 = 513 
        seq_view_unaligned_mode            = 514 
        seq_view_unaligned_color           = 515 
        seq_view_fill_char                 = 516 
        seq_view_fill_color                = 517 
        seq_view_label_color               = 518 
        surface_carve_normal_cutoff        = 519
        trace_atoms_mode                   = 520
        session_changed                    = 521
        ray_clip_shadows                   = 522
        mouse_wheel_scale                  = 523
        nonbonded_transparency             = 524
        ray_spec_local                     = 525
        line_color                         = 526
        ray_label_specular                 = 527
        mesh_skip                          = 528
        label_digits                       = 529
        label_distance_digits              = 530
        label_angle_digits                 = 531
        label_dihedral_digits              = 532
        surface_negative_visible           = 533
        surface_negative_color             = 534
        mesh_negative_visible              = 535
        mesh_negative_color                = 536
        group_auto_mode                    = 537
        group_full_member_names            = 538
        gradient_max_length                = 539
        gradient_min_length                = 540
        gradient_min_slope                 = 541
        gradient_normal_min_dot            = 542
        gradient_step_size                 = 543
        gradient_spacing                   = 544
        gradient_symmetry                  = 545
        ray_trace_color                    = 546
        group_arrow_prefix                 = 547
        suppress_hidden                    = 548
        session_compression                = 549
        movie_fps                          = 550
        ray_transparency_oblique           = 551
        ray_trace_trans_cutoff             = 552
        ray_trace_persist_cutoff           = 553
        ray_transparency_oblique_power     = 554
        ray_scatter                        = 555
        h_bond_from_proton                 = 556
        auto_copy_images                   = 557
        moe_separate_chains                = 558
        transparency_global_sort           = 559
        hide_long_bonds                    = 560
        auto_rename_duplicate_objects      = 561
        pdb_hetatm_guess_valences          = 562
        ellipsoid_quality                  = 563
        cgo_ellipsoid_quality              = 564
        movie_animate_by_frame             = 565
        ramp_blend_nearby_colors           = 566
        auto_defer_builds                  = 567
        ellipsoid_probability              = 568
        ellipsoid_scale                    = 569
        ellipsoid_color                    = 570
        ellipsoid_transparency             = 571
        movie_rock                         = 572
        cache_mode                         = 573
        dash_color                         = 574
        angle_color                        = 575
        dihedral_color                     = 576        
        grid_mode                          = 577
        cache_max                          = 578
        grid_slot                          = 579
        grid_max                           = 580
        cartoon_putty_transform            = 581
        rock                               = 582
        cone_quality                       = 583
        pdb_formal_charges                 = 584
        ati_bugs                           = 585
        geometry_export_mode               = 586
        mouse_grid                         = 587
        mesh_cutoff                        = 588
        mesh_carve_selection               = 589
        mesh_carve_state                   = 590
        mesh_carve_cutoff                  = 591
        mesh_clear_selection               = 592
        mesh_clear_state                   = 593
        mesh_clear_cutoff                  = 594
        mesh_grid_max                      = 595
        session_cache_optimize             = 596
        sdof_drag_scale                    = 597
        scene_buttons_mode                 = 598
        scene_buttons                      = 599
        map_auto_expand_sym                = 600
        image_copy_always                  = 601
        max_ups                            = 602
        auto_overlay                       = 603
        stick_ball_color                   = 604
        stick_h_scale                      = 605
        sculpt_pyra_inv_weight             = 606
        keep_alive                         = 607
        fit_kabsch                         = 608
        stereo_dynamic_strength            = 609
        dynamic_width                      = 610
        dynamic_width_factor               = 611
        dynamic_width_min                  = 612
        dynamic_width_max                  = 613
        draw_mode                          = 614
        clean_electro_mode                 = 615
        valence_mode                       = 616
        show_frame_rate                    = 617
        
    setting_sc = Shortcut(SettingIndex.__dict__.keys())
    
    index_list = []
    name_dict = {}
    name_list = SettingIndex.__dict__.keys()
    name_list = filter(lambda x:x[0]!='_',name_list)
    name_list.sort()
    tmp_list = map(lambda x:(getattr(SettingIndex,x),x),name_list)
    for a in tmp_list:
        name_dict[a[0]]=a[1]
        index_list.append(a[0])

    boolean_dict = {
        "true" : 1,
        "false": 0,
        "on"   : 1,
        "off"  : 0,
        "1"    : 1,
        "0"    : 0,
        "1.0"  : 1,
        "0.0"  : 0,
        }

    boolean_sc = Shortcut(boolean_dict.keys())

    quote_strip_re = re.compile("^\'.*\'$|^\".*\"$")

    def _get_index(name):
        # this may be called from C, so don't raise any exceptions...
        result = setting_sc.interpret(name)
        if is_string(result):
            if hasattr(SettingIndex,result):
                return getattr(SettingIndex,result)
            else:
                return -1
        else:
            return -1

    def _get_name(index): # likewise, this can be called from C so no exceptions
        return name_dict.get(index,"")

    def get_index_list():
        return index_list

    def get_name_list():
        return name_list

    def dump_tuples():
        for idx in range(581):
            print _get_name(idx), cmd.get_setting_tuple(idx)
            
    ###### API functions

    def set_bond(name, value, selection1, selection2=None,
                 state=0, updates=1, log=0, quiet=1, _self=cmd):
        ''' 
DESCRIPTION

    "set_bond" changes per-bond settings for all bonds which exist
    between two selections of atoms.

USAGE

    set_bond name, value, selection1 [, selection2 ]

ARGUMENTS

    name = string: name of the setting

    value = string: new value to use

    selection1 = string: first set of atoms

    selection2 = string: seconds set of atoms {default: (selection1)}

EXAMPLE

    set_bond stick_transparency, 0.7, */n+c+ca+o


NOTES

    The following per-bond settings are currently implemented.  Others
    may seem to be recognized but will currently have no effect when
    set at the per-bond level.
    
    * valence
    * line_width
    * line_color
    * stick_radius
    * stick_color
    * stick_transparency

    Note that if you attempt to use the "set" command with a per-bond
    setting over a selection of atoms, the setting change will appear
    to take, but no change will be observed.
    
PYMOL API

    cmd.set_bond ( string name, string value,
                   string selection1,
                   string selection2,
                   int state, int updates, log=0, quiet=1)

       '''
        r = DEFAULT_ERROR
        selection1 = str(selection1)
        if selection2 == None:
            selection2 = selection1
        if log:
            if ',' in str(value):
                value = str(cmd.safe_list_eval(str(value)))
            if len(selection2):
                cmd.log("set_bond %s,%s,%s\n"%(str(name),str(value),str(selection1),str(selection2)))
            else:
                cmd.log("set_bond %s,%s\n"%(str(name),str(value),str(selection1)))           
        index = _get_index(str(name))
        if(index<0):
            print "Error: unknown setting '%s'."%name
            raise QuietException
        else:
            try:
                _self.lock(_self)
                type = _cmd.get_setting_tuple(_self._COb,int(index),str(""),int(-1))[0]
                if type==None:
                    print "Error: unable to get setting type."
                    raise QuietException
                try:
                    if type==1: # boolean (also support non-zero float for truth)
                        handled = 0
                        if boolean_sc.interpret(str(value))==None:
                            try: # number, non-zero, then interpret as TRUE
                                if not (float(value)==0.0):
                                    handled = 1
                                    v = (1,)
                                else:
                                    handled = 1
                                    v = (0,)
                            except:
                                pass
                        if not handled:
                            v = (boolean_dict[
                                boolean_sc.auto_err(
                                str(value),"boolean")],)
                    elif type==2: # int (also supports boolean language for 0,1)
                        if boolean_sc.has_key(str(value)):
                            v = (boolean_dict[
                                boolean_sc.auto_err(
                                str(value),"boolean")],)
                        else:
                            v = (int(value),)
                    elif type==3: # float
                        if boolean_sc.has_key(str(value)):
                            v = (float(boolean_dict[
                                boolean_sc.auto_err(
                                str(value),"boolean")]),)
                        else:
                            v = (float(value),)
                    elif type==4: # float3 - some legacy handling req.
                        if is_string(value):
                            if not ',' in value:
                                v = string.split(value)
                            else:
                                v = eval(value)
                        else:
                            v = value
                        v = (float(v[0]),float(v[1]),float(v[2]))
                    elif type==5: # color
                        v = (str(value),)
                    elif type==6: # string
                        vl = str(value)
                        # strip outermost quotes (cheesy approach)
                        if quote_strip_re.search(vl)!=None:
                            vl=vl[1:-1]
                        v = (vl,)
                    v = (type,v)
                    if len(selection1):
                        selection1=selector.process(selection1)
                    if len(selection2):
                        selection2=selector.process(selection2)                        
                    r = _cmd.set_bond(_self._COb,int(index),v,
                                 "("+selection1+")","("+selection2+")",
                                 int(state)-1,int(quiet),
                                 int(updates))
                except:
                    traceback.print_exc()
                    if(_feedback(fb_module.cmd,fb_mask.debugging,_self)):
                        traceback.print_exc()
                        print "Error: unable to read setting value."
                    raise QuietException
            finally:
                _self.unlock(r,_self)
        if _self._raising(r,_self): raise QuietException            
        return r

        
    def set(name, value=1, selection='', state=0, updates=1, log=0,
            quiet=1,_self=cmd):
        
        '''
DESCRIPTION

    "set" changes global, object, object-state, or per-atom settings.

USAGE

    set name [,value [,selection [,state ]]]

ARGUMENTS

    name = string: setting name

    value = string: a setting value {default: 1}

    selection = string: name-pattern or selection-expression
    {default:'' (global)}

    state = a state number {default: 0 (per-object setting)}

EXAMPLES

    set orthoscopic

    set line_width, 3

    set surface_color, white, 1hpv
    
    set sphere_scale, 0.5, elem C

NOTES

    The default behavior (with a blank selection) is global.  If the
    selection is "all", then the setting entry in each individual
    object will be changed.  Likewise, for a given object, if state is
    zero, then the object setting will be modified.  Otherwise, the
    setting for the indicated state within the object will be
    modified.

    If a selection is provided as opposed to an object name, then the
    atomic setting entries are modified.

    The following per-atom settings are currently implemented.  Others
    may seem to be recognized but will have no effect when set on a
    per-atom basis.
    
    * sphere_color
    * surface_color
    * mesh_color
    * label_color
    * dot_color
    * cartoon_color
    * ribbon_color
    * transparency (for surfaces)
    * sphere_transparency
    
    Note that if you attempt to use the "set" command with a per-bond
    setting over a selection of atoms, the setting change will appear
    to take, but no change will be observed.  Please use the
    "set_bond" command for per-bond settings.
    

PYMOL API

    cmd.set(string name, string value, string selection, int state,
            int updates, int quiet)

SEE ALSO

    get, set_bond
    
'''
        r = DEFAULT_ERROR
        selection = str(selection)
        if log:
            if ',' in str(value):
                value = str(cmd.safe_list_eval(str(value)))
            if len(selection):
                cmd.log("set %s,%s,%s\n"%(str(name),str(value),str(selection)))
            else:
                cmd.log("set %s,%s\n"%(str(name),str(value)))            
        index = _get_index(str(name))
        if(index<0):
            print "Error: unknown setting '%s'."%name
            raise QuietException
        else:
            try:
                _self.lock(_self)
                type = _cmd.get_setting_tuple(_self._COb,int(index),str(""),int(-1))[0]
                if type==None:
                    print "Error: unable to get setting type."
                    raise QuietException
                try:
                    if type==1: # boolean (also support non-zero float for truth)
                        handled = 0
                        if boolean_sc.interpret(str(value))==None:
                            try: # number, non-zero, then interpret as TRUE
                                if not (float(value)==0.0):
                                    handled = 1
                                    v = (1,)
                                else:
                                    handled = 1
                                    v = (0,)
                            except:
                                pass
                        if not handled:
                            v = (boolean_dict[
                                boolean_sc.auto_err(
                                str(value),"boolean")],)
                    elif type==2: # int (also supports boolean language for 0,1)
                        if boolean_sc.has_key(str(value)):
                            v = (boolean_dict[
                                boolean_sc.auto_err(
                                str(value),"boolean")],)
                        else:
                            v = (int(value),)
                    elif type==3: # float
                        if boolean_sc.has_key(str(value)):
                            v = (float(boolean_dict[
                                boolean_sc.auto_err(
                                str(value),"boolean")]),)
                        else:
                            v = (float(value),)
                    elif type==4: # float3 - some legacy handling req.
                        if is_string(value):
                            if not ',' in value:
                                v = string.split(value)
                            else:
                                v = eval(value)
                        else:
                            v = value
                        v = (float(v[0]),float(v[1]),float(v[2]))
                    elif type==5: # color
                        v = (str(value),)
                    elif type==6: # string
                        vl = str(value)
                        # strip outermost quotes (cheesy approach)
                        if quote_strip_re.search(vl)!=None:
                            vl=vl[1:-1]
                        v = (vl,)
                    v = (type,v)
                    if len(selection):
                        selection=selector.process(selection)
                    r = _cmd.set(_self._COb,int(index),v,
                                     selection,
                                     int(state)-1,int(quiet),
                                     int(updates))
                except:
                    traceback.print_exc()
                    if(_feedback(fb_module.cmd,fb_mask.debugging,_self)):
                        traceback.print_exc()
                        print "Error: unable to read setting value."
                    raise QuietException
            finally:
                _self.unlock(r,_self)
        if _self._raising(r,_self): raise QuietException            
        return r

    def unset(name, selection='', state=0, updates=1, log=0, quiet=1, _self=cmd):
        '''
DESCRIPTION

    "unset" clear non-global settings and zeros out global settings.

USAGE

    unset name [,selection [,state ]]

EXAMPLE

    unset orthoscopic

    unset surface_color, 1hpv

    unset sphere_scale, elem C
    
NOTES

    If selection is not provided, unset changes the named global
    setting to a zero or off value.

    If a selection is provided, then "unset" undefines per-object,
    per-state, or per-atom settings.

PYMOL API

    cmd.unset(string name, string selection, int state, int updates,
                int log)

SEE ALSO

    set, set_bond
    
        '''
        r = DEFAULT_ERROR
        selection = str(selection)
        if log:
            if(len(selection)):
                cmd.log("unset %s,%s\n"%(str(name),str(selection)))
            else:
                cmd.log("set %s,0\n"%(str(name)))
        index = _get_index(str(name))
        if(index<0):
            print "Error: unknown setting '%s'."%name
            raise QuietException
        else:
            if not len(selection):
                r = set(name,0,'',state,updates,log=0,quiet=quiet,_self=_self)
            else:
                try:
                    _self.lock(_self)
                    try:
                        selection = selector.process(selection)   
                        r = _cmd.unset(_self._COb,int(index),selection,
                                            int(state)-1,int(quiet),
                                            int(updates))
                    except:
                        if(_feedback(fb_module.cmd,fb_mask.debugging,_self)):
                            traceback.print_exc()
                            raise QuietException
                        print "Error: unable to unset setting value."
                finally:
                    _self.unlock(r,_self)
        if _self._raising(r,_self): raise QuietException            
        return r
    
    def unset_bond(name,selection1,selection2=None,state=0,updates=1,log=0,quiet=1,_self=cmd):
        '''
DESCRIPTION

    "unset_bond" removes a per-bond setting for a given set of bonds.
    
USAGE

    unset name [,selection [, selection [,state ]]]

        '''
        r = DEFAULT_ERROR
        selection1 = str(selection1)
        if selection2 == None:
            selection2 = selection1
        selection2 = str(selection2)
        if log:
            if(len(selection2)):
                cmd.log("unset %s,%s\n"%(str(name),str(selection1),str(selection2)))
            else:
                cmd.log("set %s,%s\n"%(str(name),str(selection1)))
        index = _get_index(str(name))
        if(index<0):
            print "Error: unknown setting '%s'."%name
            raise QuietException
        else:
            try:
                _self.lock(_self)
                try:
                    selection1 = selector.process(selection1)
                    selection2 = selector.process(selection2)   
                    r = _cmd.unset_bond(_self._COb,int(index),selection1,selection2,
                                   int(state)-1,int(quiet),
                                   int(updates))
                except:
                    if(_feedback(fb_module.cmd,fb_mask.debugging,_self)):
                        traceback.print_exc()
                        raise QuietException
                    print "Error: unable to unset setting value."
            finally:
                _self.unlock(r,_self)
        if _self._raising(r,_self): raise QuietException            
        return r

    def get_setting(name,object='',state=0,_self=cmd): # INTERNAL
        r = DEFAULT_ERROR
        if is_string(name):
            i = _get_index(name)
        else:
            i = int(name)
        if i<0:
            print "Error: unknown setting"
            raise QuietException
        try:
            _self.lock(_self)
            r = _cmd.get_setting_tuple(_self._COb,i,str(object),int(state)-1)
            typ = r[0]
            if typ<3: # boolean, int
                value = int(r[1][0])
            elif typ<4: # float
                value = r[1][0]
            elif typ<5: # vector
                value = r[1]
            else:
                value = r[1] # color or string
        finally:
            _self.unlock(r,_self)
        if is_ok(r):
            return value
        elif _self._raising(r,_self):
            raise QuietException                     
        return r

    def get(name, selection='', state=0, quiet=1, _self=cmd):
        '''
DESCRIPTION

    "get" prints out the current value of a setting.

USAGE

    get name [, selection [, state ]]
    
EXAMPLE

    get line_width

ARGUMENTS

    name = string: setting name

    selection = string: object name (selections not yet supported)

    state = integer: state number
    
NOTES

    "get" currently only works with global, per-object, and per-state
    settings.  There is currently no way to retrieve per-atom settings.
    
PYMOL API

    cmd.get(string name, string object, int state, int quiet)

SEE ALSO

    set, set_bond, get_bond

    '''
        
        r = DEFAULT_ERROR
        state = int(state)
        if is_string(name):
            i = _get_index(name)
        else:
            i = int(name)
        if i<0:
            print "Error: unknown setting"
            raise QuietException
        try:
            _self.lock(_self)
            r = _cmd.get_setting_text(_self._COb,i,str(selection),state-1)
        finally:
            _self.unlock(r,_self)
        if is_ok(r) and (r!=None):
            r_str = str(r)
            if len(string.strip(r_str))==0:
                r_str = "\'"+r_str+"\'"
            if not quiet:
                if(selection==''):
                    print " get: %s = %s"%(name,r_str)
                elif state<=0:
                    print " get: %s = %s in object %s"%(name,r_str,selection)
                else:
                    print " get: %s = %s in object %s state %d"%(name,r_str,selection,state)
        if _self._raising(r,_self): raise QuietException
        return r
    
    def get_setting_tuple(name,object='',state=0,_self=cmd): # INTERNAL
        r = DEFAULT_ERROR
        if is_string(name):
            i = _get_index(name)
        else:
            i = int(name)
        if i<0:
            print "Error: unknown setting"
            raise QuietException
        try:
            _self.lock(_self)
            r = _cmd.get_setting_tuple(_self._COb,i,str(object),int(state)-1)
        finally:
            _self.unlock(r,_self)
        if _self._raising(r,_self): raise QuietException
        return r
    
    def get_setting_boolean(name,object='',state=0,_self=cmd): # INTERNAL
        r = DEFAULT_ERROR
        if is_string(name):
            i = _get_index(name)
        else:
            i = int(name)
        if i<0:
            print "Error: unknown setting"
            raise QuietException
        try:
            _self.lock(_self)
            r = _cmd.get_setting_of_type(_self._COb,i,str(object),int(state)-1,1)
        finally:
            _self.unlock(r,_self)
        if _self._raising(r,_self): raise QuietException
        return r
    
    def get_setting_int(name,object='',state=0,_self=cmd): # INTERNAL
        r = DEFAULT_ERROR
        if is_string(name):
            i = _get_index(name)
        else:
            i = int(name)
        if i<0:
            print "Error: unknown setting"
            raise QuietException
        try:
            _self.lock(_self)
            r = _cmd.get_setting_of_type(_self._COb,i,str(object),int(state)-1,2)
        finally:
            _self.unlock(r,_self)
        if _self._raising(r,_self): raise QuietException
        return r
    
    def get_setting_float(name,object='',state=0,_self=cmd): # INTERNAL
        r = DEFAULT_ERROR
        if is_string(name):
            i = _get_index(name)
        else:
            i = int(name)
        if i<0:
            print "Error: unknown setting"
            raise QuietException
        try:
            _self.lock(_self)
            r = _cmd.get_setting_of_type(_self._COb,i,str(object),int(state)-1,3)
        finally:
            _self.unlock(r,_self)
        if _self._raising(r,_self): raise QuietException
        return r

    def get_setting_text(name,object='',state=0,_self=cmd):  # INTERNAL
        r = DEFAULT_ERROR
        if is_string(name):
            i = _get_index(name)
        else:
            i = int(name)
        if i<0:
            print "Error: unknown setting"
            raise QuietException
        try:
            _self.lock(_self)
            r = _cmd.get_setting_text(_self._COb,i,str(object),int(state)-1)
        finally:
            _self.unlock(r,_self)
        if _self._raising(r,_self): raise QuietException
        return r

    def get_setting_updates(_self=cmd): # INTERNAL
        r = []
        if lock_attempt(_self):
            try:
                r = _cmd.get_setting_updates(_self._COb)
            finally:
                _self.unlock(r,_self)
        return r

    def get_setting_legacy(name,_self=cmd): # INTERNAL, DEPRECATED
        r = DEFAULT_ERROR
        try:
            _self.lock(_self)
            r = _cmd.get_setting(_self._COb,name)
        finally:
            _self.unlock(r,_self)
        if _self._raising(r,_self): raise QuietException
        return r

