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
        _feedback,fb_module,fb_mask 
   from cmd import is_string

   boolean_type = 1
   int_type     = 2
   float_type   = 3
   float3_type  = 4

   class SettingIndex:
      bonding_vdw_cutoff    =0
      min_mesh_spacing      =1
      dot_density           =2
      dot_mode              =3
      solvent_radius        =4
      sel_counter           =5
      bg_rgb                =6
      ambient               =7
      direct                =8
      reflect               =9
      light                =10
      power                =11
      antialias            =12
      cavity_cull          =13
      ambient_scale        =14
      single_image         =15
      movie_delay          =16
      ribbon_power         =17
      ribbon_power_b       =18
      ribbon_sampling      =19
      ribbon_radius        =20
      stick_radius         =21
      hash_max             =22
      orthoscopic          =23  # note, ortho still works...
      spec_reflect         =24
      spec_power           =25
      sweep_angle          =26
      sweep_speed          =27
      dot_hydrogens        =28
      dot_radius           =29
      ray_trace_frames     =30
      cache_frames         =31
      trim_dots            =32
      cull_spheres         =33
      test1                =34
      test2                =35
      surface_best         =36
      surface_normal       =37
      surface_quality      =38
      surface_proximity    =39
      normal_workaround    =40 
      stereo_angle         =41 
      stereo_shift         =42
      line_smooth          =43
      line_width           =44
      half_bonds           =45
      stick_quality        =46
      stick_overlap        =47
      stick_nub            =48
      all_states           =49
      pickable             =50
      auto_show_lines      =51
      idle_delay           =52
      no_idle              =53
      fast_idle            =54
      slow_idle            =55
      rock_delay           =56
      dist_counter         =57
      dash_length          =58
      dash_gap             =59
      auto_zoom            =60
      overlay              =61
      text                 =62
      button_mode          =63 
      valence              =64
      nonbonded_size       =65
      label_color          =66
      ray_trace_fog        =67
      spheroid_scale       =68
      ray_trace_fog_start  =69
      spheroid_smooth      =70
      spheroid_fill        =71
      auto_show_nonbonded  =72
      cache_display        =73
      mesh_radius          =74
      backface_cull        =75
      gamma                =76
      dot_width            =77
      auto_show_selections =78
      auto_hide_selections =79
      selection_width      =80
      selection_overlay    =81
      static_singletons    =82
      max_triangles        =83
      depth_cue            =84
      specular             =85
      shininess            =86
      sphere_quality       =87
      fog                  =88
      isomesh_auto_state   =89
      mesh_width           =90
      cartoon_sampling     =91
      cartoon_loop_radius  =92
      cartoon_loop_quality =93
      cartoon_power        =94
      cartoon_power_b      =95
      cartoon_rect_length  =96
      cartoon_rect_width   =97
      internal_gui_width   =98
      internal_gui         =99
      cartoon_oval_length  =100
      cartoon_oval_width   =101
      cartoon_oval_quality =102
      cartoon_tube_radius  =103
      cartoon_tube_quality =104
      cartoon_debug        =105
      ribbon_width         =106
      dash_width           =107
      dash_radius          =108
      cgo_ray_width_scale  =109
      line_radius          =110
      cartoon_round_helices     =111
      cartoon_refine_normals    =112
      cartoon_flat_sheets       =113
      cartoon_smooth_loops      =114
      cartoon_dumbbell_length   =115
      cartoon_dumbbell_width    =116
      cartoon_dumbbell_radius   =117
      cartoon_fancy_helices     =118
      cartoon_fancy_sheets      =119
      ignore_pdb_segi       =120
      ribbon_throw          =121
      cartoon_throw         =122
      cartoon_refine        =123
      cartoon_refine_tips   =124
      cartoon_discrete_colors   =125
      normalize_ccp4_maps   =126
      surface_poor          =127
      internal_feedback     =128
      cgo_line_width        =129
      cgo_line_radius       =130
      logging               =131
      robust_logs           =132
      log_box_selections    =133
      log_conformations     =134
      valence_size          =135
      valence_default       =135 # legacy
      surface_miserable     =136
      ray_opaque_background =137
      transparency          =138
      ray_texture           =139
      ray_texture_settings  =140
      suspend_updates       =141
      full_screen           =142
      surface_mode          =143
      surface_color         =144
      mesh_mode             =145
      mesh_color            =146
      auto_indicate_flags   =147
      surface_debug         =148
      ray_improve_shadows   =149
      smooth_color_triangle =150
      ray_default_renderer  =151
      field_of_view         =152
      reflect_power         =153
      preserve_chempy_ids   =154
      sphere_scale          =155
      two_sided_lighting    =156
      secondary_structure   =157
      auto_remove_hydrogens =158
      raise_exceptions      =159
      stop_on_exceptions    =160
      sculpting             =161
      auto_sculpt           =162
      sculpt_vdw_scale      =163
      sculpt_vdw_scale14    =164
      sculpt_vdw_weight     =165
      sculpt_vdw_weight14   =166
      sculpt_bond_weight    =167
      sculpt_angl_weight    =168
      sculpt_pyra_weight    =169
      sculpt_plan_weight    =170
      sculpting_cycles      =171
      sphere_transparency   =172
      sphere_color          =173
      sculpt_field_mask     =174
      sculpt_hb_overlap     =175
      sculpt_hb_overlap_base=176
      legacy_vdw_radii      =177
      sculpt_memory         =178
      connect_mode          =179
      cartoon_cylindrical_helices = 180
      cartoon_helix_radius  =181
      connect_cutoff        =182
   #   save_pdb_ss           =183
      sculpt_line_weight    =184
      fit_iterations        =185
      fit_tolerance         =186
      batch_prefix          =187
      stereo_mode           =188
      cgo_sphere_quality    =189
      pdb_literal_names     =190
      wrap_output           =191
      fog_start             =192
      state                 =193
      frame                 =194
      ray_shadows           =195 # legacy
      ray_shadow            =195 # improved ease of use..(for ray_sh <Tab> )
      ribbon_trace          =196
      security              =197
      stick_transparency    =198
      ray_transparency_shadows = 199
      session_version_check = 200
      ray_transparency_specular = 201
      stereo_double_pump_mono = 202
      sphere_solvent        = 203
      mesh_quality          = 204
      mesh_solvent          = 205
      dot_solvent           = 206
      ray_shadow_fudge      = 207
      ray_triangle_fudge    = 208
      debug_pick            = 209
      dot_color             = 210
      mouse_limit           = 211
      mouse_scale           = 212
      transparency_mode     = 213
      clamp_colors          = 214
      pymol_space_max_red   = 215
      pymol_space_max_green = 216
      pymol_space_max_blue  = 217
      pymol_space_min_factor= 218
      roving_origin         = 219
      roving_lines          = 220
      roving_sticks         = 221
      roving_spheres        = 222
      roving_labels         = 223
      roving_delay          = 224
      roving_selection      = 225
      roving_byres          = 226
      roving_ribbon         = 227
      roving_cartoon        = 228
      roving_polar_contacts = 229
      roving_polar_cutoff   = 230
      roving_nonbonded      = 231
      float_labels          = 232
      roving_detail         = 233
      roving_nb_spheres     = 234
      ribbon_color          = 235
      cartoon_color         = 236
      ribbon_smooth         = 237
      auto_color            = 238
      ray_interior_color    = 240
      cartoon_highlight_color = 241
      coulomb_units_factor  = 242
      coulomb_dielectric    = 243
      ray_interior_shadows  = 244
      ray_interior_texture  = 245

      roving_map1_name      = 246
      roving_map2_name      = 247
      roving_map3_name      = 248

      roving_map1_level     = 249
      roving_map2_level     = 250
      roving_map3_level     = 251

      roving_isomesh        = 252
      roving_isosurface     = 253
      scenes_changed        = 254
      gaussian_b_adjust     = 255
      pdb_standard_order    = 256
      cartoon_smooth_first  = 257
      cartoon_smooth_last   = 258
      cartoon_smooth_cycles = 259
      cartoon_flat_cycles   = 260
      max_threads           = 261
      show_progress         = 262
      use_display_lists     = 263
      cache_memory          = 264
      simplify_display_lists= 265
      retain_order          = 266
      pdb_hetatm_sort       = 267
      pdb_use_ter_records   = 268
      cartoon_trace         = 269
      ray_oversample_cutoff = 270
      gaussian_resolution   = 271
      gaussian_b_floor      = 272
      sculpt_nb_interval    = 273
      sculpt_tors_weight    = 274
      sculpt_tors_tolerance = 275
      stick_ball            = 276
      stick_ball_ratio      = 277
      stick_fixed_radius    = 278
      cartoon_transparency  = 279
      dash_round_ends       = 280
      h_bond_max_angle      = 281
      h_bond_cutoff_center  = 282
      h_bond_cutoff_edge    = 283
      h_bond_power_a        = 284
      h_bond_power_b        = 285
      h_bond_cone           = 286
      ss_helix_psi_target   = 287 
      ss_helix_psi_include  = 288
      ss_helix_psi_exclude  = 289

      ss_helix_phi_target   = 290
      ss_helix_phi_include  = 291
      ss_helix_phi_exclude  = 292

      ss_strand_psi_target  = 293
      ss_strand_psi_include = 294
      ss_strand_psi_exclude = 295

      ss_strand_phi_target  = 296
      ss_strand_phi_include = 297
      ss_strand_phi_exclude = 298
      movie_loop            = 299
      pdb_retain_ids        = 300
      pdb_no_end_record     = 301
      cgo_dot_width         = 302
      cgo_dot_radius        = 303
      defer_updates         = 304
      normalize_o_maps      = 305
      swap_dsn6_bytes       = 306
      pdb_insertions_go_first = 307
      roving_origin_z         = 308
      roving_origin_z_cushion = 309
      specular_intensity    = 310
      overlay_lines         = 311
      ray_transparency_spec_cut = 312
      internal_prompt       = 313
      normalize_grd_maps    = 314
      ray_blend_colors      = 315
      ray_blend_red         = 316
      ray_blend_green       = 317
      ray_blend_blue        = 318
      png_screen_gamma      = 319
      png_file_gamma        = 320
      editor_label_fragments = 321
      internal_gui_control_size = 322
      auto_dss              = 323
      transparency_picking_mode = 324
      virtual_trackball     = 325
      pdb_reformat_names_mode = 326
      ray_pixel_scale_to_window = 327
      label_font_id             = 328
      pdb_conect_all            = 329
      button_mode_name      = 330
      surface_type          = 331
      dot_normals           = 332
      session_migration     = 333
      mesh_normals          = 334
      mesh_type             = 335
      dot_lighting          = 336
      mesh_lighting         = 337
      surface_solvent       = 338
      triangle_max_passes   = 339
      ray_interior_reflect  = 340
      internal_gui_mode     = 341
      surface_carve_selection = 342
      surface_carve_state   = 343
      surface_carve_cutoff  = 344
      surface_clear_selection = 345
      surface_clear_state   = 346
      surface_clear_cutoff  = 347
      surface_trim_cutoff   = 348
      surface_trim_factor   = 349
      ray_max_passes        = 350
      active_selections     = 351
      ray_transparency_contrast = 352
      seq_view              = 353
      mouse_selection_mode  = 354
      seq_view_label_spacing = 355
      seq_view_label_start = 356
      seq_view_format      = 357
      seq_view_location    = 358
      seq_view_overlay     = 359
      auto_classify_atoms  = 360
      cartoon_nucleic_acid_mode = 361
      seq_view_color       = 362
      seq_view_label_mode  = 363
      surface_ramp_above_mode = 364
      stereo               = 365
      wizard_prompt_mode   = 366
      coulomb_cutoff       = 367
      slice_track_camera   = 368
      slice_height_scale   = 369
      slice_height_map     = 370
      slice_grid           = 371
      
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
      "false" : 0,
      "on"   : 1,
      "off"  : 0,
      "1"    : 1,
      "0"    : 0,
      "1.0"  : 1,
      "0.0"  : 0,
      }

   boolean_sc = Shortcut(boolean_dict.keys())


   def _get_index(name):
      # this may be called from C, so don't raise any exceptions...
      result = setting_sc.interpret(name)
      if type(result)==types.StringType:
         if hasattr(SettingIndex,result):
            return getattr(SettingIndex,result)
         else:
            return -1
      else:
         return -1

   def _get_name(index): # likewise, this can be called from C so no exceptions
      if name_dict.has_key(index):
         return name_dict[index]
      else:
         return ""

   def get_index_list():
      return index_list

   def get_name_list():
      return name_list

   ###### API functions

   def set(name,value=1,selection='',state=0,updates=1,log=0,quiet=1):
      '''
DESCRIPTION

   "set" changes one of the PyMOL state variables,

USAGE

   set name, [,value [,object-or-selection [,state ]]]

   set name = value  # (DEPRECATED)

PYMOL API

   cmd.set ( string name, string value=1,
             string selection='', int state=0,
              int updates=1, quiet=1)

NOTES

   The default behavior (with a blank selection) changes the global
   settings database.  If the selection is 'all', then the settings
   database in all individual objects will be changed.  Likewise, for
   a given object, if state is zero, then the object database will be
   modified.  Otherwise, the settings database for the indicated state
   within the object will be modified.

   If a selection is provided, then all objects in the selection will
   be affected. 

      '''
      r = None
      selection = str(selection)
      if log:
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
            lock()
            type = _cmd.get_setting_tuple(int(index),str(""),int(-1))[0]
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
                  v = (str(value),)

               v = (type,v)
               if len(selection):
                  selection=selector.process(selection)
               r = _cmd.set(int(index),v,
                            selection,
                            int(state)-1,int(quiet),
                            int(updates))
            except:
               if(_feedback(fb_module.cmd,fb_mask.debugging)):
                  traceback.print_exc()
                  print "Error: unable to read setting value."
               raise QuietException
         finally:
            unlock()
      return r

   def unset(name,selection='',state=0,updates=1,log=0,quiet=1):
      '''
DESCRIPTION

   "unset" behaves in two ways.

   If selection is not provided, unset changes the named global
   setting to a zero or off value.

   If a selection is provided, then "unset" undefines 
   object-specific or state-specific settings so that the global
   setting will be in effect.

USAGE

   unset name [,selection [,state ]]

PYMOL API

   cmd.unset ( string name, string selection = '',
            int state=0, int updates=1, int log=0 )

      '''
      r = None
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
            r = set(name,0,'',state,updates,log=0,quiet=quiet)
         else:
            try:
               lock()
               try:
                  selection = selector.process(selection)   
                  r = _cmd.unset(int(index),selection,
                                 int(state)-1,int(quiet),
                                 int(updates))
               except:
                  if(_feedback(fb_module.cmd,fb_mask.debugging)):
                     traceback.print_exc()
                     raise QuietException
                  print "Error: unable to unset setting value."
            finally:
               unlock()
      return r


   def get_setting(name,object='',state=0): # INTERNAL
      r = None
      if is_string(name):
         i = _get_index(name)
      else:
         i = int(name)
      if i<0:
         print "Error: unknown setting"
         raise QuietException
      try:
         lock()
         r = _cmd.get_setting_tuple(i,str(object),int(state)-1)
         typ = r[0]
         if typ<3: # boolean, int
            r = int(r[1][0])
         elif typ<4: # float
            r = r[1][0]
         elif typ<5: # vector
            r = r[1]
         else:
            r = r[1] # color or string
      finally:
         unlock()
      return r

   def get(name,object='',state=0,quiet=1):
      r = None
      state = int(state)
      if is_string(name):
         i = _get_index(name)
      else:
         i = int(name)
      if i<0:
         print "Error: unknown setting"
         raise QuietException
      try:
         lock()
         r = _cmd.get_setting_text(i,str(object),state-1)
      finally:
         unlock()
      if r!=None:
         if not quiet:
            if(object==''):
               print " get: %s = %s"%(name,r)
            elif state<=0:
               print " get: %s = %s in object %s"%(name,r,object)
            else:
               print " get: %s = %s in object %s state %d"%(name,r,object,state)
      return r
   
   def get_setting_tuple(name,object='',state=0): # INTERNAL
      r = None
      if is_string(name):      i = _get_index(name)
      else:
         i = int(name)
      if i<0:
         print "Error: unknown setting"
         raise QuietException
      try:
         lock()
         r = _cmd.get_setting_tuple(i,str(object),int(state)-1)
      finally:
         unlock()
      return r

   def get_setting_text(name,object='',state=0):  # INTERNAL
      r = None
      if is_string(name):
         i = _get_index(name)
      else:
         i = int(name)
      if i<0:
         print "Error: unknown setting"
         raise QuietException
      try:
         lock()
         r = _cmd.get_setting_text(i,str(object),int(state)-1)
      finally:
         unlock()
      return r

   def get_setting_updates(): # INTERNAL
      r = []
      if lock_attempt():
         try:
            r = _cmd.get_setting_updates()
         finally:
            unlock()
      return r


   def get_setting_legacy(name): # INTERNAL, DEPRECATED
      r = None
      try:
         lock()
         r = _cmd.get_setting(name)
      finally:
         unlock()
      return r
