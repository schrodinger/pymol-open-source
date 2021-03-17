cmd = __import__("sys").modules["pymol.cmd"]

class ExprShortcut(cmd.Shortcut):
    '''
    Expression shortcut for iterate/alter/label with "s." prefix
    setting autocompletion.
    '''
    def interpret(self, kee, mode=0):
        if not kee.startswith('s.'):
            return cmd.Shortcut.interpret(self, kee, mode)
        v = cmd.setting.setting_sc.interpret(kee[2:])
        if isinstance(v, str):
            return 's.' + v
        if isinstance(v, list):
            return ['s.' + v for v in v]
        return None

expr_sc = ExprShortcut([
    'segi', 'chain', 'resn', 'resi', 'name', 'alt', 'elem', 'text_type',
    'formal_charge', 'numeric_type', 'ID',
    'q', 'b', 'partial_charge', 'vdw',
    'p.', 's.',
    'elec_radius',
    'oneletter',
    'model', 'resv', 'type', 'stereo', 'rank', 'index', 'ss', 'color', 'reps',
    'protons', 'label', 'geom', 'valence', 'flags', 'cartoon',
])


def fragments_sc():
    import os
    import chempy
    return cmd.Shortcut([
        f[:-4] for f in os.listdir(chempy.path + 'fragments')
        if f.endswith('.pkl')
    ])


def vol_ramp_sc():
    from . import colorramping
    return cmd.Shortcut(colorramping.namedramps)

names_sc = lambda: cmd.Shortcut(cmd.get_names('public'))

aa_nam_e = [ names_sc                   , 'name'            , ''   ]
aa_nam_s = [ names_sc                   , 'name'            , ' '  ]
aa_nam_c = [ names_sc                   , 'name'            , ', ' ]
aa_exp_e = [ expr_sc                    , 'expression'      , ''   ]
aa_sel_e = [ cmd.selection_sc           , 'selection'       , ''   ]
aa_sel_c = [ cmd.selection_sc           , 'selection'       , ', ' ]
aa_obj_e = [ cmd.object_sc              , 'object'          , ''   ]
aa_obj_s = [ cmd.object_sc              , 'object'          , ' '  ]
aa_obj_c = [ cmd.object_sc              , 'object'          , ', ' ]
aa_set_c = [ cmd.setting.setting_sc     , 'setting'         , ', ' ]
aa_map_c = [ cmd.map_sc                 , 'map object'      , ', ' ]
aa_rep_c = [ cmd.repres_sc              , 'representation'  , ', ' ]
aa_rem_c = [ cmd.repmasks_sc            , 'representation'  , ', ' ]
aa_v_r_c = [ vol_ramp_sc                , 'volume ramp'     , ', ' ]
aa_ali_e = [ cmd.Shortcut(['align', 'super', 'cealign']), 'alignment method', '']

def wizard_sc():
    import os, pymol.wizard
    names_glob = [name[:-3] for p in pymol.wizard.__path__
            for name in os.listdir(p) if name.endswith('.py')]
    return cmd.Shortcut(names_glob)

def get_auto_arg_list(self_cmd=cmd):
    self_cmd = self_cmd._weakrefproxy

    aa_vol_c = [ lambda:
            cmd.Shortcut(self_cmd.get_names_of_type('object:volume')),
            'volume', '' ]
    aa_ramp_c = [ lambda:
            cmd.Shortcut(self_cmd.get_names_of_type('object:ramp')),
            'ramp', '' ]
    aa_scene_e = [lambda: cmd.Shortcut(cmd.get_scene_list()), 'scene', '']

    return [
# 1st
        {
        'align'          : aa_sel_c,
        'alignto'        : aa_obj_c,
        'alter'          : aa_sel_e,
        'alphatoall'     : aa_sel_c,
        'api'            : [ self_cmd.kwhash, 'command', '' ],
        'assign_stereo'  : aa_sel_e,
        'bond'           : aa_sel_e,
        'as'             : aa_rem_c,
        'bg_color'       : [ lambda c=self_cmd:c._get_color_sc(c), 'color'       , ''   ],
        'button'         : [ self_cmd.controlling.button_sc  , 'button'          , ', ' ],
        'cartoon'        : [ self_cmd.viewing.cartoon_sc     , 'cartoon'         , ', ' ],
        'cache'          : [ self_cmd.exporting.cache_action_sc , 'cache mode'   , ', ' ],
        'center'         : aa_sel_e,
        'cealign'        : aa_sel_e,
        'centerofmass'   : aa_sel_e,
        'color'          : [ lambda c=self_cmd:c._get_color_sc(c), 'color'       , ', ' ],
        'color_deep'     : [ lambda c=self_cmd:c._get_color_sc(c), 'color'       , ', ' ],
        'config_mouse'   : [ self_cmd.controlling.ring_dict_sc, 'mouse cycle'    , ''   ],
        'clean'          : aa_sel_c,
        'clip'           : [ self_cmd.viewing.clip_action_sc , 'clipping action' , ', ' ],
        'copy'           : aa_obj_c,
        'copy_to'        : aa_obj_c,
        'count_atoms'    : aa_sel_e,
        'count_discrete' : aa_sel_e,
        'create'         : aa_obj_c,
        'delete'         : aa_nam_s,
        'deprotect'      : aa_sel_e,
        'disable'        : aa_obj_s,
        'distance'       : aa_obj_e,
        'dss'            : aa_sel_e,
        'enable'         : aa_obj_s,
        'extra_fit'      : aa_sel_e,
        'extract'        : aa_obj_e,
        'feedback'       : [ self_cmd.fb_action_sc           , 'action'          , ', ' ],
        'fit'            : aa_sel_e,
        'flag'           : [ self_cmd.editing.flag_sc        , 'flag'            , ', ' ],
        'fragment'       : [ fragments_sc                    , 'fragment name'   , ''   ],
        'full_screen'    : [ self_cmd.toggle_sc              , 'option'          , ''   ],
        'fuse'           : aa_sel_e,
        'get'            : aa_set_c,
        'get_area'       : aa_sel_e,
        'get_bond'       : aa_set_c,
        'get_chains'     : aa_sel_e,
        'get_extent'     : aa_sel_e,
        'get_property_list' : aa_obj_c,
        'get_symmetry'   : aa_obj_c,
        'gradient'       : [ self_cmd.object_sc              , 'gradient'        , ', ' ],
        'group'          : [ self_cmd.group_sc               , 'group object'    , ', ' ],
        'help'           : [ self_cmd.help_sc                , 'selection'       , ''   ],
        'help_setting'   : [ self_cmd.setting.setting_sc     , 'setting'         , ''   ],
        'h_add'          : aa_sel_e,
        'hide'           : aa_rem_c,
        'isolevel'       : [ self_cmd.contour_sc             , 'contour'         , ', ' ],
        'iterate'        : aa_sel_e,
        'indicate'       : aa_sel_e,
        'intra_fit'      : aa_sel_e,
        'label'          : aa_sel_e,
        'map_set'        : aa_map_c,
        'mask'           : aa_sel_e,
        'mview'          : [ self_cmd.moving.mview_action_sc , 'action'          , ''   ],
        'map_double'     : aa_map_c,
        'map_halve'      : aa_map_c,
        'map_trim'       : aa_map_c,
        'matrix_copy'    : aa_obj_c,
        'matrix_reset'   : aa_obj_c,
        'mse2met'        : aa_sel_e,
        'order'          : aa_nam_s,
        'orient'         : aa_sel_e,
        'origin'         : aa_sel_e,
        'pair_fit'       : aa_sel_c,
        'pbc_unwrap'     : aa_obj_e,
        'pbc_wrap'       : aa_obj_e,
        'protect'        : aa_sel_e,
        'pi_interactions': aa_obj_e,
        'pseudoatom'     : aa_obj_c,
        'ramp_new'       : aa_ramp_c,
        'ramp_update'    : aa_ramp_c,
        'rebond'         : aa_obj_e,
        'rebuild'        : aa_sel_e,
        'reference'      : [ self_cmd.editing.ref_action_sc  , 'action'          , ', ' ],
        'remove'         : aa_sel_e,
        'reinitialize'   : [ self_cmd.commanding.reinit_sc   , 'option'          , ''   ],
        'scene'          : aa_scene_e,
        'sculpt_activate': aa_obj_e,
        'sculpt_deactivate': aa_obj_e,
        'sculpt_iterate' : aa_obj_c,
        'set'            : aa_set_c,
        'set_bond'       : aa_set_c,
        'set_key'        : [ lambda: cmd.Shortcut(cmd.key_mappings), 'key'       , ', ' ],
        'set_name'       : aa_nam_c,
        'set_title'      : aa_obj_c,
        'show'           : aa_rem_c,
        'smooth'         : aa_sel_e,
        'space'          : [ self_cmd.space_sc               , 'space'           , ''   ],
        'spectrum'       : aa_exp_e,
        'split_chains'   : aa_sel_e,
        'split_states'   : aa_obj_c,
        'super'          : aa_sel_c,
        'stereo'         : [ self_cmd.stereo_sc              , 'option'          , ''   ],
        'symmetry_copy'  : aa_obj_c,
        'toggle'         : aa_rem_c,
        'uniquify'       : [ expr_sc                         , 'identifier'      , ', ' ],
        'unmask'         : aa_sel_e,
        'unset'          : aa_set_c,
        'unset_bond'     : aa_set_c,
        'unset_deep'     : aa_set_c,
        'update'         : aa_sel_e,
        'valence'        : [ self_cmd.editing.order_sc       , 'order'           , ', ' ],
        'volume_color'   : aa_vol_c,
        'volume_panel'   : aa_vol_c,
        'view'           : [ self_cmd._pymol._view_dict_sc   , 'view'            , ''   ],
        'window'         : [ self_cmd.window_sc              , 'action'          , ', ' ],
        'wizard'         : [ wizard_sc                       , 'wizard'          , ', '   ],
        'zoom'           : aa_sel_e,
        },
# 2nd
        {
        'align'          : aa_sel_e,
        'alignto'        : aa_ali_e,
        'alter'          : aa_exp_e,
        'alter_state'    : aa_sel_e,
        'alphatoall'     : aa_exp_e,
        'as'             : aa_sel_e,
        'bond'           : aa_sel_e,
        'button'         : [ self_cmd.controlling.but_mod_sc , 'modifier'        , ', ' ],
        'cache'          : aa_scene_e,
        'cealign'        : aa_sel_e,
        'clean'          : aa_sel_e,
        'color'          : aa_sel_e,
        'color_deep'     : aa_obj_e,
        'copy'           : aa_obj_e,
        'copy_to'        : aa_sel_e,
        'create'         : aa_sel_c,
        'distance'       : aa_sel_e,
        'extra_fit'      : aa_obj_e,
        'extract'        : aa_sel_e,
        'feedback'       : [ self_cmd.fb_module_sc           , 'module'          , ', ' ],
        'flag'           : aa_sel_c,
        'fuse'           : aa_sel_e,
        'get'            : aa_obj_c,
        'get_bond'       : aa_sel_c,
        'get_property'   : aa_obj_c,
        'gradient'       : aa_map_c,
        'group'          : aa_nam_s,
        'hide'           : aa_sel_e,
        'isomesh'        : aa_map_c,
        'isosurface'     : aa_map_c,
        'iterate'        : aa_exp_e,
        'iterate_state'  : aa_sel_e,
        'join_states'    : aa_sel_e,
        'volume'         : aa_map_c,
        'select'         : aa_sel_e,
        'save'           : aa_sel_c,
        'label'          : aa_exp_e,
        'load'           : aa_sel_c,
        'load_traj'      : aa_obj_c,
        'map_set'        : [ self_cmd.editing.map_op_sc      , 'operator'        , ', ' ],
        'map_new'        : [ self_cmd.creating.map_type_sc   , 'map type'        , ', ' ],
        'map_trim'       : aa_sel_c,
        'morph'          : aa_sel_e,
        'movie.produce'  : [ self_cmd.movie.produce_mode_sc  , 'render mode'     , ', ' ],
        'matrix_copy'    : aa_obj_c,
        'multifilesave'  : aa_sel_c,
        'order'          : [ self_cmd.boolean_sc             , 'sort'            , ', ' ],
        'pair_fit'       : aa_sel_c,
        'pi_interactions': aa_sel_e,
        'rebuild'        : aa_rep_c,
        'reference'      : aa_sel_c,
        'scene'          : [ self_cmd.viewing.scene_action_sc, 'scene action'    , ', ' ],
        'set_name'       : aa_nam_e,
        'show'           : aa_sel_e,
        'slice_new'      : aa_map_c,
        'spectrum'       : [ self_cmd.palette_sc             , 'palette'         , ''   ],
        'super'          : aa_sel_e,
        'symexp'         : aa_obj_c,
        'symmetry_copy'  : aa_obj_c,
        'toggle'         : aa_sel_e,
        'view'           : [ self_cmd.viewing.view_sc        , 'view action'     , ''   ],
        'uniquify'       : aa_sel_e,
        'unset'          : aa_sel_c,
        'unset_bond'     : aa_sel_c,
        'unset_deep'     : aa_obj_e,
        'update'         : aa_sel_e,
        'ramp_new'       : aa_map_c,
        'valence'        : aa_sel_c,
        'volume_color'   : aa_v_r_c,
        },
#3rd
        {
        'alter_state'    : aa_exp_e,
        'button'         : [ self_cmd.controlling.but_act_sc , 'button action'   , ''   ],
        'callout'        : aa_sel_e,
        'distance'       : aa_sel_e,
        'extra_fit'      : aa_ali_e,
        'feedback'       : [ self_cmd.fb_mask_sc             , 'mask'            , ''   ],
        'flag'           : [ self_cmd.editing.flag_action_sc , 'flag action'     , ''   ],
        'get_bond'       : aa_sel_e,
        'group'          : [ self_cmd.creating.group_action_sc, 'group action'    , ''   ],
        'iterate_state'  : aa_exp_e,
        'map_set'        : [ self_cmd.map_sc                 , 'map'             , ' '  ],
        'morph'          : aa_sel_e,
        'order'          : [ self_cmd.controlling.location_sc, 'location'        , ', ' ],
        'pi_interactions': aa_sel_e,
        'ramp_update'    : [ self_cmd.creating.ramp_spectrum_sc , 'ramp color spectrum' , ', ' ],
        'set'            : aa_sel_c,
        'set_bond'       : aa_sel_c,
        'set_property'   : aa_obj_c,
        'set_atom_property' : aa_sel_e,
        'spectrum'       : aa_sel_e,
        'symexp'         : aa_sel_c,
        'uniquify'       : aa_sel_e,
        'unset_bond'     : aa_sel_c,
        'valence'        : aa_sel_c,
        'volume'         : aa_v_r_c,
        },
#4th
        {
        'ramp_new'       : [ self_cmd.creating.ramp_spectrum_sc , 'ramp color spectrum' , ', ' ],
        'map_new'        : aa_sel_c,
        'isosurface'     : aa_sel_c,
        'volume'         : aa_sel_c,
        'isomesh'        : aa_sel_c,
        'set_bond'       : aa_sel_c,
        'valence'        : aa_sel_c,
        }
        ]
