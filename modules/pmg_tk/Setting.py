#A*
#-------------------------------------------------------------------
#B* This file contains source code for the PyMOL computer program C*
#copyright 1998-2000 by Warren Lyford Delano of DeLano Scientific.  D*
#-------------------------------------------------------------------
#E* It is unlawful to modify or remove this copyright notice.  F*
#-------------------------------------------------------------------
#G* Please see the accompanying LICENSE file for further information.
#H*
#-------------------------------------------------------------------
#I* Additional authors of this source file include: -* -* NOTE: Based
#on code by John E. Grayson which was in turn -* based on code written
#by Doug Hellmann.  Z*
#-------------------------------------------------------------------

# this section is devoted to making sure that Tkinter variables which
# correspond to Menu-displayed settings are kept synchronized with
# PyMOL

from Tkinter import *
import Pmw
from pymol import cmd
import pymol.setting
import time

pm = cmd

class Setting:

    def __init__(self):

        while not cmd.ready(): # make sure PyMOL is ready for action...
            time.sleep(0.1)
            
        self.ray_trace_frames = IntVar()
        self.ray_trace_frames.set(int(cmd.get_setting_legacy('ray_trace_frames')))
        
        self.cache_frames = IntVar()
        self.cache_frames.set(int(cmd.get_setting_legacy('ray_trace_frames')))

        self.ortho = IntVar()
        self.ortho.set(int(cmd.get_setting_legacy('orthoscopic')))

        self.antialias = IntVar()
        self.antialias.set(int(cmd.get_setting_legacy('antialias')))

        self.all_states = IntVar()
        self.all_states.set(int(cmd.get_setting_legacy('all_states')))

        self.line_smooth = IntVar()
        self.line_smooth.set(int(cmd.get_setting_legacy('line_smooth')))

        self.overlay = IntVar()
        self.overlay.set(int(cmd.get_setting_legacy('overlay')))

        self.valence = IntVar()
        self.valence.set(int(cmd.get_setting_legacy('valence')))

        self.auto_zoom = IntVar()
        self.auto_zoom.set(int(cmd.get_setting_legacy('auto_zoom')))

        self.auto_show_selections = IntVar()
        self.auto_show_selections.set(int(cmd.get_setting_legacy('auto_show_selections')))

        self.auto_hide_selections = IntVar()
        self.auto_hide_selections.set(int(cmd.get_setting_legacy('auto_hide_selections')))

        self.static_singletons = IntVar()
        self.static_singletons.set(int(cmd.get_setting_legacy('static_singletons')))

        self.backface_cull = IntVar()
        self.backface_cull.set(int(cmd.get_setting_legacy('backface_cull')))

        self.depth_cue = IntVar()
        self.depth_cue.set(int(cmd.get_setting_legacy('depth_cue')))

        self.specular = IntVar()
        self.specular.set(int(cmd.get_setting_legacy('specular')))

        self.cartoon_round_helices = IntVar()
        self.cartoon_round_helices.set(int(cmd.get_setting_legacy('cartoon_round_helices')))

        self.cartoon_fancy_helices = IntVar()
        self.cartoon_fancy_helices.set(int(cmd.get_setting_legacy('cartoon_fancy_helices')))

        self.cartoon_flat_sheets = IntVar()
        self.cartoon_flat_sheets.set(int(cmd.get_setting_legacy('cartoon_flat_sheets')))

        self.cartoon_fancy_sheets = IntVar()
        self.cartoon_fancy_sheets.set(int(cmd.get_setting_legacy('cartoon_fancy_sheets')))

        self.cartoon_discrete_colors = IntVar()
        self.cartoon_discrete_colors.set(int(cmd.get_setting_legacy('cartoon_discrete_colors')))

        self.cartoon_smooth_loops = IntVar()
        self.cartoon_smooth_loops.set(int(cmd.get_setting_legacy('cartoon_smooth_loops')))

        self.cartoon_cylindrical_helices = IntVar()
        self.cartoon_cylindrical_helices.set(int(cmd.get_setting_legacy('cartoon_cylindrical_helices')))

        self.cartoon_side_chain_helper = IntVar()
        self.cartoon_side_chain_helper.set(int(cmd.get_setting_legacy('cartoon_side_chain_helper')))

        self.ignore_pdb_segi = IntVar()
        self.ignore_pdb_segi.set(int(cmd.get_setting_legacy('ignore_pdb_segi')))

        self.log_box_selections = IntVar()
        self.log_box_selections.set(int(cmd.get_setting_legacy('log_box_selections')))

        self.log_conformations = IntVar()
        self.log_conformations.set(int(cmd.get_setting_legacy('log_conformations')))

        self.two_sided_lighting = IntVar()
        self.two_sided_lighting.set(int(cmd.get_setting_legacy('two_sided_lighting')))

        self.auto_remove_hydrogens = IntVar()
        self.auto_remove_hydrogens.set(int(cmd.get_setting_legacy('auto_remove_hydrogens')))

        self.auto_sculpt = IntVar()
        self.auto_sculpt.set(int(cmd.get_setting_legacy('auto_sculpt')))

        self.sculpting = IntVar()
        self.sculpting.set(int(cmd.get_setting_legacy('sculpting')))

        self.roving_origin = IntVar()
        self.roving_origin.set(int(cmd.get_setting_legacy('roving_origin')))

        self.roving_detail = IntVar()
        self.roving_detail.set(int(cmd.get_setting_legacy('roving_detail')))

        self.ray_interior_color = IntVar()
        self.ray_interior_color.set(int(cmd.get_setting_legacy('ray_interior_color')!=-1))

        self.cartoon_highlight_color = IntVar()
        self.cartoon_highlight_color.set(int(cmd.get_setting_legacy('ray_interior_color')!=-1))

        self.use_display_lists = IntVar()
        self.use_display_lists.set(int(cmd.get_setting_legacy('use_display_lists')!=-1))

        self.virtual_trackball = IntVar()
        self.virtual_trackball.set(int(cmd.get_setting_legacy('virtual_trackball')))

        self.stereo = IntVar()
        
        self.seq_view = IntVar()
        self.seq_view.set(int(cmd.get_setting_legacy('seq_view')))
        
        self.texture_fonts = IntVar()
        self.texture_fonts.set(int(cmd.get_setting_legacy('texture_fonts')))

        self.animation = IntVar()
        self.animation.set(int(cmd.get_setting_legacy('animation')))

        self.opaque_background = IntVar()
        self.opaque_background.set(int(cmd.get_setting_legacy('opaque_background')))

        self.show_alpha_checker = IntVar()
        self.show_alpha_checker.set(int(cmd.get_setting_legacy('show_alpha_checker')))

        self.F=[ None,
                    IntVar(),
                    IntVar(),
                    IntVar(),
                    IntVar(),
                    IntVar(),
                    IntVar(),

                    IntVar(),
                    IntVar(),
                    IntVar(),
                    IntVar(),
                    IntVar(),
                    IntVar(),
                    ]
        self.SHFTF=[ None,
                    IntVar(),
                    IntVar(),
                    IntVar(),
                    IntVar(),
                    IntVar(),
                    IntVar(),

                    IntVar(),
                    IntVar(),
                    IntVar(),
                    IntVar(),
                    IntVar(),
                    IntVar(),
                    ]
        
        self.xref = { 
            'ray_trace_frames':
            (lambda s,a: (cmd.set(a,("%1.0f" % s.ray_trace_frames.get()),log=1),
                              s.cache_frames.set(s.ray_trace_frames.get()),
                              s.update('cache_frames'))),
            'cache_frames'  :
            (lambda s,a: (cmd.set(a,("%1.0f" % s.cache_frames.get()),log=1))),
            'ortho'         :
            (lambda s,a: (cmd.set(a,("%1.0f" % s.ortho.get()),log=1))),
            'antialias'     :
            (lambda s,a: (cmd.set(a,("%1.0f" % s.antialias.get()),log=1))),
            'valence'       :
            (lambda s,a: (cmd.set(a,("%1.2f" % s.valence.get()),log=1))),
            'all_states'    :
            (lambda s,a: (cmd.set(a,("%1.0f" % s.all_states.get()),log=1))),
            'line_smooth'    :
            (lambda s,a: (cmd.set(a,("%1.0f" % s.line_smooth.get()),log=1))),
            'overlay'       :
            (lambda s,a: (cmd.set(a,("%1.0f" %( s.overlay.get()*5)),log=1))),
            'auto_zoom'       :
            (lambda s,a: (cmd.set(a,("%1.0f" % s.auto_zoom.get()),log=1))),
            'auto_show_selections'       :
            (lambda s,a: (cmd.set(a,("%1.0f" % s.auto_show_selections.get()),log=1))),
            'auto_hide_selections'       :
            (lambda s,a: (cmd.set(a,("%1.0f" % s.auto_hide_selections.get()),log=1))),
            'static_singletons'       :
            (lambda s,a: (cmd.set(a,("%1.0f" % s.static_singletons.get()),log=1))),
            'backface_cull'       :
            (lambda s,a: (cmd.set(a,("%1.0f" % s.backface_cull.get()),log=1))),
            'depth_cue'       :
            (lambda s,a: s.depth_cue_set()),
            'specular'       :
            (lambda s,a: (cmd.set(a,("%1.0f" % (s.specular.get())),log=1))),

            'cartoon_round_helices'       :
            (lambda s,a: (cmd.set(a,("%1.0f" % (s.cartoon_round_helices.get())),log=1))),
            'cartoon_fancy_helices'       :
            (lambda s,a: (cmd.set(a,("%1.0f" % (s.cartoon_fancy_helices.get())),log=1))),
            'cartoon_flat_sheets'         :
            (lambda s,a: (cmd.set(a,("%1.0f" % (s.cartoon_flat_sheets.get())),log=1))),
            'cartoon_fancy_sheets'        :
            (lambda s,a: (cmd.set(a,("%1.0f" % (s.cartoon_fancy_sheets.get())),log=1))),
            'cartoon_discrete_colors'        :
            (lambda s,a: (cmd.set(a,("%1.0f" % (s.cartoon_discrete_colors.get())),log=1))),
            'cartoon_smooth_loops'        :
            (lambda s,a: (cmd.set(a,("%1.0f" % (s.cartoon_smooth_loops.get())),log=1))),
            'cartoon_cylindrical_helices'        :
            (lambda s,a: (cmd.set(a,("%1.0f" % (s.cartoon_cylindrical_helices.get())),log=1))),
            'cartoon_side_chain_helper'        :
            (lambda s,a: (cmd.set(a,("%1.0f" % (s.cartoon_side_chain_helper.get())),log=1))),


            'ignore_pdb_segi'        :
            (lambda s,a: (cmd.set(a,("%1.0f" % (s.ignore_pdb_segi.get())),log=1))),
            'log_box_selections'        :
            (lambda s,a: (cmd.set(a,("%1.0f" % (s.log_box_selections.get())),log=1))),
            'log_conformations'        :
            (lambda s,a: (cmd.set(a,("%1.0f" % (s.log_conformations.get())),log=1))),

            'two_sided_lighting'        :
            (lambda s,a: (cmd.set(a,("%1.0f" % (s.two_sided_lighting.get())),log=1))),

            'auto_remove_hydrogens'        :
            (lambda s,a: (cmd.set(a,("%1.0f" % (s.auto_remove_hydrogens.get())),log=1))),
            'auto_sculpt'        :
            (lambda s,a: (cmd.set(a,("%1.0f" % (s.auto_sculpt.get())),log=1))),
            'sculpting'        :
            (lambda s,a: (cmd.set(a,("%1.0f" % (s.sculpting.get())),log=1))),
            'roving_origin'        :
            (lambda s,a: (cmd.set(a,("%1.0f" % (s.roving_origin.get())),log=1))),

            'roving_detail'        :
            (lambda s,a: (cmd.set(a,("%1.0f" % (s.roving_detail.get())),log=1))),
            
            'ray_interior_color'        :
            (lambda s,a: (s.ray_interior_color_set())),
            'cartoon_highlight_color'        :
            (lambda s,a: (s.cartoon_highlight_color_set())),
            'use_display_lists'        :
            (lambda s,a: (cmd.set(a,("%1.0f" % (s.use_display_lists.get())),log=1))),

            'virtual_trackball'         :
            (lambda s,a: (cmd.set(a,("%1.0f" % s.virtual_trackball.get()),log=1))),
            
            'seq_view'         :
            (lambda s,a: (cmd.set(a,("%d" % s.seq_view.get()),log=1))),
            'texture_fonts'         :
            (lambda s,a: (cmd.set(a,("%d" % s.texture_fonts.get()),log=1))),

            'animation'        :
            (lambda s,a: (cmd.set(a,("%1.0f" % (s.animation.get())),log=1))),

            'opaque_background'         :
            (lambda s,a: (cmd.set(a,("%d" % s.opaque_background.get()),log=1))),

            'show_alpha_checker'         :
            (lambda s,a: (cmd.set(a,("%d" % s.show_alpha_checker.get()),log=1))),

            }

        self.update_code = {
            'ray_trace_frames':
            (lambda s,t: (s.ray_trace_frames.set(int(t[1][0])))),
            'cache_frames':
            (lambda s,t: (s.cache_frames.set(int(t[1][0])))),
            'orthoscopic':
            (lambda s,t: (s.ortho.set(int(t[1][0])))),
            'all_states':
            (lambda s,t: (s.all_states.set(int(t[1][0])))),
            'line_smooth':
            (lambda s,t: (s.line_smooth.set(int(t[1][0])))),
            'overlay':
            (lambda s,t: (s.overlay.set(t[1][0]!=0))),
            'valence':
            (lambda s,t: (s.valence.set(t[1][0]>0.0))),
            'auto_zoom':
            (lambda s,t: (s.auto_zoom.set(int(t[1][0])!=0))), 
            'auto_show_selections':
            (lambda s,t: (s.auto_show_selections.set(int(t[1][0])))), 
            'auto_hide_selections':
            (lambda s,t: (s.auto_hide_selections.set(int(t[1][0])))), 
            'static_singletons':
            (lambda s,t: (s.static_singletons.set(int(t[1][0])))), 
            'backface_cull':
            (lambda s,t: (s.backface_cull.set(int(t[1][0])))), 
            'depth_cue'       :
            (lambda s,t: (s.depth_cue.set(int(t[1][0])))),
            'specular'       :
            (lambda s,t: (s.specular.set(t[1][0]>0.0))), 
            'cartoon_round_helices':
            (lambda s,t: (s.cartoon_round_helices.set(t[1][0]!=0))),
            'cartoon_fancy_helices':
            (lambda s,t: (s.cartoon_fancy_helices.set(t[1][0]!=0))),
            'cartoon_flat_sheets':
            (lambda s,t: (s.cartoon_flat_sheets.set(t[1][0]!=0))),
            'cartoon_fancy_sheets':
            (lambda s,t: (s.cartoon_fancy_sheets.set(t[1][0]!=0))),
            'cartoon_discrete_colors':
            (lambda s,t: (s.cartoon_discrete_colors.set(t[1][0]!=0))),
            'cartoon_smooth_loops':
            (lambda s,t: (s.cartoon_smooth_loops.set(t[1][0]!=0))),
            'cartoon_cylindrical_helices':
            (lambda s,t: (s.cartoon_cylindrical_helices.set(t[1][0]!=0))),
            'cartoon_side_chain_helper':
            (lambda s,t: (s.cartoon_side_chain_helper.set(t[1][0]!=0))),
            'ignore_pdb_segi':
            (lambda s,t: (s.ignore_pdb_segi.set(t[1][0]!=0))),
            'log_box_selections':
            (lambda s,t: (s.log_box_selections.set(t[1][0]!=0))),
            'log_conformations':
            (lambda s,t: (s.log_conformations.set(t[1][0]!=0))),
            'two_sided_lighting':
            (lambda s,t: (s.two_sided_lighting.set(t[1][0]!=0))),
            'auto_remove_hydrogens':
            (lambda s,t: (s.auto_remove_hydrogens.set(t[1][0]!=0))),
            'auto_sculpt':
            (lambda s,t: (s.auto_sculpt.set(t[1][0]!=0))),
            'sculpting':
            (lambda s,t: (s.sculpting.set(t[1][0]!=0))),
            'roving_origin':
            (lambda s,t: (s.roving_origin.set(t[1][0]!=0))),
            'roving_detail':
            (lambda s,t: (s.roving_detail.set(t[1][0]!=0))),
            'ray_interior_color':
            (lambda s,t: (s.ray_interior_color.set(t[1][0]!=-1))),
            'cartoon_highlight_color':
            (lambda s,t: (s.cartoon_highlight_color.set(t[1][0]!=-1))),
            'scenes_changed' :
            (lambda s,t: (s.update_scenes())),
            'use_display_lists':
            (lambda s,t: (s.use_display_lists.set(t[1][0]!=0))),
            'virtual_trackball':
            (lambda s,t: (s.virtual_trackball.set(t[1][0]!=0))),
            'seq_view':
            (lambda s,t: (s.seq_view.set(t[1][0]!=0))),
            'texture_fonts':
            (lambda s,t: (s.texture_fonts.set(t[1][0]!=0))),
            'stereo':
            (lambda s,t: (s.stereo.set(t[1][0]!=0))),
            'animation':
            (lambda s,t: (s.animation.set(t[1][0]!=0))),         
            'opaque_background':
            (lambda s,t: (s.animation.set(t[1][0]!=0))),
            'show_alpha_checker':
            (lambda s,t: (s.animation.set(t[1][0]!=0))),       
          }
        self.active_list = [
            pymol.setting._get_index("ray_trace_frames"),
            pymol.setting._get_index("cache_frames"),
            pymol.setting._get_index("orthoscopic"),
            pymol.setting._get_index("antialias"),
            pymol.setting._get_index("all_states"),
            pymol.setting._get_index("overlay"),
            pymol.setting._get_index("line_smooth"),
            pymol.setting._get_index("valence"),
            pymol.setting._get_index("auto_zoom"),
            pymol.setting._get_index("auto_hide_selections"),
            pymol.setting._get_index("auto_show_selections"),
            pymol.setting._get_index("static_singletons"),
            pymol.setting._get_index("backface_cull"),
            pymol.setting._get_index("depth_cue"),
            pymol.setting._get_index("specular"),
            pymol.setting._get_index("cartoon_round_helices"),
            pymol.setting._get_index("cartoon_fancy_helices"),
            pymol.setting._get_index("cartoon_flat_sheets"),
            pymol.setting._get_index("cartoon_fancy_sheets"),
            pymol.setting._get_index("cartoon_discrete_colors"),
            pymol.setting._get_index("cartoon_smooth_loops"),
            pymol.setting._get_index("cartoon_cylindrical_helices"),
            pymol.setting._get_index("cartoon_side_chain_helper"),                  
            pymol.setting._get_index("ignore_pdb_segi"),         
            pymol.setting._get_index("log_box_selections"),
            pymol.setting._get_index("log_conformations"),
            pymol.setting._get_index("two_sided_lighting"),
            pymol.setting._get_index("auto_remove_hydrogens"),
            pymol.setting._get_index("auto_sculpt"),
            pymol.setting._get_index("sculpting"),
            pymol.setting._get_index("roving_origin"),
            pymol.setting._get_index("roving_detail"),                           
            pymol.setting._get_index("ray_interior_color"),
            pymol.setting._get_index("cartoon_highlight_color"),
            pymol.setting._get_index("scenes_changed"),
            pymol.setting._get_index("use_display_lists"),
            pymol.setting._get_index("virtual_trackball"),
            pymol.setting._get_index("seq_view"),
            pymol.setting._get_index("texture_fonts"),
            pymol.setting._get_index("stereo"),
            pymol.setting._get_index("animation"),
            pymol.setting._get_index("opaque_background"),
            pymol.setting._get_index("show_alpha_checker"),            
            ]

        self.active_dict = {}
        for a in self.active_list:
            self.active_dict[a] = pymol.setting._get_name(a)

    def depth_cue_set(self):
        cmd.set("depth_cue",("%1.0f" % self.depth_cue.get()),log=1)
        if pymol.setting.get_setting("ray_trace_fog")>=0:
            cmd.set("ray_trace_fog",("%1.0f" % self.depth_cue.get()),log=1) 

    def ray_interior_color_set(self):
        if(self.ray_interior_color.get()):
            cmd.set("ray_interior_color","grey20",log=1)
        else:
            cmd.set("ray_interior_color","-1",log=1)

    def cartoon_highlight_color_set(self):
        if(self.cartoon_highlight_color.get()):
            cmd.set("cartoon_highlight_color","grey50",log=1)
        else:
            cmd.set("cartoon_highlight_color","-1",log=1)
        
    def update(self,sttng):
        set_fn = self.xref[sttng]
        set_fn(self,sttng)

    def update_scenes(self):
        dict = cmd.get_scene_dict()
        for x in range(1,13):
            if dict.has_key('F%d'%x):
                self.F[x].set(1)
            else:
                self.F[x].set(0)
            if dict.has_key('SHFT-F%d'%x):
                self.SHFTF[x].set(1)
            else:
                self.SHFTF[x].set(0)
                
    def refresh(self): # get any settings changes from PyMOL and update menus
        lst = cmd.get_setting_updates()
        for a in lst:
            if a in self.active_list:
                name = self.active_dict[a]
                if self.update_code.has_key(name):
                    code = self.update_code[name]
                    if code!=None:
                        tup = cmd.get_setting_tuple(name,'',0)
                        if tup!=None:
                            apply(code,(self,tup))
        
