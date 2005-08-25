# -c

print "BEGIN-LOG"

import pymol
import cmd

# need to do this for all supported file types because the load routines diverge...

def x(st):
   print st
   cmd.set("raise_exceptions",0)
   valu = eval(st)
   print valu
   if cmd.is_error(valu):
      cmd.set("raise_exceptions",1)
      try:
         print eval(st)
      except pymol.CmdException:
         print "CmdException raised."
      except cmd.QuietException:
         print "QuietException raised."

map(x,[
'cmd.load("dat/pept.pdb")' ,
'cmd.load("dat/nonexistent.pdb")' ,

'cmd.load("dat/small01.mol")' ,
'cmd.load("dat/nonexistent.mol")' ,
'cmd.load("dat/ligs3d.sdf")' ,
'cmd.load("dat/nonexistent.sdf")',
   
'cmd.load("dat/small03.mol2")' ,
'cmd.load("dat/nonexistent.mol2")' ,

'cmd.load("dat/pept.r3d")' ,
'cmd.load("dat/nonexistent.r3d")' ,

'cmd.load("dat/pept.pkl")' ,
'cmd.load("dat/nonexistent.pkl")' ,

'cmd.get_names()' ,
'cmd.delete("all")' ,
'cmd.get_names()' ,

'cmd.get_model("nonexistent")' ,
   
'cmd.space("cmyk")',
'cmd.space("unknown")',
'cmd.space("rgb")',
'cmd.space()',

])

cmd.delete("all")
cmd.load("dat/pept.pdb")
mdl = cmd.get_model("pept")
print mdl.__class__
print len(mdl.atom)

mdl = cmd.get_model("none")
print mdl.__class__
print len(mdl.atom)

map( x, [
'cmd.get_model("nonexistent")',
'cmd.create("test","none",quiet=0)',
'cmd.create("test2","nonexistent")',
'cmd.create("test3","?allowed",quiet=0)',
'cmd.fragment("arg")',
'cmd.fragment("nonexistent")',
])

cmd.delete("all")
cmd.load("dat/pept.pdb")

map( x, [
'cmd.distance("none","none","none")',
'cmd.distance("none","invalid1","none")',
'cmd.distance("none","none","invalid2")',
'cmd.distance("none","invalid1","invalid2")',
'cmd.index("none")',
'cmd.index("first pept")',
'cmd.index("invalid")',
'cmd.identify("none")',
'cmd.identify("first pept")',
'cmd.identify("invalid")',
'cmd.identify("none",mode=1)',
'cmd.identify("first pept",mode=1)',
'cmd.identify("invalid",mode=1)',
'cmd.count_atoms("none")',
'cmd.count_atoms()',
'cmd.count_atoms("pept")',
'cmd.count_atoms("invalid")',
'cmd.count_frames()',
'cmd.count_states()',
'cmd.count_states("all")',
'cmd.count_states("none")',
'cmd.count_states("invalid")',
'cmd.get_area()',
'cmd.get_area("resi 1")',
'cmd.get_area("none")',
'cmd.get_area("invalid")',
'cmd.get_chains()',
'cmd.get_chains("pept")',
'cmd.get_chains("resi 1")',
'cmd.get_chains("none")',
'cmd.get_chains("invalid")',
'cmd.zoom()',
'cmd.zoom("pept")',
'cmd.zoom("none")',
'cmd.zoom("invalid")',
'cmd.center()',
'cmd.center("pept")',
'cmd.center("none")',
'cmd.center("invalid")',

])



print "END-LOG"

TODO="""

      #--------------------------------------------------------------------
      from querying import \
           angle,              \
           dist,               \
           dihedral,           \
           distance,           \
           export_dots,        \
           find_pairs,         \
           get_angle,          \
           get_area,           \
           get_chains,         \
           get_color_index,    \
           get_color_indices,  \
           get_object_color_index, \
           get_object_list,    \
           get_color_tuple,    \
           get_atom_coords,    \
           get_dihedral,       \
           get_distance,       \
           get_extent,         \
           get_movie_locked,   \
           get_names_of_type,  \
           get_object_matrix, \
           get_phipsi,         \
           get_position,       \
           get_povray,         \
           get_renderer,       \
           get_symmetry,       \
           get_title,          \
           get_type,           \
           id_atom,            \
           overlap,            \
           phi_psi

      #--------------------------------------------------------------------
      from selecting import \
           deselect,           \
           indicate,           \
           select,             \
           select_list,        \
           pop

      #--------------------------------------------------------------------
      from exporting import \
           png,                \
           export_coords,      \
           get_pdbstr,         \
           get_session,        \
           multisave,          \
           save               

      #--------------------------------------------------------------------
      import editing
      from editing import \
           alter,              \
           alter_list,         \
           alter_state,        \
           attach,             \
           bond,               \
           cycle_valence,      \
           deprotect,          \
           dss,                \
           edit,               \
           fix_chemistry,      \
           flag,               \
           fuse,               \
           h_add,              \
           h_fill,             \
           invert,             \
           iterate,            \
           iterate_state,      \
           map_set_border,     \
           map_double,         \
           matrix_transfer,    \
           matrix_reset,       \
           protect,            \
           push_undo,          \
           redo,               \
           remove,             \
           remove_picked,      \
           rename,             \
           replace,            \
           rotate,             \
           sculpt_purge,       \
           sculpt_deactivate,  \
           sculpt_activate,    \
           sculpt_iterate,     \
           set_dihedral,       \
           set_name,           \
           set_geometry,       \
           set_object_ttt,     \
           set_symmetry,       \
           set_title,          \
           smooth,             \
           sort,               \
           split_states,       \
           torsion,            \
           transform_object,   \
           transform_selection,\
           translate,          \
           translate_atom,     \
           unbond,             \
           undo,               \
           unpick,             \
           update

      #--------------------------------------------------------------------

      from externing import \
           cd,                 \
           ls,                 \
           paste,              \
           pwd,                \
           system

      #--------------------------------------------------------------------
      from wizarding import \
           get_wizard,         \
           get_wizard_stack,   \
           refresh_wizard,     \
           replace_wizard,     \
           set_wizard,         \
           set_wizard_stack,   \
           dirty_wizard,       \
           wizard

      #--------------------------------------------------------------------
      from fitting import \
           align,             \
           fit,               \
           rms,               \
           rms_cur,           \
           intra_fit,         \
           intra_rms,         \
           intra_rms_cur,     \
           pair_fit          

      #--------------------------------------------------------------------
      from preset import \
           simple,            \
           technical,         \
           pretty,         \
           publication

      #--------------------------------------------------------------------
      import moving
      from moving import \
           madd,              \
           mset,              \
           mclear,            \
           mdo,               \
           mappend,           \
           mmatrix,           \
           mdump,             \
           accept,            \
           decline,           \
           mpng,              \
           mview,             \
           forward,           \
           backward,          \
           rewind,            \
           middle,            \
           ending,            \
           mplay,             \
           mtoggle,           \
           mstop,             \
           mpng,              \
           mray,              \
           frame,             \
           get_movie_playing, \
           get_state,         \
           get_frame         

      #--------------------------------------------------------------------
      import viewing
      from viewing import \
           as,                 \
           bg_color,           \
           bg_colour,          \
           cartoon,            \
           clip,               \
           color,              \
           colour,             \
           del_colorection,    \
           dirty,              \
           disable,            \
           enable,             \
           full_screen,        \
           get_colorection,    \
           get_view,           \
           get_vis,            \
           get_scene_dict,     \
           hide,               \
           label,              \
           load_png,           \
           meter_reset,        \
           move,               \
           orient,             \
           origin,             \
           center,             \
           ray,                \
           rebuild,            \
           recolor,            \
           recolour,           \
           refresh,            \
           reset,              \
           rock,               \
           scene,              \
           set_color,          \
           set_colour,         \
           set_colorection,    \
           set_colorection_name,\
           set_vis,            \
           set_view,           \
           show,               \
           spectrum,           \
           stereo,             \
           toggle,             \
           turn,               \
           view,               \
           viewport,           \
           window,             \
           zoom
      #--------------------------------------------------------------------
      import setting
      from setting import \
           set,                 \
           get,                 \
           unset,               \
           get_setting_legacy,  \
           get_setting_tuple,   \
           get_setting_updates, \
           get_setting_text

      #--------------------------------------------------------------------
      import helping
      from helping import \
           show_help,           \
           help,                \
           commands

      #--------------------------------------------------------------------
      from experimenting import \
           check,              \
           dump,               \
           expfit,             \
           get_bond_print,     \
           fast_minimize,      \
           import_coords,      \
           load_coords,        \
           mem,                \
           minimize,           \
           spheroid,           \
           test
           
      from importing import \
           finish_object,      \
           load_brick,         \
           load_callback,      \
           load_cgo,           \
           load_embedded,      \
           load_map,           \
           load_model,         \
           load_object,        \
           load_traj,          \
           loadable,           \
           read_mmodstr,       \
           read_molstr,        \
           read_pdbstr,        \
           read_xplorstr,      \
           set_session,        \


      #--------------------------------------------------------------------
      import creating
      from creating import \
           copy,               \
           isodot,             \
           isolevel,           \
           isomesh,            \
           isosurface,         \
           slice_new,          \
           symexp,             \
           map_new,            \
           ramp_new

      #--------------------------------------------------------------------
      from commanding import \
           cls,                \
           do,                 \
           log,                \
           log_close,          \
           log_open,           \
           quit,               \
           resume,             \
           splash,             \
           reinitialize,       \
           sync

      #--------------------------------------------------------------------
      import controlling
      from controlling import \
           button,             \
           config_mouse,       \
           mouse,              \
           mask,               \
           order,              \
           set_key,            \
           unmask,             \
           edit_mode

"""



