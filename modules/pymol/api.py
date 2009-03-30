

#--------------------------------------------------------------------
from importing import \
      finish_object,      \
      load,               \
      load_brick,         \
      load_callback,      \
      load_cgo,           \
      load_embedded,      \
      load_map,           \
      load_model,         \
      load_object,        \
      load_traj,          \
      load_raw,           \
      loadable,           \
      read_mmodstr,       \
      read_molstr,        \
      read_sdfstr,        \
      read_pdbstr,        \
      read_xplorstr,      \
      fetch,              \
      set_session,        \
      space              

#--------------------------------------------------------------------
import creating
from creating import \
      copy,               \
      create,             \
      extract,            \
      fragment,           \
      group,              \
      gradient,           \
      isodot,             \
      isolevel,           \
      isomesh,            \
      isosurface,         \
      map_new,            \
      pseudoatom,         \
      slice_new,          \
      symexp,             \
      ramp_new,           \
      ungroup

#--------------------------------------------------------------------
import commanding
from commanding import \
      cls,                \
      delete,             \
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

#--------------------------------------------------------------------
from querying import \
      angle,              \
      auto_measure,       \
      count_atoms,        \
      count_frames,       \
      count_states,       \
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
      get_idtf,           \
      get_modal_draw,     \
      get_model,          \
      get_movie_locked,   \
      get_movie_length,   \
      get_names,          \
      get_names_of_type,  \
      get_legal_name,     \
      get_unused_name,    \
      get_object_matrix,  \
      get_mtl_obj,        \
      get_phipsi,         \
      get_position,       \
      get_povray,         \
      get_raw_alignment,  \
      get_renderer,       \
      get_symmetry,       \
      get_title,          \
      get_type,           \
      get_version,        \
      get_vrml,           \
      id_atom,            \
      identify,           \
      index,              \
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
import exporting
from exporting import \
      copy_image,         \
      cache,              \
      export_coords,      \
      get_pdbstr,         \
      get_session,        \
      get_fastastr,       \
      multisave,          \
      png,                \
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
      drag,               \
      dss,                \
      edit,               \
      fix_chemistry,      \
      flag,               \
      fuse,               \
      get_editor_scheme,  \
      h_add,              \
      h_fill,             \
      h_fix,              \
      invert,             \
      iterate,            \
      iterate_state,      \
      map_set,            \
      map_set_border,     \
      map_double,         \
      map_halve,          \
      map_trim,           \
      matrix_copy,        \
      matrix_reset,       \
      protect,            \
      push_undo,          \
      reference,          \
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
      set_object_color,   \
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
      update,             \
      valence,            \
      vdw_fit 

from editor import \
      fab
      
from computing import \
      clean              

matrix_transfer = matrix_copy # legacy

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
      super,             \
      rms,               \
      rms_cur,           \
      intra_fit,         \
      intra_rms,         \
      intra_rms_cur,     \
      pair_fit          

#--------------------------------------------------------------------
# ARE ALL OF THESE UNUSED AND/OR DEPRECATED (?)
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
      show_as,            \
      bg_color,           \
      bg_colour,          \
      cartoon,            \
      capture,            \
      clip,               \
      color,              \
      colour,             \
      del_colorection,    \
      dirty,              \
      disable,            \
      draw,               \
      enable,             \
      full_screen,        \
      get_colorection,    \
      get_view,           \
      get_vis,            \
      get_scene_dict,     \
      get_scene_list,     \
      hide,               \
      label,              \
      label2,             \
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
      scene_order,        \
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
      set_bond,            \
      get,                 \
      unset,               \
      unset_bond,          \
      get_setting_boolean, \
      get_setting_int,     \
      get_setting_float,   \
      get_setting_legacy,  \
      get_setting_tuple,   \
      get_setting_updates, \
      get_setting_text

#--------------------------------------------------------------------
import helping
from helping import \
      abort,               \
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

#--------------------------------------------------------------------
#from m4x import \
#     metaphorics

#--------------------------------------------------------------------
# Modules which contain programs used explicity as "module.xxx"

import util
import movie

# dang! Python 2.6 will break PyMOL's "as" method. 
# Proposal:
#  1. stick with Python <=2.5 for as long as possible
#  2. convert API method to cmd.show_as() and leave "as" in the scripting langauge
#  3. allow "show_as" in the scripting language
globals()['as'] = show_as
