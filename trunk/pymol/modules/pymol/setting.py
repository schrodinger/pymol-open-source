
import types
from shortcut import Shortcut

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
   auto_show_lines       =51
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
   auto_show_nonbonded   =72
   cache_display        =73
   mesh_radius          =74
   backface_cull        =75
   gamma                =76
   dot_width            =77
   auto_show_selections  =78
   auto_hide_selections  =79
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
   cartoon_round_helices =111
   cartoon_refine_normals = 112
   cartoon_flat_sheets  =113
   cartoon_smooth_loops =114
   cartoon_dumbbell_length   =  115
   cartoon_dumbbell_width    =  116
   cartoon_dumbbell_radius   =  117
   cartoon_fancy_helices    =  118
   cartoon_fancy_sheets     =  119
   ignore_pdb_segi       =120
   ribbon_throw          =121
   cartoon_throw         =122
   cartoon_refine        =123
   cartoon_refine_tips   =124
   cartoon_discrete_colors  =125
   normalize_ccp4_maps   =126
   surface_poor          =127
   internal_feedback     =128
   cgo_line_width        =129
   cgo_line_radius       =130
   logging               =131
   
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
   }
   
boolean_sc = Shortcut(name_list)
  
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
