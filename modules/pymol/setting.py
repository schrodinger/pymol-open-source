
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
   dot_size             =29
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
   autoshow_lines       =51
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
   autoshow_nonbonded   =72
   cache_display        =73
   mesh_radius          =74
   backface_cull        =75
   gamma                =76
   dot_width            =77
   autoshow_selections  =78
   autohide_selections  =79
   selection_width      =80
   selection_overlay    =81
   static_singletons    =82

setting_sc = Shortcut(SettingIndex.__dict__.keys())

name_dict = {}
lst = SettingIndex.__dict__.keys()
lst = filter(lambda x:x[0]!='_',lst)
lst = map(lambda x:(getattr(SettingIndex,x),x),lst)
for a in lst:
   name_dict[a[0]]=a[1]
del lst

boolean_dict = {
   "true" : 1,
   "false" : 0,
   "on"   : 1,
   "off"  : 0,
   "1"    : 1,
   "0"    : 0,
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

def _get_name(index):
   if name_dict.has_key(index):
      return name_dict[index]
   else:
      return ""

