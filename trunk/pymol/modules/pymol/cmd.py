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

# cmd.py 
# Python interface module for PyMol
#
# **This is the only module which should be/need be imported by 
# ** PyMol API Based Programs

# NOTE: this cmd module has grown absurdly huge, therefore expect this
# file to be broken down sometime after the 0.51 release into more
# manageable pieces.  Note that cmd.<whatever> will still work -- the
# symbols will be mapped into this namespace even after the code modules
# are separated.

# NEW CALL RETURN CONVENTIONS for _cmd.so C-layer
#
# (1) Calls into C (_cmd) should return results/status and print errors
#     and feedback (according to mask) BUT NEVER RAISE EXCEPTIONS.
## (2) In the absence of an expected return value, truth applies:
#        Success is 1, true 
#        Failure is 0, false. None, NULL
#
#     NOTE: if you need more informative failure, then code and
#           document an exception to this rule for your functions!
#
# (3) If _cmd produces a specific return result, be sure to include an
#     error result as one of the possibilities outside the range of the
#     expected return value.  For example, a negative distance or count
#
# (4) cmd.py API wrappers can then raise exceptions based on truth
#     and should return truth for success or None for failure
#     (if no exception was raised)
#

import re
import _cmd
import string
import thread
import threading
import types
import pymol
import os
import parsing
import sys

from shortcut import Shortcut

from chempy import io

#######################################################################
# symbols for early export
#######################################################################

file_ext_re= re.compile(string.join([
   "\.pdb$|\.ent$|\.mol$|",
   r"\.PDB$|\.ENT$|\.MOL$|",
   r"\.mmod$|\.mmd$|\.dat$|\.out$|",
   r"\.MMOD$|\.MMD$|\.DAT$|\.OUT$|",
   r"\.xplor$|\.pkl$|\.sdf$|", 
   r"\.XPLOR$|\.PKL$|\.SDF$|",                        
   r"\.r3d$|\.xyz$|\.xyz_[0-9]*$|", 
   r"\.R3D$|\.XYZ$|\.XYZ_[0-9]*$|",
   r"\.cc1$|\.cc2$|", # ChemDraw 3D
   r"\.CC1$|\.CC2$|",
   r"\.pmo$|", # Experimental molecular object format
   r"\.PMO$|",
   r"\.ccp4$|\.CCP4$" # CCP4   
   ],''))

safe_oname_re = re.compile(r"\+|\(|\)|\||\&|\!|\,")  # quash reserved characters
   
QuietException = parsing.QuietException

#--------------------------------------------------------------------
# shortcuts...

toggle_dict = {'on':1,'off':0,'1':1,'0':0}
toggle_sc = Shortcut(toggle_dict.keys())

stereo_dict = {'on':1,'off':0,'1':1,'0':0,'swap':-1,
               'crosseye':2,'quadbuffer':3} #,'walleye':3}
stereo_sc = Shortcut(stereo_dict.keys())

repres = {
   'everything'    : -1,
   'sticks'        : 0,
   'spheres'       : 1,
   'surface'       : 2,
   'labels'        : 3,
   'nb_spheres'    : 4,
   'cartoon'       : 5,
   'ribbon'        : 6,
   'lines'         : 7,
   'mesh'          : 8,
   'dots'          : 9,
   'dashes'        :10,
   'nonbonded'     :11,
   'cell'          :12,
   'cgo'           :13,
   'callback'      :14,
   'extent'        :15,   
}
repres_sc = Shortcut(repres.keys())

def null_function():
   pass

#--------------------------------------------------------------------
# convenient type-checking

def is_string(obj):
   return isinstance(obj,types.StringType)

def is_list(obj):
   return isinstance(obj,types.ListType)

#--------------------------------------------------------------------
# locks and threading

# the following lock is used by both C and Python to insure that no more than
# one active thread enters PyMOL at a given time. 

lock_api = pymol.lock_api
lock_api_c = pymol.lock_api_c

def lock_c(): # INTERNAL
   lock_api_c.acquire(1)

def unlock_c(): # INTERNAL
   lock_api_c.release()

def lock(): # INTERNAL
   lock_api.acquire(1) 

def lock_attempt(): # INTERNAL
   return lock_api.acquire(blocking=0)

def unlock(): # INTERNAL
   if (thread.get_ident() == pymol.glutThread):
      lock_api.release()
      _cmd.flush_now()
   else:
      lock_api.release()
      if _cmd.wait_queue(): # commands waiting to be executed?
         w = 0.0025 # NOTE: affects API perf. for "do" and delayed-exec
         while 1:
            e = threading.Event() # using this for portable delay
            e.wait(w)
            del e
            if not _cmd.wait_queue():
               break
            if w > 0.1: # wait up 0.2 sec max for PyMOL to flush queue
               if _feedback(fb_module.cmd,fb_mask.debugging):
                  fb_debug.write("Debug: avoiding possible dead-lock?\n")
               break
            w = w * 2 # wait twice as long each time until flushed

def is_glut_thread(): # internal
   if thread.get_ident() == pymol.glutThread:
      return 1
   else:
      return 0

#--------------------------------------------------------------------
# Feedback

class fb_action:
   set = 0
   enable = 1
   disable = 2
   push = 3
   pop = 4
   
fb_action_sc = Shortcut(fb_action.__dict__.keys())

class fb_module:

# This first set represents internal C systems

   all                       =0
   isomesh                   =1
   map                       =2
   matrix                    =3
   mypng                     =4
   triangle                  =5
   match                     =6
   raw                       =7
   isosurface                =8
   
   feedback                  =12
   scene                     =13
   threads                   =14  
   symmetry                  =15
   ray                       =16
   setting                   =17
   object                    =18
   ortho                     =19
   movie                     =20
   python                    =21
   extrude                   =22
   rep                       =23
   shaker                    =24
   
   coordset                  =25
   distset                   =26

   objectmolecule            =30
   objectmap                 =31
   objectmesh                =32
   objectdist                =33 
   objectcgo                 =34
   objectcallback            =35
   objectsurface             =36
   
   repwirebond               =45
   repcylbond                =46
   replabel                  =47
   repsphere                 =49
   repsurface                =50
   repmesh                   =51
   repdot                    =52
   repnonbonded              =53
   repnonbondedsphere        =54
   repdistdash               =55
   repdistlabel              =56
   repribbon                 =57
   repcartoon                =58
   sculpt                    =59
   
   executive                 =70
   selector                  =71
   editor                    =72

   export                    =75
   ccmd                      =76
   api                       =77   
   
   main                      =80  

# This second set, with negative indices
# represent python level systems

   parser                    =-1
   cmd                       =-2
   
fb_module_sc = Shortcut(fb_module.__dict__.keys())

class fb_mask:
   results =             0x01
   errors =              0x02
   actions =             0x04
   warnings =            0x08
   details =             0x10
   blather =             0x20
   debugging =           0x80
   everything =          0xFF
   
fb_mask_sc = Shortcut(fb_mask.__dict__.keys())

fb_dict ={}

for a in fb_module.__dict__.keys():
   vl = getattr(fb_module,a)
   if vl<0:
      fb_dict[vl] = 0x1F # default mask

fb_debug = sys.stderr # can redirect python debugging output elsewhere if desred...
   
def feedback(action="?",module="?",mask="?"):
   '''
DESCRIPTION

   "feedback" allows you to control what and how much text is output
   from PyMOL.

USAGE
   
   feedback action,module,mask

   action is one of ['set','enable','disable']
   module is a space-separated list of strings or simply "all"
   mask is a space-separated list of strings or simply "everything"

NOTES:

   "feedback" alone will print a list of the available choices
   
PYMOL API

   cmd.feedback(string action,string module,string mask)
   
EXAMPLES

   feedback enable, all , debugging
   feedback disable, selector, warnings actions
   feedback enable, main, blather
'''
   r = None

   # validate action

   if action=="?":
      print " feedback: available actions: set, enable, disable"
      act_kee = 0
   else:
      act_kee = fb_action_sc.interpret(action)
      if act_kee == None:
         print "Error: invalid feedback action '%s'."%action
         if _raising():
            raise QuietException
         else:
            return None
      elif not is_string(act_kee):
         print "Error: ambiguous feedback action '%s'."%action
         print action_amb
         if _raising():
            raise QuietException
         else:
            return None
      act_int = int(getattr(fb_action,act_kee))

   if (act_int<3) and ("?" in [action,module,mask]):
      if module=="?":
         print " feedback: available modules:"
         for a in fb_module.__dict__.keys():
            if a[0]!='_':
               print "   ",a
      if mask=="?":
         print " feedback: available masks:"
         for a in fb_mask.__dict__.keys():
            if a[0]!='_':
               print "   ",a
   else:
      if (act_int>=3):
         module='all'
         mask='everything'
            
      # validate and combine masks
      
      mask_int = 0
      mask_lst = string.split(mask)
      for mask in mask_lst:
         mask_kee = fb_mask_sc.interpret(mask)
         if mask_kee == None:
            print "Error: invalid feedback mask '%s'."%mask
            if _raising(): raise QuietException
            else: return None
         elif not is_string(mask_kee):
            print "Error: ambiguous feedback mask '%s'."%mask
            if _raising(): raise QuietException
            else: return None
         mask_int = int(getattr(fb_mask,mask_kee))

      # validate and iterate modules
      
      mod_lst = string.split(module)
      for module in mod_lst:
         mod_kee = fb_module_sc.interpret(module)
         if mod_kee == None:
            print "Error: invalid feedback module '%s'."%module
            if _raising(): raise QuietException
            else: return None
         elif not is_string(mod_kee):
            print "Error: ambiguous feedback module '%s'."%module
            if _raising(): raise QuietException
            else: return None
         mod_int = int(getattr(fb_module,mod_kee))
         if mod_int>=0:
            try:
               lock()
               r = _cmd.set_feedback(act_int,mod_int,mask_int)
            finally:
               unlock()
         if mod_int<=0:
            if mod_int:
               if act_int==0:
                  fb_dict[mod_int] = mask_int
               elif act_int==1:
                  fb_dict[mod_int] = fb_dict[mod_int] | mask_int
               elif act_int==2:
                  fb_dict[mod_int] = fb_dict[mod_int] & ( 0xFF - mask_int )
            else:
               for mod_int in fb_dict.keys():
                  if act_int==0:
                     fb_dict[mod_int] = mask_int
                  elif act_int==1:
                     fb_dict[mod_int] = fb_dict[mod_int] | mask_int
                  elif act_int==2:
                     fb_dict[mod_int] = fb_dict[mod_int] & ( 0xFF - mask_int )
            if _feedback(fb_module.feedback,fb_mask.debugging):
                sys.stderr.write(" feedback: mode %d on %d mask %d\n"%(
                   act_int,mod_int,mask_int))
   return r

#--------------------------------------------------------------------
# internal API routines

# status reporting

def _feedback(module,mask): # feedback test routine
   r = 0
   module = int(module)
   mask = int(mask)
   if module>0:
      try:
         lock()
         r = _cmd.feedback(module,mask)
      finally:
         unlock()
   else:
      if fb_dict.has_key(module):
         r = fb_dict[module]&mask
   return r

# movie rendering

def _mpng(*arg): # INTERNAL
   try:
      lock()   
      fname = arg[0]
      if re.search("\.png$",fname):
         fname = re.sub("\.png$","",fname)
      fname = os.path.expanduser(fname)
      fname = os.path.expandvars(fname)
      r = _cmd.mpng_(str(fname),arg[1],arg[2])
   finally:
      unlock()
   return r

# loading

def _load(oname,finfo,state,ftype,finish,discrete):
   # caller must already hold API lock
   # NOTE: state index assumes 1-based state
   r = 1
   if ftype not in (loadable.model,loadable.brick):
      if ftype == loadable.r3d:
         import cgo
         obj = cgo.from_r3d(finfo)
         if obj:
            _cmd.load_object(str(oname),obj,int(state)-1,loadable.cgo,
                             int(finish),int(discrete))
         else:
            print " load: couldn't load raster3d file."
      elif ftype == loadable.cc1: # ChemDraw 3D
         obj = io.cc1.fromFile(finfo)
         if obj:
            _cmd.load_object(str(oname),obj,int(state)-1,loadable.model,
                             int(finish),int(discrete))            
      else:
         r = _cmd.load(str(oname),finfo,int(state)-1,int(ftype),
                       int(finish),int(discrete))
   else:
      try:
         x = io.pkl.fromFile(finfo)
         if isinstance(x,types.ListType) or isinstance(x,types.TupleType):
            for a in x:
               r = _cmd.load_object(str(oname),a,int(state)-1,
                                    int(ftype),0,int(discrete))
               if(state>0):
                  state = state + 1
            _cmd.finish_object(str(oname))
         else:
            r = _cmd.load_object(str(oname),x,
                                 int(state)-1,int(ftype),
                                 int(finish),int(discrete))
      except:
         print 'Error: can not load file "%s"' % finfo
   return r

# function keys and other specials

def _special(k,x,y): # INTERNAL (invoked when special key is pressed)
   k=int(k)
   if special.has_key(k):
      if special[k][1]:
         apply(special[k][1],special[k][2],special[k][3])
   return None

# control keys

def _ctrl(k):
   if ctrl.has_key(k):
      ck = ctrl[k]
      if ck[0]!=None:
         apply(ck[0],ck[1],ck[2])
   return None

# alt keys

def _alt(k):
   if alt.has_key(k):
      ak=alt[k]
      if ak[0]!=None:
         apply(ak[0],ak[1],ak[2])
   return None

# writing PNG files (thread-unsafe)

def _png(a): # INTERNAL - can only be safely called by GLUT thread 
   try:
      lock()   
      fname = a
      if not re.search("\.png$",fname):
         fname = fname +".png"
      fname = os.path.expanduser(fname)
      fname = os.path.expandvars(fname)         
      r = _cmd.png(str(fname))
   finally:
      unlock()
   return r

# quitting (thread-specific)

def _quit():
   try:
      lock()
      r = _cmd.quit()
   finally:
      unlock()
   return r

# screen redraws (thread-specific)

def _refresh(swap_buffers=1):  # Only call with GLUT thread!
   try:
      lock()
      if thread.get_ident() == pymol.glutThread:
         if swap_buffers:
            r = _cmd.refresh_now()
         else:
            r = _cmd.refresh()
      else:
         print "Error: Ignoring an unsafe call to cmd._refresh"
   finally:      unlock()
   return r

# stereo (platform dependent )

def _sgi_stereo(flag): # SGI-SPECIFIC - bad bad bad
   
   if os.path.exists("/usr/gfx/setmon"):
      if flag:
         os.system("/usr/gfx/setmon -n 1024x768_96s")
      else:
         os.system("/usr/gfx/setmon -n 72hz")
      
# color alias interpretation

def _interpret_color(color):
   _validate_color_sc()
   new_color = color_sc.interpret(color)
   if new_color:
      if is_string(new_color):
         return new_color
      else:
         color_sc.auto_err(color,'color')
   else:
      return color

def _validate_color_sc():
   global color_sc
   if color_sc == None: # update color shortcuts if needed
      lst = get_color_indices()
      color_sc = Shortcut(map(lambda x:x[0],lst))
      color_dict = {}
      for a in lst: color_dict[a[0]]=a[1]

def _invalidate_color_sc():
   global color_sc
   color_sc = None

def _get_color_sc():
   _validate_color_sc()
   return color_sc

def _get_feedback(): # INTERNAL
   l = []
   if lock_attempt():
      try:
         r = _cmd.get_feedback()
         while r:
            l.append(r)
            r = _cmd.get_feedback()
      finally:
         unlock()
   return l
get_feedback = _get_feedback # for legacy compatibility

# testing tools

# for comparing floating point numbers calculated using
# different FPUs

def _dump_floats(lst,format="%7.3f",cnt=9):
   c = cnt
   for a in lst:
      print format%a,
      c = c -1
      if c<=0:
         print
         c=cnt
   if c!=cnt:
      print
   
# HUH?
def _adjust_coord(a,i,x):
   a.coord[i]=a.coord[i]+x
   return None

def _raising():
   return get_setting_legacy("raise_exceptions")

#######################################################################
# now import modules which depend on the above
#######################################################################

import editor

#######################################################################
# cmd module functions...
#######################################################################

def ready(): # INTERNAL
   return _cmd.ready()

def setup_global_locks(): # INTERNAL, OBSOLETE?
   pass
   
def null():
   pass

# for extending the language

def extend(name,function):
   '''
DESCRIPTION

   "extend" is an API-only function which binds a new external
   function as a command into the PyMOL scripting language.

PYMOL API

   cmd.extend(string name,function function)
   
PYTHON EXAMPLE

   def foo(moo=2): print moo
   cmd.extend('foo',foo)

   The following would now be valid within PyMOL:

   foo
   foo 3
   foo moo=5

SEE ALSO

   alias, api
   '''
   keyword[name] = [function, 0,0,',',parsing.STRICT]
   kwhash.append(name)
   

# for aliasing compound commands to a single keyword

def alias(name,command):
   '''
DESCRIPTION

   "alias" allows you to bind a commonly used command to a single word

USAGE
   
   alias name, command-sequence

PYMOL API

   cmd.alias(string name,string command)
   
EXAMPLES

   alias go,load "test.pdb"; zoom (i;500); show sticks,(i;500 a;4)
   go

SEE ALSO

   extend
   '''
   keyword[name] = [eval("lambda :do('''%s ''')"%command), 0,0,',',parsing.STRICT]
   kwhash.append(name)

def write_html_ref(file):
   lst = cmd.globals()
   f=open(file,'w')
   head = 'H2'
   f.write("<HTML><BODY><H1>Reference</H1>")
   kees = lst.keys()
   kees.sort()
   for a in kees:
      if hasattr(lst[a],'__doc__'):
         if a[0:1]!='_' and (a not in ['string','thread',
                                       'setup_global_locks',
                                       'real_system','sys','imp','glob'
                                       ]):
            doc = lst[a].__doc__
            if is_string(doc):
               if len(doc):
                  doc = string.strip(doc)
                  doc = string.replace(doc,"<","&lt;")
                  f.write("<HR SIZE=1><%s>"%head+a+"</%s>\n"%head)
                  f.write("<PRE>"+string.strip(doc)+"\n\n</PRE>")
   f.write("</BODY></HTML>")
   f.close()

def dummy(*arg):
   '''
DESCRIPTION

   This is a dummy function which simply returns None.
   Don't you wish your life was so easy?
   '''
   #'
   return None

#####################################################################
# Here is where the PyMOL Command Language and API are built.
#####################################################################


# first we need to import a set of symbols into this module's local
# namespace

#--------------------------------------------------------------------
from importing import \
     load,               \
     load_brick,         \
     load_callback,      \
     load_cgo,           \
     load_map,           \
     load_model,         \
     load_object,        \
     read_mmodstr,       \
     read_molstr,        \
     read_pdbstr,        \
     loadable    

#--------------------------------------------------------------------
from creating import \
     copy,               \
     create,             \
     fragment,           \
     isodot,             \
     isomesh,            \
     isosurface,         \
     symexp,             \
     map_new

#--------------------------------------------------------------------
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
     sync

#--------------------------------------------------------------------
import controlling
from controlling import \
     button,             \
     config_mouse,       \
     mouse,              \
     mask,               \
     set_key,            \
     unmask,             \
     edit_mode

#--------------------------------------------------------------------
from querying import \
     count_atoms,        \
     count_frames,       \
     count_states,       \
     dist,               \
     distance,           \
     export_dots,        \
     find_pairs,         \
     get_area,           \
     get_color_indices,  \
     get_color_tuple,    \
     get_dihedral,       \
     get_extent,         \
     get_model,          \
     get_names,          \
     get_names_of_type,  \
     get_phipsi,         \
     get_position,       \
     get_povray,         \
     get_renderer,       \
     get_title,          \
     get_type,           \
     id_atom,            \
     identify,           \
     index,              \
     overlap,            \
     phi_psi

#--------------------------------------------------------------------
from selecting import \
     deselect,           \
     indicate,           \
     select             

#--------------------------------------------------------------------
from exporting import \
     png,                \
     export_coords,      \
     multisave,          \
     save               

#--------------------------------------------------------------------
import editing
from editing import \
     alter,              \
     alter_state,        \
     attach,             \
     bond,               \
     cycle_valence,      \
     deprotect,          \
     edit,               \
     flag,               \
     fuse,               \
     h_add,              \
     h_fill,             \
     invert,             \
     iterate,            \
     iterate_state,      \
     map_set_border,     \
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
     set_geometry,       \
     set_title,          \
     smooth,             \
     sort,               \
     torsion,            \
     transform_object,   \
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
     refresh_wizard,     \
     set_wizard,         \
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
from moving import \
     mset,              \
     mclear,            \
     mdo,               \
     mappend,           \
     mmatrix,           \
     mpng,              \
     forward,           \
     backward,          \
     rewind,            \
     middle,            \
     ending,            \
     mplay,             \
     mstop,             \
     mpng,              \
     mray,              \
     frame,             \
     get_state,         \
     get_frame         

#--------------------------------------------------------------------
import viewing
from viewing import \
     bg_color,           \
     cartoon,            \
     clip,               \
     color,              \
     dirty,              \
     disable,            \
     enable,             \
     full_screen,        \
     get_view,           \
     hide,               \
     label,              \
     load_png,           \
     meter_reset,        \
     move,               \
     orient,             \
     origin,             \
     ray,                \
     rebuild,            \
     recolor,            \
     refresh,            \
     reset,              \
     rock,               \
     set_color,          \
     set_view,           \
     show,               \
     stereo,             \
     turn,               \
     view,               \
     viewport,           \
     zoom


#--------------------------------------------------------------------
import setting
from setting import \
     set,                 \
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
     focus,              \
     import_coords,      \
     load_coords,        \
     mem,                \
     minimize,           \
     spheroid,           \
     test

#--------------------------------------------------------------------
from m4x import \
     metaphorics

#--------------------------------------------------------------------
# Modules which contain programs used explicity as "module.xxx"

import util
import movie

# This is the main dictionary

keyword = {

   # keyword : [ command, # min_arg, max_arg, separator, mode ]

   # NOTE: min_arg, max_arg, and separator, are hold-overs from the
   #       original PyMOL parser which will eventually be removed.
   #       all new commands should use NO_CHECK or STRICT modes
   #       which make much better use of built-in python features.
   'abort'         : [ dummy             , 0 , 0 , ''  , parsing.ABORT ],
   'alias'         : [ alias             , 0 , 0 , ''  , parsing.LITERAL1 ],
   'align'         : [ align             , 0 , 0 , ''  , parsing.STRICT ],
   'alter'         : [ alter             , 0 , 0 , ''  , parsing.LITERAL1 ],
   'alter_state'   : [ alter_state       , 0 , 0 , ''  , parsing.LITERAL2 ],
   'attach'        : [ attach            , 0 , 0 , ''  , parsing.LITERAL2 ],   
   'backward'      : [ backward          , 0 , 0 , ''  , parsing.STRICT ],
   'bg_color'      : [ bg_color          , 0 , 0 , ''  , parsing.STRICT ],
   'bond'          : [ bond              , 0 , 0 , ''  , parsing.STRICT ],
   'button'        : [ button            , 0 , 0 , ''  , parsing.STRICT ],
   'cartoon'       : [ cartoon           , 0 , 0 , ''  , parsing.STRICT ],
   'cd'            : [ cd                , 0 , 0 , ''  , parsing.STRICT ],  
   'check'         : [ check             , 0 , 0 , ''  , parsing.STRICT ],
   'clip'          : [ clip              , 0 , 0 , ''  , parsing.STRICT ],
   'cls'           : [ cls               , 0 , 0 , ''  , parsing.STRICT ],
   'color'         : [ color             , 0 , 0 , ''  , parsing.STRICT ],
   'commands'      : [ helping.commands  , 0 , 0 , ''  , parsing.STRICT ],
   'config_mouse'  : [ config_mouse      , 0 , 0 , ''  , parsing.STRICT ],
   'copy'          : [ copy              , 0 , 0 , ''  , parsing.LEGACY ],
   'count_atoms'   : [ count_atoms       , 0 , 0 , ''  , parsing.STRICT ],
   'count_frames'  : [ count_frames      , 0 , 0 , ''  , parsing.STRICT ],   
   'count_states'  : [ count_states      , 0 , 0 , ''  , parsing.STRICT ],
   'cycle_valence' : [ cycle_valence     , 0 , 0 , ''  , parsing.STRICT ],
   'create'        : [ create            , 0 , 0 , ''  , parsing.LEGACY ],   
   'delete'        : [ delete            , 0 , 0 , ''  , parsing.STRICT ],
   'deprotect'     : [ deprotect         , 0 , 0 , ''  , parsing.STRICT ],
   'deselect'      : [ deselect          , 0 , 0 , ''  , parsing.STRICT ],
   'dir'           : [ ls                , 0 , 0 , ''  , parsing.STRICT ],  
   'disable'       : [ disable           , 0 , 0 , ''  , parsing.STRICT ],
   'distance'      : [ distance          , 0 , 0 , ''  , parsing.LEGACY ],   
   'dump'          : [ dump              , 0 , 0 , ''  , parsing.STRICT ],
   'dummy'         : [ dummy             , 0 , 0 , ''  , parsing.STRICT ],   
   'edit'          : [ edit              , 0 , 0 , ''  , parsing.STRICT ],
   'edit_mode'     : [ edit_mode         , 0 , 0 , ''  , parsing.STRICT ],
   'enable'        : [ enable            , 0 , 0 , ''  , parsing.STRICT ],
   'ending'        : [ ending            , 0 , 0 , ''  , parsing.STRICT ],
   'export_dots'   : [ export_dots       , 0 , 0 , ''  , parsing.STRICT ],
   'extend'        : [ extend            , 0 , 0 , ''  , parsing.STRICT ],
   'fast_minimize' : [ fast_minimize     , 1,  4 , ',' , parsing.SIMPLE ], # TO REMOVE
   'feedback'      : [ feedback          , 0,  0 , ''  , parsing.STRICT ],
   'fit'           : [ fit               , 0 , 0 , ''  , parsing.STRICT ],
   'flag'          : [ flag              , 0 , 0 , ''  , parsing.LEGACY ],
   'fork'          : [ helping.spawn    , 1 , 2 , ',' , parsing.SPAWN  ],
   'forward'       : [ forward           , 0 , 0 , ''  , parsing.STRICT ],
   'fragment'      : [ fragment          , 0 , 0 , ''  , parsing.STRICT ],
   'full_screen'   : [ full_screen       , 0 , 0 , ''  , parsing.STRICT ],
   'fuse'          : [ fuse              , 0 , 0 , ''  , parsing.STRICT ],
   'frame'         : [ frame             , 0 , 0 , ''  , parsing.STRICT ],
   'get_area'      : [ get_area          , 0 , 0 , ''  , parsing.STRICT ],
   'get_dihedral'  : [ get_dihedral      , 0 , 0 , ''  , parsing.STRICT ],
   'get_extent'    : [ get_extent        , 0 , 0 , ''  , parsing.STRICT ],   
   'get_position'  : [ get_position      , 0 , 0 , ''  , parsing.STRICT ],
   'get_title'     : [ get_title         , 0 , 0 , ''  , parsing.STRICT ],   
   'get_type'      : [ get_type          , 0 , 0 , ''  , parsing.STRICT ],
   'get_view'      : [ get_view          , 0 , 0 , ''  , parsing.STRICT ],
   'h_add'         : [ h_add             , 0 , 0 , ''  , parsing.STRICT ],
   'h_fill'        : [ h_fill            , 0 , 0 , ''  , parsing.STRICT ],
   'help'          : [ help              , 0 , 0 , ''  , parsing.STRICT ],
   'hide'          : [ hide              , 0 , 0 , ''  , parsing.STRICT ],
   'id_atom'       : [ id_atom           , 0 , 0 , ''  , parsing.STRICT ],
   'identify'      : [ identify          , 0 , 0 , ''  , parsing.STRICT ],
   'index'         : [ index             , 0 , 0 , ''  , parsing.STRICT ],
   'indicate'      : [ indicate          , 0 , 0 , ''  , parsing.STRICT ],   
   'intra_fit'     : [ intra_fit         , 0 , 0 , ''  , parsing.STRICT ],
   'intra_rms'     : [ intra_rms         , 0 , 0 , ''  , parsing.STRICT ],
   'intra_rms_cur' : [ intra_rms_cur     , 0 , 0 , ''  , parsing.STRICT ],
   'invert'        : [ invert            , 0 , 0 , ''  , parsing.STRICT ],
   'isodot'        : [ isodot            , 0 , 0 , ''  , parsing.LEGACY ],   
   'isomesh'       : [ isomesh           , 0 , 0 , ''  , parsing.LEGACY ],
   'isosurface'    : [ isosurface        , 0 , 0 , ''  , parsing.LEGACY ],   
   'iterate'       : [ iterate           , 0 , 0 , ''  , parsing.LITERAL1 ],
   'iterate_state' : [ iterate_state     , 0 , 0 , ''  , parsing.LITERAL2 ],
   'label'         : [ label             , 0 , 0 , ''  , parsing.LITERAL1 ],
   'load'          : [ load              , 0 , 0 , ''  , parsing.STRICT ],
   'load_png'      : [ load_png          , 0 , 0 , ''  , parsing.STRICT ],
   'log'           : [ log               , 0 , 0 , ''  , parsing.STRICT ],
   'log_close'     : [ log_close         , 0 , 0 , ''  , parsing.STRICT ],
   'log_open'      : [ log_open          , 0 , 0 , ''  , parsing.STRICT ],
   'ls'            : [ ls                , 0 , 0 , ''  , parsing.STRICT ],  
   'mask'          : [ mask              , 0 , 0 , ''  , parsing.STRICT ],
   'map_set_border': [ map_set_border    , 0 , 0 , ''  , parsing.STRICT ],
   'map_new'       : [ map_new           , 0 , 0 , ''  , parsing.STRICT ],    
   'mappend'       : [ mappend           , 2 , 2 , ':' , parsing.SINGLE ], 
   'mem'           : [ mem               , 0 , 0 , ''  , parsing.STRICT ],
   'meter_reset'   : [ meter_reset       , 0 , 0 , ''  , parsing.STRICT ],
   'move'          : [ move              , 0 , 0 , ''  , parsing.STRICT ],
   'mset'          : [ mset              , 0 , 0 , ''  , parsing.STRICT ],
   'mdo'           : [ mdo               , 2 , 2 , ':' , parsing.SINGLE ],
   'mpng'          : [ mpng              , 0 , 0 , ''  , parsing.STRICT ],
   'mplay'         : [ mplay             , 0 , 0 , ''  , parsing.STRICT ],
   'mray'          : [ mray              , 0 , 0 , ''  , parsing.STRICT ],
   'mstop'         : [ mstop             , 0 , 0 , ''  , parsing.STRICT ],
   'mclear'        : [ mclear            , 0 , 0 , ''  , parsing.STRICT ],
   'middle'        : [ middle            , 0 , 0 , ''  , parsing.STRICT ],
   'minimize'      : [ minimize          , 0 , 4 , ',' , parsing.SIMPLE ], # TO REMOVE
   'mmatrix'       : [ mmatrix           , 0 , 0 , ''  , parsing.STRICT ],
   'mouse'         : [ mouse             , 0 , 0 , ''  , parsing.STRICT ],
   'multisave'     : [ multisave         , 0 , 0 , ''  , parsing.STRICT ],   
   'origin'        : [ origin            , 0 , 0 , ''  , parsing.STRICT ],
   'orient'        : [ orient            , 0 , 0 , ''  , parsing.STRICT ],
   'overlap'       : [ overlap           , 0 , 0 , ''  , parsing.STRICT ],
   'pair_fit'      : [ pair_fit          , 2 ,98 , ',' , parsing.SIMPLE ],
   'phi_psi'       : [ phi_psi           , 0 , 0 , ''  , parsing.STRICT ],
   'protect'       : [ protect           , 0 , 0 , ''  , parsing.STRICT ],
   'push_undo'     : [ push_undo         , 0 , 0 , ''  , parsing.STRICT ],   
   'pwd'           : [ pwd               , 0 , 0 , ''  , parsing.STRICT ],
   'ray'           : [ ray               , 0 , 0 , ''  , parsing.STRICT ],
   'rebuild'       : [ rebuild           , 0 , 0 , ''  , parsing.STRICT ],
   'recolor'       : [ recolor           , 0 , 0 , ''  , parsing.STRICT ],   
   'redo'          : [ redo              , 0 , 0 , ''  , parsing.STRICT ],
   'refresh'       : [ refresh           , 0 , 0 , ''  , parsing.STRICT ],
   'remove'        : [ remove            , 0 , 0 , ''  , parsing.STRICT ],
   'remove_picked' : [ remove_picked     , 0 , 0 , ''  , parsing.STRICT ],
   'rename'        : [ rename            , 0 , 0 , ''  , parsing.STRICT ],
   'replace'       : [ replace           , 0 , 0 , ''  , parsing.STRICT ],
   'reset'         : [ reset             , 0 , 0 , ''  , parsing.STRICT ],
   'resume'        : [ resume            , 0 , 0 , ''  , parsing.STRICT ],
   'rewind'        : [ rewind            , 0 , 0 , ''  , parsing.STRICT ],
   'rock'          : [ rock              , 0 , 0 , ''  , parsing.STRICT ],
   'rotate'        : [ rotate            , 0 , 0 , ''  , parsing.STRICT ],   
   'run'           : [ helping.run       , 1 , 2 , ',' , parsing.RUN    ],
   'rms'           : [ rms               , 0 , 0 , ''  , parsing.STRICT ],
   'rms_cur'       : [ rms_cur           , 0 , 0 , ''  , parsing.STRICT ],
   'save'          : [ save              , 0 , 0 , ''  , parsing.STRICT ],
   'sculpt_purge'  : [ sculpt_purge      , 0 , 0 , ''  , parsing.STRICT ],   
   'sculpt_deactivate': [ sculpt_deactivate , 0 , 0 , ''  , parsing.STRICT ],
   'sculpt_activate': [ sculpt_activate  , 0 , 0 , ''  , parsing.STRICT ],
   'sculpt_iterate': [ sculpt_iterate    , 0 , 0 , ''  , parsing.STRICT ],
   'select'        : [ select            , 0 , 0 , ''  , parsing.LEGACY ],
   'set'           : [ set               , 0 , 0 , ''  , parsing.LEGACY ],
   'set_color'     : [ set_color         , 0 , 0 , ''  , parsing.LEGACY ],
   'set_dihedral'  : [ set_dihedral      , 0 , 0 , ''  , parsing.STRICT ],   
   'set_geometry'  : [ set_geometry      , 0 , 0 , ''  , parsing.STRICT ],      
   'set_title'     : [ set_title         , 0 , 0 , ''  , parsing.STRICT ],   
   'set_key'       : [ set_key           , 0 , 0 , ''  , parsing.STRICT ], # API only
   'set_view'      : [ set_view          , 0 , 0 , ''  , parsing.STRICT ],   
   'show'          : [ show              , 0 , 0 , ''  , parsing.STRICT ],
   'smooth'        : [ smooth            , 0 , 0 , ''  , parsing.STRICT ],
   'sort'          : [ sort              , 0 , 0 , ''  , parsing.STRICT ],
   'spawn'         : [ helping.spawn     , 1 , 2 , ',' , parsing.SPAWN  ],
   'spheroid'      : [ spheroid          , 0 , 0 , ''  , parsing.STRICT ],
   'splash'        : [ splash            , 0 , 0 , ''  , parsing.STRICT ],
   '_special'      : [ _special          , 0 , 0 , ''  , parsing.STRICT ],
   'stereo'        : [ stereo            , 0 , 0 , ''  , parsing.STRICT ],
   'symexp'        : [ symexp            , 0 , 0 , ''  , parsing.LEGACY ],
   'system'        : [ system            , 0 , 0 , ''  , parsing.LITERAL ],
   'test'          : [ test              , 0 , 0 , ''  , parsing.STRICT ],
   'torsion'       : [ torsion           , 0 , 0 , ''  , parsing.STRICT ],
   'translate'     : [ translate         , 0 , 0 , ''  , parsing.STRICT ],
   'turn'          : [ turn              , 0 , 0 , ''  , parsing.STRICT ],
   'quit'          : [ quit              , 0 , 0 , ''  , parsing.STRICT ],
   '_quit'         : [ _quit             , 0 , 0 , ''  , parsing.STRICT ],
   'png'           : [ png               , 0 , 0 , ''  , parsing.STRICT ],
   'unbond'        : [ unbond            , 0 , 0 , ''  , parsing.STRICT ],
   'unpick'        : [ unpick            , 0 , 0 , ''  , parsing.STRICT ],
   'undo'          : [ undo              , 0 , 0 , ''  , parsing.STRICT ],
   'unmask'        : [ unmask            , 0 , 0 , ''  , parsing.STRICT ],
   'unprotect'     : [ deprotect         , 0 , 0 , ''  , parsing.STRICT ],
   'update'        : [ update            , 0 , 0 , ''  , parsing.STRICT ],
   'view'          : [ view              , 0 , 0 , ''  , parsing.STRICT ],   
   'viewport'      : [ viewport          , 0 , 0 , ''  , parsing.STRICT ],
   'wizard'        : [ wizard            , 0 , 0 , ''  , parsing.STRICT ],
   'zoom'          : [ zoom              , 0 , 0 , ''  , parsing.STRICT ],
# utility programs
   'util.cbag'     : [ util.cbag         , 0 , 0 , ''  , parsing.STRICT ],
   'util.cbac'     : [ util.cbac         , 0 , 0 , ''  , parsing.STRICT ],
   'util.cbay'     : [ util.cbay         , 0 , 0 , ''  , parsing.STRICT ],
   'util.cbas'     : [ util.cbas         , 0 , 0 , ''  , parsing.STRICT ],
   'util.cbap'     : [ util.cbap         , 0 , 0 , ''  , parsing.STRICT ],
   'util.cbak'     : [ util.cbak         , 0 , 0 , ''  , parsing.STRICT ],
   'util.cbaw'     : [ util.cbaw         , 0 , 0 , ''  , parsing.STRICT ],
   'util.cbab'     : [ util.cbab         , 0 , 0 , ''  , parsing.STRICT ],
   'util.cbc'      : [ util.cbc          , 0 , 0 , ''  , parsing.STRICT ],
   'util.mrock'    : [ util.mrock        , 0 , 0 , ''  , parsing.STRICT ], # LEGACY
   'util.mroll'    : [ util.mroll        , 0 , 0 , ''  , parsing.STRICT ], # LEGACY
   'util.ss'       : [ util.ss           , 0 , 0 , ''  , parsing.STRICT ],# secondary structure
   'util.rainbow'  : [ util.rainbow      , 0 , 0 , ''  , parsing.STRICT ],# secondary structure
# movie programs
   'movie.rock'    : [ movie.rock        , 0 , 0 , ''  , parsing.STRICT ],
   'movie.roll'    : [ movie.roll        , 0 , 0 , ''  , parsing.STRICT ],
   'movie.load'    : [ movie.load        , 0 , 0 , ''  , parsing.STRICT ],
   'movie.zoom'    : [ movie.zoom        , 0 , 0 , ''  , parsing.STRICT ],
   'movie.screw'   : [ movie.screw       , 0 , 0 , ''  , parsing.STRICT ],
   'movie.nutate'  : [ movie.nutate      , 0 , 0 , ''  , parsing.STRICT ],
   'movie.tdroll'  : [ movie.tdroll      , 0 , 0 , ''  , parsing.STRICT ],
# activate metaphorics extensions
   'metaphorics'   : [ metaphorics       , 0 , 0 , ''  , parsing.STRICT ],
   }

kwhash = Shortcut(keyword.keys())

# informational, or API-only functions which don't exist in the
# PyMOL command language namespace

help_only = {  
   'api'           : [ helping.api          , 0 , 0 , '' , 0 ],
   'selections'    : [ helping.selections   , 0 , 0 , '' , 0 ],
   'keyboard'      : [ helping.keyboard     , 0 , 0 , '' , 0 ],
   'mouse'         : [ helping.mouse        , 0 , 0 , '' , 0 ],
   'examples'      : [ helping.examples     , 0 , 0 , '' , 0 ],
   'read_molstr'   : [ read_molstr           , 0 , 0 , '' , 0 ],
   'release'       : [ helping.release      , 0 , 0 , '' , 0 ],   
   'launching'     : [ helping.launching    , 0 , 0 , '' , 0 ],
   'load_model'    : [ load_model            , 0 , 0 , '' , 0 ],
   'movies'        : [ helping.movies       , 0 , 0 , '' , 0 ],
   'editing'       : [ helping.editing      , 0 , 0 , '' , 0 ],  
   'edit_keys'     : [ helping.edit_keys    , 0 , 0 , '' , 0 ],
   'povray'        : [ helping.povray       , 0 , 0 , '' , 0 ],
   'get_names'     : [ get_names             , 0 , 0 , '' , 0 ],
   'get_type'      : [ get_type              , 0 , 0 , '' , 0 ],
   'faster'        : [ helping.faster       , 0 , 0 , '' , 0 ],
   'sync'          : [ sync                 , 0 , 0 , '' , 0],
   'transparency'  : [ helping.transparency , 0 , 0 , '' , 0 ],  
   '@'             : [ helping.at_sign      , 0 , 0 , '' , 0 ],  
}

help_sc = Shortcut(keyword.keys()+help_only.keys())


special = {
   1        : [ 'F1'        , None                   , () , {} ],
   2        : [ 'F2'        , None                   , () , {} ],
   3        : [ 'F3'        , None                   , () , {} ],
   4        : [ 'F4'        , None                   , () , {} ],
   5        : [ 'F5'        , None                   , () , {} ],
   6        : [ 'F6'        , None                   , () , {} ],
   7        : [ 'F7'        , None                   , () , {} ],
   8        : [ 'F8'        , None                   , () , {} ],
   9        : [ 'F9'        , None                   , () , {} ],
   10       : [ 'F10'       , None                   , () , {} ],
   11       : [ 'F11'       , redo                   , () , {} ],
   12       : [ 'F12'       , undo                   , () , {} ],
   100      : [ 'left'      , backward               , () , {} ],
   101      : [ 'up'        , None                   , () , {} ],
   102      : [ 'right'     , forward                , () , {} ],
   103      : [ 'down'      , None                   , () , {} ],
   104      : [ 'pgup'      , rewind                 , () , {} ],
   105      : [ 'pgdn'      , ending                 , () , {} ],
   106      : [ 'home'      , zoom                   , () , {} ],
   107      : [ 'end'       , ending                 , () , {} ],
   108      : [ 'insert'    , rock                   , () , {} ]   
}

ctrl = {
   'A' : [ redo                   , () , {}],
   'B' : [ replace                , ('Br',1,1), {} ],
   'C' : [ replace                , ('C',4,4), {} ],
   'D' : [ remove_picked          , () , {} ],
   'E' : [ invert                 , () , {} ],   
   'F' : [ replace                , ('F',1,1), {} ],   
   'G' : [ replace                , ('H',1,1), {} ],
   'I' : [ replace                , ('I',1,1), {} ],
   'J' : [ alter                  , ('pk1','formal_charge=-1.0'), {} ],
   'K' : [ alter                  , ('pk1','formal_charge =1.0'), {} ],
   'L' : [ replace                , ('Cl',1,1) , {}],
   'N' : [ replace                , ('N',4,3) , {}],
   'O' : [ replace                , ('O',4,2) , {}],   
   'P' : [ replace                , ('P',4,1) , {}],
   'Q' : [ dist                   , () , {}],   
   'R' : [ h_fill                 , () , {} ],   
   'S' : [ replace                , ('S',4,2) , {}],
   'T' : [ bond                   , () , {} ],   
   'U' : [ alter                  , ('pk1','formal_charge =0.0') , {}],
   'W' : [ cycle_valence          , () , {}],   
   'X' : [ remove                 , ('pkfrag1',) , {} ],
   'Y' : [ attach                 , ('H',1,1) , {} ],
   'Z' : [ undo                   , () , {} ],   
   }

alt = {
   '1' : [ editor.attach_fragment  , ("pk1","formamide",5,0), {}],
   '2' : [ editor.attach_fragment  , ("pk1","formamide",4,0), {}],
   '3' : [ editor.attach_fragment  , ("pk1","sulfone",3,1), {}],
   '4' : [ editor.attach_fragment  , ("pk1","cyclobutane",4,0), {}],
   '5' : [ editor.attach_fragment  , ("pk1","cyclopentane",5,0), {}],
   '6' : [ editor.attach_fragment  , ("pk1","cyclohexane",7,0), {}],
   '7' : [ editor.attach_fragment  , ("pk1","cycloheptane",8,0), {}],
   '8' : [ editor.attach_fragment  , ("pk1","cyclopentadiene",5,0), {}],
   '9' : [ editor.attach_fragment  , ("pk1","benzene",6,0), {}],
   '0' : [ editor.attach_fragment  , ("pk1","formaldehyde",2,0), {}],
   'a' : [ editor.attach_amino_acid, ("pk1","ala"), {}],
   'b' : [ editor.attach_amino_acid, ("pk1","ace"), {}],                                 
   'c' : [ editor.attach_amino_acid, ("pk1","cys"), {}],
   'd' : [ editor.attach_amino_acid, ("pk1","asp"), {}],
   'e' : [ editor.attach_amino_acid, ("pk1","glu"), {}],
   'f' : [ editor.attach_amino_acid, ("pk1","phe"), {}],
   
   'g' : [ editor.attach_amino_acid, ("pk1","gly"), {}],
   'h' : [ editor.attach_amino_acid, ("pk1","his"), {}],
   'i' : [ editor.attach_amino_acid, ("pk1","ile"), {}],

   'j' : [ editor.attach_fragment,   ("pk1","acetylene",2,0), {}],
   'k' : [ editor.attach_amino_acid, ("pk1","lys"), {}],
   'l' : [ editor.attach_amino_acid, ("pk1","leu"), {}],
   
   'm' : [ editor.attach_amino_acid, ("pk1","met"), {}],
   'n' : [ editor.attach_amino_acid, ("pk1","asn"), {}],
   'p' : [ editor.attach_amino_acid, ("pk1","pro"), {}],
   'q' : [ editor.attach_amino_acid, ("pk1","gln"), {}],
   'r' : [ editor.attach_amino_acid, ("pk1","arg"), {}],
   
   's' : [ editor.attach_amino_acid, ("pk1","ser"), {}],
   't' : [ editor.attach_amino_acid, ("pk1","thr"), {}],
   'v' : [ editor.attach_amino_acid, ("pk1","val"), {}],
   'w' : [ editor.attach_amino_acid, ("pk1","trp"), {}],
   'y' : [ editor.attach_amino_acid, ("pk1","tyr"), {}],
   'z' : [ editor.attach_amino_acid, ("pk1","nme"), {}],                              
   }

selection_sc = lambda sc=Shortcut,gn=get_names:sc(gn('public')+['all'])
object_sc = lambda sc=Shortcut,gn=get_names:sc(gn('objects'))
map_sc = lambda sc=Shortcut,gnot=get_names_of_type:sc(gnot('object:map'))

# Table for argument autocompletion

auto_arg =[
   {
   'color'          : [ _get_color_sc          , 'color'           , ', ' ],
   'cartoon'        : [ viewing.cartoon_sc     , 'cartoon'         , ', ' ],      
   'set'            : [ setting.setting_sc     , 'setting'         , '='  ],
   'flag'           : [ editing.flag_sc        , 'flag'            , ', ' ],
   'show'           : [ repres_sc              , 'representation'  , ', ' ],
   'hide'           : [ repres_sc              , 'representation'  , ', ' ],
   'stereo'         : [ stereo_sc              , 'option'          , ''   ],
   'full_screen'    : [ toggle_sc              , 'option'          , ''   ],
   'clip'           : [ viewing.clip_action_sc , 'clipping action' , ', ' ],
   'feedback'       : [ fb_action_sc           , 'action'          , ', ' ],
   'button'         : [ controlling.button_sc  , 'button'          , ', ' ],
   'align'          : [ selection_sc           , 'selection'       , ','  ],
   'zoom'           : [ selection_sc           , 'selection'       , ''   ],
   'origin'         : [ selection_sc           , 'selection'       , ''   ],
   'protect'        : [ selection_sc           , 'selection'       , ''   ],
   'deprotect'      : [ selection_sc           , 'selection'       , ''   ],   
   'mask'           : [ selection_sc           , 'selection'       , ''   ],
   'unmask'         : [ selection_sc           , 'selection'       , ''   ],
   'delete'         : [ selection_sc           , 'selection'       , ''   ],
   'alter'          : [ selection_sc           , 'selection'       , ''   ],
   'iterate'        : [ selection_sc           , 'selection'       , ''   ],
   'iterate_state'  : [ selection_sc           , 'selection'       , ''   ],
   'help'           : [ help_sc                , 'selection'       , ''   ],         
   },
   {
   'align'          : [ selection_sc           , 'selection'       , ''   ],
   'feedback'       : [ fb_module_sc           , 'module'          , ', ' ],
   'button'         : [ controlling.but_mod_sc , 'modifier'        , ', ' ],
   'show'           : [ selection_sc           , 'selection'       , ''   ],
   'hide'           : [ selection_sc           , 'selection'       , ''   ],
   'color'          : [ selection_sc           , 'selection'       , ''   ],
   'select'         : [ selection_sc           , 'selection'       , ''   ],
   'save'           : [ selection_sc           , 'selection'       , ', ' ],
   'flag'           : [ selection_sc           , 'selection'       , ', ' ],   
   'load'           : [ selection_sc           , 'selection'       , ', ' ],
   'create'         : [ selection_sc           , 'selection'       , ', ' ],
   'symexp'         : [ object_sc              , 'object'          , ', ' ],   
   'isomesh'        : [ map_sc                 , 'map object'      , ', ' ],
   'view'           : [ viewing.view_sc        , 'view action'     , ''   ],            
   },
   {
   'feedback'       : [ fb_mask_sc             , 'mask'            , ''   ],
   'button'         : [ controlling.but_act_sc , 'button action'   , ''   ],
   'flag'           : [ editing.flag_action_sc , 'flag action'     , ''   ],      
   }
   ]
   
color_sc = None
