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

import thread
import string
import types

import pymol
import selector

import cmd
from cmd import _cmd,lock,unlock,Shortcut,QuietException,_raising, \
     _feedback,fb_module,fb_mask, \
     repres,repres_sc, \
     toggle_dict,toggle_sc


view_dict = {}
view_sc = Shortcut(['store','recall'])
view_dict_sc = Shortcut([])

def zoom(selection="all",buffer=0.0):
   '''
DESCRIPTION
  
   "zoom" scales and translates the window and the origin to cover the
   atom selection.
      
USAGE
 
   zoom object-or-selection [,buffer]
   zoom (selection) [,buffer]
 
PYMOL API

   cmd.zoom( string object-or-selection [,float buffer] )

SEE ALSO

   origin, orient
   '''
   # preprocess selection
   selection = selector.process(selection)
   #   
   try:
      lock()   
      r = _cmd.zoom(str(selection),float(buffer))
   finally:
      unlock()
   return r


clip_action_sc = Shortcut([ 'near','far','move','slab' ])

def clip(mode,offset):
   '''
DESCRIPTION
  
   "clip" alterss the near and far clipping planes according to 
      
USAGE
  
   clip {near|far|move|slab}, distance
 
EXAMPLES
 
   clip near, -5  # moves near plane away from you by 5 A
   clip far, 10   # moves far plane towards you by 10 A
   clip slab, 20  # sets slab thickness to 20 A
   clip move, -5  # moves the slab away from you by 5 A
   
PYMOL API

   cmd.clip( string mode, float distance )

SEE ALSO

   zoom, reset
   '''
   mode = clip_action_sc.auto_err(str(mode),'mode')
   try:
      lock()   
      r = _cmd.clip(str(mode),float(offset))
   finally:
      unlock()
   return r

def origin(selection="(all)",object=None,position=None):
   '''
DESCRIPTION
  
   "origin" sets the center of rotation about a selection.
   If an object name is specified, it can be used to set
   the center of rotation for the object's TTT matrix.
      
USAGE
 
   origin selection [, object [,position]]
   origin (selection)
   origin position=[1.0,2.0,3.0]

PYMOL API
 
   cmd.origin( string object-or-selection )

SEE ALSO

   zoom, orient, reset
   '''
   #'
   # preprocess selection
   selection = selector.process(selection)
   #   
   try:
      lock()
      if object==None: object=''
      if position==None: position=(0.0,0.0,0.0)
      else:
         if cmd.is_string(position):
            position = eval(position)
         selection = ''
      r = _cmd.origin(selection,str(object),
                      (float(position[0]),
                       float(position[1]),
                       float(position[2])
                       ))
   finally:
      unlock()
   return r

def orient(selection="(all)"):
   '''
DESCRIPTION
  
   "orient" aligns the principal components of the atoms in the
   selection with the XYZ axes.  The function is similar to the
   orient command in X-PLOR.
      
USAGE
 
   orient object-or-selection
   orient (selection)
 
PYMOL API
 
   cmd.orient( string object-or-selection )

SEE ALSO

   zoom, origin, reset
   '''
   # preprocess selection
   selection = selector.process(selection)
   #   
   try:
      lock()
      r = _cmd.orient(selection)
   finally:
      unlock()
   return r

   
def move(axis,distance):
   '''
DESCRIPTION
  
   "move" translates the world about one of the three primary axes.
      
USAGE
  
   move axis,distance
 
EXAMPLES
 
   move x,3
   move y,-1
    
PYMOL API
 
   cmd.move( string axis, float distance )

SEE ALSO

   turn
   '''
   try:
      lock()   
      r = _cmd.move(str(axis),float(distance))
   finally:
      unlock()
   return r

def enable(name='all'):
   '''
DESCRIPTION
  
   "enable" enable display of an object and all currently visible representations.
   
USAGE
 
   enable name
   enable all

   name = object or selection name
   
PYMOL API
 
   cmd.enable( string object-name )
 
EXAMPLE
 
   enable my_object

SEE ALSO

   show, hide, disable
   '''
   try:
      lock()   
      r = _cmd.onoff(str(name),1);
   finally:
      unlock()
   return r

def disable(name='all'):
   '''
DESCRIPTION
  
   "disable" disables display of an object and all currently visible
   representations.
   
USAGE
 
   disable name
   disable all 

   "name" is the name of an object or a named selection
   
PYMOL API
 
   cmd.disable( string name ) 
 
EXAMPLE
 
   disable my_object

SEE ALSO

   show, hide, enable   
   '''
   try:
      lock()   
      r = _cmd.onoff(str(name),0);
   finally:
      unlock()
   return r


def show(representation="",selection=""):
   '''
DESCRIPTION
  
   "show" turns on atom and bond representations.
 
   The available representations are:
    
      lines     spheres   mesh      ribbon     cartoon
      sticks    dots      surface   labels
      nonbonded nb_spheres 
   
USAGE
 
   show
   show reprentation [,object]
   show reprentation [,(selection)]
   show (selection)
   
PYMOL API
 
   cmd.show( string representation="", string selection="" )
 
EXAMPLES
 
   show lines,(name ca or name c or name n)
   show ribbon
 
NOTES

   "selection" can be an object name
   "show" alone will turn on lines for all bonds.

SEE ALSO

   hide, enable, disable
   '''
   r=1
   try:
      lock()
      if (representation=="") and (selection==""):
         r = _cmd.showhide("(all)",repres['lines'],1); # show lines by default       
      elif (representation!="") and (selection!=""):
         rep = representation
         rep = repres_sc.auto_err(rep,'representation')
         repn = repres[rep];
         # preprocess selection 
         selection = selector.process(selection)
         #   
         r = _cmd.showhide(str(selection),int(repn),1);
      elif representation=='all':
         r = _cmd.showhide("all",repres['lines'],1); # show lines by default 
      elif (representation[0:1]=='(') or (string.find(representation,'/')>=0):
         # preprocess selection
         selection = selector.process(representation)
         #                  
         r = _cmd.showhide(str(selection),repres['lines'],1);
      else: # selection==""
         rep = representation
         rep = repres_sc.auto_err(rep,'representation')
         repn = repres[rep];
         r = _cmd.showhide("all",int(repn),1);
   finally:
      unlock()
   return r

def hide(representation="",selection=""):
   '''
DESCRIPTION
  
   "hide" turns of atom and bond representations.
 
   The available representations are:
    
      lines     spheres   mesh      ribbon     cartoon
      sticks    dots      surface   labels
      nonbonded nb_spheres
   
USAGE
 
   hide reprentation [,object]
   hide reprentation [,(selection)]
   hide (selection)
   
PYMOL API
 
   cmd.hide( string representation="", string selection="")
 
EXAMPLES
 
   hide lines,all
   hide ribbon

SEE ALSO

   show, enable, disable
   '''
   r = 1
   try:
      lock()
      if (representation=="") and (selection==""):
         r = _cmd.showhide("@",0,0);      
      elif (representation!="") and (selection!=""):
         rep = representation
         rep = repres_sc.auto_err(rep,'representation')
         repn = repres[rep];
         selection = selector.process(selection)
         r = _cmd.showhide(str(selection),int(repn),0);
      elif (representation=='all'):
         r = _cmd.showhide("@",0,0);
      elif (representation[0:1]=='(') or (string.find(representation,'/')>=0):
         selection = selector.process(representation)
         r = _cmd.showhide(str(selection),-1,0);
      else: # selection == ""
         rep = representation
         rep = repres_sc.auto_err(rep,'representation')
         repn = repres[rep];
         r = _cmd.showhide("all",int(repn),0);
   finally:
      unlock()
   return r


def get_view(output=1):
   '''
DESCRIPTION
 
   "get_view" returns and optionally prints out the current view
   information in a format which can be embedded into a command
   script and used in subsequent calls to "set_view"
 
USAGE

   get_view
   
PYMOL API

   cmd.get_view(output=1)  
 
API USAGE

   cmd.get_view(0) # zero option suppresses output
'''
   
   r = None
   try:
      lock()
      r = _cmd.get_view()
   finally:
      unlock()
   if len(r):
      if cmd.get_setting_legacy("logging")!=0.0:
         print " get_view: matrix written to log file."
         cmd.log("_ set_view (\\\n","cmd.set_view((\\\n")
         cmd.log("_  %14.9f, %14.9f, %14.9f,\\\n"%r[0:3]  , "  %14.9f, %14.9f, %14.9f,\\\n"%r[0:3])
         cmd.log("_  %14.9f, %14.9f, %14.9f,\\\n"%r[4:7]  , "  %14.9f, %14.9f, %14.9f,\\\n"%r[4:7])
         cmd.log("_  %14.9f, %14.9f, %14.9f,\\\n"%r[8:11] , "  %14.9f, %14.9f, %14.9f,\\\n"%r[8:11])
         cmd.log("_  %14.9f, %14.9f, %14.9f,\\\n"%r[16:19], "  %14.9f, %14.9f, %14.9f,\\\n"%r[16:19])
         cmd.log("_  %14.9f, %14.9f, %14.9f,\\\n"%r[19:22], "  %14.9f, %14.9f, %14.9f,\\\n"%r[19:22]) 
         cmd.log("_  %14.9f, %14.9f, %14.9f )\n"%r[22:25] , "  %14.9f, %14.9f, %14.9f ))\n"%r[22:25])
         if output<2: # suppress if we have a log file open
            output=0
      if output:
         print "### cut below here and paste into script ###"
         print "set_view (\\"
         print "  %14.9f, %14.9f, %14.9f,\\"%r[0:3]
         print "  %14.9f, %14.9f, %14.9f,\\"%r[4:7]
         print "  %14.9f, %14.9f, %14.9f,\\"%r[8:11]
         print "  %14.9f, %14.9f, %14.9f,\\"%r[16:19]
         print "  %14.9f, %14.9f, %14.9f,\\"%r[19:22]
         print "  %14.9f, %14.9f, %14.9f )"%r[22:25]
         print "### cut above here and paste into script ###"
      
   r = r[0:3]+r[4:7]+r[8:11]+r[16:25]
   return r

def set_view(view):
   '''
DESCRIPTION
 
   "set_view" sets viewing information for the current scene,
   including the rotation matrix, position, origin of rotation,
   clipping planes, and the orthoscopic flag.
 
USAGE

   set_view (...)  where ... is 18 floating point numbers
   
PYMOL API

   cmd.set_view(string-or-sequence view)  
 
'''
   
   r = None
   if cmd.is_string(view):
      try:
         view = eval(view)
      except:
         print "Error: bad view argument; should be a sequence of 18 floats."
         raise QuietException
   if len(view)!=18:
      print "Error: bad view argument; should be a sequence of 18 floats."
      raise QuietException
   else:
      try:
         lock()
         r = _cmd.set_view((
            float(view[ 0]),float(view[ 1]),float(view[ 2]),0.0,
            float(view[ 3]),float(view[ 4]),float(view[ 5]),0.0,
            float(view[ 6]),float(view[ 7]),float(view[ 8]),0.0,
            0.0,0.0,0.0,1.0,
            float(view[ 9]),float(view[10]),float(view[11]),
            float(view[12]),float(view[13]),float(view[14]),
            float(view[15]),float(view[16]),float(view[17])))
      finally:
         unlock()
   return r

def view(key,action='recall'):
   '''
DESCRIPTION
 
   "view" makes it possible to save and restore viewpoints on a given
   scene within a single session.
 
USAGE
 
   view key[,action]
   view ?

   key can be any string
   action should be 'store' or 'recall' (default: 'recall')
   
PYMOL API

   cmd.view(string key,string action)
   
EXAMPLES

   view 0,store
   view 0

SEE ALSO

   get_view
   '''
   if key=='?':
      print " view: stored views:"
      parsing.dump_str_list(view_dict.keys())
   else:
      action = view_sc.auto_err(action,'action')
      if action=='recall':
         key = view_dict_sc.auto_err(key,'view')
         set_view(view_dict[key])
         if _feedback(fb_module.scene,fb_mask.blather): # redundant
            print " view: recalled."
      elif action=='store':
         view_dict_sc.append(key)
         view_dict[key]=cmd.get_view(0)
         if _feedback(fb_module.scene,fb_mask.actions):
            print " view: view stored as '%s'."%key

def stereo(state='on'):
   '''
DESCRIPTION
 
   "stereo" activates or deactives stereo mode.  Currently only
   high-end stereo graphics are supported on the SGI (stereo in a
   window).
 
USAGE
 
   stereo on
   stereo off

PYMOL API

   cmd.stereo(string state="on")
   '''
   state = toggle_dict[toggle_sc.auto_err(str(state),'toggle')]
   r = None
   if state:
      try:
         lock()
         if _cmd.stereo(1):
            r = cmd._stereo(1)
         else:
            print "Error: stereo not available"
      finally:
         unlock();
   else:
      try:
         lock()
         if _cmd.stereo(0):
            r = cmd._stereo(0)
         else:
            print "Error: stereo not available"
      finally:
         unlock();
   return r


def turn(axis,angle):
   '''
DESCRIPTION
  
   "turn" rotates the world about one of the three primary axes
      
USAGE
  
   turn axis, angle
 
EXAMPLES
 
   turn x,90
   turn y,45
    
PYMOL API
 
   cmd.turn( string axis, float angle )

SEE ALSO

   move
   '''
   try:
      lock()
      r = _cmd.turn(str(axis),float(angle))
   finally:
      unlock()
   return r


def full_screen(toggle=1):
   '''
DESCRIPTION

   "full_screen" enables or disables PyMOL's full_screen mode.  This
   is only functions well on PC's.

USAGE

   full_screen on
   full_screen off
   
'''
   toggle = toggle_dict[toggle_sc.auto_err(str(toggle),'toggle')]
   if thread.get_ident() == pymol.glutThread:
      try: 
         lock()
         r = _cmd.full_screen(int(toggle))
      finally:
         unlock()
   else:
      try:
         lock()
         r = _cmd.do("_cmd.full_screen(%d)"%int(toggle))
      finally:
         unlock()
   return r


def rock(mode=-1):
   '''
DESCRIPTION
  
   "rock" toggles Y axis rocking.
 
USAGE
 
   rock
 
PYMOL API
  
   cmd.rock()
   '''
   try:
      lock()   
      r = _cmd.rock(int(mode))
   finally:
      unlock()
   return r

def label(selection="(all)",expression=""):
   '''
DESCRIPTION
 
   "label" labels one or more atoms properties over a selection using
   the python evaluator with a separate name space for each atom.  The
   symbols defined in the name space are:
 
      name, resn, resi, chain, q, b, segi, type (ATOM,HETATM) 
      formal_charge, partial_charge, numeric_type, text_type
   
   All strings in the expression must be explicitly quoted.  This
   operation typically takes several seconds per thousand atoms
   altered.

   To clear labels, simply omit the expression or set it to ''.
   
USAGE

   label (selection),expression
   
EXAMPLES
  
   label (chain A),chain
   label (n;ca),"%s-%s" % (resn,resi)
   label (resi 200),"%1.3f" % partial_charge
   '''
   # preprocess selection
   selection = selector.process(selection)
   #
   r = 1
   try:
      lock()
      if len(str(expression))==0:
         r= _cmd.label("("+str(selection)+")",'')
      else:
         r = _cmd.label("("+str(selection)+")",'label='+str(expression))
   finally:
      unlock()   
   return r


def viewport(width=-1,height=-1):
   '''
DESCRIPTION
  
   "viewport" changes the size of the viewing port (and thus the size
   of all png files subsequently output)
      
USAGE
 
   viewport width, height
 
PYMOL API
  
   cmd.viewport(int width, int height)
   '''
   r = None
   if not cmd.is_glut_thread():
      cmd.do("viewport %d,%d"%(int(width),int(height)))
   else:
      try:
         lock()
         r = _cmd.viewport(int(width),int(height))
      finally:
         unlock()
   return r


def bg_color(color="black"):
   '''
DESCRIPTION
 
   "bg_color" sets the background color

USAGE

   bg_color [color]

PYMOL API

   cmd.color(string color="black")
   
   '''
   r = None
   color = cmd._interpret_color(color)
   try:
      lock()
      r = _cmd.bg_color(str(color))
   finally:
      unlock()
   if not r:
      if _raising(): raise QuietException
   return r

cartoon_dict = {
   'skip'        : -1,
   'automatic'   : 0,
   'loop'        : 1,
   'rectangle'   : 2,
   'oval'        : 3,
   'tube'        : 4,
   'arrow'       : 5,
   'dumbbell'    : 6,   
}

cartoon_sc = Shortcut(cartoon_dict.keys())

def cartoon(type,selection="(all)"):
   '''
DESCRIPTION
  
   "cartoon" changes the default cartoon for a set of atoms.
      
USAGE
 
   cartoon type, (selection)

   type = skip | automatic | loop | rectangle | oval | tube | arrow | dumbbell

PYMOL API
  
   cmd.cartoon(string type, string selection )

EXAMPLES

   cartoon rectangle,(chain A)
   cartoon skip,(resi 145:156)
   
NOTES

   the "automatic" mode utilizes ribbons according to the
   information in the PDB HELIX and SHEET records.

'''
   # preprocess selection
   selection = selector.process(selection)
   #
   type = cartoon_dict[cartoon_sc.auto_err(str(type),'type')];
   r = 1
   try:
      lock()   
      r = _cmd.cartoon("("+str(selection)+")",int(type))
   finally:
      unlock()
   return r

def ray(width=0,height=0,renderer=-1):
   '''
DESCRIPTION
  
   "ray" creates a ray-traced image of the current frame. This
   can take some time (up to several minutes, depending on image
   complexity).
      
USAGE
 
   ray [width,height [,renderer]]

EXAMPLES

   ray
   ray 1024,768
   ray renderer=0
   
PYMOL API
  
   cmd.ray(int width,int height,int renderer=-1)

NOTES

   renderer = -1 is default (use value in ray_default_renderer)
   renderer =  0 uses PyMOL's internal renderer
   renderer =  1 uses PovRay's renderer.  This is Unix-only
      and you must have "x-povray" in your path.  It utilizes two
      two temporary files: "tmp_pymol.pov" and "tmp_pymol.png".
      
SEE ALSO

   "help faster" for optimization tips with the builtin renderer.
   "help povray" for how to use PovRay instead of PyMOL's built-in
      ray-tracing engine.
   
   '''
   try:
      lock()   
      r = _cmd.render(int(width),int(height),int(renderer))
   finally:
      unlock()
   return r

def refresh():
   '''
DESCRIPTION
  
   "refresh" causes the scene to be refresh as soon as it is safe to
   do so.

USAGE

   refresh

PYMOL API
 
   cmd.refresh()

SEE ALSO

   rebuild
   '''
   if thread.get_ident() == pymol.glutThread:
      r = _cmd.refresh_now()
   else:
      try:
         lock()
         r = _cmd.do("cmd._refresh()")
      finally:
         unlock()
   return r

def reset(object=''):
   '''
DESCRIPTION
  
   "reset" restores the rotation matrix to identity, sets the origin
   to the center of mass (approx.) and zooms the window and clipping
   planes to cover all objects.
   
USAGE
 
   reset 
 
PYMOL API
 
   cmd.reset ( )
   '''
   try:
      lock()   
      r = _cmd.reset(0,str(object))
   finally:
      unlock()
   return r


def dirty(): # OBSOLETE?
   try:
      lock()
      r = _cmd.dirty()
   finally:
      unlock()
   return r

def meter_reset():
   '''
DESCRIPTION
  
   "meter_reset" resets the frames per secound counter
   
USAGE

   meter_reset
   '''
   try:
      lock()   
      r = _cmd.reset_rate()
   finally:
      unlock()
   return r

def load_png(filename):
   r=None
   try:
      lock()   
      r = _cmd.load_png(str(filename))
   finally:
      unlock()
   return r


def rebuild(selection='all',representation='everything'):
   '''
DESCRIPTION

   "rebuild" forces PyMOL to recreate geometric objects in
   case any of them have gone out of sync.

USAGE
   
   rebuild [selection [, representation ]]
   
PYMOL API

   cmd.rebuild(string selection = 'all', string representation = 'everything')

SEE ALSO

   refresh
'''
   r = 1
   selection = selector.process(selection)
   representation = repres_sc.auto_err(representation,'representation')
   repn = repres[representation];
   try:
      lock()
      r = _cmd.rebuild(selection,repn)
   finally:
      unlock()

def recolor(selection='all',representation='everything'):
   '''
DESCRIPTION

   "rebuild" forces PyMOL to reapply colors to geometric objects in
   case any of them have gone out of sync.

USAGE
   
   recolor [selection [, representation ]]
   
PYMOL API

   cmd.recolor(string selection = 'all', string representation = 'everything')

SEE ALSO

   recolor
'''
   r = 1
   selection = selector.process(selection)
   representation = repres_sc.auto_err(representation,'representation')
   repn = repres[representation];
   try:
      lock()
      r = _cmd.recolor(selection,repn)
   finally:
      unlock()
         
def color(color,selection="(all)"):
   '''
DESCRIPTION
  
   "color" changes the color of an object or an atom selection.
   
USAGE
 
   color color-name
   color color-name, object-name
   color color-name, (selection)
 
PYMOL API
  
   cmd.color( string color, string color-name )
 
EXAMPLES 
 
   color yellow, (name C*)
   '''
   # preprocess selection
   selection = selector.process(selection)
   color = cmd._interpret_color(color)
   #
   try:
      lock()
      r = _cmd.color(str(color),str(selection),0)
   finally:
      unlock()
   if not r:
      if _raising(): raise QuietException
   return r


def set_color(name,rgb):
   '''
DESCRIPTION
  
   "set_color" defines a new color with color indices (0.0-1.0)
   
USAGE

   set_color name, [ red-float, green-float, blue-float ]
   
   set_color name = [ red-float, green-float, blue-float ]
     # (DEPRECATED)
     
PYMOL API
  
   cmd.set_color( string name, float-list rgb )
 
EXAMPLES 
 
   set_color red = [ 1.0, 0.0, 0.0 ]
   '''
   r = 1
   if cmd.is_string(rgb):
      rgb = eval(rgb)
   if not (isinstance(rgb,types.ListType) or isinstance(rgb,types.TupleType)):
      print "Error: color specification must be a list such as [ 1.0, 0.0, 0.0 ]"
   elif len(rgb)!=3:
      print "Error: color specification must be a list such as [ 1.0, 0.0, 0.0 ]"
   else:
      try:
         lock()

         if len(rgb)==3:
            r = _cmd.colordef(str(name),float(rgb[0]),float(rgb[1]),float(rgb[2]))
            cmd._invalidate_color_sc()
         else:
            print "Error: invalid color."
      finally:
         unlock()
   return r
