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

if __name__=='pymol.viewing':
   
   import thread
   import string
   import types

   import pymol
   import selector
   import copy
   import parsing

   import cmd
   from cmd import _cmd,lock,unlock,Shortcut,QuietException,_raising, \
        _feedback,fb_module,fb_mask, \
        repres,repres_sc, \
        toggle_dict,toggle_sc,stereo_dict,stereo_sc, \
        palette_dict ,palette_sc

   rep_list = [ "lines","sticks","spheres",
                "dots","surface","mesh",
                "nonbonded", "nb_spheres",
                "cartoon","ribbon","labels"]

   view_sc = Shortcut(['store','recall','delete'])
   view_dict = {}
   view_dict_sc = Shortcut([])

   scene_dict = {}
   scene_dict_sc = Shortcut([])

   def zoom(selection="all",buffer=0.0,state=0,complete=0):
      '''
DESCRIPTION

   "zoom" scales and translates the window and the origin to cover the
   atom selection.


USAGE

   zoom [ selection [,buffer [, state [, complete ]]]]

EXAMPLES

   zoom
   zoom complete=1
   zoom (chain A)
   zoom 142/

PYMOL API

   cmd.zoom( string selection, float buffer=0.0,
             int state=0, int complete=0 )

NOTES

   state = 0 (default) use all coordinate states
   state = -1 use only coordinates for the current state
   state > 0  use coordinates for a specific state

   complete = 0 or 1:
      Normally the zoom command tries to guess an optimal zoom level
   for visualization, balancing closeness against occasional clipping
   of atoms out of the field of view.  You can change this behavior by
   setting the complete option to 1, which will guarantee that the
   atom positions for the entire selection will fit in the field of an
   orthoscopic view.  To absolutely prevent clipping, you may also
   need to add a buffer (typically 2 A) to account for the perpective
   transformation and for graphical representations which extend
   beyond the atom coordinates.

SEE ALSO

   origin, orient, center
      '''
      # preprocess selection
      selection = selector.process(selection)
      #   
      try:
         lock()   
         r = _cmd.zoom(str(selection),float(buffer),int(state)-1,int(complete))
      finally:
         unlock()
      return r

   def center(selection="all",state=0,origin=1):
      '''
DESCRIPTION

   "center" translates the window, the clipping slab, and the
   origin to point centered within the atom selection.

USAGE

   center [ selection [,state [, origin]]]

EXAMPLES

   center 145/

PYMOL API

   cmd.center( string selection, int state = 0, int origin = 1 )

NOTES

   state = 0 (default) use all coordinate states
   state = -1 use only coordinates for the current state
   state > 0  use coordinates for a specific state

   origin = 1 (default) move the origin
   origin = 0 leave the origin unchanged

SEE ALSO

   origin, orient, zoom
      '''
      # preprocess selection
      selection = selector.process(selection)
      #   
      try:
         lock()   
         r = _cmd.center(str(selection),int(state)-1,int(origin))
      finally:
         unlock()
      return r

   clip_action_sc = Shortcut([ 'near','far','move','slab','atoms' ])

   def clip(mode,offset,selection=None,state=0):
      '''
DESCRIPTION

   "clip" alters the near and far clipping planes

USAGE

   clip {near|far|move|slab|atoms}, distance [,selection [,state ]]

EXAMPLES

   clip near, -5           # moves near plane away from you by 5 A
   clip far, 10            # moves far plane towards you by 10 A
   clip move, -5           # moves the slab away from you by 5 A
   clip slab, 20           # sets slab thickness to 20 A
   clip slab, 10, resi 11  # clip 10 A slab about residue 11

   clip atoms, 5, pept     # clip atoms in "pept" with a 5 A buffer
                           # about their current camera positions

PYMOL API

   cmd.clip( string mode, float distance, string selection = None)

SEE ALSO

   zoom, reset
      '''
      mode = clip_action_sc.auto_err(str(mode),'mode')
      if selection!=None:
         selection = selector.process(selection)
      else:
         selection = ''
      try:
         lock()   
         r = _cmd.clip(str(mode),float(offset),str(selection),int(state)-1)
      finally:
         unlock()
      return r

   def origin(selection="(all)",object=None,position=None,state=0):
      '''
DESCRIPTION

   "origin" sets the center of rotation about a selection.
   If an object name is specified, it can be used to set
   the center of rotation for the object's TTT matrix.

USAGE

   origin selection [, object [,position, [, state]]]
   origin (selection)
   origin position=[1.0,2.0,3.0]

PYMOL API

   cmd.origin( string object-or-selection )

NOTES

   state = 0 (default) use all coordinate states
   state = -1 use only coordinates for the current state
   state > 0  use coordinates for a specific state

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
                          ),int(state)-1)
      finally:
         unlock()
      return r

   def orient(selection="(all)",state=0):
      '''
DESCRIPTION

   "orient" aligns the principal components of the atoms in the
   selection with the XYZ axes.  The function is similar to the
   orient command in X-PLOR.

USAGE

   orient object-or-selection [, state]
   orient (selection)

PYMOL API

   cmd.orient( string object-or-selection [, state = 0] )

NOTES

   state = 0 (default) use all coordinate states
   state = -1 use only coordinates for the current state
   state > 0  use coordinates for a specific state

SEE ALSO

   zoom, origin, reset
      '''
      # preprocess selection
      selection = selector.process(selection)
      #   
      try:
         lock()
         r = _cmd.orient(selection,int(state)-1)
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
      sticks    dots      surface   labels     extent
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

   def set_view(view,quiet=1):
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
               float(view[15]),float(view[16]),float(view[17])),quiet)
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
   view *

   key can be any string
   action should be 'store' or 'recall' (default: 'recall')

PYMOL API

   cmd.view(string key,string action)

EXAMPLES

   view 0,store
   view 0

SEE ALSO

   set_view, get_view
      '''
      global view_dict,view_dict_sc
   
      if key=='*':
         action = view_sc.auto_err(action,'action')
         if action=='delete':
            view_dict = {}
            view_dict_sc = Shortcut(view_dict.keys())                        
         else:
            print " view: stored views:"
            lst = view_dict.keys()
            lst.sort()
            parsing.dump_str_list(lst)
            
      else:
         action = view_sc.auto_err(action,'action')
         if action=='recall':
            key = view_dict_sc.auto_err(key,'view')
            set_view(view_dict[key])
            if _feedback(fb_module.scene,fb_mask.actions): # redundant
               print " view: \"%s\" recalled."%key
         elif action=='store':
            view_dict_sc.append(key)
            view_dict[key]=cmd.get_view(0)
            if _feedback(fb_module.scene,fb_mask.actions):
               print " view: view stored as \"%s\"."%key
         elif action=='delete':
            key = view_dict_sc.auto_err(key,'view')
            if view_dict.has_key(key):
               del view_dict[key]
               view_dict_sc = Shortcut(view_dict.keys())            
               if _feedback(fb_module.scene,fb_mask.actions): # redundant
                  print " view: '%s' deleted."%key


   def get_vis():
      try:
         lock()
         r = _cmd.get_vis()
      finally:
         unlock()
      return r

   def set_vis(dict):
      try:
         lock()         
         r = _cmd.set_vis(dict)
      finally:
         unlock()
      return r
     
   def scene(key,action='recall'):
      '''
DESCRIPTION

   "scene" makes it possible to save and restore scenes
   scene within a single session.

USAGE

   scene key[,action]
   scene *

   key can be any string
   action should be 'store' or 'recall' (default: 'recall')

PYMOL API

   cmd.scene(string key,string action)

EXAMPLES

   scene 0,store
   scene 0

SEE ALSO

   view, set_view, get_view
      '''
      global scene_dict,scene_dict_sc
   
      if key=='*':
         action = view_sc.auto_err(action,'action')
         if action=='delete':
            scene_dict = {}
            scene_dict_sc = Shortcut(scene_dict.keys())                        
         else:
            print " scene: stored scenes:"
            lst = scene_dict.keys()
            lst.sort()
            parsing.dump_str_list(lst)
            
      else:
         action = view_sc.auto_err(action,'action')
         if action=='recall':
            key = scene_dict_sc.auto_err(key,'scene')
            set_view(scene_dict[key][0])
            cmd.hide()
            cmd.disable()
            cmd.set_vis(scene_dict[key][1])
            cmd.frame(scene_dict[key][2])
            for rep in rep_list:
               name = "_scene_"+key+"_"+rep
               cmd.show(rep,name)
            if _feedback(fb_module.scene,fb_mask.actions): # redundant
               print " scene: \"%s\" recalled."%key
         elif action=='store':
            scene_dict_sc.append(key)
            scene_dict[key]=[cmd.get_view(0),
                             cmd.get_vis(),
                             cmd.get_frame()]
            for rep in rep_list:
               name = "_scene_"+key+"_"+rep
               cmd.select(name,"rep "+rep)
            if _feedback(fb_module.scene,fb_mask.actions):
               print " scene: scene stored as \"%s\"."%key
         elif action=='delete':
            key = scene_dict_sc.auto_err(key,'view')
            if scene_dict.has_key(key):
               del scene_dict[key]
               scene_dict_sc = Shortcut(scene_dict.keys())            
               name = "_scene_"+key+"_*"
               cmd.delete(name)
               if _feedback(fb_module.scene,fb_mask.actions):
                  print " scene: '%s' deleted."%key

               
   def session_save_views(session):
      session['view_dict']=copy.deepcopy(view_dict)
      return 1

   def session_restore_views(session):
      global view_dict,view_dict_sc
      if session.has_key('view_dict'):
         view_dict=copy.deepcopy(session['view_dict'])
         view_dict_sc = Shortcut(view_dict.keys())
      return 1

   if session_restore_views not in pymol._session_restore_tasks:
      pymol._session_restore_tasks.append(session_restore_views)

   if session_save_views not in pymol._session_save_tasks:
      pymol._session_save_tasks.append(session_save_views)

   def session_save_scenes(session):
      session['scene_dict']=copy.deepcopy(scene_dict)
      return 1

   def session_restore_scenes(session):
      global scene_dict,scene_dict_sc
      if session.has_key('scene_dict'):
         scene_dict=copy.deepcopy(session['scene_dict'])
         scene_dict_sc = Shortcut(scene_dict.keys())
      return 1

   if session_restore_scenes not in pymol._session_restore_tasks:
      pymol._session_restore_tasks.append(session_restore_scenes)

   if session_save_scenes not in pymol._session_save_tasks:
      pymol._session_save_tasks.append(session_save_scenes)

   def stereo(state='on',quiet=1):
      '''
DESCRIPTION

   "stereo" activates or deactives stereo mode.

USAGE

   stereo on
   stereo off
   stereo swap
   stereo crosseye 
   stereo quadbuffer

NOTES

   quadbuffer is the default stereo mode if hardware stereo is available
   otherwise, crosseye is the default

PYMOL API

   cmd.stereo(string state="on")
      '''
      state = stereo_dict[stereo_sc.auto_err(str(state),'state')]
      r = None
      try:
         lock()
         if state>1:
            if state==2: # cross-eye
               cmd.set("stereo_mode","2",quiet=quiet)
            elif state==3: # quad
               cmd.set("stereo_mode","1",quiet=quiet)
            elif state==4: # wall-eye
               cmd.set("stereo_mode","3",quiet=quiet)
            state=1
         if not _cmd.stereo(state):
            print "Error: Selected stereo mode is not available."
            if _raising(): raise QuietException
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

   def label(selection="(all)",expression="",quiet=1):
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
            r= _cmd.label("("+str(selection)+")",'',quiet)
         else:
            r = _cmd.label("("+str(selection)+")",'label='+str(expression),quiet)
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

   def ray(width=0,height=0,renderer=-1,angle=0.0,shift=0.0):
      '''
DESCRIPTION

   "ray" creates a ray-traced image of the current frame. This
   can take some time (up to several minutes, depending on image
   complexity).

USAGE

   ray [width,height [,renderer [,angle [,shift ]]]

   angle and shift can be used to generate matched stereo pairs
   
EXAMPLES

   ray
   ray 1024,768
   ray renderer=0

PYMOL API

   cmd.ray(int width,int height,int renderer=-1,float shift=0)

NOTES

   renderer = -1 is default (use value in ray_default_renderer)
   renderer =  0 uses PyMOL's internal renderer
   renderer =  1 uses PovRay's renderer.  This is Unix-only
      and you must have "x-povray" in your path.  It utilizes two
      two temporary files: "tmp_pymol.pov" and "tmp_pymol.png".

SEE ALSO

   "help faster" for optimization tips with the builtin renderer.
   "help povray" for how to use PovRay instead of PyMOL\'s built-in
   ray-tracing engine.

      '''
      try:
         lock()   
         r = _cmd.render(int(width),int(height),
                         int(renderer),float(angle),float(shift))
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

   def load_png(filename,movie=1,quiet=0):
      r=None
      try:
         lock()
         r = _cmd.load_png(str(filename),int(movie),int(quiet))
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

   def color(color,selection="(all)",quiet=1):
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
         r = _cmd.color(str(color),str(selection),0,int(quiet))
      finally:
         unlock()
      if not r:
         if _raising(): raise QuietException
      return r

   def spectrum(expression="count",
                palette="rainbow",
                selection="(all)",
                minimum=None,
                maximum=None,
                byres=0,quiet=1):
      '''
DESCRIPTION

   "spectrum" colors atoms using a spectrum
   
USAGE

PYMOL API

EXAMPLES 

      '''

      palette = palette_sc.auto_err(palette,'palette')
      
      (prefix,digits,first,last) = palette_dict[str(palette)]

      if (maximum==None) or (minimum==None):
         minimum = 0 # signal to auto-adjust levels
         maximum = -1

      # preprocess selection
      selection = selector.process(selection)
      #
      try:
         lock()
         r = _cmd.spectrum(str(selection),str(expression),
                           float(minimum),float(maximum),
                           int(first),int(last),str(prefix),
                           int(digits),int(byres),int(quiet))
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





