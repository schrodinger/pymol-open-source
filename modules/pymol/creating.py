#A* -------------------------------------------------------------------
#B* This file contains source code for the PyMOL computer program
#C* copyright 1998-2000 by Warren Lyford Delano of DeLano Scientific. 
#D* -------------------------------------------------------------------
#E* It is unlawful to modify or remove this copyright notice.
#F* -------------------------------------------------------------------
#G* Please see the accompanying LICENSE file for further information. 
#H* -------------------------------------------------------------------
#I* Additional authors of this source file include:
#-* Filipe Maia (slicing code)
#-* 
#-*
#Z* -------------------------------------------------------------------

if __name__=='pymol.creating':
   
   import selector
   import traceback
   import operator
   import cmd
   import string
   from cmd import _cmd,lock,unlock,Shortcut,QuietException,is_list,is_string
   from cmd import file_ext_re, sanitize_list_re
   from chempy import fragments

   map_type_dict = {
      'vdw' : 0,
      'coulomb' : 1,
      'gaussian' : 2,
      'coulomb_neutral' : 3,
      'coulomb_local' : 4,
      }

   map_type_sc = Shortcut(map_type_dict.keys())

   ramp_spectrum_dict = {
      "traditional" : 1,
      "sludge" : 2,
      "ocean" : 3,
      "hot" : 4,
      "grayable" : 5,
      "rainbow" : 6,
      "afmhot" : 7,
      "grayscale" : 8,
      }
   
   ramp_spectrum_sc = Shortcut(ramp_spectrum_dict.keys())
   
   def map_new(name,type='gaussian',grid=None,selection="(all)",buffer=None,box=None,state=0):
      '''
      state > 0: do indicated state
      state = 0: independent states in independent extents
      state = -1: current state
      state = -2: independent states in unified extent
      state = -3: all states in one map
      '''
      # preprocess selection
      r = None
      selection = selector.process(selection)
      if box!=None: # box should be [[x1,y1,z1],[x2,y2,z2]]
         if cmd.is_string(box):
            box = eval(box)
         box = (float(box[0][0]),
                float(box[0][1]),
                float(box[0][2]),
                float(box[1][0]),
                float(box[1][1]),
                float(box[1][2]))
      else:
         box = (0.0,0.0,0.0,1.0,1.0,1.0)
      if grid==None:
         grid = cmd.get_setting_legacy('gaussian_resolution')/3.0
      if buffer==None:
         buffer = cmd.get_setting_legacy('gaussian_resolution')
      grid = float(grid) # for now, uniform xyz; later (x,y,z)
      type = map_type_dict[map_type_sc.auto_err(str(type),'map type')]
      try:
         lock()
         r = _cmd.map_new(str(name),int(type),grid,str(selection),
                          float(buffer),box,int(state)-1)
      finally:
         unlock()

   def ramp_new(name,map_name,range=[-1.0,0.0,1.0],
                color=['red',[1.0,1.0,1.0],'blue'],
                map_state=1,selection='',
                beyond=2.0,within=6.0,
                sigma=2.0,zero=1):
      safe_color = string.strip(str(color))
      if(safe_color[0:1]=="["): # looks like a list
         color = eval(str(safe_color))
      else: # looks like a literal
         color = str(color)
      new_color = []
      # preprocess selection
      if selection!='':
         selection = selector.process(selection)
      if is_list(color):
         for a in color:
            if not is_list(a):
               new_color.append(list(cmd.get_color_tuple(a)))
            else:
               new_color.append(a)
      elif is_string(color):
         new_color = ramp_spectrum_dict[ramp_spectrum_sc.auto_err(str(color),'ramp color spectrum')]
      else:
         new_color=int(color)
      try:
         lock()
         r = _cmd.ramp_new(str(name),str(map_name),list(eval(str(range))),new_color,
                           int(map_state)-1,str(selection),float(beyond),float(within),
                           float(sigma),int(zero))
      finally:
         unlock()

   def isomesh(name,map,level=1.0,selection='',buffer=0.0,state=1,carve=None,source_state=0):
      '''
DESCRIPTION

   "isomesh" creates a mesh isosurface object from a map object.

USAGE

   isomesh name, map, level [,(selection) [,buffer [,state [,carve ]]]]

   name = the name for the new mesh isosurface object.

   map = the name of the map object to use for computing the mesh.

   level = the contour level.

   selection = an atom selection about which to display the mesh with
      an additional "buffer" (if provided).

   state = the state into which the object should be loaded (default=1)
      (set state=0 to append new mesh as a new state)

   carve = a radius about each atom in the selection for which to
      include density. If "carve" is not provided, then the whole
      brick is displayed.

NOTES

   If the mesh object already exists, then the new mesh will be
   appended onto the object as a new state (unless you indicate a state).

   state > 0: specific state
   state = 0: all states
   state = -1: current state
   
   source_state > 0: specific state
   source_state = 0: include all states starting with 0
   source_state = -1: current state
   source_state = -2: last state in map

SEE ALSO

   isodot, load
'''

      if selection!='':
         mopt = 1 # about a selection
      else:
         mopt = 0 # render the whole map
      # preprocess selection
      selection = selector.process(selection)
      if selection not in [ 'center', 'origin' ]:
         selection = "("+selection+")"
      #
      if carve==None:
         carve=0.0
      try:
         lock()
         r = _cmd.isomesh(str(name),0,str(map),int(mopt),
                          selection,float(buffer),
                          float(level),0,int(state)-1,float(carve),
                          int(source_state)-1)
      finally:
         unlock()
      return r

   def slice_new(name,map,state=1,source_state=0):
      '''
DESCRIPTION

   "slice_map" creates a slice object from a map object.

USAGE

   slice_map name, map, [opacity, [resolution, [state, [source_state]]]]

   name = the name for the new slice object.

   map = the name of the map object to use for computing the slice.

   opacity = opacity of the new slice (default=1)

   resolution = the number of pixels per sampling (default=5)

   state = the state into which the object should be loaded (default=1)
      (set state=0 to append new mesh as a new state)

   source_state = the state of the map from which the object should be loaded (default=0)
   
SEE ALSO

   isomesh, isodot, load
'''

      try:
         lock()
         r = _cmd.slice_new(str(name),str(map),int(state)-1,int(source_state)-1)
      finally:
         unlock()
      return r


   def isosurface(name,map,level=1.0,selection='',buffer=0.0,state=1,carve=None,
                  source_state=0,side=1,mode=3):
      '''
DESCRIPTION

   "isosurface" creates a new surface object from a map object.

USAGE

   isosurface name, map, level [,(selection) [,buffer [,state [,carve ]]]]

   name = the name for the new mesh isosurface object.

   map = the name of the map object to use for computing the mesh.

   level = the contour level.

   selection = an atom selection about which to display the mesh with
      an additional "buffer" (if provided).

   state = the state into which the object should be loaded (default=1)
      (set state=0 to append new surface as a new state)

   carve = a radius about each atom in the selection for which to
      include density. If "carve= not provided, then the whole
      brick is displayed.

NOTES

   If the surface object already exists, then the new surface will be
   appended onto the object as a new state (unless you indicate a state).

SEE ALSO

   isodot, isomesh, load
      '''
      if selection!='':
         mopt = 1 # about a selection
      else:
         mopt = 0 # render the whole map
      # preprocess selection
      selection = selector.process(selection)
      if selection not in [ 'center', 'origin' ]:
         selection = "("+selection+")"
     #
      if carve==None:
         carve=0.0
      try:
         lock()
         r = _cmd.isosurface(str(name),0,str(map),int(mopt),
                             selection,float(buffer),
                             float(level),int(mode),int(state)-1,float(carve),
                             int(source_state)-1,int(side))
      finally:
         unlock()
      return r

   def isodot(name,map,level=1.0,selection='',buffer=0.0,state=0,
              carve=None,source=0):
      '''
DESCRIPTION

   "isodot" creates a dot isosurface object from a map object.

USAGE

   isodot name = map, level [,(selection) [,buffer [, state ] ] ] 

   map = the name of the map object to use.

   level = the contour level.

   selection = an atom selection about which to display the mesh with
      an additional "buffer" (if provided).

NOTES

   If the dot isosurface object already exists, then the new dots will
   be appended onto the object as a new state.

SEE ALSO

   load, isomesh
      '''
      if selection!='':
         mopt = 1 # about a selection
      else:
         mopt = 0 # render the whole map
      # preprocess selections
      selection = selector.process(selection)
      if selection not in [ 'center', 'origin' ]:
         selection = "("+selection+")"
      #
      if carve==None:
         carve=0.0
      try:
         lock()
         r = _cmd.isomesh(str(name),0,str(map),int(mopt),
                          selection,float(buffer),
                          float(level),1,int(state)-1,
                          float(carve),int(source)-1)
      finally:
         unlock()
      return r

   def isolevel(name,level=1.0,state=0):
      '''
DESCRIPTION

   "isolevel" changes the contour level of a isodot, isosurface, or isomesh object.

USAGE

   isolevel name, level, state

      '''
      try:
         lock()
         r = _cmd.isolevel(str(name),float(level),int(state)-1)
      finally:
         unlock()
      return r

   def copy(target,source):
      '''
DESCRIPTION

   "copy" creates a new object that is an identical copy of an
   existing object

USAGE

   copy target, source

   copy target = source         # (DEPRECATED)

PYMOL API

   cmd.copy(string target,string source)

SEE ALSO

   create
      '''
      try:
         lock()
         r = _cmd.copy(str(source),str(target))
      finally:
         unlock()
      return r

   def symexp(prefix,object,selection,cutoff):
      '''
DESCRIPTION

   "symexp" creates all symmetry related objects for the specified object
   that occurs within a cutoff about an atom selection.  The new objects
   are labeled using the prefix provided along with their crystallographic
   symmetry operation and translation.

USAGE

   symexp prefix = object, (selection), cutoff

PYMOL API

   cmd.symexp( string prefix, string object, string selection, float cutoff) 

SEE ALSO

   load
      '''
      # preprocess selection
      selection=selector.process(selection)
      #
      try:
         lock()
         r = _cmd.symexp(str(prefix),str(object),"("+str(selection)+")",float(cutoff))
      finally:
         unlock()
      return r

   def fragment(name,object=None,origin=1,zoom=0,quiet=1):
      '''
DESCRIPTION

   "fragment" retrieves a 3D structure from the fragment library, which is currently
   pretty meager (just amino acids).

USAGE

   fragment name

   '''
      r = 1
      try:
         save=cmd.get_setting_legacy('auto_zoom')
         cmd.set('auto_zoom',zoom,quiet=1)
         try:
            if object==None:
               object=name
            model = fragments.get(str(name))
            la = len(model.atom)
            if la:
               mean = map(lambda x,la=la:x/la,[
                  reduce(operator.__add__,map(lambda a:a.coord[0],model.atom)),

                  reduce(operator.__add__,map(lambda a:a.coord[1],model.atom)),
                  reduce(operator.__add__,map(lambda a:a.coord[2],model.atom))])
               position = cmd.get_position()
               for c in range(0,3):
                  mean[c]=position[c]-mean[c]
                  map(lambda a,x=mean[c],c=c:cmd._adjust_coord(a,c,x),model.atom)
               mean = map(lambda x,la=la:x/la,[
                  reduce(operator.__add__,map(lambda a:a.coord[0],model.atom)),
                  reduce(operator.__add__,map(lambda a:a.coord[1],model.atom)),
                  reduce(operator.__add__,map(lambda a:a.coord[2],model.atom))])
            cmd.load_model(model,str(object),quiet=quiet)
         finally:
            cmd.set('auto_zoom',save,quiet=1)
      except:
         traceback.print_exc()
         print "Error: unable to load fragment %s" % name
         r = 0
      if not r:
         if cmd._raising(): raise QuietException
      return r

   def create(name,selection,source_state=0,target_state=0,discrete=0):
      '''
DESCRIPTION

   "create" creates a new molecule object from a selection.  It can
   also be used to create states in an existing object.

   NOTE: this command has not yet been throughly tested.

USAGE

   create name, (selection) [,source_state [,target_state ] ]

   create name = (selection) [,source_state [,target_state ] ]
     # (DEPRECATED)

   name = object to create (or modify)
   selection = atoms to include in the new object
   source_state (default: 0 - copy all states)
   target_state (default: 0)

PYMOL API

   cmd.create(string name, string selection, int state, int target_state,int discrete)

NOTES

   If the source and target states are zero (default), all states will
   be copied.  Otherwise, only the indicated states will be copied.

SEE ALSO

   load, copy
      '''
      # preprocess selection
      selection = selector.process(selection)
      #      
      try:
         lock()
         if name==None:
            sel_cnt = _cmd.get("sel_counter") + 1.0
            _cmd.legacy_set("sel_counter","%1.0f" % sel_cnt)
            name = "obj%02.0f" % sel_cnt
         _cmd.create(str(name),"("+str(selection)+")",
                     int(source_state)-1,int(target_state)-1,int(discrete))
      finally:
         unlock()
      return None





