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

    import pymol
    import selector
    import traceback
    import operator
    import cmd
    import string
    import re
    from cmd import _cmd, Shortcut, is_list, is_string, \
          file_ext_re, safe_list_eval, safe_alpha_list_eval, \
          DEFAULT_ERROR, DEFAULT_SUCCESS, is_ok, is_error, \
          is_tuple

    from chempy import fragments

    map_type_dict = {
        'vdw' : 0,
        'coulomb' : 1,
        'gaussian' : 2, # gaussian summation
        'coulomb_neutral' : 3,
        'coulomb_local' : 4,
        'gaussian_max' : 5, # gaussian maximum contributor
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
        "object" : [[-1.0,-1.0,-1.0]]
        }
    
    ramp_spectrum_sc = Shortcut(ramp_spectrum_dict.keys())

    group_action_dict = {
        "add" : 1,
        "remove" : 2,
        "open" : 3,
        "close" : 4,
        "toggle" : 5,
        "auto" : 6,
        "ungroup" : 7,
        "empty" : 8,
        "purge" : 9,
        "excise" : 10,
        }

    group_action_sc =  Shortcut(group_action_dict.keys())
    
    def group(name, members="", action='auto', quiet=1,_self=cmd):
        '''

DESCRIPTION

    "group" creates or updates a group object: a container for
    organizing objects into a hierarchy.
    
USAGE

    group name [, members [, action ]]

ARGUMENTS

    name = string: name of the group

    members = string: space-separated list of objects to include in
              the group

    action = add, remove, open, close, toggle, auto, ungroup, empty,
             purge, excise

EXAMPLE

    group kinases, 1oky 1pkg 1t46 1uwh 1z5m
    group kinases, open
    group kinases, close

NOTES

    Group objects can typically be used as arguments to commands.  In
    such cases, the command should be applied to all members of the
    group.  If the group is used as a selection, then all atoms in all
    objects in the group should be included in the selection.

    When a group objects is open, objects can be added or removed from
    the group by right-clicking and dragging in the control panel.

SEE ALSO

    ungroup, order
    
'''
        
        r = DEFAULT_ERROR        
        action = group_action_dict[group_action_sc.auto_err(str(action),'group action')]
        if name=='all': name='*'
        if action==6:
            if len(members):
                action=1
            elif (name in _self.get_names()) or ('*' in name):
                action=5
            else:
                action=1
        try:
            _self.lock(_self)
            r = _cmd.group(_self._COb,str(name),str(members),int(action),int(quiet))
        finally:
            _self.unlock(r,_self)
        if _self._raising(r,_self): raise pymol.CmdException         
        return r
    
    def ungroup(name, members="", quiet=1, _self=cmd):
        '''

DESCRIPTION

    "ungroup" removes an object from a group object, returning it to
    the top level.

USAGE

    ungroup name


SEE ALSO

    group
    
    '''
        
        r = DEFAULT_ERROR
        try:
            _self.lock(_self)
            r = _cmd.group(_self._COb,str(name),str(members),7,int(quiet))
        finally:
            _self.unlock(r,_self)
        if _self._raising(r,_self): raise pymol.CmdException         
        return r
    
    def map_new(name, type='gaussian', grid=None, selection="(all)",
                buffer=None, box=None, state=0, quiet=1, zoom=0,
                normalize=-1, clamp=[1.0,-1.0], _self=cmd):

        '''

DESCRIPTION

    "map_new" creates a map object using one of the built-in map
    generation routines.  This command not yet fully supported.

USAGE

    map_new name [, type [, grid [, selection [, buffer [, box [, state ]]]]]]

ARGUMENTS

    name = string: name of the map object to create or modify
	
	type = vdw, gaussian, gaussian_max, coulomb, coulomb_neutral, coulomb_local

    grid = float: grid spacing

    selection = string: atoms about which to generate the map

    buffer = float: cutoff 
    
    state > 0: use the indicated state
    
    state = 0: use all states independently with independent extents
    
    state = -1: use current global state
    
    state = -2: use effective object state(s)
    
    state = -3: use all states in one map
    
    state = -4: use all states independent states by with a unified extent

NOTES

    This command can be used to create low-resolution surfaces of
    protein structures.
    
    '''
        # preprocess selection
        r = DEFAULT_ERROR
        selection = selector.process(selection)
        if box!=None: # box should be [[x1,y1,z1],[x2,y2,z2]]
            if _self.is_string(box):
                box = safe_list_eval(box)
            box = (float(box[0][0]),
                   float(box[0][1]),
                   float(box[0][2]),
                   float(box[1][0]),
                   float(box[1][1]),
                   float(box[1][2]))
            box_flag = 1
        else:
            box = (0.0,0.0,0.0,1.0,1.0,1.0)
            box_flag = 0
        if grid==None:
            grid = _self.get_setting_legacy('gaussian_resolution')/3.0
        if buffer==None:
            buffer = _self.get_setting_legacy('gaussian_resolution')
        grid = float(grid) # for now, uniform xyz; later (x,y,z)

        if not is_list(clamp):
            clamp = safe_list_eval(str(clamp))
        if len(clamp)<2:
            clamp = [1.0,-1.0]
        type = map_type_dict[map_type_sc.auto_err(str(type),'map type')]
        try:
            _self.lock(_self)
            r = _cmd.map_new(_self._COb,str(name),int(type),grid,str(selection),
                             float(buffer),box,int(state)-1,
                             int(box_flag),int(quiet),int(zoom),int(normalize),
                             float(clamp[0]),float(clamp[1]))
        finally:
            _self.unlock(r,_self)
        if _self._raising(r,_self): raise pymol.CmdException         
        return r

    def ramp_new(name, map_name, range=[-1.0,0.0,1.0],
                 color=['red',[1.0,1.0,1.0],'blue'], state=0,
                 selection='', beyond=2.0, within=6.0, sigma=2.0,
                 zero=1, quiet=1, _self=cmd):

        '''
DESCRIPTION

    "ramp_new" creates a color ramp based on a map potential value or
    based on proximity to a molecular object.
    
USAGE

    ramp_new name, map_name [, range [, color [, state [, selection [,
        beyond [, within [, sigma [, zero ]]]]]]]]

ARGUMENTS

    name = string: name of the ramp object

    map_name = string: name of the map (for potential) or molecular
    object (for proximity)
    
    range = list: values corresponding to slots in the ramp

    color = list: colors corresponding to slots in the ramp

    state = integer: state identifier

    selection = selection: for automatic ranging
    
    beyond = number: with automatic ranging, are we excluding
    values beyond a certain distance from the selection?

    within = number: with automatic ranging, are we only including
    valuess within a certain distance from the selection?

    sigma = number: with automatic ranging, how many standard
    deviations from the mean do we go?

    zero = integer: with automatic ranging, do we force the central
    value to be zero?

EXAMPLES

    ramp_new e_pot_color, e_pot_map, [-10,0,10], [red,white,blue]

NOTES

    Color ramps are extremely powerful but complicated to use.

    In the simplest case, they can be used to color representations
    based on the potential values found in a map object at the
    corresponding positions in space.

    In another simple case, representations can be colored based on
    proximity to a target.  Note that since ramp targets must
    themselves be real objects (not merely selections), the "create"
    command may be needed in order to generate an appropriate target.
    
    In more complicated cases, they can be used to color
    representations on one object based atoms found in another.

    Ramps can operate recursively.  In other words, the output color
    from one ramp can be used as the input color for another.  For
    example, you could color by map potential within a certain
    distance of the target object, beyond which, a uniform color is applied.
    
    
PYMOL API

    def ramp_new(string name, string map_name, list range, list color,
                 int state, string selection, float beyond, float
                 within, float sigma, int zero, int quiet)

SEE ALSO

    load, color, create, slice, gradient
    
    '''
        r = DEFAULT_ERROR
        safe_color = string.strip(str(color))
        if(safe_color[0:1]=="["): # looks like a list
            color = safe_alpha_list_eval(str(safe_color))
        else: # looks like a literal
            color = str(color)
        new_color = []
        # preprocess selection
        if selection!='':
            selection = selector.process(selection)
        # coerce range
        range = list(safe_list_eval(str(range)))
        range = map(lambda x:float(x),range)
        if is_list(color):
            for a in color:
                if not is_list(a):
                    new_color.append(list(_self.get_color_tuple(a,4))) # incl negative RGB special colors
                else:
                    new_color.append(a)
        elif is_string(color):
            new_color = ramp_spectrum_dict[ramp_spectrum_sc.auto_err(str(color),'ramp color spectrum')]
        else:
            new_color=int(color)
        try:
            _self.lock(_self)
            r = _cmd.ramp_new(_self._COb,str(name),str(map_name),range,new_color,
                                    int(state)-1,str(selection),float(beyond),float(within),
                                    float(sigma),int(zero),int(quiet))
        finally:
            _self.unlock(r,_self)
        if _self._raising(r,_self): raise pymol.CmdException         
        return r

    def isomesh(name, map, level=1.0, selection='', buffer=0.0,
                state=1, carve=None, source_state=0, quiet=1, _self=cmd):
        '''
DESCRIPTION

    "isomesh" creates a mesh isosurface object from a map object.

USAGE

    isomesh name, map, level [, selection [, buffer [, state [, carve ]]]]

ARGUMENTS

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
        r = DEFAULT_ERROR
        if selection!='':
            region = 1 # about a selection
        else:
            region = 0 # render the whole map
        # preprocess selection
        selection = selector.process(selection)
        if selection not in [ 'center', 'origin' ]:
            selection = "("+selection+")"
        #
        if carve==None:
            carve=0.0
        try:
            _self.lock(_self)
            r = _cmd.isomesh(_self._COb,str(name),str(map),int(region),
                             selection,float(buffer),
                             float(level),0,int(state)-1,float(carve),
                             int(source_state)-1,int(quiet),
                             float(level))
        finally:
            _self.unlock(r,_self)
        if _self._raising(r,_self): raise pymol.CmdException         
        return r

    def slice_new(name, map, state=1, source_state=0, _self=cmd):
        '''
DESCRIPTION

    "slice_map" creates a slice object from a map object.

USAGE

    slice_map name, map, [opacity, [resolution, [state, [source_state]]]]

ARGUMENTS

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
        r = DEFAULT_ERROR
        try:
            _self.lock(_self)
            r = _cmd.slice_new(_self._COb,str(name),str(map),int(state)-1,int(source_state)-1)
        finally:
            _self.unlock(r,_self)
        if _self._raising(r,_self): raise pymol.CmdException         
        return r


    def isosurface(name, map, level=1.0, selection='', buffer=0.0, state=1,
                   carve=None, source_state=0, side=1, mode=3, quiet=1,
                   _self=cmd):
        '''
DESCRIPTION

    "isosurface" creates a new surface object from a map object.

USAGE

    isosurface name, map, level [, selection [, buffer [, state [, carve ]]]]

ARGUMENTS

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
        r = DEFAULT_ERROR
        if selection!='':
            region = 1 # about a selection
        else:
            region = 0 # render the whole map
        # preprocess selection
        selection = selector.process(selection)
        if selection not in [ 'center', 'origin' ]:
            selection = "("+selection+")"
      #
        if carve==None:
            carve=0.0
        try:
            _self.lock(_self)
            r = _cmd.isosurface(_self._COb,str(name),str(map),int(region),
                                      selection,float(buffer),
                                      float(level),int(mode),int(state)-1,float(carve),
                                      int(source_state)-1,int(side),int(quiet))

        finally:
            _self.unlock(r,_self)
        if _self._raising(r,_self): raise pymol.CmdException         
        return r

    def isodot(name,map,level=1.0,selection='',buffer=0.0,state=0,
                  carve=None,source_state=0,quiet=1,_self=cmd):
        '''
DESCRIPTION

    "isodot" creates a dot isosurface object from a map object.

USAGE

    isodot name = map, level [,(selection) [,buffer [, state ] ] ] 

ARGUMENTS

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
        r = DEFAULT_ERROR
        if selection!='':
            region = 1 # about a selection
        else:
            region = 0 # render the whole map
        # preprocess selections
        selection = selector.process(selection)
        if selection not in [ 'center', 'origin' ]:
            selection = "("+selection+")"
        #
        if carve==None:
            carve=0.0
        try:
            _self.lock(_self)
            r = _cmd.isomesh(_self._COb,str(name),str(map),int(region),
                             selection,float(buffer),
                             float(level),1,int(state)-1,
                             float(carve),int(source_state)-1,int(quiet),
                             float(level))
        finally:
            _self.unlock(r,_self)
        if _self._raising(r,_self): raise pymol.CmdException                  
        return r



    def isolevel(name,level=1.0,state=0,query=0,quiet=1,_self=cmd):
        '''
DESCRIPTION

    "isolevel" changes the contour level of a isodot, isosurface, or isomesh object.

USAGE

    isolevel name, level, state

        '''
        r = DEFAULT_ERROR
        try:
            _self.lock(_self)
            r = _cmd.isolevel(_self._COb,str(name),float(level),int(state)-1,int(query),int(quiet))
        finally:
            _self.unlock(r,_self)
        if not int(query):
            if _self._raising(r,_self): raise pymol.CmdException                  
        return r

    def gradient(name, map, minimum=1.0, maximum=-1.0,
                 selection='', buffer=0.0, state=0,
                 carve=None, source_state=0, quiet=1, _self=cmd):
        '''
DESCRIPTION

    "gradient" creates a gradient object from a map object.

USAGE

    gradient name = map, [ minimum, [, maximum [, selection [, buffer [, state ]]]]]

ARGUMENTS

    map = the name of the map object to use.

    minimum, maximum = minimum and maximum levels (default: full map range)

    selection = an atom selection about which to display the mesh with
        an additional "buffer" (if provided).

SEE ALSO

    load, isomesh
        '''
        r = DEFAULT_ERROR
        if selection!='':
            region = 1 # about a selection
        else:
            region = 0 # render the whole map
        # preprocess selections
        selection = selector.process(selection)
        if selection not in [ 'center', 'origin' ]:
            selection = "("+selection+")"
        #
        if carve==None:
            carve=0.0
        try:
            _self.lock(_self)
            r = _cmd.isomesh(_self._COb,str(name),str(map),int(region),
                             selection,float(buffer),
                             float(minimum),3,int(state)-1,
                             float(carve),int(source_state)-1,int(quiet),
                             float(maximum))
        finally:
            _self.unlock(r,_self)
        if _self._raising(r,_self): raise pymol.CmdException                  
        return r

    def copy(target,source,zoom=-1,_self=cmd):
        '''
DESCRIPTION

    "copy" creates a new object that is an identical copy of an
    existing object.

USAGE

    copy target, source

NOTES

    Currently, this command only works for molecular objects.

PYMOL API

    cmd.copy(string target, string source)

SEE ALSO

    create
        '''
        r = DEFAULT_ERROR      
        try:
            _self.lock(_self)
            r = _cmd.copy(_self._COb,str(source),str(target),int(zoom))
        finally:
            _self.unlock(r,_self)
        if _self._raising(r,_self): raise pymol.CmdException                  
        return r

    def symexp(prefix, object, selection, cutoff, segi=0, quiet=1,_self=cmd):
        '''
DESCRIPTION

    "symexp" creates all symmetry-related objects for the specified
    object that occur within a cutoff about an atom selection.

USAGE

    symexp prefix, object, selection, cutoff

NOTES

    Thee newly objects are labeled using the prefix provided along with
    their crystallographic symmetry operation and translation.

SEE ALSO

    load
        '''
        r = DEFAULT_ERROR
        # preprocess selection
        selection=selector.process(selection)
        #
        try:
            _self.lock(_self)
            r = _cmd.symexp(_self._COb,str(prefix),str(object),
                            "("+str(selection)+")",float(cutoff),
                            int(segi),int(quiet))
        finally:
            _self.unlock(r,_self)
        if _self._raising(r,_self): raise pymol.CmdException                           
        return r

    def fragment(name, object=None, origin=1, zoom=0, quiet=1, _self=cmd):
        '''
DESCRIPTION

    "fragment" retrieves a 3D structure from the fragment library,
    which is currently pretty meager (just amino acids).

USAGE

    fragment name

    '''
        r = DEFAULT_ERROR
        try:
            save=_self.get_setting_legacy('auto_zoom')
            if object==None:
                object=name
            model = fragments.get(str(name))
            la = len(model.atom)
            if la:
                mean = map(lambda x,la=la:x/la,[
                    reduce(operator.__add__,map(lambda a:a.coord[0],model.atom)),

                    reduce(operator.__add__,map(lambda a:a.coord[1],model.atom)),
                    reduce(operator.__add__,map(lambda a:a.coord[2],model.atom))])
                position = _self.get_position()
                for c in range(0,3):
                    mean[c]=position[c]-mean[c]
                    map(lambda a,x=mean[c],c=c:_self._adjust_coord(a,c,x),model.atom)
                mean = map(lambda x,la=la:x/la,[
                    reduce(operator.__add__,map(lambda a:a.coord[0],model.atom)),
                    reduce(operator.__add__,map(lambda a:a.coord[1],model.atom)),
                    reduce(operator.__add__,map(lambda a:a.coord[2],model.atom))])
            r = _self.load_model(model,str(object),quiet=quiet,zoom=zoom)
        except IOError:
            print "Error: unable to load fragment '%s'." % name
        except:
            traceback.print_exc()
            print "Error: unable to load fragment '%s'." % name         
        if _self._raising(r,_self): raise pymol.CmdException                                    
        return r

    def create(name, selection, source_state=0,
               target_state=0, discrete=0, zoom=-1, quiet=1,
               singletons=0, extract=None, _self=cmd):
        '''
DESCRIPTION

    "create" creates a new molecule object from a selection.  It can
    also be used to create states in an existing object.

USAGE

    create name, (selection) [,source_state [,target_state ] ]

ARGUMENTS

    name = string: name of object to create or modify

    selection = string: atoms to include in the new object

    source_state = integer: {default: 0 -- copy all states}

    target_state = integer: {default: 0}

PYMOL API

    cmd.create(string name, string selection, int state,
               int target_state, int discrete)

NOTES

    If the source and target states are zero (default), then all
    states will be copied.  Otherwise, only the indicated states will
    be copied.

SEE ALSO

    load, copy, extract
        '''
        r = DEFAULT_ERROR      
        # preprocess selection
        selection = selector.process(selection)
        #      
        try:
            _self.lock(_self)
            if name==None:
                sel_cnt = _cmd.get(_self._COb,"sel_counter") + 1.0
                _cmd.legacy_set(_self._COb,"sel_counter","%1.0f" % sel_cnt)
                name = "obj%02.0f" % sel_cnt
            r = _cmd.create(_self._COb,str(name),"("+str(selection)+")",
                            int(source_state)-1,int(target_state)-1,
                            int(discrete),int(zoom),int(quiet),int(singletons))
        finally:
            _self.unlock(r,_self)
        if not is_error(r): # temporary inefficient implementation
            if extract not in (None, 0, '0'):
                if extract not in (1, '1'):
                    extract = selector.process(extract)
                else:
                    extract = selection
                _self.remove("(("+extract+") in (%s)) and not (%s)"%(name,name))
        if _self._raising(r,_self): raise pymol.CmdException                                    
        return r

    def extract(*arg,**kw):
        '''
DESCRIPTION

    "extract" is simply a shorthand way calling the "create" command
    with the extract argument activated, so that atoms in the new
    object are removed from the source object.

USAGE

    extract name, selection [, source_state [, target_state ]]

SEE ALSO

    create
    
    '''
        
        kw['extract'] = 1
        return apply(create,arg,kw)

    pseudoatom_mode_dict = {
        "unit" : 0, # radius 0.5
        "extent" : 1,
        "rms" : 2,
#        "ellipse" : 2,  for anisotropic b-factors? 
        }

    pseudoatom_mode_sc =  Shortcut(pseudoatom_mode_dict.keys())

    unquote_re = re.compile(r"'[^']*'|\"[^\"]*\"")
    
    def unquote(s):
        s = str(s)
        if unquote_re.search(s):
            return s[1:-1]
        return s
    
    def pseudoatom(object, selection='', name='PS1', resn='PSD', resi='1', chain='P',
                   segi='PSDO', elem='PS', vdw=-1.0, hetatm=1, b=0.0, q=0.0, color='',
                   label='', pos=None, state=0, mode='rms', quiet=1,_self=cmd):
        '''
        
DESCRIPTION

    "pseudoatom" adds a pseudoatom to a molecular object, and will
    creating the molecular object if it does not yet exist.
    
USAGE

    pseudoatom object [, selection [, name [, resn [, resi [, chain
        [, segi [, elem [, vdw [, hetatm [, b [, q [, color [, label
        [, pos [, state [, mode [, quiet ]]]]]]]]]]]]]]]]]

NOTES

    "pseudoatom" can be used for a wide variety of random tasks where
    on must place an atom or a label in 3D space.
    
    '''
        
        r = DEFAULT_ERROR      
        # preprocess selection
        if len(color):
            color = _self.get_color_index(str(color))
        else:
            color = -1 # default
        selection = selector.process(selection)
        mode = pseudoatom_mode_dict[pseudoatom_mode_sc.auto_err(str(mode),'pseudoatom mode')]
        
        (name,resn,resi,chain,segi,elem,label) = map(unquote,(name,resn,resi,chain,segi,elem,label))
        #      
        try:
            _self.lock(_self)
            if pos!=None:
                if not (is_list(pos) or is_tuple(pos)):
                    pos = safe_list_eval(pos)
                pos = (float(pos[0]), # tuple-ize
                       float(pos[1]),
                       float(pos[2]))
            if len(selection.split())>1:
                selection = "("+str(selection)+")"
            r = _cmd.pseudoatom(_self._COb,str(object), str(selection),
                                str(name), str(resn), str(resi), str(chain),
                                str(segi), str(elem), float(vdw), int(hetatm),
                                float(b), float(q), str(label), pos, int(color),
                                int(state)-1, int(mode), int(quiet))
        finally:
            _self.unlock(r,_self)
        if _self._raising(r,_self): raise pymol.CmdException                                    
        return r
        

