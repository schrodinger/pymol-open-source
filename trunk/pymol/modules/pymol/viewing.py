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

if __name__=='pymol.viewing':

    import thread
    import threading
    import string
    import types
    import traceback
    import pymol
    import selector
    import copy
    import parsing
    import re
    import cmd
    
    from cmd import _cmd,lock,unlock,Shortcut,QuietException,_raising, \
          _feedback,fb_module,fb_mask, \
          repres,repres_sc, is_string, is_list, is_ok, is_error, \
          toggle_dict,toggle_sc,stereo_dict,stereo_sc, \
          palette_dict ,palette_sc, window_dict, window_sc, \
          safe_list_eval,  safe_alpha_list_eval, \
          location_code, location_sc, boolean_dict, boolean_sc, \
          DEFAULT_ERROR, DEFAULT_SUCCESS
        
    import thread
    
    rep_list = [ "lines", "sticks", "spheres", "dots", "surface",
                 "mesh", "nonbonded", "nb_spheres", "cartoon",
                 "ribbon", "labels", "slice", "ellipsoids"]

    scene_action_sc = Shortcut(['store','recall','clear','insert_before',
                                'insert_after','next','previous',
                                'start', 'update','rename','delete',
                                'order', 'sort', 'first',
                                'append'])
    scene_action_dict = {}
    scene_action_dict_sc = Shortcut([])

    view_sc = Shortcut(['store','recall','clear'])

    def zoom(selection="all", buffer=0.0, state=0, complete=0, animate=0,_self=cmd):
        '''
DESCRIPTION

    "zoom" scales and translates the window and the origin to cover the
    atom selection.

USAGE

    zoom [ selection [, buffer [, state [, complete [, animate ]]]]]
    
EXAMPLES

    zoom 
    zoom complete=1
    zoom 142/, animate=3
    zoom (chain A)

ARGUMENTS

    selection = string: selection-expression or name pattern {default: all}

    buffer = float: distance  {default: 0}
    
    state = 0: uses all coordinate states {default}
    
    state = -1: uses only coordinates for the current state
    
    state > 0: uses coordinates for a specific state

    complete = 0 or 1: will insure no atoms centers are clipped

    animate < 0: uses the default animation duration
    
    animate = 0: no animation
    
    animate > 0: animates using the provided duration in seconds
    
PYMOL API

    cmd.zoom(string selection, float buffer, int state, int complete,
             int animate)

NOTES

    The zoom command normally tries to guess an optimal zoom level for
    visualization, balancing closeness against occasional clipping of
    atoms out of the field of view.  You can change this behavior by
    setting the complete option to 1, which will guarantee that the
    atom positions for the entire selection will fit in the field of
    an orthoscopic view.

    To absolutely prevent clipping, you may also need to add an
    additional buffer (typically 2 A) to account for graphical
    representations which extend beyond the atom coordinates.  

SEE ALSO

    origin, orient, center
        '''
        # preprocess selection
        selection = selector.process(selection)
        #
        r = DEFAULT_ERROR
        try:
            _self.lock(_self)   
            r = _cmd.zoom(_self._COb,str(selection),float(buffer),
                              int(state)-1,int(complete),int(animate))
        finally:
            _self.unlock(r,_self)
        if _self._raising(r,_self): raise QuietException
        return r

    def center(selection="all", state=0, origin=1, animate=0, _self=cmd):
        '''
DESCRIPTION

    "center" translates the window, the clipping slab, and the
    origin to a point centered within the atom selection.

USAGE

    center [ selection [, state [, origin [, animate ]]]]

EXAMPLES

    chain chain B
    center 145/

ARGUMENTS

    selection = string: selection-expression or name pattern (default: "all").
    
    state = 0 (default) use all coordinate states
    
    state = -1 use only coordinates for the current state
    
    state > 0  use coordinates for a specific state

    origin = 1 (default) move the origin
    
    origin = 0 leave the origin unchanged

PYMOL API

    cmd.center(string selection, int state, int origin)

SEE ALSO

    origin, orient, zoom
        '''
        r = DEFAULT_ERROR
        # preprocess selection
        selection = selector.process(selection)
        #
        try:
            _self.lock(_self)   
            r = _cmd.center(_self._COb,str(selection),int(state)-1,int(origin),int(animate))
        finally:
            _self.unlock(r,_self)
        if _self._raising(r,_self): raise QuietException
        return r

    clip_action_sc = Shortcut([ 'near','far','move','slab','atoms' ])

    def clip(mode, distance, selection=None, state=0, _self=cmd):
        '''
DESCRIPTION

    "clip" alters the positions of the clipping planes.

USAGE

    clip mode, distance [, selection [, state ]]

ARGUMENTS 

    mode = near, far, move, slab, or atoms

    distance is a floating point value

    selection = atom selection (for mode=atoms only)

EXAMPLES

    clip near, -5           # moves near plane away from you by 5 A
    clip far, 10            # moves far plane towards you by 10 A
    clip move, -5           # moves the slab away from you by 5 A
    clip slab, 20           # sets slab thickness to 20 A
    clip slab, 10, resi 11  # clip 10 A slab about residue 11

    clip atoms, 5, pept     # clip atoms in "pept" with a 5 A buffer
                            # about their current camera positions

PYMOL API

    cmd.clip(string mode, float distance, string selection, int state)

SEE ALSO

    zoom, orient, reset
        '''
        r = DEFAULT_ERROR      
        mode = clip_action_sc.auto_err(str(mode),'mode')
        if selection!=None:
            selection = selector.process(selection)
        else:
            selection = ''
        try:
            _self.lock(_self)   
            r = _cmd.clip(_self._COb,str(mode),float(distance),
                          str(selection),int(state)-1)
        finally:
            _self.unlock(r,_self)
        if _self._raising(r,_self): raise QuietException         
        return r

    def origin(selection="(all)", object=None, position=None, state=0, _self=cmd):
        '''
DESCRIPTION

    "origin" sets the center of rotation about a selection.  If an
    object name is specified, it can be used to set the center of
    rotation for the object (for use in animation and editing).

USAGE

    origin [ selection [, object [,position, [, state ]]]]

ARGUMENTS

    selection = string: selection-expression or name-list {default: (all)}
    
    state = 0 (default) use all coordinate states
    
    state = -1 use only coordinates for the current state
    
    state > 0  use coordinates for a specific state

EXAMPLES

    origin chain A
    
    origin position=[1.0,2.0,3.0]

PYMOL API

    cmd.origin(string object-or-selection)

SEE ALSO

    zoom, orient, reset
        '''
        #'
        r = DEFAULT_ERROR      
        # preprocess selection
        selection = selector.process(selection)
        #   
        try:
            _self.lock(_self)
            if object==None: object=''
            if position==None: position=(0.0,0.0,0.0)
            else:
                if _self.is_string(position):
                    position = safe_list_eval(position)
                selection = ''
            r = _cmd.origin(_self._COb,selection,str(object),
                                 (float(position[0]),
                                  float(position[1]),
                                  float(position[2])
                                  ),int(state)-1)
        finally:
            _self.unlock(r,_self)
        if _self._raising(r,_self): raise QuietException         
        return r

    def orient(selection="(all)", state=0, animate=0, _self=cmd):
        '''
DESCRIPTION

    "orient" aligns the principal components of the atoms in the
    selection with the XYZ axes.  

USAGE

    orient [ selection [, state [, animate ]]]

ARGUMENTS

    selection = a selection-expression or name-pattern {default: (all)}

    state = 0: use all coordinate states {default}
    
    state = -1: uses only coordinates for the current state
    
    state > 0: uses coordinates for a specific state

EXAMPLES

    orient organic

NOTES

    The function is similar to the orient command in X-PLOR.

PYMOL API

    cmd.orient(string object-or-selection, int state, float animate)

SEE ALSO

    zoom, origin, reset
        '''
        r = DEFAULT_ERROR      
        # preprocess selection
        selection = selector.process(selection)
        #   
        try:
            _self.lock(_self)
            r = _cmd.orient(_self._COb,"("+selection+")",int(state)-1,float(animate))
        finally:
            _self.unlock(r,_self)
        if _self._raising(r,_self): raise QuietException         
        return r

    def move(axis, distance, _self=cmd):
        '''
DESCRIPTION

    "move" translates the camera about one of the three primary axes.

USAGE

    move axis, distance

EXAMPLES

    move x, 3
    move y, -1

PYMOL API

    cmd.move(string axis, float distance)

SEE ALSO

    turn, rotate, translate, zoom, center, clip
        '''
        r = DEFAULT_ERROR      
        try:
            _self.lock(_self)   
            r = _cmd.move(_self._COb,str(axis),float(distance))
        finally:
            _self.unlock(r,_self)
        if _self._raising(r,_self): raise QuietException         
        return r

    def enable(name='all', parents=0, _self=cmd):
        '''
DESCRIPTION

    "enable" turns on display of one or more objects and/or selections.

USAGE

    enable name

ARGUMENTS    

    name = name-pattern or selection. 

NOTES

    If name matches a selection name, then selection indicator dots
    are shown for atoms in that selection.  If name is a
    selection-expression, then all objects with atoms in that
    selection are enabled.

    For an object\'s content to be displayed in the 3D viewer, the
    object must be enabled AND at least one of the available
    representations must be shown.
    
PYMOL API

    cmd.enable(string object-name)

EXAMPLES

    enable target_protein  # enables the target_protein object

    enable 1dn2.*   # enables all entities starting with 1dn2.
    
    enable *lig     # enables all entities ending with lig
    
SEE ALSO

    show, hide, disable
        '''
        r = DEFAULT_ERROR      
        if name[0]=='(':
            selection = selector.process(name)
            try:
                _self.lock(_self)
                r = _cmd.onoff_by_sele(_self._COb,selection,1)
            finally:
                _self.unlock(r,_self)
        else:
            try:
                _self.lock(_self)   
                r = _cmd.onoff(_self._COb,str(name),1,int(parents));
            finally:
                _self.unlock(r,_self)
        if _self._raising(r,_self): raise QuietException            
        return r

    def disable(name='all',_self=cmd):
        '''
DESCRIPTION

    "disable" turns off display of one or more objects and/or selections.

USAGE

    disable name

ARGUMENTS    

    name = name-pattern or selection.

PYMOL API

    cmd.disable(string name) 

SEE ALSO

    show, hide, enable   
        '''
        r = DEFAULT_ERROR      
        if name[0]=='(':
            selection = selector.process(name)
            try:
                _self.lock(_self)
                r = _cmd.onoff_by_sele(_self._COb,selection,0)
            finally:
                _self.unlock(r,_self)
        else:
            try:
                _self.lock(_self)   
                r = _cmd.onoff(_self._COb,str(name),0,0);
            finally:
                _self.unlock(r,_self)
        if _self._raising(r,_self): raise QuietException            
        return r

    def toggle(representation="", selection="", _self=cmd):
        '''
DESCRIPTION

    "toggle" toggles the visibility of a representation within a
    selection.
    
USAGE

    toggle representation, selection

ARGUMENTS

    representation = string

    selection = string

NOTES

    If no arguments are provided, then lines are toggled for all
    objects in the aggregate.
    
PYMOL API

    cmd.toggle(string representation, string selection)

        '''
        r = DEFAULT_ERROR
        try:
            _self.lock(_self)
            if (representation=="") and (selection==""):
                r = _cmd.toggle(_self._COb,"(all)",repres['lines']); # show lines by default       
            elif representation=='object':
                if selection in cmd.get_names('all',enabled_only=1):
                    r = cmd.disable(selection)
                else:
                    r = cmd.enable(selection)
            elif (representation!="") and (selection!=""):
                rep = representation
                rep = repres_sc.auto_err(rep,'representation')
                repn = repres[rep];
                # preprocess selection 
                selection = selector.process(selection)
                #   
                r = _cmd.toggle(_self._COb,str(selection),int(repn));
            elif representation=='all':
                r = _cmd.toggle(_self._COb,"all",repres['lines']); # toggle lines by default 
            elif (representation[0:1]=='(') or (string.find(representation,'/')>=0):
                # preprocess selection
                selection = selector.process(representation)
                #                  
                r = _cmd.toggle(_self._COb,str(selection),repres['lines']);
            else: # selection==""
                rep = representation
                rep = repres_sc.auto_err(rep,'representation')
                repn = repres[rep];
                r = _cmd.toggle(_self._COb,"all",int(repn));
        finally:
            _self.unlock(r,_self)
        if _self._raising(r,_self): raise QuietException
        return r

    def show(representation="", selection="", _self=cmd):
        '''
DESCRIPTION

    "show" turns on representations for objects and selections.

USAGE

    show [ representation [, selection ]]

ARGUMENTS

    representation = lines, spheres, mesh, ribbon, cartoon, sticks,
       dots, surface, labels, extent, nonbonded, nb_spheres, slice,
       extent, slice, dashes, angles, dihedrals, cgo, cell, callback,
       or everything

    selection = string: a selection-expression or name-pattern

NOTES

    With no arguments, "show" alone turns on lines for all bonds and
    nonbonded for all atoms in all molecular objects.

EXAMPLES

    show
    show ribbon
    show lines, (name ca or name c or name n)

SEE ALSO

    hide, enable, disable

'''
        r = DEFAULT_ERROR
        try:
            _self.lock(_self)
            if (representation=="") and (selection==""):
                if is_ok(_cmd.showhide(_self._COb,"(all)",repres['lines'],1)): # show lines by default
                    r = _cmd.showhide(_self._COb,"(all)",repres['nonbonded'],2)
            elif (representation!="") and (selection!=""):
                rep = representation
                rep = repres_sc.auto_err(rep,'representation')
                repn = repres[rep];
                # preprocess selection 
                selection = selector.process(selection)
                #   
                r = _cmd.showhide(_self._COb,str(selection),int(repn),1);
            elif representation=='all':
                if is_ok(_cmd.showhide(_self._COb,"all",repres['lines'],1)): # show lines by default
                    r = _cmd.showhide(_self._COb,"all",repres['nonbonded'], 1) # nonbonded
            elif (representation[0:1]=='(') or (string.find(representation,'/')>=0):
                # preprocess selection
                selection = selector.process(representation)
                #                  
                if is_ok(_cmd.showhide(_self._COb,str(selection),repres['lines'],1)):
                    r = _cmd.showhide(_self._COb,str(selection),repres['nonbonded'],2);
            else: # selection==""
                rep = representation
                rep = repres_sc.auto_err(rep,'representation')
                repn = repres[rep]
                r = _cmd.showhide(_self._COb,"all",int(repn),1);
        finally:
            _self.unlock(r,_self)
        if _self._raising(r,_self): raise QuietException         
        return r

    def show_as(representation="", selection="", _self=cmd):
        '''
DESCRIPTION

    "as" turns on and off atom and bond representations.

    The available representations are:

        lines     spheres    mesh      ribbon     cartoon
        sticks    dots       surface   labels     extent
        nonbonded nb_spheres slice

USAGE

    as representation [, selection ]

ARGUMENTS
    
    representation = lines, spheres, mesh, ribbon, cartoon, sticks,
        dots, surface, labels, extent, nonbonded, nb_spheres, slice,
        extent, slice, dashes, angles, dihedrals, cgo, cell, callback,
        or everything

    selection = string {default: all}

EXAMPLES

    as lines, name ca or name c or name n

    as ribbon

PYMOL API

    cmd.show_as(string representation, string selection)

NOTES

    "selection" can be an object name
    "as" alone will turn on lines and nonbonded and hide everything else.

SEE ALSO

    show, hide, enable, disable
        '''
        r = DEFAULT_ERROR
        try:
            _self.lock(_self)
            if (representation=="") and (selection==""):
                if is_ok(_cmd.showhide(_self._COb,str(selection),-1,0)):
                    if is_ok(_cmd.showhide(_self._COb,"(all)",repres['lines'],1)):
                        r = _cmd.showhide(_self._COb,"(all)",repres['nonbonded'],1)
            elif (representation!="") and (selection!=""):
                rep = representation
                rep = repres_sc.auto_err(rep,'representation')
                repn = repres[rep]
                # preprocess selection 
                selection = selector.process(selection)
                #
                if is_ok(_cmd.showhide(_self._COb,str(selection),-1,0)):
                    r = _cmd.showhide(_self._COb,str(selection),int(repn),1)
            elif representation=='all':
                if is_ok(_cmd.showhide(_self._COb,str(selection),-1,0)):            
                    if if_ok(_cmd.showhide(_self._COb,"all",repres['lines'],1)): # show lines by default
                        r = _cmd.showhide(_self._COb,"all",repres['nonbonded'],1) # show nonbonded by default
            elif (representation[0:1]=='(') or (string.find(representation,'/')>=0):
                # preprocess selection
                selection = selector.process(representation)
                if is_ok(_cmd.showhide(_self._COb,str(selection),-1,0)):
                    r = _cmd.showhide(_self._COb,str(selection),repres['lines'],1)
            else: # selection==""
                rep = representation
                rep = repres_sc.auto_err(rep,'representation')
                repn = repres[rep];
                if is_ok(_cmd.showhide(_self._COb,"all",-1,0)):
                    r = _cmd.showhide(_self._COb,"all",int(repn),1);
        finally:
            _self.unlock(r,_self)
        if _self._raising(r,_self): raise QuietException
        return r

    def hide(representation="", selection="",_self=cmd):
        '''
DESCRIPTION

    "hide" turns of atom and bond representations.


USAGE

    hide [ representation [, selection ]]

ARGUMENTS

    representation = lines, spheres, mesh, ribbon, cartoon,
       sticks, dots, surface, labels, extent, nonbonded, nb_spheres,
       slice, extent, slice, dashes, angles, dihedrals, cgo, cell, callback, 
       or everything

    selection = string: a selection-expression or name-pattern

EXAMPLES

    hide lines, all
    hide ribbon

PYMOL API

    cmd.hide(string representation, string selection)

SEE ALSO

    show, enable, disable

        '''
        r = DEFAULT_ERROR
        try:
            _self.lock(_self)
            if (representation=="") and (selection==""):
                r = _cmd.showhide(_self._COb,"@",0,0);      
            elif (representation!="") and (selection!=""):
                rep = representation
                rep = repres_sc.auto_err(rep,'representation')
                repn = repres[rep];
                selection = selector.process(selection)
                r = _cmd.showhide(_self._COb,str(selection),int(repn),0);
            elif (representation=='all'):
                r = _cmd.showhide(_self._COb,"@",0,0);
            elif (representation[0:1]=='(') or (string.find(representation,'/')>=0):
                selection = selector.process(representation)
                r = _cmd.showhide(_self._COb,str(selection),-1,0);
            else: # selection == ""
                rep = representation
                rep = repres_sc.auto_err(rep,'representation')
                repn = repres[rep];
                r = _cmd.showhide(_self._COb,"all",int(repn),0);
        finally:
            _self.unlock(r,_self)
        if _self._raising(r,_self): raise QuietException         
        return r


    def get_view(output=1, quiet=1, _self=cmd):
        '''
DESCRIPTION

    "get_view" returns and optionally prints out the current view
    information in a format which can be embedded into a command
    script and can be used in subsequent calls to "set_view".

    If a log file is currently open, get_view will not write the view
    matrix to the screen unless the "output" parameter is 2.

USAGE

    get_view [output]

ARGUMENTS

    output = 0: output matrix to screen

    output = 1: do not Output matrix to screen

    output = 2: force output to screen even if log file is open

    output = 3: return formatted string instead of a list

NOTES

    Contents of the view matrix:
    
    * 0  -  8: column-major 3x3 matrix which rotates model space to camera space

    * 9  - 11: origin of rotation relative to camera (in camera space)

    * 12 - 14: origin of rotation (in model space)

    * 15: front plane distance from the camera

    * 16: rear plane distance from the camera

    * 17: orthoscopic flag (+/-) and field of view (if abs(value) > 1)

    The camera always looks down -Z with its +X left and its +Y down.

    Therefore, in the default view, model +X is to the observer\'s
    right, +Y is upward, and +Z points toward the observer.

PYMOL API

    cmd.get_view(output=1, quiet=1)

SEE ALSO

    set_view
    '''

        r = DEFAULT_ERROR
        try:
            _self.lock(_self)
            r = _cmd.get_view(_self._COb)
        finally:
            _self.unlock(r,_self)
        if is_ok(r):
            if len(r):
                if (_self.get_setting_legacy("logging")!=0.0) and (output<3):
                    if not quiet:
                        print " get_view: matrix written to log file."
                    _self.log("_ set_view (\\\n","cmd.set_view((\\\n")
                    _self.log("_  %14.9f, %14.9f, %14.9f,\\\n"%r[0:3]  ,
                              "  %14.9f, %14.9f, %14.9f,\\\n"%r[0:3])
                    _self.log("_  %14.9f, %14.9f, %14.9f,\\\n"%r[4:7]  ,
                              "  %14.9f, %14.9f, %14.9f,\\\n"%r[4:7])
                    _self.log("_  %14.9f, %14.9f, %14.9f,\\\n"%r[8:11] ,
                              "  %14.9f, %14.9f, %14.9f,\\\n"%r[8:11])
                    _self.log("_  %14.9f, %14.9f, %14.9f,\\\n"%r[16:19],
                              "  %14.9f, %14.9f, %14.9f,\\\n"%r[16:19])
                    _self.log("_  %14.9f, %14.9f, %14.9f,\\\n"%r[19:22],
                              "  %14.9f, %14.9f, %14.9f,\\\n"%r[19:22]) 
                    _self.log("_  %14.9f, %14.9f, %14.9f )\n"%r[22:25] ,
                              "  %14.9f, %14.9f, %14.9f ))\n"%r[22:25])
                    if output<2: # suppress if we have a log file open
                        output=0
                if output and (not quiet) and (output<3):
                    print "### cut below here and paste into script ###"
                    print "set_view (\\"
                    print "  %14.9f, %14.9f, %14.9f,\\"%r[0:3]
                    print "  %14.9f, %14.9f, %14.9f,\\"%r[4:7]
                    print "  %14.9f, %14.9f, %14.9f,\\"%r[8:11]
                    print "  %14.9f, %14.9f, %14.9f,\\"%r[16:19]
                    print "  %14.9f, %14.9f, %14.9f,\\"%r[19:22]
                    print "  %14.9f, %14.9f, %14.9f )"%r[22:25]
                    print "### cut above here and paste into script ###"
            if output==3:
                return ("set_view (\\\n"+
                  "  %14.9f, %14.9f, %14.9f,\\\n"%r[0:3] +
                  "  %14.9f, %14.9f, %14.9f,\\\n"%r[4:7] +
                  "  %14.9f, %14.9f, %14.9f,\\\n"%r[8:11] +
                  "  %14.9f, %14.9f, %14.9f,\\\n"%r[16:19] +
                  "  %14.9f, %14.9f, %14.9f,\\\n"%r[19:22] +
                  "  %14.9f, %14.9f, %14.9f )\n"%r[22:25])
            r = r[0:3]+r[4:7]+r[8:11]+r[16:25]
        elif _self._raising(r,_self):
            raise QuietException
        return r

    def set_view(view,animate=0,quiet=1,hand=1,_self=cmd):
        '''
DESCRIPTION

    "set_view" sets viewing information for the current scene,
    including the rotation matrix, position, origin of rotation,
    clipping planes, and the orthoscopic flag.

USAGE

    set_view [ view ] 

EXAMPLE

    set_view (\
        0.999876618,   -0.000452542,   -0.015699286,\
        0.000446742,    0.999999821,   -0.000372844,\
        0.015699454,    0.000365782,    0.999876678,\
        0.000000000,    0.000000000, -150.258514404,\
        11.842411041,   20.648729324,    8.775371552,\
        118.464958191,  182.052062988,    0.000000000 )

PYMOL API

    cmd.set_view(string-or-sequence view)  

SEE ALSO

    get_view
    '''
        r = DEFAULT_ERROR
        if cmd.is_string(view):
            try:
                view = eval(re.sub(r"[^0-9,\-\)\(\.]","",view))
            except:
                traceback.print_exc()
                print "Error: bad view argument; should be a sequence of 18 floats."
                raise QuietException
        if len(view)!=18:
            print "Error: bad view argument; should be a sequence of 18 floats."
            raise QuietException
        else:
            try:
                _self.lock(_self)
                r = _cmd.set_view(_self._COb,(
                    float(view[ 0]),float(view[ 1]),float(view[ 2]),0.0,
                    float(view[ 3]),float(view[ 4]),float(view[ 5]),0.0,
                    float(view[ 6]),float(view[ 7]),float(view[ 8]),0.0,
                    0.0,0.0,0.0,1.0,
                    float(view[ 9]),float(view[10]),float(view[11]),
                    float(view[12]),float(view[13]),float(view[14]),
                    float(view[15]),float(view[16]),float(view[17])),
                    int(quiet),float(animate),int(hand))
            finally:
                _self.unlock(r,_self)
        if _self._raising(r,_self): raise QuietException
        return r

    def view(key, action='recall', animate=-1,_self=cmd):
        '''
DESCRIPTION

    "view" saves and restore camera views.

USAGE

    view key [, action [, animate]]
    
ARGUMENTS

    key = string or *

    action = store or recall: {default: recall}

NOTES

    Views F1 through F12 are automatically bound to function keys
    provided that "set_key" has not been used to redefine the
    behaviour of the respective key, and that a "scene" has not been
    defined for that key.

EXAMPLES

    view 0, store
    view 0

PYMOL API

    cmd.view(string key, string action)

SEE ALSO

    scene, set_view, get_view
        '''
        pymol=_self._pymol
        
        if key=='*':
            action = view_sc.auto_err(action,'action')
            if action=='clear':
                pymol._view_dict = {}
                pymol._view_dict_sc = Shortcut(pymol._view_dict.keys())                        
            else:
                print " view: stored views:"
                lst = pymol._view_dict.keys()
                lst.sort()
                parsing.dump_str_list(lst)
                
        else:
            action = view_sc.auto_err(action,'action')
            if action=='recall':
                key = pymol._view_dict_sc.auto_err(key,'view')
                _self.set_view(pymol._view_dict[key],animate=animate)
                if _feedback(fb_module.scene,fb_mask.actions,_self): # redundant
                    print " view: \"%s\" recalled."%key
            elif (action=='store') or (action=='update'):
                pymol._view_dict_sc.append(key)
                pymol._view_dict[key]=_self.get_view(0)
                if _feedback(fb_module.scene,fb_mask.actions,_self):
                    print " view: view "+action+"d as \"%s\"."%key
            elif action=='clear':
                key = pymol._view_dict_sc.auto_err(key,'view')
                if pymol._view_dict.has_key(key):
                    del pymol._view_dict[key]
                    pymol._view_dict_sc = Shortcut(pymol._view_dict.keys())            
                    if _feedback(fb_module.scene,fb_mask.actions,_self): # redundant
                        print " view: '%s' deleted."%key


    def get_vis(_self=cmd):
        r = DEFAULT_ERROR      
        try:
            _self.lock(_self)
            r = _cmd.get_vis(_self._COb)
        finally:
            _self.unlock(r,_self)
        if _self._raising(r,_self): raise QuietException
        return r

    def set_vis(dict,_self=cmd):
        r = DEFAULT_ERROR      
        try:
            _self.lock(_self)         
            r = _cmd.set_vis(_self._COb,dict)
        finally:
            _self.unlock(r,_self)
        if _self._raising(r,_self): raise QuietException
        return r

    def get_colorection(key,_self=cmd):
        r = DEFAULT_ERROR      
        try:
            _self.lock(_self)         
            r = _cmd.get_colorection(_self._COb,key)
        finally:
            _self.unlock(r,_self)
        if _self._raising(r,_self): raise QuietException
        return r

    def set_colorection(dict,key,_self=cmd):
        r = DEFAULT_ERROR      
        try:
            _self.lock(_self)         
            r = _cmd.set_colorection(_self._COb,dict,key)
        finally:
            _self.unlock(r,_self)
        if _self._raising(r,_self): raise QuietException
        return r

    def set_colorection_name(dict,key,new_key,_self=cmd):
        r = DEFAULT_ERROR      
        try:
            _self.lock(_self)         
            r = _cmd.set_colorection_name(_self._COb,dict,key,new_key)
        finally:
            _self.unlock(r,_self)
        if _self._raising(r,_self): raise QuietException
        return r

    def del_colorection(dict,key,_self=cmd):
        r = DEFAULT_ERROR      
        try:
            _self.lock(_self)         
            r = _cmd.del_colorection(_self._COb,dict,key)
        finally:
            _self.unlock(r,_self)
        if _self._raising(r,_self): raise QuietException
        return r

    def get_scene_dict(_self=cmd):
        pymol=_self._pymol        
        r = DEFAULT_ERROR      
        try:
            _self.lock(_self)
            cpy = copy.deepcopy(pymol._scene_dict)
        finally:
            _self.unlock(r,_self)
        return cpy

    def get_scene_list(_self=cmd):
        return copy.deepcopy(_scene_validate_list(_self))

    scene_sort_dict = {
        'F1' : 'F01',
        'F2' : 'F02',
        'F3' : 'F03',
        'F4' : 'F04',
        'F5' : 'F05',
        'F6' : 'F06',
        'F7' : 'F07',
        'F8' : 'F08',
        'F9' : 'F09',
        'SHFT-F1' : 'SHFT-F01',
        'SHFT-F2' : 'SHFT-F02',
        'SHFT-F3' : 'SHFT-F03',
        'SHFT-F4' : 'SHFT-F04',
        'SHFT-F5' : 'SHFT-F05',
        'SHFT-F6' : 'SHFT-F06',
        'SHFT-F7' : 'SHFT-F07',
        'SHFT-F8' : 'SHFT-F08',
        'SHFT-F9' : 'SHFT-F09',
        }

    def _scene_get_unique_key(_self=cmd):
        pymol=_self._pymol
        keys = pymol._scene_dict.keys()
        while 1:
            key = "%03d"%pymol._scene_counter
            if pymol._scene_dict.has_key(key):
                pymol._scene_counter = pymol._scene_counter + 1
            else:
                break;
        return key
    
    def _scene_validate_list(_self=cmd):
        pymol=_self._pymol
        new_list = []
        new_dict = {}
        
        for a in pymol._scene_order:
            if pymol._scene_dict.has_key(a) and not new_dict.has_key(a):
                new_list.append(a)
                new_dict[a] = 1

        lst = map(lambda x:(scene_sort_dict.get(x,x),x), pymol._scene_dict.keys())
        lst.sort()
        lst = map(lambda x:x[1],lst)
        for a in lst:
            if not new_dict.has_key(a):
                new_list.append(a)
        pymol._scene_order = new_list
        # update shortcuts
        pymol._scene_dict_sc.rebuild( pymol._scene_dict.keys() )
        # update PyMOL internally
        r = DEFAULT_ERROR
        try:
            _self.lock(_self)
            r = _cmd._set_scene_names(_self._COb,pymol._scene_order)
        finally:
            _self.unlock(r,_self);
        return pymol._scene_order

    def chain_session(_self=cmd):
        import os
        # assumes locked interpreter
        r = 0
        session_file = str(_self.get("session_file"))
        re_pat = re.compile("[0-9]+\.")
        if len(session_file): # find next session file, if it exists
            mo = re_pat.search(session_file)
            if mo!=None:
                pat = mo.group(0)
                if len(pat):
                    file_no = int(float(pat)) + 1
                    new_form = r"%0"+str(len(pat)-1)+"d."
                    for new_num in range(file_no, file_no+11):
                        new_pat = new_form % new_num
                        new_file = re_pat.sub(new_pat, session_file)
                        # try both PSE and PSW
                        if not os.path.exists(new_file):
                            new_file = re.sub("\.pse$",".psw",new_file,re.I)
                        if not os.path.exists(new_file):
                            new_file = re.sub("\.psw$",".pse",new_file,re.I)
                        if os.path.exists(new_file):
                            _self.do("_ cmd.load(r'''"+new_file+"''',format='psw')")
                            return 1
        return 0

    def scene_order(names,sort=0,location='current',quiet=1,_self=cmd):
        '''
DESCRIPTION

    "scene_order" changes the ordering of scenes.

USAGE

    scene_order names, sort, location

ARGUMENTS

    names = string: a space-separated list of names

    sort = yes or no {default: no}

    location = top, current, or bottom {default: current}

EXAMPLES

    scene_order *,yes             
    scene_order F6 F4 F3 
    scene_order 003 006 004, location=top

PYMOL API

    cmd.scene_order(string names, string sort, string location)

SEE ALSO

    scene
    '''
        pymol=_self._pymol

        r = DEFAULT_ERROR
        location=location_code[location_sc.auto_err(location,'location')]
        if is_string(sort):
            sort=boolean_dict[boolean_sc.auto_err(sort,'sort option')]
        if is_string(names):
            if names=='*':
                names = _scene_validate_list(_self)
            else:
                names = names.split()
        r = DEFAULT_SUCCESS
        if len(names):
            name_dict = {}
            for name in names:
                if not pymol._scene_dict.has_key(name):
                    print " Error: scene '%s' is not defined."%name
                    r = DEFAULT_ERROR
                name_dict[name] = 1
            if sort: # special F-key-name aware sort
                sort_list = [] 
                fkey_re = re.compile("[fF]([1-9][0-9]*).*")
                for name in names:
                    mo = fkey_re.match(name)
                    if mo != None:
                        sort_list.append( ("F%02d\n"%int(mo.group(1)),name) )
                    else:
                        sort_list.append( (name,name) )
                sort_list.sort()
                names = []
                for name in sort_list:
                    names.append(name[1])
            if not is_error(r):
                old_list = _scene_validate_list(_self)
                if location < 0: # top
                    new_list = names
                    for name in old_list:
                        if not name_dict.has_key(name):
                            new_list.append(name)
                    pymol._scene_order = new_list
                elif location == 0: # current
                    start = old_list.index(names[0])
                    if start >= 0: # guaranteed?
                        new_list = []
                        for name in old_list[:start]:
                            if not name_dict.has_key(name):
                                new_list.append(name)
                        new_list.extend(names)
                        for name in old_list[start:]:
                            if not name_dict.has_key(name):
                                new_list.append(name)
                        pymol._scene_order = new_list
                else: # bottom
                    new_list = []
                    for name in old_list:
                        if not name_dict.has_key(name):
                            new_list.append(name)
                    new_list.extend(names)
                    pymol._scene_order = new_list
        if _self._raising(r,_self):
            raise pymol.CmdException
        _scene_validate_list(_self) # force scene buttons to be updated
        if not quiet and not is_error(r):
            _self.scene('*')
        return r
        
    def scene(key='auto', action='recall', message=None, view=1,
              color=1, active=1, rep=1, frame=1, animate=-1,
              new_key=None, hand=1, quiet=1, _self=cmd):

        '''
DESCRIPTION

    "scene" saves and restores scenes.  A scene consists of the camera
    view, all object activity information, all atom-wise visibilities,
    all atom-wise colors, all representations, the global frame index,
    and may contain a text message to display on playback.

USAGE

    scene [key [,action [, message, [ new_key=new-key-value ]]]]

ARGUMENTS

    key = string, new, auto, or *: use new for an automatically
    numbered new scene, use auto for the current scene (if one
    exists), and use * for all scenes (clear and recall actions only).
    
    action = store, recall, insert_after, insert_before, next,
    previous, update, rename, or clear: (default = recall).  If
    rename, then a new_key argument must be explicitly defined.

    message = string: a text message to display with the scene.

    new_key = string: the new name for the scene
    
EXAMPLES

    scene *

    scene F1, store
    scene F2, store, Please note the critical hydrogen bond shown in yellow.

    scene F1
    scene F2

    scene F1, rename, new_key=F5

NOTES

    Scenes F1 through F12 are automatically bound to function keys
    provided that "set_key" has not been used to redefine the behaviour
    of the respective key.
    
SEE ALSO

    view, set_view, get_view

        '''
        pymol=_self._pymol
        r = DEFAULT_SUCCESS
        view = int(view)
        rep = int(rep)
        color = int(color)
        active = int(active)
        frame = int(frame)
        quiet = int(quiet)

        if not view:
            animate=0
        elif animate<0:
            scene_animation = int(_self.get_setting_legacy("scene_animation"))
            if scene_animation<0:
                scene_animation = int(_self.get_setting_legacy("animation"))
            if scene_animation!=0:
                animate = _self.get_setting_legacy("scene_animation_duration")
            else:
                animate = 0
        try:
            _self.lock(_self) # manipulating global data, so need lock
            if key=='auto':
                action = scene_action_sc.auto_err(action,'action')
                if action=='recall':
                    action='next'
                if action=='start':
                    lst = _scene_validate_list(_self)
                    if len(lst):
                        key = lst[0]
                        action='recall'
            if action == 'delete':
                action='clear'             
            if action == 'append':
                action='store'
            if key=='*':
                action = scene_action_sc.auto_err(action,'action')
                if action=='clear':
                    for key in pymol._scene_dict.keys():
                        # free selections
                        list = pymol._scene_dict[key]
                        if len(list)>3:
                            colorection = list[3]
                            if colorection!=None:
                                _self.del_colorection(colorection,key) 
                        name = "_scene_"+key+"_*"
                        _self.delete(name)
                    pymol._scene_dict = {}
                    pymol._scene_order = []
                    _scene_validate_list(_self)
                elif action in ['sort']:
                    _self.scene_order(key,1,quiet=quiet)                    
                else:
                    print " scene: current order:"
                    lst = _scene_validate_list(_self)
                    parsing.dump_str_list(lst)
            else:
                action = scene_action_sc.auto_err(action,'action')
                
                if action=='insert_before':
                    key = _scene_get_unique_key(_self=_self)
                    cur_scene = setting.get("scene_current_name",_self=_self)
                    _scene_validate_list(_self)
                    if cur_scene in pymol._scene_order:
                        ix = pymol._scene_order.index(cur_scene)
                    else:
                        ix = 0
                    pymol._scene_order.insert(ix,key)
                    action='store'
                elif action=='rename':
                    if not pymol._scene_dict.has_key(key):
                        print "Error: scene '%s' not found."%key
                    elif new_key==None:
                        print "Error: must provide the 'new_key' argument"
                    else:
                        new_scene_order = []
                        for a in pymol._scene_order:
                            if a==key:
                                new_scene_order.append(new_key)
                            else:
                                new_scene_order.append(a)
                        pymol._scene_order=new_scene_order
                        pymol._scene_dict[new_key] = pymol._scene_dict[key]
                        del pymol._scene_dict[key]
                        valid_names = _self.get_names("all")
                        for rep_name in rep_list:
                            name = "_scene_"+key+"_"+rep_name
                            if name in valid_names:
                                new_name = "_scene_"+new_key+"_"+rep_name
                                _self.set_name(name,new_name)
                        list = pymol._scene_dict[new_key]
                        if len(list)>3:
                            _self.set_colorection_name(list[3],key,new_key)
                        print" scene: '%s' renamed to '%s'."%(key,new_key)
                        pymol._scene_dict_sc.rebuild( pymol._scene_dict.keys())
                        _self.set("session_changed",1,quiet=1)
                        _scene_validate_list(_self) # force update of scene buttons                        
                elif action=='insert_after':
                    key = _scene_get_unique_key(_self=_self)            
                    cur_scene = setting.get("scene_current_name",_self=_self)
                    _scene_validate_list(_self)                    
                    if cur_scene in pymol._scene_order:
                        ix = pymol._scene_order.index(cur_scene) + 1
                    else:
                        ix = len(pymol._scene_order)
                    pymol._scene_order.insert(ix,key)
                    action='store'
                if action=='recall':
                    _self.set("scenes_changed",1,quiet=1);
                    key = pymol._scene_dict_sc.auto_err(key,'scene')
                    _self.set('scene_current_name', key, quiet=1)               
                    list = pymol._scene_dict[key]
                    ll = len(list)
                    if (ll>1) and (active):
                        if list[1]!=None:
                            _self.disable()
                            _self.deselect()
                            _self.set_vis(list[1])
                    if (ll>2) and (frame):
                        if list[2]!=None:
                            if not _self.get_movie_playing(): # don't set frame when movie is already playing
                                if _self.get_frame()!=list[2]: # only set the frame when it isn't already correct
                                    _self.frame(list[2])
                    if (ll>3) and (color):
                        if list[3]!=None:
                            _self.set_colorection(list[3],key)
                    if (ll>4) and (rep):
                        if list[4]==None:
                            rep = 0
                    if (ll>5) and (message==None):
                        if list[5]!=None:
                            message=list[5]
                    if rep!=0:
                        _self.hide("(all)")
                        valid_names = _self.get_names("all")
                        for rep_name in rep_list:
                            name = "_scene_"+key+"_"+rep_name
                            if name in valid_names:
                                _self.show(rep_name,name)
                    replace_flag = 0
                    wiz = _self.get_wizard()
                    if wiz!=None:
                        if str(wiz.__class__) == 'pymol.wizard.message.Message':
                            if hasattr(wiz,'from_scene'):
                                replace_flag = 1
                    mess_flag = 0
                    if message!=None:
                        if is_string(message):
                            if len(message):
                                if(replace_flag):
                                    _self.replace_wizard("message",message)
                                else:
                                    _self.wizard("message",message)
                                _self.get_wizard().from_scene = 1
                                mess_flag = 1
                        if is_list(message):
                            if len(message):
                                if(replace_flag):
                                    apply(_self.replace_wizard,("message",)+tuple(message))
                                else:
                                    apply(_self.wizard,("message",)+tuple(message))
                                _self.get_wizard().from_scene = 1
                                mess_flag = 1
                    if replace_flag and not mess_flag:
                        _self.wizard()
                    if (ll>0) and (view):
                        if list[0]!=None:
                            _self.set_view(list[0],animate,quiet,hand)
                    if not quiet and _feedback(fb_module.scene,fb_mask.actions,_self): # redundant
                        print " scene: \"%s\" recalled."%key
                elif (action=='store') or (action=='update'):
                    if key =='new':
                        key=_scene_get_unique_key(_self=_self)               
                    if key =='auto':
                        key = setting.get("scene_current_name",_self=_self)
                        if key=='':
                            key=_scene_get_unique_key(_self=_self)
                    if not pymol._scene_dict.has_key(key):
                        pymol._scene_dict_sc.append(key)
                    else: # get rid of existing one (if exists)
                        list = pymol._scene_dict[key]
                        if (action=='update') and (message==None) and len(list)>5:
                            message = pymol._scene_dict[key][5]
                        if len(list)>3:
                            colorection = list[3]
                            if colorection!=None:
                                _self.del_colorection(colorection,key) # important -- free RAM
                        name = "_scene_"+key+"_*"
                        _self.delete(name)
                    if key not in pymol._scene_order:
                        pymol._scene_order.append(key)
                    entry = []
                    if view:
                        entry.append(_self.get_view(0))
                    else:
                        entry.append(None);
                    if active:
                        entry.append(_self.get_vis())
                    else:
                        entry.append(None)
                    if frame:
                        entry.append(_self.get_frame())
                    else:
                        entry.append(None)
                    if color:
                        entry.append(_self.get_colorection(key))
                    else:
                        entry.append(None)
                    if rep:
                        entry.append(1)
                        for rep_name in rep_list:
                            name = "_scene_"+key+"_"+rep_name
                            _self.select(name,"rep "+rep_name)
                    else:
                        entry.append(None)
                    if is_string(message):
                        if len(message)>1:
                            if (message[0:1] in [ '"',"'"] and
                                 message[-1:] in [ '"',"'"]):
                                message=message[1:-1]
                            else:
                                message = string.split(message,"\n")
                    entry.append(message)
                    pymol._scene_dict[key]=entry
                    if _feedback(fb_module.scene,fb_mask.actions,_self):
                        print " scene: scene "+action+"d as \"%s\"."%key
                    _scene_validate_list(_self)                        
                    _self.set("scenes_changed",1,quiet=1);
                    _self.set('scene_current_name',key,quiet=1)
                    _self.set("session_changed",1,quiet=1)
                elif action=='clear':
                    if key=='auto':
                        key = setting.get("scene_current_name",_self=_self)
                    key = pymol._scene_dict_sc.auto_err(key,'scene')
                    if pymol._scene_dict.has_key(key):
                        list = pymol._scene_dict[key]
                        if len(list)>3:
                            colorection = list[3]
                            if colorection!=None:
                                _self.del_colorection(colorection,key) # important -- free RAM
                        lst = _scene_validate_list(_self)
                        if key == setting.get("scene_current_name",_self=_self):
                            ix = lst.index(key) - 1
                            if ix>=0:
                                setting.set("scene_current_name",lst[ix],quiet=1,_self=_self)
                            else:
                                setting.set("scene_current_name","",quiet=1,_self=_self)                                
                        _self.set("scenes_changed",1,quiet=1);               
                        del pymol._scene_dict[key]
                        name = "_scene_"+key+"_*"
                        _self.delete(name)
                        _scene_validate_list(_self)
                        if _feedback(fb_module.scene,fb_mask.actions,_self):
                            print " scene: '%s' deleted."%key
                    _self.set("session_changed",1,quiet=1)                                                                    
                elif action=='next':
                    lst = _scene_validate_list(_self)
                    cur_scene = setting.get('scene_current_name',quiet=1,_self=_self)
                    if cur_scene in lst:
                        ix = lst.index(cur_scene) + 1
                    else:
                        ix = 0
                        animate = 0
                        if ((pymol._scene_quit_on_action==action) and
                             (setting.get("presentation",_self=_self)=="on") and 
                             (setting.get("presentation_auto_quit",_self=_self)=="on")):
                            _self.quit()
                    if ix<len(lst):
                        scene_name = lst[ix]
                        _self.set('scene_current_name', scene_name, quiet=1)
                        scene(scene_name,'recall',animate=animate,_self=_self)
                    elif setting.get("scene_loop",_self=_self)=="on": # loop back to the beginning
                        if len(lst):
                            scene_name = lst[0]
                            _self.set('scene_current_name', scene_name, quiet=1)
                            scene(scene_name,'recall',animate=animate,_self=_self)
                    else: # otherwise put up blank screen
                        _self.set('scene_current_name','',quiet=1)
                        chained = 0
                        if (setting.get("presentation",_self=_self)=="on"):
                            chained = chain_session(_self)
                            if (not chained) and (setting.get("presentation_auto_quit",_self=_self)=="on"):
                                pymol._scene_quit_on_action = action
                        if not chained: # and len(lst):
                            _self.disable() # just hide everything
                            _self.wizard()
                            
                elif action=='previous':
                    lst = _scene_validate_list(_self)            
                    cur_scene = setting.get('scene_current_name',quiet=1,_self=_self)
                    if cur_scene in lst:
                        ix = lst.index(cur_scene) - 1
                    else:
                        ix = len(lst)-1
                        animate = 0
                        if ((pymol._scene_quit_on_action==action) and
                             (setting.get("presentation",_self=_self)=="on") and 
                             (setting.get("presentation_auto_quit",_self=_self)=="on")):
                            _self.quit()
                    if ix>=0:
                        scene_name = lst[ix]
                        scene(scene_name,'recall',animate=animate,hand=-1,_self=_self)
                    elif setting.get("scene_loop",_self=_self)=="on": # loop back to the end
                        print setting.get("scene_loop",_self=_self)
                        if len(lst):
                            scene_name = lst[-1]
                            _self.set('scene_current_name', scene_name, quiet=1)
                            scene(scene_name,'recall',animate=animate,hand=-1,_self=_self)
                    else: # otherwise put up blank screen
                        _self.set('scene_current_name','',quiet=1)
                        if ((setting.get("presentation",_self=_self)=="on") and 
                             (setting.get("presentation_auto_quit",_self=_self)=="on")):
                            pymol._scene_quit_on_action = action
                        if len(lst):
                            _self.disable() # just hide everything
                            _self.wizard()
                elif action=='order':
                    _self.scene_order(key,quiet=quiet)
                elif action=='sort':
                    _self.scene_order(key,1,quiet=quiet)
                elif action=='first':
                    _self.scene_order(key,location='top',quiet=quiet)
                
        finally:
            _self.unlock(r,_self)
        return r
                        
    def session_save_views(session,_self=cmd):
        pymol=_self._pymol        
        session['view_dict']=copy.deepcopy(pymol._view_dict)
        return 1

    def session_restore_views(session,_self=cmd):
        pymol=_self._pymol
        if session.has_key('view_dict'):
            pymol._view_dict=copy.deepcopy(session['view_dict'])
            pymol._view_dict_sc.rebuild(pymol._view_dict.keys())
        return 1


    def session_save_scenes(session,_self=cmd):
        pymol=_self._pymol        
        session['scene_dict']=copy.deepcopy(pymol._scene_dict)
        session['scene_order']=copy.deepcopy(pymol._scene_order)
        return 1

    def session_restore_scenes(session,_self=cmd):
        pymol=_self._pymol        
        if session.has_key('scene_dict'):
            pymol._scene_dict = copy.deepcopy(session['scene_dict'])
            pymol._scene_dict_sc.rebuild(pymol._scene_dict.keys())
        if session.has_key('scene_order'):
            pymol._scene_order = copy.deepcopy(session['scene_order'])
        else:
            pymol._scene_order = []
        _scene_validate_list(_self=_self)         
        return 1


    def stereo(toggle='on', quiet=1, _self=cmd):
        '''
DESCRIPTION

    "stereo" activates or deactives stereo mode.

USAGE

    stereo [toggle]

ARGUMENTS

    toggle = on, off, crosseye, walleye, quadbuffer, sidebyside, or geowall
    
EXAMPLES

    stereo on
    stereo off
    stereo crosseye

NOTES

    "quadbuffer" is the default stereo mode if hardware stereo is available.
    otherwise, "crosseye" is the default.

PYMOL API

    cmd.stereo(string toggle)
        '''
        toggle = stereo_dict[stereo_sc.auto_err(str(toggle),'toggle')]
        r = DEFAULT_ERROR      
        try:
            _self.lock(_self)
            if(toggle>0) : # stereo mode code 
                _self.set("stereo_mode",str(toggle),quiet=quiet)
                toggle=1
            r = _cmd.stereo(_self._COb,toggle)
            if is_error(r):
                print "Error: Selected stereo mode is not available."
        finally:
            _self.unlock(r,_self);
        if _self._raising(r,_self): raise QuietException
        return r


    def turn(axis, angle, _self=cmd):
        '''
DESCRIPTION

    "turn" rotates the camera about one of the three primary axes,
    centered at the origin.

USAGE

    turn axis, angle

EXAMPLES

    turn x, 90
    turn y, 45

PYMOL API

    cmd.turn(string axis, float angle)

SEE ALSO

    move, rotate, translate, zoom, center, clip
        '''
        r = DEFAULT_ERROR      
        try:
            _self.lock(_self)
            r = _cmd.turn(_self._COb,str(axis),float(angle))
        finally:
            _self.unlock(r,_self)
        if _self._raising(r,_self): raise QuietException
        return r


    def full_screen(toggle=-1, _self=cmd):
        '''
DESCRIPTION

    "full_screen" enables or disables full screen mode.  

USAGE

    full_screen [toggle]

EXAMPLES


    full_screen
    full_screen on
    full_screen off

NOTES

    This does not work correctly on all platforms.  If you encounter
    trouble, try using the maximize button on the viewer window
    instead.
    
    '''
        toggle = toggle_dict[toggle_sc.auto_err(str(toggle),'toggle')]
        if thread.get_ident() == pymol.glutThread:
            try: 
                _self.lock(_self)
                r = _cmd.full_screen(_self._COb,int(toggle))
            finally:
                _self.unlock(r,_self)
        else:
            try:
                _self.lock(_self)
                r = _self._do("_cmd.full_screen(_self._COb,%d)"%int(toggle))
            finally:
                _self.unlock(r,_self)
        if _self._raising(r,_self): raise QuietException
        return r


    def rock(mode=-1,_self=cmd):
        '''
DESCRIPTION

    "rock" toggles Y axis rocking.

USAGE

    rock

PYMOL API

    cmd.rock()
        '''
        r = DEFAULT_ERROR      
        try:
            _self.lock(_self)   
            r = _cmd.rock(_self._COb,int(mode))
        finally:
            _self.unlock(r,_self)
        if _self._raising(r,_self): raise QuietException
        return r

    def label(selection="(all)", expression="", quiet=1,_self=cmd):
        '''
DESCRIPTION

    "label" labels one or more atoms in a selection by evaluating an
    Python expression referencing properties for each atom.

USAGE

    label [ selection [, expression ]]

ARGUMENTS

    selection = string: a selection-expression

    expression = string: a Python expression that can be converted to a string
    
EXAMPLES

    label chain A, chain
    label name ca,"%s-%s" % (resn,resi)
    label resi 200,"%1.3f" % partial_charge

NOTES

    The symbols defined in the label name space for each atom are:

        name, resi, resn, resv, chain, segi, model, alt, q, b, type,
        index, rank, ID, ss, vdw, elec_radius, label, elem, geom,
        flags, color, cartoon, valence, formal_charge, partial_charge,
        numeric_type, text_type

    All strings in the expression must be explicitly quoted.

    This operation typically takes several seconds per thousand atoms
    labelled.

    To clear labels, simply omit the expression or set it to ''.

        '''
        # preprocess selection
        selection = selector.process(selection)
        #
        r = DEFAULT_ERROR      
        try:
            _self.lock(_self)
            if len(str(expression))==0:
                r= _cmd.label(_self._COb,"("+str(selection)+")",'',quiet)
            else:
                r = _cmd.label(_self._COb,"("+str(selection)+")",'label='+str(expression),quiet)
        finally:
            _self.unlock(r,_self)   
        if _self._raising(r,_self): raise QuietException
        return r

    def label2(selection="(all)", expression="", quiet=1,_self=cmd):
        # preprocess selection
        selection = selector.process(selection)
        #
        r = DEFAULT_ERROR      
        try:
            _self.lock(_self)
            if len(str(expression))==0:
                r= _cmd.label2(_self._COb,"("+str(selection)+")",'',quiet)
            else:
                r = _cmd.label2(_self._COb,"("+str(selection)+")",str(expression),quiet)
        finally:
            _self.unlock(r,_self)   
        if _self._raising(r,_self): raise QuietException
        return r

    def window(action='show', x=0, y=0, width=0, height=0, _self=cmd):
        '''
DESCRIPTION

    "window" controls the visibility of PyMOL\'s output window

USAGE

    window [ action [, x [, y [, width [, height ]]]]]

NOTES
    
    This command is not fully implemented in MacPyMOL.

PYMOL API

    cmd.window(string action, int x, int y, int width, int height)

        '''
        action = window_sc.auto_err(action,'action')
        action = window_dict[str(action)]

        r = DEFAULT_ERROR      
        try:
            _self.lock(_self)
            r = _cmd.window(_self._COb,action,int(x),int(y),int(width),int(height))
        finally:
            _self.unlock(r,_self)
        if _self._raising(r,_self): raise QuietException
        return r
        
    def viewport(width=-1,height=-1,_self=cmd):
        '''
DESCRIPTION

    "viewport" changes the size of the graphics display area.

USAGE

    viewport width, height

PYMOL API

    cmd.viewport(int width, int height)
        '''
        r = None
        if not cmd.is_glut_thread():
            _self.do("viewport %d,%d"%(int(width),int(height)),0)
        else:
            try:
                _self.lock(_self)
                r = _cmd.viewport(_self._COb,int(width),int(height))
            finally:
                _self.unlock(r,_self)
        if _self._raising(r,_self): raise QuietException
        return r


    def bg_color(color="black", _self=cmd):
        '''
DESCRIPTION

    "bg_color" sets the background color.

USAGE

    bg_color [ color ]

ARGUMENTS

    color = string: color name or number {default: black}

EXAMPLES

    bg_color grey30

    bg_color
    
NOTES

    To obtain a transparent background, "unset opaque_background", and
    then use "ray".
    
SEE ALSO

    set_color, ray
    
PYMOL API

    cmd.bg_color(string color)

        '''
        color = _self._interpret_color(_self,color)
        r = DEFAULT_ERROR      
        try:
            _self.lock(_self)
            r = _cmd.bg_color(_self._COb,str(color))
        finally:
            _self.unlock(r,_self)
        if _self._raising(r,_self): raise QuietException
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
        'putty'       : 7,
    }

    cartoon_sc = Shortcut(cartoon_dict.keys())

    def cartoon(type, selection="(all)", _self=cmd):
        '''
DESCRIPTION

    "cartoon" changes the default cartoon representation for a set of atoms.

USAGE

    cartoon type, selection

ARGUMENTS

    type = automatic, skip, loop, rectangle, oval, tube, arrow, dumbbell
    
PYMOL API

    cmd.cartoon(string type, string selection)

EXAMPLES

    cartoon rectangle, chain A

    cartoon skip, resi 145-156

NOTES

    This command is rarely required since the default "automatic" mode
    chooses cartoons according to the information in the PDB HELIX and
    SHEET records.
    
    '''
        # preprocess selection
        selection = selector.process(selection)
        #
        type = cartoon_dict[cartoon_sc.auto_err(str(type),'type')];
        r = DEFAULT_ERROR      
        try:
            _self.lock(_self)   
            r = _cmd.cartoon(_self._COb,"("+str(selection)+")",int(type))
        finally:
            _self.unlock(r,_self)
        if _self._raising(r,_self): raise QuietException
        return r

    def _ray(width,height,antialias,angle,shift,renderer,quiet,_self=cmd):
        r = DEFAULT_ERROR
        try:
            _self.lock_without_glut(_self)
            try:
                _cmd.set_busy(_self._COb,1)
                r = _cmd.render(_self._COb,int(width),int(height),
                                int(antialias),
                                float(angle),
                                float(shift),int(renderer),
                                int(quiet))
            finally:
                _cmd.set_busy(_self._COb,0)
        finally:
            _self.unlock(r,_self)
        if _self._raising(r,_self): raise QuietException
        return r

    def capture(quiet=1, _self=cmd):
        _self.draw(antialias=-2,quiet=quiet)
        
    def draw(width=0, height=0, antialias=-1, quiet=1, _self=cmd):
        '''
DESCRIPTION

    "draw" creates an OpenGL-based image of the current frame.  

USAGE

    draw [width [,height [,antialias ]]]

ARGUMENTS

    width = integer {default: 0 (current)}

    height = integer {default: 0 (current)}

    antialias = integer {default: -1 (use antialias setting)}
    
EXAMPLES

    draw
    draw 1600

NOTES

    Default width and height are taken from the current viewpoint. If
    one is specified but not the other, then the missing value is
    scaled so as to preserve the current aspect ratio.

    Because this feature uses the OpenGL rendering context to piece
    together the image, it does not work when running in the
    command-line only mode.

    On certain graphics hardware, "unset opaque_background" followed
    by "draw" will produce an image with a transparent background.
    However, better results can usually be obtained using "ray".
    
PYMOL API

    cmd.draw(int width, int height, int antialias, int quiet)

SEE ALSO

    ray, png, save
'''
        # stop movies and sculpting if they're on...
        if _self.get_movie_playing():
            _self.mstop()
        if int(_self.get_setting_legacy("sculpting"))!=0:
            _self.set("sculpting","off",quiet=1)
        # make sure that there aren't any pending display events
        _self.refresh()
        #
        r = DEFAULT_ERROR
        try:
            _self.lock(_self)
            r = _cmd.draw(_self._COb,int(width),int(height),
                          int(antialias),int(quiet))
        finally:
            _self.unlock(r,_self)      
        if _self._raising(r,_self): raise QuietException
        return r

    def ray(width=0, height=0, antialias=-1, angle=0.0, shift=0.0,
            renderer=-1, quiet=1, async=0, _self=cmd):
        '''
DESCRIPTION

    "ray" creates a ray-traced image of the current frame. This
    can take some time (up to several minutes, depending on image
    complexity).

USAGE

    ray [width [,height [,antialias [,angle [,shift [,renderer [,quiet
        [,async ]]]]]]]]]

ARGUMENTS

    width = integer {default: 0 (current)}

    height = integer {default: 0 (current)}

    antialias = integer {default: -1 (use antialias setting)}

    angle = float: y-axis rotation for stereo image generation
    {default: 0.0}

    shift = float: x-axis translation for stereo image generation
    {default: 0.0}

    renderer = -1, 0, 1, or 2: respectively, default, built-in,
    pov-ray, or dry-run {default: 0}
    
    async = 0 or 1: should rendering be done in a background thread?
    
EXAMPLES

    ray
    ray 1024,768
    ray renderer=2

NOTES

    Default width and height are taken from the current viewpoint. If
    one is specified but not the other, then the missing value is
    scaled so as to preserve the current aspect ratio.
    
    angle and shift can be used to generate matched stereo pairs

    renderer = 1 uses PovRay.  This is Unix-only and you must have
        "povray" in your path.  It utilizes two two temporary files:
        "tmp_pymol.pov" and "tmp_pymol.png".

    See "help faster" for optimization tips with the builtin renderer.
    See "help povray" for how to use PovRay instead of PyMOL\'s
    built-in ray-tracing engine.

PYMOL API

    cmd.ray(int width, int height, int antialias, float angle,
            float shift, int renderer, int quiet, int async)
SEE ALSO

    draw, png, save
        '''
        arg_tup = (int(width),int(height),
                   int(antialias),float(angle),
                   float(shift),int(renderer),int(quiet),_self)
        # stop movies, rocking, and sculpting if they're on...
        if _self.get_movie_playing():
            _self.mstop()
        if int(_self.get_setting_legacy("sculpting"))!=0:
            _self.set("sculpting","off",quiet=1)
        if _self.rock(-2)>0:
            _self.rock(0)
        #
        r = DEFAULT_ERROR
        if not async:
            r = apply(_ray, arg_tup)
        else:
            render_thread = threading.Thread(target=_ray, args=arg_tup)
            render_thread.setDaemon(1)
            render_thread.start()
            r = DEFAULT_SUCCESS
        if _self._raising(r,_self): raise QuietException
        return r

    def refresh(_self=cmd):
        '''
DESCRIPTION

    "refresh" causes the scene to be redrawn as soon as the operating
    system allows it to be done.

USAGE

    refresh

PYMOL API

    cmd.refresh()

SEE ALSO

    rebuild
        '''
        r = None
        if thread.get_ident() == pymol.glutThread:
            try:
                _self.lock(_self)
                r = _cmd.refresh_now(_self._COb)
            finally:
                _self.unlock(r,_self)
        else:
            try:
                _self.lock(_self)
                r = _self._do("_ cmd._refresh(_self=cmd)",_self=_self)
            finally:
                _self.unlock(r,_self)
        if _self._raising(r,_self): raise QuietException
        return r

    def reset(object='',_self=cmd):
        '''
DESCRIPTION

    "reset" restores the rotation matrix to identity, sets the origin
    to the center of mass (approx.) and zooms the window and clipping
    planes to cover all objects.  Alternatively, it can reset object
    matrices.

USAGE

    reset [ object ]

PYMOL API

    cmd.reset()
        '''
        r = DEFAULT_ERROR      
        try:
            _self.lock(_self)   
            r = _cmd.reset(_self._COb,0,str(object))
        finally:
            _self.unlock(r,_self)
        if _self._raising(r,_self): raise QuietException
        return r


    def dirty(_self=cmd): # OBSOLETE?
        r = DEFAULT_ERROR      
        try:
            _self.lock(_self)
            r = _cmd.dirty(_self._COb)
        finally:
            _self.unlock(r,_self)
        if _self._raising(r,_self): raise QuietException
        return r

    def meter_reset(_self=cmd):
        '''
DESCRIPTION

    "meter_reset" resets the frames per secound counter.

USAGE

    meter_reset
        '''
        r = DEFAULT_ERROR      
        try:
            _self.lock(_self)   
            r = _cmd.reset_rate(_self._COb)
        finally:
            _self.unlock(r,_self)
        if _self._raising(r,_self): raise QuietException
        return r

    def load_png(filename, movie=1, stereo=-1, quiet=0, _self=cmd):
        '''
DESCRIPTION

    "load_png" loads and displays a PNG file from disk.

USAGE

    load_png filename

NOTES

    If the displayed image is too big for the window, it will be
    reduced 2-fold repeatedly until it fits.
    
    '''
        
        r = DEFAULT_ERROR      
        filename = _self.exp_path(str(filename))
        try:
            _self.lock(_self)
            r = _cmd.load_png(_self._COb,str(filename),int(movie),int(stereo),int(quiet))
        finally:
            _self.unlock(r,_self)
        if _self._raising(r,_self): raise QuietException
        return r


    def rebuild(selection='all',representation='everything',_self=cmd):
        '''
DESCRIPTION

    "rebuild" forces PyMOL to recreate geometric objects in
    case any of them have gone out of sync.

USAGE

    rebuild [selection [, representation ]]

ARGUMENTS

    selection = string {default: all}

    representation = string: {default: everything}

PYMOL API

    cmd.rebuild(string selection, string representation)

SEE ALSO

    refresh
    '''
        selection = selector.process(selection)
        representation = repres_sc.auto_err(representation,'representation')
        repn = repres[representation];
        r = DEFAULT_ERROR      
        try:
            _self.lock(_self)
            r = _cmd.rebuild(_self._COb,selection,repn)
        finally:
            _self.unlock(r,_self)
        if _self._raising(r,_self): raise QuietException
        return r
    
    def recolor(selection='all', representation='everything', _self=cmd):
        '''
DESCRIPTION

    "recolor" forces reapplication of colors to existing objects.
    
USAGE

    recolor [selection [, representation ]]

ARGUMENTS

    selection = string {default: all}

    representation = string {default: everything}
    
NOTES

    This command often needs to be executed after "set_color" has been
    used to redefine one or more existing colors.
    
PYMOL API

    cmd.recolor(string selection = 'all', string representation = 'everything')

SEE ALSO

    color, set_color
    '''
        selection = selector.process(selection)
        representation = repres_sc.auto_err(representation,'representation')
        repn = repres[representation];
        r = DEFAULT_ERROR      
        try:
            _self.lock(_self)
            r = _cmd.recolor(_self._COb,selection,repn)
        finally:
            _self.unlock(r,_self)
        if _self._raising(r,_self): raise QuietException
        return r


    def color(color, selection="(all)", quiet=1, flags=0, _self=cmd):
        '''
DESCRIPTION

    "color" changes the color of objects or atoms.

USAGE

    color color [, selection ]

ARGUMENTS

    color = string: color name or number

    selection = string: selection-expression or name-pattern
    corresponding to the atoms or objects to be colored
    {default: (all)}.

NOTES

    When using color ramps, the ramp can be used as a color.
    
PYMOL API

    cmd.color(string color, string selection, int quiet)

SEE ALSO

    set_color, recolor
    
EXAMPLE 

    color cyan
    color yellow, chain A
    '''
        # preprocess selection
        selection = selector.process(selection)
        color = _self._interpret_color(_self,str(color))
        #
        r = DEFAULT_ERROR      
        try:
            _self.lock(_self)
            r = _cmd.color(_self._COb,str(color),str(selection),int(flags),int(quiet))
        finally:
            _self.unlock(r,_self)
        if _self._raising(r,_self): raise QuietException
        return r

    def spectrum(expression="count", palette="rainbow",
                 selection="(all)", minimum=None, maximum=None,
                 byres=0, quiet=1, _self=cmd):
        
        '''
DESCRIPTION

    "spectrum" colors atoms with a spectrum of colors based on an atomic
    property.
    
USAGE

    spectrum [expression [, palette [, selection [, minimum [, maximum [, byres ]]]]]]

ARGUMENTS

    expression = count, b, q, or pc: respectively, atom count, temperature factor,
    occupancy, or partial charge {default: count}
    
    palette = string: palette name {default: rainbow}

    selection = string: atoms to color {default: (all)}

    minimum = float: {default: None (automatic)}

    maximum = float: {default: None (automatic)}

    byres = integer: controls whether coloring is applied per-residue {default: 0}

EXAMPLES

    spectrum b, blue_red, minimum=10, maximum=50

    spectrum count, rainbow_rev, chain A, byres=1

NOTES

    Available palettes include:

       blue_green blue_magenta blue_red blue_white_green
       blue_white_magenta blue_white_red blue_white_yellow blue_yellow
       cbmr cyan_magenta cyan_red cyan_white_magenta cyan_white_red
       cyan_white_yellow cyan_yellow gcbmry green_blue green_magenta
       green_red green_white_blue green_white_magenta green_white_red
       green_white_yellow green_yellow green_yellow_red magenta_blue
       magenta_cyan magenta_green magenta_white_blue
       magenta_white_cyan magenta_white_green magenta_white_yellow
       magenta_yellow rainbow rainbow2 rainbow2_rev rainbow_cycle
       rainbow_cycle_rev rainbow_rev red_blue red_cyan red_green
       red_white_blue red_white_cyan red_white_green red_white_yellow
       red_yellow red_yellow_green rmbc yellow_blue yellow_cyan
       yellow_cyan_white yellow_green yellow_magenta yellow_red
       yellow_white_blue yellow_white_green yellow_white_magenta
       yellow_white_red yrmbcg

PYMOL API

    def spectrum(string expression, string palette,
                 string selection, float minimum, float maximum,
                 int byres, int quiet)


        '''

        palette = palette_sc.auto_err(palette,'palette')
        
        (prefix,digits,first,last) = palette_dict[str(palette)]

        if (maximum==None) or (minimum==None):
            minimum = 0 # signal to auto-adjust levels
            maximum = -1

        # preprocess selection
        selection = selector.process(selection)
        #
        r = DEFAULT_ERROR      
        try:
            _self.lock(_self)
            r = _cmd.spectrum(_self._COb,str(selection),str(expression),
                                    float(minimum),float(maximum),
                                    int(first),int(last),str(prefix),
                                    int(digits),int(byres),int(quiet))
        finally:
            _self.unlock(r,_self)
        if _self._raising(r,_self): raise QuietException
        return r
    
    def set_color(name, rgb, mode=0, quiet=1, _self=cmd):
        '''
DESCRIPTION

    "set_color" defines a new color using the red, green, and blue
    (RGB) color components.

USAGE

    set_color name, rgb

ARGUMENTS

    name = string: name for the new or existing color

    rgb = list of numbers: [red, green, blue] each and all in the range
    (0.0, 1.0) or (0, 255)
    
EXAMPLES 

    set_color red, [ 1.0, 0.0, 0.0 ]

    set_color yellow, [ 255, 255, 0 ]

NOTES

    PyMOL automatically infers the range based on the input arguments.

    It may be necessary to issue "recolor" command in order to force
    recoloring of existing objects.
    
SEE ALSO

    recolor
    
PYMOL API

    cmd.set_color(string name, list-of-numbers rgb, int mode )

        '''
        r = DEFAULT_ERROR
        if _self.is_string(rgb):
            rgb = safe_list_eval(rgb)
        if not (isinstance(rgb,types.ListType) or isinstance(rgb,types.TupleType)):
            print "Error: color specification must be a list such as [ 1.0, 0.0, 0.0 ]"
        elif len(rgb)!=3:
            print "Error: color specification must be a list such as [ 1.0, 0.0, 0.0 ]"
        else:
            rgb = [float(rgb[0]),float(rgb[1]),float(rgb[2])]
            if (rgb[0]>1.0) or (rgb[1]>1.0) or (rgb[2]>1.0):
                # these days, we'll accept 0-1 or 0-255, so long as [1,1,1] is white
                rgb[0] = rgb[0]/255.0
                rgb[1] = rgb[1]/255.0
                rgb[2] = rgb[2]/255.0            
            try:
                _self.lock(_self)
                if len(rgb)==3:
                    r = _cmd.colordef(_self._COb,str(name),rgb[0],rgb[1],rgb[2],int(mode),int(quiet))
                    _self._invalidate_color_sc(_self)
                else:
                    print "Error: invalid color."
            finally:
                _self.unlock(r,_self)
        if _self._raising(r,_self): raise QuietException
        return r

# Aliases for Mother England.

    colour = color
    set_colour = set_color
    bg_colour = bg_color
    recolour = recolor
    

    import setting
    
