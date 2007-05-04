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
          safe_list_eval, lock_without_glut, DEFAULT_ERROR, DEFAULT_SUCCESS
        
    import thread
    
    rep_list = [ "lines","sticks","spheres",
                     "dots","surface","mesh",
                     "nonbonded", "nb_spheres",
                     "cartoon","ribbon","labels","slice"]

    scene_action_sc = Shortcut(['store','recall','clear','insert_before',
                                'insert_after','next','previous',
                                'start', 'update','rename','delete', 'append'])
    scene_action_dict = {}
    scene_action_dict_sc = Shortcut([])

    view_sc = Shortcut(['store','recall','clear'])
    view_dict = {}
    view_dict_sc = Shortcut([])

    scene_dict = {}
    scene_dict_sc = Shortcut([])
    scene_order = []
    scene_counter = 1
    scene_quit_on_action = ''

    def zoom(selection="all", buffer=0.0, state=0, complete=0, animate=0):
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

    selection = a selection-expression or name pattern (default: "all").

    buffer = a floating point number (default: 0): distance
    
    state = 0 (default): uses all coordinate states (default)
    
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
    additional buffer (typically 2 A) to account for the perpective
    transformation and for graphical representations which extend
    beyond the atom coordinates.

SEE ALSO

    origin, orient, center
        '''
        # preprocess selection
        selection = selector.process(selection)
        #
        r = DEFAULT_ERROR
        try:
            lock()   
            r = _cmd.zoom(str(selection),float(buffer),
                              int(state)-1,int(complete),int(animate))
        finally:
            unlock(r)
        if _raising(r): raise QuietException
        return r

    def center(selection="all",state=0,origin=1,animate=0):
        '''
DESCRIPTION

    "center" translates the window, the clipping slab, and the
    origin to a point centered within the atom selection.

USAGE

    center selection=all, state=0, origin=1, animate=0

EXAMPLES

    center 145/

PYMOL API

    cmd.center(string selection, int state, int origin)

ARGUMENTS

    state = 0 (default) use all coordinate states
    state = -1 use only coordinates for the current state
    state > 0  use coordinates for a specific state

    origin = 1 (default) move the origin
    origin = 0 leave the origin unchanged

SEE ALSO

    origin, orient, zoom
        '''
        r = DEFAULT_ERROR
        # preprocess selection
        selection = selector.process(selection)
        #
        try:
            lock()   
            r = _cmd.center(str(selection),int(state)-1,int(origin),int(animate))
        finally:
            unlock(r)
        if _raising(r): raise QuietException
        return r

    clip_action_sc = Shortcut([ 'near','far','move','slab','atoms' ])

    def clip(mode,distance,selection=None,state=0):
        '''
DESCRIPTION

    "clip" alters the near and far clipping planes

USAGE

    clip mode, distance, selection='', state=0

    mode is one of near, far, move, slab, or atoms

EXAMPLES

    clip near, -5           # moves near plane away from you by 5 A
    clip far, 10            # moves far plane towards you by 10 A
    clip move, -5           # moves the slab away from you by 5 A
    clip slab, 20           # sets slab thickness to 20 A
    clip slab, 10, resi 11  # clip 10 A slab about residue 11

    clip atoms, 5, pept     # clip atoms in "pept" with a 5 A buffer
                            # about their current camera positions

PYMOL API

    cmd.clip(string mode, float distance, string selection)

SEE ALSO

    zoom, reset
        '''
        r = DEFAULT_ERROR      
        mode = clip_action_sc.auto_err(str(mode),'mode')
        if selection!=None:
            selection = selector.process(selection)
        else:
            selection = ''
        try:
            lock()   
            r = _cmd.clip(str(mode),float(distance),
                          str(selection),int(state)-1)
        finally:
            unlock(r)
        if _raising(r): raise QuietException         
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

    cmd.origin(string object-or-selection)

ARGUMENTS

    state = 0 (default) use all coordinate states
    state = -1 use only coordinates for the current state
    state > 0  use coordinates for a specific state

SEE ALSO

    zoom, orient, reset
        '''
        #'
        r = DEFAULT_ERROR      
        # preprocess selection
        selection = selector.process(selection)
        #   
        try:
            lock()
            if object==None: object=''
            if position==None: position=(0.0,0.0,0.0)
            else:
                if cmd.is_string(position):
                    position = safe_list_eval(position)
                selection = ''
            r = _cmd.origin(selection,str(object),
                                 (float(position[0]),
                                  float(position[1]),
                                  float(position[2])
                                  ),int(state)-1)
        finally:
            unlock(r)
        if _raising(r): raise QuietException         
        return r

    def orient(selection="(all)",state=0,animate=0):
        '''
DESCRIPTION

    "orient" aligns the principal components of the atoms in the
    selection with the XYZ axes.  The function is similar to the
    orient command in X-PLOR.

USAGE

    orient object-or-selection [, state]
    orient (selection)

PYMOL API

    cmd.orient(string object-or-selection, int state, float animate)

ARGUMENTS

    state = 0 (default) use all coordinate states
    state = -1 use only coordinates for the current state
    state > 0  use coordinates for a specific state

SEE ALSO

    zoom, origin, reset
        '''
        r = DEFAULT_ERROR      
        # preprocess selection
        selection = selector.process(selection)
        #   
        try:
            lock()
            r = _cmd.orient("("+selection+")",int(state)-1,float(animate))
        finally:
            unlock(r)
        if _raising(r): raise QuietException         
        return r

    def move(axis,distance):
        '''
DESCRIPTION

    "move" translates the camera about one of the three primary axes.

USAGE

    move axis,distance

EXAMPLES

    move x,3
    move y,-1

PYMOL API

    cmd.move(string axis, float distance)

SEE ALSO

    turn, rotate, translate, zoom, center, clip
        '''
        r = DEFAULT_ERROR      
        try:
            lock()   
            r = _cmd.move(str(axis),float(distance))
        finally:
            unlock(r)
        if _raising(r): raise QuietException         
        return r

    def enable(name='all'):
        '''
DESCRIPTION

    "enable" turns on display of one or more objects and/or selections.

USAGE

    enable name

ARGUMENTS    

    name = name-pattern or selection.  If name matches a selection
    name, then selection indicator dots are shown for atoms in that
    selection.

NOTES

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
                lock()
                r = _cmd.onoff_by_sele(selection,1)
            finally:
                unlock(r)
        else:
            try:
                lock()   
                r = _cmd.onoff(str(name),1);
            finally:
                unlock(r)
        if _raising(r): raise QuietException            
        return r

    def disable(name='all'):
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
                lock()
                r = _cmd.onoff_by_sele(selection,0)
            finally:
                unlock(r)
        else:
            try:
                lock()   
                r = _cmd.onoff(str(name),0);
            finally:
                unlock(r)
        if _raising(r): raise QuietException            
        return r

    def toggle(representation="", selection=""):
        '''
DESCRIPTION

    "toggle" toggles the visibility of a representation within a
    selection.
    
USAGE

    toggle representation, selection

ARGUMENTS

    representation = na    

NOTES

    If no arguments are provided, then lines are toggled for all
    objects in the aggregate.
    
PYMOL API

    cmd.toggle(string representation, string selection)

        '''
        r = DEFAULT_ERROR
        try:
            lock()
            if (representation=="") and (selection==""):
                r = _cmd.toggle("(all)",repres['lines']); # show lines by default       
            elif (representation!="") and (selection!=""):
                rep = representation
                rep = repres_sc.auto_err(rep,'representation')
                repn = repres[rep];
                # preprocess selection 
                selection = selector.process(selection)
                #   
                r = _cmd.toggle(str(selection),int(repn));
            elif representation=='all':
                r = _cmd.toggle("all",repres['lines']); # toggle lines by default 
            elif (representation[0:1]=='(') or (string.find(representation,'/')>=0):
                # preprocess selection
                selection = selector.process(representation)
                #                  
                r = _cmd.toggle(str(selection),repres['lines']);
            else: # selection==""
                rep = representation
                rep = repres_sc.auto_err(rep,'representation')
                repn = repres[rep];
                r = _cmd.toggle("all",int(repn));
        finally:
            unlock(r)
        if _raising(r): raise QuietException
        return r

    def show(representation="", selection=""):
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

    selection = a selection-expression or name-pattern

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
            lock()
            if (representation=="") and (selection==""):
                if is_ok(_cmd.showhide("(all)",repres['lines'],1)): # show lines by default
                    r = _cmd.showhide("(all)",repres['nonbonded'],2)
            elif (representation!="") and (selection!=""):
                rep = representation
                rep = repres_sc.auto_err(rep,'representation')
                repn = repres[rep];
                # preprocess selection 
                selection = selector.process(selection)
                #   
                r = _cmd.showhide(str(selection),int(repn),1);
            elif representation=='all':
                if is_ok(_cmd.showhide("all",repres['lines'],1)): # show lines by default
                    r = _cmd.showhide("all",repres['nonbonded'], 1) # nonbonded
            elif (representation[0:1]=='(') or (string.find(representation,'/')>=0):
                # preprocess selection
                selection = selector.process(representation)
                #                  
                if is_ok(_cmd.showhide(str(selection),repres['lines'],1)):
                    r = _cmd.showhide(str(selection),repres['nonbonded'],2);
            else: # selection==""
                rep = representation
                rep = repres_sc.auto_err(rep,'representation')
                repn = repres[rep]
                r = _cmd.showhide("all",int(repn),1);
        finally:
            unlock(r)
        if _raising(r): raise QuietException         
        return r

    def show_as(representation="",selection=""):
        '''
DESCRIPTION

    "as" turns on and off atom and bond representations.

    The available representations are:

        lines     spheres    mesh      ribbon     cartoon
        sticks    dots       surface   labels     extent
        nonbonded nb_spheres slice

USAGE

    as
    as reprentation [,object]
    as reprentation [,(selection)]
    as (selection)

PYMOL API

<<<<<<< .mine
    cmd.show_as( string representation="", string selection="" )
=======
    cmd.as(string representation, string selection)
>>>>>>> .r2806

EXAMPLES

    as lines,(name ca or name c or name n)
    as ribbon

NOTES

    "selection" can be an object name
    "as" alone will turn on lines and nonbonded and hide everything else.

SEE ALSO

    show, hide, enable, disable
        '''
        r = DEFAULT_ERROR
        try:
            lock()
            if (representation=="") and (selection==""):
                if is_ok(_cmd.showhide(str(selection),-1,0)):
                    if is_ok(_cmd.showhide("(all)",repres['lines'],1)):
                        r = _cmd.showhide("(all)",repres['nonbonded'],1)
            elif (representation!="") and (selection!=""):
                rep = representation
                rep = repres_sc.auto_err(rep,'representation')
                repn = repres[rep]
                # preprocess selection 
                selection = selector.process(selection)
                #
                if is_ok(_cmd.showhide(str(selection),-1,0)):
                    r = _cmd.showhide(str(selection),int(repn),1)
            elif representation=='all':
                if is_ok(_cmd.showhide(str(selection),-1,0)):            
                    if if_ok(_cmd.showhide("all",repres['lines'],1)): # show lines by default
                        r = _cmd.showhide("all",repres['nonbonded'],1) # show nonbonded by default
            elif (representation[0:1]=='(') or (string.find(representation,'/')>=0):
                # preprocess selection
                selection = selector.process(representation)
                if is_ok(_cmd.showhide(str(selection),-1,0)):
                    r = _cmd.showhide(str(selection),repres['lines'],1)
            else: # selection==""
                rep = representation
                rep = repres_sc.auto_err(rep,'representation')
                repn = repres[rep];
                if is_ok(_cmd.showhide("all",-1,0)):
                    r = _cmd.showhide("all",int(repn),1);
        finally:
            unlock(r)
        if _raising(r): raise QuietException
        return r

    def hide(representation="", selection=""):
        '''
DESCRIPTION

    "hide" turns of atom and bond representations.


USAGE

    hide [ representation [, selection ]]

ARGUMENTS

    representation =  lines, spheres, mesh, ribbon, cartoon,
       sticks, dots, surface, labels, extent, nonbonded, nb_spheres,
       slice, extent, slice, dashes, angles, dihedrals, cgo, cell, callback, 
       or everything

    selection = a selection-expression or name-pattern

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
            unlock(r)
        if _raising(r): raise QuietException         
        return r


    def get_view(output=1,quiet=1):
        '''
DESCRIPTION

    "get_view" returns and optionally prints out the current view
    information in a format which can be embedded into a command
    script and can be used in subsequent calls to "set_view".

    If a log file is currently open, get_view will not write the view
    matrix to the screen unless the "output" parameter is 2.

USAGE

    get_view

PYMOL API

    cmd.get_view(output=1, quiet=1)
    
    my_view= cmd.get_view()

    output control:
    
        0 = output matrix to screen
        1 = do not Output matrix to screen
        2 = force output to screen even if log file is open
        3 = return formatted string instead of a list
        
API USAGE

    cmd.get_view(0) # zero option suppresses output (LEGACY approach)
    cmd.get_view(quiet=1) # suppresses output using PyMOL\'s normal "quiet" parameter.

NOTES

    Contents of the view matrix
        0  -  8 = column-major 3x3 matrix which rotates model axes to camera axes
        9  - 11 = origin or rotation relative to the camera in camera space
        12 - 14 = origin of rotation in model space
        15      = front plane distance from the camera
        16      = rear plane distance from the camera
        17      = orthoscopic flag 

SEE ALSO

    set_view
    '''

        r = DEFAULT_ERROR
        try:
            lock()
            r = _cmd.get_view()
        finally:
            unlock(r)
        if is_ok(r):
            if len(r):
                if (cmd.get_setting_legacy("logging")!=0.0) and (output!=3):
                    if not quiet:
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
                if output and not quiet and (output!=3):
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
        elif _raising(r):
            raise QuietException
        return r

    def set_view(view,animate=0,quiet=1,hand=1):
        '''
DESCRIPTION

    "set_view" sets viewing information for the current scene,
    including the rotation matrix, position, origin of rotation,
    clipping planes, and the orthoscopic flag.

USAGE

    set_view (...)  where ... is 18 floating point numbers

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
                    float(view[15]),float(view[16]),float(view[17])),
		    int(quiet),float(animate),int(hand))
            finally:
                unlock(r)
        if _raising(r): raise QuietException
        return r

    def view(key,action='recall',animate=-1):
        '''
DESCRIPTION

    "view" makes it possible to save and restore viewpoints on a given
    scene within a single session.

USAGE

    view key[,action]
    view *

    key can be any string
    action should be 'store' or 'recall' (default: 'recall')

VIEWS

    Views F1 through F12 are automatically bound to function keys
    provided that "set_key" has not been used to redefine the behaviour
    of the respective key, and that a "scene" has not been defined for
    that key.

EXAMPLES

    view 0,store
    view 0

PYMOL API

    cmd.view(string key, string action)

SEE ALSO

    scene, set_view, get_view
        '''
        global view_dict,view_dict_sc
    
        if key=='*':
            action = view_sc.auto_err(action,'action')
            if action=='clear':
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
                set_view(view_dict[key],animate=animate)
                if _feedback(fb_module.scene,fb_mask.actions): # redundant
                    print " view: \"%s\" recalled."%key
            elif (action=='store') or (action=='update'):
                view_dict_sc.append(key)
                view_dict[key]=cmd.get_view(0)
                if _feedback(fb_module.scene,fb_mask.actions):
                    print " view: view stored as \"%s\"."%key
            elif action=='clear':
                key = view_dict_sc.auto_err(key,'view')
                if view_dict.has_key(key):
                    del view_dict[key]
                    view_dict_sc = Shortcut(view_dict.keys())            
                    if _feedback(fb_module.scene,fb_mask.actions): # redundant
                        print " view: '%s' deleted."%key


    def get_vis():
        r = DEFAULT_ERROR      
        try:
            lock()
            r = _cmd.get_vis()
        finally:
            unlock(r)
        if _raising(r): raise QuietException
        return r

    def set_vis(dict):
        r = DEFAULT_ERROR      
        try:
            lock()         
            r = _cmd.set_vis(dict)
        finally:
            unlock(r)
        if _raising(r): raise QuietException
        return r

    def get_colorection(key):
        r = DEFAULT_ERROR      
        try:
            lock()         
            r = _cmd.get_colorection(key)
        finally:
            unlock(r)
        if _raising(r): raise QuietException
        return r

    def set_colorection(dict,key):
        r = DEFAULT_ERROR      
        try:
            lock()         
            r = _cmd.set_colorection(dict,key)
        finally:
            unlock(r)
        if _raising(r): raise QuietException
        return r

    def set_colorection_name(dict,key,new_key):
        r = DEFAULT_ERROR      
        try:
            lock()         
            r = _cmd.set_colorection_name(dict,key,new_key)
        finally:
            unlock(r)
        if _raising(r): raise QuietException
        return r

    def del_colorection(dict,key):
        r = DEFAULT_ERROR      
        try:
            lock()         
            r = _cmd.del_colorection(dict,key)
        finally:
            unlock(r)
        if _raising(r): raise QuietException
        return r

    def get_scene_dict():
        r = DEFAULT_ERROR      
        try:
            lock()
            cpy = copy.deepcopy(scene_dict)
        finally:
            unlock()
        return cpy

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

    def _scene_get_unique_key():
        global scene_dict,scene_order      
        global scene_counter
        keys = scene_dict.keys()
        while 1:
            key = "%03d"%scene_counter
            if scene_dict.has_key(key):
                scene_counter = scene_counter + 1
            else:
                break;
        return key
    
    def _scene_validate_list():
        global scene_dict,scene_order
        new_list = []
        new_dict = {}
        for a in scene_order:
            if scene_dict.has_key(a) and not new_dict.has_key(a):
                new_list.append(a)
                new_dict[a] = 1
        lst = map(lambda x:(scene_sort_dict.get(x,x),x), scene_dict.keys())
        lst.sort()
        lst = map(lambda x:x[1],lst)
        for a in lst:
            if not new_dict.has_key(a):
                new_list.append(a)
        scene_order = new_list
        r = DEFAULT_ERROR
        try:
            lock()
            r = _cmd._set_scene_names(scene_order)
        finally:
            unlock(r);
        return scene_order

    def chain_session():
        import os
        # assumes locked interpreter
        r = 0
        session_file = str(cmd.get("session_file"))
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
                            cmd.do("_ cmd.load(r'''"+new_file+"''',format='psw')")
                            return 1

                
        return 0
    
    def scene(key='auto', action='recall', message=None, view=1,
              color=1, active=1, rep=1, frame=1, animate=-1,
              new_key=None, hand=1, quiet=1):

        '''
DESCRIPTION

    "scene" makes it possible to save and restore multiple scenes
    scene within a single session.  A scene consists of the view, all
    object activity information, all atom-wise visibility, colors,
    representations, and the global frame index.

USAGE

    scene key [,action [, message, [ new_key=new-key-value ]]]

ARGUMENTS

    key = any-string, new, auto, or *: use new for an automatically
    numbered new scene, use auto for the current scene (if one
    exists), and use * for all scenes (clear and recall actions only).
    
    action = store, recall, insert_after, insert_before, next,
    previous, update, rename, or clear: (default = recall).  If
    rename, then a new_key argument must be explicitly defined.

    message: a text message to display with the scene.

    new_key:  the new name for the scene
    
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
        global scene_dict,scene_dict_sc,scene_order
        global scene_quit_on_action
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
            scene_animation = int(cmd.get_setting_legacy("scene_animation"))
            if scene_animation<0:
                scene_animation = int(cmd.get_setting_legacy("animation"))
            if scene_animation!=0:
                animate = cmd.get_setting_legacy("scene_animation_duration")
            else:
                animate = 0
        try:
            lock() # manipulating global data, so need lock
            if key=='auto':
                action = scene_action_sc.auto_err(action,'action')
                if action=='recall':
                    action='next'
                if action=='start':
                    lst = _scene_validate_list()
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
                    for key in scene_dict.keys():
                        # free selections
                        list = scene_dict[key]
                        if len(list)>3:
                            colorection = list[3]
                            if colorection!=None:
                                cmd.del_colorection(colorection,key) 
                        name = "_scene_"+key+"_*"
                        cmd.delete(name)
                    scene_dict = {}
                    scene_dict_sc = Shortcut(scene_dict.keys())
                    scene_order = []
                    _scene_validate_list()
                else:
                    print " scene: stored scenes:"
                    lst = _scene_validate_list()
                    parsing.dump_str_list(lst)
            else:
                action = scene_action_sc.auto_err(action,'action')
                
                if action=='insert_before':
                    key = _scene_get_unique_key()
                    cur_scene = setting.get("scene_current_name")
                    if cur_scene in scene_order:
                        ix = scene_order.index(cur_scene)
                    else:
                        ix = 0
                    scene_order.insert(ix,key)
                    _scene_validate_list()                                            
                    action='store'
                elif action=='rename':
                    if not scene_dict.has_key(key):
                        print "Error: scene '%s' not found."%key
                    elif new_key==None:
                        print "Error: must provide the 'new_key' argument"
                    else:
                        new_scene_order = []
                        for a in scene_order:
                            if a==key:
                                new_scene_order.append(new_key)
                            else:
                                new_scene_order.append(a)
                        scene_order=new_scene_order
                        scene_dict[new_key] = scene_dict[key]
                        del scene_dict[key]
                        valid_names = cmd.get_names("all")
                        for rep_name in rep_list:
                            name = "_scene_"+key+"_"+rep_name
                            if name in valid_names:
                                new_name = "_scene_"+new_key+"_"+rep_name
                                cmd.set_name(name,new_name)
                        list = scene_dict[new_key]
                        if len(list)>3:
                            cmd.set_colorection_name(list[3],key,new_key)
                        print" scene: '%s' renamed to '%s'."%(key,new_key)
                        scene_dict_sc = Shortcut(scene_dict.keys())
                        cmd.set("session_changed",1,quiet=1)
                elif action=='insert_after':
                    key = _scene_get_unique_key()            
                    cur_scene = setting.get("scene_current_name")
                    if cur_scene in scene_order:
                        ix = scene_order.index(cur_scene) + 1
                    else:
                        ix = len(scene_order)
                    scene_order.insert(ix,key)
                    _scene_validate_list()                    
                    action='store'
                if action=='recall':
                    cmd.set("scenes_changed",1,quiet=1);
                    key = scene_dict_sc.auto_err(key,'scene')
                    cmd.set('scene_current_name', key, quiet=1)               
                    list = scene_dict[key]
                    ll = len(list)
                    if (ll>1) and (active):
                        if list[1]!=None:
                            cmd.disable()
                            cmd.deselect()
                            cmd.set_vis(list[1])
                    if (ll>2) and (frame):
                        if list[2]!=None:
                            if not cmd.get_movie_playing(): # don't set frame when movie is already playing
				if cmd.get_frame()!=list[2]: # only set the frame when it isn't already correct
				    cmd.frame(list[2])
                    if (ll>3) and (color):
                        if list[3]!=None:
                            cmd.set_colorection(list[3],key)
                    if (ll>4) and (rep):
                        if list[4]==None:
                            rep = 0
                    if (ll>5) and (message==None):
                        if list[5]!=None:
                            message=list[5]
                    if rep!=0:
                        cmd.hide("(all)")
                        valid_names = cmd.get_names("all")
                        for rep_name in rep_list:
                            name = "_scene_"+key+"_"+rep_name
                            if name in valid_names:
                                cmd.show(rep_name,name)
                    replace_flag = 0
                    wiz = cmd.get_wizard()
                    if wiz!=None:
                        if str(wiz.__class__) == 'pymol.wizard.message.Message':
                            if hasattr(wiz,'from_scene'):
                                replace_flag = 1
                    mess_flag = 0
                    if message!=None:
                        if is_string(message):
                            if len(message):
                                if(replace_flag):
                                    cmd.replace_wizard("message",message)
                                else:
                                    cmd.wizard("message",message)
                                cmd.get_wizard().from_scene = 1
                                mess_flag = 1
                        if is_list(message):
                            if len(message):
                                if(replace_flag):
                                    apply(cmd.replace_wizard,("message",)+tuple(message))
                                else:
                                    apply(cmd.wizard,("message",)+tuple(message))
                                cmd.get_wizard().from_scene = 1
                                mess_flag = 1
                    if replace_flag and not mess_flag:
                        cmd.wizard()
                    if (ll>0) and (view):
                        if list[0]!=None:
                            set_view(list[0],animate,quiet,hand)
                    if not quiet and _feedback(fb_module.scene,fb_mask.actions): # redundant
                        print " scene: \"%s\" recalled."%key
                elif (action=='store') or (action=='update'):
                    if key =='new':
                        key=_scene_get_unique_key()               
                    if key =='auto':
                        key = setting.get("scene_current_name")
                        if key=='':
                            key=_scene_get_unique_key()
                    if not scene_dict.has_key(key):
                        scene_dict_sc.append(key)
                    else: # get rid of existing one (if exists)
                        list = scene_dict[key]
                        if (action=='update') and (message==None) and len(list)>5:
                            message = scene_dict[key][5]
                        if len(list)>3:
                            colorection = list[3]
                            if colorection!=None:
                                cmd.del_colorection(colorection,key) # important -- free RAM
                        name = "_scene_"+key+"_*"
                        cmd.delete(name)
                    if key not in scene_order:
                        scene_order.append(key)
                    entry = []
                    if view:
                        entry.append(cmd.get_view(0))
                    else:
                        entry.append(None);
                    if active:
                        entry.append(cmd.get_vis())
                    else:
                        entry.append(None)
                    if frame:
                        entry.append(cmd.get_frame())
                    else:
                        entry.append(None)
                    if color:
                        entry.append(cmd.get_colorection(key))
                    else:
                        entry.append(None)
                    if rep:
                        entry.append(1)
                        for rep_name in rep_list:
                            name = "_scene_"+key+"_"+rep_name
                            cmd.select(name,"rep "+rep_name)
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
                    scene_dict[key]=entry
                    if _feedback(fb_module.scene,fb_mask.actions):
                        print " scene: scene stored as \"%s\"."%key
                    _scene_validate_list()                        
                    cmd.set("scenes_changed",1,quiet=1);
                    cmd.set('scene_current_name',key,quiet=1)
                    cmd.set("session_changed",1,quiet=1)
                elif action=='clear':
                    if key=='auto':
                        key = setting.get("scene_current_name")
                    key = scene_dict_sc.auto_err(key,'scene')
                    if scene_dict.has_key(key):
                        list = scene_dict[key]
                        if len(list)>3:
                            colorection = list[3]
                            if colorection!=None:
                                cmd.del_colorection(colorection,key) # important -- free RAM
                        lst = _scene_validate_list()
                        if key == setting.get("scene_current_name"):
                            ix = lst.index(key) - 1
                            if ix>=0:
                                setting.set("scene_current_name",lst[ix],quiet=1)
                        cmd.set("scenes_changed",1,quiet=1);               
                        del scene_dict[key]
                        name = "_scene_"+key+"_*"
                        cmd.delete(name)
                        scene_dict_sc = Shortcut(scene_dict.keys())
                        _scene_validate_list()
                        if _feedback(fb_module.scene,fb_mask.actions):
                            print " scene: '%s' deleted."%key
                    cmd.set("session_changed",1,quiet=1)                                                                    
                elif action=='next':
                    lst = _scene_validate_list()
                    cur_scene = setting.get('scene_current_name',quiet=1)
                    if cur_scene in lst:
                        ix = lst.index(cur_scene) + 1
                    else:
                        ix = 0
                        animate = 0
                        if ((scene_quit_on_action==action) and
                             (setting.get("presentation")=="on") and 
                             (setting.get("presentation_auto_quit")=="on")):
                            cmd.quit()
                    if ix<len(lst):
                        scene_name = lst[ix]
                        cmd.set('scene_current_name', scene_name, quiet=1)
                        scene(scene_name,'recall',animate=animate)
                    elif setting.get("scene_loop")=="on": # loop back to the beginning
                        if len(lst):
                            scene_name = lst[0]
                            cmd.set('scene_current_name', scene_name, quiet=1)
                            scene(scene_name,'recall',animate=animate)
                    else: # otherwise put up blank screen
                        cmd.set('scene_current_name','',quiet=1)
                        chained = 0
                        if (setting.get("presentation")=="on"):
                            chained = chain_session()
                            if (not chained) and (setting.get("presentation_auto_quit")=="on"):
                                scene_quit_on_action = action
                        if not chained: # and len(lst):
                            cmd.disable() # just hide everything
                            cmd.wizard()
                            
                elif action=='previous':
                    lst = _scene_validate_list()            
                    cur_scene = setting.get('scene_current_name',quiet=1)
                    if cur_scene in lst:
                        ix = lst.index(cur_scene) - 1
                    else:
                        ix = len(lst)-1
                        animate = 0
                        if ((scene_quit_on_action==action) and
                             (setting.get("presentation")=="on") and 
                             (setting.get("presentation_auto_quit")=="on")):
                            cmd.quit()
                    if ix>=0:
                        scene_name = lst[ix]
                        scene(scene_name,'recall',animate=animate,hand=-1)
                    elif setting.get("scene_loop")=="on": # loop back to the end
                        print setting.get("scene_loop")
                        if len(lst):
                            scene_name = lst[-1]
                            cmd.set('scene_current_name', scene_name, quiet=1)
                            scene(scene_name,'recall',animate=animate,hand=-1)
                    else: # otherwise put up blank screen
                        cmd.set('scene_current_name','',quiet=1)
                        if ((setting.get("presentation")=="on") and 
                             (setting.get("presentation_auto_quit")=="on")):
                            scene_quit_on_action = action
                        if len(lst):
                            cmd.disable() # just hide everything
                            cmd.wizard()
        finally:
            unlock(r)
        return r
                        
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
        session['scene_order']=copy.deepcopy(scene_order)
        return 1

    def session_restore_scenes(session):
        global scene_dict,scene_dict_sc,scene_order
        if session.has_key('scene_dict'):
            scene_dict = copy.deepcopy(session['scene_dict'])
            scene_dict_sc = Shortcut(scene_dict.keys())
        if session.has_key('scene_order'):
            scene_order = copy.deepcopy(session['scene_order'])
        else:
            scene_order = []
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
    stereo walleye
    stereo quadbuffer

NOTES

    quadbuffer is the default stereo mode if hardware stereo is available
    otherwise, crosseye is the default.

PYMOL API

    cmd.stereo(string state)
        '''
        state = stereo_dict[stereo_sc.auto_err(str(state),'state')]
        r = DEFAULT_ERROR      
        try:
            lock()
            if state>1:
                if state==2: # cross-eye
                    cmd.set("stereo_mode","2",quiet=quiet)
                elif state==3: # quad
                    cmd.set("stereo_mode","1",quiet=quiet)
                elif state==4: # wall-eye
                    cmd.set("stereo_mode","3",quiet=quiet)
                elif state==5: # geowall
                    cmd.set("stereo_mode","4",quiet=quiet)
                elif state==6:
                    cmd.set("stereo_mode","5",quiet=quiet)
                state=1
            r = _cmd.stereo(state)
            if is_error(r):
                print "Error: Selected stereo mode is not available."
        finally:
            unlock(r);
        if _raising(r): raise QuietException
        return r


    def turn(axis,angle):
        '''
DESCRIPTION

    "turn" rotates the camera about one of the three primary axes,
    centered at the origin.

USAGE

    turn axis, angle

EXAMPLES

    turn x,90
    turn y,45

PYMOL API

    cmd.turn(string axis, float angle)

SEE ALSO

    move, rotate, translate, zoom, center, clip
        '''
        r = DEFAULT_ERROR      
        try:
            lock()
            r = _cmd.turn(str(axis),float(angle))
        finally:
            unlock(r)
        if _raising(r): raise QuietException
        return r


    def full_screen(toggle=-1):
        '''
DESCRIPTION

    "full_screen" enables or disables PyMOL\'s full screen mode.  This
    does not work well on all platforms.  

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
                unlock(r)
        else:
            try:
                lock()
                r = cmd._do("_cmd.full_screen(%d)"%int(toggle))
            finally:
                unlock(r)
        if _raising(r): raise QuietException
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
        r = DEFAULT_ERROR      
        try:
            lock()   
            r = _cmd.rock(int(mode))
        finally:
            unlock(r)
        if _raising(r): raise QuietException
        return r

    def label(selection="(all)", expression="", quiet=1):
        '''
DESCRIPTION

    "label" labels one or more atoms in a selection by evaluating an
    Python expression referencing properties for each atom.

USAGE

    label [ selection [, expression ]]

ARGUMENTS

    selection = a selection-expression

    expression = a Python expression that can be converted to a string
    
EXAMPLES

    label chain A, chain
    label name ca,"%s-%s" % (resn,resi)
    label resi 200,"%1.3f" % partial_charge

NOTES

    The symbols defined in the name space are:

        name, resi, resn, resv, chain, segi, model, alt, q, b, type,
        index, rank, ID, ss, vdw, elec_radius, label, elem, geom,
        flags, color, cartoon, valence, formal_charge, partial_charge,
        numeric_type, text_type

    All strings in the expression must be explicitly quoted.  This
    operation typically takes several seconds per thousand atoms
    altered.

    To clear labels, simply omit the expression or set it to ''.

        '''
        # preprocess selection
        selection = selector.process(selection)
        #
        r = DEFAULT_ERROR      
        try:
            lock()
            if len(str(expression))==0:
                r= _cmd.label("("+str(selection)+")",'',quiet)
            else:
                r = _cmd.label("("+str(selection)+")",'label='+str(expression),quiet)
        finally:
            unlock(r)   
        if _raising(r): raise QuietException
        return r

    def window(action='show',x=0,y=0,width=0,height=0):
        '''
DESCRIPTION

    "window" controls the visibility of PyMOL\'s output window

USAGE

    window action

PYMOL API

    cmd.window(string action)

    action = \'show\' or \'hide\'
    
        '''
        action = window_sc.auto_err(action,'action')
        action = window_dict[str(action)]

        r = DEFAULT_ERROR      
        try:
            lock()
            r = _cmd.window(action,int(x),int(y),int(width),int(height))
        finally:
            unlock(r)
        if _raising(r): raise QuietException
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
            cmd.do("viewport %d,%d"%(int(width),int(height)),0)
        else:
            try:
                lock()
                r = _cmd.viewport(int(width),int(height))
            finally:
                unlock(r)
        if _raising(r): raise QuietException
        return r


    def bg_color(color="black",_self=cmd):
        '''
DESCRIPTION

    "bg_color" sets the background color

USAGE

    bg_color [color]

PYMOL API

    cmd.bg_color(string color="black")

        '''
        color = cmd._interpret_color(_self,color)
        r = DEFAULT_ERROR      
        try:
            lock()
            r = _cmd.bg_color(str(color))
        finally:
            unlock(r)
        if _raising(r): raise QuietException
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

    def cartoon(type,selection="(all)"):
        '''
DESCRIPTION

    "cartoon" changes the default cartoon for a set of atoms.

USAGE

    cartoon type, (selection)

    type = skip | automatic | loop | rectangle | oval | tube | arrow | dumbbell

PYMOL API

    cmd.cartoon(string type, string selection)

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
        r = DEFAULT_ERROR      
        try:
            lock()   
            r = _cmd.cartoon("("+str(selection)+")",int(type))
        finally:
            unlock(r)
        if _raising(r): raise QuietException
        return r

    def _ray(width,height,antialias,angle,shift,renderer,quiet):
        r = DEFAULT_ERROR
        try:
            lock_without_glut()
            try:
                _cmd.set_busy(1)
                r = _cmd.render(int(width),int(height),
                                int(antialias),
                                float(angle),
                                float(shift),int(renderer),
                                int(quiet))
            finally:
                _cmd.set_busy(0)
        finally:
            unlock(r)
        if _raising(r): raise QuietException
        return r

    def draw(width=0,height=0,antialias=-1,quiet=1):
        # stop movies and sculpting if they're on...
        if cmd.get_movie_playing():
            cmd.mstop()
        if int(cmd.get_setting_legacy("sculpting"))!=0:
            cmd.set("sculpting","off",quiet=1)
        # make sure that there aren't any pending display events
        cmd.refresh()
        #
        r = DEFAULT_ERROR
        try:
            lock()
            r = _cmd.draw(int(width),int(height),
                          int(antialias),int(quiet))
        finally:
            unlock(r)      
        if _raising(r): raise QuietException
        return r

    def ray(width=0, height=0, antialias=-1, angle=0.0, shift=0.0,
            renderer=-1, quiet=1, async=0):
        '''
DESCRIPTION

    "ray" creates a ray-traced image of the current frame. This
    can take some time (up to several minutes, depending on image
    complexity).

USAGE

    ray [width [,height [,renderer [,antialias [,angle [,shift 
        [,renderer [,quiet [,async ]]]]]]]]]

PYMOL API

    cmd.ray(int width, int height, int antialias, float angle,
            float shift, int renderer, int quiet, int async)

EXAMPLES

    ray
    ray 1024,768
    ray renderer=0

NOTES

    Default width and height are taken from the current viewpoint. If
    one is specified, but not the other, then the missing value is
    scaled so as to preserve the current aspect ratio.
    
    angle and shift can be used to generate matched stereo pairs

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
        arg_tup = (int(width),int(height),
                   int(antialias),float(angle),
                   float(shift),int(renderer),int(quiet))
        # stop movies, rocking, and sculpting if they're on...
        if cmd.get_movie_playing():
            cmd.mstop()
        if int(cmd.get_setting_legacy("sculpting"))!=0:
            cmd.set("sculpting","off",quiet=1)
        if cmd.rock(-2)>0:
            cmd.rock(0)
        #
        r = DEFAULT_ERROR
        if not async:
            r = apply(_ray, arg_tup)
        else:
            render_thread = threading.Thread(target=_ray, args=arg_tup)
            render_thread.setDaemon(1)
            render_thread.start()
            r = DEFAULT_SUCCESS
        if _raising(r): raise QuietException
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
        r = None
        if thread.get_ident() == pymol.glutThread:
	    try:
		lock()
		r = _cmd.refresh_now()
	    finally:
		unlock(r)
        else:
            try:
                lock()
                r = cmd._do("_ cmd._refresh()")
            finally:
                unlock(r)
        if _raising(r): raise QuietException
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
        r = DEFAULT_ERROR      
        try:
            lock()   
            r = _cmd.reset(0,str(object))
        finally:
            unlock(r)
        if _raising(r): raise QuietException
        return r


    def dirty(): # OBSOLETE?
        r = DEFAULT_ERROR      
        try:
            lock()
            r = _cmd.dirty()
        finally:
            unlock(r)
        if _raising(r): raise QuietException
        return r

    def meter_reset():
        '''
DESCRIPTION

    "meter_reset" resets the frames per secound counter

USAGE

    meter_reset
        '''
        r = DEFAULT_ERROR      
        try:
            lock()   
            r = _cmd.reset_rate()
        finally:
            unlock(r)
        if _raising(r): raise QuietException
        return r

    def load_png(filename,movie=1,stereo=-1,quiet=0):
        r = DEFAULT_ERROR      
	filename = cmd.exp_path(str(filename))
        try:
            lock()
            r = _cmd.load_png(str(filename),int(movie),int(stereo),int(quiet))
        finally:
            unlock(r)
        if _raising(r): raise QuietException
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
        selection = selector.process(selection)
        representation = repres_sc.auto_err(representation,'representation')
        repn = repres[representation];
        r = DEFAULT_ERROR      
        try:
            lock()
            r = _cmd.rebuild(selection,repn)
        finally:
            unlock(r)
        if _raising(r): raise QuietException
        return r
    
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
        selection = selector.process(selection)
        representation = repres_sc.auto_err(representation,'representation')
        repn = repres[representation];
        r = DEFAULT_ERROR      
        try:
            lock()
            r = _cmd.recolor(selection,repn)
        finally:
            unlock(r)
        if _raising(r): raise QuietException
        return r

    def color(color, selection="(all)", quiet=1, flags=0, _self=cmd):
        '''
DESCRIPTION

    "color" changes the color of objects or atoms.

USAGE

    color color [, selection]

ARGUMENTS

    color = a color or ramp name

    selection = a selection-expression or name-pattern
    
PYMOL API

    cmd.color(string color, string selection, int quiet)

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
            lock()
            r = _cmd.color(_self._COb,str(color),str(selection),int(flags),int(quiet))
        finally:
            unlock(r)
        if _raising(r): raise QuietException
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
        r = DEFAULT_ERROR      
        try:
            lock()
            r = _cmd.spectrum(str(selection),str(expression),
                                    float(minimum),float(maximum),
                                    int(first),int(last),str(prefix),
                                    int(digits),int(byres),int(quiet))
        finally:
            unlock(r)
        if _raising(r): raise QuietException
        return r
    
    def set_color(name,rgb,mode=0,quiet=1,_self=cmd):
        '''
DESCRIPTION

    "set_color" defines a new color with color indices (0.0-1.0)

USAGE

    set_color name, [ red-float, green-float, blue-float ]

PYMOL API

    cmd.set_color( string name, float-list rgb, int mode )

EXAMPLES 

    set_color red = [ 1.0, 0.0, 0.0 ]
        '''
        r = DEFAULT_ERROR
        if cmd.is_string(rgb):
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
                lock()

                if len(rgb)==3:
                    r = _cmd.colordef(str(name),rgb[0],rgb[1],rgb[2],int(mode),int(quiet))
                    _self._invalidate_color_sc(_self)
                else:
                    print "Error: invalid color."
            finally:
                unlock(r)
        if _raising(r): raise QuietException
        return r

# Aliases for Mother England.

    colour = color
    set_colour = set_color
    bg_colour = bg_color
    recolour = recolor
    

    import setting
    
