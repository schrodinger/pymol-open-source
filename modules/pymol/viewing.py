#A* -------------------------------------------------------------------
#B* This file contains source code for the PyMOL computer program
#C* Copyright (c) Schrodinger, LLC.
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

from . import colorprinting

if True:

    import sys
    import threading
    import pymol
    from . import selector
    import copy
    from . import parsing
    import re
    cmd = sys.modules["pymol.cmd"]

    from .cmd import _cmd, Shortcut, \
          _feedback,fb_module,fb_mask, \
          repres,repres_sc, is_string, is_list, \
          repmasks,repmasks_sc, \
          toggle_dict,toggle_sc,stereo_dict,stereo_sc, \
          palette_dict, palette_sc, window_dict, window_sc, \
          safe_list_eval, safe_alpha_list_eval, \
          location_code, location_sc, boolean_dict, boolean_sc, \
          DEFAULT_ERROR, DEFAULT_SUCCESS

    palette_colors_dict = {
        'rainbow_cycle'     : 'magenta blue cyan green yellow orange red magenta',
        'rainbow_cycle_rev' : 'magenta red orange yellow green cyan blue magenta',
        'rainbow'           : 'blue cyan green yellow orange red',
        'rainbow_rev'       : 'red orange yellow green cyan blue',
        'rainbow2'          : 'blue cyan green yellow orange red',
        'rainbow2_rev'      : 'red orange yellow green cyan blue',
        'gcbmry'            : 'green cyan blue magenta red yellow',
        'yrmbcg'            : 'yellow red magenta blue cyan green',
        'cbmr'              : 'cyan blue magenta red',
        'rmbc'              : 'red magenta blue cyan',
    }

    rep_list = [ "lines", "sticks", "spheres", "dots", "surface",
                 "mesh", "nonbonded", "nb_spheres", "cartoon",
                 "ribbon", "labels", "slice", "ellipsoids", "volume" ]

    scene_action_sc = Shortcut(['store','recall','clear','insert_before',
                                'insert_after','next','previous',
                                'start', 'update','rename','delete',
                                'order', 'sort', 'first',
                                'append'])
    scene_action_dict = {}
    scene_action_dict_sc = Shortcut([])

    view_sc = Shortcut(['store','recall','clear'])

    def zoom(selection="all", buffer=0.0, state=0, complete=0, animate=0, *, _self=cmd):
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
        with _self.lockcm:
            r = _cmd.zoom(_self._COb,str(selection),float(buffer),
                              int(state)-1,int(complete),float(animate))
        return r

    def center(selection="all", state=0, origin=1, animate=0, *, _self=cmd):
        '''
DESCRIPTION

    "center" translates the window, the clipping slab, and the
    origin to a point centered within the atom selection.

USAGE

    center [ selection [, state [, origin [, animate ]]]]

EXAMPLES

    center chain B
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
        # preprocess selection
        selection = selector.process(selection)
        #
        with _self.lockcm:
            r = _cmd.center(_self._COb,str(selection),int(state)-1,int(origin),float(animate))
        return r

    clip_action_sc = Shortcut([ 'near','far','move','slab','atoms' ])

    def clip(mode, distance, selection=None, state=0, *, _self=cmd):
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
        mode = clip_action_sc.auto_err(str(mode),'mode')
        if selection is not None:
            selection = selector.process(selection)
        else:
            selection = ''
        with _self.lockcm:
            r = _cmd.clip(_self._COb,str(mode),float(distance),
                          str(selection),int(state)-1)
        return r

    def origin(selection="(all)", object=None, position=None, state=0, *, _self=cmd):
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
        # preprocess selection
        selection = selector.process(selection)
        #
        with _self.lockcm:
            if object is None: object=''
            if position is None: position=(0.0,0.0,0.0)
            else:
                if _self.is_string(position):
                    position = safe_list_eval(position)
                selection = ''
            r = _cmd.origin(_self._COb,selection,str(object),
                                 (float(position[0]),
                                  float(position[1]),
                                  float(position[2])
                                  ),int(state)-1)
        return r

    def orient(selection="(all)", state=0, animate=0, *, _self=cmd):
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
        # preprocess selection
        selection = selector.process(selection)
        with _self.lockcm:
            return _cmd.orient(_self._COb,"("+selection+")",int(state)-1,float(animate))

    def move(axis, distance, *, _self=cmd):
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
        with _self.lockcm:
            return _cmd.move(_self._COb,str(axis),float(distance))

    def enable(name='all', parents=0, *, _self=cmd):
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
        if name[0]=='(':
            selection = selector.process(name)
            with _self.lockcm:
                r = _cmd.onoff_by_sele(_self._COb,selection,1)
        else:
            with _self.lockcm:
                r = _cmd.onoff(_self._COb,str(name),1,int(parents));
        return r

    def disable(name='all', *, _self=cmd):
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
        if name[0]=='(':
            selection = selector.process(name)
            with _self.lockcm:
                r = _cmd.onoff_by_sele(_self._COb,selection,0)
        else:
            with _self.lockcm:
                r = _cmd.onoff(_self._COb,str(name),0,0);
        return r

    def _rep_to_repmask(rep):
        repn = 0
        for rep in rep.split():
            rep = repmasks_sc.auto_err(rep, 'representation')
            repn |= repmasks[rep]
        return repn

    def toggle(representation="lines", selection="all", *, _self=cmd):
        '''
DESCRIPTION

    "toggle" toggles the visibility of a representation within a
    selection.
    
USAGE

    toggle [ representation [, selection ]]

ARGUMENTS

    representation = string: named representation {default: lines}

    selection = string: atom selection {default: all}

NOTES

    If the representation is enabled for any atom in the selection, it will
    be turned off.

PYMOL API

    cmd.toggle(string representation, string selection)

SEE ALSO

    show, hide
        '''
        with _self.lockcm:
            if representation == 'object':
                repn = -2
            else:
                repn = _rep_to_repmask(representation)
                # preprocess selection
                selection = selector.process(selection)
            r = _cmd.toggle(_self._COb,str(selection),int(repn));
        return r

    def _showhide(rep, selection, value, _self):
        if not selection and (rep in ("", "all") or '(' in rep or '/' in rep):
            # rep looks like a selection
            selection = rep
            rep = "wire" if value else "everything"

        selection = selector.process(selection) or "all"
        repn = _rep_to_repmask(rep)

        with _self.lockcm:
            r = _cmd.showhide(_self._COb, str(selection), int(repn), value)

        return r

    def show(representation="wire", selection="", *, _self=cmd):
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
    show lines, (name CA+C+N)

SEE ALSO

    hide, enable, disable

'''
        return _showhide(representation, selection, 1, _self)

    def show_as(representation="wire", selection="", *, _self=cmd):
        '''
DESCRIPTION

    "as" turns on and off atom and bond representations.

USAGE

    as representation [, selection ]

ARGUMENTS
    
    representation = lines, spheres, mesh, ribbon, cartoon, sticks,
        dots, surface, labels, extent, nonbonded, nb_spheres, slice,
        extent, slice, dashes, angles, dihedrals, cgo, cell, callback,
        volume or everything

    selection = string {default: all}

EXAMPLES

    as lines, name CA+C+N

    as ribbon

PYMOL API

    cmd.show_as(string representation, string selection)

NOTES

    "selection" can be an object name
    "as" alone will turn on lines and nonbonded and hide everything else.

SEE ALSO

    show, hide, enable, disable
        '''
        return _showhide(representation, selection, 2, _self)

    def hide(representation="everything", selection="", *, _self=cmd):
        '''
DESCRIPTION

    "hide" turns off atom and bond representations.


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
        return _showhide(representation, selection, 0, _self)


    def get_view(output=1, quiet=1, *, _self=cmd):
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

        with _self.lockcm:
            r = _cmd.get_view(_self._COb)

        if True:
            output = int(output)
            if True:
                if (_self.get_setting_int("logging") != 0) and (output<3):
                    if not quiet:
                        print(" get_view: matrix written to log file.")
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
                    print("### cut below here and paste into script ###")
                    print("set_view (\\")
                    print("  %14.9f, %14.9f, %14.9f,\\"%r[0:3])
                    print("  %14.9f, %14.9f, %14.9f,\\"%r[4:7])
                    print("  %14.9f, %14.9f, %14.9f,\\"%r[8:11])
                    print("  %14.9f, %14.9f, %14.9f,\\"%r[16:19])
                    print("  %14.9f, %14.9f, %14.9f,\\"%r[19:22])
                    print("  %14.9f, %14.9f, %14.9f )"%r[22:25])
                    print("### cut above here and paste into script ###")
            if output==3:
                return ("set_view (\\\n"+
                  "  %14.9f, %14.9f, %14.9f,\\\n"%r[0:3] +
                  "  %14.9f, %14.9f, %14.9f,\\\n"%r[4:7] +
                  "  %14.9f, %14.9f, %14.9f,\\\n"%r[8:11] +
                  "  %14.9f, %14.9f, %14.9f,\\\n"%r[16:19] +
                  "  %14.9f, %14.9f, %14.9f,\\\n"%r[19:22] +
                  "  %14.9f, %14.9f, %14.9f )\n"%r[22:25])
            r = r[0:3]+r[4:7]+r[8:11]+r[16:25]
        return r

    def set_view(view,animate=0,quiet=1,hand=1, *, _self=cmd):
        r'''
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
        if isinstance(view, (str, bytes)):
            view = safe_list_eval(view)

        if len(view)!=18:
            raise pymol.CmdException(
                "bad view argument; should be a sequence of 18 floats")

        with _self.lockcm:
            r = _cmd.set_view(_self._COb,(
                    float(view[ 0]),float(view[ 1]),float(view[ 2]),0.0,
                    float(view[ 3]),float(view[ 4]),float(view[ 5]),0.0,
                    float(view[ 6]),float(view[ 7]),float(view[ 8]),0.0,
                    0.0,0.0,0.0,1.0,
                    float(view[ 9]),float(view[10]),float(view[11]),
                    float(view[12]),float(view[13]),float(view[14]),
                    float(view[15]),float(view[16]),float(view[17])),
                    int(quiet),float(animate),int(hand))
        return r

    def view(key, action='recall', animate=-1, *, _self=cmd):
        '''
DESCRIPTION

    "view" saves and restore camera views.

USAGE

    view key [, action [, animate]]
    
ARGUMENTS

    key = string or *

    action = store, recall, clear: {default: recall}

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
                print(" view: stored views:")
                lst = list(pymol._view_dict.keys())
                lst.sort()
                parsing.dump_str_list(lst)

        else:
            action = view_sc.auto_err(action,'action')
            if action=='recall':
                key = pymol._view_dict_sc.auto_err(key,'view')
                _self.set_view(pymol._view_dict[key],animate=animate)
                if _feedback(fb_module.scene,fb_mask.actions,_self): # redundant
                    print(" view: \"%s\" recalled."%key)
            elif (action=='store') or (action=='update'):
                pymol._view_dict_sc.append(key)
                pymol._view_dict[key]=_self.get_view(0)
                if _feedback(fb_module.scene,fb_mask.actions,_self):
                    print(" view: view "+action+"d as \"%s\"."%key)
            elif action=='clear':
                key = pymol._view_dict_sc.auto_err(key,'view')
                if key in pymol._view_dict:
                    del pymol._view_dict[key]
                    pymol._view_dict_sc = Shortcut(pymol._view_dict.keys())
                    if _feedback(fb_module.scene,fb_mask.actions,_self): # redundant
                        print(" view: '%s' deleted."%key)


    def get_viewport(output=1, quiet=1, *, _self=cmd):
        '''
DESCRIPTION

    "get_viewport" returns and optionally prints out the screen viewport size

USAGE

    get_viewport [output]

ARGUMENTS

    output = 0: do not print to screen

    output = 1 {default}: print to screen if not logging and not quiet

    output = 2: force output to screen even if log file is open

PYMOL API

    cmd.get_viewport(output=1, quiet=1)

    '''
        output = int(output)

        with _self.lockcm:
            r = _cmd.get_viewport(_self._COb)

        if _self.get_setting_int("logging") and output < 3:
            _self.log(f"_ viewport {r[0]}, {r[1]}\n", f"cmd.viewport{r}\n")
            if not quiet:
                print(" get_viewport: data written to log file.")
            if output < 2:  # suppress if we have a log file open
                output = 0

        if (0 < output < 3) and not quiet:
            print("### cut below here and paste into script ###")
            print("viewport %4d, %4d" % r)
            print("### cut above here and paste into script ###")

        if output == 3:
            colorprinting.warning(" Warning: get_viewport(3) is deprecated")
            return "viewport ( %4d, %4d )\n" % r

        return r

    def get_vis(_self=cmd):
        with _self.lockcm:
            return _cmd.get_vis(_self._COb)

    def set_vis(dict, *, _self=cmd):
        with _self.lockcm:
            return _cmd.set_vis(_self._COb, dict)

    def get_colorection(key, *, _self=cmd):
        with _self.lockcm:
            return _cmd.get_colorection(_self._COb, key)

    def set_colorection(dict,key, *, _self=cmd):
        with _self.lockcm:
            return _cmd.set_colorection(_self._COb, dict, key)

    def del_colorection(dict,key, *, _self=cmd):
        with _self.lockcm:
            return _cmd.del_colorection(_self._COb, dict, key)

    def get_scene_list(_self=cmd):
        with _self.lockcm:
            return _cmd.get_scene_order(_self._COb)

    def chain_session(_self=cmd):
        import os
        # assumes locked interpreter
        r = 0
        session_file = str(_self.get("session_file"))
        re_pat = re.compile("[0-9]+\.")
        if len(session_file): # find next session file, if it exists
            mo = re_pat.search(session_file)
            if mo is not None:
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

    def scene_order(names,sort=0,location='current',quiet=1, *, _self=cmd):
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

    # if names have spaces
    cmd.scene_order(["name 1", "name 2"])

PYMOL API

    cmd.scene_order(names: Union[list, str], sort: str, location: str)

SEE ALSO

    scene
    '''

        location = location_sc.auto_err(location,'location')
        if is_string(sort):
            sort=boolean_dict[boolean_sc.auto_err(sort,'sort option')]

        if isinstance(names, str):
            names = names.split()

        with _self.lockcm:
            return _cmd.scene_order(_self._COb, names, sort, location)

    def _scene_get_current_message(_self=cmd):
        wiz = _self.get_wizard()
        return '\n'.join(wiz.message) if (wiz is not None
                and wiz.__class__.__name__ == 'Message'
                and hasattr(wiz, 'from_scene')) else None

    def scene_recall_message(message, *, _self=cmd):
        '''
        INTERNAL, DO NOT USE.
        Display a scene message.
        '''
        wiz = _self.get_wizard()
        replace_flag = (wiz is not None
                and wiz.__class__.__name__ == 'Message'
                and hasattr(wiz, 'from_scene'))

        if message:
            if is_string(message):
                message = message.splitlines()
            elif not is_list(message):
                raise TypeError("message %s" % (type(message)))
            wizard_func = _self.replace_wizard if replace_flag else _self.wizard
            wizard_func("message", *message)
            _self.get_wizard().from_scene = 1
        elif replace_flag:
            _self.wizard()

    def scene(key='auto', action='recall', message=None, view=1,
              color=1, active=1, rep=1, frame=1, animate=-1,
              new_key=None, hand=1, quiet=1, sele="all", *, _self=cmd):
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
        action = scene_action_sc.auto_err(action, 'action')

        if is_list(message):
            message = '\n'.join(message)

        # default when called with no arguments
        if key == 'auto':
            if action == 'recall':
                action = 'next'

        # preserve message on update
        if action == 'update':
            if message is None:
                message = _scene_get_current_message(_self)

        # aliases (DEPRECATED)
        if action == 'clear':
            action = 'delete'
        elif action == 'append' or action == 'update':
            action = 'store'

        # presentation auto quit
        if (pymol._scene_quit_on_action == action and
                action in ('next', 'previous') and
                _self.get_setting_boolean("presentation") and
                _self.get_setting_boolean("presentation_auto_quit") and
                _self.get("scene_current_name") == ""):
            if not chain_session(_self):
                _self.quit()

        # call C function
        with _self.lockcm:
            r = _cmd.scene(_self._COb, key, action, message, int(view),
                    int(color), int(active), int(rep), int(frame),
                    float(animate), new_key, int(hand), sele)

        # for presentation auto quit
        pymol._scene_quit_on_action = action

        return r

    def _legacy_scene(key='auto', action='recall', message=None, view=1,
              color=1, active=1, rep=1, frame=1, animate=-1,
              new_key=None, hand=1, quiet=1, *, _self=cmd):
        ''' FOR INTERNAL USE ONLY. Stores and deletes <=1.7.4 compatible scenes. '''
        pymol=_self._pymol
        view = int(view)
        rep = int(rep)
        color = int(color)
        active = int(active)
        frame = int(frame)
        quiet = int(quiet)
        animate = 0

        with _self.lockcm:
            if key=='*':
                if action=='clear':
                    for key in pymol._scene_dict:
                        # free selections
                        scene_list = pymol._scene_dict[key]
                        if len(scene_list)>3:
                            colorection = scene_list[3]
                            if colorection is not None:
                                _self.del_colorection(colorection,key)
                        name = "_scene_"+key+"_*"
                        _self.delete(name)
                else:
                    raise ValueError('action=' + action)
            else:
                if action == 'store':
                    if key in ('new', 'auto'):
                        raise ValueError('key=' + key)
                    if key in pymol._scene_dict:
                        raise RuntimeError('update not supported')
                    if rep:
                        for rep_name in rep_list:
                            name = "_scene_"+key+"_"+rep_name
                            _self.select(name,"rep "+rep_name)
                    if is_string(message):
                        if message:
                            if (message[0:1] in [ '"',"'"] and
                                 message[-1:] in [ '"',"'"]):
                                message=message[1:-1]
                            else:
                                message = message.splitlines()
                    pymol._scene_dict[key] = [
                        _self.get_view(0) if view else None,
                        _self.get_vis() if active else None,
                        _self.get_frame() if frame else None,
                        _self.get_colorection(key) if color else None,
                        1 if rep else None,
                        message,
                    ]
                else:
                    raise ValueError('action=' + action)

    def session_save_views(session, *, _self=cmd):
        pymol=_self._pymol
        session['view_dict']=copy.deepcopy(pymol._view_dict)
        return 1

    def session_restore_views(session, *, _self=cmd):
        pymol=_self._pymol
        if 'view_dict' in session:
            pymol._view_dict=copy.deepcopy(session['view_dict'])
            pymol._view_dict_sc.rebuild(list(pymol._view_dict.keys()))
        return 1

    def session_restore_scenes(session, *, _self=cmd):
        # Restore scenes from old session files (<= 1.7.4)

        if 'scene_dict' in session:
            _self.scene('*', 'clear')

            # save initial scene
            tempname = '_initial_scene'
            while tempname in session['scene_dict']:
                tempname += '_'
            _self.scene(tempname, 'store')

            frame = 0
            if _self.get_movie_playing():
                _self.mstop()
                frame = _self.get_frame()

            for key, data in list(session['scene_dict'].items()):
                _convert_legacy_scene(key, data, _self)

            if frame:
                _self.frame(frame)
                _self.mplay()

            # restore initial scene
            _self.scene(tempname, 'recall', animate=0)
            _self.scene(tempname, 'clear')

        if 'scene_order' in session:
            _self.scene_order(' '.join(session['scene_order']))

        return 1

    def _convert_legacy_scene(key, scene_list, _self=cmd):
        # Create a scene from the given legacy scene list and finally delete
        # the colorection and rep selections.

        scene_list += [None] * 5

        view, active, frame, color, rep = [(0 if x is None else 1)
                for x in scene_list[:5]]

        if frame:
            _self.frame(scene_list[2])

        if view:
            _self.set_view(scene_list[0], 0.0)

        if active:
            _self.disable()
            _self.deselect()
            _self.set_vis(scene_list[1])

        if color:
            _self.set_colorection(scene_list[3], key)
            _self.del_colorection(scene_list[3], key)

        if rep:
            # only atomic representations
            _self.hide('everything', '(*)')
            sele_prefix = _self.get_legal_name('_scene_' + key + '_')
            for rep_name in rep_list:
                _self.show(rep_name, "?" + sele_prefix + rep_name)
            _self.delete(sele_prefix + "*")

        _self.scene(key, 'store', scene_list[5], view, color, active, rep, frame)

    def stereo(toggle='on', quiet=1, *, _self=cmd):
        '''
DESCRIPTION

    "stereo" activates or deactives stereo mode.

USAGE

    stereo [toggle]

ARGUMENTS

    toggle = on, off, crosseye, walleye, quadbuffer, sidebyside, geowall, or openvr
    
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
        with _self.lockcm:
            return _cmd.stereo(_self._COb, toggle)


    def turn(axis, angle, *, _self=cmd):
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
        with _self.lockcm:
            r = _cmd.turn(_self._COb,str(axis),float(angle))
        return r


    def full_screen(toggle=-1, *, _self=cmd):
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
        with _self.lockcm:
            if _self.is_gui_thread():
                return _cmd.full_screen(_self._COb,int(toggle))
            return _self._do("full_screen %s" % (toggle), echo=0)


    def rock(mode=-1, *, _self=cmd):
        '''
DESCRIPTION

    "rock" toggles Y axis rocking.

USAGE

    rock

PYMOL API

    cmd.rock()
        '''
        with _self.lockcm:
            r = _cmd.rock(_self._COb,int(mode))
        return r

    def label(selection="(all)", expression="", quiet=1, *, _self=cmd):
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
    label name CA,"%s-%s" % (resn,resi)
    label resi 200,"%1.3f" % partial_charge

NOTES

    The symbols defined in the label name space for each atom are:

        name, resi, resn, resv, chain, segi, model, alt, q, b, type,
        index, rank, ID, ss, vdw, elec_radius, label, elem, geom,
        flags, color, cartoon, valence, formal_charge, partial_charge,
        numeric_type, text_type, stereo

    All strings in the expression must be explicitly quoted.

    This operation typically takes several seconds per thousand atoms
    labelled.

    To clear labels, simply omit the expression or set it to ''.

        '''
        # preprocess selection
        selection = selector.process(selection)
        #
        with _self.lockcm:
            return _cmd.label(_self._COb, selection, expression, quiet)

    def label2(selection="(all)", expression="", quiet=1, *, _self=cmd):
        # preprocess selection
        selection = selector.process(selection)
        #
        with _self.lockcm:
            return _cmd.label2(_self._COb, selection, expression, quiet)

    def window(action='show', x=0, y=0, width=0, height=0, *, _self=cmd):
        '''
DESCRIPTION

    "window" controls the visibility of PyMOL\'s output window

USAGE

    window [ action [, x [, y [, width [, height ]]]]]

PYMOL API

    cmd.window(string action, int x, int y, int width, int height)

        '''
        action = window_sc.auto_err(action,'action')
        action = window_dict[str(action)]

        with _self.lockcm:
            from pymol.gui import get_qtwindow as getPyMOLWindow
            qt_window = getPyMOLWindow()
            if qt_window:
                r = DEFAULT_SUCCESS
                qt_window.window_cmd(action, int(x),int(y),int(width),int(height))
            else:
                r = _cmd.window(_self._COb,action,int(x),int(y),int(width),int(height))
        return r

    def viewport(width=-1,height=-1, *, _self=cmd):
        '''
DESCRIPTION

    "viewport" changes the size of the graphics display area.

USAGE

    viewport width, height

PYMOL API

    cmd.viewport(int width, int height)
        '''
        if cmd.is_string(width) and height == -1:
            width = _self.safe_eval(width)
            if _self.is_sequence(width):
                colorprinting.warning(" Warning: Tuple-syntax (parentheses) "
                                      "for viewport is deprecated")
                width, height = width

        if not cmd.is_gui_thread():
            _self.do("viewport %d,%d"%(int(width),int(height)),0)
            return None

        with _self.lockcm:
            return _cmd.viewport(_self._COb, int(width), int(height))


    def bg_color(color="black", *, _self=cmd):
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
        with _self.lockcm:
            r = _cmd.bg_color(_self._COb,str(color))
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
        'dash'        : 8,
        'cylinder'    : 9,
    }

    cartoon_sc = Shortcut(cartoon_dict.keys())

    def cartoon(type, selection="(all)", *, _self=cmd):
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
        with _self.lockcm:
            return _cmd.cartoon(_self._COb, selection, int(type))

    def _ray(width,height,antialias,angle,shift,renderer,quiet,_self=cmd):
        r = DEFAULT_ERROR
        try:
            _self.lock_without_glut()
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
            _self.unlock(r)
        return r

    def capture(quiet=1, *, _self=cmd):
        _self.draw(antialias=-2,quiet=quiet)

    def draw(width=0, height=0, antialias=-1, quiet=1, *, _self=cmd):
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
        if _self.get_setting_boolean("sculpting"):
            _self.set("sculpting","off",quiet=1)
        #
        def func():
            with _self.lockcm:
                # make sure that there aren't any pending display events
                # TODO could this be fixed with PYMOL-3328 (SceneUpdate)?
                _cmd.refresh_now(_self._COb)

                return _cmd.draw(_self._COb,int(width),int(height),
                          int(antialias),int(quiet))
        return _self._call_with_opengl_context(func)

    def ray(width=0, height=0, antialias=-1, angle=0.0, shift=0.0,
            renderer=-1, quiet=1, async_=0, _self=cmd, **kwargs):
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
        async_ = int(kwargs.pop('async', async_))

        if kwargs:
            raise pymol.CmdException('unknown argument: ' + ', '.join(kwargs))

        arg_tup = (int(width),int(height),
                   int(antialias),float(angle),
                   float(shift),int(renderer),int(quiet),_self)
        # stop movies, rocking, and sculpting if they're on...
        if _self.get_movie_playing():
            _self.mstop()
        if _self.get_setting_boolean("sculpting"):
            _self.set("sculpting","off",quiet=1)
        if _self.rock(-2)>0:
            _self.rock(0)
        #
        if not async_:
            r = _ray(*arg_tup)
        else:
            render_thread = threading.Thread(target=_ray, args=arg_tup)
            render_thread.setDaemon(1)
            render_thread.start()
            r = DEFAULT_SUCCESS
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
        if _self.is_gui_thread():
            return _self._refresh()
        with _self.lockcm:
            return _self._do("_ cmd._refresh()")

    def reset(object='', *, _self=cmd):
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
        with _self.lockcm:
            return _cmd.reset(_self._COb, str(object))


    def dirty(_self=cmd): # OBSOLETE?
        with _self.lockcm:
            r = _cmd.dirty(_self._COb)
        return r

    def meter_reset(_self=cmd):
        '''
DESCRIPTION

    "meter_reset" resets the frames per secound counter.

USAGE

    meter_reset
        '''
        with _self.lockcm:
            r = _cmd.reset_rate(_self._COb)
        return r

    def load_png(filename, movie=1, stereo=-1, quiet=0, *, _self=cmd):
        '''
DESCRIPTION

    "load_png" loads and displays a PNG file from disk.

USAGE

    load_png filename

NOTES

    If the displayed image is too big for the window, it will be
    reduced 2-fold repeatedly until it fits.
    
    '''

        filename = _self.exp_path(str(filename))
        with _self.lockcm:
            return _cmd.load_png(_self._COb, filename, int(movie), int(stereo),
                                 int(quiet))


    def rebuild(selection='all',representation='everything', *, _self=cmd):
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
        with _self.lockcm:
            return _cmd.rebuild(_self._COb, selection, repn)

    def recolor(selection='all', representation='everything', *, _self=cmd):
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
        with _self.lockcm:
            return _cmd.recolor(_self._COb, selection, repn)


    def color(color, selection="(all)", quiet=1, flags=0, *, _self=cmd):
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

    color_deep, set_color, recolor
    
EXAMPLE 

    color cyan
    color yellow, chain A
    '''
        # preprocess selection
        selection = selector.process(selection)
        color = _self._interpret_color(_self,str(color))

        with _self.lockcm:
            return _cmd.color(_self._COb, str(color), str(selection),
                              int(flags), int(quiet))


    def color_deep(color, name='all', quiet=1, *, _self=cmd):
        '''
DESCRIPTION

    Unset all object and atom level (not global) color settings and
    apply given color.

ARGUMENTS

    color = str: color name or number

    name = str: object name or pattern {default: all}

SEE ALSO

    color, unset_deep
        '''
        from pymol.menu import rep_setting_lists
        _self.unset_deep([s for L in rep_setting_lists for (r, s) in L if s],
                name, updates=0, quiet=quiet)
        _self.color(color, name, quiet=quiet)


    import colorsys
    _spectrumany_interpolations = {
        'hls': (colorsys.rgb_to_hls, colorsys.hls_to_rgb),
        'hsv': (colorsys.rgb_to_hsv, colorsys.hsv_to_rgb),
        'rgb': ((lambda *rgb: rgb), (lambda *rgb: rgb)),
    }

    def spectrumany(expression, colors, selection='(all)', minimum=None,
            maximum=None, quiet=1, interpolation='rgb', *, _self=cmd):
        '''
DESCRIPTION

    Pure python implementation of the spectrum command. Supports arbitrary
    color lists instead of palettes and any numerical atom property which
    works in iterate as expression.

    Non-numeric values (like resn) will be enumerated.

    This is not a separate PyMOL command but is used as a fallback in "spectrum".
        '''
        from . import CmdException

        try:
            from_rgb, to_rgb = _spectrumany_interpolations[interpolation]
        except KeyError:
            raise CmdException('interpolation must be one of {}'.format(
                list(_spectrumany_interpolations)))

        if ' ' not in colors:
            colors = palette_colors_dict.get(colors) or colors.replace('_', ' ')

        quiet, colors = int(quiet), colors.split()

        n_colors = len(colors)
        if n_colors < 2:
            raise CmdException('please provide at least 2 colors')

        col_tuples = [_self.get_color_tuple(i) for i in colors]
        if None in col_tuples:
            raise CmdException('unknown color')

        col_tuples = [from_rgb(*c) for c in col_tuples]

        expression = {'pc': 'partial_charge', 'fc': 'formal_charge',
                'resi': 'resv'}.get(expression, expression)

        if expression == 'count':
            e_list = list(range(_self.count_atoms(selection)))
        else:
            e_list = []
            _self.iterate(selection, 'e_list.append(%s)' % (expression), space=locals())

        try:
            v_list = [float(v) for v in e_list if v is not None]
        except (TypeError, ValueError):
            if not quiet:
                print(' Spectrum: Expression is non-numeric, enumerating values')
            v_list = e_list = list(map(sorted(set(e_list)).index, e_list))

        if not v_list:
            return (0., 0.)

        if minimum is None: minimum = min(v_list)
        if maximum is None: maximum = max(v_list)
        r = minimum, maximum = float(minimum), float(maximum)
        if not quiet:
            print(' Spectrum: range (%.5f to %.5f)' % r)

        val_range = maximum - minimum
        if not val_range:
            _self.color(colors[0], selection)
            return r

        e_it = iter(e_list)
        def next_color():
            v = next(e_it)
            if v is None:
                return False
            v = min(1.0, max(0.0, (float(v) - minimum) / val_range)) * (n_colors - 1)
            i = min(int(v), n_colors - 2)
            p = v - i

            col = [(col_tuples[i+1][j] * p + col_tuples[i][j] * (1.0 - p))
                    for j in range(3)]

            rgb = [int(0xFF * v) for v in to_rgb(*col)]

            return 0x40000000 + rgb[0] * 0x10000 + rgb[1] * 0x100 + rgb[2]

        _self.alter(selection, 'color = next_color() or color', space=locals())
        _self.recolor(selection)

        return r

    def spectrum(expression="count", palette="rainbow",
                 selection="(all)", minimum=None, maximum=None,
                 byres=0, quiet=1, interpolation='rgb', *, _self=cmd):

        '''
DESCRIPTION

    "spectrum" colors atoms with a spectrum of colors based on an atomic
    property.
    
USAGE

    spectrum [expression [, palette [, selection [, minimum [, maximum [, byres ]]]]]]

ARGUMENTS

    expression = count, b, q, or pc: respectively, atom count, temperature factor,
    occupancy, or partial charge {default: count}
    
    palette = string: palette name or space separated list of colors
    {default: rainbow}

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
        palette_hit = palette_sc.shortcut.get(palette)
        if palette_hit:
            palette = palette_hit

        if not expression.replace('_', '').isalpha() or not palette_hit:
            return spectrumany(expression, palette, selection,
                    minimum, maximum, quiet, interpolation, _self=_self)

        (prefix,digits,first,last) = palette_dict[str(palette)]

        if (maximum is None) or (minimum is None):
            minimum = 0 # signal to auto-adjust levels
            maximum = -1

        # preprocess selection
        selection = selector.process(selection)
        #
        with _self.lockcm:
            r = _cmd.spectrum(_self._COb,str(selection),str(expression),
                                    float(minimum),float(maximum),
                                    int(first),int(last),str(prefix),
                                    int(digits),int(byres),int(quiet))
        return r

    def set_color(name, rgb, mode=0, quiet=1, *, _self=cmd):
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
        if isinstance(rgb, (str, bytes)):
            rgb = safe_list_eval(rgb)

        if not isinstance(rgb, (list, tuple)) or len(rgb) != 3:
            raise pymol.CmdException(
                "color specification must be a list such as [ 1.0, 0.0, 0.0 ]")

        rgb = [float(c) for c in rgb]
        if rgb[0] > 1.0 or rgb[1] > 1.0 or rgb[2] > 1.0:
            rgb = [c / 0xFF for c in rgb]

        with _self.lockcm:
            r = _cmd.colordef(_self._COb, str(name), rgb[0], rgb[1], rgb[2],
                              int(mode), int(quiet))
            _self._invalidate_color_sc()
        return r

# Aliases for Mother England.

    colour = color
    set_colour = set_color
    bg_colour = bg_color
    recolour = recolor


def ipython_image(*args, _self=cmd, **kwargs):
    """Render the scene and return the image as an IPython.display.Image.

    All arguments are forwarded to cmd.png().

    @rtype IPython.display.Image
    """
    import os, tempfile
    from IPython.display import Image
    filename = tempfile.mktemp(".png")
    _self.png(filename, *args, **kwargs)
    try:
        return Image(filename)
    finally:
        os.unlink(filename)
