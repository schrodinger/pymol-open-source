#A* -------------------------------------------------------------------
#B* This file contains source code for the PyMOL computer program
#C* Copyright (c) Schrodinger, LLC.
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

from .constants import CURRENT_STATE, ALL_STATES

if True:

    import time
    from . import selector
    import pymol
    cmd = __import__("sys").modules["pymol.cmd"]
    from .cmd import _cmd,lock,unlock,Shortcut, \
          _feedback,fb_module,fb_mask,is_list, \
          DEFAULT_ERROR, DEFAULT_SUCCESS, _raising, is_ok, is_error

    def auto_measure(*, _self=cmd):
        lst = _self.get_names("selections")
        if "pk1" in lst:
            if "pk2" in lst:
                if "pk3" in lst:
                    if "pk4" in lst:
                        _self.dihedral()
                    else:
                        _self.angle()
                else:
                    _self.distance()
        _self.unpick()

    def get_volume_field(objName, state=1, copy=1, *, _self=cmd):
        '''
DESCRIPTION

    EXPERIMENTAL AND SUBJECT TO CHANGE - DO NOT USE
    API only. Get the raw data of a map or volume object.

ARGUMENTS

    objName = str: object name

    state = int: state index {default: 1}

    copy = 0/1: {default: 1} WARNING: only use copy=0 if you know what you're
    doing. copy=0 will return a numpy array which is a wrapper of the internal
    memory. If the internal memory gets freed or reallocated, this wrapper
    will become invalid.
        '''
        with _self.lockcm:
            r = _self._cmd.get_volume_field(_self._COb, objName, int(state) - 1, int(copy))
        return r

    def get_volume_histogram(objName, bins=64, range=None, *, _self=cmd):
        '''
DESCRIPTION

    API ONLY.  Get min, max, mean, stdev and histogram of a map or volume
    object as a list of length bins + 4.
        '''
        with _self.lockcm:
            r = _self._cmd.get_volume_histogram(_self._COb,objName,
                    int(bins), range or (0., 0.))
        return r

    def get_unused_name(prefix="tmp", alwaysnumber=1, *, _self=cmd):
        with _self.lockcm:
            r = _cmd.get_unused_name(_self._COb, prefix, alwaysnumber)
        return r

    def get_modal_draw(*, _self=cmd, quiet=1):
        with _self.lockcm:
            r = _cmd.get_modal_draw(_self._COb)
        return r

    def get_drag_object_name(*, _self=cmd):
        with _self.lockcm:
            r = _cmd.get_drag_object_name(_self._COb)
        return r

    def get_object_matrix(object,state=1, incl_ttt=1, *, _self=cmd):
        '''
DESCRIPTION

    "get_object_matrix" is an unsupported command that may have
    something to do with querying the transformation matrices
    associated with an object
        '''
        object = str(object)
        with _self.lockcm:
            r = _cmd.get_object_matrix(_self._COb,str(object), int(state)-1, int(incl_ttt))
        return r

    def get_object_ttt(object, quiet=1, *, _self=cmd):
        '''
DESCRIPTION

    "get_object_ttt" is an unsupported command
        '''
        quiet = int(quiet)
        with _self.lockcm:
            r = _cmd.get_object_ttt(_self._COb, str(object), -1, quiet)
        if not quiet:
            if r is None:
                print('TTT is None')
            else:
                for i in range(4):
                    if i == 3:
                        print('TTT ---------------------------+---------')
                    print('TTT %8.2f %8.2f %8.2f | %8.2f' % tuple(r[i * 4:i * 4 + 4]))
        return r

    def get_object_settings(object, state=ALL_STATES, quiet=1, *, _self=cmd):
        '''
DESCRIPTION

    "get_object_settings" is an unsupported command
        '''
        with _self.lockcm:
            r = _cmd.get_object_settings(_self._COb, str(object), int(state) - 1)
        return r

    def get_object_list(selection="(all)", quiet=1, *, _self=cmd):
        '''
        
DESCRIPTION

    "get_object_list" is an unsupported command that may have
    something to do with querying the objects covered by a selection.
    '''
        selection = selector.process(selection)
        with _self.lockcm:
            r = _cmd.get_object_list(_self._COb,str(selection))
            if not quiet:
                if(is_list(r)):
                    print(" get_object_list: ",str(r))
        return r

    def get_symmetry(selection="(all)", state=CURRENT_STATE, quiet=1, *, _self=cmd):
        '''
DESCRIPTION

    "get_symmetry" can be used to obtain the crystal
    and spacegroup parameters for a molecule or map.

USAGE

    get_symmetry object-name-or-selection

PYMOL API

    cmd.get_symmetry(string selection, int state, int quiet)


        '''
        selection = selector.process(selection)
        with _self.lockcm:
            r = _cmd.get_symmetry(_self._COb,str(selection),int(state) - 1)
        if not quiet:
            if r:
                        print(" get_symmetry: A     = %7.3f B    = %7.3f C     = %7.3f"%tuple(r[0:3]))
                        print(" get_symmetry: Alpha = %7.3f Beta = %7.3f Gamma = %7.3f"%tuple(r[3:6]))
                        print(" get_symmetry: SpaceGroup = %s"%r[6])
            else:
                        print(" get_symmetry: No symmetry defined.")
        return r

    def get_title(object, state, quiet=1, *, _self=cmd):
        '''
DESCRIPTION

    "get_title" retrieves a text string to the state of a particular
    object which will be displayed when the state is active.

USAGE

    set_title object, state

PYMOL API

    cmd.set_title(string object, int state, string text)

    '''
        with _self.lockcm:
            r = _cmd.get_title(_self._COb,str(object),int(state)-1)
        if not quiet:
                if r is not None:
                    print(" get_title: %s"%r)
        return r


    def angle(name=None, selection1="(pk1)", selection2="(pk2)",
	      selection3="(pk3)", mode=None, label=1, reset=0,
	      zoom=0, state=ALL_STATES, quiet=1, *,
              state1=-3, state2=-3, state3=-3, _self=cmd):

        '''
DESCRIPTION

    "angle" shows the angle formed between any three atoms.

USAGE

    angle [ name [, selection1 [, selection2 [, selection3 ]]]]

NOTES

    "angle" alone will show the angle angle formed by selections (pk1),
    (pk2), (pk3) which can be set using the "PkAt" mouse action
    (typically, Ctrl-middle-click)

PYMOL API

    cmd.angle(string name, string selection1, string selection2,
              string selection3)

SEE ALSO

    distance, dihedral
    '''

        r = DEFAULT_SUCCESS
        if selection1=="(pk1)":
            if "pk1" not in _self.get_names('selections'):
                if _feedback(fb_module.cmd,fb_mask.errors,_self):
                    print("cmd-Error: The 'pk1' selection is undefined.")
                r = DEFAULT_ERROR
        if selection2=="(pk2)":
            if "pk2" not in _self.get_names('selections'):
                if _feedback(fb_module.cmd,fb_mask.errors,_self):
                    print("cmd-Error: The 'pk2' selection is undefined.")
                r = DEFAULT_ERROR
        if selection3=="(pk3)":
            if "pk3" not in _self.get_names('selections'):
                if _feedback(fb_module.cmd,fb_mask.errors,_self):
                    print("cmd-Error: The 'pk3' selection is undefined.")
                r = DEFAULT_ERROR
        if is_ok(r):
            r = DEFAULT_ERROR

            # if unlabeled, then get next name in series

            if name is not None:
                nam=name
            else:
                try:
                    _self.lock(_self)
                    cnt = _self.get_setting_int("dist_counter") + 1
                    r = _self.set("dist_counter", cnt)
                    nam = "angle%02.0f" % cnt
                finally:
                    _self.unlock(r,_self)

            # defaults
            if mode is None:
                mode = 0
            # preprocess selections
            selection1 = selector.process(selection1)
            selection2 = selector.process(selection2)
            selection3 = selector.process(selection3)
            # now do the deed
            try:
                _self.lock(_self)
                if selection2!="same":
                    selection2 = "("+selection2+")"
                if selection3!="same":
                    selection3 = "("+selection3+")"
                r = _cmd.angle(_self._COb,str(nam),"("+str(selection1)+")",
                               str(selection2),
                               str(selection3),
                               int(mode),int(label),int(reset),
                               int(zoom),int(quiet),int(state)-1,
                               int(state1)-1, int(state2)-1, int(state3)-1)
            finally:
                _self.unlock(r,_self)
        if _raising(r,_self): raise pymol.CmdException
        return r

    def dihedral(name=None, selection1="(pk1)", selection2="(pk2)",
		 selection3="(pk3)", selection4="(pk4)", mode=None,
		 label=1, reset=0, zoom=0, state=ALL_STATES, quiet=1, *, _self=cmd):
        '''
DESCRIPTION

    "dihedral" shows dihedral angles formed between any four atoms.

USAGE

    dihedral [ name [, selection1 [, selection2 [, selection3 [, selection4 ]]]]]


NOTES

    "dihedral" alone will show the dihedral angle formed by selections
    (pk1), (pk2), (pk3), and (pk4), which can be set using the "PkAt"
    mouse action (typically, Ctrl-middle-click)

PYMOL API

    cmd.dihedral(string name, string selection1, string selection2,
                 string selection3, string selection4)

SEE ALSO

    distance, angle
    '''
        r = DEFAULT_SUCCESS
        if selection1=="(pk1)":
            if "pk1" not in _self.get_names('selections'):
                if _feedback(fb_module.cmd,fb_mask.errors,_self):
                    print("cmd-Error: The 'pk1' selection is undefined.")
                r = DEFAULT_ERROR
        if selection2=="(pk2)":
            if "pk2" not in _self.get_names('selections'):
                if _feedback(fb_module.cmd,fb_mask.errors,_self):
                    print("cmd-Error: The 'pk2' selection is undefined.")
                r = DEFAULT_ERROR
        if selection3=="(pk3)":
            if "pk3" not in _self.get_names('selections'):
                if _feedback(fb_module.cmd,fb_mask.errors,_self):
                    print("cmd-Error: The 'pk3' selection is undefined.")
                r = DEFAULT_ERROR
        if selection3=="(pk4)":
            if "pk4" not in _self.get_names('selections'):
                if _feedback(fb_module.cmd,fb_mask.errors,_self):
                    print("cmd-Error: The 'pk4' selection is undefined.")
                r = DEFAULT_ERROR
        if is_ok(r):
            r = DEFAULT_ERROR
            # if unlabeled, then get next name in series

            if name is not None:
                nam=name
            else:
                try:
                    _self.lock(_self)
                    cnt = _self.get_setting_int("dist_counter") + 1
                    r = _self.set("dist_counter", cnt)
                    nam = "dihedral%02.0f" % cnt
                finally:
                    _self.unlock(r,_self)

            # defaults
            if mode is None:
                mode = 0
            # preprocess selections
            selection1 = selector.process(selection1)
            selection2 = selector.process(selection2)
            selection3 = selector.process(selection3)
            selection4 = selector.process(selection4)
            # now do the deed
            try:
                _self.lock(_self)
                if selection2!="same":
                    selection2 = "("+selection2+")"
                if selection3!="same":
                    selection3 = "("+selection3+")"
                if selection4!="same":
                    selection4 = "("+selection4+")"
                r = _cmd.dihedral(_self._COb,str(nam),"("+str(selection1)+")",
                                  str(selection2),
                                  str(selection3),
                                  str(selection4),
                                  int(mode),int(label),
                                  int(reset),int(zoom),
                                  int(quiet),int(state)-1)
            finally:
                _self.unlock(r,_self)
        if _raising(r,_self): raise pymol.CmdException
        return r

    def distance(name=None, selection1="(pk1)", selection2="(pk2)",
		 cutoff=None, mode=None, zoom=0, width=None, length=None,
                 gap=None, label=1, quiet=1, reset=0, state=ALL_STATES,
                 state1=-3, state2=-3, *, _self=cmd):

        '''
DESCRIPTION

    "distance" creates a new distance object between two selections.

USAGE
    
    distance [name [, selection1 [, selection2 [, cutoff [, mode ]]]]]

ARGUMENTS

    name = string: name of the distance object to create

    selection1 = string: first atom selection

    selection2 = string: second atom selection

    cutoff = float: longest distance to show 
    
    mode = 0: all interatomic distances

    mode = 1: only bond distances

    mode = 2: only show polar contact distances

    mode = 3: like mode=0, but use distance_exclusion setting

    mode = 4: distance between centroids (does not support
              dynamic_measures; new in PyMOL 1.8.2)

    mode = 5: pi-pi and pi-cation interactions

    mode = 6: pi-pi interactions

    mode = 7: pi-cation interactions

    mode = 8: like mode=3, but cutoff is the ratio between
              distance and sum of VDW radii

    state = int: object state to create the measurement object in
    and to get coordinates from {default: 0 (all states)}

    state1, state2 = int: overrule 'state' argument to measure distances
    between different states {default: use state}

EXAMPLES

    distance mydist, 14/CA, 29/CA

    distance hbonds, all, all, 3.2, mode=2

NOTES

    The distance wizard makes measuring distances easier than using
    the "dist" command for real-time operations.

    "dist" alone will show distances between selections (pk1) and (pk1),
    which can be set using the PkAt mouse action (usually CTRL-middle-click).

PYMOL API

    cmd.distance(string name, string selection1, string selection2,
                 string cutoff, string mode )

    '''
        # handle unnamed distance
        r = DEFAULT_SUCCESS
        if name is not None:
            if len(name):
                if name[0]=='(' or ' ' in name or '/' in name: # we're one argument off...
                    if cutoff is not None:
                        mode = cutoff
                    if selection2!="(pk2)":
                        cutoff = selection2
                    if selection1!="(pk1)":
                        selection2 = selection1
                    selection1=name
                    name = None

        if selection1=="(pk1)":
            if "pk1" not in _self.get_names('selections'):
                if _feedback(fb_module.cmd,fb_mask.errors,_self):
                    print("cmd-Error: The 'pk1' selection is undefined.")
                r = DEFAULT_ERROR
        if selection2=="(pk2)":
            if "pk2" not in _self.get_names('selections'):
                if _feedback(fb_module.cmd,fb_mask.errors,_self):
                    print("cmd-Error: The 'pk2' selection is undefined.")
                r = DEFAULT_ERROR
        if is_ok(r):
            r = DEFAULT_ERROR

            # if unlabeled, then get next name in series

            if name is not None:
                nam=name
            else:
                try:
                    _self.lock(_self)
                    cnt = _self.get_setting_int("dist_counter") + 1
                    r = _self.set("dist_counter", cnt)
                    nam = "dist%02.0f" % cnt
                finally:
                    _self.unlock(r,_self)

            # defaults
            if mode is None:
                mode = 0
            if cutoff is None:
                cutoff = -1.0
            # preprocess selections
            selection1 = selector.process(selection1)
            selection2 = selector.process(selection2)
            # now do the deed
            try:
                _self.lock(_self)
                if selection2!="same":
                    selection2 = "("+selection2+")"
                r = _cmd.dist(_self._COb,str(nam),"("+str(selection1)+")",
                              str(selection2),int(mode),float(cutoff),
                              int(label),int(quiet),int(reset),
                              int(state)-1,int(zoom),
                              int(state1)-1, int(state2)-1)
                if width is not None:
                    _self.set("dash_width",width,nam)
                if length is not None:
                    _self.set("dash_length",length,nam)
                if gap is not None:
                    _self.set("dash_gap",gap,nam)
            finally:
                _self.unlock(r,_self)
        if (r<0.0) and (not quiet):
            # a negative value is an warning signal from PyMOL...
            r = DEFAULT_ERROR
        if _raising(r,_self): raise pymol.CmdException
        return r

    # LEGACY support for cmd.dist
    dist = distance

    def pi_interactions(name="",
                        selection1="all",
                        selection2="same",
                        state=ALL_STATES,
                        state1=-3,
                        state2=-3,
                        quiet=1,
                        reset=0,
                        *, _self=cmd):
        '''
DESCRIPTION

    Find pi-pi and pi-cation interactions.

    Identical to cmd.distance(..., mode=5, label=0)

SEE ALSO

    distance
        '''
        raise pymol.IncentiveOnlyException()

    def get_povray(*, _self=cmd):
        '''
DESCRIPTION

    "get_povray" returns a tuple corresponding to strings for a PovRay
    input file.

PYMOL API

    cmd.get_povray()

        '''
        with _self.lockcm:
            r = _cmd.get_povray(_self._COb)
        return r

    def get_idtf(quiet=1, *, _self=cmd):
        '''
DESCRIPTION

    "get_idft" is under development, but should eventually return an
    idtf file containing multiple objects and scenes.

PYMOL API

    cmd.idtf()

        '''
        with _self.lockcm:
            r = _cmd.get_idtf(_self._COb)

        if not quiet:
            fov = _self.get_setting_float("field_of_view")
            dist = _self.get_view()[11]
            print(" 3Daac=%3.1f, 3Droll=0, 3Dc2c=0 0 1, 3Droo=%1.2f, 3Dcoo=0 0 %1.2f" % (fov, -dist, dist))

        return r

    def get_mtl_obj(*, _self=cmd):
        '''
DESCRIPTION

    NOTE: this is an incomplete and unsupported feature.

    "get_mtl_obj" returns a tuple containing mtl and obj input files
    for use with Maya.

PYMOL API

    cmd.get_obj_mtl()

        '''
        with _self.lockcm:
            r = _cmd.get_mtl_obj(_self._COb)
        return r

    def get_version(quiet=1, *, _self=cmd):
        '''
DESCRIPTION
 
    "get_version" returns a tuple of length six containing text,
    floating point, and integer representations of the current PyMOL
    version number, build date as unix timestamp, GIT SHA and SVN
    code revision so far available.
   
PYMOL API

    cmd.get_version(int quiet)

        '''
        # get_version doesn't need the _COb and doesn't require a lock
        r = _cmd.get_version()
        if True:
            quiet = int(quiet)
            if quiet < 1 and _feedback(fb_module.cmd, fb_mask.results, _self):
                import re
                p = pymol.get_version_message(r)
                print(re.sub(r'^', ' ', p, re.M))
                if quiet < 0:
                    if r[3]:
                        print(' build date:', time.strftime('%c %Z', time.localtime(r[3])))
                    if r[4]:
                        print(' git sha:', r[4])
        return r

    def get_vrml(version=2, *, _self=cmd):
        '''
DESCRIPTION

    "get_vrml" returns a VRML2 string representing the content
    currently displayed.

PYMOL API

    cmd.get_vrml()

        '''
        with _self.lockcm:
            r = _cmd.get_vrml(_self._COb,int(version))
        return r

    def get_collada(version=2, *, _self=cmd):
        '''
DESCRIPTION

    "get_collada" returns a COLLADA string representing the content
    currently displayed.

PYMOL API

    cmd.get_collada()

        '''
        with _self.lockcm:
            r = _cmd.get_collada(_self._COb,int(version))
        return r

    def get_gltf(filename, quiet=1, *, _self=cmd):
        '''
DESCRIPTION

    "get_gltf" saves a gltf file representing the content
    currently displayed.

PYMOL API

    cmd.get_gltf()

        '''
        import shutil
        exe = shutil.which('collada2gltf') or shutil.which('COLLADA2GLTF-bin')
        if exe is None:
            raise pymol.CmdException('could not find collada2gltf')

        # https://github.com/schrodinger/pymol-open-source/issues/107
        _self.set('collada_geometry_mode', 1, quiet=quiet)

        r = _self.get_collada()

        # write collada file
        with open(filename, 'w') as handle:
            handle.write(r)

        import subprocess

        result = subprocess.call([exe, '-i', filename, '-o', filename])
                # convert collada file to gltf by using collada2gltf binary

        if not quiet:
            if result == 0:
                print(' Save: wrote "' + filename + '".')
            else:
                print(' Save-Error: no file written')

        return result

    def count_states(selection="(all)", quiet=1, *, _self=cmd):
        '''
DESCRIPTION

    "count_states" returns the number of states in the selection.

USAGE

    count_states
    
PYMOL API

    cmd.count_states(string selection)

SEE ALSO

    frame
    '''
        # preprocess selection
        selection = selector.process(selection)
        #
        with _self.lockcm:
            r = _cmd.count_states(_self._COb,selection)
        if not quiet:
                print(" cmd.count_states: %d states."%r)
        return r

    def get_movie_length(quiet=1, images=-1, *, _self=cmd):
        '''
DESCRIPTION

    "get_movie_length" returns the number of frames explicitly defined
    in the movie, not including molecular states.

PYMOL API

    cmd.count_frames()

SEE ALSO

    frame, count_states, count_frames
    '''
        with _self.lockcm:
            r = _cmd.get_movie_length(_self._COb)
            if r<0:
                if images==0:
                    r = 0
                elif images<0:
                    r = -r
            if images == 1:
                if r>0:
                    r = 0
            if r>=0 and not quiet:
                print(" cmd.get_movie_length: %d frames"%r)
        return r

    def count_frames(quiet=1, *, _self=cmd):
        '''
DESCRIPTION

    "count_frames" returns the number of frames defined for the PyMOL
    movie.

USAGE

    count_frames
    
PYMOL API

    cmd.count_frames()

SEE ALSO

    frame, count_states
    '''
        with _self.lockcm:
            r = _cmd.count_frames(_self._COb)
            if not quiet: print(" cmd.count_frames: %d frames"%r)
        return r

    def overlap(selection1, selection2, state1=1, state2=1, adjust=0.0, quiet=1, *, _self=cmd):
        '''
DESCRIPTION

    "overlap" is an unsupported command that sums up
    [(VDWi + VDWj) - distance_ij]/2 between pairs of
    selected atoms.

PYMOL API

    cmd.overlap(string selection1, string selection2 [, int state1=1, int state2=1, float adjust=0.0])

NOTES

    It does not compute the volume overlap,
    selections with more atoms will have a larger
    sum.

    '''

        # preprocess selections
        selection1 = selector.process(selection1)
        selection2 = selector.process(selection2)
        #
        with _self.lockcm:
            r = _cmd.overlap(_self._COb,str(selection1),str(selection2),
                                  int(state1)-1,int(state2)-1,
                                  float(adjust))
            if not quiet: print(" cmd.overlap: %5.3f Angstroms."%r)
        return r

    def get_movie_locked(*, _self=cmd):
        with _self.lockcm:
            r = _cmd.get_movie_locked(_self._COb)
        return r

    def get_object_color_index(name, *, _self=cmd):
        name = str(name)
        with _self.lockcm:
            r = _cmd.get_object_color_index(_self._COb,name)
        return r

    def get_color_tuple(name, mode=0, *, _self=cmd):
        '''Get the RGB color tuple for a color (range 0.0 to 1.0)
        :param name: color name or index
        :param mode: don't use (mode 4 returns negative R for special colors)
        '''
        name=str(name)
        mod = int(mode)
        if mode in (1, 2):
            print(' Warning: use get_color_indices instead of get_color_tuple(mode={})'.format(mode))
        elif mode == 3:
            print(' Warning: use get_color_index instead of get_color_tuple(mode=3)')
        with _self.lockcm:
            r = _cmd.get_color(_self._COb,name,mode)
            if r is None:
                if _feedback(fb_module.cmd,fb_mask.errors,_self):
                    print("cmd-Error: Unknown color '%s'."%name)
        return r

    def get_color_indices(all=0, *, _self=cmd):
        with _self.lockcm:
            if all:
                r = _cmd.get_color(_self._COb,'',2)
            else:
                r = _cmd.get_color(_self._COb,'',1)
        return r

    def get_color_index(color,*, _self=cmd):
        with _self.lockcm:
            r = _cmd.get_color(_self._COb,str(color),3)
        return r

    def get_renderer(quiet=1, *, _self=cmd):
        '''
DESCRIPTION

    Prints OpenGL renderer information.
        '''
        with _self.lockcm:
            r = _cmd.get_renderer(_self._COb)

        if not int(quiet):
            print(" OpenGL graphics engine:")
            print("  GL_VENDOR:   ", r[0])
            print("  GL_RENDERER: ", r[1])
            print("  GL_VERSION:  ", r[2])

        return r

    def get_phipsi(selection="(name CA)", state=CURRENT_STATE, *, _self=cmd):
        # preprocess selections
        selection = selector.process(selection)
        #
        with _self.lockcm:
            r = _cmd.get_phipsi(_self._COb,"("+str(selection)+")",int(state)-1)
        return r

    def get_atom_coords(selection, state=ALL_STATES, quiet=1, *, _self=cmd):
        '''
DESCRIPTION

    "get_atom_coords" returns the 3D coordinates of a single atom.
    
    '''
        # low performance way to get coords for a single atom
        selection = selector.process(selection)
        with _self.lockcm:
            r = _cmd.get_atom_coords(_self._COb,str(selection),int(state)-1,int(quiet))
        if not quiet:
            for a in r:
                print(" cmd.get_atom_coords: [%8.3f,%8.3f,%8.3f]"%(a))
        return r

    def get_coords(selection='all', state=1, quiet=1, *, _self=cmd):
        '''
DESCRIPTION

    API only. Get selection coordinates as numpy array.

ARGUMENTS

    selection = str: atom selection {default: all}

    state = int: state index or all states if state=0 {default: 1}
        '''
        selection = selector.process(selection)
        with _self.lockcm:
            r = _cmd.get_coords(_self._COb, selection, int(state) - 1)
            return r

    def get_coordset(name, state=1, copy=1, quiet=1, *, _self=cmd):
        '''
DESCRIPTION

    API only. Get object coordinates as numpy array.

ARGUMENTS

    selection = str: atom selection {default: all}

    state = int: state index {default: 1}

    copy = 0/1: {default: 1} WARNING: only use copy=0 if you know what you're
    doing. copy=0 will return a numpy array which is a wrapper of the internal
    coordinate set memory. If the internal memory gets freed or reallocated,
    this wrapper will become invalid.
        '''
        with _self.lockcm:
            r = _cmd.get_coordset(_self._COb, name, int(state) - 1, int(copy))
            return r


    def get_position(quiet=1, *, _self=cmd):
        '''
DESCRIPTION

    "get_position" returns the 3D coordinates of the center of the
    viewer window.

    '''

        with _self.lockcm:
            r = _cmd.get_position(_self._COb)
        if not quiet:
            print(" cmd.get_position: [%8.3f,%8.3f,%8.3f]"%(r[0],r[1],r[2]))
        return r

    def get_distance(atom1="pk1", atom2="pk2", state=CURRENT_STATE, quiet=1, *, _self=cmd):
        '''
DESCRIPTION

    "get_distance" returns the distance between two atoms.  By default, the
    coordinates used are from the current state, however an alternate
    state identifier can be provided.

USAGE

    get_distance atom1, atom2, [,state ]

EXAMPLES

    get_distance 4/n,4/c
    get_distance 4/n,4/c,state=4
    
PYMOL API

    cmd.get_distance(atom1="pk1",atom2="pk2",state=-1)

        '''
        # preprocess selections
        atom1 = selector.process(atom1)
        atom2 = selector.process(atom2)
        #
        with _self.lockcm:
            r = _cmd.get_distance(_self._COb,str(atom1),str(atom2),int(state)-1)
        if not quiet:
            print(" cmd.get_distance: %5.3f Angstroms."%r)
        return r

    def get_angle(atom1="pk1", atom2="pk2", atom3="pk3", state=CURRENT_STATE, quiet=1, *, _self=cmd):
        '''
DESCRIPTION

    "get_angle" returns the angle between three atoms.  By default, the
    coordinates used are from the current state, however an alternate
    state identifier can be provided.

USAGE

    get_angle atom1, atom2, atom3, [,state ]

EXAMPLES

    get_angle 4/n,4/c,4/ca
    get_angle 4/n,4/c,4/ca,state=4

PYMOL API

    cmd.get_angle(atom1="pk1",atom2="pk2",atom3="pk3",state=-1)

        '''
        # preprocess selections
        atom1 = selector.process(atom1)
        atom2 = selector.process(atom2)
        atom3 = selector.process(atom3)
        #
        with _self.lockcm:
            r = _cmd.get_angle(_self._COb,str(atom1),str(atom2),str(atom3),int(state)-1)
        if not quiet:
            print(" cmd.get_angle: %5.3f degrees."%r)
        return r

    def get_dihedral(atom1="pk1", atom2="pk2", atom3="pk3", atom4="pk4", state=CURRENT_STATE, quiet=1, *, _self=cmd):
        '''
DESCRIPTION

    "get_dihedral" returns the dihedral angle between four atoms.  By
    default, the coordinates used are from the current state, however
    an alternate state identifier can be provided.

    By convention, positive dihedral angles are right-handed
    (looking down the atom2-atom3 axis).

USAGE

    get_dihedral atom1, atom2, atom3, atom4 [,state ]

EXAMPLES

    get_dihedral 4/n,4/c,4/ca,4/cb
    get_dihedral 4/n,4/c,4/ca,4/cb,state=4

PYMOL API

    cmd.get_dihedral(atom1,atom2,atom3,atom4,state=-1)

        '''
        # preprocess selections
        atom1 = selector.process(atom1)
        atom2 = selector.process(atom2)
        atom3 = selector.process(atom3)
        atom4 = selector.process(atom4)
        #
        with _self.lockcm:
            r = _cmd.get_dihe(_self._COb,str(atom1),str(atom2),str(atom3),str(atom4),int(state)-1)
        if not quiet:
            print(" cmd.get_dihedral: %5.3f degrees."%r)
        return r

    def get_model(selection="(all)", state=1, ref='', ref_state=0, *, _self=cmd):
        '''
DESCRIPTION

    "get_model" returns a ChemPy "Indexed" format model from a selection.

PYMOL API

    cmd.get_model(string selection [,int state] )

        '''
        # preprocess selection
        selection = selector.process(selection)
        #
        with _self.lockcm:
            r = _cmd.get_model(_self._COb,"("+str(selection)+")",int(state)-1,str(ref),int(ref_state)-1)
        return r

    def get_bonds(selection="(all)", state=CURRENT_STATE, *, _self=cmd):
        '''
DESCRIPTION

    Get a list of (atm1, atm2, order) tuples for bonds with coordinates in
    the given state (same logic as cmd.get_model()).

    WARNING: atm1/atm2 are 0-based indices! They enumerate the atoms in the
    selection and do not necessarily correspond to the "index" atom property.

    To get a atm1/atm2 to index mapping, you can do:
    >>> stored.indices = []
    >>> cmd.iterate_state(state, selection, "stored.indices.append(index)")

SEE ALSO

    cmd.get_model().bond
        '''
        with _self.lockcm:
            return _cmd.get_bonds(_self._COb, selection, int(state) - 1)

    def get_area(selection="(all)", state=1, load_b=0, quiet=1, *, _self=cmd):
        '''
DESCRIPTION

    Get the surface area of an selection. Depends on the "dot_solvent"
    setting. With "dot_solvent=off" (default) it calculates the solvent
    excluded surface area, else the surface accessible surface.

USAGE

    get_area [ selection [, state [, load_b ]]]

ARGUMENTS

    load_b = bool: store per-atom surface area in b-factors {default: 0}

SEE ALSO

    "dot_solvent" setting, "dots" representation (show dots)
        '''
        # preprocess selection
        selection = selector.process(selection)
        #
        with _self.lockcm:
            r = _cmd.get_area(_self._COb,"("+str(selection)+")",int(state)-1,int(load_b))
        if not quiet:
            print(" cmd.get_area: %5.3f Angstroms^2."%r)
        return r

    def get_chains(selection="(all)", state=ALL_STATES, quiet=1, *, _self=cmd):
        '''
DESCRIPTION

    Print the list of chain identifiers in the given selection.

USAGE

    get_chains [ selection [, state ]]

ARGUMENTS

    selection = str: atom selection {default: all}

    state = int: CURRENTLY IGNORED
        '''
        # preprocess selection
        selection = selector.process(selection)
        #
        with _self.lockcm:
            r = _cmd.get_chains(_self._COb,"("+str(selection)+")",int(state)-1)
        if r is None:
            return []
        if not quiet:
            print(" cmd.get_chains: ",str(r))
        return r

    def get_names(type='public_objects', enabled_only=0, selection="", *, _self=cmd):
        '''
DESCRIPTION

    "get_names" returns a list of object and/or selection names.

PYMOL API

    cmd.get_names( [string: "objects"|"selections"|"all"|"public_objects"|"public_selections"] )

NOTES

    The default behavior is to return only object names.

SEE ALSO

    get_type, count_atoms, count_states
        '''
        selection = selector.process(selection)
        # this needs to be replaced with a flag & masking scheme...
        if type=='objects':
            mode = 1
        elif type=='selections':
            mode = 2
        elif type=='all':
            mode = 0
        elif type=='public':
            mode = 3
        elif type=='public_objects':
            mode = 4
        elif type=='public_selections':
            mode = 5
        elif type=='public_nongroup_objects':
            mode = 6
        elif type=='public_group_objects':
            mode = 7
        elif type=='nongroup_objects':
            mode = 8
        elif type=='group_objects':
            mode = 9
        else:
            raise pymol.CmdException("unknown type: '{}'".format(type))
        with _self.lockcm:
            r = _cmd.get_names(_self._COb,int(mode),int(enabled_only),str(selection))
        return r

    def get_legal_name(name, *, _self=cmd):
        with _self.lockcm:
            r = _cmd.get_legal_name(_self._COb,str(name))
        return r

    def get_type(name, quiet=1, *, _self=cmd):
        '''
DESCRIPTION

    "get_type" returns a string describing the named object or
     selection.

PYMOL API

    cmd.get_type(string object-name)

NOTES

    Possible return values are

    "object:molecule"
    "object:map"
    "object:mesh"
    "object:slice"
    "object:surface"
    "object:measurement"
    "object:cgo"
    "object:group"
    "object:volume"
    "selection"

SEE ALSO

    get_names
        '''
        with _self.lockcm:
            r = _cmd.get_type(_self._COb,str(name))
        if not quiet:
            print(r)
        return r

    def id_atom(selection, mode=0, quiet=1, *, _self=cmd):
        '''
DESCRIPTION

    "id_atom" returns the original source id of a single atom, or
    raises and exception if the atom does not exist or if the selection
    corresponds to multiple atoms.

PYMOL API

    list = cmd.id_atom(string selection)
        '''
        r = DEFAULT_ERROR
        selection = str(selection)
        l = identify(selection, mode)
        ll = len(l)
        if not ll:
            if _feedback(fb_module.cmd,fb_mask.errors,_self):
                print("cmd-Error: atom %s not found by id_atom." % selection)
            if _self._raising(_self=_self): raise pymol.CmdException
        elif ll>1:
            if _feedback(fb_module.cmd,fb_mask.errors,_self):
                print("cmd-Error: multiple atoms %s found by id_atom." % selection)
            if _self._raising(_self=_self): raise pymol.CmdException
        else:
            r = l[0]
            if not quiet:
                if mode:
                    print(" cmd.id_atom: (%s and id %d)"%(r[0],r[1]))
                else:
                    print(" cmd.id_atom: (id %d)"%r)
        if _raising(r,_self): raise pymol.CmdException
        return r

    def identify(selection="(all)", mode=0, quiet=1, *, _self=cmd):
        '''
DESCRIPTION

    "identify" returns a list of atom IDs corresponding to the ID code
    of atoms in the selection.

PYMOL API

    list = cmd.identify(string selection="(all)",int mode=0)

NOTES

    mode 0: only return a list of identifiers (default)
    mode 1: return a list of tuples of the object name and the identifier

        '''
        # preprocess selection
        selection = selector.process(selection)
        #
        with _self.lockcm:
            r = _cmd.identify(_self._COb,"("+str(selection)+")",int(mode)) # 0 = default mode
        if is_list(r):
            if len(r):
                if not quiet:
                    if mode:
                        for a in r:
                            print(" cmd.identify: (%s and id %d)"%(a[0],a[1]))
                    else:
                        for a in r:
                            print(" cmd.identify: (id %d)"%a)
        return r

    def index(selection="(all)", quiet=1, *, _self=cmd):
        '''
DESCRIPTION

    "index" returns a list of tuples corresponding to the
    object name and index of the atoms in the selection.

PYMOL API

    list = cmd.index(string selection="(all)")

NOTE

  Atom indices are fragile and will change as atoms are added
  or deleted.  Whenever possible, use integral atom identifiers
  instead of indices.

        '''
        # preprocess selection
        selection = selector.process(selection)
        #
        with _self.lockcm:
            r = _cmd.index(_self._COb,"("+str(selection)+")",0) # 0 = default mode
        if not quiet:
            if is_list(r):
                if len(r):
                    for a in r:
                        print(" cmd.index: (%s`%d)"%(a[0],a[1]))
        return r

    def find_pairs(selection1, selection2, state1=1, state2=1, cutoff=3.5, mode=0, angle=45, *, _self=cmd):
        '''
DESCRIPTION

    API only function. Returns a list of atom pairs. Atoms are represented as
    (model,index) tuples.

    Can be restricted to hydrogen-bonding-like contacts. WARNING: Only checks
    atom orientation, not atom type (so would hydrogen bond between carbons for
    example), so make sure to provide appropriate atom selections.

ARGUMENTS

    selection1, selection2 = string: atom selections

    state1, state2 = integer: state-index (only positive values allowed) {default: 1}

    cutoff = float: distance cutoff {default: 3.5}

    mode = integer: if mode=1, do coarse hydrogen bonding assessment {default: 0}

    angle = float: hydrogen bond angle cutoff, only if mode=1 {default: 45.0}

NOTE

    Although this does a similar job like "distance", it uses a completely
    different routine and the "mode" argument has different meanings!
        '''
        # preprocess selection
        selection1 = selector.process(selection1)
        selection2 = selector.process(selection2)
        #
        with _self.lockcm:
            r = _cmd.find_pairs(_self._COb,"("+str(selection1)+")",
                                      "("+str(selection2)+")",
                                      int(state1)-1,int(state2)-1,
                                      int(mode),float(cutoff),float(angle))
        return r

    def get_extent(selection="(all)", state=ALL_STATES, quiet=1, *, _self=cmd):
        '''
DESCRIPTION

    "get_extent" returns the minimum and maximum XYZ coordinates of a
    selection as an array:
     [ [ min-X , min-Y , min-Z ],[ max-X, max-Y , max-Z ]]

PYMOL API

    cmd.get_extent(string selection="(all)", state=0 )

        '''
        # preprocess selection
        selection = selector.process(selection)
        #
        with _self.lockcm:
            r = _cmd.get_min_max(_self._COb,str(selection),int(state)-1)
        if not quiet:
            print(" cmd.extent: min: [%8.3f,%8.3f,%8.3f]"%(r[0][0],r[0][1],r[0][2]))
            print(" cmd.extent: max: [%8.3f,%8.3f,%8.3f]"%(r[1][0],r[1][1],r[1][2]))
        return r

    def phi_psi(selection="(byres pk1)", quiet=1, *, _self=cmd):
        '''
DESCRIPTION

    "phi_psi" return the phi and psi angles for a protein atom
    selection.
    
USAGE
        '''

        r = cmd.get_phipsi(selection)
        if not quiet and isinstance(r, dict):
            for a in sorted(r):
                    _self.iterate("(%s`%d)"%a,"print(' %-9s "+
                                ("( %6.1f, %6.1f )"%r[a])+
                                "'%(resn+'-'+resi+':'))")
        return r

    def count_atoms(selection="(all)", quiet=1, state=ALL_STATES, domain='', *, _self=cmd):
        '''
DESCRIPTION

    "count_atoms" returns a count of atoms in a selection.

USAGE

    count_atoms [ selection [, quiet [, state ]]]
        '''
        r = DEFAULT_ERROR
        # preprocess selection
        selection = selector.process(selection)
        #
        try:
            _self.lock(_self)
            r = _cmd.select(_self._COb,"_count_tmp","("+str(selection)+")",1,int(state)-1,str(domain))
            _cmd.delete(_self._COb,"_count_tmp")
        finally:
            _self.unlock(r,_self)
        if not quiet: print(" count_atoms: %d atoms"%r)
        if _raising(r,_self): raise pymol.CmdException
        return r

    def count_discrete(selection, quiet=1, *, _self=cmd):
        '''
DESCRIPTION

    Count the number of discrete objects in selection.

USAGE

    count_discrete selection
        '''
        with _self.lockcm:
            r = _cmd.count_discrete(_self._COb, str(selection))
            if not int(quiet):
                print(' count_discrete: %d' % r)
            return r

    def get_names_of_type(type, public=1, *, _self=cmd):
        """
DESCRIPTION

    "get_names_of_type" will return a list of names for the given type.

        """
        obj_type = "public_objects" if public==1 else "objects"
        types = []
        mix = []
        obj = None
        try:
            obj = _self.get_names(obj_type)
        except:
            pass

        if obj:
            try:
                types = list(map(_self.get_type,obj))
                mix = list(zip(obj,types))
            except:
                pass
        lst = []
        for a in mix:
            if a[1]==type:
                lst.append(a[0])
        return lst

    def get_raw_alignment(name='', active_only=0, *, _self=cmd):
        '''
DESCRIPTION

    "get_raw_alignment" returns a list of lists of (object,index)
    tuples containing the raw per-atom alignment relationships

PYMOL API

    cmd.get_raw_alignment(name)

        '''
        with _self.lockcm:
            r = _cmd.get_raw_alignment(_self._COb,str(name),int(active_only))
        return r

    def get_object_state(name, *, _self=cmd):
        '''
DESCRIPTION

    Returns the effective object state.
        '''
        states = _self.count_states('%' + name)
        if states < 2 and _self.get_setting_boolean('static_singletons'):
            return 1
        if _self.get_setting_int('all_states', name):
            return 0
        state = _self.get_setting_int('state', name)
        if state > states:
            raise pymol.CmdException('Invalid state %d for object %s' % (state, name))
        return state

    def get_selection_state(selection, *, _self=cmd):
        '''
DESCRIPTION

    Returns the effective object state for all objects in given selection.
    Raises exception if objects are in different states.
        '''
        state_set = set(map(_self.get_object_state,
            _self.get_object_list('(' + selection + ')')))
        if len(state_set) != 1:
            if len(state_set) == 0:
                return 1
            raise pymol.CmdException('Selection spans multiple object states')
        return state_set.pop()

    def centerofmass(selection='(all)', state=CURRENT_STATE, quiet=1, *, _self=cmd):
        '''
DESCRIPTION

    Calculates the center of mass. Considers atom mass and occupancy.

ARGUMENTS

    selection = string: atom selection {default: all}

    state = integer: object state, -1 for current state, 0 for all states
    {default: -1}

NOTES

    If occupancy is 0.0 for an atom, set it to 1.0 for the calculation
    (assume it was loaded from a file without occupancy information).

SEE ALSO

    get_extent
        '''
        from chempy import cpv
        state, quiet = int(state), int(quiet)

        if state < 0:
            state = _self.get_selection_state(selection)
        if state == 0:
            states = list(range(1, _self.count_states(selection)+1))
        else:
            states = [state]

        com = cpv.get_null()
        totmass = 0.0

        for state in states:
            model = _self.get_model(selection, state)
            for a in model.atom:
                m = a.get_mass() * (a.q or 1.0)
                com = cpv.add(com, cpv.scale(a.coord, m))
                totmass += m

        if not totmass:
            raise pymol.CmdException('mass is zero')

        com = cpv.scale(com, 1./totmass)
        if not quiet:
            print(' Center of Mass: [%8.3f,%8.3f,%8.3f]' % tuple(com))
        return com

    def cif_get_array(name, key, dtype="s", quiet=1, *, _self=cmd):
        '''
DESCRIPTION

    EXPERIMENTAL AND SUBJECT TO CHANGE!

ARGUMENTS

    name = string: object name

    key = CIF data item name in lower case

    dtype = str: "s" (str), "i" (int) or "f" (float)
        '''
        with _self.lockcm:
            r = _cmd.cif_get_array(_self._COb, name, key, dtype)
        if r and not int(quiet):
            n = len(r)
            r_print = r if n < 10 else (r[:9] + ['... (%d more items)' % (n - 9)])
            print(" %s:" % (key), ', '.join(map(str, r_print)))
        return r

    def get_assembly_ids(name, quiet=1, *, _self=cmd):
        '''
DESCRIPTION

    EXPERIMENTAL AND SUBJECT TO CHANGE!
    Get the list of assembly ids for an object loaded from mmCIF.
        '''
        r = cif_get_array(name, "_pdbx_struct_assembly.id", _self=_self)
        if r and not int(quiet):
            print(" Assembly IDs: %s" % ', '.join(r))
        return r
