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

if __name__=='pymol.querying':

    import selector
    import pymol
    import cmd
    from cmd import _cmd,lock,unlock,Shortcut, \
          _feedback,fb_module,fb_mask,is_list, \
          DEFAULT_ERROR, DEFAULT_SUCCESS, _raising, is_ok, is_error

    def auto_measure(_self=cmd):
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

    def get_volume_field(objName,_self=cmd):
        r = DEFAULT_ERROR
        try:
            _self.lock(_self)
            r = _self._cmd.get_volume_field(_self._COb,objName)
        finally:
            _self.unlock(r,_self)
        return r

    def get_volume_histogram(objName,_self=cmd):
        r = DEFAULT_ERROR
        try:
            _self.lock(_self)
            r = _self._cmd.get_volume_histogram(_self._COb,objName)
        finally:
            _self.unlock(r,_self)
        return r

    def get_volume_ramp(objName,_self=cmd):
        r = DEFAULT_ERROR
        try:
            _self.lock(_self)
            r = _self._cmd.get_volume_ramp(_self._COb,objName)
        finally:
            _self.unlock(r,_self)
        return r

    def get_volume_is_updated(objName,_self=cmd):
        r = DEFAULT_ERROR
        try:
            _self.lock(_self)
            r = _self._cmd.get_volume_is_updated(_self._COb,objName)
        finally:
            _self.unlock(r,_self)
        return r
 
    # SHOULD BE MOVED SOMWHERE ELSE
    def set_volume_ramp(objName, ramp, _self=cmd):
        r = DEFAULT_ERROR
        try:
            _self.lock(_self)
            r = _self._cmd.set_volume_ramp(_self._COb,objName,ramp)
        finally:
            _self.unlock(r,_self)
        return r

    def get_unused_name(prefix="tmp",alwaysnumber=1,_self=cmd):
        r = DEFAULT_ERROR        
        # should replace this with a C function
        try:
            _self.lock(_self)  
            avoid_dict = {}
            for name in _self.get_names('all'):
                avoid_dict[name] = None
            if alwaysnumber or avoid_dict.has_key(prefix):
                counter = 1
                while 1:
                    r = prefix + "%02d"%counter
                    if not avoid_dict.has_key(r):
                        break
                    counter = counter + 1
            else:
                r=prefix

            # safe names....
            safe_chars = range(48,58)
            safe_chars.extend(range(65,91))
            safe_chars.extend(range(97,123))
            safe_chars.extend([ord('+'),ord('-'),ord('_'),ord('.')])
            # remove non-safe chars
            for c in range(len(r)):
                if ord(r[c]) not in safe_chars:
                    r = r[:c] + "_" + r[c+1:]

        finally:
            _self.unlock(r,_self)
        return r
    
    def get_modal_draw(_self=cmd,quiet=1):
        r = DEFAULT_ERROR
        try:
            _self.lock(_self)   
            r = _cmd.get_modal_draw(_self._COb)
        finally:
            _self.unlock(r,_self)
        return r

    def get_drag_object_name(_self=cmd):
        r = DEFAULT_ERROR
        try:
            _self.lock(_self)   
            r = _cmd.get_drag_object_name(_self._COb)
        finally:
            _self.unlock(r,_self)
        return r

    def get_object_matrix(object,state=1,_self=cmd):
        '''
DESCRIPTION

    "get_object_matrix" is an unsupported command that may have
    something to do with querying the transformation matrices
    associated with an object
        '''
        r = DEFAULT_ERROR
        object = str(object)
        try:
            _self.lock(_self)   
            r = _cmd.get_object_matrix(_self._COb,str(object), int(state)-1)
        finally:
            _self.unlock(r,_self)
        if _raising(r,_self): raise pymol.CmdException
        return r

    def get_object_list(selection="(all)", quiet=1, _self=cmd):
        '''
        
DESCRIPTION

    "get_object_list" is an unsupported command that may have
    something to do with querying the objects covered by a selection.
    '''
        
        r = DEFAULT_ERROR
        selection = selector.process(selection)
        try:
            _self.lock(_self)
            r = _cmd.get_object_list(_self._COb,str(selection))
            if not quiet:
                if(is_list(r)):
                    print " get_object_list: ",str(r)
        finally:
            _self.unlock(r,_self)
        if _raising(r,_self): raise pymol.CmdException
        return r

    def get_symmetry(selection="(all)",quiet=1,_self=cmd):
        '''
DESCRIPTION

    "get_symmetry" can be used to obtain the crystal
    and spacegroup parameters for a molecule
    (FUTURE - map object - FUTURE)

USAGE

    get_symmetry object-name-or-selection

PYMOL API

    cmd.get_symmetry(string selection)


        '''
        r = DEFAULT_ERROR
        selection = selector.process(selection)
        try:
            _self.lock(_self)
            r = _cmd.get_symmetry(_self._COb,str(selection))
            if not quiet:
                if(is_list(r)):
                    if(len(r)):
                        print " get_symmetry: A     = %7.3f B    = %7.3f C     = %7.3f"%tuple(r[0:3])
                        print " get_symmetry: Alpha = %7.3f Beta = %7.3f Gamma = %7.3f"%tuple(r[3:6])
                        print " get_symmetry: SpaceGroup = %s"%r[6]
                    else:
                        print " get_symmetry: No symmetry defined."
        finally:
            _self.unlock(r,_self)
        if _raising(r,_self): raise pymol.CmdException
        return r

    def get_title(object, state, quiet=1, _self=cmd):
        '''
DESCRIPTION

    "get_title" retrieves a text string to the state of a particular
    object which will be displayed when the state is active.

USAGE

    set_title object, state

PYMOL API

    cmd.set_title(string object, int state, string text)

    '''
        r = DEFAULT_ERROR
        try:
            _self.lock(_self)
            r = _cmd.get_title(_self._COb,str(object),int(state)-1)
            if not quiet:
                if r!=None:
                    print " get_title: %s"%r      
        finally:
            _self.unlock(r,_self)
        if _raising(r,_self): raise pymol.CmdException
        return r


    def angle(name=None, selection1="(pk1)", selection2="(pk2)", 
	      selection3="(pk3)", mode=None, label=1, reset=0, 
	      zoom=0, state=0, quiet=1, _self=cmd):

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
                    print "cmd-Error: The 'pk1' selection is undefined."
                r = DEFAULT_ERROR
        if selection2=="(pk2)":
            if "pk2" not in _self.get_names('selections'):
                if _feedback(fb_module.cmd,fb_mask.errors,_self):         
                    print "cmd-Error: The 'pk2' selection is undefined."
                r = DEFAULT_ERROR
        if selection3=="(pk3)":
            if "pk3" not in _self.get_names('selections'):
                if _feedback(fb_module.cmd,fb_mask.errors,_self):         
                    print "cmd-Error: The 'pk3' selection is undefined."
                r = DEFAULT_ERROR
        if is_ok(r):
            r = DEFAULT_ERROR
            
            # if unlabeled, then get next name in series

            if name!=None:
                nam=name
            else:
                try:
                    _self.lock(_self)
                    cnt = _cmd.get(_self._COb,"dist_counter") + 1.0
                    r = _cmd.legacy_set(_self._COb,"dist_counter","%1.0f" % cnt)
                    nam = "angle%02.0f" % cnt
                finally:
                    _self.unlock(r,_self)

            # defaults
            if mode == None:
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
                               int(zoom),int(quiet),int(state)-1)
            finally:
                _self.unlock(r,_self)
        if _raising(r,_self): raise pymol.CmdException
        return r

    def dihedral(name=None, selection1="(pk1)", selection2="(pk2)",
		 selection3="(pk3)", selection4="(pk4)", mode=None,
		 label=1, reset=0, zoom=0, state=0, quiet=1, _self=cmd):
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
                    print "cmd-Error: The 'pk1' selection is undefined."
                r = DEFAULT_ERROR
        if selection2=="(pk2)":
            if "pk2" not in _self.get_names('selections'):
                if _feedback(fb_module.cmd,fb_mask.errors,_self):         
                    print "cmd-Error: The 'pk2' selection is undefined."
                r = DEFAULT_ERROR
        if selection3=="(pk3)":
            if "pk3" not in _self.get_names('selections'):
                if _feedback(fb_module.cmd,fb_mask.errors,_self):         
                    print "cmd-Error: The 'pk3' selection is undefined."
                r = DEFAULT_ERROR
        if selection3=="(pk4)":
            if "pk4" not in _self.get_names('selections'):
                if _feedback(fb_module.cmd,fb_mask.errors,_self):         
                    print "cmd-Error: The 'pk4' selection is undefined."
                r = DEFAULT_ERROR
        if is_ok(r):
            r = DEFAULT_ERROR
            # if unlabeled, then get next name in series

            if name!=None:
                nam=name
            else:
                try:
                    _self.lock(_self)
                    cnt = _cmd.get(_self._COb,"dist_counter") + 1.0
                    r = _cmd.legacy_set(_self._COb,"dist_counter","%1.0f" % cnt)
                    nam = "dihedral%02.0f" % cnt
                finally:
                    _self.unlock(r,_self)

            # defaults
            if mode == None:
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
                 gap=None, label=1, quiet=1, reset=0, state=0, _self=cmd):
	
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
        if name!=None:
            if len(name):
                if name[0]=='(' or ' ' in name or '/' in name: # we're one argument off...
                    if cutoff!=None:
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
                    print "cmd-Error: The 'pk1' selection is undefined."
                r = DEFAULT_ERROR
        if selection2=="(pk2)":
            if "pk2" not in _self.get_names('selections'):
                if _feedback(fb_module.cmd,fb_mask.errors,_self):         
                    print "cmd-Error: The 'pk2' selection is undefined."
                r = DEFAULT_ERROR
        if is_ok(r):
            r = DEFAULT_ERROR

            # if unlabeled, then get next name in series

            if name!=None:
                nam=name
            else:
                try:
                    _self.lock(_self)
                    cnt = _cmd.get(_self._COb,"dist_counter") + 1.0
                    r = _cmd.legacy_set(_self._COb,"dist_counter","%1.0f" % cnt)
                    nam = "dist%02.0f" % cnt
                finally:
                    _self.unlock(r,_self)

            # defaults
            if mode == None:
                mode = 0
            if cutoff == None:
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
                              int(state)-1,int(zoom))
                if width!=None:
                    _self.set("dash_width",width,nam)
                if length!=None:
                    _self.set("dash_length",length,nam)
                if gap!=None:
                    _self.set("dash_gap",gap,nam)
            finally:
                _self.unlock(r,_self)
        if (r<0.0) and (not quiet):
            # a negative value is an warning signal from PyMOL...
            r = DEFAULT_ERROR
        if _raising(r,_self): raise pymol.CmdException
        return r

    # LEGACY support for cmd.dist
    def dist(*arg,**kw):
        return apply(distance,arg,kw)

    def get_povray(_self=cmd):
        '''
DESCRIPTION

    "get_povray" returns a tuple corresponding to strings for a PovRay
    input file.

PYMOL API

    cmd.get_povray()

        '''
        r = DEFAULT_ERROR
        try:
            _self.lock(_self)   
            r = _cmd.get_povray(_self._COb)
        finally:
            _self.unlock(r,_self)
        if _raising(r,_self): raise pymol.CmdException
        return r

    def get_idtf(_self=cmd):
        '''
DESCRIPTION

    "get_idft" is under development, but should eventually return an
    idtf file containing multiple objects and scenes.

PYMOL API

    cmd.idtf()

        '''
        r = DEFAULT_ERROR
        try:
            _self.lock(_self)   
            r = _cmd.get_idtf(_self._COb)
        finally:
            _self.unlock(r,_self)
        if _raising(r,_self): raise pymol.CmdException
        return r

    def get_mtl_obj(_self=cmd):
        '''
DESCRIPTION

    NOTE: this is an incomplete and unsupported feature.

    "get_mtl_obj" returns a tuple containing mtl and obj input files
    for use with Maya.

PYMOL API

    cmd.get_obj_mtl()

        '''
        r = DEFAULT_ERROR
        try:
            _self.lock(_self)   
            r = _cmd.get_mtl_obj(_self._COb)
        finally:
            _self.unlock(r,_self)
        if _raising(r,_self): raise pymol.CmdException
        return r
    
    def get_version(quiet=1,_self=cmd):
	'''
DESCRIPTION
 
    "get_version" returns a tuple of length three containing text,
    floating point, and integer representations of the current PyMOL
    version number.
   
PYMOL API

    cmd.get_version(int quiet)

	'''
        r = DEFAULT_ERROR
        try:
            _self.lock(_self)   
            r = _cmd.get_version(_self._COb)
        finally:
            _self.unlock(r,_self)
        if _raising(r,_self):
            raise pymol.CmdException
        else:
            if not quiet:
                if _feedback(fb_module.cmd,fb_mask.results,_self):
                    print " version: %s (%8.6f) %d"%r
        return r
    
    def get_vrml(version=2,_self=cmd): 
        '''
DESCRIPTION

    "get_vrml" returns a VRML2 string representing the content
    currently displayed.

PYMOL API

    cmd.get_vrml()

        '''
        r = DEFAULT_ERROR
        try:
            _self.lock(_self)   
            r = _cmd.get_vrml(_self._COb,int(version))
        finally:
            _self.unlock(r,_self)
        if _raising(r,_self): raise pymol.CmdException
        return r

    def count_states(selection="(all)", quiet=1, _self=cmd):
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
        r = DEFAULT_ERROR
        try:
            _self.lock(_self)
            r = _cmd.count_states(_self._COb,selection)
        finally:
            _self.unlock(r,_self)
        if is_ok(r):
            if not quiet:
                print " cmd.count_states: %d states."%r            
        if _raising(r,_self): raise pymol.CmdException
        return r
    
    def get_movie_length(quiet=1, images=-1, _self=cmd):
        '''
DESCRIPTION

    "get_movie_length" returns the number of frames explicitly defined
    in the movie, not including molecular states.

PYMOL API

    cmd.count_frames()

SEE ALSO

    frame, count_states, count_frames
    '''
        r = DEFAULT_ERROR
        try:
            _self.lock(_self)
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
                print " cmd.get_movie_length: %d frames"%r      
        finally:
            _self.unlock(r,_self)
        if _raising(r,_self): raise pymol.CmdException
        return r

    def count_frames(quiet=1, _self=cmd):
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
        r = DEFAULT_ERROR
        try:
            _self.lock(_self)
            r = _cmd.count_frames(_self._COb)
            if not quiet: print " cmd.count_frames: %d frames"%r      
        finally:
            _self.unlock(r,_self)
        if _raising(r,_self): raise pymol.CmdException
        return r

    def export_dots(object,state,_self=cmd):  
        '''
DESCRIPTION

    "export_dots" is an old unsupported command that may have
    something to do with returning the coordinates of the dot
    representation back to the Python layer.

    '''
        r = DEFAULT_ERROR
        try:
            _self.lock(_self)
            r = _cmd.export_dots(_self._COb,object,int(state)-1)
        finally:
            _self.unlock(r,_self)
        if _raising(r,_self): raise pymol.CmdException
        return r

    def overlap(selection1, selection2, state1=1, state2=1, adjust=0.0, quiet=1, _self=cmd):
        '''
DESCRIPTION

    "overlap" is an unsupported command that may have something to do
    with measuring the total amount of van der Waals overlap between
    two selections.

NOTES

    For best performance, use the smaller molecule as selection 1.
    
    '''
        
        # preprocess selections
        selection1 = selector.process(selection1)
        selection2 = selector.process(selection2)
        #
        r = DEFAULT_ERROR
        try:
            _self.lock(_self)
            r = _cmd.overlap(_self._COb,str(selection1),str(selection2),
                                  int(state1)-1,int(state2)-1,
                                  float(adjust))
            if not quiet: print " cmd.overlap: %5.3f Angstroms."%r
        finally:
            _self.unlock(r,_self)
        if _raising(r,_self): raise pymol.CmdException
        return r

    def get_movie_locked(_self=cmd):
        r = DEFAULT_ERROR
        try:
            _self.lock(_self)
            r = _cmd.get_movie_locked(_self._COb)
        finally:
            _self.unlock(r,_self)
        if _raising(r,_self): raise pymol.CmdException
        return r

    def get_object_color_index(name,_self=cmd):
        name = str(name)
        r = DEFAULT_ERROR
        try:
            _self.lock(_self)
            r = _cmd.get_object_color_index(_self._COb,name)
        finally:
            _self.unlock(r,_self)      
#        if _raising(r,_self): raise pymol.CmdException
        return r
    
    def get_color_tuple(name,mode=0,_self=cmd):
        name=str(name)
        r = DEFAULT_ERROR
        try:
            _self.lock(_self)
            r = _cmd.get_color(_self._COb,name,mode)
            if r==None:
                if _feedback(fb_module.cmd,fb_mask.errors,_self):         
                    print "cmd-Error: Unknown color '%s'."%name
        finally:
            _self.unlock(r,_self)
        if _raising(r,_self): raise pymol.CmdException
        return r

    def get_color_indices(all=0,_self=cmd):
        r = DEFAULT_ERROR
        try:
            _self.lock(_self)
            if all:
                r = _cmd.get_color(_self._COb,'',2)
            else:
                r = _cmd.get_color(_self._COb,'',1)            
        finally:
            _self.unlock(r,_self)
        if _raising(r,_self): raise pymol.CmdException
        return r

    def get_color_index(color,_self=cmd):
        r = DEFAULT_ERROR
        try:
            _self.lock(_self)
            r = _cmd.get_color(_self._COb,str(color),3)
        finally:
            _self.unlock(r,_self)
        if _raising(r,_self): raise pymol.CmdException
        return r
            
    def get_renderer(_self=cmd):  # 
        r = DEFAULT_ERROR
        try:
            _self.lock(_self)
            r = _cmd.get_renderer(_self._COb)
        finally:
            _self.unlock(r,_self)
        if _raising(r,_self): raise pymol.CmdException
        return r

    def get_phipsi(selection="(name ca)",state=-1,_self=cmd):
        # preprocess selections
        selection = selector.process(selection)
        #   
        r = DEFAULT_ERROR
        try:
            _self.lock(_self)
            r = _cmd.get_phipsi(_self._COb,"("+str(selection)+")",int(state)-1)
        finally:
            _self.unlock(r,_self)
        if _raising(r,_self): raise pymol.CmdException
        return r

    def get_atom_coords(selection, state=0,quiet=1,_self=cmd):
        '''
DESCRIPTION

    "get_atom_coords" returns the 3D coordinates of a single atom.
    
    '''
        # low performance way to get coords for a single atom
        r = []
        selection = selector.process(selection)
        try:
            _self.lock(_self)
            r = _cmd.get_atom_coords(_self._COb,str(selection),int(state)-1,int(quiet))
        finally:
            _self.unlock(r,_self)
        if r==None:
            if _self._raising(_self=_self): raise pymol.CmdException
        elif not quiet:
            for a in r:
                print " cmd.get_coords: [%8.3f,%8.3f,%8.3f]"%(a)
        if _raising(r,_self): raise pymol.CmdException
        return r

    
    def get_position(quiet=1, _self=cmd):
        '''
DESCRIPTION

    "get_position" returns the 3D coordinates of the center of the
    viewer window.

    '''
        
        try:
            _self.lock(_self)
            r = _cmd.get_position(_self._COb)
        finally:
            _self.unlock(r,_self)
        if _raising(r,_self):
            raise pymol.CmdException
        elif not quiet:
            print " cmd.get_position: [%8.3f,%8.3f,%8.3f]"%(r[0],r[1],r[2])
        return r

    def get_distance(atom1="pk1",atom2="pk2",state=-1,quiet=1,_self=cmd):
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
        r = DEFAULT_ERROR
        # preprocess selections
        atom1 = selector.process(atom1)
        atom2 = selector.process(atom2)
        #   
        r = None
        try:
            _self.lock(_self)
            r = _cmd.get_distance(_self._COb,str(atom1),str(atom2),int(state)-1)
        finally:
            _self.unlock(r,_self)
        if _raising(r,_self):
            raise pymol.CmdException
        elif not quiet:
            print " cmd.get_distance: %5.3f Angstroms."%r
        return r

    def get_angle(atom1="pk1",atom2="pk2",atom3="pk3",state=-1,quiet=1,_self=cmd):
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
        r = DEFAULT_ERROR
        try:
            _self.lock(_self)
            r = _cmd.get_angle(_self._COb,str(atom1),str(atom2),str(atom3),int(state)-1)
        finally:
            _self.unlock(r,_self)
        if _raising(r,_self):
            raise pymol.CmdException
        elif not quiet:
            print " cmd.get_angle: %5.3f degrees."%r
        return r
        
    def get_dihedral(atom1="pk1",atom2="pk2",atom3="pk3",atom4="pk4",state=-1,quiet=1,_self=cmd):
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
        r = DEFAULT_ERROR
        try:
            _self.lock(_self)
            r = _cmd.get_dihe(_self._COb,str(atom1),str(atom2),str(atom3),str(atom4),int(state)-1)
        finally:
            _self.unlock(r,_self)
        if _raising(r,_self):
            raise pymol.CmdException
        elif not quiet:
            print " cmd.get_dihedral: %5.3f degrees."%r
        return r

    def get_model(selection="(all)",state=1,ref='',ref_state=0,_self=cmd):
        '''
DESCRIPTION

    "get_model" returns a ChemPy "Indexed" format model from a selection.

PYMOL API

    cmd.get_model(string selection [,int state] )

        '''
        # preprocess selection
        selection = selector.process(selection)
        #   
        r = DEFAULT_ERROR
        try:
            _self.lock(_self)
            r = _cmd.get_model(_self._COb,"("+str(selection)+")",int(state)-1,str(ref),int(ref_state)-1)
            if r==None:
                r = DEFAULT_ERROR
        finally:
            _self.unlock(r,_self)
        if _raising(r,_self): raise pymol.CmdException
        return r

    def get_area(selection="(all)",state=1,load_b=0,quiet=1,_self=cmd):
        '''
        PRE-RELEASE functionality - API will change
        '''
        # preprocess selection
        selection = selector.process(selection)
        #
        r = DEFAULT_ERROR
        try:
            _self.lock(_self)
            r = _cmd.get_area(_self._COb,"("+str(selection)+")",int(state)-1,int(load_b))
        finally:
            _self.unlock(r,_self)
        if r<0.0: # negative area signals error condition
            if _self._raising(): raise pymol.CmdException
        elif not quiet:
            print " cmd.get_area: %5.3f Angstroms^2."%r
        return r

    def get_chains(selection="(all)",state=0,quiet=1,_self=cmd):
        '''
        PRE-RELEASE functionality - API will change

        state is currently ignored
        '''
        # preprocess selection
        selection = selector.process(selection)
        #
        r = DEFAULT_ERROR
        try:
            _self.lock(_self)
            r = _cmd.get_chains(_self._COb,"("+str(selection)+")",int(state)-1)
        finally:
            _self.unlock(r,_self)
        if r==None:
            return []
        if _raising(r,_self):
            raise pymol.CmdException
        elif not quiet:
            print " cmd.get_chains: ",str(r)
        return r


    def get_names(type='public_objects',enabled_only=0,selection="",_self=cmd):
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
        r = DEFAULT_ERROR
        # this needs to be replaced with a flag & masking scheme...
        mode = 1
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
            print "Error: unknown type: '%s'"%str(type)
            if _raising(-1,_self): raise pymol.CmdException            
        try:
            _self.lock(_self)
            r = _cmd.get_names(_self._COb,int(mode),int(enabled_only),str(selection))
        finally:
            _self.unlock(r,_self)
        if _raising(r,_self): raise pymol.CmdException
        return r

    def get_legal_name(name,_self=cmd):
        r = DEFAULT_ERROR
        try:
            _self.lock(_self)
            r = _cmd.get_legal_name(_self._COb,str(name))
        finally:
            _self.unlock(r,_self)
        if _raising(r,_self): raise pymol.CmdException
        return r

    def get_type(name,quiet=1,_self=cmd):
        '''
DESCRIPTION

    "get_type" returns a string describing the named object or
     selection or the string "nonexistent" if the name in unknown.

PYMOL API

    cmd.get_type(string object-name)

NOTES

    Possible return values are

    "object:molecule"
    "object:map"
    "object:mesh"
    "object:distance"
    "selection"

SEE ALSO

    get_names
        '''
        r = DEFAULT_ERROR
        try:
            _self.lock(_self)
            r = _cmd.get_type(_self._COb,str(name))
        finally:
            _self.unlock(r,_self)
        if is_error(r):
            if _feedback(fb_module.cmd,fb_mask.errors,_self):      
                print "cmd-Error: unrecognized name."
        elif not quiet:
            print r
        if _raising(r,_self): raise pymol.CmdException
        return r

    def id_atom(selection,mode=0,quiet=1,_self=cmd):
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
        l = apply(identify,(selection,mode,1))
        ll = len(l)
        if not ll:
            if _feedback(fb_module.cmd,fb_mask.errors,_self):
                print "cmd-Error: atom %s not found by id_atom." % selection
            if _self._raising(_self=_self): raise pymol.CmdException
        elif ll>1:
            if _feedback(fb_module.cmd,fb_mask.errors,_self):
                print "cmd-Error: multiple atoms %s found by id_atom." % selection
            if _self._raising(_self=_self): raise pymol.CmdException
        else:
            r = l[0]
            if not quiet:
                if mode:
                    print " cmd.id_atom: (%s and id %d)"%(r[0],r[1])
                else:
                    print " cmd.id_atom: (id %d)"%r
        if _raising(r,_self): raise pymol.CmdException
        return r

    def identify(selection="(all)",mode=0,quiet=1,_self=cmd):
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
        r = DEFAULT_ERROR
        try:
            _self.lock(_self)
            r = _cmd.identify(_self._COb,"("+str(selection)+")",int(mode)) # 0 = default mode
        finally:
            _self.unlock(r,_self)
        if is_list(r):
            if len(r):
                if not quiet:
                    if mode:
                        for a in r:
                            print " cmd.identify: (%s and id %d)"%(a[0],a[1])
                    else:
                        for a in r:
                            print " cmd.identify: (id %d)"%a
        if _raising(r,_self): raise pymol.CmdException
        return r

    def index(selection="(all)",quiet=1,_self=cmd):
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
        r = DEFAULT_ERROR
        try:
            _self.lock(_self)
            r = _cmd.index(_self._COb,"("+str(selection)+")",0) # 0 = default mode
        finally:
            _self.unlock(r,_self)
        if not quiet:
            if is_list(r):
                if len(r):
                    for a in r:
                        print " cmd.index: (%s`%d)"%(a[0],a[1])
        if _raising(r,_self): raise pymol.CmdException
        return r

    def find_pairs(selection1,selection2,state1=1,state2=1,cutoff=3.5,mode=0,angle=45,_self=cmd):
        '''
DESCRIPTION

    "find_pairs" is currently undocumented.

        '''
        # preprocess selection
        selection1 = selector.process(selection1)
        selection2 = selector.process(selection2)
        #      
        r = DEFAULT_ERROR
        try:
            _self.lock(_self)
            r = _cmd.find_pairs(_self._COb,"("+str(selection1)+")",
                                      "("+str(selection2)+")",
                                      int(state1)-1,int(state2)-1,
                                      int(mode),float(cutoff),float(angle))
            # 0 = default mode
        finally:
            _self.unlock(r,_self)
        if _raising(r,_self): raise pymol.CmdException
        return r

    def get_extent(selection="(all)",state=0,quiet=1,_self=cmd):
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
        r = DEFAULT_ERROR
        try:
            _self.lock(_self)
            r = _cmd.get_min_max(_self._COb,str(selection),int(state)-1)
        finally:
            _self.unlock(r,_self)
        if not r:
            if _self._raising(_self=_self): raise pymol.CmdException
        elif not quiet:
            print " cmd.extent: min: [%8.3f,%8.3f,%8.3f]"%(r[0][0],r[0][1],r[0][2])
            print " cmd.extent: max: [%8.3f,%8.3f,%8.3f]"%(r[1][0],r[1][1],r[1][2])      
        if _raising(r,_self): raise pymol.CmdException
        return r

    def phi_psi(selection="(byres pk1)", quiet=1, _self=cmd):
        '''
DESCRIPTION

    "phi_psi" return the phi and psi angles for a protein atom
    selection.
    
USAGE
        '''
        
        r = cmd.get_phipsi(selection)
        if r!=None:
            kees = r.keys()
            kees.sort()
            if not quiet:
                _self.feedback('push')
                _self.feedback('disable','executive','actions')
                for a in kees:
                    _self.iterate("(%s`%d)"%a,"print ' %-9s "+
                                ("( %6.1f, %6.1f )"%r[a])+
                                "'%(resn+'-'+resi+':')")
                _self.feedback('pop')
        elif _feedback(fb_module.cmd,fb_mask.errors,_self):      
            print "cmd-Error: can't compute phi_psi"
        if _raising(r,_self): raise pymol.CmdException
        return r


    def count_atoms(selection="(all)",quiet=1,state=0,domain='',_self=cmd):
        '''
DESCRIPTION

    "count_atoms" returns a count of atoms in a selection.

USAGE

    count_atoms (selection)

PYMOL API

    cmd.count_atoms(string selection)

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
        if not quiet: print " count_atoms: %d atoms"%r
        if _raising(r,_self): raise pymol.CmdException
        return r

    def get_names_of_type(type,_self=cmd):
        obj = _self.get_names('objects')
        types = map(get_type,obj)
        mix = map(None,obj,types)
        lst = []
        for a in mix:
            if a[1]==type:
                lst.append(a[0])
        return lst


    def get_raw_alignment(name='',active_only=0,_self=cmd):
        '''
DESCRIPTION

    "get_raw_alignment" returns a list of lists of (object,index)
    tuples containing the raw per-atom alignment relationships

PYMOL API

    cmd.get_raw_alignment(name)

        '''
        try:
            _self.lock(_self)
            r = _cmd.get_raw_alignment(_self._COb,str(name),int(active_only))
        finally:
            _self.unlock(r,_self)
        if _raising(r,_self): raise pymol.CmdException
        return r



