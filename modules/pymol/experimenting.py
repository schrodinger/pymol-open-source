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

if __name__=='pymol.experimenting':
    
    import selector
    from cmd import _cmd,lock,unlock,Shortcut,QuietException, \
          DEFAULT_ERROR, DEFAULT_SUCCESS, _raising, is_ok, is_error        
    import cmd
    import threading
    import pymol
    import string
    
    def get_bond_print(obj,max_bond,max_type,_self=cmd):
        r = DEFAULT_ERROR
        try:
            _self.lock(_self)
            r = _cmd.get_bond_print(_self._COb,str(obj),int(max_bond),int(max_type))
        finally:
            _self.unlock(r,_self)
        if _self._raising(r,_self): raise pymol.CmdException                  
        return r

    def expfit(a,b,_self=cmd): # Huh?
        r = DEFAULT_ERROR
        try:
            _self.lock(_self)   
            r = _cmd.fit(_self._COb,a,b,2)
        finally:
            _self.unlock(r,_self)
        if _self._raising(r,_self): raise pymol.CmdException                  
        return r

    def spheroid(object="",average=0,_self=cmd):  # EXPERIMENTAL
        '''
DESCRIPTION

    "spheroid" averages trajectory frames together to create
    an ellipsoid-like approximation of the actual anisotropic
    motion exhibited by the atom over a series of trajectory frames.

USAGE

    spheroid object,average

    average = number of states to average for each resulting spheroid state

    '''
        r = DEFAULT_ERROR
        try:
            print "Warning: 'spheroid' is experimental, incomplete, and unstable."
            _self.lock(_self)
            r = _cmd.spheroid(_self._COb,str(object),int(average))
        finally:
            _self.unlock(r,_self)
        if _self._raising(r,_self): raise pymol.CmdException                  
        return r

    def mem(_self=cmd):
        '''
DESCRIPTION

    "mem" Dumps current memory state to standard output. This is a
    debugging feature, not an official part of the API.

    '''
        r = DEFAULT_ERROR
        try:
            _self.lock(_self)
            r = _cmd.mem(_self._COb)
        finally:
            _self.unlock(r,_self)
        if _self._raising(r,_self): raise pymol.CmdException                  
        return r


    def check(selection=None, preserve=0):
        '''
DESCRIPTION

    "check" is unsupported command that may eventually have something
    to do with assigning forcefield parameters to a selection of
    atoms.
    
'''        
        # This function relies on code that is not currently part of PyMOL/ChemPy
        # NOTE: the realtime module relies on code that is not yet part of PyMOL/ChemPy
        from chempy.tinker import realtime
        if selection==None:
            arg = cmd.get_names("objects")
            arg = arg[0:1]
            if arg:
                if len(arg):
                    selection = arg
        if selection!=None:
            selection = selector.process(selection)
            realtime.assign("("+selection+")",int(preserve))
            realtime.setup("("+selection+")")

    def fast_minimize(*arg):
        '''
DESCRIPTION

    "fast_minimize" is an unsupported nonfunctional command that may
    eventually have something to do with doing a quick clean up of the
    molecular structure.
    
'''        
        from chempy.tinker import realtime  
        grad  = 0.01
        iter = 500
        interval = 50
        la = len(arg)
        if not la:
            arg = cmd.get_names("objects")
            arg = arg[0:1]
            la = len(arg)
        if la:
            sele  = "("+arg[0]+")"
            if la>1:
                iter = int(arg[1])
            if la>2:
                grad = float(arg[2])
            if la>3:
                interval = int(arg[3])
            t = threading.Thread(target=realtime.mini,args=(iter,grad,interval,arg[0]))
            t.setDaemon(1)
            t.start()

    def minimize(*arg):
        '''
DESCRIPTION

    "fast_minimize" is an unsupported nonfunctional command that may
    eventually have something to do with minimization.
    
'''        
        from chempy.tinker import realtime  
        grad  = 0.01
        iter = 500
        interval = 50
        la = len(arg)
        if not la:
            arg = cmd.get_names("objects")
            arg = arg[0:1]
            la = len(arg)
        if la:
            sele  = "("+arg[0]+")"
            if la>1:
                iter = int(arg[1])
            if la>2:
                grad = float(arg[2])
            if la>3:
                interval = int(arg[3])
            if realtime.setup(sele):
                t = threading.Thread(target=realtime.mini,args=(iter,grad,interval,arg[0]))
                t.setDaemon(1)
                t.start()
            else:
                print " minimize: missing parameters, can't continue"


    def dump(fnam,obj,_self=cmd):
        '''
DESCRIPTION

    "dump" is an unsupported command which may have something to do
    with outputing isosurface meshes and surface objects to a file.

    '''
        r = DEFAULT_ERROR
        try:
            _self.lock(_self)
            r = _cmd.dump(_self._COb,str(fnam),obj)
        finally:
            _self.unlock(r,_self)
        if _self._raising(r,_self): raise pymol.CmdException                  
        return r


    def dummy(*arg):
        return None

    def test(group=0,index=0,_self=cmd): # generic test routine for development
        '''
DESCRIPTION

    "dump" is an unsupported internal command.

    '''
        r = DEFAULT_ERROR
        try:
            _self.lock(_self)   
            r=_cmd.test(_self._COb,int(group),int(index))
        finally:
            _self.unlock(r,_self)
        if _self._raising(r,_self): raise pymol.CmdException                  
        return r

    def import_coords(coords,name,state,_self=cmd): # experimental
        r = DEFAULT_ERROR      
        try:
            _self.lock(_self)   
            r = _cmd.import_coords(_self._COb,str(name),int(state)-1,coords)
        finally:
            _self.unlock(r,_self)
        if _self._raising(r,_self): raise pymol.CmdException                  
        return r

    def load_coords(*arg): # UNSUPPORTED
        _self = cmd
        r = DEFAULT_ERROR
        try:
            _self.lock(_self)
            ok = 1
            ftype = loadable.model
            state = 0
            model = arg[0];
            if len(arg)<2:
                ok=0
            if len(arg)>=2:
                oname = string.strip(arg[1])
            if len(arg)>=3:
                state = int(arg[2])-1
            if ok:
                r = _cmd.load_coords(_self._COb,str(oname),model,
                                            int(state)-1,int(ftype))
            else:
                print "Error: invalid arguments."
        finally:
            _self.unlock(r,_self)
        if _self._raising(r,_self): raise pymol.CmdException                  
        return r
