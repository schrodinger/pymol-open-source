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

if __name__=='pymol.experimenting':
    
    import selector
    from cmd import _cmd,lock,unlock,Shortcut,QuietException, \
          DEFAULT_ERROR, DEFAULT_SUCCESS, _raising, is_ok, is_error        
    import cmd
    import threading
    import pymol
    
    def get_bond_print(obj,max_bond,max_type):
        r = DEFAULT_ERROR
        try:
            lock()
            r = _cmd.get_bond_print(str(obj),int(max_bond),int(max_type))
        finally:
            unlock()
        if _raising(r): raise pymol.CmdException                  
        return r

    def expfit(a,b): # Huh?
        r = DEFAULT_ERROR
        try:
            lock()   
            r = _cmd.fit(a,b,2)
        finally:
            unlock()
        if _raising(r): raise pymol.CmdException                  
        return r

    def spheroid(object="",average=0):  # EXPERIMENTAL
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
            lock()
            r = _cmd.spheroid(str(object),int(average))
        finally:
            unlock()
        if _raising(r): raise pymol.CmdException                  
        return r

    def mem():
        '''
DESCRIPTION

    "mem" Dumps current memory state to standard output. This is a
    debugging feature, not an official part of the API.

    '''
        r = DEFAULT_ERROR
        try:
            lock()
            r = _cmd.mem()
        finally:
            unlock()
        if _raising(r): raise pymol.CmdException                  
        return r


    def check(selection=None,preserve=0):
    # UNSUPPORTED
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
    # OBSOLETE, TO BE REMOVED
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
    # OBSOLETE, TO BE REMOVED
        # NOTE: the realtime module relies on code that is not yet part of PyMOL/ChemPy
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


    def dump(fnam,obj):
        r = DEFAULT_ERROR
        try:
            lock()
            r = _cmd.dump(str(fnam),obj)
        finally:
            unlock()
        if _raising(r): raise pymol.CmdException                  
        return r


    def dummy(*arg):
        return None

    def test(group=0,index=0): # generic test routine for development
        r = DEFAULT_ERROR
        try:
            lock()   
            r=_cmd.test(int(group),int(index))
        finally:
            unlock()
        if _raising(r): raise pymol.CmdException                  
        return r

    def import_coords(obj,state,mechio): # experimental
        r = DEFAULT_ERROR      
        try:
            lock()   
            r = _cmd.import_coords(str(obj),int(state)-1,mechio)
        finally:
            unlock()
        if _raising(r): raise pymol.CmdException                  
        return r

    def load_coords(*arg): # UNSUPPORTED
        _self = cmd
        r = DEFAULT_ERROR
        try:
            lock()
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
                r = _cmd.load_coords(str(oname),model,
                                            int(state)-1,int(ftype))
            else:
                print "Error: invalid arguments."
        finally:
            unlock()
        if _raising(r): raise pymol.CmdException                  
        return r
