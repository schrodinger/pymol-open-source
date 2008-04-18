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

if __name__=='pymol.exporting':
    import os
    import thread
    import selector
    import string
    import re
    import copy
    
    import pymol
    import cmd
    from cmd import _cmd,lock,unlock,Shortcut,QuietException
    from chempy import io
    from chempy.sdf import SDF,SDFRec
    from cmd import _feedback,fb_module,fb_mask, \
                     DEFAULT_ERROR, DEFAULT_SUCCESS, _raising, is_ok, is_error, \
                     is_list, is_dict, is_tuple
    import traceback

    def copy_image(quiet=1,_self=cmd): # incentive feature / proprietary
        r = DEFAULT_ERROR
        if thread.get_ident() == pymol.glutThread:
            r = _self._copy_image(_self,int(quiet))
        else:
            r = _self.do('cmd._copy_image(quiet=%d,_self=cmd)'%int(quiet))
        if _self._raising(r,_self): raise QuietException         
        return r

    def get_pdbstr(selection="all", state=-1, ref='', ref_state=-1, quiet=1,_self=cmd):
        '''
DESCRIPTION

    "get_pdbstr" in an API-only function which returns a pdb
    corresponding to the atoms in the selection provided and that are
    present in the indicated state

PYMOL API

    cmd.get_pdbstr(string selection, int state)

NOTES

    "state" is a 1-based state index for the object.

    if state is -1, then current state is used.

    if state is 0, then all states are saved.
    
    '''
        r = DEFAULT_ERROR
        try:
            _self.lock(_self)   
            r = _cmd.get_pdb(_self._COb,str(selection),int(state)-1,0,str(ref),int(ref_state),int(quiet))
        finally:
            _self.unlock(r,_self)
        if _self._raising(r,_self): raise QuietException         
        return r

    def _get_dump_str(obj):
        if is_list(obj):
            list = map(_get_dump_str,obj)
            result = "[ "+string.join(list,",\n")+" ] "
        elif is_dict(obj):
            list = []
            for key in obj.keys():
                list.append( _get_dump_str(key)+" : "+_get_dump_str(obj[key]) )
            result = "{ "+string.join(list,",\n")+" } "
        elif is_tuple(obj):
            list = map(_get_dump_str,obj)
            result = "( "+string.join(list,",\n")+" ) "
        else:
            result = str(obj)
        return result
    
    def get_session(names='', partial=0, quiet=1, compress=-1, cache=-1, _self=cmd):
        session = {}
        r = DEFAULT_SUCCESS
        for a in _self._pymol._session_save_tasks:
            if a==None:
                try:
                    _self.lock(_self)
                    r = _cmd.get_session(_self._COb,session,str(names),
                                         int(partial),int(quiet))
                finally:
                    _self.unlock(r,_self)
                try:
                    session['session'] = copy.deepcopy(_self._pymol.session)
                    if cache<0:
                        cache = _self.get_setting_int('cache_mode')
                    if cache:
                        session['cache'] = _self._pymol._cache
                except:
                    traceback.print_exc()
            else:
                try:
                    if is_error(apply(a,(session,),{'_self':_self})):
                        r = DEFAULT_ERROR
                except:
                    traceback.print_exc()
                    print "Error: An error occurred when trying to generate session."
                    print "Error: The resulting session file may be incomplete."
        if is_ok(r):
            if(compress<0):
                compress = _self.get_setting_boolean('session_compression')
            if(compress):
                import zlib
                session = zlib.compress(io.pkl.toString(session))
            return session
        elif _self._raising(r,_self):
            raise QuietException                  
        return r
        

    def png(filename, width=0, height=0, dpi=-1.0, ray=0, quiet=1, _self=cmd):
        '''
DESCRIPTION

    "png" saves a PNG format image file of the current display.

USAGE

    png filename [, width [, height [, dpi [, ray]]]]

ARGUMENTS

    filename = string: file path to be written
    
    width = integer: width in pixels {default: 0 (current)}

    height = integer: height in pixels {default: 0 (current)}

    dpi = float: dots-per-inch {default -1.0 (unspecified)}

    ray = 0 or 1: should ray be run first {default: 0 (no)}

EXAMPLES

    png image.png
    png image.png, dpi=300

NOTES

    PNG is the only image format supported by PyMOL.

SEE ALSO

    mpng, save
    
PYMOL API

    cmd.png(string filename, int width, int width, float dpi,
            int ray, int quiet)
        '''
        r = DEFAULT_ERROR
        if thread.get_ident() == pymol.glutThread:
            r = _self._png(str(filename),int(width),int(height),float(dpi),int(ray),int(quiet),_self)
        else:
            r = _self._do("cmd._png('%s',%d,%d,%1.6f,%d,%d)"%(filename,width,height,dpi,ray,quiet),_self=_self)
        if _self._raising(r,_self): raise QuietException
        return r

    def export_coords(obj,state,_self=cmd): # experimental
        r = DEFAULT_ERROR
        try:
            _self.lock(_self)   
            r = _cmd.export_coords(_self._COb,str(obj),int(state)-1)
        finally:
            _self.unlock(r,_self)
        if _self._raising(r,_self): raise QuietException
        return r

    def multisave(filename, object, state=0, _self=cmd): # experimental -- deprecated
        '''
DESCRIPTION

    "multisave" is an unsupported command that may have something to
    do with writing multistate coordinate data to fragile version &
    platform-specific binary (pmo) file.
    
    '''
        
        r = DEFAULT_ERROR
        try:
            _self.lock(_self)
            r = _cmd.multisave(_self._COb,str(filename),str(object),int(state)-1,0)
        finally:
            _self.unlock(r,_self)
        if _self._raising(r,_self): raise QuietException
        return r

    def save(filename, selection='(all)', state=-1, format='', ref='',
             ref_state=-1, quiet=1, partial=0,_self=cmd):
        '''
DESCRIPTION

    "save" writes content to a file.
    
USAGE

    save filename [, selection [, state [, format ]]]

ARGUMENTS

    filename = string: file path to be written

    selection = string: atoms to save {default: (all)}

    state = integer: state to save {default: -1 (current state)}
    
PYMOL API

    cmd.save(string file, string selection, int state, string format)

NOTES

    The file format is automatically chosen if the extesion is one of
    the supported output formats: pdb, pqr, mol, pkl, mmd, mmod, pov,
    png, pse, aln, obj, mtl, or wrl.

    If the file format is not recognized, then a PDB file is written
    by default.

    For molecular files and where applicable and supported:
    
    * if state = -1 (default), then only the current state is written.

    * if state = 0, then a multi-state output file is written.
    
SEE ALSO

    load, get_model
        '''
        # preprocess selection
        input_selection = selection
        selection = selector.process(input_selection)
        #   
        r = DEFAULT_ERROR
        lc_filename=string.lower(filename)
        if format=='':
            format = 'unknown'
            if re.search("\.pdb$|\.ent$",lc_filename):
                format = 'pdb'
            elif re.search("\.pqr$",lc_filename):
                format = 'pqr'
            elif re.search("\.mol$",lc_filename):
                format = 'mol'
            elif re.search("\.sdf$",lc_filename):
                format = 'sdf'
            elif re.search("\.pkl$",lc_filename):
                format = 'pkl'
            elif re.search("\.pkl$",lc_filename):
                format = 'pkla'
            elif re.search("\.mmd$",lc_filename):
                format = 'mmod'
            elif re.search("\.mmod$",lc_filename):
                format = 'mmod'
            elif re.search("\.pmo$",lc_filename):
                format = 'pmo'
            elif re.search("\.pov$",lc_filename):
                format = 'pov'
            elif re.search("\.png$",lc_filename):
                format = 'png'
            elif re.search("\.pse$|\.psw$",lc_filename):
                format = 'pse'
            elif re.search("\.aln$",lc_filename):
                format = 'aln'
            elif re.search("\.obj$",lc_filename):
                format = 'obj'
            elif re.search("\.mtl$",lc_filename):
                format = 'mtl'
            elif re.search("\.wrl$",lc_filename):
                format = 'wrl'
            else:
                format = str(format)
        if format=='unknown':
            if not quiet:
                print " Save-Warning: Unrecognized file type -- defaulting to PDB format."
            format='pdb'
        filename = _self.exp_path(filename)
        if format=='pdb': # standard PDB file 
            f=open(filename,"w")
            if f:
                st = ''
                try:
                    _self.lock(_self)
                    st = _cmd.get_pdb(_self._COb,"("+str(selection)+")",int(state)-1,0,
                                      str(ref),int(ref_state)-1,int(quiet))
                    if st != None:
                        r = DEFAULT_SUCCESS                    
                finally:
                    _self.unlock(r,_self=_self)
                f.write(st)
                f.close()
                if not quiet:
                    print " Save: wrote \""+filename+"\"."
        elif format=='aln':
            st = ''
            try:
                _self.lock(_self)
                st = _cmd.get_seq_align_str(_self._COb,str(selection),int(state)-1,0,int(quiet))
                if st != None:
                    r = DEFAULT_SUCCESS
            finally:
                _self.unlock(r,_self)
            if st!=None:
                f=open(filename,"w")
                f.write(st)
                f.close()
                if not quiet:
                    print " Save: wrote \""+filename+"\"."
            else:
                r = DEFAULT_ERROR
        elif format=='pqr': # PQR (modified PDB file)
            f=open(filename,"w")
            if f:
                st = ''
                try:
                    _self.lock(_self)
                    st = _cmd.get_pdb(_self._COb,"("+str(selection)+")",int(state)-1,1,
                                      str(ref),int(ref_state)-1,int(quiet))
                    if st != None:
                        r = DEFAULT_SUCCESS
                finally:
                    _self.unlock(r,_self)
                f.write(st)
                f.close()
                if not quiet:
                    print " Save: wrote \""+filename+"\"."
        elif format=='pkl': # python binary
            io.pkl.toFile(_self.get_model(selection,state,ref,ref_state),filename)
            r = DEFAULT_SUCCESS
            if not quiet:
                print " Save: wrote \""+filename+"\"."
        elif format=='pkla': # ascii override
            io.pkl.toFile(_self.get_model(selection,state,ref,ref_state),filename,bin=0)
            r = DEFAULT_SUCCESS
            if not quiet:
                print " Save: wrote \""+filename+"\"."
        elif format=='pse': # PyMOL session
            _self.set("session_file",filename,quiet=1)
            if '(' in input_selection: # ignore selections 
                input_selection=''
            io.pkl.toFile(_self.get_session(str(input_selection),int(partial),int(quiet)),filename)
            r = DEFAULT_SUCCESS
            if not quiet:
                print " Save: wrote \""+filename+"\"."
        elif format=='mmod': # macromodel
            io.mmd.toFile(_self.get_model(selection,state,ref,ref_state),filename)
            r = DEFAULT_SUCCESS
            if not quiet:
                print " Save: wrote \""+filename+"\"."
        elif format=='sdf':
            state = int(state)
            if state==0:
                first_state = 1
                last_state = cmd.count_states(selection)
            else:
                first_state = state
                last_state = state
            sdf = SDF(filename,'w')
            for state in range(first_state, last_state+1):
                rec = SDFRec(io.mol.toList(_self.get_model(selection,state,ref,ref_state))
                             + ["$$$$\n"])
                sdf.write(rec)
            r = DEFAULT_SUCCESS
            if not quiet:
                print " Save: wrote %d states to \"%s\"."%(1+last_state-first_state,filename)
        elif format=='mol': 
            io.mol.toFile(_self.get_model(selection,state,ref,ref_state),filename)
            r = DEFAULT_SUCCESS
            if not quiet:
                print " Save: wrote \""+filename+"\"."
        elif format=='png':
            r = _self.png(filename,quiet=quiet)
        elif format=='pov':
            tup = _self.get_povray()
            f=open(filename,"w")
            f.write(tup[0])
            f.write(tup[1])
            f.flush()
            f.close()
            if not quiet:
                print " Save: wrote \""+filename+"\"."
            r = DEFAULT_SUCCESS
        elif format=='obj':
            tup = _self.get_mtl_obj()
            f=open(filename,"w")
            f.write(tup[1])
            f.flush()
            f.close()
            r = DEFAULT_SUCCESS
        elif format=='mtl':
            tup = _self.get_mtl_obj()
            f=open(filename,"w")
            f.write(tup[0])
            f.flush()
            f.close()
            r = DEFAULT_SUCCESS
        elif format=='wrl':
            txt = _self.get_vrml()
            f=open(filename,"w")
            f.write(txt)
            f.flush()
            f.close()
            if not quiet:
                print " Save: wrote \""+filename+"\"."
            r = DEFAULT_SUCCESS
        if _self._raising(r,_self): raise QuietException
        return r

