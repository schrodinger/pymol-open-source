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
    from cmd import _feedback,fb_module,fb_mask, \
                     DEFAULT_ERROR, DEFAULT_SUCCESS, _raising, is_ok, is_error
    import traceback

    def get_pdbstr(selection="all", state=-1, ref='', ref_state=-1, quiet=1):
        '''
DESCRIPTION

    "get_pdbstr" in an API-only function which returns a pdb
    corresponding to the atoms in the selection provided and that are
    present in the indicated state

PYMOL API ONLY

    cmd.get_pdbstr( string selection="all", int state=0 )

NOTES

    "state" is a 1-based state index for the object.

    if state is -1, then current state is used.

    if state is 0, then all states are saved.
    
    '''
        r = DEFAULT_ERROR
        try:
            lock()   
            r = _cmd.get_pdb(str(selection),int(state)-1,0,str(ref),int(ref_state),int(quiet))
        finally:
            unlock(r)
        if _raising(r): raise QuietException         
        return r
    
    def get_session(names='', partial=0, quiet=1, compress=-1):
        session = {}
        r = DEFAULT_SUCCESS
        for a in pymol._session_save_tasks:
            if a==None:
                try:
                    lock()
                    r = _cmd.get_session(session,str(names),
                                         int(partial),int(quiet))
                finally:
                    unlock(r)
                try:
                    session['session'] = copy.deepcopy(pymol.session)
                except:
                    traceback.print_exc()
            else:
                try:
                    if is_error(apply(a,(session,))):
                        r = DEFAULT_ERROR
                except:
                    traceback.print_exc()
                    print "Error: An error occurred when trying to generate session."
                    print "Error: The resulting session file may be incomplete."
        if is_ok(r):
            if(compress<0):
                compress = cmd.get_setting_boolean('session_compression')
            if(compress):
                import zlib
                session = zlib.compress(io.pkl.toString(session))
            return session
        elif _raising(r):
            raise QuietException                  
        return r
        

    def png(filename, width=0, height=0, dpi=-1.0, ray=0, quiet=1):
        '''
DESCRIPTION

    "png" saves a PNG format image file of the current display.

USAGE

    png filename [, width [, height [, dpi [, ray]]]]

ARGUMENTS

    filename is the file to create
    
    width of image in pixels (default 0 -> viewport width)

    height of image in pixels (default 0 -> viewport height)

    dpi in dots-per-inch (default -1.0 -- unspecified)

    ray = 0 or 1 (default 0): should ray be run first?

EXAMPLES

    png image.png
    png image.png, dpi=300
       
PYMOL API

    cmd.png(string file, int width, int width, float dpi,
            int ray, int quiet)
        '''
        r = DEFAULT_ERROR
        if thread.get_ident() ==pymol.glutThread:
            r = cmd._png(str(filename),int(width),int(height),float(dpi),int(ray),int(quiet))
        else:
            r = cmd._do("cmd._png('%s',%d,%d,%1.6f,%d,%d)"%(filename,width,height,dpi,ray,quiet))
        if _raising(r): raise QuietException
        return r

    def export_coords(obj,state): # experimental
        r = DEFAULT_ERROR
        try:
            lock()   
            r = _cmd.export_coords(str(obj),int(state)-1)
        finally:
            unlock(r)
        if _raising(r): raise QuietException
        return r

    def multisave(filename,object,state=0): # experimental -- deprecated
        r = DEFAULT_ERROR
        try:
            lock()
            r = _cmd.multisave(str(filename),str(object),int(state)-1,0)
        finally:
            unlock(r)
        if _raising(r): raise QuietException
        return r

    def save(filename, selection='(all)', state=-1, format='', ref='',
             ref_state=-1, quiet=1, partial=0):
        '''
DESCRIPTION

    "save" writes content to a file.
    
USAGE

    save file [,selection [,state [,format ]]]

PYMOL API

    cmd.save(string file, string selection, int state, string format)

NOTES

    The file format is autodetected if the extesion is one of the
    supported output formats: pdb, pqr, mol, pkl, mmd, mmod, pov, png,
    pse, aln, obj, mtl, or wrl.

    Currently, if the file format is not recognized, then a PDB file
    is written by default

    For molecular files and where applicable:
    
    when state = -1 (default) only the current state is written.
    when state = 0, then a multi-state output PDB file is written.
    
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
        filename = cmd.exp_path(filename)
        if format=='pdb': # standard PDB file 
            f=open(filename,"w")
            if f:
                st = ''
                try:
                    lock()
                    st = _cmd.get_pdb("("+str(selection)+")",int(state)-1,0,
                                      str(ref),int(ref_state)-1,int(quiet))
                finally:
                    unlock()
                f.write(st)
                f.close()
                r = DEFAULT_SUCCESS
                if not quiet:
                    print " Save: wrote \""+filename+"\"."
        elif format=='aln':
            st = ''
            try:
                lock()
                st = _cmd.get_seq_align_str(str(selection),int(state)-1,0,int(quiet))
            finally:
                unlock()
            if st!=None:
                f=open(filename,"w")
                f.write(st)
                f.close()
                r = DEFAULT_SUCCESS
                if not quiet:
                    print " Save: wrote \""+filename+"\"."
            else:
                r = DEFAULT_ERROR
        elif format=='pqr': # PQR (modified PDB file)
            f=open(filename,"w")
            if f:
                st = ''
                try:
                    lock()
                    st = _cmd.get_pdb("("+str(selection)+")",int(state)-1,1,
                                      str(ref),int(ref_state)-1,int(quiet))
                finally:
                    unlock()
                f.write(st)
                f.close()
                r = DEFAULT_SUCCESS
                if not quiet:
                    print " Save: wrote \""+filename+"\"."
        elif format=='pkl': # python binary
            io.pkl.toFile(cmd.get_model(selection,state,ref,ref_state),filename)
            r = DEFAULT_SUCCESS
            if not quiet:
                print " Save: wrote \""+filename+"\"."
        elif format=='pkla': # ascii override
            io.pkl.toFile(cmd.get_model(selection,state,ref,ref_state),filename,bin=0)
            r = DEFAULT_SUCCESS
            if not quiet:
                print " Save: wrote \""+filename+"\"."
        elif format=='pse': # PyMOL session
            cmd.set("session_file",filename,quiet=1)
            if '(' in input_selection: # ignore selections 
                input_selection=''
            io.pkl.toFile(cmd.get_session(str(input_selection),int(partial),int(quiet)),filename)
            r = DEFAULT_SUCCESS
            if not quiet:
                print " Save: wrote \""+filename+"\"."
        elif format=='mmod': # macromodel
            io.mmd.toFile(cmd.get_model(selection,state,ref,ref_state),filename)
            r = DEFAULT_SUCCESS
            if not quiet:
                print " Save: wrote \""+filename+"\"."
        elif format=='mol': 
            io.mol.toFile(cmd.get_model(selection,state,ref,ref_state),filename)
            r = DEFAULT_SUCCESS
            if not quiet:
                print " Save: wrote \""+filename+"\"."
        elif format=='png':
            r = cmd.png(filename,quiet=quiet)
        elif format=='pov':
            tup = cmd.get_povray()
            f=open(filename,"w")
            f.write(tup[0])
            f.write(tup[1])
            f.flush()
            f.close()
            if not quiet:
                print " Save: wrote \""+filename+"\"."
            r = DEFAULT_SUCCESS
        elif format=='obj':
            tup = cmd.get_mtl_obj()
            f=open(filename,"w")
            f.write(tup[1])
            f.flush()
            f.close()
            r = DEFAULT_SUCCESS
        elif format=='mtl':
            tup = cmd.get_mtl_obj()
            f=open(filename,"w")
            f.write(tup[0])
            f.flush()
            f.close()
            r = DEFAULT_SUCCESS
        elif format=='wrl':
            txt = cmd.get_vrml()
            f=open(filename,"w")
            f.write(txt)
            f.flush()
            f.close()
            if not quiet:
                print " Save: wrote \""+filename+"\"."
            r = DEFAULT_SUCCESS
        if _raising(r): raise QuietException
        return r

