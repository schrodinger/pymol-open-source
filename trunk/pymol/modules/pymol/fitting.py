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

if __name__=='pymol.fitting':
    
    import cmd
    from cmd import _cmd,lock,unlock
    import selector
    import os
    import pymol
    
    from cmd import _cmd,lock,unlock,Shortcut, \
          DEFAULT_ERROR, DEFAULT_SUCCESS, _raising, is_ok, is_error
    
    def align(source, target, cutoff=2.0, cycles=2, gap=-10.0,
              extend=-0.5, max_gap=-1, object=None,
              matrix="BLOSUM62", source_state=0, target_state=0,
              quiet=1, max_skip=0, transform=1,reset=0):
        
        '''
DESCRIPTION

    "align" performs a sequence alignment followed by a structural
    alignment, and then carries out zero or more cycles of refinement
    in order to reject structural outliers found during the fit.

USAGE 

    align source, target [, object=alignment-object-name ]

NOTES

    If object is specified, then align will create an object which
    indicates paired atoms and supports visualization of the alignment
    in the sequence viewer.

EXAMPLE

    align protA////CA, protB////CA, object=alnAB

SEE ALSO

    pair_fit, fit, rms, rms_cur, intra_rms, intra_rms_cur
        '''
        r = DEFAULT_ERROR
        source = selector.process(source)
        target = selector.process(target)
        mfile = cmd.exp_path("$PYMOL_PATH/data/pymol/matrices/"+matrix)
        if object==None: object=''
        # delete existing alignment object (if asked to reset it)
        try:
            lock()
            r = _cmd.align(source,"("+target+")",float(cutoff),int(cycles),float(gap),
                                float(extend),int(max_gap),str(object),str(mfile),
                                int(source_state)-1,int(target_state)-1,
                                int(quiet),int(max_skip),int(transform),int(reset))
        finally:
            unlock(r)
        if _raising(r): raise pymol.CmdException         
        return r

    def intra_fit(selection,state=1,quiet=1,mix=0):
        '''
DESCRIPTION

    "intra_fit" fits all states of an object to an atom selection
    in the specified state.  It returns the rms values to python
    as an array.

USAGE 

    intra_fit (selection),state

PYMOL API

    cmd.intra_fit( string selection, int state )

EXAMPLES

    intra_fit ( name ca )

PYTHON EXAMPLE

    from pymol import cmd
    rms = cmd.intra_fit("(name ca)",1)

SEE ALSO

    fit, rms, rms_cur, intra_rms, intra_rms_cur, pair_fit
        '''
        # preprocess selection
        selection = selector.process(selection)
        #   
        r = DEFAULT_ERROR
        state = int(state)
        mix = int(mix)
        try:
            lock()
            r = _cmd.intrafit("("+str(selection)+")",int(state)-1,2,int(quiet),int(mix))
        finally:
            unlock(r)
        if r<0.0:
            r = DEFAULT_ERROR
        elif not quiet:
            st = 1
            for a in r:
                if a>=0.0:
                    if mix:
                        print " cmd.intra_fit: %5.3f in state %d vs mixed target"%(a,st)
                    else:
                        print " cmd.intra_fit: %5.3f in state %d vs state %d"%(a,st,state)
                st = st + 1
        if _raising(r): raise pymol.CmdException         
        return r

    def intra_rms(selection,state=0,quiet=1):
        '''
DESCRIPTION

    "intra_rms" calculates rms fit values for all states of an object
    over an atom selection relative to the indicated state.  
    Coordinates are left unchanged.  The rms values are returned
    as a python array.

PYMOL API

    cmd.intra_rms( string selection, int state)

PYTHON EXAMPLE

    from pymol import cmd
    rms = cmd.intra_rms("(name ca)",1)

SEE ALSO

    fit, rms, rms_cur, intra_fit, intra_rms_cur, pair_fit
        '''
        # preprocess selection
        selection = selector.process(selection)
        #   
        r = DEFAULT_ERROR
        state = int(state)
        try:
            lock()
            r = _cmd.intrafit("("+str(selection)+")",int(state)-1,1,int(quiet),int(0))
        finally:
            unlock(r)
        if r<0.0:
            r = DEFAULT_ERROR
        elif not quiet:
            st = 1
            for a in r:
                if a>=0.0:
                    print " cmd.intra_rms: %5.3f in state %d vs state %d"%(a,st,state)
                st = st + 1
        if _raising(r): raise pymol.CmdException         
        return r

    def intra_rms_cur(selection,state=0,quiet=1):
        '''
DESCRIPTION

    "intra_rms_cur" calculates rms values for all states of an object
    over an atom selection relative to the indicated state without
    performing any fitting.  The rms values are returned
    as a python array.

PYMOL API

    cmd.intra_rms_cur( string selection, int state)

PYTHON EXAMPLE

    from pymol import cmd
    rms = cmd.intra_rms_cur("(name ca)",1)

SEE ALSO

    fit, rms, rms_cur, intra_fit, intra_rms, pair_fit
        '''
        # preprocess selection
        selection = selector.process(selection)
        #   
        r = DEFAULT_ERROR
        state = int(state)
        try:
            lock()
            r = _cmd.intrafit("("+str(selection)+")",int(state)-1,0,int(quiet),int(0))
        finally:
            unlock(r)
        if r<0.0:
            r = DEFAULT_ERROR
        elif not quiet:
            st = 1
            for a in r:
                if a>=0.0:
                    print " cmd.intra_rms_cur: %5.3f in state %d vs state %d"%(a,st,state)
                st = st + 1
        if _raising(r): raise pymol.CmdException         
        return r

    def fit(selection,target,source_state=0,target_state=0,
              quiet=1,matchmaker=0,cutoff=2.0,cycles=0,object=None):
        '''
DESCRIPTION

    "fit" superimposes the model in the first selection on to the model
    in the second selection.  Only matching atoms in both selections
    will be used for the fit.

USAGE

    fit (selection), (target-selection)

EXAMPLES

    fit ( mutant and name ca ), ( wildtype and name ca )

SEE ALSO

    rms, rms_cur, intra_fit, intra_rms, intra_rms_cur
        '''
        r = DEFAULT_ERROR      
        a=str(selection)
        b=str(target)
        # preprocess selections
        a = selector.process(a)
        b = selector.process(b)
        #
        if object==None: object=''
        if matchmaker==0:
            sele1 = "((%s) in (%s))" % (str(a),str(b))
            sele2 = "((%s) in (%s))" % (str(b),str(a))
        else:
            sele1 = str(a)
            sele2 = str(b)
        try:
            lock()
            r = _cmd.fit(sele1,sele2,2,
                             int(source_state)-1,int(target_state)-1,
                             int(quiet),int(matchmaker),float(cutoff),
                             int(cycles),str(object))
        finally:
            unlock(r)
        if _raising(r): raise pymol.CmdException         
        return r

    def rms(selection,target,source_state=0,target_state=0,quiet=1,
              matchmaker=0,cutoff=2.0,cycles=0,object=None):
        '''
DESCRIPTION

    "rms" computes a RMS fit between two atom selections, but does not
    tranform the models after performing the fit.

USAGE

    rms (selection), (target-selection)

EXAMPLES

    fit ( mutant and name ca ), ( wildtype and name ca )

SEE ALSO

    fit, rms_cur, intra_fit, intra_rms, intra_rms_cur, pair_fit   
        '''
        r = DEFAULT_ERROR      
        a=str(selection)
        b=str(target)
        # preprocess selections
        a = selector.process(a)
        b = selector.process(b)
        #
        if object==None: object=''      
        if matchmaker==0:
            sele1 = "((%s) in (%s))" % (str(a),str(b))
            sele2 = "((%s) in (%s))" % (str(b),str(a))
        else:
            sele1 = str(a)
            sele2 = str(b)
        try:
            lock()   
            r = _cmd.fit(sele1,sele2,1,
                             int(source_state)-1,int(target_state)-1,
                             int(quiet),int(matchmaker),float(cutoff),
                             int(cycles),str(object))
        finally:
            unlock(r)
        if _raising(r): raise pymol.CmdException         
        return r

    def rms_cur(selection,target,source_state=0,target_state=0,
                    quiet=1,matchmaker=0,cutoff=2.0,cycles=0,object=None):
        '''
DESCRIPTION

    "rms_cur" computes the RMS difference between two atom
    selections without performing any fitting.

USAGE

    rms_cur (selection), (selection)

SEE ALSO

    fit, rms, intra_fit, intra_rms, intra_rms_cur, pair_fit   
        '''
        r = DEFAULT_ERROR      
        a=str(selection)
        b=str(target)
        # preprocess selections
        a = selector.process(a)
        b = selector.process(b)
        #
        if object==None: object=''            
        if matchmaker==0:
            sele1 = "((%s) in (%s))" % (str(a),str(b))
            sele2 = "((%s) in (%s))" % (str(b),str(a))
        else:
            sele1 = str(a)
            sele2 = str(b)
        try:
            lock()
            r = _cmd.fit(sele1,sele2,0,
                             int(source_state)-1,int(target_state)-1,
                             int(quiet),int(matchmaker),float(cutoff),
                             int(cycles),str(object))
        finally:
            unlock(r)
        if _raising(r): raise pymol.CmdException         
        return r

    def pair_fit(*arg):
        '''
DESCRIPTION

    "pair_fit" fits a set of atom pairs between two models.  Each atom
    in each pair must be specified individually, which can be tedious
    to enter manually.  Script files are recommended when using this
    command.

USAGE

    pair_fit (selection), (selection), [ (selection), (selection) [ ...] ]

SEE ALSO

    fit, rms, rms_cur, intra_fit, intra_rms, intra_rms_cur
        '''
        r = DEFAULT_ERROR      
        new_arg = []
        for a in arg:
            new_arg.append(selector.process(a))
        try:
            lock()   
            r = _cmd.fit_pairs(new_arg)
        finally:
            unlock(r)
        if _raising(r): raise pymol.CmdException         
        return r






