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

if True:

        cmd = __import__("sys").modules["pymol.cmd"]
        from .cmd import _cmd,lock,unlock
        from . import selector
        import os
        import pymol

        from .cmd import _cmd,lock,unlock,Shortcut, \
            DEFAULT_ERROR, DEFAULT_SUCCESS, _raising, is_ok, is_error


        def cealign(target, mobile, target_state=1, mobile_state=1, quiet=1,
                    guide=1, d0=3.0, d1=4.0, window=8, gap_max=30, transform=1,
                    object=None, *, _self=cmd):
                '''
DESCRIPTION

    "cealign" aligns two proteins using the CE algorithm.

USAGE 

    cealign target, mobile [, target_state [, mobile_state [, quiet [,
        guide [, d0 [, d1 [, window [, gap_max, [, transform ]]]]]]]]]

NOTES

    If "guide" is set PyMOL will align using only alpha carbons, which is the
    default behavior. Otherwise, PyMOL will use all atoms. If "quiet" is set
    to -1, PyMOL will print the rotation matrix as well.

    Reference: Shindyalov IN, Bourne PE (1998) Protein structure alignment by
    incremental combinatorial extension (CE) of the optimal path.  Protein
    Engineering 11(9) 739-747.

EXAMPLES

    cealign protA////CA, protB////CA

    # fetch two proteins and align them
    fetch 1rlw 1rsy, async=0
    cealign 1rlw, 1rsy

SEE ALSO

    align, pair_fit, fit, rms, rms_cur, intra_rms, intra_rms_cur, super
                '''
                quiet = int(quiet)
                window = int(window)
                guide = "" if int(guide)==0 else "and guide"

                mobile = "(%s) %s" % (mobile,guide)
                target = "(%s) %s" % (target,guide)

                # handle PyMOL's macro /// notation
                mobile = selector.process(mobile)
                target = selector.process(target)

                # make the lists for holding coordinates and IDs
                mod1 = _self.get_model(target, state=target_state)
                mod2 = _self.get_model(mobile, state=mobile_state)
                sel1 = mod1.get_coord_list()
                sel2 = mod2.get_coord_list()
                ids1 = [a.id for a in mod1.atom]
                ids2 = [a.id for a in mod2.atom]

                if len(sel1) < 2 * window:
                        print("CEalign-Error: Your target selection is too short.")
                        raise pymol.CmdException
                if len(sel2) < 2 * window:
                        print("CEalign-Error: Your mobile selection is too short.")
                        raise pymol.CmdException
                if window < 3:
                        print("CEalign-Error: window size must be an integer greater than 2.")
                        raise pymol.CmdException
                if int(gap_max) < 0:
                        print("CEalign-Error: gap_max must be a positive integer.")
                        raise pymol.CmdException

                r = DEFAULT_ERROR

                try:
                        _self.lock(_self)

                        # call the C function
                        r = _cmd.cealign( _self._COb, sel1, sel2, float(d0), float(d1), int(window), int(gap_max) )

                        (aliLen, RMSD, rotMat, i1, i2) = r
                        if quiet==-1:
                                import pprint
                                print("RMSD %f over %i residues" % (float(RMSD), int(aliLen)))
                                print("TTT Matrix:")
                                pprint.pprint(rotMat)
                        elif quiet==0:
                                print("RMSD %f over %i residues" % (float(RMSD), int(aliLen)))

                        if int(transform):
                            for model in _self.get_object_list("(" + mobile + ")"):
                                _self.transform_object(model, rotMat, state=0)

                        if object is not None:
                            obj1 = _self.get_object_list("(" + target + ")")
                            obj2 = _self.get_object_list("(" + mobile + ")")
                            if len(obj1) > 1 or len(obj2) > 1:
                                print(' CEalign-Error: selection spans multiple' + \
                                        ' objects, cannot create alignment object')
                                raise pymol.CmdException
                            tmp1 = _self.get_unused_name('_1')
                            tmp2 = _self.get_unused_name('_2')
                            _self.select_list(tmp1, obj1[0], [ids1[j] for i in i1 for j in range(i, i+window)])
                            _self.select_list(tmp2, obj2[0], [ids2[j] for i in i2 for j in range(i, i+window)])
                            _self.rms_cur(tmp2, tmp1, cycles=0, matchmaker=4, object=object)
                            _self.delete(tmp1)
                            _self.delete(tmp2)
                except SystemError:
                    # findBest might return NULL, which raises SystemError
                    print(" CEalign-Error: alignment failed")
                finally:
                        _self.unlock(r,_self)
                if _self._raising(r,_self): raise pymol.CmdException
                return ( {"alignment_length": aliLen, "RMSD" : RMSD, "rotation_matrix" : rotMat } )

        def extra_fit(selection='(all)', reference='', method='align', zoom=1,
                quiet=0, *, _self=cmd, **kwargs):
            '''
DESCRIPTION

    Like "intra_fit", but for multiple objects instead of
    multiple states.

ARGUMENTS

    selection = string: atom selection of multiple objects {default: all}

    reference = string: reference object name {default: first object in selection}

    method = string: alignment method (command that takes "mobile" and "target"
    arguments, like "align", "super", "cealign" {default: align}

    ... extra arguments are passed to "method"

SEE ALSO

    align, super, cealign, intra_fit, util.mass_align
            '''
            zoom, quiet = int(zoom), int(quiet)
            sele_name = _self.get_unused_name('_')
            _self.select(sele_name, selection, 0)
            models = _self.get_object_list(sele_name)

            if not reference:
                reference = models[0]
                models = models[1:]
            elif reference in models:
                models.remove(reference)
            else:
                _self.select(sele_name, reference, merge=1)

            if _self.is_string(method):
                if method in _self.keyword:
                    method = _self.keyword[method][0]
                else:
                    raise pymol.CmdException(method, 'Unknown method')

            for model in models:
                x = method(mobile='?%s & ?%s' % (sele_name, model),
                        target='?%s & ?%s' % (sele_name, reference), **kwargs)
                if not quiet:
                    if _self.is_sequence(x):
                        print('%-20s RMSD = %8.3f (%d atoms)' % (model, x[0], x[1]))
                    elif isinstance(x, float):
                        print('%-20s RMSD = %8.3f' % (model, x))
                    elif isinstance(x, dict) and 'RMSD' in x:
                        natoms = x.get('alignment_length', 0)
                        suffix = (' (%s atoms)' % natoms) if natoms else ''
                        print('%-20s RMSD = %8.3f' % (model, x['RMSD']) + suffix)
                    else:
                        print('%-20s' % (model,))

            if zoom:
                _self.zoom(sele_name)
            _self.delete(sele_name)


        def alignto(target='', method="cealign", selection='', quiet=1, *, _self=cmd, **kwargs):
                """
DESCRIPTION

    "alignto" aligns all other loaded objects to the target
    using the specified alignment algorithm.

    This is a wrapper for "extra_fit".

USAGE

    alignto [ target [, method [, selection [, quiet ]]]]

ARGUMENTS

    target = str: reference object name {default: first object in selection}

    method = str: (see "extra_fit") {default: cealign}

    selection = str: (see "extra_fit") {default: public objects}

EXAMPLE

    # fetch some calmodulins
    fetch 1cll 1sra 1ggz 1k95, async=0

    # align them to 1cll using cealign
    alignto 1cll, method=cealign
    alignto 1cll, object=all_to_1cll

SEE ALSO

    extra_fit, align, super, cealign, fit, rms, rms_cur, intra_fit
                """
                if not selection:
                    names = _self.get_names("public_objects", 1)
                    if not names:
                        raise pymol.CmdException('no public objects')
                    selection = '%' + ' %'.join(names)
                return extra_fit(selection, target, method, 0, quiet, _self=_self, **kwargs)



        def super(mobile, target, cutoff=2.0, cycles=5,
                          gap=-1.5, extend=-0.7, max_gap=50, object=None,
                          matrix="BLOSUM62", mobile_state=0, target_state=0,
                          quiet=1, max_skip=0, transform=1, reset=0,
                          seq=0.0, radius=12.0, scale=17.0, base=0.65,
                          coord=0.0, expect=6.0, window=3, ante=-1.0,
                          *, _self=cmd):

                '''
DESCRIPTION

    "super" performs a residue-based pairwise alignment followed by a
    structural superposition, and then carries out zero or more cycles
    of refinement in order to reject outliers.

USAGE 

    super mobile, target [, object=name ]

NOTES

    By adjusting various parameters, the nature of the initial
    alignment can be modified to include or exclude various factors
    including sequence similarity, main chain path, secondary &
    tertiary structure, and current coordinates.

EXAMPLE

    super protA////CA, protB////CA, object=supeAB

SEE ALSO

    align, pair_fit, fit, rms, rms_cur, intra_rms, intra_rms_cur
        '''
                r = DEFAULT_ERROR
                mobile = selector.process(mobile)
                target = selector.process(target)
                if object is None: object=''
                matrix = str(matrix)
                if matrix.lower() in ['none', '']:
                        mfile = ''
                elif os.path.exists(matrix):
                        mfile = matrix
                else:
                        mfile = cmd.exp_path("$PYMOL_DATA/pymol/matrices/"+matrix)
                # delete existing alignment object (if asked to reset it)
                try:
                        _self.lock(_self)

                        r = _cmd.align(_self._COb,mobile,"("+target+")",float(cutoff),
                                                   int(cycles),float(gap),float(extend),int(max_gap),
                                                   str(object),str(mfile),
                                                   int(mobile_state)-1,int(target_state)-1,
                                                   int(quiet),int(max_skip),int(transform),
                                                   int(reset),float(seq),
                                                   float(radius),float(scale),float(base),
                                                   float(coord),float(expect),int(window),
                                                   float(ante))

                finally:
                        _self.unlock(r,_self)
                if _self._raising(r,_self): raise pymol.CmdException
                return r

        def align(mobile, target, cutoff=2.0, cycles=5, gap=-10.0,
                          extend=-0.5, max_gap=50, object=None,
                          matrix="BLOSUM62", mobile_state=0, target_state=0,
                          quiet=1, max_skip=0, transform=1, reset=0, *, _self=cmd):

                '''
DESCRIPTION

    "align" performs a sequence alignment followed by a structural
    superposition, and then carries out zero or more cycles of
    refinement in order to reject structural outliers found during
    the fit. "align" does a good job on proteins with decent sequence
    similarity (identity >30%). For comparing proteins with lower
    sequence identity, the "super" and "cealign" commands perform
    better.

USAGE 

    align mobile, target [, cutoff [, cycles
        [, gap [, extend [, max_gap [, object
        [, matrix [, mobile_state [, target_state
        [, quiet [, max_skip [, transform [, reset ]]]]]]]]]]]]]

ARGUMENTS

    mobile = string: atom selection of mobile object

    target = string: atom selection of target object

    cutoff = float: outlier rejection cutoff in Angstrom {default: 2.0}

    cycles = int: maximum number of outlier rejection cycles {default: 5}

    gap, extend, max_gap: sequence alignment parameters

    object = string: name of alignment object to create {default: (no
    alignment object)}

    matrix = string: file name of substitution matrix for sequence
    alignment {default: BLOSUM62}

    mobile_state = int: object state of mobile selection {default: 0 = all states}

    target_state = int: object state of target selection {default: 0 = all states}

    transform = 0/1: do superposition {default: 1}
        
NOTES

    If object is specified, then align will create an object which
    indicates paired atoms and supports visualization of the alignment
    in the sequence viewer.

    The RMSD of the aligned atoms (after outlier rejection!) is reported
    in the text output. The all-atom RMSD can be obtained by setting
    cycles=0 and thus not doing any outlier rejection.

EXAMPLE

    align protA////CA, protB////CA, object=alnAB

SEE ALSO

    super, cealign, pair_fit, fit, rms, rms_cur, intra_rms, intra_rms_cur
                '''
                r = DEFAULT_ERROR
                mobile = selector.process(mobile)
                target = selector.process(target)
                matrix = str(matrix)
                if matrix.lower() in ['none', '']:
                        mfile = ''
                elif os.path.exists(matrix):
                        mfile = matrix
                else:
                        mfile = cmd.exp_path("$PYMOL_DATA/pymol/matrices/"+matrix)
                if object is None: object=''
                # delete existing alignment object (if asked to reset it)
                try:
                        _self.lock(_self)
                        r = _cmd.align(_self._COb,mobile,"("+target+")",
                                                   float(cutoff),int(cycles),float(gap),
                                                   float(extend),int(max_gap),str(object),str(mfile),
                                                   int(mobile_state)-1,int(target_state)-1,
                                                   int(quiet),int(max_skip),int(transform),int(reset),
                                                   -1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0, 0.0)
                finally:
                        _self.unlock(r,_self)
                if _self._raising(r,_self): raise pymol.CmdException
                return r

        def intra_fit(selection, state=1, quiet=1, mix=0, *, pbc=1, _self=cmd):
                '''
DESCRIPTION

    "intra_fit" fits all states of an object to an atom selection
    in the specified state.  It returns the rms values to python
    as an array.

USAGE 

    intra_fit selection [, state]

ARGUMENTS

    selection = string: atoms to fit

    state = integer: target state

    pbc = 0/1: Consider periodic boundary conditions {default: 1}

PYMOL API

    cmd.intra_fit( string selection, int state )

EXAMPLES

    intra_fit ( name CA )

PYTHON EXAMPLE

    from pymol import cmd
    rms = cmd.intra_fit("(name CA)",1)

SEE ALSO

    fit, rms, rms_cur, intra_rms, intra_rms_cur, pair_fit
                '''
                # preprocess selection
                selection = selector.process(selection)
                #
                r = DEFAULT_ERROR
                state = int(state)
                mix = int(mix)
                with _self.lockcm:
                    r = _cmd.intrafit(_self._COb, selection,
                              int(state) - 1, 2, int(quiet), int(mix), int(pbc))
                if not isinstance(r, list):
                        r = DEFAULT_ERROR
                elif not quiet:
                        st = 1
                        for a in r:
                                if a>=0.0:
                                        if mix:
                                                print(" cmd.intra_fit: %5.3f in state %d vs mixed target"%(a,st))
                                        else:
                                                print(" cmd.intra_fit: %5.3f in state %d vs state %d"%(a,st,state))
                                st = st + 1
                if _self._raising(r,_self): raise pymol.CmdException
                return r

        def intra_rms(selection, state=0, quiet=1, *, _self=cmd):
                '''
DESCRIPTION

    "intra_rms" calculates rms fit values for all states of an object
    over an atom selection relative to the indicated state.
    Coordinates are left unchanged.  The rms values are returned as a
    python array.

EXAMPLE

    from pymol import cmd
    rms = cmd.intra_rms("(name CA)",1)
    print rms

PYMOL API

    cmd.intra_rms(string selection, int state)

SEE ALSO

    fit, rms, rms_cur, intra_fit, intra_rms_cur, pair_fit
                '''
                # preprocess selection
                selection = selector.process(selection)
                #
                r = DEFAULT_ERROR
                state = int(state)
                try:
                        _self.lock(_self)
                        r = _cmd.intrafit(_self._COb,"("+str(selection)+")",int(state)-1,1,int(quiet),int(0))
                finally:
                        _self.unlock(r,_self)
                if not isinstance(r, list):
                        r = DEFAULT_ERROR
                elif not quiet:
                        st = 1
                        for a in r:
                                if a>=0.0:
                                        print(" cmd.intra_rms: %5.3f in state %d vs state %d"%(a,st,state))
                                st = st + 1
                if _self._raising(r,_self): raise pymol.CmdException
                return r

        def intra_rms_cur(selection, state=0, quiet=1, *, _self=cmd):
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
    rms = cmd.intra_rms_cur("(name CA)",1)

SEE ALSO

    fit, rms, rms_cur, intra_fit, intra_rms, pair_fit
                '''
                # preprocess selection
                selection = selector.process(selection)
                #
                r = DEFAULT_ERROR
                state = int(state)
                try:
                        _self.lock(_self)
                        r = _cmd.intrafit(_self._COb,"("+str(selection)+")",int(state)-1,0,int(quiet),int(0))
                finally:
                        _self.unlock(r,_self)
                if not isinstance(r, list):
                    r = DEFAULT_ERROR
                elif not quiet:
                    st = 1
                    for a in r:
                        if a>=0.0:
                            print(" cmd.intra_rms_cur: %5.3f in state %d vs state %d"%(a,st,state))
                        st = st + 1
                if _self._raising(r,_self): raise pymol.CmdException
                return r

        def fit(mobile, target, mobile_state=0, target_state=0,
		quiet=1, matchmaker=0, cutoff=2.0, cycles=0, object=None, *, _self=cmd):
            '''
DESCRIPTION

    "fit" superimposes the model in the first selection on to the model
    in the second selection. Only matching atoms in both selections will
    be used for the fit.
	
USAGE

    fit mobile, target [, mobile_state [, target_state [, quiet
        [, matchmaker [, cutoff [, cycles [, object ]]]]]]]

ARGUMENTS

    mobile = string: atom selection

    target = string: atom selection

    mobile_state = integer: object state {default=0, all states)

    target_state = integer: object state {default=0, all states)

    matchmaker = integer: how to match atom pairs {default: 0}
        -1:  assume that atoms are stored in the identical order
        0/1: match based on all atom identifiers (segi,chain,resn,resi,name,alt)
        2:   match based on ID
        3:   match based on rank
        4:   match based on index (same as -1 ?)

    cutoff = float: outlier rejection cutoff (only if cycles>0) {default: 2.0}

    cycles = integer: number of cycles in outlier rejection refinement {default: 0}

    object = string: name of alignment object to create {default: None}

EXAMPLES

	fit protA, protB

NOTES

	Since atoms are matched based on all of their identifiers
	(including segment and chain identifiers), this command is only
	helpful when comparing very similar structures.

SEE ALSO

	align, super, pair_fit, rms, rms_cur, intra_fit, intra_rms, intra_rms_cur
            '''
            r = DEFAULT_ERROR
            a=str(mobile)
            b=str(target)
            # preprocess selections
            a = selector.process(a)
            b = selector.process(b)
            #
            if object is None: object=''
            if int(matchmaker)==0:
                sele1 = "((%s) in (%s))" % (str(a),str(b))
                sele2 = "((%s) in (%s))" % (str(b),str(a))
            else:
                sele1 = str(a)
                sele2 = str(b)
            try:
                _self.lock(_self)
                r = _cmd.fit(_self._COb,sele1,sele2,2,
                        int(mobile_state)-1,int(target_state)-1,
                        int(quiet),int(matchmaker),float(cutoff),
                        int(cycles),str(object))
            finally:
                _self.unlock(r,_self)
            if r < -0.5:
                raise pymol.CmdException
            return r

        def rms(mobile, target, mobile_state=0, target_state=0, quiet=1,
			  matchmaker=0, cutoff=2.0, cycles=0, object=None, *, _self=cmd):
            '''
DESCRIPTION

    "rms" computes a RMS fit between two atom selections, but does not
    tranform the models after performing the fit.

USAGE

    rms (selection), (target-selection)

EXAMPLES

    fit ( mutant and name CA ), ( wildtype and name CA )

SEE ALSO

    fit, rms_cur, intra_fit, intra_rms, intra_rms_cur, pair_fit
            '''
            r = DEFAULT_ERROR
            a=str(mobile)
            b=str(target)
            # preprocess selections
            a = selector.process(a)
            b = selector.process(b)
            #
            if object is None: object=''
            if int(matchmaker)==0:
                sele1 = "((%s) in (%s))" % (str(a),str(b))
                sele2 = "((%s) in (%s))" % (str(b),str(a))
            else:
                sele1 = str(a)
                sele2 = str(b)
            try:
                _self.lock(_self)
                r = _cmd.fit(_self._COb,sele1,sele2,1,
                        int(mobile_state)-1,int(target_state)-1,
                        int(quiet),int(matchmaker),float(cutoff),
                        int(cycles),str(object))
            finally:
                _self.unlock(r,_self)
            if r < -0.5:
                raise pymol.CmdException
            return r

        def rms_cur(mobile, target, mobile_state=0, target_state=0,
				quiet=1, matchmaker=0, cutoff=2.0, cycles=0,
				object=None, *, _self=cmd):

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
            a=str(mobile)
            b=str(target)
            # preprocess selections
            a = selector.process(a)
            b = selector.process(b)
            #
            if object is None: object=''
            if int(matchmaker)==0:
                sele1 = "((%s) in (%s))" % (str(a),str(b))
                sele2 = "((%s) in (%s))" % (str(b),str(a))
            else:
                sele1 = str(a)
                sele2 = str(b)
            try:
                _self.lock(_self)
                r = _cmd.fit(_self._COb,sele1,sele2,0,
                        int(mobile_state)-1,int(target_state)-1,
                        int(quiet),int(matchmaker),float(cutoff),
                        int(cycles),str(object))
            finally:
                _self.unlock(r,_self)
            if r < -0.5:
                raise pymol.CmdException
            return r

        def pair_fit(*arg, quiet=0, _self=cmd):
            '''
DESCRIPTION

    "pair_fit" fits matched sets of atom pairs between two objects.

USAGE

    pair_fit selection, selection, [ selection, selection [ ... ]]

EXAMPLES

    # superimpose protA residues 10-25 and 33-46 to protB residues 22-37 and 41-54:

    pair_fit protA/10-25+33-46/CA, protB/22-37+41-54/CA

    # superimpose ligA atoms C1, C2, and C4 to ligB atoms C8, C4, and C10, respectively:

    pair_fit ligA////C1, ligB////C8, ligA////C2, ligB////C4, ligA////C3, ligB////C10

NOTES

    So long as the atoms are stored in PyMOL with the same order
    internally, you can provide just two selections.  Otherwise, you
    may need to specify each pair of atoms separately, two by two, as
    additional arguments to pair_fit.

    Script files are usually recommended when using this command.

SEE ALSO

    fit, rms, rms_cur, intra_fit, intra_rms, intra_rms_cur
            '''
            r = DEFAULT_ERROR
            if len(arg) < 2:
                raise pymol.CmdException('need at least 2 selection')
            if len(arg) % 2:
                raise pymol.CmdException('need even number of selections')
            new_arg = list(map(selector.process, arg))
            try:
                _self.lock(_self)
                r = _cmd.fit_pairs(_self._COb,new_arg, quiet)
            finally:
                _self.unlock(r,_self)
            if r < -0.5:
                raise pymol.CmdException
            return r
