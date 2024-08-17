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

import pymol
from pymol.shortcut import Shortcut
from .constants import CURRENT_STATE, ALL_STATES


class _AtomProxy:
    """
    Proxy for the "iterate" atom namespace
    """

    def __init__(self, ns):
        self.__dict__['_ns'] = ns

    def __getattr__(self, key):
        return self._ns[key]

    def __setattr__(self, key, value):
        self._ns[key] = value

    def __repr__(self):
        if self.state != 0:
            tail = f' ({self.x:.2f}, {self.y:.2f}, {self.z:.2f}) state={self.state}'
        else:
            tail = ''
        return (
            f'<{self.__class__.__name__} '
            f'/{self.model}/{self.segi}/{self.chain}/{self.resn}`{self.resi}/{self.name}`{self.alt}{tail}>'
        )

    def __dir__(self):
        from .completing import expr_sc
        return [k.split(".")[0] for k in expr_sc.keywords]

def _iterate_prepare_args(expression, space, _self):
    if not expression:
        raise pymol.CmdException('missing expression')

    if not isinstance(expression, str):
        assert callable(expression)
        assert space is None
        space = {
            "_callback": (lambda ns, f=expression, p=_AtomProxy: f(p(ns))),
            "locals": locals
        }
        expression = "_callback(locals())"
    elif space is None:
        space = _self._pymol.__dict__

    return expression, space


if True:

    import math
    from . import selector
    cmd = __import__("sys").modules["pymol.cmd"]
    from .cmd import _cmd,lock,unlock,is_string, \
          boolean_sc,boolean_dict,safe_list_eval, is_sequence, \
          DEFAULT_ERROR, DEFAULT_SUCCESS, _raising, is_ok, is_error
    from chempy import cpv

    ref_action_dict = {
        'store'     : 1,
        'recall'    : 2,
        'validate'  : 3,
        'swap'      : 4,
    }

    ref_action_sc = Shortcut(ref_action_dict.keys())

    def reference(action='validate', selection='(all)',
                  state=0, quiet=1, _self=cmd):
        r = DEFAULT_ERROR
        if is_string(action):
            action = ref_action_sc.auto_err(action,"action")
            action = ref_action_dict[action]
        else:
            action = int(action)
        selection = selector.process(selection)
        try:
            _self.lock(_self)
            r = _cmd.reference( _self._COb, int(action), str(selection),
                               int(state)-1, int(quiet))
        finally:
            _self.unlock(r,_self)
        if _self._raising(r,_self): raise pymol.CmdException
        return r

    def sculpt_purge(_self=cmd):
        '''
DESCRIPTION

    "sculpt_purge" is an unsupported feature.
    
    '''
        r = DEFAULT_ERROR
        try:
            _self.lock(_self)
            r = _cmd.sculpt_purge(_self._COb)
        finally:
            _self.unlock(r,_self)
        if _self._raising(r,_self): raise pymol.CmdException
        return r

    def sculpt_deactivate(object, _self=cmd):
        '''
DESCRIPTION

    "sculpt_deactivate" deactivates sculpting for the given object and
    clears the stored restraints.

ARGUMENTS

    object = str: name of a single object or "all"

SEE ALSO

    sculpt_activate
    '''
        r = 0
        try:
            _self.lock(_self)
            r = _cmd.sculpt_deactivate(_self._COb,str(object))
        finally:
            _self.unlock(r,_self)
        if _self._raising(r,_self): raise pymol.CmdException
        return r

    def sculpt_activate(object, state=0, match_state=-1,
                        match_by_segment=0, _self=cmd):
        '''
DESCRIPTION

    "sculpt_activate" enables sculpting for the given object. The current
    geometry (bond lengths, angles, etc.) of the given state is remembered as
    the reference geometry.

ARGUMENTS

    object = str: name of a single object or "all"

    state = int: object state or 0 for current state {default: 0}

SEE ALSO

    sculpt_iterate, sculpt_deactivate
    '''
        r = DEFAULT_ERROR
        try:
            _self.lock(_self)
            r = _cmd.sculpt_activate(_self._COb,str(object),
                                     int(state)-1,
                                     int(match_state)-1,
                                     int(match_by_segment))
        finally:
            _self.unlock(r,_self)
        if _self._raising(r,_self): raise pymol.CmdException
        return r

    def split_states(object, first=1, last=0, prefix=None, _self=cmd):
        '''
DESCRIPTION

    "split_states" separates a multi-state molecular object into a set
    of single-state molecular objects.

USAGE

    split_states object [, first [, last [, prefix ]]]
    
EXAMPLE

    load docking_hits.sdf
    split_states docking_hits, prefix=hit
    delete docking_hits
    
SEE ALSO

    join_states
        '''
        r = DEFAULT_SUCCESS
        prefix_set = bool(prefix)
        first=int(first)
        last=int(last)
        prefix_a = first - 1

        # support selections, not only objects
        sele = _self.get_unused_name('_sele_Abnqh5s5VS')
        _self.select(sele, object, 0)

        # all names to check for name conflicts
        names = set(_self.get_names('all'))

        # iterate over objects
        for object in _self.get_object_list(sele):
            olast = _self.count_states('%' + object)

            if 0 < last < olast:
                olast = last

            for a in range(first, olast + 1):
                name = None
                if prefix_set:
                    prefix_a += 1
                else:
                    prefix_a = a
                    prefix = object + "_"
                    name = _self.get_title(object,a)

                if not name:
                    name = prefix + "%04d" % prefix_a
                elif name in names:
                    name = _self.get_unused_name(name)
                names.add(name)

                r = _self.create(name, "?%s & ?%s" % (sele, object), a, 1)

                if is_error(r):
                    break

        _self.delete(sele)
        if _self._raising(r,_self): raise pymol.CmdException
        return r

    def sculpt_iterate(object, state=CURRENT_STATE, cycles=10, _self=cmd):
        '''
DESCRIPTION

    "sculpt_iterate" performs a simple energy minimization of atomic
    coordinates based on the geometry restraints which were defined with
    the "sculpt_activate" invocation and which are selected in the
    "sculpt_field_mask" setting. Sculpting currently supports local
    geometry restraints and vdw repulsion, but no solvation or
    electrostatic effects.

ARGUMENTS

    object = str: name of a single object or "all"

    state = int: object state or -1 for current state, 0 for all states
    {default: -1} (changed in PyMOL 2.5: 0 used to be "current state" as well)

    cycles = int: number of iterations {default: 10}

SEE ALSO

    commands: sculpt_activate, sculpt_deactivate
    settings: "sculpting" setting, all "sculpt_*" settings
    '''
        r = DEFAULT_ERROR
        try:
            _self.lock(_self)
            r = _cmd.sculpt_iterate(_self._COb,str(object),int(state)-1,int(cycles))
        finally:
            _self.unlock(r,_self)
        if _self._raising(r,_self): raise pymol.CmdException
        return r

    def smooth(selection="all", passes=1, window=5, first=1,
               last=0, ends=0, quiet=1, *, cutoff=-1, pbc=1, _self=cmd):

        '''
DESCRIPTION

    "smooth" performs a window average of coordinate states.  

USAGE

    smooth [ selection [, passes [, window [, first [, last [, ends]]]]]]

ARGUMENTS

    ends = 0 or 1: controls whether or not the end states are also smoothed
    using a weighted asymmetric window

    cutoff = float: Distance cutoff for atom movement between two frames.

    pbc = 0/1: Consider periodic boundary conditions {default: 1}

NOTES

    This type of averaging is often used to suppress high-frequency
    vibrations in a molecular dynamics trajectory.

SEE ALSO

    load_traj

    '''
        selection = selector.process(selection)
        with _self.lockcm:
            return _cmd.smooth(_self._COb, selection, int(passes), int(window),
                               int(first) - 1,
                               int(last) - 1, int(ends), int(quiet),
                               float(cutoff), int(pbc))

    def pbc_unwrap(oname, bymol=True, *, _self=cmd):
        '''
DESCRIPTION

    Unwrap molecules or atoms from PBC box so that they don't jump
    across periodic boundaries.

ARGUMENTS

    oname = str: object name

    bymol = 0/1: Unwrap by molecule, not by atom {default: 1}
        '''
        with _self.lockcm:
            return _cmd.pbc_unwrap(_self._COb, oname, int(bymol))

    def pbc_wrap(oname, center=None, *, _self=cmd):
        '''
DESCRIPTION

    Wrap molecules into PBC box.

ARGUMENTS

    oname = str: object name

    center = list or None: Center position in model space, or None
    to use average of first coordinate state.

EXAMPLE

    pbc_wrap trajectory, center=[0, 0, 0]
        '''
        if isinstance(center, str):
            center = _self.safe_list_eval(center)
        with _self.lockcm:
            return _cmd.pbc_wrap(_self._COb, oname, center)

    def set_state_order(name, order, quiet=1, _self=cmd):
        '''
DESCRIPTION

    API only. Set the order of states for an object.

ARGUMENTS

    name = str: object name

    order = list of int: index array (1-based state indices)

EXAMPLE

    # reverse the order of a 20 model object
    cmd.set_state_order('1nmr', range(20, 0, -1))
        '''
        with _self.lockcm:
            return _cmd.set_state_order(_self._COb, name, [i - 1 for i in order])

    def set_discrete(name, discrete=1, quiet=1, _self=cmd):
        '''
DESCRIPTION

    Convert discrete to non-discrete object or vice versa.
        '''
        with _self.lockcm:
            return _cmd.set_discrete(_self._COb, name, int(discrete))

    def set_symmetry(selection,
            a, b, c, alpha, beta, gamma, spacegroup="P1",
            state=-1, quiet=1,
            _self=cmd):

        '''
DESCRIPTION

    "set_symmetry" defines or redefines the crystal and spacegroup
    parameters for a molecule or map object.

USAGE

    set_symmetry selection, a, b, c, alpha, beta, gamma, spacegroup

ARGUMENTS

    selection = str: object name pattern

PYMOL API

    cmd.set_symmetry(string selection, float a, float b, float c,
          float alpha, float beta, float gamma, string spacegroup)

        '''
        with _self.lockcm:
            r = _cmd.set_symmetry(_self._COb,str(selection), int(state) - 1,
                                         float(a),float(b),float(c),
                                         float(alpha),float(beta),float(gamma),
                                         str(spacegroup),
                                         int(quiet))
        return r

    def symmetry_copy(source_name,
                      target_name,
                      source_state=1,
                      target_state=1,
                      quiet=1,
                      _self=cmd):
        """
DESCRIPTION

    "symmetry_copy" copies symmetry information from one object to another.

USAGE

    symmetry_copy source_name, target_name, source_state, target_state

ARGUMENTS

    source_name = str: object name
    target_name = str: object name pattern
    source_state = int: object state (maps only)
    target_state = int: object state (maps only)

NOTES

    Molecular objects don't support individual states yet.
        """
        with _self.lockcm:
            return _cmd.symmetry_copy(_self._COb,
                                   str(source_name),
                                   str(target_name),
                                   int(source_state)-1,
                                   int(target_state)-1,
                                   int(quiet))


    def set_name(old_name, new_name,_self=cmd):
        '''
DESCRIPTION

    "set_name" changes the name of an object or selection.
    
USAGE

    set_name old_name, new_name
    
PYMOL API

    cmd.set_name(string old_name, string new_name)

        '''
        r = DEFAULT_ERROR
        try:
            _self.lock(_self)
            r = _cmd.set_name(_self._COb,str(old_name),
                                    str(new_name))
        finally:
            _self.unlock(r,_self)
        if _self._raising(r,_self): raise pymol.CmdException
        return r


    def set_geometry(selection, geometry, valence, _self=cmd):
        '''
DESCRIPTION

    "set_geometry" changes PyMOL\'s assumptions about the proper valence
    and geometry of atoms in the selection.

USAGE

    set_geometry selection, geometry, valence

NOTES

    Immature functionality. See code for details.

PYMOL API

    cmd.set_geometry(string selection, int geometry, int valence)

SEE ALSO

    remove, attach, fuse, bond, unbond
    '''
        r = DEFAULT_ERROR
        # preprocess selection
        selection = selector.process(selection)
        try:
            _self.lock(_self)
            r = _cmd.set_geometry(_self._COb,str(selection),int(geometry),int(valence))
        finally:
            _self.unlock(r,_self)
        if _self._raising(r,_self): raise pymol.CmdException
        return r

    def undo(_self=cmd):
        '''
DESCRIPTION

    "undo" restores the previous conformation of the object currently
    being edited.

USAGE

    undo

SEE ALSO

    redo, push_undo
    '''
        r = DEFAULT_ERROR
        try:
            _self.lock(_self)
            r = _cmd.undo(_self._COb,-1)
        finally:
            _self.unlock(r,_self)
        if _self._raising(r,_self): raise pymol.CmdException
        return r

    def push_undo(selection, just_coordinates=1, finish_undo=0, add_objects=0, delete_objects=0, state=0, _self=cmd):
        '''
DESCRIPTION

    "push_undo" stores the current conformations of objects in the
    selection onto their individual undo rings.

    Notice: This command is only partly implemented in open-source PyMOL.

USAGE

    push_undo (all)

SEE ALSO

    undo, redo
    '''
        # preprocess selections
        selection = selector.process(selection)
        #
        r = DEFAULT_ERROR
        try:
            _self.lock(_self)
            r = _cmd.push_undo(_self._COb,"("+str(selection)+")",int(state)-1)
        finally:
            _self.unlock(r,_self)
        if _self._raising(r,_self): raise pymol.CmdException
        return r

    def redo(_self=cmd):
        '''
DESCRIPTION

    "redo" reapplies the conformational change of the object currently
    being edited.

USAGE

    redo

SEE ALSO

    undo, push_undo
    '''
        r = DEFAULT_ERROR
        try:
            _self.lock(_self)
            r = _cmd.undo(_self._COb,1)
        finally:
            _self.unlock(r,_self)
        if _self._raising(r,_self): raise pymol.CmdException
        return r

    order_dict = {
    # simulation
        '0'         : 0,
        '1'         : 1,
        '2'         : 2,
        '3'         : 3,
        '4'         : 4,
        'aromatic'  : 4,
        'guess'     : -1,
        'copy'      : -2
    }

    order_sc = Shortcut(order_dict.keys())

    def valence(order, selection1=None, selection2=None, source='',
                target_state=0, source_state=0, reset=1,
                quiet=1, *, symop="", _self=cmd):
        '''
DESCRIPTION

    "valence" modifies the valences of all existing bonds formed
    between two atom selections.
    
USAGE

    valence 2, (name C), (name O)

PYMOL API

    cmd.valence(string selection1, selection2)

SEE ALSO

    unbond, fuse, attach, replace, remove_picked
    '''
        if is_string(order):
            order = order_sc.auto_err(order,"order")
            order = order_dict[order]
        else:
            order = int(order)
        r = DEFAULT_ERROR
        # preprocess selections
        if selection1 is None:
            selection1="(pk1)"
            if selection2 is None:
                selection2="(pk2)"
        if selection2 is None:
            selection2 = selection1
        selection1 = selector.process(selection1)
        selection2 = selector.process(selection2)
        if source!='':
            source = selector.process(source)
        try:
            _self.lock(_self)
            if order>=0:
                r = _cmd.bond(_self._COb, selection1, selection2, int(order), 2, int(quiet), symop)
            else:
                r = _cmd.revalence(_self._COb,
                                   "("+selection1+")",
                                   "("+selection2+")",
                                   str(source),
                                   int(target_state)-1, int(source_state)-1,
                                   int(reset), int(quiet))
        finally:
            _self.unlock(r,_self)
        if _self._raising(r,_self): raise pymol.CmdException
        return r

    def add_bond(oname, index1, index2, order=1, _self=cmd):
        '''
DESCRIPTION

    API-only function to add a bond by atom indices (1-based, same as "index"
    in cmd.iterate()).

    To add bonds by atom selection, use cmd.bond()

ARGUMENTS

    oname = str: object name

    index1 = int: first atom index

    index2 = int: second atom index

    order = int: bond order {default: 1}

SEE ALSO

    cmd.get_bonds()
        '''
        with _self.lockcm:
            return _cmd.add_bond(_self._COb, oname, index1 - 1, index2 - 1, order)

    def rebond(oname, state=CURRENT_STATE, *, pbc=1, _self=cmd):
        '''
DESCRIPTION

    Discard all bonds and do distance based bonding.

ARGUMENTS

    oname = str: object name

    state = int: object state {default: -1 (current state)}

    pbc = 0/1: Use periodic boundary conditions (only if symmetry
    is defined for the object) {default: 1}
        '''
        with _self.lockcm:
            return _cmd.rebond(_self._COb, oname, int(state) - 1, int(pbc))

    def bond(atom1="pk1", atom2="pk2", order=1, *, quiet=1, symop="", _self=cmd):
        '''
DESCRIPTION

    "bond" creates a new bond between two selections, each of which
    should contain one atom.

USAGE

    bond [atom1, atom2 [,order]]

ARGUMENTS

    atom1 = str: Atom selection of first atom {default: pk1}

    atom2 = str: Atom selection of second atom {default: pk2}

    order = int: Bond order {default: 1}

    symop = str: Symmetry operation code for second atom (e.g. "1_555")

NOTES

    The atoms must both be within the same object.

SEE ALSO

    unbond, fuse, attach, replace, remove_picked
    '''
        r = DEFAULT_ERROR
        # preprocess selections
        atom1 = selector.process(atom1)
        atom2 = selector.process(atom2)
        try:
            _self.lock(_self)
            r = _cmd.bond(_self._COb,
                          atom1, atom2,
                          int(order),1,int(quiet), symop)
        finally:
            _self.unlock(r,_self)
        if _self._raising(r,_self): raise pymol.CmdException
        return r

    def invert(quiet=1, _self=cmd):
        '''
DESCRIPTION

    "invert" inverts the stereo-chemistry of atom (pk1), holding attached atoms
    (pk2) and (pk3) immobile.

USAGE

    invert 

NOTES

    The invert function is usually bound to CTRL-E in Editing Mode.

PYMOL API

    cmd.invert( )

    '''
        #
        with _self.lockcm:
            return _cmd.invert(_self._COb, int(quiet))

    def unbond(atom1="(pk1)", atom2="(pk2)", quiet=1, _self=cmd):
        '''
DESCRIPTION

    "unbond" removes all bonds between two selections.

USAGE

    unbond atom1,atom2

ARGUMENTS

    atom1 = string {default: (pk1)}

    atom2 = string {default: (pk2)}

PYMOL API

    cmd.unbond(selection atom1, selection atom2)

SEE ALSO

    bond, fuse, remove_picked, attach, detach, replace

    '''
        r = DEFAULT_ERROR
        # preprocess selections
        atom1 = selector.process(atom1)
        atom2 = selector.process(atom2)
        try:
            _self.lock(_self)
            r = _cmd.bond(_self._COb,
                          atom1, atom2,
                          0,0,int(quiet))
        finally:
            _self.unlock(r,_self)
        if _self._raising(r,_self): raise pymol.CmdException
        return r

    def remove(selection, quiet=1, _self=cmd):
        '''
DESCRIPTION

    "remove" eleminates the atoms in a selection from their respective
    molecular objects.

USAGE

    remove selection

EXAMPLES

    remove resi 124 

PYMOL API

    cmd.remove( string selection )

SEE ALSO

    delete
    '''
        r = DEFAULT_ERROR
        # preprocess selection
        selection = selector.process(selection)
        #
        r = 1
        try:
            _self.lock(_self)
            r = _cmd.remove(_self._COb,"("+selection+")",int(quiet))
        finally:
            _self.unlock(r,_self)
        if _self._raising(r,_self): raise pymol.CmdException
        return r


    def remove_picked(hydrogens=1,quiet=1,_self=cmd):
        '''
DESCRIPTION

    "remove_picked" removes the atom or bond currently picked for
    editing.

USAGE

    remove_picked [ hydrogens ]

NOTES

    This function is usually connected to the
    DELETE key and "CTRL-D".

    By default, attached hydrogens will also be deleted unless
    hydrogen-flag is zero.

PYMOL API

    cmd.remove_picked(integer hydrogens)

SEE ALSO

    attach, replace
    '''
        r = DEFAULT_ERROR
        try:
            _self.lock(_self)
            r = _cmd.remove_picked(_self._COb,int(hydrogens),int(quiet))
        finally:
            _self.unlock(r,_self)
        if _self._raising(r,_self): raise pymol.CmdException
        return r


    def cycle_valence(h_fill=1, quiet=1, _self=cmd):
        '''
DESCRIPTION

    "cycle_valence" cycles the valence on the currently selected bond.

USAGE

    cycle_valence [ h_fill ]

ARGUMENTS

    h_fill = 0 or 1: updated hydrogens too? {default: 1 (yes)}
    
EXAMPLE

    cycle_valence

NOTES

    If the h_fill flag is true, hydrogens will be added or removed to
    satisfy valence requirements.

    This function is usually connected to the DELETE key and "CTRL-W".

PYMOL API

    cmd.cycle_valence(int h_fill)

SEE ALSO

    remove_picked, attach, replace, fuse, h_fill
    '''
        r = DEFAULT_ERROR
        try:
            _self.lock(_self)
            r = _cmd.cycle_valence(_self._COb,quiet)
        finally:
            _self.unlock(r,_self)
        if h_fill:
            globals()['h_fill'](quiet)
        if _self._raising(r,_self): raise pymol.CmdException
        return r


    def attach(element,geometry,valence,name='',quiet=1,_self=cmd):
        '''
DESCRIPTION

    "attach" adds a single atom on to the picked atom.

USAGE

    attach element, geometry, valence

PYMOL API

    cmd.attach( element, geometry, valence )

    '''
        with _self.lockcm:
            return _cmd.attach(_self._COb,str(element),int(geometry),int(valence),str(name))

    def fuse(selection1="(pk1)", selection2="(pk2)",
             mode=0, recolor=1, move=1, _self=cmd):
        '''
DESCRIPTION

    "fuse" joins two objects into one by forming a bond.  A copy of
    the object containing the first atom is moved so as to form an
    approximately resonable bond with the second, and that copy is
    then merged with the first object.

USAGE

    fuse [ selection1 [, selection2 [, mode [, recolor [, move ]]]]]

ARGUMENTS

    selection1 = str: single atom selection (will be copied to object 2)

    selection2 = str: single atom selection

    mode = int: {default: 0}
      3: don't move and don't create a bond, just combine into single object

    recolor = bool: recolor C atoms to match target {default: 1}

    move = bool: {default: 1}

NOTES

    Each selection must include a single atom in each object.
    The atoms can both be hydrogens, in which case they are
    eliminated, or they can both be non-hydrogens, in which
    case a bond is formed between the two atoms.

SEE ALSO

    bond, unbond, attach, replace, fuse, remove_picked
    '''
        r = DEFAULT_ERROR
        # preprocess selections
        selection1 = selector.process(selection1)
        selection2 = selector.process(selection2)
        #
        try:
            _self.lock(_self)
            r = _cmd.fuse(_self._COb,str(selection1),str(selection2),
                              int(mode),int(recolor),int(move))
        finally:
            _self.unlock(r,_self)
        if _self._raising(r,_self): raise pymol.CmdException
        return r

    def unpick(_self=cmd):
        '''
DESCRIPTION

    "unpick" deletes the special "pk" atom selections (pk1, pk2, etc.)
    used in atom picking and molecular editing.

USAGE

    unpick

PYMOL API

    cmd.unpick()

SEE ALSO

    edit
        '''

        r = DEFAULT_ERROR
        try:
            _self.lock(_self)
            r = _cmd.unpick(_self._COb)
        finally:
            _self.unlock(r,_self)
        if _self._raising(r,_self): raise pymol.CmdException
        return r

    def drag(selection=None, wizard=1, edit=1, quiet=1, mode=-1, _self=cmd):
        '''
DESCRIPTION

    "drag" activates dragging for a selection, enabling the user to
    manipulate the atom coordinates of the atoms using mouse controls
    similar to those for controlling the camera.

USAGE

    drag [ selection ]

ARGUMENTS
  
    selection = string: atoms to drag.  If not provided, and dragging
    is active, then dragging is instead deactivated.

NOTES

    Currently, the selection of atom to drag must all reside in a
    single molecular object.

'''
        import pymol.wizard.dragging

        quiet = int(quiet)
        if (selection is not None) and (selection!=""):
            selection = selector.process(selection)
            if is_string(edit):
                edit=boolean_dict[boolean_sc.auto_err(edit,'boolean')]
            if is_string(wizard):
                wizard=boolean_dict[boolean_sc.auto_err(wizard,'boolean')]
            edit = int(edit)
            wizard = int(wizard)
            old_button_mode = _self.get('button_mode')
        else:
            wizard = 0
            edit = 0
            selection = ""
        #
        r = DEFAULT_ERROR
        try:
            _self.lock(_self)
            r = _cmd.drag(_self._COb,str(selection),int(quiet),int(mode))
        finally:
            _self.unlock(r,_self)
        if not is_error(r):
            if edit:
                _self.edit_mode(edit)
            if wizard:
                wiz = _self.get_wizard()
                if (wiz is None):
                    _self.wizard("dragging",old_button_mode)
                elif not isinstance(wiz, pymol.wizard.dragging.Dragging):
                    _self.wizard("dragging",old_button_mode)
                else:
                    wiz.recount()
        if _self._raising(r,_self): raise pymol.CmdException
        return r

    def edit(selection1='', selection2='none', selection3='none',
             selection4='none', pkresi=0, pkbond=1, quiet=1, _self=cmd):
        '''
DESCRIPTION

    "edit" picks atoms or a bond for editing.

USAGE

    edit selection1 [, selection2 [, selection3 [, selection4 [, pkresi [, pkbond ]]]]] 

NOTES

    If only one selection is provided, an atom is picked.

    If two selections are provided, the bond between them
    is picked (by default, if one exists).

PYMOL API

    cmd.edit(string selection1, string selection2,
             string selection3, string selection4,
             int pkresi, int pkbond, int quiet)

SEE ALSO

    unpick, remove_picked, cycle_valence, torsion
    '''
        # preprocess selections
        selection1 = selector.process(selection1)
        selection2 = selector.process(selection2)
        selection3 = selector.process(selection3)
        selection4 = selector.process(selection4)
        #
        r = DEFAULT_ERROR
        try:
            _self.lock(_self)
            r = _cmd.edit(_self._COb,str(selection1),str(selection2),
                              str(selection3),str(selection4),
                              int(pkresi),int(pkbond),int(quiet))
        finally:
            _self.unlock(r,_self)
        if _self._raising(r,_self): raise pymol.CmdException
        return r

    def get_editor_scheme(_self=cmd):
        r = DEFAULT_ERROR
        try:
            _self.lock(_self)
            r = _cmd.get_editor_scheme(_self._COb)
        finally:
            _self.unlock(r,_self)
        if _self._raising(r,_self): raise pymol.CmdException
        return r

    def torsion(angle,_self=cmd):
        '''
DESCRIPTION

    "torsion" rotates the torsion on the bond currently
    picked for editing.  The rotated fragment will correspond
    to the first atom specified when picking the bond (or the
    nearest atom, if picked using the mouse).

USAGE

    torsion angle

PYMOL API

    cmd.torsion( float angle )

SEE ALSO

    edit, unpick, remove_picked, cycle_valence
    '''
        r = DEFAULT_ERROR
        try:
            _self.lock(_self)
            r = _cmd.torsion(_self._COb,float(angle))
        finally:
            _self.unlock(r,_self)
        if _self._raising(r,_self): raise pymol.CmdException
        return r

    def h_fill(quiet=1, _self=cmd):
        '''
DESCRIPTION

    "h_fill" removes and replaces hydrogens on the atom or bond picked
    for editing.

USAGE

    h_fill

NOTES

    This is useful for fixing hydrogens after changing bond valences.

PYMOL API

    cmd.h_fill()

SEE ALSO

    edit, cycle_valence, h_add
    '''
        r = DEFAULT_ERROR
        try:
            _self.lock(_self)
            r = _cmd.h_fill(_self._COb,int(quiet))
        finally:
            _self.unlock(r,_self)
        if _self._raising(r,_self): raise pymol.CmdException
        return r

    def h_fix(selection="",quiet=1,_self=cmd):
        '''
DESCRIPTION

    "h_fix" is an unsupported command that may have something to do
    with repositioning hydrogen atoms.
    
    '''
        # preprocess selection
        selection = selector.process(selection)

        r = DEFAULT_ERROR
        try:
            _self.lock(_self)
            r = _cmd.h_fix(_self._COb,str(selection),int(quiet))
        finally:
            _self.unlock(r,_self)
        if _self._raising(r,_self): raise pymol.CmdException
        return r


    def h_add(selection="(all)", quiet=1, state=0, legacy=0, _self=cmd):
        '''
DESCRIPTION

    "h_add" adds hydrogens onto a molecule based on current valences.

USAGE

    h_add [ selection [, state ]]

ARGUMENTS

    selection = string {default: (all)}

    state = int {default: 0 (all states)}

NOTES

    Because PDB files do not normally contain bond valences for
    ligands and other nonstandard components, it may be necessary to
    manually correct ligand conformations before adding hydrogens.

SEE ALSO

    h_fill
    '''
        # preprocess selection
        selection = selector.process(selection)
        #
        r = DEFAULT_ERROR
        try:
            _self.lock(_self)
            r = _cmd.h_add(_self._COb,selection,int(quiet),
                    int(state) - 1, int(legacy))
        finally:
            _self.unlock(r,_self)
        if _self._raising(r,_self): raise pymol.CmdException
        return r



    def sort(object="",_self=cmd):
        '''
DESCRIPTION

    "sort" reorders atoms in the structure.  It usually only necessary
    to run this routine after an "alter" command which has modified the
    names of atom properties.  Without an argument, sort will resort
    all atoms in all objects.

USAGE

    sort [object]

PYMOL API

    cmd.sort(string object)

SEE ALSO

    alter
    '''
        r = DEFAULT_ERROR
        try:
            _self.lock(_self)
            r = _cmd.sort(_self._COb,str(object))
        finally:
            _self.unlock(r,_self)
        if _self._raising(r,_self): raise pymol.CmdException
        return r


    def replace(element, geometry, valence, h_fill=1, name="",
		quiet=1, _self=cmd):
        '''
DESCRIPTION

    "replace" replaces the picked atom with a new atom.

USAGE

    replace element, geometry, valence [, h_fill [, name ]]

NOTES

    Immature functionality. See code for details.

PYMOL API

    cmd.replace(string element, int geometry, int valence, int h_fill,
                string name)

SEE ALSO

    remove, attach, fuse, bond, unbond
    '''
        r = DEFAULT_ERROR
        if "pk1" not in _self.get_names("selections"):
            print(" Error: you must first pick an atom to replace.")
            raise pymol.CmdException
        try:
            if h_fill: # strip off existing hydrogens
                remove("((neighbor pk1) and elem H)",quiet=quiet)
            _self.lock(_self)
            r = _cmd.replace(_self._COb,str(element),int(geometry),int(valence),str(name),quiet)
        finally:
            _self.unlock(r,_self)
        if _self._raising(r,_self): raise pymol.CmdException
        return r

    def rename(selection="all", force=0, quiet=1, _self=cmd):
        '''
DESCRIPTION

    "rename" creates new atom names which are unique within residues.

USAGE

    rename selection [, force ]

PYMOL API

    cmd.rename(string selection, int force )

SEE ALSO

    alter
    '''
        r = DEFAULT_ERROR
        selection = "("+selector.process(selection)+")"
        try:
            _self.lock(_self)
            r = _cmd.rename(_self._COb,str(selection),int(force),int(quiet))
        finally:
            _self.unlock(r,_self)
        if _self._raising(r,_self): raise pymol.CmdException
        return r

    def dss(selection="(all)", state=0, context=None, preserve=0,
            quiet=1, _self=cmd):

        '''
DESCRIPTION

    "dss" defines secondary structure based on backbone geometry
    and hydrogen bonding patterns.
    
USAGE

    dss selection, state

ARGUMENT

    selection = string: {default: (all)}

    state = integer: {default: 0 -- all states}
    
EXAMPLE

    dss

NOTES

    With PyMOL, heavy emphasis is placed on cartoon aesthetics, and so
    both hydrogen bonding patterns and backbone geometry are used in
    the assignment process.  Depending upon the local context, helix
    and strand assignments are made based on geometry, hydrogen
    bonding, or both.

    This command will generate results which differ slightly from DSSP
    and other programs.  Most deviations occur in borderline or
    transition regions.  Generally speaking, PyMOL is more strict, thus
    assigning fewer helix/sheet residues, except for partially
    distorted helices, which PyMOL tends to tolerate.
    
    WARNING: This algorithm has not yet been rigorously validated.

    If you dislike one or more of the assignments made by dss, you can
    use the alter command to make changes (followed by "rebuild").
    For example:
    
        alter 123-125/, ss=\'L\'
        alter pk1, ss=\'S\'
        alter 90/, ss=\'H\'
        rebuild

PYMOL API

    cmd.dss(string selection, int state)
        
        '''
        # preprocess selections
        selection = selector.process(selection)
        r = DEFAULT_ERROR
        if context is None:
            context = ""
        else:
            context = selector.process(context)
        #
        try:
            _self.lock(_self)
            r = _cmd.dss(_self._COb,str(selection),int(state)-1,str(context),
                            int(preserve),int(quiet))
        finally:
            _self.unlock(r,_self)
        if _self._raising(r,_self): raise pymol.CmdException
        return r

    def alter(selection, expression, quiet=1, space=None, _self=cmd):

        '''
DESCRIPTION

    "alter" changes atomic properties using an expression evaluated
    within a temporary namespace for each atom.

USAGE

    alter selection, expression

EXAMPLES

    alter chain A, chain='B'
    alter all, resi=str(int(resi)+100)
    sort

NOTES

    Symbols defined (* = read only):

    name, resn, resi, resv, chain, segi, elem, alt, q, b, vdw, type,
    partial_charge, formal_charge, elec_radius, text_type, label, 
    numeric_type, model*, state*, index*, ID, rank, color, ss,
    cartoon, flags

    All strings must be explicitly quoted.  This operation typically
    takes several seconds per thousand atoms altered.  

    You may need to issue a "rebuild" in order to update associated
    representations.
    
    WARNING: You should always issue a "sort" command on an object
    after modifying any property which might affect canonical atom
    ordering (names, chains, etc.).  Failure to do so will confound
    subsequent "create" and "byres" operations.  

SEE ALSO

    alter_state, iterate, iterate_state, sort
        '''
        expression, space = _iterate_prepare_args(expression, space, _self)

        # preprocess selections
        selection = selector.process(selection)

        with _self.lockcm:
            return _cmd.alter(_self._COb, selection, expression, False,
                              int(quiet), dict(space))

    def alter_list(object, expr_list, quiet=1, space=None, _self=cmd):
        '''
DESCRIPTION

    "alter_list" is an unsupported feature.
    
        '''
        if space is None:
            space = _self._pymol.__dict__

        with _self.lockcm:
            return _cmd.alter_list(_self._COb, object, list(expr_list),
                                   int(quiet), dict(space))


    def iterate(selection, expression, quiet=1, space=None, _self=cmd):

        '''
DESCRIPTION

    "iterate" iterates over an expression within a temporary namespace
    for each atom.

USAGE

    iterate selection, expression

EXAMPLES

    stored.net_charge = 0
    iterate all, stored.net_charge = stored.net_charge + partial_charge
    print(stored.net_charge)
    
    stored.names = []
    iterate all, stored.names.append(name)
    print(stored.names)
    
    # Using a Python callback (new in PyMOL 2.5)
    names = []
    cmd.iterate("all", lambda atom: names.append(atom.name))
    print(names)

NOTES

    Unlike with the "alter" command, atomic properties cannot be
    altered.  Other than that, the commands are identical.

SEE ALSO

    iterate_state, alter, alter_state
        '''
        expression, space = _iterate_prepare_args(expression, space, _self)

        # preprocess selection
        selection = selector.process(selection)

        with _self.lockcm:
            return _cmd.alter(_self._COb, selection, expression, True,
                              int(quiet), dict(space))

    def alter_state(state, selection, expression, quiet=1,
                    space=None, atomic=1, _self=cmd):

        '''
DESCRIPTION

    "alter_state" changes atom coordinates and flags over a particular
    state and selection using the Python evaluator with a temporary
    namespace for each atomic coordinate.

USAGE

    alter_state state, selection, expression

EXAMPLES

    alter_state 1, all, x=x+5
    rebuild
    
NOTES

    By default, most of the symbols from "alter" are available for use
    on a read-only basis.  

    It is usually necessary to "rebuild" representations once your
    alterations are complete.
    
SEE ALSO

    iterate_state, alter, iterate
        '''
        expression, space = _iterate_prepare_args(expression, space, _self)

        # preprocess selection
        selection = selector.process(selection)
        #
        state = int(state)

        with _self.lockcm:
            return _cmd.alter_state(_self._COb,
                                    int(state) - 1, selection, expression,
                                    False, int(quiet), dict(space))

    def iterate_state(state, selection, expression, quiet=1,
                      space=None, atomic=1, _self=cmd):

        '''
DESCRIPTION

    "iterate_state" is to "alter_state" as "iterate" is to "alter"

USAGE

    iterate_state state, selection, expression

EXAMPLES

    stored.sum_x = 0.0
    iterate_state 1, all, stored.sum_x = stored.sum_x + x
    print(stored.sum_x)
    
SEE ALSO

    iterate, alter, alter_state
        '''
        expression, space = _iterate_prepare_args(expression, space, _self)

        # preprocess selection
        selection = selector.process(selection)

        with _self.lockcm:
            return _cmd.alter_state(_self._COb,
                                    int(state) - 1, selection, expression,
                                    True, int(quiet), dict(space))

    def translate(vector=[0.0,0.0,0.0], selection="all", state=-1,
                  camera=1, object=None, object_mode=0, _self=cmd):

        '''
DESCRIPTION

    "translate" translates the atomic coordinates of atoms in a
    selection.  Alternatively, is modifies the matrix associated with
    a particular object or object-state.

USAGE

    translate vector [, selection [, state [, camera [, object ]]]]

ARGUMENTS

    vector = float vector: translation vector

    selection = string: atoms whose coordinates should be modified {default: all}

    state > 0: only the indicated state is modified

    state = 0: all states are modified

    state = -1: only current state is modified {default}

    camera = 0 or 1: is the vector in camera coordinates? {default: 1 (yes)}

    object = string: object name (only if rotating object matrix) {default: None}


PYMOL API

    cmd.translate(list vector, string selection, int state, int
                  camera, string object)

EXAMPLES

    translate [1,0,0], name CA

NOTES

    "translate" can be used to translate the atomic coordinates of a
    molecular object.  Behavior differs depending on whether or not
    the "object" parameter is specified.

    If object is None, then translate translates atomic coordinates
    according to the vector provided for the selection and in the state
    provided.  All representation geometries will need to be
    regenerated to reflect the new atomic coordinates.

    If object is set to an object name, then selection is ignored and
    instead of translating the atomic coordinates, the object\'s
    overall representation display matrix is modified.  This option is
    for use in animations only.

    The "camera" option controls whether the camera or the model\'s
    axes are used to interpret the translation vector.

        '''
        r = DEFAULT_ERROR
        object_mode = int(object_mode)
        if _self.is_string(vector):
            vector = safe_list_eval(vector)
        if not _self.is_sequence(vector):
            print("Error: bad vector.")
            raise pymol.CmdException
        else:
            vector = [float(vector[0]),float(vector[1]),float(vector[2])]
            selection = selector.process(selection)
            camera=int(camera)
            view = _self.get_view(0)
            if camera:
                mat = [ view[0:3],view[3:6],view[6:9] ]
                shift = cpv.transform(mat,vector)
            else:
                shift = vector
            if object is None:
                ttt = [1.0,0.0,0.0,shift[0],
                         0.0,1.0,0.0,shift[1],
                         0.0,0.0,1.0,shift[2],
                         0.0,0.0,0.0,1.0]
                r=_self.transform_selection(selection,ttt,state=state)
            elif object_mode==0: # update the TTT display matrix
                try:
                    _self.lock(_self)
                    r=_cmd.translate_object_ttt(_self._COb,str(object),shift)
                finally:
                    _self.unlock(r,_self)
            elif object_mode==1: # either updates TTT or coordinates & history
                     # depending on the current matrix mode
                matrix = [1.0, 0.0, 0.0, shift[0],
                          0.0, 1.0, 0.0, shift[1],
                          0.0, 0.0, 1.0, shift[2],
                          0.0, 0.0, 0.0, 1.0]
                try:
                    _self.lock(_self)
                    r = _cmd.transform_object(_self._COb,str(object),int(state)-1,
                                                      list(matrix),0,'',1)
                finally:
                    _self.unlock(r,_self)
            else:
                print(" Error: translate: unrecognized object_mode")
        if _self._raising(r,_self): raise pymol.CmdException
        return r

    def rotate(axis='x', angle=0.0, selection="all", state=-1, camera=1,
               object=None, origin=None, object_mode=0, _self=cmd):

        '''
DESCRIPTION

    "rotate" rotates the atomic coordinates of atoms in a selection
    about an axis.  Alternatively, it modifies the matrix associated
    with a particular object or object state.

USAGE

    rotate axis, angle [, selection [, state [, camera [, object 
        [, origin ]]]]]

ARGUMENTS

    axis = x, y, z, or float vector: axis about which to rotate

    angle = float: degrees of rotation

    selection = string: atoms whose coordinates should be modified {default: all}
    
    state > 0: only the indicated state is modified

    state = 0: all states are modified

    state = -1: only the current state is modified {default}

    camera = 0 or 1: is the axis specific in camera coordinates? {default: 1}

    object = string: object name (only if rotating object matrix) {default: None}

    origin = float vector: origin of rotateion {default: None}
        
EXAMPLES

    rotate x, 45, pept

    rotate [1,1,1], 10, chain A

NOTES

    Behavior differs depending on whether or not the "object"
    parameter is specified.

    If object is None, then the atomic coordinates are modified
    directly, and all representation geometries will need to be
    regenerated to reflect the new atomic coordinates.

    If object is set to an object name, then the selection field is
    ignored and instead of translating the atomic coordinates, the
    object matrix is modified.  This option is only intended for use
    in animations and is not yet fully supported.

PYMOL API

    cmd.rotate(list-or-string axis, float angle, string selection, 
               int state, int camera, string object)

        '''
        r = DEFAULT_ERROR
        object_mode = int(object_mode)
        have_origin = 0
        if axis in ['x','X']:
            axis = [1.0,0.0,0.0]
        elif axis in ['y','Y']:
            axis = [0.0,1.0,0.0]
        elif axis in ['z','Z']:
            axis = [0.0,0.0,1.0]
        else:
            axis = safe_list_eval(str(axis))
        if not _self.is_list(axis):
            print("Error: bad axis.")
            raise pymol.CmdException
        else:
            axis = [float(axis[0]),float(axis[1]),float(axis[2])]
            angle = math.pi*float(angle)/180.0
            view = _self.get_view(0)
            if origin is not None:
                have_origin = 1
                if _self.is_string(origin):
                    if ',' in origin:
                        origin = safe_list_eval(origin) # should be a sequence of floats
                    else:
                        _self.lock(_self)
                        try:
                            origin = _cmd.get_origin(_self._COb,str(origin))
                        finally:
                            unlock(-1)
                origin = [float(origin[0]),float(origin[1]),float(origin[2])]
            else:
                origin = [view[12],view[13],view[14]]
            camera=int(camera)
            if camera:
                vmat = [ view[0:3],view[3:6],view[6:9] ]
                axis = cpv.transform(vmat,axis)
            mat = cpv.rotation_matrix(angle,axis)
            if object is None:
                ttt = [mat[0][0],mat[0][1],mat[0][2],origin[0],
                       mat[1][0],mat[1][1],mat[1][2],origin[1],
                       mat[2][0],mat[2][1],mat[2][2],origin[2],
                       -origin[0],-origin[1],-origin[2], 1.0]
                r=_self.transform_selection(selection,ttt,state=state)
            elif object_mode==0:
                _self.lock(_self)
                try:
                    if not have_origin:
                        origin = _cmd.get_origin(_self._COb,str(object))
                    if is_sequence(origin):
                        ttt = [mat[0][0],mat[0][1],mat[0][2], origin[0],
                               mat[1][0],mat[1][1],mat[1][2], origin[1],
                               mat[2][0],mat[2][1],mat[2][2], origin[2],
                               -origin[0], -origin[1], -origin[2], 1.0]
                        r=_cmd.combine_object_ttt(_self._COb,str(object),ttt)
                finally:
                    _self.unlock(r,_self)
                if not is_sequence(origin):
                    print(" Error: rotate: unknown object '%s'."%object)
                    if _self._raising(r,_self):
                        raise pymol.CmdException
            elif object_mode==1:

                matrix = [mat[0][0],mat[0][1],mat[0][2], origin[0],
                          mat[1][0],mat[1][1],mat[1][2], origin[1],
                          mat[2][0],mat[2][1],mat[2][2], origin[2],
                          -origin[0],-origin[1],-origin[2], 1.0]
                try:
                    _self.lock(_self)
                    r = _cmd.transform_object(_self._COb,str(object),int(state)-1,
                                                      list(matrix),0,'',0)
                finally:
                    _self.unlock(r,_self)

            else:
                print(" Error: rotate: unrecognized object_mode")
        if _self._raising(r,_self): raise pymol.CmdException
        return r

    def look_at(target_obj: str, mobile_obj: str = "_Camera", *, _self=cmd):
        """
DESCRIPTION

    "look_at" modifies a rotation of an object (or view) so that its forward (z axis)
    faces the center of a target object.

ARGUMENTS

    mobile_obj the object to rotate

    target_obj the object to look at

USAGE

    look_at target_obj

PYMOL API

    cmd.look_at(string target_obj)

    """
        with _self.lockcm:
            return _cmd.look_at(_self._COb, target_obj, mobile_obj)


    def move_on_curve(mobile_obj: str,
                      curve_obj: str,
                      t: float,
                      *,
                      _self=cmd) -> None:
        '''
DESCRIPTION

    "move_on_curve" moves an object along a curve.

ARGUMENTS

    mobile_obj the object to move
    curve_obj the curve to move the object along

USAGE

    move_on_curve mobile_obj, curve_obj, t

PYMOL API

    cmd.move_on_curve(string mobile_obj, string curve_obj, float t)
        '''
        with _self.lockcm:
            return _cmd.move_on_curve(_self._COb, mobile_obj, curve_obj, t)


    def set_title(object, state, text, _self=cmd):
        '''
DESCRIPTION

    "set_title" attaches a text string to the state of a particular
    object which can be displayed next to the object name when that
    state is active.  This is useful for display the energies of a set
    of conformers.

USAGE

    set_title object, state, text

PYMOL API

    cmd.set_title(string object, int state, string text)

    '''
        r = DEFAULT_ERROR
        try:
            _self.lock(_self)
            r = _cmd.set_title(_self._COb,str(object),int(state)-1,str(text))
        finally:
            _self.unlock(r,_self)

    def set_object_ttt(object, ttt, state=0, quiet=1, homogenous=0,
                       _self=cmd):
        '''
DESCRIPTION

    "set_object_ttt" is an API-only function which sets the TTT matrix
    (view transformation) for an object.

    When a movie is defined and the object has key frames for object
    motions, then the key frames take priority and update the TTT matrix
    while the movie is playing.

    Unlike a homogenous matrix where the last row is always [0,0,0,1],
    a TTT matrix may have a pre-translation vector in the last row.

ARGUMENTS

    object = str: object name

    ttt = list of 16 floats: TTT matrix

    state = int: UNUSED, TTT matrices are not state specific

    homogenous = 0/1: NAME IS MISLEADING AND IMPLEMENTATION
    POSSIBLY WRONG! If 1, then transpose the input matrix and
    set the last column (post-translation) to [0,0,0,1].

SEE ALSO

    cmd.transform_object, cmd.matrix_reset
        '''
        r = None
        if _self.is_string(ttt):
            ttt = safe_list_eval(str(ttt))
        if homogenous: # passed a homogenous matrix, so do the best we can
            ttt = [ # NOTE: this appears to be incorrect...
                ttt[ 0], ttt[ 4], ttt[ 8], 0.0,
                ttt[ 1], ttt[ 5], ttt[ 9], 0.0,
                ttt[ 2], ttt[ 6], ttt[10], 0.0,
                ttt[ 3], ttt[ 7], ttt[11], 1.0]
        try:
            _self.lock(_self)
            r = _cmd.set_object_ttt(_self._COb,str(object),
				    (float(ttt[ 0]),
				     float(ttt[ 1]),
				     float(ttt[ 2]),
				     float(ttt[ 3]),
				     float(ttt[ 4]),
				     float(ttt[ 5]),
				     float(ttt[ 6]),
				     float(ttt[ 7]),
				     float(ttt[ 8]),
				     float(ttt[ 9]),
				     float(ttt[10]),
				     float(ttt[11]),
				     float(ttt[12]),
				     float(ttt[13]),
				     float(ttt[14]),
				     float(ttt[15])),
				    int(state)-1,int(quiet))
        finally:
            _self.unlock(r,_self)
        if _self._raising(r,_self): raise pymol.CmdException
        return r

    def transform_selection(selection, matrix, state=-1, log=0,
                            homogenous=0, transpose=0, _self=cmd):
        '''

DESCRIPTION

    "transform_selection" transforms the atomic coordinates of a
    selection.

PYMOL API
   
    cmd.transform_selection(string selection, list matrix, int state,
                            int log, int homogenous, int transpose):

NOTES

    Note that when homogenous is zero, the input matrix is NOT a
    standard homogenous 4x4 transformation matrix.  Instead it is
    something PyMOL-specific which consists of the following:

    1) a 3x3 matrix containing the rotation in the upper-left quadrant

    2) a 1x3 translation to be applied *before* rotation in the bottom row
        (matrix[12],matrix[13],matrix[14]).

    3) a 3x1 translation to be applied *after* rotation in the right-hand
        column (matrix[3],matrix[7],matrix[11])

    In other words, if the matrix is:

    [  m0  m1  m2  m3 \\
       m4  m5  m6  m7 \\
       m8  m9 m10 m11 \\
      m12 m13 m14 m15 ] 

    Atoms will be transformed as follows

    Y = M X

    y0 = m0*(x0+m12) + m1*(x1+m13) +  m2*(x2+m14) + m3 \\
    y1 = m4*(x0+m12) + m5*(x1+m13) +  m6*(x2+m14) + m7 \\
    y2 = m8*(x0+m12) + m9*(x1+m13) + m10*(x2+m14) + m11 

        '''
        r = DEFAULT_ERROR
        selection = selector.process(selection)
        if int(transpose):
            matrix = [ matrix[0], matrix[4], matrix[8 ], matrix[12],
                       matrix[1], matrix[5], matrix[9 ], matrix[13],
                       matrix[2], matrix[6], matrix[10], matrix[14],
                       matrix[3], matrix[7], matrix[11], matrix[15]]
        try:
            _self.lock(_self)
            r = _cmd.transform_selection(_self._COb,str(selection),int(state)-1,
                                                  list(matrix),int(log),int(homogenous))
        finally:
            _self.unlock(r,_self)
        if _self._raising(r,_self): raise pymol.CmdException
        return r

    def transform_object(name, matrix, state=-1, log=0, selection='',
                         homogenous=0, transpose=0, _self=cmd):
        '''
DESCRIPTION

    "transform_object" in an API-only function which applies a
    transformation matrix to an object.

    If setting "matrix_mode" > 0 and selection is empty, then this
    function operates on the TTT (movie) matrix.

ARGUMENTS

    name = str: object name

    matrix = list of 16 floats: transformation matrix

    state = int: object state {default: -1}

    log = 0/1: write action to log file (only applies if object is a
    molecular object) {default: 0}

    selection = str: atom selection (only applies to molecular objects,
    if empty then the whole object state is transformed).

    homogenous = 0/1: if 0, then matrix[12:15] may contain a pre-translation,
    otherwise those values must be zeros (see also cmd.transform_selection)
    {default: 0}

    transpose = 0/1: matrix is 0=row-major, 1=column-major {default: 0}

SEE ALSO

    cmd.transform_selection, cmd.set_object_ttt, cmd.matrix_reset
        '''
        r = DEFAULT_ERROR
        if int(transpose):
            matrix = [ matrix[0], matrix[4], matrix[8 ], matrix[12],
                          matrix[1], matrix[5], matrix[9 ], matrix[13],
                          matrix[2], matrix[6], matrix[10], matrix[14],
                          matrix[3], matrix[7], matrix[11], matrix[15]]
        try:
            _self.lock(_self)
            r = _cmd.transform_object(_self._COb,str(name),int(state)-1,
                                              list(matrix),int(log),
                                              str(selection),
                                              int(homogenous))
        finally:
            _self.unlock(r,_self)
        if _self._raising(r,_self): raise pymol.CmdException
        return r

    def matrix_copy(source_name='', target_name='',
                    source_mode=-1, target_mode=-1,
                    source_state=1, target_state=1,
                    target_undo=1, log=0, quiet=1,
                    _self=cmd):
        '''

DESCRIPTION
        
    "matrix_copy" copies a transformation matrix from one object to
    another.

USAGE

    matrix_copy source_name, target_name

NOTES

    This command is often used after a protein structure alignment to
    bring other related objects into the same frame of reference.

SEE ALSO

    matrix_reset, align, fit, pair_fit
'''

        r = DEFAULT_ERROR
        if source_name is None:
            source_name = ''
        target_name = str(target_name).strip()
        source_name = str(source_name).strip()
        if (target_name == '' and source_name != ''): # tentative -- create a new command instead?
            mat = _self.get_object_matrix(source_name,source_state)
            view = _self.get_view()
            new_view = (
                mat[ 0]*view[ 0] + mat[ 4]*view[ 3] + mat[ 8]*view[ 6],
                mat[ 0]*view[ 1] + mat[ 4]*view[ 4] + mat[ 8]*view[ 7],
                mat[ 0]*view[ 2] + mat[ 4]*view[ 5] + mat[ 8]*view[ 8],
                mat[ 1]*view[ 0] + mat[ 5]*view[ 3] + mat[ 9]*view[ 6],
                mat[ 1]*view[ 1] + mat[ 5]*view[ 4] + mat[ 9]*view[ 7],
                mat[ 1]*view[ 2] + mat[ 5]*view[ 5] + mat[ 9]*view[ 8],
                mat[ 2]*view[ 0] + mat[ 6]*view[ 3] + mat[10]*view[ 6],
                mat[ 2]*view[ 1] + mat[ 6]*view[ 4] + mat[10]*view[ 7],
                mat[ 2]*view[ 2] + mat[ 6]*view[ 5] + mat[10]*view[ 8],
                view[ 9] ,         view[10] ,         view[11],
                mat[ 0]*view[12] + mat[ 1]*view[13] + mat[ 2]*view[14] -
                mat[ 0]* mat[ 3] - mat[ 4]* mat[ 7] - mat[ 8]* mat[11],
                mat[ 4]*view[12] + mat[ 5]*view[13] + mat[ 6]*view[14] -
                mat[ 1]* mat[ 3] - mat[ 5]* mat[ 7] - mat[ 9]* mat[11],
                mat[ 8]*view[12] + mat[ 9]*view[13] + mat[10]*view[14] -
                mat[ 2]* mat[ 3] - mat[ 6]* mat[ 7] - mat[10]* mat[11],
                view[15] ,         view[16] ,         view[17] )
            r = _self.set_view(new_view)
        else:
            with _self.lockcm:
                r = _cmd.matrix_copy(_self._COb,str(source_name),
                                                 str(target_name),
                                                 int(source_mode),
                                                 int(target_mode),
                                                 int(source_state)-1,
                                                 int(target_state)-1,
                                                 int(target_undo),
                                                 int(log),
                                                 int(quiet))
            if _self._raising(r,_self):
                raise pymol.CmdException
        return r

    def matrix_reset(name, state=1, mode=-1, log=0, quiet=1,_self=cmd):
        '''

DESCRIPTION
        
    "matrix_reset" resets the transformation for an object.

USAGE

    matrix_reset name [, state [, mode ]]

ARGUMENTS

    name = str: object name

    state = int: object state {default: 1}

    mode = int: {defualt: -1 = matrix_mode or 0}
      0: transformation was applied to coordinates
      1: reset TTT matrix (movie transformation)
      2: reset state matrix

SEE ALSO

    matrix_copy, align, super, fit, pair_fit
    
'''

        r = DEFAULT_ERROR
        try:
            _self.lock(_self)
            r = _cmd.reset_matrix(_self._COb,str(name),
                                         int(mode),
                                         int(state)-1,
                                         int(log),
                                         int(quiet))
        finally:
            _self.unlock(r,_self)
        if _self._raising(r,_self):
            raise pymol.CmdException
        return r


    def translate_atom(sele1, v0, v1, v2, state=0, mode=0,
                       log=0, _self=cmd):

        r = DEFAULT_ERROR
        sele1 = selector.process(sele1)
        try:
            _self.lock(_self)
            r = _cmd.translate_atom(_self._COb,str(sele1),float(v0),float(v1),
                                            float(v2),int(state)-1,int(mode),int(log))
        finally:
            _self.unlock(r,_self)
        if _self._raising(r,_self): raise pymol.CmdException
        return r

    def update(target, source, target_state=0, source_state=0,
               matchmaker=1, quiet=1, _self=cmd):

        '''
DESCRIPTION

    "update" transfers coordinates from one selection to another.
USAGE

    update (target-selection),(source-selection)

EXAMPLES

    update target,(variant)

NOTES

    Currently, this applies across all pairs of states.  Fine
    control will be added later.

SEE ALSO

    load
    '''
        r = DEFAULT_ERROR
        a=target
        b=source
        # preprocess selections
        a = selector.process(a)
        b = selector.process(b)
        #
        if a[0]!='(': a="("+str(a)+")"
        if b[0]!='(': b="("+str(b)+")"
        try:
            _self.lock(_self)
            r = _cmd.update(_self._COb,str(a),str(b),int(target_state)-1,
                                 int(source_state)-1,int(matchmaker),int(quiet))
        finally:
            _self.unlock(r,_self)
        if _self._raising(r,_self): raise pymol.CmdException
        return r


    def set_dihedral(atom1, atom2, atom3, atom4, angle, state=1, quiet=1, _self=cmd):
        '''
DESCRIPTION

    "set_dihedral" changes the dihedral angle formed between the four
    bonded atoms provided.  The atoms must be acyclic.
    
USAGE

    set_dihedral atom1, atom2, atom3, atom4, angle [, state [, quiet ]]

NOTES

    Because set_dihedral uses the molecular editing capability,
    numbered "pk" atom selections (if any) will be redefined by this
    operation.
    
PYMOL API

    cmd.set_dihedral(string atom1, string atom2, string atom3, string atom4,
                     float angle, int state, int quiet)

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
            r = _cmd.set_dihe(_self._COb,atom1,atom2,atom3,atom4,
                                    float(angle),int(state)-1,int(quiet))
        finally:
            _self.unlock(r,_self)
        if _self._raising(r,_self): raise pymol.CmdException
        return r

    map_op_dict = {
        'minimum'       : 0,
        'maximum'       : 1,
        'sum'           : 2,
        'average'       : 3,
        'difference'    : 4,
        'copy'          : 5,
        'unique'        : 6,
        }

    map_op_sc = Shortcut(map_op_dict.keys())

    def map_set(name, operator, operands='', target_state=0,
                source_state=0, zoom=0, quiet=1, _self=cmd):

        '''
DESCRIPTION

    "map_set" provides a number of common operations on and between maps.

USAGE

    map_set name, operator, operands, target_state, source_state

    operator may be "minimum, maximum, average, sum, or difference"

EXAMPLES

    map my_sum, add, map1 map2 map3
    map my_avg, average, map1 map2 map3
    
NOTES

    source_state = 0 means all states
    target_state = -1 means current state
    
    experimental
    
SEE ALSO

    map_new
        '''
        r = DEFAULT_ERROR
        operator_index = map_op_dict[map_op_sc.auto_err(operator,'operator')]
        try:
            _self.lock(_self)
            r = _cmd.map_set(_self._COb,str(name), int(operator_index), str(operands),
                         int(target_state)-1, int(source_state)-1, int(zoom), int(quiet))
        finally:
            _self.unlock(r,_self)
        if _self._raising(r,_self): raise pymol.CmdException
        return r

    def map_set_border(name, level=0.0, state=0, _self=cmd):

        '''
DESCRIPTION

    "map_set_border" is a function (reqd by PDA) which allows you to set the
    level on the edge points of a map

USAGE

    map_set_border name, level

NOTES

    unsupported.
    
SEE ALSO

    load
        '''
        r = DEFAULT_ERROR
        try:
            _self.lock(_self)
            r = _cmd.map_set_border(_self._COb,str(name),float(level),int(state))
        finally:
            _self.unlock(r,_self)
        if _self._raising(r,_self): raise pymol.CmdException
        return r

    def map_double(name, state=0, _self=cmd):
        '''
DESCRIPTION

    "map_double" resamples a map at twice the current resolution.

NOTES

     The amount of memory required to store the map will increase
     eight-fold.

USAGE

    map_double map_name, state
        '''
        r = DEFAULT_ERROR
        try:
            _self.lock(_self)
            r = _cmd.map_double(_self._COb,str(name),int(state)-1)
        finally:
            _self.unlock(r,_self)
        if _self._raising(r,_self): raise pymol.CmdException
        return r

    def map_halve(name, state=0, smooth=1, _self=cmd):

        '''
DESCRIPTION

    "map_halve" resamples a map at half the current resolution.  

USAGE

    map_halve map_name, state
    
NOTES

    The amount of memory required to store the map will decrease
    eight-fold.

SEE ALSO

    map_double
        '''

        r = DEFAULT_ERROR
        try:
            _self.lock(_self)
            r = _cmd.map_halve(_self._COb,str(name),int(state)-1,int(smooth))
        finally:
            _self.unlock(r,_self)
        if _self._raising(r,_self): raise pymol.CmdException
        return r

    def map_trim(name, selection, buffer=0.0, map_state=0, sele_state=0, quiet=1, _self=cmd):
        '''
DESCRIPTION

    "map_trim" is an unsupported command that may have something to do with
    reducing the extent of a map to cover just a single selection of atoms.

    '''

        r = DEFAULT_ERROR
        # preprocess selection
        selection = selector.process(selection)
        #
        try:
            _self.lock(_self)
            r = _cmd.map_trim(_self._COb,str(name),str(selection),
                                  float(buffer),int(map_state)-1,
                                  int(sele_state)-1,int(quiet))
        finally:
            _self.unlock(r,_self)
        if _self._raising(r,_self): raise pymol.CmdException
        return r


    def protect(selection="(all)", quiet=1, _self=cmd):
        '''
DESCRIPTION

    "protect" protects a set of atoms from tranformations performed
    using the editing features.  This is most useful when you are
    modifying an internal portion of a chain or cycle and do not wish
    to affect the rest of the molecule.

USAGE

    protect (selection)

PYMOL API

    cmd.protect(string selection)

SEE ALSO

    deprotect, mask, unmask, mouse, editing
    '''
        r = DEFAULT_ERROR
        # preprocess selection
        selection = selector.process(selection)
        #
        try:
            _self.lock(_self)
            r = _cmd.protect(_self._COb,"("+str(selection)+")",1,int(quiet))
        finally:
            _self.unlock(r,_self)
        if _self._raising(r,_self): raise pymol.CmdException
        return r

    def deprotect(selection="(all)", quiet=1, _self=cmd):
        '''
DESCRIPTION

    "deprotect" reverses the effect of the "protect" command.

USAGE

    deprotect (selection)

PYMOL API

    cmd.deprotect(string selection)

SEE ALSO

    protect, mask, unmask, mouse, editing
    '''
        r = DEFAULT_ERROR
        # preprocess selection
        selection = selector.process(selection)
        #
        try:
            _self.lock(_self)
            r = _cmd.protect(_self._COb,"("+str(selection)+")",0,int(quiet))
        finally:
            _self.unlock(r,_self)
        if _self._raising(r,_self): raise pymol.CmdException
        return r

    flag_dict = {
    # simulation
        'focus'         : 0,
        'free'          : 1,
        'restrain'      : 2,
        'fix'           : 3,
        'exclude'       : 4,
        'study'         : 5,

    # rendering
        'exfoliate'     : 24,
        'ignore'        : 25,
        'no_smooth'     : 26,
    }

    flag_sc = Shortcut(flag_dict.keys())

    flag_action_dict = {
        'reset' : 0,
        'set'   : 1,
        'clear' : 2,
        }

    flag_action_sc = Shortcut(flag_action_dict.keys())

    def fix_chemistry(selection1="all", selection2="all",
		      invalidate=1, quiet=1, _self=cmd):
        '''
DESCRIPTION
   
    "fix chemistry" is an unsupported feature.

'''
        selection1 = selector.process(selection1)
        selection2 = selector.process(selection2)
        with _self.lockcm:
            return _cmd.fix_chemistry(_self._COb,str(selection1),
				   str(selection2),int(invalidate),
				   int(quiet))

    def set_object_color(name, color, quiet=1, _self=cmd):
        r = DEFAULT_ERROR
        try:
            _self.lock(_self)
            r = _cmd.set_object_color(_self._COb,str(name),str(color),int(quiet))
        finally:
            _self.unlock(r,_self)
        if _raising(r,_self): raise pymol.CmdException
        return r

    def flag(flag, selection, action="reset", quiet=1, _self=cmd):
        '''
DESCRIPTION

    "flag" sets the indicated flag for atoms in the selection and
     clears the indicated flag for atoms not in the selection.  

USAGE

    flag flag, selection [, action ]

ARGUMENTS

    flag = str or int: Flag name or number

    selection = str: atom selection

    action = reset: {default} set flag for atoms in selection and clear it for all others

    action = set: set the flag for atoms in selection, leaving other atoms unchanged

    action = clear: clear the flag for selected atoms, leaving other atoms unchanged

EXAMPLES  

    flag free, (resi 45 x; 6)

NOTES

    This is primarily useful for passing selection information into
    Chempy models, which have a 32 bit attribute "flag" which holds
    this information.

    If the 'auto_indicate_flags' setting is true, then PyMOL will automatically
    create a selection called "indicate" which contains all atoms with that flag
    after applying the command.

    SPECIAL FLAGS

    * Flags 0-5 are reserved for molecular modeling

        focus      0 = Atoms of Interest (i.e. a ligand in an active site) \\
        free       1 = Free Atoms (free to move subject to a force-field) \\
        restrain   2 = Restrained Atoms (typically harmonically contrained) \\
        fix        3 = Fixed Atoms (no movement allowed) \\
        exclude    4 = Atoms which should not be part of any simulation
        study      5

    * Flags 6-7 are for protein and nucleic acid classification

    * Flags 8-15 are free for end users to manipulate

    * Flags 16-21 are reserved for external GUIs and linked applications

    * Flags 22-23 are for temporary use only (flag 23 used for coverage
      tracking when assigning parameters in chempy.champ.assign)

    * Flags 24-31 are reserved for PyMOL internal usage

        exfoliate 24 = DEPRECATED: Use "hide surface, sele" instead \\
        ignore    25 = Ignore atoms altogether when surfacing \\
        no_smooth 26 = Do not smooth atom position

PYMOL API

    cmd.flag(int flag, string selection, string action="reset",
             int indicate=0)

        '''
        if isinstance(flag, str) and not flag.isdigit():
            flag = flag_dict[flag_sc.auto_err(flag, "flag")]
        else:
            flag = int(flag)

        # preprocess selection
        selection = selector.process(selection)
        action = flag_action_dict[flag_action_sc.auto_err(action,'action')]

        with _self.lockcm:
            return _cmd.flag(_self._COb, flag, selection, action, int(quiet))

    def vdw_fit(selection1, selection2, state1=1,state2=1,buffer=0.24,quiet=1,_self=cmd):
        '''
DESCRIPTION

    "vdw_fit" is an unsupported feature.
    
    '''
        # preprocess selections
        selection1 = selector.process(selection1)
        selection2 = selector.process(selection2)
        #
        r = DEFAULT_ERROR
        try:
            _self.lock(_self)
            r = _cmd.vdw_fit(_self._COb,str(selection1),int(state1)-1,
                             str(selection2),int(state2)-1,
                             float(buffer),int(quiet))
        finally:
            _self.unlock(r,_self)
        if _self._raising(r,_self): raise pymol.CmdException
        return r

    def split_chains(selection='(all)', prefix=None, group=None, quiet=1, _self=cmd):
        '''
DESCRIPTION

    Create a single object for each chain in selection

SEE ALSO

    split_states
        '''
        count = 0
        models = _self.get_object_list('(' + selection + ')')
        names_list = []
        for model in models:
            for chain in _self.get_chains('(%s) and model %s' % (selection, model)):
                count += 1
                if not prefix:
                    name = '%s_%s' % (model, chain)
                else:
                    name = '%s%04d' % (prefix, count)
                _self.create(name, '(%s) and model %s and chain "%s"' %
                        (selection, model, chain), zoom=0)
                names_list.append(name)
            _self.disable(model)

        if group:
            _self.group(group, ' '.join(names_list), 'add')

    def alphatoall(selection='polymer', properties='b', operator='byca',
            quiet=1, _self=cmd):
        '''
DESCRIPTION

    Expand any given property of the CA atoms to all atoms in the residue

ARGUMENTS

    selection = string: atom selection {default: polymer}

    properties = string: space separated list of atom properties {default: b}
        '''
        properties = '(' + ','.join(properties.split()) + ')'
        space = {'props': {}}
        _self.iterate('%s (%s)' % (operator, selection),
                'props[model,segi,chain,resi] = ' + properties,
                space=space)
        _self.alter(selection,
                properties + ' = props.get((model,segi,chain,resi), ' + properties + ')',
                space=space)
        if not int(quiet):
            print(' Modified %d residues' % (len(space['props'])))

    def mse2met(selection='all', quiet=1, _self=cmd):
        '''
DESCRIPTION

    Mutate selenomethionine to methionine
        '''
        x = _self.alter('(%s) and MSE/SE' % selection, 'name="SD";elem="S"')
        _self.flag('ignore', '(%s) and resn MSE' % (selection), 'clear')
        _self.alter('(%s) and resn MSE' % selection, 'resn="MET";type="ATOM"')
        if not int(quiet):
            print(' Altered %d MSE residues to MET' % (x))


    def _base(i, numerals, _emptyzero=False):
        if i == 0:
            return '' if _emptyzero else numerals[0]

        b = len(numerals)
        return _base(i // b, numerals, True) + numerals[i % b]


    def uniquify(identifier, selection, reference='', quiet=1, _self=cmd):
        '''
DESCRIPTION

    Make `identifier` unique with respect to reference selection.

ARGUMENTS

    identifier = str: atom identifier (chain, segi, etc.)

    selection = str: atom selection to modify

    reference = str: atom selection whose identifiers must not be
    present in the first selection {default: !selection}

EXAMPLE

    fetch 1a00 1hbb, async=0
    uniquify chain, 1hbb
    # 1hbb now has chains E,F,G,H
        '''
        import itertools
        import string

        if not reference:
            reference = '!(' + selection + ')'

        p = identifier

        set_ref = set()
        set_sel = set()
        mapping = {}
        space = {'set_ref': set_ref, 'set_sel': set_sel, 'mapping': mapping}

        _self.iterate(reference, 'set_ref.add(' + p + ')', space=space)
        _self.iterate(selection, 'set_sel.add(' + p + ')', space=space)

        set_union = set_ref | set_sel
        set_inter = set_ref & set_sel

        if not set_inter:
            return

        baseargs = ()
        i_iter = itertools.count(1)

        if isinstance(next(iter(set_inter)), int):
            basefunc = int
        elif p in ('resi',):
            basefunc = str
        else:
            basefunc = _base
            baseargs = (string.ascii_uppercase + '123456789',)
            i_iter = itertools.count(0)

        for name in set_inter:
            for i in i_iter:
                newname = basefunc(i, *baseargs)
                if newname not in set_union:
                    mapping[name] = newname
                    set_union.add(newname)
                    break

        _self.alter(selection, '%s = mapping.get(%s, %s)' % (p, p, p),
                space=space)

        if not int(quiet):
            print(' Uniquify: renamed %d %s identifier(s)' % (len(mapping), p))


    def copy_to(name, selection, rename='chain segi ID', zoom=-1, quiet=1, _self=cmd):
        '''
DESCRIPTION

    Copies selection to object `name` (all states) and by default
    renames chain, segi and ID identifiers to avoid naming conflicts.

ARGUMENTS

    name = str: object name to modify

    selection = str: atom selection (will be copied to `name`)

    rename = str: space separated list of identifiers to rename
    {default: chain segi ID}

SEE ALSO

    create, fuse
        '''
        temp = _self.get_unused_name('_tmp')

        try:
            _self.create(temp, selection, zoom=0)
            _self.disable(' '.join(_self.get_object_list(selection)))

            for prop in rename.split():
                if prop.upper() == 'ID':
                    _self.alter(temp, 'ID = -1')
                else:
                    uniquify(prop, temp, '?' + name, quiet=quiet)

            _self.create(name, '?' + name + ' ' + temp, zoom=zoom)
            _self.unpick()

            if not int(quiet):
                n = _self.count_atoms(temp)
                print(' Copied %d atoms to object %s' % (n, name))
        finally:
            _self.delete(temp)
