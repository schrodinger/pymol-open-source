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

    from . import selector

    cmd = __import__("sys").modules["pymol.cmd"]

    from .cmd import _cmd,Shortcut, \
          DEFAULT_ERROR, DEFAULT_SUCCESS, _raising, is_ok, is_error

    import pymol

    def deselect(_self=cmd):
        '''
DESCRIPTION

    "deselect" disables any and all visible selections

USAGE

    deselect

PYMOL API

    cmd.deselect()
        '''
        r = DEFAULT_SUCCESS
        arg = _self.get_names("selections",enabled_only=1)
        for a in arg:
            _self.disable(a)
        if _self._raising(r,_self): raise pymol.CmdException
        return r


    def select(name, selection="", enable=-1, quiet=1, merge=0, state=0, domain='',_self=cmd):
        '''
DESCRIPTION

    "select" creates a named atom selection from a
    selection-expression.

USAGE

    select name, selection [, enable [, quiet [, merge [, state [, domain ]]]]]

ARGUMENTS

    name = a unique name for the selection

    selection = a selection-expression

NOTES

    If a selection-expression with explicit surrounding parethenses is
    provided as the first argument, then the default selection name
    is used as the name argument.

EXAMPLES 

    select chA, chain A
    select ( resn HIS )
    select near142, resi 142 around 5

PYMOL API

    cmd.select(string name, string selection)

SEE ALSO

    delete
        '''
        with _self.lockcm:
            return _cmd.select(
                _self._COb,  #
                "" if name is None else str(name),
                str(selector.process(selection)),
                int(quiet),
                int(state) - 1,
                str(domain),
                int(enable),
                int(merge))


    def pop(name, source, enable=-1, quiet=1, _self=cmd):
        '''
DESCRIPTION

    "pop" provides a mechanism of iterating through an atom selection
    atom by atom, where each atom is sequentially assigned to the
    named selection.
    
USAGE

    pop name, source
    
EXAMPLE

    select src, name CA

    python
    while cmd.pop("tmp","src"):
        cmd.zoom("tmp",2, animate=1)
        for a in range(30):
           cmd.refresh()
           time.sleep(0.05)
    python end
    
PYMOL API

    cmd.deselect()
        '''
        r = DEFAULT_ERROR
        try:
            _self.lock(_self)
            r = _cmd.pop(_self._COb,str(name),str(source),int(quiet))
            if is_ok(r):
                enable = int(enable)
                if enable>0:
                    r = _cmd.onoff(_self._COb,str(name),1,0);
                elif enable == 0:
                    r = _cmd.onoff(_self._COb,str(name),0,0)
        finally:
            _self.unlock(r,_self)
        if _self._raising(r,_self): raise pymol.CmdException
        return r

    id_type_dict = {
        'index' : 0,
        'id'    : 1,
        'rank'  : 2,
        }

    id_type_sc = Shortcut(id_type_dict.keys())

    def select_list(name,object,id_list,state=0,mode='id',quiet=1,_self=cmd):
        '''
DESCRIPTION

    API only. Select by atom indices within a single object.

    Returns the number of selected atoms.

ARGUMENTS

    name = str: a unique name for the selection

    object = str: object name

    id_list = list of integers: ID, index, or rank list.

    state = int: object state, to limit selection to atoms which have
    coordinates in that state (-1 = current, 0 = ignore) {default: 0}

    mode = id|index|rank: {default: id}
        '''
        #
        mode = id_type_dict[id_type_sc.auto_err(mode,'identifier type')]
        with _self.lockcm:
            return _cmd.select_list(_self._COb, name, object, id_list,
                                    int(state) - 1, int(mode), int(quiet))

    def indicate(selection="(all)",_self=cmd):
        '''
DESCRIPTION

    "indicate" shows a visual representation of an atom selection.

USAGE

    indicate (selection)

PYMOL API

    cmd.count(string selection)

        '''
        r = DEFAULT_ERROR
        # preprocess selection
        selection = selector.process(selection)
        #
        try:
            _self.lock(_self)
            r = _cmd.select(_self._COb,"indicate","("+str(selection)+")",1,-1,'')
            if is_error(r):
                _self.delete("indicate")
            else:
                _self.enable("indicate")
        finally:
            _self.unlock(r,_self)
        if _self._raising(r,_self): raise pymol.CmdException
        return r

    def objsele_state_iter(selection, state=0, _self=cmd):
        '''
DESCRIPTION

    API only. Get (object-specific-selection, object-state) tuples for all
    objects in selection.
        '''
        for oname in _self.get_object_list('(' + selection + ')'):
            osele = '(%s) & ?%s' % (selection, oname)
            if state < 0:
                first = last = _self.get_object_state(oname)
            else:
                first = last = state
            if first == 0:
                first = 1
                last = _self.count_states('%' + oname)
            for ostate in range(first, last + 1):
                yield osele, ostate
