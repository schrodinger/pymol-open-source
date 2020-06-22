from . import editor
cmd = __import__("sys").modules["pymol.cmd"]
from .cmd import DEFAULT_SUCCESS, DEFAULT_ERROR

# persistent storage between copy/paste
class _PersistentEditing:

    def __init__(self,self_cmd=cmd):

        # private-like

        self._cmd = self_cmd

        # init == 0 : uninitialized cmd and no data
        # init == 1 : initialized cmd and no data
        # init == 2 : initialized cmd and data
        self._init = 0

        # string name of the temporary persistent object

        self._obj = None

        # last active selection

        self._sel = None

    def deferred_init(self):
        # ecah of these methods will have to start with
        # if not self.init ...
        # because this is imported while 'cmd' is being created
        # and we need cmd functionality; this should be resolved
        # in pymol2
        if self._init == 0:
            self._obj = self._cmd.get_unused_name("_persistent_obj")
            self._init = 1

    def set_sel(self,sel):
        self._sel = sel

    def clean_obj(self):
        self._cmd.delete(self._obj)
        self._init = 1

    def create_tmp(self,sel,extract):

        if sel is None: return

        # creates an invisible temporary persistent
        # object for pasting; it is created via CTRL-C
        # for copy or CTRL-X for cut, differing in the
        # extract flag
        self.deferred_init()

        self.set_sel(sel)

        # in case something goes wrong, try & finally
        try:
            if self._init == 2:
                self.clean_obj()

            self._cmd.set("suspend_updates", 1)

            self._cmd.create(self._obj, self._sel, extract=extract, zoom=0)
            self._cmd.disable(self._obj)
            self._cmd.enable(self._sel)

            # success, now we have data

            self._init = 2

        finally:
            self._cmd.set("suspend_updates", 0)

# in pymol2 this needs to be abstracted better
# by making _persistent a member of Cmd so we can
# call the cmd functions as necessary; for now, it's
# global to this module
_persistent = _PersistentEditing(cmd)

# user commands via keyboard

_kCopy  = 0
_kPaste = 1
_kCut   = 2

def editing_ring(action, space=_persistent, self_cmd=cmd):

    sel = get_active_selection_name(self_cmd)

    space.set_sel(sel)

    # COPY current selection into a new hidden object
    if action==_kCopy:
        if sel:
            space.create_tmp(sel,extract=0)

    # CUT current selection into a new hidden object
    elif action==_kCut:
        if sel:
            space.create_tmp(sel,extract=1)

    # PASTE current hidden object into a new object
    elif action==_kPaste:

        if space._init < 2: return

        try:
            self_cmd.set("suspend_updates", 1)

            # copying or pasting; enable that object
            # and create the duplicate

            self_cmd.enable(space._obj)

            self_cmd.copy(self_cmd.get_unused_name("obj"),space._obj,zoom=0)

            # re-hide the temporary

            self_cmd.disable(space._obj)

        finally:
            self_cmd.set("suspend_updates", 0)


#
# collapse down to one function in a switch
#
def get_selection_copy_name(pfx,self_cmd=cmd):
    # obj => obj_copy01
    # obj_copy01 => obj_copy02
    # obj_copy02 => obj_copy03

    # WRONG -- get_unused_name doesn't work as expected

    import re
    if len(re.findall("_copy\d+$",pfx))!=0:
        return self_cmd.get_unused_name(pfx)
    else:
        return self_cmd.get_unused_name(pfx+"_copy")

def invert_active_selection(self_cmd=cmd):
    sel = get_active_selection_name(self_cmd)
    if sel:
        self_cmd.select(sel, 'not %s' % sel)

def get_active_selection_name(self_cmd=cmd):
    # find the name of the active selection

    for sel in self_cmd.get_names("public_selections", enabled_only=1):
        return sel

    return None

def get_default_keys(_self=cmd):
    _self = _self._weakrefproxy

    emptydict = {}
    keys = {
        'left'          : ( _self.backward              , (), emptydict ),
        'right'         : ( _self.forward               , (), emptydict ),
        'pgup'          : 'scene action=previous',
        'pgdn'          : 'scene action=next',
        'home'          : 'zoom animate=-1',
        'end'           : 'mtoggle',
        'insert'        : 'rock',
        'SHFT-left'     : 'backward',
        'SHFT-right'    : 'forward',
        'SHFT-pgup'     : 'scene action=previous',
        'SHFT-pgdn'     : 'scene action=next',
        'SHFT-home'     : 'rewind',
        'SHFT-end'      : 'ending',
        'SHFT-insert'   : 'rock',
        # F1-F2: scene Fn, store
        'CTRL-left'     : 'backward',
        'CTRL-right'    : 'forward',
        'CTRL-pgup'     : ( _self.scene                 , ('', 'insert_before'), emptydict ),
        'CTRL-pgdn'     : ( _self.scene                 , ('', 'insert_after'), emptydict ),
        'CTRL-home'     : 'zoom animate=-1',
        'CTRL-end'      : 'scene new, store',
        'CTRL-insert'   : 'scene auto, store',
        'CTRL-A'        : 'select sele, all, 1',
        'CTRL-C'        : ( editing_ring                , (),  {'action': _kCopy, 'space':_persistent,'self_cmd':_self}),
        'CTRL-F'        : 'wizard find',
        'CTRL-H'        : 'help edit_keys',
        'CTRL-I'        : ( invert_active_selection     , (),  emptydict ),
        'CTRL-L'        : 'util.ligand_zoom()',
        'CTRL-T'        : 'bond;unpick',
        'CTRL-V'        : ( editing_ring                , (),  {'action': _kPaste, 'space':_persistent,'self_cmd':_self}),
        'CTRL-X'        : ( editing_ring                , (),  {'action': _kCut, 'space':_persistent,'self_cmd':_self}),
        'CTRL-Y'        : 'redo',
        'CTRL-Z'        : 'undo',
        'ALT-left'      : 'backward',
        'ALT-right'     : 'forward',
        'ALT-pgup'      : 'rewind',
        'ALT-pgdn'      : 'ending',
        'ALT-home'      : 'zoom animate=-1',
        'ALT-end'       : 'ending',
        'ALT-insert'    : 'rock',
        'ALT-1'         : ( editor.attach_fragment      , ("pk1","formamide",5,0), emptydict),
        'ALT-2'         : ( editor.attach_fragment      , ("pk1","formamide",4,0), emptydict),
        'ALT-3'         : ( editor.attach_fragment      , ("pk1","sulfone",3,1), emptydict),
        'ALT-4'         : ( editor.attach_fragment      , ("pk1","cyclobutane",4,0), emptydict),
        'ALT-5'         : ( editor.attach_fragment      , ("pk1","cyclopentane",5,0), emptydict),
        'ALT-6'         : ( editor.attach_fragment      , ("pk1","cyclohexane",7,0), emptydict),
        'ALT-7'         : ( editor.attach_fragment      , ("pk1","cycloheptane",8,0), emptydict),
        'ALT-8'         : ( editor.attach_fragment      , ("pk1","cyclopentadiene",5,0), emptydict),
        'ALT-9'         : ( editor.attach_fragment      , ("pk1","benzene",6,0), emptydict),
        'ALT-0'         : ( editor.attach_fragment      , ("pk1","formaldehyde",2,0), emptydict),
        'ALT-A'         : ( editor.attach_amino_acid    , ("pk1","ala"), emptydict),
        'ALT-B'         : ( editor.attach_amino_acid    , ("pk1","ace"), emptydict),
        'ALT-C'         : ( editor.attach_amino_acid    , ("pk1","cys"), emptydict),
        'ALT-D'         : ( editor.attach_amino_acid    , ("pk1","asp"), emptydict),
        'ALT-E'         : ( editor.attach_amino_acid    , ("pk1","glu"), emptydict),
        'ALT-F'         : ( editor.attach_amino_acid    , ("pk1","phe"), emptydict),
        'ALT-G'         : ( editor.attach_amino_acid    , ("pk1","gly"), emptydict),
        'ALT-H'         : ( editor.attach_amino_acid    , ("pk1","his"), emptydict),
        'ALT-I'         : ( editor.attach_amino_acid    , ("pk1","ile"), emptydict),
        'ALT-J'         : ( editor.attach_fragment      , ("pk1","acetylene",2,0), emptydict),
        'ALT-K'         : ( editor.attach_amino_acid    , ("pk1","lys"), emptydict),
        'ALT-L'         : ( editor.attach_amino_acid    , ("pk1","leu"), emptydict),
        'ALT-M'         : ( editor.attach_amino_acid    , ("pk1","met"), emptydict),
        'ALT-N'         : ( editor.attach_amino_acid    , ("pk1","asn"), emptydict),
        'ALT-P'         : ( editor.attach_amino_acid    , ("pk1","pro"), emptydict),
        'ALT-Q'         : ( editor.attach_amino_acid    , ("pk1","gln"), emptydict),
        'ALT-R'         : ( editor.attach_amino_acid    , ("pk1","arg"), emptydict),
        'ALT-S'         : ( editor.attach_amino_acid    , ("pk1","ser"), emptydict),
        'ALT-T'         : ( editor.attach_amino_acid    , ("pk1","thr"), emptydict),
        'ALT-V'         : ( editor.attach_amino_acid    , ("pk1","val"), emptydict),
        'ALT-W'         : ( editor.attach_amino_acid    , ("pk1","trp"), emptydict),
        'ALT-Y'         : ( editor.attach_amino_acid    , ("pk1","tyr"), emptydict),
        'ALT-Z'         : ( editor.attach_amino_acid    , ("pk1","nme"), emptydict),
        # F1-F2: scene SHFT-Fn, store
        'CTSH-left'     : 'backward',
        'CTSH-right'    : 'forward',
        'CTSH-pgup'     : ( _self.scene                 , ('', 'insert_before'), emptydict ),
        'CTSH-pgdn'     : ( _self.scene                 , ('', 'insert_after'), emptydict ),
        'CTSH-home'     : 'zoom animate=-1',
        'CTSH-end'      : 'mtoggle',
        'CTSH-insert'   : 'rock',
        'CTSH-A'        : 'redo',
        'CTSH-B'        : 'replace Br,1,1',
        'CTSH-C'        : 'replace C,4,4',
        'CTSH-D'        : 'remove_picked',
        'CTSH-E'        : 'invert',
        'CTSH-F'        : 'replace F,1,1',
        'CTSH-G'        : 'replace H,1,1',
        'CTSH-I'        : 'replace I,1,1',
        'CTSH-J'        : 'alter pk1,formal_charge=-1.',
        'CTSH-K'        : 'alter pk1,formal_charge=1.',
        'CTSH-L'        : 'replace Cl,1,1',
        'CTSH-N'        : 'replace N,4,3',
        'CTSH-O'        : 'replace O,4,2',
        'CTSH-P'        : 'replace P,4,1',
        'CTSH-R'        : 'h_fill',
        'CTSH-S'        : 'replace S,4,2',
        'CTSH-T'        : 'bond;unpick',
        'CTSH-U'        : 'alter pk1,formal_charge=0.',
        'CTSH-W'        : 'cycle_valence',
        'CTSH-X'        : ( _self.auto_measure          , (), emptydict ),
        'CTSH-Y'        : 'attach H,1,1',
        'CTSH-Z'        : 'undo',
    }

    # F1-F12 scene ..., store
    for i in range(1, 13):
        Fn = 'F' + str(i)
        keys['CTRL-' + Fn] = ( _self.scene, (          Fn, 'store'), emptydict )
        keys['CTSH-' + Fn] = ( _self.scene, ('SHFT-' + Fn, 'store'), emptydict )

    return keys
