from . import editor
cmd = __import__("sys").modules["pymol.cmd"]
from .cmd import DEFAULT_SUCCESS, DEFAULT_ERROR

# persistent storage between copy/paste
class _PersistentEditing:

    def __init__(self,self_cmd=cmd):
        self._cmd = self_cmd
        self._obj = None

    def get_clipboard_object(self):
        '''name of the temporary persistent object'''
        return self._obj

    def create_tmp(self,sel,extract):
        '''copy or cut selection to clipboard'''
        if self._obj:
            self._cmd.delete(self._obj)
        else:
            self._obj = self._cmd.get_unused_name("_persistent_obj")

        auto_hide_sele = self._cmd.get_setting_boolean("auto_hide_selections")
        self._cmd.set("auto_hide_selections", 0, updates=0)

        try:
            self._cmd.create(self._obj, sel, extract=extract, zoom=0)
            self._cmd.disable(self._obj)
        finally:
            self._cmd.set("auto_hide_selections", auto_hide_sele, updates=0)

# user commands via keyboard

_kCopy  = 0
_kPaste = 1
_kCut   = 2

def editing_ring(action, *, _self=cmd):
    """
DESCRIPTION

    Helper function for copy/cut/paste of molecular selections.

ARGUMENTS

    action = cut/copy/paste/invert
    """

    space = getattr(_self, "_editing_ring_space", None)
    if space is None:
        space = _PersistentEditing(_self)
        _self._editing_ring_space = space

    # PASTE current hidden object into a new object
    if action in (_kPaste, "paste"):
        clipobj = space.get_clipboard_object()
        if not clipobj:
            print("Nothing on clipboard")
        else:
            _self.copy(_self.get_unused_name("obj"), clipobj, zoom=0)
        return

    sels = _self.get_names("public_selections", enabled_only=1)
    if not sels:
        print("No active selection")
        return

    sel = sels[0]
    assert sel

    # COPY current selection into a new hidden object
    if action in (_kCopy, "copy"):
        space.create_tmp(sel, extract=0)

    # CUT current selection into a new hidden object
    elif action in (_kCut, "cut"):
        space.create_tmp(sel, extract=1)

    # INVERT the current selection
    elif action == "invert":
        _self.select(sel, 'not ' + sel)

    else:
        raise AttributeError(action)


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
        'CTRL-C'        : 'editing_ring copy',
        'CTRL-F'        : 'wizard find',
        'CTRL-H'        : 'help edit_keys',
        'CTRL-I'        : 'editing_ring invert',
        'CTRL-L'        : 'util.ligand_zoom()',
        'CTRL-T'        : 'bond;unpick',
        'CTRL-V'        : 'editing_ring paste',
        'CTRL-X'        : 'editing_ring cut',
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
