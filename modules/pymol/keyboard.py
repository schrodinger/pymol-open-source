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
    from pymol.shortcut_manager import shortcut_dict_ref, ShortcutIndex
    keys = {}
    for key in shortcut_dict_ref:
        keys[key] = shortcut_dict_ref[key][ShortcutIndex.COMMAND]

    return keys
