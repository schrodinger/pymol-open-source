from enum import IntEnum

from pymol import setting
from pymol import save_shortcut as shortcut_saver
from pymol.keyboard import get_default_keys
from pymol.shortcut_dict import shortcut_dict_ref

class ShortcutIndex(IntEnum):
    COMMAND = 0
    DESCRIPT = 1
    USER_DEF = 2

class ShortcutManager():
    def __init__(self, saved_shortcuts, cmd):
        self.cmd = cmd
        self.saved_shortcuts = saved_shortcuts
        self.default_bindings = get_default_keys(self.cmd)
        self.cmd.shortcut_dict = {key: list(value) for key,value in shortcut_dict_ref.items()}

        # Tuple of keys that are reserved for the system
        self.reserved_keys = ('CTRL-S','CTRL-E','CTRL-O','CTRL-M','up','down')

    def check_saved_dict(self):
        '''
        Checks if there is a saved file. Updates shortcut_dict if needed.
        '''
        for key, saved_shortcut_list in self.saved_shortcuts.items():
            if key in self.cmd.shortcut_dict:
                self.cmd.shortcut_dict[key][ShortcutIndex.USER_DEF] = saved_shortcut_list[ShortcutIndex.USER_DEF]
            else:
                self.cmd.shortcut_dict.update({key:saved_shortcut_list})

    def check_key_mappings(self):
        '''
        Corrects mismatched key mappings and adds missing keys.
        There are three locations of shortcuts:
        1) cmd.key_mappings: Current state of key mappings for PyMOL
        2) cmd.default_bindings: Default values of key mappings
        3) cmd.shortcut_dict: Belongs to this class instance, for table
        '''
        current_mappings = self.cmd.key_mappings
        missing_keys = []
        mismatch_keys = []
        empty_keys = []
        manual_reset_keys = []

        for key in current_mappings:
            if not current_mappings[key]:
                empty_keys.append(key)
            if key in self.default_bindings:
                if current_mappings[key] != self.default_bindings[key]:
                    mismatch_keys.append(key)
                elif current_mappings[key] == self.default_bindings[key] and self.cmd.shortcut_dict[key][2]:
                    manual_reset_keys.append(key)
            else:
                missing_keys.append(key)
        for key in mismatch_keys:
            self.cmd.shortcut_dict[key][ShortcutIndex.USER_DEF] = str(current_mappings[key])
        for key in manual_reset_keys:
            self.cmd.shortcut_dict[key][ShortcutIndex.USER_DEF] = ''
        for key in missing_keys:
            if key in self.cmd.shortcut_dict:
                self.cmd.shortcut_dict[key][ShortcutIndex.USER_DEF] = str(current_mappings[key])
            else:
                # User defined keys using set_key outside of menu
                self.cmd.shortcut_dict[key] = ['','',str(current_mappings[key])]
        for key in empty_keys:
            self.cmd.shortcut_dict[key][ShortcutIndex.USER_DEF] = 'Deleted'

    def save_shortcuts(self):
        '''
        Removes unused keys before calling save_shortcuts.
        These unused keys represent keybinding defined by the user but then assigned to nothing.
        '''
        self.remove_unused(self.cmd.shortcut_dict)
        shortcut_saver.save_shortcuts(self.cmd.shortcut_dict)

    def remove_unused(self, full_dict):
        '''
        Removes items from dictionary if first and third element are missing from the value list. 
        '''
        unused_keys = []
        for key in full_dict:
            if not full_dict[key][0] and not full_dict[key][2]:
                unused_keys.append(key)
        for key in unused_keys:
            del full_dict[key]
        return full_dict

    def reset_all_default(self):
        '''
        Iterates over all values to restore their default commands.
        This will restore them to the values from keyboard.py
        '''
        delete_keys = []

        for key in self.cmd.shortcut_dict:
            self.cmd.shortcut_dict[key][ShortcutIndex.USER_DEF] = ''

            if key in self.default_bindings:
                binding = self.default_bindings[key]
                if isinstance(binding,str):
                    self.cmd.set_key(key,binding)
                elif isinstance(binding,tuple):
                    self.cmd.set_key(key,binding[0],binding[1],binding[2])
                else:
                    print("Incorrect type found when resetting defaults")
            else:
                self.cmd.set_key(key,'')
                delete_keys.append(key)

        for key in delete_keys:
            del self.cmd.shortcut_dict[key]
            del self.cmd.key_mappings[key]

        print("Restored default keybindings")

    def create_new_shortcut(self, new_key, new_binding):
        '''
        Creates a new shortcut after checking existing and reserved keys.
        '''
        if new_key not in self.cmd.shortcut_dict and new_key not in self.reserved_keys:
            try:
                self.cmd.set_key(new_key,new_binding)
                print('Assigning ',new_key, ' to ',new_binding)
            except Exception:
                print("This cannot be bound.")
            else:
                self.cmd.shortcut_dict.update({new_key:['','',new_binding]})
        elif new_key in self.cmd.shortcut_dict:
            try:
                self.cmd.set_key(new_key,new_binding)
                print('Assigning ',new_key, ' to ',new_binding)
            except Exception:
                print("This cannot be bound.")
            else:
                self.cmd.shortcut_dict[new_key][2] = new_binding
        else:
            print("This key is reserved.")
