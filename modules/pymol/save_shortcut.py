import sys
import os
import threading
import json

_SHORTCUTS_SAVE_FILE = u'~/.pymol/shortcuts_save.json'

from os.path import expanduser, expandvars

def get_shortcut_save_filename():
    '''
    Returns the users shortcuts_save filename.
    '''
    filename = expandvars(_SHORTCUTS_SAVE_FILE)
    filename = expanduser(filename)
    return filename

def save_shortcuts(shortcuts_dict):
    '''
    Saves the contents of shortcuts_dict to the filename supplied by
    get_shortcut_save_filename().
    '''
    save_file = get_shortcut_save_filename()

    try:
        os.makedirs(os.path.dirname(save_file), 0o750)
    except OSError:
        # This will trigger if the directory has already been made for license.
        pass

    try:
        with open(save_file, 'w') as savefile:
            json.dump(shortcuts_dict, savefile)
            print("Saved shortcuts to file ", save_file)
    except Exception as e:
        print("Unable to save to file.")

def load_shortcuts_dict():
    '''
    Attempts to load the filename from get_shortcut_save_filename.
    Quietly passes file doesn't exist (i.e. nothing has been saved yet)
    Dictionary is returned through save_dict.
    '''
    save_dict = {}
    save_file = get_shortcut_save_filename()
    try:
        with open(save_file) as save_dict_file:
            save_dict = json.load(save_dict_file)
    except FileNotFoundError as file_error:
        # This represents the case where no save file has been created yet and should pass quietly. 
        pass
    except Exception as e:
        print("No shortcut save file has been loaded.")
        print(e)
    return save_dict

def setkey_from_dict(save_dict,cmd):
    '''
    Iterates over a dictionary and sets the corresponding keys with cmd.set_key.
    '''
    for key in save_dict:
        if save_dict[key][2]:
            cmd.set_key(key,save_dict[key][2])

def load_and_set(cmd):
    '''
    Loads save_dict from file, sets using setkey_from_dict. Returns save_dict.
    '''
    save_dict = load_shortcuts_dict()
    setkey_from_dict(save_dict,cmd)
    return save_dict
