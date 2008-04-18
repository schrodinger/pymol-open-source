
from pymol import cmd

class Wizard:

    event_mask_pick    = 1
    event_mask_select  = 2
    event_mask_key     = 4
    event_mask_special = 8 
    event_mask_scene   = 16 # scene changed
    event_mask_state   = 32 # state changed
    event_mask_frame   = 64 # frame changed
    event_mask_dirty   = 128 # anything changed (via OrthoDirty)
    
    def __init__(self,_self=cmd):
        self.menu = {}
        self.prompt = None
        self.panel = None
        self.cmd = _self
        self._validate_instance() 
        
    def _validate_instance(self):
        _pymol = self.cmd._pymol
        if not hasattr(_pymol.session, 'wizard_storage'):
            _pymol.session.wizard_storage = {}
        wizard_storage = _pymol.session.wizard_storage
        class_key = str(self.__class__)
        self.session = wizard_storage.get(class_key,{})
        # create a dictionary for per-class storage in session
        if not self.session:
            wizard_storage[class_key] = self.session
        
    def get_prompt(self):
        return self.prompt

    def get_panel(self):
        return self.panel

    def get_event_mask(self):
        return Wizard.event_mask_pick + Wizard.event_mask_select

    def do_scene(self):
        return None

    def do_state(self,state):
        return None

    def do_frame(self,frame):
        return None

    def do_dirty(self):
        return None
    
    def do_pick(self,bondFlag):
        return None

    def do_select(self,name):
        return None

    def do_key(self,k,x,y,mod):
        return None

    def do_special(self,k,x,y,mod):
        return None

    def cleanup(self):
        pass

    def get_menu(self,tag):
        result = None
        if self.menu.has_key(tag):
            result = self.menu[tag]
        return result

