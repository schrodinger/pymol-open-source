


class Wizard:

    event_mask_pick    = 1
    event_mask_select  = 2
    event_mask_key     = 4
    event_mask_special = 8
    event_mask_scene   = 16
    event_mask_state   = 32
    event_mask_frame   = 64
    
    def __init__(self):
        self.menu = {}
        self.prompt = None
        self.panel = None
        
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

