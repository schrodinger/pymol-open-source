from pymol.wizard import Wizard
from pymol import cmd
import pymol
import types

class Toggle(Wizard):

    def __init__(self,*arg,**kw):
        self.message = []
        for a in arg:
            if not isinstance(a,types.ListType):
                self.message.append(a)
            else:
                self.message.extend(a)
        for a in self.message:
            print a
        self.message_visible = 1
            
    def toggle(self):
        self.message_visible = not self.message_visible
        cmd.refresh_wizard()

    def get_prompt(self):
        if self.message_visible:
            self.prompt = self.message
        else:
            self.prompt = None
        return self.prompt

    def get_panel(self):
        panel = [
            [ 2, 'Toggle Fullscreen', 
              'cmd.full_screen(apply(lambda x:{ "off":"on", "on":"off"}[x],(cmd.get("full_screen"),)))'],
            [ 2, 'Toggle Stereo 3D', 
              'cmd.stereo(apply(lambda x:{ "off":"on", "on":"off"}[x],(cmd.get("stereo"),)))'],
            [ 2, 'Toggle Message', 
              'cmd.get_wizard().toggle()'],
            [ 2, 'Dismiss', 'cmd.set_wizard()' ]
            ]
        if len(self.message)==0:
            del panel[2]
        return panel



