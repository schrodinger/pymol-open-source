from pymol.wizard import Wizard
from pymol import cmd
import pymol
import types

class Fedora(Wizard):

    def __init__(self,*arg):
        self.message = []
        for a in arg:
            if not isinstance(a,types.ListType):
                self.message.append(a)
            else:
                self.message.extend(a)
            
    def get_prompt(self):
        self.prompt = self.message
        return self.prompt

    def do_pick(self,bondFlag):
        if(bondFlag):
            print " " # clear out misleading bond pick pk2 information...
        cmd.unpick()

    def do_select(self,name):
#      cmd.deselect()
        pass
        
    def get_panel(self):
        return [
            [ 1, 'Message', '' ],
            [ 2, 'Dismiss', 'cmd.set_wizard()' ]
            ]



