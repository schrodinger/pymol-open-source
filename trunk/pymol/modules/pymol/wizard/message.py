from pymol.wizard import Wizard
from pymol import cmd
import pymol
import types

class Message(Wizard):

    def __init__(self,*arg,**kw):
        _self = kw.get('_self',cmd)
        Wizard.__init__(self,_self)        
        self.message = []
        for a in arg:
            if not isinstance(a,types.ListType):
                self.message.append(a)
            else:
                self.message.extend(a)
        for a in self.message:
            print " "+a
        self.dismiss = int(kw.get("dismiss",1))

    def get_prompt(self):
        self.prompt = self.message
        return self.prompt

    def get_panel(self):
        if not hasattr(self,'dismiss'):
            self.dismiss=1
        if self.dismiss==1:
            return [
                [ 1, 'Message', '' ],
                [ 2, 'Dismiss', 'cmd.set_wizard()' ]
                ]
        else:
            return []





