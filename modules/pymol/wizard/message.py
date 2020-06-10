from pymol.wizard import Wizard
from pymol import cmd
import pymol
import types

import re
_nuke_color_re = re.compile(r"\\[0-9][0-9][0-9]")

class Message(Wizard):

    def __init__(self, *arg, dismiss=1, _self=cmd):
        Wizard.__init__(self,_self)
        self.message = []
        for a in arg:
            if not isinstance(a,list):
                self.message.append(a)
            else:
                self.message.extend(a)
        for a in self.message:
            print(" " + _nuke_color_re.sub('',a))
        self.dismiss = int(dismiss)

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
