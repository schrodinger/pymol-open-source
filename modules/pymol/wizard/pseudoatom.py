from pymol.wizard import Wizard

from pymol import cmd
import pymol
import types

class Pseudoatom(Wizard):

    def __init__(self,mode='label',pos='[0.0,0.0,0.0]',_self=cmd):
        Wizard.__init__(self,_self)
        self.mode = mode
        if mode == 'label':
            self.prefix = 'Label text: \888'
        self.text = ''
        self.pos = pos

    def get_event_mask(self):
        return Wizard.event_mask_key

    def do_key(self,k,x,y,m):
        if k in [8,127]:
            self.text = self.text[:-1]
        elif k==27:
            self.cmd.set_wizard()
        elif k==32:
            self.text = self.text + " "
        elif k>32:
            self.text = self.text + chr(k)
        elif k==10 or k==13:
            self.text = self.text.strip()
            if self.mode=='label':
                obj_name = self.cmd.get_unused_name(self.text[0:14].lower(),0)
                self.cmd.pseudoatom(obj_name,pos=self.pos,label=self.text)
            self.cmd.set_wizard()
        self.cmd.refresh_wizard()
        return 1

    def get_prompt(self):
        self.prompt = [ self.prefix + self.text + "_" ]
        return self.prompt

    def get_panel(self):
        return [
            [ 2, 'Cancel', 'cmd.set_wizard()' ]
            ]
