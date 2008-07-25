from pymol.wizard import Wizard

from pymol import cmd
import pymol
import types
import string

class Renaming(Wizard):

    def __init__(self,old_name,mode='object',_self=cmd):
        Wizard.__init__(self,_self)
        
        self.prefix = 'Renaming \\999%s\\--- to: \\999'%old_name
        self.old_name = old_name
        self.new_name = old_name
        self.mode = mode
        
    def get_event_mask(self):
        return Wizard.event_mask_key

    def do_key(self,k,x,y,m):
        if k in [8,127]:
            self.new_name = self.new_name[:-1]
        elif k==27:
            cmd.set_wizard()
            cmd.refresh()
        elif k==32:
            self.new_name = self.new_name + "_"
        elif k>32:
            self.new_name= self.new_name + chr(k)
        elif k==10 or k==13:
            self.new_name = string.strip(self.new_name)
            if self.mode=='object':
                cmd.do("set_name %s,%s"%(self.old_name,self.new_name),log=0)
            elif self.mode=='scene':
                cmd.do("scene %s,rename,new_key=%s"%(self.old_name,self.new_name),log=0)                
            cmd.set_wizard()
            cmd.refresh()
            return 1
        cmd.refresh_wizard()
        return 1
        
    def get_prompt(self):
        self.prompt = [ self.prefix + self.new_name + "_" ]
        return self.prompt

    def get_panel(self):
        return [
            [ 1, 'Renaming', '' ],
            [ 2, 'Cancel', 'cmd.set_wizard()' ]
            ]


