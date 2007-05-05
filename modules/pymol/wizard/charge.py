
from pymol.wizard import Wizard
from pymol import cmd
import pymol

default_mode = 'labchg'

class Charge(Wizard):

    def __init__(self,_self=cmd):
        Wizard.__init__(self,_self)
        
        self.modes = [
            'labchg',
            'addchg',
            'cpychg',
            'zrochg',
            'mzochg',
            'cbachg',
            'movchg',
            'rbachg',
            'sumchg',
            ]

        self.mode = default_mode
        self.status = 0
        
        self.mode_name = {
            'labchg':'Show',
            'cpychg':'Copy',
            'addchg':'Add',
            'zrochg':'Zero',
            'mzochg':'Move & Zero (atom)',
            'cbachg':'Move & Zero (resi)',
            'movchg':'Move & Remove (atom)',
            'rbachg':'Move & Remove (resi)',
            'sumchg':'Get Total Charge'
            }
        
        # initialize mode menu
        
        smm = []
        smm.append([ 2, 'Atom Charge Mode', '' ])
        for a in self.modes:
            smm.append([ 1, self.mode_name[a], 'cmd.get_wizard().set_mode("'+a+'")'])

        self.menu['mode']=smm
        
        self.memory = 0
            

    def get_panel(self):
        return [
            [ 1, 'Charge Wizard',''],
            [ 3, self.mode_name[self.mode],'mode'],
            [ 2, 'Clear','cmd.get_wizard().clear()'],
            [ 2, 'Done','cmd.set_wizard()'],
            ]

    def cleanup(self):
        global default_mode
        default_mode = self.mode
        self.clear()
        
    def clear(self):
        self.set_status(0)
        if 'wcharge' in self.cmd.get_names('selections'):
            if self.mode!='sumchg':
                self.cmd.edit("wcharge")
                self.cmd.label("pkmol",'') # fastest clear command
            else:
                self.cmd.label("wcharge",'') # fastest clear command            
            self.cmd.delete("wcharge")
            self.cmd.unpick()
        self.cmd.unpick()
        self.cmd.refresh_wizard()
        
    def get_prompt(self):
        self.prompt = None
        if self.mode == 'cpychg':
            if self.status==0:
                self.prompt = [ 'Pick source atom...' ]
            elif self.status==1:
                self.prompt = [ 'Pick destination atom on which to assign charge %6.4f'%self.partial_charge ]
        if self.mode == 'addchg':
            if self.status==0:
                self.prompt = [ 'Pick source atom...' ]
            elif self.status==1:
                self.prompt = [ 'Pick destination atom on which to add charge %6.4f'%self.partial_charge ]
        if self.mode == 'mzochg':
            if self.status==0:
                self.prompt = [ 'Pick source atom to copy and zero...' ]
            elif self.status==1:
                self.prompt = [ 'Pick destination atom on which to add charge %6.4f'%self.partial_charge ]
        if self.mode == 'movchg':
            if self.status==0:
                self.prompt = [ 'Pick source atom to copy and destroy...' ]
            elif self.status==1:
                self.prompt = [ 'Pick destination atom on which to add charge %6.4f'%self.partial_charge ]
        if self.mode == 'sumchg':
            if self.status==0:
                self.prompt = [ 'Pick an atom on the chain...' ]
            if self.status==1:
                self.prompt = [ 'Total charge on the chain is %6.4f'%self.partial_charge,
                                     'Pick an atom on the chain...' ]
        if self.mode == 'cbachg':
            if self.status==0:
                self.prompt = [ 'Pick source residue to copy and zero...' ]
            elif self.status==1:
                self.prompt = [ 'Pick destination residue on which to add charges.']
        if self.mode == 'rbachg':
            if self.status==0:
                self.prompt = [ 'Pick source residue to copy and remove...' ]
            elif self.status==1:
                self.prompt = [ 'Pick destination residue on which to add charges.']
        if self.mode == 'zrochg':
                self.prompt = [ 'Pick atom on which to zero charge...' ]


        if self.mode == 'labchg':
                self.prompt = [ 'Pick atom on which to show charge...' ]

        if "wcharge" in self.cmd.get_names('selections'):
            pymol.stored.charge = 0
            if self.cmd.iterate("(byres wcharge)",
                                "stored.charge = stored.charge + partial_charge"):
                self.prompt.insert(0,"Total charge on the residue is %6.4f"%pymol.stored.charge)
                
        return self.prompt
    
    def set_mode(self,mode):
        if mode in self.modes:
            self.mode = mode
        self.status = 0
        self.cmd.refresh_wizard()
        
    def set_status(self,status):
        self.status = status
        self.cmd.refresh_wizard()

    def do_pick(self,bondFlag):
        if bondFlag:
            print " Error: please select a single atom"
            
        if self.mode == 'cpychg':
            # picking up
            if self.status==0:
                if self.cmd.iterate("(pk1)","stored.charge = partial_charge"):
                    self.partial_charge = pymol.stored.charge
                    self.status = 1
                    self.cmd.label("(pk1)","'%6.4f'%partial_charge")
                    self.cmd.select("wcharge","(pk1)")
                    self.cmd.unpick()
                    self.cmd.enable("wcharge")
                    
            # dropping off
            elif self.status==1:
                pymol.stored.charge=self.partial_charge
                if self.cmd.alter("(pk1)","partial_charge = stored.charge"):
                    self.status = 0
                    self.cmd.label("(pk1)","'%6.4f'%partial_charge")
                    self.cmd.select("wcharge","(pk1)")
                    self.cmd.unpick()
                    self.cmd.enable("wcharge")

        if self.mode == 'cbachg' or self.mode == 'rbachg':
            # picking up
            if self.status==0:
                pymol.stored.chg_dict = {}
                if self.cmd.iterate("(byres pk1)","stored.chg_dict[name] = partial_charge"):
                    self.charge_dict = pymol.stored.chg_dict
                    self.status = 1
                    self.cmd.label("(byres pk1)","'%6.4f'%partial_charge")
                    self.cmd.select("wcharge","(byres pk1)")
                    self.cmd.unpick()
                    self.cmd.enable("wcharge")
                    
            # dropping off
            elif self.status==1:
                pymol.stored.valid_atoms = []
                if self.cmd.iterate("(byres pk1)","stored.valid_atoms.append(name)"):
                    kees = self.charge_dict.keys()
                    valid_dict = {}
                    for a in pymol.stored.valid_atoms:
                        if a in kees:
                            valid_dict[a] = 1
                    pymol.stored.chg_dict = self.charge_dict
                    # copy/add charges
                    for a in valid_dict.keys():
                        self.cmd.alter("((byres pk1) and name %s)"%a,
                                     "partial_charge = partial_charge + stored.chg_dict[name]")
                        if self.mode == 'rbachg':
                            self.cmd.remove("((wcharge) and name %s)"%a)
                        else:
                            self.cmd.alter("((wcharge) and name %s)"%a,
                                         "partial_charge = 0")
                    self.status = 0
                    # update labels
                    self.cmd.label("(wcharge or (byres pk1))","'%6.4f'%partial_charge")
                    # show which atoms had charges moved
                    self.cmd.select("wcharge","(none)")
                    for a in valid_dict.keys():
                        self.cmd.select("wcharge","(wcharge or ((byres pk1) and name %s))"%a)
                    self.cmd.unpick()
                    self.cmd.enable("wcharge")

        if self.mode == 'addchg':
            # picking up
            if self.status==0:
                if self.cmd.iterate("(pk1)","stored.charge = partial_charge"):
                    self.partial_charge = pymol.stored.charge
                    self.status = 1
                    self.cmd.label("(pk1)","'%6.4f'%partial_charge")
                    self.cmd.select("wcharge","(pk1)")
                    self.cmd.unpick()
                    self.cmd.enable("wcharge")
                    
            # dropping off
            elif self.status==1:
                pymol.stored.charge=self.partial_charge
                if self.cmd.alter("(pk1)","partial_charge = partial_charge + stored.charge"):
                    self.status = 0
                    self.cmd.label("(pk1)","'%6.4f'%partial_charge")
                    self.cmd.select("wcharge","(pk1)")
                    self.cmd.unpick()
                    self.cmd.enable("wcharge")
                    
        if self.mode == 'mzochg':
            # picking up
            if self.status==0:
                if self.cmd.iterate("(pk1)","stored.charge = partial_charge"):
                    self.partial_charge = pymol.stored.charge
                    self.status = 1
                    self.cmd.label("(pk1)","'%6.4f'%partial_charge")
                    self.cmd.select("wcharge","(pk1)")
                    self.cmd.unpick()
                    self.cmd.enable("wcharge")
                    
            # dropping off
            elif self.status==1:
                pymol.stored.charge=self.partial_charge
                if self.cmd.alter("(pk1)","partial_charge = partial_charge + stored.charge"):
                    self.status = 0
                    self.cmd.label("(pk1)","'%6.4f'%partial_charge")
                    self.cmd.alter("(wcharge)","partial_charge=0")
                    self.cmd.label("(wcharge)","'%6.4f'%partial_charge")               
                    self.cmd.select("wcharge","(pk1)")
                    self.cmd.unpick()
                    self.cmd.enable("wcharge")

        if self.mode == 'movchg':
            # picking up
            if self.status==0:
                if self.cmd.iterate("(pk1)","stored.charge = partial_charge"):
                    self.partial_charge = pymol.stored.charge
                    self.status = 1
                    self.cmd.label("(pk1)","'%6.4f'%partial_charge")
                    self.cmd.select("wcharge","(pk1)")
                    self.cmd.unpick()
                    self.cmd.enable("wcharge")
                    
            # dropping off
            elif self.status==1:
                pymol.stored.charge=self.partial_charge
                if self.cmd.alter("(pk1)","partial_charge = partial_charge + stored.charge"):
                    self.status = 0
                    self.cmd.label("(pk1)","'%6.4f'%partial_charge")
                    self.cmd.remove("wcharge")
                    self.cmd.select("wcharge","(pk1)")
                    self.cmd.unpick()
                    self.cmd.enable("wcharge")

        if self.mode == 'zrochg':
            if self.cmd.alter("(pk1)","partial_charge = 0.0"):
                self.cmd.label("(pk1)","'%6.4f'%partial_charge")
                self.cmd.select("wcharge","(pk1)")
                self.cmd.unpick()
                self.cmd.enable("wcharge")
                    
        if self.mode == 'labchg':
            self.cmd.label("(pk1)","'%6.4f'%partial_charge")
            self.cmd.select("wcharge","(pk1)")
            self.cmd.unpick()
            self.cmd.enable("wcharge")
                    
        if self.mode == 'sumchg':
            pymol.stored.charge = 0.0
            if self.cmd.iterate("(pkmol)","stored.charge = stored.charge + partial_charge"):
                self.partial_charge = pymol.stored.charge
                self.status = 1
                self.cmd.select("wcharge","(pkmol)")
                self.cmd.unpick()
                self.cmd.enable("wcharge")
            
        self.cmd.refresh_wizard()
