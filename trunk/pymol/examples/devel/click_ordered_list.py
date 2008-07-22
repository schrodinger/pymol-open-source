# example showing how one might assembled an ordered list of atoms
# from mouse clicked through use of a PyMOL Wizard.

from pymol.wizard import Wizard
from pymol import cmd, util
import pymol
import traceback

click_sele = "_clicked"
disp_sele = "_displayed"

# Wizard class definition 

def quote_blanks(st):
    st = str(st)
    if not len(st):
        return "''"
    else:
        return st
    
class Clicker(Wizard):

    def __init__(self,_self=cmd):
        Wizard.__init__(self,_self)
        self.cmd.unpick()
        self.click_list = []
        cmd.select(click_sele,"none",enable=0)
        cmd.select(disp_sele,"none",enable=1)

    def get_prompt(self):
        self.prompt = [ 'Please click a set of atoms one by one (%d clicked so far)...' % len(self.click_list) ]
        return self.prompt

    def get_panel(self):
        return [
            [ 2, 'Print Atom List','cmd.get_wizard().print_list()'],
            [ 2, 'Remove Last','cmd.get_wizard().remove_last()'],
            [ 2, 'Reset Atom List','cmd.get_wizard().reset_list()'],                        
            [ 2, 'Done','cmd.set_wizard()'],
            ]
    
    def do_select(self,name):
        if name!=disp_sele:
            self.cmd.select(disp_sele,name + " or ?" + disp_sele,enable=1)
        dict = { 'x' : [] }
        self.cmd.iterate(disp_sele+" and not ?"+click_sele,
                    'x.append( (model,segi,chain,resn,resi,name,alt) )',
                    space=dict)
        if len(dict['x']): # toggled on
            entry = tuple(map(quote_blanks,dict['x'][0]))
            sele_str = "/%s/%s/%s/%s`%s/%s`%s"%entry
            self.cmd.select(click_sele,"?"+click_sele+" or "+sele_str,enable=0)
            self.click_list.append(sele_str)
            self.cmd.refresh_wizard()
        else: # toggled off
            new_list = []
            for entry in self.click_list:
                if self.cmd.count_atoms(entry+" and "+disp_sele):
                    new_list.append(entry)
            self.click_list = new_list
            self.cmd.select(click_sele,disp_sele,enable=0)
            cmd.refresh_wizard()
        self.cmd.delete(name)
        self.cmd.select(disp_sele,click_sele,enable=1)            
    
    def remove_last(self):
        self.click_list = self.click_list[:-1]
        self.cmd.select(click_sele,"none")
        for entry in self.click_list:
            self.cmd.select(click_sele,click_sele+" or "+entry)
        self.cmd.select(disp_sele,click_sele,enable=1)
        self.cmd.refresh_wizard()
                                
    def reset_list(self):
        self.click_list = []
        cmd.select(click_sele,"none",enable=0)
        cmd.select(disp_sele,"none",enable=1)
        self.cmd.refresh_wizard()
        
    def print_list(self):
        cnt = 1
        print "Atoms Clicked (in order):"
        for entry in self.click_list:
            print "Atom %d: %s"%(cnt,entry)
            cnt = cnt + 1

    def cleanup(self):
        self.cmd.delete(disp_sele)        
        self.cmd.delete(click_sele)
        
# load an example structure

cmd.load("$TUT/1hpv.pdb")

# color by chain (aesthetics)

util.cbc()

# put the mouse into single-atom selection mode

cmd.set('mouse_selection_mode',0)

# activate the wizard

cmd.set_wizard(Clicker()) 

