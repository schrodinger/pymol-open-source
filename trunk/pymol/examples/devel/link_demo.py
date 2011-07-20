# experimental demonstration of firing up a web browser based on
# clicks in the 3D viewer

# simply "run" this script from within PyMOL itself (file menu) or
# ./pymol link_demo.py

import webbrowser
import pymol
from pymol.wizard import Wizard
from pymol import cmd, util
import pymol

# Wizard class definition 

class Clickurl(Wizard):

    def __init__(self,_self=cmd):
        Wizard.__init__(self,_self)
        self.cmd.unpick()

    def get_prompt(self):
        self.prompt = [ 'Please click a labelled atom...' ]
        return self.prompt

    def get_panel(self):
        return [
            [ 2, 'Done','cmd.set_wizard()'],
            ]
    
    def do_select(self,name):
        pymol.stored.link = ''
        cmd.iterate(name,"stored.link=text_type") 
        if len(pymol.stored.link):
            webbrowser.open(pymol.stored.link)
        cmd.delete(name)

# load an example structure

cmd.load("$TUT/1hpv.pdb")

# color by chain (aesthetics)

util.cbc()

# store the links as atom text_types

cmd.alter("name ca",r"text_type='http://delsci.info/cgi-bin/click.cgi?residue=%s%s%s'%(resn,resi,chain)")

# put the mouse into single-atom selection mode

cmd.set('mouse_selection_mode',0)

# just show ribbon (means we can only select labelled C-alphas)

cmd.show_as("cartoon") 

# set up the labels

cmd.label("name ca","'Link'") 

# color the labels white

cmd.set("label_color", 'white')

# activate the wizard

cmd.set_wizard(Clickurl()) 

