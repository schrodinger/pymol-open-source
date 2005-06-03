# Warren L. DeLano

from pymol.wizard import Wizard
from pymol import cmd
import pymol

saved_mode = 2 # toggle
saved_scope = 1 # byres
saved_what = 3 # sticks
saved_color = 5 # magenta

class Appearance(Wizard):


    mode_dict = {
        0 : ['Color', '_ cmd.color' ],
        1 : ['Color (elem c)', '_ util.color_carbon' ],
        2 : ['Toggle','_ cmd.toggle'],
        3 : ['Show', '_ cmd.show'],
        4 : ['Hide', '_ cmd.hide'],
        5 : ['Select', '_ cmd.select'],
        }

    what_dict = {
        1 : ['Lines', 'lines'],
        2 : ['Nonbonded', 'nonbonded' ],
        3 : ['Sticks', 'sticks' ],
        4 : ['Ribbon', 'ribbon' ],
        5 : ['Cartoon', 'cartoon' ],
        6 : ['Labels', 'labels' ],
        7 : ['Dots', 'dots' ],
        8 : ['Spheres', 'spheres' ],
        9 : ['NB Spheres', 'nb_spheres' ],
        10: ['Mesh', 'mesh' ],
        11: ['Surface', 'surface' ],
        }

    color_dict = {
        1 : ['\\900red', 'red' ],
        2 : ['\\090green', 'green' ],
        3 : ['\\009blue'  , 'blue' ],
        4 : ['\\990yellow' , 'yellow' ],
        5 : ['\\909magenta', 'magenta' ],
        6 : ['\\099cyan'    ,  'cyan' ],
        7 : ['\\955salmon'   , 'salmon' ],
        8 : ['\\595lime'      , 'lime' ],
        9 : ['\\967pink'  , 'pink' ],
        10 : ['\\559slate' ,  'slate' ],
        11 : ['\\949violet' , 'violet' ],    
        12 : ['\\950orange'  , 'orange' ],
        13 : ['\\059marine'   ,  'marine' ],
        14 : ['\\905hotpink' , 'hotpink' ],
        }

    scope_dict = {
        0: ['By Atom', ''],
        1: ['By Residue', 'byres' ],
        2: ['By Chain', 'bychain' ],
        3: ['By Segment', 'bysegment' ],
        4: ['By Object', 'byobject' ],
        5: ['By Molecule', 'bymol' ],
        }

    def __init__(self):

        cmd.deselect()
        cmd.unpick()
        Wizard.__init__(self)
        self.selection_mode = cmd.get_setting_legacy("mouse_selection_mode")
        cmd.set("mouse_selection_mode",0) # set selection mode to atomic      
        self.current_mode = saved_mode
        self.current_what = saved_what
        self.current_scope = saved_scope
        self.current_color = saved_color
        self.menu['mode'] = [
            [2, 'Mode', '' ],
            [1, self.mode_dict[0][0], 'cmd.get_wizard().set_mode(0)' ],
            [1, self.mode_dict[1][0], 'cmd.get_wizard().set_mode(1)' ],
            [ 0, ''           , '' ],         
            [1, self.mode_dict[2][0], 'cmd.get_wizard().set_mode(2)' ],
            [1, self.mode_dict[3][0], 'cmd.get_wizard().set_mode(3)' ],
            [1, self.mode_dict[4][0], 'cmd.get_wizard().set_mode(4)' ],
#         [1, self.mode_dict[5][0], 'cmd.get_wizard().set_mode(4)' ],         
            ]

        self.menu['what'] = [
            [2, 'What', '' ],
            [1, self.what_dict[1][0], 'cmd.get_wizard().set_what(1)' ],
            [1, self.what_dict[2][0], 'cmd.get_wizard().set_what(2)' ],
            [1, self.what_dict[3][0], 'cmd.get_wizard().set_what(3)' ],
            [1, self.what_dict[4][0], 'cmd.get_wizard().set_what(4)' ],
            [1, self.what_dict[5][0], 'cmd.get_wizard().set_what(5)' ],
            [0, '', '' ],
            [1, self.what_dict[6][0], 'cmd.get_wizard().set_what(6)' ],
            [0, '', '' ],         
            [1, self.what_dict[7][0], 'cmd.get_wizard().set_what(7)' ],
            [1, self.what_dict[8][0], 'cmd.get_wizard().set_what(8)' ],
            [1, self.what_dict[9][0], 'cmd.get_wizard().set_what(9)' ],
            [0, '', '' ],
            [1, self.what_dict[10][0], 'cmd.get_wizard().set_what(10)' ],
            [1, self.what_dict[11][0], 'cmd.get_wizard().set_what(11)' ],
            ]

        self.menu['color'] = [
            [2, 'Color', '' ],
            [1, self.color_dict[1][0], 'cmd.get_wizard().set_color(1)' ],
            [1, self.color_dict[2][0], 'cmd.get_wizard().set_color(2)' ],
            [1, self.color_dict[3][0], 'cmd.get_wizard().set_color(3)' ],
            [1, self.color_dict[4][0], 'cmd.get_wizard().set_color(4)' ],
            [1, self.color_dict[5][0], 'cmd.get_wizard().set_color(5)' ],
            [1, self.color_dict[6][0], 'cmd.get_wizard().set_color(6)' ],
            [1, self.color_dict[7][0], 'cmd.get_wizard().set_color(7)' ],
            [1, self.color_dict[8][0], 'cmd.get_wizard().set_color(8)' ],
            [1, self.color_dict[9][0], 'cmd.get_wizard().set_color(9)' ],
            [1, self.color_dict[10][0], 'cmd.get_wizard().set_color(10)' ],
            [1, self.color_dict[11][0], 'cmd.get_wizard().set_color(11)' ],
            [1, self.color_dict[12][0], 'cmd.get_wizard().set_color(12)' ],
            [1, self.color_dict[13][0], 'cmd.get_wizard().set_color(13)' ],
            [1, self.color_dict[14][0], 'cmd.get_wizard().set_color(14)' ],         
            ]

        self.menu['scope'] = [
            [ 2, 'Scope', ''],
            [ 1, self.scope_dict[0][0], 'cmd.get_wizard().set_scope(0)' ],
            [ 1, self.scope_dict[1][0], 'cmd.get_wizard().set_scope(1)' ],
            [ 1, self.scope_dict[2][0], 'cmd.get_wizard().set_scope(2)' ],
            [ 1, self.scope_dict[3][0], 'cmd.get_wizard().set_scope(3)' ],
            [ 1, self.scope_dict[4][0], 'cmd.get_wizard().set_scope(4)' ],        
            [ 0, ''           , '' ],
            [ 1, self.scope_dict[5][0], 'cmd.get_wizard().set_scope(5)' ],                 
            ]

    def set_scope(self,scope):
        scope = int(scope)
        if scope in self.scope_dict:
            self.current_scope = scope
        cmd.refresh_wizard()

    def set_what(self,what):
        what = int(what)
        if what in self.what_dict:
            self.current_what = what
        cmd.refresh_wizard()

    def set_color(self,color):
        color = int(color)
        if color in self.color_dict:
            self.current_color = color
        cmd.refresh_wizard()

    def set_mode(self,mode):
        mode = int(mode)
        if mode in self.mode_dict:
            self.current_mode = mode
        cmd.refresh_wizard()

    def undo(self):
        print "no undo!"
        
    def get_prompt(self):
        self.prompt = []
        return self.prompt

    def get_panel(self):
        panel = [
                [ 1, 'Appearance Wizard',''],
                [ 3, self.mode_dict[self.current_mode][0], 'mode' ],
                ]
        if self.current_mode in [0,1]: # color mode
            panel.append(
                [ 3, self.color_dict[self.current_color][0], 'color' ])
        elif self.current_mode in [2,3,4]: # what mode
            panel.append(
            [ 3, self.what_dict[self.current_what][0], 'what' ])
        else: # select mode
            panel.append(
                [ 1, 'Atoms' , '' ])
            
        panel.extend([
                [ 3, self.scope_dict[self.current_scope][0], 'scope' ],
                #         [ 2, 'Undo', 'cmd.get_wizard().undo()' ],
                [ 2, 'Done','cmd.set_wizard()'],
                ])
        return panel
            

    def do_pick(self,bondFlag):
        if self.current_mode in [0,1]: # color
            sele = "(%s pk1)"%self.scope_dict[self.current_scope][1]
            color = self.color_dict[self.current_color][1]
            mode = self.mode_dict[self.current_mode][1]
            cmmd = mode+'("%s","%s")'%(color,sele)
            cmd.do(cmmd,log=0)
        elif self.current_mode in [2,3,4]: # show/hide/toggle
            sele = "(%s pk1)"%self.scope_dict[self.current_scope][1]
            what = self.what_dict[self.current_what][1]
            mode = self.mode_dict[self.current_mode][1]
            cmmd = mode+'("%s","%s")'%(what,sele)
            cmd.do(cmmd,log=0)
        else: # select
            pass
        cmd.unpick()
        cmd.refresh_wizard()
        return 1

    def do_select(self,selection):
        if self.current_mode in [0,1]: # color
            sele = "(%s %s)"%(self.scope_dict[self.current_scope][1],selection)
            color = self.color_dict[self.current_color][1]
            mode = self.mode_dict[self.current_mode][1]
            cmmd = mode+'("%s","%s")'%(color,sele)
            cmd.do(cmmd,log=0)
        elif self.current_mode in [2,3,4]: # show/hide/toggle
            sele = "(%s %s)"%(self.scope_dict[self.current_scope][1],selection)
            what = self.what_dict[self.current_what][1]
            mode = self.mode_dict[self.current_mode][1]
            cmmd = mode+'("%s","%s")'%(what,sele)
            cmd.do(cmmd,log=0)
        cmd.delete(selection)
        cmd.deselect()
        cmd.unpick()
        cmd.refresh_wizard()
        return 1

    def cleanup(self):
        cmd.set("mouse_selection_mode",self.selection_mode) # restore selection mode      
        global saved_scope
        saved_scope = self.current_scope
        global saved_mode
        saved_mode = self.current_mode
        global saved_color
        saved_color = self.current_color
        global saved_what
        saved_what = self.current_what
