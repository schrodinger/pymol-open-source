
from pymol.wizard import Wizard
from pymol import cmd
from pymol import util

import pymol

sele_pre = "sclpt_wz_"
sele_prefix_len = len(sele_pre)

cent_sele = sele_pre+ "cent"
free_sele = sele_pre+ "free"
fix_sele = sele_pre+ "fix"
excl_sele = sele_pre+ "excl"

default_mode = 'by_resi'
default_radius = 6.0
default_cushion = 6.0

# wizard status

HAVE_SELECTIONS = 1
NO_SELECTIONS = 0

class Sculpting(Wizard):

    cutoff = 3.5
    
    def __init__(self,_self=cmd):

        Wizard.__init__(self,_self)
        
        self.status = NO_SELECTIONS
        self.error = None
        self.object_name = None
        self.radius = default_radius
        self.cushion = default_cushion

        # activate editing mode if not already active

        if "Editing" not in self.cmd.get("button_mode_name"):
            self.cmd.edit_mode(1)
            self.restore_edit_mode = 0
        else:
            self.restore_edit_mode = None

        # activate bumps

        self.restore_sculpt_vdw_vis_mode = self.cmd.get('sculpt_vdw_vis_mode')
        self.cmd.set("sculpt_vdw_vis_mode")

        # deactivate all objects

        self.cmd.sculpt_deactivate("all")

        # turn on sculpting

        self.cmd.set("sculpting")

        # unmask all atoms so that they can be picked

        self.cmd.unmask("all")

        # mode selection subsystem
        
        self.mode = default_mode
        self.modes = [
            'ligand_rx',
#         'ligand_re',         
#         'by_atom',
            'by_resi',
            ]

        self.mode_name = {
            'by_atom':'Atom Shells',
            'by_resi':'Residue Shells',
            'ligand_rx' :'One Residue',
            'ligand_re' :'Residue v. Free',                  
            'ligand_cx' :'Chain v. Fixed',
            'ligand_ce' :'Chain v. Free',                  
            }

        smm = []
        smm.append([ 2, 'Selection Mode', '' ])
        for a in self.modes:
            smm.append([ 1, self.mode_name[a], 'cmd.get_wizard().set_mode("'+a+'")'])
        self.menu['mode']=smm

        self.menu['radius'] = [[ 2, 'Mobile Atom Sphere', '' ],
                                      [1, '4.0 A Radius','cmd.get_wizard().set_radius(4)'],
                                      [1, '5.0 A Radius','cmd.get_wizard().set_radius(5)'],
                                      [1, '6.0 A Radius','cmd.get_wizard().set_radius(6)'],
                                      [1, '8.0 A Radius','cmd.get_wizard().set_radius(8)'],
                                      [1, '10.0 A Radius','cmd.get_wizard().set_radius(10)'],
                                      [1, '15.0 A Radius','cmd.get_wizard().set_radius(15)'],
                                      [1, '20.0 A Radius','cmd.get_wizard().set_radius(20)'],
                                      ]
                                        
        self.menu['cushion'] = [[ 2, 'Fixed Atom Cushion', '' ],
            [1, '2.0 A Cushion','cmd.get_wizard().set_cushion(2)'],
            [1, '3.0 A Cushion','cmd.get_wizard().set_cushion(3)'],
            [1, '4.0 A Cushion','cmd.get_wizard().set_cushion(4)'],
            [1, '6.0 A Cushion','cmd.get_wizard().set_cushion(6)'],
            [1, '8.0 A Cushion','cmd.get_wizard().set_cushion(8)'],
            [1, '10.0 A Cushion','cmd.get_wizard().set_cushion(10)'],
            [1, '12.0 A Cushion','cmd.get_wizard().set_cushion(12)'],                              
            ]

# generic set routines

    def set_mode(self,mode):
        if mode in self.modes:
            self.mode = mode
        self.cmd.refresh_wizard()

    def set_radius(self,radius):
        self.radius=radius
        self.update_selections()
        self.cmd.refresh_wizard()
        
    def set_cushion(self,cushion):
        self.cushion=cushion
        self.update_selections()
        self.cmd.refresh_wizard()
        
    def update_selections(self):
        if self.status == HAVE_SELECTIONS:
            if self.mode=='by_resi':
                self.cmd.select(free_sele,"((byobj %s) and byres ((%s) x; %8.3f))"%
                              (cent_sele,cent_sele,self.radius))
                self.cmd.select(fix_sele,"((byobj %s) and (not %s) and byres ((%s) x; %8.3f))"%
                              (cent_sele,free_sele,cent_sele,self.radius+self.cushion))
                self.cmd.select(excl_sele,"((byobj %s) and (not (%s|%s)))"%
                              (cent_sele,free_sele,fix_sele))
            elif self.mode=='ligand_rx':
                self.cmd.select(free_sele,"((byobj %s) and byres (%s))"%
                              (cent_sele,cent_sele))
                self.cmd.select(fix_sele,"((byobj %s) and (not %s) and byres ((%s) x; %8.3f))"%
                              (cent_sele,free_sele,free_sele,self.cushion))
                self.cmd.select(excl_sele,"((byobj %s) and (not (%s|%s)))"%
                              (cent_sele,free_sele,fix_sele))
            self.cmd.protect("(byobj %s)"%cent_sele)
            self.cmd.deprotect(free_sele)
            self.cmd.flag('exclude',"(byobj %s)"%cent_sele,"clear")
            self.cmd.flag('exclude',excl_sele,"set")
            self.cmd.color('grey',excl_sele)
            self.cmd.unmask("(byobj %s)"%cent_sele)
            self.cmd.mask(excl_sele)
            self.cmd.mask(fix_sele)
            self.cmd.zoom(free_sele,3,animate=1)
            util.cbac(fix_sele,_self=self.cmd)
            util.cbag(free_sele,_self=self.cmd)
            self.cmd.disable('indicate')
            self.cmd.disable(cent_sele)            
            self.cmd.disable(free_sele)
            self.cmd.disable(fix_sele)
            self.cmd.disable(excl_sele)
            self.cmd.unpick()
            for obj in self.cmd.get_names(selection=cent_sele):
                self.cmd.push_undo(obj)
                self.cmd.sculpt_activate(obj)
                        
    def set_object_mode(self,mode):
        if mode in self.object_modes:
            self.object_mode = mode
        self.status = NO_SELECTIONS
        self.cmd.refresh_wizard()
        
    def get_panel(self):
        return [
            [ 1, 'Sculpting',''],
            [ 3, self.mode_name[self.mode],'mode'],
            [ 3, "Radius: %3.1f A"%self.radius,'radius'],
            [ 3, "Cushion: %3.1f A"%self.cushion,'cushion'],
            [ 2, 'Toggle Sculpting', 'cmd.set("sculpting",{"off":"on"}.get(cmd.get("sculpting"),0))' ],
            [ 2, 'Toggle Bumps', 'cmd.set("sculpt_vdw_vis_mode",{0:1}.get(int(cmd.get("sculpt_vdw_vis_mode")),0))'],
#            [ 2, 'Undo Last Change', 'cmd.undo()'],
            [ 2, 'Relocate','cmd.get_wizard().free_all()'],
            [ 2, 'Done','cmd.set_wizard()'],
            ]

    def free_all(self):
        self.clear()
        
    def clear(self):
        self.cmd.unmask("all")
        self.cmd.deprotect("all")
        self.cmd.sculpt_deactivate("all")
        if self.status == HAVE_SELECTIONS:
            util.cbag("(byobj %s)"%cent_sele, _self=self.cmd)
            self.cmd.flag("exclude","(byobj %s)"%cent_sele,"clear")
        self.status = NO_SELECTIONS
        self.cmd.delete(sele_pre+"*")
        self.cmd.refresh_wizard()
        
    def cleanup(self):
        global default_mode, default_radius, default_cushion
        default_mode = self.mode
        default_radius = self.radius
        default_cushion = self.cushion
        if self.restore_edit_mode != None:
            self.cmd.edit_mode(self.restore_edit_mode)
        self.cmd.set("sculpt_vdw_vis_mode",self.restore_sculpt_vdw_vis_mode)
        self.cmd.set("sculpting",0)
        self.clear()
        
    def get_prompt(self):
        self.prompt = None
        if self.status == NO_SELECTIONS:
            self.prompt = [ 'Please pick the center atom...']
        if self.error!=None:
            self.prompt.append(self.error)
        return self.prompt
    
    def do_pick(self,bondFlag):
        global dist_count
        if self.status == NO_SELECTIONS:
            if self.cmd.select(cent_sele,"pk1"):
                self.status = HAVE_SELECTIONS
                # save current coordinates 
                for obj in self.cmd.get_names(selection=cent_sele):
                    self.cmd.push_undo(obj)
            self.update_selections()
            self.cmd.refresh_wizard()
            return 1
        else:
            return 0
        
