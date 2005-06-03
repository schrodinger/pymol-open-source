
from pymol.wizard import Wizard
from pymol import cmd
from pymol import util

import pymol

sele_pre = "sw_"
sele_prefix_len = len(sele_pre)

cent_sele = sele_pre+ "cent"
free_sele = sele_pre+ "free"
fix_sele = sele_pre+ "fix"
excl_sele = sele_pre+ "excl"

default_mode = 'by_resi'
default_radius = 6.0
default_cushion = 6.0

class Sculpting(Wizard):

    cutoff = 3.5
    
    def __init__(self):

        Wizard.__init__(self)
        
        self.status = 0 # 0 no atoms selections, 1 atom selected
        self.error = None
        self.object_name = None
        self.radius = default_radius
        self.cushion = default_cushion
        
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
            'ligand_rx' :'Residue v. Fixed',
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
        cmd.refresh_wizard()

    def set_radius(self,radius):
        self.radius=radius
        self.update_selections()
        cmd.refresh_wizard()
        
    def set_cushion(self,cushion):
        self.cushion=cushion
        self.update_selections()
        cmd.refresh_wizard()
        
    def update_selections(self):
        if self.status==1:
            if self.mode=='by_resi':
                cmd.select(free_sele,"((byobj %s) and byres ((%s) x; %8.3f))"%
                              (cent_sele,cent_sele,self.radius))
                cmd.select(fix_sele,"((byobj %s) and (not %s) and byres ((%s) x; %8.3f))"%
                              (cent_sele,free_sele,cent_sele,self.radius+self.cushion))
                cmd.select(excl_sele,"((byobj %s) and (not (%s|%s)))"%
                              (cent_sele,free_sele,fix_sele))
            elif self.mode=='ligand_rx':
                cmd.select(free_sele,"((byobj %s) and byres (%s))"%
                              (cent_sele,cent_sele))
                cmd.select(fix_sele,"((byobj %s) and (not %s) and byres ((%s) x; %8.3f))"%
                              (cent_sele,free_sele,free_sele,self.cushion))
                cmd.select(excl_sele,"((byobj %s) and (not (%s|%s)))"%
                              (cent_sele,free_sele,fix_sele))
            cmd.protect("(byobj %s)"%cent_sele)
            cmd.deprotect(free_sele)
            cmd.flag('exclude',"(byobj %s)"%cent_sele,"clear")
            cmd.flag('exclude',excl_sele,"set")
            cmd.color('grey',excl_sele)
            cmd.zoom(free_sele,3.0)
            util.cbac(fix_sele)
            util.cbag(free_sele)
            cmd.disable('indicate')
            cmd.disable(cent_sele)            
            cmd.disable(free_sele)
            cmd.disable(fix_sele)
            cmd.disable(excl_sele)
            cmd.unpick()
            
            self.status = 1
            
    def set_object_mode(self,mode):
        if mode in self.object_modes:
            self.object_mode = mode
        self.status = 0
        cmd.refresh_wizard()
        
    def get_panel(self):
        return [
            [ 1, 'Sculpting',''],
            [ 3, self.mode_name[self.mode],'mode'],
            [ 3, "Radius: %3.1f A"%self.radius,'radius'],
            [ 3, "Cushion: %3.1f A"%self.cushion,'cushion'],
            [ 2, 'Relocate','cmd.get_wizard().free_all()'],
            [ 2, 'Done','cmd.set_wizard()'],
            ]

    def free_all(self):
        if self.status==1:
            cmd.deprotect("(byobj %s)"%cent_sele)
            cmd.protect(fix_sele)
            cmd.flag('exclude',"(byobj %s)"%cent_sele,"clear")
            util.cbag("(byobj %s)"%cent_sele)
        self.clear()
        cmd.refresh_wizard()
        
    def clear(self):
        cmd.delete(sele_pre+"*")
        self.status = 0
        cmd.refresh_wizard()
        
    def cleanup(self):
        global default_mode, default_radius, default_cushion
        default_mode = self.mode
        default_radius = self.radius
        default_cushion = self.cushion
        self.clear()
        
    def get_prompt(self):
        self.prompt = None
        if self.status==0:
            self.prompt = [ 'Please pick the center atom...']
        if self.error!=None:
            self.prompt.append(self.error)
        return self.prompt
    
    def do_pick(self,bondFlag):
        global dist_count
        if not self.status:
            if cmd.select(cent_sele,"pk1"):
                self.status = 1
                self.update_selections()
            cmd.refresh_wizard()
            return 1
        else:
            return 0
        
