from pymol.wizard import Wizard
from pymol import cmd
import pymol

import traceback

sele_prefix = "_dw"
sele_prefix_len = len(sele_prefix)

dist_prefix = "dist"

default_mode = 'pairs'
default_object_mode  = 'append'
dist_count = 0

class Distance(Wizard):

    cutoff = 3.5

    def __init__(self,_self=cmd):

        cmd.unpick();
        Wizard.__init__(self,_self)

        self.status = 0 # 0 no atoms selections, 1 atom selected
        self.error = None
        self.object_name = None

        # mode selection subsystem

        self.mode = default_mode
        self.modes = [
            'polar',
            'heavy',
            'neigh',
            'pairs',
            ]

        self.mode_name = {
            'polar':'Polar Neighbors',
            'heavy':'Heavy Neighbors',
            'neigh':'Neighbors',
            'pairs':'Pairwise Distances',
            }

        smm = []
        smm.append([ 2, 'Measurement Mode', '' ])
        for a in self.modes:
            smm.append([ 1, self.mode_name[a], 'cmd.get_wizard().set_mode("'+a+'")'])
        self.menu['mode']=smm

        # overwrite mode selection subsystem

        self.object_mode=default_object_mode
        self.object_modes = [
            'overwr',
            'append',
            ]
        self.object_mode_name = {
            'overwr':'Replace Previous',
            'append':'Create New',
            }

        smm = []
        smm.append([ 2, 'New Distances?', '' ])
        for a in self.object_modes:
            smm.append([ 1, self.object_mode_name[a], 'cmd.get_wizard().set_object_mode("'+a+'")'])
        self.menu['object_mode']=smm
        self.selection_mode = cmd.get_setting_int("mouse_selection_mode")
        cmd.set("mouse_selection_mode",0) # set selection mode to atomic
        cmd.deselect() # disable the active selection (if any)

# generic set routines

    def set_mode(self,mode):
        if mode in self.modes:
            self.mode = mode
        self.status = 0
        self.clear()
        cmd.refresh_wizard()

    def set_object_mode(self,mode):
        if mode in self.object_modes:
            self.object_mode = mode
        self.status = 0
        cmd.refresh_wizard()


    def get_panel(self):
        return [
            [ 1, 'Distance Measurement',''],
            [ 3, self.mode_name[self.mode],'mode'],
            [ 3, self.object_mode_name[self.object_mode],'object_mode'],
            [ 2, 'Delete Last' , 'cmd.get_wizard().delete_last()'],
            [ 2, 'Delete All' , 'cmd.get_wizard().delete_all()'],
            [ 2, 'Done','cmd.set_wizard()'],
            ]

    def cleanup(self):
        global default_mode, default_object_mode
        default_mode = self.mode
        default_object_mode = self.object_mode
        self.clear()
        cmd.set("mouse_selection_mode",self.selection_mode) # restore selection mode

    def clear(self):
        cmd.delete(sele_prefix+"*")

    def get_prompt(self):
        self.prompt = None
        if self.mode == 'pairs':
            if self.status==0:
                self.prompt = [ 'Please click on the first atom...']
            elif self.status==1:
                self.prompt = [ 'Please click on the second atom...' ]
        elif self.mode in [ 'polar', 'neigh' ]:
            self.prompt = [ 'Please click an atom...']
        if self.error is not None:
            self.prompt.append(self.error)
        return self.prompt

    def delete_last(self):
        global dist_count
        if self.status==0:
            if dist_count>0:
                cmd.delete(dist_prefix+"%02d"%dist_count)
                dist_count = dist_count - 1
        self.status=0
        self.error = None
        self.clear()
        cmd.refresh_wizard()

    def delete_all(self):
        global dist_count
        dist_count = 0
        cmd.delete(dist_prefix+"*")
        self.status=0
        self.error = None
        self.clear()
        cmd.refresh_wizard()

    def do_select(self,name): # map selects into picks
        cmd.unpick()
        cmd.edit(name)
        cmd.delete(name)
        self.do_pick(0)

    def do_pick(self,bondFlag):
        global dist_count

        if bondFlag:
            self.error = "Error: please select an atom, not a bond."
            print(self.error)
        else:
            if self.mode == 'pairs':
                if self.status==0:
                    name = sele_prefix
                    cmd.select(name,"(pk1)")
                    self.status = 1
                    self.error = None
                elif self.status==1:
                    if ((self.object_mode=='append') or (not dist_count)):
                        dist_count = dist_count + 1
                    else:
                        cmd.delete(dist_prefix+"%2d"%dist_count)
                    name = dist_prefix + "%02d"%dist_count
                    cmd.dist(name,sele_prefix,"(pk1)")
                    cmd.delete(sele_prefix)
                    cmd.unpick()
                    cmd.enable(name)
                    self.status = 0
            elif self.mode in ['neigh','polar','heavy']:
                if ((self.object_mode=='append') or (not dist_count)):
                    dist_count = dist_count + 1
                else:
                    cmd.delete(dist_prefix+"%2d"%dist_count)
                name = dist_prefix + "%02d"%dist_count
                cnt = 0
                if self.mode == 'neigh':
                    cnt = cmd.select(sele_prefix,"(v. and (pk1 a; %f) and (not (nbr. pk1)) and (not (nbr. (nbr. pk1))) and (not (nbr. (nbr. (nbr. pk1)))))"%self.__class__.cutoff)
                elif self.mode == 'polar':
                    cnt = cmd.select(sele_prefix,"(v. and (pk1 a; %f) and (e. n,o) and (not (nbr. pk1)) and (not (nbr. (nbr. pk1))) and (not (nbr. (nbr. (nbr. pk1)))))"%self.__class__.cutoff)
                elif self.mode == 'heavy':
                    cnt = cmd.select(sele_prefix,"(v. and (pk1 a; %f) and (not h.) and (not (nbr. pk1)) and (not (nbr. (nbr. pk1))) and (not (nbr. (nbr. (nbr. pk1)))))"%self.__class__.cutoff)
                cmd.delete(name)
                if cnt:
                    cmd.dist(name,"(pk1)",sele_prefix)
                else:
                    print(" Wizard: No neighbors found.")
                cmd.delete(sele_prefix)
                cmd.unpick()
                cmd.enable(name)
        cmd.refresh_wizard()
