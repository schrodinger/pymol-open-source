
from pymol.wizard import Wizard
from pymol import cmd
import pymol

import traceback

sele_prefix = "_mw"
sele_prefix_len = len(sele_prefix)

indi_sele = "_indicate_mw"

obj_prefix = "measure"

class Measurement(Wizard):

    cutoff = 3.5
    
    def __init__(self,_self=cmd):
        Wizard.__init__(self,_self)

        self.cmd.unpick()
        
        self.status = 0 # 0 no atoms selections, 1 atom selected, 2 atoms selected, 3 atoms selected
        self.error = None
        self.object_name = None

        # mode selection subsystem

        self.mode = self.session.get('default_mode','pairs')
        self.modes = [
            'pairs',
            'angle',
            'dihed', 
            'polar',
            'heavy',
            'neigh',
            ]
        
        self.mode_name = {
            'polar':'Polar Neighbors',
            'heavy':'Heavy Neighbors',
            'neigh':'Neighbors',
            'pairs':'Distances',
            'angle':'Angles',
            'dihed':'Dihedrals',
            }

        smm = []
        smm.append([ 2, 'Measurement Mode', '' ])
        for a in self.modes:
            smm.append([ 1, self.mode_name[a], 'cmd.get_wizard().set_mode("'+a+'")'])
        self.menu['mode']=smm

        # overwrite mode selection subsystem
        
        self.object_mode = self.session.get('default_object_mode','append')
        self.object_modes = [
            'merge', 
            'overwr',
            'append',
            ]
        self.object_mode_name = {
            'merge':'Merge With Previous',
            'overwr':'Replace Previous',
            'append':'Create New Object',         
            }

        smm = []
        smm.append([ 2, 'New Measurements?', '' ])
        for a in self.object_modes:
            smm.append([ 1, self.object_mode_name[a], 'cmd.get_wizard().set_object_mode("'+a+'")'])
        self.menu['object_mode']=smm
        self.selection_mode = self.cmd.get_setting_legacy("mouse_selection_mode")
        self.cmd.set("mouse_selection_mode",0) # set selection mode to atomic
        self.cmd.deselect() # disable the active selection (if any)

    def _validate_instance(self):
        Wizard._validate_instance(self)
        if not hasattr(self,'meas_count'):
            self.meas_count = self.session.get('meas_count',0)            
            
    def get_name(self,untaken=1,increment=1):
        self._validate_instance()
        if increment or self.meas_count<1:
            self.meas_count = self.meas_count + 1
        obj_name = obj_prefix+"%02d"%self.meas_count
        if untaken:
            name_dict = {}
            for tmp_name in cmd.get_names("all"):
                name_dict[tmp_name] = None
            while obj_name in name_dict:
                self.meas_count = self.meas_count + 1
                obj_name = obj_prefix+"%02d"%self.meas_count
        return obj_name
    
# generic set routines

    def set_mode(self,mode):
        if mode in self.modes:
            self.mode = mode
        self.status = 0
        self.clear_input()
        self.cmd.refresh_wizard()

    def set_object_mode(self,mode):
        if mode in self.object_modes:
            self.object_mode = mode
        self.status = 0
        self.cmd.refresh_wizard()

        
    def get_panel(self):
        if int(self.cmd.get("mouse_selection_mode")!=0):
            self.cmd.set("mouse_selection_mode",0)
        return [
            [ 1, 'Measurement',''],
            [ 3, self.mode_name[self.mode],'mode'],
            [ 3, self.object_mode_name[self.object_mode],'object_mode'],
            [ 2, 'Delete Last Object' , 'cmd.get_wizard().delete_last()'],
            [ 2, 'Delete All Measurements' , 'cmd.get_wizard().delete_all()'],
            [ 2, 'Done','cmd.set_wizard()'],
            ]

    def cleanup(self):
        self.session['default_mode'] = self.mode
        self.session['default_object_mode'] = self.object_mode
        self.clear_input()
        self.cmd.set("mouse_selection_mode",self.selection_mode) # restore selection mode
        
    def clear_input(self):
        self.cmd.delete(sele_prefix+"*") 
        self.cmd.delete(indi_sele)
        self.status = 0
        
    def get_prompt(self):
        self.prompt = None
        if self.mode in ['pairs', 'angle', 'dihed' ]:
            if self.status==0:
                self.prompt = [ 'Please click on the first atom...']
            elif self.status==1:
                self.prompt = [ 'Please click on the second atom...' ]
            elif self.status==2:
                self.prompt = [ 'Please click on the third atom...' ]
            elif self.status==3:
                self.prompt = [ 'Please click on the fourth atom...' ]
        elif self.mode in [ 'polar', 'neigh', 'heavy' ]:
            self.prompt = [ 'Please click an atom...']
        if self.error!=None:
            self.prompt.append(self.error)
        return self.prompt
    
    def delete_last(self):
        self._validate_instance()
        if self.status==0:
            if self.meas_count>0:
                name = self.get_name(0,0)
                self.cmd.delete(name)
                self.meas_count = self.meas_count - 1
        self.status=0
        self.error = None
        self.clear_input()
        self.cmd.refresh_wizard()

    def delete_all(self):
        
        self.meas_count = 0
        self.cmd.delete(obj_prefix+"*")
        self.status=0
        self.error = None
        self.clear_input()
        self.cmd.refresh_wizard()

    def do_select(self,name): # map selects into picks
        self.cmd.unpick()
        try:
            self.cmd.edit(name + " and not " + sele_prefix + "*") # note, using new object name wildcards
            self.cmd.delete(name)
            self.do_pick(0)
        except pymol.CmdException:
            if self.status:
                sele_name = sele_prefix + str(self.status-1)         
                self.cmd.select(indi_sele, sele_name)
                self.cmd.enable(indi_sele)

    def do_pick(self,bondFlag):
        if bondFlag:
            self.error = "Error: please select an atom, not a bond."
            print self.error
        else:
            reset = 1
            sele_name = sele_prefix + str(self.status)
            if self.mode == 'pairs':
                if self.status==0:
                    self.cmd.select(sele_name,"(pk1)")
                    self.cmd.select(indi_sele, sele_name)
                    self.cmd.enable(indi_sele)
                    self.status = 1
                    self.error = None
                elif self.status==1:
                    obj_name  = self.get_name((self.object_mode=='append'),
                                              (self.object_mode=='append'))
                    if self.object_mode=='merge':
                        reset = 0
                    self.cmd.dist(obj_name,sele_prefix+"0","(pk1)",reset=reset)
                    self.cmd.enable(obj_name)
                    self.clear_input()
                    self.status = 0
                self.cmd.unpick()
            elif self.mode == 'angle':
                if self.status<2:
                    self.cmd.select(sele_name,"(pk1)")
                    self.cmd.unpick()
                    self.cmd.select(indi_sele, sele_name)
                    self.cmd.enable(indi_sele)
                    self.status = self.status + 1
                    self.error = None
                else:
                    obj_name = self.get_name((self.object_mode=='append'),
                                             (self.object_mode=='append'))
                    if self.object_mode=='merge':
                        reset = 0
                    self.cmd.angle(obj_name, sele_prefix+"0", sele_prefix+"1",
                                 "(pk1)", reset=reset)
                    self.cmd.enable(obj_name)
                    self.clear_input()
                    self.status = 0
                self.cmd.unpick()
            elif self.mode == 'dihed':
                if self.status<3:
                    self.cmd.select(sele_name,"(pk1)")
                    self.cmd.unpick()
                    self.cmd.select(indi_sele, sele_name)
                    self.cmd.enable(indi_sele)
                    self.status = self.status + 1
                    self.error = None
                else:
                    obj_name = self.get_name((self.object_mode=='append'),
                                             (self.object_mode=='append'))
                    if self.object_mode=='merge':
                        reset = 0
                    self.cmd.dihedral(obj_name, sele_prefix+"0", sele_prefix+"1",
                                     sele_prefix+"2", "(pk1)", reset=reset)
                    self.cmd.enable(obj_name)
                    self.clear_input()
                    self.status = 0
                self.cmd.unpick()
            elif self.mode in ['neigh','polar','heavy']:
                reset = 1
                obj_name = self.get_name((self.object_mode=='append'),
                                             (self.object_mode=='append'))
                if self.object_mode=='merge':
                    reset = 0
                cnt = 0
                if self.mode == 'neigh':
                    cnt = self.cmd.select(sele_prefix,
"(v. and (pk1 a; %f) and (not (nbr. pk1)) and (not (nbr. (nbr. pk1))) and (not (nbr. (nbr. (nbr. pk1)))))"
                                          %self.__class__.cutoff)
                elif self.mode == 'polar':
                    cnt = self.cmd.select(sele_prefix,
"(v. and (pk1 a; %f) and (e. n,o) and (not (nbr. pk1)) and (not (nbr. (nbr. pk1))) and (not (nbr. (nbr. (nbr. pk1)))))"
                    %self.__class__.cutoff)            
                elif self.mode == 'heavy':
                    cnt = self.cmd.select(sele_prefix,
"(v. and (pk1 a; %f) and (not h.) and (not (nbr. pk1)) and (not (nbr. (nbr. pk1))) and (not (nbr. (nbr. (nbr. pk1)))))"
                    %self.__class__.cutoff)            
                if cnt:
                    self.cmd.dist(obj_name,"(pk1)",sele_prefix,reset=reset)
                else:
                    print " Wizard: No neighbors found."
                self.clear_input()
                self.cmd.unpick()
                self.cmd.enable(obj_name)
        self.cmd.refresh_wizard()
