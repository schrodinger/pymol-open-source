from pymol.wizard import Wizard
from pymol import cmd
import pymol

import traceback

sele_prefix = "_mw"
sele_prefix_len = len(sele_prefix)

indi_sele = "_indicate_mw"
obj_prefix = "measure"

class Measurement(Wizard):
    modes = [
        'pairs',
        'rings',
        'angle',
        'dihed',
        'polar',
        'heavy',
        'neigh',
        'hbond',
        ]
    mode_name = {
        'polar':'Polar Neighbors',
        'heavy':'Heavy Neighbors',
        'neigh':'Neighbors',
        'pairs':'Distances',
        'angle':'Angles',
        'dihed':'Dihedrals',
        'hbond':'Polar Contacts',
        'rings':'Distances to Rings',
        }
    neighbor_modes = [
        'same',
        'other',
        'enabled',
        'all',
        'in_object',
        'in_selection'
        ]
    object_modes = [
        'merge',
        'overwr',
        'append',
        ]
    object_mode_name = {
        'merge':'Merge With Previous',
        'overwr':'Replace Previous',
        'append':'Create New Object',
        }

    def __init__(self,_self=cmd):
        Wizard.__init__(self,_self)

        self.pick_state = [0, 0, 0, 0]

        self.cmd.unpick()

        self.cutoff = self.cmd.get_setting_float("neighbor_cutoff")
        self.heavy_neighbor_cutoff = self.cmd.get_setting_float("heavy_neighbor_cutoff")
        self.polar_neighbor_cutoff = self.cmd.get_setting_float("polar_neighbor_cutoff")
        self.hbond_cutoff = self.cmd.get_setting_float("h_bond_cutoff_center")

        self.status = 0 # 0 no atoms selections, 1 atom selected, 2 atoms selected, 3 atoms selected
        self.error = None
        self.object_name = None

        # mode selection subsystem
        self.mode = self.session.get('default_mode','pairs')

        self.neighbor_target = ""

        # TODO:
        # make this a function, and call it when we call refresh wizard
        # to update the object/selection list
        smm = []
        smm.append([ 2, 'Measurement Mode', '' ])
        for a in self.modes:
            if a in ("neigh", "polar", "heavy"):
                smm.append([ 1, self.mode_name[a], self.neighbor_submenu(a)])
            else:
                smm.append([ 1, self.mode_name[a], 'cmd.get_wizard().set_mode("'+a+'")'])
        self.menu['mode']=smm

        # overwrite mode selection subsystem
        self.object_mode = self.session.get('default_object_mode','append')

        smm = []
        smm.append([ 2, 'New Measurements?', '' ])
        for a in self.object_modes:
            smm.append([ 1, self.object_mode_name[a], 'cmd.get_wizard().set_object_mode("'+a+'")'])

        self.menu['object_mode']=smm
        # initially select atoms, but now users can change this
        self.selection_mode = self.cmd.get_setting_int("mouse_selection_mode")
        self.cmd.set("mouse_selection_mode",0) # set selection mode to atomic
        self.cmd.deselect() # disable the active selection (if any)
        self.mouse_mode = 0

    def get_event_mask(self):
        """
        Sets what this Wizard listens for.  event_mask_dirty is too coarse an event,
        so we should consider something more fine grain for mouse mode updates
        """
        return Wizard.event_mask_pick + Wizard.event_mask_select + Wizard.event_mask_dirty

    def neighbor_objects(self,a):
        """
        for neighbor selecting, populate the menu with these names
        """
        list = self.cmd.get_names("public_objects",1)[0:25] # keep this practical
        list = [x for x in list if self.cmd.get_type(x)=="object:molecule"]
        result = [[ 2, 'Object: ', '']]
        for b in list:
            result.append( [ 1, b,  'cmd.get_wizard().set_neighbor_target("%s","%s")' % (a,b)])
        return result

    def neighbor_selections(self,a):
        """
        get list of public selections for populating the menu
        """
        list = self.cmd.get_names("public_selections",1)[0:25] # keep this practical
        list = [x for x in list if self.cmd.get_type(x)=="selection"]
        result = [[ 2, 'Selections: ', '']]
        for b in list:
            result.append( [ 1, b,  'cmd.get_wizard().set_neighbor_target("%s","%s")' % (a,b)])
        return result

    def neighbor_submenu(self,a,_self=cmd):
        return [ [2, self.mode_name[a]+": ", ''],
                 [1, "in all objects",  'cmd.get_wizard().set_neighbor_target("'+a+'","all")'],
                 [1, "in object", self.neighbor_objects(a) ],               # while this is a submenu this is also a command, so [1, ...]
                 [1, "in selection", self.neighbor_selections(a) ],         # not [2, ...]
                 [1, "in other objects",  'cmd.get_wizard().set_neighbor_target("'+a+'","other")'],
                 [1, "in same object",  'cmd.get_wizard().set_neighbor_target("'+a+'", "same")'],
            ]

    def _validate_instance(self):
        Wizard._validate_instance(self)
        if not hasattr(self,'meas_count'):
            self.meas_count = self.session.get('meas_count',0)

    def get_name(self,untaken=1,increment=1):
        """
        get a name for the next measurement object
        """
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

    def set_neighbor_target(self,mode,target):
        """
        sets the neighbor target in the menu
        """
        self.set_mode(mode)
        self.neighbor_target=target
        self.status = 0
        self.clear_input()
        self.cmd.refresh_wizard()

    def set_mode(self,mode):
        """
        sets what we're measuring, distance, angle, dihedral, etc.
        """
        if mode in self.modes:
            self.mode = mode
        # if setting mode, we're restarting the selection process
        self.status = 0
        self.clear_input()

        if self.mode=='hbond':
            self.cmd.set("mouse_selection_mode", 5)

        self.cmd.refresh_wizard()

    def set_object_mode(self,mode):
        if mode in self.object_modes:
            self.object_mode = mode
        self.status = 0
        self.cmd.refresh_wizard()


    def get_panel(self):
        return [
            [ 1, 'Measurement',''],
            [ 3, self.mode_name[self.mode],'mode'],
            [ 3, self.object_mode_name[self.object_mode],'object_mode'],
            [ 2, 'Delete Last Object' , 'cmd.get_wizard().delete_last()'],
            [ 2, 'Delete All Measurements' , 'cmd.get_wizard().delete_all()'],
            [ 2, 'Done','cmd.set_wizard()'],
            ]

    def cleanup(self):
        """
        restore user session how we found it
        """
        self.session['default_mode'] = self.mode
        self.session['default_object_mode'] = self.object_mode
        self.clear_input()
        self.cmd.set("mouse_selection_mode",self.selection_mode) # restore selection mode

    def clear_input(self):
        """
        delete our user selections for this wizard
        """
        self.cmd.delete(sele_prefix+"*")
        self.cmd.delete(indi_sele)
        self.cmd.delete("pk1")
        self.status = 0

    def get_selection_name(self):
        """
        Return a textual description of the mouse mode
        """
        if self.cmd.get("mouse_selection_mode",   quiet=1)=="0":
            return ("atom","")
        elif self.cmd.get("mouse_selection_mode", quiet=1)=="1":
            return ("residue"," br. ")
        elif self.cmd.get("mouse_selection_mode", quiet=1)=="2":
            return ("chain", " bc. ")
        elif self.cmd.get("mouse_selection_mode", quiet=1)=="3":
            return ("segment", " bs. ")
        elif self.cmd.get("mouse_selection_mode", quiet=1)=="4":
            return ("object", " bo. ")
        elif self.cmd.get("mouse_selection_mode", quiet=1)=="5":
            return ("molecule", " bm. ")
        elif self.cmd.get("mouse_selection_mode", quiet=1)=="6":
            return ("C-alpha", " bca. ")

    def get_prompt(self):
        (what, code) = self.get_selection_name()
        self.prompt = None
        if self.mode in ['pairs', 'angle', 'dihed', 'hbond', 'rings' ]:
            if self.status==0:
                self.prompt = [ 'Please click on the first %s...' % what]
            elif self.status==1:
                self.prompt = [ 'Please click on the second %s...' % what ]
            elif self.status==2:
                self.prompt = [ 'Please click on the third %s...' % what]
            elif self.status==3:
                self.prompt = [ 'Please click on the fourth %s...' % what]
        elif self.mode in [ 'polar', 'neigh', 'heavy', 'surf' ]:
            (what, code) = self.get_selection_name()
            letterN=""
            if what[0] in ('a', 'e', 'i', 'o', 'u'):
                letterN = "n"
            else:
                letterN = ""
            self.prompt = [ 'Please click a%s %s...' % (letterN, what)]
        if self.error is not None:
            self.prompt.append(self.error)
        return self.prompt

    def delete_last(self):
        """
        Corresponds to the "Delete Last Object" menu button
        """
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
        """
        Corresponds to the "Delete All Measurements" menu button
        """
        self.meas_count = 0
        self.cmd.delete(obj_prefix+"*")
        self.status=0
        self.error = None
        self.clear_input()
        self.cmd.refresh_wizard()

    def do_pick_state(self, state):
        self.pick_state[self.status] = state

    def do_select(self,name): # map selects into picks
        self.cmd.unpick()
        try:
            self.cmd.select("pk1", name + " and not " + sele_prefix + "*") # note, using new object name wildcards
            self.cmd.delete(name)
            self.do_pick(0)
        except pymol.CmdException:
            if self.status:
                sele_name = sele_prefix + str(self.status-1)
                self.cmd.select(indi_sele, sele_name)
                self.cmd.enable(indi_sele)

    def do_pick(self,bondFlag):
        # update pk1 based on current mouse mode
        (what,code) = self.get_selection_name()
        self.cmd.select("pk1", code + "(pk1)")
        if bondFlag:
            self.error = "Error: please select an atom, not a bond."
            print(self.error)
        else:
            reset = 1
            sele_name = sele_prefix + str(self.status)
            state1, state2, state3, state4 = self.pick_state

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

                    if state1 == state2 and self.cmd.get_selection_state('?pk1') != 0:
                        state1 = state2 = -3

                    self.cmd.dist(obj_name,"(v. and " + sele_prefix+"0)","(v. and (pk1))",reset=reset,
                            state1=state1, state2=state2)

                    self.cmd.enable(obj_name)
                    self.clear_input()
                    self.status = 0
                self.cmd.unpick()
            elif self.mode == 'rings':
                if self.status==0:
                    self.cmd.select(sele_name,"(?pk1 | byring ?pk1)")
                    self.cmd.select(indi_sele, sele_name)
                    self.cmd.enable(indi_sele)
                    self.status = 1
                    self.error = None
                elif self.status==1:
                    obj_name  = self.get_name((self.object_mode=='append'),
                                              (self.object_mode=='append'))
                    if self.object_mode=='merge':
                        reset = 0

                    if state1 == state2 and self.cmd.get_selection_state('?pk1') != 0:
                        state1 = state2 = -3

                    self.cmd.distance(obj_name,
                            "?" + sele_prefix + "0",
                            "(?pk1 | byring ?pk1)",
                            mode=4, reset=reset,
                            state1=state1, state2=state2)
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

                    if state1 == state2 == state3 and self.cmd.get_selection_state('?pk1') != 0:
                        state1 = state2 = state3 = -3

                    self.cmd.angle(obj_name, "(v. and " + sele_prefix+"0)", "(v. and " + sele_prefix+"1)",
                                   "(v. and (pk1))", reset=reset,
                                   state1=state1, state2=state2, state3=state3)
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
                    self.cmd.dihedral(obj_name, "(v. and " + sele_prefix+"0)", "(v. and " + sele_prefix+"1)",
                                      "(v. and " + sele_prefix+"2)", "(v. and (pk1))", reset=reset)
                    self.cmd.enable(obj_name)
                    self.clear_input()
                    self.status = 0
                self.cmd.unpick()
            if self.mode == 'hbond':
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
                    self.cmd.dist(obj_name,"(v. and " + sele_prefix+"0)","(v. and (pk1))",mode=2,cutoff=self.hbond_cutoff,reset=reset)
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
                sel_mod = ""
                if self.mode in ('neigh', 'polar', 'heavy'):
                    if self.neighbor_target=="same":
                        sel_mod = "bm. pk1"
                    elif self.neighbor_target=="other":
                        sel_mod = "(not bm. pk1)"
                    elif self.neighbor_target=="enabled":
                        sel_mod = "(enabled)"
                    elif self.neighbor_target=="all":
                        sel_mod = "all"
                    else:
                        sel_mod = "(%s)" % self.neighbor_target

                cutoffType=0
                if self.mode == 'neigh':
                    cnt = self.cmd.select(sele_prefix,
"(v. and (pk1 a; %f) and (not (nbr. pk1)) and (not (nbr. (nbr. pk1))) and (not (nbr. (nbr. (nbr. pk1)))) and (%s))"
                                          %(self.cutoff, sel_mod))
                    cutoffType=self.cutoff
                elif self.mode == 'polar':
                    cnt = self.cmd.select(sele_prefix,
"(v. and (pk1 a; %f) and (e. n,o) and (not (nbr. pk1)) and (not (nbr. (nbr. pk1))) and (not (nbr. (nbr. (nbr. pk1)))) and (%s))"
                    %(self.polar_neighbor_cutoff, sel_mod))
                    cutoffType = self.polar_neighbor_cutoff
                elif self.mode == 'heavy':
                    cutoffType = self.heavy_neighbor_cutoff
                    cnt = self.cmd.select(sele_prefix,
"(v. and (pk1 a; %f) and (not h.) and (not (nbr. pk1)) and (not (nbr. (nbr. pk1))) and (not (nbr. (nbr. (nbr. pk1)))) and (%s))"
                    %(self.heavy_neighbor_cutoff, sel_mod))
                if cnt:
                    self.cmd.dist(obj_name,"(pk1)",sele_prefix,cutoff=cutoffType,reset=reset)
                else:
                    print(" Wizard: No neighbors found.")
                self.clear_input()
                self.cmd.unpick()
                self.cmd.enable(obj_name)
        self.cmd.refresh_wizard()

    def do_dirty(self):
        if self.mouse_mode != self.cmd.get("mouse_selection_mode"):
            self.mouse_mode = self.cmd.get("mouse_selection_mode")
            self.cmd.refresh_wizard()
