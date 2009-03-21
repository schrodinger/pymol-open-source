
import sys
import os
from glob import glob
import traceback
import string
import types

from Tkinter import *
import tkMessageBox
import Pmw

from pymol import editor
from pymol import computing
from pymol.wizard import Wizard
from pymol.parsing import QuietException

from chempy import cpv

import pymol

WID = 5
MAX_COL = 12
imgDict = {}
active_sele = "_builder_active" # object we're working on...
newest_sele = "_builder_added" # last atom added?

display_sele = "_build_display"


#############################################################
# action-directing Wizards

class ActionWizard(Wizard):

    def __init__(self,_self):
        Wizard.__init__(self,_self)
        self.actionHash = str(self.__class__)
        
    def setActionHash(self, action_hash):
        self.actionHash = action_hash

    def activateOrDismiss(self):
        activate_flag = 1
        cur_wiz = self.cmd.get_wizard()
        if cur_wiz != None:
            if cur_wiz.__class__ == self.__class__:
                if cur_wiz.actionHash == self.actionHash:
                    activate_flag = 0
        if activate_flag:
            self.cmd.set_wizard(self,replace=1)
            self.cmd.refresh_wizard()
        else:
            self.actionWizardDone()
        return activate_flag

    def actionWizardDone(self):
        self.cmd.delete(active_sele)
        self.cmd.unpick()
        self.cmd.set_wizard()
        self.cmd.refresh_wizard()

    def activeSeleValid(self):
        if active_sele in self.cmd.get_names("selections"):
            if self.cmd.select(active_sele, "byobj "+active_sele)<1:
                self.cmd.delete(active_sele)
        elif "pk1" in self.cmd.get_names("selections"):
            self.cmd.select(active_sele,"byobj pk1")
        else:
            obj_list = self.cmd.get_names("objects",enabled_only=1)
            if len(obj_list)==1:
                if self.cmd.select(active_sele, obj_list[0])<1:
                    self.cmd.delete(active_sele)
        return active_sele in self.cmd.get_names("selections")

class CleanWizard(ActionWizard):

    def __init__(self,_self):
        self.clean_obj = None
        ActionWizard.__init__(self,_self)

    def run_job(self):
        if active_sele in self.cmd.get_names("selections"):
            obj_list = self.cmd.get_object_list(active_sele)
            if len(obj_list)==1:
                self.cmd.unpick()
                self.cmd.set_wizard()
                self.cmd.refresh_wizard()
                self.cmd.do("_ cmd.clean('%s',message='''Cleaning %s...''',async=1)"%(active_sele,obj_list[0]))
#                computing.CleanJob(self.cmd,active_sele,
#                                   message="Cleaning %s..."%obj_list[0])

    def do_pick(self, bondFlag):
        if active_sele in self.cmd.get_names("selections"):
            obj_list = self.cmd.get_object_list(active_sele)
            if len(obj_list)!=1:
                self.cmd.delete(active_sele)
        else:
            self.cmd.select(active_sele, "byobj pk1")
        self.cmd.unpick()
        self.cmd.deselect()
        obj_list = self.cmd.get_object_list(active_sele)
        if isinstance(obj_list,types.ListType) and (len(obj_list)==1):
            self.run_job()
        else:
            print "Error: can only clean one object at a time"

    def toggle(self):
        if self.activateOrDismiss():
            if self.activeSeleValid():
                self.run_job()

    def get_prompt(self):
        return ["Pick object to clean..."]

    
    def get_panel(self):
        return [
            [ 1, 'Clean', ''],
            [ 2, 'Done','cmd.set_wizard()'],
            ]

class SculptWizard(ActionWizard):

    def __init__(self,_self):
        ActionWizard.__init__(self,_self)
        self.sculpt_object = None
        
    def sculpt_activate(self):
        if active_sele in self.cmd.get_names("selections"):
            obj_list = self.cmd.get_object_list(active_sele)
            if len(obj_list)==1:
                obj_name = obj_list[0]
                self.cmd.push_undo(obj_name)
                self.cmd.sculpt_activate(obj_name)
                self.cmd.set("sculpting",1)
                self.sculpt_object = obj_name
                self.cmd.sculpt_activate(obj_name)
                if int(self.cmd.get("sculpt_vdw_vis_mode")):
                    self.cmd.show("cgo",obj_name)
                self.cmd.unpick()
                self.cmd.refresh_wizard()
            else:
                print "Error: cannot sculpt more than one object at a time"
        
    def sculpt_deactivate(self):
        if self.sculpt_object != None:
            self.cmd.set("sculpt_vdw_vis_mode","0",self.sculpt_object)
            self.cmd.sculpt_iterate(self.sculpt_object,self.cmd.get_state(),0)
            self.cmd.unset("sculpt_vdw_vis_mode",self.sculpt_object)
            self.cmd.sculpt_deactivate(self.sculpt_object)
            self.sculpt_object = None
            self.cmd.refresh_wizard()            
        
    def do_pick(self, bondFlag):
        if self.sculpt_object == None:
            self.cmd.select(active_sele, "byobj pk1")
            self.sculpt_activate()
        else:
            return 0 # already sculpting, so handle like a normal edit
        
    def toggle(self):
        if self.activateOrDismiss():
            if self.activeSeleValid():
                self.sculpt_activate()

    def get_prompt(self):
        if self.sculpt_object == None:
            return ["Pick object to sculpt..."]
        else:
            return ["Sculpting %s..."%self.sculpt_object]

    def finish_sculpting(self):
        if self.sculpt_object:
            self.sculpt_deactivate()
        self.cmd.set("sculpting",0)
        self.cmd.delete(active_sele)
        self.cmd.set_wizard()
        self.cmd.refresh_wizard()
        
    def scramble(self,mode):
        if self.cmd.count_atoms(self.sculpt_object):
            sc_tmp = "_scramble_tmp"
            if mode == 0:
                self.cmd.select(sc_tmp,self.sculpt_object+
                            " and not (fixed or restrained)")
            if mode == 1:
                self.cmd.select(sc_tmp,self.sculpt_object+
                            " and not (fixed)")
            extent = self.cmd.get_extent(sc_tmp)
            center = self.cmd.get_position(sc_tmp)
            radius = 1.25*cpv.length(cpv.sub(extent[0],extent[1]))
            self.cmd.alter_state(self.cmd.get_state(), sc_tmp,
                                 "(x,y,z)=rsp(pos,rds)",
                space= { 'rsp' :  cpv.random_displacement,
                         'pos' : center,
                         'rds' : radius })
            self.cmd.delete(sc_tmp)

    def get_panel(self):
        return [
            [ 1, 'Sculpt', ''],
            [ 2, 'Undo', 'cmd.undo()'],

            [ 2, 'Switch Object', 'cmd.get_wizard().sculpt_deactivate()'],
            [ 2, 'Scramble Unrestrained Coords.', 'cmd.get_wizard().scramble(0)'],
            [ 2, 'Scramble Unfixed Coords.', 'cmd.get_wizard().scramble(1)'],
            [ 2, 'Done','cmd.get_wizard().finish_sculpting()'],
            ]

    def cleanup(self):
        self.sculpt_deactivate()
        ActionWizard.cleanup(self)

class RepeatableActionWizard(ActionWizard):

    def __init__(self,_self):
        ActionWizard.__init__(self,_self)
        self.repeating = 0

    def repeat(self):
        self.repeating = 1
        self.cmd.refresh_wizard()

    def getRepeating(self):
        return self.repeating

    def activateRepeatOrDismiss(self):
        activate_flag = 1
        cur_wiz = self.cmd.get_wizard()
        if cur_wiz != None:
            if cur_wiz.__class__ == self.__class__:
                if cur_wiz.actionHash == self.actionHash:
                    if cur_wiz.getRepeating():
                        activate_flag = 0
                    else:
                        self.repeat()
                elif cur_wiz.getRepeating():
                    self.repeat()
        if activate_flag:
            self.cmd.set_wizard(self,replace=1)
            self.repeat() # always repeating for now...
            self.cmd.refresh_wizard()
        else:
            self.actionWizardDone()
        return activate_flag

    def cleanup(self):
        self.cmd.unpick()
        ActionWizard.cleanup(self)

class ReplaceWizard(RepeatableActionWizard):

    def do_pick(self, bondFlag):
        self.cmd.select(active_sele, "bymol pk1") 
        self.cmd.replace(self.symbol, self.geometry, self.valence)
        if not self.getRepeating():
            self.actionWizardDone()

    def toggle(self, symbol, geometry, valence, text):
        self.symbol = symbol
        self.geometry = geometry
        self.valence = valence
        self.text = text
        self.setActionHash( (symbol,geometry,valence,text) )
        self.activateRepeatOrDismiss()
                                    
    def get_prompt(self):
        if self.getRepeating():
            return ["Pick atoms to replace with %s..."%self.text]
        else:
            return ["Pick atom to replace with %s..."%self.text]            

    def get_panel(self):
        if self.getRepeating():
            return [
                [ 1, 'Replacing Multiple Atoms',''],
                [ 2, 'Done','cmd.set_wizard()'],
                ]
        else:
            return [
                [ 1, 'Replacing an Atom',''],
                [ 2, 'Replace Multiple Atoms','cmd.get_wizard().repeat()'],
                [ 2, 'Done','cmd.set_wizard()'],
                ]

class AttachWizard(RepeatableActionWizard):

    def __init__(self,_self):
        RepeatableActionWizard.__init__(self,_self)
        self.mode = 0

    def do_pick(self, bondFlag):
        if self.mode == 0:
            self.cmd.select(active_sele, "bymol pk1") 
            editor.attach_fragment("pk1", self.fragment, self.position, self.geometry, _self=self.cmd)
        elif self.mode == 1:            
            self.cmd.select(active_sele, "bymol pk1") 
            editor.combine_fragment("pk1", self.fragment, self.position,
                                    self.geometry, _self=self.cmd)
            self.mode = 0
            self.cmd.refresh_wizard()
        self.cmd.unpick()
        if not self.getRepeating():
            self.actionWizardDone()
            
    def toggle(self, fragment, position, geometry, text):
        self.fragment = fragment
        self.position = position
        self.geometry = geometry
        self.text = text
        self.setActionHash( (fragment, position, geometry, text) )
        self.activateRepeatOrDismiss()

    def create_new(self):
        names = self.cmd.get_names("objects")
        num = 1
        while 1:
            name = "obj%02d"%num
            if name not in names:
                break
            num = num + 1
        self.cmd.fragment(self.fragment, name)
        if not self.getRepeating():
            self.actionWizardDone()
    
    def get_prompt(self):
        if self.mode == 0:
            if self.getRepeating():
                return ["Pick locations to attach %s..."%self.text]
            else:
                return ["Pick location to attach %s..."%self.text]            
        else:
            return ["Pick object to combine %s into..."%self.text]

    def combine(self):
        self.mode = 1
        self.cmd.refresh_wizard()

    def get_panel(self):
        if self.getRepeating():
            return [
                [ 1, 'Attaching Multiple Fragmnts',''],
                [ 2, 'Create As New Object','cmd.get_wizard().create_new()'],                
                [ 2, 'Combine w/ Existing Object','cmd.get_wizard().combine()'],                
                [ 2, 'Done','cmd.set_wizard()'],
                ]
        else:
            return [
                [ 1, 'Attaching One Fragment',''],
                [ 2, 'Create As New Object','cmd.get_wizard().create_new()'],
                [ 2, 'Combine w/ Existing Object','cmd.get_wizard().combine()'],                
                [ 2, 'Attach Multiple Fragments','cmd.get_wizard().repeat()'],
                [ 2, 'Done','cmd.set_wizard()'],
                ]

class AminoAcidWizard(RepeatableActionWizard):

    def __init__(self,_self):
        RepeatableActionWizard.__init__(self,_self)
        self.mode = 0

    def do_pick(self, bondFlag):
        if self.mode == 0:
            self.cmd.select(active_sele, "bymol pk1") 
            try:
                editor.attach_amino_acid("pk1", self.aminoAcid, _self=self.cmd)        
            except QuietException:
                pass
        elif self.mode == 1:
            self.cmd.select(active_sele, "bymol pk1") 
            editor.combine_fragment("pk1", self.aminoAcid, 0, 1, _self=self.cmd)
            self.mode = 0
            self.cmd.refresh_wizard()

        self.cmd.unpick()
        if not self.getRepeating():
            self.actionWizardDone()
            
    def toggle(self, amino_acid):
        self.aminoAcid = amino_acid
        self.setActionHash( (amino_acid,) )
        self.activateRepeatOrDismiss()

    def create_new(self):
        names = self.cmd.get_names("objects")
        num = 1
        while 1:
            name = "obj%02d"%num
            if name not in names:
                break
            num = num + 1
        editor.attach_amino_acid("pk1", self.aminoAcid, object=name, _self=self.cmd)        
        if not self.getRepeating():
            self.actionWizardDone()

    def combine(self):
        self.mode = 1
        self.cmd.refresh_wizard()
    
    def get_prompt(self):
        if self.mode == 0:
            if self.getRepeating():
                return ["Pick locations to attach %s..."%self.aminoAcid]
            else:
                return ["Pick location to attach %s..."%self.aminoAcid]            
        else:
            return ["Pick object to combine %s into..."%self.aminoAcid]

    def get_panel(self):
        if self.getRepeating():
            return [
                [ 1, 'Attaching Multiple Residues',''],
                [ 2, 'Create As New Object','cmd.get_wizard().create_new()'],                
                [ 2, 'Combine w/ Existing Object','cmd.get_wizard().combine()'],                
                [ 2, 'Done','cmd.set_wizard()'],
                ]
        else:
            return [
                [ 1, 'Attaching Amino Acid',''],
                [ 2, 'Create As New Object','cmd.get_wizard().create_new()'],
                [ 2, 'Combine w/ Existing Object','cmd.get_wizard().combine()'],                
                [ 2, 'Attach Multiple...','cmd.get_wizard().repeat()'],
                [ 2, 'Done','cmd.set_wizard()'],
                ]
            
class ValenceWizard(RepeatableActionWizard):

    def cleanup(self):
        self.cmd.button('single_left','none','PkAt')
        self.cmd.button('double_left','none','MovA')
        
    def do_pick(self, bondFlag):
        self.cmd.select(active_sele, "bymol pk1") 
        if bondFlag:
            if self.order>=0:
                self.cmd.valence(self.order, "pk1", "pk2")
                self.cmd.h_fill()
            else:
                self.cmd.cycle_valence()
            self.cmd.unpick()
        else:
            self.cmd.button('double_left','none','PkBd')
            self.cmd.button('single_left','none','PkBd')
            self.cmd.unpick()
        if not self.getRepeating():
            self.actionWizardDone()

            
    def toggle(self, order, text):
        self.order = order
        self.text = text
        self.setActionHash( (order,text) )
        self.activateRepeatOrDismiss()
        if self.cmd.get_wizard() == self:
            # get us into bond picking mode...
            self.cmd.button('double_left','none','PkBd')
            self.cmd.button('single_left','none','PkBd') 
                                    
    def get_prompt(self):
        if self.getRepeating():
            return ["Pick bonds to set as %s..."%self.text]
        else:
            return ["Pick bond to set as %s..."%self.text]            

    def get_panel(self):
        if self.getRepeating():
            return [
                [ 1, 'Setting Multiple Valences',''],
                [ 2, 'Done','cmd.set_wizard()'],
                ]
        else:
            return [
                [ 1, 'Set a Bond Valence',''],
                [ 2, 'Set Multiple Valences','cmd.get_wizard().repeat()'],
                [ 2, 'Done','cmd.set_wizard()'],
                ]

class ChargeWizard(RepeatableActionWizard):

    def do_pick(self, bondFlag):
        self.cmd.select(active_sele, "bymol pk1") 
        self.cmd.alter("pk1","formal_charge=%s" % self.charge)
        self.cmd.h_fill()
        if abs(float(self.charge))>0.0001:
            self.cmd.label("pk1","'''"+self.text+"'''")
        else:
            self.cmd.label("pk1")
        self.cmd.unpick()
        if not self.getRepeating():
            self.actionWizardDone()
            
    def toggle(self, charge, text):
        self.charge = charge
        self.text = text
        self.setActionHash( (charge,text) )
        self.activateRepeatOrDismiss()
                                    
    def get_prompt(self):
        if self.getRepeating():
            return ["Pick atoms to set charge = %s..."%self.text]
        else:
            return ["Pick atom to set charge = %s..."%self.text]            

    def get_panel(self):
        if self.getRepeating():
            return [
                [ 1, 'Setting Multiple Charges',''],
                [ 2, 'Done','cmd.set_wizard()'],
                ]
        else:
            return [
                [ 1, 'Setting Atom Charge',''],
                [ 2, 'Modify Multiple Atoms','cmd.get_wizard().repeat()'],
                [ 2, 'Done','cmd.set_wizard()'],
                ]

class InvertWizard(RepeatableActionWizard):

    def do_pick(self, bondFlag):
        self.cmd.select(active_sele, "bymol pk1") 
        picked = collectPicked(self.cmd)
        if picked == ["pk1","pk2","pk3"]:
            self.cmd.invert()
            self.cmd.unpick()
            if not self.getRepeating():
                self.actionWizardDone()
        self.cmd.refresh_wizard()
        
    def toggle(self):
        self.activateRepeatOrDismiss()
                                    
    def get_prompt(self):
        if "pk1" in self.cmd.get_names("selections"):
            if "pk2" in self.cmd.get_names("selections"):
                return ["Pick the second stationary atom..."]
            else:
                return ["Pick the first stationary atom..."]                
        else:
            return ["Pick origin atom for inversion..."]

    def get_panel(self):
        if self.getRepeating():
            return [
                [ 1, 'Inverting Multiple',''],
                [ 2, 'Done','cmd.set_wizard()'],
                ]
        else:
            return [
                [ 1, 'Inverting Stereocenter',''],
                [ 2, 'Invert Multiple','cmd.get_wizard().repeat()'],
                [ 2, 'Done','cmd.set_wizard()'],
                ]
        
class BondWizard(RepeatableActionWizard):

    def do_pick(self, bondFlag):
        self.cmd.select(active_sele, "bymol pk1") 
        picked = collectPicked(self.cmd)
        if picked == ["pk1","pk2"]:
            if ((self.cmd.count_atoms("pk1 and hydro") == self.cmd.count_atoms("pk1")) and
                (self.cmd.count_atoms("(pk1 extend 1) and not hydro")==1) and
                (self.cmd.count_atoms("pk2 and hydro") == self.cmd.count_atoms("pk2")) and
                (self.cmd.count_atoms("(pk2 extend 1) and not hydro")==1)):
                self.cmd.select("pk1","(pk1 extend 1) and not hydro")
                self.cmd.select("pk2","(pk2 extend 1) and not hydro")
            self.cmd.bond("pk1", "pk2")
            self.cmd.h_fill()
            self.cmd.unpick()
            if not self.getRepeating():
                self.actionWizardDone()
        self.cmd.refresh_wizard()

    def toggle(self):
        self.activateRepeatOrDismiss()
                                    
    def get_prompt(self):
        if "pk1" in self.cmd.get_names("selections"):
            return ["Pick second atom for bond..."]
        else:
            return ["Pick first atom for bond..."]

    def get_panel(self):
        if self.getRepeating():
            return [
                [ 1, 'Creating Multiple Bonds',''],
                [ 2, 'Done','cmd.set_wizard()'],
                ]
        else:
            return [
                [ 1, 'Creating Bond',''],
                [ 2, 'Create Multiple Bonds','cmd.get_wizard().repeat()'],
                [ 2, 'Done','cmd.set_wizard()'],
                ]

class UnbondWizard(RepeatableActionWizard):

    def cleanup(self):
        self.cmd.button('single_left','none','PkAt')
        
    def do_pick(self, bondFlag):
        self.cmd.select(active_sele, "bymol pk1") 
        if bondFlag:
            self.cmd.unbond("pk1", "pk2")
            self.cmd.h_fill()
            self.cmd.unpick()
        else:
            self.cmd.button('single_left','none','PkBd')
            self.cmd.unpick()
        if not self.getRepeating():
            self.actionWizardDone()

            
    def toggle(self):
        self.activateRepeatOrDismiss()
        if self.cmd.get_wizard() == self:
            self.cmd.button('single_left','none','PkBd') # get us into bond picking mode...
                                    
    def get_prompt(self):
        if self.getRepeating():
            return ["Pick bonds to delete..."]
        else:
            return ["Pick bond to delete..."]

    def get_panel(self):
        if self.getRepeating():
            return [
                [ 1, 'Deleting Multiple Bonds',''],
                [ 2, 'Done','cmd.set_wizard()'],
                ]
        else:
            return [
                [ 1, 'Deleting a Bond',''],
                [ 2, 'Delete Multiple Bonds','cmd.get_wizard().repeat()'],
                [ 2, 'Done','cmd.set_wizard()'],
                ]

class HydrogenWizard(RepeatableActionWizard):

    def run_add(self):
        if self.mode == 'add':
            self.cmd.h_add(active_sele)
            self.cmd.delete(active_sele)

    def do_pick(self, bondFlag):
        self.cmd.select(active_sele, "bymol pk1") 
        if self.mode == 'fix':
            self.cmd.h_fill()
            self.cmd.unpick()
        elif self.mode == 'add':
            self.cmd.unpick()
            self.run_add()
        if not self.getRepeating():
            self.actionWizardDone()
            
    def toggle(self,mode):
        self.mode = mode
        self.setActionHash( (mode,) )

        if self.mode == 'add':
            if self.activateOrDismiss():
                if self.activeSeleValid():
                    self.run_add()
        else:
            self.activateRepeatOrDismiss()
                                    
    def get_prompt(self):
        if self.mode == 'fix':
            if self.getRepeating():           
                return ["Pick atom upon which to fix hydrogens..."]
            else:
                return ["Pick atoms upon which to fix hydrogens..."]
        else:
            return ["Pick molecule upon which to add hydrogens..."]

    def get_panel(self):
        if self.getRepeating():
            if self.mode == 'fix':
                return [
                    [ 1, 'Fixing Hydrogens',''],
                    [ 2, 'Done','cmd.set_wizard()'],
                    ]
            else:
                return [
                    [ 1, 'Adding Hydrogens',''],
                    [ 2, 'Done','cmd.set_wizard()'],
                    ]
        else:
            if self.mode == 'fix':
                return [
                    [ 1, 'Fixing Hydrogens',''],
                    [ 2, 'Fix Multiple Atoms','cmd.get_wizard().repeat()'],
                    [ 2, 'Done','cmd.set_wizard()'],
                    ]
            else:
                return [
                    [ 1, 'Adding Hydrogens',''],
                    [ 2, 'Add To Multiple...','cmd.get_wizard().repeat()'],
                    [ 2, 'Done','cmd.set_wizard()'],
                    ]

class RemoveWizard(RepeatableActionWizard):

    def do_pick(self, bondFlag):
        cnt = self.cmd.select(active_sele,
                              "((pk1 and not hydro) extend 1) and not hydro")
        self.cmd.remove_picked()
        if cnt:
            self.cmd.fix_chemistry(active_sele)
            self.cmd.h_add(active_sele)
        self.cmd.delete(active_sele)
        if not self.getRepeating():
            self.actionWizardDone()
            
    def toggle(self):
        self.activateRepeatOrDismiss()
                                    
    def get_prompt(self):
        if self.getRepeating():
            return ["Pick atoms to delete..."]
        else:
            return ["Pick atom to delete..."]

    def get_panel(self):
        if self.getRepeating():
            return [
                [ 1, 'Deleting Atoms',''],
                [ 2, 'Done','cmd.set_wizard()'],
                ]
        else:
            return [
                [ 1, 'Deleting an Atom',''],
                [ 2, 'Delete Multiple Atoms','cmd.get_wizard().repeat()'],
                [ 2, 'Done','cmd.set_wizard()'],
                ]

class AtomFlagWizard(ActionWizard):

    def update_display(self):
        if active_sele in self.cmd.get_names("selections"):
            self.cmd.select(display_sele,active_sele+
                            " and flag %d"%self.flag)
            self.cmd.enable(display_sele)
        else:
            self.cmd.delete(display_sele)
        self.cmd.refresh_wizard()
        
    def do_pick(self, bondFlag):
        if not(active_sele not in self.cmd.get_names("selections")):
            if self.cmd.count_atoms("pk1 and flag %d"%self.flag):
                self.cmd.flag(self.flag,"pk1","clear")
            else:
                self.cmd.flag(self.flag,"pk1","set")            
        self.cmd.select(active_sele, "byobj pk1")
        self.cmd.unpick()
        self.update_display()
       
    def do_select(self,selection):
        if selection == display_sele:
            self.cmd.flag(self.flag,active_sele+" and "+display_sele,"set")
            self.cmd.flag(self.flag,active_sele+" and not "+display_sele,"clear")
        self.cmd.refresh_wizard()
        self.update_display()

    def get_prompt(self):
        if active_sele not in self.cmd.get_names("selections"):
            return ["Pick object to operate on..."]
        else:
            self.cmd.reference("validate",active_sele) # overbroad
            if self.flag == 2:
                return ["Toggle restrained atoms..."]
            elif self.flag ==3:
                return ["Toggle fixed atoms..."]
            else:
                return ["Toggle unknown atom flag..."]
        
    def toggle(self,flag=0):
        self.flag = flag
        if self.activateOrDismiss():
            if self.activeSeleValid():
                self.update_display()
            else:
                self.cmd.deselect()
            self.cmd.unpick()
                            
    def do_all(self):
        if active_sele in self.cmd.get_names("selections"):
            self.cmd.flag(self.flag,active_sele,"set")
            self.update_display()
        
    def do_less(self,mode):
        if active_sele in self.cmd.get_names("selections"):
            if mode == 0:
                self.cmd.flag(self.flag,"(( byobj " + active_sele + 
                              " ) and not flag %d) extend 1"%self.flag,"clear")
            elif mode == 1:
                self.cmd.flag(self.flag,"byres ((( byobj " + active_sele + 
                              " ) and not flag %d) extend 1)"%self.flag,"clear")
#            elif mode == 2:
#                self.cmd.flag(self.flag,"name ca and byres (byres (( byobj " + active_sele + 
#                              " ) and not flag %d) extend 1)"%self.flag,"clear")
            self.update_display()

    def do_cas(self,mode):
        if active_sele in self.cmd.get_names("selections"):
            if mode == 1:
                self.cmd.flag(self.flag,active_sele,"clear")                
                self.cmd.flag(self.flag,active_sele+" and polymer and name ca","set")                
            elif mode == 0:
                self.cmd.flag(self.flag,active_sele+
                              " and flag %d and polymer and name ca"%self.flag,"set")
                self.cmd.flag(self.flag,active_sele+
                              " and not (polymer and name ca)","clear")                
                
            self.update_display()
            
    def do_more(self,mode):
        if active_sele in self.cmd.get_names("selections"):
            if mode == 0:
                self.cmd.flag(self.flag,active_sele + 
                              " and (flag %d extend 1)"%self.flag,"set")
            elif mode == 1:
                self.cmd.flag(self.flag,"byres ("+ active_sele + 
                              " and (byres flag %d) extend 1)"%self.flag,"set")
                
            elif mode == 2:
                self.cmd.flag(self.flag,"byres ("+ active_sele + 
                              " and flag %d )"%self.flag,"set")
                
            self.update_display()

    def do_none(self):
        if active_sele in self.cmd.get_names("selections"):
            self.cmd.flag(self.flag,active_sele,"clear")
            self.update_display()
    
    def do_store(self):
        if active_sele in self.cmd.get_names("selections"):
            self.cmd.reference("store",active_sele)

    def do_recall(self):
        if active_sele in self.cmd.get_names("selections"):
            self.cmd.reference("recall",active_sele)

    def do_swap(self):
        if active_sele in self.cmd.get_names("selections"):
            self.cmd.reference("swap",active_sele)

    def get_panel(self):
        title = {2:"Restrained Atoms",
                 3:"Fixed Atoms"}.get(self.flag)
        verb = {2:"Restrain", 3:"Fix"}.get(self.flag)

        result = [
            [ 1, title, ''],
            [ 2, "All",'cmd.get_wizard().do_all()'],
            [ 2, "All C-alphas",'cmd.get_wizard().do_cas(1)'],            
            [ 2, "More (byres)",'cmd.get_wizard().do_more(1)'],
            [ 2, "More",'cmd.get_wizard().do_more(0)'],
            [ 2, "Byresidue",'cmd.get_wizard().do_more(2)'],            
            [ 2, "Less", 'cmd.get_wizard().do_less(0)'],
            [ 2, "Less (by residue)", 'cmd.get_wizard().do_less(1)'],
            [ 2, "Only C-alphas",'cmd.get_wizard().do_cas(0)'],                        
            [ 2, "None", 'cmd.get_wizard().do_none()'],
            [ 2, 'Done','cmd.set_wizard()'],
            ]

        if self.flag == 2:
            result[-1:-1] = [
            [ 2, "Store Reference Coords.", 'cmd.get_wizard().do_store()'],
            [ 2, "Recall Reference Coords.", 'cmd.get_wizard().do_recall()'],
            [ 2, "Swap Reference Coords.", 'cmd.get_wizard().do_swap()']]

        return result
 
    def cleanup(self):
        self.cmd.delete(display_sele)
        Wizard.cleanup(self)

class FixAtomWizard(AtomFlagWizard):
    pass

class RestAtomWizard(AtomFlagWizard):
    pass

##############################################################
# base widgets

class GuiFrame(Frame):
    def __init__(self, parent):
        Frame.__init__(self, parent, width=WID*MAX_COL, height=1, borderwidth=0)
        self.grid(sticky=W)
        self.row = 0
        self.col = 1

    def nextRow(self):
        self.row += 1
        self.col = 1

    def nextColumn(self, columnspan=1):
        self.col += 1
        if self.col > MAX_COL:
            self.nextRow()

class GuiLabel:
    def __init__(self, frame, text, width=WID*2):
        frame.nextRow()
        l = Label(frame, text=text, padx=1, pady=1, justify=RIGHT, width=width)
        l.grid(row=frame.row, column=frame.col, sticky=W+E)
        frame.nextColumn()

class GuiButton:
    def __init__(self,frame, text, command, comment="", colspan=1):
        if(text[0:1]=='|'):
            width=3*colspan
        else:
            width=WID*colspan
        b = Button(frame, text=text, padx=1, pady=1, # note pady=0 better with older Xft due to scaling
            width=width, borderwidth=2, relief=RIDGE,
            command = command)
        b.grid(row=frame.row, column=frame.col, columnspan=colspan,
            sticky=W+E)
        balloon = Pmw.Balloon()
        balloon.bind(b, comment)
        frame.nextColumn(colspan)

class GuiImgButton:
    def __init__(self, frame, img, command, comment="", colspan=1):
        b = Button(frame, image=imgDict[img], padx=1, pady=0, height=17,
            width=WID*colspan, borderwidth=2, relief=RIDGE,
            command = command)
        b.img = img
        b.grid(row=frame.row, column=frame.col, columnspan=colspan,
            sticky=S+W+E)
        balloon = Pmw.Balloon()
        balloon.bind(b, comment)
        frame.nextColumn(colspan)

class GuiRadiobutton:
    def __init__(self, frame, text, var, val, command, colspan=1):
        b = Radiobutton(frame, text=text,
            width=WID*colspan, borderwidth=2, indicatoron=0,
            variable=var, value=val, command=command)
        b.grid(row=frame.row, column=frame.col, columnspan=colspan,
            sticky=W+E)
        frame.nextColumn(colspan)


############################################################


class AtomFrame(GuiFrame):
    def __init__(self, parent):
        self.builder = parent.builder
        self.cmd = self.builder.cmd
        GuiFrame.__init__(self, parent)
#        GuiLabel(self, "Atoms")
        GuiButton(self, "H", lambda s=self: s.replace("H",1,1, "Hydrogen"), "Hydrogen")
        GuiButton(self, "C", lambda s=self: s.replace("C",4,4, "Carbon"), "Carbon")
        GuiButton(self, "N", lambda s=self: s.replace("N",4,3, "Nitrogen"), "Nitrogen")
        GuiButton(self, "O", lambda s=self: s.replace("O",4,2, "Oxygen"), "Oxygen")
        GuiButton(self, "P", lambda s=self: s.replace("P",4,3, "Phosphorous"), "Phosphorous")
        GuiButton(self, "S", lambda s=self: s.replace("S",2,2, "Sulfur"), "Sulfur")
        GuiButton(self, "F", lambda s=self: s.replace("F",1,1, "Fluorine"), "Fluorine")
        GuiButton(self, "Cl",lambda s=self: s.replace("Cl",1,1, "Chlorine"), "Chlorine")
        GuiButton(self, "Br",lambda s=self: s.replace("Br",1,1, "Bromine"), "Bromine")
        GuiButton(self, "I", lambda s=self: s.replace("I",1,1, "Iodine"), "Iodine")
        GuiButton(self, "-CF3", lambda s=self: s.grow("trifluoromethane",4,0,"trifluoro"),
                  "Trifluoromethane")
        GuiButton(self, "-OMe", lambda s=self: s.grow("methanol",5,0,"methoxy"), "Methanol")

    def replace(self, atom, geometry, valence, text):
        picked = collectPicked(self.cmd)
        if len(picked):
            self.cmd.select(active_sele,"byobj "+picked[0])
            self.cmd.replace(atom, geometry, valence)
            self.builder.doAutoPick()
        else:
            ReplaceWizard(_self=self.cmd).toggle(atom,geometry,valence,text)

    def grow(self, name, pos, geom, text):
        if "pk1" in self.cmd.get_names("selections"):
            self.cmd.select(active_sele,"byobj pk1")            
            editor.attach_fragment("pk1", name, pos, geom, _self=self.cmd)
            self.builder.doAutoPick()
        else:
            self.cmd.unpick()
            AttachWizard(self.cmd).toggle(name, pos, geom, text)
            
class FragmentFrame(GuiFrame):
    def __init__(self, parent):
        self.builder = parent.builder
        self.cmd = self.builder.cmd
        GuiFrame.__init__(self, parent)

#        GuiLabel(self, "Fragments")

        GuiButton(self, "CH4", lambda s=self: s.grow("methane",1,0,"methyl"), "Methane")
        GuiButton(self, "C=C", lambda s=self: s.grow("ethylene",4,0,"vinyl"), "Enthylene")
        GuiButton(self, "C#C", lambda s=self: s.grow("acetylene",2,0,"alkynl"), "Acetylene")
        GuiButton(self, "C#N", lambda s=self: s.grow("cyanide",2,0,"cyano"), "Cyanide")
        GuiButton(self, "C=O", lambda s=self: s.grow("formaldehyde",2,0,"carbonyl",), "Aldehyde")
        GuiButton(self, "C=OO", lambda s=self: s.grow("formic",4,0,"carboxyl"), "Formic Acid")
        GuiButton(self, "C=ON", lambda s=self: s.grow("formamide",5,0,"C->N amide"), "C->N amide")
        GuiButton(self, "NC=O", lambda s=self: s.grow("formamide",3,1,"N->C amide"), "N->C amide")
        GuiButton(self, "S=O2", lambda s=self: s.grow("sulfone",3,1,"sulfonyl"), "Sulfone")
        GuiButton(self, "P=O3", lambda s=self: s.grow("phosphite",4,0,"phosphoryl"), "Phosphite")
        GuiButton(self, "N=O2", lambda s=self: s.grow("nitro",3,0,"nitro"), "Nitro")
        GuiButton(self, "",None,"Unassigned")

#        GuiLabel(self, "Rings")
        self.nextRow()
        GuiImgButton(self, "cyc3", lambda s=self: s.grow("cyclopropane",4,0,"cyclopropyl"), "Cyclopropane")
        GuiImgButton(self, "cyc4", lambda s=self: s.grow("cyclobutane",4,0,"cyclobutyl"), "Cyclobutane")
        GuiImgButton(self, "cyc5", lambda s=self: s.grow("cyclopentane",5,0,"cyclopentyl"), "Cyclopentane")
        GuiImgButton(self, "cyc6", lambda s=self: s.grow("cyclohexane",7,0,"cyclohexyl"), "Cyclohexane")
        GuiImgButton(self, "cyc7", lambda s=self: s.grow("cycloheptane",8,0,"cycloheptyl"), "Cycloheptane")
        #self.nextRow()
        GuiImgButton(self, "aro5", lambda s=self: s.grow("cyclopentadiene",5,0,"cyclopentadienyl"), "Cyclopentadiene")
        GuiImgButton(self, "aro6", lambda s=self: s.grow("benzene",6,0,"phenyl"), "Benzene")
        GuiImgButton(self, "aro65",lambda s=self: s.grow("indane",12,0,"indanyl"), "Indane")
        GuiImgButton(self, "aro66",lambda s=self: s.grow("napthylene",13,0,"napthyl"), "Napthylene")
        GuiImgButton(self, "aro67",lambda s=self: s.grow("benzocycloheptane",13,0,
                                                         "benzocycloheptyl"), "Benzocycloheptane")
        GuiButton(self, "",None,"Unassigned")
        GuiButton(self, "",None,"Unassigned")

    def grow(self, name, pos, geom, text):
        if "pk1" in self.cmd.get_names("selections"):
            self.cmd.select(active_sele,"byobj pk1")            
            editor.attach_fragment("pk1", name, pos, geom, _self=self.cmd)
            self.builder.doAutoPick()
        else:
            self.cmd.unpick()
            AttachWizard(self.cmd).toggle(name, pos, geom, text)
#            editor.attach_fragment("", name, pos, geom, _self=self.cmd)
#            self.cmd.unpick()
#            self.cmd.select(active_sele,name)
#            self.builder.doAutoPick()            
        
##############################################################

class ModifyFrame(GuiFrame):
    def __init__(self, parent):
        self.builder = parent.builder
        self.cmd = self.builder.cmd
        GuiFrame.__init__(self, parent)
        GuiLabel(self, " Charge", width=6)
        GuiButton(self, "+1", lambda s=self: s.setCharge("1.0","+1"), "Positive Charge")
        GuiButton(self, "0", lambda s=self: s.setCharge("0.0","neutral"), "Neutral Charge")
        GuiButton(self, "-1", lambda s=self: s.setCharge("-1.0","-1"), "Negative Charge")

        l = Label(self, text="Bonds", width=7)
        l.grid(row=self.row, column=self.col, sticky=E)
        self.nextColumn()

        GuiButton(self, "Create", self.createBond, "Create bond between pk1 and pk2")
        GuiButton(self, "Delete", self.deleteBond, "Delete bond between pk1 and pk2")
        GuiButton(self, "Cycle", self.cycleBond, "Cycle bond valence");
        #self.nextRow()
        GuiButton(self, "|", lambda s=self: s.setOrder("1", "single"), "Create single bond")
        GuiButton(self, "||", lambda s=self: s.setOrder("2", "double", ), "Create double bond")
        GuiButton(self, "|||", lambda s=self: s.setOrder("3", "triple"), "Create triple bond")
#        GuiButton(self, "Cycle", lambda s=self: s.cmd.cycle_valence(1), "Cycle valence (single -> double -> triple")
        GuiButton(self, "Arom", lambda s=self: s.setOrder("4", "aromatic"), "Create aromatic bond")

    def setCharge(self, charge, text):
        picked = collectPicked(self.cmd)
        if len(picked)>0:
            for sele in picked:
                self.cmd.alter(sele,"formal_charge=%s" % charge)
                self.cmd.h_fill()
                if abs(float(charge))>0.0001:
                    self.cmd.label(sele,"'''"+text+"'''")
                else:
                    self.cmd.label(sele)
            self.cmd.unpick()
        else:
            ChargeWizard(self.cmd).toggle(charge, text)
                
    def createBond(self):
        picked = collectPicked(self.cmd)
        if picked == ["pk1","pk2"]:
            if ((self.cmd.count_atoms("pk1 and hydro") == self.cmd.count_atoms("pk1")) and
                (self.cmd.count_atoms("(pk1 extend 1) and not hydro")==1) and
                (self.cmd.count_atoms("pk2 and hydro") == self.cmd.count_atoms("pk2")) and
                (self.cmd.count_atoms("(pk2 extend 1) and not hydro")==1)):
                self.cmd.select("pk1","(pk1 extend 1) and not hydro")
                self.cmd.select("pk2","(pk2 extend 1) and not hydro")
            self.cmd.bond("pk1", "pk2")
            self.cmd.h_fill()
            self.cmd.unpick()
        else:
            BondWizard(self.cmd).toggle()


    def deleteBond(self):
        picked = collectPicked(self.cmd)
        if picked == ["pk1","pk2"]:
            self.cmd.unbond("pk1", "pk2")
            self.cmd.h_fill()
            self.cmd.unpick()
        else:
            self.cmd.unpick()
            UnbondWizard(self.cmd).toggle()

    def cycleBond(self):
        picked = collectPicked(self.cmd)
        if picked == ["pk1","pk2"]:
            self.cmd.cycle_valence()
            self.cmd.unpick()
        else:
            ValenceWizard(_self=self.cmd).toggle(-1,"Cycle bond")        

    def setOrder(self, order, text):
        picked = collectPicked(self.cmd)
        if picked == ["pk1","pk2"]:
            self.cmd.unbond("pk1", "pk2")
            self.cmd.bond("pk1", "pk2", order)
            self.cmd.h_fill()
            self.cmd.unpick()
        else:
            self.cmd.unpick()
            ValenceWizard(_self=self.cmd).toggle(order,text)        
class EditFrame(GuiFrame):
    def __init__(self, parent):
        self.builder = parent.builder
        self.cmd = self.builder.cmd
        GuiFrame.__init__(self, parent)

        GuiLabel(self, " Atoms", width=6)
        GuiButton(self, "Fix H", self.fixH, "Fix hydrogens on picked atoms.")
        GuiButton(self, "Add H", self.addH, "Add hydrogens to entire molecule")
        GuiButton(self, "Invert", self.invert, "Invert stereochemistry around pk1 (pk2 and pk3 will remain fixed)")
        #self.nextRow()
#        GuiButton(self, "Center", self.center, "Center pk1 (or molecule)")
        GuiButton(self, "Delete", self.removeAtom, "Remove atoms")
        #GuiButton(self, "Reset", self.reset)
        GuiButton(self, "Clear", self.clear, "Delete everything")
        l = Label(self, text="Model", width=5)
        l.grid(row=self.row, column=self.col, sticky=E)
        self.nextColumn()
        GuiButton(self, "Clean", self.clean, "Cleanup Structure")
        GuiButton(self, "Sculpt", self.sculpt, "Molecular Sculpting")
        GuiButton(self, "Fix", self.fix, "Fix Atom Positions")
        GuiButton(self, "Rest", self.rest, "Restrain Atom Positions")
        GuiButton(self, "Undo", self.undo, "Undo Changes")

    def fixH(self):
        picked = collectPicked(self.cmd)
        if len(picked):
            self.cmd.h_fill()
            self.cmd.unpick()
        else:
            HydrogenWizard(_self=self.cmd).toggle('fix')
        
    def addH(self):
        picked = collectPicked(self.cmd)
        if len(picked):
            self.cmd.h_add("pkmol")
            self.cmd.unpick()            
        else:
            HydrogenWizard(_self=self.cmd).toggle('add')

    def invert(self):
        picked = collectPicked(self.cmd)
        if picked == ["pk1","pk2","pk3"]:
            self.cmd.invert()
            self.cmd.unpick()
        else:
            self.cmd.unpick()
            InvertWizard(self.cmd).toggle()

    def center(self):
        if "pk1" in self.cmd.get_names("selections"):
            self.cmd.zoom("pk1", 5.0, animate=-1)
        else:
            self.cmd.zoom("all", 3.0, animate=-1)

    def removeAtom(self):
        picked = collectPicked(self.cmd)
        if len(picked):
            if self.cmd.count_atoms("?pkbond"):
                self.cmd.edit("(pk1)","(pk2)",pkbond=0)
            cnt = self.cmd.select(active_sele,
                   "(((?pkset or ?pk1) and not hydro) extend 1) and not hydro")
            self.cmd.remove_picked()
            if cnt:
                self.cmd.fix_chemistry(active_sele)
                self.cmd.h_add(active_sele)
            self.cmd.delete(active_sele)
            self.cmd.unpick()
        else:
            RemoveWizard(self.cmd).toggle()
            
    def reset(self):
        self.cmd.unpick()

    def clear(self):
        check = tkMessageBox.askokcancel("Confirm", "Really delete everything?")
        if check:
            self.cmd.delete("all")

    def sculpt(self):
        picked = collectPicked(self.cmd)
        if len(picked):
            self.cmd.select(active_sele, string.join(picked," or "))
        SculptWizard(_self=self.cmd).toggle()

    def clean(self):
        picked = collectPicked(self.cmd)
        if len(picked):
            self.cmd.select(active_sele, "pkmol")
            self.cmd.unpick()
        CleanWizard(_self=self.cmd).toggle()
        
    def undo(self):
        warn("Sorry, the undo button is not yet implemented.")
        
    def fix(self):
        picked = collectPicked(self.cmd)
        if len(picked):
            self.cmd.select(active_sele,"pk1")
            self.cmd.deselect()
        else:
            self.cmd.delete(active_sele)
        FixAtomWizard(_self=self.cmd).toggle(3)

    def rest(self):
        picked = collectPicked(self.cmd)
        if len(picked):
            self.cmd.select(active_sele,"byobj ("+string.join(picked," or ")+")")
            self.cmd.deselect()
        else:
            self.cmd.delete(active_sele)
        RestAtomWizard(_self=self.cmd).toggle(2)
            
############################################################

class AminoAcidFrame(GuiFrame):
    def attach(self, aa):
        picked = collectPicked(self.cmd)
        if len(picked)==1:
            try:
                editor.attach_amino_acid(picked[0], aa, _self=self.cmd)
            except:
                traceback.print_exc()
            self.builder.doZoom()
        else:
            self.cmd.unpick()
            AminoAcidWizard(_self=self.cmd).toggle(aa)    
        
    def __init__(self, parent):
        self.builder = parent.builder
        self.cmd = self.builder.cmd
        GuiFrame.__init__(self, parent)
        #GuiLabel(self, "Amino Acids")

        aaList = ["Ace", "Ala", "Arg", "Asn", "Asp", "Cys", "Gln", "Glu", "Gly", 
                  "His", "Ile", "Leu", "Lys", "Met", "Phe", "Pro", "Ser", 
                  "Thr", "Trp", "Tyr", "Val", "NMe", "NHH", ]
        for aa in aaList:
            r = aa.lower()
            GuiButton(self, aa, lambda s=self,r=r: s.attach(r), "Build %s residue" % aa)


class SecStructFrame(GuiFrame):
    def __init__(self, parent):
        self.builder = parent.builder
        self.cmd = self.builder.cmd
        GuiFrame.__init__(self, parent)
        self.secStructL = ["alpha helix", "parallel beta sheet",
                           "antiparallel beta sheet"]
        self.ss = int(float(self.cmd.get("secondary_structure")))-1
        ss = self.secStructL[self.ss]

        GuiLabel(self, "Secondary Structure:", WID*4)
        GuiButton(self, "Next...", self.toggleSS, "Toggle secondary structure for building")
        self.ssText = Label(self, text = ss, justify = LEFT)
        self.ssText.grid(row=self.row, column=self.col, columnspan=4, sticky=W)

    def toggleSS(self):
        self.ss += 1
        if self.ss == len(self.secStructL):
            self.ss = 0
        self.ssText.configure(text=self.secStructL[self.ss])
        self.cmd.set("secondary_structure", "%d" % (self.ss+1))


############################################################

class ChemFrame(Frame):
    def __init__(self, parent):
        """module for constructing chemical gui & widgets for molecule
        building and editing with PyMol
        """
        self.builder = parent.builder
        Frame.__init__(self, parent, width=WID*MAX_COL, relief=SUNKEN, bd=2)
        AtomFrame(self)
        FragmentFrame(self)

class ProteinFrame(Frame):
    def __init__(self, parent):
        """module for constructing tkinter gui & widgets for molecule
        building and editing with PyMol
        """
        self.builder = parent.builder
        Frame.__init__(self, parent, width=WID*MAX_COL, relief=SUNKEN, bd=2)
        AminoAcidFrame(self)
        SecStructFrame(self)

class EditorFrame(Frame):
    def __init__(self, parent):
        """module for constructing chemical gui & widgets for molecule
        building and editing with PyMol
        """
        self.builder = parent.builder
        Frame.__init__(self, parent, width=WID*MAX_COL)  #, relief=GROOVE, bd=1)
        ModifyFrame(self)
        EditFrame(self)

############################################################

class Builder(Frame):
    def __init__(self, app, parent, *kw, **args):
        """
        module for constructing tkinter gui & widgets for molecule
        building and editing with PyMol

        Builder constructs the top level gui that will contain
        other widgets for working with chemicals, proteins, etc.
        """
        self.builder = self
        self.app = app
        self.cmd = app.pymol.cmd
        
        Frame.__init__(self, parent, *kw, **args)
        self.deferred = 1
        
    def deferred_activate(self):

        # avoids incurring the overhead of launching the builder unless
        # we need it
        
        if self.deferred:
            self.deferred = 0
            
            # get bitmaps for buttons
            #        imgDir = os.path.split(__file__)[0] + "/bitmaps/"
            
            imgDir = os.path.join(os.environ['PYMOL_DATA'], "pmg_tk/bitmaps/builder")
#            print imgDir
            imgList = glob("%s/aro*.gif" % imgDir) + glob("%s/cyc*.gif" % imgDir)
            for imgFile in imgList:
                imgName = os.path.splitext(os.path.split(imgFile)[1])[0]
                if imgName not in imgDict.keys():
                    imgDict[imgName] = PhotoImage(file=imgFile)

            # construct everything
            self.constructMain()

        # unsafe approach disabled
        # self.bind_all("<ButtonRelease-1>", self.doSelect)

    def doAutoPick(self, old_atoms=None):
        if self.autoPik.get():
            self.cmd.unpick()
            if self.cmd.select(newest_sele,"(byobj "+active_sele+") and not "+active_sele)==0:
                self.cmd.select(newest_sele, active_sele)
            new_list = self.cmd.index(newest_sele+" and hydro")
            if len(new_list)==0:
                new_list = self.cmd.index(newest_sele)                
            if new_list:
                index = new_list.pop()
                try:
                    self.cmd.edit("%s`%d" % index)
                    if self.cmd.get_wizard()!=None:
                        self.cmd.do("_ cmd.get_wizard().do_pick(0)")
                except pymol.CmdException:
                    print " doAutoPick-Error: exception"
            self.doZoom()

    def doZoom(self, *ignore):
#        print "zoom",self.autoZoom.get()
        if self.autoZoom.get():
            if "pk1" in self.cmd.get_names("selections"):
                self.cmd.zoom("((neighbor pk1) extend 4)", 4.0, animate=-1)

#    def doValence(self, *ignore):
#        self.cmd.set("valence", self.showValence.get())

    def doVdw(self, *ignore):
        self.cmd.set("sculpt_vdw_vis_mode", self.showVdw.get())

    def doElectro(self, *ignore):
        self.cmd.set("clean_electro_mode", self.electroMode.get())


##############################################################
# main constructor methods

    def constructMain(self):
        self.chemB = Button(self, text="Chemical", 
            relief=SUNKEN, command = self.toggleChemProtein)
        self.chemB.grid(row=1, column=1, sticky=N+S+W+E)

        self.protB = Button(self, text="Protein", 
            command = self.toggleChemProtein)
        self.protB.grid(row=2, column=1, sticky=N+S+W+E)

        self.protFrame = ProteinFrame(self)
        self.chemFrame = ChemFrame(self)
        self.editFrame = EditorFrame(self)
        self.editFrame.grid(row=5, column=2, rowspan=2, sticky=N+S+W+E)
        self.toggleChemProtein()

        self.autoPik = IntVar()
        self.autoPik.set(0)
#        autoB = Checkbutton(self, text="autopick", 
#            borderwidth=1, pady=0, justify=LEFT, variable=self.autoPik, 
#            onvalue=1, offvalue=0, command=self.cmd.unpick)
#        autoB.grid(row=4, column=1, sticky=W)

        self.autoZoom = IntVar()
        self.autoZoom.set(0)
#        Checkbutton(self, text="autozoom", 
#            borderwidth=1, pady=0, justify=LEFT, variable=self.autoZoom, 
#            onvalue=1, offvalue=0, command=self.doZoom).grid(row=5, column=1,
#            sticky=W)

        self.showValence = StringVar()
        self.showValence.set(self.cmd.get("valence"))
#        Checkbutton(self, text="valence", 
#            borderwidth=1, pady=0, justify=LEFT, variable=self.showValence, 
#            onvalue="on", offvalue="off", command=self.doValence).grid(row=6, 
#            column=1, sticky=W)

        self.showVdw = StringVar()
        self.showVdw.set(self.cmd.get("sculpt_vdw_vis_mode"))
        Checkbutton(self, text="Bumps", 
            borderwidth=1, pady=0, justify=LEFT, variable=self.showVdw, 
            onvalue="on", offvalue="off", command=self.doVdw).grid(row=5, 
            column=1, sticky=W)

        self.electroMode = StringVar()
        self.electroMode.set(int(self.cmd.get("clean_electro_mode")))
        Checkbutton(self, text="El-stat", 
            borderwidth=1, pady=0, justify=LEFT, variable=self.electroMode, 
            onvalue=1, offvalue=0, command=self.doElectro).grid(row=6, 
            column=1,sticky=W)

    def toggleChemProtein(self):
        if self.chemFrame.grid_info():
            self.chemB.configure(relief=RAISED)
            self.chemFrame.grid_forget()
            self.protB.configure(relief=SUNKEN)
            self.protFrame.grid(row=1, column=2, rowspan=4, sticky=W)
        else:
            self.chemB.configure(relief=SUNKEN)
            self.chemFrame.grid(row=1, column=2, rowspan=4, sticky=W)
            self.protB.configure(relief=RAISED)
            self.protFrame.grid_forget()


############################################################

def getSeleDict(self_cmd):
    result = {}
    for sele in self_cmd.get_names("selections"):
        result[sele] = 1
    return result

def collectPicked(self_cmd):
    result = []
    sele_dict = getSeleDict(self_cmd)
    for sele in ["pk1","pk2","pk3","pk4"]:
        if sele_dict.has_key(sele):
            result.append(sele)
    return result

def getAtoms(self_cmd, nAtom):
    """counts how many atoms are selected by (pk) selections
    @param nAtom: number of atoms to select
    @type nAtom: integer
    @returns: 1 if number selected == nAtoms, 0 if not
    @rtype: integer
    """
    count = 0
    for i in range(1,5):
        if "pk%d" % i in self_cmd.get_names("selections"):
            count += 1

    if count == nAtom:
        return 1
    else:
        msg = "Please select (pk1) first..."
        for w in range(2, nAtom+1):
            msg += " + (pk%d)" % w
        warn(msg)
        return 0

def warn(text):
    check = tkMessageBox.showerror("Error", text)


