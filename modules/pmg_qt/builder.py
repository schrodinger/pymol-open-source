
import itertools, os
from glob import glob

import pymol
import pymol._gui

from .pymol_gl_widget import PyMOLGLWidget
from pymol.Qt import QtGui, QtCore
from pymol.Qt import QtWidgets
from pymol.Qt.utils import PopupOnException

Qt = QtCore.Qt

from pymol.wizard import Wizard
from pymol.parsing import QuietException

from pymol import editor
from pymol import computing

from chempy import cpv

active_sele = "_builder_active" # object we're working on...
newest_sele = "_builder_added" # last atom added?
display_sele = "_build_display"

undocontext = editor.undocontext

def undoablemethod(sele):
    def decorator(func):
        def wrapper(self, *args, **kwargs):
            with undocontext(self.cmd, sele):
                return func(self, *args, **kwargs)
        return wrapper
    return decorator

#############################################################
# Action-directing Wizards

class ActionWizard(Wizard):

    def __init__(self, _self=pymol.cmd):
        Wizard.__init__(self,_self)
        self.actionHash = str(self.__class__)

    def setActionHash(self, action_hash):
        self.actionHash = action_hash

    def activateOrDismiss(self):
        activate_flag = 1
        cur_wiz = self.cmd.get_wizard()
        if cur_wiz is not None:
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
            else:
                enabled_list = self.cmd.get_names("objects",enabled_only=1)
                active_obj_list = self.cmd.get_object_list(active_sele)
                if len(active_obj_list) != 1:
                    self.cmd.delete(active_sele)
                elif active_obj_list[0] not in enabled_list:
                    self.cmd.delete(active_sele)
        if "pk1" in self.cmd.get_names("selections"):
            self.cmd.select(active_sele,"byobj pk1")
        else:
            enabled_list = self.cmd.get_names("objects",enabled_only=1)
            if len(enabled_list)==1:
                if self.cmd.select(active_sele, enabled_list[0])<1:
                    self.cmd.delete(active_sele)
        return active_sele in self.cmd.get_names("selections")


class CleanWizard(ActionWizard):

    def __init__(self, _self=pymol.cmd):
        self.clean_obj = None
        ActionWizard.__init__(self,_self)

    def run_job(self):
        if active_sele in self.cmd.get_names("selections"):
            obj_list = self.cmd.get_object_list(active_sele)
            if len(obj_list)==1:
                self.cmd.unpick()
                self.cmd.set_wizard()
                self.cmd.refresh_wizard()
                self.cmd.do("_ cmd.clean('%s',message='''Cleaning %s...''',async_=1)"%(active_sele,obj_list[0]))

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
        if isinstance(obj_list,list) and (len(obj_list)==1):
            self.run_job()
        else:
            print("Error: can only clean one object at a time")

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

    def __init__(self, _self=pymol.cmd):
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
                print("Error: cannot sculpt more than one object at a time")

    def sculpt_deactivate(self):
        if ((self.sculpt_object is not None) and
            self.sculpt_object in self.cmd.get_names()):
            self.cmd.set("sculpt_vdw_vis_mode","0",self.sculpt_object)
            self.cmd.sculpt_iterate(self.sculpt_object,self.cmd.get_state(),0)
            self.cmd.unset("sculpt_vdw_vis_mode",self.sculpt_object)
            self.cmd.sculpt_deactivate(self.sculpt_object)
            self.sculpt_object = None
            self.cmd.refresh_wizard()

    def do_pick(self, bondFlag):
        if self.sculpt_object is None:
            self.cmd.select(active_sele, "byobj pk1")
            self.sculpt_activate()
        else:
            return 0 # already sculpting, so handle like a normal edit

    def toggle(self):
        if self.activateOrDismiss():
            if self.activeSeleValid():
                self.sculpt_activate()

    def get_prompt(self):
        if self.sculpt_object is None:
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

    def __init__(self, _self=pymol.cmd):
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
        if cur_wiz is not None:
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

    def __init__(self, _self=pymol.cmd):
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

    def __init__(self, _self=pymol.cmd, ss=-1):
        RepeatableActionWizard.__init__(self,_self)
        self.mode = 0
        self.setSecondaryStructure(ss)

    def setSecondaryStructure(self, ss):
        self._secondary_structure = ss

    def attach_monomer(self, objectname=""):
         editor.attach_amino_acid("?pk1", self.aminoAcid, object=objectname,
                ss=self._secondary_structure,
                 _self=self.cmd)

    def do_pick(self, bondFlag):
        # since this function can change any position of atoms in a related
        # molecule, bymol is used
        if self.mode == 0:
            self.cmd.select(active_sele, "bymol pk1")
            try:
                with undocontext(self.cmd, "bymol ?pk1"):
                    self.attach_monomer(self.aminoAcid)
            except QuietException:
                fin = -1
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
        self.attach_monomer(self.aminoAcid)

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

    @undoablemethod("(?pk1 ?pk2) extend 1")
    def do_pick(self, bondFlag):
        self.cmd.select(active_sele, "bymol pk1")
        if bondFlag:
            if int(self.order)>=0:
                self.cmd.valence(self.order, "pk1", "pk2")
                self.cmd.h_fill()
            else:
                self.cmd.cycle_valence()
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

    @undoablemethod("bymol ?pk1")
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

    @PopupOnException.decorator
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

    @staticmethod
    def staticaction(cmd):
        picked = collectPicked(cmd)
        if picked != ["pk1","pk2"]:
            return False

        cmd.select(active_sele, "bymol ?pk1")

        if (    cmd.count_atoms("?pk1&hydro") and
                cmd.count_atoms("?pk2&hydro") and
                cmd.count_atoms("(?pk1 extend 1)&!hydro") and
                cmd.count_atoms("(?pk2 extend 1)&!hydro")):
            # two hydrogens picked -> pick their heavy neighbors instead
            cmd.select("pk1","(pk1 extend 1) and not hydro")
            cmd.select("pk2","(pk2 extend 1) and not hydro")

        with undocontext(cmd, "(?pk1 ?pk2) extend 1"):
            cmd.bond("pk1", "pk2")
            cmd.h_fill()

        cmd.unpick()
        return True

    def do_pick(self, bondFlag):
        if self.staticaction(self.cmd):
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

    @undoablemethod("(?pk1 ?pk2) extend 1")
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

    @undoablemethod("?pk1 extend 1")
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


#############################################################
### Helper functions

def getSeleDict(self_cmd):
    result = {}
    for sele in self_cmd.get_names("selections"):
        result[sele] = 1
    return result


def collectPicked(self_cmd):
    result = []
    sele_dict = getSeleDict(self_cmd)
    for sele in ["pk1","pk2","pk3","pk4"]:
        if sele in sele_dict:
            result.append(sele)
    return result


#############################################################
### Actual GUI


def makeFragmentButton():
    btn = QtWidgets.QPushButton()
    btn.setAttribute(Qt.WA_LayoutUsesWidgetRect) # OS X workaround
    btn.setSizePolicy(
            QtWidgets.QSizePolicy.Minimum,
            QtWidgets.QSizePolicy.MinimumExpanding)
    btn.setAutoDefault(False)
    return btn


class _BuilderPanel(QtWidgets.QWidget):

    def __init__(self, parent=None, app=None):
        super(_BuilderPanel, self).__init__(parent)

        self.setWindowTitle("Builder")
        self.setObjectName("builder")
        self.cmd = app.pymol.cmd

        self.layout = QtWidgets.QVBoxLayout()
        self.setLayout(self.layout)
        self.buttons_layout = QtWidgets.QVBoxLayout()

        self.tabs = QtWidgets.QTabWidget(self)
        self.layout.setContentsMargins(5, 5, 5, 5);
        self.layout.setSpacing(5);
        self.layout.addWidget(self.tabs)
        self.layout.addLayout(self.buttons_layout)
        self.layout.addStretch()

        self.fragments_layout = QtWidgets.QGridLayout()
        self.fragments_layout.setContentsMargins(5, 5, 5, 5);
        self.fragments_layout.setSpacing(5);
        self.fragments_tab = QtWidgets.QWidget()
        self.fragments_tab.setLayout(self.fragments_layout)
        self.protein_layout = QtWidgets.QGridLayout()
        self.protein_layout.setContentsMargins(5, 5, 5, 5);
        self.protein_layout.setSpacing(5);
        self.protein_tab = QtWidgets.QWidget()
        self.protein_tab.setLayout(self.protein_layout)

        self.tabs.addTab(self.fragments_tab, "Chemical")
        self.tabs.addTab(self.protein_tab, "Protein")

        self.getIcons()

        buttons = [
            [ ("H", "Hydrogen", lambda: self.replace("H", 1, 1, "Hydrogen")),
              ("C", "Carbon", lambda: self.replace("C", 4, 4, "Carbon")),
              ("N", "Nitrogen", lambda: self.replace("N", 4, 3, "Nitrogen")),
              ("O", "Oxygen", lambda: self.replace("O", 4, 2, "Oxygen")),
              ("P", "Phosphorus", lambda: self.replace("P",4,3, "Phosphorous")),
              ("S", "Sulfur", lambda: self.replace("S",2,2, "Sulfur")),
              ("F", "Fluorine", lambda: self.replace("F",1,1, "Fluorine")),
              ("Cl", "Chlorrine", lambda: self.replace("Cl",1,1, "Chlorine")),
              ("Br", "Bromine", lambda: self.replace("Br",1,1, "Bromine")),
              ("I", "Iodine", lambda: self.replace("I",1,1, "Iodine")),
              ("-CF3", "Trifluoromethane", lambda: self.grow("trifluoromethane",4,0, "trifluoro")),
              ("-OMe", "Methanol", lambda: self.grow("methanol",5,0, "methoxy")),
            ],
            [ ("CH4", "Methyl", lambda: self.grow("methane",1,0,"methyl")),
              ("C=C", "Ethylene", lambda: self.grow("ethylene",4,0,"vinyl")),
              ("C#C", "Acetylene", lambda: self.grow("acetylene",2,0,"alkynl")),
              ("C#N", "Cyanide", lambda: self.grow("cyanide",2,0,"cyano")),
              ("C=O", "Aldehyde", lambda: self.grow("formaldehyde",2,0,"carbonyl",)),
              ("C=OO", "Formic Acid", lambda: self.grow("formic",4,0,"carboxyl")),
              ("C=ON", "C->N amide", lambda: self.grow("formamide",5,0,"C->N amide")),
              ("NC=O", "N->C amide", lambda: self.grow("formamide",3,1,"N->C amide")),
              ("S=O2", "Sulfone", lambda: self.grow("sulfone",3,1,"sulfonyl")),
              ("P=O3", "Phosphite", lambda: self.grow("phosphite",4,0,"phosphoryl")),
              ("N=O2", "Nitro", lambda: self.grow("nitro",3,0,"nitro")),
            ],
            [
              ("#cyc3", "Cyclopropane", lambda: self.grow("cyclopropane",4,0,"cyclopropyl")),
              ("#cyc4", "Cyclobutane", lambda: self.grow("cyclobutane",4,0,"cyclobutyl")),
              ("#cyc5", "Cyclopentane", lambda: self.grow("cyclopentane",5,0,"cyclopentyl")),
              ("#cyc6", "Cyclohexane", lambda: self.grow("cyclohexane",7,0,"cyclohexyl")),
              ("#cyc7", "Cycloheptane", lambda: self.grow("cycloheptane",8,0,"cycloheptyl")),
              ("#aro5", "Cyclopentadiene", lambda: self.grow("cyclopentadiene",5,0,"cyclopentadienyl")),
              ("#aro6", "Benzene", lambda: self.grow("benzene",6,0,"phenyl")),
              ("#aro65", "Indane", lambda: self.grow("indane",12,0,"indanyl")),
              ("#aro66", "Napthylene", lambda: self.grow("napthylene",13,0,"napthyl")),
              ("#aro67", "Benzocycloheptane", lambda: self.grow("benzocycloheptane",13,0, "benzocycloheptyl")),
            ]
        ]

        self.btn_icons = {}

        requestsize = QtCore.QSize(48, 48)
        for row, btn_row in enumerate(buttons):
            for col, bb in enumerate(btn_row):
                btn_label, btn_tooltip, btn_command = bb
                btn = makeFragmentButton()
                if btn_label.startswith('#'):
                    icons = self.icons[btn_label[1:]]
                    btn.setIcon(icons[0])
                    btn.setIconSize(icons[1].actualSize(requestsize))
                    self.btn_icons[btn] = icons
                else:
                    btn.setText(btn_label)
                btn.setToolTip(btn_tooltip)
                btn.clicked.connect(btn_command)
                self.fragments_layout.addWidget(btn, row, col)

        buttons = [
            [ 'Ace', 'Ala', 'Arg', 'Asn', 'Asp', 'Cys', 'Gln', 'Glu', 'Gly', 'His', 'Ile', 'Leu' ],
            [ 'Lys', 'Met', 'Phe', 'Pro', 'Ser', 'Thr', 'Trp', 'Tyr', 'Val', 'NMe', 'NHH' ]
        ]
        for row, btn_row in enumerate(buttons):
            for col, btn_label in enumerate(btn_row):
                btn = makeFragmentButton()
                btn.setText(btn_label)
                btn.setToolTip("Build %s residue" % btn_label)
                res = btn_label.lower()
                slot = lambda val=None, s=self,r=res: s.attach(r)
                btn.clicked.connect(slot)
                self.protein_layout.addWidget(btn, row, col)

        lab = QtWidgets.QLabel('Secondary Structure:')
        lab_cols = 3
        self.ss_cbox = QtWidgets.QComboBox()
        self.ss_cbox.addItem("Alpha Helix")
        self.ss_cbox.addItem("Beta Sheet (Anti-Parallel)")
        self.ss_cbox.addItem("Beta Sheet (Parallel)")
        self.protein_layout.addWidget(lab, 2, 0, 1, lab_cols)
        self.protein_layout.addWidget(self.ss_cbox, 2, lab_cols, 1, 4)
        self.ss_cbox.currentIndexChanged[int].connect(self.ssIndexChanged)

        buttons = [
            [
              ( "@Atoms:", None, None),
              ( "Fix H", "Fix hydrogens on picked atoms", self.fixH),
              ( "Add H", "Add hydrogens to entire molecule", self.addH),
              ( "Invert", "Invert stereochemistry around pk1 (pk2 and pk3 will remain fixed)", self.invert),
              ( "Delete", "Remove atoms", self.removeAtom),
              ( "Clear", "Delete everything", self.clear),
              ( "@   Charge:", None, None),
              ( " +1 ", "Positive Charge", lambda: self.setCharge(1,"+1")),
              ( "  0 ", "Neutral Charge", lambda: self.setCharge(0,"neutral")),
              ( " -1 ", "Negative Charge", lambda: self.setCharge(-1,"-1")),
            ],
            [
              ( "@Bonds:", None, None),
              ( "Create", "Create bond between pk1 and pk2", self.createBond),
              ( "Delete", "Delete bond between pk1 and pk2", self.deleteBond),
              ( "Cycle", "Cycle bond valence", self.cycleBond),
              ( "  |  ", "Create single bond", lambda: self.setOrder("1", "single")),
              ( " || ", "Create double bond", lambda: self.setOrder("2", "double")),
              ( " ||| ", "Create triple bond", lambda: self.setOrder("3", "triple")),
              ( "Arom", "Create aromatic bond", lambda: self.setOrder("4", "aromatic")),
              ( "@   Model:", None, None),
              ( "Clean", "Cleanup structure", self.clean),
              ( "Sculpt", "Molecular sculpting", self.sculpt),
              ( "Fix", "Fix atom positions", self.fix),
              ( "Rest", "Restrain atom positions", self.rest),
            ],
            [
              ( "$El-stat", "Electrostatics term for 'Clean' action", "clean_electro_mode"),
              ( "@   ", None, None),
              ( "$Bumps", "Show VDW contacts during sculpting", "sculpt_vdw_vis_mode"),
              ( "@   ", None, None),
              ( "#Undo Enabled", "", "suspend_undo"),
              ( "Undo", "Undo last change", self.undo),
              ( "Redo", "Redo last change", self.redo),
            ]
        ]

        for row, btn_row in enumerate(buttons):
            btn_row_layout = QtWidgets.QHBoxLayout()
            self.buttons_layout.addLayout(btn_row_layout)
            for col, bb in enumerate(btn_row):
                btn_label, btn_tooltip, btn_command = bb
                if btn_label[0] == '@':
                    btn = QtWidgets.QLabel(btn_label[1:])
                elif btn_label[0] in ('#', '$'):
                    btn = QtWidgets.QCheckBox(btn_label[1:])
                    setting = btn_command
                    value = self.cmd.get_setting_int(setting)
                    if btn_label[0] == '$':
                        btn.setChecked(bool(value))
                        @btn.toggled.connect
                        def _(checked, n=setting):
                            self.cmd.set(n, checked, quiet=0)
                    else:
                        btn.setChecked(not value)
                        @btn.toggled.connect
                        def _(checked, n=setting):
                            self.cmd.set(n, not checked, quiet=0)
                else:
                    btn = makeFragmentButton()
                    btn.setText(btn_label)
                    btn.clicked.connect(btn_command)
                if btn_tooltip:
                    btn.setToolTip(btn_tooltip)
                btn_row_layout.addWidget(btn)
            btn_row_layout.addStretch()

    def getIcons(self):
        self.icons = {}
        # use old Tk icons
        imgDir = os.path.join(os.environ['PYMOL_DATA'], "pmg_tk/bitmaps/builder")
        imgList = glob("%s/aro*.gif" % imgDir) + glob("%s/cyc*.gif" % imgDir)
        for imgFile in imgList:
            imgName = os.path.splitext(os.path.split(imgFile)[1])[0]
            if imgName not in list(self.icons.keys()):
                image = QtGui.QImage(imgFile)
                pixmap = QtGui.QPixmap.fromImage(image)
                image.invertPixels()
                inv_pixmap = QtGui.QPixmap.fromImage(image)
                self.icons[imgName] = (QtGui.QIcon(pixmap), QtGui.QIcon(inv_pixmap))

    def showEvent(self, event):
        self.cmd.set("editor_auto_measure", 0)
        self.cmd.set("auto_overlay")
        self.cmd.set("valence")
        self.cmd.edit_mode(1)

    def grow(self, name, pos, geom, text):
        if "pk1" in self.cmd.get_names("selections"):
            self.cmd.select(active_sele,"byobj pk1")
            editor.attach_fragment("pk1", name, pos, geom, _self=self.cmd)
            self.doAutoPick()
        else:
            self.cmd.unpick()
            AttachWizard(self.cmd).toggle(name, pos, geom, text)

    def replace(self, atom, geometry, valence, text):
        picked = collectPicked(self.cmd)
        if len(picked):
            self.cmd.select(active_sele,"byobj "+picked[0])
            self.cmd.replace(atom, geometry, valence)
            self.doAutoPick()
        else:
            ReplaceWizard(_self=self.cmd).toggle(atom,geometry,valence,text)

    def attach(self, aa):
        ss = self.ss_cbox.currentIndex() + 1
        picked = collectPicked(self.cmd)
        if len(picked)==1:
            try:
                with undocontext(self.cmd, "bymol %s" % picked[0]):
                    editor.attach_amino_acid(picked[0], aa,
                            ss=ss, _self=self.cmd)
            except:
                fin = -1
            self.doZoom()
        else:
            self.cmd.unpick()
            AminoAcidWizard(_self=self.cmd, ss=ss).toggle(aa)

    def ssIndexChanged(self, index):
        w = self.cmd.get_wizard()
        if isinstance(w, AminoAcidWizard):
            w.setSecondaryStructure(index + 1)

    def doAutoPick(self, old_atoms=None):
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
                if self.cmd.get_wizard() is not None:
                    self.cmd.do("_ cmd.get_wizard().do_pick(0)")
            except pymol.CmdException:
                print(" doAutoPick-Error: exception")
        self.doZoom()

    def doZoom(self, *ignore):
        if "pk1" in self.cmd.get_names("selections"):
            self.cmd.zoom("((neighbor pk1) extend 4)", 4.0, animate=-1)

    def setCharge(self, charge, text):
        picked = collectPicked(self.cmd)
        if len(picked)>0:
            sele = "?pk1 ?pk2 ?pk3 ?pk4"
            with undocontext(self.cmd, sele):
                self.cmd.alter(sele,"formal_charge=%s" % charge)
                self.cmd.h_fill()
                self.cmd.label(sele,'"%+d" % formal_charge if formal_charge else ""')
            self.cmd.unpick()
        else:
            ChargeWizard(self.cmd).toggle(charge, text)

    def createBond(self):
        if not BondWizard.staticaction(self.cmd):
            BondWizard(self.cmd).toggle()

    def deleteBond(self):
        picked = collectPicked(self.cmd)
        if picked == ["pk1","pk2"]:
            with undocontext(self.cmd, "(?pk1 ?pk2) extend 1"):
                self.cmd.unbond("pk1", "pk2")
                self.cmd.h_fill()
            self.cmd.unpick()
        else:
            self.cmd.unpick()
            UnbondWizard(self.cmd).toggle()

    def cycleBond(self):
        picked = collectPicked(self.cmd)
        if picked == ["pk1","pk2"]:
            with undocontext(self.cmd, "(?pk1 ?pk2) extend 1"):
                self.cmd.cycle_valence()
                self.cmd.unpick()
        else:
            ValenceWizard(_self=self.cmd).toggle(-1,"Cycle bond")

    @undoablemethod("(?pk1 ?pk2) extend 1")
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

    @PopupOnException.decorator
    def invert(self, _=None):
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
            with undocontext(self.cmd, "(?pkset ?pk1) extend 1"):
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
        QMB = QtWidgets.QMessageBox
        check = QMB.question(None, "Confirm",
            "Really delete everything?", QMB.Yes | QMB.No)
        if check == QMB.Yes:
            self.cmd.delete("all")
            self.cmd.refresh_wizard()

    def sculpt(self):
        picked = collectPicked(self.cmd)
        if len(picked):
            self.cmd.select(active_sele, ' or '.join(picked))
        SculptWizard(_self=self.cmd).toggle()

    def clean(self):
        picked = collectPicked(self.cmd)
        if len(picked):
            self.cmd.select(active_sele, "pkmol")
            self.cmd.unpick()
        CleanWizard(_self=self.cmd).toggle()

    def undo(self):
        self.cmd.undo()

    def redo(self):
        self.cmd.redo()

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
            self.cmd.select(active_sele,"byobj ("+' or '.join(picked)+")")
            self.cmd.deselect()
        else:
            self.cmd.delete(active_sele)
        RestAtomWizard(_self=self.cmd).toggle(2)


def BuilderPanelDocked(parent, *args, **kwargs):
    widget = _BuilderPanel(parent, *args, **kwargs)
    window = QtWidgets.QDockWidget(parent)
    window.setWindowTitle("Builder")
    window.setWidget(widget)
    window.setFloating(True)
    return window
