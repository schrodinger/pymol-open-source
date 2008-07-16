
import sys
import os
from glob import glob
import traceback
import string

from Tkinter import *
import tkMessageBox
import Pmw

from pymol import editor
import pymol

WID = 5
MAX_COL = 12
imgDict = {}
active_sele = "_builder_active"
newest_sele = "_builder_added"

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
    def __init__(self, frame, text, command, comment="", colspan=1):
        b = Button(frame, text=text, padx=1, pady=1, 
            width=WID*colspan, borderwidth=2, relief=RIDGE,
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
        GuiLabel(self, "Atoms")
        GuiButton(self, "H", lambda s=self: s.replace("H",1,1), "Hydrogen")
        GuiButton(self, "C", lambda s=self: s.replace("C",4,4), "Carbon")
        GuiButton(self, "N", lambda s=self: s.replace("N",4,3), "Nitrogen")
        GuiButton(self, "O", lambda s=self: s.replace("O",4,2), "Oxygen")
        GuiButton(self, "P", lambda s=self: s.replace("P",4,3), "Phosphorous")
        GuiButton(self, "S", lambda s=self: s.replace("S",2,2), "Sulfur")
        GuiButton(self, "F", lambda s=self: s.replace("F",1,1), "Fluorine")
        GuiButton(self, "Cl", lambda s=self: s.replace("Cl",1,1), "Chlorine")
        GuiButton(self, "Br", lambda s=self: s.replace("Br",1,1), "Bromine")
        GuiButton(self, "I", lambda s=self: s.replace("I",1,1), "Iodine")

    def replace(self, atom, geometry, valence):
        if getAtoms(self.cmd, 1):
            if "pk1" in self.cmd.get_names("selections"):
                self.cmd.select(active_sele,"byobj pk1")
                self.cmd.replace(atom, geometry, valence)
                self.builder.doAutoPick()
            else:
                warn("Please select (pk1) first...")

class FragmentFrame(GuiFrame):
    def __init__(self, parent):
        self.builder = parent.builder
        self.cmd = self.builder.cmd
        GuiFrame.__init__(self, parent)

        GuiLabel(self, "Fragments")
        GuiButton(self, "CH4", lambda s=self: s.grow("methane",1,0), "Methane")
        GuiButton(self, "C=C", lambda s=self: s.grow("ethylene",4,0))
        GuiButton(self, "C#C", lambda s=self: s.grow("acetylene",2,0), "Acetylene")
        GuiButton(self, "NC=O", lambda s=self: s.grow("formamide",3,1), "N->C amide")
        GuiButton(self, "C=ON", lambda s=self: s.grow("formamide",5,0), "C->N amide")
        GuiButton(self, "C=O", lambda s=self: s.grow("formaldehyde",2,0), "Aldehyde")
        GuiButton(self, "S=O2", lambda s=self: s.grow("sulfone",3,1), "Sulfone")

        GuiLabel(self, "Rings")
        GuiImgButton(self, "cyc4", lambda s=self: s.grow("cyclobutane",4,0), "Cyclobutane")
        GuiImgButton(self, "cyc5", lambda s=self: s.grow("cyclopentane",5,0), "Cyclopentane")
        GuiImgButton(self, "cyc6", lambda s=self: s.grow("cyclohexane",7,0), "Cyclohexane")
        GuiImgButton(self, "cyc7", lambda s=self: s.grow("cycloheptane",8,0), "Cycloheptane")
        #self.nextRow()
        GuiImgButton(self, "aro5", lambda s=self: s.grow("cyclopentadiene",5,0), "Cyclopentadiene")
        GuiImgButton(self, "aro6", lambda s=self: s.grow("benzene",6,0), "Phenyl")

    def grow(self, name, fill, geom):
        if "pk1" in self.cmd.get_names("selections"):
            self.cmd.select(active_sele,"byobj pk1")            
            editor.attach_fragment("pk1", name, fill, geom, _self=self.cmd)
            self.builder.doAutoPick()
        else:
            editor.attach_fragment("", name, fill, geom, _self=self.cmd)
            self.cmd.unpick()
            self.cmd.select(active_sele,name)
            self.builder.doAutoPick()            
        
##############################################################

class ModifyFrame(GuiFrame):
    def __init__(self, parent):
        self.builder = parent.builder
        self.cmd = self.builder.cmd
        GuiFrame.__init__(self, parent)
        GuiLabel(self, "Charge")
        GuiButton(self, "+1", lambda s=self: s.setCharge("1.0"), "Positive Charge")
        GuiButton(self, "0", lambda s=self: s.setCharge("0.0"), "Neutral Charge")
        GuiButton(self, "-1", lambda s=self: s.setCharge("-1.0"), "Negative Charge")

        l = Label(self, text="Bond", width=10)
        l.grid(row=self.row, column=self.col, sticky=E)
        self.nextColumn()

        GuiButton(self, "Create", self.createBond, "Create bond between pk1 and pk2")
        GuiButton(self, "Delete", self.deleteBond, "Delete bond between pk1 and pk2")
        #self.nextRow()
        GuiButton(self, "|", lambda s=self: s.setOrder("1"), "Create single bond")
        GuiButton(self, "||", lambda s=self: s.setOrder("2"), "Create double bond")
        GuiButton(self, "|||", lambda s=self: s.setOrder("3"), "Create triple bond")
#        GuiButton(self, "Cycle", lambda s=self: s.cmd.cycle_valence(1), "Cycle valence (single -> double -> triple")
        GuiButton(self, "Arom", lambda s=self: s.setOrder("4"), "Create aromatic bond")

    def setCharge(self, charge):
        self.cmd.alter("pk1","formal_charge=%s" % charge)
        self.cmd.h_fill()

    def createBond(self):
        if getAtoms(self.cmd, 2):
            if self.cmd.count_atoms("(pk1 or pk2) and (elem h,f,cl,br)") > 0:
                warn("Can not create bond")
            else:
                self.cmd.bond("pk1", "pk2")
                self.cmd.h_fill()
                self.cmd.unpick()

    def deleteBond(self):
        if getAtoms(self.cmd, 2):
            self.cmd.unbond("pk1", "pk2")
            self.cmd.h_fill()
            self.cmd.unpick()

    def setOrder(self, order):
        if getAtoms(self.cmd, 2):
            self.cmd.unbond("pk1", "pk2")
            self.cmd.bond("pk1", "pk2", order)
            self.cmd.h_fill()
            self.cmd.unpick()

class CleanJob:
    def __init__(self,self_cmd,sele):
        self.cmd = self_cmd
        # this code will moved elsewhere
        ok = 1
        try:
            from freemol import mengine
        except:
            ok = 0
            print "Error: Unable to import module freemol.mengine"
        if ok:
            if not mengine.validate():
                ok = 0
                print "Error: Unable to validate freemol.mengine"
        if not ok:
            warn("Please be sure that FreeMOL is correctly installed.")
        else:
            from chempy import io
            obj_list = self_cmd.get_object_list("bymol ("+sele+")")
            ok = 0
            result = None
            if len(obj_list)==1:
                obj_name = obj_list[0]
                state = self_cmd.get_state()
                sdf_list = io.mol.toList(self_cmd.get_model(obj_name,state=state)) + ["$$$$\n"]
                result = mengine.run(string.join(sdf_list,''))
                if result != None:
                    if len(result):
                        clean_sdf = result[0]
                        clean_mol = clean_sdf.split("$$$$")[0]
                        if len(clean_mol):
                            clean_name = "builder_clean_tmp"
                            self_cmd.read_molstr(clean_mol, clean_name, zoom=0)
                            self_cmd.set("retain_order","1",clean_name)
                            self_cmd.fit(clean_name, obj_name, matchmaker=4, quiet=0,
                                         mobile_state=1, target_state=state)
                            self_cmd.update(obj_name, clean_name, matchmaker=0,
                                            source_state=1, target_state=state)
                            self_cmd.delete(clean_name)
                            ok = 1
            if not ok:
                warn("Cleanup failed.  Invalid input or software malfuction?")
                if result != None:
                    if len(result)>1:
                        print result[1]
                            
class EditFrame(GuiFrame):
    def __init__(self, parent):
        self.builder = parent.builder
        self.cmd = self.builder.cmd
        GuiFrame.__init__(self, parent)

        GuiLabel(self, "Editing")
        GuiButton(self, "Fix H", lambda s=self: s.cmd.h_fill(), "Fix hydrogens on pk1")
        GuiButton(self, "Add H", lambda s=self: s.cmd.h_add("pkmol"), "Add hydrogens to molecule")
        GuiButton(self, "Invert", self.invert, "Invert stereochemistry around pk1 (pk2 and pk3 will remain fixed)")
        #self.nextRow()
        GuiButton(self, "Center", self.center, "Center pk1 (or molecule)")
        GuiButton(self, "Delete", self.deleteAtom, "Delete pk1")
        #GuiButton(self, "Reset", self.reset)
        GuiButton(self, "Clear", self.clear, "Delete everything")
        l = Label(self, text="Structure", width=10)
        l.grid(row=self.row, column=self.col, sticky=E)
        self.nextColumn()
        GuiButton(self, "Clean", self.clean, "Cleanup Structure")
        GuiButton(self, "Sculpt", self.sculpt, "Molecular Sculpting")
        GuiButton(self, "Undo", self.undo, "Undo Changes")

    def invert(self):
        if getAtoms(self.cmd, 3):
            self.cmd.invert()

    def center(self):
        if "pk1" in self.cmd.get_names("selections"):
            self.cmd.zoom("pk1", 5.0, animate=-1)
        else:
            self.cmd.zoom("all", 3.0, animate=-1)

    def deleteAtom(self):
        if getAtoms(self.cmd, 1):
            self.cmd.remove_picked()
            self.cmd.unpick()

    def reset(self):
        self.cmd.unpick()

    def clear(self):
        check = tkMessageBox.askokcancel("Confirm", "Really delete everything?")
        if check:
            self.cmd.delete("all")

    def sculpt(self):
        print "sculpt: to come"

    def clean(self):
        if getAtoms(self.cmd, 1):
            CleanJob(self.cmd,"pk1")
        
    def undo(self):
        print "undo: to come"        
        
############################################################

class AminoAcidFrame(GuiFrame):
    def attach(self, aa):
        try:
            editor.attach_amino_acid("pk1", aa, _self=self.cmd)
        except:
            traceback.print_exc()
        self.builder.doZoom()
        
    def __init__(self, parent):
        self.builder = parent.builder
        self.cmd = self.builder.cmd
        GuiFrame.__init__(self, parent)
        #GuiLabel(self, "Amino Acids")

        aaList = ["Ala", "Asp", "Asn", "Arg", "Cys", "Glu", "Gln", "Gly", 
                  "His", "Ile", "Leu", "Lys", "Met", "Phe", "Pro", "Ser", 
                  "Thr", "Trp", "Tyr", "Val", "Ace", "NMe"]
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

    def doValence(self, *ignore):
        self.cmd.set("valence", self.showValence.get())

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
        self.autoPik.set(1)
        autoB = Checkbutton(self, text="autopick", 
            borderwidth=1, pady=0, justify=LEFT, variable=self.autoPik, 
            onvalue=1, offvalue=0, command=self.cmd.unpick)
        autoB.grid(row=4, column=1, sticky=W)

        self.autoZoom = IntVar()
        self.autoZoom.set(0)
        Checkbutton(self, text="autozoom", 
            borderwidth=1, pady=0, justify=LEFT, variable=self.autoZoom, 
            onvalue=1, offvalue=0, command=self.doZoom).grid(row=5, column=1,
            sticky=W)

        self.showValence = StringVar()
        self.showValence.set(self.cmd.get("valence"))
        Checkbutton(self, text="valence", 
            borderwidth=1, pady=0, justify=LEFT, variable=self.showValence, 
            onvalue="on", offvalue="off", command=self.doValence).grid(row=6, 
            column=1, sticky=W)


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

def getAtoms(self_cmd, nAtom):
    """counts how many atoms are selected by (pk) objects
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


