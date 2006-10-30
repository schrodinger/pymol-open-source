
import sys
import os
from glob import glob
import traceback

from Tkinter import *
import tkMessageBox
import Pmw

# NOTE: these all need to be converted over to instance.cmd, etc.
# instead of using global modules

from pymol import cmd
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
        GuiFrame.__init__(self, parent)
        GuiLabel(self, "Atoms")
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
        if getAtoms(1):
            if "pk1" in cmd.get_names("selections"):
                cmd.select(active_sele,"byobj pk1")
                cmd.replace(atom, geometry, valence)
                self.builder.doAutoPick()
            else:
                warn("Please select (pk1) first...")

class FragmentFrame(GuiFrame):
    def __init__(self, parent):
        self.builder = parent.builder
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
        if "pk1" in cmd.get_names("selections"):
            cmd.select(active_sele,"byobj pk1")            
            editor.attach_fragment("pk1", name, fill, geom)
            self.builder.doAutoPick()
        else:
            editor.attach_fragment("", name, fill, geom)
            cmd.unpick()
            cmd.select(active_sele,name)
            self.builder.doAutoPick()            
        
##############################################################

class ModifyFrame(GuiFrame):
    def __init__(self, parent):
        self.builder = parent.builder        
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
        GuiButton(self, "Cycle", lambda: cmd.cycle_valence(1), "Cycle valence (single -> double -> triple")

    def setCharge(self, charge):
        cmd.alter("pk1","formal_charge=%s" % charge)
        cmd.h_fill()

    def createBond(self):
        if getAtoms(2):
            if cmd.count_atoms("(pk1 or pk2) and (elem h,f,cl,br)") > 0:
                warn("Can not create bond")
            else:
                cmd.bond("pk1", "pk2")
                cmd.h_fill()
                cmd.unpick()

    def deleteBond(self):
        if getAtoms(2):
            cmd.unbond("pk1", "pk2")
            cmd.h_fill()
            cmd.unpick()

    def setOrder(self, order):
        if getAtoms(2):
            cmd.unbond("pk1", "pk2")
            cmd.bond("pk1", "pk2", order)
            cmd.h_fill()
            cmd.unpick()

class EditFrame(GuiFrame):
    def __init__(self, parent):
        self.builder = parent.builder
        GuiFrame.__init__(self, parent)

        GuiLabel(self, "Editing")
        GuiButton(self, "Fix H", lambda: cmd.h_fill(), "Fix hydrogens on pk1")
        GuiButton(self, "Add H", lambda: cmd.h_add("pkmol"), "Add hydrogens to molecule")
        GuiButton(self, "Invert", self.invert, "Invert stereochemistry around pk1 (pk2 and pk3 will remain fixed)")
        #self.nextRow()
        GuiButton(self, "Center", self.center, "Center pk1 (or molecule)")
        GuiButton(self, "Delete", self.deleteAtom, "Delete pk1")
        #GuiButton(self, "Reset", self.reset)
        GuiButton(self, "Clear", self.clear, "Delete everything")

    def invert(self):
        if getAtoms(3):
            cmd.invert()

    def center(self):
        if "pk1" in cmd.get_names("selections"):
            cmd.zoom("pk1", 5.0)
        else:
            cmd.zoom("all", 3.0)

    def deleteAtom(self):
        if getAtoms(1):
            cmd.remove_picked()
            cmd.unpick()

    def reset(self):
        cmd.unpick()

    def clear(self):
        check = tkMessageBox.askokcancel("Confirm", "Really delete everything?")
        if check:
            cmd.delete("all")


############################################################

class AminoAcidFrame(GuiFrame):
    def attach(self, aa):
        try:
            editor.attach_amino_acid("pk1", aa)
        except:
            traceback.print_exc()
        self.builder.doZoom()
        
    def __init__(self, parent):
        self.builder = parent.builder        
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
        GuiFrame.__init__(self, parent)
        self.secStructL = ["alpha helix", "parallel beta sheet",
                           "antiparallel beta sheet"]
        self.ss = int(float(cmd.get("secondary_structure")))-1
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
        cmd.set("secondary_structure", "%d" % (self.ss+1))


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
    def __init__(self, parent, *kw, **args):
        """
        module for constructing tkinter gui & widgets for molecule
        building and editing with PyMol

        Builder constructs the top level gui that will contain
        other widgets for working with chemicals, proteins, etc.
        """
        self.builder = self

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
            cmd.unpick()
            if cmd.select(newest_sele,"(byobj "+active_sele+") and not "+active_sele)==0:
                cmd.select(newest_sele, active_sele)
            new_list = cmd.index(newest_sele+" and hydro")
            if len(new_list)==0:
                new_list = cmd.index(newest_sele)                
            if new_list:
                index = new_list.pop()
                try:
                    cmd.edit("%s`%d" % index)
                    if cmd.get_wizard()!=None:
                        cmd.do("_ cmd.get_wizard().do_pick(0)")
                except pymol.CmdException:
                    print " doAutoPick-Error: exception"
            self.doZoom()

    def doZoom(self, *ignore):
#        print "zoom",self.autoZoom.get()
        if self.autoZoom.get():
            if "pk1" in cmd.get_names("selections"):
                cmd.zoom("((neighbor pk1) expand 5)", 5.0)

    def doValence(self, *ignore):
        cmd.set("valence", self.showValence.get())

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
            onvalue=1, offvalue=0, command=cmd.unpick)
        autoB.grid(row=4, column=1, sticky=W)

        self.autoZoom = IntVar()
        self.autoZoom.set(0)
        Checkbutton(self, text="autozoom", 
            borderwidth=1, pady=0, justify=LEFT, variable=self.autoZoom, 
            onvalue=1, offvalue=0, command=self.doZoom).grid(row=5, column=1,
            sticky=W)

        self.showValence = StringVar()
        self.showValence.set(cmd.get("valence"))
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

def getAtoms(nAtom):
    """counts how many atoms are selected by (pk) objects
    @param nAtom: number of atoms to select
    @type nAtom: integer
    @returns: 1 if number selected == nAtoms, 0 if not
    @rtype: integer
    """
    count = 0
    for i in range(1,5):
        if "pk%d" % i in cmd.get_names("selections"):
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


