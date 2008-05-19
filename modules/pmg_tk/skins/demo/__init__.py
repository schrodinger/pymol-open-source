
# Demonstration of how you would create a new skin
# from scratch

import os

from Tkinter import *
from tkFileDialog import *

import Pmw

from pmg_tk.skins import PMGSkin

class Demo(PMGSkin):

    def fileOpenDialog(self):
        if not hasattr(self,'fileOpenPath'):
            self.fileOpenPath = os.getcwd()
        ftypes =  [("All Readable","*.pdb"),
                   ("All Readable","*.ccp4"),
                   ("All Readable","*.xplor"),
                   ("All Readable","*.mol"),
                   ("All Readable","*.mol2"),
                   ("All Readable","*.sdf"),
                   ("All Readable","*.xyz"),                                         
                   ("All Readable","*.r3d"),
                   ("All Readable","*.cc1"),
                   ("All Readable","*.cc2"),                                         
                   ("All Readable","*.ent"),
                   ("All Readable","*.dat"),
                   ("All Readable","*.out"),
                   ("All Readable","*.mmd"),
                   ("All Readable","*.mmod"),
                   ("All Readable","*.pse"),
                   ("All Readable","*.phi"),
                   ("All Readable","*.fld"),
                   ("All Readable","*.grd"),
                   ("All Readable","*.o"),
                   ("All Readable","*.omap"),                                         
                   ("All Readable","*.brix"),
                   ("All Readable","*.dx"),
                   ("All Readable","*.pqr"),
                   ("All Readable","*.p5m"),
                   ("All Readable","*.p1m"),
                   ("All Readable","*.cube"),
                   ("All Readable","*.cif"),
                   ("All Readable","*.moe"), # proprietary format
                   ("All Readable","*.mae"), # proprietary format
                   ("PDB File","*.pdb"),
                   ("All Files","*.*"),
                   ("All Files","*"),                                         
                   ("PDB File","*.ent"),
                   ("PyMOL Session","*.pse"),
                   ("CCP4 Map","*.ccp4"),                                         
                   ("XPLOR Map","*.xplor"),
                   ("MOL2/Multi-MOL2","*.mol2"),
                   ("Macromodel File","*.dat"),
                   ("Macromodel File","*.out"),
                   ("Macromodel File","*.mmd"),
                   ("Macromodel File","*.mmod"),
                   ("BRIX/O Map","*.o"),
                   ("BRIX/O Map","*.omap"),
                   ("BRIX/O Map","*.brix"),
                   ("CIF","*.cif"),
                   ("Gaussian Cube Map","*.cube"),
                   ("DX Map","*.dx"),                                         
                   ("AVS (MEAD) Field","*.fld"),                                         
                   ("MOL File","*.mol"),
                   ("MOE File","*.moe"), # proprietary format
                   ("MOE File","*.mae"), # proprietary format                       
                   ("ChemPy Model","*.pkl"),
                   ("Raster3D Scene","*.r3d"),
                   ("SDF File","*.sdf"),
                   ("ChemDraw3D File","*.cc1"),
                   ("ChemDraw3D File","*.cc2"),
                   ("Tinker XYZ File","*.xyz")
                   ]
        ofile_list = askopenfilename(initialdir = self.fileOpenPath,
                                     filetypes=ftypes,
                                     multiple=1) 
        for ofile in ofile_list:
            if len(ofile):
                self.fileOpenPath = os.path.dirname(ofile)
                self.cmd.load(ofile,quiet=0)

    def createMenuBar(self):
        self.menuBar = Pmw.MenuBar(self.root, # balloon=self.balloon,
                                   hull_relief=RAISED, hull_borderwidth=1)
        self.menuBar.pack(fill=X)

        self.menuBar.addmenu('File', 'File Menu',tearoff=TRUE)

        self.menuBar.addmenuitem('File', 'command', 'Open file.',
                                label='Open...',
                                command=self.fileOpenDialog)

    def createInterface(self):

        # create the menu bar
        self.createMenuBar()
        
    def setup(self):

        # call the parent method
        PMGSkin.setup(self)
        
        # name the application
        self.root.title("Demonstration PyMOL Skin")

        # create the user interface
        self.createInterface()

        # pack the root window
        self.app._hull.pack(side=LEFT, fill=BOTH, expand=YES)

        
    def __init__(self,app):
        PMGSkin.__init__(self,app)
        self.app = app
        self.pymol = app.pymol
        self.cmd = app.pymol.cmd
                
def __init__(app):
    app.set_skin(Demo(app))

    
