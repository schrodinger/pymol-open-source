#A* -------------------------------------------------------------------
#B* This file contains source code for the PyMOL computer program
#C* copyright 1998-2000 by Warren Lyford Delano of DeLano Scientific. 
#D* -------------------------------------------------------------------
#E* It is unlawful to modify or remove this copyright notice.
#F* -------------------------------------------------------------------
#G* Please see the accompanying LICENSE file for further information. 
#H* -------------------------------------------------------------------
#I* Additional authors of this source file include:
#-*
#-* NOTE: Based on code by John E. Grayson which was in turn 
#-* based on code written by Doug Hellmann. 
#Z* -------------------------------------------------------------------

import sys
import os
from glob import glob
import re
import traceback
import Pmw

if True:
    import queue as Queue
    from tkinter import *
    from tkinter.filedialog import *
    import tkinter.messagebox as tkMessageBox

try:
    # monkey patch Pmw's error message box
    def _reporterror(func, args):
        import pymol

        exc_type, exc_value, exc_traceback = sys.exc_info()
        msg = "Sorry this command was not successful at this time"
        if issubclass(exc_type, (pymol.cmd.QuietException)):
            tkMessageBox.showerror('Error', msg)
        elif issubclass(exc_type, (pymol.CmdException)):
            tkMessageBox.showerror(getattr(exc_value, 'label', 'Error'), str(exc_value.message).strip() or msg)
        else:
            _reporterror.orig(func, args)

    PmwBase = sys.modules[Pmw.MegaWidget.__module__]
    _reporterror.orig = PmwBase._reporterror
    PmwBase._reporterror = _reporterror

    del PmwBase
except Exception as e:
    print('monkey patching _reporterror failed:', e)

class PMGApp(Pmw.MegaWidget):

    def initOS(self):
         # Initialize platform-specific options
         if sys.platform == 'darwin':
             self.initializeTk_mac()
         elif sys.platform[:3] == 'win':
             self.initializeTk_win32()
         elif sys.platform[:5] == 'linux':
             self.initializeTk_unix()
         else:
             self.initializeTk_unix()
#       self.root.tk.call('tk','scaling',1)

         # try to get the windows properly aligned...
         osFrame = { 'win32' : (4,60), 'irix'   : (0,41),
                     'darwin': (0,51), 'cygwin' : (0,60),
                     'linux' : (0,31), 'linux2' : (0,31) }

         self.frameXAdjust, self.frameYAdjust = osFrame.get(sys.platform, (0, 51))
         
    def initializeTk_win32(self):
        pass
        
    def initializeTk_mac(self):
        pass
        
    def initializeTk_unix(self):
        pass

    def getLoadableFileTypes(self):
        return [("All Readable","*.pdb"),
                ("All Readable","*.pdb1"),                
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
                ("All Readable","*.psw"),
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
                ("All Readable","*.maegz"), # proprietary format
                ("All Readable","*.cms"), # proprietary format
                ("All Readable","*.idx"), # proprietary format
                ("All Readable","*.fasta"),
                ("All Readable","*.aln"),
                ("All Readable","*.acnt"),
                ("All Readable","*.mtz"),                 
                ("All Readable","*.vis"),
                ("All Readable","*.psf"),
                ("All Readable","*.pdbml"),
                ("All Readable","*.xml"),
                ("All Readable","*.xml.gz"),
                ("All Readable","*.pdbqt"),
                ("All Readable","*.cml"),
                ("All Readable","*.mmtf"),
                ("PDB File","*.pdb"),
                ("PDB1 File","*.pdb1"),                
                ("All Files","*.*"),
                ("All Files","*"),  
                ("PDB File","*.ent"),
                ("PyMOL Session","*.pse"),
                ("PyMOL Show","*.psw"),
                ("CCP4 Map","*.ccp4"), 
                ("XPLOR Map","*.xplor"),
                ("MOL2/Multi-MOL2","*.mol2"),
                ("Macromodel File","*.dat"),
                ("Macromodel File","*.out"),
                ("Macromodel File","*.mmd"),
                ("Macromodel File","*.mmod"),
#                ("MTZ Reflection File","*.mtz"),
                ("BRIX/O Map","*.o"),
                ("BRIX/O Map","*.omap"),
                ("BRIX/O Map","*.brix"),
                ("CIF","*.cif"),
                ("Gaussian Cube Map","*.cube"),
                ("DX Map","*.dx"),     
                ("AVS (MEAD) Field","*.fld"),    
                ("MOL File","*.mol"),
                ("MOE File","*.moe"), # proprietary format
                ("MAE File","*.mae"), # proprietary format  
                ("MAE File","*.maegz"), # proprietary format
                ("MAE File","*.cms"), # proprietary format
                ("Desmond Trajectory","*.idx"), # proprietary format
                ("VIS File","*.vis"),
                ("ChemPy Model","*.pkl"),
                ("Raster3D Scene","*.r3d"),
                ("SDF File","*.sdf"),
                ("ChemDraw3D File","*.cc1"),
                ("ChemDraw3D File","*.cc2"),
                ("XYZ File","*.xyz"),
                ("Fasta File","*.fasta"),
                ("CLUSTAL file","*.aln"),
                ("ACNT Map","*.acnt"),                 
                ("Protein Structure File","*.psf"),
                ("PDBML","*.pdbml"),
                ("PDBML","*.xml"),
                ("PDBML","*.xml.gz"),
                ("PDBQT","*.pdbqt"),
                ("Chemical Markup Language","*.cml"),
                ("MMTF","*.mmtf"),
                ("MMTF","*.mmtf.gz"),
                ]

    def initializeTk_colors_common(self):
        #self.root.option_add('*background', 'grey')   #let system decide
        self.root.option_add('*foreground', 'black')
        self.root.option_add('*EntryField.Entry.background', 'white')
        self.root.option_add('*Entry.background', 'white')      
        self.root.option_add('*MessageBar.Entry.background', 'gray85')
        self.root.option_add('*Listbox*background', 'white')
        self.root.option_add('*Listbox*selectBackground', 'dark slate blue')
        self.root.option_add('*Listbox*selectForeground', 'white')
        
    def quit_app(self):
        self.pymol.cmd.log_close()
        self.pymol.cmd.quit()  # avoid logging this - it's inconvenient...

    def flush_fifo_once(self):
        # flush the external GUI fifo command queue
        while not self.fifo.empty():
            try:
                cmmd = self.fifo.get(0)
                if isinstance(cmmd, str):
                    exec(cmmd)
                else:
                    cmmd()
            except:
                traceback.print_exc()
        
    def flush_fifo(self):
        self.flush_fifo_once()
        if self.allow_after:
            self.root.after(20,self.flush_fifo) # 50X a second
        
    def run(self,poll=0):
        # this call to mainloop needs to be replaced with something revocable
        self.flush_fifo_once()
        keep_alive = 1
        if poll:
            import time
            while keep_alive:
                self.root.update()
                time.sleep(0.05)
        else:
            self.root.mainloop()
            
        self.quit_app()

    def execute(self,cmmd):  
        self.fifo.put(cmmd)

    def my_show(self,win,center=1):
        if sys.platform!='linux2':
            win.show()
        else: # autocenter, deiconify, and run mainloop
            # this is a workaround for a bug in the
            # interaction between Tcl/Tk and common Linux
            # window managers (namely KDE/Gnome) which causes
            # an annoying 1-2 second delay in opening windows!
            if center:
                tw = win.winfo_reqwidth()+100
                th = win.winfo_reqheight()+100
                vw = win.winfo_vrootwidth()
                vh = win.winfo_vrootheight()
                x = max(0,(vw-tw)/2)
                y = max(0,(vh-th)/2)
                win.geometry(newGeometry="+%d+%d"%(x,y))
            win.deiconify()

    def my_withdraw(self,win):
        if sys.platform!='linux2':
            win.withdraw()
        else: 
            win.destroy()

    def _initializePlugins(self):
        from pymol.plugins import legacysupport
        return legacysupport.initializePlugins(self)

    def addSkinMenuItems(self,menuBar,menuName):
        if not hasattr(self,'skinNameList'):
            # find installed skins
            skin_pattern = re.sub(r"[\/\\][^\/\\]*$","/skins/*/__init__.py*",__file__)
            raw_list = glob(skin_pattern)
            unique = {}
            for a in raw_list:
                key = re.sub(r"[\/\\]__init__\.py.*$","",a)
                unique[re.sub(r".*[\/\\]","",key)] = 1
            name_list = list(unique.keys())
            name_list.sort()
            self.skinNameList = name_list
        if hasattr(self,'skinNameList'):
            for name in self.skinNameList:
                caps_name = name.capitalize()
                menuBar.addmenuitem(menuName, 'command', name,
                                    label=caps_name,
                                    command=lambda s=self,n=name:s.setSkin(n))

        
    def setSkin(self,skin,run=1):
        if isinstance(skin,str):
            inv = sys.modules.get("pymol.invocation",None)
            if inv!=None:
                module_path = inv.options.gui +".skins."+ skin
                __import__(inv.options.gui +".skins."+ skin)
                skin = sys.modules[module_path].__init__(self)
        if skin != self.skin:
            if self.skin != None:
                self.skin.takedown()
            self.skin = skin

        if run:
            self.runSkin()
            
    def runSkin(self):
        if self.skin != None:
            self.skin.setup()
    
    def __init__(self, pymol_instance, skin):

        # prevent overloading
        self.initializePlugins = self._initializePlugins

        self.allow_after = 1 # easy switch for troubleshooting threads

        self.pymol = pymol_instance
        
        if self.pymol._ext_gui != None:
        
            raise RuntimeError  # only one PMGApp should ever be instantiated
        
        else:

            # create a FIFO so that PyMOL can send code to be executed by the GUI thread
            
            self.fifo = Queue.Queue(0)

            # create a pymol global so that PyMOL can find the external GUI

            self.pymol._ext_gui = self

            self.skin = None
            
            # initialize Tcl/Tk

            self.root = Tk() # creates the root window for the application

            # color scheme

            self.initializeTk_colors_common()

             # operating-system dependencies

            self.initOS()

            # Python megawigit initialization

            Pmw.initialise(self.root)
            
            # Initialize the base class
            
            Pmw.MegaWidget.__init__(self, parent=self.root)

            # read the command line arguments regarding:
            # - the size of the root window
            # - the skin to use

            inv = sys.modules.get("pymol.invocation",None)
            if inv != None:
                if skin == None:
                    skin = inv.options.skin
                self.frameWidth = inv.options.win_x + 220
                self.frameXPos = inv.options.win_px - self.frameXAdjust
                self.frameHeight = inv.options.ext_y
                self.frameYPos = inv.options.win_py - (
                         self.frameHeight + self.frameYAdjust)
                self.setSkin(skin,run=0)
                
            # define the size of the root window
            
            import platform
            if sys.platform == 'darwin' and platform.mac_ver()[0] >= '10.9':
                # let OS X Maverics place the window automatically, to avoid
                # off-screen placement in multi-monitor setup
                self.root.geometry('%dx%d' % (self.frameWidth, self.frameHeight))
            else:
                self.root.geometry('%dx%d+%d+%d' % (
                self.frameWidth, self.frameHeight, self.frameXPos, self.frameYPos))
            
            # activate polling on the fifo

            if self.allow_after:
                self.root.after(1000,self.flush_fifo)

            # and let 'er rip

            self.runSkin()

