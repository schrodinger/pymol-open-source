
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
        ofile_list = askopenfilename(initialdir = self.fileOpenPath,
                                     filetypes=self.app.getLoadableFileTypes(),
                                     multiple=1) 
        for ofile in ofile_list:
            if len(ofile):
                self.fileOpenPath = os.path.dirname(ofile)
                self.cmd.load(ofile,quiet=0)

    def quit(self):
        self.cmd.log_close()
        self.cmd.quit()

    def createMenuBar(self):
        self.menuBar = Pmw.MenuBar(self.root, # balloon=self.balloon,
                                   hull_relief=RAISED, hull_borderwidth=1)
        self.menuBar.pack(fill=X)        

        self.menuBar.addmenu('File', 'File Menu',tearoff=TRUE)

        self.menuBar.addmenuitem('File', 'command', 'Open file',
                                label='Open...',
                                command=self.fileOpenDialog)

        self.menuBar.addcascademenu('File', 'Skin', 'Skin',
                                             label='Skin',tearoff=TRUE)

        self.app.addSkinMenuItems(self.menuBar,'Skin')
        
        self.menuBar.addmenuitem('File', 'command', 'Quit PyMOL',
                                label='Quit',
                                command=self.quit)

    def destroyMenuBar(self):
        self.menuBar.pack_forget()
        self.menuBar.destroy()
        del self.menuBar
        
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
        self.app._hull.pack(side=TOP, anchor=CENTER, expand=NO, fill=NONE,
                            padx=0, pady=0, ipadx=0, ipady=0)
        
    def takedown(self):
        self.destroyMenuBar()
        PMGSkin.takedown(self)        
        
    def __init__(self,app):

        PMGSkin.__init__(self,app)
        self.app = app
        self.pymol = app.pymol
        self.cmd = app.pymol.cmd
                
def __init__(app):
    return Demo(app)
    
