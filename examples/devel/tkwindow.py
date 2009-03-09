# demonstration plugin which opens up a window on startup
#
# you must copy this file into modules/pmg_tk/startup/
#
# pymol will then run the script on startup

from Tkinter import *
from tkFileDialog import *
from pymol import cmd

class myTkWindow:

    def __init__(self,app):
        self.root = Tk()
        self.app = app
        self.root.geometry(newGeometry="+900+50")
        Button(self.root,
               text="Open File",
               command=lambda s=self:s.open()).pack()

    def open(self):
        file_list = askopenfilename(multiple=1)
        for file in file_list:
            self.app.pymol.cmd.load(file)
        
def __init__(app):
    myTkWindow(app)
    
        
