
# Demonstration of how you would create a new skin
# from scratch

from Tkinter import *
from tkFileDialog import *

from pmg_tk.skins import PMGSkin

class Demo(PMGSkin):
    
    def setup(self):

        # call the parent method
        PMGSkin.setup(self)
        
        # name the application
        self.root.title("Demonstration PyMOL Skin")

        # pack the root window
        self.app._hull.pack(side=LEFT, fill=BOTH, expand=YES)

    def __init__(self,app):
        PMGSkin.__init__(self,app)
        self.app = app
        
def __init__(app):
    app.set_skin(Demo(app))

    
