

class PMGSkin:

    def setup(self):
        # start with a clean packer
        self.app._hull.pack_forget()
    
    def takedown(self):
        pass
    
    def __init__(self,app):
        self.app = app
        self.root = app.root
        self.pymol = app.pymol
            
