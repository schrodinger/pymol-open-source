

class PMGSkin:

    def setup(self):
        # if a "look" is active, take it down first before starting new one
        if self.app.skin != self: 
            self.app.skin.takedown()
        self.app.skin = self;
        # (should replace above with accessor calls)
        
    def takedown(self):
        pass

    def __init__(self,app):
        self.app = app
        self.root = app.root
        self.pymol = app.pymol
            
