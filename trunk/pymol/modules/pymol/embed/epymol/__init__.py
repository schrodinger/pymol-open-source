
from pymol.embed import EmbeddedPyMOL

class ePyMOL(EmbeddedPyMOL):
    def __init__(self):
        self.ep_init()
        # initial mouse position
        self.lastx = self.x = 30
        self.lasty = self.y = 30
        
    def SetSize(self, width, height):
        self.ep_reshape(width,height)

    def OnChar(self,code):
        self.ep_char(0,0,code,0,0,0)

    def OnSpecial(self,code):
        self.ep_special(0,0,code,0,0,0)
        
    def OnPaint(self):
        self.OnDraw()

    def OnMouseDown(self,*arg):
        self.ep_mouse_down(*arg)

    def OnMouseUp(self,*arg):
        self.ep_mouse_up(*arg[0:2])
        
    def OnMouseMotion(self,*arg):
        self.ep_motion(*arg)                  

    def OnDraw(self):
        self.ep_draw()
        
    def OnIdle(self):
        self.ep_idle()

    def GetRedisplay(self):
        return self.ep_get_redisplay()

    
    def CheckPyMOL(self):
        pass
        #if self.ep_get_redisplay():
        #   self.Repaint()
            
    def Repaint(self):
        pass



