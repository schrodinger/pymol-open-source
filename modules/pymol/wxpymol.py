
from wxPython.wx       import *
from wxPython.glcanvas import *
from OpenGL.GL import *
from OpenGL.GLUT import *
import _cmd
import threading
import pymol

if 1:

    class MyCanvasBase(wxGLCanvas):
        def __init__(self, parent):
            wxGLCanvas.__init__(self, parent, -1)
            self.init = false
            # initial mouse position
            self.lastx = self.x = 30
            self.lasty = self.y = 30
            EVT_ERASE_BACKGROUND(self, self.OnEraseBackground)
            EVT_SIZE(self, self.OnSize)
            EVT_PAINT(self, self.OnPaint)
            EVT_LEFT_DOWN(self, self.OnMouseDown)  # needs fixing...
            EVT_LEFT_UP(self, self.OnMouseUp)
            EVT_MOTION(self, self.OnMouseMotion)
            EVT_IDLE(self,self.OnIdle)
            
        def OnEraseBackground(self, event):
            pass # Do nothing, to avoid flashing on MSW.

        def OnSize(self, event):
            size = self.GetClientSize()
            if self.GetContext():
                self.SetCurrent()
                glViewport(0, 0, size.width, size.height)
                _cmd.p_glut_event(2,size.width,size.height,0,0,0)
                
        def OnPaint(self, event):
            dc = wxPaintDC(self)
            self.SetCurrent()
            if not self.init:
                self.InitGL()
                self.init = true
            self.OnDraw()

        def OnMouseDown(self, evt):
            self.CaptureMouse()
            x,y = evt.GetPosition()
            if evt.LeftIsDown():
                _cmd.p_glut_event(3,x,y,15,9,0)
                self.Repaint()
                self.mouse='l'
            elif evt.MiddleIsDown():
                _cmd.p_glut_event(3,x,y,16,9,0)
                self.mouse='m'
                self.Repaint()                
            elif evt.RightIsDown():
                _cmd.p_glut_event(3,x,y,17,9,0)
                self.mouse='r'
                self.Repaint()
            self.CheckPyMOL()                
#            print "mouse down"
            
        def OnMouseUp(self, evt):
            x,y = evt.GetPosition()
            if self.mouse=='l':
                _cmd.p_glut_event(3,x,y,15,10,0)
                self.Repaint()
            elif self.mouse=='m':
                _cmd.p_glut_event(3,x,y,16,10,0)
                self.Repaint()
            elif self.mouse=='r':
                _cmd.p_glut_event(3,x,y,17,10,0)
                self.Repaint()                
            self.CheckPyMOL()                
#            print "mouse up"
            self.ReleaseMouse()
            
        def OnMouseMotion(self, evt):
            if evt.Dragging():
                x,y = evt.GetPosition()
                if evt.LeftIsDown():
                    _cmd.p_glut_event(4,x,y,15,0,0)
                    self.Repaint()
                elif evt.MiddleIsDown():
                    _cmd.p_glut_event(4,x,y,16,0,0)                
                    self.Repaint()
                elif evt.RightIsDown():
                    _cmd.p_glut_event(4,x,y,17,0,0)
                    self.Repaint()
                self.CheckPyMOL()                
#                print "mouse motion"
            
        def OnIdle(self,evt):
            _cmd.p_glut_event(0,0,0,0,0,0) # idle event
#            self.AddPendingEvent(wxIdleEvent())
#            evt.RequestMore(true)
            self.CheckPyMOL()

        def Repaint(self):
            self.AddPendingEvent(wxPaintEvent())
            
    class CubeCanvas(MyCanvasBase):
        def InitGL(self):

            _cmd.runwxpymol()
            pymol._swap_buffers = lambda s=self: s.swap()
            _cmd.p_glut_event(2,640,480,0,0,0) # initial reshape event
            self.CheckPyMOL()

        def swap(self):
            if self.GetContext():
                self.SetCurrent()
#                print "swapping"
                self.SwapBuffers()
            
        def OnDraw(self):
            _cmd.p_glut_event(1,0,0,0,0,0) # draw event

        def CheckPyMOL(self):
            if _cmd.p_glut_get_redisplay():
                self.Repaint()

#----------------------------------------------------------------------

overview = """\
"""
#----------------------------------------------------------------------

class MyApp(wxApp):
  def OnInit(self):
      frame = wxFrame(None, -1, "CubeCanvas", size=(640,480))
      canvas = CubeCanvas(frame)
      frame.Show(true)
      self.SetTopWindow(frame)

      import threading
      return TRUE

app = MyApp(0)
app.MainLoop()

