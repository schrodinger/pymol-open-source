
from wxPython.wx       import *
from wxPython.glcanvas import *
from OpenGL.GL import *
from OpenGL.GLUT import *
import _cmd
import threading
import pymol

P_GLUT_IDLE_EVENT          =  0
P_GLUT_DISPLAY_EVENT       =  1
P_GLUT_RESHAPE_EVENT       =  2
P_GLUT_MOUSE_EVENT         =  3
P_GLUT_MOTION_EVENT        =  4

P_GLUT_ACTIVE_ALT             =  32
P_GLUT_ACTIVE_CTRL            =  64
P_GLUT_ACTIVE_SHIFT           =  128

P_GLUT_DOWN                   =  9
P_GLUT_UP                     =  10
P_GLUT_KEY_DOWN               =  11
P_GLUT_KEY_LEFT               =  12
P_GLUT_KEY_RIGHT              =  13
P_GLUT_KEY_UP                 =  14
P_GLUT_LEFT_BUTTON            =  15
P_GLUT_MIDDLE_BUTTON          =  16
P_GLUT_RIGHT_BUTTON           =  17

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
            EVT_LEFT_DOWN(self, self.OnMouseDown)  
            EVT_LEFT_UP(self, self.OnMouseUp)
            EVT_MIDDLE_DOWN(self, self.OnMouseDown)
            EVT_MIDDLE_UP(self, self.OnMouseUp)
            EVT_RIGHT_DOWN(self, self.OnMouseDown) 
            EVT_RIGHT_UP(self, self.OnMouseUp)
            EVT_MOTION(self, self.OnMouseMotion)
            EVT_IDLE(self,self.OnIdle)
            self.mouse = ''
            self.mod = 0
            
        def OnEraseBackground(self, event):
            pass # Do nothing, to avoid flashing on MSW.

        def OnSize(self, event):
            size = self.GetClientSize()
            if self.GetContext():
                self.SetCurrent()
                glViewport(0, 0, size.width, size.height)
                _cmd.p_glut_event(P_GLUT_RESHAPE_EVENT,size.width,size.height,0,0,0)
                
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
            self.mod = 0
            if evt.ShiftDown():
                self.mod = self.mod + P_GLUT_ACTIVE_SHIFT
            if evt.ControlDown():
                self.mod = self.mod + P_GLUT_ACTIVE_CTRL
            if evt.LeftIsDown():
                _cmd.p_glut_event(P_GLUT_MOUSE_EVENT,x,y,P_GLUT_LEFT_BUTTON,P_GLUT_DOWN,self.mod)
                self.Repaint()
                self.mouse='l'
            elif evt.MiddleIsDown():
                _cmd.p_glut_event(P_GLUT_MOUSE_EVENT,x,y,P_GLUT_MIDDLE_BUTTON,P_GLUT_DOWN,self.mod)
                self.mouse='m'
                self.Repaint()                
            elif evt.RightIsDown():
                _cmd.p_glut_event(P_GLUT_MOUSE_EVENT,x,y,P_GLUT_RIGHT_BUTTON,P_GLUT_DOWN,self.mod)
                self.mouse='r'
                self.Repaint()
            self.CheckPyMOL()                
#            print "mouse down"
            
        def OnMouseUp(self, evt):
            x,y = evt.GetPosition()
            if self.mouse=='l':
                _cmd.p_glut_event(P_GLUT_MOUSE_EVENT,x,y,P_GLUT_LEFT_BUTTON,P_GLUT_UP,self.mod)
                self.Repaint()
            elif self.mouse=='m':
                _cmd.p_glut_event(P_GLUT_MOUSE_EVENT,x,y,P_GLUT_MIDDLE_BUTTON,P_GLUT_UP,self.mod)
                self.Repaint()
            elif self.mouse=='r':
                _cmd.p_glut_event(P_GLUT_MOUSE_EVENT,x,y,P_GLUT_RIGHT_BUTTON,P_GLUT_UP,self.mod)
                self.Repaint()                
            self.mouse = ''
            self.CheckPyMOL()
#            print "mouse up"
            self.ReleaseMouse()
            
        def OnMouseMotion(self, evt):
            if evt.Dragging():
                x,y = evt.GetPosition()
                if evt.LeftIsDown():
                    _cmd.p_glut_event(P_GLUT_MOTION_EVENT,x,y,P_GLUT_LEFT_BUTTON,0,self.mod)
                    self.Repaint()
                elif evt.MiddleIsDown():
                    _cmd.p_glut_event(P_GLUT_MOTION_EVENT,x,y,P_GLUT_MIDDLE_BUTTON,0,self.mod)                
                    self.Repaint()
                elif evt.RightIsDown():
                    _cmd.p_glut_event(P_GLUT_MOTION_EVENT,x,y,P_GLUT_RIGHT_BUTTON,0,self.mod)
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
            _cmd.p_glut_event(P_GLUT_RESHAPE_EVENT,800,480,0,0,0) # initial reshape event
            self.CheckPyMOL()

        def swap(self):
            if self.GetContext():
                self.SetCurrent()
#                print "swapping"
                self.SwapBuffers()
            
        def OnDraw(self):
            _cmd.p_glut_event(P_GLUT_DISPLAY_EVENT,0,0,0,0,0) # draw event

        def CheckPyMOL(self):
            if _cmd.p_glut_get_redisplay():
                self.Repaint()

#----------------------------------------------------------------------

overview = """\
"""
#----------------------------------------------------------------------

class MyApp(wxApp):
  def OnInit(self):
      frame = wxFrame(None, -1, "PyMOL", size=(800,480))
      canvas = CubeCanvas(frame)
      frame.Show(true)
      self.SetTopWindow(frame)

      import threading
      return TRUE

app = MyApp(0)
app.MainLoop()

