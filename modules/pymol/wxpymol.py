
from wxPython.wx       import *
from wxPython.glcanvas import *
from OpenGL.GL import *
from OpenGL.GLUT import *
import _cmd

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
            #EVT_LEFT_DOWN(self, self.OnMouseDown)  # needs fixing...
            #EVT_LEFT_UP(self, self.OnMouseUp)
            #EVT_MOTION(self, self.OnMouseMotion)

        def OnEraseBackground(self, event):
            pass # Do nothing, to avoid flashing on MSW.

        def OnSize(self, event):
            size = self.GetClientSize()
            if self.GetContext():
                self.SetCurrent()
                glViewport(0, 0, size.width, size.height)

        def OnPaint(self, event):
            dc = wxPaintDC(self)
            self.SetCurrent()
            if not self.init:
                self.InitGL()
                self.init = true
            self.OnDraw()

        def OnMouseDown(self, evt):
            self.CaptureMouse()

        def OnMouseUp(self, evt):
            self.ReleaseMouse()

        def OnMouseMotion(self, evt):
            if evt.Dragging() and evt.LeftIsDown():
                self.x, self.y = self.lastx, self.lasty
                self.x, self.y = evt.GetPosition()
                self.Refresh()

    class CubeCanvas(MyCanvasBase):
        def InitGL(self):
            # set viewing projection
            glMatrixMode(GL_PROJECTION);
            glFrustum(-0.5, 0.5, -0.5, 0.5, 1.0, 3.0);

            # position viewer
            glMatrixMode(GL_MODELVIEW);
            glTranslatef(0.0, 0.0, -2.0);

            # position object
            glRotatef(self.y, 1.0, 0.0, 0.0);
            glRotatef(self.x, 0.0, 1.0, 0.0);

            glEnable(GL_DEPTH_TEST);
            glEnable(GL_LIGHTING);
            glEnable(GL_LIGHT0);
            _cmd.runwxpymol()
            _cmd.p_glut_event(2)
            
            
        def OnDraw(self):

            # clear color and depth buffers
            glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

            # draw six faces of a cube
            glBegin(GL_QUADS)
            glNormal3f( 0.0, 0.0, 1.0)
            glVertex3f( 0.5, 0.5, 0.5)
            glVertex3f(-0.5, 0.5, 0.5)
            glVertex3f(-0.5,-0.5, 0.5)
            glVertex3f( 0.5,-0.5, 0.5)

            glNormal3f( 0.0, 0.0,-1.0)
            glVertex3f(-0.5,-0.5,-0.5)
            glVertex3f(-0.5, 0.5,-0.5)
            glVertex3f( 0.5, 0.5,-0.5)
            glVertex3f( 0.5,-0.5,-0.5)

            glNormal3f( 0.0, 1.0, 0.0)
            glVertex3f( 0.5, 0.5, 0.5)
            glVertex3f( 0.5, 0.5,-0.5)
            glVertex3f(-0.5, 0.5,-0.5)
            glVertex3f(-0.5, 0.5, 0.5)

            glNormal3f( 0.0,-1.0, 0.0)
            glVertex3f(-0.5,-0.5,-0.5)
            glVertex3f( 0.5,-0.5,-0.5)
            glVertex3f( 0.5,-0.5, 0.5)
            glVertex3f(-0.5,-0.5, 0.5)

            glNormal3f( 1.0, 0.0, 0.0)
            glVertex3f( 0.5, 0.5, 0.5)
            glVertex3f( 0.5,-0.5, 0.5)
            glVertex3f( 0.5,-0.5,-0.5)
            glVertex3f( 0.5, 0.5,-0.5)

            glNormal3f(-1.0, 0.0, 0.0)
            glVertex3f(-0.5,-0.5,-0.5)
            glVertex3f(-0.5,-0.5, 0.5)
            glVertex3f(-0.5, 0.5, 0.5)
            glVertex3f(-0.5, 0.5,-0.5)
            glEnd()

            glRotatef(self.lasty - self.y, 1.0, 0.0, 0.0);
            glRotatef(self.lastx - self.x, 0.0, 1.0, 0.0);


            _cmd.p_glut_event(0)
            _cmd.p_glut_event(0)
            _cmd.p_glut_event(1)            
            
            self.SwapBuffers()


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
      return TRUE

app = MyApp(0)
app.MainLoop()

