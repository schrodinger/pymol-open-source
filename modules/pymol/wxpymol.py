
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
P_GLUT_CHAR_EVENT          =  5

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
            EVT_CHAR(self,self.OnChar)
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

        def OnChar(self, evt):
            self.mod = 0
            if evt.ShiftDown():
                self.mod = self.mod + P_GLUT_ACTIVE_SHIFT
            if evt.ControlDown():
                self.mod = self.mod + P_GLUT_ACTIVE_CTRL
            if evt.MetaDown():
                self.mod = self.mod + P_GLUT_ACTIVE_ALT
            _cmd.p_glut_event(P_GLUT_CHAR_EVENT,evt.GetX(),evt.GetY(),evt.GetKeyCode(),0,self.mod)
                              
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
            _cmd.p_glut_event(P_GLUT_RESHAPE_EVENT,800,500,0,0,0) # initial reshape event
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

ID_ABOUT=101  
ID_OPEN=102 
ID_BUTTON1=110 
ID_EXIT=200
        
class MainWindow(wxFrame):
   def __init__(self,parent,id,title):  
      self.dirname='' 
      wxFrame.__init__(self,parent,-4, title, style=wxDEFAULT_FRAME_STYLE|  
                              wxNO_FULL_REPAINT_ON_RESIZE)  
      self.control = wxTextCtrl(self, 1, style=wxTE_MULTILINE)  
      self.CreateStatusBar() # A Statusbar in the bottom of the window  
      # Setting up the menu.  
      filemenu= wxMenu()  
      filemenu.Append(ID_OPEN, "&Open"," Open a file to edit")  
      filemenu.AppendSeparator()  
      filemenu.Append(ID_ABOUT, "&About"," Information about this program")  
      filemenu.AppendSeparator()  
      filemenu.Append(ID_EXIT,"E&xit"," Terminate the program")  
      # Creating the menubar.  
      menuBar = wxMenuBar()  
      menuBar.Append(filemenu,"&File") # Adding the "filemenu" to the MenuBar  
      self.SetMenuBar(menuBar)  # Adding the MenuBar to the Frame content.  
      EVT_MENU(self, ID_ABOUT, self.OnAbout) 
      EVT_MENU(self, ID_EXIT, self.OnExit) 
      EVT_MENU(self, ID_OPEN, self.OnOpen) 

      self.sizer2 = wxBoxSizer(wxHORIZONTAL) 
      self.buttons=[] 
      for i in range(0,6): 
         self.buttons.append(wxButton(self, ID_BUTTON1+i, "Button &"+`i`)) 
         self.sizer2.Add(self.buttons[i],1,wxEXPAND) 

      self.canvas = CubeCanvas(self)

      # Use some sizers to see layout options 
      self.sizer=wxBoxSizer(wxVERTICAL) 
      self.sizer.Add(self.control,1,wxEXPAND) 
      self.sizer.Add(self.sizer2,0,wxEXPAND) 
      self.sizer.Add(self.canvas,1,wxEXPAND)
      
      #Layout sizers 
      self.SetSizer(self.sizer) 
      self.SetAutoLayout(1) 
      self.sizer.Fit(self) 

      self.Show(1)  

   def OnAbout(self,e):  
      d= wxMessageDialog( self, " A sample editor \n"  
                     " in wxPython","About Sample Editor", wxOK)  
                     # Create a message dialog box  
      d.ShowModal() # Shows it  
      d.Destroy() # finally destroy it when finished.  

   def OnExit(self,e):  
      self.Close(true)  # Close the frame.  

   def OnOpen(self,e):  
      """ Open a file""" 
      dlg = wxFileDialog(self, "Choose a file", self.dirname, "", "*.*", wxOPEN)  
      if dlg.ShowModal() == wxID_OK:  
         self.filename=dlg.GetFilename()  
         self.dirname=dlg.GetDirectory()  
         f=open(self.dirname+'\\'+self.filename,'r')  
         self.control.SetValue(f.read())  
         f.close()  
      dlg.Destroy()  
       
   
class MyApp(wxApp):
  def OnInit(self):
      mw = MainWindow(None, -1,"PyMOL")
      import threading
      return TRUE

#      frame = wxFrame(None, -1, "PyMOL", size=(1024,768))
#      frame.canvas = CubeCanvas(frame)
#      frame.control = wxTextCtrl(frame,1,style=wxTE_MULTILINE)
#      frame.sizer = wxBoxSizer(wxVERTICAL)
#      frame.sizer.Add(frame.canvas,1,wxEXPAND)
#      frame.sizer.Add(frame.control,1,wxEXPAND)
#      frame.sizer.Fit(frame)
#      frame.Show(true)
#      self.SetTopWindow(frame)

app = MyApp(0)
app.MainLoop()

