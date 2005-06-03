
from wxPython.wx       import *
from wxPython.glcanvas import *
import _cmd
import threading

import __main__

__main__.pymol_argv = [ 'pymol', '-xq' ]

import pymol
from pymol.embed import EmbeddedPyMOL

code_dict = {
    342 : 'F1',
    343 : 'F2',
    344 : 'F3',
    345 : 'F4',
    346 : 'F5',
    347 : 'F6',
    348 : 'F7',
    349 : 'F8',
    350 : 'F9',
    351 : 'F10',
    352 : 'F11',
    353 : 'F12',
    324 : 'INSERT',
    315 : 'HOME',
    312 : 'PAGE_UP',
    314 : 'END',
    313 : 'PAGE_DOWN',
    316 : 'LEFT',
    317 : 'UP',
    319 : 'DOWN',
    318 : 'RIGHT',
    }

def pst(st):
    print st
    return 1

class MyCanvasBase(wxGLCanvas,EmbeddedPyMOL):
    def __init__(self, parent):

        # wxWANTS_CHARS so we get TAB, etc.
        wxGLCanvas.__init__(self, parent, -1, style=wxWANTS_CHARS)
        self.init = false
        # initial mouse position
        self.lastx = self.x = 30
        self.lasty = self.y = 30
        self.have_focus = 0
        EVT_KILL_FOCUS(self, self.OnKillFocus)
        EVT_SET_FOCUS(self, self.OnSetFocus)
        EVT_CHAR(self, lambda s:pst("char"))
        EVT_ERASE_BACKGROUND(self, self.OnEraseBackground)
        EVT_SIZE(self, self.OnSize)
        EVT_PAINT(self, self.OnPaint)
        EVT_LEFT_DOWN(self, self.OnMouseDown)  
        EVT_LEFT_UP(self, self.OnMouseUp)
        EVT_MIDDLE_DOWN(self, self.OnMouseDown)
        EVT_MIDDLE_UP(self, self.OnMouseUp)
        EVT_MOUSEWHEEL(self, self.OnMouseWheel)
        EVT_RIGHT_DOWN(self, self.OnMouseDown) 
        EVT_RIGHT_UP(self, self.OnMouseUp)
        EVT_MOTION(self, self.OnMouseMotion)
        EVT_IDLE(self,self.OnIdle)
        EVT_CHAR(self,self.OnChar)

    def OnSetFocus(self,event):
        self.have_focus = 1
#      wxGLCanvas.OnSetFocus(event)
        
    def OnKillFocus(self,event):
        self.have_focus = 0
#      wxGLCanvas.OnKillFocus(event)      
        
    def ProcessEvent(self,event):
        pass
        #print "event"
        #wxGLCanvas.ProcessEvent(event)
        
    def OnEraseBackground(self, event):
        pass # Do nothing, to avoid flashing on MSW.

    def OnSize(self, event):
        size = self.GetClientSize()
        if self.GetContext():
            self.SetCurrent()
            self.ep_reshape(int(size.width),int(size.height))

    def OnChar(self, evt):
        code = evt.GetKeyCode()
        if code<256:
            self.ep_char(evt.GetX(),evt.GetY(),code,
                             evt.ShiftDown(),evt.ControlDown(),evt.MetaDown())
        else:
            code = code_dict.get(code)
            if code!=None:
                self.ep_special(evt.GetX(),evt.GetY(),code,
                             evt.ShiftDown(),evt.ControlDown(),evt.MetaDown())
        self.CheckPyMOL()
        
    def OnPaint(self, event):
        #print "paint"
        dc = wxPaintDC(self)
        self.SetCurrent()
        if not self.init:
            self.InitGL()
            self.init = true
        self.OnDraw()

    def OnMouseDown(self, evt):
        if not self.have_focus: # restore keyboard focus
            self.SetFocus()
        self.CaptureMouse()
        x,y = evt.GetPosition()
        self.ep_mouse_down(x,y,
                            evt.LeftIsDown(),evt.MiddleIsDown(),evt.RightIsDown(),
                            evt.ShiftDown(),evt.ControlDown(),evt.MetaDown())
        self.CheckPyMOL()            
      #         print "mouse down"

    def OnMouseUp(self, evt):
        x,y = evt.GetPosition()
        self.ep_mouse_up(x,y)
        self.CheckPyMOL()
      #         print "mouse up"
        self.ReleaseMouse()

    def OnMouseMotion(self, evt):
        if evt.Dragging():
            x,y = evt.GetPosition()
            self.ep_motion(x,y,
                            evt.LeftIsDown(),evt.MiddleIsDown(),evt.RightIsDown(),
                            evt.ShiftDown(),evt.ControlDown(),evt.MetaDown())
            self.CheckPyMOL()
        else:
            x,y = evt.GetPosition()
            self.ep_passive_motion(x,y,
                            evt.ShiftDown(),evt.ControlDown(),evt.MetaDown())
            self.CheckPyMOL()
            
            #            print "mouse motion"

    def OnMouseWheel(self, evt):
        x,y = evt.GetPosition()
        self.ep_wheel(x,y,
                             evt.GetWheelRotation(),
                             evt.ShiftDown(),evt.ControlDown(),evt.MetaDown())
          
    def OnIdle(self,evt):
        self.ep_idle()
        evt.RequestMore(true)
        self.CheckPyMOL()

    def Repaint(self):
        self.AddPendingEvent(wxPaintEvent())

class PyMOLCanvas(MyCanvasBase):
    def InitGL(self):
        self.ep_init()
        self.ep_set_swap_callback(self.swap)
        self.CheckPyMOL()

    def swap(self):
        if self.GetContext():
            self.SetCurrent()
            self.SwapBuffers()

    def OnDraw(self):
        self.ep_draw()

    def CheckPyMOL(self):
        if self.ep_get_redisplay():
            self.Refresh()
            #self.Repaint()

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
        wxFrame.__init__(self,parent,-4, title, size=(900,700),style=wxDEFAULT_FRAME_STYLE|  
                                        wxNO_FULL_REPAINT_ON_RESIZE)

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
        #
        
        self.splitterH = wxSplitterWindow(self, -1, style=wxNO_3D|wxSP_3D)
        self.splitterV = wxSplitterWindow(self.splitterH, -1, style=wxNO_3D|wxSP_3D)
                
#      self.sizer2 = wxBoxSizer(wxHORIZONTAL) 
#      self.buttons=[] 
#      for i in range(0,6): 
#         self.buttons.append(wxButton(self, ID_BUTTON1+i, "Button &"+`i`)) 
#         self.sizer2.Add(self.buttons[i],1,wxEXPAND) 

        def EmptyHandler(evt): pass
        EVT_ERASE_BACKGROUND(self.splitterH, EmptyHandler)
        EVT_ERASE_BACKGROUND(self.splitterV, EmptyHandler)

        self.nb = wxNotebook(self.splitterH, -1, style=wxCLIP_CHILDREN)

        self.splitterH.SplitVertically(self.splitterV,self.nb)

        self.right_box_size = 200
        self.splitterH.SetSashPosition(-self.right_box_size, true)
        self.splitterH.SetMinimumPaneSize(1)

        self.pymol_canvas = PyMOLCanvas(self.splitterV)
        self.pymol = self.pymol_canvas.ep_get_pymol()

        EVT_MOVE(self, self.OnMove)
        EVT_SIZE(self, self.OnSize)
        
        self.control = wxTextCtrl(self.splitterV, -1, style=wxTE_MULTILINE|wxTE_PROCESS_TAB)
        
        self.splitterV.SplitHorizontally(self.control,self.pymol_canvas)
        self.splitterV.SetSashPosition(120, true)
        self.splitterV.SetMinimumPaneSize(20)
        
        # Use some sizers to see layout options 
#      self.sizer=wxBoxSizer(wxVERTICAL) 
#      self.sizer.Add(self.control,1,wxEXPAND) 
#      self.sizer.Add(self.sizer2,0,wxEXPAND) 
#      self.sizer.Add(self.canvas,1,wxEXPAND)
        
        #Layout sizers 
#      self.SetSizer(self.sizer) 
#      self.SetAutoLayout(1) 
#      self.sizer.Fit(self) 

        self.Show(1)  

    def OnSize(self,e):
        size = e.GetSize()
        self.SetSize(e.GetSize())
        self.splitterH.SetSashPosition(size.width-self.right_box_size)
        self.pymol_canvas.Repaint()
        e.Skip()
        
    def OnMove(self,e):
        self.pymol_canvas.Repaint()
        
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
            self.pymol.cmd.load(self.dirname+'\\'+self.filename)
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

