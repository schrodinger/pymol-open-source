#A* -------------------------------------------------------------------
#B* This file contains source code for the PyMOL computer program
#C* copyright 1998-2000 by Warren Lyford Delano of DeLano Scientific. 
#D* -------------------------------------------------------------------
#E* It is unlawful to modify or remove this copyright notice.
#F* -------------------------------------------------------------------
#G* Please see the accompanying LICENSE file for further information. 
#H* -------------------------------------------------------------------
#I* Additional authors of this source file include:
#-*
#-* NOTE: Based on code by John E. Grayson which was in turn 
#-* based on code written by Doug Hellmann. 
#Z* -------------------------------------------------------------------

from Tkinter import *
import Pmw
import sys, string

class AbstractApp(Pmw.MegaWidget):        
   appversion      = ''
   appname         = ''
   copyright       = ''
   contactweb      = ''
   contactemail    = ''
          
   frameWidth      = 860 # these days, overriden by values in pymol.invocation
   frameHeight     = 120
   frameXPos       = 0
   frameYPos       = 0
   frameAdjust     = 51
   
   padx         = 1
   pady         = 1
   balloonhelp    = 1
   
   busyCursor = 'watch'
   
   def __init__(self, **kw):
      optiondefs = (
         ('padx',         1,               Pmw.INITOPT),
         ('pady',         1,               Pmw.INITOPT))
      self.defineoptions(kw, optiondefs)

      self.root = Tk()
      self.initializeTk(self.root)
      Pmw.initialise(self.root)
      self.root.title(self.appname)

      if sys.platform=='win32':
         self.frameAdjust = 60 # 51 sufficient for Win2k, but 60 needed for XP
      elif sys.platform[0:4]=='irix':
         self.frameAdjust = 41
      elif sys.platform!='linux':
         self.frameAdjust = 31
      else:
         self.frameAdjust = 51
         
      inv = sys.modules.get("pymol.invocation",None)
      if inv!=None:
         self.frameWidth = inv.options.win_x+220
         self.frameXPos = inv.options.win_px
         self.frameHeight = inv.options.ext_y
         self.frameYPos = inv.options.win_py - (
            self.frameHeight + self.frameAdjust)
      self.root.geometry('%dx%d+%d+%d' % (
         self.frameWidth, self.frameHeight,
         self.frameXPos, self.frameYPos))

      # Initialize the base class
      Pmw.MegaWidget.__init__(self, parent=self.root)
      
      # initialize the application
      self.appInit()
      
      # create the interface
      self.__createInterface()
      
      # create a table to hold the cursors for
      # widgets which get changed when we go busy
      self.preBusyCursors = None
      
      # pack the container and set focus
      # to ourselves
      self._hull.pack(side=LEFT, fill=BOTH, expand=YES)
      self.focus_set()

      # initialize our options
      self.initialiseoptions(AbstractApp)
      
   def appInit(self):
      # Called before interface is created (should be overridden).
      pass
      
   def initializeTk(self, root):
      self.pad = ' ' 
      # Initialize platform-specific options
      if sys.platform == 'mac':
         self.__initializeTk_mac(root)
      elif sys.platform[:3] == 'win':
         self.__initializeTk_win32(root)
      elif sys.platform[:5] == 'linux':
         self.__initializeTk_unix(root)
      else:
         self.__initializeTk_unix(root)

   def __initializeTk_colors_common(self, root):
      root.option_add('*background', 'grey')
      root.option_add('*foreground', 'black')
      root.option_add('*EntryField.Entry.background', 'white')
      root.option_add('*Entry.background', 'white')      
      root.option_add('*MessageBar.Entry.background', 'gray85')
      root.option_add('*Listbox*background', 'white')
      root.option_add('*Listbox*selectBackground', 'dark slate blue')
      root.option_add('*Listbox*selectForeground', 'white')
                  
   def __initializeTk_win32(self, root):
      self.__initializeTk_colors_common(root)
      root.option_add('*Font', 'Tahoma 8')
      self.pad = ' '
#      root.option_add('*EntryField.Entry.Font', 'Courier 10')
#      root.option_add('*Listbox*Font', 'Courier 10')
      
   def __initializeTk_mac(self, root):
      self.__initializeTk_colors_common(root)
      
   def __initializeTk_unix(self, root):
      self.__initializeTk_colors_common(root)

   def busyStart(self, newcursor=None):
      if not newcursor:
         newcursor = self.busyCursor
      newPreBusyCursors = {}
      for component in self.busyWidgets:
         newPreBusyCursors[component] = component['cursor']
         component.configure(cursor=newcursor)
         component.update_idletasks()
      self.preBusyCursors = (newPreBusyCursors, self.preBusyCursors)
      
   def busyEnd(self):
      if not self.preBusyCursors:
         return
      oldPreBusyCursors = self.preBusyCursors[0]
      self.preBusyCursors = self.preBusyCursors[1]
      for component in self.busyWidgets:
         try:
            component.configure(cursor=oldPreBusyCursors[component])
         except KeyError:
            pass
         component.update_idletasks()
           
   def __createAboutBox(self):
      Pmw.aboutversion(self.appversion)
      Pmw.aboutcopyright(self.copyright)
      Pmw.aboutcontact(
        'For more information, browse to: %s\n or send email to: %s' %\
                 (self.contactweb, self.contactemail))
      self.about = Pmw.AboutDialog(self._hull, 
                            applicationname=self.appname)
      self.about.withdraw()
      return None
   
   def showAbout(self):
      # Create the dialog to display about and contact information.
      self.about.activate(geometry='centerscreenfirst')
#      self.about.show()
#      self.about.focus_set()
      
   def toggleBalloon(self):
      if self.toggleBalloonVar.get():
         self.__balloon.configure(state = 'both')
      else:
         self.__balloon.configure(state = 'status')

   def __createMenuBar(self):
      self.menuBar = self.createcomponent('menubar', (), None,
                                 Pmw.MenuBar,
                                 (self._hull,),
#                                 hull_relief=RAISED,
#                                 hull_borderwidth=0,
                                 balloon=self.balloon())

      self.menuBar.pack(fill=X)
      self.menuBar.addmenu('Help', 'About %s' % self.appname, side='right')
      self.menuBar.addmenu('File', 'File Input',tearoff=TRUE)
      self.menuBar.addmenu('Edit', 'Text Editing',tearoff=TRUE)
                     
   def createMenuBar(self):
#override this
      pass

   def __createBalloon(self):
      # Create the balloon help manager for the frame.
      # Create the manager for the balloon help
      self.__balloon = self.createcomponent('balloon', (), None,
                                   Pmw.Balloon, (self._hull,))

   def balloon(self):
      return self.__balloon

   def __createDataArea(self):
      # Create data area where data entry widgets are placed.
      self.dataArea = self.createcomponent('dataarea',
                                  (), None,
                                  Frame, (self._hull,), 
                                  relief=SUNKEN, 
                                  bd=1)
      self.dataArea.pack(side=LEFT, fill=BOTH, expand=YES,
                     padx=1, pady=1)

   def __createCommandArea(self):
      # Create a command area for application-wide buttons.
      self.commandFrame = self.createcomponent('commandframe', (), None,
         Frame,(self._hull,),relief=SUNKEN,bd=1)
      self.commandFrame.place(width=500)
      self.commandFrame.pack(side=TOP, 
                   expand=NO, 
                   fill=BOTH,
                   padx=1,
                   pady=1)

   def __createMessageBar(self):
      # Create the message bar area for help and status messages.
      frame = self.createcomponent('bottomtray', (), None,
                            Frame,(self._hull,), relief=SUNKEN)
      self.__messageBar = self.createcomponent('messagebar',
                                      (), None,
                                     Pmw.MessageBar, 
                                     (frame,),
                                     #entry_width = 40,
                                     entry_relief=SUNKEN,
                                     entry_bd=1,
                                     labelpos=None)
      self.__messageBar.pack(side=LEFT, expand=YES, fill=X)

#      self.__progressBar = ProgressBar.ProgressBar(frame,
#                                    fillColor='slateblue',
#                                    doLabel=1,
#                                    width=150)
#      self.__progressBar.frame.pack(side=LEFT, expand=NO, fill=NONE)

#      self.updateProgress(0)
      frame.pack(side=BOTTOM, expand=NO, fill=X)
               
      self.__balloon.configure(statuscommand = \
                         self.__messageBar.helpmessage)

   def messageBar(self):
      return self.__messageBar

#   def updateProgress(self, newValue=0, newMax=0):
#      self.__progressBar.updateProgress(newValue, newMax)

   def bind(self, child, balloonHelpMsg, statusHelpMsg=None):
      # Bind a help message and/or status message to a widget.
      self.__balloon.bind(child, balloonHelpMsg, statusHelpMsg)

   def get_dataArea(self):
      # Retrieve the interior site where widgets should go.
      return self.dataArea

   def get_commandFrame(self):
      # Retrieve the command frame where buttons go.
      return self.commandFrame

   def __createInterface(self):
      self.__createBalloon()
      self.__createMenuBar()
      self.__createDataArea()
      self.__createCommandArea()
      self.__createMessageBar()
      self.__createAboutBox()
      #
      # Create the parts of the interface
      # which can be modified by subclasses
      #
      self.busyWidgets = ( self.root, )
      self.createMenuBar()
      self.createInterface()

   def createInterface(self):
      # Override this method to create the interface for the app.
      pass
      
   def quit_app(self):
      # Override this method to create the interface for the app.
      pass
      
   def main(self):
      # This method should be left intact!
      self.pack()
      self.mainloop()
      self.quit_app()
      
   def run(self):
      self.main()
