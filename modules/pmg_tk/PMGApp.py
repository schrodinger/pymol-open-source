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
from tkFileDialog import *

from AbstractApp import AbstractApp
from Setting import Setting

import Pmw
import sys, string
from pymol import cmd
import re
import thread
import threading
import os

class PMGApp(AbstractApp):

   appversion     = '0.43'
   appname       = 'PyMOL Molecular Graphics System'
   copyright      = 'Copyright (C) 1998-2001 by Warren DeLano of\nDeLano Scientific. All rights reserved.'
   contactweb     = 'http://www.pymol.org'
   contactemail   = 'warren@delanoscientific.com'

   def buttonAdd(self,frame,text,cmd):
      newBtn=self.createcomponent('button', (), None,
         Button,frame,text=text,highlightthickness=0,
         command=cmd,padx=0,pady=0)
      newBtn.pack(side=LEFT,fill=BOTH,expand=YES)
      
   def createButtons(self):
      row2 = self.createcomponent('row2', (), None,
         Frame,self.get_commandFrame(),bd=0)
      row2.pack(side=TOP,fill=BOTH,expand=YES)
      btn_reset = self.buttonAdd(row2,'Reset',cmd.reset)
      btn_rtrace = self.buttonAdd(row2,'Ray Trace',lambda : cmd.do("ray"))
      btn_reset = self.buttonAdd(row2,'Rock',cmd.rock)

#      row3 = self.createcomponent('row3', (), None,
#         Frame,self.get_commandFrame(),bd=0)
#      row3.pack(side=TOP,fill=BOTH,expand=YES)
#      btn_reset = self.buttonAdd(row3,'Barf',self.quit)
#      btn_reset = self.buttonAdd(row3,'Now',self.quit)
#      btn_reset = self.buttonAdd(row3,'Eat',self.quit)
#      btn_reset = self.buttonAdd(row3,'Shrimp',self.quit)

      row1 = self.createcomponent('row1', (), None,
         Frame,self.get_commandFrame(),bd=0)
      row1.pack(side=TOP,fill=BOTH,expand=YES)
      btn_rewind = self.buttonAdd(row1,'|<',cmd.rewind)
      btn_back = self.buttonAdd(row1,'<',cmd.backward)
      btn_play = self.buttonAdd(row1,'Play',cmd.mplay)
      btn_stop = self.buttonAdd(row1,'Stop',cmd.mstop)
      btn_forward = self.buttonAdd(row1,'>',cmd.forward)
      btn_last = self.buttonAdd(row1,'>|',cmd.ending)
      btn_ccache = self.buttonAdd(row1,'MClear',cmd.mclear)


   def createMain(self):
      self.command = StringVar()
      self.entry = self.createcomponent('entry', (), None,
                           Entry,
                           (self.get_dataArea(),),
                           justify=LEFT,
                           width=50,
                           textvariable=self.command)
      self.entry.pack(side=BOTTOM,expand=NO,fill=X)
      self.entry.bind('<Return>',lambda event,w=self.command:
         (cmd.do(w.get()),cmd.dirty(),w.set('')))

      self.output = self.createcomponent('output', (), None,
                           Pmw.ScrolledText,
                           (self.get_dataArea(),))

      text = self.output.component('text')
      if sys.platform!='win32':
         text.configure(font=('Courier',12))         
      else:
         text.configure(font=('Courier',9))
      text.configure(width=72)
      self.output.after(1000,self.update_feedback)
      self.output.pack(side=BOTTOM,expand=YES,fill=BOTH)
      self.bind(self.entry, 'Command Input Area')

   def update_feedback(self):
      for a in cmd.get_feedback():
         self.output.insert(END,"\n")
         self.output.insert(END,a)
         self.output.see(END)
         self.lineCount = self.lineCount + 1
         if self.lineCount > 10000:
            self.output.delete('0.0','%i.%i' % (self.lineCount-5000,0))
            self.lineCount=5000
      self.output.after(50,self.update_feedback)

   def createInterface(self):
         AbstractApp.createInterface(self)
         self.createButtons()
         self.createMain()
         self.lineCount = 0
         
   def quit_app(self):
      cmd.quit()

   def file_open(self):
      ofile = askopenfilename(filetypes=[("PDB File","*.pdb"),
                    ("MOL File","*.mol"),("XPLOR Map","*.xplor"),
                    ("ChemPy Model","*.pkl"),("All Files","*.*")])
      if len(ofile):
         cmd.do("cmd.load('%s')"%ofile)

   def file_run(self):
      ofile = askopenfilename(filetypes=[("PyMOL Script","*.pml"),("Python Program","*.py")])
      if re.search("\.py$",ofile):
         cmd.do("run "+ofile);      
      else:
         cmd.do("@"+ofile);

   def file_save(self):
      sfile = asksaveasfilename(filetypes=[("PNG File","*.png")])
      if len(sfile):
         cmd.png(sfile)

   def file_savemovie(self):
      sfile = asksaveasfilename(filetypes=[("Numbered PNG Files","*.png")])
      if len(sfile):
         cmd.mpng(sfile)

   def createMenuBar(self):
      self.menuBar.addmenuitem('Help', 'command',
                               'Get information on application', 
                               label='About', command = cmd.splash)

      self.menuBar.addmenuitem('Help', 'separator', '')
      
      self.menuBar.addmenuitem('Help', 'command', 'Help on the API',
                               label='API',
                               command = lambda: cmd.show_help("api"))      

      self.menuBar.addmenuitem('Help', 'command', 'Help on Launching',
                               label='Launching',
                               command = lambda: cmd.show_help("launching"))      

      self.menuBar.addmenuitem('Help', 'command', 'Help on Commands',
                               label='Commands',
                               command = lambda: cmd.show_help("commands"))      

      self.menuBar.addmenuitem('Help', 'separator', '')

      self.menuBar.addmenuitem('Help', 'command', 'Help on the API',
                               label='API',
                               command = lambda: cmd.show_help("api"))      

      self.menuBar.addmenuitem('Help', 'command', 'Help on Selections',
                               label='Selections',
                               command = lambda: cmd.show_help("selections"))      

      self.menuBar.addmenuitem('Help', 'command', 'Example Selections',
                               label='Example Selections',
                               command = lambda: cmd.show_help("examples"))      

      self.menuBar.addmenuitem('Help', 'command', 'Help on the Mouse',
                               label='Mouse',
                               command = lambda: cmd.show_help("mouse"))      

      self.menuBar.addmenuitem('Help', 'command', 'Help on the Keyboard',
                               label='Keyboard',
                               command = lambda: cmd.show_help("keyboard"))      

      self.menuBar.addmenuitem('Help', 'command', 'Help on Stereo',
                               label='Stereo',
                               command = lambda: cmd.show_help("stereo"))      


      self.toggleBalloonVar = IntVar()
      self.toggleBalloonVar.set(1)
      self.setting = Setting()

      self.menuBar.addmenuitem('Help', 'separator', '')
      
      self.menuBar.addmenuitem('Help', 'checkbutton',
                         'Toggle balloon help',
                         label='Balloon help',
                        variable = self.toggleBalloonVar,
                        command=self.toggleBalloon)

      self.menuBar.addmenuitem('File', 'command', 'Open structure file.',
                        label='Open...',
                        command=self.file_open)

#      self.menuBar.addmenuitem('File', 'command', 'Open sequential files.',
#                        label='Open Sequence...',
#                        command=self.file_open)

      self.menuBar.addmenuitem('File', 'separator', '')

      self.menuBar.addmenuitem('File', 'command', 'Run program or script.',
                        label='Run...',
                        command=self.file_run)

      self.menuBar.addmenuitem('File', 'separator', '')

      self.menuBar.addmenuitem('File', 'command', 'Save current image.',
                        label='Save Image...',
                        command=self.file_save)

      self.menuBar.addmenuitem('File', 'command', 'Save all frames.',
                        label='Save Movie...',
                        command=self.file_savemovie)

      self.menuBar.addmenuitem('File', 'separator', '')

      self.menuBar.addmenuitem('File', 'command', 'Quit PyMOL',
                        label='Quit',
                        command=cmd.quit)

      self.menuBar.addmenu('Movie', 'Movie Control')

      self.menuBar.addmenuitem('Movie', 'checkbutton',
                         'Photorealistic images.',
                         label='Ray Trace Frames',
                        variable = self.setting.ray_trace_frames,
                        command = lambda s=self: s.setting.update('ray_trace_frames'))

      self.menuBar.addmenuitem('Movie', 'checkbutton',
                         'Save images in memory.',
                         label='Cache Frames',
                        variable = self.setting.cache_frames,
                        command = lambda s=self: s.setting.update('cache_frames'))

      self.menuBar.addmenuitem('Movie', 'command', 'Flush Cache',
                               label='Flush Cache',
                               command = lambda: cmd.mclear())

      self.menuBar.addmenuitem('Movie', 'separator', '')

      self.menuBar.addmenuitem('Movie', 'command', 'Maximum Speed',
                               label='Maximum Speed',
                               command = lambda: cmd.set("movie_delay","0"))

      self.menuBar.addmenuitem('Movie', 'command', '30 FPS',
                               label='30 FPS',
                               command = lambda: cmd.set("movie_delay","33"))

      self.menuBar.addmenuitem('Movie', 'command', '15 FPS',
                               label='15 FPS',
                               command = lambda: cmd.set("movie_delay","66"))

      self.menuBar.addmenuitem('Movie', 'command', '5 FPS',
                               label='5 FPS',
                               command = lambda: cmd.set("movie_delay","200"))

      self.menuBar.addmenuitem('Movie', 'command', '1 FPS',
                               label='1 FPS',
                               command = lambda: cmd.set("movie_delay","1000"))

      self.menuBar.addmenuitem('Movie', 'command', '0.3 FPS',
                               label='0.3 FPS',
                               command = lambda: cmd.set("movie_delay","3000"))

      self.menuBar.addmenuitem('Movie', 'separator', '')

      self.menuBar.addmenuitem('Movie', 'command', 'Reset Meter',
                               label='Reset Meter',
                               command = lambda: cmd.meter_reset())

      self.menuBar.addmenu('Display', 'Display Control')


      self.menuBar.addmenuitem('Display', 'command', 'Clear Text Output',
                               label='Clear Text',
                               command = lambda: cmd.cls())

      self.menuBar.addmenuitem('Display', 'command', 'Hide Text Output',
                               label='Hide Text',
                               command = lambda: cmd.set("text","0"))

      self.menuBar.addmenuitem('Display', 'command', 'Show Text Output',
                               label='Show Text',
                               command = lambda: cmd.set("text","1"))

      self.menuBar.addmenuitem('Display', 'separator', '')
      
      self.menuBar.addmenuitem('Display', 'command', 'Stereo On',
                               label='Stereo On',
                               command = lambda: cmd.stereo("on"))

      self.menuBar.addmenuitem('Display', 'command', 'Stereo Off',
                               label='Stereo Off',
                               command = lambda: cmd.stereo("off"))

      self.menuBar.addmenu('Options', 'Configuration Control')

      self.menuBar.addmenuitem('Options', 'checkbutton',
                         'Show bond valences.',
                         label='Show bond valences',
                        variable = self.setting.valence,
                        command = lambda s=self: s.setting.update('valence'))

      self.menuBar.addmenuitem('Options', 'checkbutton',
                         'Superimpose all molecular states.',
                         label='Show All States',
                        variable = self.setting.all_states,
                        command = lambda s=self: s.setting.update('all_states'))

      self.menuBar.addmenuitem('Options', 'checkbutton',
                         'Disable perspective.',
                         label='Orthoscopic',
                        variable = self.setting.ortho,
                        command = lambda s=self: s.setting.update('ortho'))

      self.menuBar.addmenuitem('Options', 'checkbutton',
                         'Smooth raytracing.',
                         label='Antialias',
                        variable = self.setting.antialias,
                        command = lambda s=self: s.setting.update('antialias'))

      self.menuBar.addmenuitem('Options', 'checkbutton',
                         'Auto Zoom.',
                         label='Auto Zoom',
                        variable = self.setting.auto_zoom,
                        command = lambda s=self: s.setting.update('auto_zoom'))

      self.menuBar.addmenuitem('Options', 'checkbutton',
                         'Overlay',
                         label='Overlay Text on Graphics',
                        variable = self.setting.overlay,
                        command = lambda s=self: s.setting.update('overlay'))

      self.menuBar.addmenuitem('Options', 'checkbutton',
                         'Normals Bug Workaround',
                         label='Normals Bug Workaround',
                        variable = self.setting.normal_workaround,
                        command = lambda s=self: s.setting.update('normal_workaround'))

      self.menuBar.addmenu('Mouse', 'Mouse Configuration')

      self.menuBar.addmenuitem('Mouse', 'command', 'Visualization',
                               label='Visualization',
                               command = lambda: cmd.edit_mode("off"))

      self.menuBar.addmenuitem('Mouse', 'command', 'Editing',
                               label='Editing',
                               command = lambda: cmd.edit_mode("on"))

      self.menuBar.addmenu('Wizards', 'Task Wizards')
      

      self.menuBar.addmenuitem('Wizards', 'command', 'Distances',
                               label='Distances',
                               command = lambda: cmd.wizard("distance"))

      self.menuBar.addmenuitem('Wizards', 'command', 'Pair Fitting',
                               label='Pair Fitting',
                               command = lambda: cmd.wizard("pair_fit"))

      self.menuBar.addmenuitem('Wizards', 'command', 'Charge',
                               label='Charge',
                               command = lambda: cmd.wizard("charge"))



