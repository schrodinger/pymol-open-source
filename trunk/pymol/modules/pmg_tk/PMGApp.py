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
from SetEditor import SetEditor

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
      btn_stop = self.buttonAdd(row1,'Stop',cmd.mstop)
      btn_play = self.buttonAdd(row1,'Play',cmd.mplay)
      btn_forward = self.buttonAdd(row1,'>',cmd.forward)
      btn_last = self.buttonAdd(row1,'>|',cmd.ending)
      btn_ccache = self.buttonAdd(row1,'MClear',cmd.mclear)

   def my_show(self,win,center=1):
      if sys.platform!='linux2':
         win.show()
      else: # autocenter, deiconify, and run mainloop
         # this is a workaround for a bug in the
         # interaction between Tcl/Tk and common Linux
         # window managers (namely KDE/Gnome) which causes
         # an annoying 1-2 second delay in opening windows!
         if center:
            tw = win.winfo_reqwidth()+100
            th = win.winfo_reqheight()+100
            vw = win.winfo_vrootwidth()
            vh = win.winfo_vrootheight()
            x = max(0,(vw-tw)/2)
            y = max(0,(vh-tw)/2)
            win.geometry(newGeometry="+%d+%d"%(x,y))
         win.deiconify()
#         win.show()
         
   def my_withdraw(self,win):
      if sys.platform!='linux2':
         win.withdraw()
      else: # autocenter, deiconify, and run mainloop
         win.destroy()

   def my_activate(self,win,center=1,focus=None):
      if sys.platform!='linux2':
         win.activate()
      else: # autocenter, deiconify, and run mainloop
         # this is a workaround for a bug in the
         # interaction between Tcl/Tk and common Linux
         # window managers (namely KDE/Gnome) which causes
         # an annoying 1-2 second delay in opening windows!
         if center:
            tw = win.winfo_reqwidth()+100
            th = win.winfo_reqheight()+100
            vw = win.winfo_vrootwidth()
            vh = win.winfo_vrootheight()
            x = max(0,(vw-tw)/2)
            y = max(0,(vh-tw)/2)
            win.geometry(newGeometry="+%d+%d"%(x,y))
         win.deiconify()
         if focus!=None:
            focus.focus_set()
         win.mainloop()
         
   def my_deactivate(self,win):
      if sys.platform!='linux2':
         win.deactivate()
      else: # autocenter, deiconify, and run mainloop
         win.destroy()
      
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
      if sys.platform=='linux2':
         self.my_fw_font=('Courier',12)
      elif sys.platform=='win32':
         self.my_fw_font=('Courier',9)
      else:
         self.my_fw_font=('Courier',9)
         
      text.configure(font = self.my_fw_font)
      text.configure(width=72)
      self.output.after(1000,self.update_feedback)
      self.output.after(1000,self.update_menus)
      self.output.pack(side=BOTTOM,expand=YES,fill=BOTH)
      self.bind(self.entry, 'Command Input Area')
      self.initialdir = os.getcwd()

   def update_feedback(self):
      for a in cmd.get_feedback():
         self.output.insert(END,"\n")
         self.output.insert(END,a)
         self.output.see(END)
         self.lineCount = self.lineCount + 1
         if self.lineCount > 10000:
            self.output.delete('0.0','%i.%i' % (self.lineCount-5000,0))
            self.lineCount=5000
      self.output.after(50,self.update_feedback) # update feedback window 20X a second

   def update_menus(self):
      self.setting.refresh()
      self.output.after(500,self.update_menus) # update menus twice a second
      
   def createInterface(self):
         AbstractApp.createInterface(self)
         self.createButtons()
         self.createMain()
         self.lineCount = 0
         
   def quit_app(self):
      cmd.quit()

   def file_open(self):
      ofile = askopenfilename(initialdir = self.initialdir,
                              filetypes=[("PDB File","*.pdb"),
                                         ("All Files","*.*"),
                                         ("PDB File","*.ent"),
                                         ("Macromodel File","*.dat"),
                                         ("Macromodel File","*.out"),
                                         ("Macromodel File","*.mmd"),
                                         ("Macromodel File","*.mmod"),
                                         ("MOL File","*.mol"),
                                         ("ChemPy Model","*.pkl"),
                                         ("Raster3D Scene","*.r3d"),
                                         ("SDF File","*.sdf"),
                                         ("XPLOR Map","*.xplor"),
                                         ("Tinker XYZ File","*.xyz")
                                         ])
      if len(ofile):
         cmd.load(ofile)

   def file_save(self):
      lst = cmd.get_names('all')
      lst = filter(lambda x:x[0]!="_",lst)
      self.dialog = Pmw.SelectionDialog(self.root,title="Save",
                          buttons = ('OK', 'Cancel'),
                                   defaultbutton='OK',
                          scrolledlist_labelpos=N,
                          label_text='Which object or selection would you like to save?',
                          scrolledlist_items = lst,
                          command = self.file_save2)
      if len(lst):
         listbox = self.dialog.component('scrolledlist')      
         listbox.selection_set(0)
      self.my_show(self.dialog)
      
   def file_save2(self,result):
      if result!='OK':
         self.my_withdraw(self.dialog)
         del self.dialog
      else:
         sels = self.dialog.getcurselection()
         if len(sels)!=0:
            sfile = sels[0]+".pdb"
            self.my_withdraw(self.dialog)
            del self.dialog
            if result=='OK':
               sfile = asksaveasfilename(initialfile = sfile,
                                         initialdir = self.initialdir,
                                         filetypes=[("PDB File","*.pdb"),
                                                    ("MOL File","*.mol"),
                                                    ("MMD File","*.mmd"),
                                                    ("PKL File","*.pkl"),
                                                    ])
               if len(sfile):
                  self.initialdir = re.sub(r"[^\/\\]*$","",sfile)
                  cmd.save(sfile,"(%s)"%sels[0])
         
   def file_run(self):
      ofile = askopenfilename(initialdir = os.getcwd(),
                   filetypes=[("PyMOL Script","*.pml"),
                              ("Python Program","*.py"),
                              ("Python Program","*.pyc"),
                              ("PyMOL Program","*.pym")])
      if len(ofile):
         dir = re.sub(r"[^\/\\]*$","",ofile)
         os.chdir(dir)	
         if re.search("\.py$",ofile):
            cmd.do("run "+ofile);      
         else:
            cmd.do("@"+ofile);

   def file_savepng(self):
      sfile = asksaveasfilename(initialdir = self.initialdir,
             filetypes=[("PNG File","*.png")])
      if len(sfile):
         self.initialdir = re.sub(r"[^\/\\]*$","",sfile)
         cmd.png(sfile)
         
      
   def file_savemovie(self):
      sfile = asksaveasfilename(filetypes=[("Numbered PNG Files","*.png")])
      if len(sfile):
         self.initialdir = re.sub(r"[^\/\\]*$","",sfile)
         cmd.mpng(sfile)

   def demo1(self):
      cmd.disable()
      cmd.do("cd $PYMOL_PATH")
      cmd.delete("pept")
      cmd.delete("pept_dist")
      cmd.load("test/dat/pept.pdb")
      cmd.show("sticks","(pept and not i;5:7)")
      cmd.show("surface","(pept and i;5,6)")
      cmd.show("mesh","(pept and i;1,11,12,13)")
      cmd.show("spheres","(pept and i;2,12,9,4 and not n;c,n,o,ca)")
      cmd.show("dots","(i;8)")
      cmd.dist("pept_dist","(i;1&n;OD2)","(i;13&n;OG1)")
      cmd.set("dot_width","2");

   def demo2(self):
      cmd.disable()
      cmd.do("cd $PYMOL_PATH")
      cmd.delete("cgo1")
      cmd.delete("cgo2")
      cmd.do("cd $PYMOL_PATH")
      cmd.load("test/dat/pept.r3d","cgo1")
      cmd.load("test/dat/3al1.r3d","cgo2")
      cmd.zoom()

   def demo3(self):
      cmd.disable()
      cmd.do("cd $PYMOL_PATH")
      cmd.do("run examples/devel/cgo03.py")

   def demo4(self):
      cmd.disable()
      cmd.delete("arg")
      cmd.fragment("arg")
      cmd.zoom("arg",2)
      cmd.show("sticks","arg")
      cmd.feedback('dis','sel','res')
      for a in xrange(1,181):
         cmd.edit("(arg and n;cd)","(arg and n;cg)")
         cmd.torsion("6")
         cmd.unpick()
         cmd.edit("(arg and n;cb)","(arg and n;ca)")
         cmd.torsion("2")
         cmd.unpick()
         cmd.refresh()
      cmd.feedback('ena','sel','res')

   def performance(self,mode):
      if mode==0: # maximum quality
         cmd.set('line_smooth',1)
         cmd.set('depth_cue',1)         
         cmd.set('specular',1)
         cmd.set('surface_quality',1)
         cmd.set('stick_quality',15)
         cmd.set('sphere_quality',2)
         cmd.do("rebuild")
      elif mode==33:
         cmd.set('line_smooth',1)         
         cmd.set('depth_cue',1)         
         cmd.set('specular',1)
         cmd.set('surface_quality',0)
         cmd.set('stick_quality',8)
         cmd.set('sphere_quality',1)
         cmd.do("rebuild")
      elif mode==66: # good perfomance
         cmd.set('line_smooth',0)
         cmd.set('depth_cue',0)         
         cmd.set('specular',1)
         cmd.set('surface_quality',0)
         cmd.set('stick_quality',8)
         cmd.set('sphere_quality',1)
         cmd.do("rebuild")         
      else: # maximum performance
         cmd.set('line_smooth',0)
         cmd.set('depth_cue',0)
         cmd.set('specular',0)
         cmd.set('surface_quality',0)
         cmd.set('stick_quality',5)
         cmd.set('sphere_quality',0)
         cmd.do("rebuild")         
         
   def createMenuBar(self):
      self.menuBar.addmenuitem('Help', 'command',
                               'Get information on application', 
                               label='About', command = cmd.splash)

      self.menuBar.addmenuitem('Help', 'command', 'Release Notes',
                               label='Release Notes',
                               command = lambda: cmd.show_help("release"))      

      self.menuBar.addmenuitem('Help', 'separator', '')
      

      self.menuBar.addmenuitem('Help', 'command', 'Help on Commands',
                               label='Commands',
                               command = lambda: cmd.show_help("commands"))      

      self.menuBar.addmenuitem('Help', 'command', 'Help on Launching',
                               label='Launching',
                               command = lambda: cmd.show_help("launching"))      

      self.menuBar.addmenuitem('Help', 'separator', '')

      self.menuBar.addmenuitem('Help', 'command', 'Help on Selections',
                               label='Select Command',
                               command = lambda: cmd.show_help("select"))      

      self.menuBar.addmenuitem('Help', 'command', 'Help on Selections',
                               label='Selection Syntax',
                               command = lambda: cmd.show_help("selections"))      

      self.menuBar.addmenuitem('Help', 'command', 'Example Selections',
                               label='Selection Examples',
                               command = lambda: cmd.show_help("examples"))      

      self.menuBar.addmenuitem('Help', 'separator', '')
      

      self.menuBar.addmenuitem('Help', 'command', 'Help on the Mouse',
                               label='Mouse',
                               command = lambda: cmd.show_help("mouse"))      

      self.menuBar.addmenuitem('Help', 'command', 'Help on the Keyboard',
                               label='Keyboard',
                               command = lambda: cmd.show_help("keyboard"))      

      self.menuBar.addmenuitem('Help', 'command', 'Help on Molecular Editing',
                               label='Molecular Editing',
                               command = lambda: cmd.show_help("editing"))      

      self.menuBar.addmenuitem('Help', 'command', 'Help on Molecular Editing',
                               label='Molecular Editing Keys',
                               command = lambda: cmd.show_help("edit_keys"))      

      self.menuBar.addmenuitem('Help', 'command', 'Help on Stereo',
                               label='Stereo',
                               command = lambda: cmd.show_help("stereo"))      

      self.menuBar.addmenuitem('Help', 'separator', '')
      

      self.menuBar.addmenuitem('Help', 'command', 'Help on the API',
                               label='API',
                               command = lambda: cmd.show_help("api"))      

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

      self.menuBar.addmenuitem('File', 'command', 'Save structure file.',
                        label='Save Molecule...',
                        command=self.file_save)

#      self.menuBar.addmenuitem('File', 'command', 'Open sequential files.',
#                        label='Open Sequence...',
#                        command=self.file_open)

      self.menuBar.addmenuitem('File', 'command', 'Save current image.',
                        label='Save Image...',
                        command=self.file_savepng)

      self.menuBar.addmenuitem('File', 'command', 'Save all frames.',
                        label='Save Movie...',
                        command=self.file_savemovie)

      self.menuBar.addmenuitem('File', 'separator', '')
      
      self.menuBar.addmenuitem('File', 'command', 'Run program or script.',
                        label='Run...',
                        command=self.file_run)

      self.menuBar.addmenuitem('File', 'separator', '')

      self.menuBar.addmenuitem('File', 'command', 'Quit PyMOL',
                        label='Quit',
                        command=cmd.quit)

      self.menuBar.addmenuitem('Edit', 'command',
                         'To Copy: Use Ctrl-C',
                         label='To copy text use Ctrl-C',
                               state='disabled',
                        command =  None)

      self.menuBar.addmenuitem('Edit', 'command',
                         'To Paste, Use Ctrl-V',
                         label='To paste text use Ctrl-V',
                               state='disabled',                               
                        command =  None)

      self.menuBar.addmenu('Movies', 'Movie Control')

      self.menuBar.addmenuitem('Movies', 'checkbutton',
                         'Photorealistic images.',
                         label='Ray Trace Frames',
                        variable = self.setting.ray_trace_frames,
                        command = lambda s=self: s.setting.update('ray_trace_frames'))

      self.menuBar.addmenuitem('Movies', 'checkbutton',
                         'Save images in memory.',
                         label='Cache Frames',
                        variable = self.setting.cache_frames,
                        command = lambda s=self: s.setting.update('cache_frames'))

      self.menuBar.addmenuitem('Movies', 'command', 'Flush Image Cache',
                               label='Flush Image Cache',
                               command = lambda: cmd.mclear())

      self.menuBar.addmenuitem('Movies', 'separator', '')

      self.menuBar.addmenuitem('Movies', 'checkbutton',
                         'Static Singleton Objects.',
                         label='Static Singleton Objects',
                        variable = self.setting.static_singletons,
                        command = lambda s=self: s.setting.update('static_singletons'))

      self.menuBar.addmenuitem('Movies', 'checkbutton',
                         'Superimpose all molecular states.',
                         label='Show All States',
                        variable = self.setting.all_states,
                        command = lambda s=self: s.setting.update('all_states'))

      self.menuBar.addmenuitem('Movies', 'separator', '')
      
      self.menuBar.addmenuitem('Movies', 'command', 'Maximum Speed',
                               label='Maximum Speed',
                               command = lambda: cmd.set("movie_delay","0"))

      self.menuBar.addmenuitem('Movies', 'command', '30 FPS',
                               label='30 FPS',
                               command = lambda: cmd.set("movie_delay","33"))

      self.menuBar.addmenuitem('Movies', 'command', '15 FPS',
                               label='15 FPS',
                               command = lambda: cmd.set("movie_delay","66"))

      self.menuBar.addmenuitem('Movies', 'command', '5 FPS',
                               label='5 FPS',
                               command = lambda: cmd.set("movie_delay","200"))

      self.menuBar.addmenuitem('Movies', 'command', '1 FPS',
                               label='1 FPS',
                               command = lambda: cmd.set("movie_delay","1000"))

      self.menuBar.addmenuitem('Movies', 'command', '0.3 FPS',
                               label='0.3 FPS',
                               command = lambda: cmd.set("movie_delay","3000"))

      self.menuBar.addmenuitem('Movies', 'separator', '')

      self.menuBar.addmenuitem('Movies', 'command', 'Reset Meter',
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

      self.menuBar.addmenuitem('Display', 'separator', '')

      self.menuBar.addmenuitem('Display', 'command', 'Maximum Performance',
                               label='Maximum Performance',
                               command = lambda s=self: s.performance(100))

      self.menuBar.addmenuitem('Display', 'command', 'Reasonable Performance',
                               label='Reasonable Performance',
                               command = lambda s=self: s.performance(66))

      self.menuBar.addmenuitem('Display', 'command', 'Reasonable Quality',
                               label='Reasonable Quality',
                               command = lambda s=self: s.performance(33))

      self.menuBar.addmenuitem('Display', 'command', 'Maximum Quality',
                               label='Maximum Quality',
                               command = lambda s=self: s.performance(0))

      self.menuBar.addmenu('Settings', 'Configuration Control')


      self.menuBar.addmenuitem('Settings', 'command',
                         'Edit PyMOL Settings',
                         label='Edit All...',
                               command = lambda s=self: SetEditor(s))

      self.menuBar.addmenuitem('Settings', 'separator', '')
      
      self.menuBar.addmenuitem('Settings', 'checkbutton',
                         'Show Valences.',
                         label='Show Valences',
                        variable = self.setting.valence,
                        command = lambda s=self: s.setting.update('valence'))


      self.menuBar.addmenuitem('Settings', 'checkbutton',
                         'Disable perspective.',
                         label='Orthoscopic View',
                        variable = self.setting.ortho,
                        command = lambda s=self: s.setting.update('ortho'))

      self.menuBar.addmenuitem('Settings', 'checkbutton',
                         'Smooth Lines.',
                         label='Smooth Lines',
                        variable = self.setting.line_smooth,
                        command = lambda s=self: s.setting.update('line_smooth'))

      self.menuBar.addmenuitem('Settings', 'checkbutton',
                         'Depth Cue.',
                         label='Depth Cue',
                        variable = self.setting.depth_cue,
                        command = lambda s=self: s.setting.update('depth_cue'))

      self.menuBar.addmenuitem('Settings', 'checkbutton',
                         'Specular Reflections.',
                         label='Specular Reflections',
                        variable = self.setting.specular,
                        command = lambda s=self: s.setting.update('specular'))

      self.menuBar.addmenuitem('Settings', 'checkbutton',
                         'Auto Zoom.',
                         label='Auto Zoom New Objects',
                        variable = self.setting.auto_zoom,
                        command = lambda s=self: s.setting.update('auto_zoom'))


      self.menuBar.addmenuitem('Settings', 'checkbutton',
                         'Overlay',
                         label='Overlay Text on Graphics',
                        variable = self.setting.overlay,
                        command = lambda s=self: s.setting.update('overlay'))

      self.menuBar.addmenuitem('Settings', 'checkbutton',
                         'Smooth raytracing.',
                         label='Antialiased Rendering',
                        variable = self.setting.antialias,
                        command = lambda s=self: s.setting.update('antialias'))

      self.menuBar.addmenuitem('Settings', 'checkbutton',
                         'Cull Backfaces when Rendering',
                         label='Cull Backfaces when Rendering',
                        variable = self.setting.backface_cull,
                        command = lambda s=self: s.setting.update('backface_cull'))

      self.menuBar.addmenu('Mouse', 'Mouse Configuration')

      self.menuBar.addmenuitem('Mouse', 'command', 'Visualization',
                               label='Visualization',
                               command = lambda: cmd.edit_mode("off"))

      self.menuBar.addmenuitem('Mouse', 'command', 'Editing',
                               label='Editing',
                               command = lambda: cmd.edit_mode("on"))

      self.menuBar.addmenu('Wizards', 'Task Wizards')
      

      self.menuBar.addmenuitem('Wizards', 'command', 'Distance',
                               label='Distance',
                               command = lambda: cmd.wizard("distance"))

      self.menuBar.addmenuitem('Wizards', 'command', 'Mutagenesis',
                               label='Mutagenesis',
                               command = lambda: cmd.wizard("mutagenesis"))

      self.menuBar.addmenuitem('Wizards', 'command', 'Pair Fitting',
                               label='Pair Fitting',
                               command = lambda: cmd.wizard("pair_fit"))

      self.menuBar.addmenuitem('Wizards', 'command', 'Charge',
                               label='Charge',
                               command = lambda: cmd.wizard("charge"))

      self.menuBar.addmenu('Demos', 'Demonstrations')

      self.menuBar.addmenuitem('Demos', 'command', 'Representations',
                               label='Representations',
                               command = self.demo1)

      self.menuBar.addmenuitem('Demos', 'command', 'Scripted Animation',
                               label='Scripted Animation',
                               command = self.demo4)

      self.menuBar.addmenuitem('Demos', 'command', 'Molscript/Raster3D Input',
                               label='Molscript/Raster3D Input',
                               command = self.demo2)

      self.menuBar.addmenuitem('Demos', 'command', 'Compiled Graphics Objects',
                               label='Compiled Graphics Objects',
                               command = self.demo3)
