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
from ColorEditor import ColorEditor
from Demo import Demo

import Pmw
import sys, string
import pymol
from pymol import cmd
from pymol import util
from pymol import parser

import re
import thread
import threading
import os
from glob import glob
import __builtin__
import traceback
import Queue

def complete(event,str,widget,self):
   st = parser.complete(str.get())
   if st:
      str.set(st)
      widget.icursor(len(st))
   self.focus_entry = 1
   return 1

class PMGApp(AbstractApp):

   appversion     = '0.86'
   appname       = 'PyMOL Molecular Graphics System'
   copyright      = 'Copyright (C) 1998-2002 by Warren DeLano of\nDeLano Scientific LLC. All rights reserved.'
   contactweb     = 'http://www.pymol.org'
   contactemail   = 'warren@delanoscientific.com'

   def appInit(self): # create a global variable for the external gui
      pymol._ext_gui = self
      self.fifo = Queue.Queue(0)
      self.save_file = ''

   def execute(self,cmmd): 
      self.fifo.put(cmmd)

   def buttonAdd(self,frame,text,cmd):
      newBtn=Button(frame,
                    text=text,highlightthickness=0,
                    command=cmd,padx=0,pady=0)
      newBtn.pack(side=LEFT,fill=BOTH,expand=YES)
      
   def createButtons(self):
      
      row2 = self.createcomponent('row2', (), None,
         Frame,self.get_commandFrame(),bd=0)
      row2.pack(side=TOP,fill=BOTH,expand=YES)
      btn_reset = self.buttonAdd(row2,'Reset',lambda: cmd.do("_ reset"))
      btn_reset = self.buttonAdd(row2,'Zoom',lambda: cmd.do("_ zoom"))      
      btn_rtrace = self.buttonAdd(row2,'Ray',lambda : cmd.do("_ ray"))
      btn_reset = self.buttonAdd(row2,'Rock',lambda :cmd.do("_ rock"))

      row3 = self.createcomponent('row3', (), None,
         Frame,self.get_commandFrame(),bd=0)
      row3.pack(side=TOP,fill=BOTH,expand=YES)
      btn_unpick = self.buttonAdd(row3,'Unpick',lambda : cmd.do("_ unpick"))
      btn_hidesele = self.buttonAdd(row3,'Hide Sele',self.hide_sele)
      btn_getview = self.buttonAdd(row3,'Get View',lambda : cmd.get_view()) # doesn't get logged

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
      btn_rewind = self.buttonAdd(row1,'|<',lambda : cmd.do("_ rewind"))
      btn_back = self.buttonAdd(row1,'<',lambda : cmd.do("_ backward"))
      btn_stop = self.buttonAdd(row1,'Stop',lambda : cmd.do("_ mstop"))
      btn_play = self.buttonAdd(row1,'Play',lambda : cmd.do("_ mplay"))
      btn_forward = self.buttonAdd(row1,'>',lambda : cmd.do("_ forward"))
      btn_last = self.buttonAdd(row1,'>|',lambda : cmd.do("_ ending"))
      btn_ccache = self.buttonAdd(row1,'MClear',lambda : cmd.do("_ mclear"))

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

      self.entry.bind('<Tab>',lambda event,w=self.command,
                      e=self.entry,s=self:
                      complete(event,w,e,s))

#      self.entry.bind('<Up>',lambda event,w=self.command,
#                      e=self.entry,s=self:cmd.do("print 1"))

#      self.entry.bind('<Down>',lambda event,w=self.command,
#                      e=self.entry,s=self:
#                      cmd.do("print 2"))

      self.output = self.createcomponent('output', (), None,
                           Pmw.ScrolledText,
                           (self.get_dataArea(),))

      text = self.output.component('text')
      if sys.platform[:5]=='linux':
         self.my_fw_font=('lucida console',12)
      elif sys.platform[:3]=='win':
         self.my_fw_font=('lucida console',8) # Courier 9
      else:
         self.my_fw_font=('lucida console',10)
                                                                               
      text.configure(font = self.my_fw_font)
      text.configure(width=72)
      self.focus_entry=0
      self.output.after(1000,self.flush_commands)
      self.output.after(1000,self.update_feedback)
      self.output.after(1000,self.update_menus)
      self.output.pack(side=BOTTOM,expand=YES,fill=BOTH)
      self.bind(self.entry, 'Command Input Area')
      self.initialdir = os.getcwd()
      self.log_file = "log.pml"

   def flush_commands(self):
      # flush the external GUI fifo command queue
      while not self.fifo.empty():
         try:
            cmmd = self.fifo.get(0)
            exec cmmd
         except:
            traceback.print_exc()
      self.output.after(20,self.flush_commands) # 50X a second
      
   def update_feedback(self):
      if self.focus_entry:
         self.focus_entry=0
         self.entry.focus_set()
      for a in cmd.get_feedback():
         self.output.insert(END,"\n")
         self.output.insert(END,a)
         self.output.see(END)
         self.lineCount = self.lineCount + 1
         if self.lineCount > 10000:
            self.output.delete('0.0','%i.%i' % (self.lineCount-5000,0))
            self.lineCount=5000
         self.entry.focus_set()
      self.output.after(100,self.update_feedback) # 10X a second

   def update_menus(self):
      self.setting.refresh()
      self.output.after(500,self.update_menus) # twice a second
      
   def createInterface(self):
         AbstractApp.createInterface(self)
         self.createButtons()
         self.createMain()
         self.lineCount = 0
         startup_pattern = re.sub(r"\/[^\/]*$","/startup/*.py*",__file__)
         # startup_pattern = os.environ['PYMOL_PATH']+"/modules/pmg_tk/startup/*.py*"
         raw_list = glob(startup_pattern)
         unique = {}
         for a in raw_list:
            unique[re.sub(r".*\/|\.py.*$","",a)] = 1
         for name in unique.keys():
            if name != "__init__":
               mod_name = "pmg_tk.startup."+name
               __builtin__.__import__(mod_name)
               mod = sys.modules[mod_name]
               mod.__init__(self)

   def quit_app(self):
      cmd.log_close()
      cmd.quit()  # avoid logging this - it is inconvenient...

   def file_open(self):
      ofile = askopenfilename(initialdir = self.initialdir,
                              filetypes=[("All Readable","*.pdb"),
                                         ("All Readable","*.ccp4"),
                                         ("All Readable","*.xplor"),
                                         ("All Readable","*.mol"),                                         
                                         ("All Readable","*.sdf"),
                                         ("All Readable","*.xyz"),                                         
                                         ("All Readable","*.r3d"),
                                         ("All Readable","*.cc1"),
                                         ("All Readable","*.cc2"),                                         
                                         ("All Readable","*.ent"),
                                         ("All Readable","*.dat"),
                                         ("All Readable","*.out"),
                                         ("All Readable","*.mmd"),
                                         ("All Readable","*.mmod"),
                                         ("All Readable","*.pse"),                                         
                                         ("PDB File","*.pdb"),
                                         ("All Files","*.*"),
                                         ("All Files","*"),                                         
                                         ("PDB File","*.ent"),
                                         ("PyMOL Session","*.pse"),
                                         ("CCP4 Map","*.ccp4"),                                         
                                         ("XPLOR Map","*.xplor"),
                                         ("Macromodel File","*.dat"),
                                         ("Macromodel File","*.out"),
                                         ("Macromodel File","*.mmd"),
                                         ("Macromodel File","*.mmod"),
                                         ("MOL File","*.mol"),
                                         ("ChemPy Model","*.pkl"),
                                         ("Raster3D Scene","*.r3d"),
                                         ("SDF File","*.sdf"),
                                         ("ChemDraw3D File","*.cc1"),
                                         ("ChemDraw3D File","*.cc2"),
                                         ("Tinker XYZ File","*.xyz")
                                         ])
      if len(ofile):
         self.initialdir = re.sub(r"[^\/\\]*$","",ofile)         
         cmd.log("load %s\n"%ofile,"cmd.load('%s')\n"%ofile)
         cmd.load(ofile)

   def log_open(self):
      sfile = asksaveasfilename(initialfile = self.log_file,
                                initialdir = self.initialdir,
                                filetypes=[
                                           ("PyMOL Script","*.pml"),
                                           ("PyMOL Program","*.pym"),
                                           ("Python Program","*.py"),
                                           ("All Files","*.*"),
                                           ("All Files","*"),
                                           ])
      if len(sfile):
         self.initialdir = re.sub(r"[^\/\\]*$","",sfile)
         self.log_file = re.sub(r"^.*[^\/\\]","",sfile)
         cmd.log_open(sfile)

   def log_resume(self,append_only=0):
      ofile = askopenfilename(initialdir = os.getcwd(),
                   filetypes=[("All Resumable","*.pml"),
                              ("All Resumable","*.pym"),
                              ("All Resumable","*.py"),
                              ("PyMOL Script","*.pml"),
                              ("PyMOL Program","*.pym"),
                              ("Python Program","*.py"),
                              ("All Files","*.*"),                                           
                              ("All Files","*"),
                              ])
      if len(ofile):
         self.initialdir = re.sub(r"[^\/\\]*$","",ofile)
         self.log_file = re.sub(r"^.*[^\/\\]","",ofile)
         os.chdir(self.initialdir)	         
         cmd.resume(ofile)

   def log_append(self,append_only=0):
      ofile = askopenfilename(initialdir = os.getcwd(),
                   filetypes=[("All Appendable","*.pml"),
                              ("All Appendable","*.pym"),
                              ("All Appendable","*.py"),
                              ("PyMOL Script","*.pml"),
                              ("PyMOL Program","*.pym"),
                              ("Python Program","*.py"),
                              ("All Files","*.*"),                                           
                              ("All Files","*"),
                              ])
      if len(ofile):
         self.initialdir = re.sub(r"[^\/\\]*$","",ofile)
         self.log_file = re.sub(r"^.*[^\/\\]","",ofile)
         os.chdir(self.initialdir)	         
         cmd.log_open(ofile,'a')

   def session_save(self):
      if self.save_file!='':
         cmd.log("save %s,format=pse\n"%(self.save_file),
                 "cmd.save('%s',format='pse')\n"%(self.save_file))
         cmd.save(self.save_file,"","pse")
      else:
         self.session_save_as()

   def session_save_as(self):
      sfile = asksaveasfilename(initialfile = self.save_file,
                                initialdir = self.initialdir,
                                filetypes=[
         ("PyMOL Session File","*.pse"),
         ])
      if len(sfile):
         if re.search(r"\.pse$|\.PSE$",sfile)==None:
            sfile=sfile+".pse"
         self.initialdir = re.sub(r"[^\/\\]*$","",sfile)
         cmd.log("save %s,format=pse\n"%(sfile),
                 "cmd.save('%s',format='pse')\n"%(sfile))
         cmd.save(sfile,"",format='pse')
         self.save_file = sfile
   
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
                                         filetypes=[
                                                    ("PDB File","*.pdb"),
                                                    ("MOL File","*.mol"),
                                                    ("MMD File","*.mmd"),
                                                    ("PKL File","*.pkl"),
                                                    ])
               if len(sfile):
                  self.initialdir = re.sub(r"[^\/\\]*$","",sfile)
                  cmd.log("save %s,(%s)\n"%(sfile,sels[0]),
                          "cmd.save('%s','(%s)')\n"%(sfile,sels[0]))
                  cmd.save(sfile,"(%s)"%sels[0])

   def hide_sele(self):
      cmd.log("util.hide_sele()\n","util.hide_sele()\n")
      util.hide_sele()
         
   def file_run(self):
      ofile = askopenfilename(initialdir = os.getcwd(),
                   filetypes=[("All Runnable","*.pml"),
                              ("All Runnable","*.pym"),
                              ("All Runnable","*.py"),
                              ("All Runnable","*.pyc"),
                              ("PyMOL Script","*.pml"),
                              ("Python Program","*.py"),
                              ("Python Program","*.pyc"),
                              ("PyMOL Program","*.pym"),
                              ("All Files","*.*"),                                           
                              ("All Files","*"),
                              ])
      if len(ofile):
         dir = re.sub(r"[^\/\\]*$","",ofile)
         os.chdir(dir)	
         if re.search("\.pym*$|\.PYM*$",ofile):
            cmd.do("run "+ofile);      
         else:
            cmd.do("@"+ofile);

   def file_savepng(self):
      sfile = asksaveasfilename(initialdir = self.initialdir,
             filetypes=[("PNG File","*.png")])
      if len(sfile):
         self.initialdir = re.sub(r"[^\/\\]*$","",sfile)
         cmd.log("png %s\n"%sfile,"cmd.png('%s')\n"%sfile)
         cmd.png(sfile)
         
      
   def file_savemovie(self):
      sfile = asksaveasfilename(filetypes=[("Numbered PNG Files","*.png")])
      if len(sfile):
         self.initialdir = re.sub(r"[^\/\\]*$","",sfile)
         cmd.log("mpng %s\n"%sfile,"cmd.mpng('%s')\n"%sfile)         
         cmd.mpng(sfile)
   
   def createMenuBar(self):
      self.menuBar.addmenuitem('Help', 'command',
                               'Get information on application', 
                               label='About', command = lambda : cmd.do("_ splash"))

      self.menuBar.addmenuitem('Help', 'command', 'Release Notes',
                               label='Release Notes',
                               command = lambda: cmd.do("_ cmd.show_help('release')"))

      self.menuBar.addmenuitem('Help', 'separator', '')

      self.menuBar.addmenuitem('Help', 'command', 'Help on Commands',
                               label='Commands',
                               command = lambda: cmd.do("_ cmd.show_help('commands')"))

      self.menuBar.addmenuitem('Help', 'command', 'Help on Launching',
                               label='Launching',
                               command = lambda: cmd.do("_ cmd.show_help('launching')"))      

      self.menuBar.addmenuitem('Help', 'separator', '')

      self.menuBar.addmenuitem('Help', 'command', 'Help on Selections',
                               label='Select Command',
                               command = lambda: cmd.do("_ cmd.show_help('select')"))      

      self.menuBar.addmenuitem('Help', 'command', 'Help on Selections',
                               label='Selection Syntax',
                               command = lambda: cmd.do("_ cmd.show_help('selections')"))      

      self.menuBar.addmenuitem('Help', 'command', 'Example Selections',
                               label='Selection Examples',
                               command = lambda: cmd.do("_ cmd.show_help('examples')"))      

      self.menuBar.addmenuitem('Help', 'separator', '')
      

      self.menuBar.addmenuitem('Help', 'command', 'Help on the Mouse',
                               label='Mouse',
                               command = lambda: cmd.do("_ cmd.show_help('mouse')"))      

      self.menuBar.addmenuitem('Help', 'command', 'Help on the Keyboard',
                               label='Keyboard',
                               command = lambda: cmd.do("_ cmd.show_help('keyboard')"))      

      self.menuBar.addmenuitem('Help', 'command', 'Help on Molecular Editing',
                               label='Molecular Editing',
                               command = lambda: cmd.do("_ cmd.show_help('editing')"))      

      self.menuBar.addmenuitem('Help', 'command', 'Help on Molecular Editing',
                               label='Molecular Editing Keys',
                               command = lambda: cmd.do("_ cmd.show_help('edit_keys')"))      

      self.menuBar.addmenuitem('Help', 'command', 'Help on Stereo',
                               label='Stereo',
                               command = lambda: cmd.do("_ cmd.show_help('stereo')"))      

      self.menuBar.addmenuitem('Help', 'separator', '')
      

      self.menuBar.addmenuitem('Help', 'command', 'Help on the API',
                               label='API',
                               command = lambda: cmd.do("_ cmd.show_help('api')"))      

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
                        label=self.pad+'Open...',
                        command=self.file_open)


      self.menuBar.addmenuitem('File', 'command', 'Save session.',
                        label=self.pad+'Save Session',
                        command=self.session_save)

      self.menuBar.addmenuitem('File', 'command', 'Save session.',
                        label=self.pad+'Save Session As...',
                        command=self.session_save_as)

      self.menuBar.addmenuitem('File', 'command', 'Save structure file.',
                        label=self.pad+'Save Molecule...',
                        command=self.file_save)

#      self.menuBar.addmenuitem('File', 'command', 'Open sequential files.',
#                        label=self.pad+'Open Sequence...',
#                        command=self.file_open)

      self.menuBar.addmenuitem('File', 'command', 'Save current image.',
                        label=self.pad+'Save Image...',
                        command=self.file_savepng)

      self.menuBar.addmenuitem('File', 'command', 'Save all frames.',
                        label=self.pad+'Save Movie...',
                        command=self.file_savemovie)

      self.menuBar.addmenuitem('File', 'separator', '')
      
      self.menuBar.addmenuitem('File', 'command', 'Open log file.',
                        label=self.pad+'Log...',
                        command=self.log_open)

      self.menuBar.addmenuitem('File', 'command', 'Resume log file.',
                        label=self.pad+'Resume...',
                        command=self.log_resume)

      self.menuBar.addmenuitem('File', 'command', 'Append log file.',
                        label=self.pad+'Append...',
                        command=self.log_append)

      self.menuBar.addmenuitem('File', 'command', 'Close log file.',
                        label=self.pad+'Close Log',
                        command=cmd.log_close)

      self.menuBar.addmenuitem('File', 'command', 'Run program or script.',
                        label=self.pad+'Run...',
                        command=self.file_run)


      self.menuBar.addmenuitem('File', 'separator', '')

      self.menuBar.addmenuitem('File', 'command', 'Quit PyMOL',
                        label=self.pad+'Quit',
                        command=self.quit)

      self.menuBar.addmenuitem('File', 'separator', '')
      
      self.menuBar.addmenuitem('File', 'checkbutton',
                         'Log Conformations.',
                         label=self.pad+'Log Conformations',
                        variable = self.setting.log_conformations,
                        command = lambda s=self: s.setting.update('log_conformations'))

      self.menuBar.addmenuitem('File', 'checkbutton',
                         'Log Box Selections.',
                         label=self.pad+'Log Box Selections',
                        variable = self.setting.log_box_selections,
                        command = lambda s=self: s.setting.update('log_box_selections'))

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

      self.menuBar.addmenuitem('Edit', 'separator', '')

      self.menuBar.addmenuitem('Edit', 'command', 'Undo Conformation',
                               label='Undo Conformation [Ctrl-Z]',
                               command = lambda: cmd.do("_ undo"))

      self.menuBar.addmenuitem('Edit', 'command', 'Redo Conformation',
                               label='Redo Conformation [Ctrl-A]',
                               command = lambda: cmd.do("_ redo"))

      self.menuBar.addmenuitem('Edit', 'separator', '')
      
      self.menuBar.addmenuitem('Edit', 'command', 'Cycle Bond Valence',
                               label='Cycle Bond Valence [Ctrl-W]',
                               command = lambda: cmd.do("_ cycle_valence"))

      self.menuBar.addmenuitem('Edit', 'command', 'Fill Hydrogens',
                               label='Fill Hydrogens on (pk1) [Ctrl-R]',
                               command = lambda: cmd.do("_ h_fill"))

      self.menuBar.addmenuitem('Edit', 'command', 'Invert',
                               label='Invert (lb)-(pk1)-(rb) [Ctrl-E]',
                               command = lambda: cmd.do("_ invert"))

      self.menuBar.addmenuitem('Edit', 'command', 'Form Bond',
                               label='Create Bond (lb)-(rb) [Ctrl-T]',
                               command = lambda: cmd.do("_ bond"))


      self.menuBar.addmenuitem('Edit', 'separator', '')

      
      self.menuBar.addmenuitem('Edit', 'command', 'Remove (pk1)',
                               label='Remove (pk1) [Ctrl-D]',
                               command = lambda: cmd.do("_ remove pk1"))

      self.menuBar.addmenuitem('Edit', 'command', 'Remove (pkfrag1)',
                               label='Remove (pkfrag1) [Ctrl-X]',
                               command = lambda: cmd.do("_ remove pkfrag1"))

      self.menuBar.addmenuitem('Edit', 'command', 'Remove (pkchain)',
                               label='Remove (pkchain)',
                               command = lambda: cmd.do("_ remove pkchain"))

      self.menuBar.addmenuitem('Edit', 'separator', '')
      
      self.menuBar.addmenuitem('Edit', 'command', 'Make Positive',
                               label='Make (pk1) Positive [Ctrl-K]',
                               command = lambda: cmd.do("_ alter pk1,formal_charge=1.0"))

      self.menuBar.addmenuitem('Edit', 'command', 'Make Negative',
                               label='Make (pk1) Negative [Ctrl-J]',
                               command = lambda: cmd.do("_ alter pk1,formal_charge=-1.0"))

      self.menuBar.addmenuitem('Edit', 'command', 'Make Neutral',
                               label='Make (pk1) Neutral',
                               command = lambda: cmd.do("_ alter pk1,formal_charge=-0.0"))


      self.menuBar.addmenu('Fragment', 'Fragment')

      self.menuBar.addmenuitem('Fragment', 'command', 'Acetylene',
                               label='Acetylene [Alt-J]',
                               command = lambda: cmd.do(
         "_ editor.attach_fragment('pk1','acetylene',2,0)"))

      self.menuBar.addmenuitem('Fragment', 'command', 'Amide N->C',
                               label='Amide N->C [Alt-1]',
                               command = lambda: cmd.do(
         "_ editor.attach_fragment('pk1','formamide',3,1)"))

      self.menuBar.addmenuitem('Fragment', 'command', 'Amide C->N',
                               label='Amide C->N [Alt-2]',
                               command = lambda: cmd.do(
         "_ editor.attach_fragment('pk1','formamide',5,0)"))

      self.menuBar.addmenuitem('Fragment', 'command', 'Bromine',
                               label='Bromine [Ctrl-B]',
                               command = lambda: cmd.do("_ replace Br,1,1"))

      self.menuBar.addmenuitem('Fragment', 'command', 'Carbon',
                               label='Carbon [Ctrl-C]',
                               command = lambda: cmd.do("_ replace C,4,4"))

      self.menuBar.addmenuitem('Fragment', 'command', 'Carbonyl',
                               label='Carbonyl [Alt-0]',
                               command = lambda: cmd.do(
         "_ editor.attach_fragment('pk1','formaldehyde',2,0)"))

      self.menuBar.addmenuitem('Fragment', 'command', 'Chlorine',
                               label='Chlorine [Ctrl-L]',
                               command = lambda: cmd.do("_ replace Cl,1,1"))

      self.menuBar.addmenuitem('Fragment', 'command', 'Cyclobutyl',
                               label='Cyclobutyl [Alt-4]',
                               command = lambda: cmd.do(
         "_ editor.attach_fragment('pk1','cyclobutane',4,0)"))

      self.menuBar.addmenuitem('Fragment', 'command', 'Cyclopentyl',
                               label='Cyclopentyl [Alt-5]',
                               command = lambda: cmd.do(
         "_ editor.attach_fragment('pk1','cyclopentane',5,0)"))

      self.menuBar.addmenuitem('Fragment', 'command', 'Cyclopentadiene',
                               label='Cyclopentadiene [Alt-8]',
                               command = lambda: cmd.do(
           "_ editor.attach_fragment('pk1','cyclopentadiene',5,0)"))

      self.menuBar.addmenuitem('Fragment', 'command', 'Cyclohexyl',
                               label='Cyclohexyl [Alt-6]',
                               command = lambda: cmd.do(
         "_ editor.attach_fragment('pk1','cyclohexane',7,0)"))

      self.menuBar.addmenuitem('Fragment', 'command', 'Cycloheptyl',
                               label='Cycloheptyl [Alt-7]',
                               command = lambda: cmd.do(
         "_ editor.attach_fragment('pk1','cycloheptane',8,0)"))

      self.menuBar.addmenuitem('Fragment', 'command', 'Fluorine',
                               label='Fluorine [Ctrl-F]',
                               command = lambda: cmd.do("_ replace F,1,1"))

      self.menuBar.addmenuitem('Fragment', 'command', 'Iodine',
                               label='Iodine [Ctrl-I]',
                               command = lambda: cmd.do("_ replace I,1,1"))

      self.menuBar.addmenuitem('Fragment', 'command', 'Methane',
                               label='Methane',
                               command = lambda: cmd.do(
         "_ editor.attach_fragment('pk1','methane',1,0)"))

      self.menuBar.addmenuitem('Fragment', 'command', 'Nitrogen',
                               label='Nitrogen [Ctrl-N]',
                               command = lambda: cmd.do("_ replace N,4,3"))

      self.menuBar.addmenuitem('Fragment', 'command', 'Oxygen',
                               label='Oxygen [Ctrl-O]',
                               command = lambda: cmd.do("_ replace O,4,2"))

      self.menuBar.addmenuitem('Fragment', 'command', 'Phenyl',
                               label='Phenyl [Alt-9]',
                               command = lambda: cmd.do(
         "_ editor.attach_fragment('pk1','benzene',6,0)"))

      self.menuBar.addmenuitem('Fragment', 'command', 'Sulfer',
                               label='Sulfer [Ctrl-S]',
                               command = lambda: cmd.do("_ replace S,2,2"))

      self.menuBar.addmenuitem('Fragment', 'command', 'Sulfonyl',
                               label='Sulfonyl [Alt-3]',
                               command = lambda: cmd.do(
         "_ editor.attach_fragment('pk1','sulfone',3,1)"))

      self.menuBar.addmenuitem('Fragment', 'command', 'Phosphorus',
                               label='Phosphorus [Ctrl-P]',
                               command = lambda: cmd.do("_ replace P,4,3"))

      self.menuBar.addmenu('Residue', 'Residue')

      self.menuBar.addmenuitem('Residue', 'command', 'Acetyl',
                               label='Acetyl [Alt-B]',
                               command = lambda: cmd.do(
         "_ editor.attach_amino_acid('pk1','ace')"))

      self.menuBar.addmenuitem('Residue', 'command', 'Alanine',
                               label='Alanine [Alt-A]',
                               command = lambda: cmd.do(
         "_ editor.attach_amino_acid('pk1','ala')"))

      self.menuBar.addmenuitem('Residue', 'command', 'Aspartate',
                               label='Aspartate [Alt-D]',
                               command = lambda: cmd.do(
         "_ editor.attach_amino_acid('pk1','asp')"))

      self.menuBar.addmenuitem('Residue', 'command', 'Asparagine',
                               label='Asparagine [Alt-N]',
                               command = lambda: cmd.do(
         "_ editor.attach_amino_acid('pk1','asn')"))

      self.menuBar.addmenuitem('Residue', 'command', 'Arginine',
                               label='Arginine [Alt-R]',
                               command = lambda: cmd.do(
         "_ editor.attach_amino_acid('pk1','arg')"))

      self.menuBar.addmenuitem('Residue', 'command', 'Cysteine',
                               label='Cysteine [Alt-C]',
                               command = lambda: cmd.do(
         "_ editor.attach_amino_acid('pk1','cys')"))

      self.menuBar.addmenuitem('Residue', 'command', 'Glutamate',
                               label='Glutamate [Alt-E]',
                               command = lambda: cmd.do(
         "_ editor.attach_amino_acid('pk1','glu')"))

      self.menuBar.addmenuitem('Residue', 'command', 'Glutamine',
                               label='Glutamine [Alt-N]',
                               command = lambda: cmd.do(
         "_ editor.attach_amino_acid('pk1','gln')"))

      self.menuBar.addmenuitem('Residue', 'command', 'Glycine',
                               label='Glycine [Alt-G]',
                               command = lambda: cmd.do(
         "_ editor.attach_amino_acid('pk1','gly')"))

      self.menuBar.addmenuitem('Residue', 'command', 'Histidine',
                               label='Histidine [Alt-H]',
                               command = lambda: cmd.do(
         "_ editor.attach_amino_acid('pk1','his')"))

      self.menuBar.addmenuitem('Residue', 'command', 'Isoleucine',
                               label='Isoleucine [Alt-I]',
                               command = lambda: cmd.do(
         "_ editor.attach_amino_acid('pk1','ile')"))

      self.menuBar.addmenuitem('Residue', 'command', 'Leucine',
                               label='Leucine [Alt-L]',
                               command = lambda: cmd.do(
         "_ editor.attach_amino_acid('pk1','leu')"))

      self.menuBar.addmenuitem('Residue', 'command', 'Lysine',
                               label='Lysine [Alt-K]',
                               command = lambda: cmd.do(
         "_ editor.attach_amino_acid('pk1','lys')"))

      self.menuBar.addmenuitem('Residue', 'command', 'Methionine',
                               label='Methionine [Alt-M]',
                               command = lambda: cmd.do(
         "_ editor.attach_amino_acid('pk1','met')"))

      self.menuBar.addmenuitem('Residue', 'command', 'N-Methyl',
                               label='N-Methyl [Alt-Z]',
                               command = lambda: cmd.do(
         "_ editor.attach_amino_acid('pk1','nme')"))

      self.menuBar.addmenuitem('Residue', 'command', 'Phenylalanine',
                               label='Phenylalanine [Alt-F]',
                               command = lambda: cmd.do(
         "_ editor.attach_amino_acid('pk1','phe')"))

      self.menuBar.addmenuitem('Residue', 'command', 'Proline',
                               label='Proline [Alt-P]',
                               command = lambda: cmd.do(
         "_ editor.attach_amino_acid('pk1','pro')"))

      self.menuBar.addmenuitem('Residue', 'command', 'Serine',
                               label='Serine [Alt-S]',
                               command = lambda: cmd.do(
         "_ editor.attach_amino_acid('pk1','ser')"))

      self.menuBar.addmenuitem('Residue', 'command', 'Threonine',
                               label='Threonine [Alt-T]',
                               command = lambda: cmd.do(
         "_ editor.attach_amino_acid('pk1','thr')"))

      self.menuBar.addmenuitem('Residue', 'command', 'Tryptophan',
                               label='Tryptophan [Alt-W]',
                               command = lambda: cmd.do(
         "_ editor.attach_amino_acid('pk1','trp')"))

      self.menuBar.addmenuitem('Residue', 'command', 'Tyrosine',
                               label='Tyrosine [Alt-Y]',
                               command = lambda: cmd.do(
         "_ editor.attach_amino_acid('pk1','tyr')"))

      self.menuBar.addmenuitem('Residue', 'command', 'Valine',
                               label='Valine [Alt-V]',
                               command = lambda: cmd.do(
         "_ editor.attach_amino_acid('pk1','val')"))


      self.menuBar.addmenuitem('Residue', 'separator', '')
      self.menuBar.addmenuitem('Residue', 'command', 'Helix',
                               label='Helix',
                               command = lambda: cmd.do("_ set secondary_structure,1"))

      self.menuBar.addmenuitem('Residue', 'command', 'Antiparallel Beta Sheet',
                               label='Antiparallel Beta Sheet',
                               command = lambda: cmd.do("_ set secondary_structure,2"))

      self.menuBar.addmenuitem('Residue', 'command', 'Parallel Beta Sheet',
                               label='Parallel Beta Sheet',
                               command = lambda: cmd.do("_ set secondary_structure,3"))

      self.menuBar.addmenu('Movie', 'Movie Control')

      
      self.menuBar.addcascademenu('Movie', 'Speed', 'Playback Speed',
                                  label=self.pad+'Speed')

      self.menuBar.addmenuitem('Speed', 'command', 'Maximum',
                               label=self.pad+'Maximum',
                               command = lambda: cmd.set("movie_delay","0",log=1))

      self.menuBar.addmenuitem('Speed', 'command', '30 FPS',
                               label=self.pad+'30 FPS',
                               command = lambda: cmd.set("movie_delay","33",log=1))

      self.menuBar.addmenuitem('Speed', 'command', '15 FPS',
                               label=self.pad+'15 FPS',
                               command = lambda: cmd.set("movie_delay","66",log=1))

      self.menuBar.addmenuitem('Speed', 'command', '5 FPS',
                               label=self.pad+'5 FPS',
                               command = lambda: cmd.set("movie_delay","200",log=1))

      self.menuBar.addmenuitem('Speed', 'command', '1 FPS',
                               label=self.pad+'1 FPS',
                               command = lambda: cmd.set("movie_delay","1000",log=1))

      self.menuBar.addmenuitem('Speed', 'command', '0.3 FPS',
                               label=self.pad+'0.3 FPS',
                               command = lambda: cmd.set("movie_delay","3000",log=1))

      self.menuBar.addmenuitem('Movie', 'command', 'Reset Meter',
                               label=self.pad+'Reset Meter',
                               command = lambda: cmd.do("_ meter_reset"))

      self.menuBar.addmenuitem('Movie', 'separator', '')

      self.menuBar.addmenuitem('Movie', 'checkbutton',
                         'Photorealistic images.',
                         label=self.pad+'Render Frames',
                        variable = self.setting.ray_trace_frames,
                        command = lambda s=self: s.setting.update('ray_trace_frames'))

      self.menuBar.addmenuitem('Movie', 'checkbutton',
                         'Save images in memory.',
                         label=self.pad+'Cache Frames',
                        variable = self.setting.cache_frames,
                        command = lambda s=self: s.setting.update('cache_frames'))

      self.menuBar.addmenuitem('Movie', 'command', 'Flush Cache',
                               label=self.pad+'Flush Cache',
                               command = lambda: cmd.mclear())

      self.menuBar.addmenuitem('Movie', 'separator', '')

      self.menuBar.addmenuitem('Movie', 'checkbutton',
                         'Static Singletons Objects',
                         label=self.pad+'Static Singletons',
                        variable = self.setting.static_singletons,
                        command = lambda s=self: s.setting.update('static_singletons'))

      self.menuBar.addmenuitem('Movie', 'checkbutton',
                         'Superimpose all molecular states.',
                         label=self.pad+'Show All States',
                        variable = self.setting.all_states,
                        command = lambda s=self: s.setting.update('all_states'))


      self.menuBar.addmenu('Rendering', 'Rendering Control')

      self.menuBar.addmenuitem('Rendering', 'checkbutton',
                         'Smooth raytracing.',
                         label=self.pad+'Antialias',
                        variable = self.setting.antialias,
                        command = lambda s=self: s.setting.update('antialias'))

      self.menuBar.addmenuitem('Rendering', 'separator', '')
      
      self.menuBar.addcascademenu('Rendering', 'Shadows', 'Shadows',
                               label=self.pad+'Shadows')

      self.menuBar.addmenuitem('Shadows', 'command', 'None',
                               label='None',
                               command = lambda : cmd.do("_ util.ray_shadows('none')"))

      self.menuBar.addmenuitem('Shadows', 'command', 'Light',
                               label='Light',
                               command = lambda : cmd.do("_ util.ray_shadows('light')"))

      self.menuBar.addmenuitem('Shadows', 'command', 'Matte',
                               label='Matte',
                               command = lambda : cmd.do("_ util.ray_shadows('matte')"))

      self.menuBar.addmenuitem('Shadows', 'command', 'Medium',
                               label='Medium',
                               command = lambda : cmd.do("_ util.ray_shadows('medium')"))

      self.menuBar.addmenuitem('Shadows', 'command', 'Heavy',
                               label='Heavy',
                               command = lambda : cmd.do("_ util.ray_shadows('heavy')"))

      self.menuBar.addmenuitem('Shadows', 'command', 'Black',
                               label='Black',
                               command = lambda : cmd.do("_ util.ray_shadows('black')"))

      self.menuBar.addcascademenu('Rendering', 'Texture', 'Texture',
                               label=self.pad+'Texture')

      self.menuBar.addmenuitem('Texture', 'command', 'None',
                               label='None',
                               command = lambda : cmd.do("_ cmd.set('ray_texture',0)"))

      self.menuBar.addmenuitem('Texture', 'command', 'Matte 1',
                               label='Matte 1',
                               command = lambda : cmd.do("_ cmd.set('ray_texture',1)"))

      self.menuBar.addmenuitem('Texture', 'command', 'Matte 2',
                               label='Matte 2',
                               command = lambda : cmd.do("_ cmd.set('ray_texture',4)"))

      self.menuBar.addmenuitem('Texture', 'command', 'Swirl 1',
                               label='Swirl 1',
                               command = lambda : cmd.do("_ cmd.set('ray_texture',2)"))

      self.menuBar.addmenuitem('Texture', 'command', 'Swirl 2',
                               label='Swirl 2',
                               command = lambda : cmd.do("_ cmd.set('ray_texture',3)"))

      self.menuBar.addmenuitem('Texture', 'command', 'Fiber',
                               label='Fiber',
                               command = lambda : cmd.do("_ cmd.set('ray_texture',5)"))



      self.menuBar.addmenuitem('Rendering', 'separator', '')

      self.menuBar.addmenuitem('Rendering', 'checkbutton',
                         'Cull Backfaces when Rendering',
                         label=self.pad+'Cull Backfaces',
                        variable = self.setting.backface_cull,
                        command = lambda s=self: s.setting.update('backface_cull'))

      self.menuBar.addmenu('Display', 'Display Control')


      
      self.menuBar.addmenuitem('Display', 'command', 'Stereo On',
                               label='Stereo On',
                               command = lambda: cmd.do("_ stereo on"))

      self.menuBar.addmenuitem('Display', 'command', 'Stereo Off',
                               label='Stereo Off',
                               command = lambda: cmd.do("_ stereo off"))

      self.menuBar.addcascademenu('Display', 'Stereo', 'Stereo',
                               label='Stereo')

      self.menuBar.addmenuitem('Stereo', 'command', 'Cross-Eye Stereo',
                               label='Cross-Eye Stereo',
                               command = lambda: cmd.do("_ stereo crosseye"))

      self.menuBar.addmenuitem('Stereo', 'command', 'Quad-Buffered Stereo',
                               label='Quad-Buffered Stereo',
                               command = lambda: cmd.do("_ stereo quadbuffer"))

      self.menuBar.addmenuitem('Stereo', 'separator', '')

      self.menuBar.addmenuitem('Stereo', 'command', 'Swap Sides',
                               label='Swap Sides',
                               command = lambda: cmd.do("_ stereo swap"))

      self.menuBar.addmenuitem('Display', 'separator', '')

      self.menuBar.addcascademenu('Display', 'Performance', 'Performance',
                                  label='Performance')

      self.menuBar.addmenuitem('Performance', 'command', 'Maximum Performance',
                               label='Maximum Performance',
                               command = lambda : cmd.do("_ util.performance(100)"))

      self.menuBar.addmenuitem('Performance', 'command', 'Reasonable Performance',
                               label='Reasonable Performance',
                               command = lambda : cmd.do("_ util.performance(66)"))
      
      self.menuBar.addmenuitem('Performance', 'command', 'Reasonable Quality',
                               label='Reasonable Quality',
                               command = lambda : cmd.do("_ util.performance(33)"))

      self.menuBar.addmenuitem('Performance', 'command', 'Maximum Quality',
                               label='Maximum Quality',
                               command = lambda : cmd.do("_ util.performance(0)"))

      self.menuBar.addcascademenu('Display', 'Background', 'Background',
                                  label='Background')

      self.menuBar.addmenuitem('Background', 'command', 'White Background',
                               label='White',
                               command = lambda : cmd.do("_ cmd.bg_color('white')"))

      self.menuBar.addmenuitem('Background', 'command', 'Light Grey',
                               label='Light Grey',
                               command = lambda : cmd.do("_ cmd.bg_color('grey80')"))

      self.menuBar.addmenuitem('Background', 'command', 'Grey',
                               label='Grey',
                               command = lambda : cmd.do("_ cmd.bg_color('grey50')"))


      self.menuBar.addmenuitem('Background', 'command', 'Black Background',
                               label='Black',
                               command = lambda : cmd.do("_ cmd.bg_color('black')"))

      self.menuBar.addmenuitem('Display', 'separator', '')
      
      self.menuBar.addmenuitem('Display', 'checkbutton',
                         'Disable perspective.',
                         label=self.pad+'Orthoscopic View',
                        variable = self.setting.ortho,
                        command = lambda s=self: s.setting.update('ortho'))

      self.menuBar.addmenuitem('Display', 'checkbutton',
                         'Show Valences.',
                         label=self.pad+'Show Valences',
                        variable = self.setting.valence,
                        command = lambda s=self: s.setting.update('valence'))


      self.menuBar.addmenuitem('Display', 'checkbutton',
                         'Smooth Lines.',
                         label=self.pad+'Smooth Lines',
                        variable = self.setting.line_smooth,
                        command = lambda s=self: s.setting.update('line_smooth'))

      self.menuBar.addmenuitem('Display', 'checkbutton',
                         'Depth Cue Fog.',
                         label=self.pad+'Depth Cue & Fog',
                        variable = self.setting.depth_cue,
                        command = lambda s=self: s.setting.update('depth_cue'))

      self.menuBar.addmenuitem('Display', 'checkbutton',
                         'Two Sided Lighting.',
                         label=self.pad+'Two Sided Lighting',
                        variable = self.setting.two_sided_lighting,
                        command = lambda s=self: s.setting.update('two_sided_lighting'))

      self.menuBar.addmenuitem('Display', 'checkbutton',
                         'Specular Reflections.',
                         label=self.pad+'Specular Reflections',
                        variable = self.setting.specular,
                        command = lambda s=self: s.setting.update('specular'))



      self.menuBar.addmenu('Setting', 'Configuration Control')

      self.menuBar.addmenuitem('Setting', 'command',
                         'Edit PyMOL Settings',
                         label=self.pad+'Edit All...',
                               command = lambda s=self: SetEditor(s))

      self.menuBar.addmenuitem('Setting', 'command',
                         'Edit PyMOL Colors',
                         label=self.pad+'Colors...',
                               command = lambda s=self: ColorEditor(s))


      self.menuBar.addmenuitem('Setting', 'separator', '')
      

      self.menuBar.addmenuitem('Setting', 'checkbutton',
                               'Ignore PDB segi.',
                               label=self.pad+'Ignore PDB Segment Identifier',
                               variable = self.setting.ignore_pdb_segi,
                               command = lambda s=self: s.setting.update('ignore_pdb_segi'))

      self.menuBar.addmenuitem('Setting', 'checkbutton',
                         'Auto-Zoom.',
                         label=self.pad+'Auto-Zoom New Objects',
                        variable = self.setting.auto_zoom,
                        command = lambda s=self: s.setting.update('auto_zoom'))

      self.menuBar.addmenuitem('Setting', 'checkbutton',
                         'Auto-Show Selections.',
                         label=self.pad+'Auto-Show New Selections',
                        variable = self.setting.auto_show_selections,
                        command = lambda s=self: s.setting.update('auto_show_selections'))

      self.menuBar.addmenuitem('Setting', 'checkbutton',
                         'Auto-Hide Selections.',
                         label=self.pad+'Auto-Hide Selections',
                        variable = self.setting.auto_hide_selections,
                        command = lambda s=self: s.setting.update('auto_hide_selections'))

      self.menuBar.addmenuitem('Setting', 'checkbutton',
                         'Auto-Remove Hydrogens.',
                         label=self.pad+'Auto-Remove Hydrogens',
                        variable = self.setting.auto_remove_hydrogens,
                        command = lambda s=self: s.setting.update('auto_remove_hydrogens'))

      self.menuBar.addmenuitem('Setting', 'separator', '')


      self.menuBar.addmenuitem('Setting', 'command', 'Show Text Output',
                               label='Show Text',
                               command = lambda: cmd.set("text","1",log=1))

      self.menuBar.addmenuitem('Setting', 'command', 'Hide Text Output',
                               label='Hide Text',
                               command = lambda: cmd.set("text","0",log=1))

      self.menuBar.addmenuitem('Setting', 'checkbutton',
                         'Overlay Text Output on Graphics',
                         label=self.pad+'Overlay Text',
                        variable = self.setting.overlay,
                        command = lambda s=self: s.setting.update('overlay'))

      
      self.menuBar.addmenu('Mouse', 'Mouse Configuration')

      self.menuBar.addmenuitem('Mouse', 'command', '3 Button Universal Cycle',
                               label='3 Button Universal Cycle',
                               command = lambda: cmd.config_mouse('three_button'))

      self.menuBar.addmenuitem('Mouse', 'command', '2 Button Viewing Cycle',
                               label='2 Button Viewing Cycle',
                               command = lambda: cmd.config_mouse('two_button'))

      self.menuBar.addmenuitem('Mouse', 'command', 'Setup 2 Button Editing Cycle',
                               label='2 Button Editing Cycle',
                               command = lambda: cmd.config_mouse('two_button_editing'))

      self.menuBar.addmenuitem('Mouse', 'separator', '')

      self.menuBar.addmenuitem('Mouse', 'command', '3 Button Viewing Mode',
                               label='3 Button Viewing Mode',
                               command = lambda: cmd.mouse('three_button_viewing'))

      self.menuBar.addmenuitem('Mouse', 'command', '3 Button Editing Mode',
                               label='3 Button Editing Mode',
                               command = lambda: cmd.mouse('three_button_editing'))

      self.menuBar.addmenuitem('Mouse', 'command', '2 Button Viewing Mode',
                               label='2 Button Viewing Mode',
                               command = lambda: cmd.mouse('two_button_viewing'))

      self.menuBar.addmenuitem('Mouse', 'command', '2 Button Selecting Mode',
                               label='2 Button Selecting Mode',
                               command = lambda: cmd.mouse('two_button_selecting'))

      self.menuBar.addmenuitem('Mouse', 'command', '2 Button Editing Mode',
                               label='2 Button Editing Mode',
                               command = lambda: cmd.mouse('two_button_editing'))

      self.menuBar.addmenu('Cartoon', 'Cartoon Properties')

      self.menuBar.addmenuitem('Cartoon', 'checkbutton',
                         'Round Helices',
                         label=self.pad+'Round Helices',
                        variable = self.setting.cartoon_round_helices,
                        command = lambda s=self: s.setting.update('cartoon_round_helices'))

      self.menuBar.addmenuitem('Cartoon', 'checkbutton',
                         'Fancy Helices',
                         label=self.pad+'Fancy Helices',
                        variable = self.setting.cartoon_fancy_helices,
                        command = lambda s=self: s.setting.update('cartoon_fancy_helices'))

      self.menuBar.addmenuitem('Cartoon', 'checkbutton',
                         'Cylindrical Helices',
                         label=self.pad+'Cylindrical Helices',
                        variable = self.setting.cartoon_cylindrical_helices,
                        command = lambda s=self: s.setting.update('cartoon_cylindrical_helices'))

      self.menuBar.addmenuitem('Cartoon', 'checkbutton',
                         'Flat Sheets',
                         label=self.pad+'Flat Sheets',
                        variable = self.setting.cartoon_flat_sheets,
                        command = lambda s=self: s.setting.update('cartoon_flat_sheets'))


      self.menuBar.addmenuitem('Cartoon', 'checkbutton',
                         'Fancy Sheets',
                         label=self.pad+'Fancy Sheets',
                        variable = self.setting.cartoon_fancy_sheets,
                        command = lambda s=self: s.setting.update('cartoon_fancy_sheets'))

      self.menuBar.addmenuitem('Cartoon', 'checkbutton',
                         'Smooth Loops',
                         label=self.pad+'Smooth Loops',
                        variable = self.setting.cartoon_smooth_loops,
                        command = lambda s=self: s.setting.update('cartoon_smooth_loops'))

      self.menuBar.addmenuitem('Cartoon', 'checkbutton',
                         'Discrete Colors',
                         label=self.pad+'Discrete Colors',
                        variable = self.setting.cartoon_discrete_colors,
                        command = lambda s=self: s.setting.update('cartoon_discrete_colors'))


      self.menuBar.addmenu('Wizard', 'Task Wizards')
      
      self.menuBar.addmenuitem('Wizard', 'command', 'Density Map Wizard',
                               label='Density',
                               command = lambda: cmd.do("_ wizard density"))

      self.menuBar.addmenuitem('Wizard', 'command', 'Distance',
                               label='Distance',
                               command = lambda: cmd.do("_ wizard distance"))

      self.menuBar.addmenuitem('Wizard', 'command', 'Filter',
                               label='Filter',
                               command = lambda: cmd.do("_ wizard filter"))

      self.menuBar.addmenuitem('Wizard', 'command', 'Mutagenesis',
                               label='Mutagenesis',
                               command = lambda: cmd.do("_ wizard mutagenesis"))

      self.menuBar.addmenuitem('Wizard', 'command', 'Pair Fitting',
                               label='Pair Fitting',
                               command = lambda: cmd.do("_ wizard pair_fit"))

      self.menuBar.addmenuitem('Wizard', 'command', 'Sculpting',
                               label='Sculpting',
                               command = lambda: cmd.do("_ wizard sculpting"))

      self.menuBar.addmenuitem('Wizard', 'separator', '')
      
      self.menuBar.addmenuitem('Wizard', 'command', 'Label',
                               label='Label',
                               command = lambda: cmd.do("_ wizard label"))

      self.menuBar.addmenuitem('Wizard', 'command', 'Charge',
                               label='Charge',
                               command = lambda: cmd.do("_ wizard charge"))


      self.menuBar.addmenu('Sculpting', 'Sculpting')

      self.menuBar.addmenuitem('Sculpting', 'checkbutton',
                         'Auto-Sculpt.',
                         label=self.pad+'Auto-Sculpting',
                        variable = self.setting.auto_sculpt,
                        command = lambda s=self: s.setting.update('auto_sculpt'))

      self.menuBar.addmenuitem('Sculpting', 'checkbutton',
                         'Sculpting.',
                         label=self.pad+'Sculpting',
                        variable = self.setting.sculpting,
                        command = lambda s=self: s.setting.update('sculpting'))

      self.menuBar.addmenuitem('Sculpting', 'separator', '')
      

      self.menuBar.addmenuitem('Sculpting', 'command', 'Activate',
                               label='Activate',
                               command = lambda: cmd.do("_ sculpt_activate all"))

      self.menuBar.addmenuitem('Sculpting', 'command', 'Deactivate',
                               label='Deactivate',
                               command = lambda: cmd.do("_ sculpt_deactivate all"))

      self.menuBar.addmenuitem('Sculpting', 'command', 'Clear Memory',
                               label='Clear Memory',
                               command = lambda: cmd.do("_ sculpt_purge"))

      self.menuBar.addmenuitem('Sculpting', 'separator', '')

      self.menuBar.addmenuitem('Sculpting', 'command', '1 Cycle/Update',
                               label='1 Cycle per Update',
                               command = lambda: cmd.do("_ set sculpting_cycles=1"))

      self.menuBar.addmenuitem('Sculpting', 'command', '3 Cycles/Update',
                               label='3 Cycles per Update',
                               command = lambda: cmd.do("_ set sculpting_cycles=3"))

      self.menuBar.addmenuitem('Sculpting', 'command', '10 Cycles/Update',
                               label='10 Cycles per Update',
                               command = lambda: cmd.do("_ set sculpting_cycles=10"))

      self.menuBar.addmenuitem('Sculpting', 'command', '33 Cycles/Update',
                               label='33 Cycles per Update',
                               command = lambda: cmd.do("_ set sculpting_cycles=33"))

      self.menuBar.addmenuitem('Sculpting', 'command', '100 Cycles/Update',
                               label='100 Cycles per Update',
                               command = lambda: cmd.do("_ set sculpting_cycles=100"))

      self.menuBar.addmenuitem('Sculpting', 'command', '333 Cycles/Update',
                               label='333 Cycles per Update',
                               command = lambda: cmd.do("_ set sculpting_cycles=333"))

      self.menuBar.addmenuitem('Sculpting', 'command', '1000 Cycles/Update',
                               label='1000 Cycles per Update',
                               command = lambda: cmd.do("_ set sculpting_cycles=1000"))

      self.menuBar.addmenuitem('Sculpting', 'separator', '')

      self.menuBar.addmenuitem('Sculpting', 'command', 'Bonds Only',
                               label='Bonds Only',
                               command = lambda: cmd.do("_ set sculpt_field_mask=1"))

      self.menuBar.addmenuitem('Sculpting', 'command', 'Bonds & Angles Only',
                               label='Bonds & Angles Only',
                               command = lambda: cmd.do("_ set sculpt_field_mask=3"))

      self.menuBar.addmenuitem('Sculpting', 'command', 'All Except VDW',
                               label='All Except VDW',
                               command = lambda: cmd.do("_ set sculpt_field_mask=31"))

      self.menuBar.addmenuitem('Sculpting', 'command', 'All Except 1-4 VDW',
                               label='All Except 1-4 VDW',
                               command = lambda: cmd.do("_ set sculpt_field_mask=63"))

      self.menuBar.addmenuitem('Sculpting', 'command', 'All Terms',
                               label='All Terms',
                               command = lambda: cmd.do("_ set sculpt_field_mask=127"))

      self.menuBar.addmenu('Demo', 'Demonstrations')

      self.demo = Demo()
      
      self.menuBar.addmenuitem('Demo', 'command', 'Representations',
                               label='Representations',
                               command = lambda s=self:s.demo('rep'))

      self.menuBar.addmenuitem('Demo', 'command', 'Cartoon Ribbons',
                               label='Cartoon Ribbons',
                               command = lambda s=self:s.demo('cartoon'))

      self.menuBar.addmenuitem('Demo', 'command', 'Transparency',
                               label='Transparency',
                               command = lambda s=self:s.demo('trans'))

      self.menuBar.addmenuitem('Demo', 'command', 'Ray Tracing',
                               label='Ray Tracing',
                               command = lambda s=self:s.demo('ray'))

      self.menuBar.addmenuitem('Demo', 'command', 'Sculpting',
                               label='Sculpting',
                               command = lambda s=self:s.demo('sculpt'))

      self.menuBar.addmenuitem('Demo', 'command', 'Scripted Animation',
                               label='Scripted Animation',
                               command = lambda s=self:s.demo('anime'))

      self.menuBar.addmenuitem('Demo', 'command', 'Molscript/Raster3D Input',
                               label='Molscript/Raster3D Input',
                               command = lambda s=self:s.demo('raster3d'))

      self.menuBar.addmenuitem('Demo', 'command', 'Compiled Graphics Objects',
                               label='Compiled Graphics Objects',
                               command = lambda s=self:s.demo('cgo'))
