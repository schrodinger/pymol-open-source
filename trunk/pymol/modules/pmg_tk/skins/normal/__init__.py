

import sys, string
import re
import thread
import os

from Tkinter import *
from tkFileDialog import *
import tkMessageBox

import Pmw

from pymol.wizard import cleanup

from pmg_tk.Setting import Setting
from pmg_tk.SetEditor import SetEditor
from pmg_tk.ColorEditor import ColorEditor

from pmg_tk.skins import PMGSkin
from builder import Builder


class Normal(PMGSkin):

    pad = ' ' # extra space in menus
    
    appname        = 'PyMOL Tcl/Tk GUI'
    appversion     = '1.00'
    copyright      = ('Copyright (C) 1998-2006 by Warren DeLano and \n'+
                            'DeLano Scientific LLC. All rights reserved.')
    contactweb     = 'http://www.pymol.org'
    contactemail   = 'warren@delanoscientific.com'
    
    # responsible for setup and takedown of the normal skin

    def complete(self,event):
        st = self.cmd._parser.complete(self.command.get())
        if st:
            self.command.set(st)
            self.entry.icursor(len(st))
        self.focus_entry = 1
        return 1

    def createDataArea(self):
        # Create data area where data entry widgets are placed.
        self.dataArea = self.app.createcomponent('dataarea',
                                             (), None,
                                             Frame, (self.app._hull,), 
                                             relief=SUNKEN, 
                                             bd=1)
        self.dataArea.pack(side=LEFT, fill=BOTH, expand=YES,
                            padx=1, pady=1)

    def createCommandArea(self):
        # Create a command area for application-wide buttons.
        self.commandFrame = self.app.createcomponent('commandframe', (), None,
            Frame,(self.app._hull,),relief=SUNKEN,bd=1)
        self.commandFrame.place(width=500)
        self.commandFrame.pack(side=TOP, 
                         expand=NO, 
                         fill=BOTH,
                         padx=1,
                         pady=1)

    def createMessageBar(self):
        # Create the message bar area for help and status messages.
        frame = self.app.createcomponent('bottomtray', (), None,
                                     Frame,(self.app._hull,), relief=SUNKEN)
        self.__messageBar = self.app.createcomponent('messagebar',
                                                  (), None,
                                                 Pmw.MessageBar, 
                                                 (frame,),
                                                 #entry_width = 40,
                                                 entry_relief=SUNKEN,
                                                 entry_bd=1,
                                                 labelpos=None)
        self.__messageBar.pack(side=LEFT, expand=NO, fill=X)


    def confirm_quit(self,e=None):
        if int(self.cmd.get_setting_legacy("session_changed")):
            session_file = self.cmd.get_setting_text("session_file")
            if session_file != '':
                message = "Save the current session '%s'?"%os.path.split(session_file)[1]
            else:
                message = "Save the current session?"
            check = tkMessageBox._show("Save Session", message,
                                       tkMessageBox.QUESTION, tkMessageBox.YESNOCANCEL)
            if check==tkMessageBox.YES:
                if self.session_save():
                    self.quit_app()
            elif check==tkMessageBox.NO:
                self.quit_app()
        else:
            self.quit_app()

    def quit_app(self):
        self.cmd.log_close()
        self.cmd.quit()  # avoid logging this - it is inconvenient...


    def buttonAdd(self,frame,text,cmd):
        newBtn=Button(frame,
                          text=text,highlightthickness=0,
                          command=cmd,padx=0,pady=0)
        newBtn.pack(side=LEFT,fill=BOTH,expand=YES)

    def get_view(self):
        self.cmd.get_view(quiet=0)
        try:
            str = self.cmd.get_view(3,quiet=1)
            self.root.clipboard_clear()
            self.root.clipboard_append(str)
            self.last_view = str
            self.app.selection_clear()
            self.app.selection_own()
            self.app.selection_handle(lambda a,b,s=self:s.last_view)
            print " PyMOL: Viewing matrix copied to clipboard."
        except:
            traceback.print_exc()
        
    def createButtons(self):
        
        row1 = self.app.createcomponent('row1', (), None,
            Frame,self.commandFrame,bd=0)
        row1.pack(side=TOP,fill=BOTH,expand=YES)
        btn_reset = self.buttonAdd(row1,'Reset',lambda s=self: s.cmd.do("_ reset"))
        btn_reset = self.buttonAdd(row1,'Zoom',lambda s=self: s.cmd.do("_ zoom animate=1"))
        btn_rtrace = self.buttonAdd(row1,'Draw',lambda s=self: s.cmd.do("_ draw"))        
        btn_rtrace = self.buttonAdd(row1,'Ray',lambda s=self: s.cmd.do("_ ray async=1"))
        btn_reset = self.buttonAdd(row1,'Rock',lambda s=self: s.cmd.do("_ rock"))

        row2 = self.app.createcomponent('row2', (), None,
            Frame,self.commandFrame,bd=0)
        row2.pack(side=TOP,fill=BOTH,expand=YES)
        btn_unpick = self.buttonAdd(row2,'Unpick',lambda s=self: s.cmd.do("_ unpick"))
        btn_hidesele = self.buttonAdd(row2,'Deselect',self.hide_sele)
        btn_getview = self.buttonAdd(row2,'Get View',lambda s=self: s.get_view()) # doesn't get logged


        row3 = self.app.createcomponent('row3', (), None,
            Frame,self.commandFrame,bd=0)
        row3.pack(side=TOP,fill=BOTH,expand=YES)
        btn_rewind = self.buttonAdd(row3,'|<',lambda s=self: s.cmd.do("_ rewind"))
        btn_back = self.buttonAdd(row3,'<',lambda s=self: s.cmd.do("_ backward"))
        btn_stop = self.buttonAdd(row3,'Stop',lambda s=self: s.cmd.do("_ mstop"))
        btn_play = self.buttonAdd(row3,'Play',lambda s=self: s.cmd.do("_ mplay"))
        btn_forward = self.buttonAdd(row3,'>',lambda s=self: s.cmd.do("_ forward"))
        btn_last = self.buttonAdd(row3,'>|',lambda s=self: s.cmd.do("_ ending"))
        btn_ccache = self.buttonAdd(row3,'MClear',lambda s=self: s.cmd.do("_ mclear"))

        row4 = self.app.createcomponent('row4', (), None,
            Frame,self.commandFrame,bd=0)
        row4.pack(side=TOP,fill=BOTH,expand=YES)
        self.cmdB = self.buttonAdd(row4,'Command',
                                            lambda s=self:
                                            s.toggleFrame(s.cmdFrame))
        self.buildB = self.buttonAdd(row4,'Builder',
                                              lambda s=self:
                                              s.toggleFrame(s.buildFrame))
        
        
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
                y = max(0,(vh-th)/2)
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

    def back(self):
        if not self.history_cur:
            self.history[0] = self.command.get()
        self.history_cur = (self.history_cur + 1) & self.history_mask
        self.command.set(self.history[self.history_cur])
        l = len(self.history[self.history_cur])
        self.entry.icursor(l)
    
    def forward(self):
        if not self.history_cur:
            self.history[0] = self.command.get()
        self.history_cur = (self.history_cur - 1) & self.history_mask
        self.command.set(self.history[self.history_cur])
        l = len(self.history[self.history_cur])
        self.entry.icursor(l)
        
    def do(self,cmmd):
        self.history[0]=cmmd
        self.history.insert(0,'') # always leave blank at 0
        self.history.pop(self.history_mask+1)
        self.history_cur = 0
        self.cmd.do(cmmd)
                          
    def createConsole(self):
        self.command = StringVar()      
        self.lineCount = 0
        self.history_mask = 0xFF
        self.history = [''] * (self.history_mask+1)
        self.history_cur = 0

        self.cmdFrame = Frame(self.dataArea)
        self.buildFrame = Builder(self.dataArea)
        
        self.toggleFrame(self.cmdFrame,startup=1)

        self.entry = Entry(self.cmdFrame, justify=LEFT, width=50,
             textvariable=self.command)
        self.entry.pack(side=BOTTOM,expand=NO,fill=X)
        self.output = Pmw.ScrolledText(self.cmdFrame)
        self.output.pack(side=TOP, fill=BOTH, expand=YES)      

        self.entry.bind('<Return>', lambda e, s=self:
             (s.do(s.command.get()), s.cmd.dirty(), s.command.set('')))
        self.entry.bind('<Tab>', lambda e, s=self: s.complete(e))
        self.entry.bind('<Up>', lambda e, s=self: s.back())
        self.entry.bind('<Down>', lambda e, s=self: s.forward())
        self.root.protocol("WM_DELETE_WINDOW", lambda s=self: s.confirm_quit())
        
        self.initialdir = os.getcwd()
        self.log_file = "log.pml"      

#      self.entry = self.app.createcomponent('entry', (), None,
#                           Entry,
#                           (self.dataArea,),
#                           justify=LEFT,
#                           width=50,
###                           textvariable=self.command)

        text = self.output.component('text')
        self.text = text      
        if sys.platform[:5]=='linux':
            text.tk.call('tk','scaling',1)
            self.font = 'fixed' # should be available on any X11-based platform
            self.my_fw_font=(self.font,10)
        elif sys.platform[:3]=='win': 
            self.font = 'lucida console' # only available on windows
            self.my_fw_font=(self.font,8) 
        else:
            text.tk.call('tk','scaling',1)
            self.font = 'fixed' # should be available on any X11-based platform
            self.my_fw_font=(self.font,10)
                                                                                                         
        text.configure(font = self.my_fw_font)
        text.configure(width=72)


        self.balloon.bind(self.entry, 'Command Input Area')
        
        self.focus_entry=0
        self.refocus_entry=0
        if self.app.allow_after:
            self.output.after(1000,self.update_feedback)
            self.output.after(1000,self.update_menus)
            
        self.output.pack(side=BOTTOM,expand=YES,fill=BOTH)
        self.app.bind(self.entry, 'Command Input Area')

        self.app.bind_all('<F1>',lambda e,s=self: s.cmd.do("cmd._special(1,0,0)"))
        self.app.bind_all('<F2>',lambda e,s=self: s.cmd.do("cmd._special(2,0,0)"))
        self.app.bind_all('<F3>',lambda e,s=self: s.cmd.do("cmd._special(3,0,0)"))
        self.app.bind_all('<F4>',lambda e,s=self: s.cmd.do("cmd._special(4,0,0)"))
        self.app.bind_all('<F5>',lambda e,s=self: s.cmd.do("cmd._special(5,0,0)"))
        self.app.bind_all('<F6>',lambda e,s=self: s.cmd.do("cmd._special(6,0,0)"))
        self.app.bind_all('<F7>',lambda e,s=self: s.cmd.do("cmd._special(7,0,0)"))
        self.app.bind_all('<F8>',lambda e,s=self: s.cmd.do("cmd._special(8,0,0)"))
        self.app.bind_all('<F9>',lambda e,s=self: s.cmd.do("cmd._special(9,0,0)"))
        self.app.bind_all('<F10>',lambda e,s=self: s.cmd.do("cmd._special(10,0,0)"))
        self.app.bind_all('<F11>',lambda e,s=self: s.cmd.do("cmd._special(11,0,0)"))
        self.app.bind_all('<F12>',lambda e,s=self: s.cmd.do("cmd._special(12,0,0)"))
        self.entry.bind('<Prior>',lambda e,s=self: s.cmd.do("cmd._special(104,0,0)"))
        self.entry.bind('<Next>',lambda e,s=self: s.cmd.do("cmd._special(105,0,0)"))
        self.entry.bind('<Home>',lambda e,s=self: s.cmd.do("cmd._special(106,0,0)"))
        self.entry.bind('<End>',lambda e,s=self: s.cmd.do("cmd._special(107,0,0)"))
        if sys.platform=='darwin':
            if self.app.pymol.invocation.options.external_gui==3:
                self.root.bind_all('<Leave>',lambda e,s=self: s.left(e)) # for MacPyMOLX11Hybrid
                self.root.bind_all('<Enter>',lambda e,s=self: s.focus_in(e))

    def focus_in(self,event):
        if self.refocus_entry:
            self.cmd.do("_ cmd.window('defocus')")
            self.refocus_entry = 0
            self.entry.focus_set()

    def left(self,event):
        if id(event.widget) == id(self.root):
            if ((event.y>event.widget.winfo_height())):
                self.cmd.do("_ cmd.window('focus')")
                self.root.focus_set()
                self.refocus_entry = 1

    def update_feedback(self):
        if self.focus_entry:
            self.focus_entry=0
            self.entry.focus_set()
        for a in self.cmd.get_feedback():
            self.output.insert(END,"\n")
            self.output.insert(END,a)
            self.output.see(END)
            self.lineCount = self.lineCount + 1
            if self.lineCount > 10000:
                self.output.delete('0.0','%i.%i' % (self.lineCount-5000,0))
                self.lineCount=5000
#            self.entry.focus_set()
        self.updating = 1
        progress = self.cmd.get_progress()
        if progress>=0.0:
            self.messageBar.message("busy","Progress %d%%..."%int(progress*100))
        else:
            self.messageBar.resetmessages("busy")
        if self.app.allow_after:
            self.output.after(100,self.update_feedback) # 10X a second

    def toggleFrame(self, frame, startup=0):
        if frame not in self.dataArea.slaves():
            # clear all frames in dataArea
            for f in self.dataArea.slaves():
                f.pack_forget()
            # add requested frame to data area
            frame.pack(side=BOTTOM, fill=BOTH, expand=YES)
        else:
            # clear frame from data area
            if frame != self.cmdFrame:
                frame.pack_forget()
                # next command will cause command frame to be turned on if
                # nothing else is visible... might not want this behavior
                self.cmdFrame.pack(side=BOTTOM, fill=BOTH, expand=YES)
        if not startup:
            if frame==self.cmdFrame:
                self.cmd.edit_mode(0)
            elif frame==self.buildFrame:
                frame.deferred_activate()
                self.cmd.edit_mode(1)
                self.cmd.set("valence","1")
            
    def update_menus(self):
        self.setting.refresh()
        if self.app.allow_after:
            self.output.after(500,self.update_menus) # twice a second

    def file_open(self,tutorial=0):
        if not tutorial:
            initdir = self.initialdir
            ftypes =  [("All Readable","*.pdb"),
                       ("All Readable","*.ccp4"),
                       ("All Readable","*.xplor"),
                       ("All Readable","*.mol"),
                       ("All Readable","*.mol2"),
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
                       ("All Readable","*.phi"),
                       ("All Readable","*.fld"),
                       ("All Readable","*.grd"),
                       ("All Readable","*.o"),
                       ("All Readable","*.omap"),                                         
                       ("All Readable","*.brix"),
                       ("All Readable","*.dx"),
                       ("All Readable","*.pqr"),
                       ("All Readable","*.p5m"),
                       ("All Readable","*.p1m"),
                       ("All Readable","*.cube"),
                       ("All Readable","*.moe"), # proprietary format
                       ("PDB File","*.pdb"),
                       ("All Files","*.*"),
                       ("All Files","*"),                                         
                       ("PDB File","*.ent"),
                       ("PyMOL Session","*.pse"),
                       ("CCP4 Map","*.ccp4"),                                         
                       ("XPLOR Map","*.xplor"),
                       ("MOL2/Multi-MOL2","*.mol2"),
                       ("Macromodel File","*.dat"),
                       ("Macromodel File","*.out"),
                       ("Macromodel File","*.mmd"),
                       ("Macromodel File","*.mmod"),
                       ("BRIX/O Map","*.o"),
                       ("BRIX/O Map","*.omap"),
                       ("BRIX/O Map","*.brix"),
                       ("Gaussian Cube Map","*.cube"),
                       ("DX Map","*.dx"),                                         
                       ("AVS (MEAD) Field","*.fld"),                                         
                       ("MOL File","*.mol"),
                       ("MOE File","*.moe"), # proprietary format
                       ("ChemPy Model","*.pkl"),
                       ("Raster3D Scene","*.r3d"),
                       ("SDF File","*.sdf"),
                       ("ChemDraw3D File","*.cc1"),
                       ("ChemDraw3D File","*.cc2"),
                       ("Tinker XYZ File","*.xyz")
                       ]
        else:
            initdir = os.environ['TUT']
            # only list file extensions that are used for tutorial data
            ftypes = [("Tutorial Data","*.pdb"),]
        ofile = askopenfilename(initialdir = initdir,
                                        filetypes=ftypes)
        if len(ofile):
            if not tutorial:
                self.initialdir = re.sub(r"[^\/\\]*$","",ofile)
            try:
                self.cmd.log("load %s\n"%ofile,"cmd.load('%s',quiet=0)\n"%ofile)
                if (string.lower(ofile[-4:])=='.pse') and (ofile!=self.save_file):
                    self.save_file = '' # remove ambiguous default
                self.cmd.load(ofile,quiet=0)
            except self.pymol.CmdException:
                print "Error: unable to open file '%s'"%ofile

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
            self.cmd.log_open(sfile)

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
#            os.chdir(self.initialdir)                 
            self.cmd.resume(ofile)

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
#            os.chdir(self.initialdir)                 
            self.cmd.log_open(ofile,'a')

    def session_save(self):
        self.save_file = self.cmd.get_setting_text("session_file")
        if self.save_file!='':
            self.cmd.log("save %s,format=pse\n"%(self.save_file),
                      "cmd.save('%s',format='pse')\n"%(self.save_file))
            self.cmd.save(self.save_file,"","pse",quiet=0)
            return 0
        else:
            return self.session_save_as()

    def session_save_as(self):
        (self.initialdir, self.save_file) = os.path.split(self.cmd.get_setting_text("session_file"))
        sfile = asksaveasfilename(initialfile = self.save_file,
                                          initialdir = self.initialdir,
                                          filetypes=[
            ("PyMOL Session File","*.pse"),
            ("PyMOL Show File","*.psw"),
            ])
        if len(sfile):
            if re.search(r"\.pse$|\.PSE$|\.psw$|\.PSW$",sfile)==None:
                sfile=sfile+".pse"
            self.initialdir = re.sub(r"[^\/\\]*$","",sfile)
            self.cmd.log("save %s,format=pse\n"%(sfile),
                      "cmd.save('%s',format='pse')\n"%(sfile))
            self.cmd.save(sfile,"",format='pse',quiet=0)
            self.save_file = sfile
            self.cmd.set("session_file",self.save_file,quiet=1)
            return 1
        else:
            return 0
    
    def file_save(self):
        lst = self.cmd.get_names('all')
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
                sfile = sels[0] +".pdb"
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
                        self.cmd.log("save %s,(%s)\n"%(sfile,sels[0]),
                                  "cmd.save('%s','(%s)')\n"%(sfile,sels[0]))
                        self.cmd.save(sfile,"(%s)"%sels[0],quiet=0)

    def hide_sele(self):
        self.cmd.log("util.hide_sele()\n","util.hide_sele()\n")
        self.util.hide_sele()
            
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
            self.__script__ = ofile
#            os.chdir(dir)        
            if re.search("\.pym*$|\.PYM*$",ofile):
                self.cmd.do("run "+ofile);      
            else:
                self.cmd.do("@"+ofile);

    def file_save_png(self):
        sfile = asksaveasfilename(initialdir = self.initialdir,
                 filetypes=[("PNG File","*.png")])
        if len(sfile):
            self.initialdir = re.sub(r"[^\/\\]*$","",sfile)
            self.cmd.log("png %s\n"%sfile,"cmd.png('%s')\n"%sfile)
            self.cmd.png(sfile,quiet=0)

    def file_save_wrl(self):
        sfile = asksaveasfilename(initialdir = self.initialdir,
                 filetypes=[("VRML 2 WRL File","*.wrl")])
        if len(sfile):
            self.initialdir = re.sub(r"[^\/\\]*$","",sfile)
            self.cmd.log("save %s\n"%sfile,"cmd.save('%s')\n"%sfile)
            self.cmd.save(sfile,quiet=0)
            
    def file_save_pov(self):
        sfile = asksaveasfilename(initialdir = self.initialdir,
                 filetypes=[("POV File","*.pov")])
        if len(sfile):
            self.initialdir = re.sub(r"[^\/\\]*$","",sfile)
            self.cmd.log("save %s\n"%sfile,"cmd.save('%s')\n"%sfile)
            self.cmd.save(sfile,quiet=0)
        
    def file_savemovie(self):
        sfile = asksaveasfilename(filetypes=[("Numbered PNG Files","*.png")])
        if len(sfile):
            self.initialdir = re.sub(r"[^\/\\]*$","",sfile)
            self.cmd.log("mpng %s\n"%sfile,"cmd.mpng('%s')\n"%sfile)         
            self.cmd.mpng(sfile)

    def about_plugins(self):
        about = Pmw.MessageDialog((self.app._hull),
                                          title = 'About Plugins',
                                          message_text =
     'Plugins are external modules which extend PyMOL\'s capabilities.\n\n Available plugins (if any) are shown in the Plugin menu.\n\nIf no plugins are listed, then either none have been installed, \nor those that are installed are not yet functional.')
        about.activate(geometry='centerscreenfirst')      

    def transparency_menu(self,name,label,setting_name):
        
        self.menuBar.addcascademenu('Transparency', name, label, label=label)
        
        for lab, val in [ ('Off', 0.0), ('20%', 0.2), ('40%', 0.4), 
                                ('50%', 0.5), ('60%', 0.6), ('80%', 0.8) ]:
            self.menuBar.addmenuitem(name,  'command', lab, label=lab,
                                             command = lambda v=val,s=self,sn=setting_name: s.cmd.set(sn, v))

    def createMenuBar(self):
        self.menuBar = Pmw.MenuBar(self.root, balloon=self.balloon,
                                            hull_relief=RAISED, hull_borderwidth=1) 
        self.menuBar.pack(fill=X)

        self.menuBar.addmenu('Tutorial', 'Tutorial', side='right')      

        self.menuBar.addmenuitem('Tutorial', 'command', 'Open tutorial data file.',
                                label='Open File...',
                                command=lambda s=self: s.file_open(tutorial=1))

# to come
#        self.menuBar.addmenuitem('Tutorial', 'separator', '')
#
#        self.menuBar.addmenuitem('Tutorial', 'command', 'Beginners',
#                                         label='Beginners',
#                                         command = lambda s=self: None)

        self.menuBar.addmenu('Help', 'About %s' % self.appname, side='right')      
        self.menuBar.addmenuitem('Help', 'command',
                                         'Get information on application', 
                                         label='About', command = lambda s=self: s.cmd.do("_ splash"))

        self.menuBar.addmenuitem('Help', 'separator', '')
        
        self.menuBar.addmenuitem('Help', 'command', 'Demo',
                                         label='Demo',
                                         command = lambda s=self: s.cmd.do(
            "_ replace_wizard demo,cartoon"))


#      self.menuBar.addmenuitem('Help', 'command', 'Release Notes',
#                               label='Release Notes',
#                               command = lambda s=self: s.cmd.do("_ cmd.show_help('release')"))

#      self.menuBar.addmenuitem('Help', 'separator', '')
        
#      self.menuBar.addmenuitem('Help', 'command', 'Help on Commands',
#                               label='Commands',
#                               command = lambda s=self: s.cmd.do("_ cmd.show_help('commands')"))

#      self.menuBar.addmenuitem('Help', 'command', 'Help on Launching',
#                               label='Launching',
#                               command = lambda s=self: s.cmd.do("_ cmd.show_help('launching')"))      

#      self.menuBar.addmenuitem('Help', 'separator', '')

#      self.menuBar.addmenuitem('Help', 'command', 'Help on Selections',
#                               label='Select Command',
#                               command = lambda s=self: s.cmd.do("_ cmd.show_help('select')"))      

#      self.menuBar.addmenuitem('Help', 'command', 'Help on Selections',
#                               label='Selection Syntax',
#                               command = lambda s=self: s.cmd.do("_ cmd.show_help('selections')"))      

#      self.menuBar.addmenuitem('Help', 'command', 'Example Selections',
#                               label='Selection Examples',
#                               command = lambda s=self: s.cmd.do("_ cmd.show_help('examples')"))      

#      self.menuBar.addmenuitem('Help', 'separator', '')
        

#      self.menuBar.addmenuitem('Help', 'command', 'Help on the Mouse',
#                               label='Mouse',
#                               command = lambda s=self: s.cmd.do("_ cmd.show_help('mouse')"))      

#      self.menuBar.addmenuitem('Help', 'command', 'Help on the Keyboard',
#                               label='Keyboard',
#                               command = lambda s=self: s.cmd.do("_ cmd.show_help('keyboard')"))      

#      self.menuBar.addmenuitem('Help', 'command', 'Help on Molecular Editing',
#                               label='Molecular Editing',
#                               command = lambda s=self: s.cmd.do("_ cmd.show_help('editing')"))      

#      self.menuBar.addmenuitem('Help', 'command', 'Help on Molecular Editing',
#                               label='Molecular Editing Keys',
#                               command = lambda s=self: s.cmd.do("_ cmd.show_help('edit_keys')"))      

#      self.menuBar.addmenuitem('Help', 'command', 'Help on Stereo',
#                               label='Stereo',
#                               command = lambda s=self: s.cmd.do("_ cmd.show_help('stereo')"))      

#      self.menuBar.addmenuitem('Help', 'separator', '')
        

#      self.menuBar.addmenuitem('Help', 'command', 'Help on the API',
#                               label='API',
#                               command = lambda s=self: s.cmd.do("_ cmd.show_help('api')"))      

        self.toggleBalloonVar = IntVar()
        self.toggleBalloonVar.set(0)
        self.setting = Setting()

#      self.menuBar.addmenuitem('Help', 'separator', '')
        
#      self.menuBar.addmenuitem('Help', 'checkbutton',
#                         'Toggle balloon help',
#                         label='Balloon help',
#                        variable = self.toggleBalloonVar,
#                        command=self.toggleBalloon)

        self.menuBar.addmenu('File', 'File Input',tearoff=TRUE)

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

        self.menuBar.addcascademenu('File', 'SaveImageAs', 'Save Image As',
                                             label=self.pad+'Save Image As',tearoff=FALSE)

        self.menuBar.addmenuitem('SaveImageAs', 'command', 'Save current image as PNG Image.',
                                label='PNG...',
                                command=self.file_save_png)

        self.menuBar.addmenuitem('SaveImageAs', 'separator', '')
        
        self.menuBar.addmenuitem('SaveImageAs', 'command', 'Save current image as VRML.',
                                label='VRML 2...',
                                command=self.file_save_wrl)
        
        self.menuBar.addmenuitem('SaveImageAs', 'command', 'Save current image as PovRay input.',
                                label='POV-Ray...',
                                command=self.file_save_pov)

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
                                command=self.cmd.log_close)

        self.menuBar.addmenuitem('File', 'command', 'Run program or script.',
                                label=self.pad+'Run...',
                                command=self.file_run)


        self.menuBar.addmenuitem('File', 'separator', '')

        self.menuBar.addmenuitem('File', 'command', 'Quit PyMOL',
                                label=self.pad+'Quit',
                                command=self.confirm_quit)

        self.menuBar.addmenuitem('File', 'command', 'Reinitialize PyMOL',
                                label=self.pad+'Reinitialize',
                                command=self.cmd.reinitialize)

#      self.menuBar.addmenuitem('File', 'separator', '')
        
#      self.menuBar.addmenuitem('File', 'checkbutton',
#                         'Log Conformations.',
#                         label=self.pad+'Log Conformations',
#                        variable = self.setting.log_conformations,
#                        command = lambda s=self: s.setting.update('log_conformations'))

#      self.menuBar.addmenuitem('File', 'checkbutton',
#                         'Log Box Selections.',
#                         label=self.pad+'Log Box Selections',
#                        variable = self.setting.log_box_selections,
#                        command = lambda s=self: s.setting.update('log_box_selections'))

        self.menuBar.addmenu('Edit', 'Text Editing',tearoff=TRUE)
        self.menuBar.addmenuitem('Edit', 'command',
                                 'To Copy: Use Ctrl-C in TclTk GUI',
                                 label='To copy text use Ctrl-C in the TclTk GUI',
                                         state='disabled',
                                command =  None)

        self.menuBar.addmenuitem('Edit', 'command',
                                 'To Paste, Use Ctrl-V in TclTk GUI',
                                 label='To paste text use Ctrl-V in the TckTk GUI',
                                         state='disabled',                               
                                command =  None)

        self.menuBar.addmenuitem('Edit', 'separator', '')

        self.menuBar.addmenuitem('Edit', 'command', 'Undo Conformation',
                                         label='Undo Conformation [Ctrl-Z]',
                                         command = lambda s=self: s.cmd.do("_ undo"))

        self.menuBar.addmenuitem('Edit', 'command', 'Redo Conformation',
                                         label='Redo Conformation [Ctrl-A]',
                                         command = lambda s=self: s.cmd.do("_ redo"))

        self.menuBar.addmenu('Build', 'Build',tearoff=TRUE)

        self.menuBar.addcascademenu('Build', 'Fragment', 'Fragment',
                                             label='Fragment',tearoff=TRUE)
        
#      self.menuBar.addmenu('Fragment', 'Fragment')

        self.menuBar.addmenuitem('Fragment', 'command', 'Acetylene',
                                         label='Acetylene [Alt-J]',
                                         command = lambda s=self: s.cmd.do(
            "_ editor.attach_fragment('pk1','acetylene',2,0)"))

        self.menuBar.addmenuitem('Fragment', 'command', 'Amide N->C',
                                         label='Amide N->C [Alt-1]',
                                         command = lambda s=self: s.cmd.do(
            "_ editor.attach_fragment('pk1','formamide',3,1)"))

        self.menuBar.addmenuitem('Fragment', 'command', 'Amide C->N',
                                         label='Amide C->N [Alt-2]',
                                         command = lambda s=self: s.cmd.do(
            "_ editor.attach_fragment('pk1','formamide',5,0)"))

        self.menuBar.addmenuitem('Fragment', 'command', 'Bromine',
                                         label='Bromine [Ctrl-B]',
                                         command = lambda s=self: s.cmd.do("_ replace Br,1,1"))

        self.menuBar.addmenuitem('Fragment', 'command', 'Carbon',
                                         label='Carbon [Ctrl-C]',
                                         command = lambda s=self: s.cmd.do("_ replace C,4,4"))

        self.menuBar.addmenuitem('Fragment', 'command', 'Carbonyl',
                                         label='Carbonyl [Alt-0]',
                                         command = lambda s=self: s.cmd.do(
            "_ editor.attach_fragment('pk1','formaldehyde',2,0)"))

        self.menuBar.addmenuitem('Fragment', 'command', 'Chlorine',
                                         label='Chlorine [Ctrl-L]',
                                         command = lambda s=self: s.cmd.do("_ replace Cl,1,1"))

        self.menuBar.addmenuitem('Fragment', 'command', 'Cyclobutyl',
                                         label='Cyclobutyl [Alt-4]',
                                         command = lambda s=self: s.cmd.do(
            "_ editor.attach_fragment('pk1','cyclobutane',4,0)"))

        self.menuBar.addmenuitem('Fragment', 'command', 'Cyclopentyl',
                                         label='Cyclopentyl [Alt-5]',
                                         command = lambda s=self: s.cmd.do(
            "_ editor.attach_fragment('pk1','cyclopentane',5,0)"))

        self.menuBar.addmenuitem('Fragment', 'command', 'Cyclopentadiene',
                                         label='Cyclopentadiene [Alt-8]',
                                         command = lambda s=self: s.cmd.do(
              "_ editor.attach_fragment('pk1','cyclopentadiene',5,0)"))

        self.menuBar.addmenuitem('Fragment', 'command', 'Cyclohexyl',
                                         label='Cyclohexyl [Alt-6]',
                                         command = lambda s=self: s.cmd.do(
            "_ editor.attach_fragment('pk1','cyclohexane',7,0)"))

        self.menuBar.addmenuitem('Fragment', 'command', 'Cycloheptyl',
                                         label='Cycloheptyl [Alt-7]',
                                         command = lambda s=self: s.cmd.do(
            "_ editor.attach_fragment('pk1','cycloheptane',8,0)"))

        self.menuBar.addmenuitem('Fragment', 'command', 'Fluorine',
                                         label='Fluorine [Ctrl-F]',
                                         command = lambda s=self: s.cmd.do("_ replace F,1,1"))

        self.menuBar.addmenuitem('Fragment', 'command', 'Iodine',
                                         label='Iodine [Ctrl-I]',
                                         command = lambda s=self: s.cmd.do("_ replace I,1,1"))

        self.menuBar.addmenuitem('Fragment', 'command', 'Methane',
                                         label='Methane',
                                         command = lambda s=self: s.cmd.do(
            "_ editor.attach_fragment('pk1','methane',1,0)"))

        self.menuBar.addmenuitem('Fragment', 'command', 'Nitrogen',
                                         label='Nitrogen [Ctrl-N]',
                                         command = lambda s=self: s.cmd.do("_ replace N,4,3"))

        self.menuBar.addmenuitem('Fragment', 'command', 'Oxygen',
                                         label='Oxygen [Ctrl-O]',
                                         command = lambda s=self: s.cmd.do("_ replace O,4,2"))

        self.menuBar.addmenuitem('Fragment', 'command', 'Phenyl',
                                         label='Phenyl [Alt-9]',
                                         command = lambda s=self: s.cmd.do(
            "_ editor.attach_fragment('pk1','benzene',6,0)"))

        self.menuBar.addmenuitem('Fragment', 'command', 'Sulfer',
                                         label='Sulfer [Ctrl-S]',
                                         command = lambda s=self: s.cmd.do("_ replace S,2,2"))

        self.menuBar.addmenuitem('Fragment', 'command', 'Sulfonyl',
                                         label='Sulfonyl [Alt-3]',
                                         command = lambda s=self: s.cmd.do(
            "_ editor.attach_fragment('pk1','sulfone',3,1)"))

        self.menuBar.addmenuitem('Fragment', 'command', 'Phosphorus',
                                         label='Phosphorus [Ctrl-P]',
                                         command = lambda s=self: s.cmd.do("_ replace P,4,3"))

#      self.menuBar.addmenu('Residue', 'Residue')

        self.menuBar.addcascademenu('Build', 'Residue', 'Residue',
                                             label='Residue',tearoff=TRUE)

 
        self.menuBar.addmenuitem('Residue', 'command', 'Acetyl',
                                         label='Acetyl [Alt-B]',
                                         command = lambda s=self: s.cmd.do(
            "_ editor.attach_amino_acid('pk1','ace')"))

        self.menuBar.addmenuitem('Residue', 'command', 'Alanine',
                                         label='Alanine [Alt-A]',
                                         command = lambda s=self: s.cmd.do(
            "_ editor.attach_amino_acid('pk1','ala')"))

        self.menuBar.addmenuitem('Residue', 'command', 'Amine',
                                         label='Amine',
                                         command = lambda s=self: s.cmd.do(
            "_ editor.attach_amino_acid('pk1','nhh')"))

        self.menuBar.addmenuitem('Residue', 'command', 'Aspartate',
                                         label='Aspartate [Alt-D]',
                                         command = lambda s=self: s.cmd.do(
            "_ editor.attach_amino_acid('pk1','asp')"))

        self.menuBar.addmenuitem('Residue', 'command', 'Asparagine',
                                         label='Asparagine [Alt-N]',
                                         command = lambda s=self: s.cmd.do(
            "_ editor.attach_amino_acid('pk1','asn')"))

        self.menuBar.addmenuitem('Residue', 'command', 'Arginine',
                                         label='Arginine [Alt-R]',
                                         command = lambda s=self: s.cmd.do(
            "_ editor.attach_amino_acid('pk1','arg')"))

        self.menuBar.addmenuitem('Residue', 'command', 'Cysteine',
                                         label='Cysteine [Alt-C]',
                                         command = lambda s=self: s.cmd.do(
            "_ editor.attach_amino_acid('pk1','cys')"))

        self.menuBar.addmenuitem('Residue', 'command', 'Glutamate',
                                         label='Glutamate [Alt-E]',
                                         command = lambda s=self: s.cmd.do(
            "_ editor.attach_amino_acid('pk1','glu')"))

        self.menuBar.addmenuitem('Residue', 'command', 'Glutamine',
                                         label='Glutamine [Alt-N]',
                                         command = lambda s=self: s.cmd.do(
            "_ editor.attach_amino_acid('pk1','gln')"))

        self.menuBar.addmenuitem('Residue', 'command', 'Glycine',
                                         label='Glycine [Alt-G]',
                                         command = lambda s=self: s.cmd.do(
            "_ editor.attach_amino_acid('pk1','gly')"))

        self.menuBar.addmenuitem('Residue', 'command', 'Histidine',
                                         label='Histidine [Alt-H]',
                                         command = lambda s=self: s.cmd.do(
            "_ editor.attach_amino_acid('pk1','his')"))

        self.menuBar.addmenuitem('Residue', 'command', 'Isoleucine',
                                         label='Isoleucine [Alt-I]',
                                         command = lambda s=self: s.cmd.do(
            "_ editor.attach_amino_acid('pk1','ile')"))

        self.menuBar.addmenuitem('Residue', 'command', 'Leucine',
                                         label='Leucine [Alt-L]',
                                         command = lambda s=self: s.cmd.do(
            "_ editor.attach_amino_acid('pk1','leu')"))

        self.menuBar.addmenuitem('Residue', 'command', 'Lysine',
                                         label='Lysine [Alt-K]',
                                         command = lambda s=self: s.cmd.do(
            "_ editor.attach_amino_acid('pk1','lys')"))

        self.menuBar.addmenuitem('Residue', 'command', 'Methionine',
                                         label='Methionine [Alt-M]',
                                         command = lambda s=self: s.cmd.do(
            "_ editor.attach_amino_acid('pk1','met')"))

        self.menuBar.addmenuitem('Residue', 'command', 'N-Methyl',
                                         label='N-Methyl [Alt-Z]',
                                         command = lambda s=self: s.cmd.do(
            "_ editor.attach_amino_acid('pk1','nme')"))

        self.menuBar.addmenuitem('Residue', 'command', 'Phenylalanine',
                                         label='Phenylalanine [Alt-F]',
                                         command = lambda s=self: s.cmd.do(
            "_ editor.attach_amino_acid('pk1','phe')"))

        self.menuBar.addmenuitem('Residue', 'command', 'Proline',
                                         label='Proline [Alt-P]',
                                         command = lambda s=self: s.cmd.do(
            "_ editor.attach_amino_acid('pk1','pro')"))

        self.menuBar.addmenuitem('Residue', 'command', 'Serine',
                                         label='Serine [Alt-S]',
                                         command = lambda s=self: s.cmd.do(
            "_ editor.attach_amino_acid('pk1','ser')"))

        self.menuBar.addmenuitem('Residue', 'command', 'Threonine',
                                         label='Threonine [Alt-T]',
                                         command = lambda s=self: s.cmd.do(
            "_ editor.attach_amino_acid('pk1','thr')"))

        self.menuBar.addmenuitem('Residue', 'command', 'Tryptophan',
                                         label='Tryptophan [Alt-W]',
                                         command = lambda s=self: s.cmd.do(
            "_ editor.attach_amino_acid('pk1','trp')"))

        self.menuBar.addmenuitem('Residue', 'command', 'Tyrosine',
                                         label='Tyrosine [Alt-Y]',
                                         command = lambda s=self: s.cmd.do(
            "_ editor.attach_amino_acid('pk1','tyr')"))

        self.menuBar.addmenuitem('Residue', 'command', 'Valine',
                                         label='Valine [Alt-V]',
                                         command = lambda s=self: s.cmd.do(
            "_ editor.attach_amino_acid('pk1','val')"))


        self.menuBar.addmenuitem('Residue', 'separator', '')
        self.menuBar.addmenuitem('Residue', 'command', 'Helix',
                                         label='Helix',
                                         command = lambda s=self: s.cmd.do("_ set secondary_structure,1"))

        self.menuBar.addmenuitem('Residue', 'command', 'Antiparallel Beta Sheet',
                                         label='Antiparallel Beta Sheet',
                                         command = lambda s=self: s.cmd.do("_ set secondary_structure,2"))

        self.menuBar.addmenuitem('Residue', 'command', 'Parallel Beta Sheet',
                                         label='Parallel Beta Sheet',
                                         command = lambda s=self: s.cmd.do("_ set secondary_structure,3"))

        self.menuBar.addmenuitem('Build', 'separator', '')

        self.menuBar.addcascademenu('Build', 'Sculpting', 'Sculpting',
                                             label='Sculpting',tearoff=TRUE)

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
                                         command = lambda s=self: s.cmd.do("_ sculpt_activate all"))

        self.menuBar.addmenuitem('Sculpting', 'command', 'Deactivate',
                                         label='Deactivate',
                                         command = lambda s=self: s.cmd.do("_ sculpt_deactivate all"))

        self.menuBar.addmenuitem('Sculpting', 'command', 'Clear Memory',
                                         label='Clear Memory',
                                         command = lambda s=self: s.cmd.do("_ sculpt_purge"))

        self.menuBar.addmenuitem('Sculpting', 'separator', '')

        self.menuBar.addmenuitem('Sculpting', 'command', '1 Cycle/Update',
                                         label='1 Cycle per Update',
                                         command = lambda s=self: s.cmd.do("_ set sculpting_cycles=1"))

        self.menuBar.addmenuitem('Sculpting', 'command', '3 Cycles/Update',
                                         label='3 Cycles per Update',
                                         command = lambda s=self: s.cmd.do("_ set sculpting_cycles=3"))

        self.menuBar.addmenuitem('Sculpting', 'command', '10 Cycles/Update',
                                         label='10 Cycles per Update',
                                         command = lambda s=self: s.cmd.do("_ set sculpting_cycles=10"))

        self.menuBar.addmenuitem('Sculpting', 'command', '33 Cycles/Update',
                                         label='33 Cycles per Update',
                                         command = lambda s=self: s.cmd.do("_ set sculpting_cycles=33"))

        self.menuBar.addmenuitem('Sculpting', 'command', '100 Cycles/Update',
                                         label='100 Cycles per Update',
                                         command = lambda s=self: s.cmd.do("_ set sculpting_cycles=100"))

        self.menuBar.addmenuitem('Sculpting', 'command', '333 Cycles/Update',
                                         label='333 Cycles per Update',
                                         command = lambda s=self: s.cmd.do("_ set sculpting_cycles=333"))

        self.menuBar.addmenuitem('Sculpting', 'command', '1000 Cycles/Update',
                                         label='1000 Cycles per Update',
                                         command = lambda s=self: s.cmd.do("_ set sculpting_cycles=1000"))

        self.menuBar.addmenuitem('Sculpting', 'separator', '')

#define cSculptBond  0x01
#define cSculptAngl  0x02
#define cSculptPyra  0x04
#define cSculptPlan  0x08
#define cSculptLine  0x10
#define cSculptVDW   0x20
#define cSculptVDW14 0x40
#define cSculptTors  0x80

        self.menuBar.addmenuitem('Sculpting', 'command', 'Bonds Only',
                                         label='Bonds Only',
                                         command = lambda s=self: s.cmd.do("_ set sculpt_field_mask=%d"%(
            0x01)))

        self.menuBar.addmenuitem('Sculpting', 'command', 'Bonds & Angles Only',
                                         label='Bonds & Angles Only',
                                         command = lambda s=self: s.cmd.do("_ set sculpt_field_mask=%d"%(
            0x01+0x02)))
        
        self.menuBar.addmenuitem('Sculpting', 'command', 'Local Geometry Only',
                                         label='Local Geometry Only',
                                         command = lambda s=self: s.cmd.do("_ set sculpt_field_mask=%d"%(
            0x01+0x02+0x04+0x08+0x10)))

        self.menuBar.addmenuitem('Sculpting', 'command', 'All Except VDW',
                                         label='All Except VDW',
                                         command = lambda s=self: s.cmd.do("_ set sculpt_field_mask=%d"%(
            0x01+0x02+0x04+0x08+0x10+0x80)))

        self.menuBar.addmenuitem('Sculpting', 'command', 'All Except 1-4 VDW & Torsions',
                                         label='All Except 1-4 VDW & Torsions',
                                         command = lambda s=self: s.cmd.do("_ set sculpt_field_mask=%d"%(
            0x01+0x02+0x04+0x08+0x10+0x20)))

        self.menuBar.addmenuitem('Sculpting', 'command', 'All Terms',
                                         label='All Terms',
                                         command = lambda s=self: s.cmd.do("_ set sculpt_field_mask=%d"%(
            0xFF)))


        self.menuBar.addmenuitem('Build', 'separator', '')
        
        self.menuBar.addmenuitem('Build', 'command', 'Cycle Bond Valence',
                                         label='Cycle Bond Valence [Ctrl-W]',
                                         command = lambda s=self: s.cmd.do("_ cycle_valence"))

        self.menuBar.addmenuitem('Build', 'command', 'Fill Hydrogens',
                                         label='Fill Hydrogens on (pk1) [Ctrl-R]',
                                         command = lambda s=self: s.cmd.do("_ h_fill"))

        self.menuBar.addmenuitem('Build', 'command', 'Invert',
                                         label='Invert (pk2)-(pk1)-(pk3) [Ctrl-E]',
                                         command = lambda s=self: s.cmd.do("_ invert"))

        self.menuBar.addmenuitem('Build', 'command', 'Form Bond',
                                         label='Create Bond (pk1)-(pk2) [Ctrl-T]',
                                         command = lambda s=self: s.cmd.do("_ bond"))


        self.menuBar.addmenuitem('Build', 'separator', '')

        
        self.menuBar.addmenuitem('Build', 'command', 'Remove (pk1)',
                                         label='Remove (pk1) [Ctrl-D]',
                                         command = lambda s=self: s.cmd.do("_ remove pk1"))

        self.menuBar.addmenuitem('Build', 'separator', '')
        
        self.menuBar.addmenuitem('Build', 'command', 'Make Positive',
                                         label='Make (pk1) Positive [Ctrl-K]',
                                         command = lambda s=self: s.cmd.do("_ alter pk1,formal_charge=1.0"))

        self.menuBar.addmenuitem('Build', 'command', 'Make Negative',
                                         label='Make (pk1) Negative [Ctrl-J]',
                                         command = lambda s=self: s.cmd.do("_ alter pk1,formal_charge=-1.0"))

        self.menuBar.addmenuitem('Build', 'command', 'Make Neutral',
                                         label='Make (pk1) Neutral [Ctrl-U]',
                                         command = lambda s=self: s.cmd.do("_ alter pk1,formal_charge=-0.0"))


        self.menuBar.addmenu('Movie', 'Movie Control',tearoff=TRUE)

        
        self.menuBar.addcascademenu('Movie', 'Speed', 'Playback Speed',
                                             label=self.pad+'Speed')

        self.menuBar.addmenuitem('Speed', 'command', 'Maximum',
                                         label=self.pad+'Maximum',
                                         command = lambda s=self: s.cmd.set("movie_delay","0",log=1))

        self.menuBar.addmenuitem('Speed', 'command', '30 FPS',
                                         label=self.pad+'30 FPS',
                                         command = lambda s=self: s.cmd.set("movie_delay","33",log=1))

        self.menuBar.addmenuitem('Speed', 'command', '15 FPS',
                                         label=self.pad+'15 FPS',
                                         command = lambda s=self: s.cmd.set("movie_delay","66",log=1))

        self.menuBar.addmenuitem('Speed', 'command', '5 FPS',
                                         label=self.pad+'5 FPS',
                                         command = lambda s=self: s.cmd.set("movie_delay","200",log=1))

        self.menuBar.addmenuitem('Speed', 'command', '1 FPS',
                                         label=self.pad+'1 FPS',
                                         command = lambda s=self: s.cmd.set("movie_delay","1000",log=1))

        self.menuBar.addmenuitem('Speed', 'command', '0.3 FPS',
                                         label=self.pad+'0.3 FPS',
                                         command = lambda s=self: s.cmd.set("movie_delay","3000",log=1))

        self.menuBar.addmenuitem('Movie', 'command', 'Reset Meter',
                                         label=self.pad+'Reset Meter',
                                         command = lambda s=self: s.cmd.do("_ meter_reset"))

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
                                         command = lambda s=self: s.cmd.mclear())

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

        self.menuBar.addmenu('Display', 'Display Control',tearoff=TRUE)

        self.menuBar.addmenuitem('Display', 'checkbutton',
                                 'Sequence',
                                 label='Sequence',
                                variable = self.setting.seq_view,
                                command = lambda s=self: s.setting.update('seq_view'))      

        self.menuBar.addcascademenu('Display', 'Sequence', 'Sequence Mode',
                                         label='Sequence Mode')

        self.menuBar.addmenuitem('Sequence', 'command', 'All Residue Numbers',
                                         label='All Residue Numbers',
                                         command = lambda s=self: s.cmd.do("_ set seq_view_label_mode,2"))

        self.menuBar.addmenuitem('Sequence', 'command', 'Top Sequence Only',
                                 label='Top Sequence Only',
                                 command = lambda s=self: s.cmd.do("_ set seq_view_label_mode,1"))

        self.menuBar.addmenuitem('Sequence', 'command', 'Object Names Only',
                                 label='Object Names Only',
                                 command = lambda s=self: s.cmd.do("_ set seq_view_label_mode,0"))

        self.menuBar.addmenuitem('Sequence', 'command', 'No Labels',
                                 label='No Labels',
                                 command = lambda s=self: s.cmd.do("_ set seq_view_label_mode,3"))

        self.menuBar.addmenuitem('Sequence', 'separator', '')

        self.menuBar.addmenuitem('Sequence', 'command', 'Residue Codes',
                                         label='Residue Codes',
                                         command = lambda s=self: s.cmd.do("_ set seq_view_format,0"))

        self.menuBar.addmenuitem('Sequence', 'command', 'Residue Names',
                                         label='Residues',
                                         command = lambda s=self: s.cmd.do("_ set seq_view_format,1"))

        self.menuBar.addmenuitem('Sequence', 'command', 'Chain Identifiers',
                                         label='Chains',
                                         command = lambda s=self: s.cmd.do("_ set seq_view_format,3"))

        self.menuBar.addmenuitem('Sequence', 'command', 'Atom Names',
                                         label='Atoms',
                                         command = lambda s=self: s.cmd.do("_ set seq_view_format,2"))

        self.menuBar.addmenuitem('Sequence', 'command', 'States',
                                         label='States',
                                         command = lambda s=self: s.cmd.do("_ set seq_view_format,4"))


        self.menuBar.addmenuitem('Display', 'checkbutton',
                                 'Stereo',
                                 label='Stereo',
                                variable = self.setting.stereo,
                                command = lambda s=self: s.cmd.do("_ stereo "+
                                                                          ('off','on')[s.setting.stereo.get()]))
        
#      self.menuBar.addmenuitem('Display', 'command', 'Stereo On',
#                               label='Stereo On',
#                              command = lambda s=self: s.cmd.do("_ stereo on"))

#      self.menuBar.addmenuitem('Display', 'command', 'Stereo Off',
#                               label='Stereo Off',
#                               command = lambda s=self: s.cmd.do("_ stereo off"))

        self.menuBar.addcascademenu('Display', 'Stereo', 'Stereo Mode',
                                         label='Stereo Mode')

        self.menuBar.addmenuitem('Stereo', 'command', 'Cross-Eye Stereo',
                                         label='Cross-Eye Stereo',
                                         command = lambda s=self: s.cmd.do("_ stereo crosseye"))

        self.menuBar.addmenuitem('Stereo', 'command', 'Wall-Eye Stereo',
                                         label='Wall-Eye Stereo',
                                         command = lambda s=self: s.cmd.do("_ stereo walleye"))

        self.menuBar.addmenuitem('Stereo', 'command', 'Quad-Buffered Stereo',
                                         label='Quad-Buffered Stereo',
                                         command = lambda s=self: s.cmd.do("_ stereo quadbuffer"))

        self.menuBar.addmenuitem('Stereo', 'separator', '')


        self.menuBar.addmenuitem('Stereo', 'command', 'Swap Sides',
                                         label='Swap Sides',
                                         command = lambda s=self: s.cmd.do("_ stereo swap"))

        self.menuBar.addmenuitem('Display', 'separator', '')
        
        self.menuBar.addcascademenu('Display', 'Zoom', 'Zoom',
                                             label='Zoom')

        self.menuBar.addmenuitem('Zoom', 'command', '4 Angstrom Sphere',
                                         label='4 Angstrom Sphere',
                                         command = lambda s=self: s.cmd.do("_ zoom center,4"))

        self.menuBar.addmenuitem('Zoom', 'command', '6 Angstrom Sphere',
                                         label='6 Angstrom Sphere',
                                         command = lambda s=self: s.cmd.do("_ zoom center,6"))

        self.menuBar.addmenuitem('Zoom', 'command', '8 Angstrom Sphere',
                                         label='8 Angstrom Sphere',
                                         command = lambda s=self: s.cmd.do("_ zoom center,8"))

        self.menuBar.addmenuitem('Zoom', 'command', '12 Angstrom Sphere',
                                         label='12 Angstrom Sphere',
                                         command = lambda s=self: s.cmd.do("_ zoom center,12"))

        self.menuBar.addmenuitem('Zoom', 'command', '20 Angstrom Sphere',
                                         label='20 Angstrom Sphere',
                                         command = lambda s=self: s.cmd.do("_ zoom center,20"))

        self.menuBar.addmenuitem('Zoom', 'command', 'All',
                                         label='All',
                                         command = lambda s=self: s.cmd.do("_ zoom all"))

        self.menuBar.addmenuitem('Zoom', 'command', 'Complete',
                                         label='Complete',
                                         command = lambda s=self: s.cmd.do("_ zoom all,complete=1"))

        self.menuBar.addcascademenu('Display', 'Clip', 'Clip',
                                             label='Clip')

        self.menuBar.addmenuitem('Clip', 'command', 'Nothing',
                                         label='Nothing',
                                         command = lambda s=self: s.cmd.do("_ clip atoms,5,all"))

        self.menuBar.addmenuitem('Clip', 'command', '8 Angstrom Slab',
                                         label='8 Angstrom Slab',
                                         command = lambda s=self: s.cmd.do("_ clip slab,8"))

        self.menuBar.addmenuitem('Clip', 'command', '12 Angstrom Slab',
                                         label='12 Angstrom Slab',
                                         command = lambda s=self: s.cmd.do("_ clip slab,10"))

        self.menuBar.addmenuitem('Clip', 'command', '16 Angstrom Slab',
                                         label='16 Angstrom Slab',
                                         command = lambda s=self: s.cmd.do("_ clip slab,15"))

        self.menuBar.addmenuitem('Clip', 'command', '20 Angstrom Slab',
                                         label='20 Angstrom Slab',
                                         command = lambda s=self: s.cmd.do("_ clip slab,20"))

        self.menuBar.addmenuitem('Clip', 'command', '30 Angstrom Slab',
                                         label='30 Angstrom Slab',
                                         command = lambda s=self: s.cmd.do("_ clip slab,30"))


        self.menuBar.addmenuitem('Display', 'separator', '')

        self.menuBar.addcascademenu('Display', 'Background', 'Background',
                                             label='Background')

        self.menuBar.addmenuitem('Background', 'command', 'White Background',
                                         label='White',
                                         command = lambda s=self: s.cmd.do("_ cmd.bg_color('white')"))

        self.menuBar.addmenuitem('Background', 'command', 'Light Grey',
                                         label='Light Grey',
                                         command = lambda s=self: s.cmd.do("_ cmd.bg_color('grey80')"))

        self.menuBar.addmenuitem('Background', 'command', 'Grey',
                                         label='Grey',
                                         command = lambda s=self: s.cmd.do("_ cmd.bg_color('grey50')"))


        self.menuBar.addmenuitem('Background', 'command', 'Black Background',
                                         label='Black',
                                         command = lambda s=self: s.cmd.do("_ cmd.bg_color('black')"))

        self.menuBar.addmenuitem('Background', 'separator', '')
        
        self.menuBar.addmenuitem('Background', 'checkbutton',
                                 'Opaque Background Color',
                                 label=self.pad+'Opaque',
                                variable = self.setting.opaque_background,
                                command = lambda s=self: s.setting.update('opaque_background'))

        self.menuBar.addmenuitem('Background', 'checkbutton',
                                 'Show Alpha Checker',
                                 label=self.pad+'Show Alpha Checker',
                                variable = self.setting.show_alpha_checker,
                                command = lambda s=self: s.setting.update('show_alpha_checker'))

        self.menuBar.addcascademenu('Display', 'Color Space', 'Color Space',
                                             label='Color Space')

        self.menuBar.addmenuitem('Color Space', 'command', 'CMYK (for publications)',
                                         label='CMYK (for publications)',
                                         command = lambda s=self: s.cmd.do("_ cmd.space('cmyk')"))

        self.menuBar.addmenuitem('Color Space', 'command', 'PyMOL (for video & web)',
                                         label='PyMOL (for video & web)',
                                         command = lambda s=self: s.cmd.do("_ cmd.space('pymol')"))

        self.menuBar.addmenuitem('Color Space', 'command', 'RGB (default)',
                                         label='RGB (default)',
                                         command = lambda s=self: s.cmd.do("_ cmd.space('rgb')"))

        self.menuBar.addcascademenu('Display', 'Performance', 'Quality',
                                             label='Quality')

        self.menuBar.addmenuitem('Performance', 'command', 'Maximum Performance',
                                         label='Maximum Performance',
                                         command = lambda s=self: s.cmd.do("_ util.performance(100)"))

        self.menuBar.addmenuitem('Performance', 'command', 'Reasonable Performance',
                                         label='Reasonable Performance',
                                         command = lambda s=self: s.cmd.do("_ util.performance(66)"))
        
        self.menuBar.addmenuitem('Performance', 'command', 'Reasonable Quality',
                                         label='Reasonable Quality',
                                         command = lambda s=self: s.cmd.do("_ util.performance(33)"))

        self.menuBar.addmenuitem('Performance', 'command', 'Maximum Quality',
                                         label='Maximum Quality',
                                         command = lambda s=self: s.cmd.do("_ util.performance(0)"))

        
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
                                 'Depth Cue (Fogging).',
                                 label=self.pad+'Depth Cue',
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

        self.menuBar.addmenuitem('Display', 'checkbutton',
                                 'Use Display Lists.',
                                 label=self.pad+'Use Display Lists',
                                variable = self.setting.use_display_lists,
                                command = lambda s=self: s.setting.update('use_display_lists'))

        self.menuBar.addmenuitem('Display', 'checkbutton',
                                 'Texture Fonts',
                                 label=self.pad+'Texture Fonts',
                                variable = self.setting.texture_fonts,
                                command = lambda s=self: s.setting.update('texture_fonts'))

        self.menuBar.addmenuitem('Display', 'checkbutton',
                                 'Animation',
                                 label=self.pad+'Animation',
                                variable = self.setting.animation,
                                command = lambda s=self: s.setting.update('animation'))

        self.menuBar.addmenu('Setting', 'Configuration Control',tearoff=TRUE)

        self.menuBar.addmenuitem('Setting', 'command',
                                 'Edit PyMOL Settings',
                                 label=self.pad+'Edit All...',
                                         command = lambda s=self: SetEditor(s))

        self.menuBar.addmenuitem('Setting', 'command',
                                 'Edit PyMOL Colors',
                                 label=self.pad+'Colors...',
                                         command = lambda s=self: ColorEditor(s))

        self.menuBar.addcascademenu('Setting', 'Cartoon', 'Cartoon',
                                             label=self.pad+'Cartoon')

        self.menuBar.addmenuitem('Cartoon', 'checkbutton',
                                 'Side Chain Helper',
                                 label=self.pad+'Side Chain Helper',
                                variable = self.setting.cartoon_side_chain_helper,
                                command = lambda s=self: s.setting.update('cartoon_side_chain_helper'))

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

        self.menuBar.addmenuitem('Cartoon', 'checkbutton',
                                 'Highlight Color',
                                 label=self.pad+'Highlight Color',
                                variable = self.setting.cartoon_highlight_color,
                                command = lambda s=self: s.setting.update('cartoon_highlight_color'))


        self.menuBar.addcascademenu('Setting', 'Transparency', 'Transparency',
                                             label=self.pad+'Transparency')

        self.transparency_menu('CartoonTransparency','Cartoon','cartoon_transparency')
        self.transparency_menu('SurfaceTransparency','Surface','transparency')
        self.transparency_menu('StickTransparency','Stick','stick_transparency')
        self.transparency_menu('SphereTransparency','Sphere','sphere_transparency')      

        self.menuBar.addmenuitem('Transparency', 'separator', '')
        

        self.menuBar.addmenuitem('Transparency', 'command', 'Uni-Layer',
                                         label='Uni-Layer',
                                         command = lambda s=self: s.cmd.do(
            "_ cmd.set('transparency_mode',2);cmd.set('backface_cull',1);cmd.set('two_sided_lighting',0)"))

        self.menuBar.addmenuitem('Transparency', 'command', 'Multi-Layer',
                                         label='Multi-Layer',
                                         command = lambda s=self: s.cmd.do(
            "_ cmd.set('transparency_mode',1);cmd.set('backface_cull',0);cmd.set('two_sided_lighting',1)"))

        self.menuBar.addmenuitem('Transparency', 'command', 'Fast and Ugly',
                                         label='Fast and Ugly',
                                         command = lambda s=self: s.cmd.do(
            "_ cmd.set('transparency_mode',0);cmd.set('backface_cull',1);cmd.set('two_sided_lighting',0)"))

                
        self.menuBar.addcascademenu('Setting', 'Rendering', 'Rendering',
                                             label=self.pad+'Rendering')

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
                                         command = lambda s=self: s.cmd.do("_ util.ray_shadows('none')"))

        self.menuBar.addmenuitem('Shadows', 'command', 'Light',
                                         label='Light',
                                         command = lambda s=self: s.cmd.do("_ util.ray_shadows('light')"))

        self.menuBar.addmenuitem('Shadows', 'command', 'Medium',
                                         label='Medium',
                                         command = lambda s=self: s.cmd.do("_ util.ray_shadows('medium')"))

        self.menuBar.addmenuitem('Shadows', 'command', 'Heavy',
                                         label='Heavy',
                                         command = lambda s=self: s.cmd.do("_ util.ray_shadows('heavy')"))

        self.menuBar.addmenuitem('Shadows', 'command', 'Black',
                                         label='Black',
                                         command = lambda s=self: s.cmd.do("_ util.ray_shadows('black')"))

        self.menuBar.addmenuitem('Shadows', 'separator', '')
        
        self.menuBar.addmenuitem('Shadows', 'command', 'Matte',
                                         label='Matte',
                                         command = lambda s=self: s.cmd.do("_ util.ray_shadows('matte')"))
        self.menuBar.addmenuitem('Shadows', 'command', 'Soft',
                                         label='Soft',
                                         command = lambda s=self: s.cmd.do("_ util.ray_shadows('soft')"))

        self.menuBar.addmenuitem('Shadows', 'command', 'Occlusion',
                                        label='Occlusion',
                                         command = lambda s=self: s.cmd.do("_ util.ray_shadows('occlusion')"))



        self.menuBar.addcascademenu('Rendering', 'Texture', 'Texture',
                                         label=self.pad+'Texture')

        self.menuBar.addmenuitem('Texture', 'command', 'None',
                                         label='None',
                                         command = lambda s=self: s.cmd.do("_ cmd.set('ray_texture',0)"))

        self.menuBar.addmenuitem('Texture', 'command', 'Matte 1',
                                         label='Matte 1',
                                         command = lambda s=self: s.cmd.do("_ cmd.set('ray_texture',1)"))

        self.menuBar.addmenuitem('Texture', 'command', 'Matte 2',
                                         label='Matte 2',
                                         command = lambda s=self: s.cmd.do("_ cmd.set('ray_texture',4)"))

        self.menuBar.addmenuitem('Texture', 'command', 'Swirl 1',
                                         label='Swirl 1',
                                         command = lambda s=self: s.cmd.do("_ cmd.set('ray_texture',2)"))

        self.menuBar.addmenuitem('Texture', 'command', 'Swirl 2',
                                         label='Swirl 2',
                                         command = lambda s=self: s.cmd.do("_ cmd.set('ray_texture',3)"))

        self.menuBar.addmenuitem('Texture', 'command', 'Fiber',
                                         label='Fiber',
                                         command = lambda s=self: s.cmd.do("_ cmd.set('ray_texture',5)"))

        self.menuBar.addcascademenu('Rendering', 'Interior Texture', 'Interior Texture',
                                         label=self.pad+'Interior Texture')

        self.menuBar.addmenuitem('Interior Texture', 'command', 'None',
                                         label='None',
                                         command = lambda s=self: s.cmd.do("_ cmd.set('ray_interior_texture',0)"))

        self.menuBar.addmenuitem('Interior Texture', 'command', 'Matte 1',
                                         label='Matte 1',
                                         command = lambda s=self: s.cmd.do("_ cmd.set('ray_interior_texture',1)"))

        self.menuBar.addmenuitem('Interior Texture', 'command', 'Matte 2',
                                         label='Matte 2',
                                         command = lambda s=self: s.cmd.do("_ cmd.set('ray_interior_texture',4)"))

        self.menuBar.addmenuitem('Interior Texture', 'command', 'Swirl 1',
                                         label='Swirl 1',
                                         command = lambda s=self: s.cmd.do("_ cmd.set('ray_interior_texture',2)"))

        self.menuBar.addmenuitem('Interior Texture', 'command', 'Swirl 2',
                                         label='Swirl 2',
                                         command = lambda s=self: s.cmd.do("_ cmd.set('ray_interior_texture',3)"))

        self.menuBar.addmenuitem('Interior Texture', 'command', 'Fiber',
                                         label='Fiber',
                                         command = lambda s=self: s.cmd.do("_ cmd.set('ray_interior_texture',5)"))


        self.menuBar.addcascademenu('Rendering', 'Memory', 'Memory',
                                         label=self.pad+'Memory')

        self.menuBar.addmenuitem('Memory', 'command', 'Use Less (slower)',
                                         label='Use Less (slower)',
                                         command = lambda s=self: s.cmd.do("_ cmd.set('hash_max',70,quiet=0)"))

        self.menuBar.addmenuitem('Memory', 'command', 'Use Standard Amount',
                                         label='Use Standard Amount',
                                         command = lambda s=self: s.cmd.do("_ cmd.set('hash_max',100,quiet=0)"))

        self.menuBar.addmenuitem('Memory', 'command', 'Use More (faster)',
                                         label='Use More (faster)',
                                         command = lambda s=self: s.cmd.do("_ cmd.set('hash_max',170,quiet=0)"))

        self.menuBar.addmenuitem('Memory', 'command', 'Use Even More',
                                         label='Use Even More',
                                         command = lambda s=self: s.cmd.do("_ cmd.set('hash_max',230,quiet=0)"))

        self.menuBar.addmenuitem('Memory', 'command', 'Use Most',
                                         label='Use Most',
                                         command = lambda s=self: s.cmd.do("_ cmd.set('hash_max',300,quiet=0)"))

        self.menuBar.addmenuitem('Rendering', 'separator', '')

        self.menuBar.addmenuitem('Rendering', 'checkbutton',
                                 'Cull Backfaces when Rendering',
                                 label=self.pad+'Cull Backfaces',
                                variable = self.setting.backface_cull,
                                command = lambda s=self: s.setting.update('backface_cull'))


        self.menuBar.addmenuitem('Rendering', 'checkbutton',
                                 'Opaque Interior Colors',
                                 label=self.pad+'Opaque Interiors',
                                variable = self.setting.ray_interior_color,
                                command = lambda s=self: s.setting.update('ray_interior_color'))

        self.menuBar.addmenuitem('Setting', 'separator', '')

        self.menuBar.addcascademenu('Setting', 'Output', 'Output Size',
                                             label=self.pad+'Output Size')

        self.menuBar.addmenuitem('Output', 'command', '8',
                                         label='8 Point',
                                         command = lambda s=self:
                                         s.text.configure(font=(s.font,8)))

        self.menuBar.addmenuitem('Output', 'command', '9',
                                         label='9 Point',
                                         command = lambda s=self:
                                         s.text.configure(font=(s.font,9)))

        self.menuBar.addmenuitem('Output', 'command', '10',
                                         label='10 Point',
                                         command = lambda s=self:
                                         s.text.configure(font=(s.font,10)))

        self.menuBar.addmenuitem('Output', 'command', '11',
                                         label='11 Point',
                                         command = lambda s=self:
                                         s.text.configure(font=(s.font,11)))

        self.menuBar.addmenuitem('Output', 'command', '12',
                                         label='12 Point',
                                         command = lambda s=self:
                                         s.text.configure(font=(s.font,12)))

        self.menuBar.addcascademenu('Setting', 'Control', 'Control Size',
                                             label=self.pad+'Control Size')

        self.menuBar.addmenuitem('Control', 'command', '12',
                                         label='12',
                                         command = lambda s=self:
                                         s.cmd.do("_ set internal_gui_control_size,12,quiet=1"))

        self.menuBar.addmenuitem('Control', 'command', '14',
                                         label='14',
                                         command = lambda s=self:
                                         s.cmd.do("_ set internal_gui_control_size,14,quiet=1"))

        self.menuBar.addmenuitem('Control', 'command', '16',
                                         label='16',
                                         command = lambda s=self:
                                         s.cmd.do("_ set internal_gui_control_size,16,quiet=1"))

        self.menuBar.addmenuitem('Control', 'command', '18',
                                         label='18 (default)',
                                         command = lambda s=self:
                                         s.cmd.do("_ set internal_gui_control_size,18,quiet=1"))

        self.menuBar.addmenuitem('Control', 'command', '20',
                                         label='20',
                                         command = lambda s=self:
                                         s.cmd.do("_ set internal_gui_control_size,20,quiet=1"))

        self.menuBar.addmenuitem('Control', 'command', '24',
                                         label='24',
                                         command = lambda s=self:
                                         s.cmd.do("_ set internal_gui_control_size,24,quiet=1"))

        self.menuBar.addmenuitem('Control', 'command', '30',
                                         label='30',
                                         command = lambda s=self:
                                         s.cmd.do("_ set internal_gui_control_size,30,quiet=1"))

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
                                         label=self.pad+'Show Text',
                                         command = lambda s=self: s.cmd.set("text","1",log=1))

        self.menuBar.addmenuitem('Setting', 'command', 'Hide Text Output',
                                         label=self.pad+'Hide Text',
                                         command = lambda s=self: s.cmd.set("text","0",log=1))

        self.menuBar.addmenuitem('Setting', 'checkbutton',
                                 'Overlay Text Output on Graphics',
                                 label=self.pad+'Overlay Text',
                                variable = self.setting.overlay,
                                command = lambda s=self: s.setting.update('overlay'))

        self.menuBar.addmenu('Scene', 'Scene Storage',tearoff=TRUE)

        self.menuBar.addmenuitem('Scene', 'command', 'Next',
                                         label=self.pad+'Next [PgDn]',
                                         command = lambda s=self: s.cmd.scene('auto','next'))

        self.menuBar.addmenuitem('Scene', 'command', 'Previous',
                                         label=self.pad+'Previous [PgUp]',
                                         command = lambda s=self: s.cmd.scene('auto','previous'))

        self.menuBar.addmenuitem('Scene', 'separator', '')
        
        self.menuBar.addmenuitem('Scene', 'command', 'Append',
                                         label=self.pad+'Append',
                                         command = lambda s=self: s.cmd.scene('new','store'))

        self.menuBar.addmenuitem('Scene', 'command', 'Insert Before',
                                         label=self.pad+'Insert (before)',
                                         command = lambda s=self: s.cmd.scene('','insert_before'))

        self.menuBar.addmenuitem('Scene', 'command', 'Insert After',
                                         label=self.pad+'Insert (after)',
                                         command = lambda s=self: s.cmd.scene('','insert_after'))

        self.menuBar.addmenuitem('Scene', 'command', 'Update',
                                         label=self.pad+'Update',
                                         command = lambda s=self: s.cmd.scene('auto','update'))

#      self.menuBar.addmenuitem('Scene', 'command', 'Annotate',
#                               label=self.pad+'Append',
#                               command = lambda s=self: s.cmd.scene('new','store'))

        self.menuBar.addmenuitem('Scene', 'separator', '')

        self.menuBar.addmenuitem('Scene', 'command', 'Delete',
                                         label=self.pad+'Delete',
                                         command = lambda s=self: s.cmd.scene('auto','clear'))

        self.menuBar.addmenuitem('Scene', 'separator', '')

        self.menuBar.addcascademenu('Scene', 'Recall', 'Recall',
                                             label=self.pad+'Recall')

        self.menuBar.addcascademenu('Scene', 'Store', 'Store',
                                             label=self.pad+'Store')

#      self.menuBar.addcascademenu('Store', 'StoreSHFT', 'StoreSHFT',
#                                  label=self.pad+'Shift')

        self.menuBar.addcascademenu('Scene', 'Clear', 'Clear',
                                             label=self.pad+'Clear')

#      self.menuBar.addcascademenu('Scene', 'SceneSHFT', 'SceneSHFT',
#                                  label=self.pad+'Shift')

        for x in range(1,13):
            self.menuBar.addmenuitem('Store', 'checkbutton', 'F%d'%x,
                                             label=self.pad+'F%d'%x,
                                             variable = self.setting.F[x],
                                             command = lambda x=x,s=self: s.cmd.do("scene F%d,store"%x))

            self.menuBar.addmenuitem('Recall', 'checkbutton', 'Recall F%d'%x,
                                             label=self.pad+'F%d'%x,
                                             variable = self.setting.F[x],
                                             command = lambda x=x,s=self: s.cmd.do("scene F%d"%x))

            self.menuBar.addmenuitem('Clear', 'checkbutton', 'F%d'%x,
                                         label=self.pad+'F%d'%x,
                                         variable = self.setting.F[x],
                                         command = lambda x=x,s=self: s.cmd.do("scene F%d,clear"%x))

#         self.menuBar.addmenuitem('ClearSHFT', 'checkbutton', 'SHFT-F%d'%x,
#                                  label=self.pad+'SHFT-F%d'%x,
#                                  variable = self.setting.SHFTF[x],
#                                  command = lambda x=x,s=self: s.cmd.do("scene SHFT-F%d,clear"%x))
            
        self.menuBar.addmenu('Mouse', 'Mouse Configuration',tearoff=TRUE)

        self.menuBar.addcascademenu('Mouse', 'SelectionMode', 'Selection Mode',
                                             label='Selection Mode')

        self.menuBar.addmenuitem('SelectionMode', 'command', 'Atoms',
                                         label='Atoms',
                                         command = lambda s=self:
                                         s.cmd.do("_ set mouse_selection_mode,0,quiet=1"))

        self.menuBar.addmenuitem('SelectionMode', 'command', 'Residues',
                                         label='Residues',
                                         command = lambda s=self:
                                         s.cmd.do("_ set mouse_selection_mode,1,quiet=1"))

        self.menuBar.addmenuitem('SelectionMode', 'command', 'Chains',
                                         label='Chains',
                                         command = lambda s=self:
                                         s.cmd.do("_ set mouse_selection_mode,2,quiet=1"))

        self.menuBar.addmenuitem('SelectionMode', 'command', 'Segments',
                                         label='Segments',
                                         command = lambda s=self:
                                         s.cmd.do("_ set mouse_selection_mode,3,quiet=1"))

        self.menuBar.addmenuitem('SelectionMode', 'command', 'Objects',
                                         label='Objects',
                                         command = lambda s=self:
                                         s.cmd.do("_ set mouse_selection_mode,4,quiet=1"))

        self.menuBar.addmenuitem('SelectionMode', 'separator', '')

        self.menuBar.addmenuitem('SelectionMode', 'command', 'Molecules',
                                         label='Molecules',
                                         command = lambda s=self:
                                         s.cmd.do("_ set mouse_selection_mode,5,quiet=1"))

        self.menuBar.addmenuitem('SelectionMode', 'separator', '')

        self.menuBar.addmenuitem('SelectionMode', 'command', 'C-alphas',
                                         label='C-alphas',
                                         command = lambda s=self:
                                         s.cmd.do("_ set mouse_selection_mode,6,quiet=1"))

        self.menuBar.addmenuitem('Mouse', 'separator', '')
        self.menuBar.addmenuitem('Mouse', 'command', '3 Button Viewing Mode',
                                         label='3 Button Viewing Mode',
                                         command = lambda s=self: s.cmd.mouse('three_button_viewing'))

        self.menuBar.addmenuitem('Mouse', 'command', '3 Button Editing Mode',
                                         label='3 Button Editing Mode',
                                         command = lambda s=self: s.cmd.mouse('three_button_editing'))

        self.menuBar.addmenuitem('Mouse', 'command', '2 Button Viewing Mode',
                                         label='2 Button Viewing Mode',
                                         command = lambda s=self: s.cmd.mouse('two_button_viewing'))

        self.menuBar.addmenuitem('Mouse', 'command', '2 Button Selecting Mode',
                                         label='2 Button Selecting Mode',
                                         command = lambda s=self: s.cmd.mouse('two_button_selecting'))

        self.menuBar.addmenuitem('Mouse', 'command', '2 Button Editing Mode',
                                         label='2 Button Editing Mode',
                                         command = lambda s=self: s.cmd.mouse('two_button_editing'))

        self.menuBar.addmenuitem('Mouse', 'command', '1 Button Viewing Mode',
                                         label='1 Button Viewing Mode',
                                         command = lambda s=self: s.cmd.mouse('one_button_viewing'))

        self.menuBar.addmenuitem('Mouse', 'separator', '')

        self.menuBar.addmenuitem('Mouse', 'checkbutton',
                                 'Virtual Trackball.',
                                 label=self.pad+'Virtual Trackball',
                                variable = self.setting.virtual_trackball,
                                command = lambda s=self: s.setting.update('virtual_trackball'))

        self.menuBar.addmenuitem('Mouse', 'checkbutton',
                                 'Roving Origin.',
                                 label=self.pad+'Roving Origin',
                                variable = self.setting.roving_origin,
                                command = lambda s=self: s.setting.update('roving_origin'))

        self.menuBar.addmenuitem('Mouse', 'checkbutton',
                                 'Roving Detail.',
                                 label=self.pad+'Roving Detail',
                                variable = self.setting.roving_detail,
                                command = lambda s=self: s.setting.update('roving_detail'))

        self.menuBar.addmenuitem('Mouse', 'separator', '')
        

        self.menuBar.addmenuitem('Mouse', 'command', '3 Button Editing Cycle',
                                         label='3 Button Editing Cycle',
                                         command = lambda s=self: s.cmd.config_mouse('three_button'))

#        self.menuBar.addmenuitem('Mouse', 'command', '3 Button Motions Cycle',
#                                         label='3 Button Motions Cycle',
#                                         command = lambda s=self: s.cmd.config_mouse('three_button_motions'))

        self.menuBar.addmenuitem('Mouse', 'command', '2 Button Viewing Cycle',
                                         label='2 Button Viewing Cycle',
                                         command = lambda s=self: s.cmd.config_mouse('two_button'))

        self.menuBar.addmenuitem('Mouse', 'command', 'Setup 2 Button Editing Cycle',
                                         label='2 Button Editing Cycle',
                                         command = lambda s=self: s.cmd.config_mouse('two_button_editing'))

        self.menuBar.addmenu('Wizard', 'Task Wizards',tearoff=TRUE)

        self.menuBar.addmenuitem('Wizard', 'command', 'Appearance',
                                         label='Appearance',
                                         command = lambda s=self: s.cmd.do("_ wizard appearance"))

        self.menuBar.addmenuitem('Wizard', 'command', 'Measurement',
                                         label='Measurement',
                                         command = lambda s=self: s.cmd.do("_ wizard measurement"))

        self.menuBar.addmenuitem('Wizard', 'command', 'Mutagenesis',
                                         label='Mutagenesis',
                                         command = lambda s=self: s.cmd.do("_ wizard mutagenesis"))

        self.menuBar.addmenuitem('Wizard', 'command', 'Pair Fitting',
                                         label='Pair Fitting',
                                         command = lambda s=self: s.cmd.do("_ wizard pair_fit"))

        self.menuBar.addmenuitem('Wizard', 'separator', '')
            
        self.menuBar.addmenuitem('Wizard', 'command', 'Density Map Wizard',
                                         label='Density',
                                         command = lambda s=self: s.cmd.do("_ wizard density"))

        self.menuBar.addmenuitem('Wizard', 'command', 'Filter',
                                         label='Filter',
                                         command = lambda s=self: s.cmd.do("_ wizard filter"))


        self.menuBar.addmenuitem('Wizard', 'command', 'Sculpting',
                                         label='Sculpting',
                                         command = lambda s=self: s.cmd.do("_ wizard sculpting"))

        if cleanup.auto_configure()>0:
            self.menuBar.addmenuitem('Wizard', 'separator', '')
        
            self.menuBar.addmenuitem('Wizard', 'command', 'Cleanup',
                                             label='Cleanup',
                                             command = lambda s=self: s.cmd.do("_ wizard cleanup"))

        self.menuBar.addmenuitem('Wizard', 'separator', '')
        
        self.menuBar.addmenuitem('Wizard', 'command', 'Label',
                                         label='Label',
                                         command = lambda s=self: s.cmd.do("_ wizard label"))

        self.menuBar.addmenuitem('Wizard', 'command', 'Charge',
                                         label='Charge',
                                         command = lambda s=self: s.cmd.do("_ wizard charge"))

        self.menuBar.addmenuitem('Wizard', 'separator', '')
        
        self.menuBar.addcascademenu('Wizard', 'Demo', 'Demo',
                                             label='Demo',tearoff=TRUE)

        self.menuBar.addmenuitem('Demo', 'command', 'Representations',
                                         label='Representations',
                                         command = lambda s=self: s.cmd.do(
            "_ replace_wizard demo,reps"))

        self.menuBar.addmenuitem('Demo', 'command', 'Cartoon Ribbons',
                                         label='Cartoon Ribbons',
                                         command = lambda s=self: s.cmd.do(
            "_ replace_wizard demo,cartoon"))

        self.menuBar.addmenuitem('Demo', 'command', 'Roving Detail',
                                         label='Roving Detail',
                                         command = lambda s=self: s.cmd.do(
            "_ replace_wizard demo,roving"))
        
        self.menuBar.addmenuitem('Demo', 'command', 'Roving Density',
                                         label='Roving Density',
                                         command = lambda s=self: s.cmd.do(
            "_ replace_wizard demo,roving_density"))

        self.menuBar.addmenuitem('Demo', 'command', 'Transparency',
                                         label='Transparency',
                                         command = lambda s=self: s.cmd.do(
            "_ replace_wizard demo,trans"))

        self.menuBar.addmenuitem('Demo', 'command', 'Ray Tracing',
                                         label='Ray Tracing',
                                         command = lambda s=self: s.cmd.do(
            "_ replace_wizard demo,ray"))

        self.menuBar.addmenuitem('Demo', 'command', 'Sculpting',
                                         label='Sculpting',
                                         command = lambda s=self: s.cmd.do(
            "_ replace_wizard demo,sculpt"))

        self.menuBar.addmenuitem('Demo', 'command', 'Scripted Animation',
                                         label='Scripted Animation',
                                         command = lambda s=self: s.cmd.do(
            "_ replace_wizard demo,anime"))


        self.menuBar.addmenuitem('Demo', 'command', 'Electrostatics',
                                         label='Electrostatics',
                                         command = lambda s=self: s.cmd.do(
            "_ replace_wizard demo,elec"))


        self.menuBar.addmenuitem('Demo', 'command', 'Compiled Graphics Objects',
                                         label='Compiled Graphics Objects',
                                         command = lambda s=self: s.cmd.do(
            "_ replace_wizard demo,cgo"))

        self.menuBar.addmenuitem('Demo', 'command', 'MolScript/Raster3D Input',
                                         label='Molscript/Raster3D Input',
                                         command = lambda s=self: s.cmd.do(
            "_ replace_wizard demo,raster3d"))

        self.menuBar.addmenuitem('Demo', 'separator', '')
        
        self.menuBar.addmenuitem('Demo', 'command', 'End Demonstration',
                                         label='End Demonstration',
                                         command = lambda s=self: s.cmd.do(
            '_ replace_wizard demo,finish'))

        self.menuBar.addmenu('Plugin', 'Plugin',tearoff=TRUE)      

        self.menuBar.addmenuitem('Plugin', 'command', 'About',
                                         label='About Plugins',
                                         command = lambda s=self: s.about_plugins())

        self.menuBar.addmenuitem('Plugin', 'command', 'Install Plugin',
                                         label='Install Plugin...',
                                         command = lambda s=self: s.app.install_plugin())
        
        self.menuBar.addmenuitem('Plugin', 'separator', '')



    def createInterface(self):

        self.balloon = Pmw.Balloon(self.root)

        self.createMenuBar()

        self.app.menuBar = self.menuBar # to support legacy plugins
        
        self.app.initialize_plugins()
        
        self.createDataArea()

        self.createCommandArea()

#      self.__createAboutBox()
        Pmw.aboutversion(self.appversion)
        Pmw.aboutcopyright(self.copyright)
        Pmw.aboutcontact(
             'For more information, browse to: %s\n or send email to: %s' %\
             (self.contactweb, self.contactemail))
        self.about = Pmw.AboutDialog(self.root, applicationname=self.appname)
        self.about.withdraw()
        
        # Create the parts of the interface

#      self.initPlugins()

        self.buttonArea = Frame(self.root)
        self.buttonArea.pack(side=TOP, anchor=W)
        self.createButtons()

        self.createMessageBar()
        self.messageBar = Pmw.MessageBar(self.commandFrame, entry_width = 40,
             entry_relief='sunken', entry_borderwidth=1) #, labelpos = 'w')
        self.messageBar.pack(side=BOTTOM, anchor=W, fill=X, expand=1)
#        btn_interrupt = self.buttonAdd(self.commandFrame,'Interrupt',lambda s=self: s.cmd.interrupt())
                
        self.balloon.configure(statuscommand = self.messageBar.helpmessage)

        self.createConsole()


    def setup(self):

        # call the parent method
        PMGSkin.setup(self)
        
        # name the application
        self.root.title(self.appname)

        # create the user interface
        self.createInterface()

        # pack the root window
        self.app._hull.pack(side=LEFT, fill=BOTH, expand=YES)

        # and set focus
        self.root.focus_set()

    def __init__(self,app):

        PMGSkin.__init__(self,app)
        self.save_file = ''
        self.cmd = app.pymol.cmd
        self.util = app.pymol.util

def __init__(app):
    app.set_skin(Normal(app))

    
