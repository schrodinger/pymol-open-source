DEBUG = False

import sys
import re
import threading
import os
import time

if True:
    from tkinter import *
    import tkinter.filedialog as tkFileDialog
    import tkinter.messagebox as tkMessageBox
    import tkinter.font as tkFont
    _next_method_name = '__next__'

import Pmw

from pymol.wizard import cleanup

from pmg_tk.Setting import Setting
from pmg_tk.SetEditor import SetEditor
from pmg_tk.ColorEditor import ColorEditor

from pmg_tk.skins import PMGSkin
from .builder import Builder

import pymol._gui

import traceback

root = None

def encode(s):
    # obsolete since removal of py2 
    return s

def split_tk_file_list(pattern):
    filenames = []
    while True:
        pattern = pattern.strip()
        if not pattern:
            break
        sep = None
        if pattern[0] == '{':
            pattern = pattern[1:]
            sep = '}'
        a = pattern.split(sep, 1)
        filenames.append(a[0])
        pattern = a[1] if len(a) == 2 else ''
    return filenames

def asksaveasfilename(*args, **kwargs):
    filename = tkFileDialog.asksaveasfilename(*args, **kwargs)
    return encode(filename)

def askopenfilename(*args, **kwargs):
    filename = tkFileDialog.askopenfilename(*args, **kwargs)
    if not filename:
        return filename
    multiple = kwargs.get('multiple', 0)
    if not multiple:
        filename = [filename]
    elif not isinstance(filename, (list, tuple)):
        filename = split_tk_file_list(filename)
    filename = map(os.path.normpath, filename)
    filename = map(encode, filename)
    filename = list(filename)
    if not multiple:
        return filename[0]
    return filename

def _darwin_browser_open(url):
    os.popen("open "+url,'r').read()
    
def darwin_browser_open(url):
    t = threading.Thread(target=_darwin_browser_open,args=(url,))
    t.setDaemon(1)
    t.start()

def _doAsync(self_cmd,cmmd,dirty=0):
    self_cmd.do(cmmd) # force strict ordering of commands
    if dirty:
        self_cmd.dirty()

def _def_ext(ext): # platform-specific default extension handling
    if sys.platform != 'win32': 
        ext = None # default extensions don't work right under X11/Tcl/Tk
    return ext

class Normal(PMGSkin, pymol._gui.PyMOLDesktopGUI):

    pad = ' ' # extra space in menus
    
    appname        = 'The PyMOL Molecular Graphics System'
    appversion     = '0.0.0.0' # will be set in __init__
    copyright      = ('Copyright (C) 2003-%d\n' % (time.localtime().tm_year,) +
                      'Schrodinger LLC.\n'+
                      'All rights reserved.')
    contactweb     = 'http://www.pymol.org'
    contactemail   = 'sales@schrodinger.com'
    
    # responsible for setup and takedown of the normal skin

    def _inc_fontsize(self, delta, font):
        size = font.cget('size')
        sign = -1 if size < 0 else 1
        size = max(5, abs(size) + delta)
        font.configure(size=size * sign)

    def inc_fontsize(self, delta=1):
        for name in tkFont.names():
            self._inc_fontsize(delta, tkFont.nametofont(name))

    def inc_fontsize_dialog(self):
        dialog = Toplevel(self.root)
        grid = dialog
        kw = {'row': 0, 'sticky': 'w', 'padx': 5, 'pady': 5}
        col = getattr(iter(range(5)), _next_method_name)
        Button(grid, text=' - ', command=lambda: self.inc_fontsize(-1)).grid(column=col(), **kw)
        Button(grid, text=' + ', command=lambda: self.inc_fontsize( 1)).grid(column=col(), **kw)
        Label(grid, text='All GUI Font Sizes').grid(column=col(), **kw)
        kw['row'] = 1
        col = getattr(iter(range(5)), _next_method_name)
        Button(grid, text=' - ', command=lambda: self._inc_fontsize(-1, self.fixedfont)).grid(column=col(), **kw)
        Button(grid, text=' + ', command=lambda: self._inc_fontsize( 1, self.fixedfont)).grid(column=col(), **kw)
        Label(grid, text='Output Font Size').grid(column=col(), **kw)
        dialog.title('GUI Font Size')

    @property
    def initialdir(self):
        '''
        Be in sync with cd/pwd on the console until the first file has been
        browsed, then remember the last directory.
        '''
        return self._initialdir or os.getcwd()

    @initialdir.setter
    def initialdir(self, value):
        self._initialdir = value

    def cd_dialog(self):
        self.cmd.cd(encode(tkFileDialog.askdirectory(
            title="Change Working Directory",
            initialdir=self.initialdir)) or '.', quiet=0)

    def createDataArea(self):
        # Create data area where data entry widgets are placed.
        self.dataArea = self.app.createcomponent('dataarea',
                                             (), None,
                                             Frame, (self.app._hull,), 
                                             relief=SUNKEN, 
                                             bd=1)
        self.dataArea.pack(side=LEFT, fill=BOTH, expand=YES,
                            padx=1, pady=1)

    def destroyDataArea(self):
        self.app.destroycomponent('dataarea')
        
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

    def destroyCommandArea(self):
        self.app.destroycomponent('commandframe')
        
    def createMessageBar(self):
        self.messageBar = Pmw.MessageBar(self.commandFrame, entry_width = 25,
             entry_relief='sunken', entry_borderwidth=1) #, labelpos = 'w')


        self.abortButton=Button(self.commandFrame,
                      text='Abort',highlightthickness=0,
                      command=lambda s=self:self.abort(),padx=0,pady=0)
        self.abortButton.pack(side=RIGHT,fill=BOTH,expand=YES)

        self.abortButton=Button(self.commandFrame,
                                text='Rebuild',highlightthickness=0,
                                #                                state=DISABLED,
                                command=lambda s=self:self.rebuild(),padx=0,pady=0)
        self.abortButton.pack(side=RIGHT,fill=BOTH,expand=YES)

        self.messageBar.pack(side=BOTTOM, anchor=W, fill=X, expand=1)
        self.balloon.configure(statuscommand = self.messageBar.helpmessage)

    def destroyMessageBar(self):
        self.messageBar.destroy()
        
    def get_current_session_file(self):
        session_file = self.cmd.get_setting_text("session_file")        
        session_file = session_file.replace("\\","/") # always use unix-like path separators
        return session_file

    def set_current_session_file(self, session_file):
        session_file = session_file.replace("\\","/") # always use unix-like path separators
        self.cmd.set("session_file",session_file)

    def confirm_quit(self,e=None):
        if self.cmd.get_setting_boolean("session_changed"):
            session_file = self.get_current_session_file()
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


    def buttonAdd(self,frame,text,cmmd):
        newBtn=Button(frame,
                      text=text,highlightthickness=0,
                      command=cmmd,padx=0,pady=0)
        newBtn.pack(side=LEFT,fill=BOTH,expand=YES)
        return newBtn

    def get_view(self):
        self.cmd.get_view(2, quiet=0)
        try:
            str = self.cmd.get_view(3,quiet=1)
            self.root.clipboard_clear()
            self.root.clipboard_append(str)
            self.last_view = str
            self.app.selection_clear()
            self.app.selection_own()
            self.app.selection_handle(lambda a,b,s=self:s.last_view)
            print(" PyMOL: Viewing matrix copied to clipboard.")
        except:
            traceback.print_exc()
        
    def createButtons(self):
        self.buttonArea = Frame(self.root)
        self.buttonArea.pack(side=TOP, anchor=W)
        
        row1 = self.app.createcomponent('row1', (), None,
            Frame,self.commandFrame,bd=0)
        row1.pack(side=TOP,fill=BOTH,expand=YES)
        btn_reset = self.buttonAdd(row1,'Reset',lambda s=self: s.cmd.do("_ reset"))
        btn_reset = self.buttonAdd(row1,'Zoom',lambda s=self: s.cmd.do("_ zoom animate=-1"))
        btn_orient = self.buttonAdd(row1,'Orient',lambda s=self: s.cmd.do("_ orient animate=1"))        
        btn_rtrace = self.buttonAdd(row1,'Draw',lambda s=self: s.cmd.do("_ draw"))        
        btn_rtrace = self.buttonAdd(row1,'Ray',lambda s=self: s.cmd.do("_ ray async=1"))

        row2 = self.app.createcomponent('row2', (), None,
            Frame,self.commandFrame,bd=0)
        row2.pack(side=TOP,fill=BOTH,expand=YES)
        btn_unpick = self.buttonAdd(row2,'Unpick',lambda s=self: s.cmd.do("_ unpick"))
        btn_hidesele = self.buttonAdd(row2,'Deselect', lambda: self.cmd.do("_ deselect"))
        btn_reset = self.buttonAdd(row2,'Rock',lambda s=self: s.cmd.do("_ rock"))
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
        self.volB = self.buttonAdd(row4, 'Volume',
                                    self.newVolumeFrame)
        # initialize disabled
        # self.volB.config(state=DISABLED)

    def newVolumeFrame(self):
        volumes = self.cmd.get_names_of_type("object:volume", public=1)
        if not volumes:
            return
        if len(volumes) == 1:
            self.cmd.volume_panel(volumes[0])
            return
        def callback():
            sels = listbox.getcurselection()
            if sels:
                self.cmd.volume_panel(sels[0])
            window.destroy()
        title = 'Select a volume object'
        window = Toplevel(self.app.root)
        window.title(title)
        listbox = Pmw.ScrolledListBox(window,
                labelpos='nw',
                label_text=title,
                items=volumes,
                selectioncommand=callback)
        listbox.pack(padx=5, pady=5)
        x, y = window.winfo_pointerxy()
        window.geometry('+%d+%d' % (x - 20, y - 20))

    def destroyButtonArea(self):
        self.app.destroycomponent('row1')
        self.app.destroycomponent('row2')
        self.app.destroycomponent('row3')
        self.app.destroycomponent('row4')
        self.buttonArea.destroy()

    def my_show(self,win,center=1):
        win.show()
            
    def my_withdraw(self,win):
        if sys.platform!='linux2':
            win.withdraw()
        else: 
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

    def doAsync(self,cmmd):
        t = threading.Thread(target=_doAsync,args=(self.cmd,cmmd))
        t.setDaemon(1)
        t.start()
        
    def command_get(self):
        return self.command.get()

    def command_set(self, v):
        return self.command.set(v)

    def command_set_cursor(self, i):
        self.entry.icursor(i)

    def dump(self,event):
        print(dir(event))
        print(event.keysym, event.keycode)
        
    def createConsole(self):
        self.command = StringVar()      
        self.lineCount = 0
        self._setup_history()

        self.cmdFrame = Frame(self.dataArea)
        self.buildFrame = Builder(self.app, self.dataArea)
        
        self.toggleFrame(self.cmdFrame,startup=1)

        self.entryFrame = Frame(self.cmdFrame)
        self.entryFrame.pack(side=BOTTOM,expand=NO,fill=X)
        self.entry_label = Label(self.entryFrame, text="PyMOL>", padx=1, pady=1, justify=RIGHT)
        self.entry_label.pack(side=LEFT,expand=NO,fill=X)        
        self.entry = Entry(self.entryFrame, justify=LEFT, width=70,
             textvariable=self.command)
        self.entry.pack(side=LEFT,expand=YES,fill=X)
        self.output = Pmw.ScrolledText(self.cmdFrame)
        self.output.pack(side=TOP, fill=BOTH, expand=YES)      

        self.entry.bind('<Return>', lambda e, s=self:
             (s.doTypedCommand(s.command.get()), s.command.set('')))
        self.entry.bind('<Tab>', lambda e, s=self: s.complete(e))
        self.entry.bind('<Up>', lambda e, s=self: s.back())
        self.entry.bind('<Down>', lambda e, s=self: s.forward())
        self.entry.bind('<Control-Up>', lambda e: self.back_search())
        self.root.protocol("WM_DELETE_WINDOW", lambda s=self: s.confirm_quit())
        
        self.log_file = "log.pml"      

#      self.entry = self.app.createcomponent('entry', (), None,
#                           Entry,
#                           (self.dataArea,),
#                           justify=LEFT,
#                           width=50,
###                           textvariable=self.command)

        text = self.output.component('text')
        self.text = text      

        if sys.platform.startswith('win'):
            self.font = 'lucida console' # only available on windows
            self.my_fw_font=(self.font,8) 
            self.fixedfont.configure(family=self.font, size=self.my_fw_font[1])
        else:
            text.tk.call('tk','scaling',1)
            self.font = 'fixed' # should be available on any X11-based platform
            self.my_fw_font=(self.font,10)
            if sys.platform == 'darwin':
                self.fixedfont.configure(size=11)

        text.configure(width=74)

        self.balloon.bind(self.entry, '''Command Input Area

Get the list of commands by hitting <TAB>

Get the list of arguments for one command with a question mark:
PyMOL> color ?

Read the online help for a command with "help":
PyMOL> help color

Get autocompletion for many arguments by hitting <TAB>
PyMOL> color ye<TAB>    (will autocomplete "yellow")
''')
        
        if self.app.allow_after:
            self.output.after(100,self.update_feedback)
            self.output.after(100,self.update_menus)
            
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
        self.app.bind_all('<Control-F1>',lambda e,s=self: s.cmd.do("cmd._special(1,0,0,2)"))
        self.app.bind_all('<Control-F2>',lambda e,s=self: s.cmd.do("cmd._special(2,0,0,2)"))
        self.app.bind_all('<Control-F3>',lambda e,s=self: s.cmd.do("cmd._special(3,0,0,2)"))
        self.app.bind_all('<Control-F4>',lambda e,s=self: s.cmd.do("cmd._special(4,0,0,2)"))
        self.app.bind_all('<Control-F5>',lambda e,s=self: s.cmd.do("cmd._special(5,0,0,2)"))
        self.app.bind_all('<Control-F6>',lambda e,s=self: s.cmd.do("cmd._special(6,0,0,2)"))
        self.app.bind_all('<Control-F7>',lambda e,s=self: s.cmd.do("cmd._special(7,0,0,2)"))
        self.app.bind_all('<Control-F8>',lambda e,s=self: s.cmd.do("cmd._special(8,0,0,2)"))
        self.app.bind_all('<Control-F9>',lambda e,s=self: s.cmd.do("cmd._special(9,0,0,2)"))
        self.app.bind_all('<Control-F10>',lambda e,s=self: s.cmd.do("cmd._special(10,0,0,2)"))
        self.app.bind_all('<Control-F11>',lambda e,s=self: s.cmd.do("cmd._special(11,0,0,2)"))
        self.app.bind_all('<Control-F12>',lambda e,s=self: s.cmd.do("cmd._special(12,0,0,2)"))                
        
        self.entry.bind('<Prior>',lambda e,s=self: s.cmd.do("cmd._special(104,0,0)"))
        self.entry.bind('<Next>',lambda e,s=self: s.cmd.do("cmd._special(105,0,0)"))
        self.entry.bind('<Control-Prior>',lambda e,s=self: s.cmd.do("cmd._special(104,0,0,2)"))
        self.entry.bind('<Control-Next>',lambda e,s=self: s.cmd.do("cmd._special(105,0,0,2)"))
        self.entry.bind('<Home>',lambda e,s=self: s.cmd.do("cmd._special(106,0,0)"))
        self.entry.bind('<End>',lambda e,s=self: s.cmd.do("cmd._special(107,0,0)"))

    def update_feedback(self):
        feedback = self.cmd._get_feedback(self.cmd)
        if feedback!=None:
            self.text.configure(state='normal')
            for a in feedback:
                self.output.insert(END,"\n")
                self.output.insert(END,a)
                self.output.see(END)
                self.lineCount = self.lineCount + 1
                if self.lineCount > 10000:
                    self.output.delete('0.0','%i.%i' % (self.lineCount-5000,0))
                    self.lineCount=5000
            self.text.configure(state='disabled')
        progress = self.cmd.get_progress()
        if progress>=0.0:
#            self.abortButton.config(state=NORMAL)
            self.messageBar.message("busy","Progress %d%%..."%int(progress*100))
        else:
#            self.abortButton.config(state=DISABLED)            
            self.messageBar.resetmessages("busy")
        if self.app.allow_after:
            if feedback == None: # PyMOL busy, so try more aggressively to get lock
                self.output.after(10,self.update_feedback) # 100X a second                
            else:
                self.output.after(100,self.update_feedback) # 10X a second

    def abort(self):
        self.cmd.interrupt()
#        self.abortButton.config(state=DISABLED)

    def rebuild(self):
        self.doAsync("_ rebuild")
        
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
                frame = self.cmdFrame
        if not startup:
            if frame == self.cmdFrame:
                if self.edit_mode != None:
                    self.cmd.edit_mode(self.edit_mode)
                    self.edit_mode = None
                if self.auto_overlay != None:
                    self.cmd.set("auto_overlay",self.auto_overlay)
                    self.auto_overlay = None
                if self.valence != None:
                    self.cmd.set("valence",self.valence)
            elif frame == self.buildFrame:
                frame.deferred_activate()
                if "Editing" not in self.cmd.get("button_mode_name"):
                    self.cmd.edit_mode(1)
                    self.edit_mode = 0
                self.valence = self.cmd.get("valence")
                self.cmd.set("valence","1")
                self.auto_overlay = self.cmd.get("auto_overlay")
                self.cmd.set("auto_overlay",1)
            
    def update_menus(self):
        self.setting.refresh()

        # 2019-06-25 Disabled because cmd.get_names_of_type() is a blocking
        # command if the API is locked, blocks progress display.
        if False:
            # volume frame is closed, update the button
            if len(self.cmd.get_names_of_type("object:volume",public=1))>0:
                self.volB.config(state=NORMAL)
            else:
                self.volB.config(state=DISABLED)
        # keep calling
        if self.app.allow_after:
            self.output.after(500,self.update_menus) # twice a second

    def file_open(self,tutorial=0):

        if not tutorial:
            initdir = self.initialdir
            ftypes = self.app.getLoadableFileTypes()
        else:
            initdir = os.environ['TUT']
            # only list file extensions that are used for tutorial data
            ftypes = [("Tutorial Data","*.pdb"),]
        if TkVersion>8.3:
            ofile_list = askopenfilename(initialdir = initdir,
                                         filetypes=ftypes,
                                         multiple=1) # new option in Tk 8.4
        else:
            ofile_list = [ askopenfilename(initialdir = initdir,
                                         filetypes=ftypes) ]

        for ofile in ofile_list:
            if len(ofile):
                if not tutorial:
                    self.initialdir = os.path.dirname(ofile)
                try:
                    if ofile[-4:].lower() == '.pse' and ofile != self.save_file:
                        self.save_file = '' # remove ambiguous default
                    self.cmd.do('_ /cmd.load(%s, quiet=0)' % repr(ofile))
                except self.pymol.CmdException:
                    print("Error: unable to open file '%s'"%ofile)

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
            self.initialdir = os.path.dirname(sfile)
            self.log_file = os.path.basename(sfile)
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
            self.initialdir = os.path.dirname(ofile)
            self.log_file = os.path.basename(ofile)
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
            self.initialdir = os.path.dirname(ofile)
            self.log_file = os.path.basename(ofile)
#            os.chdir(self.initialdir)                 
            self.cmd.log_open(ofile,'a')

    def session_save(self):
        self.save_file = self.get_current_session_file()
        if self.save_file!='':
            self.cmd.log("save %s,format=pse\n"%(self.save_file),
                      "cmd.save('%s',format='pse')\n"%(self.save_file))
#            self.cmd.save(self.save_file,"","pse",quiet=0)
#            self.cmd.set("session_changed",0)
            self.cmd.do("_ cmd.save('''%s''','','pse',quiet=0)"%self.save_file) # do this in the main thread to block cmd.quit, etc.
            self.cmd.do("_ cmd.set('session_changed',0)")
            return 1
        else:
            return self.session_save_as()

    def session_save_as(self):
        (self.initialdir, self.save_file) = os.path.split(self.get_current_session_file())
        (save_file, def_ext) = os.path.splitext(self.save_file)
        sfile = asksaveasfilename(defaultextension = _def_ext(def_ext),
                                  initialfile = save_file,  
                                  initialdir = self.initialdir,
                                  filetypes=[
            ("PyMOL Session File","*.pse"),
            ("PyMOL Show File","*.psw"),
            ])
        if len(sfile):
            if re.search(r"\.pse$|\.PSE$|\.psw$|\.PSW$",sfile)==None:
                sfile=sfile+".pse"
            self.initialdir = os.path.dirname(sfile)
            self.cmd.log("save %s,format=pse\n"%(sfile),
                      "cmd.save('%s',format='pse')\n"%(sfile))
#            self.cmd.save(sfile,"",format='pse',quiet=0)
#            self.cmd.set("session_changed",0)
            self.save_file = sfile
#            self.cmd.set("session_file",self.save_file)
            self.set_current_session_file(self.save_file)
            # do this in the main thread to block cmd.quit, etc.
            self.cmd.do("_ cmd.save('''%s''','','pse',quiet=0)"%self.save_file) 
            self.cmd.do("_ cmd.set('session_changed',0)")
            return 1
        else:
            return 0


    def file_save(self):
        """
        File->Save Molecule, now with filtering
        """
        def command(result):
            if result == 'OK':
                self.file_save2(
                        dialog.getcurselection(),
                        multiple_files_option.getvalue(),
                        states_option.getvalue())
            self.my_withdraw(dialog)

        def update_save_listbox():
            lst = self.cmd.get_names('public')
            searchstr = filter_entry.getvalue()

            if searchstr:
                lst = [x for x in lst if searchstr in x]

            dialog.component("scrolledlist").setlist(lst)

        dialog = Pmw.SelectionDialog(self.root,
                                          title="Save",
                                          buttons = ('OK', 'Cancel'),
                                          defaultbutton='OK',
                                          scrolledlist_labelpos=N,
                                          scrolledlist_listbox_selectmode=EXTENDED,
                                          label_text='Which object or selection would you like to save?',
                                          scrolledlist_items = (),  # used to be 'lst'
                                          command = command)

        filter_entry = Pmw.EntryField(dialog.interior(),
                                     labelpos='w',
                                     modifiedcommand=update_save_listbox,
                                     validate=None,
                                     value="",
                                     label_text="Filter:")
        filter_entry.pack(pady=6, fill='x', expand=0, padx=10)

        multiple_files_option = Pmw.RadioSelect( dialog.interior(),
                                                      labelpos='w',
                                                      orient='vertical',
                                                      selectmode='single',
                                                      label_text="Save to...",
                                                      buttontype="radiobutton",
                                                      )
        multiple_files_option.add("one file")
        multiple_files_option.add("multiple files")
        multiple_files_option.invoke("one file")
        multiple_files_option.pack(side='left', pady=8)
                                                      
                                                 
        states_option = Pmw.RadioSelect( dialog.interior(),
                                              labelpos='w',
                                              orient='vertical',
                                              selectmode='single',
                                              label_text='Saved state...',
                                              buttontype="radiobutton"
                                              )
        states_option.add("all")
        states_option.add("global")
        states_option.add("object's current")
        states_option.invoke("global")
        states_option.pack(side='right', pady=8)
                                               
                                                

        # The listbox is created empty.  Fill it now.
        update_save_listbox()

        if len(dialog.component('scrolledlist').get()):
            # set focus on the first item
            listbox = dialog.component('scrolledlist')
            listbox.selection_set(0)

        self.my_show(dialog)
        
    def file_save2(self, sels, multiple_files_flag, state_flag):
        filetypes_save = [
            ("PDB File","*.pdb"),
            ("MOL File","*.mol"),
            ("MOL2 File","*.mol2"),
            ("MMD File","*.mmd"),
            ("PKL File","*.pkl"),
            ("SDF File","*.sdf"),
            ("PDBx/mmCIF","*.cif"),
            ("PQR","*.pqr"),
            ("Maestro","*.mae"),
            ("XYZ","*.xyz"),
        ]
        if True:
            # save N>1 objects to ONE file
            if multiple_files_flag == "one file" and len(sels)>=1:
                        sfile = '_'.join(sels) if len(sels) < 3 else \
                            sels[0] + '-and-%d-more' % (len(sels) - 1)
                        sfile = asksaveasfilename(defaultextension = _def_ext(".pdb"),
                                                  initialfile = sfile,
                                                  initialdir = self.initialdir,
                                                  filetypes=filetypes_save)
                        if len(sfile):
                            # maybe use PDBSTRs for saving multiple files to multiple states
                            self.initialdir = os.path.dirname(sfile)
                            save_sele = ' or '.join(["("+str(x)+")" for x in sels])
                            self.cmd.log("save %s,(%s)\n"%(sfile,save_sele),
                                         "cmd.save('%s','(%s)')\n"%(sfile,save_sele))
                            if state_flag == "all":
                                self.cmd.save(sfile,"(%s)"%save_sele,state=0,quiet=0)
                            elif state_flag == "object's current":
                                ap = 0
                                for sel in sels:
                                    s = int(self.cmd.get("state", str(sel)))
                                    self.cmd.multisave(sfile,str(sel),state=s, quiet=0, append=ap)
                                    ap = 1
                            else:
                                self.cmd.save(sfile,"(%s)"%save_sele,quiet=0)
                            return
            else:
                # save to many files

                for curName in sels:
                    ## print "Result is: ", result
                    ## print "Sels is: ", sels
                    ## print "CurName is: ", curName
                    ## print "State flag is: ", state_flag

                    # The only special case for saving files is when the user selects a multi-state object
                    # and wants to save that to multiple files, each state in one file.
                    doSplit=False
                    if state_flag=='all':
                        stateSave = "0"
                        if len(sels)==1:
#                            print "User wants to split a file"
                            doSplit=True
                    elif state_flag=='global':
                        stateSave = self.cmd.get_state()
                    elif state_flag=="object's current":
                        stateSave = int(self.cmd.get("state",curName))
#                        print "Saving curren't object's state as: ", stateSave
                    else: # default to current global
                        stateSave = "state=", self.cmd.get_state()

                    if True:
                        sfile = asksaveasfilename(defaultextension = _def_ext(".pdb"),
                                                  initialfile = curName,
                                                  initialdir = self.initialdir,
                                                  filetypes = filetypes_save)
                        # now save the file (customizing states as necessary)
#                        print "sfile is: ", sfile

                        if len(sfile):
                            # maybe use PDBSTRs for saving multiple files to multiple states
                            self.initialdir = os.path.dirname(sfile)
                            save_sele = str("("+curName+")")

                            if doSplit:
                                # save each state in "save_sele" to file "sfile" as 'sfile_stateXYZ.pdb'
                                s = self.cmd.count_states(save_sele)
                                for stateSave in range(1,int(s)+1):
                                    save_file = sfile
                                    # _state004
                                    inter = "_state" + str(stateSave).zfill(len(str(s))+1)
                                    # g either MATCHES *.pdb or not.  If so, save, name_stateXYZ.pdb
                                    g = re.search("(.*)(\..*)$", save_file)
                                    if g!=None:
                                        # 1PDB_state004.pdb
                                        save_file = g.groups()[0] + inter + g.groups()[1]
                                    else:
                                        # user entered a file w/o an extension name: eg, '1abc'
                                        # this saves to, '1abc_state00XYZ'
                                        save_file = save_file + inter

                                    self.cmd.log("save %s,(%s)\n"%(save_file,save_sele),
                                                 "cmd.save('%s','(%s)', state='%s')\n"%(save_file,save_sele,stateSave))
                                    self.cmd.save(save_file,"(%s)"%save_sele,state=stateSave,quiet=0)
                            else:
                                save_file = sfile

                                # just save current selection to one file
                                self.cmd.log("save %s,(%s)\n"%(save_file,save_sele),
                                             "cmd.save('%s','(%s)', state='%s')\n"%(save_file,save_sele,stateSave))
                                self.cmd.save(save_file,"(%s)"%save_sele,state=stateSave,quiet=0)


    def edit_pymolrc(self):
        from pmg_tk import TextEditor
        TextEditor.edit_pymolrc(self)

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
            self.__script__ = ofile
            if re.search("\.pym*$|\.PYM*$",ofile):
                self.cmd.do("run "+ofile);      
            else:
                self.cmd.do("@"+ofile);

    def file_save_png(self):
        sfile = asksaveasfilename(defaultextension = _def_ext(".png"),
                                  initialdir = self.initialdir,
                 filetypes=[("PNG File","*.png")])
        if len(sfile):
            self.initialdir = os.path.dirname(sfile)
            self.cmd.log("png %s\n"%sfile,"cmd.png('%s')\n"%sfile)
            self.cmd.png(sfile,quiet=0)

    def file_save_wrl(self):
        sfile = asksaveasfilename(defaultextension = _def_ext(".wrl"),
                                  initialdir = self.initialdir,
                 filetypes=[("VRML 2 WRL File","*.wrl")])
        if len(sfile):
            self.initialdir = os.path.dirname(sfile)
            self.cmd.log("save %s\n"%sfile,"cmd.save('%s')\n"%sfile)
            self.cmd.save(sfile,quiet=0)
            
    def file_save_dae(self):
        sfile = asksaveasfilename(defaultextension = _def_ext(".dae"),
                                  initialdir = self.initialdir,
                 filetypes=[("COLLADA File","*.dae")])
        if len(sfile):
            self.initialdir = os.path.dirname(sfile)
            self.cmd.log("save %s\n"%sfile,"cmd.save('%s')\n"%sfile)
            self.cmd.save(sfile,quiet=0)

    def file_save_pov(self):
        sfile = asksaveasfilename(defaultextension = _def_ext(".pov"),
                                  initialdir = self.initialdir,
                                  filetypes=[("POV File","*.pov")])
        if len(sfile):
            self.initialdir = os.path.dirname(sfile)
            self.cmd.log("save %s\n"%sfile,"cmd.save('%s')\n"%sfile)
            self.cmd.save(sfile,quiet=0)

    def file_save_mpeg(self):
        try:
            from pymol import mpeg_encode
            if not mpeg_encode.validate():
                print("produce-error: Unable to validate pymol.mpeg_encode")
                raise
        except:
            tkMessageBox.showerror("Error",
                "MPEG encoder missing.\nThe FreeMOL add-ons may not be installed.")
            return

        def command(value):
            mQual = int(w_quality.get())
            mode = 'ray' if w_ray.get() else 'draw'
            viewport = int(w_viewport[0].get()), int(w_viewport[1].get())
            dialog.destroy()

            if value != 'OK':
                return

            sfile = asksaveasfilename(defaultextension = _def_ext(".mpg"),
                                      initialdir = self.initialdir,
                                      filetypes=[("MPEG movie file","*.mpg")])
            if len(sfile):
                self.initialdir = os.path.dirname(sfile)
                mQual = self.cmd.get_setting_int("movie_quality")
                self.cmd.log("movie.produce %s,quality=%d,quiet=0\n"%(sfile,mQual),
                             "cmd.movie.produce('''%s''',quality=%d,quiet=0)\n"%(sfile,mQual))
                self.cmd.movie.produce(sfile, mode, quality=mQual, quiet=0,
                        width=viewport[0], height=viewport[1])

        dialog = Pmw.Dialog(title='Movie Settings', buttons=('OK', 'Cancel'),
                defaultbutton='OK', command=command)
        parent = dialog.interior()
        gridkw = {'padx': 5, 'pady': 5, 'sticky': W, 'row': 0}

        Label(parent, text='Encoding Quality (0-100)',).grid(column=0, **gridkw)
        w_quality = Pmw.Counter(parent,
                entryfield_value=self.cmd.get_setting_int("movie_quality"),
                entryfield_validate={'validator': 'integer', 'min': 0, 'max': 100})
        w_quality.grid(column=1, **gridkw)

        gridkw['row'] += 1
        Label(parent, text='Ray Trace Frames').grid(column=0, **gridkw)
        w_ray = BooleanVar(parent, self.cmd.get_setting_boolean('ray_trace_frames'))
        Checkbutton(parent, variable=w_ray).grid(column=1, **gridkw)

        w_viewport = []
        for text, value in zip(('Width', 'Height'), self.cmd.get_viewport()):
            gridkw['row'] += 1
            Label(parent, text=text + ' (pixels)').grid(column=0, **gridkw)
            w = Pmw.Counter(parent, entryfield_value=value, entryfield_validate={'validator': 'integer', 'min': 0})
            w.grid(column=1, **gridkw)
            w_viewport.append(w)

    def file_save_mpng(self):
        sfile = asksaveasfilename(initialdir = self.initialdir,
                                  filetypes=[("Numbered PNG Files","*.png")])
        if len(sfile):
            self.initialdir = os.path.dirname(sfile)
            self.cmd.log("mpng %s\n"%sfile,"cmd.mpng('%s')\n"%sfile)         
            self.cmd.mpng(sfile,modal=-1)

    def cat_terms(self):
        for path in [ "$PYMOL_PATH/LICENSE.txt", "$PYMOL_PATH/LICENSE.TXT", "$PYMOL_PATH/LICENSE" ]:
            path = self.pymol.cmd.exp_path(path)
            if os.path.exists(path):
                print(open(path).read().strip())
                return 
        print(" Error: no license terms found.")

    def toggleClickThrough(self, toggle):
        if toggle:
            os.system(
            "defaults write com.apple.x11 wm_click_through -bool true")
            os.system(
            "defaults write org.x.X11 wm_click_through -bool true")
            os.system(
            "defaults write com.apple.x11 wm_ffm -bool true")
            os.system(
            "defaults write org.x.X11 wm_ffm -bool true")
            print("Enabled wm_click_through and wm_ffm.", end=' ')
        else:
            os.system(
            "defaults write com.apple.x11 wm_click_through -bool false")
            os.system(
            "defaults write org.x.X11 wm_click_through -bool false")
            os.system(
            "defaults write com.apple.x11 wm_ffm -bool false")
            os.system(
            "defaults write org.x.X11 wm_ffm -bool false")
            print("Disabled wm_click_through and wm_ffm.", end=' ')
        print("Please restart X11.")

    def createMenuBar(self):
        self.menuBar = Pmw.MenuBar(self.root, balloon=self.balloon,
                                   hull_relief=RAISED, hull_borderwidth=1) 
        self.menuBar.pack(fill=X)

        addmenuitem = self.menuBar.addmenuitem
        addcascademenu = self.menuBar.addcascademenu

        self.setting = Setting(self.app)

        def _addmenu(data, parent):
            for item in data:
                if item[0] == 'separator':
                    addmenuitem(parent, 'separator', '')
                elif item[0] == 'menu':
                    label = item[1]
                    menulabel = parent + '/' + label
                    self.menuBar.addcascademenu(parent, menulabel,
                        label, label=label, tearoff=FALSE)
                    _addmenu(item[2], menulabel)
                elif item[0] == 'command':
                    label = item[1]
                    command = item[2]
                    if command is None:
                        if DEBUG:
                            print('warning: skipping', label, parent)
                    else:
                        if isinstance(command, str):
                            command = lambda c=command: self.cmd.do(c)
                        addmenuitem(parent, 'command', label, label=label, command=command)
                elif item[0] == 'check':
                    label = item[1]
                    var = getattr(self.setting, item[2])
                    if len(item) > 4:
                        addmenuitem(parent, 'checkbutton', label, label=label, variable=var, onvalue=item[3], offvalue=item[4])
                    else:
                        addmenuitem(parent, 'checkbutton', label, label=label, variable=var)
                elif item[0] == 'radio':
                    label = item[1]
                    var = getattr(self.setting, item[2])
                    value = item[3]
                    addmenuitem(parent, 'radiobutton', label=label, value=value, variable=var)
                elif DEBUG:
                    print('error:', item)

        for _, label, data in self.get_menudata():
            assert _ == 'menu'
            self.menuBar.addmenu(label, label, tearoff=TRUE)
            _addmenu(data, label)

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

#        self.toggleBalloonVar = IntVar()
#        self.toggleBalloonVar.set(0)

#      self.menuBar.addmenuitem('Help', 'separator', '')
        
#      self.menuBar.addmenuitem('Help', 'checkbutton',
#                         'Toggle balloon help',
#                         label='Balloon help',
#                        variable = self.toggleBalloonVar,
#                        command=self.toggleBalloon)

        if sys.platform == 'win32' and self.app.pymol.invocation.options.incentive_product:
                self.menuBar.addmenuitem('Edit', 'separator', '')

                self.menuBar.addmenuitem('Edit', 'command',
                                     'Copy Image',
                                     label='Copy Image to Clipboard',
                                     command = lambda s=self:s.cmd.copy_image(quiet=0))

                self.menuBar.addmenuitem('Edit', 'checkbutton',
                                 'Auto-Copy Images',
                                 label='Auto-Copy Images',
                                 variable = self.setting.auto_copy_images,
                                 )

        self.menuBar.addmenuitem('Edit', 'separator', '')

        self.menuBar.addmenuitem('Edit', 'command',
                                 'To Copy Text: Use Ctrl-C in TclTk GUI',
                                 label='To copy text use Ctrl-C in the TclTk GUI',
                                         state='disabled',
                                command =  None)

        self.menuBar.addmenuitem('Edit', 'command',
                                 'To Paste Text, Use Ctrl-V in TclTk GUI',
                                 label='To paste text use Ctrl-V in the TckTk GUI',
                                         state='disabled',                               
                                command =  None)

        if sys.platform == 'darwin':
            self.menuBar.addmenuitem('Mouse', 'separator', '')

            self.menuBar.addcascademenu('Mouse', 'MacX11Focus', 'Mac OS X11',
                                        label='Mac OS X11')
            self.menuBar.addmenuitem('MacX11Focus', 'command',
                                     'Enable Click Through',
                                     label='Enable Click Through',
                                     command = lambda s=self: 
                                     s.toggleClickThrough(1))

            self.menuBar.addmenuitem('MacX11Focus', 'command',
                                     'Disable Click Through',
                                     label='Disable Click Through',
                                     command = lambda s=self: 
                                     s.toggleClickThrough(0))

        # hook up scene menu updates
        index = self.pymol.setting.index_dict.get('scenes_changed')
        self.setting.active_dict[index] = self.update_scene_menu

    def settings_edit_all_dialog(self):
        SetEditor(self)

    def edit_colors_dialog(self):
        ColorEditor(self)

    def update_scene_menu(self):
        scene_list = self.cmd.get_scene_list()
        for action in ['recall', 'clear']:
            parent = 'Scene/' + action.capitalize()
            self.menuBar.deletemenuitems(parent, 0, 999)
            for k in scene_list:
                self.menuBar.addmenuitem(parent, 'command', k, label=k,
                        command=lambda k=k, a=action: self.cmd.scene(k, a))
        parent = 'Scene/Store'
        self.menuBar.deletemenuitems(parent, 0, 11)
        for i in range(12):
            k = 'F' + str(i + 1)
            self.scene_F_keys[i].set(1 if k in scene_list else 0)
            self.menuBar.addmenuitem(parent, 'checkbutton', k, label=k,
                    variable=self.scene_F_keys[i],
                    command=lambda k=k: self.cmd.scene(k, 'store'))

    def show_about(self):
        Pmw.aboutversion(self.appversion)
        Pmw.aboutcopyright(self.copyright)
        Pmw.aboutcontact(
             'For more information, browse to: %s\n or send email to: %s' %\
             (self.contactweb, self.contactemail))
        self.about = Pmw.AboutDialog(self.root, applicationname=self.appname)
        self.my_activate(self.about)
        self.about.withdraw()
        
    def createInterface(self):

        self.balloon = Pmw.Balloon(self.root)

        self.createMenuBar()

        self.app.menuBar = self.menuBar # to support legacy plugins    
        
        self.app.initializePlugins()
        
        self.createDataArea()

        self.createCommandArea()

        self.createButtons()

        self.createMessageBar()

        self.createConsole()

    def setup(self):

        # call the parent method
        PMGSkin.setup(self)
        
        # name the application
        self.root.title(self.appname)

        # create the user interface
        self.createInterface()

        # pack the root window
        self.app._hull.pack(side=LEFT, fill=BOTH, expand=YES, anchor=CENTER)

        # and set focus
        if hasattr(self,'entry'): self.entry.focus_set()

    def takedown(self):
        self.destroyMessageBar()
        self.destroyDataArea()
        self.destroyCommandArea()
        self.destroyButtonArea()
        self.balloon.destroy()
        self.menuBar.destroy()
        
    def __init__(self,app):
        global root
        root = app.root

        PMGSkin.__init__(self,app)
        Normal.appversion = app.pymol.cmd.get_version()[0]
        Normal.appversion += " Incentive Product" \
                if app.pymol.invocation.options.incentive_product else \
                " Open-Source"
        self.app = app
        self.save_file = ''
        self.cmd = app.pymol.cmd
        self.util = app.pymol.util
        self.movie_command = None
        self.movie_start = 1
        self.auto_overlay = None
        self.edit_mode = None
        self.valence = None
        self._initialdir = ''
        self.fixedfont = tkFont.nametofont('TkFixedFont')
        self.scene_F_keys = [IntVar(root) for _ in range(12)]

def __init__(app):
    return Normal(app)

    
