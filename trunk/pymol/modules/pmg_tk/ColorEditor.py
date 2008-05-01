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

# this section is devoted to making sure that Tkinter variables which
# correspond to Menu-displayed settings are kept synchronized with
# PyMOL

from Tkinter import *
import tkColorChooser

import Pmw
import string
import copy

class NewColor:
    def __init__(self,app,parent):
        self.app = app
        self.cmd = app.pymol.cmd
        self.parent = parent
        items = []

        self.dialog = Pmw.PromptDialog(self.app.root,title='Create color named...',
                                  buttons = ('Create', 'Cancel'),
                                              defaultbutton='Set',
                                  buttonboxpos=S,
                                  command = self.command)
        self.dialog.geometry("300x120")
        self.entryfield = self.dialog.component('entryfield')
        self.entry = self.entryfield.component('entry')

        self.entryfield.setentry('')
        self.dialog.protocol('WM_DELETE_WINDOW',self.cancel)

        app.my_activate(self.dialog,focus=self.entry)

    def cancel(self,event=None):
        self.command(result='Done')

    def command(self,result=None):
        if result=='Create':
            st = string.strip(self.entry.get())
            if len(st):
                self.parent.update(st)
            self.app.my_deactivate(self.dialog)
            if len(st):
                ColorEdit(self.app,st,self.parent,[1.0,1.0,1.0])
        else:
            self.app.my_deactivate(self.dialog)
        
class ColorEdit:
    def __init__(self,app,name,parent,rgb):
        self.app = app
        self.parent = parent
        self.name = name
        self.cmd = app.pymol.cmd
        
        color = tkColorChooser.Chooser(
            initialcolor='#%02x%02x%02x'%(
            int(rgb[0]*255),int(rgb[1]*255),int(rgb[2]*255)),
            title="Modify color").show()

        if color:
            if color[0]!=None:
                rgb = color[0]
                self.cmd.do("set_color %s,[%5.2f,%5.2f,%5.2f]"%(name,
                         float(rgb[0]/255.0),
                         float(rgb[1]/255.0),
                         float(rgb[2]/255.0)))
                self.cmd.do("recolor")
                
class ColorEditor:

    def __init__(self,app):

        self.app = app
        self.list = []
        self.cmd = app.pymol.cmd
        lst = self.cmd.get_color_indices()
        if lst == None:
            lst = []
        lst.sort()
        for a in lst:
            self.list.append("%-30s"%(a[0]))

        self.index = {}
        c = 0
        for a in lst:
            self.index[a[0]] = a[1]
            c = c + 1
            
        self.dialog = Pmw.SelectionDialog(self.app.root,title="Settings",
                                  buttons = ('New', 'Edit', 'Done'),
                                              defaultbutton='Edit',
                                  scrolledlist_labelpos=N,
                                  label_text='Double click to edit',
                                  scrolledlist_items = self.list,
                                  command = self.command)
        self.dialog.geometry("500x400")
        self.listbox = self.dialog.component('scrolledlist')
        self.listbox.component('listbox').configure(font=app.my_fw_font)
        self.dialog.protocol('WM_DELETE_WINDOW',self.cancel)
        app.my_show(self.dialog)

    def cancel(self,event=None):
        self.command(result='Done')
        
    def update(self,name):
        if not self.index.has_key(name):
            self.listbox.insert(0,"%s"%name)
            self.listbox.selection_clear()         
            self.listbox.selection_set(0)
            for a in self.index.keys():
                self.index[a]=self.index[a]+1
            self.index[name]=0
        else:
            self.listbox.selection_clear()         
            self.listbox.selection_set(self.index[name])
            
    def command(self,result):
        if result=='Done':
            self.app.my_withdraw(self.dialog)
        elif result=='Edit':
            sels = self.dialog.getcurselection()
            if len(sels)!=0:
                color = string.strip(sels[0])
                ColorEdit(self.app,color,self,self.cmd.get_color_tuple(color))
        else:
            NewColor(self.app,self)
            
