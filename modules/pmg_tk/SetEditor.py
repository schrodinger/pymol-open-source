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
import Pmw
from pymol import cmd
import pymol.setting
import string
import copy

class SingleEdit:
   def __init__(self,app,name,parent):
      self.app = app
      self.parent = parent
      self.name = name
      items = []

      prompt = name
      
      self.dialog = Pmw.PromptDialog(self.app.root,title=prompt,
                          buttons = ('Set', 'Cancel'),
                                   defaultbutton='Set',
                          buttonboxpos=S,
                          command = self.command)
      self.dialog.geometry("300x120")
      self.entryfield = self.dialog.component('entryfield')
      self.entry = self.entryfield.component('entry')

      txt = cmd.get_setting_text(name)
      self.entryfield.setentry(txt)
      self.entry.selection_range(0,len(txt))

      app.my_activate(self.dialog,focus=self.entry)

   def command(self,result):
      if result=='Set':
         st = string.strip(self.entry.get())
         if len(st):
            cmd.set(self.name,st,log=1)
         self.parent.update(self.name)
      self.app.my_deactivate(self.dialog)
            
class SetEditor:

   def __init__(self,app):

      self.app = app
      self.list = []
      for a in  pymol.setting.get_index_list():
         self.list.append("%-30s %s"%(pymol.setting._get_name(a),
                    cmd.get_setting_text(a,'',-1)))

      self.index = {}
      c = 0
      for a in pymol.setting.get_name_list():
         self.index[a] = c
         c = c + 1
         
      self.dialog = Pmw.SelectionDialog(self.app.root,title="Settings",
                          buttons = ('Edit', 'Done'),
                                   defaultbutton='Edit',
                          scrolledlist_labelpos=N,
                          label_text='Double click to edit',
                          scrolledlist_items = self.list,
                          command = self.command)
      self.dialog.geometry("500x400")
      self.listbox = self.dialog.component('scrolledlist')
      self.listbox.component('listbox').configure(font=app.my_fw_font)
      
      app.my_show(self.dialog)

   def update(self,name):
      if self.index.has_key(name):
         idx = self.index[name]
         self.listbox.delete(idx,idx)
         self.listbox.insert(idx,"%-30s %s"%(name,
                    cmd.get_setting_text(pymol.setting._get_index(name))))
         self.listbox.selection_set(idx)
         
   def command(self,result):
      if result=='Done':
         self.app.my_withdraw(self.dialog)
      else:
         sels = self.dialog.getcurselection()
         if len(sels)!=0:
            setting = string.strip(sels[0][0:30])
            SingleEdit(self.app,setting,self)

               
