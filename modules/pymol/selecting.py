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
#-* 
#-*
#Z* -------------------------------------------------------------------

import selector

import cmd

from cmd import _cmd,lock,unlock,Shortcut,QuietException
from cmd import _feedback,fb_module,fb_mask

def deselect():
   '''
DESCRIPTION
  
   "deselect" disables any and all visible selections
 
USAGE
 
   deselect
 
PYMOL API
  
   cmd.deselect()
   '''
   arg = cmd.get_names("selections")
   for a in arg:
      cmd.disable(a)

def select(name,selection="",quiet=0,show=0):
   '''
DESCRIPTION
  
   "select" creates a named selection from an atom selection.
 
USAGE
 
   select (selection)
   select name, (selection)
   select name = (selection)            # (DEPRECATED)
 
PYMOL API
  
   cmd.select(string name, string selection)
 
EXAMPLES 

   select near , (ll expand 8)
   select near , (ll expand 8)
   select bb, (name ca,n,c,o )

NOTES

   'help selections' for more information about selections.
   '''   
   try:
      lock()
      if selection=="":
         sel_cnt = _cmd.get("sel_counter") + 1.0
         _cmd.legacy_set("sel_counter","%1.0f" % sel_cnt)
         selection = name
         name = "sel%02.0f" % sel_cnt
      else:
         name = name
      # preprocess selection (note: inside TRY)
      selection = selector.process(selection)
      #
      r = _cmd.select(str(name),str(selection),int(quiet))
      if r and show:
         r = _cmd.onoff(str(name),1);
   finally:
      unlock()
   return r

def indicate(selection="(all)"):
   '''
DESCRIPTION
  
   "indicate" shows a visual representation of an atom selection.
 
USAGE
 
   indicate (selection)
 
PYMOL API
  
   cmd.count(string selection)
 
   '''
   # preprocess selection
   selection = selector.process(selection)
   #      
   try:
      lock()   
      r = _cmd.select("_indicate","("+str(selection)+")",1)
      cmd.enable("_indicate")
   finally:
      unlock()
   return r
