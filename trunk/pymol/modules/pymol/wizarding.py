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

if __name__=='pymol.wizarding':
   
   import imp
   import sys
   import string
   import cmd
   from cmd import _cmd,lock,unlock,Shortcut,QuietException,_raising, \
        _feedback,fb_module,fb_mask

   def wizard(name):
      '''
   DESCRIPTION

      "wizard" launches on of the built-in wizards.  There are special
      Python scripts which work with PyMOL in order to obtain direct user
      interaction and easily peform complicated tasks.

   USAGE

      wizard name

   PYMOL API

      cmd.wizard(string name)

   EXAMPLE

      wizard distance  # launches the distance measurement wizard
   '''
      import wizard
      try:
         if not sys.modules.has_key(name):
            mod_tup = imp.find_module(name,wizard.__path__)
            mod_obj = imp.load_module(name,mod_tup[0],mod_tup[1],mod_tup[2])
         else:
            mod_obj = sys.modules[name]
         if mod_obj:
            oname = string.capitalize(name)
            if hasattr(mod_obj,oname):
               wiz = apply(getattr(mod_obj,oname))
               if wiz:
                  set_wizard(wiz)
                  cmd.do("refresh")
         else:
            print "Error: Sorry, couldn't import the '"+name+"' wizard."         
      except ImportError:
         print "Error: Sorry, couldn't import the '"+name+"' wizard."


   def set_wizard(*arg): # INTERNAL
      r = None
      wiz = None
      if len(arg):
         wiz=arg[0]
      try:
         lock()
         r = _cmd.set_wizard(wiz)
      finally:
         unlock()
      return r

   def refresh_wizard(*arg): # INTERNAL
      r = None
      wiz = None
      if len(arg):
         wiz=arg[0]
      try:
         lock()
         r = _cmd.refresh_wizard()
      finally:
         unlock()
      return r

   def get_wizard(*arg): # INTERNAL
      r = None
      try:
         lock()
         r = _cmd.get_wizard()
      finally:
         unlock()
      return r



