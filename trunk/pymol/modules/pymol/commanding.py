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

if __name__=='pymol.commanding':

   import thread
   import string
   import re
   import os
   import time
   import threading
   import traceback

   import cmd
   import pymol
   from cmd import _cmd,lock,unlock,Shortcut,QuietException
   from cmd import _feedback,fb_module,fb_mask,is_list

   def resume(fname):
      if os.path.exists(fname):
         if(re.search(r"\.py$|\.PY$|\.pym$|.PYM$",fname)):
            cmd.do("run %s"%fname)
         else:
            cmd.do("@%s"%fname)
      cmd.do("log_open %s,a"%fname)

   def log_open(fname='log.pml',mode='w'):
      try:
         try:
            if pymol._log_file!=None:
               pymol._log_file.close()
         except:
            pass
         pymol._log_file = open(fname,mode)
         if _feedback(fb_module.cmd,fb_mask.details): # redundant
            if mode!='a':
               print " Cmd: logging to '%s'."%fname
            else:
               print " Cmd: appending to '%s'."%fname            
         if mode=='a':
            pymol._log_file.write("\n") # always start on a new line
         if(re.search(r"\.py$|\.PY$|\.pym$|.PYM$",fname)):
            cmd.set("logging",2,quiet=1)
         else:
            cmd.set("logging",1,quiet=1)
      except:
         print"Error: unable to open log file '%s'"%fname
         pymol._log_file = None
         cmd.set("logging",0,quiet=1)
         traceback.print_exc()

   def log(text,alt_text=None):
      if pymol._log_file!=None:
         mode = cmd.get_setting_legacy("logging")
         if mode:
            if mode==1:
               pymol._log_file.write(text)
            elif mode==2:
               if alt_text!=None:
                  pymol._log_file.write(alt_text)
               else:
                  pymol._log_file.write("cmd.do('''%s''')\n"%string.strip(text))
            pymol._log_file.flush()

   def log_close():
      if pymol._log_file!=None:
         pymol._log_file.close()
         cmd.set("logging",0,quiet=1)
         if _feedback(fb_module.cmd,fb_mask.details): # redundant
            print " Cmd: log closed."

   def cls(): 
      '''
DESCRIPTION

   "cls" clears the output buffer.

USAGE

   cls
   '''
      r = None
      try:
         lock()
         r = _cmd.cls()
      finally:
         unlock()
      return r

   def splash(mode=0):
      '''
DESCRIPTION

   "splash" shows the splash screen information.

USAGE

   splash
      '''
      r = None
      mode = int(mode)
      if mode == 1: # just show PNG
         png_path = os.path.expandvars("$PYMOL_PATH/data/pymol/splash.png")
         if os.path.exists(png_path):
            cmd.do("_ cmd.load_png('%s',0,quiet=1)"%png_path)
            r = 1
      else:
         if cmd.get_setting_legacy("internal_feedback")>0.1:
            cmd.set("text","1",quiet=1)
         print
         try:
            lock()
            r = _cmd.splash()
         finally:
            unlock()
      return r

   def reinitialize():
      '''
DESCRIPTION

   "reinitialize" reinitializes PyMOL

USAGE

   reinitialize
      '''
      try:
         lock()
         r = _cmd.reinitialize()
      finally:
         unlock()
      return r


   def sync(timeout=1.0,poll=0.05):
      '''
DESCRIPTION

   "sync" is an API-only function which waits until all current
   commmands have been executed before returning.  A timeout
   can be used to insure that this command eventually returns.

PYMOL API

   cmd.sync(float timeout=1.0,float poll=0.05)

SEE ALSO

   frame
   '''
      now = time.time()
      timeout = float(timeout)
      poll = float(poll)
      if _cmd.wait_queue(): # commands waiting to be executed?
         while 1:
            e = threading.Event() # using this for portable delay
            e.wait(poll)
            del e
            if not _cmd.wait_queue():
               break
            if (time.time()-now)>timeout:
               break


   def do(commands):
      # WARNING: don't call this routine if you already have the API lock
      '''
DESCRIPTION

   "do" makes it possible for python programs to issue simple PyMOL
   commands as if they were entered on the command line.

PYMOL API

   cmd.do( commands )

USAGE (PYTHON)

   from pymol import cmd
   cmd.do("load file.pdb")
      '''
      r = None
      if is_list(commands):
         cmmd_list = commands
      else:
         cmmd_list = [ commands ]
      for cmmd in cmmd_list:
         lst = string.split(cmmd,"\n")
         if len(lst)<2:
            for a in lst:
               if(len(a)):
                  try:
                     lock()
                     r = None
                     r = _cmd.do(a)
                  finally:
                     unlock()
         else:
            defer = cmd.get_setting_legacy("defer_updates")
            cmd.set('defer_updates',1)
            for a in lst:
               if(len(a)):
                  try:
                     lock()
                     r = None
                     r = _cmd.do(a)
                  finally:
                     unlock()
            cmd.set('defer_updates',defer)
      return r


   def quit():
      '''
DESCRIPTION

   "quit" terminates the program. 

USAGE

   quit

PYMOL API

   cmd.quit()
      '''
      if thread.get_ident() == pymol.glutThread:
         try: 
            lock()
            r = _cmd.quit()
         finally:
            unlock()
      else:
         try:
            lock()
            r = _cmd.do("_ cmd.quit()")
            thread.exit()
         finally:
            unlock()
      return r

   def delete(name):
      '''
DESCRIPTION

   "delete" removes an object or a selection. 

USAGE

   delete name  
   delete all   # deletes all objects

   name = name of object or selection

PYMOL API

   cmd.delete (string name = object-or-selection-name )

SEE ALSO

   remove
      '''
      try:
         lock()   
         r = _cmd.delete(str(name))
      finally:
         unlock()
      return r
