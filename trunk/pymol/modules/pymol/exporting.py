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

if __name__=='pymol.exporting':
   import os
   import thread
   import selector
   import string
   import re

   import pymol
   import cmd
   from cmd import _cmd,lock,unlock,Shortcut,QuietException
   from chempy import io
   from cmd import _feedback,fb_module,fb_mask
   import traceback

   def get_pdbstr(selection="all", state=0):
      '''
DESCRIPTION

   "get_pdbstr" in an API-only function which returns a pdb
   corresponding to the atoms in the selection provided and that are
   present in the indicated state

PYMOL API ONLY

   cmd.get_pdbstr( string selection="all", int state=0 )

NOTES

   "state" is a 1-based state index for the object.

   if state is zero, then current state is used.
   
   '''
      r = None
      try:
         lock()   
         r = _cmd.get_pdb(str(selection),int(state)-1)
      finally:
         unlock()
      return r
   
   def get_session():
      session = {}
      r = 1
      for a in pymol._session_save_tasks:
         if a==None:
            try:
               lock()
               r = _cmd.get_session(session)
            finally:
               unlock()
         else:
            try:
               r = apply(a,(session,))
            except:
               traceback.print_exc()
               print "Error: An error occurred when trying to generate session."
               print "Error: The resulting session file may be incomplete."
      return session

   def png(filename,quiet=1):
      '''
DESCRIPTION

   "png" writes a png format image file of the current image to disk.

USAGE

   png filename

PYMOL API

   cmd.png( string file )
      '''
      if thread.get_ident() ==pymol.glutThread:
         r = cmd._png(str(filename),int(quiet))
      else:
         r = _cmd.do("cmd._png('"+str(filename)+"')",0)
      return r

   def export_coords(obj,state): # experimental
      r = None
      try:
         lock()   
         r = _cmd.export_coords(str(obj),int(state)-1)
      finally:
         unlock()
      return r

   def multisave(filename,object,state=0):
      r = 1
      try:
         lock()
         _cmd.multisave(str(filename),str(object),int(state)-1,0)
      finally:
         unlock()
      return r

   def save(filename,selection='(all)',state=0,format='',quiet=1):
      '''
DESCRIPTION

   "save" writes selected atoms to a file.  The file format is
   autodetected if the extesion is ".pdb", ".pse", ".mol", ".mmod", or
   ".pkl"

   Note that if the file extension ends in ".pse" (PyMOL Session), the
   complete PyMOL state is always saved to the file (the selection and
   state parameters are thus ignored).

USAGE

   save file [,(selection) [,state [,format]] ]

PYMOL API

   cmd.save(file, selection, state, format)

NOTES

   When saving a session file, then "state" has no effect.
   When state = 0 (default), only the current state is written.
   When state = -1, then a multi-state output file is written (PDB only).
   
SEE ALSO

   load, get_model
      '''
      # preprocess selection
      selection = selector.process(selection)
      #   
      r = 1
      if format=='':
         format = 'pdb'
         lc_filename=string.lower(filename)
         if re.search("\.pdb$|\.ent$",lc_filename):
            format = 'pdb'
         elif re.search("\.mol$",lc_filename):
            format = 'mol'
         elif re.search("\.pkl$",lc_filename):
            format = 'pkl'
         elif re.search("\.pkl$",lc_filename):
            format = 'pkla'
         elif re.search("\.mmd$",lc_filename):
            format = 'mmod'
         elif re.search("\.mmod$",lc_filename):
            format = 'mmod'
         elif re.search("\.pmo$",lc_filename):
            format = 'pmo'
         elif re.search("\.png$",lc_filename):
            format = 'png'
         elif re.search("\.pse$",lc_filename):
            format = 'pse'
      else:
         format = str(format)
      filename = os.path.expanduser(filename)
      filename = os.path.expandvars(filename)
      if format=='pdb':
         f=open(filename,"w")
         if f:
            try:
               lock()
               st = _cmd.get_pdb("("+str(selection)+")",int(state)-1)
            finally:
               unlock()
               f.write(st)
               f.close()
            r = None
            if not quiet:
               print " Save: wrote \""+filename+"\"."
      elif format=='pkl': # python binary
         io.pkl.toFile(cmd.get_model(selection,state),filename)
         if not quiet:
            print " Save: wrote \""+filename+"\"."
      elif format=='pkla': # ascii override
         io.pkl.toFile(cmd.get_model(selection),filename,bin=0)
         if not quiet:
            print " Save: wrote \""+filename+"\"."
      elif format=='pse': # PyMOL session
         io.pkl.toFile(cmd.get_session(),filename)
         if not quiet:
            print " Save: wrote \""+filename+"\"."
      elif format=='mmod': # macromodel
         io.mmd.toFile(cmd.get_model(selection),filename)
         if not quiet:
            print " Save: wrote \""+filename+"\"."
      elif format=='mol': 
         io.mol.toFile(cmd.get_model(selection),filename)
         if not quiet:
            print " Save: wrote \""+filename+"\"."
      elif format=='png':
         cmd.png(filename,quiet=quiet)
      return r

