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

def get_session():
   '''
NOT YET IMPLEMENTED
'''
   r = None
   try:
      lock()   
      r = _cmd.get_session(1)
   finally:
      unlock()
   return r

def png(filename):
   '''
DESCRIPTION
  
   "png" writes a png format image file of the current image to disk.
   
USAGE
 
   png filename
 
PYMOL API
 
   cmd.png( string file )
   '''
   if thread.get_ident() ==pymol.glutThread:
      r = cmd._png(str(filename))
   else:
      r = _cmd.do("cmd._png('"+str(filename)+"')")
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

def save(filename,selection='(all)',state=0,format=''):
   '''
DESCRIPTION
  
   "save" writes selected atoms to a file.  The file format is
   autodetected if the extesion is ".pdb" or ".pkl"
 
USAGE
 
   save file [,(selection) [,state [,format]] ]
 
PYMOL API
  
   cmd.save(file, selection, state, format)

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
         print " Save: wrote \""+filename+"\"."
   elif format=='pkl': # python binary
      io.pkl.toFile(cmd.get_model(selection,state),filename)
      print " Save: wrote \""+filename+"\"."
   elif format=='pkla': # ascii override
      io.pkl.toFile(cmd.get_model(selection),filename,bin=0)
      print " Save: wrote \""+filename+"\"."
   elif format=='pse': # PyMOL session
      io.pkl.toFile(cmd.get_session(),filename)
      print " Save: wrote \""+filename+"\"."
   elif format=='mmod': # macromodel
      io.mmd.toFile(cmd.get_model(selection),filename)
      print " Save: wrote \""+filename+"\"."
   elif format=='mol': 
      io.mol.toFile(cmd.get_model(selection),filename)
      print " Save: wrote \""+filename+"\"."
   return r

