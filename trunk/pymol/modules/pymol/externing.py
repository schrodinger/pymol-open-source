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
import pymol
import string
import parsing
import thread

from glob import glob
from cmd import _cmd,lock,unlock,Shortcut,QuietException
from cmd import _feedback,fb_module,fb_mask

def cd(dir):
   '''
DESCRIPTION

   "cd" changes the current working directory.

USAGE
   
   cd <path>

SEE ALSO

   pwd, ls, system
   '''
   dir = os.path.expanduser(dir)
   dir = os.path.expandvars(dir)
   os.chdir(dir)  # raises on error
   return 1

def pwd():
   '''
DESCRIPTION

   Print current working directory.

USAGE
   
   pwd

SEE ALSO

   cd, ls, system
   '''
   print os.getcwd()
   return 1

def ls(pattern=None):
   '''
DESCRIPTION

   List contents of the current working directory.

USAGE
   
   ls [pattern]
   dir [pattern]

EXAMPLES

   ls
   ls *.pml

SEE ALSO

   cd, pwd, system   
   '''
   if pattern==None:
      pattern = "*"
   else:
      pattern = os.path.expanduser(pattern)
      pattern = os.path.expandvars(pattern)
   if string.find("*",pattern)<0:
      lst = glob(pattern+"/*")
   else:
      lst = []
   if not len(lst):
      lst = glob(pattern)
   if len(lst):
      lst.sort()
      lst = parsing.list_to_str_list(lst)
      for a in lst:
         print a
   else:
      print " ls: Nothing found.  Is that a valid path?"
   return 1

def system(command,sync=1):
   '''
DESCRIPTION

   "system" executes a command in a subshell under Unix or Windows.

USAGE

   system command 

PYMOL API

   cmd.system(string command,int sync=1)

NOTES

   sync can only be specified from the Python level (not the command language)
   
   if sync is 0, then the command is run in a separate thread (default 1)
   whose identifier is returned in r

   if sync is 1, then the result code from "system" is returned in r
   
SEE ALSO

   ls, cd, pwd
   '''
   if sync:
      r = _cmd.system(str(command),int(sync))
   else:
      r = thread.start_new(_cmd.system,(command,0))
      
   return r # special meaning

def paste(): # INTERNAL
   r=1
   lst = []
   if hasattr(pymol,"machine_get_clipboard"):
      lst = pymol.machine_get_clipboard()
   if len(lst):
      while 1:
         if len(lst[-1]):
            if ord(lst[-1][-1])>32: # trim off final CR
               break;
         else:
            break;
         lst[-1]=lst[-1][:-1]
      _cmd.paste(lst)      
   return r 
