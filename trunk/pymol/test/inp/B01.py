# 
# basic run through the PDB with wire and cartoon display
#
      
from glob import glob

import threading
import time
import whrandom
from pymol import cmd
import sys, os, os.path

ent_dir = "pdb"

def load():
   list = glob("pdb/*/*")
   l = len(list)
   c = 0
   for file in list:
      c = c + 1
      cmd.delete('pdb')
      cmd.load(file,'pdb')
      cmd.refresh()
      cmd.hide()
      cmd.show('cartoon')
      cmd.color('red','ss h')
      cmd.color('yellow','ss s')
      cmd.orient('pdb')
      sys.__stderr__.write(".")
      sys.__stderr__.flush()
      n = cmd.count_states()
      if n>1:
         for a in range(1,n+1):
            cmd.refresh()
            cmd.frame(a)
            time.sleep(0.025)
         sys.__stderr__.write(" %d of %d"%(c,l))
         sys.__stderr__.write("\n")
         sys.__stderr__.flush()
      else:
         cmd.refresh()
load()


