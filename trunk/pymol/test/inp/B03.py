# 
# basic run through the PDB with ray-tracing
#
      
from glob import glob

import threading
import time
import whrandom
from pymol import cmd
import sys, os, os.path


ent_dir = "pdb"


def load():
   cmd.set("valence")
   r = 0
   list = glob("pdb/*/*")
#   while list[0]!="pdb/f8/pdb1f8u":
#      list.pop(0)
   for file in list:
      print file
      cmd.delete('pdb')
      cmd.load(file,'pdb')
      cmd.orient('pdb')
      cmd.refresh()
      cmd.hide()
      cmd.show("sticks")
      sys.__stderr__.write(".")
      sys.__stderr__.flush()
      n = cmd.count_states()
      cmd.refresh()
      if n>1:
         sys.__stderr__.write(file+"\n")
         sys.__stderr__.flush()
         for a in range(1,n+1):
            cmd.forward()
            cmd.refresh()
load()


