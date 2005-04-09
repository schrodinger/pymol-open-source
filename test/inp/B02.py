# 
# basic run through the PDB with ray-tracing
#
      
from glob import glob

import threading
import time
from pymol import cmd
import sys, os, os.path


ent_dir = "pdb"


def load():
   r = 0
   rep = [ "lines","sticks","spheres","dots","ribbon","cartoon" ]
   list = glob("pdb/*/*")
#   while list[0]!="pdb/f8/pdb1f8u":
#      list.pop(0)
   for file in list:
      print file
      cmd.delete('pdb')
      cmd.load(file,'pdb')
      cmd.orient('pdb')
      cmd.color('red','ss h')
      cmd.color('yellow','ss s')
      cmd.hide()
      if cmd.count_atoms()<15000:
         cmd.show(rep[r],"all")
      elif cmd.count_atoms()<50000:
         cmd.show("cartoon","all")	
      else:
         cmd.show("lines","all")	

      r = r + 1
      if r>=len(rep): r=0;
      sys.__stderr__.write(".")
      sys.__stderr__.flush()
      n = cmd.count_states()
      print file
      cmd.ray(160,120)
      if n>1:
         sys.__stderr__.write("\n")
         sys.__stderr__.flush()
      cmd.refresh()
load()


