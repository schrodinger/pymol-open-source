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
   r = 0
   rep = [ "lines","sticks","spheres","dots","ribbon","cartoon" ]
   list = glob("pdb/*/*")
   for file in list:
      print file
      cmd.delete('pdb')
      cmd.load(file,'pdb')
      cmd.orient('pdb')
      cmd.color('red','ss h')
      cmd.color('yellow','ss s')
      cmd.hide()
      cmd.show(rep[r],"all")
      r = r + 1
      if r>=len(rep): r=0;
      sys.__stderr__.write(".")
      sys.__stderr__.flush()
      n = cmd.count_states()
      cmd.ray(100,75)
      if n>1:
         sys.__stderr__.write("\n")
         sys.__stderr__.flush()
      cmd.refresh()
load()


