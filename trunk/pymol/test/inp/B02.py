# 
# basic run through the PDB with ray-tracing
#
      
from glob import glob

import threading
import time
import whrandom
from pymol import cmd
import sys, os, os.path


cmd.viewport(100,75)

ent_dir = "pdb"

r = 0
rep = [ "lines","sticks","spheres","ribbon","cartoon" ]

def load():
   list = glob("pdb/*/*")
   for file in list:
      print file
      cmd.delete('pdb')
      cmd.load(file,'pdb')
      cmd.refresh()
      cmd.hide()
      cmd.show(rep[r])
      r = r + 1
      if r>=len(rep): r=0;
      cmd.color('red','ss h')
      cmd.color('yellow','ss s')
      cmd.orient('pdb')
      sys.__stderr__.write(".")
      sys.__stderr__.flush()
      n = cmd.count_states()
      cmd.ray()
      if n>1:
         sys.__stderr__.write("\n")
         sys.__stderr__.flush()
load()


