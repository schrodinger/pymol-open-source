# 
# multi-object run through PDB directories
#
      
from glob import glob

import threading
import time
import whrandom
from pymol import cmd
import sys, os, os.path

ent_dir = "pdb/*"

def load():
   r = 0
   reps = [ "lines", "sticks", "cartoon" ]
   list = glob(ent_dir)
   for dir in list:
      print dir
      cmd.delete("all")
      for file in glob(dir+"/pdb*"):
         cmd.load(file)
         cmd.refresh()
         time.sleep(0.5)
      cmd.zoom()
      for rep in reps:
         cmd.show(rep)
         cmd.refresh()
         time.sleep(0.5)         
load()


