# 
# multi-object run through PDB directories
#
      
from glob import glob
import traceback
import threading
import time
import whrandom
from pymol import cmd
import sys, os, os.path

ent_dir = "pdb/*"

def load():
   try:
      cmd.mplay()
      r = 0
      reps = [ "lines", "sticks", "cartoon" ]
      list = glob(ent_dir)
      for dir in list:
         print dir
         cmd.delete("all")
         for file in glob(dir+"/pdb*"):
            cmd.load(file)
            cmd.refresh()
            cmd.sync()
            time.sleep(0.1)
         cmd.zoom()
         for rep in reps:
            cmd.show(rep)
            cmd.refresh()
            cmd.sync()
            time.sleep(0.1)
         time.sleep(1)
   except:
      traceback.print_exc()
      
#time.sleep(4)
load()


