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
      r = 0
      reps = [ "lines"]
      list = glob(ent_dir)
      for dir in list:
         print dir
         print "deleting"
         cmd.delete("all")
         for file in glob(dir+"/pdb*"):
            print "loading"
            cmd.load(file)
            print "refreshing"
            cmd.refresh()
            print "sleeping"
            time.sleep(0.05)
         cmd.zoom()
         for rep in reps:
            print "showing"
            cmd.show(rep)
            print "refreshing"
            cmd.refresh()
            print "sleeping"
            time.sleep(0.05)
         print "sleeping longer"
         time.sleep(0.3)
         print "next!"
   except:
      traceback.print_exc()

# allow the Tcl/Tk GUI time to get out of the way...
time.sleep(3)

# now start
load()


