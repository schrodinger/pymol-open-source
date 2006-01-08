# -c 
# Check to make sure PyMOL's C secondary structure assignment code
# works on all PDB structures 
      
from glob import glob
import traceback
import threading
import time
from pymol import cmd
import sys, os, os.path
import string

ent_dir = "pdb/*"

def load():
   try:
      r = 0
      list = glob(ent_dir)
      list.sort()
#      list = [ "pdb/vq" ]
      for dir in list:
         sys.__stdout__.write("\n"+dir)
         sys.__stdout__.flush()
         for file in glob(dir+"/pdb*"):
            name = os.path.split(file)[-1]
            name = string.split(name,'.')[0]
            cmd.disable()
            cmd.load(file,name)
            cmd.as("cartoon",name)
            cmd.refresh()
            cmd.dss(name)
            cmd.refresh()
            time.sleep(0.1)
            sys.__stdout__.write(".")
            sys.__stdout__.flush()
         sys.__stdout__.write("("+str(cmd.count_atoms())+")")
         sys.__stdout__.flush()         
         cmd.dss()
         cmd.delete('all')
   except:
      traceback.print_exc()

# allow the Tcl/Tk GUI time to get out of the way...
time.sleep(3)

cmd.set("cartoon_smooth_loops",0)
cmd.set("cartoon_ring_mode",3)
cmd.set("defer_builds_mode",3) 

cmd.feedback('disable','symmetry objectmolecule executive','everything')
# now start
load()


