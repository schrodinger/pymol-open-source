# -c 
# Check to make sure PyMOL's C secondary structure assignment code
# works on all PDB structures 
      
from glob import glob
import traceback
import threading
import time
from pymol import cmd
import sys, os, os.path

ent_dir = "pdb/*"

def load():
   try:
      r = 0
      list = glob(ent_dir)
      list.sort()
#      list = [ "pdb/dw" ]
      for dir in list:
         sys.__stdout__.write("\n"+dir)
         sys.__stdout__.flush()
         for file in glob(dir+"/pdb*"):
            name = os.path.split(file)[-1][:-4]
            cmd.load(file,name)
            cmd.hide("everything",name)
            cmd.disable()
            cmd.show("cartoon",name)
            cmd.enable(name)
            cmd.refresh()
            cmd.dss(name)
            cmd.refresh()
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
cmd.feedback('disable','symmetry objectmolecule executive','everything')
# now start
load()


