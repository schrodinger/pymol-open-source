# -c 
# I/O check for ability to read & write PDB files in order and with atom names preserved
#
      
from glob import glob
import traceback
import threading
import time
from pymol import cmd
import sys, os, os.path

ent_dir = "pdb/*"

ref_file = "tmp/ref.pdb"
out_file = "tmp/out.pdb"
cmp_file = "tmp/cmp.pdb"
dif_file = "tmp/dif.pdb"

def load():
   try:
      r = 0
      list = glob(ent_dir)
      list.sort()
      list.reverse()
#      list = ["pdb/bn/"]
      for dir in list:
         sys.__stdout__.write("\n"+dir)
         sys.__stdout__.flush()
         for file in glob(dir+"/pdb*"):
            cmd.delete('all')
            cmd.load(file)
            if cmd.count_states()>1: # skip multi-model objects
               break
            if cmd.count_atoms()<0:
               break
            cmd.save(out_file)
            os.system("awk '/^ATOM  /{print substr($0,12,16),\" \"};{next;}' < %s | sed 's/  *//g' > %s"%(out_file,cmp_file))
            os.system("awk '/^ATOM  /{print substr($0,12,16),\" \"};{next;}' < %s | sed 's/  *//g' > %s"%(file,ref_file))
            os.system("diff %s %s > %s"%(ref_file,cmp_file,dif_file))
                      
            f=open(dif_file)
            l=f.readlines()
            f.close()
            if len(l):
               print
               for a in l:
                  print a,
               # save it so we have something to look at...
               os.system("/bin/cp -f %s %s_s"%(cmp_file,cmp_file))
               os.system("/bin/cp -f %s %s_s"%(ref_file,ref_file))               
               os.system("/bin/cp -f %s %s_src"%(out_file,cmp_file))
               os.system("/bin/cp -f %s %s_src"%(file,ref_file))
               print file
            else:
               sys.__stdout__.write(".")
               sys.__stdout__.flush()
   except:
      traceback.print_exc()

# allow the Tcl/Tk GUI time to get out of the way...
time.sleep(3)

cmd.set("retain_order",1)
cmd.feedback('disable','symmetry objectmolecule executive','everything')

# now start
load()


