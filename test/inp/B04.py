# -c 
# I/O check for ability to read & write PDB files in order and with atom names preserved
#
      
from glob import glob
import traceback
import threading
import time
import whrandom
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
      for dir in list:
         for file in glob(dir+"/pdb*"):
            cmd.delete('all')
            cmd.load(file)
            if cmd.count_states()>1: # skip multi-model objects
               break
            if cmd.count_atoms()<0:
               break
            cmd.save(out_file)
            os.system("awk '/^ATOM  /{print substr($0,12,40)};{next;}' < %s > %s"%(out_file,cmp_file))
            os.system("awk '/^ATOM  /{print substr($0,12,40)};{next;}' < %s  > %s"%(file,ref_file))
            os.system("diff %s %s > %s"%(ref_file,cmp_file,dif_file))
                      
            f=open(dif_file)
            l=f.readlines()
            f.close()
            print file
            if len(l):
               for a in l:
                  print a,
               # save it so we have something to look at...
               os.system("/bin/cp -f %s %s_s"%(cmp_file,cmp_file))
               os.system("/bin/cp -f %s %s_s"%(ref_file,ref_file))               
            
   except:
      traceback.print_exc()

# allow the Tcl/Tk GUI time to get out of the way...
time.sleep(3)

cmd.set("retain_order",1)
cmd.feedback('disable','symmetry objectmolecule executive','everything')

# now start
load()


