#!/usr/bin/python 

import os
import sys
import re
import glob
import string

# list of tests to run ...

standard_tests = ['cmd', 'chempy']

#-------------

# uniform parameter handling

argv = sys.argv

if len(argv):
   if re.search("\.py$",argv[0]):
      argv = argv[1:]

# ---

pymol = "pymol"
cmp = "cmp"
ref = "ref"
inp = "inp"
if not os.path.exists(cmp):
   os.mkdir(cmp)

if len(argv):
   tests = argv
else:
   tests = standard_tests

for test in tests:
   for ifil in glob.glob( inp+"/"+test+"*" ):
      f = open(ifil) # get command for running test
      cmd = f.readline()
      f.close()
      cmd = string.strip(re.sub(r"^#*","",cmd))
      tst = re.sub(r".*/","",ifil) # get exact test name
      tst = re.sub(r"\..*","",tst)
      print " run_tests: "+tst+"..."
#      print "    "+cmd+" "+ifil+" >& "+cmp+"/"+tst+".log"
      os.system(cmd+" "+ifil+" > "+cmp+"/"+tst+".log")
#      os.system(cmd+" "+ifil)
      for ref_fil in ( glob.glob(ref+"/"+tst+".*")):
         postfx = re.sub(r".*/","",ref_fil)
         cmp_fil = cmp+"/"+postfx
         if os.path.exists(cmp_fil):
            if os.system("diff -q "+ref_fil+" "+cmp_fil+" >/dev/null"):
               print " run_tests: DIFFERS: '"+postfx+"'"
         else:
            print " run_tests: MISSING: '"+postfx+"'"
      for cmp_fil in ( glob.glob(cmp+"/"+tst+".*")):
         postfx = re.sub(r".*/","",cmp_fil)
         ref_fil = ref+"/"+postfx
         if not os.path.exists(ref_fil):
            print " run_tests: no reference version for '"+postfx+"'"

