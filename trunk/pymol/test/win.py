#!/usr/bin/python

import os
import sys
import re
import glob
import string
import traceback
import time

# list of tests to run ...

# prefixes
# C = command line scripts
# G = graphics and internal GUI

standard_prefixes = ['C', 'G' ]

# set module path for generic python stuff

ppath = ''
if os.environ.has_key('PYTHONPATH'):
   ppath = os.environ['PYTHONPATH'] + ":"
ppath = ppath + os.getcwd() + "/../modules"

os.environ['PYTHONPATH'] = ppath

#print ppath

# uniform parameter handling

argv = sys.argv

if len(argv):
   if re.search("\.py$",argv[0]):
      argv = argv[1:]

# ---

pymol = "pymol"
cmmd = "c:\pymolws\pymol.bat "
cmp = "cmp"
ref = "ref"
inp = "inp"
tmp = "tmp"

if not os.path.exists(cmp):
   os.mkdir(cmp)

if not os.path.exists(tmp):
   os.mkdir(tmp)

if len(argv)>1:
   tests = argv
else:
   tests = standard_prefixes

for test in tests:
   flist = glob.glob( inp+"/"+test+"*" )
   cvs = inp+"/CVS"
   if cvs in flist:
      flist.remove(cvs)
   flist.sort()
   for ifil in flist:
      # get options
      f = open(ifil)
      opt = string.strip(f.readline())
      opt = re.sub("^\s*\#","",opt)
      f.close()
      
      tst = re.sub(r".*/|.*\\","",ifil) # get exact test name without suffix
      tst = re.sub(r"\..*","",tst)
      
      print " run_tests: "+tst+"..."
      
      syscmd = cmmd+" -x "+opt+" "+ifil+" > tmp.txt"
      print syscmd
      os.system(syscmd)

      # generate log file
 
      f = open("tmp.txt")
      g = open("cmp\\"+tst+".log","w")

      echo = 0 
      while 1:
         l = f.readline()
         if not l: break
         ll=string.strip(l)
         if ll=='BEGIN-LOG':
            echo = 1
         elif ll=='END-LOG':
            echo = 0
         elif echo:
            g.write(l)

      f.close()
      g.close()

# now compare

      f = open("cmp\\"+tst+".log","r")
      g = open("ref\\"+tst+".log","r")
      while 1:
         lf = f.readline()
         lg = g.readline()
         if (not lf) and not (lg):
            break
         if string.strip(lf)!=string.strip(lg):
            print "<",lf
            print ">",lg

print "done"
time.sleep(60)

# 