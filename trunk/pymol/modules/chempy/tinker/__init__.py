
import os
import shutil
import glob
import re
import string

prefix = "tinker_run"

def run(command,in_prefix,out_prefix,tokens):
   for a in glob.glob(prefix+".*"):
      os.unlink(a)
   for src in glob.glob(in_prefix+".*"):
      dst = string.split(src,'.')
      dst = prefix+'.'+dst[len(dst)-1]
      shutil.copyfile(src,dst)
   pipe = os.popen(bin_path+command,"w")
   if not pipe:
      print "Error: can't run tinker!!!"
      raise RunError
   for a in tokens:
      pipe.write(a+"\n")
   pipe.close()
   for src in glob.glob(prefix+".*_2"):
      dst = string.replace(src,'_2','')
      if os.path.exists(dst):
         os.unlink(dst)
      os.rename(src,dst)      
   for src in glob.glob(prefix+".*"):
      dst = string.split(src,'.')
      dst = out_prefix+'.'+dst[len(dst)-1]
      if os.path.exists(dst):
         os.unlink(dst)
      os.rename(src,dst)
   
if os.environ.has_key('TINKER_PATH'):
   base = os.environ['TINKER_PATH']
   bin_path = base + '/bin/'
   params_path = base + '/params/'
   if os.environ.has_key('PYMOL_PATH'):
      pymol_path = os.environ['PYMOL_PATH']
      test_path = pymol_path + '/modules/chempy/tinker/'
      if os.path.exists(test_path):
         params_path = test_path
else:
   base = ''
   bin_path = ''
   params_path = ''

