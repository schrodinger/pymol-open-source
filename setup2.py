#
# This script only applies if you are performing a Python Distutils-based
# installation of PyMOL.
#
# Please run this script after using Distutils "setup.py":
#   python setup.py build
#   python setup.py install
#   python install.py

import os
import re
import sys
from distutils import dir_util,file_util

if sys.platform=='win32':
   launch_script = "pymol.bat"
elif sys.platform=='cygwin':
   launch_script = "pymol"
else:
   launch_script = "pymol.com"
   
try:
   pymol_launch = 3
   import pymol
   
   if not (os.path.exists("data") and os.path.exists("LICENSE")):
      print ' Please run "python setup2.py" from the PyMOL source directory.'
   else:
      pymol_file = sys.modules['pymol'].__file__
      if pymol_file[-4:]==".pyc":
         pymol_file = pymol_file[0:-1]
      pymol_path = re.sub(r"[\/\\][^\/\\]*$","/pymol_path",pymol_file)

      print pymol_path
      # Create PYMOL_PATH directory
      dir_util.mkpath(pymol_path)
      
      # Copy everything we need into it
      dir_util.copy_tree("data",pymol_path+"/data",1,1,0,1,1,0)
      dir_util.copy_tree("test",pymol_path+"/test",1,1,0,1,1,0)
      dir_util.copy_tree("scripts",pymol_path+"/scripts",1,1,0,1,1,0)
      dir_util.copy_tree("examples",pymol_path+"/examples",1,1,0,1,1,0)
      file_util.copy_file("LICENSE",pymol_path+"/",1,1,1,None,1,0)      
      
      # Now build a startup script for PyMOL
      f = open(launch_script,'w')
      if sys.platform!='win32':
         f.write("#!/bin/sh\n")
      pymol_init = re.sub(r"[\/\\][^\/\\]*$","/__pymol_path",sys.modules['pymol'].__file__)
      python_exe = sys.executable
      if python_exe[0:2]=="./":
         python_exe=os.getcwd()+"/"+python_exe[2:]
      if sys.platform!='win32':
         f.write(python_exe+" "+pymol_file+" $*\n")
      else:
         f.write('"%s" "%s" %%1 %%2 %%3 %%4 %%5 %%6 %%7 %%8 %%9\n'%(python_exe,pymol_file))
      f.close()
      os.chmod(launch_script,0755)
      print '''

Created "%s" which can be used to launch PyMOL.  You may wish to copy
this file into a standard location such as /usr/bin or /usr/local/bin.

'''%launch_script
      
except ImportError:
   print '''
   
    The scripts "setup.py" and "setup2.py" are for performing a 
 Distutils-based installation of PyMOL into an external Python
 instance (such as the shared python of a Linux distribution).

 Please run

    python setup.py install
    
 before running this script

    python setup2.py
    
 '''



