#
# This script only applies if you are performing a Python Distutils-based
# installation of PyMOL.
#
# Please run this script after using Distutils "setup.py":
#   python setup.py build
#   python setup.py install
#   python setup2.py

import os
import re
import sys
from distutils import dir_util,file_util

if sys.platform=='win32':
   launch_script = "pymol.bat"
elif sys.platform=='cygwin':
   launch_script = "pymol"
else:
   launch_script = "pymol"

help_str = '''
The scripts "setup.py" and "setup2.py" are for performing a
Distutils-based installation of PyMOL into an external Python instance,
such as the shared python of a Linux distribution.
'''

if ((('uninstall' not in sys.argv) and ('install' not in sys.argv)) or
    ((('uninstall' in sys.argv) and ('install' in sys.argv)))):
   print help_str
   print "usage:\n    python setup2.py install \n    python setup2.py uninstall\n"
   sys.exit(0)
   
if 'uninstall' in sys.argv:
   uninstall = 1
else:
   uninstall = 0
   
try:
   pymol_launch = 3
   import pymol
   
   if not (os.path.exists("data") and os.path.exists("LICENSE")):
      print ' Please run "python setup2.py" from the PyMOL source directory.'
   else:
      pymol_file = sys.modules['pymol'].__file__
      if pymol_file[-4:]==".pyc":
         pymol_file = pymol_file[0:-1]

      site_packages = os.path.split(os.path.split(pymol_file)[0])[0]

      if uninstall:
         print "\nUninstalling from %s ..."%site_packages
         pymol_path = os.path.join(site_packages,'pymol')
         if os.path.exists(pymol_path):
            print " removing %s..."%pymol_path
            dir_util.remove_tree(pymol_path)
         pmg_tk_path = os.path.join(site_packages,'pmg_tk')
         if os.path.exists(pmg_tk_path):
            print " removing %s..."%pmg_tk_path
            dir_util.remove_tree(pmg_tk_path)
         pmg_wx_path = os.path.join(site_packages,'pmg_wx')
         if os.path.exists(pmg_wx_path):
            print " removing %s..."%pmg_wx_path
            dir_util.remove_tree(pmg_wx_path)
         chempy_path = os.path.join(site_packages,'chempy')
         if os.path.exists(chempy_path):
            print " removing %s..."%chempy_path
            dir_util.remove_tree(chempy_path)

         if 'pmw' in sys.argv:
            pmw_path = os.path.join(site_packages,'Pmw')
            if os.path.exists(pmw_path):
               print " removing %s..."%pmw_path
               dir_util.remove_tree(pmw_path)

         if os.path.exists(launch_script):
            print " removing %s..."%launch_script
            os.unlink(launch_script)
               
      else:
         print "\n Installing into %s ..."%site_packages

#         pymol_path = re.sub(r"[\/\\][^\/\\]*$","/pymol_path",pymol_file)
         pymol_path = os.path.join(os.path.split(pymol_file)[0],"pymol_path")

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

         if sys.platform=='darwin':
            f.write('if [ "$DISPLAY" == "" ]; then\nDISPLAY=":0.0"\nexport DISPLAY\nfi\n')

         pymol_init = re.sub(r"[\/\\][^\/\\]*$","/__pymol_path",sys.modules['pymol'].__file__)
         python_exe = sys.executable
         if python_exe[0:2]=="./":
            python_exe=os.getcwd()+"/"+python_exe[2:]
         if sys.platform!='win32':
            f.write(python_exe+" "+pymol_file+" \"$@\"\n")
         else:
            f.write('"%s" "%s" %%1 %%2 %%3 %%4 %%5 %%6 %%7 %%8 %%9\n'%(python_exe,pymol_file))
         f.close()
         os.chmod(launch_script,0755)

         if 'pmw' in sys.argv: # only working under unix shell
            os.system("tar -C %s -zxvf modules/pmg_tk/pmw.tgz"%site_packages)
            
         print '''
 Created "./%s" which can be used to launch PyMOL.  You may wish to copy
 this file into a standard location such as /usr/bin or /usr/local/bin.
   '''%launch_script
      
except ImportError:
   print help_str
   print ''' Please run:

     python setup.py install
    
 before this script:

     python setup2.py install
     python setup2.py uninstall    
 '''
   print "Error: Unable to import pymol.  Have you run setup.py?"
   


