
# Before using this script, you must built a monolithic
# PyMOL DLL which gets deposited into "pymol/_cmd.dll"

# (monolithic means that it also includes the C code for
# Numeric, PyOpenGL, ExtensionClass, and sglite)

build = "c:/pymolws"
product = build+"/win32module"

import os
import shutil

#
# copy all .py and .dll files from modules
#

def rec_copy(source,target,path):
    list = os.listdir(source+"/"+path)
    for entry in list:
        src = source+"/"+path+"/"+entry
        if entry not in ['pmg_wx','CVS','launch_pymol.py','compile_pymol.py']:
            if os.path.isfile(src):
                if entry[-4:] not in ['.ilk','.pyc']:
                    trg = target+"/"+path
                    if not os.path.exists(trg):
                        os.mkdir(trg)
                    trg = trg + "/" +entry    
                    shutil.copyfile(src,trg)
            elif os.path.isdir(src):
                trg = target+"/"+path
                if not os.path.exists(trg):
                    os.mkdir(trg)
                rec_copy(source,target,path+"/"+entry)

if not os.path.exists(product):
    os.mkdir(product)
rec_copy(build+"/modules",product,"pymol")
shutil.copyfile(build+"/glut32.dll",product+"/pymol/glut32.dll")
pymol_path = product+"/pymol/pymol_path"
if not os.path.exists(pymol_path):
    os.mkdir(pymol_path)
rec_copy(build,pymol_path,"data")
rec_copy(build,pymol_path,"test")
modules = pymol_path + "/modules"
if not os.path.exists(modules):
    os.mkdir(modules)    
rec_copy(build+"/modules",modules,"Pmw")
rec_copy(build+"/modules",modules,"chempy")
rec_copy(build+"/modules",modules,"pmg_tk")
rec_copy(build+"/modules",modules,"Numeric")

examples = build + "/examples/launching"
shutil.copyfile(examples+"/launch.py",product+"/launch.py")
shutil.copyfile(examples+"/launch_demo.py",product+"/launch_demo.py")
shutil.copyfile(examples+"/launch_no_tk.py",product+"/launch_no_tk.py")
shutil.copyfile(examples+"/launch_viewer1.py",product+"/launch_viewer1.py")
shutil.copyfile(examples+"/launch_viewer2.py",product+"/launch_viewer2.py")

f=open(product+"/README.TXT",'w')
f.write('''
This is an experimental win32 module build of PyMOL.

You can either run PyMOL directly from this folder, with "pymol" being
implicitly added to the module search path, or you can or install
PyMOL globally by moving the "pymol" directory into the main Python folder.

Please direct any questions to warren@delanoscientific.com

Several launch scripts are provided to get you started.
''')
f.close()
print "Win32 Module Product Built"
print "Hit return to close"
import sys
print sys.stdin.readline()

                
