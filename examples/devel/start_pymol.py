# ================================================================
# NOTE: THIS APPROACH IS NOT RECOMMENDED SINCE THE PyMOL LAUNCH
#       PROCESS IS SUBJECT TO CHANGE WITHOUT WARNING
#
# The recommended approach is to just use PyMOL as your python
# interpreter:
#
#     "pymol -r <script.py>" or
#     "pymol -qcr <script.py>" instead of
#     "python <script.py>"
# 
# However, If must invoke PyMOL from an external interpreter, then
# you willl need to:
#
#  (1) prepare the environment first by sourcing "pymol.csh"
#      or an equivalent shell script which defines PYMOL_PATH,
#      LD_LIBRARY_PATH, TCL_LIBRARY, etc.
#
#  (2) include the following code block at the *top* level of your
#      program script.
#
# ================================================================
# PyMOL launch code
#
# === Provide arguments to PyMOL (first one must be "pymol")
pymol_argv = [ "pymol", "-q" ]
#
# === Launch the PyMOL thread(s)
import os,threading,__main__,__builtin__
threading.Thread(target=__builtin__.execfile,
     args=(os.environ['PYMOL_PATH']+"/modules/launch_pymol.py",
           __main__.__dict__,__main__.__dict__)).start()
#
# === Wait until PyMOL is ready to receive commands
e=threading.Event()
while not hasattr(__main__,'pymol'):
   e.wait(0.01)
while not pymol._cmd.ready():
   e.wait(0.01)
#
# PyMOL is now launched, you can now import "pymol" modules.
# ===============================================================

from pymol import cmd

cmd.load("$PYMOL_PATH/test/dat/pept.pdb")
print " The surface area is: %8.3f"%cmd.get_area()


