# ================================================================
# NOTE: THIS APPROACH IS NOT RECOMMENDED, BUT THE PyMOL LAUNCH
#       PROCESS IS SUBJECT TO CHANGE
#
# The recommended approach is to just use:
#     "pymol -r <script.py>" or
#     "pymol -q -c -r <script.py>" instead of
#     "python <script.py>"
# 
# However, If must invoke PyMOL from an external interpreter,
# you need to:
#
#  (1) prepare the environment by sourcing "pymol.csh"
#      or an equivalent (to define PYMOL_PATH,
#      LD_LIBRARY_PATH, TCL_LIBRARY, etc.)
#
#  (2) include the following code at the *top* level of your
#      program script.
#
# === Define any argument to PyMOL you might want to indicate
pymol_argv = [ "pymol", "-c", "-q" ]
# (Note: first argument must be "pymol")
#
# === Launch the PyMOL thread
import os,threading,__main__,__builtin__
t = threading.Thread(target=__builtin__.execfile,
    args=(os.environ['PYMOL_PATH']+"/modules/launch_pymol.py",
          __main__.__dict__,__main__.__dict__))
t.start()
#
# === Wait until PyMOL is ready to receive commands
e=threading.Event()
while not hasattr(__main__,'pymol'):
   e.wait(0.01)
while not pymol._cmd.ready():
   e.wait(0.01)
# You can now safely make use of pymol
# ===============================================================

#from pymol import cmd

#cmd.load("$PYMOL_PATH/test/dat/pept.pdb")
#print " The surface area is: %8.3f"%cmd.get_area()


