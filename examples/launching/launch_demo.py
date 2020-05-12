# Importing the PyMOL module will create the window.

import pymol

# Call the function below before using any PyMOL modules.

# THIS DOES NOT WORK ON macOS
pymol.finish_launching()

# Now we can import cmd

from pymol import cmd

cmd.load("$PYMOL_PATH/test/dat/pept.pdb")
cmd.show("sticks")

