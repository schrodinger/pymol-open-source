
# Tell PyMOL we don't want any GUI features and
# that we don't even want the command line.

import __main__
__main__.pymol_argv = [ 'pymol', '-qxif', '0' ]

# Importing the PyMOL module will create the window.

import pymol

# Call the function below before using any PyMOL modules.

pymol.finish_launching()
