
# Tell PyMOL we don't want its external Tcl/Tk GUI.

import __main__
__main__.pymol_argv = [ 'pymol', '-qx' ]

# Importing the PyMOL module will create the window.

import pymol

# Call the function below before using any PyMOL modules.

pymol.finish_launching()
