
# Tell PyMOL we don't want any GUI features.

import __main__
__main__.pymol_argv = [ 'pymol', '-qxi' ]

# Importing the PyMOL module will create the window.

import pymol

# Call the function below before using any PyMOL modules.

pymol.finish_launching()
