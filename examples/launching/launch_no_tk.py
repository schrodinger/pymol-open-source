
# Tell PyMOL we don't want its external Tcl/Tk GUI.

import pymol

# Call the function below before using any PyMOL modules.

pymol.finish_launching([ 'pymol', '-qx' ])
