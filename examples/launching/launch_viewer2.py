
# Tell PyMOL we don't want any GUI features and
# that we don't even want the command line.

import pymol

# THIS DOES NOT WORK ON macOS
pymol.finish_launching([ 'pymol', '-qxif', '0' ])
