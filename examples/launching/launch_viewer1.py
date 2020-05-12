
# Tell PyMOL we don't want any GUI features.

import pymol

# THIS DOES NOT WORK ON macOS
pymol.finish_launching([ 'pymol', '-qxi' ])
