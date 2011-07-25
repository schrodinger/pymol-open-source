## Automatically adapted for numpy.oldnumeric Jul 09, 2010 by -c

# PyOpenGL: modified for usage inside of PyMOL

try:
    import numpy.oldnumeric as multiarray
    _numeric = 1
except ImportError:
    _numeric = 0

__version__ = "1.5.6b1"
