# PyOpenGL: modified for usage inside of PyMOL

try:
    import multiarray
    _numeric = 1
except ImportError:
    _numeric = 0

__version__ = "1.5.6b1"
