from pymol import opengl

from _glut import *
try:
    from _glut import *
    
except ImportError:
    print "Importing the GLUT library failed."
    import sys
    if sys.platform == 'win32':
        print "The GLUT32.DLL file should be in your PATH"
    print "Consult the documentation for opengl.glut"
    print "for installation suggestions."
	
from glutconst import *
