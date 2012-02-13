# $Id$
import sys
from pymol import opengl

if opengl._numeric:
    from numpy.oldnumeric import ArrayType
    try:
        import _glu_num
        _glu = _glu_num
    except ImportError:
        import _glu
    except SystemError:
        import _glu
else:
    import _glu

from gluconst import *
from pymol.opengl.gl import Error, CarefulFunction

origdict = _glu.__dict__.copy()
sys.modules['pymol.opengl.glu'].__dict__.update(origdict)

carefuldict = {}
for name, func in origdict.items():
    if callable(func):
        carefuldict[name] = CarefulFunction(name, func)
# These do the same sorts of things that the C versions would

def careful():
    cd = carefuldict.copy()
    if cd.has_key('error'): del cd['error']
    sys.modules['pymol.opengl.glu'].__dict__.update(cd)

import string
def fast():
    cd = origdict.copy()
    if cd.has_key('error'): del cd['error']
    if cd.has_key('glconst'): del cd['glconst']
    sys.modules['pymol.opengl.glu'].__dict__.update(cd)
