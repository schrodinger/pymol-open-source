## Automatically adapted for numpy.oldnumeric Jul 09, 2010 by -c

import sys
from pymol import opengl

"""
This module performs some nasty trickery in order to automatically
detect OpenGL errors, by wrapping each callable function with a
function which calls the function, but afterwards checks to see if an
error code was set, in which case it raises an exception

This is highly experimental.  To disable it, call the fast()
function.  To reenable it, call the careful() function.

"""

if opengl._numeric:
    from numpy.oldnumeric import ArrayType
    try:
        import _opengl_num
        _opengl = _opengl_num
    except ImportError:
        import _opengl
    except SystemError:
        import _opengl
    try:
        import openglutil_num
        openglutil =openglutil_num
    except ImportError:
        import openglutil
    except SystemError:
        import openglutil      
else:
    import _opengl
    import openglutil

origdict = _opengl.__dict__.copy()
origdict.update(openglutil.__dict__)

from glconst import *

_safetylevel = 0

# Constants
import glconst
constant_names = {}
for k, v in glconst.__dict__.items():
    if type(v) == type(1):
        constant_names[v] = k
    
error = origdict['glGetError']
class Error:
    def __init__(self, num):
        self.num = num
        self.str = constant_names.get(num, 'Unknown Error')
    def __str__(self):
        return '%s (%d)' % (self.str, self.num)

class CarefulFunction:
    def __init__(self, name, function):
        self.name = name
        self.func = function
    def __repr__(self):
        return "<careful function wrapper around %s>" % self.function
    def __call__(self, *args, **kw):
        # did someone check in a debug (non-production) version of this module?
        # this was _not_ done in the 1.1 release (causes huge gobs of printing
        # when the demos are run).
##        print "calling %s with %s, %s" % (self.func, args, kw)
        retval = apply(self.func, args, kw)
        err = error()
        if err:
            if _verysafe:
                raise Error(err)
            else:
                print Error(err)
        return retval

carefuldict = {}
for name, func in origdict.items():
    if callable(func):
        carefuldict[name] = CarefulFunction(name, func)
# These do the same sorts of things that the C versions would

def careful():
    cd = carefuldict.copy()
    if cd.has_key('error'): del cd['error']
    sys.modules['pymol.opengl.gl'].__dict__.update(cd)

import string
def fast():
    cd = origdict.copy()
    if cd.has_key('error'): del cd['error']
    if cd.has_key('glconst'): del cd['glconst']
    sys.modules['pymol.opengl.gl'].__dict__.update(cd)

# by default, be careful
if _safetylevel >= 1:
    careful()
else:
    fast()
if _safetylevel >= 2:
    _verysafe = 1
else:
    _verysafe = 0
    
def glGetIntegerv(*args):
    r = apply(glGetDoublev, args)
    if isinstance(r, type(())):
        return tuple(map(int, r))
    elif opengl._numeric and isinstance(r, ArrayType): 
        return tuple(r.astype('i').tolist())
    else:
        return int(r)
glGetInteger = glGetIntegerv
glGetBooleanv = glGetBoolean = glGetFloatv = glGetFloat = glGetDouble

# for b/w compatibility w/ 1.0??
GL_LEVELS = 0x0010

__version__ = '1.5.4'
