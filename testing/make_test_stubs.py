import sys, os
from pymol import cmd

def make_test_stubs(module, out=sys.stdout):
    '''
    Utility function to generate a test class for a PyMOL module,
    filled with skipped stub methods.

    Example:

    PyMOL> run make_test_stubs.py
    PyMOL> make_test_stubs pymol.fitting, /tmp/fitting.py
    '''
    if isinstance(module, str):
        __import__(module)
        module = sys.modules[module]

    if isinstance(out, str):
        out = open(out, 'w')

    print >> out, '''
from pymol import cmd, testing, stored

class Test%s(testing.PyMOLTestCase):
''' % (module.__name__.split('.')[-1].capitalize())

    type_func = type(cmd.align)
    ns = vars(module)
    for name in sorted(ns):
        value = ns[name]
        if isinstance(value, type_func) and \
                value.__module__ == module.__name__ and \
                value == getattr(cmd, name, None):
            nameCap = name.replace('_', ' ').title().replace(' ', '')
            print >> out, '    def test%s(self):' % nameCap
            print >> out, '        cmd.%s' % value.__name__
            print >> out, '        self.skipTest("TODO")'
            print >> out, ''

cmd.extend('make_test_stubs', make_test_stubs)
