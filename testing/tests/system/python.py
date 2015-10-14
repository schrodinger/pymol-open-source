'''
Python installation tests
'''

import sys
from pymol import testing

required_modules = [
    'Image',
    'numpy',
]

# when building on OSX with system python, some modules are not
# available by default. Only test for them if we build with our
# own python distribution which should ship all these modules.
if (not sys.executable.startswith('/System/Library/Frameworks/Python.framework/Versions/2.7')
        and not testing.PYMOL_EDU
        ):
    required_modules += [
#        'OpenGL', # we don't have it on Windows for PyMOL > 1.7.6
        'matplotlib',
    ]

class TestSystem(testing.PyMOLTestCase):

    def testHasModules(self):
        failed = []
        for name in required_modules:
            try:
                __import__(name, level=0)
            except ImportError:
                failed.append(name)
        self.assertEqual(failed, [])
