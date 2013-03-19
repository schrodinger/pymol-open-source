'''
Python installation tests
'''

from pymol import testing

required_modules = [
    'Image',
    'numpy',
    'matplotlib',
    'OpenGL',
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
