'''
vis file support
'''

import unittest
from pymol import cmd, testing, stored

try:
    import h5py
except ImportError:
    h5py = None

@unittest.skipIf(not h5py, 'TODO: h5py required')
@testing.requires('incentive')
class Test1441(testing.PyMOLTestCase):

    def test_map(self):
        cmd.load(self.datafile("Structure_potential.vis"))

        self.assertEqual(1, len(cmd.get_names_of_type('object:map')))
        self.assertEqual(2, len(cmd.get_names_of_type('object:surface')))

    def test_surface(self):
        cmd.load(self.datafile("surf2.vis"))

        self.assertEqual(1, len(cmd.get_names_of_type('object:cgo')))
