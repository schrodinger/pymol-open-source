'''
vis file support
'''

from pymol import cmd, testing, stored

@testing.requires('incentive')
class Test1441(testing.PyMOLTestCase):

    def test_map(self):
        cmd.load(self.datafile("Structure_potential.vis"))

        self.assertEqual(1, len(cmd.get_names_of_type('object:map')))
        self.assertEqual(2, len(cmd.get_names_of_type('object:surface')))

    def test_surface(self):
        cmd.load(self.datafile("surf2.vis"))

        self.assertEqual(1, len(cmd.get_names_of_type('object:cgo')))
