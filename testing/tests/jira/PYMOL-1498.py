'''
MAE sticks representation
'''

from pymol import cmd, testing, stored

@testing.requires('incentive')
class Test1498(testing.PyMOLTestCase):

    def test(self):
        cmd.load(self.datafile("repr.mae"))

        self.assertEqual(60, cmd.count_atoms('repr.x_lines and rep lines'))
        self.assertEqual(60, cmd.count_atoms('repr.x_sticks and rep sticks'))
        
        self.assertEqual(60, cmd.count_atoms('repr.x_cartoon and rep cartoon'))
        self.assertEqual(60, cmd.count_atoms('repr.x_cart_lines and rep cartoon'))
        self.assertEqual(60, cmd.count_atoms('repr.x_cart_lines and rep lines'))

        self.assertEqual(60, cmd.count_atoms('repr.x_spheres and rep spheres'))
        self.assertEqual(60, cmd.count_atoms('repr.x_thinsticks and rep sticks'))
        
        self.assertEqual(60, cmd.count_atoms('repr.x_ballstick and rep sticks'))
        self.assertEqual(60, cmd.count_atoms('repr.x_ballstick and rep spheres'))
