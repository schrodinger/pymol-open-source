'''
Test ramp color interpolation along a cylinder
'''

from pymol import cmd, testing, stored

class TestCylInterpol(testing.PyMOLTestCase):

    @testing.requires_version('1.8.0.2')
    def test(self):
        cmd.viewport(200, 100)

        cmd.fragment('ethylene', 'm1')
        cmd.remove('hydro')
        cmd.pseudoatom('p1', 'index 1')
        cmd.disable('p1')
        cmd.ramp_new('r1', 'p1', [0.1, 1.2], ['red', 'blue'])
        cmd.disable('r1')
        cmd.color('r1')
        cmd.show_as('sticks', 'm1')

        cmd.orient('m1')
        cmd.move('z', 11)

        self.ambientOnly()

        img = self.get_imagearray()

        self.assertImageHasColor('red', img)
        self.assertImageHasColor('blue', img)
        self.assertImageHasColor([.5, .0, .5], img, delta=1)
