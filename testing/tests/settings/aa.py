'''
test antialias shaders
'''

import numpy
from pymol import cmd, testing

class TestAA(testing.PyMOLTestCase):

    @testing.requires_version('1.7.4')
    @testing.requires('no_edu')
    @testing.requires('incentive')
    @testing.requires('gui')
    def testAA(self):
        cmd.viewport(100, 100)

        self.ambientOnly()

        cmd.fragment('gly')
        cmd.show_as('spheres')
        cmd.color('white')
        cmd.zoom()

        # b/w image, we expect only two color values
        img = self.get_imagearray()
        self.assertTrue(len(numpy.unique(img[...,:3])) == 2)

        for aa in (1, 2):
            cmd.set('antialias_shader', aa)

            # smoothed edges, we expect more than two color values
            img = self.get_imagearray()
            self.assertTrue(len(numpy.unique(img[...,:3])) > 2)
