'''
test various settings
'''

import numpy
from pymol import cmd, testing

class TestSettings(testing.PyMOLTestCase):
    @testing.foreach(True, False)
    def testStickBall(self, use_shader):
        '''
        Test some stick_ball* settings
        '''
        cmd.viewport(100, 100)

        cmd.set('use_shaders', use_shader)

        self.ambientOnly()

        cmd.set('stick_ball')
        cmd.set('stick_ball_ratio', 2.0)
        cmd.set('stick_ball_color', 'blue')
        cmd.set('stick_color', 'red')

        cmd.fragment('gly')
        cmd.orient()

        cmd.color('green')
        cmd.show_as('sticks')

        img = self.get_imagearray()
        self.assertImageHasColor('blue', img)
        self.assertImageHasColor('red', img)
        self.assertImageHasNotColor('green', img)

    @testing.requires_version('1.7.5')
    @testing.requires('incentive')
    @testing.requires('no_edu')
    @testing.requires('gui')
    def testTransparencyMode3(self):
        '''
        Test if something gets rendered in transparency_mode=3
        '''
        cmd.viewport(100, 100)

        cmd.fragment('gly')
        cmd.zoom()

        cmd.set('transparency_mode', 3)
        cmd.set('transparency', 0.5)
        cmd.set('stick_transparency', 0.5)
        cmd.set('sphere_transparency', 0.5)
        cmd.set('cartoon_transparency', 0.5)

        for rep in ['sphere', 'stick', 'surface']:
            cmd.show_as(rep)

            # check on black screen to see if any shader failed
            img = self.get_imagearray()
            self.assertTrue(img[...,:3].any())

    @testing.requires_version('1.7.4')
    @testing.requires('incentive')
    @testing.requires('no_edu')
    @testing.requires('gui')
    def testAA(self):
        '''
        Make a black/white image and check if gray pixels are found with
        antialias_shader=1/2
        '''
        cmd.viewport(100, 100)

        cmd.set('use_shaders', True)

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

    @testing.foreach(0, 1)
    @testing.requires_version('1.7.5')
    @testing.requires('incentive')
    @testing.requires('gui')
    def testTrilines(self, trilines):
        cmd.viewport(100, 100)

        cmd.set('use_shaders', True)

        self.ambientOnly()

        cmd.set('dynamic_width', 0)
        cmd.set('line_width', 5)
        cmd.set('line_smooth', 1)

        cmd.fragment('ethylene')
        cmd.show_as('lines')
        cmd.color('white')
        cmd.orient()

        cmd.set('trilines', trilines)

        # check percentage of covered pixels
        img = self.get_imagearray()
        npixels = img.shape[0] * img.shape[1]
        covered = numpy.count_nonzero(img[...,:3])

        print "covered=", covered, " npixels=" , npixels, " ratio=" , covered / float(npixels)
        self.assertTrue(0.14 < covered / float(npixels) < 0.165)
