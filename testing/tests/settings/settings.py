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
    @testing.requires('gui')
    def testTrilines(self, trilines):
        cmd.viewport(100, 100)

        cmd.set('use_shaders', True)

        self.ambientOnly()

        cmd.set('dynamic_width', 0)
        cmd.set('line_width', 5)
        cmd.set('line_smooth', 0)

        cmd.fragment('ethylene')
        cmd.show_as('lines')
        cmd.color('white')
        cmd.orient()

        cmd.set('trilines', trilines)

        # check percentage of covered pixels
        img = self.get_imagearray()
        npixels = img.shape[0] * img.shape[1]
        covered = numpy.count_nonzero(img[...,:3])
        ratio = covered / float(npixels)
        msg = "covered=%d npixels=%d ratio=%f" % (covered, npixels, ratio)
        self.assertTrue(0.14 < ratio < 0.165, msg)

    @testing.requires_version('1.7.5')
    @testing.requires('gui')
    def testChromadepth(self):
        cmd.viewport(100, 100)

        self.ambientOnly()

        cmd.set('use_shaders')
        cmd.set('chromadepth')
        cmd.set('orthoscopic')

        cmd.pseudoatom(pos=(0, 0, -5))
        cmd.pseudoatom(pos=(0, 0, 5))
        cmd.color('gray')
        cmd.show_as('spheres')
        cmd.zoom()
        cmd.turn('y', 20)

        img = self.get_imagearray()
        self.assertImageHasColor('blue', img)
        self.assertImageHasColor('red', img)

    @testing.requires_version('2.3')
    def testAtomCartoonTransparency(self):
        '''
        Test atom-level cartoon_transparency
        '''
        self.ambientOnly()
        cmd.viewport(100, 100)
        cmd.set('opaque_background')
        cmd.set('ray_shadow', 0)

        cmd.load(self.datafile('1oky-frag.pdb'))
        cmd.orient()

        cmd.show_as('cartoon')
        cmd.select('sele', 'resi 106-111', 0)

        # hide blended parts between residues
        cmd.cartoon('skip', 'resi 105+112')

        cmd.color('0xFF0000')
        cmd.color('0x0000FF', 'sele')

        # all opaque
        img = self.get_imagearray()
        self.assertImageHasColor('0xFF0000', img)
        self.assertImageHasColor('0x0000FF', img)

        # object-level transparency
        cmd.set('cartoon_transparency', 0.4)
        img = self.get_imagearray()
        self.assertImageHasColor('0x990000', img) # 40%
        self.assertImageHasColor('0x000099', img) # 40%

        # atom-level full-opaque
        cmd.set('cartoon_transparency', 0.0, 'sele')
        img = self.get_imagearray()
        self.assertImageHasColor('0x990000', img) # 40%
        self.assertImageHasColor('0x0000FF', img) #  0%

        # atom-level semi-transparent
        cmd.set('cartoon_transparency', 0.6, 'sele')
        img = self.get_imagearray()
        self.assertImageHasColor('0x990000', img) # 40%
        self.assertImageHasColor('0x000066', img) # 60%

        # atom-level full-transparent (expect only two color values)
        cmd.set('cartoon_transparency', 0.0)
        cmd.set('cartoon_transparency', 1.0, 'sele')
        img = self.get_imagearray()
        self.assertImageHasColor('0xFF0000', img)
        self.assertTrue(len(numpy.unique(img[...,:3])) == 2)
