'''
test transparency settings
'''

from pymol import cmd, testing

class TestTransparency(testing.PyMOLTestCase):

    @testing.requires('incentive')
    @testing.requires('no_edu')
    @testing.requires('gui')
    def testTransparencyMode3(self):
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
