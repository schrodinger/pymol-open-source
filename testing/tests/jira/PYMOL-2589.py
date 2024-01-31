'''
some stick_ball caps missing
'''

import pymol
from pymol import cmd, util, testing, stored

class Test2589(testing.PyMOLTestCase):

    def test(self):
        if not (pymol.invocation.options.no_gui or cmd.get_setting_int('use_shaders')):
            self.skipTest('no ray or shaders')

        self.ambientOnly()
        cmd.set('valence', 0)

        cmd.viewport(350, 200)
        cmd.fragment('ile')
        cmd.remove('hydro')
        cmd.show_as('sticks')
        cmd.orient()

        cmd.set('stick_ball', 0)
        img1 = self.get_imagearray()

        cmd.set('stick_ball', 1)
        img2 = self.get_imagearray()

        self.assertImageEqual(img1, img2, count=10)
