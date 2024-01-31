'''
background during stereo doesn't change
'''

import os
from pymol import cmd, CmdException, testing, stored

@testing.requires('gui')
class Test(testing.PyMOLTestCase):

    def test(self):
        cmd.set('bg_rgb', 'red')
        cmd.stereo('anaglyph')
        self.assertImageHasColor('red')
