'''
PYMOL-802
Ray tracing is broken in 1.6.0

Solution: for edge sampling, if opaque background and result pixel
          is not opaque, fill it in with the background color by interpolating
          rgb values and setting alpha to 1

Test: as many combinations of settings for ray tracing to make sure alpha is correctly
      set in all pixels
'''

import os
from pymol import cmd, CmdException, testing, stored

class TestPYMOL802(testing.PyMOLTestCase):

    @testing.foreach((-1,'on',True), (1,'on',True), (0,'on',False), (-1,'off',False))
    def test802(self, ray_opaque_background, opaque_background, trans_check):
        cmd.load(self.datafile('1oky.pdb.gz'))
        cmd.hide()
        cmd.show('cartoon')
        cmd.set('ray_opaque_background', ray_opaque_background)
        cmd.set('opaque_background', opaque_background)
        if trans_check:
            self.assertImageHasNoTransparency(self.get_imagearray(width=100, height=100, ray=1))
        else:
            self.assertImageHasTransparency(self.get_imagearray(width=100, height=100, ray=1))
            
