'''
PYMOL-2095
background in ray tracing not drawn properly

'''

import os
from pymol import cmd, CmdException, testing, stored

class TestPYMOL834(testing.PyMOLTestCase):

    @testing.foreach((1,0,False), (1,1,True), (0,1,False))
    
    def test2095(self, bg_gradient, opaque_background, hasColor):
        cmd.load(self.datafile("1oky-frag.pdb"))
        cmd.zoom(buffer=20)
        cmd.set("bg_gradient", bg_gradient)
        cmd.set("bg_rgb_top", "blue")
        cmd.set("bg_rgb_bottom", "blue")
        cmd.set("opaque_background", opaque_background)
        cmd.color("red")
        if hasColor:
            self.assertImageHasColor('blue', self.get_imagearray(width=100, height=100, ray=1))
        else:
            self.assertImageHasNotColor('blue', self.get_imagearray(width=100, height=100, ray=1))

