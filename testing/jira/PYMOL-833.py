'''
PYMOL-833
png file gets not written instantly if width/height arguments given 
'''

import os
from pymol import cmd, testing, stored

class TestPYMOL833(testing.PyMOLTestCase):

    @testing.requires('gui')
    def testPngExists(self):
        '''
        Save a PNG image with width/height specified and
        check if the file exists.
        '''
        cmd.pseudoatom('m1')
        cmd.show('spheres')

        with testing.mktemp('.png') as filename:
            cmd.png(filename, width=100, height=100, ray=0)
            cmd.draw()
            self.assertTrue(os.path.exists(filename), 'png file not written')
