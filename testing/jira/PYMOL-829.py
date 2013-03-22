'''
Regression test for PYMOL-829
Certain mouse actions in edit mode cause objects to be ungrouped
'''

import unittest
from pymol import cmd, testing

@testing.requires('gui')
class Test829(testing.PyMOLTestCase):

    def testAngleColor(self):
        cmd.fragment('ala')
        cmd.color('blue')
        cmd.angle('ang01', 'index 8', 'index 5', 'index 2')
        cmd.hide('label')
        cmd.set('fog', '0')
        cmd.set('ambient', '1')
        cmd.set('angle_color', 'red')
        cmd.zoom()
        #need to check to make sure screen image has red in it
        img_array = self.get_imagearray(width=100, height=100, ray=0)
        self.assertImageHasColor('red', img_array, delta=[10,0,0])

    def testDihedralColor(self):
        cmd.fragment('ala')
        cmd.color('blue')
        cmd.dihedral('dih01', 'index 8', 'index 5', 'index 2', 'index 1')
        cmd.hide('label')
        cmd.set('fog', '0')
        cmd.set('ambient', '1')
        cmd.set('dihedral_color', 'red')
        cmd.zoom()

        #need to check to make sure screen image has red in it
        img_array = self.get_imagearray(width=100, height=100, ray=0)
        self.assertImageHasColor('red', img_array)

