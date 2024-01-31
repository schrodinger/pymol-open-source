'''
bg_image_filename in session files
'''

import shutil
import os, unittest
from pymol import cmd, testing

@testing.requires('incentive', 'gui')
class Test1571(testing.PyMOLTestCase):

    def test(self):
        cmd.set('opaque_background')
        with testing.mktemp('tmp.pse') as session_filename:
            with testing.mktemp('bg.png') as bg_image_filename:
                cmd.bg_color('red')
                cmd.png(bg_image_filename)
                cmd.bg_color('blue')
                cmd.set('bg_image_filename', bg_image_filename)
                cmd.save(session_filename)

            cmd.load(session_filename)
        x = cmd.get('bg_image_filename')
        self.assertTrue(x.startswith('data:'))
        self.assertImageHasColor('red', delta=4,
                msg='bg_image_filename=data: not rendered')

