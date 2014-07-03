'''
PYMOL-833
png file gets not written instantly if width/height arguments given 
'''

import os
from pymol import cmd, testing, stored
import threading
import unittest

class TestPYMOL833(testing.PyMOLTestCase):

    def _pngExistsAsync(self, filename, sync_in_thread):
        cmd.png(filename, width=100, height=100, ray=0)
        cmd.draw()
        if sync_in_thread:
            cmd.sync()

    @testing.foreach( (True, False), (False, True), (False, False) )
    @unittest.skip("PYMOL-833 is outstanding")
    def testPngExists(self, sync_in_main_thread, sync_in_thread):
        '''
        Save a PNG image with width/height specified and
        check if the file exists. a valid workaround would
        be to call sync(), but if sync() is called inside the
        thread, then it still doesnt work.  ideally we shouldnt
        need to call sync() at all, the png call should just do
        it
        '''
        cmd.pseudoatom('m1')
        cmd.show('spheres')
        with testing.mktemp('.png') as filename:
            th = threading.Thread(target=self._pngExistsAsync, args=[filename, sync_in_thread])
            th.start()
            th.join()
            if sync_in_main_thread:
                cmd.sync()
            self.assertTrue(os.path.exists(filename), 'png file not written')
