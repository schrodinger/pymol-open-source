'''
PYMOL-834
Invalid list evaluation causes stack trace and function doesn't run

Solution: raise CmdException, which is handled quiet by PyMOL
'''

import os
from pymol import cmd, CmdException, testing, stored

class TestPYMOL834(testing.PyMOLTestCase):

    def _load_data(self):
        cmd.fragment('ala', 'm1')
        cmd.map_new('foo_map', 'gaussian', 1.0)

    def _ramp_new(self):
        cmd.ramp_new('foo_ramp', 'foo_map', '[red, white, blue, green]')

    def test834(self):
        '''
        Save a PNG image with width/height specified and
        check if the file exists.
        '''
        self._load_data()
        self.assertRaises(CmdException, self._ramp_new)
