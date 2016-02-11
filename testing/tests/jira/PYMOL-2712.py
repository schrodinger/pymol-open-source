'''
reinitialize has deferred effects
'''

from pymol import cmd, CmdException, testing, stored

class TestPYMOL2712(testing.PyMOLTestCase):

    @testing.requires_version('1.8.1.0')
    def test_reinit(self):
        # bug was: reinitialization happens after 'frag gly'
        cmd.do('reinitialize;fragment gly')
        self.assertEqual(cmd.count_atoms(), 7)
