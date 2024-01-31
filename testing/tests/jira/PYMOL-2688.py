'''
PYMOL-2688 1.8.0 crashes with 1.5.0 session file
'''

from pymol import cmd, CmdException, testing, stored

class TestPYMOL2688(testing.PyMOLTestCase):

    @testing.requires_version('1.8.0.4')
    def test(self):
        cmd.load('PYMOL-2688.pse')
        self.assertEqual(['myramp'], cmd.get_names())
