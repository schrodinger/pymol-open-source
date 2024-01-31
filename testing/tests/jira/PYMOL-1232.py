'''
PYMOL-1232
Cannot load session file after "wizard message"
'''

from pymol import cmd, CmdException, testing, stored

class TestPYMOL1232(testing.PyMOLTestCase):

    def test1232(self):
        with testing.mktemp('.pse') as filename:
            cmd.save(filename)
            cmd.wizard('message', 'foo')
            cmd.load(filename)
