'''
PYMOL-2668
crash with cylinder shader and dist dash < 1.0 Angstrom
'''

from pymol import cmd, CmdException, testing, stored

class TestPYMOL2668(testing.PyMOLTestCase):

    @testing.requires_version('1.8.0.0') # fixed in 1.8.0.1
    @testing.requires('gui')
    def test2668(self):
        cmd.pseudoatom('m1', pos=(-.4, 0, 0))
        cmd.pseudoatom('m2', pos=( .4, 0, 0))
        cmd.distance('d0', 'm1', 'm2')
        cmd.draw()
