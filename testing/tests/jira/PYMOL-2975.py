'''
some elements assigned wrong from SDF and XYZ files
'''

import os
from pymol import cmd, CmdException, testing, stored

@testing.requires_version('1.9')
class TestElem(testing.PyMOLTestCase):

    @testing.foreach('mol', 'xyz')
    def test(self, ext):
        cmd.load(self.datafile('elements.' + ext))
        for elem in 'Ne Si Cd Sn Ce Nd Sm Dy Ho Hf Hg Np Cm Cf No Db Sg Hs Ds Cn'.split():
            self.assertEqual(cmd.count_atoms('elem ' + elem), 1, elem)
