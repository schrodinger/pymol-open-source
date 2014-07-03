'''
atom numbering following h_add
'''

from pymol import cmd, testing, stored

class Test1553(testing.PyMOLTestCase):

    def test(self):
        cmd.fab('AG')
        cmd.remove('hydro')
        cmd.h_add('resn ALA')

        # with this bug: count = 6
        self.assertEqual(11, cmd.count_atoms('byres (last hydro)'))
