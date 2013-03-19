'''
PYMOL-1240
PyMOL should save ANISOU records in PDBs -- and other records
'''

from pymol import cmd, testing, stored

v_pdbstr_anisou = (
    'ATOM      1  N   ARG A 197       5.287   1.830   7.079  1.00 30.85           N  \n'
    'ANISOU    1  N   ARG A 197     3652   4226   3841   -164    -54    -45       N  \n'
    'END\n')

class TestPYMOL1241(testing.PyMOLTestCase):

    def testSaveANISO(self):
        cmd.read_pdbstr(v_pdbstr_anisou, 'm1')
        v = cmd.get_pdbstr()
        self.assertEqual(v, v_pdbstr_anisou, 'ANISOU records missing:\n' + v)
