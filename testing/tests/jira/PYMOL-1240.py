'''
PYMOL-1240
PyMOL should save ANISOU records in PDBs -- and other records
'''

from pymol import cmd, testing, stored

v_pdbstr_anisou = (
    'ATOM      1  N   ARG A 197       5.287   1.830   7.079  1.00 30.85           N  \n'
    'ANISOU    1  N   ARG A 197      652   4226   3841   -164    -54    -45       N  \n'
    'END\n')

v_pdbstr_rotated = (
    'ATOM      1  N   ARG A 197      24.081   8.670   7.079  1.00 30.85           N  \n'
    'ANISOU    1  N   ARG A 197     1839   3790   3091  -1033   1258    503       N  \n'
    'END\n')

class TestPYMOL1240(testing.PyMOLTestCase):

    def testSaveANISO(self):
        cmd.set('pdb_use_ter_records', 0);

        cmd.read_pdbstr(v_pdbstr_anisou, 'm1')
        v = cmd.get_pdbstr()
        self.assertEqual(v, v_pdbstr_anisou, 'ANISOU records missing:\n' + v)

        cmd.rotate('y', 30, object='m1')
        cmd.translate([20, 0, 0], object='m1')
        cmd.rotate('z', 20, object='m1')

        v = cmd.get_pdbstr()
        self.assertEqual(v, v_pdbstr_rotated, 'ANISOU not rotated in PDB string' + v)

        v = cmd.get_model('m1').atom[0].u_aniso
        self.assertArrayEqual(v, [0.183853, 0.378995, 0.309052, -0.103350, 0.125751, 0.050349],
                1e-5, 'ANISOU not rotated in chempy model')
