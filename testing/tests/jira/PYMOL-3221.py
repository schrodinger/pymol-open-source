'''
PYMOL-3221
MAE ANISOU support
'''

from pymol import cmd, testing

v_pdbstr_anisou = (
    'ATOM      1  N   GLU A 114      24.832  -7.270  -5.728  1.00 33.91           N  \n'
    'ANISOU    1  N   GLU A 114     6968   3709   2207   -518    495    146       N  \n'
    'ATOM      2  CA  GLU A 114      25.839  -6.416  -5.102  1.00 33.68           C  \n'
    'ATOM      3  CB  GLU A 114      25.893  -5.047  -5.787  1.00 34.17           C  \n'
    'ATOM      4  CG  GLU A 114      26.342  -5.083  -7.241  1.00 38.31           C  \n'
    'ANISOU    4  CG  GLU A 114     7707   4148   2699   -531    590    255       C  \n'
    'ATOM      5  CD  GLU A 114      25.207  -5.382  -8.202  1.00 40.38           C  \n'
    'ANISOU    5  CD  GLU A 114     7970   4550   2823   -523    543    278       C  \n'
    'TER   \n'
    'END\n')

v_pdbstr_rotated = (
    'ATOM      1  N   GLU A 114      26.964  -6.052  -5.274  1.00 33.91           N  \n'
    'ANISOU    1  N   GLU A 114     3700   2782   6402   -635  -1525   -432       N  \n'
    'ATOM      2  CA  GLU A 114      26.822  -5.036  -6.315  1.00 33.68           C  \n'
    'ATOM      3  CB  GLU A 114      25.355  -4.880  -6.726  1.00 34.17           C  \n'
    'ATOM      4  CG  GLU A 114      24.745  -6.113  -7.378  1.00 38.31           C  \n'
    'ANISOU    4  CG  GLU A 114     4106   3362   7086   -636  -1669   -513       C  \n'
    'ATOM      5  CD  GLU A 114      24.227  -7.119  -6.367  1.00 40.38           C  \n'
    'ANISOU    5  CD  GLU A 114     4393   3584   7366   -777  -1636   -513       C  \n'
    'TER   \n'
    'END\n')

@testing.requires('incentive')
@testing.requires_version('2.3.1')
class Test3221(testing.PyMOLTestCase):

    def test(self):
        cmd.read_pdbstr(v_pdbstr_anisou, 'm1')

        with testing.mktemp('.mae') as maefilename:
            cmd.save(maefilename, 'm1')
            cmd.delete('*')
            cmd.load(maefilename, 'm1')

        self.assertEqual(cmd.get_pdbstr('m1'), v_pdbstr_anisou)

        cmd.read_pdbstr(v_pdbstr_rotated, 'm2')
        cmd.align('m1', 'm2')

        with testing.mktemp('.mae') as maefilename:
            cmd.save(maefilename, 'm1')
            cmd.delete('*')
            cmd.load(maefilename, 'm1')

        # Attention: Not sure if rotation is numerically stable across
        # platforms. Fuzzy comparison might be needed.
        self.assertEqual(cmd.get_pdbstr('m1'), v_pdbstr_rotated)
