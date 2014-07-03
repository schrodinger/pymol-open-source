'''
PYMOL-1241
PyMOL doesn't read last ANISOU record
'''

from pymol import cmd, testing, stored

v_pdbstr_anisou = (
    'ATOM      1  N   ARG A 197       5.287   1.830   7.079  1.00 30.85           N  \n'
    'ANISOU    1  N   ARG A 197     3652   4226   3841   -164    -54    -45       N  \n'
    'ATOM      2  CA  ARG A 197       6.677   2.258   6.809  1.00 28.67           C  \n'
    'ANISOU    2  CA  ARG A 197     3650   3727   3513   -123    -43    -41       C  \n'
    'ATOM      3  C   ARG A 197       7.070   3.179   7.964  1.00 26.05           C  \n'
    'ANISOU    3  C   ARG A 197     3219   3434   3242   -225    -95     22       C  \n'
    'END\n')

class TestPYMOL1241(testing.PyMOLTestCase):

    def test(self):
        cmd.read_pdbstr(v_pdbstr_anisou, 'm1')
        for a in cmd.get_model().atom:
            self.assertNotEqual(a.u_aniso, [0.0, 0.0, 0.0, 0.0, 0.0, 0.0], 'empty ANISOU')
