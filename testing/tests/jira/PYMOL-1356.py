'''
PYMOL-1356
rotate ANISOU when aligning molecules
'''

from pymol import cmd, testing, stored

v_pdbstr_anisou = (
    'ATOM      1  N   GLU A 114      24.832  -7.270  -5.728  1.00 33.91           N  \n'
    'ANISOU    1  N   GLU A 114     6968   3709   2207   -518    495    146       N  \n'
    'ATOM      2  CA  GLU A 114      25.839  -6.416  -5.102  1.00 33.68           C  \n'
    'ANISOU    2  CA  GLU A 114     6957   3558   2282   -477    534    166       C  \n'
    'ATOM      3  CB  GLU A 114      25.893  -5.047  -5.787  1.00 34.17           C  \n'
    'ANISOU    3  CB  GLU A 114     7123   3591   2269   -451    534    251       C  \n'
    'ATOM      4  CG  GLU A 114      26.342  -5.083  -7.241  1.00 38.31           C  \n'
    'ANISOU    4  CG  GLU A 114     7707   4148   2699   -531    590    255       C  \n'
    'ATOM      5  CD  GLU A 114      25.207  -5.382  -8.202  1.00 40.38           C  \n'
    'ANISOU    5  CD  GLU A 114     7970   4550   2823   -523    543    278       C  \n'
    'END\n')

v_pdbstr_rotated = (
    'ATOM      1  N   GLU A 114      26.964  -6.052  -5.274  1.00 33.91           N  \n'
    'ANISOU    1  N   GLU A 114     3700   2782   6402   -635  -1526   -432       N  \n'
    'ATOM      2  CA  GLU A 114      26.822  -5.036  -6.315  1.00 33.68           C  \n'
    'ANISOU    2  CA  GLU A 114     3593   2814   6390   -532  -1532   -461       C  \n'
    'ATOM      3  CB  GLU A 114      25.355  -4.880  -6.726  1.00 34.17           C  \n'
    'ANISOU    3  CB  GLU A 114     3539   2892   6552   -584  -1566   -491       C  \n'
    'ATOM      4  CG  GLU A 114      24.745  -6.113  -7.378  1.00 38.31           C  \n'
    'ANISOU    4  CG  GLU A 114     4106   3362   7086   -637  -1669   -512       C  \n'
    'ATOM      5  CD  GLU A 114      24.227  -7.118  -6.367  1.00 40.38           C  \n'
    'ANISOU    5  CD  GLU A 114     4393   3584   7366   -777  -1636   -512       C  \n'
    'END\n')

views = [
  ( -0.119232267,    0.595452785,   -0.794493258,\
    -0.545022011,    0.629605889,    0.553667188,\
     0.829900384,    0.499031276,    0.249465555,\
     0.000000000,    0.000000000,  -14.178204536,\
    25.622600555,   -5.839600086,   -6.412000179,\
    11.178204536,   17.178203583,  -20.000000000 ),
  (  0.847085655,   -0.119858012,   -0.517764270,\
     0.405014724,    0.776385665,    0.482895762,\
     0.344105840,   -0.618756235,    0.706209421,\
    -0.000000000,    0.000000000,  -14.178204536,\
    25.622600555,   -5.839799881,   -6.412000179,\
    11.178204536,   17.178203583,  -20.000000000 ),
]

class Test1356(testing.PyMOLTestCase):

    def _rep(self, name='m1'):
        cmd.disable()
        cmd.color('yellow')
        cmd.show_as('ellipsoids')
        cmd.enable(name)

    @testing.foreach(0, 1)
    def test(self, matrix_mode):
        cmd.viewport(100,100)
        cmd.set('matrix_mode', matrix_mode)
        cmd.set('ambient', 1.0)

        cmd.read_pdbstr(v_pdbstr_anisou, 'm1', zoom=0)
        cmd.read_pdbstr(v_pdbstr_rotated, 'm2', zoom=0)
        self._rep()

        cmd.set_view(views[0])
        ref = self.get_imagearray()
        self.assertImageHasColor('yellow')

        cmd.align('m1', 'm2')
        cmd.set_view(views[1])
        self.assertImageEqual(ref, count=100, delta=1, msg='ANISOU not aligned on display')

        pdbstr = cmd.get_pdbstr('m1')
        cmd.delete('*')
        cmd.read_pdbstr(pdbstr, 'm1', zoom=0)
        self._rep()
        self.assertImageEqual(ref, count=100, msg='ANISOU not aligned in PDB')
