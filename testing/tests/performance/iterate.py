'''
Stress testing for iterate/alter
'''

from pymol import cmd, testing, stored

class StressIterateAlter(testing.PyMOLTestCase):

    def load_big_example_multistate(self):
        # 4434 atoms in 61 states (270474 atoms)
        cmd.load(self.datafile('2cas.pdb.gz'))
        cmd.load(self.datafile('2cas.dcd'))

    def testAlterState(self):
        self.load_big_example_multistate()

        v_count = cmd.count_atoms() * cmd.count_states()
        assert v_count > 10**5

        xyz = []
        with self.timing('i', 5.0):
            cmd.iterate_state(0, 'all', 'xyz.append((x,y,z))', space=locals())

        self.assertEqual(v_count, len(xyz))

        with self.timing('a', 5.0):
            cmd.alter_state(0, 'all', '(x,y,z) = next(xyz_rev)',
                    space={'xyz_rev': reversed(xyz), 'next': next})

        cmd.iterate_state(0, 'last all', 'stored.xyz = (x,y,z)')
        self.assertEqual(stored.xyz, xyz[0])
