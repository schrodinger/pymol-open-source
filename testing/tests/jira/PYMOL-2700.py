'''
mol2 partial charge export
'''

from pymol import cmd, CmdException, testing, stored

class TestPYMOL2700(testing.PyMOLTestCase):

    @testing.requires_version('1.8.0.4')
    def testMol2PartialChargeExport(self):
        cmd.fragment('gly')

        charges_pre = []
        cmd.iterate('*', 'charges_pre.append(partial_charge)', space=locals())

        # fragment library has charges, but better check...
        self.assertTrue(sum(abs(c) for c in charges_pre) > 0.0)

        with testing.mktemp('.mol2') as filename:
            cmd.save(filename)
            cmd.delete('*')
            cmd.load(filename)

        charges_post = []
        cmd.iterate('*', 'charges_post.append(partial_charge)', space=locals())

        self.assertArrayEqual(charges_pre, charges_post, delta=1e-2)
