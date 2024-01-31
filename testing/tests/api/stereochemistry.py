import os
from pymol import movie, cmd, testing, stored

@testing.requires('incentive')
@testing.requires_version('2.2')
class TestStereochemistry(testing.PyMOLTestCase):

    def test_assign_stereo_schrodinger(self):
        if not os.getenv("SCHRODINGER"):
            self.skipTest("SCHRODINGER not set")

        self.assertEqual(cmd.count_atoms("stereo R"), 0)

        cmd.load(self.datafile("1ehz-5.pdb"), "m1")
        cmd.assign_stereo(method="schrodinger")

        self.assertEqual(cmd.count_atoms("stereo R"), 15)
        self.assertEqual(cmd.count_atoms("stereo S"), 5)

    def test_assign_stereo_rdkit(self):
        try:
            import rdkit
        except ImportError:
            self.skipTest("rdkit not available")

        self.assertEqual(cmd.count_atoms("stereo R"), 0)

        cmd.load(self.datafile("1ehz-5.pdb"), "m1")
        cmd.assign_stereo(method="rdkit")

        self.assertEqual(cmd.count_atoms("stereo R & not elem P"), 15)
        self.assertEqual(cmd.count_atoms("stereo S"), 5)
