from pymol import cmd, testing, stored
from chempy import cpv

class TestFitting(testing.PyMOLTestCase):

    def testAlign(self):
        cmd.load(self.datafile("1oky-frag.pdb"), "m1")
        cmd.load(self.datafile("1t46-frag.pdb"), "m2")
        cycles = 2
        r = cmd.align("m1", "m2", object="aln", cycles=cycles, transform=0)
        self.assertAlmostEqual(r[0], 1.52563, delta=1e-4)
        self.assertEqual(r[1], 177)
        self.assertEqual(r[1], cmd.count_atoms("aln") / 2)
        self.assertEqual(r[2], cycles)

    def testAlignto(self):
        cmd.fragment("gly", "m1")
        cmd.copy("m2", "m1")
        cmd.alignto(method="fit", mobile_state=1, target_state=1)

    def testCealign(self):
        cmd.load(self.datafile("1oky-frag.pdb"), "m1")
        cmd.load(self.datafile("1t46-frag.pdb"), "m2")
        r = cmd.cealign("m2", "m1", object="aln")
        self.assertAlmostEqual(r["RMSD"], 1.90375, delta=1e-4)
        alen = r["alignment_length"]
        self.assertEqual(alen, 40)
        self.assertEqual(alen, cmd.count_atoms("aln") / 2)

    def testFit(self):
        cmd.fragment("gly", "m1")
        cmd.create("m2", "m1")
        rms = cmd.fit("m1", "m2")
        self.assertEqual(rms, 0.0)
        rms = cmd.rms("m1", "m2")
        self.assertEqual(rms, 0.0)
        rms = cmd.rms_cur("m1", "m2")
        self.assertEqual(rms, 0.0)

    def testIntraFit(self):
        cmd.fragment("gly", "m1")
        cmd.create("m1", "m1", 1, 2)
        rms_list = cmd.intra_fit("m1")
        self.assertArrayEqual(rms_list, [-1.0, 0.0])
        rms_list = cmd.intra_rms("m1")
        self.assertArrayEqual(rms_list, [-1.0, 0.0])
        rms_list = cmd.intra_rms_cur("m1")
        self.assertArrayEqual(rms_list, [-1.0, 0.0])

    def testIntraRms(self):
        # see intra_fit
        pass

    def testIntraRmsCur(self):
        # see intra_fit
        pass

    def testPairFit(self):
        cmd.fragment('trp')
        cmd.fragment('his')

        # 1 atom
        sele = ('trp and guide', 'his and guide')
        pos = map(cmd.get_atom_coords, sele)
        vec = cpv.sub(*pos)
        mat_ref = [
            1.0, 0.0, 0.0, -vec[0],
            0.0, 1.0, 0.0, -vec[1],
            0.0, 0.0, 1.0, -vec[2],
            0.0, 0.0, 0.0, 1.0]
        rms = cmd.pair_fit(*sele)
        self.assertEqual(rms, 0.0)
        mat = cmd.get_object_matrix('trp')
        self.assertArrayEqual(mat, mat_ref, 1e-4)

        # 2 atoms
        sele += ('trp & name CB', 'his & name CB')
        rms = cmd.pair_fit(*sele)
        self.assertAlmostEqual(rms, 0.0082, delta=1e-4)

        # 4 atoms
        sele += ('trp & name CG', 'his & name CG',
                 'trp & name CD1', 'his & name CD2')
        rms = cmd.pair_fit(*sele)
        self.assertAlmostEqual(rms, 0.0713, delta=1e-4)

    def testRms(self):
        # see fit
        pass

    def testRmsCur(self):
        # see fit
        pass

    def testSuper(self):
        cmd.load(self.datafile("1oky-frag.pdb"), "m1")
        cmd.load(self.datafile("1t46-frag.pdb"), "m2")
        r = cmd.super("m1", "m2", object="aln")
        self.assertAlmostEqual(r[0], 0.9667, delta=1e-4)
        self.assertEqual(r[1], 172)
        self.assertEqual(r[1], cmd.count_atoms("aln") / 2)

