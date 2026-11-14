"""
search_sequence (modules/pymol/finding.py), cluster_ligands and
color_group_with_gradient (modules/pymol/grouping.py).
"""

from pymol import CmdException, cmd, grouping, testing


def residues(selection):
    found = []
    cmd.iterate(
        "(%s) and name CA" % selection,
        "found.append((model, chain, int(resi)))",
        space={"found": found},
    )
    return found


class TestFinding(testing.PyMOLTestCase):
    def test_search_sequence(self):
        cmd.fab("ACDNGSACD", "p1", chain="A")
        cmd.fab("GSACD", "p2", chain="B")
        cmd.alter("p2", "resi=str(int(resi) + 100)")
        self.assertEqual(cmd.search_sequence("acd", "hit"), 3)
        self.assertEqual(
            residues("hit"),
            [("p1", "A", n) for n in (1, 2, 3, 7, 8, 9)]
            + [("p2", "B", n) for n in (103, 104, 105)],
        )
        # whole residues are selected, not only the CA atoms
        self.assertEqual(
            cmd.count_atoms("hit and p2"), cmd.count_atoms("p2 and resi 103-105")
        )
        # a regular expression, the default name, a selection to look in
        self.assertEqual(cmd.search_sequence("N[^P][ST]", selection="p1"), 1)
        self.assertEqual(residues("seq_NPST"), [("p1", "A", n) for n in (4, 5, 6)])
        # no match across two chains, an empty selection when nothing is found
        self.assertEqual(cmd.search_sequence("CDGS", "none_found"), 0)
        self.assertEqual(cmd.count_atoms("none_found"), 0)
        self.assertRaises(CmdException, cmd.search_sequence, "A[")
        self.assertIn("search_sequence", cmd.keyword)

    def ligands(self):
        """Six one-atom objects: three near the origin, two near x=50, one
        alone at x=200."""
        for name, x in zip("abcdef", (0, 3, 6, 50, 52, 200)):
            cmd.pseudoatom("lig_" + name, pos=[x, 0, 0])
        cmd.group("ligs", "lig_*")

    def test_cluster_ligands(self):
        self.ligands()
        groups = cmd.cluster_ligands("lig_*")
        self.assertEqual(groups, ["cluster_1", "cluster_2", "cluster_3"])
        self.assertEqual(grouping.members("cluster_1"), ["lig_a", "lig_b", "lig_c"])
        self.assertEqual(grouping.members("cluster_2"), ["lig_d", "lig_e"])
        self.assertEqual(grouping.members("cluster_3"), ["lig_f"])
        # a chain of small steps is one cluster: 0, 3, 6 with a cutoff of 4
        self.assertEqual(len(cmd.cluster_ligands("lig_*", cutoff=4)), 3)
        # a number of groups instead of a cutoff; emptied groups go away
        self.assertEqual(cmd.cluster_ligands("lig_*", 2), ["cluster_1", "cluster_2"])
        self.assertEqual(len(grouping.members("cluster_1")), 5)
        self.assertNotIn("cluster_3", cmd.get_names())
        cmd.cluster_ligands("lig_*", prefix="site", cutoff=1000)
        self.assertEqual(len(grouping.members("site_1")), 6)
        self.assertRaises(CmdException, cmd.cluster_ligands, "lig_a")
        self.assertIn("cluster_ligands", cmd.keyword)

    def test_color_group_with_gradient(self):
        self.ligands()
        cmd.color_group_with_gradient("ligs", "red")
        colors = []
        cmd.iterate("ligs", "colors.append(color)", space={"colors": colors})
        rgb = [cmd.get_color_tuple(c) for c in colors]
        reds = [c[0] for c in rgb]
        self.assertEqual(reds, sorted(reds))  # dark to light
        self.assertAlmostEqual(reds[0], 0.3, places=2)
        self.assertAlmostEqual(reds[-1], 1.0, places=2)
        self.assertTrue(all(c[1] == 0 and c[2] == 0 for c in rgb))
        # a second gradient leaves the first one alone
        cmd.color_group_with_gradient("lig_a", "blue")
        self.assertAlmostEqual(
            cmd.get_color_tuple(self.color_of("lig_f"))[0], 1.0, places=2
        )
        self.assertRaises(CmdException, cmd.color_group_with_gradient, "none")
        self.assertIn("color_group_with_gradient", cmd.keyword)

    def color_of(self, name):
        colors = []
        cmd.iterate(name, "colors.append(color)", space={"colors": colors})
        return colors[0]
