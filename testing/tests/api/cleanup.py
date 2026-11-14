"""
Structure cleanup: the "clean" submenu of the action menus
(pymol.menu.clean) removes water, ions and ligands; "clean" minimizes a
small molecule with rdkit.
"""

import sys
from unittest import mock

from pymol import CmdException, cmd, computing, menu, testing

try:
    import rdkit  # noqa: F401
except ImportError:
    rdkit = None

needs_rdkit = testing.unittest.skipUnless(rdkit, "rdkit not installed")

# a dipeptide, two waters, a sodium ion and an ethanol ligand
PDB = """\
ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00  0.00           N
ATOM      2  CA  ALA A   1       1.458   0.000   0.000  1.00  0.00           C
ATOM      3  C   ALA A   1       2.009   1.420   0.000  1.00  0.00           C
ATOM      4  O   ALA A   1       1.251   2.390   0.000  1.00  0.00           O
ATOM      5  CB  ALA A   1       1.988  -0.773  -1.199  1.00  0.00           C
ATOM      6  N   GLY A   2       3.332   1.536   0.000  1.00  0.00           N
ATOM      7  CA  GLY A   2       3.988   2.839   0.000  1.00  0.00           C
ATOM      8  C   GLY A   2       5.504   2.693   0.000  1.00  0.00           C
ATOM      9  O   GLY A   2       6.057   1.593   0.000  1.00  0.00           O
HETATM   10  O   HOH A 101      10.000  10.000  10.000  1.00  0.00           O
HETATM   11  O   HOH A 102      12.000  10.000  10.000  1.00  0.00           O
HETATM   12 NA    NA A 201      15.000  10.000  10.000  1.00  0.00          NA
HETATM   13  C1  EOH A 301      20.000  10.000  10.000  1.00  0.00           C
HETATM   14  C2  EOH A 301      21.500  10.000  10.000  1.00  0.00           C
HETATM   15  O   EOH A 301      22.000  11.300  10.000  1.00  0.00           O
END
"""
PARTS = {"polymer": 9, "solvent": 2, "inorganic": 1, "organic": 3}


class TestCleanup(testing.PyMOLTestCase):
    def counts(self):
        return {part: cmd.count_atoms("m and " + part) for part in PARTS}

    def test_clean_menu(self):
        gone = {
            "remove waters": ["solvent"],
            "remove ions": ["inorganic"],
            "remove waters and ions": ["solvent", "inorganic"],
            "remove ligands": ["organic"],
        }
        for label, parts in gone.items():
            cmd.delete("*")
            cmd.read_pdbstr(PDB, "m")
            cmd.fragment("ala", "other")  # outside the menu's target
            self.assertEqual(self.counts(), PARTS)
            rows = {row[1][4:]: row[2] for row in menu.clean(cmd, "m")[1:5]}
            cmd.do(rows[label], echo=0)
            expected = {p: 0 if p in parts else n for p, n in PARTS.items()}
            self.assertEqual(self.counts(), expected, label)
            self.assertEqual(cmd.count_atoms("other"), 10)

    def test_action_menus(self):
        cmd.read_pdbstr(PDB, "m")
        cmd.group("g", "m")
        for entries in (
            menu.mol_action(cmd, "m"),
            menu.sele_action(cmd, "m"),
            menu.group_action(cmd, "g"),
            menu.all_action(cmd, "all"),
        ):
            labels = [row[1] for row in entries]
            self.assertEqual(labels.count("clean"), 1)
            # "remove waters" lives in the submenu only
            self.assertFalse([text for text in labels if "waters" in text])
            submenu = entries[labels.index("clean")][2]
            self.assertEqual(len(submenu), 7)


class TestClean(testing.PyMOLTestCase):
    def test_without_rdkit(self):
        cmd.fragment("methane", "m")
        with mock.patch.dict(sys.modules, {"rdkit": None}):
            with self.assertRaisesRegex(CmdException, r"pymolish\[chem\]"):
                cmd.clean("m")

    @needs_rdkit
    def test_ethane(self):
        cmd.load_smiles("CC", "ethane")
        coords = cmd.get_coords("ethane")
        # pull one methyl group 0.5 A away along the bond
        axis = coords[1] - coords[0]
        cmd.translate(
            list(0.5 * axis / (axis @ axis) ** 0.5),
            "ethane and (index 2 or (hydro and neighbor index 2))",
            camera=0,
        )
        self.assertGreater(cmd.get_distance("index 1", "index 2"), 1.9)
        cmd.clean("ethane")
        self.assertAlmostEqual(cmd.get_distance("index 1", "index 2"), 1.53, delta=0.02)
        self.assertEqual(cmd.count_atoms("ethane"), 8)

    @needs_rdkit
    def test_without_hydrogens(self):
        # heavy atoms only: hydrogens are added for the run, not to the object
        cmd.load_smiles("CCO", "m")
        cmd.remove("hydro")
        cmd.alter_state(1, "m and index 1", "x += 0.4")
        before = cmd.get_coords("m")
        cmd.clean("m")
        self.assertEqual(cmd.count_atoms("m"), 3)
        self.assertGreater(abs(cmd.get_coords("m") - before).max(), 0.05)
        self.assertAlmostEqual(cmd.get_distance("index 1", "index 2"), 1.52, delta=0.03)

    @needs_rdkit
    def test_states(self):
        cmd.load_smiles("CCCCO", "m", 3)
        self.assertGreater(cmd.count_states("m"), 1)
        first = cmd.get_coords("m", 1).copy()
        cmd.alter_state(2, "m and index 1", "x += 0.4")
        moved = cmd.get_coords("m", 2).copy()
        cmd.clean("m", state=2)
        self.assertTrue((cmd.get_coords("m", 1) == first).all())
        self.assertFalse((cmd.get_coords("m", 2) == moved).all())

    @needs_rdkit
    def test_refused(self):
        cmd.read_pdbstr(PDB, "m")
        cmd.fragment("methane", "other")
        before = cmd.get_coords("all")
        for args, kwargs in (
            (("m and resi 1",), {}),  # a part of the peptide
            (("organic or other",), {}),  # two objects
            (("none",), {}),
            (("other",), {"present": "m"}),
            (("other",), {"fix": "elem C"}),
            (("other",), {"restrain": "elem C"}),
            (("other",), {"method": "uff"}),
        ):
            self.assertRaises(CmdException, cmd.clean, *args, **kwargs)
        with mock.patch.object(computing, "CLEAN_MAX_ATOMS", 4):
            with self.assertRaisesRegex(CmdException, "minimize"):
                cmd.clean("other")
        self.assertTrue((cmd.get_coords("all") == before).all())
        cmd.clean("m and organic")  # the ligand alone is a whole molecule
