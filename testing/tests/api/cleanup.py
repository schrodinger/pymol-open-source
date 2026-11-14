"""
Structure cleanup: the "clean" submenu of the action menus
(pymol.menu.clean) removes water, ions and ligands.
"""

from pymol import cmd, menu, testing

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
            rows = {row[1][4:]: row[2] for row in menu.clean(cmd, "m")[1:]}
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
            self.assertEqual(len(submenu), 5)
