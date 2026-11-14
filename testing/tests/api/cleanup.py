"""
Structure cleanup: the "clean" submenu of the action menus
(pymol.menu.clean) removes water, ions and ligands; "clean" minimizes a
small molecule with rdkit; "fix_structure" repairs a protein with PDBFixer.
"""

import sys
from unittest import mock

from pymol import CmdException, cmd, computing, menu, testing

try:
    import rdkit  # noqa: F401
except ImportError:
    rdkit = None

needs_rdkit = testing.unittest.skipUnless(rdkit, "rdkit not installed")

try:
    import pdbfixer  # noqa: F401
except ImportError:
    pdbfixer = None

needs_pdbfixer = testing.unittest.skipUnless(pdbfixer, "pdbfixer not installed")

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
            self.assertEqual(len(submenu), 8)


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


class TestFixStructure(testing.PyMOLTestCase):
    def test_without_pdbfixer(self):
        cmd.fragment("ala", "m")
        with mock.patch.dict(sys.modules, {"pdbfixer": None}):
            with self.assertRaisesRegex(CmdException, r"pymolish\[fix\]"):
                cmd.fix_structure("m")
        self.assertIn("fix_structure", cmd.keyword)

    @needs_pdbfixer
    def test_missing_atoms(self):
        cmd.fab("AKA", "p")
        cmd.remove("hydro")
        heavy = cmd.count_atoms("p")
        cmd.remove("p and resi 2 and name NZ+CE")
        self.assertEqual(cmd.fix_structure("p"), "p_fixed")
        # the two side-chain atoms are back, and the terminal oxygen is new
        self.assertEqual(cmd.count_atoms("p_fixed and not hydro"), heavy + 1)
        self.assertEqual(cmd.count_atoms("p_fixed and resi 2 and name NZ+CE"), 2)
        self.assertGreater(cmd.count_atoms("p_fixed and hydro"), 20)
        self.assertEqual(cmd.count_atoms("p"), heavy - 2)  # the original stays
        # nothing is overwritten; options off leave the atoms as they are
        name = cmd.fix_structure("p", missing_atoms=0, hydrogens=0)
        self.assertEqual(name, "p_fixed01")
        self.assertEqual(cmd.count_atoms(name), heavy - 2)
        self.assertEqual(cmd.fix_structure("p", name="mine", hydrogens=0), "mine")
        self.assertRaises(CmdException, cmd.fix_structure, "p or mine")
        self.assertRaises(CmdException, cmd.fix_structure, "none")

    @needs_pdbfixer
    def test_hydrogens_by_ph(self):
        cmd.fab("AHA", "p")
        cmd.remove("hydro")
        low = cmd.fix_structure("p", ph=4.0, name="low")
        high = cmd.fix_structure("p", ph=9.0, name="high")
        self.assertEqual(
            cmd.count_atoms(low + " and resn HIS and hydro"),
            cmd.count_atoms(high + " and resn HIS and hydro") + 1,
        )

    @needs_pdbfixer
    def test_nonstandard_and_heterogens(self):
        cmd.read_pdbstr(PDB, "m")
        cmd.alter("m and resi 1", "resn='DAL'")  # D-alanine becomes ALA
        cmd.fix_structure("m", hydrogens=0, name="full")
        self.assertEqual(cmd.count_atoms("full and resn DAL"), 0)
        self.assertEqual(cmd.count_atoms("full and resn ALA and name CA"), 1)
        self.assertEqual(cmd.count_atoms("full and solvent"), 2)
        self.assertEqual(cmd.count_atoms("full and resn EOH+NA"), 4)
        cmd.fix_structure("m", hydrogens=0, nonstandard=0, keep_water=0, name="dry")
        self.assertEqual(cmd.count_atoms("dry and resn DAL and name CA"), 1)
        self.assertEqual(cmd.count_atoms("dry and solvent"), 0)
        self.assertEqual(cmd.count_atoms("dry and resn EOH+NA"), 4)
        cmd.fix_structure("m", hydrogens=0, keep_heterogens=0, name="bare")
        self.assertEqual(cmd.count_atoms("bare and resn EOH+NA"), 0)
        self.assertEqual(cmd.count_atoms("bare and solvent"), 2)

    @needs_pdbfixer
    def test_dialog(self):
        import os

        os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
        from pymol.Qt import QtWidgets

        app = QtWidgets.QApplication.instance() or QtWidgets.QApplication([])
        from pmg_qt import file_dialogs
        from pmg_qt.pymol_qt_gui import PyMOLQtGUI

        class Window(QtWidgets.QWidget):
            load_form = PyMOLQtGUI.load_form

        window = Window()
        window.cmd = cmd
        cmd.fab("AKA", "p")
        cmd.remove("hydro or (resi 2 and name NZ)")
        row = [r for r in menu.clean(cmd, "p") if r[1] == "fix structure..."][0]
        self.assertIn('fix_structure_dialog("p")', row[2])

        form = file_dialogs.fix_structure_dialog(window, "p")
        ok = form.buttons.button(QtWidgets.QDialogButtonBox.StandardButton.Ok)
        self.assertEqual(form.output_command.text(), "fix_structure p, async_=1")
        form.input_object.setEditText("nothing here")
        self.assertFalse(ok.isEnabled())
        form.input_object.setEditText("p")
        form.input_ph.setValue(5.5)
        form.input_water.setChecked(False)
        self.assertEqual(
            form.output_command.text(),
            "fix_structure p, keep_water=0, ph=5.5, async_=1",
        )
        form.input_hydrogens.setChecked(False)
        self.assertFalse(form.input_ph.isEnabled())
        self.assertEqual(
            form.output_command.text(),
            "fix_structure p, hydrogens=0, keep_water=0, async_=1",
        )
        # the background thread needs a running main loop: call directly
        with mock.patch.object(cmd, "async_", lambda func, *a, **k: func(*a, **k)):
            ok.click()
        self.assertEqual(cmd.count_atoms("p_fixed and resi 2 and name NZ"), 1)
        self.assertEqual(cmd.count_atoms("p_fixed and hydro"), 0)
        del app
