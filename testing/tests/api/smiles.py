"""
load_smiles (modules/pymol/smiles.py) and Build > From SMILES. Skipped
without rdkit, except the message that names the extra.
"""

import os
import sys
import tempfile
from unittest import mock

from pymol import CmdException, cmd, smiles, testing

try:
    import rdkit  # noqa: F401
except ImportError:
    rdkit = None

needs_rdkit = testing.unittest.skipUnless(rdkit, "rdkit not installed")


class TestSmiles(testing.PyMOLTestCase):
    def test_without_rdkit(self):
        with mock.patch.dict(sys.modules, {"rdkit": None}):
            with self.assertRaisesRegex(CmdException, r"pymolish\[chem\]"):
                cmd.load_smiles("CCO")
        self.assertIn("load_smiles", cmd.keyword)

    @needs_rdkit
    def test_one(self):
        self.assertEqual(cmd.load_smiles('"CC(=O)O"', "acid"), ["acid"])
        self.assertEqual(cmd.count_atoms("acid"), 8)  # with hydrogens
        self.assertEqual(cmd.count_atoms("acid and elem H"), 4)
        extent = cmd.get_extent("acid")
        self.assertGreater(extent[1][2] - extent[0][2], 0.5)  # 3D, not flat
        # the same input gives the same coordinates; a used name gets a number
        self.assertEqual(cmd.load_smiles("CC(=O)O", "acid"), ["acid01"])
        self.assertAlmostEqual(cmd.rms_cur("acid", "acid01"), 0, places=3)
        self.assertEqual(cmd.load_smiles("C#N"), ["mol01"])
        self.assertRaises(CmdException, cmd.load_smiles, "not a smiles")
        self.assertEqual(len(cmd.get_names()), 3)

    @needs_rdkit
    def test_conformers(self):
        cmd.load_smiles("CCCCCCO", "hexanol", 8, "g")
        self.assertGreater(cmd.count_states("hexanol"), 1)
        self.assertLessEqual(cmd.count_states("hexanol"), 8)
        self.assertEqual(cmd.count_atoms("hexanol"), 21)  # states, not copies
        self.assertEqual(cmd.get_object_list("g"), ["hexanol"])

    @needs_rdkit
    def test_file(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = os.path.join(tmp, "ligands.csv")
            with open(path, "w") as handle:
                handle.write(
                    "smiles,name\n# comment\nCCO,ethanol\n\n"
                    '"c1ccccc1";"benzene"\nXYZ,bad\nCO\nN\tb\n'
                )
            names = cmd.load_smiles(path, "lig", group="ligands")
        # "b" is a reserved word: PyMOL names the object "b_"
        self.assertEqual(names, ["ethanol", "benzene", "lig", "b_"])
        self.assertEqual(cmd.get_object_list("ligands"), names)
        self.assertEqual(cmd.count_atoms("benzene"), 12)

    @needs_rdkit
    def test_menu(self):
        os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
        from pymol.Qt import QtWidgets

        app = QtWidgets.QApplication.instance() or QtWidgets.QApplication([])
        from pmg_qt import file_dialogs

        app.clipboard().setText(" C=O ")
        window = mock.Mock(cmd=cmd)
        asked = {}

        def get_text(parent, title, label, text=""):
            asked["text"] = text
            return "C=C", True

        with mock.patch.object(QtWidgets.QInputDialog, "getText", get_text):
            file_dialogs.build_from_smiles(window)
        self.assertEqual(asked["text"], "C=O")  # from the clipboard
        self.assertEqual(cmd.count_atoms("mol01"), 6)  # ethene, typed
        app.clipboard().setText("two words")
        with mock.patch.object(QtWidgets.QInputDialog, "getText", get_text):
            file_dialogs.build_from_smiles(window)
        self.assertEqual(asked["text"], "")
