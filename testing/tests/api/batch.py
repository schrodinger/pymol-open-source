"""
load_files and export_objects (modules/pymol/batch.py) and their two File
dialogs, with real forms on an offscreen QApplication.
"""

import os
import pathlib
import tempfile

from pymol import CmdException, cmd, testing


class TestBatch(testing.PyMOLTestCase):
    def setUp(self):
        super().setUp()
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = pathlib.Path(self.tmp.name)
        (self.root / "sub").mkdir()
        cmd.fragment("ala", "tmp")
        for name in ("b.pdb", "a.pdb", "sub/a.pdb", "other.sdf", "UP.PDB"):
            cmd.save(str(self.root / name), "tmp")
        cmd.delete("tmp")
        (self.root / "junk.pdb").write_text("not a structure\n")

    def test_load_files(self):
        self.assertIn("load_files", cmd.keyword)
        self.assertIn("export_objects", cmd.keyword)
        names = cmd.load_files(str(self.root))
        # sorted, junk.pdb skipped; "b" is a reserved word and gets "_"
        self.assertEqual(names, ["UP", "a", "b_"])
        self.assertEqual(cmd.get_names(), names)
        # a name in use gets a number, nothing is replaced or appended to
        self.assertEqual(cmd.load_files(str(self.root), pattern="a"), ["a01"])
        self.assertEqual(cmd.count_states("a"), 1)
        cmd.delete("*")
        names = cmd.load_files(str(self.root), "pdb", "g", "[ab]", 1, "x_")
        self.assertEqual(names, ["x_a", "x_b", "x_sub_a"])
        self.assertEqual(cmd.get_object_list("g"), names)
        self.assertEqual(cmd.load_files(str(self.root), ".SDF"), ["other"])
        self.assertRaises(CmdException, cmd.load_files, str(self.root), "cif")
        self.assertRaises(CmdException, cmd.load_files, str(self.root / "no"))

    def test_export_objects(self):
        cmd.load_files(str(self.root), group="g")
        out = self.root / "out" / "deep"
        written = cmd.export_objects("g", str(out), "pdb")
        self.assertEqual(sorted(os.listdir(out)), ["UP.pdb", "a.pdb", "b_.pdb"])
        self.assertEqual(len(written), 3)
        # files that exist are kept unless overwrite=1
        (out / "a.pdb").write_text("mine")
        self.assertEqual(len(cmd.export_objects("g", str(out), "pdb")), 0)
        self.assertEqual((out / "a.pdb").read_text(), "mine")
        cmd.export_objects("a", str(out), "pdb", overwrite=1)
        self.assertIn("ATOM", (out / "a.pdb").read_text())
        # of each object only the selected atoms; any format "save" knows
        cmd.export_objects("(a or b_) and name CA", str(out), "xyz")
        self.assertEqual((out / "b_.xyz").read_text().split()[0], "1")
        self.assertFalse((out / "UP.xyz").exists())
        self.assertRaises(CmdException, cmd.export_objects, "none", str(out))

    def test_dialogs(self):
        os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
        from pymol.Qt import QtWidgets

        app = QtWidgets.QApplication.instance() or QtWidgets.QApplication([])
        from pmg_qt import file_dialogs
        from pmg_qt.pymol_qt_gui import PyMOLQtGUI

        class Window(QtWidgets.QWidget):
            load_form = PyMOLQtGUI.load_form
            initialdir = str(self.root / "sub")

        window = Window()
        window.cmd = cmd
        ok = QtWidgets.QDialogButtonBox.StandardButton.Ok

        form = file_dialogs.file_open_folder(window)
        self.assertEqual(form.output_count.text(), "1 files match")
        form.input_dir.setText(str(self.root / "missing"))
        self.assertEqual(form.output_count.text(), "Not a folder")
        self.assertFalse(form.buttons.button(ok).isEnabled())
        form.input_dir.setText(str(self.root))
        form.input_recursive.setChecked(True)
        form.input_filter.setText("a")
        form.input_group.setText("g")
        self.assertEqual(form.output_count.text(), "2 files match")
        self.assertEqual(
            form.output_command.text(),
            "load_files %s, pdb, g, pattern=a, recursive=1" % self.root,
        )
        form.buttons.button(ok).click()
        self.assertEqual(cmd.get_object_list("g"), ["a", "sub_a"])
        self.assertEqual(window.initialdir, str(self.root))

        form = file_dialogs.file_export_objects(window)
        self.assertEqual(form.input_selection.currentText(), "g")
        self.assertEqual(form.output_count.text(), "2 objects")
        form.input_selection.setEditText("no such thing")
        self.assertFalse(form.buttons.button(ok).isEnabled())
        form.input_selection.setEditText("sub_*")
        form.input_dir.setText(str(self.root / "out"))
        form.input_format.setEditText("pdb")
        self.assertEqual(
            form.output_command.text(),
            "export_objects sub_*, %s, pdb" % (self.root / "out"),
        )
        form.buttons.button(ok).click()
        self.assertEqual(os.listdir(self.root / "out"), ["sub_a.pdb"])
        del app
