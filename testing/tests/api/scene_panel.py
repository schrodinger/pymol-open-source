"""
Scenes panel: rows follow the scene list and the current scene; recall, move,
rename (bad names refused), add, delete. Headless (offscreen QApplication),
where scenes have no thumbnail.
"""

import os

from pymol import cmd, testing


class TestScenePanel(testing.PyMOLTestCase):
    def setUp(self):
        super().setUp()
        os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
        from pymol.Qt import QtWidgets

        self.app = QtWidgets.QApplication.instance() or QtWidgets.QApplication([])

    def test_panel(self):
        from pmg_qt.scene_panel import ScenePanel

        cmd.fragment("gly")
        for name in ("a", "b", "c"):
            cmd.turn("x", 30)
            cmd.scene(name, "store")
        panel = ScenePanel(cmd)
        # without a GUI, quit crashes on a Qt widget that is still alive (PExit)
        self.addCleanup(lambda: (panel.deleteLater(), self.app.processEvents()))
        rows = panel.list

        def names():
            panel.refresh()
            return [rows.item(i).text() for i in range(rows.count())]

        self.assertEqual(names(), ["a", "b", "c"])
        self.assertEqual(rows.currentItem().text(), "c")

        rows.itemActivated.emit(rows.item(0))
        self.assertEqual(cmd.get("scene_current_name"), "a")
        names()
        self.assertEqual(rows.currentRow(), 0)

        panel._move(1)  # Ctrl+Down
        self.assertEqual(cmd.get_scene_list(), ["b", "a", "c"])
        self.assertEqual(rows.currentItem().text(), "a")
        panel._move(1)
        panel._move(1)  # already last
        panel._move(-1)
        self.assertEqual(cmd.get_scene_list(), ["b", "a", "c"])
        self.assertEqual(names(), ["b", "a", "c"])

        rows.item(0).setText("first")
        self.app.processEvents()
        self.assertEqual(cmd.get_scene_list(), ["first", "a", "c"])
        # "" would delete the scene, "a" would overwrite scene a
        for bad in ("", "a", "two words"):
            rows.item(0).setText(bad)
            self.app.processEvents()
            self.assertEqual(cmd.get_scene_list(), ["first", "a", "c"])
            self.assertEqual(names(), ["first", "a", "c"])

        panel._add()
        self.assertEqual(len(names()), 4)
        rows.setCurrentRow(1)
        panel._scene("clear")
        self.assertEqual(names()[:2], ["first", "c"])

        cmd.scene("*", "clear")  # a command, not the panel
        self.assertEqual(names(), [])
