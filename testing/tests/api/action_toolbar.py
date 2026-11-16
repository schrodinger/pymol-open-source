"""
Action toolbar: selection mode menu and label, rock state, zoom, presets
menu, panel and More entries, an icon on every button, the palette menu, the
UI scale menu.
Headless (offscreen QApplication).
"""

import os
from collections import defaultdict

from pymol import cmd, testing


class TestActionToolbar(testing.PyMOLTestCase):
    def setUp(self):
        super().setUp()
        os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
        from pymol.Qt import QtWidgets

        self.app = QtWidgets.QApplication.instance() or QtWidgets.QApplication([])

    def test_toolbar(self):
        from pmg_qt.action_toolbar import ActionToolbar

        hits = []
        callbacks = defaultdict(list)
        bar = ActionToolbar(
            cmd,
            callbacks,
            panels=[("Builder", lambda: hits.append("builder"))],
            more=[("Get View", lambda: hits.append("view"))],
        )
        texts = {a.text(): a for a in bar.actions()}

        # selection mode: menu sets the setting, the setting sets the label
        self.assertEqual(bar.selecting.text(), "Residues")
        bar.selecting.menu().actions()[2].trigger()
        self.assertEqual(cmd.get_setting_int("mouse_selection_mode"), 2)
        self.assertEqual(bar.selecting.text(), "Chains")
        for callback in callbacks[cmd.setting._get_index("mouse_selection_mode")]:
            callback(0)
        self.assertEqual(bar.selecting.text(), "Atoms")

        # rock follows the setting
        self.assertFalse(bar.rock.isChecked())
        for callback in callbacks[cmd.setting._get_index("rock")]:
            callback(1)
        self.assertTrue(bar.rock.isChecked())

        cmd.fragment("ala")
        view = cmd.get_view()
        cmd.turn("x", 90)
        cmd.move("z", -30)
        bar.zoom.click()
        texts["Orient"].trigger()
        self.assertNotEqual(cmd.get_view(), view)

        bar.presets.menu().aboutToShow.emit()
        self.assertIn("simple", [a.text() for a in bar.presets.menu().actions()])

        texts["Builder"].trigger()
        more = {a.text(): a for a in bar.more.menu().actions()}
        more["Get View"].trigger()
        self.assertEqual(hits, ["builder", "view"])

    def test_icons(self):
        from pymol.Qt import QtCore, QtWidgets
        from pmg_qt.action_toolbar import ActionToolbar

        bar = ActionToolbar(
            cmd,
            defaultdict(list),
            panels=[
                ("Builder", lambda: None),
                ("Scenes", lambda: None),
                ("Draw/Ray", QtWidgets.QMenu()),
            ],
        )
        buttons = [
            b
            for b in bar.findChildren(QtWidgets.QToolButton)
            if b.objectName() != "qt_toolbar_ext_button"  # Qt's overflow arrow
        ]
        self.assertEqual(len(buttons), 13)
        self.assertEqual([b.text() for b in buttons if b.icon().isNull()], [])
        IconOnly = QtCore.Qt.ToolButtonStyle.ToolButtonIconOnly
        self.assertEqual(
            sorted(b.toolTip() for b in buttons if b.toolButtonStyle() == IconOnly),
            ["Colors that new objects get", "More view actions", "Redo", "Undo"],
        )

    def test_palette(self):
        from unittest import mock
        from pymol import palettes
        from pmg_qt import action_toolbar, theme

        # -k (how the tests run): the saved palette is not applied
        with mock.patch.object(theme, "object_palette", return_value="deep"):
            bar = action_toolbar.ActionToolbar(cmd, defaultdict(list))
        self.assertEqual(palettes.active, "default")
        menu = bar.colors.menu()
        actions = {a.data(): a for a in menu.actions() if a.data()}
        self.assertEqual(list(actions), ["default", *palettes.PALETTES])
        self.assertEqual(actions["deep"].text(), "Deep (white background)")
        self.assertFalse(actions["deep"].icon().isNull())
        try:
            with mock.patch.object(theme, "set_object_palette") as save:
                actions["deep"].trigger()
            save.assert_called_once_with("deep")
            self.assertEqual(palettes.active, "deep")

            cmd.palette("soft")  # typed: the menu follows when it opens
            menu.aboutToShow.emit()
            self.assertEqual([n for n, a in actions.items() if a.isChecked()], ["soft"])

            cmd.fragment("ala")
            menu.actions()[-1].trigger()  # Recolor objects now
            self.assertEqual(cmd.count_atoms("elem C & color soft_01"), 3)

            # without -k the saved palette is back at start
            cmd.palette("default")
            with (
                mock.patch.object(theme, "object_palette", return_value="deep"),
                mock.patch.object(
                    action_toolbar.pymol.invocation.options, "pymolrc", []
                ),
            ):
                action_toolbar.ActionToolbar(cmd, defaultdict(list))
            self.assertEqual(palettes.active, "deep")
        finally:
            cmd.palette("default")

    def test_scale(self):
        from unittest import mock
        from pmg_qt import action_toolbar, theme

        with mock.patch.object(theme, "ui_scale", return_value=1.8):
            bar = action_toolbar.ActionToolbar(cmd, defaultdict(list))
        self.assertEqual(bar.scale.text(), "180%")
        actions = bar.scale.menu().actions()
        self.assertEqual(
            [a.text() for a in actions], ["100%", "120%", "150%", "180%", "200%"]
        )
        self.assertEqual([a.text() for a in actions if a.isChecked()], ["180%"])
        with (
            mock.patch.object(theme, "set_ui_scale") as save,
            mock.patch.object(action_toolbar.QtWidgets.QMessageBox, "information"),
        ):
            actions[1].trigger()
        save.assert_called_once_with(1.2)
