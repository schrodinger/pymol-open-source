"""
Objects panel: rows and group nesting, the menu table against pymol.menu,
menu building, visibility toggle, filter, refresh after delete, the states
panel and the wizard buttons. Headless (offscreen QApplication).
"""

import os
import re

from pymol import cmd, menu, testing


class TestObjectPanel(testing.PyMOLTestCase):
    def setUp(self):
        super().setUp()
        os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
        from pymol.Qt import QtWidgets

        self.app = QtWidgets.QApplication.instance() or QtWidgets.QApplication([])

    def test_snapshot(self):
        from pmg_qt.object_panel import snapshot

        cmd.fragment("gly", "g1")
        cmd.fragment("ala", "a2")
        cmd.disable("a2")
        cmd.group("grp", "a2")
        cmd.select("sele", "name CA")
        rows = {r[0]: r for r in snapshot(cmd)}
        self.assertEqual(next(iter(rows)), "all")
        self.assertEqual(rows["a2"], ("a2", "object:molecule", False, "grp"))
        self.assertEqual(rows["g1"], ("g1", "object:molecule", True, ""))
        self.assertEqual(rows["grp"][1], "object:group")
        self.assertEqual(rows["sele"][1], "selection")

    def test_menu_table(self):
        from pmg_qt.object_panel import MENUS, menu_entries

        for functions in MENUS.values():
            self.assertEqual(len(functions), 5)
            for function in filter(None, functions):
                self.assertTrue(callable(getattr(menu, function)), function)
        cmd.fragment("ala")
        for button in range(5):
            self.assertTrue(menu_entries(cmd, "ala", "object:molecule", button))
        self.assertIsNone(menu_entries(cmd, "x", "object:map", 3))  # no labels

    def test_fill_menu(self):
        from pymol.Qt import QtWidgets
        from pmg_qt.object_panel import fill_menu

        cmd.fragment("ala")
        popup = QtWidgets.QMenu()
        entries = [
            [2, "Color:", ""],
            [1, "\\900red", 'cmd.color("red","ala")'],
            [0, "", ""],
            [1, "more", lambda: [[1, "hide", 'cmd.disable("ala")']]],
        ]
        entries.append([1, r"\933delete", 'cmd.delete("none")'])
        fill_menu(popup, entries, cmd)
        title = popup.actions()[0]  # a caption, skipped by the keyboard
        self.assertEqual(title.defaultWidget().text(), "Color:")
        self.assertFalse(title.isEnabled())
        delete = popup.actions()[-1]  # a known command gets its line icon
        self.assertFalse(delete.icon().isNull())
        red, more = (a for a in popup.actions() if a.text() in ("red", "more"))
        self.assertFalse(red.icon().isNull())  # swatch from the \900 markup
        red.trigger()
        colors = set()
        cmd.iterate("ala", "colors.add(color)", space={"colors": colors})
        self.assertEqual(colors, {cmd.get_color_index("red")})
        self.assertTrue(more.menu().isEmpty())  # filled when first opened
        more.menu().aboutToShow.emit()
        self.assertEqual([a.text() for a in more.menu().actions()], ["hide"])

    def test_icons(self):
        from pymol.Qt import QtGui
        from pmg_qt import icons

        white = QtGui.QColor("white")
        for label in icons.LABELS:
            image = icons.menu_icon(label, white).pixmap(16, 16).toImage()
            # drawn in the color asked for, on a clear background
            drawn = [
                image.pixelColor(x, y)
                for x in range(image.width())
                for y in range(image.height())
                if image.pixelColor(x, y).alpha() > 100
            ]
            self.assertTrue(drawn, label)
            self.assertTrue(all(c.lightness() > 240 for c in drawn), label)
        self.assertIsNone(icons.menu_icon("no such command", white))
        # destructive rows get the trash icon, not a red square
        source = open(menu.__file__).read()
        risky = re.findall(r"(?:del|rem)_col \+ '([^']+)'", source)
        self.assertGreater(len(risky), 20)
        for label in risky:
            self.assertIsNotNone(icons.menu_icon(label, white), label)
        # every top-level command of the viewport's menu has an icon
        for entry in menu.main_menu(cmd, (0, 0, 0), (0, 0)):
            if entry[0] == 1:
                self.assertIn(entry[1].replace(r"\933", ""), icons.LABELS)

    def test_panel(self):
        from pmg_qt.object_panel import Checked, ObjectPanel, Unchecked

        cmd.fragment("gly", "g1")
        cmd.fragment("ala", "a2")
        cmd.group("grp", "a2")
        panel = ObjectPanel(cmd)
        panel.refresh()
        self.assertIs(panel.items["a2"].parent(), panel.items["grp"])
        self.assertEqual(
            [panel.items["g1"].text(c) for c in range(1, 6)], list("ASHLC")
        )

        # a click on the name hides the object, the next one shows it again
        panel._clicked(panel.items["g1"], 0)
        self.assertNotIn("g1", cmd.get_names("objects", enabled_only=1))
        self.assertEqual(panel.items["g1"].checkState(0), Unchecked)
        panel._clicked(panel.items["g1"], 0)
        self.assertIn("g1", cmd.get_names("objects", enabled_only=1))
        self.assertEqual(panel.items["g1"].checkState(0), Checked)
        # a double click arrives after a click hid the row: shown again
        panel._clicked(panel.items["g1"], 0)
        panel._double_clicked(panel.items["g1"], 0)
        self.assertIn("g1", cmd.get_names("objects", enabled_only=1))

        panel.filter.setText("a2")
        self.assertTrue(panel.items["g1"].isHidden())
        self.assertFalse(panel.items["grp"].isHidden())  # its member matches

        cmd.delete("g1")
        panel.refresh()
        self.assertNotIn("g1", panel.items)

    def test_states(self):
        from pmg_qt.object_panel import StatesPanel

        panel = StatesPanel(cmd)
        cmd.fragment("ala")
        panel.refresh()
        self.assertTrue(panel.isHidden())  # nothing to step through

        for state in (2, 3, 4):
            cmd.create("ala", "ala", 1, state)
        panel.refresh()
        self.assertFalse(panel.isHidden())
        self.assertEqual(panel.title.text(), "States")
        self.assertEqual(panel.counter.text(), "1/4")
        self.assertEqual(panel.counter.maximum(), 4)

        for button, frame in (("last", 4), ("first", 1), ("forward", 2)):
            panel.buttons[button].click()
            self.assertEqual(cmd.get_frame(), frame)
            self.assertEqual(panel.counter.value(), frame)
        panel.buttons["previous"].click()
        self.assertEqual(cmd.get_frame(), 1)
        panel.counter.setValue(3)  # typed or stepped by the user
        self.assertEqual(cmd.get_state(), 3)
        cmd.frame(2)  # changed by command: the next poll shows it
        panel.refresh()
        self.assertEqual(panel.counter.text(), "2/4")

        play = panel.buttons["play"]
        play.click()
        self.assertTrue(cmd.get_movie_playing())
        self.assertEqual(play.toolTip(), "Stop")
        play.click()
        self.assertFalse(cmd.get_movie_playing())
        self.assertEqual(play.toolTip(), "Play")

        speeds = {a.text(): a for a in panel.speed.menu().actions()}
        speeds["5 fps"].trigger()
        self.assertEqual(cmd.get_setting_float("movie_fps"), 5)
        self.assertEqual(panel.speed.text(), "5 fps")
        cmd.set("movie_fps", 0.3)
        panel.refresh()
        self.assertEqual(panel.speed.text(), "0.3 fps")

        cmd.mset("1 x8")  # a movie: frames, not states
        panel.refresh()
        self.assertEqual(panel.title.text(), "Frames")
        self.assertEqual(panel.counter.maximum(), 8)

        cmd.delete("ala")
        cmd.mclear()
        cmd.mset()
        panel.refresh()
        self.assertTrue(panel.isHidden())

    def test_wizard(self):
        from pmg_qt.object_panel import ObjectPanel

        panel = ObjectPanel(cmd)
        cmd.wizard("measurement")
        panel.wizard.refresh()
        buttons = [
            panel.wizard.buttons.itemAt(i).widget()
            for i in range(panel.wizard.buttons.count())
        ]
        self.assertFalse(panel.wizard.isHidden())
        done = [b for b in buttons if b.text() == "Done"][0]
        self.assertTrue(any(b.menu() for b in buttons))  # mode popups
        done.click()
        self.assertIsNone(cmd.get_wizard())
        panel.wizard.refresh()
        self.assertTrue(panel.wizard.isHidden())
