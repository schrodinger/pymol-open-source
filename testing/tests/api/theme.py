"""
Dark palette has readable text on all groups and apply() switches back.
The saved UI scale reaches Qt unless the user exported their own.
Headless (offscreen QApplication).
"""

import os

from pymol import testing


class TestTheme(testing.PyMOLTestCase):
    def test_palette(self):
        os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
        from pymol.Qt import QtWidgets
        from pmg_qt import theme

        app = QtWidgets.QApplication.instance() or QtWidgets.QApplication([])
        p = theme.dark_palette()
        for group in (theme.Group.Active, theme.Group.Inactive, theme.Group.Disabled):
            bg = p.color(group, theme.Role.Window).lightness()
            fg = p.color(group, theme.Role.WindowText).lightness()
            self.assertGreater(fg - bg, 60)
        theme.apply(app, True)
        self.assertLess(app.palette().window().color().lightness(), 100)
        menu = QtWidgets.QMenu()
        self.assertLess(menu.palette().window().color().lightness(), 100)
        theme.apply(app, False)
        self.assertGreater(app.palette().window().color().lightness(), 150)
        self.assertEqual(app.style().objectName(), "fusion")

    def test_ui_scale(self):
        from unittest import mock
        from pmg_qt import theme

        with mock.patch.object(theme, "ui_scale", return_value=1.5):
            env = {}
            theme.apply_ui_scale(env)
            self.assertEqual(env, {"QT_SCALE_FACTOR": "1.5"})
            env = {"QT_SCALE_FACTOR": "2"}
            theme.apply_ui_scale(env)
            self.assertEqual(env["QT_SCALE_FACTOR"], "2")
        with mock.patch.object(theme, "ui_scale", return_value=1.0):
            env = {}
            theme.apply_ui_scale(env)
            self.assertEqual(env, {})

    def test_prefer_x11_on_wsl(self):
        from pmg_qt import theme

        env = {"WSL_DISTRO_NAME": "Ubuntu", "DISPLAY": ":0"}
        theme.prefer_x11_on_wsl(env)
        self.assertEqual(env["QT_QPA_PLATFORM"], "xcb;wayland")
        env = {"WSL_DISTRO_NAME": "Ubuntu", "DISPLAY": ":0", "QT_QPA_PLATFORM": "w"}
        theme.prefer_x11_on_wsl(env)
        self.assertEqual(env["QT_QPA_PLATFORM"], "w")  # the user's choice wins
        for env in ({"DISPLAY": ":0"}, {"WSL_DISTRO_NAME": "Ubuntu"}):
            theme.prefer_x11_on_wsl(env)
            self.assertNotIn("QT_QPA_PLATFORM", env)

    def test_syntax_colors_readable_on_dark(self):
        from pmg_qt.syntax import readable

        for name in ("blue", "darkBlue", "black", "brown", "darkGreen"):
            self.assertGreaterEqual(readable(name, True).lightness(), 170)
            self.assertLess(readable(name, False).lightness(), 170)

    def test_default_ui_scale(self):
        from unittest import mock
        from pmg_qt import theme

        with mock.patch.object(theme, "_settings") as settings:
            settings.return_value.value.side_effect = lambda key, default: default
            self.assertEqual(theme.ui_scale(), 1.5)
            settings.return_value.value.side_effect = lambda key, default: "bad"
            self.assertEqual(theme.ui_scale(), 1.5)
            settings.return_value.value.side_effect = lambda key, default: "1.2"
            self.assertEqual(theme.ui_scale(), 1.2)

    def test_float_docks_as_dialogs(self):
        os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
        from pymol.Qt import QtCore, QtWidgets
        from pmg_qt import theme

        app = QtWidgets.QApplication.instance() or QtWidgets.QApplication([])
        mask = QtCore.Qt.WindowType.WindowType_Mask
        f = theme._DialogDocks()
        app.installEventFilter(f)
        try:
            win = QtWidgets.QMainWindow()
            dock = QtWidgets.QDockWidget("d", win)
            dock.setFloating(True)
            dock.show()
            self.assertEqual(dock.windowFlags() & mask, QtCore.Qt.WindowType.Dialog)
            self.assertTrue(dock.isVisible() and dock.isFloating())
            win.addDockWidget(QtCore.Qt.DockWidgetArea.TopDockWidgetArea, dock)
            dock.setFloating(False)  # docked again: no longer a window
            self.assertFalse(dock.isWindow())
        finally:
            app.removeEventFilter(f)

        # only where WSLg draws X11 windows
        with_filter = []
        fake = type("A", (), {"platformName": lambda s: "xcb"})()
        fake.installEventFilter = with_filter.append
        theme.float_docks_as_dialogs(fake, {})
        self.assertEqual(with_filter, [])
