"""
Dark palette has readable text on all groups and apply() switches back.
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
        theme.apply(app, False)
        self.assertGreater(app.palette().window().color().lightness(), 150)
        self.assertEqual(app.style().objectName(), "fusion")
