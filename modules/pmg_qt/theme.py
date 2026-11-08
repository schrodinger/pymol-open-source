"""
Dark/light palette for the Qt GUI. The 3D viewport and its internal menu are
drawn by PyMOL itself and are dark already; this makes the Qt panels match.
"""

from pymol.Qt import QtCore, QtGui, QtWidgets

Role = QtGui.QPalette.ColorRole
Group = QtGui.QPalette.ColorGroup

_DARK = {
    Role.Window: (53, 53, 53),
    Role.WindowText: (220, 220, 220),
    Role.Base: (35, 35, 35),
    Role.AlternateBase: (45, 45, 45),
    Role.ToolTipBase: (35, 35, 35),
    Role.ToolTipText: (220, 220, 220),
    Role.Text: (220, 220, 220),
    Role.Button: (53, 53, 53),
    Role.ButtonText: (220, 220, 220),
    Role.BrightText: (255, 80, 80),
    Role.Link: (90, 160, 240),
    Role.Highlight: (42, 130, 218),
    Role.HighlightedText: (255, 255, 255),
    Role.PlaceholderText: (140, 140, 140),
}


def dark_palette():
    p = QtGui.QPalette()
    for role, rgb in _DARK.items():
        p.setColor(role, QtGui.QColor(*rgb))
    grey = QtGui.QColor(127, 127, 127)
    for role in (Role.WindowText, Role.Text, Role.ButtonText, Role.HighlightedText):
        p.setColor(Group.Disabled, role, grey)
    p.setColor(Group.Disabled, Role.Highlight, QtGui.QColor(80, 80, 80))
    return p


def is_dark():
    v = QtCore.QSettings("pymolish", "pymolish-studio").value("theme/dark", True)
    return v in (True, "true", 1, "1")


def apply(app, dark):
    # ponytail: Fusion + palette works on every platform; Qt's
    # QStyleHints.setColorScheme is a no-op on Linux (checked with 6.12)
    app.setStyle("Fusion")
    app.setPalette(dark_palette() if dark else app.style().standardPalette())


def set_dark(dark):
    QtCore.QSettings("pymolish", "pymolish-studio").setValue("theme/dark", bool(dark))
    app = QtWidgets.QApplication.instance()
    apply(app, dark)
    for w in (
        app.topLevelWidgets()
    ):  # styled widgets ignore palette changes until repolished
        w.setStyleSheet(w.styleSheet())
