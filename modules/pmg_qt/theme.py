"""
Dark/light palette for the Qt GUI. The 3D viewport and its internal menu are
drawn by PyMOL itself and are dark already; this makes the Qt panels match.
"""

import os

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


def _settings():
    return QtCore.QSettings("pymolish", "pymolish-studio")


def is_dark():
    return _settings().value("theme/dark", True) in (True, "true", 1, "1")


# Multiplies whatever scale the OS reports. Needed where the OS reports none
# (WSLg, bare X11: a 2560x1600 laptop panel shows up as 96 dpi, scale 1).
SCALES = (1, 1.2, 1.5, 1.8, 2)
DEFAULT_SCALE = 1.5


def ui_scale():
    try:
        return float(_settings().value("ui/scale", DEFAULT_SCALE))
    except (TypeError, ValueError):
        return DEFAULT_SCALE


def set_ui_scale(scale):
    _settings().setValue("ui/scale", float(scale))


def object_palette():
    """The palette picked in the toolbar: see pymol.palettes."""
    return str(_settings().value("colors/palette", "default"))


def set_object_palette(name):
    _settings().setValue("colors/palette", name)


def apply_ui_scale(environ=os.environ):
    """Call before QApplication exists: Qt reads the factor once, at startup.
    A factor the user exported wins over the saved one."""
    scale = ui_scale()
    if scale != 1 and not any(
        v in environ for v in ("QT_SCALE_FACTOR", "QT_SCREEN_SCALE_FACTORS")
    ):
        environ["QT_SCALE_FACTOR"] = "%g" % scale


def prefer_x11_on_wsl(environ=os.environ):
    """Call before QApplication exists. WSLg's Wayland compositor does not
    keep popups on the screen, so a menu at the window edge runs off it and
    its submenus cannot be reached. Through X11 Qt places menus itself: they
    stay inside the screen and submenus open to the left when they must."""
    if "WSL_DISTRO_NAME" in environ and environ.get("DISPLAY"):
        environ.setdefault("QT_QPA_PLATFORM", "xcb;wayland")


class _RoundMenus(QtCore.QObject):
    def eventFilter(self, watched, event):
        if (
            event.type() == QtCore.QEvent.Type.Polish
            and isinstance(watched, QtWidgets.QMenu)
            and not watched.testAttribute(_TRANSLUCENT)
        ):
            watched.setAttribute(_TRANSLUCENT)
            watched.setWindowFlags(
                watched.windowFlags()
                | QtCore.Qt.WindowType.FramelessWindowHint
                | QtCore.Qt.WindowType.NoDropShadowWindowHint
            )
        return False


_TRANSLUCENT = QtCore.Qt.WidgetAttribute.WA_TranslucentBackground


def round_menus(app):
    """The stylesheet rounds the corners of every menu; this makes the
    window behind those corners see-through, for each menu before it is
    first shown, wherever it was created."""
    app.installEventFilter(_RoundMenus(app))


def apply(app, dark):
    # Fusion + palette darkens every Qt-drawn menu and widget. The color
    # scheme hint darkens what the OS draws (macOS menu bar, Windows title
    # bars and system menus); it is a no-op on Linux (checked with 6.12).
    scheme = QtCore.Qt.ColorScheme
    app.styleHints().setColorScheme(scheme.Dark if dark else scheme.Light)
    app.setStyle("Fusion")
    app.setPalette(dark_palette() if dark else app.style().standardPalette())


def set_dark(dark):
    _settings().setValue("theme/dark", bool(dark))
    app = QtWidgets.QApplication.instance()
    apply(app, dark)
    for w in (
        app.topLevelWidgets()
    ):  # styled widgets ignore palette changes until repolished
        w.setStyleSheet(w.styleSheet())
