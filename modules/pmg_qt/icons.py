"""
Thin-line icons for menu entries and toolbar buttons. The drawings are Lucide icons (ISC license,
data/pmg_qt/icons/lucide/LICENSE), tinted with the color of the text they sit
next to.
"""

import functools
import os

from pymol.Qt import QtCore, QtGui, QtWidgets

# menu label (pymol.menu, color markup removed) -> icon file
LABELS = {
    "new": "circle-plus",
    "pseudoatom": "circle-dot",
    "zoom": "zoom-in",
    "zoom (vis)": "zoom-in",
    "orient": "compass",
    "orient (vis)": "compass",
    "center": "crosshair",
    "center (vis)": "crosshair",
    "origin": "locate-fixed",
    "reset": "rotate-ccw",
    "movie": "clapperboard",
    "scene": "image",
    "enable": "eye",
    "show": "eye",
    "disable": "eye-off",
    "hide": "eye-off",
    "(all)": "layers",
    "(visible)": "scan-eye",
    "ray": "aperture",
    "delete": "trash-2",
    "delete all": "trash-2",
    "delete selection": "trash-2",
    "remove": "trash-2",
    "remove atoms": "trash-2",
    "reinitialize": "refresh-cw",
    "update": "refresh-cw",
    "quit": "log-out",
    "drag": "move",
    "drag object matrix": "move",
    "drag object coords": "move",
    "atom": "atom",
    "residue": "hexagon",
    "chain": "link",
    "segment": "spline",
    "object": "box",
    "molecule": "share-2",
    "fragment": "puzzle",
    "fragment+joint(s)": "puzzle",
    "actions": "zap",
    "color": "palette",
    "preset": "sliders-horizontal",
    "label": "tag",
    "ss": "dna",
    "clean": "sparkles",
    "select": "mouse-pointer-click",
    "create object": "copy-plus",
    "copy to object": "copy",
    "duplicate": "copy",
    "extract object": "scissors",
    "rename": "pencil",
    "rename selection": "pencil",
    "find": "search",
    "compute": "calculator",
    "clear": "eraser",
    "play": "play",
    "stop": "square",
    "rewind": "skip-back",
    "next": "skip-forward",
    # states panel
    "first": "skip-back",
    "previous": "rewind",
    "forward": "fast-forward",
    "last": "skip-forward",
    # action toolbar
    "undo": "undo-2",
    "redo": "redo-2",
    "rock": "rotate-3d",
    "presets": "sliders-horizontal",
    "builder": "hammer",
    "scenes": "image",
    "draw/ray": "aperture",
    "ui scale": "scaling",
    "more actions": "ellipsis",
}


@functools.lru_cache(maxsize=None)
def _pixmap(name, color, ratio):
    path = os.path.expandvars("$PYMOL_DATA/pmg_qt/icons/lucide/%s.svg" % name)
    with open(path) as handle:
        svg = handle.read()
    # the files draw 2 px lines on a 24 px grid in "currentColor"
    svg = svg.replace("currentColor", color).replace(
        'stroke-width="2"', 'stroke-width="1.5"'
    )
    data = QtCore.QBuffer()
    data.setData(svg.encode())
    reader = QtGui.QImageReader(data, b"svg")
    side = round(16 * ratio)
    reader.setScaledSize(QtCore.QSize(side, side))
    pixmap = QtGui.QPixmap.fromImage(reader.read())
    pixmap.setDevicePixelRatio(ratio)
    return pixmap


def menu_icon(label, color):
    """The icon for a menu label, drawn in `color` (a QColor); None if the
    label has none. "delete object" falls back to the icon of "delete"."""
    name = LABELS.get(label) or LABELS.get(label.split(" ")[0])
    if name:
        ratio = QtWidgets.QApplication.instance().devicePixelRatio()
        return QtGui.QIcon(_pixmap(name, color.name(), ratio))
