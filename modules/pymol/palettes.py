"""
Object color palettes: the order in which new objects get their carbon
color. "palette" swaps the table the loaders walk (AutoColor in
layer1/Color.cpp); the colors are named colors, <palette>_01 and up.
"""

import pymol

cmd = __import__("sys").modules["pymol.cmd"]
from .cmd import _cmd  # noqa: E402

# name -> (background the colors are picked for, hex colors in cycle order).
# testing/tests/api/palettes.py checks contrast and neighbor distance.
PALETTES = {
    # Material Design accent colors, near-duplicates removed
    "vivid": (
        "dark",
        "00e676 40c4ff ff4081 ffd740 b388ff ff6e40 18ffff 69f0ae e040fb 448aff"
        " b2ff59 ff80ab ffff00 7c4dff ff9100 00b0ff ea80fc 76ff03 ff1744 ffc400"
        " aeea00",
    ),
    # Material Design 800 and 900 tones, near-duplicates removed
    "deep": (
        "white",
        "1565c0 c62828 2e7d32 6a1b9a e65100 00838f ad1457 283593 00695c 827717"
        " 4e342e d81b60 b8860b 37474f 8d6e63",
    ),
    # PASTELS of the pymolish plugin
    "pastel": (
        "dark",
        "ffcccc ccffcc ccccff ffffcc ffccff ccffff ffe6cc e6ccff ccffe6 ffcce6"
        " e6e6cc cce6ff f2ccd9 d9f2d9 e6d9ff ffe6e6 ccf2f2 f2f2cc e6ccf2 d9ffd9"
        " ffd9cc ccd9f2 f2ffcc d9ccff fff2d9",
    ),
    # 12 hues at half saturation and equal luminance, then 10 darker ones
    "soft": (
        "white",
        "6192c3 da6d6d 4e9d76 cd67cd 73994d 7f7fff b28559 4d9a9a d56aa0 4f9f4f"
        " ad73e7 929249 4b68bd a15840 307866 ac4592 427930 7555d4 7d6a32 367488"
        " b34762 9a49b6",
    ),
    # Okabe-Ito, then Paul Tol's bright and muted sets, values as published.
    # Left out: black, the greys, and Tol colors close to an earlier entry.
    # Only the first seven are designed as one set.
    "colorblind": (
        "any",
        "e69f00 56b4e9 009e73 f0e442 0072b2 d55e00 cc79a7 4477aa 228833 ccbb44"
        " ee6677 aa3377 332288 88ccee 44aa99 999933 ddcc77 882255",
    ),
}
PALETTES = {name: (bg, colors.split()) for name, (bg, colors) in PALETTES.items()}

active = "default"  # process-wide, like the table in Color.cpp


def apply(_self=cmd):
    """Register the colors of the active palette and hand them to the
    loaders. Also runs after reinitialize and a session load: both rebuild
    the color list, so the indices change."""
    names = []
    for i, color in enumerate(PALETTES.get(active, ("", ()))[1], 1):
        names.append("%s_%02d" % (active, i))
        _self.set_color(names[-1], [int(color[j : j + 2], 16) for j in (0, 2, 4)])
    indices = [_self.get_color_index(name) for name in names]
    with _self.lockcm:
        _cmd.set_auto_colors(_self._COb, indices)  # empty: the PyMOL order


def session_restore(session, _self=cmd):
    apply(_self)
    return True


def palette(name="", recolor=0, quiet=1, *, _self=cmd):
    """
    DESCRIPTION

        "palette" sets the colors that new objects cycle through. Without
        arguments it lists the palettes and marks the active one.

    USAGE

        palette [ name [, recolor ]]

    ARGUMENTS

        name = default, vivid, deep, pastel, soft or colorblind.
        "default" is the PyMOL order (green, cyan, light magenta, ...).

        recolor = 0/1: also color the carbons of the loaded objects, from the
        first palette color on {default: 0}

    EXAMPLE

        palette vivid
        palette soft, 1

    NOTES

        vivid and pastel are made for a dark background, deep and soft for a
        white one. The colors are named <palette>_01, <palette>_02, ... and can
        be used like any color. A "palette" line in pymolrc applies at every
        start.

    SEE ALSO

        color, set_color, util.color_objs
    """
    global active
    if name:
        if name != "default" and name not in PALETTES:
            raise pymol.CmdException(
                "unknown palette '%s', choose from: default %s"
                % (name, " ".join(PALETTES))
            )
        active = name
        apply(_self)
        _self.set("auto_color_next", 0)
        if not int(quiet):
            print(" palette: new objects use '%s'" % name)
    elif not int(recolor):
        for name in ("default", *PALETTES):
            print(" %s %s" % ("*" if name == active else " ", name))
    if int(recolor):
        _self.set("auto_color_next", 0)
        for name in _self.get_names("public_nongroup_objects"):
            if _self.get_type(name) == "object:molecule":
                color = str(_self.get_color_index("auto"))
                _self.set_object_color(name, color)
                _self.color(color, "?%s & elem C" % name)
