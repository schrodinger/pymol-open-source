"""
Object color palettes: the data (size, contrast against the background each
palette is made for, distance between neighbors in the cycle) and the
"palette" command (new objects, recolor, default, reinitialize, sessions).
"""

import math

import pymol

from pymol import cmd, palettes, testing

BACKGROUND = {"dark": (0, 0, 0), "white": (1, 1, 1)}


def rgb(color):
    return [int(color[i : i + 2], 16) / 255 for i in (0, 2, 4)]


def luminance(color):  # WCAG 2 relative luminance
    r, g, b = (
        v / 12.92 if v <= 0.03928 else ((v + 0.055) / 1.055) ** 2.4 for v in color
    )
    return 0.2126 * r + 0.7152 * g + 0.0722 * b


def contrast(a, b):
    low, high = sorted((luminance(a), luminance(b)))
    return (high + 0.05) / (low + 0.05)


class TestPalettes(testing.PyMOLTestCase):
    def tearDown(self):
        cmd.palette("default")
        super().tearDown()

    def carbon_colors(self):
        """Color index of the first carbon of each object."""
        colors = []
        for name in cmd.get_names():
            cmd.iterate(
                "first (%s & elem C)" % name, "colors.append(color)", space=locals()
            )
        return colors

    def test_data(self):
        self.assertEqual(
            list(palettes.PALETTES), ["vivid", "deep", "pastel", "soft", "colorblind"]
        )
        for name, (background, colors) in palettes.PALETTES.items():
            with self.subTest(name):
                self.assertTrue(1 <= len(colors) <= 30)
                self.assertEqual(len(set(colors)), len(colors))
                values = [rgb(c) for c in colors]
                for a, b in zip(values, values[1:] + values[:1]):
                    self.assertGreater(math.dist(a, b), 0.14, (name, a, b))
                # colorblind: published sets for any background, not judged
                for value in values if background in BACKGROUND else ():
                    self.assertGreaterEqual(contrast(value, BACKGROUND[background]), 3)

    def test_new_objects(self):
        cmd.palette("vivid")
        for name in ("ala", "gly", "ser"):
            cmd.fragment(name)
        vivid = [cmd.get_color_index("vivid_%02d" % i) for i in (1, 2, 3)]
        self.assertEqual(self.carbon_colors(), vivid)
        self.assertEqual(cmd.get_object_color_index("gly"), vivid[1])
        self.assertArrayEqual(cmd.get_color_tuple("vivid_01"), rgb("00e676"), 1e-3)
        # only carbons
        self.assertEqual(cmd.count_atoms("ala & not elem C & color vivid_01"), 0)

        cmd.delete("*")
        cmd.palette("default")
        for name in ("ala", "gly", "ser"):
            cmd.fragment(name)
        self.assertEqual(self.carbon_colors(), [26, 5, 154])

    def test_recolor(self):
        for name in ("ala", "gly"):
            cmd.fragment(name)
        cmd.palette("deep")
        self.assertEqual(self.carbon_colors(), [26, 5])  # untouched
        cmd.palette("deep", recolor=1)
        deep = [cmd.get_color_index("deep_%02d" % i) for i in (1, 2, 3)]
        self.assertEqual(self.carbon_colors(), deep[:2])
        self.assertEqual(cmd.count_atoms("ala & elem O & color deep_01"), 0)
        cmd.fragment("ser")  # continues the cycle
        self.assertEqual(self.carbon_colors(), deep)

        cmd.palette("soft")
        cmd.util.color_objs()
        soft = cmd.get_color_index("soft_02")
        self.assertEqual(
            cmd.count_atoms("gly & color %d" % soft), cmd.count_atoms("gly")
        )

    def test_wraps_and_unknown(self):
        cmd.palette("deep")
        n = len(palettes.PALETTES["deep"][1])
        cmd.set("auto_color_next", n - 1)
        cmd.fragment("ala")
        cmd.fragment("gly")
        self.assertEqual(
            self.carbon_colors(),
            [cmd.get_color_index("deep_%02d" % i) for i in (n, 1)],
        )
        self.assertRaises(pymol.CmdException, cmd.palette, "nosuch")
        self.assertEqual(palettes.active, "deep")

    def test_survives_reinitialize_and_session(self):
        cmd.palette("soft")
        cmd.reinitialize()
        cmd.fragment("ala")
        self.assertEqual(self.carbon_colors(), [cmd.get_color_index("soft_01")])
        self.assertArrayEqual(cmd.get_color_tuple("soft_01"), rgb("6192c3"), 1e-3)

        # a session saved without the palette colors drops them from the
        # color list; loading it must not leave the table pointing at nothing
        cmd.palette("default")
        cmd.reinitialize()
        cmd.fragment("gly")
        with testing.mktemp(".pse") as filename:
            cmd.save(filename)
            cmd.palette("soft")
            cmd.load(filename)
        cmd.fragment("ala")
        # the session brings its own auto_color_next: 1
        self.assertEqual(self.carbon_colors(), [26, cmd.get_color_index("soft_02")])
        self.assertArrayEqual(cmd.get_color_tuple("soft_01"), rgb("6192c3"), 1e-3)
