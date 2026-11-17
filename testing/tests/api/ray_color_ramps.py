'''
ray_color_ramps with molecular surfaces and isomeshes

Per-vertex ramp colors (RC) must be ignored when ray_color_ramps is off.
Stale RC data after std::vector::clear() made the ray tracer evaluate the
ramp per pixel anyway.
https://github.com/schrodinger/pymol-open-source/issues/532
'''

from pymol import cmd, testing


@testing.requires_version('3.2')
class TestRayColorRamps(testing.PyMOLTestCase):

    def _atomic_ramp(self):
        cmd.fragment("gly", "m1")
        cmd.color("blue")
        cmd.color("yellow", "elem C+O")
        cmd.ramp_new("ramp", "m1", [0, 5], ["atomic", "atomic"])

    def _assert_ray_color_ramps(self, max_colors_on, min_colors_off):
        self.ambientOnly()

        # ray_color_ramps=1: ramp evaluated per pixel -> nearest atom color,
        # a few discrete colors
        cmd.set("ray_color_ramps", 1)
        img = self.get_imagearray(width=100, height=100, ray=1)
        self.assertLess(self.imageCountColors(img), max_colors_on)

        # ray_color_ramps=0: per-vertex colors interpolated -> gradient
        # (with the bug: same few discrete colors as ray_color_ramps=1)
        cmd.set("ray_color_ramps", 0)
        img = self.get_imagearray(width=100, height=100, ray=1)
        self.assertGreater(self.imageCountColors(img), min_colors_off)

    def testSurface(self):
        self._atomic_ramp()
        cmd.set("surface_color", "ramp", "m1")
        cmd.show_as("surface", "m1")
        cmd.orient("m1")

        # observed: 5 colors (on), ~260 colors (off)
        self._assert_ray_color_ramps(10, 100)

    def testIsomesh(self):
        self._atomic_ramp()
        cmd.set("gaussian_resolution", 3.0)
        cmd.set("gaussian_b_floor", 40)
        cmd.map_new("map", "gaussian", 2.0, "m1", 5.0)
        cmd.isomesh("mesh", "map", 1.0)
        cmd.color("ramp", "mesh")
        cmd.disable("*")
        cmd.enable("mesh")
        cmd.orient("m1")

        # observed: 3 colors (on), ~58 colors (off)
        self._assert_ray_color_ramps(10, 20)
