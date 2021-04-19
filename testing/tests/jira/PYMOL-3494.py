'''
ray_color_ramps regression in 2.4
'''

from pymol import cmd, testing

class Test3494(testing.PyMOLTestCase):
    def test(self):
        cmd.fragment("gly", "m1")
        cmd.color('blue')
        cmd.color('yellow', 'elem C+O')
        cmd.set("gaussian_resolution", 3.0)
        cmd.set("gaussian_b_floor", 40)
        cmd.map_new("map", "gaussian", 2.0, "m1", 5.0)
        cmd.isosurface("isosurf", "map", 1.0)
        cmd.ramp_new("ramp", "m1", [0, 5], ["atomic", "atomic"])
        cmd.color("ramp", "isosurf")
        cmd.disable("*")
        cmd.enable("isosurf")

        cmd.set_view((\
             0.696023464,    0.704662204,    0.137847140,\
            -0.434485614,    0.260496944,    0.862185299,\
             0.571640551,   -0.659994185,    0.487477839,\
            -0.000000065,   -0.000000011,  -22.690608978,\
             0.078630894,   -0.295513719,    0.061826922,\
            17.590608597,   27.790609360,  -20.000000000 ))

        self.ambientOnly()

        # PyMOL 2.4 has >3000 mismatch pixels
        # The 312 pixels are outline pixels (ffast-math differences?)
        # Delta=1 needed for older builds to pass (ffast-math differences?)

        cmd.set("ray_color_ramps", 0)
        img = self.get_imagearray(width=100, height=100, ray=1)

        # Test is not compatible with different isosurface_algorithm settings
        # self.assertImageEqual("PYMOL-3494-ray_color_ramps-0.png", img, delta=1)

        # 258 colors observed in gradient
        self.assertTrue(self.imageCountColors(img) > 250)

        cmd.set("ray_color_ramps", 1)
        img = self.get_imagearray(width=100, height=100, ray=1)

        # Test is not compatible with different isosurface_algorithm settings
        # self.assertImageEqual("PYMOL-3494-ray_color_ramps-1.png", img, delta=1, count=312)

        # 2 discrete colors + 1 background color
        self.assertEqual(self.imageCountColors(img), 3)
