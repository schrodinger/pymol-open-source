"""
PYMOL-3487 Partial load of ramps
"""

from os.path import join
from pymol import cmd, testing

data_dir = "PYMOL-3487-data"

@testing.requires_version("2.5")
class Test3487(testing.PyMOLTestCase):

    def _assert_colors(self, colors=['red', 'yellow']):
        self.ambientOnly()
        cmd.orient()
        img = self.get_imagearray(width=100, height=100, ray=1)
        for color in colors:
            self.assertImageHasColor(color, img)
        return img

    def test_order1(self):
        cmd.load(join(data_dir, "session1.pze"), partial=1)
        ramp1_index = cmd.get_color_index("ramp1")
        cmd.load(join(data_dir, "session2.pze"), partial=1)
        self.assertEqual(cmd.get_color_index("ramp1"), ramp1_index)
        cmd.load(join(data_dir, "session3.pze"), partial=1)
        self.assertEqual(cmd.get_color_index("ramp1"), ramp1_index)
        self.assertEqual(cmd.get_session()["color_ext"], [
            ["ramp1", 1],
            ["ramp2", 1],
        ])
        self._assert_colors()
        self.assertEqual(ramp1_index, -10)
        ramp2_index = cmd.get_color_index("ramp2")
        self.assertEqual(ramp2_index, -11)

    def test_order2(self):
        cmd.load(join(data_dir, "session2.pze"), partial=1)
        ramp1_index = cmd.get_color_index("ramp1")
        cmd.load(join(data_dir, "session3.pze"), partial=1)
        self.assertEqual(cmd.get_color_index("ramp1"), ramp1_index)
        ramp2_index = cmd.get_color_index("ramp2")
        cmd.load(join(data_dir, "session1.pze"), partial=1)
        self.assertEqual(cmd.get_color_index("ramp1"), ramp1_index)
        self.assertEqual(cmd.get_session()["color_ext"], [
            ["ramp1", 1],
            ["ramp2", 1],
        ])
        self._assert_colors()
        self.assertEqual(ramp1_index, -10)
        self.assertEqual(ramp2_index, -11)

    def test_order3(self):
        cmd.load(join(data_dir, "session3.pze"), partial=1)
        ramp2_index = cmd.get_color_index("ramp2")
        cmd.load(join(data_dir, "session1.pze"), partial=1)
        self.assertEqual(cmd.get_color_index("ramp2"), ramp2_index)
        ramp1_index = cmd.get_color_index("ramp1")
        cmd.load(join(data_dir, "session2.pze"), partial=1)
        self.assertEqual(cmd.get_color_index("ramp1"), ramp1_index)
        self.assertEqual(cmd.get_session()["color_ext"], [
            ["ramp1", 1],
            ["ramp2", 1],
        ])
        self._assert_colors()
        self.assertEqual(ramp1_index, -10)
        self.assertEqual(ramp2_index, -11)

    def test_overwrite_existing(self):
        cmd.pseudoatom('existing1', pos=(1.4, 3.8, 0.0))
        cmd.pseudoatom('existing2', pos=(-1.0, 0.5, -0.1))
        cmd.ramp_new("ramp1", "existing1", [0, 10], ["green", "green"])
        ramp1_index = cmd.get_color_index("ramp1")
        cmd.ramp_new("ramp3", "existing2", [0, 10], ["magenta", "magenta"])
        ramp3_index = cmd.get_color_index("ramp3")
        self.assertEqual(cmd.get_color_index("ramp1"), ramp1_index)
        self.assertNotEqual(ramp3_index, ramp1_index)
        self.assertEqual(cmd.get_session()["color_ext"], [
            ["ramp1", 1],
            ["ramp3", 1],
        ])
        self.assertEqual(ramp1_index, -10)
        self.assertEqual(ramp3_index, -11)

        cmd.show_as('surface', 'existing1 existing2')
        cmd.set('surface_color', 'ramp1', 'existing2')
        cmd.set('surface_color', 'ramp3', 'existing1')
        cmd.disable('ramp*')
        self._assert_colors(['green', 'magenta'])

        cmd.load(join(data_dir, "session1.pze"), partial=1)
        cmd.load(join(data_dir, "session2.pze"), partial=1)
        cmd.load(join(data_dir, "session3.pze"), partial=1)
        self.assertEqual(cmd.get_color_index("ramp1"), ramp1_index)
        self.assertEqual(cmd.get_color_index("ramp3"), ramp3_index)
        self.assertEqual(cmd.get_session()["color_ext"], [
            ["ramp1", 1],
            ["ramp3", 1],
            ["ramp2", 1],
        ])
        ramp2_index = cmd.get_color_index("ramp2")
        self.assertEqual(ramp1_index, -10)
        self.assertEqual(ramp3_index, -11)
        self.assertEqual(ramp2_index, -12)

        cmd.recolor()
        img = self._assert_colors(['red', 'yellow', 'magenta'])
        self.assertImageHasNotColor('green', img)

        self.assertEqual(cmd.get_setting_int("surface_color", "mol1"), ramp1_index)
        self.assertEqual(cmd.get_setting_int("surface_color", "mol2"), ramp2_index)
        self.assertEqual(cmd.get_setting_int("surface_color", "existing1"), ramp3_index)
        self.assertEqual(cmd.get_setting_int("surface_color", "existing2"), ramp1_index)
