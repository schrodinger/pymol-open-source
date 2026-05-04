from pymol import cmd
from pymol import test_utils


def _setup_cell_scene(color: str) -> None:
    cmd.load(test_utils.datafile("1rx1.pdb"), "m1")
    cmd.hide("everything")
    cmd.show("cell")
    cmd.set("cell_color", color)
    cmd.bg_color("white")
    cmd.set("ray_opaque_background", 1)
    cmd.orient("m1")
    cmd.zoom("m1", buffer=20)
    test_utils.ambientOnly(cmd)


def test_cell_color_ray():
    _setup_cell_scene("red")
    img = test_utils.get_imagearray(cmd, width=200, height=200, ray=1)
    assert test_utils._imageHasColor(cmd, "red", img)
