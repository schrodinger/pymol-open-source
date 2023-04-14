from pymol import cmd
import pytest


def test_look_at():
    ori_view = cmd.get_view()
    cmd.pseudoatom("M1", pos=[10, 0, 0])
    cmd.look_at("M1")
    new_view = cmd.get_view()
    assert ori_view != new_view

    ref_new_view = (0.980580688,    0.000000000,   -0.196116135,
                    0.000000000,    1.000000000,    0.000000000,
                    0.196116135,    0.000000000,    0.980580688,
                    -9.805807114,    0.000000000,  -49.029033661,
                    0.000000000,    0.000000000,    0.000000000,
                    40.000000000,  100.000000000,  -20.000000000)
    assert ref_new_view == pytest.approx(new_view)
