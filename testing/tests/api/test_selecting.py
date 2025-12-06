from pymol import cmd


def test_select_operators():
    cmd.reinitialize()
    cmd.pseudoatom(pos=[1, 2, 3], b=5)

    assert cmd.count_atoms("b > 4") == 1
    assert cmd.count_atoms("b < 6") == 1
    assert cmd.count_atoms("b = 5") == 1
    assert cmd.count_atoms("b >= 5") == 1
    assert cmd.count_atoms("b <= 5") == 1

    assert cmd.count_atoms("b > 4 & x == 1") == 1
    assert cmd.count_atoms("b > 5 & x == 1") == 0
    assert cmd.count_atoms("b < 6 & y <= 3") == 1
    assert cmd.count_atoms("b = 5 & z >= 2") == 1
