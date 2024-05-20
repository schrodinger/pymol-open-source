from pymol import cmd
from pymol import test_utils


@test_utils.requires_version("3.0")
def test_bcif():
    cmd.load(test_utils.datafile("115d.bcif.gz"))
    assert cmd.count_atoms() == 407
