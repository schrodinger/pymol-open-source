'''
Testing: pymol.importing
'''

import pytest

from pymol import cmd


def test_undo_clean():
    cmd.fragment('his')
    xyz1 = cmd.get_model('his').get_coord_list()
    cmd.clean('his')
    xyz2 = cmd.get_model('his').get_coord_list()
    assert len(xyz1) == len(xyz2)
    # coords have changed
    assert xyz1 != xyz2
    cmd.undo()
    xyz2 = cmd.get_model('his').get_coord_list()
    assert xyz1 == xyz2
    cmd.redo()
    xyz2 = cmd.get_model('his').get_coord_list()
    assert xyz1 != xyz2
