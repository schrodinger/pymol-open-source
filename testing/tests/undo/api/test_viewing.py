'''
Testing: pymol.importing
'''

import pytest

from pymol import cmd
from chempy import cpv


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


def test_undo_origin():
    cmd.pseudoatom('m1')
    cmd.pseudoatom('m2')
    cmd.pseudoatom('m3', pos=[1, 0, 0])
    # by selection
    cmd.origin('m3')
    cmd.rotate('y', 90, 'm1')
    # by position
    cmd.origin(position=[-1, 0, 0])
    cmd.undo()
    cmd.rotate('y', 90, 'm2')
    coords = []
    cmd.iterate_state(1, 'm1 m2', 'coords.append([x,y,z])', space=locals())
    d = cpv.distance(*coords)
    assert d == pytest.approx(0)

    cmd.delete("*")
    cmd.pseudoatom('m1')
    cmd.pseudoatom('m2')
    cmd.pseudoatom('m3', pos=[1, 0, 0])
    # by selection
    cmd.origin('m3')
    cmd.rotate('y', 90, 'm1')
    # by position
    cmd.origin(position=[-1, 0, 0])
    cmd.undo()
    cmd.redo()
    cmd.rotate('y', 90, 'm2')
    coords = []
    cmd.iterate_state(1, 'm1 m2', 'coords.append([x,y,z])', space=locals())
    d = cpv.distance(*coords)
    assert d == pytest.approx(2 * 2**0.5)
