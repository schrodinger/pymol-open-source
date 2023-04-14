'''
Testing: pymol.importing
'''

import pytest

from pymol import cmd
from pymol import test_utils


def _undo_assert(fn, current_value, original_value):
    assert current_value == fn()
    cmd.undo()
    assert original_value == fn()
    cmd.redo()
    assert current_value == fn()


def test_undo_load_raw():
    cmd.load_raw(open(test_utils.datafile(
        "sampletrajectory.pdb")).read(), "pdb")
    _undo_assert(cmd.count_atoms, 115, 0)


@pytest.mark.parametrize("topext,trjext",
                         [('.pdb', '.dcd'),
                          ('.pdb', '.crd'),
                          ('.gro', '.xtc'),
                          ])
def test_undo_load_traj(topext, trjext):
    pdbfile = test_utils.datafile("sampletrajectory" + topext)
    dcdfile = test_utils.datafile("sampletrajectory" + trjext)
    cmd.load(pdbfile)
    cmd.load_traj(dcdfile, state=0)
    _undo_assert(cmd.count_states, 11, 1)
    cmd.delete('*')
    cmd.load(pdbfile)
    cmd.load_traj(dcdfile, state=1, interval=2, max=3)
    _undo_assert(cmd.count_states, 3, 1)
    cmd.delete('*')
    cmd.load(pdbfile)
    cmd.load_traj(dcdfile, state=1, start=3, stop=5)
    _undo_assert(cmd.count_states, 3, 1)
    cmd.delete('*')
    cmd.load(pdbfile)
    cmd.load_traj(dcdfile, state=1, stop=9, average=3)
    _undo_assert(cmd.count_states, 3, 1)
    cmd.delete('*')
