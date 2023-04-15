'''
Testing: pymol.editing
'''

import pytest

from pymol import cmd


def test_undo_protect():
    cmd.pseudoatom('m1', pos=[0., 0., 0.])
    cmd.pseudoatom('m1', pos=[1., 0., 0.])
    cmd.protect('m1`1')
    cmd.undo()
    cmd.translate([0., 0., 1.])
    assert [0., 0., 1.] == cmd.get_atom_coords('m1`1')
    assert [1., 0., 1.] == cmd.get_atom_coords('m1`2')
    cmd.protect('m1`1')
    cmd.undo()
    cmd.redo()
    cmd.translate([0., 1., 0.])
    assert [0., 0., 1.] == cmd.get_atom_coords('m1`1')
    assert [1., 1., 1.] == cmd.get_atom_coords('m1`2')


def test_undo_deprotect():
    cmd.pseudoatom('m1', pos=[0., 0., 0.])
    cmd.pseudoatom('m1', pos=[1., 0., 0.])
    cmd.protect('m1`1')
    cmd.deprotect()
    cmd.undo()
    cmd.translate([0., 0., 1.])
    assert [0., 0., 0.] == cmd.get_atom_coords('m1`1')
    assert [1., 0., 1.] == cmd.get_atom_coords('m1`2')
    cmd.protect('m1`1')
    cmd.deprotect()
    cmd.undo()
    cmd.redo()
    cmd.translate([0., 1., 0.])
    assert [0., 1., 0.] == cmd.get_atom_coords('m1`1')
    assert [1., 1., 1.] == cmd.get_atom_coords('m1`2')


def assert_array_equal(arr1, arr2, _not=False):
    import numpy

    a1 = numpy.asarray(arr1)
    a2 = numpy.asarray(arr2)

    assert a1.shape == a2.shape
    assert not _not == numpy.allclose(a1, a2, 0, 0)


def assert_array_not_equal(arr1, arr2):
    assert_array_equal(arr1, arr2, True)


def test_undo_update():
    # 3 states
    cmd.fragment('gly', 'm1')
    cmd.create('m1', 'm1', 1, 2)
    cmd.create('m1', 'm1', 1, 3)
    # second object, 90 degree rotates
    cmd.copy('m2', 'm1')
    cmd.rotate('x', 90, '(m2)', state=0)
    # reference coordsets
    cs = cmd.get_coordset
    cs1 = cs('m1', 1)
    cs2 = cs('m2', 1)
    # m2/3 will change (pre-check)
    assert_array_equal(cs2, cs('m2', 3))
    assert_array_not_equal(cs1, cs('m2', 3))
    # update explicit state
    cmd.update('m2', 'm1', 3, 2)
    # m2/3 has changed
    assert_array_equal(cs1, cs('m2', 3))
    assert_array_not_equal(cs2, cs('m2', 3))
    cmd.undo()
    assert_array_not_equal(cs1, cs('m2', 3))
    assert_array_equal(cs2, cs('m2', 3))
    cmd.redo()
    assert_array_equal(cs1, cs('m2', 3))
    assert_array_not_equal(cs2, cs('m2', 3))


def test_undo_fit():
    cmd.fragment("gly", "m1")
    cmd.create("m2", "m1")
    cmd.rotate("y", "90", "m2")
    rms = cmd.rms_cur("m1", "m2")
    assert rms != pytest.approx(0.0)
    cmd.fit("m1", "m2")
    rms = cmd.rms_cur("m1", "m2")
    assert rms == pytest.approx(0.0)
    cmd.undo()
    rms = cmd.rms_cur("m1", "m2")
    assert rms != pytest.approx(0.0)
    cmd.redo()
    rms = cmd.rms_cur("m1", "m2")
    assert rms == pytest.approx(0.0)
