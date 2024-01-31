'''
Testing: pymol.creating
'''

import pytest

from pymol import cmd
from pymol import test_utils


def test_undo_fragment():
    frag_name = "ala"
    cmd.fragment(frag_name)
    assert frag_name in cmd.get_names()
    cmd.undo()
    assert frag_name not in cmd.get_names()
    cmd.redo()
    assert frag_name in cmd.get_names()


def test_undo_create():
    cmd.fragment("ala")
    obj_name = "foo"
    cmd.create(obj_name, "ala", 1, 1)
    assert obj_name in cmd.get_names()
    cmd.undo()
    assert obj_name not in cmd.get_names()
    cmd.redo()
    assert obj_name in cmd.get_names()


def test_undo_pseudoatom():
    cmd.pseudoatom('m1')
    assert cmd.count_atoms() == 1
    cmd.undo()
    assert cmd.count_atoms() == 0
    cmd.redo()
    assert cmd.count_atoms() == 1


def test_undo_set_raw_alignment():
    cmd.fab('ACDEF', 'm1')
    cmd.fab('CDE', 'm2')
    index_m1 = [('m1', 12), ('m1', 23), ('m1', 35)]
    index_m2 = [('m2',  2), ('m2', 13), ('m2', 25)]
    raw = [list(t) for t in zip(index_m1, index_m2)]
    cmd.set_raw_alignment('aln', raw)
    assert cmd.index('m1 & aln') == index_m1
    assert cmd.index('m2 & aln') == index_m2
    cmd.undo()
    assert 'aln' not in cmd.get_names()
    cmd.redo()
    assert cmd.index('m1 & aln') == index_m1
    assert cmd.index('m2 & aln') == index_m2


def test_undo_isosurface():
    cmd.fragment('gly', 'm1')
    cmd.set('gaussian_b_floor', 30)
    cmd.set('mesh_width', 5)
    cmd.map_new('map')
    cmd.delete('m1')
    surf_name = 'isoSurf'
    cmd.isosurface(surf_name, 'map')
    assert surf_name in cmd.get_names()
    cmd.undo()
    assert surf_name not in cmd.get_names()
    cmd.redo()
    assert surf_name in cmd.get_names()


def test_undo_isodot():
    cmd.fragment('gly', 'm1')
    cmd.set('gaussian_b_floor', 30)
    cmd.set('mesh_width', 5)
    cmd.map_new('map')
    cmd.delete('m1')
    dot_name = 'isoDot'
    cmd.isosurface(dot_name, 'map')
    assert dot_name in cmd.get_names()
    cmd.undo()
    assert dot_name not in cmd.get_names()
    cmd.redo()
    assert dot_name in cmd.get_names()


def test_undo_isomesh():
    cmd.fragment('gly', 'm1')
    cmd.set('gaussian_b_floor', 30)
    cmd.set('mesh_width', 5)
    cmd.map_new('map')
    cmd.delete('m1')

    # make mesh
    cmd.isomesh('mesh', 'map')
    assert 'mesh' in cmd.get_names()
    cmd.undo()
    assert 'mesh' not in cmd.get_names()
    cmd.redo()
    assert 'mesh' in cmd.get_names()


def test_undo_copy():
    cmd.fragment('ala', 'm1')
    cmd.copy('m3', 'm1')
    cmd.undo()
    assert 'm3' not in cmd.get_names()
    cmd.redo()
    assert 'm3' in cmd.get_names()
    cmd.fragment('his', 'm2')
    cmd.copy('m3', 'm2')
    cmd.undo()
    assert cmd.count_atoms('m1') == cmd.count_atoms('m3')
    cmd.redo()
    assert cmd.count_atoms('m2') == cmd.count_atoms('m3')


def get_objs_in_group(group):
    obj_list = []
    cmd.iterate('g1', 'obj_list.append(model)', space=locals())
    return obj_list


def test_undo_group():
    cmd.pseudoatom('m1')
    cmd.pseudoatom('m2')
    cmd.pseudoatom('m3')
    cmd.group('g1', 'm1 m2')
    assert get_objs_in_group('g1') == ['m1', 'm2']
    cmd.undo()
    assert get_objs_in_group('g1') == []
    cmd.redo()
    assert get_objs_in_group('g1') == ['m1', 'm2']


def test_undo_ungroup():
    cmd.pseudoatom('m1')
    cmd.pseudoatom('m2')
    cmd.pseudoatom('m3')
    cmd.group('g1', 'm1 m2')
    cmd.ungroup('m2')
    assert get_objs_in_group('g1') == ['m1']
    cmd.undo()
    assert get_objs_in_group('g1') == ['m1', 'm2']
    cmd.redo()
    assert get_objs_in_group('g1') == ['m1']


def test_undo_isolevel():
    cmd.fragment('gly', 'm1')
    cmd.set('gaussian_b_floor', 30)
    cmd.set('mesh_width', 5)
    cmd.map_new('map')
    cmd.delete('m1')
    # make mesh
    cmd.isodot('dot', 'map')
    cmd.isodot('dot', 'map', source_state=1, state=-2)
    # check mesh presence by color
    meshcolor = 'red'
    cmd.color(meshcolor, 'dot')
    test_utils.ambientOnly(cmd)
    assert test_utils.imageHasColor(cmd, meshcolor)
    cmd.frame(2)
    assert cmd.get_state() == 2
    assert test_utils.imageHasColor(cmd, meshcolor)
    cmd.isolevel('dot', 5)
    assert test_utils.imageHasColor(cmd, meshcolor)
    cmd.isolevel('dot', 10)
    assert not test_utils.imageHasColor(cmd, meshcolor)
    cmd.undo()
    assert test_utils.imageHasColor(cmd, meshcolor)
    cmd.redo()
    assert not test_utils.imageHasColor(cmd, meshcolor)


def test_undo_gradient():
    pass


def test_undo_symexp():
    pass


def test_undo_extract():
    pass
