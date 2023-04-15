'''
Testing: pymol.importing
'''

import pytest

from pymol import cmd

from pymol.constants import ALL_STATES


def _undo_sele(sele_name,
               current_sele_atom_count,
               current_sele_list,
               original_sele_atom_count,
               original_sele_list):
    assert current_sele_atom_count == cmd.count_atoms(sele_name)
    assert current_sele_list == cmd.get_names("selections")
    cmd.undo()
    assert original_sele_list == cmd.get_names("selections")
    cmd.redo()
    assert current_sele_atom_count == cmd.count_atoms(sele_name)
    assert current_sele_list == cmd.get_names("selections")


def test_undo_select():
    cmd.fragment("gly", "m1")
    NC = 2
    # auto_number_selections=0
    cmd.select("elem C")
    _undo_sele("sele", NC, ["sele"], 0, [])

    assert 0 == cmd.get_setting_int("sel_counter")

    cmd.delete("sele")


def test_undo_auto_number_select():
    cmd.fragment("gly", "m1")
    NC = 2
    # auto_number_selections=1
    cmd.set('auto_number_selections', 1)
    cmd.set('sel_counter', 3)

    cmd.select("elem C")
    _undo_sele("sel04", NC, ["sel04"], 0, [])
    assert 4 == cmd.get_setting_int("sel_counter")

    cmd.set('auto_number_selections', 1)  # ??

    cmd.select(None, "elem C")
    _undo_sele("sel05", NC, ["sel04", "sel05"], 0, ["sel04"])
    assert 5 == cmd.get_setting_int("sel_counter")

    cmd.delete("sel*")

    # name=None always numbers the selection
    cmd.set('auto_number_selections', 0)
    cmd.select(None, "elem C")

    _undo_sele("sel06", NC, ["sel06"], 0, [])
    assert 6 == cmd.get_setting_int("sel_counter")

    cmd.delete("sel*")


def _undo_assert_selections(
    target_sele,
    prev_atom_count,
    prev_all_names,
    prev_enabled_names,
    curr_atom_count,
    curr_all_names,
    curr_enabled_names,
    state=ALL_STATES
):
    assert curr_atom_count == cmd.count_atoms(target_sele, state=state)
    assert curr_all_names == cmd.get_names("selections", enabled_only=0)
    assert curr_enabled_names == cmd.get_names("selections", enabled_only=1)
    cmd.undo()
    assert prev_atom_count == cmd.count_atoms(target_sele, state=state)
    assert prev_all_names == cmd.get_names("selections", enabled_only=0)
    assert prev_enabled_names == cmd.get_names("selections", enabled_only=1)
    cmd.redo()
    assert curr_atom_count == cmd.count_atoms(target_sele, state=state)
    assert curr_all_names == cmd.get_names("selections", enabled_only=0)
    assert curr_enabled_names == cmd.get_names("selections", enabled_only=1)


def test_undo_select_merge():
    cmd.fragment("gly", "m1")
    NC = 2

    # merge with non-existing, enable=0
    cmd.select("foo", "elem C", 0, merge=1)
    _undo_assert_selections("?foo", 0, [], [], NC, ["foo"], [])

    # merge, enable=1
    cmd.select("foo", "elem N", 1)
    _undo_assert_selections("?foo", NC, ["foo"], [], 1, ["foo"], ["foo"])

    cmd.select("foo", "elem C", -1, merge=1)
    _undo_assert_selections(
        "?foo", NC - 1, ["foo"], ["foo"], NC + 1, ["foo"], ["foo"])

    cmd.select("foo", "elem O", -1, merge=2)
    assert NC + 2 == cmd.count_atoms("foo")
    assert ["foo"] == cmd.get_names("selections", enabled_only=0)
    assert ["foo"] == cmd.get_names("selections", enabled_only=1)

    # merge, enable=0
    cmd.select("foo", "elem N", 1)
    cmd.select("foo", "elem N", 0)
    _undo_assert_selections("?foo", 1, ["foo"], ["foo"], 1, ["foo"], [])

    cmd.select("foo", "elem C", -1, merge=1)
    _undo_assert_selections("?foo", NC - 1, ["foo"], [], NC + 1, ["foo"], [])

    cmd.select("foo", "elem O", -1, merge=2)
    _undo_assert_selections("?foo", NC + 1, ["foo"], [], NC - 1, ["foo"], [])

    # state
    cmd.delete("sele")
    cmd.create('m1', 'm1 & elem C', 1, 2)
    cmd.select('present')
    _undo_assert_selections("present", 7, ["foo"], [], 7, [
                            "foo", "sele"], ["sele"])

    cmd.delete("sele")

    cmd.select('present', state=1)
    _undo_assert_selections("present", 7, ["foo"], [], 7, [
                            "foo", "sele"], ["sele"])

    cmd.delete("sele")

    cmd.select('present', state=2)
    _undo_assert_selections("present", 2, ["foo"], [], 2, [
                            "foo", "sele"], ["sele"], state=2)

    cmd.delete("sele")
    cmd.select('present', state=3)
    _undo_assert_selections("present", 0, ["foo"], [], 0, [
                            "foo", "sele"], ["sele"], state=3)

    # domain
    cmd.delete("sele")
    cmd.delete("foo")
    cmd.select("foo", "elem C")

    # foo is disabled; must test enabled state
    cmd.select("bar", "name CA+N+O", domain="foo")
    _undo_assert_selections("?bar", 0, ["foo"], ["foo"], 1, [
                            "foo", "bar"], ["bar"])
