"""
A set interrupt flag (``cmd.interrupt()``) must not make neighbor-based
commands silently return wrong results, and a refused ``quit`` must not leave
the flag set.

https://github.com/schrodinger/pymol-open-source/issues/533
"""

import contextlib

import pytest

from pymol import cmd, _cmd, CmdException

# poly-Ala helix: has backbone H-bonds and distance-based bonds only
SEQ = "AAAAAAAAAAAAAAAA"


@contextlib.contextmanager
def interrupt_flag(_self=cmd):
    """
    Set the interrupt flag like cmd.interrupt() does. Without a GUI idle loop
    nothing clears it, so always reset it to not leak into other tests.
    """
    _self.interrupt()
    try:
        yield
    finally:
        _cmd.interrupt(_self._COb, 0)


def make_helix(name="m", _self=cmd):
    _self.fab(SEQ, name, ss=1)


def count_atoms_or_error(selection):
    """
    Some selection operators may honor the interrupt by failing with an error.
    That's fine, silently returning a wrong count is not.
    """
    try:
        return cmd.count_atoms(selection)
    except CmdException:
        return None


def test_within_with_interrupt_flag():
    make_helix()
    within = "m within 4 of (m and resi 5)"
    around = "(m and resi 5) around 6"
    n_within = cmd.count_atoms(within)
    n_around = cmd.count_atoms(around)
    assert n_within > 0
    assert n_around > 0

    with interrupt_flag():
        assert cmd.count_atoms(within) == n_within
        assert count_atoms_or_error(around) in (n_around, None)


def test_load_bonds_with_interrupt_flag(tmp_path):
    make_helix()
    # no CONECT records for ATOM records, bonds come from distance search
    filename = str(tmp_path / "helix.pdb")
    cmd.save(filename, "m")
    cmd.delete("*")

    cmd.load(filename, "ref")
    n_bonds = len(cmd.get_model("ref").bond)
    assert n_bonds > 0

    with interrupt_flag():
        cmd.load(filename, "m")
    assert len(cmd.get_model("m").bond) == n_bonds


def test_dss_with_interrupt_flag():
    make_helix()
    cmd.dss("m")
    n_helix = cmd.count_atoms("m and ss H")
    assert n_helix > 0

    cmd.alter("m", "ss='L'")
    with interrupt_flag():
        cmd.dss("m")
    assert cmd.count_atoms("m and ss H") == n_helix


def test_polar_contacts_with_interrupt_flag():
    make_helix()
    # average distance, 0.0 if no contacts were found
    avg = cmd.distance("ref", "m", "m", mode=2)
    assert avg > 0

    with interrupt_flag():
        assert cmd.distance("hb", "m", "m", mode=2) == pytest.approx(avg)


def test_refused_quit_does_not_leave_interrupt_flag():
    import pymol2

    p = pymol2.PyMOL(scheme="widget")  # sets no_quit
    p.start()
    try:
        make_helix(_self=p.cmd)
        area = p.cmd.get_area("m")
        assert area > 0

        p.cmd.quit()  # refused

        # dot surface computation aborts if the interrupt flag is set
        assert p.cmd.get_area("m") == pytest.approx(area)
    finally:
        _cmd.interrupt(p.cmd._COb, 0)
        p.stop()
