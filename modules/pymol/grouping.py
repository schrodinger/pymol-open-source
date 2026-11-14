"""
Group commands: step through the members of a group, show or hide them,
copy or sort a group, split objects into chains. Ported from the pymolish
plugin (MIT, same author); new code, MIT.
"""

import pymol

cmd = __import__("sys").modules["pymol.cmd"]
from .cmd import _cmd  # noqa: E402


def members(group, _self=cmd):
    """Names of the direct members of a group, in panel order."""
    if group not in _self.get_names("group_objects"):
        raise pymol.CmdException("'%s' is not a group" % group)
    with _self.lockcm:
        return _cmd.get_group_members(_self._COb, group).split()


def _active_group(_self):
    """The group to step through when none is named: the only group, or the
    only one with an enabled member."""
    groups = _self.get_names("group_objects")
    if len(groups) > 1:
        enabled = set(_self.get_names("objects", enabled_only=1))
        groups = [g for g in groups if enabled.intersection(members(g, _self))]
    if len(groups) != 1:
        raise pymol.CmdException("name the group to step through")
    return groups[0]


def _step(group, direction, quiet, _self):
    group = group or _active_group(_self)
    names = members(group, _self)
    if not names:
        raise pymol.CmdException("group '%s' is empty" % group)
    enabled = set(_self.get_names("objects", enabled_only=1))
    shown = [i for i, name in enumerate(names) if name in enabled]
    # nothing shown: forward starts at the first member, backward at the last
    start = shown[0] if shown else (-1 if direction > 0 else 0)
    name = names[(start + direction) % len(names)]
    for i in shown:
        _self.disable(names[i])
    _self.enable(name, parents=1)
    if not int(quiet):
        print(" %s: %s" % (group, name))
    return name


def group_next(group="", quiet=1, *, _self=cmd):
    """
    DESCRIPTION

        "group_next" shows the next member of a group and hides the others,
        so that a group can be viewed one object at a time. After the last
        member comes the first.

    USAGE

        group_next [ group ]

    ARGUMENTS

        group = string: name of the group {default: the only group, or the
        only group with a visible member}

    NOTES

        The key is CTRL-N; CTRL-P is "group_previous".

    SEE ALSO

        group_previous, group_enable, group_disable
    """
    return _step(group, 1, quiet, _self)


def group_previous(group="", quiet=1, *, _self=cmd):
    """
    DESCRIPTION

        "group_previous" shows the previous member of a group and hides the
        others. Before the first member comes the last.

    USAGE

        group_previous [ group ]

    SEE ALSO

        group_next, group_enable, group_disable
    """
    return _step(group, -1, quiet, _self)


def _set_enabled(group, each, on, _self):
    names = members(group, _self)[:: max(1, int(each))]
    for name in names:
        if on:
            _self.enable(name, parents=1)
        else:
            _self.disable(name)
    return names


def group_enable(group, each=1, *, _self=cmd):
    """
    DESCRIPTION

        "group_enable" shows the members of a group. "enable" on a group
        only shows the group itself.

    USAGE

        group_enable group [, each ]

    ARGUMENTS

        group = string: name of the group

        each = integer: every n-th member only {default: 1, all}

    EXAMPLE

        group_enable poses
        group_enable poses, 10

    SEE ALSO

        group_disable, group_next, enable
    """
    return _set_enabled(group, each, True, _self)


def group_disable(group, each=1, *, _self=cmd):
    """
    DESCRIPTION

        "group_disable" hides the members of a group; the group stays
        enabled.

    USAGE

        group_disable group [, each ]

    SEE ALSO

        group_enable, group_next, disable
    """
    return _set_enabled(group, each, False, _self)


def _copy_rows(group, new_group, suffix, _self):
    """(copy, source, group of the copy) per member, nested groups included."""
    rows = []
    for name in members(group, _self):
        rows.append((name + suffix, name, new_group))
        if _self.get_type(name) == "object:group":
            rows += _copy_rows(name, name + suffix, suffix, _self)
    return rows


def copy_group(group, new_group="", suffix="_cpy", quiet=1, *, _self=cmd):
    """
    DESCRIPTION

        "copy_group" copies a group with all its members, groups inside it
        included. Hidden members stay hidden in the copy.

    USAGE

        copy_group group [, new_group [, suffix ]]

    ARGUMENTS

        group = string: name of the group

        new_group = string: name of the copy {default: group + suffix}

        suffix = string: added to the name of every member {default: _cpy}

    EXAMPLE

        copy_group poses
        copy_group poses, backup, _bak

    SEE ALSO

        copy, group, sort_group
    """
    new_group = new_group or group + suffix
    rows = _copy_rows(group, new_group, suffix, _self)
    taken = set(_self.get_names("all"))
    clash = [n for n in [new_group] + [row[0] for row in rows] if n in taken]
    if clash:
        raise pymol.CmdException("already exists: " + " ".join(clash))
    enabled = set(_self.get_names("objects", enabled_only=1))
    _self.group(new_group)
    for new, source, parent in rows:
        if _self.get_type(source) == "object:group":
            _self.group(new)
        else:
            _self.copy(new, source, zoom=0)
        if source not in enabled:
            _self.disable(new)
        _self.group(parent, new, "add")
    if not int(quiet):
        print(" copy_group: %d objects in '%s'" % (len(rows), new_group))
    return [row[0] for row in rows]


def sort_group(group, reverse=0, *, _self=cmd):
    """
    DESCRIPTION

        "sort_group" sorts the members of a group by name, ignoring case.

    USAGE

        sort_group group [, reverse ]

    ARGUMENTS

        group = string: name of the group

        reverse = 0/1: Z to A {default: 0}

    SEE ALSO

        order, group
    """
    names = sorted(members(group, _self), key=str.lower, reverse=bool(int(reverse)))
    if names:
        _self.order(" ".join(names))
    return names


def is_sorted(group, _self=cmd):
    """True if the members are in A to Z order; the menu row then sorts Z
    to A."""
    names = members(group, _self)
    return len(names) > 1 and names == sorted(names, key=str.lower)


def extract_chain(selection, chains="", group="", quiet=1, *, _self=cmd):
    """
    DESCRIPTION

        "extract_chain" makes one new object per chain, named
        <object>_<chain>. The source objects are not changed.

    USAGE

        extract_chain selection [, chains [, group ]]

    ARGUMENTS

        selection = string: atoms to take; a group name takes all its
        members

        chains = string: chain identifiers joined with "+" {default: every
        chain}

        group = string: group to put the new objects in {default: none}

    EXAMPLE

        extract_chain 1abc
        extract_chain poses, A+B, poses_AB

    SEE ALSO

        split_chains, create, extract
    """
    wanted = chains.replace("+", " ").split()
    rows = []
    for obj in _self.get_object_list(selection):
        sele = "(%s) and %%%s" % (selection, obj)
        for chain in wanted or _self.get_chains(sele):
            chain_sele = '%s and chain "%s"' % (sele, chain)
            if _self.count_atoms(chain_sele):
                name = _self.get_legal_name("%s_%s" % (obj, chain))
                rows.append((name, chain_sele))
    if not rows:
        raise pymol.CmdException("no atoms in the chains asked for")
    taken = set(_self.get_names("all"))
    clash = [name for name, _ in rows if name in taken]
    if clash:
        raise pymol.CmdException("already exists: " + " ".join(clash))
    for name, chain_sele in rows:
        _self.create(name, chain_sele, zoom=0)
        if group:
            _self.group(group, name, "add")
    if not int(quiet):
        print(" extract_chain: created %d objects" % len(rows))
    return [name for name, _ in rows]
