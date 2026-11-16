"""
Load every structure file of a folder; save many objects to one file each.
Ported from the pymolish plugin (MIT, same author); new code, MIT.
"""

import fnmatch
import os
import pathlib

import pymol

cmd = __import__("sys").modules["pymol.cmd"]


def find_files(directory=".", format="pdb", pattern="*", recursive=0):
    """Sorted paths of the files "load_files" would load."""
    root = pathlib.Path(cmd.exp_path(directory))
    if not root.is_dir():
        raise pymol.CmdException("not a directory: %s" % directory)
    suffix = "." + format.lower().lstrip(".")
    return sorted(
        (
            path
            for path in (root.rglob("*") if int(recursive) else root.iterdir())
            if path.suffix.lower() == suffix
            and fnmatch.fnmatch(path.stem, pattern)
            and path.is_file()
        ),
        # Windows paths compare without case: same order on every system
        key=lambda path: path.parts,
    )


def load_files(
    directory=".",
    format="pdb",
    group="",
    pattern="*",
    recursive=0,
    prefix="",
    quiet=1,
    *,
    _self=cmd,
):
    """
    DESCRIPTION

        "load_files" loads every file of one format from a directory, one
        object per file.

    USAGE

        load_files [ directory [, format [, group [, pattern [, recursive
            [, prefix ]]]]]]

    ARGUMENTS

        directory = string: where to look {default: working directory}

        format = string: file extension, any that "load" reads {default: pdb}

        group = string: group to put the objects in {default: none}

        pattern = string: wildcard for the file name without extension
        {default: *}

        recursive = 0/1: also look in subdirectories; their names become
        part of the object name {default: 0}

        prefix = string: put in front of every object name

    EXAMPLE

        load_files ~/poses, sdf, poses
        load_files ., cif, pattern=model_*, recursive=1

    NOTES

        A file that cannot be read is reported and skipped. An object is
        never replaced: a name in use gets a number.

    SEE ALSO

        load, loadall, export_objects, group
    """
    paths = find_files(directory, format, pattern, recursive)
    if not paths:
        raise pymol.CmdException("no %s files match in %s" % (format, directory))
    root = pathlib.Path(cmd.exp_path(directory))
    loaded = []
    known = set(_self.get_names("objects"))
    for path in paths:
        stem = str(path.relative_to(root).with_suffix("")).replace(os.sep, "_")
        name = _self.get_unused_name(_self.get_legal_name(prefix + stem), 0)
        try:
            _self.load(str(path), name, zoom=0, quiet=1)
        except pymol.CmdException as error:
            print(" load_files: skipped %s: %s" % (path, error))
            continue
        # ask for the names: "b.pdb" is loaded as "b_", "b" is a reserved word
        for name in [n for n in _self.get_names("objects") if n not in known]:
            if _self.get_type(name) == "object:molecule" and not (
                _self.count_atoms("%" + name)
            ):
                _self.delete(name)
                print(" load_files: skipped %s: no atoms" % path)
            else:
                known.add(name)
                loaded.append(name)
    if group and loaded:
        _self.group(group, " ".join(loaded))
    if _self.get_setting_boolean("auto_zoom"):
        _self.zoom()
    if not int(quiet):
        print(" load_files: loaded %d of %d files" % (len(loaded), len(paths)))
    return loaded


def export_objects(
    selection, directory=".", format="cif", overwrite=0, quiet=1, *, _self=cmd
):
    """
    DESCRIPTION

        "export_objects" saves every object of a selection to its own file,
        <directory>/<object>.<format>.

    USAGE

        export_objects selection [, directory [, format [, overwrite ]]]

    ARGUMENTS

        selection = string: a group name, a name pattern or any atom
        selection; of each object the selected atoms are saved

        directory = string: created when missing {default: working
        directory}

        format = string: file extension, any that "save" writes {default: cif}

        overwrite = 0/1: replace files that exist {default: 0, skip them}

    EXAMPLE

        export_objects poses, ~/out, pdb
        export_objects model_* and polymer, format=fasta

    SEE ALSO

        save, multisave, load_files
    """
    names = _self.get_object_list(selection)
    if not names:
        raise pymol.CmdException("no molecular objects in '%s'" % selection)
    root = pathlib.Path(cmd.exp_path(directory))
    root.mkdir(parents=True, exist_ok=True)
    written = []
    for name in names:
        path = root / ("%s.%s" % (name, format.lower().lstrip(".")))
        if path.exists() and not int(overwrite):
            print(" export_objects: %s exists, skipped" % path)
            continue
        _self.save(str(path), "(%s) and %%%s" % (selection, name), quiet=1)
        written.append(str(path))
    if not int(quiet):
        print(" export_objects: wrote %d of %d files" % (len(written), len(names)))
    return written
