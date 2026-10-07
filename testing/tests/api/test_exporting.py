import errno
import os
import stat
import sys
import tempfile

import pytest

from pymol import cmd
from pymol import exporting
from pymol import test_utils


@test_utils.requires_version("3.2")
def test_bcif_export():
    """Test BCIF export and round-trip"""
    # Create a simple structure
    cmd.fragment("ala")
    orig_count = cmd.count_atoms("ala")
    assert orig_count == 10

    # Export to BCIF
    with tempfile.NamedTemporaryFile(suffix='.bcif', delete=False) as f:
        bcif_file = f.name

    try:
        cmd.save(bcif_file, "ala")
        assert os.path.exists(bcif_file)
        assert os.path.getsize(bcif_file) > 0

        # Load back and verify
        cmd.delete("all")
        cmd.load(bcif_file, "test_loaded")
        loaded_count = cmd.count_atoms("test_loaded")
        assert loaded_count == orig_count, f"Atom count mismatch: {loaded_count} != {orig_count}"
    finally:
        if os.path.exists(bcif_file):
            os.unlink(bcif_file)


@test_utils.requires_version("3.2")
def test_bcif_export_multi_object():
    """Test BCIF export with multiple objects"""
    cmd.fragment("ala")
    cmd.fragment("gly")
    ala_count = cmd.count_atoms("ala")
    gly_count = cmd.count_atoms("gly")

    with tempfile.NamedTemporaryFile(suffix='.bcif', delete=False) as f:
        bcif_file = f.name

    try:
        cmd.save(bcif_file, "all")
        assert os.path.getsize(bcif_file) > 0

        cmd.delete("all")
        cmd.load(bcif_file)

        names = cmd.get_object_list()
        assert len(names) == 2, f"Expected 2 objects, got {len(names)}: {names}"
        assert cmd.count_atoms(names[0]) == ala_count
        assert cmd.count_atoms(names[1]) == gly_count
    finally:
        if os.path.exists(bcif_file):
            os.unlink(bcif_file)


@pytest.mark.parametrize("ext", ["pse", "pse.gz", "pdb"])
def test_save_failure_keeps_existing_file(ext, tmp_path, monkeypatch):
    """A failed write must not truncate an existing file (#520)"""
    filename = str(tmp_path / f"model.{ext}")
    cmd.fragment("ala")
    cmd.save(filename)
    with open(filename, "rb") as handle:
        original = handle.read()

    class DiskFull:
        def __init__(self, handle):
            self.handle = handle

        def __enter__(self):
            return self

        def __exit__(self, *exc):
            self.handle.close()

        def write(self, data):
            self.handle.write(data[:10])
            raise OSError(errno.ENOSPC, os.strerror(errno.ENOSPC))

    def failing_open(file, mode="r", *args, **kwargs):
        handle = open(file, mode, *args, **kwargs)
        return handle if "r" in mode else DiskFull(handle)

    monkeypatch.setattr(exporting, "open", failing_open, raising=False)

    cmd.fragment("gly")
    with pytest.raises(OSError):
        cmd.save(filename)

    with open(filename, "rb") as handle:
        assert handle.read() == original
    assert os.listdir(tmp_path) == [f"model.{ext}"]


@pytest.mark.skipif(sys.platform == "win32", reason="POSIX file modes")
def test_save_keeps_file_mode(tmp_path):
    filename = str(tmp_path / "model.pdb")
    cmd.fragment("ala")
    cmd.save(filename)
    os.chmod(filename, 0o640)
    cmd.save(filename)
    assert stat.S_IMODE(os.stat(filename).st_mode) == 0o640


@pytest.mark.skipif(sys.platform == "win32", reason="POSIX file modes")
@pytest.mark.skipif(hasattr(os, "geteuid") and os.geteuid() == 0,
                    reason="root can write read-only files")
def test_save_refuses_read_only_file(tmp_path):
    filename = str(tmp_path / "model.pdb")
    cmd.fragment("ala")
    cmd.save(filename)
    os.chmod(filename, 0o444)
    with pytest.raises(PermissionError):
        cmd.save(filename)


@pytest.mark.skipif(sys.platform == "win32", reason="symlinks need privileges")
def test_save_writes_through_symlink(tmp_path):
    target = tmp_path / "target.pdb"
    link = tmp_path / "link.pdb"
    target.write_text("")
    link.symlink_to(target)
    cmd.fragment("ala")
    cmd.save(str(link))
    assert link.is_symlink()
    assert target.stat().st_size > 0
