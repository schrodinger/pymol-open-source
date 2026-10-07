import base64
import errno
import json
import os
import stat
import struct
import sys
import tempfile

import numpy
import pytest

import pymol
from pymol import cmd
from pymol import exporting
from pymol import test_utils

requires_gltf = test_utils.requires_capability(
    "gltf", "Requires PyMOL built with --json=true (native glTF/GLB export)")


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


def _srgb_to_linear(c):
    """Reference implementation of the glTF sRGB to linear transfer."""
    if c <= 0.0:
        return 0.0
    if c >= 1.0:
        return 1.0
    if c <= 0.04045:
        return c / 12.92
    return ((c + 0.055) / 1.055) ** 2.4


def _read_accessor(gltf, bin_data, accessor_index):
    """Return an accessor's values as an array with one row per element."""
    acc = gltf['accessors'][accessor_index]
    view = gltf['bufferViews'][acc['bufferView']]
    offset = view.get('byteOffset', 0) + acc.get('byteOffset', 0)
    ncomp = {'SCALAR': 1, 'VEC3': 3}[acc['type']]
    dtype = {5126: '<f4', 5125: '<u4'}[acc['componentType']]
    values = numpy.frombuffer(bin_data, dtype, acc['count'] * ncomp, offset)
    return values.reshape(acc['count'], ncomp)


def _read_attribute(gltf, bin_data, name, mesh_index=0):
    """Return a mesh's vertex attribute (e.g. COLOR_0) as an array."""
    prim = gltf['meshes'][mesh_index]['primitives'][0]
    return _read_accessor(gltf, bin_data, prim['attributes'][name])


def _parse_glb(filepath):
    """Parse a GLB file and return (gltf_json, bin_data)."""
    with open(filepath, 'rb') as f:
        return _parse_glb_bytes(f.read())


def _parse_glb_bytes(data):
    """Parse GLB data and return (gltf_json, bin_data)."""
    magic, version, length = struct.unpack_from('<III', data, 0)
    assert magic == 0x46546C67, f"Bad GLB magic: {hex(magic)}"
    assert version == 2
    assert length == len(data)

    json_len, json_type = struct.unpack_from('<II', data, 12)
    assert json_type == 0x4E4F534A  # "JSON"
    gltf = json.loads(data[20:20 + json_len])

    bin_offset = 20 + json_len
    bin_len, bin_type = struct.unpack_from('<II', data, bin_offset)
    assert bin_type == 0x004E4942  # "BIN\0"
    bin_data = data[bin_offset + 8:bin_offset + 8 + bin_len]

    return gltf, bin_data


@test_utils.requires_version("3.2")
@requires_gltf
def test_glb_export_sticks():
    """Test GLB export with stick representation"""
    cmd.fragment("ala")
    cmd.show_as("sticks")

    with tempfile.NamedTemporaryFile(suffix='.glb', delete=False) as f:
        glb_file = f.name

    try:
        cmd.save(glb_file)
        assert os.path.exists(glb_file)
        assert os.path.getsize(glb_file) > 0

        gltf, bin_data = _parse_glb(glb_file)

        # Validate glTF structure
        assert gltf['asset']['version'] == '2.0'
        assert 'PyMOL' in gltf['asset']['generator']
        assert len(gltf['scenes']) == 1
        assert len(gltf['meshes']) >= 1
        assert len(gltf['materials']) >= 1
        assert len(gltf['buffers']) == 1
        assert len(bin_data) > 0

        # PyMOL renders without back-face culling
        assert all(mat['doubleSided'] for mat in gltf['materials'])

        # Check mesh has required attributes
        prim = gltf['meshes'][0]['primitives'][0]
        assert 'POSITION' in prim['attributes']
        assert 'NORMAL' in prim['attributes']
        assert 'COLOR_0' in prim['attributes']
        assert 'indices' in prim
        assert prim['mode'] == 4  # TRIANGLES

        # Check accessor types
        pos_acc = gltf['accessors'][prim['attributes']['POSITION']]
        assert pos_acc['type'] == 'VEC3'
        assert pos_acc['componentType'] == 5126  # FLOAT
        assert pos_acc['count'] > 0
        assert 'min' in pos_acc
        assert 'max' in pos_acc
    finally:
        if os.path.exists(glb_file):
            os.unlink(glb_file)


@test_utils.requires_version("3.2")
@requires_gltf
def test_gltf_export_sticks():
    """Test native, self-contained glTF 2.0 export"""
    cmd.fragment("ala")
    cmd.show_as("sticks")

    with tempfile.NamedTemporaryFile(suffix='.gltf', delete=False) as f:
        gltf_file = f.name

    try:
        cmd.save(gltf_file)

        with open(gltf_file, encoding='utf-8') as handle:
            gltf = json.load(handle)

        assert gltf['asset']['version'] == '2.0'
        assert len(gltf['meshes']) >= 1

        buffer = gltf['buffers'][0]
        prefix = 'data:application/octet-stream;base64,'
        assert buffer['uri'].startswith(prefix)
        binary = base64.b64decode(buffer['uri'][len(prefix):])
        assert len(binary) == buffer['byteLength']
    finally:
        if os.path.exists(gltf_file):
            os.unlink(gltf_file)


@test_utils.requires_version("3.2")
@requires_gltf
def test_glb_export_spheres():
    """Test GLB export API with sphere representation and vertex colors"""
    cmd.pseudoatom("atoms", pos=[0.0, 0.0, 0.0], color="red")
    cmd.pseudoatom("atoms", pos=[4.0, 0.0, 0.0], color="blue")
    cmd.show_as("spheres")

    gltf, bin_data = _parse_glb_bytes(cmd.get_glb())
    assert len(gltf['meshes']) >= 1

    prim = gltf['meshes'][0]['primitives'][0]
    pos_acc = gltf['accessors'][prim['attributes']['POSITION']]
    assert pos_acc['count'] > 0

    colors = _read_attribute(gltf, bin_data, 'COLOR_0')
    assert any(r > 0.99 and g < 0.01 and b < 0.01 for r, g, b in colors)
    assert any(r < 0.01 and g < 0.01 and b > 0.99 for r, g, b in colors)


@test_utils.requires_version("3.2")
@requires_gltf
def test_glb_export_midtone_color_is_linear():
    """Mid-tone colors are converted from PyMOL display space to linear"""
    cmd.pseudoatom("atoms", pos=[0.0, 0.0, 0.0], color="grey50")
    cmd.show_as("spheres")

    display = cmd.get_color_tuple("grey50")
    expected = tuple(_srgb_to_linear(c) for c in display)

    # only meaningful if the color actually is a mid-tone
    assert all(0.01 < c < 0.99 for c in display)
    assert all(abs(e - d) > 0.05 for e, d in zip(expected, display))

    with tempfile.NamedTemporaryFile(suffix='.glb', delete=False) as f:
        glb_file = f.name

    try:
        cmd.save(glb_file)
        gltf, bin_data = _parse_glb(glb_file)

        colors = _read_attribute(gltf, bin_data, 'COLOR_0')
        assert len(colors)
        for color in colors:
            assert all(abs(c - e) < 1e-5 for c, e in zip(color, expected)), \
                f"{color} != {expected}"
    finally:
        if os.path.exists(glb_file):
            os.unlink(glb_file)


@test_utils.requires_version("3.2")
@requires_gltf
def test_glb_export_surface():
    """Test GLB export with surface representation"""
    cmd.fragment("gly")
    cmd.show_as("surface")

    with tempfile.NamedTemporaryFile(suffix='.glb', delete=False) as f:
        glb_file = f.name

    try:
        cmd.save(glb_file)
        gltf, bin_data = _parse_glb(glb_file)
        prim = gltf['meshes'][0]['primitives'][0]
        pos_acc = gltf['accessors'][prim['attributes']['POSITION']]
        assert pos_acc['count'] > 0

        # adjacent triangles share vertices, a closed surface has about
        # half as many vertices as triangles (unshared would be 3x)
        tri = _read_accessor(gltf, bin_data, prim['indices']).reshape(-1, 3)
        assert pos_acc['count'] < len(tri)
        assert len(numpy.unique(tri)) == pos_acc['count']
        assert (tri[:, 0] != tri[:, 1]).all()
        assert (tri[:, 1] != tri[:, 2]).all()
        assert (tri[:, 0] != tri[:, 2]).all()
    finally:
        if os.path.exists(glb_file):
            os.unlink(glb_file)


@test_utils.requires_version("3.2")
@requires_gltf
def test_glb_export_cartoon():
    """Test GLB export with cartoon representation (needs enough residues)"""
    cmd.load(test_utils.datafile('1bna.cif'))
    cmd.show_as("cartoon")

    with tempfile.NamedTemporaryFile(suffix='.glb', delete=False) as f:
        glb_file = f.name

    try:
        cmd.save(glb_file)
        assert os.path.getsize(glb_file) > 0

        gltf, bin_data = _parse_glb(glb_file)
        assert len(gltf['meshes']) >= 1

        prim = gltf['meshes'][0]['primitives'][0]
        pos_acc = gltf['accessors'][prim['attributes']['POSITION']]
        assert pos_acc['count'] > 100  # cartoon should have many vertices

        normals = _read_attribute(gltf, bin_data, 'NORMAL')
        assert len(normals) == pos_acc['count']
        for normal in normals:
            length_squared = sum(component * component for component in normal)
            assert abs(length_squared - 1.0) < 1e-5, normal
    finally:
        if os.path.exists(glb_file):
            os.unlink(glb_file)


@test_utils.requires_version("3.2")
@requires_gltf
def test_glb_export_transparent():
    """Test GLB export with transparency creates BLEND material"""
    cmd.fragment("ala")
    cmd.show_as("spheres")
    cmd.set("sphere_transparency", 0.5)

    with tempfile.NamedTemporaryFile(suffix='.glb', delete=False) as f:
        glb_file = f.name

    try:
        cmd.save(glb_file)
        assert os.path.getsize(glb_file) > 0

        gltf, _ = _parse_glb(glb_file)

        # Should have a material with BLEND alpha mode
        has_blend = any(
            mat.get('alphaMode') == 'BLEND'
            for mat in gltf['materials']
        )
        assert has_blend, "Transparent geometry should produce BLEND material"
    finally:
        if os.path.exists(glb_file):
            os.unlink(glb_file)


@test_utils.requires_version("3.2")
@requires_gltf
def test_glb_export_winding():
    """Triangles face the same way as their vertex normals"""
    # viewers flip the normals of back facing triangles on double sided
    # materials, so a wrong winding renders dark
    cmd.fragment("trp")
    cmd.show_as("sticks")
    cmd.show("spheres")
    cmd.set("sphere_scale", 0.3)

    gltf, bin_data = _parse_glb_bytes(cmd.get_glb())

    for mesh_index, mesh in enumerate(gltf['meshes']):
        pos = _read_attribute(gltf, bin_data, 'POSITION', mesh_index)
        nrm = _read_attribute(gltf, bin_data, 'NORMAL', mesh_index)
        tri = _read_accessor(gltf, bin_data,
                             mesh['primitives'][0]['indices']).reshape(-1, 3)
        face = numpy.cross(pos[tri[:, 1]] - pos[tri[:, 0]],
                           pos[tri[:, 2]] - pos[tri[:, 0]])
        vertex = nrm[tri[:, 0]] + nrm[tri[:, 1]] + nrm[tri[:, 2]]
        nondegenerate = numpy.linalg.norm(face, axis=1) > 1e-6
        dots = numpy.einsum('ij,ij->i', face, vertex)[nondegenerate]
        assert len(dots) > 0
        assert (dots > 0).all(), f"{(dots <= 0).sum()} inverted triangles"


@test_utils.requires_version("3.2")
@requires_gltf
def test_glb_export_view_orientation():
    """Scene is oriented like the view, unless geometry_export_mode=1"""
    cmd.pseudoatom("atoms", pos=[0.0, 0.0, 0.0])
    cmd.pseudoatom("atoms", pos=[10.0, 0.0, 0.0])
    cmd.show_as("spheres")
    cmd.zoom()
    cmd.turn("y", 90)

    def get_bounds():
        gltf, _ = _parse_glb_bytes(cmd.get_glb())
        prim = gltf['meshes'][0]['primitives'][0]
        pos_acc = gltf['accessors'][prim['attributes']['POSITION']]
        return pos_acc['min'], pos_acc['max']

    # model x axis points along view z, centered on the origin of rotation
    lo, hi = get_bounds()
    assert hi[0] - lo[0] < 5
    assert hi[2] - lo[2] > 10
    assert abs(hi[2] + lo[2]) < 1e-3

    cmd.set("geometry_export_mode", 1)
    lo, hi = get_bounds()
    assert hi[0] - lo[0] > 10
    assert hi[2] - lo[2] < 5
    assert abs(hi[0] + lo[0] - 10.0) < 1e-3


@test_utils.requires_version("3.2")
@requires_gltf
@pytest.mark.parametrize("suffix", [".glb", ".gltf"])
def test_glb_export_empty(suffix):
    """Export with no exportable geometry raises and writes no file"""
    cmd.fragment("ala")
    cmd.show_as("cartoon")  # too few residues for cartoon

    with tempfile.NamedTemporaryFile(suffix=suffix, delete=False) as f:
        out_file = f.name

    # Remove the temp file so we can check if save creates one
    os.unlink(out_file)

    try:
        with pytest.raises(pymol.CmdException):
            cmd.save(out_file)
        assert not os.path.exists(out_file)
    finally:
        if os.path.exists(out_file):
            os.unlink(out_file)


@pytest.mark.parametrize("ext", [
    "pse", "pse.gz", "pdb",
    pytest.param("glb", marks=requires_gltf),
    pytest.param("gltf", marks=requires_gltf),
])
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
