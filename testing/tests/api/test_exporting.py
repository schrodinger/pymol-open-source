import base64
import json
import os
import struct
import tempfile

from pymol import cmd
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


def _parse_glb(filepath):
    """Parse a GLB file and return (gltf_json, bin_data)."""
    with open(filepath, 'rb') as f:
        data = f.read()
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
def test_glb_export_spheres():
    """Test GLB export API with sphere representation and vertex colors"""
    cmd.pseudoatom("atoms", pos=[0.0, 0.0, 0.0], color="red")
    cmd.pseudoatom("atoms", pos=[4.0, 0.0, 0.0], color="blue")
    cmd.show_as("spheres")

    with tempfile.NamedTemporaryFile(suffix='.glb', delete=False) as f:
        glb_file = f.name

    try:
        cmd.get_glb(glb_file)
        assert os.path.getsize(glb_file) > 0

        gltf, bin_data = _parse_glb(glb_file)
        assert len(gltf['meshes']) >= 1

        prim = gltf['meshes'][0]['primitives'][0]
        pos_acc = gltf['accessors'][prim['attributes']['POSITION']]
        assert pos_acc['count'] > 0

        color_acc = gltf['accessors'][prim['attributes']['COLOR_0']]
        color_view = gltf['bufferViews'][color_acc['bufferView']]
        color_offset = color_view.get('byteOffset', 0)
        color_offset += color_acc.get('byteOffset', 0)
        colors = struct.unpack_from(
            '<%df' % (color_acc['count'] * 3), bin_data, color_offset)
        colors = list(zip(colors[0::3], colors[1::3], colors[2::3]))
        assert any(r > 0.99 and g < 0.01 and b < 0.01 for r, g, b in colors)
        assert any(r < 0.01 and g < 0.01 and b > 0.99 for r, g, b in colors)
    finally:
        if os.path.exists(glb_file):
            os.unlink(glb_file)


@test_utils.requires_version("3.2")
def test_glb_export_surface():
    """Test GLB export with surface representation"""
    cmd.fragment("gly")
    cmd.show_as("surface")

    with tempfile.NamedTemporaryFile(suffix='.glb', delete=False) as f:
        glb_file = f.name

    try:
        cmd.save(glb_file)
        gltf, _ = _parse_glb(glb_file)
        prim = gltf['meshes'][0]['primitives'][0]
        pos_acc = gltf['accessors'][prim['attributes']['POSITION']]
        assert pos_acc['count'] > 0
    finally:
        if os.path.exists(glb_file):
            os.unlink(glb_file)


@test_utils.requires_version("3.2")
def test_glb_export_cartoon():
    """Test GLB export with cartoon representation (needs enough residues)"""
    cmd.load(test_utils.datafile('1bna.cif'))
    cmd.show_as("cartoon")

    with tempfile.NamedTemporaryFile(suffix='.glb', delete=False) as f:
        glb_file = f.name

    try:
        cmd.save(glb_file)
        assert os.path.getsize(glb_file) > 0

        gltf, _ = _parse_glb(glb_file)
        assert len(gltf['meshes']) >= 1

        prim = gltf['meshes'][0]['primitives'][0]
        pos_acc = gltf['accessors'][prim['attributes']['POSITION']]
        assert pos_acc['count'] > 100  # cartoon should have many vertices
    finally:
        if os.path.exists(glb_file):
            os.unlink(glb_file)


@test_utils.requires_version("3.2")
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
def test_glb_export_empty():
    """Test GLB export with no exportable geometry produces no file"""
    cmd.fragment("ala")
    cmd.show_as("cartoon")  # too few residues for cartoon

    with tempfile.NamedTemporaryFile(suffix='.glb', delete=False) as f:
        glb_file = f.name

    # Remove the temp file so we can check if save creates one
    os.unlink(glb_file)

    cmd.save(glb_file)
    assert not os.path.exists(glb_file)
