import shutil

from pymol import cmd
from pymol import test_utils
import pytest

_has_pdb2pqr = bool(
    shutil.which('pdb2pqr') or
    shutil.which('pdb2pqr30') or
    shutil.which('pdb2pqr_cli')
)


@test_utils.requires_version("3.0")
def test_look_at():
    ori_view = cmd.get_view()
    cmd.pseudoatom("M1", pos=[10, 0, 0])
    cmd.look_at("M1")
    new_view = cmd.get_view()
    assert ori_view != new_view

    ref_new_view = (0.980580688,    0.000000000,   -0.196116135,
                    0.000000000,    1.000000000,    0.000000000,
                    0.196116135,    0.000000000,    0.980580688,
                    -9.805807114,    0.000000000,  -49.029033661,
                    0.000000000,    0.000000000,    0.000000000,
                    40.000000000,  100.000000000,  -20.000000000)
    assert ref_new_view == pytest.approx(new_view)


@pytest.mark.skipif(not _has_pdb2pqr, reason="pdb2pqr not installed")
def test_protonate():
    cmd.load("testing/data/1rx1.pdb")

    # Apply visual settings to verify preservation
    cmd.color("green", "1rx1")
    cmd.show("sticks", "1rx1")
    heavy_count = cmd.count_atoms("1rx1 and not hydro")

    cmd.protonate("1rx1", pH=7.4)

    # Heavy atoms preserved
    assert cmd.count_atoms("1rx1 and not hydro") == heavy_count

    # Hydrogens were added
    assert cmd.count_atoms("1rx1 and hydro") > 0

    # Visual settings preserved on heavy atoms
    colors = set()
    cmd.iterate("1rx1 and not hydro", "colors.add(color)", space=locals())
    assert 3 in colors  # 3 = green color index

    # At pH 7.4, most Asp carboxylates should be deprotonated
    # (PROPKA may predict borderline pKa for some residues)
    n_asp = cmd.count_atoms("1rx1 and not hydro and resn ASP and name OD1+OD2")
    asp_carboxyl_h = cmd.count_atoms(
        "1rx1 and hydro and neighbor (resn ASP and name OD1+OD2)")
    assert asp_carboxyl_h < n_asp, \
        "Most Asp carboxylates should be deprotonated at pH 7.4"


@pytest.mark.skipif(not _has_pdb2pqr, reason="pdb2pqr not installed")
def test_protonate_low_pH():
    cmd.load("testing/data/1rx1.pdb")

    cmd.protonate("1rx1", pH=2.0)

    # At pH 2.0, Asp carboxylates should be protonated
    asp_carboxyl_h = cmd.count_atoms(
        "1rx1 and hydro and neighbor (resn ASP and name OD1+OD2)")
    assert asp_carboxyl_h > 0, \
        "Asp carboxylates should be protonated at pH 2.0"


def test_protonate_no_object():
    with pytest.raises(Exception):
        cmd.protonate("nonexistent_object")


def test_protonate_invalid_pH():
    cmd.fragment("gly")
    with pytest.raises(Exception):
        cmd.protonate("gly", pH=-1.0)
    with pytest.raises(Exception):
        cmd.protonate("gly", pH=15.0)


def test_protonate_fallback():
    """Test textbook pKa fallback (no pdb2pqr needed)."""
    from pymol.editing import _protonate_fallback

    cmd.load("testing/data/1rx1.pdb")
    cmd.color("green", "1rx1")
    heavy_count = cmd.count_atoms("1rx1 and not hydro")

    # Use fallback directly at pH 7.4
    _protonate_fallback("all", "1rx1", 7.4, 0, 1, _self=cmd)

    # Heavy atoms preserved
    assert cmd.count_atoms("1rx1 and not hydro") == heavy_count

    # Hydrogens were added
    assert cmd.count_atoms("1rx1 and hydro") > 0

    # Colors preserved
    colors = set()
    cmd.iterate("1rx1 and not hydro", "colors.add(color)", space=locals())
    assert 3 in colors

    # At pH 7.4, Asp carboxylates deprotonated (pKa 3.65)
    asp_h = cmd.count_atoms(
        "1rx1 and hydro and neighbor (resn ASP and name OD1+OD2)")
    assert asp_h == 0, \
        "Asp carboxylates should be deprotonated at pH 7.4"

    # At pH 7.4, His deprotonated (pKa 6.0) — no H on ND1
    his_nd1_h = cmd.count_atoms(
        "1rx1 and hydro and neighbor (resn HIS and name ND1)")
    assert his_nd1_h == 0, \
        "His ND1 should be deprotonated at pH 7.4"

    # At pH 7.4, Lys protonated (pKa 10.53) — 3H on NZ
    lys_h = cmd.count_atoms(
        "1rx1 and hydro and neighbor (resn LYS and name NZ)")
    n_lys = cmd.count_atoms("1rx1 and not hydro and resn LYS and name NZ")
    assert lys_h == n_lys * 3, \
        "Lys NZ should have 3H (protonated) at pH 7.4"


def test_protonate_fallback_low_pH():
    """Test fallback at low pH — carboxylates should be protonated."""
    from pymol.editing import _protonate_fallback

    cmd.load("testing/data/1rx1.pdb")

    _protonate_fallback("all", "1rx1", 2.0, 0, 1, _self=cmd)

    # Asp carboxylates protonated (pH 2.0 < pKa 3.65)
    asp_h = cmd.count_atoms(
        "1rx1 and hydro and neighbor (resn ASP and name OD1+OD2)")
    assert asp_h > 0, \
        "Asp carboxylates should be protonated at pH 2.0"

    # His protonated (pH 2.0 < pKa 6.0) — H on ND1
    his_nd1_h = cmd.count_atoms(
        "1rx1 and hydro and neighbor (resn HIS and name ND1)")
    n_his = cmd.count_atoms("1rx1 and not hydro and resn HIS and name ND1")
    assert his_nd1_h == n_his, \
        "His ND1 should be protonated at pH 2.0"


def test_protonate_fallback_pH_transitions():
    """Formal charges should drive h_add across textbook pKa thresholds."""
    from pymol.editing import _protonate_fallback

    cmd.fab("EHK", "m1")

    expected = [
        (3.5, 30),
        (5.9, 29),
        (6.1, 28),
        (13.0, 27),
    ]

    for pH, hydrogen_count in expected:
        _protonate_fallback("m1", "m1", pH, 0, 1, _self=cmd)
        assert cmd.count_atoms("m1 and hydro") == hydrogen_count


def _formal_charges(selection):
    charges = []
    cmd.iterate(selection, "charges.append(formal_charge)",
                space={"charges": charges})
    return charges


def _resns(selection):
    names = set()
    cmd.iterate(selection, "names.add(resn)", space={"names": names})
    return names


def test_protonate_fallback_free_cysteine():
    """A free thiol loses its H and turns into a thiolate above pKa 8.18."""
    from pymol.editing import _protonate_fallback

    cmd.fab("ACA", "m1")

    _protonate_fallback("m1", "m1", 7.0, 0, 1, _self=cmd)
    assert cmd.count_atoms("m1 and hydro and neighbor (name SG)") == 1
    assert _formal_charges("m1 and name SG") == [0]

    _protonate_fallback("m1", "m1", 9.0, 0, 1, _self=cmd)
    assert cmd.count_atoms("m1 and hydro and neighbor (name SG)") == 0
    assert _formal_charges("m1 and name SG") == [-1]

    # and back down again
    _protonate_fallback("m1", "m1", 7.0, 0, 1, _self=cmd)
    assert cmd.count_atoms("m1 and hydro and neighbor (name SG)") == 1
    assert _formal_charges("m1 and name SG") == [0]


def test_protonate_fallback_disulfide_stays_neutral():
    """A bridged cysteine has no ionizable H, so it must not pick up -1."""
    from pymol.editing import _protonate_fallback

    cmd.fab("CGC", "m1")
    cmd.bond("m1 and resi 1 and name SG", "m1 and resi 3 and name SG")

    _protonate_fallback("m1", "m1", 13.0, 0, 1, _self=cmd)

    assert cmd.count_atoms("m1 and hydro and neighbor (name SG)") == 0
    assert _formal_charges("m1 and name SG") == [0, 0]


def test_protonate_fallback_disulfide_in_discrete_object():
    """A discrete object keeps each atom in one state only, so the bridged
    cysteines outside the current state must still be recognized."""
    from pymol.editing import _protonate_fallback

    cmd.fab("CGC", "m1")
    cmd.bond("m1 and resi 1 and name SG", "m1 and resi 3 and name SG")
    cmd.create("m2", "m1", 1, 1, discrete=1)
    cmd.create("m2", "m2", 1, 2, discrete=1)

    _protonate_fallback("m2", "m2", 13.0, 0, 1, _self=cmd)

    assert _formal_charges("m2 and name SG") == [0, 0, 0, 0]


def test_protonate_fallback_tyrosine_and_arginine():
    """Tyr and Arg bracket the pH range covered by the fallback table."""
    from pymol.editing import _protonate_fallback

    cmd.fab("YR", "m1")

    _protonate_fallback("m1", "m1", 7.0, 0, 1, _self=cmd)
    assert cmd.count_atoms("m1 and hydro and neighbor (name OH)") == 1
    assert _formal_charges("m1 and name OH") == [0]
    assert cmd.count_atoms(
        "m1 and hydro and neighbor (name NE+NH1+NH2)") == 5
    assert _formal_charges("m1 and name NH1") == [1]

    _protonate_fallback("m1", "m1", 13.0, 0, 1, _self=cmd)
    assert cmd.count_atoms("m1 and hydro and neighbor (name OH)") == 0
    assert _formal_charges("m1 and name OH") == [-1]
    assert cmd.count_atoms(
        "m1 and hydro and neighbor (name NE+NH1+NH2)") == 4
    assert _formal_charges("m1 and name NH1") == [0]


def _build_pdb_residue(seq, resn, obj_name):
    """Build `seq`, rename it to `resn` and round-trip it through PDB so that
    the reader assigns the bond orders and formal charges that name implies."""
    cmd.fab(seq, "_prot_tmp")
    cmd.alter("_prot_tmp", f"resn = {resn!r}")
    pdbstr = cmd.get_pdbstr("_prot_tmp and not hydro")
    cmd.delete("_prot_tmp")
    cmd.read_pdbstr(pdbstr, obj_name)


def _build_his_tautomer(resn, obj_name):
    _build_pdb_residue("H", resn, obj_name)


@pytest.mark.parametrize("resn,neutral_h_on", [
    ("HIS", "NE2"),
    ("HID", "ND1"),
    ("HIE", "NE2"),
])
def test_protonate_fallback_his_tautomers(resn, neutral_h_on):
    """Both imidazole tautomers must titrate on their own free nitrogen."""
    from pymol.editing import _protonate_fallback

    other = "ND1" if neutral_h_on == "NE2" else "NE2"
    _build_his_tautomer(resn, "m1")

    # above pKa 6.0: neutral, single ring H on the tautomer's nitrogen
    _protonate_fallback("m1", "m1", 7.0, 0, 1, _self=cmd)
    assert cmd.count_atoms(
        f"m1 and hydro and neighbor (name {neutral_h_on})") == 1
    assert cmd.count_atoms(f"m1 and hydro and neighbor (name {other})") == 0
    assert _formal_charges("m1 and name ND1+NE2") == [0, 0]

    # below pKa 6.0: imidazolium, one H on each ring nitrogen
    _protonate_fallback("m1", "m1", 5.0, 0, 1, _self=cmd)
    assert cmd.count_atoms(
        f"m1 and hydro and neighbor (name {neutral_h_on})") == 1
    assert cmd.count_atoms(f"m1 and hydro and neighbor (name {other})") == 1
    assert sum(_formal_charges("m1 and name ND1+NE2")) == 1


@pytest.mark.parametrize("resn", ["HIS", "HISB", "HISE", "HIE", "HIP",
                                  "HISH", "HISP"])
def test_protonate_his_epsilon_aliases(resn):
    """Names the PDB reader double bonds CE1=ND1 titrate on ND1."""
    from pymol.editing import _get_formal_charge

    assert _get_formal_charge(resn, 'ND1', 0, 5.0) == 1
    assert _get_formal_charge(resn, 'NE2', 0, 5.0) == 0
    assert _get_formal_charge(resn, 'ND1', 0, 7.0) == 0
    assert _get_formal_charge(resn, 'NE2', 0, 7.0) == 0


@pytest.mark.parametrize("resn", ["HID", "HISA", "HISD"])
def test_protonate_his_delta_aliases(resn):
    """Names the PDB reader double bonds CE1=NE2 titrate on NE2."""
    from pymol.editing import _get_formal_charge

    assert _get_formal_charge(resn, 'NE2', 0, 5.0) == 1
    assert _get_formal_charge(resn, 'ND1', 0, 5.0) == 0
    assert _get_formal_charge(resn, 'NE2', 0, 7.0) == 0
    assert _get_formal_charge(resn, 'ND1', 0, 7.0) == 0


def test_protonate_titratable_selection_matches_table():
    """The alter selection must match an atom exactly when
    _get_formal_charge has a pKa for it - missing atoms are never titrated,
    extra atoms lose their chemFlag for nothing."""
    from pymol.editing import (_RESN_ALIAS, _TITRATABLE_SELECTION,
                               _get_formal_charge)

    sentinel = object()

    cmd.fab("EHKRYCDA", "m1")
    # the residue of "EHKRYCDA" whose atom names each canonical name expects
    resi_of = {'GLU': 1, 'HIE': 2, 'HID': 2, 'LYS': 3,
               'ARG': 4, 'TYR': 5, 'CYS': 6, 'ASP': 7}
    seen, matched = set(), set()

    def snapshot():
        cmd.iterate("m1", "acc.add((resn, name))", space={"acc": seen})
        cmd.iterate(f"(m1) and ({_TITRATABLE_SELECTION})",
                    "acc.add((resn, name))", space={"acc": matched})

    snapshot()
    for alias, canonical in sorted(_RESN_ALIAS.items()):
        cmd.alter(f"m1 and resi {resi_of[canonical]}", f"resn = {alias!r}")
        snapshot()

    expected = {
        (resn, name) for resn, name in seen
        if _get_formal_charge(resn, name, sentinel, 7.0) is not sentinel
    }

    assert matched
    assert matched == expected


def test_protonate_fallback_lowercase_resn():
    """PyMOL selections ignore case, so the pKa lookup must too."""
    from pymol.editing import _protonate_fallback, _get_formal_charge

    assert _get_formal_charge('his', 'nd1', 0, 5.0) == 1
    assert _get_formal_charge('asp', 'od2', 0, 7.0) == -1

    cmd.fab("D", "m1")
    cmd.alter("m1", "resn = resn.lower()")
    cmd.alter("m1", "name = name.lower()")

    _protonate_fallback("m1", "m1", 7.0, 0, 1, _self=cmd)

    assert cmd.count_atoms("m1 and hydro and neighbor (name OD1+OD2)") == 0
    assert _formal_charges("m1 and name OD2") == [-1]


def test_protonate_fallback_preserves_untitrated_chemistry():
    """Writing formal_charge clears chemFlag, so only titratable atoms may
    be written - otherwise set_geometry corrections are silently discarded."""
    from pymol.editing import _protonate_fallback

    cmd.fab("GD", "m1")
    cmd.remove("m1 and hydro")
    # planar with valence 3: CA has N and C as neighbours, so exactly one H
    cmd.set_geometry("m1 and resi 1 and name CA", 2, 3)

    _protonate_fallback("m1", "m1", 7.0, 0, 1, _self=cmd)

    assert cmd.count_atoms(
        "m1 and hydro and neighbor (resi 1 and name CA)") == 1
    assert cmd.count_atoms("m1 and hydro and neighbor (name OD1+OD2)") == 0
    assert _formal_charges("m1 and name OD2") == [-1]


@pytest.mark.parametrize("seq,alias,name,acidic,basic", [
    ("R", "ARGP", "NH1", 1, 0),
    ("K", "LYSP", "NZ", 1, 0),
    ("D", "ASPM", "OD2", 0, -1),
    ("E", "GLUM", "OE2", 0, -1),
])
def test_protonate_fallback_reader_charged_aliases(seq, alias, name,
                                                   acidic, basic):
    """The PDB reader pre-charges the ARGP/ASPM/GLUM/LYSP spellings, so the
    fallback has to titrate them instead of leaving the reader's charge."""
    from pymol.editing import _protonate_fallback

    _build_pdb_residue(seq, alias, "m1")
    assert _resns("m1") == {alias}

    _protonate_fallback("m1", "m1", 1.0, 0, 1, _self=cmd)
    assert _formal_charges(f"m1 and name {name}") == [acidic]

    _protonate_fallback("m1", "m1", 13.0, 0, 1, _self=cmd)
    assert _formal_charges(f"m1 and name {name}") == [basic]


@pytest.mark.parametrize("resn", ["HIP", "HISH", "HISP"])
def test_protonate_fallback_demoted_cation_is_renamed(resn):
    """The PDB reader restores ND1 = +1 for these names whatever the file
    says, so a demoted residue has to be renamed to survive a round trip."""
    from pymol.editing import _protonate_fallback

    _build_his_tautomer(resn, "m1")

    # below pKa 6.0 the name still describes the residue, so it is kept
    _protonate_fallback("m1", "m1", 5.0, 0, 1, _self=cmd)
    assert _resns("m1") == {resn}
    assert sum(_formal_charges("m1 and name ND1+NE2")) == 1
    assert cmd.count_atoms("m1 and hydro and neighbor (name ND1+NE2)") == 2

    # above it, the name would resurrect the proton the fallback removed
    _protonate_fallback("m1", "m1", 7.4, 0, 1, _self=cmd)
    assert _resns("m1") == {"HIS"}
    assert _formal_charges("m1 and name ND1+NE2") == [0, 0]

    cmd.read_pdbstr(cmd.get_pdbstr("m1"), "m2")
    assert _formal_charges("m2 and name ND1+NE2") == [0, 0]
    assert cmd.count_atoms("m2 and hydro and neighbor (name ND1)") == 0


def test_protonate_fallback_thioether_crosslink_stays_neutral():
    """The free-thiol test has to be connectivity based: an SG bridged to
    another residue's CB is no more ionizable than a disulfide."""
    from pymol.editing import _protonate_fallback

    cmd.fab("CGC", "m1")
    cmd.bond("m1 and resi 1 and name SG", "m1 and resi 3 and name CB")

    _protonate_fallback("m1", "m1", 13.0, 0, 1, _self=cmd)

    assert _formal_charges("m1 and resi 1 and name SG") == [0]
    assert cmd.count_atoms(
        "m1 and hydro and neighbor (resi 1 and name SG)") == 0
    # the cysteine which is still free keeps titrating
    assert _formal_charges("m1 and resi 3 and name SG") == [-1]


def test_protonate_fallback_keeps_unselected_hydrogens():
    """Only the selection is rebuilt, so only its hydrogens may be stripped."""
    from pymol.editing import _protonate_fallback

    cmd.fab("AD", "m1")
    cmd.h_add("m1")
    untouched = cmd.count_atoms("m1 and hydro and resi 1")
    assert untouched > 0

    _protonate_fallback("m1 and resi 2", "m1", 13.0, 0, 1, _self=cmd)

    assert cmd.count_atoms("m1 and hydro and resi 1") == untouched
    assert _formal_charges("m1 and name OD2") == [-1]
    assert cmd.count_atoms("m1 and hydro and neighbor (name OD1+OD2)") == 0


def _build_single_nuc(nuc_acid, nuc_type, obj_name, chain='A'):
    """Helper: build a single nucleotide and return the object name.

    Creates a nucleotide from scratch via attach_nuc_acid, then optionally
    alters the chain to the requested value.
    """
    from pymol import editor
    # Build from scratch (0 atoms selected) — always creates chain 'A'
    cmd.select("sele", "none")
    editor.attach_nuc_acid("sele", nuc_acid, nuc_type, object=obj_name, dbl_helix=False)
    if chain != 'A':
        cmd.alter(obj_name, f"chain='{chain}';segi='{chain}'")
        cmd.rebuild()
    return obj_name


def test_attach_nuc_acid_extend_with_chain():
    """attach_nuc_acid can extend a nucleotide that has a chain ID."""
    from pymol import editor

    obj = _build_single_nuc("atp", "RNA", "test_nuc1", chain="A")
    initial_count = cmd.count_atoms(obj)

    # Select O3' on chain A to extend
    cmd.select("sele", f"{obj} & name O3' & chain A")
    assert cmd.count_atoms("sele") == 1

    editor.attach_nuc_acid("sele", "utp", "RNA", dbl_helix=False)
    assert cmd.count_atoms(obj) > initial_count


def test_attach_nuc_acid_extend_empty_chain():
    """attach_nuc_acid should work when chain ID is empty (GH #502)."""
    from pymol import editor

    obj = _build_single_nuc("atp", "RNA", "test_nuc2", chain="")
    initial_count = cmd.count_atoms(obj)

    # Select O3' on empty chain to extend
    cmd.select("sele", f"{obj} & name O3' & chain ''")
    assert cmd.count_atoms("sele") == 1

    editor.attach_nuc_acid("sele", "utp", "RNA", dbl_helix=False)
    assert cmd.count_atoms(obj) > initial_count
