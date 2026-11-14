"""
Molecules from SMILES, with 3D coordinates from rdkit (optional dependency,
extra "chem"). Ported from the pymolish plugin (MIT, same author); new
code, MIT.
"""

import os
import re

import pymol

cmd = __import__("sys").modules["pymol.cmd"]


def rdkit():
    """(Chem, AllChem), or a CmdException that names the extra."""
    try:
        from rdkit import Chem, RDLogger
        from rdkit.Chem import AllChem
    except ImportError:
        raise pymol.CmdException(
            "this command needs rdkit: pip install 'pymolish[chem]'"
        ) from None
    RDLogger.DisableLog("rdApp.*")  # parse errors are reported here instead
    return Chem, AllChem


def _entries(text):
    """(SMILES, name) pairs: one from the text itself, or one per line of a
    file with the SMILES in the first column and a name in the second."""
    path = cmd.exp_path(text)
    if not os.path.isfile(path):
        return [(text, "")]
    rows = []
    with open(path, encoding="utf-8-sig") as handle:
        for line in handle:
            # a SMILES has no space, comma or semicolon: any of them is a
            # column separator, also in CSV files
            parts = [p.strip("\"'") for p in re.split(r"[\s,;]+", line.strip())]
            if parts[0] and not parts[0].startswith("#"):
                rows.append((parts[0], parts[1] if len(parts) > 1 else ""))
    if rows and rows[0][0].lower() == "smiles":  # header
        del rows[0]
    return rows


def _embed(text, conformers):
    """MOL blocks of the 3D conformers of one SMILES; [] if rdkit cannot."""
    Chem, AllChem = rdkit()
    mol = Chem.MolFromSmiles(text)
    if mol is None:
        return []
    mol = Chem.AddHs(mol)
    params = AllChem.ETKDGv3()
    params.randomSeed = 42  # the same input gives the same coordinates
    params.pruneRmsThresh = 0.5
    ids = list(AllChem.EmbedMultipleConfs(mol, conformers, params))
    if not ids:  # macrocycles and strained rings need random start coordinates
        params.useRandomCoords = True
        ids = list(AllChem.EmbedMultipleConfs(mol, conformers, params))
    if not ids:
        return []
    if AllChem.MMFFHasAllMoleculeParams(mol):
        AllChem.MMFFOptimizeMoleculeConfs(mol, maxIters=500)
    else:
        AllChem.UFFOptimizeMoleculeConfs(mol, maxIters=500)
    return [Chem.MolToMolBlock(mol, confId=i) for i in ids]


def load_smiles(smiles, name="", conformers=1, group="", quiet=1, *, _self=cmd):
    """
    DESCRIPTION

        "load_smiles" builds a molecule from a SMILES string, with
        hydrogens and minimized 3D coordinates. Given a file, it builds one
        molecule per line.

    USAGE

        load_smiles smiles [, name [, conformers [, group ]]]

    ARGUMENTS

        smiles = string: a SMILES in quotes, or a text or CSV file with the
        SMILES in the first column and, if present, a name in the second

        name = string: object name; for a file, the name of the lines that
        have none {default: mol01, mol02, ...}

        conformers = integer: how many conformers to try for, loaded as
        states of the object; similar ones are dropped {default: 1}

        group = string: group to put the objects in {default: none}

    EXAMPLE

        load_smiles "CC(=O)Oc1ccccc1C(=O)O", aspirin
        load_smiles "c1ccccc1O", phenol, 20
        load_smiles ligands.smi, group=ligands

    NOTES

        Needs rdkit, which the downloadable builds include; for a pip
        install: pip install 'pymolish[chem]'. Coordinates come from
        ETKDGv3 and an MMFF94 minimization (UFF where MMFF94 has no
        parameters). A line that is not a SMILES is reported and skipped.
        "split_states" turns conformers into objects.

    SEE ALSO

        fab, fragment, split_states
    """
    rdkit()
    entries = _entries(smiles.strip().strip("\"'"))
    known = set(_self.get_names("objects"))
    loaded = []
    for text, label in entries:
        blocks = _embed(text, max(1, int(conformers)))
        if not blocks:
            print(" load_smiles: skipped, no 3D structure for '%s'" % text)
            continue
        wanted = _self.get_legal_name(label or name)
        obj = _self.get_unused_name(wanted or "mol", 0 if wanted else 1)
        for state, block in enumerate(blocks, 1):
            _self.read_molstr(block, obj, state, discrete=0, zoom=0)
        # ask for the name: "b" is loaded as "b_", "b" is a reserved word
        new = [n for n in _self.get_names("objects") if n not in known]
        known.update(new)
        loaded += new
    if not loaded:
        raise pymol.CmdException("no molecule built from '%s'" % smiles)
    if group:
        _self.group(group, " ".join(loaded))
    if _self.get_setting_boolean("auto_zoom"):
        _self.zoom(" or ".join("%" + n for n in loaded))
    if not int(quiet):
        print(" load_smiles: built %d of %d molecules" % (len(loaded), len(entries)))
    return loaded
