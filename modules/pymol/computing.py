#A* -------------------------------------------------------------------
#B* This file contains source code for the PyMOL computer program
#C* Copyright (c) Schrodinger, LLC
#D* -------------------------------------------------------------------
#E* It is unlawful to modify or remove this copyright notice.
#F* -------------------------------------------------------------------
#G* Please see the accompanying LICENSE file for further information.
#H* -------------------------------------------------------------------
#I* Additional authors of this source file include:
#-*
#-*
#-*
#Z* -------------------------------------------------------------------

import sys
cmd_module = __import__("sys").modules["pymol.cmd"]

import pymol

# shortcut: no pocket and no ligand parameters beyond MMFF94; above this many
# atoms the time and the missing protein terms make the result useless
CLEAN_MAX_ATOMS = 1000


def clean(selection, present='', state=-1, fix='', restrain='',
          method='mmff', async_=0, save_undo=1, message=None,
          _self=cmd_module, **kwargs):
    '''
DESCRIPTION

    "clean" minimizes the geometry of a small molecule with the MMFF94
    force field of rdkit: bond lengths, angles and clashes inside the
    molecule are corrected, atoms stay close to where they were.

USAGE

    clean selection [, state ]

ARGUMENTS

    selection = string: one whole molecule of one object, up to 1000 atoms

    state = integer: state to change {default: -1, the current state}

NOTES

    Needs rdkit: pip install 'pymolish[chem]'; the downloadable builds
    include it.

    Missing hydrogens are added for the calculation only. Bond orders are
    used as they are: a ligand from a PDB file has only single bonds, so
    set its double bonds first (Build menu, or "valence"), otherwise
    planar rings pucker.

    The atoms around the molecule are not seen: a ligand can move into
    the protein. For a protein use "minimize".

    "present", "fix" and "restrain" are not supported and raise an error;
    the only method is mmff. "save_undo" and "message" are accepted and
    not used.

SEE ALSO

    minimize, fix_structure, sculpt_activate, h_add
    '''
    from .smiles import rdkit

    if present or fix or restrain:
        raise pymol.CmdException(
            'clean: "present", "fix" and "restrain" are not supported')
    if method != 'mmff':
        raise pymol.CmdException('clean: the only method is mmff')
    if int(kwargs.pop('async', async_)):
        return _self.async_(clean, selection, state=state, _self=_self)
    Chem, AllChem = rdkit()
    state = int(state) if int(state) > 0 else _self.get_state()
    count = _self.count_atoms(selection)
    if len(_self.get_object_list(selection)) != 1:
        raise pymol.CmdException('clean: give atoms of exactly one object')
    if count != _self.count_atoms('bymolecule (%s)' % selection):
        raise pymol.CmdException(
            'clean: the selection is a part of a molecule; give all of it')
    if count > CLEAN_MAX_ATOMS:
        raise pymol.CmdException(
            'clean is for small molecules: %d atoms, the limit is %d. '
            'For a protein use "minimize".' % (count, CLEAN_MAX_ATOMS))
    mol = Chem.MolFromMolBlock(
        _self.get_str('mol', selection, state), removeHs=False)
    if mol is None:
        raise pymol.CmdException(
            'clean: rdkit cannot read the molecule; check valences and '
            'formal charges')
    if all(a.GetAtomicNum() > 1 for a in mol.GetAtoms()) and all(
            b.GetBondTypeAsDouble() == 1 for b in mol.GetBonds()):
        print(' clean: no hydrogens and only single bonds. If double '
              'bonds are missing, set them and run clean again.')
    mol = Chem.AddHs(mol, addCoords=True)  # appended after the atoms we have
    if not AllChem.MMFFHasAllMoleculeParams(mol):
        raise pymol.CmdException(
            'clean: MMFF94 has no parameters for an atom of this molecule')
    AllChem.MMFFOptimizeMolecule(mol, maxIters=1000)
    _self.load_coords(
        mol.GetConformer().GetPositions()[:count], selection, state)
