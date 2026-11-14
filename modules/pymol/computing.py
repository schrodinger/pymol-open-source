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


def _openmm():
    '''(PDBFixer class, openmm.app module), or a CmdException that names
    the extra.'''
    try:
        import openmm.app
        from pdbfixer import PDBFixer
    except ImportError:
        raise pymol.CmdException(
            "this command needs pdbfixer and openmm: "
            "pip install 'pymolish[fix]'") from None
    return PDBFixer, openmm.app


def _new_object(pdbstr, name, _self):
    '''Load PDB text into a new object; a name in use gets a number.'''
    known = set(_self.get_names('objects'))
    _self.read_pdbstr(pdbstr, _self.get_unused_name(name, 0), zoom=0)
    # ask for the name: a reserved word like "all" is loaded as "all_"
    return [n for n in _self.get_names('objects') if n not in known][0]


def fix_structure(selection, missing_atoms=1, hydrogens=1, ph=7.0,
                  nonstandard=1, keep_water=1, keep_heterogens=1, name='',
                  async_=0, quiet=1, *, _self=cmd_module):
    '''
DESCRIPTION

    "fix_structure" repairs a protein or nucleic acid with PDBFixer and
    loads the result as a new object. The original is not changed.

USAGE

    fix_structure selection [, missing_atoms [, hydrogens [, ph
        [, nonstandard [, keep_water [, keep_heterogens [, name ]]]]]]]

ARGUMENTS

    selection = string: atoms of one object

    missing_atoms = 0/1: build heavy atoms that are missing in residues
    that are present, also the terminal oxygen {default: 1}

    hydrogens = 0/1: add hydrogens {default: 1}

    ph = float: pH for the protonation of standard residues {default: 7.0}

    nonstandard = 0/1: replace nonstandard residues with the standard one
    they come from, for example MSE with MET {default: 1}

    keep_water = 0/1 {default: 1}

    keep_heterogens = 0/1: keep ligands and ions {default: 1}

    name = string: name of the new object {default: <object>_fixed}

    async_ = 0/1: work in the background {default: 0}

EXAMPLE

    fix_structure 1abc
    fix_structure 1abc and chain A, ph=5.5, keep_water=0

NOTES

    Needs pdbfixer and openmm: pip install 'pymolish[fix]'; the
    downloadable builds include them.

    Residues that are missing as a whole (gaps in the chain) are not
    built: PDBFixer finds them from SEQRES records, which a loaded object
    does not keep.

    For the hydrogens of a ligand PDBFixer downloads its definition from
    files.rcsb.org; without a connection the ligand gets none. The result
    comes back as PDB text: bond orders of ligands are single, colors and
    representations are the defaults.

SEE ALSO

    minimize, clean, h_add
    '''
    import io

    if int(async_):
        return _self.async_(
            fix_structure, selection, missing_atoms, hydrogens, ph,
            nonstandard, keep_water, keep_heterogens, name, 0, quiet,
            _self=_self)
    PDBFixer, app = _openmm()
    objects = _self.get_object_list(selection)
    if len(objects) != 1:
        raise pymol.CmdException('fix_structure: give atoms of one object')
    if not int(keep_water):
        selection = '(%s) and not solvent' % selection
    before = _self.count_atoms(selection)
    try:
        fixer = PDBFixer(pdbfile=io.StringIO(_self.get_pdbstr(selection)))
        fixer.findMissingResidues()  # none without SEQRES; the later steps ask
        if int(nonstandard):
            fixer.findNonstandardResidues()
            fixer.replaceNonstandardResidues()
        if not int(keep_heterogens):
            fixer.removeHeterogens(keepWater=True)
        if int(missing_atoms):
            fixer.findMissingAtoms()
            fixer.addMissingAtoms(seed=0)
        if int(hydrogens):
            fixer.addMissingHydrogens(float(ph))
        out = io.StringIO()
        app.PDBFile.writeFile(
            fixer.topology, fixer.positions, out, keepIds=True)
    except pymol.CmdException:
        raise
    except Exception as error:  # PDBFixer and OpenMM raise many types
        raise pymol.CmdException('fix_structure: %s' % error) from None
    name = _new_object(out.getvalue(), name or objects[0] + '_fixed', _self)
    if not int(quiet):
        print(' fix_structure: "%s" has %d atoms, the selection had %d'
              % (name, _self.count_atoms('%' + name), before))
    return name
