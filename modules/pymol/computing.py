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

    class Fixer(PDBFixer):
        # ponytail: PDBFixer 1.12 stops on a ligand definition without
        # ideal coordinates ("?", for example PID). Such a ligand gets no
        # hydrogens; delete this class when PDBFixer reads these files.
        def _downloadCCDDefinition(self, residueName):
            try:
                return super()._downloadCCDDefinition(residueName)
            except ValueError:
                print(' fix_structure: no hydrogens for %s, its definition '
                      'at files.rcsb.org has no ideal coordinates'
                      % residueName)
                self._ccdCache[residueName.upper()] = None

    return Fixer, openmm.app


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

    hydrogens = 0/1: build all hydrogens for the pH; the ones that are
    there are replaced {default: 1}

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
    files.rcsb.org; without a connection, or when the definition has no
    ideal coordinates, the ligand gets none. The result
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
    if int(hydrogens):
        # all of them are built for the pH; hydrogens from "h_add" have
        # other names and PDBFixer would add a second set next to them
        selection = '(%s) and not hydro' % selection
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


# heavy atoms are held near their start: 10 kcal/mol/A^2, in kJ/mol/nm^2
MINIMIZE_RESTRAINT = 4184.0


def minimize(selection, maxiter=500, name='', async_=0, quiet=1, *,
             _self=cmd_module):
    '''
DESCRIPTION

    "minimize" relaxes a protein with OpenMM: Amber14 force field, implicit
    solvent (GBn2, 10 A cutoff), heavy atoms held near where they are. Clashes and
    strained geometry go away, the fold does not move. The result is a new
    object; the original is not changed.

USAGE

    minimize selection [, maxiter [, name ]]

ARGUMENTS

    selection = string: atoms of one object: protein or nucleic acid
    with all heavy atoms, water and simple ions

    maxiter = integer: most minimization steps, 0 for until converged
    {default: 500}

    name = string: name of the new object {default: <object>_min}

    async_ = 0/1: work in the background {default: 0}

EXAMPLE

    fix_structure 1abc and polymer, name=prot
    minimize prot

NOTES

    Needs openmm: pip install 'pymolish[fix]'; the downloadable builds
    include it. The CPU platform is used.

    The hydrogens of "fix_structure" are kept. When hydrogens are missing
    or have names the force field does not know (from "h_add" or "fab"),
    all of them are built again for pH 7. Missing heavy atoms, also the
    terminal oxygen, stop the command: run "fix_structure" first. There are no
    parameters for ligands and modified residues: remove them or give
    "polymer". A chain with a gap is refused, because the force field
    would bond the two ends.

    Prints the energy before and after in kJ/mol and the RMSD of the
    heavy atoms.

SEE ALSO

    fix_structure, clean, sculpt_activate
    '''
    import io

    if int(async_):
        return _self.async_(
            minimize, selection, maxiter, name, 0, quiet, _self=_self)
    _, app = _openmm()
    import numpy
    import openmm
    from openmm import unit

    objects = _self.get_object_list(selection)
    if len(objects) != 1:
        raise pymol.CmdException('minimize: give atoms of one object')
    try:
        pdb = app.PDBFile(io.StringIO(_self.get_pdbstr(selection)))
        forcefield = app.ForceField(
            'amber14-all.xml', 'amber14/tip3p.xml', 'implicit/gbn2.xml')
        topology, positions = pdb.topology, pdb.positions
        # measured on 260 residues: 6 s with this cutoff, 270 s without,
        # the same heavy-atom RMSD (docs/baseline.md)
        options = dict(nonbondedMethod=app.CutoffNonPeriodic,
                       nonbondedCutoff=1.0 * unit.nanometer)
        try:
            system = forcefield.createSystem(topology, **options)
        except ValueError:
            # hydrogens are missing, or come from "h_add" or "fab" with
            # names the force field does not know: build them all
            modeller = app.Modeller(topology, positions)
            modeller.delete([a for a in topology.atoms()
                             if a.element is not None
                             and a.element.symbol == 'H'])
            modeller.addHydrogens(forcefield)
            topology, positions = modeller.topology, modeller.positions
            system = forcefield.createSystem(topology, **options)
            print(' minimize: hydrogens built for pH 7')
        start = positions.value_in_unit(unit.nanometer)
        for a, b in topology.bonds():
            if a.residue != b.residue and {a.name, b.name} == {'C', 'N'}:
                # 0.13 nm is a peptide bond; poor models reach 0.25
                if numpy.linalg.norm(
                        numpy.subtract(start[a.index], start[b.index])) > 0.3:
                    raise pymol.CmdException(
                        'minimize: chain %s has a gap between residues %s '
                        'and %s; give one piece at a time'
                        % (a.residue.chain.id, a.residue.id, b.residue.id))
        hold = openmm.CustomExternalForce(
            '0.5*k*((x-x0)^2+(y-y0)^2+(z-z0)^2)')
        hold.addGlobalParameter('k', MINIMIZE_RESTRAINT)
        for xyz in 'x0', 'y0', 'z0':
            hold.addPerParticleParameter(xyz)
        heavy = [a.index for a in topology.atoms()
                 if a.element is not None and a.element.symbol != 'H']
        for index in heavy:
            hold.addParticle(index, start[index])
        hold.setForceGroup(1)  # so that the energy is reported without it
        system.addForce(hold)
        context = openmm.Context(
            system, openmm.VerletIntegrator(0.001),
            openmm.Platform.getPlatformByName('CPU'))
        context.setPositions(start)

        def energy():
            state = context.getState(getEnergy=True, groups={0})
            return state.getPotentialEnergy().value_in_unit(
                unit.kilojoule_per_mole)

        before = energy()
        openmm.LocalEnergyMinimizer.minimize(
            context, maxIterations=int(maxiter))
        after = energy()
        end = context.getState(getPositions=True).getPositions(asNumpy=True)
        out = io.StringIO()
        app.PDBFile.writeFile(topology, end, out, keepIds=True)
    except pymol.CmdException:
        raise
    except ValueError as error:  # no template: unknown residue or atoms
        raise pymol.CmdException(
            'minimize: %s\n The force field covers standard residues with '
            'all their atoms, water and simple ions. Run "fix_structure" '
            'first, or leave out ligands with "polymer".' % error) from None
    except Exception as error:  # OpenMM raises many types
        raise pymol.CmdException('minimize: %s' % error) from None
    name = _new_object(out.getvalue(), name or objects[0] + '_min', _self)
    moved = numpy.asarray(end.value_in_unit(unit.angstrom))[heavy]
    moved -= 10 * numpy.asarray(start)[heavy]
    rmsd = float(numpy.sqrt((moved ** 2).sum(axis=1).mean()))
    if not int(quiet):
        print(' minimize: energy %.0f -> %.0f kJ/mol, heavy atoms moved '
              '%.2f A (RMSD), new object "%s"' % (before, after, rmsd, name))
    return name, before, after, rmsd
