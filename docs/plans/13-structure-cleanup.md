# 13 — Clean, fix and optimize structures

[Roadmap](../../PLAN.md) · Request 2026-11-14 · **Done 2026-11-14, 4.12.0 to
4.12.3.** The text below is the plan as written before the work; the notes
at each slice say what changed.

**Four layers, cheapest first: remove water and ions (no dependency), make the
dead `clean` work for ligands (rdkit), repair proteins (PDBFixer), relax
proteins (OpenMM). Each layer ships alone.**

## What exists

- **`clean` is dead.** `cmd.clean` raises `IncentiveOnlyException`
  (`modules/pymol/computing.py:20`), yet A > clean is in the object and
  selection menus (`modules/pymol/menu.py`, `mol_action`, `sele_action`).
- **Removal today:** A > remove waters. Nothing for ions. Selectors `solvent`,
  `inorganic`, `organic`, `polymer`, `metals` exist.
- **Hydrogens:** `h_add`, `remove hydro`, native.
- **Build > Sculpting:** native, but geometry regularization, not a force field.
- **No code rebuilds missing atoms or residues.**
- **`chempy/bmin`, `chempy/tinker`:** wrappers for external programs; only
  `experimenting.py:84` and `:116` use tinker. Removal candidates, separate PR.

| Package (PyPI, 2026-11-14) | License | Size, needs |
| --- | --- | --- |
| pdbfixer 1.12.0 | MIT | source only, pure Python; needs numpy, openmm >= 8.2, `legacy-cgi` on Python 3.13+ |
| openmm 8.6.1 | MIT and LGPL parts; confirm from its license files | 14 MB Linux wheel |
| rdkit 2026.9.2 | BSD-3 | arrives with [step 12](12-native-commands.md) |

## Decisions

- **PDBFixer over the alternatives.** Nothing in the tree builds missing heavy
  atoms. pdb2pqr can, but brings no minimizer; OpenMM comes with PDBFixer and
  gives layer 13d for free. Fallback if OpenMM cannot be bundled: pdb2pqr for
  13c, no 13d.
- **Selectors, not residue-name lists.** Ions are `inorganic`. Ceiling:
  carbon-containing additives (GOL, EDO, ACT, PEG) count as ligands and stay;
  add one short `resn` list if users ask.
- **Never overwrite.** Fix and optimize write `<object>_fixed` / `<object>_min`;
  the original stays. Removal edits in place and is undoable like `remove`.
- **PDBFixer and OpenMM: approved by the user on 2026-11-14**, recorded in
  `AGENTS.md` §2.2. Extra `fix`, imported lazily; without it the command raises
  `CmdException` naming the extra. Release bundles include both; record the
  bundle size before and after in `docs/baseline.md`.
- **CPU platform only** in bundles: no CUDA, OpenCL or HIP packages
  (`AGENTS.md` §2.3).
- **Long work off the GUI thread** (`AGENTS.md` §2.3): 13b to 13d run through
  the existing `async_` path; the busy indicator shows; no cancel button.

### A > clean becomes a submenu

| Entry | Runs | Slice |
| --- | --- | --- |
| remove waters | `remove solvent and (sele)` | 13a |
| remove ions | `remove inorganic and (sele)` | 13a |
| remove waters and ions | both | 13a |
| remove ligands | `remove organic and (sele)` | 13a |
| optimize ligand geometry | `clean sele` | 13b |
| fix structure... | dialog, then `fix_structure` | 13c |
| optimize protein geometry | `minimize` | 13d |

The top-level "remove waters" row goes away (it moved): one less duplicate for
[step 14](14-menus.md).

## Slices

**13a, remove. Done in 4.12.0.** Menu rows only; no new command, the selectors are the API.
Test: PDB string with protein, water, Na and a ligand; after each entry the
matching `count_atoms` is 0 and the rest is unchanged.

**13b, `clean` for small molecules. Done in 4.12.1.** Fill the existing `cmd.clean` signature:
selection to MOL (`cmd.get_str`), `MMFFOptimizeMolecule`, coordinates back with
`cmd.load_coords`. Method name stays `mmff`.
- Cap at about 1,000 atoms; above it the message points to 13d.
- `fix` and `restrain` arguments: raise "not supported" rather than ignore.
- Ceiling: the pocket around a ligand is not seen.
- Test: stretched C-C bond in ethane returns to 1.53 Å ± 0.02; skipped without rdkit.
- Public behavior change (`clean` works): `ChangeLog`.

**13c, `fix_structure`. Done in 4.12.2, without `missing_residues`:** PDBFixer
finds missing residues from SEQRES records, and a PyMOL object does not keep
them, so the option could never do anything on this path. The bundle ships
only the Reference and CPU platforms of OpenMM.
`fix_structure object [, missing_atoms=1, hydrogens=1, ph=7.0, nonstandard=1, missing_residues=0, keep_water=1, keep_heterogens=1, name=]`
- Path: `cmd.get_pdbstr` → `PDBFixer(pdbfile=StringIO)` → `PDBFile.writeFile(keepIds=True)` → `cmd.read_pdbstr`.
- Dialog `fix_structure.ui`: those checkboxes and the pH. Missing residues is
  off by default: built loops are unrefined; the dialog says so.
- Test: build `fab AKA`, delete one side-chain atom, fix, atom count is back;
  skipped without pdbfixer.

**13d, `minimize`. Done in 4.12.3**, with a 1.0 nm cutoff (45 times faster,
same result, `docs/baseline.md`); the nonfunctional tinker-based
`cmd.minimize` and `cmd.fast_minimize` gave up the name. OpenMM, Amber14, implicit solvent, heavy atoms restrained
so clashes relax without drift; `maxiter` argument.
- Standard residues, water and ions only. A ligand without parameters gives a
  clear error: remove it or run on the protein alone.
- Ceiling: no ligand force field (openmmforcefields and OpenFF are hundreds of
  MB); revisit only with a measured need.
- Test: energy goes down, heavy-atom RMSD to the start stays under 0.5 Å.
- Record time for a 300-residue protein in the PR.

Every slice adds its section to `docs/features.md` and rows to
`docs/licenses.md`.

## Deleted

- The `IncentiveOnlyException` stub in `clean`.
- The duplicated "remove waters" rows in three action menus.
- Not ported: the plugin's `remove_ions`, `remove_metals`, `water`, `ions`.

## Skipped

- Membrane systems, explicit solvent, dynamics, GPU platforms.
- Ligand protonation and parameters.
- Cancel button, progress percentage.
- Removing `chempy/bmin` and `chempy/tinker`: its own `remove:` PR.
