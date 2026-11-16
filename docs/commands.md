# Commands added by pymolish-studio

Commands that PyMOL does not have. Each one works on the command line, in
scripts as `cmd.<name>`, and has `help <name>`. Since version 4.11.0 they
replace the commands of the pymolish plugin; remove the plugin from
`~/.pymol/startup/`, or its commands of the same name win.

## Groups

| Command | What it does | Example |
| --- | --- | --- |
| `group_next [group]` | Shows the next member of a group and hides the others. Key: Ctrl-N. | `group_next poses` |
| `group_previous [group]` | The same, backward. Key: Ctrl-P. | `group_previous` |
| `group_enable group [, each]` | Shows the members of a group, or every n-th one. | `group_enable poses, 10` |
| `group_disable group [, each]` | Hides the members; the group stays enabled. | `group_disable poses` |
| `copy_group group [, new_group [, suffix]]` | Copies a group with all members, groups inside it included. | `copy_group poses, backup, _bak` |
| `sort_group group [, reverse]` | Sorts the members by name, ignoring case. | `sort_group poses` |
| `extract_chain selection [, chains [, group]]` | Makes one new object per chain, `<object>_<chain>`. The source is not changed. | `extract_chain poses, A+B, poses_AB` |

| `cluster_ligands selection [, count [, cutoff [, prefix]]]` | Sorts objects into groups by where they are: by a distance cutoff, or into a given number of groups. | `cluster_ligands poses` |
| `color_group_with_gradient selection [, color]` | Gives every object its own shade of one color, dark to light. | `color_group_with_gradient poses, orange` |

The group row in the Objects panel has four of them in its A menu:
duplicate group, sort members (a second click sorts Z to A), enable members,
disable members.

## Many files at once

| Command | What it does | Example |
| --- | --- | --- |
| `load_files [directory [, format [, group [, pattern [, recursive [, prefix]]]]]]` | Loads every file of one format from a folder, one object per file. Unreadable files are reported and skipped; no object is replaced. | `load_files ~/poses, sdf, poses` |
| `export_objects selection [, directory [, format [, overwrite]]]` | Saves every object of a selection (a group, a name pattern, any selection) to `<object>.<format>`. Files that exist are skipped unless `overwrite=1`. | `export_objects poses, ~/out, pdb` |

In the window: File > Open Folder... and File > Export Objects.... Both show
the command they will run.

## Structures from the internet

| Command | What it does | Example |
| --- | --- | --- |
| `list_uniprot_structures accession` | Prints the PDB entries of a protein with method, resolution and chains. Loads nothing. | `list_uniprot_structures P04637` |
| `fetch_by_uniprot accession [, limit [, group [, async_]]]` | Downloads the PDB entries of a protein into one group. Give a limit: a well studied protein has hundreds. | `fetch_by_uniprot P00918, 5` |
| `fetch_uniprot_ligands accession [, min_weight [, limit [, group [, async_]]]]` | Loads every different small molecule found in the PDB entries of a protein into the group `<accession>_ligands`. Each is the ideal model, alone, not in its binding place. Ions and additives lighter than `min_weight` (100 Da) are left out. | `fetch_uniprot_ligands P00918, limit=20` |
| `fetch_fold query [, name [, color [, async_]]]` | Loads a predicted structure: from AlphaFold DB for a UniProt accession or an AlphaFold id, from the ESMFold server for a sequence. | `fetch_fold P04637` |
| `color_plddt [selection [, legend]]` | Colors a predicted model by confidence with the four AlphaFold DB colors and adds a legend. | `color_plddt AF_P04637_F1` |

In the window: File > Get by UniProt... and File > Get Predicted
Structure.... These commands ask rest.uniprot.org, data.rcsb.org,
files.rcsb.org, alphafold.ebi.ac.uk and api.esmatlas.com; a sequence given to `fetch_fold` is sent to the ESMFold
server.

## Molecules from SMILES

| Command | What it does | Example |
| --- | --- | --- |
| `load_smiles smiles [, name [, conformers [, group]]]` | Builds a molecule from a SMILES in quotes, with hydrogens and minimized 3D coordinates; conformers become states. A text or CSV file gives one molecule per line (SMILES, then an optional name). | `load_smiles "CC(=O)Oc1ccccc1C(=O)O", aspirin` |

In the window: Build > From SMILES..., filled from the clipboard. Needs
rdkit: the downloadable builds have it, a pip install needs
`pip install 'pymolish[chem]'`.

## Sequence

| Command | What it does | Example |
| --- | --- | --- |
| `search_sequence pattern [, name [, selection]]` | Selects the residues of every protein chain that match a sequence or a regular expression. | `search_sequence N[^P][ST], sequon` |

## Repair

| Command | What it does | Example |
| --- | --- | --- |
| `fix_structure selection [, missing_atoms [, hydrogens [, ph [, nonstandard [, keep_water [, keep_heterogens [, name]]]]]]]` | Repairs a protein with PDBFixer and loads the result as `<object>_fixed`: missing atoms, hydrogens for a pH, nonstandard residues. The original stays. | `fix_structure 1abc, ph=5.5` |

| `minimize selection [, maxiter [, name]]` | Relaxes a protein with OpenMM (Amber14, implicit solvent), heavy atoms held near their place; the result is `<object>_min`. Needs complete residues and no ligands. | `minimize prot` |

`clean` is a PyMOL command that works here: it minimizes one small molecule
with rdkit. In the window: A > clean on an object, a group or a selection.
`fix_structure` and `minimize` need pdbfixer and openmm: the downloadable builds have
them, a pip install needs `pip install 'pymolish[fix]'`.

## Other

| Command | What it does | Example |
| --- | --- | --- |
| `palette [name [, recolor]]` | Sets the colors that new objects cycle through. | `palette vivid` |
