# 12 — pymolish plugin commands, native

[Roadmap](../../PLAN.md) · Request 2026-11-14 · **Done 2026-11-14, 4.11.0 to
4.11.4.** 18 commands (19 since 4.12.9: `fetch_uniprot_ligands` came back on
the user's request of 2026-11-16) in five new files (`grouping.py`, `batch.py`,
`fetching.py`, `smiles.py`, `finding.py`), listed in
[docs/commands.md](../commands.md). The text below is the plan as written
before the work (plugin 0.1.27, `b970e67`).

**Port the useful command bodies of the
[pymolish plugin](https://github.com/adrian-jasinski/pymolish) into the
application, with menu entries and two dialogs. Leave the plugin's framework
and its meta commands behind: 58 commands become about 17.**

## What exists

- **Plugin:** MIT, same author. 7,600 lines of source, 58 commands.

  | Part | Lines | Native equivalent |
  | --- | --- | --- |
  | `core/` registry, completion, logging, validators | 700 | `cmd.extend`, `help`, `cmd.auto_arg`, `CmdException` |
  | `core/deps.py`, `core/upgrade.py` | 370 | none needed once built in |
  | `__init__.py` bootstrap, `gui/` | 480 | `pymol.Qt`, `.ui` forms |
  | `extensions/` | 6,100 | port, see triage |

- **Its hooks:** wraps `pymol.menu.group_action`; binds Ctrl-N and Ctrl-P.
- **Its optional dependencies:** numpy, matplotlib, biopython, rdkit, PyQt5.
- **Already native here:** File > Open takes several files
  (`pymol_qt_gui.py:662`); `extra_fit` aligns many objects
  (`modules/pymol/fitting.py:203`); UI scale; selectors `solvent`,
  `inorganic`, `metals`, `pepseq`; dialog pattern with a live command preview
  in `modules/pmg_qt/file_dialogs.py`.

## Decisions

- **Port bodies, not the framework.** New files, so upstream merges do not
  conflict and the MIT boundary is clear: `modules/pymol/grouping.py`,
  `batch.py`, `fetching.py`. Registered like every command family (`api.py`,
  `keywords.py`), so `cmd.<name>`, `help <name>` and tab completion work.
  Docstrings rewritten in PyMOL's DESCRIPTION / USAGE form.
- **No new top-level menu.**

### Command triage

| Family | Keep | Drop, and why |
| --- | --- | --- |
| Groups, 18 → 7 | `group_next`, `group_previous`, `group_enable`, `group_disable`, `copy_group`, `sort_group`, `extract_chain` (one command for one or all chains) | `group_info`, `group_status`, `list_*`: listings; `group_toggle`, `sort_all_groups`, `create_sorted_group`, `copy_group_objects`, `merge_to_*`: variants. Say which you use and it stays. |
| I/O, 9 → 3 | `load_files` (with `recursive=`), `export_objects` (group or name pattern), `load_smiles` | `load_recursive`, `list_loadable_files`, `export_group`, `export_by_pattern`: folded in; `paste_smiles`: becomes a menu entry; `export_byres_to_csv`: niche |
| Selection, 3 → 1 | `search_sequence`, rebuilt on `pepseq` and `re` | `water`, `ions`: selectors `solvent` and `inorganic` exist |
| Structure, 10 → 4 | `fetch_by_uniprot`, `fetch_fold`, `cluster_ligands`, `list_uniprot_structures` | `multialign`: `extra_fit` (wrap only if it cannot take a group); `list_objects_with_prefix`: `cmd.get_names`; `remove_ions`, `remove_metals`: [step 13](13-structure-cleanup.md); `fetch_uniprot_ligands`: brought back in 4.12.9 |
| Visualization, 9 → 2 | `color_plddt`, `color_group_with_gradient` | palettes: [step 11](11-color-palettes.md); `color_secondary_structure`: C > by ss exists; `ui_scale`: toolbar |
| Meta, 9 → 0 | | `help` is the help; version is `get_version`; logo is the About box; no upgrade or installer |

### Dependencies

| Package | Used by | Decision |
| --- | --- | --- |
| numpy | `cluster_ligands` | already a runtime dependency |
| matplotlib | declared, imported nowhere | not added |
| biopython | `search_sequence`, UniProt metadata | not added; `pepseq`, `re`, `urllib`, `json` do it. Ceiling: no fuzzy sequence match. |
| PyQt5 | plugin GUI | not added; PySide6 through `pymol.Qt` |
| rdkit | `load_smiles` (SMILES to 3D needs a real embedder) and ligand clean in step 13 | **add**: BSD-3, 30 MB wheel, pulls Pillow |

**rdkit: approved by the user on 2026-11-14**, recorded in `AGENTS.md` §2.2.

- Release bundles include it: that is the build users download.
- For `pip install` it is the extra `chem`, imported lazily; without it the
  command raises `CmdException` naming the extra.
- Not a hard dependency: that would cost every install 30 MB and break on a
  Python version rdkit has no wheel for. Check 3.14 wheels first; a missing
  wheel drops rdkit from that Python's bundle, it does not block the release.
- Record bundle size before and after in `docs/baseline.md`.

### GUI

| Where | Entry | Opens |
| --- | --- | --- |
| File | Open Folder... | `load_folder.ui`: directory picker, format, recursive, name filter, group name, live count of matching files, command preview |
| File | Export Objects... | small form: group or pattern, directory, format |
| File | Get by UniProt..., Get Predicted Structure... | one text field each |
| Build | From SMILES... | one text field, pre-filled from the clipboard |
| Group row, A menu | duplicate, sort, enable members, disable members | written into `menu.group_action`, no wrapper |
| Keys | Ctrl-N, Ctrl-P | next and previous member; both are free in `modules/pymol/shortcut_dict.py` |

Every dialog runs the command through `cmd.do`, so it lands in the log.

## Slices, one PR each

1. **12a groups:** `grouping.py`, group menu rows, keys. Done in 4.11.0.
2. **12b batch I/O:** `batch.py`, the two File dialogs. Time 500 small files;
   over one second, move the load loop off the GUI thread. Done in 4.11.1:
   0.4 s, the loop stays on the GUI thread.
3. **12c fetch:** `fetching.py`, `color_plddt`, two File entries. Network tests
   stay out of `--run all`. Done in 4.11.2.
4. **12d rdkit:** `load_smiles`, Build entry, extra `chem`, spec, CI, licenses row. Done in 4.11.3.
5. **12e:** `search_sequence`, `cluster_ligands`. Done in 4.11.4, with
   `color_group_with_gradient`, which no slice had named.

Each slice: one real test per kept command in `testing/tests/api/` (the
plugin's mock-based tests are not ported), a `ChangeLog` line, a row in
`docs/commands.md` (hand-written table: command, what, example).

## Deleted or not carried over

- Plugin framework, meta commands, installer, upgrader, docs generator:
  about 1,500 lines never enter this repository.
- 41 commands, per the triage.
- Residue-name lists for ions and metals.

## Skipped

- Guard against the plugin also being installed in `~/.pymol/startup` (it
  would shadow the native names): document "remove the plugin"; add a notice
  only if someone hits it.
- Fuzzy sequence search, UniProt metadata through biopython.
- The plugin repository itself: unchanged, it still serves upstream PyMOL.
