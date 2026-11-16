# pymolish-studio — reviewed roadmap

Reviewed 2026-11-08 against revision `3250e7a3`, the working tree, and primary
documentation linked below. This is a draft implementation plan: source inspection
and UI research are complete; builds, benchmarks, and prototypes have not run.

**Keep Qt Widgets with PySide6, simplify the supported product, measure performance,
modernize the existing window, and deliver downloadable builds early.** Changing
bindings does not move molecular computation to the GPU. Follow
[AGENTS.md](AGENTS.md): smaller, faster, simpler, prettier.

Finished steps are in [PLAN-DONE.md](PLAN-DONE.md). A step moves there when its
last milestone ships; its detailed plan stays in `docs/plans/`.

## Open steps

| Step | Detailed plan | Recommendation |
| --- | --- | --- |
| 1 | [Performance and GPU](docs/plans/01-performance-gpu.md) | Baseline first; optimize measured CPU, rendering, and interaction bottlenecks. |
| 3 | [Modern workspace and metadata](docs/plans/03-modern-ui.md) | Improve the existing window; ship docks and an object metadata table first. |
| 8 | [Modern C++ in the core](docs/plans/08-cpp-modernization.md) | No rewrite. One idiom per PR behind sanitizers and tests; vendored C stays; files go in profile order so cleanup and step 1 meet. |
| 12 | [pymolish plugin commands, native](docs/plans/12-native-commands.md) | Port command bodies, not the framework: 58 commands become about 17, two dialogs, no new top menu. rdkit approved. |
| 13 | [Clean, fix and optimize structures](docs/plans/13-structure-cleanup.md) | Four layers: remove water and ions, ligand `clean` with rdkit, PDBFixer repair, OpenMM relax. PDBFixer and OpenMM approved. |
| 14 | [Action and right-click menus](docs/plans/14-menus.md) | Red squares are a rendering artifact: one-line fix now. Merge duplicated menus after 12 and 13, from a list the user approves. |

## Validation of the original draft

| Original idea or claim | Verdict and correction |
| --- | --- |
| PySide6 is the best fit | **Keep.** Qt 6 is already supported. PyQt6 is another binding to the same toolkit, not a different renderer. See step 0 for the sourced comparison. |
| Qt already provides docking | **Confirmed, partly implemented.** `PyMOLQtGUI` enables nested/tabbed docks. Workspace layout persistence was not found in the inspected GUI code. |
| Require a Core profile during binding cleanup | **Reject as an early change.** `Scene.cpp` and `Ortho.cpp` still use fixed-function calls. Audit and migrate required rendering paths separately. |
| Remove Tk, its adapters, GLUT and other bindings in one PR | **Split.** Qt plugin initialization still uses legacy helpers, and `legacysupport.py` imports `pmg_tk.startup`. Preserve needed functions before removing their implementation. |
| Enable OpenMP by default | **Already largely true.** `setup.py` defaults it on outside macOS. Marching cubes uses it; `Ray.cpp` has its own workers. Check races and oversubscription before adding parallelism. |
| Move PDB/mmCIF parsing from Python to C++ | **Correct the premise.** Native readers already exist in `ObjectMolecule.cpp`, `ObjectMolecule2.cpp`, and `CifMoleculeReader.cpp`. Profile the actual load path. |
| GPU ray tracing and marching cubes are next | **Defer.** OpenGL 4.3 compute cannot be the common requirement for native macOS OpenGL. Step 1 gives the platform constraint and experiment gate. |
| Add five performance scripts from scratch | **Reuse first.** Existing loading, shader, iteration and asynchronous-build checks provide a starting point. Some are excluded from `--run all`; headless tests do not establish GPU performance. |
| Stop updating `ChangeLog` | **Reject.** Public `cmd.*` behavior changes must still update it. Preserve historical entries. |
| Optional libraries are immediately removable | **Conditional.** COLLADA, MMTF, trajectories, and glTF have real dependencies. CI explicitly enables `json=True`; retained native glTF needs nlohmann-json. |
| PNG or glTF replaces movie encoding | **Incorrect.** Image sequences, encoded movies, and 3D scenes are different outputs. Preserve an explicit supported movie workflow. |
| Replace the whole window beside the old one | **Prefer incremental changes.** Two windows duplicate command, dialog, plugin, and lifecycle integration. Extract only the panel being implemented. |
| Object and sequence panels can simply become docks | **Underestimated.** Existing internal controls involve C++ drawing; a Qt replacement needs separate behavior and performance checks. |
| `cmd.get_coords()` is a NumPy view | **Correct.** It returns selection coordinates as an array; the explicit borrowed-memory option is `get_coordset(copy=0)`, with documented lifetime hazards. |
| `cmd.add_dock()` is a 20-line task | **Defer the API.** Ownership, thread access, headless behavior, and layout restoration matter. Use existing window access first. |
| Cross-platform binaries are a small mechanical task | **Underestimated.** Resources, shared libraries, Qt plugins, architecture, and clean-machine testing need milestones. |
| Test PRs only on Linux/Python 3.12 | **Reject for this migration.** Preserve three-OS and Python 3.11–3.14 coverage; optimize caching and artifact jobs first. |
| Latest pins prove compatibility; LTS defines the runtime floor | **Not established.** The PySide6 6.12.0 pin is real, but a minimum version needs tests. LTS does not imply unrestricted access to all maintenance packages. |
| Replace GLEW, collapse directories, cap setup.py at 300 lines | **Skip without evidence.** These can add coupling and churn without improving a measured workflow. |

Evidence and implementation boundaries are linked in each step. “Keep” means a
recommended direction, not a claim that an implementation has passed tests.

## Open work, in execution order

The topic numbers remain the user's original numbers. The letters identify
milestones in the detailed plans, not new topics.

| Order | Work | Exit condition | Relative size |
| --- | --- | --- | --- |
| 5 | **1b — measured optimization** and **3 — workspace increments** | Performance changes have before/after evidence; docks and table ship in reviewable test builds. | 3 slices 1-3 done 2026-11-08 (4.5.0: objects dock, layout persistence, Panels menu); slice 4 done (4.5.11 Qt Objects panel replaces the viewport's object list; 4.5.14 workspace layout: action toolbar on top, command area at the bottom, toggle bar); slice 5 partly done (4.5.2 dark theme, 4.5.7 UI Scale for high-resolution screens and dark-theme color audit); 1b started (4.5.10 session load 6.6x faster, 4.5.12 surface build 1.24x, 4.7.7 to 4.7.9 surface build another 1.5x; GUI frame time measured 2026-11-14, software rendering under WSL is the main cost); renderer default on WSL, sticks and surface draw cost, and slice 6 open |
| 9 | **14a — red squares** | Destructive menu rows show a red trash icon. | Small. No inputs, no dependencies: first. 10a and 10b, the rest of orders 9 and 10, are done (4.9.0). |
| 12 | **12a to 12e — native commands** | Kept plugin commands work without the plugin, from the command line and from File, Build and the group menu. | Large; five PRs. |
| 13 | **13a to 13d — cleanup** | Water and ions go in one click; `clean` works on ligands; a protein is repaired and relaxed into a new object. | Medium to large. 13a needs nothing; 13b needs 12d. |
| 14 | **14b, 14c — menu consolidation** | Approved keep, merge, drop list applied; fewer rows and fewer lines in `menu.py`. | Medium. Last, so the rows from 12 and 13 are audited once. |

Orders 9 to 14 were requested on 2026-11-14 and do not block order 5.

### Working rules for orders 9 to 14

- One slice per session and per PR. Each plan names the files and lines to
  open, so a session starts there instead of exploring again.
- [ponytail](https://github.com/dietrichgebert/ponytail) is already the rule
  (`AGENTS.md` §1): each plan has a "Skipped" list and, where code goes away,
  a "Deleted" list. Nothing on a "Skipped" list is built without a new request.
- [caveman](https://github.com/JuliusBrussee/caveman) shortens replies only;
  optional. The larger saving is the narrow slice and a fresh session per slice.
- Decided on 2026-11-14: rdkit, PDBFixer and OpenMM are approved optional
  dependencies; the palette is remembered in `QSettings`; `LICENSE` gets a
  header above the unchanged notice. Still the user's: the license line for
  the two images (9, 10), which dropped plugin commands to keep (12), the
  menu list (14b).

Cleanup precedes speculative optimization. Recording a baseline is preparation.
UI work needs the binding decision and a usable test build, not completion of
every unrelated removal. GPU compute is not on the critical path. Plugin
compatibility work starts before Tk removal; final extension documentation can
follow the first UI increments. Artifact work can start after the binding change
while the remaining legacy cleanup proceeds.

## Shared completion criteria

- Preserve read/write PDB, mmCIF, SDF/MOL, PSE; PNG export; basic selection,
  representations, command entry, and supported plugin loading.
- Each non-trivial implementation change leaves a runnable check in
  `testing/tests/` or `layerCTest/`; removals run the full Python suite.
- Rendering/UI changes need real GUI and platform checks. Report skipped tests.
- Keep CI coverage, document public API changes in `ChangeLog`, and preserve
  `LICENSE`, `AUTHORS`, and copyright headers.
- New dependencies, data-file deletions, and file-format output changes follow
  the decisions required by `AGENTS.md` during implementation. These plans do not
  install, delete, or change them.
- Measure runtime, memory, artifact size, and build time before fixing hard
  targets. No speedup, binary size, or delivery date is promised by this draft.
