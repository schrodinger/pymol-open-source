# pymolish-studio — reviewed roadmap

Reviewed 2026-11-08 against revision `3250e7a3`, the working tree, and primary
documentation linked below. This is a draft implementation plan: source inspection
and UI research are complete; builds, benchmarks, and prototypes have not run.

**Keep Qt Widgets with PySide6, simplify the supported product, measure performance,
modernize the existing window, and deliver downloadable builds early.** Changing
bindings does not move molecular computation to the GPU. Follow
[AGENTS.md](AGENTS.md): smaller, faster, simpler, prettier.

## Plans for all six steps

| Step | Detailed plan | Recommendation |
| --- | --- | --- |
| 0 | [UI engine research](docs/plans/00-ui-engine.md) | PySide6 + Qt Widgets; preserve the renderer while simplifying bindings. |
| 1 | [Performance and GPU](docs/plans/01-performance-gpu.md) | Baseline first; optimize measured CPU, rendering, and interaction bottlenecks. |
| 2 | [Code and dependencies](docs/plans/02-cleanup-dependencies.md) | One feature removal per PR; trace dependencies before deleting. |
| 3 | [Modern workspace and metadata](docs/plans/03-modern-ui.md) | Improve the existing window; ship docks and an object metadata table first. |
| 4 | [Features and external libraries](docs/plans/04-extensibility.md) | Preserve `cmd` and the existing loader; document small integration patterns. |
| 5 | [Builds and test binaries](docs/plans/05-builds-binaries.md) | Improve incremental builds and deliver test artifacts before large changes. |

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

## Recommended execution order

The topic numbers remain the user's original numbers. The letters identify
milestones in the detailed plans, not new topics.

| Order | Work | Exit condition | Relative size |
| --- | --- | --- | --- |
| 1 | **0 — UI decision** | PySide6/Widgets recommendation, supported platforms and GL smoke checks are explicit. | Small research task; implemented in working tree 2026-11-08, see step 0 Status |
| 2 | **5a + 1a — feedback and baseline** | Repeatable build/test commands and recorded build/performance timings before removals. | Done 2026-11-08: [docs/baseline.md](docs/baseline.md), `testing/bench.py` |
| 3 | **2a — binding and legacy GUI cleanup**, with **4a — plugin contract** | Explicit runtime binding; alternative bindings removed in one PR; plugin migration precedes separate Tk and GLUT removals. | Done 2026-11-08: 4.1.0 bindings, 4.3.0 Tk + plugin contract, 4.3.1 GLUT |
| 4 | **5b — downloadable application**, alongside **2b — independent removals** | Tester launches without installing Python; each removal passes its own checks. | 5b done 2026-11-08 (4.4.0, Linux verified locally, CI artifacts for 3 OS); 2b done 2026-11-08 (4.4.1 to 4.4.8): chempy.fast, VTK-m, mpeg_encode, RPC/HTTP/web, OpenVR, COLLADA+libxml2, MMTF, dev utilities. Kept: VMD plugins/netcdf, xwin paste, POV-Ray renderer |
| 5 | **1b — measured optimization** and **3 — workspace increments** | Performance changes have before/after evidence; docks and table ship in reviewable test builds. | 3 slices 1-3 done 2026-11-08 (4.5.0: objects dock, layout persistence, Panels menu); 1b and slices 4-6 open |
| 6 | **4b — extension examples and proven helpers** | Command/menu integration documented; any dock helper has real consumers. | Small–medium |

Cleanup precedes speculative optimization. Recording a baseline is preparation.
UI work needs the binding decision and a usable test build, not completion of
every unrelated removal. GPU compute is not on the critical path. Plugin
compatibility work starts before Tk removal; final extension documentation can
follow the first UI increments. Artifact work can start after the binding change
while the remaining legacy cleanup proceeds.

### First concrete changes

1. Record an optimized baseline and current headless/GUI test status. Resolve
   existing failures before using those results to judge a change.
2. Promote the tested PySide6 binding to runtime metadata and establish a coherent
   Python minimum with dependencies and CI. Keep dev pins exact.
3. PR title `remove: alternative Qt bindings`; commit subject
   `remove(gui): standardize on PySide6`. Preserve used import aliases; keep Tk,
   GLUT, and renderer changes separate.
4. Deliver a first test artifact, then complete the Tk/plugin migration and other
   removals one at a time. Docks and the object table are the first visible UI work.

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
