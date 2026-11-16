# pymolish-studio — reviewed roadmap

First reviewed 2026-11-08 against revision `3250e7a3`. Checked again 2026-11-16
at 4.12.6 (`f30d3525`): the open steps were compared with the code, the detailed
plans and the git log. Nothing was built or measured for that check.

**Keep Qt Widgets with PySide6, simplify the supported product, measure performance,
modernize the existing window, and deliver downloadable builds early.** Changing
bindings does not move molecular computation to the GPU. Follow
[AGENTS.md](AGENTS.md): smaller, faster, simpler, prettier.

Finished steps are in [PLAN-DONE.md](PLAN-DONE.md), and so are the shipped
milestones of the steps that are still open. A step moves there when its last
milestone ships; its detailed plan stays in `docs/plans/`.

## Open steps

| Step | Detailed plan | What is left |
| --- | --- | --- |
| 1 | [Performance and GPU](docs/plans/01-performance-gpu.md) | Renderer default on WSL (waits for the user), draw cost of sticks and surfaces, the two triangulation passes, the ray tracer's voxel table. One measured bottleneck per change. |
| 3 | [Modern workspace and metadata](docs/plans/03-modern-ui.md) | Scenes as a dock, rows moved by dragging, a picked atom selects its row, sequence viewer and scene buttons still drawn by the core, spacing of the remaining dialogs, removal of the old C++ panel (waits for the user). Slice 6, menus, is step 14. |
| 8 | [Modern C++ in the core](docs/plans/08-cpp-modernization.md) | No rewrite. One idiom per PR behind sanitizers and tests; vendored C stays; files go in profile order so cleanup and step 1 meet. Next: 8h (run Catch2 in CI), 8i (dead `#ifdef` branches), then 8c to 8g. |
| 14 | [Action and right-click menus](docs/plans/14-menus.md) | 14c: apply the table of 14b. Rows 1 to 4 are dead rows and one bug fix; rows 5 to 8 move rows people know. Both wait for the user. |

## Validation of the original draft

Rows that still steer an open step. The rows about finished steps are in
[PLAN-DONE.md](PLAN-DONE.md).

| Original idea or claim | Verdict and correction |
| --- | --- |
| Require a Core profile during binding cleanup | **Reject as an early change.** `Scene.cpp` and `Ortho.cpp` still use fixed-function calls. Audit and migrate required rendering paths separately. |
| Enable OpenMP by default | **Already largely true.** `setup.py` defaults it on outside macOS. Marching cubes uses it; `Ray.cpp` has its own workers. Check races and oversubscription before adding parallelism. |
| Move PDB/mmCIF parsing from Python to C++ | **Correct the premise.** Native readers already exist in `ObjectMolecule.cpp`, `ObjectMolecule2.cpp`, and `CifMoleculeReader.cpp`. Profile the actual load path. |
| GPU ray tracing and marching cubes are next | **Defer.** OpenGL 4.3 compute cannot be the common requirement for native macOS OpenGL. Step 1 gives the platform constraint and experiment gate. |
| Replace the whole window beside the old one | **Prefer incremental changes.** Two windows duplicate command, dialog, plugin, and lifecycle integration. Extract only the panel being implemented. |
| Object and sequence panels can simply become docks | **Underestimated.** Existing internal controls involve C++ drawing; a Qt replacement needs separate behavior and performance checks. The Objects panel was rebuilt in Qt in 4.5.11; the sequence viewer is still drawn by the core. |
| Test PRs only on Linux/Python 3.12 | **Reject for this migration.** Preserve three-OS and Python 3.11–3.14 coverage; optimize caching and artifact jobs first. **Not what CI does:** it runs 3.13 and 3.14 on three systems while `pyproject.toml` allows 3.11. See "Waiting for the user". |
| Replace GLEW, collapse directories, cap setup.py at 300 lines | **Skip without evidence.** These can add coupling and churn without improving a measured workflow. |

Evidence and implementation boundaries are linked in each step. “Keep” means a
recommended direction, not a claim that an implementation has passed tests.

## Open work, in execution order

The topic numbers remain the user's original numbers. The letters identify
milestones in the detailed plans, not new topics.

| Order | Work | Exit condition | State |
| --- | --- | --- | --- |
| 5 | **1b — measured optimization** and **3 — workspace increments** | Performance changes have before/after evidence; docks and table ship in reviewable test builds. | Open in 1b: renderer default on WSL, sticks and surface draw cost. Open in 3: the list in the table above. Shipped parts are in [PLAN-DONE.md](PLAN-DONE.md). |
| 14 | **14c — menu consolidation** | Approved keep, merge, drop list applied; fewer rows and fewer lines in `menu.py`. | Medium. Waits for the user's answer on the table in `docs/plans/14-menus.md`. `menu.py` is at 1,809 lines after 4.12.5. |
| none | **8h, 8i, then 8c to 8g — modern C++** | Each slice passes the suite, the sanitizer job, `bench.py` and a PSE round trip. | Step 8 never got an order number. 8h is one CI line. 8i needs `unifdef` on the development machine. Many PRs, no date. |

Order 14 does not block order 5. One overlap: slice 6 of step 3 (menus) is 14c
and is tracked only there. GPU compute is not on the critical path.

### Working rules

- One slice per session and per PR. Each plan names the files and lines to
  open, so a session starts there instead of exploring again.
- [ponytail](https://github.com/dietrichgebert/ponytail) is already the rule
  (`AGENTS.md` §1): each plan has a "Skipped" list and, where code goes away,
  a "Deleted" list. Nothing on a "Skipped" list is built without a new request.
- [caveman](https://github.com/JuliusBrussee/caveman) shortens replies only;
  optional. The larger saving is the narrow slice and a fresh session per slice.

## Waiting for the user

| # | Decision | Where it is described | What it holds up |
| --- | --- | --- | --- |
| 1 | Menu table: approve rows 1 to 4 (dead builders, a row that always raises, a copy-paste bug, a row that labels nothing); say which of rows 5 to 8 (two merges, the regrouped A menu, the flat right-click rows) to apply. | [14-menus.md](docs/plans/14-menus.md), "14b result" | 14c, all of it |
| 2 | On WSL, should the app select the `d3d12` driver itself? Frames go from 50 to 120 ms down to 0.7 to 11 ms; lines are then 1 pixel wide, so `trilines` goes on with it and mesh lines stay thin. | [01-performance-gpu.md](docs/plans/01-performance-gpu.md), Status 2026-11-14 | The first open item of 1b |
| 3 | Merge to `master`. `master` is at 4.8.0; 18 commits (steps 10 to 13, 14a, 14b, 4.12.4 to 4.12.6) are on `feature5` and `feature6` with no pull request, and the last two are not pushed. | `git log origin/master..feature6` | Release bundles, and every check marked "first CI run": the icon of the Windows program, the release notes, the `chem` and `fix` extras on three systems |
| 4 | Python floor: `requires-python = ">=3.11"`, but CI tests 3.13 and 3.14 only. Raise the floor to 3.13, or add 3.11 and 3.12 jobs. | `pyproject.toml`, `.github/workflows/build.yml`, the row above | Nothing; 3.11 and 3.12 are promised and never tested |
| 5 | Are merges from upstream PyMOL still wanted? No upstream remote is configured. If yes, step 8 is limited to files this fork already changed heavily. If no, the reason for "no whole-tree `ruff format`" in `AGENTS.md` §2.4 is gone. | [08-cpp-modernization.md](docs/plans/08-cpp-modernization.md), Risk | The scope of 8i and 8c to 8g |
| 6 | Delete the old C++ object panel (Display > Internal GUI)? The Qt panel replaced it in 4.5.11. `pymolish -x` still uses it. | [03-modern-ui.md](docs/plans/03-modern-ui.md), Objects panel | One removal PR |
| 7 | Rights line for the icon and the logo: who made them, which license. | [09-icons-logo.md](docs/plans/09-icons-logo.md), `docs/licenses.md` | A complete license list; no code |
| 8 | The product name contains the mark "PyMOL": keep it and accept the risk, ask Schrodinger in writing, or rename. Kept so far. | [10-readme-license.md](docs/plans/10-readme-license.md), Wording | Nothing |
| 9 | Dropped plugin commands to bring back: the group variants and `fetch_uniprot_ligands`. | [12-native-commands.md](docs/plans/12-native-commands.md), the keep and drop table | Nothing |

Not decisions, but the user's to do: remove the pymolish plugin from
`~/.pymol/startup` (it still wins over the native commands), and look at the
Windows artifact once item 3 has run.

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
