# 2 — Code, feature, and dependency cleanup

[Roadmap](../../PLAN.md) · Draft; candidates are not deletion approvals

**Remove one unsupported feature completely, then measure the result.** Identify
every retained consumer before removing a dependency. Preserve PDB, mmCIF,
SDF/MOL, PSE and PNG workflows throughout.

## 2a — GUI foundation

1. Follow [step 0](00-ui-engine.md) for explicit PySide6 runtime metadata and a
   separate alternative-binding removal. Preserve aliases with live users.
2. Trace plugin startup/menu registration through
   [legacysupport.py](../../modules/pymol/plugins/legacysupport.py),
   [plugins/__init__.py](../../modules/pymol/plugins/__init__.py),
   [Qt initialization](../../modules/pmg_qt/pymol_qt_gui.py), and
   [installation.py](../../modules/pymol/plugins/installation.py).
   `legacysupport.py` imports `pmg_tk.startup`; this is a live dependency.
3. Move needed discovery/menu/dialog behavior into the existing Qt path,
   preserving the contracts documented in [step 4](04-extensibility.md).
   Explicitly list Tk-only plugins that will cease to work.
4. Remove Tk and its now-unneeded adapters as one feature removal after the
   preparation passes tests. Migrate retained startup content first. Remove GLUT
   in a separate PR with its own audit and GUI/headless checks.
5. Raise Python's floor to 3.11 coherently with metadata, dependencies, installers
   and CI. Delete obsolete shims after that floor takes effect; old-looking syntax
   alone is not evidence that a branch is unused.

## 2b — Candidate removal queue

This is a starting inventory, not a completed reference audit. Reorder by proven
lack of use and maintenance cost; optional does not mean dead.

| Candidate | Feature that disappears | Check before removal |
| --- | --- | --- |
| OpenVR / `contrib/vr` | VR interaction | Invocation/build flags, input/scene hooks, shaders and resources |
| VTK-m | Alternate contouring path | Native contouring still meets retained workloads and produces correct results |
| `modules/web`, `pymolhttpd.py`, `rpc.py` | Browser/server/remote control | Separate independently removable features; inspect API tests and examples |
| `xwin.py` | External-window integration | Platform callers and launch options |
| `mpeg_encode.py` | That movie encoder | Keep PNG sequences and explicit remaining encoding routes; inspect fallbacks in `movie.py` |
| `povray.py` | External POV-Ray integration | Trace ray modes and `.pov` output separately; preserve native CPU ray tracing |
| COLLADA / libxml2 | COLLADA export | Inspect any retained glTF fallback through COLLADA |
| MMTF / msgpack / mmtf-c | MMTF files and related fetch routes | Audit local-file/API consumers, not only server format preferences |
| VMD molfile plugins / netcdf | Potential trajectory/volume readers | Inventory formats individually; retain the subset actual workflows require |
| Unused development utilities | Their maintenance burden | Search CI, docs and script references first |

Build switches live in [setup.py](../../setup.py). Native glTF/GLB needs
`json=True` and nlohmann-json, which [CI](../../.github/workflows/build.yml)
enables. [exporting.py](../../modules/pymol/exporting.py) documents an external
fallback too. Keep this capability/dependency unless a separate feature decision
changes it. glTF does not replace encoded movie output.

## Dependency policy and updates

| Class | Plan |
| --- | --- |
| Runtime Python | numpy and PySide6 with tested floors; no unconditional startup imports of optional analysis libraries |
| Dev/test Python | Keep packages with actual consumers; pin exactly, including the currently unpinned Biopython entry |
| Core native | Keep OpenGL, GLEW, GLM, libpng and freetype |
| Feature-dependent native | Retain libxml2, msgpack, netcdf, nlohmann-json and others until retained feature requirements are resolved |
| New tooling | Ruff is a proposed dev addition under the dependency rule; prefer it over separate formatter/import-sorting packages |

`pytest` and Pillow serve the test framework; requests serves API tests;
Biopython serves [sequence alignment](../../modules/pymol/seqalign.py). Absence
from runtime requirements does not make a test dependency unnecessary. Retained
optional features need an installation route and a useful missing-package error.

For each update, read publisher release notes and Python/platform support, choose
an available version, build and run appropriate checks. Do not invent future pins
or upgrade the native stack wholesale. Runtime floors and exact dev pins differ.

## Removal procedure

1. Search code, headers, `.pml`, examples, tests, data references, build scripts
   and docs for imports, symbols, macros and flags. Starting example:
   `rg -n --hidden -g '!.git/**' 'pmg_tk|mimic_tk|mimic_pmg_tk|Pmw' .`.
   Inspect dynamic imports and registration tables too.
2. List disappearing behavior and remaining callers in the PR. Preserve retained
   behavior first; do not delete its regression tests to hide breakage.
3. Remove implementation, flags, imports, menus, dedicated tests and live docs
   references together. Keep historical `ChangeLog` records and protected notices.
   Data-file deletions follow the repository's permission rule.
4. Build; run `pymol -ckqy testing/testing.py --run all`, relevant Catch2 checks,
   affected GUI smoke checks and protected file round trips.
5. Record removed files/dependencies and useful build/artifact size changes.
   PR title: `remove: <thing>`; commit: `remove(<scope>): <summary>`. Public
   `cmd.*` behavior changes still require `ChangeLog` updates.

## Additional ideas and boundaries

- **Feature/dependency table:** document the final formats and optional
  integrations so later removals have a clear contract.
- **Reproducible native inputs:** CI clones moving MMTF/msgpack branches. If
  retained, use immutable revisions or suitable system packages; do not vendor
  another copy.
- **Import/package audit:** compare imported modules with bundled contents;
  remove accidental inclusion only after import and resource checks.
- **Formatting:** add a verified exact Ruff pin, then use a separate
  formatting-only change instead of mixing churn with removals.

Skip GLEW replacement, `layerGraphics/` relocation, arbitrary setup.py line limits,
and font/icon/data deletion unless a concrete maintenance or dependency problem
disappears with the change.
