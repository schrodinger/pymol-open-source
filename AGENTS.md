# AGENTS.md — rules for AI agents working on pymolish-studio

pymolish-studio is a fork of open-source PyMOL (C++ core + Python API + Qt GUI).
The goal of the fork: **smaller, faster, simpler, prettier**. Every change is
judged against those four words, in that order of ease: delete first, speed up
second, simplify third, beautify fourth.

## 1. Mindset: ponytail

Follow https://github.com/dietrichgebert/ponytail. Short form:

1. Does it need to exist? If speculative, skip it and say so in one line.
2. Stdlib / C++17 standard library does it? Use that.
3. Platform or already-linked library does it (Qt, numpy, OpenGL, GLM)? Use that.
4. Only then write the minimum code that works.

- Deletion beats addition. Boring beats clever.
- No abstraction with one user. No config for a value that never changes. No "for later" scaffolding.
- Shortest working diff wins. Fewest files touched.
- Mark deliberate shortcuts with a `ponytail:` comment naming the ceiling and the upgrade path.
- Non-trivial logic leaves one runnable check behind (a test in `testing/tests/` or a Catch2 case). Trivial one-liners need none.
- Output: code first, then at most three lines on what was skipped and when to add it.

## 2. Priorities (what to work on, in order)

### 2.1 Clean
- Remove dead code, unused features, legacy paths. Candidates already known:
  - `modules/pmg_tk/` (Tk GUI), `--glut` legacy GUI, `modules/pmg_qt/mimic_tk.py`, `mimic_pmg_tk.py`, Pmw mentions.
  - `--openvr` / `contrib/vr`, `--use-vtkm`, `modules/web/`, `pymolhttpd.py`, `rpc.py`, `xwin.py`, `mpeg_encode.py`, `povray.py`.
  - PyQt5 / PySide2 fallbacks in `modules/pymol/Qt/` once PySide6 is the only target.
  - Python < 3.11 compatibility shims (`requires-python` will move to `>=3.11`).
- Before deleting a feature: `grep -rn` for every reference (Python, C++, `.pml`, tests, `data/`), delete all of them in the same PR, run the test suite.
- One removal per PR. Title it `remove: <thing>`. The commit body lists what was deleted and why nothing depends on it.

### 2.2 Dependencies
- Runtime Python deps: `numpy` and `PySide6`. Nothing else without an explicit reason written in the PR.
- Native deps to keep: OpenGL, GLEW, GLM, libpng, freetype. Everything optional (libxml2, msgpack, netcdf, nlohmann-json, VTK-m, OpenVR, GLUT) is a removal candidate unless a user-facing feature we keep needs it.
- Never vendor a library when the system package or a header-only include does the job. Never add a dependency for something a few lines can do.
- Pin dev tools exactly in `pyproject.toml` `[project.optional-dependencies].dev`; keep runtime deps as floors (`>=`).
- Updating a library: bump, build, run tests, read the changelog for breaking changes, done. No wrapper layers "to isolate us from the library".

### 2.3 Speed
- Measure first. Use `testing/tests/performance/` or a `timeit`/`perf` run; put the before/after numbers in the PR.
- Order of attack: algorithmic fix > fewer allocations / copies > SIMD or OpenMP (`--use-openmp` already exists) > GPU.
- GPU means the shaders in `data/shaders/` and the CGO / Scene rendering path in `layer1/`. Prefer moving work into existing shaders over adding a compute framework. No CUDA / OpenCL / Vulkan without a measured 5x win on a real workload.
- Keep the GUI thread responsive: long work goes off-thread or into the existing C++ layer, never in a Python loop on the main thread.
- Don't micro-optimize code that isn't on a profile.

### 2.4 Simple and beautiful
- **Code:** small functions, clear names, no comment that repeats the code. C++ follows `.clang-format` (2-space, 80 col, Linux braces). Python follows `ruff format` defaults (add it; no black, no isort separately).
- **UI (Qt, `modules/pmg_qt/`):** one modern style, consistent spacing, native widgets over custom-painted ones, `.ui` forms over hand-built layouts when the form is non-trivial. Fewer menus, fewer dialogs, fewer toggles. Every new UI element must be reachable by keyboard.
- Prefer Qt's built-in features (QSS, QSettings, QAction shortcuts, QDockWidget) over new UI frameworks. No web views, no Electron, no QML unless the Qt widgets path is proven inadequate.

## 3. Repo map

| Path | What |
|---|---|
| `layer0/` … `layer5/` | C++ core, lowest to highest layer. `layer1/` = rendering (Scene, CGO, Ray, shaders, fonts). `layer3/` = Executive / high-level ops. `layer4/` = Python bindings, `layer5/` = main. |
| `layerGraphics/` | small graphics helpers |
| `ov/`, `contrib/` | vendored C libs (champ, mmtf-c, pocketfft, uiuc molfile plugins, vr) |
| `include/` | third-party headers |
| `data/shaders/` | GLSL, compiled into `ShaderText.cpp` by `create_shadertext.py` |
| `modules/pymol/` | Python API (`cmd`), one file per command family |
| `modules/pymol/Qt/` | PySide/PyQt compatibility shim |
| `modules/pmg_qt/` | Qt GUI (the only GUI we keep) |
| `modules/pmg_tk/` | legacy Tk GUI (removal candidate) |
| `modules/chempy/` | molecule data model |
| `testing/` | Python test framework; `testing/tests/api`, `/performance`, … |
| `layerCTest/` | Catch2 C++ tests (`--config-settings testing=True`) |
| `setup.py`, `_custom_build/`, `CMakeLists.txt` | build (setuptools + CMake). Build flags are in `setup.py` `parser.add_argument`. |

Note: `layer1/FontTTF.h` and `FontTTF2.h` are embedded font data, not code.

## 4. Build and test

```sh
# dev build, incremental, with C++ tests
pip install --verbose --no-build-isolation --config-settings testing=True '.[dev]'

# run the whole Python suite (headless)
pymol -ckqy testing/testing.py --run all

# run one test file
pymol -ckqy testing/testing.py testing/tests/api/<file>.py

# format
clang-format -i <changed .cpp/.h>
ruff format modules testing     # once ruff is in dev deps
```

CI (`.github/workflows/build.yml`) builds Linux / Windows / macOS on Python 3.11–3.14 and runs the suite. A PR that breaks CI is not done.

## 5. Change rules

- Never change behavior of a public `cmd.*` API without saying so in the PR and updating `ChangeLog`. Removal of a command is allowed (see 2.1) but must be deliberate, not a side effect.
- Keep file formats readable and writable: PDB, mmCIF, SDF/MOL, PSE, PNG export. Everything else may be cut.
- Don't reformat files you didn't otherwise change. Formatting-only PRs are separate and say so.
- Don't touch `LICENSE`, `AUTHORS`, copyright headers.
- Don't commit build output (`build/`, `generated/`, `*.pyc`).
- Commit messages: `<type>(<scope>): <summary>` with type in `remove | perf | simplify | ui | fix | deps | ci | docs`. Example: `remove(gui): drop Tk GUI and Pmw support`.
- Commit after each completed task. No AI attribution trailers in commit messages. Never force-push `master`.
- Version lives in `layer0/Version.h` (`_PyMOL_VERSION`). Scheme `4.X.Y`: finishing a roadmap step from `PLAN.md` bumps `X`; a follow-up improvement inside the current step bumps `Y`. Keep `_PyMOL_VERSION_int` in sync (`X.Y.Z` → `XYYYZZZ`); `get_version()` and PSE headers use it.

## 6. When unsure

- Pick the smaller change and say what you skipped.
- If a removal might break a feature a user cares about, do the removal in a branch and list what disappears in the PR description. Scaling scope down is the user's call.
- Don't ask permission for reversible work that follows from these rules. Do ask before deleting data files, changing file-format output, or adding a dependency.
