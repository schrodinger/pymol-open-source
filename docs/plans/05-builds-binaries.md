# 5 — Fast builds and downloadable test binaries

[Roadmap](../../PLAN.md) · Draft; no build timings or bundle results yet

**Shorten the edit/build/test loop, then provide applications testers can launch
without Python.** An installable wheel and a standalone application are separate
deliverables. Build speed, launch speed and download size also need separate data.

## Status 2026-11-08

5a done: build-loop timings, editable install validated, `--jobs` fix. See [docs/baseline.md](../baseline.md).

5b done in 4.4.0: `packaging/pymolish.spec` (PyInstaller one-folder, data under `pymol/pymol_path`), wheel + bundle uploaded as CI artifacts from the Python 3.12 job on Linux, Windows and macOS (arm64 Miniforge now). Verified locally on Linux: bundle runs headless and with the GUI from another directory; see baseline for sizes. Windows and macOS bundles are unverified until CI runs; signing and installers not started.

4.5.1: bundles are archived (`.tar.gz` Linux/macOS, `.zip` Windows) and a push to `master` publishes them with the wheels on a GitHub pre-release tagged `v<version>`; the workflow also runs on `workflow_dispatch`. Downloads: the repository's Releases page.

## Current state and risks to verify

- [Build backend](../../_custom_build/backend.py) supports wheels and editable
  builds; [setup.py](../../setup.py) invokes CMake per extension and already uses
  parallel compilation. Do not add a second build backend to obtain these basics.
- `--jobs`/`JOBS` are declared, but `build_cmake` constructs its non-Windows build
  parallelism from `os.cpu_count()`. Trace and make the limit effective before
  assuming the documented setting limits memory use.
- Generated shader content avoids writes when unchanged. Build-info generation
  embeds the current time unless `SOURCE_DATE_EPOCH` is set; measure whether this
  invalidates incremental work before changing its behavior.
- [CI](../../.github/workflows/build.yml) defines 12 OS/Python combinations;
  Linux sets `DEBUG=1`. These tests are not release-performance measurements.
- The macOS job uses `macos-latest` but downloads an x86_64 Miniforge installer.
  Make runner, interpreter, native libraries and output architecture explicit
  before claiming an arm64 artifact; do not infer actual architecture from a label.
- Installation copies resources outside ordinary Python modules. Check
  `PYMOL_PATH`, `.ui` files, data, fonts, icons and shared libraries in the actual
  wheel/bundle, not just from a source checkout.

## 5a — Baseline and local iteration first

1. Record cold release build, unchanged rebuild, one C++ source rebuild, shader
   edit rebuild, Python-only edit and test times on a named machine. Include peak
   memory and output sizes. The draft's sub-minute rebuild is a target to evaluate,
   not an established result.
2. Document the existing build/test commands with prerequisites. For an already
   provisioned development environment, validate the editable route:

   ```sh
   python -m pip install --verbose --no-build-isolation --config-settings testing=True --editable '.[dev]'
   pymol -ckqy testing/testing.py --run all
   ```

   Confirm the custom backend's editable behavior before adopting this as the
   recommended loop. Build prerequisites must already be installed for
   no-build-isolation; unset `DEBUG` for release timing (`DEBUG=0` is still truthy
   in the current build script).

3. Reuse build directories for compatible local builds. Verify Python edits do
   not rebuild C++; make the existing job limit effective and test a bounded
   parallel build on a memory-constrained machine.
4. Evaluate ccache on supported Make/Ninja toolchains. Prefer the standard CMake
   compiler-launcher setting, with a tiny explicit pass-through only if required.
   The launcher environment variable requires CMake 3.17, while current metadata
   permits 3.13.3; do not silently assume it works at the floor. Windows Visual
   Studio needs separate generator/cache validation. [CMake launcher docs](https://cmake.org/cmake/help/latest/envvar/CMAKE_LANG_COMPILER_LAUNCHER.html)
5. Measure cache hit rate and rebuild time. If build timestamps cause avoidable
   work, use an explicit source epoch for reproducible artifacts while preserving
   accurate commit/version metadata. Do not introduce stale generated output.

No mandatory `uv`, Make wrapper, Ninja, or new compiler is required initially.
Add tooling only if measurement shows it reduces work enough to justify upkeep.

## 5b — Wheels and the first application bundle

### Wheels first

- Build a wheel with the existing backend and install it into a clean environment
  outside the checkout. Verify imports, launch, protected formats and data paths.
- Initially use one Python version, provisionally 3.12, for preview artifacts;
  keep the full supported Python test matrix. Python-extension wheels remain
  interpreter/ABI-specific; one wheel is not universal.
- Use existing native OS jobs first. Consider `cibuildwheel` only when publishing
  portable wheels across targets warrants its extra tooling; it builds/tests
  platform wheels but does not freeze a standalone GUI app.
  [cibuildwheel](https://cibuildwheel.pypa.io/en/stable/)
- For broadly distributed Linux wheels, define the compatibility baseline and
  validate/repair shared libraries. For all platforms, inspect dependent native
  libraries rather than relying on the builder's environment.

### Standalone application

1. Evaluate **PyInstaller one-folder** as the first packaging experiment, built
   from the tested wheel in a lean runtime environment. It bundles Python and
   dependencies, and output is platform-specific. One-folder is easier to inspect
   than a self-extracting one-file app. [PyInstaller operating model](https://pyinstaller.org/en/stable/operating-mode.html)
2. Compare Qt's existing `pyside6-deploy` if PyInstaller cannot package this app
   cleanly. It uses Nuitka and is not a dependency-free shortcut; compare build
   time, startup, size and plugin behavior before choosing.
   [Qt deployment tool](https://doc.qt.io/qtforpython-6/deployment/deployment-pyside6-deploy.html)
3. Prove one target end-to-end, then expand. Initial targets are Linux x86_64,
   Windows x86_64 and macOS arm64 after architecture alignment. The first local
   proof can be Linux; Windows `.exe` delivery is an explicit early milestone.
4. Ship Windows as a ZIP containing the executable and required files, Linux as
   an archive with launcher, and macOS as a `.app` archive. A Windows `.exe` in
   one-folder mode needs its neighboring files. Installers and DMGs can follow.
5. Include Qt platform/image plugins, PyMOL extension libraries, numpy, forms and
   required resources/notices. Test dynamic imports and supported user plugin
   discovery explicitly. Exclude dev tools only when import/resource checks pass.
6. Keep a wheel-based environment available for plugins requiring extra native
   packages. Define the first frozen bundle's supported external-plugin scope;
   do not claim arbitrary pip-installed extensions will work in it.

PyInstaller or other new tooling is a proposed dev dependency, not installed by
this plan. Choose an available exact pin and follow the repository dependency
rule when implementing. Keep packaging requirements out of core runtime deps.

## CI and artifact delivery

- Preserve Linux/Windows/macOS and Python 3.11–3.14 tests on source PRs. Add cache
  keys scoped to OS, architecture, compiler, Python ABI and dependency/build
  options; never reuse incompatible objects just to improve a hit rate.
- Separate fast preview artifacts from full release packaging. Build previews
  for relevant successful changes or manual requests; avoid regenerating native
  artifacts for documentation-only changes. Inspect actual cost before changing
  required test coverage.
- Name artifacts with commit, OS, architecture and build mode. Include concise
  launch instructions, checksums, supported optional features and known limits.
  Start with CI artifacts; release publishing is a separate release action.
- Pin retained native source fetches to immutable revisions. Keep test and
  performance configurations explicit; an accidental debug binary is not a
  representative performance preview.
- Before public releases, establish platform signing/notarization and installer
  handling as needed. For previews, identify unsigned status honestly; do not
  instruct testers to disable system protection globally.

## Additional ideas and acceptance

- **Version visible in the app:** reuse the existing About/version path to show
  the artifact's commit so feedback is attributable to the correct build.
- **Small smoke fixture:** bundle or otherwise supply a permitted local example
  for launch, rotation, selection, PNG export and PSE reopen without a download.
- **Artifact comparison:** report compressed size and startup time alongside
  build time before adding one-file packaging or more compiler optimization.
- **Reproducible feedback:** include the tested OS/GPU/Qt details and a short
  workflow checklist with each preview, without adding telemetry.

A milestone passes only after the downloaded artifact runs from a different
directory on a clean target machine without Python or the source checkout. Test
spaces/non-ASCII paths, file dialogs, a real GL context, supported plugins and
protected formats. Headless CI smoke checks supplement real desktop checks.
Building a wheel alone does not complete the executable-delivery request.
