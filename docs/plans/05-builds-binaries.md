# 5 — Fast builds and downloadable test binaries

[Roadmap](../../PLAN.md) · Draft; no build timings or bundle results yet

**Shorten the edit/build/test loop, then provide applications testers can launch
without Python.** An installable wheel and a standalone application are separate
deliverables. Build speed, launch speed and download size also need separate data.

## Status 2026-11-08

5a done: build-loop timings, editable install validated, `--jobs` fix. See [docs/baseline.md](../baseline.md).

5b done in 4.4.0: `packaging/pyinstaller.spec` (PyInstaller one-folder, data under `pymol/pymol_path`), wheel + bundle uploaded as CI artifacts from the Python 3.13 job on Linux, Windows and macOS (arm64 Miniforge now). Verified locally on Linux: bundle runs headless and with the GUI from another directory; see baseline for sizes. Windows and macOS bundles are unverified until CI runs; signing and installers not started.

4.5.1: bundles are archived (`.tar.gz` Linux/macOS, `.zip` Windows) and a push to `master` publishes them with the wheels on a GitHub pre-release tagged `v<version>`; the workflow also runs on `workflow_dispatch`. Downloads: the repository's Releases page.

4.5.3: CI runs only for pushes/PRs to `master`, on Python 3.13 and 3.14. Windows and macOS get their native libraries through pixi (`pixi.toml`, conda-forge only; the upstream `schrodinger` channel and the msgpack-c git clone are gone; catch2 is `2.13.*` because 2.13.3 has no osx-arm64 build). Linux keeps apt. Releases carry a `SHA256SUMS` file. Note: this fork has recorded zero Actions runs so far; Actions must be enabled once in the repository's Actions tab before any of this executes.

2026-11-16: 5c proposed below (builds without GitHub Actions). CI is red on one test, not on size; nothing in the workflow was changed.

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
   pymolish -ckqy testing/testing.py --run all
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

- Preserve Linux/Windows/macOS and Python 3.13–3.14 tests on source PRs. Add cache
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

## 5c — Release builds without GitHub Actions (proposal, 2026-11-16)

Question: can the Windows and macOS bundles be built on the development
machine (Windows 11 x64 host, Ubuntu 24.04 in WSL2, 32 threads, 332 GB free on
`C:`) instead of on GitHub Actions? Short answer: Linux yes, Windows yes after
a one-time compiler install, macOS no.

### What was measured first

- Bundle size is not what stops Actions. The `workflow_dispatch` run of
  2026-11-13 built all three bundles there: Linux 111 MB, macOS 79 MB, Windows
  75 MB. With the `chem` and `fix` extras the Linux bundle built locally today
  is 419 MB unpacked and 142 MB as `.tar.gz`. A release file may be up to
  2 GiB, and a public repository pays nothing for minutes or artifact storage.
- Every run since 2026-11-14 fails in the **Test** step on Linux and Windows,
  on one test: `testLoadEmbedded` (`Error: embedded data 'wats' not found`).
  The bundle steps come after the tests, so no bundle is made. The macOS job
  passes and uploads its bundle on every run. The same test passes locally
  (`--run all` on `feature6`: 1060 tests, OK). The cause is not found yet.
- So the cheapest way to get all three bundles is to fix that one test. The
  rest of this section is the route that does not depend on Actions.

### Where each target can be built

| Target | Where | Verdict | Why |
|---|---|---|---|
| Linux x86_64 | `podman` on this machine, image `ubuntu:22.04` | Works now | podman 4.9.3 is installed. The container keeps the glibc 2.35 baseline of the CI runner; a bundle built directly in WSL (glibc 2.39) does not start on older distributions. |
| Windows x86_64 | The Windows host of this machine, started from WSL | Works after a one-time install | WSL can run `pwsh.exe`. The host has git and winget but no MSVC, no pixi and no CMake. |
| macOS arm64 | Not on this machine | Not possible | See below. |

### Windows: build on the host, not in a container

One-time setup on the Windows side:

```powershell
# needs an administrator once; the account in use is a standard user
winget install Microsoft.VisualStudio.2022.BuildTools --override "--passive --wait --add Microsoft.VisualStudio.Workload.VCTools --includeRecommended"
# per user, no administrator
winget install prefix-dev.pixi
```

Each build, from WSL. The commands inside are the ones of the `build-Windows`
job in [build.yml](../../.github/workflows/build.yml), unchanged:

```sh
git clone --depth 1 "file://$PWD" /mnt/c/b/pymolish   # NTFS path: MSBuild is unreliable on \\wsl.localhost
cd /mnt/c/b/pymolish && pwsh.exe -NoProfile -Command '
  pixi run -e py313 pip wheel -v --no-deps -w dist --config-settings testing=True --config-settings json=True .
  $whl = (Get-ChildItem dist\*.whl | Select-Object -First 1).FullName
  pixi run -e py313 pip install "$whl[dev,chem,fix]"
  pixi run -e py313 pymolish -ckqy testing\testing.py --run all
  pixi run -e py313 pyinstaller --noconfirm packaging/pyinstaller.spec
  Compress-Archive -Path dist\pymolish -DestinationPath dist\pymolish-Windows-X64.zip'
```

Then start `dist\pymolish\pymolish.exe` from another directory: this is a real
desktop with a real GPU, which the CI runner is not. Publishing stays the
`gh release upload` line of the `release` job, run by hand from WSL.

Not verified: none of the above ran, because MSVC is not installed. Open
points for the first run: whether `create_buildinfo.py` is content with a
shallow clone, and long paths (hence the short `C:\b`).

Rejected for Windows:

- **Windows containers.** They need the Windows Containers feature and an
  administrator, a base image of several GB, and MSVC inside the image anyway.
  The host already is Windows; a container adds work and removes nothing.
  podman in WSL runs Linux containers only: a container shares the kernel of
  its host.
- **Cross-compiling in a Linux container** (MinGW or clang-cl, plus Wine).
  PyInstaller does not cross-freeze: it must run on the system it targets.
  Running it under Wine with PySide6 is unsupported and would ship a binary
  that no test ever ran natively.
- **MinGW through pixi, without an administrator.** `setup.py` has a MinGW
  branch, so this may work, but it yields a different binary from the MSVC one
  CI tests. Try it only if nobody can install the Build Tools.

### macOS: why not on this machine

1. **No macOS containers exist.** A container uses the kernel of its host;
   there is no macOS kernel on a Windows or Linux machine. The "macOS in
   Docker" projects are virtual machines in disguise.
2. **A macOS virtual machine is not an option either.** Apple's license allows
   macOS only on Apple hardware, and the target is arm64: Apple-silicon macOS
   does not boot in a VM on an x86 processor.
3. **Cross-compiling is half a solution.** The C++ extension could in theory
   be built with a cross toolchain and Apple's SDK (same license limit), but
   PyInstaller still has to run on macOS to make the bundle, and nothing built
   here could be started here. An untested bundle fails the acceptance rule at
   the end of this plan.

Ways to get the macOS bundle, cheapest first:

1. Keep only the `build-MacOS` job on GitHub Actions. It is green today, takes
   2 to 6 minutes and produces an 81 MB bundle.
2. A real Mac (owned or borrowed): install pixi, run the three commands of the
   `build-MacOS` job. Signing and notarization also need a Mac.
3. A rented Mac by the hour (for example AWS EC2 Mac, MacStadium, Scaleway).
   Costs money; only worth it if 1 and 2 are both ruled out.

### Proposed order

1. Decide: fix `testLoadEmbedded` on CI (all three bundles come back, no new
   setup), or go local for Linux and Windows.
2. If local: get the Build Tools installed, run the Windows block once by
   hand, fix what breaks, then put both blocks into one script
   `packaging/build_local.sh` (Linux in podman, Windows on the host). No script
   before the commands have run once.
3. macOS stays on Actions until a Mac is available.

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
