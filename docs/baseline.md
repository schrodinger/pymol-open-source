# Baseline: build loop and performance

Measured 2026-11-08 on commit a3105772 (code version 4.1.0). Rerun both
sections after any build or performance change and paste the new numbers into
the PR.

## Machine

| | |
| --- | --- |
| CPU | AMD Ryzen 9 8945HX, 16 cores / 32 threads |
| RAM | 17 GB |
| GPU | NVIDIA RTX 5070 Laptop, driver 621.01 (unused: all numbers below are headless) |
| OS | Ubuntu 24.04 under WSL2, Linux 6.6.87 |
| Toolchain | gcc 13.3, CMake 3.28, Python 3.12.4 in `.venv`, PySide6 6.12.0 |
| Build | release (`DEBUG` unset), `use-msgpackc=no`, VMD plugins on, OpenMP on, no ccache, no Ninja |

## Build loop (5a)

Command for every row except the test suite:

```sh
uv pip install --no-build-isolation --config-settings use-msgpackc=no -e '.[dev]'
```

| Step | Wall time |
| --- | --- |
| Cold build, `build/` removed | 54.4 s |
| Editable install from an existing `build/` (first editable switch) | 70.6 s |
| Editable reinstall, nothing changed | 30.3 s |
| Non-editable reinstall, nothing changed | 27.2 s |
| One C++ file touched (`layer2/RepSphere.cpp`), non-editable | 7.9 s |
| One shader touched (`data/shaders/sphere.fs`), non-editable | 10.4 s |
| One Python file touched, non-editable reinstall | 9.6 s |
| One Python file touched, editable | 0 s (no command needed, `pymol.__file__` points into `modules/`) |
| Headless suite `pymol -ckqy testing/testing.py --run all` | 14.2 s wall, 997 tests |
| `build/` size | 36 MB |
| Wheel `pymol-*-cp312-linux_x86_64.whl` | 8.3 MB |
| PyInstaller one-folder bundle, Linux (PySide6 pruned of Quick/Qml/Pdf/Network) | 292 MB unpacked, 114 MB tar.gz |

Notes:
- The editable route is validated: the suite passes with the same 4 MMTF-only
  errors as the wheel install (msgpack headers not installed on this machine).
- A no-op reinstall costs 27 to 30 s in both modes, while a reinstall after
  touching one file, measured back-to-back, costs 8 to 10 s. The difference is
  outside CMake (uv resolution and setuptools overhead) and was not
  investigated. For Python edits use the editable install and skip the
  reinstall entirely.
- `--config-settings jobs=N` and `JOBS=N` now limit CMake parallelism on
  Linux and macOS; before this change the option was parsed and ignored.
- ccache: not installed here, so not measured. CMake honors
  `CMAKE_CXX_COMPILER_LAUNCHER=ccache` from the environment without code changes.

## Size after the cleanup steps (4.1.0 to 4.4.8)

Lines of C++ and Python per directory, fork point `3250e7a3` versus `4.4.8`.
Overall diff since the fork point: 255 files, +1584 / -46358 lines.

| Directory | Before | After |
| --- | --- | --- |
| layer0 | 38039 | 37616 |
| layer1 | 343734 | 341395 |
| layer2 | 73405 | 73063 |
| layer3 | 42170 | 41946 |
| layer4 | 7975 | 7866 |
| layer5 | 5926 | 4469 |
| modules/pymol | 52641 | 50020 |
| modules/pmg_qt | 7247 | 7104 |
| modules/pmg_tk | 5479 | removed |
| modules/chempy | 21797 | 21219 |
| contrib | 45668 | 39295 |
| Wheel (Linux, cp312) | 8.3 MB | 8.2 MB |

Native dependencies dropped: libxml2, mmtf-cpp, OpenVR, VTK-m, GLUT, Tk/Pmw.

## Performance (1a)

Command: `pymol -ckqy testing/bench.py` (5 repeats, median, `max_threads` =
CPU count, headless, no GL). The map is in sigma units
(`normalize_ccp4_maps` on).


Commit a3105772, PyMOL 4.1.0, Python 3.12.4, Linux x86_64, 32 CPUs, PSE 16.8 MB

| Workload | median | min..max | peak RSS |
| --- | --- | --- | --- |
| load 1aon.pdb.gz (58870 atoms) | 0.167 s | 0.158..0.233 s | 104 MB |
| ray 800x600, cartoon (1aon) | 4.111 s | 4.017..4.487 s | 2408 MB |
| ray 800x600, surface (1aon) | 15.091 s | 14.822..15.950 s | 2745 MB |
| isosurface emd_1155.ccp4 level 2.0 (+ray 64x64) | 1.053 s | 1.011..1.108 s | 3070 MB |
| isosurface level change 2.0 -> 3.0 (+ray 64x64) | 0.091 s | 0.085..0.093 s | 3334 MB |
| save PSE (1aon, cartoon) | 0.187 s | 0.181..0.191 s | 3334 MB |
| load PSE (1aon, cartoon) | 1.020 s | 1.006..1.063 s | 3334 MB |
| startup `pymol -ckq -d quit` | 0.082 s | 0.079..0.085 s | nan MB |

Observations worth a profile later, in order of size:
- Isosurface at 1σ on this map takes 10.7 s and 9 GB peak RSS (one-off probe,
  not in the table). 2σ is 1.4 s. Mesh memory, not marching cubes, is the issue.
- CPU ray of a 58 k-atom surface: 15 s at 800x600. Cartoon: 4.3 s. Both at
  2.4 to 2.7 GB peak RSS.
- Loading the PSE (17 MB) takes 6x longer than parsing the gzipped PDB it was
  made from.
- Startup is 80 ms headless, not a target.

Not measured: GUI frame latency and GPU time (needs a real GL context, see
`testing/tests/performance/shaders.py`), Windows, macOS.
