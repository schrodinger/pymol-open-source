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
  not in the table). 2σ is 1.4 s. Corrected in 1b below: the cost is the ray
  tracer's setup for 3.2 M triangles, not the mesh.
- CPU ray of a 58 k-atom surface: 15 s at 800x600. Cartoon: 4.3 s. Both at
  2.4 to 2.7 GB peak RSS.
- Loading the PSE (17 MB) takes 6x longer than parsing the gzipped PDB it was
  made from. Fixed in 4.5.10, see below.
- Startup is 80 ms headless, not a target.

Not measured here: GUI frame time (see "GUI frame time" below, 4.7.9),
Windows, macOS.

## 1b results

### Session load (4.5.10)

Profile (`cProfile` around `cmd.load('bench.pse')`): all of the time was in
Python's pure-Python unpickler (`pickle._loads`, 30.9 M Python calls), used
since the Python 3 port so that Python 2 strings decode as text or bytes.
`chempy.pkl` now tries the C unpickler with the `undefined` codec, which
raises on the first Python 2 string, and falls back to the old path. Sessions
written by Python 3 hold no such strings and get the same objects back.
Check: `testing/tests/api/pkl.py`.

| Workload (`pymolish -ckqy testing/bench.py`) | 4.5.9 | 4.5.10 |
| --- | --- | --- |
| load PSE (1aon, cartoon, 16.8 MB) | 0.92 s (two runs: 0.88, 0.95) | 0.139 s (0.138..0.177) |
| load 1aon.pdb.gz, for scale | 0.161 s | 0.161 s |

Python 3.13.12, same machine as above. The 4.1.0 table above has 1.020 s.
Not faster: sessions saved by Python 2 builds or with `pse_export_version`
below 1.9, and the bundled `.pkl` data files (Python 2 pickles, 0.1 s for the
largest). The remaining 0.14 s is the C++ `set_session`.

### Surface build (4.5.12)

Where "ray 800x600, surface (1aon)" goes, from the ray tracer's own feedback
(`feedback enable, ray, blather`) and temporary timers in `SurfaceJobRun`:

| Part | Time |
| --- | --- |
| Building the surface (58870 atoms, 452340 points, 905126 triangles) | 9.5 s |
| ... solvent dots (`SolventDotNew`) | 1.5 s |
| ... dots to surface points | 2.7 s |
| ... refine and eliminate | 1.9 s |
| ... triangulation (`TrianglePointsToSurface`) | 3.5 s |
| Ray tracing the finished surface at 800x600 | 5.7 s (voxel table 1.2 s) |

All four build parts ran on one thread, and the build also blocks the GUI on
`show surface`. The "dots to surface points" loop tests every solvent dot
independently, so it now runs in fixed chunks under OpenMP and appends the
chunks in order.

| Workload | 4.5.11 | 4.5.12 |
| --- | --- | --- |
| surface build (1aon, +ray 16x16), three runs | 11.90, 11.97, 12.09 s | 9.48, 9.66, 9.66 s |

Same output: the POV-Ray export of the surface (every triangle) has the same
MD5 as the 4.5.11 build for seven cases (whole complex, chain A at
`surface_quality` 1, `surface_solvent`, `surface_type` 2 and 4,
`surface_cavity_mode` 1, partial selection) and for 1, 5 and 32 threads.
Check: `testing/tests/api/surface_threads.py`. Ray-traced PNGs are not a usable
comparison: their hashes differ between two runs of the same build.

Left serial: triangulation and refinement are advancing-front algorithms, and
`SolventDotNew` shares its counters across atoms. They are the next 6.9 s.

### Surface build again (4.7.7 to 4.7.9)

Temporary timers in `SurfaceJobRun`, `TrianglePointsToSurface` and
`MapSetupExpress`, 1aon, default settings, 32 threads. One change per commit.
Rows marked "about" were not timed again on 4.7.9: they are the 4.7.6 time
minus the measured saving of the parts inside them.

| Part of the surface build | 4.7.6 | 4.7.9 |
| --- | --- | --- |
| Solvent dots, per-atom pass | 1.07 s | 0.12 s |
| Solvent dots, cavity cull | 0.31 s | 0.31 s |
| Dots to surface points | 0.46 s | 0.46 s |
| Eliminate close points (3 passes) | 0.87 s | about 0.6 s |
| Eliminate troublesome points (4 passes) | 1.17 s | about 0.8 s |
| Triangulation | 3.55 s | about 2.4 s |
| ... of it: search for the next start point (132 searches) | 1.17 s | 0.18 s |
| Neighbor tables, 10 of about 450 000 points (inside the rows above) | 1.8 s | 0.6 s |

| Workload (median of 5) | 4.7.6 | 4.7.7 | 4.7.8 | 4.7.9 |
| --- | --- | --- | --- | --- |
| surface build (1aon, +ray 16x16) | 9.44 s | 8.63 s | 7.31 s | 6.24 s |

- 4.7.7: the triangulation's search for a start point walked the neighbors
  of vertices that were already triangulated and rejected them one neighbor
  at a time. The test now runs once per vertex.
- 4.7.8: the per-atom solvent-dot pass runs in chunks under OpenMP and
  appends the chunks in atom order.
- 4.7.9: `MapSetupExpress` (the neighbor table behind every distance search)
  counts first, then fills a vector of the exact size in parallel. Above
  10 000 points only. Per table 0.16-0.19 s to 0.04-0.07 s.

Same output: the POV-Ray export of the surface has the same MD5 as the 4.7.6
build for 13 cases (1aon whole, chain A, partial selection; 1rx1;
`surface_quality` -1, 1, 2; `surface_type` 2, 3, 4; `surface_solvent`;
`surface_cavity_mode` 1, 2; `surface_mode` 1) and for 1, 5 and 32 threads. The
old and the new `MapSetupExpress` ran side by side on every call of the
Python suite: identical tables. Checks: `testing/tests/api/surface_threads.py`,
`layerCTest/Test_Map.cpp`.

Full table, commit 3e4218d9, PyMOL 4.7.9, same machine, Python 3.13.12:

| Workload | 4.7.6 | 4.7.9 |
| --- | --- | --- |
| load 1aon.pdb.gz (58870 atoms) | 0.147 s | 0.191 s |
| ray 800x600, cartoon (1aon) | 4.286 s | 4.565 s |
| surface build (1aon, +ray 16x16) | 9.438 s | 6.510 s |
| ray 800x600, surface (1aon) | 13.065 s | 8.954 s |
| isosurface emd_1155.ccp4 level 2.0 (+ray 64x64) | 1.072 s | 1.046 s |
| isosurface level change 2.0 -> 3.0 (+ray 64x64) | 0.090 s | 0.093 s |
| save PSE (1aon, cartoon) | 0.180 s | 0.211 s |
| load PSE (1aon, cartoon) | 0.151 s | 0.175 s |
| startup `pymolish -ckq -d quit` | 0.073 s | 0.100 s |

The rows that got slower are machine noise, not the changes: the 4.7.9 run
came after an hour of builds and benchmarks on a laptop, and startup, which
none of the changes touch, moved the most. Load was checked on its own
because it builds a neighbor table: 20 loads, three runs each, old
`MapSetupExpress` 216 / 155 / 161 ms, new 225 / 160 / 157 ms.

On one thread (`OMP_NUM_THREADS=1` and `max_threads` 1, which is also what
a build without OpenMP does, the macOS default): surface build 12.36 s to
10.60 s, load 150 ms to 137 ms. The changes cost nothing where there are no
extra cores.

Left: the two advancing-front passes of the triangulation (0.85 s and
0.82 s), the point elimination passes (order-dependent, serial), the ray
tracer's own voxel table (1.0 s per ray of this surface).

### Ray tracer's voxel table (4.13.3)

Where a ray of the finished 1aon surface goes (800x600, 32 threads, 905 126
triangles; the ray tracer's own feedback plus temporary timers in
`BasisMakeMap`): collecting primitives 0.3 s, voxel tables 0.9 s, tracing
3.2 s. The two tables are built side by side, one thread each: the light's
takes 0.3 s, the camera's 0.75 s, of which 0.68 s is `MapSetupExpressPerp`
(62.5 million rows, pushed one by one into a vector that grows and is then
shrunk). It now counts first and fills a vector of the exact size in
parallel, the 4.7.9 pattern.

| | 4.13.2 | 4.13.3 |
| --- | --- | --- |
| Camera table, flags of voxel-spanning triangles included | 0.68 s | 0.28 s |
| ... with one OpenMP thread (`OMP_NUM_THREADS=1`) | 0.68 s | 0.38 s |
| Voxel phase of the ray (both tables) | 0.9 to 1.1 s | 0.5 s |
| ray 800x600 of the built surface, median of 5 (fastest) | 4.73 s (4.26 s) | 4.21 s (4.06 s) |

Same output: old and new function ran side by side and their three arrays
were compared on every call of the Python suite (209 calls) and of 10 rays
(1aon surface, cartoon, spheres, sticks with lines and labels,
`ray_trace_mode` 1, no shadows, a close-up; 1rx1 transparent surface, wide
field of view; an isosurface): identical. Check: `layerCTest/Test_Map.cpp`.

No table of `bench.py` this time: a second pymolish was running on the
machine and single rows varied 4x between repeats (cartoon ray 4.3 to 17 s
in one run). Cartoon ray and isosurface were timed alone instead, with
OpenMP on and off: the quiet repeats are the same (3.7 to 4.5 s, 0.95 s).

Left in the voxel phase: 0.17 s of the camera table's 0.28 s is the
zero-fill of the 250 MB list before the parallel fill (`std::vector` cannot
skip it; 115 users of `EList` make a type change a slice of its own), and
the light's table (`MapSetupExpressXY`, 0.2 s) now sets the time of the
phase. Tracing itself is 3.2 s of the 4.2 s: 28 s on one thread, so 32
threads give 9x.

### GUI frame time (first measurement, 4.7.9)

Command: `pymolish -qk $PWD/testing/bench_gui.py`. 1aon, viewport
952x747 px, 300 frames per row while rotating. A frame is one PyMOL draw
plus `glFinish`, outside Qt's paint cycle, so display synchronization is not
in the numbers. "First frame" includes building the representation.

Under WSLg the GUI gets Mesa's software renderer unless told otherwise. The
same Mesa has a D3D12 driver that reaches both GPUs of this laptop:

| Median frame | llvmpipe (default) | `GALLIUM_DRIVER=d3d12` (AMD Radeon 610M) | + `MESA_D3D12_DEFAULT_ADAPTER_NAME=NVIDIA` (RTX 5070) |
| --- | --- | --- | --- |
| lines | 55.1 ms | 0.9 ms | 2.3 ms |
| sticks | 51.7 ms | 23.1 ms | 11.3 ms |
| spheres | 24.6 ms | 6.3 ms | 0.8 ms |
| cartoon | 50.6 ms | 11.2 ms | 0.7 ms |
| surface | 122.7 ms | 27.1 ms | 10.4 ms |

| RTX 5070 | first frame | median | p95 | p99 | max | draw calls (median) |
| --- | --- | --- | --- | --- | --- | --- |
| lines | 0.11 s | 2.3 ms | 2.5 ms | 3.0 ms | 3.1 ms | 0.0 ms |
| sticks | 0.18 s | 11.3 ms | 13.7 ms | 36.3 ms | 85.8 ms | 6.0 ms |
| spheres | 0.08 s | 0.8 ms | 0.9 ms | 2.1 ms | 2.4 ms | 0.0 ms |
| cartoon | 0.55 s | 0.7 ms | 0.8 ms | 0.9 ms | 1.4 ms | 0.0 ms |
| surface | 5.75 s | 10.4 ms | 13.6 ms | 15.3 ms | 19.2 ms | 7.6 ms |

| llvmpipe | first frame | median | p95 | p99 | max | draw calls (median) |
| --- | --- | --- | --- | --- | --- | --- |
| lines | 0.26 s | 55.1 ms | 57.4 ms | 67.8 ms | 102.5 ms | 53.0 ms |
| sticks | 0.15 s | 51.7 ms | 54.4 ms | 55.9 ms | 58.2 ms | 49.8 ms |
| spheres | 0.11 s | 24.6 ms | 28.9 ms | 39.2 ms | 44.5 ms | 13.4 ms |
| cartoon | 0.78 s | 50.6 ms | 55.5 ms | 59.5 ms | 63.9 ms | 47.3 ms |
| surface | 6.00 s | 122.7 ms | 136.7 ms | 146.2 ms | 155.5 ms | 121.5 ms |

Idle with a cartoon shown, 5 s: 0.8 % of one core, 48 timer wakeups per
second (the 20 ms poll in `PyMOLGLWidget._pymolProcess`), no repaints.

What the numbers say:
- On this machine the renderer choice is worth 5 to 70 times; nothing in the
  drawing code comes close. The app does not select the D3D12 driver itself.
  The same scene (cartoon, sticks, spheres, lines, transparent surface)
  differs from the software image in about 1 % of the pixels, and a click
  selects the same residue on all three renderers. One visible difference:
  the D3D12 driver draws every line 1 pixel wide. `set trilines` restores
  the width for lines, nonbonded and dashes (3.1 ms instead of 2.3 ms per
  frame); mesh stays thin.
- On the GPU, sticks are the slowest representation and the only one with a
  long tail (p99 36 ms, max 86 ms); sticks and surface spend 6 to 8 ms per
  frame on the CPU before the GPU starts. These are the next two profiles.
- Idle cost is small; the 20 ms poll is not worth replacing for speed.

### Isosurface at 1σ: not the mesh

`isosurface` of `emd_1155.ccp4` at 1σ plus `ray 64, 64`: 10.8 s and 7.7 GB
peak, 0.7 GB afterwards. Building the object takes 0.00 s; a second `ray`
with the mesh already built takes 8 s again. Ray feedback: 3 221 521
triangles, 2.5 s converting them to primitives, 6.7 s filling the ray
tracer's voxel table (`MapSetupExpressPerp` lists every vertex in up to 27
neighbor cells, about 62 M rows for the 0.9 M-triangle protein surface). So
the cost belongs to CPU ray tracing of multi-million-triangle scenes, not to
marching cubes or mesh storage, and the OpenGL display path does not pay it.
Not optimized: a 1σ contour of a whole EM map is noise, 2σ takes 1.3 s.


### Bundle size with rdkit (4.11.3)

Linux bundle built on this machine with `packaging/pyinstaller.spec`,
2026-11-14, same tree, without and with the `chem` extra installed:

| | unpacked | tar.gz |
| --- | --- | --- |
| without rdkit | 287 MB | 109 MB |
| with rdkit 2026.9.2 (and the Pillow it needs) | 344 MB | 129 MB |

| with rdkit, pdbfixer 1.12.0 and openmm 8.6.1 (4.12.2) | 419 MB | 142 MB |

rdkit adds 57 MB unpacked and 20 MB to the download; pdbfixer and openmm add
75 MB and 13 MB. 64 MB of that is the `openmm` package with its force-field
files; the OpenCL plugins of the wheel are left out. `fix_structure` runs in
the bundle with the Reference and CPU platforms. `load_smiles` runs in
the bundle with no change to the spec except the two license texts. rdkit
2026.9.2 has wheels for Python 3.13 and 3.14.

### minimize (4.12.3)

`fix_structure`, then `minimize` with 500 steps, CPU platform, 32 threads,
2026-11-14. Amber14 with GBn2 implicit solvent, heavy atoms restrained.

| Protein | Atoms with hydrogens | No cutoff | 2.0 nm cutoff | 1.0 nm cutoff |
| --- | --- | --- | --- | --- |
| 12CA, 260 residues | 4009 | 126 to 272 s | 52 s | 4.3 to 5.9 s |
| AlphaFold P04637, 393 residues | 6220 | 397 s | not run | 7.3 s |

The heavy-atom RMSD to the start is the same with and without the cutoff
(0.14 A for 12CA), so the command uses 1.0 nm. `fix_structure` takes 0.7 s
on both. One force evaluation of 12CA: 1.2 s without a cutoff, 0.05 s with.

### Batch load (4.11.1)

`load_files` on a folder of 502 small PDB files (24 atoms each, 2 of them
unreadable), measured 2026-11-14: 0.30 to 0.57 s without a window, 0.42 s in
the window under WSLg. A plain `cmd.load` loop over the same files takes
0.31 s, so the command adds nothing measurable. Under the one-second limit
of step 12b: the load loop stays on the GUI thread.
