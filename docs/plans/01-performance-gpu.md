# 1 — Performance and GPU

[Roadmap](../../PLAN.md) · 1a done, 1b in progress; timings are in the Status sections and in [docs/baseline.md](../baseline.md)

**Improve the operation a profile identifies.** Interactive graphics already use
OpenGL. Loading, geometry generation, ray tracing, callbacks and GPU submission
are different workloads and need separate measurements.

## Status 2026-11-08

1a done: headless baseline in `testing/bench.py`, numbers and first observations in [docs/baseline.md](../baseline.md). GUI frame latency not measured. 1b first change in 4.5.10: session load 0.92 s to 0.14 s (C unpickler with exact fallback), numbers in the baseline. Second change in 4.5.12: the surface build, which is 9.5 of the 15 s of "ray a surface", runs its dot-to-point loop in parallel (11.97 s to 9.66 s, identical geometry). The 1 sigma isosurface cost was traced to the ray tracer's voxel table, not the mesh, and left alone. Next candidates: triangulation (3.5 s) and refinement (1.9 s) of the surface build, `SolventDotNew` (1.5 s), GUI frame latency (never measured).

## Status 2026-11-14

Three more 1b changes, numbers in [docs/baseline.md](../baseline.md): the surface build of 1aon went from 9.44 s to 6.24 s with the same triangles (4.7.7 start-point search of the triangulation, 4.7.8 parallel solvent dots, 4.7.9 neighbor table sized first and filled in parallel). Ray tracing that surface at 800x600: 13.1 s to 9.0 s.

GUI frame time is now measured: `testing/bench_gui.py`. Main finding: under WSLg the GUI runs on Mesa's software renderer (50 to 120 ms per frame for 1aon); with `GALLIUM_DRIVER=d3d12` the same frames take 0.7 to 11 ms on the RTX 5070. Open decision: whether the app selects that driver itself on WSL. It draws lines 1 pixel wide, so the default would need `trilines` with it, and mesh lines stay thin.

Found while measuring, fixed in 4.7.10: atom labels were solid blocks in the GUI on every driver without `EXT_geometry_shader4` (all Mesa drivers). An upstream refactor of `RepLabel` had dropped the switch that turns shaders off for the fallback label path. Check: `testing/tests/api/label_render.py` (needs the GUI).

Next candidates: sticks on the GPU (11 ms per frame, p99 36 ms, 6 ms of it CPU), surface draw submission (7.6 ms CPU per frame), the two advancing-front passes of the triangulation (0.85 s and 0.82 s), the ray tracer's voxel table (1.0 s per ray of the 1aon surface; the 4.7.9 count-then-fill pattern may carry over to `MapSetupExpressXYVert`).

## Status 2026-11-16

1b change in 4.13.3, numbers in [docs/baseline.md](../baseline.md): the ray tracer's voxel table for the camera view (`MapSetupExpressPerp`) is counted first and filled in parallel, 0.68 s to 0.28 s for the 1aon surface; a ray of that surface 4.7 s to 4.2 s. Tracing is now 3.2 s of those 4.2 s (9x on 32 threads).

Next candidates: sticks on the GPU and surface draw submission (GUI, wait for a quiet machine and the driver decision), the two advancing-front passes of the triangulation (0.85 s and 0.82 s), thread scaling of the trace loop (9x on 16 cores / 32 threads), the zero-fill of the voxel list (0.17 s).

## Evidence and corrections

- Reuse [performance tests](../../testing/tests/performance), especially
  `pdbloading.py`, `shaders.py`, `iterate.py` and `async_builds.py`. Loading and
  shader stress checks have `no_run_all`; the shader test also requires a GUI.
- [setup.py](../../setup.py) defaults OpenMP on outside macOS.
  [Marching cubes](../../layer0/marching_cubes.cpp) already has parallel loops;
  [Ray.cpp](../../layer1/Ray.cpp) has `max_threads` and dedicated workers.
- PDB/mmCIF have native readers in
  [ObjectMolecule.cpp](../../layer2/ObjectMolecule.cpp),
  [ObjectMolecule2.cpp](../../layer2/ObjectMolecule2.cpp), and
  [CifMoleculeReader.cpp](../../layer2/CifMoleculeReader.cpp). Python's
  [importing.py](../../modules/pymol/importing.py) is not proof of Python parsing.
- Rendering spans [CGO](../../layer1/CGO.cpp), [Scene](../../layer1/Scene.cpp),
  and [shaders](../../data/shaders). GPU rendering can still be CPU-bound.
- Inspect existing [Qt worker utilities](../../modules/pymol/Qt/utils.py) and
  movie-export change `f4766ab1`. They do not establish that all `cmd` operations
  are safe on arbitrary threads.

## 1a — Baseline before cleanup or optimization

Use local fixtures and deterministic settings. Start with existing small/medium
structures; add a reproducibly generated large case only when needed. Keep fixture
preparation and downloads outside the measurement.

| Workload | Measure separately | Correctness check |
| --- | --- | --- |
| Startup and PDB/mmCIF load | Startup, I/O/decompression, parsing, bond/representation work; peak RSS | Atom/bond counts, coordinates, chains and states |
| Cartoon/sticks/surface rotation | First build, warm frame latency, CPU submission and GPU time where available | Picking, transparency, labels and reference images |
| Fixed-resolution CPU ray | End-to-end wall time and thread scaling | Image comparison with a justified tolerance |
| Map contour/isosurface | Build time, memory and repeated threshold changes | Same geometry/appearance for the same input and level |
| PSE save/reopen | Time, peak memory and file size | Objects, states, settings and views survive |
| GUI interaction | Idle CPU and latency during load, representation changes and export | Window remains usable; completion/errors are shown |

Record commit, release/debug flags, compiler, OS, CPU/GPU/driver, Python/Qt,
workload size/hash, viewport size, thread/shader settings, and warm/cold state.
Warm up, repeat at least five times, and report median/spread. Sample enough frames
for meaningful tail latency; do not infer it from five frame samples.

Compare optimized builds: Linux CI currently sets `DEBUG=1`. Record display
synchronization and GPU/software rendering; enqueue time is not completed GPU
execution. Separate profiling from final timing. Extend existing checks using
`perf_counter`, `cProfile`, or available native profilers before adding tools.

Put a result table and runnable command in each optimization PR. Until execution,
results remain **not measured**. Enforce timing budgets only after establishing a
repeatable reference machine; shared CI results may be too noisy.

## 1b — One measured bottleneck per change

1. **Algorithms and copies:** avoid rebuilding unchanged representations; find
   repeated model/coordinate copies, Python/C++ crossings and decompression.
   Cache only with a demonstrated cost and explicit invalidation rules.
2. **Existing CPU paths:** improve allocations/locality, then evaluate
   vectorization, OpenMP or worker counts. Check determinism, memory, nested
   parallelism, the GIL and PyMOL locks. Do not add pragmas to Ray by assumption.
3. **Existing GPU paths:** inspect draw calls, uploads, buffer reuse, shader
   variants and stalls. Try batching/instancing where measured, preserving
   selection and picking passes.
4. **Interaction:** coalesce redundant updates if profiling justifies it. Move
   suitable computation into existing workers/native paths; keep widgets and GL
   context operations on their owning threads. Native code can still block the
   GUI or hold the GIL.
5. **Remeasure:** keep gains exceeding observed noise without unacceptable memory
   or compatibility regressions. Reprofile before choosing another target.

Deleting fallback rendering or forcing `use_shaders` is a compatibility decision.
Audit overlays, text, CGOs, picking, stereo, image export and supported hardware
before proposing it; a Core profile is not a prerequisite for these optimizations.

## GPU compute gate

Apple documents native OpenGL support through 4.1 on Apple silicon, with OpenGL
deprecated. An OpenGL 4.3 compute requirement would exclude that native path.
Any experiment needs capability checks and a supported fallback; a new UI toolkit
does not fix this. [Apple WWDC technical presentation](https://developer.apple.com/videos/play/wwdc2020/10631/?time=1481)

Propose a bounded ray-tracing or contouring experiment only after profiling shows
that kernel dominates a real workload. A new compute framework requires the
repository's measured **5× end-to-end win**, including transfers, compilation,
synchronization and equivalent output quality. Report total application benefit,
hardware coverage and fallback/maintenance cost. A kernel-only win is insufficient.
No CUDA/OpenCL/Vulkan dependency is planned now.

## Additional ideas and completion checks

- Measure **idle wakeups** and redundant representation invalidation before
  pursuing another renderer; less work may also improve battery use.
- Track **peak memory and first useful frame**, not only total loading time.
- Profile **metadata refresh and selection synchronization** after step 3, so
  the panel cannot add expensive work to every render frame.
- Add cancellation at existing safe boundaries for long exports; do not promise
  cancellation for an indivisible native operation.

Each optimization needs a runnable regression check, before/after evidence and
relevant image/geometry or file round-trip checks. Rendering changes need real GUI
validation. Memory mapping, universal instancing, a new ray tracer and automatic
quality reduction remain deferred until measurements justify them.
