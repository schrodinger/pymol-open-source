# 1 — Performance and GPU

[Roadmap](../../PLAN.md) · Draft; no timings measured in this review

**Improve the operation a profile identifies.** Interactive graphics already use
OpenGL. Loading, geometry generation, ray tracing, callbacks and GPU submission
are different workloads and need separate measurements.

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
