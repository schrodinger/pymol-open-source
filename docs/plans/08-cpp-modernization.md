# 8 — Modern C++ in the core

[Roadmap](../../PLAN.md) · Request 2026-11-13 · Draft; counts are `grep` counts,
nothing was built, timed or changed for this plan

**Modernize the C-style C++ of the core one idiom and one file at a time, behind
tests and sanitizers. No from-scratch rewrite. Vendored C stays C.**

## What exists

- **The core is already C++.** `layer0/` to `layer5/` and `layerGraphics/` hold
  165 `.cpp` files, 219,000 lines, and no `.c` file. MSVC builds with
  `/std:c++17`; GCC and Clang use their default. The work is "C written in
  C++", not a language port.
- **Real C is vendored only:** `contrib/champ` (10 files, 7,500 lines) and the
  VMD molfile plugins in `contrib/uiuc` (30 files, 23,500 lines).
- **Upstream started the same migration and left it half done:**

  | Old idiom | Sites | Modern form already in the tree | Sites |
  | --- | --- | --- | --- |
  | `VLAlloc`, `VLACheck`, ... macros | 902 | `pymol::vla<T>` | 222 |
  | `Alloc`, `Calloc`, `FreeP`, `malloc`, `free` | 714 | `std::vector` | 676 |
  | raw `new` / `delete` | 589 | `std::unique_ptr` | 102 |
  | `sprintf`, `strcpy`, `strcat`, `strncpy` | 517 | `std::string` | 465 |
  | fixed buffers (`OrthoLineType`, `WordType`, ...) | 221 | | |
  | `ok_assert` / `ok_raise` chains, `goto` | 246 | `pymol::Result` | 419 |
  | C pointer casts | 572 | | |
  | `#define` (includes header guards) | 1,917 | | |

  About 3,800 sites plus the macros.
- **Safety net today:** 167 Python test files, 23 Catch2 files,
  `testing/bench.py`. No sanitizer build. Test coverage of the core was never
  measured.

## Decisions

- **No rewrite from scratch.** 219,000 lines whose behavior is defined by the
  code itself, not by a specification; a rewrite loses function by
  construction. Small behavior-neutral steps do not.
- **Vendored C is not touched.** The molfile plugins track VMD; rewriting them
  gains nothing and makes every plugin update a manual port. Champ: same.
- **The CPython C API code stays as it is** (2,479 `PyObject` uses, mostly
  `layer4/Cmd.cpp`). The API is C. pybind11 would be a new dependency for no
  user-visible gain.
- **One idiom per PR, behavior-neutral,** commit type `simplify(<scope>)`. A
  modernization PR never changes an algorithm, and a speed PR never sweeps
  idioms. A regression then has one cause.
- **Speed stays step 1's job.** Modernization does not make code faster by
  itself and can make it slower (`std::vector` zero-fills where `malloc` did
  not; `std::string` allocates inside loops). The two steps share an order
  instead: files are modernized in profile order, so a hot function is cleaned
  just before or after it is optimized, and `bench.py` gates both.
- **Delete first.** A line removed under step 2 is a line not modernized.
- **Session layout is frozen.** Structs written by `pse_binary_dump`
  (`layer2/AtomInfoHistory.h`) keep their layout; PSE files must load and save
  as before.
- **No public change.** `cmd.*`, settings and file formats are untouched, so
  `ChangeLog` needs no entry.

## Order

| Slice | Work | Why here |
| --- | --- | --- |
| 8a | Safety net: one Linux CI job with `-fsanitize=address,undefined` running the Python suite and Catch2; one `gcov` run that lists core files the suite never enters | Without it "no lost function" cannot be shown. A file with no coverage gets a test before it is touched. |
| 8b | `sprintf` / `strcpy` / fixed buffers to `std::string` and `snprintf` | The only idiom that is a memory-safety bug class. Mechanical. |
| 8c | `Alloc` / `FreeP` / `malloc` and raw `new` / `delete` to `std::vector` and `std::unique_ptr` | Removes manual cleanup paths, which makes 8e possible. |
| 8d | VLA macros to `pymol::vla<T>` or `std::vector`; delete the macro layer when its last user is gone | Finishes upstream's migration; ends in a deletion. |
| 8e | `ok` chains and `goto` cleanup to `pymol::Result` and early return | Needs 8c: early return is only safe once ownership is automatic. |
| 8f | Function-like macros to `constexpr` / inline functions; C casts to named casts | Lowest value; do it only in files already open for another slice. |
| 8g | `ov/` (3,700 lines of hand-written hash maps and string interning) to `std::unordered_map`, if `bench.py` is not slower | The one slice that deletes a directory. Measure before deciding. |

Inside a slice, go file by file in profile order: `Executive.cpp`,
`Selector.cpp`, `ObjectMolecule.cpp`, `CGO.cpp`, `Ray.cpp`, `RepSurface.cpp`
first.

## Risk

- **Tests prove only what they cover.** Rendering and GUI paths do not run
  headless. Slices that touch `layer1/` need the real GUI check, and 8a's
  coverage list decides where a test comes first.
- **Upstream merges.** Every modernized file stops taking upstream patches
  cleanly. No upstream remote is configured today; if merges from
  `schrodinger/pymol-open-source` are still wanted, restrict the work to files
  this fork already changed heavily.
- **Size.** This is not one step finished in a week. Each slice is many PRs.
  No date or line-count target is promised.

## Skipped, and when to add

- **C++20 / 23** (`std::span`, `std::format`): add when a slice needs one and
  the three CI compilers support it. It is a one-line build flag.
- **Splitting the large files** (`Executive.cpp` is 17,700 lines) or
  redesigning the object class hierarchy: only when a change inside is blocked
  by it.
- **`clang-tidy` `modernize-*` automatic fixes:** decide in 8a. It is a new
  pinned dev tool, and whole-tree automatic fixes break the one-idiom,
  reviewable-diff rule; per-file use may pay off in 8f.

## Checks

Every PR: the full Python suite, Catch2, the sanitizer job, `bench.py` within
noise of the baseline in [docs/baseline.md](../baseline.md), and a PSE load and
save round trip. `layer1/` changes: real GUI.

## Status

**8a done 2026-11-13 in 4.7.2.**

- CI job `sanitize-Linux` builds with
  `-fsanitize=address,undefined -fno-sanitize-recover=undefined` and runs the
  Python suite; any report fails the job. Run it locally with the same
  `CXXFLAGS` / `CFLAGS` / `LDFLAGS` and
  `LD_PRELOAD=$(gcc -print-file-name=libasan.so) ASAN_OPTIONS=detect_leaks=0`.
- The first run found three bugs, fixed in `c2de3087`: `free()` of an object
  made with `new` (`PyMOL_Free`), left shifts of negative ints in the sculpting
  hash, and four uninitialized `bool` in `ObjectCGOState`.
- Left open on purpose: the alignment check is off. CGO keeps structs that
  hold pointers inside its `float` buffer, so every access is misaligned.
  Harmless on x86-64 and arm64; fixing it is a CGO storage change, not a
  sweep. Leak detection is off because Python does not free everything at
  exit.
- Coverage of the headless suite (`--coverage` build, `gcov`): **about half of
  the core**, 96,000 of 191,600 executable lines (68,400 of them are the
  stroke-font table in `VFont.cpp`).

  | Never entered by the suite | Lines |
  | --- | --- |
  | `layer1/CGOGL.cpp`, `SceneRender.cpp`, `ScenePicking.cpp`, `Picking.cpp` (OpenGL drawing and picking) | 2,682 |
  | `layer1/SceneMouse.cpp`, `layer3/Interactions.cpp`, `layer1/ScrollBar.cpp` (mouse) | 2,130 |
  | `layer4/PopUp.cpp`, `Menu.cpp` (viewport menus) | 696 |
  | `layer0/GenericBuffer.cpp`, `ShaderPrg.cpp`, `PostProcess.cpp`, `layerGraphics/gl/GLVertexBuffer.cpp` (GPU buffers, shaders) | 890 |
  | `layer2/ObjectCurve.cpp`, `RepSphereImmediate.cpp`, `RepSphereGenerate.cpp` | 520 |
  | `ov/src/OVOneToAny.cpp`, `OVRandom.cpp` | 300 |

  | Large and under 25 % | Covered |
  | --- | --- |
  | `layer1/CGO.cpp` (6,343 lines) | 17 % |
  | `layer5/PyMOL.cpp` (1,481) | 24 % |
  | `layer2/RepLabel.cpp` (1,296) | 15 % |
  | `layer3/Seeker.cpp` (1,181) | 7 % |
  | `layer0/ShaderMgr.cpp` (988) | 6 % |
  | `layer2/ObjectSlice.cpp` (752) | 10 % |
  | `layer1/Control.cpp`, `Seq.cpp`, `ButMode.cpp` (1,544) | 9 to 21 % |

  For scale: `Executive.cpp` 55 %, `ObjectMolecule.cpp` 70 %, `Selector.cpp`
  80 %, `Ray.cpp` 59 %, `Cmd.cpp` 76 %, `Scene.cpp` 29 %, `RepSurface.cpp` 32 %.

  Rule that follows: in the files above only compiler-checked mechanical
  changes are allowed without a new test or a real GUI check. Hand-written
  changes there wait for one.
- Not done: the Catch2 tests are compiled in CI but no step runs them
  (`pymolish --testing`); Catch2 is not installed on the development machine,
  so adding the step could not be tried.

**8b in progress.**

- 4.7.3: 227 of 251 `sprintf` are `snprintf(buf, std::size(buf), ...)`.
  `std::size` does not compile on a pointer, so the compiler proved every
  destination is an array. Left for a hand pass: 17 calls whose destination is
  a pointer parameter (the callee does not know the size) and 7 that write at
  an offset into a larger buffer.
- 4.7.4: 159 of 197 `strcpy` are `UtilNCopy(dst, src, std::size(dst))`, the
  bounded copy the tree already had. Left for a hand pass: 28 with a pointer
  destination and 10 with a computed one.
- 4.7.5: 22 of 27 `strcat` are `UtilNConcat(dst, src, std::size(dst))`. Left:
  5 with a pointer destination.
- 4.7.6, hand pass over what the sweeps left. Three real defects, with a test
  (`testing/tests/api/long_strings.py`) that fails under the sanitizer job
  without the fixes: an unterminated `strncpy` in `SettingGetTextValue`
  (reachable with `set` of a long text value), an off-by-one in the C label
  parser (`cmd.label2`), and unbounded `sprintf` of the file name in the PDB
  reader's load message (C API only).
- Open in 8b: about 60 `sprintf` / `strcpy` / `strcat` calls still write
  through a pointer. Read one by one, each is bounded by its own logic today
  (constant text, a name of known maximum length, or a size check just
  before), so they are not defects; they go away when the function gets a
  modern signature (`std::string` return or a size argument), which is 8c/8e
  work in the same files. The other 20 `strncpy` were read: each ends up
  terminated or needs no terminator, though several rely on the destination
  struct having been zero-filled. Then the fixed buffer types themselves
  (`OrthoLineType`, `WordType`, ...).

## Review 2026-11-14: is it finished, and what is next

`grep` counts and CI configuration only; nothing was built or changed.

**Not finished.** 8a is done, 8b is mostly done, 8c to 8g are not started.

| Idiom | Left today | Slice |
| --- | --- | --- |
| `sprintf` / `strcpy` / `strcat` through a pointer | 31 / 44 / 5 | 8b, closes with 8c and 8e |
| `Alloc` / `FreeP` / `malloc` / `free` | about 420 | 8c |
| raw `new` / `delete` | about 500 / 110 | 8c |
| VLA macros | about 700 | 8d |
| `ok_assert` / `ok_raise` / `goto`, `CHECKOK` | about 190, 260 | 8e |
| `ov/` users | about 110 | 8g |

`NULL` is already gone (27 left against 4,000 `nullptr`).

### Order from here: the safest work first

The working version is the constraint, so the next slices are the ones the
compiler or a byte comparison can prove neutral. Hand-written changes come
after them.

| Slice | Work | Proof that nothing changed |
| --- | --- | --- |
| 8h | Run the Catch2 tests in CI: add `pymolish --testing` to the Linux job. They are compiled today and never run. | It is a new check, not a change. |
| 8i | Delete `#ifdef` branches for macros nothing defines. Neither `setup.py` nor any source defines `PURE_OPENGL_ES_2` (211 guards), `_PYMOL_NOPY` (145), `_WEBGL` (66), `_PYMOL_NO_RAY`, `_PYMOL_LIB`, `PYMOL_EVAL`, `WITH_UNUSED_FUNCTIONS`, `_PYMOL_EMBEDDED`, `PICKLETOOLS`, `_PYMOL_ACTIVEX`, `_PYMOL_INCENTIVE`, `_MACPYMOL_XCODE`. About 460 simple guards and at least 1,100 dead lines; compound conditions were not counted. One macro per PR, with `unifdef` (a tool on the development machine, not a dependency; not installed yet). | The compiler sees the same text: `g++ -E -P` output of every touched file is identical before and after, apart from `__LINE__` values. Then the suite. |
| 8c to 8g | As planned above, file by file in profile order. | Suite, sanitizer job, `bench.py`, PSE round trip. |

8i comes before 8c because it removes lines that 8c to 8e would otherwise
modernize, and because it is the only slice with a mechanical proof.
`os_gl_glut.h` and `os_gl_glut_pretend.h` (GLUT is gone since 4.3.1) go in the
same slice if their constants have no other user.

### New mechanisms and libraries, judged

CI compilers: GCC 11 (Ubuntu 22.04), Apple Clang, MSVC. Cold build is 54 s
([baseline](../baseline.md)).

| Idea | Verdict |
| --- | --- |
| C++20 flag | **Later.** All three compilers have it. Today it would replace only `include/pymol/span.h` (58 lines, 21 uses). Turn it on when a slice needs `std::span` or designated initializers; it is one line in `CMakeLists.txt` and `setup.py`. |
| `std::format` or the fmt library | **Skip.** GCC 11 has no `<format>`; fmt is a new dependency. `snprintf` is bounded since 8b. |
| `std::expected` (C++23) | **Skip.** `pymol::Result` does the job at 419 sites; GCC 11 lacks it. |
| pybind11 | **Skip**, decided above. |
| Eigen | **Skip.** GLM is linked and used at 600 sites. `layer0/Vector.cpp` and `Matrix.cpp` (3,500 lines of hand-written math on `float*`) move to GLM only inside functions step 1 is optimizing anyway. |
| Drop `include/tnt` (10,800 vendored lines) | **Candidate, not neutral.** Two files use it for a 3x3 eigen solve, SVD and LU (`Matrix.cpp`, `ccealignmodule.cpp`). A replacement changes results in the last digits, so it needs an RMSD tolerance test for `cealign` and `fit` first. After 8g. |
| `std::from_chars` or fast_float for number parsing | **Only on a profile.** 315 `atof` / `sscanf` sites; float `from_chars` is recent in Apple's libc++. Belongs to step 1 if a loader profile shows parsing. |
| TBB or a thread pool | **Skip.** OpenMP is there (14 pragmas), `Ray.cpp` has its own workers. |
| `-Wall -Wextra` | **Cheap, once.** Add to the sanitizer job without `-Werror`, read the list once, fix what is a real defect (uninitialized values), ignore the rest. |
| `clang-tidy` | **Still open**, as above: per file in 8f, never whole-tree. Not installed on the development machine. |
| ccache, Ninja, unity build, precompiled headers | **Skip.** 54 s cold build is not a problem. |
| Fonts as `.ttf` data files instead of `FontTTF.h` (16.8 MB of source) | **Skip.** Smaller repository, but the file must then be found at run time on three platforms: a packaging risk for no user-visible gain. |

### Code structure

- **The layers are not layers.** In about 120 cases a file includes a header from a higher
  layer (`layer0` from `layer1`: 20 files; `layer2` from `layer5`: 26;
  `layer1` from `layer5`: 16). The directory names describe history, not
  dependencies.
- **Do not move or rename directories anyway.** A move changes no behavior and
  helps no user, makes every `git blame` and upstream comparison harder, and
  the roadmap already rejected collapsing directories without evidence. The
  repo map in `AGENTS.md` is the cheap fix for "where is what".
- **What does improve structure without risk is deletion:** 8i (dead
  branches), 8d (the VLA macro layer and `MemoryDebug`), 8g (`ov/`), then
  `include/tnt`. Each ends with files gone.
- **Splitting large files** (`Executive.cpp`, 17,700 lines) stays skipped
  until a change inside is blocked by the size. A split is a pure move, so
  when it is needed it gets its own PR with no other edit in it.
