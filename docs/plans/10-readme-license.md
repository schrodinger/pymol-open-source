# 10 — README and license files

[Roadmap](../../PLAN.md) · Request 2026-11-14 · Done in 4.9.0, see
[Status](#status). Not legal advice.

**New README that says what this is: a fork of Open-Source PyMOL. PyMOL code
stays under its license, word for word; new code is MIT. Text first (10a),
images when [step 9](09-icons-logo.md) delivers them (10b).**

## What exists

| File | State |
| --- | --- |
| `README.md` | 30 lines, titled "Open-Source PyMOL", shows the PyMOL logo, sends readers to pymol.org |
| `README` | 281 lines of upstream text from 2003; says "Python license", which is wrong today |
| `LICENSE` | Upstream notice, unchanged |
| `pyproject.toml` | `license = {file = "LICENSE"}`; only author listed is Schrodinger |
| `packaging/pyinstaller.spec` | copies `LICENSE` into the bundle |
| `AGENTS.md` §5 | "Don't touch `LICENSE`" |
| `docs/` | `plugins.md`, `baseline.md`, `plans/`; plus `INSTALL`, `DEVELOPERS`, `PACKAGING`, `ChangeLog` at the root |

## What the PyMOL license requires

| Requirement (from `LICENSE`) | How this step meets it |
| --- | --- |
| Copyright notice "appears in all copies" | Upstream text stays byte for byte in `LICENSE`; it ships in the repository, the wheel and the bundle. |
| Copyright and permission notice "appear in supporting documentation" | README license section links `LICENSE`; `docs/licenses.md` links it again. |
| Name of Schrodinger not used "in advertising or publicity" | README names Schrodinger only as copyright holder and trademark owner. |
| Derivative "must be plainly distinguished" from Schrodinger's PyMOL | Title is `pymolish-studio`; second line says independent fork, not affiliated, official product at pymol.org. |
| Trademark use limited to four slogans plus footnote | README uses "Based on PyMOL(TM) technology." and ends with "PyMOL is a trademark of Schrodinger, LLC." |

The PyMOL license is permissive and compatible with MIT for added code.

## Wording

The license **requires** the notices and the distinction. It **permits**, and
does not require, four slogans; whoever uses one must add the footnote. Texts
are copied from `LICENSE` and used exactly, including "(TM)" and the spelling
"Schrodinger".

| Text | Status | Goes where |
| --- | --- | --- |
| Upstream copyright and permission notice | required | `LICENSE` (repository, wheel, bundle) |
| "pymolish-studio is an independent fork. It is not a Schrodinger product and is not affiliated with or endorsed by Schrodinger, LLC." | required in substance (wording ours) | README line 2, About box, `docs/licenses.md`, release notes |
| "pymolish-studio. Based on PyMOL(TM) technology." | allowed slogan, used with the product name | README under the title, About box, release notes |
| "pymolish-studio. Contains PyMOL(TM) source code." | allowed slogan, used with the product name | README license section, `docs/licenses.md`, `pyproject.toml` description |
| "PyMOL is a trademark of Schrodinger, LLC." | required wherever a slogan appears | last line of README and `docs/licenses.md`, About box, release notes |
| "Includes PyMOL(TM).", "Built using PyMOL(TM)." | allowed, not used | two slogans are enough |

**Product name with a slogan.** The product name comes first and ends with a
period; the slogan follows as its own sentence, unchanged. Not
"pymolish-studio is based on PyMOL(TM) technology": the license permits the
four slogans as written, so they are not reworded or built into a sentence.

### Trademarks

| Mark | Owner | Rule in our documents |
| --- | --- | --- |
| PyMOL | Schrodinger, LLC | **Required by `LICENSE`.** Written "PyMOL(TM)" in every slogan and at the first mention in each user-facing document (README, `docs/features.md`, `docs/licenses.md`, About box, release notes). The footnote "PyMOL is a trademark of Schrodinger, LLC." ends each of them. |
| Schrodinger | Schrodinger, LLC | Named only as copyright holder and trademark owner. Never as author of the fork, sponsor or contact. |
| pymolish-studio, pymolish | none registered | No "(TM)" or "(R)" on our own name. Claiming a mark that contains "PyMOL" would make the name question below worse. |
| Qt, Python, Windows, macOS, Linux, GitHub, AlphaFold and other third-party names | their owners | No license in this project requires a notice. One line at the end of `docs/licenses.md`: "All other names are trademarks of their respective owners." |

Outside these lines, "PyMOL" appears only as plain reference ("fork of
Open-Source PyMOL", file formats, the upstream link), never in a logo or a
claim of endorsement.

**What the slogans do not cover.** Putting the product name next to a slogan
shows the two are different products. It does not grant the right to a
product name that contains "PyMOL": the license says every other use of the
mark needs "specific, written prior permission" from Schrodinger. Options,
user's call: keep the name and accept the risk, ask Schrodinger in writing,
or rename. This plan keeps the name and changes nothing else.

Places that need the lines and are easy to miss:

- `.github/workflows/build.yml:243`: the `--notes` text of each release.
- `pyproject.toml:7`: the description still reads "PyMOL is a ...".
- The About box and splash, handled in [step 9](09-icons-logo.md).

## Decisions

Approved by the user on 2026-11-14: the license layout below and the
`AGENTS.md` change (already made, with the dependency decisions).

- **Two license files.**
  - `LICENSE`: a short header (what is under which license, links to both
    files and to upstream), then the upstream notice unchanged.
  - `LICENSE-MIT`: MIT, "Copyright (c) 2026 Adrian Jasinski and pymolish-studio
    contributors".
- **Honest boundary.** New files are MIT. Changes inside upstream files are
  MIT too, but those files still carry the PyMOL notice. A distribution needs
  both files. The header says exactly that. No per-file headers are added;
  existing headers are not touched.
- **`AGENTS.md` §5, done:** the Schrodinger notice inside `LICENSE` stays byte
  for byte; text may be added above it. Adding a header does not change or
  remove a notice, so the license's own condition holds.
- **Third-party list:** `docs/licenses.md`, one table: component, license,
  where its text is, link. Read each license from its source, not from memory.
  Rows to cover: PyMOL, new code, Lucide, champ, molfile plugins, pocketfft,
  embedded fonts, and what bundles ship (PySide6/Qt is LGPL-3.0: the bundle
  must carry that text), numpy, GLEW, GLM, libpng, freetype, msgpack-c,
  netcdf, nlohmann-json. Steps 12 and 13 add rdkit, OpenMM, PDBFixer.
- **Old `README`:** `git mv README docs/README-upstream.txt`, content
  unchanged. Keeps its notice, ends the two-README confusion.
  `grep -rn "README'" setup.py _custom_build packaging` first.
- **`pyproject.toml`:** `license-files = ["LICENSE", "LICENSE-MIT"]` if the
  pinned setuptools supports it; add the fork maintainer to `authors`, keep
  Schrodinger. Spec copies `LICENSE-MIT` too.

## README outline (target: under 80 lines)

1. Logo (10b), name, one line: what it is. Directly under the title:
   "pymolish-studio. Based on PyMOL(TM) technology."
2. "Fork of [Open-Source PyMOL](https://github.com/schrodinger/pymol-open-source).
   Independent; not a Schrodinger product, not affiliated with or endorsed by
   Schrodinger, LLC. Official PyMOL: pymol.org."
3. CI badge, one screenshot (10b).
4. Features: 6 to 8 bullets, each one line, each linking into `docs/`.
5. Download (Releases), install from source (link `INSTALL`).
6. Documentation table: `docs/features.md`, `docs/plugins.md`,
   `docs/licenses.md`, `PLAN.md`, `ChangeLog`.
7. License: "pymolish-studio. Contains PyMOL(TM) source code."; PyMOL code →
   [`LICENSE`](../../LICENSE); new code → `LICENSE-MIT`; third parties →
   `docs/licenses.md`; upstream repository and its license.
8. Last line: "PyMOL is a trademark of Schrodinger, LLC."

New doc: only `docs/features.md` (removed, added, changed versus upstream;
source: `PLAN-DONE.md` and `ChangeLog`). Steps 11 to 14 each add their section.

## Steps

- **10a, no images:** `LICENSE` header, `LICENSE-MIT`, README,
  `docs/features.md`, `docs/licenses.md`, move old `README`, `pyproject.toml`
  (license files, authors, description), spec, release notes text.
- **10b, after step 9:** logo at the top of the README (referenced from
  `data/pymol/icons/`, no second copy), one screenshot `docs/img/window.png`
  taken from pymolish itself.
- **Check:** one test in `testing/tests/system/`: every relative link in
  `README.md` and `docs/*.md` points at an existing file, `LICENSE` ends
  with the upstream text (compare against a stored SHA-256 of that tail), and
  every file that contains a slogan also contains the trademark footnote and
  has the product name directly before the slogan; no file writes
  "pymolish(TM)" or "pymolish-studio(TM)".
- PR diff of `LICENSE` must show zero removed lines.

Commit: `docs(license): new README, MIT for new code, PyMOL notice kept`.

## Status

Done 2026-11-14 in 4.9.0, 10a and 10b together. Checked:
`testing/tests/system/license_docs.py`; `git diff LICENSE` shows 21 added
lines and none removed; the wheel metadata lists both license files.

Differences from the plan above:

- **README shows the icon, not the logo.** `logo.png` has near-black text on
  a transparent background and does not show on a dark page. The icon reads
  on both. The logo also spells the name "PyMOLish": see the name question
  under [Wording](#wording) before using it in public.
- **`pyproject.toml`:** `license` is the SPDX expression
  `LicenseRef-Open-Source-PyMOL AND MIT` with `license-files`; this needs
  setuptools 77, so the build floor went from 69.2.0 to 77.0.3. Homepage and
  Documentation point at this repository, not at pymol.org.
- **Qt license text:** `packaging/licenses/` holds LGPL-3.0 and GPL-3.0
  (copied from Debian's `/usr/share/common-licenses`); the spec puts them in
  the bundle. PySide6 wheels do not carry the text.
- **More rows in `docs/licenses.md`** than planned: TNT/JAMA, CE alignment,
  US-align and the Gentium font were found in the tree.
- **`setup.py`** installs `LICENSE-MIT` next to `LICENSE`.

Open:

- **License line for the icon and logo:** not delivered;
  `docs/licenses.md` says so.
- **Notices of native libraries in bundles:** GLEW, libpng, FreeType,
  netCDF, HDF5, curl, zlib, zstd and Python are shipped in a bundle without
  their license texts. `docs/licenses.md` lists them as "not included".
  Collect them in the spec when bundles stop being previews.
- **FreeType credit line** in `docs/licenses.md` uses the years of the
  Ubuntu package (1996-2023); a bundle built from conda-forge may ship a
  newer FreeType.
- **Bundle and release notes** were not built here: first CI run.

### Follow-up, 4.9.1

Requested 2026-11-14: documents that a person can read without knowing the
project. The legal lines above are unchanged and still tested.

- README is shorter and starts from what the program does; it has a
  three-step "Try it" and a table "I want to ... / Read ...".
- New `docs/getting-started.md`: first structure, the window, the mouse,
  looks, saving. Its commands were run headless on 4.9.0.
- `docs/licenses.md` opens with a plain summary; `docs/features.md` says who
  it is for.
- `INSTALL` became `docs/install.md` with the commands that work today
  (`pymolish`, `--config-settings`, pixi on Windows and macOS).
- `DEVELOPERS` told contributors to transfer copyright to Schrodinger. It
  moved unchanged to `docs/DEVELOPERS-upstream.txt`; the new
  `CONTRIBUTING.md` says contributions are MIT.
- Not touched: `PACKAGING` (its `python setup.py install` line is outdated;
  the replacement was not tested), `docs/plugins.md`, `docs/baseline.md`.

## Skipped

- Rewriting `INSTALL`, `DEVELOPERS`, `PACKAGING`: fix only commands that are wrong.
- SPDX headers in every file, a `NOTICE` file, REUSE compliance: add if a
  packager asks.
- A documentation site: Markdown in `docs/` is enough.
- A list of every third-party trademark: the one general line covers them.

---

PyMOL is a trademark of Schrodinger, LLC.
