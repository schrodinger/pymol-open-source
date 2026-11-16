# 9 — Own icon and logo

[Roadmap](../../PLAN.md) · Request 2026-11-14 · Done in 4.8.0, see
[Status](#status). Not legal advice.

**Replace the one PyMOL image the application still shows, and stop the window
from calling itself PyMOL.**

## What exists

- **One PyMOL image:** `data/pymol/icons/icon2.svg`, the upstream logo. Two
  users: `make_pymol_qicon` (`modules/pmg_qt/pymol_qt_gui.py:1219`) and
  `README.md:5`.
- **Toolbar and menu icons are not PyMOL's:** 51 Lucide files (ISC) in
  `data/pmg_qt/icons/lucide/`. They stay.
- **Bundles have no icon:** `packaging/pyinstaller.spec` calls `EXE(...)`
  without `icon=`.
- **The window still presents itself as PyMOL:**

  | Where | Text now |
  | --- | --- |
  | `pymol_qt_gui.py:129`, `:1264` | window title `PyMOL` |
  | `pymol_qt_gui.py:924` `show_about` | "The PyMOL Molecular Graphics System", `https://pymol.org`, `sales@schrodinger.com` |
  | `layer1/Ortho.cpp:2578` `OrthoSplash` | "PyMOL(TM) Molecular Graphics System", copyright lines |

## Does the license allow it

`LICENSE` has two parts.

| Clause | Effect on this step |
| --- | --- |
| Permission to "modify, distribute, and distribute modified versions" | Replacing images is allowed. |
| "don't change or remove any of the Copyright notices" | An icon is not a notice. The copyright lines in the About box and the splash text stay word for word. |
| Derivative software "must be plainly distinguished from any and all PyMOL products distributed by Schrodinger" | The license asks for this change. The PyMOL logo, the title `PyMOL` and the Schrodinger sales address are what break it today. |
| Only four slogans may use the mark without written permission | About box uses "Based on PyMOL(TM) technology.", the not-affiliated line and the footnote "PyMOL is a trademark of Schrodinger, LLC."; exact texts in [step 10, Wording](10-readme-license.md#wording). |

Open risk, user's call: the product name "pymolish" contains the mark. The
license does not cover product names. This plan does not rename anything.

## Needed from the user

| File | Form |
| --- | --- |
| Application icon | PNG, square, 1024 x 1024, transparent background, readable at 16 px |
| Logo | PNG, readable on light and dark pages (or two files) |
| Rights | One line: who made them, under which license they ship. Must not resemble the PyMOL logo. |

## Decisions

- **One master PNG, no size set.** `data/pymol/icons/pymolish.png` and
  `logo.png`. `QIcon` scales it. Ceiling: blurry at 16 px; then add one
  hand-tuned 32 px file.
- **Bundle icon from the same PNG.** `icon='../data/pymol/icons/pymolish.png'`
  in the spec. PyInstaller converts to `.ico` when Pillow is installed in the
  build job: one line in `.github/workflows/build.yml`, no committed `.ico` or
  `.icns`. Verify on the Windows artifact.
- **macOS:** the spec builds a folder, not an `.app`, so the Dock icon comes
  from `app.setWindowIcon` at run time. Nothing more to do.
- **Delete `icon2.svg`** in the same PR (`grep -rn icon2` first). It is a data
  file: this plan is the request to delete it.
- **Text:** title `pymolish-studio`. Splash: one added line
  ("pymolish-studio, an independent fork"); upstream lines unchanged. About
  box, in this order (wording from
  [step 10, Wording](10-readme-license.md#wording)):

  ```
  pymolish-studio <version>
  Based on PyMOL(TM) technology.

  An independent fork of Open-Source PyMOL. Not a Schrodinger product;
  not affiliated with or endorsed by Schrodinger, LLC.

  Copyright (C) Schrödinger, LLC.          <- existing lines, unchanged
  All rights reserved.
  New code: MIT, pymolish-studio contributors.
  <repository URL>

  PyMOL is a trademark of Schrodinger, LLC.
  ```

  Dropped: `https://pymol.org` as the product page and `sales@schrodinger.com`.
  No "(TM)" on our own name.
- **Rename `make_pymol_qicon`?** No. Plugins may import it.

## Steps

1. Add the two PNGs; point `make_pymol_qicon` at the icon.
2. Title, About box, splash line.
3. Spec `icon=`; Pillow in the build job.
4. Delete `icon2.svg`; README image is handled in [step 10](10-readme-license.md).
5. Check: one test in `testing/tests/system/` loads the PNG with `QImage`
   (works headless), asserts it is not null and square, that `cmd.splash`
   output still contains the Schrodinger copyright line, and that the About
   text contains the slogan and the trademark footnote.
6. Look at it: title bar, task bar, About box on Linux; the Windows artifact.

Commit: `ui(brand): own application icon, logo and window title`.
`ChangeLog`: `cmd.splash` prints one more line.

## Status

Done 2026-11-14 in 4.8.0 with the two files the user delivered: icon
1254 x 1254 (not resized), logo 2172 x 724, both with a transparent
background. Checked: `testing/tests/system/branding.py`; real GUI on WSLg
shows the title `pymolish-studio` and the icon, readable at 16 px;
PyInstaller's own converter turns the PNG into `.ico` and `.icns`.

Differences from the plan above:

- No line in `build.yml`: Pillow was already pinned in the `dev` extra, which
  the bundle jobs install.
- `README.md` points at `pymolish.png` so the deleted `icon2.svg` leaves no
  broken image. The rest of the README is still [step 10](10-readme-license.md).
- The Help menu row reads "About pymolish-studio".

Open:

- **Rights line** for the two images (who made them, which license): not
  delivered yet. `docs/licenses.md` has the row and says so.
- **Windows artifact:** icon on `pymolish.exe` not looked at; first CI run.
- **Logo on dark pages:** `logo.png` has near-black text on a transparent
  background. The README uses the icon for that reason.
- **Edge pixels:** both files have a thin halo of faint red pixels left by
  the background removal. Not visible at icon sizes; clean in the source art.

## Skipped

- Splash image, `.desktop` file, installer icons: no installer exists yet.
- Renaming the `pymol` package, `PYMOL_DATA`, or the product: separate decision.
- Replacing the Lucide icons.

---

PyMOL is a trademark of Schrodinger, LLC.
