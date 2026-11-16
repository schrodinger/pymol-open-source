# What is different

This page is for people who know Open-Source PyMOL(TM) and want to know what
pymolish-studio changed. New to both? Start with
[Getting started](getting-started.md).

The short version: the molecular engine, the commands and the file formats are
the same. The window is new, some slow operations are faster, and old
interfaces and unused features are removed.

Numbers in brackets are pymolish-studio versions; the
[ChangeLog](../ChangeLog) has the details.

## The window

| What | How it works now |
| --- | --- |
| Toolkit | One GUI, Qt 6 (4.1.0). |
| Theme | Dark by default; Display > Dark Theme switches it (4.5.2). |
| Size | The percent button in the toolbar sets the size of text and controls, 100 to 200 %. It starts at 150 % (4.7.0). |
| Layout | Toolbar with icons on top, viewport in the middle, log and command line at the bottom. The layout is remembered (4.5.0, 4.5.14). |
| Objects panel | Lists objects, selections and groups. Click a row to show or hide, double-click to zoom; each row has the Action, Show, Hide, Label and Color menus (4.5.11, 4.5.15). |
| States panel | First, previous, play, next, last, jump to a state, speed. Appears when there is more than one state or movie frame (4.6.0). |
| Menus | The right-click menu opens at the pointer as a normal menu, with icons (4.5.18). |
| Name | Title, About box and icon say pymolish-studio (4.8.0). |
| Object colors | The palette button in the toolbar picks the colors that new objects get. The choice is remembered (4.10.0). |

### Object color palettes

Each new object takes the next color of the palette for its carbons. Pick one
in the toolbar or type `palette vivid`; `palette soft, 1` also recolors what is
loaded. `palette default` is the PyMOL order.

| Palette | Made for | Colors |
| --- | --- | --- |
| `vivid` | dark background | ![vivid](img/palette-vivid.png) |
| `deep` | white background, print | ![deep](img/palette-deep.png) |
| `pastel` | dark background | ![pastel](img/palette-pastel.png) |
| `soft` | white background | ![soft](img/palette-soft.png) |
| `colorblind` | any background | ![colorblind](img/palette-colorblind.png) |

`colorblind` is Okabe-Ito followed by Paul Tol's bright and muted sets, with
the published values. Only the first seven colors are designed as one set.

## Speed

Each change was measured before and after on the same machine. Method and all
numbers: [baseline.md](baseline.md).

| Work | Before | After | Version |
| --- | --- | --- | --- |
| Load a 17 MB session (`.pse`) | 0.92 s | 0.14 s | 4.5.10 |
| Build the surface of 1aon (58870 atoms), with a small ray | 11.97 s | 9.66 s | 4.5.12 |
| The same surface build | 9.44 s | 6.24 s | 4.7.7 to 4.7.9 |

## Getting it

- **Ready-to-run bundles** for Linux, Windows and macOS (arm64) on the
  [Releases page](https://github.com/adrian-jasinski/pymolish-studio/releases).
  No Python needed. Unsigned preview builds (4.4.0, 4.5.1).
- **The program is called `pymolish`.** In Python it is still
  `import pymol`, so scripts keep working (4.5.5).
- **Python 3.11 or newer** when you build it yourself (4.1.0):
  [Install from source](install.md).

## Removed

Each row was unused, outdated or had a better replacement.

| Removed | Version | Use instead |
| --- | --- | --- |
| PyQt5, PyQt6 and PySide2 support | 4.1.0 | PySide6 |
| Tk GUI, Pmw and the Tk plugin manager | 4.3.0 | The Qt GUI. Plugins that build Tk windows no longer work. |
| GLUT GUI | 4.3.1 | The Qt GUI |
| `chempy.fast` | 4.4.1 | Nothing; it was unused |
| VTK-m isosurface backend | 4.4.2 | The internal marching cubes |
| `mpeg_encode` movie encoder | 4.4.3 | ffmpeg or convert |
| XML-RPC server (`-R`), built-in HTTP server, web examples | 4.4.4 | Nothing |
| OpenVR support | 4.4.5 | Nothing |
| COLLADA (`.dae`) export and import | 4.4.6 | glTF export |
| MMTF import and export | 4.4.7 | BinaryCIF |
| Startup image in the viewport | 4.6.1 | The notice is still printed to the log |

## Kept

- The `cmd` API and the plugin loader: [Writing a plugin](plugins.md).
- Reading and writing PDB, mmCIF, SDF/MOL and PSE sessions; PNG export.
- Ray tracing, movies, the POV-Ray renderer, the VMD molfile plugins.

What comes next is in the [roadmap](../PLAN.md).

---

PyMOL is a trademark of Schrodinger, LLC.
