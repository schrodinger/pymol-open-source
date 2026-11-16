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
| Shorter menus | The A menu of an object has 10 rows instead of 19; the less used ones are under "more". Right click on a row of the Objects panel gives Show as, Color, rename, duplicate and delete at once (4.13.0). |
| Name | Title, About box and icon say pymolish-studio (4.8.0). |
| Object colors | The palette button in the toolbar shows the palettes as color strips. Picking one recolors the loaded objects, and new objects and "by chain" follow it. The choice is remembered (4.10.0, 4.13.1). |
| Presets and looks | The Presets button lists the presets in order, with one "show solvent" / "hide solvent" row and a glass preset (cartoon inside a transparent surface). Below them are looks for the whole scene: figure (white, black outline), soft, flat, comic, sketch and x-ray. Outlines appear in ray-traced images (4.14.0). Every preset but "simple" draws the protein as cartoon; the ligand presets draw the ligand as balls on thin sticks; a "cartoon transparency" row makes the cartoon see-through until the next preset (4.14.3). Presets wire (the cartoon as a wire mesh) and illustration (shades of blue inside a faint outlined surface, a colored line drawing when ray traced); "shades" in the color menus colors a chain, or the chains, in shades of one color (4.14.4). |

### Object color palettes

Each new object takes the next color of the palette for its carbons, and
C > by chain gives each chain the next color. Picking a palette in the toolbar
also recolors what is loaded; typed, `palette vivid` leaves it alone and
`palette soft, 1` recolors. `palette default` is the PyMOL order.

| Palette | Made for | Colors |
| --- | --- | --- |
| `vivid` | dark background | ![vivid](img/palette-vivid.png) |
| `deep` | white background, print | ![deep](img/palette-deep.png) |
| `pastel` | dark background | ![pastel](img/palette-pastel.png) |
| `soft` | white background | ![soft](img/palette-soft.png) |
| `colorblind` | any background | ![colorblind](img/palette-colorblind.png) |
| `nord` | dark background | ![nord](img/palette-nord.png) |
| `dracula` | dark background | ![dracula](img/palette-dracula.png) |
| `catppuccin` | dark background | ![catppuccin](img/palette-catppuccin.png) |

`colorblind` is Okabe-Ito followed by Paul Tol's bright and muted sets, with
the published values. Only the first seven colors are designed as one set.
`nord`, `dracula` and `catppuccin` are the accent colors of the editor themes
of those names (Nord frost and aurora, Dracula, Catppuccin Mocha).

## New commands and menu entries

| What | How it works now |
| --- | --- |
| Commands from the pymolish plugin | 19 commands are built in: groups, folders, UniProt structures and ligands, predicted structures, SMILES, sequence search. [Commands](commands.md) lists them (4.11.0 to 4.11.4, 4.12.9). |
| Clean | The A menu of an object, a group, a selection and of "all" has a "clean" submenu: remove waters, remove ions, remove waters and ions, remove ligands. Ions are what the `inorganic` selector matches; additives that contain carbon (glycerol, acetate) count as ligands (4.12.0). |
| Optimize a ligand | `clean` works: it minimizes one small molecule (up to 1000 atoms) with the MMFF94 force field of rdkit. In the menu: A > clean > optimize ligand geometry; the Clean button of the Builder uses it too. In Open-Source PyMOL the command only raises an error. The atoms around the ligand are not seen, and bond orders are used as they are (4.12.1). |
| Fix a structure | `fix_structure` repairs a protein with PDBFixer and loads the result as a new object: missing atoms of residues that are present, hydrogens for a pH, nonstandard residues replaced; water, ligands and ions kept or dropped. In the menu: A > clean > fix structure... opens the options. Gaps in the chain are not built (4.12.2). |
| Relax a protein | `minimize` relaxes a protein with OpenMM (Amber14, implicit solvent) while holding the heavy atoms near their place, and loads the result as a new object. 260 residues take about 5 seconds. In the menu: A > clean > optimize protein geometry. It needs complete residues: run fix structure first. Ligands have no parameters and must be left out (4.12.3). |

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
- **Python 3.13 or newer** when you build it yourself (4.12.8; it was
  3.11 since 4.1.0):
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
