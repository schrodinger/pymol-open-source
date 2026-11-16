# 15 — Presets and looks

[Roadmap](../../PLAN.md) · Request 2026-11-16 · **15a done 2026-11-16 in
4.14.0.** 15b and 15c are open and wait for the user.

**The preset menu is ordered, has one solvent row instead of twin rows, and a
glass preset. Below it the toolbar menu has five looks: light, background and
outline for the whole scene.**

## What existed

- **Presets:** `modules/pymol/preset.py`, functions that take a selection and
  set representations, colors and per-object settings. `menu.presets`
  (`modules/pymol/menu.py`) listed 14 rows in no order; the toolbar's Presets
  button and A > preset on every row show the same list.
- **Twin rows:** simple and simple (no solvent), pretty and pretty (with
  solvent), publication and publication (with solvent). Each pair differs
  only in the water.
- **Light:** only Plugin > Lighting Settings, five buttons (Default, Metal,
  Plastic, Rubber, X-Ray) that set light and nothing else.
- **`c_peptide_marketing2.pse`,** the user's example for a white background.
  It differs from the defaults in: white background, the plugin's Rubber
  light (`ambient` 0.05, `direct` 0.2, `reflect` 0.5, `shininess` 10),
  `power` 1.67, `depth_cue` off, `ray_trace_mode` 1 with `ray_trace_color`
  black, `antialias` 3, `ray_opaque_background` 1, and a white surface at
  `transparency` 0.6.

## What was tried

Every candidate was ray traced and captured from the real viewport on 1oky
(cartoon, surface, ligand).

| Tried | Result | Verdict |
| --- | --- | --- |
| The settings of the session file | Matte colors and black outlines on white | Look **figure** |
| Plugin: X-Ray | Inverted light, bright edges; unlike anything else | Look **x-ray** |
| Plugin: Metal, Plastic, Rubber | On cartoon and on surface hard to tell from each other and from the default | Not in the menu; they stay in the plugin |
| ChimeraX `lighting soft` (no key or fill light, ambient 1.5 from 64 directions) and silhouettes | PyMOL has eight lights at most: `light_count` 8, `direct` 0, `ambient` 0.35, `reflect` 0.9, soft shadows (`ray_shadow_decay_factor` 0.1). Bright, even colors; with the glass preset close to a ChimeraX picture | Look **soft** |
| Illustration style (`ray_trace_mode` 3, `ambient` 1, no direct light, no shadows) | Unshaded color areas with outlines, after David Goodsell's pictures | Look **flat** |
| Glass: transparent white surface plus `ray_transparency_oblique` 1 | The surface gets denser toward its edges, in ray-traced images and in the viewport | Preset **glass** |
| Glass: more specular light, `ray_transparency_oblique_power` 3 and 8, `ray_transparency_contrast` 1.6 | No visible gain over the row above | Not set |
| Glass: `transparency_mode` 3 | Grey rims in the viewport when combined with the oblique setting | Not set |

Outlines come from the ray tracer (`ray_trace_mode`). The viewport has no
outline pass, so figure, soft and flat show their outline only in images made
with Draw/Ray > Ray. That is 15b.

## Decisions

- **A preset says what is shown, a look says how it is lit.** Presets work on
  a selection and are in both menus. Looks are global settings, so they are
  only in the toolbar menu, under a "Look:" caption.
- **One solvent row.** It reads "show solvent" or "hide solvent", acts at
  once, and the presets picked afterwards keep the choice. It starts hidden.
  While hidden, a preset runs on the selection without its water, so
  "ligands" and "technical" also draw no contacts to hidden water. "default"
  is the exception: it shows everything, as after loading.
- **No function is removed.** `preset.simple_no_solv`, `pretty_solv`,
  `pub_solv` and the others keep their behavior; only the menu rows changed.
  The menu calls `preset.from_menu(name, selection)`.
- **A look replaces the one before it.** `preset.look(name)` first restores
  the default of every setting that any look sets (`unset` on a global
  setting does that), then sets its own. No table of defaults.
- **The next preset undoes glass.** `_prepare` now also restores
  `surface_color` and the global `ray_transparency_oblique`, as it already
  did for `two_sided_lighting`.

## The menu after 15a

| Group | Rows |
| --- | --- |
| Reset | default · classified |
| Plain | simple · technical · ball and stick · b factor putty |
| Pictures | pretty · publication · glass |
| Binding | ligands · ligand sites ▸ · protein interface |
| Option | show solvent / hide solvent |
| Look (toolbar only) | default · figure · soft · flat · x-ray |

13 preset rows instead of 14: three twin rows left, glass and the solvent row
came.

## Slices

**15a, menu, glass and looks. Done in 4.14.0.**
- `preset.py`: `glass`, `from_menu`, `show_solvent`, `LOOKS`, `look`.
- `menu.py`: `presets` reordered, `preset_ligand_sites` through `from_menu`,
  new builder `looks`.
- `action_toolbar.py`: the Presets button appends the looks.
- Checks: `testing/tests/api/preset.py` (`testGlass`, `testMenuSolvent`,
  `testLook`), `menu.py` (`testPresets`), `action_toolbar.py`. Real window on
  WSLg: every look, glass, the menu. Ray-traced: every look, glass in three
  looks.
- Commit: `ui(presets): ordered menu, one solvent row, glass preset, looks`.

**15b, outline in the viewport. Open, the user decides.**
- A screen-space edge pass over the depth buffer after the scene is drawn,
  in `data/shaders/` and the Scene draw path of `layer1/`. New shader code,
  to be measured with `testing/bench_gui.py` before and after.

**15c, ChimeraX object colors. Open, the user decides.**
- One more entry in `palettes.PALETTES` (tan, sky blue, plum, light green,
  salmon, ...), the pattern of [step 11](11-color-palettes.md). One line and
  the docstring.

## Skipped

- Rows for Metal, Plastic and Rubber: see the table.
- Remembering the solvent choice after a restart. Add it with `QSettings`
  when asked.
- `look` as a typed command. `preset.look("flat")` works on the command line.
- A mark on the look that is active.
- The Lighting plugin reading `LOOKS`, and its sliders following a look.
- A new color for the sticks of "ball and stick". They are white, so on a
  white background they show only where a look draws shading or an outline.
- A ChimeraX representation preset. "classified" already shows cartoon,
  ligand sticks and ion spheres.

## Sources

- PyMOLWiki, [Ray](https://wiki.pymol.org/index.php/Ray): the four values of
  `ray_trace_mode`.
- ChimeraX, [lighting](https://rbvi.ucsf.edu/chimerax/docs/user/commands/lighting.html):
  what its presets soft, full and flat set.
- Oxford Protein Informatics Group,
  [Making Pretty Pictures in PyMOL v2](https://www.blopig.com/blog/?p=12018):
  white background, `ray_trace_mode` 1, `antialias` 3.
