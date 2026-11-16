# 15 — Presets and looks

[Roadmap](../../PLAN.md) · Request 2026-11-16 · **15a done 2026-11-16 in
4.14.0, 15d in 4.14.3, 15f in 4.14.4.** 15b, 15c and 15e are open and wait for
the user.

**The preset menu is ordered, has one solvent row instead of twin rows, a
glass preset and a cartoon transparency row. Every preset but "simple" draws
the protein as cartoon. Below it the toolbar menu has seven looks: light,
background and outline for the whole scene.**

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

Second round, request of 2026-11-16 (15d and 15e), same structure:

| Tried | Result | Verdict |
| --- | --- | --- |
| The gist's `BallnStick` on the ligand: `stick_radius` 0.07, `sphere_scale` 0.18 (hydrogens 0.13), carbons `grey85` | Small balls on thin sticks; the ligand stands out from the protein's lines in every look | In **ligands** and **ligand sites** |
| The gist's black sticks | Invisible on the default black background | `grey30` instead: visible on black and on white |
| `cartoon_transparency` 0.25, 0.5, 0.75 | The ligand shows through the cartoon; viewport and ray-traced | Row **cartoon transparency** |
| `ray_trace_mode` 3 with default light | Color steps, but dark | Not used |
| `ray_trace_mode` 3, bright light (`ambient` 0.5, `direct` 0.6, `reflect` 0.3), no shadow, `antialias` 1 | Bright colors in four steps inside thick outlines; `antialias` 1 makes the lines thicker | Look **comic** |
| `ray_trace_mode` 2 | Black outlines on white, a line drawing | Look **sketch** |
| Gallery, Goodsell-like: `ray_trace_gain` 0, 0.03 | The outline looks the same as with the default 0.12 | Not set |
| `ambient_occlusion_mode` 1 on a surface | Pockets get dark, more depth on an opaque surface; black blotches on the glass surface | 15e |
| Fog on a white background (`depth_cue` 1, `fog_start` 0.45), after the gallery's Geis picture | The back of the molecule fades to white | 15e |
| `cartoon_cylindrical_helices` 1 | Helices as cylinders, a schematic drawing | 15e |
| Gallery, "Complex Stylized Protein": grey balls inside a rim-only surface (`ray_transparency_oblique_power` 20) | Works; reads only on small proteins | 15e |
| `ray_texture` 2 (in the gist) on a surface | No visible difference | Not set |

Third round, request of 2026-11-16 with six example pictures (15f):

| Tried | Result | Verdict |
| --- | --- | --- |
| `cartoon_debug` 3: the cartoon builder writes line strips instead of triangle strips | Helices, sheets and loops as a wire mesh, in the viewport and ray traced; close to the two wire pictures | Preset **wire** |
| `cartoon_fancy_helices` 0 and `cartoon_sampling` 14 with the wire | Little difference; a denser mesh only gets darker | Not set |
| `spectrum` with three colors: a light tint, the color, a dark shade | One color from light to dark along the chain | `util.shades`, rows **shades** |
| The look flat on shades of blue | `ray_trace_mode` 3 rounds every color to four steps: the shades became cyan, white and teal | flat now uses `ray_trace_mode` 1 |
| Faint surface (white, `transparency` 0.85) around the cartoon, look flat | The surface has no outline: the ray tracer leaves out everything more than 5 % transparent | See the next row |
| `ray_trace_trans_cutoff` 1 | The surface gets its silhouette line and small crease marks, as in the proteinbase pictures | In the look **flat**; preset **illustration** |
| The same cutoff in the look figure with the glass preset | Heavy black lines over the glass | Only in flat |
| A pale blue surface at `transparency` 0.7 | A little closer to the examples, but wrong after recoloring in another color | White at 0.85 |
| Dark end of the shades at 40 % of the color | Lost on a black background | 50 % |

Outlines come from the ray tracer (`ray_trace_mode`). The viewport has no
outline pass, so figure, soft and flat show their outline only in images made
with Draw/Ray > Ray. That is 15b. Comic and sketch depend on it most: in
the viewport comic is a bright light on white, and sketch is the default
light on white. The preset illustration is flat blue shapes there.

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

- **Ribbon only in "simple".** `technical`, `ligands` and the ligand sites
  presets draw the protein as cartoon (15d). "ligand sites > cartoon" is now
  that preset without its surface.
- **The ligand of the ligand presets is a ball-and-stick drawing** (15d).
  `_prepare` already restored `sphere_scale`, `stick_radius` and
  `stick_color`, so the next preset undoes it.
- **Cartoon transparency is a submenu, not a switch:** off, 25%, 50%, 75%.
  It sets `cartoon_transparency` on the objects of the menu; `_prepare`
  restores it, so it lasts until the next preset.

- **"illustration" also picks the look flat** (15f). It is the one preset
  that sets a look: without the outline the picture is not there. Another
  preset does not undo the look; the Look rows do.
- **flat keeps colors as they are** (15f): `ray_trace_mode` 1. The color
  steps are comic's.
- **Shades are a color action, not a preset:** `util.shades(color,
  selection, chains)`, in the color menus next to rainbow and next to "by
  chain". Eight colors in the menu; any color name works when typed.

## The menu after 15f

| Group | Rows |
| --- | --- |
| Reset | default · classified |
| Plain | simple · technical · ball and stick · b factor putty |
| Pictures | pretty · publication · glass · wire · illustration |
| Binding | ligands · ligand sites ▸ · protein interface |
| Option | cartoon transparency ▸ · show solvent / hide solvent |
| Look (toolbar only) | default · figure · soft · flat · comic · sketch · x-ray |

16 preset rows, two more than before step 15: three twin rows left; glass,
wire, illustration, the solvent row and cartoon transparency came.

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

**15d, cartoon, ligand style, transparency, comic and sketch. Done in 4.14.3.**
- `preset.py`: `_ball_and_stick_ligand`; `technical`, `ligands`,
  `ligand_sites` show cartoon; `_prepare` restores `cartoon_transparency`;
  looks `comic` and `sketch`.
- `menu.py`: builder `cartoon_transparency`, one row in `presets`.
- Checks: `testing/tests/api/preset.py` (`testTechnical`, `testLigands`,
  `testLook`), `menu.py` (`testPresets`). Real window on WSLg: ligands on
  black and on white, ligand sites, technical, 50 % cartoon, comic, sketch,
  the menu. Ray-traced: the same, and every row of the table above.
- Commit: `ui(presets): cartoon in every preset, ligand balls and sticks,
  cartoon transparency, comic and sketch looks`.

**15f, wire, illustration, shades. Done in 4.14.4.**
- `preset.py`: `wire`, `illustration`; `_prepare` restores `cartoon_debug`;
  the look flat.
- `util.py`: `shades`. `menu.py`: builder `shades`, rows in `spectrum`,
  `by_chain` and `presets`.
- Checks: `testing/tests/api/preset.py` (`testWire`, `testIllustration`),
  `util.py` (`test_shades`), `menu.py` (`testShades`, `testPresets`). Real
  window on WSLg: wire, the preset after wire, illustration, shades by
  chain. Ray-traced: illustration in two colors, wire, publication and
  glass in the new flat.
- Commit: `ui(presets): wire and illustration presets, shades of one color`.

**15b, outline in the viewport. Open, the user decides.**
- A screen-space edge pass over the depth buffer after the scene is drawn,
  in `data/shaders/` and the Scene draw path of `layer1/`. New shader code,
  to be measured with `testing/bench_gui.py` before and after.

**15c, ChimeraX object colors. Open, the user decides.**
- One more entry in `palettes.PALETTES` (tan, sky blue, plum, light green,
  salmon, ...), the pattern of [step 11](11-color-palettes.md). One line and
  the docstring.

**15e, more effects. Open, the user picks.** Each was tried (table above)
except the first, which needs code.

| Effect | What it takes |
| --- | --- |
| Color steps in the viewport, so comic looks like comic while working | One new setting and a few lines in `data/shaders/compute_color_for_light.fs` that round the light to steps. Together with 15b it makes the viewport match the ray-traced comic picture. |
| Fog toward a white background, for depth in the white looks | Two settings in a look (`depth_cue` 1, `fog_start` 0.45), or a "depth" row. Works in the viewport too. |
| Helices as cylinders, a schematic preset | One preset: `cartoon_cylindrical_helices` 1 on the publication preset; `_prepare` restores it. |
| Ambient occlusion in the soft look | `ambient_occlusion_mode` 1 in the look; glass must switch it off for its own objects. To be timed on a large surface first. |
| "Balls in glass" preset from the gallery | One preset; only for small proteins. |
| Glow and blur with distance, as in the wire example pictures | Not tried, needs code: a blur pass over the drawn picture, the same kind of pass as 15b. |

## Skipped

- Rows for Metal, Plastic and Rubber: see the table.
- Remembering the solvent choice after a restart. Add it with `QSettings`
  when asked.
- `look` as a typed command. `preset.look("flat")` works on the command line.
- A mark on the look that is active.
- The Lighting plugin reading `LOOKS`, and its sliders following a look.
- A new color for the sticks of "ball and stick". They are white, so on a
  white background they show only where a look draws shading or an outline.
- The rest of the gist: Bondi radii, slate nitrogens, `valence_size`, its
  dash style, the orbital and density commands.
- A ball-and-stick ligand in "pretty", "publication" and "glass". They keep
  sticks.
- A mark on the active cartoon transparency.
- A cross-hatched mesh on sheets. The wire draws them as a zigzag.
- A color picker for the shades, and "shades (elem C)" rows.
- A ChimeraX representation preset. "classified" already shows cartoon,
  ligand sticks and ion spheres.

## Sources

- PyMOLWiki, [Ray](https://wiki.pymol.org/index.php/Ray): the four values of
  `ray_trace_mode`.
- ChimeraX, [lighting](https://rbvi.ucsf.edu/chimerax/docs/user/commands/lighting.html):
  what its presets soft, full and flat set.
- bobbypaton, [pymol_style.py](https://gist.github.com/bobbypaton/1cdc4784f3fc8374467bae5eb410edef):
  the `BallnStick` values.
- PyMOLWiki, [Gallery](https://pymolwiki.org/Gallery): Goodsell-like,
  Stylized Ball and Stick, Complex Stylized Protein, Tribute to Irving Geis.
- Oxford Protein Informatics Group,
  [Making Pretty Pictures in PyMOL v2](https://www.blopig.com/blog/?p=12018):
  white background, `ray_trace_mode` 1, `antialias` 3.
