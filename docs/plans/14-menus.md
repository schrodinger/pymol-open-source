# 14 — Action and right-click menus

[Roadmap](../../PLAN.md) · Request 2026-11-14 · **Done 2026-11-16 in 4.13.0**:
14a in 4.10.1, 14b on 2026-11-14, 14c in 4.13.0. What 14c did, and where it
left the draft, is under "14c result".

**The red squares are a rendering artifact, fixed in one line (14a, do it
early). The real cleanup (14b, 14c) waits until steps 12 and 13 have added
their rows, so the menus are audited once.**

## What exists

- **One source for every popup:** `modules/pymol/menu.py`, 1,854 lines, 129
  builder functions returning `[kind, label, command]` rows.
- **Three places show them:**

  | Place | Code |
  | --- | --- |
  | A / S / H / L / C buttons on a row | `MENUS` table, `modules/pmg_qt/object_panel.py:35` |
  | Right click on a row | the same five as submenus, `object_panel.py:566` |
  | Right click in the viewport | `pick_menu`, `main_menu` through `popup_at_cursor`, `object_panel.py:177` |

- **Red squares, the cause.** 26 destructive rows carry the color code `\933`
  (`del_col`, `rem_col`, `menu.py:25`). `fill_menu` (`object_panel.py:172`)
  tints the row's icon with that color when the label is in `icons.LABELS`,
  and otherwise draws a 12 px square of the color. "delete object",
  "delete group", "delete everything", "delete selections" and
  "remove waters" are not in `LABELS`, so they get a red square. They are not
  error markers. The same code correctly draws swatches in the color menus.
- **Duplicates and dead rows found so far:**

  | What | Where |
  | --- | --- |
  | `sele_action2`: 10 of its 16 rows are also in `sele_action`; one caller | `menu.py:1195`, `:1718` |
  | `group_action`: 17 of its 21 rows are also in `mol_action` | `menu.py:1221`, `:1252` |
  | `duplicate`, marked `# broken...` | `menu.py:1185` |
  | `test1`, `test2` builders, no caller found in Python or C++ | `menu.py:1487`, `:1494` |
  | zoom, orient, center, origin at the top of every A menu and in the toolbar | removed from the A and right-click menus in 4.12.5 (decision 2026-11-14); the toolbar and a double click on a row stay |

## Decisions

- **Keep `menu.py` and its function names.** Plugins wrap them (the pymolish
  plugin wraps `group_action`), and the C++ internal GUI reads the same lists.
  This step edits rows, it does not redesign the format.
- **Menu content is a user-facing change:** 14b produces the list, the user
  approves it, 14c applies it.

## Slices

**14a, red squares. Done in 4.10.1.**
- `icons.menu_icon`: fall back to the label's first word, so every
  "delete ..." and "remove ..." row gets the trash icon, tinted red.
  One line. Squares remain only where a row has a color and no icon: the
  color menus.
- Remove the broken `duplicate` row.
- Test: for every row built with `del_col` or `rem_col`, `menu_icon` returns
  an icon.
- Commit: `fix(ui): destructive menu rows show a red trash icon, not a square`.

**14b, audit. Done 2026-11-14; the table for approval is below.**
- A throwaway script calls every builder with a dummy name and prints
  menu, path, label, command. Not committed.
- From that table: rows with the same command in one menu, rows whose command
  raises `IncentiveOnlyException`, rows unused in the three places above.
- Output: one keep / merge / drop table in this file, for the user's approval.
  Also the before numbers: rows per menu, lines in `menu.py`.

### 14b result, 2026-11-14 (4.12.3): for approval

A throwaway script called every builder with real objects (a peptide, a
selection, a group, a map, a mesh, a surface, a distance) and walked the
submenus. Before numbers: `menu.py` has 1,869 lines and 127 builders
(1,809 lines after 4.12.5; the row counts below are from before it).

| Menu | Top rows | Rows at all levels |
| --- | --- | --- |
| `mol_action` (A on an object) | 23 | 120 |
| `group_action` (A on a group) | 25 | 121 |
| `sele_action` (A on a selection) | 19 | 145 |
| `sele_action2` (right click on a selected atom, "actions") | 16 | 108 |
| `all_action` (A on "all") | 12 | 69 |
| `map_action`, `mesh_action`, `surface_action`, `slice_action` | 14, 9, 9, 12 | 40, 25, 25, 14 |
| `mol_show`, `mol_hide`, `mol_labels` | 21, 22, 17 | 43, 23, 25 |
| `mol_color` | 15 | 878 |

No action menu has two rows with the same command, with one exception
(row 3 below). The color menus repeat colors on purpose (one color sits in
"greens" and in "tints").

| # | Verdict | What | Why |
| --- | --- | --- | --- |
| 1 | Drop | Builders `test1`, `test2`, `mol_toggle`, `simple_show`, `simple_hide` | No caller in Python, C++ or `data/`. |
| 2 | Drop | find > pi interactions > all (in six menus, one source row) | `cmd.pi_interactions` raises `IncentiveOnlyException`. "pi-pi" and "pi-cation" next to it work and stay. |
| 3 | Fix | modify > exclude > polymer and > inorganic | Copy-paste error: both run "not organic", like the organic row. |
| 4 | Drop | label > stereochemistry | The `stereo` property is only filled by `assign_stereo`, which raises `IncentiveOnlyException`; the row labels nothing. |
| 5 | Merge | `sele_action2` into `sele_action` with a flag | Their 16 and 19 top rows share 11. The right-click variant leaves out drag and clean, which its parent menu already has, and shows the five "modify" rows flat. |
| 6 | Merge | `group_action` and `mol_action` into one builder | 17 of the top rows are the same; the differences are labels (rename object / rename group), drag matrix against drag, and the four group rows. |
| 7 | Regroup | A menu of an object: rename, duplicate, delete · clean ▸ · preset ▸ · find ▸ · align ▸ · hydrogens ▸ · group ▸ · more ▸ | 19 top rows become 10. "more" holds drag, reset matrix, generate, assign sec. struc., state, masking, sequence, movement, compute. "copy to object" moves under "duplicate". |
| 8 | Add | Right click on a row: show as ▸, color ▸, rename, duplicate, delete on top, then the five submenus | The essentials one click away, as planned. |
| 9 | Keep | Everything else. Done in 4.12.5, against the first draft of this row: zoom, orient, center and origin left every A and right-click menu | The toolbar zooms, orients and centers; a double click on a row zooms on it. |

Rows 1 to 4 remove dead things and fix a bug: they change nothing a user
relies on. Rows 5 to 8 move rows that people know by position. Decided
2026-11-16: the user left the choice open ("select best option and apply");
all eight rows were applied, because each one removes rows or duplicated
code and none removes a function.

**14c, apply. Done 2026-11-16 in 4.13.0.**

### 14c result

| | Before 14b | Before 14c | After |
| --- | --- | --- | --- |
| `menu.py` lines | 1,869 | 1,809 | 1,745 |
| Builders | 127 | 127 | 122 |
| A on an object, top rows | 23 | 19 | 10 |
| A on a group, top rows | 25 | 21 | 13 |
| A on a selection, top rows | 19 | 15 | 16 |

- One builder for objects and groups: `mol_action(cmd, name, group=0)`.
  `group_action` stays as a one-line name, because the C++ panel asks for it
  and the pymolish plugin wraps it. The rest of the old menu is in `mol_more`.
- One builder for selections: `sele_action(cmd, name, picked=0)`.
  `sele_action2` is deleted; the right click on a selected atom calls
  `sele_action(..., 1)`, which leaves out drag and clean.
- Right click on a row of the Objects panel: `context_menu` in
  `object_panel.py`.

Where the result differs from the draft below:

- **"duplicate object" is one click** (`cmd.copy` to `<name>_cpy`), like
  "duplicate group". "copy to object" is under "more", not under "duplicate".
- **The selection menu gained a row:** "duplicate selection" was only in the
  right-click variant. Order is now rename, duplicate, delete in every menu.
- **The right-click variant of the selection menu** shows "modify" as a
  submenu instead of five flat rows, and has align and hydrogens.
- **Right click on a row has Color once**, on top: Show as, Color · rename,
  duplicate, delete · Action, Show, Hide, Label. A row type without a part
  (a measurement has no "show as" and no labels) leaves it out.
- **A group keeps its three member rows on top** (sort, enable, disable).

Checks: `testing/tests/api/menu.py` `testConsolidated` (row counts, no two
rows with one command, the commands compile, duplicate works twice) and
`testExclude`; `testing/tests/api/object_panel.py` `test_context_menu`. Real
window on WSLg at 150 %: the right-click menu, A on an object, A on a group.

The draft, as written before the work:

- One builder for selection actions (`sele_action2` becomes a flag); one for
  object and group actions.
- Target A menu, about 10 top rows:
  rename, duplicate, delete · clean ▸ · preset ▸ ·
  find ▸ · align ▸ · hydrogens ▸ · group ▸ · more ▸ (drag, reset
  matrix, generate, state, masking, sequence, movement, compute).
- Right click on a row: essentials flat on top (show as ▸, color ▸,
  rename, duplicate, delete), then the five submenus as today.
- Essentials added: a working `duplicate`, the group rows from
  [step 12](12-native-commands.md), the clean submenu from step 13. Anything
  else comes from the approved table, not from this draft.
- Delete `test1`, `test2` and every builder left without a caller
  (`grep -rn` across Python, C++, `data/`, tests first).
- Test: no menu has two rows with the same command; every command string in
  the object, group and selection menus parses.
- Report rows and `menu.py` lines before and after. Target: fewer of both.
- Commit: `simplify(ui): merge duplicated action menus`.

## Skipped

- Rewriting `menu.py` as data tables or native `QMenu` code.
- Icons for every row.
- The viewport's internal C++ menu drawing.
- New shortcuts for menu rows.
