# 3 — Modern workspace and metadata panel

[Roadmap](../../PLAN.md) · Slices 1 to 4 shipped, 5 partly; slice 6 (menus) was done as [step 14](14-menus.md) in 4.13.0

**Evolve the existing main window into a clear molecular workspace.** Deliver
movable, resizable panels and a useful metadata table before replacing every
internal control. Preserve command-line workflows and existing `.ui` dialogs.

## Status 2026-11-08

Slices 1 to 3 shipped in 4.5.0: `modules/pmg_qt/object_table.py` (model + dock, headless test in `testing/tests/api/object_table.py`), named docks (`dock_command`, `dock_objects`, `dock_builder`), `QSettings` layout save/restore with off-screen recovery, Display > Panels with Reset Layout. Viewport-to-table synchronization, object controls migration (slice 4), menu cleanup (6) is open. Slice 5 started in 4.5.2: `modules/pmg_qt/theme.py` (Fusion + dark palette, default on, Display > Dark Theme, persisted); spacing/labels not touched yet. Slice 5 continued in 4.5.7: Display > UI Scale and the dark-theme color audit, see [Display scale, typography and color](#display-scale-typography-and-color).

Slice 4 shipped in 4.5.11: `modules/pmg_qt/object_panel.py` replaces the 4.5.0 object table and the object list of the viewport's internal GUI, see [Objects panel](#objects-panel). Open: menu cleanup (6), the items listed in that section, spacing audit of the remaining dialogs.

## Status 2026-11-16

Scenes are a dock since 4.13.4 (Workspace layout, step 6). Checks: `testing/tests/api/scene_panel.py`; real GUI on WSLg at 150 %: thumbnails, rename through the editor, Ctrl+Down, Enter. Not exercised: a row dragged with the real mouse (the drop goes through the same `moveRow` as Ctrl+Up / Ctrl+Down). The scene buttons the core draws in the viewport are still there next to the dock (`scene_buttons`); they are part of the open row "Sequence viewer, scene buttons".

## Reference products and current constraints

Use commercial PyMOL for familiar molecular workflows and Maestro for its
relationship between a workspace and a project table. Maestro's published changes
include table configuration and selecting the corresponding table entry from a
picked workspace atom. These are useful interaction references, not a requirement
to reproduce its complete suite or proprietary assets. [PyMOL product](https://www.pymol.org/),
[Maestro table improvements](https://www.schrodinger.com/life-science/download/release-notes/release-2025-2/),
[workspace-to-table selection](https://www.schrodinger.com/life-science/download/release-notes/release-2025-4/)

The [current main window](../../modules/pmg_qt/pymol_qt_gui.py) already has docks,
menus and command entry. The internal object controls and sequence viewer also
involve [Executive](../../layer3/Executive.cpp) and
[Ortho](../../layer1/Ortho.cpp). They are not standalone Qt widgets that can be
reparented into docks. Keep the sequence display initially; replace it only as a
separately scoped, tested feature.

## Proposed workspace

| Region | Initial content | Behavior |
| --- | --- | --- |
| Center | Existing OpenGL viewport | Stable parent/context; largest useful area |
| Left, later increment | Objects, groups and selections | Visibility controls and existing actions through named menus |
| Bottom | New object metadata table | Sort/filter, row selection and object navigation |
| Bottom, tabbed or separate | Existing command entry and log | History/completion retained; easy keyboard focus |
| Right, when opened | Existing properties, builder or tool panels | Dock, float, resize, close and reopen |
| Within viewport initially | Existing sequence display | Preserve current selection behavior and settings |

Use `QDockWidget` and `QMainWindow` capabilities. “Unpinned” initially means
floating/detached; IDE-style edge auto-hide is a separate, deferred requirement.
Provide View → Panels, Reset Layout, and keyboard-accessible move/float actions.
Qt supports movable, floatable and closable docks directly.
[QDockWidget](https://doc.qt.io/qt-6/qdockwidget.html)

## Implementation slices

1. **Docks and persistence:** improve the existing window, name docks stably,
   save geometry/state in `QSettings`, and restore after widgets exist. Handle a
   missing monitor or invalid layout by returning tools to a visible location.
   Keep Reset Layout reachable with all tool panels closed. Store UI layout in
   settings, not by changing PSE output. [QMainWindow persistence](https://doc.qt.io/qt-6/qmainwindow.html)
2. **Object table MVP:** use `QTableView`, `QAbstractTableModel` and a
   `QSortFilterProxyModel`. Start read-only with one row per object, not per atom.
   Provide object name, type, enabled state, atom count where applicable, and
   state count. Add existing object properties only when the current API exposes
   them; show unavailable values explicitly instead of inventing metadata.
3. **Selection and navigation:** select the corresponding object from a row and
   provide Zoom/Orient actions. Distinguish a highlighted row, an atom selection,
   and object visibility; selecting a row must not silently hide other objects.
   Add viewport-to-table synchronization after the basic model is stable.
4. **Object controls:** migrate the most-used internal object actions into Qt
   only after parity for groups, names, visibility and command-driven updates.
   Keep legacy controls until their replacement passes those checks, then remove
   the replaced path deliberately. Avoid permanent duplicate object panels.
5. **Visual consistency:** standardize spacing, typography and action labels;
   follow the system palette with minimal QSS where needed. Reuse available icons
   and Qt standard icons first. Keep text labels/tooltips for ambiguous actions.
6. **Menus and discovery:** simplify actual common tasks while retaining access
   to advanced commands. Consider a small action search only after menu cleanup
   shows a discovery problem; reuse registered `QAction`s rather than treating
   every `cmd` keyword as an executable UI action with no arguments.

## Workspace layout

Request 2026-11-12, with two screenshots of commercial PyMOL 3 as the
reference: an action toolbar on top, the content panel on the right with a
toggle bar under it, the log and command line at the bottom. We take the
arrangement and the kind of control in each place; no icons, artwork or
product text are copied, and features we do not have (timeline, license
banner) are not imitated.

```
+---------------------------------------------------------------------+
| File  Edit  Build  Movie  Display  Setting  Scene  Mouse  Wizard ... |
| [Residues v] Undo Redo | Zoom v  Orient  Rock  Presets v            |
|                               (movie)   Builder  Scenes  Draw/Ray v  More v |
+------------------------------------------------+--------------------+
|                                                | o all      A S H L C|
|                 viewport                       | o 1abc     A S H L C|
|                                                |   (sele)   A S H L C|
|                                                | [wizard buttons]    |
|                                                | Mouse v Wizard v SEQ >_  State 1/1 |
+------------------------------------------------+--------------------+
| log                                                                 |
| PyMOL> [ command line                              ] [progress] [Abort] |
+---------------------------------------------------------------------+
```

| Step | What | Where | State |
| --- | --- | --- | --- |
| 1 | **Action toolbar** on top: selection mode (what a click selects), Undo, Redo, Zoom (menu: all, selection, center, reset), Orient, Rock (shows its state), Presets (same list as A > preset for all). Right side: Builder, Scenes, Draw/Ray (the existing render popup), More (Unpick, Deselect, Get View, Rebuild, Properties, Reset, movie commands). The movie buttons moved to the states panel in 4.6.0 ([step 6](06-states-panel.md)). | New `modules/pmg_qt/action_toolbar.py`, a `QToolBar` | Done 4.5.14 |
| 2 | **Command area at the bottom**: log and command line move from the top dock to the bottom; progress bar and Abort sit at the end of the command line. The grid of 18 quick buttons is deleted: every one of them is in the toolbar or its More menu. | `pymol_qt_gui.py` | Done 4.5.14 |
| 3 | **Content panel look**: no header row, visibility shown as a dot, A/S/H/L/C as rounded chips (accent colors on the "all" row, neutral on the others; since 4.5.15 every shown row, see Row state and colors), taller rows. Still the same checkable tree, so the keyboard works as before. | `object_panel.py` | Done 4.5.14 |
| 4 | **Toggle bar** under the content panel: mouse mode (the Mouse menu), Wizard (the Wizard menu), SEQ (sequence viewer), `>_` (show or hide the command area). The states panel sits above it ([step 6](06-states-panel.md)). | `object_panel.py` footer, filled by the window | Done 4.5.14 |
| 5 | Saved layouts are versioned again (2) so the command area starts at the bottom. `pymolish -x` (no external GUI) hides the toolbar and the command area as before. | `pymol_qt_gui.py` | Done 4.5.14 |
| 6 | **Scenes as a docked panel** with thumbnails, like the reference's left panel: one row per scene, the current one marked; Enter or double click recalls, F2 renames, dragging or Ctrl+Up / Ctrl+Down reorders; Add, Update, Rename, Delete. Hidden until the toolbar's Scenes button, Scene > Scenes or Display > Panels opens it; its place is saved with the layout. Polled from the window's 0.5 s timer while shown (marked `ponytail:`). The separate Scene Panel window (397 lines) is deleted. | New `modules/pmg_qt/scene_panel.py` (dock `dock_scenes`), `pymol_qt_gui.py` | Done 4.13.4 |
| 7 | Icons for the toolbar, from the Lucide set of step 11. | `action_toolbar.py`, `icons.py` | Done 4.7.0 ([step 7](07-toolbar-icons-scale.md)) |
| 11 | **Menu look.** Request 2026-11-12 with a reference picture: rounded menus, roomy rows, thin-line icons, hairline separators. One style for every menu in the window (`pymol.sty`); `theme.round_menus` makes the window corners see-through. Icons: 43 Lucide drawings (ISC license, `data/pmg_qt/icons/lucide/`), matched to menu labels in `icons.py` and tinted with the text color; delete, reinitialize and quit keep their red. Menu titles are muted captions. Tear-off rows are removed. | `icons.py`, `pymol.sty`, `theme.py`, `object_panel.py` | Done 4.5.18 |
| 8 | Timeline. Not planned: no timeline feature exists in this code base. | | Skipped |
| 10 | **Menus stay on the screen under WSLg.** Reported 2026-11-12: with the window maximized, the Color menu of the Objects panel opened its submenu to the right, off the screen. Cause: under Wayland the compositor places popups, and WSLg's does not flip or slide them. On WSL the GUI now starts through X11 (`QT_QPA_PLATFORM=xcb;wayland` unless the user set one), where Qt places menus itself. Checked maximized at 150 %: menu inside the screen, submenu to its left. | `theme.py` `prefer_x11_on_wsl` | Done 4.5.17 |
| 9 | The viewport no longer draws its own `PyMOL>` prompt and feedback line (`internal_feedback` 0 when the command area exists; `pymolish -x` keeps it). Text typed while the viewport has focus goes to the command line, which is shown and focused; Ctrl and Alt combinations and the special keys stay viewport shortcuts. `set internal_feedback, 1` in pymolrc restores the old prompt and key handling. | `pymol_qt_gui.py` | Done 4.5.16 |

Checks: `testing/tests/api/action_toolbar.py` (selection mode, rock state,
presets menu, zoom), `object_panel.py` tests extended for the
footer. Real GUI on WSLg at 100 % and 150 %.

## Objects panel

Detail for slice 4. Request 2026-11-12: the right panel (object names with the
action, show, hide, label and color buttons) could not be resized with the
mouse in any discoverable way, and its letters only grow in whole steps through
`display_scale_factor`.

**Why a rebuild.** The old panel is drawn by the C++ core with 8x13 bitmap
glyphs; its layout is hard-coded in multiples of 8 pixels across Executive,
ButMode, Control, Wizard and PopUp. A fractional size there means rewriting the
text path of all five. Its menus, though, are plain data from
`modules/pymol/menu.py`, so a Qt panel shows the same menus with native text.

| Part | Decision | State |
| --- | --- | --- |
| Object list | `QTreeWidget` in the `dock_objects` dock, right side. One row per object or selection, group members nested, "all" first. Checkbox toggles visibility, double click zooms (changed in 4.5.15, see Row state and colors). | Done 4.5.11 |
| A / S / H / L / C | Same menu functions per object type as `CExecutive::click`, table `MENUS`; submenus are filled when first opened; color markup in labels becomes a swatch. Right click or the Menu key opens the five menus for the row. | Done 4.5.11 |
| Text size | Qt text: follows the application font and Display > UI Scale at any step, including 125 %. | Done 4.5.11 |
| Resize | Drag the dock edge. The dock also floats, tabs and closes (Display > Panels). | Done 4.5.11 |
| Wizards | Title, buttons and popups of the active wizard (`get_panel`, `get_menu`) under the list; the dock reopens when a wizard starts. | Done 4.5.11 |
| Mouse mode, selection mode, state | One status line under the list. | Done 4.5.11 |
| Old internal GUI | Off by default in the GUI; Display > Internal GUI or `set internal_gui` in pymolrc brings it back. Its left edge can now be dragged, with a resize cursor. Not deleted: remove the C++ panel in its own PR after the Qt panel has been used for a release. | Kept on purpose |
| Updates | Polled from the window's 0.5 s timer: names, types, visibility. Counts are looked up for tooltips and copied rows only. Marked `ponytail:`; upgrade path is a change notification from the executive. | Done 4.5.11 |
| Reorder and regroup by dragging rows | Use `order` and `group` meanwhile. | Open |
| Group open/closed shared with the core | Expansion is view state of the tree for now. | Open |
| Picking an atom selects its row | Viewport-to-panel synchronization. | Open |
| Menus the core opens (right click in the viewport or on a sequence residue, scene buttons, timeline) | Qt menus since 4.5.18: the core hands its menu list to `pymol.menu.popup`, which the window sets (`popup_at_cursor`). Without a window (`-x` keeps the hook too; `-c` has none) the GL popup is still used. | Done 4.5.18 |
| Sequence viewer, scene buttons | Still drawn by the core in whole scale steps. | Open |
| Per-object motion menu (movie column) | Only shown by the old panel when a movie exists. | Open |

### Row state and colors

Request 2026-11-12: a shown object should look like the "all" row, a hidden
one should recede, and one click on the row should switch between the two.

| Step | What | State |
| --- | --- | --- |
| 1 | **Chip colors:** A, S and H in three blues (dark to light), L violet, C stays green. One table, `ACCENTS`. | Done 4.5.15 |
| 2 | **Shown rows** use those colors on every row, not only on "all", next to the green dot. The neutral grey chips are deleted. | Done 4.5.15 |
| 3 | **Hidden rows** draw the same chips and the name darker; the row returns to full color while the mouse is over it. "all" never dims. | Done 4.5.15 |
| 4 | **Click toggles:** a click on the dot or the name shows or hides the object in the viewport (`enable` / `disable`) and the row follows at once. Space does the same for the current row. Ctrl or Shift click still only selects rows for copying. | Done 4.5.15 |
| 5 | **Double click** shows the object and zooms on it, since its first click already toggled the row. | Done 4.5.15 |

Checks: `testing/tests/api/object_panel.py` (rows and nesting, every name in
`MENUS` exists in `pymol.menu`, menu building, toggle, filter, delete, status
line, wizard buttons, click and double click on a row). Real GUI on WSLg at 150 %: panel, color menu,
measurement wizard, edge drag of the old panel (400 to 460).

## Display scale, typography and color

Detail for slice 5. Problem reported 2026-11-12 on a 2560x1600 laptop panel
under WSLg: every letter in menus, panels and plugins is too small, and the
Lighting Settings plugin has blue and red text that disappears on the dark theme.

**Cause, measured.** Qt 6 follows the scale the OS reports. WSLg (and bare X11)
report none: the panel arrives as 96 dpi, device pixel ratio 1, default font
"Sans Serif 9pt" = 12 px. Windows, macOS and Wayland/GNOME/KDE desktops report
a real scale, so they were already correct. The physical size WSLg reports is
made up (650 mm for a laptop panel), so the size cannot be detected.

| Item | Decision | State |
| --- | --- | --- |
| Text, widgets, icons, spacing | One knob: Display > UI Scale (100 to 300 %), saved in `QSettings` `ui/scale`, handed to Qt as `QT_SCALE_FACTOR` before the application object exists. It multiplies the OS scale. An exported `QT_SCALE_FACTOR` wins. Applies on next start. | Done 4.5.7; since 4.7.0 picked from a toolbar button, 100 to 200 %, default 150 % ([step 7](07-toolbar-icons-scale.md)) |
| Viewport menus, object panel, labels | Drawn by the C++ core with the integer `display_scale_factor`. Now rounded from the window scale, not truncated: 150 % gives 2, was 1. | Done 4.5.7 |
| Text colors | No hard-coded foreground color in Qt code, forms or bundled plugins: bold or the palette instead. Text on a hard-coded light background sets its own dark color. Links use the palette link color. | Done 4.5.7: Lighting Settings, APBS, PNG dialog, Fetch, Plugin Manager, script editor syntax colors |
| First start on an unscaled high-resolution screen | Offer 150 % once when the OS reports scale 1 and the screen is 2560 px wide or more. Offer, do not guess: a 27" monitor at that width is fine at 100 %. Either answer is saved, so it is asked once. | Done 4.5.8; removed in 4.7.0, 150 % is the default |
| 125 % in the viewport | `display_scale_factor` is an integer, so 125 % still draws viewport text at 1x. The object list no longer depends on it (Objects panel, 4.5.11); what remains is the viewport popup menus and the sequence viewer. | Open, only if reported |
| Live rescale without restart | Qt reads the factor once. Skip unless restarts prove painful. | Skipped |
| Separate font-size setting | A second knob next to UI Scale. Skip unless the scale alone is not enough. | Skipped |
| Fixed point sizes in rich text | `help_shortcut.ui` (13 pt, a macOS font) and `pluginmanager.ui` (9 pt, Ubuntu) carried a designer font in their HTML. Removed, so both use the application font. | Done 4.5.13 |
| Spacing and labels | Hand-built dialogs with zero margins (Lighting Settings fixed); remaining dialogs audited one at a time. | Open |

Rule for plugins (in [docs/plugins.md](../plugins.md)): do not set text colors
or pixel font sizes; use the palette, `<b>`, and layouts with default spacing.

Checks: `testing/tests/api/theme.py` (saved scale reaches Qt, exported value
wins, syntax colors readable on dark). Real GUI checked on WSLg at 100 % and
150 %: window scale 1.5, `display_scale_factor` 2.

## Table behavior and performance

- Use the proxy model's row mapping when sorted/filtered. Do not use a displayed
  row number as object identity. Refresh identity mapping after rename, deletion,
  reinitialize or session replacement; discard stale references.
- Batch metadata queries and keep a snapshot for display. `data()` should not
  query the molecular engine or copy coordinates for every cell paint.
- Start with Refresh and refresh after known actions in this window. Command
  scripts/plugins can also mutate the scene: mark the initial limitation clearly
  and add an existing notification hook when a real need justifies it. Avoid
  polling every atom or rebuilding the table on each frame.
- Prevent selection feedback loops. If a rename or delete invalidates a row during
  synchronization, refresh or clear it rather than selecting an unrelated object.
- Keep Qt model updates on the GUI thread. A worker may prepare safe detached
  results and queue them back; model APIs are not thread-safe.
  [Qt model threading](https://doc.qt.io/qt-6/qabstractitemmodel.html#thread-safety)
- Add residue/atom rows later through bounded queries or lazy loading. A million
  Python row dictionaries and widgets are not an acceptable default.

An initial manual-refresh implementation should carry a `ponytail:` comment naming
its ceiling (external mutations need Refresh) and upgrade path (existing change
notifications). Do not introduce a general event bus just for this table.

## Additional ideas, ordered by usefulness

- **Copy selected rows:** tab-separated clipboard output for quick reporting;
  add CSV export via the stdlib only when a workflow needs a saved table.
- **Column chooser and remembered widths:** useful once real metadata creates
  more columns than fit, without a new preferences dialog.
- **Selection inspector:** summarize the active object/selection before adding
  editable atom properties or a new molecular data model.
- **Source provenance:** show file source and load time only when recorded by a
  defined load path. Existing sessions may lack them. Decide persistence
  separately; do not silently alter PSE or invent historical values.
- **Saved workspace presets:** defer until testing demonstrates a need beyond
  the last layout plus Reset Layout.

## Acceptance checks and scope

Use a real GUI test build to load multiple structures, float/resize/tab/close and
reopen panels, restart, restore the layout, and recover after monitor changes.
Repeat key tasks with keyboard only and at multiple display scales on supported
OSes. Keep the GL viewport and its context stable while moving tool panels.

Leave runnable tests for table sorting/filtering and identity after delete/rename,
plus layout recovery where practical. Verify load/save, PNG, selection, command
history/completion and plugin menus. Measure table refresh and rotation with the
panel visible/hidden using step 1's baseline.

No second complete main-window implementation, custom-painted docking system,
external icon pack, editable atom grid, or sequence-view rewrite is required for
the first release. Add them only for a demonstrated unmet requirement.
