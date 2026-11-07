# 0 — UI engine research

[Roadmap](../../PLAN.md) · Reviewed 2026-11-08 · Decision proposal

**Use Qt Widgets through PySide6 as the only supported binding.** Modernize the
existing application incrementally. This retains the native molecular renderer
and forms; it does not accelerate ray tracing or loading by itself.

## Status — implemented 2026-11-08 (uncommitted)

Done in the working tree, pending review and commit as
`remove(gui): standardize on PySide6`:

- `modules/pymol/Qt/__init__.py` imports PySide6 only. Kept for plugins:
  `PYQT_NAME`, `QT_API=pyside6`, the `QtWidgets.QAction/QActionGroup/QShortcut/
  QOpenGLWidget` names and `QtCore.pyqtSignal/pyqtSlot`. Dropped: PyQt5, PySide2,
  PyQt6 branches, `_Qt_pre`, `QtCoreModels`, `QT_VERSION`, `QT_VERSION_STR`,
  `QtOpenGL` import. Callers of `QtCoreModels` and `QT_VERSION_STR` migrated.
- `modules/pymol/Qt/utils.py`: `loadUi` uses `QUiLoader` only; `uic` subprocess
  path removed.
- `modules/pmg_qt/pymol_gl_widget.py`: `QOpenGLWidget` only; `QGLWidget`,
  accumulation-buffer and Qt-version branches removed. Requested GL profile
  unchanged (plan step 5).
- `modules/pmg_qt/TextEditor.py`: no PyQt5 import fallback.
- `pyproject.toml`: `PySide6>=6.12` is a runtime dependency; dev pin stays
  `6.12.0`; `requires-python = ">=3.11"` to match CI. Two `sys.version_info`
  guards in `commanding.py` removed.
- `INSTALL`, `ChangeLog` updated.

Verification on this machine (Linux x86_64, WSL2, Python 3.12.3, PySide6 6.12.0,
Qt 6.12.0, OpenGL 4.5 / GLSL 4.50 via WSLg):

- Shim import, alias presence and `diagnostics_qt()` output checked.
- Headless suite `pymol -ckqy testing/testing.py --run all`, built with VMD
  plugins on and `use-msgpackc=no` (mmtf-cpp headers not installed locally):
  997 tests, 277 skipped, 4 errors, all MMTF export tests needing msgpack.
  No failure touches Qt, commanding, or diagnostics.
- GUI smoke (`-kq` script): load PDB, turn, cartoon/sticks/spheres, color,
  select, `png`, `save`/reload PSE, quit. Passed (340 atoms after PSE reload).
  Behavior identical with the original shim and GL widget swapped back in:
  framebuffer `png ray=0` from a startup script is blank in both, and both
  segfault on `quit` under WSLg/xcb. Neither is a regression of this change.
  `png ray=1` in the same scenario produces a real image (5445 distinct colors).

Not verified here: Windows, macOS, stereo, high-DPI picking, plugin dialogs.
Those remain exit criteria before tagging a release.

## Repository evidence

- [Qt shim](../../modules/pymol/Qt/__init__.py): tries PyQt5, PySide2, PyQt6,
  then PySide6, unless overridden; supplies aliases used by callers.
- [Main window](../../modules/pmg_qt/pymol_qt_gui.py): `QMainWindow`, nested and
  tabbed docks, command entry, plugin menus, and reusable dialogs.
- [Viewport](../../modules/pmg_qt/pymol_gl_widget.py): Qt 6 uses
  `QOpenGLWidget`; a `QGLWidget` branch remains for older bindings.
- [Metadata](../../pyproject.toml): numpy is the only runtime dependency;
  PySide6 6.12.0 is a dev pin. Python allows 3.9+, while
  [CI](../../.github/workflows/build.yml) targets 3.11–3.14 on three OS families.
- [GUI commands](../../modules/pymol/_gui.py) and
  [plugins](../../modules/pymol/plugins/__init__.py) are shared integration
  points. The Qt UI still depends on Tk-era adapters.

## Options and tradeoffs

These recommendations are engineering judgments for this repository, not claims
that the other frameworks cannot build scientific applications.

| Option | Useful capability | Cost here | Recommendation |
| --- | --- | --- | --- |
| **PySide6 / Qt Widgets** | Existing forms, menus, model/view tables, docking and GL widget | Binding cleanup, interaction design, packaging validation | **Adopt.** It is Qt's official Python binding. [Qt for Python](https://doc.qt.io/qtforpython-6/) |
| PyQt6 / Qt Widgets | The same underlying Qt 6 toolkit | Binding/plugin migration without a demonstrated product gain | Do not switch. GPLv3/commercial licensing differs from PySide6. [Riverbank](https://www.riverbankcomputing.com/software/pyqt/intro) |
| PyQt5 / PySide2 | Older plugin compatibility | Extra bindings and legacy branches | Retire after compatibility checks; avoid a blanket Qt 5 end-of-life claim. |
| Qt Quick / QML | Declarative animation and a GPU scene graph | New UI code and native-renderer integration | Reconsider only for a requirement Widgets fails. QML does not automatically port PyMOL's renderer to Metal/Vulkan. [Qt scene graph](https://doc.qt.io/qt-6/qtquick-visualcanvas-scenegraph.html) |
| Dear PyGui / Dear ImGui | GPU-drawn interactive tools | Replace forms, interaction behavior, menus and viewport integration | No demonstrated migration benefit. Do not claim docking/persistence is categorically absent. [Dear PyGui](https://dearpygui.readthedocs.io/en/latest/about/what-why.html) |
| Tauri / Electron / Qt WebEngine | Web UI ecosystem | New web/native boundary around the renderer and Python API | Defer. They do not share one bundle-size estimate: Tauri uses platform webviews. [Tauri](https://v2.tauri.app/reference/webview-versions/) |
| wx / GTK / Tk | Other desktop widget choices | Rebuild working Qt integration without a demonstrated need | No prototype warranted in this scope. |

PySide6 has LGPLv3/GPLv3/commercial distributions; packaged modules and components
matter. Keep their required notices in binaries. An LTS label does not mean all
maintenance packages are freely available. [Qt distribution documentation](https://doc.qt.io/qtforpython-6.10/commercial/index.html)

## Draft implementation plan

1. Establish a GUI smoke scenario: load a local structure, rotate, pick atoms,
   change representation, open a dialog/plugin menu, export PNG, save/reopen PSE,
   and close. Record OS, architecture, Qt version, and GL renderer.
2. Promote PySide6 to runtime metadata with a verified floor; keep the exact dev
   pin. Start from the existing 6.12.0 pin, available on
   [PyPI](https://pypi.org/project/PySide6/6.12.0/). A lower floor such as 6.8
   requires minimum-version tests and compatible wheels across the support matrix.
3. In a separate removal PR, retain `pymol.Qt` as the import location and select
   PySide6 directly. Audit `PYQT_NAME`, `QT_API`, `QtCoreModels`, `pyqtSignal`,
   `pyqtSlot`, `QT_VERSION`, and QtWidgets aliases before removing them. Migrate
   live callers in the same change or keep the used alias.
4. Remove the `QGLWidget` branch once PySide6 is mandatory. Test stereo,
   high-DPI picking, multisampling, framebuffer export and context cleanup.
   Qt added stereo-buffer APIs in 6.5; the old no-stereo comment is not a current
   toolkit limitation. [QOpenGLWidget](https://doc.qt.io/qt-6/qopenglwidget.html)
5. Leave the GL profile unchanged. Fixed-function calls remain in
   [Scene](../../layer1/Scene.cpp) and [Ortho](../../layer1/Ortho.cpp). Changing
   the requested format alone does not remove those dependencies. [Apple's
   Core-profile migration guide](https://developer.apple.com/library/archive/documentation/GraphicsImaging/Conceptual/OpenGL-MacProgGuide/UpdatinganApplicationtoSupportOpenGL3/UpdatinganApplicationtoSupportOpenGL3.html)
6. Follow [step 2](02-cleanup-dependencies.md) for separate Tk and GLUT removals.

## Additional ideas worth testing

- **Recoverable layouts:** save dock state and geometry with Reset Layout.
  Give docks stable `objectName` values for restoration.
  [QMainWindow](https://doc.qt.io/qt-6/qmainwindow.html)
- **Floating tools, stable viewport:** detach ordinary docks first and keep the
  GL viewport central. GL reparenting can recreate its context; avoid adding that
  lifecycle problem merely to provide unpinned tools.
- **Usability before restyling:** repeat load/select/color/export using the
  keyboard at two display scales. Fix focus and clipping before adding assets.

## Exit criteria and deferrals

Binding cleanup needs the Python suite, plugin initialization, `.ui` loading and
GUI smoke checks on supported platforms. A headless import alone is insufficient.
The research recommendation is complete; this runtime verification remains work.

Skip alternate-engine prototypes, custom docking libraries and renderer rewrites
until a concrete Widgets limitation demonstrates the need.
