# 4 — Features and external libraries

[Roadmap](../../PLAN.md) · Draft; preserve existing extension contracts

**Keep the `cmd` API and existing plugin loader as the integration boundary.**
Make a small command or Qt tool easy to add without a new framework, dependency
container, plugin marketplace or speculative C++ ABI.

## Status 2026-11-08

4a done in 4.3.0: `pymol.plugins.startup` namespace, headless no-op app, Qt `PmwMenuBar`, `app.root` raises a clear error, Qt-backed dialog proxies, `testing/tests/api/plugins.py`. 4b examples not started.

## Existing capabilities and corrections

- [Plugin loader](../../modules/pymol/plugins/__init__.py): startup discovery,
  `plugin_load`, `__init_plugin__` initialization and menu registration already
  exist. Installing an arbitrary Python package does not automatically register
  it with this loader; document an explicit supported import/registration route.
- [Legacy support](../../modules/pymol/plugins/legacysupport.py) and
  [Qt initialization](../../modules/pmg_qt/pymol_qt_gui.py) remain connected.
  Removing them before preserving discovery and menus would break current tools.
- [GUI access](../../modules/pymol/gui.py) provides `get_qtwindow`; use the actual
  window and existing helpers before proposing a public `cmd.add_dock` API.
- [querying.py](../../modules/pymol/querying.py) exposes `get_model`, `get_coords`
  and `get_coordset`. `get_coords` returns selection coordinates as an array;
  `get_coordset(copy=0)` explicitly borrows internal memory and can become invalid
  after reallocation or deletion. Do not describe these as equivalent views.
- [setup.py](../../setup.py) supplies native sources/libraries to the small
  [CMake file](../../CMakeLists.txt). A universal `find_package`/`--with-*`
  convention is not the current architecture and should not be added speculatively.

## 4a — Contract and migration before Tk removal

1. Inventory built-in plugins, startup directories and supported user examples.
   Classify command-only, Qt UI and Tk-only integrations; state which continue.
2. Preserve `cmd.extend`, documented initialization behavior, menu registration,
   startup discovery and the Qt plugin manager while detaching them from Tk.
   Preserve existing startup content rather than deleting its container blindly.
3. Keep command-only plugins usable in headless mode. Qt-only tools should report
   the absence of a GUI clearly and avoid importing Qt at module import when the
   command path does not need it.
4. Keep optional library imports inside the command/tool that needs them. Show a
   useful missing-package message when invoked; ordinary startup should still work.
5. Add focused checks for command loading without a GUI, Qt menu registration,
   startup discovery and an optional dependency being absent. Test before and
   after removing the adapters; record deliberate compatibility changes.

The Tk-only retirement is a user-facing change to list in its removal PR. Do not
promise automatic conversion of arbitrary third-party dialogs.

## 4b — Small examples after the first workspace increments

Document one command-only example and one Qt panel example using existing APIs.
Include installation/discovery, expected versions, menu/shortcut registration,
parent ownership, closing/reopening the panel and where errors appear. Use
`docs/plugins.md` when implementing this documentation; it is not created by this
planning task.

Start with a local `QDockWidget` in the example. Consider a small helper in the
existing GUI/plugin module only after two real consumers need identical behavior.
Its contract must address stable dock names, lifetime, duplicate activation,
GUI-thread access, headless failure and layout restoration. Do not put a GUI-only
helper into `cmd` merely to make it easy to find.

For analysis, query only required data. Use coordinate arrays for numeric work
and `get_model` when full atom/bond metadata is needed. Workers should own copied
input and return detached results through existing thread utilities; keep PyMOL
mutation on established safe paths. Never retain a borrowed coordinate view
across deletion, state changes or concurrent mutation without a proven lifetime.

## External library policy

| Need | Smallest supported approach |
| --- | --- |
| Numeric array work | Existing numpy |
| Table data, simple export, configuration | Qt model/view and Python stdlib |
| Optional chemistry/sequence feature | Plugin-owned dependency, imported on use; an extra only for a real retained feature |
| Native functionality missing from the core | Measure the need, then use a system/header-only dependency through the existing build path |
| Qt tool panel | Existing widget, action and docking APIs |

Core runtime dependencies remain numpy and PySide6. Any additional dependency
needs the repository's explicit reason and approval when implementing it. Keep
dev pins exact; avoid adding RDKit/Biopython extras merely for future use.

Plugin environments and binary bundles need different expectations. A wheel in a
user-managed Python environment can support separately installed native packages.
A frozen app contains a fixed interpreter and native dependencies; do not promise
that arbitrary `pip install` or native wheels will work inside it. Define and test
a limited external-plugin contract for the first bundle in
[step 5](05-builds-binaries.md), while retaining the wheel route for advanced use.

## Additional ideas and acceptance checks

- **Capability-driven UI:** use existing capabilities such as
  `cmd.get_capabilities()` where applicable to explain unavailable optional
  actions instead of catching unrelated failures after activation.
- **Lazy loading:** measure startup contribution before loading every optional
  plugin; preserve explicit autoload preferences.
- **Useful bug reports:** expose already-available version, renderer and feature
  information in the existing diagnostics/help flow before inventing telemetry.
- **One documented example per real pattern:** keep examples executable in tests
  where they contain logic; avoid a broad SDK or hot-reload infrastructure.

Completion means the examples load in their stated environments, work without
unrelated optional packages, close cleanly, and preserve headless/core behavior.
Defer stable native plugin ABI, plugin dependency resolution and public dock APIs
until actual extension consumers establish the requirement.
