# Writing a plugin

A plugin is one Python file (or a package) in a startup directory. pymolish
imports it at launch as `pymol.plugins.startup.<name>` and calls its
`__init_plugin__(app=None)`.

| | |
| --- | --- |
| User directory | `~/.pymol/startup/`, on Windows `%APPDATA%\pymol\startup\`. Plugin > Plugin Manager can add search paths (restart needed). |
| Bundled plugins | `data/startup/` in this repository |
| Needs | pymolish 4.3.0 or later, Python 3.13+, PySide6 through `pymol.Qt` only |
| Without a GUI (`pymolish -c`) | Commands load. A plugin that calls `addmenuitemqt` is skipped with a one-line notice. |

Both examples below were run on 4.5.10, headless and in the GUI.

## A command

`~/.pymol/startup/centroid.py`. Adds `centroid` to the command line with
help (`help centroid`) and tab completion. No Qt, so it also works in scripts
run with `-c`.

```python
from pymol import cmd, CmdException


def centroid(selection="all", quiet=1):
    """
    DESCRIPTION

        Print and return the geometric center of a selection.

    USAGE

        centroid [ selection ]
    """
    xyz = cmd.get_coords(selection)  # numpy array, only the atoms asked for
    if xyz is None:
        raise CmdException("no atoms in selection")
    center = xyz.mean(axis=0).tolist()
    if not int(quiet):
        print(" centroid: %.3f %.3f %.3f" % tuple(center))
    return center


def __init_plugin__(app=None):
    cmd.extend("centroid", centroid)
    cmd.auto_arg[0]["centroid"] = cmd.auto_arg[0]["zoom"]  # tab completion
```

`cmd.extend` registers the command-line keyword; from Python call the function
itself. Raise `CmdException` for user errors: the command line shows the
message without a traceback.

## A panel

`~/.pymol/startup/atomcount.py`. Adds Plugin > Atom Count, which opens a dock
the user can move, float, close and reopen.

```python
from pymol import cmd, gui, plugins

dock = None


def show_panel():
    global dock
    if dock is None:
        # Qt is imported on use, so loading the plugin costs nothing
        from pymol.Qt import QtCore, QtWidgets

        window = gui.get_qtwindow()
        dock = QtWidgets.QDockWidget("Atom Count", window)
        dock.setObjectName("dock_atomcount")  # the saved layout finds it by name
        label = QtWidgets.QLabel()
        button = QtWidgets.QPushButton("&Count visible atoms")
        button.clicked.connect(
            lambda: label.setText("<b>%d</b> atoms" % cmd.count_atoms("visible"))
        )
        body = QtWidgets.QWidget()
        layout = QtWidgets.QVBoxLayout(body)
        layout.addWidget(button)
        layout.addWidget(label)
        layout.addStretch()
        dock.setWidget(body)
        window.addDockWidget(QtCore.Qt.DockWidgetArea.RightDockWidgetArea, dock)
    dock.show()
    dock.raise_()


def __init_plugin__(app=None):
    plugins.addmenuitemqt("Atom Count", show_panel)
```

- **One instance.** Keep the dock in a module variable and show it again; a
  second menu click must not build a second panel.
- **Parent and name.** The main window owns the dock (`gui.get_qtwindow()`).
  A stable `objectName` starting with `dock_` lets the saved layout put it
  back where the user left it.
- **Keyboard.** Give buttons a mnemonic (`&Count`), keep the default tab order.
- **Threads.** Callbacks run on the GUI thread. Call `cmd` directly for quick
  queries; anything slow belongs in a command the user can interrupt, not in a
  Python loop inside a button handler.

## Looks: let the application decide

The application sets the theme (dark by default) and the size (Display > UI
Scale, for high-resolution screens). A plugin that hard-codes either one breaks
on somebody's screen.

- No text colors: `<font color=blue>` is unreadable on the dark theme. Use
  `<b>` for emphasis and the palette for everything else.
- A colored background needs its own text color:
  `setStyleSheet("background: #fc6; color: black")`.
- No pixel font sizes and no fixed widget sizes; use layouts with their default
  margins and spacing.
- Import Qt from `pymol.Qt`, not from `PySide6` or `PyQt5` directly.

## Data

Ask for what you need: `cmd.get_coords(selection)` returns a numpy array of
just those atoms; `cmd.get_model` builds every atom and bond as Python objects
and is the slow choice for numeric work. Import optional libraries inside the
function that uses them, so a missing package costs one feature, not the plugin.
