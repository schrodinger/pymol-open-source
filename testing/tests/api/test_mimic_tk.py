"""
Tests for the Tk -> Qt module shims (modules/pmg_qt/mimic_tk.py), which let
legacy Tk plugins import tkinter.messagebox and tkinter.filedialog in the
Qt GUI.
"""
import importlib.util
import sys
import warnings

import pytest

pytest.importorskip('pymol.Qt')

from pmg_qt import mimic_tk

_MISSING = object()


def _restore(namespace, key, value):
    if value is _MISSING:
        namespace.pop(key, None)
    else:
        namespace[key] = value


@pytest.fixture
def uncached_tk_submodules():
    """
    Drop cached tkinter submodules so the next import goes through
    sys.meta_path, and restore the previous state afterwards.
    """
    tkinter = pytest.importorskip('tkinter')
    saved = {}
    for fullname in mimic_tk.mapping:
        attr = fullname.rpartition('.')[2]
        saved[fullname] = (sys.modules.pop(fullname, _MISSING),
                           vars(tkinter).pop(attr, _MISSING))
    yield
    for fullname, (module, value) in saved.items():
        _restore(sys.modules, fullname, module)
        _restore(vars(tkinter), fullname.rpartition('.')[2], value)


@pytest.mark.parametrize('fullname', sorted(mimic_tk.mapping))
def test_loader_creates_mapped_module(fullname):
    importer = mimic_tk.MimicTkImporter()
    spec = importer.find_spec(fullname, None)

    module = importlib.util.module_from_spec(spec)
    assert module is mimic_tk.mapping[fullname]
    spec.loader.exec_module(module)


def test_loader_resolves_module_names(monkeypatch):
    monkeypatch.setitem(mimic_tk.mapping, 'tkinter.mimic_tk_test', 'colorsys')
    importer = mimic_tk.MimicTkImporter()
    spec = importer.find_spec('tkinter.mimic_tk_test', None)

    import colorsys
    assert importlib.util.module_from_spec(spec) is colorsys


def test_find_spec_ignores_other_modules():
    importer = mimic_tk.MimicTkImporter()
    assert importer.find_spec('tkinter.simpledialog', None) is None


def test_import_tkinter_submodules(uncached_tk_submodules):
    # Python 3.15 no longer falls back to load_module(); older versions
    # emit an ImportWarning before doing so.
    with warnings.catch_warnings():
        warnings.simplefilter('error', ImportWarning)
        import tkinter.messagebox
        from tkinter import filedialog

    assert tkinter.messagebox is mimic_tk.qtMessageBox
    assert filedialog is mimic_tk.qtFileDialog
    assert sys.modules['tkinter.messagebox'] is mimic_tk.qtMessageBox
    assert sys.modules['tkinter.filedialog'] is mimic_tk.qtFileDialog
