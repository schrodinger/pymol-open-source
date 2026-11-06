"""
Tests for pymol.testing.import_from_file (testing/testing.py), which the
test driver uses to load itself as "pymol.testing" and every unittest file.
"""
import sys
import warnings

import pytest

from pymol import testing

MODULE_NAME = 'pymol_test_import_from_file_dummy'


@pytest.fixture
def module_name():
    sys.modules.pop(MODULE_NAME, None)
    yield MODULE_NAME
    sys.modules.pop(MODULE_NAME, None)


def _import(filename, name=None):
    # Python 3.15 removed SourceFileLoader.load_module(); older versions
    # emit a DeprecationWarning when it is called.
    with warnings.catch_warnings():
        warnings.simplefilter('error', DeprecationWarning)
        return testing.import_from_file(str(filename), name)


def test_returns_registered_module(tmp_path, module_name):
    path = tmp_path / 'dummy.py'
    path.write_text('import sys\n'
                    'VALUE = 42\n'
                    'SELF = sys.modules[__name__]\n')

    module = _import(path, module_name)

    assert module.VALUE == 42
    assert module.__name__ == module_name
    assert module.__file__ == str(path)
    assert sys.modules[module_name] is module
    # registered in sys.modules before the body runs, like load_module()
    assert module.SELF is module


def test_default_name(tmp_path):
    path = tmp_path / 'dummy.py'
    path.write_text('VALUE = 7\n')

    module = _import(path)
    try:
        assert module.VALUE == 7
        assert not module.__name__.endswith('.py')
        assert sys.modules[module.__name__] is module
    finally:
        sys.modules.pop(module.__name__, None)


def test_failed_import_is_not_registered(tmp_path, module_name):
    path = tmp_path / 'broken.py'
    path.write_text('raise RuntimeError("boom")\n')

    with pytest.raises(RuntimeError, match='boom'):
        _import(path, module_name)

    assert module_name not in sys.modules


def test_invalid_extension(tmp_path):
    path = tmp_path / 'dummy.txt'
    path.write_text('')

    with pytest.raises(ValueError, match='invalid extension'):
        _import(path)
