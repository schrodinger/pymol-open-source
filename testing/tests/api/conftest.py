import pytest

from pymol import cmd


@pytest.fixture(autouse=True)
def setup():
    cmd.reinitialize()
    yield
