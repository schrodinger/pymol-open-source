import os
import sys

import pytest

import pymol
from pymol import cmd

has_multi_undo = 'multi_undo' in pymol.get_capabilities()


@pytest.fixture(autouse=True)
def setup():
    if not has_multi_undo:
        pytest.skip()
    cmd.reinitialize()
    yield


@pytest.fixture(autouse=True)
def undo_enable(setup):  # setup to run first
    cmd.undo_enable()
    yield
    cmd.undo_disable()
