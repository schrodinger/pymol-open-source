from pytest import fixture

from pymol.commanding import declare_command, Selection
from pymol import cmd

@fixture(scope='function')
def reinitialize():
    cmd.reinitialize()


def test_declare_command_casting():
   from pathlib import Path

   @declare_command
   def func(a: int, b: Path):
       assert isinstance(a, int) and a == 1
       assert isinstance(b, Path) and b == Path("/tmp")
   func(1, "/tmp")
   cmd.do('func 1, "/tmp"')

def test_declare_command_bool():
   @declare_command
   def func(a: bool, b: bool):
       assert a
       assert not b
   func("True", "no")
   cmd.do('func True, no')

def test_declare_command_default():
   @declare_command
   def func(a: Selection="sele"):
      assert a == "a"
   func("a")
   cmd.do('func b')

def test_declare_command_docstring():
   @declare_command
   def func():
      """docstring"""
   assert func.__doc__ == "docstring"