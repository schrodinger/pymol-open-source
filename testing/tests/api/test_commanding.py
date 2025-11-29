from pytest import mark
from pymol import cmd
import sys
from typing import List, Union, Any, Tuple
from pathlib import Path


def test_docstring():
    @cmd.new_command
    def func():
        """docstring"""
    assert func.__doc__ == "docstring"


def test_bool(capsys):
    @cmd.new_command
    def func(a: bool, b: bool):
        assert a
        assert not b
    cmd.do("func yes, 0")
    out, err = capsys.readouterr()
    assert out == '' and err == ''


def test_generic(capsys):
    @cmd.new_command
    def func(
        nullable_point: Tuple[float, float, float],
        my_var: Union[int, float] = 10,
        my_foo: Union[int, float] = 10.0,
        extended_calculation: bool = True,
        old_style: Any = "Old behavior"
    ):
        assert nullable_point == (1., 2., 3.)
        assert extended_calculation
        assert isinstance(my_var, int)
        assert isinstance(my_foo, float)
        assert old_style == "Old behavior"
    cmd.do("func nullable_point=1 2 3, my_foo=11.0")
    out, err = capsys.readouterr()
    assert out + err == ''

def test_path(capsys):
    @cmd.new_command
    def func(dirname: Path = Path('.')):
        assert dirname.exists()
    cmd.do('func ..')
    cmd.do('func')
    out, err = capsys.readouterr()
    assert out + err == ''



@mark.skip("This function does not works as expected")
def test_any(capsys):
    @cmd.new_command
    def func(old_style: Any):
        assert old_style is RuntimeError
    func(RuntimeError)    
    cmd.do("func RuntimeError")
    out, err = capsys.readouterr()
    assert 'AssertionError' not in out+err

def test_list(capsys):
    @cmd.new_command
    def func(a: List):
        assert a[1] == "2"

    cmd.do("func 1 2 3")
    out, err = capsys.readouterr()
    assert out + err == ''

    @cmd.new_command
    def func(a: List[int]):
        assert a[1] == 2

    cmd.do("func 1 2 3")
    out, err = capsys.readouterr()
    assert out + err == ''

def test_tuple(capsys):
    @cmd.new_command
    def func(a: Tuple[str, int]):
        assert a == ("fooo", 42)

    cmd.do("func fooo 42")
    out, err = capsys.readouterr()
    assert out + err == ''


def test_default(capsys):
    @cmd.new_command
    def func(a: str="sele"):
        assert a == "sele"
        
    cmd.do('func')
    out, err = capsys.readouterr()
    assert out + err == ''

@mark.skipif(
    sys.version_info < (3, 11),
    reason="Requires StrEnum of Python 3.11+"
)
def test_str_enum(capsys):
    from enum import StrEnum
    class E(StrEnum):
        A = "a"
    @cmd.new_command
    def func(e: E):
        assert e == E.A
        assert isinstance(e, E)
    cmd.do('func a')
    out, err = capsys.readouterr()
    assert out + err == ''

