import sys
from pytest import mark
from typing import List, Union, Any, Tuple, Optional
from pathlib import Path

from pymol import cmd
from pymol.commanding import ArgumentParsingError


def test_docstring():
    @cmd.new_command
    def func():
        """
            docstring
            a
        """
    assert func.__doc__ == "docstring\na"


def test_bool(capsys):
    @cmd.new_command
    def func(a: bool, b: bool):
        assert a is True
        assert not b
    cmd.do("func yes, 0")
    out, err = capsys.readouterr()
    assert out + err == ''


def test_generic(capsys):
    @cmd.new_command
    def func(
        nullable_point: Tuple[float, float, float],
        my_var: Union[int, float] = 10,
        my_foo: int | float = 10.0,
        null_ptr: Optional[bool] = None,
        extended_calculation: bool = True,
        old_style: Any = "Old behavior"
    ):
        assert nullable_point == (1., 2., 3.)
        assert extended_calculation
        assert isinstance(my_var, int)
        assert isinstance(my_foo, float)
        assert null_ptr is None
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

    @cmd.new_command
    def func(a: List[bool]):
        assert a.pop(0) == False
        assert a.pop(0) == True
    cmd.do("func 0 yes")
    out, err = capsys.readouterr()
    assert out + err == ''


def test_tuple(capsys):
    @cmd.new_command
    def func(a: Tuple[str, int]):
        assert a == ("fooo a", 42)
    cmd.do("func 'fooo a' 42")
    out, err = capsys.readouterr()
    assert out + err == ''


def test_default(capsys):
    @cmd.new_command
    def func(a: str="sele"):
        assert a == "sele"
        
    cmd.do('func')
    out, err = capsys.readouterr()
    assert out + err == ''


def test_enum(capsys):
    from enum import Enum
    class E(Enum):
        A = 1
        B = 2
    @cmd.new_command
    def func(e: E):
        assert e == E.A
        assert isinstance(e, E)
    cmd.do('func A')
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
        B = "b"
    @cmd.new_command
    def func(e: E):
        assert e == E.A
        assert isinstance(e, E)
    cmd.do('func a')
    out, err = capsys.readouterr()
    assert out + err == ''


def test_quiet(capsys):
    @cmd.new_command
    def func(quiet: bool=True):
        assert not quiet
    cmd.do('func')
    out, err = capsys.readouterr()
    assert out + err == ''


def test_argument_error():
    err = ArgumentParsingError('my_var', "Short error message.")
    assert str(err) == "Failed at parsing 'my_var'. Short error message."


def test_call_error(capsys):
    @cmd.new_command
    def func(
        my_var: Union[int, float] = 10,
        my_foo: int | float = 10.0,
        null_ptr: Optional[bool] = None,
        extended_calculation: bool = True,
        old_style: Any = "Old behavior"
    ):
        assert isinstance(my_var, int)
        assert isinstance(my_foo, float)
        assert null_ptr is None
        assert extended_calculation
        assert old_style == "Old behavior"
    
    try:
        cmd._pymol.invocation.options.exit_on_error = 0

        cmd.do("func my_foo=a")
        out, err = capsys.readouterr()
        assert "Failed at parsing 'my_foo'." in (out + err)

        cmd.do("func extended_calculation=a")
        out, err = capsys.readouterr()
        assert (
            "Failed at parsing 'extended_calculation'."
            " Can't parse 'a' as bool."
            " Supported true values are yes, 1, true, on, y."
            " Supported false values are no, 0, false, off, n."
        ) in (out + err)

    finally:
        cmd._pymol.invocation.options.exit_on_error = 1
