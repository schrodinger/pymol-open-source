from pytest import fixture, skip

from pymol.commanding import declare_command, Selection
from pymol import cmd, testing, stored


@fixture(scope="function")
def reinitialize():
    cmd.reinitialize()
    yield


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
    cmd.do("func True, no")


def test_declare_command_default():
    @declare_command
    def func(a: Selection = "sele"):
        assert a == "a"

    func("a")
    cmd.do("func b")


def test_declare_command_docstring():
    @declare_command
    def func():
        """docstring"""

    assert func.__doc__ == "docstring"


def test_alias(reinitialize):
    stored.v = None
    cmd.alias("foo", "/stored.v = 123")
    cmd.do("_ foo", echo=0)
    assert stored.v == 123


# def testCls(reinitialize):
#     cmd.set("internal_prompt", 0)
#     cmd.set("text", 1)
#     cmd.cls()

#     # check on black screen
#     img = self.get_imagearray()
#     assert not img[..., :3].any()


def test_delete():
    breakpoint()
    cmd.pseudoatom("m1")
    cmd.pseudoatom("m2")
    cmd.pseudoatom("m3")
    cmd.delete("m1 m2")
    assert set(cmd.get_names()) == set(["m3"])


@skip("tested with other methods")
def test_do(self):
    # tested with other methods
    pass
