import sys

import __main__
from pytest import mark
from pymol import cmd, testing, stored
from typing import List
from typing import Optional, Any, Tuple, Union, List
from pathlib import Path


class TestCommanding(testing.PyMOLTestCase):

    def testAlias(self):
        stored.v = None
        cmd.alias('foo', '/stored.v = 123')
        cmd.do('_ foo', echo=0)
        self.assertEqual(stored.v, 123)

    @testing.requires('gui', 'no_run_all')
    def testCls(self):
        cmd.set('internal_prompt', 0)
        cmd.set('text', 1)
        cmd.cls()

        # check on black screen
        img = self.get_imagearray()
        self.assertFalse(img[...,:3].any())

    def testDelete(self):
        cmd.pseudoatom('m1')
        cmd.pseudoatom('m2')
        cmd.pseudoatom('m3')
        cmd.delete('m1 m2')
        self.assertEqual(cmd.get_names(), ['m3'])

    def testDeleteStates(self):
        cmd.pseudoatom('m1', state=10)
        self.assertEqual(cmd.count_states('m1'), 10)
        cmd.delete_states('m1', '4-8')
        self.assertEqual(cmd.count_states('m1'), 5)

    def testDo(self):
        # tested with other methods
        pass

    def testExtend(self):
        def check(v):
            stored.v = None
            cmd.do('foo', echo=0)
            self.assertEqual(stored.v, v)

        cmd.extend('foo', lambda: setattr(stored, 'v', 123))
        check(123)

        @cmd.extend
        def foo():
            stored.v = 456
        check(456)

    def testLog(self):
        with testing.mktemp('.pml') as logfile:
            cmd.log_open(logfile)
            cmd.do('_ color blue')
            cmd.log('hello world')
            cmd.log_close()
            lines = [_f for _f in map(str.strip, open(logfile)) if _f]
            self.assertEqual(lines, ['color blue', 'hello world'])

    @testing.requires_version('1.8.1.0')
    def testLog2(self):
        """
        Tests robustness of different logging methods:
        1) Python implementation (cmd.log) vs. C implementation (PLog, via cmd.do)
        2) pml vs. py syntax
        3) handling of quoted input (''' broken in 1.8.2)
        """
        self.ambientOnly()
        cmd.viewport(100, 100)
        cmd.fragment('gly')
        cmd.orient()
        cmd.show_as('spheres')

        for ext in ['.pml', '.py']:
            with testing.mktemp(ext) as logfile:
                cmd.log_open(logfile)
                cmd.do('_ color blue')
                cmd.do('/cmd.color("yellow", "elem O")')
                cmd.do('cmd.color("""green""",' " '''elem N''')")
                cmd.log('bg red\n')
                cmd.log('', 'cmd.color(\'magenta\', "hydro")\n')
                cmd.log_close()

                cmd.color('white')
                cmd.bg_color('white')

                if ext == '.pml':
                    cmd.do('@' + logfile)
                else:
                    cmd.do('run ' + logfile)

                img = self.get_imagearray()
                self.assertImageHasColor('blue', img)
                self.assertImageHasColor('yellow', img)
                self.assertImageHasColor('green', img)
                self.assertImageHasColor('red', img)
                self.assertImageHasColor('magenta', img)
                self.assertImageHasNotColor('white', img)

    def testLogClose(self):
        # see testLog
        pass

    def testLogOpen(self):
        # see testLog
        pass

    def testQuit(self):
        cmd.quit
        self.skipTest("cannot test quit")

    def testReinitialize(self):
        def check(v, names):
            self.assertEqual(v, cmd.get('pdb_conect_all'))
            self.assertEqual(cmd.get_names(), names)

        cmd.pseudoatom('m1')
        v = cmd.set('pdb_conect_all')
        cmd.reinitialize('store_defaults')

        cmd.reinitialize('original_settings')
        check('off', ['m1'])

        cmd.reinitialize('settings')
        check('on', ['m1'])

        cmd.reinitialize('purge_defaults')
        cmd.reinitialize('settings')
        check('off', ['m1'])

        cmd.reinitialize('everything')
        check('off', [])

    def testResume(self):
        with testing.mktemp('.pml') as logfile:
            with open(logfile, 'w') as f:
                print('bg yellow', file=f)
            cmd.resume(logfile)
            self.assertEqual('yellow', cmd.get('bg_rgb'))
            cmd.log('hello world')
            cmd.log_close()
            lines = [_f for _f in map(str.strip, open(logfile)) if _f]
            self.assertEqual(lines, ['bg yellow', 'hello world'])

    def testSplash(self):
        cmd.feedback('disable', 'all', 'output')
        cmd.splash()

    def testSync(self):
        cmd.sync()

    @testing.foreach(
        ('local', 'pymol', False),
        ('global', 'pymol', True),
        ('module', '', False),
        ('main', '__main__', True),
        ('private', '__main__', False),
    )
    def testRun(self, namespace, mod, rw):
        stored.tmp = False
        with testing.mktemp('.py') as filename:
            varname = '_tmp_' + namespace
            with open(filename, 'w') as handle:
                print('from pymol import stored', file=handle)
                if mod:
                    print('stored.tmp = (__name__ == "%s")' % (mod), file=handle)
                else:
                    print('stored.tmp = True', file=handle)
                print(varname + ' = True', file=handle)
            cmd.do('run %s, %s' % (filename, namespace), 0, 0)
            self.assertTrue(stored.tmp)
            if mod:
                self.assertEqual(rw, hasattr(sys.modules[mod], varname))


def test_declare_command_casting():
    from pathlib import Path
    @cmd.declare_command
    def func(a: int, b: Path):
        assert isinstance(a, int) and a == 1
        assert isinstance(b, Path) and "/tmp" == str(b)
    cmd.do('func 1, /tmp')


@mark.skip(reason="API not implemented yet")
def test_declare_command_optional(capsys):
    @cmd.declare_command
    def func(a: Optional[int] = None):
        assert a is None
    cmd.do("func")
    out, err = capsys.readouterr()
    assert out+err == ''

    @cmd.declare_command
    def func(a: Optional[int] = None):
        assert a is 10
    cmd.do("func 10")
    out, err = capsys.readouterr()
    assert out+err == ''

def test_declare_command_docstring():
    @cmd.declare_command
    def func():
        """docstring"""
    assert func.__doc__ == "docstring"


def test_declare_command_bool(capsys):
    @cmd.declare_command
    def func(a: bool, b: bool):
        assert a
        assert not b

    cmd.do("func yes, 0")
    out, err = capsys.readouterr()
    assert out == '' and err == ''


def test_declare_command_generic(capsys):
    @cmd.declare_command
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

def test_declare_command_path(capsys):
    @cmd.declare_command
    def func(dirname: Path = Path('.')):
        assert dirname.exists()
    cmd.do('func ..')
    cmd.do('func')
    out, err = capsys.readouterr()
    assert out + err == ''

def test_declare_command_any(capsys):
    @cmd.declare_command
    def func(old_style: Any):
        assert old_style != "RuntimeError"
    cmd.do("func RuntimeError")
    out, err = capsys.readouterr()
    assert 'AssertionError' in out+err

def test_declare_command_list(capsys):
    @cmd.declare_command
    def func(a: List):
        assert a[1] == "2"
    cmd.do("func 1 2 3")
    out, err = capsys.readouterr()
    assert out + err == ''

    @cmd.declare_command
    def func(a: List[int]):
        assert a[1] == 2
    cmd.do("func 1 2 3")
    out, err = capsys.readouterr()
    assert out + err == ''

def test_declare_command_tuple(capsys):
    @cmd.declare_command
    def func(a: Tuple[str, int]):
        assert a == ("fooo", 42)
    cmd.do("func fooo 42")
    out, err = capsys.readouterr()
    assert out + err == ''

def test_declare_command_arg_docs():
    @cmd.declare_command
    def func(
        # multiline
        # documentation works
        foo: int, # inline
        a: str,
        # bar are strings
        bar: Tuple[str, int], # continued...
        b: Any = 10, # The new old age
        # aaaa
        c: Any = 'a' # b
    ):
        "main description"
        pass

    assert func.__arg_docs['foo'] == "multiline documentation works inline"
    assert func.__arg_docs['a'] == ""
    assert func.__arg_docs['bar'] == "bar are strings continued..."
    assert func.__arg_docs['b'] == 'The new old age'
    assert func.__arg_docs['c'] == 'aaaa b'
    assert func.__annotations__['foo'] == int
    assert func.__annotations__['bar'] == Tuple[str, int]

def test_declare_command_default():
    @cmd.declare_command
    def func(a: str="sele"):
        assert a == "a"
    func("a")
    cmd.do('func a')