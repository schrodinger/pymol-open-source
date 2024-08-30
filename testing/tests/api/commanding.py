from __future__ import print_function

import sys
import pytest

import pymol
import __main__
from pymol import cmd, testing, stored

from typing import List




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
        assert isinstance(b, (Path, str)) and "/tmp" == str(b)
    func(1, "/tmp")
    cmd.do('func 1, /tmp')


def test_declare_command_default(capsys):
    from pymol.commanding import Selection
    @cmd.declare_command
    def func(a: Selection = "sele"):
        assert a == "sele"
    func()
    cmd.do("func")
    out, err = capsys.readouterr()
    assert out == ''

def test_declare_command_docstring():
    @cmd.declare_command
    def func():
        """docstring"""
    assert func.__doc__ == "docstring"

    @cmd.declare_command
    def func():
        """
        docstring
        Test:
            --foo
        """
    assert func.__doc__ == "docstring\nTest:\n    --foo"


def test_declare_command_type_return(capsys):
    @cmd.declare_command
    def func() -> int:
        return 1

    assert func() == 1
    out, err = capsys.readouterr()
    assert out == ''

    @cmd.declare_command
    def func():
        return 1
    assert func() == 1

def test_declare_command_list_str(capsys):
    @cmd.declare_command
    def func(a: List[str]):
        print(a[-1])

    func(["a", "b", "c"])
    cmd.do('func a b c')
    out, err = capsys.readouterr()
    assert out == 'c\nc\n'

def test_declare_command_list_int(capsys):
    @cmd.declare_command
    def func(a: List[int]):
        print(a[-1] ** 2)
        return a[-1] ** 2

    assert func([1, 2, 3]) == 9
    cmd.do('func 1 2 3')
    out, err = capsys.readouterr()
    assert out == '9\n9\n'


def test_declare_command_list_float(capsys):
    @cmd.declare_command
    def func(a: List[float]):
        print(a[-1]**2)
        return a[-1]**2

    assert func([1.1, 2.0, 3.0]) == 9.0
    cmd.do('func 1 2 3')
    out, err = capsys.readouterr()
    assert out == '9.0\n9.0\n'


def test_declare_command_bool(capsys):
    @cmd.declare_command
    def func(a: bool, b: bool):
        assert a
        assert not b

    func(True, False)

    cmd.do("func yes, no")
    out, err = capsys.readouterr()
    assert out == '' and err == ''