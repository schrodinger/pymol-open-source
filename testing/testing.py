'''
Infrastructure for PyMOL testing

Usage:
    pymol testing.py --run all             Run all tests
    pymol testing.py --run some/file.py    Run tests from given files

PyMOL test cases should subclass pymol.testing.PyMOLTestCase and provide
either one "runTest" method or at least one "test*" method.
'''

from __future__ import print_function

import os
import sys
import pymol
import collections
import platform
import inspect

try:
    WindowsError
except NameError:
    WindowsError = None

try:
    basestring
except NameError:
    basestring = (str, bytes)

def compareListFunction(x, y):
    return collections.Counter(x) == collections.Counter(y)

def import_from_file(filename, name=None):
    import imp
    if name is None:
        try:
             name = os.path.relpath(filename).replace('.', '_')
        except ValueError:
             name = os.path.basename(filename).replace('.', '_')
    for suffix in imp.get_suffixes():
        if filename.endswith(suffix[0]):
            break
    else:
        raise ValueError('invalid extension: "%s"' % filename)
    return imp.load_module(name, open(filename), filename, suffix)

if __name__ != 'pymol.testing':
    # pymol foo.py    -> __name__ == 'pymol'
    # pymol -r foo.py -> __name__ == '__main__'

    _file = inspect.currentframe().f_code.co_filename
    pymol.testing = import_from_file(_file, 'pymol.testing')
    pymol.testing.cli()

else:
    import uuid
    import time
    import unittest
    import itertools
    import tempfile
    import argparse

    try:
        import Image
    except ImportError:
        from PIL import Image
        sys.modules['Image'] = Image

    from pymol import cmd
    from pymol.invocation import options

    PYMOL_VERSION = cmd.get_version()
    PYMOL_EDU = 'Edu' in PYMOL_VERSION[0]
    is_win64bit = "Windows" in platform.system() and sys.maxsize > 2**32

    usage = 'pymol [pymol options] %s [test options]' % (os.path.basename(__file__))
    parser = argparse.ArgumentParser("pymol", usage=usage)
    parser.add_argument('--xml', action='store_true')
    parser.add_argument('filenames', nargs='*', default=[])
    parser.add_argument('--out', default=sys.stdout)
    parser.add_argument('--offline', action='store_true')
    parser.add_argument('--no-mmlibs', action='store_true')
    parser.add_argument('--no-undo', action='store_true')
    parser.add_argument('--verbosity', type=int, default=2)

    have_dash_dash = __file__.startswith(sys.argv[0]) or '--run' in sys.argv
    cliargs = parser.parse_known_args(None if have_dash_dash else [])[0]

    run_all = False
    max_threads = int(cmd.get('max_threads'))

    cmd.set('use_shaders')
    use_shaders = cmd.get_setting_boolean('use_shaders')

    pymol_test_dir = os.path.abspath(os.path.dirname(__file__))

    deferred_unlink = []
    deferred_rmtree = []

    class requires_version(object):
        '''
        Decorator for restricting to PyMOL version
        '''

        def __init__(self, version):
            self.version = version

        def _tupleize(self, strversion):
            r = []
            for x in strversion.split('.'):
                try:
                    r.append(int(x))
                except ValueError:
                    break
            return tuple(r)

        def __call__(self, func):
            if isinstance(self.version, int):
                test = self.version <= PYMOL_VERSION[2]
            elif isinstance(self.version, float):
                test = self.version <= PYMOL_VERSION[1]
            else:
                test = self._tupleize(self.version) <= self._tupleize(PYMOL_VERSION[0])

            if not test:
                return unittest.skip('version %s' % (self.version))(func)

            return func

    class requires(object):
        '''
        Decorator for test methods which only should be executed
        under certain conditions.

        Example:

            >>> @requires('gui')
            >>> def testSomething(self):
            >>>     do_something()
        '''
        def __init__(self, *flags):
            self.flags = flags

        def __call__(self, func):

            flags = dict.fromkeys(self.flags, True)
            flags_known = []

            def hasflag(flag):
                flags_known.append(flag)
                return flags.pop(flag, False)

            if hasflag('shaders') and not use_shaders:
                return unittest.skip('shaders')(func)

            if hasflag('gui') and options.no_gui:
                return unittest.skip('no gui')(func)

            if hasflag('incentive') and not options.incentive_product:
                return unittest.skip('no incentive')(func)

            if hasflag('no_edu') and PYMOL_EDU:
                return unittest.skip('no edu')(func)

            if hasflag('network') and cliargs.offline:
                return unittest.skip('no network')(func)

            if hasflag('mmlibs') and cliargs.no_mmlibs:
                return unittest.skip('no mmlibs')(func)

            if hasflag('undo') and cliargs.no_undo:
                return unittest.skip('no undo')(func)

            if hasflag('no_run_all') and run_all:
                return unittest.skip('skip with all')(func)

            if hasflag('multicore') and max_threads <= 1:
                return unittest.skip('no multicore')(func)

            if hasflag('properties') and not options.incentive_product:
                return unittest.skip('no pymol.properties')(func)

            if hasflag('freemol') and (not options.incentive_product or PYMOL_EDU):
                return unittest.skip('no freemol')(func)

            if hasflag('no_win64bit') and is_win64bit:
                return unittest.skip('skip 64bit')(func)

            if flags:
                raise ValueError('unknown flags: ' + ', '.join(flags)
                        + '; choices: ' + ', '.join(sorted(flags_known)))

            return func

    class mktemp(object):
        '''
        Context manager which returns a temporary filename and
        deletes the file in the end, if it exists.
        '''
        def __init__(self, suffix=''):
            self.filename = tempfile.mktemp(suffix)
        def __enter__(self):
            return self.filename
        def __exit__(self, exc_type, exc_value, traceback):
            if os.path.exists(self.filename):
                try:
                    os.remove(self.filename)
                except WindowsError:
                    deferred_unlink.append(self.filename)

    class mkdtemp(object):
        '''
        Context manager for temporary directory
        '''
        def __init__(self):
            self.name = tempfile.mkdtemp()
        def __enter__(self):
            return self.name
        def __exit__(self, exc_type, exc_value, traceback):
            import shutil
            if os.path.exists(self.name):
                try:
                    shutil.rmtree(self.name)
                except WindowsError:
                    deferred_rmtree.append(self.name)

    class foreachList(list):
        pass

    class foreach(object):
        '''
        Decorator to call a method with arguments.

        If you have multiple decorators, this one must be the first (outer
        most) one because it does not return a function and thus cannot be
        further processed by other decorators.

        Examples:

            >>> @testing.foreach(1, 2, 3)
            >>> @someotherdecorator
            >>> def testSomething(self, a):
            >>>     print a

            Will print:
            ... 1
            ... 2
            ... 3

            >>> @testing.foreach((1,'A'), (2,'B'))
            >>> def testSomething(self, a, b):
            >>>     print a, b

            Will print:
            ... 1 A
            ... 2 B

            >>> @testing.foreach.zip((1,2), ('A','B'))
            >>> def testSomething(self, a, b):
            >>>     print a, b

            Will print:
            ... 1 A
            ... 2 B

            >>> @testing.foreach.product((1,2), ('A','B'))
            >>> def testSomething(self, a, b):
            >>>     print a, b

            Will print:
            ... 1 A
            ... 1 B
            ... 2 A
            ... 2 B

        '''
        def __init__(self, *args):
            self.args = args

        def __call__(self, func):
            r = foreachList()
            for args in self.args:
                if not isinstance(args, (tuple, list)):
                    args = (args,)
                def wrapper(self, a=args):
                    return func(self, *a)
                r.append([wrapper, args])  # need to pass the arguments in to set test name
            return r

        @classmethod
        def zip(cls, *args):
            args = zip(*args)
            return cls(*args)

        @classmethod
        def product(cls, *args):
            args = itertools.product(*args)
            return cls(*args)

    class PyMOLTestCaseMeta(type):
        '''
        Metaclass for PyMOLTestCase. Plays together with the foreach decorator.
        '''
        def __init__(self, *a, **k):
            if self.__module__ == 'pymol.testing':
                return

            for k, v in list(vars(self).items()):
                if isinstance(v, foreachList):
                    for c, fargs in enumerate(v, 1):
                        f, args = fargs
                        # set test name to function name plus arguments (delimited by '_')
                        setattr(self, '%s__%s' % (k, '_'.join(str(e) for e in args)), f)
                    delattr(self, k)

    class TimingCM(object):
        '''
        Timing context manager
        '''
        def __init__(self, test, msg=None, max=None):
            self.test = test
            self.msg = msg
            self.max = max
        def __enter__(self):
            self.start = time.time()
        def __exit__(self, exc_type, exc_value, traceback):
            if exc_type:
                return
            delta = time.time() - self.start
            if self.max and delta > self.max:
                msg = 'slow: %fs > %fs' % (delta, self.max)
                if self.msg:
                    msg = self.msg + ', ' + msg
                raise AssertionError(msg)
            self.test.timings.append((self.msg, delta))

    class PyMOLTestCase(PyMOLTestCaseMeta("Base", (unittest.TestCase,), {})):
        '''
        Common PyMOL unit tests should subclass this.

        Each tests starts with a clean (reinitialized) PyMOL session and
        from the directory where the file is located.
        '''

        if sys.version_info.major > 2:
            assertEquals = unittest.TestCase.assertEqual
            assertItemsEqual = unittest.TestCase.assertCountEqual

        moddirs = {}

        def setUp(self):
            self.oldcwd = os.getcwd()
            cmd.reinitialize()
            cmd.viewport(640, 480)

            if cliargs.no_undo:
                cmd.set('suspend_undo', updates=0)

            cwd = self.moddirs[type(self).__module__]
            os.chdir(cwd)

            cmd.feedback('push')
            cmd.feedback('disable', 'all', 'details actions')
            self.timings = []

        def tearDown(self):
            cmd.feedback('pop')
            os.chdir(self.oldcwd)

        def _getColorTuple(self, color):
            if isinstance(color, (tuple, list)):
                return tuple(color)
            return cmd.get_color_tuple(color)
                
        def assertColorEqual(self, color1, color2):
            self.assertEqual(self._getColorTuple(color1), self._getColorTuple(color2))

        def assertImageEqual(self, img1, img2=None, delta=0, count=0, msg='images not equal'):
            '''
            Test if two images are the same.

            img1, img2 can be either filenames, Image (PIL) objects
            or numpy arrays.

            delta > 0 is for inexact match (image data is 0..255 int)

            count is the number of allowed pixel mismatches.
            '''
            import numpy

            if isinstance(img1, basestring) and not \
                    os.path.exists(img1):
                print(' Generating reference img:', img1)
                self.png(img1)
                return

            data1 = self.get_imagearray(img1)
            data2 = self.get_imagearray(img2)

            self.assertEqual(data1.shape, data2.shape,
                    'image shapes not equal ')

            diff = abs(data1 - data2)

            noff = numpy.sum(diff > delta)
            if noff > count * data1.shape[-1]:
                filename = tempfile.mktemp('diff.png')

                diffimg = Image.fromarray((255 - diff.reshape(data1.shape)).astype(numpy.uint8))
                diffimg.save(filename)

                self.assertTrue(False, msg + ' (%d) %s' % (noff, filename))

        def _imageHasColor(self, color, img, delta=0):
            if isinstance(color, str):
                color = [int(v*255) for v in cmd.get_color_tuple(color)]
            else:
                color = list(color)
            dim = img.shape[-1]
            if dim == len(color) + 1:
                dim -= 1
                img = img[...,:dim]
            diff = abs(img.reshape((-1, dim)) - color)
            return (diff - delta <= 0).prod(1).sum()

        def save_imagearray(self, img, filename=None):
            if not filename:
                filename = tempfile.mktemp('.png')

            img = Image.fromarray(img)
            img.save(filename)

            return filename

        def _assertImageHasColor(self, test, color, img, delta, msg):
            import numpy

            img = self.get_imagearray(img)
            has_color = self._imageHasColor(color, img, delta)

            if bool(has_color) != test:
                filename = self.save_imagearray(img)
                self.assertTrue(False, msg + ', ' + filename)

        def assertImageHasColor(self, color, img=None, delta=0, msg=''):
            if not msg:
                msg = 'no such color: ' + str(color)
            self._assertImageHasColor(True, color, img, delta, msg)

        def assertImageHasNotColor(self, color, img=None, delta=0, msg=''):
            if not msg:
                msg = 'color found: ' + str(color)
            self._assertImageHasColor(False, color, img, delta, msg)

        def assertImageHasTransparency(self, img=None):
            img = self.get_imagearray(img)
            self.assertTrue((img[:,:,3] < 255).any())

        def assertImageHasNoTransparency(self, img=None):
            img = self.get_imagearray(img)
            if img.shape[-1] == 4:
                self.assertTrue((img[:,:,3] == 255).all())

        def assertArrayEqual(self, a1, a2, delta=0, msg='arrays not equal', _not=False):
            '''
            Test if two (multi-)dimensional numeric arrays are (almost) equal.
            '''
            import numpy

            a1 = numpy.asarray(a1)
            a2 = numpy.asarray(a2)

            self.assertEqual(a1.shape, a2.shape, msg + ' (shape)')
            self.assertEqual(not _not, numpy.allclose(a1, a2, 0, delta), msg)

        def assertArrayNotEqual(self, a1, a2, delta=0, msg='arrays equal'):
            return self.assertArrayEqual(a1, a2, delta, msg, True)

        def timing(self, *args, **kwargs):
            '''
            Timing context manager for feedback and maximum runtime assertion.
            Will show the runtime in seconds next to the OK message if tests
            are run with verbose=2.

            Optional arguments:
            msg = string: short label
            max = float: maximum allowed runtime in seconds

            Example:

                >>> with self.timing():
                >>>     so_something()

                >>> # maximum runtime assertion
                >>> with self.timing(max=3.0):
                >>>     so_something()
            '''
            return TimingCM(self, *args, **kwargs)

        def datafile(self, filename):
            '''
            Return path to filename, the current directory and the data
            directory are searched for filename.
            '''
            if os.path.exists(filename):
                return filename
            return os.path.join(pymol_test_dir, 'data', filename)

        def get_imagearray(self, img=None, **kwargs):
            '''
            Get bitmap data as a numpy array.
            
            img can be either a filename or a Image (PIL) object.
            '''
            if PYMOL_EDU and (options.no_gui or 'ray' in kwargs):
                self.skipTest("edu no-ray")

            import numpy
            
            if img is None:
                with mktemp('.png') as filename:
                    self.png(filename, **kwargs)
                    return self.get_imagearray(filename)

            if isinstance(img, numpy.ndarray):
                return img

            if isinstance(img, basestring):
                img = Image.open(img)

            if not isinstance(img, Image.Image):
                raise TypeError('img must be filename or Image instance')
        
            return numpy.array(img.getdata(),
                               numpy.uint8).reshape((img.size[1], img.size[0], -1))

        def png(self, filename, *args, **kwargs):
            '''
            Save image to filename, with antialias=0.
            '''
            cmd.set('antialias', 0)
            cmd.png(filename, *args, **kwargs)
            cmd.draw()

        def ambientOnly(self):
            cmd.set('ambient', 1)
            cmd.set('antialias', 0)
            cmd.set('light_count', 1)
            cmd.set('depth_cue', 0)

            # needed for open-source
            cmd.set('reflect', 0)
            cmd.set('direct', 0)

    class PyMOLTestResult(unittest.runner.TextTestResult):
        def addSuccess(self, test):
            if not (self.showAll and test.timings):
                return super(PyMOLTestResult, self).addSuccess(test)

            unittest.result.TestResult.addSuccess(self, test)
            msg = 'ok (%s)' % ', '.join(
                    ('%s: %.3fs' % (m, t) if m else '%.3fs' % t)
                    for (m, t) in test.timings)
            self.stream.writeln(msg)

            filename = os.getenv("PYMOLTESTTIMINGS",
                    os.path.join(pymol_test_dir, "timings.tab"))
            with open(filename, "a") as handle:
                for i, (m, t) in enumerate(test.timings):
                    version = cmd.get_version()
                    buildinfo = version[3:] or [0, "", 0]
                    print('\t'.join([
                        '%f' % time.time(),
                        '%012x' % uuid.getnode(),
                        '%f' % t,
                        type(test).__name__ + '.' + test._testMethodName,
                        str(m or i),
                        version[0],
                        buildinfo[1],
                        '%d' % buildinfo[2],
                        platform.platform(),
                        platform.node(),
                    ]), file=handle)

    def run_testfiles(filenames='all', verbosity=2, out=sys.stderr, **kwargs):
        '''
DESCRIPTION

    Run one or multiple unit test files as a test suite.

USAGE

    run_testfiles file1 file2 ... [, verbosity [, out ]]
        '''
        import glob

        if filenames in ('all', ['all']):
            global run_all
            run_all = True
            filenames = os.path.join(pymol_test_dir, 'tests', '*', '*.py')

        if isinstance(filenames, basestring):
            filenames = [filename
                    for pattern in filenames.split()
                    for filename in glob.glob(cmd.exp_path(pattern))]

        suite = unittest.TestSuite()

        for filename in filenames:
            if os.path.isdir(filename):
                filenames.extend(glob.glob(os.path.join(filename, '*.py')))
                continue

            mod = import_from_file(filename)

            # hacky: register working directory with test cases
            dirname = os.path.abspath(os.path.dirname(filename))
            PyMOLTestCase.moddirs[mod.__name__] = dirname

            suite.addTest(unittest.defaultTestLoader
                    .loadTestsFromModule(mod))

        if not 'xml' in kwargs:
            kwargs['xml'] = False
        if kwargs['xml']:
            import xmlrunner
            testresult = xmlrunner.XMLTestRunner(output=out, verbosity=int(verbosity)).run(suite)
        else:
            if isinstance(out, str):
                out = open(out, 'w')
            testresult = unittest.TextTestRunner(stream=out,
                                                 resultclass=PyMOLTestResult, verbosity=int(verbosity)).run(suite)

        while deferred_unlink:
            os.unlink(deferred_unlink.pop())

        while deferred_rmtree:
            import subprocess
            subprocess.call(['rd', '/s', '/q', deferred_rmtree.pop()], shell=True)

        return len(testresult.errors) + len(testresult.failures)

    def cli():
        '''
        Test suite client application.
        '''
        if not cliargs.filenames:
            # silently do nothing
            return

        nfail = run_testfiles(**vars(cliargs))
        cmd.quit(nfail)

    cmd.extend('run_testfiles', run_testfiles)

