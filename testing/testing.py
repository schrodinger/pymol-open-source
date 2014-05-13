'''
Infrastructure for PyMOL testing

Usage:
    pymol testing.py --run all             Run all tests
    pymol testing.py --run some/file.py    Run tests from given files

PyMOL test cases should subclass pymol.testing.PyMOLTestCase and provide
either one "runTest" method or at least one "test*" method.
'''

import os
import sys
import pymol
import collections
import platform

try:
    WindowsError
except NameError:
    WindowsError = None

def compareListFunction(x, y):
    return collections.Counter(x) == collections.Counter(y)

def import_from_file(filename, name=None):
    import imp
    if name is None:
        name = os.path.relpath(filename).replace('.', '_')
    for suffix in imp.get_suffixes():
        if filename.endswith(suffix[0]):
            break
    else:
        raise ValueError('invalid extension: "%s"' % filename)
    return imp.load_module(name, open(filename), filename, suffix)

if __name__ != 'pymol.testing':
    # pymol foo.py    -> __name__ == 'pymol'
    # pymol -r foo.py -> __name__ == '__main__'

    _cli = __name__ in ('pymol', '__main__')
    _file = pymol.__script__ if _cli else __file__
    pymol.testing = import_from_file(_file, 'pymol.testing')
    pymol.testing.cli()

else:
    import uuid
    import time
    import unittest
    import itertools
    import tempfile
    import argparse

    from pymol import cmd
    from pymol.invocation import options

    usage = 'pymol [pymol options] %s [test options]' % (os.path.basename(__file__))
    parser = argparse.ArgumentParser(usage=usage)
    parser.add_argument('--run', dest='filenames', nargs='*', default=[])
    parser.add_argument('--out', default=sys.stdout)
    parser.add_argument('--offline', action='store_true')
    parser.add_argument('--verbosity', type=int, default=2)
    cliargs = parser.parse_known_args()[0]
    run_all = False
    max_threads = int(cmd.get('max_threads'))

    cmd.set('use_shaders')
    use_shaders = cmd.get_setting_boolean('use_shaders')

    pymol_test_dir = os.path.abspath(os.path.dirname(__file__))

    deferred_unlink = []

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

            if hasflag('network') and cliargs.offline:
                return unittest.skip('no network')(func)

            if hasflag('no_run_all') and run_all:
                return unittest.skip('skip with all')(func)

            if hasflag('multicore') and max_threads <= 1:
                return unittest.skip('no multicore')(func)

            if hasflag('properties') and not options.incentive_product:
                return unittest.skip('no pymol.properties')(func)

            if hasflag('freemol') and not options.incentive_product:
                return unittest.skip('no freemol')(func)

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
                os.remove(self.filename)

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
                shutil.rmtree(self.name)

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

            for k, v in vars(self).items():
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

    class PyMOLTestCase(unittest.TestCase):
        '''
        Common PyMOL unit tests should subclass this.

        Each tests starts with a clean (reinitialized) PyMOL session and
        from the directory where the file is located.
        '''
        __metaclass__ = PyMOLTestCaseMeta

        moddirs = {}

        def setUp(self):
            self.oldcwd = os.getcwd()
            cmd.reinitialize()
            cmd.viewport(640, 480)

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
                print ' Generating reference img:', img1
                self.png(img1)
                return

            data1 = self.get_imagearray(img1)
            data2 = self.get_imagearray(img2)

            self.assertEqual(data1.shape, data2.shape,
                    'image shapes not equal ')

            noff = numpy.sum(abs(data1 - data2) > delta)
            self.assertLessEqual(noff, count * data1.shape[-1], msg + ' (%d)' % noff)

        def _imageHasColor(self, color, img=None, delta=0):
            if isinstance(color, str):
                color = [int(v*255) for v in cmd.get_color_tuple(color)]
            else:
                color = list(color)
            img = self.get_imagearray(img)
            dim = img.shape[-1]
            if dim == len(color) + 1:
                color.append(255)
            if isinstance(delta, list) and dim == len(delta) + 1:
                delta.append(0)
            diff = abs(img.reshape((-1, dim)) - color)
            return (diff - delta <= 0).prod(1).sum()

        def assertImageHasColor(self, color, img=None, delta=0):
            self.assertTrue(self._imageHasColor(color, img, delta),
                    'no such color: ' + str(color))

        def assertImageHasNotColor(self, color, img=None, delta=0):
            self.assertFalse(self._imageHasColor(color, img, delta),
                    'color found: ' + str(color))

        def assertImageHasTransparency(self, img=None):
            img = self.get_imagearray(img)
            self.assertTrue((img[:,:,3] < 255).any())

        def assertImageHasNoTransparency(self, img=None):
            img = self.get_imagearray(img)
            if img.shape[-1] == 4:
                self.assertTrue((img[:,:,3] == 255).all())

        def assertArrayEqual(self, a1, a2, delta=0, msg='arrays not equal'):
            '''
            Test if two (multi-)dimensional numeric arrays are (almost) equal.
            '''
            import numpy

            a1 = numpy.asarray(a1)
            a2 = numpy.asarray(a2)

            self.assertEqual(a1.shape, a2.shape, msg + ' (shape)')
            self.assertTrue(numpy.allclose(a1, a2, 0, delta), msg)

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
            import Image, numpy
            
            if img is None:
                filename = tempfile.mktemp('.png')
                try:
                    self.png(filename, **kwargs)
                    return self.get_imagearray(filename)
                finally:
                    try:
                        os.unlink(filename)
                    except WindowsError:
                        deferred_unlink.append(filename)

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
            cmd.unset('antialias')
            cmd.png(filename, *args, **kwargs)
            cmd.draw()

        def ambientOnly(self):
            cmd.set('ambient', 1)
            cmd.set('antialias', 0)
            cmd.set('light_count', 1)
            cmd.set('depth_cue', 0)

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
                    print >> handle, '\t'.join([
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
                    ])

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
            filenames = os.path.join(pymol_test_dir, '*', '*.py')

        if isinstance(filenames, basestring):
            filenames = [filename
                    for pattern in filenames.split()
                    for filename in glob.glob(cmd.exp_path(pattern))]

        if isinstance(out, str):
            out = open(out, 'w')

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

        testresult = unittest.TextTestRunner(stream=out,
                resultclass=PyMOLTestResult,
                verbosity=int(verbosity)).run(suite)

        while deferred_unlink:
            os.unlink(deferred_unlink.pop())

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

