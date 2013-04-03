PyMOL Testing Framework
=======================

This is a python unittest based testing framework for PyMOL.

Feel free to clone this repository and contribute.
https://help.github.com/articles/fork-a-repo

How to Run Tests
----------------

From the system command line:

    pymol -q /path/to/pymol-testing/runall.pml

From the pymol command line:

    PyMOL> run /path/to/pymol-testing/testing.py
    PyMOL> run_testfiles path/to/*.py

How to write Tests
------------------

A test is a method of a `pymol.testing.PyMOLTestCase` subclass and is named  
with a "test" prefix. Test files go into subdirectories like "api/" or 
"performance/". Each test method will start

*   with a clean PyMOL session (does reinitialize)
*   from the directory of its file (does os.chdir)

Decorators
----------

*   `@testing.requires(keyword)`: Takes one or more keywords and skips the
    test, if the requirement is not met. Valid keywords include: gui,
    incentive, network, `no_run_all`, multicore
*   `@testing.foreach(args1, args2, ...)`: Decorates a test method which
    takes arguments

Context Managers
----------------

*    `PyMOLTestCase.timing(max=sec)`: Gives feedback on running time and
     makes the test fail if it takes longer than `max` seconds.
*    `PyMOLTestCase.mktemp()`: Generates a temporary file name and deletes
     the file at context exit.

Assertion Methods
-----------------

In addition to the standard assertion methods of the unittest.TestCase
class, these assertion methods are currently available:

*   `assertImageEqual(img1, img2=None, delta=0)`: Takes two images (as
    filenames, PIL objects or numpy arrays). If only one image is given,
    grab the current PyMOL scene as second imagae
*   `assertImageHasColor(color, img=None, delta=0)`: Takes a color name, RGB
    or RGBA color tuple and an image (or grab the current scene) and check
    if it contains the given color
*   `assertImageHasTransparency`
*   `assertImageHasNoTransparency`
*   `assertColorEqual`
*   `assertArrayEqual`

Code Example
------------

    from pymol import cmd, testing, stored

    class TestSomething(testing.PyMOLTestCase):

        def testFoo(self):
            cmd.fragment('ala')
            cmd.show_as('sticks')
            cmd.orient()
            self.assertImageEqual('ref-foo.png')

        @testing.requires('gui')
        def testBar(self):
            cmd.some_opengl_feature()
