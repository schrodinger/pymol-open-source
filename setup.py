#!/usr/bin/env python
#
# This script only applies if you are performing a Python Distutils-based
# installation of PyMOL.
#
# It may assume that all of PyMOL's external dependencies are
# pre-installed into the system.

from distutils.core import setup, Extension
from distutils.util import change_root
from distutils.errors import *
from distutils.command.install import install
from distutils.command.build_py import build_py
from glob import glob
import shutil
import sys, os, re

# non-empty DEBUG variable turns off optimization and adds -g flag
DEBUG = bool(os.getenv('DEBUG', ''))

WIN = sys.platform.startswith('win')
MAC = sys.platform.startswith('darwin')

# handle extra arguments
class options:
    osx_frameworks = True
    jobs = int(os.getenv('JOBS', 0))
    no_libxml = False
    no_glut = True
    use_msgpackc = 'guess'
    help_distutils = False
    testing = False

try:
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--no-glut', action="store_true")
    parser.add_argument('--glut', dest='no_glut', action="store_false",
            help="link with GLUT (legacy GUI)")
    parser.add_argument('--no-osx-frameworks', dest='osx_frameworks',
            action="store_false")
    parser.add_argument('--osx-frameworks', action="store_true",
            help="on MacOS use OpenGL and GLUT frameworks instead of shared "
            "libraries from XQuartz. Note that the GLUT framework has no "
            "mouse wheel support, so this option is generally not desired.")
    parser.add_argument('--jobs', '-j', type=int, help="for parallel builds "
            "(defaults to number of processors)")
    parser.add_argument('--no-libxml', action="store_true",
            help="skip libxml2 dependency, disables COLLADA export")
    parser.add_argument('--use-msgpackc', choices=('c++11', 'c', 'guess', 'no'),
            help="c++11: use msgpack-c header-only library; c: link against "
            "shared library; no: disable fast MMTF load support")
    parser.add_argument('--help-distutils', action="store_true",
            help="show help for distutils options and exit")
    parser.add_argument('--testing', action="store_true",
            help="Build C-level tests")
    options, sys.argv[1:] = parser.parse_known_args(namespace=options)
except ImportError:
    print("argparse not available")

if options.help_distutils:
    sys.argv.append("--help")

if True:
    import monkeypatch_distutils
    monkeypatch_distutils.set_parallel_jobs(options.jobs)


def forms_uic(build_lib='modules'):
    '''
    Convert Qt UI files in "modules/pmg_qt/forms" to Python files in place
    '''

def get_prefix_path():
    '''
    Return a list of paths which will be searched for "include",
    "include/freetype2", "lib", "lib64" etc.
    '''
    try:
        return os.environ['PREFIX_PATH'].split(os.pathsep)
    except KeyError:
        pass

    if sys.platform.startswith("freebsd"):
        return ["/usr/local"]

    X11 = ['/usr/X11'] * (not options.osx_frameworks)

    if sys.platform == 'darwin':
        for prefix in ['/sw', '/opt/local', '/usr/local']:
            if sys.executable.startswith(prefix):
                return [prefix] + X11

    if is_conda_env():
        if sys.platform.startswith('win'):
            return [os.path.join(sys.prefix, 'Library')]

        return [sys.prefix] + X11

    return ['/usr'] + X11


def is_conda_env():
    return (
        'conda' in sys.prefix or
        'conda' in sys.version or
        'Continuum' in sys.version or
        sys.prefix == os.getenv('CONDA_PREFIX'))


def guess_msgpackc():
    for prefix in prefix_path:
        f = os.path.join(prefix, 'include', 'msgpack', 'version_master.h')

        try:
            m = re.search(r'MSGPACK_VERSION_MAJOR\s+(\d+)', open(f).read())
        except EnvironmentError:
            continue

        if m is not None:
            major = int(m.group(1))
            if major > 1:
                return 'c++11'

    return 'no'


class build_py_pymol(build_py):
    def run(self):
        build_py.run(self)
        forms_uic(self.build_lib)

class install_pymol(install):
    pymol_path = None
    bundled_pmw = False
    no_launcher = False

    user_options = install.user_options + [
        ('pymol-path=', None, 'PYMOL_PATH'),
        ('bundled-pmw', None, 'install bundled Pmw module'),
        ('no-launcher', None, 'skip installation of the pymol launcher'),
        ]

    def finalize_options(self):
        install.finalize_options(self)
        if self.pymol_path is None:
            self.pymol_path = os.path.join(self.install_libbase, 'pymol', 'pymol_path')
        elif self.root is not None:
            self.pymol_path = change_root(self.root, self.pymol_path)

    def run(self):
        install.run(self)
        self.install_pymol_path()

        if not self.no_launcher:
            self.make_launch_script()

        if self.bundled_pmw:
            raise Exception('--bundled-pmw has been removed, please install Pmw from '
                    'https://github.com/schrodinger/pmw-patched')

    def unchroot(self, name):
        if self.root is not None and name.startswith(self.root):
            return name[len(self.root):]
        return name

    def copy_tree_nosvn(self, src, dst):
        ignore = lambda src, names: set([
        ]).intersection(names)
        if os.path.exists(dst):
            shutil.rmtree(dst)
        print('copying %s -> %s' % (src, dst))
        shutil.copytree(src, dst, ignore=ignore)

    def copy(self, src, dst):
        copy = self.copy_tree_nosvn if os.path.isdir(src) else self.copy_file
        copy(src, dst)

    def install_pymol_path(self):
        self.mkpath(self.pymol_path)
        for name in [ 'LICENSE', 'data', 'test', 'scripts', 'examples', ]:
            self.copy(name, os.path.join(self.pymol_path, name))

    def make_launch_script(self):
        if sys.platform.startswith('win'):
           launch_script = 'pymol.bat'
        else:
           launch_script = 'pymol'

        self.mkpath(self.install_scripts)
        launch_script = os.path.join(self.install_scripts, launch_script)

        python_exe = os.path.abspath(sys.executable)
        pymol_file = self.unchroot(os.path.join(self.install_libbase, 'pymol', '__init__.py'))
        pymol_path = self.unchroot(self.pymol_path)

        with open(launch_script, 'w') as out:
            if sys.platform.startswith('win'):
                # paths relative to launcher, if possible
                try:
                    python_exe = '%~dp0\\' + os.path.relpath(python_exe, self.install_scripts)
                except ValueError:
                    pass
                try:
                    pymol_file = '%~dp0\\' + os.path.relpath(pymol_file, self.install_scripts)
                except ValueError:
                    pymol_file = os.path.abspath(pymol_file)

                # out.write('set PYMOL_PATH=' + pymol_path + os.linesep)
                out.write('"%s" "%s"' % (python_exe, pymol_file))
                out.write(' %*' + os.linesep)
            else:
                out.write('#!/bin/sh' + os.linesep)
                out.write('export PYMOL_PATH="%s"' % pymol_path + os.linesep)
                out.write('"%s" "%s" "$@"' % (python_exe, pymol_file) + os.linesep)

        os.chmod(launch_script, 0o755)

#============================================================================

# should be something like (build_base + "/generated"), but that's only
# known to build and install instances
generated_dir = os.path.join(os.environ.get("PYMOL_BLD", "build"), "generated")

import create_shadertext
create_shadertext.create_all(generated_dir)

# can be changed with environment variable PREFIX_PATH
prefix_path = get_prefix_path()

inc_dirs = [
]

pymol_src_dirs = [
    "ov/src",
    "layer0",
    "layer1",
    "layer2",
    "layer3",
    "layer4",
    "layer5",
    "modules/cealign/src",
    generated_dir,
]

def_macros = [
    ("_PYMOL_LIBPNG", None),
    ("_PYMOL_FREETYPE", None),
    ("_PYMOL_INLINE", None),
]

libs = ["png", "freetype"]
lib_dirs = []
ext_comp_args = [
    "-Wno-narrowing",
    # legacy stuff
    '-Wno-char-subscripts',
    # optimizations
    "-ffast-math",
    "-funroll-loops",
    "-O0" if DEBUG else "-O3",
    "-fcommon",
]
ext_link_args = []
ext_objects = []
data_files = []
ext_modules = []

if True:
    # VMD plugin support
    inc_dirs += [
        'contrib/uiuc/plugins/include',
    ]
    pymol_src_dirs += [
        'contrib/uiuc/plugins/molfile_plugin/src',
    ]
    def_macros += [
        ("_PYMOL_VMD_PLUGINS", None),
    ]

if not options.no_libxml:
    # COLLADA support
    def_macros += [
        ("_HAVE_LIBXML", None)
    ]
    libs += ["xml2"]

if options.use_msgpackc == 'guess':
    options.use_msgpackc = guess_msgpackc()

if options.use_msgpackc == 'no':
    def_macros += [("_PYMOL_NO_MSGPACKC", None)]
else:
    if options.use_msgpackc == 'c++11':
        def_macros += [("MMTF_MSGPACK_USE_CPP11", None)]
    else:
        libs += ['msgpackc']

    pymol_src_dirs += ["contrib/mmtf-c"]

if options.no_glut:
    def_macros += [
        ("_PYMOL_NO_MAIN", None),
    ]

if options.testing:
    pymol_src_dirs += ["layerCTest"]
    def_macros += [("_PYMOL_CTEST", None)]

inc_dirs += pymol_src_dirs

#============================================================================
if MAC:
        libs += ["GLEW"]

        if options.osx_frameworks:
            ext_link_args += [
                "-framework", "OpenGL",
            ] + (not options.no_glut) * [
                "-framework", "GLUT",
            ]
            def_macros += [
                ("_PYMOL_OSX", None),
            ]
        else:
            libs += [
                "GL",
            ] + (not options.no_glut) * [
                "glut",
            ]

if WIN:
        # clear
        libs = []
        ext_comp_args = []

        def_macros += [
            ("WIN32", None),
        ]

        libs += [
            "Ws2_32",   # htonl
        ]

        if True:

            libs += [
                "glew32",
                "freetype",
                "libpng",
            ] + (not options.no_glut) * [
                "freeglut",
            ] + (not options.no_libxml) * [
                "libxml2",
            ]

        if DEBUG:
            ext_comp_args += ['/Z7']
            ext_link_args += ['/DEBUG']

        libs += [
            "opengl32",
        ]

if not (MAC or WIN):
        libs += [
            "GL",
            "GLEW",
        ] + (not options.no_glut) * [
            "glut",
        ]

if True:
    try:
        import numpy
        inc_dirs += [
            numpy.get_include(),
        ]
        def_macros += [
            ("_PYMOL_NUMPY", None),
        ]
    except ImportError:
        print("numpy not available")

if True:
    for prefix in prefix_path:
        for dirs, suffixes in [
                [inc_dirs, [("include",), ("include", "freetype2"), ("include", "libxml2")]],
                [lib_dirs, [("lib64",), ("lib",)]],
                ]:
            dirs.extend(filter(os.path.isdir, [os.path.join(prefix, *s) for s in suffixes]))

if True:
    if sys.platform.startswith("freebsd"):
        libs += ["execinfo"]

    # optimization currently causes a clang segfault on OS X 10.9 when
    # compiling layer2/RepCylBond.cpp
    if sys.platform == 'darwin':
        ext_comp_args += ["-fno-strict-aliasing"]

def get_pymol_version():
    return re.findall(r'_PyMOL_VERSION "(.*)"', open('layer0/Version.h').read())[0]

def get_sources(subdirs, suffixes=('.c', '.cpp')):
    return sorted([f for d in subdirs for s in suffixes for f in glob(d + '/*' + s)])

def get_packages(base, parent='', r=None):
    from os.path import join, exists
    if r is None:
        r = []
    if parent:
        r.append(parent)
    for name in os.listdir(join(base, parent)):
        if '.' not in name and exists(join(base, parent, name, '__init__.py')):
            get_packages(base, join(parent, name), r)
    return r

package_dir = dict((x, os.path.join(base, x))
        for base in ['modules']
        for x in get_packages(base))

ext_modules += [
    Extension("pymol._cmd",
              get_sources(pymol_src_dirs),
              include_dirs = inc_dirs,
              libraries = libs,
              library_dirs = lib_dirs,
              define_macros = def_macros,
              extra_link_args = ext_link_args,
              extra_compile_args = ext_comp_args,
              extra_objects = ext_objects,
    ),

    Extension("chempy.champ._champ",
        get_sources(['contrib/champ']),
        include_dirs=["contrib/champ"],
    ),
]

distribution = setup ( # Distribution meta-data
    cmdclass  = {
        'build_py': build_py_pymol,
        'install': install_pymol,
    },
    name      = "pymol",
    version   = get_pymol_version(),
    author    = "Schrodinger",
    url       = "http://pymol.org",
    contact   = "pymol-users@lists.sourceforge.net",
    description = ("PyMOL is a Python-enhanced molecular graphics tool. "
        "It excels at 3D visualization of proteins, small molecules, density, "
        "surfaces, and trajectories. It also includes molecular editing, "
        "ray tracing, and movies. Open Source PyMOL is free to everyone!"),

    package_dir = package_dir,
    packages = list(package_dir),
    package_data = {'pmg_qt': ['forms/*.ui']},

    ext_modules = ext_modules,
    data_files  = data_files,
)
