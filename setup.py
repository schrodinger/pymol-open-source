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
import platform

import multiprocessing.pool
import monkeypatch_distutils

# handle extra arguments
class options:
    osx_frameworks = False
    jobs = int(os.getenv('JOBS', 0))
    no_libxml = False
    use_msgpackc = 'c++11'
    help_distutils = False
    no_cxx11 = False

# OS X <= 10.8
if sys.platform == 'darwin' and tuple(
        map(int, platform.mac_ver()[0].split('.'))) < (10, 9):
    options.no_cxx11 = True

try:
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--osx-frameworks', action="store_true",
            help="on MacOS use OpenGL and GLUT frameworks instead of shared "
            "libraries from XQuartz. Note that the GLUT framework has no "
            "mouse wheel support, so this option is generally not desired.")
    parser.add_argument('--jobs', '-j', type=int, help="for parallel builds "
            "(defaults to number of processors)")
    parser.add_argument('--no-libxml', action="store_true",
            help="skip libxml2 dependency, disables COLLADA export")
    parser.add_argument('--use-msgpackc', choices=('c++11', 'c', 'no'),
            help="c++11: use msgpack-c header-only library; c: link against "
            "shared library; no: disable fast MMTF load support")
    parser.add_argument('--no-cxx11', action="store_true", help="Disable "
            "C++11 std library features. Will still require C++11 'auto' "
            "keyword support.")
    parser.add_argument('--help-distutils', action="store_true",
            help="show help for distutils options and exit")
    options, sys.argv[1:] = parser.parse_known_args(namespace=options)
except ImportError:
    print("argparse not available")

if options.help_distutils:
    sys.argv.append("--help")

if options.jobs != 1:
    monkeypatch_distutils.pmap = multiprocessing.pool.ThreadPool(options.jobs or None).map

def posix_find_lib(names, lib_dirs):
    # http://stackoverflow.com/questions/1376184/determine-if-c-library-is-installed-on-unix
    from subprocess import Popen, PIPE
    args = ["cc", "-shared", "-o", os.devnull] + ["-L" + d for d in lib_dirs]
    for name in names:
        p = Popen(args + ["-l" + name], stdout=PIPE, stderr=PIPE)
        p.communicate()
        if p.wait() == 0:
            return name
    raise IOError('could not find any of ' + str(names))

class build_py_pymol(build_py):
    def run(self):
        build_py.run(self)

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
            import tarfile
            pmwtgz = "modules/pmg_tk/pmw-py%d.tgz" % (sys.version_info[0])
            if not os.path.exists(pmwtgz):
                pmwtgz = "modules/pmg_tk/pmw.tgz"
            tar = tarfile.open(pmwtgz)
            tar.extractall(self.install_libbase)
            tar.close()

    def unchroot(self, name):
        if self.root is not None and name.startswith(self.root):
            return name[len(self.root):]
        return name

    def copy_tree_nosvn(self, src, dst):
        ignore = lambda src, names: set(['.svn']).intersection(names)
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
                out.write('set PYMOL_PATH=' + pymol_path + os.linesep)
                out.write('"%s" "%s"' % (python_exe, pymol_file))
                out.write(' %1 %2 %3 %4 %5 %6 %7 %8 %9' + os.linesep)
            else:
                out.write('#!/bin/sh' + os.linesep)
                if sys.platform.startswith('darwin'):
                    out.write('[ "$DISPLAY" == "" ] && export DISPLAY=":0.0"' + os.linesep)
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
prefix_path = ["/usr", "/usr/X11"]

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
    ("_PYMOL_INLINE", None),
]

libs = []
pyogl_libs = []
lib_dirs = []
ext_comp_args = [
    # legacy stuff
    '-Wno-write-strings',
    '-Wno-unused-function',
    '-Wno-char-subscripts',
]
ext_link_args = []
data_files = []
ext_modules = []

if options.no_cxx11:
    def_macros += [
        ('_PYMOL_NO_CXX11', None),
    ]
    if options.use_msgpackc == 'c++11':
        options.use_msgpackc = 'no'

if True:
    # VMD plugin support
    pymol_src_dirs += [
        'contrib/uiuc/plugins/include',
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

if options.use_msgpackc == 'no':
    def_macros += [("_PYMOL_NO_MSGPACKC", None)]
else:
    if options.use_msgpackc == 'c++11':
        def_macros += [("MMTF_MSGPACK_USE_CPP11", None)]
    else:
        libs += ['msgpackc']

    pymol_src_dirs += ["contrib/mmtf-c"]

inc_dirs = list(pymol_src_dirs)

#============================================================================
if sys.platform=='win32': 
    # NOTE: this branch not tested in years and may not work...
    inc_dirs += [
              "win32/include"]
    libs=["opengl32","glu32","glut32","libpng","zlib"]
    pyogl_libs = ["opengl32","glu32","glut32"]
    lib_dirs=["win32/lib"]
    def_macros += [
                ("WIN32",None),
                ("_PYMOL_LIBPNG",None),
                ]
    ext_link_args=['/NODEFAULTLIB:"LIBC"']
#============================================================================
elif sys.platform=='cygwin':
    # NOTE: this branch not tested in years and may not work...
    libs=["glut32","opengl32","glu32","png"]
    pyogl_libs = ["glut32","opengl32","glu32"]
    lib_dirs=["/usr/lib/w32api"]
    def_macros += [
                ("CYGWIN",None),
                ("_PYMOL_LIBPNG",None)]
#============================================================================
else: # unix style (linux, mac, ...)

    def_macros += [
            ("_PYMOL_FREETYPE",None),
            ("NO_MMLIBS",None),
            ]

    if sys.platform == 'darwin':
        for prefix in ['/sw', '/opt/local', '/usr/local']:
            if sys.executable.startswith(prefix):
                prefix_path.insert(0, prefix)
    elif sys.platform.startswith("freebsd"):
        prefix_path = ["/usr/local"]

    if 'anaconda' in sys.executable.lower():
        prefix_path.insert(0, sys.executable.split('/bin/')[0])

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

    libs += ["png", "freetype"]

    try:
        prefix_path = os.environ['PREFIX_PATH'].split(os.pathsep)
    except KeyError:
        pass

    for prefix in prefix_path:
        for dirs, suffixes in [
                [inc_dirs, [("include",), ("include", "freetype2"), ("include", "libxml2")]],
                [lib_dirs, [("lib64",), ("lib",)]],
                ]:
            dirs.extend(filter(os.path.isdir, [os.path.join(prefix, *s) for s in suffixes]))

    if sys.platform == 'darwin' and options.osx_frameworks:
        ext_link_args += [
            "-framework", "OpenGL",
            "-framework", "GLUT",
        ]
    else:
        glut = posix_find_lib(['glut', 'freeglut'], lib_dirs)
        pyogl_libs += ["GL", "GLU", glut]

    libs += ["GLEW"]
    libs += pyogl_libs

    ext_comp_args += ["-ffast-math", "-funroll-loops", "-fcommon"]

    # optimization currently causes a clang segfault on OS X 10.9 when
    # compiling layer2/RepCylBond.cpp
    if sys.platform != 'darwin':
        ext_comp_args += ["-O3"]

def get_pymol_version():
    return re.findall(r'_PyMOL_VERSION "(.*)"', open('layer0/Version.h').read())[0]

def get_sources(subdirs, suffixes=('.c', '.cpp')):
    return [f for d in subdirs for s in suffixes for f in glob(d + '/*' + s)]

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

    ext_modules = ext_modules,
    data_files  = data_files,
)
