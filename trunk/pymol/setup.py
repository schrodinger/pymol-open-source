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
from distutils.command.build import build
from glob import glob
import shutil
import sys, os, re

import distutils.ccompiler
import multiprocessing.pool

def CCompiler_compile(self, sources, output_dir=None, macros=None,
        include_dirs=None, debug=0, extra_preargs=None, extra_postargs=None,
        depends=None):
    '''
    Enable parallel and incremental build.

    To do a clean build, please remove the "build" directory.
    '''
    macros, objects, extra_postargs, pp_opts, build = self._setup_compile(
            output_dir, macros, include_dirs, sources, depends, extra_postargs)
    cc_args = self._get_cc_args(pp_opts, debug, extra_preargs)

    def _single_compile(obj):
        try:
            src, ext = build[obj]
        except KeyError:
            return
        try:
            if not self.force and \
                    os.path.getmtime(obj) > \
                    os.path.getmtime(src):
                return
        except OSError:
            pass
        self._compile(obj, src, ext, cc_args, extra_postargs, pp_opts)

    pmap(_single_compile, objects)
    return objects

jobs = int(os.getenv('JOBS', 0))
pmap = map if jobs == 1 else multiprocessing.pool.ThreadPool(jobs or None).map

distutils.ccompiler.CCompiler.compile = CCompiler_compile

def posix_find_lib(names, lib_dirs):
    # http://stackoverflow.com/questions/1376184/determine-if-c-library-is-installed-on-unix
    from subprocess import Popen, PIPE
    args = ["gcc", "-o", os.devnull, "-x", "c", "-"] + ["-L" + d for d in lib_dirs]
    for name in names:
        p = Popen(args + ["-l" + name], stdin=PIPE, stdout=PIPE, stderr=PIPE)
        p.communicate("int main(){}")
        if p.wait() == 0:
            return name
    raise IOError('could not find any of ' + str(names))

class install_pymol(install):
    pymol_path = None
    bundled_pmw = False
    user_options = install.user_options + [
        ('pymol-path=', None, 'PYMOL_PATH'),
        ('bundled-pmw', None, 'install bundled Pmw module'),
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
        self.make_launch_script()

        if self.bundled_pmw:
            os.system("tar -C %s -zxvf modules/pmg_tk/pmw.tgz" % self.install_libbase)

    def unchroot(self, name):
        if self.root is not None and name.startswith(self.root):
            return name[len(self.root):]
        return name

    def copy_tree_nosvn(self, src, dst):
        ignore = lambda src, names: set(['.svn']).intersection(names)
        if os.path.exists(dst):
            shutil.rmtree(dst)
        print 'copying', src, '->', dst
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
                    out.write('if [ "$DISPLAY" == "" ]; then DISPLAY=":0.0"; export DISPLAY; fi' + os.linesep)
                out.write('PYMOL_PATH="%s"; export PYMOL_PATH' % pymol_path + os.linesep)
                out.write('"%s" "%s" "$@"' % (python_exe, pymol_file) + os.linesep)

        os.chmod(launch_script, 0755)
        self.mkpath(self.install_scripts)
        self.copy(launch_script, self.install_scripts)

#============================================================================

import create_shadertext
create_shadertext.create_shadertext(
        "data/shaders",
        "shadertext.txt",
        "generated/include/ShaderText.h",
        "generated/src/ShaderText.c")

pymol_src_dirs = [
    "ov/src",
    "layer0",
    "layer1",
    "layer2",
    "layer3",
    "layer4",
    "layer5",
    "modules/cealign/src",
    "modules/cealign/src/tnt",
    'generated/src',
    'generated/include',
]

def_macros = [
    ("_PYMOL_MODULE", None),
]

libs = []
pyogl_libs = []
lib_dirs = []
ext_comp_args = []
ext_link_args = []

if True:
    # VMD plugin support
    pymol_src_dirs += [
        'contrib/uiuc/plugins/include',
        'contrib/uiuc/plugins/molfile_plugin/src',
    ]
    def_macros += [
        ("_PYMOL_VMD_PLUGINS", None),
    ]

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
            ("_PYMOL_LIBPNG",None),
            ("_PYMOL_FREETYPE",None),
            ("_PYMOL_INLINE",None),
            ("_PYMOL_NUMPY",None),
            ("_PYMOL_OPENGL_SHADERS",None),
            ("NO_MMLIBS",None),
            ("_PYMOL_CGO_DRAWARRAYS",None),
            ("_PYMOL_CGO_DRAWBUFFERS",None),
            ("_CGO_DRAWARRAYS",None),
            ("_PYMOL_GL_CALLLISTS",None),
            ("OPENGL_ES_2",None),
            ]

    libs += ["png", "freetype"]

    try:
        prefix_path = os.environ['PREFIX_PATH'].split(os.pathsep)
    except KeyError:
        prefix_path = ["/usr", "/usr/X11", "/opt/local", "/sw"]

    for prefix in prefix_path:
        inc_dirs += filter(os.path.isdir, [prefix + s for s in ["/include", "/include/freetype2"]])
        lib_dirs += filter(os.path.isdir, [prefix + s for s in ["/lib"]])

    glut = posix_find_lib(['glut', 'freeglut'], lib_dirs)
    for _libs in (libs, pyogl_libs):
        _libs += ["GL", "GLU", "GLEW", glut]

    ext_comp_args = ["-ffast-math", "-funroll-loops", "-O3", "-fcommon"]

def get_pymol_version():
    return re.findall(r'_PyMOL_VERSION "(.*)"', open('layer0/Version.h').read())[0]

def get_sources(subdirs, suffixes=('.c', '.cpp')):
    return [f for d in subdirs for s in suffixes for f in glob(d + '/*' + s)]

def get_packages(base, parent='', r=[]):
    from os.path import join, exists
    if parent:
        r.append(parent)
    for name in os.listdir(join(base, parent)):
        if '.' not in name and exists(join(base, parent, name, '__init__.py')):
            get_packages(base, join(parent, name))
    return r

def pyogl_extension(name, sources):
    return Extension(name, sources, inc_dirs, def_macros, None, lib_dirs, pyogl_libs,
            extra_compile_args=ext_comp_args, extra_link_args=ext_link_args)

distribution = setup ( # Distribution meta-data
    cmdclass  = {'install': install_pymol},
    name      = "pymol",
    version   = get_pymol_version(),
    author    = "Schrodinger",
    url       = "http://pymol.org",
    contact   = "pymol-users@lists.sourceforge.net",
    description = "PyMOL is a Python-enhanced molecular graphics tool. It excels at 3D visualization of proteins, small molecules, density, surfaces, and trajectories. It also includes molecular editing, ray tracing, and movies. Open Source PyMOL is free to everyone!", 

    package_dir = {'' : 'modules'},
    packages = get_packages('modules'),

    ext_modules = [
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

        pyogl_extension("pymol.opengl.glu._glu_num", ["contrib/pyopengl/_glu_nummodule.c"]),
        pyogl_extension("pymol.opengl.glu._glu", ["contrib/pyopengl/_glumodule.c"]),
        pyogl_extension("pymol.opengl.glut._glut", ["contrib/pyopengl/_glutmodule.c"]),
        pyogl_extension("pymol.opengl.gl._opengl_num", ["contrib/pyopengl/_opengl_nummodule.c"]),
        pyogl_extension("pymol.opengl.gl._opengl", ["contrib/pyopengl/_openglmodule.c"]),
        pyogl_extension("pymol.opengl.gl.openglutil", ["contrib/pyopengl/openglutil.c"]),
        pyogl_extension("pymol.opengl.gl.openglutil_num", ["contrib/pyopengl/openglutil_num.c"]),
    ],
)
