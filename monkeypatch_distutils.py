'''
Enhances distutils with:
    - better C++ support
    - incremental builds (w/o header dependency handling)
    - parallel builds
'''

import os
import sys
import distutils.sysconfig
import distutils.unixccompiler

# threaded parallel map (optional)
pmap = map

if sys.platform == 'darwin':
    import _osx_support
    _osx_support._UNIVERSAL_CONFIG_VARS += ('CXXFLAGS', 'LDCXXSHARED',)
    _osx_support._COMPILER_CONFIG_VARS += ('LDCXXSHARED',)

distutils.unixccompiler.UnixCCompiler.executables.update({
    'compiler_cxx'      : ["c++"],
    'compiler_so_cxx'   : ["c++"],
    'linker_so_cxx'     : ["c++", "-shared"],
    'linker_exe_cxx'    : ["c++"],
})

def monkeypatch(orig):
    '''
    Decorator to replace a function or class method. Makes the
    unpatched function available as <patchedfunction>._super
    '''
    def wrapper(func):
        func._super = orig
        func.__name__ = orig.__name__
        parent = getattr(orig, 'im_self', None)
        if parent is None:
            parent = getattr(orig, 'im_class', None)
        if parent is None:
            parent = sys.modules[orig.__module__]
        setattr(parent, orig.__name__, func)
        return func
    return wrapper

@monkeypatch(distutils.sysconfig.customize_compiler)
def customize_compiler(compiler):
    customize_compiler._super(compiler)

    (cxx, ccshared, ldcxxshared) = \
            distutils.sysconfig.get_config_vars('CXX', 'CCSHARED', 'LDCXXSHARED')

    cxxflags = os.environ.get('CXXFLAGS', '') + ' ' + os.environ.get('CPPFLAGS', '')
    ldcxxshared = os.environ.get('LDCXXSHARED', ldcxxshared) + \
            ' ' + os.environ.get('LDFLAGS', '') + \
            ' ' + os.environ.get('CXXFLAGS', '') + \
            ' ' + os.environ.get('CPPFLAGS', '')

    cxx_cmd = cxx + ' ' + cxxflags

    # C++11 by default (c++0x works for now and supports more compiler versions)
    if '-std=' not in cxx_cmd:
        cxx_cmd += ' -std=c++0x'

    compiler.set_executables(
            compiler_cxx=cxx_cmd,
            compiler_so_cxx=cxx_cmd + ' ' + ccshared,
            linker_so_cxx=ldcxxshared,
            linker_exe_cxx=cxx)

@monkeypatch(distutils.unixccompiler.UnixCCompiler.compile)
def compile(self, sources, output_dir=None, macros=None,
        include_dirs=None, debug=0, extra_preargs=None, extra_postargs=None,
        depends=None):
    '''
    Enable parallel and incremental build.

    To do a clean build, please remove the "build" directory.
    '''
    if os.getenv('DEBUG', ''):
        debug = 1

    macros, objects, extra_postargs, pp_opts, build = self._setup_compile(
            output_dir, macros, include_dirs, sources, depends, extra_postargs)
    cc_args = self._get_cc_args(pp_opts, debug, extra_preargs)

    compiler_so = self.compiler_so
    compiler_so_cxx = self.compiler_so_cxx

    if sys.platform == 'darwin':
        compiler_so = _osx_support.compiler_fixup(compiler_so, cc_args + extra_postargs)
        compiler_so_cxx = _osx_support.compiler_fixup(compiler_so_cxx, cc_args + extra_postargs)

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

        # _compile
        compiler = compiler_so_cxx \
                if self.detect_language(src) == 'c++' \
                else compiler_so
        try:
            self.spawn(compiler + cc_args + [src, '-o', obj] + extra_postargs)
        except distutils.errors.DistutilsExecError, msg:
            raise distutils.errors.CompileError, msg

    pmap(_single_compile, objects)
    return objects
