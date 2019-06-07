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

try:
    import _osx_support
    _osx_support._UNIVERSAL_CONFIG_VARS += ('CXXFLAGS', 'LDCXXSHARED',)
    _osx_support._COMPILER_CONFIG_VARS += ('LDCXXSHARED',)
except ImportError:
    _osx_support = None

distutils.unixccompiler.UnixCCompiler.executables.update({
    'compiler_cxx'      : ["c++"],
    'compiler_so_cxx'   : ["c++"],
    'linker_so_cxx'     : ["c++", "-shared"],
    'linker_exe_cxx'    : ["c++"],
})


def set_parallel_jobs(N):
    '''
    Set the number of parallel build jobs.
    N=1 : single threaded
    N=0 : use number of CPUs
    '''
    global pmap

    if N == 1:
        pmap = map
    else:
        from multiprocessing import pool
        pmap = pool.ThreadPool(N or None).map


def monkeypatch(parent, name):
    '''
    Decorator to replace a function or class method. Makes the
    unpatched function available as <patchedfunction>._super
    '''
    def wrapper(func):
        orig = getattr(parent, name)
        func._super = orig
        func.__name__ = name
        setattr(parent, name, func)
        return func
    return wrapper

@monkeypatch(distutils.sysconfig, 'customize_compiler')
def customize_compiler(compiler):
    # remove problematic flags
    if sys.platform == 'linux' and (
            'clang' in os.getenv('CC', '') or
            'clang' in os.getenv('LD', '')):
        import re
        re_flto = re.compile(r'-flto\S*')
        config_vars = distutils.sysconfig.get_config_vars()
        for (key, value) in config_vars.items():
            if re_flto.search(str(value)) is not None:
                config_vars[key] = re_flto.sub('', value)

    customize_compiler._super(compiler)

    if compiler.compiler_type != "unix":
        return

    (cxx, ccshared, ldcxxshared) = \
            distutils.sysconfig.get_config_vars('CXX', 'CCSHARED', 'LDCXXSHARED')

    cxx = os.environ.get('CXX') or cxx
    cxxflags = os.environ.get('CXXFLAGS', '') + ' ' + os.environ.get('CPPFLAGS', '')
    ldcxxshared = os.environ.get('LDCXXSHARED', ldcxxshared or '') + \
            ' ' + os.environ.get('LDFLAGS', '') + \
            ' ' + os.environ.get('CXXFLAGS', '') + \
            ' ' + os.environ.get('CPPFLAGS', '')

    cxx_cmd = cxx + ' ' + cxxflags

    # C++11 by default
    if '-std=' not in cxx_cmd:
        cxx_cmd += ' -std=c++11'

    compiler.set_executables(
            compiler_cxx=cxx_cmd,
            compiler_so_cxx=cxx_cmd + ' ' + ccshared,
            linker_so_cxx=ldcxxshared,
            linker_exe_cxx=cxx)

@monkeypatch(distutils.unixccompiler.UnixCCompiler, 'compile')
def compile(self, sources, output_dir=None, macros=None,
        include_dirs=None, debug=0, extra_preargs=None, extra_postargs=None,
        depends=None):
    '''
    Enable parallel and incremental build.

    To do a clean build, please remove the "build" directory.
    '''
    macros, objects, extra_postargs, pp_opts, build = self._setup_compile(
            output_dir, macros, include_dirs, sources, depends, extra_postargs)
    cc_args = self._get_cc_args(pp_opts, debug, extra_preargs)

    compiler_so = self.compiler_so
    compiler_so_cxx = self.compiler_so_cxx

    if sys.platform == 'darwin' and _osx_support is not None:
        compiler_so = _osx_support.compiler_fixup(compiler_so, cc_args + extra_postargs)
        compiler_so_cxx = _osx_support.compiler_fixup(compiler_so_cxx, cc_args + extra_postargs)

    mtimes = {}

    def need_compile(obj):
        try:
            src, ext = build[obj]
        except KeyError:
            return False

        if self.force:
            return True

        try:
            # parse .d file if it exists (Makefile syntax)
            deps = []
            with open(os.path.splitext(obj)[0] + '.d') as handle:
                contents = handle.read().split(': ', 1)[-1]
                for line in contents.splitlines():
                    deps.extend(line.rstrip('\n\r\\').split())
        except EnvironmentError:
            # no .d file, consider src the only dependency
            deps = [src]

        try:
            obj_mtime = os.path.getmtime(obj)
            for dep in deps:
                if dep not in mtimes:
                    mtimes[dep] = os.path.getmtime(dep)
                if obj_mtime < mtimes[dep]:
                    return True
        except OSError:
            return True

        return False

    def _single_compile(obj):
        src = build[obj][0]

        # _compile
        compiler = compiler_so_cxx \
                if self.detect_language(src) == 'c++' \
                else compiler_so
        try:
            self.spawn(compiler + cc_args + [src, '-o', obj] + extra_postargs)
        except distutils.errors.DistutilsExecError as msg:
            raise distutils.errors.CompileError(msg)

    for _ in pmap(_single_compile, filter(need_compile, objects)):
        pass

    return objects
