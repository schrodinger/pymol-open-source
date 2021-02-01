'''
Enhances distutils with:
    - better C++ support
    - incremental builds
    - parallel builds
'''

import os
import shlex
import subprocess
import sys
import distutils.sysconfig
import distutils.ccompiler
import distutils.unixccompiler
from distutils import log
from distutils.errors import (DistutilsExecError, DistutilsPlatformError,
                              CompileError, LibError, LinkError)

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


def incremental_parallel_compile(single_compile, objects, force):
    '''
    Call `single_compile` on each item in `objects` (in parallel) if
    any dependency (source file) is newer than the object file.
    '''
    mtimes = {}

    def deps(obj):
        # parse .d file (Makefile syntax)
        with open(os.path.splitext(obj)[0] + '.d') as handle:
            contents = handle.read().split(': ', 1)[-1]
            contents = contents.replace('\\\n', ' ')
            for dep in shlex.split(contents):
                yield dep

    def need_compile(obj):
        if force:
            return True

        try:
            obj_mtime = os.path.getmtime(obj)
            for dep in deps(obj):
                if dep not in mtimes:
                    mtimes[dep] = os.path.getmtime(dep)
                if obj_mtime < mtimes[dep]:
                    return True
        except EnvironmentError:
            return True

        return False

    for _ in pmap(single_compile, filter(need_compile, objects)):
        pass


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


def strip_broken_isysroot(args, start=0):
    '''
    Strip -isysroot which don't exist from args.

    Those can come from lib/python2.7/_sysconfigdata_x86_64_apple_darwin13_4_0.py
    '''
    assert isinstance(args, list)

    try:
        i = args.index('-isysroot', start)
    except ValueError:
        return

    strip_broken_isysroot(args, i + 2)

    if not os.path.isdir(args[i + 1]):
        args[i:i + 2] = []


@monkeypatch(distutils.sysconfig, 'customize_compiler')
def customize_compiler(compiler):
    # remove problematic flags
    if sys.platform == 'linux' and (
            'icpc' in os.getenv('CXX', '') or
            'clang' in os.getenv('CC', '') or
            'clang' in os.getenv('LD', '')):
        import re
        re_flto = re.compile(r'-flto\S*|-fno-semantic-interposition')
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

    # strips non-existing -isysroot
    strip_broken_isysroot(compiler_so)
    strip_broken_isysroot(compiler_so_cxx)

    if sys.platform == 'darwin' and _osx_support is not None:
        # strips duplicated -isysroot
        compiler_so = _osx_support.compiler_fixup(compiler_so, cc_args + extra_postargs)
        compiler_so_cxx = _osx_support.compiler_fixup(compiler_so_cxx, cc_args + extra_postargs)

    # generate dependency (.d) file
    cc_args.append('-MMD')

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

    incremental_parallel_compile(_single_compile, objects, self.force)

    return objects

if sys.platform.startswith('win'):
    try:
        from distutils import _msvccompiler as msvccompiler
    except ImportError:
        from distutils import msvccompiler

    @monkeypatch(msvccompiler.MSVCCompiler, 'compile')
    def compile(self, sources,
                output_dir=None, macros=None, include_dirs=None, debug=0,
                extra_preargs=None, extra_postargs=None, depends=None):
        '''
        Enable parallel and incremental build.
        '''
        if not self.initialized:
            self.initialize()
        compile_info = self._setup_compile(output_dir, macros, include_dirs,
                                           sources, depends, extra_postargs)
        macros, objects, extra_postargs, pp_opts, build = compile_info

        compile_opts = extra_preargs or []
        compile_opts.append ('/c')
        if debug:
            compile_opts.extend(self.compile_options_debug)
        else:
            compile_opts.extend(self.compile_options)

        env = dict(os.environ)
        try:
            env['PATH'] = self._paths
        except AttributeError:
            pass

        def _single_compile(obj):
            try:
                src, ext = build[obj]
            except KeyError:
                return

            add_cpp_opts = False

            if debug:
                # pass the full pathname to MSVC in debug mode,
                # this allows the debugger to find the source file
                # without asking the user to browse for it
                src = os.path.abspath(src)

            if ext in self._c_extensions:
                input_opt = "/Tc" + src
            elif ext in self._cpp_extensions:
                input_opt = "/Tp" + src
                add_cpp_opts = True
            elif ext in self._rc_extensions:
                # compile .RC to .RES file
                input_opt = src
                output_opt = "/fo" + obj
                try:
                    self.spawn([self.rc] + pp_opts +
                               [output_opt] + [input_opt])
                except DistutilsExecError as msg:
                    raise CompileError(msg)
                return
            elif ext in self._mc_extensions:
                # Compile .MC to .RC file to .RES file.
                h_dir = os.path.dirname(src)
                rc_dir = os.path.dirname(obj)
                try:
                    # first compile .MC to .RC and .H file
                    self.spawn([self.mc] +
                               ['-h', h_dir, '-r', rc_dir] + [src])
                    base, _ = os.path.splitext (os.path.basename (src))
                    rc_file = os.path.join (rc_dir, base + '.rc')
                    # then compile .RC to .RES file
                    self.spawn([self.rc] +
                               ["/fo" + obj] + [rc_file])

                except DistutilsExecError as msg:
                    raise CompileError(msg)
                return
            else:
                # how to handle this file?
                raise CompileError("Don't know how to compile %s to %s"
                                   % (src, obj))

            args = [self.cc] + compile_opts + pp_opts
            if add_cpp_opts:
                args.append('/EHsc')
            args.append(input_opt)
            args.append("/Fo" + obj)
            args.extend(extra_postargs)

            try:
                msvc_spawn_and_write_d_file(obj, src, args, env, self.dry_run)
            except DistutilsExecError as msg:
                raise CompileError(msg)

        incremental_parallel_compile(_single_compile, objects, self.force)

        return objects

    def msvc_spawn_and_write_d_file(obj, src, cmd, env, dry_run):
        '''
        Run command with /showIncludes and convert output to dependency (.d) file.
        '''
        log.info(' '.join(cmd))

        if dry_run:
            return

        process = subprocess.Popen(cmd + ['/showIncludes'], env=env,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                universal_newlines=True)
        out = process.communicate()[0]

        deps = set([src])

        for line in out.splitlines():
            if not line.startswith('Note: including file:'):
                sys.stderr.write(line + '\n')
                continue

            dep = line[21:].strip()
            dep_lower = dep.lower()
            if not (
                    # filter out system headers
                    'microsoft visual studio' in dep_lower or
                    'windows kits' in dep_lower):
                deps.add(dep)

        if process.returncode != 0:
            raise DistutilsExecError("command %r failed with exit status %d"
                    % (cmd, process.returncode))

        with open(os.path.splitext(obj)[0] + '.d', 'w') as handle:
            handle.write(': ')
            for dep in deps:
                handle.write(' \\\n')
                handle.write(dep.replace('\\', '\\\\').replace(' ', '\\ '))
