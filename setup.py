#!/usr/bin/env python
#
# This script only applies if you are performing a Python setuptools-based
# installation of PyMOL.
#
# It may assume that all of PyMOL's external dependencies are
# pre-installed into the system.

import argparse
import glob
import io as cStringIO
import os
import pathlib
import re
import shutil
import sys
import sysconfig
import time
from collections import defaultdict
from subprocess import PIPE, Popen

import numpy
from setuptools import Extension, setup
from setuptools.command.build_ext import build_ext
from setuptools.command.build_py import build_py
from setuptools.command.install import install

# non-empty DEBUG variable turns off optimization and adds -g flag
DEBUG = bool(os.getenv("DEBUG", ""))
WIN = sys.platform.startswith("win")
MAC = sys.platform.startswith("darwin")


# Have to copy from "create_shadertext.py" script due to the use of pyproject.toml
# Full explanation:
# https://github.com/pypa/setuptools/issues/3939
def create_all(generated_dir, pymoldir="."):
    """
    Generate various stuff
    """
    create_shadertext(
        os.path.join(pymoldir, "data", "shaders"),
        generated_dir,
        os.path.join(generated_dir, "ShaderText.h"),
        os.path.join(generated_dir, "ShaderText.cpp"),
    )
    create_buildinfo(generated_dir, pymoldir)


class openw(object):
    """
    File-like object for writing files. File is actually only
    written if the content changed.
    """

    def __init__(self, filename):
        if os.path.exists(filename):
            self.out = cStringIO.StringIO()
            self.filename = filename
        else:
            os.makedirs(os.path.dirname(filename), exist_ok=True)
            self.out = open(filename, "w")
            self.filename = None

    def close(self):
        if self.out.closed:
            return
        if self.filename:
            with open(self.filename) as handle:
                oldcontents = handle.read()
            newcontents = self.out.getvalue()
            if oldcontents != newcontents:
                self.out = open(self.filename, "w")
                self.out.write(newcontents)
        self.out.close()

    def __getattr__(self, name):
        return getattr(self.out, name)

    def __enter__(self):
        return self

    def __exit__(self, *a, **k):
        self.close()

    def __del__(self):
        self.close()


def create_shadertext(shaderdir, shaderdir2, outputheader, outputfile):
    outputheader = openw(outputheader)
    outputfile = openw(outputfile)

    include_deps = defaultdict(set)
    ifdef_deps = defaultdict(set)

    # get all *.gs *.vs *.fs *.shared from the two input directories
    shaderfiles = set()
    for sdir in [shaderdir, shaderdir2]:
        for ext in ["gs", "vs", "fs", "shared", "tsc", "tse"]:
            shaderfiles.update(
                map(os.path.basename, sorted(glob.glob(os.path.join(sdir, "*." + ext))))
            )

    varname = "_shader_cache_raw"
    outputheader.write("extern const char * %s[];\n" % varname)
    outputfile.write("const char * %s[] = {\n" % varname)

    for filename in sorted(shaderfiles):
        shaderfile = os.path.join(shaderdir, filename)
        if not os.path.exists(shaderfile):
            shaderfile = os.path.join(shaderdir2, filename)

        with open(shaderfile, "r") as handle:
            contents = handle.read()

        if True:
            outputfile.write('"%s", ""\n' % (filename))

            for line in contents.splitlines():
                line = line.strip()

                # skip blank lines and obvious comments
                if not line or line.startswith("//") and not "*/" in line:
                    continue

                # write line, quoted, escaped and with a line feed
                outputfile.write(
                    '"%s\\n"\n' % line.replace("\\", "\\\\").replace('"', r"\"")
                )

                # include and ifdef dependencies
                if line.startswith("#include"):
                    include_deps[line.split()[1]].add(filename)
                elif line.startswith("#ifdef") or line.startswith("#ifndef"):
                    ifdef_deps[line.split()[1]].add(filename)

            outputfile.write(",\n")

    outputfile.write("0};\n")

    # include and ifdef dependencies
    for varname, deps in [("_include_deps", include_deps), ("_ifdef_deps", ifdef_deps)]:
        outputheader.write("extern const char * %s[];\n" % varname)
        outputfile.write("const char * %s[] = {\n" % varname)
        for name, itemdeps in deps.items():
            outputfile.write('"%s", "%s", 0,\n' % (name, '", "'.join(sorted(itemdeps))))
        outputfile.write("0};\n")

    outputheader.close()
    outputfile.close()


def create_buildinfo(outputdir, pymoldir="."):
    try:
        sha = (
            Popen(["git", "rev-parse", "HEAD"], cwd=pymoldir, stdout=PIPE)
            .stdout.read()
            .strip()
            .decode()
        )
    except OSError:
        sha = ""

    with openw(os.path.join(outputdir, "PyMOLBuildInfo.h")) as out:
        print(
            """
#define _PyMOL_BUILD_DATE %d
#define _PYMOL_BUILD_GIT_SHA "%s"
        """
            % (time.time(), sha),
            file=out,
        )


# handle extra arguments
def str2bool(v: str) -> bool:
    if v.lower() == "true":
        return True
    elif v.lower() == "false":
        return False
    else:
        raise argparse.ArgumentTypeError("Boolean value expected.")


class options:
    osx_frameworks = True
    jobs = int(os.getenv("JOBS", 0))
    no_libxml = False
    no_glut = True
    use_msgpackc = "guess"
    testing = False
    openvr = False
    use_openmp = "no" if MAC else "yes"
    use_vtkm = "no"
    vmd_plugins = True


parser = argparse.ArgumentParser()
parser.add_argument(
    "--glut", dest="no_glut", type=str2bool, help="link with GLUT (legacy GUI)"
)
parser.add_argument(
    "--no-osx-frameworks",
    dest="osx_frameworks",
    help="on MacOS use XQuartz instead of native frameworks",
    type=str2bool,
)
parser.add_argument(
    "--jobs",
    "-j",
    type=int,
    help="for parallel builds " "(defaults to number of processors)",
)
parser.add_argument(
    "--no-libxml",
    type=str2bool,
    help="skip libxml2 dependency, disables COLLADA export",
)
parser.add_argument("--use-openmp", choices=("yes", "no"), help="Use OpenMP")
parser.add_argument(
    "--use-vtkm",
    choices=("1.5", "1.6", "1.7", "no"),
    help="Use VTK-m for isosurface generation",
)
parser.add_argument(
    "--use-msgpackc",
    choices=("c++11", "c", "guess", "no"),
    help="c++11: use msgpack-c header-only library; c: link against "
    "shared library; no: disable fast MMTF load support",
)
parser.add_argument("--testing", type=str2bool, help="Build C-level tests")
parser.add_argument("--openvr", dest="openvr", type=str2bool)
parser.add_argument(
    "--no-vmd-plugins",
    dest="vmd_plugins",
    type=str2bool,
    help="Disable VMD molfile plugins (libnetcdf dependency)",
)
options, sys.argv[1:] = parser.parse_known_args(namespace=options)


def get_prefix_path() -> list[str]:
    """
    Return a list of paths which will be searched for "include",
    "include/freetype2", "lib", "lib64" etc.
    """
    paths = []

    if (prefix_path := os.environ.get("PREFIX_PATH")) is not None:
        paths += prefix_path.split(os.pathsep)

    if sys.platform.startswith("freebsd"):
        paths += ["/usr/local"]

    if not options.osx_frameworks:
        paths += ["usr/X11"]

    if MAC:
        for prefix in ["/sw", "/opt/local", "/usr/local"]:
            if sys.base_prefix.startswith(prefix):
                paths += [prefix]

    if is_conda_env():
        if WIN:
            if "CONDA_PREFIX" in os.environ:
                paths += [os.path.join(os.environ["CONDA_PREFIX"], "Library")]
            paths += [os.path.join(sys.prefix, "Library")]

        paths += [sys.prefix] + paths

    paths += ["/usr"]

    return paths


def is_conda_env():
    return (
        "conda" in sys.prefix
        or "conda" in sys.version
        or "Continuum" in sys.version
        or sys.prefix == os.getenv("CONDA_PREFIX")
    )


def guess_msgpackc():
    for prefix in prefix_path:
        for suffix in ["h", "hpp"]:
            f = os.path.join(prefix, "include", "msgpack", f"version_master.{suffix}")

            try:
                m = re.search(r"MSGPACK_VERSION_MAJOR\s+(\d+)", open(f).read())
            except EnvironmentError:
                continue

            if m is not None:
                major = int(m.group(1))
                if major > 1:
                    return "c++11"

    return "no"


class CMakeExtension(Extension):

    def __init__(
        self,
        name,
        sources,
        include_dirs=[],
        libraries=[],
        library_dirs=[],
        define_macros=[],
        extra_link_args=[],
        extra_compile_args=[],
    ):
        # don't invoke the original build_ext for this special extension
        super().__init__(name, sources=[])
        self.sources = sources
        self.include_dirs = include_dirs
        self.libraries = libraries
        self.library_dirs = library_dirs
        self.define_macros = define_macros
        self.extra_link_args = extra_link_args
        self.extra_compile_args = extra_compile_args


class build_ext_pymol(build_ext):
    def initialize_options(self) -> None:
        super().initialize_options()
        if DEBUG and not WIN:
            self.debug = False

    def run(self):
        for ext in self.extensions:
            self.build_cmake(ext)

    def build_cmake(self, ext):
        cwd = pathlib.Path().absolute()

        # these dirs will be created in build_py, so if you don't have
        # any python sources to bundle, the dirs will be missing
        name_split = ext.name.split(".")
        target_name = name_split[-1]
        build_temp = pathlib.Path(self.build_temp) / target_name
        build_temp.mkdir(parents=True, exist_ok=True)
        extdir = pathlib.Path(self.get_ext_fullpath(ext.name))
        extdirabs = extdir.absolute()

        extdir.parent.mkdir(parents=True, exist_ok=True)

        def concat_paths(paths):
            return "".join(path.replace("\\", "/") + ";" for path in paths)

        config = "Debug" if DEBUG else "Release"
        lib_output_dir = str(extdir.parent.absolute())
        all_files = ext.sources
        all_src = concat_paths(all_files)
        all_defs = "".join(mac[0] + ";" for mac in ext.define_macros)
        all_libs = "".join(f"{lib};" for lib in ext.libraries)
        all_ext_link = " ".join(ext.extra_link_args)
        all_comp_args = "".join(f"{arg};" for arg in ext.extra_compile_args)
        all_lib_dirs = concat_paths(ext.library_dirs)
        all_inc_dirs = concat_paths(ext.include_dirs)

        lib_mode = "RUNTIME" if WIN else "LIBRARY"

        shared_suffix = sysconfig.get_config_var("EXT_SUFFIX")

        cmake_args = [
            f"-DTARGET_NAME={target_name}",
            f"-DCMAKE_{lib_mode}_OUTPUT_DIRECTORY={lib_output_dir}",
            f"-DCMAKE_BUILD_TYPE={config}",
            f"-DALL_INC_DIR={all_inc_dirs}",
            f"-DALL_SRC={all_src}",
            f"-DALL_DEF={all_defs}",
            f"-DALL_LIB_DIR={all_lib_dirs}",
            f"-DALL_LIB={all_libs}",
            f"-DALL_COMP_ARGS={all_comp_args}",
            f"-DALL_EXT_LINK={all_ext_link}",
            f"-DSHARED_SUFFIX={shared_suffix}",
        ]

        # example of build args
        build_args = ["--config", config]
        if not WIN:  # Win /MP flag on compilation level
            cpu_count = os.cpu_count() or 1
            build_args += [f"-j{cpu_count}"]

        os.chdir(str(build_temp))
        self.spawn(["cmake", str(cwd)] + cmake_args)
        if not self.dry_run:
            self.spawn(["cmake", "--build", "."] + build_args)

        if WIN:
            # Move up from VS release folder
            cmake_lib_loc = pathlib.Path(
                lib_output_dir, "Release", f"{target_name}{shared_suffix}"
            )
            if cmake_lib_loc.exists():
                shutil.move(cmake_lib_loc, extdirabs)

        # Troubleshooting: if fail on line above then delete all possible
        # temporary CMake files including "CMakeCache.txt" in top level dir.
        os.chdir(str(cwd))


class build_py_pymol(build_py):
    def run(self):
        build_py.run(self)


class install_pymol(install):
    pymol_path = None
    bundled_pmw = False
    no_launcher = False

    user_options = install.user_options + [
        ("pymol-path=", None, "PYMOL_PATH"),
        ("bundled-pmw", None, "install bundled Pmw module"),
        ("no-launcher", None, "skip installation of the pymol launcher"),
    ]

    def finalize_options(self):
        install.finalize_options(self)

        self.pymol_path_is_default = self.pymol_path is None

        if self.pymol_path is None:
            self.pymol_path = os.path.join(self.install_libbase, "pymol", "pymol_path")
        elif self.root is not None:
            self.pymol_path = install_pymol.change_root(self.root, self.pymol_path)

    def run(self):
        super().run()
        self.install_pymol_path()

        if not self.no_launcher:
            self.make_launch_script()

        if self.bundled_pmw:
            raise Exception(
                "--bundled-pmw has been removed, please install Pmw from "
                "https://github.com/schrodinger/pmw-patched"
            )

    def unchroot(self, name):
        if self.root is not None and name.startswith(self.root):
            return name[len(self.root) :]
        return name

    def copy_tree_nosvn(self, src, dst):
        def ignore(src, names):
            return set([]).intersection(names)

        if os.path.exists(dst):
            shutil.rmtree(dst)
        print("copying %s -> %s" % (src, dst))
        shutil.copytree(src, dst, ignore=ignore)

    def copy(self, src, dst):
        copy = self.copy_tree_nosvn if os.path.isdir(src) else self.copy_file
        copy(src, dst)

    def install_pymol_path(self):
        self.mkpath(self.pymol_path)
        for name in [
            "LICENSE",
            "data",
            "test",
            "examples",
        ]:
            self.copy(name, os.path.join(self.pymol_path, name))

        if options.openvr:
            self.copy(
                "contrib/vr/README.md", os.path.join(self.pymol_path, "README-VR.txt")
            )

    def make_launch_script(self):
        if sys.platform.startswith("win"):
            launch_script = "pymol.bat"
        else:
            launch_script = "pymol"

        self.mkpath(self.install_scripts)
        launch_script = os.path.join(self.install_scripts, launch_script)

        python_exe = os.path.abspath(sys.executable)
        site_packages_dir = sysconfig.get_path('purelib')
        pymol_file = self.unchroot(
            os.path.join(site_packages_dir, "pymol", "__init__.py")
        )
        pymol_path = self.unchroot(self.pymol_path)

        with open(launch_script, "w") as out:
            if WIN:
                # paths relative to launcher, if possible
                try:
                    python_exe = "%~dp0\\" + os.path.relpath(
                        python_exe, self.install_scripts
                    )
                except ValueError:
                    pass
                try:
                    pymol_file = "%~dp0\\" + os.path.relpath(
                        pymol_file, self.install_scripts
                    )
                except ValueError:
                    pymol_file = os.path.abspath(pymol_file)

                if not self.pymol_path_is_default:
                    out.write(f"set PYMOL_PATH={pymol_path}" + os.linesep)
                out.write('"%s" "%s"' % (python_exe, pymol_file))
                out.write(" %*" + os.linesep)
            else:
                out.write("#!/bin/sh" + os.linesep)
                if not self.pymol_path_is_default:
                    out.write(f'export PYMOL_PATH="{pymol_path}"' + os.linesep)
                out.write('exec "%s" "%s" "$@"' % (python_exe, pymol_file) + os.linesep)

        os.chmod(launch_script, 0o755)


# ============================================================================


# should be something like (build_base + "/generated"), but that's only
# known to build and install instances
generated_dir = os.path.join(os.environ.get("PYMOL_BLD", "build"), "generated")

create_all(generated_dir)

# can be changed with environment variable PREFIX_PATH
prefix_path = get_prefix_path()

inc_dirs = [
    "include",
]

pymol_src_dirs = [
    "ov/src",
    "layer0",
    "layer1",
    "layer2",
    "layer3",
    "layer4",
    "layer5",
    generated_dir,
]

def_macros = [
    ("_PYMOL_LIBPNG", None),
    ("_PYMOL_FREETYPE", None),
]

if DEBUG and not WIN:
    def_macros += [
        # bounds checking in STL containers
        ("_GLIBCXX_ASSERTIONS", None),
    ]

libs = ["png", "freetype"]
lib_dirs = []
ext_comp_args = (
    [
        "-Werror=return-type",
        "-Wunused-variable",
        "-Wno-switch",
        "-Wno-narrowing",
        # legacy stuff
        "-Wno-char-subscripts",
        # optimizations
        "-Og" if DEBUG else "-O3",
    ]
    if not WIN
    else ["/MP"]
)
ext_link_args = []
ext_objects = []
data_files = []
ext_modules = []

if options.use_openmp == "yes":
    def_macros += [
        ("PYMOL_OPENMP", None),
    ]
    if MAC:
        ext_comp_args += ["-Xpreprocessor", "-fopenmp"]
        libs += ["omp"]
    elif WIN:
        ext_comp_args += ["/openmp"]
    else:
        ext_comp_args += ["-fopenmp"]
        ext_link_args += ["-fopenmp"]

if options.vmd_plugins:
    # VMD plugin support
    inc_dirs += [
        "contrib/uiuc/plugins/include",
    ]
    pymol_src_dirs += [
        "contrib/uiuc/plugins/molfile_plugin/src",
    ]
    def_macros += [
        ("_PYMOL_VMD_PLUGINS", None),
    ]

if not options.no_libxml:
    # COLLADA support
    def_macros += [("_HAVE_LIBXML", None)]
    libs += ["xml2"]

if options.use_msgpackc == "guess":
    options.use_msgpackc = guess_msgpackc()

if options.use_msgpackc == "no":
    def_macros += [("_PYMOL_NO_MSGPACKC", None)]
else:
    if options.use_msgpackc == "c++11":
        def_macros += [
            ("MMTF_MSGPACK_USE_CPP11", None),
            ("MSGPACK_NO_BOOST", None),
        ]
    else:
        libs += ["msgpackc"]

    pymol_src_dirs += ["contrib/mmtf-c"]

if options.no_glut:
    def_macros += [
        ("_PYMOL_NO_MAIN", None),
    ]

if options.testing:
    pymol_src_dirs += ["layerCTest"]
    def_macros += [("_PYMOL_CTEST", None)]

if options.openvr:
    def_macros += [("_PYMOL_OPENVR", None)]
    pymol_src_dirs += [
        "contrib/vr",
    ]

inc_dirs += pymol_src_dirs

# ============================================================================
if MAC:
    libs += ["GLEW"]
    def_macros += [("PYMOL_CURVE_VALIDATE", None)]

    if options.osx_frameworks:
        ext_link_args += [
            "-framework OpenGL",
        ] + (not options.no_glut) * [
            "-framework GLUT",
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

    def_macros += [
        ("WIN32", None),
    ]

    libs += [
        "Advapi32",  # Registry (RegCloseKey etc.)
        "Ws2_32",  # htonl
    ]

    libs += (
        [
            "glew32",
            "freetype",
            "libpng",
        ]
        + (not options.no_glut)
        * [
            "freeglut",
        ]
        + (not options.no_libxml)
        * [
            "libxml2",
        ]
    )

    if DEBUG:
        ext_comp_args += ["/Z7"]
        ext_link_args += ["/DEBUG"]

    libs += [
        "opengl32",
    ]
    # TODO: Remove when we move to setup-CMake
    ext_comp_args += ["/std:c++17"]

if not (MAC or WIN):
    libs += [
        "GL",
        "GLEW",
    ] + (not options.no_glut) * [
        "glut",
    ]

if options.use_vtkm != "no":
    for prefix in prefix_path:
        vtkm_inc_dir = os.path.join(prefix, "include", f"vtkm-{options.use_vtkm}")
        if os.path.exists(vtkm_inc_dir):
            break
    else:
        raise LookupError(
            "VTK-m headers not found." f' PREFIX_PATH={":".join(prefix_path)}'
        )
    def_macros += [
        ("_PYMOL_VTKM", None),
    ]
    inc_dirs += [
        vtkm_inc_dir,
        vtkm_inc_dir + "/vtkm/thirdparty/diy/vtkmdiy/include",
        vtkm_inc_dir + "/vtkm/thirdparty/lcl/vtkmlcl",
    ] + (options.use_vtkm == "1.5") * [
        vtkm_inc_dir + "/vtkm/thirdparty/diy",
        vtkm_inc_dir + "/vtkm/thirdparty/taotuple",
    ]
    libs += [
        f"vtkm_cont-{options.use_vtkm}",
        (
            f"vtkm_filter-{options.use_vtkm}"
            if options.use_vtkm == "1.5"
            else f"vtkm_filter_contour-{options.use_vtkm}"
        ),
    ]

if options.vmd_plugins:
    libs += [
        "netcdf",
    ]

if options.openvr:
    libs += [
        "openvr_api",
    ]

inc_dirs += [
    numpy.get_include(),
]
def_macros += [
    ("_PYMOL_NUMPY", None),
]

for prefix in prefix_path:
    for dirs, suffixes in [
        [
            inc_dirs,
            [
                ("include",),
                ("include", "freetype2"),
                ("include", "libxml2"),
                ("include", "openvr"),
            ],
        ],
        [lib_dirs, [("lib64",), ("lib",)]],
    ]:
        dirs.extend(filter(os.path.isdir, [os.path.join(prefix, *s) for s in suffixes]))

# optimization currently causes a clang segfault on OS X 10.9 when
# compiling layer2/RepCylBond.cpp
if MAC:
    ext_comp_args += ["-fno-strict-aliasing"]


def get_pymol_version():
    return re.findall(r'_PyMOL_VERSION "(.*)"', open("layer0/Version.h").read())[0]


def get_sources(subdirs, suffixes=(".c", ".cpp")):
    return sorted(
        [f for d in subdirs for s in suffixes for f in glob.glob(d + "/*" + s)]
    )


def get_packages(base, parent="", r=None):
    from os.path import exists, join

    if r is None:
        r = []
    if parent:
        r.append(parent)
    for name in os.listdir(join(base, parent)):
        if "." not in name and exists(join(base, parent, name, "__init__.py")):
            get_packages(base, join(parent, name), r)
    return r


package_dir = dict(
    (x, os.path.join(base, x)) for base in ["modules"] for x in get_packages(base)
)

# Python includes
inc_dirs.append(sysconfig.get_paths()["include"])
inc_dirs.append(sysconfig.get_paths()["platinclude"])

champ_inc_dirs = ["contrib/champ"]
champ_inc_dirs.append(sysconfig.get_paths()["include"])
champ_inc_dirs.append(sysconfig.get_paths()["platinclude"])

if WIN:
    # pyconfig.py forces linking against pythonXY.lib on MSVC
    py_lib = pathlib.Path(sysconfig.get_paths()["stdlib"]).parent / "libs"
    lib_dirs.append(str(py_lib))

ext_modules += [
    CMakeExtension(
        name="pymol._cmd",
        sources=get_sources(pymol_src_dirs),
        include_dirs=inc_dirs,
        libraries=libs,
        library_dirs=lib_dirs,
        define_macros=def_macros,
        extra_link_args=ext_link_args,
        extra_compile_args=ext_comp_args,
    ),
    CMakeExtension(
        name="chempy.champ._champ",
        sources=get_sources(["contrib/champ"]),
        include_dirs=champ_inc_dirs,
        library_dirs=lib_dirs,
    ),
]

setup(
    cmdclass={
        "build_ext": build_ext_pymol,
        "build_py": build_py_pymol,
        "install": install_pymol,
    },
    version=get_pymol_version(),
    ext_modules=ext_modules,
)
