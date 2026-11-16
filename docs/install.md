# Install from source

Most people do not need this page: the
[Releases page](https://github.com/adrian-jasinski/pymolish-studio/releases)
has ready-to-run bundles. Build from source when you want to change the code
or need a system the bundles do not cover.

The build compiles the C++ core with CMake and installs the Python package
with `pip`. It takes a few minutes the first time.

## What you need

| Needed | For |
| --- | --- |
| C++17 compiler (gcc 8 or newer), Python 3.11 or newer with `pip` | the build; `pip` fetches CMake itself |
| OpenGL, GLEW, GLM, libpng, FreeType | drawing, pictures, text |
| PySide6, numpy | the window and arrays; `pip` installs them |

| Optional | For | Switch |
| --- | --- | --- |
| msgpack-c 2.1.5 or newer | BinaryCIF files | `use-msgpackc=no` builds without it |
| netCDF | trajectory formats of the VMD plugins | `vmd-plugins=no` builds without it |
| nlohmann-json 3.10 or newer | glTF and GLB export | off unless `json=True` |
| Catch2 | C++ tests | off unless `testing=True` |

Give a switch to `pip` as `--config-settings name=value`.

## Linux (Ubuntu)

These are the packages the CI build installs; add `build-essential` if you
have no compiler.

```sh
sudo apt-get install libfreetype6-dev libglew-dev libglm-dev libpng-dev \
    libmsgpack-dev libnetcdf-dev nlohmann-json3-dev catch2

python3 -m venv .venv && . .venv/bin/activate
pip install --config-settings json=True .
pymolish
```

## Windows and macOS

The native libraries come from conda-forge through
[pixi](https://pixi.sh); `pixi.toml` lists them. This is how CI builds.

```sh
pixi run -e py313 pip install --config-settings json=True .
pixi run -e py313 pymolish
```

## Working on the code

An editable install: Python edits under `modules/` take effect on the next
start. After a C++ or shader edit run the same command again; CMake rebuilds
only what changed.

```sh
pip install cmake numpy setuptools
pip install --no-build-isolation --config-settings testing=True -e '.[dev]'
```

| Task | Command |
| --- | --- |
| Run all tests, no window | `pymolish -ckqy testing/testing.py --run all` |
| Run one test file | `pymolish -ckqy testing/testing.py --run testing/tests/api/<file>.py` |
| Measure performance, see [baseline.md](baseline.md) | `pymolish -ckqy testing/bench.py` |
| Build a bundle into `dist/pymolish/` | `pyinstaller --noconfirm packaging/pyinstaller.spec` |

Useful when a build misbehaves:

- `--config-settings jobs=N` or `JOBS=N` limits parallel compile jobs.
- `PREFIX_PATH` is a colon-separated list of directories to search for
  headers and libraries, for example `$HOME/msgpack-c:/opt/local`.
- `CXX`, `CXXFLAGS` and `LDFLAGS` choose the compiler and its flags.
- `export CMAKE_CXX_COMPILER_LAUNCHER=ccache` uses ccache.

Next: [Contributing](../CONTRIBUTING.md).
