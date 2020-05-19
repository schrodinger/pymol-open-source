
----------------------------------------------------------------------
INSTALLATION VIA COMPILATION 
----------------------------------------------------------------------

See also: http://pymolwiki.org/index.php/Linux_Install

REQUIREMENTS

    - C++11 compiler (e.g. gcc 4.7+)
    - Python 3.6+
    - Pmw (Python Megawidgets) (optional, for legacy GUI/plugins)
      https://github.com/schrodinger/pmw-patched
    - OpenGL
    - GLEW
    - GLUT (freeglut) (optional, enable with --glut)
    - libpng
    - freetype
    - libxml2 (optional, for COLLADA export, disable with --no-libxml)
    - msgpack-c 2.1.5+ (optional, for fast MMTF loading and export,
        disable with --use-msgpackc=no)
    - mmtf-cpp (for fast MMTF export, disable with --use-msgpackc=no)
    - PyQt5, PyQt4, PySide2 or PySide (optional, will fall back to Tk
        interface if compiled with --glut)
    - glm
    - catch2 (optional, enable with --testing)
    - openvr 1.0.x (optional, enable with --openvr)
    - libnetcdf (optional, disable with --no-vmd-plugins)

SETUP OPTIONS

    python setup.py --help
    python setup.py --help-distutils
    python setup.py --help-distutils install

    Special install options:
      --pymol-path=       installation directory for PyMOL data ($PYMOL_PATH)
      --no-launcher       skip installation of the pymol launcher

    Environment variables:
      PREFIX_PATH   Colon-delimited list of paths to search for headers and
                    libraries, e.g. $HOME/mmtf-cpp:$HOME/msgpack-c:/opt/local
      CXX           C++ compiler command
      CC            C compiler command
      CXXFLAGS      C++ compiler flags
      CFLAGS        C compiler and linker flags
      CPPFLAGS      C/C++ preprocessor flags, e.g. -I/tmp/msgpack-c/include
      LDFLAGS       linker flags

INSTALLATION

    python setup.py install --prefix=~/someplace

RUNNING PyMOL

    ~/someplace/bin/pymol

Good luck!

