option(PYMOL_GLUT "Disable GLUT" OFF)
option(PYMOL_LIBXML2 "Enable LIBXML2 support" OFF)
option(PYMOL_OSX_FRAMEWORKS "Enable XQuartz instead of native frameworks" ON)
option(PYMOL_USE_OPENMP "Enable OpenMP support" ON)
option(PYMOL_TESTING "Build C-level tests" OFF)
option(PYMOL_OPENVR "Enable openvr support" OFF)
option(PYMOL_VMD_PLUGINS "Disable VMD molfile plugins (libnetcdf dependency)"
       ON)
option(PYMOL_USE_VTKM "Use VTK-m for isosurface generation" OFF)
option(
  PYMOL_USE_MSGPACKC
  "c++11: use msgpack-c header-only library; c: link against shared library; no: disable fast MMTF load support"
  ON)
