find_package(NetCDF MODULE QUIET)

if(NetCDF_FOUND)
  list(APPEND PYMOL_LIBS NetCDF::NetCDF)
else()
  find_library(NETCDF_LIB NAMES netcdf)
  find_path(NETCDF_INCLUDE_DIR NAMES netcdf.h)

  if(NETCDF_LIB AND NETCDF_INCLUDE_DIR)
    list(APPEND PYMOL_LIBS ${NETCDF_LIB})
    list(APPEND PYMOL_INC_DIRS ${NETCDF_INCLUDE_DIR})
  else()
    message(FATAL_ERROR "NetCDF not found!")
  endif()
endif()

include(${PROJECT_SOURCE_DIR}/contrib/uiuc/plugins/sources.cmake)
list(APPEND PYMOL_DEF_MACROS "_PYMOL_VMD_PLUGINS")
