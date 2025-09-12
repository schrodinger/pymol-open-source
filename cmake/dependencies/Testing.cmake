find_package(Catch2 REQUIRED)

include(${PROJECT_SOURCE_DIR}/layerCTest/sources.cmake)
list(APPEND PYMOL_LIBS Catch2::Catch2)
list(APPEND PYMOL_DEF_MACROS "_PYMOL_CTEST")
