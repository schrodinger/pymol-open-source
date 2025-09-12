find_package(OpenMP)
list(APPEND PYMOL_DEF_MACROS "PYMOL_OPENMP")

if(OpenMP_C_FOUND)
  list(APPEND PYMOL_COMPILE_OPTIONS_C ${OpenMP_C_FLAGS})
  list(APPEND PYMOL_LIBS ${OpenMP_C_LIBRARIES})
endif()

if(OpenMP_CXX_FOUND)
  list(APPEND PYMOL_COMPILE_OPTIONS_CXX ${OpenMP_CXX_FLAGS})
  if(OpenMP_CXX_LIBRARIES AND NOT OpenMP_CXX_LIBRARIES STREQUAL
                              OpenMP_C_LIBRARIES)
    list(APPEND PYMOL_LIBS ${OpenMP_CXX_LIBRARIES})
  endif()
endif()

if(APPLE)
  list(APPEND PYMOL_COMPILE_OPTIONS "-Xpreprocessor" "-fopenmp")
  list(APPEND PYMOL_LIBS "omp")

elseif(WIN32)
  list(APPEND PYMOL_COMPILE_OPTIONS "/openmp")

else()
  list(APPEND PYMOL_COMPILE_OPTIONS "-fopenmp")
  list(APPEND PYMOL_LINK_OPTIONS "-fopenmp")

endif()
