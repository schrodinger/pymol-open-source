if(PYMOL_USE_MSGPACKC STREQUAL "guess")

  find_path(
    MSGPACK_CXX_INCLUDE_DIR
    NAMES msgpack.hpp
    HINTS ${CONDA_INCLUDE_DIR})

  if(MSGPACK_CXX_INCLUDE_DIR)
    message(STATUS "Found msgpack C++11 headers: ${MSGPACK_CXX_INCLUDE_DIR}")
    set(PYMOL_USE_MSGPACKC "c++11")
  else()

    find_path(
      MSGPACK_INCLUDE_DIR
      NAMES msgpack.h
      HINTS ${CONDA_INCLUDE_DIR})

    if(MSGPACK_INCLUDE_DIR)
      message(STATUS "Found msgpack C headers: ${MSGPACK_INCLUDE_DIR}")
      set(PYMOL_USE_MSGPACKC "c")

    else()
      message(
        FATAL_ERROR
          "Could not auto-detect msgpack. Please specify PYMOL_USE_MSGPACKC as 'c' or 'c++11'"
      )
    endif()

  endif()
endif()

if(PYMOL_USE_MSGPACKC STREQUAL "c++11")
  list(APPEND PYMOL_DEF_MACROS "MMTF_MSGPACK_USE_CPP11" "MSGPACK_NO_BOOST")

  find_path(
    MSGPACK_CXX_INCLUDE_DIR
    NAMES msgpack.hpp
    HINTS ${CONDA_INCLUDE_DIR})

  if(NOT MSGPACK_CXX_INCLUDE_DIR)
    message(FATAL_ERROR "msgpack C++11 headers not found")
  endif()

  list(APPEND PYMOL_INC_DIRS ${MSGPACK_CXX_INCLUDE_DIR})

elseif(PYMOL_USE_MSGPACKC STREQUAL "c")
  find_path(
    MSGPACK_INCLUDE_DIR
    NAMES msgpack.h
    HINTS ${CONDA_INCLUDE_DIR})

  find_library(
    MSGPACK_LIBRARY
    NAMES msgpack msgpackc msgpack-c
    HINTS ${CONDA_LIBRARY_DIR})

  if(NOT MSGPACK_INCLUDE_DIR OR NOT MSGPACK_LIBRARY)
    message(FATAL_ERROR "msgpack C library not found")
  endif()

  list(APPEND PYMOL_INC_DIRS ${MSGPACK_INCLUDE_DIR})
  list(APPEND PYMOL_LIBS ${MSGPACK_LIBRARY})

else()
  message(ERROR "Unknow value for PYMOL_USE_MSGPACKC key")

endif()

include(${PROJECT_SOURCE_DIR}/contrib/mmtf-c/sources.cmake)
