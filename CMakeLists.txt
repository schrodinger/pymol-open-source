cmake_minimum_required(VERSION 3.15)

project(${TARGET_NAME})

set(CMAKE_VERBOSE_MAKEFILE on)
set(Python_FIND_VIRTUALENV FIRST)
set(Python_FIND_STRATEGY LOCATION)

add_library(${TARGET_NAME} SHARED ${ALL_SRC})

foreach(FOLDER_NAME IN ITEMS layer0 layer1 layer2 layer3 layer4 layer5 ov)
    include(${CMAKE_CURRENT_SOURCE_DIR}/${FOLDER_NAME}/sources.cmake)
endforeach()

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
include(Helpers)
include(Dependencies)
# Darwin, Windows, Linux
include(platform/${CMAKE_SYSTEM_NAME})
# MSVC, GCC, Clang
include(compiler/${CMAKE_CXX_COMPILER_FRONTEND_VARIANT})

target_compile_options(${TARGET_NAME} PRIVATE
    ${PYMOL_COMPILE_OPTIONS}
)

execute_process(
    COMMAND ${Python_EXECUTABLE}  -c "import sysconfig; print(sysconfig.get_config_var('EXT_SUFFIX'))"
    OUTPUT_VARIABLE PYTHON_MODULE_SUFFIX
    OUTPUT_STRIP_TRAILING_WHITESPACE
)

set_output_path(${TARGET_NAME} "${Python_SITEARCH}/pymol")

target_compile_features(${TARGET_NAME} PRIVATE cxx_std_17)

target_include_directories(${TARGET_NAME} PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/include
    ${CMAKE_CURRENT_SOURCE_DIR}/contrib/pocketfft
    ${ALL_INC_DIR}
    ${Python_INCLUDE_DIRS}
    ${Python_NumPy_INCLUDE_DIRS}
    ${PYMOL_INC_DIRS}
)

target_link_libraries(${TARGET_NAME}
    ${PYMOL_LIBS}
)

target_compile_definitions(${TARGET_NAME} PUBLIC
    ${PYMOL_DEF_MACROS}
    _PYMOL_NUMPY
)

add_subdirectory(contrib/champ)
