set(CHEMPY_PROJECT_NAME _champ)

add_library(${CHEMPY_PROJECT_NAME} MODULE
    champ.c
    champ_module.c
    chiral.c
    err2.c
    feedback2.c
    list.c
    os_memory.c
    sort.c
    strblock.c
    vla.c
)

include(Helpers)

target_link_options(${CHEMPY_PROJECT_NAME} PRIVATE
    "$<$<PLATFORM_ID:Darwin>:LINKER:-undefined,dynamic_lookup>"
)

target_include_directories(${CHEMPY_PROJECT_NAME} PUBLIC
    ${Python_INCLUDE_DIRS}
    ${CMAKE_CURRENT_SOURCE_DIR}
)
target_link_directories(${CHEMPY_PROJECT_NAME} PUBLIC
    $<$<CXX_COMPILER_ID:MSVC>:${Python_LIBRARY_DIRS}>
)

target_link_libraries(${CHEMPY_PROJECT_NAME} PUBLIC
    # pyconfig.py forces linking against pythonXY.lib on MSVC
    $<$<CXX_COMPILER_ID:MSVC>:${Python_LIBRARIES}>
)

set_output_path(${CHEMPY_PROJECT_NAME} "${Python_SITEARCH}/chempy/champ")

target_compile_features(${CHEMPY_PROJECT_NAME} PRIVATE cxx_std_17)
