target_compile_options(${TARGET_NAME} PRIVATE
    /MP
    /std:c++17
    $<$<CONFIG:Debug>:/Z7>
)

target_link_libraries(${TARGET_NAME}
    # pyconfig.py forces linking against pythonXY.lib on MSVC
    $<$<CONFIG:Debug>:/DEBUG>
)
