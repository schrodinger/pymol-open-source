target_compile_options(${TARGET_NAME} PRIVATE
    -Werror=return-type
    -Wunused-variable
    -Wno-switch
    -Wno-narrowing
    -Wno-char-subscripts
    $<$<CONFIG:Debug>:-Og>
    $<$<NOT:$<CONFIG:Debug>>:-O3>
)

target_compile_definitions(${TARGET_NAME} PUBLIC
    # bounds checking in STL containers
    $<$<CONFIG:Debug>:_GLIBCXX_ASSERTIONS>
)
