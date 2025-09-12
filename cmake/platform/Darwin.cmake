target_link_options(${TARGET_NAME} PRIVATE
    LINKER:-undefined,dynamic_lookup
)

target_compile_options(${TARGET_NAME} PRIVATE
    # optimization currently causes a clang segfault on OS X 10.9 when
    # compiling layer2/RepCylBond.cpp
    -fno-strict-aliasing
)

target_compile_definitions(${TARGET_NAME} PUBLIC
    PYMOL_CURVE_VALIDATE
)
