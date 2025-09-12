find_package(Freetype REQUIRED)

list(APPEND PYMOL_LIBS Freetype::Freetype)
list(APPEND PYMOL_DEF_MACROS "_PYMOL_FREETYPE")
