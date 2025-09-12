find_package(OpenVR REQUIRED)

include(${PROJECT_SOURCE_DIR}/contrib/vr/sources.cmake)
list(APPEND PYMOL_LIBS OpenVR::OpenVR)
list(APPEND PYMOL_DEF_MACROS "_PYMOL_OPENVR")
