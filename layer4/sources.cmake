target_sources(${TARGET_NAME} PRIVATE
    ${CMAKE_CURRENT_LIST_DIR}/Cmd.cpp
    ${CMAKE_CURRENT_LIST_DIR}/Menu.cpp
    ${CMAKE_CURRENT_LIST_DIR}/PopUp.cpp
)

target_include_directories(${TARGET_NAME} PUBLIC
    ${CMAKE_CURRENT_LIST_DIR}
)
