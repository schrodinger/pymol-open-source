target_sources(${TARGET_NAME} PRIVATE
    ${CMAKE_CURRENT_LIST_DIR}/mmtf_parser.cpp
)

target_include_directories(${TARGET_NAME} PUBLIC
    ${CMAKE_CURRENT_LIST_DIR}
)
