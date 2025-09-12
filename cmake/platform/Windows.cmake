target_compile_definitions(${TARGET_NAME} PUBLIC
    WIN32
)

target_link_directories(${TARGET_NAME} PUBLIC
    ${Python_LIBRARY_DIRS}
)

target_link_libraries(${TARGET_NAME}
    ${Python_LIBRARIES}
)

