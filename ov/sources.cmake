target_sources(${TARGET_NAME} PRIVATE
    ${CMAKE_CURRENT_LIST_DIR}/src/OVContext.cpp
    ${CMAKE_CURRENT_LIST_DIR}/src/OVHeap.cpp
    ${CMAKE_CURRENT_LIST_DIR}/src/OVHeapArray.cpp
    ${CMAKE_CURRENT_LIST_DIR}/src/OVLexicon.cpp
    ${CMAKE_CURRENT_LIST_DIR}/src/OVOneToAny.cpp
    ${CMAKE_CURRENT_LIST_DIR}/src/OVOneToOne.cpp
    ${CMAKE_CURRENT_LIST_DIR}/src/OVRandom.cpp
    ${CMAKE_CURRENT_LIST_DIR}/src/ov_utility.cpp
)

target_include_directories(${TARGET_NAME} PUBLIC
    ${CMAKE_CURRENT_LIST_DIR}/src
)
