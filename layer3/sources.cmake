target_sources(${TARGET_NAME} PRIVATE
    ${CMAKE_CURRENT_LIST_DIR}/AtomIterators.cpp
    ${CMAKE_CURRENT_LIST_DIR}/CifDataValueFormatter.cpp
    ${CMAKE_CURRENT_LIST_DIR}/Editor.cpp
    ${CMAKE_CURRENT_LIST_DIR}/Executive.cpp
    ${CMAKE_CURRENT_LIST_DIR}/ExecutivePython.cpp
    ${CMAKE_CURRENT_LIST_DIR}/Interactions.cpp
    ${CMAKE_CURRENT_LIST_DIR}/MaeExportHelpers.cpp
    ${CMAKE_CURRENT_LIST_DIR}/MoleculeExporter.cpp
    ${CMAKE_CURRENT_LIST_DIR}/MovieScene.cpp
    ${CMAKE_CURRENT_LIST_DIR}/PlugIOManager.cpp
    ${CMAKE_CURRENT_LIST_DIR}/RingFinder.cpp
    ${CMAKE_CURRENT_LIST_DIR}/Seeker.cpp
    ${CMAKE_CURRENT_LIST_DIR}/Selector.cpp
    ${CMAKE_CURRENT_LIST_DIR}/SelectorTmp.cpp
    ${CMAKE_CURRENT_LIST_DIR}/SpecRec.cpp
    ${CMAKE_CURRENT_LIST_DIR}/SpecRecSpecial.cpp
)

target_include_directories(${TARGET_NAME} PUBLIC
    ${CMAKE_CURRENT_LIST_DIR}
)
