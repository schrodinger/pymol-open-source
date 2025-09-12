target_sources(
  ${TARGET_NAME}
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/Test_Algorithm.cpp
          ${CMAKE_CURRENT_LIST_DIR}/Test_Bezier.cpp
          ${CMAKE_CURRENT_LIST_DIR}/Test_cache_ptr.cpp
          ${CMAKE_CURRENT_LIST_DIR}/Test_CCrystal.cpp
          ${CMAKE_CURRENT_LIST_DIR}/Test_CifFile.cpp
          ${CMAKE_CURRENT_LIST_DIR}/Test_Classic_VLA.cpp
          ${CMAKE_CURRENT_LIST_DIR}/Test_copyable_ptr.cpp
          ${CMAKE_CURRENT_LIST_DIR}/Test_Event.cpp
          ${CMAKE_CURRENT_LIST_DIR}/Test_Executive.cpp
          ${CMAKE_CURRENT_LIST_DIR}/Test_Image.cpp
          ${CMAKE_CURRENT_LIST_DIR}/Test_List.cpp
          ${CMAKE_CURRENT_LIST_DIR}/Test_Picking.cpp
          ${CMAKE_CURRENT_LIST_DIR}/Test_Result.cpp
          ${CMAKE_CURRENT_LIST_DIR}/Test_ScrollBar.cpp
          ${CMAKE_CURRENT_LIST_DIR}/Test_ShaderPreprocessor.cpp
          ${CMAKE_CURRENT_LIST_DIR}/Test_SymOp.cpp
          ${CMAKE_CURRENT_LIST_DIR}/Test_Test.cpp
          ${CMAKE_CURRENT_LIST_DIR}/Test_TTT.cpp
          ${CMAKE_CURRENT_LIST_DIR}/Test_type_traits.cpp
          ${CMAKE_CURRENT_LIST_DIR}/Test_Util.cpp
          ${CMAKE_CURRENT_LIST_DIR}/Test_VLA.cpp
          ${CMAKE_CURRENT_LIST_DIR}/Test_zstring_view.cpp
          ${CMAKE_CURRENT_LIST_DIR}/Test.cpp)

target_include_directories(${TARGET_NAME} PUBLIC ${CMAKE_CURRENT_LIST_DIR})
