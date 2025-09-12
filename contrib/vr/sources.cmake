target_sources(
  ${TARGET_NAME}
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/OpenVRController.cpp
          ${CMAKE_CURRENT_LIST_DIR}/OpenVRControllerModel.cpp
          ${CMAKE_CURRENT_LIST_DIR}/OpenVRLaser.cpp
          ${CMAKE_CURRENT_LIST_DIR}/OpenVRMenu.cpp
          ${CMAKE_CURRENT_LIST_DIR}/OpenVRMode.cpp
          ${CMAKE_CURRENT_LIST_DIR}/OpenVRQuad.cpp
          ${CMAKE_CURRENT_LIST_DIR}/OpenVRScenePicker.cpp
          ${CMAKE_CURRENT_LIST_DIR}/OpenVRStub.cpp
          ${CMAKE_CURRENT_LIST_DIR}/OpenVRStubDevice.cpp
          ${CMAKE_CURRENT_LIST_DIR}/OpenVRUtils.cpp)

target_include_directories(${TARGET_NAME} PUBLIC ${CMAKE_CURRENT_LIST_DIR})
