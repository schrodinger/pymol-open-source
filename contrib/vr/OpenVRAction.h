/*
OpenVR for PyMOL Copyright Notice
=====================================

The OpenVR for PyMOL source code is copyrighted, but you can freely use and
copy it as long as you don't change or remove any of the Copyright notices.
OpenVR for PyMOL is made available under the following open-source license
terms:

------------------------------------------------------------------------------
Copyright (c) 2018 EPAM Systems, Inc.

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
------------------------------------------------------------------------------

*/

#ifndef _H_OpenVRAction
#define _H_OpenVRAction

// system headers
#include "openvr.h"

struct OpenVRAction {
  char const* path;

  enum Type {
    TYPE_DIGITAL,
    TYPE_ANALOG,
    TYPE_POSE,
  } type;

  vr::VRActionHandle_t handle;
  vr::VRInputValueHandle_t restrictToDevice;

  union {
    struct {
      bool bActive;
      vr::VRInputValueHandle_t activeOrigin;
    };
    vr::InputDigitalActionData_t digital;
    vr::InputAnalogActionData_t analog;
    vr::InputPoseActionData_t pose;
  };

  vr::InputOriginInfo_t origin;

  explicit OpenVRAction(vr::IVRInput* Input, char const* path, Type type = TYPE_DIGITAL)
    : path(path)
    , type(type)
    , handle(vr::k_ulInvalidActionHandle)
    , restrictToDevice(vr::k_ulInvalidInputValueHandle)
  {
    bActive = false;
    Input->GetActionHandle(path, &handle);
  }

  bool IsPressed() const {
    return bActive && digital.bState;
  }

  bool WasPressed() const {
    return bActive && digital.bChanged && digital.bState;
  }

  bool WasReleased() const {
    return bActive && digital.bChanged && !digital.bState;
  }

  bool WasPressedOrReleased() const {
    return bActive && digital.bChanged;
  }

  bool PoseValid() const {
    return bActive && pose.pose.bPoseIsValid;
  }

  unsigned DeviceIndex() const {
    return origin.trackedDeviceIndex;
  }

  void Update(vr::IVRInput* Input) {
    vr::EVRInputError error = vr::VRInputError_WrongType;
    if (handle != vr::k_ulInvalidActionHandle) {
      switch (type) {
        case TYPE_DIGITAL:
          error = Input->GetDigitalActionData(handle, &digital, sizeof(digital), restrictToDevice);
          break;
        case TYPE_ANALOG:
          error = Input->GetAnalogActionData(handle, &analog, sizeof(analog), restrictToDevice);
          break;
        case TYPE_POSE:
          error =
#if defined(_PYMOL_OPENVR_MINOR) && _PYMOL_OPENVR_MINOR >= 4
              // OpenVR SDK 1.4.18: Split GetPoseActionData into two functions
              Input->GetPoseActionDataRelativeToNow
#else
              Input->GetPoseActionData
#endif
              (handle, vr::TrackingUniverseSeated, 0.0f, &pose, sizeof(pose), restrictToDevice);
          break;
        default:
          break;
      }
    }
    if (error != vr::VRInputError_None)
      bActive = false;
    if (bActive)
      error = Input->GetOriginTrackedDeviceInfo(activeOrigin, &origin, sizeof(origin));
    if (error != vr::VRInputError_None) {
      origin.devicePath = vr::k_ulInvalidInputValueHandle;
      origin.trackedDeviceIndex = vr::k_unTrackedDeviceIndexInvalid;
      origin.rchRenderModelComponentName[0] = '\0';
    }
  }
};

#endif // _H_OpenVRAction
