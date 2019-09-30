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

// this headers
#include "OpenVRStubDevice.h"

// system headers
#include <string.h>

namespace vr {
namespace stub {

const Device_t DeviceList_t::devices[] = {{
  TrackedDeviceClass_HMD, {
    DeviceProperty_t(Prop_ModelNumber_String, "Stub HMD"),
    DeviceProperty_t(Prop_SerialNumber_String, "0DD-F00D0000"),
  },
}};

unsigned DeviceProperty_t::Get(char* buffer, unsigned bufferSize, ETrackedPropertyError* error) const {
  // validate prop type
  if (m_type != k_unStringPropertyTag) {
    if (error)
      *error = TrackedProp_WrongDataType;
    return 0;
  }

  // validate buffer size
  unsigned size = strlen(m_asString) + 1;
  if (!buffer || bufferSize < size) {
    if (error)
      *error = TrackedProp_BufferTooSmall;
    return size;
  }

  strcpy(buffer, m_asString);
  if (error)
    *error = TrackedProp_Success;
  return size;
}

DeviceProperty_t const* Device_t::FindProperty(ETrackedDeviceProperty prop) const
{
  for (unsigned i = 0; i < MAX_DEVICE_PROPERTIES && props[i].IsValid(); ++i) {
    if (props[i].Is(prop)) {
      return &props[i];
    }
  }
  return nullptr;
}

Device_t const* DeviceList_t::GetDevice(unsigned index)
{
  return index >= MAX_DEVICES || devices[index].deviceClass == TrackedDeviceClass_Invalid ? nullptr : &devices[index];
}

} // stub
} // vr 

