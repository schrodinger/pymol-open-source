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

#ifndef _H_OpenVRStubDevice
#define _H_OpenVRStubDevice

// system headers
#include "openvr.h"

namespace vr {
namespace stub {

class DeviceProperty_t {
  ETrackedDeviceProperty m_prop;
  PropertyTypeTag_t m_type;
  union {
    char const* m_asString;
  };

public:
  DeviceProperty_t() : m_prop(Prop_Invalid) {}
  DeviceProperty_t(ETrackedDeviceProperty prop, char const* value) : m_prop(prop), m_type(k_unStringPropertyTag), m_asString(value) {}

  bool Is(ETrackedDeviceProperty prop) const {
    return m_prop == prop;
  }

  bool IsValid() const {
    return m_prop != Prop_Invalid;
  }

  unsigned Get(char* buffer, unsigned bufferSize, ETrackedPropertyError* error) const;
};

struct Device_t {
  static const unsigned MAX_DEVICE_PROPERTIES = 2;

  // should be a plain old struct w/o constructor for easier initialization
  ETrackedDeviceClass deviceClass;
  DeviceProperty_t props[MAX_DEVICE_PROPERTIES];

  DeviceProperty_t const* FindProperty(ETrackedDeviceProperty prop) const;
};

struct DeviceList_t {
  static const unsigned MAX_DEVICES = 1;

  static const Device_t devices[MAX_DEVICES];

  static Device_t const* GetDevice(unsigned index);
};

} // stub
} // vr 

#endif /* _H_OpenVRStubDevice */
