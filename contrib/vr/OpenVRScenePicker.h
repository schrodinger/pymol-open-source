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

#ifndef _H_OpenVRScenePicker
#define _H_OpenVRScenePicker

// local headers
#include "OpenVRLaserTarget.h"
#include "OpenVRMode.h"

class OpenVRScenePicker : public OpenVRLaserTarget {
  OpenVRInputHandlers* m_inputHandlers;
  unsigned m_ownerID;
  int m_clickX;
  int m_clickY;
  bool m_active;
  float m_matrix[16];

public:
  OpenVRScenePicker();

  void Init(OpenVRInputHandlers* inputHandlers);
  void Free();

  void Activate(unsigned ownerID, int x, int y);
  void Deactivate();
  bool IsActive() const;

  bool LaserShoot(float const* origin, float const* dir, float const* color, float* distance = 0);
  void LaserClick(int glutButton, int glutState);
  bool IsLaserAllowed(unsigned deviceIndex) const;
  float const* GetLaserColor() const;

  float const* GetMatrix() const;
};

#endif // _H_OpenVRScenePicker
