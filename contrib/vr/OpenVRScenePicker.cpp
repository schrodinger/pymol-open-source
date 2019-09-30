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
#include "OpenVRScenePicker.h"

// system headers
#include "os_gl.h"
#include <math.h>
#include <string.h>

// local headers
#include "OpenVRUtils.h"

OpenVRScenePicker::OpenVRScenePicker()
: m_inputHandlers(0)
, m_ownerID(~0u)
, m_clickX(0)
, m_clickY(0)
, m_active(false)
{
  memset(m_matrix, 0, sizeof(m_matrix));
  m_matrix[0] = m_matrix[5] = m_matrix[10] = m_matrix[15] = 1.0f;
}

void OpenVRScenePicker::Init(OpenVRInputHandlers* inputHandlers)
{
  m_inputHandlers = inputHandlers;
}

void OpenVRScenePicker::Free()
{
}

void OpenVRScenePicker::Activate(unsigned ownerID, int x, int y)
{
  m_ownerID = ownerID;
  m_clickX = x;
  m_clickY = y;
  m_active = true;
}

void OpenVRScenePicker::Deactivate()
{
  m_ownerID = ~0u;
  m_active = false;
}

bool OpenVRScenePicker::IsActive() const
{
  return m_active;
}

bool OpenVRScenePicker::LaserShoot(float const* origin, float const* dir, float const* color, float* distance /* = 0 */)
{
  float matrix[16];
  float* right = &matrix[0];  // head coordinates X-axis
  float* up = &matrix[4];     // head coordinates Y-axis
  float* back = &matrix[8];   // head coordinates Z-axis
  float* pivot = &matrix[12]; // head coordinates origin

  pivot[0] = origin[0];
  pivot[1] = origin[1];
  pivot[2] = origin[2];
  pivot[3] = 1.0f;

  back[0] = -dir[0];
  back[1] = -dir[1];
  back[2] = -dir[2];
  back[3] = 0.0f;

  if (fabsf(back[1]) < 0.95f) {
    up[0] = 0.0f;
    up[1] = 1.0f;
    up[2] = 0.0f;
    up[3] = 0.0f;
    OpenVRUtils::VectorCrossProduct(up, back, right);
    OpenVRUtils::VectorNormalize(right);
    OpenVRUtils::VectorCrossProduct(back, right, up);
    OpenVRUtils::VectorNormalize(up);
  } else {
    right[0] = 1.0f;
    right[1] = 0.0f;
    right[2] = 0.0f;
    right[3] = 0.0f;
    OpenVRUtils::VectorCrossProduct(back, right, up);
    OpenVRUtils::VectorNormalize(up);
    OpenVRUtils::VectorCrossProduct(up, back, right);
    OpenVRUtils::VectorNormalize(right);
  }

  OpenVRUtils::MatrixFastInverseGLGL(matrix, m_matrix);

  return false;
}

void OpenVRScenePicker::LaserClick(int glutButton, int glutState)
{
  m_inputHandlers->MouseFunc(glutButton, glutState, m_clickX, m_clickY, 0);
}

bool OpenVRScenePicker::IsLaserAllowed(unsigned deviceIndex) const
{
  return m_active && (m_ownerID == deviceIndex || m_ownerID == ~0u);
}

float const* OpenVRScenePicker::GetMatrix() const
{
  return m_matrix;
}

float const* OpenVRScenePicker::GetLaserColor() const
{
  static const float color[4] = {1.0f, 1.0f, 0.0f, 0.5f};
  return color;
}
