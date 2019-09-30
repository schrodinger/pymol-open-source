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

// this header
#include "OpenVRController.h"

// system headers
#include "os_gl.h"

OpenVRController::OpenVRController()
  : m_init(false)
  , m_deviceIndex(~0U)
  , m_bShowController(true)
  , m_gripIsPressed(false)
  , m_pRenderModel(NULL)
{
}

void OpenVRController::Init()
{
  m_hintsQuad = new OpenVRQuad();
  m_hintsQuad->SetSize(0.04f, 0.04f);
  m_hintsQuad->SetAlpha(0.25f);
  m_laser.Init();
  m_init = true;
}

void OpenVRController::Free()
{
  m_laser.Free();
  delete m_hintsQuad;
  m_hintsQuad = 0;
  m_init = false;
}

bool OpenVRController::IsInitialized()
{
  return m_init;
}

void OpenVRController::SetHintsTexture(GLuint hintsTexture, unsigned spriteCount)
{
  m_hintsQuad->SetTexture(hintsTexture, spriteCount);
  m_hintsQuad->SetMirror(true);
}

void OpenVRController::SetHintsIndex(unsigned index)
{
  m_hintsQuad->SetSprite(index);
}

void OpenVRController::Show(bool isVisible)
{
  m_bShowController = isVisible;
}

bool OpenVRController::IsVisible() const
{
  return m_bShowController;
}

void OpenVRController::LaserShow(bool isVisible)
{
  if (isVisible)
    m_laser.Show();
  else
    m_laser.Hide();
}

bool OpenVRController::IsLaserVisible() const
{
  return m_laser.IsVisible();
}

void OpenVRController::Draw()
{
  GL_DEBUG_FUN();

  glPushMatrix();
  glMultMatrixf(m_pose);

  if (m_pRenderModel) {
    m_pRenderModel->Draw();

    glPushMatrix();
    glTranslatef(0.0f, 0.005f, 0.049f);
    glRotatef(96.8f, 1.0f, 0.0f, 0.0f);
    m_hintsQuad->Draw();
    glPopMatrix();
  }

  m_laser.Draw();

  glPopMatrix();
}

bool OpenVRController::GetLaserRay(float* origin, float* dir) const
{
  if (IsLaserVisible()) {
    if (origin) {
      origin[0] = m_pose[12];
      origin[1] = m_pose[13];
      origin[2] = m_pose[14];
    }
    if (dir) {
      dir[0] = -m_pose[8];
      dir[1] = -m_pose[9];
      dir[2] = -m_pose[10];
    }
    return true;
  }
  return false;
}

unsigned OpenVRController::GetLaserDeviceIndex() const
{
  return m_deviceIndex;
}

void OpenVRController::SetLaserLength(float length)
{
  m_laser.SetLength(length);
}

void OpenVRController::SetLaserWidth(float width)
{
  m_laser.SetWidth(width);
}

void OpenVRController::SetLaserColor(float r, float g, float b, float a)
{
  m_laser.SetColor(r, g, b, a);
}

void OpenVRController::SetLaserColor(float const color[])
{
  m_laser.SetColor(color[0], color[1], color[2], color[3]);
}
