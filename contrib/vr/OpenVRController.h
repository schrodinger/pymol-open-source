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

#ifndef _H_OpenVRController
#define _H_OpenVRController

// system headers
#include "openvr.h"

// local headers
#include "OpenVRLaserSource.h"
#include "OpenVRControllerModel.h"
#include "OpenVRLaser.h"

class OpenVRController : public OpenVRLaserSource {
public:
  OpenVRController();
 
  void Init();
  void Free();
  bool IsInitialized();

  void SetHintsTexture(GLuint hintsTexture, unsigned spriteCount);
  void SetHintsIndex(unsigned index);

  void Draw();

  float *GetPose() {return m_pose;} // it's not safe =)
  float *GetWorldToControllerMatrix() {return m_worldToController;} // it's not safe =)

  void Show(bool isVisible);
  bool IsVisible() const;

  void LaserShow(bool isVisible);
  bool IsLaserVisible() const;
  bool GetLaserRay(float* origin, float* dir) const;
  unsigned GetLaserDeviceIndex() const;
  void SetLaserLength(float length);
  void SetLaserWidth(float width);
  void SetLaserColor(float r, float g, float b, float a);
  void SetLaserColor(float const color[]);

  bool isGripPressed() const { return m_gripIsPressed; }
  void pressGrip(bool press) { m_gripIsPressed = press; }

public:
  vr::TrackedDeviceIndex_t m_deviceIndex;
  OpenVRControllerModel *m_pRenderModel;
  std::string m_sRenderModelName;

private:
  bool m_init;
  bool m_bShowController;
  GLfloat m_pose[16]; // model2world matrix 
  OpenVRLaser m_laser;
  GLfloat m_worldToController[16];
  bool m_gripIsPressed;
  OpenVRQuad* m_hintsQuad;
};

#endif /* _H_OpenVRController */
