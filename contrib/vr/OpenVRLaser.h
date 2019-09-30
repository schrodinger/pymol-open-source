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

#ifndef _H_OpenVRLaser
#define _H_OpenVRLaser

// system headers
#include "os_gl.h"

class OpenVRLaser {
public:
  static double MAX_LENGTH;
  static double MIN_WIDTH;

public:
  OpenVRLaser();

  void Init();
  void Free();

  void Show();
  void Hide();
  bool IsVisible() const;

  void SetColor(float r, float g, float b, float a = 1.0f);
  void GetColor(float* color) const;

  void SetLength(float length);
  void SetWidth(float width);

  void Draw();

private:
  void InitGeometry();
  void FreeGeometry();

  bool InitShaders();
  void FreeShaders();

private:
  bool m_valid;
  bool m_visible;
  float m_length;
  float m_width;
  float m_color[4];

  // geometry
  GLuint m_vertexArrayID;
  GLuint m_vertexBufferID;
  GLuint m_vertexCount;

  // shader
  GLuint m_programID;
  GLint m_scaleUniform;
  GLint m_colorUniform;
};

inline void OpenVRLaser::Show()
{
  m_visible = true;
}

inline void OpenVRLaser::Hide()
{
  m_visible = false;
}

inline bool OpenVRLaser::IsVisible() const {
  return m_visible;
}

#endif /* _H_OpenVRLaser */
