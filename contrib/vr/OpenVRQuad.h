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

#ifndef _H_OpenVRQuad
#define _H_OpenVRQuad

// system headers
#include "os_gl.h"

class OpenVRQuad {

public:
  OpenVRQuad();
  ~OpenVRQuad();

  void SetTexture(GLuint textureID, unsigned spriteCount = 1);
  void SetSprite(unsigned index);
  void SetSize(float width, float height);
  void SetAlpha(float alpha);
  void SetMirror(bool mirror);

  void Draw();

private:
  void InitGeometry();
  void FreeGeometry();

  bool InitShaders();
  void FreeShaders();

private:
  float m_width;
  float m_height;
  unsigned m_spriteIndex;
  unsigned m_spriteCount;
  float m_alpha;
  bool m_mirror;

  // geometry
  GLuint m_vertexArrayID;
  GLuint m_vertexBufferID;
  GLuint m_vertexCount;

  // texture
  GLuint m_textureID;

  // shader
  GLuint m_programID;
  GLint m_positionMulAddUniform;
  GLint m_texcoordMulAddUniform;
  GLint m_colorMulAddUniform;
  GLint m_textureUniform;
};

#endif /* _H_OpenVRQuad */
