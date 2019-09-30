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

#ifndef _H_OpenVRMenu
#define _H_OpenVRMenu

// system headers
#include "openvr.h"
#include "os_gl.h"

// pymol headers
#include "PyMOLGlobals.h"
#include "ShaderMgr.h"

// local headers
#include "OpenVRMode.h"
#include "OpenVRLaserTarget.h"

class OpenVRMenu : public OpenVRLaserTarget {
  struct Hotspot_t {
    int x, y, radius;
    float color[4];

    Hotspot_t() : x(0), y(0), radius(2) {
      color[0] = 1.0f;
      color[1] = 1.0f;
      color[2] = 1.0f;
      color[3] = 0.25f;
    }
  };

public:
  OpenVRMenu();

  void Init(OpenVRInputHandlers* inputHandlers);
  void Free();

  void SetSize(float distance, float fovTangent);
  void SetAlpha(float alpha);
  void SetSceneColor(float grayLevel, float alpha);
  void SetBackdropColor(float grayLevel, float alpha);
  void SetBackdrop(bool enable);
  void SetOverlay(bool enable);

  void Crop(unsigned x, unsigned y, unsigned width, unsigned height);
  void Start(unsigned width, unsigned height, bool clear);
  void Finish();

  void Show(GLfloat const* headMatrix, unsigned ownerID);
  void Hide();
  bool IsVisible() const;
  unsigned GetOwnerID() const;

  void ShowtHotspot(int x, int y, float const color[] = 0);
  void HideHotspot();

  void Draw(GLuint sceneTextureID = 0);

  bool LaserShoot(float const* origin, float const* dir, float const* color, float* distance = 0);
  void LaserClick(int glutButton, int glutState);
  bool IsLaserAllowed(unsigned deviceIndex) const;
  float const* GetLaserColor() const;

private:
  void InitGeometry();
  void FreeGeometry();

  bool InitShaders();
  void FreeShaders();

  void InitBuffers(unsigned width, unsigned height);
  void FreeBuffers();

  void DrawBackdrop();

  bool IntersectRay(GLfloat const* origin, GLfloat const* dir, int* x, int* y, float* distance = 0);

private:
  OpenVRInputHandlers* m_inputHandlers;
  unsigned m_width;
  unsigned m_height;
  unsigned m_visibleX;
  unsigned m_visibleY;
  unsigned m_visibleWidth;
  unsigned m_visibleHeight;
  float m_alpha;
  float m_sceneColor;
  float m_sceneAlpha;
  float m_backdropColor;
  float m_backdropAlpha;
  float m_distance;
  float m_fovTangent;
  bool m_valid;
  bool m_visible;
  bool m_backdrop;
  bool m_overlay;
  unsigned m_ownerID;

  GLfloat m_matrix[16];
  GLfloat m_worldHalfWidth;
  GLfloat m_worldHalfHeight;

  Hotspot_t m_hotspot;

  // offscreen framebuffer
  GLuint m_frameBufferID;
  GLuint m_guiTextureID;

  // geometry
  GLuint m_vertexArrayID;
  GLuint m_vertexBufferID;
  GLuint m_vertexCount;

  // main shader
  GLuint m_programID;
  GLint m_guiTransformUniform;
  GLint m_hotspotTransformUniform;
  GLint m_hotspotColorUniform;
  GLint m_factorsUniform;
  GLint m_guiTextureUniform;
  GLint m_sceneTextureUniform;

  // backdrop shader
  GLuint m_backdropProgramID;
  GLint m_backdropColorUniform;
};

inline bool OpenVRMenu::IsVisible() const {
  return m_visible;
}

inline unsigned OpenVRMenu::GetOwnerID() const {
  return m_ownerID;
}

inline void OpenVRMenu::SetSize(float distance, float fovTangent) {
  m_distance = distance;
  m_fovTangent = fovTangent;
}

inline void OpenVRMenu::SetAlpha(float alpha) {
  m_alpha = alpha;
}

inline void OpenVRMenu::SetSceneColor(float grayLevel, float alpha) {
  m_sceneColor = grayLevel;
  m_sceneAlpha = alpha;
}

inline void OpenVRMenu::SetBackdropColor(float grayLevel, float alpha) {
  m_backdropColor = grayLevel;
  m_backdropAlpha = alpha;
}

inline void OpenVRMenu::SetBackdrop(bool enable) {
  m_backdrop = enable;
}

inline void OpenVRMenu::SetOverlay(bool enable) {
  m_overlay = enable;
}

#endif /* _H_OpenVRMenu */
