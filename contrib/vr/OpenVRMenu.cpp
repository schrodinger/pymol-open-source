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
#include "OpenVRMenu.h"

// local headers
#include "OpenVRUtils.h"

OpenVRMenu::OpenVRMenu()
  : m_width(640)
  , m_height(480)
  , m_visibleX(0)
  , m_visibleY(0)
  , m_visibleWidth(640)
  , m_visibleHeight(480)
  , m_alpha(1.0f)
  , m_sceneColor(0.2f)
  , m_sceneAlpha(0.75f)
  , m_backdropColor(0.2f)
  , m_backdropAlpha(0.75f)
  , m_distance(3.0f)
  , m_fovTangent(0.75f)
  , m_valid(false)
  , m_visible(false)
  , m_backdrop(false)
  , m_overlay(false)
  , m_worldHalfWidth(2.0f)
  , m_worldHalfHeight(1.5f)
  , m_frameBufferID(0)
  , m_guiTextureID(0)
  , m_vertexArrayID(0)
  , m_vertexBufferID(0)
  , m_vertexCount(0)
  , m_programID(0)
  , m_guiTransformUniform(-1)
  , m_hotspotTransformUniform(-1)
  , m_hotspotColorUniform(-1)
  , m_factorsUniform(-1)
  , m_guiTextureUniform(-1)
  , m_sceneTextureUniform(-1)
  , m_backdropProgramID(0)
  , m_backdropColorUniform(-1)
{
  memset(m_matrix, 0, sizeof(m_matrix));
  m_matrix[0] = m_matrix[5] = m_matrix[10] = m_matrix[15] = 1.0f;

  HideHotspot();
}

void OpenVRMenu::Init(OpenVRInputHandlers* inputHandlers)
{
  m_inputHandlers = inputHandlers;
  InitGeometry();
  m_valid = InitShaders();
}

void OpenVRMenu::Free()
{
  FreeBuffers();
  FreeShaders();
  FreeGeometry();
}

void OpenVRMenu::InitGeometry()
{
  static struct Vertex_t {
    Vector3f position;
    Vector3f texcoord;
  } vertices[] = {
    {{-1.0f, +1.0f}, {0.0f, 1.0f}},
    {{+1.0f, +1.0f}, {1.0f, 1.0f}},
    {{-1.0f, -1.0f}, {0.0f, 0.0f}},
    {{+1.0f, -1.0f}, {1.0f, 0.0f}},
  };

  m_vertexCount = sizeof(vertices) / sizeof(vertices[0]);

  // Create and bind a VAO to hold state for this model
  glGenVertexArrays(1, &m_vertexArrayID);
  glBindVertexArray(m_vertexArrayID);

  // Populate a vertex buffer
  glGenBuffers(1, &m_vertexBufferID);
  glBindBuffer(GL_ARRAY_BUFFER, m_vertexBufferID);
  glBufferData(GL_ARRAY_BUFFER, sizeof(vertices), vertices, GL_STATIC_DRAW);

  // Identify the components in the vertex buffer
  glEnableVertexAttribArray(0);
  glVertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE, sizeof(Vertex_t), (void *)offsetof(Vertex_t, position));
  glEnableVertexAttribArray(1);
  glVertexAttribPointer(1, 2, GL_FLOAT, GL_FALSE, sizeof(Vertex_t), (void *)offsetof(Vertex_t, texcoord));

  glBindBuffer(GL_ARRAY_BUFFER, 0);
  glBindVertexArray(0);
}

void OpenVRMenu::FreeGeometry()
{
  if (m_vertexBufferID) {
    glDeleteVertexArrays(1, &m_vertexArrayID);
    glDeleteBuffers(1, &m_vertexBufferID);
    m_vertexArrayID = 0;
    m_vertexBufferID = 0;
  }
}

bool OpenVRMenu::InitShaders()
{
  const char* vs =
    "attribute vec2 position;\n"
    "attribute vec2 texcoords_in;\n\n"
    "varying vec2 texcoords;\n\n"
    "void main() {\n"
      "texcoords = texcoords_in;\n"
      "gl_Position = gl_ModelViewProjectionMatrix * vec4(position, 0.0, 1.0);\n"
    "}\n";
  const char* fs =
    "uniform vec4 guiTransform;\n"
    "uniform sampler2D guiTexture;\n"
    "uniform sampler2D sceneTexture;\n"
    "uniform vec4 hotspotTransform;\n"
    "uniform vec4 hotspotColor;\n"
    "uniform vec4 factors;\n"
    "varying vec2 texcoords;\n"
    "void main() {\n"
      "vec2 delta = hotspotTransform.xy * texcoords + hotspotTransform.zw;\n"
      "float spotFactor = step(dot(delta, delta), 1.0);\n"
      "vec2 guiTexcoords = guiTransform.xy * texcoords + guiTransform.zw;\n"
      "vec4 guiColor = texture2D(guiTexture, guiTexcoords);\n"
      "vec4 sceneColor = texture2D(sceneTexture, guiTexcoords);\n"
      "sceneColor = vec4(mix(sceneColor.rgb, vec3(1.0, 1.0, 1.0), factors.x), factors.y);\n"
      "vec4 color = vec4(mix(sceneColor.rgb, guiColor.rgb, guiColor.a), mix(sceneColor.a, 1.0, guiColor.a) * factors.a);\n"
      "gl_FragColor = color + spotFactor * hotspotColor.a * vec4(hotspotColor.rgb, factors.a);\n"
    "}\n";
  const char* fsBackdrop =
    "uniform vec4 backdropColor;\n"
    "void main() {\n"
      "gl_FragColor = backdropColor;\n"
    "}\n";

  m_programID = OpenVRUtils::CompileProgram(vs, fs);
  m_guiTransformUniform = glGetUniformLocation(m_programID, "guiTransform");
  m_hotspotTransformUniform = glGetUniformLocation(m_programID, "hotspotTransform");
  m_hotspotColorUniform = glGetUniformLocation(m_programID, "hotspotColor");
  m_factorsUniform = glGetUniformLocation(m_programID, "factors");
  m_guiTextureUniform = glGetUniformLocation(m_programID, "guiTexture");
  m_sceneTextureUniform = glGetUniformLocation(m_programID, "sceneTexture");

  m_backdropProgramID = OpenVRUtils::CompileProgram(vs, fsBackdrop);
  m_backdropColorUniform = glGetUniformLocation(m_backdropProgramID, "backdropColor");

  return m_programID && m_hotspotTransformUniform != -1 && m_hotspotColorUniform != -1;
}

void OpenVRMenu::FreeShaders()
{
  glDeleteProgram(m_backdropProgramID);
  glDeleteProgram(m_programID);
}

void OpenVRMenu::InitBuffers(unsigned width, unsigned height)
{
  if (m_frameBufferID)
    FreeBuffers();

  m_visibleX = m_visibleY = 0;
  m_width = m_visibleWidth = width;
  m_height = m_visibleHeight = height;

  // framebuffer
  glGenFramebuffersEXT(1, &m_frameBufferID);
  glBindFramebufferEXT(GL_FRAMEBUFFER, m_frameBufferID);

  // - color
  glGenTextures(1, &m_guiTextureID);
  glBindTexture(GL_TEXTURE_2D, m_guiTextureID);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAX_LEVEL, 0);
  glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA8, m_width, m_height, 0, GL_RGBA, GL_UNSIGNED_BYTE, NULL);
  glFramebufferTexture2DEXT(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, m_guiTextureID, 0);

  glBindTexture(GL_TEXTURE_2D, 0);
  glBindFramebufferEXT(GL_FRAMEBUFFER, 0);
}

void OpenVRMenu::FreeBuffers()
{
  glDeleteTextures(1, &m_guiTextureID);
  glDeleteFramebuffers(1, &m_frameBufferID);
}

void OpenVRMenu::Crop(unsigned x, unsigned y, unsigned width, unsigned height)
{
  m_visibleX = x;
  m_visibleY = y;
  m_visibleWidth = width;
  m_visibleHeight = height;
}

void OpenVRMenu::Start(unsigned width, unsigned height, bool clear)
{
  if (m_width != width || m_height != height || !m_frameBufferID)
    InitBuffers(width, height);

  glBindFramebufferEXT(GL_FRAMEBUFFER, m_frameBufferID);

  if (clear) {
    glClearColor(0.0f, 0.0f, 0.0f, 0.0f);
    glClear(GL_DEPTH_BUFFER_BIT | GL_COLOR_BUFFER_BIT);
    glClearColor(0.0f, 0.0f, 0.0f, 1.0f);
  }
}

void OpenVRMenu::Finish()
{
  glBindFramebufferEXT(GL_FRAMEBUFFER, 0);
}

void OpenVRMenu::Show(GLfloat const* headMatrix, unsigned ownerID)
{
  memcpy(m_matrix, headMatrix, sizeof(m_matrix));

  // constrain matrix to a zero roll angle
  float* right = m_matrix;
  float* up = m_matrix + 4;
  float* back = m_matrix + 8;
  if (fabsf(back[1]) < 0.95f) {
    up[0] = 0.0f;
    up[1] = 1.0f;
    up[2] = 0.0f;
    OpenVRUtils::VectorCrossProduct(up, back, right);
    OpenVRUtils::VectorNormalize(right);
    OpenVRUtils::VectorCrossProduct(back, right, up);
    OpenVRUtils::VectorNormalize(up);
  }

  m_visible = true;
  m_ownerID = ownerID;
}

void OpenVRMenu::Hide()
{
  m_visible = false;
  m_ownerID = ~0U;
}

void OpenVRMenu::ShowtHotspot(int x, int y, float const color[] /* = 0 */)
{
  m_hotspot.x = x;
  m_hotspot.y = y;
  if (color) {
    m_hotspot.color[0] = color[0];
    m_hotspot.color[1] = color[1];
    m_hotspot.color[2] = color[2];
    m_hotspot.color[3] = color[3];
  }
}

void OpenVRMenu::HideHotspot()
{
  m_hotspot.x = -2 * m_hotspot.radius;
  m_hotspot.y = -2 * m_hotspot.radius;
}

void OpenVRMenu::DrawBackdrop()
{
  if (!m_backdropProgramID)
    return;

  glMatrixMode(GL_PROJECTION);
  glPushMatrix();
  glLoadIdentity();
  glMatrixMode(GL_MODELVIEW);
  glPushMatrix();
  glLoadIdentity();

  glDepthFunc(GL_ALWAYS);
  glDepthMask(GL_FALSE);
  glUseProgram(m_backdropProgramID);
  glBindVertexArray(m_vertexArrayID);

  if (m_backdropColorUniform != -1) {
    float backdropColor[4] = {m_backdropColor, m_backdropColor, m_backdropColor, m_backdropAlpha};
    glUniform4fv(m_backdropColorUniform, 1, backdropColor);
  }

  glDrawArrays(GL_TRIANGLE_STRIP, 0, m_vertexCount);

  glBindVertexArray(0);
  glUseProgram(0);
  glDepthMask(GL_TRUE);
  glDepthFunc(GL_LESS);

  glMatrixMode(GL_PROJECTION);
  glPopMatrix();
  glMatrixMode(GL_MODELVIEW);
  glPopMatrix();
}

void OpenVRMenu::Draw(GLuint sceneTextureID /* = 0 */)
{
  if (!m_valid || !m_visible)
    return;

  if (m_backdrop)
    DrawBackdrop();

  // recalculate extents of the panel
  m_worldHalfWidth = m_distance * m_fovTangent;
  m_worldHalfHeight = m_worldHalfWidth * m_height / m_width;

  // take copping into account if requested
  float m_visibleWidthScale = (float)m_visibleWidth / m_width;
  float m_visibleHeightScale = (float)m_visibleHeight / m_height;
  m_worldHalfWidth *= m_visibleWidthScale;
  m_worldHalfHeight *= m_visibleHeightScale;

  glPushMatrix();
  glMultMatrixf(m_matrix);
  glTranslatef(0.0f, 0.0f, -m_distance);
  glScalef(m_worldHalfWidth, m_worldHalfHeight, 1.0f);

  if (m_overlay)
    glDepthFunc(GL_ALWAYS);
  glUseProgram(m_programID);
  glBindVertexArray(m_vertexArrayID);
  glActiveTexture(GL_TEXTURE0);
  glBindTexture(GL_TEXTURE_2D, m_guiTextureID);
  glActiveTexture(GL_TEXTURE1);
  glBindTexture(GL_TEXTURE_2D, sceneTextureID);

  float guiTransform[4] = {
    m_visibleWidthScale, m_visibleHeightScale,
    1.0f - m_visibleWidthScale, 1.0f - m_visibleHeightScale,
  };
  float hotspotTransform[4] = {
    (float)m_visibleWidth / m_hotspot.radius,
    (float)m_visibleHeight / m_hotspot.radius,
    -(float)(m_hotspot.x - m_visibleX) / m_hotspot.radius,
    -(float)(m_hotspot.y - m_visibleY) / m_hotspot.radius,
  };
  float factors[4] = {m_sceneColor, m_sceneAlpha, 0.0f, m_alpha};

  glUniform4fv(m_guiTransformUniform, 1, guiTransform);
  glUniform4fv(m_hotspotTransformUniform, 1, hotspotTransform);
  glUniform4fv(m_hotspotColorUniform, 1, m_hotspot.color);
  glUniform4fv(m_factorsUniform, 1, factors);
  glUniform1i(m_guiTextureUniform, 0);
  glUniform1i(m_sceneTextureUniform, 1);

  glDrawArrays(GL_TRIANGLE_STRIP, 0, m_vertexCount);

  glActiveTexture(GL_TEXTURE1);
  glBindTexture(GL_TEXTURE_2D, 0);
  glActiveTexture(GL_TEXTURE0);
  glBindTexture(GL_TEXTURE_2D, 0);
  glBindVertexArray(0);
  glUseProgram(0);
  if (m_overlay)
    glDepthFunc(GL_LESS);

  glPopMatrix();
}

bool OpenVRMenu::IntersectRay(GLfloat const* origin, GLfloat const* dir, int* x, int* y, float* out_distance /* = 0 */)
{
  GLfloat const* right = &m_matrix[0];  // head coordinates X-axis
  GLfloat const* up = &m_matrix[4];     // head coordinates Y-axis
  GLfloat const* normal = &m_matrix[8]; // head coordinates Z-axis
  GLfloat const* pivot = &m_matrix[12]; // head coordinates origin

  // check ray direction
  float DdotN = dir[0] * normal[0] + dir[1] * normal[1] + dir[2] * normal[2];
  static const float COS100 = -0.173648f;
  if (DdotN > COS100)
    return false; // no chance/need to intersect

  // check ray distance
  float OtoP[3] = {pivot[0] - origin[0], pivot[1] - origin[1], pivot[2] - origin[2]};
  float OtoPdotN = OtoP[0] * normal[0] + OtoP[1] * normal[1] + OtoP[2] * normal[2];
  float distance = (OtoPdotN - m_distance) / DdotN;
  if (distance < 0.0f)
    return false;

  // get raw intersection
  float DdotX = dir[0] * right[0] + dir[1] * right[1] + dir[2] * right[2];
  float OtoPdotX = OtoP[0] * right[0] + OtoP[1] * right[1] + OtoP[2] * right[2];
  float rawX = distance * DdotX - OtoPdotX;
  float DdotY = dir[0] * up[0] + dir[1] * up[1] + dir[2] * up[2];
  float OtoPdotY = OtoP[0] * up[0] + OtoP[1] * up[1] + OtoP[2] * up[2];
  float rawY = distance * DdotY - OtoPdotY;

  // check out of bounds
  int pixelX = (int)((rawX + m_worldHalfWidth) * 0.5f * m_visibleWidth / m_worldHalfWidth + 0.5f);
  int pixelY = (int)((rawY + m_worldHalfHeight) * 0.5f * m_visibleHeight / m_worldHalfHeight + 0.5f);
  if (pixelX < 0 || pixelY < 0 || pixelX >= m_visibleWidth || pixelY >= m_visibleHeight)
    return false;

  *x = pixelX + m_visibleX;
  *y = pixelY + m_visibleY;
  if (out_distance)
    *out_distance = distance;
  return true;
}

bool OpenVRMenu::LaserShoot(float const* origin, float const* dir, float const* color, float* distance /* = 0 */)
{
  int x, y;
  bool hit = IntersectRay(origin, dir, &x, &y, distance);
  if (hit) {
    ShowtHotspot(x, y, color);
    m_inputHandlers->MotionFunc(x, y, 0);
  } else {
    HideHotspot();
  }
  return hit;
}

void OpenVRMenu::LaserClick(int glutButton, int glutState)
{
  if (m_hotspot.x >= m_visibleX && m_hotspot.y >= m_visibleY && m_hotspot.x < m_visibleX + m_visibleWidth && m_hotspot.y < m_visibleY + m_visibleHeight) {
    m_inputHandlers->MouseFunc(glutButton, glutState, m_hotspot.x, m_hotspot.y, 0);
  }
}

bool OpenVRMenu::IsLaserAllowed(unsigned deviceIndex) const
{
  return m_visible && (m_ownerID == deviceIndex || m_ownerID == ~0u);
}

float const* OpenVRMenu::GetLaserColor() const
{
  static const float color[4] = {0.0f, 1.0f, 1.0f, 1.0f};
  return color;
}
