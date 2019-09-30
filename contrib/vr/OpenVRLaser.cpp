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
#include "OpenVRLaser.h"

// pymol header
#include "Setting.h"

// local headers
#include "OpenVRUtils.h"

double OpenVRLaser::MAX_LENGTH = 100.0f;
double OpenVRLaser::MIN_WIDTH = 1.0f;

OpenVRLaser::OpenVRLaser()
  : m_valid(false)
  , m_visible(false)
  , m_length(OpenVRLaser::MAX_LENGTH)
  , m_width(OpenVRLaser::MIN_WIDTH)
  , m_vertexArrayID(0)
  , m_vertexBufferID(0)
  , m_vertexCount(0)
  , m_programID(0)
  , m_scaleUniform(-1)
  , m_colorUniform(-1)
{
  SetColor(0.0f, 1.0f, 1.0f, 0.25f);
}

void OpenVRLaser::Init()
{
  InitGeometry();
  m_valid = InitShaders();
}

void OpenVRLaser::Free()
{
  FreeShaders();
  FreeGeometry();
}

void OpenVRLaser::InitGeometry()
{
  static struct Vertex_t {
    float position[3];
  } vertices[] = {
    {{0.0f, 0.0f, 0.0f}},
    {{0.0f, 0.0f, -1.0f}},
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
  glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(Vertex_t), (void *)offsetof(Vertex_t, position));

  glBindBuffer(GL_ARRAY_BUFFER, 0);
  glBindVertexArray(0);
}

void OpenVRLaser::FreeGeometry()
{
  if (m_vertexBufferID) {
    glDeleteVertexArrays(1, &m_vertexArrayID);
    glDeleteBuffers(1, &m_vertexBufferID);
    m_vertexArrayID = 0;
    m_vertexBufferID = 0;
  }
}

bool OpenVRLaser::InitShaders()
{
  const char* vs =
    "attribute vec3 position;\n"
    "uniform vec4 scale;\n"
    "varying vec2 texcoords;\n\n"
    "void main() {\n"
      "gl_Position = gl_ModelViewProjectionMatrix * vec4(position * scale.xyz, 1.0);\n"
    "}\n";
  const char* fs =
    "uniform vec4 color;\n"
    "void main() {\n"
      "gl_FragColor = color;\n"
    "}\n";

  m_programID = OpenVRUtils::CompileProgram(vs, fs);
  m_scaleUniform = glGetUniformLocation(m_programID, "scale");
  m_colorUniform = glGetUniformLocation(m_programID, "color");
  return m_programID && m_scaleUniform != -1 && m_colorUniform != -1;
}

void OpenVRLaser::FreeShaders()
{
  glDeleteProgram(m_programID);
}

void OpenVRLaser::SetColor(float r, float g, float b, float a /* = 1.0f */)
{
  m_color[0] = r;
  m_color[1] = g;
  m_color[2] = b;
  m_color[3] = a;
}

void OpenVRLaser::GetColor(float* color) const
{
  color[0] = m_color[0];
  color[1] = m_color[1];
  color[2] = m_color[2];
  color[3] = m_color[3];
}

void OpenVRLaser::SetLength(float length)
{
  m_length = length > 0.0f ? length : OpenVRLaser::MAX_LENGTH;
}

void OpenVRLaser::SetWidth(float width)
{
  m_width = width > 1.0f ? width : OpenVRLaser::MIN_WIDTH;
}

void OpenVRLaser::Draw()
{
  if (!m_valid || !m_visible)
    return;

  glUseProgram(m_programID);
  glBindVertexArray(m_vertexArrayID);

  float scale[4] = {1.0f, 1.0f, m_length, 1.0f};

  glUniform4fv(m_scaleUniform, 1, scale);
  glUniform4fv(m_colorUniform, 1, m_color);

  glLineWidth(m_width); // set here to avoid the influence of previous rendering steps
  glDrawArrays(GL_LINES, 0, m_vertexCount);

  glBindVertexArray(0);
  glUseProgram(0);
}
