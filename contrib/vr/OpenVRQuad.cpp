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
#include "OpenVRQuad.h"

// local headers
#include "OpenVRUtils.h"

OpenVRQuad::OpenVRQuad()
  : m_width(1.0f)
  , m_height(1.0f)
  , m_spriteIndex(0)
  , m_spriteCount(1)
  , m_alpha(1.0f)
  , m_mirror(false)
  , m_vertexArrayID(0)
  , m_vertexBufferID(0)
  , m_vertexCount(0)
  , m_textureID(0)
  , m_programID(0)
  , m_positionMulAddUniform(-1)
  , m_texcoordMulAddUniform(-1)
  , m_colorMulAddUniform(-1)
  , m_textureUniform(-1)
{
  InitGeometry();
  InitShaders();
}

OpenVRQuad::~OpenVRQuad()
{
  FreeShaders();
  FreeGeometry();
}

void OpenVRQuad::InitGeometry()
{
  struct Vertex_t {
    float position[2];
    float texcoord[2];
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

void OpenVRQuad::FreeGeometry()
{
  if (m_vertexBufferID) {
    glDeleteVertexArrays(1, &m_vertexArrayID);
    glDeleteBuffers(1, &m_vertexBufferID);
    m_vertexArrayID = 0;
    m_vertexBufferID = 0;
  }
}

bool OpenVRQuad::InitShaders()
{
  const char* attributes[] = {"positionRaw", "texcoordRaw", 0};
  const char* vs =
    "uniform vec4 positionMulAdd;\n"
    "uniform vec4 texcoordMulAdd;\n"
    "attribute vec2 positionRaw;\n"
    "attribute vec2 texcoordRaw;\n\n"
    "varying vec2 texcoord;\n\n"
    "void main() {\n"
      "vec2 position;\n"
      "texcoord = texcoordMulAdd.xy * texcoordRaw + texcoordMulAdd.zw;\n"
      "position = positionMulAdd.xy * positionRaw + positionMulAdd.zw;\n"
      "gl_Position = gl_ModelViewProjectionMatrix * vec4(position, 0.0, 1.0);\n"
    "}\n";
  const char* fs =
    "uniform vec4 colorMulAdd;\n"
    "uniform sampler2D texture;\n"
    "varying vec2 texcoord;\n"
    "void main() {\n"
      "vec4 color = texture2D(texture, texcoord);\n"
      "gl_FragColor = colorMulAdd.xxxy * color + colorMulAdd.zzzw;\n"
    "}\n";

  m_programID = OpenVRUtils::CompileProgram(vs, fs, attributes);
  m_positionMulAddUniform = glGetUniformLocation(m_programID, "positionMulAdd");
  m_texcoordMulAddUniform = glGetUniformLocation(m_programID, "texcoordMulAdd");
  m_colorMulAddUniform = glGetUniformLocation(m_programID, "colorMulAdd");
  m_textureUniform = glGetUniformLocation(m_programID, "texture");
  return m_programID;
}

void OpenVRQuad::FreeShaders()
{
  glDeleteProgram(m_programID);
}

void OpenVRQuad::SetTexture(GLuint textureID, unsigned spriteCount /* = 1 */)
{
  m_textureID = textureID;
  m_spriteCount = spriteCount;
  m_spriteIndex = 0;
}

void OpenVRQuad::SetSprite(unsigned index)
{
  m_spriteIndex = index % m_spriteCount;
}

void OpenVRQuad::SetSize(float width, float height)
{
  m_width = width;
  m_height = height;
}

void OpenVRQuad::SetAlpha(float alpha)
{
  m_alpha = alpha;
}

void OpenVRQuad::SetMirror(bool mirror)
{
  m_mirror = mirror;
}

void OpenVRQuad::Draw()
{
  glUseProgram(m_programID);
  glBindVertexArray(m_vertexArrayID);
  glActiveTexture(GL_TEXTURE0);
  glBindTexture(GL_TEXTURE_2D, m_textureID);

  float positionMulAdd[4] = {0.5f * m_width, 0.5f * m_height, 0.0f, 0.0f};
  float texcoordMulAdd[4] = {1.0f / m_spriteCount, m_mirror ? -1.0f : 1.0f, (float)m_spriteIndex / m_spriteCount, m_mirror ? 1.0f : 0.0f};
  float hasTexture = m_textureID ? 1.0f : 0.0f;
  float colorMulAdd[4] = {hasTexture, hasTexture * m_alpha, 1.0f - hasTexture, (1.0f - hasTexture) * m_alpha};

  glUniform4fv(m_positionMulAddUniform, 1, positionMulAdd);
  glUniform4fv(m_texcoordMulAddUniform, 1, texcoordMulAdd);
  glUniform4fv(m_colorMulAddUniform, 1, colorMulAdd);
  glUniform1i(m_textureUniform, 0);

  glDrawArrays(GL_TRIANGLE_STRIP, 0, m_vertexCount);

  glActiveTexture(GL_TEXTURE0);
  glBindTexture(GL_TEXTURE_2D, 0);
  glBindVertexArray(0);
  glUseProgram(0);
}
