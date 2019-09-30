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
#include "OpenVRUtils.h"

// system headers
#include <stdio.h>
#include <math.h>
#include <string>

// pymol headers
#include "MyPNG.h"

namespace OpenVRUtils {

static GLuint CompileShader(GLenum shaderType, char const* shader)
{
  GLuint shaderID = glCreateShader(shaderType);
  glShaderSource(shaderID, 1, &shader, NULL);
  glCompileShader(shaderID);

  GLint success = GL_FALSE;
  glGetShaderiv(shaderID, GL_COMPILE_STATUS, &success);
  if (!success) {
    int infoLogLength = 0;
    glGetShaderiv(shaderID, GL_INFO_LOG_LENGTH, &infoLogLength);
    if (infoLogLength > 0) {
      char *infoLog = new char[infoLogLength];
      glGetShaderInfoLog(shaderID, infoLogLength, 0, infoLog);
      printf("GLSL shader compilation failed, log follows:\n"
        "=============================================================================\n"
        "%s"
        "=============================================================================\n", infoLog);
      delete[] infoLog;
    }
    glDeleteShader(shaderID);
    shaderID = 0;
  }
  return shaderID;
}

GLuint CompileProgram(char const* vertexShader,  char const* fragmentShader, char const* attributes[] /* = 0 */)
{
  GLuint vertexShaderID = CompileShader(GL_VERTEX_SHADER, vertexShader);
  GLuint fragmentShaderID = CompileShader(GL_FRAGMENT_SHADER, fragmentShader);

  GLuint programID = 0;
  if (vertexShaderID && fragmentShaderID) {
    programID = glCreateProgram();
    glAttachShader(programID, vertexShaderID);
    glAttachShader(programID, fragmentShaderID);

    if (attributes)
      for (int i = 0; attributes[i] != 0; ++i)
        glBindAttribLocation(programID, i, attributes[i]);

    glLinkProgram(programID);

    GLint success = GL_FALSE;
    glGetProgramiv(programID, GL_LINK_STATUS, &success);
    if (!success) {
      int infoLogLength = 0;
      glGetProgramiv(programID, GL_INFO_LOG_LENGTH, &infoLogLength);
      if (infoLogLength > 0) {
        char *infoLog = new char[infoLogLength];
        glGetProgramInfoLog(programID, infoLogLength, 0, infoLog);
        printf("GLSL program linking failed, log follows:\n"
          "=============================================================================\n"
          "%s"
          "=============================================================================\n", infoLog);
        delete[] infoLog;
      }

      glDeleteProgram(programID);
      programID = 0;
    }
  }

  if (vertexShaderID) {
    glDeleteShader(vertexShaderID);
  }

  if (fragmentShaderID) {
    glDeleteShader(fragmentShaderID);
  }

  if (programID) {
    glUseProgram(programID);
    glUseProgram(0);
  }

  return programID;
}

GLuint LoadTexture(unsigned width, unsigned height, unsigned char const* ptr)
{
  GLuint texture = 0;

  glGenTextures(1, &texture);
  glBindTexture(GL_TEXTURE_2D, texture);

  glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, width, height, 0, GL_RGBA, GL_UNSIGNED_BYTE, ptr);
  glGenerateMipmap(GL_TEXTURE_2D);

  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR);

  GLfloat maxAniso;
  glGetFloatv(GL_MAX_TEXTURE_MAX_ANISOTROPY_EXT, &maxAniso);
  glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MAX_ANISOTROPY_EXT, maxAniso);

  glBindTexture(GL_TEXTURE_2D, 0);
  return texture;
}

GLuint LoadTexture(char const* filename)
{
  auto filePath = std::string(getenv("PYMOL_DATA"))
                      .append(PATH_SEP "openvr" PATH_SEP)
                      .append(filename);
  auto image = MyPNGRead(filePath.c_str());
  if (image) {
    return LoadTexture(image->getWidth(), image->getHeight(), image->bits());
  }
  return 0;
}

void VectorNormalize(float v[])
{
  double len = sqrt((v[0] * v[0]) + (v[1] * v[1]) + (v[2] * v[2]));
  if (len > 1.0e-5f) {
    float a = 1.0f / len;
    v[0] *= a;
    v[1] *= a;
    v[2] *= a;
  }
}

void VectorCrossProduct(float const v1[], float const v2[], float cross[])
{
  cross[0] = (v1[1] * v2[2]) - (v1[2] * v2[1]);
  cross[1] = (v1[2] * v2[0]) - (v1[0] * v2[2]);
  cross[2] = (v1[0] * v2[1]) - (v1[1] * v2[0]);
}

void MatrixFastInverseGLGL(float const srcGL44[], float dstGL44[])
{
  float const (*src)[4] = (float const (*)[4])srcGL44;
  float (*dst)[4] = (float (*)[4])dstGL44;

  // transpose rotation
  dst[0][0] = src[0][0];
  dst[0][1] = src[1][0];
  dst[0][2] = src[2][0];
  dst[0][3] = 0.0f;
  dst[1][0] = src[0][1];
  dst[1][1] = src[1][1];
  dst[1][2] = src[2][1];
  dst[1][3] = 0.0f;
  dst[2][0] = src[0][2];
  dst[2][1] = src[1][2];
  dst[2][2] = src[2][2];
  dst[2][3] = 0.0f;
  
  // transpose-rotated negative translation
  dst[3][0] = -(src[0][0] * src[3][0] + src[0][1] * src[3][1] + src[0][2] * src[3][2]);
  dst[3][1] = -(src[1][0] * src[3][0] + src[1][1] * src[3][1] + src[1][2] * src[3][2]);
  dst[3][2] = -(src[2][0] * src[3][0] + src[2][1] * src[3][1] + src[2][2] * src[3][2]);
  dst[3][3] = 1.0f;
}

void MatrixFastInverseVRGL(float const srcVR34[], float dstGL44[])
{
    float const (*src)[4] = (float const (*)[4])srcVR34;
    float (*dst)[4] = (float (*)[4])dstGL44;

    // transpose rotation
    dst[0][0] = src[0][0];
    dst[0][1] = src[0][1];
    dst[0][2] = src[0][2];
    dst[0][3] = 0.0f;
    dst[1][0] = src[1][0];
    dst[1][1] = src[1][1];
    dst[1][2] = src[1][2];
    dst[1][3] = 0.0f;
    dst[2][0] = src[2][0];
    dst[2][1] = src[2][1];
    dst[2][2] = src[2][2];
    dst[2][3] = 0.0f;
    
    // transpose-rotated negative translation
    dst[3][0] = -(src[0][0] * src[0][3] + src[1][0] * src[1][3] + src[2][0] * src[2][3]);
    dst[3][1] = -(src[0][1] * src[0][3] + src[1][1] * src[1][3] + src[2][1] * src[2][3]);
    dst[3][2] = -(src[0][2] * src[0][3] + src[1][2] * src[1][3] + src[2][2] * src[2][3]);
    dst[3][3] = 1.0f;
}

void MatrixCopyVRGL(float const srcVR34[], float dstGL44[])
{
  float const (*src)[4] = (float const (*)[4])srcVR34;
  float (*dst)[4] = (float(*)[4])dstGL44;
  dst[0][0] = src[0][0]; dst[0][1] = src[1][0]; dst[0][2] = src[2][0]; dst[0][3] = 0.0f;
  dst[1][0] = src[0][1]; dst[1][1] = src[1][1]; dst[1][2] = src[2][1]; dst[1][3] = 0.0f;
  dst[2][0] = src[0][2]; dst[2][1] = src[1][2]; dst[2][2] = src[2][2]; dst[2][3] = 0.0f;
  dst[3][0] = src[0][3]; dst[3][1] = src[1][3]; dst[3][2] = src[2][3]; dst[3][3] = 1.0f;
}

} // namespace OpenVRUtils
