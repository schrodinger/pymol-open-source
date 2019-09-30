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

#ifndef _H_OpenVRUtils
#define _H_OpenVRUtils

// system headers
#include "os_gl.h"

namespace OpenVRUtils {

GLuint CompileProgram(char const* vertexShader,  char const* fragmentShader, char const* attributes[] = 0);

GLuint LoadTexture(unsigned width, unsigned height, unsigned char const* ptr);
GLuint LoadTexture(char const* filename);

void VectorNormalize(float v[]);
void VectorCrossProduct(float const v1[], float const v2[], float cross[]);

void MatrixFastInverseGLGL(float const srcGL44[], float dstGL44[]);
void MatrixFastInverseVRGL(float const srcVR34[], float dstGL44[]);
void MatrixCopyVRGL(float const srcVR34[], float dstGL44[]);

} // namespace OpenVRUtils

#endif // _H_OpenVRUtils
