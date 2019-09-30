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

#ifndef _H_OpenVRControllerModel
#define _H_OpenVRControllerModel

// system headers
#include "openvr.h"

// pymol headers
#include "PyMOLGlobals.h"
#include "ShaderMgr.h"

// local headers
#include "OpenVRQuad.h"

class OpenVRControllerModel
{
public:
  OpenVRControllerModel(const std::string & sRenderModelName);
  ~OpenVRControllerModel();

  bool Init(PyMOLGlobals * G, const vr::RenderModel_t & vrModel, const vr::RenderModel_TextureMap_t & vrDiffuseTexture);
  void Free();

  void Draw();

  const std::string & GetName() const { return m_sModelName; }

private:
  void InitGeometry(const vr::RenderModel_t &vrModel);
  void FreeGeometry();

  void InitTexture(const vr::RenderModel_TextureMap_t &vrDiffuseTexture);
  void FreeTexture();

  bool InitShaders(PyMOLGlobals * G);
  void FreeShaders();

private:
  GLuint m_glVertBuffer;
  GLuint m_glIndexBuffer;
  GLuint m_glVertArray;
  GLuint m_glTexture;
  GLsizei m_unVertexCount;
  std::string m_sModelName;
  CShaderPrg *m_pShader;
};

void ShutdownRenderModels();
OpenVRControllerModel *FindOrLoadRenderModel(PyMOLGlobals *G, const char *pchRenderModelName);

#endif /* _H_OpenVRController */
