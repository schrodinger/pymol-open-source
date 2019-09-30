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
#include "OpenVRControllerModel.h"

// system headers
#include <vector>
#include "os_std.h"
#include "os_gl.h"

// pymol headers
#include "Feedback.h"
#include "Matrix.h"

// local headers
#include "OpenVRController.h"

static std::vector<OpenVRControllerModel *> s_vecRenderModels;

//-----------------------------------------------------------------------------
// Purpose: Create/destroy GL Render Models
//-----------------------------------------------------------------------------
OpenVRControllerModel::OpenVRControllerModel(const std::string & sRenderModelName) : m_sModelName( sRenderModelName ) {
  m_glIndexBuffer = 0;
  m_glVertArray = 0;
  m_glVertBuffer = 0;
  m_glTexture = 0;
  m_pShader = NULL;
}

OpenVRControllerModel::~OpenVRControllerModel() {
  Free();
}

//-----------------------------------------------------------------------------
// Purpose: Allocates and populates the GL resources for a render model
//-----------------------------------------------------------------------------
bool OpenVRControllerModel::Init(PyMOLGlobals *G, const vr::RenderModel_t &vrModel, const vr::RenderModel_TextureMap_t &vrDiffuseTexture) {  
  InitGeometry(vrModel);
  InitTexture(vrDiffuseTexture);
  InitShaders(G);
  return true;
}

//-----------------------------------------------------------------------------
// Purpose: Frees the GL resources for a render model
//-----------------------------------------------------------------------------
void OpenVRControllerModel::Free() {
  FreeShaders();
  FreeTexture();
  FreeGeometry();
}

void OpenVRControllerModel::InitGeometry(const vr::RenderModel_t &vrModel) {
 // create and bind a VAO to hold state for this model
  glGenVertexArrays( 1, &m_glVertArray );
  glBindVertexArray( m_glVertArray );

  // Populate a vertex buffer
  glGenBuffers( 1, &m_glVertBuffer );
  glBindBuffer( GL_ARRAY_BUFFER, m_glVertBuffer );
  glBufferData( GL_ARRAY_BUFFER, sizeof(vr::RenderModel_Vertex_t) *vrModel.unVertexCount, vrModel.rVertexData, GL_STATIC_DRAW);

  // Identify the components in the vertex buffer
  glEnableVertexAttribArray(0);
  glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof( vr::RenderModel_Vertex_t ), (void *)offsetof( vr::RenderModel_Vertex_t, vPosition));
  glEnableVertexAttribArray(1);
  glVertexAttribPointer(1, 2, GL_FLOAT, GL_FALSE, sizeof( vr::RenderModel_Vertex_t ), (void *)offsetof(vr::RenderModel_Vertex_t, rfTextureCoord));

  // Create and populate the index buffer
  glGenBuffers( 1, &m_glIndexBuffer );
  glBindBuffer( GL_ELEMENT_ARRAY_BUFFER, m_glIndexBuffer );
  glBufferData( GL_ELEMENT_ARRAY_BUFFER, sizeof( uint16_t ) * vrModel.unTriangleCount * 3, vrModel.rIndexData, GL_STATIC_DRAW );

  glBindVertexArray( 0 );

  m_unVertexCount = vrModel.unTriangleCount * 3;
}

void OpenVRControllerModel::FreeGeometry() {
  if( m_glVertBuffer ) {
    glDeleteBuffers(1, &m_glIndexBuffer);
    glDeleteVertexArrays( 1, &m_glVertArray );
    glDeleteBuffers(1, &m_glVertBuffer);
    m_glIndexBuffer = 0;
    m_glVertArray = 0;
    m_glVertBuffer = 0;
  }
}

void OpenVRControllerModel::InitTexture(const vr::RenderModel_TextureMap_t &vrDiffuseTexture) {
// create and populate the texture
  glGenTextures(1, &m_glTexture );
  glBindTexture( GL_TEXTURE_2D, m_glTexture );

  glTexImage2D( GL_TEXTURE_2D, 0, GL_RGBA, vrDiffuseTexture.unWidth, vrDiffuseTexture.unHeight,
    0, GL_RGBA, GL_UNSIGNED_BYTE, vrDiffuseTexture.rubTextureMapData );

  // If this renders black ask McJohn what's wrong.
  glGenerateMipmap(GL_TEXTURE_2D);

  glTexParameteri( GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE );
  glTexParameteri( GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE );
  glTexParameteri( GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR );
  glTexParameteri( GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR );

  GLfloat fLargest;
  glGetFloatv( GL_MAX_TEXTURE_MAX_ANISOTROPY_EXT, &fLargest );
  glTexParameterf( GL_TEXTURE_2D, GL_TEXTURE_MAX_ANISOTROPY_EXT, fLargest );

  glBindTexture( GL_TEXTURE_2D, 0 ); 
}

void OpenVRControllerModel::FreeTexture() {
  if (m_glTexture) {
    glDeleteTextures(1, &m_glTexture);
    m_glTexture = 0;  
  }
}

bool OpenVRControllerModel::InitShaders(PyMOLGlobals * G) {
  m_pShader = new CShaderPrg(G, "CRenderModel", "vrcontroller.vs", "vrcontroller.fs");
  if (m_pShader) {
    m_pShader->reload();
    glBindAttribLocation(m_pShader->id, 0, "position");
    glBindAttribLocation(m_pShader->id, 1, "texcoords_in");
    glLinkProgram(m_pShader->id);
  }
  return m_pShader != NULL;
}

void OpenVRControllerModel::FreeShaders() {
  if (m_pShader) {
    delete m_pShader;
    m_pShader = NULL;
  }
}

//-----------------------------------------------------------------------------
// Purpose: Draws the render model
//-----------------------------------------------------------------------------
void OpenVRControllerModel::Draw() {
  GL_DEBUG_FUN();
  m_pShader->Enable();
  glBindVertexArray( m_glVertArray );

  glActiveTexture( GL_TEXTURE0 );
  glBindTexture( GL_TEXTURE_2D, m_glTexture );

  glDrawElements( GL_TRIANGLES, m_unVertexCount, GL_UNSIGNED_SHORT, 0 );

  glBindVertexArray( 0 );
  m_pShader->Disable();
}

// Purpose:
//-----------------------------------------------------------------------------
static
void ThreadSleep( unsigned long nMilliseconds )
{
#if defined(_WIN32)
  ::Sleep( nMilliseconds );
#elif defined(POSIX)
  usleep( nMilliseconds * 1000 );
#endif
}

// Purpose: Finds a render model we've already loaded or loads a new one
//-----------------------------------------------------------------------------
OpenVRControllerModel *FindOrLoadRenderModel(PyMOLGlobals * G, const char *pchRenderModelName)
{
  OpenVRControllerModel *pRenderModel = NULL;
  for(std::vector<OpenVRControllerModel *>::iterator i = s_vecRenderModels.begin(); i != s_vecRenderModels.end(); i++) {
    if(!
#ifdef _WIN32
        stricmp
#else
        strcasecmp
#endif
        ((*i)->GetName().c_str(), pchRenderModelName)){
      pRenderModel = *i;
      break;
    }
  }

  // load the model if we didn't find one
  if(!pRenderModel) {
    vr::RenderModel_t *pModel;
    vr::EVRRenderModelError error;
    while (true) {
      error = vr::VRRenderModels()->LoadRenderModel_Async(pchRenderModelName, &pModel);
      if (error != vr::VRRenderModelError_Loading)
        break;

      ThreadSleep( 1 );
    }

    if (error != vr::VRRenderModelError_None) {
      PRINTF
        " Unable to load render model %s - %s\n",
        pchRenderModelName,
        vr::VRRenderModels()->GetRenderModelErrorNameFromEnum(error)
      ENDF(G);
      return NULL; // move on to the next tracked device
    }

    vr::RenderModel_TextureMap_t *pTexture;
    while (true) {
      error = vr::VRRenderModels()->LoadTexture_Async(pModel->diffuseTextureId, &pTexture);
      if (error != vr::VRRenderModelError_Loading)
        break;
      ThreadSleep( 1 );
    }

    if (error != vr::VRRenderModelError_None) {
      PRINTF
        " Unable to load render texture id:%d for render model %s\n",
        pModel->diffuseTextureId,
        pchRenderModelName
      ENDF(G);
      vr::VRRenderModels()->FreeRenderModel(pModel);
      return NULL; // move on to the next tracked device
    }

    pRenderModel = new OpenVRControllerModel(pchRenderModelName);
    if (!pRenderModel->Init(G, *pModel, *pTexture)) {
      PRINTF
        " Unable to create GL model from render model %s\n",
        pchRenderModelName 
      ENDF(G);
      delete pRenderModel;
      pRenderModel = NULL;
    } else {
      s_vecRenderModels.push_back(pRenderModel);
    }
    vr::VRRenderModels()->FreeRenderModel(pModel);
    vr::VRRenderModels()->FreeTexture(pTexture);
  }
  return pRenderModel;
}

void ShutdownRenderModels() {
  for(std::vector<OpenVRControllerModel *>::iterator i = s_vecRenderModels.begin(); i != s_vecRenderModels.end(); i++) {
    (*i)->Free();
    delete (*i);
  }
  s_vecRenderModels.clear();
}
