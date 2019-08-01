#include "PostProcess.h"

#include "gl_def.h"
#include "GenericBuffer.h"
#include "os_gl.h"
#include "pymol/memory.h"

#ifndef _PYMOL_NO_AA_SHADERS
#endif

void PostProcess::activateRTAsTexture(std::size_t idx, GLuint textureUnit)
{
  glActiveTexture(GL_TEXTURE0 + textureUnit);
  auto& rt = m_renderTargets[idx];
  if (rt && rt->textures().front()) {
    rt->textures().front()->bind();
  }
}

void PostProcess::activateTexture(std::size_t idx, GLuint textureUnit)
{
  glActiveTexture(GL_TEXTURE0 + textureUnit);
  auto& tex = m_textures[idx];
  if (tex) {
    tex->bind();
  }
}

void PostProcess::bindRT(std::size_t idx, bool clear)
{
  if (idx < m_renderTargets.size()) {
    auto& t = m_renderTargets[idx];
    if (t) {
      t->bind(clear);
    }
  }
}

void PostProcess::bindFBORBO(std::size_t idx)
{
  if (idx < m_renderTargets.size()) {
    auto& t = m_renderTargets[idx];
    if (t) {
      t->bindFBORBO();
    }
  }
}

const renderTarget_t::shape_type PostProcess::size(std::size_t idx) const
    noexcept
{
  return m_renderTargets[idx]->size();
}

#ifndef _PYMOL_NO_AA_SHADERS
#endif

OIT_PostProcess::OIT_PostProcess(int width, int height, renderBuffer_t* rbo)
{
  if (TM3_IS_ONEBUF) {
    auto rt0 = pymol::make_unique<renderTarget_t>(width, height);
    rt0->layout({{4, rt_layout_t::FLOAT}}, rbo);
    m_renderTargets.push_back(std::move(rt0));

    auto rt1 = pymol::make_unique<renderTarget_t>(width, height);
    rt1->layout({{1, rt_layout_t::FLOAT}}, rt0->rbo());
    m_renderTargets.push_back(std::move(rt1));
  } else {
    std::vector<rt_layout_t> layouts;
    layouts.emplace_back(4, rt_layout_t::FLOAT);
    if (GLEW_VERSION_3_0) {
      layouts.emplace_back(1, rt_layout_t::FLOAT);
    } else {
      layouts.emplace_back(2, rt_layout_t::FLOAT);
    }

    auto rt0 = pymol::make_unique<renderTarget_t>(width, height);
    rt0->layout(std::move(layouts), rbo);
    m_renderTargets.push_back(std::move(rt0));
  }
}

void OIT_PostProcess::bindRT(std::size_t idx, bool clear)
{
#if !defined(PURE_OPENGL_ES_2) || defined(_WEBGL)
  if (TM3_IS_ONEBUF) {
    auto& rt = m_renderTargets[idx - 1];
    if (rt)
      rt->fbo()->bind();
  } else {
    const GLenum bufs[] = {GL_COLOR_ATTACHMENT0_EXT, GL_COLOR_ATTACHMENT1_EXT};
    if (!m_renderTargets.empty()) {
      auto& rt = m_renderTargets.front();
      if (rt) {
        rt->fbo()->bind();
      }
    }
    glDrawBuffers(2, bufs);
  }
  glClearColor(0.f, 0.f, 0.f, 0.f);
  glClear(GL_COLOR_BUFFER_BIT);
  glDepthMask(GL_FALSE);
  glEnable(GL_DEPTH_TEST);
  glEnable(GL_BLEND);
  glBlendFuncSeparate(GL_ONE, GL_ONE, GL_ONE, GL_ONE_MINUS_SRC_ALPHA);
#endif
}

void OIT_PostProcess::activateRTAsTexture(std::size_t idx, GLuint textureUnit)
{
  glActiveTexture(GL_TEXTURE0 + textureUnit);
  if (TM3_IS_ONEBUF) {
    auto& rt = m_renderTargets[idx];
    auto& tex = rt->textures().front();
    if (tex) {
      tex->bind();
    }
  } else {
    auto& rt = m_renderTargets.front();
    if (rt) {
      rt->textures()[idx]->bind();
    }
  }
}

