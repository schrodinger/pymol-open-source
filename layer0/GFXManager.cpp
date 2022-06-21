#include "GFXManager.h"

GFXManager::GFXManager(CShaderMgr* shaderMgr)
    : m_shaderMgr(shaderMgr)
{
}

GFXAPIBackend GFXManager::backend() const noexcept
{
  return m_apiBackend;
}
