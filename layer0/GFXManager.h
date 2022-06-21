#pragma once

enum class GFXAPIBackend { OPENGL };

class CShaderMgr;

class GFXManager
{
public:
  GFXManager(CShaderMgr* shaderMgr);
  GFXAPIBackend backend() const noexcept;

private:
  CShaderMgr* m_shaderMgr = nullptr;
  GFXAPIBackend m_apiBackend = GFXAPIBackend::OPENGL;
};
