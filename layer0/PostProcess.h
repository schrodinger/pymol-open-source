#pragma once

#include <memory>
#include <vector>

#include "GenericBuffer.h"

class PostProcess
{
protected:
  std::vector<std::unique_ptr<renderTarget_t>> m_renderTargets;
  std::vector<std::unique_ptr<textureBuffer_t>> m_textures;

public:
  /**
   * Activates/Binds render-target textures.
   * @param idx Index of rendertarget to bind
   * @param textureUnit offset of texture unit to assign (0 for GL_TEXTURE0, 1
   * for GL_TEXTURE1, etc...)
   * @note indices should preferably be passed in as an enum for clarity
   * Usage: Use render target from previous render pass as an input to
   * sample from (e.g. sample from normal map created in previous pass for
   * lighting in current pass.
   */
  virtual void activateRTAsTexture(std::size_t idx, GLuint textureUnit);

  /**
   * Activates/Binds non-render-target textures.
   * @param idx Index of texture to bind
   * @param textureUnit offset of texture unit to assign (0 for GL_TEXTURE0, 1
   * for GL_TEXTURE1, etc...)
   * @note indices should preferably be passed in as an enum for clarity
   * Usage: Use textures created/imported offline to sample from (e.g. Sample
   * from static image for background)
   */
  virtual void activateTexture(std::size_t idx, GLuint textureUnit);

  /**
   * Bind render target for rendering into
   * @param idx Index of render target to draw to
   * @param clear determines if the attachment texture should be cleared before
   * rendering.
   * Usage: Use for preparing to draw to a render target's FBO's texture.
   * Note: For RBO drawing, use bindFBORBO
   */
  virtual void bindRT(std::size_t idx, bool clear = true);

  /**
   * Bind render target (with RBO) for rendering into
   * @param idx Index of render target to draw to
   * @param clear determines if the attachment texture should be cleared before
   * rendering.
   * Usage: Use for preparing to draw to a render target's FBO's render buffer.
   * Note: For texture drawing, use bindRT
   */
  virtual void bindFBORBO(std::size_t idx);

  /**
   *
   * Obtain rendering dimension
   * @param idx Index of render target to query size from
   * Note: A subclass of Postprocess shall populate their render targets in
   * its constructor else results are undefined.
   */
  const renderTarget_t::shape_type size(std::size_t idx = 0) const noexcept;

  virtual ~PostProcess() = default;
};

#ifndef _PYMOL_NO_AA_SHADERS
#endif

class OIT_PostProcess : public PostProcess
{
public:
  struct OITRT {
    enum { ACCUM = 0, REVEALAGE = 1, TOTAL = 2 };
  };
  OIT_PostProcess(int width, int height, renderBuffer_t* rbo);

  void bindRT(std::size_t idx, bool clear = true) override;

  void activateRTAsTexture(std::size_t idx, GLuint textureUnit) override;
};

