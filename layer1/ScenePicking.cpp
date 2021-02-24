

#include"Base.h"
#include"Scene.h"
#include"ScenePicking.h"
#include"SceneRender.h"
#include"ShaderMgr.h"
#include"MemoryDebug.h"
#include"PyMOL.h"
#include"P.h"
#include "Err.h"
#include "Picking.h"

#define cRange 7

int SceneDoXYPick(PyMOLGlobals * G, int x, int y, int click_side)
{
  CScene *I = G->Scene;
  int defer_builds_mode = SettingGetGlobal_i(G, cSetting_defer_builds_mode);

  if(defer_builds_mode == 5)    /* force generation of a pickable version */
    SceneUpdate(G, true);
  if(OrthoGetOverlayStatus(G) || SettingGetGlobal_i(G, cSetting_text))
    SceneRender(G, NULL, 0, 0, NULL, 0, 0, 0, 0);       /* remove overlay if present */
  SceneDontCopyNext(G);

  I->LastPicked.context.object = NULL;
  SceneRender(G, &I->LastPicked, x, y, NULL, 0, 0, click_side, 0);
  return (I->LastPicked.context.object != NULL);
  /* did we pick something? */
}

/**
 * Query the number of RED, GREEN, BLUE, and ALPHA bitplanes from OpenGL
 * @param[out] pickconv picking instance to update
 */
static void PickColorConverterSetRgbaBitsFromGL(
    PyMOLGlobals* G, PickColorConverter& pickconv)
{
  GLint rgba_bits[4] = {4, 4, 4, 0};
  int max_check_bits = 0;

#ifdef _WEBGL
  // Can't turn off antialiasing in WebPyMOL, so we need to add some check bits
  // to filter out antialiased pixels.
  max_check_bits = 2;
#endif

  if (!SettingGet<bool>(G, cSetting_pick32bit)) {
    pickconv.setRgbaBits(rgba_bits, max_check_bits);
    return;
  }

  GLint currentFrameBuffer = G->ShaderMgr->default_framebuffer_id;

  if (SettingGet<bool>(G, cSetting_use_shaders)) {
    glGetIntegerv(GL_FRAMEBUFFER_BINDING, &currentFrameBuffer);
  }

  if (currentFrameBuffer != G->ShaderMgr->default_framebuffer_id) {
    glBindFramebuffer(GL_FRAMEBUFFER, G->ShaderMgr->default_framebuffer_id);
  }

  glGetIntegerv(GL_RED_BITS, rgba_bits + 0);
  glGetIntegerv(GL_GREEN_BITS, rgba_bits + 1);
  glGetIntegerv(GL_BLUE_BITS, rgba_bits + 2);
  glGetIntegerv(GL_ALPHA_BITS, rgba_bits + 3);

  PRINTFD(G, FB_Scene)
  " %s: GL RGBA BITS: (%d, %d, %d, %d)\n", __func__, rgba_bits[0], rgba_bits[1],
      rgba_bits[2], rgba_bits[3] ENDFD;

  if (currentFrameBuffer != G->ShaderMgr->default_framebuffer_id) {
    glBindFramebuffer(GL_FRAMEBUFFER, currentFrameBuffer);
  }

  pickconv.setRgbaBits(rgba_bits, max_check_bits);
}

/**
 * Get picking indices for the rectangle (x,y,w,h)
 */
static std::vector<unsigned> SceneGetPickIndices(PyMOLGlobals* G,
    SceneUnitContext* context, int x, int y, int w, int h, GLenum gl_buffer)
{
  CScene *I = G->Scene;
  auto& pickmgr = I->pickmgr;
  const auto use_shaders = SettingGet<bool>(G, cSetting_use_shaders);

  SceneGLClearColor(0.0, 0.0, 0.0, 0.);

  if (!pickmgr.m_valid) {
    PickColorConverterSetRgbaBitsFromGL(G, pickmgr);
  }

  const auto shift_per_pass = pickmgr.getTotalBits();
  const auto pass_max = use_shaders ? SHADER_PICKING_PASSES_MAX : 99;

  std::vector<unsigned> indices(w * h);

  if(I->grid.active)
    GridGetGLViewport(G, &I->grid);

  for (int pass = 0;; ++pass) {
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    pickmgr.m_pass = pass;

    if (!pickmgr.m_valid || !use_shaders) {
      pickmgr.resetCount();
    }
    {
      int slot;
      for(slot = 0; slot <= I->grid.last_slot; slot++) {
        if(I->grid.active) {
          GridSetGLViewport(&I->grid, slot);
        }
        SceneRenderAll(
            G, context, NULL, &pickmgr, RenderPass::Antialias, true, 0.0F, &I->grid, 0, SceneRenderWhich::AllObjects);
      }
    }

#ifndef _PYMOL_NO_MAIN
    const auto debug_pick = SettingGet<int>(G, cSetting_debug_pick);
    if(debug_pick) {
      PyMOL_SwapBuffers(G->PyMOL);
      PSleep(G, 1000000 * debug_pick / 4);
      PyMOL_SwapBuffers(G->PyMOL);
    }
#endif

#ifndef PURE_OPENGL_ES_2
    if (!hasFrameBufferBinding()) {
      glReadBuffer(gl_buffer);
    }
#endif

    // assume glReadPixels does not overrun the buffer - previous version of
    // PyMOL allocated a larger buffer to account for "buggy glReadPixels"
    std::vector<unsigned char> buffer(4 * indices.size());
    PyMOLReadPixels(x, y, w, h, GL_RGBA, GL_UNSIGNED_BYTE, &buffer[0]);

    for (size_t i = 0; i < indices.size(); ++i) {
      const auto* c = &buffer[4 * i];
      indices[i] |= pickmgr.indexFromColor(c) << (shift_per_pass * pass);
    }

    if (pickmgr.count() < (1ULL << shift_per_pass * (pass + 1))) {
      break;
    }

    if (pass + 1 == pass_max) {
      PRINTFB(G, FB_Scene, FB_Warnings)
        " Scene-Warning: Maximum number of picking passes exceeded\n"
        " (%u picking colors, %u color bits)\n",
        pickmgr.count(), shift_per_pass ENDFB(G);
      break;
    }
  }

  if(I->grid.active)
    GridSetGLViewport(&I->grid, -1);

  pickmgr.m_valid = true;

  return indices;
}

/**
 * Pick a single point
 *
 * @param[out] pick Copy picking result here
 * @param x X position in the window to pick
 * @param y Y position in the window to pick
 */
static void SceneRenderPickingSinglePick(PyMOLGlobals* G,
    SceneUnitContext* context, Picking* pick, int x, int y,
    GLenum render_buffer)
{
  CScene *I = G->Scene;
  const int debug_pick = SettingGet<int>(G, cSetting_debug_pick);
  const int cRangeVal = DIP2PIXEL(cRange);
  const int h = (cRangeVal * 2 + 1), w = (cRangeVal * 2 + 1);

  auto indices = SceneGetPickIndices(
      G, context, x - cRangeVal, y - cRangeVal, w, h, render_buffer);

  assert(!indices.empty());

  /* now find the correct pixel */
  unsigned int index = 0;
  for (int d = 0; (d < cRangeVal); ++d) {
    for (int a = -d; (a <= d); ++a) {
      for (int b = -d; (b <= d); ++b) {
        index = indices[a + cRangeVal + (b + cRangeVal) * w];
        if (index) {
          a = d = cRangeVal;
          break;
        }
      }
    }
  }

  auto* pik = I->pickmgr.getIdentifier(index);

  if (pik) {
    *pick = *pik;
    if(debug_pick) {
      PRINTFB(G, FB_Scene, FB_Details)
	" SceneClick-Detail: obj %p index %d bond %d\n",
	pick->context.object, pick->src.index, pick->src.bond ENDFB(G);
    }
    // if cPickableNoPick then set object to NULL since nothing picked
    if (pick->src.bond == cPickableNoPick)
      pick->context.object = NULL;
  } else {
    pick->context.object = NULL;
  }

#ifndef PURE_OPENGL_ES_2
  /* Picking changes the Shading model to GL_FLAT,
     we need to change it back to GL_SMOOTH. This is because
     bg_grad() might be called in OrthoDoDraw() before GL 
     settings are set in SceneRender() */
  //      glEnable(GL_COLOR_MATERIAL);
  glShadeModel(SettingGetGlobal_b(G, cSetting_pick_shading) ? GL_FLAT : GL_SMOOTH);
#endif
}

/**
 * Pick all items in a rectangle
 * @param[in,out] smp Defines the (x,y,w,h) rectangle (input), and takes the
 * picking result in `smp->picked` (output)
 */
static
void SceneRenderPickingMultiPick(PyMOLGlobals * G, SceneUnitContext *context, Multipick * smp, GLenum render_buffer){
  CScene *I = G->Scene;
  Picking previous;

  assert(smp->picked.empty());

  auto indices = SceneGetPickIndices(G, context, smp->x, smp->y,
      std::max(1, smp->w), std::max(1, smp->h), render_buffer);

  /* need to scissor this */ 
  for (auto index : indices) {
    const auto* pik = I->pickmgr.getIdentifier(index);
    if (pik) {
      if (pik->src.index != previous.src.index ||
          pik->context.object != previous.context.object) {
        if (pik->context.object->type == cObjectMolecule) {
          smp->picked.push_back(*pik);
        }
        previous = *pik;
      }
    }
  }
#ifndef PURE_OPENGL_ES_2
  /* Picking changes the Shading model to GL_FLAT,
     we need to change it back to GL_SMOOTH. This is because
     bg_grad() might be called in OrthoDoDraw() before GL 
     settings are set in SceneRender() */
  //      glEnable(GL_COLOR_MATERIAL);
  glShadeModel(SettingGetGlobal_b(G, cSetting_pick_shading) ? GL_FLAT : GL_SMOOTH);
#endif
}

void SceneRenderPicking(PyMOLGlobals * G, int stereo_mode, int *click_side, int stereo_double_pump_mono, 
			Picking * pick, int x, int y, Multipick * smp, SceneUnitContext *context,
			GLenum render_buffer){
  CScene *I = G->Scene;

  if (render_buffer == GL_BACK) {
    render_buffer = G->DRAW_BUFFER0;
  }

  SceneSetupGLPicking(G);

  if (!stereo_double_pump_mono){
    switch (stereo_mode) {
    case cStereo_crosseye:
    case cStereo_walleye:
    case cStereo_sidebyside:
      glViewport(I->rect.left, I->rect.bottom, I->Width / 2, I->Height);
      break;
    case cStereo_geowall:
      *click_side = OrthoGetWrapClickSide(G);
      break;
    }
  }
#ifndef PURE_OPENGL_ES_2
  glPushMatrix();           /* 1 */
#endif

  switch (stereo_mode) {
  case cStereo_crosseye:
    ScenePrepareMatrix(G, (*click_side > 0) ? 1 : 2);
    break;
  case cStereo_walleye:
  case cStereo_geowall:
  case cStereo_sidebyside:
    ScenePrepareMatrix(G, (*click_side < 0) ? 1 : 2);
    break;
#ifdef _PYMOL_OPENVR
  case cStereo_openvr:
    ScenePrepareMatrix(G, 0, cStereo_openvr);
    break;
#endif
  }
  G->ShaderMgr->SetIsPicking(true);
  if(pick) {
    SceneRenderPickingSinglePick(G, context, pick, x, y, render_buffer);
  } else if(smp) {
    SceneRenderPickingMultiPick(G, context, smp, render_buffer);
  }
  G->ShaderMgr->SetIsPicking(false);
#ifndef PURE_OPENGL_ES_2
  glPopMatrix();            /* 1 */
#endif
}

/*========================================================================*/
int SceneMultipick(PyMOLGlobals * G, Multipick * smp)
{
  CScene *I = G->Scene;
  int click_side = 0;
  int defer_builds_mode = SettingGetGlobal_i(G, cSetting_defer_builds_mode);

  if(defer_builds_mode == 5)    /* force generation of a pickable version */
    SceneUpdate(G, true);

  if(OrthoGetOverlayStatus(G) || SettingGetGlobal_i(G, cSetting_text))
    SceneRender(G, NULL, 0, 0, NULL, 0, 0, 0, 0);       /* remove overlay if present */
  SceneDontCopyNext(G);
  if(StereoIsAdjacent(G)) {
    if(smp->x > (I->Width / 2))
      click_side = 1;
    else
      click_side = -1;
    smp->x = smp->x % (I->Width / 2);
  }
  SceneRender(G, NULL, 0, 0, smp, 0, 0, click_side, 0);
  SceneDirty(G);
  return (1);
}
