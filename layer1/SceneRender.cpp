/*
 * (c) Schrodinger, Inc.
 */

#include <algorithm>

#include "CGO.h"
#include "Control.h"
#include "Editor.h"
#include "Err.h"
#include "Executive.h"
#include "Feedback.h"
#include "Matrix.h"
#include "P.h"
#include "Picking.h"
#include "PyMOLOptions.h"
#include "Scene.h"
#include "ScenePicking.h"
#include "SceneRay.h"
#include "ShaderMgr.h"
#include "Util.h"
#include "main.h"

#ifdef _PYMOL_OPENVR
#include "OpenVRMode.h"
#endif

/* EXPERIMENTAL VOLUME RAYTRACING DATA */
extern float* rayDepthPixels;
extern int rayVolume, rayWidth, rayHeight;

static void SetDrawBufferForStereo(
    PyMOLGlobals* G, CScene* I, int stereo_mode, int times, int fog_active);
static void SceneDrawStencilInBuffer(
    PyMOLGlobals* G, CScene* I, int stereo_mode);

static void SceneRenderStereoLoop(PyMOLGlobals* G, int timesArg,
    int must_render_stereo, int stereo_mode, bool render_to_texture,
    const Offset2D& pos, const Extent2D& oversizeExtent,
    int stereo_double_pump_mono, int curState, float* normal,
    SceneUnitContext* context, float width_scale, int fog_active,
    bool onlySelections, bool noAA);

static void SceneRenderAA(PyMOLGlobals* G);

static void PrepareViewPortForStereoImpl(PyMOLGlobals* G, CScene* I,
    int stereo_mode, bool offscreen, int times, const Offset2D& pos,
    const Extent2D& oversizeExtent, GLenum draw_mode,
    int position /* left=0, right=1 */);

static void PrepareViewPortForMonoInitializeViewPort(PyMOLGlobals* G, CScene* I,
    int stereo_mode, bool offscreen, int times, const Offset2D& pos,
    const Extent2D& oversizeExtent);

static void PrepareViewPortForStereo(PyMOLGlobals* G, CScene* I,
    int stereo_mode, bool offscreen, int times, const Offset2D& pos,
    const Extent2D& oversizeExtent);

static void PrepareViewPortForStereo2nd(PyMOLGlobals* G, CScene* I,
    int stereo_mode, bool offscreen, int times, const Offset2D& pos,
    const Extent2D& oversizeExtent);

static void InitializeViewPortToScreenBlock(PyMOLGlobals* G, CScene* I,
    const Offset2D& pos, const Extent2D& oversizeExtent, int* stereo_mode,
    float* width_scale);

static void SceneSetPrepareViewPortForStereo(PyMOLGlobals* G,
    PrepareViewportForStereoFuncT prepareViewportForStereo, int times,
    const Offset2D& pos, const Extent2D& oversizeExtent, int stereo_mode,
    float width_scale);

static CGO* GenerateUnitScreenCGO(PyMOLGlobals* G);

static void SceneRenderPostProcessStack(PyMOLGlobals* G);

static int stereo_via_stencil(int stereo_mode)
{
  switch (stereo_mode) {
  case cStereo_stencil_by_row:
  case cStereo_stencil_by_column:
  case cStereo_stencil_checkerboard:
  case cStereo_stencil_custom:
    return true;
  }
  return false;
}

static int render_stereo_blend_into_full_screen(int stereo_mode)
{
  switch (stereo_mode) {
  case cStereo_stencil_by_row:
  case cStereo_stencil_by_column:
  case cStereo_stencil_checkerboard:
  case cStereo_stencil_custom:
  case cStereo_anaglyph:
  case cStereo_dynamic:
  case cStereo_clone_dynamic:
    return true;
  }
  return false;
}

void GridSetViewport(PyMOLGlobals* G, GridInfo* I, int slot)
{
  if (slot)
    I->slot = slot + I->first_slot - 1;
  else
    I->slot = slot;
  /* if we are in grid mode, then prepare the grid slot viewport */
  if (slot < 0) {
    SceneSetViewport(G, I->cur_view);
  } else if (!slot) { /* slot 0 is the full screen */
    Rect2D view{};
    view.offset = Offset2D{};
    view.extent.width = I->cur_view.extent.width / I->n_col;
    view.extent.height = I->cur_view.extent.height / I->n_row;
    if (I->n_col < I->n_row) {
      view.extent.width *= I->n_col;
      view.extent.height *= I->n_col;
    } else {
      view.extent.width *= I->n_row;
      view.extent.height *= I->n_row;
    }
    view.offset.x += I->cur_view.offset.x +
                     (I->cur_view.extent.width - view.extent.width) / 2;
    view.offset.y += I->cur_view.offset.y;
    SceneSetViewport(G, view);
    I->context = ScenePrepareUnitContext(view.extent);
  } else {
    int abs_grid_slot = slot - I->first_slot;
    int grid_col = abs_grid_slot % I->n_col;
    int grid_row = (abs_grid_slot / I->n_col);
    Rect2D view{};
    view.offset.x = (grid_col * I->cur_view.extent.width) / I->n_col;
    view.extent.width =
        ((grid_col + 1) * I->cur_view.extent.width) / I->n_col - view.offset.x;
    view.offset.y = I->cur_view.extent.height -
                    ((grid_row + 1) * I->cur_view.extent.height) / I->n_row;
    view.extent.height =
        (I->cur_view.extent.height -
            ((grid_row) *I->cur_view.extent.height) / I->n_row) -
        view.offset.y;
    view.offset.x += I->cur_view.offset.x;
    view.offset.y += I->cur_view.offset.y;
    I->cur_viewport_size = view.extent;
    SceneSetViewport(G, view);
    I->context = ScenePrepareUnitContext(view.extent);
  }
}

static void glBlendFunc_default()
{
  if (glBlendFuncSeparate) {
    glBlendFuncSeparate(
        GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA, GL_ONE, GL_ONE_MINUS_SRC_ALPHA);
  } else {
    // OpenGL 1.x (e.g. remote desktop on Windows)
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
  }
}

/*========================================================================*/
/* SceneRender: Responsible for rendering the scene, whether its picking
                (SceneRenderPicking) or rendering (SceneRenderStereoLoop).
                It also takes calls anti-aliasing (SceneRenderAA) if
                necessary after rendering and before selection markers.
 */
void SceneRender(PyMOLGlobals* G, const SceneRenderInfo& renderInfo)
{
  /* think in terms of the camera's world */
  CScene* I = G->Scene;
  float normal[4] = {0.0, 0.0, 1.0, 0.0};
  auto aspRat = SceneGetAspectRatio(G);
  double start_time = 0.0;

  float width_scale = 0.0F;
  auto stereo = SettingGet<bool>(G, cSetting_stereo);
  bool use_shaders = SettingGet<bool>(G, cSetting_use_shaders);
  int last_grid_active = I->grid.active;
  I->n_texture_refreshes = 0;
#if defined(_WEBGL) && defined(PYMOL_EVAL)
  if (!OrthoEvalCheck(G))
    return;
#endif
  PRINTFD(G, FB_Scene)
  " SceneRender: entered. pick %p x %d y %d smp %p\n", (void*) renderInfo.pick,
      renderInfo.mousePos.x, renderInfo.mousePos.y,
      (void*) renderInfo.sceneMultipick ENDFD;

  G->ShaderMgr->Check_Reload();

  auto const last_grid_shape = std::array<int, 2>{I->grid.n_col, I->grid.n_row};

  auto grid_mode = SettingGet<GridMode>(G, cSetting_grid_mode);
  int grid_size = 0;
  if (grid_mode != GridMode::NoGrid) {
    grid_size = SceneGetGridSize(G, grid_mode);
    GridUpdate(&I->grid, aspRat, grid_mode, grid_size);
    if (I->grid.active)
      aspRat *= I->grid.asp_adjust;
  } else {
    I->grid.active = false;
  }

  auto const grid_shape = std::array<int, 2>{I->grid.n_col, I->grid.n_row};

  if (last_grid_shape != grid_shape &&
      SettingGet<BgGradient>(G, cSetting_bg_gradient) == BgGradient::Grid) {
    OrthoBackgroundTextureNeedsUpdate(G);
  }

  if (last_grid_active != I->grid.active || grid_size != I->last_grid_size) {
    G->ShaderMgr->ResetUniformSet();
  }
  I->last_grid_size = grid_size;
  G->ShaderMgr->FreeAllVBOs();
#ifndef _PYMOL_IOS
  SceneUpdateAnimation(G);
#endif

  auto render_buffer = SceneMustDrawBoth(G) ? GL_BACK_LEFT : G->DRAW_BUFFER0;

  int stereo_mode = I->StereoMode;
  bool postprocessOnce{false};
  switch (stereo_mode) {
  case cStereo_walleye:
  case cStereo_crosseye:
    aspRat = aspRat / 2;
  case cStereo_sidebyside:
  case cStereo_anaglyph:
    postprocessOnce = stereo;
    break;
  default:
    postprocessOnce = !stereo;
  }
  if (G->HaveGUI && G->ValidContext) {

    if (Feedback(G, FB_OpenGL, FB_Debugging))
      PyMOLCheckOpenGLErr("SceneRender checkpoint 0");

    int stereo_double_pump_mono = false;
    bool must_render_stereo =
        (stereo && stereo_mode != 0); // are we doing stereo?
    if (!must_render_stereo) {
      if (G->StereoCapable &&
          SettingGet<int>(G, NULL, NULL, cSetting_stereo_double_pump_mono)) {
        /* force stereo rendering */
        must_render_stereo = true;
        stereo_double_pump_mono = true;
      }
    }
    /* if we seem to be configured for hardware stereo,
       but can't actually do it, then fallback on mono --
       this would happen for instance if fullscreen is stereo-component
       and windowed is not */
    if (must_render_stereo && (stereo_mode < cStereo_crosseye) &&
        !(G->StereoCapable)) {
      must_render_stereo = false;
    }

    /* If we are rendering a stereo_mode that stencils, define the stencil
     * buffer */
    if (must_render_stereo && stereo_via_stencil(stereo_mode)) {
      if (!I->StencilValid) {
        SceneDrawStencilInBuffer(G, I, stereo_mode);
        I->StencilValid = true;
      }
    }

    render_buffer = G->DRAW_BUFFER0; // GL_BACK

    if (must_render_stereo) {
      switch (stereo_mode) {
      case cStereo_quadbuffer: /* hardware stereo */
      case cStereo_clone_dynamic:
      case cStereo_openvr:
        render_buffer = GL_BACK_LEFT;
        break;
      }
    }

    OrthoDrawBuffer(G, render_buffer);

    if (Feedback(G, FB_OpenGL, FB_Debugging))
      PyMOLCheckOpenGLErr("SceneRender checkpoint 1");

    auto view_save = SceneGetViewport(G);
    InitializeViewPortToScreenBlock(G, I, renderInfo.mousePos,
        renderInfo.oversizeExtent, &stereo_mode, &width_scale);

    if (!(renderInfo.pick || renderInfo.sceneMultipick))
      bg_grad(G);

#ifndef _WEBGL
    glLineWidth(SettingGet<float>(G, cSetting_line_width));
#endif
    glEnable(GL_DEPTH_TEST);

    /* get matrixes for unit objects */
#ifndef PURE_OPENGL_ES_2
    if (SettingGet<bool>(G, cSetting_line_smooth)) {
      if (!(renderInfo.pick || renderInfo.sceneMultipick)) {
        glEnable(GL_LINE_SMOOTH);
        glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);
      }
    } else {
      glDisable(GL_LINE_SMOOTH);
    }
    glPointSize(SettingGet<float>(G, cSetting_dot_width));

    if (ALWAYS_IMMEDIATE_OR(!use_shaders)) {
      glEnable(GL_NORMALIZE); /* get rid of this to boost performance */

      glMatrixMode(GL_MODELVIEW);
      glLoadIdentity();
      /* must be done with identity MODELVIEW */
      SceneProgramLighting(G);
    }
#endif
    Extent2D scene_extent{static_cast<std::uint32_t>(I->Width),
        static_cast<std::uint32_t>(I->Height)};
    auto context = ScenePrepareUnitContext(scene_extent);
    /* do standard 3D objects */
    /* Set up the clipping planes */

    int curState = -1;
    if (!SettingGet<bool>(G, cSetting_all_states)) {
      curState = std::max(-1, SettingGet<int>(G, cSetting_state) - 1);
    }
    if (!SettingGet<bool>(G, cSetting_ortho)) {
      double ymax = I->m_view.m_clipSafe().m_front * GetFovWidth(G) / 2.0;
      double ymin = -ymax;
      double xmin = ymin * aspRat;
      double xmax = ymax * aspRat;
      glFrustum44f(I->ProjectionMatrix, xmin, xmax, ymin, ymax,
          stereo_mode == cStereo_openvr ? 0.1f : I->m_view.m_clipSafe().m_front,
          I->m_view.m_clipSafe().m_back);
    } else {
      Extent2D extent;
      extent.height =
          std::max(R_SMALL4, -I->m_view.pos().z) * GetFovWidth(G) / 2.f;
      extent.width = extent.height * aspRat;
      glOrtho44f(I->ProjectionMatrix, -static_cast<GLfloat>(extent.width),
          static_cast<GLfloat>(extent.width),
          -static_cast<GLfloat>(extent.height),
          static_cast<GLfloat>(extent.height), I->m_view.m_clipSafe().m_front,
          I->m_view.m_clipSafe().m_back);
    }

#ifndef PURE_OPENGL_ES_2
    if (ALWAYS_IMMEDIATE_OR(!use_shaders)) {
      glMatrixMode(GL_PROJECTION);
      glLoadMatrixf(I->ProjectionMatrix);
      glMatrixMode(GL_MODELVIEW);
    }
#endif
    ScenePrepareMatrix(G, 0);

    /* get the Z axis vector for sorting transparent objects */

    if (SettingGet<bool>(G, cSetting_transparency_global_sort) &&
        SettingGet<bool>(G, cSetting_transparency_mode)) {
      if (!I->AlphaCGO)
        I->AlphaCGO = CGONew(G);
    } else {
      CGOFree(I->AlphaCGO);
    }

    /* make note of how large pixels are at the origin  */

    I->VertexScale = SceneGetScreenVertexScale(G, NULL);

    /* determine the direction in which we are looking relative */

    /* 2. set the normals to reflect light back at the camera */

    float zAxis[4] = {0.0, 0.0, 1.0, 0.0};
    MatrixInvTransformC44fAs33f3f(
        glm::value_ptr(I->m_view.rotMatrix()), zAxis, normal);
    copy3f(normal, I->ViewNormal);

    if (SettingGet<bool>(G, cSetting_normal_workaround)) {
      I->LinesNormal[0] = 0.0;
      I->LinesNormal[1] = 0.0;
      I->LinesNormal[2] = 1.0;
      /* for versions of GL that don't transform GL_LINES normals */
    } else {
      I->LinesNormal[0] = I->ViewNormal[0];
      I->LinesNormal[1] = I->ViewNormal[1];
      I->LinesNormal[2] = I->ViewNormal[2];
    }

    PRINTFD(G, FB_Scene)
    " SceneRender: matrices loaded. rendering objects...\n" ENDFD;

    /* 1. render all objects */
    if (renderInfo.pick || renderInfo.sceneMultipick) {
      SceneRenderPicking(G, stereo_mode, renderInfo.clickSide,
          stereo_double_pump_mono, renderInfo.pick, renderInfo.mousePos.x,
          renderInfo.mousePos.y, renderInfo.sceneMultipick, &context,
          render_buffer);
    } else {
      int times = 1;
      bool render_to_texture_for_pp{false};
      /* STANDARD RENDERING */

      start_time = UtilGetSeconds(G);

      glEnable(GL_BLEND);
      glBlendFunc_default();

      glEnable(GL_DITHER);

#ifndef PURE_OPENGL_ES_2
      if (ALWAYS_IMMEDIATE_OR(!use_shaders)) {
        glColorMaterial(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE);
        glEnable(GL_COLOR_MATERIAL);
        glShadeModel(
            SettingGet<bool>(G, cSetting_pick_shading) ? GL_FLAT : GL_SMOOTH);

        if (use_shaders) {
          glDisable(GL_ALPHA_TEST);
        } else {
          // for immediate mode labels (with shaders, this would cause the OS X
          // R9 bugs!)
          glAlphaFunc(GL_GREATER, 0.05F);
          glEnable(GL_ALPHA_TEST);
        }

        if (G->Option->multisample)
          glEnable(0x809D); /* GL_MULTISAMPLE_ARB */
        glColor4ub(255, 255, 255, 255);
        glNormal3fv(normal);
      }
#endif

      auto fog_active = SceneSetFog(G);

#ifndef _PYMOL_NO_AA_SHADERS
      if (renderInfo.oversizeExtent.width == 0 &&
          renderInfo.oversizeExtent.height == 0) {
        render_to_texture_for_pp =
            SettingGet<int>(G, cSetting_antialias_shader);
      }
      if (render_to_texture_for_pp) {
        if (!must_render_stereo || postprocessOnce) {
          G->ShaderMgr->bindOffscreen(I->Width, I->Height, &I->grid);
          bg_grad(G);
        }
      }
#endif
      /* rendering for visualization */

      /*** THIS IS AN UGLY EXPERIMENTAL
       *** VOLUME + RAYTRACING COMPOSITION CODE
       ***/
      if (rayVolume && rayDepthPixels) {
        SceneRenderRayVolume(G, I);
        rayVolume--;
      }
      /*** END OF EXPERIMENTAL CODE ***/

      switch (stereo_mode) {
      case cStereo_clone_dynamic:
      case cStereo_dynamic:
        times = 2;
        break;
      }
      PRINTFD(G, FB_Scene)
      " SceneRender: I->StereoMode %d must_render_stereo %d\n    StereoCapable "
      "%d\n",
          stereo_mode, must_render_stereo, G->StereoCapable ENDFD;

      bool onlySelections{false};
      SceneRenderStereoLoop(G, times, must_render_stereo, stereo_mode,
          render_to_texture_for_pp, renderInfo.mousePos,
          renderInfo.oversizeExtent, stereo_double_pump_mono, curState, normal,
          &context, width_scale, fog_active, onlySelections, postprocessOnce);

      if (render_to_texture_for_pp) {
        /* BEGIN rendering the selection markers, should we put all of this into
           a function, so it can be called above as well? */

#ifndef PURE_OPENGL_ES_2
        if (!must_render_stereo || postprocessOnce) {
          SceneSetPrepareViewPortForStereo(G,
              PrepareViewPortForMonoInitializeViewPort, times,
              renderInfo.mousePos, renderInfo.oversizeExtent, stereo_mode,
              width_scale);
          SceneRenderPostProcessStack(G);
        }
#endif
        bool renderToTexture{false};
        onlySelections = true;
        SceneRenderStereoLoop(G, times, must_render_stereo, stereo_mode,
            renderToTexture, renderInfo.mousePos, renderInfo.oversizeExtent,
            stereo_double_pump_mono, curState, normal, &context, width_scale,
            fog_active, onlySelections, postprocessOnce);
      }

#ifndef PURE_OPENGL_ES_2
      if (ALWAYS_IMMEDIATE_OR(!use_shaders)) {
        glDisable(GL_FOG);
        glDisable(GL_LIGHTING);
        glDisable(GL_LIGHT0);
        glDisable(GL_LIGHT1);
        glDisable(GL_COLOR_MATERIAL);
        glDisable(GL_DITHER);
      }
#endif
    }

#ifndef PURE_OPENGL_ES_2
    if (ALWAYS_IMMEDIATE_OR(!use_shaders)) {
      glLineWidth(1.0);
      glDisable(GL_LINE_SMOOTH);
      glDisable(GL_BLEND);
      glDisable(GL_NORMALIZE);
      glDisable(GL_DEPTH_TEST);
      glDisable(GL_ALPHA_TEST);
      if (G->Option->multisample)
        glDisable(0x809D); /* GL_MULTISAMPLE_ARB */
    }
#endif
    SceneSetViewport(G, view_save);

    if (Feedback(G, FB_OpenGL, FB_Debugging))
      PyMOLCheckOpenGLErr("SceneRender final checkpoint");
  }

  PRINTFD(G, FB_Scene)
  " SceneRender: rendering complete.\n" ENDFD;

  if (!(renderInfo.pick ||
          renderInfo.sceneMultipick)) { /* update frames per second field */
    I->LastRender = UtilGetSeconds(G);
    I->ApproxRenderTime = I->LastRender - start_time;

    if (I->CopyNextFlag) {
      start_time = I->LastRender - start_time;
      if ((start_time > 0.10) || (MainSavingUnderWhileIdle()))
        if (!(ControlIdling(G)))
          if (SettingGet<bool>(G, cSetting_cache_display)) {
            if (!I->CopyType) {
              SceneCopy(G, render_buffer, false, false);
            }
          }
    } else {
      I->CopyNextFlag = true;
    }
    if (renderInfo.forceCopy && !(I->CopyType)) {
      SceneCopy(G, render_buffer, true, false);
      I->CopyType = 2; /* do not display force copies */
    }
  }

#ifdef _PYMOL_OPENVR
  if (stereo_mode == cStereo_openvr &&
      (!SettingGet<bool>(G, cSetting_text) ||
          SettingGet<int>(G, cSetting_openvr_gui_text) == 2)) {
    Block* scene_block = I;
    int scene_width = scene_block->rect.right - scene_block->rect.left;
    int scene_height = scene_block->rect.top - scene_block->rect.bottom;
    OpenVRSceneFinish(G, scene_block->rect.left, scene_block->rect.bottom,
        scene_width, scene_height);
  }
#endif

  PRINTFD(G, FB_Scene)
  " SceneRender: leaving...\n" ENDFD;
}

#ifndef _PYMOL_NO_AA_SHADERS
static void AppendCopyWithChangedShader(
    PyMOLGlobals* G, CGO* destCGO, CGO* srcCGO, int frommode, int tomode)
{
  CGO* cgo = CGONew(G);
  CGOAppendNoStop(cgo, srcCGO);
  CGOChangeShadersTo(cgo, frommode, tomode);
  CGOAppendNoStop(destCGO, cgo);
  CGOFreeWithoutVBOs(cgo);
}
#endif

/* SceneRenderAA: renders Anti-aliasing from the I->offscreen_texture texture,
                  depending on the antialias_shader setting, FXAA (1 stage)
                  or SMAA (3 stages) are rendered using framebuffers
   I->offscreen2_fb, I->offscreen3_fb, and into the screen block
 */
void SceneRenderAA(PyMOLGlobals* G)
{
#ifndef _PYMOL_NO_AA_SHADERS
  CScene* I = G->Scene;
  int ok = true;
  glBindFramebufferEXT(
      GL_FRAMEBUFFER_EXT, G->ShaderMgr->default_framebuffer_id);
  if (!I->offscreenCGO) {
    CGO* unitCGO = GenerateUnitScreenCGO(G);
    ok &= unitCGO != NULL;
    if (ok) {
      auto antialias_mode = SettingGet<int>(G, cSetting_antialias_shader);

      I->offscreenCGO = CGONew(G);

      switch (antialias_mode) {
      case 0:
        break;
      case 1: // fxaa
        AppendCopyWithChangedShader(G, I->offscreenCGO, unitCGO,
            GL_DEFAULT_SHADER_WITH_SETTINGS, GL_FXAA_SHADER);
        break;
      default:
        AppendCopyWithChangedShader(G, I->offscreenCGO, unitCGO,
            GL_DEFAULT_SHADER_WITH_SETTINGS, GL_SMAA1_SHADER);
        if (antialias_mode != 3) { // not 1nd Pass as output
          CGODisable(I->offscreenCGO, GL_SMAA1_SHADER);
          AppendCopyWithChangedShader(G, I->offscreenCGO, unitCGO,
              GL_DEFAULT_SHADER_WITH_SETTINGS, GL_SMAA2_SHADER);
          CGODisable(I->offscreenCGO, GL_SMAA2_SHADER);
          if (antialias_mode != 4) { // not 2nd Pass as output
            AppendCopyWithChangedShader(G, I->offscreenCGO, unitCGO,
                GL_DEFAULT_SHADER_WITH_SETTINGS, GL_SMAA3_SHADER);
            CGODisable(I->offscreenCGO, GL_SMAA3_SHADER);
          }
        }
        break;
      }
      CGOStop(I->offscreenCGO);
      CGOFreeWithoutVBOs(unitCGO);
      I->offscreenCGO->use_shader = true;
    } else {
      I->offscreenCGO = NULL;
    }
  }
  if (ok && I->offscreenCGO) {
    CGORender(I->offscreenCGO, NULL, NULL, NULL, NULL, NULL);
    G->ShaderMgr->Disable_Current_Shader();
    glBindTexture(GL_TEXTURE_2D, 0);
    glEnable(GL_DEPTH_TEST);
    glBindFramebufferEXT(
        GL_FRAMEBUFFER_EXT, G->ShaderMgr->default_framebuffer_id);
  }
#endif
}

static void SceneRenderAllObject(PyMOLGlobals* G, CScene* I,
    SceneUnitContext* context, RenderInfo* info, float* normal, int state,
    pymol::CObject* obj, GridInfo* grid, int* slot_vla, int fat)
{
  if (!SceneGetDrawFlag(grid, slot_vla, obj->grid_slot))
    return;

  auto use_shader = info->use_shaders;

#ifndef _WEBGL
  glLineWidth(fat ? 3.0 : 1.0);
#endif

  switch (obj->getRenderContext()) {
  case pymol::RenderContext::UnitWindow:
    // e.g. Gadgets/Ramps
    {
      float projSave[16];
      copy44f(I->ProjectionMatrix, projSave);

      if (grid->active) {
        context = &grid->context;
      }

      glOrtho44f(I->ProjectionMatrix, context->unit_left, context->unit_right,
          context->unit_top, context->unit_bottom, context->unit_front,
          context->unit_back);

#ifndef PURE_OPENGL_ES_2
      if (ALWAYS_IMMEDIATE_OR(!use_shader)) {
        glPushAttrib(GL_LIGHTING_BIT);

        glMatrixMode(GL_PROJECTION);
        glLoadMatrixf(I->ProjectionMatrix);

        glMatrixMode(GL_MODELVIEW);
        glPushMatrix();
        glLoadIdentity();

        float vv[4] = {0.f, 0.f, -1.f, 0.f}, dif[4] = {1.f, 1.f, 1.f, 1.f};
        glLightfv(GL_LIGHT0, GL_POSITION, vv);
        glLightfv(GL_LIGHT0, GL_DIFFUSE, dif);

        glNormal3f(0.0F, 0.0F, 1.0F);
      }
#endif

      info->state = ObjectGetCurrentState(obj, false);
      obj->render(info);

      copy44f(projSave, I->ProjectionMatrix);

#ifndef PURE_OPENGL_ES_2
      if (ALWAYS_IMMEDIATE_OR(!use_shader)) {
        glMatrixMode(GL_PROJECTION);
        glLoadMatrixf(I->ProjectionMatrix);

        glMatrixMode(GL_MODELVIEW);
        glPopMatrix();

        glPopAttrib();
      }
#endif
    }
    break;
  case pymol::RenderContext::Camera: /* context/grid 0 is all slots */
  default:
    ScenePushModelViewMatrix(G);

#ifndef PURE_OPENGL_ES_2
    if (normal && Feedback(G, FB_OpenGL, FB_Debugging))
      glNormal3fv(normal);
#endif

    if (!grid->active ||
        grid->mode == GridMode::NoGrid ||
        grid->mode == GridMode::ByObject) {
      info->state = ObjectGetCurrentState(obj, false);
      obj->render(info);
    } else if (grid->slot) {
      if (grid->mode == GridMode::ByObjectStates) {
        if ((info->state = state + grid->slot - 1) >= 0)
          obj->render(info);
      } else if (grid->mode == GridMode::ByObjectByState) {
        info->state = grid->slot - obj->grid_slot - 1;
        if (info->state >= 0 && info->state < obj->getNFrame())
          obj->render(info);
      }
    }

    ScenePopModelViewMatrix(G, !use_shader);
    break;
  }
}

/*========================================================================
 * SceneRenderAll: Renders all CObjects in the scene
 *
 * context: context info
 * normal: initial normal (for immediate mode)
 * pass: which pass (opaque, antialias, transparent)
 * fat: wide lines (i.e., for picking)
 * width_scale: specifies width_scale and sampling
 * grid: grid information
 * dynamic_pass: for specific stereo modes dynamic and clone_dynamic
 * which: enum specifying which objects (AllObjects, OnlyGadgets,
 * OnlyNonGadgets, GadgetsLast)
 */
void SceneRenderAll(PyMOLGlobals* G, SceneUnitContext* context, float* normal,
    PickColorManager* pickmgr, RenderPass pass, int fat, float width_scale,
    GridInfo* grid, int dynamic_pass, SceneRenderWhich which_objects)
{
  CScene* I = G->Scene;
  int state = SceneGetState(G);
  RenderInfo info;
#if defined(_WEBGL) && defined(PYMOL_EVAL)
  if (!OrthoEvalCheck(G))
    return;
#endif
  info.pick = pickmgr;
  info.pass = pass;
  info.vertex_scale = I->VertexScale;
  info.fog_start = I->FogStart;
  info.fog_end = I->FogEnd;
  info.front = I->m_view.m_clipSafe().m_front;
  info.use_shaders = SettingGetGlobal_b(G, cSetting_use_shaders);
#if defined(_PYMOL_IOS) && !defined(_WEBGL)
  /* For now, on IOS, just crank up the sampling for text,
     TODO : need to figure out where this change should really be */
  info.sampling = 2;
#else
  info.sampling = 1;
#endif
  info.alpha_cgo = I->AlphaCGO;
  info.ortho = SettingGetGlobal_b(G, cSetting_ortho);
  if (I->StereoMode && dynamic_pass && (!info.pick)) {
    int stereo_mode = SettingGetGlobal_i(G, cSetting_stereo_mode);
    switch (stereo_mode) {
    case cStereo_dynamic:
    case cStereo_clone_dynamic:
      info.line_lighting = true;
      break;
    }
  }

  if (I->StereoMode) {
    float buffer;
    float stAng, stShift;
    stAng = SettingGetGlobal_f(G, cSetting_stereo_angle);
    stShift = SettingGetGlobal_f(G, cSetting_stereo_shift);
    stShift = (float) (stShift * fabs(I->m_view.pos().z) / 100.0);
    stAng =
        (float) (stAng * atan(stShift / fabs(I->m_view.pos().z)) * 90.0 / cPI);
    buffer = fabs(I->Width * I->VertexScale * tan(cPI * stAng / 180.0));
    info.stereo_front = I->m_view.m_clipSafe().m_front + buffer;
  } else {
    info.stereo_front = I->m_view.m_clipSafe().m_front;
  }

  info.back = I->m_view.m_clipSafe().m_back;
  SceneGetViewNormal(G, info.view_normal);

  if (info.alpha_cgo && (pass == RenderPass::Opaque)) {
    CGOReset(info.alpha_cgo);
    CGOSetZVector(
        info.alpha_cgo, I->ModMatrix[2], I->ModMatrix[6], I->ModMatrix[10]);
  }

  if (SettingGetGlobal_b(G, cSetting_dynamic_width)) {
    info.dynamic_width = true;
    info.dynamic_width_factor =
        SettingGetGlobal_f(G, cSetting_dynamic_width_factor);
    info.dynamic_width_min = SettingGetGlobal_f(G, cSetting_dynamic_width_min);
    info.dynamic_width_max = SettingGetGlobal_f(G, cSetting_dynamic_width_max);
  }

  if (width_scale != 0.0F) {
    info.width_scale_flag = true;
    info.width_scale = width_scale;
    info.sampling = (int) info.width_scale;
    if (info.sampling < 1)
      info.sampling = 1;
  }
  {
    auto slot_vla = I->m_slots.data();
    switch (which_objects) {
    case SceneRenderWhich::AllObjects:
      for (auto obj : I->Obj) {
        /* EXPERIMENTAL RAY-VOLUME COMPOSITION CODE */
        if (!rayVolume || obj->type == cObjectVolume) {
          SceneRenderAllObject(
              G, I, context, &info, normal, state, obj, grid, slot_vla, fat);
        }
      }
      break;
    case SceneRenderWhich::OnlyGadgets:
      for (auto obj : I->GadgetObjs) {
        SceneRenderAllObject(
            G, I, context, &info, normal, state, obj, grid, slot_vla, fat);
      }
      break;
    case SceneRenderWhich::OnlyNonGadgets:
      for (auto obj : I->NonGadgetObjs) {
        // ObjectGroup used to have fRender = NULL
        if (obj->type != cObjectGroup) {
          SceneRenderAllObject(
              G, I, context, &info, normal, state, obj, grid, slot_vla, fat);
        }
      }
      break;
    case SceneRenderWhich::GadgetsLast:
      // Gadgets Last
      for (auto obj : I->NonGadgetObjs) {
        /* EXPERIMENTAL RAY-VOLUME COMPOSITION CODE */
        if (obj->type !=
                cObjectGroup && // ObjectGroup used to have fRender = NULL
            (!rayVolume || obj->type == cObjectVolume)) {
          SceneRenderAllObject(
              G, I, context, &info, normal, state, obj, grid, slot_vla, fat);
        }
      }
      for (auto obj : I->GadgetObjs) {
        SceneRenderAllObject(
            G, I, context, &info, normal, state, obj, grid, slot_vla, fat);
      }
      break;
    }
  }

  if (info.alpha_cgo) {
    CGOStop(info.alpha_cgo);
    /* this only works when all objects are rendered in the same frame of
     * reference */
    if (pass == RenderPass::Transparent) {
      CGORenderAlpha(info.alpha_cgo, &info, 0);
    }
  }
}
/*==================================================================================*/
/* DoRendering: This is the function that is responsible for looping through
   each rendering pass (opaque, then antialiased, then transparent) for each
   grid slot (only one grid slot if full screen).  It also implements
   transparency_mode 3, (weighted, blended order-independent transparency) which
   renders the opaque/antialiased passes to the offscreen texture with a depth
   texture, renders the transparent pass to the OIT offscreen texture, calls
   OIT_copy to copy the opaque to the screen (if necessary, i.e., not already
   rendering to AA texture), and then calls the OIT rendering pass that computes
   the resulting image. This function also renders only the selections
   (onlySelections) for all grids or the full screen.
 */
static void DoRendering(PyMOLGlobals* G, CScene* I, GridInfo* grid, int times,
    int curState, float* normal, SceneUnitContext* context, float width_scale,
    bool onlySelections, bool excludeSelections)
{
  const RenderPass passes[] = {
      RenderPass::Opaque, RenderPass::Antialias, RenderPass::Transparent};
  bool use_shaders = (bool) SettingGetGlobal_b(G, cSetting_use_shaders);
  bool t_mode_3_os =
      use_shaders && SettingGetGlobal_i(G, cSetting_transparency_mode) == 3;
  bool t_mode_3 = !onlySelections && t_mode_3_os;
  GLint currentFrameBuffer;

#if !defined(PURE_OPENGL_ES_2) || defined(_WEBGL)
  if (t_mode_3) {
    glGetIntegerv(GL_FRAMEBUFFER_BINDING, &currentFrameBuffer);
    // currentFrameBuffer: 0 - rendering to screen, need to render opaque to
    // offscreen buffer
    //                     non-0 - already rendering to AA texture, need to use
    //                     I->offscreen_depth_rb
    //                             transparent (OIT) pass
    // In the case of jymol the currentFramebuffer is not 0 so we are checking
    // against the default framebuffer
    if (currentFrameBuffer == G->ShaderMgr->default_framebuffer_id) {
      G->ShaderMgr->bindOffscreen(I->Width, I->Height, &I->grid);
      bg_grad(G);
    }
    glEnable(GL_DEPTH_TEST);
  }
#endif
  if (grid->active) { // && !offscreen)
    grid->cur_view = SceneGetViewport(G);
  }
  {
    int slot;
    bool cont = true;
    bool t_first_pass = true;
    G->ShaderMgr->stereo_draw_buffer_pass = 0;
    for (auto pass :
        passes) { /* render opaque, then antialiased, then transparent... */
      if (!cont) {
        break;
      }
#if !defined(PURE_OPENGL_ES_2) || defined(_WEBGL)
      if (t_mode_3 && pass == RenderPass::Transparent) {
        G->ShaderMgr->Disable_Current_Shader();
        int drawbuf = 1;
        if (TM3_IS_ONEBUF) {
          if (!t_first_pass) {
            G->ShaderMgr->stereo_draw_buffer_pass = 1;
          }
          drawbuf = t_first_pass ? 1 : 2;
        }
        G->ShaderMgr->bindOffscreenOIT(I->Width, I->Height, drawbuf);
        G->ShaderMgr->oit_pp->bindRT(
            drawbuf); // for transparency pass, render to OIT texture
        if (currentFrameBuffer == G->ShaderMgr->default_framebuffer_id) {
          SceneInitializeViewport(G, true);
        }
      }
#endif
      for (slot = 0; slot <= grid->last_slot; slot++) {
        if (grid->active) {
          GridSetViewport(G, grid, slot);
        } else if (slot) {
          break; // if grid is off, then just get out of loop after 1st pass
                 // (full screen)
        }
        /* render picked atoms */
        /* render the debugging CGO */
#ifdef PURE_OPENGL_ES_2
        if (!onlySelections) {
          EditorRender(G, curState);
          CGORender(G->DebugCGO, NULL, NULL, NULL, NULL, NULL);
        }
#else
        if (!use_shaders)
          glPushMatrix(); /* 2 */
        if (!onlySelections && !t_mode_3)
          EditorRender(G, curState);
        if (!use_shaders) {
          glPopMatrix();  /* 1 */
          glPushMatrix(); /* 2 */
        }
        if (!onlySelections) {
          if (!use_shaders)
            glNormal3fv(normal);
          CGORender(G->DebugCGO, NULL, NULL, NULL, NULL, NULL);
        }
        if (!use_shaders) {
          glPopMatrix();  /* 1 */
          glPushMatrix(); /* 2 */
        }
#endif
        /* render all objects */
        if (!onlySelections) {
#if !defined(PURE_OPENGL_ES_2) || defined(_WEBGL)
          if (t_mode_3) {
            if (pass == RenderPass::Opaque) {
              EditorRender(G, curState);
            }
            // transparency-mode == 3 render all objects for this pass
            SceneRenderAll(G, context, normal, NULL, pass, false, width_scale,
                grid, times, SceneRenderWhich::OnlyNonGadgets); // opaque
          } else {
#else
          {
#endif
            // transparency-mode != 3 render all objects for each pass
            for (const auto pass2 :
                passes) { /* render opaque, then antialiased, then
                             transparent... */
              SceneRenderAll(G, context, normal, NULL, pass2, false,
                  width_scale, grid, times, SceneRenderWhich::GadgetsLast);
            }
            cont = false;
          }
        } else if (t_mode_3_os && pass == RenderPass::Opaque) {
          // onlySelections and t_mode_3, render only gadgets
          glEnable(GL_BLEND); // need to blend for the text onto the gadget
                              // background
          glBlendFunc_default();

          SceneRenderAll(G, context, normal, NULL,
              RenderPass::Transparent /* gadgets render in transp pass */,
              false, width_scale, grid, times, SceneRenderWhich::OnlyGadgets);
          glDisable(GL_BLEND);
        }
#ifdef PURE_OPENGL_ES_2
        if (!excludeSelections) {
          if (!grid->active ||
              slot > 0) { /* slot 0 is the full screen in grid mode, so don't
                             render selections */
            int s = grid->active && grid->mode == GridMode::ByObject ? slot : 0;
            ExecutiveRenderSelections(G, curState, s, grid);
          }
        }
#else
        if (!use_shaders) {
          glPopMatrix(); /* 1 */
          /* render selections */
          glPushMatrix(); /* 2 */
          glNormal3fv(normal);
        }
        if (!t_mode_3 && !excludeSelections) {
          if (!grid->active ||
              slot > 0) { /* slot 0 is the full screen in grid mode, so don't
                             render selections */
            int s = grid->active && grid->mode == GridMode::ByObject ? slot : 0;
            ExecutiveRenderSelections(G, curState, s, grid);
          }
        }
        if (!use_shaders) {
          glPopMatrix(); /* 1 */
        }
#endif
      } // end slot loop
      if (TM3_IS_ONEBUF) {
        if (t_mode_3 && pass == RenderPass::Transparent && t_first_pass) {
          pass = RenderPass::Antialias;
          t_first_pass = false;
          continue;
        }
      }
#if !defined(PURE_OPENGL_ES_2) || defined(_WEBGL)
      if (t_mode_3 && pass == RenderPass::Transparent) {
        glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, currentFrameBuffer);
        glBindTexture(GL_TEXTURE_2D, 0);
        if (grid->active)
          GridSetViewport(G, grid, -1);
        if (currentFrameBuffer ==
            G->ShaderMgr
                ->default_framebuffer_id) { // if rendering to screen, need to
                                            // render offscreen opaque to screen
          SceneInitializeViewport(G, false);
          if (!I->offscreenOIT_CGO_copy) {
            // TODO G->ShaderMgr->Reload_Copy_Shaders();
            I->offscreenOIT_CGO_copy = GenerateUnitScreenCGO(G);
            CGOChangeShadersTo(I->offscreenOIT_CGO_copy,
                GL_DEFAULT_SHADER_WITH_SETTINGS, GL_OIT_COPY_SHADER);
            I->offscreenOIT_CGO_copy->use_shader = true;
          }
          CGORender(I->offscreenOIT_CGO_copy, NULL, NULL, NULL, NULL, NULL);
        }
        if (!I->offscreenOIT_CGO) {
          I->offscreenOIT_CGO = GenerateUnitScreenCGO(G);
          CGOChangeShadersTo(I->offscreenOIT_CGO,
              GL_DEFAULT_SHADER_WITH_SETTINGS, GL_OIT_SHADER);
          I->offscreenOIT_CGO->use_shader = true;
        }
        CGORender(I->offscreenOIT_CGO, NULL, NULL, NULL, NULL, NULL);

        glBlendFunc_default();

        if ((currentFrameBuffer == G->ShaderMgr->default_framebuffer_id) &&
            t_mode_3) {
          // onlySelections and t_mode_3, render only gadgets
          SceneRenderAll(G, context, normal, NULL,
              RenderPass::Transparent /* gadgets render in transp pass */,
              false, width_scale, grid, times, SceneRenderWhich::OnlyGadgets);
        }

        glDisable(GL_BLEND);
        glDepthMask(GL_TRUE);

        if (!excludeSelections) {
          grid->cur_view = SceneGetViewport(G);
          for (slot = 0; slot <= grid->last_slot; slot++) {
            if (grid->active) {
              GridSetViewport(G, grid, slot);
            }
            if (!grid->active ||
                slot > 0) { /* slot 0 is the full screen in grid mode, so don't
                               render selections */
              int s =
                  grid->active && grid->mode == GridMode::ByObject ? slot : 0;
              ExecutiveRenderSelections(G, curState, s, grid);
            }
          }
        }
      }
#endif
    }
  }
  if (grid->active)
    GridSetViewport(G, grid, -1);
}

/*==================================================================================*/
/* SceneRenderStereoLoop: This is the function that is responsible for rendering
   all objects either in a monoscopic or stereo display.  It prepares the
   viewport, offscreen textures (if necessary), draws the background (if
   necessary) and calls DoRendering either once (for monoscopic) or twice (for
   stereo)
 */
void SceneRenderStereoLoop(PyMOLGlobals* G, int timesArg,
    int must_render_stereo, int stereo_mode, bool render_to_texture,
    const Offset2D& pos, const Extent2D& oversizeExtent,
    int stereo_double_pump_mono, int curState, float* normal,
    SceneUnitContext* context, float width_scale, int fog_active,
    bool onlySelections, bool offscreenPrepared)
{
  CScene* I = G->Scene;
  int times = timesArg;
  bool shouldPrepareOffscreen =
      !onlySelections && render_to_texture && !offscreenPrepared;
  bool use_shaders = (bool) SettingGetGlobal_b(G, cSetting_use_shaders);

  // only cStereo_clone_dynamic and cStereo_dynamic has times=2, otherwise
  // times=1
  while (times--) {
    if (must_render_stereo) {
      bool anaglyph = G->ShaderMgr && stereo_mode == cStereo_anaglyph;
      /* STEREO RENDERING (real or double-pumped) */
      PRINTFD(G, FB_Scene)
      " SceneRender: left hand stereo...\n" ENDFD;

      /* LEFT HAND STEREO */
      if (anaglyph) {
        G->ShaderMgr->stereo_flag = -1; // left eye
        G->ShaderMgr->stereo_blend = 0;
      }

#ifdef _PYMOL_OPENVR
      int savedWidth, savedHeight;
      if (stereo_mode == cStereo_openvr) {
        savedWidth = I->Width;
        savedHeight = I->Height;
        OpenVRGetWidthHeight(G, &I->Width, &I->Height);
      }
#endif

      SceneSetPrepareViewPortForStereo(G, PrepareViewPortForStereo, times, pos,
          oversizeExtent, stereo_mode, width_scale);

      if (!shouldPrepareOffscreen) {
        PrepareViewPortForStereo(
            G, I, stereo_mode, render_to_texture, times, pos, oversizeExtent);
      }
#ifndef PURE_OPENGL_ES_2
      if (use_shaders)
        glPushMatrix(); // 1
      if (shouldPrepareOffscreen) {
        G->ShaderMgr->bindOffscreen(I->Width, I->Height, &I->grid);
        bg_grad(G);
      }
#endif
      ScenePrepareMatrix(G, stereo_double_pump_mono ? 0 : 1, stereo_mode);
      DoRendering(G, I, &I->grid, times, curState, normal, context, width_scale,
          onlySelections, render_to_texture);

#ifndef PURE_OPENGL_ES_2
      if (use_shaders)
        glPopMatrix(); // 0
#endif

#ifdef _PYMOL_OPENVR
      // TODO Check if this is the correct place for this block. In openvr
      // branch, was last block in DoHandedStereo.
      if (stereo_mode == cStereo_openvr) {
        OpenVRDraw(G);
        OpenVREyeFinish(G);
      }
#endif

      PRINTFD(G, FB_Scene)
      " SceneRender: right hand stereo...\n" ENDFD;
      if (shouldPrepareOffscreen) {
        SceneRenderPostProcessStack(G);
      }

      /* RIGHT HAND STEREO */
      if (anaglyph) {
        G->ShaderMgr->stereo_flag = 1; // right eye
        G->ShaderMgr->stereo_blend =
            render_stereo_blend_into_full_screen(stereo_mode);
      }
      SceneSetPrepareViewPortForStereo(G, PrepareViewPortForStereo2nd, times,
          pos, oversizeExtent, stereo_mode, width_scale);
      if (!shouldPrepareOffscreen) {
        PrepareViewPortForStereo2nd(
            G, I, stereo_mode, render_to_texture, times, pos, oversizeExtent);
      }
#ifndef PURE_OPENGL_ES_2
      if (!use_shaders)
        glPushMatrix(); // 1
      if (shouldPrepareOffscreen) {
        G->ShaderMgr->bindOffscreen(I->Width, I->Height, &I->grid);
      }
      if (shouldPrepareOffscreen ||
          (stereo_mode == cStereo_quadbuffer && !onlySelections) // PYMOL-2342
      ) {
        bg_grad(G);
      }
#endif
      ScenePrepareMatrix(G, stereo_double_pump_mono ? 0 : 2, stereo_mode);
      glClear(GL_DEPTH_BUFFER_BIT);
      DoRendering(G, I, &I->grid, times, curState, normal, context, width_scale,
          onlySelections, render_to_texture);
      if (anaglyph) {
        G->ShaderMgr->stereo_flag = 0;
        G->ShaderMgr->stereo_blend = 0;
      }
#ifndef PURE_OPENGL_ES_2
      if (!use_shaders)
        glPopMatrix(); // 0
#endif

#ifdef _PYMOL_OPENVR
      if (stereo_mode == cStereo_openvr) {
        // TODO Check if this is the correct place for this block. In openvr
        // branch, was last block in DoHandedStereo.
        OpenVRDraw(G);
        OpenVREyeFinish(G);

        I->Width = savedWidth;
        I->Height = savedHeight;
      }
#endif

      /* restore draw buffer */
      if (shouldPrepareOffscreen) {
        SceneRenderPostProcessStack(G);
      }
      SetDrawBufferForStereo(G, I, stereo_mode, times, fog_active);
    } else {
      if (G->ShaderMgr) {
        G->ShaderMgr->stereo_flag = 0;
        G->ShaderMgr->stereo_blend = 0;
      }
      /* MONOSCOPING RENDERING (not double-pumped) */
      if (!I->grid.active && render_to_texture) {
        SceneSetViewport(G, 0, 0, I->Width, I->Height);
        if (!onlySelections)
          bg_grad(G);
      }
      if (Feedback(G, FB_OpenGL, FB_Debugging))
        PyMOLCheckOpenGLErr("Before mono rendering");
      SceneSetPrepareViewPortForStereo(G,
          PrepareViewPortForMonoInitializeViewPort, times, pos, oversizeExtent,
          stereo_mode, width_scale);
      DoRendering(G, I, &I->grid, times, curState, normal, context, width_scale,
          onlySelections, render_to_texture);
      if (Feedback(G, FB_OpenGL, FB_Debugging))
        PyMOLCheckOpenGLErr("during mono rendering");
    }
  }
}

void PrepareViewPortForMonoInitializeViewPort(PyMOLGlobals* G, CScene* I,
    int stereo_mode, bool offscreen, int times, const Offset2D& pos,
    const Extent2D& oversizeExtent)
{
  float width_scale;
  InitializeViewPortToScreenBlock(
      G, I, pos, oversizeExtent, &stereo_mode, &width_scale);
}

void PrepareViewPortForStereo(PyMOLGlobals* G, CScene* I, int stereo_mode,
    bool offscreen, int times, const Offset2D& pos,
    const Extent2D& oversizeExtent)
{
  PrepareViewPortForStereoImpl(G, I, stereo_mode, offscreen, times, pos,
      oversizeExtent, GL_BACK_LEFT, 0);
}

void PrepareViewPortForStereo2nd(PyMOLGlobals* G, CScene* I, int stereo_mode,
    bool offscreen, int times, const Offset2D& pos,
    const Extent2D& oversizeExtent)
{
  PrepareViewPortForStereoImpl(G, I, stereo_mode, offscreen, times, pos,
      oversizeExtent, GL_BACK_RIGHT, 1);
}

void InitializeViewPortToScreenBlock(PyMOLGlobals* G, CScene* I,
    const Offset2D& pos, const Extent2D& oversizeExtent, int* stereo_mode,
    float* width_scale)
{
  if (oversizeExtent.width != 0 && oversizeExtent.height != 0) {
    Rect2D want_view{};
    want_view.offset.x = static_cast<std::int32_t>(I->rect.left + pos.x);
    want_view.offset.y = static_cast<std::int32_t>(I->rect.bottom + pos.y);
    want_view.extent = oversizeExtent;
    SceneSetViewport(G, want_view);
    auto got_view = SceneGetViewport(G);
#ifndef _WEBGL
    if (want_view != got_view) {
      PRINTFB(G, FB_Scene, FB_Warnings)
      "Scene-Warning: glViewport failure.\n" ENDFB(G);
    }
#endif
    switch (*stereo_mode) {
    case cStereo_geowall:
      *stereo_mode = 0;
      break;
    }
    *width_scale = ((float) (oversizeExtent.width)) / I->Width;
  } else {
    Rect2D view{};
    view.offset = Offset2D{I->rect.left, I->rect.bottom};
    view.extent = SceneGetExtent(G);
    SceneSetViewport(G, view);
  }
}

void SceneSetPrepareViewPortForStereo(PyMOLGlobals* G,
    PrepareViewportForStereoFuncT prepareViewportForStereo, int times,
    const Offset2D& pos, const Extent2D& oversizeExtent, int stereo_mode,
    float width_scale)
{
  CScene* I = G->Scene;
  I->vp_prepareViewPortForStereo = prepareViewportForStereo;
  I->vp_times = times;
  I->vp_pos = pos;
  I->vp_oversize = oversizeExtent;
  I->vp_stereo_mode = stereo_mode;
  I->vp_width_scale = width_scale;
}

/* PrepareViewPortForStereoImpl : sets up viewport and GL state for stereo_modes
 */
void PrepareViewPortForStereoImpl(PyMOLGlobals* G, CScene* I, int stereo_mode,
    bool offscreen, int times, const Offset2D& pos,
    const Extent2D& oversizeExtent, GLenum draw_mode,
    int position /* left=0, right=1 */)
{
  int position_inv = 1 - position;
  switch (stereo_mode) {
  case cStereo_quadbuffer: /* hardware */
    OrthoDrawBuffer(G, draw_mode);
    SceneSetViewport(G, I->rect.left, I->rect.bottom, I->Width, I->Height);
    break;
  case cStereo_crosseye: /* side by side, crosseye */
    if (offscreen) {
      SceneSetViewport(
          G, position_inv * I->Width / 2, 0, I->Width / 2, I->Height);
    } else if (oversizeExtent.width != 0 && oversizeExtent.height != 0) {
      SceneSetViewport(G,
          I->rect.left + (position_inv * oversizeExtent.width / 2) + pos.x,
          I->rect.bottom + pos.y, oversizeExtent.width / 2,
          oversizeExtent.height);
    } else {
      SceneSetViewport(G, I->rect.left + (position_inv * I->Width / 2),
          I->rect.bottom, I->Width / 2, I->Height);
    }
    break;
  case cStereo_walleye:
  case cStereo_sidebyside:
    if (offscreen) {
      SceneSetViewport(G, position * I->Width / 2, 0, I->Width / 2, I->Height);
    } else if (oversizeExtent.width != 0 && oversizeExtent.height != 0) {
      SceneSetViewport(G,
          I->rect.left + (position * oversizeExtent.width / 2) + pos.x,
          I->rect.bottom + pos.y, oversizeExtent.width / 2,
          oversizeExtent.height);
    } else {
      SceneSetViewport(G, I->rect.left + (position * I->Width / 2),
          I->rect.bottom, I->Width / 2, I->Height);
    }
    break;
  case cStereo_geowall:
    if (offscreen) {
      SceneSetViewport(G, position * I->Width / 2, 0, I->Width / 2, I->Height);
    } else {
      SceneSetViewport(G, I->rect.left + (position * G->Option->winX / 2),
          I->rect.bottom, I->Width, I->Height);
    }
    break;
  case cStereo_stencil_by_row:
  case cStereo_stencil_by_column:
  case cStereo_stencil_checkerboard:
    if (I->StencilValid) {
      glStencilFunc(GL_EQUAL, position_inv, 1);
      glStencilOp(GL_KEEP, GL_KEEP, GL_KEEP);
      glEnable(GL_STENCIL_TEST);
    }
    break;
  case cStereo_stencil_custom:
    break;
  case cStereo_anaglyph:
    /* glClear(GL_ACCUM_BUFFER_BIT); */
#ifdef _PYMOL_IOS
    if (position)
      glClear(GL_DEPTH_BUFFER_BIT);
    glEnable(GL_BLEND);
    glBlendEquation(GL_FUNC_ADD);
    glColorMask(position_inv, position, position, true);
#else

    if (TM3_IS_ONEBUF) {
      glColorMask(position_inv, position, position, true);
    } else {
      if (GLEW_VERSION_3_0 &&
          SettingGetGlobal_i(G, cSetting_transparency_mode) == 3) {
        // if GL 3.0 is available, use glColorMaski to mask only first draw
        // buffer for anaglyph in transparency_mode 3
        glColorMaski(0, position_inv, position, position, true);
      } else {
        glColorMask(position_inv, position, position, true);
      }
    }

    if (position)
      glClear(GL_DEPTH_BUFFER_BIT);
#endif
    break;
#ifndef PURE_OPENGL_ES_2
  case cStereo_clone_dynamic:
    if (position_inv) {
      glClear(GL_ACCUM_BUFFER_BIT);
      OrthoDrawBuffer(G, GL_BACK_LEFT);
      if (times) {
        float dynamic_strength =
            SettingGetGlobal_f(G, cSetting_stereo_dynamic_strength);
        float vv[4] = {0.75F, 0.75F, 0.75F, 1.0F};
        vv[0] = dynamic_strength;
        vv[1] = dynamic_strength;
        vv[2] = dynamic_strength;
        glMaterialfv(GL_FRONT_AND_BACK, GL_EMISSION, vv);
        glAccum(GL_ADD, 0.5);
        glDisable(GL_FOG);
      }
    } else {
      GLenum err;
      if (times) {
        glAccum(GL_ACCUM, -0.5);
      } else {
        glAccum(GL_ACCUM, 0.5);
      }
      if ((err = glGetError())) {
        PRINTFB(G, FB_Scene, FB_Errors)
        "Stereo Error 0x%x: stereo_mode=12 clone_dynamic requires access to "
        "the accumulation buffer,\n"
        "you need to start PyMOL with the -t argument, setting back to "
        "default\n",
            err ENDFB(G);
        SettingSetGlobal_i(G, cSetting_stereo_mode, cStereo_crosseye);
        SceneSetStereo(G, 0);
        return;
      }
      glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    }
    break;
  case cStereo_dynamic:
    if (position_inv) {
      if (times) {
        float dynamic_strength =
            SettingGetGlobal_f(G, cSetting_stereo_dynamic_strength);
        float vv[4] = {0.75F, 0.75F, 0.75F, 1.0F};
        vv[0] = dynamic_strength;
        vv[1] = dynamic_strength;
        vv[2] = dynamic_strength;
        glClearAccum(0.5, 0.5, 0.5, 0.5);
        glClear(GL_ACCUM_BUFFER_BIT);
        glMaterialfv(GL_FRONT_AND_BACK, GL_EMISSION, vv);
        glDisable(GL_FOG);
        SceneSetViewport(G, I->rect.left + G->Option->winX / 2, I->rect.bottom,
            I->Width, I->Height);
      } else {
        glClearAccum(0.0, 0.0, 0.0, 0.0);
        glClear(GL_ACCUM_BUFFER_BIT);
        SceneSetViewport(G, I->rect.left, I->rect.bottom, I->Width, I->Height);
      }
    } else {
      GLenum err;
      if (times) {
        glAccum(GL_ACCUM, -0.5);
      } else {
        glAccum(GL_ACCUM, 0.5);
        glEnable(GL_SCISSOR_TEST);
      }
      if ((err = glGetError())) {
        if (SettingGetGlobal_i(G, cSetting_stereo_mode) != cStereo_crosseye) {
          PRINTFB(G, FB_Scene, FB_Errors)
          "Stereo Error 0x%x: stereo_mode=11 dynamic requires access to the "
          "accumulation buffer,\n"
          "you need to start PyMOL with the -t argument, setting back to "
          "default\n",
              err ENDFB(G);
          SettingSetGlobal_i(G, cSetting_stereo_mode, cStereo_crosseye);
          SceneSetStereo(G, 0);
        }
        return;
      }
      glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
      if (!times) {
        glDisable(GL_SCISSOR_TEST);
      }
    }
    break;
#ifdef _PYMOL_OPENVR
  case cStereo_openvr:
    OpenVREyeStart(G, position);
    break;
#endif
#endif
  }
}

/* SetDrawBufferForStereo : called after 2nd/right eye rendered in stereo to
   reset GL state properly based on what was changed in
   PrepareViewPortForStereoImpl per stereo_mode
 */
void SetDrawBufferForStereo(
    PyMOLGlobals* G, CScene* I, int stereo_mode, int times, int fog_active)
{
  switch (stereo_mode) {
  case cStereo_quadbuffer:
    OrthoDrawBuffer(G, GL_BACK_LEFT); /* leave us in a stereo context
                                         (avoids problems with cards than can't
                                         handle use of mono contexts) */
    break;
  case cStereo_crosseye:
  case cStereo_walleye:
  case cStereo_sidebyside:
  case cStereo_openvr:
    OrthoDrawBuffer(G, GL_BACK);
    break;
  case cStereo_geowall:
    break;
  case cStereo_stencil_by_row:
  case cStereo_stencil_by_column:
  case cStereo_stencil_checkerboard:
    glDisable(GL_STENCIL_TEST);
    break;
  case cStereo_stencil_custom:
    break;
  case cStereo_anaglyph:
    glColorMask(true, true, true, true);
    break;
  case cStereo_clone_dynamic:
#ifndef PURE_OPENGL_ES_2
    glAccum(GL_ACCUM, 0.5);
#endif
    if (times) {
      float vv[4] = {0.0F, 0.0F, 0.0F, 0.0F};
#ifndef PURE_OPENGL_ES_2
      glMaterialfv(GL_FRONT_AND_BACK, GL_EMISSION, vv);
      if (fog_active)
        glEnable(GL_FOG);
#endif
      glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
      OrthoDrawBuffer(G, GL_BACK_RIGHT);
    }
#ifndef PURE_OPENGL_ES_2
    glAccum(GL_RETURN, 1.0);
#endif
    OrthoDrawBuffer(G, GL_BACK_LEFT);
    break;
  case cStereo_dynamic:
#ifndef PURE_OPENGL_ES_2
    glAccum(GL_ACCUM, 0.5);
#endif
    if (times) {
      float vv[4] = {0.0F, 0.0F, 0.0F, 0.0F};
#ifndef PURE_OPENGL_ES_2
      glMaterialfv(GL_FRONT_AND_BACK, GL_EMISSION, vv);
      if (fog_active)
        glEnable(GL_FOG);
#endif
      glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    }
#ifndef PURE_OPENGL_ES_2
    glAccum(GL_RETURN, 1.0);
#endif
    if (times) {
      SceneSetViewport(
          G, I->rect.left, I->rect.bottom, I->Width + 2, I->Height + 2);
      glScissor(
          I->rect.left - 1, I->rect.bottom - 1, I->Width + 2, I->Height + 2);
      glEnable(GL_SCISSOR_TEST);
      glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
      glDisable(GL_SCISSOR_TEST);
    } else {
      glDisable(GL_SCISSOR_TEST);
    }
    break;
  }
}

void SceneInitializeViewport(PyMOLGlobals* G, bool offscreen)
{
  CScene* I = G->Scene;
  if (offscreen)
    SceneSetViewport(G, 0, 0, I->Width, I->Height);
  else {
    if (I->vp_prepareViewPortForStereo) {
      GLint currentFrameBuffer;
      glGetIntegerv(GL_FRAMEBUFFER_BINDING, &currentFrameBuffer);
      if (currentFrameBuffer ==
          G->ShaderMgr->default_framebuffer_id) { // if writing to screen, then
                                                  // set viewport to screen
        float width_scale;
        // this is called before preparing view port, since the prepare function
        // doesn't setup/change the viewport in all modes
        InitializeViewPortToScreenBlock(
            G, I, I->vp_pos, I->vp_oversize, &I->vp_stereo_mode, &width_scale);
      }
      I->vp_prepareViewPortForStereo(
          G, I, I->vp_stereo_mode, 0, I->vp_times, I->vp_pos, I->vp_oversize);
    } else {
      PRINTFB(G, FB_Scene, FB_Errors)
      " SceneInitializeViewport: I->vp_prepareViewPortForStereo=NULL\n" ENDFB(
          G);
    }
  }
}

void SceneDrawStencilInBuffer(PyMOLGlobals* G, CScene* I, int stereo_mode)
{
  GLint viewport[4];
#ifdef _PYMOL_IOS
  {
    int width, height;
    SceneGetWidthHeight(G, &width, &height);
    viewport[0] = viewport[1] = 0;
    viewport[2] = width;
    viewport[3] = height;
  }
#else
  glGetIntegerv(GL_VIEWPORT, viewport);
#endif

#ifndef PURE_OPENGL_ES_2
  glMatrixMode(GL_PROJECTION);
  glPushMatrix();
  glLoadIdentity();
  glOrtho(0, viewport[2], 0, viewport[3], -10.0, 10.0);
  glMatrixMode(GL_MODELVIEW);
  glPushMatrix();
  glLoadIdentity();
  glTranslatef(0.33F, 0.33F, 0.0F);

  glDisable(GL_ALPHA_TEST);
  glDisable(GL_LIGHTING);
  glDisable(GL_FOG);
  glDisable(GL_NORMALIZE);
  glDisable(GL_COLOR_MATERIAL);
  glDisable(GL_LINE_SMOOTH);
  glShadeModel(
      SettingGetGlobal_b(G, cSetting_pick_shading) ? GL_FLAT : GL_SMOOTH);
  glDisable(0x809D); /* GL_MULTISAMPLE_ARB */
#endif

  glDisable(GL_DEPTH_TEST);
  glDisable(GL_DITHER);
  glDisable(GL_BLEND);

  glDisable(GL_STENCIL_TEST);
  glClearStencil(0);
  glColorMask(false, false, false, false);
  glDepthMask(false);
  glClear(GL_STENCIL_BUFFER_BIT);

  glEnable(GL_STENCIL_TEST);
  glStencilFunc(GL_ALWAYS, 1, 1);
  glStencilOp(GL_KEEP, GL_KEEP, GL_REPLACE);

#ifndef PURE_OPENGL_ES_2
  {
    int h = viewport[3], w = viewport[2];
    glLineWidth(1.0);
    switch (stereo_mode) {
    case cStereo_stencil_by_row: {
      int parity = I->StencilParity;
      int y;
      glBegin(GL_LINES);
      for (y = 0; y < h; y += 2) {
        glVertex2i(0, y + parity);
        glVertex2i(w, y + parity);
      }
      glEnd();
    } break;
    case cStereo_stencil_by_column: {
      int x;
      glBegin(GL_LINES);
      for (x = 0; x < w; x += 2) {
        glVertex2i(x, 0);
        glVertex2i(x, h);
      }
      glEnd();
    } break;
    case cStereo_stencil_checkerboard: {
      int i, m = 2 * ((h > w) ? h : w);
      glBegin(GL_LINES);
      for (i = 0; i < m; i += 2) {
        glVertex2i(i, 0);
        glVertex2i(0, i);
      }
      glEnd();
    } break;
    }
  }
#endif

  glColorMask(true, true, true, true);
  glDepthMask(true);

#ifndef PURE_OPENGL_ES_2
  glMatrixMode(GL_MODELVIEW);
  glPopMatrix();
  glMatrixMode(GL_PROJECTION);
  glPopMatrix();
#endif
}

CGO* GenerateUnitScreenCGO(PyMOLGlobals* G)
{
  CGO cgo(G);
  CGOBegin(&cgo, GL_TRIANGLE_STRIP);
  CGOVertex(&cgo, -1.f, -1.f, 0.98f);
  CGOVertex(&cgo, 1.f, -1.f, 0.98f);
  CGOVertex(&cgo, -1.f, 1.f, 0.98f);
  CGOVertex(&cgo, 1.f, 1.f, 0.98f);
  CGOEnd(&cgo);
  assert(cgo.has_begin_end);
  return CGOOptimizeToVBONotIndexed(&cgo, 0);
}

static void SceneRenderPostProcessStack(PyMOLGlobals* G)
{
  /**
   * TODO: Postprocess rendering will eventually reside here.
   * For now we'll keep this generalized rendering function to
   * decouple antialiasing logic from framebuffer management.
   */
  SceneRenderAA(G);
}
