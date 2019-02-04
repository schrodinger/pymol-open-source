

#include"Base.h"
#include"Scene.h"
#include"ScenePicking.h"
#include"SceneRender.h"
#include"ShaderMgr.h"
#include"MemoryDebug.h"
#include"PyMOL.h"
#include"P.h"
#include "MacPyMOL.h"
#include "Err.h"

typedef unsigned char pix[4];
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
  I->invPick = false;
  return (I->LastPicked.context.object != NULL);
  /* did we pick something? */
}

static
unsigned int SceneFindTriplet(PyMOLGlobals * G, int x, int y, GLenum gl_buffer, bool bits32)
{
  unsigned int result = 0;
  /*int before_check[100];
     int *int_ptr;
   */
  pix *buffer = NULL;
  pix *extra_safe_buffer = NULL;

  /*int after_check[100]; */
  /* pix_array *array_ptr;
     char *safe_place;
   */
  int a, b, d, flag;
  float contentScaleFactor = DIP2PIXEL(1);
  int cRangeVal = contentScaleFactor < 1.5 ? 7 : 21;

  int h = (cRangeVal * 2 + 1), w = (cRangeVal * 2 + 1);

  int debug = false;
  unsigned char *c;
  int bits15 = false;
  GLint rb, gb, bb, ab;
  int bkrd_alpha = 0xFF;
  int check_alpha = false;

  if(G->HaveGUI && G->ValidContext) {   /*just in case */

    glGetIntegerv(GL_RED_BITS, &rb);
    glGetIntegerv(GL_GREEN_BITS, &gb);
    glGetIntegerv(GL_BLUE_BITS, &bb);
    glGetIntegerv(GL_ALPHA_BITS, &ab);

    bits15 = (rb == 5) && (gb == 5) && (bb == 5);
    if((rb < 4) && (gb < 4) && (bb < 4)){
      PRINTFB(G, FB_Scene, FB_Errors) "SceneFindTriplet: ERROR: not enough colors to pick: rb=%d gb=%d bb=%d\n", rb, gb, bb ENDFB(G);
      return 0;
    }
    if(Feedback(G, FB_Scene, FB_Debugging))
      debug = true;

#ifndef PURE_OPENGL_ES_2
    if (!hasFrameBufferBinding())
      glReadBuffer(gl_buffer);
#endif
    extra_safe_buffer = pymol::malloc<pix>(w * h * 21);
    buffer = extra_safe_buffer + (w * h * 10);

    PyMOLReadPixels(x - cRangeVal, y - cRangeVal, cRangeVal * 2 + 1, cRangeVal * 2 + 1, GL_RGBA,
                    GL_UNSIGNED_BYTE, &buffer[0][0]);

    if(debug) {
      for(a = 0; a <= (cRangeVal * 2); a++) {
        for(b = 0; b <= (cRangeVal * 2); b++)
          printf("%2x ",
                 (buffer[a + b * w][0] + buffer[a + b * w][1] +
                  buffer[a + b * w][2]) & 0xFF);
        printf("\n");
      }
      printf("\n");
      for(a = 0; a <= (cRangeVal * 2); a++) {
        for(b = 0; b <= (cRangeVal * 2); b++)
          printf("%02x ", (buffer[a + b * w][3]) & 0xFF);
        printf("\n");
      }
      printf("\n");
      for(a = 0; a <= (cRangeVal * 2); a++) {
        for(b = 0; b <= (cRangeVal * 2); b++)
          printf("%02x%02x%02x ", (buffer[a + b * w][0]) & 0xFF,
                 (buffer[a + b * w][1]) & 0xFF, (buffer[a + b * w][2]) & 0xFF);
        printf("\n");
      }
      printf("\n");
    }

    /* first, check to make sure bkrd_alpha is correct 
       (this is a bug for systems with broken alpha, such as Extreme 3D on Solaris 8 */

    if (bits32){
      check_alpha = false;
    } else {
      flag = true;
      for(d = 0; ab && flag && (d < cRangeVal); d++)
        for(a = -d; flag && (a <= d); a++)
          for(b = -d; flag && (b <= d); b++) {
            c = &buffer[(a + cRangeVal) + (b + cRangeVal) * w][0];
            if(c[3] == bkrd_alpha) {
              check_alpha = true;
              flag = false;
            }
          }
    }

    /* now find the correct pixel */
    flag = true;
    bool check_green_bit = !bits32;
    bool strict = true;
    if (bits32 || bits15)
      strict = false;
    for(d = 0; flag && (d < cRangeVal); d++)
      for(a = -d; flag && (a <= d); a++)
        for(b = -d; flag && (b <= d); b++) {
          c = &buffer[(a + cRangeVal) + (b + cRangeVal) * w][0];
          if(
             ((c[3] == bkrd_alpha) || (!check_alpha)) &&
             ((bits15 && c[1]) || (c[1] & 0x8) || !check_green_bit) &&
             ((!strict) || ((c[1] & 0xF) == 8 && (c[0] & 0xF) == 0 && (c[2] & 0xF) == 0))
             ) {           /* only consider intact, saturated pixels */
	    if (bits15){  /* workaround for 15 bit rendering, for some reason red/green need rounding */
	      c[0] += 0x8;
	      c[2] += 0x8;
	    }
            if (bits32){
              result = (c[0] & 0xFF) + ((c[1] & 0xFF) << 8) + ((c[2] & 0xFF) << 16) + ((c[3] & 0xFF) << 24);
	      if (result)
		flag = false;
	    }
            else {
              result = ((c[0] >> 4) & 0xF) + (c[1] & 0xF0) + ((c[2] << 4) & 0xF00);
	      flag = false;
	    }
          }
        }
    FreeP(extra_safe_buffer);
  }
  return (result);
}

bool SceneHas32BitColor(PyMOLGlobals * G) {
#ifndef PURE_OPENGL_ES_2
  GLint bits;
  GLint currentFrameBuffer;

  ok_assert(1, SettingGetGlobal_b(G, cSetting_use_shaders));
  ok_assert(1, SettingGetGlobal_b(G, cSetting_pick32bit));

  glGetIntegerv(GL_FRAMEBUFFER_BINDING, &currentFrameBuffer);

  if (currentFrameBuffer != G->ShaderMgr->default_framebuffer_id) {
    glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, G->ShaderMgr->default_framebuffer_id);
  }

  glGetIntegerv(GL_ALPHA_BITS,  &bits); ok_assert(2, bits >= 8);
  glGetIntegerv(GL_BLUE_BITS,   &bits); ok_assert(2, bits >= 8);
  glGetIntegerv(GL_GREEN_BITS,  &bits); ok_assert(2, bits >= 8);
  glGetIntegerv(GL_RED_BITS,    &bits);

ok_except2:
  if (currentFrameBuffer != G->ShaderMgr->default_framebuffer_id) {
    glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, currentFrameBuffer);
  }

  ok_assert(1, bits >= 8);

  PRINTFD(G, FB_Scene) "Scene-DEBUG: 32bit picking\n" ENDFD;
  return true;

ok_except1:
  PRINTFD(G, FB_Scene) "Scene-DEBUG: 16bit picking\n" ENDFD;
#endif
  return false;
}

static
void SceneRenderPickingSinglePick(PyMOLGlobals * G, SceneUnitContext *context, Picking * pick, int x, int y, GLenum render_buffer){
  /* atom picking HACK - obfuscative coding */
  CScene *I = G->Scene;
  int debug_pick = 0;
  unsigned int lowBits = 0, highBits = 0;
  unsigned int index;
  bool bits32 = SceneHas32BitColor(G);

  debug_pick = SettingGetGlobal_i(G, cSetting_debug_pick);
  SceneGLClearColor(0.0, 0.0, 0.0, 0.);

  if (I->pickVLA.empty()){
    I->pickVLA.resize(5000);
  }

  if(I->grid.active)
    GridGetGLViewport(G, &I->grid);

  // two passes in 16bit picking mode, if needed
  for (int pass = 0;; ++pass) {
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    if (I->invPick || !SettingGetGlobal_b(G, cSetting_use_shaders)){
      I->pickVLA.begin()->src.index = 0;
      I->pickVLA.begin()->src.bond = 2 + pass;
    } else {
      I->pickVLA.begin()->src.bond = 0 + pass;
    }
    {
      int slot;
      for(slot = 0; slot <= I->grid.last_slot; slot++) {
        if(I->grid.active) {
          GridSetGLViewport(&I->grid, slot);
        }
        SceneRenderAll(G, context, NULL, std::addressof(I->pickVLA), 0, true, 0.0F, &I->grid, 0, 0, bits32);
      }
    }

    if(debug_pick) {
      PyMOL_SwapBuffers(G->PyMOL);
      PSleep(G, 1000000 * debug_pick / 4);
      PyMOL_SwapBuffers(G->PyMOL);
    }

    if (pass == 1) {
      highBits = SceneFindTriplet(G, x, y, render_buffer, false);
      index += (highBits << 12);
      break;
    }

    index = lowBits = SceneFindTriplet(G, x, y, render_buffer, bits32);

    if (bits32 || I->pickVLA[0].src.index < (1 << 12)) {
      // no need for a second pass
      break;
    }
  }

  if(I->grid.active)
    GridSetGLViewport(&I->grid, -1);

  if(debug_pick) {
    if (bits32){
      PRINTFB(G, FB_Scene, FB_Details)
	" SceneClick-Detail: lowBits=%u index %u < %u?\n", lowBits, index, I->pickVLA.begin()->src.index ENDFB(G);
    } else {
      PRINTFB(G, FB_Scene, FB_Details)
	" SceneClick-Detail: lowBits=%u highBits=%u index %u < %u?\n", lowBits, highBits, index, I->pickVLA.begin()->src.index ENDFB(G);
    }
  }

  if(index && (index <= I->pickVLA.begin()->src.index)) {
    *pick = I->pickVLA[index];       /* return object info */
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

/*========================================================================*/
static
unsigned int *SceneReadTriplets(PyMOLGlobals * G, int x, int y, int w, int h,
                                GLenum gl_buffer, bool bits32)
{
  unsigned int *result = NULL;
  pix *buffer = NULL;
  pix *extra_safe_buffer = NULL;
  int a, b;
  unsigned char *c;
  int cc = 0;
  int bits15 = false;
  int bkrd_alpha = 0xFF;
  int check_alpha = false;

  GLint rb, gb, bb, ab;

  if(w < 1)
    w = 1;
  if(h < 1)
    h = 1;
  if(G->HaveGUI && G->ValidContext) {   /*just in case */
    glGetIntegerv(GL_RED_BITS, &rb);
    glGetIntegerv(GL_GREEN_BITS, &gb);
    glGetIntegerv(GL_BLUE_BITS, &bb);
    glGetIntegerv(GL_ALPHA_BITS, &ab);

    bits15 = (rb == 5) && (gb == 5) && (bb == 5);
    if((rb < 4) && (gb < 4) && (bb < 4)){
      PRINTFB(G, FB_Scene, FB_Errors) "SceneReadTriplet: ERROR: not enough colors to pick: rb=%d gb=%d bb=%d\n", rb, gb, bb ENDFB(G);
      return 0;
    }
    /* create some safe RAM on either side of the read buffer -- buggy
       ReadPixels implementations tend to trash RAM surrounding the
       target block */

    extra_safe_buffer = pymol::malloc<pix>(w * h * 11);
    buffer = extra_safe_buffer + (w * h * 5);

    result = VLAlloc(unsigned int, w * h);
#ifndef PURE_OPENGL_ES_2
    if (!hasFrameBufferBinding())
      glReadBuffer(gl_buffer);
#endif
    PyMOLReadPixels(x, y, w, h, GL_RGBA, GL_UNSIGNED_BYTE, &buffer[0][0]);

    /* first, check to make sure bkrd_alpha is correct 
       (this is a bug for systems with broken alpha, such as Extreme 3D on Solaris 8 */

    if (bits32)
      check_alpha = false;
    else
      for(a = 0; ab && a < w; a++)
	for(b = 0; b < h; b++) {
	  c = &buffer[a + b * w][0];
	  if(c[3] == bkrd_alpha) {
	    check_alpha = true;
	  }
        }

    bool check_green_bit = !bits32;
    bool strict = true;
    if (bits32 || bits15)
      strict = false;
    /* now read pixels */
    for(a = 0; a < w; a++)
      for(b = 0; b < h; b++) {
        c = &buffer[a + b * w][0];
        if(
           ((c[3] == bkrd_alpha) || (!check_alpha)) &&
           ((bits15 && c[1]) || (c[1] & 0x8) || !check_green_bit) &&
           ((!strict) || ((c[1] & 0xF) == 8 && (c[0] & 0xF) == 0 && (c[2] & 0xF) == 0))
           ) {             /* only consider intact, saturated pixels */
          VLACheck(result, unsigned int, cc + 1);
	  if (bits15){  /* workaround for 15 bit rendering, for some reason red/green need rounding */
	    c[0] += 0x8;
	    c[2] += 0x8;
	  }
          if (bits32)
            result[cc] = (c[0] & 0xFF) + ((c[1] & 0xFF) << 8) + ((c[2] & 0xFF) << 16) + ((c[3] & 0xFF) << 24);
          else {
            result[cc] = ((c[0] >> 4) & 0xF) + (c[1] & 0xF0) + ((c[2] << 4) & 0xF00);
          }
          result[cc + 1] = b + a * h; // resulting pixel offset
          if (result[cc] || (!bits32))  // !bits32 because it should hit the high-order bits
            cc += 2;
        }
      }
    FreeP(extra_safe_buffer);
    VLASize(result, unsigned int, cc);
  }
  return (result);
}

static
void SceneRenderPickingMultiPick(PyMOLGlobals * G, SceneUnitContext *context, Multipick * smp, GLenum render_buffer){
  /* multiple atom picking HACK - even more obfuscative coding */
  CScene *I = G->Scene;
  Picking *pik;
  unsigned int *lowBitVLA = NULL, *highBitVLA = NULL;
  int high, low;
  unsigned int lastIndex = 0;
  unsigned int index;
  void *lastPtr = NULL;
  int nPick;
  int nHighBits = 0, nLowBits;
  bool bits32 = SceneHas32BitColor(G);

  SceneGLClearColor(0.0, 0.0, 0.0, 0.0);

  if (I->pickVLA.empty()){
    I->pickVLA.resize(5000);
  }

  if(I->grid.active)
    GridGetGLViewport(G, &I->grid);

  // two passes in 16bit picking mode, if needed
  for (int pass = 0;; ++pass) {
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    if (I->invPick || !SettingGetGlobal_b(G, cSetting_use_shaders)){
      I->pickVLA.begin()->src.index = 0;
      I->pickVLA.begin()->src.bond = 2 + pass;
    } else {
      I->pickVLA.begin()->src.bond = 0 + pass;
    }

    {
      int slot;
      for(slot = 0; slot <= I->grid.last_slot; slot++) {
        if(I->grid.active) {
          GridSetGLViewport(&I->grid, slot);
        }
        SceneRenderAll(G, context, NULL, std::addressof(I->pickVLA), 0, true, 0.0F, &I->grid, 0, 0, bits32);
      }
    }

    if (pass == 1) {
      highBitVLA = SceneReadTriplets(G, smp->x, smp->y, smp->w, smp->h, render_buffer, false);
      nHighBits = VLAGetSize(highBitVLA);
      break;
    }

    lowBitVLA = SceneReadTriplets(G, smp->x, smp->y, smp->w, smp->h, render_buffer, bits32);
    nLowBits = VLAGetSize(lowBitVLA);

    if (bits32 || I->pickVLA.begin()->src.index < (1 << 12)) {
      // no need for a second pass
      bits32 = true; // continue like with 32bit picking
      break;
    }
  }

  if(I->grid.active)
    GridSetGLViewport(&I->grid, -1);

  /* need to scissor this */ 
  nPick = 0;
  if(nLowBits && (bits32 || nHighBits)) {
    low = 0;
    high = 0;
    while((low < nLowBits) && (bits32 || high < nHighBits)) {
      if(bits32 || lowBitVLA[low + 1] == highBitVLA[high + 1]) {
        if (bits32){
          // 32bit picking
          index = lowBitVLA[low];
        } else {
          index = lowBitVLA[low] + (highBitVLA[high] << 12);
        }
        if(index && (index <= I->pickVLA.begin()->src.index)) {
          pik = I->pickVLA.data() + index;  /* just using as a tmp */
          if((pik->src.index != lastIndex) || (pik->context.object != lastPtr)) {
            if(((CObject *) pik->context.object)->type == cObjectMolecule) {
              nPick++;    /* start from 1 */
              VLACheck(smp->picked, Picking, nPick);
              smp->picked[nPick] = *pik;  /* return atom/object info -- will be redundant */
            }
            lastIndex = pik->src.index;
            lastPtr = pik->context.object;
          }
        }
        low += 2;
        high += 2;
      } else if(lowBitVLA[low + 1] < highBitVLA[high + 1])
        low += 2;
      else
        high += 2;
    }
  }
  smp->picked[0].src.index = nPick;
#ifndef PURE_OPENGL_ES_2
  /* Picking changes the Shading model to GL_FLAT,
     we need to change it back to GL_SMOOTH. This is because
     bg_grad() might be called in OrthoDoDraw() before GL 
     settings are set in SceneRender() */
  //      glEnable(GL_COLOR_MATERIAL);
  glShadeModel(SettingGetGlobal_b(G, cSetting_pick_shading) ? GL_FLAT : GL_SMOOTH);
#endif
  VLAFreeP(lowBitVLA);
  VLAFreeP(highBitVLA);
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
