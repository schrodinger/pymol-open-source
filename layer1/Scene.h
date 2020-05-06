
/* 
A* -------------------------------------------------------------------
B* This file contains source code for the PyMOL computer program
C* copyright 1998-2000 by Warren Lyford Delano of DeLano Scientific. 
D* -------------------------------------------------------------------
E* It is unlawful to modify or remove this copyright notice.
F* -------------------------------------------------------------------
G* Please see the accompanying LICENSE file for further information. 
H* -------------------------------------------------------------------
I* Additional authors of this source file include:
-* 
-* 
-*
Z* -------------------------------------------------------------------
*/

#ifndef _H_Scene
#define _H_Scene

#include"Base.h"
#include"PyMOLObject.h"
#include"Ortho.h"
#include"View.h"
#include"Result.h"
#include"SceneView.h"
#include"SceneDef.h"
#include"SceneRender.h"
#include"ShaderMgr.h"
#include"pymol/zstring_view.h"

#include <vector>

#define SDOF_NORMAL_MODE 0
#define SDOF_CLIP_MODE 1
#define SDOF_DRAG_MODE 2

// TODO: define remaining cSceneRay_MODEs (VRML, COLLADA, etc.)
#define cSceneRay_MODE_IDTF 7

#define cSceneImage_Default -1
#define cSceneImage_Normal 0
#define cSceneImage_Draw 1
#define cSceneImage_Ray 2

#define SceneScrollBarMargin DIP2PIXEL(1)
#define SceneScrollBarWidth DIP2PIXEL(13)
#define SceneGetExactScreenVertexScale SceneGetScreenVertexScale

enum cSceneClip {
  cSceneClip_invalid = -1,
  cSceneClip_near = 0,
  cSceneClip_far = 1,
  cSceneClip_move = 2,
  cSceneClip_slab = 3,
  cSceneClip_atoms = 4,
};

float SceneGetDynamicLineWidth(RenderInfo * info, float line_width);
void SceneInvalidateStencil(PyMOLGlobals * G);
int SceneHasImage(PyMOLGlobals * G);
int SceneInit(PyMOLGlobals * G);
void SceneUpdate(PyMOLGlobals * G, int force);
int SceneRenderCached(PyMOLGlobals * G);
int SceneSetFog(PyMOLGlobals *G);
void SceneSetFogUniforms(PyMOLGlobals * G, CShaderPrg *);
void SceneSetFrame(PyMOLGlobals * G, int mode, int frame);
int SceneGetFrame(PyMOLGlobals * G);
int SceneGetState(PyMOLGlobals * G);

void SceneDirty(PyMOLGlobals * G);      /* scene dirty, but leave the overlay if one exists */
void SceneInvalidate(PyMOLGlobals * G); /* scene dirty and remove the overlay */
void SceneChanged(PyMOLGlobals * G);    /* update 3D objects */

int SceneCountFrames(PyMOLGlobals * G);
int SceneGetNFrame(PyMOLGlobals * G, int *has_movie=nullptr);
void SceneSetMatrix(PyMOLGlobals * G, float *);
float *SceneGetMatrix(PyMOLGlobals * G);

#define SceneGetModMatrix SceneGetModelViewMatrix

float *SceneGetPmvMatrix(PyMOLGlobals * G);

float SceneGetScreenVertexScale(PyMOLGlobals * G, const float *v1);
bool SceneGetVisible(PyMOLGlobals * G, const float *v1);
float SceneGetDepth(PyMOLGlobals * G, const float *v1);
float SceneGetRawDepth(PyMOLGlobals * G, const float *v1);

void SceneIdle(PyMOLGlobals * G);
void SceneFree(PyMOLGlobals * G);

int SceneGetDrawFlagGrid(PyMOLGlobals * G, GridInfo * grid, int slot);
int SceneDeferRay(PyMOLGlobals * G,
                  int ray_width,
                  int ray_height,
                  int mode,
                  float angle, float shift, int quiet, int show_timing, int antialias);
int SceneMakeMovieImage(PyMOLGlobals * G,
    int show_timing, int validate, int mode,
    int width=0, int height=0);
int SceneValidateImageMode(PyMOLGlobals * G, int mode, bool defaultdraw);

bool ScenePNG(PyMOLGlobals* G, const char* png, float dpi, int quiet,
    int prior_only, int format, std::vector<unsigned char>* outbuf = nullptr);
int SceneCopyExternal(PyMOLGlobals * G, int width, int height, int rowbytes,
                      unsigned char *dest, int mode);

void SceneResetMatrix(PyMOLGlobals * G);

void SceneRestartFrameTimer(PyMOLGlobals * G);

void SceneGetEyeNormal(PyMOLGlobals * G, float *v1, float *normal);

void SceneRotateAxis(PyMOLGlobals * G, float angle, char axis);
void SceneRotate(
    PyMOLGlobals* G, float angle, float x, float y, float z, bool dirty = true);
void SceneTranslate(PyMOLGlobals * G, float x, float y, float z);

void SceneTranslateScaled(PyMOLGlobals * G, float x, float y, float z, int sdof_mode);
void SceneRotateScaled(PyMOLGlobals * G, float rx, float ry, float rz, int sdof_mode);

void SceneClip(PyMOLGlobals * G, int plane, float movement, const char *sele, int state);
pymol::Result<> SceneClipFromMode(PyMOLGlobals* G, pymol::zstring_view mode, float movement,
    pymol::zstring_view sele, int state);
std::pair<int, int> SceneGetImageSize(PyMOLGlobals * G);
float SceneGetGridAspectRatio(PyMOLGlobals * G);
void SceneScale(PyMOLGlobals * G, float scale);
void SceneResetNormalCGO(PyMOLGlobals * G, CGO *cgo, int lines);
void SceneResetNormal(PyMOLGlobals * G, int lines);
void SceneResetNormalUseShader(PyMOLGlobals * G, int lines, short use_shader);
void SceneResetNormalToViewVector(PyMOLGlobals * G, short use_shader);
void SceneResetNormalUseShaderAttribute(PyMOLGlobals * G, int lines, short use_shader, int attr);
void SceneGetResetNormal(PyMOLGlobals * G, float *normal, int lines);

int SceneObjectAdd(PyMOLGlobals * G, pymol::CObject* obj);
int SceneObjectDel(PyMOLGlobals * G, pymol::CObject* obj, int allow_purge);

/**
 * @brief removes object from Scene
 * @param obj Object to be removed
 * @return if Object removed--always true if obj==nullptr
 */

bool SceneObjectRemove(PyMOLGlobals* G, pymol::CObject* obj);
int SceneObjectIsActive(PyMOLGlobals* G, pymol::CObject* obj);
void SceneOriginSet(PyMOLGlobals * G, const float *origin, int preserve);
void SceneOriginGet(PyMOLGlobals * G, float *origin);
void SceneWindowSphere(PyMOLGlobals * G, const float *location, float radius);
void SceneRelocate(PyMOLGlobals * G, const float *location);
Block *SceneGetBlock(PyMOLGlobals * G);
void SceneApplyMatrix(PyMOLGlobals * G, float *m);
void SceneSetStereo(PyMOLGlobals * G, bool flag);
int SceneGetStereo(PyMOLGlobals * G);
void SceneDontCopyNext(PyMOLGlobals * G);
void SceneGetViewNormal(PyMOLGlobals * G, float *v);
void SceneClipSet(PyMOLGlobals * G, float front, float back);
void SceneGetView(PyMOLGlobals * G, SceneViewType view);
void SceneSetView(PyMOLGlobals * G, const SceneViewType view,
                  int quiet, float animate, int hand);
void SceneRestartSweepTimer(PyMOLGlobals * G);
int SceneViewEqual(SceneViewType left, SceneViewType right);
void SceneToViewElem(PyMOLGlobals * G, CViewElem * elem, const char *scene_name);
void SceneFromViewElem(PyMOLGlobals * G, CViewElem * elem, int dirty);
void SceneGetCenter(PyMOLGlobals * G, float *pos);
void SceneGetWidthHeight(PyMOLGlobals * G, int *width, int *height);
void SceneGetWidthHeightStereo(PyMOLGlobals * G, int *width, int *height);
void SceneInvalidateCopy(PyMOLGlobals * G, int free_buffer);

void SceneSetCardInfo(PyMOLGlobals * G, const char *vendor, const char *renderer, const char *version);
void SceneGetCardInfo(PyMOLGlobals * G, char **vendor, char **renderer, char **version);
int SceneLoadPNG(PyMOLGlobals * G, const char *fname, int movie_flag, int stereo, int quiet);

void SceneSetDefaultView(PyMOLGlobals * G);
void SceneRovingDirty(PyMOLGlobals * G);
int SceneRovingCheckDirty(PyMOLGlobals * G);
void SceneRovingUpdate(PyMOLGlobals * G);
void SceneRovingChanged(PyMOLGlobals * G);
void SceneRovingPostpone(PyMOLGlobals * G);
int SceneReinitialize(PyMOLGlobals * G);
void SceneUpdateStereoMode(PyMOLGlobals * G);
void SceneSuppressMovieFrame(PyMOLGlobals * G);
int SceneDeferClick(Block * block, int button, int x, int y, int mod);
int SceneDeferDrag(Block * block, int x, int y, int mod);
int SceneDeferImage(PyMOLGlobals * G, int width, int height, const char *filename,
                    int antialias, float dpi, int format, int quiet);
const char *SceneGetSeleModeKeyword(PyMOLGlobals * G);
void SceneUpdateStereo(PyMOLGlobals * G);
float ScenePushRasterMatrix(PyMOLGlobals * G, float *v);
void ScenePopRasterMatrix(PyMOLGlobals * G);
void ScenePrimeAnimation(PyMOLGlobals * G);
void SceneLoadAnimation(PyMOLGlobals * G, double duration, int hand);
int SceneMustDrawBoth(PyMOLGlobals * G);
float SceneGetReflectScaleValue(PyMOLGlobals * G, int limit = 8);
float SceneGetSpecularValue(PyMOLGlobals * G, float spec, int limit = 8);

void SceneGetAdjustedLightValues(PyMOLGlobals * G,
    float *ptr_spec,
    float *ptr_spec_power,
    float *ptr_spec_direct,
    float *ptr_spec_direct_power,
    int limit = 8);

void SceneAbortAnimation(PyMOLGlobals * G);
void SceneObjectUpdateThread(CObjectUpdateThreadInfo * T);
int SceneCaptureWindow(PyMOLGlobals * G);

void SceneZoom(PyMOLGlobals * G, float scale);

int SceneGetTwoSidedLighting(PyMOLGlobals * G);
int SceneGetTwoSidedLightingSettings(PyMOLGlobals * G, const CSetting *set1, const CSetting *set2);

float SceneGetLineWidthForCylinders(PyMOLGlobals * G, RenderInfo * info, float line_width);
float SceneGetLineWidthForCylindersStatic(PyMOLGlobals * G, RenderInfo * info, float dynamic_line_width_arg, float line_width);

void ScenePushModelViewMatrix(PyMOLGlobals * G);
void ScenePopModelViewMatrix(PyMOLGlobals * G, bool);

float *SceneGetModelViewMatrix(PyMOLGlobals * G);
float *SceneGetProjectionMatrix(PyMOLGlobals * G);
void SceneSetBackgroundColorAlreadySet(PyMOLGlobals * G, int);
int SceneGetBackgroundColorAlreadySet(PyMOLGlobals * G);
void SceneSetDoNotClearBackground(PyMOLGlobals * G, int);
int SceneGetDoNotClearBackground(PyMOLGlobals * G);

void SceneProgramLighting(PyMOLGlobals * G, CShaderPrg * shaderPrg = NULL);
void SceneGLClear(PyMOLGlobals * G, GLbitfield mask);

#ifdef _PYMOL_IOS
void SceneTranslateSceneXYWithScale(PyMOLGlobals * G, float x, float y);
int SceneIsTwisting(PyMOLGlobals * G);
#endif

void SceneUpdateAnimation(PyMOLGlobals * G);

void SceneSetupGLPicking(PyMOLGlobals * G);

int SceneDrawImageOverlay(PyMOLGlobals * G, int override  ORTHOCGOARG);

int SceneIncrementTextureRefreshes(PyMOLGlobals * G);

void SceneResetTextureRefreshes(PyMOLGlobals * G);

void SceneGetScaledAxes(
    PyMOLGlobals* G, pymol::CObject* obj, float* xn, float* yn);
void SceneGetScaledAxesAtPoint(PyMOLGlobals * G, float *pt, float *xn, float *yn);

int SceneGetCopyType(PyMOLGlobals * G);

void SceneGenerateMatrixToAnotherZFromZ(PyMOLGlobals *G, float *convMatrix, float *curpt, float *pt);
void SceneAdjustZtoScreenZ(PyMOLGlobals *G, float *pos, float z);
float SceneGetCurrentBackSafe(PyMOLGlobals *G);
float SceneGetCurrentFrontSafe(PyMOLGlobals *G);
void SceneSetPointToWorldScreenRelative(PyMOLGlobals *G, float *pos, float *screenPt);

int StereoIsAdjacent(PyMOLGlobals * G);

int SceneGetGridSize(PyMOLGlobals * G, int grid_mode);

void GridUpdate(GridInfo * I, float asp_ratio, int mode, int size);
void GridGetRayViewport(GridInfo * I, int width, int height);
void GridSetRayViewport(GridInfo * I, int slot, int *x, int *y, int *width,
                        int *height);
int SceneGetDrawFlag(GridInfo * grid, int *slot_vla, int slot);

void SceneApplyImageGamma(PyMOLGlobals * G, unsigned int *buffer, int width,
                          int height);
void ScenePrepareUnitContext(SceneUnitContext * context, int width, int height);

float GetFovWidth(PyMOLGlobals * G);

void ScenePrepareMatrix(PyMOLGlobals * G, int mode, int stereo_mode = 0);

void SceneCopy(PyMOLGlobals * G, GLenum buffer, int force, int entire_window);

// FIXME use pymol matrices
void SceneGetModel2WorldMatrix(PyMOLGlobals * G, float *matrix);
void SceneSetModel2WorldMatrix(PyMOLGlobals * G, float const *matrix);
float SceneGetScale(PyMOLGlobals * G);

void ScenePickAtomInWorld(PyMOLGlobals * G, int x, int y, float *atomWorldPos);

void SceneInvalidatePicking(PyMOLGlobals * G);

pymol::Image* SceneImagePrepare(PyMOLGlobals * G, bool prior_only);

void SceneDoRoving(PyMOLGlobals * G, float old_front,
                   float old_back, float old_origin,
                   int adjust_flag, int zoom_flag);
void UpdateFrontBackSafe(CScene *I);
int stereo_via_adjacent_array(int stereo_mode);
int SceneDeferredClick(DeferredMouse * dm);
int SceneDeferredDrag(DeferredMouse * dm);
int SceneDeferredRelease(DeferredMouse * dm);

#endif

