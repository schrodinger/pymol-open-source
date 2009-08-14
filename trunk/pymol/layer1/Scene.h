
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

typedef struct {
  unsigned char *data;
  int size;
  int width, height;
  int stereo;                   /* indicates data actually contains two back to back full-screen images */
  int needs_alpha_reset;        /* needs alpha reset */
} ImageType;

#define cSceneViewSize 25

#define SDOF_NORMAL_MODE 0
#define SDOF_CLIP_MODE 1
#define SDOF_DRAG_MODE 2

#define cSceneRay_MODE_IDTF 7

#define cSceneImage_Default -1
#define cSceneImage_Normal 0
#define cSceneImage_Draw 1
#define cSceneImage_Ray 2

typedef float SceneViewType[cSceneViewSize];


/* all information required to define the geometry of a particular view,
   for shipping to and from python as a list of floats
   0-15 = 4x4 rotation matrix 
   16-18 = position
   19-21 = origin
   22    = front plane
   23    = rear plane
   24    = orthoscopic flag 
*/

float SceneGetDynamicLineWidth(RenderInfo * info, float line_width);
void SceneInvalidateStencil(PyMOLGlobals * G);
int SceneHasImage(PyMOLGlobals * G);
int SceneInit(PyMOLGlobals * G);
void SceneDone(PyMOLGlobals * G);
void SceneUpdate(PyMOLGlobals * G, int force);
int SceneRenderCached(PyMOLGlobals * G);
void SceneRender(PyMOLGlobals * G, Picking * pick, int x, int y,
                 Multipick * smp, int oversize_width, int oversize_height,
                 int click_side, int force_copy);
void SceneSetFrame(PyMOLGlobals * G, int mode, int frame);
int SceneSetNames(PyMOLGlobals * G, PyObject * list);
int SceneGetFrame(PyMOLGlobals * G);
int SceneGetState(PyMOLGlobals * G);
int SceneGetButtonMargin(PyMOLGlobals * G);

void SceneDirty(PyMOLGlobals * G);      /* scene dirty, but leave the overlay if one exists */
void SceneInvalidate(PyMOLGlobals * G); /* scene dirty and remove the overlay */
void SceneChanged(PyMOLGlobals * G);    /* update 3D objects */

void SceneCountFrames(PyMOLGlobals * G);
int SceneGetNFrame(PyMOLGlobals * G, int *has_movie);
void SceneSetMatrix(PyMOLGlobals * G, float *);
float *SceneGetMatrix(PyMOLGlobals * G);

void SceneReshape(Block * block, int width, int height);
float SceneGetScreenVertexScale(PyMOLGlobals * G, float *v1);

void SceneTest(PyMOLGlobals * G);
void SceneIdle(PyMOLGlobals * G);
void SceneFree(PyMOLGlobals * G);
void SceneRay(PyMOLGlobals * G, int width, int height, int mode,
              char **headerVLA, char **charVLA,
              float angle, float shift, int quiet,
              G3dPrimitive ** g3d, int show_timing, int antialias);
void SceneDoRay(PyMOLGlobals * G, int width, int height, int mode,
                char **headerVLA, char **charVLA,
                float angle, float shift, int quiet,
                G3dPrimitive ** g3d, int show_timing, int antialias);
int SceneDeferRay(PyMOLGlobals * G,
                  int ray_width,
                  int ray_height,
                  int mode,
                  float angle, float shift, int quiet, int show_timing, int antialias);
int SceneMakeMovieImage(PyMOLGlobals * G, int show_timing, int validate, int mode);

int ScenePNG(PyMOLGlobals * G, char *png, float dpi, int quiet,
             int prior_only, int format);
int SceneCopyExternal(PyMOLGlobals * G, int width, int height, int rowbytes,
                      unsigned char *dest, int mode);

void SceneResetMatrix(PyMOLGlobals * G);

void SceneRestartFrameTimer(PyMOLGlobals * G);

void ScenePerspective(PyMOLGlobals * G, int flag);
void SceneGetEyeNormal(PyMOLGlobals * G, float *v1, float *normal);

void SceneRotate(PyMOLGlobals * G, float angle, float x, float y, float z);
void SceneTranslate(PyMOLGlobals * G, float x, float y, float z);

void SceneTranslateScaled(PyMOLGlobals * G, float x, float y, float z, int sdof_mode);
void SceneRotateScaled(PyMOLGlobals * G, float rx, float ry, float rz, int sdof_mode);

void SceneClip(PyMOLGlobals * G, int plane, float movement, char *sele, int state);
void SceneGetImageSize(PyMOLGlobals * G, int *width, int *height);

void SceneScale(PyMOLGlobals * G, float scale);
void SceneResetNormal(PyMOLGlobals * G, int lines);

int SceneObjectAdd(PyMOLGlobals * G, CObject * obj);
int SceneObjectDel(PyMOLGlobals * G, CObject * obj, int allow_purge);
int SceneObjectIsActive(PyMOLGlobals * G, CObject * obj);
void SceneOriginSet(PyMOLGlobals * G, float *origin, int preserve);
void SceneOriginGet(PyMOLGlobals * G, float *origin);
void SceneWindowSphere(PyMOLGlobals * G, float *location, float radius);
void SceneRelocate(PyMOLGlobals * G, float *location);
Block *SceneGetBlock(PyMOLGlobals * G);
void SceneApplyMatrix(PyMOLGlobals * G, float *m);
void SceneSetStereo(PyMOLGlobals * G, int flag);
int SceneGetStereo(PyMOLGlobals * G);
void SceneDontCopyNext(PyMOLGlobals * G);
void ScenePrepareExit(PyMOLGlobals * G);
void SceneGetViewNormal(PyMOLGlobals * G, float *v);
void SceneClipSet(PyMOLGlobals * G, float front, float back);
void SceneGetView(PyMOLGlobals * G, SceneViewType view);
void SceneSetView(PyMOLGlobals * G, SceneViewType view,
                  int quiet, float animate, int hand);
void SceneRestartSweepTimer(PyMOLGlobals * G);
int SceneViewEqual(SceneViewType left, SceneViewType right);
void SceneToViewElem(PyMOLGlobals * G, CViewElem * elem, char *scene_name);
void SceneFromViewElem(PyMOLGlobals * G, CViewElem * elem, int dirty);
void SceneGetPos(PyMOLGlobals * G, float *pos);
void SceneGetWidthHeight(PyMOLGlobals * G, int *width, int *height);
int SceneMultipick(PyMOLGlobals * G, Multipick * smp);
void SceneInvalidateCopy(PyMOLGlobals * G, int free_buffer);

void SceneSetCardInfo(PyMOLGlobals * G, char *vendor, char *renderer, char *version);
void SceneGetCardInfo(PyMOLGlobals * G, char **vendor, char **renderer, char **version);
int SceneLoadPNG(PyMOLGlobals * G, char *fname, int movie_flag, int stereo, int quiet);

void SceneSetDefaultView(PyMOLGlobals * G);
void SceneApplyRotMatrix(PyMOLGlobals * G, float *src, float *dst);
void SceneRovingDirty(PyMOLGlobals * G);
int SceneRovingCheckDirty(PyMOLGlobals * G);
void SceneRovingUpdate(PyMOLGlobals * G);
void SceneRovingChanged(PyMOLGlobals * G);
void SceneRovingPostpone(PyMOLGlobals * G);
void SceneCleanupStereo(PyMOLGlobals * G);
int SceneReinitialize(PyMOLGlobals * G);
void SceneUpdateStereoMode(PyMOLGlobals * G);
void SceneSuppressMovieFrame(PyMOLGlobals * G);
int SceneDeferClick(Block * block, int button, int x, int y, int mod);
int SceneDeferRelease(Block * block, int button, int x, int y, int mod);
int SceneDeferDrag(Block * block, int x, int y, int mod);
int SceneDeferImage(PyMOLGlobals * G, int width, int height, char *filename,
                    int antialias, float dpi, int format, int quiet);
char *SceneGetSeleModeKeyword(PyMOLGlobals * G);
void SceneUpdateStereo(PyMOLGlobals * G);
void ScenePushRasterMatrix(PyMOLGlobals * G, float *v);
void ScenePopRasterMatrix(PyMOLGlobals * G);
void ScenePrimeAnimation(PyMOLGlobals * G);
void SceneLoadAnimation(PyMOLGlobals * G, double duration, int hand);
int SceneMustDrawBoth(PyMOLGlobals * G);
float SceneGetReflectScaleValue(PyMOLGlobals * G, int limit);
float SceneGetSpecularValue(PyMOLGlobals * G, float spec, int limit);
void SceneAbortAnimation(PyMOLGlobals * G);
void SceneObjectUpdateThread(CObjectUpdateThreadInfo * T);
int SceneCaptureWindow(PyMOLGlobals * G);

#endif
