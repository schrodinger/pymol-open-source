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

#define cSceneViewSize 25
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

int SceneInit(PyMOLGlobals *G);
void SceneDone(PyMOLGlobals *G);
void SceneUpdate(PyMOLGlobals *G);
int SceneRenderCached(PyMOLGlobals *G);
void SceneRender(PyMOLGlobals *G,Pickable *pick,int x,int y,Multipick *smp);
void SceneSetFrame(PyMOLGlobals *G,int mode,int frame);
int SceneGetFrame(PyMOLGlobals *G);
int SceneGetState(PyMOLGlobals *G);
void SceneDirty(PyMOLGlobals *G); /* must update transformation */
void SceneChanged(PyMOLGlobals *G); /* must actually update 3D objects */
void SceneCountFrames(PyMOLGlobals *G) ;
int SceneGetNFrame(PyMOLGlobals *G);
void SceneSetMatrix(PyMOLGlobals *G,float *);
float *SceneGetMatrix(PyMOLGlobals *G);

void SceneReshape(Block *block,int width,int height);
float SceneGetScreenVertexScale(PyMOLGlobals *G,float *v1);

void SceneTest(PyMOLGlobals *G);
void SceneIdle(PyMOLGlobals *G);
void SceneFree(PyMOLGlobals *G);
void SceneRay(PyMOLGlobals *G,int width,int height,int mode,
              char **headerVLA,char **charVLA,
              float angle,float shift,int quiet);
void SceneMakeMovieImage(PyMOLGlobals *G);

void ScenePNG(PyMOLGlobals *G,char *png,int quiet);
int SceneCopyExternal(PyMOLGlobals *G,int width, int height,int rowbytes,unsigned char *dest);

void SceneResetMatrix(PyMOLGlobals *G);

void SceneRestartTimers(PyMOLGlobals *G);

void ScenePerspective(PyMOLGlobals *G,int flag);

void SceneRotate(PyMOLGlobals *G,float angle,float x,float y,float z);
void SceneTranslate(PyMOLGlobals *G,float x,float y, float z);
void SceneClip(PyMOLGlobals *G,int plane,float movement,char *sele,int state);

void SceneScale(PyMOLGlobals *G,float scale);
void SceneResetNormal(PyMOLGlobals *G,int lines);

void SceneObjectAdd(PyMOLGlobals *G,CObject *obj);
void SceneObjectDel(PyMOLGlobals *G,CObject *obj);
void SceneOriginSet(PyMOLGlobals *G,float *origin,int preserve);
void SceneOriginGet(PyMOLGlobals *G,float *origin);
void SceneWindowSphere(PyMOLGlobals *G,float *location,float radius);
void SceneRelocate(PyMOLGlobals *G,float *location);
Block *SceneGetBlock(PyMOLGlobals *G);
void SceneApplyMatrix(PyMOLGlobals *G,float *m);
void SceneSetStereo(PyMOLGlobals *G,int flag);
int SceneGetStereo(PyMOLGlobals *G);
void ScenePurgeCopy(PyMOLGlobals *G);
void SceneDontCopyNext(PyMOLGlobals *G);
void ScenePrepareExit(PyMOLGlobals *G);
void SceneGetViewNormal(PyMOLGlobals *G,float *v);
void SceneClipSet(PyMOLGlobals *G,float front,float back);
void SceneGetView(PyMOLGlobals *G,SceneViewType view);
void SceneSetView(PyMOLGlobals *G,SceneViewType view,int quiet);

void SceneToViewElem(PyMOLGlobals *G,CViewElem *elem);
void SceneFromViewElem(PyMOLGlobals *G,CViewElem *elem);
void SceneGetPos(PyMOLGlobals *G,float *pos);
void SceneGetWidthHeight(PyMOLGlobals *G,int *width,int *height);
int SceneMultipick(PyMOLGlobals *G,Multipick *smp);

void SceneSetCardInfo(PyMOLGlobals *G,char *vendor,char *renderer,char *version);
void SceneGetCardInfo(PyMOLGlobals *G,char **vendor,char **renderer,char **version);
int SceneLoadPNG(PyMOLGlobals *G,char *fname,int movie_flag,int quiet);

void SceneSetDefaultView(PyMOLGlobals *G);
void SceneApplyRotMatrix(PyMOLGlobals *G,float *src,float *dst);
void SceneRovingDirty(PyMOLGlobals *G);
int SceneRovingCheckDirty(PyMOLGlobals *G);
void SceneRovingUpdate(PyMOLGlobals *G);
void SceneRovingChanged(PyMOLGlobals *G);
void SceneRovingPostpone(PyMOLGlobals *G);
void SceneCleanupStereo(PyMOLGlobals *G);
int SceneReinitialize(PyMOLGlobals *G);
void SceneUpdateStereoMode(PyMOLGlobals *G);
void SceneSuppressMovieFrame(PyMOLGlobals *G);
int SceneClick(Block *block,int button,int x,int y,int mod);
int SceneRelease(Block *block,int button,int x,int y,int mod);
int SceneDrag(Block *block,int x,int y,int mod);
char *SceneGetSeleModeKeyword(PyMOLGlobals *G);
void SceneUpdateStereo(PyMOLGlobals *G);
void ScenePushRasterMatrix(PyMOLGlobals *G,float *v);
void ScenePopRasterMatrix(PyMOLGlobals *G);

#endif



