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

void SceneInit(void);
void SceneDone(void);
void SceneUpdate(void);
int SceneRenderCached(void);
void SceneRender(Pickable *pick,int x,int y,Multipick *smp);
void SceneSetFrame(int mode,int frame);
int SceneGetFrame(void);
int SceneGetState(void);
void SceneDirty(void); /* must update transformation */
void SceneChanged(void); /* must actually update 3D objects */
void SceneCountFrames(void) ;
int SceneGetNFrame(void);
void SceneSetMatrix(float *);
float *SceneGetMatrix(void);

void SceneReshape(Block *block,int width,int height);

void SceneTest(void);
void SceneIdle(void);
void SceneFree(void);
void SceneRay(int width,int height,int mode,
              char **headerVLA,char **charVLA,
              float angle,float shift,int quiet);
void SceneMakeMovieImage(void);

void ScenePNG(char *png,int quiet);
int SceneCopyExternal(int width, int height,int rowbytes,unsigned char *dest);

void SceneResetMatrix(void);

void SceneRestartTimers(void);

void ScenePerspective(int flag);

void SceneRotate(float angle,float x,float y,float z);
void SceneTranslate(float x,float y, float z);
void SceneClip(int plane,float movement,char *sele,int state);

void SceneScale(float scale);
void SceneResetNormal(int lines);

void SceneObjectAdd(CObject *obj);
void SceneObjectDel(CObject *obj);
void SceneOriginSet(float *origin,int preserve);
void SceneOriginGet(float *origin);
void SceneWindowSphere(float *location,float radius);
void SceneRelocate(float *location);
Block *SceneGetBlock(void);
void SceneApplyMatrix(float *m);
void SceneSetStereo(int flag);
int SceneGetStereo(void);
void ScenePurgeCopy(void);
void SceneDontCopyNext(void);
void ScenePrepareExit(void);
void SceneGetViewNormal(float *v);
void SceneClipSet(float front,float back);
void SceneGetView(SceneViewType view);
void SceneSetView(SceneViewType view,int quiet);

void SceneToViewElem(CViewElem *elem);
void SceneFromViewElem(CViewElem *elem);
void SceneGetPos(float *pos);
void SceneGetWidthHeight(int *width,int *height);
int SceneMultipick(Multipick *smp);

void SceneSetCardInfo(char *vendor,char *renderer,char *version);
void SceneGetCardInfo(char **vendor,char **renderer,char **version);
int SceneLoadPNG(char *fname,int movie_flag,int quiet);

void SceneSetDefaultView(void);
void SceneApplyRotMatrix(float *src,float *dst);
void SceneRovingDirty(void);
int SceneRovingCheckDirty(void);
void SceneRovingUpdate(void);
void SceneRovingChanged(void);
void SceneRovingPostpone(void);
void SceneCleanupStereo(void);
int SceneReinitialize(void);
void SceneUpdateStereoMode(void);
void SceneSuppressMovieFrame(void);
int SceneClick(Block *block,int button,int x,int y,int mod);
int SceneRelease(Block *block,int button,int x,int y,int mod);
int SceneDrag(Block *block,int x,int y,int mod);
char *SceneGetSeleModeKeyword(void);
void SceneUpdateStereo(void);
#endif



