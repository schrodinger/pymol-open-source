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
#include"Object.h"
#include"Ortho.h"

typedef float SceneViewType[25]; 
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
void SceneRender(Pickable *pick,int x,int y);
void SceneSetFrame(int mode,int frame);
int SceneGetFrame(void);
int SceneGetState(void);
void SceneDirty(void); /* must update transformation */
void SceneChanged(void); /* must actually update 3D objects */
void SceneCountFrames(void) ;
void SceneSetMatrix(float *);
float *SceneGetMatrix(void);

void SceneReshape(Block *block,int width,int height);

void SceneTest(void);
void SceneIdle(void);
void SceneFree(void);
void SceneRay(void);
void SceneMakeMovieImage(void);
void SceneCopy(int buffer);
void SceneRay(void);

void ScenePNG(char *png);
void SceneResetMatrix(void);

void SceneRestartTimers(void);

void ScenePerspective(int flag);

void SceneRotate(float angle,float x,float y,float z);
void SceneTranslate(float x,float y, float z);
void SceneClip(int plane,float movement);

void SceneScale(float scale);
void SceneResetNormal(int lines);

void SceneObjectAdd(Object *obj);
void SceneObjectDel(Object *obj);
void SceneOriginSet(float *origin,int preserve);
void SceneOriginGet(float *origin);
void SceneWindowSphere(float *location,float radius);
Block *SceneGetBlock(void);
void SceneApplyMatrix(float *m);
void SceneSetStereo(int flag);
void ScenePurgeCopy(void);
void SceneDontCopyNext(void);
void ScenePrepareExit(void);
void SceneGetViewNormal(float *v);
void SceneClipSet(float front,float back);
void SceneGetView(SceneViewType view);
void SceneSetView(SceneViewType view);

#endif



