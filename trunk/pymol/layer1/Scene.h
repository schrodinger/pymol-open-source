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

void SceneInit(void);
void SceneDone(void);
void SceneUpdate(void);
int SceneRenderCached(void);
void SceneRender(Pickable *pick,int x,int y);
void SceneSetFrame(int mode,int frame);
int SceneGetFrame(void);
void SceneDirty(void); /* must update transformation */
void SceneChanged(void); /* must update actually 3D objects */
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
void SceneWindowSphere(float *location,float radius);
Block *SceneGetBlock(void);
void SceneApplyMatrix(float *m);

#endif



