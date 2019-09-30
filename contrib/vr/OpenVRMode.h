/*
OpenVR for PyMOL Copyright Notice
=====================================

The OpenVR for PyMOL source code is copyrighted, but you can freely use and
copy it as long as you don't change or remove any of the Copyright notices.
OpenVR for PyMOL is made available under the following open-source license
terms:

------------------------------------------------------------------------------
Copyright (c) 2018 EPAM Systems, Inc.

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
------------------------------------------------------------------------------

*/

#ifndef _H_OpenVRStereo
#define _H_OpenVRStereo

// pymol headers
#include "PyMOLGlobals.h"

bool OpenVRAvailable(PyMOLGlobals * G);
bool OpenVRReady(PyMOLGlobals * G);

int OpenVRInit(PyMOLGlobals * G);
void OpenVRFree(PyMOLGlobals * G);

// PyMOL-style enum
enum {
  cAction_scene_next,
  cAction_scene_prev,

  cAction_movie_toggle,
  cAction_movie_next,
  cAction_movie_prev,
};

class OpenVRInputHandlers {
public:
  virtual void KeyboardFunc(unsigned char k, int x, int y, int mod) {}
  virtual void SpecialFunc(int k, int x, int y, int mod) {}
  virtual int MouseFunc(int button, int state, int x, int y, int mod) { return 0; }
  virtual int MotionFunc(int x, int y, int mod) { return 0; }
  virtual void ActionFunc(int a) {}
};

void OpenVRSetInputHandlers(PyMOLGlobals * G, OpenVRInputHandlers* handlers);

void OpenVRFeedback(PyMOLGlobals * G);

void OpenVRFrameStart(PyMOLGlobals * G);
void OpenVREyeStart(PyMOLGlobals * G, int eye);
void OpenVREyeFinish(PyMOLGlobals * G);
void OpenVRSceneFinish(PyMOLGlobals * G, unsigned sceneX, unsigned sceneY, unsigned sceneWidth, unsigned sceneHeight);
void OpenVRFrameFinish(PyMOLGlobals * G);

void OpenVRGetWidthHeight(PyMOLGlobals * G, int* width, int* height);

void OpenVRMenuBufferStart(PyMOLGlobals * G, unsigned width, unsigned height, bool clear = true);
void OpenVRMenuBufferFinish(PyMOLGlobals * G);
void OpenVRMenuToggle(PyMOLGlobals * G, unsigned deviceIndex = ~0U);
void OpenVRMenuCrop(PyMOLGlobals * G, unsigned x, unsigned y, unsigned width, unsigned height);
void OpenVRMenuSettingsChanged(PyMOLGlobals * G);

float* OpenVRGetWorldToHead(PyMOLGlobals * G);
float* OpenVRGetHeadToEye(PyMOLGlobals * G);
float* OpenVRGetControllerPose(PyMOLGlobals * G);
float* OpenVRGetEyeProjection(PyMOLGlobals * G, float near_plane, float far_plane);
float const* OpenVRGetPickingMatrix(PyMOLGlobals * G);
void OpenVRGetPickingProjection(PyMOLGlobals * G, float near_plane, float far_plane, float *matrix);
float const *OpenVRGetMolecule2WorldMatrix(PyMOLGlobals * G, float *scaler);

void OpenVRLoadProjectionMatrix(PyMOLGlobals * G, float near_plane, float far_plane);
void OpenVRLoadPickingProjectionMatrix(PyMOLGlobals * G, float near_plane, float far_plane);
void OpenVRLoadWorld2EyeMatrix(PyMOLGlobals * G);

bool OpenVRIsMoleculeCaptured(PyMOLGlobals * G);

void OpenVRHandleInput(PyMOLGlobals * G, int SceneX, int SceneY, int SceneWidth, int SceneHeight, float *model2World);

void OpenVRDraw(PyMOLGlobals * G);

void OpenVRClippingChanged(PyMOLGlobals * G);
void OpenVRLaserWidthChanged(PyMOLGlobals * G);

void OpenVRUpdateScenePickerLength(PyMOLGlobals * G, float *PickWorldPoint);
bool OpenVRIsScenePickerActive(PyMOLGlobals * G);

#endif /* _H_OpenVRStereo */
