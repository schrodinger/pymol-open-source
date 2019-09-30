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

// this header
#include "OpenVRMode.h"

// system headers
#include "os_std.h"
#include "os_gl.h"
#include "os_python.h"
#include <string>
#include <vector>
#include "openvr.h"

// pymol headers
#include "PyMOLOptions.h"
#include "Setting.h"
#include "Feedback.h"
#include "Matrix.h"
#include "Ortho.h"

// local headers
#include "OpenVRUtils.h"
#include "OpenVRStub.h"
#include "OpenVRController.h"
#include "OpenVRMenu.h"
#include "OpenVRActionList.h"
#include "OpenVRScenePicker.h"
#include "OpenVRLaserTarget.h"

struct CEye {
  vr::EVREye Eye;

  GLfloat HeadToEyeMatrix[16];
  GLfloat ProjectionMatrix[16];

  GLuint FrameBufferID;
  GLuint DepthBufferID;
  GLuint ColorBufferID;

  GLuint ResolveBufferID;
  GLuint ResolveTextureID;

  vr::Texture_t Texture;

  float Left, Right, Top, Bottom; // projection params
};

enum EHand
{
  HLeft = 0,
  HRight = 1,

  Hand_Count
};

enum EUserActionSet
{
  UserActionSet_Mouse,
  UserActionSet_Scene,
  UserActionSet_Movie,

  UserActionSet_Count
};

enum EUserAction
{
  UserAction_None,

  UserAction_Mouse_LClick,
  UserAction_Mouse_MClick,
  UserAction_Mouse_RClick,

  UserAction_Scene_Prev,
  UserAction_Scene_Next,

  UserAction_Movie_Prev,
  UserAction_Movie_Toggle,
  UserAction_Movie_Next,
};

static EUserAction s_userActionMapping[UserActionSet_Count][3] = {
  {UserAction_Mouse_LClick, UserAction_None/*UserAction_Mouse_MClick*/, UserAction_Mouse_RClick}, // UserActionSet_Mouse
  {UserAction_Scene_Prev, UserAction_None, UserAction_Scene_Next}, // UserActionSet_Scene
  {UserAction_Movie_Prev, UserAction_Movie_Toggle, UserAction_Movie_Next}, // UserActionSet_Movie
};

struct COpenVR {
  // Such structures used to be calloc-ed, this replicates that
  void *operator new(size_t size) {
    void *mem = ::operator new(size);
    memset(mem, 0, size);
    return mem;
  }

  vr::EVRInitError InitError;
  vr::IVRSystem* System;
  vr::IVRCompositor* Compositor;
  vr::IVRInput* Input;
  vr::TrackedDevicePose_t Poses[vr::k_unMaxTrackedDeviceCount]; // todo remove from globals?

  GLfloat HeadPose[16];
  GLfloat WorldToHeadMatrix[16];

  unsigned Width;
  unsigned Height;

  CEye* Eye;
  CEye Left;
  CEye Right;

  OpenVRController Hands[Hand_Count];
  GLuint ControllerHintsTexture;

  bool ForcedFront;
  bool ForcedBack;
 
  OpenVRMenu Menu;
  GLuint MenuSplashTexture;

  OpenVRScenePicker Picker;

  OpenVRInputHandlers* Handlers;
  OpenVRActionList* Actions;
  EUserActionSet UserActionSet[Hand_Count];

  // mouse cursor imitation information
  int startX, startY;
  int deltaX, deltaY;

  float moleculeToWorldMatrix[16];
  float moleculeToCapturingController[16];
  float controllersCenterToWorld[16];
  int capturingHandIdx;
  float controllersDistance;
};

static char const* deviceClassNames[] = {
  "Invalid",
  "Head-Mounted Display",
  "Controller",
  "Generic Tracker",
  "Reference Point",
  "Accessory",
};
static const int deviceClassNamesCount = sizeof(deviceClassNames) / sizeof(*deviceClassNames);

struct CMouseEvent {
  unsigned deviceIndex;
  int button;
  int state;
  CMouseEvent() : deviceIndex(0), button(0), state(0) {}
  CMouseEvent(unsigned i, int b, int s) : deviceIndex(i), button(b), state(s) {}
  CMouseEvent(unsigned i, int b, bool s) : deviceIndex(i), button(b), state(s ? P_GLUT_DOWN : P_GLUT_UP) {}
};

class OpenVROrthoInputHandlers : public OpenVRInputHandlers {
  PyMOLGlobals *G;

public:
  explicit OpenVROrthoInputHandlers(PyMOLGlobals *G) : G(G) {}

  void KeyboardFunc(unsigned char k, int x, int y, int mod) override {
    OrthoKey(G, k, x, y, mod);
  }
  void SpecialFunc(int k, int x, int y, int mod) override {
    OrthoSpecial(G, k, x, y, mod);
  }
  int MouseFunc(int button, int state, int x, int y, int mod) override {
    return OrthoButtonDefer(G, button, state, x, y, mod);
  }
  int MotionFunc(int x, int y, int mod) override {
    return OrthoDrag(G, x, y, mod);
  }
  void ActionFunc(int a) override {
    switch (a) {
    case cAction_scene_next:
      OrthoCommandIn(G, "cmd.scene('','next')");
      break;

    case cAction_scene_prev:
      OrthoCommandIn(G, "cmd.scene('','previous')");
      break;

    case cAction_movie_toggle:
      OrthoCommandIn(G, "mtoggle");
      break;

    case cAction_movie_next:
      OrthoCommandIn(G, "forward");
      break;

    case cAction_movie_prev:
      OrthoCommandIn(G, "backward");
      break;
    }

    OrthoInvalidateDoDraw(G);
  }
};

void UpdateDevicePoses(PyMOLGlobals * G);

bool OpenVRAvailable(PyMOLGlobals *)
{
  return vr::stub::VR_IsHmdPresent();
}

bool OpenVRReady(PyMOLGlobals * G) 
{
  COpenVR *I = G->OpenVR;
  return I && I->InitError == vr::VRInitError_None && I->System != NULL;
}

static bool EyeInit(CEye * I, vr::EVREye eye, int scene_width, int scene_height)
{
  I->Eye = eye;

  // framebuffer
  glGenFramebuffersEXT(1, &I->FrameBufferID);
  glBindFramebufferEXT(GL_FRAMEBUFFER, I->FrameBufferID);

  // - depth
  glGenRenderbuffersEXT(1, &I->DepthBufferID);
  glBindRenderbufferEXT(GL_RENDERBUFFER, I->DepthBufferID);
  glRenderbufferStorageMultisampleEXT(GL_RENDERBUFFER, 4, GL_DEPTH_COMPONENT, scene_width, scene_height);
  glFramebufferRenderbufferEXT(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT, GL_RENDERBUFFER, I->DepthBufferID);

  // - color
  glGenTextures(1, &I->ColorBufferID);
  glBindTexture(GL_TEXTURE_2D_MULTISAMPLE, I->ColorBufferID);
  glTexImage2DMultisample(GL_TEXTURE_2D_MULTISAMPLE, 4, GL_RGBA8, scene_width, scene_height, GL_TRUE);
  glFramebufferTexture2DEXT(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D_MULTISAMPLE, I->ColorBufferID, 0);
  glBindTexture(GL_TEXTURE_2D_MULTISAMPLE, 0);
  
  // resolve buffer
  glGenFramebuffersEXT(1, &I->ResolveBufferID);
  glBindFramebufferEXT(GL_FRAMEBUFFER, I->ResolveBufferID);
  
  // - color
  glGenTextures(1, &I->ResolveTextureID);
  glBindTexture(GL_TEXTURE_2D, I->ResolveTextureID);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAX_LEVEL, 0);
  glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA8, scene_width, scene_height, 0, GL_RGBA, GL_UNSIGNED_BYTE, NULL);
  glFramebufferTexture2DEXT(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, I->ResolveTextureID, 0);
  glBindTexture(GL_TEXTURE_2D, 0);

  // VR texture
  I->Texture.handle = (void*)(size_t)I->ResolveTextureID;
  I->Texture.eType = vr::TextureType_OpenGL;
  I->Texture.eColorSpace = vr::ColorSpace_Gamma;

  // check FBO status
  GLenum status = glCheckFramebufferStatusEXT(GL_FRAMEBUFFER);
  glBindFramebufferEXT(GL_FRAMEBUFFER, 0);
  return (status == GL_FRAMEBUFFER_COMPLETE);
}

static void EyeFree(CEye * I)
{
  glDeleteTextures(1, &I->ResolveTextureID);
  glDeleteFramebuffers(1, &I->ResolveBufferID);
  glDeleteTextures(1, &I->ColorBufferID);
  glDeleteRenderbuffers(1, &I->DepthBufferID);
  glDeleteFramebuffers(1, &I->FrameBufferID);
}

int OpenVRInit(PyMOLGlobals * G)
{
  if(G->OpenVR)
    return 1; // already initialized
  
  vr::stub::VR_StubEnable(G->Option->openvr_stub);
  if (!OpenVRAvailable(G))
    return 0; // don't bother initializing the whole system

  COpenVR *I = G->OpenVR = new COpenVR();
  if(!I) 
    return 0;

  I->InitError = vr::VRInitError_None;
  I->System = vr::stub::VR_Init(&I->InitError, vr::VRApplication_Scene);
  if (I->InitError != vr::VRInitError_None) {
    I->System = NULL;
    return 0;
  }

  I->Compositor = vr::stub::VRCompositor();
  I->ForcedFront = true;

  OpenVRSetInputHandlers(G, new OpenVROrthoInputHandlers(G));

  I->Input = vr::stub::VRInput(); 
  if (I->Input) {
    // init manifest
    auto manifestPath = std::string(getenv("PYMOL_DATA"))
                            .append(PATH_SEP "openvr" PATH_SEP "actions.json");
    I->Input->SetActionManifestPath(manifestPath.c_str());

    I->Actions = new OpenVRActionList(I->Input);
  }

  I->capturingHandIdx = -1;
  
  return 1;
}

void OpenVRFree(PyMOLGlobals * G)
{
  ShutdownRenderModels();

  if(!G->OpenVR)
    return;

  COpenVR *I = G->OpenVR;
  if(I->System) {
    vr::stub::VR_Shutdown();

    I->Picker.Free();
    I->Menu.Free();
    delete I->Handlers;

    I->Hands[HLeft].Free();
    I->Hands[HRight].Free();

    EyeFree(&I->Right);
    EyeFree(&I->Left);

    I->System = NULL;
  }

  delete G->OpenVR;
  G->OpenVR = NULL;
}

static void OpenVRInitPostponed(PyMOLGlobals * G)
{
  COpenVR *I = G->OpenVR;
  if(!OpenVRReady(G))
    return;

  if (!I->Width || !I->Height) {
    I->System->GetRecommendedRenderTargetSize(&I->Width, &I->Height);
    EyeInit(&I->Left, vr::Eye_Left, I->Width, I->Height);
    EyeInit(&I->Right, vr::Eye_Right, I->Width, I->Height);

    I->ControllerHintsTexture = OpenVRUtils::LoadTexture("hints_vive_controller.png");
    I->MenuSplashTexture = OpenVRUtils::LoadTexture("menu_splash.png");

    I->Menu.Init(I->Handlers);
    I->Picker.Init(I->Handlers);

    OpenVRMenuSettingsChanged(G);
  }

  float width = SettingGetGlobal_f(G, cSetting_openvr_laser_width);
  for (int i = HLeft; i <= HRight; ++i) {
    OpenVRController &hand = I->Hands[i];
    if (!hand.IsInitialized()) {
      hand.Init();
      hand.SetHintsTexture(I->ControllerHintsTexture, UserActionSet_Count);
      hand.SetLaserWidth(width);
    }
  }
}

void OpenVRSetInputHandlers(PyMOLGlobals * G, OpenVRInputHandlers* handlers)
{
  COpenVR *I = G->OpenVR;
  if (I) {
    if (I->Handlers)
      delete I->Handlers;
    I->Handlers = handlers;
  }
}

static std::string GetStringTrackedDeviceProperty(vr::IVRSystem *System, vr::TrackedDeviceIndex_t index, vr::TrackedDeviceProperty prop)
{
  uint32_t length = System->GetStringTrackedDeviceProperty(index, prop, NULL, 0);
  if(length != 0) {
    std::string buffer(length, 0);
    if (System->GetStringTrackedDeviceProperty(index, prop, &buffer[0], length) != 0) {
      return buffer;
    }
  }

  return std::string("<ERROR>");
}

void OpenVRFeedback(PyMOLGlobals * G)
{
  COpenVR *I = G->OpenVR;
  if(vr::stub::VR_IsStubEnabled()) {
    G->Feedback->add(" OpenVR stub is enabled.\n");
  }
  if(!OpenVRAvailable(G)) {
    G->Feedback->add(" OpenVR system is not available.\n");
  } else if(!OpenVRReady(G)) {
    PRINTF
      " OpenVR system is not ready: %s.\n",
      I ? vr::stub::VR_GetVRInitErrorAsEnglishDescription(I->InitError) : "Failed to initialize properly"
    ENDF(G);
  } else {
    G->Feedback->add(" Detected OpenVR system. Devices being currently tracked:\n");

    bool found = false;
    for(uint32_t i = 0; i < vr::k_unMaxTrackedDeviceCount; ++i) {
      vr::ETrackedDeviceClass deviceClass = I->System->GetTrackedDeviceClass(i);
      if(deviceClass != vr::TrackedDeviceClass_Invalid) {
        found = true;

        char const* className = (0 <= deviceClass && deviceClass < deviceClassNamesCount) ? deviceClassNames[deviceClass] : "<ERROR>";
        std::string model = GetStringTrackedDeviceProperty(I->System, i, vr::Prop_ModelNumber_String);
        std::string serial = GetStringTrackedDeviceProperty(I->System, i, vr::Prop_SerialNumber_String);

        PRINTF "  %02u: %s (%s %s)\n", i, className, model.c_str(), serial.c_str() ENDF(G);
      }
    }
    if(!found) {
      G->Feedback->add("  No valid devices found.\n");
    }
  }
  G->Feedback->add("\n");
}

void OpenVRFrameStart(PyMOLGlobals * G)
{
  COpenVR *I = G->OpenVR;
  if(!OpenVRReady(G))
    return;

  // create OpenGL assets on the first use
  OpenVRInitPostponed(G);

  // get matrices from tracked devices
  if (I->Compositor->WaitGetPoses(I->Poses, vr::k_unMaxTrackedDeviceCount, NULL, 0) != vr::VRCompositorError_None) {
    G->Feedback->add("  Cannot update device poses\n");
  } 
  UpdateDevicePoses(G);
}

void OpenVREyeStart(PyMOLGlobals * G, int eye)
{
  COpenVR *I = G->OpenVR;
  if(!OpenVRReady(G))
    return;

  GL_DEBUG_FUN();

  CEye *E = I->Eye = eye ? &I->Right : &I->Left;

  glBindFramebufferEXT(GL_FRAMEBUFFER, E->FrameBufferID);
  glViewport(0, 0, I->Width, I->Height);
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
}

void OpenVREyeFinish(PyMOLGlobals * G)
{
  COpenVR *I = G->OpenVR;
  if(!OpenVRReady(G))
    return;

  GL_DEBUG_FUN();

  CEye *E = I->Eye;
  if(!E)
    return;

  if(G->Option->multisample)
    glDisable(0x809D);       /* GL_MULTISAMPLE_ARB */

  glBindFramebufferEXT(GL_READ_FRAMEBUFFER, E->FrameBufferID);
  glBindFramebufferEXT(GL_DRAW_FRAMEBUFFER, E->ResolveBufferID);
  glBlitFramebufferEXT(0, 0, I->Width, I->Height, 0, 0, I->Width, I->Height, GL_COLOR_BUFFER_BIT, GL_LINEAR);
  glBindFramebufferEXT(GL_READ_FRAMEBUFFER, 0);
  glBindFramebufferEXT(GL_DRAW_FRAMEBUFFER, 0);
  glBindFramebufferEXT(GL_FRAMEBUFFER, 0);

  if(G->Option->multisample)
    glEnable(0x809D);       /* GL_MULTISAMPLE_ARB */
  
  I->Eye = NULL;
}

void OpenVRSceneFinish(PyMOLGlobals * G, unsigned sceneX, unsigned sceneY, unsigned sceneWidth, unsigned sceneHeight)
{
  COpenVR *I = G->OpenVR;
  if(!OpenVRReady(G))
    return;

  GL_DEBUG_FUN();

  // find a proper rectangle with the scene aspect ratio
  unsigned width = I->Height * sceneWidth / sceneHeight;
  unsigned height = I->Width * sceneHeight / sceneWidth;
  unsigned dx = 0, dy = 0;
  if (width < I->Width) {
    dx = (I->Width - width) / 2;
    height = I->Height;
  } else {
    dy = (I->Height - height) / 2;
    width = I->Width;
  }

  // display a copy of the VR framebuffer in the main PyMOL window
  glDrawBuffer(GL_BACK);
  glBindFramebufferEXT(GL_READ_FRAMEBUFFER, I->Left.ResolveBufferID);
  glBlitFramebufferEXT(dx, dy, dx + width, dy + height, sceneX, sceneY, sceneX + sceneWidth, sceneY + sceneHeight, GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT, GL_NEAREST);
  glBindFramebufferEXT(GL_READ_FRAMEBUFFER, 0);
}

void OpenVRFrameFinish(PyMOLGlobals * G)
{
  COpenVR *I = G->OpenVR;
  if(!OpenVRReady(G))
    return;

  GL_DEBUG_FUN();

  // send rendered pictures into the headset
  I->Compositor->Submit(vr::Eye_Left, &I->Left.Texture);
  I->Compositor->Submit(vr::Eye_Right, &I->Right.Texture);
}

void OpenVRGetWidthHeight(PyMOLGlobals * G, int* width, int* height)
{
  COpenVR *I = G->OpenVR;
  if (I) {
    *width = I->Width;
    *height = I->Height;
  }
}

void OpenVRMenuBufferStart(PyMOLGlobals * G, unsigned width, unsigned height, bool clear /* = true */)
{
  COpenVR *I = G->OpenVR;
  if(!OpenVRReady(G))
    return;

  I->Menu.Start(width, height, clear);
}

void OpenVRMenuBufferFinish(PyMOLGlobals * G)
{
  COpenVR *I = G->OpenVR;
  if(!OpenVRReady(G))
    return;

  I->Menu.Finish();
}

void OpenVRMenuToggle(PyMOLGlobals * G, unsigned deviceIndex /* = ~0U */)
{
  COpenVR *I = G->OpenVR;
  if(!OpenVRReady(G))
    return;

  if (!I->Menu.IsVisible()) {
    I->Menu.Show(I->HeadPose, deviceIndex);
    I->ForcedBack = true;
  } else {
    unsigned ownerIndex = I->Menu.GetOwnerID();
    if (deviceIndex == ownerIndex || deviceIndex == ~0U || ownerIndex == ~0U) {
      I->Menu.Hide();
      I->ForcedBack = false;
    }
  }
}

void OpenVRMenuCrop(PyMOLGlobals * G, unsigned x, unsigned y, unsigned width, unsigned height)
{
  COpenVR *I = G->OpenVR;
  if(!OpenVRReady(G))
    return;

  I->Menu.Crop(x, y, width, height);
}

void OpenVRMenuSettingsChanged(PyMOLGlobals * G)
{
  COpenVR *I = G->OpenVR;
  if(!OpenVRReady(G))
    return;

  float distance = SettingGetGlobal_f(G, cSetting_openvr_gui_distance);
  float fov = SettingGetGlobal_f(G, cSetting_openvr_gui_fov);
  I->Menu.SetSize(distance, tanf(fov * PI / 180.0f));

  float sceneColor = SettingGetGlobal_f(G, cSetting_openvr_gui_scene_color);
  float sceneAlpha = SettingGetGlobal_f(G, cSetting_openvr_gui_scene_alpha);
  I->Menu.SetSceneColor(sceneColor, sceneAlpha);

  float backColor = SettingGetGlobal_f(G, cSetting_openvr_gui_back_color);
  float backAlpha = SettingGetGlobal_f(G, cSetting_openvr_gui_back_alpha);
  I->Menu.SetBackdropColor(backColor, backAlpha);
}

float* OpenVRGetHeadToEye(PyMOLGlobals * G)
{
  COpenVR *I = G->OpenVR;
  if(!OpenVRReady(G) || !I->Eye)
    return NULL;

  CEye *E = I->Eye;
  vr::HmdMatrix34_t EyeToHeadTransform = I->System->GetEyeToHeadTransform(E->Eye);
  OpenVRUtils::MatrixFastInverseVRGL((const float *)EyeToHeadTransform.m, E->HeadToEyeMatrix);
  
  return E->HeadToEyeMatrix;
}

float* OpenVRGetWorldToHead(PyMOLGlobals * G) {
  COpenVR *I = G->OpenVR;
  if(!OpenVRReady(G))
    return NULL;

  return I->WorldToHeadMatrix;
}

float* OpenVRGetControllerPose(PyMOLGlobals * G, EHand handIdx) {
  COpenVR *I = G->OpenVR;
  if(!OpenVRReady(G))
    return NULL;

  return I->Hands[handIdx].GetPose();
}

void OpenVRGetProjection(float left, float right, float top, float bottom, float near_plane, float far_plane, float *matrix)
{
  if (!matrix)
    return;

  // fast affine inverse matrix, row major to column major, whew...
  {
    float (*dst)[4] = (float(*)[4])matrix;
    float dx = (right - left);
    float dy = (bottom - top);
    float dz = far_plane - near_plane;

    // transpose rotation
    dst[0][0] = 2.0f / dx;
    dst[0][1] = 0.0f;
    dst[0][2] = 0.0f;
    dst[0][3] = 0.0f;
    
    dst[1][0] = 0.0f;
    dst[1][1] = 2.0f / dy;
    dst[1][2] = 0.0f;
    dst[1][3] = 0.0f;
    
    dst[2][0] = (right + left) / dx;
    dst[2][1] = (top + bottom) / dy;
    dst[2][2] = -(far_plane + near_plane) / dz;
    dst[2][3] = -1.0f;
    
    dst[3][0] = 0.0f;
    dst[3][1] = 0.0f;
    dst[3][2] = -2.0f * far_plane * near_plane / dz;
    dst[3][3] = 0.0f;
  }

  return;
}

void CheckNearFarPlaneSettings(PyMOLGlobals * G, float &near_plane, float &far_plane) {
  COpenVR *I = G->OpenVR;

  if (I->ForcedFront || SettingGetGlobal_b(G, cSetting_openvr_disable_clipping)) {
    near_plane = SettingGetGlobal_f(G, cSetting_openvr_near_plane);
  }
  if (I->ForcedBack || SettingGetGlobal_b(G, cSetting_openvr_disable_clipping)) {
    far_plane = SettingGetGlobal_f(G, cSetting_openvr_far_plane);
  }
}

float* OpenVRGetEyeProjection(PyMOLGlobals * G, float near_plane, float far_plane)
{
  COpenVR *I = G->OpenVR;
  if(!OpenVRReady(G) || !I->Eye)
    return NULL;

  CheckNearFarPlaneSettings(G, near_plane, far_plane);
 
  CEye *E = I->Eye;
  I->System->GetProjectionRaw(E->Eye, &(E->Left), &(E->Right), &(E->Top), &(E->Bottom));
  OpenVRGetProjection(E->Left, E->Right, E->Top, E->Bottom, near_plane, far_plane, E->ProjectionMatrix);
  return E->ProjectionMatrix;
}

void  OpenVRGetPickingProjection(PyMOLGlobals * G, float near_plane, float far_plane, float *matrix)
{
  COpenVR *I = G->OpenVR;
  if(!OpenVRReady(G))
    return;

  CheckNearFarPlaneSettings(G, near_plane, far_plane);

  // take avarage projection params from eyes
  float left, right, top, bottom;
  CEye &LEye = I->Left, &REye = I->Right;
  left = (LEye.Left + REye.Left) * 0.5f;
  right = (LEye.Right + REye.Right) * 0.5f;
  top = (LEye.Top + REye.Top) * 0.5f;
  bottom = (LEye.Bottom + REye.Bottom) * 0.5f;
  OpenVRGetProjection(left, right, top, bottom, near_plane, far_plane, matrix);
  return;
}

float const* OpenVRGetPickingMatrix(PyMOLGlobals * G)
{
  COpenVR *I = G->OpenVR;
  if(!OpenVRReady(G))
    return NULL;

  return I->Picker.GetMatrix();
}

void OpenVRLoadProjectionMatrix(PyMOLGlobals * G, float near_plane, float far_plane)
{
  glLoadMatrixf(OpenVRGetEyeProjection(G, near_plane, far_plane));
}

void OpenVRLoadPickingProjectionMatrix(PyMOLGlobals * G, float near_plane, float far_plane)
{
  float matrix[16];
  OpenVRGetPickingProjection(G, near_plane, far_plane, matrix);
  glLoadMatrixf(matrix);
}

void OpenVRLoadWorld2EyeMatrix(PyMOLGlobals * G)
{
  glLoadMatrixf(OpenVRGetHeadToEye(G));
  glMultMatrixf(OpenVRGetWorldToHead(G));
}

std::string GetTrackedDeviceString(PyMOLGlobals * G, vr::TrackedDeviceIndex_t unDevice, vr::TrackedDeviceProperty prop, vr::TrackedPropertyError *peError = NULL )
{
  COpenVR *I = G->OpenVR;
  if (!I || !I->System) 
    return "";

  uint32_t unRequiredBufferLen = I->System->GetStringTrackedDeviceProperty( unDevice, prop, NULL, 0, peError );
  if( unRequiredBufferLen == 0 )
    return "";

  char *pchBuffer = new char[ unRequiredBufferLen ];
  unRequiredBufferLen = I->System->GetStringTrackedDeviceProperty( unDevice, prop, pchBuffer, unRequiredBufferLen, peError );
  std::string sResult = pchBuffer;
  delete [] pchBuffer;
  return sResult;
}

void ProcessButtonDragAsMouse(PyMOLGlobals * G, OpenVRAction *action, int glutButton, int screenCenterX, int screenCenterY) {
  COpenVR *I = G->OpenVR;
  if (!action || !I) return;

  OpenVRInputHandlers* Handlers = I->Handlers;
  if (!Handlers)
   return;

  // imitate mouse cursor position from controller camera position
  float (*mat)[4] = (float (*)[4])I->Hands[HRight].GetPose();
  int x = (int)(mat[3][0]* 500.0f), y = (int)(mat[3][1] * 500.0f); // magic factors

  bool nowPressed = action->IsPressed();
  if (action->WasPressedOrReleased()) {
    if (nowPressed) {
      I->startX = x;
      I->startY = y;
      Handlers->MouseFunc(glutButton, P_GLUT_DOWN, screenCenterX, screenCenterY, 0);
    } else {
      Handlers->MouseFunc(glutButton, P_GLUT_UP, I->deltaX + screenCenterX, I->deltaY + screenCenterY, 0);
    }
  }
  if (nowPressed) {
    I->deltaX = x - I->startX;
    I->deltaY = y - I->startY;
    Handlers->MotionFunc(I->deltaX + screenCenterX, I->deltaY + screenCenterY, 0);
  }
}

bool OpenVRIsMoleculeCaptured(PyMOLGlobals * G) {
  COpenVR *I = G->OpenVR;
  return I->Hands[HLeft].isGripPressed() || I->Hands[HRight].isGripPressed();
}

void CalculateScalingPivotToWorldMatrix(PyMOLGlobals * G, float *pivotToWorldMatrix) {
  COpenVR *I = G->OpenVR;
  if (!I || !pivotToWorldMatrix)
    return;

  identity44f(pivotToWorldMatrix);
  float *lPose = I->Hands[HLeft].GetPose();
  float *rPose = I->Hands[HRight].GetPose();
  // get center traslation matrix
  average3f(&lPose[12], &rPose[12], &pivotToWorldMatrix[12]); 
}

float const *OpenVRGetMolecule2WorldMatrix(PyMOLGlobals * G, float *scaler) {
  COpenVR *I = G->OpenVR;
  if (!I)
    return NULL;

  float temp[16];
  OpenVRActionList* Actions = I->Actions;
  if (Actions->LGrip->IsPressed() && Actions->RGrip->IsPressed()) {
    // translate after the scaling pivot
    float pivotToWorldMatrix[16];
    CalculateScalingPivotToWorldMatrix(G, pivotToWorldMatrix);
    memcpy(temp, pivotToWorldMatrix, sizeof(temp));
    // scale due to changing distance between controllers
    if (scaler) {
      float newDistance = diff3f(&(I->Hands[HLeft].GetPose()[12]), &(I->Hands[HRight].GetPose()[12]));
      *scaler = newDistance / I->controllersDistance;
      I->controllersDistance = newDistance;
    }
  } else {
    memcpy(temp, I->Hands[I->capturingHandIdx].GetPose(), sizeof(temp));
    *scaler = 1.0f;
  }
  MatrixMultiplyC44f(I->moleculeToCapturingController, temp);
  memcpy(I->moleculeToWorldMatrix, temp, sizeof(I->moleculeToWorldMatrix));
  return I->moleculeToWorldMatrix;
}

void AttachMoleculeToController(PyMOLGlobals * G, int handIdx) {
  COpenVR *I = G->OpenVR;
  I->capturingHandIdx = handIdx;
  // catch up
  memcpy(I->moleculeToCapturingController, I->Hands[handIdx].GetWorldToControllerMatrix(), sizeof(I->moleculeToCapturingController));
  MatrixMultiplyC44f(I->moleculeToWorldMatrix, I->moleculeToCapturingController);
}

void AttachMoleculeToCenter(PyMOLGlobals * G) {
  COpenVR *I = G->OpenVR;

  // get scaling pivot center = controllers mass center
  float worldToPivotMatrix[16];
  CalculateScalingPivotToWorldMatrix(G, worldToPivotMatrix);
  // inverse transform to be exactly WorldToPivot
  worldToPivotMatrix[12] *= -1.0;
  worldToPivotMatrix[13] *= -1.0;
  worldToPivotMatrix[14] *= -1.0;
  // attach molecule to pivot
  memcpy(I->moleculeToCapturingController, worldToPivotMatrix, sizeof(I->moleculeToCapturingController));
  MatrixMultiplyC44f(I->moleculeToWorldMatrix, I->moleculeToCapturingController);
  // save distance between controllers
  I->controllersDistance = diff3f(&(I->Hands[HLeft].GetPose()[12]), &(I->Hands[HRight].GetPose()[12]));
}

void HandleLaser(PyMOLGlobals * G, int centerX, int centerY, CMouseEvent const& mouseEvent);

void OpenVRHandleInput(PyMOLGlobals * G, int SceneX, int SceneY, int SceneWidth, int SceneHeight, float *model2World)
{
  COpenVR *I = G->OpenVR;
  if(!OpenVRReady(G))
    return;

  vr::VREvent_t event;
  while (I->System->PollNextEvent(&event, sizeof(event)))
    /* pass */;

  if (!I->Input)
    return;

  OpenVRActionList* Actions = I->Actions;
  Actions->Update(I->Input);

  OpenVRController& LeftHand = I->Hands[HLeft];
  OpenVRController& RightHand = I->Hands[HRight];
  int centerX = SceneX + SceneWidth / 2;
  int centerY = SceneY + SceneHeight / 2;

  // update VR GUI state
  if (Actions->ToggleMenu->WasPressed()) {
    OpenVRMenuToggle(G);
  }

  // update controllers visibility
  LeftHand.Show(Actions->LeftHand->PoseValid());
  RightHand.Show(Actions->RightHand->PoseValid());

  // process grips
  {
    I->Hands[HLeft].pressGrip(Actions->LGrip->IsPressed());
    I->Hands[HRight].pressGrip(Actions->RGrip->IsPressed());

    if (OpenVRIsMoleculeCaptured(G)) {
      memcpy(I->moleculeToWorldMatrix, model2World, sizeof(I->moleculeToWorldMatrix)); 
    }
    if (Actions->LGrip->WasPressed() && !Actions->RGrip->IsPressed()) {
      AttachMoleculeToController(G, HLeft);
    }
    if (Actions->RGrip->WasPressed() && !Actions->LGrip->IsPressed()) {
      AttachMoleculeToController(G, HRight);
    }
    // TODO make it being less ugly
    if ((Actions->RGrip->WasPressed() && Actions->LGrip->IsPressed()) ||
       (Actions->LGrip->WasPressed() && Actions->RGrip->IsPressed())) {
      AttachMoleculeToCenter(G);
    }
    if (Actions->LGrip->WasReleased() && Actions->RGrip->IsPressed()) {
      AttachMoleculeToController(G, HRight);
    }
    if (Actions->RGrip->WasReleased() && Actions->LGrip->IsPressed()) {
      AttachMoleculeToController(G, HLeft);
    }
    if ((Actions->LGrip->WasReleased() && !Actions->RGrip->IsPressed()) ||
        (Actions->RGrip->WasReleased() && !Actions->LGrip->IsPressed())) {
      I->capturingHandIdx = -1;
    }
  }

  // switch user action sets
  {
    int increment = 0;
    unsigned deviceIndex = 0;

    if (Actions->ActionSetNext->WasPressed()) {
      deviceIndex = Actions->ActionSetNext->DeviceIndex();
      increment = +1;
    } else if (Actions->ActionSetPrev->WasPressed()) {
      deviceIndex = Actions->ActionSetPrev->DeviceIndex();
      increment = -1;
    }

    if (increment) {
      EHand handIndex = EHand(deviceIndex == RightHand.m_deviceIndex);
      I->UserActionSet[handIndex] = EUserActionSet((I->UserActionSet[handIndex] + UserActionSet_Count + increment) % UserActionSet_Count);
      I->Hands[handIndex].SetHintsIndex(I->UserActionSet[handIndex]);
    }
  }

  // process user actions
  CMouseEvent mouseEvent;
  OpenVRAction* userActions[] = {Actions->Action1, Actions->Action2, Actions->Action3};
  for (int i = 0, n = sizeof(userActions) / sizeof(*userActions); i < n; ++i) {
    OpenVRAction* action = userActions[i];
    if (action->WasPressedOrReleased()) {
      EHand handIndex = EHand(action->DeviceIndex() == RightHand.m_deviceIndex);
      EUserActionSet userActionSet = I->UserActionSet[handIndex];
      EUserAction userAction = s_userActionMapping[userActionSet][i];

      switch (userAction) {
      case UserAction_Mouse_LClick:
        mouseEvent = CMouseEvent(action->DeviceIndex(), P_GLUT_LEFT_BUTTON, action->IsPressed());
        break;
      case UserAction_Mouse_MClick:
        mouseEvent = CMouseEvent(action->DeviceIndex(), P_GLUT_MIDDLE_BUTTON, action->IsPressed());
        break;
      case UserAction_Mouse_RClick:
        mouseEvent = CMouseEvent(action->DeviceIndex(), P_GLUT_RIGHT_BUTTON, action->IsPressed());
        break;
      case UserAction_Scene_Prev:
        if (action->IsPressed())
          I->Handlers->ActionFunc(cAction_scene_prev);
        break;
      case UserAction_Scene_Next:
        if (action->IsPressed())
          I->Handlers->ActionFunc(cAction_scene_next);
        break;
      case UserAction_Movie_Prev:
        if (action->IsPressed())
          I->Handlers->ActionFunc(cAction_movie_prev);
        break;
      case UserAction_Movie_Toggle:
        if (action->IsPressed())
          I->Handlers->ActionFunc(cAction_movie_toggle);
        break;
      case UserAction_Movie_Next:
        if (action->IsPressed())
          I->Handlers->ActionFunc(cAction_movie_next);
        break;
      }
    }
  }

  HandleLaser(G, centerX, centerY, mouseEvent);
}

void HandleLaser(PyMOLGlobals * G, int centerX, int centerY, CMouseEvent const& mouseEvent)
{
  COpenVR* I = G->OpenVR;
  OpenVRActionList* Actions = I->Actions;

  // hide all lasers
  I->Menu.HideHotspot();
  for (size_t laserIndex = 0; laserIndex < sizeof(I->Hands) / sizeof(I->Hands[0]); ++laserIndex) {
    I->Hands[laserIndex].LaserShow(false);
  }

  // detect a laser source
  OpenVRLaserSource* laserSource = 0;
  if (Actions->Laser->IsPressed()) {
    for (size_t laserIndex = 0; laserIndex < sizeof(I->Hands) / sizeof(I->Hands[0]); ++laserIndex) {
      if (Actions->Laser->DeviceIndex() == I->Hands[laserIndex].GetLaserDeviceIndex() && I->UserActionSet[laserIndex] == UserActionSet_Mouse) {
        laserSource = &I->Hands[laserIndex];
        break;
      }
    }
  }

  bool menuHit = false;

  if (laserSource) {
    I->Picker.Activate(laserSource->GetLaserDeviceIndex(), centerX, centerY);

    float origin[3], dir[3];
    laserSource->LaserShow(true);
    laserSource->GetLaserRay(origin, dir);

    // shoot the laser
    OpenVRLaserTarget* laserTarget = 0;
    OpenVRLaserTarget* targets[] = {&I->Menu, &I->Picker};
    for (int i = 0, n = sizeof(targets) / sizeof(targets[0]); i < n && !laserTarget; ++i) {
      float distance = 0.0f;
      if (targets[i]->IsLaserAllowed(laserSource->GetLaserDeviceIndex()) && 
          targets[i]->LaserShoot(origin, dir, targets[i]->GetLaserColor(), &distance)) {
        laserTarget = targets[i];
        laserSource->SetLaserLength(distance);
        laserSource->SetLaserColor(laserTarget->GetLaserColor());
        if (laserTarget == &I->Menu)
          menuHit = true;
      }
    }

    // laser missed
    float missedColor[4] = {1.0f, 1.0f, 0.0f, 0.5f};
    if (!laserTarget) {
      laserTarget = &I->Picker;
      if (!SettingGetGlobal_b(G, cSetting_openvr_cut_laser)) {
        laserSource->SetLaserLength(0.0f);
      }
      laserSource->SetLaserColor(missedColor);
    }

    if (mouseEvent.deviceIndex == laserSource->GetLaserDeviceIndex()) {
      laserTarget->LaserClick(mouseEvent.button, mouseEvent.state);
    }

  } else {
    I->Picker.Deactivate();
  }

  float alpha = SettingGetGlobal_f(G, cSetting_openvr_gui_alpha);
  int useAlpha = SettingGetGlobal_i(G, cSetting_openvr_gui_use_alpha);
  I->Menu.SetAlpha(useAlpha == 0 || useAlpha == 2 && menuHit ? 1.0f : alpha);

  int useBackdrop = SettingGetGlobal_i(G, cSetting_openvr_gui_use_backdrop);
  I->Menu.SetBackdrop(useBackdrop == 1 || useBackdrop == 2 && menuHit ? true : false);

  int overlay = SettingGetGlobal_i(G, cSetting_openvr_gui_overlay);
  I->Menu.SetOverlay(overlay == 1 || overlay == 2 && menuHit ? true : false);
}

void UpdateDevicePoses(PyMOLGlobals * G) {
  COpenVR *I = G->OpenVR;

  for (uint32_t nDevice = 0; nDevice < vr::k_unMaxTrackedDeviceCount; nDevice++) {
    vr::TrackedDevicePose_t &pose = I->Poses[nDevice];
    if (pose.bPoseIsValid) {
      vr::ETrackedDeviceClass device = I->System->GetTrackedDeviceClass(nDevice);
      switch (device) {
        case vr::TrackedDeviceClass_HMD:
          OpenVRUtils::MatrixCopyVRGL((const float *)pose.mDeviceToAbsoluteTracking.m, I->HeadPose);
          OpenVRUtils::MatrixFastInverseVRGL((const float *)pose.mDeviceToAbsoluteTracking.m, I->WorldToHeadMatrix);
          break;
        case vr::TrackedDeviceClass_Controller:
          {
            vr::ETrackedControllerRole role = I->System->GetControllerRoleForTrackedDeviceIndex(nDevice);

            OpenVRController* hand = 0;
            if (role == vr::TrackedControllerRole_LeftHand) {
              hand = &I->Hands[HLeft];
            } else if (role == vr::TrackedControllerRole_RightHand) {
              hand = &I->Hands[HRight];
            }

            if (hand) {
              OpenVRUtils::MatrixCopyVRGL((const float *)pose.mDeviceToAbsoluteTracking.m, (float *)hand->GetPose());
              OpenVRUtils::MatrixFastInverseVRGL((const float *)pose.mDeviceToAbsoluteTracking.m, (float *)hand->GetWorldToControllerMatrix());
              hand->m_deviceIndex = nDevice;
              std::string sRenderModelName = GetTrackedDeviceString(G, nDevice, vr::Prop_RenderModelName_String);
              if (sRenderModelName != hand->m_sRenderModelName) {
                hand->m_pRenderModel = FindOrLoadRenderModel(G, sRenderModelName.c_str());
                hand->m_sRenderModelName = sRenderModelName;
              }
            }
          }
          break;
        default:
          break;
      }
    }
  }
}

void OpenVRDraw(PyMOLGlobals * G)
{
  COpenVR *I = G->OpenVR;
  if(!OpenVRReady(G))
    return;

  GL_DEBUG_FUN();

  glPushMatrix();
  OpenVRLoadWorld2EyeMatrix(G);

  // render menu if present
  I->Menu.Draw(I->MenuSplashTexture);

  // render controllers
  for (int i = HLeft; i <= HRight; ++i) {
    I->Hands[i].Draw();
  }

  glPopMatrix();
}

void OpenVRClippingChanged(PyMOLGlobals * G) {
  static bool s_oldDepthCue = true;
  bool clipping = SettingGetGlobal_b(G, cSetting_openvr_disable_clipping);
  if (clipping) {
    s_oldDepthCue = SettingGetGlobal_b(G, cSetting_depth_cue);
    SettingSetGlobal_b(G, cSetting_depth_cue, false);
  } else {
    SettingSetGlobal_b(G, cSetting_depth_cue, s_oldDepthCue);
  }
}

void OpenVRLaserWidthChanged(PyMOLGlobals * G) {
  COpenVR *I = G->OpenVR;
  float width = SettingGetGlobal_f(G, cSetting_openvr_laser_width);
  for (int i = HLeft; i <= HRight; ++i) {
    OpenVRController &hand = I->Hands[i];
    hand.SetLaserWidth(width);
  }
}

void OpenVRUpdateScenePickerLength(PyMOLGlobals * G, float *PickWorldPoint)
{
  COpenVR *I = G->OpenVR;
  if(!OpenVRReady(G))
    return;
  
  // get active ray
  for (int i = HLeft; i <= HRight; ++i) {
    OpenVRController &hand = I->Hands[i];
    if (hand.IsLaserVisible()) {
      // get ray start point in world
      float laserStartPointWorld[3];
      if (hand.GetLaserRay(laserStartPointWorld, NULL)) {
        // calc distance to pointed atom in world CS
        float dist = diff3f(laserStartPointWorld, PickWorldPoint);
        // set new ray length
        hand.SetLaserLength(dist);
      }
    }
  }
}

bool OpenVRIsScenePickerActive(PyMOLGlobals * G) {
  COpenVR *I = G->OpenVR;
  return (I && I->Picker.IsActive());
}
