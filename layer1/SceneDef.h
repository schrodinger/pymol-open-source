
/* 
A* -------------------------------------------------------------------
B* This file contains source code for the PyMOL computer program
C* Copyright (c) Schrodinger, LLC. 
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

#ifndef _H_SceneDef
#define _H_SceneDef

#include"gl_def.h"
#include"Ortho.h"
#include"Util.h"
#include"View.h"
#include"Image.h"
#include "Picking.h"
#include"ScrollBar.h"
#include"SceneElem.h"
#include"SceneView.h"
#include"Rect.h"
#include "Camera.h"
#include<list>
#include<vector>

#include <glm/mat4x4.hpp>

#define TRN_BKG 0x30
#define MAX_ANI_ELEM 300

namespace pymol
{
  struct CObject;
}

typedef struct {
  float unit_left, unit_right, unit_top, unit_bottom, unit_front, unit_back;
} SceneUnitContext;


struct Offset2D
{
  std::int32_t x;
  std::int32_t y;
};

struct Extent2D
{
  std::uint32_t width;
  std::uint32_t height;
};

struct Rect2D
{
  Offset2D offset;
  Extent2D extent;
};


inline bool operator==(const Rect2D& rectA, const Rect2D& rectB) {
  return rectA.offset.x == rectB.offset.x && rectA.offset.y == rectB.offset.y &&
         rectA.extent.width == rectB.extent.width &&
         rectA.extent.height == rectB.extent.height;
}

inline bool operator!=(const Rect2D& rectA, const Rect2D& rectB) {
  return !(rectA == rectB);
}

struct GridInfo {
  int n_col;
  int n_row;
  int first_slot;
  int last_slot;
  float asp_adjust;
  int active;
  int size;
  int slot;
  GridMode mode;
  Rect2D cur_view;
  Extent2D cur_viewport_size;
  SceneUnitContext context;     /* for whole-window display */
};

using PrepareViewportForStereoFuncT = void (*)(
    PyMOLGlobals*, CScene*, int, bool, int, const Offset2D&,
    const Extent2D&);

class CScene : public Block {
 public:
  std::list<pymol::CObject*> Obj, GadgetObjs, NonGadgetObjs;
  pymol::Camera m_view{};
  float InvMatrix[16]{};          /* WARNING: column major, as per OpenGL spec */
  float PmvMatrix[16]{};
  float Scale{1.0F};
  int Width{640}, Height{480};
  int Button{};
  int LastX{}, LastY{};
  int StartX{}, StartY{};
  int LastWinX{}, LastWinY{};
  double LastClickTime{};
  int LastButton{}, LastMod{};
  int PossibleSingleClick{};
  double LastReleaseTime{};
  double SingleClickDelay{};
  float ViewNormal[3]{}, LinesNormal[3]{};
  float TextColor[3]{0.2F, 1.0F, 0.2F};
  double SweepTime{};
  bool DirtyFlag{true};
  bool ChangedFlag{};
  int CopyType{};
  bool CopyNextFlag{true}, CopyForced{};
  int NFrame { 0 };
  int HasMovie { 0 };
  std::shared_ptr<pymol::Image> Image { nullptr };
  bool MovieFrameFlag{};
  double LastRender{}, RenderTime{}, LastFrameTime{}, LastFrameAdjust{};
  double LastSweep{}, LastSweepTime{};
  float LastSweepX{}, LastSweepY{};
  int RockFrame{};
  Picking LastPicked{};
  int StereoMode{};
  OrthoLineType vendor{}, renderer{}, version{};
  bool SculptingFlag{};
  int SculptingSave{};
  bool RovingDirtyFlag{};
  bool RovingCleanupFlag{};
  double RovingLastUpdate{};
  int Threshold{}, ThresholdX{}, ThresholdY{};
  float LastPickVertex[3]{}, LastClickVertex[3]{};
  bool LastPickVertexFlag{};
  bool LoopFlag{};
  int LoopMod{};
  BlockRect LoopRect{};
  CViewElem ani_elem[MAX_ANI_ELEM + 1]{};
  int cur_ani_elem{}, n_ani_elem{};
  int LastStateBuilt{-1};
  bool AnimationStartFlag{};
  double AnimationStartTime{};
  double AnimationLagTime{};
  int AnimationStartFrame{};
  double ApproxRenderTime{};
  float VertexScale{0.01F};
  float FogStart{};
  float FogEnd{};

  /* Scene Names */
  int ButtonsShown{}, ButtonDrag{}, ButtonMargin{}, ButtonsValid{};
  int Over{-1}, Pressed{-1}, PressMode{}, HowFarDown{}, NSkip{};
  int ScrollBarActive{};
  OrthoLineType ReorderLog{};
  ScrollBar m_ScrollBar;
  std::vector<SceneElem> SceneVec;
  CGO *AlphaCGO { nullptr };

  std::vector<int> m_slots;

  int StencilValid{}, StencilParity{};
  bool ReinterpolateFlag{};
  pymol::CObject* ReinterpolateObj{nullptr};
  pymol::CObject* MotionGrabbedObj{nullptr};

  short prev_no_z_rotation1{}, prev_no_z_rotation2{};
  int orig_x_rotation{}, orig_y_rotation{};

  std::vector<glm::mat4> m_ModelViewMatrixStack;
  glm::mat4 modelViewMatrix{};

  glm::mat4 projectionMatrix{};
  int background_color_already_set{};
  int do_not_clear{};
  GridInfo grid{};
  int last_grid_size{};
  int n_texture_refreshes { 0 };
  CGO *offscreenCGO { nullptr };
  CGO *offscreenOIT_CGO { nullptr };
  CGO *offscreenOIT_CGO_copy { nullptr };
  PrepareViewportForStereoFuncT vp_prepareViewPortForStereo = nullptr;
  Offset2D vp_pos{};
  Extent2D vp_oversize{};
  int vp_times{}, vp_stereo_mode{};
  float vp_width_scale{};
  PickColorManager pickmgr;

  CScene(PyMOLGlobals * G) : Block(G), m_ScrollBar(G, false) {}

  virtual int click(int button, int x, int y, int mod) override;
  virtual int release(int button, int x, int y, int mod) override;
  virtual int drag(int x, int y, int mod) override;
  virtual void draw(CGO* orthoCGO) override;
  virtual void reshape(int width, int height) override;

  SceneView getSceneView() const { return m_view.getView(); }
  void setSceneView(const SceneView& view);
};

#endif
