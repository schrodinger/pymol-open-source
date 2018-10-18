
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

#include"Base.h"
#include"PyMOLObject.h"
#include"Ortho.h"
#include"View.h"
#include"ScrollBar.h"
#include<list>
#include<vector>

#ifdef PURE_OPENGL_ES_2
# define GLEW_EXT_gpu_shader4 false
# define GLEW_EXT_geometry_shader4 false
# ifdef _WEBGL
#  include <emscripten/val.h>
#  define GLEW_EXT_draw_buffers2 !emscripten::val::module_property("ONEBUFFER").as<bool>()
# else
#  define GLEW_EXT_draw_buffers2 false
# endif
#endif

#define TM3_IS_ONEBUF !GLEW_EXT_draw_buffers2

#define TRN_BKG 0x30
#define MAX_ANI_ELEM 300

struct DeferredMouse : public CDeferred {
  Block *block { nullptr };
  int button { 0 };
  int x { 0 };
  int y { 0 };
  int mod { 0 };
  double when { 0.0 };
  int mode_override { 0 };
  DeferredMouse(PyMOLGlobals * G) : CDeferred(G){}
};

struct DeferredImage : public CDeferred {
  int width { 0 };
  int height { 0 };
  std::string filename;
  int quiet { 0 };
  int antialias { 0 };
  float dpi { 0.0f };
  int entire_window { 0 };
  int format { 0 };
  DeferredImage(PyMOLGlobals * G) : CDeferred(G){}
};

struct DeferredRay : public CDeferred {
  int ray_width { 0 };
  int ray_height { 0 };
  int mode { 0 };
  float angle { 0.0f };
  float shift { 0.0f };
  int quiet { 0 };
  int show_timing { 0 };
  int antialias { 0 };
  DeferredRay(PyMOLGlobals * G) : CDeferred(G){}
};

typedef struct {
  int len;
  const char *name;
  int x1, y1, x2, y2, drawn;
} SceneElem;

typedef struct {
  unsigned char *data;
  int size;
  int width, height;
  int stereo;                   /* indicates data actually contains two back to back full-screen images */
  int needs_alpha_reset;        /* needs alpha reset */
} ImageType;

typedef struct {
  float unit_left, unit_right, unit_top, unit_bottom, unit_front, unit_back;
} SceneUnitContext;

typedef struct {
  int n_col;
  int n_row;
  int first_slot;
  int last_slot;
  float asp_adjust;
  int active;
  int size;
  int slot;
  int mode;
  GLint cur_view[4];
  GLint cur_viewport_size[2];
  SceneUnitContext context;     /* for whole-window display */
} GridInfo;

class CScene : public Block {
 public:
  std::list<CObject*> Obj, GadgetObjs, NonGadgetObjs;
  float RotMatrix[16];          /* WARNING: column major, as per OpenGL spec */
  float InvMatrix[16];          /* WARNING: column major, as per OpenGL spec */
  float PmvMatrix[16];
  float Scale;
  int Width, Height;
  int Button;
  int LastX, LastY;
  int StartX, StartY;
  int LastWinX, LastWinY;
  double LastClickTime;
  int LastButton, LastMod;
  int PossibleSingleClick;
  double LastReleaseTime;
  double SingleClickDelay;
  float ViewNormal[3], LinesNormal[3];
  float Pos[3], Origin[3];
  float H;
  float Front, Back, FrontSafe, BackSafe;
  float TextColor[3];
  double SweepTime;
  int DirtyFlag;
  int ChangedFlag;
  int CopyType, CopyNextFlag, CopyForced;
  int NFrame { 0 };
  int HasMovie { 0 };
  ImageType *Image { nullptr };
  int MovieOwnsImageFlag;
  int MovieFrameFlag;
  double LastRender, RenderTime, LastFrameTime, LastFrameAdjust;
  double LastSweep, LastSweepTime;
  float LastSweepX, LastSweepY;
  int RockFrame;
  Picking LastPicked;
  int StereoMode;
  OrthoLineType vendor, renderer, version;
  int SculptingFlag, SculptingSave;
  int RovingDirtyFlag;
  int RovingCleanupFlag;
  double RovingLastUpdate;
  int Threshold, ThresholdX, ThresholdY;
  CView *View { nullptr };
  float LastPickVertex[3], LastClickVertex[3];
  int LastPickVertexFlag;
  int LoopFlag;
  int LoopMod;
  BlockRect LoopRect;
  CViewElem ani_elem[MAX_ANI_ELEM + 1];
  int cur_ani_elem, n_ani_elem;
  int LastStateBuilt;
  int AnimationStartFlag;
  double AnimationStartTime;
  double AnimationLagTime;
  int AnimationStartFrame;
  double ApproxRenderTime;
  float VertexScale;
  float FogStart;
  float FogEnd;

  /* Scene Names */
  int ButtonsShown, ButtonDrag, ButtonMargin, ButtonsValid;
  int Over, Pressed, PressMode, HowFarDown, NSkip;
  int ScrollBarActive;
  int ReorderFlag;
  OrthoLineType ReorderLog;
  ScrollBar m_ScrollBar;
  char *SceneNameVLA { nullptr };
  SceneElem *SceneVLA { nullptr };
  int NScene;
  CGO *AlphaCGO { nullptr };

  int *SlotVLA { nullptr };

  int StencilValid, StencilParity;
  int ReinterpolateFlag;
  CObject *ReinterpolateObj { nullptr };
  CObject *MotionGrabbedObj { nullptr };

  short prev_no_z_rotation1, prev_no_z_rotation2;
  int orig_x_rotation, orig_y_rotation;

  std::vector<float> m_ModelViewMatrixStack;
  int m_ModelViewMatrixStackDepth { 0 };

  union {
    float ModelViewMatrix[16];
    float ModMatrix[16];          // old alias, deprecated
  };

  float ProjectionMatrix[16];
  int background_color_already_set;
  int do_not_clear;
  GridInfo grid;
  int last_grid_size;
  int n_texture_refreshes { 0 };
  CGO *offscreenCGO { nullptr };
  CGO *offscreenOIT_CGO { nullptr };
  CGO *offscreenOIT_CGO_copy { nullptr };
  void (*vp_prepareViewPortForStereo)(PyMOLGlobals *, CScene *, int, short, int, int, int, int, int);
  int vp_times, vp_x, vp_y, vp_owidth, vp_oheight, vp_stereo_mode;
  float vp_width_scale;
  Picking *pickVLA { nullptr };
  bool invPick; // if set, picking should be re-built

  CScene(PyMOLGlobals * G) : Block(G), m_ScrollBar(G, false) {}

  virtual int click(int button, int x, int y, int mod) override;
  virtual int release(int button, int x, int y, int mod) override;
  virtual int drag(int x, int y, int mod) override;
  virtual void draw(CGO* orthoCGO) override;
  virtual void reshape(int width, int height) override;

  // PYMOL-2561, PYMOL-2711 : This structure used to be calloc-ed, this replicates that
  void *operator new(size_t size) {
    void *mem = ::operator new(size);
    memset(mem, 0, size);
    return mem;
  }
};

#endif
