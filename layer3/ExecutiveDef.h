#pragma once

#include "Ortho.h"
#include "ScrollBar.h"

class SpecRec;
class CGO;
struct _CObject;
struct CTracker;
struct _OVLexicon;
struct _OVOneToOne;
struct ExecutiveObjectOffset;



struct PanelRec {
  SpecRec *spec;
  unsigned nest_level;
  bool is_group = false;
  bool is_open = false;

  PanelRec(SpecRec* spec_, unsigned nest_level_)
      : spec(spec_)
      , nest_level(nest_level_)
  {
  }
};

struct ListMember{
  int list_id;
  int next;
};

enum class ExecutiveDragMode {
  Off,
  Visibility,
  Reorder,
  VisibilityWithCamera
};

enum class ExecutiveToggleMode {
  DeferVisibility,
  ImmediateVisibility,
  HoverActivate,
  CenterActivateDeactivatePrevious,
  ZoomActivateDeactivatePrevious,
  ZoomExclusiveActivate
};

struct CExecutive : public Block {
  SpecRec *Spec {};
  CTracker *Tracker {};
  int Width {}, Height {}, HowFarDown { 0 };
  int ScrollBarActive { 0 };
  int NSkip { 0 };
  ScrollBar m_ScrollBar;
  pymol::CObject *LastEdited { nullptr };
  ExecutiveDragMode DragMode = ExecutiveDragMode::Off;
  ExecutiveToggleMode ToggleMode = ExecutiveToggleMode::DeferVisibility;
  int Pressed { -1 }, Over { -1 }, LastOver {}, OldVisibility {}, PressedWhat {}, OverWhat {};
  SpecRec *LastChanged { nullptr }, *LastZoomed { nullptr }, *RecoverPressed { nullptr };
  int ReorderFlag { false };
  OrthoLineType ReorderLog {};
#ifndef GLUT_FULL_SCREEN
  // freeglut has glutLeaveFullScreen, no need to remember window dimensions
  int oldPX {}, oldPY {}, oldWidth {}, oldHeight {};
#endif
  int all_names_list_id {}, all_obj_list_id {}, all_sel_list_id {};
  OVLexicon *Lex {};
  OVOneToOne *Key {};
  bool ValidGroups { false };
  bool ValidSceneMembers { false };
  int ValidGridSlots {};

  std::vector<PanelRec> Panel{};

#ifdef _WEBGL
#endif
  int CaptureFlag {};
  int LastMotionCount {};
  CGO *selIndicatorsCGO { nullptr };
  int selectorTexturePosX { 0 }, selectorTexturePosY { 0 }, selectorTextureAllocatedSize { 0 }, selectorTextureSize { 0 };
  short selectorIsRound { 0 };

  // AtomInfoType::unique_id -> (object, atom-index)
  ExecutiveObjectOffset *m_eoo {}; // VLA of (object, atom-index)
  OVOneToOne *m_id2eoo {}; // unique_id -> m_eoo-index

  CExecutive(PyMOLGlobals * G) : Block(G), m_ScrollBar(G, false) {};

  int release(int button, int x, int y, int mod) override;
  int click(int button, int x, int y, int mod) override;
  int drag(int x, int y, int mod) override;
  void draw(CGO* orthoCGO) override;
  void reshape(int width, int height) override;
};
