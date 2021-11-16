#pragma once

struct RenderInfo;
struct Rep;
struct PyMOLGlobals;
struct CSetting;

struct CCGORenderer {
  PyMOLGlobals* G = nullptr;
  RenderInfo* info = nullptr;
  Rep* rep = nullptr;
  const float* color = nullptr;
  float alpha{};
  short sphere_quality{};
  bool isPicking{};
  unsigned pick_pass() const noexcept;
  bool use_shader{}; // OpenGL 1.4+, e.g., glEnableVertexAttribArray() (on) vs.
                     // glEnableClientState() (off)
  bool debug{};
  CSetting* set1 = nullptr;
  CSetting* set2 = nullptr;
};

bool CGORendererInit(PyMOLGlobals* G);
void CGORendererFree(PyMOLGlobals* G);
