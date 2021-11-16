#include "CGORenderer.h"

#include "Base.h"
#include "MemoryDebug.h"
#include "PyMOLGlobals.h"
#include "Rep.h"
#include "Setting.h"

struct RenderInfo;
struct Rep;
struct PyMOLGlobals;
struct CSetting;

unsigned CCGORenderer::pick_pass() const noexcept
{
  return info->pick->m_pass;
}

bool CGORendererInit(PyMOLGlobals* G)
{
  G->CGORenderer = new CCGORenderer();
  if (!G->CGORenderer) {
    return false;
  }
  G->CGORenderer->G = G;
  G->CGORenderer->isPicking = false;
  G->CGORenderer->alpha = 1.0F;
  return true;
}

void CGORendererFree(PyMOLGlobals* G)
{
  DeleteP(G->CGORenderer);
}
