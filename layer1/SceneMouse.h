#pragma once

#include "Picking.h"
#include "pymol/zstring_view.h"

struct Block;
struct PyMOLGlobals;

namespace pymol
{
struct CObject;
}

struct NamedPickContext {
  std::string name;
  int state;
};

struct NamedPicking {
  Pickable src;
  NamedPickContext context;
  NamedPicking(const Picking& pick);
};

void SceneClickObject(PyMOLGlobals* G, pymol::CObject* obj, const NamedPicking& LastPicked,
    int mode, pymol::zstring_view sel_mode_kw);
void SceneClickTransformObject(
    PyMOLGlobals* G, pymol::CObject* obj, const NamedPicking& LastPicked, int mode, bool is_single_click);
void SceneClickPickBond(PyMOLGlobals* G, int x, int y, int mode, const NamedPicking& LastPicked);

/**
 * Handles mouse clicks in the scene
 * @param block Scene block
 * @param button button pressed
 * @param x x coordinate
 * @param y y coordinate
 * @param mod modifier key
 * @param when time of event
 *
 * @pre Active rendering context
 */
void SceneClick(Block* block, int button, int x, int y, int mod, double when);

/**
 * Handles mouse release in the scene
 * @param block Scene block
 * @param button button released
 * @param x x coordinate
 * @param y y coordinate
 * @param mod modifier key
 * @param when time of event
 *
 * @pre Active rendering context
 */
void SceneRelease(Block* block, int button, int x, int y, int mod, double when);

/**
 * Handles mouse drag in the scene
 * @param block Scene block
 * @param button button released
 * @param x x coordinate
 * @param y y coordinate
 * @param mod modifier key
 * @param when time of event
 *
 * @pre Active rendering context
 */
void SceneDrag(Block* block, int x, int y, int mod, double when);

/**
 * Handles mouse drag in the scene
 * @param block Scene block
 * @param button button released
 * @param x x coordinate
 * @param y y coordinate
 * @param mod modifier key
 * @param when time of event
 *
 * @pre Active rendering context
 */
void SceneMouseMove(Block* block, int x, int y, int mod, double when);

/**
 * Deselect or toggle off entries
 * @param button button pressed
 * @param mod modifier key
 * @param mode button mode
 * @return if something has been deselected or toggled
 */

void SceneClickPickNothing(PyMOLGlobals* G, int button, int mod, int mode);

