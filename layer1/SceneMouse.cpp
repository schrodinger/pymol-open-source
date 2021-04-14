#include "SceneMouse.h"
#include "ButMode.h"
#include "Editor.h"
#include "Executive.h"
#include "Matrix.h"
#include "Menu.h"
#include "ObjectDist.h"
#include "ObjectGadget.h"
#include "ObjectMolecule.h"
#include "ObjectSlice.h"
#include "P.h"
#include "PyMOLGlobals.h"
#include "PyMOLObject.h"
#include "SceneDef.h"
#include "ScenePicking.h"
#include "Selector.h"
#include "Seq.h"
#include "Wizard.h"

#include <glm/glm.hpp>

#define cDoubleTime 0.35

NamedPicking::NamedPicking(const Picking& pick)
  : src(pick.src)
{
  if (pick.context.object != nullptr) {
    context.name = pick.context.object->Name;
  }
  context.state = pick.context.state;
}

static void SceneNoteMouseInteraction(PyMOLGlobals* G)
{
  SceneAbortAnimation(G);
  if (SettingGet_b(G, NULL, NULL, cSetting_mouse_restart_movie_delay)) {
    SceneRestartFrameTimer(G);
  }
}

static int SceneLoopClick(Block* block, int button, int x, int y, int mod)
{
  PyMOLGlobals* G = block->m_G;
  CScene* I = G->Scene;
  I->LoopRect.left = x;
  I->LoopRect.top = y;
  I->LoopRect.right = x;
  I->LoopRect.bottom = y;
  I->LoopFlag = true;
  I->LoopMod = mod;
  OrthoSetLoopRect(G, true, &I->LoopRect);
  OrthoGrab(G, block);
  return 1;
}

static int SceneLoopDrag(Block* block, int x, int y, int mod)
{
  PyMOLGlobals* G = block->m_G;
  CScene* I = G->Scene;
  I->LoopRect.right = x;
  I->LoopRect.bottom = y;
  OrthoSetLoopRect(G, true, &I->LoopRect);
  return 1;
}

static int SceneLoopRelease(Block* block, int button, int x, int y, int mod)
{
  PyMOLGlobals* G = block->m_G;
  CScene* I = G->Scene;
  int tmp;
  int mode;
  mode = ButModeTranslate(G, button, I->LoopMod);

  if (I->LoopRect.top < I->LoopRect.bottom) {
    tmp = I->LoopRect.top;
    I->LoopRect.top = I->LoopRect.bottom;
    I->LoopRect.bottom = tmp;
  }
  if (I->LoopRect.right < I->LoopRect.left) {
    tmp = I->LoopRect.right;
    I->LoopRect.right = I->LoopRect.left;
    I->LoopRect.left = tmp;
  }
  OrthoSetLoopRect(G, false, &I->LoopRect);
  ExecutiveSelectRect(G, &I->LoopRect, mode);
  I->LoopFlag = false;
  OrthoUngrab(G);
  OrthoDirty(G);
  return 1;
}

void SceneClickButtonAddTo(PyMOLGlobals* G, pymol::CObject* obj,
    pymol::zstring_view selName, pymol::zstring_view buffer,
    pymol::zstring_view sel_mode_kw)
{
  CScene* I = G->Scene;
  if (SelectorIndexByName(G, selName.c_str()) >= 0) {
    auto buf2 = pymol::string_format(
        "(((%s) or %s(%s)) and not ((%s(%s)) and %s(%s)))", selName,
        sel_mode_kw, buffer, sel_mode_kw, buffer, sel_mode_kw, selName);
    SelectorCreate(G, selName.c_str(), buf2.c_str(), nullptr, false, nullptr);
    if (obj->type == cObjectMolecule) {
      if (SettingGetGlobal_i(G, cSetting_logging)) {
        auto objMol = (ObjectMolecule*) obj;
        auto atomSele = ObjectMoleculeGetAtomSeleLog(
            objMol, I->LastPicked.src.index, false);
        auto select = pymol::string_format(
            "(((%s) or %s(%s)) and not ((%s(%s)) and %s(%s)))", selName,
            sel_mode_kw, atomSele, sel_mode_kw, atomSele, sel_mode_kw, selName);
        auto pLogBuffer =
            pymol::string_format("cmd.select('%s',\"%s(%s)\",enable=1)",
                selName, sel_mode_kw, select);
        PLog(G, pLogBuffer, cPLog_pym);
      }
    }
  } else {
    auto buf2 = pymol::string_format("%s(%s)", sel_mode_kw, buffer);
    SelectorCreate(G, selName.c_str(), buf2.c_str(), NULL, false, NULL);
    if (obj->type == cObjectMolecule) {
      if (SettingGet<int>(G, cSetting_logging)) {
        auto objMol = (ObjectMolecule*) obj;
        auto atomSele = ObjectMoleculeGetAtomSeleLog(
            objMol, I->LastPicked.src.index, false);
        auto select = pymol::string_format(
            "cmd.select('%s',\"%s(%s)\")", selName, sel_mode_kw, atomSele);
        PLog(G, select, cPLog_pym);
      }
    }
  }
  if (SettingGet<bool>(G, cSetting_auto_hide_selections))
    ExecutiveHideSelections(G);
  if (SettingGet<bool>(G, cSetting_auto_show_selections))
    ExecutiveSetObjVisib(G, selName, 1, false);
  WizardDoSelect(G, selName.c_str(), I->LastPicked.context.state);
}

/**
 * Check for Double click
 * @param x cursor X position
 * @param Y cursor Y position
 * @param mod modifier key (Shift, ctrl, etc...)
 * @param when Time which this click happened
 * @return new button
 */
static int SceneClickCheckDoubleClick(
    PyMOLGlobals* G, int button, int x, int y, int mod, double when)
{
  auto I = G->Scene;
  // check for double click (within 0.35s and 10sq. pixels
  if (((ButModeCheckPossibleSingleClick(G, button, mod) || (!mod)) &&
          ((when - I->LastClickTime) < cDoubleTime))) {
    int dx = abs(I->LastWinX - x);
    int dy = abs(I->LastWinY - y);
    if ((dx < 10) && (dy < 10) && (I->LastButton == button)) {
      switch (button) {
      case P_GLUT_LEFT_BUTTON:
        button = P_GLUT_DOUBLE_LEFT;
        break;
      case P_GLUT_MIDDLE_BUTTON:
        button = P_GLUT_DOUBLE_MIDDLE;
        break;
      case P_GLUT_RIGHT_BUTTON:
        button = P_GLUT_DOUBLE_RIGHT;
        break;
      }
    }
  }
  return button;
}

/**
 * @param button which button pressed
 * @param x Cursor X position
 * @param Y Cursor Y position
 * @param mod modifier key (Shift, ctrl, etc...)
 * @return true if click was handled (resolved with action)
 */
static int SceneClickSceneButton(
    PyMOLGlobals* G, int button, int x, int y, int mod)
{
  auto I = G->Scene;
  for (std::size_t i = 0u; i < I->SceneVec.size(); i++) {
    auto& elem = I->SceneVec[i];
    if (elem.drawn && elem.rect.contains(x, y)) {
      switch (button) {
      case P_GLUT_LEFT_BUTTON: /* normal activate (with interpolation) */
        I->Pressed = i;
        I->Over = i;
        I->PressMode = 1;
        SceneDirty(G);
        return true;
      case P_GLUT_MIDDLE_BUTTON: /* rapid browse mode */
        I->Pressed = i;
        I->PressMode = 2;
        I->Over = i;

        {
          const char* cur_name =
              SettingGetGlobal_s(G, cSetting_scene_current_name);
          int animate = -1;
          if (mod & cOrthoCTRL)
            animate = 0;
          if (cur_name && elem.name != cur_name) {
            auto buffer = pymol::string_format(
                "cmd.scene('''%s''',animate=%d)", elem.name, animate);
            PParse(G, buffer);
            PFlush(G);
            PLog(G, buffer, cPLog_pym);
          }
        }
        return true;
      case P_GLUT_RIGHT_BUTTON: /* drag or menu... */
        I->Pressed = i;
        I->PressMode = 3;
        I->Over = i;
        return true;
      }
      break;
    }
  }
  return false;
}

/**
 * Handles object-related click events
 * @param obj clicked object
 * @param LastPicked last picked info
 * @param mode mouse button mode
 * @param sel_mode_kw Current selection mode operator for the mouse
 *
 */
void SceneClickObject(PyMOLGlobals* G, pymol::CObject* obj, const NamedPicking& LastPicked,
    int mode, pymol::zstring_view sel_mode_kw)
{
  std::string selName;
  switch (obj->type) {
  case cObjectMolecule: {
    if (Feedback(G, FB_Scene, FB_Results)) {
      auto buffer = obj->describeElement(LastPicked.src.index);
      PRINTF " You clicked %s", buffer.c_str() ENDF(G);
      OrthoRestorePrompt(G);
    }
    auto buffer =
        pymol::string_format("%s`%d", obj->Name, LastPicked.src.index + 1);
    switch (mode) { // Populate SelName/buffers
    case cButModeLB:
    case cButModeAddToLB:
      selName = "lb";
      break;
    case cButModeMB:
    case cButModeAddToMB:
      selName = "mb";
      break;
    case cButModeRB:
    case cButModeAddToRB:
      selName = "rb";
      break;
    case cButModeSeleSet:
    case cButModeSeleToggle:
      ExecutiveGetActiveSeleName(
          G, selName, true, SettingGet<int>(G, cSetting_logging));
      break;
    case cButModeDragMol: {
      auto objMol = (ObjectMolecule*) obj;
      auto atomSele =
          ObjectMoleculeGetAtomSeleLog(objMol, LastPicked.src.index, false);
      auto pLogBuffer =
          pymol::string_format("cmd.drag(\"bymol (%s)\")", atomSele);
      PParse(G, pLogBuffer);
      PLog(G, buffer, cPLog_pym);
    } break;
    case cButModeDragObj: {
      auto objMol = (ObjectMolecule*) obj;
      auto atomSele =
          ObjectMoleculeGetAtomSeleLog(objMol, LastPicked.src.index, false);
      auto pLogBuffer =
          pymol::string_format("cmd.drag(\"byobject (%s)\")", atomSele);
      PParse(G, pLogBuffer);
      PLog(G, pLogBuffer, cPLog_pym);
    } break;
    case cButModeOrigAt:
      SceneNoteMouseInteraction(G);
      {
        float v1[3];

        if (ObjectMoleculeGetAtomTxfVertex((ObjectMolecule*) obj,
                LastPicked.context.state, LastPicked.src.index, v1)) {
          EditorFavorOrigin(G, v1);
          ExecutiveOrigin(G, NULL, true, NULL, v1, 0);
        }
      }
      if (obj->type == cObjectMolecule) {
        if (SettingGetGlobal_i(G, cSetting_logging)) {
          auto objMol = (ObjectMolecule*) obj;
          auto atomSele =
              ObjectMoleculeGetAtomSeleLog(objMol, LastPicked.src.index, false);
          auto pLogBuffer =
              pymol::string_format("cmd.origin(\"%s\")", atomSele);
          PLog(G, pLogBuffer, cPLog_pym);
        }
        if (Feedback(G, FB_Scene, FB_Results)) {
          auto buffer = obj->describeElement(LastPicked.src.index);
          PRINTF " You clicked %s", buffer.c_str() ENDF(G);
          OrthoRestorePrompt(G);
        }
      }
      PRINTFB(G, FB_Scene, FB_Actions)
      " Scene: Origin set.\n" ENDFB(G);
      break;
    case cButModeCent:
      SceneNoteMouseInteraction(G);
      {
        float v1[3];

        if (ObjectMoleculeGetAtomTxfVertex((ObjectMolecule*) obj,
                LastPicked.context.state, LastPicked.src.index, v1)) {
          ExecutiveCenter(G, nullptr, 0, true, -1, v1, true);
        }
      }

      if (SettingGet<int>(G, cSetting_logging)) {
        auto objMol = (ObjectMolecule*) obj;
        auto atomSele =
            ObjectMoleculeGetAtomSeleLog(objMol, LastPicked.src.index, false);
        auto pLogBuffer = pymol::string_format("cmd.center(\"%s\",state=%d)",
            atomSele, LastPicked.context.state + 1);
        PLog(G, pLogBuffer, cPLog_pym);
      }
      break;
    }
    switch (mode) {
    case cButModeSimpleClick:
      break;
    case cButModeLB:
    case cButModeMB:
    case cButModeRB:
    case cButModeSeleSet: {
      auto buf2 = pymol::string_format("(%s(%s))", sel_mode_kw, buffer);
      SelectorCreate(G, selName.c_str(), buf2.c_str(), nullptr, false, nullptr);
    }
      if (SettingGet<bool>(G, cSetting_auto_hide_selections))
        ExecutiveHideSelections(G);
      if (SettingGet<bool>(G, cSetting_auto_show_selections))
        ExecutiveSetObjVisib(G, selName, 1, false);
      if (obj->type == cObjectMolecule) {
        if (SettingGet<int>(G, cSetting_logging)) {
          auto objMol = (ObjectMolecule*) obj;
          auto buf1 =
              ObjectMoleculeGetAtomSeleLog(objMol, LastPicked.src.index, false);
          auto pLogBuffer =
              pymol::string_format("cmd.select('%s',\"%s(%s)\",enable=1)",
                  selName, sel_mode_kw, buf1);
          PLog(G, pLogBuffer, cPLog_pym);
        }
      }
      WizardDoSelect(G, selName.c_str(), LastPicked.context.state);
      break;
    case cButModeAddToLB:
    case cButModeAddToMB:
    case cButModeAddToRB:
    case cButModeSeleToggle:
      SceneClickButtonAddTo(G, obj, selName, buffer.c_str(), sel_mode_kw);
      break;
    }
  }
  case cObjectGadget:
    break;
  default:
    EditorInactivate(G);
    break;
  }
}

void SceneClickTransformObject(
    PyMOLGlobals* G, pymol::CObject* obj, const NamedPicking& LastPicked, int mode, bool is_single_click)
{
  auto I = G->Scene;
  switch (obj->type) {
  case cObjectMolecule:
    switch (mode) {
    case cButModeMenu: {
      ObjectMolecule* objMol = (ObjectMolecule*) obj;
      int active_sele = ExecutiveGetActiveSele(G);
      if (active_sele && SelectorIsMember(G,
                             objMol->AtomInfo[LastPicked.src.index].selEntry,
                             active_sele)) {
        /* user clicked on a selected atom */
        ObjectNameType name;
        ExecutiveGetActiveSeleName(
            G, name, false, SettingGetGlobal_i(G, cSetting_logging));
        MenuActivate2Arg(G, I->LastWinX, I->LastWinY + 20, /* selection menu */
            I->LastWinX, I->LastWinY, is_single_click, "pick_sele", name, name);
      } else {
        /* user clicked on an atom not in a selection */
        auto buffer = obj->describeElement(LastPicked.src.index);
        auto atomSele = ObjectMoleculeGetAtomSeleLog(
            (ObjectMolecule*) obj, LastPicked.src.index, false);
        MenuActivate2Arg(G, I->LastWinX, I->LastWinY + 20, I->LastWinX,
            I->LastWinY, is_single_click, "pick_menu", buffer.c_str(),
            atomSele.c_str());
      }
    } break;
    case cButModePickAtom1:
      if (obj && obj->type == cObjectMolecule) {
        if (Feedback(G, FB_Scene, FB_Results)) {
          auto buffer = obj->describeElement(LastPicked.src.index);
          PRINTF " You clicked %s -> (%s)\n", buffer.c_str(), cEditorSele1 ENDF(G);
        }
        if (SettingGet<int>(G, cSetting_logging)) {
          auto objMol = (ObjectMolecule*) obj;
          auto atomSele = ObjectMoleculeGetAtomSeleLog(
              objMol, LastPicked.src.index, false);
          auto pLogBuffer =
              pymol::string_format("cmd.edit(\"%s\",pkresi=1)", atomSele);
          PLog(G, pLogBuffer, cPLog_pym);
        }
        OrthoRestorePrompt(G);
        auto buffer = pymol::string_format(
            "%s`%d", obj->Name, LastPicked.src.index + 1);
        EditorInactivate(G);
        SelectorCreate(G, cEditorSele1, buffer.c_str(), nullptr, true, nullptr);
        EditorActivate(G, SettingGet<int>(G, cSetting_state) - 1, false);
        if (EditorActive(G)) {
          EditorDefineExtraPks(G);
        }
        WizardDoPick(G, 0, LastPicked.context.state);
      }
      break;
    case cButModePickAtom:
      if (obj && obj->type == cObjectMolecule) {
        auto buffer = obj->describeElement(LastPicked.src.index);
        if (EditorIsBondMode(G)
            /* &&!(EditorIsAnActiveObject(G,(ObjectMolecule*)obj)) */
        ) {
          EditorInactivate(G);
          EditorLogState(G, false);
        }
        if ((!EditorIsBondMode(G)) &&
            EditorDeselectIfSelected(
                G, (ObjectMolecule*) obj, LastPicked.src.index, true)) {
          PRINTF " You unpicked %s.", buffer.c_str() ENDF(G);
          if (EditorActive(G))
            EditorDefineExtraPks(G);
          EditorLogState(G, false);
        } else {
          if (EditorIsBondMode(G) &&
              EditorDeselectIfSelected(
                  G, (ObjectMolecule*) obj, LastPicked.src.index, false)) {
            EditorInactivate(G);
          }
          WordType name;
          EditorGetNextMultiatom(G, name);

          PRINTFB(G, FB_Scene, FB_Results)
          " You clicked %s -> (%s)\n", buffer.c_str(), name ENDFB(G);
          /* TODO: logging */

          auto buffer = pymol::string_format(
              "%s`%d", obj->Name, LastPicked.src.index + 1);
          ExecutiveDelete(G, name);
          SelectorCreate(G, name, buffer.c_str(), nullptr, true, nullptr);
          EditorActivate(G, SettingGet<int>(G, cSetting_state) - 1, false);
          if (EditorActive(G)) {
            EditorDefineExtraPks(G);
          }
          EditorLogState(G, false);
          WizardDoPick(G, 0, LastPicked.context.state);
        }
      }
      break;
    }
    break;
  case cObjectGadget:
    break;
  default:
    EditorInactivate(G);
    break;
  }
}

/**
 * Handle logic for objects with pickable bonds
 * @param x cursor X
 * @param y cursor Y
 * @param mode mouse button mode
 * @param lastPicked last picked info
 */

void SceneClickPickBond(PyMOLGlobals* G, int x, int y, int mode, const NamedPicking& LastPicked)
{
  auto I = G->Scene;
  auto obj = ExecutiveFindObject<ObjectMolecule>(G, LastPicked.context.name);
  EditorInactivate(G);
  if (!obj) {
    return;
  }
  if (Feedback(G, FB_Scene, FB_Results)) {
    auto buffer = obj->describeElement(I->LastPicked.src.index);
    PRINTF " You clicked %s -> (%s)", buffer.c_str(), cEditorSele1 ENDF(G);
    OrthoRestorePrompt(G);
  }

  /*        ObjectMoleculeChooseBondDir(objMol,I->LastPicked.bond,
     &I->LastPicked.src.index,&atIndex); */

  {
    auto buffer = pymol::string_format("%s`%d", obj->Name, I->LastPicked.src.index + 1);
    SelectorCreate(G, cEditorSele1, buffer.c_str(), NULL, true, NULL);
  }

  if (LastPicked.src.bond >= 0) {
    int atIndex = obj->Bond[LastPicked.src.bond].index[0];
    if (atIndex == LastPicked.src.index)
      atIndex = obj->Bond[LastPicked.src.bond].index[1];
    if (Feedback(G, FB_Scene, FB_Results)) {
      auto buffer = obj->describeElement(atIndex);
      PRINTF " You clicked %s -> (%s)", buffer.c_str(), cEditorSele2 ENDF(G);
      OrthoRestorePrompt(G);
    }

    if (SettingGetGlobal_i(G, cSetting_logging)) {
      auto buf1 = ObjectMoleculeGetAtomSeleLog(
          obj, LastPicked.src.index, false);
      auto buf2 = ObjectMoleculeGetAtomSeleLog(obj, atIndex, false);
      auto buffer = pymol::string_format("cmd.edit(\"%s\",\"%s\")", buf1, buf2);
      PLog(G, buffer, cPLog_pym);
    }
    auto buffer = pymol::string_format("%s`%d", obj->Name, atIndex + 1);
    SelectorCreate(G, cEditorSele2, buffer.c_str(), NULL, true, NULL);
    EditorActivate(G, SettingGetGlobal_i(G, cSetting_state) - 1, true);

    if (mode == cButModePkTorBnd) {
      /* get ready to drag */
      SceneDontCopyNext(G);
      EditorPrepareDrag(G, obj, -1, LastPicked.src.index,
          SettingGetGlobal_i(G, cSetting_state) - 1, mode);
      I->SculptingFlag = 1;
      I->SculptingSave =
          obj->AtomInfo[LastPicked.src.index].protekted;
      obj->AtomInfo[LastPicked.src.index].protekted = 2;
    }
    WizardDoPick(G, 1, LastPicked.context.state);
  } else {
    WizardDoPick(G, 0, LastPicked.context.state);
  }
  if (SettingGetGlobal_b(G, cSetting_auto_hide_selections))
    ExecutiveHideSelections(G);
}

/**
 * Handle logic for picking nothing
 * @param button which button pressed
 * @param mod modifier key (Shift, ctrl, etc...)
 * @param mode mouse button mode
 */

void SceneClickPickNothing(PyMOLGlobals* G, int button, int mod, int mode)
{
  auto I = G->Scene;
  switch (mode) {
  case cButModeSeleSet: {
    ObjectNameType name;
    if (ExecutiveGetActiveSeleName(
            G, name, false, SettingGet<int>(G, cSetting_logging))) {
      SelectorCreate(G, name, "none", nullptr, true, nullptr);
      if (SettingGet<int>(G, cSetting_logging)) {
        auto buf2 = pymol::string_format("cmd.select('%s','none')\n", name);
        PLog(G, buf2, cPLog_no_flush);
      }
      SeqDirty(G);
    }
  }
  case cButModeSeleToggle: {
    ObjectNameType name;

    if (ExecutiveGetActiveSeleName(
            G, name, false, SettingGet<int>(G, cSetting_logging))) {
      ExecutiveSetObjVisib(G, name, 0, false);
      if (SettingGet<int>(G, cSetting_logging)) {
        auto buf2 = pymol::string_format("cmd.disable('%s')\n", name);
        PLog(G, buf2, cPLog_no_flush);
      }
    }
  } break;
  case cButModeSimpleClick:
    PyMOL_SetClickReady(G->PyMOL, "", -1, button, mod, I->LastWinX,
        I->Height - (I->LastWinY + 1), nullptr, 0);
    break;
  }
  PRINTFB(G, FB_Scene, FB_Blather)
  " %s: no atom found nearby.\n", __func__ ENDFB(G);
  /* this here to prevent display weirdness after
     an unsuccessful picking pass... not sure it helps though */
  SceneInvalidate(G);
  OrthoRestorePrompt(G);
}

static int get_stereo_x(int x, int* last_x, int width, int* click_side)
{
  int width_2 = width / 2;
  int width_3 = width / 3;
  if (!last_x) {
    if (x > width_2) {
      x -= width_2;
      if (click_side)
        *click_side = 1; /* right */
    } else {
      if (click_side)
        *click_side = -1; /* left */
    }
  } else {
    if ((x - (*last_x)) > width_3) {
      x -= width_2;
      if (click_side)
        *click_side = 1; /* right */
    } else if (((*last_x) - x) > width_3) {
      x += width_2;
      if (click_side)
        *click_side = 1; /* right */
    } else {
      if (click_side)
        *click_side = -1; /* left */
    }
  }
  return x;
}

/*========================================================================*/
static int SceneClick(
    Block* block, int button, int x, int y, int mod, double when)
{
  PyMOLGlobals* G = block->m_G;
  CScene* I = G->Scene;
  const char* sel_mode_kw = "";
  int is_single_click =
      ((button == P_GLUT_SINGLE_LEFT) || (button == P_GLUT_SINGLE_MIDDLE) ||
          (button == P_GLUT_SINGLE_RIGHT));

  if (!is_single_click) {
    int click_handled = false;

    if (I->ButtonsShown) {
      /* check & handle a click on the scrollbar */
      if (I->ScrollBarActive) {
        if ((x - I->rect.left) < (SceneScrollBarWidth + SceneScrollBarMargin)) {
          click_handled = true;
          I->m_ScrollBar.click(button, x, y, mod);
        }
      }
      if (!click_handled) {
        for (const auto& elem : I->SceneVec) {
          if (elem.drawn && elem.rect.contains(x, y)) {
            click_handled = true;
            break;
          }
        }
      }
    }

    if (!click_handled) {
      button = SceneClickCheckDoubleClick(G, button, x, y, mod, when);
    } // end not click handled

    if (ButModeCheckPossibleSingleClick(G, button, mod) || (!mod)) {
      I->PossibleSingleClick = 1;
    } else {
      const char* but_mode_name =
          SettingGetGlobal_s(G, cSetting_button_mode_name);
      if (but_mode_name && but_mode_name[0] == '1') {
        I->PossibleSingleClick = 1;
      } else {
        I->PossibleSingleClick = 0;
      }
    }
  } // end not single-click

  int click_handled = false;
  int click_side = 0;

  I->LastWinX = x;
  I->LastWinY = y;
  I->LastClickTime = when;
  I->LastButton = button;
  I->LastMod = mod;
  I->Threshold = 0;

  SceneGetCenter(G, I->LastClickVertex);
  {
    float vScale = SceneGetExactScreenVertexScale(G, I->LastClickVertex);
    float v[3];
    v[0] = -(I->Width / 2 - (x - I->rect.left)) * vScale;
    v[1] = -(I->Height / 2 - (y - I->rect.bottom)) * vScale;
    v[2] = 0;
    MatrixInvTransformC44fAs33f3f(I->m_view.m_rotMatrix, v, v);
    add3f(v, I->LastClickVertex, I->LastClickVertex);
  }

  if (I->ButtonsShown) {

    if (I->ScrollBarActive) {
      if ((x - I->rect.left) < (SceneScrollBarWidth + SceneScrollBarMargin)) {
        click_handled = true;
        I->m_ScrollBar.click(button, x, y, mod);
      }
    }
    if (!click_handled) {
      click_handled = SceneClickSceneButton(G, button, x, y, mod);
    }
  }
  if (!click_handled) {
    int mode = ButModeTranslate(
        G, button, mod); /* trying to work around something... */

    I->Button = button;
    I->SculptingSave = 0;
    switch (mode) {
    case cButModeScaleSlabExpand:
      SceneNoteMouseInteraction(G);
      SceneClip(G, 5,
          1.0F + (0.2 * SettingGetGlobal_f(G, cSetting_mouse_wheel_scale)),
          NULL, 0);
      break;
    case cButModeScaleSlabShrink:
      SceneNoteMouseInteraction(G);
      SceneClip(G, 5,
          1.0F - (0.2 * SettingGetGlobal_f(G, cSetting_mouse_wheel_scale)),
          NULL, 0);
      break;
    case cButModeMoveSlabForward:
      SceneNoteMouseInteraction(G);
      {
        float old_front = I->m_view.m_clip.m_front;
        float old_back = I->m_view.m_clip.m_back;
        float old_origin = -I->m_view.m_pos[2];
        SceneClip(G, 6,
            0.1F * SettingGetGlobal_f(G, cSetting_mouse_wheel_scale), NULL, 0);
        SceneDoRoving(G, old_front, old_back, old_origin, true, false);
      }
      break;
    case cButModeMoveSlabBackward:
      SceneNoteMouseInteraction(G);
      {
        float old_front = I->m_view.m_clip.m_front;
        float old_back = I->m_view.m_clip.m_back;
        float old_origin = -I->m_view.m_pos[2];

        SceneClip(G, 6,
            -0.1F * SettingGetGlobal_f(G, cSetting_mouse_wheel_scale), NULL, 0);
        SceneDoRoving(G, old_front, old_back, old_origin, true, false);
      }
      break;
    case cButModeZoomForward:
      SceneNoteMouseInteraction(G);
      {
        float factor =
            -((I->m_view.m_clipSafe.m_front + I->m_view.m_clipSafe.m_back) /
                2) *
            0.1 * SettingGetGlobal_f(G, cSetting_mouse_wheel_scale);
        if (factor <= 0.0F) {
          I->m_view.m_pos[2] += factor;
          I->m_view.m_clip.m_front -= factor;
          I->m_view.m_clip.m_back -= factor;
          UpdateFrontBackSafe(I);
        }
      }
      break;
    case cButModeZoomBackward:
      SceneNoteMouseInteraction(G);
      {
        float factor =
            ((I->m_view.m_clipSafe.m_front + I->m_view.m_clipSafe.m_back) / 2) *
            0.1F * SettingGetGlobal_f(G, cSetting_mouse_wheel_scale);
        if (factor >= 0.0F) {
          I->m_view.m_pos[2] += factor;
          I->m_view.m_clip.m_front -= factor;
          I->m_view.m_clip.m_back -= factor;
          UpdateFrontBackSafe(I);
        }
      }
      break;
    case cButModeMoveSlabAndZoomForward:
      SceneNoteMouseInteraction(G);
      {
        float old_front = I->m_view.m_clip.m_front;
        float old_back = I->m_view.m_clip.m_back;
        float old_origin = -I->m_view.m_pos[2];
        SceneClip(G, 6,
            0.1F * SettingGetGlobal_f(G, cSetting_mouse_wheel_scale), NULL, 0);
        SceneDoRoving(G, old_front, old_back, old_origin, true, true);
      }
      break;
    case cButModeMoveSlabAndZoomBackward:
      SceneNoteMouseInteraction(G);
      {
        float old_front = I->m_view.m_clip.m_front;
        float old_back = I->m_view.m_clip.m_back;
        float old_origin = -I->m_view.m_pos[2];
        SceneClip(G, 6,
            -0.1F * SettingGetGlobal_f(G, cSetting_mouse_wheel_scale), NULL, 0);
        SceneDoRoving(G, old_front, old_back, old_origin, true, true);
      }
      break;
    case cButModeRectAdd: /* deprecated */
    case cButModeRectSub: /* deprecated */
    case cButModeRect:    /* deprecated */
    case cButModeSeleAddBox:
    case cButModeSeleSetBox:
    case cButModeSeleSubBox:
      return SceneLoopClick(block, button, x, y, mod);
      break;
    case cButModeRotDrag:
    case cButModeMovDrag:
    case cButModeMovDragZ:
      SceneNoteMouseInteraction(G);
      SceneDontCopyNext(G);

      y = y - I->margin.bottom;
      x = x - I->margin.left;

      if (stereo_via_adjacent_array(I->StereoMode))
        x = get_stereo_x(x, NULL, I->Width, NULL);

      I->LastX = x;
      I->LastY = y;
      EditorReadyDrag(G, SettingGetGlobal_i(G, cSetting_state) - 1);

      if (EditorDraggingObjectMatrix(G)) {
        auto obj = EditorDragObject(G);
        if (obj) {
          if (SettingGetGlobal_b(G, cSetting_movie_auto_store)) {
            ObjectTranslateTTT(obj, NULL, true);
            I->MotionGrabbedObj = obj;
            obj->Grabbed = true;
            if (SettingGetGlobal_i(G, cSetting_movie_auto_interpolate)) {
              I->ReinterpolateFlag = true;
              I->ReinterpolateObj = obj;
            }
          } else {
            ObjectTranslateTTT(obj, NULL, false);
          }
        }
      }
      break;
    case cButModeRotXYZ:
    case cButModeTransXY:
    case cButModeTransZ:
    case cButModeClipNF:
    case cButModeClipN:
    case cButModeClipF:
    case cButModeRotZ:
    case cButModeInvRotZ:
    case cButModeRotL:
    case cButModeMovL:
    case cButModeMvzL:
      SceneNoteMouseInteraction(G);
      SceneDontCopyNext(G);

      y = y - I->margin.bottom;
      x = x - I->margin.left;

      if (stereo_via_adjacent_array(I->StereoMode))
        x = get_stereo_x(x, NULL, I->Width, NULL);

      I->LastX = x;
      I->LastY = y;
      break;
    case cButModePickAtom1:
    case cButModePickAtom:
    case cButModeMenu:
      if (stereo_via_adjacent_array(I->StereoMode))
        x = get_stereo_x(x, NULL, I->Width, NULL);

      if (SceneDoXYPick(G, x, y, click_side)) {
        auto obj = I->LastPicked.context.object;
        y = y - I->margin.bottom;
        x = x - I->margin.left;
        I->LastX = x;
        I->LastY = y;
        SceneClickTransformObject(G, obj, I->LastPicked, mode, is_single_click);
      } else { /* no atom picked */
        switch (mode) {
        case cButModeMenu:

          MenuActivate3fv(G, I->LastWinX, I->LastWinY, I->LastWinX, I->LastWinY,
              is_single_click, "main_menu", I->LastClickVertex);
          break;
        default:
          EditorInactivate(G);
          if (SettingGetGlobal_i(G, cSetting_logging)) {
            PLog(G, "cmd.edit()", cPLog_pym);
          }
          break;
        }
      }
      SceneDirty(G);
      break;
    case cButModePickBond:
    case cButModePkTorBnd:
      if (stereo_via_adjacent_array(I->StereoMode))
        x = get_stereo_x(x, NULL, I->Width, &click_side);

      if (SceneDoXYPick(G, x, y, click_side)) {
        y = y - I->margin.bottom;
        x = x - I->margin.left;
        I->LastX = x;
        I->LastY = y;

        if (mode == cButModePkTorBnd) {
          I->Threshold = 3;
          I->ThresholdX = x;
          I->ThresholdY = y;
        }
        SceneClickPickBond(G, x, y, mode, I->LastPicked);
      } else {
        EditorInactivate(G);
        EditorLogState(G, false);
      }
      SceneInvalidate(G);
      break;
    case cButModeRotObj:
    case cButModeMovObj:
    case cButModeMovObjZ:
    case cButModeRotView:
    case cButModeMovView:
    case cButModeMovViewZ:
    case cButModeRotFrag:
    case cButModeMovFrag:
    case cButModeMovFragZ:
    case cButModeTorFrag:
    case cButModeMoveAtom:
    case cButModeMoveAtomZ:
      if (stereo_via_adjacent_array(I->StereoMode))
        x = get_stereo_x(x, NULL, I->Width, &click_side);

      if (SceneDoXYPick(G, x, y, click_side)) {
        auto obj = I->LastPicked.context.object;
        y = y - I->margin.bottom;
        x = x - I->margin.left;
        I->LastX = x;
        I->LastY = y;
        switch (obj->type) {
        case cObjectMolecule:

          if (I->LastPicked.src.bond == cPickableLabel) {
            /* if user picks a label with move object/move fragment,
               then move the object/fragment, not the label */

            switch (mode) {
            case cButModeRotObj:
            case cButModeMovObj:
            case cButModeMovObjZ:
            case cButModeRotFrag:
            case cButModeMovFrag:
            case cButModeMovFragZ:
            case cButModeMovViewZ:
            case cButModeRotView:
            case cButModeMovView:
              I->LastPicked.src.bond = cPickableAtom;
              break;
            }
          }

          switch (mode) {
          case cButModeMovViewZ:
          case cButModeRotView:
          case cButModeMovView: {
            if (SettingGetGlobal_b(G, cSetting_movie_auto_store)) {
              ObjectTranslateTTT(obj, NULL, true);
              I->MotionGrabbedObj = obj;
              obj->Grabbed = true;
              if (SettingGetGlobal_i(G, cSetting_movie_auto_interpolate)) {
                I->ReinterpolateFlag = true;
                I->ReinterpolateObj = obj;
              }
            } else {
              ObjectTranslateTTT(obj, NULL, false);
            }
          } break;
          }

          if (I->LastPicked.src.bond >= cPickableAtom) {
            if (Feedback(G, FB_Scene, FB_Results)) {
              auto buffer = obj->describeElement(I->LastPicked.src.index);
              PRINTF " You clicked %s", buffer.c_str() ENDF(G);
              OrthoRestorePrompt(G);
            }
          }
          {
            auto objMol = (ObjectMolecule*) obj;
            EditorPrepareDrag(G, obj, -1, I->LastPicked.src.index,
                SettingGetGlobal_i(G, cSetting_state) - 1, mode);

            if (I->LastPicked.src.bond >= cPickableAtom) {
              I->SculptingFlag = 1;
              I->SculptingSave =
                  objMol->AtomInfo[I->LastPicked.src.index].protekted;
              objMol->AtomInfo[I->LastPicked.src.index].protekted = 2;
            }
          }
          break;
        case cObjectSlice:
          if (ObjectSliceGetVertex((ObjectSlice*) obj, I->LastPicked.src.index,
                  I->LastPicked.src.bond, I->LastPickVertex)) {
            I->LastPickVertexFlag = true;
          }
          break;
        case cObjectMeasurement:
          break;
        case cObjectGadget:
          break;
        default:
          EditorInactivate(G);
          break;
        }
      }
      break;

    case cButModeSeleSet:
    case cButModeSeleToggle:
      sel_mode_kw = SceneGetSeleModeKeyword(G);

      /* intentional pass through */

    case cButModeLB:
    case cButModeMB:
    case cButModeRB:
    case cButModeAddToLB:
    case cButModeAddToMB:
    case cButModeAddToRB:
    case cButModeSimpleClick:
    case cButModeOrigAt:
    case cButModeCent:
    case cButModeDragMol:
    case cButModeDragObj:
      if (stereo_via_adjacent_array(I->StereoMode))
        x = get_stereo_x(x, NULL, I->Width, &click_side);

      if (SceneDoXYPick(G, x, y, click_side)) {
        auto obj = I->LastPicked.context.object;

        assert(I->LastPicked.src.bond != cPickableNoPick);

        if (mode == cButModeSimpleClick) {
          float pos_store[3];
          const float* pos = nullptr;
          int const index = I->LastPicked.src.index; /* 1-based */
          int const state = I->LastPicked.context.state;
          auto objmol = dynamic_cast<ObjectMolecule const*>(obj);

          if (objmol &&
              ObjectMoleculeGetAtomTxfVertex(objmol, state, index, pos_store))
            pos = pos_store;

          PyMOL_SetClickReady(G->PyMOL, obj->Name, index, button, mod,
              I->LastWinX, I->Height - (I->LastWinY + 1), pos,
              state + 1 /* 1-based */, I->LastPicked.src.bond);
        }
        SceneClickObject(G, obj, I->LastPicked, mode, sel_mode_kw);
      } else { // Nothing picked
        SceneClickPickNothing(G, button, mod, mode);
      }
    }

    I->StartX = I->LastX;
    I->StartY = I->LastY;
  }
  return (1);
}

/*========================================================================*/
static int SceneRelease(
    Block* block, int button, int x, int y, int mod, double when)
{
  PyMOLGlobals* G = block->m_G;
  CScene* I = G->Scene;
  int release_handled = false;
  if (I->ButtonsShown && I->PressMode) {
    if (I->ScrollBarActive) {
      if ((x - I->rect.left) < (SceneScrollBarWidth + SceneScrollBarMargin)) {
        I->m_ScrollBar.release(button, x, y, mod);
        release_handled = true;
      }
    }
    if (!release_handled) {
      int ungrab = true;
      if (I->PressMode) {
        I->Over = -1;
        auto elem = &I->SceneVec.front();
        for (std::size_t i = 0u; i < I->SceneVec.size(); i++) {
          auto& elem_ = I->SceneVec[i];
          if (elem_.drawn && elem_.rect.contains(x, y)) {
            elem = &elem_;
            I->Over = i;
            break;
          }
        }
        if (I->Over >= 0) {
          release_handled = true;
          switch (I->PressMode) {
          case 1:
            if (I->Over == I->Pressed) {
              auto buffer = pymol::string_format("cmd.scene('''%s''')", elem->name);
              PParse(G, buffer);
              PFlush(G);
              PLog(G, buffer, cPLog_pym);
            }
            break;
          case 2: {
            const char* cur_name =
                SettingGetGlobal_s(G, cSetting_scene_current_name);
            if (cur_name && elem->name != cur_name) {
              auto buffer = pymol::string_format("cmd.scene('''%s''')", elem->name);
              PParse(G, buffer);
              PFlush(G);
              PLog(G, buffer, cPLog_pym);
            }
          } break;
          case 3:
            if (I->Pressed == I->Over) {
              Block* block = MenuActivate1Arg(G, I->LastWinX,
                  I->LastWinY + 20, /* scene menu */
                  I->LastWinX, I->LastWinY, true, "scene_menu", elem->name.c_str());
              if (block)
                block->drag(x, y, mod);
              ungrab = false;
            }
            break;
          }
        }
      }
      I->LastPickVertexFlag = false;
      I->Pressed = -1;
      I->Over = -1;
      I->PressMode = 0;
      if (ungrab)
        OrthoUngrab(G);
    }
  }
  if (!release_handled) {
    ObjectMolecule* obj;
    I->LastReleaseTime = when;
    if (I->PossibleSingleClick == 1) {
      double slowest_single_click = 0.25F;
      double diff = when - I->LastClickTime;

      slowest_single_click += I->ApproxRenderTime;

      if ((diff < 0.0) || (diff > slowest_single_click))
        I->PossibleSingleClick = 0;
      else {
        int but = -1;
        I->PossibleSingleClick = 2;
        I->SingleClickDelay = 0.15;

        switch (I->LastButton) {
        case P_GLUT_LEFT_BUTTON:
          but = P_GLUT_DOUBLE_LEFT;
          break;
        case P_GLUT_MIDDLE_BUTTON:
          but = P_GLUT_DOUBLE_MIDDLE;
          break;
        case P_GLUT_RIGHT_BUTTON:
          but = P_GLUT_DOUBLE_RIGHT;
          break;
        }
        if (but > 0) {
          int mode = ButModeTranslate(G, but, mod);
          if (mode == cButModeNone)
            I->SingleClickDelay =
                0.0; /* no double-click set? force immediate single click */
        }
      }
    }
    if (I->LoopFlag) {
      I->PossibleSingleClick = 0;
      return SceneLoopRelease(block, button, x, y, mod);
    }
    OrthoUngrab(G);
    I->LoopFlag = false;
    if (I->SculptingFlag) {
      /* SettingSet(G,cSetting_sculpting,1); */
      obj = (ObjectMolecule*) I->LastPicked.context.object;
      if (obj) {
        obj->AtomInfo[I->LastPicked.src.index].protekted = I->SculptingSave;
      }
      I->SculptingFlag = 0;
    }
  }
  if (I->ReinterpolateFlag && I->ReinterpolateObj) {
    if (ExecutiveValidateObjectPtr(G, I->ReinterpolateObj, 0)) {
      ObjectMotionReinterpolate(I->ReinterpolateObj);
    }
    I->ReinterpolateFlag = true;
    I->ReinterpolateObj = NULL;
  }
  if (I->MotionGrabbedObj) {
    if (ExecutiveValidateObjectPtr(G, I->MotionGrabbedObj, 0)) {
      I->MotionGrabbedObj->Grabbed = false;
      I->MotionGrabbedObj = NULL;
    }
  }
  return 1;
}

int SceneDeferredClick(DeferredMouse* dm)
{
  SceneClick(dm->block, dm->button, dm->x, dm->y, dm->mod, dm->when);
  return 1;
}

int SceneDeferredRelease(DeferredMouse* dm)
{
  SceneRelease(dm->block, dm->button, dm->x, dm->y, dm->mod, dm->when);
  return 1;
}

/*========================================================================*/
static int SceneDrag(Block* block, int x, int y, int mod, double when)
{
  PyMOLGlobals* G = block->m_G;
  CScene* I = G->Scene;
  float scale, vScale;
  float v1[3], v2[3], n1[3], n2[3], r1, r2, cp[3], v3[3];
  float dx, dy, dt;
  float axis[3], axis2[3], theta, omega;
  float old_front, old_back, old_origin;
  int mode;
  int eff_width;
  int moved_flag;
  int adjust_flag;
  int drag_handled = false;
  int virtual_trackball;
  pymol::CObject* obj;

  if (I->PossibleSingleClick) {
    double slowest_single_click_drag = 0.15;
    if ((when - I->LastClickTime) > slowest_single_click_drag) {
      I->PossibleSingleClick = 0;
    }
  }

  if (I->LoopFlag) {
    return SceneLoopDrag(block, x, y, mod);
  }
  if (I->ButtonsShown && I->PressMode) {
    if (I->ButtonsValid) {
      drag_handled = true;
      I->Over = -1;
      auto elem = &I->SceneVec.front();
      for (std::size_t i = 0u; i < I->SceneVec.size(); i++) {
        auto& elem_ = I->SceneVec[i];
        if (elem_.drawn && elem_.rect.contains(x, y)) {
          elem = &elem_;
          I->Over = i;
          OrthoDirty(G);
          break;
        }
      }
      switch (I->PressMode) {
      case 2:
        if (I->Over >= 0) {
          if (I->Pressed != I->Over) {
            const char* cur_name =
                SettingGetGlobal_s(G, cSetting_scene_current_name);
            if (cur_name && elem->name != cur_name) {
              int animate = -1;
              if (mod & cOrthoCTRL)
                animate = 0;
              auto buffer = pymol::string_format("cmd.scene('''%s''',animate=%d)", elem->name,
                  animate);
              PParse(G, buffer);
              PFlush(G);
              PLog(G, buffer, cPLog_pym);
            }
            I->Pressed = I->Over;
          }
        } else {
          I->Pressed = -1;
        }
      case 3:
        if ((I->Over >= 0) && (I->Pressed != I->Over))
          I->PressMode = 4; /* activate dragging */
        break;
      }

      if (I->PressMode == 4) { /* dragging */
        if ((I->Over >= 0) && (I->Pressed != I->Over) && (I->Pressed >= 0)) {

          auto pressed = &I->SceneVec[I->Pressed];
          std::string buffer;

          if (I->Over > 0) { /* not over the first scene in list */
            SceneElem* first = elem - 1;
            SceneElem* second = pressed;
            if (first >= pressed) {
              first = elem;
              second = pressed;
            }
            buffer = pymol::string_format("cmd.scene_order(['''%s''','''%s'''])", first->name,
                second->name);
          } else {
            buffer = pymol::string_format("cmd.scene_order(['''%s'''],location='top')",
                pressed->name);
          }
          PParse(G, buffer);
          PFlush(G);
          PLog(G, buffer, cPLog_pym);
          I->Pressed = I->Over;
          I->ButtonsValid = false;
          if (SettingGetGlobal_b(G, cSetting_scene_buttons)) {
            OrthoInvalidateDoDraw(G);
          }
        }
      }
    }
  }

  if (!drag_handled) {

    mode = ButModeTranslate(G, I->Button, mod);

    y = y - I->margin.bottom;

    scale = (float) I->Height;
    if (scale > I->Width)
      scale = (float) I->Width;
    scale = 0.45F * scale;
    SceneInvalidateCopy(G, false);
    SceneDontCopyNext(G);
    switch (mode) {
    case cButModePickAtom:
      obj = I->LastPicked.context.object;
      if (obj)
        switch (obj->type) {
        case cObjectGadget: {
          auto gad = static_cast<ObjectGadget*>(obj);

          ObjectGadgetGetVertex(
              gad, I->LastPicked.src.index, I->LastPicked.src.bond, v1);

          vScale = SceneGetExactScreenVertexScale(G, v1);
          if (stereo_via_adjacent_array(I->StereoMode)) {
            x = get_stereo_x(x, &I->LastX, I->Width, NULL);
          }

          /* transform into model coodinate space */
          switch (obj->getRenderContext()) {
          case pymol::RenderContext::UnitWindow: {
            float divisor;
            divisor = (float) I->Width;
            if (I->Height < I->Width)
              divisor = (float) I->Height;
            v2[0] = (x - I->LastX) / divisor;
            v2[1] = (y - I->LastY) / divisor;
            v2[2] = 0;
          } break;
          default:
          case pymol::RenderContext::Camera:
            v2[0] = (x - I->LastX) * vScale;
            v2[1] = (y - I->LastY) * vScale;
            v2[2] = 0;
            MatrixInvTransformC44fAs33f3f(I->m_view.m_rotMatrix, v2, v2);
            break;
          }
          add3f(v1, v2, v2);
          ObjectGadgetSetVertex(
              gad, I->LastPicked.src.index, I->LastPicked.src.bond, v2);
          if (I->LastPicked.src.index) { // pick id on gadget is 1 for band (to
                                         // change values), 0 for everything
                                         // else (to move it around)
            SceneChanged(G); // changing values, need to update gadget text
          } else {
            PyMOL_NeedRedisplay(G->PyMOL); // moving gadget, just re-draw
          }
        } break;
        }
      I->LastX = x;
      I->LastY = y;
      break;
    case cButModeRotDrag:
      eff_width = I->Width;
      if (stereo_via_adjacent_array(I->StereoMode)) {
        eff_width = I->Width / 2;
        x = get_stereo_x(x, &I->LastX, I->Width, NULL);
      }

      virtual_trackball =
          SettingGet_i(G, NULL, NULL, cSetting_virtual_trackball);

      if (virtual_trackball == 2 &&
          (!I->prev_no_z_rotation1 || !I->prev_no_z_rotation2)) {
        /* when virtual_trackball=2 and twisting, need to set v1,v2 relative to
         * orig */
        v2[0] = (float) (eff_width / 2) - I->orig_x_rotation;
        v2[1] = (float) (I->Height / 2) - I->orig_y_rotation;

        v1[0] = v2[0] - (float) (x - I->LastX);
        v1[1] = v2[1] - (float) (y - I->LastY);

      } else if (virtual_trackball == 1) {
        v1[0] = (float) (eff_width / 2) - x;
        v1[1] = (float) (I->Height / 2) - y;

        v2[0] = (float) (eff_width / 2) - I->LastX;
        v2[1] = (float) (I->Height / 2) - I->LastY;

      } else {
        v1[0] = (float) (I->LastX) - x;
        v1[1] = (float) (I->LastY) - y;

        v2[0] = 0;
        v2[1] = 0;
      }

      r1 = (float) sqrt1f(v1[0] * v1[0] + v1[1] * v1[1]);
      r2 = (float) sqrt1f(v2[0] * v2[0] + v2[1] * v2[1]);

      {
        short r1lt, r2lt;
        if (virtual_trackball == 2) {
          r1lt = I->prev_no_z_rotation1;
          r2lt = I->prev_no_z_rotation2;
        } else {
          r1lt = r1 < scale;
          r2lt = r2 < scale;
          I->prev_no_z_rotation1 = r1 < scale;
          I->prev_no_z_rotation2 = r2 < scale;
          I->orig_x_rotation = x;
          I->orig_y_rotation = y;
        }
        if (r1lt) {
          v1[2] = (float) sqrt1f(scale * scale - r1 * r1);
        } else {
          v1[2] = 0.0;
        }
        if (r2lt) {
          v2[2] = (float) sqrt1f(scale * scale - r2 * r2);
        } else {
          v2[2] = 0.0;
        }
      }
      normalize23f(v1, n1);
      normalize23f(v2, n2);
      cross_product3f(n1, n2, cp);
      theta =
          (float) (SettingGet_f(G, NULL, NULL, cSetting_mouse_scale) * 2 * 180 *
                   asin(sqrt1f(cp[0] * cp[0] + cp[1] * cp[1] + cp[2] * cp[2])) /
                   cPI);
      dx = (v1[0] - v2[0]);
      dy = (v1[1] - v2[1]);
      dt = (float) (SettingGet_f(G, NULL, NULL, cSetting_mouse_limit) *
                    sqrt1f(dx * dx + dy * dy) / scale);

      if (theta > dt)
        theta = dt;

      normalize23f(cp, axis);

      axis[2] = -axis[2];

      theta = theta / (1.0F + (float) fabs(axis[2]));

      /* transform into model coodinate space */
      MatrixInvTransformC44fAs33f3f(I->m_view.m_rotMatrix, axis, v2);
      v1[0] = (float) (cPI * theta / 180.0);
      EditorDrag(G, NULL, -1, mode, SettingGetGlobal_i(G, cSetting_state) - 1,
          v1, v2, v3);
      I->LastX = x;
      I->LastY = y;
      break;
    case cButModeMovDrag:
    case cButModeMovDragZ:
      if (I->Threshold) {
        if ((abs(x - I->ThresholdX) > I->Threshold) ||
            (abs(y - I->ThresholdY) > I->Threshold)) {
          I->Threshold = 0;
        }
      }
      if (!I->Threshold) {

        copy3f(I->m_view.m_origin, v1);
        vScale = SceneGetExactScreenVertexScale(G, v1);
        if (stereo_via_adjacent_array(I->StereoMode)) {
          x = get_stereo_x(x, &I->LastX, I->Width, NULL);
        }

        if (mode == cButModeMovDragZ) {
          v2[0] = 0;
          v2[1] = 0;
          v2[2] = -(y - I->LastY) * vScale;
        } else {
          v2[0] = (x - I->LastX) * vScale;
          v2[1] = (y - I->LastY) * vScale;
          v2[2] = 0;
        }

        v3[0] = 0.0F;
        v3[1] = 0.0F;
        v3[2] = 1.0F;

        /* transform into model coodinate space */
        MatrixInvTransformC44fAs33f3f(I->m_view.m_rotMatrix, v2, v2);
        MatrixInvTransformC44fAs33f3f(I->m_view.m_rotMatrix, v3, v3);

        EditorDrag(G, NULL, -1, mode, SettingGetGlobal_i(G, cSetting_state) - 1,
            v1, v2, v3);
      }
      I->LastX = x;
      I->LastY = y;
      break;
    case cButModeRotObj:
    case cButModeMovObj:
    case cButModeMovObjZ:
    case cButModeRotView:
    case cButModeMovView:
    case cButModeMovViewZ:
    case cButModeRotFrag:
    case cButModeMovFrag:
    case cButModeMovFragZ:
    case cButModeTorFrag:
    case cButModeMoveAtom:
    case cButModeMoveAtomZ:
    case cButModePkTorBnd:
      obj = I->LastPicked.context.object;
      if (obj) {
        if (I->Threshold) {
          if ((abs(x - I->ThresholdX) > I->Threshold) ||
              (abs(y - I->ThresholdY) > I->Threshold)) {
            I->Threshold = 0;
          }
        }
        if (!I->Threshold)
          switch (obj->type) {
          case cObjectGadget: /* note repeated above */
          {
            auto gad = static_cast<ObjectGadget*>(obj);

            ObjectGadgetGetVertex(
                gad, I->LastPicked.src.index, I->LastPicked.src.bond, v1);

            vScale = SceneGetExactScreenVertexScale(G, v1);
            if (stereo_via_adjacent_array(I->StereoMode)) {
              x = get_stereo_x(x, &I->LastX, I->Width, NULL);
            }

            /* transform into model coodinate space */
            switch (obj->getRenderContext()) {
            case pymol::RenderContext::UnitWindow: {
              float divisor;
              divisor = (float) I->Width;
              if (I->Height < I->Width)
                divisor = (float) I->Height;
              v2[0] = (x - I->LastX) / divisor;
              v2[1] = (y - I->LastY) / divisor;
              v2[2] = 0;
            } break;
            default:
            case pymol::RenderContext::Camera:
              v2[0] = (x - I->LastX) * vScale;
              v2[1] = (y - I->LastY) * vScale;
              v2[2] = 0;
              MatrixInvTransformC44fAs33f3f(I->m_view.m_rotMatrix, v2, v2);
              break;
            }
            add3f(v1, v2, v2);
            ObjectGadgetSetVertex(
                gad, I->LastPicked.src.index, I->LastPicked.src.bond, v2);
            if (I->LastPicked.src
                    .index) {  // pick id on gadget is 1 for band (to change
                               // values), 0 for everything else (to move it
                               // around)
              SceneChanged(G); // changing values, need to update gadget text
            } else {
              PyMOL_NeedRedisplay(G->PyMOL); // moving gadget, just re-draw
            }
          } break;
          case cObjectMolecule:
            if (ObjectMoleculeGetAtomTxfVertex((ObjectMolecule*) obj,
                    I->LastPicked.context.state, I->LastPicked.src.index, v1)) {
              /* scale properly given the current projection matrix */
              vScale = SceneGetExactScreenVertexScale(G, v1);

              if (stereo_via_adjacent_array(I->StereoMode)) {
                x = get_stereo_x(x, &I->LastX, I->Width, NULL);
              }

              switch (mode) {
              case cButModeMovFragZ:
              case cButModeMovObjZ:
              case cButModeMovViewZ:
              case cButModeMoveAtomZ:
                v2[0] = 0;
                v2[1] = 0;
                v2[2] = -(y - I->LastY) * vScale;
                break;
              default:
                v2[0] = (x - I->LastX) * vScale;
                v2[1] = (y - I->LastY) * vScale;
                v2[2] = 0;
                break;
              }

              v3[0] = 0.0F;
              v3[1] = 0.0F;
              v3[2] = 1.0F;

              /* transform into model coodinate space */
              MatrixInvTransformC44fAs33f3f(I->m_view.m_rotMatrix, v2, v2);
              MatrixInvTransformC44fAs33f3f(I->m_view.m_rotMatrix, v3, v3);

              if (I->LastPicked.src.bond >= cPickableAtom) {
                if ((mode != cButModeMoveAtom) && (mode != cButModeMoveAtomZ)) {
                  EditorDrag(G, obj, I->LastPicked.src.index, mode,
                      SettingGetGlobal_i(G, cSetting_state) - 1, v1, v2, v3);
                  switch (mode) {
                  case cButModeMovViewZ:
                  case cButModeRotView:
                  case cButModeMovView:
                    if (SettingGetGlobal_i(G, cSetting_movie_auto_store) &&
                        SettingGetGlobal_i(
                            G, cSetting_movie_auto_interpolate)) {
                      I->ReinterpolateFlag = true;
                      I->ReinterpolateObj = obj;
                    }
                    break;
                  }
                } else {
                  int log_trans =
                      SettingGetGlobal_b(G, cSetting_log_conformations);
                  ObjectMolecule* objOM = (ObjectMolecule*) obj;
                  ObjectMoleculeMoveAtom(objOM, I->LastPicked.context.state,
                      I->LastPicked.src.index, v2, 1, log_trans);
                  /* -- JV - if this object knows about distances, then move
                   * them if necessary */
                  /* check the dynamic_measures setting and make sure the object
                   * has a distance measure, first  */
                  /* obviated by new method
                  if (SettingGetGlobal_i(G, cSetting_dynamic_measures))
                    ObjectMoleculeMoveDist( (ObjectMolecule *) obj,
                  SettingGetGlobal_i(G, cSetting_state)-1,
                  I->LastPicked.src.index, v2, 1, log_trans);
                  */
                  SceneInvalidate(G);
                }
              } else {
                int log_trans =
                    SettingGetGlobal_b(G, cSetting_log_conformations);
                ObjectMolecule* objOM = (ObjectMolecule*) obj;
                float diffInPixels[2];
                diffInPixels[0] = (x - I->LastX);
                diffInPixels[1] = (y - I->LastY);
                ObjectMoleculeMoveAtomLabel(objOM, I->LastPicked.context.state,
                    I->LastPicked.src.index, v2, log_trans, diffInPixels);
                SceneInvalidate(G);
              }
            }
            break;
          case cObjectSlice: {
            ObjectSlice* slice = (ObjectSlice*) obj;

            if (I->LastPickVertexFlag) {
              auto state = SceneGetState(G);

              copy3f(I->LastPickVertex, v1);

              vScale = SceneGetExactScreenVertexScale(G, v1);

              if (stereo_via_adjacent_array(I->StereoMode)) {
                x = get_stereo_x(x, &I->LastX, I->Width, NULL);
              }

              v2[0] = (x - I->LastX) * vScale;
              v2[1] = (y - I->LastY) * vScale;
              v2[2] = 0;

              v3[0] = 0.0F;
              v3[1] = 0.0F;
              v3[2] = 1.0F;

              /* transform into model coodinate space */
              MatrixInvTransformC44fAs33f3f(I->m_view.m_rotMatrix, v2, v2);
              MatrixInvTransformC44fAs33f3f(I->m_view.m_rotMatrix, v3, v3);

              ObjectSliceDrag(slice, state, mode, v1, v2, v3);
            }
          } break;
          case cObjectMeasurement:
            if (ObjectDistGetLabelTxfVertex((ObjectDist*) obj,
                    I->LastPicked.context.state, I->LastPicked.src.index, v1)) {
              /* scale properly given the current projection matrix */
              vScale = SceneGetExactScreenVertexScale(G, v1);
              if (stereo_via_adjacent_array(I->StereoMode)) {
                x = get_stereo_x(x, &I->LastX, I->Width, NULL);
              }

              switch (mode) {
              case cButModeMovFragZ:
              case cButModeMovObjZ:
              case cButModeMovViewZ:
              case cButModeMoveAtomZ:
                v2[0] = 0;
                v2[1] = 0;
                v2[2] = -(y - I->LastY) * vScale;
                break;
              default:
                v2[0] = (x - I->LastX) * vScale;
                v2[1] = (y - I->LastY) * vScale;
                v2[2] = 0;
                break;
              }

              v3[0] = 0.0F;
              v3[1] = 0.0F;
              v3[2] = 1.0F;

              /* transform into model coodinate space */
              MatrixInvTransformC44fAs33f3f(I->m_view.m_rotMatrix, v2, v2);
              MatrixInvTransformC44fAs33f3f(I->m_view.m_rotMatrix, v3, v3);

              if (I->LastPicked.src.bond == cPickableLabel) {
                int log_trans =
                    SettingGetGlobal_b(G, cSetting_log_conformations);
                ObjectDistMoveLabel((ObjectDist*) obj,
                    I->LastPicked.context.state, I->LastPicked.src.index, v2, 1,
                    log_trans);
                SceneInvalidate(G);
              }
            }
            break;
          default:
            break;
          }
      }
      I->LastX = x;
      I->LastY = y;
      break;
    case cButModeTransXY:

      SceneNoteMouseInteraction(G);

      vScale = SceneGetExactScreenVertexScale(G, I->m_view.m_origin);
      if (stereo_via_adjacent_array(I->StereoMode)) {

        x = get_stereo_x(x, &I->LastX, I->Width, NULL);
      }

      v2[0] = (x - I->LastX) * vScale;
      v2[1] = (y - I->LastY) * vScale;
      v2[2] = 0.0F;

      moved_flag = false;
      if (I->LastX != x) {
        I->m_view.m_pos[0] += v2[0];
        I->LastX = x;
        SceneInvalidate(G);
        moved_flag = true;
      }
      if (I->LastY != y) {
        I->m_view.m_pos[1] += v2[1];
        I->LastY = y;
        SceneInvalidate(G);
        moved_flag = true;
      }

      EditorFavorOrigin(G, NULL);
      if (moved_flag && SettingGetGlobal_b(G, cSetting_roving_origin)) {
        SceneGetCenter(G, v2); /* gets position of center of screen */
        SceneOriginSet(G, v2, true);
      }
      if (moved_flag && SettingGetGlobal_b(G, cSetting_roving_detail)) {
        SceneRovingDirty(G);
      }
      break;
    case cButModeRotXYZ:
    case cButModeRotZ:
    case cButModeInvRotZ:
    case cButModeTransZ:
    case cButModeClipNF:
    case cButModeClipN:
    case cButModeClipF:
    case cButModeRotL:
    case cButModeMovL:
    case cButModeMvzL:

      SceneNoteMouseInteraction(G);
      virtual_trackball =
          SettingGet_i(G, NULL, NULL, cSetting_virtual_trackball);
      eff_width = I->Width;
      if (stereo_via_adjacent_array(I->StereoMode)) {
        eff_width = I->Width / 2;
        x = get_stereo_x(x, &I->LastX, I->Width, NULL);
      }
      if (virtual_trackball == 2 &&
          (!I->prev_no_z_rotation1 || !I->prev_no_z_rotation2)) {
        /* when virtual_trackball=2 and twisting, need to set v1,v2 relative to
         * orig */
        v2[0] = (float) (eff_width / 2) - I->orig_x_rotation;
        v2[1] = (float) (I->Height / 2) - I->orig_y_rotation;

        v1[0] = v2[0] - (float) (x - I->LastX);
        v1[1] = v2[1] - (float) (y - I->LastY);

      } else if (virtual_trackball == 1) {
        v1[0] = (float) (eff_width / 2) - x;
        v1[1] = (float) (I->Height / 2) - y;

        v2[0] = (float) (eff_width / 2) - I->LastX;
        v2[1] = (float) (I->Height / 2) - I->LastY;

      } else {
        v1[0] = (float) (I->LastX) - x;
        v1[1] = (float) (I->LastY) - y;

        v2[0] = 0;
        v2[1] = 0;
      }

      r1 = (float) sqrt1f(v1[0] * v1[0] + v1[1] * v1[1]);
      r2 = (float) sqrt1f(v2[0] * v2[0] + v2[1] * v2[1]);

      {
        short r1lt, r2lt;
        if (SettingGet_i(G, NULL, NULL, cSetting_virtual_trackball) == 2) {
          r1lt = I->prev_no_z_rotation1;
          r2lt = I->prev_no_z_rotation2;
        } else {
          r1lt = r1 < scale;
          r2lt = r2 < scale;
          I->prev_no_z_rotation1 = r1 < scale;
          I->prev_no_z_rotation2 = r2 < scale;
          I->orig_x_rotation = x;
          I->orig_y_rotation = y;
        }
        if (r1lt) {
          float val = scale * scale - r1 * r1;
          short isNeg = val < 0.f;
          v1[2] = (float) sqrt(fabs(val)) * (isNeg ? -1.f : 1.f);
        } else {
          v1[2] = 0.0;
        }
        if (r2lt) {
          float val = scale * scale - r2 * r2;
          short isNeg = val < 0.f;
          v2[2] = (float) sqrt(fabs(val)) * (isNeg ? -1.f : 1.f);
        } else {
          v2[2] = 0.0;
        }
      }
      normalize23f(v1, n1);
      normalize23f(v2, n2);
      cross_product3f(n1, n2, cp);
      theta =
          (float) (SettingGet_f(G, NULL, NULL, cSetting_mouse_scale) * 2 * 180 *
                   asin(sqrt1f(cp[0] * cp[0] + cp[1] * cp[1] + cp[2] * cp[2])) /
                   cPI);

      dx = (v1[0] - v2[0]);
      dy = (v1[1] - v2[1]);
      dt = (float) (SettingGet_f(G, NULL, NULL, cSetting_mouse_limit) *
                    sqrt1f(dx * dx + dy * dy) / scale);
      if (theta > dt)
        theta = dt;

      normalize23f(cp, axis);

      theta = theta / (1.0F + (float) fabs(axis[2]));

      v1[2] = 0.0;
      v2[2] = 0.0;
      normalize23f(v1, n1);
      normalize23f(v2, n2);
      cross_product3f(n1, n2, cp);
      omega =
          (float) (2 * 180 *
                   asin(sqrt1f(cp[0] * cp[0] + cp[1] * cp[1] + cp[2] * cp[2])) /
                   cPI);
      normalize23f(cp, axis2);

      old_front = I->m_view.m_clip.m_front;
      old_back = I->m_view.m_clip.m_back;
      old_origin = -I->m_view.m_pos[2];

      moved_flag = false;
      adjust_flag = false;
      switch (mode) {
      case cButModeRotXYZ:
        if (I->LastX != x) {
          SceneRotate(G, theta, axis[0], axis[1], -axis[2]);
          I->LastX = x;
          adjust_flag = true;
        }
        if (I->LastY != y) {
          SceneRotate(G, theta, axis[0], axis[1], -axis[2]);
          I->LastY = y;
          adjust_flag = true;
        }
        break;
      case cButModeRotZ:
        if (I->LastX != x) {
          SceneRotate(G, omega, axis2[0], axis2[1], -axis2[2]);
          I->LastX = x;
          adjust_flag = true;
        }
        if (I->LastY != y) {
          SceneRotate(G, omega, axis2[0], axis2[1], -axis2[2]);
          I->LastY = y;
          adjust_flag = true;
        }
        break;
      case cButModeInvRotZ:
        if (I->LastX != x) {
          SceneRotate(G, (I->LastX - x) / 2.0F, 0.0F, 0.0F, 1.0F);
          I->LastX = x;
          adjust_flag = true;
        }
        break;
      case cButModeTransZ:
        // zoom
        if (I->LastY != y) {
          {
            // the 5.0 min depth below is an empirical value
            float factor = SettingGetGlobal_f(G, cSetting_mouse_z_scale) *
                           (y - I->LastY) / 400.0 *
                           std::max(5.f, -I->m_view.m_pos[2]);
            if (!SettingGetGlobal_b(G, cSetting_legacy_mouse_zoom))
              factor = -factor;
            I->m_view.m_pos[2] += factor;
            I->m_view.m_clip.m_front -= factor;
            I->m_view.m_clip.m_back -= factor;
            UpdateFrontBackSafe(I);
          }
          I->LastY = y;
          SceneInvalidate(G);
          adjust_flag = true;
        }
        break;
      case cButModeClipNF:
        if (I->LastX != x) {
          I->m_view.m_clip.m_back -= (((float) x) - I->LastX) / 10;
          I->LastX = x;
          moved_flag = true;
        }
        if (I->LastY != y) {
          I->m_view.m_clip.m_front -= (((float) y) - I->LastY) / 10;
          I->LastY = y;
          moved_flag = true;
        }
        if (moved_flag) {
          SceneClipSet(G, I->m_view.m_clip.m_front, I->m_view.m_clip.m_back);
        }
        break;
      case cButModeClipN:
        if (I->LastX != x || I->LastY != y) {
          I->m_view.m_clip.m_front -= (x - I->LastX + y - I->LastY) / 10.f;
          I->LastX = x;
          I->LastY = y;
          SceneClipSet(G, I->m_view.m_clip.m_front, I->m_view.m_clip.m_back);
          moved_flag = true;
        }
        break;
      case cButModeClipF:
        if (I->LastX != x || I->LastY != y) {
          I->m_view.m_clip.m_back -= (x - I->LastX + y - I->LastY) / 10.f;
          I->LastX = x;
          I->LastY = y;
          SceneClipSet(G, I->m_view.m_clip.m_front, I->m_view.m_clip.m_back);
          moved_flag = true;
        }
        break;
      case cButModeRotL:
      case cButModeMovL: {
        /* when light_count == 1, there is an ambient light;
         * when light_count == 2, there are two lights, the ambient and
         * a directional, called "light". When there are three, it's the
         * ambient, light and the third is light2, and so on.
         *
         * Should we use an off-by-one here to make this easier for the
         * user to understand? light1 is ambient, light2 is first directional
         *
         */
        float pos[3];
        int which_light;
        float ms = 0.01;

        which_light = light_setting_indices
            [glm::clamp(SettingGetGlobal_i(G, cSetting_edit_light), 1, 9) - 1];

        copy3f(SettingGet<const float*>(G, which_light), pos);

        pos[0] += (float) dx * ms;
        pos[1] += (float) dy * ms;

        SettingSet_3fv(G->Setting, which_light, pos);
        SettingGenerateSideEffects(G, which_light, NULL, 0, 1);
        I->LastX = x;
        I->LastY = y;
      } break;
      case cButModeMvzL: {
        float pos[3];
        int which_light;
        float ms = 0.01;
        float factor = 0.f;

        /* when light_count == 1, there is an ambient light;
         * when light_count == 2, there are two lights, the ambient and
         * a directional, called "light". When there are three, it's the
         * ambient, light and the third is light2, and so on.
         */

        which_light = light_setting_indices
            [glm::clamp(SettingGetGlobal_i(G, cSetting_edit_light), 1, 9) - 1];

        if (I->LastY != y) {
          factor =
              400 /
              ((I->m_view.m_clipSafe.m_front + I->m_view.m_clipSafe.m_back) /
                  2);
          if (factor >= 0.0F) {
            factor = SettingGetGlobal_f(G, cSetting_mouse_z_scale) *
                     (((float) y) - I->LastY) / factor;
            if (!SettingGetGlobal_b(G, cSetting_legacy_mouse_zoom))
              factor = -factor;
          }
        }

        copy3f(SettingGet<const float*>(G, which_light), pos);
        SettingGenerateSideEffects(G, which_light, NULL, 0, 1);

        pos[2] -= factor * ms;

        SettingSet_3fv(G->Setting, which_light, pos);
        I->LastX = x;
        I->LastY = y;
      } break;
      }
      if (moved_flag)
        SceneDoRoving(G, old_front, old_back, old_origin, adjust_flag, false);
    }
  }
  if (I->PossibleSingleClick) {
    int max_single_click_drag = 4;
    int dx = abs(I->StartX - I->LastX);
    int dy = abs(I->StartY - I->LastY);
    if ((dx > max_single_click_drag) || (dy > max_single_click_drag)) {
      I->PossibleSingleClick = false;
    }
  }
  return (1);
}

int SceneDeferredDrag(DeferredMouse* dm)
{
  SceneDrag(dm->block, dm->x, dm->y, dm->mod, dm->when);
  return 1;
}
