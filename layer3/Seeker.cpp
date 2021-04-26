
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

#include <algorithm>

#include"os_python.h"
#include"os_predef.h"
#include"os_std.h"
#include"os_gl.h"
#include"Err.h"
#include"Util.h"
#include"Seq.h"
#include"Seeker.h"
#include"MemoryDebug.h"
#include"Executive.h"
#include"P.h"
#include"Selector.h"
#include"Wizard.h"
#include"Scene.h"
#include"Menu.h"
#include "Lex.h"

struct CSeeker : public CSeqHandler {
  bool dragging;
  SeekerDragInfo dragInfo{};
  double LastClickTime;


  CSeqRow* click(PyMOLGlobals * G, std::vector<CSeqRow>& rowVLA, int button, int row, int col,
                      int mod, int x, int y) override;
  CSeqRow* drag(PyMOLGlobals * G, std::vector<CSeqRow>& rowVLA, int row, int col, int mod) override;
  CSeqRow* release(PyMOLGlobals * G, std::vector<CSeqRow>& rowVLA, int button, int row, int col,
                        int mod) override;
  void refresh(PyMOLGlobals * G, std::vector<CSeqRow>& rowVLA) override;
};

static CSeqRow *SeekerDrag(PyMOLGlobals * G, std::vector<CSeqRow>& rowVLA, int row, int col, int mod);

static void SeekerBuildSeleFromAtomList(PyMOLGlobals * G, const char *obj_name, int *atom_list,
                                        const char *sele_name, int start_fresh)
{
  ObjectMolecule *obj = ExecutiveFindObjectMoleculeByName(G, obj_name);

  if(start_fresh) {
    SelectorCreateFromObjectIndices(G, sele_name, obj, atom_list, -1);
  } else {
    SelectorCreateFromObjectIndices(G, cTempSeekerSele2, obj, atom_list, -1);

    auto buf1 = pymol::string_format("?%s|?%s", sele_name, cTempSeekerSele2);
    SelectorCreate(G, sele_name, buf1.c_str(), nullptr, true, nullptr);
    ExecutiveDelete(G, cTempSeekerSele2);
  }
}

static void SeekerSelectionToggleRange(PyMOLGlobals * G, std::vector<CSeqRow>& rowVLA, int row_num,
                                       int col_first, int col_last, int inc_or_excl,
                                       int start_over)
{
  char selName[WordLength];

  if(row_num >= 0) {
    CSeqCol *col;
    char prefix[3] = "";
    int logging = SettingGetGlobal_i(G, cSetting_logging);
    int col_num;
    int *atom_vla = NULL;
    int n_at = 0;
    int at_idx;
    int *atom_list;

    ObjectMolecule *obj;
    if(logging == cPLog_pml)
      strcpy(prefix, "_ ");
    auto row = &rowVLA[row_num];
    if((obj = ExecutiveFindObjectMoleculeByName(G, row->name))) {
      atom_vla = VLAlloc(int, obj->NAtom / 10);
      for(col_num = col_first; col_num <= col_last; col_num++) {
        col = row->col + col_num;
        if(!col->spacer) {
          if(!start_over) {
            if(inc_or_excl)
              col->inverse = true;
            else
              col->inverse = false;
          } else {
            col->inverse = true;
          }
          atom_list = row->atom_lists + col->atom_at;
          while((at_idx = (*(atom_list++))) >= 0) {     /* build one extra long list 
                                                           so that we only call selector once */
            VLACheck(atom_vla, int, n_at);
            atom_vla[n_at++] = at_idx;
          }
        }
      }
      VLACheck(atom_vla, int, n_at);
      atom_vla[n_at] = -1;
      SeekerBuildSeleFromAtomList(G, row->name, atom_vla, cTempSeekerSele, true);
      VLAFreeP(atom_vla);

      {
        const char *sele_mode_kw;
        sele_mode_kw = SceneGetSeleModeKeyword(G);

        if(logging)
          SelectorLogSele(G, cTempSeekerSele);

        {

          std::string buf1;
          ExecutiveGetActiveSeleName(G, selName, true, logging);

          /* selection or deselecting? */

          if(!start_over) {
            if(inc_or_excl) {
              buf1 = pymol::string_format("((%s(?%s)) or %s(?%s))",
                      sele_mode_kw, selName, sele_mode_kw, cTempSeekerSele);
            } else {
              buf1 = pymol::string_format("((%s(?%s)) and not %s(?%s))",
                      sele_mode_kw, selName, sele_mode_kw, cTempSeekerSele);
            }
          } else {
            buf1 = pymol::string_format("%s(?%s)", sele_mode_kw, cTempSeekerSele);
          }

          /* create the new active selection */

          SelectorCreate(G, selName, buf1.c_str(), nullptr, true, nullptr);
          {
            auto buf2 = pymol::string_format("%scmd.select(\"%s\",\"%s\",enable=1)\n", prefix, selName,
                    buf1);
            PLog(G, buf2, cPLog_no_flush);
          }
          WizardDoSelect(G, selName);

        }

        ExecutiveDelete(G, cTempSeekerSele);
        if(logging) {
          auto buf2 = pymol::string_format("%scmd.delete(\"%s\")\n", prefix, cTempSeekerSele);
          PLog(G, buf2, cPLog_no_flush);
          PLogFlush(G);
        }

        if(SettingGetGlobal_b(G, cSetting_auto_show_selections))
          ExecutiveSetObjVisib(G, selName, 1, false);
        SceneInvalidate(G);
      }
    }
  }
}

static void SeekerSelectionToggle(PyMOLGlobals * G, std::vector<CSeqRow>& rowVLA, int row_num,
                                  int col_num, int inc_or_excl, int start_over)
{
  char selName[WordLength];

  if(row_num >= 0) {
    int *atom_list;
    char prefix[3] = "";
    int logging = SettingGetGlobal_i(G, cSetting_logging);

    if(logging == cPLog_pml)
      strcpy(prefix, "_ ");
    auto row = &rowVLA[row_num];
    auto col = &row->col[col_num];
    if(!col->spacer)
      if(ExecutiveFindObjectByName(G, row->name)) {
        const char *sele_mode_kw;
        atom_list = row->atom_lists + col->atom_at;

        /* build up a selection consisting of residue atoms */

        SeekerBuildSeleFromAtomList(G, row->name, atom_list, cTempSeekerSele, true);
        sele_mode_kw = SceneGetSeleModeKeyword(G);

        if(logging)
          SelectorLogSele(G, cTempSeekerSele);

        {
          std::string buf1;
          ExecutiveGetActiveSeleName(G, selName, true, logging);

          /* selection or deselecting? */

          if(!start_over) {
            if(inc_or_excl) {
              if(!col->spacer) {
                col->inverse = true;
                buf1 = pymol::string_format("((%s(?%s)) or %s(%s))",
                        sele_mode_kw, selName, sele_mode_kw, cTempSeekerSele);
              }
            } else {
              if(!col->spacer) {
                col->inverse = false;
                buf1 = pymol::string_format("((%s(?%s)) and not %s(%s))",
                        sele_mode_kw, selName, sele_mode_kw, cTempSeekerSele);
              }
            }
          } else {
            if(!col->spacer) {
              col->inverse = true;
              buf1 = pymol::string_format("%s(%s)", sele_mode_kw, cTempSeekerSele);
            }
          }

          /* create the new active selection */

          SelectorCreate(G, selName, buf1.c_str(), nullptr, true, nullptr);
          {
            auto buf2 = pymol::string_format("%scmd.select(\"%s\",\"%s\",enable=1)\n", prefix, selName,
                    buf1);
            PLog(G, buf2, cPLog_no_flush);
          }
          WizardDoSelect(G, selName);
        }

        ExecutiveDelete(G, cTempSeekerSele);
        if(logging) {
          auto buf2 = pymol::string_format("%scmd.delete(\"%s\")\n", prefix, cTempSeekerSele);
          PLog(G, buf2, cPLog_no_flush);
          PLogFlush(G);
        }

        if(SettingGetGlobal_b(G, cSetting_auto_show_selections))
          ExecutiveSetObjVisib(G, selName, 1, false);
        SceneInvalidate(G);
      }
  }
}

void SeekerSelectionUpdateCenter(PyMOLGlobals * G, std::vector<CSeqRow>& rowVLA, int row_num,
                                        int col_num, int start_over)
{
  if (row_num < 0) {
     return;
  }
  auto row = &rowVLA[row_num];
  auto col = &row->col[col_num];

  if(col->spacer) {
    return;
  }
  if(auto obj = ExecutiveFindObjectByName(G, row->name)) {

    if(col->state)
      SettingSetSmart_i(G, obj->Setting.get(), nullptr, cSetting_state, col->state);

    auto* atom_list = &row->atom_lists[col->atom_at];

    SeekerBuildSeleFromAtomList(G, row->name, atom_list, cTempCenterSele,
                                start_over);
    auto logging = SettingGet<bool>(G, cSetting_logging);
    if(logging)
      SelectorLogSele(G, cTempCenterSele);
  }
}

void SeekerSelectionCenter(PyMOLGlobals * G, int action)
{
  char prefix[3] = "";
  int logging = SettingGetGlobal_i(G, cSetting_logging);
  if(logging == cPLog_pml)
    strcpy(prefix, "_ ");

  switch (action) {
  case 0:                      /* center cumulative */
    ExecutiveCenter(G, cTempCenterSele, -1, true, -1, NULL, true);
    if(logging) {
      auto buf2 = pymol::string_format("%scmd.center(\"%s\")\n", prefix, cTempCenterSele);
      PLog(G, buf2, cPLog_no_flush);
      PLogFlush(G);
    }
    break;
  case 1:                      /* zoom */
    ExecutiveWindowZoom(G, cTempCenterSele, 0.0, -1, false, -1, true);
    if(logging) {
      auto buf2 = pymol::string_format("%scmd.zoom(\"%s\")\n", prefix, cTempCenterSele);
      PLog(G, buf2, cPLog_no_flush);
      PLogFlush(G);
    }
    break;
  case 2:                      /* center seeker */
    {
      char selName[WordLength];
      if(ExecutiveGetActiveSeleName(G, selName, true, logging)) {
        ExecutiveCenter(G, selName, -1, true, -1, NULL, true);
        if(logging) {
          auto buf2 = pymol::string_format("%scmd.center(\"%s\")\n", prefix, selName);
          PLog(G, buf2, cPLog_no_flush);
          PLogFlush(G);
        }
      }
    }
    break;
  }
}

#define cDoubleTime 0.35

static CSeqRow *SeekerClick(PyMOLGlobals * G, std::vector<CSeqRow>& rowVLA, int button, int row_num,
                            int col_num, int mod, int x, int y)
{
  CSeqRow *row;
  CSeqCol *col;
  /*  char selName[WordLength]; */
  CSeeker *I = G->Seeker;
  int logging = SettingGetGlobal_i(G, cSetting_logging);
  int continuation = false;
  if((row_num < 0) || (col_num < 0)) {
    switch (button) {
    case P_GLUT_LEFT_BUTTON:
      if((UtilGetSeconds(G) - I->LastClickTime) < cDoubleTime) {
        char name[WordLength];
        if(ExecutiveGetActiveSeleName(G, name, false, false)) {
          SelectorCreate(G, name, "none", NULL, true, NULL);
          if(logging) {
            auto buf2 = pymol::string_format("cmd.select('%s','none', enable=1)", name);
            PLog(G, buf2, cPLog_no_flush);
          }
          SeqDirty(G);
        }
      }
      I->LastClickTime = UtilGetSeconds(G);
      break;
    }
  } else {
    row = rowVLA.data() + row_num;
    col = row->col + col_num;
    I->dragging = false;
    I->dragInfo.button = button;
    I->box_row = row_num;
    I->box_stop_col = col_num;
    if((I->dragInfo.row == row_num) && (button == P_GLUT_LEFT_BUTTON) && (mod & cOrthoSHIFT)) {
      continuation = true;
    } else {
      I->dragInfo.row = -1;         /* invalidate */
      I->box_start_col = col_num;
    }

    switch (button) {
    case P_GLUT_RIGHT_BUTTON:
      {
        ObjectMolecule *obj;
        char name[WordLength];

        if(ExecutiveGetActiveSeleName(G, name, false, logging) && col->inverse) {
          MenuActivate2Arg(G, x, y + 16, x, y, false, "pick_sele", name, name);
        } else if((obj = ExecutiveFindObjectMoleculeByName(G, row->name))) {
          {
            int *atom_list;
            char prefix[3] = "";
            int logging = SettingGetGlobal_i(G, cSetting_logging);

            if(logging == cPLog_pml)
              strcpy(prefix, "_ ");

            if(ExecutiveFindObjectByName(G, row->name)) {
              atom_list = row->atom_lists + col->atom_at;

              /* build up a selection consisting of residue atoms */

              if((*atom_list) >= 0) {

                auto buffer = ObjectMoleculeGetAtomSele(obj, *atom_list);

                SeekerBuildSeleFromAtomList(G, row->name, atom_list, cTempSeekerSele,
                                            true);
                if(logging)
                  SelectorLogSele(G, cTempSeekerSele);

                MenuActivate2Arg(G, x, y + 16, x, y, false, "seq_option", cTempSeekerSele,
                                 buffer.c_str());

              }
            }
          }
        }
      }
      break;
    case P_GLUT_MIDDLE_BUTTON:
      if(!col->spacer) {
        ObjectMolecule *obj;
        I->dragInfo.start_col = col_num;
        I->dragInfo.last_col = col_num;
        I->dragInfo.row = row_num;
        I->dragging = true;
        SeekerSelectionUpdateCenter(G, rowVLA, row_num, col_num, true);
        int action = mod & cOrthoCTRL ? 1 : 0;
        SeekerSelectionCenter(G, action);
        I->box_active = true;
        if(col->state && (obj = ExecutiveFindObjectMoleculeByName(G, row->name))) {
          SettingSetSmart_i(G, obj->Setting.get(), NULL, cSetting_state, col->state);
          SceneChanged(G);
        }
      }
      break;
    case P_GLUT_LEFT_BUTTON:
      if(!col->spacer) {
        int start_over = false;
        int center = 0;
        ObjectMolecule *obj;
        if(mod & cOrthoCTRL) {
          center = 2;
        }
        int codes = SettingGet_i(G, row->obj->Setting.get(), NULL, cSetting_seq_view_format);
        if(row->obj->DiscreteFlag && SettingGet_b(G,
                                           row->obj->Setting.get(),
                                           NULL, cSetting_seq_view_discrete_by_state))
          codes = 4;
        if (codes != 4 || row->obj->DiscreteFlag) { // keep only non-discrete states selectable
          if(!continuation) {
            I->dragInfo.start_col = col_num;
            I->dragInfo.last_col = col_num;
            I->dragInfo.row = row_num;
            I->dragInfo.dir = 0;
            I->dragInfo.start_toggle = true;
          } else {
            int tmp;
            if(((col_num < I->dragInfo.start_col) && (I->dragInfo.last_col > I->dragInfo.start_col)) ||
               ((col_num > I->dragInfo.start_col) && (I->dragInfo.last_col < I->dragInfo.start_col))) {
              tmp = I->dragInfo.last_col;
              I->dragInfo.last_col = I->dragInfo.start_col;
              I->dragInfo.start_col = tmp;
              I->dragInfo.dir = -I->dragInfo.dir;
            }
          }
          I->dragging = true;

          I->box_active = true;
          if(continuation) {
            SeekerDrag(G, rowVLA, row_num, col_num, mod);
          } else {
            if(col->inverse && !start_over) {
              SeekerSelectionToggle(G, rowVLA, row_num, col_num, false, false);
              I->dragInfo.setting = false;
            } else {
              SeekerSelectionToggle(G, rowVLA, row_num, col_num, true, start_over);
              I->dragInfo.setting = true;
            }
          }
        }

        if(center)
          SeekerSelectionCenter(G, 2);

        if(col->state && (obj = ExecutiveFindObjectMoleculeByName(G, row->name))) {
          SettingSetSmart_i(G, obj->Setting.get(), NULL, cSetting_state, col->state);
          SceneChanged(G);
        }
      }
      break;
    }
  }

  return NULL;
}

static void SeekerRefresh(PyMOLGlobals * G, std::vector<CSeqRow>& rowVLA)
{
  if(rowVLA.empty()) {
    return;
  }
  auto nRow = rowVLA.size();
  int sele = ExecutiveGetActiveSele(G);

  if(sele < 0)
    sele = SelectorIndexByName(G, "_seeker_hilight");

  for(int b = 0; b < nRow; b++) {
    auto row = &rowVLA[b];

    ObjectMolecule *obj;
    if((obj = ExecutiveFindObjectMoleculeByName(G, row->name))) {
      const AtomInfoType *atInfo = obj->AtomInfo.data();
      int selected;

      if(sele < 0) {
        for(int a = 0; a < row->nCol; a++) {
          auto col = &row->col[a];
          col->inverse = false;
        }
      } else {
        for(int a = 0; a < row->nCol; a++) {

          auto col = &row->col[a];
          if(!col->spacer) {
            selected = false;
            auto atom_list = &row->atom_lists[col->atom_at];

            int at;
            while((at = (*atom_list)) >= 0) {
              atom_list++;
              if(SelectorIsMember(G, atInfo[at].selEntry, sele)) {
                selected = true;
              }
            }

            if(selected)
              col->inverse = true;
            else
              col->inverse = false;
          } else
            col->inverse = false;
        }
      }
    }
  }
}

static CSeqRow *SeekerDrag(PyMOLGlobals * G, std::vector<CSeqRow>& rowVLA, int row, int col, int mod)
{
  CSeeker *I = G->Seeker;
  int a;

  if((row >= 0) && (col >= 0) && (I->dragging)) {
    I->box_stop_col = col;

    switch (I->dragInfo.button) {
    case P_GLUT_LEFT_BUTTON:
      if(col != I->dragInfo.last_col) {

        if(I->dragInfo.dir) {
          if(I->dragInfo.dir > 0) {
            if(col <= I->dragInfo.start_col) {
              col = I->dragInfo.start_col;
              if(I->dragInfo.start_toggle) {
                SeekerSelectionToggle(G, rowVLA, I->dragInfo.row, I->dragInfo.start_col,
                                      !I->dragInfo.setting, false);
                I->dragInfo.start_toggle = false;
              }
            } else if(col > I->dragInfo.start_col) {
              if(!I->dragInfo.start_toggle) {
                SeekerSelectionToggle(G, rowVLA, I->dragInfo.row, I->dragInfo.start_col,
                                      I->dragInfo.setting, false);
                I->dragInfo.start_toggle = true;
              }
            }
          } else if(I->dragInfo.dir < 0) {
            if(col >= I->dragInfo.start_col) {
              col = I->dragInfo.start_col;
              if(I->dragInfo.start_toggle) {
                SeekerSelectionToggle(G, rowVLA, I->dragInfo.row, I->dragInfo.start_col,
                                      !I->dragInfo.setting, false);
                I->dragInfo.start_toggle = false;
              }
            } else if(col < I->dragInfo.start_col) {
              if(!I->dragInfo.start_toggle) {
                SeekerSelectionToggle(G, rowVLA, I->dragInfo.row, I->dragInfo.start_col,
                                      I->dragInfo.setting, false);
                I->dragInfo.start_toggle = true;
              }
            }
          }
        }

        if((I->dragInfo.last_col < I->dragInfo.start_col) && (col > I->dragInfo.start_col)) {
          SeekerSelectionToggleRange(G, rowVLA, I->dragInfo.row, I->dragInfo.last_col,
                                     I->dragInfo.start_col - 1, !I->dragInfo.setting, false);
          I->dragInfo.last_col = I->dragInfo.start_col;
        }
        if((I->dragInfo.last_col > I->dragInfo.start_col) && (col < I->dragInfo.start_col)) {
          SeekerSelectionToggleRange(G, rowVLA, I->dragInfo.row, I->dragInfo.start_col + 1,
                                     I->dragInfo.last_col, !I->dragInfo.setting, false);
          I->dragInfo.last_col = I->dragInfo.start_col;
        }
        if(I->dragInfo.start_col == I->dragInfo.last_col) {
          if(col > I->dragInfo.start_col) {
            if(!I->dragInfo.dir)
              I->dragInfo.dir = 1;
            I->dragInfo.last_col = I->dragInfo.start_col + 1;
            SeekerSelectionToggle(G, rowVLA, I->dragInfo.row, I->dragInfo.last_col,
                                  I->dragInfo.setting, false);
          } else if(col < I->dragInfo.start_col) {
            if(!I->dragInfo.dir)
              I->dragInfo.dir = -1;
            I->dragInfo.last_col = I->dragInfo.start_col - 1;
            SeekerSelectionToggle(G, rowVLA, I->dragInfo.row, I->dragInfo.last_col,
                                  I->dragInfo.setting, false);
          }
        }
        if(I->dragInfo.start_col < I->dragInfo.last_col) {

          if(col > I->dragInfo.last_col) {
            SeekerSelectionToggleRange(G, rowVLA, I->dragInfo.row, I->dragInfo.last_col + 1, col,
                                       I->dragInfo.setting, false);
          } else {
            SeekerSelectionToggleRange(G, rowVLA, I->dragInfo.row, col + 1, I->dragInfo.last_col,
                                       !I->dragInfo.setting, false);
          }
        } else {

          if(col < I->dragInfo.last_col) {
            SeekerSelectionToggleRange(G, rowVLA, I->dragInfo.row, col, I->dragInfo.last_col - 1,
                                       I->dragInfo.setting, false);
          } else {
            SeekerSelectionToggleRange(G, rowVLA, I->dragInfo.row, I->dragInfo.last_col, col - 1,
                                       !I->dragInfo.setting, false);
          }
        }
        I->dragInfo.last_col = col;

        if(mod & cOrthoCTRL) {
          SeekerSelectionCenter(G, 2);
        }

      }
      break;
    case P_GLUT_MIDDLE_BUTTON:
      if(col != I->dragInfo.last_col) {
        int action = 0;
        int start_over = false;

        if(mod & cOrthoCTRL) {
          action = 1;
        }
        if(!(mod & cOrthoSHIFT)) {
          start_over = true;
          I->box_start_col = col;
          SeekerSelectionUpdateCenter(G, rowVLA, I->dragInfo.row, col, start_over);
        } else {
          if(I->dragInfo.start_col == I->dragInfo.last_col) {
            if(col > I->dragInfo.start_col) {
              I->dragInfo.last_col = I->dragInfo.start_col + 1;
              SeekerSelectionUpdateCenter(G, rowVLA, I->dragInfo.row, I->dragInfo.last_col,
                                          start_over);
            } else if(col < I->dragInfo.start_col) {
              I->dragInfo.last_col = I->dragInfo.start_col - 1;
              SeekerSelectionUpdateCenter(G, rowVLA, I->dragInfo.row, I->dragInfo.last_col,
                                          start_over);
            }
          }
          if(I->dragInfo.start_col < I->dragInfo.last_col) {

            if(col > I->dragInfo.last_col) {
              for(a = I->dragInfo.last_col + 1; a <= col; a++) {
                SeekerSelectionUpdateCenter(G, rowVLA, I->dragInfo.row, a, start_over);
              }
            }
          } else {

            if(col < I->dragInfo.last_col) {
              for(a = I->dragInfo.last_col - 1; a >= col; a--) {
                SeekerSelectionUpdateCenter(G, rowVLA, I->dragInfo.row, a, start_over);
              }
            }
          }
        }
        I->dragInfo.last_col = col;

        SeekerSelectionCenter(G, action);
      }
      break;
    }
  }
  return NULL;
}

static CSeqRow *SeekerRelease(PyMOLGlobals * G, std::vector<CSeqRow>& rowVLA, int button,
                              int row, int col, int mod)
{
  CSeeker *I = G->Seeker;
  I->dragging = false;

  I->box_active = false;
  return NULL;
}

char SeekerGetAbbr(PyMOLGlobals * G, const char *abbr, char water, char unknown)
{

  switch (abbr[0]) {
  case 'A':
    switch (abbr[1]) {
    case 'L':
      if(abbr[2] == 'A')
        return 'A';
      break;
    case 'R':
      if(abbr[2] == 'G')
        return 'R';
      break;
    case 'S':
      switch (abbr[2]) {
      case 'P':
        return 'D';
        break;
      case 'N':
        return 'N';
        break;
      }
      break;
    }
    break;
  case 'C':
    switch (abbr[1]) {
    case 'Y':
      switch (abbr[2]) {
      case 'S':
      case 'X':
        return 'C';
        break;
      }
      break;
    }
    break;
  case 'G':
    switch (abbr[1]) {
    case 'L':
      switch (abbr[2]) {
      case 'N':
        return 'Q';
        break;
      case 'U':
        return 'E';
        break;
      case 'Y':
        return 'G';
        break;
      }
    }
    break;
  case 'H':
    switch (abbr[1]) {
    case 'I':
      switch (abbr[2]) {
      case 'S':
      case 'D':
      case 'E':
        return 'H';
        break;
      }
      break;
    case 'O':
      switch (abbr[2]) {
      case 'H':
        return water;
        break;
      }
      break;
    case '2':
      switch (abbr[2]) {
      case 'O':
        return water;
        break;
      }
      break;
    }
  case 'I':
    switch (abbr[1]) {
    case 'L':
      switch (abbr[2]) {
      case 'E':
        return 'I';
        break;
      }
    }
    break;
  case 'L':
    switch (abbr[1]) {
    case 'E':
      switch (abbr[2]) {
      case 'U':
        return 'L';
        break;
      }
      break;
    case 'Y':
      switch (abbr[2]) {
      case 'S':
        return 'K';
        break;
      }
      break;
    }
    break;
  case 'M':
    switch (abbr[1]) {
    case 'E':
      switch (abbr[2]) {
      case 'T':
        return 'M';
        break;
      }
    break;
    case 'S':
      switch (abbr[2]) {
      case 'E': // MSE (SELENOMETHIONINE)
        return 'M';
      }
    }
    break;
  case 'P':
    switch (abbr[1]) {
    case 'H':
      switch (abbr[2]) {
      case 'E':
        return 'F';
        break;
      }
      break;
    case 'R':
      switch (abbr[2]) {
      case 'O':
        return 'P';
        break;
      }
      break;
    }
    break;
  case 'S':
    switch (abbr[1]) {
    case 'E':
      switch (abbr[2]) {
      case 'R':
        return 'S';
        break;
      case 'C': // SEC (SELENOCYSTEINE)
        return 'U';
        break;
      }
      break;
    case 'O':                  /* SOL -- gromacs solvent residue */
      switch (abbr[2]) {
      case 'L':
        return water;
        break;
      }
      break;
    }
    break;
  case 'T':
    switch (abbr[1]) {
    case 'H':
      switch (abbr[2]) {
      case 'R':
        return 'T';
        break;
      }
      break;
    case 'I':
      switch (abbr[2]) {
      case 'P':
        return water;
        break;
      }
      break;
    case 'R':
      switch (abbr[2]) {
      case 'P':
        return 'W';
        break;
      }
      break;
    case 'Y':
      switch (abbr[2]) {
      case 'R':
        return 'Y';
        break;
      }
      break;
    }
    break;
  case 'V':
    switch (abbr[1]) {
    case 'A':
      switch (abbr[2]) {
      case 'L':
        return 'V';
        break;
      }
      break;
    }
    break;
  case 'W':
    switch (abbr[1]) {
    case 'A':
      switch (abbr[2]) {
      case 'T':
        return water;
        break;
      }
      break;
    }
    break;

  }

  return unknown;
}

static int SeekerFindColor(PyMOLGlobals * G, const AtomInfoType * ai, int n_more_plus_one)
{
  int result = ai->color;      /* default -- use first atom color */
  const AtomInfoType *ai0 = ai;
  while(1) {
    if(ai0->flags & cAtomFlag_guide)    /* best use guide color */
      return ai0->color;
    if(ai0->protons == cAN_C)   /* or use carbon color */
      result = ai0->color;
    n_more_plus_one--;
    if(n_more_plus_one > 0) {
      ai0++;
      if(!AtomInfoSameResidueP(G, ai, ai0))
        break;
    } else
      break;
  }
  return result;
}

static int SeekerFindTag(PyMOLGlobals * G, const AtomInfoType * ai, int sele, int codes,
                         int n_more_plus_one)
{
  int result = 0;      /* default -- no tag */
  const AtomInfoType *ai0 = ai;
  while(1) {
    int tag = SelectorIsMember(G, ai0->selEntry, sele);
    if(tag && (codes < 2) && (ai0->flags & cAtomFlag_guide))    /* use guide atom if present */
      return tag;
    if(result < tag) {
      if(!result)
        result = tag;
      else if((codes < 2) && (ai0->flags & cAtomFlag_guide))    /* residue based and on guide atom */
        result = tag;
    }
    n_more_plus_one--;
    if(n_more_plus_one > 0) {
      int do_break = false;
      ai0++;
      switch (codes) {
      case 0:
      case 1:
        if(!AtomInfoSameResidueP(G, ai, ai0))
          do_break = true;
        break;
      case 2:                  /* atoms */
        do_break = true;
        break;
      case 3:                  /* chains */
        if(!AtomInfoSameChainP(G, ai, ai0))
          do_break = true;
        break;
      }
      if(do_break)
        break;
    } else
      break;
  }
  return result;
}

void SeekerUpdate(PyMOLGlobals * G)
{
  /*  pymol::CObject *o = NULL;
     int s; */

  void *hidden = NULL;
  const AtomInfoType *ai;
  ObjectMolecule *obj;
  int nRow = 0;
  int label_mode = 0;
  int codes = 0;
  int max_row = 50;
  int default_color = 0;
  int align_sele = -1;          /* alignment selection */
  const int MAXCONSECUTIVEGAPS = 9;
  CSeqRow *row, *lab = NULL;
  std::vector<CSeqRow> row_vla;
  /* FIRST PASS: get all the residues represented properly */
  label_mode = SettingGetGlobal_i(G, cSetting_seq_view_label_mode);

  align_sele = ExecutiveGetActiveAlignmentSele(G);

  while (ExecutiveIterateObjectMolecule(G, &obj, &hidden)) {
    bool isHidden = (SettingGet<bool>(G, cSetting_hide_underscore_names) && (obj->Name[0] == '_'));
    if (obj->Enabled && (SettingGet<bool>(*obj, cSetting_seq_view)) && !isHidden) {
      int a;
      const AtomInfoType *last = NULL, *last_segi = NULL, *last_chain = NULL;
      const CoordSet *last_disc = NULL;
      int last_state;
      int last_abbr = true;
      int last_spacer = false;
      int nCol = 0;
      int nListEntries = 1;     /* first list starts at 1 always... */
      int est_col = obj->NAtom / 5 + 1;
      int est_char = obj->NAtom * 4;
      int first_atom_in_label;
      int missing_color = SettingGet_i(G, obj->Setting.get(), NULL, cSetting_seq_view_fill_color);
      const CoordSet *cs = obj->DiscreteFlag ? NULL : obj->getCoordSet(-2 /* current */);
      bool atom_in_state;

      int gapMode = SettingGet_i(G, obj->Setting.get(), nullptr, cSetting_seq_view_gap_mode);
      int min_pad = -1;
      CSeqCol *r1 = NULL, *l1 = NULL;   /* *col */

      if(nRow >= max_row)
        break;

      codes = SettingGet_i(G, obj->Setting.get(), NULL, cSetting_seq_view_format);
      if(obj->DiscreteFlag && SettingGet_b(G,
                                           obj->Setting.get(),
                                           NULL, cSetting_seq_view_discrete_by_state))
        codes = 4;
      default_color = SettingGet_i(G, obj->Setting.get(), NULL, cSetting_seq_view_color);

      /* allocate a row for labels, if present
         the text for the labels and the residues will line up exactly 
       */

      VecCheck(row_vla, nRow);
      if((label_mode == 2) || ((label_mode == 1) && (!nRow))) {
        lab = row_vla.data() + nRow++;
        lab->txt = pymol::vla<char>(est_char);
        lab->col = pymol::vla<CSeqCol>(est_col);
        lab->label_flag = true;
      } else {
        lab = NULL;
      }

      VecCheck(row_vla, nRow);

      row = row_vla.data() + nRow;
      if(lab)
        lab = row - 1;          /* critical! */
      row->txt = pymol::vla<char>(est_char);
      row->col = pymol::vla<CSeqCol>(est_col);
      row->fill = pymol::vla<CSeqCol>(est_col / 8);
      row->atom_lists = pymol::vla<int>(obj->NAtom + est_col + 1);
      row->atom_lists[0] = -1;  /* terminate the blank listQ (IMPORTANT!) */
      row->char2col = pymol::vla<int>(est_char);
      row->obj = obj;
      strcpy(row->name, obj->Name);
      row->color = obj->Color;
      ai = obj->AtomInfo;

      /* copy object name onto label row */

      if(lab) {

        int st_len;
        /* copy label text */

        VLACheck(lab->col, CSeqCol, nCol);
        l1 = lab->col + nCol;
        l1->start = lab->len;
        UtilConcatVLA(&lab->txt, &lab->len, "/");
        UtilConcatVLA(&lab->txt, &lab->len, obj->Name);
        l1->stop = lab->len;
        st_len = l1->stop - l1->start;

        if(label_mode == 2) {
          /* blank equivalent text for sequence row below the fixed label */
          VLACheck(row->col, CSeqCol, nCol);
          r1 = row->col + nCol;
          r1->start = row->len;
          UtilFillVLA(&row->txt, &row->len, ' ', st_len);
          r1->stop = row->len;
          r1->spacer = true;
          nCol++;
        }
      }
      if(label_mode < 2) {      /* no label rows, so put object name into left-hand column */

        /* copy label text */

        VLACheck(row->col, CSeqCol, nCol);
        r1 = row->col + nCol;
        r1->start = row->len;
        UtilConcatVLA(&row->txt, &row->len, "/");
        UtilConcatVLA(&row->txt, &row->len, obj->Name);
        r1->stop = row->len;
        r1->spacer = true;
        row->column_label_flag = true;
        row->title_width = row->len;
        nCol++;
      } else if(label_mode == 3) {      /* otherwise just insert a blank zero-length column */
        VLACheck(row->col, CSeqCol, nCol);
        r1 = row->col + nCol;
        r1->start = row->len;
        UtilConcatVLA(&row->txt, &row->len, "");
        r1->stop = row->len;
        r1->spacer = true;
        nCol++;
      }

      if(lab) {

        int st_len;
        /* copy label text */

        VLACheck(lab->col, CSeqCol, nCol);
        l1 = lab->col + nCol;
        l1->start = lab->len;
        if(obj->NAtom) {
          UtilConcatVLA(&lab->txt, &lab->len, "/");
          UtilConcatVLA(&lab->txt, &lab->len, LexStr(G, ai->segi));
          UtilConcatVLA(&lab->txt, &lab->len, "/");
          UtilConcatVLA(&lab->txt, &lab->len, LexStr(G, ai->chain));
          UtilConcatVLA(&lab->txt, &lab->len, "/");
        } else {
          UtilConcatVLA(&lab->txt, &lab->len, "///");
        }

        l1->stop = lab->len;
        st_len = l1->stop - l1->start;

        last_segi = ai;
        last_chain = ai;
        /* blank equivalent text for sequence row below the fixed label */
        VLACheck(row->col, CSeqCol, nCol);
        r1 = row->col + nCol;
        r1->start = row->len;
        UtilFillVLA(&row->txt, &row->len, ' ', st_len);
        r1->stop = row->len;
        r1->spacer = true;
        nCol++;
      } else {                  /* if no labels, just insert a space in row below */
        VLACheck(row->col, CSeqCol, nCol);
        r1 = row->col + nCol;
        r1->start = row->len;
        UtilConcatVLA(&row->txt, &row->len, " ");
        r1->stop = row->len;
        r1->spacer = true;
        nCol++;
      }

      last_state = -1;
      for(a = 0; a < obj->NAtom; a++) {
        first_atom_in_label = false;
        if(lab && !AtomInfoSameSegmentP(G, last_segi, ai)) {

          int st_len;

          if(row->len < min_pad) {
            row->len = min_pad;
          }
          min_pad = -1;

          /* copy label text */

          VLACheck(lab->col, CSeqCol, nCol);
          l1 = lab->col + nCol;
          l1->start = lab->len;
          UtilConcatVLA(&lab->txt, &lab->len, "/");
          UtilConcatVLA(&lab->txt, &lab->len, LexStr(G, ai->segi));
          UtilConcatVLA(&lab->txt, &lab->len, "/");
          UtilConcatVLA(&lab->txt, &lab->len, LexStr(G, ai->chain));
          UtilConcatVLA(&lab->txt, &lab->len, "/");
          l1->stop = lab->len;
          st_len = l1->stop - l1->start;

          /* blank equivalent text for sequence row */
          VLACheck(row->col, CSeqCol, nCol);
          r1 = row->col + nCol;
          r1->start = row->len;
          UtilFillVLA(&row->txt, &row->len, ' ', st_len);
          r1->stop = row->len;
          r1->spacer = true;
          nCol++;

          last_abbr = false;
          last_spacer = true;
          last_segi = ai;
          last_chain = ai;

        } else if(lab && !AtomInfoSameChainP(G, last_chain, ai)) {

          int st_len;

          if(row->len < min_pad) {
            row->len = min_pad;
          }
          min_pad = -1;

          /* copy label text */

          VLACheck(lab->col, CSeqCol, nCol);
          l1 = lab->col + nCol;
          l1->start = lab->len;
          UtilConcatVLA(&lab->txt, &lab->len, "/");
          UtilConcatVLA(&lab->txt, &lab->len, LexStr(G, ai->chain));
          UtilConcatVLA(&lab->txt, &lab->len, "/");
          l1->stop = lab->len;
          st_len = l1->stop - l1->start;

          /* blank equivalent text for sequence row */
          VLACheck(row->col, CSeqCol, nCol);
          r1 = row->col + nCol;
          r1->start = row->len;
          UtilFillVLA(&row->txt, &row->len, ' ', st_len);
          r1->stop = row->len;
          r1->spacer = true;
          nCol++;

          last_abbr = false;
          last_spacer = true;
          last_chain = ai;
        }

        if(min_pad < 0)
        {
          min_pad = row->len + 1; // + strlen(resi)
          for (int v = ai->resv; v; v /= 10) {
            min_pad++;
          }
          if (ai->inscode) {
            min_pad++;
          }
        }

        atom_in_state = (cs && cs->atmToIdx(a) >= 0);

        int gapsNeeded{0};

        if(gapMode != GapMode::NONE
           && AtomInfoSameChainP(G, ai, last)
           && (ai->flags & last->flags & cAtomFlag_polymer)
           && align_sele < 0){
            gapsNeeded = ai->resv - last->resv - 1;

            if(gapsNeeded > 1 && gapMode == GapMode::SINGLE){
              gapsNeeded = 1;
            }
        }

        auto push_gap = [&](const char* str)
          {
            auto str_size = strlen(str);
            UtilConcatVLA(&row->txt, &row->len, str);
            VLACheck(row->col, CSeqCol, nCol + str_size);
            r1 = row->col + nCol;
            for(int i = 0; i < str_size; i++){
                r1->color = missing_color;
                r1->spacer = true;
                r1->stop = r1->start + 1;
                auto lastStop = r1->stop;
                nCol++;
                r1 = row->col + nCol;
                r1->start = lastStop;
            }
          };
        switch (codes) {
        case 0:                /* one letter residue codes */
          if(!AtomInfoSameResidueP(G, last, ai)) {
            char abbr[2] = "1";
            last = ai;

            VLACheck(row->col, CSeqCol, nCol);
            r1 = row->col + nCol;
            r1->start = row->len;

            //Only include non-consecutive gaps when not doing alignment
            if(gapsNeeded > 0 && gapsNeeded <= MAXCONSECUTIVEGAPS){
              for(int g = 0; g < gapsNeeded; ++g){
                push_gap("-");
              }
            }
            else if(gapsNeeded > MAXCONSECUTIVEGAPS){
                push_gap("---...---");
            }

            if(obj->DiscreteFlag)
              r1->state = ai->discrete_state;

            first_atom_in_label = true;

            // single letter codes for polymer/solvent
            if (!(ai->flags & (cAtomFlag_organic | cAtomFlag_inorganic))) {
              abbr[0] = SeekerGetAbbr(G, LexStr(G, ai->resn), 'O', 0);
            } else {
              abbr[0] = 0;
            }

            r1->hint_no_space = last_abbr || last_spacer;

            if(!abbr[0]) {
              if(last_abbr) {
                UtilConcatVLA(&row->txt, &row->len, " ");
                r1->start = row->len;
              }

              if(ai->resn)
                UtilConcatVLA(&row->txt, &row->len, LexStr(G, ai->resn));
              else
                UtilConcatVLA(&row->txt, &row->len, "''");

              r1->stop = row->len;

              UtilConcatVLA(&row->txt, &row->len, " ");
            } else {
              UtilConcatVLA(&row->txt, &row->len, abbr);
              r1->is_abbr = true;
              r1->stop = row->len;
            }

            if(!atom_in_state)
              r1->color = missing_color;
            else if(default_color < 0)
              r1->color = SeekerFindColor(G, ai, obj->NAtom - a);
            else
              r1->color = default_color;
            if(align_sele >= 0) {
              r1->tag = SeekerFindTag(G, ai, align_sele, codes, obj->NAtom - a);
            } else {
              r1->tag = 0;
            }
            nCol++;
            last_abbr = abbr[0];
          }

          break;
        case 1:                /* explicit residue codes */
          if(!AtomInfoSameResidueP(G, last, ai)) {
            last = ai;

            VLACheck(row->col, CSeqCol, nCol);
            r1 = row->col + nCol;
            r1->start = row->len;
            if(obj->DiscreteFlag)
              r1->state = ai->discrete_state;
            first_atom_in_label = true;

            //Only include non-consecutive gaps when not doing alignment
            if(gapsNeeded > 0 && gapsNeeded <= MAXCONSECUTIVEGAPS){
              for(int g = 0; g < gapsNeeded; ++g){
                  push_gap("--- ");
              }
            }
            else if(gapsNeeded > MAXCONSECUTIVEGAPS){
                push_gap("---...--- ");
            }

            if(ai->resn)
              UtilConcatVLA(&row->txt, &row->len, LexStr(G, ai->resn));
            else
              UtilConcatVLA(&row->txt, &row->len, "''");
            r1->stop = row->len;

            if(!atom_in_state)
              r1->color = missing_color;
            else if(default_color < 0)
              r1->color = SeekerFindColor(G, ai, obj->NAtom - a);
            else
              r1->color = default_color;
            if(align_sele >= 0) {
              r1->tag = SeekerFindTag(G, ai, align_sele, codes, obj->NAtom - a);
            } else {
              r1->tag = 0;
            }
            UtilConcatVLA(&row->txt, &row->len, " ");
            nCol++;
          }
          break;
        case 2:                /* atom names */
          VLACheck(row->col, CSeqCol, nCol);
          r1 = row->col + nCol;
          r1->start = row->len;
          first_atom_in_label = true;
          if(ai->name)
            UtilConcatVLA(&row->txt, &row->len, LexStr(G, ai->name));
          else
            UtilConcatVLA(&row->txt, &row->len, "''");
          r1->stop = row->len;

          if(!atom_in_state)
            r1->color = missing_color;
          else if(default_color < 0)
            r1->color = ai->color;
          else
            r1->color = default_color;
          if(align_sele >= 0) {
            r1->tag = SeekerFindTag(G, ai, align_sele, codes, obj->NAtom - a);
          } else {
            r1->tag = 0;
          }
          UtilConcatVLA(&row->txt, &row->len, " ");
          nCol++;
          break;
        case 3:                /* chains */
          if(!AtomInfoSameChainP(G, last, ai)) {
            last = ai;

            VLACheck(row->col, CSeqCol, nCol);
            r1 = row->col + nCol;
            r1->start = row->len;
            first_atom_in_label = true;

            if(ai->chain)
              UtilConcatVLA(&row->txt, &row->len, LexStr(G, ai->chain));
            else
              UtilConcatVLA(&row->txt, &row->len, "''");
            r1->stop = row->len;
            if(default_color < 0)
              r1->color = SeekerFindColor(G, ai, obj->NAtom - a);
            else
              r1->color = default_color;
            if(align_sele >= 0) {
              r1->tag = SeekerFindTag(G, ai, align_sele, codes, obj->NAtom - a);
            } else {
              r1->tag = 0;
            }
            UtilConcatVLA(&row->txt, &row->len, " ");
            nCol++;
          }
          break;
        case 4:                /* state names */
          if(obj->DiscreteFlag) {
            CoordSet *cs;
            if((cs = obj->DiscreteCSet[a]) != last_disc) {
              last_disc = cs;
              if(cs) {
                default_color = SettingGet_i(G, cs->Setting.get(), obj->Setting.get(),
                                             cSetting_seq_view_color);
                VLACheck(row->col, CSeqCol, nCol);
                r1 = row->col + nCol;
                r1->start = row->len;
                r1->color = default_color;
                first_atom_in_label = true;
                if(cs->Name[0])
                  UtilConcatVLA(&row->txt, &row->len, cs->Name);
                else {
                  auto buf1 = pymol::string_format("%d", ai->discrete_state);
                  UtilConcatVLA(&row->txt, &row->len, buf1.c_str());
                }
                r1->stop = row->len;
                r1->state = ai->discrete_state;
                UtilConcatVLA(&row->txt, &row->len, " ");
                nCol++;
              }
            }
          } else {
            /* non-discrete objects simply get their states enumerated
               without selections */

            if(last_state < 0) {
              int b;
              const CoordSet *cs;
              last_state = 1;
              first_atom_in_label = true;
              for(b = 0; b < obj->NCSet; b++) {
                cs = obj->CSet[b];
                if(cs) {
                  default_color = SettingGet_i(G, cs->Setting.get(), obj->Setting.get(),
                                               cSetting_seq_view_color);

                  VLACheck(row->col, CSeqCol, nCol);
                  r1 = row->col + nCol;
                  r1->state = b + 1;
                  r1->start = row->len;
                  r1->atom_at = nListEntries + 1;       /* tricky & dangerous */
                  r1->color = default_color;
                  if(cs->Name[0])
                    UtilConcatVLA(&row->txt, &row->len, cs->Name);
                  else {
                    auto buf1 = pymol::string_format("%d", b + 1);
                    UtilConcatVLA(&row->txt, &row->len, buf1.c_str());
                  }
                  r1->stop = row->len;
                  UtilConcatVLA(&row->txt, &row->len, " ");
                  nCol++;
                }
              }
            }
          }
          break;
        case 5:                /* movie frames */
          break;
        }

        if(first_atom_in_label) {
          if(nCol > 1) {        /* terminate current list, if any */
            VLACheck(row->atom_lists, int, nListEntries);
            row->atom_lists[nListEntries] = -1;
            nListEntries++;
          }
          if(r1) {
            r1->atom_at = nListEntries;
          }
        }

        VLACheck(row->atom_lists, int, nListEntries);
        row->atom_lists[nListEntries] = a;
        nListEntries++;
        ai++;
      }

      if(lab) {
        /*        if(lab->len<row->len) {
           lab->len = row->len;
           } */
        VLASize(lab->txt, char, lab->len + 1);
        lab->txt[lab->len] = 0;
        VLACheck(lab->col, CSeqCol, nCol);      /* make sure we've got column records for labels too */
        lab->nCol = nCol;

        /*if(row->len<lab->len) {
           row->len = lab->len;
           } */
      }

      VLASize(row->txt, char, row->len + 1);
      row->txt[row->len] = 0;

      row->nCol = nCol;

      /* terminate last atom list */
      VLACheck(row->atom_lists, int, nListEntries);
      row->atom_lists[nListEntries] = -1;
      nListEntries++;
      nRow++;
    }
  }

  /* SECOND PASS: align columns to reflect current alignment and fixed labels */
  if(nRow) {
    int a, b;
    int nCol;
    int maxCol = 0;
    int done_flag = false;
    /* find out the maximum number of columns */

    for(a = 0; a < nRow; a++) {
      row = row_vla.data() + a;
      nCol = row->nCol;
      row->accum = 0;           /* initialize the accumulators */
      row->current = 0;
      if(maxCol < nCol)
        maxCol = nCol;
    }

    if(align_sele < 0) {
      /* in the simplest mode, just start each sequence in the same column */

      b = 0;
      while(!done_flag) {
        int max_offset = 0;
        done_flag = true;
        for(a = 0; a < nRow; a++) {
          row = row_vla.data() + a;
          if(!row->label_flag) {
            if(b < row->nCol) {
              CSeqCol *r1 = row->col + b;
              done_flag = false;

              r1->offset = r1->start + row->accum;
              if(max_offset < r1->offset)
                max_offset = r1->offset;
            }
          }
        }
        for(a = 0; a < nRow; a++) {
          row = row_vla.data() + a;
          if(!row->label_flag) {
            if(b < row->nCol) {
              CSeqCol *r1 = row->col + b;
              if(b < 3) {
                if(r1->offset < max_offset) {
                  row->accum += max_offset - r1->offset;
                }
              }
              r1->offset = r1->start + row->accum;
            }
          }
        }
        b++;
      }
    } else {
      /* in alignment mode, line up the tags */
      int stagger = false;
      /* intialize current columns and get the starting character */
      int current = 0;
      int first = true;
      switch (SettingGetGlobal_i(G, cSetting_seq_view_unaligned_mode)) {
      case 0:
      case 1:
      case 2:
        stagger = false;
        break;
      default:
        stagger = true;
        break;
      }
      for(a = 0; a < nRow; a++) {
        row = row_vla.data() + a;
        row->cCol = 0;
        if((!row->label_flag) && ((row->cCol < row->nCol))) {
          if(current < row->accum)
            current = row->accum;
        }
      }
      done_flag = false;
      while(!done_flag) {
        int hint_tagged_no_space = true;
        done_flag = true;
        {
          {
            /* insert untagged entries into their own columns */
            int untagged_flag = true;
            int saw_untagged_no_abbr = false;
            int hint_untagged_space = false;
            while(untagged_flag) {
              int space_added = false;
              int max_width = 0;
              untagged_flag = false;

              /* first get the spaces in... */

              for(a = 0; a < nRow; a++) {
                row = row_vla.data() + a;
                if((!row->label_flag) && (row->cCol < row->nCol)) {
                  CSeqCol *r1 = row->col + row->cCol;
                  if(!r1->tag) {        /* not aligned */
                    int text_len = (r1->stop - r1->start);
                    if((!first) && (!space_added) && (row->cCol > 2) &&
                       (codes || (((!r1->is_abbr) && (!r1->spacer))) ||
                        hint_untagged_space || (r1->is_abbr && (!r1->hint_no_space)))) {
                      /* insert space */
                      current++;
                      space_added = true;
                    }
                    if(max_width < text_len)
                      max_width = text_len;
                  }
                }
              }

              /* then do the rest */

              for(a = 0; a < nRow; a++) {
                row = row_vla.data() + a;
                if((!row->label_flag) && (row->cCol < row->nCol)) {
                  CSeqCol *r1 = row->col + row->cCol;
                  if(!r1->tag) {        /* not aligned */
                    int text_len = (r1->stop - r1->start);
                    untagged_flag = true;
                    done_flag = false;
                    saw_untagged_no_abbr |= (!r1->is_abbr) && (!r1->spacer);

                    first = false;
                    r1->offset = current;
                    r1->unaligned = true;

                    if(!r1->spacer) {
                      int aa;

                      for(aa = 0; aa < nRow; aa++) {    /* infill populate other rows with dashes */
                        if(aa != a) {
                          CSeqRow *row2 = row_vla.data() + aa;
                          if(!row2->label_flag) {
                            if(row2->cCol < row2->nCol) {
                              CSeqCol *r2 = row2->col + row2->cCol;
                              if(stagger || r2->tag || r2->spacer) {
                                VLACheck(row2->fill, CSeqCol, row2->nFill);
                                r2 = row2->fill + row2->nFill;
                                r2->stop = text_len;
                                r2->offset = current;
                                row2->nFill++;
                              }
                            } else {
                              CSeqCol *r2 = row2->col + row2->cCol;
                              VLACheck(row2->fill, CSeqCol, row2->nFill);
                              r2 = row2->fill + row2->nFill;
                              r2->stop = text_len;
                              r2->offset = current;
                              row2->nFill++;
                            }
                          }
                        }
                      }
                    }
                    if(stagger)
                      current += text_len;
                    else if(max_width < text_len)
                      max_width = text_len;
                  }
                }
              }
              if(!stagger)
                current += max_width;
              if(saw_untagged_no_abbr) {
                hint_untagged_space = true;
                hint_tagged_no_space = false;
              } else {
                hint_untagged_space = false;
                hint_tagged_no_space = true;
              }
              saw_untagged_no_abbr = false;
              for(a = 0; a < nRow; a++) {
                row = row_vla.data() + a;
                if((!row->label_flag) && (row->cCol < row->nCol)) {
                  CSeqCol *r1 = row->col + row->cCol;
                  if(!r1->tag) {
                    row->cCol++;
                  }
                }
              }
            }
          }
        }

        {
          /* next insert match-tagged entries into the same column */
          int min_tag = 0;
          for(a = 0; a < nRow; a++) {
            row = row_vla.data() + a;
            if((!row->label_flag) && (row->cCol < row->nCol)) {
              CSeqCol *r1 = row->col + row->cCol;
              if(r1->tag && ((min_tag > r1->tag) || (!min_tag))) {
                min_tag = r1->tag;
              }
            }
          }
          if(min_tag) {
            int width, max_width = 0;
            int space_added = false;
            int rep;
            for(rep = 0; rep < 2; rep++)
              for(a = 0; a < nRow; a++) {
                row = row_vla.data() + a;
                if((!row->label_flag) && (row->cCol < row->nCol)) {
                  CSeqCol *r1 = row->col + row->cCol;
                  if(r1->tag == min_tag) {
                    if((!first) && (!space_added) &&
                       (codes || (((!r1->is_abbr) && (!r1->spacer))) ||
                        (r1->is_abbr
                         && (!(r1->hint_no_space || hint_tagged_no_space))))) {
                      /* insert space */
                      current++;
                      space_added = true;
                    }
                    done_flag = false;
                    first = false;
                    r1->offset = current;
                    width = (r1->stop - r1->start);
                    if(max_width < width)
                      max_width = width;
                    /*   row->cCol++;  */
                  }
                }
              }
            {
              int aa;
              for(aa = 0; aa < nRow; aa++) {    /* infill populate other rows with dashes */
                CSeqRow *row2 = row_vla.data() + aa;
                if(!row2->label_flag) {
                  if(row2->cCol < row2->nCol) {
                    CSeqCol *r1 = row2->col + row2->cCol;
                    if(r1->tag != min_tag) {
                      CSeqCol *r2;
                      VLACheck(row2->fill, CSeqCol, row2->nFill);
                      r2 = row2->fill + row2->nFill;
                      r2->stop = max_width;
                      r2->offset = current;
                      row2->nFill++;
                    }
                  } else {
                    CSeqCol *r2;
                    VLACheck(row2->fill, CSeqCol, row2->nFill);
                    r2 = row2->fill + row2->nFill;
                    r2->stop = max_width;
                    r2->offset = current;
                    row2->nFill++;
                  }
                }
              }
            }
            for(a = 0; a < nRow; a++) {
              row = row_vla.data() + a;
              if((!row->label_flag) && (row->cCol < row->nCol)) {
                CSeqCol *r1 = row->col + row->cCol;
                if(r1->tag == min_tag) {
                  row->cCol++;
                }
              }
            }
            current += max_width;
          }
        }
      }
    }

    for(a = 0; a < nRow; a++) {
      row = row_vla.data() + a;
      nCol = row->nCol;
      if(row->label_flag)
        lab = row;
      else {
        for(b = 0; b < nCol; b++) {
          CSeqCol *r1 = row->col + b, *l1 = NULL;
          if(lab) {
            l1 = lab->col + b;  /* if a fixed label is present, 
                                   get the final offset from the residue line */
            if(l1->stop)
              l1->offset = r1->offset;
          }
        }
        lab = NULL;
      }
    }

  }

  /* THIRD PASS: fill in labels, based on actual residue spacing */

  if(nRow && (codes != 4)) {
    int a, b, c;
    int nCol;
    for(a = 0; a < nRow; a++) {
      lab = row_vla.data() + a;

      if(lab->label_flag) {
        int next_open = 0;
        int *atom_list;
        int st_len;
        int div, sub;
        int draw_it;
        int n_skipped = 0;
        int last_resv = -1;
        const AtomInfoType *last_ai = NULL;
        const ObjectMolecule *obj;
        const AtomInfoType *ai;
        row = lab + 1;
        nCol = row->nCol;
        obj = row->obj;
        div = SettingGet_i(G, obj->Setting.get(), NULL, cSetting_seq_view_label_spacing);
        sub = SettingGet_i(G, obj->Setting.get(), NULL, cSetting_seq_view_label_start);

        for(b = 0; b < nCol; b++) {
          CSeqCol *r1 = row->col + b;
          CSeqCol *l1 = lab->col + b;

          ai = NULL;
          if(r1->atom_at) {
            atom_list = row->atom_lists + r1->atom_at;
            if(*atom_list >= 0)
              ai = obj->AtomInfo + (*atom_list);        /* get first atom in list */
          }
          if(l1->stop) {        /* if label is already present, just line it up */
            l1->offset = r1->offset;
          } else if((r1->offset >= next_open) && ai) {
            if((div > 1) && (codes != 2)) {
              if(!((ai->resv - sub) % div))
                draw_it = true;
              else
                draw_it = false;
            } else {
              draw_it = true;
            }
            if(ai->resv != (last_resv + 1))     /* gap in sequence?  then draw label ASAP */
              draw_it = true;
            if(n_skipped >= (div + div))        /* don't skip too many without a label! */
              draw_it = true;

            if(AtomInfoSameResidueP(G, last_ai, ai))    /* don't ever draw a residue label twice */
              draw_it = false;

            if(draw_it) {
              n_skipped = 0;
              last_ai = ai;
              l1->start = lab->len;
              if(codes == 2) {
                UtilConcatVLA(&lab->txt, &lab->len, LexStr(G, ai->resn));
                UtilConcatVLA(&lab->txt, &lab->len, "`");
              }

              {
                char resi[8];
                AtomResiFromResv(resi, sizeof(resi), ai);
                UtilConcatVLA(&lab->txt, &lab->len, resi);
              }
              l1->stop = lab->len;
              st_len = l1->stop - l1->start + 1;
              l1->offset = r1->offset;
              next_open = r1->offset + st_len;

              /* make sure this label doesn't conflict with a fixed label */

              for(c = b + 1; c < nCol; c++) {
                CSeqCol *l2 = lab->col + c;
                if(l2->offset && (l2->offset < next_open)) {
                  l1->start = 0;
                  l1->stop = 0;
                  break;
                }
                if((c - b) > st_len)    /* only search as many columns as characters */
                  break;
              }
            } else
              n_skipped++;
          }

          if(ai)
            last_resv = ai->resv;
        }
      }
    }
  }

  /* FOURTH PASS: simply fill in character offsets */
  if(nRow) {
    int a, b;
    int nCol;
    int start, stop;
    for(a = 0; a < nRow; a++) {
      row = row_vla.data() + a;
      row->ext_len = 0;

      if(!row->label_flag) {
        nCol = row->nCol;

        for(b = 0; b < nCol; b++) {
          CSeqCol *r1 = row->col + b;
          stop = r1->offset + (r1->stop - r1->start);
          if(row->ext_len < stop)
            row->ext_len = stop;
        }
        VLACheck(row->char2col, int, row->ext_len);
        UtilZeroMem(row->char2col.data(), row->ext_len);
        for(b = 0; b < nCol; b++) {
          CSeqCol *r1 = row->col + b;
          int c;
          start = r1->offset;
          stop = r1->offset + (r1->stop - r1->start);
          for(c = start; c < stop; c++)
            row->char2col[c] = b + 1;
        }
      }
    }
  }
  SeqSetRow(G, std::move(row_vla), nRow);
  SeqSetHandler(G, G->Seeker);
}

int SeekerInit(PyMOLGlobals * G)
{
  CSeeker *I = NULL;
  if((I = (G->Seeker = new CSeeker()))) {
    I->dragInfo.row = -1;
    I->LastClickTime = UtilGetSeconds(G) - 1.0F;
    return 1;
  } else {
    return 0;
  }
}

void SeekerFree(PyMOLGlobals * G)
{
  DeleteP(G->Seeker);
}

CSeqRow* CSeeker::click(PyMOLGlobals * G, std::vector<CSeqRow>& rowVLA, int button, int row, int col,
                    int mod, int x, int y)
{
  return SeekerClick(G, rowVLA, button, row, col, mod, x, y);
}

CSeqRow* CSeeker::drag(PyMOLGlobals * G, std::vector<CSeqRow>& rowVLA, int row, int col, int mod)
{
  return SeekerDrag(G, rowVLA, row, col, mod);
}

CSeqRow* CSeeker::release(PyMOLGlobals * G, std::vector<CSeqRow>& rowVLA, int button, int row, int col,
                      int mod)
{
  return SeekerRelease(G, rowVLA, button, row, col, mod);
}

void CSeeker::refresh(PyMOLGlobals * G, std::vector<CSeqRow>& rowVLA)
{
  return SeekerRefresh(G, rowVLA);
}

void SeekerSetDragInfo(PyMOLGlobals* G, const SeekerDragInfo& dragInfo)
{
  G->Seeker->dragInfo = dragInfo;
}

SeekerDragInfo SeekerGetDragInfo(PyMOLGlobals* G)
{
  return G->Seeker->dragInfo;
}
