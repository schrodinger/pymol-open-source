
/* 
A* -------------------------------------------------------------------
B* This file contains source code for the PyMOL computer program
C* copyright 1998-2006 by Warren Lyford Delano of DeLano Scientific. 
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

#include <set>

#include"os_python.h"

#include"os_predef.h"
#include"os_std.h"
#include"os_gl.h"

#include"OOMac.h"
#include"ObjectAlignment.h"
#include"Base.h"
#include"MemoryDebug.h"
#include"CGO.h"
#include"Scene.h"
#include"Setting.h"
#include"PConv.h"
#include"main.h"
#include"Color.h"
#include"Executive.h"
#include"OVContext.h"
#include"Util.h"
#include"Selector.h"
#include"Seq.h"
#include"Seeker.h"
#include"ShaderMgr.h"
#include "Lex.h"

/* 
   just so you don't forget...

   alignment vlas are zero-separated lists of groups of unique atom identifiers 

*/

static int GroupOrderKnown(PyMOLGlobals * G,
                           int *curVLA,
                           const int *newVLA,
                           int cur_start,
                           int new_start, ObjectMolecule * guide, int *action)
{
  int order_known = false;
  if(guide) {
    int c, id;
    int cur_offset = -1;
    int new_offset = -1;

    /* find lowest offset within the cur group */
    c = cur_start;
    while((id = curVLA[c++])) {
      auto eoo = ExecutiveUniqueIDAtomDictGet(G, id);
      if (eoo && eoo->obj == guide) {
        if((cur_offset < 0) || (eoo->atm < cur_offset))
          cur_offset = eoo->atm;
      }
    }

    /* find lowest offset within the new group */
    c = new_start;
    while((id = newVLA[c++])) {
      auto eoo = ExecutiveUniqueIDAtomDictGet(G, id);
      if (eoo && eoo->obj == guide) {
        if((new_offset < 0) || (eoo->atm < new_offset))
          new_offset = eoo->atm;
      }
    }

    if((new_offset >= 0) && (cur_offset >= 0)) {
      if(new_offset < cur_offset) {
        order_known = true;
        *action = -1;
      } else if(new_offset > cur_offset) {
        order_known = true;
        *action = 1;
      }
    }
  }
  return order_known;
}

static int AlignmentFindTag(PyMOLGlobals * G, AtomInfoType * ai, int sele,
                            int n_more_plus_one)
{
  int result = 0;      /* default -- no tag */
  AtomInfoType *ai0 = ai;
  while(1) {
    int tag = SelectorIsMember(G, ai0->selEntry, sele);
    if(tag && (ai0->flags & cAtomFlag_guide))   /* use guide atom if present */
      return tag;
    if(result < tag) {
      if(!result)
        result = tag;
      else if(ai0->flags & cAtomFlag_guide)     /* residue based and on guide atom */
        result = tag;
    }
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

/**
 * Get single letter abbreviation for CLUSTAL output.
 *
 * See also:
 * pymol.exporting._resn_to_aa
 * AtomInfoKnownNucleicResName()
 * AtomInfoKnownProteinResName()
 * SeekerGetAbbr()
 */
static char get_abbr(PyMOLGlobals * G, const AtomInfoType * ai) {
  const char * resn = LexStr(G, ai->resn);
  const char unknown = ((ai->flags & cAtomFlag_polymer)) ? '?' : 0;

  if ((ai->flags & cAtomFlag_nucleic)) {
    if (resn[0] == 'D') {
      ++resn;
    }

    if (strlen(resn) != 1) {
      return unknown;
    }

    return resn[0];
  }

  return SeekerGetAbbr(G, resn, 0, unknown);
}

int ObjectAlignmentAsStrVLA(PyMOLGlobals * G, ObjectAlignment * I, int state, int format,
                            char **str_vla)
{
  int ok = true;
  ov_size len = 0;
  char *vla = VLAlloc(char, 1000);
  int force_update = false;
  int active_only = false;
  int max_name_len = 12;        /* default indentation */

  if(state < 0)
    state = I->getCurrentState();
  if(state < 0)
    state = SceneGetState(G);
  if(state >= 0 && state < I->getNFrame()) {
    ObjectAlignmentState *oas = I->State.data() + state;
    if(oas->alignVLA) {
      if(state != I->SelectionState) {  /* get us a selection for the current state */
        I->ForceState = state;
        force_update = true;
        I->update();
      }

      switch (format) {
      case 0:                  /* aln */
        UtilConcatVLA(&vla, &len, "CLUSTAL\n\n");
        break;
      }

      {
        int align_sele = SelectorIndexByName(G, I->Name);
        if(align_sele >= 0) {
          int nRow = 0;
          ov_size nCol = 0;
          CSeqRow *row_vla = NULL, *row;
          char *cons_str = NULL;
          void *hidden = NULL;

          ObjectMolecule *obj;

          {
            row_vla = VLACalloc(CSeqRow, 10);

            /* first, find out which objects are included in the
               alignment and count the name length */

            while(ExecutiveIterateObjectMolecule(G, &obj, &hidden)) {
              if((obj->Enabled || !active_only) && (obj->Name[0] != '_')) {
                int a;
                const AtomInfoType *ai = obj->AtomInfo.data();
                for(a = 0; a < obj->NAtom; a++) {
                  if(SelectorIsMember(G, ai->selEntry, align_sele)) {
                    int name_len = strlen(obj->Name);
                    if(max_name_len < name_len)
                      max_name_len = name_len;
                    VLACheck(row_vla, CSeqRow, nRow);
                    row = row_vla + nRow;
                    row->obj = obj;
                    row->nCol = obj->NAtom;
                    nRow++;
                    break;
                  }
                  ai++;
                }
              }
            }

            /* next, figure out how many total columns exist */

            {
              int done = false;
              while(!done) {
                int a;
                int min_tag = -1;
                int untagged_col = false;
                done = true;
                for(a = 0; a < nRow; a++) {
                  row = row_vla + a;
                  while(row->cCol < row->nCol) {        /* advance to next tag in each row & find lowest */
                    AtomInfoType *ai = row->obj->AtomInfo + row->cCol;
                    done = false;
                    if(AtomInfoSameResidueP(G, row->last_ai, ai)) {
                      row->cCol++;
                    } else if(!get_abbr(G, ai)) {      /* not a polymer residue */
                      row->cCol++;
                    } else {
                      int tag =
                        AlignmentFindTag(G, ai, align_sele, row->nCol - row->cCol);
                      if(tag) { /* we're at a tagged atom... */
                        if(min_tag > tag)
                          min_tag = tag;
                        else if(min_tag < 0)
                          min_tag = tag;
                        break;
                      } else {
                        untagged_col = true;
                        break;
                      }
                    }
                  }
                  if(untagged_col)
                    break;
                }
                if(untagged_col) {
                  nCol++;
                  /* increment all untagged atoms */
                  for(a = 0; a < nRow; a++) {
                    row = row_vla + a;
                    if(row->cCol < row->nCol) {
                      AtomInfoType *ai = row->obj->AtomInfo + row->cCol;
                      int tag =
                        AlignmentFindTag(G, ai, align_sele, row->nCol - row->cCol);
                      if(!tag) {
                        row->last_ai = ai;
                        row->cCol++;
                      }
                    }
                  }
                } else if(min_tag >= 0) {
                  /* increment all matching tagged atoms */
                  nCol++;
                  for(a = 0; a < nRow; a++) {
                    row = row_vla + a;
                    if(row->cCol < row->nCol) {
                      AtomInfoType *ai = row->obj->AtomInfo + row->cCol;
                      int tag =
                        AlignmentFindTag(G, ai, align_sele, row->nCol - row->cCol);
                      if(tag == min_tag) {      /* advance past this tag */
                        row->cCol++;
                        row->last_ai = ai;
                      }
                    }
                  }
                }
              }
            }
            /* allocate storage for the sequence alignment */

            cons_str = pymol::calloc<char>(nCol + 1);  /* conservation string */

            {
              int a;
              for(a = 0; a < nRow; a++) {
                row = row_vla + a;
                row->txt = pymol::vla<char>(nCol + 1);
                row->len = 0;
                row->last_ai = NULL;
                row->cCol = 0;
              }
            }

            nCol = 0;

            {
              int done = false;
              while(!done) {
                int a;
                int min_tag = -1;
                int untagged_col = false;
                done = true;
                for(a = 0; a < nRow; a++) {
                  row = row_vla + a;
                  while(row->cCol < row->nCol) {        /* advance to next tag in each row & find lowest */
                    AtomInfoType *ai = row->obj->AtomInfo + row->cCol;
                    done = false;
                    if(AtomInfoSameResidueP(G, row->last_ai, ai)) {
                      row->cCol++;
                    } else if(!get_abbr(G, ai)) {      /* not a polymer residue */
                      row->cCol++;
                    } else {
                      int tag =
                        AlignmentFindTag(G, ai, align_sele, row->nCol - row->cCol);
                      if(tag) { /* we're at a tagged atom... */
                        if(min_tag > tag)
                          min_tag = tag;
                        else if(min_tag < 0)
                          min_tag = tag;
                        break;
                      } else {
                        untagged_col = true;
                        break;
                      }
                    }
                  }
                }
                if(untagged_col) {
                  /* increment all untagged atoms */
                  for(a = 0; a < nRow; a++) {
                    row = row_vla + a;
                    if(row->cCol < row->nCol) {
                      AtomInfoType *ai = row->obj->AtomInfo + row->cCol;
                      int tag =
                        AlignmentFindTag(G, ai, align_sele, row->nCol - row->cCol);
                      if(!tag) {
                        if(!AtomInfoSameResidueP(G, row->last_ai, ai)) {
                          row->last_ai = ai;
                          row->txt[row->len] = get_abbr(G, ai);
                        } else {
                          row->txt[row->len] = '-';
                        }
                        row->len++;
                        row->cCol++;
                      } else {
                        row->txt[row->len] = '-';
                        row->len++;
                      }
                    } else {
                      row->txt[row->len] = '-';
                      row->len++;
                    }
                  }
                  cons_str[nCol] = ' ';
                  nCol++;
                } else if(min_tag >= 0) {
                  char cons_abbr = ' ';
                  int abbr_cnt = 0;
                  /* increment all matching tagged atoms */
                  for(a = 0; a < nRow; a++) {
                    row = row_vla + a;
                    if(row->cCol < row->nCol) {
                      AtomInfoType *ai = row->obj->AtomInfo + row->cCol;
                      int tag =
                        AlignmentFindTag(G, ai, align_sele, row->nCol - row->cCol);
                      if(tag == min_tag) {      /* advance past this tag */
                        char abbr;
                        if(!AtomInfoSameResidueP(G, row->last_ai, ai)) {
                          row->last_ai = ai;
                          abbr = get_abbr(G, ai);
                          if(cons_abbr == ' ')
                            cons_abbr = abbr;
                          else if(cons_abbr != abbr)
                            cons_abbr = 0;
                          abbr_cnt++;
                        } else {
                          abbr = '-';
                          cons_abbr = 0;
                        }
                        row->txt[row->len] = abbr;
                        row->len++;
                        row->cCol++;
                      } else {
                        cons_abbr = 0;
                        row->txt[row->len] = '-';
                        row->len++;
                      }
                    } else {
                      cons_abbr = 0;
                      row->txt[row->len] = '-';
                      row->len++;
                    }
                  }
                  if(abbr_cnt > 1) {
                    if(cons_abbr)
                      cons_str[nCol] = '*';     /* aligned and identical */
                    else
                      cons_str[nCol] = '.';     /* aligned but not identical */
                  } else
                    cons_str[nCol] = ' ';
                  nCol++;
                }
              }
            }

            {
              int block_width = 76 - (max_name_len + 1);
              int done = false;
              ov_size seq_len = 0;
              int a;
              while(!done) {
                done = true;
                for(a = 0; a < nRow; a++) {
                  row = row_vla + a;
                  UtilNPadVLA(&vla, &len, row->obj->Name, max_name_len + 1);
                  if(seq_len < row->len) {
                    UtilNPadVLA(&vla, &len, row->txt + seq_len, block_width);
                  }
                  UtilConcatVLA(&vla, &len, "\n");
                }
                if(seq_len < nCol) {
                  UtilNPadVLA(&vla, &len, "", max_name_len + 1);
                  UtilNPadVLA(&vla, &len, cons_str + seq_len, block_width);
                  UtilConcatVLA(&vla, &len, "\n");
                }
                seq_len += block_width;
                for(a = 0; a < nRow; a++) {
                  row = row_vla + a;
                  if(seq_len < row->len) {
                    done = false;
                    break;
                  }
                }
                UtilConcatVLA(&vla, &len, "\n");
              }
            }
          }

          /* free up resources */
          if(row_vla) {
            int a;
            for(a = 0; a < nRow; a++) {
              row = row_vla + a;
              row->txt.freeP();
            }
          }
          FreeP(cons_str);
          VLAFreeP(row_vla);
        }
      }
    }
  }

  if(force_update) {
    I->update();
  }

  VLASize(vla, char, len + 1);
  vla[len] = 0;
  *str_vla = vla;
  return ok;
}

static int *AlignmentMerge(PyMOLGlobals * G, int *curVLA, const int *newVLA,
                           ObjectMolecule * guide, ObjectMolecule * flush)
{
  /* curVLA and newVLA must be properly sized and zero terminated... */
  int *result = NULL;
  int n_result = 0;

  {
    {

      int n_cur = VLAGetSize(curVLA);
      int n_new = VLAGetSize(newVLA);

      // get the set of non-guide objects in the new alignment, they need
      // to be flushed (removed) from the current alignment
      std::set<const ObjectMolecule*> flushobjects;
      if (flush) {
        flushobjects.insert(flush);
      } else {
        for (const int *it = newVLA, *it_end = it + n_new; it != it_end; ++it) {
          if (*it) {
            auto eoo = ExecutiveUniqueIDAtomDictGet(G, *it);
            if (eoo && eoo->obj != guide) {
              flushobjects.insert(eoo->obj);
            }
          }
        }
      }

      /* first, go through and eliminate old matching atoms between guide and flush (if any) */
      {
        int cur_start = 0;
        while(cur_start < n_cur) {

          while((cur_start < n_cur) && !curVLA[cur_start]) {
            cur_start++;
          }

          {
            int other_seen = 0;
            int flush_seen = false;
            ObjectMolecule *obj;

            {
              int cur = cur_start;
              int id;
              while((id = curVLA[cur])) {
                auto eoo = ExecutiveUniqueIDAtomDictGet(G, id);
                if (eoo) {
                  obj = eoo->obj;
                  if(flushobjects.count(obj)) {
                    flush_seen = true;
                  } else {
                    other_seen++;
                  }
                }
                cur++;
              }
            }

            if(flush_seen) {    /* eliminate flush atoms */
              int cur = cur_start;
              int id;
              while((id = curVLA[cur])) {
                auto eoo = ExecutiveUniqueIDAtomDictGet(G, id);
                if (eoo) {
                  obj = eoo->obj;
                  if(flushobjects.count(obj)) {
                    int tmp = cur;
                    while(curVLA[tmp]) {
                      curVLA[tmp] = curVLA[tmp + 1];
                      tmp++;
                    }
                  }
                }
                cur++;
              }
            }

            if(other_seen < 2) {     /* eliminate orphaned atoms */
              int cur = cur_start;
              while(curVLA[cur])
                curVLA[cur++] = 0;
            }

            while(curVLA[cur_start])
              cur_start++;
            while((cur_start < n_cur) && !curVLA[cur_start]) {
              cur_start++;
            }
          }
        }
      }

      /* now combine the alignments */

      {
        OVOneToAny *used = OVOneToAny_New(G->Context->heap);
        OVOneToAny *active = OVOneToAny_New(G->Context->heap);
        int cur_start = 0;
        int new_start = 0;

        result = VLAlloc(int, ((n_cur < n_new) ? n_new : n_cur));

        while((cur_start < n_cur) || (new_start < n_new)) {

          int action;           /* -1 = insert new, 0 = merge, 1 = insert cur */

          /* make sure both lists are queued up on the next identifier */

          while((cur_start < n_cur) && !curVLA[cur_start])
            cur_start++;
          while((new_start < n_new) && !newVLA[new_start])
            new_start++;

          if(newVLA[new_start]) {       /* default is to insert new first... */
            action = -1;
          } else {
            action = 1;
          }

          if((cur_start < n_cur) && (new_start < n_new) &&
             curVLA[cur_start] && newVLA[new_start]) {
            /* both lists active */

            int c, id;
            int overlapping = false;

            OVOneToAny_Reset(active);
            c = cur_start;
            while((id = curVLA[c++])) { /* record active atoms */
              OVOneToAny_SetKey(active, id, 1);
            }

            c = new_start;
            while((id = newVLA[c++])) { /* see if there are any matches */
              if(OVreturn_IS_OK(OVOneToAny_GetKey(active, id))) {
                overlapping = true;
                break;
              }
            }

            if(overlapping) {
              /* overlapping, so merge */
              action = 0;

            } else {
              /* non-overlapping, so we need to figure out which goes first... */
              if(!GroupOrderKnown(G, curVLA, newVLA,
                                  cur_start, new_start, guide, &action)) {
                int c, id;
                ObjectMolecule *obj, *last_obj = NULL;
                c = cur_start;
                while((id = curVLA[c++])) {
                  auto eoo = ExecutiveUniqueIDAtomDictGet(G, id);
                  if (eoo) {
                    obj = eoo->obj;
                    if(obj != last_obj) {
                      if(GroupOrderKnown(G, curVLA, newVLA,
                                         cur_start, new_start, obj, &action))
                        break;
                      else
                        last_obj = obj;
                    }
                  }
                }

                /* if order isn't set by now, then it doesn't matter...
                   so group new will go in first */
              }
            }
          }

          /* check assumptions */

          if((action < 1) && !(new_start < n_new))
            action = 1;
          else if((action > (-1)) && !(cur_start < n_cur))
            action = -1;

          /* take action */

          {
            int id;

            switch (action) {
            case -1:           /* insert new */
              if(new_start < n_new) {
                while((id = newVLA[new_start])) {
                  if(OVOneToAny_GetKey(used, id).status == OVstatus_NOT_FOUND) {
                    if(OVreturn_IS_OK(OVOneToAny_SetKey(used, id, 1))) {
                      VLACheck(result, int, n_result);
                      result[n_result] = id;
                      n_result++;
                    }
                  }
                  new_start++;
                }
                while((new_start < n_new) && (!newVLA[new_start]))
                  new_start++;
              }
              VLACheck(result, int, n_result);
              result[n_result] = 0;
              n_result++;
              break;
            case 0:            /* merge, with cur going first */
              if(new_start < n_new) {
                while((id = newVLA[new_start])) {
                  if(OVOneToAny_GetKey(used, id).status == OVstatus_NOT_FOUND) {
                    if(OVreturn_IS_OK(OVOneToAny_SetKey(used, id, 1))) {
                      VLACheck(result, int, n_result);
                      result[n_result] = id;
                      n_result++;
                    }
                  }
                  new_start++;
                }
                while((new_start < n_new) && (!newVLA[new_start]))
                  new_start++;
              }
              if(cur_start < n_cur) {
                while((id = curVLA[cur_start])) {
                  if(OVOneToAny_GetKey(used, id).status == OVstatus_NOT_FOUND) {
                    if(OVreturn_IS_OK(OVOneToAny_SetKey(used, id, 1))) {
                      VLACheck(result, int, n_result);
                      result[n_result] = id;
                      n_result++;
                    }
                  }
                  cur_start++;
                }
                while((cur_start < n_cur) && (!curVLA[cur_start]))
                  cur_start++;
              }
              VLACheck(result, int, n_result);
              result[n_result] = 0;
              n_result++;
              break;
            case 1:            /* insert cur */
              if(cur_start < n_cur) {
                while((id = curVLA[cur_start])) {
                  if(OVOneToAny_GetKey(used, id).status == OVstatus_NOT_FOUND) {
                    if(OVreturn_IS_OK(OVOneToAny_SetKey(used, id, 1))) {
                      VLACheck(result, int, n_result);
                      result[n_result] = id;
                      n_result++;
                    }
                  }
                  cur_start++;
                }
                while((cur_start < n_cur) && (!curVLA[cur_start]))
                  cur_start++;
                VLACheck(result, int, n_result);
                result[n_result] = 0;
                n_result++;
              }
              break;
            }
          }
        }
        OVOneToAny_DEL_AUTO_NULL(active);
        OVOneToAny_DEL_AUTO_NULL(used);
      }
    }
  }
  if(result && n_result && (!result[n_result - 1])) {
    VLACheck(result, int, n_result);
    result[n_result] = 0;
    n_result++;
  }
  VLASize(result, int, n_result);
  return result;
}

static PyObject *ObjectAlignmentStateAsPyList(ObjectAlignmentState * I)
{
  PyObject *result = NULL;

  result = PyList_New(2);
  if(I->alignVLA) {
    PyList_SetItem(result, 0, PConvIntVLAToPyList(I->alignVLA));
  } else {
    PyList_SetItem(result, 0, PConvAutoNone(NULL));
  }
  PyList_SetItem(result, 1, PyString_FromString(I->guide));
  return (PConvAutoNone(result));
}

static PyObject *ObjectAlignmentAllStatesAsPyList(ObjectAlignment * I)
{

  PyObject *result = NULL;
  int a;
  result = PyList_New(I->getNFrame());
  for(a = 0; a < I->getNFrame(); a++) {
    PyList_SetItem(result, a, ObjectAlignmentStateAsPyList(I->State.data() + a));
  }
  return (PConvAutoNone(result));

}

static int ObjectAlignmentStateFromPyList(PyMOLGlobals * G, ObjectAlignmentState * I,
                                          PyObject * list, int version)
{
  int ok = true;
  int ll = 0;
  if(ok)
    ok = (list != NULL);
  if(ok)
    ok = PyList_Check(list);
  if(ok)
    ll = PyList_Size(list);
  /* TO SUPPORT BACKWARDS COMPATIBILITY...
     Always check ll when adding new PyList_GetItem's */
  if(ok && (ll > 1)) {
    PConvPyListToIntVLA(PyList_GetItem(list, 0), &I->alignVLA);
    strcpy(I->guide, PyString_AsString(PyList_GetItem(list, 1)));

    for (auto& align : I->alignVLA) {
      if (align) {
        align = SettingUniqueConvertOldSessionID(G, align);
      }
    }
  }
  return (ok);
}

static int ObjectAlignmentAllStatesFromPyList(ObjectAlignment * I, PyObject * list,
                                              int version)
{
  int ok = true;
  int a;
  if(ok)
    ok = PyList_Check(list);
  if(ok) {
    int nstates = PyList_Size(list);
    I->State.resize(nstates);
    for(a = 0; a < nstates; a++) {
      auto *val = PyList_GetItem(list, a);
      ok =
        ObjectAlignmentStateFromPyList(I->G, I->State.data() + a, val,
                                       version);
      if(!ok)
        break;
    }
  }
  return (ok);
}

int ObjectAlignmentNewFromPyList(PyMOLGlobals * G, PyObject * list,
                                 ObjectAlignment ** result, int version)
{
  int ok = true;
  ObjectAlignment *I = NULL;
  (*result) = NULL;
  if(ok)
    ok = (list != NULL);
  if(ok)
    ok = PyList_Check(list);

  I = new ObjectAlignment(G);
  if(ok)
    ok = (I != NULL);

  if(ok){
    auto *val = PyList_GetItem(list, 0);
    ok = ObjectFromPyList(G, val, I);
  }
  if(ok)
    ok = ObjectAlignmentAllStatesFromPyList(I, PyList_GetItem(list, 2), version);
  if(ok) {
    (*result) = I;
    ObjectAlignmentRecomputeExtent(I);
  } else {
    /* cleanup? */
  }
  return (ok);
}

PyObject *ObjectAlignmentAsPyList(ObjectAlignment * I)
{
  PyObject *result = NULL;

  result = PyList_New(3);
  PyList_SetItem(result, 0, ObjectAsPyList(I));
  PyList_SetItem(result, 1, PyInt_FromLong(I->getNFrame()));
  PyList_SetItem(result, 2, ObjectAlignmentAllStatesAsPyList(I));

  return (PConvAutoNone(result));
}


/*========================================================================*/

void ObjectAlignmentRecomputeExtent(ObjectAlignment * I)
{
  float mx[3], mn[3];
  int extent_flag = false;
  int a;
  for(a = 0; a < I->getNFrame(); a++)
    if(I->State[a].primitiveCGO) {
      if(CGOGetExtent(I->State[a].primitiveCGO.get(), mn, mx)) {
        if(!extent_flag) {
          extent_flag = true;
          copy3f(mx, I->ExtentMax);
          copy3f(mn, I->ExtentMin);
        } else {
          max3f(mx, I->ExtentMax, I->ExtentMax);
          min3f(mn, I->ExtentMin, I->ExtentMin);
        }
      }
    }
  I->ExtentFlag = extent_flag;
}


/*========================================================================*/
void ObjectAlignment::update()
{
  auto I = this;
  int update_needed = false;
  {
    int a;
    for(a = 0; a < getNFrame(); a++) {
      ObjectAlignmentState *oas = I->State.data() + a;
      if(!oas->valid){
        update_needed = true;
      }
    }
  }
  if(update_needed) {
    {
      int a;
      for(a = 0; a < getNFrame(); a++) {
        ObjectAlignmentState *oas = I->State.data() + a;
	if(!oas->valid){
          ObjectMolecule *guide_obj = NULL;
          if(oas->guide[0]) {
            guide_obj = ExecutiveFindObjectMoleculeByName(G, oas->guide);
          }
          if(I->SelectionState == a)
            I->SelectionState = -1;

          oas->primitiveCGO.reset();

          oas->id2tag.clear();

          {
            CGO *cgo = CGONew(G);

            if(oas->alignVLA) {
              int id, b = 0, c;
              auto& vla = oas->alignVLA;
              int n_id = vla.size();
              float mean[3], vert[3], gvert[3];
              int n_coord = 0;
              int tag = SELECTOR_BASE_TAG + 1;
              auto& id2tag = oas->id2tag;

              CGOBegin(cgo, GL_LINES);

              while(b < n_id) {

                int gvert_valid;
                while((b < n_id) && (!vla[b]))
                  b++;

                if(!(b < n_id))
                  break;

                c = b;
                n_coord = 0;
                gvert_valid = false;
                zero3f(mean);
                while((id = vla[c++])) {
                  auto eoo = ExecutiveUniqueIDAtomDictGet(G, id);
                  if (eoo) {
                    if(ObjectMoleculeGetAtomVertex(eoo->obj, a,
                                                   eoo->atm, vert)) {
                      n_coord++;
                      add3f(vert, mean, mean);
                      if(eoo->obj == guide_obj) {
                        copy3f(vert, gvert);
                        gvert_valid = true;
                      }
                    }
                  }
                }

                if(n_coord > 2) {       /* >2 points, then draw to mean or guide vertex */
                  float scale = 1.0F / n_coord;

                  scale3f(mean, scale, mean);

                  c = b;
                  while((id = vla[c++])) {
                    auto eoo = ExecutiveUniqueIDAtomDictGet(G, id);
                    if (eoo) {
                      if(ObjectMoleculeGetAtomVertex(eoo->obj, a,
                                                     eoo->atm, vert)) {
                        if(gvert_valid) {
                          if(eoo->obj != guide_obj) {
                            cgo->add<cgo::draw::line>(gvert, vert);
                          }
                        } else {
                          cgo->add<cgo::draw::line>(mean, vert);
                        }
                      }
                    }
                  }
                } else if(n_coord) {    /* if 2 points, then simply draw a line */
                  float first[3];
                  int first_flag = true;
                  c = b;
                  while((id = vla[c++])) {
                    auto eoo = ExecutiveUniqueIDAtomDictGet(G, id);
                    if (eoo) {
                      if(ObjectMoleculeGetAtomVertex(eoo->obj, a,
                                                     eoo->atm, vert)) {
                        if(first_flag) {
                          copy3f(vert, first);
                          first_flag = false;
                        } else {
                          cgo->add<cgo::draw::line>(first, vert);
                        }
                      }
                    }
                  }
                }
                /* update the it2tag dictionary */

                tag++;

                while((b < n_id) && vla[b]) {
                  id2tag[vla[b]] = tag;
                  b++;
                }
              }
              CGOEnd(cgo);
            }

            CGOStop(cgo);
            oas->primitiveCGO.reset(cgo);
            if (!CGOHasOperationsOfType(oas->primitiveCGO.get(), cgo::draw::line::op_code)){
              oas->primitiveCGO.reset();
            }
          }
          oas->valid = true;
        }
      }
    }
  }
  if(I->SelectionState < 0) {
    int state = -1;
    if(I->ForceState >= 0) {
      state = I->ForceState;
      I->ForceState = 0;
    } else {
      state = I->getCurrentState();
    }
    // TODO do these fallbacks make any sense?
    if(state < 0)
      state = SceneGetState(G);
    if(state >= I->getNFrame())
      state = I->getNFrame() - 1;
    if(state < 0)
      state = 0;
    if(state < I->getNFrame()) {
      ObjectAlignmentState *oas = I->State.data() + state;
      if(!oas->id2tag.empty()) {
        SelectorDelete(G, I->Name);
        SelectorCreateFromTagDict(G, I->Name, oas->id2tag, false);
        I->SelectionState = state;
      }
    }
  }
  SceneInvalidate(I->G);
}


/*========================================================================*/

int ObjectAlignment::getNFrame() const
{
  return State.size();
}


/*========================================================================*/

void ObjectAlignment::render(RenderInfo * info)
{
  auto I = this;
  int state = info->state;
  CRay *ray = info->ray;
  auto pick = info->pick;
  const RenderPass pass = info->pass;
  ObjectAlignmentState *sobj = NULL;
  const float *color;

  ObjectPrepareContext(I, info);

  color = ColorGet(G, I->Color);

  if (pick)
    return;

  if(pass == RenderPass::Opaque || ray) {
    if((I->visRep & cRepCGOBit)) {

      for(StateIterator iter(G, I->Setting.get(), state, I->getNFrame()); iter.next();) {
        sobj = I->State.data() + iter.state;

        if (!sobj->primitiveCGO)
          continue;

	if(ray) {
	    CGORenderRay(sobj->primitiveCGO.get(), ray, info, color, NULL, I->Setting.get(), NULL);
	} else if(G->HaveGUI && G->ValidContext) {
#ifndef PURE_OPENGL_ES_2
	  if(!info->line_lighting)
	    glDisable(GL_LIGHTING);
#endif
	  SceneResetNormal(G, true);
          bool use_shader = SettingGetGlobal_b(G, cSetting_use_shaders);

          CGO * cgo = NULL;

          if (use_shader) {
            bool as_cylinders =
              SettingGetGlobal_b(G, cSetting_alignment_as_cylinders) &&
              SettingGetGlobal_b(G, cSetting_render_as_cylinders);

            bool trilines = !as_cylinders && SettingGetGlobal_b(G, cSetting_trilines);

            if (sobj->renderCGO && (
                  (as_cylinders ^ sobj->renderCGO_has_cylinders) ||
                  (trilines ^ sobj->renderCGO_has_trilines))){
              sobj->renderCGO.reset();
            }

            if (!sobj->renderCGO) {
              int shader =
                as_cylinders    ? GL_CYLINDER_SHADER :
                trilines        ? GL_TRILINES_SHADER : GL_LINE_SHADER;

              CGO *tmpCGO = CGONew(G), *tmp2CGO = NULL;
              CGOEnable(tmpCGO, shader);
              CGOSpecial(tmpCGO, SET_ALIGNMENT_UNIFORMS_ATTRIBS);

              if (as_cylinders) {
                tmp2CGO = CGOConvertLinesToCylinderShader(sobj->primitiveCGO.get(), tmpCGO, false);
              } else if (trilines) {
                tmp2CGO = CGOConvertToTrilinesShader(sobj->primitiveCGO.get(), tmpCGO, false);
              } else {
                tmp2CGO = CGOConvertToLinesShader(sobj->primitiveCGO.get(), tmpCGO, false);
              }

              tmpCGO->free_append(tmp2CGO);

              CGODisable(tmpCGO, shader);

              sobj->renderCGO.reset(tmpCGO);
              sobj->renderCGO_has_cylinders = as_cylinders;
              sobj->renderCGO_has_trilines = trilines;
            }

            cgo = sobj->renderCGO.get();
          } else {
            cgo = sobj->primitiveCGO.get();
          }

          if (cgo) {
            CGORenderGL(cgo, color, I->Setting.get(), NULL, info, NULL);
          }

#ifndef PURE_OPENGL_ES_2
	  glEnable(GL_LIGHTING);
#endif
	}
      }
    }
  }
}

void ObjectAlignment::invalidate(cRep_t rep, cRepInv_t level, int state)
{
  if((rep == cRepAll) || (rep == cRepCGO)) {
    for(StateIterator iter(G, Setting.get(), state, getNFrame()); iter.next();) {
      ObjectAlignmentState& sobj = State[iter.state];
      sobj.valid = false;
      sobj.renderCGO.reset();
    }
  }
}


/*========================================================================*/
ObjectAlignment::ObjectAlignment(PyMOLGlobals * G) : pymol::CObject(G)
{
  type = cObjectAlignment;
}


/*========================================================================*/
ObjectAlignment *ObjectAlignmentDefine(PyMOLGlobals * G,
                                       ObjectAlignment * obj,
                                       const pymol::vla<int>& align_vla,
                                       int state,
                                       int merge,
                                       ObjectMolecule * guide, ObjectMolecule * flush)
{
  ObjectAlignment *I = NULL;

  if(obj) {
    if(obj->type != cObjectAlignment)       /* TODO: handle this */
      obj = NULL;
  }
  if(!obj) {
    I = new ObjectAlignment(G);
  } else {
    I = obj;
    I->invalidate(cRepAll, cRepInvRep, state);
  }
  if(state < 0)
    state = I->getNFrame();

  VecCheck(I->State, state);

  {
    ObjectAlignmentState *oas = I->State.data() + state;
    oas->valid = false;
    if(guide) {
      strcpy(oas->guide, guide->Name);
    }
    if(align_vla.data()) {
      if(merge && oas->alignVLA) {
        int *new_vla = AlignmentMerge(G, oas->alignVLA.data(), align_vla.data(), guide, flush);
        if(new_vla) {
          oas->alignVLA = pymol::vla_take_ownership(new_vla);
        }
      } else {
        oas->alignVLA = align_vla;
      }
    } else {
      VLAFreeP(oas->alignVLA);
    }
  }
  if(I) {
    ObjectAlignmentRecomputeExtent(I);
  }
  SceneChanged(G);
  SceneCountFrames(G);
  return (I);
}

pymol::CObject* ObjectAlignment::clone() const
{
  return new ObjectAlignment(*this);
}
