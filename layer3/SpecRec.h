#pragma once

#include"PyMOLObject.h"
#include"CGO.h"

/**
 * Specification record (a row in the object menu panel)
 */
class SpecRec {
public:
  /* NOTE: must zero-init with CALLOC */
  int type;
  WordType name;                /*only used for selections */
  pymol::CObject* obj;
  SpecRec *next;
  int visible;            /* This is actually when object is "Enabled", not visible */

  ObjectNameType group_name;

  /* not pickled */
  int sele_color;
  int hilight;                  /* 0 = none, 1 = name, 2 = group control (if any) */
  int previous;
  int cand_id;
  SpecRec *group;
  int group_member_list_id;
  int in_scene;
  int in_panel;
  int grid_slot;
  CGO *gridSlotSelIndicatorsCGO;

  void setEnabled(PyMOLGlobals* G, bool enabled);
  const char* baseName() const;
  bool isHidden(bool hide_underscore_names) const;
  bool isHiddenNotRecursive(bool hide_underscore_names) const;
  bool isChildOf(SpecRec const*) const;
};

