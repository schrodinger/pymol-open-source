/*
 * MAE format export helper functions
 *
 * (c) 2018 Schrodinger, Inc.
 */

#include <string>

#include "os_std.h"

#include "MaeExportHelpers.h"
#include "Color.h"
#include "Lex.h"
#include "Setting.h"
#include "Executive.h"
#include "SpecRec.h"

/**
 * Setting getter
 */
template <typename V>
V SettingGet(PyMOLGlobals * G, const SeleCoordIterator& iter, int index) {
  V value = SettingGet<V>(G, iter.cs->Setting.get(), iter.obj->Setting.get(), index);
  return AtomSettingGetWD(G, iter.getAtomInfo(), index, value);
}

/*
 * enums from
 * mmshare/include/mmct.h
 * mmshare/include/mmctg.h
 */

enum MM_CTAtomStyle {
  MMCT_ATOM_NOSTYLE = 0,
  MMCT_ATOM_CIRCLE,
  MMCT_ATOM_CPK,
  MMCT_ATOM_BALLNSTICK
};

enum MM_CTBondStyle {
  MMCT_BOND_NOSTYLE = 0,
  MMCT_BOND_WIRE,
  MMCT_BOND_TUBE,
  MMCT_BOND_BALLNSTICK
};

enum MM_CTRibbonStyle {
  MMCT_RIBBON_STYLE_NONE      = 0,  //
  MMCT_RIBBON_STYLE_CARTOON   = 1,  // auto
  MMCT_RIBBON_STYLE_RIBBON    = 2,  // auto
  MMCT_RIBBON_STYLE_TUBE      = 3,  // tube, cartoon_tube_radius=0.30
  MMCT_RIBBON_STYLE_THINTUBE  = 4,  // tube, cartoon_tube_radius=0.15
  MMCT_RIBBON_STYLE_CURVELINE = 5,  // ribbon, ribbon_sampling=5
  MMCT_RIBBON_STYLE_CALINE    = 6,  // ribbon, ribbon_as_cylinders=0
  MMCT_RIBBON_STYLE_CATUBE    = 7,  // ribbon, ribbon_as_cylinders=1
};

/**
 * Get the MM_CTAtomStyle for the current atom of the iterator
 */
int MaeExportGetAtomStyle(PyMOLGlobals * G,
    const SeleCoordIterator& iter)
{
  auto ai = iter.getAtomInfo();

  if (ai->visRep & cRepSphereBit)
    return MMCT_ATOM_CPK;

  if ((ai->visRep & cRepNonbondedSphereBit) && !ai->bonded)
    return MMCT_ATOM_BALLNSTICK;

  if ((ai->visRep & cRepCylBit) && ai->bonded
      && SettingGet<bool>(G, iter, cSetting_stick_ball)
      && SettingGet<float>(G, iter, cSetting_stick_ball_ratio) > 1.0f)
    return MMCT_ATOM_BALLNSTICK;

  return MMCT_ATOM_NOSTYLE;
}

/**
 * Get the MM_CTBondStyle for the bond between two atoms
 */
int MaeExportGetBondStyle(const AtomInfoType * ai1, const AtomInfoType * ai2) {
  if (ai1->visRep & ai2->visRep & cRepCylBit)
    return MMCT_BOND_TUBE;

  // assuming line_stick_helper=on
  if ((ai1->visRep & (cRepLineBit | cRepCylBit)) &&
      (ai2->visRep & (cRepLineBit | cRepCylBit)))
    return MMCT_BOND_WIRE;

  return MMCT_BOND_NOSTYLE;
}

/**
 * Get the MM_CTRibbonStyle for an atom
 */
int MaeExportGetRibbonStyle(const AtomInfoType * ai) {
  if (ai->visRep & cRepCartoonBit) {
    switch (ai->cartoon) {
      case cCartoon_skip:
        return MMCT_RIBBON_STYLE_NONE;
      case cCartoon_loop:
      case cCartoon_tube:
      case cCartoon_putty:
        return MMCT_RIBBON_STYLE_TUBE;
    }
    return MMCT_RIBBON_STYLE_CARTOON;
  }

  if (ai->visRep & cRepRibbonBit) {
    return MMCT_RIBBON_STYLE_CALINE;
  }

  return MMCT_RIBBON_STYLE_NONE;
}

/**
 * Get the MAE ribbon color index and the hex formatted RGB ribbon color for
 * the current atom of the iterator
 */
void MaeExportGetRibbonColor(PyMOLGlobals * G,
    const SeleCoordIterator& iter,
    char * ribbon_color_rgb)
{
  const auto ai = iter.getAtomInfo();

  if (!(ai->flags & cAtomFlag_guide))
    return;

  int setting_index;
  if (ai->visRep & cRepCartoonBit) {
    setting_index = cSetting_cartoon_color;
  } else if (ai->visRep & cRepRibbonBit) {
    setting_index = cSetting_ribbon_color;
  } else {
    return;
  }

  int color = SettingGet<int>(G, iter, setting_index);

  if (color > 0) {
    const float * rgb = ColorGet(G, color);
    sprintf(ribbon_color_rgb, "%02X%02X%02X",
        int(rgb[0] * 255),
        int(rgb[1] * 255),
        int(rgb[2] * 255));
  }
}

/**
 * Get the MAE label user text for an atom.
 *
 * Quotes and backslashes are escaped.
 */
std::string MaeExportGetLabelUserText(PyMOLGlobals * G,
    const AtomInfoType * ai)
{
  std::string label_user_text;

  if (ai->label) {
    const char * label = LexStr(G, ai->label);
    for (const char * p = label; *p; ++p) {
      if (*p == '"' || *p == '\\')
        label_user_text += '\\';
      label_user_text += *p;
    }
  }

  return label_user_text;
}

/**
 * Get the MAE group title/id
 */
std::string MaeExportGetSubGroupId(PyMOLGlobals * G,
    const pymol::CObject * obj)
{
  std::string subgroupid;
  const SpecRec * rec = nullptr;

  // obj -> spec rec
  for (ObjectIterator iter(G); iter.next();) {
    if (iter.getObject() == obj) {
      rec = iter.getSpecRec();
      break;
    }
  }

  // "->".join(grouphierarchy)
  for (; rec && rec->group_name[0]; rec = rec->group) {
    if (!subgroupid.empty()) {
      subgroupid.insert(0, "->");
    }
    subgroupid.insert(0, rec->group_name);
  }

  return subgroupid;
}

/**
 * Get parsable string representation, with quotes and escaped
 * quotes/backslashes if needed.
 */
std::string MaeExportStrRepr(const char * text)
{
  if (text[0] /* not empty string */) {
    bool needquotes = false;

    // check accepted ascii characters
    for (const char * p = text; *p; ++p) {
      if (*p < '$' || *p > 'z' || *p == '\\') {
        needquotes = true;
        break;
      }
    }

    if (!needquotes) {
      return text;
    }
  }

  std::string quoted_text;
  quoted_text.reserve(strlen(text) + 2);
  quoted_text += '"';

  for (const char * p = text; *p; ++p) {
    if (*p == '"' || *p == '\\')
      quoted_text += '\\';
    quoted_text += *p;
  }

  quoted_text += '"';
  return quoted_text;
}
