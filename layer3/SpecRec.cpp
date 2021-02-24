/*
 * (c) Schrodinger, Inc.
 */

#include "SpecRec.h"

#include <cstring>

/**
 * If name is "<group_name>.<base_name>" then return <base_name>, otherwise
 * just return name.
 */
const char* SpecRec::baseName() const
{
  if (auto const pos = std::strlen(group_name)) {
    if (std::strncmp(name, group_name, pos) == 0 && name[pos] == '.') {
      return name + pos + 1;
    }
  }

  return name;
}

/**
 * True if this record is a direct or indirect child of `other`.
 *
 * @pre ExecutiveUpdateGroups()
 */
bool SpecRec::isChildOf(SpecRec const* other) const
{
  return group == other || (group && group->isChildOf(other));
}

/**
 * True if this record is hidden from the object-menu-panel.
 *
 * @param hide_underscore_names Whether names that start with an underscore
 * should be hidden
 *
 * @pre ExecutiveUpdateGroups()
 */
bool SpecRec::isHidden(bool hide_underscore_names) const
{
  if (hide_underscore_names) {
    for (auto rec = this; rec; rec = rec->group) {
      if (rec->baseName()[0] == '_')
        return true;
    }
  }

  return false;
}

/**
 * Like isHidden() but assume that the parent group is not hidden.
 *
 * @pre Parent group is not hidden
 */
bool SpecRec::isHiddenNotRecursive(bool hide_underscore_names) const
{
  assert(!group || !group->isHidden(hide_underscore_names));

  return hide_underscore_names && baseName()[0] == '_';
}
