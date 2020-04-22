/*
 * RAII for temporary selections
 *
 * (c) 2020 Schrodinger, Inc.
 */

#include "Selector.h"

#include <utility>

SelectorTmp::SelectorTmp(SelectorTmp&& other)
{
  *this = std::move(other);
  assert(!other.m_name[0]);
  assert(other.m_count == -1);
}

/**
 * Factory which propagates errors.
 * @param sele Selection expression
 * @param empty_is_error If true, then an empty expression is an error,
 * otherwise not (but getIndex() returns cSelectionInvalid)
 */
pymol::Result<SelectorTmp> SelectorTmp::make(
    PyMOLGlobals* G, const char* sele, bool empty_is_error)
{
  if (empty_is_error && !sele[0]) {
    return pymol::Error("Empty expression");
  }

  SelectorTmp self;
  self.m_G = G;
  auto res = SelectorGetTmpResult(G, sele, self.m_name);
  if (res) {
    assert(!empty_is_error || self.m_name[0]);
    self.m_count = res.result();
    return std::move(self);
  }
  return res.error_move();
}

/**
 * Factory which propagates errors.
 * @param sele Object name pattern or selection expression
 */
pymol::Result<SelectorTmp2> SelectorTmp2::make(
    PyMOLGlobals* G, const char* sele, bool empty_is_error)
{
  if (empty_is_error && !sele[0]) {
    return pymol::Error("Empty expression");
  }

  SelectorTmp2 self;
  self.m_G = G;
  auto res = SelectorGetTmp2Result(G, sele, self.m_name);
  if (res) {
    assert(!empty_is_error || self.m_name[0]);
    self.m_count = res.result();
    return std::move(self);
  }
  return res.error_move();
}
