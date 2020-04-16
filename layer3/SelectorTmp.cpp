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

pymol::Result<SelectorTmp> SelectorTmp::make(PyMOLGlobals* G, const char* sele)
{
  SelectorTmp self;
  self.m_G = G;
  auto res = SelectorGetTmpResult(G, sele, self.m_name);
  if (res) {
    self.m_count = res.result();
    return std::move(self);
  }
  return res.error_move();
}

pymol::Result<SelectorTmp2> SelectorTmp2::make(
    PyMOLGlobals* G, const char* sele)
{
  SelectorTmp2 self;
  self.m_G = G;
  auto res = SelectorGetTmp2Result(G, sele, self.m_name);
  if (res) {
    self.m_count = res.result();
    return std::move(self);
  }
  return res.error_move();
}
