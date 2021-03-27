/*
 * Selector data types
 *
 * (c) 2016 Schrodinger, Inc.
 */

#pragma once

#include "os_std.h"
#include "pymol/memory.h"

#include "AtomIterators.h"
#include "Selector.h"
#include "ObjectMolecule.h"

#include "OVLexicon.h"
#include "OVOneToAny.h"

#define cNDummyModels 2
#define cNDummyAtoms 2

struct TableRec {
  int model;
  int atom;
  int index;
  float f1;
};

struct SelectionInfoRec {
  SelectorID_t ID = 0;
  std::string name;

  ObjectMolecule* theOneObject = nullptr;
  int theOneAtom = -1;

  bool justOneObject() const { return theOneObject != nullptr; }
  bool justOneAtom() const { return justOneObject() && theOneAtom >= 0; }

  SelectionInfoRec() = default;
  SelectionInfoRec(SelectorID_t id, std::string name_)
      : ID(id)
      , name(std::move(name_))
  {
  }
};


struct CSelectorManager
{
  std::vector<MemberType> Member;
  SelectorMemberOffset_t FreeMember = 0;
  std::vector<SelectionInfoRec> Info;
  SelectorID_t NSelection = 0;
  std::unordered_map<std::string, int> Key;
  CSelectorManager();
};

struct CSelector {
  PyMOLGlobals* G = nullptr;
  CSelectorManager* mgr = nullptr;
  std::vector<ObjectMolecule*> Obj;
  std::vector<TableRec> Table;
  pymol::cache_ptr<ObjectMolecule> Origin;
  pymol::cache_ptr<ObjectMolecule> Center;
  int NCSet = 0; // Seems to hold the largest NCSet in Obj
  bool SeleBaseOffsetsValid = false;
  CSelector(PyMOLGlobals* G, CSelectorManager* mgr);
  CSelector(const CSelector&) = default;
  CSelector& operator=(const CSelector&) = default;
  CSelector(CSelector&&) = default;
  CSelector& operator=(CSelector&&) = default;
  ~CSelector();
};

/**
 * Iterator over the selector table. If `SelectorUpdateTable(G,
 * cSelectorUpdateTableAllStates, -1)` was called, this would be all atoms.
 *
 * Does NOT provide coord or coordset access
 *
 * @pre Selector table is up-to-date
 */
class SelectorAtomIterator : public AbstractAtomIterator
{
  CSelector* selector;

public:
  int a; //!< index in selector table

  SelectorAtomIterator(CSelector* I)
      : selector(I)
  {
    reset();
  }

  void reset() override { a = cNDummyAtoms - 1; }

  bool next() override;
};
