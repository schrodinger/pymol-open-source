/*
 * Selector data types
 *
 * (c) 2016 Schrodinger, Inc.
 */

#pragma once

#include "os_std.h"
#include "pymol/memory.h"

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
  pymol::copyable_ptr<ObjectMolecule> Origin;
  pymol::copyable_ptr<ObjectMolecule> Center;
  int NCSet = 0; // Seems to hold the largest NCSet in Obj
  bool SeleBaseOffsetsValid = false;
  CSelector(PyMOLGlobals* G, CSelectorManager* mgr);
  CSelector(CSelector&&) = default;
  CSelector& operator=(CSelector&&) = default;
  ~CSelector();
};

