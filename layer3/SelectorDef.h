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
  int ID = 0;
  bool justOneObjectFlag = false;
  ObjectMolecule* theOneObject = nullptr;
  bool justOneAtomFlag = false;
  int theOneAtom = 0;
  SelectionInfoRec() = default;
  SelectionInfoRec(int id) : ID(id) {}
};


struct CSelectorManager
{
  std::vector<MemberType> Member;
  int NMember = 0;
  int FreeMember = 0;
  std::vector<std::string> Name;
  std::vector<SelectionInfoRec> Info;
  static int TmpCounter;
  int NSelection = 0;
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

