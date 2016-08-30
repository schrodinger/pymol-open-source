/*
 * Selector data types
 *
 * (c) 2016 Schrodinger, Inc.
 */

#pragma once

#include "os_std.h"

#include "Selector.h"
#include "ObjectMolecule.h"

#include "OVLexicon.h"
#include "OVOneToAny.h"

#define SelectorWordLength 1024
typedef char SelectorWordType[SelectorWordLength];

#define cNDummyModels 2
#define cNDummyAtoms 2

struct TableRec {
  int model;
  int atom;
  int index;
  float f1;
};

struct SelectionInfoRec {
  int ID;
  int justOneObjectFlag;
  ObjectMolecule *theOneObject;
  int justOneAtomFlag;
  int theOneAtom;
};

struct CSelector {
  MemberType *Member;           /* Must be first in structure, so that we can get this w/o knowing the struct */
  SelectorWordType *Name;       /* this seems rather excessive, since name len < ObjNameMax ... */
  SelectionInfoRec *Info;
  int NSelection, NActive;
  int TmpCounter;
  int NMember;
  int FreeMember;
  ObjectMolecule **Obj;
  TableRec *Table;
  float *Vertex;
  int *Flag1, *Flag2;
  ov_size NAtom;
  ov_size NModel;
  int NCSet;
  int SeleBaseOffsetsValid;
  ObjectMolecule *Origin, *Center;
  OVLexicon *Lex;
  OVOneToAny *Key;
  OVOneToOne *NameOffset;
};
