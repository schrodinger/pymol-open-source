/*
 * This file contains source code for the PyMOL computer program
 * Copyright (c) Schrodinger, LLC.
 *
 * Handling of "interned" strings in PyMOLGlobals::Lexicon
 */

#pragma once

#include "PyMOLGlobals.h"
#include "OVLexicon.h"
#include "pymol/zstring_view.h"

#define LexDec(G, i) OVLexicon_DecRef(G->Lexicon, i)
#define LexInc(G, i) OVLexicon_IncRef(G->Lexicon, i)
#define LexNumeric(i) (i)

/**
 * Get the pointer to the internal string storage for reference `i`.
 */
inline const char * LexStr(PyMOLGlobals * G, const lexidx_t & i) {
  return (i) ? OVLexicon_FetchCString(G->Lexicon, i) : "";
}

/**
 * Lookup or insert new string `s` into the global lexicon and return
 * it's numerical reference. Will always return 0 for the empty string.
 */
inline lexidx_t LexIdx(PyMOLGlobals * G, pymol::zstring_view s) {
  return s && !s.empty() ? OVLexicon_GetFromCString(G->Lexicon, s.c_str()).word : 0;
}

/**
 * Assignment (i = j) with reference count update for old and new value.
 */
inline void LexAssign(PyMOLGlobals * G, lexidx_t& i, const lexidx_t& j) {
  if (i != j) {
    LexDec(G, i);
    i = j;
    LexInc(G, i);
  }
}

inline void LexAssign(PyMOLGlobals * G, lexidx_t& i, const char * s) {
  LexDec(G, i);
  i = LexIdx(G, s);
}

#define LEX_BORROW_NOTFOUND -1

/**
 * Lookup string `s` without inserting or incrementing the ref count. If
 * `s` is not in the lexicon, return -1.
 */
inline lexidx_t LexBorrow(PyMOLGlobals * G, const char * s) {
  auto result = OVLexicon_BorrowFromCString(G->Lexicon, s);
  return (result.status == OVstatus_SUCCESS) ? result.word : LEX_BORROW_NOTFOUND;
}
