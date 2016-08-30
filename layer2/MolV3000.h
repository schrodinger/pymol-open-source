/*
 * MOL/SDF V3000 input support for PyMOL.
 *
 * (c) Schrodinger, Inc.
 */

#include "PyMOLGlobals.h"
#include "AtomInfo.h"

const char * MOLV3000Parse(PyMOLGlobals * G,
    const char * buffer,
    AtomInfoType *& atInfo,
    BondType     *& bond,
    float        *& coord,
    int & nAtom,
    int & nBond)
  ;
