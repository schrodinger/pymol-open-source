/*
 * MOL/SDF V3000 input support for PyMOL.
 *
 * (c) Schrodinger, Inc.
 */

struct PyMOLGlobals;
struct AtomInfoType;
struct BondType;

const char * MOLV3000Parse(PyMOLGlobals * G,
    const char * buffer,
    AtomInfoType *& atInfo,
    BondType     *& bond,
    float        *& coord,
    int & nAtom,
    int & nBond)
  ;
