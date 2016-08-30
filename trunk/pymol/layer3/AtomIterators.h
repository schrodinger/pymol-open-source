/*
 * Atom iterators
 *
 * (c) 2014 Schrodinger, Inc.
 */

#ifndef _H_ATOMITERATORS
#define _H_ATOMITERATORS

#include "PyMOLGlobals.h"
#include "ObjectMolecule.h"
#include "CoordSet.h"
#include "AtomInfo.h"

/*
 * Atom iterator base class
 */
class AbstractAtomIterator {
protected:
  int atm;      // atom index in object molecule
  int idx;      // atom index in coordset

public:
  ObjectMolecule * obj;
  CoordSet * cs;

  // virtual destructor for dynamic types
  virtual ~AbstractAtomIterator() {}

  // resets the internal state to start over
  virtual void reset() = 0;

  // advance the internal state to the next atom, return false if there is no
  // next atom
  virtual bool next() = 0;

  // get current atom
  AtomInfoType * getAtomInfo() {
    return obj->AtomInfo + atm;
  };

  // get current atom's coordinates
  float * getCoord() {
    return cs->Coord + (3 * idx);
  };

  // get current atom index in object molecule
  int getAtm() {
    return atm;
  }

  // get current coordinate index (atom index in coordset)
  int getIdx() {
    return idx;
  }
};

/*
 * State specific iterator over an atom selection. Similar to cmd.iterate_state.
 * If state is -1 then iterate over all states.
 *
 * SeleCoordIterator iter(G, sele, state);
 * while(iter.next()) {
 *   dump3f(iter.getCoord(), "coords");
 * }
 */
class SeleCoordIterator : public AbstractAtomIterator {
  PyMOLGlobals * G;
  int statearg; // state argument, can be -1 (all), -2 (current), -3 (effective)
  int statemax; // largest state in selection
  bool per_object;              // whether to iterate over object states or global states
  ObjectMolecule * prev_obj;    // for per_object=true

public:
  int a;        // index in selection
  int state;    // current state

  void init(PyMOLGlobals * G_, int sele_, int state_);

  SeleCoordIterator() {} // undefined state until "init()" called
  SeleCoordIterator(PyMOLGlobals * G_, int sele_, int state_) {
    init(G_, sele_, state_);
  };

  void reset();
  bool next();

  // return true if iterating over all states
  bool isMultistate() const {
    return statearg == -1;
  }

  // return true if iterating over all states of an object before advancing
  // to the next object, rather than iterating over global states
  bool isPerObject() const {
    return per_object;
  }

  void setPerObject(bool per_object_) {
    per_object = per_object_ && isMultistate();
  }

private:
  bool nextStateInPrevObject() {
    if (prev_obj && (++state) < prev_obj->NCSet) {
      a = prev_obj->SeleBase - 1;
      return true;
    }
    return false;
  }
};

/*
 * Iterator over an atom selection. Similar to cmd.iterate.
 *
 * Does NOT provide coord or coordset access
 *
 * SeleAtomIterator iter(G, sele);
 * while(iter.next()) {
 *   ai = iter.getAtomInfo();
 * }
 */
class SeleAtomIterator : public AbstractAtomIterator {
  PyMOLGlobals * G;
  int sele;     // selection
  char * stmp;  // temporary named selection

public:
  int a;        // index in selection

  SeleAtomIterator(PyMOLGlobals * G_, int sele_) {
    G = G_;
    sele = sele_;
    stmp = NULL;
    reset();
  };

  SeleAtomIterator(PyMOLGlobals * G_, const char * sele_);
  ~SeleAtomIterator();

  void reset();
  bool next();
};

/*
 * CoordSet atom iterator, iterates over sorted atoms
 *
 * CoordSetAtomIterator iter(G, cs);
 * while(iter.next()) {
 *   dump3f(iter.getCoord(), "coords");
 * }
 */
class CoordSetAtomIterator : public AbstractAtomIterator {
public:

  CoordSetAtomIterator(CoordSet * cs_) {
    cs = cs_;
    obj = cs->Obj;
    reset();
  }

  void reset() {
    atm = -1;
  }

  bool next();
};


#endif

// vi:sw=2:expandtab
