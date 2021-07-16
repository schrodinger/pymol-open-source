/*
 * Atom iterators
 *
 * (c) 2014 Schrodinger, Inc.
 */

#ifndef _H_ATOMITERATORS
#define _H_ATOMITERATORS

#include "PyMOLGlobals.h"

struct CoordSet;
struct ObjectMolecule;
struct AtomInfoType;

/**
 * Atom iterator base class
 *
 * The initial state of the iterator points "before" the first atom, there is no
 * current atom until next() is called the first time.
 */
class AbstractAtomIterator {
protected:
  int atm;      //!< Current atom index in object molecule
  int idx = -1; //!< Current coordinate index (atom index in coordset)

public:
  ObjectMolecule* obj;    //!< Current object molecule
  CoordSet* cs = nullptr; //!< Current coordiante set

  // virtual destructor for dynamic types
  virtual ~AbstractAtomIterator() = default;

  /// Resets the internal state to start over
  virtual void reset() = 0;

  /// Advance the internal state to the next atom.
  /// @return False if there is no next atom
  virtual bool next() = 0;

  /// Current atom
  AtomInfoType * getAtomInfo();
  const AtomInfoType * getAtomInfo() const;

  /// Current atom's coordinates
  float* getCoord();

  /// Current atom index in object molecule
  int getAtm() const {
    return atm;
  }

  /// Current coordinate index (atom index in coordset)
  int getIdx() const {
    return idx;
  }
};

/**
 * State specific iterator over an atom selection. Similar to `cmd.iterate_state`.
 * If state is -1 then iterate over all states.
 *
 * @verbatim
   SeleCoordIterator iter(G, sele, state);
   while(iter.next()) {
     dump3f(iter.getCoord(), "coords");
   }
   @endverbatim
 */
class SeleCoordIterator : public AbstractAtomIterator {
  PyMOLGlobals* G = nullptr;
  StateIndex_t statearg; // state argument, can be -1 (all), -2 (current), -3 (effective)
  StateIndex_t statemax; // largest state in selection
  bool per_object;              // whether to iterate over object states or global states
  ObjectMolecule * prev_obj;    // for per_object=true
  SelectorID_t sele = -1 /* cSelectionInvalid */;

public:
  int a;        //!< index in selection
  StateIndex_t state; //!< current state

  /// Undefined state until copy assigned from a valid state
  SeleCoordIterator() = default;
  SeleCoordIterator(PyMOLGlobals*, SelectorID_t sele_, StateIndex_t state_,
      bool update_table = true);

  void reset();
  bool next();

  /// Return true if iterating over all states
  bool isMultistate() const {
    return statearg == cStateAll;
  }

  /// Return true if iterating over all states of an object before advancing
  /// to the next object, rather than iterating over global states
  bool isPerObject() const {
    return per_object;
  }

  /// @see isPerObject()
  void setPerObject(bool per_object_) {
    per_object = per_object_ && isMultistate();
  }

private:
  bool nextStateInPrevObject();
};

/**
 * Iterator over an atom selection. Similar to `cmd.iterate`.
 *
 * Does NOT provide coord or coordset access
 *
 * @verbatim
   SeleAtomIterator iter(G, sele);
   while(iter.next()) {
     ai = iter.getAtomInfo();
   }
   @endverbatim
 */
class SeleAtomIterator : public AbstractAtomIterator {
  PyMOLGlobals * G;
  SelectorID_t sele; //!< selection
  char* stmp = nullptr; //!< temporary named selection

public:
  int a;        //!< index in selection

  /// Iterate over an existing atom selection
  /// @pre ::SelectorUpdateTable was called
  SeleAtomIterator(PyMOLGlobals* G_, SelectorID_t sele_)
      : G(G_)
      , sele(sele_)
  {
    reset();
  };

  SeleAtomIterator(PyMOLGlobals * G_, const char * sele_);
  ~SeleAtomIterator();

  void reset();
  bool next();
};

/**
 * CoordSet atom iterator, iterates over sorted atoms
 *
 * @verbatim
   CoordSetAtomIterator iter(G, cs);
   while(iter.next()) {
     dump3f(iter.getCoord(), "coords");
   }
   @endverbatim
 */
class CoordSetAtomIterator : public AbstractAtomIterator {
public:

  CoordSetAtomIterator(CoordSet * cs_);

  void reset() {
    atm = -1;
  }

  bool next();
};


#endif

// vi:sw=2:expandtab
