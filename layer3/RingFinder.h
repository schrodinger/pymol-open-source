#pragma once

#include <cstddef>
#include <vector>

struct AtomInfoType;
struct ObjectMolecule;

/**
 * Ring finder
 */
class AbstractRingFinder
{
  ObjectMolecule* m_obj = nullptr;
  std::vector<int> m_indices;

  void recursion(int atm, int depth);

protected:
  AbstractRingFinder() = delete;
  AbstractRingFinder(unsigned maxringsize)
      : m_indices(maxringsize)
  {
  }

  /**
   * Optional object preparation.
   */
  virtual void prepareObject(ObjectMolecule* obj) {}

  /**
   * Optional atom filter.
   * @return True to exclude an atom, false to include it.
   */
  virtual bool atomIsExcluded(const AtomInfoType&) const { return false; }

  /**
   * Is called when a ring is found.
   *
   * @param obj Molecule object
   * @param indices Ring atom indices
   * @param size Ring size
   */
  virtual void onRingFound(
      ObjectMolecule* obj, const int* indices, size_t size) = 0;

public:
  /*
   * Does a depth-first search for all paths of length in
   * range [3, maxringsize], which lead back to `atm` and
   * don't visit any atom twice.
   */
  void apply(ObjectMolecule* obj, int atm);
};
