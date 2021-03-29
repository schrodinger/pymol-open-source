/*
 * (c) Schrodinger, Inc.
 */

#pragma once

struct ObjectMolecule;

/**
 * Proxy for ObjectMolecule::Neighbor entries of a given atom.
 */
class AtomNeighbors
{
  const int* m_neighbor;

public:
  struct value_type {
    // Member layout must match ObjectMolecule::Neighbor.
    // The data layout is described in ObjectMolecule::getNeighborArray.
    int atm;  ///< Atom index
    int bond; ///< Bond index
  };

  AtomNeighbors(const ObjectMolecule* I, int atm);

  /// Get the i'th neighbor
  const value_type& operator[](size_t i) const { return *(begin() + i); }

  /// Number of neighbors
  unsigned size() const { return m_neighbor[0]; }

  const value_type* begin() const
  {
    return reinterpret_cast<value_type const*>(m_neighbor + 1);
  }

  const value_type* end() const { return begin() + size(); }
};
