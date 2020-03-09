#include "RingFinder.h"
#include "ObjectMolecule.h"

void AbstractRingFinder::recursion(int atm, int depth)
{
  int atm_neighbor, j;

  m_indices[depth] = atm;

  ITERNEIGHBORATOMS(m_obj->Neighbor, atm, atm_neighbor, j)
  {
    // check bond order
    if (m_obj->Bond[m_obj->Neighbor[j + 1]].order < 1)
      continue;

    // check atom filter
    if (atomIsExcluded(m_obj->AtomInfo[atm_neighbor]))
      continue;

    // check if closing a ring of size >= 3
    if (depth > 1 && atm_neighbor == m_indices[0]) {
      // found ring, add it to the selection
      onRingFound(m_obj, m_indices.data(), depth + 1);
    } else if (depth < m_indices.size() - 1) {
      // check for undesired ring with start != 0
      int i = depth;
      while ((--i) >= 0)
        if (atm_neighbor == m_indices[i])
          break; // stop recursion
      if (i == -1) {
        recursion(atm_neighbor, depth + 1);
      }
    }
  }
}

void AbstractRingFinder::apply(ObjectMolecule* obj, int atm)
{
  if (m_obj != obj) {
    m_obj = obj;
    ObjectMoleculeUpdateNeighbors(m_obj);
    prepareObject(m_obj);
  }

  recursion(atm, 0);
}
