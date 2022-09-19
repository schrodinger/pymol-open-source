#pragma once

#include "PyMOLGlobals.h"

struct DistSet;

namespace pymol
{

enum InteractionDir {
  cInteractionNone,
  cInteractionBoth,
  cInteractionForward,
};

/**
 * Halogen bond criteria
 */
class HalogenBondCriteria
{

public:
  // Cutoff distance
  float m_distance = 3.5f;

  // halogen as a donor minimum donor angle
  float m_as_donor_min_donor_angle = 140.0f;

  // halogen as a donor minimum acceptor angle
  float m_as_donor_min_acceptor_angle = 90.0f;

  // halogen as an acceptor minimum donor angle
  float m_as_acceptor_min_donor_angle = 120.0f;

  // halogen as an acceptor minimum acceptor angle
  float m_as_acceptor_min_acceptor_angle = 90.0f;

  // halogen as an acceptor maximum acceptor angle
  float m_as_acceptor_max_acceptor_angle = 170.0f;

public:
  /**
   * Construct a new Halogen Bond Criteria object
   * and initialize it's member variables
   */
  HalogenBondCriteria(PyMOLGlobals* G);
};

/**
 * Salt Bridge criteria
 */
class SaltBridgeCriteria
{

public:
  // Cutoff distance
  float m_distance = 5.0f;

public:
  /**
   * Construct a new Salt Bridge Criteria object and initialize it's member variables
   */
  SaltBridgeCriteria(PyMOLGlobals* G);
};

DistSet* FindPiInteractions(PyMOLGlobals* G,
    DistSet* ds,           //
    int sele1, int state1, //
    int sele2, int state2, //
    bool pipi = true,      //
    InteractionDir picat = cInteractionBoth);

/**
 * Find Halogen bond interactions
 *
 * @param ds - DistSet
 * @param sele1 - selections index
 * @param state1 - state index
 * @param sele2 - selection index
 * @param state2 - state index
 * @param cutoff - cutoff distance
 * @param result - result average distance
 *
 * @return DistSet - distance set
 */
DistSet* FindHalogenBondInteractions(PyMOLGlobals* G, DistSet* ds, int sele1,
    int state1, int sele2, int state2, float cutoff, float* result);

/**
 * Find Salt-bridge interactions
 *
 * @param ds - DistSet
 * @param sele1 - selections index
 * @param state1 - state index
 * @param sele2 - selection index
 * @param state2 - state index
 * @param cutoff - cutoff distance
 * @param result - result average distance
 *
 * @return DistSet - distance set
 */
DistSet* FindSaltBridgeInteractions(PyMOLGlobals* G, DistSet* ds, int sele1,
    int state1, int sele2, int state2, float cutoff, float* result);
}
