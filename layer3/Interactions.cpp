/**
 * Atomic interactions
 *
 * (c) 2020 Schrodinger, Inc.
 */

#include "Interactions.h"
#include "AtomInfo.h"
#include "AtomIterators.h"
#include "DistSet.h"
#include "Feedback.h"
#include "Map.h"
#include "ObjectMolecule.h"
#include "RingFinder.h"
#include "Selector.h"
#include "SelectorDef.h"
#include "AtomIterators.h"

#include <glm/vec3.hpp>
#include <glm/gtc/type_ptr.hpp>

#include <map>
#include <set>
#include <vector>

using AtomIndices = std::vector<int>;
using ObjRings = std::map<const ObjectMolecule*, std::set<AtomIndices>>;
using ObjAtoms = std::map<const ObjectMolecule*, AtomIndices>;
using Coords = std::vector<glm::vec3>;

namespace
{

/**
 * Ring finder which exposes founds rings in the `m_rings` member variable.
 */
class RingSetFinder : public AbstractRingFinder
{
  bool m_planar = false;

  /// The collected rings
  ObjRings m_rings;

  /**
   * @param planar If true, then only find planar rings
   */
  RingSetFinder(bool planar = false, int maxringsize = 7)
      : AbstractRingFinder(maxringsize)
      , m_planar(planar)
  {
  }

protected:
  void prepareObject(ObjectMolecule* obj) override
  {
    if (m_planar) {
      ObjectMoleculeVerifyChemistry(obj, cSelectorUpdateTableAllStates);
    }
  }

  bool atomIsExcluded(const AtomInfoType& atom) const override
  {
    return m_planar && atom.geom != cAtomInfoPlanar;
  }

  void onRingFound(
      ObjectMolecule* obj, const int* indices, size_t size) override
  {
    AtomIndices ring(indices, indices + size);
    std::sort(ring.begin(), ring.end());
    m_rings[obj].insert(std::move(ring));
  }

public:
  friend ObjRings FindRings(PyMOLGlobals*, int, bool);
};

/**
 * Find rings
 */
ObjRings FindRings(PyMOLGlobals* G, int sele, bool planar)
{
  RingSetFinder ringfinder(planar);

  for (SeleAtomIterator iter(G, sele); iter.next();) {
    ringfinder.apply(iter.obj, iter.getAtm());
  }

  return std::move(ringfinder.m_rings);
}

/**
 * Find cations
 */
ObjAtoms FindCations(PyMOLGlobals* G, int sele)
{
  ObjectMolecule* obj = nullptr;
  ObjAtoms cations;

  for (SeleAtomIterator iter(G, sele); iter.next();) {
    if (iter.obj != obj) {
      obj = iter.obj;
      ObjectMoleculeFixChemistry(obj, sele, sele, false);
      ObjectMoleculeVerifyChemistry(obj, cSelectorUpdateTableAllStates);
    }

    if (iter.getAtomInfo()->formalCharge > 0) {
      cations[obj].push_back(iter.getAtm());
    }
  }

  return cations;
}

/**
 * Simplified ring representation by center and normal.
 */
struct CNRing {
  glm::vec3 center{};
  glm::vec3 normal{};

  CNRing(const Coords& ringcoords)
  {
    for (const auto& xyz : ringcoords) {
      center += xyz;
    }

    if (!ringcoords.empty()) {
      center /= ringcoords.size();
    }

    if (ringcoords.size() >= 3) {
      auto v01 = ringcoords[1] - ringcoords[0];
      auto v12 = ringcoords[2] - ringcoords[1];
      normal = glm::normalize(glm::cross(v01, v12));
    }
  }
};

/**
 * Helper function to convert atom indices to simplified ring structures.
 */
std::vector<CNRing> cnrings_from_objrings(const ObjRings& objrings, int state)
{
  std::vector<CNRing> cnrings;

  for (auto& objitem : objrings) {
    const auto* obj = objitem.first;
    const auto* cs = obj->getCoordSet(state);

    if (!cs) {
      continue;
    }

    for (auto& ring : objitem.second) {
      Coords ringcoords;

      for (int atm : ring) {
        auto idx = cs->atmToIdx(atm);
        if (idx >= 0) {
          const float* v = cs->coordPtr(idx);
          ringcoords.emplace_back(v[0], v[1], v[2]);
        }
      }

      cnrings.emplace_back(ringcoords);
    }
  }

  return cnrings;
}

/**
 * Helper function to to convert atom indices to coordinates.
 */
Coords coords_from_objatoms(const ObjAtoms& objatoms, int state)
{
  Coords coords;

  for (auto& objitem : objatoms) {
    const auto* obj = objitem.first;
    const auto* cs = obj->getCoordSet(state);

    if (!cs) {
      continue;
    }

    for (int atm : objitem.second) {
      auto idx = cs->atmToIdx(atm);
      if (idx >= 0) {
        const float* v = cs->coordPtr(idx);
        coords.emplace_back(v[0], v[1], v[2]);
      }
    }
  }

  return coords;
}

/**
 * Helper function to convert to flat memory layout
 */
template <typename Range>
std::vector<float> flatten_ring_centers(const Range& rings)
{
  std::vector<float> flat;
  for (const auto& ring : rings) {
    flat.push_back(ring.center[0]);
    flat.push_back(ring.center[1]);
    flat.push_back(ring.center[2]);
  }
  return flat;
}

/**
 * Smaller angle between two vectors
 */
float angle_acute_degrees(glm::vec3 const& v1, glm::vec3 const& v2)
{
  return glm::degrees(std::acos(
      std::abs(glm::dot(v1, v2)) / (glm::length(v1) * glm::length(v2))));
}

void DistSetAddDistance(DistSet* ds,  //
    const float* v1, const float* v2, //
    int state1, int state2,           //
    AtomInfoType* ai1 = nullptr, AtomInfoType* ai2 = nullptr)
{
  auto G = ds->G;

  ds->MeasureInfo.emplace_front();
  auto& info = ds->MeasureInfo.front();
  info.offset = ds->NIndex;
  info.id[0] = ai1 ? AtomInfoCheckUniqueID(G, ai1) : 0;
  info.id[1] = ai2 ? AtomInfoCheckUniqueID(G, ai2) : 0;
  info.state[0] = state1;
  info.state[1] = state2;
  info.measureType = cRepDash;

  ds->Coord.reserve((ds->NIndex + 2) * 3);
  for (size_t i = 0; i != 3; ++i) {
    ds->Coord[(ds->NIndex + 0) * 3 + i] = v1[i];
    ds->Coord[(ds->NIndex + 1) * 3 + i] = v2[i];
  }

  ds->NIndex += 2;
}

} // namespace

namespace pymol
{

/**
 * Find pi-pi and/or pi-cation interactions.
 *
 * @param ds Measurement state to update, or nullptr to create a new one.
 * @param sele1 Selection index
 * @param state1 Object state
 * @param sele2 Selection index
 * @param state2 Object state
 * @param pipi If true, then search for pi-pi interactions
 * @param picat If "both" then search for pi-cation interactions in both
 * directions. If "forward" then search for rings in selection 1 and for cations
 * in selection 2.
 */
DistSet* FindPiInteractions(PyMOLGlobals* G,
    DistSet* ds,           //
    int sele1, int state1, //
    int sele2, int state2, //
    bool pipi,             //
    InteractionDir picat)
{
  // These constants are borrowed from mmshare/include/structureinteraction.h
  constexpr auto RING_ALIGNMENT_MAX_ANGLE = 40.0;
  constexpr auto DEFAULT_PI_CATION_MAXIMUM_DISTANCE = 6.6;
  constexpr auto DEFAULT_PI_CATION_MAXIMUM_ANGLE = 30.0;
  constexpr auto PIPI_FACE_TO_FACE_MAXIMUM_DISTANCE = 4.4;
  constexpr auto PIPI_FACE_TO_FACE_MAXIMUM_ANGLE = 30.0;
  constexpr auto PIPI_EDGE_TO_FACE_MAXIMUM_DISTANCE = 5.5;
  constexpr auto PIPI_EDGE_TO_FACE_MINIMUM_ANGLE = 60.0;

  const bool sele_is_same = sele1 == sele2 && state1 == state2;

  auto rings1 = cnrings_from_objrings(FindRings(G, sele1, true), state1);

  // data structure for fast neighbor lookup
  std::unique_ptr<MapType> centers1map(
      MapNew(G, -DEFAULT_PI_CATION_MAXIMUM_DISTANCE,
          flatten_ring_centers(rings1).data(), rings1.size(), nullptr));

  if (!ds) {
    ds = DistSetNew(G);
  }

  if (pipi) {
    auto rings2 =
        sele_is_same ? rings1
                     : cnrings_from_objrings(FindRings(G, sele2, true), state2);

    int i = -1;
    for (const auto& ring2 : rings2) {
      ++i;
      for (int j : MapEIter(*centers1map, glm::value_ptr(ring2.center))) {
        if (sele_is_same && j >= i) {
          continue;
        }

        auto const& ring1 = rings1[j];
        auto const v = ring2.center - ring1.center;
        auto const distance = glm::length(v);

        if (distance < 1e-2 || //
            distance > PIPI_EDGE_TO_FACE_MAXIMUM_DISTANCE) {
          continue;
        }

        auto const normal_to_v_angle_i = angle_acute_degrees(ring2.normal, v);
        auto const normal_to_v_angle_j = angle_acute_degrees(ring1.normal, v);

        if (normal_to_v_angle_i > RING_ALIGNMENT_MAX_ANGLE &&
            normal_to_v_angle_j > RING_ALIGNMENT_MAX_ANGLE) {
          // collinear
          continue;
        }

        auto const angle = angle_acute_degrees(ring2.normal, ring1.normal);

        if (angle < PIPI_FACE_TO_FACE_MAXIMUM_ANGLE &&
            distance < PIPI_FACE_TO_FACE_MAXIMUM_DISTANCE) {
          PRINTFB(G, FB_DistSet, FB_Blather)
            "face-to-face %d %d\n", i, j ENDFB(G);
        } else if (angle > PIPI_EDGE_TO_FACE_MINIMUM_ANGLE) {
          PRINTFB(G, FB_DistSet, FB_Blather)
            "edge-to-face %d %d\n", i, j ENDFB(G);
        } else {
          continue;
        }

        DistSetAddDistance(ds,
            glm::value_ptr(ring1.center), //
            glm::value_ptr(ring2.center), state1, state2);
      }
    }
  }

  if (picat != cInteractionNone) {
    auto cations2 = coords_from_objatoms(FindCations(G, sele2), state2);

    int i = -1;
    for (const auto& cation2 : cations2) {
      ++i;
      for (int j : MapEIter(*centers1map, glm::value_ptr(cation2))) {
        auto const v = cation2 - rings1[j].center;
        auto const distance = glm::length(v);

        if (distance > DEFAULT_PI_CATION_MAXIMUM_DISTANCE) {
          continue;
        }

        if (angle_acute_degrees(rings1[j].normal, v) >
            DEFAULT_PI_CATION_MAXIMUM_ANGLE) {
          // collinear
          continue;
        }

        DistSetAddDistance(ds,
            glm::value_ptr(rings1[j].center), //
            glm::value_ptr(cation2), state1, state2);

        PRINTFB(G, FB_DistSet, FB_Blather)
          "pi-cat %d %d\n", i, j ENDFB(G);
      }
    }

    // vice-versa
    if (!sele_is_same && picat == cInteractionBoth) {
      FindPiInteractions(G, ds, sele2, state2, sele1, state1, cInteractionNone,
          cInteractionForward);
    }
  }

  return ds;
}

/**
 * Test halogen-bond criteria when halogen is an acceptor.
 * Halogen acceptor bonds have relations between the donor H atom (H) and the
 * atom bonded to it (D), the halogen acceptor (X) and the atom bonded to it
 * (B), D–H...X–B. The following must be true for a valid acceptor halogen bond:
 *
 * The H...X distance must be less than a specified maximum distance.
 * The D–H...X angle must be greater than a specified minimum value.
 * The H...X–B angle must be greater than a specified minimum value.
 * The H...X–B angle must be less than a specified maximum value.
 *
 * @param v_x_h - X-H vector
 * @param v_d_h - D-H vector
 * @param v_x_b - X-B vector
 * @param hbc - halogen-bond criteria
 *
 * @return true - if it meets a halogen-bond interaction criteria
 * @return false - if it does not meet a halogen-bond interaction criteria
 */
static bool TestHalogenBondAcceptor(
    float* v_x_h, float* v_d_h, float* v_x_b, HalogenBondCriteria* hbc)
{
  float n_v_x_h[3];
  normalize23f(v_x_h, n_v_x_h);

  float n_v_d_h[3];
  normalize23f(v_d_h, n_v_d_h);

  float n_v_x_b[3];
  normalize23f(v_x_b, n_v_x_b);

  float dp = dot_product3f(n_v_x_h, n_v_d_h);
  float angle = 180.0 * acos(dp) / PI;
  if (angle < hbc->m_as_acceptor_min_donor_angle) {
    return false;
  }

  dp = dot_product3f(n_v_x_h, n_v_x_b);
  angle = 180.0 * acos(dp) / PI;
  if (angle < hbc->m_as_acceptor_min_acceptor_angle ||
      angle > hbc->m_as_acceptor_max_acceptor_angle) {
    return false;
  }

  float dist = length3f(v_x_h);
  return dist <= hbc->m_distance;
}

/**
 * Test halogen-bond criteria when halogen is a donor.
 * Halogen bonds in which the halogen acts as donor are defined in a similar way
 * to hydrogen bonds, by relations between four atoms: the donor halogen atom
 * (X), the donor atom (D) bonded to it, the acceptor atom (A), and another
 * neighbor atom (B) bonded to A, represented as D–X...A–B. The following must
 * be true for a valid donor halogen bond:
 *
 * The X...A distance must be less than a specified maximum distance.
 * The D–X...A angle must be greater than a specified minimum value.
 * The X...A–B angle must be greater than a specified minimum value.
 *
 * @param v_a_x - A-X vector
 * @param v_d_x - D-X vector
 * @param v_a_b - A-B vector
 * @param hbc - halogen-bond criteria
 *
 * @return true - if it meets a halogen-bond interaction criteria
 * @return false - if it does not meet a halogen-bond interaction criteria
 */
static bool TestHalogenBondDonor(
    float* v_a_x, float* v_d_x, float* v_a_b, HalogenBondCriteria* hbc)
{
  float n_v_a_x[3];
  normalize23f(v_a_x, n_v_a_x);

  float n_v_d_x[3];
  normalize23f(v_d_x, n_v_d_x);

  float n_v_a_b[3];
  normalize23f(v_a_b, n_v_a_b);

  float angle = dot_product3f(n_v_a_x, n_v_d_x);
  angle = 180.0 * acos(angle) / PI;
  if (angle < hbc->m_as_donor_min_donor_angle) {
    return false;
  }

  angle = dot_product3f(n_v_a_x, n_v_a_b);
  angle = 180.0 * acos(angle) / PI;
  if (angle < hbc->m_as_donor_min_acceptor_angle) {
    return false;
  }

  float dist = length3f(v_a_x);
  return dist <= hbc->m_distance;
}

/**
 * Check halogen-bond as a donor
 *
 * @param don_obj - donor object molecule
 * @param don_atom - donor atom id
 * @param don_state - donor state index
 * @param acc_obj - acceptor object molecule
 * @param acc_atom - acceptor atom id
 * @param acc_state - acceptor state index
 * @param hbc - halogen-bond criteria
 *
 * @return true - if it meets a halogen-bond interaction criteria
 * @return false - if it does not meet a halogen-bond interaction criteria
 */
static bool CheckHalogenBondAsDonor(ObjectMolecule* don_obj, int don_atom,
    int don_state, ObjectMolecule* acc_obj, int acc_atom, int acc_state,
    HalogenBondCriteria* hbc)
{
  bool result = false;
  const CoordSet* csD = nullptr;
  const CoordSet* csA = nullptr;

  // first, check for existence of coordinate sets
  if (don_state >= 0 && don_state < don_obj->NCSet) {
    csD = don_obj->CSet[don_state];
  }

  if (acc_state >= 0 && acc_state < acc_obj->NCSet) {
    csA = acc_obj->CSet[acc_state];
  }

  if (csD == nullptr) {
    return false;
  }

  if (csA == nullptr) {
    return false;
  }

  if (don_atom >= don_obj->NAtom) {
    return false;
  }

  if (acc_atom >= acc_obj->NAtom) {
    return false;
  }

  AtomInfoType* don_atom_info = don_obj->AtomInfo + don_atom;

  if (don_atom_info->protons == cAN_Cl || don_atom_info->protons == cAN_Br ||
      don_atom_info->protons == cAN_I) {

    // now check for coordinates of these actual atoms
    auto idxD = csD->atmToIdx(don_atom);
    auto idxA = csA->atmToIdx(acc_atom);

    if (idxA >= 0 && idxD >= 0) {

      const float* vDon = csD->coordPtr(idxD);
      const float* vAcc = csA->coordPtr(idxA);

      float v_d[3];
      if (ObjectMoleculeGetNeighborVector(don_obj, don_atom, don_state, v_d)) {

        float v_b[3];
        if (ObjectMoleculeGetNeighborVector(
                acc_obj, acc_atom, acc_state, v_b)) {

          float v_d_x[3];
          subtract3f(v_d, vDon, v_d_x);

          float v_a_b[3];
          subtract3f(vAcc, v_b, v_a_b);

          float v_a_x[3];
          subtract3f(vAcc, vDon, v_a_x);

          int result = TestHalogenBondDonor(v_a_x, v_d_x, v_a_b, hbc);
          if (result) {
            return true;
          }
        }
      }
    }
  }
  return result;
}

/**
 * Check halogen-bond as acceptor
 *
 * @param don_obj - donor object molecule
 * @param don_atom - donor atom id
 * @param don_state - donor state index
 * @param acc_obj - acceptor object molecule
 * @param acc_atom - acceptor atom id
 * @param acc_state - acceptor state index
 * @param hbc - halogen-bond criteria
 *
 * @return true - if it meets a halogen-bond interaction criteria
 * @return false - if it does not meet a halogen-bond interaction criteria
 */
static bool CheckHalogenBondAsAcceptor(ObjectMolecule* don_obj, int don_atom,
    int don_state, ObjectMolecule* acc_obj, int acc_atom, int acc_state,
    HalogenBondCriteria* hbc)
{
  bool result = false;
  const CoordSet* csD = nullptr;
  const CoordSet* csA = nullptr;

  // first, check for existence of coordinate sets
  if (don_state >= 0 && don_state < don_obj->NCSet) {
    csD = don_obj->CSet[don_state];
  }

  if (acc_state >= 0 && acc_state < acc_obj->NCSet) {
    csA = acc_obj->CSet[acc_state];
  }

  if (csD == nullptr) {
    return false;
  }

  if (csA == nullptr) {
    return false;
  }

  if (don_atom >= don_obj->NAtom) {
    return false;
  }

  if (acc_atom >= acc_obj->NAtom) {
    return false;
  }

  AtomInfoType* don_atom_info = don_obj->AtomInfo + don_atom;
  AtomInfoType* acc_atom_info = acc_obj->AtomInfo + acc_atom;

  if ((don_atom_info->protons == cAN_H || don_atom_info->protons == cAN_N ||
          don_atom_info->protons == cAN_O || don_atom_info->protons == cAN_S) &&
      (acc_atom_info->protons == cAN_Cl || acc_atom_info->protons == cAN_Br ||
          acc_atom_info->protons == cAN_I)) {

    auto idxD = csD->atmToIdx(don_atom);
    auto idxA = csA->atmToIdx(acc_atom);

    if ((idxA >= 0) && (idxD >= 0)) {
      const float* vDon = csD->coordPtr(idxD);
      const float* vAcc = csA->coordPtr(idxA);

      float v_d[3];
      if (ObjectMoleculeGetNeighborVector(don_obj, don_atom, don_state, v_d)) {

        float v_b[3];
        if (ObjectMoleculeGetNeighborVector(
                acc_obj, acc_atom, acc_state, v_b)) {
          float v_d_h[3];
          subtract3f(v_d, vDon, v_d_h);

          float v_x_b[3];
          subtract3f(vAcc, v_b, v_x_b);

          float v_x_h[3];
          subtract3f(vAcc, vDon, v_x_h);

          bool result = TestHalogenBondAcceptor(v_x_h, v_d_h, v_x_b, hbc);
          if (result) {
            return true;
          }
        }
      }
    }
  }

  return result;
}

SaltBridgeCriteria::SaltBridgeCriteria(PyMOLGlobals* G)
{
  m_distance =
      SettingGet<float>(G, nullptr, nullptr, cSetting_salt_bridge_distance);
}

HalogenBondCriteria::HalogenBondCriteria(PyMOLGlobals* G)
{
  m_distance =
      SettingGet<float>(G, nullptr, nullptr, cSetting_halogen_bond_distance);
  m_as_donor_min_donor_angle = SettingGet<float>(
      G, nullptr, nullptr, cSetting_halogen_bond_as_donor_min_donor_angle);
  m_as_donor_min_acceptor_angle = SettingGet<float>(
      G, nullptr, nullptr, cSetting_halogen_bond_as_donor_min_acceptor_angle);
  m_as_acceptor_min_donor_angle = SettingGet<float>(
      G, nullptr, nullptr, cSetting_halogen_bond_as_acceptor_min_donor_angle);
  m_as_acceptor_min_acceptor_angle = SettingGet<float>(G, nullptr, nullptr,
      cSetting_halogen_bond_as_acceptor_min_acceptor_angle);
  m_as_acceptor_max_acceptor_angle = SettingGet<float>(G, nullptr, nullptr,
      cSetting_halogen_bond_as_acceptor_max_acceptor_angle);
}

/**
 * Prepare neighbor tables
 *
 * @param sele1 - selections index
 * @param state1 - state index
 * @param sele2 - selection index
 * @param state2 - state index
 *
 * @return maximum number of atoms
 */
int PrepareNeighborTables(
    PyMOLGlobals* G, int sele1, int state1, int sele2, int state2)
{
  CSelector* I = G->Selector;

  // update states: if the two are the same, update that one state, else update
  // all states
  if (state1 < 0 || state2 < 0 || state1 != state2) {
    SelectorUpdateTable(G, cSelectorUpdateTableAllStates, -1);
  } else {
    SelectorUpdateTable(G, state1, -1);
  }

  // find and prepare (neighbortables) in any participating Molecular objects
  // fill in all the neighbor tables
  int max_n_atom = I->Table.size();
  ObjectMolecule* lastObj = nullptr;
  for (int a = cNDummyAtoms; a < I->Table.size(); a++) {
    // foreach atom in the session, get its identifier and ObjectMolecule to
    // which it belongs
    int at = I->Table[a].atom; // grab the atom ID from the Selectors->Table
    ObjectMolecule* obj =
        I->Obj[I->Table[a].model]; // -- JV -- quick way to get an object from
                                   // an atom
    int s = obj->AtomInfo[at]
                .selEntry; // grab the selection entry# from this Atoms Info
    if (obj != lastObj) {
      if (max_n_atom < obj->NAtom) {
        max_n_atom = obj->NAtom;
      }
      // if the current atom is in sele1 or sele2 then update it's object's
      // neighbor table
      if (SelectorIsMember(G, s, sele1) || SelectorIsMember(G, s, sele2)) {
        // if hbonds (so, more than just distance)
        ObjectMoleculeVerifyChemistry(obj, -1);
        lastObj = obj;
      }
    }
  }

  return max_n_atom;
}

/**
 * Insert distance info into DistSet and copy atom's coordinates
 *
 * @param ds - distance set
 * @param state1 - state index
 * @param state2 - state index
 * @param ai1 - atom info
 * @param ai2 - atom info
 * @param atom1_vv - atom coordinates
 * @param atom2_vv - atom coordinates
 * @param numVerts - number of vertices
 */
void InsertDistanceInfo(PyMOLGlobals* G, DistSet* ds, int state1, int state2,
    AtomInfoType* ai1, AtomInfoType* ai2, float* atom1_vv, float* atom2_vv,
    int numVerts)
{
  // Insert DistInfo records for updating distances
  // Init/Add the elem to the DistInfo list
  ds->MeasureInfo.emplace_front();
  auto* atom1Info = &ds->MeasureInfo.front();

  // TH
  atom1Info->id[0] = AtomInfoCheckUniqueID(G, ai1);
  atom1Info->id[1] = AtomInfoCheckUniqueID(G, ai2);

  atom1Info->offset = numVerts; // offset into this DSet's Coord
  atom1Info->state[0] = state1; // state1 of sel1
  atom1Info->state[1] = state2;
  atom1Info->measureType = cRepDash; // DISTANCE-dash

  auto& coords = ds->Coord;

  // see if coords has room at another 6 floats
  VLACheck(coords, float, (numVerts * 3) + 6);
  float* vv0 = coords + (numVerts * 3);

  if (atom1_vv != nullptr && atom2_vv != nullptr) {
    const size_t count_to_copy = 3;

    std::copy_n(atom1_vv, count_to_copy, vv0);
    vv0 += count_to_copy;
    atom1_vv += count_to_copy;

    std::copy_n(atom2_vv, count_to_copy, vv0);
    vv0 += count_to_copy;
    atom2_vv += count_to_copy;
  }
}

/**
 * Create coverage vector
 *
 * @param sele1 - selection index
 * @param sele2 - selection index
 * @return vector of booleans which determines if a given atom appears in sele1
 * and sele2
 */
static std::vector<bool> CreateCoverage(PyMOLGlobals* G, int sele1, int sele2)
{
  CSelector* I = G->Selector;
  std::vector<bool> result(I->Table.size());

  // coverage determines if a given atom appears in sele1 and sele2
  for (SelectorAtomIterator iter(I); iter.next();) {
    int s = iter.getAtomInfo()->selEntry;
    if (SelectorIsMember(G, s, sele1) && SelectorIsMember(G, s, sele2)) {
      result[iter.a] = true;
    }
  }

  return result;
}

DistSet* FindHalogenBondInteractions(PyMOLGlobals* G, DistSet* ds, int sele1,
    int state1, int sele2, int state2, float cutoff, float* result)
{
  CSelector* I = G->Selector;

  int numVerts = 0;
  *result = 0.0f;
  // if the dist set exists, get info from it, otherwise get a new one
  if (ds == nullptr) {
    ds = DistSetNew(G);
  } else {
    numVerts = ds->NIndex; // number of vertices
  }

  auto& coords = ds->Coord;
  coords.reserve(10);

  int max_n_atom = PrepareNeighborTables(G, sele1, state1, sele2, state2);

  HalogenBondCriteria halogenbcRec(G);
  cutoff = halogenbcRec.m_distance;
  if (cutoff < 0.0f) {
    const float max_cutoff = 1000.0f;
    cutoff = max_cutoff;
  }

  // coverage determines if a given atom appears in sel1 and sel2
  std::vector<bool> coverage = CreateCoverage(G, sele1, sele2);

  // this creates an interleaved list of ints for mapping ids to states within a
  // given neighborhood
  std::vector<int> interstate_vector =
      SelectorGetInterstateVector(G, sele1, state1, sele2, state2, cutoff);
  int cnt = interstate_vector.size() / 2;

  std::vector<int> zero(max_n_atom);
  std::vector<int> scratch(max_n_atom);

  float dist_sum = 0.0f;
  int dist_cnt = 0;

  // for each state
  for (int a = 0; a < cnt; a++) {

    // get the interstate atom identifier for the two atoms to distance
    int a1 = interstate_vector[a * 2];
    int a2 = interstate_vector[a * 2 + 1];

    // check their coverage to avoid duplicates
    if (a1 < a2 || (a1 != a2 && !(coverage[a1] && coverage[a2])) ||
        (state1 != state2)) {
      // eliminate reverse duplicates

      // get the object-local atom ID
      int at1 = I->Table[a1].atom;
      int at2 = I->Table[a2].atom;

      if (sele1 == sele2 && at1 > at2) {
        continue;
      }

      // get the object for this global atom ID
      ObjectMolecule* obj1 = I->Obj[I->Table[a1].model];
      ObjectMolecule* obj2 = I->Obj[I->Table[a2].model];

      // the states are valid for these two atoms
      if (state1 < obj1->NCSet && state2 < obj2->NCSet) {
        // get the coordinate sets for both atoms
        CoordSet* cs1 = obj1->CSet[state1];
        CoordSet* cs2 = obj2->CSet[state2];
        if (cs1 != nullptr && cs2 != nullptr) {
          // for bonding
          float* don_vv = nullptr;
          float* acc_vv = nullptr;

          // grab the appropriate atom information for this object-local atom
          AtomInfoType* ai1 = obj1->AtomInfo + at1;
          AtomInfoType* ai2 = obj2->AtomInfo + at2;

          int idx1 = cs1->atmToIdx(at1);
          int idx2 = cs2->atmToIdx(at2);

          if (idx1 >= 0 && idx2 >= 0) {
            // actual distance calculation from ptA to ptB
            float dist =
                (float) diff3f(cs1->coordPtr(idx1), cs2->coordPtr(idx2));

            // if we pass the bonding cutoff
            if (dist < cutoff) {

              bool a_keeper = false;

              if (ai1->hb_donor) {
                a_keeper = CheckHalogenBondAsAcceptor(
                    obj1, at1, state1, obj2, at2, state2, &halogenbcRec);
                if (a_keeper) {
                  don_vv = cs1->coordPtr(idx1);
                  acc_vv = cs2->coordPtr(idx2);
                }
              } else if (ai2->hb_donor) {
                a_keeper = CheckHalogenBondAsAcceptor(
                    obj2, at2, state2, obj1, at1, state1, &halogenbcRec);
                if (a_keeper) {
                  don_vv = cs2->coordPtr(idx2);
                  acc_vv = cs1->coordPtr(idx1);
                }
              }

              if (a_keeper == false) {
                if (ai2->hb_acceptor) {
                  a_keeper = CheckHalogenBondAsDonor(
                      obj1, at1, state1, obj2, at2, state2, &halogenbcRec);
                  if (a_keeper) {
                    don_vv = cs1->coordPtr(idx1);
                    acc_vv = cs2->coordPtr(idx2);
                  }
                } else if (ai1->hb_acceptor) {
                  a_keeper = CheckHalogenBondAsDonor(
                      obj2, at2, state2, obj1, at1, state1, &halogenbcRec);
                  if (a_keeper) {
                    don_vv = cs2->coordPtr(idx2);
                    acc_vv = cs1->coordPtr(idx1);
                  }
                }
              }

              if (a_keeper) {
                InsertDistanceInfo(
                    G, ds, state1, state2, ai1, ai2, don_vv, acc_vv, numVerts);
                dist_cnt++;
                dist_sum += dist;
                numVerts += 2;
              }
            }
          }
        }
      }
    }
  }

  if (dist_cnt > 0) {
    (*result) = dist_sum / dist_cnt;
  }

  if (coords) {
    coords.resize((numVerts + 1) * 3);
  }

  ds->NIndex = numVerts;

  return ds;
}

DistSet* FindSaltBridgeInteractions(PyMOLGlobals* G, DistSet* ds, int sele1,
    int state1, int sele2, int state2, float cutoff, float* result)
{
  CSelector* I = G->Selector;

  int numVerts = 0;
  *result = 0.0f;
  // if the dist set exists, get info from it, otherwise get a new one
  if (ds == nullptr) {
    ds = DistSetNew(G);
  } else {
    numVerts = ds->NIndex; // number of vertices
  }

  auto& coords = ds->Coord;
  coords.reserve(10);

  int max_n_atom = PrepareNeighborTables(G, sele1, state1, sele2, state2);

  SaltBridgeCriteria saltbcRec(G);
  cutoff = saltbcRec.m_distance;
  if (cutoff < 0.0f) {
    const float max_cutoff = 1000.0f;
    cutoff = max_cutoff;
  }

  // coverage determines if a given atom appears in sel1 and sel2
  std::vector<bool> coverage = CreateCoverage(G, sele1, sele2);

  // this creates an interleaved list of ints for mapping ids to states within a
  // given neighborhood
  std::vector<int> interstate_vector =
      SelectorGetInterstateVector(G, sele1, state1, sele2, state2, cutoff);
  int cnt = interstate_vector.size() / 2;

  std::vector<int> zero(max_n_atom);
  std::vector<int> scratch(max_n_atom);

  float dist_sum = 0.0f;
  int dist_cnt = 0;

  // for each state
  for (int a = 0; a < cnt; a++) {

    // get the interstate atom identifier for the two atoms to distance
    int a1 = interstate_vector[a * 2];
    int a2 = interstate_vector[a * 2 + 1];

    // check their coverage to avoid duplicates
    if (a1 < a2 || (a1 != a2 && !(coverage[a1] && coverage[a2])) ||
        (state1 != state2)) {
      // eliminate reverse duplicates

      // get the object-local atom ID
      int at1 = I->Table[a1].atom;
      int at2 = I->Table[a2].atom;

      if (sele1 == sele2 && at1 > at2) {
        continue;
      }

      // get the object for this global atom ID
      ObjectMolecule* obj1 = I->Obj[I->Table[a1].model];
      ObjectMolecule* obj2 = I->Obj[I->Table[a2].model];

      // the states are valid for these two atoms
      if (state1 < obj1->NCSet && state2 < obj2->NCSet) {
        // get the coordinate sets for both atoms
        CoordSet* cs1 = obj1->CSet[state1];
        CoordSet* cs2 = obj2->CSet[state2];
        if (cs1 != nullptr && cs2 != nullptr) {
          // for bonding
          float* anion_vv = nullptr;
          float* cation_vv = nullptr;

          // grab the appropriate atom information for this object-local atom
          AtomInfoType* ai1 = obj1->AtomInfo + at1;
          AtomInfoType* ai2 = obj2->AtomInfo + at2;

          int atom1_charge = ai1->formalCharge;
          int atom2_charge = ai2->formalCharge;
          if (atom1_charge * atom2_charge >= 0) {
            continue;
          }

          if (ai1->isHydrogen() || ai2->isHydrogen()) {
            continue;
          }

          int idx1 = cs1->atmToIdx(at1);
          int idx2 = cs2->atmToIdx(at2);

          if (idx1 >= 0 && idx2 >= 0) {
            // actual distance calculation from ptA to ptB
            float dist =
                (float) diff3f(cs1->coordPtr(idx1), cs2->coordPtr(idx2));

            // if we pass the bonding cutoff
            if (dist < cutoff) {

              if (atom1_charge < 0) {
                anion_vv = cs1->coordPtr(idx1);
                cation_vv = cs2->coordPtr(idx2);
              } else {
                anion_vv = cs2->coordPtr(idx2);
                cation_vv = cs1->coordPtr(idx1);
              }

              InsertDistanceInfo(G, ds, state1, state2, ai1, ai2, anion_vv,
                  cation_vv, numVerts);
              dist_cnt++;
              dist_sum += dist;
              numVerts += 2;
            }
          }
        }
      }
    }
  }

  if (dist_cnt > 0) {
    (*result) = dist_sum / dist_cnt;
  }

  if (coords) {
    coords.resize((numVerts + 1) * 3);
  }

  ds->NIndex = numVerts;

  return ds;
}

} // namespace pymol
