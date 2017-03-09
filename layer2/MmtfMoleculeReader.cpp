/*
 * MMTF Molecule Reader
 *
 * (c) 2016 Schrodinger, Inc.
 */

#include "os_std.h"

#include "ObjectMolecule.h"
#include "Feedback.h"

#ifndef _PYMOL_NO_MSGPACKC

#include <mmtf_parser.h>

#include <algorithm>

#include "AssemblyHelpers.h"
#include "AtomInfo.h"
#include "Err.h"
#include "Lex.h"
#include "MemoryDebug.h"
#include "Rep.h"

const int ss_map[] = {
    'H', // pi helix
    'L', // bend
    'H', // alpha helix
    'S', // extended
    'H', // 3-10 helix
    'S', // bridge
    'L', // turn
    'L', // coil
    0
};

/*
 * Get the assembly coordinate sets
 *
 * See also: read_pdbx_struct_assembly (CifMoleculeReader.cpp)
 */
static
CoordSet ** get_assembly_csets(PyMOLGlobals * G,
    const MMTF_container * container,
    const AtomInfoType * atInfo,
    const CoordSet * cset)
{
  const char * assembly_id = SettingGetGlobal_s(G, cSetting_assembly);

  if (!assembly_id || !assembly_id[0])
    return NULL;

  const MMTF_BioAssembly * assembly = NULL;

  // look up assembly by name
  for (auto a = container->bioAssemblyList,
      a_end = a + container->bioAssemblyListCount;
      a != a_end; ++a) {
    if (strcmp(a->name, assembly_id) == 0) {
      assembly = a;
      break;
    }
  }

  if (!assembly) {
    PRINTFB(G, FB_Executive, FB_Details)
      " ExecutiveLoad-Detail: No such assembly: '%s'\n", assembly_id ENDFB(G);
    return NULL;
  }

  PRINTFB(G, FB_Executive, FB_Details)
    " ExecutiveLoad-Detail: Creating assembly '%s'\n", assembly_id ENDFB(G);

  int ncsets = assembly->transformListCount;
  CoordSet ** csets = VLACalloc(CoordSet *, ncsets);

  for (int state = 0; state < ncsets; ++state) {
    auto trans = assembly->transformList + state;

    // get set of chains for this transformation
    std::set<lexidx_t> chains_set;
    for (auto ci_it = trans->chainIndexList,
        ci_it_end = ci_it + trans->chainIndexListCount;
        ci_it != ci_it_end; ++ci_it) {
      const char * chain = container->chainIdList[*ci_it];
      auto result = OVLexicon_BorrowFromCString(G->Lexicon, chain);
      if (OVreturn_IS_OK(result)) {
        chains_set.insert(result.word);
      }
    }

    // copy and transform
    csets[state] = CoordSetCopyFilterChains(cset, atInfo, chains_set);
    CoordSetTransform44f(csets[state], trans->matrix);
  }

  return csets;
}
#endif

ObjectMolecule * ObjectMoleculeReadMmtfStr(PyMOLGlobals * G, ObjectMolecule * I,
    const char * st, int st_len, int frame, int discrete, int quiet, int multiplex, int zoom)
{
#ifdef _PYMOL_NO_MSGPACKC
  PRINTFB(G, FB_ObjectMolecule, FB_Errors)
    " Error: This build has no fast MMTF support.\n" ENDFB(G);
  return NULL;
#else

  if (I) {
    PRINTFB(G, FB_ObjectMolecule, FB_Errors)
      " Error: loading MMTF into existing object not supported, please use 'create'\n"
      "        to append to an existing object.\n" ENDFB(G);
    return NULL;
  }

  if (multiplex > 0) {
    PRINTFB(G, FB_ObjectMolecule, FB_Errors)
      " Error: multiplex not supported for MMTF\n" ENDFB(G);
    return NULL;
  }

  MMTF_container * container = MMTF_container_new();

  if (!MMTF_unpack_from_string(st, st_len, container)) {
    PRINTFB(G, FB_ObjectMolecule, FB_Errors)
      " Error: Failed to load MMTF file\n" ENDFB(G);
    MMTF_container_free(container);
    return NULL;
  }

  PRINTFB(G, FB_ObjectMolecule, FB_Details)
    " MMTF structureId: '%s', mmtfVersion: '%s'\n",
    (container->structureId ? container->structureId : ""),
    (container->mmtfVersion ? container->mmtfVersion : "") ENDFB(G);

  int chainIndex = 0;
  int groupIndex = 0;
  int atomIndex = 0;

  // template atom
  AtomInfoType tai;
  memset(&tai, 0, sizeof(AtomInfoType));
  tai.visRep = RepGetAutoShowMask(G);

  I = ObjectMoleculeNew(G, /* discrete */ 1);
  I->Obj.Color = AtomInfoUpdateAutoColor(G);
  I->NAtom = container->numAtoms;
  I->NCSet = container->numModels;
  I->Bond = VLAlloc(BondType, container->numBonds);
  VLASize(I->AtomInfo, AtomInfoType, I->NAtom);
  VLASize(I->CSet, CoordSet *, I->NCSet);

  int nindexEstimate = container->numAtoms / std::max(1, container->numModels);
  bool use_auth = SettingGetGlobal_b(G, cSetting_cif_use_auth);

  // symmetry
  if (container->unitCell &&
      container->spaceGroup &&
      container->spaceGroup[0]) {
    CSymmetry * symmetry = I->Symmetry = SymmetryNew(G);

    for (int i = 0; i < 3; i++) {
      symmetry->Crystal->Dim[i] = container->unitCell[i];
      symmetry->Crystal->Angle[i] = container->unitCell[i + 3];
    }

    strncpy(symmetry->SpaceGroup, container->spaceGroup, WordLength - 1);

    SymmetryUpdate(symmetry);
    CrystalUpdate(symmetry->Crystal);
  }

  // models (states)
  for (int modelIndex = 0; modelIndex < container->numModels; ++modelIndex) {
    int modelChainCount = container->chainsPerModel[modelIndex];

    CoordSet * cset = CoordSetNew(G);
    cset->Coord = VLAlloc(float, 3 * nindexEstimate);
    cset->IdxToAtm = VLAlloc(int, nindexEstimate);
    cset->Obj = I;
    I->CSet[modelIndex] = cset;

    // chains
    for (int j = 0; j < modelChainCount; ++j, ++chainIndex) {
      if (container->chainNameList)
        LexAssign(G, tai.chain, container->chainNameList[chainIndex]);
      if (use_auth)
        LexAssign(G, tai.segi, container->chainIdList[chainIndex]);

      int chainGroupCount = container->groupsPerChain[chainIndex];

      // groups (residues)
      for (int k = 0; k < chainGroupCount; ++k, ++groupIndex) {
        const MMTF_GroupType * group =
          container->groupList +
          container->groupTypeList[groupIndex];

        if (container->secStructList)
          tai.ssType[0] = ss_map[container->secStructList[groupIndex]];

        LexAssign(G, tai.resn, group->groupName);
        tai.hetatm = group->singleLetterCode == '?';
        tai.flags = tai.hetatm ? cAtomFlag_ignore : 0;

        if (use_auth) {
          tai.resv = container->groupIdList[groupIndex];
          if (container->insCodeList)
            tai.setInscode(container->insCodeList[groupIndex]);
        } else if (container->sequenceIndexList) {
          tai.resv = container->sequenceIndexList[groupIndex];
        }

        // bonds
        for (int l = 0, offset = atomIndex; l < group->bondAtomListCount / 2; ++l) {
          BondTypeInit2(I->Bond + (I->NBond++),
              offset + group->bondAtomList[l * 2],
              offset + group->bondAtomList[l * 2 + 1],
              group->bondOrderListCount ? group->bondOrderList[l] : 1);
        }

        int groupAtomCount = group->atomNameListCount;

        // atoms
        for (int l = 0; l < groupAtomCount; ++l, ++atomIndex) {
          AtomInfoType * ai = I->AtomInfo + atomIndex;
          AtomInfoCopy(G, &tai, ai);

          ai->rank = atomIndex;
          ai->formalCharge = group->formalChargeList[l];
          LexAssign(G, ai->name, group->atomNameList[l]);
          strncpy(ai->elem, group->elementList[l], cElemNameLen);

          if (container->atomIdList)
            ai->id = container->atomIdList[atomIndex];
          if (container->bFactorList)
            ai->b = container->bFactorList[atomIndex];
          if (container->occupancyList)
            ai->q = container->occupancyList[atomIndex];
          if (container->altLocList)
            ai->alt[0] = container->altLocList[atomIndex];

          AtomInfoAssignParameters(G, ai);
          AtomInfoAssignColors(G, ai);

          VLACheck(cset->IdxToAtm, int, cset->NIndex);
          VLACheck(cset->Coord, float, 3 * cset->NIndex + 2);
          cset->IdxToAtm[cset->NIndex] = atomIndex;
          float * coord = cset->coordPtr(cset->NIndex);
          cset->NIndex++;

          coord[0] = container->xCoordList[atomIndex];
          coord[1] = container->yCoordList[atomIndex];
          coord[2] = container->zCoordList[atomIndex];
        }
      }
    }

    VLASize(cset->Coord, float, 3 * cset->NIndex);
    VLASize(cset->IdxToAtm, int, cset->NIndex);
  }

  if (atomIndex != I->NAtom) {
    PRINTFB(G, FB_ObjectMolecule, FB_Warnings)
      " MMTF-Inconsistency-Warning: atom count mismatch (%d != numAtoms=%d)\n",
      atomIndex, I->NAtom ENDFB(G);
  }

  AtomInfoPurge(G, &tai);

  for (int l = 0; l < container->bondAtomListCount / 2; ++l) {
    BondTypeInit2(I->Bond + (I->NBond++),
        container->bondAtomList[l * 2],
        container->bondAtomList[l * 2 + 1],
        container->bondOrderListCount ? container->bondOrderList[l] : 1);
  }

  if (container->numBonds != I->NBond) {
    PRINTFB(G, FB_ObjectMolecule, FB_Warnings)
      " MMTF-Inconsistency-Warning: bond count mismatch (%d != numBonds=%d)\n",
      I->NBond, container->numBonds ENDFB(G);
  }

  if (discrete < 1) {
    ObjectMoleculeSetDiscrete(G, I, 0);
  } else {
    I->updateAtmToIdx();
  }

  // assemblies
  if (!I->DiscreteFlag) {
    CoordSet ** assembly_csets = get_assembly_csets(G, container,
        I->AtomInfo, I->CSet[0]);

    ObjectMoleculeSetAssemblyCSets(I, assembly_csets);
  }

  MMTF_container_free(container);

  ObjectMoleculeInvalidate(I, cRepAll, cRepInvAll, -1);
  ObjectMoleculeUpdateNonbonded(I);
  ObjectMoleculeAutoDisableAtomNameWildcard(I);

  return I;
#endif
}
