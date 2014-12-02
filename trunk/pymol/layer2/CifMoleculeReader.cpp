/*
 * Read a molecule from CIF
 *
 * (c) 2014 Schrodinger, Inc.
 */

#include "os_predef.h"
#include "os_std.h"

#include "MemoryDebug.h"
#include "Err.h"

#include "Base.h"
#include "Util.h"
#include "Scene.h"
#include "ObjectMolecule.h"
#include "CifFile.h"

#include <string>
#include <iostream>
#include <stdexcept>

/*
 * Add one bond without checking if it already exists
 */
static void ObjectMoleculeAddBond2(ObjectMolecule * I, int i1, int i2, int order) {
  VLACheck(I->Bond, BondType, I->NBond);
  BondTypeInit2(I->Bond + I->NBond, i1, i2, order);
  I->NBond++;
}

/*
 * Get the distance between two atoms in ObjectMolecule
 */
static float GetDistance(ObjectMolecule * I, int i1, int i2) {
  CoordSet *cset1, *cset2;
  int *AtmToIdx;

  if (I->DiscreteFlag) {
    cset1 = I->DiscreteCSet[i1];
    cset2 = I->DiscreteCSet[i2];
    AtmToIdx = I->DiscreteAtmToIdx;
  } else {
    cset1 = cset2 = I->CSet[0];
    AtmToIdx = cset1->AtmToIdx;
  }

  float v[3];
  subtract3f(
      cset1->Coord + AtmToIdx[i1] * 3,
      cset2->Coord + AtmToIdx[i2] * 3, v);
  return length3f(v);
}

/*
 * Bond order string to int
 */
static int bondOrderLookup(const char * order) {
  switch (order[0]) {
    case 'a': case 'A': // arom
      return 4;
    case 't': case 'T': // triple
      return 3;
    case 'd': case 'D':
      switch (order[1]) {
        case 'e': case 'E': // deloc
          return 4;
      }
      // double
      return 2;
  }
  // single
  return 1;
}

/*
 * Datastructure for bond_dict
 * bond_dict[resn][name1][name2] = order
 */
typedef std::map<std::string,
        std::map<std::string,
        std::map<std::string, int> > > bond_dict_t;
bond_dict_t bond_dict;

/*
 * parse components.cif into dictionary
 */
void update_components_bond_dict() {
  const char *name1, *name2;
  int order_value;
  const cif_array *arr_id_1, *arr_id_2, *arr_order;

  if (bond_dict.size())
    return;

  const char *filename = getenv("COMPONENTS_CIF");
  if (!filename || !filename[0])
    filename = "components.cif";

  cif_file * cif = new cif_file(filename);

  for (m_str_cifdatap_t::iterator data_it = cif->datablocks.begin(),
      data_it_end = cif->datablocks.end(); data_it != data_it_end; ++data_it) {

    const std::string &resn = data_it->first;
    cif_data * data = data_it->second;

    if( !(arr_id_1  = data->get_arr("_chem_comp_bond.atom_id_1")) ||
        !(arr_id_2  = data->get_arr("_chem_comp_bond.atom_id_2")) ||
        !(arr_order = data->get_arr("_chem_comp_bond.value_order")))
      continue;

    int nrows = arr_id_1->get_nrows();

    for (int i = 0; i < nrows; i++) {
      name1 = arr_id_1->as_s(i);
      name2 = arr_id_2->as_s(i);

      // make sure name1 < name2
      if (strcmp(name1, name2) < 0) {
        const char *tmp = name1; name1 = name2; name2 = tmp;
      }

      const char *order = arr_order->as_s(i);
      order_value = bondOrderLookup(order);

      bond_dict[resn][name1][name2] = order_value;
    }
  }
}

/*
 * Add bonds for one residue, with atoms spanning from i_start to i_end-1,
 * based on components.cif
 */
void ConnectComponent(ObjectMolecule * I, int i_start, int i_end) {
  if (i_end - i_start < 2)
    return;

  AtomInfoType *a1, *a2, *ai = I->AtomInfo;
  int order;
  bond_dict_t::const_iterator d_it;

  // get residue bond dictionary
  d_it = bond_dict.find(ai[i_start].resn);
  if (d_it == bond_dict.end())
    return;

  // for all pairs of atoms in given set
  for (int i1 = i_start + 1; i1 < i_end; i1++) {
    for (int i2 = i_start; i2 < i1; i2++) {
      a1 = ai + i1;
      a2 = ai + i2;

      // don't connect different alt codes
      if (a1->alt[0] && a2->alt[0] && strcmp(a1->alt, a2->alt) != 0) {
        continue;
      }

      // make sure name1 < name2
      if (strcmp(a1->name, a2->name) < 0) {
        AtomInfoType *atmp = a1; a1 = a2; a2 = atmp;
      }

      // lookup if atoms are bonded
      try {
        order = d_it->second.at(a1->name).at(a2->name);
      } catch (const std::out_of_range& e) {
        continue;
      }

      // make bond
      ObjectMoleculeAddBond2(I, i1, i2, order);
    }
  }
}

/*
 * Add intra residue bonds based on components.cif, and common polymer
 * connecting bonds (C->N, O3*->P)
 */
int ObjectMoleculeConnectComponents(ObjectMolecule * I)
{
  PyMOLGlobals * G = I->Obj.G;
  int i_start = 0, i_prev_c = 0, i_prev_o3 = 0;

  // read components.cif
  update_components_bond_dict();

  // user feedback if components.cif not available
  if (bond_dict.empty()) {
    PRINTFB(G, FB_ObjectMolecule, FB_Errors)
      " Error: Please download 'components.cif' from http://www.wwpdb.org/ccd.html\n"
      " and place it in the current directory or set the COMPONENTS_CIF environment"
      " variable.\n"
      ENDFB(G);
    return false;
  }

  // reserve some memory for new bonds
  if (!I->Bond) {
    I->Bond = VLACalloc(BondType, I->NAtom * 4);
  } else {
    VLACheck(I->Bond, BondType, I->NAtom * 4);
  }

  for (int i = 0; i < I->NAtom; i++) {
    const char *name = I->AtomInfo[i].name;

    // intra-residue
    if(!AtomInfoSameResidue(G, I->AtomInfo + i_start, I->AtomInfo + i)) {
      ConnectComponent(I, i_start, i);
      i_start = i;
    }

    // inter-residue polymer bonds
    if (strcmp("C", name) == 0) {
      i_prev_c = i;
    } else if (strncmp("O3", name, 2) == 0 && (name[2] == '*' || name[2] == '\'')) {
      // name in ('O3*', "O3'")
      i_prev_o3 = i;
    } else {
      int i_prev =
        (strcmp("N", name) == 0) ? i_prev_c :
        (strcmp("P", name) == 0) ? i_prev_o3 : -1;

      if (i_prev >= 0 && !AtomInfoSameResidue(G,
            I->AtomInfo + i_prev, I->AtomInfo + i)
          && GetDistance(I, i_prev, i) < 1.8) {
        // make bond
        ObjectMoleculeAddBond2(I, i_prev, i, 1);
      }
    }
  }

  // intra-residue (last residue)
  ConnectComponent(I, i_start, I->NAtom);

  // clean up
  VLASize(I->Bond, BondType, I->NBond);

  return true;
}

/*
 * secondary structure hash
 */
class sshashkey {
public:
  int asym_id;
  std::string resi;
  sshashkey() {};
  sshashkey(int asym_id_, const char *resi_, const char *ins_code = NULL)
    : asym_id(asym_id_), resi(resi_) {
    if (ins_code)
      resi.append(ins_code);
  }
  int compare(const sshashkey &other) const {
    int test = resi.compare(other.resi);
    if (test == 0)
      test = (asym_id - other.asym_id);
    return test;
  }
  bool operator<(const sshashkey &other) const { return compare(other) < 0; }
  bool operator>(const sshashkey &other) const { return compare(other) > 0; }
};
class sshashvalue {
public:
  char ss;
  sshashkey end;
};
typedef std::map<sshashkey, sshashvalue> sshashmap;
void sshashmap_clear(PyMOLGlobals * G, sshashmap &ssrecords) {
  // decrement Lexicon references (should go into ~sshashkey(), but
  // the PyMOLGlobals is not known there)
  for (sshashmap::iterator it = ssrecords.begin(),
      it_end = ssrecords.end(); it != it_end; ++it) {
    LexDec(G, it->first.asym_id);
    LexDec(G, it->second.end.asym_id);
  }
  ssrecords.clear();
}

/*
 * Read CELL and SYMMETRY
 */
CSymmetry * read_symmetry(PyMOLGlobals * G, cif_data * data) {
  cif_array * cell[6] = {
    data->get_arr("_cell?length_a"),
    data->get_arr("_cell?length_b"),
    data->get_arr("_cell?length_c"),
    data->get_arr("_cell?angle_alpha"),
    data->get_arr("_cell?angle_beta"),
    data->get_arr("_cell?angle_gamma")
  };

  for (int i = 0; i < 6; i++)
    if (cell[i] == NULL)
      return NULL;

  CSymmetry * symmetry = SymmetryNew(G);
  if (!symmetry)
    return NULL;

  for (int i = 0; i < 3; i++) {
    symmetry->Crystal->Dim[i] = cell[i]->as_d();
    symmetry->Crystal->Angle[i] = cell[i + 3]->as_d();
  }

  strncpy(symmetry->SpaceGroup,
      data->get_opt("_symmetry?space_group_name_h-m")->as_s(),
      WordLength - 1);

  symmetry->PDBZValue = data->get_opt("_cell.z_pdb")->as_i(1);

  return symmetry;
}

/*
 * Read CHEM_COMP_ATOM
 */
CoordSet ** read_chem_comp_atom_model(PyMOLGlobals * G, cif_data * data,
    AtomInfoType ** atInfoPtr) {

  cif_array * arr_x = data->get_arr("_chem_comp_atom.model_cartn_x");
  cif_array * arr_y = data->get_arr("_chem_comp_atom.model_cartn_y");
  cif_array * arr_z = data->get_arr("_chem_comp_atom.model_cartn_z");

  if (!arr_x || !arr_y || !arr_z) {
    arr_x = data->get_arr("_chem_comp_atom.pdbx_model_cartn_x_ideal");
    arr_y = data->get_arr("_chem_comp_atom.pdbx_model_cartn_y_ideal");
    arr_z = data->get_arr("_chem_comp_atom.pdbx_model_cartn_z_ideal");

    if (!arr_x || !arr_y || !arr_z) {
      return false;
    }
  }

  cif_array * arr_name            = data->get_opt("_chem_comp_atom.atom_id");
  cif_array * arr_symbol          = data->get_opt("_chem_comp_atom.type_symbol");
  cif_array * arr_resn            = data->get_opt("_chem_comp_atom.comp_id");
  cif_array * arr_partial_charge  = data->get_opt("_chem_comp_atom.partial_charge");
  cif_array * arr_formal_charge   = data->get_opt("_chem_comp_atom.charge");

  int nrows = arr_x->get_nrows();
  AtomInfoType *ai;
  int atomCount = 0, nAtom = nrows;
  float * coord = VLAlloc(float, 3 * nAtom);

  for (int i = 0; i < nrows; i++) {
    VLACheck(*atInfoPtr, AtomInfoType, atomCount);
    ai = *atInfoPtr + atomCount;
    memset((void*) ai, 0, sizeof(AtomInfoType));

    ai->rank = atomCount;
    ai->id = atomCount + 1;

    strncpy(ai->name, arr_name->as_s(i), cAtomNameLen);
    strncpy(ai->resn, arr_resn->as_s(i), cResnLen);
    strncpy(ai->elem, arr_symbol->as_s(i), cElemNameLen);

    ai->partialCharge = arr_partial_charge->as_d(i);
    ai->formalCharge = arr_formal_charge->as_i(i);

    ai->hetatm = 1;

    memset((void*) ai->visRep, 0, sizeof(ai->visRep));
    ai->visRep[cRepLine] = true;
    ai->visRep[cRepNonbonded] = true;

    AtomInfoAssignParameters(G, ai);
    AtomInfoAssignColors(G, ai);

    coord[atomCount * 3 + 0] = arr_x->as_d(i);
    coord[atomCount * 3 + 1] = arr_y->as_d(i);
    coord[atomCount * 3 + 2] = arr_z->as_d(i);

    atomCount++;
  }

  VLASize(coord, float, 3 * atomCount);
  VLASize(*atInfoPtr, AtomInfoType, atomCount);

  CoordSet ** csets = VLACalloc(CoordSet*, 1);
  csets[0] = CoordSetNew(G);
  csets[0]->NIndex = atomCount;
  csets[0]->Coord = coord;

  return csets;
}

/*
 * Read ATOM_SITE
 */
CoordSet ** read_atom_site(PyMOLGlobals * G, cif_data * data,
    AtomInfoType ** atInfoPtr, short * fractional) {

  cif_array *arr_x, *arr_y, *arr_z;
  cif_array *arr_name, *arr_resn, *arr_resi, *arr_chain, *arr_symbol,
            *arr_group_pdb, *arr_alt, *arr_ins_code, *arr_b, *arr_u,
            *arr_q, *arr_ID, *arr_mod_num, *arr_entity_id, *arr_segi;

  if ((arr_x = data->get_arr("_atom_site?cartn_x")) &&
      (arr_y = data->get_arr("_atom_site?cartn_y")) &&
      (arr_z = data->get_arr("_atom_site?cartn_z"))) {
    *fractional = 0;
  } else if (
      (arr_x = data->get_arr("_atom_site?fract_x")) &&
      (arr_y = data->get_arr("_atom_site?fract_y")) &&
      (arr_z = data->get_arr("_atom_site?fract_z"))) {
    *fractional = 1;
  } else {
    return NULL;
  }

  arr_name        = data->get_opt("_atom_site.auth_atom_id",
                                  "_atom_site.label_atom_id",
                                  "_atom_site_label");
  arr_resn        = data->get_opt("_atom_site.auth_comp_id",
                                  "_atom_site.label_comp_id");
  arr_resi        = data->get_opt("_atom_site.auth_seq_id",
                                  "_atom_site.label_seq_id");
  arr_chain       = data->get_arr("_atom_site.auth_asym_id");
  arr_segi        = data->get_opt("_atom_site.label_asym_id");
  arr_symbol      = data->get_opt("_atom_site?type_symbol");
  arr_group_pdb   = data->get_opt("_atom_site.group_pdb");
  arr_alt         = data->get_opt("_atom_site.label_alt_id");
  arr_ins_code    = data->get_opt("_atom_site.pdbx_pdb_ins_code");
  arr_b           = data->get_opt("_atom_site?b_iso_or_equiv");
  arr_u           = data->get_arr("_atom_site?u_iso_or_equiv"); // NULL
  arr_q           = data->get_opt("_atom_site?occupancy");
  arr_ID          = data->get_opt("_atom_site.id",
                                  "_atom_site_label");
  arr_mod_num     = data->get_opt("_atom_site.pdbx_pdb_model_num");
  arr_entity_id   = data->get_arr("_atom_site.label_entity_id"); // NULL

  if (!arr_chain)
    arr_chain = arr_segi;

  int nrows = arr_x->get_nrows();
  const char * resi;
  AtomInfoType *ai;
  int atomCount = 0;
  int first_model_num = arr_mod_num->as_i(0);

  for (int i = 0, n = nrows; i < n; i++) {
    if (arr_mod_num->as_i(i) != first_model_num)
      continue;

    VLACheck(*atInfoPtr, AtomInfoType, atomCount);
    ai = *atInfoPtr + atomCount;

    ai->rank = atomCount;
    ai->alt[0] = arr_alt->as_s(i)[0];

    ai->id = arr_ID->as_i(i);
    ai->b = (arr_u != NULL) ?
             arr_u->as_d(i) * 78.95683520871486 : // B = U * 8 * pi^2
             arr_b->as_d(i);
    ai->q = arr_q->as_d(i, 1.0);

    strncpy(ai->name, arr_name->as_s(i), cAtomNameLen);
    strncpy(ai->resn, arr_resn->as_s(i), cResnLen);
    strncpy(ai->elem, arr_symbol->as_s(i), cElemNameLen);
    strncpy(ai->segi, arr_segi->as_s(i), cSegiLen);

    ai->chain = LexIdx(G, arr_chain->as_s(i));

    ai->hetatm = 'H' == arr_group_pdb->as_s(i)[0];

    resi = arr_resi->as_s(i);
    ai->resv = atoi(resi);
    strncpy(ai->resi, resi, cResnLen);
    UtilNConcat(ai->resi, arr_ins_code->as_s(i), sizeof(ResIdent));

    memset((void*) ai->visRep, 0, sizeof(ai->visRep));
    ai->visRep[cRepLine] = true;
    ai->visRep[cRepNonbonded] = true;

    AtomInfoAssignParameters(G, ai);
    AtomInfoAssignColors(G, ai);

    if (arr_entity_id != NULL) {
      /* BEGIN PROPRIETARY CODE SEGMENT (see disclaimer in "os_proprietary.h") */
      /* END PROPRIETARY CODE SEGMENT */
    }

    atomCount++;
  }

  VLASize(*atInfoPtr, AtomInfoType, atomCount);

  float * coord = NULL;
  int ncsets = 0, mod_num, current_mod_num = 0;
  CoordSet ** csets = VLACalloc(CoordSet*, nrows / atomCount);

  for (int i = 0; i < nrows; i++) {
    mod_num = arr_mod_num->as_i(i, 1);

    if (i % atomCount == 0) {
      VLACheck(csets, CoordSet*, ncsets);
      csets[ncsets] = CoordSetNew(G);
      csets[ncsets]->NIndex = atomCount;
      csets[ncsets]->Coord = coord = VLAlloc(float, 3 * atomCount);
      ncsets++;
    } else if (current_mod_num != mod_num) {
      std::cout << "PDBX_PDB_MODEL_NUM problem" << std::endl;
    }

    current_mod_num = mod_num;

    *(coord++) = arr_x->as_d(i);
    *(coord++) = arr_y->as_d(i);
    *(coord++) = arr_z->as_d(i);
  }

  VLASize(csets, CoordSet*, ncsets);

  return csets;
}

/*
 * Read missing residues.
 * Append CA atoms to atInfoPtr, with no modification to coord sets. Sorting
 * will be necesarry to move those atoms to the correct place in the sequence.
 */
bool read_pdbx_unobs_or_zero_occ_residues(PyMOLGlobals * G, cif_data * data,
    AtomInfoType ** atInfoPtr) {

  cif_array *arr_resn, *arr_resi, *arr_chain, *arr_segi,
            *arr_poly_flag, *arr_ins_code, *arr_mod_num;

  if((arr_resn    = data->get_arr("_pdbx_unobs_or_zero_occ_residues.auth_comp_id",
                                  "_pdbx_unobs_or_zero_occ_residues.label_comp_id")) == NULL ||
     (arr_resi    = data->get_arr("_pdbx_unobs_or_zero_occ_residues.auth_seq_id",
                                  "_pdbx_unobs_or_zero_occ_residues.label_seq_id")) == NULL)
    return false;

  arr_poly_flag   = data->get_opt("_pdbx_unobs_or_zero_occ_residues.polymer_flag");
  arr_ins_code    = data->get_opt("_pdbx_unobs_or_zero_occ_residues.pdb_ins_code");
  arr_mod_num     = data->get_opt("_pdbx_unobs_or_zero_occ_residues.pdb_model_num");
  arr_segi        = data->get_opt("_pdbx_unobs_or_zero_occ_residues.label_asym_id");
  arr_chain       = data->get_arr("_pdbx_unobs_or_zero_occ_residues.auth_asym_id");

  if (!arr_chain)
    arr_chain = arr_segi;

  int nrows = arr_resn->get_nrows();
  const char * resi;
  AtomInfoType *ai;
  int atomCount = VLAGetSize(*atInfoPtr);
  int fake_id = 0;

  if (atomCount > 0)
    fake_id = (*atInfoPtr + atomCount - 1)->id;

  for (int i = 0, n = nrows; i < n; i++) {
    if (arr_mod_num->as_i(i, 1) != 1)
      continue;

    if ('N' == arr_poly_flag->as_s(i)[0])
      continue;

    VLACheck(*atInfoPtr, AtomInfoType, atomCount); // auto-zero
    ai = *atInfoPtr + atomCount;

    ai->rank = atomCount;
    ai->id = (++fake_id);

    strncpy(ai->name, "CA", cAtomNameLen);
    strncpy(ai->resn, arr_resn->as_s(i), cResnLen);
    ai->elem[0] = 'C';
    strncpy(ai->segi, arr_segi->as_s(i), cSegiLen);

    ai->chain = LexIdx(G, arr_chain->as_s(i));

    resi = arr_resi->as_s(i);
    ai->resv = atoi(resi);
    strncpy(ai->resi, resi, cResnLen);
    UtilNConcat(ai->resi, arr_ins_code->as_s(i), sizeof(ResIdent));

    AtomInfoAssignParameters(G, ai);
    AtomInfoAssignColors(G, ai);

    atomCount++;
  }

  VLASize(*atInfoPtr, AtomInfoType, atomCount);

  return true;
}

/*
 * Read secondary structure from STRUCT_CONF or STRUCT_SHEET_RANGE
 */
bool read_ss_(PyMOLGlobals * G, cif_data * data, char ss, sshashmap &ssrecords) {
  cif_array *arr_beg_chain, *arr_beg_resi,
            *arr_end_chain, *arr_end_resi;

  std::string prefix = "_struct_conf.";
  if (ss == 'S')
    prefix = "_struct_sheet_range.";

  if (!(arr_beg_chain = data->get_arr((prefix + "beg_auth_asym_id").c_str(),
                                      (prefix + "beg_label_asym_id").c_str())) ||
      !(arr_beg_resi  = data->get_arr((prefix + "beg_auth_seq_id").c_str(),
                                      (prefix + "beg_label_seq_id").c_str())) ||
      !(arr_end_chain = data->get_arr((prefix + "end_auth_asym_id").c_str(),
                                      (prefix + "end_label_asym_id").c_str())) ||
      !(arr_end_resi  = data->get_arr((prefix + "end_auth_seq_id").c_str(),
                                      (prefix + "end_label_seq_id").c_str())))
    return false;

  cif_array *arr_beg_ins_code = data->get_opt((prefix + "pdbx_beg_pdb_ins_code").c_str());
  cif_array *arr_end_ins_code = data->get_opt((prefix + "pdbx_end_pdb_ins_code").c_str());

  int nrows = arr_beg_chain->get_nrows();

  for (int i = 0; i < nrows; i++) {
    sshashkey key(
      LexIdx(G, arr_beg_chain->as_s(i)),
      arr_beg_resi->as_s(i),
      arr_beg_ins_code->as_s(i));

    sshashvalue &value = ssrecords[key];
    value.ss = ss;
    new (&value.end) sshashkey(
        LexIdx(G, arr_end_chain->as_s(i)),
        arr_end_resi->as_s(i),
        arr_end_ins_code->as_s(i));
  }

  return true;
}

/*
 * Read secondary structure
 */
bool read_ss(PyMOLGlobals * G, cif_data * datablock, AtomInfoType * atInfo) {
  sshashmap ssrecords;

  read_ss_(G, datablock, 'H', ssrecords);
  read_ss_(G, datablock, 'S', ssrecords);

  if (ssrecords.empty())
    return false;

  AtomInfoType *aj, *ai, *atoms_end = atInfo + VLAGetSize(atInfo);

  for (ai = atInfo; ai < atoms_end; ai++) {
    if (strcmp(ai->name, "CA"))
      continue;

    sshashkey key(ai->chain, ai->resi);
    sshashmap::iterator it = ssrecords.find(key);

    if (it == ssrecords.end())
      continue;

    sshashvalue &value = it->second;

    for (aj = ai; aj < atoms_end; aj++) {
      if (strcmp(aj->name, "CA"))
        continue;

      aj->ssType[0] = value.ss;

      if (value.end.resi == aj->resi && value.end.asym_id == aj->chain)
        break;
    }
  }

  // decrement Lexicon references (should go into ~sshashkey())
  sshashmap_clear(G, ssrecords);

  return true;
}

/*
 * Read anisotropic temperature factors from ATOM_SITE or ATOM_SITE_ANISOTROP
 */
bool read_atom_site_aniso(PyMOLGlobals * G, cif_data * data,
    AtomInfoType * atInfo) {

  cif_array *arr_label, *arr_u11, *arr_u22, *arr_u33, *arr_u12, *arr_u13, *arr_u23;
  bool mmcif = true;
  float factor = 1.0;

  if ((arr_label = data->get_arr("_atom_site_anisotrop.id", "_atom_site.id"))) {
    // mmCIF, assume _atom_site_id is numeric and look up by atom ID
    // Warning: according to mmCIF spec, id can be any alphanumeric string
  } else if ((arr_label = data->get_arr("_atom_site_aniso_label"))) {
    // small molecule CIF, lookup by atom name
    mmcif = false;
  } else {
    return false;
  }

  if ((arr_u11 = data->get_arr("_atom_site_anisotrop.u[1][1]", "_atom_site_aniso_u_11", "_atom_site.aniso_u[1][1]"))) {
    // U
    arr_u22 = data->get_opt("_atom_site_anisotrop.u[2][2]", "_atom_site_aniso_u_22", "_atom_site.aniso_u[2][2]");
    arr_u33 = data->get_opt("_atom_site_anisotrop.u[3][3]", "_atom_site_aniso_u_33", "_atom_site.aniso_u[3][3]");
    arr_u12 = data->get_opt("_atom_site_anisotrop.u[1][2]", "_atom_site_aniso_u_12", "_atom_site.aniso_u[1][2]");
    arr_u13 = data->get_opt("_atom_site_anisotrop.u[1][3]", "_atom_site_aniso_u_13", "_atom_site.aniso_u[1][3]");
    arr_u23 = data->get_opt("_atom_site_anisotrop.u[2][3]", "_atom_site_aniso_u_23", "_atom_site.aniso_u[2][3]");
  } else if (
      (arr_u11 = data->get_arr("_atom_site_anisotrop.b[1][1]", "_atom_site_aniso_b_11", "_atom_site.aniso_b[1][1]"))) {
    // B
    factor = 0.012665147955292222; // U = B / (8 * pi^2)
    arr_u22 = data->get_opt("_atom_site_anisotrop.b[2][2]", "_atom_site_aniso_b_22", "_atom_site.aniso_b[2][2]");
    arr_u33 = data->get_opt("_atom_site_anisotrop.b[3][3]", "_atom_site_aniso_b_33", "_atom_site.aniso_b[3][3]");
    arr_u12 = data->get_opt("_atom_site_anisotrop.b[1][2]", "_atom_site_aniso_b_12", "_atom_site.aniso_b[1][2]");
    arr_u13 = data->get_opt("_atom_site_anisotrop.b[1][3]", "_atom_site_aniso_b_13", "_atom_site.aniso_b[1][3]");
    arr_u23 = data->get_opt("_atom_site_anisotrop.b[2][3]", "_atom_site_aniso_b_23", "_atom_site.aniso_b[2][3]");
  } else {
    return false;
  }

  AtomInfoType *ai;
  int nAtom = VLAGetSize(atInfo);

  std::map<int, AtomInfoType*> id_dict;
  std::map<std::string, AtomInfoType*> name_dict;

  // build dictionary
  for (int i = 0; i < nAtom; i++) {
    ai = atInfo + i;
    if (mmcif) {
      id_dict[ai->id] = ai;
    } else {
      std::string key(ai->name);
      name_dict[key] = ai;
    }
  }

  // read aniso table
  for (int i = 0; i < arr_u11->get_nrows(); i++) {
    try {
      if (mmcif) {
        int key = arr_label->as_i(i);
        ai = id_dict.at(key);
      } else {
        std::string key(arr_label->as_s(i));
        ai = name_dict.at(key);
      }
    } catch (const std::out_of_range& e) {
      std::cout << "atom lookup failed" << std::endl;
      continue;
    }

    ai->U11 = arr_u11->as_d(i) * factor;
    ai->U22 = arr_u22->as_d(i) * factor;
    ai->U33 = arr_u33->as_d(i) * factor;
    ai->U12 = arr_u12->as_d(i) * factor;
    ai->U13 = arr_u13->as_d(i) * factor;
    ai->U23 = arr_u23->as_d(i) * factor;
  }

  return true;
}

/*
 * Read GEOM_BOND
 */
bool read_geom_bond_atom_site_labels(PyMOLGlobals * G, cif_data * data,
    AtomInfoType * atInfo, CoordSet * cset) {

  cif_array *arr_ID_1, *arr_ID_2;
  if ((arr_ID_1 = data->get_arr("_geom_bond.atom_site_id_1",
                                "_geom_bond_atom_site_label_1")) == NULL ||
      (arr_ID_2 = data->get_arr("_geom_bond.atom_site_id_2",
                                "_geom_bond_atom_site_label_2")) == NULL)
    return false;

  cif_array *arr_symm_1 = data->get_opt("_geom_bond?site_symmetry_1");
  cif_array *arr_symm_2 = data->get_opt("_geom_bond?site_symmetry_2");

  int nrows = arr_ID_1->get_nrows();
  int nAtom = VLAGetSize(atInfo);
  int nBond = 0;

  BondType *bond = cset->TmpBond = VLACalloc(BondType, 6 * nAtom);

  std::map<std::string, int> name_dict;

  // build dictionary
  for (int i = 0; i < nAtom; i++) {
    std::string key(atInfo[i].name);
    name_dict[key] = i;
  }

  // read table
  for (int i = 0; i < nrows; i++) {
    if (strcmp(arr_symm_1->as_s(i),
               arr_symm_2->as_s(i)))
      // don't bond to symmetry mates
      continue;

    std::string key1(arr_ID_1->as_s(i));
    std::string key2(arr_ID_2->as_s(i));

    try {
      int i1 = name_dict.at(key1);
      int i2 = name_dict.at(key2);

      nBond++;
      BondTypeInit2(bond++, i1, i2, 1);

    } catch (const std::out_of_range& e) {
      std::cout << "name lookup failed" << std::endl;
    }
  }

  if (nBond) {
    VLASize(cset->TmpBond, BondType, nBond);
    cset->NTmpBond = nBond;
  } else {
    VLAFreeP(cset->TmpBond);
  }

  return true;
}

/*
 * Read bonds from STRUCT_CONN
 */
bool read_struct_conn_(PyMOLGlobals * G, cif_data * data,
    AtomInfoType * atInfo, CoordSet * cset) {

  cif_array *col_type_id;
  cif_array *col_asym_id[2], *col_comp_id[2], *col_seq_id[2], *col_atom_id[2];
  cif_array *col_alt_id[2], *col_ins_code[2], *col_symm[2];

  if ((col_type_id    = data->get_arr("_struct_conn.conn_type_id")) == NULL ||
      (col_asym_id[0] = data->get_arr("_struct_conn.ptnr1_auth_asym_id",
                                      "_struct_conn.ptnr1_label_asym_id")) == NULL ||
      (col_comp_id[0] = data->get_arr("_struct_conn.ptnr1_auth_comp_id",
                                      "_struct_conn.ptnr1_label_comp_id")) == NULL ||
      (col_seq_id[0]  = data->get_arr("_struct_conn.ptnr1_auth_seq_id",
                                      "_struct_conn.ptnr1_label_seq_id")) == NULL ||
      (col_atom_id[0] = data->get_arr("_struct_conn.ptnr1_label_atom_id")) == NULL ||
      (col_asym_id[1] = data->get_arr("_struct_conn.ptnr2_auth_asym_id",
                                      "_struct_conn.ptnr2_label_asym_id")) == NULL ||
      (col_comp_id[1] = data->get_arr("_struct_conn.ptnr2_auth_comp_id",
                                      "_struct_conn.ptnr2_label_comp_id")) == NULL ||
      (col_seq_id[1]  = data->get_arr("_struct_conn.ptnr2_auth_seq_id",
                                      "_struct_conn.ptnr2_label_seq_id")) == NULL ||
      (col_atom_id[1] = data->get_arr("_struct_conn.ptnr2_label_atom_id")) == NULL)
    return false;

  col_alt_id[0]   = data->get_opt("_struct_conn.pdbx_ptnr1_label_alt_id");
  col_ins_code[0] = data->get_opt("_struct_conn.pdbx_ptnr1_pdb_ins_code");
  col_symm[0]     = data->get_opt("_struct_conn.ptnr1_symmetry");
  col_alt_id[1]   = data->get_opt("_struct_conn.pdbx_ptnr2_label_alt_id");
  col_ins_code[1] = data->get_opt("_struct_conn.pdbx_ptnr2_pdb_ins_code");
  col_symm[1]     = data->get_opt("_struct_conn.ptnr2_symmetry");

  int nrows = col_type_id->get_nrows();
  int nAtom = VLAGetSize(atInfo);
  int nBond = 0;

  BondType *bond = cset->TmpBond = VLACalloc(BondType, 6 * nAtom);

  std::map<std::string, int> name_dict;

  for (int i = 0; i < nAtom; i++) {
    std::string key(LexStr(G, atInfo[i].chain));
    key += '/';
    key += atInfo[i].resn;
    key += '/';
    key += atInfo[i].resi;
    key += '/';
    key += atInfo[i].name;
    key += '/';
    key += atInfo[i].alt;
    name_dict[key] = i;
  }

  for (int i = 0; i < nrows; i++) {
    const char * type_id = col_type_id->as_s(i);
    if (strncasecmp(type_id, "covale", 6) && strcasecmp(type_id, "modres"))
      // ignore non-covalent bonds (metalc, hydrog)
      continue;
    if (strcmp(col_symm[0]->as_s(i),
               col_symm[1]->as_s(i)))
      // don't bond to symmetry mates
      continue;

    std::string key[2];
    for (int j = 0; j < 2; j++) {
      key[j] += col_asym_id[j]->as_s(i);
      key[j] += '/';
      key[j] += col_comp_id[j]->as_s(i);
      key[j] += '/';
      key[j] += col_seq_id[j]->as_s(i);
      key[j] += col_ins_code[j]->as_s(i);
      key[j] += '/';
      key[j] += col_atom_id[j]->as_s(i);
      key[j] += '/';
      key[j] += col_alt_id[j]->as_s(i);
    }

    try {
      int i1 = name_dict.at(key[0]);
      int i2 = name_dict.at(key[1]);

      nBond++;
      BondTypeInit2(bond++, i1, i2, 1);

    } catch (const std::out_of_range& e) {
      std::cout << "name lookup failed" << std::endl;
    }
  }

  if (nBond) {
    VLASize(cset->TmpBond, BondType, nBond);
    cset->NTmpBond = nBond;
  } else {
    VLAFreeP(cset->TmpBond);
  }

  return true;
}

/*
 * Read bonds from CHEM_COMP_BOND
 */
bool read_chem_comp_bond_atom_ids(PyMOLGlobals * G, cif_data * data,
    AtomInfoType * atInfo, CoordSet * cset) {

  cif_array *col_ID_1, *col_ID_2, *col_comp_id;

  if ((col_ID_1    = data->get_arr("_chem_comp_bond.atom_id_1")) == NULL ||
      (col_ID_2    = data->get_arr("_chem_comp_bond.atom_id_2")) == NULL ||
      (col_comp_id = data->get_arr("_chem_comp_bond.comp_id")) == NULL)
    return false;

  cif_array *col_order = data->get_opt("_chem_comp_bond.value_order");

  int nrows = col_ID_1->get_nrows();
  int nAtom = VLAGetSize(atInfo);
  int nBond = 0;

  BondType *bond = cset->TmpBond = VLACalloc(BondType, 6 * nAtom);

  std::map<std::string, int> name_dict;

  for (int i = 0; i < nAtom; i++) {
    std::string key(atInfo[i].name);
    name_dict[key] = i;
  }

  for (int i = 0; i < nrows; i++) {
    std::string key1(col_ID_1->as_s(i));
    std::string key2(col_ID_2->as_s(i));
    const char * order = col_order->as_s(i);

    try {
      int i1 = name_dict.at(key1);
      int i2 = name_dict.at(key2);
      int order_value = bondOrderLookup(order);

      nBond++;
      BondTypeInit2(bond++, i1, i2, order_value);

    } catch (const std::out_of_range& e) {
      std::cout << "name lookup failed" << std::endl;
    }
  }

  if (nBond) {
    VLASize(cset->TmpBond, BondType, nBond);
    cset->NTmpBond = nBond;
  } else {
    VLAFreeP(cset->TmpBond);
  }

  return true;
}

/*
 * Read a molecule from a data block
 *
 * States are returned as a coord set VLA, atoms are written to atInfoPtr, and
 * bonds are stored with first coordinate set.
 */
CoordSet ** ObjectMoleculeCifData2CoordSets(PyMOLGlobals * G, cif_data * datablock,
    AtomInfoType ** atInfoPtr, short * fractional) {

  CoordSet ** csets = NULL;

  if ((csets = read_atom_site(G, datablock, atInfoPtr, fractional))) {
    // bonds
    if (
        read_geom_bond_atom_site_labels(G, datablock, *atInfoPtr, csets[0]) ||
        read_struct_conn_(G, datablock, *atInfoPtr, csets[0]) ||
        0);

    // anisou
    read_atom_site_aniso(G, datablock, *atInfoPtr);

    // secondary structure
    read_ss(G, datablock, *atInfoPtr);

    // missing residues
    read_pdbx_unobs_or_zero_occ_residues(G, datablock, atInfoPtr);

  } else if ((csets = read_chem_comp_atom_model(G, datablock, atInfoPtr))) {
    // bonds
    read_chem_comp_bond_atom_ids(G, datablock, *atInfoPtr, csets[0]);
  }

  return csets;
}

ObjectMolecule *ObjectMoleculeReadCifStr(PyMOLGlobals * G, ObjectMolecule * I,
                                      char *st, int frame,
                                      int discrete, int quiet, int multiplex,
                                      char *new_name)
{
  CoordSet ** csets = NULL;
  AtomInfoType *atInfo;
  int isNew = !I;
  int nAtom, ncsets;
  short fractional = false;

  cif_data * datablock;
  cif_file * cif = new cif_file(st, 2);
  ok_assert(1, cif);
  ok_assert(2, cif->datablocks.size());
  datablock = cif->datablocks.begin()->second;

  if (discrete > 0 || multiplex > 0) {
    PRINTFB(G, FB_ObjectMolecule, FB_Errors)
      " Error: discrete and multiplex not yet supported for CIF\n" ENDFB(G);
    return NULL;
  }

  // allocate either ObjectMolecule or AtomInfo
  if(isNew) {
    I = ObjectMoleculeNew(G, (discrete > 0));
    I->Obj.Color = AtomInfoUpdateAutoColor(G);
    atInfo = I->AtomInfo;
  } else {
    atInfo = VLACalloc(AtomInfoType, 10);
  }

  // read coordsets from datablock
  csets = ObjectMoleculeCifData2CoordSets(G, datablock, &atInfo, &fractional);
  ok_assert(3, csets);

  // get number of atoms and coordinate sets
  nAtom = csets[0]->NIndex;
  ncsets = VLAGetSize(csets);

  // get frame where to insert coordsets
  if(frame < 0)
    frame = I->NCSet;

  // initialize the new coordsets (not data, but indices, etc.)
  for (int i = 0; i < ncsets; i++) {
    csets[i]->Obj = I;
    csets[i]->enumIndices();
  }

  // get atoms into ObjectMolecule
  if(isNew) {
    I->AtomInfo = atInfo;   /* IMPORTANT to reassign: this VLA may have moved! */
    I->NAtom = VLAGetSize(atInfo);
  } else {
    // NOTE:
    // - will sort cset
    // - will sort atInfo
    // - will release atInfo
    ObjectMoleculeMerge(I, atInfo, csets[0], false, cAIC_MOLMask, false);
    atInfo = NULL;

    // transfer sorting from cset[0] to cset[1..n]
    for (int i = 1; i < ncsets; i++) {
      for(int j = 0; j < nAtom; j++) {
        csets[i]->IdxToAtm[j] = csets[0]->IdxToAtm[j];
        csets[i]->AtmToIdx[j] = csets[0]->AtmToIdx[j];
      }
    }
  }

  // get coordinate sets into ObjectMolecule
  for (int i = 0; i < ncsets; i++) {
    int j = frame + i;

    if(I->NCSet <= j)
      I->NCSet = j + 1;

    VLACheck(I->CSet, CoordSet *, j);

    // if this state already exists, free it
    if(I->CSet[j])
      I->CSet[j]->fFree();

    I->CSet[j] = csets[i];
  }

  // handle symmetry and update fractional -> cartesian
  I->Symmetry = read_symmetry(G, datablock);
  if (I->Symmetry) {
    SymmetryAttemptGeneration(I->Symmetry, false);

    if(fractional && I->Symmetry->Crystal) {
      CrystalUpdate(I->Symmetry->Crystal);
      for (int i = 0; i < ncsets; i++) {
        CoordSetFracToReal(csets[i], I->Symmetry->Crystal);
      }
    }
  }

  // grow IdxToAtm and AtmToIdx if more atoms than coordinates
  ObjectMoleculeExtendIndices(I, frame);

  // sort atoms internally
  ObjectMoleculeSort(I);

  // create bonds
  if(isNew) {
    if(SettingGetGlobal_i(G, cSetting_connect_mode) == 4) {
      ObjectMoleculeConnectComponents(I);
    } else {
      ObjectMoleculeConnect(I, &I->NBond, &I->Bond, I->AtomInfo, csets[0], true, -1);
    }
  }

  // computationally intense update tasks
  SceneCountFrames(G);
  ObjectMoleculeInvalidate(I, cRepAll, cRepInvAll, -1);
  if (isNew) {
    ObjectMoleculeUpdateIDNumbers(I);
    ObjectMoleculeUpdateNonbonded(I);
  }

  VLAFreeP(csets);
  delete cif;
  return I;
ok_except3:
  PRINTFB(G, FB_ObjectMolecule, FB_Errors)
    " Error: no coordinates found in CIF\n" ENDFB(G);
  if(isNew) {
    ObjectMoleculeFree(I);
  } else {
    VLAFreeP(atInfo);
  }
  return NULL;
ok_except2:
  PRINTFB(G, FB_ObjectMolecule, FB_Errors)
    " Error: empty CIF\n" ENDFB(G);
  delete cif;
ok_except1:
  return NULL;
}

// vi:sw=2:ts=2:expandtab
