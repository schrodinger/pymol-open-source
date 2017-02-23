/*
 * Read a molecule from CIF
 *
 * (c) 2014 Schrodinger, Inc.
 */

#include <algorithm>
#include <string>
#include <map>
#include <set>
#include <vector>
#include <memory>

#include "os_predef.h"
#include "os_std.h"

#include "MemoryDebug.h"
#include "Err.h"

#include "AssemblyHelpers.h"
#include "AtomInfo.h"
#include "Base.h"
#include "Executive.h"
#include "P.h"
#include "Util.h"
#include "Scene.h"
#include "Rep.h"
#include "ObjectMolecule.h"
#include "CifFile.h"
#include "CifBondDict.h"
#include "Util2.h"
#include "Vector.h"
#include "Lex.h"

// canonical amino acid three letter codes
const char * aa_three_letter[] = {
  "ALA", // A
  "ASX", // B  for ambiguous asparagine/aspartic-acid
  "CYS", // C
  "ASP", // D
  "GLU", // E
  "PHE", // F
  "GLY", // G
  "HIS", // H
  "ILE", // I
  NULL,  // J
  "LYS", // K
  "LEU", // L
  "MET", // M
  "ASN", // N
  "HOH", // O  for water
  "PRO", // P
  "GLN", // Q
  "ARG", // R
  "SER", // S
  "THR", // T
  NULL,  // U
  "VAL", // V
  "TRP", // W
  NULL,  // X  for other
  "TYR", // Y
  "GLX", // Z  for ambiguous glutamine/glutamic acid
};

// amino acid one-to-three letter code translation
static const char * aa_get_three_letter(char aa) {
  if (aa < 'A' || aa > 'Z')
    return "UNK";
  const char * three = aa_three_letter[aa - 'A'];
  return (three) ? three : "UNK";
}

// dictionary content types
enum CifDataType {
  CIF_UNKNOWN,
  CIF_CORE,       // small molecule
  CIF_MMCIF,      // macromolecular structure
  CIF_CHEM_COMP   // chemical component
};

// simple 1-indexed string storage
class seqvec_t : public std::vector<std::string> {
public:
  void set(int i, const char * mon_id) {
    if (i < 1) {
      printf("error: i(%d) < 1\n", i);
      return;
    }
    if (i > size())
      resize(i);
    (*this)[i - 1] = mon_id;
  }

  const char * get(int i) const {
    if (i < 1 || i > size())
      return NULL;
    return (*this)[i - 1].c_str();
  }
};

// structure to collect information about a data block
struct CifContentInfo {
  PyMOLGlobals * G;
  CifDataType type;
  bool fractional;
  bool use_auth;
  std::set<lexidx_t> chains_filter;
  std::set<std::string> polypeptide_entities; // entity ids
  std::map<std::string, seqvec_t> sequences;  // entity_id -> [resn1, resn2, ...]

  bool is_excluded_chain(const char * chain) {
    if (chains_filter.empty())
      return false;

    auto result = OVLexicon_BorrowFromCString(G->Lexicon, chain);
    if (OVreturn_IS_OK(result))
      return is_excluded_chain(result.word);

    return false;
  }

  bool is_excluded_chain(lexidx_t chain) {
    return (!chains_filter.empty() &&
        chains_filter.count(chain) == 0);
  }

  bool is_polypeptide(const char * entity_id) {
    return polypeptide_entities.count(entity_id);
  }

  CifContentInfo(PyMOLGlobals * G, bool use_auth=true) :
    G(G),
    type(CIF_UNKNOWN),
    fractional(false),
    use_auth(use_auth) {}
};

/*
 * Make a string key that represents the collection of alt id, asym id,
 * atom id, comp id and seq id components of the label for a macromolecular
 * atom site.
 */
static std::string make_mm_atom_site_label(PyMOLGlobals * G, AtomInfoType * a) {
  char resi[8];
  AtomResiFromResv(resi, sizeof(resi), a);

  std::string key(LexStr(G, a->chain));
  key += '/';
  key += LexStr(G, a->resn);
  key += '/';
  key += resi;
  key += '/';
  key += LexStr(G, a->name);
  key += '/';
  key += a->alt;
  return key;
}
static std::string make_mm_atom_site_label(PyMOLGlobals * G, const char * asym_id,
    const char * comp_id, const char * seq_id, const char * ins_code,
    const char * atom_id, const char * alt_id) {
  std::string key(asym_id);
  key += '/';
  key += comp_id;
  key += '/';
  key += seq_id;
  key += ins_code;
  key += '/';
  key += atom_id;
  key += '/';
  key += alt_id;
  return key;
}

/*
 * Get first non-NULL element
 */
template <typename T>
static T VLAGetFirstNonNULL(T * vla) {
  int n = VLAGetSize(vla);
  for (int i = 0; i < n; ++i)
    if (vla[i])
      return vla[i];
  return NULL;
}

/*
 * Lookup one key in a map, return true if found and
 * assign output reference `value1`
 */
template <typename Map, typename Key, typename T>
inline bool find1(Map& dict, T& value1, const Key& key1) {
  auto it = dict.find(key1);
  if (it == dict.end())
    return false;
  value1 = it->second;
  return true;
}

/*
 * Lookup two keys in a map, return true if both found and
 * assign output references `value1` and `value2`.
 */
template <typename Map, typename Key, typename T>
inline bool find2(Map& dict,
    T& value1, const Key& key1,
    T& value2, const Key& key2) {
  if (!find1(dict, value1, key1))
    return false;
  if (!find1(dict, value2, key2))
    return false;
  return true;
}

static void AtomInfoSetEntityId(PyMOLGlobals * G, AtomInfoType * ai, const char * entity_id) {
  ai->custom = LexIdx(G, entity_id);

#ifdef _PYMOL_IP_EXTRAS
  PropertySet(G, ai, "entity_id", entity_id);
#endif
}

/*
 * Add one bond without checking if it already exists
 */
static void ObjectMoleculeAddBond2(ObjectMolecule * I, int i1, int i2, int order) {
  VLACheck(I->Bond, BondType, I->NBond);
  BondTypeInit2(I->Bond + I->NBond, i1, i2, order);
  I->NBond++;
}

/*
 * Distance based connectivity for discrete objects
 */
static void ObjectMoleculeConnectDiscrete(ObjectMolecule * I) {
  for (int i = 0; i < I->NCSet; i++) {
    if (!I->CSet[i])
      continue;

    int nbond = 0;
    BondType * bond = NULL;

    ObjectMoleculeConnect(I, &nbond, &bond, I->AtomInfo, I->CSet[i], true, 3);

    if (!bond)
      continue;

    if (!I->Bond) {
      I->Bond = bond;
    } else {
      VLASize(I->Bond, BondType, I->NBond + nbond);
      memcpy(I->Bond + I->NBond, bond, nbond * sizeof(*bond));
      VLAFreeP(bond);
    }

    I->NBond += nbond;
  }
}

/*
 * Get the distance between two atoms in ObjectMolecule
 */
static float GetDistance(ObjectMolecule * I, int i1, int i2) {
  const CoordSet *cset;
  int idx1 = -1, idx2 = -1;

  // find first coordset which contains both atoms
  if (I->DiscreteFlag) {
    cset = I->DiscreteCSet[i1];
    if (cset == I->DiscreteCSet[i2]) {
      idx1 = I->DiscreteAtmToIdx[i1];
      idx2 = I->DiscreteAtmToIdx[i2];
    }
  } else {
    for (int i = 0; i < I->NCSet; ++i) {
      if ((cset = I->CSet[i])) {
        if ((idx1 = cset->AtmToIdx[i1]) != -1 &&
            (idx2 = cset->AtmToIdx[i2]) != -1) {
          break;
        }
      }
    }
  }

  if (idx1 == -1 || idx2 == -1)
    return 999.f;

  float v[3];
  subtract3f(
      cset->coordPtr(idx1),
      cset->coordPtr(idx2), v);
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
 * Read bonds from CHEM_COMP_BOND in `bond_dict` dictionary
 */
static bool read_chem_comp_bond_dict(const cif_data * data, bond_dict_t &bond_dict) {
  const cif_array *arr_id_1, *arr_id_2, *arr_order, *arr_comp_id;

  if( !(arr_id_1  = data->get_arr("_chem_comp_bond.atom_id_1")) ||
      !(arr_id_2  = data->get_arr("_chem_comp_bond.atom_id_2")) ||
      !(arr_order = data->get_arr("_chem_comp_bond.value_order")) ||
      !(arr_comp_id = data->get_arr("_chem_comp_bond.comp_id"))) {

    if ((arr_comp_id = data->get_arr("_chem_comp_atom.comp_id"))) {
      // atom(s) but no bonds (e.g. metals)
      bond_dict.set_unknown(arr_comp_id->as_s());
      return true;
    }

    return false;
  }

  const char *name1, *name2, *resn;
  int order_value;
  int nrows = arr_id_1->get_nrows();

  for (int i = 0; i < nrows; i++) {
    resn = arr_comp_id->as_s(i);
    name1 = arr_id_1->as_s(i);
    name2 = arr_id_2->as_s(i);

    const char *order = arr_order->as_s(i);
    order_value = bondOrderLookup(order);

    bond_dict.set(resn, name1, name2, order_value);
  }

  return true;
}

/*
 * parse $PYMOL_DATA/chem_comp_bond-top100.cif (subset of components.cif) into
 * a static (global) dictionary.
 */
static bond_dict_t * get_global_components_bond_dict(PyMOLGlobals * G) {
  static bond_dict_t bond_dict;

  if (bond_dict.empty()) {
    const char * pymol_data = getenv("PYMOL_DATA");
    if (!pymol_data || !pymol_data[0])
      return NULL;

    std::string path(pymol_data);
    path.append(PATH_SEP).append("chem_comp_bond-top100.cif");
    cif_file cif(path.c_str());

    for (m_str_cifdatap_t::iterator data_it = cif.datablocks.begin(),
        data_it_end = cif.datablocks.end(); data_it != data_it_end; ++data_it) {
      read_chem_comp_bond_dict(data_it->second, bond_dict);
    }
  }

  return &bond_dict;
}

/*
 * True for N-H1 and N-H3, those are not in the chemical components dictionary.
 */
static bool is_N_H1_or_H3(PyMOLGlobals * G,
    const AtomInfoType * a1,
    const AtomInfoType * a2) {
  if (a2->name == G->lex_const.N) {
    a2 = a1;
  } else if (a1->name != G->lex_const.N) {
    return false;
  }

  return (a2->name == G->lex_const.H1 || a2->name == G->lex_const.H3);
}

/*
 * Add bonds for one residue, with atoms spanning from i_start to i_end-1,
 * based on components.cif
 */
static void ConnectComponent(ObjectMolecule * I, int i_start, int i_end,
    bond_dict_t * bond_dict) {

  if (i_end - i_start < 2)
    return;

  auto G = I->Obj.G;
  AtomInfoType *a1, *a2, *ai = I->AtomInfo;
  int order;

  // get residue bond dictionary
  auto res_dict = bond_dict->get(G, LexStr(G, ai[i_start].resn));
  if (res_dict == NULL)
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

      // restart if we hit the next residue in bulk solvent (atoms must
      // not be sorted for this)
      // TODO artoms are sorted at this point
      if (a1->name == a2->name) {
        i_start = i1;
        break;
      }

      // lookup if atoms are bonded
      order = res_dict->get(LexStr(G, a1->name), LexStr(G, a2->name));
      if (order < 0) {
        if (!is_N_H1_or_H3(G, a1, a2) || GetDistance(I, i1, i2) > 1.2)
          continue;

        order = 1;
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
static int ObjectMoleculeConnectComponents(ObjectMolecule * I,
    bond_dict_t * bond_dict=NULL) {

  PyMOLGlobals * G = I->Obj.G;
  int i_start = 0, i_prev_c = 0, i_prev_o3 = 0;

  if (!bond_dict) {
    // read components.cif
    if (!(bond_dict = get_global_components_bond_dict(G)))
      return false;
  }

  // reserve some memory for new bonds
  if (!I->Bond) {
    I->Bond = VLACalloc(BondType, I->NAtom * 4);
  } else {
    VLACheck(I->Bond, BondType, I->NAtom * 4);
  }

  for (int i = 0;; ++i) {
    // intra-residue
    if(!AtomInfoSameResidue(G, I->AtomInfo + i_start, I->AtomInfo + i)) {
      ConnectComponent(I, i_start, i, bond_dict);
      i_start = i;
    }

    if (i == I->NAtom)
      break;

    // ignore alt coords for inter-residue bonding
    if (I->AtomInfo[i].alt[0] && I->AtomInfo[i].alt[0] != 'A')
      continue;

    const char *name = LexStr(G, I->AtomInfo[i].name);

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

  // clean up
  VLASize(I->Bond, BondType, I->NBond);

  return true;
}

/*
 * secondary structure hash
 */
class sshashkey {
public:
  lexidx_t chain; // borrowed ref
  int resv;
  char inscode;

  void assign(int asym_id_, int resv_, char ins_code_ = '\0') {
    chain = asym_id_;
    resv = resv_;
    inscode = ins_code_;
  }

  // comparable to sshashkey and AtomInfoType
  template <typename T> int compare(const T &other) const {
    int test = resv - other.resv;
    if (test == 0) {
      test = (chain - other.chain);
      if (test == 0)
        test = inscode - other.inscode;
    }
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

// std::array for pre-C++11
template <typename T, size_t N>
class myarray {
  T m_data[N];
  public:
  T * data() { return m_data; }
};

// PDBX_STRUCT_OPER_LIST type
typedef std::map<std::string, myarray<float, 16> > oper_list_t;

// type for parsed PDBX_STRUCT_OPER_LIST
typedef std::vector<std::vector<std::string> > oper_collection_t;

/*
 * Parse operation expressions like (1,2)(3-6)
 */
static oper_collection_t parse_oper_expression(const std::string &expr) {
  using namespace std;

  oper_collection_t collection;

  // first step to split parenthesized chunks
  vector<string> a_vec = strsplit(expr, ')');

  // loop over chunks (still include leading '(')
  for (auto a_it = a_vec.begin(); a_it != a_vec.end(); ++a_it) {
    const char * a_chunk = a_it->c_str();

    // finish chunk
    while (*a_chunk == '(')
      ++a_chunk;

    // skip empty chunks
    if (!*a_chunk)
      continue;

    collection.resize(collection.size() + 1);
    oper_collection_t::reference ids = collection.back();

    // split chunk by commas
    vector<string> b_vec = strsplit(a_chunk, ',');

    // look for ranges
    for (vector<string>::iterator
        b_it = b_vec.begin();
        b_it != b_vec.end(); ++b_it) {
      // "c_d" will have either one (no range) or two items
      vector<string> c_d = strsplit(*b_it, '-');

      ids.push_back(c_d[0]);

      if (c_d.size() == 2)
        for (int i = atoi(c_d[0].c_str()) + 1,
                 j = atoi(c_d[1].c_str()) + 1; i < j; ++i)
        {
          char i_str[16];
          snprintf(i_str, sizeof(i_str), "%d", i);
          ids.push_back(i_str);
        }
    }
  }

  return collection;
}

/*
 * Get chains which are part of the assembly
 *
 * assembly_chains: output set
 * assembly_id: ID of the assembly or NULL to use first assembly
 */
static bool get_assembly_chains(PyMOLGlobals * G,
    const cif_data * data,
    std::set<lexidx_t> &assembly_chains,
    const char * assembly_id) {

  const cif_array *arr_id, *arr_asym_id_list;

  if ((arr_id           = data->get_arr("_pdbx_struct_assembly_gen.assembly_id")) == NULL ||
      (arr_asym_id_list = data->get_arr("_pdbx_struct_assembly_gen.asym_id_list")) == NULL)
    return false;

  for (int i = 0, nrows = arr_id->get_nrows(); i < nrows; ++i) {
    if (strcmp(assembly_id, arr_id->as_s(i)))
      continue;

    const char * asym_id_list = arr_asym_id_list->as_s(i);
    std::vector<std::string> chains = strsplit(asym_id_list, ',');
    for (auto it = chains.begin(); it != chains.end(); ++it) {
      assembly_chains.insert(LexIdx(G, it->c_str()));
    }
  }

  return !assembly_chains.empty();
}

/*
 * Read assembly
 *
 * atInfo: atom info array to use for chain check
 * cset: template coordinate set to create assembly coordsets from
 * assembly_id: assembly identifier
 *
 * return: assembly coordinates as VLA of coordinate sets
 */
static
CoordSet ** read_pdbx_struct_assembly(PyMOLGlobals * G,
    const cif_data * data,
    const AtomInfoType * atInfo,
    const CoordSet * cset,
    const char * assembly_id) {

  const cif_array *arr_id, *arr_assembly_id, *arr_oper_expr, *arr_asym_id_list;


  if ((arr_id           = data->get_arr("_pdbx_struct_oper_list.id")) == NULL ||
      (arr_assembly_id  = data->get_arr("_pdbx_struct_assembly_gen.assembly_id")) == NULL ||
      (arr_oper_expr    = data->get_arr("_pdbx_struct_assembly_gen.oper_expression")) == NULL ||
      (arr_asym_id_list = data->get_arr("_pdbx_struct_assembly_gen.asym_id_list")) == NULL)
    return NULL;

  const cif_array * arr_matrix[] = {
    data->get_opt("_pdbx_struct_oper_list.matrix[1][1]"),
    data->get_opt("_pdbx_struct_oper_list.matrix[1][2]"),
    data->get_opt("_pdbx_struct_oper_list.matrix[1][3]"),
    data->get_opt("_pdbx_struct_oper_list.vector[1]"),
    data->get_opt("_pdbx_struct_oper_list.matrix[2][1]"),
    data->get_opt("_pdbx_struct_oper_list.matrix[2][2]"),
    data->get_opt("_pdbx_struct_oper_list.matrix[2][3]"),
    data->get_opt("_pdbx_struct_oper_list.vector[2]"),
    data->get_opt("_pdbx_struct_oper_list.matrix[3][1]"),
    data->get_opt("_pdbx_struct_oper_list.matrix[3][2]"),
    data->get_opt("_pdbx_struct_oper_list.matrix[3][3]"),
    data->get_opt("_pdbx_struct_oper_list.vector[3]")
  };

  // build oper_list from _pdbx_struct_oper_list

  oper_list_t oper_list;

  for (int i = 0, nrows = arr_id->get_nrows(); i < nrows; ++i) {
    float * matrix = oper_list[arr_id->as_s(i)].data();

    identity44f(matrix);

    for (int j = 0; j < 12; ++j) {
      matrix[j] = arr_matrix[j]->as_d(i);
    }
  }

  CoordSet ** csets = NULL;
  int csetbeginidx = 0;

  // assembly
  for (int i = 0, nrows = arr_oper_expr->get_nrows(); i < nrows; ++i) {
    if (strcmp(assembly_id, arr_assembly_id->as_s(i)))
      continue;

    const char * oper_expr    = arr_oper_expr->as_s(i);
    const char * asym_id_list = arr_asym_id_list->as_s(i);

    oper_collection_t collection = parse_oper_expression(oper_expr);
    std::vector<std::string> chains = strsplit(asym_id_list, ',');
    std::set<lexidx_t> chains_set;
    for (auto it = chains.begin(); it != chains.end(); ++it) {
      auto result = OVLexicon_BorrowFromCString(G->Lexicon, it->c_str());
      if (OVreturn_IS_OK(result)) {
        chains_set.insert(result.word);
      }
    }

    // new coord set VLA
    int ncsets = 1;
    for (auto c_it = collection.begin(); c_it != collection.end(); ++c_it) {
      ncsets *= c_it->size();
    }

    if (!csets) {
      csets = VLACalloc(CoordSet*, ncsets);
    } else {
      csetbeginidx = VLAGetSize(csets);
      VLASize(csets, CoordSet*, csetbeginidx + ncsets);
    }

    // for cartesian product
    int c_src_len = 1;

    // coord set for subset of atoms
    CoordSet ** c_csets = csets + csetbeginidx;
    c_csets[0] = CoordSetCopyFilterChains(cset, atInfo, chains_set);

    // build new coord sets
    for (auto c_it = collection.rbegin(); c_it != collection.rend(); ++c_it) {
      // copy
      int j = c_src_len;
      while (j < c_src_len * c_it->size()) {
        // cartesian product
        for (int k = 0; k < c_src_len; ++k, ++j) {
          c_csets[j] = CoordSetCopy(c_csets[k]);
        }
      }

      // transform
      j = 0;
      for (auto s_it = c_it->begin(); s_it != c_it->end(); ++s_it) {
        const float * matrix = oper_list[*s_it].data();

        // cartesian product
        for (int k = 0; k < c_src_len; ++k, ++j) {
          CoordSetTransform44f(c_csets[j], matrix);
        }
      }

      // cartesian product
      // Note: currently, "1m4x" seems to be the only structure in the PDB
      // which uses a cartesian product expression
      c_src_len *= c_it->size();
    }
  }

  // return assembly coordsets
  return csets;
}

/*
 * Set ribbon_trace_atoms and cartoon_trace_atoms for CA/P only models
 */
static bool read_pdbx_coordinate_model(PyMOLGlobals * G, cif_data * data, ObjectMolecule * mol) {
  const cif_array * arr_type = data->get_arr("_pdbx_coordinate_model.type");
  const cif_array * arr_asym = data->get_arr("_pdbx_coordinate_model.asym_id");

  if (!arr_type || !arr_asym)
    return false;

  // affected chains
  std::set<const char*, strless2_t> asyms;

  // collect CA/P-only chain identifiers
  for (int i = 0, nrows = arr_type->get_nrows(); i < nrows; ++i) {
    const char * type = arr_type->as_s(i);
    // no need anymore to check "CA ATOMS ONLY", since nonbonded CA are
    // now (v1.8.2) detected automatically in RepCartoon and RepRibbon
    if (strcmp(type, "P ATOMS ONLY") == 0) {
      asyms.insert(arr_asym->as_s(i));
    }
  }

  if (asyms.empty())
    return false;

  // set on atom-level
  for (int i = 0, nrows = VLAGetSize(mol->AtomInfo); i < nrows; ++i) {
    AtomInfoType * ai = mol->AtomInfo + i;
    if (asyms.count(LexStr(G, ai->segi))) {
      SettingSet(G, cSetting_cartoon_trace_atoms, true, ai);
      SettingSet(G, cSetting_ribbon_trace_atoms,  true, ai);
    }
  }

  return true;
}

/*
 * Read CELL and SYMMETRY
 */
static CSymmetry * read_symmetry(PyMOLGlobals * G, cif_data * data) {
  const cif_array * cell[6] = {
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

  symmetry->PDBZValue = data->get_opt("_cell.z_pdb")->as_i(0, 1);

  // register symmetry operations if given
  const cif_array * arr_as_xyz = data->get_arr(
      "_symmetry_equiv?pos_as_xyz",
      "_space_group_symop?operation_xyz");
  if (arr_as_xyz) {
    std::vector<std::string> sym_op;
    for (int i = 0, n = arr_as_xyz->get_nrows(); i < n; ++i) {
      sym_op.push_back(arr_as_xyz->as_s(i));
    }
    SymmetrySpaceGroupRegister(G, symmetry->SpaceGroup, sym_op);
  }

  return symmetry;
}

/*
 * Read CHEM_COMP_ATOM
 */
static CoordSet ** read_chem_comp_atom_model(PyMOLGlobals * G, cif_data * data,
    AtomInfoType ** atInfoPtr) {

  const cif_array *arr_x, *arr_y = NULL, *arr_z = NULL;

  if ((arr_x = data->get_arr("_chem_comp_atom.pdbx_model_cartn_x_ideal"))) {
    arr_y = data->get_arr("_chem_comp_atom.pdbx_model_cartn_y_ideal");
    arr_z = data->get_arr("_chem_comp_atom.pdbx_model_cartn_z_ideal");
  } else if ((arr_x = data->get_arr("_chem_comp_atom.x"))) {
    arr_y = data->get_arr("_chem_comp_atom.y");
    arr_z = data->get_arr("_chem_comp_atom.z");
  } else if ((arr_x = data->get_arr("_chem_comp_atom.model_cartn_x"))) {
    arr_y = data->get_arr("_chem_comp_atom.model_cartn_y");
    arr_z = data->get_arr("_chem_comp_atom.model_cartn_z");
  }

  if (!arr_x || !arr_y || !arr_z) {
    return NULL;
  }

  PRINTFB(G, FB_Executive, FB_Details)
    " ExecutiveLoad-Detail: Detected chem_comp CIF\n" ENDFB(G);

  const cif_array * arr_name            = data->get_opt("_chem_comp_atom.atom_id");
  const cif_array * arr_symbol          = data->get_opt("_chem_comp_atom.type_symbol");
  const cif_array * arr_resn            = data->get_opt("_chem_comp_atom.comp_id");
  const cif_array * arr_partial_charge  = data->get_opt("_chem_comp_atom.partial_charge");
  const cif_array * arr_formal_charge   = data->get_opt("_chem_comp_atom.charge");
  const cif_array * arr_stereo          = data->get_opt("_chem_comp_atom.pdbx_stereo_config");

  int nrows = arr_x->get_nrows();
  AtomInfoType *ai;
  int atomCount = 0, nAtom = nrows;
  float * coord = VLAlloc(float, 3 * nAtom);
  int auto_show = RepGetAutoShowMask(G);

  for (int i = 0; i < nrows; i++) {
    if (arr_x->is_missing(i))
      continue;

    VLACheck(*atInfoPtr, AtomInfoType, atomCount);
    ai = *atInfoPtr + atomCount;
    memset((void*) ai, 0, sizeof(AtomInfoType));

    ai->rank = atomCount;
    ai->id = atomCount + 1;

    LexAssign(G, ai->name, arr_name->as_s(i));
    LexAssign(G, ai->resn, arr_resn->as_s(i));
    strncpy(ai->elem, arr_symbol->as_s(i), cElemNameLen);

    ai->partialCharge = arr_partial_charge->as_d(i);
    ai->formalCharge = arr_formal_charge->as_i(i);

    ai->hetatm = 1;

    ai->visRep = auto_show;
    ai->stereo = convertCharToStereo(arr_stereo->as_s(i)[0]);

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
 * Map model number to state (1-based)
 */
class ModelStateMapper {
  bool remap;
  std::map<int, int> mapping;
public:
  ModelStateMapper(bool remap) : remap(remap) {}

  int operator()(int model) {
    if (!remap)
      return model;

    int state = mapping[model];

    if (!state) {
      state = mapping.size();
      mapping[model] = state;
    }

    return state;
  }
};

/*
 * Read ATOM_SITE
 *
 * atInfoPtr: atom info array to fill
 * info: data content configuration to populate with collected information
 *
 * return: models as VLA of coordinate sets
 */
static CoordSet ** read_atom_site(PyMOLGlobals * G, cif_data * data,
    AtomInfoType ** atInfoPtr, CifContentInfo &info, bool discrete) {

  const cif_array *arr_x, *arr_y, *arr_z;
  const cif_array *arr_name = NULL, *arr_resn = NULL, *arr_resi = NULL,
            *arr_chain = NULL, *arr_symbol,
            *arr_group_pdb, *arr_alt, *arr_ins_code = NULL, *arr_b, *arr_u,
            *arr_q, *arr_ID, *arr_mod_num, *arr_entity_id, *arr_segi;

  if ((arr_x = data->get_arr("_atom_site?cartn_x")) &&
      (arr_y = data->get_arr("_atom_site?cartn_y")) &&
      (arr_z = data->get_arr("_atom_site?cartn_z"))) {
  } else if (
      (arr_x = data->get_arr("_atom_site?fract_x")) &&
      (arr_y = data->get_arr("_atom_site?fract_y")) &&
      (arr_z = data->get_arr("_atom_site?fract_z"))) {
    info.fractional = true;
  } else {
    return NULL;
  }

  if (info.use_auth) {
    arr_name      = data->get_arr("_atom_site.auth_atom_id");
    arr_resn      = data->get_arr("_atom_site.auth_comp_id");
    arr_resi      = data->get_arr("_atom_site.auth_seq_id");
    arr_chain     = data->get_arr("_atom_site.auth_asym_id");
    arr_ins_code  = data->get_arr("_atom_site.pdbx_pdb_ins_code");
  }

  if (!arr_name) arr_name = data->get_arr("_atom_site.label_atom_id");
  if (!arr_resn) arr_resn = data->get_opt("_atom_site.label_comp_id");

  const cif_array *arr_label_seq_id = data->get_opt("_atom_site.label_seq_id");

  // PDBe provides unique seq_ids for bulk het groups
  if (!arr_resi) arr_resi = data->get_arr("_atom_site.pdbe_label_seq_id");
  if (!arr_resi) arr_resi = arr_label_seq_id;

  if (arr_name) {
    info.type = CIF_MMCIF;
    PRINTFB(G, FB_Executive, FB_Details)
      " ExecutiveLoad-Detail: Detected mmCIF\n" ENDFB(G);
  } else {
    arr_name      = data->get_opt("_atom_site_label");
    info.type = CIF_CORE;
    PRINTFB(G, FB_Executive, FB_Details)
      " ExecutiveLoad-Detail: Detected small molecule CIF\n" ENDFB(G);
  }

  arr_segi        = data->get_opt("_atom_site.label_asym_id");
  arr_symbol      = data->get_opt("_atom_site?type_symbol");
  arr_group_pdb   = data->get_opt("_atom_site.group_pdb");
  arr_alt         = data->get_opt("_atom_site.label_alt_id");
  arr_b           = data->get_opt("_atom_site?b_iso_or_equiv");
  arr_u           = data->get_arr("_atom_site?u_iso_or_equiv"); // NULL
  arr_q           = data->get_opt("_atom_site?occupancy");
  arr_ID          = data->get_opt("_atom_site.id",
                                  "_atom_site_label");
  arr_mod_num     = data->get_opt("_atom_site.pdbx_pdb_model_num");
  arr_entity_id   = data->get_arr("_atom_site.label_entity_id"); // NULL

  const cif_array * arr_color = data->get_arr("_atom_site.pymol_color");
  const cif_array * arr_reps  = data->get_arr("_atom_site.pymol_reps");
  const cif_array * arr_ss    = data->get_opt("_atom_site.pymol_ss");

  if (!arr_chain)
    arr_chain = arr_segi;

  ModelStateMapper model_to_state(!SettingGetGlobal_i(G, cSetting_pdb_honor_model_number));
  int nrows = arr_x->get_nrows();
  AtomInfoType *ai;
  int atomCount = 0;
  int auto_show = RepGetAutoShowMask(G);
  int first_model_num = model_to_state(arr_mod_num->as_i(0, 1));
  CoordSet * cset;
  int mod_num, ncsets = 0;

  // collect number of atoms per model and number of coord sets
  std::map<int, int> atoms_per_model;
  for (int i = 0, n = nrows; i < n; i++) {
    mod_num = model_to_state(arr_mod_num->as_i(i, 1));

    if (mod_num < 1) {
      PRINTFB(G, FB_ObjectMolecule, FB_Errors)
        " Error: model numbers < 1 not supported: %d\n", mod_num ENDFB(G);
      return NULL;
    }

    atoms_per_model[mod_num - 1] += 1;

    if (ncsets < mod_num)
      ncsets = mod_num;
  }

  // set up coordinate sets
  CoordSet ** csets = VLACalloc(CoordSet*, ncsets);
  for (auto it = atoms_per_model.begin(); it != atoms_per_model.end(); ++it) {
    csets[it->first] = cset = CoordSetNew(G);
    cset->Coord = VLAlloc(float, 3 * it->second);
    cset->IdxToAtm = VLAlloc(int, it->second);
  }

  // mm_atom_site_label -> atom index (1-indexed)
  std::map<std::string, int> name_dict;

  for (int i = 0, n = nrows; i < n; i++) {
    lexidx_t segi = LexIdx(G, arr_segi->as_s(i));

    if (info.is_excluded_chain(segi)) {
      LexDec(G, segi);
      continue;
    }

    mod_num = model_to_state(arr_mod_num->as_i(i, 1));

    // copy coordinates into coord set
    cset = csets[mod_num - 1];
    int idx = cset->NIndex++;
    float * coord = cset->coordPtr(idx);
    coord[0] = arr_x->as_d(i);
    coord[1] = arr_y->as_d(i);
    coord[2] = arr_z->as_d(i);

    if (!discrete && ncsets > 1) {
      // mm_atom_site_label aggregate
      std::string key = make_mm_atom_site_label(G,
          arr_chain->as_s(i),
          arr_resn->as_s(i),
          arr_resi->as_s(i),
          arr_ins_code ? arr_ins_code->as_s(i) : "",
          arr_name->as_s(i),
          arr_alt->as_s(i));

      // check if this is not a new atom
      if (mod_num != first_model_num) {
        int atm = name_dict[key] - 1;
        if (atm >= 0) {
          cset->IdxToAtm[idx] = atm;
          continue;
        }
      }

      name_dict[key] = atomCount + 1;
    }

    cset->IdxToAtm[idx] = atomCount;

    VLACheck(*atInfoPtr, AtomInfoType, atomCount);
    ai = *atInfoPtr + atomCount;

    ai->rank = atomCount;
    ai->alt[0] = arr_alt->as_s(i)[0];

    ai->id = arr_ID->as_i(i);
    ai->b = (arr_u != NULL) ?
             arr_u->as_d(i) * 78.95683520871486 : // B = U * 8 * pi^2
             arr_b->as_d(i);
    ai->q = arr_q->as_d(i, 1.0);

    strncpy(ai->elem, arr_symbol->as_s(i), cElemNameLen);

    ai->chain = LexIdx(G, arr_chain->as_s(i));
    ai->name = LexIdx(G, arr_name->as_s(i));
    ai->resn = LexIdx(G, arr_resn->as_s(i));
    ai->segi = segi; // steal reference

    if ('H' == arr_group_pdb->as_s(i)[0]) {
      ai->hetatm = 1;
      ai->flags = cAtomFlag_ignore;
    }

    ai->resv = arr_resi->as_i(i);
    ai->temp1 = arr_label_seq_id->as_i(i); // for add_missing_ca

    if (arr_ins_code) {
      ai->setInscode(arr_ins_code->as_s(i)[0]);
    }

    if (arr_reps) {
      ai->visRep = arr_reps->as_i(i, auto_show);
      ai->flags |= cAtomFlag_inorganic; // suppress auto_show_classified
    } else {
      ai->visRep = auto_show;
    }

    ai->ssType[0] = arr_ss->as_s(i)[0];

    AtomInfoAssignParameters(G, ai);

    if (arr_color) {
      ai->color = arr_color->as_i(i);
    } else {
      AtomInfoAssignColors(G, ai);
    }

    if (arr_entity_id != NULL) {
      AtomInfoSetEntityId(G, ai, arr_entity_id->as_s(i));
    }

    atomCount++;
  }

  VLASize(*atInfoPtr, AtomInfoType, atomCount);

  return csets;
}

/*
 * Update `info` with entity polymer information
 */
static bool read_entity_poly(PyMOLGlobals * G, const cif_data * data, CifContentInfo &info) {
  const cif_array *arr_entity_id = NULL, *arr_type = NULL,
        *arr_num = NULL, *arr_mon_id = NULL;

  if (!(arr_entity_id     = data->get_arr("_entity_poly.entity_id")) ||
      !(arr_type          = data->get_arr("_entity_poly.type")))
    return false;

  const cif_array * arr_seq_one_letter = data->get_arr("_entity_poly.pdbx_seq_one_letter_code");

  // polypeptides
  for (int i = 0, n = arr_entity_id->get_nrows(); i < n; i++) {
    if (!strncasecmp("polypeptide", arr_type->as_s(i), 11)) {
      const char * entity_id = arr_entity_id->as_s(i);
      info.polypeptide_entities.insert(entity_id);

      if (arr_seq_one_letter) {
        // sequences
        auto& entity_sequence = info.sequences[entity_id];
        const char * one = arr_seq_one_letter->as_s(i);
        for (int i = 0; *one; ++one) {
          if (strchr(" \t\r\n", *one)) // skip whitespace
            continue;

          if (*one == '(') {
            const char * end = strchr(one, ')');
            if (!end)
              break;

            std::string three(one + 1, end - one - 1);
            entity_sequence.set(++i, three.c_str());
            one = end;
          } else {
            entity_sequence.set(++i, aa_get_three_letter(*one));
          }
        }
      }
    }
  }

  if (!arr_seq_one_letter) {
    // sequences
    if ((arr_entity_id     = data->get_arr("_entity_poly_seq.entity_id")) &&
        (arr_num           = data->get_arr("_entity_poly_seq.num")) &&
        (arr_mon_id        = data->get_arr("_entity_poly_seq.mon_id"))) {
      for (int i = 0, n = arr_entity_id->get_nrows(); i < n; i++) {
        info.sequences[arr_entity_id->as_s(i)].set(
            arr_num->as_i(i),
            arr_mon_id->as_s(i));
      }
    }
  }

  return true;
}

#if 0
/*
 * Read missing residues.
 *
 * This is under the assumption that _pdbx_poly_seq_scheme has complete
 * mon_id, asym_id, and seq_id arrays, but auth_seq_num only for residues with
 * coordinates (present in _atom_site) and ? values for missing residues.
 *
 * Append CA atoms to atInfoPtr, with no modification to coord sets. Sorting
 * will be necesarry to move those atoms to the correct place in the sequence.
 */
static bool read_pdbx_poly_seq_scheme(PyMOLGlobals * G, const cif_data * data,
    AtomInfoType ** atInfoPtr, CifContentInfo &info) {

  const cif_array *arr_resn = NULL, *arr_resi = NULL, *arr_chain = NULL,
    *arr_segi = NULL, *arr_ins_code = NULL, *arr_auth_seq_num = NULL,
    *arr_entity_id = NULL;

  if (!(arr_resn          = data->get_arr("_pdbx_poly_seq_scheme.mon_id")) ||
      !(arr_segi          = data->get_arr("_pdbx_poly_seq_scheme.asym_id")) ||
      !(arr_entity_id     = data->get_arr("_pdbx_poly_seq_scheme.entity_id")) ||
      !(arr_auth_seq_num  = data->get_arr("_pdbx_poly_seq_scheme.auth_seq_num")))
    return false;

  if (info.use_auth) {
    arr_resi      = data->get_arr("_pdbx_poly_seq_scheme.pdb_seq_num");
    arr_chain     = data->get_arr("_pdbx_poly_seq_scheme.pdb_strand_id");
    arr_ins_code  = data->get_arr("_pdbx_poly_seq_scheme.pdb_ins_code");
  }

  if (!arr_resi && !(arr_resi = data->get_arr("_pdbx_poly_seq_scheme.seq_id")))
    return false;

  if (!arr_chain)
    arr_chain = arr_segi;

  int nrows = arr_resn->get_nrows();
  const char * resi;
  AtomInfoType *ai;
  int atomCount = VLAGetSize(*atInfoPtr);

  for (int i = 0, n = nrows; i < n; i++) {
    if (!arr_auth_seq_num->is_missing(i))
      continue;

    const char * segi = arr_segi->as_s(i);
    if (info.is_excluded_chain(segi))
      continue;

    const char * entity_id = arr_entity_id->as_s(i);
    if (!info.is_polypeptide(entity_id))
      continue;

    VLACheck(*atInfoPtr, AtomInfoType, atomCount); // auto-zero
    ai = *atInfoPtr + atomCount;

    ai->rank = atomCount;
    ai->id = -1;

    ai->elem[0] = 'C';

    ai->name = LexIdx(G, "CA");
    ai->resn = LexIdx(G, arr_resn->as_s(i));
    ai->segi = LexIdx(G, segi);
    ai->chain = LexIdx(G, arr_chain->as_s(i));

    ai->resv = arr_resi->as_i(i);

    if (arr_ins_code) {
      ai->setInscode(arr_ins_code->as_s(i)[0]);
    }

    AtomInfoAssignParameters(G, ai);
    AtomInfoAssignColors(G, ai);
    AtomInfoSetEntityId(G, ai, entity_id);

    atomCount++;
  }

  VLASize(*atInfoPtr, AtomInfoType, atomCount);

  return true;
}
#endif

/*
 * Sub-routine for `add_missing_ca`
 */
static void add_missing_ca_sub(PyMOLGlobals * G,
    AtomInfoType *& atInfo,
    int& current_resv,
    int& atomCount,
    const int i_ref, int resv,
    const seqvec_t * current_seq,
    const char * entity_id)
{
  if (!atInfo[i_ref].temp1)
    return;

  for (++current_resv; current_resv < resv; ++current_resv) {
    const char * resn = current_seq->get(current_resv);
    if (!resn)
      continue;

    VLACheck(atInfo, AtomInfoType, atomCount);

    AtomInfoType *ai = atInfo + atomCount;

    ai->rank = atomCount;
    ai->id = -1;

    ai->elem[0] = 'C';
    LexAssign(G, ai->name, "CA");
    LexAssign(G, ai->resn, resn);
    LexAssign(G, ai->segi, atInfo[i_ref].segi);
    LexAssign(G, ai->chain, atInfo[i_ref].chain);

    ai->temp1 = current_resv;
    ai->resv = current_resv + (atInfo[i_ref].resv - atInfo[i_ref].temp1);

    AtomInfoAssignParameters(G, ai);
    AtomInfoAssignColors(G, ai);
    AtomInfoSetEntityId(G, ai, entity_id);

    ++atomCount;
  }
}

/*
 * Read missing residues / full sequence
 *
 * This function relies on the label_seq_id numbering which must be available
 * in the `temp1` kludge field.
 *
 * Use the _entity_poly and _entity_poly_seq information to identify
 * missing residues in partially present chains. Add CA atoms for those
 * to present complete sequences in the sequence viewer.
 */
static bool add_missing_ca(PyMOLGlobals * G,
    AtomInfoType *& atInfo, CifContentInfo &info) {

  int oldAtomCount = VLAGetSize(atInfo);
  int atomCount = oldAtomCount;
  int current_resv = 0;
  const seqvec_t * current_seq = NULL;
  const char * current_entity_id = "";

  for (int i = 0; i < oldAtomCount; ++i) {
    const char * entity_id = LexStr(G, atInfo[i].custom);

    if (i == 0
        || atInfo[i].chain != atInfo[i - 1].chain
        || strcmp(entity_id, current_entity_id)) {
      // finish prev seq
      if (current_seq && i > 0) {
         add_missing_ca_sub(G,
             atInfo, current_resv, atomCount,
             i - 1, current_seq->size() + 1,
             current_seq, current_entity_id);
      }

      current_resv = 0;
      current_seq = NULL;
      current_entity_id = entity_id;

      if (info.is_polypeptide(entity_id) && !info.is_excluded_chain(atInfo[i].segi)) {
        // get new sequence
        auto it = info.sequences.find(entity_id);
        if (it != info.sequences.end()) {
          current_seq = &it->second;
        }
      }

    } else if (i > 0 && atInfo[i].temp1 == atInfo[i - 1].temp1) {
      continue;
    }

    if (current_seq) {
      add_missing_ca_sub(G,
          atInfo, current_resv, atomCount,
          i, atInfo[i].temp1,
          current_seq, entity_id);
    }
  }

  // finish last seq
  if (current_seq) {
    add_missing_ca_sub(G,
        atInfo, current_resv, atomCount,
        oldAtomCount - 1, current_seq->size() + 1,
        current_seq, current_entity_id);
  }

  VLASize(atInfo, AtomInfoType, atomCount);

  return true;
}

/*
 * Read secondary structure from STRUCT_CONF or STRUCT_SHEET_RANGE
 */
static bool read_ss_(PyMOLGlobals * G, cif_data * data, char ss,
    sshashmap &ssrecords, CifContentInfo &info)
{
  const cif_array *arr_beg_chain = NULL, *arr_beg_resi = NULL,
                  *arr_end_chain = NULL, *arr_end_resi = NULL,
                  *arr_beg_ins_code = NULL, *arr_end_ins_code = NULL;

  std::string prefix = "_struct_conf.";
  if (ss == 'S')
    prefix = "_struct_sheet_range.";

  if (info.use_auth &&
      (arr_beg_chain = data->get_arr((prefix + "beg_auth_asym_id").c_str())) &&
      (arr_beg_resi  = data->get_arr((prefix + "beg_auth_seq_id").c_str())) &&
      (arr_end_chain = data->get_arr((prefix + "end_auth_asym_id").c_str())) &&
      (arr_end_resi  = data->get_arr((prefix + "end_auth_seq_id").c_str()))) {
    // auth only
    arr_beg_ins_code = data->get_arr((prefix + "pdbx_beg_pdb_ins_code").c_str());
    arr_end_ins_code = data->get_arr((prefix + "pdbx_end_pdb_ins_code").c_str());
  } else if (
      !(arr_beg_chain = data->get_arr((prefix + "beg_label_asym_id").c_str())) ||
      !(arr_beg_resi  = data->get_arr((prefix + "beg_label_seq_id").c_str())) ||
      !(arr_end_chain = data->get_arr((prefix + "end_label_asym_id").c_str())) ||
      !(arr_end_resi  = data->get_arr((prefix + "end_label_seq_id").c_str()))) {
    return false;
  }

  const cif_array *arr_conf_type_id = (ss == 'S') ? NULL :
    data->get_arr("_struct_conf.conf_type_id");

  int nrows = arr_beg_chain->get_nrows();
  sshashkey key;

  for (int i = 0; i < nrows; i++) {
    // first character of conf_type_id (one of H, S, T)
    char ss_i = arr_conf_type_id ? arr_conf_type_id->as_s(i)[0] : ss;

    // exclude TURN_* (include HELX_* and STRN)
    if (ss_i == 'T')
      continue;

    key.assign(
      LexBorrow(G, arr_beg_chain->as_s(i)),
      arr_beg_resi->as_i(i),
      arr_beg_ins_code ? arr_beg_ins_code->as_s(i)[0] : '\0');

    sshashvalue &value = ssrecords[key];
    value.ss = ss_i;
    value.end.assign(
        LexBorrow(G, arr_end_chain->as_s(i)),
        arr_end_resi->as_i(i),
        arr_end_ins_code ? arr_end_ins_code->as_s(i)[0] : '\0');
  }

  return true;
}

/*
 * Read secondary structure
 */
static bool read_ss(PyMOLGlobals * G, cif_data * datablock,
    AtomInfoType * atInfo, CifContentInfo &info)
{
  sshashmap ssrecords;

  read_ss_(G, datablock, 'H', ssrecords, info);
  read_ss_(G, datablock, 'S', ssrecords, info);

  if (ssrecords.empty())
    return false;

  AtomInfoType *aj, *ai, *atoms_end = atInfo + VLAGetSize(atInfo);
  sshashkey key;

  for (ai = atInfo; ai < atoms_end;) {
    // advance to the next residue
    aj = ai;
    while (++ai < atoms_end &&
        AtomInfoSameResidue(G, aj, ai)) {}

    // check if residue is the beginning of a secondary structure element
    key.assign(aj->chain, aj->resv, aj->inscode);
    sshashmap::iterator it = ssrecords.find(key);

    if (it == ssrecords.end())
      continue;

    sshashvalue &value = it->second;

    // assign ss type to all atoms in the segment
    bool hit_end_residue = false;
    for (; aj < atoms_end; aj++) {
      if (value.end.compare(*aj) == 0) {
        hit_end_residue = true;
      } else if (hit_end_residue) {
        break;
      }
      aj->ssType[0] = value.ss;
    }
  }

  return true;
}

/*
 * Read the SCALEn matrix into 4x4 `matrix`
 */
static bool read_atom_site_fract_transf(PyMOLGlobals * G, const cif_data * data, float * matrix) {
  const cif_array *arr_transf[12];

  if (!(arr_transf[0] = data->get_arr("_atom_sites.fract_transf_matrix[1][1]", "_atom_sites_fract_tran_matrix_11")))
    return false;

  arr_transf[1]  = data->get_opt("_atom_sites.fract_transf_matrix[1][2]", "_atom_sites_fract_tran_matrix_12");
  arr_transf[2]  = data->get_opt("_atom_sites.fract_transf_matrix[1][3]", "_atom_sites_fract_tran_matrix_13");
  arr_transf[3]  = data->get_opt("_atom_sites.fract_transf_vector[1]", "_atom_sites_fract_tran_vector_1");
  arr_transf[4]  = data->get_opt("_atom_sites.fract_transf_matrix[2][1]", "_atom_sites_fract_tran_matrix_21");
  arr_transf[5]  = data->get_opt("_atom_sites.fract_transf_matrix[2][2]", "_atom_sites_fract_tran_matrix_22");
  arr_transf[6]  = data->get_opt("_atom_sites.fract_transf_matrix[2][3]", "_atom_sites_fract_tran_matrix_23");
  arr_transf[7]  = data->get_opt("_atom_sites.fract_transf_vector[2]", "_atom_sites_fract_tran_vector_2");
  arr_transf[8]  = data->get_opt("_atom_sites.fract_transf_matrix[3][1]", "_atom_sites_fract_tran_matrix_31");
  arr_transf[9]  = data->get_opt("_atom_sites.fract_transf_matrix[3][2]", "_atom_sites_fract_tran_matrix_32");
  arr_transf[10] = data->get_opt("_atom_sites.fract_transf_matrix[3][3]", "_atom_sites_fract_tran_matrix_33");
  arr_transf[11] = data->get_opt("_atom_sites.fract_transf_vector[3]", "_atom_sites_fract_tran_vector_3");

  for (int i = 0; i < 12; ++i)
    matrix[i] = arr_transf[i]->as_d(0);

  zero3f(matrix + 12);
  matrix[15] = 1.f;

  return true;
}

/*
 * Read anisotropic temperature factors from ATOM_SITE or ATOM_SITE_ANISOTROP
 */
static bool read_atom_site_aniso(PyMOLGlobals * G, cif_data * data,
    AtomInfoType * atInfo) {

  const cif_array *arr_label, *arr_u11, *arr_u22, *arr_u33, *arr_u12, *arr_u13, *arr_u23;
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
      std::string key(LexStr(G, ai->name));
      name_dict[key] = ai;
    }
  }

  // read aniso table
  for (int i = 0; i < arr_u11->get_nrows(); i++) {
    ai = NULL;

    if (mmcif) {
      find1(id_dict, ai, arr_label->as_i(i));
    } else {
      find1(name_dict, ai, arr_label->as_s(i));
    }

    if (!ai) {
      // expected for multi-models
      continue;
    }

    float * anisou = ai->get_anisou();
    anisou[0] = arr_u11->as_d(i) * factor;
    anisou[1] = arr_u22->as_d(i) * factor;
    anisou[2] = arr_u33->as_d(i) * factor;
    anisou[3] = arr_u12->as_d(i) * factor;
    anisou[4] = arr_u13->as_d(i) * factor;
    anisou[5] = arr_u23->as_d(i) * factor;
  }

  return true;
}

/*
 * Read GEOM_BOND
 *
 * return: BondType VLA
 */
static BondType * read_geom_bond(PyMOLGlobals * G, cif_data * data,
    AtomInfoType * atInfo) {

  const cif_array *arr_ID_1, *arr_ID_2;
  if ((arr_ID_1 = data->get_arr("_geom_bond.atom_site_id_1",
                                "_geom_bond_atom_site_label_1")) == NULL ||
      (arr_ID_2 = data->get_arr("_geom_bond.atom_site_id_2",
                                "_geom_bond_atom_site_label_2")) == NULL)
    return NULL;

  const cif_array *arr_symm_1 = data->get_opt("_geom_bond?site_symmetry_1");
  const cif_array *arr_symm_2 = data->get_opt("_geom_bond?site_symmetry_2");

  int nrows = arr_ID_1->get_nrows();
  int nAtom = VLAGetSize(atInfo);
  int nBond = 0;

  BondType *bondvla, *bond;
  bondvla = bond = VLACalloc(BondType, 6 * nAtom);

  // name -> atom index
  std::map<std::string, int> name_dict;

  // build dictionary
  for (int i = 0; i < nAtom; i++) {
    std::string key(LexStr(G, atInfo[i].name));
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

    int i1, i2;
    if (find2(name_dict, i1, key1, i2, key2)) {

      nBond++;
      BondTypeInit2(bond++, i1, i2, 1);

    } else {
      PRINTFB(G, FB_Executive, FB_Details)
        " Executive-Detail: _geom_bond name lookup failed: %s %s\n",
        key1.c_str(), key2.c_str() ENDFB(G);
    }
  }

  if (nBond) {
    VLASize(bondvla, BondType, nBond);
  } else {
    VLAFreeP(bondvla);
  }

  return bondvla;
}

/*
 * Read CHEMICAL_CONN_BOND
 *
 * return: BondType VLA
 */
static BondType * read_chemical_conn_bond(PyMOLGlobals * G, cif_data * data) {

  const cif_array *arr_number, *arr_atom_1, *arr_atom_2, *arr_type;

  if ((arr_number = data->get_arr("_atom_site?chemical_conn_number")) == NULL ||
      (arr_atom_1 = data->get_arr("_chemical_conn_bond?atom_1")) == NULL ||
      (arr_atom_2 = data->get_arr("_chemical_conn_bond?atom_2")) == NULL ||
      (arr_type   = data->get_arr("_chemical_conn_bond?type")) == NULL)
    return NULL;

  int nAtom = arr_number->get_nrows();
  int nBond = arr_atom_1->get_nrows();

  BondType *bondvla, *bond;
  bondvla = bond = VLACalloc(BondType, nBond);

  // chemical_conn_number -> atom index
  std::map<int, int> number_dict;

  // build dictionary
  for (int i = 0; i < nAtom; i++) {
    number_dict[arr_number->as_i(i)] = i;
  }

  // read table
  int i1, i2;
  for (int i = 0; i < nBond; i++) {
    if (find2(number_dict,
          i1, arr_atom_1->as_i(i),
          i2, arr_atom_2->as_i(i))) {
      BondTypeInit2(bond++, i1, i2,
          bondOrderLookup(arr_type->as_s(i)));
    } else {
      PRINTFB(G, FB_Executive, FB_Details)
        " Executive-Detail: _chemical_conn_bond name lookup failed\n" ENDFB(G);
    }
  }

  return bondvla;
}

/*
 * Read bonds from STRUCT_CONN
 *
 * Output:
 *   cset->TmpBond
 *   cset->NTmpBond
 */
static bool read_struct_conn_(PyMOLGlobals * G, cif_data * data,
    AtomInfoType * atInfo, CoordSet * cset,
    CifContentInfo &info) {

  const cif_array *col_type_id = data->get_arr("_struct_conn.conn_type_id");

  if (!col_type_id)
    return false;

  const cif_array
    *col_asym_id[2] = {NULL, NULL},
    *col_comp_id[2] = {NULL, NULL},
    *col_seq_id[2] = {NULL, NULL},
    *col_atom_id[2] = {NULL, NULL},
    *col_alt_id[2] = {NULL, NULL},
    *col_ins_code[2] = {NULL, NULL},
    *col_symm[2] = {NULL, NULL};

  if (info.use_auth) {
    col_asym_id[0] = data->get_arr("_struct_conn.ptnr1_auth_asym_id");
    col_comp_id[0] = data->get_arr("_struct_conn.ptnr1_auth_comp_id");
    col_seq_id[0]  = data->get_arr("_struct_conn.ptnr1_auth_seq_id");
    col_atom_id[0] = data->get_arr("_struct_conn.ptnr1_auth_atom_id");
    col_asym_id[1] = data->get_arr("_struct_conn.ptnr2_auth_asym_id");
    col_comp_id[1] = data->get_arr("_struct_conn.ptnr2_auth_comp_id");
    col_seq_id[1]  = data->get_arr("_struct_conn.ptnr2_auth_seq_id");
    col_atom_id[1] = data->get_arr("_struct_conn.ptnr2_auth_atom_id");

    col_alt_id[0]  = data->get_arr("_struct_conn.pdbx_ptnr1_auth_alt_id");
    col_alt_id[1]  = data->get_arr("_struct_conn.pdbx_ptnr2_auth_alt_id");

    // auth only
    col_ins_code[0] = data->get_arr("_struct_conn.pdbx_ptnr1_pdb_ins_code");
    col_ins_code[1] = data->get_arr("_struct_conn.pdbx_ptnr2_pdb_ins_code");
  }

  // for assembly chain filtering
  const cif_array *col_label_asym_id[2] = {
    data->get_arr("_struct_conn.ptnr1_label_asym_id"),
    data->get_arr("_struct_conn.ptnr2_label_asym_id")
  };

  if ((!col_asym_id[0] && !(col_asym_id[0] = col_label_asym_id[0])) ||
      (!col_comp_id[0] && !(col_comp_id[0] = data->get_arr("_struct_conn.ptnr1_label_comp_id"))) ||
      (!col_seq_id[0]  && !(col_seq_id[0]  = data->get_arr("_struct_conn.ptnr1_label_seq_id"))) ||
      (!col_atom_id[0] && !(col_atom_id[0] = data->get_arr("_struct_conn.ptnr1_label_atom_id"))) ||
      (!col_asym_id[1] && !(col_asym_id[1] = col_label_asym_id[1])) ||
      (!col_comp_id[1] && !(col_comp_id[1] = data->get_arr("_struct_conn.ptnr2_label_comp_id"))) ||
      (!col_seq_id[1]  && !(col_seq_id[1]  = data->get_arr("_struct_conn.ptnr2_label_seq_id"))) ||
      (!col_atom_id[1] && !(col_atom_id[1] = data->get_arr("_struct_conn.ptnr2_label_atom_id"))))
    return false;

  if (!col_alt_id[0]) col_alt_id[0] = data->get_opt("_struct_conn.pdbx_ptnr1_label_alt_id");
  if (!col_alt_id[1]) col_alt_id[1] = data->get_opt("_struct_conn.pdbx_ptnr2_label_alt_id");

  col_symm[0]     = data->get_opt("_struct_conn.ptnr1_symmetry");
  col_symm[1]     = data->get_opt("_struct_conn.ptnr2_symmetry");

  const cif_array *col_order = data->get_opt("_struct_conn.pdbx_value_order");

  int nrows = col_type_id->get_nrows();
  int nAtom = VLAGetSize(atInfo);
  int nBond = 0;

  BondType *bond = cset->TmpBond = VLACalloc(BondType, 6 * nAtom);

  // identifiers -> coord set index
  std::map<std::string, int> name_dict;

  for (int i = 0; i < nAtom; i++) {
    int idx = cset->atmToIdx(i);
    if (idx != -1)
      name_dict[make_mm_atom_site_label(G, atInfo + i)] = idx;
  }

  for (int i = 0; i < nrows; i++) {
    const char * type_id = col_type_id->as_s(i);
    if (strncasecmp(type_id, "covale", 6) &&
        strcasecmp(type_id, "modres") &&
#ifdef _PYMOL_IP_EXTRAS
        strcasecmp(type_id, "metalc") &&
#endif
        strcasecmp(type_id, "disulf"))
      // ignore non-covalent bonds (saltbr, hydrog)
      continue;
    if (strcmp(col_symm[0]->as_s(i),
               col_symm[1]->as_s(i)))
      // don't bond to symmetry mates
      continue;

    std::string key[2];
    for (int j = 0; j < 2; j++) {
      const char * asym_id = col_asym_id[j]->as_s(i);

      if (col_label_asym_id[j] &&
          info.is_excluded_chain(col_label_asym_id[j]->as_s(i)))
        goto next_row;

      // doen't work with label_seq_id and bulk solvent
      const char * seq_id = col_seq_id[j]->as_s(i);
      if (!seq_id[0])
        goto next_row;

      key[j] = make_mm_atom_site_label(G,
          asym_id,
          col_comp_id[j]->as_s(i),
          seq_id,
          col_ins_code[j] ? col_ins_code[j]->as_s(i) : "",
          col_atom_id[j]->as_s(i),
          col_alt_id[j]->as_s(i));
    }

    int i1, i2;
    if (find2(name_dict, i1, key[0], i2, key[1])) {
      // zero-order bond for metal coordination
      int order = strcasecmp(type_id, "metalc") ? 1 : 0;

      if (order) {
        order = bondOrderLookup(col_order->as_s(i));
      }

      nBond++;
      BondTypeInit2(bond++, i1, i2, order);

    } else {
      PRINTFB(G, FB_Executive, FB_Details)
        " Executive-Detail: _struct_conn name lookup failed: %s %s\n",
        key[0].c_str(), key[1].c_str() ENDFB(G);
    }

    // label to "continue" from inner for-loop
next_row:;
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
 *
 * return: BondType VLA
 */
static BondType * read_chem_comp_bond(PyMOLGlobals * G, cif_data * data,
    AtomInfoType * atInfo) {

  const cif_array *col_ID_1, *col_ID_2, *col_comp_id;

  if ((col_ID_1    = data->get_arr("_chem_comp_bond.atom_id_1")) == NULL ||
      (col_ID_2    = data->get_arr("_chem_comp_bond.atom_id_2")) == NULL ||
      (col_comp_id = data->get_arr("_chem_comp_bond.comp_id")) == NULL)
    return NULL;

  const cif_array *col_order = data->get_opt("_chem_comp_bond.value_order");

  int nrows = col_ID_1->get_nrows();
  int nAtom = VLAGetSize(atInfo);
  int nBond = 0;

  BondType *bondvla, *bond;
  bondvla = bond = VLACalloc(BondType, 6 * nAtom);

  // name -> atom index
  std::map<std::string, int> name_dict;

  for (int i = 0; i < nAtom; i++) {
    std::string key(LexStr(G, atInfo[i].name));
    name_dict[key] = i;
  }

  for (int i = 0; i < nrows; i++) {
    std::string key1(col_ID_1->as_s(i));
    std::string key2(col_ID_2->as_s(i));
    const char * order = col_order->as_s(i);

    int i1, i2;
    if (find2(name_dict, i1, key1, i2, key2)) {
      int order_value = bondOrderLookup(order);

      nBond++;
      BondTypeInit2(bond++, i1, i2, order_value);

    } else {
      PRINTFB(G, FB_Executive, FB_Details)
        " Executive-Detail: _chem_comp_bond name lookup failed: %s %s\n",
        key1.c_str(), key2.c_str() ENDFB(G);
    }
  }

  if (nBond) {
    VLASize(bondvla, BondType, nBond);
  } else {
    VLAFreeP(bondvla);
  }

  return bondvla;
}

/*
 * Create a new (multi-state) object-molecule from datablock
 */
static ObjectMolecule *ObjectMoleculeReadCifData(PyMOLGlobals * G, cif_data * datablock, int discrete)
{
  CoordSet ** csets = NULL;
  int ncsets;
  CifContentInfo info(G, SettingGetGlobal_b(G, cSetting_cif_use_auth));
  const char * assembly_id = SettingGetGlobal_s(G, cSetting_assembly);

  if (assembly_id && assembly_id[0]) {
    if (!get_assembly_chains(G, datablock, info.chains_filter, assembly_id))
      PRINTFB(G, FB_Executive, FB_Details)
        " ExecutiveLoad-Detail: No such assembly: '%s'\n", assembly_id ENDFB(G);
  }

  // allocate ObjectMolecule
  ObjectMolecule * I = ObjectMoleculeNew(G, (discrete > 0));
  I->Obj.Color = AtomInfoUpdateAutoColor(G);

  // read coordsets from datablock
  if ((csets = read_atom_site(G, datablock, &I->AtomInfo, info, I->DiscreteFlag))) {
    // anisou
    read_atom_site_aniso(G, datablock, I->AtomInfo);

    // secondary structure
    read_ss(G, datablock, I->AtomInfo, info);

    // trace atoms
    read_pdbx_coordinate_model(G, datablock, I);

    // polymer information
    read_entity_poly(G, datablock, info);

    // missing residues
    if (!I->DiscreteFlag && !SettingGetGlobal_i(G, cSetting_retain_order)) {
#if 0
      read_pdbx_poly_seq_scheme(G, datablock, &I->AtomInfo, info);
#else
      add_missing_ca(G, I->AtomInfo, info);
#endif
    }
  } else if ((csets = read_chem_comp_atom_model(G, datablock, &I->AtomInfo))) {
    info.type = CIF_CHEM_COMP;
  } else {
    ObjectMoleculeFree(I);
    return NULL;
  }

  // get number of atoms and coordinate sets
  I->NAtom = VLAGetSize(I->AtomInfo);
  ncsets = VLAGetSize(csets);

  // initialize the new coordsets (not data, but indices, etc.)
  for (int i = 0; i < ncsets; i++) {
    if (csets[i]) {
      csets[i]->Obj = I;
      if (!csets[i]->IdxToAtm)
        csets[i]->enumIndices();
    }
  }

  // get coordinate sets into ObjectMolecule
  VLAFreeP(I->CSet);
  I->CSet = csets;
  I->NCSet = ncsets;
  I->updateAtmToIdx();

  // handle symmetry and update fractional -> cartesian
  I->Symmetry = read_symmetry(G, datablock);
  if (I->Symmetry) {
    SymmetryUpdate(I->Symmetry);

    if(I->Symmetry->Crystal) {
      float sca[16];

      CrystalUpdate(I->Symmetry->Crystal);

      if(info.fractional) {
        for (int i = 0; i < ncsets; i++) {
          if (csets[i])
            CoordSetFracToReal(csets[i], I->Symmetry->Crystal);
        }
      } else if (info.chains_filter.empty() &&
          read_atom_site_fract_transf(G, datablock, sca)) {
        // don't do this for assemblies
        for (int i = 0; i < ncsets; i++) {
          if (csets[i])
            CoordSetInsureOrthogonal(G, csets[i], sca, I->Symmetry->Crystal);
        }
      }
    }
  }

  // coord set to use for distance based bonding and for attaching TmpBond
  CoordSet * cset = VLAGetFirstNonNULL(csets);

  // create bonds
  switch (info.type) {
    case CIF_CHEM_COMP:
      I->Bond = read_chem_comp_bond(G, datablock, I->AtomInfo);
      break;
    case CIF_CORE:
      I->Bond = read_geom_bond(G, datablock, I->AtomInfo);

      if (!I->Bond)
        I->Bond = read_chemical_conn_bond(G, datablock);

      break;
    case CIF_MMCIF:
      if (cset) {
        // sort atoms internally
        ObjectMoleculeSort(I);

        // bonds from file, goes to cset->TmpBond
        read_struct_conn_(G, datablock, I->AtomInfo, cset, info);

        // macromolecular bonding
        bond_dict_t bond_dict_local;
        if (read_chem_comp_bond_dict(datablock, bond_dict_local)) {
          ObjectMoleculeConnectComponents(I, &bond_dict_local);
        } else if(SettingGetGlobal_i(G, cSetting_connect_mode) == 4) {
          // read components.cif
          ObjectMoleculeConnectComponents(I);
        }
      }
      break;
    case CIF_UNKNOWN:
      printf("coding error...\n");
  }

  // if non of the above created I->Bond, then do distance based bonding
  if (!I->Bond) {
    if (I->DiscreteFlag) {
      ObjectMoleculeConnectDiscrete(I);
    } else if (cset) {
      ObjectMoleculeConnect(I, &I->NBond, &I->Bond, I->AtomInfo, cset, true, 3);
    }

    // guess valences for distance based bonding
    if (SettingGetGlobal_b(G, cSetting_pdb_hetatm_guess_valences)) {
      ObjectMoleculeGuessValences(I, 0, NULL, NULL, false);
    }
  } else {
    if (!I->NBond)
      I->NBond = VLAGetSize(I->Bond);

    // bonds from coordset
    if (cset && cset->TmpBond && cset->NTmpBond) {
      for (int i = 0; i < cset->NTmpBond; ++i) {
        ObjectMoleculeAddBond2(I,
            cset->IdxToAtm[cset->TmpBond[i].index[0]],
            cset->IdxToAtm[cset->TmpBond[i].index[1]],
            cset->TmpBond[i].order);
      }
      VLASize(I->Bond, BondType, I->NBond);
      VLAFreeP(cset->TmpBond);
    }
  }

  // assemblies
  if (cset && !info.chains_filter.empty()) {
    PRINTFB(G, FB_Executive, FB_Details)
      " ExecutiveLoad-Detail: Creating assembly '%s'\n", assembly_id ENDFB(G);

    CoordSet **assembly_csets = read_pdbx_struct_assembly(G, datablock,
        I->AtomInfo, cset, assembly_id);

    ObjectMoleculeSetAssemblyCSets(I, assembly_csets);
  }

  // computationally intense update tasks
  SceneCountFrames(G);
  ObjectMoleculeInvalidate(I, cRepAll, cRepInvAll, -1);
  ObjectMoleculeUpdateIDNumbers(I);
  ObjectMoleculeUpdateNonbonded(I);
  ObjectMoleculeAutoDisableAtomNameWildcard(I);

  // hetatm classification if `group_PDB` record missing
  if (info.type == CIF_MMCIF && !datablock->get_arr("_atom_site.group_pdb")) {
    I->need_hetatm_classification = true;
  }

  return I;
}

/*
 * Read one or multiple object-molecules from a CIF file. If there is only one
 * or multiplex=0, then return the object-molecule. Otherwise, create each
 * object - named by its data block name - and return NULL.
 */
ObjectMolecule *ObjectMoleculeReadCifStr(PyMOLGlobals * G, ObjectMolecule * I,
                                      const char *st, int frame,
                                      int discrete, int quiet, int multiplex,
                                      int zoom)
{
  if (I) {
    PRINTFB(G, FB_ObjectMolecule, FB_Errors)
      " Error: loading mmCIF into existing object not supported, please use 'create'\n"
      "        to append to an existing object.\n" ENDFB(G);
    return NULL;
  }

  if (multiplex > 0) {
    PRINTFB(G, FB_ObjectMolecule, FB_Errors)
      " Error: loading mmCIF with multiplex=1 not supported, please use 'split_states'.\n"
      "        after loading the object." ENDFB(G);
    return NULL;
  }

  const char * filename = NULL;
#ifndef _PYMOL_NO_CXX11
  auto cif = std::make_shared<cif_file>(filename, st);
#else
  cif_file _cif_stack(filename, st);
  auto cif = &_cif_stack;
#endif

  for (auto it = cif->datablocks.begin(); it != cif->datablocks.end(); ++it) {
    ObjectMolecule * obj = ObjectMoleculeReadCifData(G, it->second, discrete);

    if (!obj) {
      PRINTFB(G, FB_ObjectMolecule, FB_Errors)
        " mmCIF-Error: no coordinates found in data_%s\n", it->first ENDFB(G);
      continue;
    }

#if !defined(_PYMOL_NOPY) && !defined(_PYMOL_NO_CXX11)
    // we only provide access from the Python API so far
    if (SettingGetGlobal_b(G, cSetting_cif_keepinmemory)) {
      obj->m_cifdata = it->second;
      obj->m_ciffile = cif;
    }
#endif

    if (cif->datablocks.size() == 1 || multiplex == 0)
      return obj;

    // multiplexing
    ObjectSetName((CObject*) obj, it->first);
    ExecutiveDelete(G, obj->Obj.Name);
    ExecutiveManageObject(G, (CObject*) obj, zoom, true);
  }

  return NULL;
}

/*
 * Bond dictionary getter, with on-demand download of residue dictionaries
 */
const bond_dict_t::mapped_type * bond_dict_t::get(PyMOLGlobals * G, const char * resn, bool try_download) {
  auto key = make_key(resn);
  auto it = find(key);

  if (it != end())
    return &it->second;

  if (unknown_resn.count(key))
    return NULL;

#ifndef _PYMOL_NOPY
  if (try_download) {
    int blocked = PAutoBlock(G);
    bool downloaded = false;

    // call into Python
    PyObject * pyfilename = PYOBJECT_CALLMETHOD(G->P_inst->cmd,
        "download_chem_comp", "siO", resn,
        !Feedback(G, FB_Executive, FB_Details),
        G->P_inst->cmd);

    if (pyfilename) {
      const char * filename = PyString_AsString(pyfilename);

      // update
      if ((downloaded = (filename && filename[0]))) {
        cif_file cif(filename);
        for (auto it = cif.datablocks.begin(); it != cif.datablocks.end(); ++it)
          read_chem_comp_bond_dict(it->second, *this);
      }

      Py_DECREF(pyfilename);
    }

    PAutoUnblock(G, blocked);

    if (downloaded) {
      // second attempt to look up, from eventually updated dictionary
      return get(G, resn, false);
    }
  }
#endif

  PRINTFB(G, FB_Executive, FB_Warnings)
    " ExecutiveLoad-Warning: No _chem_comp_bond data for residue '%s'\n", resn
    ENDFB(G);

  // don't try downloading again
  unknown_resn.insert(key);

  return NULL;
}

// vi:sw=2:ts=2:expandtab
