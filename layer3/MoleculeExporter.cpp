/*
 * Molecular file formats export
 *
 * (c) 2016 Schrodinger, Inc.
 */

#include <vector>
#include <map>
#include <algorithm>
#include <cstdarg>
#include <clocale>
#include <memory>

#ifndef _PYMOL_NO_MSGPACKC
#include <mmtf.hpp>
#include <mmtf/export_helpers.hpp>
#endif

#include "os_std.h"

#include "MoleculeExporter.h"
#include "Selector.h"
#include "SelectorDef.h"
#include "Executive.h"
#include "Vector.h"
#include "ObjectMolecule.h"
#include "CoordSet.h"
#include "Mol2Typing.h"
#include "MmodTyping.h"
#include "Color.h"
#include "Util.h"
#include "Lex.h"
#include "P.h"
#include "PConv.h"
#include "CifDataValueFormatter.h"
#include "MaeExportHelpers.h"

#ifdef _PYMOL_IP_PROPERTIES
#include "Property.h"
#endif

/*
 * Get the capitalized element symbol. Assume that ai->elem is either all
 * caps (e.g. loaded from PDB) or capitalized.
 */
class ElemCanonicalizer {
  ElemName m_buffer;
public:
  const char * operator() (const AtomInfoType * ai) {
    const char * elem = ai->elem;
    if (ai->protons < 1 || !elem[0] || !elem[1] || islower(elem[1]))
      return elem;
    m_buffer[0] = elem[0];
    UtilNCopyToLower(m_buffer + 1, elem + 1, sizeof(ElemName) - 1);
    return m_buffer;
  }
};

/*
 * "sprintf" into VLA string buffer. Will grow (and reallocate) the VLA if
 * needed.
 *
 * Returns the number of characters written (excluding the null byte).
 */
static
int VLAprintf(pymol::vla<char>& vla, int offset, const char * format, ...) {
  int n, size = vla.size() - offset;
  va_list ap;

  va_start(ap, format);
  n = vsnprintf(vla + offset, std::max(0, size), format, ap);
  va_end(ap);

#ifdef WIN32
  // Windows API: If len > count, then count characters are stored in buffer,
  // no null-terminator is appended, and a negative value is returned.
  if (n < 0) {
    va_start(ap, format);
    n = _vscprintf(format, ap);
    va_end(ap);
  }
#endif

  // C99 standard: a return value of size or more means that the output was
  // truncated (glibc since 2.1; not on Windows)
  if (n >= size) {
    vla.check(offset + n);
    va_start(ap, format);
    vsprintf(vla + offset, format, ap);
    va_end(ap);
  }

  return n;
}

// for "multisave" behavior
enum {
  cMolExportGlobal     = 0,
  cMolExportByObject   = 1,
  cMolExportByCoordSet = 2,
};

// for re-indexed bonds
struct BondRef {
  BondType * ref;
  int id1, id2;
};

// for re-indexed atoms
struct AtomRef {
  AtomInfoType * ref;
  float coord[3];
  int id;
};

/*
 * Abstract base class for exporting molecular selections
 */
struct MoleculeExporter {
  pymol::vla<char> m_buffer; // out

protected:
  int m_offset;
  CoordSet        * m_last_cs;
  ObjectMolecule  * m_last_obj;
  int m_last_state;

  PyMOLGlobals * G;
  SeleCoordIterator m_iter;

  bool m_retain_ids;
  int m_id;

  struct matrix_t {
    double storage[16];
    const double * ptr;
  } m_mat_ref,  // inverse of reference object transformation
    m_mat_full, // for ANISOU
    m_mat_move; // for coordinates (TTT)

  // atom coordinate
private:
  float m_coord_tmp[3];
protected:
  const float *m_coord;

  // how to restart models
  int m_multi;

  std::vector<BondRef> m_bonds;
  std::vector<int> m_tmpids;

public:
  // quasi constructor (easier to inherit than a real constructor)
  virtual void init(PyMOLGlobals * G_) {
    G = G_;

    m_buffer.resize(1280);
    m_buffer[0] = '\0';

    m_mat_ref.ptr = nullptr;
    m_offset = 0;
    m_last_cs  = nullptr;
    m_last_obj = nullptr;
    m_last_state = -1;
    m_retain_ids = false;
    m_id = 0;

    setMulti(getMultiDefault());
  }

  // destructor
  virtual ~MoleculeExporter() = default;

  /*
   * Do the export (e.g. populate "m_buffer" with file contents)
   */
  void execute(int sele, int state);

  /*
   * how to handle selections which span multiple objects
   */
  void setMulti(int multi) {
    if (multi != -1)
      m_multi = multi;
  }

private:
  /*
   * Reset the "index" fields in the selector table
   */
  void resetTmpIDs() {
    m_tmpids.resize(m_iter.obj->NAtom, 0);
    std::fill(m_tmpids.begin(), m_tmpids.end(), 0);
  }

protected:
  /*
   * Get the current atom's running index (1-based). Most formats use
   * this index to reference atoms from the bonds table. If "retain_ids"
   * is true, then this is AtomInfoType::id.
   */
  int getTmpID() {
    return m_tmpids[m_iter.getAtm()];
  }

  /*
   * Get the running index for the given atom.
   */
  int getTmpID(int atm) {
    return m_tmpids[atm];
  }

  /*
   * Get the current state title if not empty, or the object name otherwise.
   */
  const char * getTitleOrName() {
    if (!m_iter.cs)
      return "untitled";
    return m_iter.cs->Name[0] ? m_iter.cs->Name : m_iter.obj->Obj.Name;
  }

public:
  /*
   * Set the reference coordinate frame to the inverse of the given
   * object's transformation matrix.
   */
  void setRefObject(const char * ref_object, int ref_state);

private:
  /*
   * Update a transformation matrix for the current state
   */
  void updateMatrix(matrix_t& matrix, bool history);

  /*
   * Populates the "m_bonds" vector
   */
  void populateBondRefs();

protected:
  // functions to be implemented by derived classes
  virtual int getMultiDefault() const = 0;
  virtual bool isExcludedBond(int atm1, int atm2);
  virtual bool isExcludedBond(const BondType * bond);
  virtual void writeAtom() = 0;
  virtual void writeBonds() = 0;
  virtual void beginObject();
  virtual void beginCoordSet();
  virtual void endObject();
  virtual void endCoordSet();
  virtual void beginMolecule() {}
  virtual void beginFile() {}
};

/*
 * Return true if this bond should not be exported
 */
bool MoleculeExporter::isExcludedBond(int atm1, int atm2) {
  return false;
}
bool MoleculeExporter::isExcludedBond(const BondType * bond) {
  return isExcludedBond(bond->index[0], bond->index[1]);
}

void MoleculeExporter::beginObject() {
  if (m_multi != cMolExportByCoordSet) {
    resetTmpIDs();
  }

  if (m_multi == cMolExportByObject) {
    beginMolecule();
  }
}

void MoleculeExporter::beginCoordSet() {
  if (m_multi == cMolExportByCoordSet) {
    resetTmpIDs();

    beginMolecule();
  }
}

void MoleculeExporter::endObject() {
  if (m_multi != cMolExportByCoordSet) {
    populateBondRefs();
  }

  if (m_multi == cMolExportByObject) {
    writeBonds();
    m_id = 0;
  }
}

void MoleculeExporter::endCoordSet() {
  if (m_multi == cMolExportByCoordSet) {
    populateBondRefs();
    writeBonds();
    m_id = 0;
  }
}

void MoleculeExporter::execute(int sele, int state) {
  m_iter.init(G, sele, state);
  m_iter.setPerObject(m_multi != cMolExportGlobal);

  beginFile();

  while (m_iter.next()) {
    if (m_last_cs != m_iter.cs) {
      if (m_last_cs) {
        endCoordSet();
      } else if (m_multi == cMolExportGlobal) {
        beginMolecule();
      }

      if (m_last_obj != m_iter.obj) {
        if (m_last_obj)
          endObject();

        beginObject();
        m_last_obj = m_iter.obj;
      }

      // update transformation matrices
      updateMatrix(m_mat_full, true);
      updateMatrix(m_mat_move, false);

      beginCoordSet();
      m_last_cs = m_iter.cs;
    }

    // for bonds
    if (!m_tmpids[m_iter.getAtm()]) {
      m_id = m_retain_ids ? m_iter.getAtomInfo()->id : (m_id + 1);
      m_tmpids[m_iter.getAtm()] = m_id;
    }

    // atom coordinate
    m_coord = m_iter.getCoord();
    if (m_mat_move.ptr) {
      transform44d3f(m_mat_move.ptr, m_coord, m_coord_tmp);
      m_coord = m_coord_tmp;
    }

    writeAtom();
  }

  if (m_last_cs)
    endCoordSet();
  if (m_last_obj)
    endObject();
  else if (m_multi == cMolExportGlobal)
    // empty selection
    beginMolecule();

  if (m_multi == cMolExportGlobal) {
    writeBonds();
  }

  m_buffer.resize(m_offset);
}

void MoleculeExporter::setRefObject(const char * ref_object, int ref_state) {
  double matrix[16];

  m_mat_ref.ptr = nullptr;

  if (!ref_object || !ref_object[0])
    return;

  auto base = ExecutiveFindObjectByName(G, ref_object);
  if (!base)
    return;

  if(ref_state < 0) {
    ref_state = ObjectGetCurrentState(base, true);
  }

  if(ObjectGetTotalMatrix(base, ref_state, true, matrix)) {
    invert_special44d44d(matrix, m_mat_ref.storage);
    m_mat_ref.ptr = m_mat_ref.storage;
  }
}

void MoleculeExporter::updateMatrix(matrix_t& matrix, bool history) {
  const auto& ref = m_mat_ref.ptr;
  if (ObjectGetTotalMatrix(reinterpret_cast<CObject*>(m_iter.obj),
        m_iter.state, history, matrix.storage)) {
    if (ref) {
      left_multiply44d44d(ref, matrix.storage);
    }
    matrix.ptr = matrix.storage;
  } else {
    matrix.ptr = ref;
  }
}

void MoleculeExporter::populateBondRefs() {
  auto& obj = m_last_obj;
  int id1, id2;

  for (auto bond = obj->Bond, bond_end = obj->Bond + obj->NBond;
      bond != bond_end; ++bond) {

    auto atm1 = bond->index[0];
    auto atm2 = bond->index[1];

    if (!(id1 = getTmpID(atm1)) ||
        !(id2 = getTmpID(atm2)))
      continue;

    if (isExcludedBond(bond))
      continue;

    if (id1 > id2)
      std::swap(id1, id2);

    // emit bond
    m_bonds.emplace_back(BondRef { bond, id1, id2 });
  }
}

// ---------------------------------------------------------------------------------- //

struct MoleculeExporterPDB : public MoleculeExporter {
  bool m_conect_all;
  bool m_conect_nodup;
  bool m_mdl_written;
  bool m_use_ter_records;
  const AtomInfoType * m_pre_ter = nullptr;
  PDBInfoRec m_pdb_info;

  // quasi constructor
  void init(PyMOLGlobals * G_) override {
    MoleculeExporter::init(G_);

    UtilZeroMem((void *) &m_pdb_info, sizeof(PDBInfoRec));

    m_conect_all    = false;
    m_mdl_written   = false;
    m_conect_nodup  = SettingGetGlobal_b(G, cSetting_pdb_conect_nodup);
    m_retain_ids    = SettingGetGlobal_b(G, cSetting_pdb_retain_ids);
    m_use_ter_records = SettingGetGlobal_b(G, cSetting_pdb_use_ter_records);
  }

  int getMultiDefault() const override {
    // single entry format by default, but we also support writing multiple
    // entries by writing a HEADER record for every object (1) or state (2)
    return cMolExportGlobal;
  }

  void writeENDMDL() {
    if (m_mdl_written) {
      m_offset += VLAprintf(m_buffer, m_offset, "ENDMDL\n");
      m_mdl_written = false;
    }
  }

  void writeEND() {
    if (!SettingGetGlobal_b(G, cSetting_pdb_no_end_record)) {
      m_offset += VLAprintf(m_buffer, m_offset, "END\n");
    }
  }

  /*
   * Write a TER record if the previous atom was polymer and `ai`
   * is NULL, non-polymer, or has a different chain identifier.
   */
  void writeTER(const AtomInfoType* ai) {
    if (!m_use_ter_records)
      return;

    if (ai && !(ai->flags & cAtomFlag_polymer)) {
      ai = nullptr;
    }

    if (m_pre_ter && !(ai && ai->chain == m_pre_ter->chain)) {
      m_offset += VLAprintf(m_buffer, m_offset, "TER   \n");
    }

    m_pre_ter = ai;
  }

  void writeAtom() override {
    writeTER(m_iter.getAtomInfo());

    CoordSetAtomToPDBStrVLA(G, &m_buffer, &m_offset, m_iter.getAtomInfo(),
        m_coord, getTmpID() - 1, &m_pdb_info, m_mat_full.ptr);
  }

  void writeBonds() override {
    writeENDMDL();

    std::map<int, std::vector<int>> conect;

    for (auto& bond : m_bonds) {
      int order = m_conect_nodup ? 1 : bond.ref->order;
      for (int i = 0; i < 2; ++i) {
        for (int d = 0; d < order; ++d) {
          conect[bond.id1].push_back(bond.id2);
        }
        std::swap(bond.id1, bond.id2);
      }
    }

    m_bonds.clear();

    for (auto& rec : conect) {
      for (int i = 0, i_end = rec.second.size(); i != i_end;) {
        m_offset += VLAprintf(m_buffer, m_offset, "CONECT%5d", rec.first);
        // up to 4 bonds per record
        for (int j = std::min(i + 4, i_end); i != j; ++i) {
          m_offset += VLAprintf(m_buffer, m_offset, "%5d", rec.second[i]);
        }
        m_offset += VLAprintf(m_buffer, m_offset, "\n");
      }
    }

    writeEND();
  }

  void writeCryst1() {
    const auto& sym = m_iter.cs->Symmetry ? m_iter.cs->Symmetry : m_iter.obj->Symmetry;

    if (sym && sym->Crystal) {
      const auto& dim   = sym->Crystal->Dim;
      const auto& angle = sym->Crystal->Angle;
      m_offset += VLAprintf(m_buffer, m_offset,
          "CRYST1%9.3f%9.3f%9.3f%7.2f%7.2f%7.2f %-11s%4d\n",
          dim[0], dim[1], dim[2], angle[0], angle[1], angle[2],
          sym->SpaceGroup, sym->PDBZValue);
    }
  }

  void beginObject() override {
    MoleculeExporter::beginObject();

    m_conect_all = SettingGet_b(G, m_iter.obj->Obj.Setting, nullptr, cSetting_pdb_conect_all);

    if (m_multi == cMolExportByObject) {
      m_offset += VLAprintf(m_buffer, m_offset, "HEADER    %.40s\n", m_iter.obj->Obj.Name);
      writeCryst1();
    }
  }

  void beginCoordSet() override {
    MoleculeExporter::beginCoordSet();

    if (m_multi == cMolExportByCoordSet) {
      m_offset += VLAprintf(m_buffer, m_offset, "HEADER    %.40s\n", getTitleOrName());
      writeCryst1();
    }

    if (m_iter.isMultistate()
        && (m_iter.isPerObject() || m_iter.state != m_last_state)) {
      m_offset += VLAprintf(m_buffer, m_offset, "MODEL     %4d\n", m_iter.state + 1);
      m_last_state = m_iter.state;
      m_mdl_written = true;
    }
  }

  void endCoordSet() override {
    writeTER(nullptr);

    MoleculeExporter::endCoordSet();

    if (m_iter.isPerObject() || m_iter.state != m_last_state) {
      writeENDMDL();
    }
  }

  bool isExcludedBond(int atm1, int atm2) override {
    const auto& atInfo = m_last_obj->AtomInfo;
    return !(m_conect_all || atInfo[atm1].hetatm || atInfo[atm2].hetatm);
  }
};

// ---------------------------------------------------------------------------------- //

struct MoleculeExporterPQR: public MoleculeExporterPDB {
  // quasi constructor
  void init(PyMOLGlobals * G_) override {
    MoleculeExporterPDB::init(G_);

    m_pdb_info.variant = PDB_VARIANT_PQR;
    m_pdb_info.pqr_workarounds = SettingGetGlobal_b(G, cSetting_pqr_workarounds);
  }

  bool isExcludedBond(int atm1, int atm2) override {
    // no bonds for PQR format
    return true;
  }
};

// ---------------------------------------------------------------------------------- //

struct MoleculeExporterCIF : public MoleculeExporter {
  const char * m_molecule_name;
  CifDataValueFormatter cifrepr;

  // quasi constructor
  void init(PyMOLGlobals * G_) override {
    MoleculeExporter::init(G_);

    // The formatter uses a circular buffer of the given size. The size
    // must be at least equal to the maximum number of invocations of `cifrepr`
    // in any member function (count two for each call with type `char`).
    cifrepr.m_buf.resize(10);

    m_retain_ids    = SettingGetGlobal_b(G, cSetting_pdb_retain_ids);
    m_molecule_name = "multi";

    m_offset += VLAprintf(m_buffer, m_offset,
        "# generated by PyMOL " _PyMOL_VERSION "\n");
  }

  int getMultiDefault() const override {
    return cMolExportByObject;
  }

  void writeCellSymmetry() {
    const auto& sym = m_iter.cs->Symmetry ? m_iter.cs->Symmetry : m_iter.obj->Symmetry;

    if (sym && sym->Crystal) {
      const auto& dim   = sym->Crystal->Dim;
      const auto& angle = sym->Crystal->Angle;
      m_offset += VLAprintf(m_buffer, m_offset, "#\n"
          "_cell.entry_id %s\n"
          "_cell.length_a %.3f\n"
          "_cell.length_b %.3f\n"
          "_cell.length_c %.3f\n"
          "_cell.angle_alpha %.2f\n"
          "_cell.angle_beta  %.2f\n"
          "_cell.angle_gamma %.2f\n"
          "_symmetry.entry_id %s\n"
          "_symmetry.space_group_name_H-M %s\n",
          cifrepr(m_molecule_name),
          dim[0], dim[1], dim[2], angle[0], angle[1], angle[2],
          cifrepr(m_molecule_name),
          cifrepr(sym->SpaceGroup));
    }
  }

  void beginMolecule() override {
    MoleculeExporter::beginMolecule();

    switch (m_multi) {
      case cMolExportByObject:   m_molecule_name = m_iter.obj->Obj.Name; break;
      case cMolExportByCoordSet: m_molecule_name = getTitleOrName(); break;
    }

    // data header
    m_offset += VLAprintf(m_buffer, m_offset, "#\n"
        "data_%s\n_entry.id %s\n", m_molecule_name,
        cifrepr(m_molecule_name));

    // symmetry
    writeCellSymmetry();

    // atom table header
    m_offset += VLAprintf(m_buffer, m_offset, "#\n"
        "loop_\n"
        "_atom_site.group_PDB\n"
        "_atom_site.id\n"
        "_atom_site.type_symbol\n"
        "_atom_site.label_atom_id\n"
        "_atom_site.label_alt_id\n"
        "_atom_site.label_comp_id\n"
        "_atom_site.label_asym_id\n"
        "_atom_site.label_entity_id\n"
        "_atom_site.label_seq_id\n"
        "_atom_site.pdbx_PDB_ins_code\n"
        "_atom_site.Cartn_x\n"
        "_atom_site.Cartn_y\n"
        "_atom_site.Cartn_z\n"
        "_atom_site.occupancy\n"
        "_atom_site.B_iso_or_equiv\n"
        "_atom_site.pdbx_formal_charge\n"
        "_atom_site.auth_asym_id\n"
        "_atom_site.pdbx_PDB_model_num\n");
  }

  void writeAtom() override {
    const AtomInfoType * ai = m_iter.getAtomInfo();
    const char * entity_id = nullptr;

#ifdef _PYMOL_IP_PROPERTIES
    char entity_id_buf[16];
    if (ai->prop_id) {
      entity_id = PropertyGetAsString(G, ai->prop_id, "entity_id", entity_id_buf);
    }
#endif

    if (!entity_id) {
      entity_id = LexStr(G, ai->custom);
    }

    m_offset += VLAprintf(m_buffer, m_offset,
        "%-6s %-3d %s %-3s " // type .. name
        "%s %-3s %s %s " // alt .. entity_id
        "%d %s %6.3f %6.3f %6.3f " // resv .. z
        "%4.2f %6.2f %d %s %d\n",  // q .. state
        ai->hetatm ? "HETATM" : "ATOM",
        getTmpID(),
        cifrepr(ai->elem),
        cifrepr(LexStr(G, ai->name)),
        cifrepr(ai->alt),
        cifrepr(LexStr(G, ai->resn)),
        cifrepr(LexStr(G, ai->segi)),
        cifrepr(entity_id),
        ai->resv,
        cifrepr(ai->inscode, "?"),
        m_coord[0], m_coord[1], m_coord[2],
        ai->q, ai->b, ai->formalCharge,
        cifrepr(LexStr(G, ai->chain)),
        m_iter.state + 1);
  }

  void writeBonds() override {
    // TODO
    // Bond categories:
    // 1) within residue -> _chem_comp_bond
    // 2) polymer links -> implicit
    // 3) others -> _struct_conn
    m_bonds.clear();
  }

  bool isExcludedBond(int atm1, int atm2) override {
    // TODO
    // polymer links
    return true;
  }
};

// ---------------------------------------------------------------------------------- //

/*
 * Experimental PyMOL style annotated mmCIF. Exports colors, representation,
 * and secondary structure.
 */
struct MoleculeExporterPMCIF : public MoleculeExporterCIF {
  void beginMolecule() override {
    MoleculeExporterCIF::beginMolecule();

    m_offset += VLAprintf(m_buffer, m_offset, "#\n"
        "_atom_site.pymol_color\n"
        "_atom_site.pymol_reps\n"
        "_atom_site.pymol_ss\n");
  }

  void writeAtom() override {
    MoleculeExporterCIF::writeAtom();

    const AtomInfoType * ai = m_iter.getAtomInfo();

    m_offset += VLAprintf(m_buffer, m_offset,
        "%d %d %s\n",
        ai->color,
        ai->visRep,
        cifrepr(ai->ssType));
  }

  void writeBonds() override {
    if (m_bonds.empty())
      return;

    m_offset += VLAprintf(m_buffer, m_offset, "#\n"
        "loop_\n"
        "_pymol_bond.atom_site_id_1\n"
        "_pymol_bond.atom_site_id_2\n"
        "_pymol_bond.order\n");

    for (auto& bond : m_bonds) {
      m_offset += VLAprintf(m_buffer, m_offset, "%d %d %d\n",
          bond.id1, bond.id2, bond.ref->order);
    }

    m_bonds.clear();
  }

  bool isExcludedBond(int atm1, int atm2) override {
    return false;
  }
};

// ---------------------------------------------------------------------------------- //

struct MoleculeExporterMOL : public MoleculeExporter {
  int m_chiral_flag;
  std::vector<AtomRef> m_atoms;
  ElemCanonicalizer elemGetter;

  int getMultiDefault() const override {
    // single entry format
    return cMolExportGlobal;
  }

  void writeAtom() override {
    const auto ai = m_iter.getAtomInfo();

    if (ai->stereo)
      m_chiral_flag = 1;

    m_atoms.emplace_back(AtomRef { ai, { m_coord[0], m_coord[1], m_coord[2] }, getTmpID() });
  }

  void writeCTabV3000() {
    m_offset += VLAprintf(m_buffer, m_offset,
        "  0  0  0  0  0  0  0  0  0  0999 V3000\n"
        "M  V30 BEGIN CTAB\n"
        "M  V30 COUNTS %d %d 0 0 %d\n"
        "M  V30 BEGIN ATOM\n",
        m_atoms.size(), m_bonds.size(), m_chiral_flag);

    // write atoms
    for (auto& atom : m_atoms) {
      auto ai = atom.ref;

      m_offset += VLAprintf(m_buffer, m_offset, "M  V30 %d %s %.4f %.4f %.4f 0",
          atom.id, elemGetter(ai), atom.coord[0], atom.coord[1], atom.coord[2]);

      if (ai->formalCharge)
        m_offset += VLAprintf(m_buffer, m_offset, " CHG=%d", ai->formalCharge);

      if (ai->stereo)
        m_offset += VLAprintf(m_buffer, m_offset, " CFG=%d", ai->stereo);

      m_offset += VLAprintf(m_buffer, m_offset, "\n");
    }

    m_atoms.clear();

    m_offset += VLAprintf(m_buffer, m_offset,
        "M  V30 END ATOM\n"
        "M  V30 BEGIN BOND\n");

    // write bonds
    int n_bonds = 0;
    for (auto& bond : m_bonds) {
      m_offset += VLAprintf(m_buffer, m_offset, "M  V30 %d %d %d %d\n",
          ++n_bonds, bond.ref->order, bond.id1, bond.id2);
    }

    m_bonds.clear();

    // end
    m_offset += VLAprintf(m_buffer, m_offset,
        "M  V30 END BOND\n"
        "M  V30 END CTAB\n"
        "M  END\n");
  }

  void writeCTabV2000() {
    // counts line
    m_offset += VLAprintf(m_buffer, m_offset,
        "%3d%3d  0  0%3d  0  0  0  0  0999 V2000\n",
        (int) m_atoms.size(), (int) m_bonds.size(), m_chiral_flag);

    // write atoms
    for (auto& atom : m_atoms) {
      auto ai = atom.ref;
      int chg = ai->formalCharge;
      m_offset += VLAprintf(m_buffer, m_offset,
          "%10.4f%10.4f%10.4f %-3s 0  %1d  %1d  0  0  0  0  0  0  0  0  0\n",
          atom.coord[0], atom.coord[1], atom.coord[2],
          elemGetter(ai), chg ? (4 - chg) : 0, (int) ai->stereo);
    }

    m_atoms.clear();

    // write bonds
    for (auto& bond : m_bonds) {
      m_offset += VLAprintf(m_buffer, m_offset, "%3d%3d%3d%3d  0  0  0\n",
          bond.id1, bond.id2, bond.ref->order, (int) bond.ref->stereo);
    }

    m_bonds.clear();

    // end
    m_offset += VLAprintf(m_buffer, m_offset, "M  END\n");
  }

  void writeBonds() override {
    if (m_atoms.size() > 999 || m_bonds.size() > 999) {
      PRINTFB(G, FB_ObjectMolecule, FB_Warnings)
        " Warning: MOL/SDF 999 atom/bond limit reached, using V3000\n" ENDFB(G);
      writeCTabV3000();
    } else {
      writeCTabV2000();
    }
  }

  void beginMolecule() override {
    MoleculeExporter::beginMolecule();

    m_offset += VLAprintf(m_buffer, m_offset,
        "%s\n  PyMOL%3.3s          3D                             0\n\n",
        getTitleOrName(), _PyMOL_VERSION);

    m_chiral_flag = 0;
  }

  bool isExcludedBond(const BondType * bond) override {
    // MOL format doesn't know zero order bonds. Writing them as order "0"
    // will produce a nonstandard file which may be rejected by other
    // applications (e.g. RDKit).
    if (bond->order == 0) {
      return !SettingGetGlobal_b(G, cSetting_sdf_write_zero_order_bonds);
    }

    return false;
  }
};

// ---------------------------------------------------------------------------------- //

struct MoleculeExporterSDF : public MoleculeExporterMOL {

  int getMultiDefault() const override {
    // multi-entry format
    return cMolExportByCoordSet;
  }

  void writeProperties() {
#ifdef _PYMOL_IP_PROPERTIES
#endif
  }

  void writeBonds() override {
    MoleculeExporterMOL::writeBonds();

    writeProperties();
    m_offset += VLAprintf(m_buffer, m_offset, "$$$$\n");
  }
};

// ---------------------------------------------------------------------------------- //

struct MOL2_SubSt {
  AtomInfoType * ai;
  int root_id;
  const char * resn;
};

static const char MOL2_bondTypes[][3] = {
  "nc", "1", "2", "3", "ar"
};

struct MoleculeExporterMOL2 : public MoleculeExporter {
  int m_n_atoms; // atom count
  int m_counts_offset; // offset for deferred counts writing
  std::vector<MOL2_SubSt> m_substs; // substructures

  int getMultiDefault() const override {
    // multi-entry format
    return cMolExportByCoordSet;
  }

  void beginFile() override {
    m_offset += VLAprintf(m_buffer, m_offset,
        "# created with PyMOL " _PyMOL_VERSION "\n");
  }

  void beginMolecule() override {
    MoleculeExporter::beginMolecule();

    // RTI MOLECULE
    m_offset += VLAprintf(m_buffer, m_offset, "@<TRIPOS>MOLECULE\n"
        "%s\n", getTitleOrName());

    // defer until number of substructures known
    m_counts_offset = m_offset;
    m_offset += VLAprintf(m_buffer, m_offset,
        "X X X                   \n" // deferred
        "SMALL\n"
        "USER_CHARGES\n"
        "@<TRIPOS>ATOM\n");

    m_n_atoms = 0;
  }

  void beginObject() override {
    MoleculeExporter::beginObject();
    ObjectMoleculeVerifyChemistry(m_iter.obj, m_iter.state);
  }

  void writeAtom() override {
    const auto ai = m_iter.getAtomInfo();

    if (m_substs.empty() ||
        !AtomInfoSameResidue(G, ai, m_substs.back().ai)) {
      m_substs.emplace_back(MOL2_SubSt { ai, getTmpID(),
          ai->resn ? LexStr(G, ai->resn) : "UNK" });
    }

    // RTI ATOM
    // atom_id atom_name x y z atom_type [subst_id
    //   [subst_name [charge [status_bit]]]]
    m_offset += VLAprintf(m_buffer, m_offset,
        "%d\t%4s\t%.3f\t%.3f\t%.3f\t%2s\t%d\t%s%d%.1s\t%.3f\t%s\n",
        getTmpID(),
        ai->name ? LexStr(G, ai->name) : ai->elem[0] ? ai->elem : "X",
        m_coord[0], m_coord[1], m_coord[2],
        getMOL2Type(m_iter.obj, m_iter.getAtm()),
        m_substs.size(),
        m_substs.back().resn, ai->resv, &ai->inscode, // subst_name
        ai->partialCharge,
        (ai->flags & cAtomFlag_solvent) ? "WATER" : "");

    ++m_n_atoms;
  }

  void writeBonds() override {
    // atom count
    m_counts_offset += sprintf(m_buffer + m_counts_offset, "%d %d %d",
        m_n_atoms, (int) m_bonds.size(), (int) m_substs.size());
    m_buffer[m_counts_offset] = ' '; // overwrite terminator

    // RTI BOND
    // bond_id origin_atom_id target_atom_id bond_type [status_bits]

    m_offset += VLAprintf(m_buffer, m_offset, "@<TRIPOS>BOND\n");

    int bond_id = 0;
    for (auto& bond : m_bonds) {
      m_offset += VLAprintf(m_buffer, m_offset, "%d %d %d %s\n",
          ++bond_id,
          bond.id1,
          bond.id2,
          MOL2_bondTypes[bond.ref->order % 5]);
    }

    m_bonds.clear();

    // RTI SUBSTRUCTURE
    // subst_id subst_name root_atom [subst_type [dict_type
    //   [chain [sub_type [inter_bonds [status [comment]]]]]]]

    m_offset += VLAprintf(m_buffer, m_offset, "@<TRIPOS>SUBSTRUCTURE\n");

    int subst_id = 0;
    for (auto& subst : m_substs) {
      const auto& ai = subst.ai;
      m_offset += VLAprintf(m_buffer, m_offset, "%d\t%s%d%.1s\t%d\t%s\t1 %s\t%s\n",
          ++subst_id,
          subst.resn, ai->resv, &ai->inscode,     // subst_name
          subst.root_id,                          // root_atom
          (ai->flags & cAtomFlag_polymer) ? "RESIDUE" : "GROUP",
          ai->chain ? LexStr(G, ai->chain) : ai->segi ? LexStr(G, ai->segi) : "****",
          subst.resn);
    }

    m_substs.clear();
  }

};

// ---------------------------------------------------------------------------------- //

struct MoleculeExporterMAE : public MoleculeExporter {
  int m_n_atoms;
  int m_n_atoms_offset;
  int m_n_arom_bonds;
  std::map<int, const AtomInfoType *> m_atoms;

  // quasi constructor
  void init(PyMOLGlobals * G_) override {
    MoleculeExporter::init(G_);

    m_n_arom_bonds = 0;
  }

  int getMultiDefault() const override {
    // multi-entry format
    return cMolExportByCoordSet;
  }

  void beginFile() override {
    m_offset += VLAprintf(m_buffer, m_offset,
        "{ s_m_m2io_version ::: 2.0.0 }\n"
        "# created with PyMOL " _PyMOL_VERSION " #\n");
  }

  void beginMolecule() override {
    MoleculeExporter::beginMolecule();

    std::string groupid = MaeExportGetSubGroupId(G, reinterpret_cast<CObject*>(m_iter.obj));

    m_offset += VLAprintf(m_buffer, m_offset,
        "\nf_m_ct {\n"
        "s_m_subgroupid\n"
        "s_m_title\n"
        ":::\n"
        "\"%s\"\n"
        "\"%s\"\n",
        groupid.c_str(),
        getTitleOrName());

    // defer until number of atoms known
    m_n_atoms_offset = m_offset;

    m_offset += VLAprintf(m_buffer, m_offset,
        "m_atom[X]            {\n" // place holder
        "# First column is atom index #\n"
        "i_m_mmod_type\n"
        "r_m_x_coord\n"
        "r_m_y_coord\n"
        "r_m_z_coord\n"
        "i_m_residue_number\n"
        "s_m_insertion_code\n"
        "s_m_chain_name\n"
        "s_m_pdb_residue_name\n"
        "s_m_pdb_atom_name\n"
        "i_m_atomic_number\n"
        "i_m_formal_charge\n"
        "s_m_color_rgb\n"
        "i_m_secondary_structure\n"
        "r_m_pdb_occupancy\n"
        "i_pdb_PDB_serial\n"
        "i_m_visibility\n"
        "i_m_representation\n"
        "i_m_ribbon_style\n"
        "i_m_ribbon_color\n"
        "s_m_ribbon_color_rgb\n"
        "s_m_label_format\n"
        "i_m_label_color\n"
        "s_m_label_user_text\n"
        ":::\n");

    m_n_atoms = 0;
  }

  void beginObject() override {
    MoleculeExporter::beginObject();
    ObjectMoleculeVerifyChemistry(m_iter.obj, m_iter.state);
  }

  void writeAtom() override {
    const auto ai = m_iter.getAtomInfo();
    const float * rgb = ColorGet(G, ai->color);

    char inscode[3] = { ai->inscode, 0 };
    if (!inscode[0]) {
      strcpy(inscode, "<>");
    }

    ResName resn = "";
    AtomName name = "X";
    if (ai->resn)
      AtomInfoGetAlignedPDBResidueName(G, ai, resn);
    if (ai->name)
      AtomInfoGetAlignedPDBAtomName(G, ai, resn, name);

    m_offset += VLAprintf(m_buffer, m_offset,
        "%d %d %.3f %.3f %.3f %d %s %s \"%-4s\" \"%-4s\" %d %d %02X%02X%02X %d %.2f %d\n",
        getTmpID(),
        getMacroModelAtomType(ai),
        m_coord[0], m_coord[1], m_coord[2],
        ai->resv,
        inscode,
        ai->chain ? LexStr(G, ai->chain) : "\" \"",
        resn,
        name,
        ai->protons,
        ai->formalCharge,
        int(rgb[0] * 255),
        int(rgb[1] * 255),
        int(rgb[2] * 255),
        ai->ssType[0] == 'H' ? 1 : ai->ssType[0] == 'S' ? 2 : 0,
        ai->q,
        ai->id);

    char ribbon_color_rgb[7] = "<>";
    MaeExportGetRibbonColor(G, m_iter, ribbon_color_rgb);
    std::string label_user_text = MaeExportGetLabelUserText(G, ai);

    m_offset += VLAprintf(m_buffer, m_offset,
        "%d %d %d %d %s \"%s\" 2 \"%s\"\n",
        (ai->visRep & ~(cRepCartoonBit | cRepRibbonBit)) ? 1 : 0,
        MaeExportGetAtomStyle(G, m_iter),
        MaeExportGetRibbonStyle(ai),
        ribbon_color_rgb[0] == '<' ? 3 /* calphaatom */ : 0 /* constant */,
        ribbon_color_rgb,
        label_user_text.empty() ? "" : "%UT",
        label_user_text.c_str());

    m_atoms[getTmpID()] = ai;

    ++m_n_atoms;
  }

  void writeBonds() override {
    // atom count
    m_n_atoms_offset += sprintf(m_buffer + m_n_atoms_offset, "m_atom[%d]", m_n_atoms);
    m_buffer[m_n_atoms_offset] = ' '; // overwrite terminator

    if (!m_bonds.empty()) {
      // table with zero rows not allowed

      m_offset += VLAprintf(m_buffer, m_offset,
          ":::\n"
          "}\n" // end atoms table
          "m_bond[%d] {\n"
          "# First column is bond index #\n"
          "i_m_from\n"
          "i_m_to\n"
          "i_m_order\n"
          "i_m_from_rep\n"
          "i_m_to_rep\n"
          ":::\n", (int) m_bonds.size());

      int b = 0;
      for (auto& bond : m_bonds) {
        int order = bond.ref->order;
        if (order > 3) {
          ++m_n_arom_bonds;
          order = 1; // aromatic bonds not supported
        }

        m_offset += VLAprintf(m_buffer, m_offset, "%d %d %d %d\n", ++b,
            bond.id1,
            bond.id2,
            order);

        int style = MaeExportGetBondStyle(m_atoms[bond.id1], m_atoms[bond.id2]);
        m_offset += VLAprintf(m_buffer, m_offset, "%d %d\n", style, style);
      }

      m_bonds.clear();
    }

    // end molecule
    m_offset += VLAprintf(m_buffer, m_offset,
        ":::\n"
        "}\n" // end bonds (or atoms) table
        "}\n");

    if (m_n_arom_bonds > 0) {
      PRINTFB(G, FB_ObjectMolecule, FB_Warnings)
        " Warning: aromatic bonds not supported by MAE format, "
        "exporting as single bonds\n" ENDFB(G);
      m_n_arom_bonds = 0;
    }
  }
};

// ---------------------------------------------------------------------------------- //

struct MoleculeExporterXYZ : public MoleculeExporter {
  int m_n_atoms;
  int m_n_atoms_offset;

  int getMultiDefault() const override {
    // multi-entry format
    return cMolExportByCoordSet;
  }

  void beginMolecule() override {
    MoleculeExporter::beginMolecule();

    // defer until number of atoms known
    m_n_atoms = 0;
    m_n_atoms_offset = m_offset;

    m_offset += VLAprintf(m_buffer, m_offset,
        "X         \n" // natoms (deferred)
        "%s\n",        // comment
        getTitleOrName());
  }

  void writeAtom() override {
    const auto ai = m_iter.getAtomInfo();

    m_offset += VLAprintf(m_buffer, m_offset,
        "%s %f %f %f\n", ai->elem, m_coord[0], m_coord[1], m_coord[2]);

    ++m_n_atoms;
  }

  void writeBonds() override {
    // atom count
    m_n_atoms_offset += sprintf(m_buffer + m_n_atoms_offset, "%d", m_n_atoms);
    m_buffer[m_n_atoms_offset] = ' '; // overwrite terminator
  }

  bool isExcludedBond(int atm1, int atm2) override {
    // no bonds for XYZ format
    return true;
  }
};

// ---------------------------------------------------------------------------------- //

#ifndef _PYMOL_NO_MSGPACKC
class MoleculeExporterMMTF : public MoleculeExporter {
  mmtf::StructureData m_raw;
  mmtf::GroupType * m_residue = nullptr;
  const AtomInfoType * m_last_ai = nullptr;
  ElemCanonicalizer m_elemGetter;

  std::vector<int32_t> colorList;
  std::vector<int32_t> repsList;

public:
  int getMultiDefault() const override {
    return cMolExportGlobal;
  }

  void beginCoordSet() override {
    m_raw.chainsPerModel.emplace_back(0);
    m_last_ai = nullptr;
  }

  void writeAtom() override {
    const AtomInfoType * ai = m_iter.getAtomInfo();

    m_raw.xCoordList.emplace_back(m_coord[0]);
    m_raw.yCoordList.emplace_back(m_coord[1]);
    m_raw.zCoordList.emplace_back(m_coord[2]);

    // PyMOL customs
    colorList.emplace_back(ai->color);
    repsList.emplace_back(ai->visRep);

    bool is_same_residue = false;
    bool is_same_chain = AtomInfoSameChainP(G, ai, m_last_ai);

    if (!is_same_chain) {
      m_raw.chainsPerModel.back() += 1;
      m_raw.groupsPerChain.emplace_back(0); // increment with every group
      m_raw.chainIdList.emplace_back(LexStr(G, ai->segi));
      m_raw.chainNameList.emplace_back(LexStr(G, ai->chain));
    } else {
      is_same_residue = AtomInfoSameResidueP(G, ai, m_last_ai);
    }

    if (!is_same_residue) {
      m_raw.groupsPerChain.back() += 1;
      m_raw.groupTypeList.emplace_back(m_raw.groupList.size());
      m_raw.groupIdList.emplace_back(ai->resv);
      m_raw.insCodeList.emplace_back(ai->inscode);
      m_raw.secStructList.emplace_back(
          ai->ssType[0] == 'H' ? 2 :
          ai->ssType[0] == 'S' ? 3 : -1);

      m_raw.groupList.emplace_back();
      m_residue = &m_raw.groupList.back();
      m_residue->groupName = LexStr(G, ai->resn);
    }

    m_residue->formalChargeList.emplace_back(ai->formalCharge);
    m_residue->atomNameList.emplace_back(LexStr(G, ai->name));
    m_residue->elementList.emplace_back(m_elemGetter(ai));

    m_raw.bFactorList.emplace_back(ai->b);
    m_raw.occupancyList.emplace_back(ai->q);
    m_raw.altLocList.emplace_back(ai->alt[0]);

    // don't assign to 'm_retain_ids', it messes up bonds
    if (SettingGetGlobal_b(G, cSetting_pdb_retain_ids)) {
      m_raw.atomIdList.emplace_back(ai->id);
    }

    m_last_ai = ai;
  }

  void writeBonds() override {
    m_raw.numAtoms = m_raw.xCoordList.size();
    m_raw.numGroups = m_raw.groupIdList.size();
    m_raw.numChains = m_raw.chainIdList.size();
    m_raw.numModels = m_raw.chainsPerModel.size();

    mmtf::BondAdder bondadder(m_raw);

    for (const auto &bond : m_bonds) {
      bondadder(bond.id1 - 1, bond.id2 - 1, bond.ref->order);
    }

    mmtf::compressGroupList(m_raw);

    packMsgpack();
  }

  void packMsgpack() {
    msgpack::zone _zone;

    auto data_map = mmtf::encodeToMap(m_raw, _zone);
    data_map["pymolColorList"] = msgpack::object(colorList, _zone);
    data_map["pymolRepsList"] = msgpack::object(repsList, _zone);

    std::stringstream stream;
    msgpack::pack(stream, data_map);
    auto buffer = stream.str();
    auto bufferSize = buffer.size();
    m_buffer.resize(bufferSize);
    std::memcpy(m_buffer, buffer.data(), bufferSize);
    m_offset = bufferSize;
  }
};
#endif

/*========================================================================*/

/*
 * Export the given selection to a molecular file format.
 *
 * Return the file contents as a VLA (must be VLAFree'd by caller) or
 * NULL if the format is not known.
 *
 * format:      pdb, sdf, ...
 * selection:   atom selection expression
 * state:       object state (-1 for all, -2/-3 for current)
 * ref_object:  name of a reference object which defines the frame of
 *              reference for exported coordinates
 * ref_state:   reference object state
 * multi:       defines how to handle selections which span multiple objects
 *              -1: use format-specific default
 *               0: one global "molecule" (default for PDB)
 *               1: molecules per objects
 *               2: molecules per states (default for sdf, mol2)
 */
pymol::vla<char> MoleculeExporterGetStr(PyMOLGlobals * G,
    const char *format,
    const char *selection,
    int state,
    const char *ref_object,
    int ref_state,
    int multi,
    bool quiet)
{
  SelectorTmp tmpsele1(G, selection);
  int sele = tmpsele1.getIndex();

  if (sele < 0)
    return nullptr;

  std::unique_ptr<MoleculeExporter> exporter;

  if (ref_state < -1)
   ref_state = state;

  // do "effective" current states
  if (state == -2)
    state = -3;

  if (strcmp(format, "pdb") == 0) {
    exporter.reset(new MoleculeExporterPDB);
  } else if (strcmp(format, "pmcif") == 0) {
    exporter.reset(new MoleculeExporterPMCIF);
  } else if (strcmp(format, "cif") == 0) {
    exporter.reset(new MoleculeExporterCIF);
  } else if (strcmp(format, "sdf") == 0) {
    exporter.reset(new MoleculeExporterSDF);
  } else if (strcmp(format, "pqr") == 0) {
    exporter.reset(new MoleculeExporterPQR);
  } else if (strcmp(format, "mol2") == 0) {
    exporter.reset(new MoleculeExporterMOL2);
  } else if (strcmp(format, "mol") == 0) {
    exporter.reset(new MoleculeExporterMOL);
  } else if (strcmp(format, "xyz") == 0) {
    exporter.reset(new MoleculeExporterXYZ);
  } else if (strcmp(format, "mae") == 0) {
    exporter.reset(new MoleculeExporterMAE);
  } else if (strcmp(format, "mmtf") == 0) {
#ifndef _PYMOL_NO_MSGPACKC
    exporter.reset(new MoleculeExporterMMTF);
#else
    PRINTFB(G, FB_ObjectMolecule, FB_Errors)
      " Error: This build has no fast MMTF support.\n" ENDFB(G);
    return nullptr;
#endif
  } else {
    PRINTFB(G, FB_ObjectMolecule, FB_Errors)
      " Error: unknown format: '%s'\n", format ENDFB(G);
    return nullptr;
  }

  // Ensure "." decimal point in printf. It's possible to change this from
  // Python, so don't rely on a persistent global value.
  std::setlocale(LC_NUMERIC, "C");

  exporter->init(G);
  exporter->setMulti(multi);
  exporter->setRefObject(ref_object, ref_state);
  exporter->execute(sele, state);

  return std::move(exporter->m_buffer);
}

/*========================================================================*/

#ifndef _PYMOL_NOPY
/*
 * See MoleculeExporterGetPyBonds
 */
struct MoleculeExporterPyBonds : public MoleculeExporter {
  PyObject *m_bond_list; // out

protected:
  int getMultiDefault() const override {
    // single-entry format
    return cMolExportGlobal;
  }

  void writeAtom() override {}

  void writeBonds() override {
    size_t nBond = m_bonds.size();
    m_bond_list = PyList_New(nBond);

    for (size_t b = 0; b < nBond; ++b) {
      const auto& bond = m_bonds[b];
      PyList_SetItem(m_bond_list, b, Py_BuildValue("iii",
            bond.id1 - 1, bond.id2 - 1, bond.ref->order));
    }

    m_bonds.clear();
  }
};

/*========================================================================*/

/*
 * Get all bonds within the given selection.
 *
 * Returns a list of (atm1, atm2, order) tuples, where atm1 and atm2 are
 * the 0-based atom indices within the selection.
 *
 * Return value: New reference
 */
PyObject *MoleculeExporterGetPyBonds(PyMOLGlobals * G,
    const char *selection, int state)
{
  SelectorTmp tmpsele1(G, selection);
  int sele = tmpsele1.getIndex();
  if (sele < 0)
    return nullptr;

  int unblock = PAutoBlock(G);

  MoleculeExporterPyBonds exporter;
  exporter.init(G);
  exporter.execute(sele, state);

  if (PyErr_Occurred())
    PyErr_Print();

  PAutoUnblock(G, unblock);

  return exporter.m_bond_list;
}

/*========================================================================*/

/*
 * Creates a `chempy.models.Indexed` instance from a selection.
 *
 * Changes from the old implementation (SelectorGetChemPyModel):
 * - drops support for `cs->Spheroid`
 *
 */
struct MoleculeExporterChemPy : public MoleculeExporter {
  PyObject *m_model; // out

protected:
  int m_n_cs; // number of coordinate sets
  float m_ref_tmp[3];
  PyObject *m_atom_list;

public:
  // quasi constructor
  void init(PyMOLGlobals * G_) override {
    MoleculeExporter::init(G_);

    m_model = nullptr;
    m_n_cs = 0;
    m_atom_list = nullptr;
  }

protected:
  int getMultiDefault() const override {
    // single-entry format
    return cMolExportGlobal;
  }

  void beginMolecule() override {
    MoleculeExporter::beginMolecule();

    m_model = PYOBJECT_CALLMETHOD(P_models, "Indexed", "");
    if (m_model) {
      m_atom_list = PyList_New(0);
      PyObject_SetAttrString(m_model, "atom", m_atom_list);
      Py_DECREF(m_atom_list);
    }
  }

  void beginCoordSet() override {
    MoleculeExporter::beginCoordSet();
    ++m_n_cs;
  }

  const float * getRefPtr() {
    RefPosType  * ref_pos = m_iter.cs->RefPos;
    const float * ref_ptr = nullptr;

    if (ref_pos) {
      ref_pos += m_iter.getIdx();
      if (ref_pos->specified) {
        ref_ptr = ref_pos->coord;
        if (m_mat_move.ptr) {
          transform44d3f(m_mat_move.ptr, ref_ptr, m_ref_tmp);
          ref_ptr = m_ref_tmp;
        }
      }
    }

    return ref_ptr;
  }

  void writeAtom() override {
    PyObject *atom = CoordSetAtomToChemPyAtom(G,
        m_iter.getAtomInfo(), m_coord, getRefPtr(),
        m_iter.getAtm(), m_mat_full.ptr);

    if (atom) {
      PyList_Append(m_atom_list, atom);
      Py_DECREF(atom);
    }
  }

  void writeProperties() {
    if (!m_last_cs)
      return;

    // only support properties if export doesn't span multiple coordsets
    if (m_n_cs != 1)
      return;

    // title
    if (m_last_cs->Name[0]) {
      PyObject *molecule = PyObject_GetAttrString(m_model, "molecule");
      if (molecule) {
        PyObject_SetAttrString(molecule, "title", PyString_FromString(m_last_cs->Name));
        Py_DECREF(molecule);
      }
    }

    // properties
#ifdef _PYMOL_IP_PROPERTIES
#endif
  }

  void writeBonds() override {
    if (!m_model)
      return;

    bool error = false;
    size_t nBond = m_bonds.size();
    PyObject *bond_list = PyList_New(nBond);

    for (size_t b = 0; b < nBond; ++b) {
      PyObject *bnd = PYOBJECT_CALLMETHOD(P_chempy, "Bond", "");
      if (!bnd) {
        error = true;
        break;
      }

      const auto& bond = m_bonds[b];
      int index[] = { bond.id1 - 1, bond.id2 - 1 };
      PConvInt2ToPyObjAttr(bnd, "index", index);
      PConvIntToPyObjAttr(bnd, "order",     bond.ref->order);
      PConvIntToPyObjAttr(bnd, "id",        bond.ref->id);
      PConvIntToPyObjAttr(bnd, "stereo",    bond.ref->stereo);
      PyList_SetItem(bond_list, b, bnd);    /* steals bnd reference */
    }

    if (!error) {
      PyObject_SetAttrString(m_model, "bond", bond_list);
    }

    Py_DECREF(bond_list);
    m_bonds.clear();

    writeProperties();
  }
};

/*========================================================================*/

PyObject *ExecutiveSeleToChemPyModel(PyMOLGlobals * G,
    const char *s1, int state,
    const char *ref_object, int ref_state)
{
  if (state == -1)
    state = 0; // no multi-state support

  if (ref_state < -1)
    ref_state = state;

  int sele = SelectorIndexByName(G, s1);
  if (sele < 0)
    return nullptr;

  int unblock = PAutoBlock(G);

  MoleculeExporterChemPy exporter;
  exporter.init(G);
  exporter.setRefObject(ref_object, ref_state);
  exporter.execute(sele, state);

  if (PyErr_Occurred())
    PyErr_Print();

  PAutoUnblock(G, unblock);

  return exporter.m_model;
}
#endif

// vi:sw=2
