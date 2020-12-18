/*
 * MOL/SDF V3000 input support for PyMOL.
 *
 * According to:
 * http://c4.cabrillo.edu/404/ctfile.pdf
 * http://infochim.u-strasbg.fr/recherche/Download/Fragmentor/MDL_SDF.pdf
 *
 * (c) Schrodinger, Inc.
 */

#include <string>

#include "os_std.h"

#include "MolV3000.h"

#include "MemoryDebug.h"
#include "Feedback.h"
#include "Parse.h"
#include "Vector.h"
#include "Rep.h"
#include "Lex.h"
#include "strcasecmp.h"

/**
 * Read a MOL/SDF V3000 line into the `out` buffer, stripping the leading
 * "M  V30 " and taking line continuation (trailing hyphen) into account.
 *
 * Moves the input pointer `p` to the next line.
 *
 * Return false if line doesn't start with "M  V30 ".
 */
static
bool MOLV3000ReadLine(const char * &p, std::string &out) {
  out.clear();

  for (bool continued = true; continued;) {
    if (strncmp(p, "M  V30 ", 7))
      return false;

    // find beginning of next line
    p += 7;
    auto next = ParseNextLine(p);

    // find end of line
    auto last = next;
    if (last > p && last[-1] == '\n') --last;
    if (last > p && last[-1] == '\r') --last;

    // check for line continuation (hyphen)
    if ((continued = last > p && last[-1] == '-')) {
      --last;
    }

    // append to out buffer
    out.append(p, last);

    p = next;
  }

  return true;
}

/**
 * Parse the next keyword=value pair from `p` and advance `p`.
 */
static
bool MOLV3000ReadKeyValue(const char *& p,
    std::string &key,
    std::string &value)
{
  // skip whitespace
  while (*p && (*p == ' ' || *p == '\t')) ++p;

  auto begin = p;
  auto termchars = " \t";

  // parse key
  for (;; ++p) {
    if (!*p)
      return false;

    if (*p == '=') {
      key.assign(begin, p);
      begin = ++p;
      if (*begin == '(') {
        termchars = ")";
      }
      break;
    }
  }

  // parse value
  while (!strchr(termchars, *p)) ++p;
  if (*begin == '(' && *p == ')') ++p;
  value.assign(begin, p);

  return true;
}

/**
 * Parse the "Extended Connection Table" (V3000) from a MOL/SDF file.
 *
 * Returns a pointer to the end of the parsed content on success, or NULL on
 * failure.
 *
 * buffer: input buffer, should begin with "M  V30 "
 * atInfo: atom VLA
 * bond:   bond VLA
 * coord:  coordinates VLA
 * nAtom:  output variable for number of atoms
 * nBond:  output variable for number of bonds
 */
const char * MOLV3000Parse(PyMOLGlobals * G,
    const char * buffer,
    AtomInfoType *& atInfo,
    BondType     *& bond,
    float        *& coord,
    int & nAtom,
    int & nBond)
{
  char cc[16]; // buffer for short words like "BEGIN" or "COUNTS"
  bool inside_atom = false;
  bool inside_bond = false;
  bool inside_any  = false;
  int auto_show = RepGetAutoShowMask(G);
  const char * error = nullptr;
  AtomInfoType * ai = nullptr;
  std::string line;
  std::string key, value;

  while (MOLV3000ReadLine(buffer, line)) {
    auto p = line.c_str();

    // save position after "V30"
    auto p_data = p;

    p = ParseWordCopy(cc, p, sizeof(cc));
    bool is_end = (strcasecmp(cc, "END") == 0);

    if (inside_any) {
      // skip
      if (is_end) {
        inside_any = false;
      }
    } else if (inside_atom) {
      if (is_end) {
        inside_atom = false;
        continue;
      }

      int index, p_offset;
      char type[4];
      float xyz[3];
      auto n = sscanf(p_data, "%d %3s %f %f %f%n %*d%n", &index, type,
            xyz, xyz + 1, xyz + 2, &p_offset, &p_offset);

      if (n != 5) {
        error = "failed to parse atom line";
        break;
      }

      p = p_data + p_offset;

      if (index < 1 || index > nAtom) {
        error = "atom index out of range";
        break;
      }

      if (atInfo) {
        ai = atInfo + index - 1;
        ai->name = LexIdx(G, type);
        ai->visRep = auto_show;
        ai->hetatm = true;
        ai->id = index;
        ai->rank = index - 1;

        copy3(xyz, coord + 3 * (index - 1));

        AtomInfoAssignParameters(G, ai);
        AtomInfoAssignColors(G, ai);

        while (MOLV3000ReadKeyValue(p, key, value)) {
          if (key == "CHG") {
            ai->formalCharge = atoi(value.c_str());
          } else if (key == "CFG") {
            ai->stereo = atoi(value.c_str());
          }
        }
      }
    } else if (inside_bond) {
      if (is_end) {
        inside_bond = false;
        continue;
      }

      int index, type, atom1, atom2, p_offset;
      auto n = sscanf(p_data, "%d %d %d %d%n", &index, &type,
          &atom1, &atom2, &p_offset);

      if (n != 4) {
        error = "failed to parse bond line";
        break;
      }

      if (bond) {
        if (index < 1 || index > nBond) {
          error = "bond index out of range";
          break;
        }

        if (type == 7) {
          type = 2; // double or aromatic
        } else if (type > 4) {
          type = 1; // 5: single or double, 6: single or aromatic, 8: any
        }

        BondTypeInit2(bond + index - 1, atom1 - 1, atom2 - 1, type);
      }

      p = p_data + p_offset;
    } else if (strcasecmp(cc, "BEGIN") == 0) {
      p = ParseWordCopy(cc, p, sizeof(cc));
      if (strcasecmp(cc, "CTAB") == 0) {
        // ignore
      } else if (strcasecmp(cc, "ATOM") == 0) {
        inside_atom = true;
      } else if (strcasecmp(cc, "BOND") == 0) {
        inside_bond = true;
      } else {
        inside_any = true;
      }
    } else if (strcasecmp(cc, "COUNTS") == 0) {
      if (sscanf(p, "%d %d", &nAtom, &nBond) != 2) {
        error = "COUNTS parsing failed";
        break;
      } else {
        if(atInfo)
          VLACheck(atInfo, AtomInfoType, nAtom);
        if(coord)
          VLACheck(coord, float, 3 * nAtom);
        if(bond)
          VLACheck(bond, BondType, nBond);
      }
    } else {
      // ignore
    }
  }

  if (!error && (inside_atom || inside_bond)) {
    error = "expected 'M  V30'";
  }

  if (error) {
    PRINTFB(G, FB_ObjectMolecule, FB_Errors)
      " MOL-V3000-Error: %s.\n", error ENDFB(G);
    return nullptr;
  }

  return buffer;
}
