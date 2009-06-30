/* MACHINE GENERATED FILE, DO NOT EDIT! */

#define VMDPLUGIN molfile_maeffplugin
#define STATIC_PLUGIN 1

//
// Version info for VMD plugin tree:
//   $Id: maeffplugin.cxx,v 1.19 2009/05/18 16:09:21 johns Exp $
//
// Version info for last sync with D. E. Shaw Research:
//  //depot/desrad/main/sw/libs/vmd_plugins,DIST/maeffplugin.cxx#6
//

/*
Copyright 2009, D. E. Shaw Research, LLC
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are
met:

* Redistributions of source code must retain the above copyright
  notice, this list of conditions, and the following disclaimer.

* Redistributions in binary form must reproduce the above copyright
  notice, this list of conditions, and the following disclaimer in the
  documentation and/or other materials provided with the distribution.

* Neither the name of D. E. Shaw Research, LLC nor the names of its
  contributors may be used to endorse or promote products derived from
  this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
"AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

*/

#if defined(_MSC_VER)
#ifndef DESRES_WIN32
#define DESRES_WIN32
#endif
#endif

#include <molfile_plugin.h>

#include <stdlib.h>
#include <string.h>

#include <fstream>
#include <vector>
#include <map>
#include <cmath>
#include <set>

#include <iostream>
#include <sstream>
#include <string>
#include <stdexcept>
#include <vector>

#ifdef DESRES_WIN32
#define ssize_t int
#define M_PI (3.1415926535897932385)
#define M_PI_2 (1.5707963267948966192)

#if defined(_MSC_VER)
#ifndef snprintf
#define snprintf _snprintf
#endif
#endif

#endif

namespace {

    /*!
     * \brief Takes a stream and returns maestro tokens.
     * This tokenizer is built on streams and uses a small, tight
     * finite state automata to construct a token
     */
    class Tokenizer {

      /*! \brief Actions for the DFA token builder
       *
       */
      typedef enum {
        DONE = 0,
        SKIPWHITE,  // 1 
        INCOMMENT,  //  2
        CHOOSEKIND,  // 3
        SINGLECHAR,  // 4
        STARTSTRING,  // 5 
        INSTRING,  // 6
        ESCAPE,  // 7 
        STARTOTHER,  // 8 
        CONTINUEOTHER  // 9
      } action_t;

      //! \brief The current character
      char m_c;

      //! \brief the stream for the file we're parsing
      std::ifstream &m_input;

      //! \brief The current token
      char * m_token;

      //! \brief number of malloc'ed bytes in m_token
      ssize_t max_token_size;

      //! \brief True iff the token is already read
      bool m_isfresh;

      //! \brief Current line in file
      unsigned m_line;

      //! \brief Line where token starts
      unsigned m_tokenline;

      //! \brief Get current character
      /*!
      * Returns the current character in the file
      * @return current char
      */
      inline char peek() { return m_c; }

      //! \brief Read a new character
      /*!
      * Read into the current character and update line
      * and point information.
      * @return (new) current char
      */
      inline char read() {
        m_c = m_input.get();
        if (m_c == '\n') m_line++;
        return m_c;
      }

      //! \brief True if character is a token
      /*!
      * A few Maestro tokens are just 1 char long.  We use
      * this test to short-circuit the tokenizer when at
      * [, ], {, or }.
      */
      static bool issingle(char c) {
        return (c == '[' || c == ']' || c == '{' || c == '}');
      }

      Tokenizer(const Tokenizer&); // No copy c'tor

    public:

      //! \brief Tokenizer based on a file stream
      Tokenizer(std::ifstream &in);

      //! \brief Clean up and release buffers.
      ~Tokenizer();

      //! \brief Current token under cursor
      const char * token();

      //! \brief Advance to next token
      void next();

      //! \brief Line number associated with current token
      unsigned line() const;

      //! \brief File seek point for current token
      size_t point() const;

      //! \brief Predict a particular token
      const char * predict(const char * match="");

      //! \brief For while(not_a(match)) loops
      bool not_a(const char * match=END_OF_FILE);

      static const char * END_OF_FILE;
    };
}

/*!
 * Build from an istream
 * @param input The stream to parse
 */
Tokenizer::Tokenizer(std::ifstream &in)
  : m_c(0),
    m_input(in),
    m_token(NULL),
    max_token_size(0),
    m_isfresh(false),
    m_line(1),
    m_tokenline(1)
{
  max_token_size = 16;
  m_token = (char *)malloc(max_token_size);

  // grab 1st token.
  read();
}

/*!
 * The destructor cleans up any heap allocated temporaries created
 * during construction.
 */
Tokenizer::~Tokenizer() {
  if (m_token) free(m_token);
}

/*!
 * This routine assembles a token character-by-character.
 * At its heart is a little DFA that looks at a character
 * to determine the next state.  For instance, the DFA
 * starts in a SKIPWHITE state and stays there unless
 * it finds a comment (#) or other character.  In the
 * INCOMMENT state, it looks for end-of-line before
 * returning to SKIPWHITE.  Similarly, all tokens are
 * defined using this one-character lookahead.
 * @return The current (possibly new) token 
 */
const char * Tokenizer::token() {
  // -----------------------------------------------
  // Keep returning the same token until next()
  // is called.
  // -----------------------------------------------
  if (m_isfresh) return m_token;
  
  // -----------------------------------------------
  // End of file simply returns an empty string
  // -----------------------------------------------

  // begin at start of token space
  char * ptr = m_token;
  m_isfresh = true;
  
  unsigned state = SKIPWHITE;
  char c = peek();
  bool good = false;
  ssize_t diff;
  while(state != DONE && c >= 0) {
    // make sure we have space in m_token for 2 more characters
    if ((diff = ptr-m_token) >= max_token_size-1) {
      m_token = (char *)realloc( m_token, 2*max_token_size );
      ptr = m_token + diff;
      max_token_size *= 2;
    }
    //std::cout << int(c) << " " << state << std::endl;
    switch(state) {
    case SKIPWHITE:
      // -----------------------------------------------
      // Skip whitespace and see if its a token or a
      // comment
      // -----------------------------------------------
      if (::isspace(c)) {
        c = read();
      } else if (c == '#') {
        state = INCOMMENT;
        c = read();
      } else {
        state = CHOOSEKIND;
      }
      break;
    case INCOMMENT:
      // -----------------------------------------------
      // On a comment, read until end of line
      // -----------------------------------------------
      if (c == '\n' || c == '#') state = SKIPWHITE;
      c = read();
      break;
    case CHOOSEKIND:
      // -----------------------------------------------
      // We []{} are single character tokens,
      // Strings start with "
      // Everything else starts with some other character
      // -----------------------------------------------
      if (issingle(c)) {
        state = SINGLECHAR;
      } else if (c == '"') {
        state = STARTSTRING;
      } else {
        state = STARTOTHER;
      }
      break;
    case SINGLECHAR:
      good = true;
      m_tokenline = m_line;
      *ptr++ = c;
      *ptr++ = '\0';
      read();
      state = DONE;
      break;
    case STARTSTRING:
      good = true;
      m_tokenline = m_line;
      *ptr++ = c;
      read(); // Skip opening quote
      c = peek();
      state = INSTRING;
      break;
    case INSTRING:
      if ( c == '"' ) {
        *ptr++ = c;
        *ptr++ = '\0';
        state = DONE;
      } else if ( c == '\\' ) {
        state = ESCAPE;
      } else {
        *ptr++ = c;
      }
      c = read();
      break;
    case ESCAPE:
      *ptr++ = c;
      state = INSTRING;
      c = read();
      break;
    case STARTOTHER:
      good = true;
      m_tokenline = m_line;
      state = CONTINUEOTHER;
      break;
    case CONTINUEOTHER:
      if ( issingle(c) || ::isspace(c) || c == '#' || c == '"' ) {
        *ptr++ = '\0';
        state = DONE;
      } else {
        *ptr++ = c;
        c = read();
      }
      break;
    }
  }
  
  // -----------------------------------------------
  // Maybe we just read trailing whitespace...
  // -----------------------------------------------
  if (!good) *m_token = '\0';
  
  return m_token;
}

/*!
 * Set state to read a new token on the next request.
 */
void Tokenizer::next() {
  m_isfresh = false;
}

/*!
 * Line associated with current token
 * @return The line number
 */
unsigned Tokenizer::line() const {
  return m_tokenline;
}

/*!
 * The predictive, recursive-descent parsers I use do a lot
 * of "I expect the next token to look like this" calls.
 * (e.g. I expect a "{" here).  This simplifies the logic
 * for that.
 * @param match
 * @return The matching token body
 */
const char *
Tokenizer::predict(const char * match) {
  const char * tok = token();
  if (strcmp(match, "") && strcmp(tok, match)) {
    std::stringstream str;
    str << "Line " <<line()<< " predicted '" <<std::string(match)<< "' have '"
        << (isprint(tok[0])?tok:"<unprintable>")
        << "'" << std::endl;
    throw std::runtime_error(str.str());
  }
  next();
  return tok;
}

/*!
 * Another common pattern is replication.  So, while (not_a("}")) { ... }
 * This function makes that easy.
 * @param match Token body to try to match
 * @return True on a match, False on EOF or a non-match
 */
bool Tokenizer::not_a(const char * match) {
  const char * tok = token();
  if (!strcmp(tok, END_OF_FILE)) return false; // EOF always quits
  return strcmp(tok, match);
}

/*!
 * The special end of file token (string of length 1 containing a null byte)
 * Normal parsing will not create this character sequence, so it makes a
 * good special token.
 */
const char * Tokenizer::END_OF_FILE = "";

namespace {
  struct element {
    double daltons;
    const char* abbreviation;
    const char* name;
  };
}

// url = "http://physics.nist.gov/cgi-bin/Elements/elInfo.pl?element=%d&context=noframes"%element
static struct element amu[] = {
  {1.00794,"H","Hydrogen"},
  {4.002602,"He","Helium"},
  {6.941,"Li","Lithium"},
  {9.012182,"Be","Beryllium"},
  {10.811,"B","Boron"},
  {12.0107,"C","Carbon"},
  {14.0067,"N","Nitrogen"},
  {15.9994,"O","Oxygen"},
  {18.9984032,"F","Fluorine"},
  {20.1797,"Ne","Neon"},
  {22.989770,"Na","Sodium"},
  {24.3050,"Mg","Magnesium"},
  {26.981538,"Al","Aluminum"},
  {28.0855,"Si","Silicon"},
  {30.973761,"P","Phosphorus"},
  {32.065,"S","Sulfur"},
  {35.453,"Cl","Chlorine"},
  {39.0983,"K","Potassium"},
  {39.948,"Ar","Argon"},
  {40.078,"Ca","Calcium"},
  {44.955910,"Sc","Scandium"},
  {47.867,"Ti","Titanium"},
  {50.9415,"V","Vanadium"},
  {51.9961,"Cr","Chromium"},
  {54.938049,"Mn","Manganese"},
  {55.845,"Fe","Iron"},
  {58.6934,"Ni","Nickel"},
  {58.933200,"Co","Cobalt"},
  {63.546,"Cu","Copper"},
  {65.409,"Zn","Zinc"},
  {69.723,"Ga","Gallium"},
  {72.64,"Ge","Germanium"},
  {74.92160,"As","Arsenic"},
  {78.96,"Se","Selenium"},
  {79.904,"Br","Bromine"},
  {83.798,"Kr","Krypton"},
  {85.4678,"Rb","Rubidium"},
  {87.62,"Sr","Strontium"},
  {88.90585,"Y","Yttrium"},
  {91.224,"Zr","Zirconium"},
  {92.90638,"Nb","Niobium"},
  {95.94,"Mo","Molybdenum"},
  {101.07,"Ru","Ruthenium"},
  {102.90550,"Rh","Rhodium"},
  {106.42,"Pd","Palladium"},
  {107.8682,"Ag","Silver"},
  {112.411,"Cd","Cadmium"},
  {114.818,"In","Indium"},
  {118.710,"Sn","Tin"},
  {121.760,"Sb","Antimony"},
  {126.90447,"I","Iodine"},
  {127.60,"Te","Tellurium"},
  {131.293,"Xe","Xenon"},
  {132.90545,"Cs","Cesium"},
  {137.327,"Ba","Barium"},
  {138.9055,"La","Lanthanum"},
  {140.116,"Ce","Cerium"},
  {140.90765,"Pr","Praseodymium"},
  {144.24,"Nd","Neodymium"},
  {150.36,"Sm","Samarium"},
  {151.964,"Eu","Europium"},
  {157.25,"Gd","Gadolinium"},
  {158.92534,"Tb","Terbium"},
  {162.500,"Dy","Dysprosium"},
  {164.93032,"Ho","Holmium"},
  {167.259,"Er","Erbium"},
  {168.93421,"Tm","Thulium"},
  {173.04,"Yb","Ytterbium"},
  {174.967,"Lu","Lutetium"},
  {178.49,"Hf","Hafnium"},
  {180.9479,"Ta","Tantalum"},
  {183.84,"W","Tungsten"},
  {186.207,"Re","Rhenium"},
  {190.23,"Os","Osmium"},
  {192.217,"Ir","Iridium"},
  {195.078,"Pt","Platinum"},
  {196.96655,"Au","Gold"},
  {200.59,"Hg","Mercury"},
  {204.3833,"Tl","Thallium"},
  {207.2,"Pb","Lead"},
  {208.98038,"Bi","Bismuth"},
  {231.03588,"Pa","Protactinium"},
  {232.0381,"Th","Thorium"},
  {238.02891,"U","Uranium"}
};

static const int nelements = sizeof(amu)/sizeof(amu[0]);

static std::pair<int, const char *> 
find_element_by_amu(double target) {
  int left = 0;
  int right = nelements-1;

  // -----------------------------------------------
  // Knuth's binary search
  // -----------------------------------------------
  while(left <= right) {
    int mid = (left+right)/2;
    if (target> amu[mid].daltons) {
      left = mid + 1;
    } else if (target< amu[mid].daltons) {
      right = mid - 1;
    } else {
      /* Exact match (unlikely) */
      left = right = mid;
      return std::pair<int, const char *>( left+1, amu[left].abbreviation );
    }
  }

  // -----------------------------------------------
  // CAUTION: at this point, the meanings of 
  // left and right are switched (i.e. left >= right,
  // see the while() loop above if you don't believe me!
  // -----------------------------------------------
  int swap = left;
  left = right;
  right = swap;

  if (left < 0) left = right;
  if (right > nelements-1) right = left;

  if (target - amu[left].daltons < amu[right].daltons - target) {
    return std::pair<int, const char *>( left+1, amu[left].abbreviation );
  }

  return std::pair<int, const char *>( right+1, amu[left].abbreviation );
}

static std::pair<double, const char *> 
find_element_by_atomic_number(int target) {
  if (target < 1) target=1;
  if (target >= nelements) target = nelements-1;
  return std::pair<double,const char *> ( amu[target-1].daltons,
                                          amu[target-1].abbreviation );
}

namespace {

  struct schema_t { /*GCOV-IGNORE*/
    char type;
    std::string attr;
  };

  struct site {
    float mass;
    float charge;
    bool  pseudo;
    site() : mass(0), charge(0), pseudo(0) {}
  };

  struct vsite {
    int ai;             // parent atom of virtual
    std::string funct;  // ffio_funct
  };
  // mapping from index in sites to entry in ffio_virtuals
  typedef std::map<int,vsite> VirtualsMap;

  struct bond_t {
    int from, to;
    float order;
    bond_t() {}
    bond_t(int f, int t, float o) : from(f), to(t), order(o) {}
  };
  
  struct pos_t {
    float x, y, z;
    pos_t() {}
    pos_t(float x_, float y_, float z_) : x(x_), y(y_), z(z_) {}
  };
  struct vel_t {
    float x, y, z;
    vel_t() {}
    vel_t(float x_, float y_, float z_) : x(x_), y(y_), z(z_) {}
  };

  struct ct_data {
    int natoms;    // number of physical atoms from m_atom block
    int npseudos;  // number of pseudoatoms from ffio_pseudo

    std::vector<molfile_atom_t> particles;
    std::vector<pos_t>          position;
    std::vector<vel_t>          velocity;
    std::vector<site>           sites;
    std::vector<bond_t>         bonds;

    std::map<size_t,int> atommap, pseudomap;

    VirtualsMap virtuals;
    ct_data() : natoms(0), npseudos(0) {}
  };

  typedef std::vector<schema_t> Schema;
  typedef std::vector<std::string> Row;
  typedef std::map<std::string, std::string> AttrMap;
  typedef std::map<int, ct_data> CtMap;

  /*! entries in fepio_fep table */
  struct fep_elem {
    //! variable names correspond to fepio_fep field names
    //!@{
    int ti, tj;
    int ai, aj, ak, al;
    int am, an, ao, ap; // Allow treatment of torsion-torsion terms
    int moiety; //!< The moiety that each mapped term belongs to.
    fep_elem() : ti(-1), tj(-1), ai(-1), aj(-1), ak(-1), al(-1),
                 am(-1), an(-1), ao(-1), ap(-1), moiety(-1) {}
    //!@}
  };
  /*! Represnts a single table in fepio_fep */
  typedef std::vector<fep_elem> FepList;
  /*! Represnts the entire fepio_fep section */
  typedef std::map<std::string, FepList > FepioMapping;

  double dotprod(const double *x, const double *y) {
    return x[0]*y[0] + x[1]*y[1] + x[2]*y[2];
  }

  struct Handle {
    std::ofstream output;
    bool eof;
    double A[3], B[3], C[3];        // global cell vectors
    int optflags;

    int stage1, stage2;  // alchemical stages
    FepioMapping fepmap;

    int nparticles;
    std::vector<int> bond_from, bond_to;
    std::vector<float> bond_order;
    std::vector<molfile_atom_t> particles; // for writing
    CtMap ctmap;

    Handle() 
    : eof(false), stage1(0), stage2(0), nparticles(0) {
      for (int i=0; i<3; i++) A[i] = B[i] = C[i] = 0;
      A[0] = 1;
      B[1] = 1;
      C[2] = 1;
    }

    void set_box( AttrMap &attrs ) {
      std::string abox("chorus_box_a_");
      std::string bbox("chorus_box_b_");
      std::string cbox("chorus_box_c_");
      for (int i=0; i<3; i++) {
        char x='x'+i;
        abox[abox.size()-1] = x;
        bbox[bbox.size()-1] = x;
        cbox[cbox.size()-1] = x;
        A[i] = atof(attrs[abox].c_str());
        B[i] = atof(attrs[bbox].c_str());
        C[i] = atof(attrs[cbox].c_str());
      }
    }
    void set_box( const molfile_timestep_t *ts) {
      // Convert VMD's unit cell information
      double cosBC = sin( ((90 - ts->alpha ) / 180) * M_PI );
      double cosAC = sin( ((90 - ts->beta  ) / 180) * M_PI );
      double cosAB = sin( ((90 - ts->gamma ) / 180) * M_PI );
      double sinAB = cos( ((90 - ts->gamma ) / 180) * M_PI );

      double Ax = ts->A;
      double Ay = 0;
      double Az = 0;
      double Bx = ts->B * cosAB;
      double By = ts->B * sinAB;
      double Bz = 0;
      double Cx,Cy,Cz;
      if (sinAB != 0) {
        Cx = cosAC;
        Cy = (cosBC - cosAC*cosAB) / sinAB;
        Cz = sqrt(1-Cx*Cx-Cy*Cy);
        Cx *= ts->C;
        Cy *= ts->C;
        Cz *= ts->C;
      } else {
        Cx=Cy=Cz=0;
      }
      A[0] = Ax; A[1] = Ay; A[2] = Az;
      B[0] = Bx; B[1] = By; B[2] = Bz;
      C[0] = Cx; C[1] = Cy; C[2] = Cz;
    }
    void get_box( molfile_timestep_t *ts) const {
      ts->A = sqrt(dotprod(A,A));
      ts->B = sqrt(dotprod(B,B));
      ts->C = sqrt(dotprod(C,C));

      if (ts->A == 0 || ts->B == 0 || ts->C == 0) {
        // bogosity in the unit cell.  Warn the user and
        // provide 90 degree angles for all.
        fprintf(stderr, "WARNING: Some unit cell dimensions were zero; all unit cell angles set to 90.\n");
        ts->alpha = ts->beta = ts->gamma = 90;
      } else {

        // compute angles
        double cosAB = dotprod(A,B)/(ts->A * ts->B);
        double cosAC = dotprod(A,C)/(ts->A * ts->C);
        double cosBC = dotprod(B,C)/(ts->B * ts->C);

        // clamp
        if (cosAB > 1.0) cosAB = 1.0; else if (cosAB < -1.0) cosAB = -1.0;
        if (cosAC > 1.0) cosAC = 1.0; else if (cosAC < -1.0) cosAC = -1.0;
        if (cosBC > 1.0) cosBC = 1.0; else if (cosBC < -1.0) cosBC = -1.0;

        // convert to angles using asin to avoid nasty rounding when we are
        // close to 90 degree angles.
        ts->alpha = 90.0 - asin(cosBC) * 90.0 / M_PI_2; /* cosBC */
        ts->beta  = 90.0 - asin(cosAC) * 90.0 / M_PI_2; /* cosAC */
        ts->gamma = 90.0 - asin(cosAB) * 90.0 / M_PI_2; /* cosAB */
      }
    }
  };

  class Array {
  protected:
    Handle *h;
    const int m_ct;

  public:
    Array(Handle *h_, int ct) : h(h_), m_ct(ct) {}
    virtual ~Array() {}

    virtual void set_schema( const Schema &schema ) {}
    virtual void insert_row(const Row &row) {}

    static void get_str(const std::string &value, char *arr, int N) {
      if (value == "<>") return;
      if (value.size() && value[0] == '"' && value[value.size()-1]) {
        strncpy(arr, value.substr(1,value.size()-2).c_str(), N);
      } else {
        strncpy(arr, value.c_str(), N);
      }
    }

    static void get_int(const std::string &value, int &ival) {
      ival=atoi(value.c_str());
    }
    static void get_float(const std::string &value, float &fval) {
      fval=atof(value.c_str());
    }
  };
#define GET_STR( val_, arr_) do { get_str(val_, arr_, sizeof(arr_)); } while(0)

  class AtomArray : public Array {
    int i_name, i_resname, i_resid, i_x, i_y, i_z, i_vx, i_vy, i_vz, 
        i_anum, i_chain, i_segid;
    std::vector<molfile_atom_t> &atoms;
    std::vector<pos_t> &pos;
    std::vector<vel_t> &vel;
    int &natoms;

  public:
    AtomArray(Handle *h_, int ct) 
    : Array(h_, ct),
      i_name(-1), i_resname(-1), i_resid(-1), 
      i_x(-1), i_y(-1), i_z(-1), 
      i_vx(-1), i_vy(-1), i_vz(-1),
      i_anum(-1), i_chain(-1), i_segid(-1),
      atoms( h->ctmap[m_ct].particles ),
      pos( h->ctmap[m_ct].position),
      vel( h->ctmap[m_ct].velocity),
      natoms( h->ctmap[m_ct].natoms )
    {
      h->optflags = MOLFILE_INSERTION;
    }
    virtual void set_schema( const Schema &schema ) {
      for (unsigned i=0; i<schema.size(); i++) {
        const std::string &attr=schema[i].attr;
        if      (attr=="m_pdb_atom_name")    i_name=i;
        else if (attr=="m_pdb_residue_name") i_resname=i;
        else if (attr=="m_residue_number")   i_resid=i;
        else if (attr=="m_x_coord")          i_x=i;
        else if (attr=="m_y_coord")          i_y=i;
        else if (attr=="m_z_coord")          i_z=i;
        else if (attr=="ffio_x_vel")         i_vx=i;
        else if (attr=="ffio_y_vel")         i_vy=i;
        else if (attr=="ffio_z_vel")         i_vz=i;
        else if (attr=="m_atomic_number")  { i_anum=i; h->optflags |= MOLFILE_ATOMICNUMBER; }
        else if (attr=="m_chain_name")       i_chain=i;
        else if (attr=="m_pdb_segment_name") i_segid=i;
      }
    }

    virtual void insert_row(const Row &row) {
      molfile_atom_t a;
      memset(&a, 0, sizeof(molfile_atom_t));
      if (i_name>=0)    GET_STR(row[i_name], a.name);
      if (i_name>=0)    GET_STR(row[i_name], a.type);
      if (i_resname>=0) GET_STR(row[i_resname], a.resname);
      if (i_resid>=0)   get_int(row[i_resid], a.resid);
      if (i_segid>=0)   GET_STR(row[i_segid], a.segid);
      if (i_chain>=0)   GET_STR(row[i_chain], a.chain);

      a.insertion[0] = 'A' + m_ct-1;
      if (i_anum>=0)    get_int(row[i_anum], a.atomicnumber);

      // if we didn't get an atom name, try to get one from the atomic number.
      bool bad_name=true;
      for (const char *p=a.name; *p; ++p) {
        if (!isspace(*p)) {
          bad_name = false;
          break;
        }
      }
      if (bad_name && a.atomicnumber>0) {
        strncpy( a.name, 
                 find_element_by_atomic_number(a.atomicnumber).second,
                 sizeof(a.name) );
      }
      
      // if we didn't get a segid, encode the ct index.
      if (!strlen(a.segid)) {
        snprintf(a.segid, 4, "C%d", m_ct);
      }
      atoms.push_back(a);
      natoms += 1;

      pos_t pnt(0,0,0);
      vel_t v(0,0,0);
      if (i_x>=0 && i_y>=0 && i_z>=0) {
        get_float(row[i_x], pnt.x);
        get_float(row[i_y], pnt.y);
        get_float(row[i_z], pnt.z);
      }
      if (i_vx>=0 && i_vy>=0 && i_vz>=0) {
        get_float(row[i_vx], v.x);
        get_float(row[i_vy], v.y);
        get_float(row[i_vz], v.z);
      }
      pos.push_back(pnt);
      vel.push_back(v);
    }
  };

  class SitesArray : public Array {
    int i_mass, i_charge, i_type;
    std::vector<site> &sites;
  public:
    SitesArray(Handle *h_, int ct) 
    : Array(h_, ct), i_mass(-1), i_charge(-1), i_type(-1),
      sites(h->ctmap[ct].sites) 
    {}

    virtual void set_schema( const Schema &schema ) {
      for (unsigned i=0; i<schema.size(); i++) {
        const std::string &attr=schema[i].attr;
        if      (attr=="ffio_mass")   { i_mass=i; h->optflags |= MOLFILE_MASS;}
        else if (attr=="ffio_charge") { i_charge=i; h->optflags |= MOLFILE_CHARGE; }
        else if (attr=="ffio_type")   i_type=i;
      }
    }
    virtual void insert_row(const Row &row) {
      site s;
      if (i_mass>=0) get_float(row[i_mass], s.mass);
      if (i_charge>=0) get_float(row[i_charge], s.charge);
      if (i_type>=0) {
        char type[32];
        GET_STR(row[i_type], type);
        s.pseudo = !strcmp(type, "pseudo");
      }
      sites.push_back(s);
    }
  };
  
  struct BondArray : public Array {
    int i_from, i_to, i_order;
    std::vector<bond_t> &bonds;

  public:
    BondArray(Handle *h_, int ct)
    : Array(h_, ct), i_from(-1), i_to(-1), i_order(-1),
      bonds(h->ctmap[m_ct].bonds)
    {}
    virtual void set_schema( const Schema &schema ) {
      for (unsigned i=0; i<schema.size(); i++) {
        const std::string &attr=schema[i].attr;
        if      (attr=="m_from")  i_from=i;
        else if (attr=="m_to")    i_to=i;
        else if (attr=="m_order") i_order=i;
      }
    }
    virtual void insert_row(const Row &row) {
      if (i_from>=0 && i_to>=0) {
        int from, to, order;
        get_int(row[i_from], from);
        get_int(row[i_to], to);
        if (from < to) {
          if (i_order>=0) get_int(row[i_order], order);
          else            order=1;
          bonds.push_back(bond_t( from, to, order ));
        }
      }
    }
  };

  class VirtualsArray : public Array {
    int i_index, i_ai, i_funct;
    std::string default_funct;

  public:
    VirtualsArray(Handle *h_, int ct, const std::string &def_funct)
    : Array(h_, ct), i_index(-1), i_ai(-1), i_funct(-1),
      default_funct(def_funct)
    {}
    virtual void set_schema( const Schema &schema ) {
      for (unsigned i=0; i<schema.size(); i++) {
        const std::string &attr=schema[i].attr;
        if      (attr=="ffio_index")   i_index=i;
        else if (attr=="ffio_ai")      i_ai=i;
        else if (attr=="ffio_funct")   i_funct=i;
      }
    }
    virtual void insert_row( const Row &row ) {
      if (i_index<0 || i_ai<0) return;
      vsite v;
      int pseudo;
      get_int( row[i_ai], v.ai);
      get_int( row[i_index], pseudo );
      v.funct = (i_funct>=0) ? row[i_funct] : default_funct;
      h->ctmap[m_ct].virtuals[pseudo] = v;
    }
  };

  class PseudoArray : public Array {
    int i_x, i_y, i_z, i_vx, i_vy, i_vz;
    int i_resname, i_chain, i_segid, i_resid;
    std::vector<molfile_atom_t> &atoms;
    std::vector<pos_t> &pos;
    std::vector<vel_t> &vel;
    int &npseudos;

  public:
    PseudoArray(Handle *h_, int ct) 
    : Array(h_, ct),
      i_x(-1), i_y(-1), i_z(-1), 
      i_vx(-1), i_vy(-1), i_vz(-1),
      i_resname(-1), i_chain(-1), i_segid(-1), i_resid(-1),
      atoms( h->ctmap[m_ct].particles ),
      pos( h->ctmap[m_ct].position),
      vel( h->ctmap[m_ct].velocity),
      npseudos( h->ctmap[m_ct].npseudos)
    {}
    virtual void set_schema( const Schema &schema ) {
      for (unsigned i=0; i<schema.size(); i++) {
        const std::string &attr=schema[i].attr;
        if      (attr=="ffio_x_coord")          i_x=i;
        else if (attr=="ffio_y_coord")          i_y=i;
        else if (attr=="ffio_z_coord")          i_z=i;
        else if (attr=="ffio_x_vel")         i_vx=i;
        else if (attr=="ffio_y_vel")         i_vy=i;
        else if (attr=="ffio_z_vel")         i_vz=i;

        else if (attr=="ffio_pdb_residue_name") i_resname=i;
        else if (attr=="ffio_chain_name")       i_chain=i;
        else if (attr=="ffio_pdb_segment_name") i_segid=i;
        else if (attr=="ffio_residue_number")   i_resid=i;
      }
    }

    virtual void insert_row( const Row &row ) {
      molfile_atom_t a;
      memset(&a, 0, sizeof(molfile_atom_t));

      a.insertion[0] = 'A' + m_ct-1;
      strcpy(a.name, "pseudo");
      strcpy(a.type, "pseudo");
      if (i_resname>=0) GET_STR(row[i_resname], a.resname);
      if (i_chain>=0) GET_STR(row[i_chain], a.chain);
      if (i_segid>=0) GET_STR(row[i_segid], a.segid);
      if (i_resid>=0) get_int(row[i_resid], a.resid);

      atoms.push_back(a);
      npseudos += 1;

      pos_t p(0,0,0);
      vel_t v(0,0,0);
      if (i_x>=0 && i_y>=0 && i_z>=0) {
        get_float(row[i_x], p.x);
        get_float(row[i_y], p.y);
        get_float(row[i_z], p.z);
      }
      if (i_vx>=0 && i_vy>=0 && i_vz>=0) {
        get_float(row[i_vx], v.x);
        get_float(row[i_vy], v.y);
        get_float(row[i_vz], v.z);
      }
      pos.push_back(p);
      vel.push_back(v);
    }
  };

  class FepioArray : public Array {
    std::string m_name;
    int i_ai, i_aj;
  public:
    FepioArray(Handle *h_, int ct, const std::string &name) 
    : Array(h_, ct), m_name(name), i_ai(-1), i_aj(-1) {}

    void set_schema( const Schema &schema ) {
      for (unsigned i=0; i<schema.size(); i++) {
        const std::string &attr=schema[i].attr;
        if      (attr=="fepio_ai") i_ai=i;
        else if (attr=="fepio_aj") i_aj=i;
      }
    }
    void insert_row( const Row &row ) {
      if (i_ai<0 || i_aj<0) return;
      fep_elem elem;
      get_int( row[i_ai], elem.ai );
      get_int( row[i_aj], elem.aj );
      h->fepmap[m_name].push_back(elem);
    }
  };

  struct Block {
    Handle *h;
    const std::string m_name;
    const int m_ct;
    bool m_full_system;
    std::vector<Array *> m_arrays;

    Block(Handle *h_, const std::string &name_, int ct) 
    : h(h_), m_name(name_), m_ct(ct), m_full_system(false) {
    }

    virtual ~Block() {
      for (unsigned i=0; i<m_arrays.size(); i++) delete m_arrays[i];
    }

    void set_attrs( AttrMap &attrs) {
      if (m_name=="f_m_ct") {
        // check for full system
        if (attrs["ffio_ct_type"]=="full_system") {
          m_full_system=true;
          return;
        }
        if (attrs.find("chorus_box_ax")!=attrs.end()) {
          h->set_box(attrs);
        }
        if (attrs.find("fepio_stage")!=attrs.end()) {
          int stage = atoi(attrs["fepio_stage"].c_str());
          if      (stage==1) h->stage1 = m_ct;
          else if (stage==2) h->stage2 = m_ct;
        }
      }
    }

    Block new_block(const std::string &name) {
      Block block(h, m_name + "_" + name, m_ct);
      block.m_full_system = m_full_system;
      return block;
    }
    Array& new_array(const std::string &name) {

      // create array subclass based on the name
      Array *arr=NULL;
      if (m_full_system) {
        arr = new Array(h, m_ct);

      } else if (name=="m_atom") {
        arr=new AtomArray(h, m_ct);

      } else if (name=="ffio_pseudo") {
        arr=new PseudoArray(h, m_ct);

      } else if (name=="ffio_virtuals") {
        arr=new VirtualsArray(h, m_ct, "virtual");

      } else if (name=="ffio_polarizable") {
        arr=new VirtualsArray(h, m_ct, "polar");

      } else if (name=="ffio_sites") {
        arr=new SitesArray(h, m_ct);

      } else if (m_name=="f_m_ct_fepio_fep" && name=="fepio_atommaps") {
        arr=new FepioArray(h, m_ct, name);

      } else if (name=="m_bond") {
        arr=new BondArray(h, m_ct);

      } else {
        arr = new Array(h, m_ct);
      }
      m_arrays.push_back(arr);
      return *arr;
    }
  };

  /// parser callbacks

  Schema predict_schema(Tokenizer& tokenizer) {
    Schema schemas;
    while(tokenizer.not_a(":::")) {
      schema_t schema;
      std::string token = tokenizer.token();
      if (token[0] != 'b' && token[0] != 'i' && token[0] != 'r' && token[0] != 's') {
        std::stringstream str;
        str << "Line " << tokenizer.line() << " predicted a schema, but " 
            << token << " didn't start b_ i_ r_ or s_ ";
        throw std::runtime_error(str.str());
      }
      schema.type = token[0];
      schema.attr = token.substr(2);
      //schema.doc  = tokenizer.optional_comment();
      schemas.push_back(schema);
      tokenizer.next();
    }
    return schemas;
  }

  void
  predict_schema_and_values(Block& M, Tokenizer& tokenizer) {
    Schema schema = predict_schema(tokenizer);
    AttrMap attrs;
    tokenizer.predict(":::");
    for (unsigned i=0;i<schema.size();++i) {
      std::string value = tokenizer.predict();
      if (value == "<>" || value == "") continue; // use default element
      // Strip quotes if present
      if (value[0] == '"' && value[value.size()-1]) {
        value = value.substr(1,value.size()-2);
      }
      attrs[schema[i].attr] = value;
    }
    M.set_attrs(attrs);
  }

  void
  check_name(const Tokenizer& tokenizer,const std::string& name) {
    if (name.size() > 0 && !(isalpha(name[0]) || name[0] == '_')) {
      std::stringstream str;
      str << "Line " << tokenizer.line() << " predicted a block name have " << name << std::endl;
      throw std::runtime_error(str.str());
    }
  }

  // forward declaration
  void predict_block(Block& M, Tokenizer& tokenizer);

  void predict_blockbody(Block& subblock, Tokenizer& tokenizer) {
    tokenizer.predict("{");
    predict_schema_and_values(subblock,tokenizer);
    while(tokenizer.not_a("}")) {
      predict_block(subblock,tokenizer);
    }
    tokenizer.predict("}");
  }

  void predict_arraybody(Array& subarray, Tokenizer& tokenizer) {
  
    // Read header
    tokenizer.predict("[");
    tokenizer.predict();
    tokenizer.predict("]");
    tokenizer.predict("{");
  
    // Read schema
    Schema schema = predict_schema(tokenizer);
    subarray.set_schema(schema);
    size_t width = schema.size();
    Row row(width);
    tokenizer.predict(":::");
  
    // Read rows
    while(tokenizer.not_a(":::")) {
      // throw away row index
      tokenizer.predict();
      for(unsigned i=0;i<width;++i) {
        row[i] = tokenizer.predict();
      }
      subarray.insert_row(row);
    }
  
    tokenizer.predict(":::");
  
    tokenizer.predict("}");
  }
  

  void predict_nameless_block(std::string name,Block& M, Tokenizer& tokenizer) {
    // -----------------------------------------------
    // May be an array
    // -----------------------------------------------
    std::string tok = tokenizer.token();
    if (tok == "[") {
      predict_arraybody(M.new_array(name),tokenizer);
    } 
  
    // -----------------------------------------------
    // Otherwise just a block
    // -----------------------------------------------
    else {
      Block subblock = M.new_block(name);
      predict_blockbody(subblock,tokenizer);
    }
  }

/**************************************************************************/
/* LOCAL  **************         check_name        ************************/


  void predict_block(Block& M, Tokenizer& tokenizer) {
    std::string name = tokenizer.predict();
    check_name(tokenizer,name);
    predict_nameless_block(name,M,tokenizer);
  }

  void fill_nameless( Block &block, Tokenizer& tokenizer) {
    predict_blockbody(block,tokenizer);
  }

}


namespace {

  std::string quotify( const std::string &s ) {
    // empty string --> quoted ""
    if (s == "") return "\"\"";

    std::string raw(s);

    // Check for non-printable characters and "
    for(std::string::iterator p=raw.begin(), en=raw.end();
        p != en; ++p) {
      if (isspace(*p) || !isprint(*p) || *p == '"' || *p == '<' || *p == '\\') {
        std::string escaped(raw.begin(),p);
        for(;p!=en;++p) {
          // We only support space and tab
          if (isspace(*p) && !(*p == ' ' || *p == '\t')) {
            throw std::invalid_argument("unprintable whitespace in '" + raw + '\'');
          }

          if (*p == '"') {
            escaped += "\\\"";
          } else if (*p == '\\') {
            escaped += "\\\\";
          } else {
            escaped += *p;
          }

        }
        raw = '"' + escaped + '"';
        break;
      }
    }
    return raw;
  }

  void write_meta(std::ofstream &output) {
    output << "{" << std::endl
           << "  s_m_m2io_version" << std::endl
           << "  :::" << std::endl
           << "  2.0.0" << std::endl
           << "}" << std::endl;
  }
  void write_ct_header( std::ofstream &output, 
                        const double *A, const double *B, const double *C ) {
    output << "f_m_ct {\n"
           << "  r_chorus_box_ax\n"
           << "  r_chorus_box_ay\n"
           << "  r_chorus_box_az\n"
           << "  r_chorus_box_bx\n"
           << "  r_chorus_box_by\n"
           << "  r_chorus_box_bz\n"
           << "  r_chorus_box_cx\n"
           << "  r_chorus_box_cy\n"
           << "  r_chorus_box_cz\n"
           << "  :::\n";
    int i;
    for (i=0; i<3; i++) output << "  " << A[i] << std::endl;
    for (i=0; i<3; i++) output << "  " << B[i] << std::endl;
    for (i=0; i<3; i++) output << "  " << C[i] << std::endl;
  }
  void write_ct_atoms(  std::ofstream &output,
                        const std::map<size_t,int>& atommap,
                        const std::vector<molfile_atom_t> &atoms,
                        const float *pos, const float *vel ) {

    output << "  m_atom[" << atommap.size() << "] {\n";
    output << "    s_m_pdb_atom_name\n";
    output << "    s_m_pdb_residue_name\n";
    output << "    s_m_chain_name\n";
    output << "    s_m_pdb_segment_name\n";
    output << "    i_m_residue_number\n";
    output << "    r_m_x_coord\n";
    output << "    r_m_y_coord\n";
    output << "    r_m_z_coord\n";
    if (vel) {
      output << "    r_ffio_x_vel\n";
      output << "    r_ffio_y_vel\n";
      output << "    r_ffio_z_vel\n";
    }
    output << "    i_m_atomic_number\n";
    output << "    i_m_mmod_type\n";
    output << "    i_m_color\n";
    output << "    i_m_visibility\n";
    output << "    i_m_formal_charge\n";
    output << "    r_m_charge1\n";
    output << "    r_m_charge2\n";
    output << "    s_m_mmod_res\n";
    output << "    s_m_grow_name\n";
    output << "    s_m_insertion_code\n";
    output << "    :::\n";

    for (std::map<size_t,int>::const_iterator i=atommap.begin(); 
        i!=atommap.end(); ++i) {
      const molfile_atom_t &a = atoms[i->first];
      output << "    " << i->second << ' '
             << quotify(a.name) << ' '
             << quotify(a.resname) << ' '
             << quotify(a.chain) << ' '
             << quotify(a.segid) << ' '
             << a.resid << ' '
             << pos[0+3*i->first] << ' '
             << pos[1+3*i->first] << ' '
             << pos[2+3*i->first] << ' ';
      if (vel) {
        output << vel[0+3*i->first] << ' '
               << vel[1+3*i->first] << ' '
               << vel[2+3*i->first] << ' ';
      }
      // get atomic number from mass if not provided explicitly.
      int anum = a.atomicnumber;
      if (anum < 1) anum = find_element_by_amu(a.mass).first;
      // the other setting make Maestro happy
      int color=2; // gray
      int mmod=64; // mmod_type; 64="any atom"
      switch  (anum) {
        case 1:  color=21; mmod=48; break;  // H
        case 3:  color=4;  mmod=11; break;  // Li+ ion
        case 6:  color=2 ; mmod=14; break;  // C
        case 7:  color=43; mmod=40; break;  // N
        case 8:  color=70; mmod=23; break;  // O
        case 9:  color=8;  mmod=56; break;  // F
        case 11: color=4;  mmod=66; break;  // Na+ ion
        case 12: color=4;  mmod=72; break;  // Mg2+ ion
        case 14: color=14; mmod=60; break;  // Si
        case 15: color=15; mmod=53; break;  // P
        case 16: color=13; mmod=52; break;  // S
        case 17: color=13; mmod=102; break;  // Cl- ion
        case 19: color=4;  mmod=67; break;  // K+ ion
        case 20: color=4;  mmod=70; break;  // Ca2+ ion
        default: ;
      }
      static const std::string blank("\" \"");
      output << anum << ' '
             << mmod << ' '       // mmod_type
             << color << ' '   // m_color
             << 1 << ' '       // m_visibility
             << 0   << ' '     // formal charge
             << 0.0 << ' '     // charge1
             << 0.0 << ' '     // charge2
             << blank << ' '  // mmod_res
             << blank << ' '  // m_grow_name
             << quotify(a.insertion) << ' ' // m_insertion_code
             << std::endl;    
    }
    output << "    :::\n";
    output << "  }\n";
  }

  void write_ct_bonds( std::ofstream &output, 
                       const std::vector<bond_t> &bonds ) {

    // don't write 0-element m_bond (ev85392)
    if (!bonds.size()) return;
    output << "  m_bond[" << bonds.size() << "] {\n"
           << "    i_m_from\n"
           << "    i_m_to\n"
           << "    i_m_order\n"
           << "    :::\n";
    for (unsigned i=0; i<bonds.size(); i++) {
      output << "    " 
             << i+1 << ' '
             << bonds[i].from << ' '
             << bonds[i].to << ' '
             << (int)floorf(0.5+bonds[i].order) << "\n";
    }
    output << "    :::\n"
           << "  }\n";
  }

  void write_ct_ffio_header( std::ofstream &output ) {
    output << "  ffio_ff {\n"
           << "    :::\n";
  }
  void write_ct_ffio_footer( std::ofstream &output ) {
    output << "  }\n";
  }

  void write_ct_sites( std::ofstream &output, 
                       const std::vector<site>& sites ) {

    output << "    ffio_sites[" << sites.size() << "] {\n"
           << "      s_ffio_type\n"
           << "      r_ffio_charge\n"
           << "      r_ffio_mass\n"
           << "      :::\n";
    for (size_t i=0; i<sites.size(); i++) {
      output << "      " << i+1 << ' '
             << (sites[i].pseudo ? "pseudo " : "atom ")
             << sites[i].charge << ' '
             << sites[i].mass << "\n";
    }
    output << "      :::\n";
    output << "    }\n";
  }

  void write_ct_pseudos( std::ofstream &output,
                         const std::map<size_t,int> &pseudos,
                         const std::vector<molfile_atom_t> &particles,
                         const float *pos, const float *vel ) {
    if (!pseudos.size()) return;
    output << "    ffio_pseudo[" << pseudos.size() << "] {\n"
           << "      r_ffio_x_coord\n"
           << "      r_ffio_y_coord\n"
           << "      r_ffio_z_coord\n"
           << "      s_ffio_pdb_residue_name\n"
           << "      s_ffio_chain_name\n"
           << "      s_ffio_pdb_segment_name\n"
           << "      i_ffio_residue_number\n";
    if (vel) output << "      r_ffio_x_vel\n"
                    << "      r_ffio_y_vel\n"
                    << "      r_ffio_z_vel\n";
    output << "      :::\n";

    for (std::map<size_t,int>::const_iterator i=pseudos.begin(); 
        i!=pseudos.end(); ++i) {
      const molfile_atom_t &a = particles[i->first];
      output << "      " 
             << i->second << ' '
             << pos[0+3*i->first] << ' ' 
             << pos[1+3*i->first] << ' ' 
             << pos[2+3*i->first] << ' ' 
             << quotify(a.name) << ' '
             << quotify(a.chain) << ' '
             << quotify(a.segid) << ' '
             << a.resid;
      if (vel) { 
        output << ' ' 
               << vel[0+3*i->first] << ' ' 
               << vel[1+3*i->first] << ' ' 
               << vel[2+3*i->first];
      }
      output << "\n";
    }
    output << "      :::\n";
    output << "    }\n";
  }

  void write_ct_footer( std::ofstream &output ) {
    output << "}\n";
  }
}

////// support for alchemical mae files

namespace {
#if 1
  /*! Map m_bond records from stage2 into stage1
   *  @param a2inv mapping from ct2 atom index to 
   *  @param ct1 stage 1 ct
   *  @param ct2 stage 2 ct
   */
  void fixup_m_bond(const std::map<int,int>& a2inv, 
                    ct_data &ct1, ct_data &ct2) {

    // go through every m_bond record in ct2 and add it to ct1 after
    // mapping the atom indices.  
    if (!ct2.bonds.size()) return; // nothing to do

    // Keep track of the bonds we already have
    typedef std::set<std::pair<int,int> > BondSet;
    BondSet bondset;
    unsigned i;
    for (i=0; i<ct1.bonds.size(); i++) {
      const bond_t &b = ct1.bonds[i];
      bondset.insert(std::make_pair(b.from,b.to));
    }

    // add the others
    for (i=0; i<ct2.bonds.size(); i++) {
      const bond_t &b = ct2.bonds[i];
      int from = b.from;
      int to   = b.to;
      std::map<int,int>::const_iterator from_iter = a2inv.find(from);
      std::map<int,int>::const_iterator to_iter   = a2inv.find(to);
      if (from_iter == a2inv.end() ||
          to_iter == a2inv.end())  {
        fprintf(stderr, "Missing entry in fepio_atommap for %d %d\n", from, to );
        throw std::runtime_error("Bad fepio_atommap");
      }
      BondSet::value_type p(from_iter->second, to_iter->second);
      if (bondset.find(p) != bondset.end()) continue;
      bondset.insert(p);
      ct1.bonds.push_back( bond_t( p.first, p.second, 1 ));
    }
  }
#endif

  /*! append alchemical forms of atoms and groups onto mapped original forms
   * @param fm1 stage 1 ct
   * @param fm2 stage 2 ct
   * @param map fepio_fep map
   */
  void alchemical_combine( Handle *h ) {
    // check if alchemical stages were assigned
    if (h->stage1<1 || h->stage2<1) return;
    fprintf(stderr, "alchemical system detected\n");
    ct_data &stage1 = h->ctmap[h->stage1];
    ct_data &stage2 = h->ctmap[h->stage2];

    // Map perturbed ct atom numbers into combined ct numbers
    std::map<int,int> a2_inv_map;
    for (int i=1; i<=stage2.natoms; i++) a2_inv_map[i]=i;

    // find the atom mapping
    FepioMapping::const_iterator atoms = h->fepmap.find("fepio_atommaps");
    if (atoms != h->fepmap.end()) {
      for (FepList::const_iterator i=atoms->second.begin();
          i!=atoms->second.end(); ++i) {
        int const ai = i->ai;
        int const aj = i->aj;
        if (ai > 0 && aj > 0) {
          a2_inv_map[aj] = ai;
        } else if ( ai > 0 && aj < 0 ) {
          // nothing to do
	} else if ( ai < 0 && aj > 0 ) {
	  a2_inv_map[aj] = ::abs(ai);

          // copy atom to stage1
          molfile_atom_t &atom = stage2.particles.at( aj-1 );
          atom.insertion[0] = 'A' + h->stage1-1;
          stage1.particles.push_back(atom);
          stage1.natoms += 1;

          // copy site as well
          stage1.sites.push_back( stage2.sites.at(aj-1));
          stage1.sites[stage1.sites.size()-1].charge = 0.0;
          
          // copy position and velocity
          stage1.position.push_back( stage2.position.at(aj-1) );
          stage1.velocity.push_back( stage2.velocity.at(aj-1) );

	} else {
          fprintf(stderr, "ai(%d) and aj(%d) < 0 in atommap\n", ai, aj);
        }
      }
    }
    fixup_m_bond(a2_inv_map, stage1, stage2);
    h->ctmap.erase(h->stage2);
  }
}

namespace {

  void *open_file_read( const char *fname, const char *ftype, int *vmdatoms) {

    std::ifstream in(fname, std::ifstream::in | std::ifstream::binary );
    if (!in) return NULL;

    Handle *h = new Handle;
    *vmdatoms = 0;

    try {
      // parse the mae file contents
      Tokenizer tokenizer(in);

      Block meta(h, "meta", 0);
      fill_nameless( meta, tokenizer);
      int ct=1;
      while (tokenizer.not_a()) {
        std::string name = tokenizer.predict();
        Block block(h, name, ct++);
        fill_nameless( block, tokenizer );
      }
      alchemical_combine(h);
    }
    catch (std::exception &e) {
      fprintf(stderr, "Reading mae file failed: %s\n", e.what());
      delete h;
      return NULL;
    }

    // post-processing once all arrays have been read
    for (CtMap::const_iterator i=h->ctmap.begin(); i!=h->ctmap.end(); ++i) {
      int natoms = i->second.natoms;
      int npseudos = i->second.npseudos;
      int nparticles = natoms + npseudos;
      int nsites = i->second.sites.size();

      *vmdatoms += nparticles;

      if (nsites>0) {
        // check particle/site consistency
        if (nsites > nparticles ) {
          fprintf(stderr, "ERROR: Too many ffio_sites records in ct %d\n", 
              i->first);
          delete h;
          return NULL;
        }
        int nblks = nparticles / nsites;
        int Np = npseudos / nblks;
        int Na = natoms / nblks;
        if (Np + Na != nsites) {
          fprintf(stderr, "ERROR: Number of particles in ct %d not a multiple of the number of ffio_sites\n", i->first);
          delete h;
          return NULL;
        }
      }

    }
    return h;
  }

  int read_structure(void* v, int* optflags, molfile_atom_t *atoms) {
    Handle *h = reinterpret_cast<Handle *>(v);

    // apply sites information
    for (CtMap::iterator i=h->ctmap.begin(); i!=h->ctmap.end(); ++i) {
      ct_data &ct = i->second;
      int natoms = ct.natoms;
      int npseudos = ct.npseudos;
      int nsites = ct.sites.size();
      int nparticles = natoms + npseudos;

      if (nsites) {

        int nblks = nparticles / nsites;
        int Np = npseudos / nblks;
        int Na = natoms / nblks;
        int atom_number = 0;
        int pseudo_number = natoms;

        // mapping from index in site block to relative order in site block
        std::map<int,int> pseudo_index;

        for (int j=0; j<nsites; ++j) {
          float mass = ct.sites[j].mass;
          float charge = ct.sites[j].charge;

          if (!ct.sites[j].pseudo) {
            // regular atom
            for(int k=0; k<nblks; ++k) {
              int ia = atom_number + k*Na;
              ct.particles[ia].mass = mass;
              ct.particles[ia].charge = charge;
            }
            ++atom_number;

          } else {
            // pseudo
            int n = pseudo_index.size();
            pseudo_index[j+1] = n;
            for(int k=0; k<nblks; ++k) {
              int ip = pseudo_number + k*Np;
              ct.particles[ip].mass = mass;
              ct.particles[ip].charge = charge;
            }
            ++pseudo_number;
          }
        }

        // check ffio_virtuals for pseudobonds
        for (VirtualsMap::const_iterator v=ct.virtuals.begin();
            v!=ct.virtuals.end(); ++v) {
          int index = v->first;
          int ai = v->second.ai;
          if (ai && index) {
            atom_number = ai - 1;
            pseudo_number = natoms + pseudo_index[index];

            for (int k=0; k<nblks; k++) {
              int from = atom_number + k*Na;
              int to = pseudo_number + k*Np;

              // copy residue information from parent to the pseudo
              strcpy( ct.particles[to].resname, ct.particles[from].resname );
              strcpy( ct.particles[to].chain,   ct.particles[from].chain);
              strcpy( ct.particles[to].segid,   ct.particles[from].segid);
              ct.particles[to].resid = ct.particles[from].resid;
              strncpy( ct.particles[to].type, v->second.funct.c_str(), 
                       sizeof(ct.particles[to].type ) );

              // add a pseudobond between the pseudo and the parent atom
              ct.bonds.push_back(bond_t( from+1, to+1, 1 ));
            }
          }
        }
      }
      memcpy( atoms, &ct.particles[0], nparticles*sizeof(molfile_atom_t) );
      atoms += nparticles;
    }

    *optflags = h->optflags;
    return MOLFILE_SUCCESS;
  }

#if vmdplugin_ABIVERSION > 14
  int read_bonds(void *v, int *nbonds, int **from, int **to, float** order,
                 int **bondtype, int *nbondtypes, char ***bondtypename) {
#else
  int read_bonds(void *v, int *nbonds, int **from, int **to, float** order) {
#endif
    Handle *h = reinterpret_cast<Handle *>(v);

    // map bonds from ct arrays to from/to/order arrays
    int offset=0;
    for (CtMap::const_iterator i=h->ctmap.begin(); i!=h->ctmap.end(); ++i) {
      const ct_data &ct = i->second;
      for (std::vector<bond_t>::const_iterator b=i->second.bonds.begin(),
           e=i->second.bonds.end(); b!=e; ++b) {
        h->bond_from.push_back( b->from + offset );
        h->bond_to.push_back( b->to + offset );
        h->bond_order.push_back( b->order );
      }
      offset += ct.particles.size();
    }
    *nbonds = h->bond_from.size();
    *from = &h->bond_from[0];
    *to = &h->bond_to[0];
    *order = &h->bond_order[0];
#if vmdplugin_ABIVERSION > 14
    *bondtype = NULL;
    *nbondtypes = 0;
    *bondtypename = NULL;
#endif
    return MOLFILE_SUCCESS;
  }

  int read_timestep_metadata(void *v, molfile_timestep_metadata *m) {
    m->has_velocities = true;
    m->count = 1;
    return MOLFILE_SUCCESS;
  }

  int read_next_timestep(void *v, int natoms, molfile_timestep_t *ts) {
    Handle *h = reinterpret_cast<Handle *>(v);
    if (h->eof) return MOLFILE_EOF;

    float *pos = ts->coords;
    float *vel = ts->velocities;
    for (CtMap::const_iterator i=h->ctmap.begin(); i!=h->ctmap.end(); ++i) {
      const ct_data &ct = i->second;
      unsigned nparticles = ct.position.size();
      memcpy( pos, &ct.position[0], 3*nparticles*sizeof(float) );
      pos += 3*nparticles;
      if (vel) {
        memcpy( vel, &ct.velocity[0], 3*nparticles*sizeof(float) );
        vel += 3*nparticles;
      }
    }

    h->get_box(ts);
    h->eof = true;
    return MOLFILE_SUCCESS;
  }

  void close_file_read( void *v) {
    delete reinterpret_cast<Handle *>(v);
  }

  void *open_file_write(const char *path, const char *type, int natoms) {
    Handle *h = new Handle;
    h->output.open(path);
    if (!h->output) {
      fprintf(stderr, "Could not open '%s' for writing.\n", path);
      delete h;
      return NULL;
    }
    h->nparticles = natoms;
    h->particles.resize(natoms);
    return h;
  }

  int write_structure(void *v, int optflags, const molfile_atom_t *atoms) {
    Handle *h = reinterpret_cast<Handle *>(v);
    h->optflags = optflags;
    memcpy(&h->particles[0], atoms, h->particles.size()*sizeof(molfile_atom_t));

    // assign ct for each particle, and count ct atoms and pseudos
    std::vector<int> atom_ct(h->nparticles);
    char last_insertion = 0;
    int ct = 1;
    for (int i=0; i<h->nparticles; i++) {
      const molfile_atom_t &a = atoms[i];
      if (i==0) last_insertion = a.insertion[0];
      else if (a.insertion[0] != last_insertion) {
        last_insertion = a.insertion[0];
        ++ct;
      }
      atom_ct[i] = ct;
      ct_data &data = h->ctmap[ct];
      site s;
      s.charge = a.charge;
      s.mass = a.mass;

      if ((optflags & MOLFILE_ATOMICNUMBER) && a.atomicnumber<1) {
        data.pseudomap[i] = ++data.npseudos;
        s.pseudo = true;
      } else {
        data.atommap[i] = ++data.natoms;
        s.pseudo = false;
      }
      data.sites.push_back(s);
    }

    // add bonds to their proper ct
    int badbonds=0, skipped=0;
    for (unsigned b=0; b<h->bond_from.size(); b++) {
      int from = h->bond_from[b]-1;
      int to   = h->bond_to[b]-1;
      float order = h->bond_order[b];
      if (from > to) continue;
      int ct = atom_ct[from];
      if (ct != atom_ct[to]) {
        ++badbonds;
        continue;
      }
      // map 0-based bond index to 1-based ct index
      ct_data &data = h->ctmap[ct];
      // add the bond only if it's between atoms, not pseudos
      std::map<size_t,int>::const_iterator ifrom=data.atommap.find(from);
      std::map<size_t,int>::const_iterator ito=data.atommap.find(to);
      if (ifrom != data.atommap.end() && ito != data.atommap.end())
        data.bonds.push_back(bond_t(ifrom->second, ito->second, order));
      else 
        ++skipped;
    }
    if (badbonds) {
      fprintf(stderr, "Could not store all bonds in mae file\n");
      fprintf(stderr, "Check that no bonded atoms have different insertion\n");
      return MOLFILE_ERROR;
    }
    if (skipped) {
      fprintf(stderr, "Info) Skipped %d pseudobonds.\n", skipped);
    }

    return MOLFILE_SUCCESS;
  }

#if vmdplugin_ABIVERSION > 14
  int write_bonds(void *v, int nbonds, int *from, int *to, float *order,
                  int *bondtype, int nbondtypes, char **bondtypename) {
#else
  int write_bonds(void *v, int nbonds, int *from, int *to, float *order) {
#endif
    Handle *h = reinterpret_cast<Handle *>(v);
    h->bond_from.resize(nbonds);
    h->bond_to.resize(nbonds);
    h->bond_order.resize(nbonds);
    memcpy( &h->bond_from[0], from, nbonds*sizeof(int) );
    memcpy( &h->bond_to[0],   to,   nbonds*sizeof(int) );
    for (int i=0; i<nbonds; i++) h->bond_order[i] = order ? order[i] : 1;

    return MOLFILE_SUCCESS;
  }

  int write_timestep(void *v, const molfile_timestep_t *ts) {
    Handle *h = reinterpret_cast<Handle *>(v);
    if (h->eof) {
      fprintf(stderr, "Cannot write multiple frames to mae file\n");
      return MOLFILE_EOF;
    }
    try {
      h->set_box(ts);

      write_meta(h->output);
      for (CtMap::const_iterator i=h->ctmap.begin(); i!=h->ctmap.end(); ++i) {
        const ct_data &data = i->second;
        write_ct_header( h->output, h->A, h->B, h->C );
        write_ct_atoms(  h->output, data.atommap, h->particles, ts->coords, ts->velocities );
        write_ct_bonds(  h->output, data.bonds );

        write_ct_ffio_header( h->output );
          write_ct_sites( h->output, data.sites );
          write_ct_pseudos( h->output, data.pseudomap, h->particles, ts->coords, ts->velocities );
        write_ct_ffio_footer( h->output );

        write_ct_footer( h->output );
      }
    }
    catch (std::exception &e) {
      fprintf(stderr, e.what());
      return MOLFILE_ERROR;
    }
    return MOLFILE_SUCCESS;
  }

  void close_file_write( void *v ) {
    Handle *h = reinterpret_cast<Handle *>(v);
    h->output.close();
    delete h;
  }
}

///////////////////////////////////////////////////////////////////
//
// Plugin Interface
//
// ////////////////////////////////////////////////////////////////

static molfile_plugin_t maeff;

VMDPLUGIN_EXTERN int VMDPLUGIN_init (void) {
  /* Plugin for maeff trajectory files */
  ::memset(&maeff,0,sizeof(maeff));
  maeff.abiversion = vmdplugin_ABIVERSION;
  maeff.type = MOLFILE_PLUGIN_TYPE;
  maeff.name = "mae";
  maeff.prettyname = "Maestro File";
  maeff.author = "D. E. Shaw Research";
  maeff.majorv = 3;
  maeff.minorv = 3;
  maeff.is_reentrant = VMDPLUGIN_THREADUNSAFE;

  maeff.filename_extension = "mae,maeff,cms";
  maeff.open_file_read = open_file_read;
  maeff.read_structure = read_structure;
  maeff.read_bonds = read_bonds;
  maeff.read_timestep_metadata = read_timestep_metadata;
  maeff.read_next_timestep = read_next_timestep;
  maeff.close_file_read = close_file_read;

  maeff.open_file_write = open_file_write;
  maeff.write_structure = write_structure;
  maeff.write_bonds = write_bonds;
  maeff.write_timestep = write_timestep;
  maeff.close_file_write = close_file_write;

  return VMDPLUGIN_SUCCESS;
}

VMDPLUGIN_EXTERN int VMDPLUGIN_register(void *v, vmdplugin_register_cb cb) {
  cb(v,reinterpret_cast<vmdplugin_t*>(&maeff));
  return VMDPLUGIN_SUCCESS;
}

VMDPLUGIN_EXTERN int VMDPLUGIN_fini(void) {
  return VMDPLUGIN_SUCCESS;
}

#if defined(TEST_MAEFFPLUGIN)

int main(int argc, char *argv[]) {
  if (argc != 3) {
    fprintf(stderr, "Usage: %s input.mae output.mae\n", argv[0]);
    return 1;
  }
  int natoms=0;
  void *v = open_file_read( argv[1], "mae", &natoms);
  if (!v) return 1;
  // read an existing mae file
  printf("%d atoms\n", natoms);
  std::vector<molfile_atom_t> atoms(natoms);
  memset(&atoms[0], 0, natoms*sizeof(molfile_atom_t));
  int optflags=0;
  if (read_structure(v, &optflags, &atoms[0])!=MOLFILE_SUCCESS) return 1;
  int nbonds;

  // read bonds
  int *from, *to;
  float *order;
  if (read_bonds(v, &nbonds, &from, &to, &order) !=MOLFILE_SUCCESS) return 1;
  std::vector<int> bond_from(nbonds), bond_to(nbonds);
  std::vector<float> bond_order(nbonds);
  memcpy(&bond_from[0], from, nbonds*sizeof(int));
  memcpy(&bond_to[0], to, nbonds*sizeof(int));
  memcpy(&bond_order[0], order, nbonds*sizeof(float));

  // read coordinates
  molfile_timestep_t ts[1];
  ts->coords = new float[3*natoms];
  ts->velocities = new float[3*natoms];
  if (read_next_timestep(v, natoms, ts) != MOLFILE_SUCCESS) return 1;

  // close
  close_file_read(v);

  // write file
  v = open_file_write( argv[2], "mae", natoms );
  if (!v) return 1;
  if (write_bonds(v, nbonds, &bond_from[0], &bond_to[0], &bond_order[0]) 
      != MOLFILE_SUCCESS) return 1;
  if (write_structure( v, optflags, &atoms[0] ) != MOLFILE_SUCCESS) return 1;
  if (write_timestep( v, ts ) != MOLFILE_SUCCESS ) return 1;
  close_file_write(v);

  delete [] ts->coords;
  delete [] ts->velocities;

  // make sure we can still parse the file we wrote
  int new_natoms;
  v = open_file_read( argv[2], "mae", &new_natoms );
  if (!v) return 1;
  if (new_natoms != natoms) {
    fprintf(stderr, "number of atoms changed: %d -> %d\n", natoms, new_natoms);
    return 1;
  }
  close_file_read(v);
  return 0;
}

#endif
