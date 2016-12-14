/* MACHINE GENERATED FILE, DO NOT EDIT! */

#define VMDPLUGIN molfile_pdbxplugin
#define STATIC_PLUGIN 1

/***************************************************************************
 *cr
 *cr            (C) Copyright 1995-2016 The Board of Trustees of the
 *cr                        University of Illinois
 *cr                         All Rights Reserved
 *cr
 ***************************************************************************/
/***************************************************************************
 * RCS INFORMATION:
 *
 *      $RCSfile: pdbxplugin.C,v $
 *      $Author: johns $       $Locker:  $             $State: Exp $
 *      $Revision: 1.17 $       $Date: 2016/11/28 05:01:22 $
 *
 ***************************************************************************/

#include "largefiles.h"   /* platform dependent 64-bit file I/O defines */

#include "molfile_plugin.h"
#include "periodic_table.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#if !(defined(WIN32) || defined(WIN64))
#include <sys/time.h>
#endif

#include "inthash.h"

//#define PDBX_DEBUG 1

//used for reading author defined values needed when reading bond info
#define CHAIN_SIZE 4
#define TYPE_SIZE 8
#define BUFFER_SIZE 1024

struct list_node { 
  unsigned int next;
  unsigned int index;
};

// class definition
typedef struct pdbxParser {
  FILE *file;
  int natoms;
  int nbonds;
  int* resid_auth;
  char * chain_auth;
  char * type_auth;
  float* xyz;
  int* bondsTo;
  int* bondsFrom;
  bool error;
  int table[64];
  unsigned int tableSize;
  inthash_t bondHash;
  list_node * hashMem;
} pdbxParser;

// XXX yuck, this needs to go!
static unsigned char charToNum[128];

enum TableColums {
  COLUMN_NUMBER,
  COLUMN_NAME,
  COLUMN_TYPE,
  COLUMN_TYPE_AUTH,
  COLUMN_RESNAME,
  COLUMN_RESID,
  COLUMN_RESID_AUTH,
  COLUMN_INSERTION,
  COLUMN_X,
  COLUMN_Y,
  COLUMN_Z,
  COLUMN_OCCUPANCY,
  COLUMN_BFACTOR,
  COLUMN_CHARGE,
  COLUMN_CHAIN,
  COLUMN_CHAIN_AUTH,
  COLUMN_JUNK
};

/* Opens the file, finds the number of atoms, and allocates arrays */
/* Reads in and stores information from the file */
static pdbxParser* create_pdbxParser(const char* filepath);

/*Reads through the file and stores data */
static int parseStructureFaster(molfile_atom_t * atoms, int * optflags, pdbxParser* parser);

static int setCoordinatesFast(float* coords, pdbxParser* parser);

static bool readRMSDBonds(molfile_atom_t * atoms, pdbxParser* parser);
static bool readAngleBonds(molfile_atom_t * atoms, pdbxParser* parser);
static bool readBonds(molfile_atom_t * atoms, pdbxParser* parser);

/* Parse through file and return the total number of atoms */
/* Will rewind the file to the start */
/* Returns -1 if the number of atoms cannot be found */
static int parseNumberAtoms(pdbxParser* parser);

/* returns true if str starts with "_atom_site." */
static inline bool isAtomSite(char * str);

static inline bool isValidateRMSDBond(char * str);

/* Assumes that str contains a single floating point number and */
/* returns it as a float. NO ERROR CHECKING */
/* Must be passed a null terminating string */
/* Wrote specifically to parse strings returned from getNextWord */
static float stringToFloat(char * str);

/* Takes a string str and finds the next word starting from pos*/
/* word must be allocated and suffiently large, does NO ERROR CHECKING */
/* After returning, word will contain the next word and pos will be updated */
/* to point to the current position in str */
static void getNextWord(char * str, void * word, int& pos);

/* Takes a string str and finds the next word starting from pos*/
/* word must be allocated and suffiently large, does NO ERROR CHECKING */
/* After returning, word will contain the next word and pos will be updated */
/* to point to the current position in str */
static void skipNextWord(char * str, void * word, int& pos);

/* Returns a unique int id for an atom based on the chain and resid */
static inline int getUniqueResID(char * chainstr, int resid);

static void initCharToNum();

#define WB_SIZE 1024

#if 0 
static const char atomSiteHeader[] =
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
  "_atom_site.Cartn_x_esd\n"
  "_atom_site.Cartn_y_esd\n"
  "_atom_site.Cartn_z_esd\n"
  "_atom_site.occupancy_esd\n"
  "_atom_site.B_iso_or_equiv_esd\n"
  "_atom_site.pdbx_formal_charge\n"
  "_atom_site.auth_seq_id\n"
  "_atom_site.auth_comp_id\n"
  "_atom_site.auth_asym_id\n"
  "_atom_site.auth_atom_id\n"
  "_atom_site.pdbx_PDB_model_num\n";
#endif


static const char atomSiteHeader[] =
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
  "_atom_site.pdbx_formal_charge\n"
  "_atom_site.auth_asym_id\n";


typedef struct pdbxWriter {
  FILE* fd;
  char writeBuf[WB_SIZE];
  char pdbName[256];
  int bufferCount;
  molfile_atom_t* atoms;
  const float* coordinates;
  int numatoms;
} pdbxWriter;

static void writeBuffer(pdbxWriter* writer);
static void writeIntro(pdbxWriter* writer);
static void write(const char* str, pdbxWriter* writer);
static void writeAtomSite(pdbxWriter* writer);
static void close(pdbxWriter* writer);
static pdbxWriter* create_pdbxWriter(const char* filename, int numAtoms);
static void addAtoms(const molfile_atom_t* atoms, int optflags, pdbxWriter* writer);
static void addCoordinates(const float* coords, pdbxWriter* writer);
static void writeFile(pdbxWriter* writer);

// class implementation 

static pdbxParser* create_pdbxParser(const char* filepath) {
  pdbxParser* parser = new pdbxParser;
  char buffer[BUFFER_SIZE];
  int numberAtoms;
  parser->xyz = NULL;
  parser->hashMem = NULL;
  parser->chain_auth = NULL;
  parser->resid_auth = NULL;
  parser->type_auth = NULL;
  parser->error = false;
  parser->bondsTo = NULL;
  parser->bondsFrom = NULL;
  parser->file = fopen(filepath, "r");
  if (!parser->file) {
    printf("pdbxplugin) cannot open file %s\n", filepath);
    return NULL;
  }
  if (NULL == fgets(buffer, BUFFER_SIZE, parser->file)) {
    printf("pdbxplugin) cannot read file %s\n", filepath);
    return NULL;
  }

  /* Find the number of atoms */
  parser->natoms = parseNumberAtoms(parser);
  numberAtoms = parser->natoms;
  if (parser->natoms <= 0) {
    printf("pdbxplugin) Could not get atom number\n");
    return NULL;
  }
  initCharToNum();
  parser->xyz = new float[numberAtoms*3];
  parser->hashMem = new list_node[numberAtoms+1];
  parser->chain_auth = new char[numberAtoms*CHAIN_SIZE];
  parser->resid_auth = new int [numberAtoms];
  parser->type_auth = new char[numberAtoms * TYPE_SIZE];
  return parser;
}

void delete_pdbxParser(pdbxParser* parser) {
  fclose(parser->file);
  if (parser->xyz != NULL) {
    delete [] parser->xyz;
    parser->xyz = NULL;
  }
  if (parser->type_auth != NULL) {
    delete [] parser->type_auth;
    parser->type_auth = NULL;
  }
  if (parser->resid_auth != NULL) {
    delete [] parser->resid_auth;
    parser->resid_auth = NULL;
  }
  if (parser->hashMem != NULL) {
    delete [] parser->hashMem;
    parser->hashMem = NULL;
  }
  if (parser->chain_auth != NULL) {
    delete [] parser->chain_auth;
    parser->chain_auth = NULL;
  }
  if (parser->type_auth != NULL) {
    inthash_destroy(&parser->bondHash);
  }
}

static void skipNextWord(char * str, void * word, int& pos) {
  /* Handle case if we start at end of line */
  if (str[pos] == '\0' || str[pos] == '\n') {
    return;
  }
  /* move forward until we hit non-whitespace */
  while(str[pos] == ' ') {
    ++pos;
  }
  /* increment pos until we hit a whitespace */
  while (str[pos++] != ' ') {}
}

static void getNextWord(char * str, void * word, int& pos) {
  char * w = (char*) word;
  int wordpos = 0;
  /* Handle case if we start at end of line */
  if (str[pos] == '\0' || str[pos] == '\n') {
    return;
  }
  /* move forward until we hit non-whitespace */
  while(str[pos] == ' ') {
    ++pos;
  }
  /* increment pos until we hit a whitespace */
  while (str[pos] != ' ') {
    w[wordpos++] = str[pos++];
  }
  w[wordpos] = '\0';
  /* Increment pos to point to first char that has not been read */
  ++pos;
}

static float stringToFloat(char * str) {
  bool neg = (str[0] == '-');
  unsigned int total = 0;
  unsigned pos = neg ? 1 : 0;
  unsigned int num = 0;
  unsigned int denom = 1;
  float retval;
  /* calculate integer before the decimal */
  while (str[pos] != '.') {
    total = (total*10) + str[pos] - '0';
    ++pos;
  }
  ++pos;
  /* Find the fraction representing the decimal */
  while (str[pos] != '\0') {
    num = (num * 10) + str[pos] - '0';
    denom *= 10;
    ++pos;
  }
  retval = (float)total + (double)num/(double)denom;
  if (neg)
    retval *= -1;
  return retval;
}

static int parseNumberAtoms(pdbxParser* parser) {
  char buffer[BUFFER_SIZE];
  char wordbuffer[64];
  int numatoms = 0;
  int i;
  int tableSize = 0;

  // skip past junk at start of file, stop when we get to atomSite data
  while (NULL == strstr(buffer, "_atom_site.")) {
    // if this is true then we couldnt find the numatoms
    if (NULL == fgets(buffer, BUFFER_SIZE, parser->file))
      return -1;
  }

  while (!(NULL == strstr(buffer, "_atom_site."))) {
    sscanf(buffer+11, "%s", wordbuffer); // table is used in parseStructure too
    /* assign integer values to each column */
    if (0 == strcmp(wordbuffer, "id")) {
      parser->table[tableSize] = COLUMN_NUMBER;
    } else if (0 == strcmp(wordbuffer, "type_symbol")) {
      parser->table[tableSize] = COLUMN_NAME;
    } else if (0 == strcmp(wordbuffer, "label_comp_id")) {
      parser->table[tableSize] = COLUMN_RESNAME;
    } else if (0 == strcmp(wordbuffer, "label_asym_id")) {
      parser->table[tableSize] = COLUMN_CHAIN;
    } else if (0 == strcmp(wordbuffer, "auth_asym_id")) {
      parser->table[tableSize] = COLUMN_CHAIN_AUTH;
    } else if (0 == strcmp(wordbuffer, "Cartn_x")) {
      parser->table[tableSize] = COLUMN_X;
    } else if (0 == strcmp(wordbuffer, "Cartn_y")) {
      parser->table[tableSize] = COLUMN_Y;
    } else if (0 == strcmp(wordbuffer, "Cartn_z")) {
      parser->table[tableSize] = COLUMN_Z;
    } else if (0 == strcmp(wordbuffer, "label_seq_id")) {
      parser->table[tableSize] = COLUMN_RESID;
    } else if (0 == strcmp(wordbuffer, "auth_seq_id")) {
      parser->table[tableSize] = COLUMN_RESID_AUTH;
    } else if (0 == strcmp(wordbuffer, "pdbx_PDB_ins_code")) {
      parser->table[tableSize] = COLUMN_INSERTION;
    } else if (0 == strcmp(wordbuffer, "B_iso_or_equiv")) {
      parser->table[tableSize] = COLUMN_BFACTOR;
    } else if (0 == strcmp(wordbuffer, "occupancy")) {
      parser->table[tableSize] = COLUMN_OCCUPANCY;
    } else if (0 == strcmp(wordbuffer, "label_atom_id")) {
      parser->table[tableSize] = COLUMN_TYPE;
    } else if (0 == strcmp(wordbuffer, "auth_atom_id")) {
      parser->table[tableSize] = COLUMN_TYPE_AUTH;
    } else if (0 == strcmp(wordbuffer, "pdbx_formal_charge")) {
      parser->table[tableSize] = COLUMN_CHARGE;
    } else {
      parser->table[tableSize] = COLUMN_JUNK;
    }

    // if this is true then we couldnt find the numatoms
    if (NULL == fgets(buffer, BUFFER_SIZE, parser->file))
      return -1;

    tableSize++;
  }

  // increment numatoms until we get to the end of the file
  while (buffer[0] != '#') {
    // if this is true then we couldnt find the numatoms
    if (NULL == fgets(buffer, BUFFER_SIZE, parser->file))
      return -1;
    ++numatoms;
  }

  rewind(parser->file);
  /* Cut off any junk columns from the end of table */

  i = tableSize;
  while(parser->table[--i] == COLUMN_JUNK){}
  tableSize = i+1;
  parser->tableSize = tableSize;

  if (numatoms == 0) {
    printf("pdbxplugin) Could not parse atom number from file\n");
    return -1;
  }

  return numatoms;
}


static void initCharToNum() {
  int i;
  int j = 1;

  for (i=0; i<128; i++)
    charToNum[i] = -1;

  i = 'A';
  while (i <= 'Z')
    charToNum[i++] = j++;
  i = 'a';
  while (i <= 'z')
    charToNum[i++] = j++;
  i = '0';
  while (i <= '9')
    charToNum[i++] = j++;
}


static inline int getUniqueResID(char * chainstr, int resid) {
  int i, uid;
  int length = strlen(chainstr);
  // Assuming max length of chainstr is 3 chars
  //Each char can be respresented by <= 6 bits since only a-z, A-Z, and 0-9 are valid values (62 possible values)
  uid = 1 + charToNum[(int)chainstr[0]];
  uid <<= 6;

  if (length == 1) {
    uid <<= 12;
  } else if (length == 2) {
    uid += charToNum[(int)chainstr[1]];
    uid <<= 12;
  } else if (length == 3) {
    uid += charToNum[(int)chainstr[1]];
    uid = (uid << 6) + charToNum[(int) chainstr[2]];
    uid <<= 6;
  }

  // First 18 bits of uid dedicated to 3 letters of chainstr
  uid <<= 12;
  uid += (0xFFF & resid); //add 12 least significant bits of resid to fill the remaining 10 bits of uid

  return uid;
}



#define ATOM_TYPE      0
#define ATOM_RESNAME   1
#define ATOM_INSERTION 2
#define ATOM_CHAIN     3
#define MAX_COLUMNS    64
#define MAX_OPTIONAL_AUTH_FIELDS  2

#define FLAG_CHAIN_LENGTH 0x01
#define FLAG_CHARGE       0x02
#define FLAG_INSERTION    0x04
#define FLAG_BFACTOR      0x08
#define FLAG_OCCUPANCY    0x10

static int parseStructureFaster(molfile_atom_t * atoms, int * optflags, pdbxParser* parser) {
  int i, count, atomdata, pos, idx, xyzcount, j;
  char buffer[BUFFER_SIZE];
  char namebuffer[8];
  char occupancybuffer[16];
  char bfactorbuffer[16];
  char chargebuffer[16];
  char residbuffer[8];
  char residAuthbuffer[8];
  char chainbuffer[16];
  char trash[16];
  char xbuffer[16];
  char ybuffer[16];
  char zbuffer[16];
  void * pillars[MAX_COLUMNS];
  molfile_atom_t * atom;
  int badptecount = 0;
  int chargecount = 0;
  int occupancycount = 0;
  int bfactorcount = 0;
  char oldChain[8];
  unsigned char parseFlags = 0;
  chainbuffer[1] = '\0';
  chainbuffer[2] = '\0';
  int hashTemp;
  int hashCount = 1;
  int head;
  int chainAuthIdx = MAX_COLUMNS-1, typeIdx = MAX_COLUMNS-1, resnameIdx = MAX_COLUMNS-1;
  int insertionIdx = MAX_COLUMNS-1, typeAuthIdx = MAX_COLUMNS-1;
#if (vmdplugin_ABIVERSION >= 20)
  int chainIdx = MAX_COLUMNS-1;
#endif
  char * chainAuth = parser->chain_auth;
  char * typeAuth = parser->type_auth;
  int tableSize = parser->tableSize;
  int* table = parser->table;
  unsigned char doBonds = 0;

  /* Initialize hash table used later when reading the special bonds */
  inthash_init(&parser->bondHash, parser->natoms);

  for (i=0; i<tableSize; i++) {
    switch (table[i]) {
      case COLUMN_NUMBER:
        pillars[i] = trash;
      break;
      case COLUMN_NAME:
        pillars[i] = namebuffer;
      break;
      case COLUMN_TYPE:
        pillars[i] = atoms->type;
        typeIdx = i;
      break;
      case COLUMN_TYPE_AUTH:
        pillars[i] = typeAuth;
        typeAuthIdx = i;
      break;
      case COLUMN_RESNAME:
        pillars[i] = atoms->resname;
        resnameIdx = i;
      break;
      case COLUMN_RESID:
        pillars[i] = residbuffer;
      break;
      case COLUMN_RESID_AUTH:
        pillars[i] = residAuthbuffer;
        doBonds++;
      break;
      case COLUMN_INSERTION:
        pillars[i] = atoms->insertion;
        insertionIdx = i;
      break;
      case COLUMN_X:
        pillars[i] = xbuffer;
      break;
      case COLUMN_Y:
        pillars[i] = ybuffer;
      break;
      case COLUMN_Z:
        pillars[i] = zbuffer;
      break;
      case COLUMN_OCCUPANCY:
        pillars[i] = occupancybuffer;
      break;
      case COLUMN_BFACTOR:
        pillars[i] = bfactorbuffer;
      break;
      case COLUMN_CHARGE:
        pillars[i] = chargebuffer;
      break;
      case COLUMN_CHAIN:
#if (vmdplugin_ABIVERSION < 20)
        pillars[i] = chainbuffer;
#else
        pillars[i] = atoms->chain;
        chainIdx = i;
#endif
      break;
      case COLUMN_CHAIN_AUTH:
        pillars[i] = chainAuth;
        chainAuthIdx = i;
        doBonds++;
      break;
      default:
        pillars[i] = trash;
      break;
    }
  }

  /* If the two optional auth fields are not present, don't look for extra bonds */
  if (doBonds != MAX_OPTIONAL_AUTH_FIELDS) {
    doBonds = 0;
  }

  atomdata = 0;
  /* Start parsing     */
  /* Skip through junk */
  while (!atomdata) {
    if (NULL == fgets(buffer, BUFFER_SIZE, parser->file)) {
      printf("pdbxplugin) failure while reading file\n");
      parser->error = true;
      return -1;
    }
    if (!(NULL == strstr(buffer, "_atom_site.")))
      atomdata = 1;
  }

  /* Skip through the atomdata table declaration*/
  while (atomdata) {
    if (NULL == fgets(buffer, BUFFER_SIZE, parser->file)) {
      printf("pdbxplugin) failure while reading file\n");
      parser->error = true;
      return -1;
    }
    if (NULL == strstr(buffer, "_atom_site."))
      atomdata = 0;
  }

  count = 0;
  atom = atoms;
  do {
    pos = 0;
    for (i=0; i<tableSize; ++i) {
      if (table[i] == COLUMN_JUNK) {
        /* if we don't want this column, update pos to point to the next column */
        skipNextWord(buffer, buffer, pos);
      }
      else {
        /* will copy each column string into the atom struct */
        /* or save the string if we need to convert it */
        getNextWord(buffer, pillars[i], pos);
      }
    }
    /* Coordinates must be saved until timestep is called */
    xyzcount = count*3;
    /*replacing atof with stringToFloat will increase performance */
    parser->xyz[xyzcount] = atof(xbuffer);
    parser->xyz[xyzcount+1] = atof(ybuffer);
    parser->xyz[xyzcount+2] = atof(zbuffer);

    atom->resid = atoi(residbuffer);
    if (doBonds && residAuthbuffer[0] != '.' && residAuthbuffer[0] != '?') {
      parser->resid_auth[count] = atoi(residAuthbuffer);

      /* add atom to hash table */
      hashTemp = getUniqueResID(chainAuth, parser->resid_auth[count]);

      if (-1 != (head = inthash_insert(&parser->bondHash, hashTemp, hashCount))) {
        /* key already exists, so we have to "add" a node to the linked list for this residue */
        /* Since we can't change the pointer in the hash table, we insert the node at the second
         * position in the list
         */
        parser->hashMem[hashCount].next = parser->hashMem[head].next;
        parser->hashMem[head].next = hashCount;
      }
      /* "add" node to list */
      parser->hashMem[hashCount++].index = count;
    }

    // XXX replace '?' or '.' insertion codes with a NUL char 
    // indicating an empty insertion code.
    if (atom->insertion[0] == '?' || atom->insertion[0] == '.') {
      atom->insertion[0] = '\0';
    }

//TODO: figure out what this conditional should be
#if (vmdplugin_ABIVERSION < 20)
    /* check to see if the chain length is greater than 2 */
    if (chainbuffer[2] != '\0' && chainbuffer[1] != '\0') {
      chainbuffer[2] = '\0';
      parseFlags |= FLAG_CHAIN_LENGTH;
    }
    atom->chain[0] = chainbuffer[0];
    atom->chain[1] = chainbuffer[1];
#endif

    /* Assign these to the pdbx_data struct */
    if (bfactorbuffer[0] != '.' && bfactorbuffer[0] != '.') {
      atom->bfactor = atof(bfactorbuffer);
      ++bfactorcount;
      parseFlags |= FLAG_BFACTOR;
    }
    else
      atom->bfactor = 0.0;

    if (occupancybuffer[0] != '.' && occupancybuffer[0] != '?') {
      atom->occupancy = atof(occupancybuffer);
      ++occupancycount;
      parseFlags |= FLAG_OCCUPANCY;
    }
    else 
      atom->occupancy = 0.0;

    if (chargebuffer[0] != '.' && chargebuffer[0] != '?') {
      atom->charge = atof(chargebuffer);
      ++chargecount;
      parseFlags |= FLAG_CHARGE;
    }
    else
      atom->charge = 0.0;

    idx = get_pte_idx_from_string(namebuffer);

    /* check for parenthesis in atom type */
    if (atom->type[0] == '"') {
      /* only save what is inside the parenthesis */
      i = 1;
      while (atom->type[i] != '"') {
        atom->type[i-1] = atom->type[i];
        ++i;
      }
      atom->type[i-1] = '\0';
    }
    /* atom->name and atom-> are the same */
    strcpy(atom->name, atom->type);

    /* Set periodic table values */
    if (idx) {
      atom->atomicnumber = idx;
      atom->mass = get_pte_mass(idx);
      atom->radius = get_pte_vdw_radius(idx);
    }
    else {
      ++badptecount;
    }
    if (NULL == fgets(buffer, BUFFER_SIZE, parser->file)) {
      printf("pdbxplugin) failure while reading file\n");
      parser->error = true;
      return -1;
    }

    ++count;
    ++atom;
    typeAuth += TYPE_SIZE;
    if (doBonds) {
      chainAuth += CHAIN_SIZE;
      pillars[chainAuthIdx] = chainAuth;
    }
    pillars[typeAuthIdx] = typeAuth;
    pillars[typeIdx] = atom->type;
    pillars[resnameIdx] = atom->resname;
    pillars[insertionIdx] = atom->insertion;
#if (vmdplugin_ABIVERSION >= 20)
    pillars[chainIdx] = atom->chain;
#endif
  } while (buffer[0] != '#'); //do until all the atoms have been read

  /* after we finish parsing, set optflags */
#if (vmdplugin_ABIVERSION < 20)
  if (parseFlags & FLAG_CHAIN_LENGTH) {
    printf("pdbxplugin) WARNING: This plugin ABI does not support chain names longer than two characters. Some chain names have been truncated.\n");
  }
#endif

  if (badptecount == 0)
    *optflags |= MOLFILE_MASS | MOLFILE_RADIUS | MOLFILE_ATOMICNUMBER;

  if (parseFlags & FLAG_CHARGE)
    *optflags |= MOLFILE_CHARGE;

  if (parseFlags & FLAG_BFACTOR)
    *optflags |= MOLFILE_BFACTOR;

  if (parseFlags & FLAG_OCCUPANCY)
    *optflags |= MOLFILE_OCCUPANCY;


  if (badptecount > 0) {
    printf("pdbxplugin) encountered %d bad element indices!\n", badptecount);
    return -1;
  }

  return 0;
}

static int setCoordinatesFast(float* coords, pdbxParser* parser) {
  int i, j;
  j = 0;
  for(i=0; i<parser->natoms; i++) {
    coords[0] = parser->xyz[j];
    coords[1] = parser->xyz[j+1];
    coords[2] = parser->xyz[j+2];
    coords += 3;
    j += 3;
  }
  return 0;
}

#define BOND_JUNK      0
#define BOND_NAME_1    1
#define BOND_CHAIN_1   2
#define BOND_RESNAME_1 3
#define BOND_RESID_1   4
#define BOND_NAME_2    5
#define BOND_CHAIN_2   6
#define BOND_RESNAME_2 7
#define BOND_RESID_2   8

static bool readAngleBonds(molfile_atom_t * atoms, pdbxParser* parser) {
  char buffer[BUFFER_SIZE];
  int bondTable[32];
  char* columns[32];
  int bondTableSize = 0;
  int bnum = 0;
  int i, pos, j, k;
  int* newBondsTo, *newBondsFrom;
  fpos_t filePos;
  char junk[16];
  char modelNum[8];
  char name1[8];
  char wordbuffer[16];
  char name2[8];
  char chain1[8];
  char chain2[8];
  char resid1buffer[8];
  char resid2buffer[8];
  int resid1, resid2;
  int uid1, uid2;
  int aIdx1, aIdx2;
  molfile_atom_t * atom;

  /* skip through the file until we find the bond information */
  do {
    if (NULL == fgets(buffer, BUFFER_SIZE, parser->file)) {
      return false;
    }
  } while (NULL == strstr(buffer,"_pdbx_validate_rmsd_angle."));

  fgetpos(parser->file, &filePos);

  //if (sscanf(  return if two words in one table definition line

  /* Parse table header data */
  while (NULL != strstr(buffer,"_pdbx_validate_rmsd_angle.")) {
    sscanf(buffer+26, "%s", wordbuffer); // table is used in parseStructure too
    /* assign integer values to each column */
    if (0 == strcmp(wordbuffer, "auth_atom_id_1")) {
      columns[bondTableSize] = (char*)name1;
    }
    else if (0 == strcmp(wordbuffer, "auth_asym_id_1")) {
      columns[bondTableSize] = (char*)chain1;
    }
    else if (0 == strcmp(wordbuffer, "auth_comp_id_1")) {
      columns[bondTableSize] = (char*)junk;
    }
    else if (0 == strcmp(wordbuffer, "auth_seq_id_1")) {
      columns[bondTableSize] = (char*)resid1buffer;
    }
    else if (0 == strcmp(wordbuffer, "auth_atom_id_2")) {
      columns[bondTableSize] = (char*)name2;
    }
    else if (0 == strcmp(wordbuffer, "auth_asym_id_2")) {
      columns[bondTableSize] = (char*)chain2;
    }
    else if (0 == strcmp(wordbuffer, "auth_comp_id_2")) {
      columns[bondTableSize] = (char*)junk;
    }
    else if (0 == strcmp(wordbuffer, "auth_seq_id_2")) {
      columns[bondTableSize] = (char*)resid2buffer;
    }
    else {
      columns[bondTableSize] = (char*)junk;
    }
    ++bondTableSize;

    if (NULL == fgets(buffer, BUFFER_SIZE, parser->file)) {
      printf("pdbxplugin) could not read bond information.\n");
      return false;
    }

  }

  /* figure out how many bonds are being defined */
  while (buffer[0] != '#') {
    if (NULL == fgets(buffer, BUFFER_SIZE, parser->file)) {
      printf("pdbxplugin) could not read bond information.\n");
      return false;
    }
    ++bnum;
  }

  int test = parser->nbonds + bnum;
  if ((newBondsTo = (int*)realloc((void*)parser->bondsTo, (parser->nbonds + bnum) * sizeof(int))) == NULL)
    return false;
  if ((newBondsFrom = (int*)realloc((void*)parser->bondsFrom, (parser->nbonds + bnum) * sizeof(int))) == NULL)
    return false;
  parser->bondsTo = newBondsTo;
  parser->bondsFrom = newBondsFrom;

  /* Skip back to the start of the bond info */
  fsetpos(parser->file, &filePos);
  if (NULL == fgets(buffer, BUFFER_SIZE, parser->file)) {
    printf("pdbxplugin) could not read bond information.\n");
    return false;
  }
  /* Skip through the header */
  while (NULL != strstr(buffer,"_pdbx_validate_rmsd_angle.")) {
    if (NULL == fgets(buffer, BUFFER_SIZE, parser->file)) {
      printf("pdbxplugin) could not read bond information.\n");
      return false;
    }
  }

  bnum = 0;
  while (buffer[0] != '#') {
    pos = 0;
    /* copy each column of the table into the appropriate columns index */
    for (i=0; i<bondTableSize; ++i) {
        getNextWord(buffer, columns[i], pos);
    }
    resid1 = atoi(resid1buffer);
    resid2 = atoi(resid2buffer);
    /* get unique res ID for hash table lookup */
    uid1 = getUniqueResID(chain1, resid1);
    uid2 = getUniqueResID(chain2, resid2);
    k = 0;

    /* find the atoms in the hash table */
    if ( ((uid1 = inthash_lookup(&parser->bondHash, uid1)) != -1) && ((uid2 = inthash_lookup(&parser->bondHash, uid2)) != -1) ) {
      // because the hashtable is residue specifc, loop through 
      // all atoms in the residue to find the correct one
      // Find atom 1 
      do {
        aIdx1 = parser->hashMem[uid1].index;
        if (strcmp(name1, parser->type_auth + aIdx1 * TYPE_SIZE) == 0 && 
            parser->resid_auth[aIdx1] == resid1 &&
            strcmp(chain1, parser->chain_auth + aIdx1 * CHAIN_SIZE) == 0) {
          k++;
          break;
        } else {
          uid1 = parser->hashMem[uid1].next;
        }
      } while (uid1 != 0); //0 indicates end of "list"

      // Find atom 2
      do {
        aIdx2 = parser->hashMem[uid2].index;
        if (strcmp(name2, parser->type_auth + aIdx2 * TYPE_SIZE) == 0 && 
            parser->resid_auth[aIdx2] == resid2 &&
            strcmp(chain2, parser->chain_auth + aIdx2 * CHAIN_SIZE) == 0) {
          k++;
          break;
        } else {
          uid2 = parser->hashMem[uid2].next;
        }
      } while (uid2 != 0); // 0 indicates end of "list"

      if (k == 2) {
        parser->bondsFrom[parser->nbonds + bnum] = aIdx1+1;  // vmd doesn't use 0 based index for bond info?
        parser->bondsTo[parser->nbonds + bnum] = aIdx2+1;
        ++bnum;
      }
    }
#ifdef PDBX_DEBUG
    else {
      printf("^^^^Could locate bond^^^^^, %s %d\n", chain1, resid1);
      printf("Error finding atom ");
      if (uid1 == 0)
        printf("1 ");
      if (uid2 == 0)
        printf("2");
      printf("\n");
    }
#endif

    if (NULL == fgets(buffer, BUFFER_SIZE, parser->file)) {
      printf("pdbxplugin) could not read RMSD bond deviation information.\n");
      return false;
    }
  }

  parser->nbonds += bnum;
#ifdef PDBX_DEBUG
  printf("pdbxplugin) nbonds defined: %d\n", nbonds);
#endif
  return (bnum == 0) ? false : true;
}
 

static bool readRMSDBonds(molfile_atom_t * atoms, pdbxParser* parser) {
  char buffer[BUFFER_SIZE];
  int bondTable[32];
  char* columns[32];
  int bondTableSize = 0;
  int bnum = 0;
  int i, pos, j, k;
  fpos_t filePos;
  char junk[16];
  char modelNum[8];
  char name1[8];
  char wordbuffer[16];
  char name2[8];
  char chain1[8];
  char chain2[8];
  char resid1buffer[8];
  char resid2buffer[8];
  int resid1, resid2;
  int uid1, uid2;
  int aIdx1, aIdx2;
  molfile_atom_t * atom;

  /* skip through the file until we find the bond information */
  do {
    if (NULL == fgets(buffer, BUFFER_SIZE, parser->file)) {
      parser->nbonds = 0;
      return false;
    }
  } while (!isValidateRMSDBond(buffer));

  fgetpos(parser->file, &filePos);

  //if (sscanf(  return if two words in one table definition line

  /* Parse table header data */
  while (isValidateRMSDBond(buffer)) {
    sscanf(buffer+25, "%s", wordbuffer); // table is used in parseStructure too
    /* assign integer values to each column */
    if (0 == strcmp(wordbuffer, "auth_atom_id_1")) {
      columns[bondTableSize] = (char*)name1;
    }
    else if (0 == strcmp(wordbuffer, "auth_asym_id_1")) {
      columns[bondTableSize] = (char*)chain1;
    }
    else if (0 == strcmp(wordbuffer, "auth_comp_id_1")) {
      columns[bondTableSize] = (char*)junk;
    }
    else if (0 == strcmp(wordbuffer, "auth_seq_id_1")) {
      columns[bondTableSize] = (char*)resid1buffer;
    }
    else if (0 == strcmp(wordbuffer, "auth_atom_id_2")) {
      columns[bondTableSize] = (char*)name2;
    }
    else if (0 == strcmp(wordbuffer, "auth_asym_id_2")) {
      columns[bondTableSize] = (char*)chain2;
    }
    else if (0 == strcmp(wordbuffer, "auth_comp_id_2")) {
      columns[bondTableSize] = (char*)junk;
    }
    else if (0 == strcmp(wordbuffer, "auth_seq_id_2")) {
      columns[bondTableSize] = (char*)resid2buffer;
    }
    else {
      columns[bondTableSize] = (char*)junk;
    }
    ++bondTableSize;

    if (NULL == fgets(buffer, BUFFER_SIZE, parser->file)) {
      printf("pdbxplugin) could not read bond information.\n");
      return false;
    }
    
  }

  /* figure out how many bonds are being defined */
  while (buffer[0] != '#') {
    if (NULL == fgets(buffer, BUFFER_SIZE, parser->file)) {
      printf("pdbxplugin) could not read bond information.\n");
      return false;
    }
    ++bnum;
  }

  parser->nbonds = bnum;
  parser->bondsTo = (int*)malloc(bnum * sizeof(int));
  parser->bondsFrom = (int*)malloc(bnum * sizeof(int));

  /* Skip back to the start of the bond info */
  fsetpos(parser->file, &filePos);
  if (NULL == fgets(buffer, BUFFER_SIZE, parser->file)) {
    printf("pdbxplugin) could not read bond information.\n");
    return false;
  }
  /* Skip through the header */
  while (isValidateRMSDBond(buffer)) {
    if (NULL == fgets(buffer, BUFFER_SIZE, parser->file)) {
      printf("pdbxplugin) could not read bond information.\n");
      return false;
    }
  }

  bnum = 0;
  while (buffer[0] != '#') {
    pos = 0;
    /* copy each column of the table into the appropriate columns index */
    for (i=0; i<bondTableSize; ++i) {
        getNextWord(buffer, columns[i], pos);
    }
    resid1 = atoi(resid1buffer);
    resid2 = atoi(resid2buffer);
    /* get unique res ID for hash table lookup */
    uid1 = getUniqueResID(chain1, resid1);
    uid2 = getUniqueResID(chain2, resid2);
    k = 0;

    /* find the atoms in the hash table */
    if ( ((uid1 = inthash_lookup(&parser->bondHash, uid1)) != -1) && ((uid2 = inthash_lookup(&parser->bondHash, uid2)) != -1) ) {
      /* because the hashtable is residue specifc, loop through all atoms in the residue to find the correct one */
      /* Find atom 1 */
      do {
        aIdx1 = parser->hashMem[uid1].index;
        if (strcmp(name1, parser->type_auth + aIdx1 * TYPE_SIZE) == 0 && parser->resid_auth[aIdx1] == resid1 &&
            strcmp(chain1, parser->chain_auth + aIdx1 * CHAIN_SIZE) == 0){
          k++;
          break;
        }
        else
          uid1 = parser->hashMem[uid1].next;
      } while (uid1 != 0); //0 indicates end of "list"

      /* find atom 2 */
      do {
        aIdx2 = parser->hashMem[uid2].index;
        if (strcmp(name2, parser->type_auth + aIdx2 * TYPE_SIZE) == 0 && parser->resid_auth[aIdx2] == resid2 &&
            strcmp(chain2, parser->chain_auth + aIdx2 * CHAIN_SIZE) == 0){
          k++;
          break;
        }
        else
          uid2 = parser->hashMem[uid2].next;
      } while (uid2 != 0); // 0 indicates end of "list"

      if (k == 2) {
        parser->bondsFrom[bnum] = aIdx1+1;  // vmd doesn't use 0 based index for bond info?
        parser->bondsTo[bnum] = aIdx2+1;
        ++bnum;
      }
    }
#ifdef PDBX_DEBUG
    else {
      printf("^^^^Could locate bond^^^^^, %s %d\n", chain1, resid1);
      printf("Error finding atom ");
      if (uid1 == 0)
        printf("1 ");
      if (uid2 == 0)
        printf("2");
      printf("\n");
    }
#endif

    if (NULL == fgets(buffer, BUFFER_SIZE, parser->file)) {
      printf("pdbxplugin) could not read RMSD bond deviation information.\n");
      return false;
    }
  }
  parser->nbonds = bnum;
#ifdef PDBX_DEBUG
  printf("pdbxplugin) nbonds defined: %d\n", nbonds);
#endif
  return (bnum > 0);
}


static bool readBonds(molfile_atom_t * atoms, pdbxParser* parser) {
  bool retval = false;
  retval |= readRMSDBonds(atoms, parser);
  retval |= readAngleBonds(atoms, parser);
  return retval;
}


static inline bool isValidateRMSDBond(char * str) {
  /* return str[0-24] == "_pdbx_validate_rmsd_bond." */
  return (str[0] == '_' && str[1] == 'p' && str[2] == 'd' && str[3] == 'b' &&
          str[4] == 'x' && str[5] == '_' && str[6] == 'v' && str[7] == 'a' &&
          str[8] == 'l' && str[9] == 'i' && str[10] == 'd' && str[11] == 'a' &&
          str[12] == 't' && str[13] == 'e' && str[14] == '_' && str[15] == 'r' &&
          str[16] == 'm' && str[17] == 's' && str[18] == 'd' && str[19] == '_' &&
          str[20] == 'b' && str[21] == 'o' && str[22] == 'n' && str[23] == 'd' &&
          str[24] == '.');
}

static inline bool isAtomSite(char * str) {
  return (str[0]=='_' && str[1]=='a' && str[2]=='t' && str[3]=='o' && str[4]=='m' &&
          str[5]=='_' && str[6]=='s' && str[7]=='i' && str[8]=='t' && str[9]=='e' &&
          str[10]=='.');
}

/* start of pdbxWriter implementation */

static pdbxWriter* create_pdbxWriter(const char* filename, int numAtoms) {
  pdbxWriter* writer = new pdbxWriter;
  int length = strlen(filename);
  int start = 0;
  int end = length;
  int i;
  writer->numatoms = numAtoms;
  writer->bufferCount = 0;

  writer->fd = fopen(filename, "w");
  /* get name of pdb file */
  for (i=0; i<length; ++i) {
    if (filename[i] == '/' || filename[i] == '\\') {
      if (i+1 < length)
        start = i+1;
    }
    if (filename[i] == '.')
      end = i;
  }
  strncpy(writer->pdbName, filename + start, end - start);
  writer->pdbName[end-start] = '\0';
  return writer;
}

static void addCoordinates(const float* coords, pdbxWriter* writer) {
  writer->coordinates = coords;
}

static void addAtoms(const molfile_atom_t* atomlist, int optflags, pdbxWriter* writer) {
  int i;
  writer->atoms = new molfile_atom_t[writer->numatoms];
  molfile_atom_t* atoms = writer->atoms;
  
  memcpy(atoms, atomlist, writer->numatoms * sizeof(molfile_atom_t));

  /* If occ, bfactor, and insertion aren't given, we assign defaultvalues. */
  if (!(optflags & MOLFILE_OCCUPANCY)) {
    for (i=0; i<writer->numatoms; i++)
      atoms[i].occupancy = 0.0f;
  }
  if (!(optflags & MOLFILE_BFACTOR)) {
    for (i=0; i<writer->numatoms; i++)
      atoms[i].bfactor= 0.0f;
  }
  if (!(optflags & MOLFILE_INSERTION)) {
    for (i=0; i<writer->numatoms; i++) {
      atoms[i].insertion[0] =' ';
      atoms[i].insertion[1] ='\0';
    }
  }
  if (!(optflags & MOLFILE_ALTLOC)) {
    for (i=0; i<writer->numatoms; i++) {
      atoms[i].altloc[0]=' ';
      atoms[i].altloc[1]='\0';
    }
  }
  if (!(optflags & MOLFILE_ATOMICNUMBER)) {
    for (i=0; i<writer->numatoms; i++)
      atoms[i].atomicnumber = 0;
  }
}

static void writeAtomSite(pdbxWriter* writer) {
  char lineBuffer[BUFFER_SIZE];
  int i;
  const float* x, *y, *z;
  molfile_atom_t* atoms = writer->atoms;
  x = writer->coordinates;
  y = x+1;
  z = x+2;

  for (i=0; i<writer->numatoms; ++i) {
    sprintf(lineBuffer,"ATOM %d %s %s . %s %s . %d ? %f %f %f %f %f %s\n",
            i+1, atoms[i].name, atoms[i].type, atoms[i].resname, atoms[i].chain,
            atoms[i].resid, *x, *y, *z, atoms[i].occupancy,
            atoms[i].charge, atoms[i].chain);
    x += 3;
    y += 3;
    z += 3;
    write(lineBuffer, writer);
  }
}


static void writeFile(pdbxWriter* writer) {
  /* write PDBx header */
  writeIntro(writer);
  write(atomSiteHeader, writer);
  writeAtomSite(writer);
  write("#\n", writer);
  close(writer);
}

static void writeIntro(pdbxWriter* writer) {
  write("data_", writer);
  write(writer->pdbName, writer);
  write("\n", writer);
}

static void close(pdbxWriter* writer) {
  writeBuffer(writer);
  fclose(writer->fd);
}

static void write(const char* str, pdbxWriter* writer) {
  int length = strlen(str);
  int copy_size;
  int num_copied = 0;

  if (length + writer->bufferCount < WB_SIZE) {
    memcpy(writer->writeBuf + writer->bufferCount, str, length);
    writer->bufferCount += length;
  }
  else do {
    copy_size = WB_SIZE - writer->bufferCount;
    if (copy_size + num_copied > length) {
      copy_size = length - num_copied;
    }
    memcpy(writer->writeBuf + writer->bufferCount, str + num_copied, copy_size);
    writer->bufferCount += copy_size;
    num_copied += copy_size;
    if (writer->bufferCount == WB_SIZE) {
      writeBuffer(writer);
    }
  } while (num_copied < length);
}

static void writeBuffer(pdbxWriter* writer) {
  if (writer->bufferCount == 0)
    return;
  fwrite(writer->writeBuf, sizeof(char), writer->bufferCount, writer->fd);
  writer->bufferCount = 0;
}


/*
 * API functions start here
 */

typedef struct {
  pdbxParser * parser;
  pdbxWriter * writer;
  int natoms;
  molfile_atom_t *atomlist;
  molfile_metadata_t *meta;
  int readTS;
} pdbx_data;


static void * open_pdbx_read(const char *filepath, const char *filetype,
                             int *natoms) {
  pdbx_data *data;
  data = new pdbx_data;
  data->readTS = 0;
  data->parser = create_pdbxParser(filepath);
  data->natoms = data->parser->natoms;
  *natoms = data->natoms;
  if (*natoms == 0) //If no atoms were found this is not a pdb file
    return NULL;
  if (data->parser->error)
    return NULL;
  return data;
}

static int read_pdbx_structure(void * mydata, int *optflags, molfile_atom_t *atoms) {
  pdbx_data * data = (pdbx_data *)mydata;
  *optflags = MOLFILE_NOOPTIONS;
#if 0
  // XXX the current pdbx code doesn't actually throw any exceptions so the try
  //     block causes linkage problems on Android if it is used as-is.
  try {
    parseStructureFaster(atoms, optflags, data->parser);
  } catch (...) {
    printf("pdbxplugin) Error while trying to read file\n");
    return MOLFILE_ERROR;
  }
#else
  if (parseStructureFaster(atoms, optflags, data->parser)) {
    printf("pdbxplugin) Error while trying to parse pdbx structure\n");
    return MOLFILE_ERROR;
  }
#endif

  printf("pdbxplugin) Starting to read bonds...\n");
  readBonds(atoms, data->parser);
  *optflags |= MOLFILE_BONDSSPECIAL;
  return MOLFILE_SUCCESS;
}

static int read_bonds(void *v, int *nbonds, int **fromptr, int **toptr,
                      float ** bondorder,int **bondtype,
                      int *nbondtypes, char ***bondtypename) {
  pdbx_data * data = (pdbx_data *)v;
  if (data->parser->nbonds == 0) {
    *nbonds = 0;
    *fromptr = NULL;
    *toptr = NULL;
  } else {
    *nbonds = data->parser->nbonds;
    *fromptr = data->parser->bondsFrom;
    *toptr = data->parser->bondsTo;
  }
  *bondorder = NULL;
  *bondtype = NULL;
  *nbondtypes = 0;
  *bondtypename = NULL;

  return MOLFILE_SUCCESS;
}

static int read_pdbx_timestep(void * mydata, int natoms, molfile_timestep_t *ts) {
  pdbx_data * data = (pdbx_data *)mydata;
  if (data->readTS)
    return MOLFILE_ERROR;
  data->readTS = 1;
  setCoordinatesFast(ts->coords, data->parser);

  return MOLFILE_SUCCESS;
}

static void close_pdbx_read(void *v) {
  pdbx_data * data = (pdbx_data *)v;
  delete_pdbxParser(data->parser);
  delete data;
}

static void* open_file_write(const char *path, const char *filetypye, int natoms) {
  FILE* fd;
  pdbx_data * data = new pdbx_data;
  data->writer = create_pdbxWriter(path, natoms);
  return data;
}

static int write_structure(void* v, int optflags, const molfile_atom_t *atoms) {
  pdbx_data * data = (pdbx_data *)v;
  addAtoms(atoms, optflags, data->writer);
  return MOLFILE_SUCCESS;
}

static int write_timestep(void *v, const molfile_timestep_t* ts) {
  pdbx_data * data = (pdbx_data *)v;
  addCoordinates(ts->coords, data->writer);
  writeFile(data->writer);
  return MOLFILE_SUCCESS;
}

static void close_file_write(void* v) {
  pdbx_data* data = (pdbx_data*)v;
  delete data->writer;
  delete data;
}

/*
 * Initialization stuff down here
 */

static molfile_plugin_t plugin;

VMDPLUGIN_API int VMDPLUGIN_init() {
  memset(&plugin, 0, sizeof(molfile_plugin_t));
  plugin.abiversion = vmdplugin_ABIVERSION;
  plugin.type = MOLFILE_PLUGIN_TYPE;
  plugin.name = "pdbx";
  plugin.prettyname = "mmCIF/PDBX";
  plugin.author = "Brendan McMorrow";
  plugin.majorv = 0;
  plugin.minorv = 9;
  plugin.is_reentrant = VMDPLUGIN_THREADSAFE;
  plugin.filename_extension = "cif";
  plugin.open_file_read = open_pdbx_read;
  plugin.read_structure = read_pdbx_structure;
  plugin.read_next_timestep = read_pdbx_timestep;
  plugin.read_bonds = read_bonds;
  plugin.open_file_write = open_file_write;
  plugin.write_structure = write_structure;
  plugin.write_timestep = write_timestep;
  plugin.close_file_write = close_file_write;
  plugin.close_file_read = close_pdbx_read;
  return VMDPLUGIN_SUCCESS;
}

VMDPLUGIN_API int VMDPLUGIN_register(void *v, vmdplugin_register_cb cb) {
  (*cb)(v, (vmdplugin_t *)&plugin);
  return VMDPLUGIN_SUCCESS;
}

VMDPLUGIN_API int VMDPLUGIN_fini() {
  return VMDPLUGIN_SUCCESS;
}



#ifdef TEST_PLUGIN

int main(int argc, char *argv[]) {
  molfile_timestep_t timestep;
  pdbx_data *v;
  int natoms, bnum, nbtypes;
  float *border[1];
  float *x, *y, *z;
  int *btype[1];
  char **btypenames[1];
  int i, set;

//  while (--argc) {
    ++argv;

  struct timeval  tot1, tot2;
  gettimeofday(&tot1, NULL);
    if (*argv != NULL)
      v = (pdbx_data*)open_pdbx_read(*argv, "pdbx", &natoms);
    else
   //   v = (pdbx_data*)open_pdbx_read("/Users/Brendan/pdbx/3j3q.cif", "pdbx", &natoms);
      v = (pdbx_data*)open_pdbx_read("/home/brendanbc1/Downloads/3j3q.cif", "pdbx", &natoms);
    if (!v) {
      fprintf(stderr, "main) open_pdbx_read failed for file %s\n", *argv);
      return 1;
    }
    fprintf(stderr, "main) open_pdbx_read succeeded for file %s\n", *argv);
    fprintf(stderr, "main) number of atoms: %d\n", natoms);

    set = 0;
    molfile_atom_t * atoms = new molfile_atom_t[natoms];
    if (MOLFILE_SUCCESS == read_pdbx_structure(v, &set, atoms)){
      printf("xyz structure successfully read.\n");
    } else {
      fprintf(stderr, "main) error reading pdbx file\n");
    }

    i = 0;
    timestep.coords = (float *)malloc(3*sizeof(float)*natoms);
    if (!read_pdbx_timestep(v, natoms, &timestep)) {
      fprintf(stderr, "main) open_pdbx_read succeeded for file %s\n", *argv);
    } else {
      fprintf(stderr, "main) Failed to read timestep\n");
    }

    gettimeofday(&tot2, NULL);
    printf ("Total time to read file: %f seconds\n",
           (double) (tot2.tv_usec - tot1.tv_usec) / 1000000 +
           (double) (tot2.tv_sec - tot1.tv_sec));
    close_pdbx_read(v);


    printf("Writing file...\n");
    v = (pdbx_data*)open_file_write("/home/brendanbc1/test.cif", 0, natoms);

    //v = (pdbx_data*)open_file_write("/Users/Brendan/test.cif", 0, natoms);
    printf("File opened for writing...\n");
    //printf("%d\n", v->writer->numatoms);
    write_structure(v, set, (const molfile_atom_t*)atoms);
    printf("Structure information gathered...\n");
    write_timestep(v, &timestep);
    printf("File written...\n");
    close_file_write(v);
    printf("File closed.\n");
    delete [] atoms;
         /*
    printf("Writing pdbx.txt\n");
    x = timestep.coords; y = x+1;
    z = x+2;
    FILE *f;
    f = fopen("pdbx.txt", "w");
    for(i=0; i<natoms; i++) {
      fprintf(f, "%i %d %s %s %s %f %f %f\n", atoms[i].atomicnumber, atoms[i].resid, atoms[i].chain, atoms[i].resname, atoms[i].type, *x,*y,*z);
      //fprintf(stderr, "%i\t%s  %s\t%s  %s  %i  %f\t%f\t%f\t%f\t%f\t%f\n", i+1, atoms[i].name, atoms[i].type,
        //      atoms[i].chain, atoms[i].resname, atoms[i].resid, *x, *y, *z, atoms[i].occupancy, atoms[i].bfactor, atoms[i].charge);
      x+=3;
      y+=3;
      z+=3;
    }
    fclose(f);
    printf("main) pdbx.txt written\n");
 // }
 */
  return 0;
}
#endif
