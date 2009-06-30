/* MACHINE GENERATED FILE, DO NOT EDIT! */

#define VMDPLUGIN molfile_pdbplugin
#define STATIC_PLUGIN 1

/***************************************************************************
 *cr
 *cr            (C) Copyright 1995-2009 The Board of Trustees of the
 *cr                        University of Illinois
 *cr                         All Rights Reserved
 *cr
 ***************************************************************************/

/***************************************************************************
 * RCS INFORMATION:
 *
 *      $RCSfile: pdbplugin.c,v $
 *      $Author: johns $       $Locker:  $             $State: Exp $
 *      $Revision: 1.72 $       $Date: 2009/04/29 15:45:32 $
 *
 ***************************************************************************/

/*
 * PDB file format specifications:
 *   http:

 */

#include "largefiles.h"   /* platform dependent 64-bit file I/O defines */

#include "molfile_plugin.h"
#include "readpdb.h"
#include "periodic_table.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/*
 * API functions start here
 */

typedef struct {
  FILE *fd;
  int first_frame;
  int natoms;
  molfile_atom_t *atomlist;
  molfile_metadata_t *meta;
  int nconect;
  int nbonds, maxbnum;
  int *from, *to, *idxmap;
} pdbdata;

static void *open_pdb_read(const char *filepath, const char *filetype, 
    int *natoms) {
  FILE *fd;
  pdbdata *pdb;
  char pdbstr[PDB_BUFFER_LENGTH];
  int indx, nconect;

  fd = fopen(filepath, "r");
  if (!fd) 
    return NULL;
  pdb = (pdbdata *)malloc(sizeof(pdbdata));
  pdb->fd = fd;
  pdb->meta = (molfile_metadata_t *) malloc(sizeof(molfile_metadata_t));
  memset(pdb->meta, 0, sizeof(molfile_metadata_t));

  pdb->meta->remarklen = 0;
  pdb->meta->remarks = NULL;

  *natoms=0;
  nconect=0;
  do {
    indx = read_pdb_record(pdb->fd, pdbstr);
    if (indx == PDB_ATOM) {
      *natoms += 1;
    } else if (indx == PDB_CONECT) {
      nconect++;
    } else if (indx == PDB_HEADER) {
      get_pdb_header(pdbstr, pdb->meta->accession, pdb->meta->date, NULL);
      if (strlen(pdb->meta->accession) > 0) 
        strcpy(pdb->meta->database, "PDB");
    } else if (indx == PDB_REMARK || indx == PDB_CONECT || indx == PDB_UNKNOWN) {
      int len=strlen(pdbstr);
      int newlen = len + pdb->meta->remarklen;

      char *newstr=realloc(pdb->meta->remarks, newlen + 1);
      if (newstr != NULL) {
        pdb->meta->remarks = newstr;
        pdb->meta->remarks[pdb->meta->remarklen] = '\0';
        memcpy(pdb->meta->remarks + pdb->meta->remarklen, pdbstr, len);
        pdb->meta->remarks[newlen] = '\0';
        pdb->meta->remarklen = newlen;
      }
    }
 
  } while (indx != PDB_END && indx != PDB_EOF);

  /* If no atoms were found, this is probably not a PDB file! */
  if (!*natoms) {
    fprintf(stderr, "PDB file '%s' contains no atoms.\n", filepath);
    if (pdb->meta->remarks != NULL)
      free(pdb->meta->remarks);
    if (pdb->meta != NULL)
      free(pdb->meta);
    free(pdb);
    return NULL;
  }

  rewind(pdb->fd); /* if ok, rewind file and prepare to parse it for real */
  pdb->natoms = *natoms;
  pdb->nconect = nconect;
  pdb->nbonds = 0;
  pdb->maxbnum = 0;
  pdb->from = NULL;
  pdb->to = NULL;
  pdb->idxmap = NULL;
  pdb->atomlist = NULL;

#if defined(VMDUSECONECTRECORDS)
  /* allocate atom index translation table if we have 99,999 atoms or less */
  /* and we have conect records to process                                 */
  if (pdb->natoms < 100000 && pdb->nconect > 0) {
    pdb->idxmap = (int *) malloc(100000 * sizeof(int));
    memset(pdb->idxmap, 0, 100000 * sizeof(int));
  }
#endif
 
  return pdb; 
}

static int read_pdb_structure(void *mydata, int *optflags, 
    molfile_atom_t *atoms) { 
  pdbdata *pdb = (pdbdata *)mydata;
  molfile_atom_t *atom;
  char pdbrec[PDB_BUFFER_LENGTH];
  int i, rectype, atomserial, pteidx;
  char ridstr[8];
  char elementsymbol[3];
  int badptecount = 0;
  long fpos = ftell(pdb->fd);

  *optflags = MOLFILE_INSERTION | MOLFILE_OCCUPANCY | MOLFILE_BFACTOR |
              MOLFILE_ALTLOC | MOLFILE_ATOMICNUMBER | MOLFILE_BONDSSPECIAL;

  i = 0;
  do {
    rectype = read_pdb_record(pdb->fd, pdbrec);
    switch (rectype) {
    case PDB_ATOM:
      atom = atoms+i;
      get_pdb_fields(pdbrec, strlen(pdbrec), &atomserial, 
          atom->name, atom->resname, atom->chain, atom->segid, 
          ridstr, atom->insertion, atom->altloc, elementsymbol,
          NULL, NULL, NULL, &atom->occupancy, &atom->bfactor);

      if (pdb->idxmap != NULL && atomserial < 100000) {
        pdb->idxmap[atomserial] = i; /* record new serial number translation */ 
      }
 
      atom->resid = atoi(ridstr);

      /* determine atomic number from the element symbol */
      pteidx = get_pte_idx_from_string(elementsymbol);
      atom->atomicnumber = pteidx;
      if (pteidx != 0) {
        atom->mass = get_pte_mass(pteidx);
        atom->radius = get_pte_vdw_radius(pteidx);
      } else {
        badptecount++; /* unrecognized element */
      }
 
      strcpy(atom->type, atom->name);
      i++;
      break;

    case PDB_CONECT:
      /* only read CONECT records for structures where we know they can */
      /* be valid for all of the atoms in the structure                 */
      if (pdb->idxmap != NULL) {
        get_pdb_conect(pdbrec, pdb->natoms, pdb->idxmap, 
                       &pdb->maxbnum, &pdb->nbonds, &pdb->from, &pdb->to);
      }
      break;

    default:
      /* other record types are ignored in the structure callback */
      /* and are dealt with in the timestep callback or elsewhere */
      break;
    }
  } while (rectype != PDB_END && rectype != PDB_EOF);

  fseek(pdb->fd, fpos, SEEK_SET);

  /* if all atoms are recognized, set the mass and radius flags too,  */
  /* otherwise let VMD guess these for itself using it's own methods  */
  if (badptecount == 0) {
    *optflags |= MOLFILE_MASS | MOLFILE_RADIUS;
  }

  return MOLFILE_SUCCESS;
}

static int read_bonds(void *v, int *nbonds, int **fromptr, int **toptr, 
                      float ** bondorder,int **bondtype, 
                      int *nbondtypes, char ***bondtypename) {
  pdbdata *pdb = (pdbdata *)v;
  
  *nbonds = 0;
  *fromptr = NULL;
  *toptr = NULL;
  *bondorder = NULL; /* PDB files don't have bond order information */
  *bondtype = NULL;
  *nbondtypes = 0;
  *bondtypename = NULL;

/* The newest plugin API allows us to return CONECT records as 
 * additional bonds above and beyond what the distance search returns.
 * Without that feature, we otherwise have to check completeness and
 * ignore them if they don't look to be fully specified for this molecule */
#if !defined(MOLFILE_BONDSSPECIAL)
  if (pdb->natoms >= 100000) {
    printf("pdbplugin) Warning: more than 99,999 atoms, ignored CONECT records\n");
    return MOLFILE_SUCCESS;
  } else if (((float) pdb->nconect / (float) pdb->natoms) <= 0.85) {
    printf("pdbplugin) Warning: Probable incomplete bond structure specified,\n");
    printf("pdbplugin)          ignoring CONECT records\n");
    return MOLFILE_SUCCESS;
  } else if (pdb->nconect == 0) {
    return MOLFILE_SUCCESS;
  }
#endif

  *nbonds = pdb->nbonds;
  *fromptr = pdb->from;
  *toptr = pdb->to;

  return MOLFILE_SUCCESS;
}


/* 
 * 
 */
static int read_next_timestep(void *v, int natoms, molfile_timestep_t *ts) {
  pdbdata *pdb = (pdbdata *)v;
  char pdbstr[PDB_BUFFER_LENGTH];
  int indx, i;
  float *x, *y, *z;
  float occup, bfac;
  if (pdb->natoms == 0) 
    return MOLFILE_ERROR; /* EOF */
  if (ts) {
    x = ts->coords;
    y = x+1;
    z = x+2;
  } else {
    x = y = z = 0;
  } 
  i = 0;
  do {
    indx = read_pdb_record(pdb->fd, pdbstr);
    if((indx == PDB_END || indx == PDB_EOF) && (i < pdb->natoms)) {
      return MOLFILE_ERROR;
    } else if(indx == PDB_ATOM) {
      if(i++ >= pdb->natoms) {
        break;      
      }
      /* just get the coordinates, and store them */
      if (ts) {
        get_pdb_coordinates(pdbstr, x, y, z, &occup, &bfac);
        x += 3;
        y += 3;
        z += 3;
      } 
    } else if (indx == PDB_CRYST1) {
      if (ts) {
        get_pdb_cryst1(pdbstr, &ts->alpha, &ts->beta, &ts->gamma,
                               &ts->A, &ts->B, &ts->C);
      }
    }
  } while(!(indx == PDB_END || indx == PDB_EOF));

  return MOLFILE_SUCCESS;
}

static void close_pdb_read(void *v) { 
  pdbdata *pdb = (pdbdata *)v;
  if (pdb->fd != NULL)
    fclose(pdb->fd);
  if (pdb->idxmap != NULL)
    free(pdb->idxmap);
  if (pdb->meta->remarks != NULL)
    free(pdb->meta->remarks);
  if (pdb->meta != NULL) 
    free(pdb->meta);
  free(pdb);
}

static void *open_file_write(const char *path, const char *filetype, 
    int natoms) {

  FILE *fd;
  pdbdata *pdb;
  fd = fopen(path, "w");
  if (!fd) {
    fprintf(stderr, "Unable to open file %s for writing\n", path);
    return NULL;
  }
  pdb = (pdbdata *)malloc(sizeof(pdbdata));
  pdb->fd = fd;
  pdb->natoms = natoms; 
  pdb->atomlist = NULL;
  pdb->first_frame = 1;
  return pdb;
}
 
static int write_structure(void *v, int optflags, 
    const molfile_atom_t *atoms) {

  int i;
  pdbdata *pdb = (pdbdata *)v;
  int natoms = pdb->natoms;
  pdb->atomlist = (molfile_atom_t *)malloc(natoms*sizeof(molfile_atom_t));
  memcpy(pdb->atomlist, atoms, natoms*sizeof(molfile_atom_t));

  /* If occ, bfactor, and insertion aren't given, we assign defaultvalues. */
  if (!(optflags & MOLFILE_OCCUPANCY)) {
    for (i=0; i<natoms; i++) pdb->atomlist[i].occupancy = 0.0f;
  }
  if (!(optflags & MOLFILE_BFACTOR)) {
    for (i=0; i<natoms; i++) pdb->atomlist[i].bfactor= 0.0f;
  }
  if (!(optflags & MOLFILE_INSERTION)) {
    for (i=0; i<natoms; i++) {
      pdb->atomlist[i].insertion[0] =' ';
      pdb->atomlist[i].insertion[1] ='\0';
    }
  }
  if (!(optflags & MOLFILE_ALTLOC)) {
    for (i=0; i<natoms; i++) {
      pdb->atomlist[i].altloc[0]=' ';
      pdb->atomlist[i].altloc[1]='\0';
    }
  }
  if (!(optflags & MOLFILE_ATOMICNUMBER)) {
    for (i=0; i<natoms; i++) pdb->atomlist[i].atomicnumber = 0;
  }

  /* TODO: put bonds into CONECT records? */
  return MOLFILE_SUCCESS;
}

/* SEQRES records look like this:

COLUMNS        DATA TYPE       FIELD         DEFINITION
---------------------------------------------------------------------------------
 1 -  6        Record name     "SEQRES"

 9 - 10        Integer         serNum        Serial number of the SEQRES record
                                             for the current chain.  Starts at 1
                                             and increments by one each line.
                                             Reset to 1 for each chain.

12             Character       chainID       Chain identifier.  This may be any
                                             single legal character, including a
                                             blank which is used if there is
                                             only one chain.

14 - 17        Integer         numRes        Number of residues in the chain.
                                             This value is repeated on every
                                             record.

20 - 22        Residue name    resName       Residue name.

24 - 26        Residue name    resName       Residue name.

... and so forth out to 68-70, for a total of 13 in each line (except possibly
the last.

source:
http:

*/

/*
 * However, we don't use them right now because of several issues that
 * can't presently be resolved satisfactorily in VMD:

According to the RCSB, SEQRES records have to contain all residues, not
just those in the structure, which means VMD will usually produce incorrect
output and there's nothing we can do about it.  The RCSB actually specifies
that all residues in the chain have to present in the SEQRES records, even
if they're not in the structure.
  
We can never know which residues to output.  Our current system of outputting   
everything is just terrible when you have 20,000 waters in your system; we
have to fix this immediately.  We could almost get away with making a hash
table of the names of protein and nucleic acid residues and only write chains
containing those residues.  However, there's this little snippet from the
specification:
  
* Heterogens which are integrated into the backbone of the chain are listed
  as being part of the chain and are included in the SEQRES records for
  that chain.
  
That means that we can never know what might appear in the sequence unless we
also read HET records and keep track of them in VMD as well.  We shouldn't 
get people depending on such fallible SEQRES records.
  
And of course, there's the fact that no other program that we know of besides   
CE needs these SEQRES records.

 * Uncomment the write_seqres line in write_timestep to turn them back on.
 */


#if 0
static void write_seqres(FILE * fd, int natoms, const molfile_atom_t *atomlist) {
  int i=0;
  while (i < natoms) {
    int k, serNum;
    int j = i;
    int ires, nres = 1;
    int resid = atomlist[i].resid;
    /* Count up the number of residues in the chain */
    const char *chain = atomlist[i].chain;
    while (j < natoms && !strcmp(chain, atomlist[j].chain)) {
      if (resid != atomlist[j].resid) {
        nres++;
        resid = atomlist[j].resid;
      }
      j++;
    }
    /* There are nres residues in the chain, from atoms i to j. */
    serNum = 1;
    ires = 1;
    resid = atomlist[i].resid;
    fprintf(fd, "SEQRES  %2d %c %4d  ",  serNum, chain[0], nres);
    serNum = 2;
    fprintf(fd, "%3s ", atomlist[i].resname);
    for (k=i; k<j; k++) {
      if (resid != atomlist[k].resid) {
        resid = atomlist[k].resid;
        if (!(ires % 13)) {
          fprintf(fd, "\nSEQRES  %2d %c %4d  ",  serNum, chain[0], nres);
          serNum++;
        }
        fprintf(fd, "%3s ", atomlist[k].resname);
        ires++;
      }
    }
    i = j;
    fprintf(fd, "\n");
  }
}
#endif

/*
CRYST1 records look like this:
The CRYST1 record presents the unit cell parameters, space group, and Z value. If the structure was not determined by crystallographic means, CRYST1 simply defines a unit cube. 


Record Format 

COLUMNS       DATA TYPE      FIELD         DEFINITION
-------------------------------------------------------------
 1 -  6       Record name    "CRYST1"

 7 - 15       Real(9.3)      a             a (Angstroms).

16 - 24       Real(9.3)      b             b (Angstroms).

25 - 33       Real(9.3)      c             c (Angstroms).

34 - 40       Real(7.2)      alpha         alpha (degrees).

41 - 47       Real(7.2)      beta          beta (degrees).

48 - 54       Real(7.2)      gamma         gamma (degrees).

56 - 66       LString        sGroup        Space group.

67 - 70       Integer        z             Z value.

* If the coordinate entry describes a structure determined by a technique
other than crystallography, CRYST1 contains a = b = c = 1.0, alpha =
beta = gamma = 90 degrees, space group = P 1, and Z = 1.

We will use "P 1" and "1" for space group and z value, as recommended, but
we'll populate the other fields with the unit cell information we do have.

*/
  
static void write_cryst1(FILE *fd, const molfile_timestep_t *ts) {
  fprintf(fd, "CRYST1%9.3f%9.3f%9.3f%7.2f%7.2f%7.2f P 1           1\n", 
    ts->A, ts->B, ts->C, ts->alpha, ts->beta, ts->gamma);
}


static int write_timestep(void *v, const molfile_timestep_t *ts) {
  pdbdata *pdb = (pdbdata *)v; 
  const molfile_atom_t *atom;
  const float *pos;
  int i;
  char elementsymbol[3];

  if (pdb->natoms == 0)
    return MOLFILE_SUCCESS;

  if (pdb->first_frame) {
    /* Turn off SEQRES writing for now; see comments above.
    write_seqres(pdb->fd, pdb->natoms, pdb->atomlist);
    */
    write_cryst1(pdb->fd, ts);
    pdb->first_frame = 0;
  }
  atom = pdb->atomlist;
  pos = ts->coords;
  for (i=0; i<pdb->natoms; i++) {
    /*
     * The 8.3 format for position, occupancy, and bfactor permits values 
     * only in the range of -999.9994 to 9999.9994 (so that they round
     * to the range [-999.999, 9999.999]).  If values fall outside of that
     * range, fail and emit an error message rather than generate a
     * misformatted PDB file.
     */
#define PDBBAD(x) ((x) < -999.9994f || (x) > 9999.9994f)
    if (PDBBAD(pos[0]) || PDBBAD(pos[1]) || PDBBAD(pos[2]) ||
		PDBBAD(atom->occupancy) || PDBBAD(atom->bfactor)) {
	    fprintf(stderr, "PDB WRITE ERROR: Position, occupancy, or b-factor (beta) for atom %d\n", i);
      fprintf(stderr, "                 cannot be written in PDB format.\n");
      fprintf(stderr, "                 File will be truncated.\n");
      return MOLFILE_ERROR;
    }

    /* check the atomicnumber and format the atomic element symbol string */
    strcpy(elementsymbol, (atom->atomicnumber < 1) ? "  " : get_pte_label(atom->atomicnumber));
    elementsymbol[0] = toupper(elementsymbol[0]);
    elementsymbol[1] = toupper(elementsymbol[1]);
 
    if (!write_raw_pdb_record(pdb->fd,  
        "ATOM  ", i+1, atom->name, atom->resname, atom->resid, 
        atom->insertion, atom->altloc, elementsymbol,
        pos[0], pos[1], pos[2], 
        atom->occupancy, atom->bfactor, atom->chain, atom->segid)) {
      fprintf(stderr, 
          "PDB: Error encoutered writing atom %d; file may be incomplete.\n", 
          i+1);
      return MOLFILE_ERROR;
    }
    ++atom;
    pos += 3;
  }
  fprintf(pdb->fd, "END\n");

  return MOLFILE_SUCCESS;
}
 
static void close_file_write(void *v) {
  pdbdata *pdb = (pdbdata *)v; 
  fclose(pdb->fd);
  free(pdb->atomlist);
  free(pdb);
}

static int read_molecule_metadata(void *v, molfile_metadata_t **metadata) {
  pdbdata *pdb = (pdbdata *)v; 
  *metadata = pdb->meta;
  return MOLFILE_SUCCESS;
}

/*
 * Initialization stuff down here
 */

static molfile_plugin_t plugin;
 
VMDPLUGIN_API int VMDPLUGIN_init() {
  memset(&plugin, 0, sizeof(molfile_plugin_t));
  plugin.abiversion = vmdplugin_ABIVERSION;
  plugin.type = MOLFILE_PLUGIN_TYPE;
  plugin.name = "pdb";
  plugin.prettyname = "PDB";
  plugin.author = "Justin Gullingsrud, John Stone";
  plugin.majorv = 1;
  plugin.minorv = 16;
  plugin.is_reentrant = VMDPLUGIN_THREADSAFE;
  plugin.filename_extension = "pdb,ent";
  plugin.open_file_read = open_pdb_read;
  plugin.read_structure = read_pdb_structure;
  plugin.read_bonds = read_bonds;
  plugin.read_next_timestep = read_next_timestep;
  plugin.close_file_read = close_pdb_read;
  plugin.open_file_write = open_file_write;
  plugin.write_structure = write_structure;
  plugin.write_timestep = write_timestep;
  plugin.close_file_write = close_file_write;
  plugin.read_molecule_metadata = read_molecule_metadata;
  return VMDPLUGIN_SUCCESS;
}

VMDPLUGIN_API int VMDPLUGIN_register(void *v, vmdplugin_register_cb cb) {
  (*cb)(v, (vmdplugin_t *)&plugin);
  return VMDPLUGIN_SUCCESS;
}

VMDPLUGIN_API int VMDPLUGIN_fini() {
  return VMDPLUGIN_SUCCESS;
}

