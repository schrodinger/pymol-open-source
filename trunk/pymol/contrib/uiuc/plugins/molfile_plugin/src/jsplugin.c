/* MACHINE GENERATED FILE, DO NOT EDIT! */

#define VMDPLUGIN molfile_jsplugin
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
 *      $RCSfile: jsplugin.c,v $
 *      $Author: johns $       $Locker:  $             $State: Exp $
 *      $Revision: 1.34 $       $Date: 2009/04/29 15:45:30 $
 *
 ***************************************************************************
 * DESCRIPTION:
 *   Code for reading and writing John's molecular dynamics trajectory files
 *   for I/O performance testing and as a simple example plugin.
 *
 *   No programs actually use this format, so it's only a test/example code.
 *
 *   This plugin reads and writes mostly the same data that a DCD file 
 *   contains, but the ordering of I/O operations avoids having to do
 *   scatter/gather passes on the buffers, so it ends up performing 
 *   two to three times faster than the DCD format for the same data set.
 *
 *   Best measured I/O performance in VMD so far is 631 MB/sec on a Sun V880z 
 *   Best standalone I/O performance so far is 688 MB/sec.
 *
 *  Standalone test binary compilation flags:
 *  cc -fast -xarch=v9a -I../../include -DTEST_JSPLUGIN jsplugin.c \
 *    -o ~/bin/readjs -lm
 *
 *  Profiling flags:
 *  cc -xpg -fast -xarch=v9a -g -I../../include -DTEST_JSPLUGIN jsplugin.c \
 *    -o ~/bin/readjs -lm
 *
 ***************************************************************************/

#include "largefiles.h"   /* platform dependent 64-bit file I/O defines */

#include <sys/stat.h>
#include <sys/types.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define VMDPLUGIN_STATIC
#include "hash.h"
#include "fastio.h"
#include "endianswap.h"
#include "molfile_plugin.h"

#ifndef M_PI_2
#define M_PI_2 1.57079632679489661922
#endif

#define JSHEADERSTRING   "JS Binary Structure and Trajectory File Format"                
#define JSMAGICNUMBER    0x31337
#define JSENDIANISM      0x12345678

#define JSMAJORVERSION   2
#define JSMINORVERSION   3

#define JSNFRAMESOFFSET  (strlen(JSHEADERSTRING) + 20)

#define JSNOERR             0
#define JSBADFILE           1
#define JSBADFORMAT         2


#define JSOPT_NOOPTIONS     0x00000000

#define JSOPT_STRUCTURE     0x00000001
#define JSOPT_BONDS         0x00000002
#define JSOPT_BONDORDERS    0x00000004
#define JSOPT_ANGLES        0x00000008
#define JSOPT_CTERMS        0x00000010

#define JSOPT_OCCUPANCY     0x00000100
#define JSOPT_BFACTOR       0x00000200
#define JSOPT_MASS          0x00000400
#define JSOPT_CHARGE        0x00000800
#define JSOPT_RADIUS        0x00001000
#define JSOPT_ATOMICNUMBER  0x00002000

typedef struct {
  fio_fd fd;
  int natoms;

#if JSMAJORVERSION > 1
  /* structure info */
  int optflags;
  molfile_atom_t *atomlist;
  molfile_metadata_t *meta;

  /* bond info */
  int nbonds;
  int *bondfrom;
  int *bondto;
  float *bondorders;

  /* angle/dihedral/improper/cross-term info */
  int numangles, *angles;
  int numdihedrals, *dihedrals;
  int numimpropers, *impropers;
  int numcterms, *cterms;
#endif

  /* trajectory info */
  int nframes;
  double tsdelta;
  int reverseendian;
  int with_unitcell;

} jshandle;

static void *open_js_read(const char *path, const char *filetype, int *natoms) {
  jshandle *js;
  int jsmagicnumber, jsendianism, jsmajorversion, jsminorversion;
  struct stat stbuf;
  char strbuf[1024];
  int rc = 0;

  if (!path) return NULL;

  /* See if the file exists, and get its size */
  memset(&stbuf, 0, sizeof(struct stat));
  if (stat(path, &stbuf)) {
    printf("jsplugin) Could not access file '%s'.\n", path);
    return NULL;
  }

  js = (jshandle *)malloc(sizeof(jshandle));
  memset(js, 0, sizeof(jshandle));
  if (fio_open(path, FIO_READ, &js->fd) < 0) {
    printf("jsplugin) Could not open file '%s' for reading.\n", path);
    free(js);
    return NULL;
  }

  /* emit header information */
  fio_fread(strbuf, strlen(JSHEADERSTRING), 1, js->fd);
  strbuf[strlen(JSHEADERSTRING)] = '\0';
  if (strcmp(strbuf, JSHEADERSTRING)) {
    printf("jsplugin) Bad trajectory header!\n");
    printf("jsplugin) Read string: %s\n", strbuf);
    return NULL;
  }

  fio_read_int32(js->fd, &jsmagicnumber);
  fio_read_int32(js->fd, &jsendianism);
  fio_read_int32(js->fd, &jsmajorversion);
  fio_read_int32(js->fd, &jsminorversion);
  fio_read_int32(js->fd, &js->natoms);
  fio_read_int32(js->fd, &js->nframes);
  if ((jsmagicnumber != JSMAGICNUMBER) || (jsendianism != JSENDIANISM)) {
    printf("jsplugin) opposite endianism file, enabling byte swapping\n");
    js->reverseendian = 1;
    swap4_aligned(&jsmagicnumber, 1);
    swap4_aligned(&jsendianism, 1);
    swap4_aligned(&jsmajorversion, 1);
    swap4_aligned(&jsminorversion, 1);
    swap4_aligned(&js->natoms, 1);
    swap4_aligned(&js->nframes, 1);
  } else {
    printf("jsplugin) native endianism file\n");
  }

  if ((jsmagicnumber != JSMAGICNUMBER) || (jsendianism != JSENDIANISM)) {
    printf("jsplugin) read_jsreader returned %d\n", rc);
    fio_fclose(js->fd);
    free(js);
    return NULL;
  }
 
  if (jsmajorversion != JSMAJORVERSION) {
    printf("jsplugin) major version mismatch\n");
    printf("jsplugin)   file version: %d\n", jsmajorversion);
    printf("jsplugin)   plugin version: %d\n", JSMAJORVERSION);
    fio_fclose(js->fd);
    free(js);
    return NULL;
  }
 
  *natoms = js->natoms;
  return js;
}


#if JSMAJORVERSION > 1

static int read_js_structure(void *mydata, int *optflags,
                             molfile_atom_t *atoms) {
  jshandle *js = (jshandle *) mydata;
  int i;

  *optflags = MOLFILE_NOOPTIONS; /* set to no options until we read them */

  /* write flags data to the file */
  fio_read_int32(js->fd, &js->optflags); 
  if (js->reverseendian)
    swap4_aligned(&js->optflags, 1);
printf("jsplugin) read option flags: %0x08x\n", js->optflags);

  /* determine whether or not this file contains structure info or not */
  if (js->optflags & JSOPT_STRUCTURE) {
    int numatomnames, numatomtypes, numresnames, numsegids, numchains;
    char **atomnames = NULL;
    char **atomtypes = NULL;
    char **resnames = NULL;
    char **segids = NULL;
    char **chains = NULL;
    short *shortbuf = NULL; /* temp buf for decoding atom records */
    int *intbuf = NULL;     /* temp buf for decoding atom records */
    float *fltbuf = NULL;   /* temp buf for decoding atom records */
 
    /* read in block of name string table sizes */
    fio_read_int32(js->fd, &numatomnames); 
    fio_read_int32(js->fd, &numatomtypes); 
    fio_read_int32(js->fd, &numresnames);
    fio_read_int32(js->fd, &numsegids);
    fio_read_int32(js->fd, &numchains); 
    if (js->reverseendian) {
      swap4_aligned(&numatomnames, js->natoms);
      swap4_aligned(&numatomtypes, js->natoms);
      swap4_aligned(&numresnames, js->natoms);
      swap4_aligned(&numsegids, js->natoms);
      swap4_aligned(&numchains, js->natoms);
    }

printf("jsplugin) reading string tables...\n");
printf("jsplugin) %d %d %d %d %d\n",
       numatomnames, numatomtypes, numresnames, numsegids, numchains);

    /* allocate string tables */
    atomnames = (char **) malloc(numatomnames * sizeof(char *));
    atomtypes = (char **) malloc(numatomtypes * sizeof(char *));
    resnames  = (char **) malloc(numresnames  * sizeof(char *));
    segids    = (char **) malloc(numsegids    * sizeof(char *));
    chains    = (char **) malloc(numchains    * sizeof(char *));

printf("jsplugin)   atom names...\n");
    /* read in the string tables */
    for (i=0; i<numatomnames; i++) {
      atomnames[i] = (char *) malloc(16 * sizeof(char));
      fio_fread(atomnames[i], 16 * sizeof(char), 1, js->fd);
    }

printf("jsplugin)   atom types...\n");
    for (i=0; i<numatomtypes; i++) {
      atomtypes[i] = (char *) malloc(16 * sizeof(char));
      fio_fread(atomtypes[i], 16 * sizeof(char), 1, js->fd);
    }

printf("jsplugin)   residue names...\n");
    for (i=0; i<numresnames; i++) {
      resnames[i] = (char *) malloc(8 * sizeof(char));
      fio_fread(resnames[i], 8 * sizeof(char), 1, js->fd);
    }

printf("jsplugin)   segment names...\n");
    for (i=0; i<numsegids; i++) {
      segids[i] = (char *) malloc(8 * sizeof(char));
      fio_fread(segids[i], 8 * sizeof(char), 1, js->fd);
    }

printf("jsplugin)   chain names...\n");
    for (i=0; i<numchains; i++) {
      chains[i] = (char *) malloc(2 * sizeof(char));
      fio_fread(chains[i], 2 * sizeof(char), 1, js->fd);
    }

printf("jsplugin) reading numeric field tables...\n");
    /* read in all of the atom fields */
    shortbuf = (void *) malloc(js->natoms * sizeof(short));

printf("jsplugin)   atom name indices...\n");
    /* read in atom names */
    fio_fread(shortbuf, js->natoms * sizeof(short), 1, js->fd);
    if (js->reverseendian)
      swap2_aligned(shortbuf, js->natoms);
    for (i=0; i<js->natoms; i++) {
      strcpy(atoms[i].name, atomnames[shortbuf[i]]);
    }    

printf("jsplugin)   atom type indices...\n");
    /* read in atom types */
    fio_fread(shortbuf, js->natoms * sizeof(short), 1, js->fd);
    if (js->reverseendian)
      swap2_aligned(shortbuf, js->natoms);
    for (i=0; i<js->natoms; i++) {
      strcpy(atoms[i].type, atomtypes[shortbuf[i]]);
    }    

printf("jsplugin)   residue name indices...\n");
    /* read in resnames */
    fio_fread(shortbuf, js->natoms * sizeof(short), 1, js->fd);
    if (js->reverseendian)
      swap2_aligned(shortbuf, js->natoms);
    for (i=0; i<js->natoms; i++) {
      strcpy(atoms[i].resname, resnames[shortbuf[i]]);
    }    
    
printf("jsplugin)   segment name indices...\n");
    /* read in segids */
    fio_fread(shortbuf, js->natoms * sizeof(short), 1, js->fd);
    if (js->reverseendian)
      swap2_aligned(shortbuf, js->natoms);
    for (i=0; i<js->natoms; i++) {
      strcpy(atoms[i].segid, segids[shortbuf[i]]);
    }    

printf("jsplugin)   chain name indices...\n");
    /* read in chains */
    fio_fread(shortbuf, js->natoms * sizeof(short), 1, js->fd);
    if (js->reverseendian)
      swap2_aligned(shortbuf, js->natoms);
    for (i=0; i<js->natoms; i++) {
      strcpy(atoms[i].chain, chains[shortbuf[i]]);
    }    

    if (shortbuf != NULL) {
      free(shortbuf);
      shortbuf=NULL;
    }

    /* 
     * read in integer data blocks 
     */
    intbuf = (int *) malloc(js->natoms * sizeof(int));

printf("jsplugin)   residue indices...\n");
    /* read in resid */
    fio_fread(intbuf, js->natoms * sizeof(int), 1, js->fd);
    if (js->reverseendian)
      swap4_aligned(intbuf, js->natoms);
    for (i=0; i<js->natoms; i++) {
      atoms[i].resid = intbuf[i];
    }    
     
    if (intbuf != NULL) {
      free(intbuf);
      intbuf = NULL;
    }


printf("jsplugin) reading optional per-atom tables...\n");
    /*
     * read in optional single-precision float data blocks
     */ 
    if (js->optflags & (JSOPT_OCCUPANCY | JSOPT_BFACTOR | 
        JSOPT_MASS | JSOPT_RADIUS | JSOPT_CHARGE)) 
      fltbuf = (void *) malloc(js->natoms * sizeof(float));

    /* read in optional data if it exists */
    if (js->optflags & JSOPT_OCCUPANCY) {
printf("jsplugin)   occupancy...\n");
      *optflags |= MOLFILE_OCCUPANCY;
      fio_fread(fltbuf, js->natoms * sizeof(float), 1, js->fd);
      if (js->reverseendian)
        swap4_aligned(fltbuf, js->natoms);
      for (i=0; i<js->natoms; i++) {
        atoms[i].occupancy = fltbuf[i];
      }    
    }

    if (js->optflags & JSOPT_BFACTOR) {
printf("jsplugin)   bfactor...\n");
      *optflags |= MOLFILE_BFACTOR;
      fio_fread(fltbuf, js->natoms * sizeof(float), 1, js->fd);
      if (js->reverseendian)
        swap4_aligned(fltbuf, js->natoms);
      for (i=0; i<js->natoms; i++) {
        atoms[i].bfactor = fltbuf[i];
      }    
    }

    if (js->optflags & JSOPT_MASS) { 
printf("jsplugin)   mass...\n");
      *optflags |= MOLFILE_MASS;
      fio_fread(fltbuf, js->natoms * sizeof(float), 1, js->fd);
      if (js->reverseendian)
        swap4_aligned(fltbuf, js->natoms);
      for (i=0; i<js->natoms; i++) {
        atoms[i].mass = fltbuf[i];
      }    
    }

    if (js->optflags & JSOPT_CHARGE) { 
printf("jsplugin)   charge...\n");
      *optflags |= MOLFILE_CHARGE;
      fio_fread(fltbuf, js->natoms * sizeof(float), 1, js->fd);
      if (js->reverseendian)
        swap4_aligned(fltbuf, js->natoms);
      for (i=0; i<js->natoms; i++) {
        atoms[i].charge = fltbuf[i];
      }    
    }

    if (js->optflags & JSOPT_RADIUS) { 
printf("jsplugin)   radius...\n");
      *optflags |= MOLFILE_RADIUS;
      fio_fread(fltbuf, js->natoms * sizeof(float), 1, js->fd);
      if (js->reverseendian)
        swap4_aligned(fltbuf, js->natoms);
      for (i=0; i<js->natoms; i++) {
        atoms[i].radius = fltbuf[i];
      }    
    }

    if (fltbuf != NULL) {
      free(fltbuf);
      fltbuf=NULL;
    }

    /*
     * read in optional integer data blocks
     */ 
    if (js->optflags & JSOPT_ATOMICNUMBER)
      intbuf = (void *) malloc(js->natoms * sizeof(int));

    if (js->optflags & JSOPT_ATOMICNUMBER) { 
printf("jsplugin)   atomic number...\n");
      *optflags |= MOLFILE_ATOMICNUMBER;
      fio_fread(intbuf, js->natoms * sizeof(int), 1, js->fd);
      if (js->reverseendian)
        swap4_aligned(intbuf, js->natoms);
      for (i=0; i<js->natoms; i++) {
        atoms[i].atomicnumber = intbuf[i];
      }    
    }

    if (intbuf != NULL) {
      free(intbuf);
      intbuf = NULL;
    }


    /*
     * read in bonds and fractional bond orders
     */ 
    if (js->optflags & JSOPT_BONDS) {
      fio_fread(&js->nbonds, sizeof(int), 1, js->fd);
      if (js->reverseendian)
        swap4_aligned(&js->nbonds, 1);
printf("jsplugin)   %d bonds...\n", js->nbonds);

      js->bondfrom = (int *) malloc(js->nbonds * sizeof(int));
      js->bondto = (int *) malloc(js->nbonds * sizeof(int));
      fio_fread(js->bondfrom, js->nbonds * sizeof(int), 1, js->fd);
      fio_fread(js->bondto, js->nbonds * sizeof(int), 1, js->fd);
      if (js->reverseendian) {
        swap4_aligned(js->bondfrom, js->nbonds);
        swap4_aligned(js->bondto, js->nbonds);
      }

      if (js->optflags & JSOPT_BONDORDERS) {
printf("jsplugin)   bond orders...\n");
        js->bondorders = (void *) malloc(js->nbonds * sizeof(float));
        fio_fread(js->bondorders, js->nbonds * sizeof(float), 1, js->fd);
        if (js->reverseendian)
          swap4_aligned(js->bondorders, js->nbonds);
      }
    }

    if (js->optflags & JSOPT_ANGLES) {
      fio_fread(&js->numangles, sizeof(int), 1, js->fd);
      if (js->reverseendian)
        swap4_aligned(&js->numangles, 1);
printf("jsplugin)   %d angles...\n", js->numangles);
      js->angles = (int *) malloc(3 * js->numangles * sizeof(int));
      fio_fread(js->angles, sizeof(int)*3*js->numangles, 1, js->fd);
      if (js->reverseendian)
        swap4_aligned(&js->angles, 3*js->numangles);

      fio_fread(&js->numdihedrals, sizeof(int), 1, js->fd);
      if (js->reverseendian)
        swap4_aligned(&js->numdihedrals, 1);
printf("jsplugin)   %d dihedrals...\n", js->numdihedrals);
      js->dihedrals = (int *) malloc(4 * js->numdihedrals * sizeof(int));
      fio_fread(js->dihedrals, sizeof(int)*4*js->numdihedrals, 1, js->fd);
      if (js->reverseendian)
        swap4_aligned(&js->dihedrals, 4*js->numdihedrals);

      fio_fread(&js->numimpropers, sizeof(int), 1, js->fd);
      if (js->reverseendian)
        swap4_aligned(&js->numimpropers, 1);
      js->impropers = (int *) malloc(4 * js->numimpropers * sizeof(int));
printf("jsplugin)   %d impropers...\n", js->numimpropers);
      fio_fread(js->impropers, sizeof(int)*4*js->numimpropers, 1, js->fd);
      if (js->reverseendian)
        swap4_aligned(&js->impropers, 4*js->numimpropers);
    }
    if (js->optflags & JSOPT_CTERMS) {
      fio_fread(&js->numcterms, sizeof(int), 1, js->fd);
      if (js->reverseendian)
        swap4_aligned(&js->numcterms, 1);
      js->cterms = (int *) malloc(8 * js->numcterms * sizeof(int));
printf("jsplugin)   %d cterms...\n", js->numcterms);
      fio_fread(js->cterms, sizeof(int)*8*js->numcterms, 1, js->fd);
      if (js->reverseendian)
        swap4_aligned(&js->cterms, 8*js->numcterms);
    }

printf("jsplugin) final optflags: %08x\n", *optflags);
printf("jsplugin) structure information complete\n");

    return MOLFILE_SUCCESS;
  }

printf("jsplugin) no structure information available\n");

  /* else, we have no structure information */
  return MOLFILE_NOSTRUCTUREDATA;
}


static int read_js_bonds(void *v, int *nbonds, int **fromptr, int **toptr, 
                         float **bondorder, int **bondtype, 
                         int *nbondtypes, char ***bondtypename) {
  jshandle *js = (jshandle *)v;

  *nbonds = 0;
  *fromptr = NULL;
  *toptr = NULL;
  *bondorder = NULL;
  *bondtype = NULL;
  *nbondtypes = 0;
  *bondtypename = NULL;

  if (js->optflags & JSOPT_BONDS) {
    /* save bond info until we actually write out the structure file */
    *nbonds = js->nbonds;
    *fromptr = js->bondfrom;
    *toptr = js->bondto;

    if (js->optflags & JSOPT_BONDORDERS) {
      *bondorder = js->bondorders;
    }
  }

  return MOLFILE_SUCCESS;
}

#if vmdplugin_ABIVERSION > 14
static int read_js_angles(void *v, int *numangles, int **angles, 
                          int **angletypes, int *numangletypes, 
                          char ***angletypenames, int *numdihedrals,
                          int **dihedrals, int **dihedraltypes, 
                          int *numdihedraltypes, char ***dihedraltypenames,
                          int *numimpropers, int **impropers, 
                          int **impropertypes, int *numimpropertypes, 
                          char ***impropertypenames, int *numcterms, 
                          int **cterms, int *ctermcols, int *ctermrows) {
  jshandle *js = (jshandle *)v;

  /* initialize data to zero */
  *numangles         = 0;
  *angles            = NULL;
  *angletypes        = NULL;
  *numangletypes     = 0;
  *angletypenames    = NULL;
  *numdihedrals      = 0;
  *dihedrals         = NULL;
  *dihedraltypes     = NULL;
  *numdihedraltypes  = 0;
  *dihedraltypenames = NULL;
  *numimpropers      = 0;
  *impropers         = NULL;
  *impropertypes     = NULL;
  *numimpropertypes  = 0;
  *impropertypenames = NULL;
  *numcterms         = 0;
  *cterms            = NULL;
  *ctermrows         = 0;
  *ctermcols         = 0;

  *numangles = js->numangles;
  *angles = js->angles;

  *numdihedrals = js->numdihedrals;
  *dihedrals = js->dihedrals;

  *numimpropers = js->numimpropers;
  *impropers = js->impropers;

  *numcterms = js->numcterms;
  *cterms = js->cterms;
  *ctermcols = 0;
  *ctermrows = 0;

  return MOLFILE_SUCCESS;
}
#else
static int read_js_angles(void *v,
               int *numangles,    int **angles,    double **angleforces,
               int *numdihedrals, int **dihedrals, double **dihedralforces,
               int *numimpropers, int **impropers, double **improperforces,
               int *numcterms,    int **cterms,
               int *ctermcols,    int *ctermrows,  double **ctermforces) {
  jshandle *js = (jshandle *)v;

  *numangles = js->numangles;
  *angles = js->angles;
  *angleforces = NULL;

  *numdihedrals = js->numdihedrals;
  *dihedrals = js->dihedrals;
  *dihedralforces = NULL;

  *numimpropers = js->numimpropers;
  *impropers = js->impropers;
  *improperforces = NULL;

  *numcterms = js->numcterms;
  *cterms = js->cterms;
  *ctermcols = 0;
  *ctermrows = 0;
  *ctermforces = NULL;

  return MOLFILE_SUCCESS;
}
#endif

#endif


static int read_js_timestep(void *v, int natoms, molfile_timestep_t *ts) {
  jshandle *js = (jshandle *)v;
  int framelen = (natoms * 3 * sizeof(float)) + (6 * sizeof(double));;

  /* if we have a valid ts pointer, read the timestep, otherwise skip it */ 
  if (ts != NULL) {
    int readlen; 
    fio_iovec iov[2];
    double unitcell[6];
    unitcell[0] = unitcell[2] = unitcell[5] = 1.0f;
    unitcell[1] = unitcell[3] = unitcell[4] = 90.0f;
  
    /* setup the I/O vector */
    iov[0].iov_base = (fio_caddr_t) ts->coords;   /* read coordinates    */
    iov[0].iov_len  = natoms * 3 * sizeof(float);
    iov[1].iov_base = (fio_caddr_t) &unitcell[0]; /* read PBC unit cell  */
    iov[1].iov_len  = 6 * sizeof(double);

    /* Do all of the reads with a single syscall, for peak efficiency. */
    /* On smart kernels, readv() causes only one context switch, and   */
    /* can effeciently scatter the reads to the various buffers.       */
    readlen = fio_readv(js->fd, &iov[0], 2); 
  
    /* check the number of read bytes versus what we expected */
    if (readlen != framelen)
      return MOLFILE_EOF;

    /* perform byte swapping if necessary */
    if (js->reverseendian) {
      swap4_aligned(ts->coords, natoms * 3);
      swap8_aligned(unitcell, 6);
    }

    /* copy unit cell values into VMD */
    ts->A = unitcell[0];
    ts->B = unitcell[1];
    ts->C = unitcell[2];
    ts->alpha = 90.0 - asin(unitcell[3]) * 90.0 / M_PI_2;
    ts->beta  = 90.0 - asin(unitcell[4]) * 90.0 / M_PI_2;
    ts->gamma = 90.0 - asin(unitcell[5]) * 90.0 / M_PI_2;
  } else {
    /* skip this frame, seek to the next frame */
    if (fio_fseek(js->fd, framelen, FIO_SEEK_CUR)) 
      return MOLFILE_EOF;
  }
 
  return MOLFILE_SUCCESS;
}
 

static void close_js_read(void *v) {
  jshandle *js = (jshandle *)v;
  fio_fclose(js->fd);

#if JSMAJORVERSION > 1
  if (js->bondfrom)
    free(js->bondfrom);
  if (js->bondto)
    free(js->bondto);
  if (js->bondorders)
    free(js->bondorders);

  /* free angle data */
  if (js->angles != NULL)
    free(js->angles);
  if (js->dihedrals != NULL)
    free(js->dihedrals);
  if (js->impropers != NULL)
    free(js->impropers);
  if (js->cterms)
    free(js->cterms);
#endif

  free(js);
}


static void *open_js_write(const char *path, const char *filetype, int natoms) {
  jshandle *js;

  js = (jshandle *)malloc(sizeof(jshandle));
  memset(js, 0, sizeof(jshandle));
  if (fio_open(path, FIO_WRITE, &js->fd) < 0) {
    printf("jsplugin) Could not open file %s for writing\n", path);
    free(js);
    return NULL;
  }

  js->natoms = natoms;
  js->with_unitcell = 1;

  /* emit header information */
  fio_write_str(js->fd, JSHEADERSTRING);
  fio_write_int32(js->fd, JSMAGICNUMBER);
  fio_write_int32(js->fd, JSENDIANISM);
  fio_write_int32(js->fd, JSMAJORVERSION);
  fio_write_int32(js->fd, JSMINORVERSION);

  /* write number of atoms */
  fio_write_int32(js->fd, natoms);

  /* write number of frames, to be updated later */
  js->nframes = 0;
  fio_write_int32(js->fd, js->nframes);

  return js;
}


#if JSMAJORVERSION > 1

static int write_js_structure(void *mydata, int optflags,
                              const molfile_atom_t *atoms) {
  jshandle *js = (jshandle *) mydata;
  int i;

  js->optflags |= JSOPT_STRUCTURE;

  if (optflags & MOLFILE_OCCUPANCY)
    js->optflags |= JSOPT_OCCUPANCY;

  if (optflags & MOLFILE_BFACTOR)
    js->optflags |= JSOPT_BFACTOR;

  if (optflags & MOLFILE_BFACTOR)
    js->optflags |= JSOPT_BFACTOR;

  if (optflags & MOLFILE_MASS)
    js->optflags |= JSOPT_MASS;

  if (optflags & MOLFILE_CHARGE)
    js->optflags |= JSOPT_CHARGE;
 
  if (optflags & MOLFILE_RADIUS)
    js->optflags |= JSOPT_RADIUS;

  if (optflags & MOLFILE_ATOMICNUMBER)
    js->optflags |= JSOPT_ATOMICNUMBER;

  /* write flags data to the file */
  fio_write_int32(js->fd, js->optflags); 
printf("jsplugin) writing option flags: %0x08x\n", js->optflags);

printf("jsplugin) writing structure...\n");
  /* determine whether or not this file contains structure info or not */
  if (js->optflags & JSOPT_STRUCTURE) {
    int numatomnames, numatomtypes, numresnames, numsegids, numchains;
    char **atomnames = NULL;
    char **atomtypes = NULL;
    char **resnames = NULL;
    char **segids = NULL;
    char **chains = NULL;
    short *shortbuf = NULL; /* temp buf for encoding atom records */
    int *intbuf = NULL;     /* temp buf for encoding atom records */
    float *fltbuf = NULL;   /* temp buf for encoding atom records */

    hash_t tmphash;         /* temporary hash table */
    hash_t atomnamehash;
    hash_t atomtypehash;
    hash_t resnamehash;
    hash_t segidhash;
    hash_t chainhash;
    int hashcnt;


printf("jsplugin) counting atom names, types, etc...\n");
    /* generate hash tables to count the number of unique strings */
    hash_init(&tmphash, 127);
    for (i=0; i<js->natoms; i++)
      hash_insert(&tmphash, atoms[i].name, 0);
    numatomnames = tmphash.entries; /* XXX need a query API for this... */
    hash_destroy(&tmphash);

    hash_init(&tmphash, 127);
    for (i=0; i<js->natoms; i++)
      hash_insert(&tmphash, atoms[i].type, 0);
    numatomtypes = tmphash.entries; /* XXX need a query API for this... */
    hash_destroy(&tmphash);

    hash_init(&tmphash, 127);
    for (i=0; i<js->natoms; i++)
      hash_insert(&tmphash, atoms[i].resname, 0);
    numresnames = tmphash.entries; /* XXX need a query API for this... */
    hash_destroy(&tmphash);

    hash_init(&tmphash, 127);
    for (i=0; i<js->natoms; i++)
      hash_insert(&tmphash, atoms[i].segid, 0);
    numsegids = tmphash.entries; /* XXX need a query API for this... */
    hash_destroy(&tmphash);

    hash_init(&tmphash, 127);
    for (i=0; i<js->natoms; i++)
      hash_insert(&tmphash, atoms[i].chain, 0);
    numchains = tmphash.entries; /* XXX need a query API for this... */
    hash_destroy(&tmphash);
 
printf("jsplugin) writing unique string counts...\n");
printf("jsplugin) %d %d %d %d %d\n",
       numatomnames, numatomtypes, numresnames, numsegids, numchains);

    /* write block of name string table sizes */
    fio_write_int32(js->fd, numatomnames); 
    fio_write_int32(js->fd, numatomtypes); 
    fio_write_int32(js->fd, numresnames);
    fio_write_int32(js->fd, numsegids);
    fio_write_int32(js->fd, numchains); 

printf("jsplugin) writing string tables...\n");

    atomnames = (char **) malloc(numatomnames * sizeof(char *));
    atomtypes = (char **) malloc(numatomtypes * sizeof(char *));
    resnames = (char **) malloc(numresnames * sizeof(char *));
    segids = (char **) malloc(numsegids * sizeof(char *));
    chains = (char **) malloc(numchains * sizeof(char *));

printf("jsplugin)   atom names...\n");
    /* generate and write out the string tables */
    hash_init(&atomnamehash, 127);
    for (hashcnt=0,i=0; i<js->natoms; i++) {
      /* add a new string table entry for hash inserts that don't yet exist */
      if (hash_insert(&atomnamehash, atoms[i].name, hashcnt) == HASH_FAIL) {
        atomnames[hashcnt] = (char *) malloc(16 * sizeof(char));
        strcpy(atomnames[hashcnt], atoms[i].name);
        hashcnt++;
      }
    }
    for (i=0; i<numatomnames; i++) {
      fio_fwrite(atomnames[i], 16 * sizeof(char), 1, js->fd);
    }


printf("jsplugin)   atom types...\n");
    hash_init(&atomtypehash, 127);
    for (hashcnt=0,i=0; i<js->natoms; i++) {
      /* add a new string table entry for hash inserts that don't yet exist */
      if (hash_insert(&atomtypehash, atoms[i].type, hashcnt) == HASH_FAIL) {
        atomtypes[hashcnt] = (char *) malloc(16 * sizeof(char));
        strcpy(atomtypes[hashcnt], atoms[i].type);
        hashcnt++;
      }
    }
    for (i=0; i<numatomtypes; i++) {
      fio_fwrite(atomtypes[i], 16 * sizeof(char), 1, js->fd);
    }


printf("jsplugin)   residue names...\n");
    hash_init(&resnamehash, 127);
    for (hashcnt=0,i=0; i<js->natoms; i++) {
      /* add a new string table entry for hash inserts that don't yet exist */
      if (hash_insert(&resnamehash, atoms[i].resname, hashcnt) == HASH_FAIL) {
        resnames[hashcnt] = (char *) malloc(8 * sizeof(char));
        strcpy(resnames[hashcnt], atoms[i].resname);
        hashcnt++;
      }
    }
    for (i=0; i<numresnames; i++) {
      fio_fwrite(resnames[i], 8 * sizeof(char), 1, js->fd);
    }


printf("jsplugin)   segment names...\n");
    hash_init(&segidhash, 127);
    for (hashcnt=0,i=0; i<js->natoms; i++) {
      /* add a new string table entry for hash inserts that don't yet exist */
      if (hash_insert(&segidhash, atoms[i].segid, hashcnt) == HASH_FAIL) {
        segids[hashcnt] = (char *) malloc(8 * sizeof(char));
        strcpy(segids[hashcnt], atoms[i].segid);
        hashcnt++;
      }
    }
    for (i=0; i<numsegids; i++) {
      fio_fwrite(segids[i], 8 * sizeof(char), 1, js->fd);
    }


printf("jsplugin)   chain names...\n");
    hash_init(&chainhash, 127);
    for (hashcnt=0,i=0; i<js->natoms; i++) {
      /* add a new string table entry for hash inserts that don't yet exist */
      if (hash_insert(&chainhash, atoms[i].chain, hashcnt) == HASH_FAIL) {
        chains[hashcnt] = (char *) malloc(2 * sizeof(char));
        strcpy(chains[hashcnt], atoms[i].chain);
        hashcnt++;
      }
    }
    for (i=0; i<numchains; i++) {
      fio_fwrite(chains[i], 2 * sizeof(char), 1, js->fd);
    }


printf("jsplugin) writing numeric field tables...\n");
    /* write out all of the atom fields */
    shortbuf = (void *) malloc(js->natoms * sizeof(short));

    /* write out atom names */
    for (i=0; i<js->natoms; i++) {
      shortbuf[i] = hash_lookup(&atomnamehash, atoms[i].name);
    }    
    fio_fwrite(shortbuf, js->natoms * sizeof(short), 1, js->fd);

    /* write out atom types */
    for (i=0; i<js->natoms; i++) {
      shortbuf[i] = hash_lookup(&atomtypehash, atoms[i].type);
    }    
    fio_fwrite(shortbuf, js->natoms * sizeof(short), 1, js->fd);

    /* write out resnames */
    for (i=0; i<js->natoms; i++) {
      shortbuf[i] = hash_lookup(&resnamehash, atoms[i].resname);
    }    
    fio_fwrite(shortbuf, js->natoms * sizeof(short), 1, js->fd);
    
    /* write out segids */
    for (i=0; i<js->natoms; i++) {
      shortbuf[i] = hash_lookup(&segidhash, atoms[i].segid);
    }    
    fio_fwrite(shortbuf, js->natoms * sizeof(short), 1, js->fd);

    /* write out chains */
    for (i=0; i<js->natoms; i++) {
      shortbuf[i] = hash_lookup(&chainhash, atoms[i].chain);
    }    
    fio_fwrite(shortbuf, js->natoms * sizeof(short), 1, js->fd);

    if (shortbuf != NULL) {
      free(shortbuf);
      shortbuf=NULL;
    }

    /* done with hash tables */
    hash_destroy(&atomnamehash);
    hash_destroy(&atomtypehash);
    hash_destroy(&resnamehash);
    hash_destroy(&segidhash);
    hash_destroy(&chainhash);


    /* 
     * write out integer data blocks 
     */
    intbuf = (int *) malloc(js->natoms * sizeof(int));

printf("jsplugin)   residue indices...\n");
    /* write out resid */
    for (i=0; i<js->natoms; i++) {
      intbuf[i] = atoms[i].resid;
    }    
    fio_fwrite(intbuf, js->natoms * sizeof(int), 1, js->fd);
     
    if (intbuf != NULL) {
      free(intbuf);
      intbuf = NULL;
    }

printf("jsplugin) writing optional per-atom tables...\n");
    /*
     * write out optional single-precision float data blocks
     */ 
    if (js->optflags & (JSOPT_OCCUPANCY | JSOPT_BFACTOR | 
        JSOPT_MASS | JSOPT_RADIUS | JSOPT_CHARGE)) 
      fltbuf = (void *) malloc(js->natoms * sizeof(float));

    /* write out optional data if it exists */

    if (js->optflags & JSOPT_OCCUPANCY) {
printf("jsplugin)   writing occupancy...\n");
      for (i=0; i<js->natoms; i++) {
        fltbuf[i] = atoms[i].occupancy;
      }    
      fio_fwrite(fltbuf, js->natoms * sizeof(float), 1, js->fd);
    }

    if (js->optflags & JSOPT_BFACTOR) {
printf("jsplugin)   writing bfactor...\n");
      for (i=0; i<js->natoms; i++) {
        fltbuf[i] = atoms[i].bfactor;
      }    
      fio_fwrite(fltbuf, js->natoms * sizeof(float), 1, js->fd);
    }

    if (js->optflags & JSOPT_MASS) { 
printf("jsplugin)   writing mass...\n");
      for (i=0; i<js->natoms; i++) {
        fltbuf[i] = atoms[i].mass;
      }    
      fio_fwrite(fltbuf, js->natoms * sizeof(float), 1, js->fd);
    }

    if (js->optflags & JSOPT_CHARGE) { 
printf("jsplugin)   writing charge...\n");
      for (i=0; i<js->natoms; i++) {
        fltbuf[i] = atoms[i].charge;
      }    
      fio_fwrite(fltbuf, js->natoms * sizeof(float), 1, js->fd);
    }

    if (js->optflags & JSOPT_RADIUS) { 
printf("jsplugin)   writing radius...\n");
      for (i=0; i<js->natoms; i++) {
        fltbuf[i] = atoms[i].radius;
      }    
      fio_fwrite(fltbuf, js->natoms * sizeof(float), 1, js->fd);
    }

    if (fltbuf != NULL) {
      free(fltbuf);
      fltbuf=NULL;
    }


    /*
     * write out optional integer data blocks
     */ 
    if (js->optflags & JSOPT_ATOMICNUMBER)
      intbuf = (void *) malloc(js->natoms * sizeof(int));

    if (js->optflags & JSOPT_ATOMICNUMBER) { 
printf("jsplugin)   writing atomic number...\n");
      for (i=0; i<js->natoms; i++) {
        intbuf[i] = atoms[i].atomicnumber;
      }    
      fio_fwrite(intbuf, js->natoms * sizeof(int), 1, js->fd);
    }

    if (intbuf != NULL) {
      free(intbuf);
      intbuf = NULL;
    }


    /*
     * write out bonds and fractional bond orders
     */ 
    if (js->optflags & JSOPT_BONDS) {
printf("jsplugin) writing bonds...\n");
      fio_fwrite(&js->nbonds, sizeof(int), 1, js->fd);
      fio_fwrite(js->bondfrom, js->nbonds * sizeof(int), 1, js->fd);
      fio_fwrite(js->bondto, js->nbonds * sizeof(int), 1, js->fd);

      if (js->optflags & JSOPT_BONDORDERS) {
printf("jsplugin) writing bond orders...\n");
        fio_fwrite(js->bondorders, js->nbonds * sizeof(float), 1, js->fd);
      }
    }

    /*
     * write out angles/dihedrals/impropers/cross-terms
     */
    if (js->optflags & JSOPT_ANGLES) {
printf("jsplugin) writing angles/dihedrals/impropers...\n");
      fio_fwrite(&js->numangles, sizeof(int), 1, js->fd);
      fio_fwrite(js->angles, sizeof(int)*3*js->numangles, 1, js->fd);

      fio_fwrite(&js->numdihedrals, sizeof(int), 1, js->fd);
      fio_fwrite(js->dihedrals, sizeof(int)*4*js->numdihedrals, 1, js->fd);

      fio_fwrite(&js->numimpropers, sizeof(int), 1, js->fd);
      fio_fwrite(js->impropers, sizeof(int)*4*js->numimpropers, 1, js->fd);
    }
    if (js->optflags & JSOPT_CTERMS) {
printf("jsplugin) writing cross-terms\n");
      fio_fwrite(&js->numcterms, sizeof(int), 1, js->fd);
      fio_fwrite(js->cterms, sizeof(int)*8*js->numcterms, 1, js->fd);
    }

    return MOLFILE_SUCCESS;
  }

  /* else, we have no structure information */
  return MOLFILE_NOSTRUCTUREDATA;
}


static int write_js_bonds(void *mydata, int nbonds, int *fromptr, int *toptr, 
                          float *bondorder,  int *bondtype, 
                          int nbondtypes, char **bondtypename) {
  jshandle *js = (jshandle *) mydata;

  if (nbonds > 0 && fromptr != NULL && toptr != NULL) {
    js->optflags |= JSOPT_BONDS; 

    /* save bond info until we actually write out the structure file */
    js->nbonds = nbonds;
    js->bondfrom = (int *) malloc(nbonds * sizeof(int));
    memcpy(js->bondfrom, fromptr, nbonds * sizeof(int));
    js->bondto = (int *) malloc(nbonds * sizeof(int));
    memcpy(js->bondto, toptr, nbonds * sizeof(int));

    if (bondorder != NULL) {
      js->optflags |= JSOPT_BONDORDERS;
      js->bondorders = (float *) malloc(nbonds * sizeof(float));
      memcpy(js->bondorders, bondorder, nbonds * sizeof(float));
    }
  }

  return MOLFILE_SUCCESS;
}

#if vmdplugin_ABIVERSION > 14
static int write_js_angles(void * v, int numangles, const int *angles,
                           const int *angletypes, int numangletypes,
                           const char **angletypenames, int numdihedrals, 
                           const int *dihedrals, const int *dihedraltype,
                           int numdihedraltypes, const char **dihedraltypenames,
                           int numimpropers, const int *impropers, 
                           const int *impropertypes, int numimpropertypes, 
                           const char **impropertypenames, int numcterms, 
                           const int *cterms, int ctermcols, int ctermrows) {
  jshandle *js = (jshandle *) v;

  /* save info until we actually write out the structure file */
  js->numangles = numangles;
  js->numdihedrals = numdihedrals;
  js->numimpropers = numimpropers;
  js->numcterms = numcterms;

  if (js->numangles > 0 || js->numdihedrals > 0 || js->numimpropers > 0) {
    js->optflags |= JSOPT_ANGLES;

    js->angles = (int *) malloc(3*js->numangles*sizeof(int));
    memcpy(js->angles, angles, 3*js->numangles*sizeof(int));
    js->dihedrals = (int *) malloc(4*js->numdihedrals*sizeof(int));
    memcpy(js->dihedrals, dihedrals, 4*js->numdihedrals*sizeof(int));
    js->impropers = (int *) malloc(4*js->numimpropers*sizeof(int));
    memcpy(js->impropers, impropers, 4*js->numimpropers*sizeof(int));
  }
  if (js->numcterms > 0) {
    js->optflags |= JSOPT_CTERMS;

    js->cterms = (int *) malloc(8*js->numcterms*sizeof(int));
    memcpy(js->cterms, cterms, 8*js->numcterms*sizeof(int));
  }

  return MOLFILE_SUCCESS;
}
#else
static int write_js_angles(void * v,
        int numangles,    const int *angles,    const double *angleforces,
        int numdihedrals, const int *dihedrals, const double *dihedralforces,
        int numimpropers, const int *impropers, const double *improperforces,
        int numcterms,   const int *cterms,
        int ctermcols, int ctermrows, const double *ctermforces) {
  jshandle *js = (jshandle *) v;

  /* save info until we actually write out the structure file */
  js->numangles = numangles;
  js->numdihedrals = numdihedrals;
  js->numimpropers = numimpropers;
  js->numcterms = numcterms;

  if (js->numangles > 0 || js->numdihedrals > 0 || js->numimpropers > 0) {
    js->optflags |= JSOPT_ANGLES;

    js->angles = (int *) malloc(3*js->numangles*sizeof(int));
    memcpy(js->angles, angles, 3*js->numangles*sizeof(int));
    js->dihedrals = (int *) malloc(4*js->numdihedrals*sizeof(int));
    memcpy(js->dihedrals, dihedrals, 4*js->numdihedrals*sizeof(int));
    js->impropers = (int *) malloc(4*js->numimpropers*sizeof(int));
    memcpy(js->impropers, impropers, 4*js->numimpropers*sizeof(int));
  }
  if (js->numcterms > 0) {
    js->optflags |= JSOPT_CTERMS;

    js->cterms = (int *) malloc(8*js->numcterms*sizeof(int));
    memcpy(js->cterms, cterms, 8*js->numcterms*sizeof(int));
  }

  return MOLFILE_SUCCESS;
}
#endif
#endif


static int write_js_timestep(void *v, const molfile_timestep_t *ts) { 
  jshandle *js = (jshandle *)v;
  double unitcell[6];

  js->nframes++; /* increment frame count written to the file so far */

  unitcell[0] = ts->A;
  unitcell[1] = ts->B;
  unitcell[2] = ts->C;
  unitcell[3] = sin((M_PI_2 / 90.0) * (90.0 - ts->alpha));
  unitcell[4] = sin((M_PI_2 / 90.0) * (90.0 - ts->beta));
  unitcell[5] = sin((M_PI_2 / 90.0) * (90.0 - ts->gamma));

  /* coordinates for all atoms */
  fio_fwrite(ts->coords, js->natoms * 3 * sizeof(float), 1, js->fd);

  /* PBC unit cell info */ 
  fio_fwrite(&unitcell[0], 6 * sizeof(double), 1, js->fd);

  return MOLFILE_SUCCESS;
}

static void close_js_write(void *v) {
  jshandle *js = (jshandle *)v;

  /* update the trajectory header information */
  fio_fseek(js->fd, JSNFRAMESOFFSET, FIO_SEEK_SET);
  fio_write_int32(js->fd, js->nframes);
  fio_fseek(js->fd, 0, FIO_SEEK_END);

  fio_fclose(js->fd);

#if JSMAJORVERSION > 1
  if (js->bondfrom)
    free(js->bondfrom);
  if (js->bondto)
    free(js->bondto);
  if (js->bondorders)
    free(js->bondorders);

  if (js->angles)
    free(js->angles);
  if (js->dihedrals)
    free(js->dihedrals);
  if (js->impropers)
    free(js->impropers);
  if (js->cterms)
    free(js->cterms);
#endif

  free(js);
}


/*
 * Initialization stuff here
 */
static molfile_plugin_t plugin;

VMDPLUGIN_API int VMDPLUGIN_init() {
  memset(&plugin, 0, sizeof(molfile_plugin_t));
  plugin.abiversion = vmdplugin_ABIVERSION;
  plugin.type = MOLFILE_PLUGIN_TYPE;
  plugin.name = "js";
  plugin.prettyname = "js";
  plugin.author = "John Stone";
  plugin.majorv = JSMAJORVERSION;
  plugin.minorv = JSMINORVERSION;
  plugin.is_reentrant = VMDPLUGIN_THREADSAFE;
  plugin.filename_extension = "js";
  plugin.open_file_read = open_js_read;
#if JSMAJORVERSION > 1
  plugin.read_structure = read_js_structure;
  plugin.read_bonds = read_js_bonds;
  plugin.read_angles = read_js_angles;
#endif
  plugin.read_next_timestep = read_js_timestep;
  plugin.close_file_read = close_js_read;
  plugin.open_file_write = open_js_write;
#if JSMAJORVERSION > 1
  plugin.write_structure = write_js_structure;
  plugin.write_bonds = write_js_bonds;
  plugin.write_angles = write_js_angles;
#endif
  plugin.write_timestep = write_js_timestep;
  plugin.close_file_write = close_js_write;
  return VMDPLUGIN_SUCCESS;
}

VMDPLUGIN_API int VMDPLUGIN_register(void *v, vmdplugin_register_cb cb) {
  (*cb)(v, (vmdplugin_t *)&plugin);
  return VMDPLUGIN_SUCCESS;
}

VMDPLUGIN_API int VMDPLUGIN_fini() {
  return VMDPLUGIN_SUCCESS;
}

  
#ifdef TEST_JSPLUGIN

#include <sys/time.h>

/* get the time of day from the system clock, and store it (in seconds) */
double time_of_day(void) {
#if defined(_MSC_VER)
  double t;

  t = GetTickCount();
  t = t / 1000.0;

  return t;
#else
  struct timeval tm;
  struct timezone tz;

  gettimeofday(&tm, &tz);
  return((double)(tm.tv_sec) + (double)(tm.tv_usec)/1000000.0);
#endif
}

int main(int argc, char *argv[]) {
  molfile_timestep_t timestep;
  void *v;
  jshandle *js;
  int i, natoms;
  float sizeMB =0.0, totalMB = 0.0;
  double starttime, endtime, totaltime = 0.0;

  while (--argc) {
    ++argv; 
    natoms = 0;
    v = open_js_read(*argv, "js", &natoms);
    if (!v) {
      printf("jsplugin) open_js_read failed for file %s\n", *argv);
      return 1;
    }
    js = (jshandle *)v;
    sizeMB = ((natoms * 3.0) * js->nframes * 4.0) / (1024.0 * 1024.0);
    totalMB += sizeMB; 
    printf("jsplugin) file: %s\n", *argv);
    printf("jsplugin)   %d atoms, %d frames, size: %6.1fMB\n", natoms, js->nframes, sizeMB);

    starttime = time_of_day();
    timestep.coords = (float *)malloc(3*sizeof(float)*natoms);
    for (i=0; i<js->nframes; i++) {
      int rc = read_js_timestep(v, natoms, &timestep);
      if (rc) {
        printf("jsplugin) error in read_js_timestep on frame %d\n", i);
        return 1;
      }
    }
    endtime = time_of_day();
    close_js_read(v);
    totaltime += endtime - starttime;
    printf("jsplugin)  Time: %5.1f seconds\n", endtime - starttime);
    printf("jsplugin)  Speed: %5.1f MB/sec, %5.1f timesteps/sec\n", sizeMB / (endtime - starttime), (js->nframes / (endtime - starttime)));
  }
  printf("jsplugin) Overall Size: %6.1f MB\n", totalMB);
  printf("jsplugin) Overall Time: %6.1f seconds\n", totaltime);
  printf("jsplugin) Overall Speed: %5.1f MB/sec\n", totalMB / totaltime);
  return 0;
}
      
#endif

