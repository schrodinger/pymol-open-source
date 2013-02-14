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
 *      $Revision: 1.56 $       $Date: 2011/07/19 14:55:46 $
 *
 ***************************************************************************
 * DESCRIPTION:
 *   This plugin implements a high-performance binary molecular structure
 *   and trajectory storage format.  This file format currently uses a simple 
 *   non-redundant hash table approach for compression of per-atom
 *   character strings, properties, and tags.  Trajectory data is stored 
 *   with a file structure that avoids the need to transpose or convert
 *   dense blocks of cartesian coordinates into the most commonly used
 *   interleaved  x/y/z ordering.  The file structure also enables zero-copy 
 *   vectorized I/O methods to be used high-performance reads for 
 *   visualization and analysis with reduced operating system overhead.
 *   The plugin optionally supports the use of a block-based file structure
 *   and block-aligned memory buffers for direct I/O that bypasses the 
 *   OS filesystem buffer caches for multi-gigabyte-per-second read rates
 *   from SSD RAID arrays.
 *
 *   At present, only VMD, NAMD, and psfgen make use of this format.
 *   It started out as a test/example code and is slowly becoming
 *   more robust and featureful.
 *
 *   We should be able to implement a selective read approach that gathers
 *   discontiguous blocks of the trajectory using the POSIX lio_listio()
 *   APIs.  On Unix we currently use I/O wrappers that are based on the 
 *   lseek() and readv() APIs.  By using lio_listio() we could eliminate
 *   the separate lseek calls and service multiple timestep reads in a 
 *   single request, even included cases with discontiguous requests.
 *
 *   VMD test results for Linux host with an 8-way RAID0 of commodity 
 *   Intel 510 SSDs with SATA III 6Gbit/sec interfaces:
 *     Non-direct I/O using standard readv(): 1203 MB/sec
 *     Direct I/O, readv(), 4KB blocked file: 2130 MB/sec
 *
 *  Standalone test binary compilation flags for 64-bit Linux:
 *  cc -O3 -m64 -I../../include -DTEST_JSPLUGIN jsplugin.c \
 *    -o ~/bin/readjs -lm
 *  
 *  Standalone test binary compilation flags for 64-bit Linux w/ CUDA:
 *  cc -O3 -m64 -I../../include -I/usr/local/encap/cuda-4.0/include \
 *    -DTEST_JSPLUGIN -DENABLECUDATESTS jsplugin.c \
 *    -o ~/bin/readjs -L/usr/local/encap/cuda-4.0/lib64 -lcudart -lm
 *
 *  Standalone test binary compilation flags for Solaris:
 *  cc -fast -xarch=v9a -I../../include -DTEST_JSPLUGIN jsplugin.c \
 *    -o ~/bin/readjs -lm
 *
 *  Profiling flags for Solaris:
 *  cc -xpg -fast -xarch=v9a -g -I../../include -DTEST_JSPLUGIN jsplugin.c \
 *    -o ~/bin/readjs -lm
 *
 ***************************************************************************/

#if 1
#define ENABLEJSSHORTREADS 1
#endif

#define VMDPLUGIN_STATIC
#include "largefiles.h"   /* platform dependent 64-bit file I/O defines */
#include "fastio.h"       /* must come before others, for O_DIRECT...   */

#include <sys/stat.h>
#include <sys/types.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "hash.h"
#include "endianswap.h"
#include "molfile_plugin.h"


/* allocate memory and return a pointer that is aligned on a given   */
/* byte boundary, to be used for page- or sector-aligned I/O buffers */
/* We use this if posix_memalign() is not available...               */
#if 1 /* sizeof(unsigned long) == sizeof(void*) */
#define myintptrtype unsigned long
#elif 1   /* sizeof(size_t) == sizeof(void*) */
#define myintptrtype size_t
#else
#define myintptrtype uintptr_t  /* C99 */
#endif
static void *alloc_aligned_ptr(size_t sz, size_t blocksz, void **unalignedptr) {
  /* pad the allocation to an even multiple of the block size */
  size_t padsz = (sz + (blocksz - 1)) & (~(blocksz - 1));
  void * ptr = malloc(padsz + blocksz);
  *unalignedptr = ptr;
  return (void *) ((((myintptrtype) ptr) + (blocksz-1)) & (~(blocksz-1)));
}


#ifndef M_PI_2
#define M_PI_2 1.57079632679489661922
#endif

#define JSHEADERSTRING   "JS Binary Structure and Trajectory File Format"                
#define JSMAGICNUMBER    0x31337
#define JSENDIANISM      0x12345678

#define JSMAJORVERSION   2
#define JSMINORVERSION   9

#define JSNFRAMESOFFSET  (strlen(JSHEADERSTRING) + 20)

#define JSNOERR             0
#define JSBADFILE           1
#define JSBADFORMAT         2


#define JSOPT_NOOPTIONS     0x00000000  /* no structure, only coords    */

/* Timesteps are block-size padded and page- or sector-aligned for      */
/* direct I/O, using  OS-specific APIs that completely bypass the OS    */
/* kernel filesystem buffer cache.                                      */
/* The use of direct I/O APIs can raise performance tremendously on     */
/* high-end RAIDs.  Tests on an 8-way RAID0 of Intel 510 SSDs raise the */
/* peak I/O rate from 1100 MB/sec up to 2020 MB/sec with direct I/O.    */
#define JSOPT_TS_BLOCKIO    0x10000000

/* large data blocks */
#define JSOPT_STRUCTURE     0x00000001  /* file contains structure data */
#define JSOPT_BONDS         0x00000002  /* file contains bond info      */
#define JSOPT_BONDORDERS    0x00000004  /* file contains bond orders    */
#define JSOPT_ANGLES        0x00000008  /* file contains angle info     */
#define JSOPT_CTERMS        0x00000010  /* file contains cross-terms    */

/* optional per-atom fields */
#define JSOPT_OCCUPANCY     0x00000100  /* file contains occupancy      */
#define JSOPT_BFACTOR       0x00000200  /* file contains b-factor       */
#define JSOPT_MASS          0x00000400  /* file contains masses         */
#define JSOPT_CHARGE        0x00000800  /* file contains charges        */
#define JSOPT_RADIUS        0x00001000  /* file contains radii          */
#define JSOPT_ATOMICNUMBER  0x00002000  /* file contains atomic numbers */

typedef struct {
  fio_fd fd;
  int natoms;

#if JSMAJORVERSION > 1
  int parsed_structure;        /* flag indicating structure is parsed   */
  char *path;                  /* path to file                          */

  /* info for block-based direct I/O */ 
  int directio_enabled;        /* block-based direct I/O is available   */
  fio_fd directio_fd;          /* block-based direct I/O using O_DIRECT */
  int directio_block_size;     /* block size to use for direct ts I/O   */
  void *directio_ucell_ptr;    /* unaligned unit cell buffer ptr        */
  void *directio_ucell_blkbuf; /* block-aligned unit cell buffer pt r   */

  /* timestep file offset, block padding, and stride information */
  fio_size_t ts_file_offset;   /* file offset to first timestep         */
  fio_size_t ts_crd_sz;        /* size of TS coordinates                */
  fio_size_t ts_crd_padsz;     /* size of TS block-padded coordinates   */
  fio_size_t ts_ucell_sz;      /* size of TS unit cell                  */
  fio_size_t ts_ucell_padsz;   /* size of TS block-padded unit cell     */
  
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
    perror("jsplugin) stat: ");
/*    return NULL; */
  }

  js = (jshandle *)malloc(sizeof(jshandle));
  memset(js, 0, sizeof(jshandle));
#if JSMAJORVERSION > 1
  js->parsed_structure=0;
  js->directio_block_size=1;
  js->directio_ucell_ptr = NULL;
  js->directio_ucell_blkbuf = NULL;

  js->directio_enabled=0;
  js->ts_file_offset=0;
  js->ts_crd_sz=0;
  js->ts_ucell_sz=0;
  js->ts_crd_padsz=0;
  js->ts_ucell_padsz=0;
#endif

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
    fio_fclose(js->fd);
    free(js);
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

  /* copy path if we succeeded in opening the file */
  js->path = (char *) calloc(strlen(path)+1, 1);
  strcpy(js->path, path);

  return js;
}


#if JSMAJORVERSION > 1

/* Compute the file offset for the first timestep and move */
/* the file pointer to the correct position to read/write  */
/* the first timestep.  Takes care of block alignment when */
/* needed.                                                 */ 
static int js_calc_timestep_blocking_info(void *mydata) {
  fio_size_t ts_block_offset;
  fio_size_t bszmask;
  jshandle *js = (jshandle *) mydata;
  int iorc=0;

  /* Record the current file offset so we can use it to */
  /* compute the absolute offset to the first timestep. */
  js->ts_file_offset = fio_ftell(js->fd);

  /* pad current offset to the start of the next block  */ 
  bszmask = js->directio_block_size - 1;
  ts_block_offset = (js->ts_file_offset + bszmask) & (~bszmask);

  printf("jsplugin) TS block size %d  curpos: %d  blockpos: %d\n", 
         (int) js->directio_block_size, 
         (int) js->ts_file_offset, 
         (int) ts_block_offset);

  /* seek to the first block of the first timestep */
  js->ts_file_offset = ts_block_offset;
  if (js->directio_enabled)
    iorc = fio_fseek(js->directio_fd, js->ts_file_offset, FIO_SEEK_SET);
  else
    iorc = fio_fseek(js->fd, js->ts_file_offset, FIO_SEEK_SET);
  if (iorc < 0) {
    perror("jsplugin) fseek(): ");
  }

  /* compute timestep block padding/skipping for both */
  /* coordinate blocks and unit cell blocks           */
  js->ts_crd_sz = js->natoms * 3 * sizeof(float);
  js->ts_crd_padsz = (js->ts_crd_sz + bszmask) & (~bszmask);

  js->ts_ucell_sz = 6 * sizeof(double);
  js->ts_ucell_padsz = (js->ts_ucell_sz + bszmask) & (~bszmask);

  /* allocate TS unit cell buffer in an aligned, block-size-multiple buffer */
  /* unaligned unit cell buffer ptr */
#if defined(USE_POSIX_MEMALIGN)
  if (posix_memalign((void**) &js->directio_ucell_ptr, 
      js->directio_block_size, js->ts_ucell_padsz)) {
    printf("jsplugin) Couldn't allocate aligned unit cell block buffer!\n");
  }
  /* the returned pointer is already block-aligned, and can free() */
  js->directio_ucell_blkbuf = js->directio_ucell_ptr;
#else
  js->directio_ucell_blkbuf = (float *) 
    alloc_aligned_ptr(js->ts_ucell_padsz, js->directio_block_size, 
                      (void**) &js->directio_ucell_ptr);
#endif

  printf("jsplugin) TS crds sz: %ld psz: %ld  ucell sz: %ld psz: %ld\n",
         (long) js->ts_crd_sz,
         (long) js->ts_crd_padsz, 
         (long) js->ts_ucell_sz, 
         (long) js->ts_ucell_padsz);

  return MOLFILE_SUCCESS;
}


static int read_js_structure(void *mydata, int *optflags,
                             molfile_atom_t *atoms) {
  jshandle *js = (jshandle *) mydata;
  int i;

  if (optflags != NULL)
    *optflags = MOLFILE_NOOPTIONS; /* set to no options until we read them */

  /* read flags data from the file */
  fio_read_int32(js->fd, &js->optflags); 
  if (js->reverseendian)
    swap4_aligned(&js->optflags, 1);
printf("jsplugin) read option flags: %0x08x\n", js->optflags);

  /* Check to see if block-based trajectory I/O is used  */
  /* and read in the block size for this file.           */
  if (js->optflags & JSOPT_TS_BLOCKIO) {
    fio_fread(&js->directio_block_size, sizeof(int), 1, js->fd);
    if (js->reverseendian)
      swap4_aligned(&js->directio_block_size, 1);

    printf("jsplugin) Block-based I/O enabled: block size %d bytes\n", 
           js->directio_block_size);

    if (fio_open(js->path, FIO_READ | FIO_DIRECT, &js->directio_fd) < 0) {
      printf("jsplugin) Direct I/O unavailable for file '%s'\n", js->path);
    } else {
      js->directio_enabled = 1;
      printf("jsplugin) Direct I/O enabled for file '%s'\n", js->path);
    } 
  }

#if defined(ENABLEJSSHORTREADS)
  /* test code for an implementation that does short reads that */
  /* skip bulk solvent, useful for faster loading of very large */
  /* structures                                                 */
  if (getenv("VMDJSMAXATOMIDX") != NULL) {
    fio_size_t bszmask;

    int maxatomidx = atoi(getenv("VMDJSMAXATOMIDX"));
    if (maxatomidx < 0)
      maxatomidx = 0;
    if (maxatomidx >= js->natoms)
      maxatomidx = js->natoms - 1;

    printf("jsplugin) Short-reads of timesteps enabled: %d / %d atoms (%.2f%%)\n",
           maxatomidx, js->natoms, 100.0*(maxatomidx+1) / ((float) js->natoms));
  }
#endif

  /* Mark the handle to indicate we've parsed the structure.             */
  /* If any errors occur after this point, they are likely fatal anyway. */
  js->parsed_structure = 1;

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
      swap4_aligned(&numatomnames, 1);
      swap4_aligned(&numatomtypes, 1);
      swap4_aligned(&numresnames, 1);
      swap4_aligned(&numsegids, 1);
      swap4_aligned(&numchains, 1);
    }

printf("jsplugin) reading string tables...\n");
printf("jsplugin) %d %d %d %d %d\n",
       numatomnames, numatomtypes, numresnames, numsegids, numchains);

    /* skip forward to first TS if the caller gives us NULL ptrs */
    if (optflags == NULL && atoms == NULL) {
      size_t offset=0;
      offset += numatomnames * (16 * sizeof(char));
      offset += numatomtypes * (16 * sizeof(char));
      offset += numresnames  * (8 * sizeof(char));
      offset += numsegids    * (8 * sizeof(char));
      offset += numchains    * (2 * sizeof(char));
      offset += js->natoms * sizeof(short); /* atom name indices    */
      offset += js->natoms * sizeof(short); /* atom type indices    */
      offset += js->natoms * sizeof(short); /* residue name indices */
      offset += js->natoms * sizeof(short); /* segment name indices */
      offset += js->natoms * sizeof(short); /* chain name indices   */
      offset += js->natoms * sizeof(int);   /* residue indices      */
      
      /* optional per-atom fields */
      if (js->optflags & JSOPT_OCCUPANCY)
        offset += js->natoms * sizeof(float); 
      if (js->optflags & JSOPT_BFACTOR)
        offset += js->natoms * sizeof(float); 
      if (js->optflags & JSOPT_MASS)
        offset += js->natoms * sizeof(float); 
      if (js->optflags & JSOPT_CHARGE)
        offset += js->natoms * sizeof(float); 
      if (js->optflags & JSOPT_RADIUS)
        offset += js->natoms * sizeof(float); 
      if (js->optflags & JSOPT_ATOMICNUMBER)
        offset += js->natoms * sizeof(int);

      fio_fseek(js->fd, offset, FIO_SEEK_CUR);
      offset=0;

      /* these require actually seeking as we process... */
      if (js->optflags & JSOPT_BONDS) {
        fio_fread(&js->nbonds, sizeof(int), 1, js->fd);
        if (js->reverseendian)
          swap4_aligned(&js->nbonds, 1);
printf("jsplugin)   %d bonds...\n", js->nbonds);

        offset += 2 * js->nbonds * sizeof(int);
        if (js->optflags & JSOPT_BONDORDERS)
          offset += js->nbonds * sizeof(float);

        fio_fseek(js->fd, offset, FIO_SEEK_CUR);
        offset=0;
      }

      if (js->optflags & JSOPT_ANGLES) {
        fio_fread(&js->numangles, sizeof(int), 1, js->fd);
        if (js->reverseendian)
          swap4_aligned(&js->numangles, 1);
printf("jsplugin)   %d angles...\n", js->numangles);
        fio_fseek(js->fd, sizeof(int)*3*js->numangles, FIO_SEEK_CUR);

        fio_fread(&js->numdihedrals, sizeof(int), 1, js->fd);
        if (js->reverseendian)
          swap4_aligned(&js->numdihedrals, 1);
printf("jsplugin)   %d dihedrals...\n", js->numdihedrals);
        fio_fseek(js->fd, sizeof(int)*4*js->numdihedrals, FIO_SEEK_CUR);

        fio_fread(&js->numimpropers, sizeof(int), 1, js->fd);
        if (js->reverseendian)
          swap4_aligned(&js->numimpropers, 1);
printf("jsplugin)   %d impropers...\n", js->numimpropers);
        fio_fseek(js->fd, sizeof(int)*4*js->numimpropers, FIO_SEEK_CUR);
      }

      if (js->optflags & JSOPT_CTERMS) {
        fio_fread(&js->numcterms, sizeof(int), 1, js->fd);
        if (js->reverseendian)
          swap4_aligned(&js->numcterms, 1);
 printf("jsplugin)   %d cterms...\n", js->numcterms);
        fio_fseek(js->fd, sizeof(int)*8*js->numcterms, FIO_SEEK_CUR);
      }
  
      /* record the file offset for the first timestep */
      js_calc_timestep_blocking_info(js);

      return MOLFILE_SUCCESS;
    }


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
    shortbuf = (short *) malloc(js->natoms * sizeof(short));

printf("jsplugin)   atom name indices...\n");
    /* read in atom names */
    fio_fread(shortbuf, js->natoms * sizeof(short), 1, js->fd);
    if (js->reverseendian)
      swap2_aligned(shortbuf, js->natoms);
    for (i=0; i<js->natoms; i++) {
      strcpy(atoms[i].name, atomnames[shortbuf[i]]);
    }
    for (i=0; i<numatomnames; i++)
      free(atomnames[i]);
    free(atomnames);

printf("jsplugin)   atom type indices...\n");
    /* read in atom types */
    fio_fread(shortbuf, js->natoms * sizeof(short), 1, js->fd);
    if (js->reverseendian)
      swap2_aligned(shortbuf, js->natoms);
    for (i=0; i<js->natoms; i++) {
      strcpy(atoms[i].type, atomtypes[shortbuf[i]]);
    }
    for (i=0; i<numatomtypes; i++)
      free(atomtypes[i]);
    free(atomtypes);

printf("jsplugin)   residue name indices...\n");
    /* read in resnames */
    fio_fread(shortbuf, js->natoms * sizeof(short), 1, js->fd);
    if (js->reverseendian)
      swap2_aligned(shortbuf, js->natoms);
    for (i=0; i<js->natoms; i++) {
      strcpy(atoms[i].resname, resnames[shortbuf[i]]);
    }
    for (i=0; i<numresnames; i++)
      free(resnames[i]);
    free(resnames);
    
printf("jsplugin)   segment name indices...\n");
    /* read in segids */
    fio_fread(shortbuf, js->natoms * sizeof(short), 1, js->fd);
    if (js->reverseendian)
      swap2_aligned(shortbuf, js->natoms);
    for (i=0; i<js->natoms; i++) {
      strcpy(atoms[i].segid, segids[shortbuf[i]]);
    }
    for (i=0; i<numsegids; i++)
      free(segids[i]);
    free(segids);

printf("jsplugin)   chain name indices...\n");
    /* read in chains */
    fio_fread(shortbuf, js->natoms * sizeof(short), 1, js->fd);
    if (js->reverseendian)
      swap2_aligned(shortbuf, js->natoms);
    for (i=0; i<js->natoms; i++) {
      strcpy(atoms[i].chain, chains[shortbuf[i]]);
    }
    for (i=0; i<numchains; i++)
      free(chains[i]);
    free(chains);

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
      fltbuf = (float *) malloc(js->natoms * sizeof(float));

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
      intbuf = (int *) malloc(js->natoms * sizeof(int));

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
        js->bondorders = (float *) malloc(js->nbonds * sizeof(float));
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
        swap4_aligned(js->angles, 3*js->numangles);

      fio_fread(&js->numdihedrals, sizeof(int), 1, js->fd);
      if (js->reverseendian)
        swap4_aligned(&js->numdihedrals, 1);
printf("jsplugin)   %d dihedrals...\n", js->numdihedrals);
      js->dihedrals = (int *) malloc(4 * js->numdihedrals * sizeof(int));
      fio_fread(js->dihedrals, sizeof(int)*4*js->numdihedrals, 1, js->fd);
      if (js->reverseendian)
        swap4_aligned(js->dihedrals, 4*js->numdihedrals);

      fio_fread(&js->numimpropers, sizeof(int), 1, js->fd);
      if (js->reverseendian)
        swap4_aligned(&js->numimpropers, 1);
      js->impropers = (int *) malloc(4 * js->numimpropers * sizeof(int));
printf("jsplugin)   %d impropers...\n", js->numimpropers);
      fio_fread(js->impropers, sizeof(int)*4*js->numimpropers, 1, js->fd);
      if (js->reverseendian)
        swap4_aligned(js->impropers, 4*js->numimpropers);
    }
    if (js->optflags & JSOPT_CTERMS) {
      fio_fread(&js->numcterms, sizeof(int), 1, js->fd);
      if (js->reverseendian)
        swap4_aligned(&js->numcterms, 1);
      js->cterms = (int *) malloc(8 * js->numcterms * sizeof(int));
printf("jsplugin)   %d cterms...\n", js->numcterms);
      fio_fread(js->cterms, sizeof(int)*8*js->numcterms, 1, js->fd);
      if (js->reverseendian)
        swap4_aligned(js->cterms, 8*js->numcterms);
    }

printf("jsplugin) final optflags: %08x\n", *optflags);
printf("jsplugin) structure information complete\n");

    /* record the file offset for the first timestep */
    js_calc_timestep_blocking_info(js);

    return MOLFILE_SUCCESS;
  }

printf("jsplugin) no structure information available\n");

  /* record the file offset for the first timestep */
  js_calc_timestep_blocking_info(js);

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
  fio_size_t framelen;

#if JSMAJORVERSION > 1
  /* If we haven't yet read (or skipped) the structure data, then we    */
  /* need to begin by skipping past it before we try to read the        */
  /* first timestep.  In the case of files with block-aligned timesteps,*/
  /* this will also get our file pointer to the right block-aligned     */
  /* location.                                                          */
  if (!js->parsed_structure)
    read_js_structure(v, NULL, NULL);
#endif

  /* compute total read/seek size of timestep */
  framelen = js->ts_crd_padsz + js->ts_ucell_padsz;

  /* if we have a valid ts pointer, read the timestep, otherwise skip it */ 
  if (ts != NULL) {
    fio_size_t readlen; 
    fio_iovec iov[2];

    /* set unit cell pointer to the TS block-aligned buffer area */
    double *unitcell = (double *) js->directio_ucell_blkbuf;

    unitcell[0] = unitcell[2] = unitcell[5] = 1.0f;
    unitcell[1] = unitcell[3] = unitcell[4] = 90.0f;

#if defined(ENABLEJSSHORTREADS)
    /* test code for an implementation that does short reads that */
    /* skip bulk solvent, useful for faster loading of very large */
    /* structures                                                 */
    if (getenv("VMDJSMAXATOMIDX") != NULL) {
      fio_size_t bszmask;
      int maxatompadsz, skipatompadsz;

      int maxatomidx = atoi(getenv("VMDJSMAXATOMIDX"));
      if (maxatomidx < 0)
        maxatomidx = 0;
      if (maxatomidx >= js->natoms)
        maxatomidx = js->natoms - 1;

      /* pad max read to the start of the next block  */
      bszmask = js->directio_block_size - 1;
      maxatompadsz = ((maxatomidx*3*sizeof(float)) + bszmask) & (~bszmask);
      skipatompadsz = js->ts_crd_padsz - maxatompadsz;

      readlen=0;
      if (js->directio_enabled) {
        if (fio_fread(ts->coords, maxatompadsz, 1, js->directio_fd) == 1)
          readlen = maxatompadsz;
        if (fio_fseek(js->directio_fd, skipatompadsz, FIO_SEEK_CUR) == 0)
          readlen += skipatompadsz;
        if (fio_fread(unitcell, js->ts_ucell_padsz, 1, js->directio_fd) == 1)
          readlen += js->ts_ucell_padsz;
      } else {
        if (fio_fread(ts->coords, maxatompadsz, 1, js->fd) == 1)
          readlen = maxatompadsz;
        if (fio_fseek(js->fd, skipatompadsz, FIO_SEEK_CUR) == 0)
          readlen += skipatompadsz;
        if (fio_fread(unitcell, js->ts_ucell_padsz, 1, js->fd) == 1)
          readlen += js->ts_ucell_padsz;
      }

#if 0
      /* clear all non-read atom coords to zeros */
      memset(ts->coords+3*maxatomidx,0,3*sizeof(float)*(js->natoms-maxatomidx));
#endif

    }  else {
#endif
 
    /* setup the I/O vector */
    iov[0].iov_base = (fio_caddr_t) ts->coords;   /* read coordinates    */
    iov[0].iov_len  = js->ts_crd_padsz;
    iov[1].iov_base = (fio_caddr_t) unitcell;     /* read PBC unit cell  */
    iov[1].iov_len  = js->ts_ucell_padsz;

    /* Do all of the reads with a single syscall, for peak efficiency. */
    /* On smart kernels, readv() causes only one context switch, and   */
    /* can effeciently scatter the reads to the various buffers.       */
    if (js->directio_enabled)
      readlen = fio_readv(js->directio_fd, &iov[0], 2); 
    else
      readlen = fio_readv(js->fd, &iov[0], 2); 

#if defined(ENABLEJSSHORTREADS)
   }
#endif 
 
    /* check the number of read bytes versus what we expected */
    if (readlen != framelen) {
      if (readlen < 0) {
        perror("jsplugin) fio_readv(): ");
      } else {
        printf("jsplugin) mismatched read: %ld, expected %ld\n", 
               (long) readlen, (long) framelen);
      }

      return MOLFILE_EOF;
    }

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
    if (js->directio_enabled) {
      if (fio_fseek(js->directio_fd, framelen, FIO_SEEK_CUR)) 
        return MOLFILE_EOF;
    } else {
      if (fio_fseek(js->fd, framelen, FIO_SEEK_CUR)) 
        return MOLFILE_EOF;
    }
  }
 
  return MOLFILE_SUCCESS;
}
 

static void close_js_read(void *v) {
  jshandle *js = (jshandle *)v;
  fio_fclose(js->fd);

#if JSMAJORVERSION > 1
  if (js->path)
    free(js->path);

  if (js->directio_enabled)
    fio_fclose(js->directio_fd);

  if (js->directio_ucell_ptr)
    free(js->directio_ucell_ptr);

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
#if JSMAJORVERSION > 1
  js->parsed_structure=0;
  js->directio_block_size=1;
  js->directio_ucell_ptr = NULL;
  js->directio_ucell_blkbuf = NULL;

  js->directio_enabled=0;
  js->ts_file_offset=0;
  js->ts_crd_sz=0;
  js->ts_ucell_sz=0;
  js->ts_crd_padsz=0;
  js->ts_ucell_padsz=0;
#endif

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

#if 1
  if (getenv("VMDJSBLOCKIO")) {
    js->optflags |= JSOPT_TS_BLOCKIO;
    js->directio_block_size = 4096; 
  }
#endif

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

  /* Check to see if block-based trajectory I/O is used  */
  /* and read in the block size for this file.           */
  if (js->optflags & JSOPT_TS_BLOCKIO) {
    fio_fwrite(&js->directio_block_size, sizeof(int), 1, js->fd);
    printf("jsplugin) Block-based I/O enabled: block size %d bytes\n", 
           js->directio_block_size);
  }

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
        atomnames[hashcnt] = (char *) calloc(1, 16 * sizeof(char));
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
        atomtypes[hashcnt] = (char *) calloc(1, 16 * sizeof(char));
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
        resnames[hashcnt] = (char *) calloc(1, 8 * sizeof(char));
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
        segids[hashcnt] = (char *) calloc(1, 8 * sizeof(char));
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
        chains[hashcnt] = (char *) calloc(1, 2 * sizeof(char));
        strcpy(chains[hashcnt], atoms[i].chain);
        hashcnt++;
      }
    }
    for (i=0; i<numchains; i++) {
      fio_fwrite(chains[i], 2 * sizeof(char), 1, js->fd);
    }


printf("jsplugin) writing numeric field tables...\n");
    /* write out all of the atom fields */
    shortbuf = (short *) malloc(js->natoms * sizeof(short));

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
      fltbuf = (float *) malloc(js->natoms * sizeof(float));

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
      intbuf = (int *) malloc(js->natoms * sizeof(int));

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

    /* update the file offset for the first timestep */
    js_calc_timestep_blocking_info(js);

    return MOLFILE_SUCCESS;
  }

  /* update the file offset for the first timestep */
  js_calc_timestep_blocking_info(js);

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

  /* set unit cell pointer to the TS block-aligned buffer area */
  double *unitcell = (double *) js->directio_ucell_blkbuf;

  js->nframes++; /* increment frame count written to the file so far */

  unitcell[0] = ts->A;
  unitcell[1] = ts->B;
  unitcell[2] = ts->C;
  unitcell[3] = sin((M_PI_2 / 90.0) * (90.0 - ts->alpha));
  unitcell[4] = sin((M_PI_2 / 90.0) * (90.0 - ts->beta));
  unitcell[5] = sin((M_PI_2 / 90.0) * (90.0 - ts->gamma));

  /* coordinates for all atoms */
  fio_fwrite(ts->coords, js->ts_crd_padsz, 1, js->fd);

  /* PBC unit cell info */ 
  fio_fwrite(unitcell, js->ts_ucell_padsz, 1, js->fd);

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
  if (js->directio_ucell_ptr)
    free(js->directio_ucell_ptr);

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

#if defined(ENABLECUDATESTS)
#include <cuda_runtime.h>
#endif

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
  float *coords=NULL;
  float *aligncoords=NULL;
  void *v;
  jshandle *js;
  int i, natoms, sz, blocksz;
  float sizeMB =0.0, totalMB = 0.0;
  double starttime, endtime, totaltime = 0.0;
  int do_io = 1;

  printf("Standalone tests for JS plugin:\n");
  
  if (getenv("VMDJSNOIO") != NULL)
    do_io = 0;

  if (do_io)
    printf("  Timestep disk I/O enabled.\n");
  else
    printf("  Timestep disk I/O DISABLED.\n");

#if defined(ENABLECUDATESTS)
  printf("  CUDA GPU support compiled in.\n");

  

  

  cudaError_t crc;
  int maxatomidx=-1;
  int devcount;
  float *devptr=NULL;
  crc = cudaGetDeviceCount(&devcount);

  printf("  GPU device count: %d\n", devcount);
  if (devcount==0)
    printf("  No GPU devices, continuing with host only...\n");

  

  if (getenv("VMDJSCUDATESTS") == NULL) {
    devcount = 0;
    printf("  GPU tests disabled.\n");
    printf("  Enable GPU tests with VMDJSCUDATESTS env variable\n");
  } else {
    printf("  Disable GPU tests by unsetting VMDJSCUDATESTS env variable\n");
  }

#if defined(ENABLEJSSHORTREADS)
  /* test code for an implementation that does short reads that */
  /* skip bulk solvent, useful for faster loading of very large */
  /* structures                                                 */
  if (getenv("VMDJSMAXATOMIDX") != NULL) {
    fio_size_t bszmask;

    maxatomidx = atoi(getenv("VMDJSMAXATOMIDX"));
    if (maxatomidx < 0)
      maxatomidx = 0;
    if (maxatomidx >= js->natoms)
      maxatomidx = js->natoms - 1;

    printf("jsplugin) Short-copies of GPU timesteps enabled: %d / %d atoms (%.2f%%)\n",
           maxatomidx, js->natoms, 100.0*(maxatomidx+1) / ((float) js->natoms));
  }
#endif
#endif


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

    /* ensure we have a large enough allocation so we can align */
    /* the starting pointer to a blocksz page boundary          */
    blocksz = 4096;
    sz = 3*sizeof(float)*natoms + blocksz;
    aligncoords = (float *) alloc_aligned_ptr(sz, blocksz, (void**) &coords);
    timestep.coords = aligncoords;

#if defined(ENABLECUDATESTS)
    if (crc == cudaSuccess && devcount > 0) {
      printf("jsplugin) allocating GPU memory buffer for CUDA tests...\n");
      crc = cudaMalloc((void**) &devptr, sz);
      if (crc != cudaSuccess) {
        printf("Failed to allocate GPU buffer!\n");
        return -1;
      }
    }
#endif

    /* loop over all timesteps ... */
    for (i=0; i<js->nframes; i++) {
      if (do_io) {
        int rc = read_js_timestep(v, natoms, &timestep);
        if (rc) {
          printf("jsplugin) error in read_js_timestep on frame %d\n", i);
          /* return 1; */
        }
      }

#if defined(ENABLECUDATESTS)
      if (crc == cudaSuccess && devcount > 0) {
        size_t bsz = (maxatomidx >= 0) ? (maxatomidx+1) : natoms;
        bsz *= 3*sizeof(float);
        crc = cudaMemcpy(devptr, timestep.coords, bsz, cudaMemcpyHostToDevice);
      }
#endif
    }

#if defined(ENABLECUDATESTS)
    /* wait for any pending GPU calls to complete */
    cudaThreadSynchronize();
#endif

    endtime = time_of_day();
    close_js_read(v);

#if defined(ENABLECUDATESTS)
    if (crc == cudaSuccess && devcount > 0) {
      cudaFree(devptr);
    }
#endif
    free(coords);

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

