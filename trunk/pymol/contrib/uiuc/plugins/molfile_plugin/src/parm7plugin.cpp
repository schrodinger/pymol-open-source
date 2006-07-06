/* MACHINE GENERATED FILE, DO NOT EDIT! */

#define VMDPLUGIN molfile_parm7plugin
#define STATIC_PLUGIN 1

/***************************************************************************
 *cr
 *cr            (C) Copyright 1995-2006 The Board of Trustees of the
 *cr                        University of Illinois
 *cr                         All Rights Reserved
 *cr
 ***************************************************************************/

/***************************************************************************
 * RCS INFORMATION:
 *
 *      $RCSfile: parm7plugin.C,v $
 *      $Author: johns $       $Locker:  $             $State: Exp $
 *      $Revision: 1.22 $       $Date: 2006/02/23 19:36:45 $
 *
 ***************************************************************************/

#include <string.h>
#include "molfile_plugin.h"
#include "ReadPARM7.h"

typedef struct {
  parmstruct *prm;
  int popn;
  FILE *fd;
  int nbonds;
  int *from, *to;
} parmdata;

static void *open_parm7_read(const char *filename, const char *,int *natoms) {
  FILE *fd;
  int popn = 0;
  if(!(fd = open_parm7_file(filename, &popn))) {
    fprintf(stderr, "Cannot open parm file '%s'\n", filename);
    return NULL;
  }
  parmstruct *prm = read_parm7_header(fd);
  if (!prm) {
    close_parm7_file(fd, popn);
    return NULL; 
  }

  *natoms = prm->Natom;
  parmdata *p = new parmdata;
  memset(p, 0, sizeof(parmdata));
  p->prm = prm;
  p->popn = popn;
  p->fd = fd;
  p->from = new int[prm->Nbonh + prm->Nbona];
  p->to   = new int[prm->Nbonh + prm->Nbona];
  return p;
}

static int read_parm7_structure(void *mydata, int *optflags, molfile_atom_t *atoms) {
  parmdata *p = (parmdata *)mydata;
  const parmstruct *prm = p->prm;
  FILE *file = p->fd;
  char buf[85];
  char field[85];
  char *resnames = NULL;
  *optflags = 0;

  while (fgets(buf, 85, file)) {
    // find the next line starting with %FLAG, indicating a new section
    if (strncmp(buf, "%FLAG ", 6)) continue;
    // parse field and format indicators
    sscanf(buf+6, "%s\n", field); // type of record
    fscanf(file, "%s\n", buf);    // format
    if (!strcmp(field, "ATOM_NAME")) {
      if (!parse_parm7_atoms(buf, prm->Natom, atoms, file)) break;
    } else if (!strcmp(field, "CHARGE")) {
      *optflags |= MOLFILE_CHARGE;
      if (!parse_parm7_charge(buf, prm->Natom, atoms, file)) break;
    } else if (!strcmp(field, "MASS")) {
      *optflags |= MOLFILE_MASS;
      if (!parse_parm7_mass(buf, prm->Natom, atoms, file)) break;
    } else if (!strcmp(field, "AMBER_ATOM_TYPE")) {
      if (!parse_parm7_atype(buf, prm->Natom, atoms, file)) break;
    } else if (!strcmp(field, "RESIDUE_LABEL")) {
      resnames = new char[4*prm->Nres];
      if (!parse_parm7_resnames(buf, prm->Nres, resnames, file)) break;
    } else if (!strcmp(field, "RESIDUE_POINTER")) {
      if (!resnames) {
        fprintf(stderr, 
            "PARM7: cannot parse RESIDUE_POINTER before RESIDUE_LABEL\n");
        continue;
      }
      if (!parse_parm7_respointers(buf, prm->Natom, atoms, 
                                   prm->Nres, resnames, file)) 
        break;
    } else if (!strcmp(field, "BONDS_WITHOUT_HYDROGEN")) {
      if (!parse_parm7_bonds(buf, prm->Nbona, p->from+p->nbonds,
            p->to+p->nbonds, file)) break;
      p->nbonds += prm->Nbona;
    } else if (!strcmp(field, "BONDS_INC_HYDROGEN")) {
      if (!parse_parm7_bonds(buf, prm->Nbonh, p->from+p->nbonds,
            p->to+p->nbonds, file)) break;
      p->nbonds += prm->Nbonh;
    }
  }

  // unused items
  for (int i=0; i<prm->Natom; i++) {
    atoms[i].chain[0] = '\0';
    atoms[i].segid[0] = '\0';
  }

  delete [] resnames;
  return MOLFILE_SUCCESS;
}

static int read_parm7_bonds(void *v, int *nbonds, int **fromptr, int **toptr, float **bondorderptr){
  parmdata *p = (parmdata *)v;
  *nbonds = p->nbonds;
  *fromptr = p->from;
  *toptr = p->to;
  *bondorderptr = NULL; // parm files don't contain bond order information
  return MOLFILE_SUCCESS;
}

static void close_parm7_read(void *mydata) {
  parmdata *p = (parmdata *)mydata;
  close_parm7_file(p->fd, p->popn);
  delete p->prm;
  delete [] p->from;
  delete [] p->to;
  delete p;
}
 
/*
 * Initialization stuff down here
 */

static molfile_plugin_t parm7plugin = {
  vmdplugin_ABIVERSION,          // ABI Version
  MOLFILE_PLUGIN_TYPE,           // type of plugin
  "parm7",                       // short name of plugin
  "AMBER7 Parm",                 // pretty name of plugin
  "Brian Bennion, Justin Gullingsrud, John E. Stone", // authors
  0,                             // major version
  7,                             // minor version
  VMDPLUGIN_THREADUNSAFE,        // is not reentrant
  "prmtop,parm7",                // filename extensions
  open_parm7_read,        
  read_parm7_structure,  
  read_parm7_bonds,
  0,                             // read_next_timestep
  close_parm7_read,     
  0,
  0,
  0,
  0,
};

VMDPLUGIN_EXTERN int VMDPLUGIN_init(){
  return VMDPLUGIN_SUCCESS;
}
VMDPLUGIN_EXTERN int VMDPLUGIN_fini(){
  return VMDPLUGIN_SUCCESS;
}
VMDPLUGIN_EXTERN int VMDPLUGIN_register(void *v, vmdplugin_register_cb cb) {
  (*cb)(v,(vmdplugin_t *)&parm7plugin);
  return VMDPLUGIN_SUCCESS;
}
