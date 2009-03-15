/* MACHINE GENERATED FILE, DO NOT EDIT! */

#define VMDPLUGIN molfile_parmplugin
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
 *      $RCSfile: parmplugin.C,v $
 *      $Author: johns $       $Locker:  $             $State: Exp $
 *      $Revision: 1.32 $       $Date: 2009/02/20 22:28:41 $
 *
 ***************************************************************************/

#include <string.h>
#include "ReadPARM.h"
#include "molfile_plugin.h"

typedef struct {
  ReadPARM *rp;
  FILE *parm;
  int natoms;
  int *from, *to;
} parmdata;

static void *open_parm_read(const char *filename, const char *, 
    int *natoms) {
 
  FILE *parm;
  ReadPARM *rp = new ReadPARM;
  if(!(parm = rp->open_parm_file(filename))) {
    fprintf(stderr, "parmplugin) Cannot open parm file '%s'\n", filename);
    delete rp;
    return NULL;
  }
  
  if (rp->readparm(parm) != 0) {
    delete rp;
    // XXX should we call close_parm_file???
    return NULL; 
  }
  *natoms = rp->get_parm_natoms();
  
  parmdata *p = new parmdata;
  memset(p, 0, sizeof(parmdata));
  p->rp = rp;
  p->parm = parm;
  p->natoms = *natoms;
  return p;
}

static int read_parm_structure(void *mydata, int *optflags,
    molfile_atom_t *atoms) {
  
  parmdata *p = (parmdata *)mydata;
  ReadPARM *rp = p->rp;
  rp->get_parm_boxInfo();
  int i;
 
  *optflags = MOLFILE_CHARGE | MOLFILE_MASS;

  for (i=0; i<p->natoms; i++) {
    molfile_atom_t *atom = atoms+i;
    // XXX Why isn't there a return code for error on read????
    rp->get_parm_atom(i, atom->name, atom->type, atom->resname, atom->segid, 
        &atom->resid, &atom->charge, &atom->mass);
    atom->chain[0] = '\0';
  }
  // XXX amber box info not supported in the API
  return MOLFILE_SUCCESS;
}

static int read_parm_bonds(void *v, int *nbonds, int **fromptr, int **toptr, 
                           float **bondorderptr,  int **bondtype, 
                           int *nbondtypes, char ***bondtypename) {
  parmdata *p = (parmdata *)v;
  ReadPARM *rp = p->rp;
  int i, j;
  int numbonds = rp->get_parm_nbonds();
  p->from = (int *)malloc(numbonds*sizeof(int));
  p->to = (int *)malloc(numbonds*sizeof(int));
  j = 0;
  for (i=0; i<numbonds; i++) {
    int a1, a2;
    rp->get_parm_bond(i, (&a1)-i, (&a2)-i);
    if ( a1 <= p->natoms && a2 <= p->natoms) {
      p->from[j] = a1;
      p->to[j] = a2;
      j++;
    } else {
      printf("parmplugin) skipping bond (%d %d)\n", a1, a2); 
    }
  }
  *nbonds = j;
  *fromptr = p->from;
  *toptr = p->to;
  *bondorderptr = NULL; // PARM files don't have bond order information
  *bondtype = NULL;
  *nbondtypes = 0;
  *bondtypename = NULL;

  return MOLFILE_SUCCESS;
}


static void close_parm_read(void *mydata) {
  parmdata *p = (parmdata *)mydata;
  p->rp->close_parm_file(p->parm);
  if (p->from) free(p->from);
  if (p->to) free(p->to);
  delete p->rp;
}
 
/*
 * Initialization stuff down here
 */

static molfile_plugin_t plugin;

VMDPLUGIN_EXTERN int VMDPLUGIN_init(void) {
  memset(&plugin, 0, sizeof(molfile_plugin_t));
  plugin.abiversion = vmdplugin_ABIVERSION;
  plugin.type = MOLFILE_PLUGIN_TYPE;
  plugin.name = "parm";
  plugin.prettyname = "AMBER Parm";
  plugin.author = "Justin Gullingsrud, John Stone";
  plugin.majorv = 4;
  plugin.minorv = 3;
  plugin.is_reentrant = VMDPLUGIN_THREADSAFE;
  plugin.filename_extension = "parm";
  plugin.open_file_read = open_parm_read;
  plugin.read_structure = read_parm_structure;
  plugin.read_bonds = read_parm_bonds;
  plugin.close_file_read = close_parm_read;
  return VMDPLUGIN_SUCCESS;
}

VMDPLUGIN_EXTERN int VMDPLUGIN_register(void *v, vmdplugin_register_cb cb) {
  (*cb)(v, (vmdplugin_t *)(void *)&plugin);
  return 0;
}

VMDPLUGIN_EXTERN int VMDPLUGIN_fini(void) { return VMDPLUGIN_SUCCESS; }

