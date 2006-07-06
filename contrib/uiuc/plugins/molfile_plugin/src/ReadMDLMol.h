/***************************************************************************
 *cr                                                                       
 *cr            (C) Copyright 1995-2006 The Board of Trustees of the           
 *cr                        University of Illinois                       
 *cr                         All Rights Reserved                        
 *cr                                                                   
 *cr    Portions contributed and copyright (C) 1998 by Andrew Dalke and
 *cr    Bioreason, Inc.                                                
 *cr                                                                   
 *cr    Some information comes from the Daylight Information Systems'
 *cr    contrib program "mol2smi" which was placed into the public domain.
 *cr
 ***************************************************************************/

/***************************************************************************
 * RCS INFORMATION:
 *
 *	$RCSfile: ReadMDLMol.h,v $
 *	$Author: johns $	$Locker:  $		$State: Exp $
 *	$Revision: 1.7 $	$Date: 2003/12/31 20:14:01 $
 *
 ***************************************************************************
 * DESCRIPTION:
 *
 *  Read and write MDL ".mol" files 
 *
 ***************************************************************************/
#ifndef READ_MDL_MOL_H
#define READ_MDL_MOL_H

// Yeah, I know, I should use streams instead of FILE*s
#include <stdio.h>

typedef struct mdl_atom {
  mdl_atom(void) {
    x = 0.0; y = 0.0; z = 0.0; symbol[0] = 0;
    mass_difference = 0; charge = 0; stereo_parity = 0; hcount = 0;
    stereo_care_box = 0;
    valence = 0;
  }
  
  float x, y, z;
  char symbol[4];
  int mass_difference;
  int charge;           // actually, 4 - charge!
  int stereo_parity;
  int hcount;           // actually, hcount + 1 !

  // terms past here exist in some Mol files (according to the Daylight
  // code) but I didn't find examples of use.  These terms are
  // placeholders for future use.
  int stereo_care_box;
  int valence;
} mdl_atom_struct;

typedef struct mdl_bond {
  mdl_bond(void) {
    bond_from = 0; bond_to = 0; bond_type = 0; bond_stereo = 0;
  }
  int bond_from;    // uses base 1
  int bond_to;      // uses base 1
  int bond_type;    // eg, 1, 2, 3, 4
  int bond_stereo;
} mdl_bond_struct;

int read_mdl_header(FILE *infile, int *natoms, int *nbonds);

int read_mdl_atom(FILE *infile, mdl_atom *atom);
int write_mdl_atom(FILE *outfile, const mdl_atom& atom);

int read_mdl_bond(FILE *infile, mdl_bond *bond);
int write_mdl_bond(FILE *outfile, const mdl_bond& bond);

int write_mdl_trailer(FILE *outfile);

#endif

