/***************************************************************************
 *cr
 *cr            (C) Copyright 1995-2002 The Board of Trustees of the
 *cr                        University of Illinois
 *cr                         All Rights Reserved
 *cr
 ***************************************************************************/

/***************************************************************************
 * RCS INFORMATION:
 *
 *      $RCSfile: molfile_plugin.h,v $
 *      $Author: johns $       $Locker:  $             $State: Exp $
 *      $Revision: 1.34 $       $Date: 2006/02/16 20:32:40 $
 *
 ***************************************************************************/

/** @file 
 * API for C extensions to define a way to load structure, coordinate,
 * trajectory, and volumetric data files  
 */ 

#ifndef MOL_FILE_PLUGIN_H
#define MOL_FILE_PLUGIN_H

#include "vmdplugin.h"

/**
 * Define a common plugin type to be used when registering the plugin.
 */
#define MOLFILE_PLUGIN_TYPE "mol file reader"

/**
 * File converter plugins use the same API  but register under a different
 * type so that regular file readers can have priority.
 */
#define MOLFILE_CONVERTER_PLUGIN_TYPE "mol file converter"

/* File plugin symbolic constants for better code readability */
#define MOLFILE_SUCCESS           0   /**< succeeded in reading file      */
#define MOLFILE_EOF              -1   /**< end of file                    */
#define MOLFILE_ERROR            -1   /**< error reading/opening a file   */
#define MOLFILE_NOSTRUCTUREDATA  -2   /**< no structure data in this file */

#define MOLFILE_NUMATOMS_UNKNOWN -1   /**< unknown number of atoms       */
#define MOLFILE_NUMATOMS_NONE     0   /**< no atoms in this file type    */


/**
 * File level comments, origin information, and annotations.
 */
typedef struct {
  char database[81];   /**< database of origin, if any        */
  char accession[81];  /**< database accession code, if any   */
  char date[81];       /**< date/time stamp for this data     */
  char title[81];      /**< brief title for this data         */
  int remarklen;       /**< length of remarks string          */
  char *remarks;       /**< free-form remarks about data      */
} molfile_metadata_t;


/* 
 * Struct for specifying atoms in a molecular structure.  The first 
 * six components are required, the rest are optional and their presence is 
 * indicating by setting the corresponding bit in optsflag.  When omitted,
 * the application (for read_structure) or plugin (for write_structure) 
 * must be able to supply default values if the missing parameters are 
 * part of its internal data structure.
 * Note that it is not possible to specify coordinates with this structure.
 * This is intentional; all coordinate I/O is done with the read_timestep and 
 * write_timestep functions. 
 */

/**
 * Per-atom attributes and information.
 */
typedef struct {
  /* these fields absolutely must be set or initialized to empty */
  char name[16];      /**< required atom name string             */
  char type[16];      /**< required atom type string             */
  char resname[8];    /**< required residue name string          */
  int resid;          /**< required integer residue ID           */
  char segid[8];      /**< required segment name string, or ""   */
  char chain[2];      /**< required chain name, or ""            */

  /* rest are optional; use optflags to specify what's present   */
  char altloc[2];     /**< optional PDB alternate location code  */
  char insertion[2];  /**< optional PDB insertion code           */
  float occupancy;    /**< optional occupancy value              */
  float bfactor;      /**< optional B-factor value               */
  float mass;         /**< optional mass value                   */
  float charge;       /**< optional charge value                 */
  float radius;       /**< optional radius value                 */
  int atomicnumber;   /**< optional element atomic number        */
} molfile_atom_t;

/*@{*/
/** Plugin optional data field availability flag */
#define MOLFILE_NOOPTIONS     0x0000 /**< no optional data                 */
#define MOLFILE_INSERTION     0x0001 /**< insertion codes provided         */
#define MOLFILE_OCCUPANCY     0x0002 /**< occupancy data provided          */
#define MOLFILE_BFACTOR       0x0004 /**< B-factor data provided           */
#define MOLFILE_MASS          0x0008 /**< Atomic mass provided             */
#define MOLFILE_CHARGE        0x0010 /**< Atomic charge provided           */
#define MOLFILE_RADIUS        0x0020 /**< Atomic VDW radius provided       */
#define MOLFILE_ALTLOC        0x0040 /**< Multiple conformations present   */
#define MOLFILE_ATOMICNUMBER  0x0080 /**< Atomic element number provided   */
#define MOLFILE_BONDSSPECIAL  0x0100 /**< Only non-standard bonds provided */
/*@}*/


/*
 * Per-timestep atom coordinates and periodic cell information
 */ 
typedef struct {
  float *coords;  /**< space for coordinates of all atoms, arranged xyzxyzxyz */

  /*@{*/   
  /**
   * Unit cell specification of the form A, B, C, alpha, beta, gamma.
   * notes: A, B, C are side lengths of the unit cell
   * alpha = angle between b and c
   *  beta = angle between a and c
   * gamma = angle between a and b
   */ 
  float A, B, C, alpha, beta, gamma; 
  /*@}*/   
} molfile_timestep_t;


/**
 * Metadata for volumetric datasets, read initially and used for subsequent
 * memory allocations and file loading.  
 */
typedef struct {
  char dataname[256]; /**< name of volumetric data set                    */

  float origin[3];    /**< origin: origin of volume (x=0, y=0, z=0 corner */

  /*
   * x/y/z axis:
   * These the three cell sides, providing both direction and length
   * (not unit vectors) for the x, y, and z axes.  In the simplest
   * case, these would be <size,0,0> <0,size,0> and <0,0,size) for 
   * an orthogonal cubic volume set.  For other cell shapes these
   * axes can be oriented non-orthogonally, and the parallelpiped
   * may have different side lengths, not just a cube/rhombus.
   */
  float xaxis[3];     /**< direction (and length) for X axis              */ 
  float yaxis[3];     /**< direction (and length) for Y axis              */
  float zaxis[3];     /**< direction (and length) for Z axis              */

  /*
   * x/y/z size: 
   * Number of grid cells along each axis.  This is _not_ the
   * physical size of the box, this is the number of voxels in each
   * direction, independent of the shape of the volume set. 
   */
  int xsize;          /**< number of grid cells along the X axis          */
  int ysize;          /**< number of grid cells along the Y axis          */
  int zsize;          /**< number of grid cells along the Z axis          */

  int has_color;      /**< flag indicating presence of voxel color data   */
} molfile_volumetric_t;


/**
 *  Enumeration of all of the supported graphics objects that can be read
 *  from graphics file reader plugins.
 */
enum molfile_graphics_type {
  MOLFILE_POINT,  MOLFILE_TRIANGLE, MOLFILE_TRINORM, MOLFILE_NORMS, 
  MOLFILE_LINE,   MOLFILE_CYLINDER, MOLFILE_CAPCYL,  MOLFILE_CONE,    
  MOLFILE_SPHERE, MOLFILE_TEXT,     MOLFILE_COLOR,   MOLFILE_TRICOLOR
};

/**
 *  Individual graphics object/element data
 */ 
typedef struct {
  int type;             /* One of molfile_graphics_type */
  int style;            /* A general style parameter    */
  float size;           /* A general size parameter     */
  float data[9];        /* All data for the element     */
} molfile_graphics_t;


/*
 * Types for raw graphics elements stored in files.  Data for each type
 * should be stored by the plugin as follows:

type        data                                     style       size
----        ----                                     -----       ----
point       x, y, z                                              pixel size
triangle    x1,y1,z1,x2,y2,z2,x3,y3,z3                 
trinorm     x1,y1,z1,x2,y2,z2,x3,y3,z3                 
            the next array element must be NORMS
tricolor    x1,y1,z1,x2,y2,z2,x3,y3,z3                 
            the next array elements must be NORMS
            the following element must be COLOR, with three RGB triples
norms       x1,y1,z1,x2,y2,z2,x3,y3,z3                 
line        x1,y1,z1,x2,y2,z2                        0=solid     pixel width
                                                     1=stippled
cylinder    x1,y1,z1,x2,y2,z2                        resolution  radius
capcyl      x1,y1,z1,x2,y2,z2                        resolution  radius
sphere      x1,y1,z1                                 resolution  radius
text        x, y, z, up to 24 bytes of text                      pixel size
color       r, g, b
*/


/**
 * Main file reader API.  Any function in this struct may be NULL
 * if not implemented by the plugin; the application checks this to determine
 * what functionality is present in the plugin. 
 */ 
typedef struct {
  /**
   * Required header 
   */
  vmdplugin_HEAD

  /**
   * Filename extension for this file type.  May be NULL if no filename 
   * extension exists and/or is known.  For file types that match several
   * common extensions, list them in a comma separated list such as:
   *  "pdb,ent,foo,bar,baz,ban"
   * The comma separated list will be expanded when filename extension matching
   * is performed.  If multiple plugins solicit the same filename extensions,
   * the first one tried/checked "wins".
   */
  const char *filename_extension;

  /**
   * Try to open the file for reading.  Return an opaque handle, or NULL on
   * failure. Set the number of atoms; if the number of atoms cannot be 
   * determined, set natoms to MOLFILE_NUMATOMS_UNKNOWN. 
   * Filetype should be the name under which this plugin was registered;
   * this is provided so that plugins can provide the same function pointer
   * to handle multiple file types.
   */
  void *(* open_file_read)(const char *filepath, const char *filetype, 
      int *natoms);
  
  /**
   * Read molecular structure from the given file handle.  atoms is allocated
   * by the caller and points to space for natoms.
   * On success, place atom information in the passed-in pointer.  
   * optflags specifies which optional fields in the atoms will be set by
   * the plugin.
   */
  int (*read_structure)(void *, int *optflags, molfile_atom_t *atoms);

  /**
   * Read bond information for the molecule.  On success the arrays from
   * and to should point to the (one-based) indices of bonded atoms.
   * Each unique bond should be specified only once, so file formats that list
   * bonds twice will need post-processing before the results are returned to
   * the caller.
   * If the plugin provides bond information, but the file loaded doesn't 
   * actually contain any bond info, the nbonds parameter should be
   * set to 0 and from/to should be set to NULL to indicate that no bond
   * information was actually present, and automatic bond search should be
   * performed.  
   * If the plugin provides bond order information, the bondorder array
   * will contain the bond order for each from/to pair.  If not, the bondorder
   * pointer should be set to NULL, in which case the caller will provide a 
   * default bond order value of 1.0.
   * These arrays must be freed by the plugin in the close_file_read function.
   * This function can be called only after read_structure().  
   * Return MOLFILE_SUCCESS if no errors occur. 
   */
  int (*read_bonds)(void *, int *nbonds, int **from, int **to, float **bondorder);

  /**
   * Read the next timestep from the file.  Return MOLFILE_SUCCESS, or 
   * MOLFILE_EOF on EOF.  If the molfile_timestep_t argument is NULL, then 
   * the frame should be skipped.  Otherwise, the application must prepare 
   * molfile_timestep_t by allocating space in coords for the corresponding 
   * number of coordinates.  
   * The natoms parameter exists because some coordinate file formats 
   * (like CRD) cannot determine for themselves how many atoms are in a 
   * timestep; the app must therefore obtain this information elsewhere
   * and provide it to the plugin.
   */
  int (* read_next_timestep)(void *, int natoms, molfile_timestep_t *);

  /**
   * Close the file and release all data.  The handle cannot be reused.
   */
  void (* close_file_read)(void *);
   
  /**
   * Open a coordinate file for writing using the given header information.
   * Return an opaque handle, or NULL on failure.  The application must
   * specify the number of atoms to be written. 
   * filetype should be the name under which this plugin was registered.
   */
  void *(* open_file_write)(const char *filepath, const char *filetype, 
      int natoms);
  
  /**
   * Write structure information.  Return success.
   */
  int (* write_structure)(void *, int optflags, const molfile_atom_t *atoms);

  /**
   * Write a timestep to the coordinate file.  Return MOLFILE_SUCCESS if no
   * errors occur.  If the file contains structure information in each 
   * timestep (like a multi-entry PDB), it will have to cache the information 
   * from the initial calls from write_structure.
   */
  int (* write_timestep)(void *, const molfile_timestep_t *);
  
  /**
   * Close the file and release all data.  The handle cannot be reused.
   */
  void (* close_file_write)(void *);

  /**
   * Retrieve metadata pertaining to volumetric datasets in this file.
   * Set nsets to the number of volumetric data sets, and set *metadata
   * to point to an array of molfile_volumetric_t.  The array is owned by
   * the plugin and should be freed by close_file_read().  The application
   * may call this function any number of times.
   */
  int (* read_volumetric_metadata)(void *, int *nsets, 
        molfile_volumetric_t **metadata);

  /** 
   * Read the specified volumetric data set into the space pointed to by 
   * datablock.  The set is specified with a zero-based index.  The space 
   * allocated for the datablock must be equal to
   * xsize * ysize * zsize.  No space will be allocated for colorblock 
   * unless has_color is nonzero; in that case, colorblock should be
   * filled in with three RGB floats per datapoint.
   */
  int (* read_volumetric_data)(void *, int set, float *datablock, 
        float *colorblock);

  /**
   * Read raw graphics data stored in this file.   Return the number of data
   * elements and the data itself as an array of molfile_graphics_t in the 
   * pointer provided by the application.  The plugin is responsible for 
   * freeing the data when the file is closed.
   */
  int (* read_rawgraphics)(void *, int *nelem, const molfile_graphics_t **data);

  /**
   * Read molecule metadata such as what database (if any) this file/data
   * came from, what the accession code for the database is, textual remarks
   * and other notes pertaining to the contained structure/trajectory/volume
   * and anything else that's informative at the whole file level.
   */ 
  int (* read_molecule_metadata)(void *, molfile_metadata_t **metadata);
  
  /**
   * Write bond information for the molecule.  The arrays from
   * and to point to the (one-based) indices of bonded atoms.
   * Each unique bond will be specified only once by the caller. 
   * File formats that list bonds twice will need to emit both the 
   * from/to and to/from versions of each.
   * This function must be called before write_structure().  
   * Unlike the read_bonds() routine, the bondorder information
   * will always be provided (each bondorder can be set to 1 if unknown).
   * Return MOLFILE_SUCCESS if no errors occur. 
   */
  int (* write_bonds)(void *, int nbonds, int *from, int *to, float *bondorder);

} molfile_plugin_t;

#endif

