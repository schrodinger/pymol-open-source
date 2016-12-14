/* MACHINE GENERATED FILE, DO NOT EDIT! */

#define VMDPLUGIN molfile_abinitplugin
#define STATIC_PLUGIN 1

/*
 *  ABINIT plugin for VMD
 *  Rob Lahaye, Sungkyunkwan University, Korea
 *  August 2010
 *
 *  LICENSE:
 *    You can include my code for as long as it is part of an 
 *  open source project. Hence: I do not allow my code to be 
 *  part of any closed source project.
 *
 *  ABINIT manual   
 *  http://www.abinit.org/
 * 
 *  LINUX
 *  gcc -O2 -Wall -fPIC -I. -I$VMDBASEDIR/plugins/include -c abinitplugin.c
 *  gcc -shared -o abinitplugin.so abinitplugin.o
 *
 *  MACOSX
 *  c++ -O2 -Wall -I. -I$VMDBASEDIR/plugins/include -c abinitplugin.c
 *  c++ -bundle -o abinitplugin.so abinitplugin.o
 *
 *  Install
 *  copy abinitplugin.so $VMDBASEDIR/plugins/$ARCH/molfile
 */

 /*
  * This plugin does NOT read the general input (*.in) and output (*.out) files of abinit;
  * their syntax is far too flexible and complex to be read in all its varieties.
  * Most (soon hopefully: all) other output files can be read by this plugin. 
  *
  * The VMD generic plugin routines jump to the appropriate routine, depending on the
  * selected file type:
  *
  * open_file_read / read_structure / read_next_timestep / 
  *            |__ GEO_open_file_read/read_structure/read_next_timestep
  *            |__ DEN_POT_WFK_open_file_read/read_structure/read_next_timestep
  *
  * read_volumetric_metadata / read_volumetric_data
  *            |__ DEN_read_volumetric_metadata/read_volumetric_data
  *            |__ POT_read_volumetric_metadata/read_volumetric_data
  *            |__ WFK_read_volumetric_metadata/read_volumetric_data (NOT YET IMPLEMENTED)
  *
  * open_file_write / write_structure / write_timestep / close_file_write
  *            |
  *            |__ the output format is in a sloppy abinit input file format; however,
  *                this file cannot serve as an input file to abinit, although its
  *                syntax is supposed to be easy to understand for a human :).
  */

#include <stdio.h>
#include <stdlib.h>
#if defined(_MSC_VER)
#include <io.h>
#define F_OK 0
// Visual C++ 2005 incorrectly displays a warning about the use of POSIX APIs
// on Windows, which is supposed to be POSIX compliant...
#define access _access
#else
#include <unistd.h>
#endif
#include <math.h>
#include <string.h>

#include "molfile_plugin.h"
#include "periodic_table.h"
#include "unit_conversion.h"

#define LINESIZE 2048  /* maximum length of a line */
#define NATOM_MAX 300  /* maximum number of atoms */

#define DBGPRINT if(1) fprintf

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* Structure to store a file's binary type */
enum Endianness { little_endian, big_endian };
typedef struct {
  enum Endianness endian; /* little or big endian*/
  int recordmarker;       /* 4 or 8 bytes */
} binary_t;


/* Structure to store the header data of an abinit binary file
 * See for example
 *   http://www.abinit.org/documentation/helpfiles/for-v6.0/users/abinit_help.html#6
 */
typedef struct {
  binary_t bintype;
  char codvsn[7];
  int headform,fform;
  int bantot, date, intxc, ixc, natom, ngfft[3], nkpt, npsp,
      nspden, nspinor, nsppol, nsym, ntypat, occopt, pertcase, usepaw;
  double ecut, ecutdg, ecutsm, ecut_eff, qptn[3], rprimd[3][3], stmbias, tphysel, tsmear;
  int usewvl, *istwfk, *nband, *npwarr, *so_psp, *symafm, *symrel[3][3], *typat;
  double *kpt[3], *occ, *tnons[3], *znucltypat, *wtk;

  char title[133];
  double znuclpsp, zionpsp;
  int pspso, pspdat, pspcod, pspxc, lmn_size;

  double residm, *xred[3], etotal, fermie;

  int cplex;
} abinit_binary_header_t;


/* Structure to store important abinit data so that it can be easily shared among the routines */
typedef struct {
  FILE *file;                  /* file pointer for reading or writing */
  char *filename;              /* file name for reading or writing */
  char filetype[4];            /* type of data file: e.g. GEOmetric or DENsity file */

  float rotmat[3][3];          /* rotation matrix, stored for periodic display hack */

  /* variables for atomic structure */
  float rprimd[3][3];          /* Real space PRIMitive translations of the unit cell */
  int natom;                   /* total number of atoms */
  int typat[NATOM_MAX];        /* number of atoms per atom type */
  molfile_atom_t *atomlist;

  /* volumetric variables for data of charge density, wavefunction, etc. */
  int nvolsets;                /* number of volumetric datasets */
  molfile_volumetric_t *vol;   /* volume set metadata */

  abinit_binary_header_t *hdr; /* header info of an abinit binary file */
} abinit_plugindata_t;


static int binread(void *, size_t, FILE *, binary_t);
static abinit_binary_header_t *abinit_header(FILE *);


/* Allocate memory for header */
static abinit_binary_header_t *abinit_header_malloc()
{
  abinit_binary_header_t *hdr = (abinit_binary_header_t *)malloc(sizeof(abinit_binary_header_t));

  /* zero all bytes in the structure */
  if (hdr) memset(hdr, 0, sizeof(abinit_binary_header_t));
  else fprintf(stderr, "\n\nABINIT plugin) ERROR: cannot allocate memory for header.\n");

  return hdr;
}


/* Allocate memory for plugin data */
static abinit_plugindata_t *abinit_plugindata_malloc()
{
  abinit_plugindata_t *data = (abinit_plugindata_t *)malloc(sizeof(abinit_plugindata_t));

  /* zero all bytes in the structure */
  if (data) memset(data, 0, sizeof(abinit_plugindata_t));
  else fprintf(stderr, "\n\nABINIT plugin) ERROR: cannot allocate memory for plugin data.\n");

  return data;
}


/* Free up the header data */
static void abinit_header_free(abinit_binary_header_t *hdr)
{
  int i;

  if (!hdr) return;

  if (hdr->istwfk) free(hdr->istwfk);
  if (hdr->nband) free(hdr->nband);
  if (hdr->npwarr) free(hdr->npwarr);
  if (hdr->so_psp) free(hdr->so_psp);
  if (hdr->symafm) free(hdr->symafm);
  for (i = 0; i < 3; ++i) {
    int j;
    for (j = 0; j < 3; ++j) if (hdr->symrel[i][j]) free(hdr->symrel[i][j]);
    if (hdr->kpt[i]) free(hdr->kpt[i]);
    if (hdr->tnons[i]) free(hdr->tnons[i]);
    if (hdr->xred[i]) free(hdr->xred[i]);
  }
  if (hdr->typat) free(hdr->typat);
  if (hdr->occ) free(hdr->occ);
  if (hdr->znucltypat) free(hdr->znucltypat);
  if (hdr->wtk) free(hdr->wtk);

  free(hdr);
  hdr = NULL;
}


/* Free up the plugin data */
static void abinit_plugindata_free(abinit_plugindata_t *data)
{
  if (!data) return;

  if (data->file) fclose(data->file);
  if (data->filename) free(data->filename);
  if (data->atomlist) free(data->atomlist);
  if (data->vol) free(data->vol);

  abinit_header_free(data->hdr);

  free(data);
  data = NULL;
}


/* This rotation matrix is determined by two vectors.
 * The matrix aligns the first vector along the x-axis
 * and places the second vector in the xy-plane.
 *
 * Input: data-structure, which contains the matrix
 */
static void abinit_buildrotmat(abinit_plugindata_t *data)
{
  float const *const a = data->rprimd[0];
  float const *const b = data->rprimd[1];

  /* Rotate first about y-axis and z-axis to align vector a along the x-axis.
   * phi  : angle between vector a and its projection on the xy-plane
   * theta: angle between projection of vector a on the xy-plane and the x-axis
   */
  const double len   = sqrt(a[0]*a[0] + a[1]*a[1]);
  const double phi   = atan2((double) a[2], (double) len);
  const double theta = atan2((double) a[1], (double) a[0]);

  const double cph = cos(phi);
  const double cth = cos(theta);
  const double sph = sin(phi);
  const double sth = sin(theta);

  /* Rotate about x-axis to place b in the xy-plane. */
  const double psi = atan2(-sph*cth*b[0] - sph*sth*b[1] + cph*b[2],-sth*b[0] + cth*b[1]);
  const double cps = cos(psi);
  const double sps = sin(psi);

  data->rotmat[0][0] =  cph * cth;
  data->rotmat[0][1] =  cph * sth;
  data->rotmat[0][2] =  sph;
  data->rotmat[1][0] = -sth * cps - sph * cth * sps;
  data->rotmat[1][1] =  cth * cps - sph * sth * sps;
  data->rotmat[1][2] =  cph * sps; 
  data->rotmat[2][0] =  sth * sps - sph * cth * cps;
  data->rotmat[2][1] = -cth * sps - sph * sth * cps; 
  data->rotmat[2][2] =  cph * cps;

  DBGPRINT(stderr, "   ROTATION MATRIX: %f   %f   %f\n", data->rotmat[0][0], data->rotmat[0][1], data->rotmat[0][2]);
  DBGPRINT(stderr, "                    %f   %f   %f\n", data->rotmat[1][0], data->rotmat[1][1], data->rotmat[1][2]);
  DBGPRINT(stderr, "                    %f   %f   %f\n", data->rotmat[2][0], data->rotmat[2][1], data->rotmat[2][2]);
}


/* Read a non-empty line from stream, remove comments (#... or !...)
 * and strip redundant whitespaces.
 *
 * Input: string, which can contain the line
 *        input stream
 *
 * Return: the stripped line, or NULL when EOF is reached
 */
static char *abinit_readline(char *line, FILE *stream)
{
  char *lineptr;

  if (!line || !stream) return NULL;

  do {
    int i;
    char *cptr;

    /* read one line from the stream */
    lineptr = fgets(line, LINESIZE, stream);

    /* first remove comment from the line */
    for (i = 0; i < strlen(line); ++i) {
        if (line[i] == '#' || line[i] == '!') {line[i] = '\0'; break;}
    }

    /* next remove redundant white spaces at the end of the line */
    for (cptr = &line[strlen(line) - 1]; isspace(*cptr); --cptr) *cptr = '\0';

    /* continue for as long as EOF is not reached and the line is empty */
  } while (lineptr != NULL && strlen(line) == 0);

  return lineptr;
}


/* Abinit uses and generates several types of files:
 *   foobar_GEO: geometry file (formatted/ascii text)
 *   foobar_DEN: charge density file (unformatted/binary)
 *   foobar_WFK: wavefunction file (unformatted/binary)
 *   foobar_POT: potential file (unformatted/binary)
 *
 * If the filetype is already set, then we only compare the
 * filetype with the given string, otherwise we first find
 * out about the filetype and then do the compare.
 *
 * Input: abinit-data-structure (may or may not have the filetype set)
          string to compare with
 *
 * Return: comparison result: 1 (equal), 0 (not equal or error)
 */
static int abinit_filetype(abinit_plugindata_t *data, char const *cmp)
{
  char lineptr[LINESIZE];

  if (!data || !cmp) return 0;

  /* if filetype is already set, then only compare */
  if (strlen(data->filetype) != 0) return (strncmp(data->filetype, cmp, 3) == 0);

  /* first try to read the abinit binary header */
  data->hdr = abinit_header(data->file);
  if (data->hdr) {
    /* header is read successfully,
     * therefore it must be an abinit unformatted (binary) file
     */

    switch(data->hdr->fform) {
    case   2: /* WFK Wave Function file */
              strcpy(data->filetype, "WFK");
              break;
    case  52: /* DEN Charge Density file */
              strcpy(data->filetype, "DEN");
              break;
    case 102: /* POT Potential file */
              strcpy(data->filetype, "POT");
              break;
    default:  /* Error */
              strcpy(data->filetype, "ERR");
              break;
    }

  } else {
    /* reading the abinit binary header failed;
     * we now resort to formatted (ascii text) file.
     */

    /* read the first non-empty line of the file */
    rewind(data->file);
    abinit_readline(lineptr, data->file);

    if (strstr(lineptr, " GEO file"))
      strcpy(data->filetype, "GEO"); /* GEO geometry file */
    else
      strcpy(data->filetype, "ERR"); /* Error */

    /* set the file pointer back to the beginning of the file */
    rewind(data->file);
  }

  return (strncmp(data->filetype, cmp, 3) == 0);
}


/* Construct a new file name by finding the first integer number
 * in the existing file name, and create the same file name but
 * with this integer number incremented.
 *
 * Note: the string 'filename' is changed into the new filename
 * ONLY if the file with the new filename exists! 
 *
 * Input: filename
 *
 * Return: 0 (success), 1 (no more files), 2 (error)
 */
static int increment_filename(char *filename)
{
  int i;
  char *newfilename = NULL, *endpart = NULL;

  DBGPRINT(stderr, "Enter increment_filename\n");

  DBGPRINT(stderr, "increment_filename: filename = %s \n", filename);
  /* search for integer in the filename starting from the end of the string */
  for (i = strlen(filename) - 1; i >= 0 && !newfilename; --i) {

    /* endpart points to the string AFTER the integer */
    if (!endpart && isdigit(filename[i])) endpart = strdup(filename + i + 1);

    /* allocate newfilename when integer is found */
    if (endpart && !newfilename && !isdigit(filename[i])) {
       newfilename = (char *)malloc(sizeof(char) * (2 + strlen(filename)));
       if (!newfilename) {
         free(endpart);
         return 2;
       }

       /* first copy part of the string BEFORE the integer number */
       strncpy(newfilename, filename, i + 1);

       /* second append the incremented integer number and add 'endpart' of the filename */
       sprintf(newfilename + i + 1, "%d%s", 1 + atoi(filename + i + 1), endpart);
    }
  }

  /* if the endpart is not found, then the file name does not have a number in it,
   * and therefore this filename is not part of a series of timesteps.
   * note that if this occurs, then 'newfilename' must also still be NULL !
   */
  if (!endpart) {
    DBGPRINT(stderr, "Exit increment_filename\n");
    return 1;
  }
  
  /* clean up */
  free(endpart);

  /* check if we can access the new file */
  if (access(newfilename, F_OK) != 0) {
    /* the new file cannot be accessed, so it does not exist */
    free(newfilename);
    DBGPRINT(stderr, "Exit increment_filename\n");
    return 1;
  } else {
    /* the new filename exists! Replace the old file name. */
    strcpy(filename, newfilename);
    free(newfilename);
    DBGPRINT(stderr, "increment_filename: filename = %s \n", filename);
    DBGPRINT(stderr, "Exit increment_filename\n");
    return 0;
  }
}


/* Geometry files are generated by ABINIT when "prtgeo 1" is given in
 * the input file. Note that this is NOT the abinit default!
 *
 * Input: data-structure
 *        reference to the number of atoms
 *
 * Return: data-structure, or NULL on error
 */
static void *GEO_open_file_read(abinit_plugindata_t *data, int *natoms)
{
  char lineptr[LINESIZE], atomname[NATOM_MAX][10];
  int i, idx;

  DBGPRINT(stderr, "Enter GEO_open_file_read\n");

  /* go to the line with the text 'XMOL data' */
  while (abinit_readline(lineptr, data->file) != NULL) {
    if (strstr(lineptr, "XMOL data")) break;
  }
  if (!strstr(lineptr, "XMOL data")) {
    fprintf(stderr, "\n\nABINIT read) ERROR: '%s' has no 'XMOL data...' lines.\n", data->filename);
    return NULL;
  }

  /* next non-empty line has the number of atoms */
  if (abinit_readline(lineptr, data->file) == NULL) {
    fprintf(stderr, "\n\nABINIT read) ERROR: cannot find the number of atoms in file '%s'.\n", data->filename);
    return NULL;
  }
  data->natom = atoi(lineptr);
  if (data->natom <= 0 || data->natom > NATOM_MAX) {
    fprintf(stderr, "\n\nABINIT read) ERROR: file '%s' has %d number of atoms.\n", data->filename, data->natom);
    return NULL;
  }

  /* read through the list of atom names and ignore their positions */
  for (i = 0; i < NATOM_MAX; ++i) data->typat[i] = atomname[i][0] = '\0';
  for (idx = i = 0; i < data->natom; ++i) {
    int n;
    char name[10];
    if (1 != fscanf(data->file, "%s %*f %*f %*f", name)) {
      fprintf(stderr, "\n\nABINIT read) ERROR: file '%s' does not have the atom list.\n", data->filename);
      return NULL;
    }

    /* compare current atom name with previous read list and set typat accordingly */
    for (n = 0; n < idx; ++n) if (strcmp(atomname[n], name) == 0) break;
    if (n == idx) strcpy(atomname[idx++], name);
    data->typat[i] = n + 1;

    DBGPRINT(stderr, "   \"%s\": name = %s : data->typat[%d] = %d\n", data->filetype, atomname[n], i, data->typat[i]);
  }

  rewind(data->file);

  *natoms = data->natom;

  DBGPRINT(stderr, "Exit GEO_open_file_read\n");
  return data;
}


static int GEO_read_structure(abinit_plugindata_t *data, int *optflags, molfile_atom_t *atomlist)
{
  char lineptr[LINESIZE];
  int i, status;

  DBGPRINT(stderr, "Enter GEO_read_structure\n");

  /* go to the line with the 'XMOL data' */
  do {
    char *line = abinit_readline(lineptr, data->file);
    status = line != NULL && !strstr(lineptr, "XMOL data");
  } while (status);

  /* skip line with atom numbers */
  abinit_readline(lineptr, data->file);

  /* find atom types in XMOL list */
  for (i = 0; i < data->natom; ++i) {
    molfile_atom_t *const atom = &(atomlist[i]);

    /* required fields */
    if (1 != fscanf(data->file, "%s %*f %*f %*f", atom->name)) {
      fprintf(stderr, "\n\nABINIT read) ERROR: file '%s' does not have the atom list.\n", data->filename);
      return MOLFILE_ERROR;
    }
    strncpy(atom->type, atom->name, sizeof(atom->type));
    atom->resname[0] = '\0';
    atom->resid = 1;
    atom->segid[0]='\0';
    atom->chain[0]='\0';

    /* Optional fields (defined in *optflags) */
    atom->atomicnumber = get_pte_idx(atom->name);
    atom->mass = get_pte_mass(atom->atomicnumber);
    atom->radius = get_pte_vdw_radius(atom->atomicnumber);

    DBGPRINT(stderr, "   atom %d : %d (%s)\n", i, atom->atomicnumber, atom->name);
  }

  /* tell which of the optional fields in the molfile_atom_t structure are provided */
  *optflags = MOLFILE_ATOMICNUMBER | MOLFILE_MASS | MOLFILE_RADIUS; 

  rewind(data->file);

  DBGPRINT(stderr, "Exit GEO_read_structure\n");
  return MOLFILE_SUCCESS;
}


static int GEO_read_next_timestep(abinit_plugindata_t *data, int natoms, molfile_timestep_t *ts)
{
  char lineptr[LINESIZE];
  float *a, *b, *c;
  int i , status;

  DBGPRINT(stderr, "Enter GEO_read_next_timestep\n");

  /* At the very first call the file pointer will be non-NULL.
   * All consecutive calls will have a NULL file pointer (see at the end
   * of this routine) and in that case "increment" the filename for the
   * data of the next timestep; if this file does not exist, then there
   * are no more timesteps.
   */
  if (!data->file) {
    if (increment_filename(data->filename) != 0) return MOLFILE_EOF;

    data->file = fopen(data->filename, "r");
    if (!data->file) return MOLFILE_EOF;
  }

  DBGPRINT(stderr, "GEO_read_next_timestep: filename = %s \n", data->filename);

  /* go to the line with 'Primitive vectors...' */
  do {
    char *line = abinit_readline(lineptr, data->file);
    status = ( line != NULL && !strstr(lineptr, "Primitive vectors") );
  } while (status);

  /* read unit cell vectors from file */
  for (i = 0; i < 3; ++i) {
    float length, *r = data->rprimd[i];
    if (3 != fscanf(data->file, "%*s %f %f %f", &r[0], &r[1], &r[2])) return MOLFILE_EOF;

    /* convert length units from Bohr to Angstrom */
    r[0] *= BOHR_TO_ANGS;
    r[1] *= BOHR_TO_ANGS;
    r[2] *= BOHR_TO_ANGS;

    /* lengths of the respective unit cell vector */
    length = sqrt(r[0]*r[0] + r[1]*r[1] + r[2]*r[2]);
    switch (i) {
      case 0: ts->A = length; break;
      case 1: ts->B = length; break;
      case 2: ts->C = length; break;
    }
  }
  abinit_buildrotmat(data);

  /* determine angles between the vectors of the unit cell */
  a = data->rprimd[0];
  b = data->rprimd[1];
  c = data->rprimd[2];

  /* alpha: angle (in degrees) between b and c */
  ts->alpha = (180.0/M_PI) * acos( (b[0]*c[0] + b[1]*c[1] + b[2]*c[2]) / (ts->B*ts->C) );

  /* beta: angle (in degrees) between a and c */
  ts->beta  = (180.0/M_PI) * acos( (a[0]*c[0] + a[1]*c[1] + a[2]*c[2]) / (ts->A*ts->C) );

  /* gamma: angle (in degrees) between a and b */
  ts->gamma = (180.0/M_PI) * acos( (a[0]*b[0] + a[1]*b[1] + a[2]*b[2]) / (ts->A*ts->B) );

  for (i = 0; i < 9; ++i) DBGPRINT(stderr, "   data->rprimd[%d][%d] = %f %s", i%3, i/3, data->rprimd[i%3][i/3], ((i+1)%3 == 0 ? "\n" : ""));

  /* go to the line with 'XMOL data' */
  do {
    char *line = abinit_readline(lineptr, data->file);
    status = line != NULL && !strstr(lineptr, "XMOL data");
  } while (status);

  /* skip line with atom numbers */
  abinit_readline(lineptr, data->file);

  /* read coordinates (in Angstrom) of each atom */
  for (i = 0; i < data->natom; ++i) {
    float x, y, z, *coords = &ts->coords[3*i];
    fscanf(data->file, "%*s %f %f %f", &x, &y, &z);
    coords[0] = data->rotmat[0][0]*x + data->rotmat[0][1]*y + data->rotmat[0][2]*z;
    coords[1] = data->rotmat[1][0]*x + data->rotmat[1][1]*y + data->rotmat[1][2]*z;
    coords[2] = data->rotmat[2][0]*x + data->rotmat[2][1]*y + data->rotmat[2][2]*z;
  }

  /* close the file and NULLify as preparation for the next timestep call */
  fclose(data->file);
  data->file = NULL;

  DBGPRINT(stderr, "Exit GEO_read_next_timestep\n");
  return MOLFILE_SUCCESS;
}


/* Electron density (DEN), potential (POT), and wavefunction (WFK)
 * files have the structure information in the generic header,
 * followed by the volumetric data. For all three we can determine
 * the structure with the same reading routine.
 * The volumetric data is the similar for density and potential (DEN/POT)
 * files, but is different for wavefunction (WFK) files.
 */
static void *DEN_POT_WFK_open_file_read(abinit_plugindata_t *data, int *natoms)
{
  int i;

  DBGPRINT(stderr, "Enter DEN_POT_WFK_open_file_read\n");

  data->natom = data->hdr->natom;

  if (data->natom <= 0 || data->natom > NATOM_MAX) return NULL;

  for (i = 0; i < data->natom; ++i) data->typat[i] = data->hdr->typat[i];

  for (i = 0; i < data->natom; ++i) DBGPRINT(stderr, "   \"%s\": data->typat[%d] = %d\n", data->filetype, i, data->typat[i]);

  *natoms = data->natom;

  DBGPRINT(stderr, "Exit DEN_POT_WFK_open_file_read\n");
  return data;
}


static int DEN_POT_WFK_read_structure(abinit_plugindata_t *data, int *optflags, molfile_atom_t *atomlist)
{
  int i;

  DBGPRINT(stderr, "Enter DEN_POT_WFK_read_structure\n");

  /* set the atom types, names, etc. */
  for (i = 0; i < data->natom; ++i) {
    molfile_atom_t *const atom = &(atomlist[i]);

    /* optional fields (defined in *optflags) */
    atom->atomicnumber = (int)floor(0.5 + data->hdr->znucltypat[data->hdr->typat[i] - 1]);
    atom->mass = get_pte_mass(atom->atomicnumber);
    atom->radius = get_pte_vdw_radius(atom->atomicnumber);

    /* required fields */
    strncpy(atom->name, get_pte_label(atom->atomicnumber), sizeof(atom->name));    
    strncpy(atom->type, atom->name, sizeof(atom->type));
    atom->resname[0] = '\0';
    atom->resid = 1;
    atom->segid[0]='\0';
    atom->chain[0]='\0';

    DBGPRINT(stderr, "   atom %d : %d (%s)\n", i, atom->atomicnumber, atom->name);
  }

  /* tell which of the optional fields in the molfile_atom_t structure are provided */
  *optflags = MOLFILE_ATOMICNUMBER | MOLFILE_MASS | MOLFILE_RADIUS; 

  DBGPRINT(stderr, "Exit DEN_POT_WFK_read_structure\n");
  return MOLFILE_SUCCESS;
}


static int DEN_POT_WFK_read_next_timestep(abinit_plugindata_t *data, int natoms, molfile_timestep_t *ts)
{
  float *a, *b, *c;
  int i;

  DBGPRINT(stderr, "Enter DEN_POT_WFK_read_next_timestep\n");

  /* At the very first call the file pointer will be non-NULL.
   * All consecutive calls will have a NULL file pointer (see at the end
   * of this routine) and in that case "increment" the filename for the
   * data of the next timestep; if this file does not exist, then there
   * are no more timesteps.
   */
  if (!data->file) {
    /* alas, volumetric data cannot (yet ?) be read as timesteps */
    return MOLFILE_EOF;

    if (increment_filename(data->filename) != 0) return MOLFILE_EOF;

    data->file = fopen(data->filename, "r");
    if (!data->file) return MOLFILE_EOF;

    /* read the new header info from the new file */
    abinit_header_free(data->hdr);
    data->hdr = abinit_header(data->file);
    if (!data->hdr) return MOLFILE_EOF;
  }

  /* get unit cell vectors and convert them to Angstrom */
  for (i = 0; i < 3; ++i) {
    float length;
    int k;
    for (k = 0; k < 3; ++k) data->rprimd[i][k] = data->hdr->rprimd[i][k] * BOHR_TO_ANGS;
    length = sqrt(pow(data->rprimd[i][0], 2) + pow(data->rprimd[i][1], 2) + pow(data->rprimd[i][2], 2));
    switch (i) {
      case 0: ts->A = length; break;
      case 1: ts->B = length; break;
      case 2: ts->C = length; break;
    }    
  }
  abinit_buildrotmat(data);

  for (i = 0; i < 9; ++i) DBGPRINT(stderr, "   data->rprimd[%d][%d] = %f %s", i%3, i/3, data->rprimd[i%3][i/3], ((i+1)%3 == 0 ? "\n" : ""));

  /* determine angles between the vectors of the unit cell */
  a = data->rprimd[0];
  b = data->rprimd[1];
  c = data->rprimd[2];

  /* alpha: angle (in degrees) between b and c */
  ts->alpha = (180.0/M_PI) * acos( (b[0]*c[0] + b[1]*c[1] + b[2]*c[2]) / (ts->B*ts->C) );

  /* beta: angle (in degrees) between a and c */
  ts->beta  = (180.0/M_PI) * acos( (a[0]*c[0] + a[1]*c[1] + a[2]*c[2]) / (ts->A*ts->C) );

  /* gamma: angle (in degrees) between a and b */
  ts->gamma = (180.0/M_PI) * acos( (a[0]*b[0] + a[1]*b[1] + a[2]*b[2]) / (ts->A*ts->B) );

  /* get coordinates (in Angstrom) of each atom */
  for (i = 0; i < data->natom; ++i) {
    double **xred = data->hdr->xred;
    float *coords = &ts->coords[3*i];
    float const x = xred[0][i] * a[0] + xred[1][i] * b[0] + xred[2][i] * c[0],
                y = xred[0][i] * a[1] + xred[1][i] * b[1] + xred[2][i] * c[1],
		z = xred[0][i] * a[2] + xred[1][i] * b[2] + xred[2][i] * c[2];

    coords[0] = data->rotmat[0][0]*x + data->rotmat[0][1]*y + data->rotmat[0][2]*z;
    coords[1] = data->rotmat[1][0]*x + data->rotmat[1][1]*y + data->rotmat[1][2]*z;
    coords[2] = data->rotmat[2][0]*x + data->rotmat[2][1]*y + data->rotmat[2][2]*z;
  }

  /* close the file and NULLify as preparation for the next timestep call */
  fclose(data->file);
  data->file = NULL;

  DBGPRINT(stderr, "Exit DEN_POT_WFK_read_next_timestep\n");
  return MOLFILE_SUCCESS;
}


static int DEN_read_volumetric_metadata(abinit_plugindata_t *data, int *nvolsets, molfile_volumetric_t **metadata)
{
  char const spintext[3][30] = { "Charge density spin up",
                                 "Charge density spin down",
                                 "Charge density spin up - down" };
  char const magntext[3][40] = { "X-projection of local magnetization",
                                 "Y-projection of local magnetization",
                                 "Z-projection of local magnetization" };
  int i;

  DBGPRINT(stderr, "Enter DEN_read_volumetric_metadata\n");

  /* Abinit provides the electron density, which is the negative of the
   * charge density. For norm-conserving pseudopotentials (usepaw = 0)
   * this is the pseudo-valence electron density; for a PAW calculation
   * (usepaw = 1) this is the pseudo-valence electron density plus the
   * compensation charge density.
   * To visualize the electron density from a PAW calculation, create
   * the density file with the pawprtden key word, not prtden (however,
   * if you want to chain this density into later calculations, you have
   * to use prtden.
   */
  if (data->hdr->usepaw) {
    fprintf(stderr, "\n\nABINIT read) WARNING: be sure that you have used \"pawprtden 1\"\n");
    fprintf(stderr,     "                      in order to visualize the electron density!\n\n");
  }

  /* Initialize the volume set list
   * Total charge density (spin up+down) : always present
   *    spin up / spin down / spin up-down \ only there for
   *    X / Y / Z local magnetization      /  spin-polarized calculations
   *                    (the latter two are absent for non-spin-polarized calculations)
   */
  data->nvolsets = (data->hdr->nspden == 1 ? 1 : 4);
  data->vol = (molfile_volumetric_t *)malloc(data->nvolsets * sizeof(molfile_volumetric_t));
  if (!data->vol) {
    fprintf(stderr, "\n\nABINIT read) ERROR: cannot allocate space for volumetric data.\n");
    return MOLFILE_ERROR;
  }

  /* get unit cell vectors and convert to Angstrom */
  for (i = 0; i < 3; ++i) {
    int k;
    for (k = 0; k < 3; ++k) data->rprimd[i][k] = data->hdr->rprimd[i][k] * BOHR_TO_ANGS;
  }
  abinit_buildrotmat(data);

  for (i = 0; i < 9; ++i) DBGPRINT(stderr, "   data->rprimd[%d][%d] = %f %s", i%3, i/3, data->rprimd[i%3][i/3], ((i+1)%3 == 0 ? "\n" : ""));

  for (i = 0; i < data->nvolsets; ++i) {

    int k;

    /* handle to the current volume set meta data */
    molfile_volumetric_t *const set = &(data->vol[i]);

    /* set the volume data name */
    if (i == 0) {
        /* the first data set always is the total charge density */
        sprintf(set->dataname, "Total charge density");
    } else if (data->hdr->nspden <= 2) {
        /* datasets contain spin polarized (up/down) densities */
        sprintf(set->dataname, "%s", spintext[i-1]);
    } else if (data->hdr->nspden == 4) {
        /* subsequent datasets contain magnetization data */
        sprintf(set->dataname, "%s", magntext[i-1]);
    } else {
        /* we should never get here.... */
        sprintf(set->dataname, "%s", "ERROR: no datasets available");
    }

    /* rotate unit cell vectors */
    for (k = 0 ; k < 3; ++k) {
      set->xaxis[k] = data->rotmat[k][0] * data->rprimd[0][0]
		    + data->rotmat[k][1] * data->rprimd[0][1]
		    + data->rotmat[k][2] * data->rprimd[0][2];

      set->yaxis[k] = data->rotmat[k][0] * data->rprimd[1][0] 
		    + data->rotmat[k][1] * data->rprimd[1][1]
		    + data->rotmat[k][2] * data->rprimd[1][2];

      set->zaxis[k] = data->rotmat[k][0] * data->rprimd[2][0] 
		    + data->rotmat[k][1] * data->rprimd[2][1]
		    + data->rotmat[k][2] * data->rprimd[2][2];
    }
    DBGPRINT(stderr, "   set->xaxis[%d] set->yaxis[%d] set->zaxis[%d]\n", k, k, k);
    for (k = 0 ; k < 3; ++k) DBGPRINT(stderr, "   %f         %f        %f\n", set->xaxis[k], set->yaxis[k], set->zaxis[k]);

    /* Add one more to the grid size and later fill the extra voxel
     * with the same value from the beginning of the row to make
     * the volumetric data smooth across the cell boundaries.
     */
    set->xsize = data->hdr->ngfft[0] + 1;
    set->ysize = data->hdr->ngfft[1] + 1;
    set->zsize = data->hdr->ngfft[2] + 1;

    set->has_color = 0;
    set->origin[0] = set->origin[1] = set->origin[2] = 0;
  }

  *nvolsets = data->nvolsets;
  *metadata = data->vol;  

  DBGPRINT(stderr, "Exit DEN_read_volumetric_metadata.\n");
  return MOLFILE_SUCCESS;
}


static int DEN_read_volumetric_data(abinit_plugindata_t *data, int set, float *datablock, float *colorblock)
{
  double const density_conversion = 1.0/pow(BOHR_TO_ANGS, 3);
  int iset;

  DBGPRINT(stderr, "Enter DEN_read_volumetric_data\n");

  /* this should never happen, nevertheless check it to prevent disasters */
  if (set >= data->nvolsets) return MOLFILE_ERROR;

  /* The density output file contains the header plus:
   * do ispden=1,nspden
   * write(unit) (rhor(ir),ir=1,cplex*ngfft(1)*ngfft(2)*ngfft(3))
   * enddo
   */
  for (iset = 0; iset <= set && iset < data->hdr->nspden; ++iset) {
    int const xsize = data->vol[iset].xsize; 
    int const ysize = data->vol[iset].ysize;
    int const zsize = data->vol[iset].zsize;

    char recordmarker[10];
    int n, ix, iy, iz;

    /* Fill the datablock with the density or magnetization data.
     * Note that for each 'iset'-loop the datablock is overwritten,
     * so that only the last datablock readings remain.
     */
    for (n = iz = 0; iz < zsize; ++iz) {
      for (iy = 0; iy < ysize; ++iy) {
        for (ix = 0; ix < xsize; ++ix, ++n) {
          double value;

          /* The datablock grid is one voxel larger in each direction in order
	   * to smoothen the density at the cell boundaries; hence fill the extra
	   * voxel with the same value as the one at the beginning of that row.
	   */
          if (ix == xsize - 1) value = datablock[n - ix];
	  else if (iy == ysize - 1) value = datablock[n - iy*xsize];
	  else if (iz == zsize - 1) value = datablock[n - iz*ysize*xsize];
	  else if (data->hdr->cplex == 1) {
            /* get the volumetric data and convert */
            binread(&value, 8, data->file, data->hdr->bintype);
	    value *= density_conversion;
	  } else if (data->hdr->cplex == 2) {
            /* get volumetric data as a complex number and convert */
	    double a, b;
            binread(&a, 8, data->file, data->hdr->bintype);
            binread(&b, 8, data->file, data->hdr->bintype);
	    value = sqrt(a*a + b*b) * density_conversion;
	  } else {
            /* we should never get here */
            return MOLFILE_ERROR;
          }

          /* fill the datablock according to the data provided */
          if (data->hdr->nspden <= 2) {
            /* data is charge density as one or two sets: total and spin up */

            if (set == 0 || set == 1) {
              /* first two sets are always directly provided as total and spin up */
              datablock[n] = value;
            } else if (set == 2) {
              /* third set we calculate from first two sets: down = total - up */
              datablock[n] = (iset == 0 ? value : datablock[n] - value);
            } else if (set == 3) {
              /* fourth set we calculate from first two sets: up - down = 2*up - total */
              datablock[n] = (iset == 0 ? -value : datablock[n] + 2*value);
            } else {
              /* we should never get here */
              return MOLFILE_ERROR;
            }

          } else if (data->hdr->nspden == 4) {
            /* data is magnetization as four sets: total charge density and X-Y-Z-magnetization */
            datablock[n] = value;

          } else {
            /* we should never get here...*/
            return MOLFILE_ERROR;
          }

	}
      }
    }

    /* skip the recordmarker bytes twice */
    fread(recordmarker, 1, data->hdr->bintype.recordmarker, data->file);
    fread(recordmarker, 1, data->hdr->bintype.recordmarker, data->file);
  }

  DBGPRINT(stderr, "Exit DEN_read_volumetric_data\n");
  return MOLFILE_SUCCESS;
}


static int POT_read_volumetric_metadata(abinit_plugindata_t *data, int *nvolsets, molfile_volumetric_t **metadata)
{
  int i;

  DBGPRINT(stderr, "Enter POT_read_volumetric_metadata\n");

  data->nvolsets = data->hdr->nspden;

  data->vol = (molfile_volumetric_t *)malloc(data->nvolsets * sizeof(molfile_volumetric_t));
  if (!data->vol) {
    fprintf(stderr, "\n\nABINIT read) ERROR: cannot allocate space for volumetric data.\n");
    return MOLFILE_ERROR;
  }

  /* get unit cell vectors and convert to Angstrom */
  for (i = 0; i < 3; ++i) {
    int k;
    for (k = 0; k < 3; ++k) data->rprimd[i][k] = data->hdr->rprimd[i][k] * BOHR_TO_ANGS;
  }
  abinit_buildrotmat(data);

  for (i = 0; i < 9; ++i) DBGPRINT(stderr, "   data->rprimd[%d][%d] = %f %s", i%3, i/3, data->rprimd[i%3][i/3], ((i+1)%3 == 0 ? "\n" : ""));

  for (i = 0; i < data->nvolsets; ++i) {

    int k;

    /* handle to the current volume set meta data */
    molfile_volumetric_t *const set = &(data->vol[i]);

    if (data->nvolsets == 1) strcpy(set->dataname, "Total potential");
    else if (data->nvolsets == 2) {
        if (i == 0) strcpy(set->dataname, "Spin up potential");
        if (i == 1) strcpy(set->dataname, "Spin down potential");
    } else if (data->nvolsets == 4) {
        if (i == 0) strcpy(set->dataname, "Spin up-up potential");
        if (i == 1) strcpy(set->dataname, "Spin down-down potential");
        if (i == 2) strcpy(set->dataname, "Real part of spin up-down potential");
        if (i == 3) strcpy(set->dataname, "Imaginary part of spin up-down potential");
    }

    /* rotate unit cell vectors */
    for (k = 0 ; k < 3; ++k) {
      set->xaxis[k] = data->rotmat[k][0] * data->rprimd[0][0]
		    + data->rotmat[k][1] * data->rprimd[0][1]
		    + data->rotmat[k][2] * data->rprimd[0][2];

      set->yaxis[k] = data->rotmat[k][0] * data->rprimd[1][0] 
		    + data->rotmat[k][1] * data->rprimd[1][1]
		    + data->rotmat[k][2] * data->rprimd[1][2];

      set->zaxis[k] = data->rotmat[k][0] * data->rprimd[2][0] 
		    + data->rotmat[k][1] * data->rprimd[2][1]
		    + data->rotmat[k][2] * data->rprimd[2][2];
    }
    DBGPRINT(stderr, "   set->xaxis[%d] set->yaxis[%d] set->zaxis[%d]\n", k, k, k);
    for (k = 0 ; k < 3; ++k) DBGPRINT(stderr, "   %f         %f        %f\n", set->xaxis[k], set->yaxis[k], set->zaxis[k]);

    /* Add one more to the grid size and fill the extra voxel
     * with the same value at the beginning of the row to make
     * the volumetric data smooth across the cell boundaries.
     */
    set->xsize = data->hdr->ngfft[0] + 1;
    set->ysize = data->hdr->ngfft[1] + 1;
    set->zsize = data->hdr->ngfft[2] + 1;

    set->has_color = 0;
    set->origin[0] = set->origin[1] = set->origin[2] = 0;
  }

  *nvolsets = data->nvolsets;
  *metadata = data->vol;  

  DBGPRINT(stderr, "Exit POT_read_volumetric_metadata.\n");
  return MOLFILE_SUCCESS;
}


static int POT_read_volumetric_data(abinit_plugindata_t *data, int set, float *datablock, float *colorblock)
{
  int n, iset;

  DBGPRINT(stderr, "Enter POT_read_volumetric_data\n");

  /* this should never happen, nevertheless check it to prevent disasters */
  if (set >= data->nvolsets) return MOLFILE_ERROR;

  /* The potential files contain the header plus:
   * do ispden=1,nspden
   * write(unit) (potential(ir),ir=1,cplex*ngfft(1)*ngfft(2)*ngfft(3))
   * enddo
   */
  for (n = iset = 0; iset <= set; ++iset) {
    int const xsize = data->vol[iset].xsize; 
    int const ysize = data->vol[iset].ysize;
    int const zsize = data->vol[iset].zsize;

    char recordmarker[10];
    int ix, iy, iz;

    /* Fill the datablock with the density data.
     * Note that for each 'iset'-loop the datablock is overwritten,
     * so that only the last datablock readings remain.
     */
    for (n = iz = 0; iz < zsize; ++iz) {
      for (iy = 0; iy < ysize; ++iy) {
        for (ix = 0; ix < xsize; ++ix, ++n) {
          double value;

          /* The datablock grid is one voxel larger in each direction in order
	   * to smoothen the density at the cell boundaries; hence fill the extra
	   * voxel with the same value as the one at the beginning of that row.
	   */
          if (ix == xsize - 1) value = datablock[n - ix];
	  else if (iy == ysize - 1) value = datablock[n - iy*xsize];
	  else if (iz == zsize - 1) value = datablock[n - iz*ysize*xsize];
	  else if (data->hdr->cplex == 1) {
            /* get the volumetric data and convert */
            binread(&value, 8, data->file, data->hdr->bintype);
	    value *= HARTREE_TO_EV;
	  } else if (data->hdr->cplex == 2) {
            /* get volumetric data as a complex number and convert */
	    double a, b;
            binread(&a, 8, data->file, data->hdr->bintype);
            binread(&b, 8, data->file, data->hdr->bintype);
	    value = sqrt(a*a + b*b) * HARTREE_TO_EV;
	  } else return MOLFILE_ERROR;

          datablock[n] = value;
	}
      }
    }

    /* THIS HAS NOT BEEN VERIFIED YET: skip the recordmarker bytes twice */
    fread(recordmarker, 1, data->hdr->bintype.recordmarker, data->file);
    fread(recordmarker, 1, data->hdr->bintype.recordmarker, data->file);
  }

  DBGPRINT(stderr, "Exit POT_read_volumetric_data\n");
  return MOLFILE_SUCCESS;
}


static int WFK_read_volumetric_metadata(abinit_plugindata_t *data, int *nvolsets, molfile_volumetric_t **metadata)
{
  /* net yet implemented */
  DBGPRINT(stderr, "Enter/Exit WFK_read_volumetric_metadata\n");
  fprintf(stderr, "\n\nABINIT read) WARNING: loading WFK is NOT YET IMPLEMENTED!\n");
  return MOLFILE_ERROR;
}


static int WFK_read_volumetric_data(abinit_plugindata_t *data, int set, float *datablock, float *colorblock)
{
  /* net yet implemented */
  DBGPRINT(stderr, "Enter/Exit WFK_read_volumetric_data: NOT YET IMPLEMENTED!\n");
  fprintf(stderr, "\n\nABINIT read) WARNING: loading WFK is NOT YET IMPLEMENTED!\n");
  return MOLFILE_ERROR;
}

/* ===================================
 * Generic vmd routines
 * ===================================
 */


/* VMD calls this one just once to verify access to the file. */
static void *open_file_read(const char *filename, const char *filetype, int *natoms)
{
  void *result = NULL;
  abinit_plugindata_t *data;

  DBGPRINT(stderr, "Enter open_file_read\n");

  /* verify that the input variables are OK */
  if (!filename || !natoms) return NULL;

  /* start with undefined value and set it after successful read */
  *natoms = MOLFILE_NUMATOMS_UNKNOWN;

  /* allocate memory for the abinit data structure */
  data = abinit_plugindata_malloc();
  if (!data) return NULL;

  /* allocate memory for filename (add an extra 10 bytes for some more flexibility) */
  data->filename = (char *)malloc( sizeof(char) * (strlen(filename) + 10) );

  /* open the file for reading */
  data->file = fopen(filename, "rb");

  if (!data->file || !data->filename) {
    abinit_plugindata_free(data);
    return NULL;
  }
  strcpy(data->filename, filename);

  if (abinit_filetype(data, "GEO"))
    result = GEO_open_file_read(data, natoms);
  else if (abinit_filetype(data, "DEN") || abinit_filetype(data, "POT") || abinit_filetype(data, "WFK"))
    result = DEN_POT_WFK_open_file_read(data, natoms);

  if (result == NULL) abinit_plugindata_free(data);

  DBGPRINT(stderr, "Exit open_file_read\n");
  return result;
}


/* VMD calls this once to find out about the atom types in the structure */
static int read_structure(void *mydata, int *optflags, molfile_atom_t *atomlist)
{
  int result = MOLFILE_ERROR;
  abinit_plugindata_t *data = (abinit_plugindata_t *)mydata;

  DBGPRINT(stderr, "Enter read_structure\n");

  if (!data || !optflags || !atomlist) return MOLFILE_ERROR;

  if (abinit_filetype(data, "GEO"))
    result = GEO_read_structure(data, optflags, atomlist);
  else if (abinit_filetype(data, "DEN") || abinit_filetype(data, "POT") || abinit_filetype(data, "WFK"))
    result = DEN_POT_WFK_read_structure(data, optflags, atomlist);

  DBGPRINT(stderr, "Exit read_structure\n");
  return result;
}


/* VMD keeps calling for a next timestep, until it gets End-Of-File here */
static int read_next_timestep(void *mydata, int natoms, molfile_timestep_t *ts)
{
  int result = MOLFILE_EOF;
  abinit_plugindata_t *data = (abinit_plugindata_t *)mydata;

  DBGPRINT(stderr, "Enter read_next_timestep\n");

  /* Save coordinatess only if we are given a timestep pointer.
   * Otherwise assume that VMD wants us to skip past it.
   */
  if (!ts || !data) return MOLFILE_EOF;

  /* Double check that the number of atoms are correct */
  if (natoms != data->natom) return MOLFILE_EOF;

  if (abinit_filetype(data, "GEO"))
    result = GEO_read_next_timestep(data, natoms, ts);
  else if (abinit_filetype(data, "DEN") || abinit_filetype(data, "POT") || abinit_filetype(data, "WFK"))
    result = DEN_POT_WFK_read_next_timestep(data, natoms, ts);

  DBGPRINT(stderr, "Exit read_next_timestep\n");
  return result;
}


static void close_file_read(void *mydata)
{
  abinit_plugindata_t *data = (abinit_plugindata_t *)mydata;

  DBGPRINT(stderr, "Enter close_read\n");

  abinit_plugindata_free(data);

  DBGPRINT(stderr, "Exit close_read\n");
}



/* This writes the basic structure data in the syntax of the abinit input format. */
static void *open_file_write(const char *filename, const char *filetype, int natoms)
{
  abinit_plugindata_t *data = abinit_plugindata_malloc();

  DBGPRINT(stderr, "Enter open_file_write\n");

  if (!data) return NULL;

  /* allocate memory for filename (add an extra 10 bytes for some more flexibility) */
  data->filename = (char *)malloc( sizeof(char) * (strlen(filename) + 10) );

  /* open the file for writing */
  data->file = fopen(filename, "w");

  if (!data->filename || !data->file) {
    abinit_plugindata_free(data);
    fprintf(stderr, "ABINIT write) ERROR: unable to open file '%s' for writing\n", filename);
    return NULL;
  }
  strcpy(data->filename, filename);

  data->natom = natoms;

  DBGPRINT(stderr, "Exit open_file_write\n");
  return data;
}


static int write_structure(void *mydata, int optflags, const molfile_atom_t *atoms)
{
  abinit_plugindata_t *data = (abinit_plugindata_t *)mydata;
  int i, znucl[NATOM_MAX], ntypat;

  DBGPRINT(stderr, "Enter write_structure\n");

  if (!data || !atoms) return MOLFILE_ERROR;

  for (i = 0; i < NATOM_MAX; ++i) znucl[i] = 0;


  for (ntypat = i = 0; i < data->natom; ++i) {
    int const idx = get_pte_idx(atoms[i].type);
    int k;

    /* check if we already have this atom's idx in the list */
    for (k = 0; k < ntypat; ++k) if (idx == znucl[k]) break;

    /* if it is not yet in the list, increment ntypat */
    if (k == ntypat) ntypat++;

    znucl[k] = idx;
    data->typat[i] = k + 1;
  }

  /* write header with info */
  fprintf(data->file, "# Format below is in a sloppy ABINIT style.\n");
  fprintf(data->file, "# See http://www.abinit.org/ for the meaning of the keywords used here.\n\n");

  /* write the atom data to file */
  fprintf(data->file, "# Definition of the atom types\nntypat %d\nznucl ", ntypat);
  for (i = 0; i < ntypat; ++i) fprintf(data->file, " %d", znucl[i]);
  fprintf(data->file, "\n\n");

  fprintf(data->file, "# Definition of the atoms\nnatom %d\ntypat ", data->natom);
  for (i = 0; i < data->natom; ++i) fprintf(data->file, " %d", data->typat[i]);
  fprintf(data->file, "\n\n");

  DBGPRINT(stderr, "Exit write_structure\n");
  return MOLFILE_SUCCESS;
}


static int write_timestep(void *mydata, const molfile_timestep_t *ts)
{
  abinit_plugindata_t *data = (abinit_plugindata_t *)mydata; 
  int i;

  DBGPRINT(stderr, "Enter write_timestep\n");

  if (!data || !ts) return MOLFILE_ERROR;

  fprintf(data->file, "# Definition of the unit cell in Bohr\n");
  fprintf(data->file, "acell %f %f %f\n", ts->A * ANGS_TO_BOHR, ts->B * ANGS_TO_BOHR, ts->C * ANGS_TO_BOHR);
  fprintf(data->file, "angdeg %f %f %f\n\n", ts->alpha, ts->beta, ts->gamma);

  fprintf(data->file, "# location of the atoms in Bohr\nxcart ");
  for (i = 0; i < data->natom; ++i) {
    float const rx = ts->coords[3*i    ]  * ANGS_TO_BOHR,
                ry = ts->coords[3*i + 1]  * ANGS_TO_BOHR,
                rz = ts->coords[3*i + 2]  * ANGS_TO_BOHR;
    fprintf(data->file, "%s%17.12f %17.12f %17.12f\n", (i != 0 ? "      " : ""), rx, ry, rz);
  }
  fprintf(data->file, "\n\n");

  DBGPRINT(stderr, "Exit write_timestep\n");
  return MOLFILE_SUCCESS;
}


static void close_file_write(void *mydata)
{
  abinit_plugindata_t *data = (abinit_plugindata_t *)mydata;
  DBGPRINT(stderr, "Enter close_file_write\n");

  abinit_plugindata_free(data);

  DBGPRINT(stderr, "Exit close_file_write\n");
}


static int read_volumetric_metadata(void *mydata, int *nvolsets, molfile_volumetric_t **metadata)
{
  int result = MOLFILE_ERROR;
  abinit_plugindata_t *data = (abinit_plugindata_t *)mydata;

  DBGPRINT(stderr, "Enter read_volumetric_metadata\n");

  if (!data || !nvolsets || !metadata) return MOLFILE_ERROR;

  if (abinit_filetype(data, "DEN"))
    result = DEN_read_volumetric_metadata(data, nvolsets, metadata);
  else if (abinit_filetype(data, "POT"))
    result = POT_read_volumetric_metadata(data, nvolsets, metadata);
  else if (abinit_filetype(data, "WFK"))
    result = WFK_read_volumetric_metadata(data, nvolsets, metadata);

  DBGPRINT(stderr, "Exit read_volumetric_metadata\n");
  return result;
}

static int read_volumetric_data(void *mydata, int set, float *datablock, float *colorblock)
{
  int result = MOLFILE_ERROR;
  abinit_plugindata_t *data = (abinit_plugindata_t *)mydata;

  DBGPRINT(stderr, "Enter read_volumetric_data\n");

  if (!data || !datablock) return MOLFILE_ERROR;

  if (abinit_filetype(data, "DEN"))
    result = DEN_read_volumetric_data(data, set, datablock, colorblock);
  else if (abinit_filetype(data, "POT"))
    result = POT_read_volumetric_data(data, set, datablock, colorblock);
  else if (abinit_filetype(data, "WFK"))
    result = WFK_read_volumetric_data(data, set, datablock, colorblock);

  DBGPRINT(stderr, "Exit read_volumetric_data\n");
  return result;
}


/* ===================================
 * Registration stuff
 * ===================================
 */

static molfile_plugin_t abinitplugin;

VMDPLUGIN_API int VMDPLUGIN_init() {
  /* zero all bytes in the plugin prior to setting the fields */
  memset(&abinitplugin, 0, sizeof(molfile_plugin_t));

  /* see molfile_plugin.h for details on the fields */
  /* header  */
  abinitplugin.abiversion   = vmdplugin_ABIVERSION; /* the ABI for the base plugin type */
  abinitplugin.type         = MOLFILE_PLUGIN_TYPE;  /* string descriptor of the plugin type */
  abinitplugin.name         = "ABINIT";             /* name for the plugin */
  abinitplugin.prettyname   = "ABINIT";             /* name in filetype list */
  abinitplugin.author       = "Rob Lahaye";         /* string identifier */
  abinitplugin.majorv       = 0;                    /* major version */
  abinitplugin.minorv       = 4;                    /* minor version */
  abinitplugin.is_reentrant = VMDPLUGIN_THREADSAFE; /* can this library be run concurrently with itself? */

  /* the rest of the plugin */
  abinitplugin.filename_extension       = "*|*_GEO|*_DEN|*_WFK|*_POT|*_VHA|*_VHXC|*_VXC"; /* filename extension for this file type */
  abinitplugin.open_file_read           = open_file_read;           /* try to open the file for reading */
  abinitplugin.read_structure           = read_structure;           /* read molecular structure from the given file handle */
  abinitplugin.read_bonds               = 0;                        /* read bond information for the molecule */
  abinitplugin.read_next_timestep       = read_next_timestep;       /* read the next timestep from the file */
  abinitplugin.close_file_read          = close_file_read;          /* close the file and release all data */
  abinitplugin.open_file_write          = open_file_write;          /* open a coordinate file for writing */
  abinitplugin.write_structure          = write_structure;          /* write a timestep to the coordinate file */
  abinitplugin.write_timestep           = write_timestep;           /* write a timestep to the coordinate file */
  abinitplugin.close_file_write         = close_file_write;         /* close the file and release all data */
  abinitplugin.read_volumetric_metadata = read_volumetric_metadata; /* retrieve metadata pertaining to volumetric datasets in this file */
  abinitplugin.read_volumetric_data     = read_volumetric_data;     /* read the specified volumetric data set into the space pointed to by datablock */
  abinitplugin.read_rawgraphics         = 0;                        /* read raw graphics data stored in this file */
  abinitplugin.read_molecule_metadata   = 0;                        /* read molecule metadata */
  abinitplugin.write_bonds              = 0;                        /* write bond information for the molecule */

  return VMDPLUGIN_SUCCESS;
}

VMDPLUGIN_API int VMDPLUGIN_register(void *v, vmdplugin_register_cb cb) {

  (*cb)(v, (vmdplugin_t *)&abinitplugin);
  return VMDPLUGIN_SUCCESS;
}

VMDPLUGIN_API int VMDPLUGIN_fini() {
  return VMDPLUGIN_SUCCESS;
}


/* ===================================
 * Procedures to read the binary data
 * ===================================
 */


/* Do a binary read of 'size' bytes from stream.
 * Store the bytes into the pointer according to the match
 * between the endianness of the file and the system reading it.
 * Return the number of bytes succesfully read.
 */
static int binread(void *ptr, size_t size, FILE *stream, binary_t bintype)
{
  char *data = (char *)ptr;
  char *storage = malloc(size * sizeof(char));
  int const result = fread(storage, 1, size, stream);
  unsigned int i = 1;
  enum Endianness myEndian = ( *(char *)&i == 1 ? little_endian : big_endian );

  for (i = 0; i < size; ++i) {
    int const index = (myEndian == bintype.endian ? i : size - i - 1);
    data[i] = storage[index];
  }

  free(storage);
  return result;
}


/* Read the full header of a binary ABINIT file until the point
 * where the data starts.
 * Return a pointer to the header struct, or NULL if the reading fails.
 */ 
static abinit_binary_header_t *abinit_header(FILE *fp)
{
  int const debug = 0;

  /* the pattern of the first 4 or 8 bytes of a binary abinit file determines the
   * endianness and recordmarker length of the file.
   */
  char const Big_Endian_4_pattern[]    = {'\x00','\x00','\x00','\x0e'},
             Big_Endian_8_pattern[]    = {'\x00','\x00','\x00','\x00','\x00','\x00','\x00','\x0e'}, 
             Little_Endian_4_pattern[] = {'\x0e','\x00','\x00','\x00'},
             Little_Endian_8_pattern[] = {'\x0e','\x00','\x00','\x00','\x00','\x00','\x00','\x00'};

  char skip[1024];
  int i, bc = 0;

  abinit_binary_header_t *hdr = abinit_header_malloc();
  if (!hdr) return NULL;

  /* be sure to start from the beginning of the file */
  rewind(fp);

  if (debug) fprintf(stderr, "START OF BINARY FILE DEBUG INFO\n");

  /* Determine endianness and record length of the binary file.
   * Note that this pattern is followed by a version string, which has no null termination.
   * Therefore the pattern is always followed by a non-null byte.
   */
  fread(skip, 1, 8, fp);
  if      (memcmp(skip, Big_Endian_4_pattern, 4) == 0)    {hdr->bintype.endian = big_endian;    hdr->bintype.recordmarker = 4;}
  else if (memcmp(skip, Big_Endian_8_pattern, 8) == 0)    {hdr->bintype.endian = big_endian;    hdr->bintype.recordmarker = 8;}
  else if (memcmp(skip, Little_Endian_8_pattern, 8) == 0) {hdr->bintype.endian = little_endian; hdr->bintype.recordmarker = 8;}
  else if (memcmp(skip, Little_Endian_4_pattern, 4) == 0) {hdr->bintype.endian = little_endian; hdr->bintype.recordmarker = 4;}
  else {abinit_header_free(hdr); return NULL;}

  if (debug) fprintf(stderr, "Binary file is in %s with a record-marker of %d bytes\n",
	                       (hdr->bintype.endian == big_endian ? "Big Endian" : "Little Endian"), hdr->bintype.recordmarker);

  /* start reading again from the beginning of the file */
  rewind(fp);

  /*The header
   * The wavefunction files, density files, and potential files all begin with the same records,
   * called the "header".
   *
   *character*6 :: codvsn
   *integer :: headform,fform
   *integer :: bantot,date,intxc,ixc,natom,ngfft(3),nkpt,npsp,
   *nspden,nspinor,nsppol,nsym,ntypat,occopt,pertcase,usepaw
   *integer :: usewvl, cplex, nspden
   *double precision :: acell(3),ecut,ecutdg,ecutsm,ecut_eff,qptn(3),rprimd(3,3),stmbias,tphysel,tsmear
   *integer :: istwfk(nkpt),nband(nkpt*nsppol),npwarr(nkpt),so_psp(npsp),&
   *& symafm(nsym),symrel(3,3,nsym),typat(natom),nrhoijsel(nspden),rhoijselect(*,nspden)
   *double precision :: kpt(3,nkpt),occ(bantot),tnons(3,nsym),znucltypat(ntypat),wtk(nkpt)
   *character*132 :: title
   *double precision :: znuclpsp,zionpsp
   *integer :: pspso,pspdat,pspcod,pspxc,lmax,lloc,mmax=integers
   *double precision :: residm,xred(3,natom),etotal,fermie,rhoij(*,nspden)
   *
   * write(unit=header) codvsn,headform,fform
   * write(unit=header) bantot,date,intxc,ixc,natom,ngfft(1:3),&
   *& nkpt,nspden,nspinor,nsppol,nsym,npsp,ntypat,occopt,pertcase,usepaw,&
   *& ecut,ecutdg,ecutsm,ecut_eff,qptn(1:3),rprimd(1:3,1:3),stmbias,tphysel,tsmear,usewvl

   * write(unit=header) istwfk(1:nkpt),nband(1:nkpt*nsppol),&
   *& npwarr(1:nkpt),so_psp(1:npsp),symafm(1:nsym),symrel(1:3,1:3,1:nsym),typat(1:natom),&
   *& kpt(1:3,1:nkpt),occ(1:bantot),tnons(1:3,1:nsym),znucltypat(1:ntypat),wtk(1:nkpt)
   * do ipsp=1,npsp
   *! (npsp lines, 1 for each pseudopotential ; npsp=ntypat, except if alchemical pseudo-atoms)
   *  write(unit=unit) title,znuclpsp,zionpsp,pspso,pspdat,pspcod,pspxc,lmn_size
   * enddo
   *!(in case of usepaw==0, final record: residm, coordinates, total energy, Fermi energy)
   * write(unit=unit) residm,xred(1:3,1:natom),etotal,fermie
   *!(in case of usepaw==1, there are some additional records)
   * if (usepaw==1)then
   *  write(unit=unit)( pawrhoij(iatom)%nrhoijsel(1:nspden),iatom=1,natom), cplex, nspden
   *  write(unit=unit)((pawrhoij(iatom)%rhoijselect(1:      nrhoijsel(ispden),ispden),ispden=1,nspden),iatom=1,natom),&
   *&                 ((pawrhoij(iatom)%rhoijp     (1:cplex*nrhoijsel(ispden),ispden),ispden=1,nspden),iatom=1,natom)
   * endif
   */

  /* skip first the recordmarker bytes */
  bc += fread(skip, 1, hdr->bintype.recordmarker, fp);

  /* code version */
  bc += fread(hdr->codvsn, sizeof(char), 6, fp); hdr->codvsn[6] = '\0';
  if (debug) fprintf(stderr, "codvsn = '%s' (code version)\n", hdr->codvsn); 

  /* format of the header */
  bc += binread(&hdr->headform, 4, fp, hdr->bintype);
  if (debug) fprintf(stderr, "headform = '%d' (format of the header)\n", hdr->headform);

  /* specification for data type: 2 for wf; 52 for density; 102 for potential */
  bc += binread(&hdr->fform, 4, fp, hdr->bintype);
  if (debug) fprintf(stderr, "fform = '%d' (2 for wf; 52 for density; 102 for potential)\n", hdr->fform);

  /* skip the recordmarker bytes two times */
  bc += fread(skip, 1, 2 * hdr->bintype.recordmarker, fp);

  /* total number of bands (sum of nband on all kpts and spins) */
  bc += binread(&hdr->bantot, 4, fp, hdr->bintype);
  if (debug) fprintf(stderr, "bantot = '%d' (sum of nband on all kpts and spins)\n", hdr->bantot);

  /* starting date */
  bc += binread(&hdr->date, 4, fp, hdr->bintype);
  if (debug) fprintf(stderr, "date = '%d' (starting date)\n", hdr->date);

  bc += binread(&hdr->intxc, 4, fp, hdr->bintype);
  if (debug) fprintf(stderr, "intxc = '%d'\n", hdr->intxc);

  bc += binread(&hdr->ixc, 4, fp, hdr->bintype);
  if (debug) fprintf(stderr, "ixc = '%d'\n", hdr->ixc);

  bc += binread(&hdr->natom, 4, fp, hdr->bintype);
  if (debug) fprintf(stderr, "natom = '%d'\n", hdr->natom);

  /* failsafe */
  if (hdr->natom <= 0) {
    fprintf(stderr, "ABINIT read) ERROR Binary Header: natom = %d is wrong!",  hdr->natom);
    abinit_header_free(hdr);
    return NULL;
  }

  bc += binread(&hdr->ngfft[0], 4, fp, hdr->bintype);
  if (debug) fprintf(stderr, "ngfft[0] = '%d'\n", hdr->ngfft[0]);
  bc += binread(&hdr->ngfft[1], 4, fp, hdr->bintype);
  if (debug) fprintf(stderr, "ngfft[1] = '%d'\n", hdr->ngfft[1]);
  bc += binread(&hdr->ngfft[2], 4, fp, hdr->bintype);
  if (debug) fprintf(stderr, "ngfft[2] = '%d'\n", hdr->ngfft[2]);

  bc += binread(&hdr->nkpt, 4, fp, hdr->bintype);
  if (debug) fprintf(stderr, "nkpt = '%d'\n", hdr->nkpt);

  bc += binread(&hdr->nspden, 4, fp, hdr->bintype);
  if (debug) fprintf(stderr, "nspden = '%d'\n", hdr->nspden);

  /* failsafe */
  if (hdr->nspden != 1 && hdr->nspden != 2 && hdr->nspden != 4) {
    fprintf(stderr, "ABINIT read) ERROR Binary Header: nspden = %d is wrong!",  hdr->nspden);
    abinit_header_free(hdr);
    return NULL;
  }

  bc += binread(&hdr->nspinor, 4, fp, hdr->bintype);
  if (debug) fprintf(stderr, "nspinor = '%d'\n", hdr->nspinor);

  bc += binread(&hdr->nsppol, 4, fp, hdr->bintype);
  if (debug) fprintf(stderr, "nsppol = '%d'\n", hdr->nsppol);

  bc += binread(&hdr->nsym, 4, fp, hdr->bintype);
  if (debug) fprintf(stderr, "nsym = '%d'\n", hdr->nsym);

  bc += binread(&hdr->npsp, 4, fp, hdr->bintype);
  if (debug) fprintf(stderr, "npsp = '%d'\n", hdr->npsp);

  bc += binread(&hdr->ntypat, 4, fp, hdr->bintype);
  if (debug) fprintf(stderr, "ntypat = '%d'\n", hdr->ntypat);

  bc += binread(&hdr->occopt, 4, fp, hdr->bintype);
  if (debug) fprintf(stderr, "occopt = '%d'\n", hdr->occopt);

  /* the index of the perturbation, 0 if GS calculation */
  bc += binread(&hdr->pertcase, 4, fp, hdr->bintype);
  if (debug) fprintf(stderr, "pertcase = '%d' (the index of the perturbation, 0 if GS calculation)\n", hdr->pertcase);

  bc += binread(&hdr->usepaw, 4, fp, hdr->bintype);
  if (debug) fprintf(stderr, "usepaw = '%d' (0=norm-conserving psps, 1=paw)\n", hdr->usepaw);

  /* failsafe */
  if (hdr->usepaw != 0 && hdr->usepaw != 1) {
    fprintf(stderr, "ABINIT read) ERROR Binary Header: usepaw = %d is wrong!", hdr->usepaw);
    abinit_header_free(hdr);
    return NULL;
  }

  bc += binread(&hdr->ecut, 8, fp, hdr->bintype);
  if (debug) fprintf(stderr, "ecut = '%g'\n", hdr->ecut);

  /* input variable (ecut for NC psps, pawecutdg for paw) */
  bc += binread(&hdr->ecutdg, 8, fp, hdr->bintype);
  if (debug) fprintf(stderr, "ecutdg = '%g' (ecut for NC psps, pawecutdg for paw)\n", hdr->ecutdg);

  bc += binread(&hdr->ecutsm, 8, fp, hdr->bintype);
  if (debug) fprintf(stderr, "ecutsm = '%g'\n", hdr->ecutsm);

  /* ecut*dilatmx**2 (dilatmx is an input variable) */
  bc += binread(&hdr->ecut_eff, 8, fp, hdr->bintype);
  if (debug) fprintf(stderr, "ecut_eff = '%g' (ecut*dilatmx**2 [dilatmx is an input variable])\n", hdr->ecut_eff);

  bc += binread(&hdr->qptn[0], 8, fp, hdr->bintype);
  if (debug) fprintf(stderr, "qptn[0] = '%g'\n", hdr->qptn[0]);
  bc += binread(&hdr->qptn[1], 8, fp, hdr->bintype);
  if (debug) fprintf(stderr, "qptn[1] = '%g'\n", hdr->qptn[1]);
  bc += binread(&hdr->qptn[2], 8, fp, hdr->bintype);
  if (debug) fprintf(stderr, "qptn[2] = '%g'\n", hdr->qptn[2]);

  for (i = 0; i < 3; ++i) {
    bc += binread(&hdr->rprimd[i][0], 8, fp, hdr->bintype);
    if (debug) fprintf(stderr, "rprimd[%d][0] = '%g'\n", i, hdr->rprimd[i][0]);
    bc += binread(&hdr->rprimd[i][1], 8, fp, hdr->bintype);
    if (debug) fprintf(stderr, "rprimd[%d][1] = '%g'\n", i, hdr->rprimd[i][1]);
    bc += binread(&hdr->rprimd[i][2], 8, fp, hdr->bintype);
    if (debug) fprintf(stderr, "rprimd[%d][2] = '%g'\n", i, hdr->rprimd[i][2]);
  }

  bc += binread(&hdr->stmbias, 8, fp, hdr->bintype);
  if (debug) fprintf(stderr, "stmbias = '%g'\n", hdr->stmbias);

  bc += binread(&hdr->tphysel, 8, fp, hdr->bintype);
  if (debug) fprintf(stderr, "tphysel = '%g'\n", hdr->tphysel);

  bc += binread(&hdr->tsmear, 8, fp, hdr->bintype);
  if (debug) fprintf(stderr, "tsmear = '%g'\n", hdr->tsmear);

  bc += binread(&hdr->usewvl, 4, fp, hdr->bintype);
  if (debug) fprintf(stderr, "usewvl = '%d'\n", hdr->usewvl);

  /* failsafe */
  if (hdr->usewvl != 0 && hdr->usewvl != 1) {
    fprintf(stderr, "ABINIT read) ERROR Binary Header: usewvl = %d is wrong!",  hdr->usewvl);
    abinit_header_free(hdr);
    return NULL;
  }

  /* skip the recordmarker bytes two times */
  bc += fread(skip, 1, 2 * hdr->bintype.recordmarker, fp);

  hdr->istwfk = (int *)malloc(sizeof(int) * hdr->nkpt);
  hdr->nband  = (int *)malloc(sizeof(int) * hdr->nkpt * hdr->nsppol);
  hdr->npwarr = (int *)malloc(sizeof(int) * hdr->nkpt);
  hdr->so_psp = (int *)malloc(sizeof(int) * hdr->npsp);
  hdr->symafm = (int *)malloc(sizeof(int) * hdr->nsym);
  hdr->typat  = (int *)malloc(sizeof(int) * hdr->natom);

  if (!hdr->istwfk || !hdr->nband || !hdr->npwarr || !hdr->so_psp || !hdr->symafm || !hdr->typat) {
    abinit_header_free(hdr);
    return NULL;
  }
  for (i = 0; i < 3; ++i) {
    int j;
    for (j = 0; j < 3; ++j) {
      hdr->symrel[i][j] = (int *)malloc(sizeof(int) * hdr->nsym);
      if (!hdr->symrel[i][j]) {abinit_header_free(hdr); return NULL;}
    }
  }

  for (i = 0; i < hdr->nkpt; ++i) {
    bc += binread(&hdr->istwfk[i], 4, fp, hdr->bintype);
    if (debug) fprintf(stderr, "istwfk[%d] = '%d'\n", i, hdr->istwfk[i]);
  }

  for (i = 0; i < hdr->nkpt * hdr->nsppol; ++i) {
    bc += binread(&hdr->nband[i], 4, fp, hdr->bintype);
    if (debug) fprintf(stderr, "nband[%d] = '%d'\n", i, hdr->nband[i]);
  }

  for (i = 0; i < hdr->nkpt; ++i) {
    bc += binread(&hdr->npwarr[i], 4, fp, hdr->bintype);
    if (debug) fprintf(stderr, "npwarr[%d] = '%d'\n", i, hdr->npwarr[i]);
  }

  for (i = 0; i < hdr->npsp; ++i) {
    bc += binread(&hdr->so_psp[i], 4, fp, hdr->bintype);
    if (debug) fprintf(stderr, "so_psp[%d] = '%d'\n", i, hdr->so_psp[i]);
  }

  for (i = 0; i < hdr->nsym; ++i) {
    bc += binread(&hdr->symafm[i], 4, fp, hdr->bintype);
    if (debug) fprintf(stderr, "symafm[%d] = '%d'\n", i, hdr->symafm[i]);
  }

  for (i = 0; i < hdr->nsym; ++i) {
    int j;
    for (j = 0; j < 3; ++j) {
      int k;     
      for (k = 0; k < 3; ++k) {
        bc += binread(&hdr->symrel[k][j][i], 4, fp, hdr->bintype);
        if (debug) fprintf(stderr, "symrel[%d][%d][%2d]= '%2d'", k, j, i, hdr->symrel[k][j][i]);
      }
      if (debug) fprintf(stderr, "\n");
    }
    if (debug) fprintf(stderr, "\n");
  }

  for (i = 0; i < hdr->natom; ++i) {
    bc += binread(&hdr->typat[i], 4, fp, hdr->bintype);
    if (debug) fprintf(stderr, "typat[%d] = '%d'\n", i, hdr->typat[i]);
  }

  for (i = 0; i < 3; ++i) {
    hdr->kpt[i] = (double *)malloc(sizeof(double) * hdr->nkpt);
    hdr->tnons[i] = (double *)malloc(sizeof(double) * hdr->nsym);
    if (!hdr->kpt[i] || !hdr->tnons[i]) {abinit_header_free(hdr); return NULL;}
  }

  hdr->occ = (double *)malloc(sizeof(double) * hdr->bantot);
  hdr->znucltypat = (double *)malloc(sizeof(double) * hdr->ntypat);
  hdr->wtk = (double *)malloc(sizeof(double) * hdr->nkpt);
  if (!hdr->occ || !hdr->znucltypat || !hdr->wtk) {abinit_header_free(hdr); return NULL;}

  for (i = 0; i < hdr->nkpt; ++i) {
    int j;
    for (j = 0; j < 3; ++j) {
      bc += binread(&hdr->kpt[j][i], 8, fp, hdr->bintype);
      if (debug) fprintf(stderr, "kpt[%d][%2d] = '%g'\n", j, i, hdr->kpt[j][i]);
    }
  }

  for (i = 0; i < hdr->bantot; ++i) {
    bc += binread(&hdr->occ[i], 8, fp, hdr->bintype);
    if (debug) fprintf(stderr, "occ[%d] = '%g'\n", i, hdr->occ[i]);
  }

  for (i = 0; i < hdr->nsym; ++i) {
    int j;
    for (j = 0; j < 3; ++j) {
      bc += binread(&hdr->tnons[j][i], 8, fp, hdr->bintype);
      if (debug) fprintf(stderr, "tnons[%d][%2d] = '%g' ", j, i, hdr->tnons[j][i]);
    }
    if (debug) fprintf(stderr, "\n");
  }

  for (i = 0; i < hdr->ntypat; ++i) {
    bc += binread(&hdr->znucltypat[i], 8, fp, hdr->bintype);
    if (debug) fprintf(stderr, "znucltypat[%d] = '%g'\n", i, hdr->znucltypat[i]);
  }

  for (i = 0; i < hdr->nkpt; ++i) {
    bc += binread(&hdr->wtk[i], 8, fp, hdr->bintype);
    if (debug) fprintf(stderr, "wtk[%d] = '%g'\n", i, hdr->wtk[i]);
  }

  for (i = 0; i < hdr->npsp; ++i) {

    /* skip the recordmarker bytes two times */
    bc += fread(skip, 1, 2 * hdr->bintype.recordmarker, fp);

    bc += fread(hdr->title, sizeof(char), 132, fp); hdr->title[132] = '\0';
    if (debug) fprintf(stderr, "title = '%s'\n", hdr->title);

    bc += binread(&hdr->znuclpsp, 8, fp, hdr->bintype);
    if (debug) fprintf(stderr, "znuclpsp = '%g'\n", hdr->znuclpsp);

    bc += binread(&hdr->zionpsp, 8, fp, hdr->bintype);
    if (debug) fprintf(stderr, "zionpsp = '%g'\n", hdr->zionpsp);

    bc += binread(&hdr->pspso, 4, fp, hdr->bintype);
    if (debug) fprintf(stderr, "pspso = '%d'\n", hdr->pspso);

    bc += binread(&hdr->pspdat, 4, fp, hdr->bintype);
    if (debug) fprintf(stderr, "pspdat = '%d'\n", hdr->pspdat);

    bc += binread(&hdr->pspcod, 4, fp, hdr->bintype);
    if (debug) fprintf(stderr, "pspcod = '%d'\n", hdr->pspcod);

    bc += binread(&hdr->pspxc, 4, fp, hdr->bintype);
    if (debug) fprintf(stderr, "pspxc = '%d'\n", hdr->pspxc);

    bc += binread(&hdr->lmn_size, 4, fp, hdr->bintype);
    if (debug) fprintf(stderr, "lmn_size = '%d'\n", hdr->lmn_size);
  }

  /* skip the recordmarker bytes two times */
  bc += fread(skip, 1, 2 * hdr->bintype.recordmarker, fp);

  for (i = 0; i < 3; ++i) {
    hdr->xred[i] = (double *)malloc(sizeof(double) * hdr->natom);
    if (!hdr->xred[i]) {abinit_header_free(hdr); return NULL;}
  }

  bc += binread(&hdr->residm, 8, fp, hdr->bintype);
  if (debug) fprintf(stderr, "residm = '%g'\n", hdr->residm);

  for (i = 0; i < hdr->natom; ++i) {
    int j;
    for (j = 0; j < 3; ++j) {
      bc += binread(&hdr->xred[j][i], 8, fp, hdr->bintype);
      if (debug) fprintf(stderr, "xred[%d][%d] = '%g'\n", j, i, hdr->xred[j][i]);

      /* failsafe */
      if (hdr->xred[j][i] < -1 || hdr->xred[j][i] > 1) {
        fprintf(stderr, "Binary Header Error: hdr->xred[%d][%d] = %g; something must be wrong!", j, i, hdr->xred[j][i]);
        {abinit_header_free(hdr); return NULL;}
      }
    }
  }

  bc += binread(&hdr->etotal, 8, fp, hdr->bintype);
  if (debug) fprintf(stderr, "etotal = '%g'\n", hdr->etotal);

  bc += binread(&hdr->fermie, 8, fp, hdr->bintype);
  if (debug) fprintf(stderr, "fermie = '%g'\n", hdr->fermie);

  if (hdr->usepaw == 1) {
    struct pawrhoij_t {
      int *nrhoijsel;
      int **rhoijselect;
      double **rhoij;
    };
    struct pawrhoij_t *pawrhoij = malloc(hdr->natom * sizeof(struct pawrhoij_t));

    for (i = 0; i < hdr->natom; ++i) {
      int j;
      pawrhoij[i].nrhoijsel = (int *)malloc(sizeof(int) * hdr->nspden);
      if (!pawrhoij[i].nrhoijsel) {abinit_header_free(hdr); return NULL;}
      for (j = 0; j < hdr->nspden; ++j){
        bc += binread(&pawrhoij[i].nrhoijsel[j], 4, fp, hdr->bintype);
        if (debug) fprintf(stderr, "pawrhoij[%d].nrhoijsel[%d] = '%d'\n", i, j, pawrhoij[i].nrhoijsel[j]);
      }
      pawrhoij[i].rhoijselect = (int **)malloc(sizeof(int) * hdr->nspden);
      pawrhoij[i].rhoij = (double **)malloc(sizeof(double) * hdr->nspden);
      if (!pawrhoij[i].rhoijselect || !pawrhoij[i].rhoij) {abinit_header_free(hdr); return NULL;}
      for (j = 0; j < hdr->nspden; ++j) {
        int k;
        pawrhoij[i].rhoijselect[j] = (int *)malloc(sizeof(int) * pawrhoij[i].nrhoijsel[j]);
        pawrhoij[i].rhoij[j] = (double *)malloc(sizeof(double) * pawrhoij[i].nrhoijsel[j]);
        if (!pawrhoij[i].rhoijselect[j]) {abinit_header_free(hdr); return NULL;}
        for (k = 0; k < pawrhoij[i].nrhoijsel[j]; ++k) {
          bc += binread(&pawrhoij[i].rhoijselect[j][k], 4, fp, hdr->bintype);
	 if (debug) fprintf(stderr, "pawrhoij[%d].rhoijselect[%d][%d] = '%d'\n", i, j, k, pawrhoij[i].rhoijselect[j][k]);
        }
        for (k = 0; k < pawrhoij[i].nrhoijsel[j]; ++k) {
          bc += binread(&pawrhoij[i].rhoij[j][k], 8, fp, hdr->bintype);
	 if (debug) fprintf(stderr, "pawrhoij[%d].rhoij[%d][%d] = '%g'\n", i, j, k, pawrhoij[i].rhoij[j][k]);
        }
      }
      for (j = 0; j < hdr->nspden; ++j) {
        free(pawrhoij[i].rhoijselect[j]);
        free(pawrhoij[i].rhoij[j]);
      }
      free(pawrhoij[i].rhoijselect);
      free(pawrhoij[i].rhoij);
      free(pawrhoij[i].nrhoijsel);
    }
    free(pawrhoij);
  } /* end of "if (usepaw == 1)" */

  /* cplex depends on other variables:
   * In GS calculations, cplex = 1.
   * In response function calculations (non-zero pertcase), cplex = 1 at the gamma point (qpt=(0,0,0))
   * and cplex = 2 if qpt/=(0,0,0).
   */
  hdr->cplex = 1;
  if (hdr->pertcase != 0 && fabs(hdr->qptn[0]) > 10e-6 && fabs(hdr->qptn[1]) > 10e-6 && fabs(hdr->qptn[2]) > 10e-6) hdr->cplex = 2;


  /* skip the recordmarker bytes two times */
  bc += fread(skip, 1, 2 * hdr->bintype.recordmarker, fp);

  if (debug) fprintf(stderr, "END OF BINARY FILE DEBUG INFO\n");

  return hdr;
}
