/* MACHINE GENERATED FILE, DO NOT EDIT! */

#define VMDPLUGIN molfile_graspplugin
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
 *      $RCSfile: graspplugin.C,v $
 *      $Author: johns $       $:  $             $State: Exp $
 *      $Revision: 1.21 $       $Date: 2009/04/29 15:45:30 $
 *
 ***************************************************************************/

/* 
 * Reader for GRASP binary surface files
 * The file format is briefly described at these two sites:
 *   http://honiglab.cpmc.columbia.edu/grasp/grasp_contents.html#A.1
 *   http://www.msg.ucsf.edu/local/programs/grasp/html/Appendix%20A.html
 */

//Modified by Biol. Angel H. Jiménez Pardo and Luis Rosales Leòn
//Visualisation Department
//DGSCA
//Universidad Nacional Autònoma de Mèxico UNAM
//2006

#include <stdio.h>
#include <math.h>
#include <string.h>
#include "molfile_plugin.h"
#include "endianswap.h"

/// Define variables by the name of Grasp property  

#define POTENTIALS    1
#define CURVATURES    2
#define DISTANCES     4
#define PROPERTY_1    8
#define PROPERTY_2   16
#define VERTEXCOLOR  32
#define GLOBALCOLOR  64

/// Define other variables
typedef float COLOUR [3]; 
typedef float COLOUR [3];
typedef unsigned int SHORTWORD;
typedef float VECTOR [3];

/// Grasp structure
typedef struct {
  SHORTWORD flag, nvert;
  COLOUR clo, cmd, chi, colors;
  VECTOR range;
} GRASSP;

/// plugin data handle
typedef struct {
  FILE *fd;
  molfile_graphics_t *graphics;
} grasp_t;


/// Sets colour
void Set_Colour(float *c, float r, float g, float b) {
  *(c+0)=r;
  *(c+1)=g;
  *(c+2)=b;
} //Set_Colour end


/// Combines colour
void ClinComb2(float c[], float k1, float c1[], float k2, float c2[]) {
  (c)[0]=(k1)*(c1)[0] + (k2)*(c2)[0];
  (c)[1]=(k1)*(c1)[1] + (k2)*(c2)[1];
  (c)[2]=(k1)*(c1)[2] + (k2)*(c2)[2];
} //ClinComb2 end


/// Gets properties
void Get_Property_Values(GRASSP *grassp, float *properties, float *colores, int nvert) {
  long  i;
  const char *name[] = {"potential","curvature","distance","property1","property2"};
  SHORTWORD  index;
  int k=0, j=0;
  float  val, weight, min, mid, max, midmin, maxmid;
  val = weight = min = mid = max = midmin = maxmid = 0.0;

  //Range values
  grassp->range[0]=-1000.0;
  grassp->range[1]= 0.0;
  grassp->range[2]= 1000.0;

  ///Checks
  index = (SHORTWORD)( (log((double) grassp->flag) / log(2.0)) + 0.5 );

  ///ojo aqui le quite el  !
  if ((grassp->flag)!=POTENTIALS) {
    if (index >=0 && index <= 4) 
      printf("graspplugin) No data available for '%s' option\n", name[index]);
    else 
      printf("graspplugin) out of range property, flag: %d index: %d\n", grassp->flag, index);
    printf("graspplugin) Will use white color instead\n");
    grassp->flag = GLOBALCOLOR;
    Set_Colour(grassp->clo, 1, 1, 1);
  } else {
    printf("graspplugin) Getting %s values.\n", name[index]);
  }

  /// init max values 
  max=  0.01;  /* should be > 0 */
  min= -0.01;  /* should be < 0 */

  // Get values
  for(i=0; i < nvert; i++) {
    if (properties[i] < min)
      min= properties[i];
    else if(properties[i] > max)
      max= properties[i];
  }

  // cut properties that are out of specified range
  if (min < grassp->range[0] || max > grassp->range[2]) {
    for(i=0; i < nvert; i++) {
      if(properties[i] < grassp->range[0])
        properties[i]= grassp->range[0];
      else if(properties[i] > grassp->range[2])
        properties[i]= grassp->range[2];
    }
  } else {
    // or reset range
    grassp->range[0] = min;
    grassp->range[2] = max;
  }

  // check mid value 
  if (grassp->range[1] <= grassp->range[0] || 
      grassp->range[1] >= grassp->range[2])
   grassp->range[1] = (grassp->range[0] + grassp->range[2]) / 2;

  printf("graspplugin) Computing colors for range %g,%g,%g\n", 
         grassp->range[0], grassp->range[1], grassp->range[2]);

  // Prepare color interpolation parameters
  min = grassp->range[0];
  mid = grassp->range[1];
  max = grassp->range[2];
  midmin = mid-min;
  maxmid = max-mid;

  // Create color for each vertex and copies to a vector
  k=0;
  for (i=0; i < nvert; i++) {
    val = properties[i];
    if (val <= mid) {
      weight = (midmin) ? (val-min)/midmin : 0;  
      ClinComb2(grassp->colors, 1-weight, grassp->clo, weight, grassp->cmd);
      for(j=0; j<=2; j++) {
        *(colores+k)=grassp->colors[j];
        k++;
      }
    } else {
      weight = (maxmid) ? (val-mid)/maxmid : 0;
      ClinComb2(grassp->colors, 1-weight, grassp->cmd, weight, grassp->chi);
      for (j=0; j<=2; j++) {
        *(colores+k)=grassp->colors[j];
        k++;
      }       
    }
    // clean up 
  }
} // Get_Property_Values end


/// Reads line3
void line3 (FILE * infile, GRASSP * grassp) {
  char line3 [81];
  fread(line3, 1, 80, infile);

  grassp->flag=0;
  int i=0;
  if (line3 [0]==',') i++;
  while (i<80 && line3[i]!=' ') {
#if 0
    // XXX the rest of the code doesn't properly process flags yet, so
    // there's no point in setting them currently.
    if(!strncmp(line3+i, "potentials", 10))   grassp->flag |= POTENTIALS;
    if(!strncmp(line3+i, "curvature",   9))   grassp->flag |= CURVATURES;
    if(!strncmp(line3+i, "distances",   9))   grassp->flag |= DISTANCES;
    if(!strncmp(line3+i, "gproperty",   9))   grassp->flag |= PROPERTY_1;
    if(!strncmp(line3+i, "g2property", 10))   grassp->flag |= PROPERTY_2;
    if(!strncmp(line3+i, "vertexcolor",11))   grassp->flag |= VERTEXCOLOR;
#endif
    i++;
  }

  // Assign default property colors
  if (grassp->flag > 0 && grassp->flag < VERTEXCOLOR) {
    switch (grassp->flag ) {
      case POTENTIALS:
        Set_Colour(grassp->clo, 1.0, 0.0, 0.0 );
        Set_Colour(grassp->cmd, 1.0, 1.0, 1.0 );
        Set_Colour(grassp->chi, 0.0, 0.0, 1.0 );
        break;

      case CURVATURES:
        Set_Colour(grassp->clo, 0.5, 0.5, 0.5 );
        Set_Colour(grassp->cmd, 1.0, 1.0, 1.0 );
        Set_Colour(grassp->chi, 0.0, 1.0, 0.0 );
        break;

     case DISTANCES:
       Set_Colour(grassp->clo, 1.0, 1.0, 1.0 );
       Set_Colour(grassp->cmd, 0.0, 0.0, 1.0 );
       Set_Colour(grassp->chi, 1.0, 0.0, 0.0 );
       break;

     default: // Global color
       Set_Colour(grassp->clo, 1.0, 0.0, 0.0 );
       Set_Colour(grassp->cmd, 0.5, 0.0, 0.5 );
       Set_Colour(grassp->chi, 0.0, 0.0, 1.0 );
       break;
   }
 }

 if (!grassp->flag)
   grassp->flag = GLOBALCOLOR; 
} //line3 end


// check endianness
static int is_little_endian(void) {
  int x=1;
  return *((char *)&x);
}   


static void *open_file_read(const char *filepath, const char *filetype,
    int *natoms) {
  FILE *fd;
  grasp_t *grasp;
  
  fd = fopen(filepath, "rb");
  if (!fd) 
    return NULL;
  grasp = new grasp_t;
  grasp->fd = fd;
  grasp->graphics = NULL;
  *natoms = 0;
  return grasp;
}


static int read_rawgraphics(void *v, int *nelem, 
    const molfile_graphics_t **data) {
  grasp_t *grasp = (grasp_t *)v;
  FILE *infile = grasp->fd;

  // Reverse engineering is your friend, and combined with FORTRAN code, voila!
  // od -c shows the header starts off:
  // \0  \0  \0   P   f   o   r   m   a   t   =   2
  // and according to ungrasp, this is a 1 for grasp versions 1.0
  // and 1.1, and 2 for grasp version 1.2
  // Also, the header lines are of length 80 characters + 4 header chars
  // + 4 trailer chars
  // The 4 bytes at the beginning/end are standard Fortran array trash

  /// Pointers grassp type
  GRASSP datax;
   
  char trash[4];
#define TRASH fread(trash, 4, 1, infile)
  char line[81];

  // FIRST LINE OF HEADER; contains format type
  TRASH; 
  fread(line, 1, 80, infile); 
  // make sure it says 'format='
  if (strncmp(line, "format=", 7) != 0) {
    printf("graspplugin) First characters of file don't look like a GRASP file\n");
    return MOLFILE_ERROR;
  }
  TRASH;

  // next char should be a 0 or 1
  char gfiletype = line[7];
  if (gfiletype == '1') {
    gfiletype = 1;
  } else if (gfiletype == '2') {
    gfiletype = 2;
  } else {
    printf("graspplugin) GRASP file is in format %c, but only '1' or '2' is supported\n", gfiletype);
    return MOLFILE_ERROR;
  }

  // SECOND LINE: contains "vertices,accessibles,normals,triangles"
  TRASH; 
  fread(line, 1, 80, infile); 
  TRASH;

  // THIRD LINE: contains (0 or more of)?
  //  "potentials,curvature,distances,gproperty,g2property,vertexcolor
  TRASH; 
  line3(infile, &datax);/// Reads line 3
  TRASH;

  // FOURTH LINE stores vertices, triangles, gridsize, lattice spacing
  int nvert, ntriangles, gridsize;
  float lattice;
  TRASH; 
  fread(line, 1, 80, infile); 
  TRASH;
  sscanf(line, "%d%d%d%f", &nvert, &ntriangles, &gridsize, &lattice);

  /// Stores color
  float *colores = new float [3*nvert];

  // FIFTH LINE stores the center (x,y,z) position
  float center[3];
  TRASH; 
  fread(line, 1, 80, infile); 
  TRASH;
  sscanf(line, "%f%f%f", center, center+1, center+2);

  float *vertex = new float[3 * nvert];
  float *access = new float[3 * nvert];
  float *normal = new float[3 * nvert];
  int *triangle = new int[3 * ntriangles];
  float *properties = new float[3* nvert];

  if (!vertex || !access || !normal || !triangle || !properties) {
    delete [] vertex;
    delete [] access;
    delete [] normal;
    delete [] triangle;
    delete [] properties;
    printf("graspplugin) Failed vertex/access/normal/triangle allocations.\n");
    return MOLFILE_ERROR;
  }

  // ungrasp says:
  //    if (filetype.eq.1) then integer*2
  //    if (filetype.eq.2) then integer*4

  // And read them in.  Who needs error checking?
  TRASH; 
  fread(vertex, 3 * sizeof(float), nvert, infile); 
  TRASH;
  TRASH; 
  fread(access, 3 * sizeof(float), nvert, infile); 
  TRASH;
  TRASH; 
  fread(normal, 3 * sizeof(float), nvert, infile); 
  TRASH;
 
  if (is_little_endian()) {
    swap4_aligned(vertex, 3*nvert);
    swap4_aligned(access, 3*nvert);
    swap4_aligned(normal, 3*nvert);
  }

  if (gfiletype == 2) {
    TRASH; 
    fread(triangle, 3 * sizeof(int), ntriangles, infile); 
    TRASH;
    TRASH;
    fread(properties, 3 * sizeof(float), nvert, infile);
    if (is_little_endian()) {
      swap4_aligned(triangle, 3*ntriangles);
      swap4_aligned(properties, 3*nvert);
    }
  } else {
#if 1
    int i;
    short *striangle = new short[3 * ntriangles];
    if (!striangle) {
      delete [] vertex;
      delete [] access;
      delete [] normal;
      delete [] triangle;
      delete [] properties;
      printf("graspplugin) Failed short triangle allocation.\n");
      return MOLFILE_ERROR;
    }

    TRASH;
    fread(striangle, sizeof(short), 3 * ntriangles, infile);
    TRASH;
    TRASH;
    fread(properties, sizeof(float), 3 * nvert, infile);
    
    if (is_little_endian()) {
    swap2_aligned(striangle, 3 * ntriangles);
    swap4_aligned(properties, 3*nvert);}
    
    for (i=0; i<3*ntriangles; i++) {
      triangle[i] = striangle[i];
    }
    delete [] striangle;  
    
#else
    // do it the slow way (converting from short to int)
    int i;
    short tmp[3];
    TRASH;
    for (i=0; i<ntriangles; i++) {
      fread(tmp, sizeof(short), 3, infile);
      if (is_little_endian()) swap2_aligned(tmp, 3);
      triangle[3*i+0] = tmp[0];
      triangle[3*i+1] = tmp[1];
      triangle[3*i+2] = tmp[2];
    }
    TRASH;
    TRASH;
    fread(properties, sizeof(float), 3 * nvert, infile);
      if (is_little_endian())
      swap4_aligned(properties, 3*nvert);

#endif
  }   
  
  /// Gets properties:  potentials, curvature, distances, gproperty, g2property or vertexcolor 
  Get_Property_Values(&datax, properties, colores, nvert);
  
  // And draw things
  grasp->graphics = new molfile_graphics_t[3*ntriangles];
  int vert1, vert2, vert3;

  for (int tri_count = 0; tri_count < ntriangles; tri_count++) {
    vert1 = triangle[3*tri_count+0] - 1;  // from 1-based FORTRAN
    vert2 = triangle[3*tri_count+1] - 1;  // to 0-based C++
    vert3 = triangle[3*tri_count+2] - 1;

    if (vert1 <      0 || vert2 <      0 || vert3 <      0 ||
        vert1 >= nvert || vert2 >= nvert || vert3 >= nvert) {
      printf("graspplugin) Error, out-of-range vertex index, aborting.\n"); 
      delete [] vertex;
      delete [] access;
      delete [] normal;
      delete [] triangle;
      delete [] properties;
      return MOLFILE_ERROR;
    }

    grasp->graphics[2*tri_count  ].type = MOLFILE_TRINORM;
    grasp->graphics[2*tri_count+1].type = MOLFILE_NORMS;
    grasp->graphics[2*tri_count+2].type = MOLFILE_COLOR;
    
    float *tridata =  grasp->graphics[2*tri_count  ].data;
    float *normdata = grasp->graphics[2*tri_count+1].data;
    float *colordata = grasp->graphics[2*tri_count+2].data;
        
    memcpy(tridata  , vertex+3*vert1, 3*sizeof(float));
    memcpy(tridata+3, vertex+3*vert2, 3*sizeof(float));
    memcpy(tridata+6, vertex+3*vert3, 3*sizeof(float));
    
    memcpy(normdata  , normal+3*vert1, 3*sizeof(float));
    memcpy(normdata+3, normal+3*vert2, 3*sizeof(float));
    memcpy(normdata+6, normal+3*vert3, 3*sizeof(float));
    
    memcpy(colordata  , properties+3*vert1, 3*sizeof(float));
    memcpy(colordata+3, properties+3*vert2, 3*sizeof(float));
    memcpy(colordata+6, properties+3*vert3, 3*sizeof(float));
  } 

  *nelem = 2*ntriangles;
  *data = grasp->graphics;

  delete [] triangle;
  delete [] normal;
  delete [] access;
  delete [] vertex;
  delete [] properties;

  return MOLFILE_SUCCESS;
}


static void close_file_read(void *v) {
  grasp_t *grasp = (grasp_t *)v;
  fclose(grasp->fd);
  delete [] grasp->graphics;
  delete grasp;
}


/*
 * Initialization stuff here
 */
static molfile_plugin_t plugin;

VMDPLUGIN_EXTERN int VMDPLUGIN_init(void) {
  memset(&plugin, 0, sizeof(molfile_plugin_t));
  plugin.abiversion = vmdplugin_ABIVERSION;
  plugin.type = MOLFILE_PLUGIN_TYPE;
  plugin.name = "grasp";
  plugin.prettyname = "GRASP";
  plugin.author = "Justin Gullingsrud, John Stone";
  plugin.majorv = 0;
  plugin.minorv = 7;
  plugin.is_reentrant = VMDPLUGIN_THREADSAFE;
  plugin.filename_extension = "srf,SRF,grasp";
  plugin.open_file_read = open_file_read;
  plugin.read_rawgraphics = read_rawgraphics;
  plugin.close_file_read = close_file_read;
  return VMDPLUGIN_SUCCESS;
}

VMDPLUGIN_EXTERN int VMDPLUGIN_register(void *v, vmdplugin_register_cb cb) {
  (*cb)(v, (vmdplugin_t *)&plugin);
  return VMDPLUGIN_SUCCESS;
}

VMDPLUGIN_EXTERN int VMDPLUGIN_fini(void) { return VMDPLUGIN_SUCCESS; }



