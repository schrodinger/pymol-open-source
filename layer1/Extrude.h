
/* 
A* -------------------------------------------------------------------
B* This file contains source code for the PyMOL computer program
C* copyright 1998-2000 by Warren Lyford Delano of DeLano Scientific. 
D* -------------------------------------------------------------------
E* It is unlawful to modify or remove this copyright notice.
F* -------------------------------------------------------------------
G* Please see the accompanying LICENSE file for further information. 
H* -------------------------------------------------------------------
I* Additional authors of this source file include:
-* 
-* 
-*
Z* -------------------------------------------------------------------
*/
#ifndef _H_Extrude
#define _H_Extrude

#include"Ray.h"
#include"CGO.h"
#include"ObjectMolecule.h"
#include"CoordSet.h"

typedef struct {
  PyMOLGlobals *G;
  int N;                        /* number of points in the extrusion segment */

  float *p;                     /* points */
  float *n;                     /* normals (3x3f) at each point */
  float *c;                     /* colors */
  int *i;                       /* atom indices */

  float r;
  float *sf;                    /* scale factors for variable-width extrusions (single point) */

  float *sv, *tv;
  float *sn, *tn;
  int Ns;

} CExtrude;

CExtrude *ExtrudeNew(PyMOLGlobals * G);

CExtrude *ExtrudeCopyPointsNormalsColors(CExtrude * orig);

int ExtrudeAllocPointsNormalsColors(CExtrude * I, int n);
void ExtrudeTruncate(CExtrude * I, int n);

/* SAUSAGE: look-up table -- nearest neighbors to extrusion points */
void ExtrudeMakePuttyLUT(CExtrude * I, ObjectMolecule * Sauce_obj);

void ExtrudeFree(CExtrude * I);

int ExtrudeCircle(CExtrude * I, int n, float size);
int ExtrudeRectangle(CExtrude * I, float width, float length, int mode);
int ExtrudeOval(CExtrude * I, int n, float width, float length);

int ExtrudeComputePuttyScaleFactors(CExtrude * I, ObjectMolecule * obj,
				    int transform,
				    float mean, float stdev, float min, float max,
				    float power, float range,
				    float min_scale, float max_scale, int smooth_window);

void ExtrudeBuildNormals1f(CExtrude * I);
void ExtrudeBuildNormals2f(CExtrude * I);
int ExtrudeComputeTangents(CExtrude * I);
int ExtrudeCylindersToCGO(CExtrude * I, CGO *cgo, float tube_radius, short is_picking);
int ExtrudeCGOSurfaceTube(CExtrude * I, CGO * cgo, int cap, float *color_override, short use_spheres);
void ExtrudeCGOSurfaceTubeToCylinders(CExtrude * I, CGO * cgo, int cap, float *color_override);
int ExtrudeCGOSurfaceVariableTube(CExtrude * I, CGO * cgo, int cap);

int ExtrudeCGOSurfacePolygon(CExtrude * I, CGO * cgo, int cap, float *color_override);
int ExtrudeCGOSurfacePolygonTaper(CExtrude * I, CGO * cgo,
                                   int sampling, float *color_override);
int ExtrudeCGOSurfaceStrand(CExtrude * I, CGO * cgo, int sampling,
                             float *color_override);
void ExtrudeCGOTraceFrame(CExtrude * I, CGO * cgo);
void ExtrudeCGOTrace(CExtrude * I, CGO * cgo);
void ExtrudeCGOTraceAxes(CExtrude * I, CGO * cgo);
int ExtrudeDumbbell1(CExtrude * I, float width, float length, int mode);
int ExtrudeDumbbell2(CExtrude * I, int n, int sign, float length, float size);
void ExtrudeDumbbellEdge(CExtrude * I, int samp, int sign, float length);

#endif
