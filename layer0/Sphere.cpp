

/* 
A* -------------------------------------------------------------------
B* This file contains source code for the PyMOL computer program
C* Copyright (c) Schrodinger, LLC. 
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

#include"os_predef.h"
#include"os_std.h"
#include"os_gl.h"

#include"Base.h"
#include"Sphere.h"
#include"Vector.h"
#include"Err.h"

#include"MemoryDebug.h"

#define FAST_SPHERE_INIT

#ifndef FAST_SPHERE_INIT
/* Twelve vertices of icosahedron on unit sphere */
#define tau 0.8506508084F       /* t=(1+sqrt(5))/2, tau=t/sqrt(1+t^2)  */
#define one 0.5257311121F       /* one=1/sqrt(1+t^2) , unit sphere     */

static const float start_points[13][3] = {
  {tau, one, 0},
  {-tau, one, 0},
  {-tau, -one, 0},
  {tau, -one, 0},
  {one, 0, tau},
  {one, 0, -tau},
  {-one, 0, -tau},
  {-one, 0, tau},
  {0, tau, one},
  {0, -tau, one},
  {0, -tau, -one},
  {0, tau, -one}
};

static const int icosahedron[21][3] = {
  {4, 8, 7},
  {4, 7, 9},
  {5, 6, 11},
  {5, 10, 6},
  {0, 4, 3},
  {0, 3, 5},
  {2, 7, 1},
  {2, 1, 6},
  {8, 0, 11},
  {8, 11, 1},
  {9, 10, 3},
  {9, 2, 10},
  {8, 4, 0},
  {11, 0, 5},
  {4, 9, 3},
  {5, 3, 10},
  {7, 8, 1},
  {6, 1, 11},
  {7, 2, 9},
  {6, 10, 2}
};
#endif

static const int mesh[30][2] = {
  {0, 3},
  {0, 4},
  {0, 5},
  {0, 8},
  {0, 11},
  {1, 2},
  {1, 6},
  {1, 7},
  {1, 8},
  {1, 11},
  {2, 6},
  {2, 7},
  {2, 9},
  {2, 10},
  {3, 4},
  {3, 5},
  {3, 9},
  {3, 10},
  {4, 7},
  {4, 8},
  {4, 9},
  {5, 6},
  {5, 10},
  {5, 11},
  {6, 10},
  {6, 11},
  {7, 8},
  {7, 9},
  {8, 11},
  {9, 10}
};

#ifdef FAST_SPHERE_INIT
#include"SphereData.h"

#else

static SphereRec *MakeDotSphere(PyMOLGlobals * G, int level);

#endif

#ifndef FAST_SPHERE_INIT
static void SphereDumpAll(CSphere *I)
{
  FILE *f;
  int i, dot_total, a, c, strip_total, seq_total, tri_total;
  SphereRec *sp;
  f = fopen("SphereData.h", "w");

  fprintf(f, "static int Sphere_NSpheres = %d;\n", NUMBER_OF_SPHERE_LEVELS);

  fprintf(f, "static int Sphere_NTri[%d] = {\n", NUMBER_OF_SPHERE_LEVELS);
  for (i=0; i < NUMBER_OF_SPHERE_LEVELS; i++){
    fprintf(f, " %d, ", I->Sphere[i]->NTri);
  }  
  fprintf(f, "\n};\n");

  fprintf(f, "static int Sphere_NStrip[%d] = {\n", NUMBER_OF_SPHERE_LEVELS);
  for (i=0; i < NUMBER_OF_SPHERE_LEVELS; i++){
    fprintf(f, " %d, ", I->Sphere[i]->NStrip);
  }  
  fprintf(f, "\n};\n");

  fprintf(f, "static int Sphere_NVertTot[%d] = {\n", NUMBER_OF_SPHERE_LEVELS);
  for (i=0; i < NUMBER_OF_SPHERE_LEVELS; i++){
    fprintf(f, " %d, ", I->Sphere[i]->NVertTot);
  }  
  fprintf(f, "\n};\n");

  fprintf(f, "static int Sphere_nDot[%d] = {\n", NUMBER_OF_SPHERE_LEVELS);
  for (i=0; i < NUMBER_OF_SPHERE_LEVELS; i++){
    fprintf(f, " %d, ", I->Sphere[i]->nDot);
  }  
  fprintf(f, "\n};\n");

  fprintf(f, "static int Sphere_dot_start[%d] = {\n", NUMBER_OF_SPHERE_LEVELS);
  dot_total = 0;
  for (i=0; i < NUMBER_OF_SPHERE_LEVELS; i++){
    fprintf(f, " %d, ", dot_total);
    dot_total += I->Sphere[i]->nDot;
  }  
  fprintf(f, "\n};\n");

  fprintf(f, "static float Sphere_dot[][3] = {\n");
  for (i=0; i < NUMBER_OF_SPHERE_LEVELS; i++){
    sp = I->Sphere[i];
    fprintf(f, "/* dots for Sphere #%d */\n", i);
    for(a = 0; a < sp->nDot; a++) {
      fprintf(f, "{ %15.12fF, %15.12fF, %15.12fF },\n",
	      sp->dot[a][0], sp->dot[a][1], sp->dot[a][2]);
    }
  }
  fprintf(f, "};\n");

  fprintf(f, "static float Sphere_area[] = {\n");
  for (i=0; i < NUMBER_OF_SPHERE_LEVELS; i++){
    sp = I->Sphere[i];
    fprintf(f, "/* areas for Sphere #%d */\n", i);
    c = 0;
    for(a = 0; a < sp->nDot; a++) {
      fprintf(f, "%15.12fF,", sp->area[a]);
    c = (c + 1) % 4;
      if (!c)
	fprintf(f, "\n");	
    }
    if (c)
      fprintf(f, "\n");
  }
  fprintf(f, "};\n");


  fprintf(f, "static int Sphere_StripLen_start[%d] = {\n", NUMBER_OF_SPHERE_LEVELS);
  strip_total = 0;
  for (i=0; i < NUMBER_OF_SPHERE_LEVELS; i++){
    fprintf(f, " %d, ", strip_total);
    strip_total += I->Sphere[i]->NStrip;
  }  
  fprintf(f, "\n};\n");

  fprintf(f, "static int Sphere_StripLen[] = {\n");
  for (i=0; i < NUMBER_OF_SPHERE_LEVELS; i++){
    sp = I->Sphere[i];
    fprintf(f, "/* StripLen for Sphere #%d */\n", i);
    c = 0;
    for(a = 0; a < sp->NStrip; a++) {
      fprintf(f, "%6d,", sp->StripLen[a]);
      c = (c + 1) % 10;
      if(!c)
	fprintf(f, "\n");
    }
    if (c)
      fprintf(f, "\n");
  }
  fprintf(f, "};\n");

  fprintf(f, "static int Sphere_Sequence_start[%d] = {\n", NUMBER_OF_SPHERE_LEVELS);
  seq_total = 0;
  for (i=0; i < NUMBER_OF_SPHERE_LEVELS; i++){
    fprintf(f, " %d, ", seq_total);
    seq_total += I->Sphere[i]->NVertTot;
  }  
  fprintf(f, "\n};\n");

  fprintf(f, "static int Sphere_Sequence[] = {\n");
  for (i=0; i < NUMBER_OF_SPHERE_LEVELS; i++){
    sp = I->Sphere[i];
    fprintf(f, "/* Sequence for Sphere #%d */\n", i);
    c = 0;
    for(a = 0; a < sp->NVertTot; a++) {
      fprintf(f, "%6d,", sp->Sequence[a]);
      c = (c + 1) % 10;
      if(!c)
	fprintf(f, "\n");
    }
    if (c)
      fprintf(f, "\n");
  }
  fprintf(f, "};\n");

  fprintf(f, "static int Sphere_Tri_start[%d] = {\n", NUMBER_OF_SPHERE_LEVELS);
  tri_total = 0;
  for (i=0; i < NUMBER_OF_SPHERE_LEVELS; i++){
    fprintf(f, " %d, ", tri_total);
    tri_total += 3 * I->Sphere[i]->NTri;
  }  
  fprintf(f, "\n};\n");

  fprintf(f, "static int Sphere_Tri[] = {\n");
  for (i=0; i < NUMBER_OF_SPHERE_LEVELS; i++){
    sp = I->Sphere[i];
    fprintf(f, "/* Tri for Sphere #%d */\n", i);
    c = 0;
    for(a = 0; a < 3* sp->NTri; a++) {
      fprintf(f, "%6d,", sp->Tri[a]);
      c = (c + 1) % 10;
      if(!c)
	fprintf(f, "\n");
    }
    if (c)
      fprintf(f, "\n");
  }
  fprintf(f, "};\n");

  fclose(f);
}
#endif

void SphereInit(PyMOLGlobals * G)
{
  CSphere *I = (G->Sphere = pymol::calloc<CSphere>(1));

#ifdef FAST_SPHERE_INIT
  I->Array = pymol::malloc<SphereRec>(Sphere_NSpheres);

  {
    int i;
    for (i=0; i<Sphere_NSpheres; i++){
      I->Array[i].area = &Sphere_area[Sphere_dot_start[i]];
      I->Array[i].dot = &Sphere_dot[Sphere_dot_start[i]];
      I->Array[i].StripLen = &Sphere_StripLen[Sphere_StripLen_start[i]];
      I->Array[i].Sequence = &Sphere_Sequence[Sphere_Sequence_start[i]];
      I->Array[i].NStrip = Sphere_NStrip[i];
      I->Array[i].NVertTot = Sphere_NVertTot[i];
      I->Array[i].nDot = Sphere_nDot[i];
      I->Array[i].Tri = &Sphere_Tri[Sphere_Tri_start[i]];
      I->Array[i].NTri = Sphere_NTri[i];
      
      if (i){
	I->Array[i].Mesh = NULL;
	I->Array[i].NMesh = 0;
      } else {
	I->Array[i].Mesh = (int *) (void *) mesh;
	I->Array[i].NMesh = 30;
      }
      I->Sphere[i] = &I->Array[i];
    }
  }
#else
  {
    int i;
    for (i=0; i<NUMBER_OF_SPHERE_LEVELS; i++){
      I->Sphere[i] = MakeDotSphere(G, i);
    }
    SphereDumpAll(I);
  }
#endif

}

#ifndef FAST_SPHERE_INIT
static void SpherePurge(SphereRec * I)
{
  /* NOTE: S->Mesh is not currently a pointer */
  mfree(I->dot);
  mfree(I->area);
  mfree(I->StripLen);
  mfree(I->Sequence);
  mfree(I->Tri);
  FreeP(I);
}
#endif

void SphereFree(PyMOLGlobals * G)
{
  CSphere *I = G->Sphere;

#ifndef FAST_SPHERE_INIT
  SpherePurge(I->Sphere[0]);
  SpherePurge(I->Sphere[1]);
  SpherePurge(I->Sphere[2]);
  SpherePurge(I->Sphere[3]);
  SpherePurge(I->Sphere[4]);
#else
  FreeP(I->Array);
#endif
  FreeP(I);
}


/* private stuff */

#ifndef FAST_SPHERE_INIT

// MAXDOT : 12, 42, 162, 642, 2562 ... :: 12 + (30 + 120 + 480 + 1920 + ... ) :: 12 + ( 30 + (30*4) + (30*4*4) + (30*4*4*4) + ... )
// MAXTRI : 80, 320, 1280, 5120, ... :: 20*1 + 20*4 + 20*4*4 + 20*4*4*4 + 20*4*4*4*4 :: 

//For NUMBER_OF_SPHERE_LEVELS=6
//#define MAXDOT 12900   // 12800
//#define MAXTRI 20500   // 20480

//For NUMBER_OF_SPHERE_LEVELS=5
#define MAXDOT 2600      // 2562
#define MAXTRI 5200      // 5120

typedef int EdgeCol[MAXDOT];    /* should move these into dynamic storage to save 3MB  mem */
typedef EdgeCol EdgeArray[MAXDOT];

typedef int Triangle[3];

typedef struct {

  float *Dot;
  EdgeArray *EdgeRef;
  Triangle *Tri;
  int NDot, NTri;

} SphereBuilderRec;

static void MakeVertex(SphereBuilderRec * S, int d1, int d2)
{
  if((*S->EdgeRef)[d1][d2] < 0) {
    average3f(S->Dot + (3 * d1), S->Dot + (3 * d2), S->Dot + (3 * S->NDot));
    (*S->EdgeRef)[d1][d2] = S->NDot;
    (*S->EdgeRef)[d2][d1] = S->NDot;
    normalize3f(S->Dot + (3 * S->NDot));
    S->NDot++;
  }
}

static float SphericalAngle(SphereBuilderRec * S, int d0, int d1, int d2)
{
  Vector3f v1, v2, s1, s2;

  /* map vector onto surface of sphere and measure angle */
  subtract3f(S->Dot + (3 * d1), S->Dot + (3 * d0), v1);
  subtract3f(S->Dot + (3 * d2), S->Dot + (3 * d0), v2);

  remove_component3f(v1, S->Dot + (3 * d0), s1);
  remove_component3f(v2, S->Dot + (3 * d0), s2);
  return (get_angle3f(s1, s2));

}

static SphereRec *MakeDotSphere(PyMOLGlobals * G, int level)
{
  SphereRec *result;
  int *TriFlag;
  int a, b, c, h, k, l, curTri, n, it;
  float area, sumArea = 0.0;
  int nStrip, *q, *s;
  int nVertTot;
  int flag;
  float vt1[3], vt2[3], vt[3];
  SphereBuilderRec SBuild, *S;
  S = &SBuild;

  S->Dot = pymol::malloc<float>(3 * MAXDOT);
  ErrChkPtr(G, S->Dot);
  S->EdgeRef = pymol::malloc<EdgeArray>(1);
  ErrChkPtr(G, S->EdgeRef);
  S->Tri = pymol::malloc<Triangle>(MAXTRI);
  ErrChkPtr(G, S->Tri);
  TriFlag = pymol::malloc<int>(MAXTRI);
  ErrChkPtr(G, TriFlag);

  S->NDot = 12;
  for(a = 0; a < S->NDot; a++) {
    for(c = 0; c < 3; c++)
      S->Dot[3 * a + c] = start_points[a][c];
    normalize3f(S->Dot + (3 * a));
  }

  S->NTri = 20;
  for(a = 0; a < S->NTri; a++)
    for(c = 0; c < 3; c++)
      S->Tri[a][c] = icosahedron[a][c];

  for(a = 0; a < MAXDOT; a++)
    for(b = 0; b < MAXDOT; b++)
      (*S->EdgeRef)[a][b] = -1;

  if(level > (NUMBER_OF_SPHERE_LEVELS-1))
    level = (NUMBER_OF_SPHERE_LEVELS-1);

  for(c = 0; c < level; c++) {
    /* create new vertices */
    for(a = 0; a < S->NTri; a++) {
      MakeVertex(S, S->Tri[a][0], S->Tri[a][1]);
      MakeVertex(S, S->Tri[a][1], S->Tri[a][2]);
      MakeVertex(S, S->Tri[a][0], S->Tri[a][2]);
    }
    /* create new triangles */
    curTri = S->NTri;
    for(a = 0; a < curTri; a++) {
      h = S->Tri[a][0];
      k = S->Tri[a][1];
      l = S->Tri[a][2];

      S->Tri[a][0] = h;
      S->Tri[a][1] = (*S->EdgeRef)[h][k];
      S->Tri[a][2] = (*S->EdgeRef)[h][l];

      S->Tri[S->NTri][0] = k;
      S->Tri[S->NTri][1] = (*S->EdgeRef)[k][h];
      S->Tri[S->NTri][2] = (*S->EdgeRef)[k][l];
      S->NTri++;

      S->Tri[S->NTri][0] = l;
      S->Tri[S->NTri][1] = (*S->EdgeRef)[l][h];
      S->Tri[S->NTri][2] = (*S->EdgeRef)[l][k];
      S->NTri++;

      S->Tri[S->NTri][0] = (*S->EdgeRef)[h][k];
      S->Tri[S->NTri][1] = (*S->EdgeRef)[k][l];
      S->Tri[S->NTri][2] = (*S->EdgeRef)[l][h];
      S->NTri++;
    }
    //    printf( "MakeDotSphere: Level: %i  S->NTri: %i\n",c, S->NTri); 
  }
  //  printf(" MakeDotSphere: NDot %i S->NTri %i\n",S->NDot,S->NTri);
  result = pymol::malloc<SphereRec>(1);
  ErrChkPtr(G, result);
  result->dot = pymol::malloc<Vector3f>(S->NDot);
  ErrChkPtr(G, result->dot);
  result->area = pymol::malloc<float>(S->NDot);
  ErrChkPtr(G, result->area);
  result->StripLen = pymol::malloc<int>(S->NTri * 3);
  ErrChkPtr(G, result->StripLen);
  result->Sequence = pymol::malloc<int>(S->NTri * 3);
  ErrChkPtr(G, result->Sequence);

  for(a = 0; a < S->NDot; a++) {
    for(c = 0; c < 3; c++)
      result->dot[a][c] = *(S->Dot + (3 * a + c));
    result->area[a] = 0.0;
  }

  /* fix normals so that v1-v0 x v2-v0 is the correct normal */

  for(a = 0; a < S->NTri; a++) {
    subtract3f(result->dot[S->Tri[a][1]], result->dot[S->Tri[a][0]], vt1);
    subtract3f(result->dot[S->Tri[a][2]], result->dot[S->Tri[a][0]], vt2);
    cross_product3f(vt1, vt2, vt);
    if(dot_product3f(vt, result->dot[S->Tri[a][0]]) < 0.0) {    /* if wrong, then interchange */
      it = S->Tri[a][2];
      S->Tri[a][2] = S->Tri[a][1];
      S->Tri[a][1] = it;
    }
  }

  for(a = 0; a < S->NTri; a++) {
    area = (float) (SphericalAngle(S, S->Tri[a][0], S->Tri[a][1], S->Tri[a][2]) +
                    SphericalAngle(S, S->Tri[a][1], S->Tri[a][0], S->Tri[a][2]) +
                    SphericalAngle(S, S->Tri[a][2], S->Tri[a][0], S->Tri[a][1]) - cPI);
    /* multiply by r^2 to get area */
    sumArea += area;
    area /= 3.0;
    result->area[S->Tri[a][0]] += area;
    result->area[S->Tri[a][1]] += area;
    result->area[S->Tri[a][2]] += area;
  }

  if(fabs(sumArea - (4 * cPI)) > 0.001) {
    printf(" MakeDotSphere: sumArea: %8.6f which is %8.6f Pi\n", sumArea, sumArea / cPI);
    ErrFatal(G, "MakeDotSphere", "Area of sphere does not sum to 4*pi!\n");
  }

  for(a = 0; a < S->NTri; a++)
    TriFlag[a] = false;

  nStrip = 0;
  nVertTot = 0;
  s = result->StripLen;
  q = result->Sequence;

  /* tesselate the sphere in a semi-efficient fashion...this could definitely be improved */

  flag = true;
  while(flag) {
    flag = false;
    a = 0;
    while(a < S->NTri) {
      if(!TriFlag[a]) {
        flag = true;

        TriFlag[a] = true;
        *(q++) = S->Tri[a][0];
        *(q++) = S->Tri[a][1];
        *(q++) = S->Tri[a][2];
        n = 3;

        b = 0;
        while(b < S->NTri) {
          if(!TriFlag[b]) {
            if(((S->Tri[b][0] == q[-2]) && (S->Tri[b][1] == q[-1])) ||
               ((S->Tri[b][0] == q[-1]) && (S->Tri[b][1] == q[-2]))) {
              *(q++) = S->Tri[b][2];
              TriFlag[b] = true;
              b = 0;
              n++;
            } else if(((S->Tri[b][0] == q[-2]) && (S->Tri[b][2] == q[-1])) ||
                      ((S->Tri[b][0] == q[-1]) && (S->Tri[b][2] == q[-2]))) {
              *(q++) = S->Tri[b][1];
              TriFlag[b] = true;
              b = 0;
              n++;
            } else if(((S->Tri[b][2] == q[-2]) && (S->Tri[b][1] == q[-1])) ||
                      ((S->Tri[b][2] == q[-1]) && (S->Tri[b][1] == q[-2]))) {
              *(q++) = S->Tri[b][0];
              TriFlag[b] = true;
              b = 0;
              n++;
            }
          }
          b++;
        }
        if(n == 3) {
          q[-3] = S->Tri[a][1];
          q[-2] = S->Tri[a][2];
          q[-1] = S->Tri[a][0];

          b = 0;
          while(b < S->NTri) {
            if(!TriFlag[b]) {
              if(((S->Tri[b][0] == q[-2]) && (S->Tri[b][1] == q[-1])) ||
                 ((S->Tri[b][0] == q[-1]) && (S->Tri[b][1] == q[-2]))) {
                *(q++) = S->Tri[b][2];
                TriFlag[b] = true;
                b = 0;
                n++;
              } else if(((S->Tri[b][0] == q[-2]) && (S->Tri[b][2] == q[-1])) ||
                        ((S->Tri[b][0] == q[-1]) && (S->Tri[b][2] == q[-2]))) {
                *(q++) = S->Tri[b][1];
                TriFlag[b] = true;
                b = 0;
                n++;
              } else if(((S->Tri[b][2] == q[-2]) && (S->Tri[b][1] == q[-1])) ||
                        ((S->Tri[b][2] == q[-1]) && (S->Tri[b][1] == q[-2]))) {
                *(q++) = S->Tri[b][0];
                TriFlag[b] = true;
                b = 0;
                n++;
              }
            }
            b++;
          }
        }
        if(n == 3) {
          q[-3] = S->Tri[a][2];
          q[-2] = S->Tri[a][0];
          q[-1] = S->Tri[a][1];
          b = 0;
          while(b < S->NTri) {
            if(!TriFlag[b]) {
              if(((S->Tri[b][0] == q[-2]) && (S->Tri[b][1] == q[-1])) ||
                 ((S->Tri[b][0] == q[-1]) && (S->Tri[b][1] == q[-2]))) {
                *(q++) = S->Tri[b][2];
                TriFlag[b] = true;
                b = 0;
                n++;
              } else if(((S->Tri[b][0] == q[-2]) && (S->Tri[b][2] == q[-1])) ||
                        ((S->Tri[b][0] == q[-1]) && (S->Tri[b][2] == q[-2]))) {
                *(q++) = S->Tri[b][1];
                TriFlag[b] = true;
                b = 0;
                n++;
              } else if(((S->Tri[b][2] == q[-2]) && (S->Tri[b][1] == q[-1])) ||
                        ((S->Tri[b][2] == q[-1]) && (S->Tri[b][1] == q[-2]))) {
                *(q++) = S->Tri[b][0];
                TriFlag[b] = true;
                b = 0;
                n++;
              }
            }
            b++;
          }
        }
        *(s++) = n;
        nVertTot += n;
        nStrip++;
      }
      a++;
    }
  }
  mfree(S->Dot);
  mfree(S->EdgeRef);
  mfree(TriFlag);
  result->Tri = (int *) S->Tri;
  result->Tri = pymol::realloc(result->Tri, S->NTri * 3);
  result->NTri = S->NTri;
  result->StripLen = pymol::realloc(result->StripLen, nStrip);
  result->Sequence = pymol::realloc(result->Sequence, nVertTot);
  result->dot = pymol::realloc(result->dot, S->NDot);
  result->area = pymol::realloc(result->area, S->NDot);
  result->nDot = S->NDot;
  result->NStrip = nStrip;
  result->NVertTot = nVertTot;
  result->Mesh = NULL;
  result->NMesh = 0;
  if(!level) {                  /* provide mesh for S->Sphere[0] only...rest, to do. */
    result->Mesh = (int *) mesh;
    result->NMesh = 30;
  }
  /*
     q=result->Sequence;
     for(a=0;a<result->NStrip;a++)
     {
     printf("%d:",result->StripLen[a]); 
     for(b=0;b<result->StripLen[a];b++)
     {
     printf("%d ",*(q++));
     }
     printf("\n");
     }
   */

  return (result);
}

#endif

void SphereRender(PyMOLGlobals * G, int level, const float *centroid, const float *color, float alpha, float radius){
#ifndef PURE_OPENGL_ES_2
  SphereRec *sp = G->Sphere->Sphere[level];
  int a, cc;
  int *q = sp->Sequence;
  float pt[3];
  if (color)
    glColor4f(color[0], color[1], color[2], alpha);
  for(a = 0; a < sp->NStrip; a++) {
    glBegin(GL_TRIANGLE_STRIP);
    cc = sp->StripLen[a];
    while(cc--) {
      glNormal3fv(sp->dot[*q]);
      mult3f(sp->dot[*q], radius, pt);
      add3f(centroid, pt, pt);
      glVertex3fv(pt);
      q++;
    }
    glEnd();
  }
#endif
}

SphereRec* GetSpheroidSphereRec(PyMOLGlobals* G)
{
  return G->Sphere->Sphere[2];
}
