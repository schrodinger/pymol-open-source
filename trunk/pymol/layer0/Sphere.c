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

#include"os_predef.h"
#include"os_std.h"

#include"Base.h"
#include"Sphere.h"
#include"Vector.h"
#include"Err.h"

#include"MemoryDebug.h"

/* Twelve vertices of icosahedron on unit sphere */
#define tau 0.8506508084F      /* t=(1+sqrt(5))/2, tau=t/sqrt(1+t^2)  */
#define one 0.5257311121F      /* one=1/sqrt(1+t^2) , unit sphere     */

static const float start_points[13][3] = {
{  tau,  one,    0 },
{ -tau,  one,    0 },
{ -tau, -one,    0 },
{  tau, -one,    0 },
{  one,   0 ,  tau },
{  one,   0 , -tau },
{ -one,   0 , -tau },
{ -one,   0 ,  tau },
{   0 ,  tau,  one },
{   0 , -tau,  one },
{   0 , -tau, -one },
{   0 ,  tau, -one }};

static const int icosahedron[21][3] = {
    { 4, 8, 7 },
    { 4, 7, 9 },
    { 5, 6, 11  },
    { 5, 10, 6 },
    { 0, 4, 3 },
    { 0, 3, 5 },
    { 2, 7, 1 },
    { 2, 1, 6 },
    { 8, 0, 11 },
    { 8, 11, 1 },
    { 9, 10, 3 },
    { 9, 2, 10 },
    { 8, 4, 0 },
    { 11, 0, 5 },
    { 4, 9, 3 },
    { 5, 3, 10 },
    { 7, 8, 1 },
    { 6, 1, 11 },
    { 7, 2, 9 },
    { 6, 10, 2 }
};

static const int mesh[30][2] = {
{ 0 , 3 },
{ 0 , 4 },
{ 0 , 5 },
{ 0 , 8 },
{ 0 , 11 },
{ 1 , 2 },
{ 1 , 6 },
{ 1 , 7 },
{ 1 , 8 },
{ 1 , 11 },
{ 2 , 6 },
{ 2 , 7 },
{ 2 , 9 },
{ 2 , 10 },
{ 3 , 4 },
{ 3 , 5 },
{ 3 , 9 },
{ 3 , 10 },
{ 4 , 7 },
{ 4 , 8 },
{ 4 , 9 },
{ 5 , 6 },
{ 5 , 10 },
{ 5 , 11 },
{ 6 , 10 },
{ 6 , 11 },
{ 7 , 8 },
{ 7 , 9 },
{ 8 , 11 },
{ 9 , 10 }                                                                        
};


#define FAST_SPHERE_INIT


#ifdef FAST_SPHERE_INIT
#include"SphereData.h"

#else

static SphereRec *MakeDotSphere(PyMOLGlobals *G,int level);

#endif

#if 0
static void SphereDump(FILE *f,char *prefix,SphereRec *sp)
{
  int a;
  int c;
  fprintf(f,"static int %s_NTri = %d;\n",prefix,sp->NTri);
  fprintf(f,"static int %s_NStrip = %d;\n",prefix,sp->NStrip);
  fprintf(f,"static int %s_NVertTot = %d;\n",prefix,sp->NVertTot);
  fprintf(f,"static int %s_nDot = %d;\n",prefix,sp->nDot);
  fprintf(f,"static float %s_dot[][3] = {\n",prefix);
  for(a=0;a<sp->nDot;a++) {
    fprintf(f,"{ %15.12fF, %15.12fF, %15.12fF },\n",
            sp->dot[a][0],
            sp->dot[a][1],
            sp->dot[a][2]);
  }
  fprintf(f,"};\n");
  
  fprintf(f,"static float %s_area[] = {\n",prefix);

  c = 0;
  for(a=0;a<sp->nDot;a++) {
    fprintf(f,"%15.12fF,",
            sp->area[a]);
    c = (c+1)%4;
    if(!c) fprintf(f,"\n");
  }
  fprintf(f,"};\n");

  fprintf(f,"static int %s_StripLen[] = {\n",prefix);
  c = 0;
  for(a=0;a<sp->NStrip;a++) {
    fprintf(f,"%6d,",
            sp->StripLen[a]);
    c = (c+1)%10;
    if(!c) fprintf(f,"\n");
  }
  fprintf(f,"};\n");

  fprintf(f,"static int %s_Sequence[] = {\n",prefix);
  c = 0;
  for(a=0;a<sp->NVertTot;a++) {
    fprintf(f,"%6d,",
            sp->Sequence[a]);
    c = (c+1)%10;
    if(!c) fprintf(f,"\n");
  }
  fprintf(f,"};\n");

  fprintf(f,"static int %s_Tri[] = {\n",prefix);
  c = 0;
  for(a=0;a<3*sp->NTri;a++) {
    fprintf(f,"%6d,",
            sp->Tri[a]);
    c = (c+1)%10;
    if(!c) fprintf(f,"\n");
  }
  fprintf(f,"};\n");


}

static void SphereDumpAll(void)
{
  FILE *f;
  f = fopen("SphereData.h","w");
  SphereDump(f,"Sphere0",I->Sphere[0]);
  SphereDump(f,"Sphere1",I->Sphere[1]);
  SphereDump(f,"Sphere2",I->Sphere[2]);
  SphereDump(f,"Sphere3",I->Sphere[3]);
  SphereDump(f,"Sphere4",I->Sphere[4]);
  fclose(f);
}
#endif


void SphereInit(PyMOLGlobals *G)
{
  register CSphere *I = (G->Sphere = Calloc(CSphere,1));

#ifdef FAST_SPHERE_INIT
  I->Array = Alloc(SphereRec,5);

  I->Array[0].area = Sphere0_area;
  I->Array[0].dot = Sphere0_dot;
  I->Array[0].StripLen = Sphere0_StripLen;
  I->Array[0].Sequence = Sphere0_Sequence;
  I->Array[0].NStrip = Sphere0_NStrip;
  I->Array[0].NVertTot = Sphere0_NVertTot;
  I->Array[0].nDot = Sphere0_nDot;
  I->Array[0].Tri = Sphere0_Tri;
  I->Array[0].NTri = Sphere0_NTri;
  I->Array[0].Mesh = (int*)mesh;
  I->Array[0].NMesh = 30;

  I->Array[1].area = Sphere1_area;
  I->Array[1].dot = Sphere1_dot;
  I->Array[1].StripLen = Sphere1_StripLen;
  I->Array[1].Sequence = Sphere1_Sequence;
  I->Array[1].NStrip = Sphere1_NStrip;
  I->Array[1].NVertTot = Sphere1_NVertTot;
  I->Array[1].nDot = Sphere1_nDot;
  I->Array[1].Tri = Sphere1_Tri;
  I->Array[1].NTri = Sphere1_NTri;
  I->Array[1].Mesh = NULL;
  I->Array[1].NMesh = 0;

  I->Array[2].area = Sphere2_area;
  I->Array[2].dot = Sphere2_dot;
  I->Array[2].StripLen = Sphere2_StripLen;
  I->Array[2].Sequence = Sphere2_Sequence;
  I->Array[2].NStrip = Sphere2_NStrip;
  I->Array[2].NVertTot = Sphere2_NVertTot;
  I->Array[2].nDot = Sphere2_nDot;
  I->Array[2].Tri = Sphere2_Tri;
  I->Array[2].NTri = Sphere2_NTri;
  I->Array[2].Mesh = NULL;
  I->Array[2].NMesh = 0;

  I->Array[3].area = Sphere3_area;
  I->Array[3].dot = Sphere3_dot;
  I->Array[3].StripLen = Sphere3_StripLen;
  I->Array[3].Sequence = Sphere3_Sequence;
  I->Array[3].NStrip = Sphere3_NStrip;
  I->Array[3].NVertTot = Sphere3_NVertTot;
  I->Array[3].nDot = Sphere3_nDot;
  I->Array[3].Tri = Sphere3_Tri;
  I->Array[3].NTri = Sphere3_NTri;
  I->Array[3].Mesh = NULL;
  I->Array[3].NMesh = 0;

  I->Array[4].area = Sphere4_area;
  I->Array[4].dot = Sphere4_dot;
  I->Array[4].StripLen = Sphere4_StripLen;
  I->Array[4].Sequence = Sphere4_Sequence;
  I->Array[4].NStrip = Sphere4_NStrip;
  I->Array[4].NVertTot = Sphere4_NVertTot;
  I->Array[4].nDot = Sphere4_nDot;
  I->Array[4].Tri = Sphere4_Tri;
  I->Array[4].NTri = Sphere4_NTri;
  I->Array[4].Mesh = NULL;
  I->Array[4].NMesh = 0;

  I->Sphere[0] = &I->Array[0];
  I->Sphere[1] = &I->Array[1];
  I->Sphere[2] = &I->Array[2];
  I->Sphere[3] = &I->Array[3];
  I->Sphere[4] = &I->Array[4];
#else
  I->Sphere[0] = MakeDotSphere(G,0);
  I->Sphere[1] = MakeDotSphere(G,1);
  I->Sphere[2] = MakeDotSphere(G,2);
  I->Sphere[3] = MakeDotSphere(G,3);
  I->Sphere[4] = MakeDotSphere(G,4);
  /*
  SphereDumpAll();
  */
  
#endif


}

#ifndef FAST_SPHERE_INIT
static void SpherePurge(SphereRec *I)
{
  /* NOTE: S->Mesh is not currently a pointer*/
  mfree(I->dot);
  mfree(I->area);
  mfree(I->StripLen);
  mfree(I->Sequence);
  mfree(I->Tri);
  FreeP(I);
}
#endif


void SphereFree(PyMOLGlobals *G)
{
  register CSphere *I = G->Sphere;

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

#define MAXDOT 2600
#define MAXTRI 6200

typedef int EdgeCol[MAXDOT]; /* should move these into dynamic storage to save 3MB  mem */
typedef EdgeCol EdgeArray[MAXDOT]; 

typedef int Triangle[3];

typedef struct {

  float *Dot;
  EdgeArray *EdgeRef;
  Triangle *Tri;
  int NDot,NTri;
  
} SphereBuilderRec;


static void MakeVertex(SphereBuilderRec *S,int d1,int d2)
{
  if((*S->EdgeRef)[d1][d2]<0)
	 {
		average3f(S->Dot+(3*d1),S->Dot+(3*d2),S->Dot+(3*S->NDot));
		(*S->EdgeRef)[d1][d2]=S->NDot;
		(*S->EdgeRef)[d2][d1]=S->NDot;
		normalize3f(S->Dot+(3*S->NDot));
		S->NDot++;
	 }
}

static float SphericalAngle(SphereBuilderRec *S,int d0,int d1,int d2)
{
  Vector3f v1,v2,s1,s2;

  /* map vector onto surface of sphere and measure angle */
  subtract3f(S->Dot+(3*d1),S->Dot+(3*d0),v1);
  subtract3f(S->Dot+(3*d2),S->Dot+(3*d0),v2);
  
  remove_component3f(v1,S->Dot+(3*d0),s1);
  remove_component3f(v2,S->Dot+(3*d0),s2);
  return(get_angle3f(s1,s2));
  
}

static SphereRec *MakeDotSphere(PyMOLGlobals *G,int level)
{
  SphereRec *result;
  int *TriFlag;
  int a,b,c,h,k,l,curTri,n,it;
  float area,sumArea=0.0;
  int nStrip,*q,*s;
  int nVertTot;
  int flag;
  float vt1[3],vt2[3],vt[3];
  SphereBuilderRec SBuild,*S;
  S = &SBuild;

  S->Dot=(float*)mmalloc(sizeof(float)*3*MAXDOT);
  ErrChkPtr(G,S->Dot);
  S->EdgeRef=(EdgeArray*)mmalloc(sizeof(EdgeArray));
  ErrChkPtr(G,S->EdgeRef);
  S->Tri=Alloc(Triangle,MAXTRI);
  ErrChkPtr(G,S->Tri);
  TriFlag=Alloc(int,MAXTRI);
  ErrChkPtr(G,TriFlag);
  
  S->NDot = 12;
  for(a=0;a<S->NDot;a++)
	 {
		for(c=0;c<3;c++)
		  S->Dot[3*a+c]=start_points[a][c];
		normalize3f(S->Dot+(3*a));
	 }

  S->NTri = 20;
  for(a=0;a<S->NTri;a++)
	 for(c=0;c<3;c++)
		S->Tri[a][c]=icosahedron[a][c];

  for(a=0;a<MAXDOT;a++)
	 for(b=0;b<MAXDOT;b++)
		(*S->EdgeRef)[a][b]=-1;

  if(level>4)
	 level=4;
  
  for(c=0;c<level;c++)
	 {
		/* create new vertices */
		for(a=0;a<S->NTri;a++)
		  {
			 MakeVertex(S,S->Tri[a][0],S->Tri[a][1]);
			 MakeVertex(S,S->Tri[a][1],S->Tri[a][2]);
			 MakeVertex(S,S->Tri[a][0],S->Tri[a][2]);
		  }		
		/* create new triangles */
		curTri=S->NTri;
		for(a=0;a<curTri;a++)
		  {
			 h=S->Tri[a][0];
			 k=S->Tri[a][1];
			 l=S->Tri[a][2];

			 S->Tri[a][0]=h;
			 S->Tri[a][1]=(*S->EdgeRef)[h][k];
			 S->Tri[a][2]=(*S->EdgeRef)[h][l];
			 
			 S->Tri[S->NTri][0]=k;
			 S->Tri[S->NTri][1]=(*S->EdgeRef)[k][h];
			 S->Tri[S->NTri][2]=(*S->EdgeRef)[k][l];
			 S->NTri++;
			 
			 S->Tri[S->NTri][0]=l;
			 S->Tri[S->NTri][1]=(*S->EdgeRef)[l][h];
			 S->Tri[S->NTri][2]=(*S->EdgeRef)[l][k];
			 S->NTri++;
			 
			 S->Tri[S->NTri][0]=(*S->EdgeRef)[h][k];
			 S->Tri[S->NTri][1]=(*S->EdgeRef)[k][l];
			 S->Tri[S->NTri][2]=(*S->EdgeRef)[l][h];
			 S->NTri++;		
		}
		/*		printf( "MakeDotSphere: Level: %i  S->NTri: %i\n",c, S->NTri); */
	 }
  /*  printf(" MakeDotSphere: NDot %i S->NTri %i\n",NDot,S->NTri);*/
  result= Alloc(SphereRec,1);
  ErrChkPtr(G,result);
  result->dot = Alloc(Vector3f,S->NDot);
  ErrChkPtr(G,result->dot);
  result->area = Alloc(float,S->NDot);
  ErrChkPtr(G,result->area);
  result->StripLen = Alloc(int,S->NTri*3);
  ErrChkPtr(G,result->StripLen);
  result->Sequence = Alloc(int,S->NTri*3);
  ErrChkPtr(G,result->Sequence);

  for(a=0;a<S->NDot;a++)
	 {
		for(c=0;c<3;c++)
		  result->dot[a][c]=*(S->Dot+(3*a+c));
		result->area[a] = 0.0;
	 }

  /* fix normals so that v1-v0 x v2-v0 is the correct normal */

  for(a=0;a<S->NTri;a++) {
    subtract3f(result->dot[S->Tri[a][1]],result->dot[S->Tri[a][0]],vt1);
    subtract3f(result->dot[S->Tri[a][2]],result->dot[S->Tri[a][0]],vt2);
    cross_product3f(vt1,vt2,vt);
    if(dot_product3f(vt,result->dot[S->Tri[a][0]])<0.0) { /* if wrong, then interchange */
      it=S->Tri[a][2]; S->Tri[a][2]=S->Tri[a][1]; S->Tri[a][1]=it;
	 }
  }    

  for(a=0;a<S->NTri;a++)
	 {
		area = (float)(SphericalAngle(S,S->Tri[a][0],S->Tri[a][1],S->Tri[a][2]) +
		  SphericalAngle(S,S->Tri[a][1],S->Tri[a][0],S->Tri[a][2]) +
				  SphericalAngle(S,S->Tri[a][2],S->Tri[a][0],S->Tri[a][1]) - cPI);
		/* multiply by r^2 to get area */
		sumArea+=area;
		area/=3.0;
		result->area[S->Tri[a][0]] += area;
		result->area[S->Tri[a][1]] += area;
		result->area[S->Tri[a][2]] += area;
	 }

  if(fabs(sumArea - (4*cPI))>0.001) {
    printf(" MakeDotSphere: sumArea: %8.6f which is %8.6f Pi\n",sumArea,sumArea/cPI);
	 ErrFatal(G,"MakeDotSphere","Area of sphere does not sum to 4*pi!\n");
  }

  
  for(a=0;a<S->NTri;a++)
	 TriFlag[a]=false;

  nStrip=0;
  nVertTot=0;
  s=result->StripLen;
  q=result->Sequence;

  /* tesselate the sphere in a semi-efficient fashion...this could definitely be improved */

  flag = true;
  while(flag) {
	 flag=false;
	 a=0;
	 while(a<S->NTri)
		{
		  if(!TriFlag[a]) {
			 flag=true;
			 
			 TriFlag[a]=true;
			 *(q++)=S->Tri[a][0];
			 *(q++)=S->Tri[a][1];
			 *(q++)=S->Tri[a][2];
			 n=3;
			 
			 b=0;
			 while(b<S->NTri)
				{
				  if(!TriFlag[b]) {
					 if(((S->Tri[b][0]==q[-2])&&(S->Tri[b][1]==q[-1]))||
						 ((S->Tri[b][0]==q[-1])&&(S->Tri[b][1]==q[-2]))) {
						*(q++)=S->Tri[b][2];
						TriFlag[b]=true;
						b=0;
						n++;
					 } else if(((S->Tri[b][0]==q[-2])&&(S->Tri[b][2]==q[-1]))||
						 ((S->Tri[b][0]==q[-1])&&(S->Tri[b][2]==q[-2]))) {
						*(q++)=S->Tri[b][1];
						TriFlag[b]=true;
						b=0;
						n++;
					 } else if(((S->Tri[b][2]==q[-2])&&(S->Tri[b][1]==q[-1]))||
						 ((S->Tri[b][2]==q[-1])&&(S->Tri[b][1]==q[-2]))) {
						*(q++)=S->Tri[b][0];
						TriFlag[b]=true;
						b=0;
						n++;
					 } 
				  }
				  b++;
				}
			 if(n==3) {
				q[-3]=S->Tri[a][1];
				q[-2]=S->Tri[a][2];
				q[-1]=S->Tri[a][0];
				
				b=0;
				while(b<S->NTri)
				  {
					 if(!TriFlag[b]) {
						if(((S->Tri[b][0]==q[-2])&&(S->Tri[b][1]==q[-1]))||
							((S->Tri[b][0]==q[-1])&&(S->Tri[b][1]==q[-2]))) {
						  *(q++)=S->Tri[b][2];
						  TriFlag[b]=true;
						  b=0;
						  n++;
						} else if(((S->Tri[b][0]==q[-2])&&(S->Tri[b][2]==q[-1]))||
									 ((S->Tri[b][0]==q[-1])&&(S->Tri[b][2]==q[-2]))) {
						  *(q++)=S->Tri[b][1];
						  TriFlag[b]=true;
						  b=0;
						  n++;
						} else if(((S->Tri[b][2]==q[-2])&&(S->Tri[b][1]==q[-1]))||
									 ((S->Tri[b][2]==q[-1])&&(S->Tri[b][1]==q[-2]))) {
						  *(q++)=S->Tri[b][0];
						  TriFlag[b]=true;
						  b=0;
						  n++;
						} 
					 }
					 b++;
				  }
			 }
			 if(n==3) {
				q[-3]=S->Tri[a][2];
				q[-2]=S->Tri[a][0];
				q[-1]=S->Tri[a][1];
				b=0;
				while(b<S->NTri)
				  {
					 if(!TriFlag[b]) {
						if(((S->Tri[b][0]==q[-2])&&(S->Tri[b][1]==q[-1]))||
							((S->Tri[b][0]==q[-1])&&(S->Tri[b][1]==q[-2]))) {
						  *(q++)=S->Tri[b][2];
						  TriFlag[b]=true;
						  b=0;
						  n++;
						} else if(((S->Tri[b][0]==q[-2])&&(S->Tri[b][2]==q[-1]))||
									 ((S->Tri[b][0]==q[-1])&&(S->Tri[b][2]==q[-2]))) {
						  *(q++)=S->Tri[b][1];
						  TriFlag[b]=true;
						  b=0;
						  n++;
						} else if(((S->Tri[b][2]==q[-2])&&(S->Tri[b][1]==q[-1]))||
									 ((S->Tri[b][2]==q[-1])&&(S->Tri[b][1]==q[-2]))) {
						  *(q++)=S->Tri[b][0];
						  TriFlag[b]=true;
						  b=0;
						  n++;
						} 
					 }
					 b++;
				  }
			 }
			 *(s++) = n;
			 nVertTot+=n;
			 nStrip++;
		  }
		  a++;
		}
  }
  mfree(S->Dot);
  mfree(S->EdgeRef);
  mfree(TriFlag);
  result->Tri = (int*)S->Tri;
  result->Tri = Realloc(result->Tri,int,S->NTri*3);
  result->NTri = S->NTri;
  result->StripLen = Realloc(result->StripLen,int,nStrip);
  result->Sequence = Realloc(result->Sequence,int,nVertTot);
  result->dot = Realloc(result->dot,Vector3f,S->NDot);
  result->area = Realloc(result->area,float,S->NDot);
  result->nDot = S->NDot;
  result->NStrip = nStrip;
  result->NVertTot = nVertTot;
  result->Mesh = NULL;
  result->NMesh = 0;
  if (!level) { /* provide mesh for S->Sphere[0] only...rest, to do.*/
    result->Mesh=(int*)mesh;
    result->NMesh=30;
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

  return(result);
}

#endif
