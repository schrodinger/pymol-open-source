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

float SphericalAngle(int d0,int d1,int d2);
void MakeVertex(int d1,int d2);

/* Twelve vertices of icosahedron on unit sphere */
#define tau 0.8506508084F      /* t=(1+sqrt(5))/2, tau=t/sqrt(1+t^2)  */
#define one 0.5257311121F      /* one=1/sqrt(1+t^2) , unit sphere     */

static float start_points[13][3] = {
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

static int icosahedron[21][3] = {
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

int mesh[30][2] = {
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

#define MAXDOT 2600
#define MAXTRI 6200

typedef int EdgeCol[MAXDOT]; /* should move these into dynamic storage to save 3MB  mem */
typedef EdgeCol EdgeArray[MAXDOT]; 

typedef int Triangle[3];

static float *Dot;
static EdgeArray *EdgeRef;
static Triangle *Tri;

static int NDot,NTri;

SphereRec *Sphere0;
SphereRec *Sphere1;
SphereRec *Sphere2;
SphereRec *Sphere3;
SphereRec *Sphere4;

#define FAST_SPHERE_INIT


#ifdef FAST_SPHERE_INIT
#include"SphereData.h"

SphereRec sSphere0,sSphere1,sSphere2,sSphere3,sSphere4,sSphere5;


#endif


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

static void SphereDumpAll()
{
  FILE *f;
  f = fopen("SphereData.h","w");
  SphereDump(f,"Sphere0",Sphere0);
  SphereDump(f,"Sphere1",Sphere1);
  SphereDump(f,"Sphere2",Sphere2);
  SphereDump(f,"Sphere3",Sphere3);
  SphereDump(f,"Sphere4",Sphere4);
  fclose(f);
}


void SphereInit(void)
{
  
#ifdef FAST_SPHERE_INIT
  sSphere0.area = Sphere0_area;
  sSphere0.dot = Sphere0_dot;
  sSphere0.StripLen = Sphere0_StripLen;
  sSphere0.Sequence = Sphere0_Sequence;
  sSphere0.NStrip = Sphere0_NStrip;
  sSphere0.NVertTot = Sphere0_NVertTot;
  sSphere0.nDot = Sphere0_nDot;
  sSphere0.Tri = Sphere0_Tri;
  sSphere0.NTri = Sphere0_NTri;
  sSphere0.Mesh = (int*)mesh;
  sSphere0.NMesh = 30;

  sSphere1.area = Sphere1_area;
  sSphere1.dot = Sphere1_dot;
  sSphere1.StripLen = Sphere1_StripLen;
  sSphere1.Sequence = Sphere1_Sequence;
  sSphere1.NStrip = Sphere1_NStrip;
  sSphere1.NVertTot = Sphere1_NVertTot;
  sSphere1.nDot = Sphere1_nDot;
  sSphere1.Tri = Sphere1_Tri;
  sSphere1.NTri = Sphere1_NTri;
  sSphere1.Mesh = NULL;
  sSphere1.NMesh = 0;

  sSphere2.area = Sphere2_area;
  sSphere2.dot = Sphere2_dot;
  sSphere2.StripLen = Sphere2_StripLen;
  sSphere2.Sequence = Sphere2_Sequence;
  sSphere2.NStrip = Sphere2_NStrip;
  sSphere2.NVertTot = Sphere2_NVertTot;
  sSphere2.nDot = Sphere2_nDot;
  sSphere2.Tri = Sphere2_Tri;
  sSphere2.NTri = Sphere2_NTri;
  sSphere2.Mesh = NULL;
  sSphere2.NMesh = 0;

  sSphere3.area = Sphere3_area;
  sSphere3.dot = Sphere3_dot;
  sSphere3.StripLen = Sphere3_StripLen;
  sSphere3.Sequence = Sphere3_Sequence;
  sSphere3.NStrip = Sphere3_NStrip;
  sSphere3.NVertTot = Sphere3_NVertTot;
  sSphere3.nDot = Sphere3_nDot;
  sSphere3.Tri = Sphere3_Tri;
  sSphere3.NTri = Sphere3_NTri;
  sSphere3.Mesh = NULL;
  sSphere3.NMesh = 0;

  sSphere4.area = Sphere4_area;
  sSphere4.dot = Sphere4_dot;
  sSphere4.StripLen = Sphere4_StripLen;
  sSphere4.Sequence = Sphere4_Sequence;
  sSphere4.NStrip = Sphere4_NStrip;
  sSphere4.NVertTot = Sphere4_NVertTot;
  sSphere4.nDot = Sphere4_nDot;
  sSphere4.Tri = Sphere4_Tri;
  sSphere4.NTri = Sphere4_NTri;
  sSphere4.Mesh = NULL;
  sSphere4.NMesh = 0;

  Sphere0 = &sSphere0;
  Sphere1 = &sSphere1;
  Sphere2 = &sSphere2;
  Sphere3 = &sSphere3;
  Sphere4 = &sSphere4;
#else
  Sphere0 = MakeDotSphere(0);
  Sphere1 = MakeDotSphere(1);
  Sphere2 = MakeDotSphere(2);
  Sphere3 = MakeDotSphere(3);
  Sphere4 = MakeDotSphere(4);
  /*
  SphereDumpAll();
  */
  
#endif


}

void SphereDone(void)
{
#ifndef FAST_SPHERE_INIT
  SphereFree(Sphere0);
  SphereFree(Sphere1);
  SphereFree(Sphere2);
  SphereFree(Sphere3);
  SphereFree(Sphere4);
#endif

}

void MakeVertex(int d1,int d2)
{
  if((*EdgeRef)[d1][d2]<0)
	 {
		average3f(Dot+(3*d1),Dot+(3*d2),Dot+(3*NDot));
		(*EdgeRef)[d1][d2]=NDot;
		(*EdgeRef)[d2][d1]=NDot;
		normalize3f(Dot+(3*NDot));
		NDot++;
	 }
}

float SphericalAngle(int d0,int d1,int d2)
{
  Vector3f v1,v2,s1,s2;

  /* map vector onto surface of sphere and measure angle */
  subtract3f(Dot+(3*d1),Dot+(3*d0),v1);
  subtract3f(Dot+(3*d2),Dot+(3*d0),v2);
  
  remove_component3f(v1,Dot+(3*d0),s1);
  remove_component3f(v2,Dot+(3*d0),s2);
  return(get_angle3f(s1,s2));
  
}

SphereRec *MakeDotSphere(int level)
{
  SphereRec *result;
  int *TriFlag;
  int a,b,c,h,k,l,curTri,n,it;
  float area,sumArea=0.0;
  int nStrip,*q,*s;
  int nVertTot;
  int flag;
  float vt1[3],vt2[3],vt[3];

  Dot=(float*)mmalloc(sizeof(float)*3*MAXDOT);
  ErrChkPtr(Dot);
  EdgeRef=(EdgeArray*)mmalloc(sizeof(EdgeArray));
  ErrChkPtr(EdgeRef);
  Tri=Alloc(Triangle,MAXTRI);
  ErrChkPtr(Tri);
  TriFlag=Alloc(int,MAXTRI);
  ErrChkPtr(TriFlag);
  
  NDot = 12;
  for(a=0;a<NDot;a++)
	 {
		for(c=0;c<3;c++)
		  Dot[3*a+c]=start_points[a][c];
		normalize3f(Dot+(3*a));
	 }

  NTri = 20;
  for(a=0;a<NTri;a++)
	 for(c=0;c<3;c++)
		Tri[a][c]=icosahedron[a][c];

  for(a=0;a<MAXDOT;a++)
	 for(b=0;b<MAXDOT;b++)
		(*EdgeRef)[a][b]=-1;

  if(level>4)
	 level=4;
  
  for(c=0;c<level;c++)
	 {
		/* create new vertices */
		for(a=0;a<NTri;a++)
		  {
			 MakeVertex(Tri[a][0],Tri[a][1]);
			 MakeVertex(Tri[a][1],Tri[a][2]);
			 MakeVertex(Tri[a][0],Tri[a][2]);
		  }		
		/* create new triangles */
		curTri=NTri;
		for(a=0;a<curTri;a++)
		  {
			 h=Tri[a][0];
			 k=Tri[a][1];
			 l=Tri[a][2];

			 Tri[a][0]=h;
			 Tri[a][1]=(*EdgeRef)[h][k];
			 Tri[a][2]=(*EdgeRef)[h][l];
			 
			 Tri[NTri][0]=k;
			 Tri[NTri][1]=(*EdgeRef)[k][h];
			 Tri[NTri][2]=(*EdgeRef)[k][l];
			 NTri++;
			 
			 Tri[NTri][0]=l;
			 Tri[NTri][1]=(*EdgeRef)[l][h];
			 Tri[NTri][2]=(*EdgeRef)[l][k];
			 NTri++;
			 
			 Tri[NTri][0]=(*EdgeRef)[h][k];
			 Tri[NTri][1]=(*EdgeRef)[k][l];
			 Tri[NTri][2]=(*EdgeRef)[l][h];
			 NTri++;		
		}
		/*		printf( "MakeDotSphere: Level: %i  NTri: %i\n",c, NTri); */
	 }
  /*  printf(" MakeDotSphere: NDot %i NTri %i\n",NDot,NTri);*/
  result= Alloc(SphereRec,1);
  ErrChkPtr(result);
  result->dot = Alloc(Vector3f,NDot);
  ErrChkPtr(result->dot);
  result->area = Alloc(float,NDot);
  ErrChkPtr(result->area);
  result->StripLen = Alloc(int,NTri*3);
  ErrChkPtr(result->StripLen);
  result->Sequence = Alloc(int,NTri*3);
  ErrChkPtr(result->Sequence);

  for(a=0;a<NDot;a++)
	 {
		for(c=0;c<3;c++)
		  result->dot[a][c]=*(Dot+(3*a+c));
		result->area[a] = 0.0;
	 }

  /* fix normals so that v1-v0 x v2-v0 is the correct normal */

  for(a=0;a<NTri;a++) {
    subtract3f(result->dot[Tri[a][1]],result->dot[Tri[a][0]],vt1);
    subtract3f(result->dot[Tri[a][2]],result->dot[Tri[a][0]],vt2);
    cross_product3f(vt1,vt2,vt);
    if(dot_product3f(vt,result->dot[Tri[a][0]])<0.0) { /* if wrong, then interchange */
      it=Tri[a][2]; Tri[a][2]=Tri[a][1]; Tri[a][1]=it;
	 }
  }    

  for(a=0;a<NTri;a++)
	 {
		area = (SphericalAngle(Tri[a][0],Tri[a][1],Tri[a][2]) +
		  SphericalAngle(Tri[a][1],Tri[a][0],Tri[a][2]) +
				  SphericalAngle(Tri[a][2],Tri[a][0],Tri[a][1]) - cPI);
		/* multiply by r^2 to get area */
		sumArea+=area;
		area/=3.0;
		result->area[Tri[a][0]] += area;
		result->area[Tri[a][1]] += area;
		result->area[Tri[a][2]] += area;
	 }

  if(fabs(sumArea - (4*cPI))>0.001) {
    printf(" MakeDotSphere: sumArea: %8.6f which is %8.6f Pi\n",sumArea,sumArea/cPI);
	 ErrFatal("MakeDotSphere","Area of sphere does not sum to 4*pi!\n");
  }

  
  for(a=0;a<NTri;a++)
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
	 while(a<NTri)
		{
		  if(!TriFlag[a]) {
			 flag=true;
			 
			 TriFlag[a]=true;
			 *(q++)=Tri[a][0];
			 *(q++)=Tri[a][1];
			 *(q++)=Tri[a][2];
			 n=3;
			 
			 b=0;
			 while(b<NTri)
				{
				  if(!TriFlag[b]) {
					 if(((Tri[b][0]==q[-2])&&(Tri[b][1]==q[-1]))||
						 ((Tri[b][0]==q[-1])&&(Tri[b][1]==q[-2]))) {
						*(q++)=Tri[b][2];
						TriFlag[b]=true;
						b=0;
						n++;
					 } else if(((Tri[b][0]==q[-2])&&(Tri[b][2]==q[-1]))||
						 ((Tri[b][0]==q[-1])&&(Tri[b][2]==q[-2]))) {
						*(q++)=Tri[b][1];
						TriFlag[b]=true;
						b=0;
						n++;
					 } else if(((Tri[b][2]==q[-2])&&(Tri[b][1]==q[-1]))||
						 ((Tri[b][2]==q[-1])&&(Tri[b][1]==q[-2]))) {
						*(q++)=Tri[b][0];
						TriFlag[b]=true;
						b=0;
						n++;
					 } 
				  }
				  b++;
				}
			 if(n==3) {
				q[-3]=Tri[a][1];
				q[-2]=Tri[a][2];
				q[-1]=Tri[a][0];
				
				b=0;
				while(b<NTri)
				  {
					 if(!TriFlag[b]) {
						if(((Tri[b][0]==q[-2])&&(Tri[b][1]==q[-1]))||
							((Tri[b][0]==q[-1])&&(Tri[b][1]==q[-2]))) {
						  *(q++)=Tri[b][2];
						  TriFlag[b]=true;
						  b=0;
						  n++;
						} else if(((Tri[b][0]==q[-2])&&(Tri[b][2]==q[-1]))||
									 ((Tri[b][0]==q[-1])&&(Tri[b][2]==q[-2]))) {
						  *(q++)=Tri[b][1];
						  TriFlag[b]=true;
						  b=0;
						  n++;
						} else if(((Tri[b][2]==q[-2])&&(Tri[b][1]==q[-1]))||
									 ((Tri[b][2]==q[-1])&&(Tri[b][1]==q[-2]))) {
						  *(q++)=Tri[b][0];
						  TriFlag[b]=true;
						  b=0;
						  n++;
						} 
					 }
					 b++;
				  }
			 }
			 if(n==3) {
				q[-3]=Tri[a][2];
				q[-2]=Tri[a][0];
				q[-1]=Tri[a][1];
				b=0;
				while(b<NTri)
				  {
					 if(!TriFlag[b]) {
						if(((Tri[b][0]==q[-2])&&(Tri[b][1]==q[-1]))||
							((Tri[b][0]==q[-1])&&(Tri[b][1]==q[-2]))) {
						  *(q++)=Tri[b][2];
						  TriFlag[b]=true;
						  b=0;
						  n++;
						} else if(((Tri[b][0]==q[-2])&&(Tri[b][2]==q[-1]))||
									 ((Tri[b][0]==q[-1])&&(Tri[b][2]==q[-2]))) {
						  *(q++)=Tri[b][1];
						  TriFlag[b]=true;
						  b=0;
						  n++;
						} else if(((Tri[b][2]==q[-2])&&(Tri[b][1]==q[-1]))||
									 ((Tri[b][2]==q[-1])&&(Tri[b][1]==q[-2]))) {
						  *(q++)=Tri[b][0];
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
  mfree(Dot);
  mfree(EdgeRef);
  mfree(TriFlag);
  result->Tri = (int*)Tri;
  result->Tri = Realloc(result->Tri,int,NTri*3);
  result->NTri = NTri;
  result->StripLen = Realloc(result->StripLen,int,nStrip);
  result->Sequence = Realloc(result->Sequence,int,nVertTot);
  result->dot = Realloc(result->dot,Vector3f,NDot);
  result->area = Realloc(result->area,float,NDot);
  result->nDot = NDot;
  result->NStrip = nStrip;
  result->NVertTot = nVertTot;
  result->Mesh = NULL;
  result->NMesh = 0;
  if (!level) { /* provide mesh for Sphere0 only...rest, to do.*/
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

void SphereFree(SphereRec *I)
{
  /* NOTE: I->Mesh is not currently a pointer*/
  mfree(I->dot);
  mfree(I->area);
  mfree(I->StripLen);
  mfree(I->Sequence);
  mfree(I->Tri);
  FreeP(I);
}



