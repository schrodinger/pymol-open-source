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

#include"os_std.h"

#include"Base.h"
#include"Sphere.h"
#include"Vector.h"
#include"Err.h"

#include"MemoryDebug.h"

float SphericalAngle(int d0,int d1,int d2);
void MakeVertex(int d1,int d2);

/* Twelve vertices of icosahedron on unit sphere */
#define tau 0.8506508084      /* t=(1+sqrt(5))/2, tau=t/sqrt(1+t^2)  */
#define one 0.5257311121      /* one=1/sqrt(1+t^2) , unit sphere     */

float start_points[13][3] = {
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

int icosahedron[21][3] = {
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

#define MAXDOT 650
#define MAXTRI 1300

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

void SphereInit(void)
{
  Sphere0 = MakeDotSphere(0);
  Sphere1 = MakeDotSphere(1);
  Sphere2 = MakeDotSphere(2);
  Sphere3 = MakeDotSphere(3);
}

void SphereDone(void)
{
  SphereFree(Sphere0);
  SphereFree(Sphere1);
  SphereFree(Sphere2);
  SphereFree(Sphere3);
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
  int a,b,c,h,k,l,curTri,n;
  float area,sumArea=0.0;
  int nStrip,*q,*s;
  int nVertTot;
  int flag;

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

  if(level>3)
	 level=3;
  
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
  result->dot = Alloc(DotRec,NDot);
  ErrChkPtr(result->dot);
  result->StripLen = Alloc(int,NTri*3);
  ErrChkPtr(result->StripLen);
  result->Sequence = Alloc(int,NTri*3);
  ErrChkPtr(result->Sequence);

  for(a=0;a<NDot;a++)
	 {
		for(c=0;c<3;c++)
		  result->dot[a].v[c]=*(Dot+(3*a+c));
		result->dot[a].area = 0.0;
	 }

  for(a=0;a<NTri;a++)
	 {
		area = (SphericalAngle(Tri[a][0],Tri[a][1],Tri[a][2]) +
		  SphericalAngle(Tri[a][1],Tri[a][0],Tri[a][2]) +
				  SphericalAngle(Tri[a][2],Tri[a][0],Tri[a][1]) - cPI);
		/* multiply by r^2 to get area */
		sumArea+=area;
		area/=3.0;
		result->dot[Tri[a][0]].area += area;
		result->dot[Tri[a][1]].area += area;
		result->dot[Tri[a][2]].area += area;
	 }

  if(fabs(sumArea - (4*cPI))>0.0001)
	 ErrFatal("MakeDotSphere","Area of sphere does not sum to 4*pi!\n");
  /*  printf(" MakeDotSphere: sumArea: %8.3f which is %8.3f Pi\n",sumArea,sumArea/cPI);*/


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
  result->nDot = NDot;
  result->NStrip = nStrip;
  result->NVertTot = nVertTot;

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
  mfree(I->dot);
  mfree(I->StripLen);
  mfree(I->Sequence);
  FreeP(I);
}



