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

#include"Isosurf.h"
#include"Tetsurf.h"
#include"MemoryDebug.h"
#include"Err.h"
#include"Crystal.h"
#include"Vector.h"
#include"Feedback.h"

#define Trace_OFF

#ifndef true
#define true 1
#endif

#ifndef false
#define false 0
#endif


#define O3(field,P1,P2,P3,offs) Ffloat3(field,P1+offs[0],P2+offs[1],P3+offs[2])

#define O3Ptr(field,P1,P2,P3,offs) Ffloat3p(field,P1+offs[0],P2+offs[1],P3+offs[2])

#define O4(field,P1,P2,P3,P4,offs) Ffloat4(field,P1+offs[0],P2+offs[1],P3+offs[2],P4)

#define O4Ptr(field,P1,P2,P3,P4,offs) Ffloat4p(field,P1+offs[0],P2+offs[1],P3+offs[2],P4)

#define I3(field,P1,P2,P3) Fint3(field,P1,P2,P3)

#define I3Ptr(field,P1,P2,P3) Fint3p(field,P1,P2,P3)

#define I4(field,P1,P2,P3,P4) Fint4(field,P1,P2,P3,P4)

#define I4Ptr(field,P1,P2,P3,P4) Fint4p(field,P1,P2,P3,P4)


typedef struct {
  float Point[3];
  float Normal[3];
  int NormalFlag;
  int  Link;
} PointType;

typedef struct {
  PointType *p[3];
  float n[3];
  int done;
} TriangleType;

typedef struct {
  int link;
  int tri;
} PointLinkType;

static TriangleType *Tri = NULL;
static PointLinkType *PtLink = NULL;

#define EdgePtPtr(field,P2,P3,P4,P5) ((PointType*)Fvoid4p(field,P2,P3,P4,P5))

#define EdgePt(field,P2,P3,P4,P5) (*((PointType*)Fvoid4p(field,P2,P3,P4,P5)))

static	CField	*VertexCodes;
static	CField	*ActiveEdges;
static	CField   *Point;

static	int	AbsDim[3],CurDim[3],CurOff[3];
static	int	Max[3];
static	CField *Coord,*Data;
static	float	Level;
static   int   Edge[6020]; /* 6017 */
static   int   EdgeStart[256];
static   int   TotPrim;

int	TetsurfInit(void);
int	TetsurfAlloc(void);
void	TetsurfFree(void);
int	TetsurfCodeVertices(void);
void	TetsurfInterpolate2(float *pt,float *v0,float l0,float *v1,float l1);
void	TetsurfInterpolate4(float *pt,float *v0,float l0,float *v1,float l1
                          ,float l2,float l3);
void	TetsurfInterpolate8(float *pt,float *v0,float l0,float *v1,float l1,
                          float l2,float l3,float l4,
                          float l5,float l6,float l7);

int	TetsurfFindActiveBoxes(int mode,int *n_strip,int n_vert,
                             int **strip_l,float **vert,
                             MapType *voxelmap,float *a_vert,
                             float carvebuffer);
void	TetsurfCode(char *bits1,char *bits2);
int	TetsurfPoints(void);

#define TetsurfSubSize		50

static void copy3fn ( float *v1,float *v2)
{
  v2[0]=v1[0];
  v2[1]=v1[1];
  v2[2]=v1[2];
}

/*===========================================================================*/
static int ProcessTetrahedron(int *edge,int nv,int v0,int v1,int v2,int v3,
                              int e01,int e02,int e03,int e12,int e13,int e23,
                              int reflect)
{
  int bits = (v3<<3)+(v2<<2)+(v1<<1)+v0;
  if(reflect)
    bits = 0xF-bits;
  switch(bits) {
  case 0x0:
    break;
  case 0x1: 
    edge[nv++]=e01;
    edge[nv++]=e02;
    edge[nv++]=e03;
    break;
  case 0x2:
    edge[nv++]=e01;
    edge[nv++]=e13;
    edge[nv++]=e12;
    break;
  case 0x3:
    edge[nv++]=e13;
    edge[nv++]=e12;
    edge[nv++]=e02;
    edge[nv++]=e03;
    edge[nv++]=e13;
    edge[nv++]=e02;
    break;
  case 0x4:
    edge[nv++]=e12;
    edge[nv++]=e23;
    edge[nv++]=e02;
    break;
  case 0x5:
    edge[nv++]=e01;
    edge[nv++]=e12;
    edge[nv++]=e03;
    edge[nv++]=e12;
    edge[nv++]=e23;
    edge[nv++]=e03;
    break;
  case 0x6:
    edge[nv++]=e01;
    edge[nv++]=e13;
    edge[nv++]=e02;
    edge[nv++]=e13;
    edge[nv++]=e23;
    edge[nv++]=e02;
    break;
  case 0x7:
    edge[nv++]=e03;
    edge[nv++]=e13;
    edge[nv++]=e23;
    break;
  case 0x8:
    edge[nv++]=e03;
    edge[nv++]=e23;
    edge[nv++]=e13;
    break;    
  case 0x9:
    edge[nv++]=e13;
    edge[nv++]=e01;
    edge[nv++]=e02;
    edge[nv++]=e02;
    edge[nv++]=e23;
    edge[nv++]=e13;
    break;    
  case 0xA:
    edge[nv++]=e01;
    edge[nv++]=e03;
    edge[nv++]=e12;
    edge[nv++]=e03;
    edge[nv++]=e23;
    edge[nv++]=e12;
    break;    
  case 0xB:
    edge[nv++]=e23;
    edge[nv++]=e12;
    edge[nv++]=e02;
    break;    
  case 0xC:
    edge[nv++]=e13;
    edge[nv++]=e02;
    edge[nv++]=e12;
    edge[nv++]=e03;
    edge[nv++]=e02;
    edge[nv++]=e13;
    break;    
  case 0xD:
    edge[nv++]=e01;
    edge[nv++]=e12;
    edge[nv++]=e13;
    break;    
  case 0xE:
    edge[nv++]=e01;
    edge[nv++]=e03;
    edge[nv++]=e02;
    break;    
  case 0xF:
    break;
  }
  return nv;
}
/*===========================================================================*/
int	TetsurfInit(void)
{
   /* there are six tetrahedrons in each cube, and there are 
      sixteen different types of tetrahedrons*/
   
   /* one will encounter 2^8 different kinds of cubes. */
   
   /* bits -> tetrahedral vertices 
         zyx
      0: 000
      1: 001
      2: 010
      3: 011
      4: 100
      5: 101
      6: 110
      7: 111
   */
   
   /* edges of the tetrahedron 

    vertex
      3210
    0:0011
    1:0101
    2:0110
    3:1001
    4:1010
    5:1100
   */

   /* edges of the cube 
      
        76543210
      0:00000011 edge 000-001
      1:00000101 edge 000-010
      2:00001001 edge 000-011
      3:00010001 edge 000-100
      4:00100001 edge 000-101
      5:01000001 edge 000-110
      6:10000001 edge 000-111
      7:00001010 edge 001-011
      8:00100010 edge 001-101
      9:10000010 edge 001-111
      0xA:00001100 edge 010-011
      0xB:01000100 edge 010-110
      0xC:10000100 edge 010-111
      0xD:00110000 edge 100-101
      0xE:01010000 edge 100-110
      0xF:10010000 edge 100-111
      0x10:10001000 edge 011-111
      0x11:10100000 edge 101-111
      0x12:11000000 edge 110-111
   */

#define cE_000_001 0x00
#define cE_000_010 0x01
#define cE_000_011 0x02
#define cE_000_100 0x03
#define cE_000_101 0x04
#define cE_000_110 0x05
#define cE_000_111 0x06
#define cE_001_011 0x07
#define cE_001_101 0x08
#define cE_001_111 0x09
#define cE_010_011 0x0A
#define cE_010_110 0x0B
#define cE_010_111 0x0C
#define cE_100_101 0x0D
#define cE_100_110 0x0E
#define cE_100_111 0x0F
#define cE_011_111 0x10
#define cE_101_111 0x11
#define cE_110_111 0x12

#define cM_000_001 0x00001
#define cM_000_010 0x00002
#define cM_000_011 0x00004
#define cM_000_100 0x00008
#define cM_000_101 0x00010
#define cM_000_110 0x00020
#define cM_000_111 0x00040
#define cM_001_011 0x00080
#define cM_001_101 0x00100
#define cM_001_111 0x00200
#define cM_010_011 0x00400
#define cM_010_110 0x00800
#define cM_010_111 0x01000
#define cM_100_101 0x02000
#define cM_100_110 0x04000
#define cM_100_111 0x08000
#define cM_011_111 0x10000
#define cM_101_111 0x20000
#define cM_110_111 0x40000

	int	ok=true;
	int	c;
   int   nv=1;
   int   last_nv;
   int v000,v100,v010,v110,v001,v101,v011,v111;

	VertexCodes=NULL;
	ActiveEdges=NULL;
	Point=NULL;
	
   last_nv = nv;
      
   for(c=0;c<256;c++) {
     
     v000=(c&0x01) ? 1 : 0;
     v001=(c&0x02) ? 1 : 0;
     v010=(c&0x04) ? 1 : 0;
     v011=(c&0x08) ? 1 : 0;
     v100=(c&0x10) ? 1 : 0;
     v101=(c&0x20) ? 1 : 0;
     v110=(c&0x40) ? 1 : 0;
     v111=(c&0x80) ? 1 : 0;

     /* tetrahedron 0: 000, 001, 011, 111 */
     nv=ProcessTetrahedron(Edge,nv,v000,v001,v011,v111,
                        cE_000_001,
                        cE_000_011,
                        cE_000_111,
                        cE_001_011,
                        cE_001_111,
                        cE_011_111,0);
     /* tetrahedron 1: 000, 001, 101, 111 */
     nv=ProcessTetrahedron(Edge,nv,v000,v001,v101,v111,
                        cE_000_001,
                        cE_000_101,
                        cE_000_111,
                        cE_001_101,
                        cE_001_111,
                        cE_101_111,1);

     /* tetrahedron 2: 000, 010, 011, 111 */
     nv=ProcessTetrahedron(Edge,nv,v000,v010,v011,v111,
                        cE_000_010,
                        cE_000_011,
                        cE_000_111,
                        cE_010_011,
                        cE_010_111,
                        cE_011_111,1);
     /* tetrahedron 3: 000, 010, 110, 111 */
     nv=ProcessTetrahedron(Edge,nv,v000,v010,v110,v111,
                        cE_000_010,
                        cE_000_110,
                        cE_000_111,
                        cE_010_110,
                        cE_010_111,
                        cE_110_111,0);

     /* tetrahedron 4: 000, 100, 101, 111 */
     nv=ProcessTetrahedron(Edge,nv,v000,v100,v101,v111,
                        cE_000_100,
                        cE_000_101,
                        cE_000_111,
                        cE_100_101,
                        cE_100_111,
                        cE_101_111,0);

     /* tetrahedron 5: 000, 100, 110, 111 */
     nv=ProcessTetrahedron(Edge,nv,v000,v100,v110,v111,
                        cE_000_100,
                        cE_000_110,
                        cE_000_111,
                        cE_100_110,
                        cE_100_111,
                        cE_110_111,1);

     Edge[nv++]=(-1); /* sentinel */
     EdgeStart[c] = last_nv;
     last_nv = nv;
   }
   /*   printf("%d\n",nv);
        for(c=1;c<nv;c++) {
        if(Edge[c]<0) {
        printf("\n");
        } else {
        printf("%02X ",Edge[c]);
        }
        }
   */
	return(ok);
}

/*===========================================================================*/
void TetsurfGetRange(Isofield *field,CCrystal *cryst,float *mn,float *mx,int *range)
{
  float fmn[3],fmx[3];
  float rmn[3],rmx[3];
  float imn[3],imx[3];
  int a;
  transform33f3f(cryst->RealToFrac,mn,fmn);
  transform33f3f(cryst->RealToFrac,mx,fmx);
  for(a=0;a<3;a++) {
    rmn[a] = F4(field->points,0,0,0,a);
    rmx[a] = F4(field->points,field->dimensions[0]-1,field->dimensions[1]-1,
                field->dimensions[2]-1,a);
  }
  
  transform33f3f(cryst->RealToFrac,rmn,imn);
  transform33f3f(cryst->RealToFrac,rmx,imx);

  for(a=0;a<3;a++) {
    range[a] = (int)((field->dimensions[a]*(fmn[a]-imn[a])/(imx[a]-imn[a])));
    if(range[a]<0) range[a]=0;
    range[a+3] = (int)((field->dimensions[a]*(fmx[a]-imn[a])/(imx[a]-imn[a]))+0.999);
    if(range[a]>field->dimensions[a])
      range[a]=field->dimensions[a];
    if(range[a+3]>field->dimensions[a])
      range[a+3]=field->dimensions[a];
  }
}
  /*===========================================================================*/
int	TetsurfVolume(Isofield *field,float level,int **num,float **vert,
                    int *range,int mode,MapType *voxelmap,float *a_vert,
                    float carvebuffer)
{
	int	ok=true;
	int	Steps[3];
	int	c,i,j,k;
   int range_store[6];
   int n_strip = 0;
   int n_vert = 0;

   TotPrim=0;
   if(range) {
     for(c=0;c<3;c++)
       {
         AbsDim[c]=field->dimensions[c];
         CurDim[c]=TetsurfSubSize+1;
         Steps[c]=1+((range[3+c]-range[c])-1)/TetsurfSubSize;
       }     
   } else {
     range=range_store;
     for(c=0;c<3;c++)
       {
         range[c]=0;
         range[3+c]=field->dimensions[c];
         AbsDim[c]=field->dimensions[c];
         CurDim[c]=TetsurfSubSize+1;
         Steps[c]=1+(AbsDim[c]-1)/TetsurfSubSize;
       }
   }
   
   /*   for(c=0;c<3;c++) {
        printf("range %d %d %d\n",c,range[c],range[c+3]);
        printf("steps %d\n",Steps[c]);
        }
   */
     
	Coord=field->points;
	Data=field->data;
	Level=level;
	if(ok) ok=TetsurfAlloc();

	if(ok)
		{

		for(i=0;i<Steps[0];i++)
		for(j=0;j<Steps[1];j++)
		for(k=0;k<Steps[2];k++)
			{
			CurOff[0]=TetsurfSubSize*i;
			CurOff[1]=TetsurfSubSize*j;
			CurOff[2]=TetsurfSubSize*k;
         for(c=0;c<3;c++)
           CurOff[c]+=range[c];
			for(c=0;c<3;c++)
           {
             Max[c]=(range[3+c]-CurOff[c]);
             if(Max[c]>(TetsurfSubSize+1))
               Max[c]=(TetsurfSubSize+1);
           }
         /*         
         for(c=0;c<3;c++)
           printf(" TetsurfVolume: c: %i CurOff[c]: %i Max[c] %i\n",c,CurOff[c],Max[c]); 
         */

			if(ok) 
           {
             if(TetsurfCodeVertices())
               n_vert=TetsurfFindActiveBoxes(mode,&n_strip,n_vert,num,vert,
                                             voxelmap,a_vert,carvebuffer);
           }
         }
      TetsurfFree();
      }
   
   if(Feedback(FB_Isosurface,FB_Actions)) { 
     if(mode!=2) {
       printf(" TetsurfVolume: Surface generated using %d vertices.\n",n_vert); 
     } else {
       printf(" TetsurfVolume: Surface generated using %d triangles.\n",TotPrim); 
     }
   }

   /* sentinel strip (0 length) */

   VLACheck(*num,int,n_strip);
   (*num)[n_strip]=0;
   (n_strip)++;
   
   /* shrinks sizes for more efficient RAM usage */

   VLASize(*vert,float,n_vert*3);
   VLASize(*num,int,n_strip);

	return(ok);
}
/*===========================================================================*/
int	TetsurfAlloc(void)
{
	int	ok=true;
	int dim4[4];
   int a;
   for(a=0;a<3;a++)
     dim4[a]=CurDim[a];
   dim4[3]=3;
   
	VertexCodes=FieldNew(CurDim,3,sizeof(int),cFieldInt);
	ErrChkPtr(VertexCodes);
	ActiveEdges=FieldNew(CurDim,3,sizeof(int),cFieldInt);
	ErrChkPtr(ActiveEdges);


   
   dim4[3]=7; /* seven different ways now... */
	Point=FieldNew(dim4,4,sizeof(PointType),cFieldOther);
	ErrChkPtr(Point);

   Tri = VLAlloc(TriangleType,50000);
   PtLink = VLAlloc(PointLinkType,50000);

	if(!(VertexCodes&&ActiveEdges&&Point))
		{
		TetsurfFree();
		ok=false;
		}
#ifdef Trace
	printf(" TetsurfAlloc: ok: %i\n",ok);
#endif
	return(ok);
}
/*===========================================================================*/
void	TetsurfFree(void)
{
  if(Tri)
    {
      VLAFreeP(Tri);
    }
  if(PtLink);
    {
      VLAFreeP(PtLink);
    }
  if(VertexCodes) 
    {
      FieldFree(VertexCodes);
      VertexCodes=NULL;
    }
  if(ActiveEdges)
    {
      FieldFree(ActiveEdges);
      ActiveEdges=NULL;
    }
  if(Point)
    {
      FieldFree(Point);
      Point=NULL;
    }
}
/*===========================================================================*/
void	TetsurfInterpolate2(float *pt,float *v0,float l0,float *v1,float l1)
{
  float	ratio;
  ratio=(Level-l0)/(l1-l0);
  pt[0]=v0[0]+(v1[0]-v0[0])*ratio;
  pt[1]=v0[1]+(v1[1]-v0[1])*ratio;
  pt[2]=v0[2]+(v1[2]-v0[2])*ratio;
}

/*===========================================================================*/
void	TetsurfInterpolate4(float *pt,float *v0,float l0,float *v1,float l1
                          ,float l2,float l3)

{
  float	ratio;
  float  v[3],l;
  average3f(v0,v1,v);
  l = (l0+l1+l2+l3)*0.25F;
  if(((l> Level)&&(l1>Level))||
     ((l<=Level)&&(l0>Level))) { /* l0 vs l */
    ratio=(Level-l0)/(l-l0);      
    pt[0]=v0[0]+(v[0]-v0[0])*ratio;
    pt[1]=v0[1]+(v[1]-v0[1])*ratio;
    pt[2]=v0[2]+(v[2]-v0[2])*ratio;
  } else {
    ratio=(Level-l1)/(l-l1);      
    pt[0]=v1[0]+(v[0]-v1[0])*ratio;
    pt[1]=v1[1]+(v[1]-v1[1])*ratio;
    pt[2]=v1[2]+(v[2]-v1[2])*ratio;
  }
}

/*===========================================================================*/
void	TetsurfInterpolate8(float *pt,float *v0,float l0,float *v1,float l1,
                          float l2,float l3,float l4,
                          float l5,float l6,float l7)
{
  float	ratio;
  float  v[3],l;
  average3f(v0,v1,v);
  l = (l0+l1+l2+l3+l4+l5+l6+l7)*0.125F;
  if(((l> Level)&&(l1>Level))||
     ((l<=Level)&&(l0>Level))) { /* l0 vs l */
    ratio=(Level-l0)/(l-l0);      
    pt[0]=v0[0]+(v[0]-v0[0])*ratio;
    pt[1]=v0[1]+(v[1]-v0[1])*ratio;
    pt[2]=v0[2]+(v[2]-v0[2])*ratio;
  } else {
    ratio=(Level-l1)/(l-l1);      
    pt[0]=v1[0]+(v[0]-v1[0])*ratio;
    pt[1]=v1[1]+(v[1]-v1[1])*ratio;
    pt[2]=v1[2]+(v[2]-v1[2])*ratio;
  }
}
/*===========================================================================*/
int	TetsurfFindActiveBoxes(int mode,int *n_strip,int n_vert,
                             int **strip_l,float **vert,
                             MapType *voxelmap,float *a_vert,float carvebuffer)
{
  int	a,b,c,i,j,k,h,l;
#ifdef Trace
	int	ECount=0;
#endif
	int i000,i001,i010,i011,i100,i101,i110,i111;
   float *c000,*c001,*c010,*c011,*c100,*c101,*c110,*c111;
   float d000,d001,d010,d011,d100,d101,d110,d111;
   int active;
   int n_active=0;
   int n_start=0;
   PointType *e[19],*p0,*p1,*p2;
   int code;
   int eidx;
   int idx;
   TriangleType *tt;
   int n_tri = 0;
   int n_link = 1;

   FieldZero(Point); /* sets initial links to zero */
   FieldZero(ActiveEdges);

   n_start = n_vert;
	for(i=0;i<(Max[0]-1);i++)
	for(j=0;j<(Max[1]-1);j++)
	for(k=0;k<(Max[2]-1);k++)
     {		
       active=0;

       i000=I3(VertexCodes,i  ,j  ,k  );
       i001=I3(VertexCodes,i  ,j  ,k+1);
       i010=I3(VertexCodes,i  ,j+1,k  );
       i011=I3(VertexCodes,i  ,j+1,k+1);
       i100=I3(VertexCodes,i+1,j  ,k  );
       i101=I3(VertexCodes,i+1,j  ,k+1);
       i110=I3(VertexCodes,i+1,j+1,k  );
       i111=I3(VertexCodes,i+1,j+1,k+1);

       if((i000!=i001)||
          (i001!=i010)||
          (i010!=i011)||
          (i011!=i100)||
          (i100!=i101)||
          (i101!=i110)||
          (i110!=i111)) { /* this is an active box */

         c000=O4Ptr(Coord,i  ,j  ,k  ,0,CurOff);
         c001=O4Ptr(Coord,i  ,j  ,k+1,0,CurOff);
         c010=O4Ptr(Coord,i  ,j+1,k  ,0,CurOff);
         c011=O4Ptr(Coord,i  ,j+1,k+1,0,CurOff);
         c100=O4Ptr(Coord,i+1,j  ,k  ,0,CurOff);
         c101=O4Ptr(Coord,i+1,j  ,k+1,0,CurOff);
         c110=O4Ptr(Coord,i+1,j+1,k  ,0,CurOff);
         c111=O4Ptr(Coord,i+1,j+1,k+1,0,CurOff);

         d000=O3(Data ,i  ,j  ,k  ,CurOff);
         d001=O3(Data ,i  ,j  ,k+1,CurOff);
         d010=O3(Data ,i  ,j+1,k  ,CurOff);
         d011=O3(Data ,i  ,j+1,k+1,CurOff);
         d100=O3(Data ,i+1,j  ,k  ,CurOff);
         d101=O3(Data ,i+1,j  ,k+1,CurOff);
         d110=O3(Data ,i+1,j+1,k  ,CurOff);
         d111=O3(Data ,i+1,j+1,k+1,CurOff);
         
         e[cE_000_001]=EdgePtPtr(Point,i  ,j  ,k  ,cE_000_001);
         e[cE_000_010]=EdgePtPtr(Point,i  ,j  ,k  ,cE_000_010);
         e[cE_000_011]=EdgePtPtr(Point,i  ,j  ,k  ,cE_000_011);
         e[cE_000_100]=EdgePtPtr(Point,i  ,j  ,k  ,cE_000_100);
         e[cE_000_101]=EdgePtPtr(Point,i  ,j  ,k  ,cE_000_101);

         e[cE_000_110]=EdgePtPtr(Point,i  ,j  ,k  ,cE_000_110);
         e[cE_000_111]=EdgePtPtr(Point,i  ,j  ,k  ,cE_000_111);
         e[cE_001_011]=EdgePtPtr(Point,i  ,j  ,k+1,cE_000_010);
         e[cE_001_101]=EdgePtPtr(Point,i  ,j  ,k+1,cE_000_100);
         e[cE_001_111]=EdgePtPtr(Point,i  ,j  ,k+1,cE_000_110);

         e[cE_010_011]=EdgePtPtr(Point,i  ,j+1,k  ,cE_000_001);
         e[cE_010_110]=EdgePtPtr(Point,i  ,j+1,k  ,cE_000_100);
         e[cE_010_111]=EdgePtPtr(Point,i  ,j+1,k  ,cE_000_101);
         e[cE_100_101]=EdgePtPtr(Point,i+1,j  ,k  ,cE_000_001);
         e[cE_100_110]=EdgePtPtr(Point,i+1,j  ,k  ,cE_000_010);

         e[cE_100_111]=EdgePtPtr(Point,i+1,j  ,k  ,cE_000_011);
         e[cE_011_111]=EdgePtPtr(Point,i  ,j+1,k+1,cE_000_100);
         e[cE_101_111]=EdgePtPtr(Point,i+1,j  ,k+1,cE_000_010);
         e[cE_110_111]=EdgePtPtr(Point,i+1,j+1,k  ,cE_000_001);

         /* Generate interpolated coordinates for all active edges */
         
         if(i000!=i001) {
           if(!(I3(ActiveEdges,i,j,k)&cM_000_001)) {
             I3(ActiveEdges,i,j,k)|=cM_000_001;
             TetsurfInterpolate2(e[cE_000_001]->Point,c000,d000,c001,d001);
           }
           active|=cM_000_001;
         }
         if(i000!=i010) {
           if(!(I3(ActiveEdges,i,j,k)&cM_000_010)) {
             I3(ActiveEdges,i,j,k)|=cM_000_010;
             TetsurfInterpolate2(e[cE_000_010]->Point,c000,d000,c010,d010);
             
           }
           active|=cM_000_010;
         }
         if(i000!=i011) {
           if(!(I3(ActiveEdges,i,j,k)&cM_000_011)) {
             I3(ActiveEdges,i,j,k)|=cM_000_011;
             TetsurfInterpolate4(e[cE_000_011]->Point,c000,d000,c011,d011,d001,d010);
           }
           active|=cM_000_011;
         }
         if(i000!=i100) {
           if(!(I3(ActiveEdges,i,j,k)&cM_000_100)) {
             I3(ActiveEdges,i,j,k)|=cM_000_100;
             TetsurfInterpolate2(e[cE_000_100]->Point,c000,d000,c100,d100);
             
           }
           active|=cM_000_100;
         }
         if(i000!=i101) {
           if(!(I3(ActiveEdges,i,j,k)&cM_000_101)) {
             I3(ActiveEdges,i,j,k)|=cM_000_101;
             TetsurfInterpolate4(e[cE_000_101]->Point,c000,d000,c101,d101,d100,d001);
             
           }
           active|=cM_000_101;
         }
         if(i000!=i110) {
           if(!(I3(ActiveEdges,i,j,k)&cM_000_110)) {
             I3(ActiveEdges,i,j,k)|=cM_000_110;
             TetsurfInterpolate4(e[cE_000_110]->Point,c000,d000,c110,d110,d100,d010);
           }
           active|=cM_000_110;
         }
         if(i000!=i111) {
           if(!(I3(ActiveEdges,i,j,k)&cM_000_111)) {
             I3(ActiveEdges,i,j,k)|=cM_000_111;
             TetsurfInterpolate8(e[cE_000_111]->Point,
                                 c000,d000,c111,d111,
                                 d001,d010,d011,d100,
                                 d101,d110);
           }
           active|=cM_000_111;
         }
         if(i001!=i011) {
           if(!(I3(ActiveEdges,i  ,j  ,k+1)&cM_000_010)) {
             I3(ActiveEdges,i  ,j  ,k+1)|=cM_000_010;
             TetsurfInterpolate2(e[cE_001_011]->Point,c001,d001,c011,d011);
           }
           active|=cM_001_011;
         }
         if(i001!=i101) {
           if(!(I3(ActiveEdges,i  ,j  ,k+1)&cM_000_100)) {
             I3(ActiveEdges,i  ,j  ,k+1)|=cM_000_100;
             TetsurfInterpolate2(e[cE_001_101]->Point,c001,d001,c101,d101);
           }
           active|=cM_001_101;
         }
         if(i001!=i111) {
           if(!(I3(ActiveEdges,i  ,j  ,k+1)&cM_000_110)) {
             I3(ActiveEdges,i  ,j  ,k+1)|=cM_000_110;
             TetsurfInterpolate4(e[cE_001_111]->Point,c001,d001,c111,d111,d101,d011);
           }
           active|=cM_001_111;
         }
         if(i010!=i011) {
           if(!(I3(ActiveEdges,i  ,j+1,k  )&cM_000_001)) {
             I3(ActiveEdges,i  ,j+1,k  )|=cM_000_001;
             TetsurfInterpolate2(e[cE_010_011]->Point,c010,d010,c011,d011);
           }
           active|=cM_010_011;
         }
         if(i010!=i110) {
           if(!(I3(ActiveEdges,i  ,j+1,k  )&cM_000_100)) {
             I3(ActiveEdges,i  ,j+1,k  )|=cM_000_100;
             TetsurfInterpolate2(e[cE_010_110]->Point,c010,d010,c110,d110);
           }
           active|=cM_010_110;
         }
         if(i010!=i111) {
           if(!(I3(ActiveEdges,i  ,j+1,k  )&cM_000_101)) {
             I3(ActiveEdges,i  ,j+1,k  )|=cM_000_101;
             TetsurfInterpolate4(e[cE_010_111]->Point,c010,d010,c111,d111,d110,d011);
           }
           active|=cM_010_111;
         }
         if(i100!=i101) {
           if(!(I3(ActiveEdges,i+1,j  ,k  )&cM_000_001)) {
             I3(ActiveEdges,i+1,j  ,k  )|=cM_000_001;
             TetsurfInterpolate2(e[cE_100_101]->Point,c100,d100,c101,d101);
           }
           active|=cM_100_101;
         }
         if(i100!=i110) {
           if(!(I3(ActiveEdges,i+1,j  ,k  )&cM_000_010)) {
             I3(ActiveEdges,i+1,j  ,k  )|=cM_000_010;
             TetsurfInterpolate2(e[cE_100_110]->Point,c100,d100,c110,d110);
           }
           active|=cM_100_110;
         }
         if(i100!=i111) {
           if(!(I3(ActiveEdges,i+1,j  ,k  )&cM_000_011)) {
             I3(ActiveEdges,i+1,j  ,k  )|=cM_000_011;
             TetsurfInterpolate4(e[cE_100_111]->Point,c100,d100,c111,d111,d101,d110);
           }
           active|=cM_100_111;
         }
         if(i011!=i111) {
           if(!(I3(ActiveEdges,i  ,j+1,k+1)&cM_000_100)) {
             I3(ActiveEdges,i  ,j+1,k+1)|=cM_000_100;
             TetsurfInterpolate2(e[cE_011_111]->Point,c011,d011,c111,d111);
           }
           active|=cM_011_111;
         }
         if(i101!=i111) {
           if(!(I3(ActiveEdges,i+1,j  ,k+1)&cM_000_010)) {
             I3(ActiveEdges,i+1,j  ,k+1)|=cM_000_010;
             TetsurfInterpolate2(e[cE_101_111]->Point,c101,d101,c111,d111);
           }
           active|=cM_101_111;
         }
         if(i110!=i111) {
           if(!(I3(ActiveEdges,i+1,j+1,k  )&cM_000_001)) {
             I3(ActiveEdges,i+1,j+1,k  )|=cM_000_001;
             TetsurfInterpolate2(e[cE_110_111]->Point,c110,d110,c111,d111);
           }
           active|=cM_110_111;
         }
                  
         if(active) {
           switch(mode) {
           case 2: 
             code=
               (i000 ? 0x01 : 0)|
               (i001 ? 0x02 : 0)|
               (i010 ? 0x04 : 0)|
               (i011 ? 0x08 : 0)|
               (i100 ? 0x10 : 0)|
               (i101 ? 0x20 : 0)|
               (i110 ? 0x40 : 0)|
               (i111 ? 0x80 : 0);
             eidx=EdgeStart[code];
             while(1) {
               idx=Edge[eidx];
               if(idx<0) break;

               /* assemble a triangle from these three points */

               VLACheck(Tri,TriangleType,n_tri);
               tt = Tri + n_tri;
               tt->p[0] = e[idx];
               tt->p[1] = e[Edge[eidx+1]];
               tt->p[2] = e[Edge[eidx+2]];

               VLACheck(PtLink,PointLinkType,n_link+3);

               /* link this triangle into the points */

               PtLink[n_link].tri = n_tri;
               PtLink[n_link].link = tt->p[0]->Link;
               tt->p[0]->Link=n_link;
               n_link++;

               PtLink[n_link].tri = n_tri;
               PtLink[n_link].link = tt->p[1]->Link;
               tt->p[1]->Link=n_link;
               n_link++;

               PtLink[n_link].tri = n_tri;
               PtLink[n_link].link = tt->p[2]->Link;
               tt->p[2]->Link=n_link;
               n_link++;
               n_tri++;
               eidx+=3;
             }
             break;
           case 1: /* lines */
             VLACheck(*vert,float,(n_vert*3)+200); 

             code=
               (i000 ? 0x01 : 0)|
               (i001 ? 0x02 : 0)|
               (i010 ? 0x04 : 0)|
               (i011 ? 0x08 : 0)|
               (i100 ? 0x10 : 0)|
               (i101 ? 0x20 : 0)|
               (i110 ? 0x40 : 0)|
               (i111 ? 0x80 : 0);
             eidx=EdgeStart[code];
             while(1) {
               idx=Edge[eidx];
               if(idx<0) break;
               copy3fn(e[idx]->Point,(*vert)+(n_vert*3));                                
               n_vert++;
               copy3fn(e[Edge[eidx+1]]->Point,(*vert)+(n_vert*3));                                
               n_vert++;
               copy3fn(e[Edge[eidx+1]]->Point,(*vert)+(n_vert*3));                                
               n_vert++;
               copy3fn(e[Edge[eidx+2]]->Point,(*vert)+(n_vert*3));               
               n_vert++;
               copy3fn(e[Edge[eidx+2]]->Point,(*vert)+(n_vert*3));               
               n_vert++;
               copy3fn(e[idx]->Point,(*vert)+(n_vert*3));                                
               n_vert++;
               eidx+=3;
             }
             break;
           case 0:
           default: /* dots */
             VLACheck(*vert,float,(n_vert*3)+200); 

             if(active&cM_000_001) {
               copy3fn(e[cE_000_001]->Point,(*vert)+(n_vert*3));
               n_vert+=1;
             }
             if(active&cM_000_010) {
               copy3fn(e[cE_000_010]->Point,(*vert)+(n_vert*3));
               n_vert+=1;
             }
             if(active&cM_000_011) {
               copy3fn(e[cE_000_011]->Point,(*vert)+(n_vert*3));
               n_vert+=1;
             }
             if(active&cM_000_100) {
               copy3fn(e[cE_000_100]->Point,(*vert)+(n_vert*3));
               n_vert+=1;
             }
             if(active&cM_000_101) {
               copy3fn(e[cE_000_101]->Point,(*vert)+(n_vert*3));
               n_vert+=1;
             }
             if(active&cM_000_110) {
               copy3fn(e[cE_000_110]->Point,(*vert)+(n_vert*3));
               n_vert+=1;
             }
             if(active&cM_000_111) {
               copy3fn(e[cE_000_111]->Point,(*vert)+(n_vert*3));
               n_vert+=1;
             }
             
             if(active&cM_001_011) {
               copy3fn(e[cE_001_011]->Point,(*vert)+(n_vert*3));
               n_vert+=1;
             }
             if(active&cM_001_101) {
               copy3fn(e[cE_001_101]->Point,(*vert)+(n_vert*3));
               n_vert+=1;
             }
             if(active&cM_001_111) {
               copy3fn(e[cE_001_111]->Point,(*vert)+(n_vert*3));
               n_vert+=1;
             }

             if(active&cM_010_011) {
               copy3fn(e[cE_010_011]->Point,(*vert)+(n_vert*3));
               n_vert+=1;
             }
             if(active&cM_010_110) {
               copy3fn(e[cE_010_011]->Point,(*vert)+(n_vert*3));
               n_vert+=1;
             }
             if(active&cM_010_111) {
               copy3fn(e[cE_010_111]->Point,(*vert)+(n_vert*3));
               n_vert+=1;
             }

             if(active&cM_100_101) {
               copy3fn(e[cE_100_101]->Point,(*vert)+(n_vert*3));
               n_vert+=1;
             }
             if(active&cM_100_110) {
               copy3fn(e[cE_100_110]->Point,(*vert)+(n_vert*3));
               n_vert+=1;
             }
             if(active&cM_100_111) {
               copy3fn(e[cE_100_111]->Point,(*vert)+(n_vert*3));
               n_vert+=1;
             }

             if(active&cM_011_111) {
               copy3fn(e[cE_011_111]->Point,(*vert)+(n_vert*3));
               n_vert+=1;
             }
             if(active&cM_101_111) {
               copy3fn(e[cE_101_111]->Point,(*vert)+(n_vert*3));
               n_vert+=1;
             }
             if(active&cM_110_111) {
               copy3fn(e[cE_101_111]->Point,(*vert)+(n_vert*3));
               n_vert+=1;
             }
             
             break;
           }
           n_active++;
         }
       }
     }
   
   switch(mode) {
   case 2:
     /* compute area-weighted normal */
     for(a=0;a<n_tri;a++) {
       float *v0,*v1,*v2;
       float vt1[3],vt2[3];

       tt = Tri + a;
       v0 = tt->p[0]->Point;
       v1 = tt->p[1]->Point;
       v2 = tt->p[2]->Point;
       tt->done=false; /* init */
       
       subtract3f(v0,v2,vt1);
       subtract3f(v1,v2,vt2);
       cross_product3f(vt2,vt1,tt->n);
     }
     /* compute normals at active points */
     for(a=0;a<n_tri;a++) {
       float v[3];

       tt = Tri + a;
       for(b=0;b<3;b++) {
         if(!tt->p[b]->NormalFlag) {
           zero3f(v);
           idx = tt->p[b]->Link;
           while(idx>0) {
             add3f(Tri[PtLink[idx].tri].n,v,v);
             idx = PtLink[idx].link;
           }
           normalize23f(v,tt->p[b]->Normal);
           tt->p[b]->NormalFlag=true;
         }
       }
     }
     /* now do an additional averaging cycle, no weighting */
     for(a=0;a<n_tri;a++) {
       tt = Tri + a;
       add3f(tt->p[0]->Normal,tt->p[1]->Normal,tt->n);
       add3f(tt->p[2]->Normal,tt->n,tt->n);
       normalize3f(tt->n);
       tt->p[0]->NormalFlag=false;
       tt->p[1]->NormalFlag=false;
       tt->p[2]->NormalFlag=false;
     }
     /* compute normals at active points */
     for(a=0;a<n_tri;a++) {
       float v[3];
       tt = Tri + a;
       for(b=0;b<3;b++) {
         if(!tt->p[b]->NormalFlag) {
           zero3f(v);
           idx = tt->p[b]->Link;
           while(idx>0) {
             add3f(Tri[PtLink[idx].tri].n,v,v);
             idx = PtLink[idx].link;
           }
           normalize23f(v,tt->p[b]->Normal);
           tt->p[b]->NormalFlag=true;
         }
       }
     }
     /* if we are carving, then exclude triangles outside region */
     if(voxelmap) {
       for(a=0;a<n_tri;a++) {
         float *v;
         tt = Tri + a;
         c=0;
         for(b=0;b<3;b++) {
           v=tt->p[b]->Point;
           MapLocus(voxelmap,v,&h,&k,&l);
           i=*(MapEStart(voxelmap,h,k,l));
           if(i) {
             j=voxelmap->EList[i++];
             while(j>=0) {
               if(within3f(a_vert+3*j,v,carvebuffer)) {
                 c++;
                 break;
               }
               j=voxelmap->EList[i++];
             }
           }
         }
         if(c<3) /* exclude this triangle from the surface */
           tt->done=true;
       }
     }

     /* now create triangle strips (not yet optimal) */
     for(a=0;a<n_tri;a++) {
       tt = Tri + a;
       n_start = n_vert;
       if(!tt->done) {

         VLACheck(*vert,float,(n_vert*3)+200); 

         /* switch order around to get "correct" triangles */
         copy3fn(tt->p[1]->Normal,(*vert)+(n_vert*3));
         n_vert++;
         copy3fn(tt->p[1]->Point,(*vert)+(n_vert*3));
         n_vert++;

         copy3fn(tt->p[0]->Normal,(*vert)+(n_vert*3));
         n_vert++;
         copy3fn(tt->p[0]->Point,(*vert)+(n_vert*3));
         n_vert++;
         
         copy3fn(tt->p[2]->Normal,(*vert)+(n_vert*3));
         n_vert++;
         copy3fn(tt->p[2]->Point,(*vert)+(n_vert*3));
         n_vert++;
         
         tt->done = true;

         p0 = tt->p[0];
         p1 = tt->p[2];
         
         while(1) {
           p2 = NULL;
           idx = p1->Link;
           while(idx>0) {
             tt=Tri+PtLink[idx].tri;

             if(!tt->done) {
               if((tt->p[0]==p0)&&(tt->p[1]==p1)) {
                 p2=tt->p[2];
                 break;
               }
               
               if((tt->p[1]==p0)&&(tt->p[2]==p1)) {
                 p2=tt->p[0];
                 break;
               }
               
               if((tt->p[2]==p0)&&(tt->p[0]==p1)) {
                 p2=tt->p[1];
                 break;
               }
               
               if((tt->p[1]==p0)&&(tt->p[0]==p1)) {
                 p2=tt->p[2];
                 break;
               }
               
               if((tt->p[2]==p0)&&(tt->p[1]==p1)) {
                 p2=tt->p[0];
                 break;
               }
               
               if((tt->p[0]==p0)&&(tt->p[2]==p1)) {
                 p2=tt->p[1];
                 break;
               }
             }
             idx = PtLink[idx].link;
           }

           if(!p2) break;
           tt->done=true;
           VLACheck(*vert,float,(n_vert*3)+200); 
           copy3fn(p2->Normal,(*vert)+(n_vert*3));
           n_vert++;
           copy3fn(p2->Point,(*vert)+(n_vert*3));
           n_vert++;
           p0=p1;
           p1=p2;
         }
       }
       if(n_vert>n_start) {
         VLACheck(*strip_l,int,*n_strip);
         (*strip_l)[*n_strip]=n_vert-n_start;
         (*n_strip)++;
       }
     }
     TotPrim+=n_tri;
     break;
     
   case 1:
     if(n_vert>n_start) {
       VLACheck(*strip_l,int,*n_strip);
       (*strip_l)[*n_strip]=n_vert-n_start;
       (*n_strip)++;
     }
     break;
   case 0: /* dots */
   default:
     if(n_vert>n_start) {
       VLACheck(*strip_l,int,*n_strip);
       (*strip_l)[*n_strip]=n_vert-n_start;
       (*n_strip)++;
     }
     break;
   }
   /*
     printf("n_strip %d\n",*n_strip);
     printf("n_active %d\n",n_active);
     printf("n_vert %d\n",n_vert);
     printf("mode %d\n",mode);
   */
	return(n_vert);
}
/*===========================================================================*/
int	TetsurfCodeVertices(void)
{
	int	i,j,k;
   int b0,b1;
   int flag1=false;
   int flag2=false;
   b0=1;
   if(Level<0.0F)
     b0=0;
   b1=1-b0;
   
	for(i=0;i<Max[0];i++)
	for(j=0;j<Max[1];j++)
	for(k=0;k<Max[2];k++)
		{
        if((O3(Data,i,j,k,CurOff)>Level)) {
          I3(VertexCodes,i,j,k)=b0;
          flag1=true;
        } else {
          I3(VertexCodes,i,j,k)=b1;
          flag2=true;
        }
		}
	return(flag1&&flag2);
}
