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

typedef struct	PointType {
	float		Point[3];
	} PointType;

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


int	TetsurfInit(void);
int	TetsurfAlloc(void);
void	TetsurfFree(void);
int	TetsurfCodeVertices(void);
void	TetsurfInterpolate(float *v1,float *l1,float *v2,float *l2,float *pt);
int	TetsurfFindActiveBoxes(int mode,int *n_strip,int n_vert,int **strip_l,float **vert);
void	TetsurfCode(char *bits1,char *bits2);
int	TetsurfPoints(void);

#define TetsurfSubSize		50

/*===========================================================================*/
static int ProcessTetrahedron(int *edge,int nv,int v0,int v1,int v2,int v3,
                              int e01,int e02,int e03,int e12,int e13,int e23)
{
  int bits = (v0<<3)+(v1<<2)+(v2<<1)+v3;
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
                        cE_011_111);
     /* tetrahedron 1: 000, 001, 101, 111 */
     nv=ProcessTetrahedron(Edge,nv,v000,v001,v101,v111,
                        cE_000_001,
                        cE_000_101,
                        cE_000_111,
                        cE_001_101,
                        cE_001_111,
                        cE_101_111);
     /* tetrahedron 2: 000, 010, 011, 111 */
     nv=ProcessTetrahedron(Edge,nv,v000,v010,v011,v111,
                        cE_000_010,
                        cE_000_011,
                        cE_000_111,
                        cE_010_011,
                        cE_010_111,
                        cE_011_111);
     /* tetrahedron 3: 000, 010, 110, 111 */
     nv=ProcessTetrahedron(Edge,nv,v000,v010,v110,v111,
                        cE_000_010,
                        cE_000_110,
                        cE_000_111,
                        cE_010_110,
                        cE_010_111,
                        cE_110_111);
     /* tetrahedron 4: 000, 100, 101, 111 */
     nv=ProcessTetrahedron(Edge,nv,v000,v100,v101,v111,
                        cE_000_100,
                        cE_000_101,
                        cE_000_111,
                        cE_100_101,
                        cE_100_111,
                        cE_101_111);
     /* tetrahedron 5: 000, 100, 110, 111 */
     nv=ProcessTetrahedron(Edge,nv,v000,v100,v110,v111,
                        cE_000_100,
                        cE_000_110,
                        cE_000_111,
                        cE_100_110,
                        cE_100_111,
                        cE_110_111);
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
int	TetsurfVolume(Isofield *field,float level,int **num,float **vert,int *range,int mode)
{
	int	ok=true;
	int	Steps[3];
	int	c,i,j,k;
   int range_store[6];
   int n_strip = 0;
   int n_vert = 0;

   if(range) {
     for(c=0;c<3;c++)
       {
         AbsDim[c]=field->dimensions[c];
         CurDim[c]=TetsurfSubSize+1;
         Steps[c]=((range[3+c]-range[c])-2)/TetsurfSubSize+1;
       }     
   } else {
     range=range_store;
     for(c=0;c<3;c++)
       {
         range[c]=0;
         range[3+c]=field->dimensions[c];
         AbsDim[c]=field->dimensions[c];
         CurDim[c]=TetsurfSubSize+1;
         Steps[c]=(AbsDim[c]-2)/TetsurfSubSize+1;
       }
   }

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
				Max[c]=range[3+c]-CurOff[c];
				if(Max[c]>(TetsurfSubSize+1))
					Max[c]=(TetsurfSubSize+1);
				}
         
#ifdef Trace
         for(c=0;c<3;c++)
           printf(" TetsurfVolume: c: %i CurOff[c]: %i Max[c] %i\n",c,CurOff[c],Max[c]); 
#endif
         
			if(ok) 
           switch(mode) { 
           case 0: /* standard mode - want lines */
             if(TetsurfCodeVertices())
               n_vert=TetsurfFindActiveBoxes(mode,&n_strip,n_vert,num,vert);
             break;
           }
			}
		TetsurfFree();
		}
   
   if(Feedback(FB_Isosurface,FB_Actions)) { 
     printf(" TetsurfVolume: Surface generated using %d vertices.\n",n_vert); 
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
   
	VertexCodes=FieldNew(CurDim,3,sizeof(int));
	ErrChkPtr(VertexCodes);
	ActiveEdges=FieldNew(CurDim,3,sizeof(int));
	ErrChkPtr(ActiveEdges);
   FieldZero(ActiveEdges);

   dim4[3]=7; /* seven different ways now... */

	Point=FieldNew(dim4,4,sizeof(PointType));
	ErrChkPtr(Point);
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
void	TetsurfInterpolate(float *v1,float *l1,float *v2,float *l2,float *pt)
{
  register float	ratio;
  ratio=(Level-*l1)/(*l2-*l1);
  pt[0]=v1[0]+(v2[0]-v1[0])*ratio;
  pt[1]=v1[1]+(v2[1]-v1[1])*ratio;
  pt[2]=v1[2]+(v2[2]-v1[2])*ratio;
}
/*===========================================================================*/
int	TetsurfFindActiveBoxes(int mode,int *n_strip,int n_vert,int **strip_l,float **vert)
{
	int	i,j,k;
#ifdef Trace
	int	ECount=0;
#endif
	int i000,i001,i010,i011,i100,i101,i110,i111;
   int active;
   int n_active=0;
   int n_start=0;


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

       /* Generate interpolated coordinates for all active edges */

       if(i000!=i001) {
         if(!(I3(ActiveEdges,i,j,k)&cM_000_001)) {
           I3(ActiveEdges,i,j,k)|=cM_000_001;
           TetsurfInterpolate(O4Ptr(Coord,i  ,j  ,k  ,0,CurOff),
                              O3Ptr(Data ,i  ,j  ,k  ,CurOff),
                              O4Ptr(Coord,i  ,j  ,k+1,0,CurOff),
                              O3Ptr(Data ,i  ,j  ,k+1,CurOff),
                              &(EdgePt(Point,i,j,k,cE_000_001).Point[0])); 
         }
         active|=cM_000_001;
       }


       if(i000!=i010) {
         if(!(I3(ActiveEdges,i,j,k)&cM_000_010)) {
           I3(ActiveEdges,i,j,k)|=cM_000_010;
           TetsurfInterpolate(O4Ptr(Coord,i  ,j  ,k  ,0,CurOff),
                              O3Ptr(Data ,i  ,j  ,k  ,CurOff),
                              O4Ptr(Coord,i  ,j+1,k  ,0,CurOff),
                              O3Ptr(Data ,i  ,j+1,k  ,CurOff),
                              &(EdgePt(Point,i,j,k,cE_000_010).Point[0])); 
         }
         active|=cM_000_010;
       }

       /*

       if(i000!=i011) {
         if(!(I3(ActiveEdges,i,j,k)&cM_000_011)) {
           I3(ActiveEdges,i,j,k)|=cM_000_011;
           TetsurfInterpolate(O4Ptr(Coord,i  ,j  ,k  ,0,CurOff),
                              O3Ptr(Data ,i  ,j  ,k  ,CurOff),
                              O4Ptr(Coord,i  ,j+1,k+1,0,CurOff),
                              O3Ptr(Data ,i  ,j+1,k+1,CurOff),
                              &(EdgePt(Point,i,j,k,cE_000_011).Point[0])); 
         }
         active|=cM_000_011;
       }
       */

       if(i000!=i100) {
         if(!(I3(ActiveEdges,i,j,k)&cM_000_100)) {
           I3(ActiveEdges,i,j,k)|=cM_000_100;
           TetsurfInterpolate(O4Ptr(Coord,i  ,j  ,k  ,0,CurOff),
                              O3Ptr(Data ,i  ,j  ,k  ,CurOff),
                              O4Ptr(Coord,i+1,j  ,k  ,0,CurOff),
                              O3Ptr(Data ,i+1,j  ,k  ,CurOff),
                              &(EdgePt(Point,i,j,k,cE_000_100).Point[0])); 
         }
         active|=cM_000_100;
       }
       /*
       if(i000!=i101) {
         if(!(I3(ActiveEdges,i,j,k)&cM_000_101)) {
           I3(ActiveEdges,i,j,k)|=cM_000_101;
           TetsurfInterpolate(O4Ptr(Coord,i  ,j  ,k  ,0,CurOff),
                              O3Ptr(Data ,i  ,j  ,k  ,CurOff),
                              O4Ptr(Coord,i+1,j  ,k+1,0,CurOff),
                              O3Ptr(Data ,i+1,j  ,k+1,CurOff),
                              &(EdgePt(Point,i,j,k,cE_000_101).Point[0])); 
         }
         active|=cM_000_101;
       }
       if(i000!=i110) {
         if(!(I3(ActiveEdges,i,j,k)&cM_000_110)) {
           I3(ActiveEdges,i,j,k)|=cM_000_110;
           TetsurfInterpolate(O4Ptr(Coord,i  ,j  ,k  ,0,CurOff),
                              O3Ptr(Data ,i  ,j  ,k  ,CurOff),
                              O4Ptr(Coord,i+1,j+1,k  ,0,CurOff),
                              O3Ptr(Data ,i+1,j+1,k  ,CurOff),
                              &(EdgePt(Point,i,j,k,cE_000_110).Point[0])); 
         }
         active|=cM_000_110;
       }
       */

       if(i000!=i111) {
         if(!(I3(ActiveEdges,i,j,k)&cM_000_111)) {
           I3(ActiveEdges,i,j,k)|=cM_000_111;
           TetsurfInterpolate(O4Ptr(Coord,i  ,j  ,k  ,0,CurOff),
                              O3Ptr(Data ,i  ,j  ,k  ,CurOff),
                              O4Ptr(Coord,i+1,j+1,k+1,0,CurOff),
                              O3Ptr(Data ,i+1,j+1,k+1,CurOff),
                              &(EdgePt(Point,i,j,k,cE_000_111).Point[0])); 
         }
         active|=cM_000_111;
       }
       /*
       if(i001!=i011) {
         if(!(I3(ActiveEdges,i  ,j  ,k+1)&cM_000_010)) {
           I3(ActiveEdges,i  ,j  ,k+1)|=cM_000_010;
           TetsurfInterpolate(O4Ptr(Coord,i  ,j  ,k+1,0,CurOff),
                              O3Ptr(Data ,i  ,j  ,k+1,CurOff),
                              O4Ptr(Coord,i  ,j+1,k+1,0,CurOff),
                              O3Ptr(Data ,i  ,j+1,k+1,CurOff),
                              &(EdgePt(Point,i  ,j  ,k+1,cE_000_010).Point[0])); 
         }
         active|=cM_001_011;
       }
       if(i001!=i101) {
         if(!(I3(ActiveEdges,i  ,j  ,k+1)&cM_000_100)) {
           I3(ActiveEdges,i  ,j  ,k+1)|=cM_000_100;
           TetsurfInterpolate(O4Ptr(Coord,i  ,j  ,k+1,0,CurOff),
                              O3Ptr(Data ,i  ,j  ,k+1,CurOff),
                              O4Ptr(Coord,i+1,j  ,k+1,0,CurOff),
                              O3Ptr(Data ,i+1,j  ,k+1,CurOff),
                              &(EdgePt(Point,i  ,j  ,k+1,cE_000_100).Point[0])); 
         }
         active|=cM_001_101;
       }
       if(i001!=i111) {
         if(!(I3(ActiveEdges,i  ,j  ,k+1)&cM_000_110)) {
           I3(ActiveEdges,i  ,j  ,k+1)|=cM_000_110;
           TetsurfInterpolate(O4Ptr(Coord,i  ,j  ,k+1,0,CurOff),
                              O3Ptr(Data ,i  ,j  ,k+1,CurOff),
                              O4Ptr(Coord,i+1,j+1,k+1,0,CurOff),
                              O3Ptr(Data ,i+1,j+1,k+1,CurOff),
                              &(EdgePt(Point,i  ,j  ,k+1,cE_000_110).Point[0])); 
         }
         active|=cM_001_111;
       }
       if(i010!=i011) {
         if(!(I3(ActiveEdges,i  ,j+1,k  )&cM_000_001)) {
           I3(ActiveEdges,i  ,j+1,k  )|=cM_000_001;
           TetsurfInterpolate(O4Ptr(Coord,i  ,j+1,k  ,0,CurOff),
                              O3Ptr(Data ,i  ,j+1,k  ,CurOff),
                              O4Ptr(Coord,i  ,j+1,k+1,0,CurOff),
                              O3Ptr(Data ,i  ,j+1,k+1,CurOff),
                              &(EdgePt(Point,i  ,j+1,k  ,cE_000_001).Point[0])); 
         }
         active|=cM_010_011;
       }
       if(i010!=i110) {
         if(!(I3(ActiveEdges,i  ,j+1,k  )&cM_000_100)) {
           I3(ActiveEdges,i  ,j+1,k  )|=cM_000_100;
           TetsurfInterpolate(O4Ptr(Coord,i  ,j+1,k  ,0,CurOff),
                              O3Ptr(Data ,i  ,j+1,k  ,CurOff),
                              O4Ptr(Coord,i+1,j+1,k  ,0,CurOff),
                              O3Ptr(Data ,i+1,j+1,k  ,CurOff),
                              &(EdgePt(Point,i  ,j+1,k  ,cE_000_100).Point[0])); 
         }
         active|=cM_010_110;
       }
       if(i010!=i111) {
         if(!(I3(ActiveEdges,i  ,j+1,k  )&cM_000_101)) {
           I3(ActiveEdges,i  ,j+1,k  )|=cM_000_101;
           TetsurfInterpolate(O4Ptr(Coord,i  ,j+1,k  ,0,CurOff),
                              O3Ptr(Data ,i  ,j+1,k  ,CurOff),
                              O4Ptr(Coord,i+1,j+1,k+1,0,CurOff),
                              O3Ptr(Data ,i+1,j+1,k+1,CurOff),
                              &(EdgePt(Point,i  ,j+1,k  ,cE_000_101).Point[0])); 
         }
         active|=cM_010_111;
       }

       if(i100!=i101) {
         if(!(I3(ActiveEdges,i+1,j  ,k  )&cM_000_001)) {
           I3(ActiveEdges,i+1,j  ,k  )|=cM_000_001;
           TetsurfInterpolate(O4Ptr(Coord,i+1,j  ,k  ,0,CurOff),
                              O3Ptr(Data ,i+1,j  ,k  ,CurOff),
                              O4Ptr(Coord,i+1,j  ,k+1,0,CurOff),
                              O3Ptr(Data ,i+1,j  ,k+1,CurOff),
                              &(EdgePt(Point,i+1,j  ,k  ,cE_000_001).Point[0])); 
         }
         active|=cM_100_101;
       }
       if(i100!=i110) {
         if(!(I3(ActiveEdges,i+1,j  ,k  )&cM_000_010)) {
           I3(ActiveEdges,i+1,j  ,k  )|=cM_000_010;
           TetsurfInterpolate(O4Ptr(Coord,i+1,j  ,k  ,0,CurOff),
                              O3Ptr(Data ,i+1,j  ,k  ,CurOff),
                              O4Ptr(Coord,i+1,j+1,k  ,0,CurOff),
                              O3Ptr(Data ,i+1,j+1,k  ,CurOff),
                              &(EdgePt(Point,i+1,j  ,k  ,cE_000_010).Point[0])); 
         }
         active|=cM_100_110;
       }
       if(i100!=i111) {
         if(!(I3(ActiveEdges,i+1,j  ,k  )&cM_000_011)) {
           I3(ActiveEdges,i+1,j  ,k  )|=cM_000_011;
           TetsurfInterpolate(O4Ptr(Coord,i+1,j  ,k  ,0,CurOff),
                              O3Ptr(Data ,i+1,j  ,k  ,CurOff),
                              O4Ptr(Coord,i+1,j+1,k+1,0,CurOff),
                              O3Ptr(Data ,i+1,j+1,k+1,CurOff),
                              &(EdgePt(Point,i+1,j  ,k  ,cE_000_011).Point[0])); 
         }
         active|=cM_100_111;
       }
       if(i011!=i111) {
         if(!(I3(ActiveEdges,i  ,j+1,k+1)&cM_000_100)) {
           I3(ActiveEdges,i  ,j+1,k+1)|=cM_000_100;
           TetsurfInterpolate(O4Ptr(Coord,i  ,j+1,k+1,0,CurOff),
                              O3Ptr(Data ,i  ,j+1,k+1,CurOff),
                              O4Ptr(Coord,i+1,j+1,k+1,0,CurOff),
                              O3Ptr(Data ,i+1,j+1,k+1,CurOff),
                              &(EdgePt(Point,i  ,j+1,k+1,cE_000_100).Point[0])); 
         }
         active|=cM_011_111;
       }
       if(i101!=i111) {
         if(!(I3(ActiveEdges,i+1,j  ,k+1)&cM_000_010)) {
           I3(ActiveEdges,i+1,j  ,k+1)|=cM_000_010;
           TetsurfInterpolate(O4Ptr(Coord,i+1,j  ,k+1,0,CurOff),
                              O3Ptr(Data ,i+1,j  ,k+1,CurOff),
                              O4Ptr(Coord,i+1,j+1,k+1,0,CurOff),
                              O3Ptr(Data ,i+1,j+1,k+1,CurOff),
                              &(EdgePt(Point,i+1,j  ,k+1,cE_000_010).Point[0])); 
         }
         active|=cM_101_111;
       }
       if(i110!=i111) {
         if(!(I3(ActiveEdges,i+1,j+1,k  )&cM_000_001)) {
           I3(ActiveEdges,i+1,j+1,k  )|=cM_000_001;
           TetsurfInterpolate(O4Ptr(Coord,i+1,j+1,k  ,0,CurOff),
                              O3Ptr(Data ,i+1,j+1,k  ,CurOff),
                              O4Ptr(Coord,i+1,j+1,k+1,0,CurOff),
                              O3Ptr(Data ,i+1,j+1,k+1,CurOff),
                              &(EdgePt(Point,i+1,j+1,k  ,cE_000_001).Point[0])); 
         }
         active|=cM_110_111;
       }
       */

       if(active) {

         VLACheck(*vert,float,(n_vert*3)+50); 
         /* make sure we have a enough storage for even the most complex box */

         switch(mode) {
         case 2:
           break;
         case 1:
           break;
         case 0:
         default: /* dots */
           if(active&cM_000_001) {
             copy3f(EdgePt(Point,i,j,k,cE_000_001).Point,
                    (*vert)+(n_vert*3));
             n_vert+=1;
           }

           if(active&cM_000_010) {
             copy3f(EdgePt(Point,i,j,k,cE_000_010).Point,
                    (*vert)+(n_vert*3));
             n_vert+=1;
           }
           if(active&cM_000_011) {
             copy3f(EdgePt(Point,i,j,k,cE_000_011).Point,
                    (*vert)+(n_vert*3));
             n_vert+=1;
           }
           if(active&cM_000_100) {
             copy3f(EdgePt(Point,i,j,k,cE_000_100).Point,
                    (*vert)+(n_vert*3));
             n_vert+=1;
           }
           if(active&cM_000_101) {
             copy3f(EdgePt(Point,i,j,k,cE_000_101).Point,
                    (*vert)+(n_vert*3));
             n_vert+=1;
           }
           if(active&cM_000_110) {
             copy3f(EdgePt(Point,i,j,k,cE_000_110).Point,
                    (*vert)+(n_vert*3));
             n_vert+=1;
           }
           if(active&cM_000_111) {
             copy3f(EdgePt(Point,i,j,k,cE_000_111).Point,
                    (*vert)+(n_vert*3));
             n_vert+=1;
           }
           break;
         }
         n_active++;
       }
     }

   switch(mode) {
   case 2:
     break;
   case 1:
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
   printf("n_strip %d\n",*n_strip);
   printf("n_active %d\n",n_active);
   printf("n_vert %d\n",n_vert);
	return(n_vert);
}
/*===========================================================================*/
int	TetsurfCodeVertices(void)
{
	int	i,j,k;
	int	VCount=0;

	for(i=0;i<Max[0];i++)
	for(j=0;j<Max[1];j++)
	for(k=0;k<Max[2];k++)
		{
		if((O3(Data,i,j,k,CurOff)>Level))
			{
			I3(VertexCodes,i,j,k)=1;
			VCount++;
			}
		else
			I3(VertexCodes,i,j,k)=0;
		}
#ifdef Trace
printf(" TetsurfCodeVertices: %i of %i vertices above level\n",VCount,
CurDim[0]*CurDim[1]*CurDim[2]);
#endif
	return(VCount);
}
