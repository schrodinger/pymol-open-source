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
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<values.h>

#include"Isosurf.h"
#include"MemoryDebug.h"
#include"Err.h"


#define Trace_OFF

#ifndef true
#define true 1
#endif

#ifndef false
#define false 0
#endif

#define O3(Name,P1,P2,P3,dims,offs) \
(*((Name)+((P1+offs[0])+((dims[0])*((P2+offs[1])+((dims[1])*(P3+offs[2])))))))

#define O3Ptr(Name,P1,P2,P3,dims,offs) \
((Name)+((P1+offs[0])+((dims[0])*((P2+offs[1])+((dims[1])*(P3+offs[2]))))))

#define O4(Name,P1,P2,P3,P4,dims,offs) \
(*((Name)+((P1+offs[0])+((dims[0])*((P2+offs[1])+((dims[1])*((P3+offs[2])+((dims[2])*(P4)))))))))

#define O4Ptr(Name,P1,P2,P3,P4,dims,offs) \
((Name)+((P1+offs[0])+((dims[0])*((P2+offs[1])+((dims[1])*((P3+offs[2])+((dims[2])*(P4))))))))

#define EdgePtPtr(Name,P2,P3,P4,P5,dims) \
((Name)+((P2)+((dims[0])*((P3)+((dims[1])*((P4)+((dims[2])*(P5))))))))

#define EdgePt(Name,P2,P3,P4,P5,dims) \
(*((Name)+((P2)+((dims[0])*((P3)+((dims[1])*((P4)+((dims[2])*(P5)))))))))

typedef struct	PointType {
	float		Point[3];
	int		NLink;
	struct PointType *(Link[4]);
	} PointType;

static	int	*VertexCodes;
static	int	*ActiveEdges;
static	PointType	*Point;
static	int	NLine;

static	int	AbsDim[3],CurDim[3],CurOff[3];
static	int	Max[3];
static	float	*Coord,*Data;
static	float	Level;
static	int	NVerts;
static	int	Code[256];

static MapType	*Map;
static int      *Num;
static int      NSeg;
static float    *Line;
static int      *Idx;

int	IsosurfInit(void);
int	IsosurfAlloc(void);
void	IsosurfFree(void);
int	IsosurfCurrent(void);
int	IsosurfCodeVertices(void);
void	IsosurfInterpolate(float *v1,float *l1,float *v2,float *l2,float *pt);
int	IsosurfFindActiveEdges(void);
int	IsosurfFindLines(void);
int	IsosurfDrawLines(void);
void	IsosurfSendLines(int Continue);
void	IsosurfCode(char *bits1,char *bits2);
void	IsosurfDotLine(float *v1, float *v2);

#define IsosurfSubSize		50

/*===========================================================================*/
/*===========================================================================*/
Isofield *IsosurfFieldAlloc(int *dims)
{
  Isofield *result;

  result=mmalloc(sizeof(Isofield));
  ErrChkPtr(result);
  result->data = mmalloc(sizeof(float)*dims[0]*dims[1]*dims[2]);
  ErrChkPtr(result->data);
  result->points = mmalloc(sizeof(float)*3*dims[0]*dims[1]*dims[2]);
  ErrChkPtr(result->points);
  result->dimensions[0]=dims[0];
  result->dimensions[1]=dims[1];
  result->dimensions[2]=dims[2];
  return(result);
}
/*===========================================================================*/
/*===========================================================================*/
void IsosurfFieldFree(Isofield *field)
{
  mfree(field->points);
  mfree(field->data);
  mfree(field);
}
/*===========================================================================*/
/*===========================================================================*/
void	IsosurfCode(char *bits1,char *bits2)
{
	int	c;
	int	b;
	int	sum1,sum2;
	
	c=0;
	while(bits1[c])
		c++;
	c--;
	sum1=0;
	b=1;
	while(c>=0)
		{
		if(bits1[c]=='1')
			sum1=sum1+b;
		b=b+b;
		c--;
		}

	c=0;
	while(bits2[c])
		c++;
	c--;
	sum2=0;
	b=1;
	while(c>=0)
		{
		if(bits2[c]=='1')
			sum2=sum2+b;
		b=b+b;
		c--;
		}

	Code[sum1]=sum2;
#ifdef Trace
	printf("IsosurfCode: %s (%i) -> %s (%i)\n",bits1,sum1,bits2, sum2);
#endif
}
/*===========================================================================*/
/*===========================================================================*/
int	IsosurfInit(void)
{
	int	ok=true;
	int	c;
	VertexCodes=NULL;
	ActiveEdges=NULL;
	Point=NULL;
	Line=NULL;
	
	for(c=0;c<255;c++)
	  Code[c]=-1;
		
/*___  
 | / |
 |/  |
 |___|
 32
*/
	IsosurfCode("10000010","100000");			
	IsosurfCode("01000001","100000");

/*___  
 | \ |
 |  \|
 |___|
 16
*/
	IsosurfCode("10010000","010000");
	IsosurfCode("01100000","010000");

/*___  
 |   |
 |  /|
 |_/_|
 8
*/
	IsosurfCode("00101000","001000");
	IsosurfCode("00010100","001000");

/*___  
 |   |
 |\  |
 |_\_|
 4
*/
	IsosurfCode("00001001","000100");
	IsosurfCode("00000110","000100");


/*___  
 | \ |
 |\ \|
 |_\_|
 16+4=20
*/

	IsosurfCode("01101001","010100");

/*___  
 | / |
 |/ /|
 |_/_|
 32+8=40
*/
	IsosurfCode("10010110","101000");


/*___  
 | | |
 | | |
 |_|_|
 2
*/
	IsosurfCode("10001000","000010");
	IsosurfCode("01000100","000010");

/*___  
 |   |
 |---|
 |___|
 1
*/
	IsosurfCode("00100010","000001");
	IsosurfCode("00010001","000001");
	
	return(ok);
}

/*===========================================================================*/
/*===========================================================================*/
int	IsosurfVolume(Isofield *field,float level,int **num,float **vert,MapType *map,int *idx)
{
	int	ok=true;
	int	Steps[3];
	int	c,i,j,k;
	int	percomp;
	int	x,y,z;

	Map = map;
	Num = *num;
	Line = *vert;
	Idx = idx;

	for(c=0;c<3;c++)
		{
		AbsDim[c]=field->dimensions[c];
		CurDim[c]=IsosurfSubSize+1;
		Steps[c]=(AbsDim[c]-2)/IsosurfSubSize+1;
		}
	Coord=field->points;
	Data=field->data;
	Level=level;
	NVerts=AbsDim[0]*AbsDim[1]*AbsDim[2];
	if(ok) ok=IsosurfAlloc();

	NLine=0;
	NSeg=0;
	VLACheck(Num,int,NSeg);
	Num[NSeg]=NLine;

	if(ok)
		{

		for(i=0;i<Steps[0];i++)
		for(j=0;j<Steps[1];j++)
		for(k=0;k<Steps[2];k++)
			{
			percomp=100*(k+(Steps[1]*(j+Steps[0]*i)))/(Steps[0]*Steps[1]*Steps[2]);
			CurOff[0]=IsosurfSubSize*i;
			CurOff[1]=IsosurfSubSize*j;
			CurOff[2]=IsosurfSubSize*k;
			for(c=0;c<3;c++)
				{
				Max[c]=AbsDim[c]-CurOff[c];
				if(Max[c]>(IsosurfSubSize+1))
					Max[c]=(IsosurfSubSize+1);
				}
			if(!(i||j||k))
				{
				for(x=0;x<Max[0];x++)
				for(y=0;y<Max[1];y++)
				for(z=0;z<Max[2];z++)
				for(c=0;c<3;c++)
					EdgePt(Point,x,y,z,c,CurDim).NLink=0;
				}

#ifdef Trace
for(c=0;c<3;c++)
	printf(" IsosurfVolume: c: %i CurOff[c]: %i Max[c] %i\n",c,CurOff[c],Max[c]); 
#endif
			if(ok) ok=IsosurfCurrent();
			}
		IsosurfFree();
		}

	Num[NSeg]=0; 

	printf(" IsosurfVolume: Surface generated using %d lines.\n",NLine); 

	*vert = Line;
	*num = Num;
	return(ok);
}
/*===========================================================================*/
/*===========================================================================*/
int	IsosurfAlloc(void)
{
	int	ok=true;
	int	NPnts;
	
	NPnts=CurDim[0]*CurDim[1]*CurDim[2];
	
	VertexCodes=(int*)mmalloc(NPnts*sizeof(int));
	ErrChkPtr(VertexCodes);
	ActiveEdges=(int*)mmalloc(NPnts*sizeof(int)*3);
	ErrChkPtr(ActiveEdges);
	Point=(PointType*)mmalloc(NPnts*sizeof(PointType)*3);
	ErrChkPtr(Point);
	if(!(VertexCodes&&ActiveEdges&&Point))
		{
		IsosurfFree();
		ok=false;
		}
#ifdef Trace
	printf(" IsosurfAlloc: ok: %i\n",ok);
#endif
	return(ok);
}
/*===========================================================================*/
/*===========================================================================*/

void	IsosurfFree(void)
{
	if(VertexCodes) 
		{
		mfree(VertexCodes);
		VertexCodes=NULL;
		}
	if(ActiveEdges)
		{
		mfree(ActiveEdges);
		ActiveEdges=NULL;
		}
	if(Point)
		{
		mfree(Point);
		Point=NULL;
		}
}
/*===========================================================================*/
/*===========================================================================*/
int	IsosurfCurrent(void)
{
	int	ok=true;
	if(IsosurfCodeVertices())
		{
		if(ok) ok=IsosurfFindActiveEdges();
		if(ok) ok=IsosurfFindLines();
		if(ok) ok=IsosurfDrawLines();
		}
	return(ok);
}
/*===========================================================================*/
/*===========================================================================*/
int	IsosurfDrawLines(void)
{
	int	c,i,j,k;
	float	*a,*b;
	int	ok=true;
	PointType	*Cur,*Start,*Next;
	int	MaxLinks,MaxL,Cnt;
	int	NLink;
#ifdef Trace
	int	LCount=0;
#endif	

	for(i=0;i<Max[0];i++)
	for(j=0;j<Max[1];j++)
	for(k=0;k<Max[2];k++)
	for(c=0;c<3;c++)
		{
		Start=EdgePtPtr(Point,i,j,k,c,CurDim);
		while(Start->NLink)
			{
			Cur=Start;
			VLACheck(Line,float,NLine*3+2);
			a=Line+(NLine*3);
			b=Cur->Point;
			*(a++)=*(b++);
			*(a++)=*(b++);
			*a=*b;
			NLine++;
			
			while(Cur)
				{
				if(Cur->NLink)
					{
					Cur->NLink--;
					NLink=Cur->NLink;
	/* Choose point which has most links */
					MaxL=NLink;
					MaxLinks=Cur->Link[MaxL]->NLink;
					Cnt=MaxL-1;
					while(Cnt>=0)
						{
						if((Cur->Link[Cnt]->NLink)>MaxLinks)
							{
							MaxL=Cnt;
							MaxLinks=Cur->Link[Cnt]->NLink;
							}
						Cnt--;
						}
					Next=Cur->Link[MaxL];
					if(MaxL!=NLink)
					  Cur->Link[MaxL]=Cur->Link[NLink];
	/* Remove double link */
					Next->NLink--;
					NLink=Next->NLink;
					Cnt=NLink;
					while(Cnt>=0)
						{
						if(Next->Link[Cnt]==Cur)
							break;
						else
							Cnt--;
						}
					if(Cnt>=0)
						{
						  if(Cnt!=NLink)
							 Next->Link[Cnt]=Next->Link[NLink];
						}
	#ifdef Trace
					else
						printf(" error: IsosurfDrawLines:  can't find double link\n");
	#endif
	
					Cur=Next;
					VLACheck(Line,float,NLine*3+2);
					a=Line+(NLine*3);
					b=Cur->Point;
					*(a++)=*(b++);
					*(a++)=*(b++);
					*a=*b;
					NLine++;
					}
				else
					{
#ifdef Trace
					LCount++;
#endif
					Cur=NULL;
					Num[NSeg]=NLine-Num[NSeg];
					NSeg++;
					VLACheck(Num,int,NSeg);
					Num[NSeg]=NLine;
					}
				}
			}
		}
#ifdef Trace
	printf(" DrawLineCount: %i\n",LCount);
#endif

	return(ok);
}
/*===========================================================================*/
/*===========================================================================*/
int	IsosurfFindLines(void)
{
	int	i,j,k,ip1,jp1,kp1;
	int	ok=true;
	int	index,cod;
	int	Max0m1,Max1m1,Max2m1;
	
#ifdef Trace
	int	LCount=0;
#endif

	PointType	*p1,*p2;
	
	Max0m1=Max[0]-1;
	Max1m1=Max[1]-1;
	Max2m1=Max[2]-1;
	for(i=0;i<Max[0];i++)
	for(j=0;j<Max[1];j++)
	for(k=0;k<Max[2];k++)
		{
		ip1=i+1;
		jp1=j+1;
		kp1=k+1;
		if((j<Max1m1)&&(k<Max2m1))
			{
			index=F4(ActiveEdges,i,j,k,1,CurDim)<<2;
			index=(index+F4(ActiveEdges,i,jp1,k,2,CurDim))<<2;
			index=(index+F4(ActiveEdges,i,j,kp1,1,CurDim))<<2;
			index=index+F4(ActiveEdges,i,j,k,2,CurDim);
			if(index)
				{
				cod=Code[index];
	#ifdef Trace
				if(index&&(cod<0))
					printf("IsosurfFindLines: bad index: %i \n",index);
	#endif
				while(cod>0)
					{
					p1=NULL;
					p2=NULL;
					switch(cod)
						{
						case 40:
						case 32:
							cod=cod-32;
							p1=EdgePtPtr(Point,i,j,k,1,CurDim);
							p2=EdgePtPtr(Point,i,j,k,2,CurDim);
							break;
						case 20:
						case 16:
							cod=cod-16;
							p1=EdgePtPtr(Point,i,j,k,1,CurDim);
							p2=EdgePtPtr(Point,i,jp1,k,2,CurDim);
							break;
						case 8:
							cod=cod-8;
							p1=EdgePtPtr(Point,i,j,kp1,1,CurDim);
							p2=EdgePtPtr(Point,i,jp1,k,2,CurDim);
							break;
						case 4:
							cod=cod-4;
							p1=EdgePtPtr(Point,i,j,kp1,1,CurDim);
							p2=EdgePtPtr(Point,i,j,k,2,CurDim);
							break;
						case 2:
							cod=cod-2;
							p1=EdgePtPtr(Point,i,j,k,1,CurDim);
							p2=EdgePtPtr(Point,i,j,kp1,1,CurDim);
							break;
						case 1:
							cod=cod-1;
							p1=EdgePtPtr(Point,i,j,k,2,CurDim);
							p2=EdgePtPtr(Point,i,jp1,k,2,CurDim);
							break;
						default:
							cod=0;
							p1=NULL;
							p2=NULL;
							break;
						}
					if(p1&&p2)
						{
							p1->Link[p1->NLink]=p2;
							p1->NLink++;
							p2->Link[p2->NLink]=p1;
							p2->NLink++;
	#ifdef Trace
	LCount++;
	#endif
						}
					}
				}
			}
		if((i<Max0m1)&&(j<Max1m1))
			{
			index=F4(ActiveEdges,i,j,k,0,CurDim)<<2;
			index=(index+F4(ActiveEdges,ip1,j,k,1,CurDim))<<2;
			index=(index+F4(ActiveEdges,i,jp1,k,0,CurDim))<<2;
			index=index+F4(ActiveEdges,i,j,k,1,CurDim);
			if(index)
				{
				cod=Code[index];
	#ifdef Trace
				if(index&&(cod<0))
					printf("IsosurfFindLines: bad index: %i \n",index);
	#endif
				while(cod>0)
					{
					switch(cod)
						{
						case 40:
						case 32:
							cod=cod-32;
							p1=EdgePtPtr(Point,i,j,k,0,CurDim);
							p2=EdgePtPtr(Point,i,j,k,1,CurDim);
							break;
						case 20:
						case 16:
							cod=cod-16;
							p1=EdgePtPtr(Point,i,j,k,0,CurDim);
							p2=EdgePtPtr(Point,ip1,j,k,1,CurDim);
							break;
						case 8:
							cod=cod-8;
							p1=EdgePtPtr(Point,i,jp1,k,0,CurDim);
							p2=EdgePtPtr(Point,ip1,j,k,1,CurDim);
							break;
						case 4:
							cod=cod-4;
							p1=EdgePtPtr(Point,i,jp1,k,0,CurDim);
							p2=EdgePtPtr(Point,i,j,k,1,CurDim);
							break;
						case 2:
							cod=cod-2;
							p1=EdgePtPtr(Point,i,j,k,0,CurDim);
							p2=EdgePtPtr(Point,i,jp1,k,0,CurDim);
							break;
						case 1:
							cod=cod-1;
							p1=EdgePtPtr(Point,i,j,k,1,CurDim);
							p2=EdgePtPtr(Point,ip1,j,k,1,CurDim);
							break;
						default:
							cod=0;
							p1=NULL;
							p2=NULL;
							break;
						}
					if(p1&&p2)
						{
							p1->Link[p1->NLink]=p2;
							p1->NLink++;
							p2->Link[p2->NLink]=p1;
							p2->NLink++;
	#ifdef Trace
	LCount++;
	#endif
						}
					}
				}
			}
		if((i<Max0m1)&&(k<Max2m1))
			{
			index=F4(ActiveEdges,i,j,k,0,CurDim)<<2;
			index=(index+F4(ActiveEdges,ip1,j,k,2,CurDim))<<2;
			index=(index+F4(ActiveEdges,i,j,kp1,0,CurDim))<<2;
			index=index+F4(ActiveEdges,i,j,k,2,CurDim);
			if(index)
				{
				cod=Code[index];
	#ifdef Trace
				if(index&&(cod<0))
					printf("IsosurfFindLines: bad index: %i \n",index);
	#endif
				while(cod>0)
					{
					switch(cod)
						{
						case 40:
						case 32:
							cod=cod-32;
							p1=EdgePtPtr(Point,i,j,k,0,CurDim);
							p2=EdgePtPtr(Point,i,j,k,2,CurDim);
							break;
						case 20:
						case 16:
							cod=cod-16;
							p1=EdgePtPtr(Point,i,j,k,0,CurDim);
							p2=EdgePtPtr(Point,ip1,j,k,2,CurDim);
							break;
						case 8:
							cod=cod-8;
							p1=EdgePtPtr(Point,i,j,k+1,0,CurDim);
							p2=EdgePtPtr(Point,ip1,j,k,2,CurDim);
							break;
						case 4:
							cod=cod-4;
							p1=EdgePtPtr(Point,i,j,kp1,0,CurDim);
							p2=EdgePtPtr(Point,i,j,k,2,CurDim);
							break;
						case 2:
							cod=cod-2;
							p1=EdgePtPtr(Point,i,j,k,0,CurDim);
							p2=EdgePtPtr(Point,i,j,kp1,0,CurDim);
							break;
						case 1:
							cod=cod-1;
							p1=EdgePtPtr(Point,i,j,k,2,CurDim);
							p2=EdgePtPtr(Point,ip1,j,k,2,CurDim);
							break;
						default:
							cod=0;
							p1=NULL;
							p2=NULL;
							break;
						}
					if(p1&&p2)
						{
							p1->Link[p1->NLink]=p2;
							p1->NLink++;
							p2->Link[p2->NLink]=p1;
							p2->NLink++;
	#ifdef Trace
	LCount++;
	#endif
						}
					}
				}
			}
		}
#ifdef Trace
printf(" IsosurfFindLines: %i lines found\n",LCount);
#endif
	return(ok);
}
/*===========================================================================*/
/*===========================================================================*/
void	IsosurfInterpolate(float *v1,float *l1,float *v2,float *l2,float *pt)
{
	int	c,offset;
	float	ratio;
	
	ratio=(Level-*l1)/(*l2-*l1);
	for(c=0;c<3;c++)
		{
		offset=c*NVerts;
		pt[c]=v1[offset]+(v2[offset]-v1[offset])*ratio;
		}
}
/*===========================================================================*/
/*===========================================================================*/
int	IsosurfFindActiveEdges(void)
{
	int	i,j,k;
	int	ok=true;
#ifdef Trace
	int	ECount=0;
#endif
	
	for(i=0;i<(Max[0]-1);i++)
	for(j=0;j<Max[1];j++)
	for(k=0;k<Max[2];k++)
		{		
		if((F3(VertexCodes,i,j,k,CurDim))&&(!F3(VertexCodes,i+1,j,k,CurDim)))
			{
#ifdef Trace
ECount++;
#endif
			F4(ActiveEdges,i,j,k,0,CurDim)=2;
			IsosurfInterpolate(
				O4Ptr(Coord,i,j,k,0,AbsDim,CurOff),
				O3Ptr(Data,i,j,k,AbsDim,CurOff),
				O4Ptr(Coord,i+1,j,k,0,AbsDim,CurOff),
				O3Ptr(Data,i+1,j,k,AbsDim,CurOff),
				&(EdgePt(Point,i,j,k,0,CurDim).Point[0]));
			}
		else if(!(F3(VertexCodes,i,j,k,CurDim))&&(F3(VertexCodes,i+1,j,k,CurDim)))
			{
#ifdef Trace
ECount++;
#endif
			F4(ActiveEdges,i,j,k,0,CurDim)=1;
			IsosurfInterpolate(
				O4Ptr(Coord,i,j,k,0,AbsDim,CurOff),
				O3Ptr(Data,i,j,k,AbsDim,CurOff),
				O4Ptr(Coord,i+1,j,k,0,AbsDim,CurOff),
				O3Ptr(Data,i+1,j,k,AbsDim,CurOff),
				&(EdgePt(Point,i,j,k,0,CurDim).Point[0]));
			}
		else
			F4(ActiveEdges,i,j,k,0,CurDim)=0;
		}

	for(i=0;i<Max[0];i++)
	for(j=0;j<(Max[1]-1);j++)
	for(k=0;k<Max[2];k++)
		{
		if((F3(VertexCodes,i,j,k,CurDim))&&(!F3(VertexCodes,i,j+1,k,CurDim)))
			{
#ifdef Trace
ECount++;
#endif
			F4(ActiveEdges,i,j,k,1,CurDim)=2;
			IsosurfInterpolate(
				O4Ptr(Coord,i,j,k,0,AbsDim,CurOff),
				O3Ptr(Data,i,j,k,AbsDim,CurOff),
				O4Ptr(Coord,i,j+1,k,0,AbsDim,CurOff),
				O3Ptr(Data,i,j+1,k,AbsDim,CurOff),
				&(EdgePt(Point,i,j,k,1,CurDim).Point[0]));
			}
		else if(!(F3(VertexCodes,i,j,k,CurDim))&&(F3(VertexCodes,i,j+1,k,CurDim)))
			{
#ifdef Trace
ECount++;
#endif
			F4(ActiveEdges,i,j,k,1,CurDim)=1;
			IsosurfInterpolate(
				O4Ptr(Coord,i,j,k,0,AbsDim,CurOff),
				O3Ptr(Data,i,j,k,AbsDim,CurOff),
				O4Ptr(Coord,i,j+1,k,0,AbsDim,CurOff),
				O3Ptr(Data,i,j+1,k,AbsDim,CurOff),
				&(EdgePt(Point,i,j,k,1,CurDim).Point[0]));
			}
		else
			F4(ActiveEdges,i,j,k,1,CurDim)=0;
		}

	for(i=0;i<Max[0];i++)
	for(j=0;j<Max[1];j++)
	for(k=0;k<(Max[2]-1);k++)
		{
		if((F3(VertexCodes,i,j,k,CurDim))&&(!F3(VertexCodes,i,j,k+1,CurDim)))
			{
#ifdef Trace
ECount++;
#endif
			F4(ActiveEdges,i,j,k,2,CurDim)=2;
			IsosurfInterpolate(
				O4Ptr(Coord,i,j,k,0,AbsDim,CurOff),
				O3Ptr(Data,i,j,k,AbsDim,CurOff),
				O4Ptr(Coord,i,j,k+1,0,AbsDim,CurOff),
				O3Ptr(Data,i,j,k+1,AbsDim,CurOff),
				&(EdgePt(Point,i,j,k,2,CurDim).Point[0]));
			}
		else if(!(F3(VertexCodes,i,j,k,CurDim))&&(F3(VertexCodes,i,j,k+1,CurDim)))
			{
#ifdef Trace
ECount++;
#endif
			F4(ActiveEdges,i,j,k,2,CurDim)=1;
			IsosurfInterpolate(
				O4Ptr(Coord,i,j,k,0,AbsDim,CurOff),
				O3Ptr(Data,i,j,k,AbsDim,CurOff),
				O4Ptr(Coord,i,j,k+1,0,AbsDim,CurOff),
				O3Ptr(Data,i,j,k+1,AbsDim,CurOff),
				&(EdgePt(Point,i,j,k,2,CurDim).Point[0]));
			}
		else
			F4(ActiveEdges,i,j,k,2,CurDim)=0;
		}
#ifdef Trace
printf(" IsosurfFindActiveEdges: %i active edges found\n",ECount);
#endif
			
	return(ok);
}
/*===========================================================================*/
/*===========================================================================*/
int	IsosurfCodeVertices(void)
{
	int	i,j,k;
	int	VCount=0;

	for(i=0;i<Max[0];i++)
	for(j=0;j<Max[1];j++)
	for(k=0;k<Max[2];k++)
		{
		if((O3(Data,i,j,k,AbsDim,CurOff)>Level))
			{
			F3(VertexCodes,i,j,k,CurDim)=1;
			VCount++;
			}
		else
			F3(VertexCodes,i,j,k,CurDim)=0;
		}
#ifdef Trace
printf(" IsosurfCodeVertices: %i of %i vertices above level\n",VCount,
CurDim[0]*CurDim[1]*CurDim[2]);
#endif
	return(VCount);
}
