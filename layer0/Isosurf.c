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
#include"MemoryDebug.h"
#include"Err.h"
#include"Crystal.h"
#include"Vector.h"
#include"Feedback.h"
#include"PConv.h"

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
	int		NLink;
	struct PointType *(Link[4]);
	} PointType;

#define EdgePtPtr(field,P2,P3,P4,P5) ((PointType*)Fvoid4p(field,P2,P3,P4,P5))

#define EdgePt(field,P2,P3,P4,P5) (*((PointType*)Fvoid4p(field,P2,P3,P4,P5)))

static	CField	*VertexCodes;
static	CField	*ActiveEdges;
static	CField   *Point;
static	int	NLine;

static	int	AbsDim[3],CurDim[3],CurOff[3];
static	int	Max[3];
static	CField *Coord,*Data;
static	float	Level;
static	int	Code[256];

static int      *Num;
static int      NSeg;
static float    *Line;

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
int	IsosurfDrawPoints(void);
int	IsosurfPoints(void);

#define IsosurfSubSize		50

PyObject *IsosurfGetPyList(Isofield *I)
{
  PyObject *result=NULL;

  result = PyList_New(3);
  PyList_SetItem(result,0,PConvIntArrayToPyList(I->dimensions,3));
  PyList_SetItem(result,1,PConvFloatArrayToPyList((float*)I->data->data,
                                                  I->dimensions[0]*
                                                  I->dimensions[1]*
                                                  I->dimensions[2]));
  PyList_SetItem(result,2,PConvFloatArrayToPyList((float*)I->points->data,
                                                  I->dimensions[0]*
                                                  I->dimensions[1]*
                                                  I->dimensions[2]*
                                                  I->dimensions[3]));
  return(PConvAutoNone(result));
}

Isofield *IsosurfNewFromPyList(PyObject *list)
{
  int ok=true;
  Isofield *result = NULL;
  int dim4[4];
  if(ok) ok=(list!=NULL);
  if(ok) ok=PyList_Check(list);
  if(ok) ok=PConvPyListToIntArrayInPlace(PyList_GetItem(list,0),dim4,3);
  dim4[3] = 3;
  
  result=mmalloc(sizeof(Isofield));
  result->data = NULL;
  result->points = NULL;

  ErrChkPtr(result);
  result->data = FieldNew(dim4,3,sizeof(float));
  ErrChkPtr(result->data);
  if(ok) ok=PConvPyListToFloatArrayInPlace(PyList_GetItem(list,1),
                                           (float*)result->data->data,
                                           dim4[0]*dim4[1]*dim4[2]);

  result->points = FieldNew(dim4,4,sizeof(float));
  ErrChkPtr(result->points);
  if(ok) ok=PConvPyListToFloatArrayInPlace(PyList_GetItem(list,2),
                                           (float*)result->points->data,
                                           dim4[0]*dim4[1]*dim4[2]*dim4[3]);

  result->dimensions[0]=dim4[0];
  result->dimensions[1]=dim4[1];
  result->dimensions[2]=dim4[2];
  if(!ok) {
    if(result) {
      if(result->dimensions)
        mfree(result->dimensions);
      if(result->points)
        mfree(result->points);
      mfree(result);
    }
  }
  return(result);
}

/*===========================================================================*/
/*===========================================================================*/
Isofield *IsosurfFieldAlloc(int *dims)
{
  int dim4[4];
  int a;
  Isofield *result;

  for(a=0;a<3;a++)
    dim4[a]=dims[a];
  dim4[3] = 3;
  
  result=mmalloc(sizeof(Isofield));
  ErrChkPtr(result);
  result->data = FieldNew(dims,3,sizeof(float));
  ErrChkPtr(result->data);
  result->points = FieldNew(dim4,4,sizeof(float));
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
  FieldFree(field->points);
  FieldFree(field->data);
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
	
	for(c=0;c<256;c++)
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
void IsosurfGetRange(Isofield *field,CCrystal *cryst,float *mn,float *mx,int *range)
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
int	IsosurfVolume(Isofield *field,float level,int **num,float **vert,int *range,int mode)
{
	int	ok=true;
	int	Steps[3];
	int	c,i,j,k;
	int	x,y,z;
   int range_store[6];
	Num = *num;
	Line = *vert;

   
   if(range) {
     for(c=0;c<3;c++)
       {
         AbsDim[c]=field->dimensions[c];
         CurDim[c]=IsosurfSubSize+1;
         Steps[c]=((range[3+c]-range[c])-2)/IsosurfSubSize+1;
       }     
   } else {
     range=range_store;
     for(c=0;c<3;c++)
       {
         range[c]=0;
         range[3+c]=field->dimensions[c];
         AbsDim[c]=field->dimensions[c];
         CurDim[c]=IsosurfSubSize+1;
         Steps[c]=(AbsDim[c]-2)/IsosurfSubSize+1;
       }
   }

	Coord=field->points;
	Data=field->data;
	Level=level;
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
			CurOff[0]=IsosurfSubSize*i;
			CurOff[1]=IsosurfSubSize*j;
			CurOff[2]=IsosurfSubSize*k;
         for(c=0;c<3;c++)
           CurOff[c]+=range[c];
			for(c=0;c<3;c++)
				{
				Max[c]=range[3+c]-CurOff[c];
				if(Max[c]>(IsosurfSubSize+1))
					Max[c]=(IsosurfSubSize+1);
				}
			if(!(i||j||k))
				{
              for(x=0;x<Max[0];x++)
                for(y=0;y<Max[1];y++)
                  for(z=0;z<Max[2];z++)
                    for(c=0;c<3;c++)
                      EdgePt(Point,x,y,z,c).NLink=0;
				}
         
#ifdef Trace
         for(c=0;c<3;c++)
           printf(" IsosurfVolume: c: %i CurOff[c]: %i Max[c] %i\n",c,CurOff[c],Max[c]); 
#endif
         
			if(ok) 
           switch(mode) { 
           case 0: /* standard mode - want lines */
             ok=IsosurfCurrent();
             break;
           case 1: /* point mode - just want points on the isosurface */
             ok=IsosurfPoints();
             break;
           }
			}
		IsosurfFree();
		}
   
   Num[NSeg]=0;  /* important - must terminate the segment list */
   
   if(Feedback(FB_Isomesh,FB_Actions)) { 
     if(mode)
       printf(" IsosurfVolume: Surface generated using %d dots.\n",NLine); 
     else
       printf(" IsosurfVolume: Surface generated using %d lines.\n",NLine); 
   }

   /* shrinks sizes for more efficient RAM usage */

   VLASize(Line,float,NLine*3);
   VLASize(Num,int,NSeg+1);

	*vert = Line;
	*num = Num;
	return(ok);
}
/*===========================================================================*/
/*===========================================================================*/
int	IsosurfAlloc(void)
{
	int	ok=true;
	int dim4[4];
   int a;
   for(a=0;a<3;a++)
     dim4[a]=CurDim[a];
   dim4[3]=3;

	VertexCodes=FieldNew(CurDim,3,sizeof(int));
	ErrChkPtr(VertexCodes);
	ActiveEdges=FieldNew(dim4,4,sizeof(int));
	ErrChkPtr(ActiveEdges);
	Point=FieldNew(dim4,4,sizeof(PointType));
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
int	IsosurfPoints(void)
{
	int	ok=true;
	if(IsosurfCodeVertices())
		{
		if(ok) ok=IsosurfFindActiveEdges();
		if(ok) ok=IsosurfDrawPoints();
		}
	return(ok);
}
/*===========================================================================*/
int	IsosurfDrawPoints(void)
{
  float *a,*b;
  int	i,j,k;
  int	ok=true;
  
 	for(i=0;i<(Max[0]-1);i++)
	for(j=0;j<Max[1];j++)
	for(k=0;k<Max[2];k++)
		{		
		if((I3(VertexCodes,i,j,k))&&(!I3(VertexCodes,i+1,j,k)))
			{
			IsosurfInterpolate(
				O4Ptr(Coord,i,j,k,0,CurOff),
				O3Ptr(Data,i,j,k,CurOff),
				O4Ptr(Coord,i+1,j,k,0,CurOff),
				O3Ptr(Data,i+1,j,k,CurOff),
				&(EdgePt(Point,i,j,k,0).Point[0]));
         
         VLACheck(Line,float,NLine*3+2);
         a=Line+(NLine*3);
         b=&(EdgePt(Point,i,j,k,0).Point[0]);
         *(a++)=*(b++);
         *(a++)=*(b++);
         *a=*b;
         NLine++;
			}
		else if(!(I3(VertexCodes,i,j,k))&&(I3(VertexCodes,i+1,j,k)))
			{
			IsosurfInterpolate(
				O4Ptr(Coord,i,j,k,0,CurOff),
				O3Ptr(Data,i,j,k,CurOff),
				O4Ptr(Coord,i+1,j,k,0,CurOff),
				O3Ptr(Data,i+1,j,k,CurOff),
				&(EdgePt(Point,i,j,k,0).Point[0]));

         VLACheck(Line,float,NLine*3+2);
         a=Line+(NLine*3);
         b=&(EdgePt(Point,i,j,k,0).Point[0]);
         *(a++)=*(b++);
         *(a++)=*(b++);
         *a=*b;
         NLine++;
			}
		else
			I4(ActiveEdges,i,j,k,0)=0;
		}

	for(i=0;i<Max[0];i++)
	for(j=0;j<(Max[1]-1);j++)
	for(k=0;k<Max[2];k++)
		{
		if((I3(VertexCodes,i,j,k))&&(!I3(VertexCodes,i,j+1,k)))
			{
			I4(ActiveEdges,i,j,k,1)=2;
			IsosurfInterpolate(
				O4Ptr(Coord,i,j,k,0,CurOff),
				O3Ptr(Data,i,j,k,CurOff),
				O4Ptr(Coord,i,j+1,k,0,CurOff),
				O3Ptr(Data,i,j+1,k,CurOff),
				&(EdgePt(Point,i,j,k,1).Point[0]));

         VLACheck(Line,float,NLine*3+2);
         a=Line+(NLine*3);
         b=&(EdgePt(Point,i,j,k,1).Point[0]);
         *(a++)=*(b++);
         *(a++)=*(b++);
         *a=*b;
         NLine++;

			}
		else if(!(I3(VertexCodes,i,j,k))&&(I3(VertexCodes,i,j+1,k)))
			{
			IsosurfInterpolate(
				O4Ptr(Coord,i,j,k,0,CurOff),
				O3Ptr(Data,i,j,k,CurOff),
				O4Ptr(Coord,i,j+1,k,0,CurOff),
				O3Ptr(Data,i,j+1,k,CurOff),
				&(EdgePt(Point,i,j,k,1).Point[0]));

         VLACheck(Line,float,NLine*3+2);
         a=Line+(NLine*3);
         b=&(EdgePt(Point,i,j,k,1).Point[0]);
         *(a++)=*(b++);
         *(a++)=*(b++);
         *a=*b;
         NLine++;

			}
		}

	for(i=0;i<Max[0];i++)
	for(j=0;j<Max[1];j++)
	for(k=0;k<(Max[2]-1);k++)
		{
		if((I3(VertexCodes,i,j,k))&&(!I3(VertexCodes,i,j,k+1)))
			{
			IsosurfInterpolate(
				O4Ptr(Coord,i,j,k,0,CurOff),
				O3Ptr(Data,i,j,k,CurOff),
				O4Ptr(Coord,i,j,k+1,0,CurOff),
				O3Ptr(Data,i,j,k+1,CurOff),
				&(EdgePt(Point,i,j,k,2).Point[0]));

         VLACheck(Line,float,NLine*3+2);
         a=Line+(NLine*3);
         b=&(EdgePt(Point,i,j,k,2).Point[0]);
         *(a++)=*(b++);
         *(a++)=*(b++);
         *a=*b;
         NLine++;

			}
		else if(!(I3(VertexCodes,i,j,k))&&(I3(VertexCodes,i,j,k+1)))
			{
			IsosurfInterpolate(
				O4Ptr(Coord,i,j,k,0,CurOff),
				O3Ptr(Data,i,j,k,CurOff),
				O4Ptr(Coord,i,j,k+1,0,CurOff),
				O3Ptr(Data,i,j,k+1,CurOff),
				&(EdgePt(Point,i,j,k,2).Point[0]));

         VLACheck(Line,float,NLine*3+2);
         a=Line+(NLine*3);
         b=&(EdgePt(Point,i,j,k,2).Point[0]);
         *(a++)=*(b++);
         *(a++)=*(b++);
         *a=*b;
         NLine++;

			}
		}

   Num[NSeg]=NLine-Num[NSeg];
   NSeg++;
   Num[NSeg]=NLine;
	return(ok);
}
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
		Start=EdgePtPtr(Point,i,j,k,c);
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
			index=I4(ActiveEdges,i,j,k,1)<<2;
			index=(index+I4(ActiveEdges,i,jp1,k,2))<<2;
			index=(index+I4(ActiveEdges,i,j,kp1,1))<<2;
			index=index+I4(ActiveEdges,i,j,k,2);
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
							p1=EdgePtPtr(Point,i,j,k,1);
							p2=EdgePtPtr(Point,i,j,k,2);
							break;
						case 20:
						case 16:
							cod=cod-16;
							p1=EdgePtPtr(Point,i,j,k,1);
							p2=EdgePtPtr(Point,i,jp1,k,2);
							break;
						case 8:
							cod=cod-8;
							p1=EdgePtPtr(Point,i,j,kp1,1);
							p2=EdgePtPtr(Point,i,jp1,k,2);
							break;
						case 4:
							cod=cod-4;
							p1=EdgePtPtr(Point,i,j,kp1,1);
							p2=EdgePtPtr(Point,i,j,k,2);
							break;
						case 2:
							cod=cod-2;
							p1=EdgePtPtr(Point,i,j,k,1);
							p2=EdgePtPtr(Point,i,j,kp1,1);
							break;
						case 1:
							cod=cod-1;
							p1=EdgePtPtr(Point,i,j,k,2);
							p2=EdgePtPtr(Point,i,jp1,k,2);
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
			index=I4(ActiveEdges,i,j,k,0)<<2;
			index=(index+I4(ActiveEdges,ip1,j,k,1))<<2;
			index=(index+I4(ActiveEdges,i,jp1,k,0))<<2;
			index=index+I4(ActiveEdges,i,j,k,1);
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
							p1=EdgePtPtr(Point,i,j,k,0);
							p2=EdgePtPtr(Point,i,j,k,1);
							break;
						case 20:
						case 16:
							cod=cod-16;
							p1=EdgePtPtr(Point,i,j,k,0);
							p2=EdgePtPtr(Point,ip1,j,k,1);
							break;
						case 8:
							cod=cod-8;
							p1=EdgePtPtr(Point,i,jp1,k,0);
							p2=EdgePtPtr(Point,ip1,j,k,1);
							break;
						case 4:
							cod=cod-4;
							p1=EdgePtPtr(Point,i,jp1,k,0);
							p2=EdgePtPtr(Point,i,j,k,1);
							break;
						case 2:
							cod=cod-2;
							p1=EdgePtPtr(Point,i,j,k,0);
							p2=EdgePtPtr(Point,i,jp1,k,0);
							break;
						case 1:
							cod=cod-1;
							p1=EdgePtPtr(Point,i,j,k,1);
							p2=EdgePtPtr(Point,ip1,j,k,1);
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
			index=I4(ActiveEdges,i,j,k,0)<<2;
			index=(index+I4(ActiveEdges,ip1,j,k,2))<<2;
			index=(index+I4(ActiveEdges,i,j,kp1,0))<<2;
			index=index+I4(ActiveEdges,i,j,k,2);
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
							p1=EdgePtPtr(Point,i,j,k,0);
							p2=EdgePtPtr(Point,i,j,k,2);
							break;
						case 20:
						case 16:
							cod=cod-16;
							p1=EdgePtPtr(Point,i,j,k,0);
							p2=EdgePtPtr(Point,ip1,j,k,2);
							break;
						case 8:
							cod=cod-8;
							p1=EdgePtPtr(Point,i,j,k+1,0);
							p2=EdgePtPtr(Point,ip1,j,k,2);
							break;
						case 4:
							cod=cod-4;
							p1=EdgePtPtr(Point,i,j,kp1,0);
							p2=EdgePtPtr(Point,i,j,k,2);
							break;
						case 2:
							cod=cod-2;
							p1=EdgePtPtr(Point,i,j,k,0);
							p2=EdgePtPtr(Point,i,j,kp1,0);
							break;
						case 1:
							cod=cod-1;
							p1=EdgePtPtr(Point,i,j,k,2);
							p2=EdgePtPtr(Point,ip1,j,k,2);
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
  float	ratio;
  ratio=(Level-*l1)/(*l2-*l1);
  pt[0]=v1[0]+(v2[0]-v1[0])*ratio;
  pt[1]=v1[1]+(v2[1]-v1[1])*ratio;
  pt[2]=v1[2]+(v2[2]-v1[2])*ratio;
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
		if((I3(VertexCodes,i,j,k))&&(!I3(VertexCodes,i+1,j,k)))
			{
#ifdef Trace
ECount++;
#endif
			I4(ActiveEdges,i,j,k,0)=2;
			IsosurfInterpolate(
				O4Ptr(Coord,i,j,k,0,CurOff),
				O3Ptr(Data,i,j,k,CurOff),
				O4Ptr(Coord,i+1,j,k,0,CurOff),
				O3Ptr(Data,i+1,j,k,CurOff),
				&(EdgePt(Point,i,j,k,0).Point[0]));
			}
		else if(!(I3(VertexCodes,i,j,k))&&(I3(VertexCodes,i+1,j,k)))
			{
#ifdef Trace
ECount++;
#endif
			I4(ActiveEdges,i,j,k,0)=1;
			IsosurfInterpolate(
				O4Ptr(Coord,i,j,k,0,CurOff),
				O3Ptr(Data,i,j,k,CurOff),
				O4Ptr(Coord,i+1,j,k,0,CurOff),
				O3Ptr(Data,i+1,j,k,CurOff),
				&(EdgePt(Point,i,j,k,0).Point[0]));
			}
		else
			I4(ActiveEdges,i,j,k,0)=0;
		}

	for(i=0;i<Max[0];i++)
	for(j=0;j<(Max[1]-1);j++)
	for(k=0;k<Max[2];k++)
		{
		if((I3(VertexCodes,i,j,k))&&(!I3(VertexCodes,i,j+1,k)))
			{
#ifdef Trace
ECount++;
#endif
			I4(ActiveEdges,i,j,k,1)=2;
			IsosurfInterpolate(
				O4Ptr(Coord,i,j,k,0,CurOff),
				O3Ptr(Data,i,j,k,CurOff),
				O4Ptr(Coord,i,j+1,k,0,CurOff),
				O3Ptr(Data,i,j+1,k,CurOff),
				&(EdgePt(Point,i,j,k,1).Point[0]));
			}
		else if(!(I3(VertexCodes,i,j,k))&&(I3(VertexCodes,i,j+1,k)))
			{
#ifdef Trace
ECount++;
#endif
			I4(ActiveEdges,i,j,k,1)=1;
			IsosurfInterpolate(
				O4Ptr(Coord,i,j,k,0,CurOff),
				O3Ptr(Data,i,j,k,CurOff),
				O4Ptr(Coord,i,j+1,k,0,CurOff),
				O3Ptr(Data,i,j+1,k,CurOff),
				&(EdgePt(Point,i,j,k,1).Point[0]));
			}
		else
			I4(ActiveEdges,i,j,k,1)=0;
		}

	for(i=0;i<Max[0];i++)
	for(j=0;j<Max[1];j++)
	for(k=0;k<(Max[2]-1);k++)
		{
		if((I3(VertexCodes,i,j,k))&&(!I3(VertexCodes,i,j,k+1)))
			{
#ifdef Trace
ECount++;
#endif
			I4(ActiveEdges,i,j,k,2)=2;
			IsosurfInterpolate(
				O4Ptr(Coord,i,j,k,0,CurOff),
				O3Ptr(Data,i,j,k,CurOff),
				O4Ptr(Coord,i,j,k+1,0,CurOff),
				O3Ptr(Data,i,j,k+1,CurOff),
				&(EdgePt(Point,i,j,k,2).Point[0]));
			}
		else if(!(I3(VertexCodes,i,j,k))&&(I3(VertexCodes,i,j,k+1)))
			{
#ifdef Trace
ECount++;
#endif
			I4(ActiveEdges,i,j,k,2)=1;
			IsosurfInterpolate(
				O4Ptr(Coord,i,j,k,0,CurOff),
				O3Ptr(Data,i,j,k,CurOff),
				O4Ptr(Coord,i,j,k+1,0,CurOff),
				O3Ptr(Data,i,j,k+1,CurOff),
				&(EdgePt(Point,i,j,k,2).Point[0]));
			}
		else
			I4(ActiveEdges,i,j,k,2)=0;
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
		if((O3(Data,i,j,k,CurOff)>Level))
			{
			I3(VertexCodes,i,j,k)=1;
			VCount++;
			}
		else
			I3(VertexCodes,i,j,k)=0;
		}
#ifdef Trace
printf(" IsosurfCodeVertices: %i of %i vertices above level\n",VCount,
CurDim[0]*CurDim[1]*CurDim[2]);
#endif
	return(VCount);
}
