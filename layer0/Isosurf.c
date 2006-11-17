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

#include"OVRandom.h"
#include"OVContext.h"

#include"Isosurf.h"
#include"MemoryDebug.h"
#include"Err.h"
#include"Crystal.h"
#include"Vector.h"
#include"Feedback.h"
#include"PConv.h"
#include"P.h"
#include"Util.h"

#define Trace_OFF

#ifndef true
#define true 1
#endif

#ifndef false
#define false 0
#endif

#define O3(field,P1,P2,P3,offs) Ffloat3(field,(P1)+offs[0],(P2)+offs[1],(P3)+offs[2])

#define O3Ptr(field,P1,P2,P3,offs) Ffloat3p(field,(P1)+offs[0],(P2)+offs[1],(P3)+offs[2])

#define O4(field,P1,P2,P3,P4,offs) Ffloat4(field,(P1)+offs[0],(P2)+offs[1],(P3)+offs[2],P4)

#define O4Ptr(field,P1,P2,P3,P4,offs) Ffloat4p(field,(P1)+offs[0],(P2)+offs[1],(P3)+offs[2],P4)

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

struct _CIsosurf {
  
  CField	*VertexCodes;
  CField	*ActiveEdges;
  CField   *Point;
  int	NLine;
  int   Skip;
  int	AbsDim[3],CurDim[3],CurOff[3];
  int	Max[3];
  CField *Coord,*Data;
  float	Level;
  int	Code[256];
  
  int      *Num;
  int      NSeg;
  float    *Line;
  
};


static int	IsosurfAlloc(PyMOLGlobals *G,CIsosurf *II);
static void	IsosurfPurge(CIsosurf *II);
static int	IsosurfCurrent(CIsosurf *II);
static int	IsosurfCodeVertices(CIsosurf *II);
static void	IsosurfInterpolate(CIsosurf *II,float *v1,float *l1,float *v2,float *l2,float *pt);
static int	IsosurfFindActiveEdges(CIsosurf *II);
static int	IsosurfFindLines(CIsosurf *II);
static int	IsosurfDrawLines(CIsosurf *II);
static void	IsosurfCode(CIsosurf *II,char *bits1,char *bits2);
static int	IsosurfDrawPoints(CIsosurf *II);
static int	IsosurfPoints(CIsosurf *II);
static int IsosurfGradients(PyMOLGlobals *G,CSetting *set1, CSetting *set2,
                            CIsosurf *II,Isofield *field,
                            int *range,float min_level, float max_level);

#define IsosurfSubSize		64

static void  _IsosurfFree(CIsosurf *I)
{
  FreeP(I);
}

void  IsosurfFree(PyMOLGlobals *G)
{
  _IsosurfFree(G->Isosurf);
  G->Isosurf = NULL;
}
/*===========================================================================*/
PyObject *IsosurfAsPyList(Isofield *field)
{
#ifdef _PYMOL_NOPY
  return NULL;
#else

  PyObject *result=NULL;

  result = PyList_New(4);

  PyList_SetItem(result,0,PConvIntArrayToPyList(field->dimensions,3));
  PyList_SetItem(result,1,PyInt_FromLong(field->save_points));
  PyList_SetItem(result,2,FieldAsPyList(field->data));
  if(field->save_points) 
    PyList_SetItem(result,3,FieldAsPyList(field->points));
  else
    PyList_SetItem(result,3,PConvAutoNone(NULL));
  return(PConvAutoNone(result));
#endif
}
/*===========================================================================*/
#ifdef _PYMOL_INLINE 
__inline__
#endif
 static void	IsosurfInterpolate(CIsosurf *I,float *v1,float *l1,float *v2,float *l2,float *pt)
{
  float	ratio;
  ratio=(I->Level-*l1)/(*l2-*l1);
  pt[0]=v1[0]+(v2[0]-v1[0])*ratio;
  pt[1]=v1[1]+(v2[1]-v1[1])*ratio;
  pt[2]=v1[2]+(v2[2]-v1[2])*ratio;
}
/*===========================================================================*/
Isofield *IsosurfNewFromPyList(PyMOLGlobals *G,PyObject *list)
{
#ifdef _PYMOL_NOPY
  return NULL;
#else

  int ok=true;
  int dim4[4];
  int a;
  int ll;

  Isofield *result = NULL;
  if(ok) ok=(list!=NULL);
  if(ok) ok=PyList_Check(list);
  if(ok) ll = PyList_Size(list);
  /* TO ENABLE BACKWARDS COMPATIBILITY...
   Always check ll when adding new PyList_GetItem's */
  if(ok) ok=((result=mmalloc(sizeof(Isofield)))!=NULL);
  if(ok) {result->data=NULL;result->points=NULL;}
  if(ok) ok=PConvPyListToIntArrayInPlace(PyList_GetItem(list,0),result->dimensions,3);
  if(ok) ok = PConvPyIntToInt(PyList_GetItem(list,1),&result->save_points);
  if(ok) ok=((result->data = FieldNewFromPyList(G,PyList_GetItem(list,2)))!=NULL);
  if(ok) {
    if(result->save_points)
      ok=((result->points = FieldNewFromPyList(G,PyList_GetItem(list,3)))!=NULL);
    else {
      for(a=0;a<3;a++)
        dim4[a]=result->dimensions[a];
      dim4[3] = 3;
      ok=((result->points = FieldNew(G,dim4,4,sizeof(float),cFieldFloat))!=NULL);
    }
  }
  if(!ok) {
    if(result) {
      if(result->data)
        FieldFree(result->data);
      if(result->points)
        FieldFree(result->points);
      mfree(result);
      result=NULL;
    }
  }
  result->gradients = NULL;
  return(result);
#endif
}

Isofield *IsosurfNewCopy(PyMOLGlobals *G,Isofield *src)
{
  int ok=true;

  Isofield *result = Calloc(Isofield,1);

  copy3f(src->dimensions, result->dimensions);
  result->save_points = src->save_points;
  
  ok = ((result->data = FieldNewCopy(G,src->data)) != NULL);
  ok = ((result->points = FieldNewCopy(G,src->points)) != NULL);

  result->gradients = NULL;
  if(!ok) {
    if(result->data) FieldFree(result->data);
    if(result->points) FieldFree(result->points);
    FreeP(result);
  }
  return result;
}

/*===========================================================================*/
void IsofieldComputeGradients(PyMOLGlobals *G,Isofield *field)
{
  int dim[4];
  int a,b,c;
  CField *data = field->data;
  CField *gradients;

  if(!field->gradients) {

    /* compute gradients relative to grid axis spacing */

    for(a=0;a<3;a++)
      dim[a]=field->dimensions[a];
    dim[3]=3;
    field->gradients = FieldNew(G,dim,4,sizeof(float),cFieldFloat);
    gradients =  field->gradients;
    dim[3] = 3;

    /* bulk internal */

    for(a=1;a<(dim[0]-1);a++) {
      for(b=1;b<(dim[1]-1);b++) {
        for(c=1;c<(dim[2]-1);c++) {
          F4(gradients,a,b,c,0) = (F3(data,a+1,b,c) - F3(data,a-1,b,c))/2.0F;
          F4(gradients,a,b,c,1) = (F3(data,a,b+1,c) - F3(data,a,b-1,c))/2.0F;
          F4(gradients,a,b,c,2) = (F3(data,a,b,c+1) - F3(data,a,b,c-1))/2.0F;
        }
      }
    }
    
    for(a=0;a<dim[0];a+=(dim[0]-1)) {
      
      /* 'a' faces */
      for(b=1;b<(dim[1]-1);b++) {
        for(c=1;c<(dim[2]-1);c++) {
          if(!a) {
            F4(gradients,a,b,c,0) = (F3(data,a+1,b,c) - F3(data,a,b,c));
          } else {
            F4(gradients,a,b,c,0) = (F3(data,a,b,c) - F3(data,a-1,b,c));
          }
          F4(gradients,a,b,c,1) = (F3(data,a,b+1,c) - F3(data,a,b-1,c))/2.0F;
          F4(gradients,a,b,c,2) = (F3(data,a,b,c+1) - F3(data,a,b,c-1))/2.0F;
        }
      }
      
      /* 'c' edges and all eight corners */
      for(b=0;b<dim[1];b+=(dim[1]-1)) {
        for(c=0;c<dim[2];c++) {
          
          if(!a) {
            F4(gradients,a,b,c,0) = (F3(data,a+1,b,c) - F3(data,a,b,c));
          } else {
            F4(gradients,a,b,c,0) = (F3(data,a,b,c) - F3(data,a-1,b,c));
          }
          
          if(!b) {
            F4(gradients,a,b,c,1) = (F3(data,a,b+1,c) - F3(data,a,b,c));
          } else {
            F4(gradients,a,b,c,1) = (F3(data,a,b,c) - F3(data,a,b-1,c));
          }
          
          if(!c) {
            F4(gradients,a,b,c,2) = (F3(data,a,b,c+1) - F3(data,a,b,c));
          } else if(c<(dim[2]-1)) {
            F4(gradients,a,b,c,2) = (F3(data,a,b,c+1) - F3(data,a,b,c-1))/2.0F;
          } else {
            F4(gradients,a,b,c,2) = (F3(data,a,b,c) - F3(data,a,b,c-1));            
          }
        }
      }
      
      /* 'b' edges  */
      for(c=0;c<dim[2];c+=(dim[2]-1)) {
        for(b=1;b<(dim[1]-1);b++) {
          if(!a) {
            F4(gradients,a,b,c,0) = (F3(data,a+1,b,c) - F3(data,a,b,c));
          } else {
            F4(gradients,a,b,c,0) = (F3(data,a,b,c) - F3(data,a-1,b,c));
          }
          
          F4(gradients,a,b,c,1) = (F3(data,a,b+1,c) - F3(data,a,b-1,c))/2.0F;          
          
          if(!c) {
            F4(gradients,a,b,c,2) = (F3(data,a,b,c+1) - F3(data,a,b,c));
          } else if(c<(dim[2]-1)) {
            F4(gradients,a,b,c,2) = (F3(data,a,b,c+1) - F3(data,a,b,c-1))/2.0F;
          } else {
            F4(gradients,a,b,c,2) = (F3(data,a,b,c) - F3(data,a,b,c-1));            
          }
        }
      }
    }
    
    for(b=0;b<dim[1];b+=(dim[1]-1)) {
      
      for(a=1;a<(dim[0]-1);a++) {
        
        /* 'b' faces */
        
        for(c=1;c<(dim[2]-1);c++) {
          F4(gradients,a,b,c,0) = (F3(data,a+1,b,c) - F3(data,a-1,b,c))/2.0F;
          if(!b) {
            F4(gradients,a,b,c,1) = (F3(data,a,b+1,c) - F3(data,a,b,c));
          } else {
            F4(gradients,a,b,c,1) = (F3(data,a,b,c) - F3(data,a,b-1,c));
          }
          F4(gradients,a,b,c,2) = (F3(data,a,b,c+1) - F3(data,a,b,c-1))/2.0F;
        }
        
        /* 'a' edges */
        
        for(c=0;c<dim[2];c+=(dim[2]-1)) {
          F4(gradients,a,b,c,0) = (F3(data,a+1,b,c) - F3(data,a-1,b,c))/2.0F;
          if(!b) {
            F4(gradients,a,b,c,1) = (F3(data,a,b+1,c) - F3(data,a,b,c));
          } else {
            F4(gradients,a,b,c,1) = (F3(data,a,b,c) - F3(data,a,b-1,c));
          }
          if(!c) {
            F4(gradients,a,b,c,2) = (F3(data,a,b,c+1) - F3(data,a,b,c));
          } else {
            F4(gradients,a,b,c,2) = (F3(data,a,b,c) - F3(data,a,b,c-1));
          }
        }
      }
    }
    
    
    /* 'c' faces */
    
    for(c=0;c<dim[2];c+=(dim[2]-1)) {
      for(a=1;a<(dim[0]-1);a++) {
        for(b=1;b<(dim[1]-1);b++) {
          F4(gradients,a,b,c,0) = (F3(data,a+1,b,c) - F3(data,a-1,b,c))/2.0F;
          F4(gradients,a,b,c,1) = (F3(data,a,b+1,c) - F3(data,a,b-1,c))/2.0F;
          if(!c) {
            F4(gradients,a,b,c,2) = (F3(data,a,b,c+1) - F3(data,a,b,c));
          } else {
            F4(gradients,a,b,c,2) = (F3(data,a,b,c) - F3(data,a,b,c-1));
          }
        }
      }
    }
  }
}

/*===========================================================================*/
Isofield *IsosurfFieldAlloc(PyMOLGlobals *G,int *dims)
{
  int dim4[4];
  int a;
  Isofield *result;

  for(a=0;a<3;a++)
    dim4[a]=dims[a];
  dim4[3] = 3;

  /* Warning: ...FromPyList also allocs and inits from the heap */

  result=mmalloc(sizeof(Isofield));
  ErrChkPtr(G,result);
  result->data = FieldNew(G,dims,3,sizeof(float),cFieldFloat);
  ErrChkPtr(G,result->data);
  result->points = FieldNew(G,dim4,4,sizeof(float),cFieldFloat);
  ErrChkPtr(G,result->points);
  result->dimensions[0]=dims[0];
  result->dimensions[1]=dims[1];
  result->dimensions[2]=dims[2];
  result->save_points=true;
  result->gradients = NULL;
  return(result);
}
/*===========================================================================*/
void IsosurfFieldFree(PyMOLGlobals *G,Isofield *field)
{
  if(field->gradients)
    FieldFree(field->gradients);
  FieldFree(field->points);
  FieldFree(field->data);
  mfree(field);
}
/*===========================================================================*/
static void	IsosurfCode(CIsosurf *II,char *bits1,char *bits2)
{
  register CIsosurf *I=II;
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

	I->Code[sum1]=sum2;
#ifdef Trace
	printf("IsosurfCode: %s (%i) -> %s (%i)\n",bits1,sum1,bits2, sum2);
#endif
}
/*===========================================================================*/
static CIsosurf *IsosurfNew(PyMOLGlobals *G)
{
	int	c;
   register CIsosurf *I = Calloc(CIsosurf,1);
              
	I->VertexCodes=NULL;
	I->ActiveEdges=NULL;
	I->Point=NULL;
	I->Line=NULL;
	I->Skip=0;
	for(c=0;c<256;c++)
	  I->Code[c]=-1;
		
/*___  
 | / |
 |/  |
 |___|
 32
*/
	IsosurfCode(I,"10000010","100000");			
	IsosurfCode(I,"01000001","100000");

/*___  
 | \ |
 |  \|
 |___|
 16
*/
	IsosurfCode(I,"10010000","010000");
	IsosurfCode(I,"01100000","010000");

/*___  
 |   |
 |  /|
 |_/_|
 8
*/
	IsosurfCode(I,"00101000","001000");
	IsosurfCode(I,"00010100","001000");

/*___  
 |   |
 |\  |
 |_\_|
 4
*/
	IsosurfCode(I,"00001001","000100");
	IsosurfCode(I,"00000110","000100");


/*___  
 | \ |
 |\ \|
 |_\_|
 16+4=20
*/

	IsosurfCode(I,"01101001","010100");

/*___  
 | / |
 |/ /|
 |_/_|
 32+8=40
*/
	IsosurfCode(I,"10010110","101000");


/*___  
 | | |
 | | |
 |_|_|
 2
*/
	IsosurfCode(I,"10001000","000010");
	IsosurfCode(I,"01000100","000010");

/*___  
 |   |
 |---|
 |___|
 1
*/
	IsosurfCode(I,"00100010","000001");
	IsosurfCode(I,"00010001","000001");
	
	return(I);
}
int IsosurfInit(PyMOLGlobals *G)
{
  G->Isosurf = IsosurfNew(G);
  return 1;
}
/*===========================================================================*/
void IsosurfGetRange(PyMOLGlobals *G,Isofield *field,
                     CCrystal *cryst,
                     float *mn,float *mx,int *range)
{
  float rmn[3],rmx[3];
  float imn[3],imx[3];
  float mix[24],imix[24];
  int a,b;
  PRINTFD(G,FB_Isosurface)
    " IsosurfGetRange: entered mn: %4.2f %4.2f %4.2f mx: %4.2f %4.2f %4.2f\n",
    mn[0],mn[1],mn[2],mx[0],mx[1],mx[2]
    ENDFD;

  for(a=0;a<3;a++) {
    rmn[a] = F4(field->points,0,0,0,a);
    rmx[a] = F4(field->points,
                field->dimensions[0]-1,
                field->dimensions[1]-1,
                field->dimensions[2]-1,a);
  }

  /* get min/max extents of map in fractional space */

  transform33f3f(cryst->RealToFrac,rmn,imn);
  transform33f3f(cryst->RealToFrac,rmx,imx);

  mix[ 0]=mn[0];
  mix[ 1]=mn[1];
  mix[ 2]=mn[2];

  mix[ 3]=mx[0];
  mix[ 4]=mn[1];
  mix[ 5]=mn[2];

  mix[ 6]=mn[0];
  mix[ 7]=mx[1];
  mix[ 8]=mn[2];

  mix[ 9]=mn[0];
  mix[10]=mn[1];
  mix[11]=mx[2];

  mix[12]=mx[0];
  mix[13]=mx[1];
  mix[14]=mn[2];

  mix[15]=mx[0];
  mix[16]=mn[1];
  mix[17]=mx[2];

  mix[18]=mn[0];
  mix[19]=mx[1];
  mix[20]=mx[2];

  mix[21]=mx[0];
  mix[22]=mx[1];
  mix[23]=mx[2];

  /* compute min/max of query in fractional space */

  for(b=0;b<8;b++) {
    transform33f3f(cryst->RealToFrac,mix+3*b,imix+3*b);
  }

  for(a=0;a<3;a++) {
    if(imx[a]!=imn[a]) { /* protect against div by zero */
      int b;
      int mini=0, maxi=0, tst_min, tst_max;
      float cur;
      for(b=0;b<8;b++) {
        cur = ((field->dimensions[a]-1)*(imix[a+3*b]-imn[a])/(imx[a]-imn[a]));
        tst_min = (int)floor(cur);
        tst_max = ((int)ceil(cur))+1;
        
        if(!b) {
          mini=tst_min;
          maxi=tst_max;
        } else {
          if(mini>tst_min)
            mini=tst_min;
          if(maxi<=tst_max)
            maxi=tst_max;
        }
      }

      range[a] = mini;

      range[a+3] = maxi;
    } else {
      range[a] = 0;
      range[a+3] = 1;
    }
    if(range[a]<0) range[a]=0;
    if(range[a]>field->dimensions[a]) range[a]=field->dimensions[a];
    if(range[a+3]<0) range[a+3]=0;
    if(range[a+3]>field->dimensions[a]) range[a+3]=field->dimensions[a];
  }
  PRINTFD(G,FB_Isosurface)
    " IsosurfGetRange: returning range: %d %d %d %d %d %d\n",
    range[0],range[1],range[2],range[3],range[4],range[5]
    ENDFD;
}
/*===========================================================================*/
int	IsosurfVolume(PyMOLGlobals *G,CSetting *set1,CSetting *set2,
                  Isofield *field,float level,int **num,
                  float **vert,int *range,int mode,int skip,float alt_level)
{
  register CIsosurf *I;
  if(PIsGlutThread()) {
    I = G->Isosurf;
  } else {
    I = IsosurfNew(G);
  }
  {
	int	ok=true;
	int	Steps[3];
	int	c,i,j,k;
	int	x,y,z;
    int range_store[6];
	I->Num = *num;
	I->Line = *vert;
    I->Skip = skip;
    if(range) {
      for(c=0;c<3;c++) {
        I->AbsDim[c]=field->dimensions[c];
        I->CurDim[c]=IsosurfSubSize+1;
        Steps[c]=((range[3+c]-range[c])-2)/IsosurfSubSize+1;
      }     
    } else {
      range=range_store;
      for(c=0;c<3;c++) {
        range[c]=0;
        range[3+c]=field->dimensions[c];
        I->AbsDim[c]=field->dimensions[c];
        I->CurDim[c]=IsosurfSubSize+1;
        Steps[c]=(I->AbsDim[c]-2)/IsosurfSubSize+1;
      }
    }
    
	I->Coord=field->points;
	I->Data=field->data;
	I->Level=level;
	if(ok) ok=IsosurfAlloc(G,I);
    
	I->NLine=0;
	I->NSeg=0;
	VLACheck(I->Num,int,I->NSeg);
	I->Num[I->NSeg]=I->NLine;

	if(ok) {
      switch(mode) {
      case 3:
        ok=IsosurfGradients(G,set1,set2,I,field,range,level,alt_level);
        break;
      default:
        for(i=0;i<Steps[0];i++)
          for(j=0;j<Steps[1];j++)
            for(k=0;k<Steps[2];k++) {
              I->CurOff[0]=IsosurfSubSize*i;
              I->CurOff[1]=IsosurfSubSize*j;
              I->CurOff[2]=IsosurfSubSize*k;
              for(c=0;c<3;c++)
                I->CurOff[c]+=range[c];
              for(c=0;c<3;c++)	{
                I->Max[c]=range[3+c]-I->CurOff[c];
                if(I->Max[c]>(IsosurfSubSize+1))
                  I->Max[c]=(IsosurfSubSize+1);
              }
              if(!(i||j||k))	{
                for(x=0;x<I->Max[0];x++)
                  for(y=0;y<I->Max[1];y++)
                    for(z=0;z<I->Max[2];z++)
                      for(c=0;c<3;c++)
                        EdgePt(I->Point,x,y,z,c).NLink=0;
              }
              
#ifdef Trace
              for(c=0;c<3;c++)
                printf(" IsosurfVolume: c: %i CurOff[c]: %i Max[c] %i\n",c,I->CurOff[c],I->Max[c]); 
#endif
              
              if(ok) 
                switch(mode) { 
                case 0: /* standard mode - want lines */
                  ok=IsosurfCurrent(I);
                  break;
                case 1: /* point mode - just want points on the isosurface */
                  ok=IsosurfPoints(I);
                  break;
                case 2:
                  /* reserved */
                break;
                }
            }
        IsosurfPurge(I);
      }
    }
   
    if(mode) {
      PRINTFB(G,FB_Isomesh,FB_Blather)
        " IsosurfVolume: Surface generated using %d dots.\n",I->NLine
        ENDFB(G);
    } else {
      PRINTFB(G,FB_Isomesh,FB_Blather)
        " IsosurfVolume: Surface generated using %d lines.\n",I->NLine
        ENDFB(G);
    }
    
    /* shrinks sizes for more efficient RAM usage */
    
    VLASize(I->Line,float,I->NLine*3);
    VLASize(I->Num,int,I->NSeg+1);
    
    I->Num[I->NSeg]=0;  /* important - must terminate the segment list */
    
	*vert = I->Line;
	*num = I->Num;  
    if(!PIsGlutThread()) {
      _IsosurfFree(I);
    }
	return(ok);
  }
}
/*===========================================================================*/
static int	IsosurfAlloc(PyMOLGlobals *G,CIsosurf *II)
{
  register CIsosurf *I = II;

	int	ok=true;
	int dim4[4];
   int a;
   for(a=0;a<3;a++)
     dim4[a]=I->CurDim[a];
   dim4[3]=3;

	I->VertexCodes=FieldNew(G,I->CurDim,3,sizeof(int),cFieldInt);
	ErrChkPtr(G,I->VertexCodes);
	I->ActiveEdges=FieldNew(G,dim4,4,sizeof(int),cFieldInt);
	ErrChkPtr(G,I->ActiveEdges);
	I->Point=FieldNew(G,dim4,4,sizeof(PointType),cFieldOther);
	ErrChkPtr(G,I->Point);
	if(!(I->VertexCodes&&I->ActiveEdges&&I->Point))
		{
		IsosurfPurge(I);
		ok=false;
		}
#ifdef Trace
	printf(" IsosurfAlloc: ok: %i\n",ok);
#endif
	return(ok);
}
/*===========================================================================*/
static void	IsosurfPurge(CIsosurf *II)
{
  register CIsosurf *I = II;
  if(I->VertexCodes) {
    FieldFree(I->VertexCodes);
    I->VertexCodes=NULL;
  }
  if(I->ActiveEdges) {
    FieldFree(I->ActiveEdges);
    I->ActiveEdges=NULL;
  }
  if(I->Point) {
    FieldFree(I->Point);
    I->Point=NULL;
  }
}
/*===========================================================================*/
static int	IsosurfCurrent(CIsosurf *II)
{
  register CIsosurf *I = II;
	int	ok=true;
	if(IsosurfCodeVertices(I)) {
      if(ok) ok=IsosurfFindActiveEdges(I);
      if(ok) ok=IsosurfFindLines(I);
      if(ok) ok=IsosurfDrawLines(I);
    }
	return(ok);
}
/*===========================================================================*/
static int	IsosurfPoints(CIsosurf *II)
{
  register CIsosurf *I = II;
  int	ok=true;
  if(IsosurfCodeVertices(I)) {
    if(ok) ok=IsosurfFindActiveEdges(I);
    if(ok) ok=IsosurfDrawPoints(I);
  }
  return(ok);
}
/*===========================================================================*/
static int IsosurfVertexInOrder(int *list,int a,int b)
{
  return(list[4*a+3]<=list[4*b+3]);
}

static int IsosurfGradients(PyMOLGlobals *G,CSetting *set1,CSetting *set2,
                            CIsosurf *II,Isofield *field,
                            int *range, float min_level, float max_level)
{
  register CIsosurf *I = II;
  int	ok=true;
  int spacing = SettingGet_f(G,set1,set2,cSetting_gradient_spacing);
  float step_size = SettingGet_f(G,set1,set2,cSetting_gradient_step_size);
  float max_walk = SettingGet_f(G,set1,set2,cSetting_gradient_max_length);
  float min_walk = SettingGet_f(G,set1,set2,cSetting_gradient_min_length);
  float min_slope = SettingGet_f(G,set1,set2,cSetting_gradient_min_slope);
  float min_dot = SettingGet_f(G,set1,set2,cSetting_gradient_min_dot);
  if(step_size<R_SMALL8)
    step_size = 0.1F;

  if(!field->gradients)
    IsofieldComputeGradients(G,field);
  if(field->gradients && I->AbsDim[0] && I->AbsDim[1] && I->AbsDim[2]) {
    int vol_size = I->AbsDim[0]*I->AbsDim[1]*I->AbsDim[2];
    int range_size = (range[3]-range[0])*(range[4]-range[1])*(range[5]-range[2]);
    int stride[3];
    int *order = Calloc(int,4*range_size);
    int *flag = Calloc(int,vol_size);
    int *active_cell = VLAlloc(int,1000);
    float cell_unit;
    stride[0] = 1;
    stride[1] = I->AbsDim[0];
    stride[2] = I->AbsDim[0]*I->AbsDim[1];

    {
      float *pos[4];
      pos[0] = Ffloat4p(I->Coord,0,0,0,0);
      pos[1] = Ffloat4p(I->Coord,1,0,0,0);
      pos[2] = Ffloat4p(I->Coord,0,1,0,0);
      pos[3] = Ffloat4p(I->Coord,0,0,1,0);
      cell_unit = (diff3f(pos[0],pos[1]) + 
                   diff3f(pos[0],pos[2]) + 
                   diff3f(pos[0],pos[3])) / 3.0F;

      max_walk /= cell_unit;
      min_walk /= cell_unit;
      step_size /= cell_unit;
      min_slope *= cell_unit;
    }

    {
      OVRandom *my_rand = OVRandom_NewBySeed(G->Context->heap,vol_size);
      int i,j,k,*p = order;
      for(k=range[2];k<range[5];k++) {
        for(j=range[1];j<range[4];j++) {
          for(i=range[0];i<range[3];i++) {
            p[0] = i;
            p[1] = j;
            p[2] = k;
            p[3] = OVRandom_Get_int31(my_rand);
            p+=4;
          }
        }
      }
      OVRandom_Del(my_rand);
    }
    UtilSortInPlace(G,order,range_size,4*sizeof(int),(UtilOrderFn*)IsosurfVertexInOrder);
    {
      int a;
      int *p = order;
      for(a=0;a<range_size;a++) {
        int abort_n_line = I->NLine;
        int abort_n_seg = I->NSeg;
        int pass;
        int n_active_cell = 0;
        float walk = max_walk;
        for(pass=0;pass<2;pass++) {
          int have_last = false;
          float last_gradient[3];
          int locus[3],*last_locus = NULL;
          float fract[3] = { 0.0F, 0.0F, 0.0F };
          int done = false;
          int n_vert = 0;
        
          locus[0] = p[0];
          locus[1] = p[1];
          locus[2] = p[2];
          
          while(!done) {
            {
              /* check array bounds violation */
              int b;
              for(b=0;b<3;b++) {
                while(fract[b]<0.0F) {
                  fract[b]+=1.0F;
                  locus[b]--;
                }
                while(fract[b]>=1.0F) { /* push everything into the lower box so the we can use it as a flag */
                  fract[b]-=1.0F;
                  locus[b]++;
                }
                while(locus[b]>=(I->AbsDim[b]-1)) {
                  if(fract[b]<=0.0F) {
                    locus[b]--;
                    fract[b]+=1.0F;
                  } else {
                    done = true;
                    break;
                  }
                }
                while(locus[b]<0) {
                  if(fract[b]>1.0F) {
                    locus[b]++;
                    fract[b]-=1.0F;
                  } else {
                    done = true;
                    break;
                  }
                }
              }
            }
            
            if(!done) {

              if((!have_last) || (have_last && ((locus[0]!=last_locus[0]) ||
                                                (locus[1]!=last_locus[1]) ||
                                                (locus[2]!=last_locus[2])))) {
                
                /* have we moved */

                /* stop if we leave target range */
                
                if((locus[0]<range[0])||(locus[1]<range[1])||(locus[2]<range[2])) 
                  done = true;

                if((locus[0]>(range[3]-2))||(locus[1]>(range[4]-2))||(locus[2]>(range[5]-2)))
                  done = true;

                /* stop if we hit an occupied square */
                
                if(!done) {
                  int *flag1 = flag + (locus[0] * stride[0]) + (locus[1] * stride[1]) + (locus[2] * stride[2]);
                  if( *flag1 ) {
                    done = true; 
                  }
                }
              }
            }

            if(!done) {
              /* stop if level exceeds desired ranges */
              float level = FieldInterpolatef(I->Data,locus[0],locus[1],locus[2],fract[0],fract[1],fract[2]);
              if((level<min_level) ||
                 (level>max_level))
                done = true;
            }

            if(!done) {
              float interp[3];
             
              /* interpolate gradient relative to grid */
              FieldInterpolate3f(field->gradients, locus, fract, interp);

              if(length3f(interp)<min_slope) { 
                /* if region is too flat, then bail */
                done = true;
              }
              
              if(!done) {
                
                /* add a line vertex at this point */
                
                float *f;
                VLACheck(I->Line,float,I->NLine*3+2);
                f=I->Line+(I->NLine*3);
                FieldInterpolate3f(field->points, locus, fract, f);
                I->NLine++;
                n_vert++;
                /* record locus for subsequent oblation */

                if((!have_last) || (have_last && ((locus[0]!=last_locus[0]) ||
                                                  (locus[1]!=last_locus[1]) ||
                                                  (locus[2]!=last_locus[2])))) {
                  
                  
                  VLACheck(active_cell, int, n_active_cell * 3 + 2);
                  {
                    int *xrd = active_cell + (n_active_cell * 3);
                    xrd[0] = locus[0];
                    xrd[1] = locus[1];
                    xrd[2] = locus[2];
                    n_active_cell++;
                    last_locus = xrd; /* warning: volatile */
                  }
                }

                /* adjust length of gradient vector */
                
                normalize3f(interp);                

                /* make sure gradient isn't too divergent to take another step */

                if(have_last) {
                  if(dot_product3f(interp,last_gradient)<min_dot)
                    done = true;
                }

                if(!done) {
                  /* take another step */

                  copy3f(interp,last_gradient);

                  if(pass) { 
                    scale3f(interp,-step_size,interp);
                  } else {
                    scale3f(interp,step_size,interp);
                  }
                  
                  walk -= step_size;
                  
                  if(walk<0.0F) 
                    done = true;
                  else {
                    /* move */
                    add3f(interp,fract,fract);
                  }
                  have_last = true;
                }

              }
            }
          }
          if(n_vert<2) {
            if(n_vert) {
              I->NLine = I->Num[I->NSeg]; /* quash isolated vertex */            
            }
          } else if(I->NLine!=I->Num[I->NSeg]) { /* any new lines? */
            VLACheck(I->Num,int,I->NSeg+1);
            I->Num[I->NSeg]=I->NLine-I->Num[I->NSeg];
            I->NSeg++;
            VLACheck(I->Num,int,I->NSeg);
            I->Num[I->NSeg]=I->NLine;
          }
        }
        if((max_walk-walk)<min_walk) { /* ignore short lines */
          I->NSeg = abort_n_seg;
          I->NLine = abort_n_line;
          I->Num[I->NSeg]=I->NLine;
        } else {
          /* keeping line, so oblate neighborhood */

          int *ac = active_cell;
          int b;
          register int cut_sq = spacing * spacing;
          for(b=0;b<n_active_cell;b++) {
            register int ii = ac[0], jj=ac[1], kk=ac[2];
            int i0 = ii - spacing;
            int j0 = jj - spacing;
            int k0 = kk - spacing;

            register int i1 = ii + spacing + 1;
            register int j1 = jj + spacing + 1;
            register int k1 = kk + spacing + 1;

            if(i0<0) i0 = 0;
            if(i1>=I->AbsDim[0]) i1 = I->AbsDim[0]-1;
            if(j0<0) j0 = 0;
            if(j1>=I->AbsDim[1]) j1 = I->AbsDim[1]-1;
            if(k0<0) k0 = 0;
            if(k1>=I->AbsDim[2]) k1 = I->AbsDim[2]-1;

            {
              register int i,j,k,a_plus_one = a+1;
              int *flag1 = flag + (i0 * stride[0]) + (j0 * stride[1]) + (k0 * stride[2]);
              for(k=k0;k<=k1;k++) {
                int *flag2 = flag1;
                int kk_sq = (kk-k);
                kk_sq = kk_sq * kk_sq;

                for(j=j0;j<=j1;j++) {
                  register int *flag3 = flag2;
                  register int jj_sq = (jj-j);
                  jj_sq = (jj_sq * jj_sq) + kk_sq;
                  
                  if( !(jj_sq>cut_sq)) {
                    for(i=i0;i<i1;i++) { 
                      if(!*flag3) {
                        register int tot_sq = (ii-i);
                        tot_sq = (tot_sq * tot_sq) + jj_sq;
                        if( !(tot_sq>cut_sq) ) { /* spherical avoidance */
                          *flag3 = a_plus_one;
                        }
                      }
                      flag3++;
                    } /* i */
                  }
                  flag2 += stride[1];
                } /* j */
                flag1 += stride[2];
              } /* k */
            }
            ac+=3;
          }
        }
        p+=4;
      }
    }
      VLAFreeP(active_cell);
    FreeP(order);
    FreeP(flag);
  }
  return(ok);
}
/*===========================================================================*/
static int	IsosurfDrawPoints(CIsosurf *II)
{
  register CIsosurf *I = II;
  float *a,*b;
  int	i,j,k;
  int	ok=true;
  
  for(i=0;i<(I->Max[0]-1);i++)
    for(j=0;j<I->Max[1];j++)
      for(k=0;k<I->Max[2];k++) {		
        if((I3(I->VertexCodes,i,j,k))&&(!I3(I->VertexCodes,i+1,j,k))) {
          IsosurfInterpolate(I,
                             O4Ptr(I->Coord,i,j,k,0,I->CurOff),
                             O3Ptr(I->Data,i,j,k,I->CurOff),
                             O4Ptr(I->Coord,i+1,j,k,0,I->CurOff),
                             O3Ptr(I->Data,i+1,j,k,I->CurOff),
                             &(EdgePt(I->Point,i,j,k,0).Point[0]));
          
          VLACheck(I->Line,float,I->NLine*3+2);
          a=I->Line+(I->NLine*3);
          b=&(EdgePt(I->Point,i,j,k,0).Point[0]);
          *(a++)=*(b++);
          *(a++)=*(b++);
          *a=*b;
          I->NLine++;
        }  else if(!(I3(I->VertexCodes,i,j,k))&&(I3(I->VertexCodes,i+1,j,k))) {
          IsosurfInterpolate(I,
                             O4Ptr(I->Coord,i,j,k,0,I->CurOff),
                             O3Ptr(I->Data,i,j,k,I->CurOff),
                             O4Ptr(I->Coord,i+1,j,k,0,I->CurOff),
                             O3Ptr(I->Data,i+1,j,k,I->CurOff),
                             &(EdgePt(I->Point,i,j,k,0).Point[0]));
          
          VLACheck(I->Line,float,I->NLine*3+2);
          a=I->Line+(I->NLine*3);
          b=&(EdgePt(I->Point,i,j,k,0).Point[0]);
          *(a++)=*(b++);
          *(a++)=*(b++);
          *a=*b;
          I->NLine++;
        } else
          I4(I->ActiveEdges,i,j,k,0)=0;
      }
  
  for(i=0;i<I->Max[0];i++)
    for(j=0;j<(I->Max[1]-1);j++)
      for(k=0;k<I->Max[2];k++)  {
        if((I3(I->VertexCodes,i,j,k))&&(!I3(I->VertexCodes,i,j+1,k))) {
          I4(I->ActiveEdges,i,j,k,1)=2;
          IsosurfInterpolate(I,
                             O4Ptr(I->Coord,i,j,k,0,I->CurOff),
                             O3Ptr(I->Data,i,j,k,I->CurOff),
                             O4Ptr(I->Coord,i,j+1,k,0,I->CurOff),
                             O3Ptr(I->Data,i,j+1,k,I->CurOff),
                             &(EdgePt(I->Point,i,j,k,1).Point[0]));
              
          VLACheck(I->Line,float,I->NLine*3+2);
          a=I->Line+(I->NLine*3);
          b=&(EdgePt(I->Point,i,j,k,1).Point[0]);
          *(a++)=*(b++);
          *(a++)=*(b++);
          *a=*b;
          I->NLine++;
              
        } else if(!(I3(I->VertexCodes,i,j,k))&&(I3(I->VertexCodes,i,j+1,k))) {
          IsosurfInterpolate(I,
                             O4Ptr(I->Coord,i,j,k,0,I->CurOff),
                             O3Ptr(I->Data,i,j,k,I->CurOff),
                             O4Ptr(I->Coord,i,j+1,k,0,I->CurOff),
                             O3Ptr(I->Data,i,j+1,k,I->CurOff),
                             &(EdgePt(I->Point,i,j,k,1).Point[0]));
              
          VLACheck(I->Line,float,I->NLine*3+2);
          a=I->Line+(I->NLine*3);
          b=&(EdgePt(I->Point,i,j,k,1).Point[0]);
          *(a++)=*(b++);
          *(a++)=*(b++);
          *a=*b;
          I->NLine++;
              
        }
      }
  
  for(i=0;i<I->Max[0];i++)
    for(j=0;j<I->Max[1];j++)
      for(k=0;k<(I->Max[2]-1);k++)        {
          if((I3(I->VertexCodes,i,j,k))&&(!I3(I->VertexCodes,i,j,k+1))) {
              IsosurfInterpolate(I,
                                 O4Ptr(I->Coord,i,j,k,0,I->CurOff),
                                 O3Ptr(I->Data,i,j,k,I->CurOff),
                                 O4Ptr(I->Coord,i,j,k+1,0,I->CurOff),
                                 O3Ptr(I->Data,i,j,k+1,I->CurOff),
                                 &(EdgePt(I->Point,i,j,k,2).Point[0]));
              
              VLACheck(I->Line,float,I->NLine*3+2);
              a=I->Line+(I->NLine*3);
              b=&(EdgePt(I->Point,i,j,k,2).Point[0]);
              *(a++)=*(b++);
              *(a++)=*(b++);
              *a=*b;
              I->NLine++;
              
          } else if(!(I3(I->VertexCodes,i,j,k))&&(I3(I->VertexCodes,i,j,k+1))) {
            IsosurfInterpolate(I,
                               O4Ptr(I->Coord,i,j,k,0,I->CurOff),
                               O3Ptr(I->Data,i,j,k,I->CurOff),
                               O4Ptr(I->Coord,i,j,k+1,0,I->CurOff),
                               O3Ptr(I->Data,i,j,k+1,I->CurOff),
                               &(EdgePt(I->Point,i,j,k,2).Point[0]));
              
            VLACheck(I->Line,float,I->NLine*3+2);
            a=I->Line+(I->NLine*3);
            b=&(EdgePt(I->Point,i,j,k,2).Point[0]);
            *(a++)=*(b++);
            *(a++)=*(b++);
            *a=*b;
            I->NLine++;
              
          }
      }

  if(I->NLine!=I->Num[I->NSeg]) { /* any new points? */
    VLACheck(I->Num,int,I->NSeg+1);
    I->Num[I->NSeg]=I->NLine-I->Num[I->NSeg];
    I->NSeg++;
    I->Num[I->NSeg]=I->NLine;
  }
  return(ok);
}
/*===========================================================================*/
static int	IsosurfDrawLines(CIsosurf *II)
{
  register CIsosurf *I = II;
	int	c,i,j,k;
	float	*a,*b;
	int	ok=true;
	PointType	*Cur,*Start,*Next;
	int	MaxLinks,MaxL,Cnt;
	int	NLink;
#ifdef Trace
	int	LCount=0;
#endif	
    
	for(i=0;i<I->Max[0];i++)
      for(j=0;j<I->Max[1];j++)
        for(k=0;k<I->Max[2];k++)
          for(c=0;c<3;c++) {
            Start=EdgePtPtr(I->Point,i,j,k,c);
            while(Start->NLink) 	{
              Cur=Start;
              VLACheck(I->Line,float,I->NLine*3+2);
              a=I->Line+(I->NLine*3);
              b=Cur->Point;
              *(a++)=*(b++);
              *(a++)=*(b++);
              *a=*b;
              I->NLine++;
			
              while(Cur)	{
				if(Cur->NLink)	{
                  Cur->NLink--;
                  NLink=Cur->NLink;
                  /* Choose point which has most links */
                  MaxL=NLink;
                  MaxLinks=Cur->Link[MaxL]->NLink;
                  Cnt=MaxL-1;
                  while(Cnt>=0)	{
                    if((Cur->Link[Cnt]->NLink)>MaxLinks)	{
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
                  while(Cnt>=0)	{
                    if(Next->Link[Cnt]==Cur)
                      break;
                    else
                      Cnt--;
                  }
                  if(Cnt>=0)	{
                    if(Cnt!=NLink)
                      Next->Link[Cnt]=Next->Link[NLink];
                  }
#ifdef Trace
                  else
                    printf(" error: IsosurfDrawLines:  can't find double link\n");
#endif
	
                  Cur=Next;
                  VLACheck(I->Line,float,I->NLine*3+2);
                  a=I->Line+(I->NLine*3);
                  b=Cur->Point;
                  *(a++)=*(b++);
                  *(a++)=*(b++);
                  *a=*b;
                  I->NLine++;
                }	else	{
#ifdef Trace
                  LCount++;
#endif
                  Cur=NULL;
                  if(I->NLine!=I->Num[I->NSeg]) { /* any new lines? */
                    VLACheck(I->Num,int,I->NSeg+1);
                    I->Num[I->NSeg]=I->NLine-I->Num[I->NSeg];
                    I->NSeg++;
                    VLACheck(I->Num,int,I->NSeg);
                    I->Num[I->NSeg]=I->NLine;
                  }
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
static int	IsosurfFindLines(CIsosurf *II)
{
  register CIsosurf *I = II;
	int	i,j,k,ip1,jp1,kp1;
	int	ok=true;
	int	index,cod;
	int	Max0m1,Max1m1,Max2m1;
    int skip = I->Skip;
#ifdef Trace
	int	LCount=0;
#endif
	PointType	*p1,*p2;
	
	Max0m1=I->Max[0]-1;
	Max1m1=I->Max[1]-1;
	Max2m1=I->Max[2]-1;
	for(i=0;i<I->Max[0];i++)
      for(j=0;j<I->Max[1];j++)
        for(k=0;k<I->Max[2];k++)	{
          ip1=i+1;
          jp1=j+1;
          kp1=k+1;
          if((j<Max1m1)&&(k<Max2m1)&&((!skip)||!(i%skip)))	{ /* i-plane */
            index=I4(I->ActiveEdges,i,j,k,1)<<2;
            index=(index+I4(I->ActiveEdges,i,jp1,k,2))<<2;
            index=(index+I4(I->ActiveEdges,i,j,kp1,1))<<2;
            index=index+I4(I->ActiveEdges,i,j,k,2);
            if(index)	{
              cod=I->Code[index];
#ifdef Trace
              if(index&&(cod<0))
                printf("IsosurfFindLines: bad index: %i \n",index);
#endif
              while(cod>0)	{
                  p1=NULL;
                  p2=NULL;
                  switch(cod)
                    {
                    case 40:
                    case 32:
                      cod=cod-32;
                      p1=EdgePtPtr(I->Point,i,j,k,1);
                      p2=EdgePtPtr(I->Point,i,j,k,2);
                      break;
                    case 20:
                    case 16:
                      cod=cod-16;
                      p1=EdgePtPtr(I->Point,i,j,k,1);
                      p2=EdgePtPtr(I->Point,i,jp1,k,2);
                      break;
                    case 8:
                      cod=cod-8;
                      p1=EdgePtPtr(I->Point,i,j,kp1,1);
                      p2=EdgePtPtr(I->Point,i,jp1,k,2);
                      break;
                    case 4:
                      cod=cod-4;
                      p1=EdgePtPtr(I->Point,i,j,kp1,1);
                      p2=EdgePtPtr(I->Point,i,j,k,2);
                      break;
                    case 2:
                      cod=cod-2;
                      p1=EdgePtPtr(I->Point,i,j,k,1);
                      p2=EdgePtPtr(I->Point,i,j,kp1,1);
                      break;
                    case 1:
                      cod=cod-1;
                      p1=EdgePtPtr(I->Point,i,j,k,2);
                      p2=EdgePtPtr(I->Point,i,jp1,k,2);
                      break;
                    default:
                      cod=0;
                      p1=NULL;
                      p2=NULL;
                      break;
                    }
                  if(p1&&p2)  {
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
            if((i<Max0m1)&&(j<Max1m1)&&((!skip)||!(k%skip)))	{ /* k-plane */
              index=I4(I->ActiveEdges,i,j,k,0)<<2;
              index=(index+I4(I->ActiveEdges,ip1,j,k,1))<<2;
              index=(index+I4(I->ActiveEdges,i,jp1,k,0))<<2;
              index=index+I4(I->ActiveEdges,i,j,k,1);
              if(index)	{
                cod=I->Code[index];
#ifdef Trace
                if(index&&(cod<0))
                  printf("IsosurfFindLines: bad index: %i \n",index);
#endif
                while(cod>0)	{
                  switch(cod)	{
                  case 40:
                  case 32:
                    cod=cod-32;
                    p1=EdgePtPtr(I->Point,i,j,k,0);
                    p2=EdgePtPtr(I->Point,i,j,k,1);
                    break;
                  case 20:
                  case 16:
                    cod=cod-16;
                    p1=EdgePtPtr(I->Point,i,j,k,0);
                    p2=EdgePtPtr(I->Point,ip1,j,k,1);
                    break;
                  case 8:
                    cod=cod-8;
                    p1=EdgePtPtr(I->Point,i,jp1,k,0);
                    p2=EdgePtPtr(I->Point,ip1,j,k,1);
                    break;
                  case 4:
                    cod=cod-4;
                    p1=EdgePtPtr(I->Point,i,jp1,k,0);
                    p2=EdgePtPtr(I->Point,i,j,k,1);
                    break;
                  case 2:
                    cod=cod-2;
                    p1=EdgePtPtr(I->Point,i,j,k,0);
                    p2=EdgePtPtr(I->Point,i,jp1,k,0);
                    break;
                  case 1:
                    cod=cod-1;
                    p1=EdgePtPtr(I->Point,i,j,k,1);
                    p2=EdgePtPtr(I->Point,ip1,j,k,1);
                    break;
                  default:
                    cod=0;
                    p1=NULL;
                    p2=NULL;
                    break;
                  }
                  if(p1&&p2)  {
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
            if((i<Max0m1)&&(k<Max2m1)&&((!skip)||!(j%skip)))	{ /* j-plane */
              index=I4(I->ActiveEdges,i,j,k,0)<<2;
              index=(index+I4(I->ActiveEdges,ip1,j,k,2))<<2;
              index=(index+I4(I->ActiveEdges,i,j,kp1,0))<<2;
              index=index+I4(I->ActiveEdges,i,j,k,2);
              if(index)      {
                cod=I->Code[index];
#ifdef Trace
                if(index&&(cod<0))
                  printf("IsosurfFindLines: bad index: %i \n",index);
#endif
                while(cod>0)	{
                  switch(cod)	{
                  case 40:
                  case 32:
                    cod=cod-32;
                    p1=EdgePtPtr(I->Point,i,j,k,0);
                    p2=EdgePtPtr(I->Point,i,j,k,2);
                    break;
                  case 20:
                  case 16:
                    cod=cod-16;
                    p1=EdgePtPtr(I->Point,i,j,k,0);
                    p2=EdgePtPtr(I->Point,ip1,j,k,2);
                    break;
                  case 8:
                    cod=cod-8;
                    p1=EdgePtPtr(I->Point,i,j,k+1,0);
                    p2=EdgePtPtr(I->Point,ip1,j,k,2);
                    break;
                  case 4:
                    cod=cod-4;
                    p1=EdgePtPtr(I->Point,i,j,kp1,0);
                    p2=EdgePtPtr(I->Point,i,j,k,2);
                    break;
                  case 2:
                    cod=cod-2;
                    p1=EdgePtPtr(I->Point,i,j,k,0);
                    p2=EdgePtPtr(I->Point,i,j,kp1,0);
                    break;
                  case 1:
                    cod=cod-1;
                    p1=EdgePtPtr(I->Point,i,j,k,2);
                    p2=EdgePtPtr(I->Point,ip1,j,k,2);
                    break;
                  default:
                    cod=0;
                    p1=NULL;
                    p2=NULL;
                    break;
                  }
                  if(p1&&p2) {
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
static int	IsosurfFindActiveEdges(CIsosurf *II)
{
  register CIsosurf *I = II;
	int	i,j,k;
	int	ok=true;
#ifdef Trace
	int	ECount=0;
#endif
	
	for(i=0;i<(I->Max[0]-1);i++)
	for(j=0;j<I->Max[1];j++)
	for(k=0;k<I->Max[2];k++)
		{		
		if((I3(I->VertexCodes,i,j,k))&&(!I3(I->VertexCodes,i+1,j,k)))
			{
#ifdef Trace
ECount++;
#endif
			I4(I->ActiveEdges,i,j,k,0)=2;
			IsosurfInterpolate(I,
				O4Ptr(I->Coord,i,j,k,0,I->CurOff),
				O3Ptr(I->Data,i,j,k,I->CurOff),
				O4Ptr(I->Coord,i+1,j,k,0,I->CurOff),
				O3Ptr(I->Data,i+1,j,k,I->CurOff),
				&(EdgePt(I->Point,i,j,k,0).Point[0]));
			}
		else if(!(I3(I->VertexCodes,i,j,k))&&(I3(I->VertexCodes,i+1,j,k)))
			{
#ifdef Trace
ECount++;
#endif
			I4(I->ActiveEdges,i,j,k,0)=1;
			IsosurfInterpolate(I,
				O4Ptr(I->Coord,i,j,k,0,I->CurOff),
				O3Ptr(I->Data,i,j,k,I->CurOff),
				O4Ptr(I->Coord,i+1,j,k,0,I->CurOff),
				O3Ptr(I->Data,i+1,j,k,I->CurOff),
				&(EdgePt(I->Point,i,j,k,0).Point[0]));
			}
		else
			I4(I->ActiveEdges,i,j,k,0)=0;
		}

	for(i=0;i<I->Max[0];i++)
	for(j=0;j<(I->Max[1]-1);j++)
	for(k=0;k<I->Max[2];k++)
		{
		if((I3(I->VertexCodes,i,j,k))&&(!I3(I->VertexCodes,i,j+1,k)))
			{
#ifdef Trace
ECount++;
#endif
			I4(I->ActiveEdges,i,j,k,1)=2;
			IsosurfInterpolate(I,
				O4Ptr(I->Coord,i,j,k,0,I->CurOff),
				O3Ptr(I->Data,i,j,k,I->CurOff),
				O4Ptr(I->Coord,i,j+1,k,0,I->CurOff),
				O3Ptr(I->Data,i,j+1,k,I->CurOff),
				&(EdgePt(I->Point,i,j,k,1).Point[0]));
			}
		else if(!(I3(I->VertexCodes,i,j,k))&&(I3(I->VertexCodes,i,j+1,k)))
			{
#ifdef Trace
ECount++;
#endif
			I4(I->ActiveEdges,i,j,k,1)=1;
			IsosurfInterpolate(I,
				O4Ptr(I->Coord,i,j,k,0,I->CurOff),
				O3Ptr(I->Data,i,j,k,I->CurOff),
				O4Ptr(I->Coord,i,j+1,k,0,I->CurOff),
				O3Ptr(I->Data,i,j+1,k,I->CurOff),
				&(EdgePt(I->Point,i,j,k,1).Point[0]));
			}
		else
			I4(I->ActiveEdges,i,j,k,1)=0;
		}

	for(i=0;i<I->Max[0];i++)
	for(j=0;j<I->Max[1];j++)
	for(k=0;k<(I->Max[2]-1);k++)
		{
		if((I3(I->VertexCodes,i,j,k))&&(!I3(I->VertexCodes,i,j,k+1)))
			{
#ifdef Trace
ECount++;
#endif
			I4(I->ActiveEdges,i,j,k,2)=2;
			IsosurfInterpolate(I,
				O4Ptr(I->Coord,i,j,k,0,I->CurOff),
				O3Ptr(I->Data,i,j,k,I->CurOff),
				O4Ptr(I->Coord,i,j,k+1,0,I->CurOff),
				O3Ptr(I->Data,i,j,k+1,I->CurOff),
				&(EdgePt(I->Point,i,j,k,2).Point[0]));
			}
		else if(!(I3(I->VertexCodes,i,j,k))&&(I3(I->VertexCodes,i,j,k+1)))
			{
#ifdef Trace
ECount++;
#endif
			I4(I->ActiveEdges,i,j,k,2)=1;
			IsosurfInterpolate(I,
				O4Ptr(I->Coord,i,j,k,0,I->CurOff),
				O3Ptr(I->Data,i,j,k,I->CurOff),
				O4Ptr(I->Coord,i,j,k+1,0,I->CurOff),
				O3Ptr(I->Data,i,j,k+1,I->CurOff),
				&(EdgePt(I->Point,i,j,k,2).Point[0]));
			}
		else
			I4(I->ActiveEdges,i,j,k,2)=0;
		}
#ifdef Trace
printf(" IsosurfFindActiveEdges: %i active edges found\n",ECount);
#endif
			
	return(ok);
}
/*===========================================================================*/
static int	IsosurfCodeVertices(CIsosurf *II)
{
  register CIsosurf *I = II;
	int	i,j,k;
	int	VCount=0;

	for(i=0;i<I->Max[0];i++)
      for(j=0;j<I->Max[1];j++)
        for(k=0;k<I->Max[2];k++) {
          if((O3(I->Data,i,j,k,I->CurOff)>I->Level)) {
            I3(I->VertexCodes,i,j,k)=1;
            VCount++;
          } else {
			I3(I->VertexCodes,i,j,k)=0;
          }
		}
#ifdef Trace
printf(" IsosurfCodeVertices: %i of %i vertices above level\n",VCount,
I->CurDim[0]*I->CurDim[1]*I->CurDim[2]);
#endif
	return(VCount);
}
