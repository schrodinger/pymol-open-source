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
-*   Larry Coopet (various optimizations)
-* 
-*
Z* ------------------------------------------------------------------- */

#include"os_predef.h"
#include"os_std.h"


#include"MemoryDebug.h"
#include"Err.h"
#include"OOMac.h"
#include"Map.h"
#include"Setting.h"
#include"Feedback.h"
#include"MemoryCache.h"
#include"Base.h"

#ifndef true
#define true 1
#endif

#ifndef false
#define false 0
#endif

static MapType *_MapNew(PyMOLGlobals *G,float range,float *vert,int nVert,
                        float *extent,int *flag,int group_id,int block_id);

float MapGetDiv(MapType *I)
{
  return I->Div;
}

void MapFree(MapType *I)
{
  if(I) {
    CacheFreeP(I->G,I->Head,I->group_id,I->block_base + cCache_map_head_offset,false);
    CacheFreeP(I->G,I->Link,I->group_id,I->block_base + cCache_map_link_offset,false);
    CacheFreeP(I->G,I->EHead,I->group_id,I->block_base + cCache_map_ehead_offset,false);
    CacheFreeP(I->G,I->EMask,I->group_id,I->block_base + cCache_map_emask_offset,false);
    VLACacheFreeP(I->G,I->EList,I->group_id,I->block_base + cCache_map_elist_offset,false);
  }
  OOFreeP(I);
}

void MapCacheInit(MapCache *M,MapType *I,int group_id,int block_base)
{
  PyMOLGlobals *G=I->G;
  M->G = G;
  M->block_base = I->block_base;
  M->Cache = CacheCalloc(G,int,I->NVert,group_id,block_base + cCache_map_cache_offset);
  M->CacheLink = CacheAlloc(G,int,I->NVert,group_id,block_base + cCache_map_cache_link_offset);
  M->CacheStart = -1;
  /*  p=M->Cache;
  for(a=0;a<I->NVert;a++)
  *(p++) = 0;*/
}

void MapCacheReset(MapCache *M)
{
	register int    i		= M->CacheStart;
	register int* 	cachep	= M->Cache;
	register int* 	clinkp	= M->CacheLink;
	register int    i1 = 0, i2 = 0, i3 = 0, i4 = 0, ii;
	while(i >= 0)  /* believe it or not, unrolling gives us almost 10%!!! */
	{
		ii = clinkp[i];
		i1 = i;
		i = ii;
		if(ii >= 0) {
			ii = clinkp[ii];
			i2 = i;
			i = ii;
		}
		cachep[i1] = 0; /* this doesn't look safe, but it is. i1-i4 are always valid indices*/
		if(ii >= 0) {
			ii = clinkp[ii];
			i3 = i;
			i = ii;
		}	
		cachep[i2] = 0; /* this doesn't look safe, but it is. i1-i4 are always valid indices*/
		if(ii >= 0) {
			ii = clinkp[ii];
			i4 = i;
			i = ii;
		}      
		cachep[i3] = 0; /* this doesn't look safe, but it is. i1-i4 are always valid indices*/
		cachep[i4] = 0; /* this doesn't look safe, but it is. i1-i4 are always valid indices*/
	}
	M->CacheStart = -1;
}

void MapCacheFree(MapCache *M,int group_id,int block_base)
{
  CacheFreeP(M->G,M->Cache,group_id,block_base + cCache_map_cache_offset,false);
  CacheFreeP(M->G,M->CacheLink,group_id,block_base + cCache_map_cache_link_offset,false);
}

#define MapSafety 0.01F

int MapInside(MapType *I,float *v,int *a,int *b,int *c) /* special version for ray-tracing */
{
  register int		atmp, btmp, ctmp;
  register float	iDiv	= I->recipDiv; 

  atmp	= (int)((v[0] - I->Min[0]) * iDiv) + MapBorder;
  btmp	= (int)((v[1] - I->Min[1]) * iDiv) + MapBorder;
  ctmp	= (int)((v[2] - I->Min[2]) * iDiv) + MapBorder;

  if(atmp < I->iMin[0])
    {
		if((I->iMin[0] - atmp) > 3)
        return(-1);
		else 
        atmp = I->iMin[0]; 
    }
  else if(atmp > I->iMax[0]) 
    { 
		if((atmp - I->iMax[0]) > 3) 
        return(-1);
		else 
        atmp = I->iMax[0]; 
    }


	if(btmp < I->iMin[1]) 
	{ 
		if((I->iMin[1] - btmp) > 3) 
			return(-1);
		else 
			btmp = I->iMin[1];
	}
	else if(btmp > I->iMax[1]) 
	{
		if((btmp - I->iMax[1]) > 3)
			return(-1);
		else 
			btmp = I->iMax[1];
	}

   /*    printf("%d %d %d\n",ctmp,I->iMin[2],I->iMax[2]);*/
	if(ctmp < I->iMin[2]) 
	{ 
		if((I->iMin[2] - ctmp) > 3) 
			return(-1);
		else 
			ctmp = I->iMin[2];
	}
	else if(ctmp > I->iMax[2]) 
	{
		if((ctmp - I->iMax[2]) > 3)
        return(0); /* just keep searching... */
		else 
			ctmp = I->iMax[2];
	}

   /*   printf("%d %d %d %d %d %d %d %d %d\n",
          atmp,I->iMin[0],I->iMax[0],
          btmp,I->iMin[1],I->iMax[1],
          ctmp,I->iMin[2],I->iMax[2]);
   */

	if(!*(MapEStart(I,atmp,btmp,ctmp)))
     return(0);
		
	*a		= atmp;
	*b		= btmp;
	*c		= ctmp;
	
	return(1);  
}

int MapInsideXY(MapType *I,float *v,int *a,int *b,int *c) /* special version for ray-tracing */
{
	int		atmp, btmp, ctmp;
	const float	iDiv	= I->recipDiv; 
		
	atmp	= (int)((v[0] - I->Min[0]) * iDiv) + MapBorder;
	btmp	= (int)((v[1] - I->Min[1]) * iDiv) + MapBorder;
	ctmp	= (int)((v[2] - I->Min[2]) * iDiv) + MapBorder + 1;

	if(atmp < I->iMin[0])
	{
		if((I->iMin[0] - atmp) > 1)
			return(false);
		else 
			atmp = I->iMin[0]; 
	}
	else if(atmp > I->iMax[0]) 
	{ 
		if((atmp - I->iMax[0]) > 1) 
			return(false);
		else 
			atmp = I->iMax[0]; 
	}

	if(btmp < I->iMin[1]) 
	{ 
		if((I->iMin[1] - btmp) > 1) 
			return(false);
		else 
			btmp = I->iMin[1];
	}
	else if(btmp > I->iMax[1]) 
	{
		if((btmp - I->iMax[1]) > 1)
			return(false);
		else 
			btmp = I->iMax[1];
	}

	if(!*(I->EMask + I->Dim[1]*atmp + btmp))
		return(false);
		
	if(ctmp < I->iMin[2])
		ctmp = I->iMin[2];
	else if(ctmp > I->iMax[2])
		ctmp = I->iMax[2];

	*a		= atmp;
	*b		= btmp;
	*c		= ctmp;
	
	return(true);  
}

#define ELIST_GROW_FACTOR 3

void MapSetupExpressXY(MapType *I,int n_vert, int negative_start) /* setup a list of XY neighbors for each square */
{
  PyMOLGlobals *G=I->G;
  int n, a,b,c,flag;
  register int d,e,i;
  unsigned int mapSize;
  int st, dim2;
  int n_alloc = n_vert * 15; /* emprical est. */
  
  PRINTFD(G,FB_Map)
	" MapSetupExpressXY-Debug: entered.\n"
	ENDFD;
  
  mapSize = I->Dim[0]*I->Dim[1]*I->Dim[2];
  I->EHead = CacheCalloc(G,int,mapSize,I->group_id,I->block_base + cCache_map_ehead_offset);
  ErrChkPtr(G,I->EHead);
  I->EList = VLACacheMalloc(G,n_alloc,sizeof(int),ELIST_GROW_FACTOR,0,
                            I->group_id,I->block_base + cCache_map_elist_offset); 
  I->EMask = CacheCalloc(G,int,I->Dim[0]*I->Dim[1],
                         I->group_id,I->block_base + cCache_map_emask_offset);
  
  n		= 1;
  dim2	= I->Dim[2];
  
  for(a = I->iMin[0]; a <= I->iMax[0]; a++) {
    for(b = I->iMin[1]; b <= I->iMax[1]; b++) {
      for(c = I->iMin[2]; c <= I->iMax[2]; c++) { /* a better alternative exists... */
        int	*iPtr1	= (I->Head + ((a-1) * I->D1D2) + ((b-1)*dim2) + c);
        
        st		= n;
        flag	= false;
        
        for(d = a-1; d <= a+1; d++) {
          /*int	*iPtr2	= (I->Head + (d * I->D1D2) + ((b-1)*dim2) + c);*/
          int	*iPtr2	= iPtr1;
          
          for(e = b-1; e <= b+1; e++) {
            /*i	= *MapFirst(I,d,e,c);*/
            i	= *iPtr2;
            if(i >= 0) {
              flag	= true;
              while(i >= 0) {
                
                VLACacheCheck(G,I->EList,int,n,I->group_id,I->block_base + cCache_map_elist_offset);
                I->EList[n]	= i;
                n++;
                i	= MapNext(I,i);
              }
            }
            
            iPtr2	+= dim2;
          }
          
          iPtr1 += I->D1D2;
        }
        
        if(flag) {
          *(I->EMask + I->Dim[1]*a + b) = true;
          *(MapEStart(I,a,b,c))= negative_start ? -st : st;
          VLACacheCheck(G,I->EList,int,n,I->group_id,I->block_base + cCache_map_elist_offset);
          I->EList[n]=-1;
          n++;
        }
      }
    }
  }
  
  PRINTFB(G,FB_Map,FB_Blather)
    " MapSetupExpressXY: %d rows in express table\n",n
    ENDFB(G);
  
  I->NEElem=n;
  VLACacheSize(G,I->EList,int,I->NEElem,I->group_id,I->block_base + cCache_map_elist_offset);
  PRINTFD(G,FB_Map)
    " MapSetupExpressXY-Debug: leaving...\n"
    ENDFD;
}

void MapSetupExpressXYVert(MapType *I,float *vert,int n_vert,int negative_start) /* setup a list of XY neighbors for each square */
{
  PyMOLGlobals *G=I->G;
  int		h, n, a,b,c;
  int		j,k,dim2;
  float	*v;
  int		*eBase, *hBase;
  int      n_alloc = n_vert * 15; /* emprical est. */
	
  PRINTFD(G,FB_Map)
    " MapSetupExpressXYVert-Debug: entered n_vert = %d negative_start = %d\n",n_vert, negative_start
	ENDFD;
	
  /*mapSize 	= I->Dim[0]*I->Dim[1]*I->Dim[2];*/
  I->EHead	= CacheCalloc(G,int, I->Dim[0]*I->Dim[1]*I->Dim[2],
                          I->group_id,I->block_base + cCache_map_ehead_offset);
  I->EMask    = CacheCalloc(G,int,I->Dim[0]*I->Dim[1],
                            I->group_id,I->block_base + cCache_map_emask_offset);
  ErrChkPtr(G,I->EHead);
  I->EList	= VLACacheMalloc(G,n_alloc,sizeof(int),ELIST_GROW_FACTOR,0,
                             I->group_id,I->block_base + cCache_map_elist_offset); /* autozero */
	
  n		= 1;
  v		= vert;
  dim2	= I->Dim[2];
		
  for(h = 0; h < n_vert; h++) {
    MapLocus(I,v,&j,&k,&c);
      
    eBase	= I->EHead + ((j-1) * I->D1D2) + ((k-1)*dim2) + c;
    hBase	= I->Head + (((j-1)-1) * I->D1D2);
      
    for(a = j-1; a <= j+1; a++) {
      int *ePtr1	= eBase;
        
      for(b = k-1; b <= k+1; b++) {
        if( *ePtr1 == 0 ) { /* if we haven't yet expanded this voxel... */
          int	st		= n;
          int	flag	= false;
          int *hPtr1 = hBase + dim2 * (b-1) + (c-1);

          register int d,e,f;

          for(d = a-1; d <= a+1; d++) {
            int *hPtr2 = hPtr1;
            for(e = b-1; e <= b+1; e++) {
              int *hPtr3 = hPtr2;
              for(f = c-1; f <= c+1; f++) {
                /*                register int i	= *MapFirst(I,d,e,f);*/
                register int i	= *hPtr3;

                if(i > -1) {
                  flag = true;
                  while(i > -1) {
                    VLACacheCheck(G,I->EList,int,n,I->group_id,I->block_base + cCache_map_elist_offset);
                    I->EList[n]	= i;
                    n++;
                    i = MapNext(I,i);
                  }
                }
                hPtr3 ++;
              }
              hPtr2 += dim2;
            }
            hPtr1 += I->D1D2;
          }
					
          if(flag) {
            *(I->EMask + I->Dim[1]*a + b) = true;
            *(MapEStart(I,a,b,c))=  negative_start ? -st : st;
            VLACacheCheck(G,I->EList,int,n,I->group_id,I->block_base + cCache_map_elist_offset);
            I->EList[n] = -1;
            n++;
          }
        }
				
        ePtr1	+= dim2;
      }
			
      eBase	+= I->D1D2;
      hBase	+= I->D1D2;
    }
		
    v += 3;
  }

  PRINTFB(G,FB_Map,FB_Blather)
    " MapSetupExpressXYVert: %d rows in express table\n",n
    ENDFB(G);

	I->NEElem = n;
    VLACacheSize(G,I->EList,int,I->NEElem,I->group_id,I->block_base + cCache_map_elist_offset);

	PRINTFD(G,FB_Map)
		" MapSetupExpressXYVert-Debug: leaving...\n"
	ENDFD;
}

void MapSetupExpressPerp(MapType *I, float *vert, float front,int nVertHint,int negative_start)
{
  PyMOLGlobals *G=I->G;
  int n=0;
  int a,b,c,i;

  unsigned int mapSize;
  int st;
  int      n_alloc = nVertHint * 15; /* emprical est. */

  register int iMin0 = I->iMin[0];
  register int iMin1 = I->iMin[1];
  register int iMax0 = I->iMax[0];
  register int iMax1 = I->iMax[1];
  register float iDiv	= I->recipDiv;
  register float min0 = I->Min[0] * iDiv;
  register float min1 = I->Min[1] * iDiv;
  register float base0, base1;
  register float perp_factor,premult,*v0;
  register int *emask, dim1, *link, *ptr1,*ptr2;
  
  PRINTFD(G,FB_Map)
    " MapSetupExpress-Debug: entered.\n"
    ENDFD;

  mapSize = I->Dim[0]*I->Dim[1]*I->Dim[2];
  I->EHead=CacheCalloc(G,int,mapSize,
                 I->group_id,I->block_base + cCache_map_ehead_offset);
  ErrChkPtr(G,I->EHead);
  I->EList=VLACacheMalloc(G,n_alloc,sizeof(int),ELIST_GROW_FACTOR,0,
                          I->group_id,I->block_base + cCache_map_elist_offset);
  
  I->EMask = CacheCalloc(G,int,I->Dim[0]*I->Dim[1],
                            I->group_id,I->block_base + cCache_map_emask_offset);
  
  emask = I->EMask;
  dim1 = I->Dim[1];
  link = I->Link;
  premult = -front*iDiv;

  n=1;
  for(a=(iMin0-1);a<=(iMax0+1);a++)
    for(b=(iMin1-1);b<=(iMax1+1);b++)
      for(c=(I->iMin[2]-1);c<=(I->iMax[2]+1);c++) {

        register int d,e,f;
        
        /* compute a "shadow" mask for all vertices */
        
        i=*MapFirst(I,a,b,c);
        while(i>=0) {
          v0 = vert + 3*i;
          perp_factor = premult/v0[2];
          base0 = v0[0] * perp_factor;
          base1 = v0[1] * perp_factor;
          
          d	= (int)(base0 - min0);
          e	= (int)(base1 - min1);
          
          d += MapBorder;
          e += MapBorder;
          
          if(d < iMin0) {
            d = iMin0; 
          } else if(d > iMax0) { 
            d = iMax0; 
          }
          
          if(e < iMin1) { 
            e = iMin1;
          } else if(e > iMax1) {
            e = iMax1;
          }
          i = link[i];
          ptr2 = (ptr1 = emask + dim1*(d-1) + (e-1));
          *(ptr2++) = true;
          *(ptr2++) = true;
          *(ptr2++) = true;
          ptr2 = (ptr1 += dim1);
          *(ptr2++) = true;
          *(ptr2++) = true;
          *(ptr2++) = true;
          ptr2 = (ptr1 += dim1);
          *(ptr2++) = true;
          *(ptr2++) = true;
          *(ptr2++) = true;
        }
        
        {
          const int am1=a-1, ap1=a+1, bm1=b-1, bp1=b+1, cm1=c-1, cp1=c+1;
          const int dim2 = I->Dim[2];
          register int flag=false;
          register int *hPtr1	= I->Head + ((am1) * I->D1D2) + ((bm1)*dim2) + cm1;
          st=n;
          for(d=am1;d<=ap1;d++) {
            register int *hPtr2 = hPtr1;
            for(e=bm1;e<=bp1;e++) {
              register int *hPtr3 = hPtr2;
              for(f=cm1;f<=cp1;f++) {
                i = *(hPtr3++);
                /*                i=*MapFirst(I,d,e,f);*/
                if(i>=0) {
                  flag=true;
                  while(i>=0) {
                    VLACacheCheck(G,I->EList,int,n,I->group_id,
                                  I->block_base + cCache_map_elist_offset);
                    I->EList[n]=i;
                    i=link[i];
                    n++;
                  }
                }
              }
              hPtr2 += dim2;
            }
            hPtr1 += I->D1D2;
          }
          if(flag) {
            *(MapEStart(I,a,b,c)) = negative_start ? -st : st;
            VLACacheCheck(G,I->EList,int,n,I->group_id,
                          I->block_base + cCache_map_elist_offset);
            I->EList[n]=-1;
            n++;
          }
        }
      }
  PRINTFB(G,FB_Map,FB_Blather)
    " MapSetupExpressPerp: %d rows in express table\n",n
    ENDFB(G);
  I->NEElem=n;
  VLACacheSize(G,I->EList,int,I->NEElem,I->group_id,I->block_base + cCache_map_elist_offset);

  PRINTFD(G,FB_Map)
    " MapSetupExpress-Debug: leaving...n=%d\n",n
    ENDFD;
  
}

void MapSetupExpress(MapType *I) /* setup a list of neighbors for each square */
{
  register PyMOLGlobals *G=I->G;
  register int n=0;
  register int c,d,e,f,i,cm1,cp2,D1D2=I->D1D2,D2=I->Dim[2];
  register int mx2=I->iMax[2];
  register int *link = I->Link;
  register int st,flag;
  register int *i_ptr3,*i_ptr4,*i_ptr5;
  register int *e_list;
#ifdef _MemoryCache_ON
  register int block_offset = I->block_base + cCache_map_elist_offset;
  register int group_id = I->group_id;
#endif
  int mx0=I->iMax[0],mx1=I->iMax[1],a,am1,ap2,*i_ptr1,b,bm1,bp2,*i_ptr2;
  unsigned int mapSize;

  PRINTFD(G,FB_Map)
    " MapSetupExpress-Debug: entered.\n"
    ENDFD;

  mapSize = I->Dim[0]*I->Dim[1]*I->Dim[2];
  I->EHead=CacheCalloc(G,int,mapSize,group_id,I->block_base + cCache_map_ehead_offset);
  ErrChkPtr(G,I->EHead);
  e_list=VLACacheMalloc(G,1000,sizeof(int),5,0,group_id,block_offset);

  n=1;
  for(a=(I->iMin[0]-1);a<=mx0;a++) {
    am1 = a-1;
    ap2 = a+2;
    i_ptr1 = I->Head + am1 * D1D2;
	 for(b=(I->iMin[1]-1);b<=mx1;b++) {
      bm1 = b-1;
      bp2 = b+2;	
      i_ptr2 = i_ptr1 + bm1 * D2;
      for(c=(I->iMin[2]-1);c<=mx2;c++) {
        st=n;
        cm1 = c-1;
        cp2 = c+2;
        flag=false;
        i_ptr5 = (i_ptr4 = (i_ptr3 = i_ptr2 + cm1));

        for(d=am1;d<ap2;d++) {
          for(e=bm1;e<bp2;e++) {
            for(f=cm1;f<cp2;f++) {
                /*i=*MapFirst(I,d,e,f);*/
              if((i= *(i_ptr5++))>=0) {
                flag=true;
                do {
                  VLACacheCheck(G,e_list,int,n,group_id,block_offset);
                  e_list[n++]=i;
                  /*i=MapNext(I,i); */
                } while( (i = link[i]) >= 0);
              }
            }
            i_ptr5 = (i_ptr4 += D2);
          }
          i_ptr5 = (i_ptr4 = (i_ptr3 += D1D2));
        }
        if(flag) {
          *(MapEStart(I,a,b,c))=st;
          VLACacheCheck(G,e_list,int,n,group_id, block_offset);
          e_list[n]=-1;
          n++;
        } else {
          *(MapEStart(I,a,b,c))=0;
        }
      }
    }
  }
  I->EList = e_list;
  I->NEElem=n;
  VLACacheSize(G,I->EList,int,I->NEElem,group_id,block_offset);
  PRINTFD(G,FB_Map)
    " MapSetupExpress-Debug: leaving...n=%d\n",n
    ENDFD;

}

void MapLocus(MapType *I,float *v,int *a,int *b,int *c)
{
	int		at, bt, ct;
	float	invDiv	= I->recipDiv; 
	
	at	= (int)((v[0] - I->Min[0]) * invDiv) + MapBorder;
	bt	= (int)((v[1] - I->Min[1]) * invDiv) + MapBorder;
	ct	= (int)((v[2] - I->Min[2]) * invDiv) + MapBorder;
	
	/* range checking...*/
	if(at < I->iMin[0])			at = I->iMin[0];
	else if(at > I->iMax[0])	at = I->iMax[0];
	
	if(bt < I->iMin[1])			bt = I->iMin[1];
	else if(bt > I->iMax[1])	bt = I->iMax[1];
	
	if(ct < I->iMin[2])			ct = I->iMin[2];
	else if(ct > I->iMax[2])	ct = I->iMax[2];
	
	*a	= at;
	*b	= bt;
	*c	= ct;
}

int *MapLocusEStart(MapType *I,float *v)
{
	register int a,b,c;
	float invDiv = I->recipDiv; 

	a=(int)(((v[0]-I->Min[0])*invDiv)+MapBorder);
	b=(int)(((v[1]-I->Min[1])*invDiv)+MapBorder);
	c=(int)(((v[2]-I->Min[2])*invDiv)+MapBorder);
	if(a<I->iMin[0]) a=I->iMin[0];
	else if(a>I->iMax[0]) a=I->iMax[0];
	if(b<I->iMin[1]) b=I->iMin[1];
	else if(b>I->iMax[1]) b=I->iMax[1];
	if(c<I->iMin[2]) c=I->iMin[2];
	else if(c>I->iMax[2]) c=I->iMax[2];
	return (I->EHead + ((a) * I->D1D2) + ((b)*I->Dim[2]) + (c));
}

int MapExclLocus(MapType *I,float *v,int *a,int *b,int *c)
{
  float invDiv = I->recipDiv; 

  *a=(int)(((v[0]-I->Min[0])*invDiv)+MapBorder);
  if(*a<I->iMin[0]) return(0);
  else if(*a>I->iMax[0]) return(0);
  *b=(int)(((v[1]-I->Min[1])*invDiv)+MapBorder);
  if(*b<I->iMin[1]) return(0);
  else if(*b>I->iMax[1]) return(0);
  *c=(int)(((v[2]-I->Min[2])*invDiv)+MapBorder);
  if(*c<I->iMin[2]) return(0);
  else if(*c>I->iMax[2]) return(0);
  return(1);
}

float MapGetSeparation(PyMOLGlobals *G,float range,float *mx,float *mn,float *diagonal)
{
  float maxSize;
  float size,maxSubDiv,divSize, subDiv[3];
  float maxCubed, subDivCubed;
  int a;
  maxSize = SettingGet(G,cSetting_hash_max);

  maxCubed = maxSize * maxSize * maxSize;

  /* find longest axis */

  subtract3f(mx,mn,diagonal);
  diagonal[0] = (float)fabs(diagonal[0]);
  diagonal[1] = (float)fabs(diagonal[1]);
  diagonal[2] = (float)fabs(diagonal[2]);
  
  size=diagonal[0];
  if(diagonal[1]>size) size=diagonal[1];
  if(diagonal[2]>size) size=diagonal[2];

  if(size==0.0) {
    diagonal[0]=1.0;
    diagonal[1]=1.0;
    diagonal[2]=1.0;
    size = 1.0;
  }

  /* compute maximum number of subdivisions */
  maxSubDiv = (float)(size/(range+MapSafety));
  if(maxSubDiv<1.0F)
    maxSubDiv = 1.0F;

  /* find resulting divSize */
  divSize = size/maxSubDiv;
  if(divSize<MapSafety)
    divSize = MapSafety;
  for(a=0;a<3;a++) {
    subDiv[a] = (float)((int)((diagonal[a] / divSize)+0.5F));
    subDiv[a] = (subDiv[a] < 1.0F) ? 1.0F : subDiv[a];
  }
  subDivCubed = subDiv[0]*subDiv[1]*subDiv[2];

  if(subDivCubed>maxCubed) {
    divSize = (float)(divSize / pow(maxCubed/subDivCubed,0.33333F));
  } else if(subDivCubed<maxCubed) {
    divSize = (float)(divSize * pow(subDivCubed/maxCubed,0.33333F));    
  }

  if(divSize<(range+MapSafety))
    divSize = range+MapSafety;

  if(Feedback(G,FB_Map,FB_Debugging)) {
    PRINTF
      " MapGetSeparation: range %8.3f divSize %8.3f size %8.3f\n",range,divSize,size
      ENDF(G);
    /*    dump3f(mx,"mx");
    dump3f(mn,"mn");
    dump3f(diagonal,"diagonal");
    printf("%8.3f\n",divSize);
  printf("divSize %8.3f\n",divSize);
    */
  }
  return(divSize);
}

MapType *MapNew(PyMOLGlobals *G,float range,float *vert,int nVert,float *extent)
{
  return(_MapNew(G,range,vert,nVert,extent,NULL,-1,0));
}

MapType *MapNewCached(PyMOLGlobals *G,float range,float *vert,int nVert,float *extent,int group_id,int block_id)
{
  return(_MapNew(G,range,vert,nVert,extent,NULL,group_id,block_id));
}

MapType *MapNewFlagged(PyMOLGlobals *G,float range,float *vert,int nVert,float *extent,int *flag)
{
  return(_MapNew(G,range,vert,nVert,extent,flag,-1,0));
}

static MapType *_MapNew(PyMOLGlobals *G,float range,float *vert,int nVert,float *extent,int *flag,
                        int group_id,int block_base)
{
  int a,c;
  int mapSize;
  int h,k,l;
  int *list;
  float *v,tmp_f;
  int firstFlag;
  Vector3f diagonal;

  OOAlloc(G,MapType);

  PRINTFD(G,FB_Map)
    " MapNew-Debug: entered.\n"
    ENDFD;
  I->G=G;
  I->group_id = group_id;
  I->block_base = block_base;
  I->Head = NULL;
  I->Link = NULL;
  I->EHead = NULL;
  I->EList = NULL;
  I->EMask = NULL;
  I->NEElem=0;
  
  I->Link=CacheAlloc(G,int,nVert,group_id,block_base + cCache_map_link_offset);
  ErrChkPtr(G,I->Link);

  for(a=0;a<nVert;a++)
	 I->Link[a] = -1;

  if(extent) {
    I->Min[0]=extent[0];
    I->Max[0]=extent[1];
    I->Min[1]=extent[2];
    I->Max[1]=extent[3];
    I->Min[2]=extent[4];
    I->Max[2]=extent[5];
  }  else	 {
    I->Min[0]=0.0F;
    I->Max[0]=0.0F;
    I->Min[1]=0.0F;
    I->Max[1]=0.0F;
    I->Min[2]=0.0F;
    I->Max[2]=0.0F;
    if(flag) {
      firstFlag=true;
      v=vert;
      for(a=0;a<nVert;a++) {
        if(flag[a]) {
          if(firstFlag) {
            for(c=0;c<3;c++) {
              I->Min[c]=v[c];
              I->Max[c]=v[c];
            }
            firstFlag=false;
          } else {
            for(c=0;c<3;c++) {
              
              if(I->Min[c]>v[c]) I->Min[c]=v[c];
              if(I->Max[c]<v[c]) I->Max[c]=v[c];
            }
          }
        }
        v+=3;
      }
    } else {
      if(nVert) {
        v=vert;
        for(c=0;c<3;c++) {
            I->Min[c] = v[c];
            I->Max[c] = v[c];
          }
        v+=3;
        for(a=1;a<nVert;a++)  {
          for(c=0;c<3;c++)  {
            if(I->Min[c]>v[c])
              I->Min[c]=v[c];
            if(I->Max[c]<v[c])
              I->Max[c]=v[c];
          }
          v+=3;
        }
      }
    }
  }
  
  /* sanity check */
  for(a=0;a<3;a++) {
    if(I->Min[a]>I->Max[a]) {
      tmp_f=I->Min[a];
      I->Max[a]=I->Min[a];
      I->Min[a]=tmp_f;
    }
  }
  
  if(Feedback(G,FB_Map,FB_Debugging)) {
    printf(" MapSetup: %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f\n",
           I->Min[0],I->Min[1],I->Min[2],
           I->Max[0],I->Max[1],I->Max[2]);
  }

  for(c=0;c<3;c++)
	 {
		I->Min[c]-=MapSafety;
		I->Max[c]+=MapSafety;
	 }

  if(range<0.0) { /* negative range is a flag to expand edges using "range".*/
	 range=-range;
	 for(c=0;c<3;c++)
		{
		  I->Min[c]-=range;
		  I->Max[c]+=range;
		}
  }

  /* compute final box size */
  I->Div = MapGetSeparation(G,range,I->Max,I->Min,diagonal);
  I->recipDiv = 1.0F/(I->Div); /* cache this */

  /* add borders to avoid special edge cases */
  I->Dim[0]=(int)((diagonal[0]/I->Div)+1+(2*MapBorder)); 
  I->Dim[1]=(int)((diagonal[1]/I->Div)+1+(2*MapBorder));
  I->Dim[2]=(int)((diagonal[2]/I->Div)+1+(2*MapBorder));

  if(Feedback(G,FB_Map,FB_Debugging)) {
    printf(" MapSetup: nVert: %d\n",nVert);
    printf(" MapSetup: I->Div: %8.3f\n",I->Div);
    printf(" MapSetup: %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f\n",
           I->Min[0],I->Min[1],I->Min[2],
           I->Max[0],I->Max[1],I->Max[2]);
    printf(" MapSetup: %8d %8d %8d\n",
           I->Dim[0],I->Dim[1],I->Dim[2]);
  }

  I->D1D2 = I->Dim[1]*I->Dim[2];

  I->iMin[0]=MapBorder;
  I->iMin[1]=MapBorder;
  I->iMin[2]=MapBorder;

  I->iMax[0]=I->Dim[0]-(1+MapBorder);
  I->iMax[1]=I->Dim[1]-(1+MapBorder);
  I->iMax[2]=I->Dim[2]-(1+MapBorder);

  /* compute size and allocate */
  mapSize = I->Dim[0]*I->Dim[1]*I->Dim[2];
  I->Head=CacheAlloc(G,int,mapSize,group_id,block_base + cCache_map_head_offset);
  /*printf("%d\n",mapSize);*/
  ErrChkPtr(G,I->Head);

  /* initialize */
  /*  for(a=0;a<I->Dim[0];a++)
	 for(b=0;b<I->Dim[1];b++)
		for(c=0;c<I->Dim[2];c++)
		*(MapFirst(I,a,b,c))=-1;*/

	/* Trick for fast clearing to -1! */
	memset( I->Head, 0xFF, mapSize * sizeof(int) );

  I->NVert = nVert;

  PRINTFD(G,FB_Map)
    " MapNew-Debug: creating 3D hash...\n"
    ENDFD;

  /* create 3-D hash of the vertices*/
  if(flag) {
	 v=vert;
	 for(a=0;a<nVert;a++)	{
       if(flag[a]) 
         if(MapExclLocus(I,v,&h,&k,&l)) {
           list = MapFirst(I,h,k,l);
           I->Link[a] = *list; 
           *list = a; /*add to top of list*/
         }
       v+=3;
     }
  } else {
    v=vert;
    for(a=0;a<nVert;a++) {
      if(MapExclLocus(I,v,&h,&k,&l)) {
        list = MapFirst(I,h,k,l);
        /*        printf("LINK %d %d %d %d %5.2f %5.2f %5.2f\n",a,h,k,l,v[0],v[1],v[2]);*/
        I->Link[a] = *list; 
        *list = a; /*add to top of list*/
      }
      v+=3;
    }
  }

  PRINTFD(G,FB_Map)
    " MapNew-Debug: leaving...\n"
    ENDFD;


  return(I);
}

