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
#include"os_gl.h"

#include"OOMac.h"
#include"ObjectMap.h"
#include"Base.h"
#include"MemoryDebug.h"
#include"Map.h"
#include"Parse.h"
#include"Isosurf.h"
#include"Vector.h"
#include"Color.h"
#include"main.h"
#include"Scene.h"
#include"PConv.h"
#include"Word.h"
#include"Vector.h"

#ifdef _PYMOL_NUMPY
typedef struct {
  PyObject_HEAD
  char *data;
  int nd;
  int *dimensions, *strides;
  PyObject *base;
  void *descr;
  int flags;
} MyArrayObject;
#endif

int ObjectMapNumPyArrayToMap(ObjectMap *I,PyObject *ary);

int ObjectMapSetBorder(ObjectMap *I,float level)
{
  int result = false;
  int a,b,c;

  c=I->FDim[2]-1;
  for(a=0;a<I->FDim[0];a++) 
    for(b=0;b<I->FDim[1];b++)
      {
        F3(I->Field->data,a,b,0) = level;
        F3(I->Field->data,a,b,c) = level;
      }

  a=I->FDim[0]-1;
  for(b=0;b<I->FDim[1];b++) 
    for(c=0;c<I->FDim[2];c++)
      {
        F3(I->Field->data,0,b,c) = level;
        F3(I->Field->data,a,b,c) = level;
      }

  b=I->FDim[1]-1;
  for(a=0;a<I->FDim[0];a++) 
    for(c=0;c<I->FDim[2];c++)
      {
        F3(I->Field->data,a,0,c) = level;
        F3(I->Field->data,a,b,c) = level;
      }
  return(result);
}

static void ObjectMapFree(ObjectMap *I);

static void ObjectMapFree(ObjectMap *I) {
  if(I->Field) {
    IsosurfFieldFree(I->Field);
    I->Field=NULL;
  }
  FreeP(I->Origin);
  FreeP(I->Dim);
  FreeP(I->Range);
  FreeP(I->Grid);
  OOFreeP(I->Crystal);
  ObjectPurge(&I->Obj);
  OOFreeP(I);
}

static void ObjectMapUpdate(ObjectMap *I) {
  SceneDirty();
}

static void ObjectMapRender(ObjectMap *I,int frame,CRay *ray,Pickable **pick,int pass)
{
  if(!pass) {
    if(I->Obj.RepVis[cRepExtent]) {
      if(ray) {
        float *vc;
        vc = ColorGet(I->Obj.Color);
        ray->fColor3fv(ray,vc);
        ray->fCylinder3fv(ray,I->Corner[0],I->Corner[1],0.20,vc,vc);
        ray->fCylinder3fv(ray,I->Corner[0],I->Corner[2],0.20,vc,vc);
        ray->fCylinder3fv(ray,I->Corner[2],I->Corner[3],0.20,vc,vc);
        ray->fCylinder3fv(ray,I->Corner[1],I->Corner[3],0.20,vc,vc);
        ray->fCylinder3fv(ray,I->Corner[0],I->Corner[4],0.20,vc,vc);
        ray->fCylinder3fv(ray,I->Corner[1],I->Corner[5],0.20,vc,vc);
        ray->fCylinder3fv(ray,I->Corner[2],I->Corner[6],0.20,vc,vc);
        ray->fCylinder3fv(ray,I->Corner[3],I->Corner[7],0.20,vc,vc);
        ray->fCylinder3fv(ray,I->Corner[4],I->Corner[5],0.20,vc,vc);
        ray->fCylinder3fv(ray,I->Corner[4],I->Corner[6],0.20,vc,vc);
        ray->fCylinder3fv(ray,I->Corner[6],I->Corner[7],0.20,vc,vc);
        ray->fCylinder3fv(ray,I->Corner[5],I->Corner[7],0.20,vc,vc);
      } else if(pick&&PMGUI) {
      } else if(PMGUI) {
        ObjectUseColor(&I->Obj);
        glDisable(GL_LIGHTING); 
        glBegin(GL_LINES);
        glVertex3fv(I->Corner[0]);
        glVertex3fv(I->Corner[1]);
        
        glVertex3fv(I->Corner[0]);
        glVertex3fv(I->Corner[2]);
        
        glVertex3fv(I->Corner[2]);
        glVertex3fv(I->Corner[3]);
        
        glVertex3fv(I->Corner[1]);
        glVertex3fv(I->Corner[3]);
        
        glVertex3fv(I->Corner[0]);
        glVertex3fv(I->Corner[4]);
        
        glVertex3fv(I->Corner[1]);
        glVertex3fv(I->Corner[5]);
        
        glVertex3fv(I->Corner[2]);
        glVertex3fv(I->Corner[6]);
        
        glVertex3fv(I->Corner[3]);
        glVertex3fv(I->Corner[7]);
        
        glVertex3fv(I->Corner[4]);
        glVertex3fv(I->Corner[5]);
        
        glVertex3fv(I->Corner[4]);
        glVertex3fv(I->Corner[6]);
        
        glVertex3fv(I->Corner[6]);
        glVertex3fv(I->Corner[7]);
        
        glVertex3fv(I->Corner[5]);
        glVertex3fv(I->Corner[7]);
        
        glEnd();
        glEnable(GL_LIGHTING);
      }
    }
  }
}

/*========================================================================*/
ObjectMap *ObjectMapNew(void)
{
OOAlloc(ObjectMap);

 ObjectInit((Object*)I);

 I->Obj.RepVis[cRepExtent]=true; 
 I->Crystal = CrystalNew();
 I->Field = NULL;
 I->Obj.type = cObjectMap;
 I->Obj.fFree = (void (*)(struct Object *))ObjectMapFree;
 I->Obj.fUpdate =  (void (*)(struct Object *)) ObjectMapUpdate;
 I->Obj.fRender =(void (*)(struct Object *, int, CRay *, Pickable **,int))ObjectMapRender;
 I->Origin = NULL;
 I->Dim = NULL;
 I->Range = NULL;
 I->Grid = NULL;
#ifdef _NOT_YET_NEEDED
  I->Obj.fGetNFrame = (int (*)(struct Object *)) ObjectMapGetNFrames;
#endif

  return(I);
}
/*========================================================================*/
int ObjectMapCCP4StrToMap(ObjectMap *I,char *CCP4Str,int bytes,int frame) {
  
  char *p;
  int *i;
  unsigned int *u;
  unsigned char *uc,c0,c1,c2,c3;
  float *f;
  float dens;
  int a,b,c,d,e;
  float v[3],vr[3],maxd,mind;
  int ok = true;
  int little_endian = 1,map_endian;
  /* CCP4 named from their docs */
  int nc,nr,ns;
  int map_mode;
  int ncstart,nrstart,nsstart;
  int nx,ny,nz;
  float xlen,ylen,zlen,alpha,beta,gamma;
  int n_skip;
  int mapc,mapr,maps;
  int cc[3],xref[3];
  int n_pts;
  double sum,sumsq;
  float mean,stdev;
  int normalize;

  normalize=(int)SettingGet(cSetting_normalize_ccp4_maps);
  maxd = FLT_MIN;
  mind = FLT_MAX;
  p=CCP4Str;
  little_endian = *((char*)&little_endian);
  map_endian = (*p||*(p+1));

  if(bytes<256*sizeof(int)) {
    PRINTFB(FB_ObjectMap,FB_Errors)
      " ObjectMapCCP4: Map appears to be truncated -- aborting."
      ENDFB;
    return(0);
  }
  if(little_endian!=map_endian) {
    PRINTFB(FB_ObjectMap,FB_Blather)
      " ObjectMapCCP4: Map appears to be reverse endian, swapping...\n"
      ENDFB;
    c = bytes;
    u = (unsigned int*)p;
    uc = (unsigned char *)u;
    while(c>3) {
      c0 = *(uc++);
      c1 = *(uc++);
      c2 = *(uc++);
      c3 = *(uc++);
      uc = (unsigned char *)u;
      *(uc++) = c3;
      *(uc++) = c2;
      *(uc++) = c1;
      *(uc++) = c0;
      u++;
      c-=4;
    }
  }
  i = (int*)p;
  nc=*(i++); /* columns */
  nr=*(i++); /* rows */
  ns=*(i++); /* sections */
  PRINTFB(FB_ObjectMap,FB_Blather)
    " ObjectMapCCP4: NC %d   NR %d   NS %d\n",
    nc,nr,ns
    ENDFB;
  map_mode = *(i++); /* mode */

  
  if(map_mode!=2) {
    PRINTFB(FB_ObjectMap,FB_Errors)
      "ObjectMapCCP4-ERR: Only map mode 2 currently supported (this map is mode %d)",map_mode
      ENDFB;
    return(0);
  }

  PRINTFB(FB_ObjectMap,FB_Blather)
    " ObjectMapCCP4: Map is mode %d.\n",map_mode
    ENDFB;

  ncstart = *(i++);
  nrstart = *(i++);
  nsstart = *(i++);

  PRINTFB(FB_ObjectMap,FB_Blather)
    " ObjectMapCCP4: NCSTART %d   NRSTART %d   NSSTART  %d\n",
    ncstart,nrstart,nsstart
    ENDFB;

  nx = *(i++);
  ny = *(i++);
  nz = *(i++);

  PRINTFB(FB_ObjectMap,FB_Blather)
    " ObjectMapCCP4: NX %d   NY %d   NZ  %d \n",
    nx,ny,nz
    ENDFB;

  xlen = *(float*)(i++);
  ylen = *(float*)(i++);
  zlen = *(float*)(i++);

  PRINTFB(FB_ObjectMap,FB_Blather)
    " ObjectMapCCP4: X %8.3f   Y %8.3f  Z  %8.3f \n",
    xlen,ylen,zlen
    ENDFB;

  alpha = *(float*)(i++);
  beta = *(float*)(i++);
  gamma = *(float*)(i++);

  PRINTFB(FB_ObjectMap,FB_Blather)
    " ObjectMapCCP4: alpha %8.3f   beta %8.3f  gamma %8.3f \n",
    alpha,beta,gamma
    ENDFB;

  mapc = *(i++);
  mapr = *(i++);
  maps = *(i++);

  PRINTFB(FB_ObjectMap,FB_Blather)
    " ObjectMapCCP4: MAPC %d   MAPR %d  MAPS  %d \n",
    mapc,mapr,maps
    ENDFB;

  i+=4;
  n_skip = *(i++);
  
  if(*(i++)) {
    PRINTFB(FB_ObjectMap,FB_Errors)
      "ObjectMapCCP4-ERR: PyMOL doesn't know how to handle skewed maps. Sorry!\n"
      ENDFB;
    return(0);
  }

  n_pts = nc*ns*nr;
  if(bytes<(n_skip + sizeof(int)*(256+n_pts))) {
    PRINTFB(FB_ObjectMap,FB_Errors)
      " ObjectMapCCP4: Map appears to be truncated -- aborting.\n"
      ENDFB;
    return(0);
  }

  if(n_pts>1) {
    f = (float*)(p+(sizeof(int)*256)+n_skip);
    c = n_pts;
    sum = 0.0;
    sumsq = 0.0;
    while(c--) {
      sumsq+=(*f)*(*f);
      sum+=*f++;
    }
    mean = sum/n_pts;
    stdev = sqrt1f((sumsq - (sum*sum/n_pts))/(n_pts-1));

    if(normalize) {
      PRINTFB(FB_ObjectMap,FB_Details)
        " ObjectMapCCP4: Normalizing with mean = %8.6f and stdev = %8.6f.\n",
        mean,stdev
        ENDFB;
    } else {
      PRINTFB(FB_ObjectMap,FB_Details)
        " ObjectMapCCP4: Map will not be normalized.\n ObjectMapCCP4: Current mean = %8.6f and stdev = %8.6f.\n",
        mean,stdev
        ENDFB;
    }

    if(stdev<0.000001)
      stdev = 1.0;
    
  } else {
    mean = 1.0;
    stdev = 1.0;
  }
  
  f = (float*)(p+(sizeof(int)*256)+n_skip);
  mapc--; /* convert to C indexing... */
  mapr--;
  maps--;

  xref[maps]=0;
  xref[mapr]=1;
  xref[maps]=2;

  xref[maps]=0;
  xref[mapr]=1;
  xref[mapc]=2;

  I->Div[0] = nx;
  I->Div[1] = ny;
  I->Div[2] = nz;

  I->FDim[mapc] = nc;
  I->Min[mapc] = ncstart;
  I->Max[mapc] = nc+ncstart;

  I->FDim[mapr] = nr;
  I->Min[mapr] = nrstart;
  I->Max[mapr] = nr+nrstart;

  I->FDim[maps] = ns;
  I->Min[maps] = nsstart;
  I->Max[maps] = ns+nsstart;
  
  I->Crystal->Dim[0] = xlen;
  I->Crystal->Dim[1] = ylen;
  I->Crystal->Dim[2] = zlen;

  I->Crystal->Angle[0] = alpha;
  I->Crystal->Angle[1] = beta;
  I->Crystal->Angle[2] = gamma;


  I->FDim[3]=3;
  if(!(I->FDim[0]&&I->FDim[1]&&I->FDim[2])) 
    ok=false;
  else {
    CrystalUpdate(I->Crystal);
    CrystalDump(I->Crystal);
    fflush(stdout);
    I->Field=IsosurfFieldAlloc(I->FDim);
    for(cc[maps]=0;cc[maps]<I->FDim[maps];cc[maps]++)
      {
        v[maps]=(cc[maps]+I->Min[maps])/((float)I->Div[maps]);

        for(cc[mapr]=0;cc[mapr]<I->FDim[mapr];cc[mapr]++) {
          v[mapr]=(cc[mapr]+I->Min[mapr])/((float)I->Div[mapr]);

          for(cc[mapc]=0;cc[mapc]<I->FDim[mapc];cc[mapc]++) {
            v[mapc]=(cc[mapc]+I->Min[mapc])/((float)I->Div[mapc]);

            if(normalize) 
              dens = (*f-mean)/stdev;
            else 
              dens = *f;
            F3(I->Field->data,cc[0],cc[1],cc[2]) = dens;
            if(maxd<*f) maxd = dens;
            if(mind>*f) mind = dens;
            f++;
            transform33f3f(I->Crystal->FracToReal,v,vr);
            for(e=0;e<3;e++) 
              F4(I->Field->points,cc[0],cc[1],cc[2],e) = vr[e];
          }
        }
      }
  }
  if(ok) {
    d = 0;
    for(c=0;c<I->FDim[2];c+=(I->FDim[2]-1))
      {
        v[2]=(c+I->Min[2])/((float)I->Div[2]);
        for(b=0;b<I->FDim[1];b+=(I->FDim[1]-1)) {
          v[1]=(b+I->Min[1])/((float)I->Div[1]);
          for(a=0;a<I->FDim[0];a+=(I->FDim[0]-1)) {
            v[0]=(a+I->Min[0])/((float)I->Div[0]);
            transform33f3f(I->Crystal->FracToReal,v,vr);
            copy3f(vr,I->Corner[d]);
            d++;
          }
        }
      }
  }
  
  if(ok) {
    v[2]=(I->Min[2])/((float)I->Div[2]);
    v[1]=(I->Min[1])/((float)I->Div[1]);
    v[0]=(I->Min[0])/((float)I->Div[0]);
    
    transform33f3f(I->Crystal->FracToReal,v,I->Obj.ExtentMin);
    
    v[2]=((I->FDim[2]-1)+I->Min[2])/((float)I->Div[2]);
    v[1]=((I->FDim[1]-1)+I->Min[1])/((float)I->Div[1]);
    v[0]=((I->FDim[0]-1)+I->Min[0])/((float)I->Div[0]);
    
    transform33f3f(I->Crystal->FracToReal,v,I->Obj.ExtentMax);
    I->Obj.ExtentFlag=true;
  }
#ifdef _UNDEFINED
  printf("%d %d %d %d %d %d %d %d %d\n",
         I->Div[0],
         I->Min[0],
         I->Max[0],
         I->Div[1],
         I->Min[1],
         I->Max[1],
         I->Div[2],
         I->Min[2],
         I->Max[2]);
  printf("Okay? %d\n",ok);
  fflush(stdout);
#endif
  if(!ok) {
    ErrMessage("ObjectMap","Error reading map");
  } else {
    printf(" ObjectMap: Map Read.  Range = %5.3f to %5.3f\n",mind,maxd);
  }
  return(ok);
}
/*========================================================================*/
int ObjectMapXPLORStrToMap(ObjectMap *I,char *XPLORStr,int frame) {
  
  char *p;
  int a,b,c,d,e;
  float v[3],vr[3],dens,maxd,mind;
  char cc[MAXLINELEN];
  int n;
  int ok = true;

  maxd = FLT_MIN;
  mind = FLT_MAX;
  p=XPLORStr;

  while(*p) {
    p = ParseNCopy(cc,p,8);
    if(!*cc) 
      p = ParseNextLine(p);
    else if(sscanf(cc,"%i",&n)==1) {
      p=ParseWordCopy(cc,p,MAXLINELEN);
      if(strstr(cc,"!NTITLE")||(!*cc)) {
        p=ParseNextLine(p);
        while(n--) {
          p=ParseNextLine(p);          
        } 
      } else if(strstr(cc,"REMARKS")) {
        p=ParseNextLine(p);          
      } else {
        break;
      }
    }
  }
  if(*p) { /* n contains first dimension */
    I->Div[0]=n;
    if(sscanf(cc,"%i",&I->Min[0])!=1) ok=false;
    p = ParseNCopy(cc,p,8); if(sscanf(cc,"%i",&I->Max[0])!=1) ok=false;
    p = ParseNCopy(cc,p,8); if(sscanf(cc,"%i",&I->Div[1])!=1) ok=false;
    p = ParseNCopy(cc,p,8); if(sscanf(cc,"%i",&I->Min[1])!=1) ok=false;    
    p = ParseNCopy(cc,p,8); if(sscanf(cc,"%i",&I->Max[1])!=1) ok=false;
    p = ParseNCopy(cc,p,8); if(sscanf(cc,"%i",&I->Div[2])!=1) ok=false;
    p = ParseNCopy(cc,p,8); if(sscanf(cc,"%i",&I->Min[2])!=1) ok=false;    
    p = ParseNCopy(cc,p,8); if(sscanf(cc,"%i",&I->Max[2])!=1) ok=false;
    p=ParseNextLine(p);
    p = ParseNCopy(cc,p,12); if(sscanf(cc,"%f",&I->Crystal->Dim[0])!=1) ok=false;
    p = ParseNCopy(cc,p,12); if(sscanf(cc,"%f",&I->Crystal->Dim[1])!=1) ok=false;
    p = ParseNCopy(cc,p,12); if(sscanf(cc,"%f",&I->Crystal->Dim[2])!=1) ok=false;
    p = ParseNCopy(cc,p,12); if(sscanf(cc,"%f",&I->Crystal->Angle[0])!=1) ok=false;
    p = ParseNCopy(cc,p,12); if(sscanf(cc,"%f",&I->Crystal->Angle[1])!=1) ok=false;
    p = ParseNCopy(cc,p,12); if(sscanf(cc,"%f",&I->Crystal->Angle[2])!=1) ok=false;
    p=ParseNextLine(p);
    p = ParseNCopy(cc,p,3);
    if(strcmp(cc,"ZYX")) ok=false;
    p=ParseNextLine(p);
    
  } else {
    ok=false;
  }
  if(ok) {
    I->FDim[0]=I->Max[0]-I->Min[0]+1;
    I->FDim[1]=I->Max[1]-I->Min[1]+1;
    I->FDim[2]=I->Max[2]-I->Min[2]+1;
    I->FDim[3]=3;
    if(!(I->FDim[0]&&I->FDim[1]&&I->FDim[2])) 
      ok=false;
    else {
      CrystalUpdate(I->Crystal);
      I->Field=IsosurfFieldAlloc(I->FDim);
      for(c=0;c<I->FDim[2];c++)
        {
          v[2]=(c+I->Min[2])/((float)I->Div[2]);
          p=ParseNextLine(p);
          for(b=0;b<I->FDim[1];b++) {
            v[1]=(b+I->Min[1])/((float)I->Div[1]);
            for(a=0;a<I->FDim[0];a++) {
              v[0]=(a+I->Min[0])/((float)I->Div[0]);
              p=ParseNCopy(cc,p,12);
              if(!cc[0]) {
                p=ParseNextLine(p);
                p=ParseNCopy(cc,p,12);                
              }
              if(sscanf(cc,"%f",&dens)!=1) {
                ok=false;
              } else {
                F3(I->Field->data,a,b,c) = dens;
                if(maxd<dens) maxd = dens;
                if(mind>dens) mind = dens;
              }
              transform33f3f(I->Crystal->FracToReal,v,vr);
              for(e=0;e<3;e++) 
                F4(I->Field->points,a,b,c,e) = vr[e];
            }
          }
          p=ParseNextLine(p);
        }
      if(ok) {
        d = 0;
        for(c=0;c<I->FDim[2];c+=(I->FDim[2]-1))
          {
            v[2]=(c+I->Min[2])/((float)I->Div[2]);
            for(b=0;b<I->FDim[1];b+=(I->FDim[1]-1)) {
              v[1]=(b+I->Min[1])/((float)I->Div[1]);
              for(a=0;a<I->FDim[0];a+=(I->FDim[0]-1)) {
                v[0]=(a+I->Min[0])/((float)I->Div[0]);
                transform33f3f(I->Crystal->FracToReal,v,vr);
                copy3f(vr,I->Corner[d]);
                d++;
              }
            }
          }
      }
    }
  }

  if(ok) {
    v[2]=(I->Min[2])/((float)I->Div[2]);
    v[1]=(I->Min[1])/((float)I->Div[1]);
    v[0]=(I->Min[0])/((float)I->Div[0]);

    transform33f3f(I->Crystal->FracToReal,v,I->Obj.ExtentMin);
    
    v[2]=((I->FDim[2]-1)+I->Min[2])/((float)I->Div[2]);
    v[1]=((I->FDim[1]-1)+I->Min[1])/((float)I->Div[1]);
    v[0]=((I->FDim[0]-1)+I->Min[0])/((float)I->Div[0]);

    transform33f3f(I->Crystal->FracToReal,v,I->Obj.ExtentMax);
    I->Obj.ExtentFlag=true;
  }
#ifdef _UNDEFINED
    printf("%d %d %d %d %d %d %d %d %d\n",
           I->Div[0],
           I->Min[0],
           I->Max[0],
           I->Div[1],
           I->Min[1],
           I->Max[1],
           I->Div[2],
           I->Min[2],
           I->Max[2]);
    printf("Okay? %d\n",ok);
    fflush(stdout);
#endif
  if(!ok) {
    ErrMessage("ObjectMap","Error reading map");
  } else {
    printf(" ObjectMap: Map Read.  Range = %5.3f to %5.3f\n",mind,maxd);
  }
    
  return(ok);
}
/*========================================================================*/
ObjectMap *ObjectMapReadXPLORStr(ObjectMap *I,char *XPLORStr,int frame)
{
  int ok=true;
  int isNew = true;

  if(!I) 
	 isNew=true;
  else 
	 isNew=false;
  if(ok) {
	 if(isNew) {
		I=(ObjectMap*)ObjectMapNew();
		isNew = true;
	 } else {
		isNew = false;
	 }
    ObjectMapXPLORStrToMap(I,XPLORStr,frame);
    SceneChanged();
  }
  return(I);
}
/*========================================================================*/
ObjectMap *ObjectMapReadCCP4Str(ObjectMap *I,char *XPLORStr,int bytes,int frame)
{
  int ok=true;
  int isNew = true;

  if(!I) 
	 isNew=true;
  else 
	 isNew=false;
  if(ok) {
	 if(isNew) {
		I=(ObjectMap*)ObjectMapNew();
		isNew = true;
	 } else {
		isNew = false;
	 }
    ObjectMapCCP4StrToMap(I,XPLORStr,bytes,frame);
    SceneChanged();
  }
  return(I);
}
/*========================================================================*/
ObjectMap *ObjectMapLoadCCP4File(ObjectMap *obj,char *fname,int frame)
{
  ObjectMap *I = NULL;
  int ok=true;
  FILE *f;
  long size;
  char *buffer,*p;
  float mat[9];

  f=fopen(fname,"r");
  if(!f)
	 ok=ErrMessage("ObjectMapLoadCCP4File","Unable to open file!");
  else
	 {
		if(Feedback(FB_ObjectMap,FB_Actions))
		  {
			printf(" ObjectMapLoadCCP4File: Loading from '%s'.\n",fname);
		  }
		
		fseek(f,0,SEEK_END);
      size=ftell(f);
		fseek(f,0,SEEK_SET);

		buffer=(char*)mmalloc(size);
		ErrChkPtr(buffer);
		p=buffer;
		fseek(f,0,SEEK_SET);
		fread(p,size,1,f);
		fclose(f);

		I=ObjectMapReadCCP4Str(obj,buffer,size,frame);

		mfree(buffer);
      CrystalDump(I->Crystal);
      multiply33f33f(I->Crystal->FracToReal,I->Crystal->RealToFrac,mat);
    }
#ifdef _UNDEFINED
  for(a=0;a<9;a++)
    printf("%10.5f\n",mat[a]);
#endif
  return(I);

}

/*========================================================================*/
ObjectMap *ObjectMapLoadXPLORFile(ObjectMap *obj,char *fname,int frame)
{
  ObjectMap *I = NULL;
  int ok=true;
  FILE *f;
  long size;
  char *buffer,*p;
  float mat[9];

  f=fopen(fname,"r");
  if(!f)
	 ok=ErrMessage("ObjectMapLoadXPLORFile","Unable to open file!");
  else
	 {
		if(Feedback(FB_ObjectMap,FB_Actions))
		  {
			printf(" ObjectMapLoadXPLORFile: Loading from '%s'.\n",fname);
		  }
		
		fseek(f,0,SEEK_END);
      size=ftell(f);
		fseek(f,0,SEEK_SET);

		buffer=(char*)mmalloc(size+255);
		ErrChkPtr(buffer);
		p=buffer;
		fseek(f,0,SEEK_SET);
		fread(p,size,1,f);
		p[size]=0;
		fclose(f);

		I=ObjectMapReadXPLORStr(obj,buffer,frame);

		mfree(buffer);
      CrystalDump(I->Crystal);
      multiply33f33f(I->Crystal->FracToReal,I->Crystal->RealToFrac,mat);
    }
#ifdef _UNDEFINED
  for(a=0;a<9;a++)
    printf("%10.5f\n",mat[a]);
#endif
  return(I);

}

/*========================================================================*/
int ObjectMapNumPyArrayToMap(ObjectMap *I,PyObject *ary) {
  
  int a,b,c,d,e;
  float v[3],dens,maxd,mind;
  int ok = true;

#ifdef _PYMOL_NUMPY
  MyArrayObject *pao;
  pao = (MyArrayObject*)ary;
#endif
  maxd = FLT_MIN;
  mind = FLT_MAX;
  if(ok) {
    I->FDim[0]=I->Dim[0];
    I->FDim[1]=I->Dim[1];
    I->FDim[2]=I->Dim[2];
    I->FDim[3]=3;

    if(!(I->FDim[0]&&I->FDim[1]&&I->FDim[2])) 
      ok=false;
    else {
      I->Field=IsosurfFieldAlloc(I->FDim);
      for(c=0;c<I->FDim[2];c++)
        {
          v[2]=I->Origin[2]+I->Grid[2]*c;
          for(b=0;b<I->FDim[1];b++) {
            v[1]=I->Origin[1]+I->Grid[1]*b;
            for(a=0;a<I->FDim[0];a++) {
              v[0]=I->Origin[0]+I->Grid[0]*a;
#ifdef _PYMOL_NUMPY
              dens = *((double*)
                (pao->data+
                 (pao->strides[0]*a)+
                 (pao->strides[1]*b)+
                 (pao->strides[2]*c)));
#else
              dens = 0.0;
#endif
              F3(I->Field->data,a,b,c) = dens;
              if(maxd<dens) maxd = dens;
              if(mind>dens) mind = dens;
              for(e=0;e<3;e++) 
                F4(I->Field->points,a,b,c,e) = v[e];
            }
          }
        }
      d = 0;
      for(c=0;c<I->FDim[2];c+=(I->FDim[2]-1))
        {
          v[2]=I->Origin[2]+I->Grid[2]*c;
          for(b=0;b<I->FDim[1];b+=(I->FDim[1]-1)) {
            v[1]=I->Origin[1]+I->Grid[1]*b;
            for(a=0;a<I->FDim[0];a+=(I->FDim[0]-1)) {
              v[0]=I->Origin[0]+I->Grid[0]*a;
              copy3f(v,I->Corner[d]);
              d++;
            }
          }
        }
    }
  }
  if(ok) {
    copy3f(I->Origin,I->Obj.ExtentMin);
    copy3f(I->Origin,I->Obj.ExtentMax);
    add3f(I->Range,I->Obj.ExtentMax,I->Obj.ExtentMax);
    I->Obj.ExtentFlag=true;
  }
  if(!ok) {
    ErrMessage("ObjectMap","Error reading map");
  } else {
    if(Feedback(FB_ObjectMap,FB_Actions)) {
      printf(" ObjectMap: Map Read.  Range = %5.3f to %5.3f\n",mind,maxd);
    }
  }
  return(ok);
}
/*========================================================================*/
ObjectMap *ObjectMapLoadChemPyBrick(ObjectMap *I,PyObject *Map,
                                           int frame,int discrete)
{
  int ok=true;
  int isNew = true;
  PyObject *tmp;


  if(!I) 
	 isNew=true;
  else 
	 isNew=false;

  if(ok) {

	 if(isNew) {
		I=(ObjectMap*)ObjectMapNew();
		isNew = true;
	 } else {
		isNew = false;
	 }
    if(PyObject_HasAttrString(Map,"origin")&&
       PyObject_HasAttrString(Map,"dim")&&
       PyObject_HasAttrString(Map,"range")&&
       PyObject_HasAttrString(Map,"grid")&&
       PyObject_HasAttrString(Map,"lvl"))
      {
        tmp = PyObject_GetAttrString(Map,"origin");
        if(tmp) {
          PConvPyListToFloatArray(tmp,&I->Origin);
          Py_DECREF(tmp);
        } else 
          ok=ErrMessage("ObjectMap","missing brick origin.");
        tmp = PyObject_GetAttrString(Map,"dim");
        if(tmp) {
          PConvPyListToIntArray(tmp,&I->Dim);
          Py_DECREF(tmp);
        } else 
          ok=ErrMessage("ObjectMap","missing brick dimension.");
        tmp = PyObject_GetAttrString(Map,"range");
        if(tmp) {
          PConvPyListToFloatArray(tmp,&I->Range);
          Py_DECREF(tmp);
        } else 
          ok=ErrMessage("ObjectMap","missing brick range.");
        tmp = PyObject_GetAttrString(Map,"grid");
        if(tmp) {
          PConvPyListToFloatArray(tmp,&I->Grid);
          Py_DECREF(tmp);
        } else
          ok=ErrMessage("ObjectMap","missing brick grid.");
        tmp = PyObject_GetAttrString(Map,"lvl");
        if(tmp) {
          ObjectMapNumPyArrayToMap(I,tmp);	 
          Py_DECREF(tmp);
        } else
          ok=ErrMessage("ObjectMap","missing brick density.");

      }
    SceneChanged();
  }
  return(I);
}

/*========================================================================*/
ObjectMap *ObjectMapLoadChemPyMap(ObjectMap *I,PyObject *Map,
                                  int frame,int discrete)
{

  int ok=true;
  int isNew = true;
  float *cobj;
  WordType format;
  float v[3],vr[3],dens,maxd,mind;
  int a,b,c,d,e;


  /*  
  double test[1000];
  for(a=0;a<1000;a++) {
      test[a]=rand()/(1.0+INT_MAX);
      }
      PyObject_SetAttrString(Map,"c_object",
      PyCObject_FromVoidPtr(test,NULL));
  */
  maxd = FLT_MIN;
  mind = FLT_MAX;

  if(!I) 
	 isNew=true;
  else 
	 isNew=false;

  if(ok) {

	 if(isNew) {
		I=(ObjectMap*)ObjectMapNew();
		isNew = true;
	 } else {
		isNew = false;
	 }

    if(!PConvAttrToStrMaxLen(Map,"format",format,sizeof(WordType)-1))
      ok=ErrMessage("LoadChemPyMap","bad 'format' parameter.");
    else if(!PConvAttrToFloatArrayInPlace(Map,"cell_dim",I->Crystal->Dim,3))
      ok=ErrMessage("LoadChemPyMap","bad 'cell_dim' parameter.");
    else if(!PConvAttrToFloatArrayInPlace(Map,"cell_ang",I->Crystal->Angle,3))
      ok=ErrMessage("LoadChemPyMap","bad 'cell_ang' parameter.");
    else if(!PConvAttrToIntArrayInPlace(Map,"cell_div",I->Div,3))
      ok=ErrMessage("LoadChemPyMap","bad 'cell_div' parameter.");
    else if(!PConvAttrToIntArrayInPlace(Map,"first",I->Min,3))
      ok=ErrMessage("LoadChemPyMap","bad 'first' parameter.");
    else if(!PConvAttrToIntArrayInPlace(Map,"last",I->Max,3))
      ok=ErrMessage("LoadChemPyMap","bad 'last' parameter.");

    if(ok) {
      if (strcmp(format,"CObjectZYXfloat")==0) {
        ok = PConvAttrToPtr(Map,"c_object",(void**)&cobj);
        if(!ok)
          ErrMessage("LoadChemPyMap","CObject unreadable.");        
      } else {
        ok=ErrMessage("LoadChemPyMap","unsupported format.");        
      }
    }
    /* good to go */

    if(ok) {
      if (strcmp(format,"CObjectZYXfloat")==0) {

        I->FDim[0]=I->Max[0]-I->Min[0]+1;
        I->FDim[1]=I->Max[1]-I->Min[1]+1;
        I->FDim[2]=I->Max[2]-I->Min[2]+1;
        if(Feedback(FB_ObjectMap,FB_Actions)) {
          printf(" LoadChemPyMap: CObjectZYXdouble %dx%dx%d\n",I->FDim[0],I->FDim[1],I->FDim[2]);        
        }
        I->FDim[3]=3;
        if(!(I->FDim[0]&&I->FDim[1]&&I->FDim[2])) 
          ok=false;
        else {
          CrystalUpdate(I->Crystal);
          I->Field=IsosurfFieldAlloc(I->FDim);
          for(c=0;c<I->FDim[2];c++)
            {
              v[2]=(c+I->Min[2])/((float)I->Div[2]);
              for(b=0;b<I->FDim[1];b++) {
                v[1]=(b+I->Min[1])/((float)I->Div[1]);
                for(a=0;a<I->FDim[0];a++) {
                  v[0]=(a+I->Min[0])/((float)I->Div[0]);
                  
                  dens = *(cobj++);

                  F3(I->Field->data,a,b,c) = dens;
                  if(maxd<dens) maxd = dens;
                  if(mind>dens) mind = dens;
                  transform33f3f(I->Crystal->FracToReal,v,vr);
                  for(e=0;e<3;e++) 
                    F4(I->Field->points,a,b,c,e) = vr[e];
                }
              }
            }

          if(ok) {
            d = 0;
            for(c=0;c<I->FDim[2];c+=(I->FDim[2]-1))
              {
                v[2]=(c+I->Min[2])/((float)I->Div[2]);
                for(b=0;b<I->FDim[1];b+=(I->FDim[1]-1)) {
                  v[1]=(b+I->Min[1])/((float)I->Div[1]);
                  for(a=0;a<I->FDim[0];a+=(I->FDim[0]-1)) {
                    v[0]=(a+I->Min[0])/((float)I->Div[0]);
                    transform33f3f(I->Crystal->FracToReal,v,vr);
                    copy3f(vr,I->Corner[d]);
                    d++;
                  }
                }
              }
          }
        }
      }
    }
    
    if(ok) {
      CrystalDump(I->Crystal);
      
      v[2]=(I->Min[2])/((float)I->Div[2]);
      v[1]=(I->Min[1])/((float)I->Div[1]);
      v[0]=(I->Min[0])/((float)I->Div[0]);
      
      transform33f3f(I->Crystal->FracToReal,v,I->Obj.ExtentMin);
      
      v[2]=((I->FDim[2]-1)+I->Min[2])/((float)I->Div[2]);
      v[1]=((I->FDim[1]-1)+I->Min[1])/((float)I->Div[1]);
      v[0]=((I->FDim[0]-1)+I->Min[0])/((float)I->Div[0]);
      
      transform33f3f(I->Crystal->FracToReal,v,I->Obj.ExtentMax);
      I->Obj.ExtentFlag=true;
    }

    if(!ok) {
      ErrMessage("ObjectMap","Error reading map");
    } else {
		if(Feedback(FB_ObjectMap,FB_Actions)) {
        printf(" ObjectMap: Map Read.  Range = %5.3f to %5.3f\n",mind,maxd);
      }
    }

    if(ok) SceneChanged();
  }
  return(I);
}

