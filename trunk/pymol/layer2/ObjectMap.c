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

#include<string.h>
#include<GL/gl.h>

#include"OOMac.h"
#include"ObjectMap.h"
#include"Base.h"
#include"MemoryDebug.h"
#include"Map.h"
#include"Debug.h"
#include"Parse.h"
#include"Isosurf.h"
#include"Vector.h"
#include"Color.h"
#include"main.h"
#include"Scene.h"

static void ObjectMapFree(ObjectMap *I);

static void ObjectMapFree(ObjectMap *I) {
  if(I->Field) {
    IsosurfFieldFree(I->Field);
    I->Field=NULL;
  }
  VLAFreeP(I->N);
  VLAFreeP(I->V);
  OOFreeP(I->Crystal);
  OOFreeP(I);
}

static void ObjectMapUpdate(ObjectMap *I) {
  if(I->ResurfaceFlag) {
    I->ResurfaceFlag=false;
    if(I->Field) 
      IsosurfVolume(I->Field,1.0,&I->N,&I->V,NULL,NULL);  
  }
}

static void ObjectMapRender(ObjectMap *I,int frame,CRay *ray,Pickable **pick)
{
  float *v=I->V;
  float *vc;/*=I->VC;*/
  int *n=I->N;
  int c;

  if(ray) {
	 if(n) {
		ray->fColor3fv(ray,ColorGet(I->Obj.Color));
		while(*n)
		  {
			 c=*(n++);
			 if(c--)
				{
				  vc+=3;
				  v+=3;
				  while(c--)
					 {
						ray->fCylinder3fv(ray,v-3,v,0.05,vc-3,vc);
						v+=3;
						vc+=3;
					 }
				}
		  }
	 }
  } else if(pick&&PMGUI) {
  } else if(PMGUI) {
	 if(n) {
		if(true) {
		  while(*n)
			 {
				c=*(n++);
				glBegin(GL_LINE_STRIP);
				SceneResetNormal(false);
				while(c--) {
				  glVertex3fv(v);
				  v+=3;
				}
				glEnd();
			 }
		} else {
		  while(*n)
			 {
				c=*(n++);
				glBegin(GL_LINE_STRIP);
				SceneResetNormal(false);
				while(c--) {
				  glColor3fv(vc);
				  vc+=3;
				  glVertex3fv(v);
				  v+=3;
				}
				glEnd();
			 }
		}
	 }
  }
}

/*========================================================================*/
ObjectMap *ObjectMapNew(void)
{
  OOAlloc(ObjectMap);

  ObjectInit((Object*)I);

  I->Crystal = CrystalNew();
  I->Field = NULL;
  I->V = VLAlloc(float,10000);
  I->N = VLAlloc(int,10000);

  I->ResurfaceFlag=true;

  I->Obj.fFree = (void (*)(struct Object *))ObjectMapFree;
  I->Obj.fUpdate =  (void (*)(struct Object *)) ObjectMapUpdate;
  I->Obj.fRender =(void (*)(struct Object *, int, CRay *, Pickable **))ObjectMapRender;

#ifdef _NOT_YET_NEEDED
  I->Obj.fGetNFrame = (int (*)(struct Object *)) ObjectMapGetNFrames;
#endif

  return(I);
}
/*========================================================================*/
int ObjectMapXPLORStrToMap(ObjectMap *I,char *XPLORStr,int frame) {
  
  unsigned char *p;
  int a,b,c,d,e;
  float v[3],vr[3];
  char cc[MAXLINELEN],tmp[MAXLINELEN];
  int n;
  int ok = true;

  p=XPLORStr;

  while(*p) {
    p = ParseNCopy(cc,p,8);
    if(!*cc) 
      p = ParseNextLine(p);
    else if(sscanf(cc,"%i",&n)==1) {
      p=ParseWordCopy(cc,p,MAXLINELEN);
      if(strstr(cc,"!NTITLE")) {
        p=ParseNextLine(p);
        while(n--) {
          p=ParseNextLine(p);          
        } 
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
#ifdef MAPDEBUG
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
#endif

  } else 
    ok=false;
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
          d=0;
          for(b=0;b<I->FDim[1];b++) {
            v[1]=(b+I->Min[1])/((float)I->Div[1]);
            for(a=0;a<I->FDim[0];a++) {
              v[0]=(a+I->Min[0])/((float)I->Div[0]);
              p=ParseNCopy(cc,p,12);
              if(sscanf(cc,"%f",F3Ptr(I->Field->data,a,b,c,I->Field->dimensions))!=1) ok=false;
              transform33f3f(I->Crystal->FracToReal,v,vr);
              for(e=0;e<3;e++) 
                F4(I->Field->points,a,b,c,e,I->Field->dimensions) = vr[e];
              d++; if(d==6) { p=ParseNextLine(p); d=0;}              
            }
        }
      }
    }
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
    I->ResurfaceFlag=true;
  }
  return(I);
}
/*========================================================================*/
ObjectMap *ObjectMapLoadXPLORFile(ObjectMap *obj,char *fname,int frame)
{
  ObjectMap *I = NULL;
  int ok=true;
  FILE *f;
  fpos_t size;
  char *buffer,*p;
  float mat[9];
  int a;

  f=fopen(fname,"r");
  if(!f)
	 ok=ErrMessage("ObjectMapLoadXPLORFile","Unable to open file!");
  else
	 {
		if(DebugState&DebugMap)
		  {
			printf(" ObjectMapLoadXPLORFile: Loading from %s.\n",fname);
		  }
		
		fseek(f,0,SEEK_END);
		fgetpos(f,&size);
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
	 }
  CrystalDump(I->Crystal);
  multiply33f33f(I->Crystal->FracToReal,I->Crystal->RealToFrac,mat);
  for(a=0;a<9;a++)
    printf("%10.5f\n",mat[a]);

  return(I);

}




