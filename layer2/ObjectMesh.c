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
#include<stdio.h>

#include"OOMac.h"
#include"ObjectMesh.h"
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

ObjectMesh *ObjectMeshNew(void);

static void ObjectMeshFree(ObjectMesh *I);

static void ObjectMeshFree(ObjectMesh *I) {
  VLAFreeP(I->N);
  VLAFreeP(I->V);
  OOFreeP(I);
}

void ObjectMeshDump(ObjectMesh *I,char *fname)
{
  float *v=I->V;
  int *n=I->N;
  int c;
  FILE *f;
  OrthoLineType buf;

  f=fopen(fname,"w");
  if(!f) 
    ErrMessage("ObjectMeshDump","can't open file for writing");
  else {
    if(n)
      while(*n)
        {
          c=*(n++);
          while(c--) {
            fprintf(f,"%10.4f%10.4f%10.4f\n",v[0],v[1],v[2]);
            v+=3;
          }
        }
    fclose(f);
    sprintf(buf,"%s written to %s\n",I->Obj.Name,fname);
    ErrOk("ObjectMeshDump",buf);
  }
}

static void ObjectMeshUpdate(ObjectMesh *I) {
  if(I->ResurfaceFlag) {
    I->ResurfaceFlag=false;
    PRINTF " ObjectMesh: updating \"%s\".\n" , I->Obj.Name ENDF
    if(I->Map->Field) IsosurfVolume(I->Map->Field,I->Level,&I->N,&I->V,I->Range,I->DotFlag); 
  }
  SceneDirty();
}

static void ObjectMeshRender(ObjectMesh *I,int frame,CRay *ray,Pickable **pick)
{
  float *v=I->V;
  float *vc;
  int *n=I->N;
  int c;

  if(ray) {
	 if(n) {
      vc = ColorGet(I->Obj.Color);
      if(I->DotFlag) {
        ray->fColor3fv(ray,vc);
        while(*n)
		  {
			 c=*(n++);
			 if(c--)
				{
				  v+=3;
				  while(c--)
					 {
						ray->fSphere3fv(ray,v,0.05);
						v+=3;
					 }
				}
		  }
      } else {
        while(*n)
          {
            c=*(n++);
            if(c--)
              {
                v+=3;
                while(c--)
                  {
                    ray->fCylinder3fv(ray,v-3,v,0.05,vc,vc);
                    v+=3;
                  }
              }
          }
      }
	 }
  } else if(pick&&PMGUI) {
  } else if(PMGUI) {
	 if(n) {
      ObjectUseColor(&I->Obj);
      while(*n)
        {
          c=*(n++);
          if(I->DotFlag) 
            glBegin(GL_POINTS);
          else 
            glBegin(GL_LINE_STRIP);
          SceneResetNormal(false);
          while(c--) {
            glVertex3fv(v);
            v+=3;
          }
          glEnd();
        }
    }
  }
}

/*========================================================================*/
ObjectMesh *ObjectMeshNew(void)
{
  OOAlloc(ObjectMesh);

  ObjectInit((Object*)I);

  I->V = VLAlloc(float,10000);
  I->N = VLAlloc(int,10000);
  I->N[0]=0;
  I->Obj.type = cObjectMesh;

  I->ResurfaceFlag=true;

  I->Obj.fFree = (void (*)(struct Object *))ObjectMeshFree;
  I->Obj.fUpdate =  (void (*)(struct Object *)) ObjectMeshUpdate;
  I->Obj.fRender =(void (*)(struct Object *, int, CRay *, Pickable **))ObjectMeshRender;

#ifdef _NOT_YET_NEEDED
  I->Obj.fGetNFrame = (int (*)(struct Object *)) ObjectMeshGetNFrames;
#endif
  
  return(I);
}


/*========================================================================*/
ObjectMesh *ObjectMeshFromBox(ObjectMap *map,float *mn,float *mx,float level,int dotFlag)
{
  ObjectMesh *I;
  I=ObjectMeshNew();
  I->Map = map;
  I->Level = level;
  I->DotFlag = dotFlag;
  IsosurfGetRange(I->Map->Field,I->Map->Crystal,mn,mx,I->Range);

  /*  printf("Brick %d %d %d %d %d %d\n",I->Range[0],I->Range[1],I->Range[2],I->Range[3],I->Range[4],I->Range[5]);*/
  SceneChanged();
  return(I);
}


