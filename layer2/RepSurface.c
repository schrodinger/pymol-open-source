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
#include<GL/gl.h>
#include<values.h>
#include<math.h>
#include"Base.h"
#include"MemoryDebug.h"
#include"OOMac.h"
#include"RepSurface.h"
#include"Map.h"
#include"Scene.h"
#include"Sphere.h"
#include"Setting.h"
#include"Color.h"
#include"ObjectMolecule.h"
#include"Triangle.h"

typedef struct RepSurface {
  Rep R;
  int N;
  int NT;
  float *V,*VN;
  int *T;
  int NDot;
  float *Dot;
  int solidFlag;
  Object *Obj;
} RepSurface;


void RepSurfaceRender(RepSurface *I,CRay *ray,Pickable **pick);
void RepSurfaceFree(RepSurface *I);

void RepSurfaceFree(RepSurface *I)
{
  FreeP(I->V);
  FreeP(I->VN);
  VLAFreeP(I->T);
  /*  VLAFreeP(I->N);*/
  OOFreeP(I);
}

void RepSurfaceGetSolventDots(RepSurface *I,CoordSet *cs,float probe_radius,SphereRec *sp);

void RepSurfaceRender(RepSurface *I,CRay *ray,Pickable **pick)
{
  float *v=I->V;
  float *vn=I->VN;
  int *t=I->T;
  int c=I->N;
  
  if(ray) {
  } else if(pick) {
  } else {
	 c=I->NT;
    if(c) {
      glBegin(GL_TRIANGLES);
      while(c--)
        {
			 
          glNormal3fv(vn+(*t)*3);
          glVertex3fv(v+(*t)*3);
			 t++;
          glNormal3fv(vn+(*t)*3);
          glVertex3fv(v+(*t)*3);
			 t++;
          glNormal3fv(vn+(*t)*3);
          glVertex3fv(v+(*t)*3);
			 t++;
        }
		glEnd();

		c=I->N;
      glBegin(GL_POINTS);
      while(c--)
        {
          glNormal3fv(vn);
          glVertex3fv(v);
          v+=3;
          vn+=3;
        }
		  glEnd();
    }
  }
}

#define solv_tole 0.02
#define minimum_sep 0.5

Rep *RepSurfaceNew(CoordSet *cs)
{
  ObjectMolecule *obj;
  int a,b,h,i,j,k,l,c;
  MapType *map,*solv_map;
  /* grid */
  float v1[3],*v0,*v,*vn,*vn0;
  float vdw;
  int SurfaceFlag = false;
  float probe_radius,probe_radius2;
  int inclH = true;
  int cullByFlag = false;
  int flag,*dot_flag,*p;
  int ds,a1,a2,c1;
  SphereRec *sp = Sphere0;
  OOAlloc(RepSurface);

 /* get current dot sampling */
  ds = (int)SettingGet(cSetting_dot_density);
  if(ds<0) ds=0;
  switch(ds) {
  case 0: sp=Sphere0; break;
  case 1: sp=Sphere1; break;
  case 2: sp=Sphere2; break;
  default: sp=Sphere3; break;
  }
  sp=Sphere1;

  cullByFlag = SettingGet(cSetting_trim_dots);
  inclH = SettingGet(cSetting_dot_hydrogens);

  probe_radius = SettingGet(cSetting_solvent_radius);
  probe_radius2 = probe_radius*probe_radius;

  I->N=0;
  I->NT=0;
  I->V=NULL;
  I->VN=NULL;
  I->T=NULL;
  I->R.fRender=(void (*)(struct Rep *, CRay *, Pickable **))RepSurfaceRender;
  I->R.fFree=(void (*)(struct Rep *))RepSurfaceFree;
  I->Obj = (Object*)(cs->Obj);
  I->Dot=NULL;

  obj = cs->Obj;

  /* don't waist time computing a Surface unless we need it!! */
  for(a=0;a<cs->NIndex;a++) {
	 if(obj->AtomInfo[cs->IdxToAtm[a]].visRep[cRepSurface])
		{
		  SurfaceFlag=true;
		  break;
		}
  }
  if(SurfaceFlag) {

	 I->V=Alloc(float,cs->NIndex*3*sp->nDot*2);
    ErrChkPtr(I->V);
	 I->VN=Alloc(float,cs->NIndex*3*sp->nDot*2);


    ErrChkPtr(I->VN);
	 I->N=0;
    v=I->V;
    vn=I->VN;
	 OrthoBusyFast(0,1);
    
	 RepSurfaceGetSolventDots(I,cs,probe_radius,Sphere2);
	 solv_map=MapNew(probe_radius,I->Dot,I->NDot,NULL);
	 map=MapNew(MAX_VDW+probe_radius,cs->Coord,cs->NIndex,NULL);
	 if(map&&solv_map)
		{
		  MapSetupExpress(solv_map);
		  MapSetupExpress(map);
        for(a=0;a<cs->NIndex;a++)
          {
            a1 = cs->IdxToAtm[a];
            if(obj->AtomInfo[a1].visRep[cRepSurface])
              if((inclH||(obj->AtomInfo[a1].name[0]!='H'))&&
                 ((!cullByFlag)||(!(obj->AtomInfo[a1].customFlag&0x2)))) {
                /* ignore if the "2" bit is set */
                c1=*(cs->Color+a);
                v0 = cs->Coord+3*a;
                vdw = cs->Obj->AtomInfo[a1].vdw;

                for(b=0;b<sp->nDot;b++)
                  {
                    v1[0]=v0[0]+vdw*sp->dot[b].v[0];
                    v1[1]=v0[1]+vdw*sp->dot[b].v[1];
                    v1[2]=v0[2]+vdw*sp->dot[b].v[2];
                    
                    MapLocus(map,v1,&h,&k,&l);
                    
                    flag=true;
                    
                    i=*(MapEStart(map,h,k,l));
                    if(i) {
                      j=map->EList[i++];
                      while(j>=0) {
                        a2 = cs->IdxToAtm[j];
                        if((inclH||(obj->AtomInfo[a2].name[0]!='H'))&&
                           ((!cullByFlag)||(!(obj->AtomInfo[a2].customFlag&0x2))))  
                          /* ignore if the "2" bit is set */
                          if(j!=a)
                            if(within3f(cs->Coord+3*j,v1,cs->Obj->AtomInfo[a2].vdw))
                              {
                                flag=false;
                                break;
                              }
                        j=map->EList[i++];
                      }
                    }
                    if(flag) /* doesn't intersect any atom */
                      {
                        MapLocus(solv_map,v1,&h,&k,&l);
                        
                        flag=true;
                        i=*(MapEStart(solv_map,h,k,l));
                        if(i) {
                          j=solv_map->EList[i++];
                          while(j>=0) {
                            if(within3f(I->Dot+3*j,v1,probe_radius+solv_tole))
                              {
                                flag=false;
                                break;
                              }
                            j=solv_map->EList[i++];
                          }
                        }
                        if(!flag) /* appropriately close to a water center */
                          {
                            *(v++)=v1[0];
                            *(v++)=v1[1];
                            *(v++)=v1[2];
                            *(vn++)=sp->dot[b].v[0];
                            *(vn++)=sp->dot[b].v[1];
                            *(vn++)=sp->dot[b].v[2];
                            I->N++;
                          }
                      }
                  }
              }
          } 
        MapFree(solv_map);
        FreeP(I->Dot);	 
        RepSurfaceGetSolventDots(I,cs,probe_radius,Sphere2);
		  if(I->NDot) {
			 solv_map=MapNew(probe_radius,I->Dot,I->NDot,NULL);
			 /* concave surfaces */
			 if(solv_map)
				{
				  MapSetupExpress(solv_map);
				  for(a=0;a<I->NDot;a++)
					 {
						v0 = I->Dot+3*a;
						vdw = probe_radius;
						for(b=0;b<sp->nDot;b++)
						  {
							 v[0]=v0[0]+vdw*sp->dot[b].v[0];
							 v[1]=v0[1]+vdw*sp->dot[b].v[1];
							 v[2]=v0[2]+vdw*sp->dot[b].v[2];
							 MapLocus(solv_map,v,&h,&k,&l);
							 flag=true;
							 i=*(MapEStart(solv_map,h,k,l));
							 if(i) {
								j=solv_map->EList[i++];
								while(j>=0) {
								  if(j!=a) 
									 {
										if(within3f(I->Dot+3*j,v,probe_radius)) {
										  flag=false;
										  break;
										}
									 }
								  j=solv_map->EList[i++];
								}
							 }
							 if(flag)
								{
								  MapLocus(map,v,&h,&k,&l);
								  flag=true;
								  i=*(MapEStart(map,h,k,l));
								  if(i) {
									 j=map->EList[i++];
									 while(j>=0) {
										a2 = cs->IdxToAtm[j];
										if((inclH||(obj->AtomInfo[a2].name[0]!='H'))&&
											((!cullByFlag)||(!(obj->AtomInfo[a2].customFlag&0x2))))  
										  /* ignore if the "2" bit is set */
										  if(j!=a)
											 if(within3f(cs->Coord+3*j,v,
															 cs->Obj->AtomInfo[a2].vdw+probe_radius+solv_tole))
												{
												  flag=false;
												  break;
												}
										j=map->EList[i++];
									 }
								  }
								  if(!flag) {
									 I->N++;
									 v+=3;
									 *(vn++)=-sp->dot[b].v[0];
									 *(vn++)=-sp->dot[b].v[1];
									 *(vn++)=-sp->dot[b].v[2];
								  }
								}
						  }
					 }
				  MapFree(solv_map);
				}
			 MapFree(map);
		  }
		  printf("%i\n",I->N);
		  FreeP(I->Dot);	 
		}
	 /* now, eliminate dots that are too close to each other*/
	 if(I->N) 
		{
		  dot_flag=Alloc(int,I->N);
		  ErrChkPtr(dot_flag);
		  for(a=0;a<I->N;a++) {
			 dot_flag[a]=1;
		  }
		  
		  map=MapNew(minimum_sep,I->V,I->N,NULL);
		  if(map)
			 {
				MapSetupExpress(map);		  
				v=I->V;
				for(a=0;a<I->N;a++) {
				  if(dot_flag[a]) {
					 MapLocus(map,v,&h,&k,&l);
					 i=*(MapEStart(map,h,k,l));
					 if(i) {
						j=map->EList[i++];
						while(j>=0) {
						  if(j!=a) 
							 {
								if(dot_flag[j]) {
								  if(within3f(I->V+(3*j),v,minimum_sep)) {
									 dot_flag[j]=0;
									 break;
								  }
								}
							 }
						  j=map->EList[i++];
						}
					 }
				  }
				  v+=3;
				}
				MapFree(map);
			 }
		  v0=I->V;
		  v=I->V;
		  vn0=I->VN;
		  vn=I->VN;
		  p=dot_flag;
		  c=I->N;
		  I->N=0;
		  for(a=0;a<c;a++)
			 {
				if(*(p++)) {
				  *(v0++)=*(v++);
				  *(v0++)=*(v++);
				  *(v0++)=*(v++);
				  *(vn0++)=*(vn++);
				  *(vn0++)=*(vn++);
				  *(vn0++)=*(vn++);
				  I->N++;
				} else {
				  v+=3;
				  vn+=3;
				}
			 }
		  FreeP(dot_flag);
		}
	 if(I->N) {	
		I->V = Realloc(I->V,float,(v0-I->V));
		I->VN = Realloc(I->VN,float,(vn0-I->VN));
		I->T=TrianglePointsToSurface(I->V,I->VN,I->N,probe_radius,&I->NT);

	 } else {
		I->V = Realloc(I->V,float,1);
		I->VN = Realloc(I->VN,float,1);
	 }

  }
  OrthoBusyFast(4,4);
  return((void*)(struct Rep*)I);
}

void RepSurfaceGetSolventDots(RepSurface *I,CoordSet *cs,float probe_radius,SphereRec *sp)
{
  ObjectMolecule *obj;
  int a,b,c,a1,a2,flag,i,h,k,l,j;
  float *v,*v0,vdw;
  MapType *map;
  int inFlag,*p,*dot_flag;
  int cavity_cull;
  float probe_radius_plus;
  int dotCnt,maxCnt,maxDot=0;
  int cnt;
  cavity_cull = (int)SettingGet(cSetting_cavity_cull);

  obj = cs->Obj;
  I->Dot=(float*)mmalloc(sizeof(float)*cs->NIndex*3*sp->nDot);
  ErrChkPtr(I->Dot);

  probe_radius_plus = probe_radius * 1.5;

  I->NDot=0;
  map=MapNew(MAX_VDW+probe_radius,cs->Coord,cs->NIndex,NULL);
  if(map)
	 {
		MapSetupExpress(map);
		maxCnt=0;
		v=I->Dot;
		for(a=0;a<cs->NIndex;a++)
		  {
			 OrthoBusyFast(a,cs->NIndex*3);
			 dotCnt=0;
			 a1 = cs->IdxToAtm[a];
			 v0 = cs->Coord+3*a;
			 vdw = cs->Obj->AtomInfo[a1].vdw+probe_radius;
          for(b=0;b<sp->nDot;b++)
            {
              v[0]=v0[0]+vdw*sp->dot[b].v[0];
              v[1]=v0[1]+vdw*sp->dot[b].v[1];
              v[2]=v0[2]+vdw*sp->dot[b].v[2];
              MapLocus(map,v,&h,&k,&l);
              flag=true;
              i=*(MapEStart(map,h,k,l));
              if(i) {
                j=map->EList[i++];
                while(j>=0) {
                  if(j!=a) 
                    {
                      a2 = cs->IdxToAtm[j];
                      if(within3f(cs->Coord+3*j,v,cs->Obj->AtomInfo[a2].vdw+probe_radius)) {
                        flag=false;
                        break;
                      }
                    }
                  j=map->EList[i++];
                }
              }
              if(flag)
                {
                  dotCnt++;
                  v+=3;
                  I->NDot++;
                }
            }
			 if(dotCnt>maxCnt)
				{
				  maxCnt=dotCnt;
				  maxDot=I->NDot-1;
				}
		  }
		MapFree(map);
	 }

  if(cavity_cull>0) {
	 dot_flag=Alloc(int,I->NDot);
	 ErrChkPtr(dot_flag);
	 for(a=0;a<I->NDot;a++) {
		dot_flag[a]=0;
	 }
	 dot_flag[maxDot]=1; /* this guarantees that we have a valid dot */

	 map=MapNew(probe_radius_plus,I->Dot,I->NDot,NULL);
	 if(map)
		{
		  MapSetupExpress(map);		  
		  flag=true;
		  while(flag) {
			 p=dot_flag;
			 v=I->Dot;
		  
			 flag=false;
			 for(a=0;a<I->NDot;a++)
				{
				  if(!dot_flag[a]) {
					 cnt=0;
					 MapLocus(map,v,&h,&k,&l);
					 i=*(MapEStart(map,h,k,l));
					 if(i) {
						j=map->EList[i++];
						while(j>=0) {
						  if(j!=a) 
							 {
								if(within3f(I->Dot+(3*j),v,probe_radius_plus)) {
								  if(dot_flag[j]) {
									 *p=true;
									 flag=true;
									 break;
								  }
								  cnt++;
								  if(cnt>cavity_cull) 
									 {
										*p=true;
										flag=true;
										break;
									 }
								}
							 }
						  j=map->EList[i++];
						}
					 }
				  }
				  v+=3;
				  p++;
				}
		  }
		  MapFree(map);
		}
	 v0=I->Dot;
	 v=I->Dot;
	 p=dot_flag;
	 c=I->NDot;
	 I->NDot=0;
	 for(a=0;a<c;a++)
		{
		  if(*(p++)) {
			 *(v0++)=*(v++);
			 *(v0++)=*(v++);
			 *(v0++)=*(v++);
			 I->NDot++;
		  } else {
			 v+=3;
		  }
		}
	 FreeP(dot_flag);
  }
  /*  printf("NDot: %d->%d\n",c,I->NDot);*/
}


