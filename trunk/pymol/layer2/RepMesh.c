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
#include"os_gl.h"

#include"Base.h"
#include"MemoryDebug.h"
#include"OOMac.h"
#include"RepMesh.h"
#include"Map.h"
#include"Isosurf.h"
#include"Scene.h"
#include"Sphere.h"
#include"Setting.h"
#include"Color.h"
#include"main.h"

typedef struct RepMesh {
  Rep R;
  int *N;
  int NTot;
  float *V,*VC;
  int NDot;
  float *Dot;
  float Radius,Width;
  int oneColorFlag;
  int oneColor;
  CObject *Obj;
  int *LastVisib;
  int *LastColor;
} RepMesh;

#include"ObjectMolecule.h"

void RepMeshRender(RepMesh *I,CRay *ray,Pickable **pick);
void RepMeshFree(RepMesh *I);
void RepMeshColor(RepMesh *I,CoordSet *cs);
int RepMeshSameVis(RepMesh *I,CoordSet *cs);

void RepMeshInit(void)
{
  IsosurfInit();
}

void RepMeshFree(RepMesh *I)
{
  FreeP(I->VC);
  VLAFreeP(I->V);
  VLAFreeP(I->N);
  FreeP(I->LastColor);
  FreeP(I->LastVisib);
  OOFreeP(I);
}

void RepMeshGetSolventDots(RepMesh *I,CoordSet *cs,float *min,float *max,float probe_radius);

void RepMeshRender(RepMesh *I,CRay *ray,Pickable **pick)
{
  float *v=I->V;
  float *vc=I->VC;
  int *n=I->N;
  int c;
  float *col=NULL;

  if(ray) {
	 if(n) {
      if(I->oneColorFlag) 
        col=ColorGet(I->oneColor);
		ray->fColor3fv(ray,ColorGet(I->Obj->Color));
		while(*n)
		  {
			 c=*(n++);
			 if(c--)
				{
				  vc+=3;
				  v+=3;
              if(I->oneColorFlag) {
                while(c--)
                  {
                    ray->fCylinder3fv(ray,v-3,v,I->Radius,col,col);
                    v+=3;
                    vc+=3;
                  }
              } else {
                while(c--)
                  {
                    ray->fCylinder3fv(ray,v-3,v,I->Radius,vc-3,vc);
                    v+=3;
                    vc+=3;
                  }
              }
				}
		  }
	 }
  } else if(pick&&PMGUI) {
  } else if(PMGUI) {
    glLineWidth(I->Width);
    if(n) {
	   glDisable(GL_LIGHTING);
		if(I->oneColorFlag) {
		  while(*n)
			 {
            glColor3fv(ColorGet(I->oneColor));
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
	   glEnable(GL_LIGHTING);
	 }
     }
  }

int RepMeshSameVis(RepMesh *I,CoordSet *cs)
{
  int same = true;
  int *lv,*lc,*cc;
  int a;
  AtomInfoType *ai;

  ai = cs->Obj->AtomInfo;
  lv = I->LastVisib;
  lc = I->LastColor;
  cc = cs->Color;

  for(a=0;a<cs->NIndex;a++)
    {
      if(*(lv++)!=(ai + cs->IdxToAtm[a])->visRep[cRepMesh] ) {
        same=false;
        break;
      }
      if(*(lc++)!=*(cc++)) {
        same=false;
        break;
      }
    }
  return(same);
}


void RepMeshColor(RepMesh *I,CoordSet *cs)
{
  MapType *map;
  int a,i0,i,j,h,k,l,c1;
  float *v0,*vc,*c0;
  int *lv,*lc,*cc;
  int first_color;
  ObjectMolecule *obj;
  float probe_radius;
  float dist,minDist;
  int inclH;
  int cullByFlag = false;
  int mesh_mode;
  int mesh_color;
  AtomInfoType *ai2;

  obj=cs->Obj;

  probe_radius = SettingGet_f(cs->Setting,obj->Obj.Setting,cSetting_solvent_radius);
  mesh_color = SettingGet_color(cs->Setting,obj->Obj.Setting,cSetting_mesh_color);
  mesh_mode = SettingGet_f(cs->Setting,obj->Obj.Setting,cSetting_mesh_mode);
  cullByFlag = (mesh_mode==cRepMesh_by_flags);
  inclH = !(mesh_mode==cRepMesh_heavy_atoms);

  if(!I->LastVisib) I->LastVisib = Alloc(int,cs->NIndex);
  if(!I->LastColor) I->LastColor = Alloc(int,cs->NIndex);
  lv = I->LastVisib;
  lc = I->LastColor;
  cc = cs->Color;
  ai2=obj->AtomInfo;
  for(a=0;a<cs->NIndex;a++)
    {
      *(lv++) = (ai2 + cs->IdxToAtm[a])->visRep[cRepMesh];
      *(lc++) = *(cc++);
    }

  I->Width = SettingGet_f(cs->Setting,obj->Obj.Setting,cSetting_mesh_width);
  I->Radius = SettingGet_f(cs->Setting,obj->Obj.Setting,cSetting_mesh_radius);

  if(I->NTot) {
	 obj=cs->Obj;
	 I->oneColorFlag=true;
	 first_color=-1;
	 if(!I->VC) I->VC = Alloc(float,3*I->NTot);
	 vc=I->VC;
	 /* now, assign colors to each point */
	 map=MapNew(MAX_VDW+probe_radius,cs->Coord,cs->NIndex,NULL);
	 if(map)
		{
		  MapSetupExpress(map);
		  for(a=0;a<I->NTot;a++)
			 {
				c1=1;
				minDist=MAXFLOAT;
				i0=-1;
				v0 = I->V+3*a;
				MapLocus(map,v0,&h,&k,&l);
				
				i=*(MapEStart(map,h,k,l));
				if(i) {
				  j=map->EList[i++];
				  while(j>=0) {
					 ai2=obj->AtomInfo+cs->IdxToAtm[j];
					 if((inclH||(!ai2->hydrogen))&&
						 ((!cullByFlag)||
                    (!(ai2->flags&cAtomFlag_ignore))))  
						{
						  dist = diff3f(v0,cs->Coord+j*3) - ai2->vdw;
						  if(dist<minDist)
							 {
								i0=j;
								minDist=dist;
							 }
						}
					 j=map->EList[i++];
				  }
				}
				if(i0>=0) {
				  c1=*(cs->Color+i0);
				  if(I->oneColorFlag) {
					 if(first_color>=0) {
						if(first_color!=c1)
						  I->oneColorFlag=false;
					 } else first_color=c1;
				  }
				}
				c0 = ColorGet(c1);
				*(vc++) = *(c0++);
				*(vc++) = *(c0++);
				*(vc++) = *(c0++);
			 }
		  MapFree(map);
		}
    if(I->oneColorFlag) {
      I->oneColor=first_color;
    }
  } 
  if(mesh_color>=0) {
    I->oneColorFlag=1;
    I->oneColor=mesh_color;
  }
  
}

Rep *RepMeshNew(CoordSet *cs)
{
  ObjectMolecule *obj;
  CoordSet *ccs;
  int a,b,c,d,h,k,l,*n;
  MapType *map,*smap;
  /* grid */
  Vector3f minE,maxE,sizeE;
  float size;
  int dims[3];
  float gridSize;
  Isofield *field;
  Vector3f point;
  float vLen,pVal,vdw;
  int aNear;
  float aLen;
  int cur;
  int meshFlag = false;
  int escFlag;
  float probe_radius,probe_radius2;
  float min_spacing;
  int visFlag;
  int mesh_mode;
  int cullByFlag;
  int inclH;
  AtomInfoType *ai1;

  OOAlloc(RepMesh);

  PRINTFD(FB_RepMesh)
	 " RepMeshNew-DEBUG: entered with coord-set %p\n",cs
	 ENDFD;
  obj = cs->Obj;

  mesh_mode = SettingGet_f(cs->Setting,obj->Obj.Setting,cSetting_mesh_mode);
  cullByFlag = (mesh_mode==cRepMesh_by_flags);
  inclH = !(mesh_mode==cRepMesh_heavy_atoms);

  visFlag=false;
  for(a=0;a<cs->NIndex;a++) {
    ai1 = obj->AtomInfo+cs->IdxToAtm[a];
	 if(ai1->visRep[cRepMesh]&&
       (inclH||(!ai1->hydrogen))&&
       ((!cullByFlag)||
        (!(ai1->flags&(cAtomFlag_exclude|cAtomFlag_ignore)))))
		{
		  visFlag=true;
		  break;
		}
  }
  if(!visFlag) {
    OOFreeP(I);
    return(NULL); /* skip if no dots are visible */
  }

  RepInit(&I->R);

  probe_radius = SettingGet_f(cs->Setting,obj->Obj.Setting,cSetting_solvent_radius);
  probe_radius2 = probe_radius*probe_radius;
  min_spacing = SettingGet_f(cs->Setting,obj->Obj.Setting,cSetting_min_mesh_spacing);

  I->N=NULL;
  I->NTot=0;
  I->V=NULL;
  I->VC=NULL;
  I->R.fRender=(void (*)(struct Rep *, CRay *, Pickable **))RepMeshRender;
  I->R.fFree=(void (*)(struct Rep *))RepMeshFree;
  I->Obj = (CObject*)(cs->Obj);
  I->R.fRecolor=(void (*)(struct Rep*, struct CoordSet*))RepMeshColor;
  I->R.fSameVis=(int (*)(struct Rep*, struct CoordSet*))RepMeshSameVis;
  I->LastVisib=NULL;
  I->LastColor=NULL;

  I->Radius = SettingGet_f(cs->Setting,obj->Obj.Setting,cSetting_mesh_radius);

  meshFlag=true;

  if(meshFlag) {

	 I->V=VLAMalloc(1000,sizeof(float),9,false);
	 I->N=VLAMalloc(100,sizeof(int),9,false);

	 I->N[0]=0;
	 
	  for(c=0;c<3;c++)
		{
		  /*		  minE[c]=(cs->Obj->Obj.extent[2*c]-MAX_VDW)-0.25;
					  maxE[c]=(cs->Obj->Obj.extent[2*c+1]+MAX_VDW)+0.25;
		  */
		  minE[c]=MAXFLOAT;
		  maxE[c]=-(MAXFLOAT);
		  } 

	  for(b=0;b<obj->NCSet;b++) {	 
		 ccs=obj->CSet[b];
		 if(ccs) {
			for(c=0;c<ccs->NIndex;c++) {
           ai1 = obj->AtomInfo+ccs->IdxToAtm[c]; /* WLD fixed 011218 */
			  if(ai1->visRep[cRepMesh]&&
              (inclH||(!ai1->hydrogen))&&
              ((!cullByFlag)||
               (!(ai1->flags&(cAtomFlag_ignore|cAtomFlag_exclude))))) {
             for(d=0;d<3;d++) {
               if(minE[d]>ccs->Coord[(3*c)+d])
                 minE[d]=ccs->Coord[(3*c)+d];
               if(maxE[d]<ccs->Coord[(3*c)+d])
                 maxE[d]=ccs->Coord[(3*c)+d];
             }
           }
			}
		 }
	  }

	  for(c=0;c<3;c++)
		{
		  minE[c]-=(MAX_VDW+0.25);
		  maxE[c]+=(MAX_VDW+0.25);
		}

	 subtract3f(maxE,minE,sizeE);
	 
	 size=sizeE[0];
	 if(sizeE[1]>size) size=sizeE[1];
	 if(sizeE[2]>size) size=sizeE[2];
	 
	 gridSize = size/80.0; /* grid size is the largest axis divided by 25 */

	 if(gridSize<min_spacing)
		gridSize=min_spacing;

	 for(c=0;c<3;c++)
		dims[c] = ((sizeE[c]/gridSize)+1.5);
	 
	 field = IsosurfFieldAlloc(dims);
	 
	 for(a=0;a<dims[0];a++)
		for(b=0;b<dims[1];b++)
		  for(c=0;c<dims[2];c++)
			 F3(field->data,a,b,c) = 2.0;

	 OrthoBusyFast(0,1);
	 RepMeshGetSolventDots(I,cs,minE,maxE,probe_radius);
	 smap=MapNew(probe_radius,I->Dot,I->NDot,NULL);
	 map=MapNew(MAX_VDW+probe_radius,cs->Coord,cs->NIndex,NULL);
	 if(map&&smap)
		{
		  MapSetupExpress(smap);
		  MapSetupExpress(map);
		  for(a=0;a<dims[0];a++)
			 {
				OrthoBusyFast(dims[0]+a,dims[0]*3);
				point[0]=minE[0]+a*gridSize;
				for(b=0;b<dims[1];b++)
				  {
					 point[1]=minE[1]+b*gridSize;
					 for(c=0;c<dims[2];c++)
						{
						  point[2]=minE[2]+c*gridSize;
                    copy3f(point,F3Ptr(field->points,a,b,c));
						  aNear = -1;
						  aLen = MAXFLOAT;
						  MapLocus(map,point,&h,&k,&l);						  
						  d=*(MapEStart(map,h,k,l));
						  if(d) {
							 cur=map->EList[d++];
							 while(cur>=0) {
                        ai1 = obj->AtomInfo+cs->IdxToAtm[cur];
                        if((inclH||(!ai1->hydrogen))&&
                           ((!cullByFlag)||
                            (!(ai1->flags&cAtomFlag_ignore)))) {
                          vLen=diffsq3f(point,cs->Coord+(cur*3));
                          if(vLen<aLen)
                            {
                              aLen=vLen;
                              aNear=cur;
                            }
                        }
								cur=map->EList[d++];
							 }
						  }						
						  if(aNear>=0)
							 {
								pVal = sqrt1f(aLen); /* pVal is the distance from atom center */
								vdw = cs->Obj->AtomInfo[cs->IdxToAtm[aNear]].vdw;
								if((pVal>=vdw)&&(pVal<(vdw+(probe_radius*1.6)))) {
								  escFlag=true;
								  /* this point lies within a water radius of the atom, so
									  lets see if it is actually near a water*/
								  aLen = MAXFLOAT;

								  MapLocus(smap,point,&h,&k,&l);						  
								  d=*(MapEStart(smap,h,k,l));
								  if(d) {
									 cur=smap->EList[d++];
									 while(cur>=0) {
                                
                              vLen=diffsq3f(point,I->Dot+(cur*3));
                              if(vLen<probe_radius2) {
                                escFlag=false;
                                break;
                              } else if(vLen<aLen) {
                                aLen=vLen;
                              }
                              
										cur=smap->EList[d++];
									 }
								  }
								  if(escFlag) {
									 if(pVal<(vdw+probe_radius)) 
										pVal=probe_radius/sqrt1f(aLen); /* yes it is - aLen is the distance*/
									 else
										pVal=0.0; /* out in bulk solvent */
								  } else {
									 pVal=pVal/vdw; 
									 if(pVal<1.0) /* no it is not near a water, so set above 1.0 to get rid of speckles*/
										pVal=1.1;
								  }
								} else { /* either too close, or too far from atom...*/
								  pVal=pVal/vdw;
								}
								if(pVal<F3(field->data,a,b,c))
								  F3(field->data,a,b,c) = pVal;
							 }
						}
				  }
			 }	 
		}
	 MapFree(smap);
	 MapFree(map);
	 FreeP(I->Dot);	 
	 OrthoBusyFast(2,3);
	 IsosurfVolume(field,1.0,&I->N,&I->V,NULL,0);
	 IsosurfFieldFree(field);
	 n=I->N;
	 I->NTot=0;
	 while(*n) I->NTot+=*(n++);
	 RepMeshColor(I,cs);
	 OrthoBusyFast(3,4);
  }
  OrthoBusyFast(4,4);
  return((void*)(struct Rep*)I);
}

void RepMeshGetSolventDots(RepMesh *I,CoordSet *cs,float *min,float *max,float probe_radius)
{
  ObjectMolecule *obj;
  int a,b,c,a1,a2,flag,i,h,k,l,j;
  float *v,*v0,vdw;
  MapType *map;
  int inFlag,*p,*dot_flag;
  SphereRec *sp = Sphere2;
  int cavity_cull;
  float probe_radius_plus;
  int dotCnt,maxCnt,maxDot=0;
  int cnt;
  int inclH,mesh_mode,cullByFlag;
  AtomInfoType *ai1,*ai2;
  obj = cs->Obj;

  cavity_cull = (int)SettingGet_f(cs->Setting,obj->Obj.Setting,cSetting_cavity_cull);

  mesh_mode = SettingGet_f(cs->Setting,obj->Obj.Setting,cSetting_mesh_mode);
  cullByFlag = (mesh_mode==cRepMesh_by_flags);
  inclH = !(mesh_mode==cRepMesh_heavy_atoms);

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

          ai1 = obj->AtomInfo+cs->IdxToAtm[a];
          if((inclH||(!ai1->hydrogen))&&
             ((!cullByFlag)||
              (!(ai1->flags&(cAtomFlag_ignore))))) {
            OrthoBusyFast(a,cs->NIndex*3);
            dotCnt=0;
            a1 = cs->IdxToAtm[a];
            v0 = cs->Coord+3*a;
            vdw = cs->Obj->AtomInfo[a1].vdw+probe_radius;
            inFlag=true;
            for(c=0;c<3;c++)
              {
                if((min[c]-v0[c])>vdw) { inFlag=false;break;};
                if((v0[c]-max[c])>vdw) { inFlag=false;break;};
              }
            if(inFlag)
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

                      ai2 = obj->AtomInfo + cs->IdxToAtm[j];
                      if((inclH||(!ai2->hydrogen))&&
                         ((!cullByFlag)||
                          (!(ai2->flags&cAtomFlag_ignore))))
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
  /* printf("NDot: %d->%d\n",c,I->NDot);*/
}


