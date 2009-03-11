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
#include"PyMOLGlobals.h"
#include"Selector.h"

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
  int *LastVisib;
  int *LastColor;
  float max_vdw;
  int mesh_type;
} RepMesh;

#include"ObjectMolecule.h"

void RepMeshFree(RepMesh *I);
void RepMeshColor(RepMesh *I,CoordSet *cs);
int RepMeshSameVis(RepMesh *I,CoordSet *cs);

void RepMeshInit(void)
{
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

int RepMeshGetSolventDots(RepMesh *I,CoordSet *cs,float *min,float *max,float probe_radius);

static void RepMeshRender(RepMesh *I,RenderInfo *info)
{
  CRay *ray = info->ray;
  Picking **pick = info->pick;
  PyMOLGlobals *G=I->R.G;
  float *v=I->V;
  float *vc=I->VC;
  int *n=I->N;
  int c;
  float *col=NULL;
  float line_width = SceneGetDynamicLineWidth(info, I->Width);

  if(ray) {
	 if(n) {
      float radius;

      if(I->Radius<=0.0F) {
        radius = ray->PixelRadius*line_width/2.0F;
      } else {
        radius = I->Radius;
      }
      /* looks like were missing some code here --

        what about mesh_type?

      */

      if(I->oneColorFlag) 
        col=ColorGet(G,I->oneColor);
      ray->fColor3fv(ray,ColorGet(G,I->R.obj->Color));
      switch(I->mesh_type) {
      case 0:
        while(*n) {
          c=*(n++);
          if(c--) {
            vc+=3;
            v+=3;
            if(I->oneColorFlag) {
              while(c--) {
                ray->fSausage3fv(ray,v-3,v,radius,col,col);
                v+=3;
                vc+=3;
              }
            } else {
              while(c--) {
                ray->fSausage3fv(ray,v-3,v,radius,vc-3,vc);
                v+=3;
                vc+=3;
              }
            }
          }
        }
      case 1:
        while(*n) {
          c=*(n++);
          if(I->oneColorFlag) {
            ray->fColor3fv(ray,col);
            while(c--) {
              ray->fSphere3fv(ray,v,radius);
              v+=3;
              vc+=3;
            }
          } else {
            while(c--) {
              ray->fColor3fv(ray,vc);
              ray->fSphere3fv(ray,v,radius);
              v+=3;
              vc+=3;
            }
          }
        }
        break;
      }
     }
  } else if(G->HaveGUI && G->ValidContext) {
    if(pick) {
    } else {
      int use_dlst;
      int lighting = SettingGet_i(G,I->R.cs->Setting,I->R.obj->Setting,cSetting_mesh_lighting);

      SceneResetNormal(G,true);
      if(!lighting) {
        if(!info->line_lighting)
          glDisable(GL_LIGHTING); 
      }
      switch(I->mesh_type) {
      case 0:
        if(info->width_scale_flag) 
          glLineWidth(line_width*info->width_scale);
        else
          glLineWidth(line_width);
        break;
      case 1: 
        if(info->width_scale_flag) 
          glPointSize(SettingGet_f(G,I->R.cs->Setting,I->R.obj->Setting,cSetting_dot_width) *
                      info->width_scale);
        else
          glPointSize(SettingGet_f(G,I->R.cs->Setting,I->R.obj->Setting,cSetting_dot_width));
        break;
      }
      use_dlst = (int)SettingGet(G,cSetting_use_display_lists);

      if(use_dlst && I->R.displayList && I->R.displayListInvalid) {
        glDeleteLists(I->R.displayList,1);
        I->R.displayList = 0;
        I->R.displayListInvalid = false;
      }
      
      if(use_dlst&&I->R.displayList) {
        glCallList(I->R.displayList);
      } else { 
        SceneResetNormal(G,false);

        if(use_dlst) {
          if(!I->R.displayList) {
            I->R.displayList = glGenLists(1);
            if(I->R.displayList) {
              glNewList(I->R.displayList,GL_COMPILE_AND_EXECUTE);
            }
          }
        }
      
        switch(I->mesh_type) {
        case 0:
          if(n) {
            if(I->oneColorFlag) {
              while(*n)
                {
                  glColor3fv(ColorGet(G,I->oneColor));
                  c=*(n++);
                  glBegin(GL_LINE_STRIP);
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
          break;
        case 1:
          glPointSize(SettingGet_f(G,I->R.cs->Setting,I->R.obj->Setting,cSetting_dot_width));
          if(n) {
            if(I->oneColorFlag) {
              while(*n) {
                glColor3fv(ColorGet(G,I->oneColor));
                c=*(n++);
                glBegin(GL_POINTS);
                while(c--) {
                  glVertex3fv(v);
                  v+=3;
                }
                glEnd();
              }
            } else {
              while(*n) {
                c=*(n++);
                glBegin(GL_POINTS);
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
          break;
        }
        if(use_dlst&&I->R.displayList) {
          glEndList();
        }
      }
      if(!lighting)
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
  PyMOLGlobals *G=cs->State.G;
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
  int state = I->R.context.state;

  obj=cs->Obj;

  probe_radius = SettingGet_f(G,cs->Setting,obj->Obj.Setting,cSetting_solvent_radius);
  mesh_color = SettingGet_color(G,cs->Setting,obj->Obj.Setting,cSetting_mesh_color);
  mesh_mode = SettingGet_i(G,cs->Setting,obj->Obj.Setting,cSetting_mesh_mode);
  cullByFlag = (mesh_mode==cRepMesh_by_flags);
  inclH = !(mesh_mode==cRepMesh_heavy_atoms);
  
  if(!I->LastVisib) I->LastVisib = Alloc(int,cs->NIndex);
  if(!I->LastColor) I->LastColor = Alloc(int,cs->NIndex);
  lv = I->LastVisib;
  lc = I->LastColor;
  cc = cs->Color;
  ai2=obj->AtomInfo;
  for(a=0;a<cs->NIndex;a++) {
      *(lv++) = (ai2 + cs->IdxToAtm[a])->visRep[cRepMesh];
      *(lc++) = *(cc++);
  }

  if(I->mesh_type!=1) {
    I->Width = SettingGet_f(G,cs->Setting,obj->Obj.Setting,cSetting_mesh_width);
    I->Radius = SettingGet_f(G,cs->Setting,obj->Obj.Setting,cSetting_mesh_radius);
  } else {
    I->Width = SettingGet_f(G,cs->Setting,obj->Obj.Setting,cSetting_dot_width);
    I->Radius = SettingGet_f(G,cs->Setting,obj->Obj.Setting,cSetting_dot_radius);
  }
  I->R.displayListInvalid=true;

  if(I->NTot) {
	 obj=cs->Obj;
    if(ColorCheckRamped(G,mesh_color)) {
      I->oneColorFlag=false;
    } else {
      I->oneColorFlag=true;
    }
	 first_color=-1;
	 if(!I->VC) I->VC = Alloc(float,3*I->NTot);
	 vc=I->VC;
	 /* now, assign colors to each point */
	 map=MapNew(G,I->max_vdw+probe_radius,cs->Coord,cs->NIndex,NULL);
	 if(map)
		{
		  MapSetupExpress(map);
		  for(a=0;a<I->NTot;a++) {
            AtomInfoType *ai0 = NULL;
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
                    dist = (float)diff3f(v0,cs->Coord+j*3) - ai2->vdw;
                    if(dist<minDist)
                      {
                        i0=j;
                        ai0=ai2;
                        minDist=dist;
                      }
                  }
                j=map->EList[i++];
              }
            }
            
            if(i0>=0) {
              int at_mesh_color;

              AtomInfoGetSetting_color(G, ai0, cSetting_mesh_color, 
                                   mesh_color, &at_mesh_color);
              
              if(at_mesh_color!=-1) {
                c1 = at_mesh_color;
              } else {
                c1 = *(cs->Color+i0);
              }
              if(I->oneColorFlag) {
                if(first_color>=0) {
                  if(first_color!=c1)
                    I->oneColorFlag=false;
                } else first_color=c1;
              }
            }
            /*
            if(ColorCheckRamped(G,mesh_color)) {
              c1 = mesh_color;
              }
            */

            if(ColorCheckRamped(G,c1)) {
              I->oneColorFlag=false;
              ColorGetRamped(G,c1,v0,vc,state);
              vc+=3;
            } else {
              c0 = ColorGet(G,c1);
              *(vc++) = *(c0++);
              *(vc++) = *(c0++);
              *(vc++) = *(c0++);
            }
			 }
		  MapFree(map);
		}
    if(I->oneColorFlag) {
      I->oneColor=first_color;
    }
  } 
  /*
    if(mesh_color>=0) {
    I->oneColorFlag=1;
    I->oneColor=mesh_color;
  }
  */
  
}

Rep *RepMeshNew(CoordSet *cs,int state)
{
  PyMOLGlobals *G=cs->State.G;
  ObjectMolecule *obj;
  CoordSet *ccs;
  int a,b,c,d,h,k,l,*n;
  MapType *map = NULL,*smap=NULL;
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
  int inSolvFlag;
  float probe_radius,probe_radius2;
  float min_spacing;
  int visFlag;
  int mesh_mode;
  int cullByFlag;
  int inclH;
  int solv_acc;
  int mesh_type ;
  int mesh_skip;
  int ok=true;

  AtomInfoType *ai1;

  OOAlloc(G,RepMesh);

  PRINTFD(G,FB_RepMesh)
	 " RepMeshNew-DEBUG: entered with coord-set %p\n",(void*)cs
	 ENDFD;
  obj = cs->Obj;
  I->R.context.object = (void*)obj;
  I->R.context.state = state;

  probe_radius = SettingGet_f(G,cs->Setting,obj->Obj.Setting,cSetting_solvent_radius);
  probe_radius2 = probe_radius*probe_radius;
  solv_acc = (SettingGet_i(G,cs->Setting,obj->Obj.Setting,cSetting_mesh_solvent));
  mesh_type = SettingGet_i(G,cs->Setting,obj->Obj.Setting,cSetting_mesh_type);
  mesh_skip = SettingGet_i(G,cs->Setting,obj->Obj.Setting,cSetting_mesh_skip);
  
  mesh_mode = SettingGet_i(G,cs->Setting,obj->Obj.Setting,cSetting_mesh_mode);
  cullByFlag = (mesh_mode==cRepMesh_by_flags);
  inclH = !(mesh_mode==cRepMesh_heavy_atoms);
  visFlag=false;
  if(obj->RepVisCache[cRepMesh])
    for(a=0;a<cs->NIndex;a++) {
      ai1 = obj->AtomInfo+cs->IdxToAtm[a];
      if(ai1->visRep[cRepMesh]&&
         (inclH||(!ai1->hydrogen))&&
         ((!cullByFlag)||
          (!(ai1->flags&(cAtomFlag_exfoliate|cAtomFlag_ignore))))) {
        visFlag=true;
        break;
      }
    }
  if(!visFlag) {
    OOFreeP(I);
    return(NULL); /* skip if no dots are visible */
  }

  I->max_vdw = ObjectMoleculeGetMaxVDW(obj) + solv_acc*probe_radius;

  RepInit(G,&I->R);

  min_spacing = SettingGet_f(G,cs->Setting,obj->Obj.Setting,cSetting_min_mesh_spacing);

  I->N=NULL;
  I->NTot=0;
  I->V=NULL;
  I->VC=NULL;
  I->NDot=0;
  I->Dot=NULL;
  I->R.fRender=(void (*)(struct Rep *, RenderInfo *))RepMeshRender;
  I->R.fFree=(void (*)(struct Rep *))RepMeshFree;
  I->R.obj = (CObject*)cs->Obj;
  I->R.cs = cs;
  I->R.fRecolor=(void (*)(struct Rep*, struct CoordSet*))RepMeshColor;
  I->R.fSameVis=(int (*)(struct Rep*, struct CoordSet*))RepMeshSameVis;
  I->LastVisib=NULL;
  I->LastColor=NULL;
  I->mesh_type = mesh_type;
  I->Radius = SettingGet_f(G,cs->Setting,obj->Obj.Setting,cSetting_mesh_radius);

  meshFlag=true;

  if(meshFlag) {
    float trim_cutoff = SettingGet_f(G,cs->Setting,obj->Obj.Setting,cSetting_mesh_cutoff);
    int trim_flag = false;
    float *trim_vla = NULL;
    MapType *trim_map = NULL;

    int carve_state = 0;
    int carve_flag = false;
    float carve_cutoff= SettingGet_f(G,cs->Setting,obj->Obj.Setting,cSetting_mesh_carve_cutoff);
    char *carve_selection = NULL;
    float *carve_vla = NULL;
    MapType *carve_map = NULL;
    
    int clear_state = 0;
    int clear_flag = false;
    float clear_cutoff = SettingGet_f(G,cs->Setting,obj->Obj.Setting,cSetting_mesh_clear_cutoff);
    char *clear_selection = NULL;
    float *clear_vla = NULL;
    MapType *clear_map = NULL;

    int mesh_max = SettingGet_i(G,cs->Setting,obj->Obj.Setting,cSetting_mesh_grid_max);
    if(mesh_max<1)
      mesh_max = 1;

    if(trim_cutoff<0.0F) {
      trim_cutoff = I->max_vdw + probe_radius;
    }
      
    if(trim_cutoff>0.0F) {
      int nc=0;
      trim_vla = VLAlloc(float,cs->NIndex*3);
      for(c=0;c<cs->NIndex;c++) {
        ai1 = obj->AtomInfo+cs->IdxToAtm[c]; 
        if(ai1->visRep[cRepMesh]&&
           (inclH||(!ai1->hydrogen))&&
           ((!cullByFlag)||
            (!(ai1->flags&(cAtomFlag_ignore|cAtomFlag_exfoliate))))) {
          VLACheck(trim_vla,float,nc*3+2);
          {
			  float *src = cs->Coord + 3*c;
          float *dst = trim_vla + 3*nc;
          *(dst++) = *(src++);
          *(dst++) = *(src++);
          *(dst++) = *(src++); 
          nc++;
		  } 
        }
      }
      if(nc) {
        trim_map = MapNew(G, trim_cutoff, trim_vla, nc, NULL);
        if(trim_map) {
          MapSetupExpress(trim_map);
          trim_flag = true;
        }
      }
    }

    if(carve_cutoff>0.0F) {
      carve_state = SettingGet_i(G,cs->Setting,obj->Obj.Setting,cSetting_mesh_carve_state) - 1;
      carve_selection = SettingGet_s(G,cs->Setting,obj->Obj.Setting,cSetting_mesh_carve_selection);
      if(carve_selection) 
        carve_map = SelectorGetSpacialMapFromSeleCoord(G,
                                                       SelectorIndexByName(G,carve_selection),
                                                       carve_state,
                                                       carve_cutoff,&carve_vla);
      if(carve_map) {
        MapSetupExpress(carve_map);
        carve_flag = true;
      }
    }
    if(clear_cutoff>0.0F) {
      clear_state = SettingGet_i(G,cs->Setting,obj->Obj.Setting,cSetting_mesh_clear_state) - 1;
      clear_selection = SettingGet_s(G,cs->Setting,obj->Obj.Setting,cSetting_mesh_clear_selection);
      if(clear_selection) 
        clear_map = SelectorGetSpacialMapFromSeleCoord(G,
                                                       SelectorIndexByName(G,clear_selection),
                                                       clear_state,
                                                       clear_cutoff,&clear_vla);
      if(clear_map) {
        MapSetupExpress(clear_map);
        clear_flag = true;
      }
    }

    I->V=VLAMalloc(1000,sizeof(float),9,false);
    I->N=VLAMalloc(100,sizeof(int),9,false);
    
    I->N[0]=0;
	
    for(c=0;c<3;c++)	{
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
               (!(ai1->flags&(cAtomFlag_ignore|cAtomFlag_exfoliate))))) {
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

     for(c=0;c<3;c++)	{
       minE[c]-=(float)(I->max_vdw+0.25F);
       maxE[c]+=(float)(I->max_vdw+0.25F);
     }

	 subtract3f(maxE,minE,sizeE);
	 
	 size=sizeE[0];
	 if(sizeE[1]>size) size=sizeE[1];
	 if(sizeE[2]>size) size=sizeE[2];
	 
	 gridSize = size/mesh_max; /* grid size is the largest axis divided by 25 */

	 if(gridSize<min_spacing)
       gridSize=min_spacing;

	 for(c=0;c<3;c++)
       dims[c] = (int)((sizeE[c]/gridSize)+1.5F);
	 
	 field = IsosurfFieldAlloc(G,dims);
	 
	 for(a=0;a<dims[0];a++)
       for(b=0;b<dims[1];b++)
         for(c=0;c<dims[2];c++)
           F3(field->data,a,b,c) = 2.0;

	 OrthoBusyFast(G,0,1);
	 if(!solv_acc)
       ok = RepMeshGetSolventDots(I,cs,minE,maxE,probe_radius);
     if(ok) {
       smap=MapNew(G,probe_radius,I->Dot,I->NDot,NULL);
       map=MapNew(G,I->max_vdw+probe_radius,cs->Coord,cs->NIndex,NULL);
     }
	 if(ok&&map&&smap)	{
       MapSetupExpress(smap);
       MapSetupExpress(map);
       for(a=0;a<dims[0];a++)	 {
         OrthoBusyFast(G,dims[0]+a,dims[0]*3);
         point[0]=minE[0]+a*gridSize;
         for(b=0;b<dims[1];b++)  {
           point[1]=minE[1]+b*gridSize;
           for(c=0;c<dims[2];c++)	{
             float bestDist;
             float dist2vdw;
             float vdw_add = solv_acc*probe_radius;
             point[2]=minE[2]+c*gridSize;
             copy3f(point,F3Ptr(field->points,a,b,c));
             aNear = -1;
             bestDist = MAXFLOAT;
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
                   vLen=(float)diff3f(point,cs->Coord+(cur*3));
                   dist2vdw = vLen - (ai1->vdw + vdw_add);
                   if(dist2vdw<bestDist) {
                     bestDist = dist2vdw;
                     aLen=vLen;
                     aNear=cur;
                   }
                 }
                 cur=map->EList[d++];
               }
             }						
             if(aNear>=0) { /* near an atom... */
               pVal = aLen; /* pVal is the distance from atom center */
               vdw = cs->Obj->AtomInfo[cs->IdxToAtm[aNear]].vdw + solv_acc*probe_radius;
               if((!solv_acc)&&(pVal>vdw)&&(pVal<(vdw+probe_radius))) { /* point outside an atom */
                 inSolvFlag=false;
                 /* this point lies within a water radius of the atom, so
                    lets see if it is actually near a water*/
                 aLen = MAXFLOAT;

                 MapLocus(smap,point,&h,&k,&l);						  
                 d=*(MapEStart(smap,h,k,l));
                 if(d) {
                   cur=smap->EList[d++];
                   while(cur>=0) {
                     vLen= diffsq3f(point,I->Dot+(cur*3));
                     if(vLen<probe_radius2) { /* within a water radius */
                       inSolvFlag=true;
                       break;
                     } else if(vLen<aLen) { /* not within water radius */
                       aLen=vLen; /* but nearest distance anyway */
                     }
                     cur=smap->EList[d++];
                   }
                 }
                 aLen = (float)sqrt1f(aLen); 
                 if(inSolvFlag) { /* point is inside a water, so a linear crossing point works fine*/
                   pVal=pVal/vdw;
                 } else { /* point is not inside a water, so we need to guestimate a value depending
                           * on where the point is */
                   if(aLen<(2*probe_radius)) {
                     pVal=1.05F*(2*probe_radius-aLen)/probe_radius; 
                     /* within a radius, so assign based on water surface, but error on the side of exclusion */
                              
                   } else {
                     pVal = 0.0; /* further than one water away... */
                   }
                 }
               } else { /* either within atom or outside solvent shell */
                 pVal=pVal/vdw;
               }
               if(pVal<F3(field->data,a,b,c))
                 F3(field->data,a,b,c) = pVal;
             }
           }
         }
         if(G->Interrupt) {
           ok=false;
           break;
         }
       }
     }
	 MapFree(smap);
     MapFree(map);
	 FreeP(I->Dot);	 
	 OrthoBusyFast(G,2,3);
     if(ok) {
       IsosurfVolume(G,NULL,NULL,field,1.0,&I->N,&I->V,NULL,mesh_type,mesh_skip,1.0F);
     }
     IsosurfFieldFree(G,field);
     if(I->N && I->V && (carve_flag||clear_flag||trim_flag)) {
       int cur_size = VLAGetSize(I->N);
       if((mesh_type==0) && cur_size) {
         int *n = I->N;
         int *new_n = VLACalloc(int,cur_size);
         int new_size = 0;
         float *new_v = I->V;
         float *v = I->V;
         while( (c=*(n++)) ) {
           int new_c = 0;
           float *last_v = NULL;
           while(c--) {
             int a_keeper = false;
             if(trim_map) {
               register int i=*(MapLocusEStart(trim_map,v));
               if(i) {
                 register int j=trim_map->EList[i++];
                 while(j>=0) {
                   register float *v_targ = trim_vla+3*j;
                   if(within3f(v_targ,v,trim_cutoff)) {
                     a_keeper = true;
                     break;
                   }
                   j=trim_map->EList[i++];
                 }
               }
             } else {
               a_keeper = true;
             }
             
             if(a_keeper && carve_map) {
               register int i=*(MapLocusEStart(carve_map,v));
               a_keeper = false;
               if(i) {
                 register int j=carve_map->EList[i++];
                 while(j>=0) {
                   register float *v_targ = carve_vla+3*j;
                   if(within3f(v_targ,v,carve_cutoff)) {
                     a_keeper = true;
                     break;
                   }
                   j=carve_map->EList[i++];
                 }
               }
             }

             if(clear_map) { 
               register int i=*(MapLocusEStart(clear_map,v));
               if(i) {
                 register int j=clear_map->EList[i++];
                 while(j>=0) {
                   if(within3f(clear_vla+3*j,v,clear_cutoff)) {
                     a_keeper = false;
                     break;
                   }
                   j=clear_map->EList[i++];
                 }
               }
             }
             if(a_keeper) {
               if(new_c) {
                 new_c++;
                 *(new_v++) = *(v++); *(new_v++) = *(v++); *(new_v++) = *(v++);
               } else {
                 if(last_v) {
                   new_c+=2;
                   *(new_v++) = *(last_v++); *(new_v++) = *(last_v++); *(new_v++) = *(last_v++);
                   *(new_v++) = *(v++); *(new_v++) = *(v++); *(new_v++) = *(v++);
                   last_v = NULL;
                 } else {
                   last_v = v;
                   v+=3;
                 }
               }
             } else {
               last_v = NULL;
               v+=3;
               if(new_c) {
                 VLACheck(new_n,int,new_size+1); /* extends the zero count sentinel */
                 new_n[new_size] = new_c;
                 new_size++;
                 new_c = 0;
               }
             }
           }
           if(new_c) {
             VLACheck(new_n,int,new_size+1); /* extends the zero count sentinel */
             new_n[new_size] = new_c;
             new_size++;
             new_c = 0;
           }
         }
         VLAFreeP(I->N);
         I->N = new_n;
       }
     }
     MapFree(trim_map);
     MapFree(carve_map);
     MapFree(clear_map);
     VLAFreeP(trim_vla);
     VLAFreeP(carve_vla);
     VLAFreeP(clear_vla);
     n=I->N;
     I->NTot=0;
	 while(*n) I->NTot+=*(n++);
	 RepMeshColor(I,cs);
	 OrthoBusyFast(G,3,4);
  }
  OrthoBusyFast(G,4,4);
  return((void*)(struct Rep*)I);
}

int RepMeshGetSolventDots(RepMesh *I,CoordSet *cs,float *min,float *max,float probe_radius)
{
  PyMOLGlobals *G=cs->State.G;
  ObjectMolecule *obj=cs->Obj;
  int a,b,c,a1,a2,flag,i,h,k,l,j;
  int ok=true;
  float *v,*v0,vdw;
  MapType *map;
  int inFlag,*p,*dot_flag;
  SphereRec *sp = G->Sphere->Sphere[0];
  int cavity_cull;
  float probe_radius_plus;
  int dotCnt,maxCnt,maxDot=0;
  int cnt;
  int inclH,mesh_mode,cullByFlag;
  AtomInfoType *ai1,*ai2;
  int ds = SettingGet_i(G,cs->Setting,obj->Obj.Setting,cSetting_mesh_quality);

  if(ds<0) ds = 0;
  if(ds>4) ds = 4;
  sp = G->Sphere->Sphere[ds];

  cavity_cull = SettingGet_i(G,cs->Setting,obj->Obj.Setting,cSetting_cavity_cull);

  mesh_mode = SettingGet_i(G,cs->Setting,obj->Obj.Setting,cSetting_mesh_mode);
  cullByFlag = (mesh_mode==cRepMesh_by_flags);
  inclH = !(mesh_mode==cRepMesh_heavy_atoms);

  I->Dot=(float*)mmalloc(sizeof(float)*cs->NIndex*3*sp->nDot);
  ErrChkPtr(G,I->Dot);

  probe_radius_plus = probe_radius * 1.5F;

  I->NDot=0;
  map=MapNew(G,I->max_vdw+probe_radius,cs->Coord,cs->NIndex,NULL);
  if(map) {
    MapSetupExpress(map);
    maxCnt=0;
    v=I->Dot;
    for(a=0;a<cs->NIndex;a++) {

      ai1 = obj->AtomInfo+cs->IdxToAtm[a];
      if((inclH||(!ai1->hydrogen))&&
         ((!cullByFlag)||
          (!(ai1->flags&(cAtomFlag_ignore))))) {
        OrthoBusyFast(G,a,cs->NIndex*3);
        dotCnt=0;
        a1 = cs->IdxToAtm[a];
        v0 = cs->Coord+3*a;
        vdw = cs->Obj->AtomInfo[a1].vdw+probe_radius;
        inFlag=true;
        for(c=0;c<3;c++) {
          if((min[c]-v0[c])>vdw) { inFlag=false;break;};
          if((v0[c]-max[c])>vdw) { inFlag=false;break;};
        }
        if(inFlag)
          for(b=0;b<sp->nDot;b++)  {
            v[0]=v0[0]+vdw*sp->dot[b][0];
            v[1]=v0[1]+vdw*sp->dot[b][1];
            v[2]=v0[2]+vdw*sp->dot[b][2];
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
                  if(j!=a)  {
                    a2 = cs->IdxToAtm[j];
                    if(within3f(cs->Coord+3*j,v,cs->Obj->AtomInfo[a2].vdw+probe_radius)) {
                      flag=false;
                      break;
                    }
                  }
                j=map->EList[i++];
              }
            }
            if(flag) {
              dotCnt++;
              v+=3;
              I->NDot++;
            }
          }
        if(dotCnt>maxCnt) {
          maxCnt=dotCnt;
          maxDot=I->NDot-1;
        }
      }
      if(G->Interrupt) {
        ok=false;
        break;
      }
    }
    MapFree(map);
  }

  if(ok&&(cavity_cull>0)) {
    dot_flag=Alloc(int,I->NDot);
    ErrChkPtr(G,dot_flag);
    for(a=0;a<I->NDot;a++) {
      dot_flag[a]=0;
    }
    dot_flag[maxDot]=1; /* this guarantees that we have a valid dot */

    map=MapNew(G,probe_radius_plus,I->Dot,I->NDot,NULL);
    if(map) {
      MapSetupExpress(map);
		  
      flag=true;
      while(flag) {
        p=dot_flag;
        v=I->Dot;
		  
        flag=false;
        for(a=0;a<I->NDot;a++)	{
          if(!dot_flag[a]) {
            cnt=0;
            MapLocus(map,v,&h,&k,&l);
            i=*(MapEStart(map,h,k,l));
            if(i) {
              j=map->EList[i++];
              while(j>=0) {
                if(j!=a) {
                  if(within3f(I->Dot+(3*j),v,probe_radius_plus)) {
                    if(dot_flag[j]) {
                      *p=true;
                      flag=true;
                      break;
                    }
                    cnt++;
                    if(cnt>cavity_cull) {
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
        if(G->Interrupt) {
          ok=false;
          break;
        }
      }
      MapFree(map);
    }
    v0=I->Dot;
    v=I->Dot;
    p=dot_flag;
    c=I->NDot;
    I->NDot=0;
    for(a=0;a<c;a++)	{
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
  if(!ok) {
    FreeP(I->Dot);
    I->NDot = 0;
  }
  return ok;
}


