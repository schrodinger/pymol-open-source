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
#include"OOMac.h"
#include"RepSphere.h"
#include"Color.h"
#include"Sphere.h"
#include"Map.h"
#include"Setting.h"
#include"main.h"

#ifdef NT
#undef NT
#endif

typedef struct RepSphere {
  Rep R;
  float *V;
  float *VC;
  SphereRec *SP;
  int *NT;
  int N,NC,NP;
  int cullFlag,spheroidFlag;
  int *LastVisib;
  int *LastColor;
} RepSphere;

#include"ObjectMolecule.h"

void RepSphereRender(RepSphere *I,CRay *ray,Pickable **pick);
void RepSphereFree(RepSphere *I);
int RepSphereSameVis(RepSphere *I,CoordSet *cs);


void RepSphereFree(RepSphere *I)
{
  FreeP(I->VC);
  FreeP(I->V);
  FreeP(I->NT);
  FreeP(I->LastColor);
  FreeP(I->LastVisib);
  RepPurge(&I->R);
  OOFreeP(I);
}

void RepSphereRender(RepSphere *I,CRay *ray,Pickable **pick)
{
  PyMOLGlobals *G=I->R.G;
  float *v=I->V,*vc;
  int c=I->N;
  int cc=0,*nt=NULL;
  int a;
  int flag;
  SphereRec *sp;
  float restart;
  float alpha;
  alpha = SettingGet_f(G,I->R.cs->Setting,I->R.obj->Setting,cSetting_sphere_transparency);

  alpha=1.0F-alpha;
  if(fabs(alpha-1.0)<R_SMALL4)
    alpha=1.0F;
  if(ray) {
    ray->fTransparentf(ray,1.0F-alpha);
    if(I->spheroidFlag) {
		sp=I->SP;
      while(c--)
        {
          vc = v;
          v+=3;
          for(a=0;a<sp->NStrip;a++) {
            cc=sp->StripLen[a];
            while((cc--)>2) {
              ray->fTriangle3fv(ray,v+3,v+9,v+15,v,v+6,v+12,vc,vc,vc);
              v+=6;
            }
            v+=12;
          }
        }
    } else {
      v=I->VC;
      c=I->NC;
      while(c--) {
        ray->fColor3fv(ray,v);
        v+=3;
        ray->fSphere3fv(ray,v,*(v+3));
        v+=4;
      }
    }
    ray->fTransparentf(ray,0.0);
  } else if(pick&&G->HaveGUI) {
    int trans_pick_mode = SettingGet_i(G,I->R.cs->Setting,
                                       I->R.obj->Setting,
                                       cSetting_transparency_picking_mode);
    if(I->R.P&&((trans_pick_mode==1)||((trans_pick_mode==2)&&(alpha>0.9F)))) {
      int i,j;
      Pickable *p;
		sp=I->SP;      
      i=(*pick)->index;
      
      p=I->R.P;
      
      if(I->spheroidFlag) {
        while(c--)
          {
            
            i++;          
            if(!(*pick)[0].ptr) {
              /* pass 1 - low order bits *            */
              glColor3ub((uchar)((i&0xF)<<4),(uchar)((i&0xF0)|0x8),(uchar)((i&0xF00)>>4)); 
              VLACheck((*pick),Pickable,i);
              p++;
              (*pick)[i] = *p; /* copy object and atom info */
            } else { 
              /* pass 2 - high order bits */           
              j=i>>12;            
              glColor3ub((uchar)((j&0xF)<<4),(uchar)((j&0xF0)|0x8),(uchar)((j&0xF00)>>4));             
            }			 
            
            v+=3;
            for(a=0;a<sp->NStrip;a++) {
              cc=sp->StripLen[a];
              glBegin(GL_TRIANGLE_STRIP);
              while((cc--)>0) {
                glNormal3fv(v);
                glVertex3fv(v+3);
                v+=6;
              }
              glEnd();
            }
          }
      } else {
        v=I->VC;
        c=I->NC;
        while(c--) {
          i++;          
          if(!(*pick)[0].ptr) {
            /* pass 1 - low order bits *            */
            glColor3ub((uchar)((i&0xF)<<4),(uchar)((i&0xF0)|0x8),(uchar)((i&0xF00)>>4)); 
            VLACheck((*pick),Pickable,i);
            p++;
            (*pick)[i] = *p; /* copy object and atom info */
          } else { 
            /* pass 2 - high order bits */           
            j=i>>12;            
            glColor3ub((uchar)((j&0xF)<<4),(uchar)((j&0xF0)|0x8),(uchar)((j&0xF00)>>4));             
          }			 
          
          {
            int *s,*q,b;
            float *v0,vdw;

            v0 = v+3;
            vdw = v[6];
            q=sp->Sequence;
            s=sp->StripLen;
            for(b=0;b<sp->NStrip;b++)
              {
                glBegin(GL_TRIANGLE_STRIP);
                for(cc=0;cc<(*s);cc++)
                  {
                    glNormal3f(sp->dot[*q][0],
                               sp->dot[*q][1],
                               sp->dot[*q][2]);
                    glVertex3f(v0[0]+vdw*sp->dot[*q][0],
                               v0[1]+vdw*sp->dot[*q][1],
                               v0[2]+vdw*sp->dot[*q][2]);
                    q++;
                  }
                glEnd();
                s++;
              }
            v+=7;
          }
        }
      }
	 (*pick)[0].index = i;
    }
  } else if(G->HaveGUI) {
    int use_dlst;
    use_dlst = (int)SettingGet(G,cSetting_use_display_lists);

    if(use_dlst&&I->R.displayList) {
      glCallList(I->R.displayList);
    } else { /* display list */

      if(use_dlst) {
        if(!I->R.displayList) {
          I->R.displayList = glGenLists(1);
          if(I->R.displayList) {
            glNewList(I->R.displayList,GL_COMPILE_AND_EXECUTE);
          }
        }
      }

      if(I->cullFlag) {
      
      if(alpha==1.0) {
        
        nt=I->NT; /* number of passes for each sphere */
        while(c--) /* iterate through all atoms */
          {
            glColor3fv(v);
            v+=3;
            cc=*(nt++);
            flag=0;
            glBegin(GL_TRIANGLE_STRIP);
            while(cc--) { /* execute loop this many times */
              restart=*(v++);
              if(restart) {
                if(flag) {
                  glEnd();
                  glBegin(GL_TRIANGLE_STRIP);
                }
                if(restart==2.0) { /* swap triangle polarity */
                  glNormal3fv(v);
                  glVertex3fv(v+3);
                }
                glNormal3fv(v);
                v+=3;
                glVertex3fv(v);
                v+=3;
                glNormal3fv(v);
                v+=3;
                glVertex3fv(v);
                v+=3;
              }
              glNormal3fv(v);
              v+=3;
              glVertex3fv(v);
              v+=3;
              flag=1;
            }
            glEnd();
          }
      } else {
    
        nt=I->NT; /* number of passes for each sphere */
        while(c--) /* iterate through all atoms */
          {
            glColor4f(v[0],v[1],v[2],alpha);
            v+=3;
            cc=*(nt++);
            flag=0;
            glBegin(GL_TRIANGLE_STRIP);
            while(cc--) { /* execute loop this many times */
              restart=*(v++);
              if(restart) {
                if(flag) {
                  glEnd();
                  glBegin(GL_TRIANGLE_STRIP);
                }
                if(restart==2.0) { /* swap triangle polarity */
                  glNormal3fv(v);
                  glVertex3fv(v+3);
                }
                glNormal3fv(v);
                v+=3;
                glVertex3fv(v);
                v+=3;
                glNormal3fv(v);
                v+=3;
                glVertex3fv(v);
                v+=3;
              }
              glNormal3fv(v);
              v+=3;
              glVertex3fv(v);
              v+=3;
              flag=1;
            }
            glEnd();
          }
        
      }
	 } else {

      if(alpha==1.0) {
        
        sp=I->SP;
        while(c--)
          {
            glColor3fv(v);
            v+=3;
            for(a=0;a<sp->NStrip;a++) {
              glBegin(GL_TRIANGLE_STRIP);
              cc=sp->StripLen[a];
              while(cc--) {
                glNormal3fv(v);
                v+=3;
                glVertex3fv(v);
                v+=3;
              }
              glEnd();
            }
          }
      } else {
        sp=I->SP;
        while(c--)
          {
            glColor4f(v[0],v[1],v[2],alpha);
            v+=3;
            for(a=0;a<sp->NStrip;a++) {
              glBegin(GL_TRIANGLE_STRIP);
              cc=sp->StripLen[a];
              while(cc--) {
                glNormal3fv(v);
                v+=3;
                glVertex3fv(v);
                v+=3;
              }
              glEnd();
            }
          }
      }
    }
      if(use_dlst&&I->R.displayList) {
        glEndList();
      }
    }
  }
}

int RepSphereSameVis(RepSphere *I,CoordSet *cs)
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
      if(*(lv++)!=(ai + cs->IdxToAtm[a])->visRep[cRepSphere] ) {
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



Rep *RepSphereNew(CoordSet *cs)
{
  PyMOLGlobals *G=cs->G;
  ObjectMolecule *obj;
  int a,b,c,a1,c1,a2,i,j,k,h,l;
  float *v,*v0,*vc,vdw,v1[3];
  float restart;
  int *q, *s,q0,q1,q2;
  int *lv,*lc,*cc;
  SphereRec *sp = G->Sphere->Sphere[0];
  int ds,*nt,flag;
  int *visFlag = NULL;
  MapType *map = NULL;
  int vFlag;
  AtomInfoType *ai2;
  int spheroidFlag = false;
  float spheroid_scale;
  float *sphLen,sphTmp,*sphNorm,*sphTmpN;
  float sphere_scale,sphere_add=0.0;
  int one_color;
  int *map_flag=NULL,*mf;

#ifdef _this_code_is_not_used
  float vv0[3],vv1[3],vv2[3];
  float tn[3],vt1[3],vt2[3],xtn[3],*tn0,*tn1,*tn2;
#endif

  OOAlloc(G,RepSphere);

  
  obj = cs->Obj;
  vFlag=false;
  if(obj->RepVisCache[cRepSphere])
    for(a=0;a<cs->NIndex;a++) {
      if(obj->AtomInfo[cs->IdxToAtm[a]].visRep[cRepSphere])
        {
          vFlag=true;
          break;
        }
    }
  if(!vFlag) {
    OOFreeP(I);
    return(NULL); /* skip if no dots are visible */
  }

  RepInit(G,&I->R);
  ds = SettingGet_i(G,cs->Setting,obj->Obj.Setting,cSetting_sphere_quality);
  if(ds<0) ds=0;
  if(ds>4) ds=4;
  sp = G->Sphere->Sphere[ds];

  one_color=SettingGet_color(G,cs->Setting,obj->Obj.Setting,cSetting_sphere_color);

  spheroid_scale=SettingGet_f(G,cs->Setting,obj->Obj.Setting,cSetting_spheroid_scale);
  if(spheroid_scale&&cs->Spheroid) 
    spheroidFlag=1;
  else
    spheroidFlag=0;

  sphere_scale=SettingGet_f(G,cs->Setting,obj->Obj.Setting,cSetting_sphere_scale);

  I->R.fRender=(void (*)(struct Rep *, CRay *, Pickable **))RepSphereRender;
  I->R.fFree=(void (*)(struct Rep *))RepSphereFree;
  I->R.fSameVis=(int (*)(struct Rep*, struct CoordSet*))RepSphereSameVis;
  I->R.fRecolor=NULL;
  I->LastVisib=NULL;
  I->LastColor=NULL;
  I->R.obj=(CObject*)obj;
  I->R.cs = cs;
    I->NP = 0;

  /* raytracing primitives */
  

  I->VC=(float*)mmalloc(sizeof(float)*cs->NIndex*7);
  ErrChkPtr(G,I->VC);
  I->NC=0;
  map_flag = Calloc(int,cs->NIndex);

  I->NT=NULL;
  
  v=I->VC; 
  mf=map_flag;

  if(SettingGet_i(G,cs->Setting,obj->Obj.Setting,cSetting_sphere_solvent)) { /* are we generating a solvent surface? */
    sphere_add = SettingGet_f(G,cs->Setting,obj->Obj.Setting,cSetting_solvent_radius); /* if so, get solvent radius */
  }
  
  if(SettingGet_f(G,cs->Setting,obj->Obj.Setting,cSetting_pickable)) {
    I->R.P=Alloc(Pickable,cs->NIndex+1);
    ErrChkPtr(G,I->R.P);
  }

  I->spheroidFlag=spheroidFlag;
  for(a=0;a<cs->NIndex;a++)
    {
      a1 = cs->IdxToAtm[a];
      if(obj->AtomInfo[a1].visRep[cRepSphere])
        {

          if(I->R.P) {
            I->NP++;
            
            I->R.P[I->NP].ptr = (void*)obj;
            I->R.P[I->NP].index = a1;
            I->R.P[I->NP].bond = -1;
          }

          *mf=true;
          I->NC++;
          if(one_color==-1)
            c1=*(cs->Color+a);
          else
            c1=one_color;
          v0 = cs->Coord+3*a;			 
          if(ColorCheckRamped(G,c1)) {
            ColorGetRamped(G,c1,v0,v);
            v+=3;
          } else {
            vc = ColorGet(G,c1); /* save new color */
            *(v++)=*(vc++);
            *(v++)=*(vc++);
            *(v++)=*(vc++);
          }
          *(v++)=*(v0++);
          *(v++)=*(v0++);
          *(v++)=*(v0++);
          *(v++)=obj->AtomInfo[a1].vdw*sphere_scale+sphere_add;
        }
      mf++;
    }

  if(I->NC) 
	 I->VC=ReallocForSure(I->VC,float,(v-I->VC));
  else
	 I->VC=ReallocForSure(I->VC,float,1);
  if(I->R.P) {
    I->R.P = Realloc(I->R.P,Pickable,I->NP+1);
    I->R.P[0].index = I->NP;
  }


  I->cullFlag = (int)SettingGet_f(G,cs->Setting,obj->Obj.Setting,cSetting_cull_spheres);
  if(spheroidFlag) I->cullFlag=false;
  if((I->cullFlag<2)&&
     (SettingGet_f(G,cs->Setting,obj->Obj.Setting,cSetting_sculpting))) 
    /* optimize build-time performance when sculpting */
    I->cullFlag=false;
  if((I->cullFlag<2)&&
     (SettingGet(G,cSetting_roving_spheres)!=0.0F))
    I->cullFlag=false;
  if(I->cullFlag) {
	 I->V=(float*)mmalloc(sizeof(float)*I->NC*(sp->NVertTot*30)); /* double check 30 */
	 ErrChkPtr(G,I->V);

	 I->NT=Alloc(int,cs->NIndex);
	 ErrChkPtr(G,I->NT);

	 visFlag = Alloc(int,sp->nDot);
	 ErrChkPtr(G,visFlag);

    map=MapNewFlagged(G,MAX_VDW*sphere_scale+sphere_add,cs->Coord,cs->NIndex,NULL,map_flag);
    if(map) MapSetupExpress(map);
  } else {
	 I->V=(float*)mmalloc(sizeof(float)*I->NC*(3+sp->NVertTot*6));
	 ErrChkPtr(G,I->V);
  }

  /* rendering primitives */

  I->N=0;
  I->SP=sp;
  v=I->V;
  nt=I->NT;

  for(a=0;a<cs->NIndex;a++)
	 {
		a1 = cs->IdxToAtm[a];
		if(obj->AtomInfo[a1].visRep[cRepSphere])
		  {
          if(one_color==-1)
            c1=*(cs->Color+a);
          else
            c1=one_color;
			 v0 = cs->Coord+3*a;
			 vdw = cs->Obj->AtomInfo[a1].vdw*sphere_scale+sphere_add;
          if(ColorCheckRamped(G,c1)) {
            ColorGetRamped(G,c1,v0,v);
            v+=3;
          } else {
            vc = ColorGet(G,c1);
            *(v++)=*(vc++);
            *(v++)=*(vc++);
            *(v++)=*(vc++);
          }

			 if(I->cullFlag&&(!spheroidFlag)) {
				for(b=0;b<sp->nDot;b++) /* Sphere culling mode - more strips, but many fewer atoms */
				  {
					 v1[0]=v0[0]+vdw*sp->dot[b][0];
					 v1[1]=v0[1]+vdw*sp->dot[b][1];
					 v1[2]=v0[2]+vdw*sp->dot[b][2];
					 
					 MapLocus(map,v1,&h,&k,&l);
					 
					 visFlag[b]=1;
					 i=*(MapEStart(map,h,k,l));
					 if(i) {
						j=map->EList[i++];
						while(j>=0) {
						  a2 = cs->IdxToAtm[j];
						  if(obj->AtomInfo[a2].visRep[cRepSphere]) {
							 if(j!=a)
								if(within3f(cs->Coord+3*j,v1,
                                    cs->Obj->AtomInfo[a2].vdw*
                                    sphere_scale + sphere_add))
								  {
									 visFlag[b]=0;
									 break;
								  }
						  }
						  j=map->EList[i++];
						}
					 }
				  }
				q=sp->Sequence;
				s=sp->StripLen;
				for(b=0;b<sp->NStrip;b++) 
				  /* this is an attempt to fill in *some* of the cracks
					* by checking to see if the center of the triangle is visible 
					* IMHO - the increase in framerates is worth missing a triangle
					* here or there, and the user can always turn off sphere culling */
				  {
					 q+=2;
					 for(c=2;c<(*s);c++) {
						  q0=*q;
						  q1=*(q-1);
						  q2=*(q-2);

						  if((!visFlag[q0])&&(!visFlag[q1])&&(!visFlag[q2]))

						  v1[0]=v0[0]+vdw*sp->dot[q0][0];
						  v1[1]=v0[1]+vdw*sp->dot[q0][1];
						  v1[2]=v0[2]+vdw*sp->dot[q0][2];

						  v1[0]+=v0[0]+vdw*sp->dot[q1][0];
						  v1[1]+=v0[1]+vdw*sp->dot[q1][1];
						  v1[2]+=v0[2]+vdw*sp->dot[q1][2];

						  v1[0]+=v0[0]+vdw*sp->dot[q2][0];
						  v1[1]+=v0[1]+vdw*sp->dot[q2][1];
						  v1[2]+=v0[2]+vdw*sp->dot[q2][2];

						  v1[0]/=3;
						  v1[1]/=3;
						  v1[2]/=3;

						  flag=true;
						  i=*(MapEStart(map,h,k,l));
						  if(i) {
							 j=map->EList[i++];
							 while(j>=0) {
								a2 = cs->IdxToAtm[j];
								if(obj->AtomInfo[a2].visRep[cRepSphere]) {
								  if(j!=a)
									 if(within3f(cs->Coord+3*j,v1,cs->Obj->AtomInfo[a2].vdw*sphere_scale+sphere_add))
										{
										  flag=false;
										  break;
										}
								}
								j=map->EList[i++];
							 }
						  }
						  if(flag)
							 {
								visFlag[q0]=1;
								visFlag[q1]=1;
								visFlag[q2]=1;
							 }
						  q++;
					 }
					 s++;
				  }
				
				*(nt)=0; /* how many passes through the triangle renderer? */
				q=sp->Sequence;
				s=sp->StripLen;

				for(b=0;b<sp->NStrip;b++)
				  {
					 restart=1.0; /* startin a new strip */
					 for(c=0;c<(*s);c++)
						{
                    if(c>1) { /* on third vertex or better */
                      q0=*q; /* get the indices of the triangle in this strip */
                      q1=*(q-1);
                      q2=*(q-2);
                      if(visFlag[q0]||(visFlag[q1])||(visFlag[q2])) /* visible? */
                        {
                          *(v++) = restart; /* store continuing string flag */
                          
                          if(restart) { /* not continuing...this is a new strip */
                            if(c&0x1) /* make sure strip starts off "right" */
                              *(v-1)=2.0;
                            *(v++)=sp->dot[q2][0]; /* normal */
                            *(v++)=sp->dot[q2][1];
                            *(v++)=sp->dot[q2][2];
                            *(v++)=v0[0]+vdw*sp->dot[q2][0]; /* point */
                            *(v++)=v0[1]+vdw*sp->dot[q2][1];
                            *(v++)=v0[2]+vdw*sp->dot[q2][2];
                            *(v++)=sp->dot[q1][0]; /* normal */
                            *(v++)=sp->dot[q1][1];
                            *(v++)=sp->dot[q1][2];
                            *(v++)=v0[0]+vdw*sp->dot[q1][0]; /* point */
                            *(v++)=v0[1]+vdw*sp->dot[q1][1];
                            *(v++)=v0[2]+vdw*sp->dot[q1][2];
                            *(v++)=sp->dot[q0][0]; /* normal */
                            *(v++)=sp->dot[q0][1];
                            *(v++)=sp->dot[q0][2];
                            *(v++)=v0[0]+vdw*sp->dot[q0][0]; /* point */
                            *(v++)=v0[1]+vdw*sp->dot[q0][1];
                            *(v++)=v0[2]+vdw*sp->dot[q0][2];
                          } else { /* continue strip */
                            *(v++)=sp->dot[q0][0]; /* normal */
                            *(v++)=sp->dot[q0][1];
                            *(v++)=sp->dot[q0][2];
                            *(v++)=v0[0]+vdw*sp->dot[q0][0]; /* point */
                            *(v++)=v0[1]+vdw*sp->dot[q0][1];
                            *(v++)=v0[2]+vdw*sp->dot[q0][2];
                          }
                          restart=0.0;
                          (*nt)++;
                        } else {
                          restart = 1.0;/* next triangle is a new strip */
                        }
                    }
                    q++;
                  }
                s++;
              }
			 } else {
				q=sp->Sequence;
				s=sp->StripLen;
            if(spheroidFlag) {
              for(b=0;b<sp->NStrip;b++)
                {
                  sphLen = cs->Spheroid+(sp->nDot*a1);
                  sphNorm = cs->SpheroidNormal+(3*sp->nDot*a1);
                  for(c=0;c<(*s);c++)
                    {
                      sphTmpN = sphNorm + 3*(*q);
                      *(v++)=*(sphTmpN++);
                      *(v++)=*(sphTmpN++);
                      *(v++)=*(sphTmpN++);
                      sphTmp = (*(sphLen+(*q)))*spheroid_scale;
                      *(v++)=v0[0]+sphTmp*sp->dot[*q][0]; /* point */
                      *(v++)=v0[1]+sphTmp*sp->dot[*q][1];
                      *(v++)=v0[2]+sphTmp*sp->dot[*q][2];
                      q++;
                    }

                  s++;
                }
            } else {
              for(b=0;b<sp->NStrip;b++)
                {
                  for(c=0;c<(*s);c++)
                    {
                      *(v++)=sp->dot[*q][0]; /* normal */
                      *(v++)=sp->dot[*q][1];
                      *(v++)=sp->dot[*q][2];
                      *(v++)=v0[0]+vdw*sp->dot[*q][0]; /* point */
                      *(v++)=v0[1]+vdw*sp->dot[*q][1];
                      *(v++)=v0[2]+vdw*sp->dot[*q][2];
                      q++;
                    }
                  s++;
                }
            }
			 }
			 I->N++;
			 if(nt) nt++;
		  }
	 }

  if(!I->LastVisib) I->LastVisib = Alloc(int,cs->NIndex);
  if(!I->LastColor) I->LastColor = Alloc(int,cs->NIndex);
  lv = I->LastVisib;
  lc = I->LastColor;
  cc = cs->Color;
  obj=cs->Obj;
  ai2=obj->AtomInfo;
  if(one_color==-1) 
    for(a=0;a<cs->NIndex;a++)
      {
        *(lv++) = (ai2 + cs->IdxToAtm[a])->visRep[cRepSphere];
        *(lc++) = *(cc++);
      }
  else 
    for(a=0;a<cs->NIndex;a++)
      {
        *(lv++) = (ai2 + cs->IdxToAtm[a])->visRep[cRepSphere];
        *(lc++) = one_color;
      }


  if(I->N) 
	 {
		I->V=ReallocForSure(I->V,float,(v-I->V));
		if(I->NT) I->NT=ReallocForSure(I->NT,int,(nt-I->NT));
	 }
  else
	 {
		I->V=ReallocForSure(I->V,float,1);
		if(I->NT) I->NT=ReallocForSure(I->NT,int,1);
	 }


  FreeP(visFlag);
  FreeP(map_flag);
  if(map)  MapFree(map);
  return((void*)(struct Rep*)I);
}


