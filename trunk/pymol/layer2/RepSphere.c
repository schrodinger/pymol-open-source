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
  int N,NC;
  int cullFlag,spheroidFlag;
  int *LastVisib;
  int *LastColor;
} RepSphere;

#include"ObjectMolecule.h"

void RepSphereRender(RepSphere *I,CRay *ray,Pickable **pick);
void RepSphereFree(RepSphere *I);
int RepSphereSameVis(RepSphere *I,CoordSet *cs);

void RepSphereInit(void)
{
}

void RepSphereFree(RepSphere *I)
{
  FreeP(I->VC);
  FreeP(I->V);
  FreeP(I->NT);
  FreeP(I->LastColor);
  FreeP(I->LastVisib);
  OOFreeP(I);
}

void RepSphereRender(RepSphere *I,CRay *ray,Pickable **pick)
{
  float *v=I->V,*vc;
  int c=I->N;
  int cc=0,*nt=NULL;
  int a;
  int flag;
  SphereRec *sp;
  float restart;

  if(ray) {
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
  } else if(pick&&PMGUI) {
  } else if(PMGUI) {
	 if(I->cullFlag) {
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
  ObjectMolecule *obj;
  int a,b,c,a1,c1,a2,i,j,k,h,l;
  float *v,*v0,*vc,vdw,v1[3],vv0[3],vv1[3],vv2[3];
  float restart;
  int *q, *s,q0,q1,q2;
  int *lv,*lc,*cc;
  SphereRec *sp = Sphere0; 
  int ds,*nt,flag;
  int *visFlag = NULL;
  MapType *map = NULL;
  int vFlag;
  AtomInfoType *ai2;
  int spheroidFlag = false;
  float spheroid_scale;
  float *sphLen,sphTmp,*sphNorm,*sphTmpN;
  float sphere_scale;
  float tn[3],vt1[3],vt2[3],xtn[3],*tn0,*tn1,*tn2;

  OOAlloc(RepSphere);

  obj = cs->Obj;
  vFlag=false;
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

  RepInit(&I->R);
  ds = (int)SettingGet_f(cs->Setting,obj->Obj.Setting,cSetting_sphere_quality);
  if(ds<0) ds=0;
  switch(ds) {
  case 0: sp=Sphere0; break;
  case 1: sp=Sphere1; break;
  case 2: sp=Sphere2; break;
  default: sp=Sphere3; break;
  }

  spheroid_scale=SettingGet_f(cs->Setting,obj->Obj.Setting,cSetting_spheroid_scale);
  if(spheroid_scale&&cs->Spheroid) 
    spheroidFlag=1;
  else
    spheroidFlag=0;

  sphere_scale=SettingGet_f(cs->Setting,obj->Obj.Setting,cSetting_sphere_scale);

  I->R.fRender=(void (*)(struct Rep *, CRay *, Pickable **))RepSphereRender;
  I->R.fFree=(void (*)(struct Rep *))RepSphereFree;
  I->R.fSameVis=(int (*)(struct Rep*, struct CoordSet*))RepSphereSameVis;
  I->R.fRecolor=NULL;
  I->LastVisib=NULL;
  I->LastColor=NULL;

  /* raytracing primitives */
  

  I->VC=(float*)mmalloc(sizeof(float)*cs->NIndex*7);
  ErrChkPtr(I->VC);
  I->NC=0;
  
  I->NT=NULL;
  
  v=I->VC; 
  
  I->spheroidFlag=spheroidFlag;
  for(a=0;a<cs->NIndex;a++)
    {
      a1 = cs->IdxToAtm[a];
      if(obj->AtomInfo[a1].visRep[cRepSphere])
        {
          I->NC++;
          c1=*(cs->Color+a);
          vc = ColorGet(c1); /* save new color */
          *(v++)=*(vc++);
          *(v++)=*(vc++);
          *(v++)=*(vc++);
          v0 = cs->Coord+3*a;			 
          *(v++)=*(v0++);
          *(v++)=*(v0++);
          *(v++)=*(v0++);
          *(v++)=obj->AtomInfo[a1].vdw*sphere_scale;
        }
    }

  if(I->NC) 
	 I->VC=(float*)mrealloc(I->VC,sizeof(float)*(v-I->VC));
  else
	 I->VC=(float*)mrealloc(I->VC,1);


  if(SettingGet_f(cs->Setting,obj->Obj.Setting,cSetting_cull_spheres)&&(!spheroidFlag)) {
	 I->V=(float*)mmalloc(sizeof(float)*cs->NIndex*(sp->NVertTot*19));
	 ErrChkPtr(I->V);

	 I->NT=Alloc(int,cs->NIndex);
	 ErrChkPtr(I->NT);

	 I->cullFlag = true;
	 visFlag = Alloc(int,sp->nDot);
	 ErrChkPtr(visFlag);

	 if(I->cullFlag) {
		map=MapNew(MAX_VDW*sphere_scale,cs->Coord,cs->NIndex,NULL);
		if(map) MapSetupExpress(map);
	 }
  } else {
	 I->V=(float*)mmalloc(sizeof(float)*cs->NIndex*(3+sp->NVertTot*6));
	 ErrChkPtr(I->V);

	 I->cullFlag = false;
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
			 c1=*(cs->Color+a);
			 v0 = cs->Coord+3*a;
			 vdw = cs->Obj->AtomInfo[a1].vdw*sphere_scale;
			 vc = ColorGet(c1);
			 *(v++)=*(vc++);
			 *(v++)=*(vc++);
			 *(v++)=*(vc++);

			 if(I->cullFlag&&(!spheroidFlag)) {
				for(b=0;b<sp->nDot;b++) /* Sphere culling mode - more strips, but many fewer atoms */
				  {
					 v1[0]=v0[0]+vdw*sp->dot[b].v[0];
					 v1[1]=v0[1]+vdw*sp->dot[b].v[1];
					 v1[2]=v0[2]+vdw*sp->dot[b].v[2];
					 
					 MapLocus(map,v1,&h,&k,&l);
					 
					 visFlag[b]=1;
					 i=*(MapEStart(map,h,k,l));
					 if(i) {
						j=map->EList[i++];
						while(j>=0) {
						  a2 = cs->IdxToAtm[j];
						  if(obj->AtomInfo[a2].visRep[cRepSphere]) {
							 if(j!=a)
								if(within3f(cs->Coord+3*j,v1,cs->Obj->AtomInfo[a2].vdw*sphere_scale))
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

						  v1[0]=v0[0]+vdw*sp->dot[q0].v[0];
						  v1[1]=v0[1]+vdw*sp->dot[q0].v[1];
						  v1[2]=v0[2]+vdw*sp->dot[q0].v[2];

						  v1[0]+=v0[0]+vdw*sp->dot[q1].v[0];
						  v1[1]+=v0[1]+vdw*sp->dot[q1].v[1];
						  v1[2]+=v0[2]+vdw*sp->dot[q1].v[2];

						  v1[0]+=v0[0]+vdw*sp->dot[q2].v[0];
						  v1[1]+=v0[1]+vdw*sp->dot[q2].v[1];
						  v1[2]+=v0[2]+vdw*sp->dot[q2].v[2];

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
									 if(within3f(cs->Coord+3*j,v1,cs->Obj->AtomInfo[a2].vdw*sphere_scale))
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
                            *(v++)=sp->dot[q2].v[0]; /* normal */
                            *(v++)=sp->dot[q2].v[1];
                            *(v++)=sp->dot[q2].v[2];
                            *(v++)=v0[0]+vdw*sp->dot[q2].v[0]; /* point */
                            *(v++)=v0[1]+vdw*sp->dot[q2].v[1];
                            *(v++)=v0[2]+vdw*sp->dot[q2].v[2];
                            *(v++)=sp->dot[q1].v[0]; /* normal */
                            *(v++)=sp->dot[q1].v[1];
                            *(v++)=sp->dot[q1].v[2];
                            *(v++)=v0[0]+vdw*sp->dot[q1].v[0]; /* point */
                            *(v++)=v0[1]+vdw*sp->dot[q1].v[1];
                            *(v++)=v0[2]+vdw*sp->dot[q1].v[2];
                            *(v++)=sp->dot[q0].v[0]; /* normal */
                            *(v++)=sp->dot[q0].v[1];
                            *(v++)=sp->dot[q0].v[2];
                            *(v++)=v0[0]+vdw*sp->dot[q0].v[0]; /* point */
                            *(v++)=v0[1]+vdw*sp->dot[q0].v[1];
                            *(v++)=v0[2]+vdw*sp->dot[q0].v[2];
                          } else { /* continue strip */
                            *(v++)=sp->dot[q0].v[0]; /* normal */
                            *(v++)=sp->dot[q0].v[1];
                            *(v++)=sp->dot[q0].v[2];
                            *(v++)=v0[0]+vdw*sp->dot[q0].v[0]; /* point */
                            *(v++)=v0[1]+vdw*sp->dot[q0].v[1];
                            *(v++)=v0[2]+vdw*sp->dot[q0].v[2];
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
                      *(v++)=v0[0]+sphTmp*sp->dot[*q].v[0]; /* point */
                      *(v++)=v0[1]+sphTmp*sp->dot[*q].v[1];
                      *(v++)=v0[2]+sphTmp*sp->dot[*q].v[2];
                      q++;
                    }

                  s++;
                }
            } else {
              for(b=0;b<sp->NStrip;b++)
                {
                  for(c=0;c<(*s);c++)
                    {
                      *(v++)=sp->dot[*q].v[0]; /* normal */
                      *(v++)=sp->dot[*q].v[1];
                      *(v++)=sp->dot[*q].v[2];
                      *(v++)=v0[0]+vdw*sp->dot[*q].v[0]; /* point */
                      *(v++)=v0[1]+vdw*sp->dot[*q].v[1];
                      *(v++)=v0[2]+vdw*sp->dot[*q].v[2];
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
  for(a=0;a<cs->NIndex;a++)
    {
      *(lv++) = (ai2 + cs->IdxToAtm[a])->visRep[cRepSphere];
      *(lc++) = *(cc++);
    }

  /*TODO: NEED TO SHRINK POINTERS HERE*/

  if(I->N) 
	 {
		I->V=(float*)mrealloc(I->V,sizeof(float)*(v-I->V));
		if(I->NT) I->NT=Realloc(I->NT,int,nt-I->NT);
	 }
  else
	 {
		I->V=(float*)mrealloc(I->V,1);
		if(I->NT) I->NT=Realloc(I->NT,int,1);
	 }
  FreeP(visFlag);
  if(map)  MapFree(map);
  return((void*)(struct Rep*)I);
}


#ifdef _this_code_is_not_used
                          /* sum normal */
                          tn0=sp->dot[q2].v; /* normal */
                          tn1=sp->dot[q1].v; /* normal */
                          tn2=sp->dot[q0].v;

                          /* sum normal */
                          add3f(tn0,tn1,tn);
                          add3f(tn2,tn,tn);

                          /* compute point */

                          vv0[0]=v0[0]+vdw*sp->dot[q2].v[0];
                          vv0[1]=v0[1]+vdw*sp->dot[q2].v[1];
                          vv0[2]=v0[2]+vdw*sp->dot[q2].v[2];

                          vv1[0]=v0[0]+vdw*sp->dot[q1].v[0];
                          vv1[1]=v0[1]+vdw*sp->dot[q1].v[1];
                          vv1[2]=v0[2]+vdw*sp->dot[q1].v[2];

                          vv2[0]=v0[0]+vdw*sp->dot[q0].v[0];
                          vv2[1]=v0[1]+vdw*sp->dot[q0].v[1];
                          vv2[2]=v0[2]+vdw*sp->dot[q0].v[2];

                          /* compute right-hand vector */

                          subtract3f(vv1,vv0,vt1);
                          subtract3f(vv2,vv1,vt2);
                          cross_product3f(vt1,vt2,xtn);
                          
                          /* reorder triangle if necessary */
                          
                          if(dot_product3f(xtn,tn)<0.0) {
                            copy3f(tn1,v);
                            v+=3;
                            copy3f(vv1,v);
                            v+=3;
                            copy3f(tn0,v);
                            v+=3;
                            copy3f(vv0,v);
                            v+=3;
                            copy3f(tn2,v);
                            v+=3;
                            copy3f(vv2,v);
                            v+=3;
                          } else {
                            copy3f(tn0,v);
                            v+=3;
                            copy3f(vv0,v);
                            v+=3;
                            copy3f(tn1,v);
                            v+=3;
                            copy3f(vv1,v);
                            v+=3;
                            copy3f(tn2,v);
                            v+=3;
                            copy3f(vv2,v);
                            v+=3;
                          }
#endif
