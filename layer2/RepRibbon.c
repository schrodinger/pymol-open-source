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

#include"Base.h"
#include"OOMac.h"
#include"RepRibbon.h"
#include"Color.h"
#include"Setting.h"
#include"Word.h"
#include"Scene.h"
#include"main.h"
#include"Feedback.h"

typedef struct RepRibbon {
  Rep R;
  float *V;
  float *VC;
  float linewidth;
  int N,NC;
  int NS;
  int NP;
} RepRibbon;

#include"ObjectMolecule.h"

void RepRibbonRender(RepRibbon *I,CRay *ray,Pickable **pick);
void RepRibbonFree(RepRibbon *I);

void RepRibbonInit(void)
{
}

void RepRibbonFree(RepRibbon *I)
{
  FreeP(I->VC);
  FreeP(I->V);
  RepFree(&I->R);
  OOFreeP(I);
}

void RepRibbonRender(RepRibbon *I,CRay *ray,Pickable **pick)
{
  float *v=I->V;
  int c=I->N;
  Pickable *p;
  int i,j,ip;
  int last;

  if(ray) {
    PRINTFD(FB_RepRibbon)
      " RepRibbonRender: rendering raytracable...\n"
      ENDFD;

	 v=I->VC;
	 c=I->NC-1;
	 if(c>0)
		while(c--) {
		  ray->fCylinder3fv(ray,v+4,v+7,*(v+3),v,v);
		  v+=10;
		}
  } else if(pick&&PMGUI) {

    PRINTFD(FB_RepRibbon)
      " RepRibbonRender: rendering pickable...\n"
      ENDFD;

	 if(c) {
      i=(*pick)->index;
      p=I->R.P;
      last=-1;
      glBegin(GL_LINES);
      while(c--)
        {
          ip=*(v++);
          if(ip!=last) {
            i++;
            last=ip;
            if(!(*pick)[0].ptr) {
              /* pass 1 - low order bits */
              
              glColor3ub((i&0xF)<<4,(i&0xF0)|0x8,(i&0xF00)>>4); /* we're encoding the index into the color */
              VLACheck((*pick),Pickable,i);
              (*pick)[i] = p[ip]; /* copy object and atom info */
            } else { 
              /* pass 2 - high order bits */
              j=i>>12;
              glColor3ub((j&0xF)<<4,(j&0xF0)|0x8,(j&0xF00)>>4); 
            }
          }	 
          v+=3;
          glVertex3fv(v);
          v+=3;
          glVertex3fv(v);
          v+=3;
        }
      glEnd();
      (*pick)[0].index = i; /* pass the count */
    }
  } else if(PMGUI) {
    

    PRINTFD(FB_RepRibbon)
      " RepRibbonRender: rendering GL...\n"
      ENDFD;

    glLineWidth(I->linewidth);
	 
    if(c) {
	   glDisable(GL_LIGHTING);
      glBegin(GL_LINES);
      SceneResetNormal(true);
	 while(c--)
		{
        v++;
		  glColor3fv(v);
		  v+=3;
		  glVertex3fv(v);
		  v+=3;
		  glVertex3fv(v);
		  v+=3;
		}
	 glEnd();
    glEnable(GL_LIGHTING);
	 }
  }
}
Rep *RepRibbonNew(CoordSet *cs)
{
  ObjectMolecule *obj;
  int a,b,a1,a2,c1,c2,*i,*s,*at,*seg,nAt,*atp;
  float *v,*v0,*v1,*v2,*v3;
  float *pv=NULL;
  float *dv=NULL;
  float *nv=NULL;
  float *tv=NULL;
  float *vc=NULL;
  float f0,f1,f2,f3,len_seg;
  float *d;
  float *dl=NULL;
  int nSeg;
  int sampling;
  float  power_a = 5;
  float power_b = 5;
  float radius,angle,ratio,dot;
  int visFlag;
  Pickable *rp;
  OOAlloc(RepRibbon);

  obj = cs->Obj;
  visFlag=false;
  for(a=0;a<cs->NIndex;a++) {
	 if(obj->AtomInfo[cs->IdxToAtm[a]].visRep[cRepRibbon])
		{
		  visFlag=true;
		  break;
		}
  }
  if(!visFlag) {
    OOFreeP(I);
    return(NULL); /* skip if not visible */
  }

  RepInit(&I->R);
  power_a=SettingGet_f(cs->Setting,obj->Obj.Setting,cSetting_ribbon_power);
  power_b=SettingGet_f(cs->Setting,obj->Obj.Setting,cSetting_ribbon_power_b);

  sampling=SettingGet_f(cs->Setting,obj->Obj.Setting,cSetting_ribbon_sampling);
  if(sampling<1) sampling=1;
  radius=SettingGet_f(cs->Setting,obj->Obj.Setting,cSetting_ribbon_radius);
  if(radius<0.01) radius=0.01;

  I->R.fRender=(void (*)(struct Rep *, CRay *, Pickable **))RepRibbonRender;
  I->R.fFree=(void (*)(struct Rep *))RepRibbonFree;
  I->R.fRecolor=NULL;

  I->linewidth = SettingGet_f(cs->Setting,obj->Obj.Setting,cSetting_ribbon_width);

  /* find all of the CA points */

  at = Alloc(int,cs->NAtIndex);
  pv = Alloc(float,cs->NAtIndex*3);
  seg = Alloc(int,cs->NAtIndex);
  
  i=at;
  v=pv;
  s=seg;

  nAt = 0;
  nSeg = 0;
  a2=-1;
  for(a1=0;a1<cs->NAtIndex;a1++)
	 {
      if(obj->DiscreteFlag) {
        if(cs==obj->DiscreteCSet[a1]) 
          a=obj->DiscreteAtmToIdx[a1];
        else 
          a=-1;
      } else 
        a=cs->AtmToIdx[a1];
		if(a>=0)
		  if(obj->AtomInfo[a1].visRep[cRepRibbon])
			 if(!obj->AtomInfo[a1].hetatm)
				if(WordMatch("CA",obj->AtomInfo[a1].name,1)<0)
				  {
					 if(a2>=0) {
                  /*						if((abs(obj->AtomInfo[a1].resv-obj->AtomInfo[a2].resv)>1)||
                                    (obj->AtomInfo[a1].chain[0]!=obj->AtomInfo[a2].chain[0])||
                                    (!WordMatch(obj->AtomInfo[a1].segi,obj->AtomInfo[a2].segi,1)))*/
                  if(!ObjectMoleculeCheckBondSep(obj,a1,a2,3)) /* CA->N->C->CA = 3 bonds */
                    a2=-1;
					 }
					 if(a2<=0)
						nSeg++;
					 *(s++) = nSeg;
					 nAt++;
					 *(i++)=a;
					 v1 = cs->Coord+3*a;		
					 *(v++)=*(v1++);
					 *(v++)=*(v1++);
					 *(v++)=*(v1++);
					 
					 a2=a1;
				  }
	 }
  if(nAt)
	 {
		/* compute differences and normals */

		s=seg;
		v=pv;
		
		dv = Alloc(float,nAt*3);
		nv = Alloc(float,nAt*3);
		dl = Alloc(float,nAt);
		v1=dv;
		v2=nv;
		d=dl;
		
		for(a=0;a<(nAt-1);a++)
		  {
			 if(*s==*(s+1))
				{
				  subtract3f(v+3,v,v1);
				  *d = length3f(v1);
				  scale3f(v1,1.0/(*d),v2);
				}
			 d++;
			 v+=3;
			 v1+=3;
			 v2+=3;
			 s++;
		  }
		
		/* compute tangents */
		
		s=seg;
		v=nv;
		
		tv = Alloc(float,nAt*3+6);
		v1=tv;
		
		*(v1++)=*(v++); /* first segment */
		*(v1++)=*(v++);
		*(v1++)=*(v++);
		s++;
		
		for(a=1;a<(nAt-1);a++)
		  {
			 if((*s==*(s-1))&&(*s==*(s+1)))
				{
				  add3f(v,(v-3),v1);
				  normalize3f(v1);			 
				}
			 else if(*s==*(s-1))
				{
				  *(v1)=*(v-3);  /* end a segment */
				  *(v1+1)=*(v-2); 
				  *(v1+2)=*(v-1); 
				}
			 else if(*s==*(s+1))
				{
				  *(v1)=*(v);   /* new segment */
				  *(v1+1)=*(v+1); 
				  *(v1+2)=*(v+2); 
				}
			 v+=3;
			 v1+=3;
			 s++;
		  }
		
		*(v1++)=*(v-3); /* last segment */
		*(v1++)=*(v-2);
		*(v1++)=*(v-1);
		
	 }

  /* okay, we now have enough info to generate smooth interpolations */

  if(nAt) {
    I->R.P=Alloc(Pickable,nAt+2);
    ErrChkPtr(I->R.P);
    I->R.P[0].index=nAt;
    rp = I->R.P + 1; /* skip first record! */
  }

  I->VC=(float*)mmalloc(sizeof(float)*cs->NIndex*10*sampling);
  ErrChkPtr(I->VC);
  I->NC=0;

  I->V=(float*)mmalloc(sizeof(float)*cs->NIndex*9*sampling);
  ErrChkPtr(I->V);

  I->N=0;
  I->NP=0;
  v=I->V;
  vc=I->VC;
  if(nAt) {
	 v1=pv; /* points */
	 v2=tv; /* tangents */
	 v3=dv; /* direction vector */
    d = dl;
	 s=seg;
	 atp=at;
    rp->ptr = (void*)obj;
    rp->index = cs->IdxToAtm[*atp];
    rp->bond = -1;
    I->NP++;
    rp++;
	 for(a=0;a<(nAt-1);a++)
		{
        rp->ptr = (void*)obj;
        rp->index = cs->IdxToAtm[*(atp+1)]; /* store pickable for n+2 */
        rp->bond = -1;
        rp++;
        I->NP++;
		  if(*s==*(s+1))
			 {
				c1=*(cs->Color+*atp);
				c2=*(cs->Color+*(atp+1));

            dot =  dot_product3f(v2,v2+3);
            angle = acos(dot);
            
            if(angle>0.001) {
              ratio=angle/sqrt1f((pow(1-cos(angle),2)+pow(sin(angle),2)));
            } else {
              ratio=1.0;
            }
            len_seg = ratio*(*d)/2.0;
				for(b=0;b<sampling;b++) /* needs optimization */
				  {

					 f0=((float)b)/sampling; /* fraction of completion */

                if(f0<0.5) /* bias sampling towards the center of the curve */
                  f0=0.5*pow(2*f0,power_a);
                else
                  f0=1.0-0.5*pow(2*(1-f0),power_a);

					 if(f0<0.5) {
						v0 = ColorGet(c1);
					 } else {
						v0 = ColorGet(c2);
					 }

                /* store index */
                if(f0<0.5)
                  *(v++)=I->NP-1;
                else
                  *(v++)=I->NP;
                
                /* store colors */

					 *(v++)=*(v0);
					 *(vc++)=*(v0++);
					 *(v++)=*(v0);
					 *(vc++)=*(v0++);
					 *(v++)=*(v0);
					 *(vc++)=*(v0);

                /* start of line/cylinder */

					 f1=1.0-f0;
                f2=pow(f0,power_b);
                f3=pow(f1,power_b);

                *(v++)=f1*(v1[0]+v2[0]*len_seg*f2)+
                  f0*(v1[3]-v2[3]*len_seg*f3);
                *(v++)=f1*(v1[1]+v2[1]*len_seg*f2)+
                  f0*(v1[4]-v2[4]*len_seg*f3);
                *(v++)=f1*(v1[2]+v2[2]*len_seg*f2)+
                  f0*(v1[5]-v2[5]*len_seg*f3);

					 *(vc++)=radius;
					 *(vc++)=*(v-3);
					 *(vc++)=*(v-2);
					 *(vc++)=*(v-1);

					 f0=((float)b+1)/sampling;

                if(f0<0.5) /* bias sampling towards the center of the curve */
                  f0=0.5*pow(2*f0,power_a);
                else
                  f0=1.0-0.5*pow(2*(1-f0),power_a);


                /* end of line/cylinder */

					 f1=1.0-f0;
                f2=pow(f0,power_b);
                f3=pow(f1,power_b);

                *(v++)=f1*(v1[0]+v2[0]*len_seg*f2)+
                  f0*(v1[3]-v2[3]*len_seg*f3);
                *(v++)=f1*(v1[1]+v2[1]*len_seg*f2)+
                  f0*(v1[4]-v2[4]*len_seg*f3);
                *(v++)=f1*(v1[2]+v2[2]*len_seg*f2)+
                  f0*(v1[5]-v2[5]*len_seg*f3);

					 *(vc++)=*(v-3);
					 *(vc++)=*(v-2);
					 *(vc++)=*(v-1);
					 I->NC++;
					 I->N++;
				  }


			 }
		  v1+=3;
		  v2+=3;
		  v3+=3;
        d++;
		  atp+=1;
		  s++;
		}


	 FreeP(dv);
	 FreeP(dl);
	 FreeP(tv);
	 FreeP(nv);
  }

  FreeP(at);
  FreeP(seg);
  FreeP(pv);
  
  if(I->N) 
	 I->V=(float*)mrealloc(I->V,sizeof(float)*(v-I->V));
  else
	 I->V=(float*)mrealloc(I->V,1);

  if(I->NC) 
	 I->VC=(float*)mrealloc(I->VC,sizeof(float)*(vc-I->VC));
  else
	 I->VC=(float*)mrealloc(I->VC,1);
  
  return((void*)(struct Rep*)I);
}




#ifdef _OLD_CODE

				for(b=0;b<sampling;b++) /* needs optimization */
				  {

					 f0=((float)b)/sampling;
					 if(f0<0.5) {
						v0 = ColorGet(c1);
					 } else {
						v0 = ColorGet(c2);
					 }

					 *(v++)=*(v0);
					 *(vc++)=*(v0++);
					 *(v++)=*(v0);
					 *(vc++)=*(v0++);
					 *(v++)=*(v0);
					 *(vc++)=*(v0);

					 f1=2*(f0-0.5);
					 if(f1<0.0) 
						f0=(-pow(-f1,1.0/(power_b))+1.0)/2.0;
					 else
						f0=(pow(f1,1.0/(power_b))+1.0)/2.0;
					 f1=1-f0;
					 f2=1-(f0*2);
					 f3=1-fabs(pow(f2,power_a));
					 f4=f0;
					 *(v++)=v1[0]+f4*v3[0]+f3*((v2[0]*f1)-(v2[3]*f0));
					 *(v++)=v1[1]+f4*v3[1]+f3*((v2[1]*f1)-(v2[4]*f0));
					 *(v++)=v1[2]+f4*v3[2]+f3*((v2[2]*f1)-(v2[5]*f0));
					 *(vc++)=radius;
					 *(vc++)=*(v-3);
					 *(vc++)=*(v-2);
					 *(vc++)=*(v-1);

					 f0=((float)b+1)/sampling;
					 f1=2*(f0-0.5);
					 if(f1<0.0) 
						f0=(-pow(-f1,1.0/(power_b))+1.0)/2.0;
					 else
						f0=(pow(f1,1.0/(power_b))+1.0)/2.0;
					 f1=1-f0;
					 f2=1-(f0*2);
					 f3=1-fabs(pow(f2,power_a));
					 f4=f0;

					 *(v++)=v1[0]+f4*v3[0]+f3*((v2[0]*f1)-(v2[3]*f0));
					 *(v++)=v1[1]+f4*v3[1]+f3*((v2[1]*f1)-(v2[4]*f0));
					 *(v++)=v1[2]+f4*v3[2]+f3*((v2[2]*f1)-(v2[5]*f0));
					 *(vc++)=*(v-3);
					 *(vc++)=*(v-2);
					 *(vc++)=*(v-1);
					 I->NC++;
					 I->N++;
				  }

#endif
