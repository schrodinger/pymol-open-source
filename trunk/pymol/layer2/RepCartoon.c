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
#include"RepCartoon.h"
#include"Color.h"
#include"Setting.h"
#include"Word.h"
#include"Scene.h"
#include"main.h"
#include"Feedback.h"
#include"CGO.h"
#include"Extrude.h"

typedef struct RepCartoon {
  Rep R;
  CGO *ray,*std;
} RepCartoon;

#include"ObjectMolecule.h"

void RepCartoonRender(RepCartoon *I,CRay *ray,Pickable **pick);
void RepCartoonFree(RepCartoon *I);

void RepCartoonFree(RepCartoon *I)
{
  if(I->ray)
    CGOFree(I->ray);
  if(I->std)
    CGOFree(I->std);
  RepFree(&I->R);
  OOFreeP(I);
}

void RepCartoonRender(RepCartoon *I,CRay *ray,Pickable **pick)
{
  if(ray) {
    PRINTFD(FB_RepCartoon)
      " RepCartoonRender: rendering raytracable...\n"
      ENDFD;
    
    if(I->ray)  
      CGORenderRay(I->ray,ray);
    else if(I->std)
      CGORenderRay(I->std,ray);    
  } else if(pick&&PMGUI) {
  } else if(PMGUI) {
    
    PRINTFD(FB_RepCartoon)
      " RepCartoonRender: rendering GL...\n"
      ENDFD;
    if(I->std)
      CGORenderGL(I->std);
  }
}

Rep *RepCartoonNew(CoordSet *cs)
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
  float loop_radius,angle,ratio,dot;
  int visFlag;
  CExtrude *ex;
  int n_p;
  int loop_quality;

  OOAlloc(RepCartoon);

  obj = cs->Obj;
  visFlag=false;
  for(a=0;a<cs->NIndex;a++) {
	 if(obj->AtomInfo[cs->IdxToAtm[a]].visRep[cRepCartoon])
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
  power_a=SettingGet_f(cs->Setting,obj->Obj.Setting,cSetting_ribbon_power);
  power_b=SettingGet_f(cs->Setting,obj->Obj.Setting,cSetting_ribbon_power_b);

  sampling = SettingGet_f(cs->Setting,obj->Obj.Setting,cSetting_cartoon_sampling);
  if(sampling<1) sampling=1;
  loop_radius = SettingGet_f(cs->Setting,obj->Obj.Setting,cSetting_cartoon_loop_radius);
  if(loop_radius<0.01) loop_radius=0.01;
  loop_quality = SettingGet_f(cs->Setting,obj->Obj.Setting,cSetting_cartoon_loop_quality);
  if(loop_quality<3) loop_quality=3;

  I->R.fRender=(void (*)(struct Rep *, CRay *, Pickable **))RepCartoonRender;
  I->R.fFree=(void (*)(struct Rep *))RepCartoonFree;
  I->R.fRecolor=NULL;

  I->ray=NULL;
  I->std=NULL;

  /* find all of the CA points */

  at = Alloc(int,cs->NIndex);
  pv = Alloc(float,cs->NIndex*3);
  seg = Alloc(int,cs->NIndex);
  
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
		  if(obj->AtomInfo[a1].visRep[cRepCartoon])
			 if(!obj->AtomInfo[a1].hetatm)
				if(WordMatch("CA",obj->AtomInfo[a1].name,1)<0)
				  {
					 if(a2>=0) {
						if((abs(obj->AtomInfo[a1].resv-obj->AtomInfo[a2].resv)>1)||
							(obj->AtomInfo[a1].chain[0]!=obj->AtomInfo[a2].chain[0])||
							(!WordMatch(obj->AtomInfo[a1].segi,obj->AtomInfo[a2].segi,1)))
						  {
							 a2=-1;
						  }
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

  I->std = CGONew();

  n_p = 0;
  if(nAt) {
    ex = ExtrudeNew();
    ExtrudeCircle(ex,loop_quality,loop_radius);
    ExtrudeAllocPointsNormalsColors(ex,cs->NIndex*(3*sampling+3));
	 v1=pv; /* points */
	 v2=tv; /* tangents */
	 v3=dv; /* direction vector */
    d = dl;
	 s=seg;
	 atp=at;
	 for(a=0;a<(nAt-1);a++)
		{
        if((!a)||((*s)!=*(s+1))) { /* new segment */
          if(n_p) {
            ExtrudeTruncate(ex,n_p);
            ExtrudeComputeTangents(ex);
            ExtrudeBuildNormals1f(ex);
            ExtrudeCGOSurfaceTube(ex,I->std,1);
          }

          ExtrudeTruncate(ex,0);
          v = ex->p;
          vc = ex->c;
          n_p = 0;
        }
		  if(ex&&(*s==*(s+1))) /* working in the same segment... */
			 {
				c1=*(cs->Color+*atp);
				c2=*(cs->Color+*(atp+1));

            dot =  dot_product3f(v2,v2+3);
            angle = acos(dot);
            
            if(angle>0.001) {
              ratio=angle/sqrt((pow(1-cos(angle),2)+pow(sin(angle),2)));
            } else {
              ratio=1.0;
            }
            len_seg = ratio*(*d)/2.0;
				for(b=0;b<sampling;b++) /* needs optimization */
				  {

                if(n_p==0) {

                  /* provide starting point on first point in segment only... */

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

                  /* store colors */
                  
                  *(vc++)=*(v0++);
                  *(vc++)=*(v0++);
                  *(vc++)=*(v0++);

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


                  /* compute tangent vector at point 
                     
                   *(vn++)=f1*(v2[0]*f2)+
                   f0*(v2[3]*f3);
                   *(vn++)=f1*(v2[1]*f2)+
                   f0*(v2[4]*f3);
                   *(vn++)=f1*(v2[2]*f2)+
                   f0*(v2[5]*f3);                 
                   vn+=6;
                  */

                  n_p++;

                }

					 f0=((float)b+1)/sampling;
                
                if(f0<0.5) {
                  v0 = ColorGet(c1);
                } else {
                  v0 = ColorGet(c2);
                }
                
                /* store colors */
                
                *(vc++)=*(v0++);
                *(vc++)=*(v0++);
                *(vc++)=*(v0++);
                
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

                /* 
                 *(vn++)=f1*(v2[0]*f2)+
                 f0*(v2[3]*f3);
                 *(vn++)=f1*(v2[1]*f2)+
                 f0*(v2[4]*f3);
                 *(vn++)=f1*(v2[2]*f2)+
                 f0*(v2[5]*f3);                 
                 vn+=6;
                 
                 if(n_p==1) { 
                                 copy3f(vn-9,vn-18);
                                 }
                */

					 n_p++;

				  }
            /*
              if(sampling>1) {
              copy3f(vn-18,vn-9);
              }
            */
			 }
		  v1+=3;
		  v2+=3;
		  v3+=3;
        d++;
		  atp+=1;
		  s++;
      }

    if(n_p) {
      ExtrudeTruncate(ex,n_p);
      ExtrudeComputeTangents(ex);
      ExtrudeBuildNormals1f(ex);
      ExtrudeCGOSurfaceTube(ex,I->std,1);
    }
    ExtrudeFree(ex);
	 FreeP(dv);
	 FreeP(dl);
	 FreeP(tv);
	 FreeP(nv);
  }

  FreeP(at);
  FreeP(seg);
  FreeP(pv);
  
  return((void*)(struct Rep*)I);
}


