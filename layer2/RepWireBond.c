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
#include"os_gl.h"

#include"OOMac.h"
#include"RepWireBond.h"
#include"Color.h"
#include"Scene.h"
#include"main.h"
#include"Setting.h"

typedef struct RepWireBond {
  Rep R;
  float *V,*VP;
  /*  Pickable *P;*/
  int N,NP;
  float Width;
  float Radius;
} RepWireBond;

#include"ObjectMolecule.h"

void RepWireBondRender(RepWireBond *I,CRay *ray,Pickable **pick);
void RepWireBondFree(RepWireBond *I);
void RepValence(float *v,float *v1,float *v2,int *other,int a1,
                int a2,float *coord,float *color,int ord,float tube_size);

void RepWireBondFree(RepWireBond *I)
{
  FreeP(I->VP);
  FreeP(I->V);
  RepFree(&I->R);
  OOFreeP(I);
}

void RepWireBondRender(RepWireBond *I,CRay *ray,Pickable **pick)
{
  float *v=I->V;
  int c=I->N;
  unsigned int i,j;
  Pickable *p;

  if(ray) {

	 v=I->V;
	 c=I->N;
	 
	 while(c--) {
      /*      printf("%8.3f %8.3f %8.3f   %8.3f %8.3f %8.3f \n",v[3],v[4],v[5],v[6],v[7],v[8]);*/
      ray->fSausage3fv(ray,v+3,v+6,I->Radius,v,v);
		v+=9;
	 }

  } else if(pick&&PMGUI) {
	 
	 i=(*pick)->index;

	 v=I->VP;
	 c=I->NP;
	 p=I->R.P;

	 glBegin(GL_LINES);
	 
	 while(c--) {

		i++;

		if(!(*pick)[0].ptr) {
		  /* pass 1 - low order bits */

          glColor3ub((uchar)((i&0xF)<<4),(uchar)((i&0xF0)|0x8),(uchar)((i&0xF00)>>4)); /* we're encoding the index into the color */
		  VLACheck((*pick),Pickable,i);
		  p++;
		  (*pick)[i] = *p; /* copy object and atom info */
		} else { 
		  /* pass 2 - high order bits */

		  j=i>>12;

          glColor3ub((uchar)((j&0xF)<<4),(uchar)((j&0xF0)|0x8),(uchar)((j&0xF00)>>4)); 

		}			 

		glVertex3fv(v);
		v+=3;
		glVertex3fv(v);
		v+=3;

	 }
	 glEnd();
	 (*pick)[0].index = i; /* pass the count */
  } else if(PMGUI) {

    glLineWidth(I->Width);
	 
	 v=I->V;
	 c=I->N;
    
    glDisable(GL_LIGHTING); 
	 glBegin(GL_LINES);	 
	 SceneResetNormal(true);
	 while(c--) {
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

Rep *RepWireBondNew(CoordSet *cs)
{
  ObjectMolecule *obj;
  int a,a1,a2,c1,c2,s1,s2,b1,b2,ord,*o;
  BondType *b;
  int half_bonds,*other=NULL;
  float valence;
  float *v,*v0,*v1,*v2,h[3];
  int visFlag;
  int maxBond;
  float tmpColor[3];

  Pickable *rp;
  AtomInfoType *ai1,*ai2;
  OOAlloc(RepWireBond);
  obj = cs->Obj;

  PRINTFD(FB_RepWireBond)
    " RepWireBondNew-Debug: entered.\n"
    ENDFD;

  visFlag=false;
  b=obj->Bond;
  for(a=0;a<obj->NBond;a++)
    {
      b1 = b->index[0];
      b2 = b->index[1];
      b++;
      if(obj->AtomInfo[b1].visRep[cRepLine]||
         obj->AtomInfo[b2].visRep[cRepLine]) {
        visFlag=true;
        break;
      }
    }
  if(!visFlag) {
    OOFreeP(I);
    return(NULL); /* skip if no dots are visible */
  }

  maxBond = 0;

  b=obj->Bond;
  for(a=0;a<obj->NBond;a++)
    {
      b1 = b->index[0];
      b2 = b->index[1];
      b++;
      
      if(obj->DiscreteFlag) {
        if((cs==obj->DiscreteCSet[b1])&&(cs==obj->DiscreteCSet[b2])) {
          a1=obj->DiscreteAtmToIdx[b1];
          a2=obj->DiscreteAtmToIdx[b2];
        } else {
          a1=-1;
          a2=-1;
        }
      } else {
        a1=cs->AtmToIdx[b1];
        a2=cs->AtmToIdx[b2];
      }
      if((a1>=0)&&(a2>=0))
        {
          maxBond++;
        }
    }


  RepInit(&I->R);

  I->R.fRender=(void (*)(struct Rep *, CRay *, Pickable **))RepWireBondRender;
  I->R.fFree=(void (*)(struct Rep *))RepWireBondFree;
  I->Width = SettingGet_f(cs->Setting,obj->Obj.Setting,cSetting_line_width);
  I->Radius = SettingGet_f(cs->Setting,obj->Obj.Setting,cSetting_line_radius);

  half_bonds = SettingGet_f(cs->Setting,obj->Obj.Setting,cSetting_half_bonds);
  valence = SettingGet_f(cs->Setting,obj->Obj.Setting,cSetting_valence);

  if(valence!=0.0) /* build list of up to 2 connected atoms for each atom */
    {
      other=Alloc(int,2*obj->NAtom);
      o=other;
      for(a=0;a<obj->NAtom;a++) {
        *(o++)=-1;
        *(o++)=-1;
      }
      b=obj->Bond;
      for(a=0;a<obj->NBond;a++)
        {
          b1 = b->index[0];
          b2 = b->index[1];
          b++;
          if(obj->DiscreteFlag) {
            if((cs==obj->DiscreteCSet[b1])&&(cs==obj->DiscreteCSet[b2])) {
              a1=obj->DiscreteAtmToIdx[b1];
              a2=obj->DiscreteAtmToIdx[b2];
            } else {
              a1=-1;
              a2=-1;
            }
          } else {
            a1=cs->AtmToIdx[b1];
            a2=cs->AtmToIdx[b2];
          }
          if((a1>=0)&&(a2>=0))
            {
              o=other+2*a1;
              if(*o!=a2) {
                if(*o>=0) o++;
                *o=a2;
              }
              o=other+2*a2;              
              if(*o!=a1) {
                if(*o>=0) o++;
                *o=a1;
              }
            }
        }
    }

  I->N=0;
  I->NP=0;
  I->V=NULL;
  I->VP=NULL;
  I->R.P=NULL;
  I->R.fRecolor=NULL;

  if(obj->NBond) {
	 I->V=(float*)mmalloc(sizeof(float)*maxBond*54);
	 ErrChkPtr(I->V);
	 	 
	 v=I->V;
	 b=obj->Bond;
	 for(a=0;a<obj->NBond;a++)
		{
        b1 = b->index[0];
        b2 = b->index[1];
        ord = b->order;
        b++;
        /*
          b1 = *(b++);
          b2 = *(b++);
          ord = (*(b++));
        */
        if(obj->DiscreteFlag) {
          if((cs==obj->DiscreteCSet[b1])&&(cs==obj->DiscreteCSet[b2])) {
            a1=obj->DiscreteAtmToIdx[b1];
            a2=obj->DiscreteAtmToIdx[b2];
          } else {
            a1=-1;
            a2=-1;
          }
        } else {
          a1=cs->AtmToIdx[b1];
          a2=cs->AtmToIdx[b2];
        }
		  if((a1>=0)&&(a2>=0))
			 {
				s1=obj->AtomInfo[b1].visRep[cRepLine];
				s2=obj->AtomInfo[b2].visRep[cRepLine];

				if(!(s1&&s2))
              if(!half_bonds) {
                s1 = 0;
                s2 = 0;
              }

				if(s1||s2)
				  {	
					 c1=*(cs->Color+a1);
					 c2=*(cs->Color+a2);
					 
					 v1 = cs->Coord+3*a1;
					 v2 = cs->Coord+3*a2;
					 
					 if((c1==c2)&&s1&&s2&&(!ColorCheckRamped(c1))) {
						

						v0 = ColorGet(c1);

                  if((valence!=0.0)&&(ord>1)&&(ord<4)) {
                    RepValence(v,v1,v2,other,a1,a2,cs->Coord,v0,ord,valence);
                    v+=ord*9;
                    I->N+=ord;
                  } else {
                    I->N++;
                    *(v++)=*(v0++);
                    *(v++)=*(v0++);
                    *(v++)=*(v0++);
                    
                    *(v++)=*(v1++);
                    *(v++)=*(v1++);
                    *(v++)=*(v1++);
                    
                    *(v++)=*(v2++);
                    *(v++)=*(v2++);
                    *(v++)=*(v2++);
                  }
                } else {
						
						h[0]=(v1[0]+v2[0])/2;
						h[1]=(v1[1]+v2[1])/2;
						h[2]=(v1[2]+v2[2])/2;
						
						if(s1)
						  {
                      
                      if(ColorCheckRamped(c1)) {
                        ColorGetRamped(c1,v1,tmpColor);
                        v0=tmpColor;
                      } else {
                        v0 = ColorGet(c1);
                      }


                      if((valence!=0.0)&&(ord>1)&&(ord<4)) {
                        RepValence(v,v1,h,other,a1,a2,cs->Coord,v0,ord,valence);
                        v+=ord*9;
                        I->N+=ord;
							 } else {

							 I->N++;
							 *(v++)=*(v0++);
							 *(v++)=*(v0++);
							 *(v++)=*(v0++);
							 
							 *(v++)=*(v1++);
							 *(v++)=*(v1++);
							 *(v++)=*(v1++);
							 
							 *(v++)=h[0];
							 *(v++)=h[1];
							 *(v++)=h[2];
                      }
                    }
						if(s2)
						  {

                      if(ColorCheckRamped(c2)) {
                        ColorGetRamped(c2,v2,tmpColor);
                        v0 = tmpColor;
                      } else {
                        v0 = ColorGet(c2);
                      }
                      if((valence!=0.0)&&(ord>1)&&(ord<4)) {
                        RepValence(v,h,v2,other,a1,a2,cs->Coord,v0,ord,valence);
                        v+=ord*9;
                        I->N+=ord;
                      } else {
                        I->N++;
                        *(v++)=*(v0++);
                        *(v++)=*(v0++);
                        *(v++)=*(v0++);
                        
                        *(v++)=h[0];
                        *(v++)=h[1];
                        *(v++)=h[2];
                        
                        *(v++)=*(v2++);
                        *(v++)=*(v2++);
                        *(v++)=*(v2++);
                      }
                      
						  }
					 }
				  }
			 }
		}

	 I->V = Realloc(I->V,float,I->N*9);

	 /* now create pickable verson */

	 if(SettingGet_f(cs->Setting,obj->Obj.Setting,cSetting_pickable)) {
		I->VP=(float*)mmalloc(sizeof(float)*maxBond*6*2);
		ErrChkPtr(I->VP);
		
		I->R.P=Alloc(Pickable,2*maxBond+1);
		ErrChkPtr(I->R.P);
		rp = I->R.P + 1; /* skip first record! */

		v=I->VP;
		b=obj->Bond;
		for(a=0;a<obj->NBond;a++)
		  {
          b1 = b->index[0];
          b2 = b->index[1];
			 b++;
          if(obj->DiscreteFlag) {
            if((cs==obj->DiscreteCSet[b1])&&(cs==obj->DiscreteCSet[b2])) {
              a1=obj->DiscreteAtmToIdx[b1];
              a2=obj->DiscreteAtmToIdx[b2];
            } else {
              a1=-1;
              a2=-1;
            }
          } else {
            a1=cs->AtmToIdx[b1];
            a2=cs->AtmToIdx[b2];
          }
			 if((a1>=0)&&(a2>=0))
				{
              ai1=obj->AtomInfo+b1;
              ai2=obj->AtomInfo+b2;
				  s1=ai1->visRep[cRepLine];
				  s2=ai2->visRep[cRepLine];
				  
              if(!(s1&&s2)) {
                if(!half_bonds) {
                  s1 = 0;
                  s2 = 0;
                }
              } 

				  if(s1||s2)
					 {	
						v1 = cs->Coord+3*a1;
						v2 = cs->Coord+3*a2;
						
						h[0]=(v1[0]+v2[0])/2;
						h[1]=(v1[1]+v2[1])/2;
						h[2]=(v1[2]+v2[2])/2;
						
						if(s1&(!ai1->masked))
						  {
							 I->NP++;
                      rp->ptr = (void*)obj;
							 rp->index = b1;
                      rp->bond = a;
                      rp++;

							 *(v++)=*(v1++);
							 *(v++)=*(v1++);
							 *(v++)=*(v1++);
							 
							 *(v++)=h[0];
							 *(v++)=h[1];
							 *(v++)=h[2];
						  }
						if(s2&(!ai2->masked))
						  {
							 I->NP++;
                      rp->ptr = (void*)obj;
							 rp->index = b2;
                      rp->bond = a;
                      rp++;
							 							 
							 *(v++)=h[0];
							 *(v++)=h[1];
							 *(v++)=h[2];
							 
							 *(v++)=*(v2++);
							 *(v++)=*(v2++);
							 *(v++)=*(v2++);
						  }
					 }
				}
		  }
		I->R.P = Realloc(I->R.P,Pickable,I->NP+1);
		I->R.P[0].index = I->NP;
		I->VP = Realloc(I->VP,float,I->NP*9);
	 }
  }
  FreeP(other);
  return((void*)(struct Rep*)I);
}




void RepValence(float *v,float *v1,float *v2,int *other,
                int a1,int a2,float *coord,float *color,int ord,
                float tube_size)
{

  float d[3],t[3],p0[3],p1[3],p2[3],*vv;
  int a3,ck;

  v[0] = color[0];
  v[1] = color[1];
  v[2] = color[2];

  v[9] = color[0];
  v[10] = color[1];
  v[11] = color[2];

  /* direction vector */

  p0[0] = (v2[0] - v1[0]);
  p0[1] = (v2[1] - v1[1]);
  p0[2] = (v2[2] - v1[2]);
  
  copy3f(p0,d);
  normalize3f(p0);
  
  /* need a third atom to get planarity*/

  a3 = -1;
  ck = other[a1*2];
  if((ck>=0)&&(ck!=a2))
    a3=ck;
  else {
    ck = other[a1*2+1];
    if((ck>=0)&&(ck!=a2))
      a3=ck;
    else {
      ck = other[a2*2];
      if((ck>=0)&&(ck!=a1))
        a3=ck;
      else {
        ck = other[a2*2+1];
        if((ck>=0)&&(ck!=a1))
          a3=ck;
        }
      }
  }

  if(a3<0) {    
    t[0] = p0[0];
    t[1] = p0[1];
    t[2] = -p0[2];
  } else {
    vv= coord+3*a3;
    t[0] = *(vv++)-v1[0];
    t[1] = *(vv++)-v1[1];
    t[2] = *(vv++)-v1[2];
    normalize3f(t);
  }
  
  cross_product3f(d,t,p1);
  
  normalize3f(p1);

  if(length3f(p1)==0.0) {
    p1[0]=p0[1];
    p1[1]=p0[2];
    p1[2]=p0[0];
    cross_product3f(p0,p1,p2);
    normalize3f(p2);
  } else {
    cross_product3f(d,p1,p2);
    
    normalize3f(p2);
  }

  /* now we have a coordinate system*/
  
  t[0] = p2[0]*tube_size;
  t[1] = p2[1]*tube_size;
  t[2] = p2[2]*tube_size;

  switch(ord) {
  case 2:
    v[0] = color[0];
    v[1] = color[1];
    v[2] = color[2];

    v[3] = v1[0] - t[0];
    v[4] = v1[1] - t[1];
    v[5] = v1[2] - t[2];
    
    v[6] = v2[0] - t[0];
    v[7] = v2[1] - t[1];
    v[8] = v2[2] - t[2];

    v[9] = color[0];
    v[10] = color[1];
    v[11] = color[2];
    
    v[12] = v1[0] + t[0];
    v[13] = v1[1] + t[1];
    v[14] = v1[2] + t[2];
    
    v[15] = v2[0] + t[0];
    v[16] = v2[1] + t[1];
    v[17] = v2[2] + t[2];
    break;
  case 3:
    t[0]=t[0]*2;
    t[1]=t[1]*2;
    t[2]=t[2]*2;

    v[0] = color[0];
    v[1] = color[1];
    v[2] = color[2];

    v[3] = v1[0] - t[0];
    v[4] = v1[1] - t[1];
    v[5] = v1[2] - t[2];
    
    v[6] = v2[0] - t[0];
    v[7] = v2[1] - t[1];
    v[8] = v2[2] - t[2];

    v[9] = color[0];
    v[10] = color[1];
    v[11] = color[2];
    
    v[12] = v1[0] + t[0];
    v[13] = v1[1] + t[1];
    v[14] = v1[2] + t[2];
    
    v[15] = v2[0] + t[0];
    v[16] = v2[1] + t[1];
    v[17] = v2[2] + t[2];


    v[18] = color[0];
    v[19] = color[1];
    v[20] = color[2];
    
    v[21] = v1[0];
    v[22] = v1[1];
    v[23] = v1[2];
    
    v[24] = v2[0];
    v[25] = v2[1];
    v[26] = v2[2];

    break;
  }
}

