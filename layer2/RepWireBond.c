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
static void RepValence(float *v,float *v1,float *v2,int *other,int a1,
                int a2,float *coord,float *color,int ord,float tube_size);


static void RepAromatic(float *v1,float *v2,int *other,
                int a1,int a2,float *coord,float *color,
                float tube_size,int half_state,float **v_ptr,int *n_ptr)
{
  float d[3],t[3],p0[3],p1[3],p2[3],*vv;
  int a3;
  float *v = *v_ptr;
  float f,f_1;
  int n=*n_ptr;
  int double_sided;

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
  
  /* need a prioritized third atom to get planarity */

  a3 = ObjectMoleculeGetPrioritizedOther(other,a1,a2,&double_sided);

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

  switch(half_state) {
  case 0: /* full bond */

    t[0] = p2[0]*tube_size*2;
    t[1] = p2[1]*tube_size*2;
    t[2] = p2[2]*tube_size*2;
    
    v[0] = color[0];
    v[1] = color[1];
    v[2] = color[2];
    
    v[3] = v1[0];
    v[4] = v1[1];
    v[5] = v1[2];
    
    v[6] = v2[0];
    v[7] = v2[1];
    v[8] = v2[2];
    

    v[9] = color[0];
    v[10] = color[1];
    v[11] = color[2];
    
    f = 0.14F;
    f_1 = 1.0F-f;

    v[12] = (f_1*v1[0]+f*v2[0]) - t[0];
    v[13] = (f_1*v1[1]+f*v2[1]) - t[1];
    v[14] = (f_1*v1[2]+f*v2[2]) - t[2];

    f = 0.4F;
    f_1 = 1.0F-f;

    v[15] = (f_1*v1[0]+f*v2[0]) - t[0];
    v[16] = (f_1*v1[1]+f*v2[1]) - t[1];
    v[17] = (f_1*v1[2]+f*v2[2]) - t[2];

    v[18] = color[0];
    v[19] = color[1];
    v[20] = color[2];

    f = 0.6F;
    f_1 = 1.0F-f;

    v[21] = (f_1*v1[0]+f*v2[0]) - t[0];
    v[22] = (f_1*v1[1]+f*v2[1]) - t[1];
    v[23] = (f_1*v1[2]+f*v2[2]) - t[2];

    f = 0.84F;
    f_1 = 1.0F-f;

    v[24] = (f_1*v1[0]+f*v2[0]) - t[0];
    v[25] = (f_1*v1[1]+f*v2[1]) - t[1];
    v[26] = (f_1*v1[2]+f*v2[2]) - t[2];
    
    v+=27;
    n+=3;

    if(double_sided) {
      
      v[ 0] = color[0];
      v[ 1] = color[1];
      v[ 2] = color[2];
      
      f = 0.14F;
      f_1 = 1.0F-f;
      
      v[ 3] = (f_1*v1[0]+f*v2[0]) + t[0];
      v[ 4] = (f_1*v1[1]+f*v2[1]) + t[1];
      v[ 5] = (f_1*v1[2]+f*v2[2]) + t[2];
      
      f = 0.4F;
      f_1 = 1.0F-f;
      
      v[ 6] = (f_1*v1[0]+f*v2[0]) + t[0];
      v[ 7] = (f_1*v1[1]+f*v2[1]) + t[1];
      v[ 8] = (f_1*v1[2]+f*v2[2]) + t[2];
      
      v[ 9] = color[0];
      v[10] = color[1];
      v[11] = color[2];
      
      f = 0.6F;
      f_1 = 1.0F-f;
      
      v[12] = (f_1*v1[0]+f*v2[0]) + t[0];
      v[13] = (f_1*v1[1]+f*v2[1]) + t[1];
      v[14] = (f_1*v1[2]+f*v2[2]) + t[2];
      
      f = 0.84F;
      f_1 = 1.0F-f;
      
      v[15] = (f_1*v1[0]+f*v2[0]) + t[0];
      v[16] = (f_1*v1[1]+f*v2[1]) + t[1];
      v[17] = (f_1*v1[2]+f*v2[2]) + t[2];
      
      v+=18;
      n+=2;

    }

    break;
  case 1:

    t[0] = p2[0]*tube_size*2;
    t[1] = p2[1]*tube_size*2;
    t[2] = p2[2]*tube_size*2;
    
    v[0] = color[0];
    v[1] = color[1];
    v[2] = color[2];
    
    v[3] = v1[0];
    v[4] = v1[1];
    v[5] = v1[2];
    
    v[6] = (v2[0]+v1[0])/2.0F;
    v[7] = (v2[1]+v1[1])/2.0F;
    v[8] = (v2[2]+v1[2])/2.0F;
    
    v[9] = color[0];
    v[10] = color[1];
    v[11] = color[2];
    
    f = 0.14F;
    f_1 = 1.0F-f;

    v[12] = (f_1*v1[0]+f*v2[0]) - t[0];
    v[13] = (f_1*v1[1]+f*v2[1]) - t[1];
    v[14] = (f_1*v1[2]+f*v2[2]) - t[2];

    f = 0.4F;
    f_1 = 1.0F-f;

    v[15] = (f_1*v1[0]+f*v2[0]) - t[0];
    v[16] = (f_1*v1[1]+f*v2[1]) - t[1];
    v[17] = (f_1*v1[2]+f*v2[2]) - t[2];

    v+=18;
    n+=2;

    if(double_sided) {

      v[ 0] = color[0];
      v[ 1] = color[1];
      v[ 2] = color[2];
      
      f = 0.14F;
      f_1 = 1.0F-f;
      
      v[ 3] = (f_1*v1[0]+f*v2[0]) + t[0];
      v[ 4] = (f_1*v1[1]+f*v2[1]) + t[1];
      v[ 5] = (f_1*v1[2]+f*v2[2]) + t[2];
      
      f = 0.4F;
      f_1 = 1.0F-f;
      
      v[ 6] = (f_1*v1[0]+f*v2[0]) + t[0];
      v[ 7] = (f_1*v1[1]+f*v2[1]) + t[1];
      v[ 8] = (f_1*v1[2]+f*v2[2]) + t[2];
      
      v+=9;
      n++;

    }
    break;
  case 2:

    t[0] = p2[0]*tube_size*2;
    t[1] = p2[1]*tube_size*2;
    t[2] = p2[2]*tube_size*2;
    
    v[0] = color[0];
    v[1] = color[1];
    v[2] = color[2];
    
    v[3] = (v2[0]+v1[0])/2.0F;
    v[4] = (v2[1]+v1[1])/2.0F;
    v[5] = (v2[2]+v1[2])/2.0F;

    v[6] = v2[0];
    v[7] = v2[1];
    v[8] = v2[2];
    
    v[9] = color[0];
    v[10] = color[1];
    v[11] = color[2];
    
    f = 0.60F;
    f_1 = 1.0F-f;

    v[12] = (f_1*v1[0]+f*v2[0]) - t[0];
    v[13] = (f_1*v1[1]+f*v2[1]) - t[1];
    v[14] = (f_1*v1[2]+f*v2[2]) - t[2];

    f = 0.84F;
    f_1 = 1.0F-f;

    v[15] = (f_1*v1[0]+f*v2[0]) - t[0];
    v[16] = (f_1*v1[1]+f*v2[1]) - t[1];
    v[17] = (f_1*v1[2]+f*v2[2]) - t[2];

    v+=18;
    n+=2;


    if(double_sided) {

      v[ 0] = color[0];
      v[ 1] = color[1];
      v[ 2] = color[2];
      
      f = 0.60F;
      f_1 = 1.0F-f;
      
      v[ 3] = (f_1*v1[0]+f*v2[0]) + t[0];
      v[ 4] = (f_1*v1[1]+f*v2[1]) + t[1];
      v[ 5] = (f_1*v1[2]+f*v2[2]) + t[2];
      
      f = 0.84F;
      f_1 = 1.0F-f;
      
      v[ 6] = (f_1*v1[0]+f*v2[0]) + t[0];
      v[ 7] = (f_1*v1[1]+f*v2[1]) + t[1];
      v[ 8] = (f_1*v1[2]+f*v2[2]) + t[2];
      
      v+=9;
      n++;

    }

    break;
  }
  *v_ptr = v;
  *n_ptr = n;

}

void RepWireBondFree(RepWireBond *I)
{
  FreeP(I->VP);
  FreeP(I->V);
  RepPurge(&I->R);
  OOFreeP(I);
}

void RepWireBondRender(RepWireBond *I,CRay *ray,Pickable **pick)
{
  PyMOLGlobals *G=I->R.G;
  float *v=I->V;
  int c=I->N;
  unsigned int i,j;
  Pickable *p;

  if(ray) {

    float radius;
    
    if(I->Radius==0.0F) {
      radius = ray->PixelRadius*I->Width/2.0F;
    } else {
        radius = I->Radius;
    }
    
	 v=I->V;
	 c=I->N;
	 
	 while(c--) {
      /*      printf("%8.3f %8.3f %8.3f   %8.3f %8.3f %8.3f \n",v[3],v[4],v[5],v[6],v[7],v[8]);*/
      ray->fSausage3fv(ray,v+3,v+6,radius,v,v);
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

    int use_dlst;
    use_dlst = (int)SettingGet(G,cSetting_use_display_lists);
    if(use_dlst&&I->R.displayList) {
      glCallList(I->R.displayList);
    } else { 

      if(use_dlst) {
        if(!I->R.displayList) {
          I->R.displayList = glGenLists(1);
          if(I->R.displayList) {
            glNewList(I->R.displayList,GL_COMPILE_AND_EXECUTE);
          }
        }
      }
      
      glLineWidth(I->Width);
      
      v=I->V;
      c=I->N;
      
      glDisable(GL_LIGHTING); 
      glBegin(GL_LINES);	 
      SceneResetNormal(G,true);
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
      if(use_dlst&&I->R.displayList) {
        glEndList();
      }
    } 
  }
}

Rep *RepWireBondNew(CoordSet *cs)
{
  PyMOLGlobals *G=cs->G;
  ObjectMolecule *obj;
  int a,a1,a2,c1,c2,s1,s2,b1,b2,ord;
  BondType *b;
  int half_bonds,*other=NULL;
  float valence;
  float *v,*v0,*v1,*v2,h[3];
  int visFlag;
  int maxSegment = 0;
  int maxBond = 0;
  float tmpColor[3];
  int valence_flag = false;
  Pickable *rp;
  AtomInfoType *ai1,*ai2;
  int cartoon_side_chain_helper = 0;
  OOAlloc(G,RepWireBond);
  obj = cs->Obj;

  PRINTFD(G,FB_RepWireBond)
    " RepWireBondNew-Debug: entered.\n"
    ENDFD;

  visFlag=false;
  b=obj->Bond;
  if(obj->RepVisCache[cRepLine])
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

  valence = SettingGet_f(G,cs->Setting,obj->Obj.Setting,cSetting_valence);
  if(valence==1.0F) /* backwards compatibility... */
    valence = SettingGet_f(G,cs->Setting,obj->Obj.Setting,cSetting_valence_size);
  valence_flag = (valence!=0.0F);
  cartoon_side_chain_helper = SettingGet_b(G,cs->Setting, obj->Obj.Setting,
                                         cSetting_cartoon_side_chain_helper);

  half_bonds = SettingGet_i(G,cs->Setting,obj->Obj.Setting,cSetting_half_bonds);

  b=obj->Bond;
  for(a=0;a<obj->NBond;a++)
    {
      b1 = b->index[0];
      b2 = b->index[1];

      
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
          if(valence_flag) {
            if((b->order>0)&&(b->order<4)) {
              maxSegment+=2*b->order;
            } else if(b->order==4) { /* aromatic */
              maxSegment+=10;
            } else {
              maxSegment+=2;
            }
          } else
            maxSegment+=2;
          maxBond++;
        }
      b++;
    }
  
  RepInit(G,&I->R);

  I->R.fRender=(void (*)(struct Rep *, CRay *, Pickable **))RepWireBondRender;
  I->R.fFree=(void (*)(struct Rep *))RepWireBondFree;
  I->Width = SettingGet_f(G,cs->Setting,obj->Obj.Setting,cSetting_line_width);
  I->Radius = SettingGet_f(G,cs->Setting,obj->Obj.Setting,cSetting_line_radius);


  I->N=0;
  I->NP=0;
  I->V=NULL;
  I->VP=NULL;
  I->R.P=NULL;
  I->R.fRecolor=NULL;

  if(obj->NBond) {

    if(valence_flag) /* build list of up to 2 connected atoms for each atom */
      other=ObjectMoleculeGetPrioritizedOtherIndexList(obj,cs);
    
	 I->V=(float*)mmalloc(sizeof(float)*maxSegment*9);
	 ErrChkPtr(G,I->V);
	 	 
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
					 

                if(cartoon_side_chain_helper) {
                  register AtomInfoType *ati1=obj->AtomInfo+b1;
                  register AtomInfoType *ati2=obj->AtomInfo+b2;
                  if(ati1->visRep[cRepCartoon]&&
                     ati2->visRep[cRepCartoon]&&
                     (!ati1->hetatm)&&(!ati2->hetatm)) {
                    register char *name1=ati1->name;
                    register int prot1=ati1->protons;
                    register char *name2=ati2->name;
                    register int prot2=ati2->protons;

                    if(prot1 == cAN_C) { 
                      if((name1[1]=='A')&&(name1[0]=='C')&&(!name1[2])) { /* CA */
                        if(prot2 == cAN_C) { 
                          if((name2[1]=='B')&&(name2[0]=='C')&&(!name2[2]))
                            c1 = c2;  /* CA-CB */
                          else if((!name2[1])&&(name2[0]=='C'))
                            s1 = s2 = 0; /* suppress CA-C */
                        } else if(prot2 == cAN_H) 
                          s1 = s2 = 0; /* suppress all CA-hydrogens */
                      }
                    } else if(prot1 == cAN_N) { 
                      if((!name1[1])&&(name1[0]=='N')) { /* N */
                        if(prot2 == cAN_C) {
                          if((name2[1]=='D')&&(name2[0]=='C')&&(!name2[2])) 
                            c1 = c2; /* N->CD in PRO */
                          else if((name2[1]=='A')&&(name2[0]=='C')&&(!name2[2]))
                            s1 = s2 = 0; /* suppress N-CA */
                          else if((!name2[1])&&(name2[0]=='C'))
                            s1 = s2 = 0; /* suppress N-C */
                        } else if(prot2 == cAN_H)
                          s1 = s2 = 0; /* suppress all N-hydrogens */
                      }
                    } else if((prot1 == cAN_O)&&(prot2 == cAN_C)) { 
                      if((!name2[1])&&(name2[0]=='C')&&
                         (((!name1[1])&&(name1[0]=='O'))||
                          ((name1[3]==0)&&(name1[2]=='T')&&(name1[1]=='X')&&(name1[0]=='O'))))
                        s1 = s2 = 0; /* suppress C-O,OXT */
                    }
                    
                    if(prot2 == cAN_C) {
                      if((name2[1]=='A')&&(name2[0]=='C')&&(!name2[2])) { /* CA */
                        if(prot1 == cAN_C) { 
                          if((name1[1]=='B')&&(name1[0]=='C')&&(!name1[2]))
                            c2 = c1; /* CA-CB */
                          else if((!name1[1])&&(name1[0]=='C'))
                            s1 = s2 = 0; /* suppress CA-C */
                        } else if(prot1 == cAN_H) 
                          s1 = s2 = 0; /* suppress all CA-hydrogens */
                      }
                    } else if(prot2 == cAN_N) {
                      if((!name2[1])&&(name2[0]=='N')) { /* N */
                        if(prot1 == cAN_C) { 
                          if((name1[1]=='D')&&(name1[0]=='C')&&(!name1[2]))
                            c2 = c1; /* N->CD in PRO */
                          else if((name1[1]=='A')&&(name1[0]=='C')&&(!name1[2]))
                            s1 = s2 = 0; /* suppress N-CA */
                          else if((!name1[1])&&(name1[0]=='C'))
                            s1 = s2 = 0; /* suppress N-C */
                        } else if(prot1 == cAN_H)
                          s1 = s2 = 0; /* suppress all N-hydrogens */
                      }
                    } else if((prot2 == cAN_O)&&(prot1 == cAN_C)) {
                      if((!name1[1])&&(name1[0]=='C')&&
                         (((!name2[1])&&(name2[0]=='O'))||
                          ((name2[3]==0)&&(name2[2]=='T')&&(name2[1]=='X')&&(name2[0]=='O'))))
                        s1 = s2 = 0; /* suppress C-O,OXT */
                    }
                  }
                }

					 if((c1==c2)&&s1&&s2&&(!ColorCheckRamped(G,c1))) {
						

						v0 = ColorGet(G,c1);

                  if((valence_flag)&&(ord>1)&&(ord<4)) {
                    RepValence(v,v1,v2,other,a1,a2,cs->Coord,v0,ord,valence);
                    v+=ord*9;
                    I->N+=ord;
                  } else if(valence_flag&&(ord==4)) { /* aromatic */
                    RepAromatic(v1,v2,other,a1,a2,cs->Coord,v0,valence,0,&v,&I->N);
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
                      
                      if(ColorCheckRamped(G,c1)) {
                        ColorGetRamped(G,c1,v1,tmpColor);
                        v0=tmpColor;
                      } else {
                        v0 = ColorGet(G,c1);
                      }

                      if((valence_flag)&&(ord>1)&&(ord<4)) {
                        RepValence(v,v1,h,other,a1,a2,cs->Coord,v0,ord,valence);
                        v+=ord*9;
                        I->N+=ord;
							 } else if(valence_flag&&(ord==4)) {
                        RepAromatic(v1,v2,other,a1,a2,cs->Coord,v0,valence,1,&v,&I->N);
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
                      if(ColorCheckRamped(G,c2)) {
                        ColorGetRamped(G,c2,v2,tmpColor);
                        v0 = tmpColor;
                      } else {
                        v0 = ColorGet(G,c2);
                      }
                      if((valence_flag)&&(ord>1)&&(ord<4)) {
                        RepValence(v,h,v2,other,a1,a2,cs->Coord,v0,ord,valence);
                        v+=ord*9;
                        I->N+=ord;
							 } else if(valence_flag&&(ord==4)) {
                        RepAromatic(v1,v2,other,a1,a2,cs->Coord,v0,valence,2,&v,&I->N);
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

	 I->V = ReallocForSure(I->V,float,(v-I->V));

	 /* now create pickable verson */

	 if(SettingGet_f(G,cs->Setting,obj->Obj.Setting,cSetting_pickable)) {
		I->VP=(float*)mmalloc(sizeof(float)*maxBond*6*2);
		ErrChkPtr(G,I->VP);
		
		I->R.P=Alloc(Pickable,2*maxBond+1);
		ErrChkPtr(G,I->R.P);
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
		I->VP = ReallocForSure(I->VP,float,(v-I->VP));
	 }
  }
  FreeP(other);
  return((void*)(struct Rep*)I);
}

static void RepValence(float *v,float *v1,float *v2,int *other,
                int a1,int a2,float *coord,float *color,int ord,
                float tube_size)
{

  float d[3],t[3],p0[3],p1[3],p2[3],*vv;
  int a3;

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
  
  /* need a prioritized third atom to get planarity */

  a3 = ObjectMoleculeGetPrioritizedOther(other,a1,a2,NULL);

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


