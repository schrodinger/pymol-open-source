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
  float linewidth;
  float radius;
  int N;
  int NS;
  int NP;
} RepRibbon;

#include"ObjectMolecule.h"

void RepRibbonFree(RepRibbon *I);

void RepRibbonInit(void)
{
}

void RepRibbonFree(RepRibbon *I)
{
  FreeP(I->V);
  RepPurge(&I->R);
  OOFreeP(I);
}

static void RepRibbonRender(RepRibbon *I,RenderInfo *info)
{
  CRay *ray = info->ray;
  Picking **pick = info->pick;
  PyMOLGlobals *G=I->R.G;
  float *v=I->V;
  int c=I->N;
  Pickable *p;
  int i,j,ip;
  int last;

  /* 

  v[0] = index1
  v[1-3] = color1
  v[4-6] = vertex1
  v[7] = index2
  v[8-10] = color2
  v[11-13] = vertex2
  v[14] = radius

  */

  if(ray) {

    float radius;
    
    if(I->radius==0.0F) {
      radius = ray->PixelRadius*I->linewidth/2.0F;
    } else {
      radius = I->radius;
    }

    PRINTFD(G,FB_RepRibbon)
      " RepRibbonRender: rendering raytracable...\n"
      ENDFD;

    if(c>0) {
      while(c--) {
        ray->fSausage3fv(ray,v+4,v+11,radius,v+1,v+8);
        v+=18;
      }
    }
  } else if(G->HaveGUI && G->ValidContext) {
    if(pick) {

      PRINTFD(G,FB_RepRibbon)
        " RepRibbonRender: rendering pickable...\n"
        ENDFD;

      if(c) {
        i=(*pick)->src.index;
        p=I->R.P;
        last=-1;
        glBegin(GL_LINES);
        while(c--)
          {
            ip=(int)*(v);
            if(ip!=last) {
              i++;
              last=ip;
              if(!(*pick)[0].src.bond) {
                /* pass 1 - low order bits */
              
                glColor3ub((uchar)((i&0xF)<<4),(uchar)((i&0xF0)|0x8),(uchar)((i&0xF00)>>4)); /* we're encoding the index into the color */
                VLACheck((*pick),Picking,i);
                (*pick)[i].src = p[ip]; /* copy object and atom info */
                (*pick)[i].context = I->R.context;
              } else { 
                /* pass 2 - high order bits */
                j=i>>12;
                glColor3ub((uchar)((j&0xF)<<4),(uchar)((j&0xF0)|0x8),(uchar)((j&0xF00)>>4)); 
              }
            }	 
            glVertex3fv(v+4);
            ip=(int)*(v+7);
            if(ip!=last) {
              glVertex3fv(v+15); /* switch colors at midpoint */
              glVertex3fv(v+15);
              i++;
              last=ip;
              if(!(*pick)[0].src.bond) {
                /* pass 1 - low order bits */
              
                glColor3ub((uchar)((i&0xF)<<4),(uchar)((i&0xF0)|0x8),(uchar)((i&0xF00)>>4)); /* we're encoding the index into the color */
                VLACheck((*pick),Picking,i);
                (*pick)[i].src = p[ip]; /* copy object and atom info */
                (*pick)[i].context = I->R.context;
              } else { 
                /* pass 2 - high order bits */
                j=i>>12;
                glColor3ub((uchar)((j&0xF)<<4),(uchar)((j&0xF0)|0x8),(uchar)((j&0xF00)>>4)); 
              }
            }	 
            glVertex3fv(v+11);
            v+=18;
          }
        glEnd();
        (*pick)[0].src.index = i; /* pass the count */
      }
    } else {
      int use_dlst;
      int ribbon_smooth;
      
      ribbon_smooth=SettingGet_i(G,NULL,I->R.obj->Setting,cSetting_ribbon_smooth);
      if(!ribbon_smooth)
        glDisable(GL_LINE_SMOOTH);

      if(info->width_scale_flag) 
        glLineWidth(I->linewidth*info->width_scale);
      else
        glLineWidth(I->linewidth);

      use_dlst = (int)SettingGet(G,cSetting_use_display_lists);
      if(use_dlst&&I->R.displayList) {
        glCallList(I->R.displayList);
      } else { 

        SceneResetNormal(G,true);
        if(use_dlst) {
          if(!I->R.displayList) {
            I->R.displayList = glGenLists(1);
            if(I->R.displayList) {
              glNewList(I->R.displayList,GL_COMPILE_AND_EXECUTE);
            }
          }
        }
      
        PRINTFD(G,FB_RepRibbon)
          " RepRibbonRender: rendering GL...\n"
          ENDFD;
      
      
        if(c) {
        
          int first = true;
        
          glDisable(GL_LIGHTING);
          glBegin(GL_LINE_STRIP);
          while(c--)
            {
              if(first) {
                glColor3fv(v+1);
                glVertex3fv(v+4);
                first=false;
              } else if(
                        (v[4]!=v[-11])||
                        (v[5]!=v[-10])||
                        (v[6]!=v[-9 ])) {
                glEnd();
                glBegin(GL_LINE_STRIP);
                glColor3fv(v+1);
                glVertex3fv(v+4);
              }
              glColor3fv(v+8);
              glVertex3fv(v+11);
              v+=18;
            }
          glEnd();
          glEnable(GL_LIGHTING);
        }
        if(use_dlst&&I->R.displayList) {
          glEndList();
        }
      }
      if(SettingGetGlobal_b(G,cSetting_line_smooth))
        glEnable(GL_LINE_SMOOTH);
    }
  }
}

static float smooth(float x,float power)
{

  if(x<=0.5) {
    if(x<=0.0) x=0.0;
    return ((float)(0.5*pow(2.0*x,power)));    
  } else {
    if(x>=1.0) x=1.0;
    return ((float)(1.0-(0.5*pow(2*(1.0-x),power))));
  }
}

Rep *RepRibbonNew(CoordSet *cs,int state)
{
  PyMOLGlobals *G=cs->State.G;
  ObjectMolecule *obj;
  int a,b,a1,a2,c1,c2,*i,*s,*at,*seg,nAt,*atp;
  float *v,*v0,*v1,*v2,*v3;
  float *pv=NULL;
  float *dv=NULL;
  float *nv=NULL;
  float *tv=NULL;

  float f0,f1,f2,f3,f4;
  float *d;
  float *dl=NULL;
  int nSeg;
  int sampling;
  float  power_a = 5;
  float power_b = 5;
  float throw;
  int visFlag;
  float dev;
  int trace,trace_mode;
  int ribbon_color;
  int na_mode;
  AtomInfoType *ai,*last_ai=NULL;
  AtomInfoType *trailing_O3p_ai=NULL, *leading_O5p_ai=NULL;
  int trailing_O3p_a = 0, leading_O5p_a = 0, leading_O5p_a1 =0;

  Pickable *rp=NULL;
  OOAlloc(G,RepRibbon);

  obj = cs->Obj;
  visFlag=false;
  if(obj->RepVisCache[cRepRibbon])
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

  RepInit(G,&I->R);
  power_a=SettingGet_f(G,cs->Setting,obj->Obj.Setting,cSetting_ribbon_power);
  power_b=SettingGet_f(G,cs->Setting,obj->Obj.Setting,cSetting_ribbon_power_b);
  throw=SettingGet_f(G,cs->Setting,obj->Obj.Setting,cSetting_ribbon_throw);
  trace=SettingGet_i(G,cs->Setting,obj->Obj.Setting,cSetting_ribbon_trace_atoms);
  trace_mode=SettingGet_i(G,cs->Setting,obj->Obj.Setting,cSetting_trace_atoms_mode);
  na_mode = SettingGet_i(G,cs->Setting,obj->Obj.Setting,cSetting_ribbon_nucleic_acid_mode);

  ribbon_color=SettingGet_color(G,cs->Setting,obj->Obj.Setting,cSetting_ribbon_color);

  sampling=SettingGet_i(G,cs->Setting,obj->Obj.Setting,cSetting_ribbon_sampling);
  if(sampling<1) sampling=1;
  I->radius=SettingGet_f(G,cs->Setting,obj->Obj.Setting,cSetting_ribbon_radius);

  I->R.fRender=(void (*)(struct Rep *, RenderInfo *))RepRibbonRender;
  I->R.fFree=(void (*)(struct Rep *))RepRibbonFree;
  I->R.fRecolor=NULL;
  I->R.obj = (CObject*)obj;
  I->linewidth = SettingGet_f(G,cs->Setting,obj->Obj.Setting,cSetting_ribbon_width);
  I->R.context.object = (void*)obj;
  I->R.context.state = state;
  
  /* find all of the CA points */

  at = Alloc(int,cs->NAtIndex*2);
  pv = Alloc(float,cs->NAtIndex*6);
  seg = Alloc(int,cs->NAtIndex*2);
  
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
      if(a>=0) {
        ai = obj->AtomInfo+a1;
        if(obj->AtomInfo[a1].visRep[cRepRibbon]) {
          if(trace||((obj->AtomInfo[a1].protons==cAN_C)&&
                     (WordMatch(G,"CA",obj->AtomInfo[a1].name,1)<0)&&
                     !AtomInfoSameResidueP(G,last_ai,ai))) {
            PRINTFD(G,FB_RepRibbon)
              " RepRibbon: found atom in %s; a1 %d a2 %d\n",obj->AtomInfo[a1].resi,a1,a2
              ENDFD;
              
            if(a2>=0) {
              /*						if((abs(obj->AtomInfo[a1].resv-obj->AtomInfo[a2].resv)>1)||
                                        (obj->AtomInfo[a1].chain[0]!=obj->AtomInfo[a2].chain[0])||
                                        (!WordMatch(G,obj->AtomInfo[a1].segi,obj->AtomInfo[a2].segi,1)))*/
              if(trace) {
                if(!AtomInfoSequential(G,obj->AtomInfo+a2,obj->AtomInfo+a1,trace_mode))
                  a2=-1;
              } else {
                if(!ObjectMoleculeCheckBondSep(obj,a1,a2,3)) /* CA->N->C->CA = 3 bonds */
                  a2=-1;
              }
            }
            PRINTFD(G,FB_RepRibbon)
              " RepRibbon: found atom in %s; a1 %d a2 %d\n",obj->AtomInfo[a1].resi,a1,a2
              ENDFD;
            last_ai = ai;
            if(a2<0) nSeg++;
            *(s++) = nSeg;
            nAt++;
            *(i++)=a;
            v1 = cs->Coord+3*a;		
            *(v++)=*(v1++);
            *(v++)=*(v1++);
            *(v++)=*(v1++);
            
            a2=a1;
          } else if(
                    (((na_mode!=1)&&(ai->protons==cAN_P) &&
                      (WordMatch(G,"P",ai->name,1)<0) ) ||
                     ((na_mode==1)&&(ai->protons==cAN_C) &&
                      (WordMatchExact(G,"C4*",ai->name,1) ||
                       WordMatchExact(G,"C4'",ai->name,1))))&&
                    !AtomInfoSameResidueP(G,last_ai,ai)) {
            if(a2>=0) {
              if(!ObjectMoleculeCheckBondSep(obj,a1,a2,6)) { /* six bonds between phosphates */
                if(trailing_O3p_ai&&((na_mode==2)||(na_mode==4))) {
                  /* 3' end of nucleic acid */
                  *(s++) = nSeg;
                  nAt++;
                  *(i++)=trailing_O3p_a;
                  v1 = cs->Coord+3*trailing_O3p_a;		
                  *(v++)=*(v1++);
                  *(v++)=*(v1++);
                  *(v++)=*(v1++);
                }
                a2=-1;
              }
            }

            trailing_O3p_ai = NULL;

            if(leading_O5p_ai&&(a2<0)&&
               ((na_mode==3)||(na_mode==4))) {
              if((!AtomInfoSameResidueP(G,ai,leading_O5p_ai))&&
                 ObjectMoleculeCheckBondSep(obj,a1,leading_O5p_a1,5)) {
                nSeg++;
                *(s++) = nSeg;
                nAt++;
                *(i++)= leading_O5p_a;
                v1 = cs->Coord+3*leading_O5p_a;		
                *(v++)=*(v1++);
                *(v++)=*(v1++);
                *(v++)=*(v1++);
                a2 = leading_O5p_a1;
              }
            }
            leading_O5p_ai = NULL;
            last_ai = ai;
            if(a2<0) nSeg++;
            *(s++) = nSeg;
            nAt++;
            *(i++)=a;
            v1 = cs->Coord+3*a;		
            *(v++)=*(v1++);
            *(v++)=*(v1++);
            *(v++)=*(v1++);
            
            a2=a1;
          } else if((a2>=0)&&
                    last_ai&&
                    (ai->protons==cAN_O)&&
                    (last_ai->protons==cAN_P)&&
                    ((na_mode==2)||(na_mode==4))&&
                    (WordMatchExact(G,"O3'",ai->name,1)||
                     WordMatchExact(G,"O3*",ai->name,1))&&
                    AtomInfoSameResidueP(G,last_ai,ai)&&
                    ObjectMoleculeCheckBondSep(obj,a1,a2,5)) {
            trailing_O3p_ai = ai;
            trailing_O3p_a = a;
          } else if((ai->protons==cAN_O)&&
                    ((na_mode==3)||(na_mode==4))&&
                    (WordMatchExact(G,"O5'",ai->name,1)||
                     WordMatchExact(G,"O5*",ai->name,1))) {
            leading_O5p_ai = ai;
            leading_O5p_a = a;
            leading_O5p_a1 = a1;
          }
        }
      }
      
    }
  if(trailing_O3p_ai&&((na_mode==2)||(na_mode==4))) {
    /* 3' end of nucleic acid */
    *(s++) = nSeg;
    nAt++;
    *(i++)=trailing_O3p_a;
    v1 = cs->Coord+3*trailing_O3p_a;		
    *(v++)=*(v1++);
    *(v++)=*(v1++);
    *(v++)=*(v1++);
  }
  PRINTFD(G,FB_RepRibbon)
    " RepRibbon: nAt %d\n",nAt
    ENDFD;
  
  if(nAt)
    {
      /* compute differences and normals */

      s=seg;
      v=pv;
		
      dv = Alloc(float,nAt*6);
      nv = Alloc(float,nAt*6);
      dl = Alloc(float,nAt*2);
      v1=dv;
      v2=nv;
      d=dl;
		
      for(a=0;a<(nAt-1);a++)
        {
          if(*s==*(s+1))
            {
              float d_1;
              subtract3f(v+3,v,v1);
              *d = (float)length3f(v1);
              if(*d>R_SMALL4) {
                d_1 = 1.0F/(*d);
                scale3f(v1,d_1,v2);
              } else if(a)  {
                copy3f(v2-3,v2); 
              } else {
                zero3f(v2);
              }
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
		
      tv = Alloc(float,nAt*6+6);
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
    I->R.P=Alloc(Pickable,2*nAt+2);
    ErrChkPtr(G,I->R.P);
    I->R.P[0].index=nAt;
    rp = I->R.P + 1; /* skip first record! */
  }

  I->V=(float*)mmalloc(sizeof(float)*2*cs->NIndex*18*sampling);
  ErrChkPtr(G,I->V);

  I->N=0;
  I->NP=0;
  v=I->V;
  if(nAt) {
    v1=pv; /* points */
    v2=tv; /* tangents */
    v3=dv; /* direction vector */
    d = dl;
    s=seg;
    atp=at;
    rp->index = cs->IdxToAtm[*atp];
    if(obj->AtomInfo[cs->IdxToAtm[*atp]].masked)
      rp->index = -1; 
    rp->bond = -1;
    I->NP++;
    rp++;
    for(a=0;a<(nAt-1);a++) {

      rp->index = cs->IdxToAtm[*(atp+1)]; /* store pickable for n+2 */
      if(obj->AtomInfo[cs->IdxToAtm[*atp]].masked)
        rp->index = -1;
      rp->bond = -1;
      rp++;
      I->NP++;
         
      PRINTFD(G,FB_RepRibbon)
        " RepRibbon: seg %d *s %d , *(s+1) %d\n",a,*s,*(s+1)
        ENDFD;

      if(*s==*(s+1)) {
        int atom_index1 = cs->IdxToAtm[*atp];
        int atom_index2 = cs->IdxToAtm[*(atp+1)];
        
        c1=*(cs->Color+*atp);
        c2=*(cs->Color+*(atp+1));
          
        if(ribbon_color>=0) {
          c1 = (c2 = ribbon_color);
        }
          
        AtomInfoGetSetting_color(G,obj->AtomInfo + atom_index1, cSetting_ribbon_color, c1, &c1);
        AtomInfoGetSetting_color(G,obj->AtomInfo + atom_index2, cSetting_ribbon_color, c2, &c2);

        dev = throw*(*d);
          
        for(b=0;b<sampling;b++) { /* needs optimization */
                  
          f0=((float)b)/sampling; /* fraction of completion */
          f0 = smooth(f0,power_a);  /* bias sampling towards the center of the curve */
                  
          if(f0<0.5) {
            v0 = ColorGet(G,c1);
          } else {
            v0 = ColorGet(G,c2);
          }

          /* store index */
          if(f0<0.5F)
            *(v++)=(float)I->NP-1;
          else
            *(v++)=(float)I->NP;
                
          /* store colors */

          *(v++)=*(v0++);
          *(v++)=*(v0++);
          *(v++)=*(v0);

          /* start of line/cylinder */

          f1=1.0F-f0;
          f2=smooth(f0,power_b);
          f3=smooth(f1,power_b);
          f4 = dev*f2*f3; /* displacement magnitude */
                
          *(v++)=f1*v1[0]+f0*v1[3]+
            f4*( f3*v2[0]-f2*v2[3] );
                
          *(v++)=f1*v1[1]+f0*v1[4]+
            f4*( f3*v2[1]-f2*v2[4] );
                
          *(v++)=f1*v1[2]+f0*v1[5]+
            f4*( f3*v2[2]-f2*v2[5] );
                
          f0=((float)b+1)/sampling;
          f0 = smooth(f0,power_a);

          if(f0<0.5) {
            v0 = ColorGet(G,c1);
          } else {
            v0 = ColorGet(G,c2);
          }

          /* store index */
          if(f0<0.5)
            *(v++)=(float)I->NP-1;
          else
            *(v++)=(float)I->NP;
                
          /* store colors */

          *(v++)=*(v0++);
          *(v++)=*(v0++);
          *(v++)=*(v0);

          /* end of line/cylinder */

          f1=1.0F-f0;
          f2=smooth(f0,power_b);
          f3=smooth(f1,power_b);
          f4 = dev*f2*f3; /* displacement magnitude */
                
          *(v++)=f1*v1[0]+f0*v1[3]+
            f4*( f3*v2[0]-f2*v2[3] );
                
          *(v++)=f1*v1[1]+f0*v1[4]+
            f4*( f3*v2[1]-f2*v2[4] );
                
          *(v++)=f1*v1[2]+f0*v1[5]+
            f4*( f3*v2[2]-f2*v2[5] );

          v++; /* radius no longer stored here... */

          average3f(v-4,v-11,v);

          v+=3;

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
	 I->V=ReallocForSure(I->V,float,(v-I->V));
  else
	 I->V=ReallocForSure(I->V,float,1);

  return((void*)(struct Rep*)I);
}


void RepRibbonRenderImmediate(CoordSet *cs, RenderInfo *info)
{
  /* performance optimized to provide a simple C-alpha trace -- no smoothing */

  PyMOLGlobals *G=cs->State.G;
  if(info->ray || info->pick || (!(G->HaveGUI && G->ValidContext)) )
    return;
  else {
    ObjectMolecule *obj = cs->Obj;
    int active = false;
    int nAtIndex = cs->NAtIndex;
    int a;
    AtomInfoType *obj_AtomInfo = obj->AtomInfo;
    AtomInfoType *ai,*last_ai = NULL;
    int trace=SettingGet_i(G,cs->Setting,obj->Obj.Setting,cSetting_ribbon_trace_atoms);
    int trace_mode=SettingGet_i(G,cs->Setting,obj->Obj.Setting,cSetting_trace_atoms_mode);
    int na_mode = SettingGet_i(G,cs->Setting,obj->Obj.Setting,cSetting_ribbon_nucleic_acid_mode);
    float linewidth = SettingGet_f(G,cs->Setting,obj->Obj.Setting,cSetting_ribbon_width);
    int a1,a2=-1;
    int color,last_color = -9;

    glLineWidth(linewidth);
    glDisable(GL_LIGHTING); 
    SceneResetNormal(G,true);      
    glBegin(GL_LINE_STRIP);
    for(a1=0;a1<nAtIndex;a1++) {
      if(obj->DiscreteFlag) {
        if(cs==obj->DiscreteCSet[a1]) 
          a=obj->DiscreteAtmToIdx[a1];
        else 
          a=-1;
      } else 
        a=cs->AtmToIdx[a1];
      if(a>=0) {
        ai = obj->AtomInfo+a1;
        if(ai->visRep[cRepRibbon]) {
          if(trace || ((ai->protons==cAN_C)&&
                       (WordMatch(G,"CA",ai->name,1)<0)&&
                       !AtomInfoSameResidueP(G,last_ai,ai))) {
            if(a2>=0) {
              if(trace) {
                if(!AtomInfoSequential(G,obj->AtomInfo+a2,obj->AtomInfo+a1,trace_mode))
                  a2=-1;
              } else {
                if(!ObjectMoleculeCheckBondSep(obj,a1,a2,3)) /* CA->N->C->CA = 3 bonds */
                  a2=-1;
              }
            }
            if(a2==-1) {
              glEnd();
              glBegin(GL_LINE_STRIP);
            }
            color = ai->color;
            if(color!=last_color) {
              last_color = color;
              glColor3fv(ColorGet(G,color));
            }
            glVertex3fv(cs->Coord+3*a);
            active = true;
            last_ai = ai;
            a2 = a1;
          } else if((((na_mode!=1)&&(ai->protons==cAN_P) &&
                      (WordMatch(G,"P",ai->name,1)<0) ) ||
                     ((na_mode==1)&&(ai->protons==cAN_C) &&
                      (WordMatchExact(G,"C4*",ai->name,1) ||
                       WordMatchExact(G,"C4'",ai->name,1))))&&
                    !AtomInfoSameResidueP(G,last_ai,ai)) {
            if(a2>=0) {
              if(!ObjectMoleculeCheckBondSep(obj,a1,a2,6)) { /* six bonds between phosphates */
                a2=-1;
              }
            }
            
            if(a2==-1) {
              glEnd();
              glBegin(GL_LINE_STRIP);
            }
            color = ai->color;
            if(color!=last_color) {
              last_color = color;
              glColor3fv(ColorGet(G,color));
            }
            glVertex3fv(cs->Coord+3*a);
            active = true;
            last_ai = ai;
            a2=a1;
          }
        }
      }
    }
    glEnd();
    glEnable(GL_LIGHTING);
    if(!active)
      cs->Active[cRepRibbon] = false;
  }
}



