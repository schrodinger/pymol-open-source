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
#include"Debug.h"
#include"Parse.h"
#include"OOMac.h"
#include"Vector.h"
#include"MemoryDebug.h"
#include"Err.h"
#include"Map.h"
#include"Selector.h"
#include"ObjectMolecule.h"
#include"Ortho.h"
#include"Util.h"
#include"Vector.h"
#include"Selector.h"
#include"Matrix.h"
#include"Scene.h"
#include"P.h"
#include"PConv.h"
#include"Executive.h"
#include"Setting.h"
#include"Sphere.h"
#include"main.h"
#include"CGO.h"
#include"Raw.h"
#include"Editor.h"
#include"Selector.h"
#include"Sculpt.h"

#define cMaxNegResi 100

#define ntrim ParseNTrim
#define nextline ParseNextLine
#define ncopy ParseNCopy
#define nskip ParseNSkip

void ObjectMoleculeRender(ObjectMolecule *I,int frame,CRay *ray,Pickable **pick,int pass);
void ObjectMoleculeCylinders(ObjectMolecule *I);
CoordSet *ObjectMoleculeMMDStr2CoordSet(char *buffer,AtomInfoType **atInfoPtr);

int ObjectMoleculeDoesAtomNeighborSele(ObjectMolecule *I, int index, int sele);

CoordSet *ObjectMoleculePMO2CoordSet(CRaw *pmo,AtomInfoType **atInfoPtr,int *restart);

void ObjectMoleculeAppendAtoms(ObjectMolecule *I,AtomInfoType *atInfo,CoordSet *cset);
CoordSet *ObjectMoleculeMOLStr2CoordSet(char *buffer,AtomInfoType **atInfoPtr);

void ObjectMoleculeUpdate(ObjectMolecule *I);
int ObjectMoleculeGetNFrames(ObjectMolecule *I);

void ObjectMoleculeDescribeElement(ObjectMolecule *I,int index,char *buffer);

void ObjectMoleculeSeleOp(ObjectMolecule *I,int sele,ObjectMoleculeOpRec *op);

void ObjectMoleculeTransformTTTf(ObjectMolecule *I,float *ttt,int state);


CoordSet *ObjectMoleculeChemPyModel2CoordSet(PyObject *model,AtomInfoType **atInfoPtr);

int ObjectMoleculeGetAtomGeometry(ObjectMolecule *I,int state,int at);
void ObjectMoleculeBracketResidue(ObjectMolecule *I,AtomInfoType *ai,int *st,int *nd);

void ObjectMoleculeAddSeleHydrogens(ObjectMolecule *I,int sele);


CoordSet *ObjectMoleculeXYZStr2CoordSet(char *buffer,AtomInfoType **atInfoPtr);
CSetting **ObjectMoleculeGetSettingHandle(ObjectMolecule *I,int state);
void ObjectMoleculeInferAmineGeomFromBonds(ObjectMolecule *I,int state);
CoordSet *ObjectMoleculeTOPStr2CoordSet(char *buffer,
                                        AtomInfoType **atInfoPtr);

ObjectMolecule *ObjectMoleculeReadTOPStr(ObjectMolecule *I,char *TOPStr,int frame,int discrete);

void ObjectMoleculeInferHBondFromChem(ObjectMolecule *I);

#define MAX_BOND_DIST 50

#if 0
static void dump_jxn(char *lab,char *q)
{
  char cc[MAXLINELEN];  
  printf("\n%s %p\n",lab,q);
  q=q-150;
  q=nextline(q);
  ncopy(cc,q,100);
  printf("0[%s]\n",cc);
  q=nextline(q);
  ncopy(cc,q,100);
  printf("1[%s]\n",cc);
  q=nextline(q);
  ncopy(cc,q,100);
  printf("2[%s]\n",cc);
  q=nextline(q);
  ncopy(cc,q,100);
  printf("3[%s]\n",cc);

}
#endif

static char *skip_fortran(int num,int per_line,char *p)
{
  int a,b;
  b=0;
  for(a=0;a<num;a++) {
    if((++b)==per_line) {
      b=0;
      p=nextline(p);
    }
  }  
  if(b) p=nextline(p);
  return(p);
}

void M4XAnnoInit(M4XAnnoType *m4x)
{
  UtilZeroMem((char*)m4x,sizeof(M4XAnnoType));
}

void M4XAnnoPurge(M4XAnnoType *m4x)
{
  int c;
  if(m4x) {
    for(c=0;c<m4x->n_context;c++) {
      VLAFreeP(m4x->context[c].hbond);
      VLAFreeP(m4x->context[c].nbond);
      VLAFreeP(m4x->context[c].site);
      VLAFreeP(m4x->context[c].ligand);
      VLAFreeP(m4x->context[c].water);
    }
    if(m4x->align) {
      M4XAlignPurge(m4x->align);
    }
    VLAFreeP(m4x->context);
  }
}

void M4XAlignInit(M4XAlignType *align)
{
  UtilZeroMem((char*)align,sizeof(M4XAlignType));
  align->id_at_point = VLACalloc(int,100);
  align->fitness = VLAlloc(float,100);
}

void M4XAlignPurge(M4XAlignType *align)
{
  VLAFreeP(align->id_at_point);
  VLAFreeP(align->fitness);
  FreeP(align);
}

void ObjectMoleculeOpRecInit(ObjectMoleculeOpRec *op)
{
  UtilZeroMem((char*)op,sizeof(ObjectMoleculeOpRec));
}

/*========================================================================*/
ObjectMolecule *ObjectMoleculeLoadTRJFile(ObjectMolecule *I,char *fname,int frame,
                                          int interval,int average,int start,
                                          int stop,int max,char *sele,int image,
                                          float *shift)
{
  int ok=true;
  FILE *f;
  char *buffer,*p;
  char cc[MAXLINELEN];  
  int n_read;
  int to_go;
  int skip_first_line = true;
  int periodic=false;
  int angles=true;
  float f0,f1,f2,*fp;
  float box[3],angle[3];
  float r_cent[3],r_trans[3];
  int r_act,r_val,r_cnt;
  float *r_fp_start=NULL,*r_fp_stop=NULL;
  int a,b,c,i;
  int *to,at_i;
  int zoom_flag=false;
  int cnt=0;
  int n_avg=0;
  int icnt;
  int ncnt=0;
  int sele0 = SelectorIndexByName(sele);
  int *xref = NULL;
  float zerovector[3]={0.0,0.0,0.0};
  CoordSet *cs = NULL;

  if(!shift)
    shift = zerovector;
  if(interval<1)
    interval=1;

  icnt=interval;
  #define BUFSIZE 4194304
  #define GETTING_LOW 10000

  f=fopen(fname,"rb");
  if(!f) {
	 ok=ErrMessage("ObjectMoleculeLoadTOPFile","Unable to open file!");
  } else
	 {
      if(!I->CSTmpl) {
        PRINTFB(FB_Errors,FB_ObjectMolecule)
          " ObjMolLoadTRJFile: Missing topology"
          ENDFB;
        return(I);
      }
      cs=CoordSetCopy(I->CSTmpl);

      if(sele0>=0) { /* build array of cross-references */
        xref = Alloc(int,I->NAtom);
        c=0;
        for(a=0;a<I->NAtom;a++) {
          if(SelectorIsMember(I->AtomInfo[a].selEntry,sele0)) {
            xref[a]=c++;
          } else {
            xref[a]=-1;
          }
        }

        for(a=0;a<I->NAtom;a++) { /* now terminate the excluded references */
          if(xref[a]<0) {
            cs->AtmToIdx[a]=-1;
          } 
        }

        to=cs->IdxToAtm;
        c=0;
        for(a=0;a<cs->NIndex;a++) { /* now fix IdxToAtm, AtmToIdx,
                                       and remap xref to coordinate space */         
          at_i = cs->IdxToAtm[a];
          if(cs->AtmToIdx[at_i]>=0) {
            *(to++)=at_i;
            cs->AtmToIdx[at_i]=c;
            xref[a]=c;
            c++;
          } else {
            xref[a]=-1;
          }
        }

        cs->NIndex=c;
        cs->IdxToAtm = Realloc(cs->IdxToAtm,int,cs->NIndex+1);
        VLASize(cs->Coord,float,cs->NIndex*3);
      }
      PRINTFB(FB_ObjectMolecule,FB_Blather) 
        " ObjMolLoadTRJFile: Loading from \"%s\".\n",fname
        ENDFB;
      buffer = (char*)mmalloc(BUFSIZE+1); /* 1 MB read buffer */
      p = buffer;
      buffer[0]=0;
      n_read = 0;
      to_go=0;
      a = 0;
      b = 0;
      c = 0;
      f1=0.0;
      f2=0.0;
      while(1)
        {
          to_go = n_read-(p-buffer);
          if(to_go<GETTING_LOW) 
            if(!feof(f)) {
              if(to_go) 
                memcpy(buffer,p,to_go);
              n_read = fread(buffer+to_go,1,BUFSIZE-to_go,f);              
              n_read = to_go + n_read;
              buffer[n_read]=0;
              p = buffer;
              if(skip_first_line) {
                p=nextline(p);
                skip_first_line=false;
              }
              to_go = n_read-(p-buffer);
            }
          if(!to_go) break;
          p=ncopy(cc,p,8);
          if((++b)==10) {
            b=0;
            p=nextline(p);
          }
          f0 = f1;
          f1 = f2;
          if(sscanf(cc,"%f",&f2)==1) {
            if((++c)==3) {
              c=0;
              if((cnt+1)>=start) {
                if(icnt<=1) {
                  if(xref) { 
                    if(xref[a]>=0)
                      fp=cs->Coord+3*xref[a];
                    else 
                      fp=NULL;
                  } else {
                    fp=cs->Coord+3*a;
                  }
                  if(fp) {
                    if(n_avg) {
                      *(fp++)+=f0;
                      *(fp++)+=f1;
                      *(fp++)+=f2;
                    } else {
                      *(fp++)=f0;
                      *(fp++)=f1;
                      *(fp++)=f2;
                    }
                  }
                }
              }
              if((++a)==I->NAtom) {
                
                cnt++;
                a=0;
                if(b) p=nextline(p);
                b=0;


                if(cs->PeriodicBoxType!=cCSet_NoPeriodicity) {
                  /* read periodic box */

                  c=0;
                  periodic=true;
                  angles=true;
                  
                  p = ncopy(cc,p,8);
                  if(sscanf(cc,"%f",&box[0])!=1) 
                    periodic=false;
                  p = ncopy(cc,p,8);
                  if(sscanf(cc,"%f",&box[1])!=1)
                    periodic=false;
                  p = ncopy(cc,p,8);
                  if(sscanf(cc,"%f",&box[2])!=1)
                    periodic=false;
                  
                  p = ncopy(cc,p,8);
                  if(sscanf(cc,"%f",&angle[0])!=1)
                    angles=false;
                  
                  p = ncopy(cc,p,8);
                  if(sscanf(cc,"%f",&angle[1])!=1)
                    angles=false;
                  
                  p = ncopy(cc,p,8);
                  if(sscanf(cc,"%f",&angle[2])!=1)
                    angles=false;
                  if(periodic) {
                    if(!cs->PeriodicBox)
                      cs->PeriodicBox=CrystalNew();
                    cs->PeriodicBox->Dim[0] = box[0];
                    cs->PeriodicBox->Dim[1] = box[1];
                    cs->PeriodicBox->Dim[2] = box[2];
                    if(angles) {
                      cs->PeriodicBox->Angle[0] = angle[0];
                      cs->PeriodicBox->Angle[1] = angle[1];
                      cs->PeriodicBox->Angle[2] = angle[2];
                    } 
                    CrystalUpdate(cs->PeriodicBox);
                    /*                    CrystalDump(cs->PeriodicBox);*/
                    p=nextline(p);
                    b=0;
                  }

                  if(cs->PeriodicBoxType==cCSet_Octahedral)
                    periodic=false; /* can't handle this yet... */
                }
                
                if((stop>0)&&(cnt>=stop))
                  break;
                if(cnt>=start) {
                  icnt--;                      
                  if(icnt>0) {
                    PRINTFB(FB_Details,FB_ObjectMolecule)
                      " ObjectMolecule: skipping set %d...\n",cnt
                      ENDFB;
                  } else {
                    icnt=interval;
                    n_avg++;
                  }
                  
                  if(icnt==interval) {
                    if(n_avg<average) {
                      PRINTFB(FB_Details,FB_ObjectMolecule)
                        " ObjectMolecule: averaging set %d...\n",cnt
                        ENDFB;
                    } else {
                      
                      /* compute average */
                      if(n_avg>1) {
                        fp=cs->Coord;
                        for(i=0;i<cs->NIndex;i++) {
                          *(fp++)/=n_avg;
                          *(fp++)/=n_avg;
                          *(fp++)/=n_avg;
                        }
                      }
                      if(periodic&&image) { /* Perform residue-based period image transformation */
                        i = 0;
                        r_cnt = 0;
                        r_act = 0; /* 0 unspec, 1=load, 2=image, 3=leave*/
                        r_val = -1;
                        while(r_act!=3) {
                          if(i>=cs->NIndex) {
                            if(r_cnt)
                              r_act = 2; 
                            else
                              r_act = 3;
                          }
                          if(r_act==0) {
                            /* start new residue */
                            r_cnt = 0;
                            r_act = 1; /* now load */
                          }
                          if(r_act==1) {
                            if(i<cs->NIndex) {
                              
                              /* is there a coordinate for atom? */
                              if(xref) { 
                                if(xref[i]>=0)
                                  fp=cs->Coord+3*xref[i];
                                else 
                                  fp=NULL;
                              } else {
                                fp=cs->Coord+3*i;
                              }
                              if(fp) { /* yes there is... */
                                if(r_cnt) {
                                  if(r_val!=I->AtomInfo[cs->IdxToAtm[i]].resv) {
                                    r_act=2; /* end of residue-> time to image */
                                  } else {
                                    r_cnt++;
                                    r_cent[0]+=*(fp++);
                                    r_cent[1]+=*(fp++);
                                    r_cent[2]+=*(fp++);
                                    r_fp_stop = fp; /* stop here */
                                    i++;
                                  }
                                } else {
                                  r_val = I->AtomInfo[cs->IdxToAtm[i]].resv;
                                  r_cnt++;
                                  r_fp_start = fp; /* start here */
                                  r_cent[0]=*(fp++);
                                  r_cent[1]=*(fp++);
                                  r_cent[2]=*(fp++);
                                  r_fp_stop = fp; /* stop here */
                                  i++;
                                }
                              } else {
                                i++;
                              }
                            } else {
                              r_act=2; /* image */
                            }
                          }

                          if(r_act==2) { /* time to image */
                            if(r_cnt) {
                              r_cent[0]/=r_cnt;
                              r_cent[1]/=r_cnt;
                              r_cent[2]/=r_cnt;
                              transform33f3f(cs->PeriodicBox->RealToFrac,r_cent,r_cent);
                              r_trans[0]=(float)fmod(1000.0+shift[0]+r_cent[0],1.0F);
                              r_trans[1]=(float)fmod(1000.0+shift[1]+r_cent[1],1.0F);
                              r_trans[2]=(float)fmod(1000.0+shift[2]+r_cent[2],1.0F);
                              r_trans[0]-=r_cent[0];
                              r_trans[1]-=r_cent[1];
                              r_trans[2]-=r_cent[2];
                              transform33f3f(cs->PeriodicBox->FracToReal,r_trans,r_trans);
                              fp=r_fp_start;
                              while(fp<r_fp_stop) {
                                *(fp++)+=r_trans[0];
                                *(fp++)+=r_trans[1];
                                *(fp++)+=r_trans[2];
                              }
                            }
                            r_act=0; /* reset */ 
                            r_cnt=0;
                          }
                        }
                      }

                      /* add new coord set */
                      if(cs->fInvalidateRep)
                        cs->fInvalidateRep(cs,cRepAll,cRepInvRep);
                      if(frame<0) frame=I->NCSet;
                      if(!I->NCSet) {
                        zoom_flag=true;
                      }
                      
                      VLACheck(I->CSet,CoordSet*,frame);
                      if(I->NCSet<=frame) I->NCSet=frame+1;
                      if(I->CSet[frame]) I->CSet[frame]->fFree(I->CSet[frame]);
                      I->CSet[frame] = cs;
                      ncnt++;
                      
                      if(average<2) {
                        PRINTFB(FB_Details,FB_ObjectMolecule)
                          " ObjectMolecule: read set %d into state %d...\n",cnt,frame+1
                          ENDFB;
                      } else {
                        PRINTFB(FB_Details,FB_ObjectMolecule)
                          " ObjectMolecule: averaging set %d...\n",cnt
                          ENDFB;
                        PRINTFB(FB_Details,FB_ObjectMolecule)
                          " ObjectMolecule: average loaded into state %d...\n",frame+1
                          ENDFB;
                      }
                      frame++;
                      cs = CoordSetCopy(cs);
                      n_avg=0;
                      if((stop>0)&&(cnt>=stop))
                        break;
                      if((max>0)&&(ncnt>=max))
                        break;
                    }
                  }
                } else {
                  PRINTFB(FB_Details,FB_ObjectMolecule)
                    " ObjectMolecule: skipping set %d...\n",cnt
                    ENDFB;
                }
              }
            }
          } else {
            PRINTFB(FB_Errors,FB_ObjectMolecule)
              " ObjMolLoadTRJFile-Error: Failed to read an expected coordinate value.\n    This trajectory does not match the loaded parameter/topology file.\n    Likely cause: either the atom count or the periodic box settings\n    are inconsistent between the two files.\n"
              ENDFB;
            break;
          }
        }
      FreeP(xref);
		mfree(buffer);
	 }
  if(cs)
    cs->fFree(cs);
  SceneChanged();
  SceneCountFrames();
  if(zoom_flag) 
    if(SettingGet(cSetting_auto_zoom)) {
      ExecutiveWindowZoom(I->Obj.Name,0.0,-1,0); /* auto zoom (all states) */
    }
  
  return(I);
}
ObjectMolecule *ObjectMoleculeLoadRSTFile(ObjectMolecule *I,char *fname,int frame)

{
  int ok=true;
  FILE *f;
  char *buffer,*p;
  char cc[MAXLINELEN];  
  float f0,f1,f2,*fp;
  int a,b,c;
  int zoom_flag=false;
  CoordSet *cs = NULL;
  int size;

  #define BUFSIZE 4194304
  #define GETTING_LOW 10000

  f=fopen(fname,"rb");
  if(!f)
	 ok=ErrMessage("ObjectMoleculeLoadRSTFile","Unable to open file!");
  else
	 {
      if(!I->CSTmpl) {
        PRINTFB(FB_Errors,FB_ObjectMolecule)
          " ObjMolLoadTRJFile: Missing topology"
          ENDFB;
        return(I);
      }
      cs=CoordSetCopy(I->CSTmpl);
      PRINTFB(FB_ObjectMolecule,FB_Blather) 
        " ObjMolLoadTRJFile: Loading from \"%s\".\n",fname
        ENDFB;


		fseek(f,0,SEEK_END);
      size=ftell(f);
		fseek(f,0,SEEK_SET);

		buffer=(char*)mmalloc(size+255);
		ErrChkPtr(buffer);
		p=buffer;
		fseek(f,0,SEEK_SET);
		fread(p,size,1,f);
		p[size]=0;
		fclose(f);

      p=nextline(p);
      p=nextline(p);

      a = 0;
      b = 0;
      c = 0;
      f1=0.0;
      f2=0.0;
      while(*p)
        {
          p=ncopy(cc,p,12);
          if((++b)==6) {
            b=0;
            p=nextline(p);
          }
          f0 = f1;
          f1 = f2;
          if(sscanf(cc,"%f",&f2)==1) {
            if((++c)==3) {
              c=0;
              fp=cs->Coord+3*a;
              *(fp++)=f0;
              *(fp++)=f1;
              *(fp++)=f2;
              
              if((++a)==I->NAtom) {
                a=0;
                if(b) p=nextline(p);
                b=0;
                /* add new coord set */
                if(cs->fInvalidateRep)
                  cs->fInvalidateRep(cs,cRepAll,cRepInvRep);
                if(frame<0) frame=I->NCSet;
                if(!I->NCSet) {
                  zoom_flag=true;
                }
                
                VLACheck(I->CSet,CoordSet*,frame);
                if(I->NCSet<=frame) I->NCSet=frame+1;
                if(I->CSet[frame]) I->CSet[frame]->fFree(I->CSet[frame]);
                I->CSet[frame] = cs;
                
                PRINTFB(FB_Details,FB_ObjectMolecule)
                  " ObjectMolecule: read coordinates into state %d...\n",frame+1
                  ENDFB;
               
                cs = CoordSetCopy(cs);
                break;
              }
            }
          } else {
            PRINTFB(FB_Errors,FB_ObjectMolecule)
              " ObjMolLoadTRJFile: atom/coordinate mismatch.\n"
              ENDFB;
            break;
          }
        }
		mfree(buffer);
	 }
  if(cs)
    cs->fFree(cs);
  
  SceneChanged();
  SceneCountFrames();
  if(zoom_flag) 
    if(SettingGet(cSetting_auto_zoom)) {
      ExecutiveWindowZoom(I->Obj.Name,0.0,-1,0); /* auto zoom (all states) */
    }
  
  return(I);
}
static char *findflag(char *p,char *flag,char *format)
{

  char cc[MAXLINELEN];
  char pat[MAXLINELEN] = "%FLAG ";
  int l;

  PRINTFD(FB_ObjectMolecule)
    " findflag: flag %s format %s\n",flag,format
    ENDFD;

  strcat(pat,flag);
  l=strlen(pat);
  while(*p) {
    p=ncopy(cc,p,l);
    if(WordMatch(cc,pat,true)<0) {
      p=nextline(p);
      break;
    }
    p=nextline(p);
    if(!*p) {
      PRINTFB(FB_ObjectMolecule,FB_Errors)
        " ObjectMolecule-Error: Unrecognized file format (can't find \"%s\").\n",pat
        ENDFB;
    }
  }

  strcpy(pat,"%FORMAT(");
  strcat(pat,format);
  strcat(pat,")");
  l=strlen(pat);
  while(*p) {
    p=ncopy(cc,p,l);
    if(WordMatch(cc,pat,true)<0) {
      p=nextline(p);
      break; 
    }
    p=nextline(p);
    if(!*p) {
      PRINTFB(FB_ObjectMolecule,FB_Errors)
        " ObjectMolecule-Error: Unrecognized file format (can't find \"%s\").\n",pat
        ENDFB;
    }
      
  }
  return(p);
}

#define nextline_top nextline

/*========================================================================*/
CoordSet *ObjectMoleculeTOPStr2CoordSet(char *buffer,
                                        AtomInfoType **atInfoPtr)
{
  char *p;
  int nAtom;
  int a,b,c,bi,last_i,at_i,aa,rc;
  float *coord = NULL;
  float *f;
  CoordSet *cset = NULL;
  AtomInfoType *atInfo = NULL,*ai;
  BondType *bond=NULL,*bd;
  int nBond=0;
  int auto_show_lines = (int)SettingGet(cSetting_auto_show_lines);
  int auto_show_nonbonded = (int)SettingGet(cSetting_auto_show_nonbonded);
  int amber7 = false;

  WordType title;
  ResName *resn;

  char cc[MAXLINELEN];
  int ok=true;
  int i0,i1,i2;

  /* trajectory parameters */

  int NTYPES,NBONH,MBONA,NTHETH,MTHETA;
  int NPHIH,MPHIA,NHPARM,NPARM,NNB,NRES;
  int NBONA,NTHETA,NPHIA,NUMBND,NUMANG,NPTRA;
  int NATYP,NPHB,IFPERT,NBPER,NGPER,NDPER;
  int MBPER,MGPER,MDPER,IFBOX,NMXRS,IFCAP;
  int NEXTRA,IPOL=0;
  int wid,col;
  float BETA;
  float BOX1,BOX2,BOX3;

  cset = CoordSetNew();  

  p=buffer;
  nAtom=0;
  if(atInfoPtr)
	 atInfo = *atInfoPtr;
  if(!atInfo)
    ErrFatal("TOPStr2CoordSet","need atom information record!");
  /* failsafe for old version..*/

  ncopy(cc,p,8);
  if(strcmp(cc,"%VERSION")==0) {
    amber7=true;
    PRINTFB(FB_ObjectMolecule,FB_Details)
      " ObjectMolecule: Attempting to read Amber7 topology file.\n"
      ENDFB;
  } else {
    PRINTFB(FB_ObjectMolecule,FB_Details)
      " ObjectMolecule: Assuming this is an Amber6 topology file.\n"
      ENDFB;
  }
  
  /* read title */
  if(amber7) {
    p = findflag(buffer,"TITLE","20a4");
  }

  p=ncopy(cc,p,20);
  title[0]=0;
  sscanf(cc,"%s",title);
  p=nextline_top(p);

  if(amber7) {

    p = findflag(buffer,"POINTERS","10I8");

    p=ncopy(cc,p,8); ok = ok && (sscanf(cc,"%d",&nAtom)==1);
    p=ncopy(cc,p,8); ok = ok && (sscanf(cc,"%d",&NTYPES)==1);
    p=ncopy(cc,p,8); ok = ok && (sscanf(cc,"%d",&NBONH)==1);
    p=ncopy(cc,p,8); ok = ok && (sscanf(cc,"%d",&MBONA)==1);
    p=ncopy(cc,p,8); ok = ok && (sscanf(cc,"%d",&NTHETH)==1);

    p=ncopy(cc,p,8); ok = ok && (sscanf(cc,"%d",&MTHETA)==1);
    p=ncopy(cc,p,8); ok = ok && (sscanf(cc,"%d",&NPHIH)==1);
    p=ncopy(cc,p,8); ok = ok && (sscanf(cc,"%d",&MPHIA)==1);
    p=ncopy(cc,p,8); ok = ok && (sscanf(cc,"%d",&NHPARM)==1);
    p=ncopy(cc,p,8); ok = ok && (sscanf(cc,"%d",&NPARM)==1);

    p=nextline_top(p);

    p=ncopy(cc,p,8); ok = ok && (sscanf(cc,"%d",&NNB)==1);
    p=ncopy(cc,p,8); ok = ok && (sscanf(cc,"%d",&NRES)==1);
    p=ncopy(cc,p,8); ok = ok && (sscanf(cc,"%d",&NBONA)==1);
    p=ncopy(cc,p,8); ok = ok && (sscanf(cc,"%d",&NTHETA)==1);
    p=ncopy(cc,p,8); ok = ok && (sscanf(cc,"%d",&NPHIA)==1);

    p=ncopy(cc,p,8); ok = ok && (sscanf(cc,"%d",&NUMBND)==1);
    p=ncopy(cc,p,8); ok = ok && (sscanf(cc,"%d",&NUMANG)==1);
    p=ncopy(cc,p,8); ok = ok && (sscanf(cc,"%d",&NPTRA)==1);
    p=ncopy(cc,p,8); ok = ok && (sscanf(cc,"%d",&NATYP)==1);
    p=ncopy(cc,p,8); ok = ok && (sscanf(cc,"%d",&NPHB)==1);

    p=nextline_top(p);

    p=ncopy(cc,p,8); ok = ok && (sscanf(cc,"%d",&IFPERT)==1);
    p=ncopy(cc,p,8); ok = ok && (sscanf(cc,"%d",&NBPER)==1);
    p=ncopy(cc,p,8); ok = ok && (sscanf(cc,"%d",&NGPER)==1);
    p=ncopy(cc,p,8); ok = ok && (sscanf(cc,"%d",&NDPER)==1);
    p=ncopy(cc,p,8); ok = ok && (sscanf(cc,"%d",&MBPER)==1);

    p=ncopy(cc,p,8); ok = ok && (sscanf(cc,"%d",&MGPER)==1);
    p=ncopy(cc,p,8); ok = ok && (sscanf(cc,"%d",&MDPER)==1);
    p=ncopy(cc,p,8); ok = ok && (sscanf(cc,"%d",&IFBOX)==1);
    p=ncopy(cc,p,8); ok = ok && (sscanf(cc,"%d",&NMXRS)==1);
    p=ncopy(cc,p,8); ok = ok && (sscanf(cc,"%d",&IFCAP)==1);
    
    p=nextline_top(p);
    p=ncopy(cc,p,8); ok = ok && (sscanf(cc,"%d",&NEXTRA)==1);

  } else {
    
    p=ncopy(cc,p,6); ok = ok && (sscanf(cc,"%d",&nAtom)==1);
    p=ncopy(cc,p,6); ok = ok && (sscanf(cc,"%d",&NTYPES)==1);
    p=ncopy(cc,p,6); ok = ok && (sscanf(cc,"%d",&NBONH)==1);
    p=ncopy(cc,p,6); ok = ok && (sscanf(cc,"%d",&MBONA)==1);
    p=ncopy(cc,p,6); ok = ok && (sscanf(cc,"%d",&NTHETH)==1);
    p=ncopy(cc,p,6); ok = ok && (sscanf(cc,"%d",&MTHETA)==1);
    
    p=ncopy(cc,p,6); ok = ok && (sscanf(cc,"%d",&NPHIH)==1);
    p=ncopy(cc,p,6); ok = ok && (sscanf(cc,"%d",&MPHIA)==1);
    p=ncopy(cc,p,6); ok = ok && (sscanf(cc,"%d",&NHPARM)==1);
    p=ncopy(cc,p,6); ok = ok && (sscanf(cc,"%d",&NPARM)==1);
    p=ncopy(cc,p,6); ok = ok && (sscanf(cc,"%d",&NNB)==1);
    p=ncopy(cc,p,6); ok = ok && (sscanf(cc,"%d",&NRES)==1);
    
    p=nextline_top(p);
    
    p=ncopy(cc,p,6); ok = ok && (sscanf(cc,"%d",&NBONA)==1);
    p=ncopy(cc,p,6); ok = ok && (sscanf(cc,"%d",&NTHETA)==1);
    p=ncopy(cc,p,6); ok = ok && (sscanf(cc,"%d",&NPHIA)==1);
    p=ncopy(cc,p,6); ok = ok && (sscanf(cc,"%d",&NUMBND)==1);
    p=ncopy(cc,p,6); ok = ok && (sscanf(cc,"%d",&NUMANG)==1);
    p=ncopy(cc,p,6); ok = ok && (sscanf(cc,"%d",&NPTRA)==1);
    
    p=ncopy(cc,p,6); ok = ok && (sscanf(cc,"%d",&NATYP)==1);
    p=ncopy(cc,p,6); ok = ok && (sscanf(cc,"%d",&NPHB)==1);
    p=ncopy(cc,p,6); ok = ok && (sscanf(cc,"%d",&IFPERT)==1);
    p=ncopy(cc,p,6); ok = ok && (sscanf(cc,"%d",&NBPER)==1);
    p=ncopy(cc,p,6); ok = ok && (sscanf(cc,"%d",&NGPER)==1);
    p=ncopy(cc,p,6); ok = ok && (sscanf(cc,"%d",&NDPER)==1);
    
    p=nextline_top(p);
    
    p=ncopy(cc,p,6); ok = ok && (sscanf(cc,"%d",&MBPER)==1);
    p=ncopy(cc,p,6); ok = ok && (sscanf(cc,"%d",&MGPER)==1);
    p=ncopy(cc,p,6); ok = ok && (sscanf(cc,"%d",&MDPER)==1);
    p=ncopy(cc,p,6); ok = ok && (sscanf(cc,"%d",&IFBOX)==1);
    p=ncopy(cc,p,6); ok = ok && (sscanf(cc,"%d",&NMXRS)==1);
    p=ncopy(cc,p,6); ok = ok && (sscanf(cc,"%d",&IFCAP)==1);
    
    p=ncopy(cc,p,6); if(sscanf(cc,"%d",&NEXTRA)!=1) NEXTRA=0;

  }

  switch(IFBOX) {
  case 2:
    cset->PeriodicBoxType = cCSet_Octahedral;
    PRINTFB(FB_ObjectMolecule,FB_Details)
      " TOPStrToCoordSet: Warning: can't currently image a truncated octahedron...\n"
      ENDFB;
    break;
  case 1:
    cset->PeriodicBoxType = cCSet_Orthogonal;
    break;
  case 0:
  default:
    cset->PeriodicBoxType = cCSet_NoPeriodicity;
    break;
  }

  p=nextline_top(p);

  if(!ok) {
    ErrMessage("TOPStrToCoordSet","Error reading counts lines");
  } else {
    PRINTFB(FB_ObjectMolecule,FB_Blather)
      " TOPStr2CoordSet: read counts line nAtom %d NBONA %d NBONH %d\n",
      nAtom,NBONA,NBONH
      ENDFB;
  }

  if(ok) {  
    VLACheck(atInfo,AtomInfoType,nAtom);

    if(amber7) {
      p = findflag(buffer,"ATOM_NAME","20a4");
    }
    /* read atoms */

    b=0;
    for(a=0;a<nAtom;a++) {
      p=ncopy(cc,p,4);
      ai=atInfo+a;
      if(!sscanf(cc,"%s",ai->name))
        ai->name[0]=0;
      if((++b)==20) {
        b=0;
        p=nextline_top(p);
      }
    }
  
    if(b) p=nextline_top(p);

    if(!ok) {
      ErrMessage("TOPStrToCoordSet","Error reading atom names");
    } else {
      PRINTFB(FB_ObjectMolecule,FB_Blather)
        " TOPStr2CoordSet: read atom names.\n"
        ENDFB;
    }

    /* read charges */

    if(amber7) {
      p = findflag(buffer,"CHARGE","5E16.8");
    }

    b=0;
    for(a=0;a<nAtom;a++) {
      p=ncopy(cc,p,16);
      ai=atInfo+a;
      if(!sscanf(cc,"%f",&ai->partialCharge))
        ok=false;
      else {
        ai->partialCharge/=18.2223F; /* convert to electron charge */
      }
      if((++b)==5) {
        b=0;
        p=nextline_top(p);
      }
    }

    if(!ok) {
      ErrMessage("TOPStrToCoordSet","Error reading charges");
    } else {
      PRINTFB(FB_ObjectMolecule,FB_Blather)
        " TOPStr2CoordSet: read charges.\n"
        ENDFB;
    }
    if(b) p=nextline_top(p);
  
    if(!amber7) {
      /* skip masses */
      
      p=skip_fortran(nAtom,5,p);
    }

    /* read LJ atom types */

    if(amber7) {
      p = findflag(buffer,"ATOM_TYPE_INDEX","10I8");
      col=10;
      wid=8;
    } else {
      col=12;
      wid=6;
    }

    b=0;
    for(a=0;a<nAtom;a++) {
      p=ncopy(cc,p,wid);
      ai=atInfo+a;
      if(!sscanf(cc,"%d",&ai->customType))
        ok=false;
      if((++b)==col) {
        b=0;
        p=nextline_top(p);
      }
    }
    if(b) p=nextline_top(p);

    if(!ok) {
      ErrMessage("TOPStrToCoordSet","Error LJ atom types");
    } else {
      PRINTFB(FB_ObjectMolecule,FB_Blather)
        " TOPStr2CoordSet: read LJ atom types.\n"
        ENDFB;
    }

    if(!amber7) {
      /* skip excluded atom counts */
      
      p=skip_fortran(nAtom,12,p);
      
      /* skip NB param arrays */
      
      p=skip_fortran(NTYPES*NTYPES,12,p);
    }

    /* read residue labels */

    if(amber7) {
      p = findflag(buffer,"RESIDUE_LABEL","20a4");
    }

    resn = Alloc(ResName,NRES);

    b=0;
    for(a=0;a<NRES;a++) {
      p=ncopy(cc,p,4);
      if(!sscanf(cc,"%s",resn[a]))
        resn[a][0]=0;
      if((++b)==20) {
        b=0;
        p=nextline_top(p);
      }
    }
    if(b) p=nextline_top(p);

    if(!ok) {
      ErrMessage("TOPStrToCoordSet","Error reading residue labels");
    } else {
      PRINTFB(FB_ObjectMolecule,FB_Blather)
        " TOPStr2CoordSet: read residue labels.\n"
        ENDFB;
    }

    /* read residue assignments */

    if(amber7) {
      p = findflag(buffer,"RESIDUE_POINTER","10I8");
      col=10;
      wid=8;
    } else {
      col=12;
      wid=6;
    }

    b=0;
    last_i=0;
    rc=0;
    for(a=0;a<NRES;a++) {
      p=ncopy(cc,p,wid);
      if(sscanf(cc,"%d",&at_i))
        {
          if(last_i)
            for(aa=(last_i-1);aa<(at_i-1);aa++) {
              ai = atInfo+aa;
              strcpy(ai->resn,resn[a-1]);
              ai->resv=rc;
              sprintf(ai->resi,"%d",rc);
            }
          rc++;
          last_i=at_i;
        }
      if((++b)==col) {
        b=0;
        p=nextline_top(p);
      }
    }
    if(b) p=nextline_top(p);
    if(last_i)
      for(aa=(last_i-1);aa<nAtom;aa++) {
        ai = atInfo+aa;
        strcpy(ai->resn,resn[NRES-1]);
        ai->resv=rc;
        sprintf(ai->resi,"%d",rc);
      }
    rc++;

    if(!ok) {
      ErrMessage("TOPStrToCoordSet","Error reading residues");
    } else {
      PRINTFB(FB_ObjectMolecule,FB_Blather)
        " TOPStr2CoordSet: read residues.\n"
        ENDFB;
    }

    FreeP(resn);

    if(!amber7) {  
      /* skip bond force constants */

      p=skip_fortran(NUMBND,5,p);
      
      /* skip bond lengths */

      p=skip_fortran(NUMBND,5,p);
      
      /* skip angle force constant */
      
      p=skip_fortran(NUMANG,5,p);
      
      /* skip angle eq */
      
      p=skip_fortran(NUMANG,5,p);
      
      /* skip dihedral force constant */
      
      p=skip_fortran(NPTRA,5,p);
      
      /* skip dihedral periodicity */
      
      p=skip_fortran(NPTRA,5,p);
      
      /* skip dihedral phases */
      
      p=skip_fortran(NPTRA,5,p);
      
      /* skip SOLTYs */
      
      p=skip_fortran(NATYP,5,p);
      
      /* skip LJ terms r12 */
      
      p=skip_fortran((NTYPES*(NTYPES+1))/2,5,p);
      
      /* skip LJ terms r6 */
      
      p=skip_fortran((NTYPES*(NTYPES+1))/2,5,p);
      
    }

    /* read bonds */

    if(amber7) {
      p = findflag(buffer,"BONDS_INC_HYDROGEN","10I8");
      col=10;
      wid=8;
    } else {
      col=12;
      wid=6;
    }
    
    nBond = NBONH + NBONA;

    bond=VLAlloc(BondType,nBond);
  
    bi = 0;
  

    b=0;
    c=0;
    i0=0;
    i1=0;
    for(a=0;a<3*NBONH;a++) {
      p=ncopy(cc,p,wid);
      i2=i1;
      i1=i0;
      if(!sscanf(cc,"%d",&i0))
        ok=false;
      if((++c)==3) {
        c=0;
        bd=bond+bi;
        bd->index[0]=(abs(i2)/3);
        bd->index[1]=(abs(i1)/3);
        bd->order=1;
        bd->stereo=0;
        bd->id = bi+1;
        bi++;
      }
      if((++b)==col) {
        b=0;
        p=nextline_top(p);
      }
    }
    if(b) p=nextline_top(p);

    if(!ok) {
      ErrMessage("TOPStrToCoordSet","Error hydrogen containing bonds");
    } else {
      PRINTFB(FB_ObjectMolecule,FB_Blather)
        " TOPStr2CoordSet: read %d hydrogen containing bonds.\n",NBONH
        ENDFB;
    }

    if(amber7) {
      p = findflag(buffer,"BONDS_WITHOUT_HYDROGEN","10I8");
      col=10;
      wid=8;
    } else {
      col=12;
      wid=6;
    }

    b=0;
    c=0;
    for(a=0;a<3*NBONA;a++) {
      p=ncopy(cc,p,wid);
      i2=i1;
      i1=i0;
      if(!sscanf(cc,"%d",&i0))
        ok=false;
      if((++c)==3) {
        c=0;
        bd=bond+bi;
        bd->index[0]=(abs(i2)/3);
        bd->index[1]=(abs(i1)/3);
        bd->order=0;
        bd->stereo=0;
        bd->id = bi+1;
        bi++;
      }
      if((++b)==col) {
        b=0;
        p=nextline_top(p);
      }
    }
    if(b) p=nextline_top(p);

    if(!ok) {
      ErrMessage("TOPStrToCoordSet","Error hydrogen free bonds");
    } else {
      PRINTFB(FB_ObjectMolecule,FB_Blather)
        " TOPStr2CoordSet: read %d hydrogen free bonds.\n",NBONA
        ENDFB;
    }

    if(!amber7) {
      /* skip hydrogen angles */
      
      p=skip_fortran(4*NTHETH,12,p);
      
      /* skip non-hydrogen angles */
      
      p=skip_fortran(4*NTHETA,12,p);
      
      /* skip hydrogen dihedrals */
      
      p=skip_fortran(5*NPHIH,12,p);
      
      /* skip non hydrogen dihedrals */
      
      p=skip_fortran(5*NPHIA,12,p);
      
      /* skip nonbonded exclusions */
      
      p=skip_fortran(NNB,12,p);
      
      /* skip hydrogen bonds ASOL */
      
      p=skip_fortran(NPHB,5,p);
      
      /* skip hydrogen bonds BSOL */
      
      p=skip_fortran(NPHB,5,p);
      
      /* skip HBCUT */
      
      p=skip_fortran(NPHB,5,p);
      
    }
    /* read AMBER atom types */

    if(amber7) {
      p = findflag(buffer,"AMBER_ATOM_TYPE","20a4");
    }

    b=0;
    for(a=0;a<nAtom;a++) {
      p=ncopy(cc,p,4);
      ai=atInfo+a;
      if(!sscanf(cc,"%s",ai->textType))
        ok=false;
      if((++b)==20) {
        b=0;
        p=nextline_top(p);
      }
    }
    if(b) p=nextline_top(p);

    if(!ok) {
      ErrMessage("TOPStrToCoordSet","Error reading atom types");
    } else {
      PRINTFB(FB_ObjectMolecule,FB_Blather)
        " TOPStr2CoordSet: read atom types.\n"
        ENDFB;
    }

    if(!amber7) {
      /* skip TREE classification */
      
      p=skip_fortran(nAtom,20,p);
      
      /* skip tree joining information */
      
      p=skip_fortran(nAtom,12,p);
      
      /* skip last atom rotated blah blah blah */
      
      p=skip_fortran(nAtom,12,p);

    }

    if(IFBOX>0) {
      
      int IPTRES,NSPM,NSPSOL;
      
      if(amber7) {
        p = findflag(buffer,"SOLVENT_POINTERS","3I8");
        wid=8;
      } else {
        wid=6;
      }
      p=ncopy(cc,p,wid); ok = ok && (sscanf(cc,"%d",&IPTRES)==1);
      p=ncopy(cc,p,wid); ok = ok && (sscanf(cc,"%d",&NSPM)==1);
      p=ncopy(cc,p,wid); ok = ok && (sscanf(cc,"%d",&NSPSOL)==1);
      
      p=nextline_top(p);
      
      if(amber7) {
        p = findflag(buffer,"ATOMS_PER_MOLECULE","10I8");
        col=10;
      } else {
        col=12;
      }
      
      /* skip num atoms per box */
      
      p=skip_fortran(NSPM,col,p);
      
      if(amber7) {
        p = findflag(buffer,"BOX_DIMENSIONS","5E16.8");
      }
      wid=16;
      
      p=ncopy(cc,p,16); ok = ok && (sscanf(cc,"%f",&BETA)==1);
      p=ncopy(cc,p,16); ok = ok && (sscanf(cc,"%f",&BOX1)==1);
      p=ncopy(cc,p,16); ok = ok && (sscanf(cc,"%f",&BOX2)==1);
      p=ncopy(cc,p,16); ok = ok && (sscanf(cc,"%f",&BOX3)==1);
      
      if(ok) {
        if(!cset->PeriodicBox) 
          cset->PeriodicBox=CrystalNew();
        cset->PeriodicBox->Dim[0] = BOX1;
        cset->PeriodicBox->Dim[1] = BOX2;
        cset->PeriodicBox->Dim[2] = BOX3;
        if((BETA > 109.47) && (BETA < 109.48)) {
          cset->PeriodicBoxType=cCSet_Octahedral;
          cset->PeriodicBox->Angle[0]=(float)(2.0*acos(1.0/sqrt(3.0))*180.0/PI);
          cset->PeriodicBox->Angle[1]=(float)(2.0*acos(1.0/sqrt(3.0))*180.0/PI);
          cset->PeriodicBox->Angle[2]=(float)(2.0*acos(1.0/sqrt(3.0))*180.0/PI);
        } else if(BETA==60.0) {
          cset->PeriodicBox->Angle[0]=60.0; /* rhombic dodecahedron (from ptraj.c) */
          cset->PeriodicBox->Angle[1]=90.0;
          cset->PeriodicBox->Angle[2]=60.0;
        } else {
          cset->PeriodicBox->Angle[0]=90.0;
          cset->PeriodicBox->Angle[1]=BETA;
          cset->PeriodicBox->Angle[2]=90.0;
        }
        
        CrystalUpdate(cset->PeriodicBox);
        /*        CrystalDump(cset->PeriodicBox);*/
      }
      /* skip periodic box */
      
      p=nextline_top(p);
      
    }
    
    if(!amber7) {

      if(IFCAP>0) {
        p=nextline_top(p);
        p=nextline_top(p);
        p=nextline_top(p);
      }

      if(IFPERT>0) {

        /* skip perturbed bond atoms */

        p=skip_fortran(2*NBPER,12,p);    

        /* skip perturbed bond atom pointers */

        p=skip_fortran(2*NBPER,12,p);    

        /* skip perturbed angles */

        p=skip_fortran(3*NGPER,12,p);    

        /* skip perturbed angle pointers */

        p=skip_fortran(2*NGPER,12,p);    

        /* skip perturbed dihedrals */

        p=skip_fortran(4*NDPER,12,p);    

        /* skip perturbed dihedral pointers */

        p=skip_fortran(2*NDPER,12,p);    

        /* skip residue names */

        p=skip_fortran(NRES,20,p);    

        /* skip atom names */

        p=skip_fortran(nAtom,20,p);    

        /* skip atom symbols */

        p=skip_fortran(nAtom,20,p);    

        /* skip unused field */

        p=skip_fortran(nAtom,5,p);    

        /* skip perturbed flags */

        p=skip_fortran(nAtom,12,p);    

        /* skip LJ atom flags */

        p=skip_fortran(nAtom,12,p);    

        /* skip perturbed charges */

        p=skip_fortran(nAtom,5,p);    

      }

      if(IPOL>0) {

        /* skip atomic polarizabilities */

        p=skip_fortran(nAtom,5,p);    

      }

      if((IPOL>0) && (IFPERT>0)) {

        /* skip atomic polarizabilities */

        p=skip_fortran(nAtom,5,p);    
    
      }
    }
    /* for future reference 

%FLAG LES_NTYP
%FORMAT(10I8)
%FLAG LES_TYPE
%FORMAT(10I8)
%FLAG LES_FAC
%FORMAT(5E16.8)
%FLAG LES_CNUM
%FORMAT(10I8)
%FLAG LES_ID
%FORMAT(10I8)

Here is the additional information for LES topology formats:
First, if NPARM ==1, LES entries are in topology (NPARM is the 10th
entry in the initial list of control parameters); otherwise the standard
format applies.
So, with NPARM=1, you just need to read a few more things at the very
end of topology file:
LES_NTYP (format: I6) ... one number, number of LES types
and four arrays:
LES_TYPE (12I6) ... NATOM integer entries
LES_FAC (E16.8) ... LES_NTYPxLES_NTYP float entries
LES_CNUM (12I6) ... NATOM integer entries
LES_ID (12I6)   ... NATOM integer entries

and that's it. Your parser must have skipped this information because it
was at the end of the file. Maybe that's good enough.



     */

    coord=VLAlloc(float,3*nAtom);

    f=coord;
    for(a=0;a<nAtom;a++) {
      *(f++)=0.0;
      *(f++)=0.0;
      *(f++)=0.0;
      ai = atInfo + a;
      ai->id = a+1; /* assign 1-based identifiers */
      AtomInfoAssignParameters(ai);
      ai->color=AtomInfoGetColor(ai);
      for(c=0;c<cRepCnt;c++) {
        ai->visRep[c] = false;
      }
      ai->visRep[cRepLine] = auto_show_lines; /* show lines by default */
      ai->visRep[cRepNonbonded] = auto_show_nonbonded; /* show lines by default */
    }
  }
  if(ok) {
    cset->NIndex=nAtom;
    cset->Coord=coord;
    cset->TmpBond=bond;
    cset->NTmpBond=nBond;
  } else {
    if(cset) 
      cset->fFree(cset);
  }
  if(atInfoPtr)
	 *atInfoPtr = atInfo;
  
  return(cset);
}

/*========================================================================*/
ObjectMolecule *ObjectMoleculeReadTOPStr(ObjectMolecule *I,char *TOPStr,int frame,int discrete)
{
  CoordSet *cset = NULL;
  AtomInfoType *atInfo;
  int ok=true;
  int isNew = true;
  unsigned int nAtom = 0;

  if(!I) 
	 isNew=true;
  else 
	 isNew=false;

  if(ok) {

	 if(isNew) {
		I=(ObjectMolecule*)ObjectMoleculeNew(discrete);
		atInfo = I->AtomInfo;
		isNew = true;
	 } else { /* never */
		atInfo=VLAMalloc(10,sizeof(AtomInfoType),2,true); /* autozero here is important */
		isNew = false;
	 }
    if(isNew) {
      AtomInfoPrimeColors();
      I->Obj.Color = AtomInfoGetCarbColor();
    }

	 cset=ObjectMoleculeTOPStr2CoordSet(TOPStr,&atInfo);	 
	 nAtom=cset->NIndex;
  }

  /* include coordinate set */
  if(ok) {
    cset->Obj = I;
    cset->fEnumIndices(cset);
    if(cset->fInvalidateRep)
      cset->fInvalidateRep(cset,cRepAll,cRepInvRep);
    if(isNew) {		
      I->AtomInfo=atInfo; /* IMPORTANT to reassign: this VLA may have moved! */
    } else {
      ObjectMoleculeMerge(I,atInfo,cset,false,cAIC_AllMask); /* NOTE: will release atInfo */
    }
    if(isNew) I->NAtom=nAtom;
    /* 
       if(frame<0) frame=I->NCSet;
       VLACheck(I->CSet,CoordSet*,frame);
       if(I->NCSet<=frame) I->NCSet=frame+1;
       if(I->CSet[frame]) I->CSet[frame]->fFree(I->CSet[frame]);
       I->CSet[frame] = cset;
    */

    if(isNew) I->NBond = ObjectMoleculeConnect(I,&I->Bond,I->AtomInfo,cset,false);
    if(cset->Symmetry&&(!I->Symmetry)) {
      I->Symmetry=SymmetryCopy(cset->Symmetry);
      SymmetryAttemptGeneration(I->Symmetry,false,false);
    }

    if(I->CSTmpl)
      if(I->CSTmpl->fFree)
        I->CSTmpl->fFree(I->CSTmpl);
    I->CSTmpl = cset; /* save template coordinate set */

    SceneCountFrames();
    ObjectMoleculeExtendIndices(I);
    ObjectMoleculeSort(I);
    ObjectMoleculeUpdateIDNumbers(I);
    ObjectMoleculeUpdateNonbonded(I);
  }
  return(I);
}

ObjectMolecule *ObjectMoleculeLoadTOPFile(ObjectMolecule *obj,char *fname,int frame,int discrete)
{
  ObjectMolecule *I=NULL;
  int ok=true;
  FILE *f;
  long size;
  char *buffer,*p;

  f=fopen(fname,"rb");
  if(!f)
	 ok=ErrMessage("ObjectMoleculeLoadTOPFile","Unable to open file!");
  else
	 {
      PRINTFB(FB_ObjectMolecule,FB_Blather) 
        " ObjectMoleculeLoadTOPFile: Loading from %s.\n",fname
        ENDFB;
		
		fseek(f,0,SEEK_END);
      size=ftell(f);
		fseek(f,0,SEEK_SET);

		buffer=(char*)mmalloc(size+255);
		ErrChkPtr(buffer);
		p=buffer;
		fseek(f,0,SEEK_SET);
		fread(p,size,1,f);
		p[size]=0;
		fclose(f);

		I=ObjectMoleculeReadTOPStr(obj,buffer,frame,discrete);

		mfree(buffer);
	 }

  return(I);
}

void ObjectMoleculeSculptClear(ObjectMolecule *I)
{
  PRINTFD(FB_ObjectMolecule)
    " ObjectMoleculeSculptClear: entered.\n"
    ENDFD;

  if(I->Sculpt) SculptFree(I->Sculpt);
  I->Sculpt=NULL;
}

void ObjectMoleculeSculptImprint(ObjectMolecule *I,int state)
{
  PRINTFD(FB_ObjectMolecule)
    " ObjectMoleculeUpdateSculpt: entered.\n"
    ENDFD;

  if(!I->Sculpt) I->Sculpt = SculptNew();
  SculptMeasureObject(I->Sculpt,I,state);
}

float ObjectMoleculeSculptIterate(ObjectMolecule *I,int state,int n_cycle)
{
  PRINTFD(FB_ObjectMolecule)
    " ObjectMoleculeIterateSculpt: entered.\n"
    ENDFD;
  if(I->Sculpt) {
    return SculptIterateObject(I->Sculpt,I,state,n_cycle);
  } else
    return 0.0F;
}

void ObjectMoleculeUpdateIDNumbers(ObjectMolecule *I)
{
  int a;
  int max;
  AtomInfoType *ai;
  BondType *b;

  if(I->AtomCounter<0) {
    max=-1;
    ai=I->AtomInfo;
    for(a=0;a<I->NAtom;a++) {
      if(ai->id>max)
        max=ai->id;
      ai++;
    }
    I->AtomCounter=max+1;
  }
  ai=I->AtomInfo;
  for(a=0;a<I->NAtom;a++) {
    if(ai->id<0) 
      ai->id=I->AtomCounter++;
    ai++;
  }

  if(I->BondCounter<0) {
    max=-1;
    b=I->Bond;
    for(a=0;a<I->NBond;a++) {
      if(b->id>max) 
        max=b->id;
      b++;
    }
    I->BondCounter=max+1;
  }
  b=I->Bond;
  for(a=0;a<I->NBond;a++) {
    if(!b->id) 
      b->id=I->BondCounter++;
    b++;
  }
}

CoordSet *ObjectMoleculePMO2CoordSet(CRaw *pmo,AtomInfoType **atInfoPtr,int *restart)
{
  int nAtom,nBond;
  int a;
  float *coord = NULL;
  CoordSet *cset = NULL;
  AtomInfoType *atInfo = NULL,*ai;
  AtomInfoType068 *atInfo068 = NULL;
  AtomInfoType076 *atInfo076 = NULL;
  AtomInfoType083 *atInfo083 = NULL;
  BondType *bond=NULL;
  BondType068 *bond068=NULL;
  BondType083 *bond083=NULL;

  int ok=true;
  int auto_show_lines;
  int auto_show_nonbonded;
  int type,size;
  float *spheroid=NULL;
  float *spheroid_normal=NULL;
  int sph_info[2];
  int version;
  auto_show_lines = (int)SettingGet(cSetting_auto_show_lines);
  auto_show_nonbonded = (int)SettingGet(cSetting_auto_show_nonbonded);

  *restart=false;
  nAtom=0;
  nBond=0;
  if(atInfoPtr)
	 atInfo = *atInfoPtr;
  
  type = RawGetNext(pmo,&size,&version);
  if(type!=cRaw_AtomInfo1) {
    ok=false;
  } else { /* read atoms */
    PRINTFD(FB_ObjectMolecule)
      " ObjectMolPMO2CoordSet: loading atom info %d bytes = %8.3f\n",size,((float)size)/sizeof(AtomInfoType)
      ENDFD;
    if(version<66) {
      PRINTFB(FB_ObjectMolecule,FB_Errors)
        " ObjectMolecule: unsupported binary file (version %d). aborting.\n",
        version
        ENDFB;
      ok=false;
    } else if(version<69) { /* legacy atom format */
      nAtom = size/sizeof(AtomInfoType068);
      atInfo068 = Alloc(AtomInfoType068,nAtom);
      ok = RawReadInto(pmo,cRaw_AtomInfo1,size,(char*)atInfo068);
      VLACheck(atInfo,AtomInfoType,nAtom);
      UtilExpandArrayElements(atInfo068,atInfo,nAtom,
                              sizeof(AtomInfoType068),sizeof(AtomInfoType));
      FreeP(atInfo068);
    } else if(version<77) { /* legacy atom format */
      nAtom = size/sizeof(AtomInfoType076);
      atInfo076 = Alloc(AtomInfoType076,nAtom);
      ok = RawReadInto(pmo,cRaw_AtomInfo1,size,(char*)atInfo076);
      VLACheck(atInfo,AtomInfoType,nAtom);
      UtilExpandArrayElements(atInfo076,atInfo,nAtom,
                              sizeof(AtomInfoType076),sizeof(AtomInfoType));
      FreeP(atInfo076);
      
    } else if(version<84) { /* legacy atom format */
      nAtom = size/sizeof(AtomInfoType083);
      atInfo083 = Alloc(AtomInfoType083,nAtom);
      ok = RawReadInto(pmo,cRaw_AtomInfo1,size,(char*)atInfo083);
      VLACheck(atInfo,AtomInfoType,nAtom);
      UtilExpandArrayElements(atInfo083,atInfo,nAtom,
                              sizeof(AtomInfoType083),sizeof(AtomInfoType));
      FreeP(atInfo083);
      
    } else {
      nAtom = size/sizeof(AtomInfoType);
      VLACheck(atInfo,AtomInfoType,nAtom);
      ok = RawReadInto(pmo,cRaw_AtomInfo1,size,(char*)atInfo);
    }
  }
  if(ok) {
    PRINTFD(FB_ObjectMolecule)
      " ObjectMolPMO2CoordSet: loading coordinates\n"
      ENDFD;
    coord = (float*)RawReadVLA(pmo,cRaw_Coords1,sizeof(float),5,false);
    if(!coord)
      ok=false;
  }
  type = RawGetNext(pmo,&size,&version);
  if(type==cRaw_SpheroidInfo1) {

    PRINTFD(FB_ObjectMolecule)
      " ObjectMolPMO2CoordSet: loading spheroid\n"
      ENDFD;

    ok = RawReadInto(pmo,cRaw_SpheroidInfo1,sizeof(int)*2,(char*)sph_info);
    if(ok) {

    PRINTFD(FB_ObjectMolecule)
      " ObjectMolPMO2CoordSet: loading spheroid size %d nsph %d\n",sph_info[0],sph_info[1]
      ENDFD;

      spheroid = (float*)RawReadPtr(pmo,cRaw_Spheroid1,&size);
      if(!spheroid)
        ok=false;

      PRINTFD(FB_ObjectMolecule)
        " ObjectMolPMO2CoordSet: loaded spheroid %p size %d \n",spheroid,size
        ENDFD;

    }
    if(ok) {
      spheroid_normal = (float*)RawReadPtr(pmo,cRaw_SpheroidNormals1,&size);
      if(!spheroid_normal)
        ok=false;
      }
      PRINTFD(FB_ObjectMolecule)
        " ObjectMolPMO2CoordSet: loaded spheroid %p size %d \n",spheroid_normal,size
        ENDFD;

    } 
  if(ok) 
      type = RawGetNext(pmo,&size,&version);    
  if(ok) {
    
    PRINTFD(FB_ObjectMolecule)
      " ObjectMolPMO2CoordSet: loading bonds\n"
      ENDFD;

    if(type!=cRaw_Bonds1) {
      ok=false;
    } else {
      if(ok) {

        /* legacy bond format */
        if(version<66) {
          PRINTFB(FB_ObjectMolecule,FB_Errors)
            " ObjectMolecule: unsupported binary file (version %d). aborting.\n",
            version
            ENDFB;
          ok=false;
        } else if(version<69) { /* legacy atom format */
          nBond = size/sizeof(BondType068);
          bond068 = Alloc(BondType068,nBond);
          ok = RawReadInto(pmo,cRaw_Bonds1,nBond*sizeof(BondType068),(char*)bond068);
          bond=VLAlloc(BondType,nBond);
          UtilExpandArrayElements(bond068,bond,nBond,
                                  sizeof(BondType068),sizeof(BondType));
          FreeP(bond068);
          for(a=0;a<nBond;a++) bond[a].id=-1; /* initialize identifiers */
        } else if(version<84) {
          nBond = size/sizeof(BondType083);
          bond083 = Alloc(BondType083,nBond);
          ok = RawReadInto(pmo,cRaw_Bonds1,nBond*sizeof(BondType083),(char*)bond083);

          bond=VLAlloc(BondType,nBond);
          UtilExpandArrayElements(bond083,bond,nBond,
                                  sizeof(BondType083),sizeof(BondType));
          FreeP(bond083);
        } else {
          bond=(BondType*)RawReadVLA(pmo,cRaw_Bonds1,sizeof(BondType),5,false);
          nBond = VLAGetSize(bond);
        }
        
        PRINTFD(FB_ObjectMolecule)
          " ObjectMolPMO2CoordSet: found %d bonds\n",nBond
          ENDFD;

        if(Feedback(FB_ObjectMolecule,FB_Debugging)) {
          for(a=0;a<nBond;a++)
            printf(" ObjectMoleculeConnect: bond %d ind0 %d ind1 %d order %d\n",
                   a,bond[a].index[0],bond[a].index[1],bond[a].order);
        }
        
      }
    }
  }

  if(ok) {
    ai=atInfo;
    for(a=0;a<nAtom;a++) {
      ai->selEntry=0;
      ai++;
    }
	 cset = CoordSetNew();
	 cset->NIndex=nAtom;
	 cset->Coord=coord;
	 cset->NTmpBond=nBond;
	 cset->TmpBond=bond;
    if(spheroid) {
      cset->Spheroid=spheroid;
      cset->SpheroidNormal=spheroid_normal;
      cset->SpheroidSphereSize=sph_info[0];
      cset->NSpheroid = sph_info[1];

    }
  } else {
	 VLAFreeP(bond);
	 VLAFreeP(coord);
    FreeP(spheroid);
    FreeP(spheroid_normal);
  }
  if(atInfoPtr)
	 *atInfoPtr = atInfo;
  if(ok) {
    type = RawGetNext(pmo,&size,&version);
    if(type==cRaw_AtomInfo1)
      *restart=true;
  }
  return(cset);
}
/*========================================================================*/
ObjectMolecule *ObjectMoleculeReadPMO(ObjectMolecule *I,CRaw *pmo,int frame,int discrete)
{

  CoordSet *cset = NULL;
  AtomInfoType *atInfo;
  int ok=true;
  int isNew = true;
  unsigned int nAtom = 0;
  int restart =false;
  int repeatFlag = true;
  int successCnt = 0;

  while(repeatFlag) {
    repeatFlag = false;
  
    if(!I) 
      isNew=true;
    else 
      isNew=false;
    
    if(ok) {
      
      if(isNew) {
        I=(ObjectMolecule*)ObjectMoleculeNew(discrete);
        atInfo = I->AtomInfo;
        isNew = true;
      } else {
        atInfo=VLAMalloc(10,sizeof(AtomInfoType),2,true); /* autozero here is important */
        isNew = false;
      }
      if(isNew) {
        AtomInfoPrimeColors();
        I->Obj.Color = AtomInfoGetCarbColor();
      }

      cset = ObjectMoleculePMO2CoordSet(pmo,&atInfo,&restart);

      if(isNew) {		
        I->AtomInfo=atInfo; /* IMPORTANT to reassign: this VLA may have moved! */
      }
      if(cset) 
        nAtom=cset->NIndex;
      else
        ok=false;
    }
    
    /* include coordinate set */
    if(ok) {
      cset->Obj = I;
      cset->fEnumIndices(cset);
      if(cset->fInvalidateRep)
        cset->fInvalidateRep(cset,cRepAll,cRepInvRep);
      if(!isNew) {		
        ObjectMoleculeMerge(I,atInfo,cset,true,cAIC_AllMask); /* NOTE: will release atInfo */
      }
      if(isNew) I->NAtom=nAtom;
      if(frame<0) frame=I->NCSet;
      VLACheck(I->CSet,CoordSet*,frame);
      if(I->NCSet<=frame) I->NCSet=frame+1;
      if(I->CSet[frame]) I->CSet[frame]->fFree(I->CSet[frame]);
      I->CSet[frame] = cset;
      if(isNew) I->NBond = ObjectMoleculeConnect(I,&I->Bond,I->AtomInfo,cset,false);

      if(cset->Symmetry&&(!I->Symmetry)) {
        I->Symmetry=SymmetryCopy(cset->Symmetry);
        SymmetryAttemptGeneration(I->Symmetry,false,false);
      }
      SceneCountFrames();
      ObjectMoleculeExtendIndices(I);
      ObjectMoleculeSort(I);
      ObjectMoleculeUpdateIDNumbers(I);
      ObjectMoleculeUpdateNonbonded(I);
      successCnt++;
      if(successCnt>1) {
        if(successCnt==2){
          PRINTFB(FB_ObjectMolecule,FB_Actions)
            " ObjectMolReadPMO: read model %d\n",1
            ENDFB;
            }
        PRINTFB(FB_ObjectMolecule,FB_Actions)
          " ObjectMolReadPMO: read model %d\n",successCnt
          ENDFB;
      }
    }
    if(restart) {
      repeatFlag=true;
      frame=frame+1;
      restart=false;
    }
  }
  return(I);
  
  }
/*========================================================================*/
ObjectMolecule *ObjectMoleculeLoadPMOFile(ObjectMolecule *obj,char *fname,int frame,int discrete)
{
  ObjectMolecule *I=NULL;
  int ok=true;
  CRaw *raw;
    
  raw = RawOpenRead(fname);
  if(!raw)
	 ok=ErrMessage("ObjectMoleculeLoadPMOFile","Unable to open file!");
  else
	 {
      PRINTFB(FB_ObjectMolecule,FB_Blather)
        " ObjectMoleculeLoadPMOFile: Loading from %s.\n",fname
        ENDFB;
		
		I=ObjectMoleculeReadPMO(obj,raw,frame,discrete);
      RawFree(raw);
	 }
  
  return(I);
}


/*========================================================================*/
int ObjectMoleculeMultiSave(ObjectMolecule *I,char *fname,int state,int append)
{
  /* version 1 writes atominfo, coords, spheroid, bonds */
  CRaw *raw = NULL;
  int ok=true;
  int a,c,a1,a2,b1,b2;
  BondType *b;
  CoordSet *cs;
  BondType *bondVLA = NULL;
  AtomInfoType *aiVLA = NULL;
  int start,stop;
  int nBond;
  int sph_info[2];
  PRINTFD(FB_ObjectMolecule)
    " ObjectMoleculeMultiSave-Debug: entered \"%s\" state=%d\n",fname,state
    ENDFD;
    
  if(append) {
    raw = RawOpenWrite(fname);
  } else {
    raw = RawOpenAppend(fname);
  }
  if(raw) {
    aiVLA = VLAMalloc(1000,sizeof(AtomInfoType),5,true);
    bondVLA = VLAlloc(BondType,4000);
    if(state<0) {
      start=0;
      stop=I->NCSet;
    } else {
      start=state;
      if(start<0)
        start=0;
      stop=state+1;
      if(stop>I->NCSet)
        stop=I->NCSet;
    }
    for(a=start;a<stop;a++) {

      PRINTFD(FB_ObjectMolecule)
        " ObjectMMSave-Debug: state %d\n",a
        ENDFD;

      cs=I->CSet[a];
      if(cs) {
        VLACheck(aiVLA,AtomInfoType,cs->NIndex);
        nBond=0;

        /* write atoms */
        for(c=0;c<cs->NIndex;c++) {
          a1 = cs->IdxToAtm[c]; /* always valid */
          aiVLA[c]=I->AtomInfo[a1];
        }
        if(ok) ok = RawWrite(raw,cRaw_AtomInfo1,sizeof(AtomInfoType)*cs->NIndex,0,(char*)aiVLA);
        
        /* write coords */
        if(ok) ok = RawWrite(raw,cRaw_Coords1,sizeof(float)*3*cs->NIndex,0,(char*)cs->Coord);

        /* write spheroid (if one exists) */
        if(cs->Spheroid&&cs->SpheroidNormal) {
          sph_info[0]=cs->SpheroidSphereSize;
          sph_info[1]=cs->NSpheroid;
          if(ok) ok = RawWrite(raw,cRaw_SpheroidInfo1,sizeof(int)*2,0,(char*)sph_info);          
          if(ok) ok = RawWrite(raw,cRaw_Spheroid1,sizeof(float)*cs->NSpheroid,0,(char*)cs->Spheroid);          
          if(ok) ok = RawWrite(raw,cRaw_SpheroidNormals1,sizeof(float)*3*cs->NSpheroid,0,(char*)cs->SpheroidNormal); 
          PRINTFD(FB_ObjectMolecule)
            " ObjectMolPMO2CoorSet: saved spheroid size %d %d\n",cs->SpheroidSphereSize,cs->NSpheroid
            ENDFD;
         
        }
        
        /* write bonds */
        b=I->Bond;
        for(c=0;c<I->NBond;c++) {
          b1 = b->index[0];
          b2 = b->index[1];
          if(I->DiscreteFlag) {
            if((cs==I->DiscreteCSet[b1])&&(cs==I->DiscreteCSet[b2])) {
              a1=I->DiscreteAtmToIdx[b1];
              a2=I->DiscreteAtmToIdx[b2];
            } else {
              a1=-1;
              a2=-1;
            }
          } else {
            a1=cs->AtmToIdx[b1];
            a2=cs->AtmToIdx[b2];
          }
          if((a1>=0)&&(a2>=0)) { 
            nBond++;
            VLACheck(bondVLA,BondType,nBond);
            bondVLA[nBond-1]=*b;
            bondVLA[nBond-1].index[0] = a1;
            bondVLA[nBond-1].index[1] = a2;
          }
          b++;

        }
        if(ok) ok = RawWrite(raw,cRaw_Bonds1,sizeof(BondType)*nBond,0,(char*)bondVLA);
      }
    }
  }
  if(raw) RawFree(raw);
  VLAFreeP(aiVLA);
  VLAFreeP(bondVLA);
  return(ok);
}
/*========================================================================*/
int ObjectMoleculeGetPhiPsi(ObjectMolecule *I,int ca,float *phi,float *psi,int state)
{
  int np=-1;
  int cm=-1;
  int c=-1;
  int n=-1;
  int result = false;
  AtomInfoType *ai;
  int n0,at;
  float v_ca[3];
  float v_n[3];
  float v_c[3];
  float v_cm[3];
  float v_np[3];

  ai=I->AtomInfo;

  if((ai[ca].name[0]=='C')&&(ai[ca].name[1]=='A'))
    {
      ObjectMoleculeUpdateNeighbors(I);
      
      /* find C */
      n0 = I->Neighbor[ca]+1;
      while(I->Neighbor[n0]>=0) {
        at = I->Neighbor[n0];
        if((ai[at].name[0]=='C')&&(ai[at].name[1]==0)) {
          c=at;
          break;
        }
        n0+=2;
      }
      
      /* find N */
      n0 = I->Neighbor[ca]+1;
      while(I->Neighbor[n0]>=0) {
        at = I->Neighbor[n0];
        if((ai[at].name[0]=='N')&&(ai[at].name[1]==0)) {
          n=at;
          break;
        }
        n0+=2;
      }
      
      /* find NP */
      if(c>=0) {
        n0 = I->Neighbor[c]+1;
        while(I->Neighbor[n0]>=0) {
          at = I->Neighbor[n0];
          if((ai[at].name[0]=='N')&&(ai[at].name[1]==0)) {
            np=at;
            break;
          }
        n0+=2;
        }
      }
      
      /* find CM */
      if(n>=0) {
        n0 = I->Neighbor[n]+1;
        while(I->Neighbor[n0]>=0) {
          at = I->Neighbor[n0];
          if((ai[at].name[0]=='C')&&(ai[at].name[1]==0)) {
            cm=at;
            break;
          }
          n0+=2;
        }
      }
      if((ca>=0)&&(np>=0)&&(c>=0)&&(n>=0)&&(cm>=0)) {
        if(ObjectMoleculeGetAtomVertex(I,state,ca,v_ca)&&
           ObjectMoleculeGetAtomVertex(I,state,n,v_n)&&
           ObjectMoleculeGetAtomVertex(I,state,c,v_c)&&
           ObjectMoleculeGetAtomVertex(I,state,cm,v_cm)&&
           ObjectMoleculeGetAtomVertex(I,state,np,v_np)) {

          (*phi)=rad_to_deg(get_dihedral3f(v_c,v_ca,v_n,v_cm));
          (*psi)=rad_to_deg(get_dihedral3f(v_np,v_c,v_ca,v_n));
          result=true;
        }
      }
    }
  return(result);
}
/*========================================================================*/
int ObjectMoleculeCheckBondSep(ObjectMolecule *I,int a0,int a1,int dist)
{
  int result = false;
  int n0;
  int stack[MAX_BOND_DIST+1];
  int history[MAX_BOND_DIST+1];
  int depth=0;
  int distinct;
  int a;
  if(dist>MAX_BOND_DIST)
    return false;
  
  ObjectMoleculeUpdateNeighbors(I);

  PRINTFD(FB_ObjectMolecule)
    " CBS-Debug: %s %d %d %d\n",I->Obj.Name,a0,a1,dist
    ENDFD;
  depth = 1;
  history[depth]=a0;
  stack[depth] = I->Neighbor[a0]+1; /* go to first neighbor */
  while(depth) { /* keep going until we've traversed tree */
    while(I->Neighbor[stack[depth]]>=0) /* end of branches? go back up one bond */
      {
        n0 = I->Neighbor[stack[depth]]; /* get current neighbor index */
        stack[depth]+=2; /* set up next neighbor */
        distinct=true; /* check to see if current candidate is distinct from ancestors */
        for(a=1;a<depth;a++) {
          if(history[a]==n0)
            distinct=false;
        }
        if(distinct) {
          if(depth<dist) { /* are not yet at the proper distance? */
            if(distinct) {
              depth++; 
              stack[depth] = I->Neighbor[n0]+1; /* then keep moving outward */
              history[depth] = n0;
            }
          } else if(n0==a1) /* otherwise, see if we have a match */
            result = true;
        }
      }
    depth--;
  }
  PRINTFD(FB_ObjectMolecule)
    " CBS-Debug: result %d\n",result
    ENDFD;
  return result;
}
/*========================================================================*/
void ObjectGotoState(ObjectMolecule *I,int state)
{
  if((I->NCSet>1)||(!SettingGet(cSetting_static_singletons))) {
    if(state>I->NCSet)
      state = I->NCSet-1;
    if(state<0)
      state = I->NCSet-1;
    SceneSetFrame(0,state);
  }
}
/*========================================================================*/
CSetting **ObjectMoleculeGetSettingHandle(ObjectMolecule *I,int state)
{
  
  if(state<0) {
    return(&I->Obj.Setting);
  } else if(state<I->NCSet) {
    if(I->CSet[state]) {
      return(&I->CSet[state]->Setting);
    } else {
      return(NULL);
    }
  } else {
    return(NULL);
  }
}
/*========================================================================*/
int ObjectMoleculeSetStateTitle(ObjectMolecule *I,int state,char *text)
{
  int result=false;
  if(state<0) state=I->NCSet-1;
  if(state>=I->NCSet) {
    PRINTFB(FB_ObjectMolecule,FB_Errors)
      "Error: invalid state %d\n",state +1
      ENDFB;
    
  } else if(!I->CSet[state]) {
    PRINTFB(FB_ObjectMolecule,FB_Errors)
      "Error: empty state %d\n",state +1
      ENDFB;
  } else {
    UtilNCopy(I->CSet[state]->Name,text,sizeof(WordType));
    result=true;
  }
  return(result);
}

/*========================================================================*/
char *ObjectMoleculeGetStateTitle(ObjectMolecule *I,int state)
{
  char *result=NULL;
  if(state<0) state=I->NCSet-1;
  if(state>=I->NCSet) {
    PRINTFB(FB_ObjectMolecule,FB_Errors)
      "Error: invalid state %d\n",state +1
      ENDFB;
  } else if(!I->CSet[state]) {
    PRINTFB(FB_ObjectMolecule,FB_Errors)
      "Error: empty state %d\n",state +1
      ENDFB;
  } else {
    result = I->CSet[state]->Name;
  }
  return(result);
}

/*========================================================================*/
void ObjectMoleculeRenderSele(ObjectMolecule *I,int curState,int sele)
{
  CoordSet *cs;
  int a,at;

  if(PMGUI) {
    if(curState>=0) {
      if(curState<I->NCSet) {
        if(I->CSet[curState]) {
          cs=I->CSet[curState];
          for(a=0;a<cs->NIndex;a++) {
            at=cs->IdxToAtm[a]; /* should work for both discrete and non-discrete objects */
            if(SelectorIsMember(I->AtomInfo[at].selEntry,sele))
              glVertex3fv(cs->Coord+3*a);
          }
        }
      } else if(SettingGet(cSetting_static_singletons)) {
        if(I->NCSet==1) {
          cs=I->CSet[0];
          if(cs) {
            for(a=0;a<cs->NIndex;a++) {
              at=cs->IdxToAtm[a]; /* should work for both discrete and non-discrete objects */
              if(SelectorIsMember(I->AtomInfo[at].selEntry,sele))
                glVertex3fv(cs->Coord+3*a);
            }
          }
        }
      }
    } else { /* all states */
      for(curState=0;curState<I->NCSet;curState++) {
        if(I->CSet[curState]) {
          cs=I->CSet[curState];
          for(a=0;a<cs->NIndex;a++) {
            at=cs->IdxToAtm[a]; /* should work for both discrete and non-discrete objects */
            if(SelectorIsMember(I->AtomInfo[at].selEntry,sele))
              glVertex3fv(cs->Coord+3*a);
          }
        }
      }
    }
  }
}

/*========================================================================*/
CoordSet *ObjectMoleculeXYZStr2CoordSet(char *buffer,AtomInfoType **atInfoPtr)
{
  char *p;
  int nAtom;
  int a,c;
  float *coord = NULL;
  CoordSet *cset = NULL;
  AtomInfoType *atInfo = NULL,*ai;
  char cc[MAXLINELEN];
  int atomCount;
  BondType *bond=NULL;
  int nBond=0;
  int b1,b2;
  WordType tmp_name;
  int auto_show_lines = (int)SettingGet(cSetting_auto_show_lines);
  int auto_show_nonbonded = (int)SettingGet(cSetting_auto_show_nonbonded);
  BondType *ii;


  p=buffer;
  nAtom=0;
  atInfo = *atInfoPtr;
  
  p=ncopy(cc,p,6);  
  if(!sscanf(cc,"%d",&nAtom)) nAtom=0;
  p=nskip(p,2);
  p=ncopy(tmp_name,p,sizeof(WordType)-1);
  p=nextline_top(p);
      
  coord=VLAlloc(float,3*nAtom);

  if(atInfo)
	 VLACheck(atInfo,AtomInfoType,nAtom);
  
  nBond=0;
  bond=VLAlloc(BondType,6*nAtom);  
  ii=bond;

  PRINTFB(FB_ObjectMolecule,FB_Blather)
	 " ObjectMoleculeReadXYZ: Found %i atoms...\n",nAtom
    ENDFB;

  a=0;
  atomCount=0;
  
  while(*p)
	 {
      ai=atInfo+atomCount;
      
      p=ncopy(cc,p,6);
      if(!sscanf(cc,"%d",&ai->id)) break;
      
      p=nskip(p,2);/* to 12 */
      p=ncopy(cc,p,3); 
      if(!sscanf(cc,"%s",ai->name)) ai->name[0]=0;
      
      ai->alt[0]=0;
      strcpy(ai->resn,"UNK");
      ai->chain[0] = 0;
      
      ai->resv=atomCount+1;
      sprintf(ai->resi,"%d",ai->resv);
      
      p=ncopy(cc,p,12);
      sscanf(cc,"%f",coord+a);
      p=ncopy(cc,p,12);
      sscanf(cc,"%f",coord+(a+1));
      p=ncopy(cc,p,12);
      sscanf(cc,"%f",coord+(a+2));
      
      ai->q=1.0;
      ai->b=0.0;
      
      ai->segi[0]=0;
      ai->elem[0]=0; /* let atom info guess/infer atom type */
      
      for(c=0;c<cRepCnt;c++) {
        ai->visRep[c] = false;
      }
      ai->visRep[cRepLine] = auto_show_lines; /* show lines by default */
      ai->visRep[cRepNonbonded] = auto_show_nonbonded; /* show lines by default */
      
      p=ncopy(cc,p,6);
      sscanf(cc,"%d",&ai->customType);
      
      /* in the absense of external tinker information, assume hetatm */
      
      ai->hetatm=1;
      
      AtomInfoAssignParameters(ai);
      ai->color=AtomInfoGetColor(ai);
      
      b1 = atomCount;
      for(c=0;c<6;c++) {
        p=ncopy(cc,p,6);
        if (!cc[0]) 
          break;
        if(!sscanf(cc,"%d",&b2))
          break;
        if(b1<(b2-1)) {
          nBond++;
          ii->index[0] = b1;
          ii->index[1] = b2-1;
          ii->order = 1; /* missing bond order information */
          ii->stereo = 0;
          ii->id = -1; /* no serial number */
        }
      }
      
      PRINTFD(FB_ObjectMolecule) 
        " ObjectMolecule-DEBUG: %s %s %s %s %8.3f %8.3f %8.3f %6.2f %6.2f %s\n",
        ai->name,ai->resn,ai->resi,ai->chain,
        *(coord+a),*(coord+a+1),*(coord+a+2),ai->b,ai->q,
        ai->segi
        ENDFD;
      
      a+=3;
      atomCount++;
      if(atomCount>=nAtom)
        break;
      p=nextline_top(p);
    }

  PRINTFB(FB_ObjectMolecule,FB_Blather) 
   " XYZStr2CoordSet: Read %d bonds.\n",nBond
    ENDFB;

  cset = CoordSetNew();
  cset->NIndex=nAtom;
  cset->Coord=coord;
  cset->TmpBond=bond;
  cset->NTmpBond=nBond;
  strcpy(cset->Name,tmp_name);
  if(atInfoPtr)
	 *atInfoPtr = atInfo;
  return(cset);
}

/*========================================================================*/
ObjectMolecule *ObjectMoleculeReadXYZStr(ObjectMolecule *I,char *PDBStr,int frame,int discrete)
{
  CoordSet *cset = NULL;
  AtomInfoType *atInfo;
  int ok=true;
  int isNew = true;
  unsigned int nAtom = 0;

  if(!I) 
	 isNew=true;
  else 
	 isNew=false;

  if(ok) {

	 if(isNew) {
		I=(ObjectMolecule*)ObjectMoleculeNew(discrete);
		atInfo = I->AtomInfo;
		isNew = true;
	 } else {
		atInfo=VLAMalloc(10,sizeof(AtomInfoType),2,true); /* autozero here is important */
		isNew = false;
	 }
    if(isNew) {
      AtomInfoPrimeColors();
      I->Obj.Color = AtomInfoGetCarbColor();
    }
    
	 cset=ObjectMoleculeXYZStr2CoordSet(PDBStr,&atInfo);	 
	 nAtom=cset->NIndex;
  }

  /* include coordinate set */
  if(ok) {
    cset->Obj = I;
    cset->fEnumIndices(cset);
    if(cset->fInvalidateRep)
      cset->fInvalidateRep(cset,cRepAll,cRepInvRep);
    if(isNew) {		
      I->AtomInfo=atInfo; /* IMPORTANT to reassign: this VLA may have moved! */
    } else {
      ObjectMoleculeMerge(I,atInfo,cset,false,cAIC_IDMask); /* NOTE: will release atInfo */
    }

    if(isNew) I->NAtom=nAtom;
    if(frame<0) frame=I->NCSet;
    VLACheck(I->CSet,CoordSet*,frame);
    if(I->NCSet<=frame) I->NCSet=frame+1;
    if(I->CSet[frame]) I->CSet[frame]->fFree(I->CSet[frame]);
    I->CSet[frame] = cset;
    if(isNew) I->NBond = ObjectMoleculeConnect(I,&I->Bond,I->AtomInfo,cset,false);
    if(cset->Symmetry&&(!I->Symmetry)) {
      I->Symmetry=SymmetryCopy(cset->Symmetry);
      SymmetryAttemptGeneration(I->Symmetry,false,false);
    }
    SceneCountFrames();
    ObjectMoleculeExtendIndices(I);
    ObjectMoleculeSort(I);
    ObjectMoleculeUpdateIDNumbers(I);
    ObjectMoleculeUpdateNonbonded(I);
  }
  return(I);
}
/*========================================================================*/
ObjectMolecule *ObjectMoleculeLoadXYZFile(ObjectMolecule *obj,char *fname,int frame,int discrete)
{
  ObjectMolecule *I=NULL;
  int ok=true;
  FILE *f;
  long size;
  char *buffer,*p;

  f=fopen(fname,"rb");
  if(!f)
	 ok=ErrMessage("ObjectMoleculeLoadXYZFile","Unable to open file!");
  else
	 {
      PRINTFB(FB_ObjectMolecule,FB_Blather) 
        " ObjectMoleculeLoadXYZFile: Loading from %s.\n",fname
        ENDFB;
		
		fseek(f,0,SEEK_END);
      size=ftell(f);
		fseek(f,0,SEEK_SET);

		buffer=(char*)mmalloc(size+255);
		ErrChkPtr(buffer);
		p=buffer;
		fseek(f,0,SEEK_SET);
		fread(p,size,1,f);
		p[size]=0;
		fclose(f);

		I=ObjectMoleculeReadXYZStr(obj,buffer,frame,discrete);

		mfree(buffer);
	 }

  return(I);
}


/*========================================================================*/
int ObjectMoleculeAreAtomsBonded(ObjectMolecule *I,int i0,int i1)
{
  int result=false;
  int a;
  BondType *b;
  b=I->Bond;
  for (a=0;a<I->NBond;a++) {
    if(i0==b->index[0]) {
      if(i1==b->index[1]) {
        result=true;
        break;
      }
    }
    if(i1==b->index[0]) {
      if(i0==b->index[1]) {
        result=true;
        break;
      }
    }
    b++;
  }
  return(result);
}
/*========================================================================*/
void ObjectMoleculeRenameAtoms(ObjectMolecule *I,int force)
{
  AtomInfoType *ai;
  int a;
  if(force) {
    ai=I->AtomInfo;
    for(a=0;a<I->NAtom;a++)
      (ai++)->name[0]=0;
  }
  AtomInfoUniquefyNames(NULL,0,I->AtomInfo,I->NAtom);  
}
/*========================================================================*/
void ObjectMoleculeAddSeleHydrogens(ObjectMolecule *I,int sele)
{
  int a,b;
  int n,nn;
  CoordSet *cs;
  CoordSet *tcs;
  int seleFlag=false;
  AtomInfoType *ai,*nai,fakeH;
  int repeatFlag=false;
  int nH;
  int *index;
  float v[3],v0[3];
  float d;

  UtilZeroMem(&fakeH,sizeof(AtomInfoType));
  fakeH.protons=1;
  ai=I->AtomInfo;
  for(a=0;a<I->NAtom;a++) {
    if(SelectorIsMember(ai->selEntry,sele)) {
      seleFlag=true;
      break;
    }
    ai++;
  }
  if(seleFlag) {
    if(!ObjectMoleculeVerifyChemistry(I)) {
      ErrMessage(" AddHydrogens","missing chemical geometry information.");
    } else if(I->DiscreteFlag) {
      ErrMessage(" AddHydrogens","can't modify a discrete object.");
    } else {

      repeatFlag=true;
      while(repeatFlag) {
        repeatFlag=false;
        nH = 0;
        ObjectMoleculeUpdateNeighbors(I);
        nai = (AtomInfoType*)VLAMalloc(1000,sizeof(AtomInfoType),1,true);        
        ai=I->AtomInfo;
        for(a=0;a<I->NAtom;a++) {
          if(SelectorIsMember(ai->selEntry,sele)) {
            n = I->Neighbor[a];
            nn = I->Neighbor[n++];
            if(nn<ai->valence) {
              VLACheck(nai,AtomInfoType,nH);
              UtilNCopy((nai+nH)->elem,"H",2);
              (nai+nH)->geom=cAtomInfoSingle;
              (nai+nH)->valence=1;
              (nai+nH)->temp1 = a; /* borrowing this field temporarily */
              ObjectMoleculePrepareAtom(I,a,nai+nH);
              nH++;
            }
          }
          ai++;
        }

        if(nH) {

          repeatFlag=true;
          cs = CoordSetNew();
          cs->Coord = VLAlloc(float,nH*3);
          cs->NIndex=nH;

          index = Alloc(int,nH);
          for(a=0;a<nH;a++) {
            index[a] = (nai+a)->temp1;
          }
          
          if(cs->fEnumIndices) cs->fEnumIndices(cs);

          cs->TmpLinkBond = VLAlloc(BondType,nH);
          for(a=0;a<nH;a++) {
            cs->TmpLinkBond[a].index[0] = (nai+a)->temp1;
            cs->TmpLinkBond[a].index[1] = a;
            cs->TmpLinkBond[a].order = 1;
            cs->TmpLinkBond[a].stereo = 0;
            cs->TmpLinkBond[a].id = -1;
          }
          cs->NTmpLinkBond = nH;

          AtomInfoUniquefyNames(I->AtomInfo,I->NAtom,nai,nH);

          ObjectMoleculeMerge(I,nai,cs,false,cAIC_AllMask); /* will free nai and cs->TmpLinkBond  */
          ObjectMoleculeExtendIndices(I);
          ObjectMoleculeUpdateNeighbors(I);

          for(b=0;b<I->NCSet;b++) { /* add coordinate into the coordinate set */
            tcs = I->CSet[b];
            if(tcs) {
              for(a=0;a<nH;a++) {
                ObjectMoleculeGetAtomVertex(I,b,index[a],v0);
                ObjectMoleculeFindOpenValenceVector(I,b,index[a],v,NULL);
                d = AtomInfoGetBondLength(I->AtomInfo+index[a],&fakeH);
                scale3f(v,d,v);
                add3f(v0,v,cs->Coord+3*a);
              }
              CoordSetMerge(tcs,cs);  
            }
          }
          FreeP(index);
          if(cs->fFree)
            cs->fFree(cs);
          ObjectMoleculeSort(I);
          ObjectMoleculeUpdateIDNumbers(I);
        } else
          VLAFreeP(nai);
      }
    }
  }
}


/*========================================================================*/
void ObjectMoleculeFuse(ObjectMolecule *I,int index0,ObjectMolecule *src,int index1,int mode)
{
  int a,b;
  AtomInfoType *ai0,*ai1,*nai;
  int n,nn;
  int at0=-1;
  int at1=-1;
  int a0,a1;
  int hydr1=-1;
  int anch1=-1;
  int ca0,ch0;
  BondType *b0,*b1;
  float *backup = NULL;
  float d,*f0,*f1;
  float va0[3],vh0[3],va1[3],vh1[3];
  float x0[3],y0[3],z0[3];
  float x1[3],y1[3],z1[3];
  float x[3],y[3],z[3];
  float t[3],t2[3];
  CoordSet *cs=NULL,*scs=NULL;
  int state1 = 0;
  CoordSet *tcs;
  int edit=1;
  OrthoLineType sele1,sele2,s1,s2;

  ObjectMoleculeUpdateNeighbors(I);
  ObjectMoleculeUpdateNeighbors(src);

  /* make sure each link point has only one neighbor */

  ai0=I->AtomInfo;
  ai1=src->AtomInfo;
  switch(mode) {
  case 0: /* fusing by replacing hydrogens */
    
    n = I->Neighbor[index0];
    nn = I->Neighbor[n++];
    if(nn==1)
      at0 = I->Neighbor[n];
    
    n = src->Neighbor[index1];
    nn = src->Neighbor[n++];
    if(nn==1)
      at1 = src->Neighbor[n];
    
    if(src->NCSet) {
      scs = src->CSet[state1];
      anch1 = scs->AtmToIdx[at1];
      hydr1 = scs->AtmToIdx[index1];
    }
    break;
  case 1: /* fuse merely by drawing a bond */
    at0 = index0;
    at1 = index1;
    
    if(src->NCSet) {
      scs = src->CSet[state1];
      anch1 = scs->AtmToIdx[at1];
    }

    break;
  }
  
  if((at0>=0)&&(at1>=0)&&scs&&(anch1>=0)) { /* have anchors and source coordinate set */

    nai=(AtomInfoType*)VLAMalloc(src->NAtom,sizeof(AtomInfoType),1,true);
    
    /* copy atoms and atom info into a 1:1 direct mapping */

    cs = CoordSetNew();
    cs->Coord = VLAlloc(float,scs->NIndex*3);
    cs->NIndex = scs->NIndex;
    for(a=0;a<scs->NIndex;a++) {
      copy3f(scs->Coord+a*3,cs->Coord+a*3);
      a1 = scs->IdxToAtm[a];
      *(nai+a) = *(ai1+a1);
      (nai+a)->selEntry=0; /* avoid duplicating selection references -> leads to hangs ! */
      (nai+a)->temp1=0; /* clear marks */
    }

    nai[at1].temp1=2; /* mark the connection point */
    
    /* copy internal bond information*/

    cs->TmpBond = VLAlloc(BondType,src->NBond);
    b1 = src->Bond;
    b0 = cs->TmpBond;
    cs->NTmpBond=0;
    for(a=0;a<src->NBond;a++) {
      a0 = scs->AtmToIdx[b1->index[0]];
      a1 = scs->AtmToIdx[b1->index[1]];
      if((a0>=0)&&(a1>=0)) {
        *b0=*b1;
        b0->index[0] = a0;
        b0->index[1] = a1;
        b0++;
        cs->NTmpBond++;
      }
      b1++;
    }

    backup = Alloc(float,cs->NIndex*3); /* make untransformed copy of coordinate set */
    for(a=0;a<cs->NIndex;a++) {
      copy3f(cs->Coord+a*3,backup+a*3);
    }
    
    switch(mode) {
    case 0:
      nai[hydr1].deleteFlag=true;
      I->AtomInfo[index0].deleteFlag=true;
      copy3f(backup+3*anch1,va1);
      copy3f(backup+3*hydr1,vh1);
      subtract3f(va1,vh1,x1); /* note reverse dir from above */
      get_system1f3f(x1,y1,z1);
      break;
    case 1:
      copy3f(backup+3*anch1,va1);
      ObjectMoleculeFindOpenValenceVector(src,state1,at1,x1,NULL);
      scale3f(x1,-1.0F,x1);
      get_system1f3f(x1,y1,z1);      
      break;
    }

    /* set up the linking bond */

    cs->TmpLinkBond = VLAlloc(BondType,1);
    cs->NTmpLinkBond = 1;
    cs->TmpLinkBond->index[0] = at0;
    cs->TmpLinkBond->index[1] = anch1;
    cs->TmpLinkBond->order = 1;
    cs->TmpLinkBond->stereo = 0;
    cs->TmpLinkBond->id = -1;
    
    if(cs->fEnumIndices) cs->fEnumIndices(cs);

    d = AtomInfoGetBondLength(ai0+at0,ai1+at1);

    AtomInfoUniquefyNames(I->AtomInfo,I->NAtom,nai,cs->NIndex);

    /* set up tags which will enable use to continue editing bond */

    if(edit) {
      for(a=0;a<I->NAtom;a++) {
        I->AtomInfo[a].temp1=0;
      }
      I->AtomInfo[at0].temp1=1;
    }

    ObjectMoleculeMerge(I,nai,cs,false,cAIC_AllMask); /* will free nai, cs->TmpBond and cs->TmpLinkBond  */

    ObjectMoleculeExtendIndices(I);
    ObjectMoleculeUpdateNeighbors(I);
    for(a=0;a<I->NCSet;a++) { /* add coordinate into the coordinate set */
      tcs = I->CSet[a];
      if(tcs) {
        switch(mode) {
        case 0:
          ca0 = tcs->AtmToIdx[at0]; /* anchor */
          ch0 = tcs->AtmToIdx[index0]; /* hydrogen */

          if((ca0>=0)&&(ch0>=0)) {
            copy3f(tcs->Coord+3*ca0,va0);
            copy3f(tcs->Coord+3*ch0,vh0);
            subtract3f(vh0,va0,x0);
            get_system1f3f(x0,y0,z0);

          }
          break;
        case 1:
          ca0 = tcs->AtmToIdx[at0]; /* anchor */

          if(ca0>=0) {
            ObjectMoleculeFindOpenValenceVector(I,a,at0,x0,NULL);
            copy3f(tcs->Coord+3*ca0,va0);
            get_system1f3f(x0,y0,z0);
            
          }
          break;
        }
        scale3f(x0,d,t2);
        add3f(va0,t2,t2);
        
        f0=backup;
        f1=cs->Coord;
        for(b=0;b<cs->NIndex;b++) { /* brute force transformation */
          subtract3f(f0,va1,t);
          scale3f(x0,dot_product3f(t,x1),x);
          scale3f(y0,dot_product3f(t,y1),y);
          scale3f(z0,dot_product3f(t,z1),z);
          add3f(x,y,y);
          add3f(y,z,f1);
          add3f(t2,f1,f1);
          f0+=3;
          f1+=3;
        }
        CoordSetMerge(tcs,cs); 
      }
    }
    switch(mode) {
    case 0:
      ObjectMoleculePurge(I);
      break;
    }
    ObjectMoleculeSort(I);
    ObjectMoleculeUpdateIDNumbers(I);
    if(edit) { /* edit the resulting bond */
      at0=-1;
      at1=-1;
      for(a=0;a<I->NAtom;a++) {
        if(I->AtomInfo[a].temp1==1)
          at0=a;
        if(I->AtomInfo[a].temp1==2)
          at1=a;
      }
      if((at0>=0)&&(at1>=0)) {
        sprintf(sele1,"%s`%d",I->Obj.Name,at1+1); /* points outward... */
        sprintf(sele2,"%s`%d",I->Obj.Name,at0+1);
        SelectorGetTmp(sele1,s1);
        SelectorGetTmp(sele2,s2);
        EditorSelect(s1,s2,NULL,NULL,false,true,true);
        SelectorFreeTmp(s1);
        SelectorFreeTmp(s2);
      }
    }
  }
  if(cs)
    if(cs->fFree)
      cs->fFree(cs);
  FreeP(backup);
}
/*========================================================================*/
int ObjectMoleculeVerifyChemistry(ObjectMolecule *I)
{
  int result=false;
  AtomInfoType *ai;
  int a;
  int flag;
  ai=I->AtomInfo;
  flag=true;
  for(a=0;a<I->NAtom;a++) {
    if(!ai->chemFlag) {
      flag=false;
    }
    ai++;
  }
  if(!flag) {
    if(I->CSet[0]) { /* right now this stuff is locked to state 0 */
      ObjectMoleculeInferChemFromBonds(I,0);
      ObjectMoleculeInferChemFromNeighGeom(I,0);
      ObjectMoleculeInferHBondFromChem(I);
      /*      ObjectMoleculeInferChemForProtein(I,0);*/
    }
    flag=true;
    ai=I->AtomInfo;
    for(a=0;a<I->NAtom;a++) {
      if(!ai->chemFlag) {
        flag=false;
        break;
      }
      ai++;
    }
  }
  if(flag)
    result=true;
  return(result);
}
/*========================================================================*/
void ObjectMoleculeAttach(ObjectMolecule *I,int index,AtomInfoType *nai)
{
  int a;
  AtomInfoType *ai;
  int n,nn;
  float v[3],v0[3],d;
  CoordSet *cs;

  ObjectMoleculeUpdateNeighbors(I);
  ai=I->AtomInfo+index;
  n = I->Neighbor[index];
  nn = I->Neighbor[n++];
  
  cs = CoordSetNew();
  cs->Coord = VLAlloc(float,3);
  cs->NIndex=1;
  cs->TmpLinkBond = VLAlloc(BondType,1);
  cs->NTmpLinkBond = 1;
  cs->TmpLinkBond->index[0]=index;
  cs->TmpLinkBond->index[1]=0;
  cs->TmpLinkBond->order=1;
  cs->TmpLinkBond->stereo=0;
  
  cs->TmpLinkBond->id = -1;
  if(cs->fEnumIndices) cs->fEnumIndices(cs);
  ObjectMoleculePrepareAtom(I,index,nai);
  d = AtomInfoGetBondLength(ai,nai);
  ObjectMoleculeMerge(I,nai,cs,false,cAIC_AllMask); /* will free nai and cs->TmpLinkBond  */
  ObjectMoleculeExtendIndices(I);
  ObjectMoleculeUpdateNeighbors(I);
  for(a=0;a<I->NCSet;a++) { /* add atom to each coordinate set */
    if(I->CSet[a]) {
      ObjectMoleculeGetAtomVertex(I,a,index,v0);
      ObjectMoleculeFindOpenValenceVector(I,a,index,v,NULL);
      scale3f(v,d,v);
      add3f(v0,v,cs->Coord);
      CoordSetMerge(I->CSet[a],cs); 
    }
  }
  ObjectMoleculeSort(I);
  ObjectMoleculeUpdateIDNumbers(I);
  if(cs->fFree)
    cs->fFree(cs);
  
}
/*========================================================================*/
int ObjectMoleculeFillOpenValences(ObjectMolecule *I,int index)
{
  int a;
  AtomInfoType *ai,*nai;
  int n,nn;
  int result=0;
  int flag = true;
  float v[3],v0[3],d;
  CoordSet *cs;

  if((index>=0)&&(index<=I->NAtom)) {  
    while(1) {
      ObjectMoleculeUpdateNeighbors(I);
      ai=I->AtomInfo+index;
      n = I->Neighbor[index];
      nn = I->Neighbor[n++];
      
      if((nn>=ai->valence)||(!flag))
        break;
      flag=false;

      cs = CoordSetNew();
      cs->Coord = VLAlloc(float,3);
      cs->NIndex=1;
      cs->TmpLinkBond = VLAlloc(BondType,1);
      cs->NTmpLinkBond = 1;
      cs->TmpLinkBond->index[0]=index;
      cs->TmpLinkBond->index[1]=0;
      cs->TmpLinkBond->order=1;
      cs->TmpLinkBond->stereo=0;
      
      cs->TmpLinkBond->id = -1;
      if(cs->fEnumIndices) cs->fEnumIndices(cs);
      nai = (AtomInfoType*)VLAMalloc(1,sizeof(AtomInfoType),1,true);
      UtilNCopy(nai->elem,"H",2);
      nai->geom=cAtomInfoSingle;
      nai->valence=1;
      ObjectMoleculePrepareAtom(I,index,nai);
      d = AtomInfoGetBondLength(ai,nai);
      ObjectMoleculeMerge(I,nai,cs,false,cAIC_AllMask); /* will free nai and cs->TmpLinkBond  */
      ObjectMoleculeExtendIndices(I);
      ObjectMoleculeUpdateNeighbors(I);
      for(a=0;a<I->NCSet;a++) { /* add atom to each coordinate set */
        if(I->CSet[a]) {
          ObjectMoleculeGetAtomVertex(I,a,index,v0);
          ObjectMoleculeFindOpenValenceVector(I,a,index,v,NULL);
          scale3f(v,d,v);
          add3f(v0,v,cs->Coord);
          CoordSetMerge(I->CSet[a],cs); 
        }
      }
      if(cs->fFree)
        cs->fFree(cs);
      result++;
      flag=true;
    }
  }
  ObjectMoleculeUpdateIDNumbers(I);
  return(result);
}

  #define MaxOcc 100

/*========================================================================*/
static int get_planer_normal(ObjectMolecule *I,int state,int index,float *normal)
{   /* NOTE assumes neighbors are defined */
  int found = false;
  int nOcc = 0;
  float occ[MaxOcc*3];
  AtomInfoType *ai = I->AtomInfo + index;
  int n,a1;
  float v0[3],v1[3],v2[3],n0[3];

  if(ObjectMoleculeGetAtomVertex(I,state,index,v0)) {        
    n = I->Neighbor[index];
    n++; /* skip count */
    while(1) { /* look for an attached non-hydrogen as a base */
      a1 = I->Neighbor[n];
      n+=2; 
      if(a1<0) break;
      if(ObjectMoleculeGetAtomVertex(I,state,a1,v1)) {        
        subtract3f(v1,v0,n0);
        normalize3f(n0); /* n0's point away from center atom */
        copy3f(n0,occ+3*nOcc);
        nOcc++; 
        if(nOcc==MaxOcc) /* safety valve */
          break;
      }
    }
    switch(ai->geom) {
    case cAtomInfoPlaner:
      if(nOcc>1) {
        cross_product3f(occ,occ+3,normal);
        if(nOcc>2) {
          cross_product3f(occ,occ+6,v2);
          if(dot_product3f(normal,v2)<0) {
            subtract3f(normal,v2,normal);
          } else {
            add3f(normal,v2,normal);
          }
          cross_product3f(occ+3,occ+6,v2);
          if(dot_product3f(normal,v2)<0) {
            subtract3f(normal,v2,normal);
          } else {
            add3f(normal,v2,normal);
          }
        }
        normalize3f(normal);
        found=true;
      }
      break;
    }
  }
  return found;
}

/*========================================================================*/
int ObjectMoleculeFindOpenValenceVector(ObjectMolecule *I,int state,
                                        int index,float *v,float *seek)
{
  CoordSet *cs;
  int nOcc = 0;
  float occ[MaxOcc*3];
  int last_occ = -1;
  int n;
  int a1;
  float v0[3],v1[3],n0[3],t[3];
  int result = false;
  AtomInfoType *ai,*ai1;
  float y[3],z[3];

  /* default is +X */
  v[0]=1.0;
  v[1]=0.0;
  v[2]=0.0;
  
  if(state<0) state=0;
  if(I->NCSet==1) state=0;
  state = state % I->NCSet;
  cs = I->CSet[state];
  if(cs) {
    if((index>=0)&&(index<=I->NAtom)) {
      ai=I->AtomInfo+index;
      if(ObjectMoleculeGetAtomVertex(I,state,index,v0)) {              
        ObjectMoleculeUpdateNeighbors(I);
        n = I->Neighbor[index];
        n++; /* skip count */
        while(1) { /* look for an attached non-hydrogen as a base */
          a1 = I->Neighbor[n];
          n+=2; 
          if(a1<0) break;
          ai1=I->AtomInfo+a1;
          if(ObjectMoleculeGetAtomVertex(I,state,a1,v1)) {        
            last_occ = a1;
            subtract3f(v1,v0,n0);
            normalize3f(n0); /* n0's point away from center atom */
            copy3f(n0,occ+3*nOcc);
            nOcc++; 
            if(nOcc==MaxOcc) /* safety valve */
              break;
          }
        }
        if((!nOcc)||(nOcc>4)||(ai->geom==cAtomInfoNone)) {
          if(!seek) 
            get_random3f(v);
          else
            copy3f(seek,v);
          result = true;
        } else {
          switch(nOcc) {
          case 1:  /* only one current occupied position */
            switch(ai->geom) {
            case cAtomInfoTetrahedral: 
              if(!seek) {
                get_system1f3f(occ,y,z);
                scale3f(occ,-0.334F,v);
                scale3f(z,  0.943F,t);
                add3f(t,v,v);              
              } else { /* point hydrogen towards sought vector */
                copy3f(seek,z);
                get_system2f3f(occ,z,y);
                scale3f(occ,-0.334F,v);
                scale3f(z,  0.943F,t);
                add3f(t,v,v);              
              }
              result = true;
              break;
            case cAtomInfoPlaner:
              {
                if(!seek) {
                  if((last_occ>=0)&&get_planer_normal(I,state,last_occ,n0)) {
                    copy3f(n0,y);
                    get_system2f3f(occ,y,z);
                  } else {
                    get_system1f3f(occ,y,z);
                  }
                  scale3f(occ,-0.500F,v);
                  scale3f(z,      0.866F,t);
                  add3f(t,v,v);
                } else {
                  copy3f(seek,z);
                  get_system2f3f(occ,z,y);
                  scale3f(occ,-0.500F,v);
                  scale3f(z,      0.866F,t);
                  add3f(t,v,v);
                }
                result = true;
              }
              break;
            case cAtomInfoLinear:
              scale3f(occ,-1.0F,v);
              result = true;
              break;
            default:
              if(!seek) 
                get_random3f(v);
              else
                copy3f(seek,v);
              result = true;
              break;
            }
            break;
          case 2:  /* only two current occupied positions */
            switch(ai->geom) {
            case cAtomInfoTetrahedral:
              add3f(occ,occ+3,t);
              get_system2f3f(t,occ,z);
              scale3f(t,-1.0F,v);
              if(seek) {
                if(dot_product3f(z,seek)<0.0F) {
                  invert3f(z);
                }
              } 
              scale3f(z,1.41F,t);
              add3f(t,v,v);              
              result = true;
              break;
            case cAtomInfoPlaner:
              add3f(occ,occ+3,t);
              scale3f(t,-1.0F,v);
              result = true;
              break;
            default:
              if(!seek) 
                get_random3f(v);
              else
                copy3f(seek,v);
              /* hypervalent */
              result = true;
              break;
            }
            break;
          case 3:  /* only three current occupied positions */
            switch(ai->geom) {
            case cAtomInfoTetrahedral:
              add3f(occ,occ+3,t);
              add3f(occ+6,t,t);
              scale3f(t,-1.0F,v);
              result = true;
              break;
            default:
              if(!seek) 
                get_random3f(v);
              else
                copy3f(seek,v);
              /* hypervalent */
              result = true;
              break;
            }
            break;
          case 4:
            if(!seek) 
              get_random3f(v);
            else
              copy3f(seek,v);
            /* hypervalent */
            result = true;
            break;
          }
        }
      }
    }
  }
  normalize3f(v);
  return(result); 
#undef MaxOcc
  
}
/*========================================================================*/
void ObjectMoleculeCreateSpheroid(ObjectMolecule *I,int average)
{
  CoordSet *cs;
  float *spheroid = NULL;
  int a,b,c,a0;
  SphereRec *sp;
  float *spl;
  float *v,*v0,*s,*f,ang,min_dist,*max_sq;
  int *i;
  float *center = NULL;
  float d0[3],n0[3],d1[3],d2[3];
  float p0[3],p1[3],p2[3];
  int t0,t1,t2,bt0,bt1,bt2;
  float dp,l,*fsum = NULL;
  float *norm = NULL;
  float spheroid_smooth;
  float spheroid_fill;
  float spheroid_ratio=0.1F; /* minimum ratio of width over length */
  float spheroid_minimum = 0.02F; /* minimum size - to insure valid normals */
  int row,*count=NULL,base;
  int nRow;
  int first=0;
  int last=0;
  int current;
  int cscount;
  int n_state=0;
  sp=Sphere1;
  
  nRow = I->NAtom*sp->nDot;


  center=Alloc(float,I->NAtom*3);
  count=Alloc(int,I->NAtom);
  fsum=Alloc(float,nRow);
  max_sq = Alloc(float,I->NAtom);

  spl=spheroid;

  spheroid_smooth=SettingGet(cSetting_spheroid_smooth);
  spheroid_fill=SettingGet(cSetting_spheroid_fill);
  /* first compute average coordinate */

  if(average<1)
    average=I->NCSet;
  current=0;
  cscount=0;
  while(current<I->NCSet) {
    if(I->CSet[current]) {
      if(!cscount)
        first=current;
      cscount++;
      last=current+1;
    }
    
    if(cscount==average)
      {
        PRINTFB(FB_ObjectMolecule,FB_Details)
          " ObjectMolecule: computing spheroid from states %d to %d.\n",
                 first+1,last
          ENDFB;

        spheroid=Alloc(float,nRow);
        
        v=center;
        i = count;
        for(a=0;a<I->NAtom;a++) {
          *(v++)=0.0;
          *(v++)=0.0;
          *(v++)=0.0;
          *(i++)=0;
        }

        for(b=first;b<last;b++) {
          cs=I->CSet[b];
          if(cs) {
            v = cs->Coord;
            for(a=0;a<cs->NIndex;a++) {
              a0=cs->IdxToAtm[a];
              v0 = center+3*a0;
              add3f(v,v0,v0);
              (*(count+a0))++;
              v+=3;
            }
          }
        }

        i=count;
        v=center;
        for(a=0;a<I->NAtom;a++) 
          if(*i) {
            (*(v++))/=(*i);
            (*(v++))/=(*i);
            (*(v++))/=(*i++);
          } else {
            v+=3;
            i++;
          }

        /* now go through and compute radial distances */

        f = fsum;
        s = spheroid;
        for(a=0;a<nRow;a++) {
          *(f++)=0.0;
          *(s++)=0.0; 
        }

        v = max_sq;
        for(a=0;a<I->NAtom;a++)
          *(v++)=0.0;

        for(b=first;b<last;b++) {
          cs=I->CSet[b];
          if(cs) {
            v = cs->Coord;
            for(a=0;a<cs->NIndex;a++) {
              a0=cs->IdxToAtm[a];
              base = (a0*sp->nDot);
              v0 = center+(3*a0);
              subtract3f(v,v0,d0); /* subtract from average */
              l = lengthsq3f(d0);
              if(l>max_sq[a0])
                max_sq[a0]=l;
              if(l>0.0) {
				float isq = (float)(1.0/sqrt1d(l));
                scale3f(d0,isq,n0);
                for(c=0;c<sp->nDot;c++) { /* average over spokes */
                  dp=dot_product3f(sp->dot[c],n0);
                  row = base + c;
                  if(dp>=0.0) {
					ang = (float)((acos(dp)/spheroid_smooth)*(cPI/2.0)); 
                    if(ang>spheroid_fill)
                      ang=spheroid_fill;
                    /* take envelop to zero over that angle */
                    if(ang<=(cPI/2.0)) {
                      dp = (float)cos(ang);
                      fsum[row] += dp*dp;
                      spheroid[row] += l*dp*dp*dp;
                    }
                  }
                }
              }
              v+=3;
            }
          }
        }

        f=fsum;
        s=spheroid;
        for(a=0;a<I->NAtom;a++) {
          min_dist = (float)(spheroid_ratio*sqrt(max_sq[a]));
          if(min_dist<spheroid_minimum)
            min_dist=spheroid_minimum;
          for(b=0;b<sp->nDot;b++) {
            if(*f>R_SMALL4) {
              (*s)=(float)(sqrt1d((*s)/(*(f++)))); /* we put the "rm" in "rms" */
            } else {
              f++;
            }
            if(*s<min_dist)
              *s=min_dist;
            s++;
          }
        }

        /* set frame 0 coordinates to the average */

         cs=I->CSet[first];
         if(cs) {
           v = cs->Coord;
           for(a=0;a<cs->NIndex;a++) {
             a0=cs->IdxToAtm[a];
             v0 = center+3*a0;
             copy3f(v0,v);
             v+=3;
           }
         }

        /* now compute surface normals */

        norm = Alloc(float,nRow*3);
        for(a=0;a<nRow;a++) {
          zero3f(norm+a*3);
        }
        for(a=0;a<I->NAtom;a++) {
          base = a*sp->nDot;
          for(b=0;b<sp->NTri;b++) {
            t0 = sp->Tri[b*3  ];
            t1 = sp->Tri[b*3+1];
            t2 = sp->Tri[b*3+2];
            bt0 = base + t0;
            bt1 = base + t1;
            bt2 = base + t2;
            copy3f(sp->dot[t0],p0);
            copy3f(sp->dot[t1],p1);
            copy3f(sp->dot[t2],p2);
            /*      scale3f(sp->dot[t0].v,spheroid[bt0],p0);
                    scale3f(sp->dot[t1].v,spheroid[bt1],p1);
                    scale3f(sp->dot[t2].v,spheroid[bt2],p2);*/
            subtract3f(p1,p0,d1);
            subtract3f(p2,p0,d2);
            cross_product3f(d1,d2,n0);
            normalize3f(n0);
            v = norm+bt0*3;
            add3f(n0,v,v);
            v = norm+bt1*3;
            add3f(n0,v,v);
            v = norm+bt2*3;
            add3f(n0,v,v);
          }
        }

        f=norm;
        for(a=0;a<I->NAtom;a++) {
          base = a*sp->nDot;
          for(b=0;b<sp->nDot;b++) {
            normalize3f(f);
            f+=3;
          }
        }
  
        if(I->CSet[first]) {
          I->CSet[first]->Spheroid=spheroid;
          I->CSet[first]->SpheroidNormal=norm;
          I->CSet[first]->NSpheroid=nRow;
        } else {
          FreeP(spheroid);
          FreeP(norm);
        }

        for(b=first+1;b<last;b++) { 
          cs=I->CSet[b];
          if(cs) {
            if(cs->fFree)
              cs->fFree(cs);
          }
          I->CSet[b]=NULL;
        }
        
        if(n_state!=first) {
          I->CSet[n_state]=I->CSet[first];
          I->CSet[first]=NULL;
        }
        n_state++;

        cscount=0;
      }
    current++;
  }
  I->NCSet=n_state;
  FreeP(center);
  FreeP(count);
  FreeP(fsum);
  FreeP(max_sq);

  ObjectMoleculeInvalidate(I,cRepSphere,cRepInvProp);
}
/*========================================================================*/
void ObjectMoleculeReplaceAtom(ObjectMolecule *I,int index,AtomInfoType *ai)
{
  if((index>=0)&&(index<=I->NAtom)) {
    memcpy(I->AtomInfo+index,ai,sizeof(AtomInfoType));
    ObjectMoleculeInvalidate(I,cRepAll,cRepInvAtoms);
    /* could we put in a refinement step here? */
  }
}
/*========================================================================*/
void ObjectMoleculePrepareAtom(ObjectMolecule *I,int index,AtomInfoType *ai)
{
  /* match existing properties of the old atom */
  int a;
  AtomInfoType *ai0;

  if((index>=0)&&(index<=I->NAtom)) {
    ai0=I->AtomInfo + index;
    ai->resv=ai0->resv;
    ai->hetatm=ai0->hetatm;
    ai->flags=ai0->flags;
    ai->geom=ai0->geom; /* ?*/
    strcpy(ai->chain,ai0->chain);
    strcpy(ai->alt,ai0->alt);
    strcpy(ai->resi,ai0->resi);
    strcpy(ai->segi,ai0->segi);
    strcpy(ai->resn,ai0->resn);    
    if((ai->elem[0]==ai0->elem[0])&&(ai->elem[1]==ai0->elem[1]))
      ai->color=ai0->color;
    else if((ai->elem[0]=='C')&&(ai->elem[1]==0)) 
      /* carbons are always colored according to the object color */
      ai->color=I->Obj.Color;
    else
      ai->color=AtomInfoGetColor(ai);
    for(a=0;a<cRepCnt;a++)
      ai->visRep[a]=ai0->visRep[a];
    ai->id=-1;
    AtomInfoUniquefyNames(I->AtomInfo,I->NAtom,ai,1);
    AtomInfoAssignParameters(ai);
  }
}
/*========================================================================*/
void ObjectMoleculePreposReplAtom(ObjectMolecule *I,int index,
                                   AtomInfoType *ai)
{
  int n;
  int a1;
  AtomInfoType *ai1;
  float v0[3],v1[3],v[3];
  float d0[3],d,n0[3];
  int cnt;
  float t[3],sum[3];
  int a;
  int ncycle;
  ObjectMoleculeUpdateNeighbors(I);
  for(a=0;a<I->NCSet;a++) {
    if(I->CSet[a]) {
      if(ObjectMoleculeGetAtomVertex(I,a,index,v0)) {
        copy3f(v0,v); /* default is direct superposition */
        ncycle=-1;
        while(ncycle) {
          cnt = 0;
          n = I->Neighbor[index];
          n++; /* skip count */
          zero3f(sum);
          while(1) { /* look for an attached non-hydrogen as a base */
            a1 = I->Neighbor[n];
            n+=2;
            if(a1<0) break;
            ai1=I->AtomInfo+a1;
            if(ai1->protons!=1) 
              if(ObjectMoleculeGetAtomVertex(I,a,a1,v1)) {        
                d = AtomInfoGetBondLength(ai,ai1);
                subtract3f(v0,v1,n0);
                normalize3f(n0);
                scale3f(n0,d,d0);
                add3f(d0,v1,t);
                add3f(t,sum,sum);
                cnt++;
              }
          }
          if(cnt) {
            scale3f(sum,1.0F/cnt,sum);
            copy3f(sum,v0);
            if((cnt>1)&&(ncycle<0))
              ncycle=5;
          }
          ncycle=abs(ncycle)-1;
        }
        if(cnt) copy3f(sum,v);
        ObjectMoleculeSetAtomVertex(I,a,index,v);            
      }
    }
  }
}
/*========================================================================*/
void ObjectMoleculeSaveUndo(ObjectMolecule *I,int state,int log)
{
  CoordSet *cs;

  FreeP(I->UndoCoord[I->UndoIter]);
  I->UndoState[I->UndoIter]=-1;
  if(state<0) state=0;
  if(I->NCSet==1) state=0;
  state = state % I->NCSet;
  cs = I->CSet[state];
  if(cs) {
    I->UndoCoord[I->UndoIter] = Alloc(float,cs->NIndex*3);
    memcpy(I->UndoCoord[I->UndoIter],cs->Coord,sizeof(float)*cs->NIndex*3);
    I->UndoState[I->UndoIter]=state;
    I->UndoNIndex[I->UndoIter] = cs->NIndex;
  }
  I->UndoIter=cUndoMask&(I->UndoIter+1);
  ExecutiveSetLastObjectEdited((CObject*)I);
  if(log) {
    OrthoLineType line;
    if(SettingGet(cSetting_logging)) {
      sprintf(line,"cmd.push_undo(\"%s\",%d)\n",I->Obj.Name,state+1);
      PLog(line,cPLog_no_flush);
    }
  }

}
/*========================================================================*/
void ObjectMoleculeUndo(ObjectMolecule *I,int dir)
{
  CoordSet *cs;
  int state;

  FreeP(I->UndoCoord[I->UndoIter]);
  I->UndoState[I->UndoIter]=-1;
  state=SceneGetState();
  if(state<0) state=0;
  if(I->NCSet==1) state=0;
  state = state % I->NCSet;
  cs = I->CSet[state];
  if(cs) {
    I->UndoCoord[I->UndoIter] = Alloc(float,cs->NIndex*3);
    memcpy(I->UndoCoord[I->UndoIter],cs->Coord,sizeof(float)*cs->NIndex*3);
    I->UndoState[I->UndoIter]=state;
    I->UndoNIndex[I->UndoIter] = cs->NIndex;
  }

  I->UndoIter=cUndoMask&(I->UndoIter+dir);
  if(!I->UndoCoord[I->UndoIter])
    I->UndoIter=cUndoMask&(I->UndoIter-dir);

  if(I->UndoState[I->UndoIter]>=0) {
    state=I->UndoState[I->UndoIter];
    if(state<0) state=0;
    
    if(I->NCSet==1) state=0;
    state = state % I->NCSet;
    cs = I->CSet[state];
    if(cs) {
      if(cs->NIndex==I->UndoNIndex[I->UndoIter]) {
        memcpy(cs->Coord,I->UndoCoord[I->UndoIter],sizeof(float)*cs->NIndex*3);
        I->UndoState[I->UndoIter]=-1;
        FreeP(I->UndoCoord[I->UndoIter]);
        if(cs->fInvalidateRep)
          cs->fInvalidateRep(cs,cRepAll,cRepInvCoord);
        SceneChanged();
      }
    }
  }
}
/*========================================================================*/
int ObjectMoleculeAddBond(ObjectMolecule *I,int sele0,int sele1,int order)
{
  int a1,a2;
  AtomInfoType *ai1,*ai2;
  int s1,s2;
  int c = 0;
  BondType *bnd;

  ai1=I->AtomInfo;
  for(a1=0;a1<I->NAtom;a1++) {
    s1=ai1->selEntry;
    if(SelectorIsMember(s1,sele0)) {
      ai2=I->AtomInfo;
      for(a2=0;a2<I->NAtom;a2++) {
        s2=ai2->selEntry;
        if(SelectorIsMember(s2,sele1)) {
          {
            VLACheck(I->Bond,BondType,I->NBond);
            bnd = I->Bond+(I->NBond);
            bnd->index[0]=a1;
            bnd->index[1]=a2;                      
            bnd->order=order;
            bnd->stereo = 0;
            bnd->id=-1;
            I->NBond++;
            c++;
            I->AtomInfo[a1].chemFlag=false;
            I->AtomInfo[a2].chemFlag=false;
          }
        }
        ai2++;
      }
    }
    ai1++;
  }
  if(c) {
    ObjectMoleculeInvalidate(I,cRepLine,cRepInvBonds);
    ObjectMoleculeInvalidate(I,cRepCyl,cRepInvBonds);
    ObjectMoleculeInvalidate(I,cRepNonbonded,cRepInvBonds);
    ObjectMoleculeInvalidate(I,cRepNonbondedSphere,cRepInvBonds);
    ObjectMoleculeInvalidate(I,cRepRibbon,cRepInvBonds);
    ObjectMoleculeInvalidate(I,cRepCartoon,cRepInvBonds);
    ObjectMoleculeUpdateIDNumbers(I);
  }
  return(c);    
}
/*========================================================================*/
int ObjectMoleculeAdjustBonds(ObjectMolecule *I,int sele0,int sele1,int mode,int order)
{
  int a0,a1;
  int offset=0;
  BondType *b0;
  int both;
  int s;
  int a;

  offset=0;
  b0=I->Bond;
  for(a=0;a<I->NBond;a++) {
    a0=b0->index[0];
    a1=b0->index[1];
    
    both=0;
    s=I->AtomInfo[a0].selEntry;
    if(SelectorIsMember(s,sele0))
      both++;
    s=I->AtomInfo[a1].selEntry;
    if(SelectorIsMember(s,sele1))
      both++;
    if(both<2) { /* reverse combo */
      both=0;
      s=I->AtomInfo[a1].selEntry;
      if(SelectorIsMember(s,sele0))
        both++;
      s=I->AtomInfo[a0].selEntry;
      if(SelectorIsMember(s,sele1))
        both++;
    }

    if(both==2) {
      switch(mode) {
      case 0: /* cycle */
        b0->order++;
        if(b0->order>3)
          b0->order=1;
        I->AtomInfo[a0].chemFlag=false;
        I->AtomInfo[a1].chemFlag=false;
        break;
      case 1: /* set */
        b0->order=order;
        I->AtomInfo[a0].chemFlag=false;
        I->AtomInfo[a1].chemFlag=false;
        break;
      }
      ObjectMoleculeInvalidate(I,cRepLine,cRepInvBonds);
      ObjectMoleculeInvalidate(I,cRepCyl,cRepInvBonds);
      ObjectMoleculeInvalidate(I,cRepNonbonded,cRepInvBonds);
      ObjectMoleculeInvalidate(I,cRepNonbondedSphere,cRepInvBonds);
      ObjectMoleculeInvalidate(I,cRepRibbon,cRepInvBonds);
      ObjectMoleculeInvalidate(I,cRepCartoon,cRepInvBonds);
    }
    b0++;
  }
  return(-offset);
}
/*========================================================================*/
int ObjectMoleculeRemoveBonds(ObjectMolecule *I,int sele0,int sele1)
{
  int a0,a1;
  int offset=0;
  BondType *b0,*b1;
  int both;
  int s;
  int a;

  offset=0;
  b0=I->Bond;
  b1=I->Bond;
  for(a=0;a<I->NBond;a++) {
    a0=b0->index[0];
    a1=b0->index[1];
    
    both=0;
    s=I->AtomInfo[a0].selEntry;
    if(SelectorIsMember(s,sele0))
      both++;
    s=I->AtomInfo[a1].selEntry;
    if(SelectorIsMember(s,sele1))
      both++;
    if(both<2) { /* reverse combo */
      both=0;
      s=I->AtomInfo[a1].selEntry;
      if(SelectorIsMember(s,sele0))
        both++;
      s=I->AtomInfo[a0].selEntry;
      if(SelectorIsMember(s,sele1))
        both++;
    }
    
    if(both==2) {
      offset--;
      b0++;
      I->AtomInfo[a0].chemFlag=false;
      I->AtomInfo[a1].chemFlag=false;
    } else if(offset) {
      *(b1++)=*(b0++); /* copy bond info */
    } else {
      *(b1++)=*(b0++); /* copy bond info */
    }
  }
  if(offset) {
    I->NBond += offset;
    VLASize(I->Bond,BondType,I->NBond);
    ObjectMoleculeInvalidate(I,cRepLine,cRepInvBonds);
    ObjectMoleculeInvalidate(I,cRepCyl,cRepInvBonds);
    ObjectMoleculeInvalidate(I,cRepNonbonded,cRepInvBonds);
    ObjectMoleculeInvalidate(I,cRepNonbondedSphere,cRepInvBonds);
    ObjectMoleculeInvalidate(I,cRepRibbon,cRepInvBonds);
    ObjectMoleculeInvalidate(I,cRepCartoon,cRepInvBonds);
  }

  return(-offset);
}
/*========================================================================*/
void ObjectMoleculePurge(ObjectMolecule *I)
{
  int a,a0,a1;
  int *oldToNew = NULL;
  int offset=0;
  BondType *b0,*b1;
  AtomInfoType *ai0,*ai1;
  
  PRINTFD(FB_ObjectMolecule)
    " ObjMolPurge-Debug: step 1, delete object selection\n"
    ENDFD;

  SelectorDelete(I->Obj.Name); /* remove the object selection and free up any selection entries*/
  /* note that we don't delete atom selection members -- those may be needed in the new object */

  PRINTFD(FB_ObjectMolecule)
    " ObjMolPurge-Debug: step 2, purge coordinate sets\n"
    ENDFD;

  for(a=0;a<I->NCSet;a++)
	 if(I->CSet[a]) 
      CoordSetPurge(I->CSet[a]);
  if(I->CSTmpl) {
    CoordSetPurge(I->CSTmpl);
  }
  PRINTFD(FB_ObjectMolecule)
    " ObjMolPurge-Debug: step 3, old-to-new mapping\n"
    ENDFD;

  oldToNew = Alloc(int,I->NAtom);
  ai0=I->AtomInfo;
  ai1=I->AtomInfo;
  for(a=0;a<I->NAtom;a++) {
    if(ai0->deleteFlag) {
      offset--;
      ai0++;
      oldToNew[a]=-1;
    } else if(offset) {
      *(ai1++)=*(ai0++);
      oldToNew[a]=a+offset;
    } else {
      oldToNew[a]=a;
      ai0++;
      ai1++;
    }
  }

  if(offset) {
    I->NAtom += offset;
    VLASize(I->AtomInfo,AtomInfoType,I->NAtom);
    for(a=0;a<I->NCSet;a++)
      if(I->CSet[a])
        CoordSetAdjustAtmIdx(I->CSet[a],oldToNew,I->NAtom);
  }

  PRINTFD(FB_ObjectMolecule)
    " ObjMolPurge-Debug: step 4, bonds\n"
    ENDFD;
  
  offset=0;
  b0=I->Bond;
  b1=I->Bond;
  for(a=0;a<I->NBond;a++) {
    a0=b0->index[0];
    a1=b0->index[1];
    if((oldToNew[a0]<0)||(oldToNew[a1]<0)) {
      offset--;
      b0++;
    } else if(offset) {
      *b1=*b0;
      b1->index[0]=oldToNew[a0]; /* copy bond info */
      b1->index[1]=oldToNew[a1];
      b0++;
      b1++;
    } else {
      *b1=*b0;
      b1->index[0]=oldToNew[a0]; /* copy bond info */
      b1->index[1]=oldToNew[a1];
      b0++;
      b1++; /* TODO check reasoning agaist above */
    }
  }
  if(offset) {
    I->NBond += offset;
    VLASize(I->Bond,BondType,I->NBond);
  }
  FreeP(oldToNew);

  PRINTFD(FB_ObjectMolecule)
    " ObjMolPurge-Debug: step 5, invalidate...\n"
    ENDFD;

  ObjectMoleculeInvalidate(I,cRepAll,cRepInvAtoms);

  PRINTFD(FB_ObjectMolecule)
    " ObjMolPurge-Debug: leaving...\n"
    ENDFD;

}
/*========================================================================*/
int ObjectMoleculeGetAtomGeometry(ObjectMolecule *I,int state,int at)
{
  /* this determines hybridization from coordinates in those few cases
   * where it is unambiguous */

  int result = -1;
  int n,nn;
  float v0[3],v1[3],v2[3],v3[3];
  float d1[3],d2[3],d3[3];
  float cp1[3],cp2[3],cp3[3];
  float avg;
  float dp;
  n  = I->Neighbor[at];
  nn = I->Neighbor[n++]; /* get count */
  if(nn==4) 
    result = cAtomInfoTetrahedral; 
  else if(nn==3) {
    /* check cross products */
    ObjectMoleculeGetAtomVertex(I,state,at,v0);    
    ObjectMoleculeGetAtomVertex(I,state,I->Neighbor[n],v1);
    ObjectMoleculeGetAtomVertex(I,state,I->Neighbor[n+2],v2);
    ObjectMoleculeGetAtomVertex(I,state,I->Neighbor[n+4],v3);
    subtract3f(v1,v0,d1);
    subtract3f(v2,v0,d2);
    subtract3f(v3,v0,d3);
    cross_product3f(d1,d2,cp1);
    cross_product3f(d2,d3,cp2);
    cross_product3f(d3,d1,cp3);
    normalize3f(cp1);
    normalize3f(cp2);
    normalize3f(cp3);
    avg=(dot_product3f(cp1,cp2)+
         dot_product3f(cp2,cp3)+
         dot_product3f(cp3,cp1))/3.0F;
    if(avg>0.75)
      result=cAtomInfoPlaner;
    else
      result=cAtomInfoTetrahedral;
  } else if(nn==2) {
    ObjectMoleculeGetAtomVertex(I,state,at,v0);    
    ObjectMoleculeGetAtomVertex(I,state,I->Neighbor[n],v1);
    ObjectMoleculeGetAtomVertex(I,state,I->Neighbor[n+2],v2);
    subtract3f(v1,v0,d1);
    subtract3f(v2,v0,d2);
    normalize3f(d1);
    normalize3f(d2);
    dp = dot_product3f(d1,d2);
    if(dp<-0.75)
      result=cAtomInfoLinear;
  }
  return(result);
}
/*========================================================================*/
void ObjectMoleculeInferChemForProtein(ObjectMolecule *I,int state)
{
  /* Infers chemical relations for a molecules under protein assumptions.
   * 
   * NOTE: this routine needs an all-atom model (with hydrogens!)
   * and it will make mistakes on non-protein atoms (if they haven't
   * already been assigned)
  */

  int a,n,a0,a1,nn;
  int changedFlag = true;
  
  AtomInfoType *ai,*ai0,*ai1=NULL;
  
  ObjectMoleculeUpdateNeighbors(I);

  /* first, try to find all amids and acids */
  while(changedFlag) {
    changedFlag=false;
    for(a=0;a<I->NAtom;a++) {
      ai=I->AtomInfo+a;
      if(ai->chemFlag) {
        if(ai->geom==cAtomInfoPlaner)
          if(ai->protons == cAN_C) {
            n = I->Neighbor[a];
            nn = I->Neighbor[n++];
            if(nn>1) {
              a1 = -1;
              while(1) {
                a0 = I->Neighbor[n];
                n+=2;
                if(a0<0) break;
                ai0 = I->AtomInfo+a0;
                if((ai0->protons==cAN_O)&&(!ai0->chemFlag)) {
                  a1=a0;
                  ai1=ai0; /* found candidate carbonyl */
                  break;
                }
              }
              if(a1>0) {
                n = I->Neighbor[a]+1;
                while(1) {
                  a0 = I->Neighbor[n];
                  if(a0<0) break;
                  n+=2;
                  if(a0!=a1) {
                    ai0 = I->AtomInfo+a0;
                    if(ai0->protons==cAN_O) {
                      if(!ai0->chemFlag) {
                        ai0->chemFlag=true; /* acid */
                        ai0->geom=cAtomInfoPlaner;
                        ai0->valence=1;
                        ai1->chemFlag=true;
                        ai1->geom=cAtomInfoPlaner;
                        ai1->valence=1;
                        changedFlag=true;
                        break;
                      }
                    } else if(ai0->protons==cAN_N) {
                      if(!ai0->chemFlag) { 
                        ai0->chemFlag=true; /* amide N */ 
                        ai0->geom=cAtomInfoPlaner;                            
                        ai0->valence=3;
                        ai1->chemFlag=true; /* amide =O */ 
                        ai1->geom=cAtomInfoPlaner;
                        ai1->valence=1;
                        changedFlag=true;
                        break;
                      } else if(ai0->geom==cAtomInfoPlaner) {
                        ai1->chemFlag=true; /* amide =O */
                        ai1->geom=cAtomInfoPlaner;
                        ai1->valence=1;
                        changedFlag=true;
                        break;
                      }
                    }
                  }
                }
              }
            }
          }
      }
    }
  }
  /* then handle aldehydes and amines (partial amides - both missing a valence) */
  
  changedFlag=true;
  while(changedFlag) {
    changedFlag=false;
    for(a=0;a<I->NAtom;a++) {
      ai=I->AtomInfo+a;
      if(!ai->chemFlag) {
        if(ai->protons==cAN_C) {
          n = I->Neighbor[a];
          nn = I->Neighbor[n++];
          if(nn>1) {
            a1 = -1;
            while(1) {
              a0 = I->Neighbor[n];
              n+=2;
              if(a0<0) break;
              ai0 = I->AtomInfo+a0;
              if((ai0->protons==cAN_O)&&(!ai0->chemFlag)) { /* =O */
                ai->chemFlag=true; 
                ai->geom=cAtomInfoPlaner;
                ai->valence=1;
                ai0->chemFlag=true;
                ai0->geom=cAtomInfoPlaner;
                ai0->valence=3;
                changedFlag=true;
                break;
              }
            }
          }
        }
        else if(ai->protons==cAN_N)
          {
            if((!ai->chemFlag)||ai->geom!=cAtomInfoLinear) {
              if(ai->formalCharge==0.0) {
                ai->chemFlag=true; 
                ai->geom=cAtomInfoPlaner;
                ai->valence=3;
              }
            }
          }
      }
    }
  }

}
/*========================================================================*/
void ObjectMoleculeInferChemFromNeighGeom(ObjectMolecule *I,int state)
{
  /* infers chemical relations from neighbors and geometry 
  * NOTE: very limited in scope */

  int a,n,a0,nn;
  int changedFlag=true;
  int geom;
  int carbonVal[10];
  
  AtomInfoType *ai,*ai2;

  carbonVal[cAtomInfoTetrahedral] = 4;
  carbonVal[cAtomInfoPlaner] = 3;
  carbonVal[cAtomInfoLinear] = 2;
  
  ObjectMoleculeUpdateNeighbors(I);
  while(changedFlag) {
    changedFlag=false;
    for(a=0;a<I->NAtom;a++) {
      ai=I->AtomInfo+a;
      if(!ai->chemFlag) {
        geom=ObjectMoleculeGetAtomGeometry(I,state,a);
        switch(ai->protons) {
        case cAN_K:
          ai->chemFlag=1;
          ai->geom=cAtomInfoNone;
          ai->valence=0;
          break;
        case cAN_H:
        case cAN_F:
        case cAN_I:
        case cAN_Br:
          ai->chemFlag=1;
          ai->geom=cAtomInfoSingle;
          ai->valence=1;
          break;
        case cAN_O:
          n = I->Neighbor[a];
          nn = I->Neighbor[n++];
          if(nn!=1) { /* water, hydroxy, ether */
            ai->chemFlag=1;
            ai->geom=cAtomInfoTetrahedral;
            ai->valence=2;
          } else { /* hydroxy or carbonyl? check carbon geometry */
            a0 = I->Neighbor[n+2];
            ai2=I->AtomInfo+a0;
            if(ai2->chemFlag) {
              if((ai2->geom==cAtomInfoTetrahedral)||
                 (ai2->geom==cAtomInfoLinear)) {
                ai->chemFlag=1; /* hydroxy */
                ai->geom=cAtomInfoTetrahedral;
                ai->valence=2;
              }
            }
          }
          break;
        case cAN_C:
          if(geom>=0) {
            ai->geom = geom;
            ai->valence = carbonVal[geom];
            ai->chemFlag=true;
          } else {
            n = I->Neighbor[a];
            nn = I->Neighbor[n++];
            if(nn==1) { /* only one neighbor */
              ai2=I->AtomInfo+I->Neighbor[n];
              if(ai2->chemFlag&&(ai2->geom==cAtomInfoTetrahedral)) {
                ai->chemFlag=true; /* singleton carbon bonded to tetC must be tetC */
                ai->geom=cAtomInfoTetrahedral;
                ai->valence=4;
              }
            }
          }
          break;
        case cAN_N:
          if(geom==cAtomInfoPlaner) {
            ai->chemFlag=true;
            ai->geom=cAtomInfoPlaner;
            ai->valence=3;
          } else if(geom==cAtomInfoTetrahedral) {
            ai->chemFlag=true;
            ai->geom=cAtomInfoTetrahedral;
            ai->valence=4;
          }
          break;
        case cAN_S:
          n = I->Neighbor[a];
          nn = I->Neighbor[n++];
          if(nn==4) { /* sulfone */
            ai->chemFlag=true;
            ai->geom=cAtomInfoTetrahedral;
            ai->valence=4;
          } else if(nn==3) { /* suloxide */
            ai->chemFlag=true;
            ai->geom=cAtomInfoTetrahedral;
            ai->valence=3;
          } else if(nn==2) { /* thioether */
            ai->chemFlag=true;
            ai->geom=cAtomInfoTetrahedral;
            ai->valence=2;
          }
          break;
        case cAN_Cl:
          ai->chemFlag=1;
          if(ai->formalCharge==0.0) {
            ai->geom=cAtomInfoSingle;
            ai->valence=1;
          } else {
            ai->geom=cAtomInfoNone;
            ai->valence=0;
          }
          break;
        }
        if(ai->chemFlag)
          changedFlag=true;
      }
    }
  }
}

/*========================================================================*/
void ObjectMoleculeInferHBondFromChem(ObjectMolecule *I)
{
  int a;
  AtomInfoType *ai;
  int a1;
  int n,nn;
  int has_hydro,may_have_lone_pair;
  /* initialize accumulators on uncategorized atoms */

  ObjectMoleculeUpdateNeighbors(I);
  ai=I->AtomInfo;
  for(a=0;a<I->NAtom;a++) {
    n = I->Neighbor[a];
    nn = I->Neighbor[n++];
    ai->hb_donor = false;
    ai->hb_acceptor = false;

    has_hydro = (nn < ai->valence); /* implicit hydrogens? */

    if(!has_hydro) {
      /* explicit hydrogens? */
      has_hydro = false;
      switch(ai->protons) {
      case cAN_N:
      case cAN_O:
        while((a1 = I->Neighbor[n])>=0) {
          n+=2; 
          if(I->AtomInfo[a1].protons==1) {
            has_hydro=true;
            break;
          }
        }
        break;
      }
    }
    
    switch(ai->protons) {
    case cAN_N:
      if(has_hydro)
        ai->hb_donor=true;
      else {
        may_have_lone_pair = false;
        
        if((!has_hydro)&&(ai->protons==cAN_N)) {
          n = I->Neighbor[a] + 1;
          while(I->Neighbor[n]>=0) {
            if(I->Neighbor[n+1]>1) { /* any double/triple/delocalized bonds? */
              may_have_lone_pair = true;
            }
            n+=2; 
          }
        }
        if((ai->formalCharge<=0)&&may_have_lone_pair) {
          ai->hb_acceptor = true;
        }
      }
      break;
    case cAN_O:
      if(has_hydro)
        ai->hb_donor = true;
      if(ai->formalCharge<=0)
        ai->hb_acceptor = true;
      break;
    }
    ai++;
  }
  
}

/*========================================================================*/
void ObjectMoleculeInferChemFromBonds(ObjectMolecule *I,int state)
{

  int a,b;
  BondType *b0;
  AtomInfoType *ai,*ai0,*ai1=NULL;
  int a0,a1;
  int expect,order;
  int n,nn;
  int changedFlag;
  /* initialize accumulators on uncategorized atoms */

  ObjectMoleculeUpdateNeighbors(I);
  ai=I->AtomInfo;
  for(a=0;a<I->NAtom;a++) {
    if(!ai->chemFlag) {
      ai->geom=0;
      ai->valence=0;
    }
    ai++;
  }
  
  /* find maximum bond order for each atom */

  b0=I->Bond;
  for(b=0;b<I->NBond;b++) {
    a0 = b0->index[0];
    a1 = b0->index[1];
    ai0=I->AtomInfo + a0;
    ai1=I->AtomInfo + a1;
    order = b0->order;
    b0++;
    if(!ai0->chemFlag) {
      if(order>ai0->geom)
        ai0->geom=order;
      ai0->valence+=order;
    }
    if(!ai1->chemFlag) {
      if(order>ai1->geom)
        ai1->geom=order;
      ai1->valence+=order;
    }
    if(order==3) { 
      /* override existing chemistry * this is a temp fix to a pressing problem...
         we need to rethink the chemisty assignment ordering (should bond
         information come first? */
      ai0->geom = cAtomInfoLinear;
      ai1->geom = cAtomInfoLinear;
      switch(ai0->protons) {
      case cAN_C:
        ai0->valence=2;
        break;
      default:
        ai0->valence=1;
      }
      switch(ai1->protons) {
      case cAN_C:
        ai1->valence=2;
        break;
      default:
        ai1->valence=1;
      }
      ai0->chemFlag=true;
      ai1->chemFlag=true;
    }
  }

  /* now set up valences and geometries */

  ai=I->AtomInfo;
  for(a=0;a<I->NAtom;a++) {
    if(!ai->chemFlag) {
      expect = AtomInfoGetExpectedValence(ai);
      n = I->Neighbor[a];
      nn = I->Neighbor[n++];
      if(ai->geom==3) {
        ai->geom = cAtomInfoLinear;
        switch(ai->protons) {
        case cAN_C:
          ai->valence=2;
          break;
        default:
          ai->valence=1;
        }
        ai->chemFlag=true;
      } else {
      if(expect<0) 
        expect = -expect; /* for now, just ignore this issue */
      /*      printf("%d %d %d %d\n",ai->geom,ai->valence,nn,expect);*/
      if(ai->valence==expect) { /* sum of bond orders equals valence */
        ai->chemFlag=true;
        ai->valence=nn;
        switch(ai->geom) { /* max bond order observed */
        case 0: ai->geom = cAtomInfoNone; break;
        case 2: ai->geom = cAtomInfoPlaner; break;
        case 3: ai->geom = cAtomInfoLinear; break;
        default: 
          if(expect==1) 
            ai->geom = cAtomInfoSingle;
          else
            ai->geom = cAtomInfoTetrahedral; 
          break;            
        }
      } else if(ai->valence<expect) { /* missing a bond */
        ai->chemFlag=true;
        ai->valence=nn+(expect-ai->valence); 
        switch(ai->geom) { 
        case 2: ai->geom = cAtomInfoPlaner; break;
        case 3: ai->geom = cAtomInfoLinear; break;
        default: 
          if(expect==1) 
            ai->geom = cAtomInfoSingle;
          else
            ai->geom = cAtomInfoTetrahedral; 
          break;
        }
      } else if(ai->valence>expect) {
        ai->chemFlag=true;
        ai->valence=nn;
        switch(ai->geom) { 
        case 2: ai->geom = cAtomInfoPlaner; break;
        case 3: ai->geom = cAtomInfoLinear; break;
        default: 
          if(expect==1) 
            ai->geom = cAtomInfoSingle;
          else
            ai->geom = cAtomInfoTetrahedral; 
          break;
        }
        if(nn>3)
          ai->geom = cAtomInfoTetrahedral;
      }
      }
    }
    ai++;
  }

  /* now go through and make sure conjugated amines are planer */
  changedFlag=true;
  while(changedFlag) {
    changedFlag=false;
    ai=I->AtomInfo;
    for(a=0;a<I->NAtom;a++) {
      if(ai->chemFlag) {
        if(ai->protons==cAN_N) 
          if(ai->formalCharge==0) 
            if(ai->geom==cAtomInfoTetrahedral) { 
              /* search for uncharged tetrahedral nitrogen */
              n = I->Neighbor[a]+1;
              while(1) {
                a0 = I->Neighbor[n];
                n+=2;
                if(a0<0) break;
                ai0 = I->AtomInfo+a0;
                if((ai0->chemFlag)&&(ai0->geom==cAtomInfoPlaner)&&
                   ((ai0->protons==cAN_C)||(ai0->protons==cAN_N))) {
                  ai->geom=cAtomInfoPlaner; /* found probable delocalization */
                  ai->valence=3; /* just in case...*/
                  changedFlag=true;
                  break;
                }
              }
            }
      }
      ai++; 
    }
  }

  /* now go through and make sure conjugated anions are planer */
  changedFlag=true;
  while(changedFlag) {
    changedFlag=false;
    ai=I->AtomInfo;
    for(a=0;a<I->NAtom;a++) {
      if(ai->chemFlag) {
        if(ai->protons==cAN_O) 
          if(ai->formalCharge==-1) 
            if((ai->geom==cAtomInfoTetrahedral)||
               (ai->geom==cAtomInfoSingle)) { /* search for anionic tetrahedral oxygen */
              n = I->Neighbor[a]+1;
              while(1) {
                a0 = I->Neighbor[n];
                n+=2;
                if(a0<0) break;
                ai0 = I->AtomInfo+a0;
                if((ai0->chemFlag)&&(ai0->geom==cAtomInfoPlaner)&&
                   ((ai0->protons==cAN_C)||(ai0->protons==cAN_N))) {
                  ai->geom=cAtomInfoPlaner; /* found probable delocalization */
                  changedFlag=true;
                  break;
                }
              }
            }
      }
      ai++;
    }
  }


}
/*========================================================================*/
int ObjectMoleculeTransformSelection(ObjectMolecule *I,int state,
                                      int sele,float *TTT,int log,char *sname) 
{
  /* if sele == -1, then the whole object state is transformed */

  int a,s;
  int flag=false;
  CoordSet *cs;
  AtomInfoType *ai;
  int logging;
  int all_states=false,inp_state;
  int ok=true;
  inp_state=state;
  if(state==-1) 
    state=ObjectGetCurrentState(&I->Obj,false);
  if(state<0) {
    all_states=true;
    state=-1;
  }
  PRINTFD(FB_ObjectMolecule)
    "ObjMolTransSele-Debug: state %d\n",state
    ENDFD;
  while(1) {
    if(all_states) {
      state++;
      if(state>=I->NCSet)
        break;
    }
    if(state<I->NCSet) {
      cs = I->CSet[state];
      if(cs) {
        if(sele>=0) {
          ai=I->AtomInfo;
          for(a=0;a<I->NAtom;a++) {
            s=ai->selEntry;
            if(!(ai->protekted==1))
              if(SelectorIsMember(s,sele))
                {
                  CoordSetTransformAtom(cs,a,TTT);
                  flag=true;
                }
            ai++;
          }
        } else {
          ai=I->AtomInfo;
          for(a=0;a<I->NAtom;a++) {
            if(!(ai->protekted==1))
              CoordSetTransformAtom(cs,a,TTT);
            ai++;
          }
          flag=true;
        }
        if(flag) 
          cs->fInvalidateRep(cs,cRepAll,cRepInvCoord);
      }
    }
    if(!all_states)
      break;
  }

  
  if(log) {
    OrthoLineType line;
    WordType sele_str = ",'";
    logging = (int)SettingGet(cSetting_logging);
    if(sele>=0) {
      strcat(sele_str,sname);
      strcat(sele_str,"'");
    }
    else
      sele_str[0]=0;
    switch(logging) {
    case cPLog_pml:
      sprintf(line,
              "_ cmd.transform_object('%s',[\\\n_ %15.9f,%15.9f,%15.9f,%15.9f,\\\n_ %15.9f,%15.9f,%15.9f,%15.9f,\\\n_ %15.9f,%15.9f,%15.9f,%15.9f,\\\n_ %15.9f,%15.9f,%15.9f,%15.9f\\\n_     ],%d,%d%s)\n",
              I->Obj.Name,
              TTT[ 0],TTT[ 1],TTT[ 2],TTT[ 3],
              TTT[ 4],TTT[ 5],TTT[ 6],TTT[ 7],
              TTT[ 8],TTT[ 9],TTT[10],TTT[11],
              TTT[12],TTT[13],TTT[14],TTT[15],inp_state+1,log,sele_str);
      PLog(line,cPLog_no_flush);
      break;
    case cPLog_pym:
      
      sprintf(line,
              "cmd.transform_object('%s',[\n%15.9f,%15.9f,%15.9f,%15.9f,\n%15.9f,%15.9f,%15.9f,%15.9f,\n%15.9f,%15.9f,%15.9f,%15.9f,\n%15.9f,%15.9f,%15.9f,%15.9f\n],%d,%d%s)\n",
              I->Obj.Name,
              TTT[ 0],TTT[ 1],TTT[ 2],TTT[ 3],
              TTT[ 4],TTT[ 5],TTT[ 6],TTT[ 7],
              TTT[ 8],TTT[ 9],TTT[10],TTT[11],
              TTT[12],TTT[13],TTT[14],TTT[15],inp_state+1,log,sele_str);
      PLog(line,cPLog_no_flush);
      break;
    default:
      break;
    }
  }
  return(ok);
}
/*========================================================================*/
int ObjectMoleculeGetAtomIndex(ObjectMolecule *I,int sele)
{
  int a,s;
  if(sele<0)
    return(-1);
  for(a=0;a<I->NAtom;a++) {
    s=I->AtomInfo[a].selEntry;
    if(SelectorIsMember(s,sele))
      return(a);
  }
  return(-1);
}
/*========================================================================*/
void ObjectMoleculeUpdateNonbonded(ObjectMolecule *I)
{
  int a;
  BondType *b;
  AtomInfoType *ai;

  if(!I->DiscreteFlag) {
    ai=I->AtomInfo;
    
    for(a=0;a<I->NAtom;a++)
      (ai++)->bonded = false;
    
    b=I->Bond;
    ai=I->AtomInfo;
    for(a=0;a<I->NBond;a++)
      {
        ai[b->index[0]].bonded=true;
        ai[b->index[1]].bonded=true;
        b++;
      }

  }
}
/*========================================================================*/
int ObjectMoleculeGetTotalAtomValence(ObjectMolecule *I,int atom)
{
  int result = 0;
  int n0;
  ObjectMoleculeUpdateNeighbors(I);
  if(atom<I->NAtom) {
    n0 = I->Neighbor[atom]+1;
    while(I->Neighbor[n0]>=0) {
      result+=I->Neighbor[n0+1];
    }
    n0+=2;
  } else {
    result = -1; /* error */
  }
  return result;
}
/*========================================================================*/
void ObjectMoleculeUpdateNeighbors(ObjectMolecule *I)
{
  /* neighbor storage structure: VERY COMPLICATED...
     
     0       list offset for atom 0 = n
     1       list offset for atom 1 = n + m + 1
     ...
     n-1     list offset for atom n-1

     n       count for atom 0 
     n+1     neighbor of atom 0
     n+2     bond index
     n+3     neighbor of atom 0
     n+4     bond index
     ...
     n+m     -1 terminator for atom 0

     n+m+1   count for atom 1
     n+m+2   neighbor of atom 1
     n+m+3   bond index
     n+m+4   neighbor of atom 1
     n+m+5   bond index
     etc.

     NOTE: all atoms have an offset and a terminator whether or not they have any bonds 
 */

  int size;
  int a,b,c,d,l0,l1,*l;
  BondType *bnd;
  if(!I->Neighbor) {
    size = (I->NAtom*3)+(I->NBond*4); 
    if(I->Neighbor) {
      VLACheck(I->Neighbor,int,size);
    } else {
      I->Neighbor=VLAlloc(int,size);
    }
    
    /* initialize */
    l = I->Neighbor;
    for(a=0;a<I->NAtom;a++)
      (*l++)=0;
    
    /* count neighbors for each atom */
    bnd = I->Bond;
    for(b=0;b<I->NBond;b++) {
      I->Neighbor[bnd->index[0]]++;
      I->Neighbor[bnd->index[1]]++;
      bnd++;
    }
    
    /* set up offsets and list terminators */
    c = I-> NAtom;
    for(a=0;a<I->NAtom;a++) {
      d = I->Neighbor[a]; /* get number of neighbors */
      I->Neighbor[c]=d; /* store neighbor count */
      I->Neighbor[a]=c+d+d+1; /* set initial position to end of list, we'll fill backwards */
      I->Neighbor[I->Neighbor[a]]=-1; /* store terminator */
      c += d + d + 2;
    }
    
    /* now load neighbors in a sequential list for each atom (reverse order) */
    bnd = I->Bond;
    for(b=0;b<I->NBond;b++) {
      l0 = bnd->index[0];
      l1 = bnd->index[1];
      bnd++;

      I->Neighbor[l0]--; 
      I->Neighbor[I->Neighbor[l0]]=b; /* store bond indices (for I->Bond) */
      I->Neighbor[l0]--;
      I->Neighbor[I->Neighbor[l0]]=l1; /* store neighbor references (I->AtomInfo, etc.)*/

      I->Neighbor[l1]--;
      I->Neighbor[I->Neighbor[l1]]=b; /* store bond indices (for I->Bond) */
      I->Neighbor[l1]--;
      I->Neighbor[I->Neighbor[l1]]=l0; /* store neighbor references (I->AtomInfo, etc.)*/
    }
    for(a=0;a<I->NAtom;a++) { /* adjust down to point to the count, not the first entry */
      if(I->Neighbor[a]>=0)
        I->Neighbor[a]--;
    }
    l=I->Neighbor;
  }
}
/*========================================================================*/
CoordSet *ObjectMoleculeChemPyModel2CoordSet(PyObject *model,AtomInfoType **atInfoPtr)
{
  int nAtom,nBond;
  int a,c;
  float *coord = NULL;
  CoordSet *cset = NULL;
  AtomInfoType *atInfo = NULL,*ai;
  float *f;
  BondType *ii,*bond=NULL;
  int ok=true;
  int auto_show_lines;
  int auto_show_nonbonded;
  int hetatm;
  int ignore_ids;

  PyObject *atomList = NULL;
  PyObject *bondList = NULL;
  PyObject *atom = NULL;
  PyObject *bnd = NULL;
  PyObject *index = NULL;
  PyObject *crd = NULL;
  PyObject *tmp = NULL;
  auto_show_lines = (int)SettingGet(cSetting_auto_show_lines);
  auto_show_nonbonded = (int)SettingGet(cSetting_auto_show_nonbonded);

  ignore_ids=!(int)SettingGet(cSetting_preserve_chempy_ids);

  nAtom=0;
  nBond=0;
  if(atInfoPtr)
	 atInfo = *atInfoPtr;

  atomList = PyObject_GetAttrString(model,"atom");
  if(atomList) 
    nAtom = PyList_Size(atomList);
  else 
    ok=ErrMessage("ObjectMoleculeChemPyModel2CoordSet","can't get atom list");


  if(ok) {
	 coord=VLAlloc(float,3*nAtom);
	 if(atInfo)
		VLACheck(atInfo,AtomInfoType,nAtom);	 
  }

  if(ok) { 
    
	 f=coord;
	 for(a=0;a<nAtom;a++)
		{
        atom = PyList_GetItem(atomList,a);
        if(!atom) 
          ok=ErrMessage("ObjectMoleculeChemPyModel2CoordSet","can't get atom");
        crd = PyObject_GetAttrString(atom,"coord");
        if(!crd) 
          ok=ErrMessage("ObjectMoleculeChemPyModel2CoordSet","can't get coordinates");
        else {
          for(c=0;c<3;c++) {
            tmp = PyList_GetItem(crd,c);
            if (tmp) 
              ok = PConvPyObjectToFloat(tmp,f++);
            if(!ok) {
              ErrMessage("ObjectMoleculeChemPyModel2CoordSet","can't read coordinates");
              break;
            }
          }
        }
        Py_XDECREF(crd);
        
        ai = atInfo+a;
        ai->id = a; /* chempy models are zero-based */
        if(!ignore_ids) { 
          if(ok) { /* get chempy atom id if extant */
            if(PTruthCallStr(atom,"has","id")) { 
              tmp = PyObject_GetAttrString(atom,"id");
              if (tmp)
                ok = PConvPyObjectToInt(tmp,&ai->id);
              if(!ok) 
                ErrMessage("ObjectMoleculeChemPyModel2CoordSet","can't read atom identifier");
              Py_XDECREF(tmp);
            } else {
              ai->id=-1;
            }
          }
        }

        if(ok) {
          tmp = PyObject_GetAttrString(atom,"name");
          if (tmp)
            ok = PConvPyObjectToStrMaxClean(tmp,ai->name,sizeof(AtomName)-1);
          if(!ok) 
            ErrMessage("ObjectMoleculeChemPyModel2CoordSet","can't read name");
          Py_XDECREF(tmp);
        }

        if(ok) {
          if(PTruthCallStr(atom,"has","text_type")) { 
            tmp = PyObject_GetAttrString(atom,"text_type");
            if (tmp)
              ok = PConvPyObjectToStrMaxClean(tmp,ai->textType,sizeof(TextType)-1);
            if(!ok) 
              ErrMessage("ObjectMoleculeChemPyModel2CoordSet","can't read text_type");
            Py_XDECREF(tmp);
          } else {
            ai->textType[0]=0;
          }
        }

        if(ok) {
          if(PTruthCallStr(atom,"has","vdw")) { 
            tmp = PyObject_GetAttrString(atom,"vdw");
            if (tmp)
              ok = PConvPyObjectToFloat(tmp,&ai->vdw);
            if(!ok) 
              ErrMessage("ObjectMoleculeChemPyModel2CoordSet","can't read vdw radius");
            Py_XDECREF(tmp);
          } else {
            ai->vdw=0.0;
          }
        }

        if(ok) {
          if(PTruthCallStr(atom,"has","stereo")) { 
            tmp = PyObject_GetAttrString(atom,"stereo");
            if (tmp)
              ok = PConvPyObjectToInt(tmp,&ai->stereo);
            if(!ok) 
              ErrMessage("ObjectMoleculeChemPyModel2CoordSet","can't read stereo");
            Py_XDECREF(tmp);
          } else {
            ai->stereo = 0;
          }
        }

        if(ok) {
          if(PTruthCallStr(atom,"has","numeric_type")) { 
            tmp = PyObject_GetAttrString(atom,"numeric_type");
            if (tmp)
              ok = PConvPyObjectToInt(tmp,&ai->customType);
            if(!ok) 
              ErrMessage("ObjectMoleculeChemPyModel2CoordSet","can't read numeric_type");
            Py_XDECREF(tmp);
          } else {
            ai->customType = cAtomInfoNoType;
          }
        }

        if(ok) {
          if(PTruthCallStr(atom,"has","formal_charge")) { 
            tmp = PyObject_GetAttrString(atom,"formal_charge");
            if (tmp)
              ok = PConvPyObjectToInt(tmp,&ai->formalCharge);
            if(!ok) 
              ErrMessage("ObjectMoleculeChemPyModel2CoordSet","can't read formal_charge");
            Py_XDECREF(tmp);
          } else {
            ai->formalCharge = 0;
          }
        }

        if(ok) {
          if(PTruthCallStr(atom,"has","partial_charge")) { 
            tmp = PyObject_GetAttrString(atom,"partial_charge");
            if (tmp)
              ok = PConvPyObjectToFloat(tmp,&ai->partialCharge);
            if(!ok) 
              ErrMessage("ObjectMoleculeChemPyModel2CoordSet","can't read partial_charge");
            Py_XDECREF(tmp);
          } else {
            ai->partialCharge = 0.0;
          }
        }

        if(ok) {
          if(PTruthCallStr(atom,"has","flags")) {         
            tmp = PyObject_GetAttrString(atom,"flags");
            if (tmp)
              ok = PConvPyObjectToInt(tmp,(int*)&ai->flags);
            if(!ok) 
              ErrMessage("ObjectMoleculeChemPyModel2CoordSet","can't read flags");
            Py_XDECREF(tmp);
          } else {
            ai->flags = 0;
          }
        }

        if(ok) {
          tmp = PyObject_GetAttrString(atom,"resn");
          if (tmp)
            ok = PConvPyObjectToStrMaxClean(tmp,ai->resn,sizeof(ResName)-1);
          if(!ok) 
            ErrMessage("ObjectMoleculeChemPyModel2CoordSet","can't read resn");
          Py_XDECREF(tmp);
        }
        
		  if(ok) {
          tmp = PyObject_GetAttrString(atom,"resi");
          if (tmp)
            ok = PConvPyObjectToStrMaxClean(tmp,ai->resi,sizeof(ResIdent)-1);
          if(!ok) 
            ErrMessage("ObjectMoleculeChemPyModel2CoordSet","can't read resi");
          else
            sscanf(ai->resi,"%d",&ai->resv);
          Py_XDECREF(tmp);
        }

        if(ok) {
          if(PTruthCallStr(atom,"has","resi_number")) {         
            tmp = PyObject_GetAttrString(atom,"resi_number");
            if (tmp)
              ok = PConvPyObjectToInt(tmp,&ai->resv);
            if(!ok) 
              ErrMessage("ObjectMoleculeChemPyModel2CoordSet","can't read resi_number");
            Py_XDECREF(tmp);
          }
        }
        
		  if(ok) {
          tmp = PyObject_GetAttrString(atom,"segi");
          if (tmp)
            ok = PConvPyObjectToStrMaxClean(tmp,ai->segi,sizeof(SegIdent)-1);
          if(!ok) 
            ErrMessage("ObjectMoleculeChemPyModel2CoordSet","can't read segi");
          Py_XDECREF(tmp);
        }

		  if(ok) {
          tmp = PyObject_GetAttrString(atom,"b");
          if (tmp)
            ok = PConvPyObjectToFloat(tmp,&ai->b);
          if(!ok) 
            ErrMessage("ObjectMoleculeChemPyModel2CoordSet","can't read b value");
          Py_XDECREF(tmp);
        }

		  if(ok) {
          tmp = PyObject_GetAttrString(atom,"q");
          if (tmp)
            ok = PConvPyObjectToFloat(tmp,&ai->q);
          if(!ok) 
            ErrMessage("ObjectMoleculeChemPyModel2CoordSet","can't read occupancy");
          Py_XDECREF(tmp);
        }

        
		  if(ok) {
          tmp = PyObject_GetAttrString(atom,"chain");
          if (tmp)
            ok = PConvPyObjectToStrMaxClean(tmp,ai->chain,sizeof(Chain)-1);
          if(!ok) 
            ErrMessage("ObjectMoleculeChemPyModel2CoordSet","can't read chain");
          Py_XDECREF(tmp);
        }
        
		  if(ok) {
          tmp = PyObject_GetAttrString(atom,"hetatm");
          if (tmp)
            ok = PConvPyObjectToInt(tmp,&hetatm);
          if(!ok) 
            ErrMessage("ObjectMoleculeChemPyModel2CoordSet","can't read hetatm");
          else
            ai->hetatm = hetatm;
          Py_XDECREF(tmp);
        }
        
		  if(ok) {
          tmp = PyObject_GetAttrString(atom,"alt");
          if (tmp)
            ok = PConvPyObjectToStrMaxClean(tmp,ai->alt,sizeof(Chain)-1);
          if(!ok) 
            ErrMessage("ObjectMoleculeChemPyModel2CoordSet","can't read alternate conformation");
          Py_XDECREF(tmp);
        }

		  if(ok) {
          tmp = PyObject_GetAttrString(atom,"symbol");
          if (tmp)
            ok = PConvPyObjectToStrMaxClean(tmp,ai->elem,sizeof(AtomName)-1);
          if(!ok) 
            ErrMessage("ObjectMoleculeChemPyModel2CoordSet","can't read symbol");
          Py_XDECREF(tmp);
        }

		  if(ok) {
          tmp = PyObject_GetAttrString(atom,"ss");
          if (tmp)
            ok = PConvPyObjectToStrMaxClean(tmp,ai->ssType,sizeof(SSType)-1);
          if(!ok) 
            ErrMessage("ObjectMoleculeChemPyModel2CoordSet","can't read secondary structure");
          Py_XDECREF(tmp);
        }

        
        for(c=0;c<cRepCnt;c++) {
          atInfo[a].visRep[c] = false;
		  }
        atInfo[a].visRep[cRepLine] = auto_show_lines; /* show lines by default */
        atInfo[a].visRep[cRepNonbonded] = auto_show_nonbonded; /* show lines by default */

		  if(ok&&atInfo) {
			 AtomInfoAssignParameters(ai);
			 atInfo[a].color=AtomInfoGetColor(ai);
		  }


		  if(!ok)
			 break;
		}
  }

  bondList = PyObject_GetAttrString(model,"bond");
  if(bondList) 
    nBond = PyList_Size(bondList);
  else
    ok=ErrMessage("ObjectMoleculeChemPyModel2CoordSet","can't get bond list");

  if(ok) {
	 bond=VLAlloc(BondType,nBond);
    ii=bond;
	 for(a=0;a<nBond;a++)
		{
        bnd = PyList_GetItem(bondList,a);
        if(!bnd) 
          ok=ErrMessage("ObjectMoleculeChemPyModel2CoordSet","can't get bond");
        index = PyObject_GetAttrString(bnd,"index");
        if(!index) 
          ok=ErrMessage("ObjectMoleculeChemPyModel2CoordSet","can't get bond indices");
        else {
          for(c=0;c<2;c++) {
            tmp = PyList_GetItem(index,c);
            if (tmp) 
              ok = PConvPyObjectToInt(tmp,&ii->index[c]);
            if(!ok) {
              ErrMessage("ObjectMoleculeChemPyModel2CoordSet","can't read coordinates");
              break;
            }
          }
        }
        if(ok) {
          tmp = PyObject_GetAttrString(bnd,"order");
          if (tmp)
            ok = PConvPyObjectToInt(tmp,&ii->order);
          if(!ok) 
            ErrMessage("ObjectMoleculeChemPyModel2CoordSet","can't read bond order");
          Py_XDECREF(tmp);
        }

        if(ok) {
          tmp = PyObject_GetAttrString(bnd,"stereo");
          if (tmp)
            ok = PConvPyObjectToInt(tmp,&ii->stereo);
          else 
            ii->stereo = 0;
          if(!ok) 
            ii->stereo = 0;
          Py_XDECREF(tmp);
        }

        ii->id=a;
        if(!ignore_ids) { 
          if(ok) { /* get unique chempy bond id if present */
            if(PTruthCallStr(bnd,"has","id")) { 
              tmp = PyObject_GetAttrString(bnd,"id");
              if (tmp)
                ok = PConvPyObjectToInt(tmp,&ii->id);
              if(!ok) 
                ErrMessage("ObjectMoleculeChemPyModel2CoordSet","can't read bond identifier");
              Py_XDECREF(tmp);
            } else {
              ii->id=-1;
            }
          }
        }
        Py_XDECREF(index);
        ii++;
      }
  }

  Py_XDECREF(atomList);
  Py_XDECREF(bondList);

  if(ok) {
	 cset = CoordSetNew();
	 cset->NIndex=nAtom;
	 cset->Coord=coord;
	 cset->NTmpBond=nBond;
	 cset->TmpBond=bond;
  } else {
	 VLAFreeP(bond);
	 VLAFreeP(coord);
  }
  if(atInfoPtr)
	 *atInfoPtr = atInfo;

  if(PyErr_Occurred())
    PyErr_Print();
  return(cset);
}


/*========================================================================*/
ObjectMolecule *ObjectMoleculeLoadChemPyModel(ObjectMolecule *I,PyObject *model,int frame,int discrete)
{
  CoordSet *cset = NULL;
  AtomInfoType *atInfo;
  int ok=true;
  int isNew = true;
  unsigned int nAtom = 0;
  PyObject *tmp,*mol;

  if(!I) 
	 isNew=true;
  else 
	 isNew=false;

  if(ok) {

	 if(isNew) {
		I=(ObjectMolecule*)ObjectMoleculeNew(discrete);
		atInfo = I->AtomInfo;
		isNew = true;
	 } else {
		atInfo=VLAMalloc(10,sizeof(AtomInfoType),2,true); /* autozero here is important */
		isNew = false;
	 }

    if(isNew) {
      AtomInfoPrimeColors();
      I->Obj.Color = AtomInfoGetCarbColor();
    }

	 cset=ObjectMoleculeChemPyModel2CoordSet(model,&atInfo);	 

    mol = PyObject_GetAttrString(model,"molecule");
    if(mol) {
      if(PyObject_HasAttrString(mol,"title")) {
        tmp = PyObject_GetAttrString(mol,"title");
        if(tmp) {
          UtilNCopy(cset->Name,PyString_AsString(tmp),sizeof(WordType));
          Py_DECREF(tmp);
          if(!strcmp(cset->Name,"untitled")) /* ignore untitled */
            cset->Name[0]=0;
        }
      }
      Py_DECREF(mol);
    }
    if(PyObject_HasAttrString(model,"spheroid")&&
       PyObject_HasAttrString(model,"spheroid_normals"))
      {
        tmp = PyObject_GetAttrString(model,"spheroid");
        if(tmp) {
          cset->NSpheroid = PConvPyListToFloatArray(tmp,&cset->Spheroid);
          if(cset->NSpheroid<0) cset->NSpheroid=0;
          Py_DECREF(tmp);
        }
        tmp = PyObject_GetAttrString(model,"spheroid_normals");
        if(tmp) {
          PConvPyListToFloatArray(tmp,&cset->SpheroidNormal);
          Py_DECREF(tmp);
        }
      }
    mol = PyObject_GetAttrString(model,"molecule");
    
	 nAtom=cset->NIndex;
  }
  /* include coordinate set */
  if(ok) {
    cset->Obj = I;
    cset->fEnumIndices(cset);
    if(cset->fInvalidateRep)
      cset->fInvalidateRep(cset,cRepAll,cRepInvRep);
    if(isNew) {	
      I->AtomInfo=atInfo; /* IMPORTANT to reassign: this VLA may have moved! */
    } else {
      ObjectMoleculeMerge(I,atInfo,cset,false,cAIC_AllMask); /* NOTE: will release atInfo */
    }
    if(isNew) I->NAtom=nAtom;
    if(frame<0) frame=I->NCSet;
    VLACheck(I->CSet,CoordSet*,frame);
    if(I->NCSet<=frame) I->NCSet=frame+1;
    if(I->CSet[frame]) I->CSet[frame]->fFree(I->CSet[frame]);
    I->CSet[frame] = cset;
    if(isNew) I->NBond = ObjectMoleculeConnect(I,&I->Bond,I->AtomInfo,cset,false);
    if(cset->Symmetry&&(!I->Symmetry)) {
      I->Symmetry=SymmetryCopy(cset->Symmetry);
      SymmetryAttemptGeneration(I->Symmetry,false,false);
    }
    SceneCountFrames();
    ObjectMoleculeExtendIndices(I);
    ObjectMoleculeSort(I);
    ObjectMoleculeUpdateIDNumbers(I);
    ObjectMoleculeUpdateNonbonded(I);
  }
  return(I);
}


/*========================================================================*/
ObjectMolecule *ObjectMoleculeLoadCoords(ObjectMolecule *I,PyObject *coords,int frame)
{
  CoordSet *cset = NULL;
  int ok=true;
  int a,l;
  PyObject *v;
  float *f;
  a=0;
  while(a<I->NCSet) {
    if(I->CSet[a]) {
      cset=I->CSet[a];
      break;
    }
    a++;
  }
  
  if(!PyList_Check(coords)) 
    ErrMessage("LoadsCoords","passed argument is not a list");
  else {
    l = PyList_Size(coords);
    if (l==cset->NIndex) {
      cset=CoordSetCopy(cset);
      f=cset->Coord;
      for(a=0;a<l;a++) {
        v=PyList_GetItem(coords,a);
/* no error checking */
        *(f++)=(float)PyFloat_AsDouble(PyList_GetItem(v,0)); 
        *(f++)=(float)PyFloat_AsDouble(PyList_GetItem(v,1));
        *(f++)=(float)PyFloat_AsDouble(PyList_GetItem(v,2));
      }
    }
  }
  /* include coordinate set */
  if(ok) {
    if(cset->fInvalidateRep)
      cset->fInvalidateRep(cset,cRepAll,cRepInvRep);

    if(frame<0) frame=I->NCSet;
    VLACheck(I->CSet,CoordSet*,frame);
    if(I->NCSet<=frame) I->NCSet=frame+1;
    if(I->CSet[frame]) I->CSet[frame]->fFree(I->CSet[frame]);
    I->CSet[frame] = cset;
    SceneCountFrames();
  }
  return(I);
}

/*========================================================================*/
void ObjectMoleculeBlindSymMovie(ObjectMolecule *I)
{
  CoordSet *frac;
  int a,c;
  int x,y,z;
  float m[16];

  if(I->NCSet!=1) {
    ErrMessage("ObjectMolecule:","SymMovie only works on objects with a single state.");
  } else if(!I->Symmetry) {
    ErrMessage("ObjectMolecule:","No symmetry loaded!");
  } else if(!I->Symmetry->NSymMat) {
    ErrMessage("ObjectMolecule:","No symmetry matrices!");    
  } else if(I->CSet[0]) {
    frac = CoordSetCopy(I->CSet[0]);
    CoordSetRealToFrac(frac,I->Symmetry->Crystal);
    for(x=-1;x<2;x++)
      for(y=-1;y<2;y++)
        for(z=-1;z<2;z++)
          for(a=0;a<I->Symmetry->NSymMat;a++) {
            if(!((!a)&&(!x)&&(!y)&&(!z))) {
              c = I->NCSet;
              VLACheck(I->CSet,CoordSet*,c);
              I->CSet[c] = CoordSetCopy(frac);
              CoordSetTransform44f(I->CSet[c],I->Symmetry->SymMatVLA+(a*16));
              identity44f(m);
              m[3] = (float)x;
              m[7] = (float)y;
              m[11] = (float)z;
              CoordSetTransform44f(I->CSet[c],m);
              CoordSetFracToReal(I->CSet[c],I->Symmetry->Crystal);
              I->NCSet++;
            }
          }
    frac->fFree(frac);
  }
  SceneChanged();
}

/*========================================================================*/
void ObjectMoleculeExtendIndices(ObjectMolecule *I)
{
  int a;
  CoordSet *cs;

  for(a=-1;a<I->NCSet;a++) {
    if(a<0) 
      cs=I->CSTmpl;
    else
      cs=I->CSet[a];
	 if(cs)
      if(cs->fExtendIndices)
        cs->fExtendIndices(cs,I->NAtom);
  }
}
/*========================================================================*/

CoordSet *ObjectMoleculeMOLStr2CoordSet(char *buffer,AtomInfoType **atInfoPtr)
{
  char *p;
  int nAtom,nBond;
  int a,c,cnt,atm,chg;
  float *coord = NULL;
  CoordSet *cset = NULL;
  AtomInfoType *atInfo = NULL;
  char cc[MAXLINELEN],cc1[MAXLINELEN],resn[MAXLINELEN] = "UNK";
  float *f;
  BondType *ii;
  BondType *bond=NULL;
  int ok=true;
  int auto_show_lines;
  int auto_show_nonbonded;
  WordType nameTmp;

  auto_show_lines = (int)SettingGet(cSetting_auto_show_lines);
  auto_show_nonbonded = (int)SettingGet(cSetting_auto_show_nonbonded);

  p=buffer;
  nAtom=0;
  if(atInfoPtr)
	 atInfo = *atInfoPtr;

  /*  p=ParseWordCopy(nameTmp,p,sizeof(WordType)-1);*/
  p=ParseNCopy(nameTmp,p,sizeof(WordType)-1);
  p=nextline(p); 
  p=nextline(p);
  p=nextline(p);

  if(ok) {
	 p=ncopy(cc,p,3);
	 if(sscanf(cc,"%d",&nAtom)!=1)
		ok=ErrMessage("ReadMOLFile","bad atom count");
  }

  if(ok) {  
	 p=ncopy(cc,p,3);
	 if(sscanf(cc,"%d",&nBond)!=1)
		ok=ErrMessage("ReadMOLFile","bad bond count");
  }

  if(ok) {
	 coord=VLAlloc(float,3*nAtom);
	 if(atInfo)
		VLACheck(atInfo,AtomInfoType,nAtom);	 
  }
  
  p=nextline(p);

  /* read coordinates and atom names */

  if(ok) { 
	 f=coord;
	 for(a=0;a<nAtom;a++)
		{
		  if(ok) {
			 p=ncopy(cc,p,10);
			 if(sscanf(cc,"%f",f++)!=1)
				ok=ErrMessage("ReadMOLFile","bad coordinate");
		  }
		  if(ok) {
			 p=ncopy(cc,p,10);
			 if(sscanf(cc,"%f",f++)!=1)
				ok=ErrMessage("ReadMOLFile","bad coordinate");
		  }
		  if(ok) {
			 p=ncopy(cc,p,10);
			 if(sscanf(cc,"%f",f++)!=1)
				ok=ErrMessage("ReadMOLFile","bad coordinate");
		  }
		  if(ok) {
          p=nskip(p,1);
			 p=ncopy(atInfo[a].name,p,3);
			 UtilCleanStr(atInfo[a].name);
			 
          for(c=0;c<cRepCnt;c++) {
            atInfo[a].visRep[c] = false;
          }
          atInfo[a].visRep[cRepLine] = auto_show_lines; /* show lines by default */
          atInfo[a].visRep[cRepNonbonded] = auto_show_nonbonded; /* show lines by default */

		  }
        if(ok) {
          p=nskip(p,2);
          p=ncopy(cc,p,3);
          if(sscanf(cc,"%d",&atInfo[a].formalCharge)==1) {
            if(atInfo[a].formalCharge) {
              atInfo[a].formalCharge = 4-atInfo[a].formalCharge;
            }
          }
          p=ncopy(cc,p,3);
          if(sscanf(cc,"%d",&atInfo[a].stereo)!=1) 
            atInfo[a].stereo=0;
        }
		  if(ok&&atInfo) {
          atInfo[a].id = a+1;
			 strcpy(atInfo[a].resn,resn);
			 atInfo[a].hetatm=true;
			 AtomInfoAssignParameters(atInfo+a);
			 atInfo[a].color=AtomInfoGetColor(atInfo+a);
          atInfo[a].alt[0]=0;
          atInfo[a].segi[0]=0;
          atInfo[a].resi[0]=0;
		  }
		  p=nextline(p);
		  if(!ok)
			 break;
		}
  }
  if(ok) {
	 bond=VLAlloc(BondType,nBond);
	 ii=bond;
	 for(a=0;a<nBond;a++)
		{
		  if(ok) {
			 p=ncopy(cc,p,3);
			 if(sscanf(cc,"%d",&ii->index[0])!=1)
				ok=ErrMessage("ReadMOLFile","bad bond atom");
		  }
		  
		  if(ok) {  
			 p=ncopy(cc,p,3);
			 if(sscanf(cc,"%d",&ii->index[1])!=1)
				ok=ErrMessage("ReadMOLFile","bad bond atom");
		  }

		  if(ok) {  
			 p=ncopy(cc,p,3);
			 if(sscanf(cc,"%d",&ii->order)!=1)
				ok=ErrMessage("ReadMOLFile","bad bond order");
		  }
        if(ok) {
			 p=ncopy(cc,p,3);
			 if(sscanf(cc,"%d",&ii->stereo)!=1)
            ii->stereo=0;

        }
        ii++;
		  if(!ok)
			 break;
		  p=nextline(p);
		}
	 ii=bond;
	 for(a=0;a<nBond;a++) {
      ii->index[0]--;/* adjust bond indexs down one */
      ii->index[1]--;
      ii++;
	 }
  }
  while(*p) { /* read M  CHG records */
    p=ncopy(cc,p,6);
    if(!strcmp(cc,"M  CHG")) {
      p=ncopy(cc,p,3);
      if(sscanf(cc,"%d",&cnt)==1) {
        while(cnt--) {
          p=ncopy(cc,p,4);
          p=ncopy(cc1,p,4);
          if(!((*cc)||(*cc1))) break;
          if((sscanf(cc,"%d",&atm)==1)&&
             (sscanf(cc1,"%d",&chg)==1)) {
            atm--;
            if((atm>=0)&&(atm<nAtom))
              atInfo[atm].formalCharge = chg;
          }
        }
      }
    }
    p=nextline(p);
  }
  if(ok) {
	 cset = CoordSetNew();
	 cset->NIndex=nAtom;
	 cset->Coord=coord;
	 cset->NTmpBond=nBond;
	 cset->TmpBond=bond;
    strcpy(cset->Name,nameTmp);
  } else {
	 VLAFreeP(bond);
	 VLAFreeP(coord);
  }
  if(atInfoPtr)
	 *atInfoPtr = atInfo;
  return(cset);
}

/*========================================================================*/
ObjectMolecule *ObjectMoleculeReadMOLStr(ObjectMolecule *I,char *MOLStr,int frame,int discrete)
{
  int ok = true;
  CoordSet *cset=NULL;
  AtomInfoType *atInfo;
  int isNew;
  int nAtom;

  if(!I) 
	 isNew=true;
  else 
	 isNew=false;

  if(isNew) {
    I=(ObjectMolecule*)ObjectMoleculeNew(discrete);
    atInfo = I->AtomInfo;
    isNew = true;
  } else {
    atInfo=VLAMalloc(10,sizeof(AtomInfoType),2,true); /* autozero here is important */
    isNew = false;
  }

  if(isNew) {
    AtomInfoPrimeColors();
    I->Obj.Color = AtomInfoGetCarbColor();
  }

  cset=ObjectMoleculeMOLStr2CoordSet(MOLStr,&atInfo);
  
  if(!cset) 
	 {
      ObjectMoleculeFree(I);
      I=NULL;
		ok=false;
	 }
  
  if(ok)
	 {
		if(frame<0)
		  frame=I->NCSet;
		if(I->NCSet<=frame)
		  I->NCSet=frame+1;
		VLACheck(I->CSet,CoordSet*,frame);
      
      nAtom=cset->NIndex;
      
      cset->Obj = I;
      cset->fEnumIndices(cset);
      if(cset->fInvalidateRep)
        cset->fInvalidateRep(cset,cRepAll,cRepInvRep);
      if(isNew) {		
        I->AtomInfo=atInfo; /* IMPORTANT to reassign: this VLA may have moved! */
      } else {
        ObjectMoleculeMerge(I,atInfo,cset,false,cAIC_MOLMask); /* NOTE: will release atInfo */
      }

      if(isNew) I->NAtom=nAtom;
      if(frame<0) frame=I->NCSet;
      VLACheck(I->CSet,CoordSet*,frame);
      if(I->NCSet<=frame) I->NCSet=frame+1;
      if(I->CSet[frame]) I->CSet[frame]->fFree(I->CSet[frame]);
      I->CSet[frame] = cset;
      
      if(isNew) I->NBond = ObjectMoleculeConnect(I,&I->Bond,I->AtomInfo,cset,false);
      
      SceneCountFrames();
      ObjectMoleculeExtendIndices(I);
      ObjectMoleculeSort(I);
      ObjectMoleculeUpdateIDNumbers(I);
      ObjectMoleculeUpdateNonbonded(I);
	 }
  return(I);
}
/*========================================================================*/
ObjectMolecule *ObjectMoleculeLoadMOLFile(ObjectMolecule *obj,char *fname,int frame,int discrete)
{
  ObjectMolecule* I=NULL;
  int ok=true;
  FILE *f;
  long size;
  char *buffer,*p;

  f=fopen(fname,"rb");
  if(!f)
	 ok=ErrMessage("ObjectMoleculeLoadMOLFile","Unable to open file!");
  else
	 {
      PRINTFB(FB_ObjectMolecule,FB_Blather)
        " ObjectMoleculeLoadMOLFile: Loading from %s.\n",fname
        ENDFB;
		
		fseek(f,0,SEEK_END);
      size=ftell(f);
		fseek(f,0,SEEK_SET);

		buffer=(char*)mmalloc(size+255);
		ErrChkPtr(buffer);
		p=buffer;
		fseek(f,0,SEEK_SET);
		fread(p,size,1,f);
		p[size]=0;
		fclose(f);
		I=ObjectMoleculeReadMOLStr(obj,buffer,frame,discrete);
		mfree(buffer);
	 }

  return(I);
}

/*========================================================================*/
void ObjectMoleculeMerge(ObjectMolecule *I,AtomInfoType *ai,
                         CoordSet *cs,int bondSearchFlag,int aic_mask)
{
  int *index,*outdex,*a2i,*i2a;
  BondType *bond=NULL;
  int a,b,c,lb=0,nb,ac,a1,a2;
  int found;
  int nAt,nBd,nBond;
  int expansionFlag = false;
  AtomInfoType *ai2;
  int oldNAtom,oldNBond;

  oldNAtom = I->NAtom;
  oldNBond = I->NBond;


  /* first, sort the coodinate set */
  
  index=AtomInfoGetSortedIndex(ai,cs->NIndex,&outdex);
  for(b=0;b<cs->NIndex;b++)
	 cs->IdxToAtm[b]=outdex[cs->IdxToAtm[b]];
  for(b=0;b<cs->NIndex;b++)
	 cs->AtmToIdx[b]=-1;
  for(b=0;b<cs->NIndex;b++)
	 cs->AtmToIdx[cs->IdxToAtm[b]]=b;
  ai2=(AtomInfoType*)VLAMalloc(cs->NIndex,sizeof(AtomInfoType),5,true); /* autozero here is important */
  for(a=0;a<cs->NIndex;a++) 
	 ai2[a]=ai[index[a]]; /* creates a sorted list of atom info records */
  VLAFreeP(ai);
  ai=ai2;

  /* now, match it up with the current object's atomic information */
	 
  for(a=0;a<cs->NIndex;a++) {
	 index[a]=-1;
	 outdex[a]=-1;
  }

  c=0;
  b=0;  
  for(a=0;a<cs->NIndex;a++) {
	 found=false;
    if(!I->DiscreteFlag) { /* don't even try matching for discrete objects */
      lb=b;
      while(b<I->NAtom) {
        ac=(AtomInfoCompare(ai+a,I->AtomInfo+b));
        if(!ac) {
          found=true;
          break;
        }
        else if(ac<0) {
          break;
        }
        b++;
      }
    }
	 if(found) {
		index[a]=b; /* store real atom index b for a in index[a] */
		b++;
	 } else {
	   index[a]=I->NAtom+c; /* otherwise, this is a new atom */
	   c++;
	   b=lb;
	 }
  }

  /* first, reassign atom info for matched atoms */

  /* allocate additional space */
  if(c)
	{
	  expansionFlag=true;
	  nAt=I->NAtom+c;
	} else {
     nAt=I->NAtom;
   }
  
  if(expansionFlag) {
	VLACheck(I->AtomInfo,AtomInfoType,nAt);
  }

  /* allocate our new x-ref tables */
  if(nAt<I->NAtom) nAt=I->NAtom;
  a2i = Alloc(int,nAt);
  i2a = Alloc(int,cs->NIndex);
  if(nAt) {
    ErrChkPtr(a2i);
  }
  if(cs->NIndex){
    ErrChkPtr(i2a);
  }
  
  for(a=0;a<cs->NIndex;a++) /* a is in original file space */
    {
		a1=cs->IdxToAtm[a]; /* a1 is in sorted atom info space */
		a2=index[a1];
		i2a[a]=a2; /* a2 is in object space */
      if(a2<oldNAtom)
        AtomInfoCombine(I->AtomInfo+a2,ai+a1,aic_mask);
      else
        *(I->AtomInfo+a2)=*(ai+a1);
    }
  
  if(I->DiscreteFlag) {
    if(I->NDiscrete<nAt) {
      VLACheck(I->DiscreteAtmToIdx,int,nAt);
      VLACheck(I->DiscreteCSet,CoordSet*,nAt);    
      for(a=I->NDiscrete;a<nAt;a++) {
        I->DiscreteAtmToIdx[a]=-1;
        I->DiscreteCSet[a]=NULL;
      }
    }
    I->NDiscrete = nAt;
  }
  
  cs->NAtIndex = nAt;
  I->NAtom = nAt;
  
  FreeP(cs->AtmToIdx);
  FreeP(cs->IdxToAtm);
  cs->AtmToIdx = a2i;
  cs->IdxToAtm = i2a;

  if(I->DiscreteFlag) {
    FreeP(cs->AtmToIdx);
    for(a=0;a<cs->NIndex;a++) {
      I->DiscreteAtmToIdx[cs->IdxToAtm[a]]=a;
      I->DiscreteCSet[cs->IdxToAtm[a]] = cs;
    }
  } else {
    for(a=0;a<cs->NAtIndex;a++)
      cs->AtmToIdx[a]=-1;
    for(a=0;a<cs->NIndex;a++)
      cs->AtmToIdx[cs->IdxToAtm[a]]=a;
  }
  
  VLAFreeP(ai);
  AtomInfoFreeSortedIndexes(index,outdex);
  
  /* now find and integrate and any new bonds */
  if(expansionFlag) { /* expansion flag means we have introduced at least 1 new atom */
    nBond = ObjectMoleculeConnect(I,&bond,I->AtomInfo,cs,bondSearchFlag);
    if(nBond) {
      index=Alloc(int,nBond);
      
      c=0;
      b=0;  
      nb=0;
      for(a=0;a<nBond;a++) { /* iterate over new bonds */
        found=false;
        b=nb; /* pick up where we left off */
        while(b<I->NBond) { 
          ac=BondCompare(bond+a,I->Bond+b);
          if(!ac) { /* zero is a match */
            found=true;
            break;
          } else if(ac<0) { /* gone past position of this bond */
            break;
          }
          b++; /* no match yet, keep looking */
        }
        if(found) {
          index[a]=b; /* existing bond...*/
          nb=b+1;
        } else { /* this is a new bond, save index and increment */
          index[a]=I->NBond+c;
          c++; 
        }
      }
      /* first, reassign atom info for matched atoms */
      if(c) {
        /* allocate additional space */
        nBd=I->NBond+c;
        
        VLACheck(I->Bond,BondType,nBd);
        
        for(a=0;a<nBond;a++) /* copy the new bonds */
          {
            a2=index[a];
            if(a2 >= I->NBond) { 
              I->Bond[a2] = bond[a];
            }
          }
        I->NBond=nBd;
      }
      FreeP(index);
    }
    VLAFreeP(bond);
  }

  if(oldNAtom) {
    if(oldNAtom==I->NAtom) {
      if(oldNBond!=I->NBond) {
        ObjectMoleculeInvalidate(I,cRepAll,cRepInvBonds);
      }
    } else {
      ObjectMoleculeInvalidate(I,cRepAll,cRepInvAtoms);
    }
  }

}
/*========================================================================*/
#if 0
ObjectMolecule *ObjectMoleculeLoadPDBFile(ObjectMolecule *obj,char *fname,
                                          int frame,int discrete,M4XAnnoType *m4x)
{
  ObjectMolecule *I=NULL;
  int ok=true;
  FILE *f;
  long size;
  char *buffer,*p;

  f=fopen(fname,"rb");
  if(!f) {
    PRINTFB(FB_ObjectMolecule,FB_Errors)
      "ObjectMolecule-ERROR: Unable to open file '%s'\n",fname
      ENDFB;
    ok=false;
  } else
	 {
      PRINTFB(FB_ObjectMolecule,FB_Blather)
        " ObjectMoleculeLoadPDBFile: Loading from %s.\n",fname
        ENDFB;
		
		fseek(f,0,SEEK_END);
      size=ftell(f);
		fseek(f,0,SEEK_SET);

		buffer=(char*)mmalloc(size+255);
		ErrChkPtr(buffer);
		p=buffer;
		fseek(f,0,SEEK_SET);
		fread(p,size,1,f);
		p[size]=0;
		fclose(f);

		I=ObjectMoleculeReadPDBStr(obj,buffer,frame,discrete,m4x);

		mfree(buffer);
	 }

  return(I);
}
#endif
/*========================================================================*/
void ObjectMoleculeAppendAtoms(ObjectMolecule *I,AtomInfoType *atInfo,CoordSet *cs)
{
  int a;
  BondType *ii;
  BondType *si;
  AtomInfoType *src,*dest;
  int nAtom,nBond;

  if(I->NAtom) {
	 nAtom = I->NAtom+cs->NIndex;
	 VLACheck(I->AtomInfo,AtomInfoType,nAtom);	 
	 dest = I->AtomInfo+I->NAtom;
	 src = atInfo;
	 for(a=0;a<cs->NIndex;a++)
		*(dest++)=*(src++);
	 I->NAtom=nAtom;
	 VLAFreeP(atInfo);
  } else {
	 if(I->AtomInfo)
		VLAFreeP(I->AtomInfo);
	 I->AtomInfo = atInfo;
	 I->NAtom=cs->NIndex;
  }
  nBond=I->NBond+cs->NTmpBond;
  if(!I->Bond)
	 I->Bond=VLAlloc(BondType,nBond);
  VLACheck(I->Bond,BondType,nBond);
  ii=I->Bond+I->NBond;
  si=cs->TmpBond;
  for(a=0;a<cs->NTmpBond;a++)
	 {
		ii->index[0]=cs->IdxToAtm[si->index[0]];
		ii->index[1]=cs->IdxToAtm[si->index[1]];
      ii->order=si->order;
      ii->stereo=si->stereo;
      ii->id=-1;
      ii++;
      si++;
	 }
  I->NBond=nBond;
}
/*========================================================================*/
CoordSet *ObjectMoleculeGetCoordSet(ObjectMolecule *I,int setIndex)
{
  if((setIndex>=0)&&(setIndex<I->NCSet))
	 return(I->CSet[setIndex]);
  else
	 return(NULL);
}
/*========================================================================*/
void ObjectMoleculeTransformTTTf(ObjectMolecule *I,float *ttt,int frame) 
{
  int b;
  CoordSet *cs;
  for(b=0;b<I->NCSet;b++)
	{
     if((frame<0)||(frame==b)) {
       cs=I->CSet[b];
       if(cs) {
         if(cs->fInvalidateRep)
           cs->fInvalidateRep(I->CSet[b],cRepAll,cRepInvCoord);
         MatrixApplyTTTfn3f(cs->NIndex,cs->Coord,ttt,cs->Coord);
       }
     }
	}
}
/*========================================================================*/
void ObjectMoleculeSeleOp(ObjectMolecule *I,int sele,ObjectMoleculeOpRec *op)
{
  int a,b,c,s,d,t_i;
  int a1,ind;
  float r,rms;
  float v1[3],v2,*vv1,*vv2,*coord,*vt,*vt1,*vt2;
  int inv_flag;
  int hit_flag = false;
  int ok = true;
  int cnt;
  int skip_flag;
  int match_flag=false;
  int offset;
  int priority;
  CoordSet *cs;
  AtomInfoType *ai,*ai0,*ai_option;
  
  PRINTFD(FB_ObjectMolecule)
    " ObjectMoleculeSeleOp-DEBUG: sele %d op->code %d\n",sele,op->code
    ENDFD;

  if(sele>=0) {
	SelectorUpdateTable();
   /* always run on entry */
	switch(op->code) {
	case OMOP_ALTR: 
   case OMOP_AlterState:
     PBlock();
     /* PBlockAndUnlockAPI() is not safe.
      * what if "v" is invalidated by another thread? */
     break;
   }
   /* */
	switch(op->code) {
	case OMOP_AddHydrogens:
     ObjectMoleculeAddSeleHydrogens(I,sele);
     break;
#ifdef _OLD_CODE
     if(!ObjectMoleculeVerifyChemistry(I)) {
       ErrMessage(" AddHydrogens","missing chemical geometry information.");
     } else {
       doneFlag=false;
       while(!doneFlag) {
         doneFlag=true;
         a=0;
         while(a<I->NAtom) {
           ai=I->AtomInfo + a;
           s=I->AtomInfo[a].selEntry;
           if(SelectorIsMember(s,sele))
             if(ObjectMoleculeFillOpenValences(I,a)) {
               hit_flag=true;
               doneFlag=false;
             }
           a++; /* realize that the atom list may have been resorted */
         }
       }
       if(hit_flag) {
         ObjectMoleculeSort(I);
         ObjectMoleculeUpdateIDNumbers(I);
       } 
     }
     break;
#endif

	case OMOP_PrepareFromTemplate:
     ai0=op->ai; /* template atom */
     for(a=0;a<I->NAtom;a++)
       {
         s=I->AtomInfo[a].selEntry;
         if(SelectorIsMember(s,sele))
           {
             ai = I->AtomInfo + a;
             ai->hetatm=ai0->hetatm;
             ai->flags=ai0->flags;
             strcpy(ai->chain,ai0->chain);
             strcpy(ai->alt,ai0->alt);
             strcpy(ai->segi,ai0->segi);
             if(op->i1==1) { /* mode 1, merge residue information */
               strcpy(ai->resi,ai0->resi);
               ai->resv=ai0->resv;
               strcpy(ai->resn,ai0->resn);    
             }
             if((ai->elem[0]==ai0->elem[0])&&(ai->elem[1]==ai0->elem[1]))
               ai->color=ai0->color;
             else
               ai->color=AtomInfoGetColor(ai);
             for(b=0;b<cRepCnt;b++)
               ai->visRep[b]=ai0->visRep[b];
             ai->id=-1;
             op->i2++;
           }
       }
     break;
     
	case OMOP_PDB1:
	  for(b=0;b<I->NCSet;b++)
       if(I->CSet[b])
		  {
			if((b==op->i1)||(op->i1<0))
			  for(a=0;a<I->NAtom;a++)
				{
				  s=I->AtomInfo[a].selEntry;
              if(SelectorIsMember(s,sele))
                {
                  if(I->DiscreteFlag) {
                    if(I->CSet[b]==I->DiscreteCSet[a])
                      ind=I->DiscreteAtmToIdx[a];
                    else
                      ind=-1;
                  } else 
                    ind=I->CSet[b]->AtmToIdx[a];
                  if(ind>=0) 
                    CoordSetAtomToPDBStrVLA(&op->charVLA,&op->i2,I->AtomInfo+a,
                                            I->CSet[b]->Coord+(3*ind),op->i3);
                  op->i3++;
                }
				}
		  }
	  break;
	case OMOP_AVRT: /* average vertex coordinate */
     for(a=0;a<I->NAtom;a++)
       {
         s=I->AtomInfo[a].selEntry;
         if((priority=SelectorIsMember(s,sele)))
           {
             cnt=0;
             for(b=0;b<I->NCSet;b++) {
               if(I->CSet[b])
                 {
                   if(I->DiscreteFlag) {
                     if(I->CSet[b]==I->DiscreteCSet[a])
                       a1=I->DiscreteAtmToIdx[a];
                     else
                       a1=-1;
                   } else 
                     a1=I->CSet[b]->AtmToIdx[a];
                   if(a1>=0) {
                     if(!cnt) {
                       VLACheck(op->vv1,float,(op->nvv1*3)+2);
                       VLACheck(op->vc1,int,op->nvv1);
                     }
                     cnt++;
                     vv2=I->CSet[b]->Coord+(3*a1);
                     vv1=op->vv1+(op->nvv1*3);
                     *(vv1++)+=*(vv2++);
                     *(vv1++)+=*(vv2++);
                     *(vv1++)+=*(vv2++);
                   }
                 }
             }
             op->vc1[op->nvv1]=cnt;
             if(cnt) {
               if(op->vp1) {
                 VLACheck(op->vp1,int,op->nvv1);
                 op->vp1[op->nvv1] = priority;
               }
               op->nvv1++;
             }
           }
       }
     break;
	case OMOP_StateVRT: /* state vertex coordinate */
     for(a=0;a<I->NAtom;a++)
       {
         s=I->AtomInfo[a].selEntry;
         if((priority=SelectorIsMember(s,sele)))
           {
             cnt=0;
             b=op->i1;
             if(b<I->NCSet)
               if(I->CSet[b])
                 {
                   if(I->DiscreteFlag) {
                     if(I->CSet[b]==I->DiscreteCSet[a])
                       a1=I->DiscreteAtmToIdx[a];
                     else
                       a1=-1;
                   } else 
                     a1=I->CSet[b]->AtmToIdx[a];
                   if(a1>=0) {
                     if(!cnt) {
                       VLACheck(op->vv1,float,(op->nvv1*3)+2);
                       VLACheck(op->vc1,int,op->nvv1);
                     }
                     cnt++;
                     vv2=I->CSet[b]->Coord+(3*a1);
                     vv1=op->vv1+(op->nvv1*3);
                     *(vv1++)+=*(vv2++);
                     *(vv1++)+=*(vv2++);
                     *(vv1++)+=*(vv2++);
                   }
                 }
             op->vc1[op->nvv1]=cnt;
             if(cnt) {
               if(op->vp1) {
                 VLACheck(op->vp1,int,op->nvv1);
                 op->vp1[op->nvv1] = priority;
               }
               op->nvv1++;
             }
           }
       }
     break;
	case OMOP_SFIT: /* state fitting within a single object */
     vt = Alloc(float,3*op->nvv2);
     cnt = 0;
     for(a=0;a<I->NAtom;a++)
       {
         s=I->AtomInfo[a].selEntry;
         if(SelectorIsMember(s,sele))
           {
             cnt++;
             break;
           }
       }
     if(cnt) { /* only perform action for selected object */
       
       for(b=0;b<I->NCSet;b++) {
         rms = -1.0;
         vt1 = vt; /* reset target vertex pointers */
         vt2 = op->vv2;
         t_i = 0; /* original target vertex index */
         if(I->CSet[b]&&(b!=op->i2))
           {
             op->nvv1=0;
             for(a=0;a<I->NAtom;a++)
               {
                 s=I->AtomInfo[a].selEntry;
                 if(SelectorIsMember(s,sele))
                   {
                     if(I->DiscreteFlag) {
                       if(I->CSet[b]==I->DiscreteCSet[a])
                         a1=I->DiscreteAtmToIdx[a];
                       else
                         a1=-1;
                     } else 
                       a1=I->CSet[b]->AtmToIdx[a];
                     if(a1>=0) {

                       match_flag=false;
                       while(t_i<op->nvv2) {
                         if(op->i1VLA[t_i]==a) {/* same atom? */
                           match_flag=true;
                           break;
                         }
                         if(op->i1VLA[t_i]<a) { /* catch up? */
                           t_i++;
                           vt2+=3;
                         } else 
                           break;
                       }
                       if(match_flag) {
                         VLACheck(op->vv1,float,(op->nvv1*3)+2);
                         vv2=I->CSet[b]->Coord+(3*a1);
                         vv1=op->vv1+(op->nvv1*3);
                         *(vv1++)=*(vv2++);
                         *(vv1++)=*(vv2++);
                         *(vv1++)=*(vv2++);
                         *(vt1++)=*(vt2);
                         *(vt1++)=*(vt2+1);
                         *(vt1++)=*(vt2+2);
                         op->nvv1++;
                       }
                     }
                   }
               }
             if(op->nvv1!=op->nvv2) {
               PRINTFB(FB_Executive,FB_Warnings)
                 "Executive-Warning: Missing atoms in state %d (%d instead of %d).\n",
                 b+1,op->nvv1,op->nvv2
                 ENDFB;
             }
             if(op->nvv1) {
               if(op->i1!=0) /* fitting flag */
                 rms = MatrixFitRMS(op->nvv1,op->vv1,vt,NULL,op->ttt);
               else 
                 rms = MatrixGetRMS(op->nvv1,op->vv1,vt,NULL);
               if(op->i1==2) 
                 ObjectMoleculeTransformTTTf(I,op->ttt,b);
             } else {
               PRINTFB(FB_Executive,FB_Warnings)
                 "Executive-Warning: No matches found for state %d.\n",b+1
                 ENDFB;
             }
           }
         VLACheck(op->f1VLA,float,b);
         op->f1VLA[b]=rms;
       }
       VLASize(op->f1VLA,float,I->NCSet);  /* NOTE this action is object-specific! */
     }
     FreeP(vt);
     break;
	case OMOP_SetGeometry: 
     for(a=0;a<I->NAtom;a++)
       {
         s=I->AtomInfo[a].selEntry;
         if(SelectorIsMember(s,sele))
           {
             ai = I->AtomInfo + a;
             ai->geom=op->i1;
             ai->valence=op->i2;
             op->i3++;
             hit_flag=true;
             break;
           }
       }
     break;
	case OMOP_OnOff:
     for(a=0;a<I->NAtom;a++)
       {
         s=I->AtomInfo[a].selEntry;
         if(SelectorIsMember(s,sele))
           {
             hit_flag=true;
             break;
           }
       }
     break;
	case OMOP_SaveUndo: /* save undo */
     for(a=0;a<I->NAtom;a++)
       {
         s=I->AtomInfo[a].selEntry;
         if(SelectorIsMember(s,sele))
           {
             hit_flag=true;
             break;
           }
       }
     break;
	case OMOP_Identify: /* identify atoms */
     for(a=0;a<I->NAtom;a++)
       {
         s=I->AtomInfo[a].selEntry;
         if(SelectorIsMember(s,sele))
           {
             VLACheck(op->i1VLA,int,op->i1);
             op->i1VLA[op->i1++]=I->AtomInfo[a].id;
           }
       }
     break;
	case OMOP_GetBFactors: 
     ai=I->AtomInfo;
     for(a=0;a<I->NAtom;a++)
       {
         s=ai->selEntry;
         if(SelectorIsMember(s,sele))
           {
             op->ff1[op->i1]=ai->b;
             op->i1++;
           }
         ai++;
       }
     break;
	case OMOP_GetOccupancies: 
     ai=I->AtomInfo;
     for(a=0;a<I->NAtom;a++)
       {
         s=ai->selEntry;
         if(SelectorIsMember(s,sele))
           {
             op->ff1[op->i1]=ai->q;
             op->i1++;
           }
         ai++;
       }
     break;
	case OMOP_GetPartialCharges: 
     ai=I->AtomInfo;
     for(a=0;a<I->NAtom;a++)
       {
         s=ai->selEntry;
         if(SelectorIsMember(s,sele))
           {
             op->ff1[op->i1]=ai->partialCharge;
             op->i1++;
           }
         ai++;
       }
     break;
	case OMOP_IdentifyObjects: /* identify atoms */
     for(a=0;a<I->NAtom;a++)
       {
         s=I->AtomInfo[a].selEntry;
         if(SelectorIsMember(s,sele))
           {
             VLACheck(op->i1VLA,int,op->i1);
             op->i1VLA[op->i1]=I->AtomInfo[a].id; 
             VLACheck(op->obj1VLA,ObjectMolecule*,op->i1);
             op->obj1VLA[op->i1]=I;
             op->i1++;
           }
       }
     break;
	case OMOP_Index: /* identify atoms */
     for(a=0;a<I->NAtom;a++)
       {
         s=I->AtomInfo[a].selEntry;
         if(SelectorIsMember(s,sele))
           {
             VLACheck(op->i1VLA,int,op->i1);
             op->i1VLA[op->i1]=a; /* NOTE: need to incr by 1 before python */
             VLACheck(op->obj1VLA,ObjectMolecule*,op->i1);
             op->obj1VLA[op->i1]=I;
             op->i1++;
           }
       }
     break;
	case OMOP_GetObjects: /* identify atoms */
     for(a=0;a<I->NAtom;a++)
       {
         s=I->AtomInfo[a].selEntry;
         if(SelectorIsMember(s,sele))
           {
             VLACheck(op->obj1VLA,ObjectMolecule*,op->i1);
             op->obj1VLA[op->i1]=I;
             op->i1++;
             break;
           }
       }
     break;
	case OMOP_CountAtoms: /* count atoms in object, in selection */
     ai=I->AtomInfo;
     for(a=0;a<I->NAtom;a++)
       {
         s=ai->selEntry;
         if(SelectorIsMember(s,sele))
           op->i1++;
         ai++;
       }
     break;
   case OMOP_PhiPsi:
     ai=I->AtomInfo;
     for(a=0;a<I->NAtom;a++)
       {
         s=ai->selEntry;
         if(SelectorIsMember(s,sele)) {
           VLACheck(op->i1VLA,int,op->i1);
           op->i1VLA[op->i1]=a;
           VLACheck(op->obj1VLA,ObjectMolecule*,op->i1);
           op->obj1VLA[op->i1]=I;
           VLACheck(op->f1VLA,float,op->i1);
           VLACheck(op->f2VLA,float,op->i1);
           if(ObjectMoleculeGetPhiPsi(I,a,op->f1VLA+op->i1,op->f2VLA+op->i1,op->i2))
             op->i1++;
         }
         ai++; 
       }
     break;
	case OMOP_Cartoon: /* adjust cartoon type */
     ai=I->AtomInfo;
     for(a=0;a<I->NAtom;a++)
       {
         s=ai->selEntry;
         if(SelectorIsMember(s,sele)) {
           ai->cartoon = op->i1;
           op->i2++;
         }
         ai++; 
       }
     break;
	case OMOP_Protect: /* protect atoms from movement */
     ai=I->AtomInfo;
     for(a=0;a<I->NAtom;a++)
       {
         s=ai->selEntry;
         if(SelectorIsMember(s,sele))
           {
             ai->protekted = op->i1;
             op->i2++;
           }
         ai++;
       }
     break;
	case OMOP_Mask: /* protect atoms from selection */
     ai=I->AtomInfo;
     for(a=0;a<I->NAtom;a++)
       {
         s=ai->selEntry;
         if(SelectorIsMember(s,sele))
           {
             ai->masked = op->i1;
             op->i2++;
           }
         ai++;
       }
     break;
	case OMOP_SetB: /* set B-value */
     ai=I->AtomInfo;
     for(a=0;a<I->NAtom;a++)
       {
         s=ai->selEntry;
         if(SelectorIsMember(s,sele))
           {
             ai->b = op->f1;
             op->i2++;
           }
         ai++;
       }
     break;
	case OMOP_Remove: /* flag atoms for deletion */
     ai=I->AtomInfo;
     if(I->DiscreteFlag) /* for now, can't remove atoms from discrete objects */
       ErrMessage("Remove","Can't remove atoms from discrete objects.");
     else
       for(a=0;a<I->NAtom;a++)
         {         
           ai->deleteFlag=false;
           s=ai->selEntry;
           if(SelectorIsMember(s,sele))
             {
               ai->deleteFlag=true;
               op->i1++;
             }
           ai++;
         }
     break;

   case OMOP_GetChains:
     ai=I->AtomInfo;
     for(a=0;a<I->NAtom;a++)
       {         
         s=ai->selEntry;
         if(SelectorIsMember(s,sele))
           {
             op->ii1[(int)ai->chain[0]]++;
             op->i1++;
           }
         ai++;
       }
     break;

   case OMOP_Spectrum:
     ai=I->AtomInfo;
     ai0=NULL;
     for(a=0;a<I->NAtom;a++)
       {         
         s=ai->selEntry;
         if(SelectorIsMember(s,sele))
           {
             skip_flag=false;
             if(op->i4&&ai0) /* byres and we've done a residue */
               if(AtomInfoSameResidue(ai,ai0))
                 skip_flag=true;
             if(!skip_flag) {
               c = (int)(0.49999+op->i1*(op->ff1[op->i3] - op->f1)/op->f2);
               if(c<0) c=0;
               if(c>=op->i1) c=op->i1-1;
               ai->color=op->ii1[c];

               /*               printf("%8.3 %8.3\n",ai->partial_charge,*/
               if(op->i4) { /* byres */
                 offset = -1;
                 while((a+offset)>=0) {
                   ai0 = I->AtomInfo + a + offset;
                   if(AtomInfoSameResidue(ai,ai0)) {
                     ai0->color = op->ii1[c];
                   } else 
                     break;
                   offset--;
                 }
                 offset = 1;
                 while((a+offset)<I->NAtom) {
                   ai0 = I->AtomInfo + a + offset;
                   if(AtomInfoSameResidue(ai,ai0)) {
                     ai0->color = op->ii1[c];
                   } else 
                     break;
                   offset++;
                 }
               }
               ai0=ai;
                      }
             op->i3++;

           }
         ai++;
       }
     break;

   case OMOP_SingleStateVertices: /* same as OMOP_VERT for a single state */
     ai=I->AtomInfo;
     if(op->cs1<I->NCSet) {
       if(I->CSet[op->cs1]) {
         b = op->cs1;
         for(a=0;a<I->NAtom;a++)
           {         
             s=ai->selEntry;
             if(SelectorIsMember(s,sele))
               {
                 op->i1++;

                 if(I->DiscreteFlag) {
                   if(I->CSet[b]==I->DiscreteCSet[a])
                     a1=I->DiscreteAtmToIdx[a];
                   else
                     a1=-1;
                 } else 
                   a1=I->CSet[b]->AtmToIdx[a];
                 if(a1>=0) {
                   VLACheck(op->vv1,float,(op->nvv1*3)+2);
                   vv2=I->CSet[b]->Coord+(3*a1);
                   vv1=op->vv1+(op->nvv1*3);
                   *(vv1++)=*(vv2++);
                   *(vv1++)=*(vv2++);
                   *(vv1++)=*(vv2++);
                   op->nvv1++;
                 }
               }
             ai++;
           }
       }
     }
     break;
   case OMOP_CSetIdxGetAndFlag:
     ai=I->AtomInfo;
     for(a=0;a<I->NAtom;a++)
       {         
         s=ai->selEntry;
         if(SelectorIsMember(s,sele))
           {
             for(b=op->cs1;b<=op->cs2;b++) {
               offset = b-op->cs1;
               if(b<I->NCSet) {
                 if(I->CSet[b]) {
                   if(I->DiscreteFlag) {
                     if(I->CSet[b]==I->DiscreteCSet[a])
                       a1=I->DiscreteAtmToIdx[a];
                     else
                       a1=-1;
                   } else 
                     a1=I->CSet[b]->AtmToIdx[a];
                   if(a1>=0) {
                     op->ii1[op->i1*offset+op->i2] = 1; /* presence flag */
                     vv1=op->vv1+3*(op->i1*offset+op->i2); /* atom-based offset */
                     vv2=I->CSet[b]->Coord+(3*a1);
                     *(vv1++)=*(vv2++);
                     *(vv1++)=*(vv2++);
                     *(vv1++)=*(vv2++);
                     op->nvv1++;
                   }
                 }
               }
             }
             op->i2++; /* atom index field for atoms within selection...*/
           }
         ai++;
       }
     break;
   case OMOP_CSetIdxSetFlagged:
     ai=I->AtomInfo;
     hit_flag=false;
     for(a=0;a<I->NAtom;a++)
       {         
         s=ai->selEntry;
         if(SelectorIsMember(s,sele))
           {
             for(b=op->cs1;b<=op->cs2;b++) {
               offset = b-op->cs1;
               if(b<I->NCSet) {
                 if(I->CSet[b]) {
                   if(I->DiscreteFlag) {
                     if(I->CSet[b]==I->DiscreteCSet[a])
                       a1=I->DiscreteAtmToIdx[a];
                     else
                       a1=-1;
                   } else 
                     a1=I->CSet[b]->AtmToIdx[a];
                   if(a1>=0) {
                     if(op->ii1[op->i1*offset+op->i2]) { /* copy flag */
                       vv1=op->vv1+3*(op->i1*offset+op->i2); /* atom-based offset */
                       vv2=I->CSet[b]->Coord+(3*a1);
                       *(vv2++)=*(vv1++);
                       *(vv2++)=*(vv1++);
                       *(vv2++)=*(vv1++);
                       op->nvv1++;
                       hit_flag=true;
                     }
                   }
                 }
               }
             }
             op->i2++; /* atom index field for atoms within selection...*/
           }
         ai++;
       }
     break;
   default:
     ai = I->AtomInfo;
     for(a=0;a<I->NAtom;a++)
		 {
		   switch(op->code) { 
         case OMOP_Flag: 
           ai->flags &= op->i2; /* clear flag using mask */
           op->i4++;
           /* no break here is intentional!  */
         case OMOP_FlagSet:
         case OMOP_FlagClear:
		   case OMOP_COLR: /* normal atom based loops */
		   case OMOP_VISI:
         case OMOP_CheckVis:
		   case OMOP_TTTF:
         case OMOP_ALTR:
         case OMOP_LABL:
         case OMOP_AlterState:
			 s=ai->selEntry;
          if(SelectorIsMember(s,sele))
            {
              switch(op->code) {
              case OMOP_Flag:
                ai->flags |= op->i1; /* set flag */
                op->i3++;
                break;
              case OMOP_FlagSet:
                ai->flags |= op->i1; /* set flag */
                op->i3++;
                break;
              case OMOP_FlagClear:
                ai->flags &= op->i2; /* clear flag */
                op->i3++;
                break;
              case OMOP_VISI:
                if(op->i1<0)
                  for(d=0;d<cRepCnt;d++) 
                    ai->visRep[d]=op->i2;                      
                else {
                  ai->visRep[op->i1]=op->i2;
                  if(op->i1==cRepCell) I->Obj.RepVis[cRepCell]=op->i2;
                }
                break;
                break;
              case OMOP_CheckVis:
                if(ai->visRep[op->i1]) {
                  op->i2 = true;
                }
                break;
              case OMOP_COLR:
                ai->color=op->i1;
                op->i2++;
                break;
              case OMOP_TTTF:
                hit_flag=true;
                break;
              case OMOP_LABL:
                if (ok) {
                  if(!op->s1[0]) {
                    ai->label[0]=0;
                    op->i1++;
                    ai->visRep[cRepLabel]=false;
                    hit_flag=true;
                  }  else {
                    if(PLabelAtom(&I->AtomInfo[a],op->s1,a)) {
                      op->i1++;
                      ai->visRep[cRepLabel]=true;
                      hit_flag=true;
                    } else
                      ok=false;
                  }
                }
                break;
              case OMOP_ALTR:
                if (ok) {
                  if(PAlterAtom(&I->AtomInfo[a],op->s1,op->i2,I->Obj.Name,a))
                    op->i1++;
                  else
                    ok=false;
                }
                break;
              case OMOP_AlterState:
                if (ok) {
                  if(op->i2<I->NCSet) {
                    cs = I->CSet[op->i2];
                    if(cs) {
                      if(I->DiscreteFlag) {
                        if(cs==I->DiscreteCSet[a])
                          a1=I->DiscreteAtmToIdx[a];
                        else
                          a1=-1;
                      } else 
                        a1=cs->AtmToIdx[a];
                      if(a1>=0) {
                        if(op->i4) 
                          ai_option = I->AtomInfo + a;
                        else
                          ai_option = NULL;
                        if(PAlterAtomState(cs->Coord+(a1*3),op->s1,op->i3,ai_option)) {
                          op->i1++;
                          hit_flag=true;
                        } else
                          ok=false;
                      }
                    }
                  }
                }
                break;
              }
              break;
            }
          break;

          /* coord-set based properties, iterating only a single coordinate set */
         case OMOP_CSetMinMax:          
         case OMOP_CSetCameraMinMax:          
         case OMOP_CSetMaxDistToPt:
         case OMOP_CSetSumVertices:
         case OMOP_CSetMoment: 
           if((op->cs1>=0)&&(op->cs1<I->NCSet)) {
             cs=I->CSet[op->cs1];
             if(cs) {
               s=ai->selEntry;
               if(SelectorIsMember(s,sele))
                 {
                   switch(op->code) {
                   case OMOP_CSetSumVertices:
                     if(I->DiscreteFlag) {
                       if(cs==I->DiscreteCSet[a])
                         a1=I->DiscreteAtmToIdx[a];
                       else
                         a1=-1;
                     } else 
                       a1=cs->AtmToIdx[a];
                     if(a1>=0)
                       {
                         coord = cs->Coord+3*a1;
                         if(op->i2) /* do we want object-transformed coordinates? */
                           if(I->Obj.TTTFlag) {
                             transformTTT44f3f(I->Obj.TTT,coord,v1);
                             coord=v1;
                           }
                         add3f(op->v1,coord,op->v1);
                         op->i1++;
                       }
                     break;
                   case OMOP_CSetMinMax:
                     if(I->DiscreteFlag) {
                       if(cs==I->DiscreteCSet[a])
                         a1=I->DiscreteAtmToIdx[a];
                       else
                         a1=-1;
                     } else 
                       a1=cs->AtmToIdx[a];
                     if(a1>=0)
                       {
                         coord = cs->Coord+3*a1;
                         if(op->i2) /* do we want object-transformed coordinates? */
                           if(I->Obj.TTTFlag) {
                             transformTTT44f3f(I->Obj.TTT,coord,v1);
                             coord=v1;
                           }
                         if(op->i1) {
                           for(c=0;c<3;c++) {
                             if(*(op->v1+c)>*(coord+c)) *(op->v1+c)=*(coord+c);
                             if(*(op->v2+c)<*(coord+c)) *(op->v2+c)=*(coord+c);
                           }
                         } else {
                           for(c=0;c<3;c++) {
                             *(op->v1+c)=*(coord+c);
                             *(op->v2+c)=*(coord+c);
                           }
                         }
                         op->i1++;
                       }
                     break;
                   case OMOP_CSetCameraMinMax:
                     if(I->DiscreteFlag) {
                       if(cs==I->DiscreteCSet[a])
                         a1=I->DiscreteAtmToIdx[a];
                       else
                         a1=-1;
                     } else 
                       a1=cs->AtmToIdx[a];
                     if(a1>=0)
                       {
                         coord = cs->Coord+3*a1;
                         if(op->i2) /* do we want object-transformed coordinates? */
                           if(I->Obj.TTTFlag) {
                             transformTTT44f3f(I->Obj.TTT,coord,v1);
                             coord=v1;
                           }
                         MatrixTransform3f(op->mat1,coord,v1); /* convert to view-space */
                         coord=v1;
                         if(op->i1) {
                           for(c=0;c<3;c++) {
                             if(*(op->v1+c)>*(coord+c)) *(op->v1+c)=*(coord+c);
                             if(*(op->v2+c)<*(coord+c)) *(op->v2+c)=*(coord+c);
                           }
                         } else {
                           for(c=0;c<3;c++) {
                             *(op->v1+c)=*(coord+c);
                             *(op->v2+c)=*(coord+c);
                           }
                         }
                         op->i1++;
                       }
                     break;
                   case OMOP_CSetMaxDistToPt:
                     if(I->DiscreteFlag) {
                       if(cs==I->DiscreteCSet[a])
                         a1=I->DiscreteAtmToIdx[a];
                       else
                         a1=-1;
                     } else 
                       a1=cs->AtmToIdx[a];
                     if(a1>=0)
                       {
                         float dist;
                         coord = cs->Coord+3*a1;
                         if(op->i2) /* do we want object-transformed coordinates? */
                           if(I->Obj.TTTFlag) {
                             transformTTT44f3f(I->Obj.TTT,coord,v1);
                             coord=v1;
                           }
                         dist = (float)diff3f(op->v1,coord);
                         if(dist>op->f1)
                           op->f1=dist;
                         op->i1++;
                       }
                     break;
                   case OMOP_CSetMoment: 
                     if(I->DiscreteFlag) {
                       if(cs==I->DiscreteCSet[a])
                         a1=I->DiscreteAtmToIdx[a];
                       else
                         a1=-1;
                     } else 
                       a1=cs->AtmToIdx[a];
                     if(a1>=0) {
                       subtract3f(cs->Coord+(3*a1),op->v1,v1);
                       v2=v1[0]*v1[0]+v1[1]*v1[1]+v1[2]*v1[2]; 
                       op->d[0][0] += v2 - v1[0] * v1[0];
                       op->d[0][1] +=    - v1[0] * v1[1];
                       op->d[0][2] +=    - v1[0] * v1[2];
                       op->d[1][0] +=    - v1[1] * v1[0];
                       op->d[1][1] += v2 - v1[1] * v1[1];
                       op->d[1][2] +=    - v1[1] * v1[2];
                       op->d[2][0] +=    - v1[2] * v1[0];
                       op->d[2][1] +=    - v1[2] * v1[1];
                       op->d[2][2] += v2 - v1[2] * v1[2];
                     }
                     break;
                     
                   }
                 }
             }
           }
           break;
		   default: /* coord-set based properties, iterating as all coordsets within atoms */
			 for(b=0;b<I->NCSet;b++)
			   if(I->CSet[b])
              {
                cs=I->CSet[b];
                inv_flag=false;
                s=ai->selEntry;
                if(SelectorIsMember(s,sele))
                  {
                    switch(op->code) {
                    case OMOP_SUMC:
                      if(I->DiscreteFlag) {
                        if(cs==I->DiscreteCSet[a])
                          a1=I->DiscreteAtmToIdx[a];
                        else
                          a1=-1;
                      } else 
                        a1=cs->AtmToIdx[a];
							 if(a1>=0)
							   {
                          coord = cs->Coord+3*a1;
                          if(op->i2) /* do we want object-transformed coordinates? */
                            if(I->Obj.TTTFlag) {
                              transformTTT44f3f(I->Obj.TTT,coord,v1);
                              coord=v1;
                            }
                          add3f(op->v1,coord,op->v1);
                          op->i1++;
							   }
							 break;
                    case OMOP_MNMX:
                      if(I->DiscreteFlag) {
                        if(cs==I->DiscreteCSet[a])
                          a1=I->DiscreteAtmToIdx[a];
                        else
                          a1=-1;
                      } else 
                        a1=cs->AtmToIdx[a];
							 if(a1>=0)
							   {
                          coord = cs->Coord+3*a1;
                          if(op->i2) /* do we want object-transformed coordinates? */
                            if(I->Obj.TTTFlag) {
                              transformTTT44f3f(I->Obj.TTT,coord,v1);
                              coord=v1;
                            }
                          if(op->i1) {
                            for(c=0;c<3;c++) {
                              if(*(op->v1+c)>*(coord+c)) *(op->v1+c)=*(coord+c);
                              if(*(op->v2+c)<*(coord+c)) *(op->v2+c)=*(coord+c);
                            }
                          } else {
                            for(c=0;c<3;c++) {
                              *(op->v1+c)=*(coord+c);
                              *(op->v2+c)=*(coord+c);
                            }
                          }
                          op->i1++;
							   }
							 break;
                    case OMOP_CameraMinMax:
                      if(I->DiscreteFlag) {
                        if(cs==I->DiscreteCSet[a])
                          a1=I->DiscreteAtmToIdx[a];
                        else
                          a1=-1;
                      } else 
                        a1=cs->AtmToIdx[a];
							 if(a1>=0)
							   {
                          coord = cs->Coord+3*a1;
                          if(op->i2) /* do we want object-transformed coordinates? */
                            if(I->Obj.TTTFlag) {
                              transformTTT44f3f(I->Obj.TTT,coord,v1);
                              coord=v1;
                            }
                          MatrixTransform3f(op->mat1,coord,v1); /* convert to view-space */
                          coord=v1;
                          if(op->i1) {
                            for(c=0;c<3;c++) {
                              if(*(op->v1+c)>*(coord+c)) *(op->v1+c)=*(coord+c);
                              if(*(op->v2+c)<*(coord+c)) *(op->v2+c)=*(coord+c);
                            }
                          } else {
                            for(c=0;c<3;c++) {
                              *(op->v1+c)=*(coord+c);
                              *(op->v2+c)=*(coord+c);
                            }
                          }
                          op->i1++;
							   }
							 break;
                    case OMOP_MaxDistToPt:
                      if(I->DiscreteFlag) {
                        if(cs==I->DiscreteCSet[a])
                          a1=I->DiscreteAtmToIdx[a];
                        else
                          a1=-1;
                      } else 
                        a1=cs->AtmToIdx[a];
							 if(a1>=0)
							   {
                          float dist;
                          coord = cs->Coord+3*a1;
                          if(op->i2) /* do we want object-transformed coordinates? */
                            if(I->Obj.TTTFlag) {
                              transformTTT44f3f(I->Obj.TTT,coord,v1);
                              coord=v1;
                            }
                          dist = (float)diff3f(coord,op->v1);
                          if(dist>op->f1)
                            op->f1=dist;
                          op->i1++;
							   }
							 break;
                    case OMOP_MDST: 
                      if(I->DiscreteFlag) {
                        if(cs==I->DiscreteCSet[a])
                          a1=I->DiscreteAtmToIdx[a];
                        else
                          a1=-1;
                      } else 
                        a1=cs->AtmToIdx[a];
							 if(a1>=0)
							   {
                          r=(float)diff3f(op->v1,cs->Coord+(3*a1));
                          if(r>op->f1)
                            op->f1=r;
							   }
							 break;
                    case OMOP_INVA:
                      if(I->DiscreteFlag) {
                        if(cs==I->DiscreteCSet[a])
                          a1=I->DiscreteAtmToIdx[a];
                        else
                          a1=-1;
                      } else 
                        a1=cs->AtmToIdx[a]; 
                      if(a1>=0)                     /* selection touches this coordinate set */ 
                        inv_flag=true;              /* so set the invalidation flag */
                      break;
                    case OMOP_VERT: 
                      if(I->DiscreteFlag) {
                        if(cs==I->DiscreteCSet[a])
                          a1=I->DiscreteAtmToIdx[a];
                        else
                          a1=-1;
                      } else 
                        a1=cs->AtmToIdx[a];
                      if(a1>=0) {
                        VLACheck(op->vv1,float,(op->nvv1*3)+2);
                        vv2=cs->Coord+(3*a1);
                        vv1=op->vv1+(op->nvv1*3);
                        *(vv1++)=*(vv2++);
                        *(vv1++)=*(vv2++);
                        *(vv1++)=*(vv2++);
                        op->nvv1++;
                      }
							 break;	
                    case OMOP_SVRT:  /* gives us only vertices for a specific coordinate set */
                      if(b==op->i1) {
                        if(I->DiscreteFlag) {
                          if(cs==I->DiscreteCSet[a])
                            a1=I->DiscreteAtmToIdx[a];
                          else
                            a1=-1;
                        } else 
                          a1=cs->AtmToIdx[a];
                        if(a1>=0) {
                          VLACheck(op->vv1,float,(op->nvv1*3)+2);
                          VLACheck(op->i1VLA,int,op->nvv1);
                          op->i1VLA[op->nvv1]=a; /* save atom index for later comparisons */
                          vv2=cs->Coord+(3*a1);
                          vv1=op->vv1+(op->nvv1*3);
                          *(vv1++)=*(vv2++);
                          *(vv1++)=*(vv2++);
                          *(vv1++)=*(vv2++);
                          op->nvv1++;
                        }
                      }
							 break;	
                      /* Moment of inertia tensor - unweighted - assumes v1 is center of molecule */
                    case OMOP_MOME: 
                      if(I->DiscreteFlag) {
                        if(cs==I->DiscreteCSet[a])
                          a1=I->DiscreteAtmToIdx[a];
                        else
                          a1=-1;
                      } else 
                        a1=cs->AtmToIdx[a];
							 if(a1>=0) {
							   subtract3f(cs->Coord+(3*a1),op->v1,v1);
							   v2=v1[0]*v1[0]+v1[1]*v1[1]+v1[2]*v1[2]; 
							   op->d[0][0] += v2 - v1[0] * v1[0];
							   op->d[0][1] +=    - v1[0] * v1[1];
							   op->d[0][2] +=    - v1[0] * v1[2];
							   op->d[1][0] +=    - v1[1] * v1[0];
							   op->d[1][1] += v2 - v1[1] * v1[1];
							   op->d[1][2] +=    - v1[1] * v1[2];
							   op->d[2][0] +=    - v1[2] * v1[0];
							   op->d[2][1] +=    - v1[2] * v1[1];
							   op->d[2][2] += v2 - v1[2] * v1[2];
							 }
							 break;
                    }
                  }
                switch(op->code) { /* full coord-set based */
                case OMOP_INVA:
                  if(inv_flag) {
                    if(op->i1<0) {
                      /* invalidate all representations */
                      for(d=0;d<cRepCnt;d++) {
                        if(cs->fInvalidateRep)
                          cs->fInvalidateRep(cs,d,op->i2);
                      }
                    } else if(cs->fInvalidateRep) 
                      /* invalidate only that particular representation */
                      cs->fInvalidateRep(cs,op->i1,op->i2);
                  }
                  break;
                }
              } /* end coordset section */
          break;
         }
         ai++;
		 }
     break;
	}
	if(hit_flag) {
	  switch(op->code) {
	  case OMOP_TTTF:
       ObjectMoleculeTransformTTTf(I,op->ttt,-1);
       break;
	  case OMOP_LABL:
       ObjectMoleculeInvalidate(I,cRepLabel,cRepInvText);
       break;
     case OMOP_AlterState: /* overly coarse - doing all states, could do just 1 */
       if(!op->i3) { /* not read_only? */
         ObjectMoleculeInvalidate(I,-1,cRepInvRep);
         SceneChanged();
       }
       break;
     case OMOP_CSetIdxSetFlagged:
       ObjectMoleculeInvalidate(I,-1,cRepInvRep);
       SceneChanged();
       break;
     case OMOP_SaveUndo:
       op->i2=true;
       ObjectMoleculeSaveUndo(I,op->i1,false);
       break;
     case OMOP_OnOff:
       ExecutiveSetObjVisib(I->Obj.Name,op->i1);
       break;
	  }
	}

   /* always run on exit...*/
	switch(op->code) {
	case OMOP_ALTR:
   case OMOP_AlterState:
     PUnblock();
     break;
   }
   /* */
  }
}
/*========================================================================*/
void ObjectMoleculeDescribeElement(ObjectMolecule *I,int index, char *buffer) 
{
  AtomInfoType *ai;
  ai=I->AtomInfo+index;
  if(ai->alt[0])
    sprintf(buffer,"%s: /%s/%s/%s/%s/%s`%s",ai->resn,I->Obj.Name,ai->segi,ai->chain,ai->resi,ai->name,ai->alt);
    else
  sprintf(buffer,"%s: /%s/%s/%s/%s/%s",ai->resn,I->Obj.Name,ai->segi,ai->chain,ai->resi,ai->name);
}

/*========================================================================*/
void ObjectMoleculeGetAtomSele(ObjectMolecule *I,int index, char *buffer) 
{
  AtomInfoType *ai;
  ai=I->AtomInfo+index;
  if(ai->alt[0]) 
    sprintf(buffer,"/%s/%s/%s/%s/%s`%s",I->Obj.Name,ai->segi,ai->chain,ai->resi,
            ai->name,ai->alt);
  else
    sprintf(buffer,"/%s/%s/%s/%s/%s`",I->Obj.Name,ai->segi,ai->chain,ai->resi,
            ai->name);   
}
/*========================================================================*/
void ObjectMoleculeGetAtomSeleLog(ObjectMolecule *I,int index, char *buffer) 
{
  AtomInfoType *ai;
  if(SettingGet(cSetting_robust_logs)) {
    ai=I->AtomInfo+index;
    
    if(ai->alt[0]) 
      sprintf(buffer,"/%s/%s/%s/%s/%s`%s",I->Obj.Name,ai->segi,ai->chain,ai->resi,
              ai->name,ai->alt);
    else
      sprintf(buffer,"/%s/%s/%s/%s/%s`",I->Obj.Name,ai->segi,ai->chain,ai->resi,
              ai->name);   
  } else {
    sprintf(buffer,"(%s`%d)",I->Obj.Name,index+1);
  }
}

void ObjectMoleculeGetAtomSeleFast(ObjectMolecule *I,int index, char *buffer) 
{
  AtomInfoType *ai;
  WordType segi,chain,resi,name,alt;
  ai=I->AtomInfo+index;
  
  if(ai->segi[0]) {
    strcpy(segi,"s;");
    strcat(segi,ai->segi);
  } else {
    strcpy(segi,"s;''");
  }
  if(ai->chain[0]) {
    strcpy(chain,"c;");
    strcat(chain,ai->chain);
  } else {
    strcpy(chain,"c;''");
  }
  if(ai->resi[0]) {
    strcpy(resi,"i;");
    strcat(resi,ai->resi);
  } else {
    strcpy(resi,"i;''");
  }
  if(ai->name[0]) {
    strcpy(name,"n;");
    strcat(name,ai->name);
  } else {
    strcpy(name,"n;''");
  }
  if(ai->alt[0]) {
    strcpy(alt,"alt ");
    strcat(alt,ai->alt);
  } else {
    strcpy(alt,"alt ''");
  }
  sprintf(buffer,"(%s&%s&%s&%s&%s&%s)",I->Obj.Name,segi,chain,resi,name,alt);
}

/*========================================================================*/
int ObjectMoleculeGetNFrames(ObjectMolecule *I)
{
  return I->NCSet;
}
/*========================================================================*/
void ObjectMoleculeUpdate(ObjectMolecule *I)
{
  int a;
  OrthoBusyPrime();
  for(a=0;a<I->NCSet;a++)
	 if(I->CSet[a]) {	
	   OrthoBusySlow(a,I->NCSet);
		PRINTFD(FB_ObjectMolecule)
		  " ObjectMolecule-DEBUG: updating state %d of \"%s\".\n" 
         , a+1, I->Obj.Name
        ENDFD;

      if(I->CSet[a]->fUpdate)
        I->CSet[a]->fUpdate(I->CSet[a]);
	 }
  if(I->Obj.RepVis[cRepCell]) {
    if(I->Symmetry) {
      if(I->Symmetry->Crystal) {
        if(I->UnitCellCGO)
          CGOFree(I->UnitCellCGO);
        I->UnitCellCGO = CrystalGetUnitCellCGO(I->Symmetry->Crystal);
      }
    }
  } 
  PRINTFD(FB_ObjectMolecule)
    " ObjectMolecule: updates complete for object %s.\n",I->Obj.Name
    ENDFD;
}
/*========================================================================*/
void ObjectMoleculeInvalidate(ObjectMolecule *I,int rep,int level)
{
  int a;
  PRINTFD(FB_ObjectMolecule)
    " ObjectMoleculeInvalidate: entered. rep: %d level: %d\n",rep,level
    ENDFD;

  if(level>=cRepInvBonds) {
    VLAFreeP(I->Neighbor); /* set I->Neighbor to NULL */
    if(I->Sculpt) {
      SculptFree(I->Sculpt);
      I->Sculpt = NULL;
    }
    ObjectMoleculeUpdateNonbonded(I);
    if(level>=cRepInvAtoms) {
      SelectorUpdateObjectSele(I);
    }
  }
  PRINTFD(FB_ObjectMolecule)
    " ObjectMoleculeInvalidate: invalidating representations...\n"
    ENDFD;

  for(a=0;a<I->NCSet;a++) 
	 if(I->CSet[a]) {	 
      if(I->CSet[a]->fInvalidateRep)
        I->CSet[a]->fInvalidateRep(I->CSet[a],rep,level);
	 }

  PRINTFD(FB_ObjectMolecule)
    " ObjectMoleculeInvalidate: leaving...\n"
    ENDFD;

}
/*========================================================================*/
int ObjectMoleculeMoveAtom(ObjectMolecule *I,int state,int index,float *v,int mode,int log)
{
  int result = 0;
  CoordSet *cs;
  if(!(I->AtomInfo[index].protekted==1)) {
    if(state<0) state=0;
    if(I->NCSet==1) state=0;
    state = state % I->NCSet;
    cs = I->CSet[state];
    if(cs) {
      result = CoordSetMoveAtom(I->CSet[state],index,v,mode);
      cs->fInvalidateRep(cs,cRepAll,cRepInvCoord);
    }
  }
  if(log) {
    OrthoLineType line,buffer;
    if(SettingGet(cSetting_logging)) {
      ObjectMoleculeGetAtomSele(I,index,buffer);
      sprintf(line,"cmd.translate_atom(\"%s\",%15.9f,%15.9f,%15.9f,%d,%d,%d)\n",
              buffer,v[0],v[1],v[2],state+1,mode,0);
      PLog(line,cPLog_no_flush);
    }
  }
  /*  if(I->Sculpt) {
      SculptIterateObject(I->Sculpt,I,state,1);
      }*/
  return(result);
}
/*========================================================================*/
int ObjectMoleculeInitBondPath(ObjectMolecule *I,ObjectMoleculeBPRec *bp )
{
  int a;
  bp->dist = Alloc(int,I->NAtom);
  bp->list = Alloc(int,I->NAtom);
  for(a=0;a<I->NAtom;a++)
    bp->dist[a]=-1;
  bp->n_atom = 0;
  return 1;
}
/*========================================================================*/
int ObjectMoleculePurgeBondPath(ObjectMolecule *I,ObjectMoleculeBPRec *bp )
{
  FreeP(bp->dist);
  FreeP(bp->list);
  return 1;
}
/*========================================================================*/
int ObjectMoleculeGetBondPaths(ObjectMolecule *I,int atom,
                               int max,ObjectMoleculeBPRec *bp)
{
  /* returns list of bond counts from atom to all others 
     dist and list must be vla array pointers or NULL */

  int a,a1,a2,n;
  int cur;
  int n_cur;
  int b_cnt = 0;

  ObjectMoleculeUpdateNeighbors(I);
  
  /* reinitialize dist array (if we've done at least one pass) */

  for(a=0;a<bp->n_atom;a++)
    bp->dist[bp->list[a]]=-1;

  bp->n_atom = 0;
  bp->dist[atom] = 0;
  bp->list[bp->n_atom] = atom;
  bp->n_atom++;
  
  cur = 0;
  while(1) {
    b_cnt++;
    if(b_cnt>max) break;

    n_cur = bp->n_atom-cur;

    /* iterate through all current atoms */

    if(!n_cur) break;
    while(n_cur--) {
      a1 = bp->list[cur++];
      n=I->Neighbor[a1]; 
      n++; /* skip cnt */
      while(1) {
        a2=I->Neighbor[n];
        n+=2;
        if(a2<0) break;
        if(bp->dist[a2]<0) { /* for each atom not yet sampled... */
          bp->dist[a2]=b_cnt;
          bp->list[bp->n_atom]=a2;
          bp->n_atom++;
        }
      }
    }
  }
  return(bp->n_atom);
}
/*========================================================================*/
int ***ObjectMoleculeGetBondPrint(ObjectMolecule *I,int max_bond,int max_type,int *dim)
{
  int a,b,i,c;
  int at1,at2;
  int ***result=NULL;
  ObjectMoleculeBPRec bp;

  dim[0]=max_type+1;
  dim[1]=max_type+1;
  dim[2]=max_bond+1;
  
  result=(int***)UtilArrayMalloc((unsigned int*)dim,3,sizeof(int));
  UtilZeroMem(**result,dim[0]*dim[1]*dim[2]*sizeof(int));
  
  ObjectMoleculeInitBondPath(I,&bp);
  for(a=0;a<I->NAtom;a++) {
    at1 = I->AtomInfo[a].customType;
    if((at1>=0)&&(at1<=max_type)) {
      ObjectMoleculeGetBondPaths(I,a,max_bond,&bp);    
      for(b=0;b<bp.n_atom;b++)
        {
          i = bp.list[b];
          at2 = I->AtomInfo[i].customType;
          if((at2>=0)&&(at2<=max_type)) {
            c=bp.dist[i];
            result[at1][at2][c]++;
          }
        }
    }
  }
  ObjectMoleculePurgeBondPath(I,&bp);
  return(result);
}
/*========================================================================*/
float ObjectMoleculeGetAvgHBondVector(ObjectMolecule *I,int atom,int state,float *v)
     /* computes average hydrogen bonding vector for an atom */
{
  float result = 0.0;
  int a1,a2,n;
  int vec_cnt = 0;
  float v_atom[3],v_neigh[3],v_diff[3],v_acc[3] = {0.0,0.0,0.0};
  CoordSet *cs;

  ObjectMoleculeUpdateNeighbors(I);

  a1 = atom;
  if(state<0) state=0;
  if(I->NCSet==1) state=0;
  state = state % I->NCSet;
  cs = I->CSet[state];
  if(cs) {
    if(CoordSetGetAtomVertex(cs,a1,v_atom)) { /* atom exists in this C-set */
      n=I->Neighbor[atom];
      n++;
      while(1) {
        a2=I->Neighbor[n];
        if(a2<0) break;
        n+=2;
        
        if(I->AtomInfo[a2].protons!=1) { /* ignore hydrogens */
          if(CoordSetGetAtomVertex(cs,a2,v_neigh)) { 
            subtract3f(v_atom,v_neigh,v_diff);
            normalize3f(v_diff);
            add3f(v_diff,v_acc,v_acc);
            vec_cnt++;
          }
        }
      }
      if(vec_cnt) {
        result = (float)length3f(v_acc);
        result = result/vec_cnt;
        normalize23f(v_acc,v);
      }
      copy3f(v_acc,v);
    }
  }
  return(result);
}
/*========================================================================*/
int ObjectMoleculeGetAtomVertex(ObjectMolecule *I,int state,int index,float *v)
{
  int result = 0;
  if(state<0) state=SettingGet_i(NULL,I->Obj.Setting,cSetting_state)-1;
  if(state<0) state=SceneGetState(); 
  if(I->NCSet==1) state=0;
  state = state % I->NCSet;
  if(I->CSet[state]) 
    result = CoordSetGetAtomVertex(I->CSet[state],index,v);
  return(result);
}
/*========================================================================*/
int ObjectMoleculeSetAtomVertex(ObjectMolecule *I,int state,int index,float *v)
{
  int result = 0;
  if(state<0) state=SettingGet_i(NULL,I->Obj.Setting,cSetting_state)-1;
  if(state<0) state=SceneGetState();
  if(I->NCSet==1) state=0;
  state = state % I->NCSet;
  if(I->CSet[state]) 
    result = CoordSetSetAtomVertex(I->CSet[state],index,v);
  return(result);
}
/*========================================================================*/
void ObjectMoleculeRender(ObjectMolecule *I,int state,CRay *ray,Pickable **pick,int pass)
{
  int a;

  PRINTFD(FB_ObjectMolecule)
    " ObjectMolecule: rendering %s...\n",I->Obj.Name
    ENDFD;

  ObjectPrepareContext(&I->Obj,ray);

  if(I->UnitCellCGO&&(I->Obj.RepVis[cRepCell])) {
    if(ray) {
      
      CGORenderRay(I->UnitCellCGO,ray,ColorGet(I->Obj.Color),
                         I->Obj.Setting,NULL);
    } else if(pick&&PMGUI) {
    } else if(PMGUI) {
      ObjectUseColor(&I->Obj);
      CGORenderGL(I->UnitCellCGO,ColorGet(I->Obj.Color),
                         I->Obj.Setting,NULL);
    }
  }

  PRINTFD(FB_ObjectMolecule)
    " ObjectMolecule: CGO's complete...\n"
    ENDFD;
  if(state<0) {
    for(a=0;a<I->NCSet;a++)
      if(I->CSet[a])
        if(I->CSet[a]->fRender)
          I->CSet[a]->fRender(I->CSet[a],ray,pick,pass);        
  } else if(state<I->NCSet) {
	 I->CurCSet=state % I->NCSet;
	 if(I->CSet[I->CurCSet]) {
      if(I->CSet[I->CurCSet]->fRender)
        I->CSet[I->CurCSet]->fRender(I->CSet[I->CurCSet],ray,pick,pass);
	 }
  } else if(I->NCSet==1) { /* if only one coordinate set, assume static */
    if(SettingGet(cSetting_static_singletons))
      if(I->CSet[0]->fRender)
        I->CSet[0]->fRender(I->CSet[0],ray,pick,pass);    
  }
  PRINTFD(FB_ObjectMolecule)
    " ObjectMolecule: rendering complete for object %s.\n",I->Obj.Name
    ENDFD;
}
/*========================================================================*/
void ObjectMoleculeDummyUpdate(ObjectMolecule *I,int mode)
{
  switch(mode) {
  case cObjectMoleculeDummyOrigin:
    SceneOriginGet(I->CSet[0]->Coord);
    break;
  case cObjectMoleculeDummyCenter:
    SceneGetPos(I->CSet[0]->Coord);
    break;
  }
}
/*========================================================================*/
ObjectMolecule *ObjectMoleculeDummyNew(int type)
{
  ObjectMolecule *I= NULL;
  
  int nAtom;
  float *coord = NULL;
  CoordSet *cset = NULL;
  AtomInfoType *atInfo = NULL;
  int frame=-1;

  I=ObjectMoleculeNew(false);

  nAtom=1;
  coord=VLAlloc(float,3*nAtom);
  zero3f(coord);
  
  atInfo=VLAMalloc(10,sizeof(AtomInfoType),2,true); /* autozero here is important */
  
  cset = CoordSetNew();
  cset->NIndex=nAtom;
  cset->Coord=coord;
  cset->TmpBond=NULL;
  cset->NTmpBond=0;
  strcpy(cset->Name,"_origin");
  
  cset->Obj = I;
  cset->fEnumIndices(cset);
  
  ObjectMoleculeMerge(I,atInfo,cset,false,cAIC_IDMask); /* NOTE: will release atInfo */
  
  if(frame<0) frame=I->NCSet;
  VLACheck(I->CSet,CoordSet*,frame);
  if(I->NCSet<=frame) I->NCSet=frame+1;
  if(I->CSet[frame]) I->CSet[frame]->fFree(I->CSet[frame]);
  I->CSet[frame] = cset;

  I->NBond = 0;
  I->Bond = VLAlloc(BondType,0);
  
  ObjectMoleculeExtendIndices(I);
  ObjectMoleculeSort(I);
  ObjectMoleculeUpdateIDNumbers(I);
  ObjectMoleculeUpdateNonbonded(I);

  return(I);
}

/*========================================================================*/

ObjectMolecule *ObjectMoleculeNew(int discreteFlag)
{
  int a;
  OOAlloc(ObjectMolecule);
  ObjectInit((CObject*)I);
  I->Obj.type=cObjectMolecule;
  I->NAtom=0;
  I->NBond=0;
  I->CSet=VLAMalloc(10,sizeof(CoordSet*),5,true); /* auto-zero */
  I->NCSet=0;
  I->Bond=NULL;
  I->AtomCounter=-1;
  I->BondCounter=-1;
  I->DiscreteFlag=discreteFlag;
  I->NDiscrete=0;
  I->UnitCellCGO=NULL;
  I->Sculpt=NULL;
  I->CSTmpl=NULL;
  if(I->DiscreteFlag) { /* discrete objects don't share atoms between states */
    I->DiscreteAtmToIdx = VLAMalloc(10,sizeof(int),6,false);
    I->DiscreteCSet = VLAMalloc(10,sizeof(CoordSet*),5,false);
    I->NDiscrete=0;
  } else {
    I->DiscreteAtmToIdx = NULL;
    I->DiscreteCSet = NULL;
  }    
  I->Obj.fRender=(void (*)(struct CObject *, int, CRay *, Pickable **,int))ObjectMoleculeRender;
  I->Obj.fFree= (void (*)(struct CObject *))ObjectMoleculeFree;
  I->Obj.fUpdate=  (void (*)(struct CObject *)) ObjectMoleculeUpdate;
  I->Obj.fGetNFrame = (int (*)(struct CObject *)) ObjectMoleculeGetNFrames;
  I->Obj.fDescribeElement = (void (*)(struct CObject *,int index,char *buffer)) ObjectMoleculeDescribeElement;
  I->Obj.fGetSettingHandle = (CSetting **(*)(struct CObject *,int state))
    ObjectMoleculeGetSettingHandle;
  I->AtomInfo=VLAMalloc(10,sizeof(AtomInfoType),2,true); /* autozero here is important */
  I->CurCSet=0;
  I->Symmetry=NULL;
  I->Neighbor=NULL;
  for(a=0;a<=cUndoMask;a++) {
    I->UndoCoord[a]=NULL;
    I->UndoState[a]=-1;
  }
  I->UndoIter=0;
  return(I);
}
/*========================================================================*/
ObjectMolecule *ObjectMoleculeCopy(ObjectMolecule *obj)
{
  int a;
  BondType *i0,*i1;
  AtomInfoType *a0,*a1;
  OOAlloc(ObjectMolecule);
  (*I)=(*obj);
  I->Symmetry=SymmetryCopy(I->Symmetry); /* null-safe */
  I->UnitCellCGO=NULL;
  I->Neighbor=NULL;
  I->Sculpt=NULL;
  for(a=0;a<=cUndoMask;a++)
    I->UndoCoord[a]=NULL;
  I->CSet=VLAMalloc(I->NCSet,sizeof(CoordSet*),5,true); /* auto-zero */
  for(a=0;a<I->NCSet;a++) {
    I->CSet[a]=CoordSetCopy(obj->CSet[a]);
    I->CSet[a]->Obj=I;
  }
  if(obj->CSTmpl)
    I->CSTmpl = CoordSetCopy(obj->CSTmpl);
  else
    I->CSTmpl=NULL;
  I->Bond=VLAlloc(BondType,I->NBond);
  i0=I->Bond;
  i1=obj->Bond;
  for(a=0;a<I->NBond;a++) {
    *(i0++)=*(i1++); /* copy structure */
  }
  
  I->AtomInfo=VLAlloc(AtomInfoType,I->NAtom);
  a0=I->AtomInfo;
  a1=obj->AtomInfo;
  for(a=0;a<I->NAtom;a++)
    *(a0++)=*(a1++);

  for(a=0;a<I->NAtom;a++) {
    I->AtomInfo[a].selEntry=0;
  }
  
  return(I);

}

/*========================================================================*/
void ObjectMoleculeFree(ObjectMolecule *I)
{
  int a;

  SceneObjectDel((CObject*)I);
  SelectorPurgeObjectMembers(I);
  for(a=0;a<I->NCSet;a++)
	 if(I->CSet[a]) {
      if(I->CSet[a]->fFree)
        I->CSet[a]->fFree(I->CSet[a]);
		I->CSet[a]=NULL;
	 }
  if(I->Symmetry) SymmetryFree(I->Symmetry);
  VLAFreeP(I->Neighbor);
  VLAFreeP(I->DiscreteAtmToIdx);
  VLAFreeP(I->DiscreteCSet);
  VLAFreeP(I->CSet);
  VLAFreeP(I->AtomInfo);
  VLAFreeP(I->Bond);
  if(I->UnitCellCGO) 
    CGOFree(I->UnitCellCGO);
  for(a=0;a<=cUndoMask;a++)
    FreeP(I->UndoCoord[a]);
  if(I->Sculpt)
    SculptFree(I->Sculpt);
  if(I->CSTmpl)
    if(I->CSTmpl->fFree)
      I->CSTmpl->fFree(I->CSTmpl);
  ObjectPurge(&I->Obj);
  OOFreeP(I);
}

/*========================================================================*/
ObjectMolecule *ObjectMoleculeReadMMDStr(ObjectMolecule *I,char *MMDStr,int frame,int discrete)
{
  int ok = true;
  CoordSet *cset=NULL;
  AtomInfoType *atInfo;
  int isNew;
  int nAtom;

  if(!I) 
	 isNew=true;
  else 
	 isNew=false;

  if(isNew) {
    I=(ObjectMolecule*)ObjectMoleculeNew(discrete);
    atInfo = I->AtomInfo;
  } else {
    atInfo=VLAMalloc(10,sizeof(AtomInfoType),2,true); /* autozero here is important */
  }
  
  if(isNew) {
    AtomInfoPrimeColors();
    I->Obj.Color = AtomInfoGetCarbColor();
  }

  cset=ObjectMoleculeMMDStr2CoordSet(MMDStr,&atInfo);  

  if(!cset) 
	 {
		VLAFreeP(atInfo);
		ok=false;
	 }
  
  if(ok)
	 {
		if(!I) 
		  I=(ObjectMolecule*)ObjectMoleculeNew(discrete);
		if(frame<0)
		  frame=I->NCSet;
		if(I->NCSet<=frame)
		  I->NCSet=frame+1;
		VLACheck(I->CSet,CoordSet*,frame);
      nAtom=cset->NIndex;
      cset->Obj = I;
      if(cset->fEnumIndices)
        cset->fEnumIndices(cset);
      if(cset->fInvalidateRep)
        cset->fInvalidateRep(cset,cRepAll,cRepInvRep);
      if(isNew) {		
        I->AtomInfo=atInfo; /* IMPORTANT to reassign: this VLA may have moved! */
        I->NAtom=nAtom;
      } else {
        ObjectMoleculeMerge(I,atInfo,cset,false,cAIC_MMDMask); /* NOTE: will release atInfo */
      }
      if(frame<0) frame=I->NCSet;
      VLACheck(I->CSet,CoordSet*,frame);
      if(I->NCSet<=frame) I->NCSet=frame+1;
      I->CSet[frame] = cset;
      if(isNew) I->NBond = ObjectMoleculeConnect(I,&I->Bond,I->AtomInfo,cset,false);
      SceneCountFrames();
      ObjectMoleculeExtendIndices(I);
      ObjectMoleculeSort(I);
      ObjectMoleculeUpdateIDNumbers(I);
      ObjectMoleculeUpdateNonbonded(I);
	 }
  return(I);
}
/*========================================================================*/
ObjectMolecule *ObjectMoleculeLoadMMDFile(ObjectMolecule *obj,char *fname,
                                          int frame,char *sepPrefix,int discrete)
{
  ObjectMolecule* I=NULL;
  int ok=true;
  FILE *f;
  int oCnt=0;
  long size;
  char *buffer,*p;
  char cc[MAXLINELEN],oName[ObjNameMax];
  int nLines;
  f=fopen(fname,"rb");
  if(!f)
	 ok=ErrMessage("ObjectMoleculeLoadMMDFile","Unable to open file!");
  else
	 {
      PRINTFB(FB_ObjectMolecule,FB_Blather)
        " ObjectMoleculeLoadMMDFile: Loading from %s.\n",fname
        ENDFB;
		fseek(f,0,SEEK_END);
      size=ftell(f);
		fseek(f,0,SEEK_SET);
		buffer=(char*)mmalloc(size+255);
		ErrChkPtr(buffer);
		p=buffer;
		fseek(f,0,SEEK_SET);
		fread(p,size,1,f);
		p[size]=0;
		fclose(f);
      p=buffer;
      while(ok) {
        ncopy(cc,p,6);
        if(sscanf(cc,"%d",&nLines)!=1)
          break;
        if(ok) {
          if(sepPrefix) {
            I=ObjectMoleculeReadMMDStr(NULL,p,frame,discrete);
            oCnt++;
            sprintf(oName,"%s-%02d",sepPrefix,oCnt);
            ObjectSetName((CObject*)I,oName);
            ExecutiveManageObject((CObject*)I,true,false);
          } else {
            I=ObjectMoleculeReadMMDStr(obj,p,frame,discrete);
            obj=I;
          }
          p=nextline(p);
          while(nLines--)
            p=nextline(p);
        }
      }
		mfree(buffer);
	 }

  return(I);
}

/*========================================================================*/
ObjectMolecule *ObjectMoleculeReadPDBStr(ObjectMolecule *I,char *PDBStr,int frame,
                                         int discrete,M4XAnnoType *m4x,char *pdb_name,
                                         char **next_pdb)
{
  CoordSet *cset = NULL;
  AtomInfoType *atInfo;
  int ok=true;
  int isNew = true;
  unsigned int nAtom = 0;
  char *start,*restart=NULL;
  int repeatFlag = true;
  int successCnt = 0;
  unsigned int aic_mask = cAIC_PDBMask;

  SegIdent segi_override=""; /* saved segi for corrupted NMR pdb files */

  start=PDBStr;
  while(repeatFlag) {
    repeatFlag = false;
  
    if(!I) 
      isNew=true;
    else 
      isNew=false;
    
    if(ok) {
      
      if(isNew) {
        I=(ObjectMolecule*)ObjectMoleculeNew(discrete);
        atInfo = I->AtomInfo;
        isNew = true;
      } else {
        atInfo=VLAMalloc(10,sizeof(AtomInfoType),2,true); /* autozero here is important */
        isNew = false;
      }
      if(isNew) {
        AtomInfoPrimeColors();
        I->Obj.Color = AtomInfoGetCarbColor();
      }

      cset=ObjectMoleculePDBStr2CoordSet(start,&atInfo,&restart,
                                         segi_override,m4x,pdb_name,next_pdb);	
      if(m4x) /* preserve original atom IDs for annotated Metaphorics files */
        if(m4x->annotated_flag)
          aic_mask = (cAIC_b|cAIC_q);
      nAtom=cset->NIndex;
    }
    if(pdb_name&&(*next_pdb)) {
      /* problematic scenario */
    }

    /* include coordinate set */
    if(ok) {

      cset->Obj = I;
      cset->fEnumIndices(cset);
      if(cset->fInvalidateRep)
        cset->fInvalidateRep(cset,cRepAll,cRepInvRep);
      if(isNew) {		
        I->AtomInfo=atInfo; /* IMPORTANT to reassign: this VLA may have moved! */
      } else {
        ObjectMoleculeMerge(I,atInfo,cset,true,aic_mask); /* NOTE: will release atInfo */
      }
      if(isNew) I->NAtom=nAtom;
      if(frame<0) frame=I->NCSet;
      VLACheck(I->CSet,CoordSet*,frame);
      if(I->NCSet<=frame) I->NCSet=frame+1;
      if(I->CSet[frame]) I->CSet[frame]->fFree(I->CSet[frame]);
      I->CSet[frame] = cset;
      if(isNew) I->NBond = ObjectMoleculeConnect(I,&I->Bond,I->AtomInfo,cset,true);
      if(cset->Symmetry&&(!I->Symmetry)) {
        I->Symmetry=SymmetryCopy(cset->Symmetry);
        SymmetryAttemptGeneration(I->Symmetry,false,false);
      }
      SceneCountFrames();
      ObjectMoleculeExtendIndices(I);
      ObjectMoleculeSort(I);
      ObjectMoleculeUpdateIDNumbers(I);
      ObjectMoleculeUpdateNonbonded(I);
      successCnt++;
      if(successCnt>1) {
        if(successCnt==2){
          PRINTFB(FB_ObjectMolecule,FB_Actions)
            " ObjectMolReadPDBStr: read MODEL %d\n",1
            ENDFB;
            }
        PRINTFB(FB_ObjectMolecule,FB_Actions)
          " ObjectMolReadPDBStr: read MODEL %d\n",successCnt
          ENDFB;
      }
    }
    if(restart) {
      repeatFlag=true;
      start=restart;
      frame=frame+1;
    }
  }
  return(I);
}
/*========================================================================*/
CoordSet *ObjectMoleculeMMDStr2CoordSet(char *buffer,AtomInfoType **atInfoPtr)
{
  char *p;
  int nAtom,nBond;
  int a,c,bPart,bOrder;
  float *coord = NULL;
  CoordSet *cset = NULL;
  AtomInfoType *atInfo = NULL,*ai;
  char cc[MAXLINELEN];
  float *f;
  BondType *ii,*bond=NULL;
  int ok=true;
  int auto_show_lines = (int)SettingGet(cSetting_auto_show_lines);
  int auto_show_nonbonded = (int)SettingGet(cSetting_auto_show_nonbonded);

  p=buffer;
  nAtom=0;
  if(atInfoPtr)
	 atInfo = *atInfoPtr;


  if(ok) {
	 p=ncopy(cc,p,6);
	 if(sscanf(cc,"%d",&nAtom)!=1)
		ok=ErrMessage("ReadMMDFile","bad atom count");
  }

  if(ok) {
	 coord=VLAlloc(float,3*nAtom);
	 if(atInfo)
		VLACheck(atInfo,AtomInfoType,nAtom);	 
  }

  if(!atInfo)
    ErrFatal("PDBStr2CoordSet","need atom information record!"); /* failsafe for old version..*/

  nBond=0;
  if(ok) {
	 bond=VLAlloc(BondType,6*nAtom);  
  }
  p=nextline(p);

  /* read coordinates and atom names */

  if(ok) { 
	 f=coord;
	 ii=bond;
	 for(a=0;a<nAtom;a++)
		{
        ai=atInfo+a;

        ai->id=a+1;
        if(ok) {
          p=ncopy(cc,p,4);
          if(sscanf(cc,"%d",&ai->customType)!=1) 
            ok=ErrMessage("ReadMMDFile","bad atom type");
        }
        if(ok) {
          if(ai->customType<=14) strcpy(ai->elem,"C");
          else if(ai->customType<=23) strcpy(ai->elem,"O");
          else if(ai->customType<=40) strcpy(ai->elem,"N");
          else if(ai->customType<=48) strcpy(ai->elem,"H");
          else if(ai->customType<=52) strcpy(ai->elem,"S");
          else if(ai->customType<=53) strcpy(ai->elem,"P");
          else if(ai->customType<=55) strcpy(ai->elem,"B");
          else if(ai->customType<=56) strcpy(ai->elem,"F");
          else if(ai->customType<=57) strcpy(ai->elem,"Cl");           
          else if(ai->customType<=58) strcpy(ai->elem,"Br");           
          else if(ai->customType<=59) strcpy(ai->elem,"I");           
          else if(ai->customType<=60) strcpy(ai->elem,"Si");           
          else if(ai->customType<=61) strcpy(ai->elem,"Du");           
          else if(ai->customType<=62) strcpy(ai->elem,"Z0");
          else if(ai->customType<=63) strcpy(ai->elem,"Lp");
          else strcpy(ai->elem,"?");
        }
        for(c=0;c<6;c++) {
          if(ok) {
            p=ncopy(cc,p,8);
            if(sscanf(cc,"%d%d",&bPart,&bOrder)!=2)
              ok=ErrMessage("ReadMMDFile","bad bond record");
            else {
              if(bPart&&bOrder&&(a<(bPart-1))) {
                nBond++;
                ii->index[0]=a;
                ii->index[1]=bPart-1;
                ii->order=bOrder;
                ii->stereo=0;
                ii->id=-1;
                ii++;
              }
            }
          }
        }
        if(ok) {
          p=ncopy(cc,p,12);
          if(sscanf(cc,"%f",f++)!=1)
            ok=ErrMessage("ReadMMDFile","bad coordinate");
        }
        if(ok) {
          p=ncopy(cc,p,12);
          if(sscanf(cc,"%f",f++)!=1)
            ok=ErrMessage("ReadMMDFile","bad coordinate");
        }
        if(ok) {
          p=ncopy(cc,p,12);
			 if(sscanf(cc,"%f",f++)!=1)
				ok=ErrMessage("ReadMMDFile","bad coordinate");
		  }
        if(ok) {
          p=nskip(p,1);
          p=ncopy(cc,p,5);
          if(sscanf(cc,"%d",&ai->resv)==1) {
            sprintf(ai->resi,"%d",ai->resv); /* range check...*/
          }
        }
        if(ok) {
          p=nskip(p,6);
          p=ncopy(cc,p,9);
			 if(sscanf(cc,"%f",&ai->partialCharge)!=1)
				ok=ErrMessage("ReadMMDFile","bad charge");
        }
        if(ok) {
          p=nskip(p,10);
          p=ncopy(cc,p,3);
          if(sscanf(cc,"%s",ai->resn)!=1)
            ai->resn[0]=0;
          ai->hetatm=true;
        }

        ai->segi[0]=0;
        ai->alt[0]=0;

        if(ok) {
          p=nskip(p,2);
          p=ncopy(ai->name,p,4);
          UtilCleanStr(ai->name);
          if(ai->name[0]==0) {
            strcpy(ai->name,ai->elem);
            sprintf(cc,"%02d",a+1);
            if((strlen(cc)+strlen(ai->name))>4)
              strcpy(ai->name,cc);
            else
              strcat(ai->name,cc);
          }

          for(c=0;c<cRepCnt;c++) {
            ai->visRep[c] = false;
          }
          ai->visRep[cRepLine] = auto_show_lines; /* show lines by default */
          ai->visRep[cRepNonbonded] = auto_show_nonbonded; /* show lines by default */

        }
        if(ok) {
          AtomInfoAssignParameters(ai);
          ai->color = AtomInfoGetColor(ai);
        }
        if(!ok)
          break;
        p=nextline(p);
      }
  }
  if(ok)
    VLASize(bond,BondType,nBond);
  if(ok) {
	 cset = CoordSetNew();
	 cset->NIndex=nAtom;
	 cset->Coord=coord;
	 cset->NTmpBond=nBond;
	 cset->TmpBond=bond;
  } else {
	 VLAFreeP(bond);
	 VLAFreeP(coord);
  }
  if(atInfoPtr)
	 *atInfoPtr = atInfo;
  return(cset);
}
