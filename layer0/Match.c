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
#include"os_limits.h"
#include"os_std.h"

#include"Base.h"
#include"MemoryDebug.h"
#include"OOMac.h"
#include"Match.h"
#include"Util.h"
#include"Feedback.h"
#include"Parse.h"

#ifndef int2 
typedef int int2[2];
#endif

CMatch *MatchNew(unsigned int na,unsigned int nb)
{
  unsigned int dim[2];
  int a,b;
  OOAlloc(CMatch);

  dim[0]=na;
  dim[1]=nb;
  I->mat = NULL;
  if(na&&nb) {
    I->mat = (float**)UtilArrayMalloc(dim,2,sizeof(float));
  }
  I->pair = NULL;
  I->na=na;
  I->nb=nb;

  /* scoring matrix is always 128^2 */
  dim[0]=128;
  dim[1]=128;
  I->smat = (float**)UtilArrayMalloc(dim,2,sizeof(float));
  for(a=0;a<128;a++)
    for(b=0;b<128;b++) 
      I->smat[a][b]=0.0;
  return(I);
}

int MatchResidueToCode(CMatch *I,int *vla,int n)
{
#define cNRES 20
  int ok=true;
  int a,b,c;
  int found;
  int rcode[cNRES],rname[cNRES];
  int *trg;
  char res[][4] = { "ALA", "A",
                     "CYS", "C",
                     "ASP", "D",
                     "GLU", "E",
                     "PHE", "F",
                     
                     "GLY", "G",
                     "HIS", "H",
                     "ILE", "I",
                     "LYS", "K",
                     "LEU", "L",
                     
                     "MET", "M",
                     "ASN", "N",
                     "PRO", "P",
                     "GLN", "Q",
                     "ARG", "R",
                     
                     "SER", "S",
                     "THR", "T",
                     "VAL", "V",
                     "TRP", "W",
                     "TYR", "Y" };
  
  /* get integral values for the residue names */
  
  for(a=0;a<cNRES;a++) {
    b = 0;
    for(c=0;c<3;c++) {
      b = (b<<8) | res[a*2][c];
    }
    rname[a]=b;
    rcode[a]=res[a*2+1][0];
  }
  
  for(b=0;b<n;b++) {
    found = 0;
    trg = vla+(b*3)+2;
    for(a=0;a<cNRES;a++)
      if(rname[a]==*trg) {
        found=true;
        *trg=rcode[a];
        break;
      }
    if(!found) {
      PRINTFB(FB_Match,FB_Warnings)
        " Match-Warning: unknown residue type %c%c%c (using X).\n",
        0xFF&(*trg>>16),0xFF&(*trg>>8),0xFF&*trg
        ENDFB;
      *trg='X';
    }
  }
  return(ok);
}

int MatchPreScore(CMatch *I,int *vla1,int n1,int *vla2,int n2)
{
  int a,b;

  PRINTFB(FB_Match,FB_Details)
    " Match: assigning %d x %d pairwise scores.\n",n1,n2
    ENDFB;

  for(a=0;a<n1;a++) {
    for(b=0;b<n2;b++) {
      I->mat[a][b] = I->smat[0x7F&vla1[a*3+2]][0x7F&vla2[b*3+2]];
    }
  }
  return 1;
}


int MatchMatrixFromFile(CMatch *I,char *fname)
{
  int ok=1;
  FILE *f;
  char *buffer;
  char *p;
  char cc[255];
  char *code=NULL;
  unsigned int x,y;
  int a;
  int n_entry;
  unsigned int size;
  
  f=fopen(fname,"r");
  if(!f) {
    PRINTFB(FB_Match,FB_Errors) 
      " Match-Error: unable to open matrix file '%s'.\n",fname
      ENDFB;
    ok=false;
  } else {
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

    /* count codes */

    p=buffer;
    n_entry = 0;
    while(*p) {
      switch(*p) {
      case '#':
        break;
      default:
        if((*p)>32) n_entry++;
        break;
      }
      p=ParseNextLine(p);
    }

    if(!n_entry) 
      ok=false;
    else {
      code=(char*)Alloc(int,n_entry);
      
      /* read codes */
      
      p=buffer;
      n_entry = 0;
      while(*p) {
        switch(*p) {
        case '#':
          break;
        default:
          if((*p)>32) {
            code[n_entry]=*p;
            n_entry++;
            break;
          }
        }
        p=ParseNextLine(p);
      }
      
      /* read values */

      p=buffer;
      while(*p) {
        switch(*p) {
        case '#':
          break;
        default:
          if((*p)>32) {
            x=*(p++);
            for(a=0;a<n_entry;a++) {
              p = ParseWordCopy(cc,p,255);
              y=(unsigned int)code[a];
              ok=sscanf(cc,"%f",&I->smat[x][y]);
              if(!ok) break;
            }
          }
          break;
        }
        if(!ok) break;
        p=ParseNextLine(p);
      }
      mfree(buffer);
    }
  }
  if(ok) {
    PRINTFB(FB_Match,FB_Details)
      " Match: read scoring matrix.\n"
      ENDFB;
  }
  FreeP(code);
  return(ok);
}

float MatchAlign(CMatch *I,float gap_penalty,float ext_penalty,int max_skip)
{
  unsigned int dim[2];
  int a,b,f,g;
  int nf,ng;
  int sf,sg;
  float **score;
  int2 **point;
  float mxv;
  int mxa,mxb;
  float tst=0.0;
  int gap=0;
  int *p;
  int cnt;

  nf = I->na+2;
  ng = I->nb+2;

  PRINTFB(FB_Match,FB_Actions)
    " MatchAlign: aligning residues (%d vs %d)...\n",I->na,I->nb
    ENDFB;

  dim[0]=nf;
  dim[1]=ng;
  VLAFreeP(I->pair);
  score = (float**)UtilArrayMalloc(dim,2,sizeof(float));
  point = (int2**)UtilArrayMalloc(dim,2,sizeof(int2));
  /* initialize the scoring matrix */
  for(f=0;f<nf;f++) {
    for(g=0;g<ng;g++) {
      score[f][g] = 0.0;
    }
  }
  /* now start walking backwards up the alignment */

  for(b=I->nb-1;b>=0;b--) {
    for(a=I->na-1;a>=0;a--) {

      /* find the maximum scoring cell accessible from this position, 
       * while taking gap penalties into account */

      mxv = FLT_MIN;
      mxa=-1;
      mxb=-1;

      /* search for asymmetric insertions and deletions */
      f = a+1;
      for(g=b+1;g<ng;g++) {
        tst = score[f][g];
        if(!((f==I->na)||(g==I->nb))) {
          gap = g-(b+1);
          if(gap) tst+=gap_penalty+ext_penalty*gap;
        }
        if(tst>mxv) {
          mxv = tst;
          mxa = f;
          mxb = g;
        }
      }
      g = b+1;
      for(f=a+1;f<nf;f++) {
        tst = score[f][g];
        if(!((f==I->na)||(g==I->nb))) {
          gap=(f-(a+1));
          if(gap) tst+=gap_penalty+ext_penalty*gap;
        }
        if(tst>mxv) {
          mxv = tst;
          mxa = f;
          mxb = g;
        }
      }

      /* search for high scoring mismatched stretches */

      sf = a+1+max_skip;
      sg = b+1+max_skip;
      if(sf>nf) sf = nf;
      if(sg>ng) sg = ng;

      for(f=a+1;f<sf;f++) {
        for(g=b+1;g<sg;g++) {
          tst = score[f][g];
          /* only penalize if we are not at the end */
          if(!((f==I->na)||(g==I->nb)))
            gap = ((f-(a+1))+(g-(b+1)));
            tst+=2*gap_penalty+ext_penalty*gap;
          }
          if(tst>mxv) {
            mxv = tst;
            mxa = f;
            mxb = g;
          }
      }
      
      /* store what the best next step is */
    
      point[a][b][0] = mxa;
      point[a][b][1] = mxb;
      
      /* and store the cumulative score for this cell */
      score[a][b] = mxv+I->mat[a][b];
      
  }
}

  if(Feedback(FB_Match,FB_Debugging)) {
    for(b=0;b<I->nb;b++) {
      for(a=0;a<I->na;a++) {
        printf("%4.1f(%2d,%2d)",score[a][b],point[a][b][0],point[a][b][1]);
      }
      printf("\n");
    }
  }

  /* find the best entry point */

  mxv = FLT_MIN;
  mxa=0;
  mxb=0;
  for(b=0;b<I->nb;b++) {
    for(a=0;a<I->na;a++) {
      tst = score[a][b];
      if(tst>mxv) {
        mxv = tst;
        mxa = a;
        mxb = b;
      }
    }
  }
  I->pair = VLAlloc(int,2*(I->na>I->nb?I->na:I->nb));
  p=I->pair;
  a = mxa;
  b = mxb;
  cnt=0;
  while((a>=0)&&(b>=0)) {
    *(p++)=a;
    *(p++)=b;
    f = point[a][b][0];
    g = point[a][b][1];
    a=f;
    b=g;
    cnt++;
  }
  PRINTFD(FB_Match)
    " MatchAlign-DEBUG: best entry %8.3f %d %d %d\n",mxv,mxa,mxb,cnt
    ENDFD;
  PRINTFB(FB_Match,FB_Results) 
    " MatchAlign: score %1.3f\n",mxv
    ENDFD;
  if(cnt)
    mxv = mxv/cnt;
  VLASize(I->pair,int,(p-I->pair));
  FreeP(score);
  FreeP(point);
  return(mxv);
}

void MatchFree(CMatch *I)
{
  FreeP(I->mat);
  FreeP(I->smat);
  VLAFreeP(I->pair);
  OOFreeP(I);
}


