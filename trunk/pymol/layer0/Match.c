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

#include"os_limits.h"
#include"os_std.h"

#include"MemoryDebug.h"
#include"OOMac.h"
#include"Match.h"
#include"Util.h"
#include"Feedback.h"

#ifndef int2 
typedef int int2[2];
#endif

CMatch *MatchNew(unsigned int na,unsigned int nb)
{
  unsigned int dim[2];
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
  return(I);
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
  float tst;
  int gap;
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

      /* search for gross insertions and deletions */
      f = a+1;
      for(g=b+1;g<ng;g++) {
        tst = score[f][g];
        /* only penalize if we are not at the end */
        if(tst!=0.0) {
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
          /* only penalize if we are not at the end */
          gap=(f-(a+1));
          if(tst!=0.0) {
            gap = f-(a+1);
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
          if(tst!=0.0) {
            gap = ((f-(a+1))+(g-(b+1)));
            tst+=2*gap_penalty+ext_penalty*gap;
          }
          if(tst>mxv) {
            mxv = tst;
            mxa = f;
            mxb = g;
          }
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
  VLAFreeP(I->pair);
  OOFreeP(I);
}


