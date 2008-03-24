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

CMatch *MatchNew(PyMOLGlobals *G,unsigned int na,unsigned int nb,int dist_mats)
{
  unsigned int dim[2];
  OOCalloc(G,CMatch);

  I->na=na;
  I->nb=nb;

  dim[0]=na;
  dim[1]=nb;
  I->G=G;
  if(na&&nb) {
    I->mat = (float**)UtilArrayCalloc(dim,2,sizeof(float));
  }

  if(dist_mats && na) {
    dim[0] = na+1;
    dim[1] = na+1;
    I->da = (float**)UtilArrayCalloc(dim,2,sizeof(float));
  }
  if(dist_mats && nb) {
    dim[0] = nb+1;
    dim[1] = nb+1;
    I->db = (float**)UtilArrayCalloc(dim,2,sizeof(float));
  }

  /* scoring matrix is always 128^2 */
  dim[0]=128;
  dim[1]=128;
  I->smat = (float**)UtilArrayCalloc(dim,2,sizeof(float));

  if(!(I->mat && I->smat && ((!dist_mats)||(I->da&&I->db)))) {
    MatchFree(I);
    I=NULL;
  }
  return(I);
}

int MatchResidueToCode(CMatch *I,int *vla,int n)
{

#define cNRES 33
  PyMOLGlobals *G=I->G;
  int ok=true;
  int a,b,c;
  int found;
  int rcode[cNRES],rname[cNRES];
  int *trg;
  char res[][4] = { 
    "A"  , "A",
    "ADE", "A",
    "C"  , "C",
    "CYT", "C",
    "G"  , "G",

    "GUA", "G",
    "T"  , "T",
    "THY", "T",
    "U"  , "T",
    "URA", "T",

    "ALA", "A",
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
    "SEP", "S", /* phosphoserine */
    "THR", "T",
    "TPO", "T", /* phosphothreonine */
    "VAL", "V",

    "TRP", "W",
    "TYR", "Y",
    "PTR", "Y"  /* phosphotyrosine */
};
  
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
      if(*trg!=0x484f48) { /* HOH */
        char st[4];
        st[0] = 0xFF&(*trg>>16);
        st[1] = 0xFF&(*trg>>8);
        st[2] = 0xFF&*trg;
        st[3] = 0;
        PRINTFB(G,FB_Match,FB_Warnings)
          " Match-Warning: unknown residue type '%s' (using X).\n",
          st
          ENDFB(G);
      }
      *trg='X';
    }
  }
  return(ok);
}

int MatchPreScore(CMatch *I,int *vla1,int n1,int *vla2,int n2,int quiet)
{
  PyMOLGlobals *G=I->G;
  int a,b;
  if(!quiet) {
    PRINTFB(G,FB_Match,FB_Details)
      " Match: assigning %d x %d pairwise scores.\n",n1,n2
      ENDFB(G);
  }

  for(a=0;a<n1;a++) {
    for(b=0;b<n2;b++) {
      I->mat[a][b] = I->smat[0x7F&vla1[a*3+2]][0x7F&vla2[b*3+2]];
    }
  }
  return 1;
}


#define BLOSUM62_ROWS 33
#define BLOSUM62_COLS 80

static char blosum62[BLOSUM62_ROWS][BLOSUM62_COLS] = {
"#  Matrix made by matblas from blosum62.iij\n",
"#  * column uses minimum score\n",
"#  BLOSUM Clustered Scoring Matrix in 1/2 Bit Units\n",
"#  Blocks Database = /data/blocks_5.0/blocks.dat\n",
"#  Cluster Percentage: >= 62\n",
"#  Entropy =   0.6979, Expected =  -0.5209\n",
"   A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V  B  Z  X  *\n",
"A  4 -1 -2 -2  0 -1 -1  0 -2 -1 -1 -1 -1 -2 -1  1  0 -3 -2  0 -2 -1  0 -4\n",
"R -1  5  0 -2 -3  1  0 -2  0 -3 -2  2 -1 -3 -2 -1 -1 -3 -2 -3 -1  0 -1 -4\n", 
"N -2  0  6  1 -3  0  0  0  1 -3 -3  0 -2 -3 -2  1  0 -4 -2 -3  3  0 -1 -4\n", 
"D -2 -2  1  6 -3  0  2 -1 -1 -3 -4 -1 -3 -3 -1  0 -1 -4 -3 -3  4  1 -1 -4\n", 
"C  0 -3 -3 -3  9 -3 -4 -3 -3 -1 -1 -3 -1 -2 -3 -1 -1 -2 -2 -1 -3 -3 -2 -4\n", 
"Q -1  1  0  0 -3  5  2 -2  0 -3 -2  1  0 -3 -1  0 -1 -2 -1 -2  0  3 -1 -4\n", 
"E -1  0  0  2 -4  2  5 -2  0 -3 -3  1 -2 -3 -1  0 -1 -3 -2 -2  1  4 -1 -4\n", 
"G  0 -2  0 -1 -3 -2 -2  6 -2 -4 -4 -2 -3 -3 -2  0 -2 -2 -3 -3 -1 -2 -1 -4\n", 
"H -2  0  1 -1 -3  0  0 -2  8 -3 -3 -1 -2 -1 -2 -1 -2 -2  2 -3  0  0 -1 -4\n", 
"I -1 -3 -3 -3 -1 -3 -3 -4 -3  4  2 -3  1  0 -3 -2 -1 -3 -1  3 -3 -3 -1 -4\n", 
"L -1 -2 -3 -4 -1 -2 -3 -4 -3  2  4 -2  2  0 -3 -2 -1 -2 -1  1 -4 -3 -1 -4\n", 
"K -1  2  0 -1 -3  1  1 -2 -1 -3 -2  5 -1 -3 -1  0 -1 -3 -2 -2  0  1 -1 -4\n", 
"M -1 -1 -2 -3 -1  0 -2 -3 -2  1  2 -1  5  0 -2 -1 -1 -1 -1  1 -3 -1 -1 -4\n", 
"F -2 -3 -3 -3 -2 -3 -3 -3 -1  0  0 -3  0  6 -4 -2 -2  1  3 -1 -3 -3 -1 -4\n", 
"P -1 -2 -2 -1 -3 -1 -1 -2 -2 -3 -3 -1 -2 -4  7 -1 -1 -4 -3 -2 -2 -1 -2 -4\n", 
"S  1 -1  1  0 -1  0  0  0 -1 -2 -2  0 -1 -2 -1  4  1 -3 -2 -2  0  0  0 -4\n", 
"T  0 -1  0 -1 -1 -1 -1 -2 -2 -1 -1 -1 -1 -2 -1  1  5 -2 -2  0 -1 -1  0 -4\n", 
"W -3 -3 -4 -4 -2 -2 -3 -2 -2 -3 -2 -3 -1  1 -4 -3 -2 11  2 -3 -4 -3 -2 -4\n", 
"Y -2 -2 -2 -3 -2 -1 -2 -3  2 -1 -1 -2 -1  3 -3 -2 -2  2  7 -1 -3 -2 -1 -4\n", 
"V  0 -3 -3 -3 -1 -2 -2 -3 -3  3  1 -2  1 -1 -2 -2  0 -3 -1  4 -3 -2 -1 -4\n", 
"B -2 -1  3  4 -3  0  1 -1  0 -3 -4  0 -3 -3 -2  0 -1 -4 -3 -3  4  1 -1 -4\n", 
"Z -1  0  0  1 -3  3  4 -2  0 -3 -3  1 -1 -3 -1  0 -1 -3 -2 -2  1  4 -1 -4\n",
"X  0 -1 -1 -1 -2 -1 -1 -1 -1 -1 -1 -1 -1 -1 -2  0  0 -2 -1 -1 -1 -1 -1 -4\n",
"* -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4  1\n",
""};

int MatchMatrixFromFile(CMatch *I,char *fname,int quiet)
{
  PyMOLGlobals *G=I->G;

  int ok=1;
  FILE *f;
  char *buffer = NULL;
  char *p;
  char cc[255];
  char *code=NULL;
  unsigned int x,y;
  int a;
  int n_entry;
  unsigned int size;

  if(fname && fname[0] 
#ifdef _PYMOL_NOPY
     /* if Python is absent, then use the hardcoded copy of BLOSUM62 */
     && !(strcmp(fname,"BLOSUM62")==0)
#endif
     ) {
    f=fopen(fname,"rb");
    if(!f) {
      PRINTFB(G,FB_Match,FB_Errors) 
        " Match-Error: unable to open matrix file '%s'.\n",fname
        ENDFB(G);
      ok=false;
    } else {
      fseek(f,0,SEEK_END);
      size=ftell(f);
      fseek(f,0,SEEK_SET);
      
      buffer=(char*)mmalloc(size+255);
      ErrChkPtr(G,buffer);
      p=buffer;
      fseek(f,0,SEEK_SET);
      fread(p,size,1,f);
      p[size]=0;
      fclose(f);
    }
  } else {
    buffer = Alloc(char, BLOSUM62_ROWS * BLOSUM62_COLS);
    if(buffer) {
      p = buffer;
      a=0;
      while(blosum62[a][0]) {
        strcpy(p,&blosum62[a][0]);
        p+=strlen(p);
        a++;
      }
    } else {
      ok=false;
    }
  }

  if(ok&&buffer) {

    /* count codes */

    p=buffer;
    n_entry = 0;
    while(*p && ok) {
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
      code=(char*)Calloc(int,n_entry);
      
      /* read codes */
      
      p=buffer;
      n_entry = 0;
      while(*p && ok) {
        switch(*p) {
        case '#':
          break;
        default:
          if((*p)>32) {
            code[n_entry]=*p;
            n_entry++;
          }
          break;
        }
        p=ParseNextLine(p);
      }
      
      /* read values */

      p=buffer;
      while((*p) && ok) {
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
            }
          }
          break;
        }
        if(!ok) break;
        p=ParseNextLine(p);
      }
    }
    mfree(buffer);
  }
  if(ok) {
    if(!quiet) {
      PRINTFB(G,FB_Match,FB_Details)
        " Match: read scoring matrix.\n"
        ENDFB(G);
    }
  }
  FreeP(code);
  return(ok);
}

int MatchAlign(CMatch *I,float gap_penalty,float ext_penalty,
               int max_gap,int max_skip,int quiet,int window,
               float ante)
{
  PyMOLGlobals *G=I->G;
  register int a,b,f,g;
  register int nf,ng;
  register int sf,sg;
  register float **score;
  register float **da,**db;
  register int na=I->na,nb=I->nb;
  unsigned int dim[2];
  register int2 **point;
  float mxv;
  int mxa,mxb;
  float tst=0.0;
  register int gap=0;
  int *p;
  int cnt;
  int ok=true;
  const float MIN_SCORE = 0.0F;
  nf = na+1;
  ng = nb+1;
  da = I->da;
  db = I->db;
  if(!quiet) {
    PRINTFB(G,FB_Match,FB_Actions)
      " MatchAlign: aligning residues (%d vs %d)...\n",na,nb
      ENDFB(G);
  }

  dim[0]=nf;
  dim[1]=ng;
  VLAFreeP(I->pair);
  score = (float**)UtilArrayCalloc(dim,2,sizeof(float));
  point = (int2**)UtilArrayCalloc(dim,2,sizeof(int2));
  if(score&&point) {

    /* initialize the scoring matrix */
    for(f=0;f<nf;f++) {
      for(g=0;g<ng;g++) {
        score[f][g] = MIN_SCORE;
        point[f][g][0] = -1;
        point[f][g][1] = -1;
      }
    }
    /* now start walking backwards up the alignment */
    {
      int second_pass = false;
      for(b=nb-1;b>=0;b--) {
        for(a=na-1;a>=0;a--) {
          /* find the maximum scoring cell accessible from this position, 
           * while taking gap penalties into account */
         
          mxv = MIN_SCORE;
          mxa=-1;
          mxb=-1;

          /* search for asymmetric insertions and deletions */
          f = a+1;
          if((max_gap>=0)&&(second_pass)) {
            sf = a+2+max_gap;
            sg = b+2+max_gap;
            if(sg>ng) sg = ng;
            if(sf>nf) sf = nf;
          } else {
            sg = ng;
            sf = nf;
          }
          for(g=b+1;g<sg;g++) {
            tst = score[f][g];

            if(window) {
              register int aa=a,bb=b,ff=f,gg=g, cc;
              tst += ante;
              for(cc=0;cc<window;cc++) {
                if((ff>=0)&&(gg>=0)&&(ff<na)&&(gg<nb)) {
                  tst -= (float)fabs(da[a][ff]-db[b][gg]);
                  aa = ff;
                  bb = gg;
                  ff = point[aa][bb][0];
                  gg = point[aa][bb][1];
                } else 
                  break;
              }
            }
            
            if(!((f==na)||(g==nb))) {
              gap = g-(b+1);
              if(gap) tst+= gap_penalty + ext_penalty*(gap-1);
            }
            if(tst>mxv) {
              mxv = tst;
              mxa = f;
              mxb = g;
            }
          }
          g = b+1;

          for(f=a+1;f<sf;f++) {
            tst = score[f][g];

            if(window) { 
              register int aa=a,bb=b,ff=f,gg=g, cc;
              tst += ante;
              for(cc=0;cc<window;cc++) {
                if((ff>=0)&&(gg>=0)&&(ff<na)&&(gg<nb)) {
                  tst -= (float)fabs(da[a][ff]-db[b][gg]);
                  aa = ff;
                  bb = gg;
                  ff = point[aa][bb][0];
                  gg = point[aa][bb][1];
                } else 
                  break;
              }
            }


            if(!((f==na)||(g==nb))) {
              gap=(f-(a+1));
              if(gap) tst+= gap_penalty + ext_penalty*(gap-1);
            }
            if(tst>mxv) {
              mxv = tst;
              mxa = f;
              mxb = g;
            }
          }
         
          if(max_skip) {
            /* search for high scoring mismatched stretches */
           
            sf = a+1+max_skip;
            sg = b+1+max_skip;
            if(sf>nf) sf = nf;
            if(sg>ng) sg = ng;
           
            for(f=a+1;f<sf;f++) {
              for(g=b+1;g<sg;g++) {
                tst = score[f][g];

                if(window) { 
                  register int aa=a,bb=b,ff=f,gg=g, cc;
                  tst += ante;
                  for(cc=0;cc<window;cc++) {
                    if((ff>=0)&&(gg>=0)&&(ff<na)&&(gg<nb)) {
                      tst -= (float)fabs(da[a][ff]-db[b][gg]);
                      aa = ff;
                      bb = gg;
                      ff = point[aa][bb][0];
                      gg = point[aa][bb][1];
                    } else 
                      break;
                  }
                }

                /* only penalize if we are not at the end */
                if(!((f==na)||(g==nb))) {
                  gap = ((f-(a+1))+(g-(b+1)));
                  if(gap>1) tst+= 2*gap_penalty + ext_penalty*(gap-2);
                }
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
        
          second_pass = true;
        }
      }
    }

    if(Feedback(G,FB_Match,FB_Debugging)) {
      for(b=0;b<nb;b++) {
        for(a=0;a<na;a++) {
          printf("%4.1f(%2d,%2d)",score[a][b],point[a][b][0],point[a][b][1]);
        }
        printf("\n");
      }
    }

    /* find the best entry point */

    mxv = MIN_SCORE;
    mxa=0;
    mxb=0;
    for(b=0;b<nb;b++) {
      for(a=0;a<na;a++) {
        tst = score[a][b];
        if(tst>mxv) {
          mxv = tst;
          mxa = a;
          mxb = b;
        }
      }
    }
    I->pair = VLAlloc(int,2*(na>nb?na:nb));
    p=I->pair;
    a = mxa;
    b = mxb;
    cnt=0;
    while((a>=0)&&(b>=0)&&(a<na)&&(b<nb)) {
      *(p++)=a;
      *(p++)=b;
      f = point[a][b][0];
      g = point[a][b][1];
      a=f;
      b=g;
      cnt++;
    }
    PRINTFD(G,FB_Match)
      " MatchAlign-DEBUG: best entry %8.3f %d %d %d\n",mxv,mxa,mxb,cnt
      ENDFD;
    if(!quiet) {
      PRINTFB(G,FB_Match,FB_Results) 
        " MatchAlign: score %1.3f\n",mxv
        ENDFD;
    }
    I->score = mxv;
    I->n_pair = cnt;
    VLASize(I->pair,int,(p-I->pair));
    FreeP(score);
    FreeP(point);
  }
  return(ok);
}

#if 0
int MatchAlign(CMatch *I,float gap_penalty,float ext_penalty,
               int max_gap,int max_skip,int quiet)
{
  PyMOLGlobals *G=I->G;
  unsigned int dim[2];
  register int a,b,f,g;
  register int nf,ng;
  register int sf,sg;
  register float **score;
  register int na=I->na,nb=I->nb;
  register int2 **point;
  float mxv;
  int mxa,mxb;
  float tst=0.0;
  register int gap=0;
  int *p;
  int cnt;
  int ok=true;
  const float MIN_SCORE = 0.0F;
  nf = na+1;
  ng = nb+1;

  if(!quiet) {
    PRINTFB(G,FB_Match,FB_Actions)
      " MatchAlign: aligning residues (%d vs %d)...\n",na,nb
      ENDFB(G);
  }

  dim[0]=nf;
  dim[1]=ng;
  VLAFreeP(I->pair);
  score = (float**)UtilArrayCalloc(dim,2,sizeof(float));
  point = (int2**)UtilArrayCalloc(dim,2,sizeof(int2));
  if(score&&point) {

    /* initialize the scoring matrix */
    for(f=0;f<nf;f++) {
      for(g=0;g<ng;g++) {
        score[f][g] = MIN_SCORE;
        point[f][g][0] = -1;
        point[f][g][1] = -1;
      }
    }
    /* now start walking backwards up the alignment */
    
    {
      int second_pass = false;
      for(b=nb-1;b>=0;b--) {
        for(a=na-1;a>=0;a--) {
         
          /* find the maximum scoring cell accessible from this position, 
           * while taking gap penalties into account */
         
          mxv = MIN_SCORE;
          mxa=-1;
          mxb=-1;
         
          /* search for asymmetric insertions and deletions */
          f = a+1;
          if((max_gap>=0)&&(second_pass)) {
            sf = a+2+max_gap;
            sg = b+2+max_gap;
            if(sg>ng) sg = ng;
            if(sf>nf) sf = nf;
          } else {
            sg = ng;
            sf = nf;
          }
          for(g=b+1;g<sg;g++) {
            tst = score[f][g];
            if(!((f==na)||(g==nb))) {
              gap = g-(b+1);
              if(gap) tst+=gap_penalty+ext_penalty*(gap-1);
            }
            if(tst>mxv) {
              mxv = tst;
              mxa = f;
              mxb = g;
            }
          }
          g = b+1;

          for(f=a+1;f<sf;f++) {
            tst = score[f][g];
            if(!((f==na)||(g==nb))) {
              gap=(f-(a+1));
              if(gap) tst+=gap_penalty+ext_penalty*(gap-1);
            }
            if(tst>mxv) {
              mxv = tst;
              mxa = f;
              mxb = g;
            }
          }
         
          if(max_skip) {
            /* search for high scoring mismatched stretches */
           
            sf = a+1+max_skip;
            sg = b+1+max_skip;
            if(sf>nf) sf = nf;
            if(sg>ng) sg = ng;
           
            for(f=a+1;f<sf;f++) {
              for(g=b+1;g<sg;g++) {
                tst = score[f][g];
                /* only penalize if we are not at the end */
                if(!((f==na)||(g==nb))) {
                  gap = ((f-(a+1))+(g-(b+1)));
                  if(gap>1) tst+=2*gap_penalty+ext_penalty*(gap-2);
                }
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
        
          second_pass = true;
        }
      }
    }

    if(Feedback(G,FB_Match,FB_Debugging)) {
      for(b=0;b<nb;b++) {
        for(a=0;a<na;a++) {
          printf("%4.1f(%2d,%2d)",score[a][b],point[a][b][0],point[a][b][1]);
        }
        printf("\n");
      }
    }

    /* find the best entry point */

    mxv = MIN_SCORE;
    mxa=0;
    mxb=0;
    for(b=0;b<nb;b++) {
      for(a=0;a<na;a++) {
        tst = score[a][b];
        if(tst>mxv) {
          mxv = tst;
          mxa = a;
          mxb = b;
        }
      }
    }
    I->pair = VLAlloc(int,2*(na>nb?na:nb));
    p=I->pair;
    a = mxa;
    b = mxb;
    cnt=0;
    while((a>=0)&&(b>=0)&&(a<na)&&(b<nb)) {
      *(p++)=a;
      *(p++)=b;
      f = point[a][b][0];
      g = point[a][b][1];
      a=f;
      b=g;
      cnt++;
    }
    PRINTFD(G,FB_Match)
      " MatchAlign-DEBUG: best entry %8.3f %d %d %d\n",mxv,mxa,mxb,cnt
      ENDFD;
    if(!quiet) {
      PRINTFB(G,FB_Match,FB_Results) 
        " MatchAlign: score %1.3f\n",mxv
        ENDFD;
    }
    I->score = mxv;
    I->n_pair = cnt;
    VLASize(I->pair,int,(p-I->pair));
    FreeP(score);
    FreeP(point);
  }
  return(ok);
}
#endif

void MatchFree(CMatch *I)
{
  FreeP(I->da);
  FreeP(I->db);
  FreeP(I->mat);
  FreeP(I->smat);
  VLAFreeP(I->pair);
  OOFreeP(I);
}


