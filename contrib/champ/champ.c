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
#include"os_python.h"

#include"os_std.h"
#include"os_memory.h"
#include"const.h"
#include"err2.h"
#include"sort.h"

#include"feedback2.h"

#include"mac.h"
#include"vla.h"
#include"champ.h"
#include"strblock.h"
#include"chiral.h"

#define _MATCHDEBUG
#define _MATCHDEBUG2
#define _RINGDEBUG
#define _AROMDEBUG

char *ChampParseAliphaticAtom(CChamp *I,char *c,int atom,int mask,int len,int imp_hyd);
char *ChampParseAromaticAtom(CChamp *I,char *c,int atom,int mask,int len,int imp_hyd);
char *ChampParseStringAtom(CChamp *I,char *c,int atom,int len);
int ChampParseAtomBlock(CChamp *I,char **c_ptr,int cur_atom);

void ChampAtomDump(CChamp *I,int index);
void ChampAtomFlagDump(CChamp *I,int index);

int ChampAddBondToAtom(CChamp *I,int atom_index,int bond_index);

int ChampMatch(CChamp *I,int template,int target,
                int unique_start,int n_wanted,int *match_start,int tag_mode);

int ChampMatch2(CChamp *I,int template,int target,
                        int start_tmpl,int start_targ,
                        int n_wanted,int *match_start,int tag_mode);

int ChampFindUniqueStart(CChamp *I,int template,int target,int *multiplicity);
int ChampUniqueListNew(CChamp *I,int atom, int unique_list);
void ChampUniqueListFree(CChamp *I,int unique_list);

void ChampMatchDump(CChamp *I,int match_idx);
void ChampMatchFree(CChamp *I,int match_idx);
void ChampMatchFreeChain(CChamp *I,int match_idx);
void ChampCountRings(CChamp *I,int index);
void ChampPrepareTarget(CChamp *I,int index);
void ChampPreparePattern(CChamp *I,int index);
char *ChampParseBlockAtom(CChamp *I,char *c,int atom,int mask,int len,int not_flag);
void ChampCountBondsEtc(CChamp *I,int index);
void ChampCheckCharge(CChamp *I,int index);
char *ChampParseTag(CChamp *I,char *c,unsigned int *map,unsigned int *not_map,int *ok);


static int num_to_ring[12] = { 
  0,
  0,
  0,
  cH_Ring3,
  cH_Ring4,
  cH_Ring5,
  cH_Ring6,

  cH_Ring7,
  cH_Ring8,
  0,
  0,
  0,};

static int num_to_valence[9] = {
  cH_0Valence,
  cH_1Valence,
  cH_2Valence,
  cH_3Valence,
  cH_4Valence,
  cH_5Valence,
  cH_6Valence,
  cH_7Valence,
  cH_8Valence };

static int num_to_degree[9] = {
  cH_0Bond,
  cH_1Bond,
  cH_2Bond,
  cH_3Bond,
  cH_4Bond,
  cH_5Bond,
  cH_6Bond,
  cH_7Bond,
  cH_8Bond
};


#ifdef _HAPPY_UT

static void diff(char *p,char *q)
{
  int same=true;
  char *s1,*s2;
  s1=p;
  s2=q;
  while((*p)&&(*q)) {
    if(*p!=*q)
      if(!((((*p)>='0')&&((*p)<='9'))&&
           (((*q)>='0')&&((*q)<='9'))))
        {
          same=false;
          break;
        }
    q++;
    p++;
  }
  if(!same)
    printf("%s\n%s\n",s1,s2);
}

int main(int argc, char *argv[])
{
  CChamp *h;
  int p;
  int pp = 0;
  char *smi;
  char buffer[1024];
  FILE *f,*g;
  int a,b,c=0;
  int len;
  int pat_list=0;
  int idx_list=0;
  int n_mat;
  int match_start;
  FeedbackInit();
  
  switch(5) {

 case 5: /* N-way cross comparison */
   h = ChampNew();
   f = fopen("test/ref/test01.smi","r");
   while(!feof(f)) {
     buffer[0]=0;
     fgets(buffer,1024,f);
      len = strlen(buffer);
      if(len) {
        if(buffer[len-1]<' ') buffer[len-1]=0;
        p = ChampSmiToPat(h,buffer);
        if(p) {
          ChampPreparePattern(h,p);
          a = h->Pat[p].unique_atom;
          printf("%d %s %d %d\n",p,buffer,a,ListLen(h->Int3,a));
          h->Pat[p].link = pat_list;
          pat_list = p;
        } else {
          printf("error %s\n",buffer);
        }
        c++;
        if(c==5000)
          break;
      }
   }
   fclose(f);
   
   c = 0;
    /* we now have a set of patterns loaded... */
    pp = pat_list;
    while(pp) {
      if(!h->Pat[pp].link)
        break;
      p = pat_list;
      n_mat = 0;
      while(p) {
        if(!h->Pat[p].link)
          break;
        
        
        n_mat += ChampMatch(h,pp,p,
                             ChampFindUniqueStart(h,pp,p,NULL),1,NULL);
        p = h->Pat[p].link;
      }
      c++;
      smi = ChampPatToSmiVLA(h,pp,NULL);
      printf("%d %d %s\n",c,n_mat,smi);
      pp = h->Pat[pp].link;
      vla_free(smi);
    }
    
    ChampFree(h);
    break;

  case 4: /* 1-way cross comparison */

    h = ChampNew();
    f = fopen("test/ref/test01.smi","r");
    while(!feof(f)) {
      buffer[0]=0;
      fgets(buffer,1024,f);
      len = strlen(buffer);
      if(len) {
        if(buffer[len-1]<' ') buffer[len-1]=0;
        p = ChampSmiToPat(h,buffer);
        if(p&&(!pp)) {
          ChampPreparePattern(h,p);
          pp = p;
          printf("%s\n",buffer);
        } else if(p) {
          /* locate unique atoms */
          ChampPreparePattern(h,p);
          b = ChampFindUniqueStart(h,pp,p,NULL);
          match_start = 0;
          n_mat = ChampMatch(h,pp,p,b,100,&match_start);
          
          /*  ChampPatDump(I,template);
              ChampMatchDump(I,match_start);*/
          ChampMatchFreeChain(h,match_start);

          if(n_mat) {
            printf("%d %d %d %d %d %d %s\n",c,n_mat,pp,p,b,h->Int3[b].value[0],buffer);
          }
          ChampPatFree(h,p);

        } else {
          printf("error %s\n",buffer);
        }
        c++;
        if(0&&(c==1000))
          break;
      }
    }
    fclose(f);
    /* we now have a set of patterns loaded... */
    ChampFree(h);
    break;

 case 3: /* self comparison */

    h = ChampNew();
    f = fopen("test/ref/test01.smi","r");
    while(!feof(f)) {
      buffer[0]=0;
      fgets(buffer,1024,f);
      len = strlen(buffer);
      if(len) {
        if(buffer[len-1]<' ') buffer[len-1]=0;
        p = ChampSmiToPat(h,buffer);
        if(p) {
          /* locate unique atoms */
          ChampPreparePattern(h,p);
          b = ChampFindUniqueStart(h,p,p,NULL);
          match_start=0;
          n_mat = ChampMatch(h,p,p,b,100,&match_start);
          if(!n_mat) {
            printf("Error!");
            exit(0);
          } else {
            printf("%d %d %s\n",c,n_mat,buffer);
          }
          ChampMatchFreeChain(h,match_start);          
          ChampPatFree(h,p);

        } else {
          printf("error %s\n",buffer);
        }
        c++;
        if(0&&(c==10))
          break;
      }
    }
    fclose(f);
    /* we now have a set of patterns loaded... */
    ChampFree(h);
    break;


  case 2: /* unique atoms */

    h = ChampNew();
    f = fopen("test/ref/test01.smi","r");
    while(!feof(f)) {
      buffer[0]=0;
      fgets(buffer,1024,f);
      len = strlen(buffer);
      if(len) {
        if(buffer[len-1]<' ') buffer[len-1]=0;
        p = ChampSmiToPat(h,buffer);
        if(p) {
          ChampPrepareUnique(h,p);
          a = h->Pat[p].unique_atom;
          printf("%d %s %d %d\n",p,buffer,a,ListLen(h->Int3,a));
          h->Pat[p].link = pat_list;
          pat_list = p;
        } else {
          printf("error %s\n",buffer);
        }
        c++;
        if(c==10)
          break;
      }
    }
    fclose(f);
    /* we now have a set of patterns loaded... */
    p = pat_list;
    while(p) {
      if(!h->Pat[p].link)
        break;
      b = ChampFindUniqueStart(h,p,h->Pat[p].link,&c);
      printf("\n%d vs %d, unique %d (mult %d)\n",p,h->Pat[p].link,b,c);
      idx_list = h->Int3[b].value[2];
      while(idx_list) {
        ChampAtomToString(h,h->Int[idx_list].value,buffer);
        printf("atom %s\n",buffer);
        idx_list = h->Int[idx_list].link;
      }
      p = h->Pat[p].link;
    }
    ChampFree(h);
    break;

  case 1: /* read & write smiles */

    /*      FeedbackEnable(FB_smiles_parsing,FB_everything);
            FeedbackEnable(FB_smiles_creation,FB_everything);*/
    h = ChampNew();
    f = fopen("test/ref/test01.smi","r");
    g = fopen("test/cmp/test01.smi","w");
    while(!feof(f)) {
      buffer[0]=0;
      fgets(buffer,1024,f);
      len = strlen(buffer);
      if(len) {
        if(buffer[len-1]<' ')
          buffer[len-1]=0;
        p = ChampSmiToPat(h,buffer);
        if(p) {
          /*      ChampPatDump(h,p);*/
          smi = ChampPatToSmiVLA(h,p,NULL);
          fprintf(g,"%s\n",smi);
          diff(buffer,smi);
          vla_free(smi);
          ChampPatFree(h,p);
        } else {
          printf("error %s\n",buffer);
        }
        if(!(c&0xFFF))
          printf(" %d\n",c);
        c++;
      }
    }
    fclose(f);
    fclose(g);
    ChampFree(h);
    break;
  }
  FeedbackFree();
  os_memory_dump();
  return 0;
}

#endif

char *ChampParseTag(CChamp *I,char *c,unsigned int *map,unsigned int *not_map,int *ok)
{
  /* parse bit masks like <1> <1,2,3> <12,3,1> etc... */

  int map_mask;
  int map_index;
  int not_flag = false;

  while(*ok) {
    if((*c)=='>') {
      c++;
      break;
    }
    if(!c) {
      *ok=false;
      break;
    }
    if(*c==';') { 
      not_flag=false;
      c++;
    } else if(*c=='!') {
      not_flag=true;
      c++;
    } else if((*c>='0')&&(*c<='9')) {
      if((*(c+1)>='0')&&(*(c+1)<='9')) {
        map_index = (*c-'0')*10+(*(c+1)-'0');
        c+=2;
      } else {
        map_index = (*c-'0');
        c++;
      }
      map_mask = 0x1;
      while(map_index) {
        map_mask = (map_mask<<1);
        map_index--;
      }
      if(not_flag) {
        *not_map|=map_mask;
      } else {
        *map|=map_mask;
      }
    } else 
      c++;
  }
  return(c);
}

static void merge_lineages(CChamp *I,int *src,int *src_mask,int *dst,int *dst_mask)
{
  int i;
  int i_src;
  int i_dst;
  i_src = *src;
  i_dst = *dst;

  while(i_src) {
    i = I->Int[i_src].value;
    if(!dst_mask[i]) {
      dst_mask[i]=1;
      *dst = ListElemPushInt(&I->Int,*dst,i);
    }
    i_src = I->Int[i_src].link;
  }

  while(i_dst) {
    i = I->Int[i_dst].value;
    if(!src_mask[i]) {
      src_mask[i]=1;
      *src = ListElemPushInt(&I->Int,*src,i);
    }
    i_dst = I->Int[i_dst].link;
  }

}

static int combine_lineage(CChamp *I,int src,int dst,int *lin_mask)
{
  int i;
  while(src) {
    i = I->Int[src].value;
    if(!lin_mask[i]) {
      lin_mask[i]=1;
      dst = ListElemPushInt(&I->Int,dst,i);
    }
    src = I->Int[src].link;
  }
  return(dst);
}

static int unrelated_lineage(CChamp *I,int index,int *lin_mask)
{
  int result = true;
#ifdef RINGDEBUG
  printf("%d %d %d %d %d\n",
         lin_mask[0],
         lin_mask[1],
         lin_mask[2],
         lin_mask[3],
         lin_mask[4]);
#endif

  while(index) {
#ifdef RINGDEBUG
    printf(" lineage: %d %d\n",I->Int[index].value,lin_mask[I->Int[index].value]);
#endif

    if(lin_mask[I->Int[index].value]) {
      result=false;
      break;
    }
    index = I->Int[index].link;
  }
  return result;
}


void ChampCountRings(CChamp *I,int index)
{
  ListPat *pat;
  ListBond *bd,*bd0,*bd1,*bd2,*bd3;
  ListAtom *at,*at0,*at1,*at2,*at3,*at4;
  int a,a0,a1,ai,bi,ai0,ai1,ai2,a2,bd_a0,bd_a1,ai3,ai4;
  int i2,i0,i1,i3;
  int bi0,bi1,bi2,bi3;
  
  int ni1;
  int n_atom = 0;
  int ring_size;
  pat = I->Pat + index;
  
  ai = pat->atom;
  while(ai) { /* count number of atoms */
    n_atom++;
    at = I->Atom + ai;
    ai = at->link;
  }

  if(n_atom) {
    ListAtom **atom = NULL;
    int *atom_idx;
    int *bonds = NULL; 
    int *neighbors = NULL;
    int lin_mask_size;
    int **lin_mask_vla = NULL;
    int *lin_mask;
    int lin_index;
    int n_expand;
    /* create a packed array of the atom pointers -- these will be our
       atom indices for this routine */

    atom = mac_malloc_array(ListAtom*,n_atom);
    atom_idx = mac_malloc_array(int,n_atom);
    a = 0;
    ai = pat->atom;
    while(ai) { /* load array and provide back-reference in mark_targ */
      at = I->Atom + ai;

      at->cycle = 0; /* initialize */
      at->class = 0;

      atom[a] = at;
      atom_idx[a] = ai;
      at->mark_targ = a; /* used for back-link */
      at->mark_read = 0; /* used to mark exclusion */
      at->mark_tmpl = 0; /* used to mark as dirty */

      at->first_targ = 0; /* used for lineage list */
      at->first_tmpl = 0; /* user for lineage mask index */

      at->first_base = 0; /* used for storing branches */
      ai = at->link;
      a++;
    }
    
    /*  build an efficient neighbor/bond list for all atoms */

    bonds = os_calloc(n_atom,sizeof(int)); /* list-starts of global bonds indices */
    neighbors = os_calloc(n_atom,sizeof(int)); /* list-starts of local atom indices */

    bi = pat->bond;
    while(bi) {
      bd = I->Bond + bi;
      
      bd->cycle = 0; /* initialize */
      bd->class = 0;
      
      ai0 = bd->atom[0]; /* global index */
      at0 = I->Atom+ai0;
      a0 = at0->mark_targ; /* local */
      ai1 = bd->atom[1];
      at1 = I->Atom+ai1;
      a1 = at1->mark_targ;

      if((bd->order)&(cH_Double|cH_Triple)) { /* record bond & atoms with Pi bonds */
        at0->class |= cH_Pi;
        at1->class |= cH_Pi;
        bd->class |= cH_Pi;
      }

      bonds[a0] = ListElemPushInt(&I->Int,bonds[a0],bi); /* enter this bond */
      bonds[a1] = ListElemPushInt(&I->Int,bonds[a1],bi);
      neighbors[a0] = ListElemPushInt(&I->Int,neighbors[a0],a1); /* cross-enter neighbors */

#ifdef RINGDEBUG
      printf(" ring: neighbor of %d is %d (%d)\n",a0,a1,neighbors[a0]);
#endif

      neighbors[a1] = ListElemPushInt(&I->Int,neighbors[a1],a0);
#ifdef RINGDEBUG
      printf(" ring: neighbor of %d is %d (%d)\n",a1,a0,neighbors[a1]);
#endif
      
      bi = bd->link;
    }
   
    /* set up data structures we'll need for ring finding */

    lin_mask_size = sizeof(int)*n_atom;
    lin_mask_vla = (int**)VLAMalloc(100,lin_mask_size,5,1); /* auto-zero */

    /* okay, repeat the ring finding process for each atom in the molecule 
       (optimize later to exclude dead-ends) */
    
    for(a=0;a<n_atom;a++) {
      int expand;
      int clean_up;
      int sz;
      int next;

      lin_index = 1;

#ifdef RINGDEBUG
      printf("==========\n");
      printf(" ring: starting with atom %d\n",a);
      ChampMemoryDump(I);
#endif
      expand = 0;
      clean_up = 0;
      sz = 1;

      at = atom[a];

      clean_up = ListElemPushInt(&I->Int,clean_up,a); /* mark this atom as dirty */
      at->mark_tmpl = true;

      at->mark_read = true; /* exclude this atom */

      expand = ListElemPushInt(&I->Int,expand,a); /* store this atom for the expansion cycle */
      vla_check(lin_mask_vla,int*,lin_index); /* assign blank lineage for this atom */
      at->first_tmpl = lin_index;
      lin_index++;
      
      while(sz<(RING_SEARCH_CUTOFF+1)) { /* while the potential loop is small...*/

        next = 0; /* start list of atoms for growth pass */

        /* find odd cycles and note growth opportunies */

        while(expand) { /* iterate over each atom to expand */
          expand = ListElemPopInt(I->Int,expand,&a0);

#ifdef RINGDEBUG
          printf(" ring: expand sz %d a0 %d\n",sz,a0);
#endif

          ni1 = neighbors[a0]; 
          while(ni1) { /* iterate over each neighbor of that atom */

#ifdef RINGDEBUG
            printf(" ring: ni1 %d\n",ni1);
#endif
            ni1 = ListElemGetInt(I->Int,ni1,&a1);
#ifdef RINGDEBUG
            printf(" ring: neighbor of %d is %d\n",a0,a1);
            printf(" ring: neighbor %d mark_read %d\n",a1,atom[a1]->mark_read);
#endif

            if(!atom[a1]->mark_read) { /* if that neighbor hasn't yet been covered... */
              if(!atom[a1]->first_base) {
                next = ListElemPushInt(&I->Int,next,a1);  /* store it for growth pass */
              }
              atom[a1]->first_base = /* and record this entry to it */
                ListElemPushInt(&I->Int,atom[a1]->first_base,a0); 
              if(!atom[a1]->mark_tmpl) {/* mark for later clean_up */
                clean_up = ListElemPushInt(&I->Int,clean_up,a1);
                atom[a1]->mark_tmpl = true;
              }

            } else if(sz==atom[a1]->mark_read) { /* ...or if the neighbor has an equal exclusion level... */

#ifdef RINGDEBUG
              printf(" ring: checking lineage between %d %d\n",a0,a1);
#endif
              if(unrelated_lineage(I,atom[a0]->first_targ, /* unrelated? */
                                   (int*)(((char*)lin_mask_vla)+(atom[a1]->first_tmpl*lin_mask_size))))
                {
                  ring_size = sz*2-1;
#ifdef RINGDEBUG
                  printf(" ring: #### found odd cycle %d\n",ring_size);
#endif
                  /* then we have a cycle of size sz*2-1 (odd) */
                  at0 = atom[a0];
                  at1 = atom[a1];
                  merge_lineages(I,
                                 &at0->first_targ,
                                 (int*)(((char*)lin_mask_vla)+(at0->first_tmpl*lin_mask_size)),
                                 &at1->first_targ,
                                 (int*)(((char*)lin_mask_vla)+(at1->first_tmpl*lin_mask_size)));
                  /* set the atom bits appropriately */
                  at0->cycle|=num_to_ring[ring_size];
                  at1->cycle|=num_to_ring[ring_size];
                  i0 = bonds[a0];
                  ai0 = atom_idx[a0];
                  ai1 = atom_idx[a1];
                  while(i0) {
                    i0 = ListElemGetInt(I->Int,i0,&bi);
                    bd = I->Bond+bi;
                    bd_a0=bd->atom[0];
                    bd_a1=bd->atom[1];
                    if(((bd_a0==ai0)&&(bd_a1==ai1))||
                       ((bd_a1==ai0)&&(bd_a0==ai1)))
                      bd->cycle|=num_to_ring[ring_size];
                  }
                }
            }
          }
        }
        
        /* find even cycles and grow pass */

        n_expand = 0;
        while(next) { /* iterate over potential branches */
          next = ListElemPopInt(I->Int,next,&a1);

#ifdef RINGDEBUG
          printf(" ring: next %d\n",a1);
#endif

          at1 = atom[a1];

          /* see if the base atoms form distinct cycles */

          i0 = atom[a1]->first_base;
          
          while(i0) {
            a0 = I->Int[i0].value;
            i2 = I->Int[i0].link;
            while(i2) {
              a2 = I->Int[i2].value;
              if(unrelated_lineage(I,atom[a0]->first_targ,
                                   (int*)(((char*)lin_mask_vla)+(atom[a2]->first_tmpl*lin_mask_size))))
                {
                  ring_size = sz*2;
                  at0 = atom[a0];
                  at2 = atom[a2];
#ifdef RINGDEBUG
                  printf(" ring: #### found even cycle %d %d %d\n",ring_size,i0,i2);
#endif
                  at0->cycle|=num_to_ring[ring_size];
                  at1->cycle|=num_to_ring[ring_size];
                  at2->cycle|=num_to_ring[ring_size];

                  i1 = bonds[a1];
                  ai0 = atom_idx[a0];
                  ai1 = atom_idx[a1];
                  ai2 = atom_idx[a2];
                  while(i1) {
                    i1 = ListElemGetInt(I->Int,i1,&bi);
                    bd = I->Bond+bi;
                    bd_a0=bd->atom[0];
                    bd_a1=bd->atom[1];
                    if(((bd_a0==ai0)&&(bd_a1==ai1))||
                       ((bd_a1==ai0)&&(bd_a0==ai1))||
                       ((bd_a0==ai2)&&(bd_a1==ai1))||
                       ((bd_a1==ai2)&&(bd_a0==ai1)))
                      bd->cycle|=num_to_ring[ring_size];
                  }

                  /* we have a cycle of size sz*2 */
                  /* set atom and bond bits appropriately */
                }
              i2 = I->Int[i2].link;
            }
            i0 = I->Int[i0].link;
          }
            
          /* allocate a new lineage for this branch atom */
          vla_check(lin_mask_vla,int*,lin_index); 
          at1->first_tmpl = lin_index;
          lin_mask = (int*)(((char*)lin_mask_vla)+(lin_mask_size*at1->first_tmpl));
          lin_index++;

          /* extend this particular lineage with each linking atom's lineage*/

          i0 = atom[a1]->first_base;
          a0 = I->Int[i0].value;

          if(a0!=a) {
            lin_mask[a0] = 1; /* add base atom */
            at1->first_targ = ListElemPushInt(&I->Int,at1->first_targ,a0);
          }
          
          while(i0) { /* and lineage of base atom */
            a0 = I->Int[i0].value;
            at0 = atom[a0];
            at1->first_targ = combine_lineage(I,at0->first_targ,at1->first_targ,lin_mask);
            i0 = I->Int[i0].link;
          }
          
          expand = ListElemPushInt(&I->Int,expand,a1);
          n_expand++;
          at1->mark_read = sz+1;

          if(!at1->mark_tmpl) {/* mark for later clean_up */
            clean_up = ListElemPushInt(&I->Int,clean_up,a1);
            at1->mark_tmpl = 1;
          }

          /* clean up */

          ListElemFreeChain(I->Int,atom[a1]->first_base); /* free base atom list */
          atom[a1]->first_base = 0;
                              
          /* copy the lineage to the new atom */

        }
        sz++;
        if(n_expand<2) {/* on a uniconnected branch, so no point in continuing */
          break;
        }
      }

      /* clean-up expansion list */
      ListElemFreeChain(I->Int,expand);
      expand = 0;

      while(clean_up) {
        clean_up=ListElemPopInt(I->Int,clean_up,&a0);
        at = atom[a0];
        
        /* reset lineage list and mask */
        lin_mask = (int*)(((char*)lin_mask_vla)+ (lin_mask_size*at->first_tmpl));
        while(at->first_targ) {
          at->first_targ = ListElemPopInt(I->Int, at->first_targ,&a1);
          lin_mask[a1] = 0;
        }
#ifdef RINGDEBUG
        for(i0=0;i0<n_atom;i0++) {
          if(lin_mask[i0])
            printf("lineage mask unclean !\n");
        }
#endif

        at->mark_read = 0; /* used to mark exclusion */
        at->mark_tmpl = 0; /* used to mark as dirty */
        
        at->first_targ = 0; /* used for lineage list */
        at->first_tmpl = 0; /* user for lineage mask index */
        
        at->first_base = 0; /* used for storing branches */


      }
    }

#ifdef RINGDEBUG
    for(a=0;a<n_atom;a++) {
      i1 = bonds[a];
      while(i1) { /* bonds of atom 1 */
        i1 = ListElemGetInt(I->Int,i1,&bi1);
        bd1 = I->Bond + bi1;
        ai1 = bd1->atom[0];
        ai2 = bd1->atom[1];
        at1 = I->Atom + ai1;
        at2 = I->Atom + ai2;
        printf("%d %d %x\n",at1->mark_targ+1,at2->mark_targ+1,bd1->cycle);
      }
    }
#endif

    /* now determine which atoms/bonds are aromatic */

    /* the following substructures are considered "aromatic" if they occur in rings of 5 or 6
     * *=*-*=*
     * *=*-*-*=*
     * NOTE  this is not a chemically accurate definition of aromaticity, just one
     * that is easy to program, and which covers most common occurences -- an 80-90% soln.
     *
     * I'll figure out something better later on...
     */
        
    ai0 = pat->atom;
    while(ai0) { /* cycle through each atom */
      
      at0 = I->Atom + ai0;
      if(at0->cycle&(cH_Ring4|cH_Ring5|cH_Ring6)) {
        /* atom 0 is cyclic */
        i0 = bonds[at0->mark_targ];
        while(i0) { /* bonds of atom 0 */
          i0 = ListElemGetInt(I->Int,i0,&bi0);
          bd0 = I->Bond + bi0;
          if((bd0->order==cH_Double)&&(bd0->cycle&(cH_Ring4|cH_Ring5|cH_Ring6))) {
            /* bond 0 is double */

            ai1 = bd0->atom[0];
            if(ai0==ai1) ai1 = bd0->atom[1];
            at1 = I->Atom + ai1;

            if(at1->cycle&(cH_Ring4|cH_Ring5|cH_Ring6)) {
              /* atom 1 is cyclic */
              
              i1 = bonds[at1->mark_targ];
              while(i1) { /* bonds of atom 1 */
                i1 = ListElemGetInt(I->Int,i1,&bi1);
                bd1 = I->Bond + bi1;
                if((bd1->order==cH_Single)&&(bd1->cycle&(cH_Ring4|cH_Ring5|cH_Ring6))) {
                  /* bond 1 is single */

                  ai2 = bd1->atom[0];
                  if(ai1==ai2) ai2 = bd1->atom[1];
                  at2 = I->Atom + ai2;

                  if(at2->cycle&(cH_Ring4|cH_Ring5|cH_Ring6)) {
                    /* atom 2 is cyclic */

                    i2 = bonds[at2->mark_targ];
                    while(i2) {
                      i2 = ListElemGetInt(I->Int,i2,&bi2);                        
                      bd2 = I->Bond + bi2;
                      if((bd2->order==cH_Double)&&(bd2->cycle&(cH_Ring4|cH_Ring5|cH_Ring6))) {
                        /* bond 2 is double */
                        
                        ai3 = bd2->atom[0];
                        if(ai2==ai3) ai3 = bd2->atom[1];
                        at3 = I->Atom + ai3;
                        
                        if(at3->cycle&(cH_Ring4|cH_Ring5|cH_Ring6)) {
                          /* atom 3 is cyclic, therefore system is aromatic */
                          
                          at3->class|=cH_Aromatic|cH_Pi;
                          at2->class|=cH_Aromatic|cH_Pi;
                          at1->class|=cH_Aromatic|cH_Pi;
                          at0->class|=cH_Aromatic|cH_Pi;
                          bd2->class|=cH_Aromatic|cH_Pi;
                          bd1->class|=cH_Aromatic|cH_Pi;
                          bd0->class|=cH_Aromatic|cH_Pi;
                        }
                      }
                    }
                    
                    i2 = bonds[at2->mark_targ];
                    while(i2) {
                      i2 = ListElemGetInt(I->Int,i2,&bi2);                        
                      bd2 = I->Bond + bi2;
                      if((bd2->order==cH_Single)&&(bd2->cycle&(cH_Ring4|cH_Ring5|cH_Ring6))) {
                        /* bond 2 is single */
                        
                        ai3 = bd2->atom[0];
                        if(ai2==ai3) ai3 = bd2->atom[1];
                        
                        if(ai3!=ai1) {
                          /* avoid backtracking... */
                          
                          at3 = I->Atom + ai3;
                          
                          if(at3->cycle&(cH_Ring4|cH_Ring5|cH_Ring6)) {
                            /* atom 3 is cyclic */

                            i3 = bonds[at3->mark_targ];
                            while(i3) {
                              i3 = ListElemGetInt(I->Int,i3,&bi3);                        
                              bd3 = I->Bond + bi3;
                              if((bd3->order==cH_Double)&&(bd3->cycle&(cH_Ring4|cH_Ring5|cH_Ring6))) {
                                /* bond 3 is double */
                                
                                ai4 = bd3->atom[0];
                                if(ai3==ai4) ai4 = bd3->atom[1];
                                at4 = I->Atom + ai4;
                                
                                if(at4->cycle&(cH_Ring4|cH_Ring5|cH_Ring6)) {
                                  
                                  /* atom 4 is cyclic, therefore system is aromatic */
                                  
                                  at4->class|=cH_Aromatic|cH_Pi;                              
                                  at3->class|=cH_Aromatic|cH_Pi;
                                  at2->class|=cH_Aromatic|cH_Pi;
                                  at1->class|=cH_Aromatic|cH_Pi;
                                  at0->class|=cH_Aromatic|cH_Pi;
                                  bd3->class|=cH_Aromatic|cH_Pi;
                                  bd2->class|=cH_Aromatic|cH_Pi;
                                  bd1->class|=cH_Aromatic|cH_Pi;
                                  bd0->class|=cH_Aromatic|cH_Pi;
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
            }
          }
        }
      }
      ai0 = at0->link;
    }

    for(a=0;a<n_atom;a++) { /* default conclusions */
      at = atom[a];
      if(!at->cycle)
        at->cycle=cH_Acyclic;
      if(!(at->class&(cH_Aromatic|cH_Aliphatic)))
        at->class|=cH_Aliphatic;
    }

    bi = pat->bond;
    while(bi) {
      bd = I->Bond + bi;
      if(!bd->cycle)
        bd->cycle=cH_Acyclic;
      if(!(bd->class&(cH_Aromatic|cH_Aliphatic)))
        bd->class=cH_Aliphatic;
      bi = bd->link;
    }
    /* now, clean up our efficient bond and neighbor lists */

    for(a=0;a<n_atom;a++) {
      ListElemFreeChain(I->Int,neighbors[a]);
      ListElemFreeChain(I->Int,bonds[a]);
    }
    mac_free(neighbors);
    mac_free(atom);
    mac_free(atom_idx);
    mac_free(bonds);
      
    /* and free the other data structures */

    vla_free(lin_mask_vla);

  }
#ifdef RINGDEBUG
  ChampPatDump(I,index);
#endif
}

void ChampCountBondsEtc(CChamp *I,int index)
{

  ListPat *pat;
  ListAtom *at,*at0,*at1;
  int ai,ai0,ai1,bi;
  ListBond *bd;
  int val_adj;

  pat = I->Pat + index;

  /* initialize */

  ai = pat->atom;
  while(ai) { 
    at = I->Atom + ai;
    at->valence=0;
    at->degree=0;
    at->tot_hydro=0;
    ai = at->link;
  }

  /* count */
  bi = pat->bond;
  while(bi) {
    bd = I->Bond + bi;
    
    ai0 = bd->atom[0];
    ai1 = bd->atom[1];
    at0 = I->Atom + ai0;
    at1 = I->Atom + ai1;
    at0->degree++;
    at1->degree++;
    if(at0->atom&cH_H) at1->tot_hydro++;
    if(at1->atom&cH_H) at0->tot_hydro++;
    switch(bd->order) {
    case cH_Single: 
      at0->valence++; 
      at1->valence++; 
      break;
    case cH_Double:
      at0->valence+=2;
      at1->valence+=2;
      break;
    case cH_Triple: 
      at0->valence+=3; 
      at1->valence+=3; 
      break;
    }
    bi = bd->link;
  }
 
  /* convert to bit masks */
  ai = pat->atom;
  while(ai) { 
    at = I->Atom + ai;
    at->degree = num_to_degree[at->degree];

    if(at->comp_imp_hydro_flag) {
      at->imp_hydro = 0;
      switch(at->charge) { /* adjust effective valence w.r.t charge */
      case cH_Neutral: val_adj=0; break;
      default: val_adj=0; break;
      case cH_Cation: val_adj=-1; break;
      case cH_Dication: val_adj=-2; break;
      case cH_Anion: val_adj=1; break;
      case cH_Dianion: val_adj=2; break;
      }
      val_adj += at->valence;
      
      switch(at->atom) { /* now compute implicit hydrogens */
      case cH_B: 
        if(val_adj<2) 
          at->imp_hydro = 2 - val_adj;
        break;
      case cH_C: 
        if(val_adj<4)
          at->imp_hydro = 4 - val_adj;
        break;
      case cH_N: 
        if(val_adj<3)
          at->imp_hydro = 3 - val_adj;
        else if(val_adj<5)
          at->imp_hydro = 5 - val_adj;
        break;
      case cH_O: 
        if(val_adj<2)
          at->imp_hydro = 2 - val_adj;
        break;
      case cH_S: 
        if(val_adj<2)
          at->imp_hydro = 2 - val_adj;
        else if(val_adj<4)
          at->imp_hydro = 4 - val_adj;
        else if(val_adj<6)
          at->imp_hydro = 6 - val_adj;
        break;
      case cH_P:
        if(val_adj<3)
          at->imp_hydro = 3 - val_adj;
        else if(val_adj<5)
          at->imp_hydro = 5 - val_adj;
        break;
      case cH_F:
      case cH_Cl:
      case cH_Br:
      case cH_I:
          if(val_adj<1)
            at->imp_hydro = 1 - val_adj;
          break;
      }
      at->valence+=at->imp_hydro; /* compute total valence */
    }
    at->tot_hydro+=at->imp_hydro;
    at->hydro_flag = true; /* we have a valid hydrogen count... */
    at->valence=num_to_valence[at->valence];
    ai = at->link;
  }
}

void ChampCheckCharge(CChamp *I,int index)
{
  ListPat *pat;
  ListAtom *at;
  int ai;

  pat = I->Pat + index;

  ai = pat->atom;
  while(ai) { 
    at = I->Atom + ai;
    if(!at->charge) {
      at->charge=cH_Neutral;
    }
    ai = at->link;
  }
}

void ChampPrepareTarget(CChamp *I,int index)
{
  ListPat *pat;
  pat = I->Pat + index;  /* NOTE: assumes I->Pat is stationary! */
  
  if(!pat->target_prep) {
    pat->target_prep=1;

    ChampCountRings(I,index);
    ChampCountBondsEtc(I,index);
    ChampCheckCharge(I,index);
    if(pat->unique_atom) 
      ChampUniqueListFree(I,pat->unique_atom);
    pat->unique_atom = ChampUniqueListNew(I,pat->atom,0);
  }
}

void ChampPreparePattern(CChamp *I,int index)
{
  ListPat *pat;
  pat = I->Pat + index;  /* NOTE: assumes I->Pat is stationary! */
  if(!pat->unique_atom)  
    pat->unique_atom = ChampUniqueListNew(I,pat->atom,0);
}

static int PTruthCallStr(PyObject *object,char *method,char *argument)
{
  int result = false;
  PyObject *tmp;
  tmp = PyObject_CallMethod(object,method,"s",argument);
  if(tmp) {
    if(PyObject_IsTrue(tmp))
      result = 1;
    Py_DECREF(tmp);
  }
  return(result);
}

static int PConvPyObjectToInt(PyObject *object,int *value)
{
  int result = true;
  PyObject *tmp;
  if(!object)
    result=false;
  else   if(PyInt_Check(object)) {
    (*value) = (int)PyInt_AsLong(object);
  } else {
    tmp = PyNumber_Int(object);
    if(tmp) {
      (*value) = (int)PyInt_AsLong(tmp);
      Py_DECREF(tmp);
    } else 
      result=false;
  }
  return(result);
}

static void UtilCleanStr(char *s) /*remove flanking white and all unprintables*/{
  char *p,*q;
  p=s;
  q=s;
  while(*p)
	 if(*p>32)
		break;
	 else 
		p++;
  while(*p)
	 if(*p>=32)
		(*q++)=(*p++);
	 else
		p++;
  *q=0;
  while(q>=s)
	 {
		if(*q>32)
		  break;
		else
		  {
			(*q)=0;
			q--;
		  }
	 }
}

static int PConvPyObjectToStrMaxClean(PyObject *object,char *value,int ln)
{
  char *st;
  PyObject *tmp;
  int result=true;
  if(!object)
    result=false;
  else if(PyString_Check(object)) {
    st = PyString_AsString(object);
    strncpy(value,st,ln);
  } else {
    tmp = PyObject_Str(object);
    if(tmp) {
      st = PyString_AsString(tmp);
      strncpy(value,st,ln);
      Py_DECREF(tmp);
    } else
      result=0;
  }
  if(ln>0)
    value[ln]=0;
  else
    value[0]=0;
  UtilCleanStr(value);
  return(result);
}

/*static int PConvPyObjectToStrMaxLen(PyObject *object,char *value,int ln)
{
  char *st;
  PyObject *tmp;
  int result=true;
  if(!object)
    result=false;
  else if(PyString_Check(object)) {
    st = PyString_AsString(object);
    strncpy(value,st,ln);
  } else {
    tmp = PyObject_Str(object);
    if(tmp) {
      st = PyString_AsString(tmp);
      strncpy(value,st,ln);
      Py_DECREF(tmp);
    } else
      result=0;
  }
  if(ln>0)
    value[ln]=0;
  else
    value[0]=0;
  return(result);
}

static int PConvAttrToStrMaxLen(PyObject *obj,char *attr,char *str,int ll)
{
  int ok=true;
  PyObject *tmp;
  if(PyObject_HasAttrString(obj,attr)) {
    tmp = PyObject_GetAttrString(obj,attr);
    ok = PConvPyObjectToStrMaxLen(tmp,str,ll);
    Py_DECREF(tmp);
  } else {
    ok=false;
  }
  return(ok);
}
*/

static int PConvPyListToFloatArrayInPlace(PyObject *obj,float *ff,int ll)
{
  int ok = true;
  int a,l;
  if(!obj) { 
    ok=false;
  } else if(!PyList_Check(obj)) {
    ok=false;
  } else {
    l=PyList_Size(obj);
    if (l!=ll) 
      ok=false;
    else {
      if(!l)
        ok=-1;
      else
        ok=l;
      for(a=0;a<l;a++)
        *(ff++) = (float)PyFloat_AsDouble(PyList_GetItem(obj,a));
    }
    /* NOTE ASSUMPTION! */
  }
  return(ok);
}

int ChampModelToPat(CChamp *I,PyObject *model) 
{

  int nAtom,nBond;
  int a;
  int ok=true;
  int cur_atom=0,last_atom = 0;
  ListAtom *at;
  int cur_bond=0,last_bond = 0;
  ListBond *bd;
  int charge=0,order=1;
  int result = 0;
  int atom1,atom2;
  int *atom_index = NULL;
  int std_flag;
  char *c;

  PyObject *atomList = NULL;
  PyObject *bondList = NULL;
  PyObject *molec = NULL;
  PyObject *atom = NULL;
  PyObject *bnd = NULL;
  PyObject *index = NULL;
  PyObject *tmp = NULL;

  nAtom=0;
  nBond=0;

  atomList = PyObject_GetAttrString(model,"atom");
  if(atomList) 
    nAtom = PyList_Size(atomList);
  else 
    ok=err_message("ChampModel2Pat","can't get atom list");

  atom_index = mac_malloc_array(int,nAtom);

  if(ok) { 
	 for(a=nAtom-1;a>=0;a--) /* reverse order */
		{
        atom = PyList_GetItem(atomList,a);
        if(!atom) 
          ok=err_message("ChampModel2Pat","can't get atom");

        cur_atom = ListElemNewZero(&I->Atom);
        at = I->Atom + cur_atom;
        at->link = last_atom;
        last_atom = cur_atom;
        at->chempy_atom = atom;
        Py_INCREF(at->chempy_atom);

        atom_index[a] = cur_atom; /* for bonds */

        if(ok) {
          tmp = PyObject_GetAttrString(atom,"name");
          if (tmp)
            ok = PConvPyObjectToStrMaxClean(tmp,at->name,NAM_SIZE-1);
          if(!ok) 
            err_message("ChampModel2Pat","can't read name");
          Py_XDECREF(tmp);
        }

        if(ok) {
          if(PTruthCallStr(atom,"has","flags")) {         
            tmp = PyObject_GetAttrString(atom,"flags");
            if (tmp)
              ok = PConvPyObjectToInt(tmp,(int*)&at->tag);
            if(!ok) 
              err_message("ChampModel2Pat","can't read flags");
            Py_XDECREF(tmp);
          } else {
            at->tag = 0;
          }
        }

        if(ok) {
          if(PTruthCallStr(atom,"has","index")) {  /* note -- chempy models have 1-based index attributes
                                                      even though the arrays are zero-based */
            tmp = PyObject_GetAttrString(atom,"index");
            if (tmp)
              ok = PConvPyObjectToInt(tmp,(int*)&at->ext_index);
            if(!ok) 
              err_message("ChampModel2Pat","can't read index");
            Py_XDECREF(tmp);
          } else {
            at->index = 0;
          }
        }

        if(ok) {
          if(PTruthCallStr(atom,"has","coord")) {         
            tmp = PyObject_GetAttrString(atom,"coord");
            if (tmp)
              ok = PConvPyListToFloatArrayInPlace(tmp,at->coord,3);
            if(!ok) 
              err_message("ChampModel2Pat","can't read coordinates");
            Py_XDECREF(tmp);
          }
        }

        if(ok) {
          if(PTruthCallStr(atom,"has","formal_charge")) { 
            tmp = PyObject_GetAttrString(atom,"formal_charge");
            if (tmp)
              ok = PConvPyObjectToInt(tmp,&charge);
            if(!ok) 
              err_message("ChampModel2Pat","can't read formal_charge");
            Py_XDECREF(tmp);
          } else {
            charge = 0;
          }
          switch(charge) {
          case 0: at->charge = cH_Neutral; break;
          case 1: at->charge = cH_Cation; break;
          case 2: at->charge = cH_Dication; break;
          case -1: at->charge = cH_Anion; break;
          case -2: at->charge = cH_Dianion; break;
          }
        }

        if(ok) {
          tmp = PyObject_GetAttrString(atom,"resn");
          if (tmp)
            ok = PConvPyObjectToStrMaxClean(tmp,at->residue,RES_SIZE-1);
          if(!ok) 
            err_message("ChampModel2Pat","can't read resn");
          Py_XDECREF(tmp);
        }
        
		  if(ok) {
          tmp = PyObject_GetAttrString(atom,"symbol");
          if (tmp)
            ok = PConvPyObjectToStrMaxClean(tmp,at->symbol,SYM_SIZE-1);
          if(!ok) 
            err_message("ChampModel2Pat","can't read symbol");
          c = at->symbol;
          std_flag = false;
          switch(*c) {
          case 'C':
            switch(*(c+1)) {
            case 'l':
            case 'L':
              ChampParseAliphaticAtom(I,c,cur_atom,cH_Cl,2,false);
              std_flag = true;
              break;
            default:
              ChampParseAliphaticAtom(I,c,cur_atom,cH_C,1,true);
              std_flag = true;
              break;
            }
            break;
          case 'H': /* nonstandard */
            ChampParseAliphaticAtom(I,c,cur_atom,cH_H,1,false);
            std_flag = true;
            break;
          case 'N':
            ChampParseAliphaticAtom(I,c,cur_atom,cH_N,1,true);
            std_flag = true;
            break;      
          case 'O':
            ChampParseAliphaticAtom(I,c,cur_atom,cH_O,1,true);
            std_flag = true;
            break;      
          case 'B':
            switch(*(c+1)) {
            case 'r':
            case 'R':
              ChampParseAliphaticAtom(I,c,cur_atom,cH_Br,2,false);
              std_flag = true;
              break;
            default:
              ChampParseAliphaticAtom(I,c,cur_atom,cH_B,1,true);
              std_flag = true;
              break;
            }
            break;
          case 'A':
            ChampParseAliphaticAtom(I,c,cur_atom,cH_A,1,true);
            std_flag = true;
            break;      
          case 'P':
            ChampParseAliphaticAtom(I,c,cur_atom,cH_P,1,true);
            std_flag = true;
            break;      
          case 'S':
            ChampParseAliphaticAtom(I,c,cur_atom,cH_S,1,true);
            std_flag = true;
            break;      
          case 'F':
            ChampParseAliphaticAtom(I,c,cur_atom,cH_F,1,false);
            std_flag = true;
            break;      
          case 'I':
            ChampParseAliphaticAtom(I,c,cur_atom,cH_I,1,false);
            std_flag = true;
            break;      
          case 'E':
            ChampParseAliphaticAtom(I,c,cur_atom,cH_E,1,false);
            std_flag = true;
            break;      
          case 'G':
            ChampParseAliphaticAtom(I,c,cur_atom,cH_G,1,false);
            std_flag = true;
            break;      
          case 'J':
            ChampParseAliphaticAtom(I,c,cur_atom,cH_J,1,false);
            std_flag = true;
            break;      
          case 'L':
            ChampParseAliphaticAtom(I,c,cur_atom,cH_L,1,false);
            std_flag = true;
            break;      
          case 'M':
            ChampParseAliphaticAtom(I,c,cur_atom,cH_M,1,false);
            std_flag = true;
            break;      
          case 'Q':
            ChampParseAliphaticAtom(I,c,cur_atom,cH_Q,1,false);
            std_flag = true;
            break;      
          case 'R':
            ChampParseAliphaticAtom(I,c,cur_atom,cH_R,1,false);
            std_flag = true;
            break;      
          case 'T':
            ChampParseAliphaticAtom(I,c,cur_atom,cH_T,1,false);
            std_flag = true;
            break;      
          case 'X':
            ChampParseAliphaticAtom(I,c,cur_atom,cH_X,1,false);
            std_flag = true;
            break;      
          case 'Z':
            ChampParseAliphaticAtom(I,c,cur_atom,cH_Z,1,false);
            std_flag = true;
            break;      
          }
          if(!std_flag) {
            ChampParseStringAtom(I,c,cur_atom,strlen(c));
          }
          Py_XDECREF(tmp);
        }
        
		  if(!ok)
			 break;
		}
  }

  bondList = PyObject_GetAttrString(model,"bond");
  if(bondList) 
    nBond = PyList_Size(bondList);
  else
    ok=err_message("ChampModel2Pat","can't get bond list");

  if(ok) {
	 for(a=nBond-1;a>=0;a--) /* reverse order */
		{
        bnd = PyList_GetItem(bondList,a);
        if(!bnd) 
          ok=err_message("ChampModel2Pat","can't get bond");

        cur_bond = ListElemNewZero(&I->Bond);
        bd = I->Bond + cur_bond;
        bd->link = last_bond;
        last_bond = cur_bond;

        bd->chempy_bond = bnd;
        Py_INCREF(bd->chempy_bond);

        if(ok) {
          tmp = PyObject_GetAttrString(bnd,"order");
          if (tmp)
            ok = PConvPyObjectToInt(tmp,&order);
          if(!ok) 
            err_message("ChampModel2Pat","can't read bond order");
          switch(order) {
          case 1: bd->order = cH_Single; break;
          case 2: bd->order = cH_Double; break;
          case 3: bd->order = cH_Triple; break;
          case 4: bd->order = cH_Single; bd->class = cH_Aromatic|cH_Pi; break;
          }
          Py_XDECREF(tmp);
        }
        
        index = PyObject_GetAttrString(bnd,"index");
        if(!index) 
          ok=err_message("ChampModel2Pat","can't get bond indices");
        else {
          if(ok) ok = PConvPyObjectToInt(PyList_GetItem(index,0),&atom1);
          if(ok) ok = PConvPyObjectToInt(PyList_GetItem(index,1),&atom2);
          if(!ok) 
            err_message("ChampModel2Pat","can't read bond atoms");
          else {
            atom1 = atom_index[atom1];
            atom2 = atom_index[atom2];

            if(order==4) {
              I->Atom[atom1].class |= cH_Aromatic|cH_Pi;
              I->Atom[atom2].class |= cH_Aromatic|cH_Pi;
            }
            bd->atom[0]=atom1;
            bd->atom[1]=atom2;
            bd->pri[0]=0;
            bd->pri[1]=0;
            ChampAddBondToAtom(I,atom1,cur_bond);
            ChampAddBondToAtom(I,atom2,cur_bond);
          }
        }
        Py_XDECREF(index);
      }
  }
  Py_XDECREF(atomList);
  Py_XDECREF(bondList);

  if(PyObject_HasAttrString(model,"molecule")) {
    molec = PyObject_GetAttrString(model,"molecule"); /* returns new reference */
  } else {
    molec = NULL;
  }

  mac_free(atom_index);

  result = ListElemNewZero(&I->Pat); 
  if(result) 
    {
      I->ActivePatList = ListElemPushInt(&I->Int,I->ActivePatList,result);
      I->Pat[result].atom = cur_atom;
      I->Pat[result].bond = cur_bond;
      I->Pat[result].chempy_molecule = molec;
      ChampPatReindex(I,result);     
  }
  return(result);
}

/* =============================================================== 
 * Class-specific Memory Management 
 * =============================================================== */

CChamp *ChampNew(void) {
  CChamp *I;

  feedback_Init(); /* only has side-effects the first time */
  ChiralInit(); 

  I = mac_calloc(CChamp); 
  I->Atom = (ListAtom*)ListNew(8000,sizeof(ListAtom));
  I->Bond = (ListBond*)ListNew(32000,sizeof(ListBond));
  I->Int  = (ListInt*)ListNew(16000,sizeof(ListInt));
  I->Int2 = (ListInt2*)ListNew(16000,sizeof(ListInt2));
  I->Int3 = (ListInt3*)ListNew(16000,sizeof(ListInt3));
  I->Tmpl = (ListTmpl*)ListNew(1000,sizeof(ListTmpl));
  I->Targ = (ListTarg*)ListNew(1000,sizeof(ListTarg));
  I->Pat = (ListPat*)ListNew(1000,sizeof(ListPat));
  I->Scope = (ListScope*)ListNew(100,sizeof(ListScope));
  I->Str = StrBlockNew(1000);
  I->Match = (ListMatch*)ListNew(1000,sizeof(ListMatch));

  I->ActivePatList = 0;

  return I;
}

void ChampMemoryDump(CChamp *I) 
{
  printf(" ChampMemory: Pat    %d\n",ListGetNAlloc(I->Pat));
  printf(" ChampMemory: Atom   %d\n",ListGetNAlloc(I->Atom));
  printf(" ChampMemory: Bond   %d\n",ListGetNAlloc(I->Bond));
  printf(" ChampMemory: Int    %d\n",ListGetNAlloc(I->Int));
  printf(" ChampMemory: Int2   %d\n",ListGetNAlloc(I->Int2));
  printf(" ChampMemory: Int3   %d\n",ListGetNAlloc(I->Int3));
  printf(" ChampMemory: Tmpl   %d\n",ListGetNAlloc(I->Tmpl));
  printf(" ChampMemory: Targ   %d\n",ListGetNAlloc(I->Targ));
  printf(" ChampMemory: Scope  %d\n",ListGetNAlloc(I->Scope));
  printf(" ChampMemory: Match  %d\n",ListGetNAlloc(I->Match));

}

int ChampMemoryUsage(CChamp *I)
{
  return(
         ListGetNAlloc(I->Pat)+
         ListGetNAlloc(I->Atom)+
         ListGetNAlloc(I->Bond)+
         ListGetNAlloc(I->Int)+
         ListGetNAlloc(I->Int2)+
         ListGetNAlloc(I->Int3)+
         ListGetNAlloc(I->Tmpl)+
         ListGetNAlloc(I->Targ)+
         ListGetNAlloc(I->Scope)+
         ListGetNAlloc(I->Match));
}

void ChampFree(CChamp *I) {

  while(I->ActivePatList) {
    ChampPatFree(I,I->ActivePatList); /* will update ActivePatList */
  }

  ListFree(I->Pat);
  ListFree(I->Atom);
  ListFree(I->Bond);
  ListFree(I->Int);
  ListFree(I->Int2);
  ListFree(I->Int3);
  ListFree(I->Tmpl);
  ListFree(I->Targ);
  ListFree(I->Scope);
  ListFree(I->Match);

  StrBlockFree(I->Str);
  mac_free(I);
}

void ChampPatReindex(CChamp *I,int index)
{
  ListPat *pt;
  ListAtom *at;
  ListBond *bd;
  int ai,bi;
  int b_idx=0;
  int a_idx=0;
  if(index) {
    pt = I->Pat + index;
    ai = pt->atom;
    while(ai) {
      at = I->Atom + ai; 
      at->index = a_idx++;
      ai = at->link;
    }
    bi = pt->bond;
    while(bi) {
      bd = I->Bond + bi; 
      bd->index = b_idx++;
      bi = bd->link;
    }
  }
}

void ChampPatFree(CChamp *I,int index)
{
  ListPat *pt;
  if(index) {
    pt = I->Pat + index;
    ChampAtomFreeChain(I,I->Pat[index].atom);
    ChampBondFreeChain(I,I->Pat[index].bond);
    if(pt->chempy_molecule) {Py_DECREF(pt->chempy_molecule);}
    ChampUniqueListFree(I,I->Pat[index].unique_atom);
    ListElemFree(I->Pat,index);
    I->ActivePatList = ListElemPurgeInt(I->Int,I->ActivePatList,index);
  }
}

void ChampAtomFree(CChamp *I,int atom)
{
  ListAtom *at;
  if(atom) {
    at = I->Atom + atom;
    if(at->chempy_atom) {Py_DECREF(at->chempy_atom);}
  }
  ListElemFree(I->Atom,atom);
}

void ChampAtomFreeChain(CChamp *I,int atom)
{
  ListAtom *at;
  int i;
  i = atom;
  while(i) {
    at = I->Atom + i;
    if(at->chempy_atom) {Py_DECREF(at->chempy_atom);}
    i = I->Atom[i].link;
  }
  ListElemFreeChain(I->Atom,atom);
}

void ChampBondFree(CChamp *I,int bond)
{
  ListBond *bd;
  if(bond) {
    bd = I->Bond + bond;
    if(bd->chempy_bond) {Py_DECREF(bd->chempy_bond);}
  }
  ListElemFree(I->Bond,bond);
}

void ChampBondFreeChain(CChamp *I,int bond)
{
  ListBond *at;
  int i;
  i = bond;
  while(i) {
    at = I->Bond + i;
    if(at->chempy_bond) {Py_DECREF(at->chempy_bond);}
    i = I->Bond[i].link;
  }
  ListElemFreeChain(I->Bond,bond);
}

/* =============================================================== 
 * Unique atoms in patterns
 * =============================================================== */

/* Unique atoms are stored in Int3 as
   0: representative atom index 
   1: total count of identical atoms
   2: list index (I->Int) of all identical atoms 
*/

int ChampUniqueListNew(CChamp *I,int atom, int unique_list)
{
  /* returns new root for unique_list */

  int cur_atom,next_atom;
  int unique,unique_atom;
  int ident;

  cur_atom = atom;
  while(cur_atom) { /* iterate through each atom in the molecule */
    next_atom = I->Atom[cur_atom].link; 
    unique = unique_list;
    while(unique) { /* check to see if it matches an existing unique atom */
      unique_atom = I->Int3[unique].value[0];
      if(ChampPatIdentical(I->Atom+cur_atom,I->Atom+unique_atom)) {
        /* if so, then just add this atom the match list for this unique atom */
        I->Int3[unique].value[1]++; /* count */
        ident = ListElemNew(&I->Int); 
        I->Int[ident].link = I->Int3[unique].value[2];
        I->Int[ident].value = cur_atom;
        I->Int3[unique].value[2] = ident;
        cur_atom = 0; /* indicate a miss */
        break;
      }
      unique = I->Int3[unique].link;
    }
    if(cur_atom) { /* found a new unique atom */
      unique_list = ListElemPush(&I->Int3,unique_list);
      I->Int3[unique_list].value[0] = cur_atom;
      I->Int3[unique_list].value[1] = 1;
#if 0
      ChampAtomDump(I,cur_atom);
#endif
      ident = ListElemNew(&I->Int); 
      I->Int[ident].value = cur_atom;
      I->Int3[unique_list].value[2] = ident; /* init list of other identical atoms */
    }
    cur_atom = next_atom;

  }
#if 0
  printf("\n");
#endif
  return(unique_list);
}

void ChampUniqueListFree(CChamp *I,int unique_list)
{
  int unique = unique_list;
  while(unique) {
    ListElemFreeChain(I->Int,I->Int3[unique].value[2]);
    unique = I->Int3[unique].link;      
  }
  ListElemFreeChain(I->Int3,unique_list);
}

/* =============================================================== 
 * Comparison
 * =============================================================== */

int ChampMatch_1V1_B(CChamp *I,int pattern,int target)
{
  ChampPreparePattern(I,pattern);
  ChampPrepareTarget(I,target);
  return(ChampMatch(I,pattern,target,
                    ChampFindUniqueStart(I,pattern,target,NULL),
                    1,NULL,false));
}

int ChampMatch_1V1_N(CChamp *I,int pattern,int target,int limit,int tag_mode)
{
  ChampPreparePattern(I,pattern);
  ChampPrepareTarget(I,target);
  return(ChampMatch(I,pattern,target,
                    ChampFindUniqueStart(I,pattern,target,NULL),
                    limit,NULL,tag_mode));
}




int ChampMatch_1V1_Map(CChamp *I,int pattern,int target,int limit,int tag_mode)
{
  int match_start = 0;

  ChampPreparePattern(I,pattern);
  ChampPrepareTarget(I,target);
  ChampMatch(I,pattern,target,
             ChampFindUniqueStart(I,pattern,target,NULL),
             limit,&match_start,tag_mode);
  return(match_start);
}

int ChampMatch_1VN_N(CChamp *I,int pattern,int list)
{
  int target;
  int c = 0;
  ChampPreparePattern(I,pattern);
  while(list) {
    target = I->Int[list].value;
    ChampPrepareTarget(I,target);    
    if(ChampMatch(I,pattern,target,
                  ChampFindUniqueStart(I,pattern,target,NULL),
                  1,NULL,false))
      
      c++;
    list = I->Int[list].link;
  }
  return(c);
}

int ChampExact_1VN_N(CChamp *I,int pattern, int list)
{
  int target;
  int c = 0;

  /* NOTE: should call generalize on molecules first!!! */

  ChampPreparePattern(I,pattern);
  while(list) {
    target = I->Int[list].value;
    if(pattern==target)
      c++;
    else {
      ChampPrepareTarget(I,target);    
      if(ChampMatch(I,pattern,target,
                    ChampFindUniqueStart(I,pattern,target,NULL),
                    1,NULL,false)) {
        if(ChampMatch(I,target,pattern,
                      ChampFindUniqueStart(I,target,pattern,NULL),
                      1,NULL,false)) {
        c++;
        }
      }
    }
    list = I->Int[list].link;
  }
  return(c);
}


int ChampMatch_NV1_N(CChamp *I,int list,int target,int limit,int tag_mode)
{
  int pattern;
  int c = 0;
  ChampPrepareTarget(I,target);    
  while(list) {
    pattern = I->Int[list].value;
    ChampPreparePattern(I,pattern);
    if(ChampMatch(I,pattern,target,
                  ChampFindUniqueStart(I,pattern,target,NULL),
                  limit,NULL,tag_mode))
      c++;
    list = I->Int[list].link;
  }
  return(c);
}

int ChampFindUniqueStart(CChamp *I,int template, int target,int *multiplicity)
{ /* returns zero for no match */

  int unique_tmpl,unique_targ;
  int tmpl_atom,targ_atom;
  int best_unique = 0;
  int score;
  int best_score = 0;

  unique_tmpl = I->Pat[template].unique_atom;
  while(unique_tmpl) {
    tmpl_atom = I->Int3[unique_tmpl].value[0];
    unique_targ = I->Pat[target].unique_atom;
    score = 0;
    while(unique_targ) {
      targ_atom = I->Int3[unique_targ].value[0];

      if(ChampAtomMatch(I->Atom + tmpl_atom,I->Atom + targ_atom)) 
        score += I->Int3[unique_targ].value[1];
      unique_targ = I->Int3[unique_targ].link;
    }
#ifdef MATCHDEBUG
    if(!score) {
      printf("unable to match: ");
      ChampAtomDump(I,tmpl_atom);
    }
#endif
    if(!score) return 0; /* no matched atom */
    
    score = score * I->Int3[unique_tmpl].value[1]; /* calculate multiplicity */
    if((!best_score)||(score<best_score)) {
      best_unique = unique_tmpl;
      best_score = score;
    }
    unique_tmpl = I->Int3[unique_tmpl].link;
  }
  if(multiplicity) *multiplicity=best_score;
  return(best_unique);
}


void ChampMatchFreeChain(CChamp *I,int match_idx)
{
  int next_idx;
  while(match_idx) {
    next_idx = I->Match[match_idx].link;
    ChampMatchFree(I,match_idx);
    match_idx = next_idx;
  }
}

void ChampMatchFree(CChamp *I,int match_idx)
{
  if(match_idx) {
    ListElemFreeChain(I->Int2,I->Match[match_idx].atom);
    ListElemFreeChain(I->Int2,I->Match[match_idx].bond);
    ListElemFree(I->Match,match_idx);
  }
}

void ChampMatchDump(CChamp *I,int match_idx)
{
  int m_atom_idx,atom_idx;
  int m_bond_idx,bond_idx;
  if(match_idx) {
    m_atom_idx = I->Match[match_idx].atom;
    m_bond_idx = I->Match[match_idx].bond;
    while(m_atom_idx) {
      atom_idx = I->Int2[m_atom_idx].value[0];
      ChampAtomDump(I,atom_idx);
      printf("(%2d,%2d)-",atom_idx,I->Atom[atom_idx].index);
      atom_idx = I->Int2[m_atom_idx].value[1];
      ChampAtomDump(I,atom_idx);
      printf("(%2d,%2d)\n",atom_idx,I->Atom[atom_idx].index);
      m_atom_idx = I->Int2[m_atom_idx].link;
    }
    while(m_bond_idx) {
      bond_idx = I->Int2[m_bond_idx].value[0];
      printf("%2d:%2d(%2d)-",I->Bond[bond_idx].atom[0],
             I->Bond[bond_idx].atom[1],bond_idx);
      bond_idx = I->Int2[m_bond_idx].value[1];
      printf("%2d:%2d(%2d)\n",I->Bond[bond_idx].atom[0],
             I->Bond[bond_idx].atom[1],bond_idx);
      m_bond_idx = I->Int2[m_bond_idx].link;
    }
  }

}

int ChampMatch(CChamp *I,int template,int target,int unique_start,
               int n_wanted,int *match_start,int tag_mode) 
{ /* returns whether or not substructure exists, but doesn't do alignment */
  int n_match = 0;
  int start_targ;
  int tmpl_atom,targ_atom;
  int rep_targ_atom;
  int unique_targ;

#ifdef MATCHDEBUG
  printf("\n\n ChampMatch: temp %d targ %d uniq %d n_want %d start %p tag %d\n",
         template,target,unique_start,n_wanted,match_start,tag_mode);

#endif

  /* we'll only need to start the search from the represenatative atom for this type
     (it isn't necc. to iterate through the others, since we'll be trying all
     combinations within the target...
  */
  if(unique_start) {
    tmpl_atom = I->Int3[unique_start].value[0]; 
    unique_targ = I->Pat[target].unique_atom;
    while(unique_targ) { /* iterate through all unique types of target atoms */
      rep_targ_atom = I->Int3[unique_targ].value[0];
      if(ChampAtomMatch(I->Atom + tmpl_atom,I->Atom + rep_targ_atom)) 
        {
          start_targ = I->Int3[unique_targ].value[2];
          while(start_targ) { /* iterate through all target atoms of this type */
            targ_atom = I->Int[start_targ].value;
            /* we now have starting atoms for each structure */
            n_match += ChampMatch2(I,template,target,
                                    tmpl_atom,targ_atom,
                                    (n_wanted-n_match),
                                   match_start,tag_mode);
            start_targ = I->Int[start_targ].link;            
            if(n_match>=n_wanted) break;
          }
        }
      if(n_match>=n_wanted) break;
      unique_targ = I->Int3[unique_targ].link;
    }
  }
  return(n_match);
}

int ChampMatch2(CChamp *I,int template,int target,
                 int start_tmpl,int start_targ,int n_wanted,
                 int *match_start,int tag_mode)

{ /* does the template covalent tree match the target? */
  
  /* basic algorithm is a multilayer mark-and-sweep traversal 
     through all atoms and bonds starting from the 
     two nucleating atoms, which are assumed to be equivalent */
  int n_match = 0;
  int stereo_template = false;


#ifdef MATCHDEBUG
  printf("\n\n ChampMatch2: %d %d\n",start_tmpl,start_targ);
#endif

  {
    /* first, initialize the marks */
    ListAtom *at;
    ListBond *bd;
    int cur_atom,cur_bond;
   
    /* template uses mark_tmpl */

    cur_atom = I->Pat[template].atom;
    while(cur_atom) {
      at = I->Atom+cur_atom;
      at->mark_tmpl=0;
      if(at->stereo) stereo_template = true;
      cur_atom = at->link;
    }
    cur_bond = I->Pat[template].bond;
    while(cur_bond) {
      bd = I->Bond+cur_bond;
      bd->mark_tmpl=0;
      cur_bond = bd->link;
    }

    /* target uses mark_targ */

    cur_atom = I->Pat[target].atom;
    while(cur_atom) {
      at = I->Atom+cur_atom;
      at->mark_targ=0;
      cur_atom = at->link;
    }
    cur_bond = I->Pat[target].bond;
    while(cur_bond) {
      bd = I->Bond+cur_bond;
      bd->mark_targ=0;
      cur_bond = bd->link;
    }
  }

  /* okay, now start our sweep... */
  {
    int tmpl_stack=0;
    int targ_stack=0;
    int tmpl_par;
    int tmpl_idx,targ_idx;
    int bond_off;
    int bond_idx,atom_idx;
    int mode=0; /* 0: iterate template atom, 1: need target atom */
    int done_flag=false;
    int match_flag;
    int bond_start;
    int atom_start;
    int stereo_match;

    ListTmpl *tmpl_ent,*parent_ent;
    ListTarg *targ_ent;
    ListAtom *tmpl_at,*targ_at,*base_at;
    ListMatch *match_ent;
    ListInt2 *int2;

    /* stack entries contain:
       0: the respective bond & atom indices
       1: current bond index in atom 
       2: backward entry index (branch point)
    */
    

    /* initialize target stack */

    targ_stack = ListElemPush(&I->Targ,targ_stack);
    targ_ent = I->Targ + targ_stack;
    targ_ent->atom = start_targ;
    targ_ent->bond = 0; /* to get here */

    /* initalize template stack */

    tmpl_stack = ListElemPush(&I->Tmpl,tmpl_stack);
    tmpl_ent = I->Tmpl + tmpl_stack;
    tmpl_ent->atom = start_tmpl;
    tmpl_ent->parent = 0;
    tmpl_ent->bond = 0; /* bond index to get here */
    tmpl_ent->match = targ_stack; /* matching target index */
    tmpl_ent->targ_start = 0; /* for matches */

    /* initialize marks */
    I->Atom[start_tmpl].mark_tmpl=1;
    I->Atom[start_tmpl].first_tmpl = tmpl_stack;
    I->Atom[start_targ].mark_targ=1;
    I->Atom[start_targ].first_targ = targ_stack;

    while(!done_flag) {


      if(!targ_stack) break;

#ifdef MATCHDEBUG
      printf("============================\n");
      printf("tmpl: ");
      tmpl_ent = I->Tmpl + tmpl_stack;
      while(1) {
        ChampAtomDump(I,tmpl_ent->atom);
        printf("(%2d) ",tmpl_ent->atom);

        if(tmpl_ent->link)
          tmpl_ent=I->Tmpl + tmpl_ent->link;
        else
          break;
      }
      printf("\ntarg: ");
      targ_ent = I->Targ + targ_stack;
      while(1) {
        ChampAtomDump(I,targ_ent->atom);
        printf("(%2d) ",targ_ent->atom); 
        if(targ_ent->link)
          targ_ent=I->Targ + targ_ent->link;
        else
          break;
      }
      printf("\n");
#endif

      switch(mode) {
      case 0: /* iterate  template atom */
        /* first, see if there is an open bond somewhere in the tree... */
        
        tmpl_par = tmpl_stack; /* start from last entry */
        while(tmpl_par) {
          tmpl_ent = I->Tmpl + tmpl_par;
          tmpl_at = I->Atom + tmpl_ent->atom; /* get atom record */
          
          bond_off = 0;  
          bond_idx = tmpl_at->bond[bond_off];
          while(bond_idx) {  /* iterate over all bonds */
            if(I->Bond[bond_idx].mark_tmpl) { /* skip over marked bonds...*/
              bond_off++;
              bond_idx = tmpl_at->bond[bond_off];
            } else
              break; /* found an open bond */
          }

          if(bond_idx) { /* there is an open bond */

#ifdef MATCHDEBUG
            printf(" tmpl: bond_idx %d is open (%2d,%2d)\n",bond_idx,
                   I->Bond[bond_idx].atom[0],I->Bond[bond_idx].atom[1]);
#endif

            atom_idx = I->Bond[bond_idx].atom[1];
            if(atom_idx==tmpl_ent->atom)
              atom_idx = I->Bond[bond_idx].atom[0];
            
            I->Bond[bond_idx].mark_tmpl++; /* mark bond "in use" */


            tmpl_stack = ListElemPush(&I->Tmpl,tmpl_stack); /* allocate new record for atom */
            tmpl_ent = I->Tmpl + tmpl_stack;
            tmpl_ent->atom = atom_idx; /* record which atom */
            tmpl_ent->parent = tmpl_par; /* record parent identity */
            tmpl_ent->bond = bond_idx; /* record bond index */
            tmpl_ent->targ_start = 0; /* launch searches from this bond in parent */
            tmpl_ent->match = 0; /* no match yet */
            
            tmpl_at = I->Atom + atom_idx;
            tmpl_at->mark_tmpl++; /* mark atom "in use" */
            if(tmpl_at->mark_tmpl==1) { /* record where it was used */
              tmpl_at->first_tmpl=tmpl_stack;
            }
      
#ifdef MATCHDEBUG      
            printf(" tmpl: created tmpl_idx %d tmpl_par %d\n",
                   tmpl_stack,tmpl_par);
            printf(" targ: tmpl atom %d:",tmpl_ent->atom);
            ChampAtomDump(I,tmpl_ent->atom);
            printf("\n");
#endif

            mode = 1; /* now we need to find matching template bond/atom */
            break; /* ...leave loop to hunt */
          } else { 
#ifdef MATCHDEBUG
            printf(" tmpl: nothing found in level %d...\n",tmpl_par);
#endif

            tmpl_par = tmpl_ent->parent; /* otherwise, proceed up tree looking for 
                                      nearest open bond */
          }
        }
        if(!tmpl_par) { /* no open bonds-> complete match */

#ifdef MATCHDEBUG2
          printf(" tmpl: EXACT MATCH DETECTED\n");
          printf(" %d %d %d %d\n",template,target,start_tmpl,start_targ);
          printf("tmpl: ");
          tmpl_ent = I->Tmpl + tmpl_stack;
          while(1) {
            ChampAtomDump(I,tmpl_ent->atom);
            printf("(%2d) ",tmpl_ent->atom);
            
            if(tmpl_ent->link)
              tmpl_ent=I->Tmpl + tmpl_ent->link;
            else
              break;
          }
          printf("\ntarg: ");
          targ_ent = I->Targ + targ_stack;
          while(1) {
            ChampAtomDump(I,targ_ent->atom);
            printf("(%2d) ",targ_ent->atom); 
            if(targ_ent->link)
              targ_ent=I->Targ + targ_ent->link;
            else
              break;
          }
          printf("\n");
#endif
          stereo_match = true;
          if(stereo_template) { /* must check stereochemistry */

            int tmpl_idx2;
            int targ_idx2;
            int bond_idx2;
            int bond_off2;
            int n_pri;
            int tmpl_pri[MAX_BOND];
            int targ_pri[MAX_BOND];
            ListTmpl *tmpl_ent2;
            ListTarg *targ_ent2;
            ListBond *tmpl_bd,*targ_bd;
            ListAtom *tmpl_at1,*targ_at1;

            tmpl_idx = tmpl_stack; /* prepare to traverse... */
            while(tmpl_idx) {
              tmpl_ent = I->Tmpl + tmpl_idx;
              I->Atom[tmpl_ent->atom].mark_read = false;
              tmpl_idx = tmpl_ent->link;
            }
            
            tmpl_idx = tmpl_stack;
            targ_idx = targ_stack;
            while(tmpl_idx&&targ_idx&&stereo_match) {
              tmpl_ent = I->Tmpl + tmpl_idx;
              targ_ent = I->Targ + targ_idx;
              tmpl_at1 = I->Atom + tmpl_ent->atom;
              targ_at1 = I->Atom + targ_ent->atom;
              if(!tmpl_at1->mark_read) { /* only compare non-virtual atoms */
                tmpl_at1->mark_read = true;
                if(tmpl_at1->stereo) {
                  
                  if(!targ_at1->stereo) {
                    /*                    
                                          printf("failure 1 %d %d %d %d \n",tmpl_at1->stereo,targ_at1->stereo,
                                          tmpl_ent->atom,targ_ent->atom);*/
                    stereo_match = false;
                  }
                  /* this is a stereo-center -- so locate the corresponding bond priorities in 
                     template and the target */

                  n_pri = 0;

                  if(tmpl_at1->imp_hydro==1) { /* deal with implicit hydrogen case */
                    tmpl_pri[n_pri]=0; /* considered first, since atom priorities always > 0 */
                    if(targ_at1->imp_hydro==1)
                      targ_pri[n_pri]=0;
                    else if(targ_at1->tot_hydro==1) {
                      /* track down the wayward hydrogen priority */
                      bond_off2 = 0;
                      bond_idx2 = targ_at1->bond[bond_off2];
                      while(bond_idx2) {  /* iterate over all bonds */
                        targ_bd = I->Bond + bond_idx2;
                        if(I->Atom[targ_bd->atom[0]].atom&cH_H) {
                          targ_pri[n_pri] = targ_bd->pri[0];
                          break;
                        } else if(I->Atom[targ_bd->atom[1]].atom&cH_H) {
                          targ_pri[n_pri] = targ_bd->pri[1];
                          break;
                        }
                        bond_off2++;
                        bond_idx2 = targ_at1->bond[bond_off2];
                        if(!bond_idx2) 
                          err_message("ChampMatch2","can't locate explicit hydrogen on target!");
                      }
                      /* failure? */

                    } else {
                      printf("failure 2\n");
                      stereo_match = false;
                    }
                    n_pri++;
                  }
                  
                  if(!stereo_match) 
                    break;
                    
                  tmpl_idx2 = tmpl_stack; /* prepare to traverse bonds... */
                  while(tmpl_idx2) {
                    tmpl_ent2 = I->Tmpl + tmpl_idx2;
                    if(tmpl_ent2->bond)
                      I->Bond[tmpl_ent2->bond].mark_read = false;
                    tmpl_idx2 = tmpl_ent2->link;
                  }

                  tmpl_idx2 = tmpl_stack;
                  targ_idx2 = targ_stack;
                  while(tmpl_idx2&&targ_idx2) {
                    tmpl_ent2 = I->Tmpl + tmpl_idx2;
                    targ_ent2 = I->Targ + targ_idx2;
                    if(tmpl_ent2->bond&&targ_ent2->bond) {
                      tmpl_bd = I->Bond + tmpl_ent2->bond;
                      if(!tmpl_bd->mark_read) {
                        tmpl_bd->mark_read = true;
                        
                        targ_bd = I->Bond + targ_ent2->bond;

                        /* if matched bond touches this atom, then 
                           get bond priority for the atom */
                        
                        if(tmpl_bd->atom[0]==tmpl_ent->atom) { 

                          tmpl_pri[n_pri]=tmpl_bd->pri[0];
                          if(targ_bd->atom[0]==targ_ent->atom) {
                            targ_pri[n_pri]=targ_bd->pri[0];
                          } else { /* assuming targ_bd->atom[1]==targ_ent->atom */
                            targ_pri[n_pri]=targ_bd->pri[1];                            
                          }
                          if(n_pri<4) 
                            n_pri++;
                          else {
                            err_message("ChampMatch2","too many connections on chiral atom!");
                          }
                        } else if(tmpl_bd->atom[1]==tmpl_ent->atom) {
                          
                          tmpl_pri[n_pri]=tmpl_bd->pri[1];
                          if(targ_bd->atom[0]==targ_ent->atom) {
                            targ_pri[n_pri]=targ_bd->pri[0];
                          } else { /* assuming targ_bd->atom[1]==targ_ent->atom */
                            targ_pri[n_pri]=targ_bd->pri[1];                            
                          }
                          if(n_pri<4) 
                            n_pri++;
                          else {
                            err_message("ChampMatch2","too many connections on chiral atom!");
                          }
                        }
                      }
                    }
                    tmpl_idx2 = tmpl_ent2->link;
                    targ_idx2 = targ_ent2->link;
                  }
                  if(n_pri<4) {
                    err_message("ChampMatch2","stereo comparison on achiral atom.");
                    stereo_match = false;
                  } else {
                    /* now test for a handedness match */
                    int tmpl_hand,targ_hand;
                    tmpl_hand = ChiralHandedness(tmpl_pri);
                    targ_hand = ChiralHandedness(targ_pri);
                    if(tmpl_at1->stereo==targ_at1->stereo) { /* same clockwizedness */
                      if(tmpl_hand!=targ_hand) { /* opposite handedness */
                        stereo_match=false;
                      }
                    } else { /* opposite clockwizedness */
                      if(tmpl_hand==targ_hand) { /* same handedness */
                        stereo_match = false;
                      }
                    }
                  }
                }
              }
              tmpl_idx = tmpl_ent->link;
              targ_idx = targ_ent->link;
            }
          }
          if(stereo_match) {
            n_match++;
            
            if(n_wanted>0) 
              if(n_match>=n_wanted) done_flag=true;
            
            if(match_start) {
              (*match_start) = ListElemPush(&I->Match,*match_start);
              match_ent = I->Match + (*match_start);
              
              atom_start=0;
              bond_start=0;
              
              tmpl_idx = tmpl_stack; /* prepare to traverse... */
              while(tmpl_idx) {
                tmpl_ent = I->Tmpl + tmpl_idx;
                I->Atom[tmpl_ent->atom].mark_read = false;
                tmpl_idx = tmpl_ent->link;
              }
              tmpl_idx = tmpl_stack;
              targ_idx = targ_stack;
              while(tmpl_idx&&targ_idx) {
                tmpl_ent = I->Tmpl + tmpl_idx;
                targ_ent = I->Targ + targ_idx;
                if(!I->Atom[tmpl_ent->atom].mark_read) { /* save non-virtual atom match */
                  I->Atom[tmpl_ent->atom].mark_read = true;
                  atom_start = ListElemPush(&I->Int2,atom_start);
                  int2 = I->Int2 + atom_start;
                  int2->value[0] = tmpl_ent->atom;
                  int2->value[1] = targ_ent->atom;
                }
                
                if(tmpl_ent->bond) { /* record bond match */
                  bond_start = ListElemPush(&I->Int2,bond_start);
                  int2 = I->Int2 + bond_start;
                  int2->value[0] = tmpl_ent->bond;
                  int2->value[1] = targ_ent->bond;
                }
                
                tmpl_idx = tmpl_ent->link;
                targ_idx = targ_ent->link;
              }
              match_ent->atom = atom_start;
              match_ent->bond = bond_start;
              
#ifdef MATCHDEBUG2
              /*            ChampPatDump(I,template);
                            ChampPatDump(I,target);*/
              ChampMatchDump(I,*match_start);
#endif
            }
            
            if(tag_mode) { /* are we using tags to mark atoms and bonds? */
              
              tmpl_idx = tmpl_stack; /* prepare to read... */
              while(tmpl_idx) {
                tmpl_ent = I->Tmpl + tmpl_idx;
                I->Atom[tmpl_ent->atom].mark_read = false;
                tmpl_idx = tmpl_ent->link;
              }
              tmpl_idx = tmpl_stack;
              targ_idx = targ_stack;
              while(tmpl_idx&&targ_idx) {
                tmpl_ent = I->Tmpl + tmpl_idx;
                targ_ent = I->Targ + targ_idx;
                if(!I->Atom[tmpl_ent->atom].mark_read) { /* save non-virtual atom match */
                  I->Atom[tmpl_ent->atom].mark_read = true;
                  if(tag_mode==cTag_merge) { /* merge */                                                                                         
                    I->Atom[targ_ent->atom].tag |= I->Atom[tmpl_ent->atom].tag;
                    I->Atom[targ_ent->atom].tag &= (0xFFFFFFFF^I->Atom[tmpl_ent->atom].not_tag);
                  } else if(tag_mode==cTag_copy) { /* copy */
                    I->Atom[targ_ent->atom].tag = I->Atom[tmpl_ent->atom].tag;
                  }
                }
                
                if(tmpl_ent->bond) { /* record bond match */
                  if(tag_mode==cTag_merge) {
                    I->Bond[targ_ent->bond].tag |= I->Bond[tmpl_ent->bond].tag;
                    I->Bond[targ_ent->bond].tag &= (0xFFFFFFFF^I->Bond[tmpl_ent->bond].not_tag);
                  } else if(tag_mode==cTag_copy) {
                    I->Bond[targ_ent->bond].tag = I->Bond[tmpl_ent->bond].tag;
                  }
                }                
                tmpl_idx = tmpl_ent->link;
                targ_idx = targ_ent->link;
              }
            }
          }
          /* back-out the last target match... */
          
          targ_ent = I->Targ + targ_stack;

          atom_idx = targ_ent->atom;
          I->Atom[atom_idx].mark_targ--; /* free-up atom */
          bond_idx = targ_ent->bond;
          I->Bond[bond_idx].mark_targ--; /* free-up bond */
          targ_stack = ListElemPop(I->Targ,targ_stack);

          mode = 1; /* now we need to find another matching template bond/atom */
        }
        break;
      case 1:
        /* try to locate matching target atom & bond */

        tmpl_ent = I->Tmpl + tmpl_stack;
        bond_off = tmpl_ent->targ_start; /* start with this bond */
        parent_ent = I->Tmpl + tmpl_ent->parent; /* get parent of current template atom */
        targ_ent = I->Targ + parent_ent->match; /* start from target atom which
                                                   matches the parent in template */
        base_at = I->Atom + targ_ent->atom;

#ifdef MATCHDEBUG
        printf(" targ: tmpl_stack %d  bond_off %d\n",tmpl_stack,bond_off);
        printf(" targ: match %d targ root atom %d:",parent_ent->match,targ_ent->atom);
        ChampAtomDump(I,targ_ent->atom);
        printf("\n");
        printf(" targ: bonds: ");
        bond_start = bond_off;
        bond_idx = base_at->bond[bond_start];
        while(bond_idx) {  /* iterate over all bonds */
          bond_start++;
          printf("%d ",bond_idx);
          bond_idx = base_at->bond[bond_start];
        }
        printf("\n");
#endif

        bond_idx = base_at->bond[bond_off];
        while(bond_idx) {  /* iterate over all bonds */

#ifdef MATCHDEBUG
          printf(" targ: trying %d bond_idx %d (%2d-%2d)\n",bond_off,bond_idx,
                 I->Bond[bond_idx].atom[0],I->Bond[bond_idx].atom[1]);
#endif

          match_flag=false;
          tmpl_ent->targ_start = bond_off + 1; /* insure we never duplicate attempt */
          
          if(!I->Bond[bond_idx].mark_targ) { /* is bond unmarked...*/

            atom_idx = I->Bond[bond_idx].atom[1]; /* find other atom */
            if(atom_idx==targ_ent->atom)
              atom_idx = I->Bond[bond_idx].atom[0];

            tmpl_at = I->Atom+tmpl_ent->atom;
            targ_at = I->Atom+atom_idx;

            if((tmpl_at->mark_tmpl-targ_at->mark_targ)==1) {/* is atom available? */
              if(ChampBondMatch(I->Bond+tmpl_ent->bond,
                                I->Bond+bond_idx)) /* does bond match? */ {
                if(ChampAtomMatch(tmpl_at,targ_at)) /* does atom match? */
                  
                  {
                    match_flag=true;
                    if(targ_at->mark_targ==1) { /* virtual */
                      /* verify that matching virtual atoms correspond to matching real atoms */
#ifdef MATCHDEBUG
                      printf(" targ: first_tmpl %d first_targ %d\n",
                             tmpl_at->first_tmpl,targ_at->first_targ);
#endif
                      if(I->Tmpl[tmpl_at->first_tmpl].match!=targ_at->first_targ)
                        match_flag=false;
                    }
                    if(match_flag) {
                      I->Bond[bond_idx].mark_targ++; /* mark bond "in use" */
                      
                      targ_stack = ListElemPush(&I->Targ,targ_stack); /* allocate new record for atom */
                      targ_ent = I->Targ + targ_stack;
                      targ_ent->atom = atom_idx;
                      targ_ent->bond = bond_idx; 

                      targ_at->mark_targ++; /* mark atom "in use" */
                      if(targ_at->mark_targ==1) { /* record where it was used */
                        targ_at->first_targ = targ_stack;
                      }
                      /* inform template entry about match */
                      
                      tmpl_ent->match = targ_stack;
                      
                      match_flag=true;
                      mode = 0; /* we have a candidate atom/bond,
                                   so return to template */
                      
#ifdef MATCHDEBUG
                      printf(" targ: MATCHED atoms %d & %d (bond %d):",
                             tmpl_ent->atom,
                             targ_ent->atom,
                             bond_idx);
                      ChampAtomDump(I,tmpl_ent->atom);
                      ChampAtomDump(I,targ_ent->atom);
                      printf("\n");
#endif
                      
                      break;
                    }
                  } else {
#ifdef MATCHDEBUG
                printf(" atom match failed %d vs %d ",tmpl_ent->atom,atom_idx);
                ChampAtomDump(I,tmpl_ent->atom);
                ChampAtomDump(I,atom_idx);
                printf("\n");
#endif
                  }
              
              } else {
#ifdef MATCHDEBUG
                printf(" bond match failed\n");
#endif

              }
            }
          }
          if(!match_flag) {
            bond_off++;
            bond_idx = base_at->bond[bond_off];
          }
        }
        if(mode) { /* unable to locate a match */ 
          /* so back-off the previous template atom... */

          tmpl_ent = I->Tmpl + tmpl_stack;
          atom_idx = tmpl_ent->atom;
          I->Atom[atom_idx].mark_tmpl--; /* free-up atom */
          bond_idx = tmpl_ent->bond;
          I->Bond[bond_idx].mark_tmpl--; /* free-up bond */
          tmpl_stack =  ListElemPop(I->Tmpl,tmpl_stack);

          /* and back-off the previous target match */
          targ_ent = I->Targ + targ_stack;
          atom_idx = targ_ent->atom;
          I->Atom[atom_idx].mark_targ--; /* free-up atom */
          bond_idx = targ_ent->bond;
          I->Bond[bond_idx].mark_targ--; /* free-up bond */
          targ_stack = ListElemPop(I->Targ,targ_stack);

        }
      }
    }
    ListElemFreeChain(I->Tmpl,tmpl_stack);
    ListElemFreeChain(I->Targ,targ_stack);
  }
#ifdef MATCHDEBUG
  printf(" ChampMatch2: returning n_match = %d\n",
         n_match);
#endif

  return n_match;
}

/* =============================================================== 
 * Smiles
 * =============================================================== */

char *ChampParseBlockAtom(CChamp *I,char *c,int atom,int mask,int len,int not_flag)
{
  ListAtom *at;
  at=I->Atom+atom;
  if(not_flag) {
    at->not_atom |= mask;
    at->neg_flag = true;
  } else {
    at->atom |= mask;
    at->pos_flag = true;
  }
  at->hydro_flag=true;
  PRINTFD(FB_smiles_parsing) 
    " ChampParseBlockAtom: called.\n"
    ENDFD;
  if(mask==cH_Sym) {
    if(len==1) {
      at->symbol[0]=*c;
      at->symbol[1]=0;
    } else if(len==2) {
      at->symbol[0]=*c;
      at->symbol[1]=*(c+1);
      at->symbol[2]=0;
    }
  }
  /* need to include code for loading symbol */
  return c+len;
}

char *ChampParseAliphaticAtom(CChamp *I,char *c,int atom,int mask,int len,int imp_hyd) 
{
  ListAtom *at;
  at=I->Atom+atom;
  at->atom |= mask;
  at->pos_flag = true;
  at->comp_imp_hydro_flag = imp_hyd;
  PRINTFD(FB_smiles_parsing) 
    " ChampParseAliphaticAtom: called.\n"
    ENDFD;
  /* need to include code for loading symbol */
  return c+len;
}

char *ChampParseAromaticAtom(CChamp *I,char *c,int atom,int mask,int len,int imp_hyd) 
{
  ListAtom *at;
  at=I->Atom+atom;
  at->atom |= mask;
  at->class |= cH_Aromatic;
  at->pos_flag = true;
  at->comp_imp_hydro_flag = imp_hyd;
  PRINTFD(FB_smiles_parsing) 
    " ChampParseAromaticAtom: called.\n"
    ENDFD;
  return c+len;
}

char *ChampParseStringAtom(CChamp *I,char *c,int atom,int len) 
{
  ListAtom *at;
  at=I->Atom+atom;
  at->atom |= cH_Any;
  at->symbol[0]=c[0];
  at->symbol[1]=c[1];
  at->pos_flag = true;
  PRINTFD(FB_smiles_parsing) 
    " ChampParseStringAtom: called.\n"
    ENDFD;
  return c+len;
}

int ChampParseNumeral(char *c);

int ChampParseNumeral(char *c) {
  switch(*c) {
  case '0':
  case '1':
  case '2':
  case '3':
  case '4':
  case '5':
  case '6':
  case '7':
  case '8':
  case '9':
    return(*c-'0');
    break;
  default:
    return(-1);
  }
}

int ChampParseAtomBlock(CChamp *I,char **c_ptr,int cur_atom) 
{
  int ok = true;
  ListAtom *at;
  char *c;
  int not_flag = false;
  int num;
  int done=false;
  int atom_seen = false;

  /*  int done;*/

  c=*c_ptr;

  at=I->Atom+cur_atom;

  at->comp_imp_hydro_flag = false;
  
  while(ok&&!done) {
    switch(*c) {
    case ']':
      done=true;
      c++;
      break;
    case 0:
      done=true;
      break;
    case '!':
      c++;
      not_flag=true;
      atom_seen = false;
      break;
    case ',':
      c++;
      atom_seen = false;
      break;
    case ';':
      not_flag=false;
      c++;
      atom_seen = false;
      break;
    case '*': /* nonstandard */
      c = ChampParseBlockAtom(I,c,cur_atom,cH_Any,1,not_flag);
      atom_seen = true;
      break;
    case '?': /* nonstandard */
      c = ChampParseBlockAtom(I,c,cur_atom,cH_NotH,1,not_flag);
      atom_seen = true;
      break;
    case '@':
      num = ChampParseNumeral(c+1);
      if(num>=0) {
        c+=2;
      } else {
        num = 1;
        c++;
        while(*c) {
          if(*c!='@')
            break;
          else {
            num++;
            c++;
          }
        }
      }
      if(num&0x1) /* odd */
        I->Atom[cur_atom].stereo = cH_Anticlock;
      else
        I->Atom[cur_atom].stereo = cH_Clockwise;
      break;
    case 'A': /* note there is no way to address the 'A' symbol ...*/
      switch(*(c+1)) {
      case 'c':
        c = ChampParseBlockAtom(I,c,cur_atom,cH_Sym,2,not_flag);
        atom_seen = true;
        break;
      case 'g':
        c = ChampParseBlockAtom(I,c,cur_atom,cH_Sym,2,not_flag);
        atom_seen = true;
        break;
      case 'l':
        c = ChampParseBlockAtom(I,c,cur_atom,cH_Sym,2,not_flag);
        atom_seen = true;
        break;
      case 'm':
        c = ChampParseBlockAtom(I,c,cur_atom,cH_Sym,2,not_flag);
        atom_seen = true;
        break;
      case 's':
        c = ChampParseBlockAtom(I,c,cur_atom,cH_Sym,2,not_flag);
        atom_seen = true;
        break;
      case 'u':
        c = ChampParseBlockAtom(I,c,cur_atom,cH_Sym,2,not_flag);
        atom_seen = true;
        break;
      default:
        if(not_flag) {
          I->Atom[cur_atom].neg_flag=true;
          I->Atom[cur_atom].not_class|=cH_Aliphatic;      
        } else {
          I->Atom[cur_atom].pos_flag=true;
          I->Atom[cur_atom].class|=cH_Aliphatic;
        }
        c++;
      }
      break;
    case 'a':
      if(not_flag) {
        I->Atom[cur_atom].neg_flag=true;
        I->Atom[cur_atom].not_class|=cH_Aromatic;
      } else {
        I->Atom[cur_atom].class|=cH_Aromatic;
        I->Atom[cur_atom].pos_flag=true;
      }
      c++;    
      break;
    case 'B':
      switch(*(c+1)) {
      case 'a':
        c = ChampParseBlockAtom(I,c,cur_atom,cH_Sym,2,not_flag);
        atom_seen = true;
        break;
      case 'e':
        c = ChampParseBlockAtom(I,c,cur_atom,cH_Sym,2,not_flag);
        atom_seen = true;
        break;
      case 'i':
        c = ChampParseBlockAtom(I,c,cur_atom,cH_Sym,2,not_flag);
        atom_seen = true;
        break;
      case 'r':
        c = ChampParseBlockAtom(I,c,cur_atom,cH_Br,2,not_flag);
        atom_seen = true;
        break;
      default:
        c = ChampParseBlockAtom(I,c,cur_atom,cH_B,1,not_flag);
        atom_seen = true;
        break;
      }
      break;
    case 'C':
      switch(*(c+1)) {
      case 'a':
        c = ChampParseBlockAtom(I,c,cur_atom,cH_Ca,2,not_flag);
        atom_seen = true;
        break;
      case 'd':
        c = ChampParseBlockAtom(I,c,cur_atom,cH_Sym,2,not_flag);
        atom_seen = true;
        break;
      case 'e':
        c = ChampParseBlockAtom(I,c,cur_atom,cH_Sym,2,not_flag);
        atom_seen = true;
        break;
      case 'o':
        c = ChampParseBlockAtom(I,c,cur_atom,cH_Sym,2,not_flag);
        atom_seen = true;
        break;
      case 'r':
        c = ChampParseBlockAtom(I,c,cur_atom,cH_Sym,2,not_flag);
        atom_seen = true;
        break;
      case 's':
        c = ChampParseBlockAtom(I,c,cur_atom,cH_Sym,2,not_flag);
        atom_seen = true;
        break;
      case 'u':
        c = ChampParseBlockAtom(I,c,cur_atom,cH_Cu,2,not_flag);
        atom_seen = true;
        break;
      case 'l':
        c = ChampParseBlockAtom(I,c,cur_atom,cH_Cl,2,not_flag);
        atom_seen = true;
        break;
      default:
        c = ChampParseBlockAtom(I,c,cur_atom,cH_C,1,not_flag);
        atom_seen = true;
        break;
      }
      break;
    case 'D':
      switch(*(c+1)) {
        case 'y':
          c = ChampParseBlockAtom(I,c,cur_atom,cH_Sym,2,not_flag);
          atom_seen = true;
          break;
      default:
        num = ChampParseNumeral(c+1);
        if(num>=0) {
          if(not_flag) {
            I->Atom[cur_atom].neg_flag=true;
            I->Atom[cur_atom].not_degree|=num_to_degree[num];
          } else {
            I->Atom[cur_atom].pos_flag=true;
            I->Atom[cur_atom].degree|=num_to_degree[num];
          }
          c+=2;
        } else
          ok=false;
        break;
      }
      break;
    case 'E':
      switch(*(c+1)) {
      case 'r':
        c = ChampParseBlockAtom(I,c,cur_atom,cH_Sym,2,not_flag);
        atom_seen = true;
        break;
      case 'u':
        c = ChampParseBlockAtom(I,c,cur_atom,cH_Sym,2,not_flag);
        atom_seen = true;
        break;
      default:
        c = ChampParseBlockAtom(I,c,cur_atom,cH_E,1,not_flag);
        atom_seen = true;
        break;
      }
      break;
    case 'F':
      switch(*(c+1)) {
      case 'e':
        c = ChampParseBlockAtom(I,c,cur_atom,cH_Fe,2,not_flag);
        atom_seen = true;
        break;
      case 'r':
        c = ChampParseBlockAtom(I,c,cur_atom,cH_Sym,2,not_flag);
        atom_seen = true;
        break;
      default:
        c = ChampParseBlockAtom(I,c,cur_atom,cH_F,1,not_flag);
        atom_seen = true;
        break;
      }
      break;      
    case 'G':
      switch(*(c+1)) {
      case 'a':
        c = ChampParseBlockAtom(I,c,cur_atom,cH_Sym,2,not_flag);
        atom_seen = true;
        break;
      case 'd':
        c = ChampParseBlockAtom(I,c,cur_atom,cH_Sym,2,not_flag);
        atom_seen = true;
        break;
      case 'e':
        c = ChampParseBlockAtom(I,c,cur_atom,cH_Sym,2,not_flag);
        atom_seen = true;
        break;
      }
      break;
    case 'H': 
      switch(*(c+1)) {
      case 'e':
        c = ChampParseBlockAtom(I,c,cur_atom,cH_Sym,2,not_flag);
        atom_seen = true;
        break;
      case 'f':
        c = ChampParseBlockAtom(I,c,cur_atom,cH_Sym,2,not_flag);
        atom_seen = true;
        break;
      case 'g':
        c = ChampParseBlockAtom(I,c,cur_atom,cH_Sym,2,not_flag);
        atom_seen = true;
        break;
      case 'o':
        c = ChampParseBlockAtom(I,c,cur_atom,cH_Sym,2,not_flag);
        atom_seen = true;
        break;
      default:
        if(!atom_seen) {
          c = ChampParseBlockAtom(I,c,cur_atom,cH_H,1,not_flag);
          atom_seen = true;
        } else {
          num = ChampParseNumeral(c+1);
          if(num>=0) {
            c+=2;
          } else {
            num = 1;
            c++;
            while(*c) {
              if(*c!='H')
                break;
              else {
                num++;
                c++;
              }
            }
          }
          I->Atom[cur_atom].imp_hydro = num;
          I->Atom[cur_atom].tot_hydro = num;
          I->Atom[cur_atom].hydro_flag = true; 
          /* turn on hydrogen count matching for this atom */
          
        }
        break;
      }
      break;
    case 'I':
      switch(*(c+1)) {
      case 'n':
        c = ChampParseBlockAtom(I,c,cur_atom,cH_Sym,2,not_flag);
        atom_seen = true;
        break;
      case 'r':
        c = ChampParseBlockAtom(I,c,cur_atom,cH_Sym,2,not_flag);
        atom_seen = true;
        break;
      default:
        c = ChampParseBlockAtom(I,c,cur_atom,cH_I,1,not_flag);
        atom_seen = true;
        break;
      }
      break;      
    case 'J':
      c = ChampParseBlockAtom(I,c,cur_atom,cH_J,1,not_flag);
        atom_seen = true;
      break;
    case 'K':
      c = ChampParseBlockAtom(I,c,cur_atom,cH_K,1,not_flag);
        atom_seen = true;
      break;      
    case 'L':
      switch(*(c+1)) {
      case 'a':
        c = ChampParseBlockAtom(I,c,cur_atom,cH_Sym,2,not_flag);
        atom_seen = true;
        break;
      case 'i':
        c = ChampParseBlockAtom(I,c,cur_atom,cH_Sym,2,not_flag);
        atom_seen = true;
        break;
      case 'u':
        c = ChampParseBlockAtom(I,c,cur_atom,cH_Sym,2,not_flag);
        atom_seen = true;
        break;
      default:
        c = ChampParseBlockAtom(I,c,cur_atom,cH_L,1,not_flag);
        atom_seen = true;
      }
      break;      
    case 'M':
      switch(*(c+1)) {
      case 'n':
        c = ChampParseBlockAtom(I,c,cur_atom,cH_Sym,2,not_flag);
        atom_seen = true;
        break;
      case 'o':
        c = ChampParseBlockAtom(I,c,cur_atom,cH_Sym,2,not_flag);
        atom_seen = true;
        break;
      case 'g':
        c = ChampParseBlockAtom(I,c,cur_atom,cH_Mg,2,not_flag);
        atom_seen = true;
        break;
      default:
        c = ChampParseBlockAtom(I,c,cur_atom,cH_M,1,not_flag);
        atom_seen = true;
        break;
      }
      break;
    case 'N':
      switch(*(c+1)) {
      case 'a':
        c = ChampParseBlockAtom(I,c,cur_atom,cH_Na,2,not_flag);      
        atom_seen = true;
        break;
      case 'b':
        c = ChampParseBlockAtom(I,c,cur_atom,cH_Sym,2,not_flag);      
        atom_seen = true;
        break;
      case 'd':
        c = ChampParseBlockAtom(I,c,cur_atom,cH_Sym,2,not_flag);      
        atom_seen = true;
        break;
      case 'i':
        c = ChampParseBlockAtom(I,c,cur_atom,cH_Sym,2,not_flag);      
        atom_seen = true;
        break;
      default:
        c = ChampParseBlockAtom(I,c,cur_atom,cH_N,1,not_flag);
        atom_seen = true;
        break;
      }
      break;      
    case 'O':
      switch(*(c+1)) {
      case 's':
        c = ChampParseBlockAtom(I,c,cur_atom,cH_Sym,2,not_flag);      
        atom_seen = true;
        break;
      default:
        c = ChampParseBlockAtom(I,c,cur_atom,cH_O,1,not_flag);
        atom_seen = true;
        break;
      }
      break;      
    case 'P':
      switch(*(c+1)) {
      case 'b':
        c = ChampParseBlockAtom(I,c,cur_atom,cH_Sym,2,not_flag);      
        atom_seen = true;
        break;
      case 'd':
        c = ChampParseBlockAtom(I,c,cur_atom,cH_Sym,2,not_flag);      
        atom_seen = true;
        break;
      case 'o':
        c = ChampParseBlockAtom(I,c,cur_atom,cH_Sym,2,not_flag);      
        atom_seen = true;
        break;
      case 'r':
        c = ChampParseBlockAtom(I,c,cur_atom,cH_Sym,2,not_flag);      
        atom_seen = true;
        break;
      case 't':
        c = ChampParseBlockAtom(I,c,cur_atom,cH_Sym,2,not_flag);      
        atom_seen = true;
        break;
      default:
        c = ChampParseBlockAtom(I,c,cur_atom,cH_P,1,not_flag);
        atom_seen = true;
      }
      break;      
    case 'p': /* Pi system */
      if(not_flag) {
        I->Atom[cur_atom].neg_flag=true;
        I->Atom[cur_atom].not_class|=cH_Pi;
      } else {
        I->Atom[cur_atom].pos_flag=true;
        I->Atom[cur_atom].class|=cH_Pi;
      }
      c++;
      break;      
    case 'R':
      switch(*(c+1)) {
      case 'b':
        c = ChampParseBlockAtom(I,c,cur_atom,cH_Sym,2,not_flag);
        atom_seen = true;
        break;
      case 'e':
        c = ChampParseBlockAtom(I,c,cur_atom,cH_Sym,2,not_flag);
        atom_seen = true;
        break;
      case 'h':
        c = ChampParseBlockAtom(I,c,cur_atom,cH_Sym,2,not_flag);
        atom_seen = true;
        break;
      case 'u':
        c = ChampParseBlockAtom(I,c,cur_atom,cH_Sym,2,not_flag);
        atom_seen = true;
        break;
      default:
        c = ChampParseBlockAtom(I,c,cur_atom,cH_R,1,not_flag);
        atom_seen = true;
        break;
      }
      break;
    case 'r':
      num = ChampParseNumeral(c+1);
      if(num>=0) {
        if(not_flag) {
          I->Atom[cur_atom].neg_flag=true;
          I->Atom[cur_atom].not_cycle|=num_to_ring[num];
        } else {
          I->Atom[cur_atom].pos_flag=true;
          I->Atom[cur_atom].cycle|=num_to_ring[num];
        }
        c+=2;
      } else {
        if(not_flag) {
          I->Atom[cur_atom].neg_flag=true;
          I->Atom[cur_atom].not_cycle|=cH_Ring3|cH_Ring4|cH_Ring5|cH_Ring6|cH_Ring7|cH_Ring8;
        } else {
          I->Atom[cur_atom].pos_flag=true;
          I->Atom[cur_atom].cycle|=cH_Ring3|cH_Ring4|cH_Ring5|cH_Ring6|cH_Ring7|cH_Ring8;
        }
        c++;
      }
      break;
    case 'Q':
      c = ChampParseBlockAtom(I,c,cur_atom,cH_Q,1,not_flag);
        atom_seen = true;
      break;      
    case 'S':
      switch(*(c+1)) {
      case 'b':
        c = ChampParseBlockAtom(I,c,cur_atom,cH_Sym,2,not_flag);
        atom_seen = true;
        break;
      case 'c':
        c = ChampParseBlockAtom(I,c,cur_atom,cH_Sym,2,not_flag);
        atom_seen = true;
        break;
      case 'e':
        c = ChampParseBlockAtom(I,c,cur_atom,cH_Se,2,not_flag);
        atom_seen = true;
        break;
      case 'i':
        c = ChampParseBlockAtom(I,c,cur_atom,cH_Sym,2,not_flag);
        atom_seen = true;
        break;
      case 'm':
        c = ChampParseBlockAtom(I,c,cur_atom,cH_Sym,2,not_flag);
        atom_seen = true;
        break;
      case 'n':
        c = ChampParseBlockAtom(I,c,cur_atom,cH_Sym,2,not_flag);
        atom_seen = true;
        break;
      case 'r':
        c = ChampParseBlockAtom(I,c,cur_atom,cH_Sym,2,not_flag);
        atom_seen = true;
        break;
      default:
        c = ChampParseBlockAtom(I,c,cur_atom,cH_S,1,not_flag);
        atom_seen = true;
        break;
      }
      break;      
    case 'T':
      switch(*(c+1)) {
      case 'a':
        c = ChampParseBlockAtom(I,c,cur_atom,cH_Sym,2,not_flag);
        atom_seen = true;
        break;
      case 'b':
        c = ChampParseBlockAtom(I,c,cur_atom,cH_Sym,2,not_flag);
        atom_seen = true;
        break;
      case 'e':
        c = ChampParseBlockAtom(I,c,cur_atom,cH_Sym,2,not_flag);
        atom_seen = true;
        break;
      case 'i':
        c = ChampParseBlockAtom(I,c,cur_atom,cH_Sym,2,not_flag);
        atom_seen = true;
        break;
      case 'h':
        c = ChampParseBlockAtom(I,c,cur_atom,cH_Sym,2,not_flag);
        atom_seen = true;
        break;
      case 'l':
        c = ChampParseBlockAtom(I,c,cur_atom,cH_Sym,2,not_flag);
        atom_seen = true;
        break;
      case 'm':
        c = ChampParseBlockAtom(I,c,cur_atom,cH_Sym,2,not_flag);
        atom_seen = true;
        break;
      default:
        c = ChampParseBlockAtom(I,c,cur_atom,cH_T,1,not_flag);
        atom_seen = true;
      }
      break;    
    case 'V':
      switch(*(c+1)) {
      default:
        c = ChampParseBlockAtom(I,c,cur_atom,cH_Sym,1,not_flag);
        atom_seen = true;
      }
      break;    
    case 'v':
      num = ChampParseNumeral(c+1);
      if(num>=0) {
        if(not_flag) {
          I->Atom[cur_atom].neg_flag=true;
          I->Atom[cur_atom].not_valence|=num_to_valence[num];
        } else {
          I->Atom[cur_atom].pos_flag=true;
          I->Atom[cur_atom].degree|=num_to_valence[num];
        }
        c+=2;
      } else 
        ok=false;
      break;
    case 'W':
      switch(*(c+1)) {
      default:
        c = ChampParseBlockAtom(I,c,cur_atom,cH_Sym,1,not_flag);
        atom_seen = true;
        break;
      }
      break;
    case 'U':
      switch(*(c+1)) {
      default:
        c = ChampParseBlockAtom(I,c,cur_atom,cH_Sym,1,not_flag);
        atom_seen = true;
        break;
      }
      break;
    case 'Y':
      switch(*(c+1)) {
      case 'b':
        c = ChampParseBlockAtom(I,c,cur_atom,cH_Sym,2,not_flag);
        atom_seen = true;
        break;
      default:
        c = ChampParseBlockAtom(I,c,cur_atom,cH_Sym,1,not_flag);
        atom_seen = true;
        break;
      }
      break;
    case 'Z':
      switch(*(c+1)) {
      case 'r':
        c = ChampParseBlockAtom(I,c,cur_atom,cH_Sym,2,not_flag);
        atom_seen = true;
        break;
      case 'n':
        c = ChampParseBlockAtom(I,c,cur_atom,cH_Zn,2,not_flag);
        atom_seen = true;
        break;
      default:
        c = ChampParseBlockAtom(I,c,cur_atom,cH_Z,1,not_flag);
        atom_seen = true;
        break;
      }
      break;
    case '+':
      num = ChampParseNumeral(c+1);
      if(num>=0) {
        c+=2;
      } else {
        num = 1;
        c++;
        while(*c) {
          if(*c!='+')
            break;
          else {
            num++;
            c++;
          }
        }
      }
      if(not_flag) {
        I->Atom[cur_atom].neg_flag=true;
        switch(num) {
        case 0: I->Atom[cur_atom].not_charge|=cH_Neutral; break;
        case 1: I->Atom[cur_atom].not_charge|=cH_Cation; break;
        case 2: I->Atom[cur_atom].not_charge|=cH_Dication; break;
        case 3: I->Atom[cur_atom].not_charge|=cH_Trication; break;
        case 4: I->Atom[cur_atom].not_charge|=cH_Tetcation; break;
        case 5: I->Atom[cur_atom].not_charge|=cH_Pentcation; break;
          
        }
      } else {
        I->Atom[cur_atom].pos_flag=true;
        switch(num) {
        case 0: I->Atom[cur_atom].charge|=cH_Neutral; break;
        case 1: I->Atom[cur_atom].charge|=cH_Cation; break;
        case 2: I->Atom[cur_atom].charge|=cH_Dication; break;
        case 3: I->Atom[cur_atom].charge|=cH_Trication; break;
        case 4: I->Atom[cur_atom].charge|=cH_Tetcation; break;
        case 5: I->Atom[cur_atom].charge|=cH_Pentcation; break;
        }
      }
      break;
    case '-':
      num = ChampParseNumeral(c+1);
      if(num>=0) {
        c+=2;
      } else {
        num = 1;
        c++;
        while(*c) {
          if(*c!='-')
            break;
          else {
            num++;
            c++;
          }
        }
      }
      if(not_flag) {
        I->Atom[cur_atom].neg_flag=true;
        switch(num) {
        case 0: I->Atom[cur_atom].not_charge|=cH_Neutral; break;
        case 1: I->Atom[cur_atom].not_charge|=cH_Anion; break;
        case 2: I->Atom[cur_atom].not_charge|=cH_Dianion; break;
        case 3: I->Atom[cur_atom].not_charge|=cH_Trianion; break;
        case 4: I->Atom[cur_atom].not_charge|=cH_Tetanion; break;
        case 5: I->Atom[cur_atom].not_charge|=cH_Pentanion; break;
        }
      } else {
        I->Atom[cur_atom].pos_flag=true;
        switch(num) {
        case 0: I->Atom[cur_atom].charge|=cH_Neutral; break;
        case 1: I->Atom[cur_atom].charge|=cH_Anion; break;
        case 2: I->Atom[cur_atom].charge|=cH_Dianion; break;
        case 3: I->Atom[cur_atom].charge|=cH_Trianion; break;
        case 4: I->Atom[cur_atom].charge|=cH_Tetanion; break;
        case 5: I->Atom[cur_atom].charge|=cH_Pentanion; break;
        }
      }
      break;
    default:
      PRINTFB(FB_smiles_parsing,FB_errors)
        " champ: error parsing atom block at '%c' in: '%s'\n",*c,*c_ptr
        ENDFB;
      c++;
      break;
    }
  }
  *c_ptr = c;
  return ok;
}


#define cSym_Null       0
#define cSym_Atom       1
#define cSym_Bond       2
#define cSym_OpenScope  3
#define cSym_CloseScope 4
#define cSym_Mark       5
#define cSym_OpenBlock  6
#define cSym_CloseBlock 7
#define cSym_Separator  8
#define cSym_Qualifier  9

char *ChampPatToSmiVLA(CChamp *I,int index,char *vla,int mode)
{
  char *result;
  int n_atom;
  int cur_scope = 0;
  int cur_atom = 0;
  int cur_bond = 0;
  int a,c,l;
  int mark[MAX_RING];
  int mark_bond[MAX_RING];
  int i;
  int left_to_do = 0;
  int next_mark =1;
  int start_atom = 0;
  AtomBuffer buf;
  ListAtom *at1;
  ListBond *bd1;
  ListScope *scp1,*scp2;

  #define concat_s(s,c) {vla_check(result,char,l+(c)+1);strcpy(result+l,(s));l+=(c);}
  #define concat_c(s) {vla_check(result,char,l+2);result[l]=(s);result[l+1]=0;l++;}


  for(a=0;a<MAX_RING;a++) {
    mark[a]=0;
  }
  
  cur_atom=I->Pat[index].atom;
  n_atom=0;
  while(cur_atom) {
    n_atom++;
    at1 = I->Atom + cur_atom;
    at1->mark_tmpl = 0;
    cur_atom = at1->link;
  }

  /*guestimate size requirements, and allocate */

  if(!vla)
    vla_malloc(result,char,n_atom*4);
  else
    result = vla;
  result[0]=0;
  l = 0;
  
  start_atom = I->Pat[index].atom;
  while(start_atom) {
    if(!I->Atom[start_atom].mark_tmpl) {
      if(l) concat_c('.')
        
      cur_scope = ListElemNewZero(&I->Scope);
      I->Scope[cur_scope].atom = start_atom;
      I->Scope[cur_scope].bond = -1; /* signals start in new scope */
      while(cur_scope) {
        scp1 = I->Scope + cur_scope;
        cur_atom = scp1->atom;
        at1=I->Atom + cur_atom;
        
        PRINTFD(FB_smiles_creation) 
          " SmiToStrVLA: scope %d cur_atom %d base_bond %d\n",cur_scope,cur_atom,scp1->base_bond
          ENDFD;
        
        if(scp1->bond<0) { /* starting new scope, so print atom and continue */
          /* print bond, if required */
          if(scp1->base_bond) {
            c = ChampBondToString(I,scp1->base_bond,buf);
            concat_s(buf,c);
          }
          /* write atom string */
          at1->mark_tmpl = 1;
          c = ChampAtomToString(I,cur_atom,buf);
          concat_s(buf,c);
          /*          sprintf(buf,":%d:",cur_atom);
                      concat_s(buf,strlen(buf));*/

          switch(mode) {
          case 0:
            break;
          case 1:
            concat_s("<",1);
            concat_s(at1->name,strlen(at1->name));
            concat_s(">",1);
            break;
          }
          /* write opening marks */
          
          for(a=0;a<MAX_BOND;a++) {
            cur_bond = at1->bond[a];
            if(!cur_bond) break;
            bd1 = I->Bond+cur_bond;
            /* write cycle indicator if necessary */
            if(bd1->atom[0]!=cur_atom) {/* opposite direction -> explicit cycle */
              if(!I->Atom[bd1->atom[0]].mark_tmpl)
                {
                  if(mark[next_mark]) {
                    for(index=0;index<9;index++) {
                      if(!mark[index]) break;
                    }
                  } else {
                    index = next_mark++;
                  }
                  if(index<MAX_RING) {
                    mark[index]=bd1->atom[0]; /* save for matching other end of cycle */
                    mark_bond[index]=cur_bond;
                    c = ChampBondToString(I,cur_bond,buf);
                    concat_s(buf,c);
                    if(index<10) {
                      concat_c('0'+index);
                    } else {
                      sprintf(buf,"%%%d",index);
                      concat_s(buf,strlen(buf));
                    }
                  }
                }
            }
          }
          
          /* write closing marks */
          
          for(index=0;index<MAX_RING;index++) {
            if(mark[index]==cur_atom) {
              c = ChampBondToString(I,mark_bond[index],buf);
              concat_s(buf,c);
              if(index<10) {
                concat_c('0'+index);
              } else {
                sprintf(buf,"%%%d",index);
                concat_s(buf,strlen(buf));
              }
              mark[index]=0;
            }
          }
          
        }
        
        /* increment bond index index counter */
        
        scp1->bond++;
        
        /* now determine whether or not we need to create a new scope, 
           and figure out which bond to work on */
        
        i = scp1->bond;
        left_to_do = 0;
        cur_bond = 0;
        while(i<MAX_BOND) {
          if(!at1->bond[i]) break;
          bd1 = I->Bond + at1->bond[i];
          if(bd1->atom[0]==cur_atom) {
            if(!I->Atom[bd1->atom[1]].mark_tmpl) { /* not yet complete */ 
              if(!cur_bond) cur_bond = at1->bond[i];
              left_to_do++;
            }
          }
          i++;
        }
        
        PRINTFD(FB_smiles_creation) 
          " SmiToStrVLA: cur_atom %d left to do %d cur_bond %d\n",cur_atom,left_to_do,cur_bond
          ENDFD;
        
        if(left_to_do>1) { /* yes, we need a new scope */
          cur_scope = ListElemPush(&I->Scope,cur_scope);
          scp2 = I->Scope + cur_scope;
          scp2->base_bond = cur_bond;
          scp2->atom = I->Bond[cur_bond].atom[1];
          scp2->bond = -1;
          concat_c('(');
          scp2->paren_flag=true;
          PRINTFD(FB_smiles_creation) 
            " SmiToStrVLA: creating new scope vs old %d\n",cur_scope
            ENDFD;

        } else if(left_to_do) { /* no we do not, so just extend current scope */
          scp1->atom = I->Bond[cur_bond].atom[1];
          scp1->base_bond = cur_bond;
          scp1->bond = -1;
          
          PRINTFD(FB_smiles_creation) 
            " SmiToStrVLA: extending scope\n"
            ENDFD;
        } else { /* nothing attached, so just close scope */
          if(scp1->paren_flag)
            concat_c(')');
          cur_scope = ListElemPop(I->Scope,cur_scope);

          PRINTFD(FB_smiles_creation) 
            " SmiToStrVLA: closing scope\n"
            ENDFD;
        }
      }
    }
    start_atom = I->Atom[start_atom].link;
  }
    
  /* trim memory usage */
  vla_set_size(result,char,strlen(result)+1);
  return result;
}

static void ChampReassignLexPri(CChamp *I,int index)
{
  /* reassigns lexical priorities based on the current tree
  NOTE: unless stereo information is stored elsewhere, this trashes it 
  */
  int n_atom;
  int cur_scope = 0;
  int cur_atom = 0;
  int cur_bond = 0;
  int a;
  int mark[MAX_RING];
  int mark_bond[MAX_RING];
  int i;
  int left_to_do = 0;
  int next_mark =1;
  int start_atom = 0;
  ListAtom *at1;
  ListBond *bd1;
  ListScope *scp1,*scp2;
  int lex_pri = 0;
  
  for(a=0;a<MAX_RING;a++) {
    mark[a]=0;
  }
  
  cur_atom=I->Pat[index].atom;
  n_atom=0;
  while(cur_atom) {
    n_atom++;
    at1 = I->Atom + cur_atom;
    at1->mark_tmpl = 0;
    cur_atom = at1->link;
  }

  start_atom = I->Pat[index].atom;
  while(start_atom) {
    if(!I->Atom[start_atom].mark_tmpl) {
      lex_pri++;
      cur_scope = ListElemNewZero(&I->Scope);
      I->Scope[cur_scope].atom = start_atom;
      I->Scope[cur_scope].bond = -1; /* signals start in new scope */
      while(cur_scope) {
        scp1 = I->Scope + cur_scope;
        cur_atom = scp1->atom;
        at1=I->Atom + cur_atom;
        
        if(scp1->bond<0) { /* starting new scope, so print atom and continue */
          /* print bond, if required */
          if(scp1->base_bond) {
            bd1 = I->Bond+scp1->base_bond;
            lex_pri++;
            bd1->pri[0] = lex_pri; 
            bd1->pri[1] = lex_pri; 
          }
          /* write atom string */
          at1->mark_tmpl = 1;
          lex_pri ++;

          /* write opening marks */
          
          for(a=0;a<MAX_BOND;a++) {
            cur_bond = at1->bond[a];
            if(!cur_bond) break;
            bd1 = I->Bond+cur_bond;
            /* write cycle indicator if necessary */
            if(bd1->atom[0]!=cur_atom) {/* opposite direction -> explicit cycle */
              if(!I->Atom[bd1->atom[0]].mark_tmpl)
                {
                  if(mark[next_mark]) {
                    for(index=0;index<9;index++) {
                      if(!mark[index]) break;
                    }
                  } else {
                    index = next_mark++;
                  }
                  if(index<MAX_RING) {
                    mark[index]=bd1->atom[0]; /* save for matching other end of cycle */
                    mark_bond[index]=cur_bond;
                    lex_pri++;
                    bd1->pri[1] = lex_pri; /* from this atom's point of view, this is the priority location*/
                  }
                }
            }
          }
          
          /* write closing marks */
          
          for(index=0;index<MAX_RING;index++) {
            if(mark[index]==cur_atom) {
              lex_pri++;
              bd1 = I->Bond + mark_bond[index];
              bd1->pri[0] = lex_pri;
              mark[index]=0;
            }
          }
        }
        
        /* increment bond index index counter */
        
        scp1->bond++;
        
        /* now determine whether or not we need to create a new scope, 
           and figure out which bond to work on */
        
        i = scp1->bond;
        left_to_do = 0;
        cur_bond = 0;
        while(i<MAX_BOND) {
          if(!at1->bond[i]) break;
          bd1 = I->Bond + at1->bond[i];
          if(bd1->atom[0]==cur_atom) {
            if(!I->Atom[bd1->atom[1]].mark_tmpl) { /* not yet complete */ 
              if(!cur_bond) cur_bond = at1->bond[i];
              left_to_do++;
            }
          }
          i++;
        }
        
        if(left_to_do>1) { /* yes, we need a new scope */
          cur_scope = ListElemPush(&I->Scope,cur_scope);
          scp2 = I->Scope + cur_scope;
          scp2->base_bond = cur_bond;
          scp2->atom = I->Bond[cur_bond].atom[1];
          scp2->bond = -1;
          lex_pri++;
          scp2->paren_flag=true;

        } else if(left_to_do) { /* no we do not, so just extend current scope */
          scp1->atom = I->Bond[cur_bond].atom[1];
          scp1->base_bond = cur_bond;
          scp1->bond = -1;

        } else { /* nothing attached, so just close scope */
          if(scp1->paren_flag)
            lex_pri++;
          cur_scope = ListElemPop(I->Scope,cur_scope);
        }
      }
    }
    start_atom = I->Atom[start_atom].link;
  }
  
}

void ChampGeneralize(CChamp *I,int index)
{
  ListBond *bd1;
  int cur_bond = 0;

  ChampPrepareTarget(I,index);

  cur_bond=I->Pat[index].bond; 
  while(cur_bond) {
    bd1 = I->Bond + cur_bond;
    if(bd1->class&cH_Aromatic) {
      bd1->order=cH_NoOrder; /* strip specific bond order
                              * from aromatic bonds */
      bd1->class=cH_Pi; /* make all bonds Pi bonds */
    }
    cur_bond = bd1->link;
  }

}

static void ChampStereoToInternal(CChamp *I, int index)
{
  
  {
     int n_atom;
     int cur_atom = 0;
     int cur_bond = 0;
     int start_atom = 0;
     int a;
     int order_handedness;
     int pri_handedness;
     ListAtom *at1;
     ListBond *bd1;
     int n_bond;
     int ati[MAX_BOND],pri[MAX_BOND];

    /* first we need to sweep the molecule in order and assign lexical
       priorities to bonds based on when atoms are written out */
    
    cur_atom=I->Pat[index].atom;
    n_atom=0;
    while(cur_atom) {
      n_atom++;
      at1 = I->Atom + cur_atom;
      at1->mark_tmpl = 0;
      at1->stereo_internal = 0; /* clear stereo orientation */
      cur_atom = at1->link;
    }

    start_atom = I->Pat[index].atom;
    while(start_atom) {
      if(!I->Atom[start_atom].mark_tmpl) {
        
        cur_atom = start_atom;
        at1=I->Atom + cur_atom;
        
        at1->mark_tmpl = 1;
        
        if(at1->stereo) {
          
          n_bond = 0;
          for(a=0;a<MAX_BOND;a++) {
            cur_bond = at1->bond[a];
            if(!cur_bond) break;
            n_bond++;
          }
          
          if(n_bond==4) {
            n_bond = 0;
            for(a=0;a<MAX_BOND;a++) {
              cur_bond = at1->bond[a];
              if(!cur_bond) break;
              bd1 = I->Bond+cur_bond;
              if(bd1->atom[0]==cur_atom) {
                pri[n_bond] = bd1->pri[0];
                ati[n_bond] = bd1->atom[1];
              } else {
                pri[n_bond] = bd1->pri[1];
                ati[n_bond] = bd1->atom[0];
              }
              n_bond++;
            }
            
            { /* get coordinates into absolute atom order */
              int idx[MAX_BOND];
              SortIntIndex(4,pri,idx);
              pri_handedness = ChiralHandedness(idx);
              SortIntIndex(4,ati,idx);
              order_handedness = ChiralHandedness(idx);
            }

            if(pri_handedness == order_handedness)
              at1->stereo_internal = at1->stereo;
            else
              at1->stereo_internal = -at1->stereo;
          }
        }
        start_atom = I->Atom[start_atom].link;
      }
    }
  }
}

static void ChampStereoFromInternal(CChamp *I, int index)
{
  
  {
     int n_atom;
     int cur_atom = 0;
     int cur_bond = 0;
     int start_atom = 0;
     int pri_handedness,order_handedness;
     int a;
     ListAtom *at1;
     ListBond *bd1;
     int n_bond;
     int ati[MAX_BOND],pri[MAX_BOND];


    /* first we need to sweep the molecule in order and assign lexical
       priorities to bonds based on when atoms are written out */
    
    cur_atom=I->Pat[index].atom;
    n_atom=0;
    while(cur_atom) {
      n_atom++;
      at1 = I->Atom + cur_atom;
      at1->mark_tmpl = 0;
      at1->stereo = 0; /* clear stereo orientation */
      cur_atom = at1->link;
    }

    start_atom = I->Pat[index].atom;
    while(start_atom) {
      if(!I->Atom[start_atom].mark_tmpl) {
        
        cur_atom = start_atom;
        at1=I->Atom + cur_atom;
        
        at1->mark_tmpl = 1;
        
        if(at1->stereo_internal) {
          
          n_bond = 0;
          for(a=0;a<MAX_BOND;a++) {
            cur_bond = at1->bond[a];
            if(!cur_bond) break;
            n_bond++;
          }
          
          if(n_bond==4) {
            n_bond = 0;
            for(a=0;a<MAX_BOND;a++) {
              cur_bond = at1->bond[a];
              if(!cur_bond) break;
              bd1 = I->Bond+cur_bond;
              if(bd1->atom[0]==cur_atom) {
                pri[n_bond] = bd1->pri[0];
                ati[n_bond] = bd1->atom[1];
              } else {
                pri[n_bond] = bd1->pri[1];
                ati[n_bond] = bd1->atom[0];
              }
              n_bond++;
            }
            
            { /* get coordinates into absolute atom order */
              int idx[MAX_BOND];
              SortIntIndex(4,pri,idx);
              pri_handedness = ChiralHandedness(idx);
              SortIntIndex(4,ati,idx);
              order_handedness = ChiralHandedness(idx);
            }
          
            if(pri_handedness==order_handedness)
              at1->stereo = at1->stereo_internal;
            else
              at1->stereo = -at1->stereo_internal;
          }
        }
        start_atom = I->Atom[start_atom].link;
      }
    }
  }
}


void ChampOrientBonds(CChamp *I,int index) 
     /* This destroys stereo information...right? */
{
  /* do a prepatory walk through the molecule to figure out how to minimize 
     explicit connections */

  int n_atom;
  int cur_scope = 0;
  int cur_atom = 0;
  int cur_bond = 0;
  int a,tmp;
  int i;
  int left_to_do = 0;
  int start_atom = 0;
  int last_atom = 0;
  ListAtom *at1;
  ListBond *bd1;
  ListScope *scp1,*scp2;

  ChampStereoToInternal(I,index); 
  /* convert stereo information into bond-order-independent internal representation */

  cur_atom=I->Pat[index].atom;
  n_atom=0;
  while(cur_atom) {
    n_atom++;
    at1 = I->Atom + cur_atom;
    at1->mark_tmpl = 0;
    cur_atom = at1->link;
  }

  cur_bond=I->Pat[index].bond; /* clear markings... */
  while(cur_bond) {
    bd1 = I->Bond + cur_bond;
    bd1->mark_tmpl = 0;
    cur_bond = bd1->link;
  }

  /* make sure that the start atom is not stereo-specified */

  last_atom = 0;
  start_atom = I->Pat[index].atom;
  while(start_atom) {
    if(!I->Atom[start_atom].stereo)
      break;
    last_atom = start_atom;
    start_atom = I->Atom[start_atom].link;
  }
  
  if(last_atom&&start_atom) { /* start atom is stereo, so we've got to change things */
    int old_first_atom = I->Pat[index].atom;
    I->Pat[index].atom = start_atom;
    I->Atom[last_atom].link = I->Atom[start_atom].link; /* excise this atom from the list */
    I->Atom[start_atom].link = old_first_atom;
  }
    

  start_atom = I->Pat[index].atom;
  while(start_atom) {
    if(!I->Atom[start_atom].mark_tmpl) {
      cur_scope = ListElemNewZero(&I->Scope);
      I->Scope[cur_scope].atom = start_atom;
      I->Scope[cur_scope].bond = -1; /* signals start in new scope */
      while(cur_scope) {
        scp1 = I->Scope + cur_scope;
        cur_atom = scp1->atom;
        at1=I->Atom + cur_atom;
        
        if(scp1->bond<0) { /* starting new scope, so print atom and continue */
          
          /* mark atom */
          at1->mark_tmpl = 1;
          
          /* reorder atom's bonds (if necessary) */
          
          for(a=0;a<MAX_BOND;a++) {
            cur_bond = at1->bond[a];
            if(!cur_bond) break;
            bd1 = I->Bond+cur_bond;
            
            if(!bd1->mark_tmpl) { /* is this the first time we've seen this bond? */
              bd1->mark_tmpl=1;
              if(bd1->atom[0]!=cur_atom) { /* reorient bond to mimize explicit cycles... */
                tmp=bd1->atom[0];
                bd1->atom[0]=bd1->atom[1];
                bd1->atom[1]=tmp;
              }
            } else { /* hmmm...not the first time...so how did we get here? */
              if(bd1->atom[0]!=scp1->base_atom) { /* not this route so... */
                tmp=bd1->atom[0];  /* reorient bond to capture explicit cycle */
                bd1->atom[0]=bd1->atom[1];
                bd1->atom[1]=tmp;
              }
            }
          }
        }
        
        /* increment bond index index counter */
        
        scp1->bond++;
        
        /* now determine whether or not we need to create a new scope, 
           and figure out which bond to work on */
        
        i = scp1->bond;
        left_to_do = 0;
        cur_bond = 0;
        while(i<MAX_BOND) {
          if(!at1->bond[i]) break;
          bd1 = I->Bond + at1->bond[i];
          
          if(!bd1->mark_tmpl) { /* is this the first time we've seen this bond? */
            bd1->mark_tmpl=1;
            if(bd1->atom[0]!=cur_atom) { /* reorient bond to mimize explicit cycles... */
              tmp=bd1->atom[0];
              bd1->atom[0]=bd1->atom[1];
              bd1->atom[1]=tmp;
            }
          }
          
          if(bd1->atom[0]==cur_atom) {
            if(!I->Atom[bd1->atom[1]].mark_tmpl) { /* not yet complete */ 
              if(!cur_bond) cur_bond = at1->bond[i];
              left_to_do++;
            }
          }
          i++;
        }
        
        if(left_to_do>1) { /* yes, we need a new scope */
          cur_scope = ListElemPush(&I->Scope,cur_scope);
          scp2 = I->Scope + cur_scope;
          scp2->base_bond = cur_bond;
          scp2->atom = I->Bond[cur_bond].atom[1];
          scp2->base_atom = cur_atom;
          scp2->bond = -1;
        } else if(left_to_do) { /* no we do not, so just extend current scope */
          scp1->atom = I->Bond[cur_bond].atom[1];
          scp1->base_bond = cur_bond;
          scp1->base_atom = cur_atom;
          scp1->bond = -1;
        } else { /* nothing attached, so just close scope */
          cur_scope = ListElemPop(I->Scope,cur_scope);
        }
      }
    }
    start_atom = I->Atom[start_atom].link;
  }
  ChampReassignLexPri(I,index);
  ChampStereoFromInternal(I,index);

}

#define R_SMALL 0.0000001F

static double sqrt1f(float f) { /* no good as a macro because f is used twice */
  if(f>0.0F)
	 return(sqrt(f));
  else
	 return 0.0;
}

static double length3f ( float *v1 )
{
  return(sqrt1f((v1[0]*v1[0]) + 
					 (v1[1]*v1[1]) + 
					 (v1[2]*v1[2])));
} 

static void normalize3f( float *v1 )
{
	double vlen = length3f(v1);
	if(vlen > R_SMALL)
	{
		float	inV	= (float)(1.0F / vlen);
		v1[0] *= inV;
		v1[1] *= inV;
		v1[2] *= inV;
	}
	else
	{
		v1[0]=v1[1]=v1[2]=0.0F;
	}
} 

static void remove_component3f ( float *v1, float *unit, float *result)
{
  float dot;

  dot = v1[0]*unit[0] + v1[1]*unit[1] + v1[2]*unit[2];
  result[0]=v1[0]-unit[0]*dot;
  result[1]=v1[1]-unit[1]*dot;
  result[2]=v1[2]-unit[2]*dot;  
}

static float dot_product3f ( float *v1, float *v2 )
{
  return( v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2]);
}

static void subtract3f ( float *v1, float *v2, float *v3 )
{
  v3[0]=v1[0]-v2[0];
  v3[1]=v1[1]-v2[1];
  v3[2]=v1[2]-v2[2];
}

static void cross_product3f ( float *v1, float *v2, float *cross )
{
  cross[0] = (v1[1]*v2[2]) - (v1[2]*v2[1]);
  cross[1] = (v1[2]*v2[0]) - (v1[0]*v2[2]);
  cross[2] = (v1[0]*v2[1]) - (v1[1]*v2[0]);
}



void ChampDetectChirality(CChamp *I,int index)
{


  ChampReassignLexPri(I,index);

  /*  ChampPatDump(I,index);*/
  
  /* now, find the chiral atoms and geometrically determine their handedness 
     based on their coordinates */

  {
     int n_atom;
     int cur_atom = 0;
     int cur_bond = 0;
     int start_atom = 0;
     int a;
     ListAtom *at1;
     ListBond *bd1;
     int pri[MAX_BOND];
     int n_bond;
     int ati[MAX_BOND];
     float *vc, *v[MAX_BOND],vd[MAX_BOND][3],vr[MAX_BOND][3],vt[3];

    /* first we need to sweep the molecule in order and assign lexical
       priorities to bonds based on when atoms are written out */
    
    cur_atom=I->Pat[index].atom;
    n_atom=0;
    while(cur_atom) {
      n_atom++;
      at1 = I->Atom + cur_atom;
      at1->mark_tmpl = 0;
      at1->stereo = 0; /* clear stereo orientation */
      cur_atom = at1->link;
    }

    start_atom = I->Pat[index].atom;
    while(start_atom) {
      if(!I->Atom[start_atom].mark_tmpl) {
        
        cur_atom = start_atom;
        at1=I->Atom + cur_atom;
        
        /*        printf("name: %s\n",at1->name);*/
        at1->mark_tmpl = 1;

        n_bond = 0;
        for(a=0;a<MAX_BOND;a++) {
          cur_bond = at1->bond[a];
          if(!cur_bond) break;
          n_bond++;
        }

        if(n_bond==4) {
        
          /* get the relative lexical priorities for each bond */

          n_bond = 0;
          for(a=0;a<MAX_BOND;a++) {
            cur_bond = at1->bond[a];
            if(!cur_bond) break;
            bd1 = I->Bond+cur_bond;
            if(bd1->atom[0]==cur_atom) {
              pri[n_bond] = bd1->pri[0];
              ati[n_bond] = bd1->atom[1];
            } else {
              pri[n_bond] = bd1->pri[1];
              ati[n_bond] = bd1->atom[0];
            }
            n_bond++;
          }

          { /* get coordinates into lexical rotational order */
            int idx[MAX_BOND];
            SortIntIndex(4,pri,idx);
            /*
            printf("pri: %d %d %d %d\n",pri[0],pri[1],pri[2],pri[3]);
            printf("idx: %d %d %d %d\n",idx[0],idx[1],idx[2],idx[3]);
            printf("ord: %d %d %d %d\n",pri[idx[0]],pri[idx[1]],pri[idx[2]],pri[idx[3]]);

            printf("%s %s %s %s\n",
                   I->Atom[ati[idx[0]]].name,
                   I->Atom[ati[idx[1]]].name,
                   I->Atom[ati[idx[2]]].name,
                   I->Atom[ati[idx[3]]].name);
            */
     
            v[0]=I->Atom[ati[idx[0]]].coord;
            v[1]=I->Atom[ati[idx[1]]].coord;
            v[2]=I->Atom[ati[idx[2]]].coord;
            v[3]=I->Atom[ati[idx[3]]].coord;
          }
            
          vc = at1->coord;
          subtract3f(v[0],vc,vd[0]);
          subtract3f(v[1],vc,vd[1]);
          subtract3f(v[2],vc,vd[2]);
          subtract3f(v[3],vc,vd[3]);
         
          normalize3f(vd[0]);
          
          remove_component3f(vd[1],vd[0],vr[1]);
          remove_component3f(vd[2],vd[0],vr[2]);
          remove_component3f(vd[3],vd[0],vr[3]);
          
          cross_product3f(vr[1],vr[2],vt);
          normalize3f(vt);

          /*          printf("%s %8.3f\n",at1->name,dot_product3f(vd[0],vt));*/
          if(dot_product3f(vd[0],vt)>0.0F) {
            at1->stereo = cH_Anticlock;
          } else {
            at1->stereo = cH_Clockwise;
          } 
        }
      start_atom = I->Atom[start_atom].link;
      }
    }
  }
}
  
  
  void ChampAtomDump(CChamp *I,int index)
{
  char buf[3];
  ChampAtomToString(I,index,buf);
  printf("%s",buf);
}

void ChampPatDump(CChamp *I,int index)
{
  int cur_atom;
  int cur_bond;
  int a;
  AtomBuffer buf;

  ListAtom *at;
  ListBond *bd;
  cur_atom = I->Pat[index].atom;
  while(cur_atom) {
    at = I->Atom+cur_atom;
    ChampAtomToString(I,cur_atom,buf);
    printf(" atom %d %3s 0x%08x nam: %s res: %s sym: %s\n",cur_atom,buf,at->atom,at->name,at->residue,at->symbol);
    printf("        cy: %x",at->cycle);
    printf(" cl: %x v: %02x D: %02x ch: %02x cy: %d st: %d ih: %d hy: %d hf: %d bo: ",
           at->class,at->valence,at->degree,at->charge,at->cycle,
           at->stereo,at->imp_hydro,at->tot_hydro,at->hydro_flag);
    for(a=0;a<MAX_BOND;a++) {
      if(!at->bond[a]) break;
      else printf("%d ",at->bond[a]);
    }
    printf("\n");
    printf("        pf: %d nf: %d !at %d !ch %d !cy %d !cl %d !D %d !v %d\n",
           at->pos_flag,at->neg_flag,at->not_atom,
           at->not_charge,at->not_cycle,at->not_class,
           at->not_degree,at->not_valence);

    cur_atom = I->Atom[cur_atom].link;
  }
  cur_bond = I->Pat[index].bond;
  while(cur_bond) {
    bd = I->Bond+cur_bond;
    printf(" bond %d 0x%01x atoms %d %d order 0x%01x cycle %x dir %d class %x pri: %d %d\n",
           cur_bond,bd->order,bd->atom[0],bd->atom[1],bd->order,bd->cycle,
           bd->direction,bd->class,bd->pri[0],bd->pri[1]);
    cur_bond = I->Bond[cur_bond].link;
  }
  fflush(stdout);
}

int ChampBondToString(CChamp *I,int index,char *buf) 
{
  ListBond *bd;

  if(index) {
    bd=I->Bond + index;
    switch(bd->order) {
    case cH_Single: buf[0]=0; break;
    case cH_Double: buf[0]='='; buf[1]=0; break;
    case cH_Triple: buf[0]='#'; buf[1]=0; break;
    }
  } else 
    buf[0]=0;
  return(strlen(buf));
}

int ChampAtomToString(CChamp *I,int index,char *buf) 
{
  /* first, determine whether this is a trivial or block atom */
  int c;
  int mask;
  int a;
  int trivial=true;
  ListAtom *at;
  
  buf[0]=0;
  if(index) {
    at=I->Atom + index;

  /* the following are non-trivial */

    trivial = trivial && (!at->neg_flag) && (!at->hydro_flag);
    trivial = trivial && !(at->charge&(cH_Cation|cH_Dication|
                                       cH_Anion|cH_Dianion|
                                       cH_Trication|cH_Trianion|
                                       cH_Tetcation|cH_Tetanion|
                                       cH_Pentcation|cH_Pentanion
                                       ));
    trivial = trivial && !(at->cycle|at->valence|at->degree);
    trivial = trivial && !(at->class&cH_Aliphatic);

    trivial = trivial && !((at->atom!=cH_Any) &&
                           at->atom&( cH_Na | cH_K  | cH_Ca | cH_Mg | 
                                      cH_Zn | cH_Fe | cH_Cu | cH_Se |
                                      cH_A  | cH_E  | cH_G  | cH_J  |
                                      cH_L  | cH_M  | cH_Q  | cH_R  |
                                      cH_T | cH_X | cH_Z ));

    trivial = trivial && (at->stereo==0);

    if(trivial&&(at->atom!=cH_Any)) {
      /* check number of atoms represented */
      c = 0;
      mask = 1;
      for(a=0;a<32;a++) {
        if(mask&at->atom) {
          c++;
          if(c>1) {
            trivial=false;
            break;
          }
        }
        mask=mask*2;
      }
    }

  
    if(trivial) {
      if(at->class&cH_Aromatic) {
        switch(at->atom) {      
        case cH_C: buf[0]='c'; buf[1]=0; break;
        case cH_N: buf[0]='n'; buf[1]=0; break;
        case cH_O: buf[0]='o'; buf[1]=0; break;
        case cH_S: buf[0]='s'; buf[1]=0; break;
        default:
          trivial=false; break;
        }
      } else {
        switch(at->atom) {
        case cH_Any: buf[0]='*'; buf[1]=0; break;
        case cH_NotH: buf[0]='?'; buf[1]=0; break;
        case cH_B: buf[0]='B'; buf[1]=0; break;
        case cH_C: buf[0]='C'; buf[1]=0; break;
        case cH_N: buf[0]='N'; buf[1]=0; break;
        case cH_O: buf[0]='O'; buf[1]=0; break;
          /*        case cH_H: buf[0]='H'; buf[1]=0; break;*/
        case cH_H: strcpy(buf,"[H]"); break;
        case cH_S: buf[0]='S'; buf[1]=0; break;
        case cH_P: buf[0]='P'; buf[1]=0; break;
        case cH_F: buf[0]='F'; buf[1]=0; break;
        case cH_Cl: buf[0]='C'; buf[1]='l'; buf[2]=0; break;
        case cH_Br: buf[0]='B'; buf[1]='r'; buf[2]=0; break;
        case cH_I:  buf[0]='I'; buf[1]=0; break;
        default: 
          trivial = false; break;
        }
      }
    }
    if(!trivial) {
      strcat(buf,"[");
      if(at->atom==cH_Sym) {
        strcat(buf,at->symbol);
      } else {
        switch(at->atom) {
        case cH_Any: strcat(buf,"*"); break;
        case cH_NotH: strcat(buf,"?"); break; 
        case cH_B: strcat(buf,"B"); break;
        case cH_C: strcat(buf,"C"); break;
        case cH_N: strcat(buf,"N"); break;
        case cH_O: strcat(buf,"O"); break;
        case cH_H: strcat(buf,"H"); break;
        case cH_S: strcat(buf,"S"); break;
        case cH_P: strcat(buf,"P"); break;
        case cH_F: strcat(buf,"F"); break;
        case cH_Cl: strcat(buf,"Cl"); break;
        case cH_Br: strcat(buf,"Br"); break;
        case cH_I: strcat(buf,"I"); break;
        case cH_Na: strcat(buf,"Na"); break;
        case cH_K: strcat(buf,"K"); break;
        case cH_Ca: strcat(buf,"Ca"); break;
        case cH_Mg: strcat(buf,"Mg"); break;
        case cH_Fe: strcat(buf,"Fe"); break;
        case cH_Zn: strcat(buf,"Zn"); break;
        case cH_Cu: strcat(buf,"Cu"); break;
        case cH_Se: strcat(buf,"Se"); break;
        case cH_A: strcat(buf,"A"); break;
        case cH_E: strcat(buf,"E"); break;
        case cH_G: strcat(buf,"G"); break;
        case cH_J: strcat(buf,"J"); break;
        case cH_L: strcat(buf,"L"); break;
        case cH_M: strcat(buf,"M"); break;
        case cH_Q: strcat(buf,"Q"); break;
        case cH_R: strcat(buf,"R"); break;
        case cH_T: strcat(buf,"T"); break;
        case cH_X: strcat(buf,"X"); break;
        case cH_Z: strcat(buf,"Z"); break;
        default: sprintf(buf,"%x",at->atom); break;
        }
      }
      switch(at->stereo) {
      case cH_Anticlock:
        strcat(buf,"@");
        break;
      case cH_Clockwise:
        strcat(buf,"@@");
        break;
      }
      if(at->imp_hydro) {
        switch(at->imp_hydro) {
        case 0:
          break;
        case 1:
          strcat(buf,"H");
          break;
        case 2:
          strcat(buf,"H2");
          break;
        case 3:
          strcat(buf,"H3");
          break;
        case 4:
          strcat(buf,"H4");
          break;
        }
      }
      if(at->charge) {
        switch(at->charge) {
        case cH_Cation:
          strcat(buf,"+");
          break;
        case cH_Anion:
          strcat(buf,"-");
          break;
        case cH_Dication:
          strcat(buf,"++");
          break;
        case cH_Dianion:
          strcat(buf,"--");
          break;
        case cH_Trication:
          strcat(buf,"+3");
          break;
        case cH_Trianion:
          strcat(buf,"-3");
          break;
        case cH_Tetcation:
          strcat(buf,"+4");
          break;
        case cH_Tetanion:
          strcat(buf,"-4");
          break;
        case cH_Pentcation:
          strcat(buf,"+5");
          break;
        case cH_Pentanion:
          strcat(buf,"-5");
          break;
        }
      }
      strcat(buf,"]");
    }
  } else 
    buf[0]=0;
  return(strlen(buf));
}


int ChampAddBondToAtom(CChamp *I,int atom_index,int bond_index) 
{
  int i;
  ListAtom *at1;
  int ok=true;
  at1 = I->Atom+atom_index;
  i=0;
  while(at1->bond[i]) i++; /* flawed */
  if(i<MAX_BOND) {
    at1->bond[i] = bond_index;
  } else {
    PRINTFB(FB_smiles_parsing,FB_errors)
      " champ: MAX_BOND exceeded...\n"
      ENDFB;
    ok=false;
  }
  return(ok);
}

int ChampSmiToPat(CChamp *I,char *c) 
{ /* returns root atom of list */
  int mark[MAX_RING]; /* ring marks 0-9 */
  int mark_pri[MAX_RING]; /* lexical priority of mark */
  int stack = 0; /* parenthetical scopes */
  int base_atom = 0;
  int last_atom = 0;
  int last_bond = 0;
  int atom_list = 0;
  int bond_list = 0;
  int bond_flag = false;
  int cur_atom = 0;
  int cur_bond = 0;
  int mark_code = 0;
  int result = 0;
  int sym;
  int ok = true;
  unsigned int bond_tags = 0;
  unsigned int bond_not_tags = 0;
  int a;
  int not_bond = false;
  int lex_pri = 0;
  char *orig_c=c;

#define save_bond() { if(last_bond) {I->Bond[last_bond].link=cur_bond;}\
          else {bond_list=cur_bond;}\
          last_bond = cur_bond;\
          cur_bond = ListElemNewZero(&I->Bond);}
  
#define save_atom() { if(last_atom) {I->Atom[last_atom].link=cur_atom;}\
          else {atom_list=cur_atom;}\
          last_atom = cur_atom;\
          cur_atom = ListElemNewZero(&I->Atom);}

  PRINTFD(FB_smiles_parsing) 
    " ChampSmiToPat: input '%s'\n",c
    ENDFD;
  
  for(a=0;a<MAX_RING;a++)
    mark[a]=0;
  cur_atom = ListElemNewZero(&I->Atom);
  cur_bond = ListElemNewZero(&I->Bond);
  
  lex_pri = 0;
  while((*c)&&ok) {
    lex_pri++;
    PRINTFD(FB_smiles_parsing) 
      " parsing: '%c' at %p\n",*c,c
      ENDFD;
    sym = cSym_Null;
    /* ============ ROOT LEVEL PARSTING ============ */
    if(((*c)>='0')&&((*c)<='9')) {
      sym = cSym_Mark;
      mark_code = (*c)-'0';
      c++;
    } else {
      switch(*c) {
      /* standard, implicit atoms, with lowest normal valences
       * B(3), C(4), N(3,5), O(2), P(3,5), S(2,4,6), F(1), Cl(1), Br(1), I(1) */
      case 'C':
        switch(*(c+1)) {
        case 'l':
        case 'L': /* be tolerate at the root level, but not withing blocks...*/
          c = ChampParseAliphaticAtom(I,c,cur_atom,cH_Cl,2,false);
          sym = cSym_Atom;
          break;
        default:
          c = ChampParseAliphaticAtom(I,c,cur_atom,cH_C,1,true);
          sym = cSym_Atom;
          PRINTFD(FB_smiles_parsing) 
            " parsed: %p\n",c
            ENDFD;
          break;
        }
        break;
      case '<': /* tag index/list */
        if(bond_flag) {
          c = ChampParseTag(I,c,&bond_tags,&bond_not_tags,&ok);
        } else {
          if(base_atom) {
            c = ChampParseTag(I,c,&I->Atom[base_atom].tag,
                              &I->Atom[base_atom].not_tag,&ok);
          } else ok=false;
        }
        sym = cSym_Qualifier;
        break;
      case '*': /* nonstandard? */
        c = ChampParseAliphaticAtom(I,c,cur_atom,cH_Any,1,false);
        sym = cSym_Atom;
        break;
      case '?': /* nonstandard */
        c = ChampParseAliphaticAtom(I,c,cur_atom,cH_NotH,1,false);
        sym = cSym_Atom;
        break;
      case 'H': /* nonstandard */
        c = ChampParseAliphaticAtom(I,c,cur_atom,cH_H,1,false);
        sym = cSym_Atom;
        break;
      case 'N':
        c = ChampParseAliphaticAtom(I,c,cur_atom,cH_N,1,true);
        sym = cSym_Atom;
        break;      
      case 'O':
        c = ChampParseAliphaticAtom(I,c,cur_atom,cH_O,1,true);
        sym = cSym_Atom;
        break;      
      case 'B':
        switch(*(c+1)) {
        case 'r':
        case 'R':
          c = ChampParseAliphaticAtom(I,c,cur_atom,cH_Br,2,true);
          sym = cSym_Atom;
          break;
        default:
          c = ChampParseAliphaticAtom(I,c,cur_atom,cH_B,1,true);
          sym = cSym_Atom;
          break;
        }
        break;
      case 'P':
        c = ChampParseAliphaticAtom(I,c,cur_atom,cH_P,1,true);
        sym = cSym_Atom;
        break;      
      case 'S':
        c = ChampParseAliphaticAtom(I,c,cur_atom,cH_S,1,true);
        sym = cSym_Atom;
        break;      
      case 'F':
        c = ChampParseAliphaticAtom(I,c,cur_atom,cH_F,1,true);
        sym = cSym_Atom;
        break;      
      case 'I':
        c = ChampParseAliphaticAtom(I,c,cur_atom,cH_I,1,true);
        sym = cSym_Atom;
        break;      
        /* standard implicit aromatic atoms */
      case 'c':
        c = ChampParseAromaticAtom(I,c,cur_atom,cH_C,1,true);
        sym = cSym_Atom;
        break;
      case 'n':
        c = ChampParseAromaticAtom(I,c,cur_atom,cH_N,1,true);
        sym = cSym_Atom;
        break;
      case 'o':
        c = ChampParseAromaticAtom(I,c,cur_atom,cH_O,1,true);
        sym = cSym_Atom;
        break;
      case 's':
        c = ChampParseAromaticAtom(I,c,cur_atom,cH_S,1,true);
        sym = cSym_Atom;
        break;
      case ';':
        c++;
        not_bond=false;
        sym = cSym_Qualifier;
        break;
      case ',':
        c++;
        sym = cSym_Qualifier;
        break;
      case '!':
        c++;
        not_bond=true;
        sym = cSym_Qualifier;
        break;
      case '-':
        c++;
        if(not_bond) 
          I->Bond[cur_bond].not_order |= cH_Single;
        else 
          I->Bond[cur_bond].order |= cH_Single;
        sym = cSym_Bond;
        break;
      case '/':
        c++;
        if(not_bond) 
          I->Bond[cur_bond].not_order |= cH_Single;
        else 
          I->Bond[cur_bond].order |= cH_Single;
        sym = cSym_Bond;
        I->Bond[cur_bond].direction = cH_Up;
        break;
      case '\\':
        c++;
        if(not_bond) 
          I->Bond[cur_bond].not_order |= cH_Single;
        else 
          I->Bond[cur_bond].order |= cH_Single;
        sym = cSym_Bond;
        I->Bond[cur_bond].direction = cH_Down;
        break;
      case '=':
        c++;
        if(not_bond)
          I->Bond[cur_bond].not_order |= cH_Double;
        else
          I->Bond[cur_bond].order |= cH_Double;
        sym = cSym_Bond;
        break;
      case '#':
        c++;
        if(not_bond)
          I->Bond[cur_bond].not_order |= cH_Triple;
        else
          I->Bond[cur_bond].order |= cH_Triple;
        sym = cSym_Bond;
        break;
      case '~':
        c++;
        if(not_bond) {
          I->Bond[cur_bond].not_order |= cH_AnyOrder;
          I->Bond[cur_bond].not_class |= cH_AnyClass;
        } else {
          I->Bond[cur_bond].order |= cH_AnyOrder;
          I->Bond[cur_bond].class |= cH_AnyClass;
        }
        sym = cSym_Bond;
        break;
      case '@':
        c++;
        if(not_bond)
          I->Bond[cur_bond].not_cycle |= cH_Cyclic;
        else
          I->Bond[cur_bond].cycle |= cH_Cyclic;
        sym = cSym_Bond;
        break;
      case ':':
        c++;
        if(not_bond)
          I->Bond[cur_bond].not_class |= cH_Aromatic;
        else
          I->Bond[cur_bond].class |= cH_Aromatic;
        sym = cSym_Bond;
        break;
      case '.': /* separator */
        c++;
        sym = cSym_Separator;
        break;
      case '%':
        c++;
        if(c) { 
          mark_code = 10*((*c)-'0');
          c++;
        } /* else error */
        if(c) {
          sym = cSym_Mark;
          mark_code += (*c)-'0';
          c++;
        } /* else error */
        break;
      case '(':
        c++;
        sym = cSym_OpenScope;
        break;
      case ')':
        c++;
        sym = cSym_CloseScope;
        break;
      case '[':
        c++;
        sym = cSym_OpenBlock;
        break;
      case ']':
        c++;
        sym = cSym_CloseBlock;
        break;
      }
    }
    if(sym==cSym_Null) {
      PRINTFB(FB_smiles_parsing,FB_errors)
        " champ: error parsing smiles string at '%c' (char %ld) in\n champ: '%s'\n",*c,c-orig_c,orig_c        
        ENDFB;
      ok=false;
    }
    if(ok) {
      /* =========== actions based on root level parsing ========== */
      switch(sym) {
      case cSym_OpenBlock:
        ok = ChampParseAtomBlock(I,&c,cur_atom);
      case cSym_Atom:
        /* was there a preceeding atom? if so, then form bond and save atom */
        if(base_atom) {
          PRINTFD(FB_smiles_parsing) 
            " ChampSmiToPtr: saving atom %d\n",last_atom
            ENDFD;
          /* backward link */
          I->Bond[cur_bond].atom[0] = base_atom;
          I->Bond[cur_bond].atom[1] = cur_atom;
          I->Bond[cur_bond].pri[0] = lex_pri;
          I->Bond[cur_bond].pri[1] = lex_pri;
          if(!bond_flag) {
            if((I->Atom[cur_atom].class&cH_Aromatic)&&
               (I->Atom[base_atom].class&cH_Aromatic))
              I->Bond[cur_bond].order = (cH_Single|cH_Aromatic); /* is this right? */
            else
              I->Bond[cur_bond].order = cH_Single;
          } 
          I->Bond[cur_bond].tag = bond_tags; /* save bond tags */
          I->Bond[cur_bond].not_tag = bond_not_tags; /* save bond tags */
          bond_tags=0;
          bond_not_tags=0;
          ok = ChampAddBondToAtom(I,cur_atom,cur_bond);
          if(ok) {
            ok = ChampAddBondToAtom(I,base_atom,cur_bond);
            save_bond();
          }
          bond_flag=false;
          not_bond=false;
        } 
        base_atom = cur_atom;
        save_atom();
        break;
      case cSym_CloseBlock: /* should never be reached */
        break;
      case cSym_OpenScope: /* push base_atom onto stack */
        stack = ListElemPushInt(&I->Int,stack,base_atom);
        break;
      case cSym_CloseScope:
        if(!stack) {
          PRINTFB(FB_smiles_parsing,FB_errors)
            " champ: stack underflow for scope...\n"
            ENDFB;
          ok=false;
        } else {
          base_atom=I->Int[stack].value;
          stack = ListElemPop(I->Int,stack);
        }
        break;
      case cSym_Bond:
        bond_flag=true;
        break;
      case cSym_Mark:
        if(base_atom) {
          if(!mark[mark_code]) { /* opening cycle */
            mark[mark_code] = base_atom;
            mark_pri[mark_code] = lex_pri;
            bond_flag = false; /* ignore the first bond valence...we'll get it from the second half of the mark*/
            not_bond = false;
          } else { /* closing cycle */
            I->Bond[cur_bond].atom[0] = base_atom;
            I->Bond[cur_bond].atom[1] = mark[mark_code];
            I->Bond[cur_bond].pri[0] = lex_pri;
            I->Bond[cur_bond].pri[1] = mark_pri[mark_code];
            if(!bond_flag) {
              I->Bond[cur_bond].order = cH_Single;
            }
            ok = ChampAddBondToAtom(I,base_atom,cur_bond);
            if(ok) {
              ok = ChampAddBondToAtom(I,mark[mark_code],cur_bond);
              save_bond();
            }
            mark[mark_code]=0;
            bond_flag=false;
            not_bond=false;
          }
        } else {
          PRINTFB(FB_smiles_parsing,FB_errors)
            " champ:  syntax error...\n"
            ENDFB;
          ok = false;
        }
        break;
      case cSym_Separator:
        base_atom = 0;
        break;
      case cSym_Qualifier:
        break;
      }
    }
  }
  if(ok&&atom_list) {
    result = ListElemNewZero(&I->Pat);
    if(result) {
      I->ActivePatList = ListElemPushInt(&I->Int,I->ActivePatList,result);
      I->Pat[result].atom = atom_list;
      I->Pat[result].bond = bond_list;
    } else
      ok=false;
  }
  if(cur_atom) ChampAtomFree(I,cur_atom);
  if(cur_bond) ChampBondFree(I,cur_bond);
  if(result) ChampPatReindex(I,result);

  PRINTFD(FB_smiles_parsing) 
    " ChampSmiToPtr: returning pattern %d atom_list %d bond_list %d\n",result,atom_list,bond_list
    ENDFD;
  
  return(result);
}  

/* =============================================================== 
 * Matching 
 * =============================================================== */

int ChampPatIdentical(ListAtom *p,ListAtom *a)
{
  if(p->pos_flag!=a->pos_flag)
    return 0;
  else {
    if(p->pos_flag) 
      if((p->atom!=a->atom)||
         (p->charge!=a->charge)||
         (p->cycle!=a->cycle)||
         (p->class!=a->class)||
         (p->degree!=a->degree)||
         (p->valence!=a->valence))
        return 0;
  }
  if(p->neg_flag!=a->neg_flag)
    return 0;
  else {
    if(p->neg_flag) 
      if((p->not_atom!=a->atom)||
         (p->not_charge!=a->charge)||
         (p->not_cycle!=a->cycle)||
         (p->not_class!=a->class)||
         (p->not_degree!=a->degree)||
         (p->not_valence!=a->valence))
        return 0;
  }
  return 1;
}

int ChampAtomMatch(ListAtom *p,ListAtom *a)
{
  if((((!p->pos_flag)||
       (((!p->atom)||(p->atom&a->atom))&&
        ((!p->charge)||(p->charge&a->charge))&&
        ((!p->cycle)||(p->cycle&a->cycle))&&   
        ((!p->class)||(p->class&a->class))&&   
        ((!p->degree)||(p->degree&a->degree))&&
        ((!p->valence)||(p->valence&a->valence))))&&
      ((!p->neg_flag)||
       (((!p->not_atom)||(!(p->not_atom&a->atom)))&&
        ((!p->not_charge)||(!(p->not_charge&a->charge)))&&
        ((!p->not_cycle)||(!(p->not_cycle&a->cycle)))&&
        ((!p->not_class)||(!(p->not_class&a->class)))&&
        ((!p->not_degree)||(!(p->not_degree&a->degree)))&&
        ((!p->not_valence)||(!(p->not_valence&a->valence)))))))
    {
      if(p->name[0])
        if(strcmp(p->name,a->name))
          return 0;
      if(p->residue[0])
        if(strcmp(p->residue,a->residue))
          return 0;
      if(p->symbol[0])
        if(strcmp(p->symbol,a->symbol))
          return 0;
      if(p->hydro_flag) {
        if(p->tot_hydro>a->tot_hydro) { /* must have at least as many hydrogens as pattern... */
          return 0;
        }
      }
      return 1;
    }
  /* what about implicit hydrogens? */

  return 0;
  
}

int ChampBondMatch(ListBond *p,ListBond *a)
{
  return((((!p->order)||(p->order&a->order))&&
          ((!p->class)||(p->class&a->class))&&
          ((!p->cycle)||(p->cycle&a->cycle))&&   
          ((!p->not_order)||(!(p->not_order&a->order)))&&   
          ((!p->not_class)||(!(p->not_class&a->class)))&&
          ((!p->not_cycle)||(!(p->not_cycle&a->cycle)))));
}
