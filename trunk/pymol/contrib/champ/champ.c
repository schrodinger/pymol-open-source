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
#include"os_memory.h"
#include"const.h"
#include"err.h"

#include"feedback.h"

#include"mac.h"
#include"vla.h"
#include"champ.h"
#include"strblock.h"

#define _MATCHDEBUG
#define _MATCHDEBUG2

char *ChampParseAliphaticAtom(CChamp *I,char *c,int atom,int mask,int len);
char *ChampParseAromaticAtom(CChamp *I,char *c,int atom,int mask,int len);
char *ChampParseStringAtom(CChamp *I,char *c,int atom,int len);
char *ChampParseBlock(CChamp *I,char *c,int atom);

void ChampPatDump(CChamp *I,int index);
void ChampAtomDump(CChamp *I,int index);

int ChampAtomToString(CChamp *I,int index,char *buf);
int ChampBondToString(CChamp *I,int index,char *buf);
int ChampAddBondToAtom(CChamp *I,int atom_index,int bond_index);

int ChampMatch(CChamp *I,int template,int target,
                int unique_start,int n_wanted,int *match_start);
int ChampMatch2(CChamp *I,int template,int target,
                        int start_tmpl,int start_targ,
                        int n_wanted,int *match_start);

int ChampFindUniqueStart(CChamp *I,int template,int target,int *multiplicity);
int ChampUniqueListNew(CChamp *I,int atom, int unique_list);
void ChampUniqueListFree(CChamp *I,int unique_list);

void ChampMatchDump(CChamp *I,int match_idx);
void ChampMatchFree(CChamp *I,int match_idx);
void ChampMatchFreeChain(CChamp *I,int match_idx);
void ChampPatReindex(CChamp *I,int index);

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
          h->Pat[p].unique_atom = ChampUniqueListNew(h,h->Pat[p].atom,0);
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
          h->Pat[p].unique_atom = ChampUniqueListNew(h,h->Pat[p].atom,0);
          pp = p;
          printf("%s\n",buffer);
        } else if(p) {
          /* locate unique atoms */
          h->Pat[p].unique_atom = ChampUniqueListNew(h,h->Pat[p].atom,0);
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
          h->Pat[p].unique_atom = ChampUniqueListNew(h,h->Pat[p].atom,0);
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
          h->Pat[p].unique_atom = ChampUniqueListNew(h,h->Pat[p].atom,0);
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
int PTruthCallStr(PyObject *object,char *method,char *argument);
int PTruthCallStr(PyObject *object,char *method,char *argument)
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

int PConvPyObjectToInt(PyObject *object,int *value);
int PConvPyObjectToInt(PyObject *object,int *value)
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

void UtilCleanStr(char *s);
void UtilCleanStr(char *s) /*remove flanking white and all unprintables*/
{
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

int PConvPyObjectToStrMaxClean(PyObject *object,char *value,int ln);
int PConvPyObjectToStrMaxClean(PyObject *object,char *value,int ln)
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

int PConvPyObjectToStrMaxLen(PyObject *object,char *value,int ln);
int PConvPyObjectToStrMaxLen(PyObject *object,char *value,int ln)
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

int PConvAttrToStrMaxLen(PyObject *obj,char *attr,char *str,int ll);
int PConvAttrToStrMaxLen(PyObject *obj,char *attr,char *str,int ll)
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

int ChampModelToPat(CChamp *I,PyObject *model) 
{

  int nAtom,nBond;
  int a;
  int ok=true;
  int cur_atom,last_atom = 0;
  ListAtom *at;
  int cur_bond,last_bond = 0;
  ListBond *bd;
  int charge,order;
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
              ChampParseAliphaticAtom(I,c,cur_atom,cH_Cl,2);
              std_flag = true;
              break;
            default:
              ChampParseAliphaticAtom(I,c,cur_atom,cH_C,1);
              std_flag = true;
              break;
            }
            break;
          case 'H': /* nonstandard */
            ChampParseAliphaticAtom(I,c,cur_atom,cH_H,1);
            std_flag = true;
            break;
          case 'N':
            ChampParseAliphaticAtom(I,c,cur_atom,cH_N,1);
            std_flag = true;
            break;      
          case 'O':
            ChampParseAliphaticAtom(I,c,cur_atom,cH_O,1);
            std_flag = true;
            break;      
          case 'B':
            switch(*(c+1)) {
            case 'r':
              ChampParseAliphaticAtom(I,c,cur_atom,cH_Br,2);
              std_flag = true;
              break;
            default:
              ChampParseAliphaticAtom(I,c,cur_atom,cH_B,1);
              std_flag = true;
              break;
            }
            break;
          case 'P':
            ChampParseAliphaticAtom(I,c,cur_atom,cH_P,1);
            std_flag = true;
            break;      
          case 'S':
            ChampParseAliphaticAtom(I,c,cur_atom,cH_S,1);
            std_flag = true;
            break;      
          case 'F':
            ChampParseAliphaticAtom(I,c,cur_atom,cH_F,1);
            std_flag = true;
            break;      
          case 'I':
            ChampParseAliphaticAtom(I,c,cur_atom,cH_I,1);
            std_flag = true;
            break;      
          }
          if(!std_flag) {
            ChampParseAliphaticAtom(I,c,cur_atom,cH_Sym,strlen(c));
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
            bd->atom[0]=atom1;
            bd->atom[1]=atom2;
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
    molec = PyObject_GetAttrString(model,"molecule");
  } else {
    molec = NULL;
  }

  mac_free(atom_index);

  result = ListElemNewZero(&I->Pat);
  I->Pat[result].atom = cur_atom;
  I->Pat[result].bond = cur_bond;
  I->Pat[result].chempy_molecule = molec;
    
  if(result) ChampPatReindex(I,result);
  return(result);
}
/* =============================================================== 
 * Class-specific Memory Management 
 * =============================================================== */

CChamp *ChampNew(void) {
  CChamp *I;

  FeedbackInit(); /* only has side-effects the first time */

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
  /* hmm...what language does this remind you of? */

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
  if(!I->Pat[pattern].unique_atom)
    I->Pat[pattern].unique_atom = ChampUniqueListNew(I,I->Pat[pattern].atom,0);
  if(!I->Pat[target].unique_atom)
    I->Pat[target].unique_atom = ChampUniqueListNew(I,I->Pat[target].atom,0);
  return(ChampMatch(I,pattern,target,
                    ChampFindUniqueStart(I,pattern,target,NULL),
                    1,NULL));
}

int ChampMatch_1V1_Map(CChamp *I,int pattern,int target,int limit)
{
  int match_start = 0;

  if(!I->Pat[pattern].unique_atom)
    I->Pat[pattern].unique_atom = ChampUniqueListNew(I,I->Pat[pattern].atom,0);
  if(!I->Pat[target].unique_atom)
    I->Pat[target].unique_atom = ChampUniqueListNew(I,I->Pat[target].atom,0);
  return(ChampMatch(I,pattern,target,
                    ChampFindUniqueStart(I,pattern,target,NULL),
                    limit,&match_start));
}

int ChampMatch_1VN_N(CChamp *I,int pattern,int list)
{
  int target;
  int c = 0;
  if(!I->Pat[pattern].unique_atom)
    I->Pat[pattern].unique_atom = ChampUniqueListNew(I,I->Pat[pattern].atom,0);
  while(list) {
    target = I->Int[list].value;
    
    if(!I->Pat[target].unique_atom)
      I->Pat[target].unique_atom = ChampUniqueListNew(I,I->Pat[target].atom,0);
    
    if(ChampMatch(I,pattern,target,
                  ChampFindUniqueStart(I,pattern,target,NULL),
                  1,NULL))
      
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
      printf("(%2d)-",atom_idx);
      atom_idx = I->Int2[m_atom_idx].value[1];
      ChampAtomDump(I,atom_idx);
      printf("(%2d)\n",atom_idx);
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

int ChampMatch(CChamp *I,int template,int target,int unique_start,int n_wanted,int *match_start) 
{ /* returns whether or not substructure exists, but doesn't do alignment */
  int n_match = 0;
  int start_targ;
  int tmpl_atom,targ_atom;
  int rep_targ_atom;
  int unique_targ;
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
                                    (n_wanted-n_match),match_start);
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
                 int *match_start)

{ /* does the template covalent tree match the target? */
  
  /* basic algorithm is a multilayer mark-and-sweep traversal 
     through all atoms and bonds starting from the 
     two nucleating atoms, which are assumed to be equivalent */
  int n_match = 0;

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
          if(tmpl_at->mark_tmpl>1) { /* virtual atom */
            bond_idx = 0; /* not qualified for futher branching... */
          } else {
            bond_off = 0;  
            bond_idx = tmpl_at->bond[bond_off];
            while(bond_idx) {  /* iterate over all bonds */
              if(I->Bond[bond_idx].mark_tmpl) { /* skip over marked bonds...*/
                bond_off++;
                bond_idx = tmpl_at->bond[bond_off];
              } else
                break; /* found an open bond */
            }
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
          n_match++;

          if(n_wanted>0) 
            if(n_match>=n_wanted) done_flag=true;
          
          if(match_start) {
            (*match_start) = ListElemPush(&I->Match,*match_start);
            match_ent = I->Match + (*match_start);
            
            atom_start=0;
            bond_start=0;

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

            if((tmpl_at->mark_tmpl-targ_at->mark_targ)==1) /* is atom available? */
              if(ChampBondMatch(I->Bond+bond_idx, /* does bond match? */
                                 I->Bond+tmpl_ent->bond))
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
  return n_match;
}

/* =============================================================== 
 * Smiles
 * =============================================================== */
char *ChampParseAliphaticAtom(CChamp *I,char *c,int atom,int mask,int len) 
{
  ListAtom *at;
  at=I->Atom+atom;
  at->atom |= mask;
  at->pos_flag = true;
  PRINTFD(FB_smiles_parsing) 
    " ChampParseAliphaticAtom: called.\n"
    ENDFD;

  return c+len;
}

char *ChampParseAromaticAtom(CChamp *I,char *c,int atom,int mask,int len) 
{
  ListAtom *at;
  at=I->Atom+atom;
  at->atom |= mask;
  at->class |= cH_Aromatic;
  at->pos_flag = true;

  PRINTFD(FB_smiles_parsing) 
    " ChampParseAromaticAtom: called.\n"
    ENDFD;
  return c+len;
}

char *ChampParseStringAtom(CChamp *I,char *c,int atom,int len) 
{
  ListAtom *at;
  at=I->Atom+atom;
  at->atom |= cH_Sym;
  at->symbol[0]=c[0];
  at->symbol[1]=c[1];
  at->pos_flag = true;

  PRINTFD(FB_smiles_parsing) 
    " ChampParseStringAtom: called.\n"
    ENDFD;
  return c+len;
}

char *ChampParseBlock(CChamp *I,char *c,int atom) 
{
  return c;
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

char *ChampPatToSmiVLA(CChamp *I,int index,char *vla)
{
  char *result;
  int n_atom;
  int cur_scope = 0;
  int cur_atom = 0;
  int cur_bond = 0;
  int a,c,l;
  int mark[MAX_RING];
  int mark_bond[MAX_RING];
  int index1;
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
            " SmiToStrVLA: creating new scope %d vs old %d\n",index1,cur_scope
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
    printf(" atom %d %s 0x%08x bonds: ",cur_atom,buf,at->atom);
    for(a=0;a<MAX_BOND;a++) {
      if(!at->bond[a]) break;
      else printf("%d ",at->bond[a]);
    }
    printf("\n");
    cur_atom = I->Atom[cur_atom].link;
  }
  cur_bond = I->Pat[index].bond;
  while(cur_bond) {
    bd = I->Bond+cur_bond;
    printf(" bond %d 0x%01x atoms %d %d order 0x%01x\n",
           cur_bond,bd->order,bd->atom[0],bd->atom[1],bd->order);
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

  if(index) {
    at=I->Atom + index;

  /* the following are non-trivial */

    trivial = trivial && !(at->neg_flag);
    trivial = trivial && !(at->charge&(cH_Cation|cH_Dication|cH_Anion|cH_Dianion));
    trivial = trivial && !(at->cycle|at->valence|at->degree);
    trivial = trivial && !(at->class&cH_Aliphatic);
    trivial = trivial && !((at->atom!=cH_Any) &&
                           at->atom&( cH_Na | cH_K  | cH_Ca | cH_Mg | cH_Zn | cH_Fe | cH_Cu | cH_Se |
                                       cH_X | cH_Y | cH_R0 | cH_R1 | cH_R2 | cH_R3 | cH_R4 |
                                      cH_R5 | cH_R6 | cH_R7 | cH_R8 | cH_R9 | cH_Sym));


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

    trivial = true;
  
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
      strcpy(buf,"?");
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
  int a;
  
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
  
  while((*c)&&ok) {
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
      /* standard, implicit atoms
       * B(3), C(4), N(3,5), O(2), P(3,5), S(2,4,6), F(1), Cl(1), Br(1), I(1) */
      case 'C':
        switch(*(c+1)) {
        case 'l':
          c = ChampParseAliphaticAtom(I,c,cur_atom,cH_Cl,2);
          sym = cSym_Atom;
          break;
        default:
          c = ChampParseAliphaticAtom(I,c,cur_atom,cH_C,1);
          sym = cSym_Atom;
          PRINTFD(FB_smiles_parsing) 
            " parsed: %p\n",c
            ENDFD;
          break;
        }
        break;
      case 'H': /* nonstandard */
        c = ChampParseAliphaticAtom(I,c,cur_atom,cH_H,1);
        sym = cSym_Atom;
        break;
      case 'N':
        c = ChampParseAliphaticAtom(I,c,cur_atom,cH_N,1);
        sym = cSym_Atom;
        break;      
      case 'O':
        c = ChampParseAliphaticAtom(I,c,cur_atom,cH_O,1);
        sym = cSym_Atom;
        break;      
      case 'B':
        switch(*(c+1)) {
        case 'r':
          c = ChampParseAliphaticAtom(I,c,cur_atom,cH_Br,2);
          sym = cSym_Atom;
          break;
        default:
          c = ChampParseAliphaticAtom(I,c,cur_atom,cH_B,1);
          sym = cSym_Atom;
          break;
        }
        break;
      case 'P':
        c = ChampParseAliphaticAtom(I,c,cur_atom,cH_P,1);
        sym = cSym_Atom;
        break;      
      case 'S':
        c = ChampParseAliphaticAtom(I,c,cur_atom,cH_S,1);
        sym = cSym_Atom;
        break;      
      case 'F':
        c = ChampParseAliphaticAtom(I,c,cur_atom,cH_F,1);
        sym = cSym_Atom;
        break;      
      case 'I':
        c = ChampParseAliphaticAtom(I,c,cur_atom,cH_I,1);
        sym = cSym_Atom;
        break;      
        /* standard implicit aromatic atoms */
      case 'c':
        c = ChampParseAromaticAtom(I,c,cur_atom,cH_C,1);
        sym = cSym_Atom;
        break;
      case 'n':
        c = ChampParseAromaticAtom(I,c,cur_atom,cH_N,1);
        sym = cSym_Atom;
        break;
      case 'o':
        c = ChampParseAromaticAtom(I,c,cur_atom,cH_O,1);
        sym = cSym_Atom;
        break;
      case 's':
        c = ChampParseAromaticAtom(I,c,cur_atom,cH_S,1);
        sym = cSym_Atom;
        break;
      case '-':
        c++;
        I->Bond[cur_bond].order = cH_Single;
        sym = cSym_Bond;
        break;
      case '=':
        c++;
        I->Bond[cur_bond].order = cH_Double;
        sym = cSym_Bond;
        break;
      case '#':
        c++;
        I->Bond[cur_bond].order = cH_Triple;
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
        }
        if(c) {
          sym = cSym_Mark;
          mark_code += (*c)-'0';
          c++;
        }
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
        " champ: error parsing smiles string\n"
        ENDFB;
      ok=false;
    }
    if(ok) {
      /* =========== actions based on root level parsing ========== */
      switch(sym) {
      case cSym_Atom:
        /* was there a preceeding atom? if so, then form bond and save atom */
        if(base_atom) {
          PRINTFD(FB_smiles_parsing) 
            " ChampSmiToPtr: saving atom %d\n",last_atom
            ENDFD;
          /* backward link */
          I->Bond[cur_bond].atom[0] = base_atom;
          I->Bond[cur_bond].atom[1] = cur_atom;
          if(!bond_flag) {
            I->Bond[cur_bond].order = cH_Single;
          }
          ok = ChampAddBondToAtom(I,cur_atom,cur_bond);
          if(ok) {
            ok = ChampAddBondToAtom(I,base_atom,cur_bond);
            save_bond();
          }
          bond_flag=false;
        } 
        base_atom = cur_atom;
        save_atom();
        break;
      case cSym_OpenScope: /* push base_atom onto stack */
        stack = ListElemPush(&I->Int,stack);
        I->Int[stack].value = base_atom;
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
            mark[mark_code]=base_atom;
          } else { /* closing cycle */
            I->Bond[cur_bond].atom[0] = base_atom;
            I->Bond[cur_bond].atom[1] = mark[mark_code];
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
      }
    }
  }
  if(ok&&atom_list) {
    result = ListElemNewZero(&I->Pat);
    I->Pat[result].atom = atom_list;
    I->Pat[result].bond = bond_list;
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
      return 1;
    }
  return 0;
  
}

int ChampBondMatch(ListBond *p,ListBond *a)
{
  return((((!p->order)||(p->order&a->order))&&
          ((!p->class)||(p->class&a->class))&&
          ((!p->cycle)||(p->cycle&a->cycle))&&   
          ((!p->not_order)||(p->not_order&a->not_order))&&   
          ((!p->not_class)||(p->not_class&a->not_class))&&
          ((!p->not_cycle)||(p->not_cycle&a->not_cycle))));
}




