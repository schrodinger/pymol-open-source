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

#include "os_gl.h"  
#include "os_std.h"

#include "Base.h"
#include "Seq.h"
#include "main.h"
#include "ScrollBar.h"
#include "MemoryDebug.h"
#include "Grap.h"
#include "PyMOLObject.h"

#include "Seeker.h"

typedef struct {
  Block *Block;
  int DragFlag;
  int ScrollBarActive;
  int NSkip;
  struct CScrollBar *ScrollBar;
  CSeqRow *Row;
  int NRow;
  int Size;
  int VisSize; 
  int Dirty;
  int LineHeight;
  int CharWidth;
  int ScrollBarWidth;
  int ScrollBarMargin;
  int CharMargin;
  CSeqHandler *Handler; /* borrowed pointer */
} CSeq;


CSeq Seq;

void SeqUpdate(void)
{
  CSeq *I=&Seq;

  if(I->Dirty) {
    SeekerUpdate();
    #if 0
    CSeqRow *row;
    I->Row = VLACalloc(CSeqRow,10);
    row = I->Row;
    
    row->txt = VLAlloc(char,200);
    strcpy(row->txt,"1         11        21        31        41        51         61        71        81");
    row->len = strlen(row->txt);
    I->NRow++;
    row = I->Row+1;
    row->txt = VLAlloc(char,200);
    strcpy(row->txt,"ACHIPENCPHIAETYCENILAAYTHENRAEWSCGLLLPQERATSWNYCVDFKLHNMFGIPLYTEWWWEHGSACKLNEMIP");
    row->len = strlen(row->txt);
    I->NRow++;

    #endif

    I->Dirty = false;
    OrthoReshape(-1,-1); /* careful, this is recursive... */
  }
}

static void SeqReshape(Block *block,int width, int height)
{
  CSeq *I=&Seq;
  BlockReshape(block,width,height);

  { /* get current sequence sizes */
    int a;
    I->Size=0;
    for(a=0;a<I->NRow;a++) {
      if(I->Row[a].len>I->Size)
        I->Size = I->Row[a].len;
    }
  }

  {
    int extra;
    I->VisSize = (I->Block->rect.right - I->Block->rect.left)/I->CharWidth;
    /*    printf("%d %d %d %d %d\n",cw,I->Block->rect.right,I->Block->rect.left,I->VisSize,I->Size);*/
    
    if(I->VisSize<1) I->VisSize = 1;
    extra = I->Size - I->VisSize;
    if(extra<=0) {
      I->ScrollBarActive = false;
    } else {
      I->ScrollBarActive = true;
      ScrollBarSetLimits(I->ScrollBar,I->Size,I->VisSize);
    }
  }
}

void SeqDirty(void)
{
  CSeq *I=&Seq;
  I->Dirty = true;
}
static int SeqDrag(Block *block,int x,int y,int mod)
{
  /*  CSeq *I=&Seq;*/
  return(1);
}

static int SeqRelease(Block *block,int button,int x,int y,int mod)
{
  CSeq *I=&Seq;  
  int pass=0;
  if(I->ScrollBarActive) {
    if((y-I->Block->rect.bottom)<I->ScrollBarWidth) {
      pass = 1;
      ScrollBarDoRelease(I->ScrollBar,button,x,y,mod);
      OrthoUngrab();
    }
  }
  I->DragFlag=false;
  return(1);
}

int SeqGetHeight(void)
{
  CSeq *I=&Seq;
  int height;

  height = 13*I->NRow;
  if(I->ScrollBarActive)
    height+=I->ScrollBarWidth;
  return(height);
}

static int FindRowCol(int x,int y,int *row_num_ptr,int *col_num_ptr)
{
  CSeq *I=&Seq;
  int result =0;
  int row_num = 0;
  int col_num;

  if(I->ScrollBarActive) {
    y-=I->ScrollBarWidth;
  } 
  row_num = (y-I->Block->rect.bottom)/I->LineHeight;
  if(row_num<I->NRow) {
    int char_num;
    CSeqRow *row;
    row_num = (I->NRow-1)-row_num;
    row = I->Row+row_num;
    char_num = (x-I->Block->rect.left-I->CharMargin)/I->CharWidth;
    if(char_num<I->VisSize) {
      char_num+=I->NSkip;
      if(char_num<0) char_num=0;
      if((char_num<row->len)&&(row->char2col)) {
        col_num = row->char2col[char_num];
        if(col_num) {
          col_num--;
          if(col_num<row->nCol) {
            result = true;
          }
        }
      }
    }
  }
  if(result) {
    *row_num_ptr = row_num;
    *col_num_ptr = col_num;
  }
  return result;
}
void SeqSetHandler(CSeqHandler *handler)
{
  CSeq *I=&Seq;
  I->Handler = handler;
}

static int SeqClick(Block *block,int button,int x,int y,int mod)
{
  CSeq *I=&Seq;
  int pass = 0;
  int row_num;
  int col_num;
  if(I->ScrollBarActive) {
    if((y-I->Block->rect.bottom)<I->ScrollBarWidth) {
      pass = 1;
      ScrollBarDoClick(I->ScrollBar,button,x,y,mod);      
    }
  } 
  if(FindRowCol(x,y,&row_num,&col_num)) {
    CSeqRow *row;
    CSeqCol *col;
    row = I->Row+row_num;
    col = row->col+col_num;
    col->inverse=!col->inverse;
    if(I->Handler)
      if(I->Handler->fClick)
        I->Handler->fClick(I->Row,button,row_num,col_num,mod);
    OrthoDirty();
  }
  return(1);
}

static void SeqDraw(Block *block)
{
  CSeq *I=&Seq;

  if(PMGUI) {
    int x = I->Block->rect.left;
    int y = I->Block->rect.bottom+I->ScrollBarMargin+1;
    
    if(I->ScrollBarActive) {
      ScrollBarSetBox(I->ScrollBar,I->Block->rect.bottom+I->ScrollBarWidth,
                      I->Block->rect.left+I->ScrollBarMargin,
                      I->Block->rect.bottom+2,
                      I->Block->rect.right-I->ScrollBarMargin);
      ScrollBarDoDraw(I->ScrollBar);
      y+=I->ScrollBarWidth;
      I->NSkip = (int)ScrollBarGetValue(I->ScrollBar);
    } else {
      I->NSkip = 0;
    }

    { /* get current sequence sizes */
      int a,b;
      float white[3] = {1,1,1};
      float black[3] = {0,0,0};
      CSeqRow *row;
      CSeqCol *col;
      glColor3fv(white);	 
      int xx,yy,ch_wid,pix_wid,tot_len;
      for(a=I->NRow-1;a>=0;a--) {
        row = I->Row+a;
        yy=y-2;
        for(b=0;b<row->nCol;b++) {
          col = row->col+b;
          if(col->start>=I->NSkip) {
            xx=x+I->CharMargin+I->CharWidth*(col->start-I->NSkip);
            ch_wid = (col->stop-col->start);
            pix_wid = I->CharWidth * ch_wid;
            tot_len = col->start+ch_wid-I->NSkip;
            if(tot_len<=I->VisSize) {
              if(col->inverse) {
                glColor3fv(white);
                glBegin(GL_POLYGON);
                glVertex2i(xx-1,yy);
                glVertex2i(xx-1,yy+I->LineHeight-1);
                glVertex2i(xx+pix_wid-1,yy+I->LineHeight-1);
                glVertex2i(xx+pix_wid-1,yy);
                glEnd();
                glColor3fv(black);
                GrapDrawSubStrFast(I->Row[a].txt,xx,y,
                                   col->start,ch_wid);
                glColor3fv(white);
              } else {
                GrapDrawSubStrFast(I->Row[a].txt,xx,y,
                                   col->start,ch_wid);
              }
            }
          }
        }
        y+=13;
      }
    }
    
  }
}

void SeqInit(void)
{
  CSeq *I=&Seq;
  I->Block = OrthoNewBlock(NULL);
  I->Block->fClick = SeqClick;
  I->Block->fDraw    = SeqDraw;
  I->Block->fDrag = SeqDrag;
  I->Block->fRelease = SeqRelease;
  I->Block->fReshape = SeqReshape;
  I->Block->active = true;
  I->Block->TextColor[0]=1.0;
  I->Block->TextColor[1]=0.75;
  I->Block->TextColor[2]=0.75;
  OrthoAttach(I->Block,cOrthoTool);
  I->DragFlag = false;
  I->ScrollBarActive = true;
  I->ScrollBar=ScrollBarNew(true);
  ScrollBarSetValue(I->ScrollBar,0);
  I->Row = NULL;
  I->NRow = 0;
  I->Dirty = true;
  I->ScrollBarWidth =16;
  I->ScrollBarMargin =2;
  I->LineHeight = 13;
  I->CharMargin = 2;
  I->CharWidth = GrapMeasureStr(" ");
}

static void SeqPurgeRowVLA(void)
{
  CSeq *I=&Seq;
  if(I->Row)
    {
      int a;
      CSeqRow *row;
      for(a=0;a<I->NRow;a++) {
        row = I->Row+a;
        VLAFreeP(row->txt);
        VLAFreeP(row->col);
        VLAFreeP(row->char2col);
        VLAFreeP(row->atom_lists);
      }
      VLAFreeP(I->Row);
    }
}

void SeqSetRowVLA(CSeqRow *row,int nRow)
{
  CSeq *I=&Seq;
  SeqPurgeRowVLA();
  I->Row = row;
  I->NRow = nRow;
}

void SeqFree(void)
{
  CSeq *I=&Seq;

  SeqPurgeRowVLA();
  if(I->ScrollBar)
    ScrollBarFree(I->ScrollBar);
  OrthoFreeBlock(I->Block);
}

Block *SeqGetBlock(void)
{
  CSeq *I=&Seq;
  {return(I->Block);}
}
