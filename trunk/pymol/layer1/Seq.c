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
#include "Scene.h"

#include "Seeker.h"

#include "Menu.h"
#include "Executive.h"

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
  int Changed;
  int Dirty;
  int LineHeight;
  int CharWidth;
  int ScrollBarWidth;
  int ScrollBarMargin;
  int CharMargin;
  int LastRow;
  CSeqHandler *Handler; /* borrowed pointer */
} CSeq;


CSeq Seq;

static int FindRowCol(int x,int y,int *row_num_ptr,int *col_num_ptr,int fixed_row)
{
  CSeq *I=&Seq;
  int result =0;
  int row_num = 0;
  int col_num = 0;

  if(I->ScrollBarActive) {
    y-=I->ScrollBarWidth;
  } 
  if(fixed_row>=0) {
    row_num = fixed_row;
  } else {
    row_num = (y-I->Block->rect.bottom)/I->LineHeight;
    row_num = (I->NRow-1)-row_num;
  }
  if((row_num>=0)&&(row_num<I->NRow)) {
    int char_num;
    CSeqRow *row;
    row = I->Row+row_num;
    char_num = (x-I->Block->rect.left-I->CharMargin)/I->CharWidth;
    if(char_num<I->VisSize) {
      char_num+=I->NSkip;
      if(char_num<0) char_num=0;
      if((char_num<row->ext_len)&&(row->char2col)) {
        col_num = row->char2col[char_num];
        if(col_num) {
          col_num--;
          if(col_num<row->nCol) {
            result = true;
          } else if(fixed_row>=0) {
            if(col_num>0) 
              col_num = row->nCol - 1;
            if(col_num<0)
              col_num = 0;
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

void SeqUpdate(void)
{
  CSeq *I=&Seq;

  if(I->Changed) {
    SeekerUpdate();
    I->Changed = false;
    I->Dirty = true;
    OrthoReshape(-1,-1); /* careful, this is recursive... */
  }
  if(I->Dirty) {
    if(I->Handler->fRefresh)
      I->Handler->fRefresh(I->Row);
    I->Dirty = false;
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
    I->VisSize = (I->Block->rect.right - I->Block->rect.left -1)/I->CharWidth;
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
  SceneDirty();
}

void SeqChanged(void)
{
  CSeq *I=&Seq;
  I->Changed = true;
}

static int SeqDrag(Block *block,int x,int y,int mod)
{
  CSeq *I=&Seq;
  int pass = 0;
  int row_num;
  int col_num;
  if(!pass) {
    if(FindRowCol(x,y,&row_num,&col_num,I->LastRow)) {
      CSeqRow *row;
      CSeqCol *col;
      row = I->Row+row_num;
      col = row->col+col_num;
      if(I->Handler)
        if(I->Handler->fDrag)
          I->Handler->fDrag(I->Row,row_num,col_num,mod);
      OrthoDirty();
    }
  }
  return(1);
}

static int SeqRelease(Block *block,int button,int x,int y,int mod)
{
  CSeq *I=&Seq;  
  int pass=0;
  /*
    if(I->ScrollBarActive) {
      if((y-I->Block->rect.bottom)<I->ScrollBarWidth) {
        pass = 1;
        ScrollBarDoRelease(I->ScrollBar,button,x,y,mod);
     OrthoUngrab();
      }
    } 
  */
  if(!pass) {
    int row_num;
    int col_num;
    if(FindRowCol(x,y,&row_num,&col_num,I->LastRow)) {
      CSeqRow *row;
      CSeqCol *col;
      row = I->Row+row_num;
      col = row->col+col_num;
      if(I->Handler)
        if(I->Handler->fRelease)
          I->Handler->fRelease(I->Row,button,row_num,col_num,mod);
      OrthoDirty();
    }
  }
  I->DragFlag=false;
  I->LastRow = -1;
  return(1);
}

int SeqGetHeight(void)
{
  CSeq *I=&Seq;
  int height = 0;

  if(I->NRow) {
    height = 13*I->NRow + 14;
    if(I->ScrollBarActive)
      height+=I->ScrollBarWidth;
  }
  return(height);
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
  if(!pass) {    
    if(FindRowCol(x,y,&row_num,&col_num,-1)) {
      CSeqRow *row;
      CSeqCol *col;
      row = I->Row+row_num;
      col = row->col+col_num;
      if(I->Handler)
        if(I->Handler->fClick)
          I->Handler->fClick(I->Row,button,row_num,col_num,mod,x,y);
      I->DragFlag=true;
      I->LastRow = row_num;
      OrthoDirty();
    } else {
      if(button == P_GLUT_RIGHT_BUTTON) {
        char name[ObjNameMax];

        if(ExecutiveGetActiveSeleName(name, false)) {
          MenuActivate2Arg(x,y+20,x,y,"pick_option",name,name);
        }
      }
    }
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
    if(I->NRow) { /* get current sequence sizes */
      int a,b;
      float white[3] = {1,1,1};
      float black[3] = {0,0,0};
      float blue[3] = {0.5,0.5,1.0};
      float *cur_color;
      CSeqRow *row;
      CSeqCol *col;
      int xx,yy,ch_wid,pix_wid,tot_len;
      int y1=y;
      int max_len = 0;
      int n_real = 0;
      for(a=I->NRow-1;a>=0;a--) {
        row = I->Row+a;
        if(row->label_flag)
          cur_color = white;
        else
          cur_color = blue;
        glColor3fv(cur_color);
        yy=y1-2;
        if(max_len<row->ext_len)
          max_len = row->ext_len;
        if(!row->label_flag)
          n_real++;
        for(b=0;b<row->nCol;b++) {
          col = row->col+b;
          if((col->offset+(col->stop-col->start))>=I->NSkip) {
            xx=x+I->CharMargin+I->CharWidth*(col->offset-I->NSkip);
            ch_wid = (col->stop-col->start);
            pix_wid = I->CharWidth * ch_wid;
            tot_len = col->offset+ch_wid-I->NSkip;
            if(tot_len<=I->VisSize) {
              if(col->inverse) {
                glColor3fv(cur_color);
                glBegin(GL_POLYGON);
                glVertex2i(xx,yy);
                glVertex2i(xx,yy+I->LineHeight-1);
                glVertex2i(xx+pix_wid,yy+I->LineHeight-1);
                glVertex2i(xx+pix_wid,yy);
                glEnd();
                glColor3fv(black);
                GrapDrawSubStrFast(row->txt,xx,y1,
                                   col->start,ch_wid);
                glColor3fv(cur_color);
              } else {
                GrapDrawSubStrFast(row->txt,xx,y1,
                                   col->start,ch_wid);
              }
            }
          }
        }
        y1+=I->LineHeight;
      }
      if(I->Handler->box_active) {
        int box_row = I->Handler->box_row;
        if((box_row>=0)&&(box_row<I->NRow)) {
          int start_col = I->Handler->box_start_col;
          int stop_col = I->Handler->box_stop_col;
          if(start_col>stop_col) {
            register int tmp = stop_col;
            stop_col=start_col;
            start_col=tmp;
          }
          row = I->Row + box_row;
          if((start_col>=0)&&(start_col<row->nCol)&&
             (stop_col>=0)&&(stop_col<row->nCol)) {
            int xx2;
            CSeqCol *col2;
            col = row->col + start_col;
            col2 = row->col + stop_col;
            
            yy=y+((I->NRow-1)-box_row)*I->LineHeight-2;
            xx=x+I->CharMargin+I->CharWidth*(col->offset-I->NSkip);
            xx2=x+I->CharMargin+I->CharWidth*(col2->offset+(col2->stop-col2->start)-I->NSkip);
            
            glBegin(GL_LINE_LOOP);
            glVertex2i(xx,yy);
            glVertex2i(xx,yy+I->LineHeight-2);
            glVertex2i(xx2,yy+I->LineHeight-2);
            glVertex2i(xx2,yy);
            glEnd();
            
          }
        }
      }
      if(I->ScrollBarActive)
        {
          int real_count = n_real;
          int mode = 0;
          float width = I->Block->rect.right - I->Block->rect.left;
          float start,stop;
          float bot,top;
          float height = (I->ScrollBarWidth - I->ScrollBarMargin);
          cur_color = blue;
          for(a=0;a<I->NRow;a++) {
            row = I->Row+a;
            if(!row->label_flag) {
              top =  I->Block->rect.bottom + I->ScrollBarMargin + (height*real_count)/n_real; 
              real_count--;
              bot = I->Block->rect.bottom + I->ScrollBarMargin + (height*real_count)/n_real;
              mode = 0;
              for(b=0;b<row->nCol;b++) {
                col = row->col+b;
                if(col->inverse&&(!mode)) {
                  start = (width*col->offset)/max_len;
                  mode=1;
                } else if((!col->inverse)&&(mode)) {
                  stop = (width*col->offset)/max_len;
                  
                  glColor3fv(cur_color);
                  glBegin(GL_POLYGON);
                  glVertex2f(start,bot);
                  glVertex2f(start,top);
                  glVertex2f(stop,top);
                  glVertex2f(stop,bot);
                  glEnd();
                  mode = 0;
                }
              }
              
              if(mode) {
                col = row->col+row->nCol-1;
                stop = width*(col->offset+(col->stop-col->start)*I->CharWidth);
                glColor3fv(cur_color);
                glBegin(GL_POLYGON);
                glVertex2f(start,bot);
                glVertex2f(start,top);
                glVertex2f(stop,top);
                glVertex2f(stop,bot);
                glEnd();
              }
              
            }
          }
          
          ScrollBarDrawHandle(I->ScrollBar,0.35F);
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
  I->LastRow=-1;
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
