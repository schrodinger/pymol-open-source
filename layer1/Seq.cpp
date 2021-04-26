
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
#include "os_gl.h"
#include "os_std.h"

#include "Base.h"
#include "Seq.h"
#include "main.h"
#include "MemoryDebug.h"
#include "PyMOLObject.h"
#include "Scene.h"
#include "Text.h"

#include "Seeker.h"

#include "Menu.h"
#include "Executive.h"
#include "Ortho.h"
#include "CGO.h"

static int SeqFindRowCol(PyMOLGlobals * G, int x, int y, int *row_num_ptr,
                         int *col_num_ptr, int fixed_row)
{
  CSeq *I = G->Seq;
  int result = 0;
  int row_num = 0;
  int col_num = 0;

  if(I->ScrollBarActive) {
    y -= DIP2PIXEL(I->ScrollBarWidth);
  }
  if(fixed_row >= 0) {
    row_num = fixed_row;
  } else {
    row_num = (y - I->rect.bottom) / DIP2PIXEL(I->LineHeight);
    row_num = (I->NRow - 1) - row_num;
  }
  if((row_num >= 0) && (row_num < I->NRow)) {
    int char_num;
    CSeqRow *row;
    row = I->Row.data() + row_num;
    char_num = (x - I->rect.left - DIP2PIXEL(I->CharMargin)) / DIP2PIXEL(I->CharWidth);
    if(row->nCol && !row->label_flag)
      if(char_num < I->VisSize) {
        char_num += I->NSkip;
        if((char_num >= 0) && (char_num < row->ext_len) && (row->char2col)) {
          col_num = row->char2col[char_num];
          if(col_num) {
            col_num--;
            if(col_num < row->nCol) {
              result = true;
            } else if(fixed_row >= 0) {
              col_num = row->nCol - 1;
              result = true;
            }
          }
        } else if(char_num == 0) {
          col_num = 0;
          result = true;
        } else {
          col_num = row->nCol - 1;
          result = true;
        }
      }
  }
  if(result) {
    *row_num_ptr = row_num;
    *col_num_ptr = col_num;
  }
  return result;
}

void SeqUpdate(PyMOLGlobals * G)
{
  CSeq *I = G->Seq;

  if(I->Changed) {
    SeekerUpdate(G);
    I->Changed = false;
    I->Dirty = true;
    OrthoReshape(G, -1, -1, false);     /* careful, this is recursive... */
  }
  if(I->Dirty) {
    I->Handler->refresh(G, I->Row);
    I->Dirty = false;
  }
}

void CSeq::reshape(int width, int height)
{
  PyMOLGlobals *G = m_G;
  CSeq *I = G->Seq;
  Block::reshape(width, height);

  {                             /* get current sequence sizes */
    int a;
    I->Size = 0;
    for(a = 0; a < I->NRow; a++) {
      if(I->Row[a].ext_len > I->Size)
        I->Size = I->Row[a].ext_len;
    }
  }

  {
    int extra;
    I->VisSize = (I->rect.right - I->rect.left - 1) / DIP2PIXEL(I->CharWidth);
    /*    printf("%d %d %d %d %d\n",cw,I->rect.right,I->rect.left,I->VisSize,I->Size); */

    if(I->VisSize < 1)
      I->VisSize = 1;
    extra = I->Size - I->VisSize;
    if(extra <= 0) {
      I->ScrollBarActive = false;
    } else {
      I->ScrollBarActive = true;
      m_ScrollBar.setLimits(I->Size, I->VisSize);
    }
  }
}

void SeqDirty(PyMOLGlobals * G)
{
  CSeq *I = G->Seq;
  I->Dirty = true;
  SceneInvalidate(G);
}

void SeqChanged(PyMOLGlobals * G)
{
  CSeq *I = G->Seq;
  I->Changed = true;
  SceneInvalidate(G);
}

int CSeq::drag(int x, int y, int mod)
{
  PyMOLGlobals *G = m_G;
  CSeq *I = G->Seq;
  int pass = 0;
  int row_num;
  int col_num;
  if(!pass) {
    if(SeqFindRowCol(G, x, y, &row_num, &col_num, I->LastRow)) {
      if(I->Handler)
        I->Handler->drag(G, I->Row, row_num, col_num, mod);
      OrthoDirty(G);
    }
  }
  return (1);
}

int CSeq::release(int button, int x, int y, int mod)
{
  PyMOLGlobals *G = m_G;
  CSeq *I = G->Seq;
  int pass = 0;
  if(!pass) {
    int row_num;
    int col_num;
    if(SeqFindRowCol(G, x, y, &row_num, &col_num, I->LastRow)) {
      if(I->Handler)
        I->Handler->release(G, I->Row, button, row_num, col_num, mod);
      OrthoDirty(G);
    } else {
      if(I->Handler)
        I->Handler->release(G, I->Row, button, -1, -1, mod);
      OrthoDirty(G);
    }
  }
  I->DragFlag = false;
  I->LastRow = -1;
  return (1);
}

int SeqGetHeight(PyMOLGlobals * G)
{
  CSeq *I = G->Seq;
  int height = 0;

  if(I->NRow) {
    height = DIP2PIXEL(I->LineHeight * I->NRow + 4);
    if(I->ScrollBarActive)
      height += DIP2PIXEL(I->ScrollBarWidth);
  }
  return (height);
}

void SeqSetHandler(PyMOLGlobals * G, CSeqHandler * handler)
{
  CSeq *I = G->Seq;
  I->Handler = handler;
}

int CSeq::click(int button, int x, int y, int mod)
{
  PyMOLGlobals *G = m_G;
  CSeq *I = G->Seq;
  int pass = 0;
  int row_num;
  int col_num;

  switch(button) {
    case P_GLUT_BUTTON_SCROLL_FORWARD:
      I->m_ScrollBar.moveBy(-1);
      return 1;
    case P_GLUT_BUTTON_SCROLL_BACKWARD:
      I->m_ScrollBar.moveBy(1);
      return 1;
  }

  if(I->ScrollBarActive) {
    if((y - rect.bottom) < DIP2PIXEL(I->ScrollBarWidth)) {
      pass = 1;
      I->m_ScrollBar.click(button, x, y, mod);
    }
  }
  if(!pass) {
    if(SeqFindRowCol(G, x, y, &row_num, &col_num, -1)) {
      if(I->Handler)
        I->Handler->click(G, I->Row, button, row_num, col_num, mod, x, y);
      I->DragFlag = true;
      I->LastRow = row_num;
      OrthoDirty(G);
    } else {
      switch (button) {
      case P_GLUT_RIGHT_BUTTON:
        {
          ObjectNameType name;
          if(ExecutiveGetActiveSeleName(G, name, false, false)) {
            MenuActivate2Arg(G, x, y + DIP2PIXEL(20), x, y, false, "pick_sele", name, name);
          }
        }
        break;
      case P_GLUT_LEFT_BUTTON:
        if(I->Handler)
          I->Handler->click(G, I->Row, button, -1, -1, mod, x, y);
        break;
      }
    }
  }
  return (1);
}

void CSeq::draw(CGO* orthoCGO)
{
  PyMOLGlobals *G = m_G;
  CSeq *I = G->Seq;

  if(G->HaveGUI && G->ValidContext) {

    int x = rect.left;
    int y = rect.bottom + DIP2PIXEL(I->ScrollBarMargin) + 1;
    float bg_color[3] = { 0.f, 0.f, 0.f }, overlay_color[3] = { 1.0F, 1.0F, 1.0F };
    int label_color_index = SettingGetGlobal_color(G, cSetting_seq_view_label_color);
    const float *label_color = ColorGet(G, label_color_index);
    copy3f(label_color, overlay_color);

    if (SettingGetGlobal_b(G, cSetting_bg_gradient)){
      if(SettingGetGlobal_b(G, cSetting_seq_view_location)) {
	copy3f(ColorGet(G, SettingGetGlobal_color(G, cSetting_bg_rgb_bottom)), bg_color);
      } else {
	copy3f(ColorGet(G, SettingGetGlobal_color(G, cSetting_bg_rgb_top)), bg_color);
      }      
    } else {
      copy3f(ColorGet(G, SettingGetGlobal_color(G, cSetting_bg_rgb)), bg_color);
    }
    if(!SettingGetGlobal_b(G, cSetting_seq_view_overlay)) {
      if (orthoCGO){
	CGOColorv(orthoCGO, bg_color);
      } else {
	glColor3fv(bg_color);
      }
      fill(orthoCGO);
    }
    if(I->ScrollBarActive) {
      I->m_ScrollBar.setBox(rect.bottom + DIP2PIXEL(I->ScrollBarWidth),
                            rect.left + DIP2PIXEL(I->ScrollBarMargin),
                            rect.bottom + DIP2PIXEL(2),
                            rect.right - DIP2PIXEL(I->ScrollBarMargin));
      I->m_ScrollBar.draw(orthoCGO);
      y += DIP2PIXEL(I->ScrollBarWidth);
      I->NSkip = static_cast<int>(I->m_ScrollBar.getValue());
    } else {
      I->NSkip = 0;
    }
    if(I->NRow) {               /* get current sequence sizes */
      int a, b;
      float black[3] = { 0, 0, 0 };
      float blue[3] = { 0.5, 0.5, 1.0 };
      const float *cur_color;
      CSeqRow *row;
      CSeqCol *col;
      int xx, yy, ch_wid, pix_wid, tot_len;
      int y1 = y;
      int max_len = 0;
      int n_real = 0;
      int vis_size = I->VisSize;
      int first_allowed;
      int max_title_width = 0;
      char fill_char;
      const char *fill_str = SettingGetGlobal_s(G, cSetting_seq_view_fill_char);
      int unaligned_color_index =
        SettingGetGlobal_color(G, cSetting_seq_view_unaligned_color);
      int fill_color_index = SettingGetGlobal_color(G, cSetting_seq_view_fill_color);
      int unaligned_mode = SettingGetGlobal_i(G, cSetting_seq_view_unaligned_mode);
      const float *unaligned_color;
      const float *fill_color;

      if(unaligned_color_index == -1) {
        switch (unaligned_mode) {
        case 3:
          unaligned_color_index = -1;
          break;
        default:
          unaligned_color_index = fill_color_index;
          break;
        }
      }
      fill_color = ColorGet(G, fill_color_index);
      if(unaligned_color_index < 0) {
        unaligned_color = NULL;
      } else {
        unaligned_color = ColorGet(G, unaligned_color_index);
      }

      if(fill_str && fill_str[0]) {
        fill_char = fill_str[0];
        if(fill_char == ' ')
          fill_char = 0;
      } else
        fill_char = 0;
      /* measure titles */

      for(a = I->NRow - 1; a >= 0; a--) {
        row = I->Row.data() + a;
        col = row->col.data();
        if(row->label_flag || row->column_label_flag) {
          row->title_width = col->offset + (col->stop - col->start);
          if(max_title_width < row->title_width)
            max_title_width = row->title_width;
        }
      }

      /* draw titles */

      cur_color = overlay_color;
      TextSetColor(G, cur_color);
      for(a = I->NRow - 1; a >= 0; a--) {
        row = I->Row.data() + a;
        yy = y1 - 2;
        col = row->col.data();
        if((row->label_flag || row->column_label_flag) && row->nCol) {
          row->title_width = col->offset + (col->stop - col->start);
          xx =
            x + DIP2PIXEL(I->CharMargin) + DIP2PIXEL(I->CharWidth) * (col->offset +
                                                (max_title_width - row->title_width));
          ch_wid = (col->stop - col->start);
          pix_wid = DIP2PIXEL(I->CharWidth * ch_wid);
          tot_len = col->offset + ch_wid - I->NSkip;
          if(tot_len <= vis_size) {
            TextDrawSubStrFast(G, row->txt, xx, y1, col->start, ch_wid ORTHOCGOARGVAR);
          }
        }
        y1 += DIP2PIXEL(I->LineHeight);
      }

      y1 = y;
      for(a = I->NRow - 1; a >= 0; a--) {
        row = I->Row.data() + a;
        cur_color = overlay_color;
	if (orthoCGO){
	  CGOColorv(orthoCGO, cur_color);
	} else {
	  glColor3fv(cur_color);
	}
        yy = y1 - 2;
        if(max_len < row->ext_len)
          max_len = row->ext_len;
        if(!row->label_flag)
          n_real++;

        if(row->label_flag) {
          first_allowed = I->NSkip + max_title_width;
        } else if(row->column_label_flag) {
          first_allowed = I->NSkip + max_title_width + 1;
        } else
          first_allowed = I->NSkip;

        for(b = 1; b < row->nCol; b++) {

          if(row->label_flag && (b > 1))
            first_allowed = I->NSkip + max_title_width + 1;

          col = row->col + b;
          if(col->offset >= first_allowed) {
            xx = x + DIP2PIXEL(I->CharMargin) + DIP2PIXEL(I->CharWidth) * (col->offset - I->NSkip);
            ch_wid = (col->stop - col->start);
            pix_wid = DIP2PIXEL(I->CharWidth * ch_wid);
            tot_len = col->offset + ch_wid - I->NSkip;
            if(tot_len <= vis_size) {
              if(row->label_flag) {
                TextSetColor(G, cur_color);
		if (orthoCGO){
		  CGOColorv(orthoCGO, cur_color);
		} else {
		  glColor3fv(cur_color);
		}
              } else if(col->unaligned && unaligned_color) {
                float tmp_color[3];
                const float *v = ColorGet(G, col->color);
                switch (unaligned_mode) {
                case 1:
                case 4:
                  average3f(v, bg_color, tmp_color);
                  TextSetColor(G, tmp_color);
		  if (orthoCGO){
		    CGOColorv(orthoCGO, tmp_color);
		  } else {
		    glColor3fv(tmp_color);
		  }
                  break;
                case 2:
                case 5:
                  average3f(v, unaligned_color, tmp_color);
                  TextSetColor(G, tmp_color);
		  if (orthoCGO){
		    CGOColorv(orthoCGO, tmp_color);
		  } else {
		    glColor3fv(tmp_color);
		  }
                  break;
                default:
                  TextSetColor(G, unaligned_color);
		  if (orthoCGO){
		    CGOColorv(orthoCGO, unaligned_color);
		  } else {
		    glColor3fv(unaligned_color);
		  }
                  break;
                }
              } else {
                const float *v = ColorGet(G, col->color);
                TextSetColor(G, v);
		if (orthoCGO){
		  CGOColorv(orthoCGO, v);
		} else {
		  glColor3fv(v);
		}
              }
              if(col->inverse) {
		if (orthoCGO){
		  CGOBegin(orthoCGO, GL_TRIANGLE_STRIP);
		  CGOVertex(orthoCGO, xx, yy, 0.f);
		  CGOVertex(orthoCGO, xx, yy + DIP2PIXEL(I->LineHeight) - 1, 0.f);
		  CGOVertex(orthoCGO, xx + pix_wid, yy, 0.f);
		  CGOVertex(orthoCGO, xx + pix_wid, yy + DIP2PIXEL(I->LineHeight) - 1, 0.f);
		  CGOEnd(orthoCGO);
		} else {
		  glBegin(GL_POLYGON);
		  glVertex2i(xx, yy);
		  glVertex2i(xx, yy + DIP2PIXEL(I->LineHeight) - 1);
		  glVertex2i(xx + pix_wid, yy + DIP2PIXEL(I->LineHeight) - 1);
		  glVertex2i(xx + pix_wid, yy);
		  glEnd();
		}
                TextSetColor(G, black);
              }
              TextDrawSubStrFast(G, row->txt, xx, y1, col->start, ch_wid ORTHOCGOARGVAR);
            }
          }
        }

        if(fill_char) {

          TextSetColor(G, fill_color);

          for(b = 0; b < row->nFill; b++) {

            col = row->fill + b;
            if(col->offset >= first_allowed) {
              xx = x + DIP2PIXEL(I->CharMargin) + DIP2PIXEL(I->CharWidth) * (col->offset - I->NSkip);
              ch_wid = (col->stop - col->start);
              pix_wid = DIP2PIXEL(I->CharWidth * ch_wid);
              tot_len = col->offset + ch_wid - I->NSkip;
              if(tot_len <= vis_size) {
                TextDrawCharRepeat(G, fill_char, xx, y1, col->start, ch_wid ORTHOCGOARGVAR);
              }
            }
          }
        }

        y1 += DIP2PIXEL(I->LineHeight);
      }

      if(I->Handler->box_active) {
        int box_row = I->Handler->box_row;
        if((box_row >= 0) && (box_row < I->NRow)) {
          int start_col = I->Handler->box_start_col;
          int stop_col = I->Handler->box_stop_col;
          if(start_col > stop_col) {
            int tmp = stop_col;
            stop_col = start_col;
            start_col = tmp;
          }
          row = I->Row.data() + box_row;
          if((start_col >= 0) && (start_col < row->nCol) &&
             (stop_col >= 0) && (stop_col < row->nCol)) {
            int xx2;
            CSeqCol *col2;
            col = row->col + start_col;
            col2 = row->col + stop_col;

            /* trim spacers (if any) */
            while(col->spacer && (start_col < stop_col)) {
              start_col++;
              col = row->col + start_col;
            }
            while(col2->spacer && (start_col < stop_col)) {
              stop_col--;
              col2 = row->col + stop_col;
            }

            yy = y + ((I->NRow - 1) - box_row) * DIP2PIXEL(I->LineHeight) - 2;
            xx = x + DIP2PIXEL(I->CharMargin) + DIP2PIXEL(I->CharWidth) * (col->offset - I->NSkip);
            xx2 =
              x + DIP2PIXEL(I->CharMargin) + DIP2PIXEL(I->CharWidth) * (col2->offset +
                                                  (col2->stop - col2->start) - I->NSkip);
	    if (orthoCGO){
	      CGOColorv(orthoCGO, overlay_color);
	      CGOLineAsTriangleStrips(orthoCGO, xx, yy, xx2, yy + DIP2PIXEL(I->LineHeight) - 2);
	      /* TODO: need to convert to triangles
	      CGOBegin(orthoCGO, GL_LINE_LOOP);
	      CGOVertex(orthoCGO, xx, yy);
	      CGOVertex(orthoCGO, xx, yy + I->LineHeight - 2);
	      CGOVertex(orthoCGO, xx2, yy + I->LineHeight - 2);
	      CGOVertex(orthoCGO, xx2, yy);
	      CGOEnd(orthoCGO);*/
	    } else {
	      glBegin(GL_LINE_LOOP);
	      glVertex2i(xx, yy);
	      glVertex2i(xx, yy + DIP2PIXEL(I->LineHeight) - 2);
	      glVertex2i(xx2, yy + DIP2PIXEL(I->LineHeight) - 2);
	      glVertex2i(xx2, yy);
	      glEnd();
	    }
	  }
	}
      }
      if(I->ScrollBarActive) {
        int real_count = n_real;
        int mode = 0;
        float width = (float) (rect.right - rect.left);
        float start = 0, stop;
        int right = 0;
        float bot, top, cent;
        float height = (float) DIP2PIXEL(I->ScrollBarWidth - I->ScrollBarMargin);
        int last_color = -1;
        cur_color = blue;
        for(a = 0; a < I->NRow; a++) {
          row = I->Row.data() + a;
          if(!row->label_flag) {
            top =
              rect.bottom + DIP2PIXEL(I->ScrollBarMargin) + (height * real_count) / n_real;
            real_count--;
            bot =
              rect.bottom + DIP2PIXEL(I->ScrollBarMargin) + (height * real_count) / n_real;
            mode = 0;
            for(b = 0; b < row->nCol; b++) {
              col = row->col + b;
              if(col->inverse && (!mode)) {
                start = (width * col->offset) / max_len;
                right = col->offset + (col->stop - col->start);
                mode = 1;
                last_color = col->color;
                if(row->label_flag)
                  cur_color = overlay_color;
                else
                  cur_color = ColorGet(G, col->color);  /* is this safe? should be for single-threading */
              } else if((!col->inverse) && (mode)) {
                if(b) {
                  stop =
                    (width * (col[-1].offset + col[-1].stop - col[-1].start)) / max_len;
                } else {
                  stop = (width * col->offset) / max_len;
                }
                if((stop - start) < 1.0F) {
                  cent = (stop + start) * 0.5F;
                  start = cent - 0.5F;
                  stop = cent + 0.5F;
                }

		if (orthoCGO){
		  CGOColorv(orthoCGO, cur_color);
		  CGOBegin(orthoCGO, GL_TRIANGLE_STRIP);
		  CGOVertex(orthoCGO, start, bot, 0.f);
		  CGOVertex(orthoCGO, start, top, 0.f);
		  CGOVertex(orthoCGO, stop, bot, 0.f);
		  CGOVertex(orthoCGO, stop, top, 0.f);
		  CGOEnd(orthoCGO);
		} else {
		  glColor3fv(cur_color);
		  glBegin(GL_POLYGON);
		  glVertex2f(start, bot);
		  glVertex2f(start, top);
		  glVertex2f(stop, top);
		  glVertex2f(stop, bot);
		  glEnd();
		}
                mode = 0;
              } else if(col->inverse && mode) {
                if(last_color != col->color) {
                  if(b) {
                    stop =
                      (width * (col[-1].offset + col[-1].stop - col[-1].start)) / max_len;
                  } else {
                    stop = (width * col->offset) / max_len;
                  }
                  if((stop - start) < 1.0F) {
                    cent = (stop + start) * 0.5F;
                    start = cent - 0.5F;
                    stop = cent + 0.5F;
                  }
		  if (orthoCGO){
		    CGOColorv(orthoCGO, cur_color);
		    CGOBegin(orthoCGO, GL_TRIANGLE_STRIP);
		    CGOVertex(orthoCGO, start, bot, 0.f);
		    CGOVertex(orthoCGO, start, top, 0.f);
		    CGOVertex(orthoCGO, stop, bot, 0.f);
		    CGOVertex(orthoCGO, stop, top, 0.f);
		    CGOEnd(orthoCGO);
		  } else {
		    glColor3fv(cur_color);
		    glBegin(GL_POLYGON);
		    glVertex2f(start, bot);
		    glVertex2f(start, top);
		    glVertex2f(stop, top);
		    glVertex2f(stop, bot);
		    glEnd();
		  }
                  start = (width * col->offset) / max_len;
                  last_color = col->color;
                  if(row->label_flag)
                    cur_color = overlay_color;
                  else
                    cur_color = ColorGet(G, col->color);        /* is this safe? should be for single-threading */
                }
                right = col->offset + (col->stop - col->start);

              }
            }

            if(mode) {
              stop = width * right / max_len;
              if((stop - start) < 1.0F) {
                cent = (stop + start) * 0.5F;
                start = cent - 0.5F;
                stop = cent + 0.5F;
              }
	      if (orthoCGO){
		CGOColorv(orthoCGO, cur_color);
		CGOBegin(orthoCGO, GL_TRIANGLE_STRIP);
		CGOVertex(orthoCGO, start, bot, 0.f);
		CGOVertex(orthoCGO, start, top, 0.f);
		CGOVertex(orthoCGO, stop, bot, 0.f);
		CGOVertex(orthoCGO, stop, top, 0.f);
		CGOEnd(orthoCGO);
	      } else {
		glColor3fv(cur_color);
		glBegin(GL_POLYGON);
		glVertex2f(start, bot);
		glVertex2f(start, top);
		glVertex2f(stop, top);
		glVertex2f(stop, bot);
		glEnd();
	      }
	    }
          }
        }

        I->m_ScrollBar.drawHandle(0.35F, orthoCGO);
      }
    }
  }
}

int SeqInit(PyMOLGlobals * G)
{
  CSeq *I = nullptr;
  if((I = (G->Seq = new CSeq(G)))) {

    I->active = true;
    I->TextColor[0] = 1.0;
    I->TextColor[1] = 0.75;
    I->TextColor[2] = 0.75;
    OrthoAttach(G, I, cOrthoTool);
    I->m_ScrollBar.setValue(0);
    return 1;
  } else
    return 0;
}

static void SeqPurgeRowVLA(PyMOLGlobals * G)
{
  G->Seq->Row.clear();
}

void SeqSetRow(PyMOLGlobals * G, std::vector<CSeqRow>&& row, int nRow)
{
  CSeq *I = G->Seq;
  I->Row = std::move(row);
  I->NRow = nRow;
}

void SeqFree(PyMOLGlobals * G)
{
  DeleteP(G->Seq);
}

Block *SeqGetBlock(PyMOLGlobals * G)
{
  CSeq *I = G->Seq;
  {
    return (I);
  }
}
