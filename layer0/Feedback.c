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
#include"MemoryDebug.h"
#include"Feedback.h"
#include"Ortho.h"

typedef struct {
  char *Stack;
  int Depth;
} CFeedback;

char *FeedbackMask;

CFeedback Feedbk;

void FeedbackInit(void)
{
  int a;
  
  CFeedback *I=&Feedbk;

  I->Stack=VLAMalloc(FB_Total,sizeof(char),5,0);
  I->Depth=0;
  FeedbackMask = I->Stack;

  for(a=0;a<FB_Total;a++) {
    FeedbackMask[a] = FB_Results | FB_Errors | FB_Warnings | FB_Actions | FB_Details;
	 /* | FB_Everything;*/
  }  
}

void FeedbackFree(void)
{
  CFeedback *I=&Feedbk;

  VLAFreeP(I->Stack);
}

/* below we'll presume that any standard feedback on the feedback
module itself will be effected at the Python level, since feedback
levels will be changed as a matter of course inside of PyMOL in order
to quietly perform complex actions.  */

void FeedbackPush(void)
{
  CFeedback *I=&Feedbk;
  int a;
  I->Depth++;
  VLACheck(I->Stack,char,(I->Depth+1)*FB_Total);
  FeedbackMask=I->Stack+(I->Depth*FB_Total);
  for(a=0;a<FB_Total;a++) {
    FeedbackMask[a] = FeedbackMask[a-FB_Total];
  }
  PRINTFD(FB_Feedback) " Feedback: push\n" ENDFD;
}

void FeedbackPop(void)
{
  CFeedback *I=&Feedbk;
  if(I->Depth) {
    I->Depth--;
    FeedbackMask=I->Stack+(I->Depth*FB_Total);
  }
  PRINTFD(FB_Feedback) " Feedback: pop\n" ENDFD;
}

void FeedbackSetMask(unsigned int sysmod,unsigned char mask)
{
  int a;
  if((sysmod>0)&&(sysmod<FB_Total)) {
    FeedbackMask[sysmod] = mask;
  } else if(!sysmod) {
    for(a=0;a<FB_Total;a++) {
      FeedbackMask[a] = mask;
    }
  }
  PRINTFD(FB_Feedback)
    " FeedbackSetMask: sysmod %d, mask 0x%02X\n",sysmod,mask 
    ENDFD;
}


void FeedbackDisable(unsigned int sysmod,unsigned char mask)
{
  int a;
  if((sysmod>0)&&(sysmod<FB_Total)) {
    FeedbackMask[sysmod] = FeedbackMask[sysmod] & (0xFF-mask);
  } else if(!sysmod) {
    for(a=0;a<FB_Total;a++) {
      FeedbackMask[a] = FeedbackMask[a] & (0xFF-mask);
    }
  }
  PRINTFD(FB_Feedback)
    " FeedbackDisable: sysmod %d, mask 0x%02X\n",sysmod,mask 
    ENDFD;

}

void FeedbackEnable(unsigned int sysmod,unsigned char mask)
{
  int a;
  if((sysmod>0)&&(sysmod<FB_Total)) {
    FeedbackMask[sysmod] = FeedbackMask[sysmod] | mask;
  } else if(!sysmod) {
    for(a=0;a<FB_Total;a++) {
      FeedbackMask[a] = FeedbackMask[a] | mask;
    }
  }
  PRINTFD(FB_Feedback)
    " FeedbackEnable: sysmod %d, mask 0x%02X\n",sysmod,mask 
    ENDFD;
  
}

void FeedbackAutoAdd(unsigned int sysmod,unsigned char mask,char *str)
{
  if(Feedback(sysmod,mask))
    OrthoAddOutput(str);
}

void FeedbackAdd(char *str)
{
  OrthoAddOutput(str);
}
