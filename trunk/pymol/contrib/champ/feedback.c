#include"os_std.h"
#include"const.h"
#include"vla.h"
#include"feedback.h"
typedef struct {
  char *Stack;
  int Depth;
} CFeedback;

char *FeedbackMask;

CFeedback Feedbk;

static int feedback_init = true;

void FeedbackInit(void)
{
  int a;
  CFeedback *I=&Feedbk;

  if(feedback_init) {
    feedback_init=false;

    
    vla_malloc(I->Stack,char,FB_total);
    I->Depth=0;
    FeedbackMask = I->Stack;
    
    for(a=0;a<FB_total;a++) {
      FeedbackMask[a] = FB_results | FB_errors | FB_warnings | FB_actions | FB_details;
      /*	  | FB_everything;*/
    }  
  }

}

void FeedbackFree(void)
{
  CFeedback *I=&Feedbk;

  vla_free(I->Stack);
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
  vla_check(I->Stack,char,(I->Depth+1)*FB_total);
  FeedbackMask=I->Stack+(I->Depth*FB_total);
  for(a=0;a<FB_total;a++) {
    FeedbackMask[a] = FeedbackMask[a-FB_total];
  }
  PRINTFD(FB_Feedback) " Feedback: push\n" ENDFD;
}

void FeedbackPop(void)
{
  CFeedback *I=&Feedbk;
  if(I->Depth) {
    I->Depth--;
    FeedbackMask=I->Stack+(I->Depth*FB_total);
  }
  PRINTFD(FB_Feedback) " Feedback: pop\n" ENDFD;
}

void FeedbackSetMask(unsigned int sysmod,unsigned char mask)
{
  int a;
  if((sysmod>0)&&(sysmod<FB_total)) {
    FeedbackMask[sysmod] = mask;
  } else if(!sysmod) {
    for(a=0;a<FB_total;a++) {
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
  if((sysmod>0)&&(sysmod<FB_total)) {
    FeedbackMask[sysmod] = FeedbackMask[sysmod] & (0xFF-mask);
  } else if(!sysmod) {
    for(a=0;a<FB_total;a++) {
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
  if((sysmod>0)&&(sysmod<FB_total)) {
    FeedbackMask[sysmod] = FeedbackMask[sysmod] | mask;
  } else if(!sysmod) {
    for(a=0;a<FB_total;a++) {
      FeedbackMask[a] = FeedbackMask[a] | mask;
    }
  }
  PRINTFD(FB_Feedback)
    " FeedbackEnable: sysmod %d, mask 0x%02X\n",sysmod,mask 
    ENDFD;
  
}


