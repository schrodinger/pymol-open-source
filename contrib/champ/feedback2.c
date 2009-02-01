#include"os_std.h"
#include"const.h"
#include"vla.h"
#include"feedback2.h"

typedef struct {
  char *Stack;
  int Depth;
} Cfeedback_;

char *feedback_Mask;

Cfeedback_ Feedbk;

static int feedback_init = true;

void feedback_Init(void)
{
  int a;
  Cfeedback_ *I=&Feedbk;

  if(feedback_init) {
    feedback_init=false;

    
    vla_malloc(I->Stack,char,FB_total);
    I->Depth=0;
    feedback_Mask = I->Stack;
    
    for(a=0;a<FB_total;a++) {
      feedback_Mask[a] = FB_results | FB_errors | FB_warnings | FB_actions | FB_details;
      /*	  | FB_everything;*/
    }  
  }

}

void feedback_Free(void)
{
  Cfeedback_ *I=&Feedbk;

  vla_free(I->Stack);
}

/* below we'll presume that any standard feedback on the feedback
module itself will be effected at the Python level, since feedback
levels will be changed as a matter of course inside of PyMOL in order
to quietly perform complex actions.  */

void feedback_Push(void)
{
  Cfeedback_ *I=&Feedbk;
  int a;
  I->Depth++;
  vla_check(I->Stack,char,(I->Depth+1)*FB_total);
  feedback_Mask=I->Stack+(I->Depth*FB_total);
  for(a=0;a<FB_total;a++) {
    feedback_Mask[a] = feedback_Mask[a-FB_total];
  }
  PRINTFD(FB_feedback_) " feedback: push\n" ENDFD;
}

void feedback_Pop(void)
{
  Cfeedback_ *I=&Feedbk;
  if(I->Depth) {
    I->Depth--;
    feedback_Mask=I->Stack+(I->Depth*FB_total);
  }
  PRINTFD(FB_feedback_) " feedback: pop\n" ENDFD;
}

void feedback_SetMask(unsigned int sysmod,unsigned char mask)
{
  int a;
  if((sysmod>0)&&(sysmod<FB_total)) {
    feedback_Mask[sysmod] = mask;
  } else if(!sysmod) {
    for(a=0;a<FB_total;a++) {
      feedback_Mask[a] = mask;
    }
  }
  PRINTFD(FB_feedback_)
    " feedbackSetMask: sysmod %d, mask 0x%02X\n",sysmod,mask 
    ENDFD;
}


void feedback_Disable(unsigned int sysmod,unsigned char mask)
{
  int a;
  if((sysmod>0)&&(sysmod<FB_total)) {
    feedback_Mask[sysmod] = feedback_Mask[sysmod] & (0xFF-mask);
  } else if(!sysmod) {
    for(a=0;a<FB_total;a++) {
      feedback_Mask[a] = feedback_Mask[a] & (0xFF-mask);
    }
  }
  PRINTFD(FB_feedback_)
    " feedbackDisable: sysmod %d, mask 0x%02X\n",sysmod,mask 
    ENDFD;

}

void feedback_Enable(unsigned int sysmod,unsigned char mask)
{
  int a;
  if((sysmod>0)&&(sysmod<FB_total)) {
    feedback_Mask[sysmod] = feedback_Mask[sysmod] | mask;
  } else if(!sysmod) {
    for(a=0;a<FB_total;a++) {
      feedback_Mask[a] = feedback_Mask[a] | mask;
    }
  }
  PRINTFD(FB_feedback_)
    " feedbackEnable: sysmod %d, mask 0x%02X\n",sysmod,mask 
    ENDFD;
  
}


