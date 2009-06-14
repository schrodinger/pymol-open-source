#ifndef _H_OVreturns
#define _H_OVreturns

typedef struct {
  ov_word status;
} OVstatus;


/* a few very common result codes */

#define OVstatus_NO_EFFECT         2
#define OVstatus_YES               1
#define OVstatus_NO                0

#define OVstatus_SUCCESS           0
#define OVstatus_FAILURE          -1
#define OVstatus_NULL_PTR         -2
#define OVstatus_OUT_OF_MEMORY    -3
#define OVstatus_NOT_FOUND        -4
#define OVstatus_DUPLICATE        -5
#define OVstatus_MISMATCH         -6
#define OVstatus_INVALID_REF_CNT  -7

#define return_OVstatus_NO_EFFECT { OVstatus _result = { OVstatus_FAILURE }; return _result; }
#define return_OVstatus_YES { OVstatus _result = { OVstatus_YES }; return _result; }
#define return_OVstatus_NO  { OVstatus _result = { OVstatus_NO }; return _result; }
#define return_OVstatus_SUCCESS { OVstatus _result = { OVstatus_SUCCESS  }; return _result; }

#define return_OVstatus_FAILURE { OVstatus _result = { OVstatus_FAILURE }; return _result; }
#define return_OVstatus_NULL_PTR { OVstatus _result = { OVstatus_NULL_PTR  }; return _result; }
#define return_OVstatus_OUT_OF_MEMORY { OVstatus _result = { OVstatus_OUT_OF_MEMORY  }; return _result; }
#define return_OVstatus_NOT_FOUND { OVstatus _result = { OVstatus_NOT_FOUND  }; return _result; }
#define return_OVstatus_DUPLICATE { OVstatus _result = { OVstatus_DUPLICATE  }; return _result; }
#define return_OVstatus_MISMATCH { OVstatus _result = { OVstatus_MISMATCH  }; return _result; }
#define return_OVstatus_INVALID_REF_CNT { OVstatus _result = { OVstatus_MISMATCH  }; return _result; }

#define OVreturn_IS_ERROR(r) ((r).status<0)
#define OVreturn_IS_OK(r) ((r).status>=0)

typedef struct {
  ov_word status;
  ov_word word;
} OVreturn_word;

typedef struct {
  ov_word status;
  ov_size size;
} OVreturn_size;

#endif
