
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

#include"os_predef.h"
#include"os_gl.h"
#include"os_std.h"

#include"main.h"
#include"PyMOLObject.h"
#include"Color.h"
#include"Ortho.h"
#include"Scene.h"
#include"Util.h"
#include"Ray.h"
#include"PConv.h"
#include"Matrix.h"
#include"MemoryDebug.h"
#include"Movie.h"
#include"View.h"
#include"Err.h"

#include"Executive.h"
#include"CGO.h"
#include"Selector.h"
#include"vla.h"
#include"pymol/type_traits.h"

void ObjectPurgeSettings(pymol::CObject * I)
{
  I->Setting.reset();
}

void ObjectMotionTrim(pymol::CObject *I, int n_frame)
{
  if(I->ViewElem) {
    VLASize(I->ViewElem,CViewElem,n_frame);
  }
}

int ObjectMotionGetLength(pymol::CObject *I)
{
  if(I->ViewElem) {
    return VLAGetSize(I->ViewElem);
  }
  return 0;
}

void ObjectMotionReinterpolate(pymol::CObject *I)
{
  float power  = SettingGet_f(I->G, NULL, I->Setting.get(), cSetting_motion_power);
  float bias   = SettingGet_f(I->G, NULL, I->Setting.get(), cSetting_motion_bias);
  int simple   = SettingGet_i(I->G, NULL, I->Setting.get(), cSetting_motion_simple);
  float linear = SettingGet_f(I->G, NULL, I->Setting.get(), cSetting_motion_linear);
  int hand     = SettingGet_i(I->G, NULL, I->Setting.get(), cSetting_motion_hand);

  /*
     int ObjectMotion(pymol::CObject * I, int action, int first,
                 int last, float power, float bias,
                 int simple, float linear, int wrap,
                 int hand, int window, int cycles, int state, int quiet);
  */
  ObjectMotion(I, 3, -1, -1, power, bias, simple, linear,
               SettingGetGlobal_b(I->G,cSetting_movie_loop) ? 1 : 0,
               hand, 5, 1, -1, 1);
}

int ObjectMotionModify(pymol::CObject *I,int action, int index, int count,int target,int freeze,int localize)
{
  int ok;

  if(I->type == cObjectGroup) { /* propagate */
    ok = ExecutiveGroupMotionModify(I->G,I,action,index,count,target,freeze);
  } else {
    ok = ViewElemModify(I->G, &I->ViewElem,action,index,count,target);
    if(ok && I->ViewElem) {
      int size = VLAGetSize(I->ViewElem);
      int n_frame = MovieGetLength(I->G);
      if(n_frame != size) { 
        /* extend entire movie */
        if(!localize)
          ExecutiveMotionExtend(I->G,true);
        if((!freeze) && SettingGetGlobal_i(I->G,cSetting_movie_auto_interpolate)) {
          ExecutiveMotionReinterpolate(I->G);
        }
      } else if((!freeze) && SettingGetGlobal_i(I->G,cSetting_movie_auto_interpolate)) {
        ObjectMotionReinterpolate(I);
      }
    }
  }
  return ok;
}

static void TTTToViewElem(float *TTT, CViewElem * elem)
{
  float *fp = TTT;
  double *dp;

  /* convert row-major TTT to column-major ViewElem */

  elem->matrix_flag = true;
  dp = elem->matrix;

  dp[0] = (double) fp[0];
  dp[1] = (double) fp[4];
  dp[2] = (double) fp[8];
  dp[3] = 0.0;

  dp[4] = (double) fp[1];
  dp[5] = (double) fp[5];
  dp[6] = (double) fp[9];
  dp[7] = 0.0;

  dp[8] = (double) fp[2];
  dp[9] = (double) fp[6];
  dp[10] = (double) fp[10];
  dp[11] = 0.0;

  dp[12] = 0.0;
  dp[13] = 0.0;
  dp[14] = 0.0;
  dp[15] = 1.0;

  /* copy inverse pre */

  elem->pre_flag = true;
  dp = elem->pre;
  *(dp++) = (double) -TTT[12];
  *(dp++) = (double) -TTT[13];
  *(dp++) = (double) -TTT[14];

  /* copy post */

  elem->post_flag = true;
  dp = elem->post;
  *(dp++) = (double) TTT[3];
  *(dp++) = (double) TTT[7];
  *(dp++) = (double) TTT[11];

}

static void TTTFromViewElem(float *TTT, CViewElem * elem)
{
  float *fp = TTT;
  double *dp;

  if(elem->matrix_flag) {
    dp = elem->matrix;

    fp[0] = (float) dp[0];
    fp[1] = (float) dp[4];
    fp[2] = (float) dp[8];
    fp[3] = 0.0;

    fp[4] = (float) dp[1];
    fp[5] = (float) dp[5];
    fp[6] = (float) dp[9];
    fp[7] = 0.0;

    fp[8] = (float) dp[2];
    fp[9] = (float) dp[6];
    fp[10] = (float) dp[10];
    fp[11] = 0.0;

    fp[12] = 0.0;
    fp[13] = 0.0;
    fp[14] = 0.0;
    fp[15] = 1.0;
  }

  if(elem->pre_flag) {
    dp = elem->pre;
    fp[12] = (float) (-*(dp++));
    fp[13] = (float) (-*(dp++));
    fp[14] = (float) (-*(dp++));
  }

  if(elem->post_flag) {
    dp = elem->post;
    fp[3] = (float) *(dp++);
    fp[7] = (float) *(dp++);
    fp[11] = (float) *(dp++);
  }
  fp[15] = 1.0F;
}

int ObjectGetSpecLevel(pymol::CObject * I, int frame)
{
  if(I->ViewElem) {
    int size = VLAGetSize(I->ViewElem);
    if(frame<0) {
      int max_level = 0;
      int i;
      for(i=0;i<size;i++) {
        if(max_level < I->ViewElem[i].specification_level)
          max_level = I->ViewElem[i].specification_level;
      }
      return max_level;
    }
    if((frame>=0) && (frame<size))
      return I->ViewElem[frame].specification_level;
    return 0;
  }
  return -1;
}

void ObjectDrawViewElem(pymol::CObject *I, BlockRect *rect,int frames ORTHOCGOARG)
{
  if(I->ViewElem) {
    ViewElemDraw(I->G,I->ViewElem,rect,frames,I->Name ORTHOCGOARGVAR);
  }
}

int ObjectMotion(pymol::CObject * I, int action, int first,
               int last, float power, float bias,
               int simple, float linear, int wrap,
               int hand, int window, int cycles, int state, int quiet)
{
  PyMOLGlobals *G = I->G;
  if(I->type == cObjectGroup) { /* propagate */
    return ExecutiveGroupMotion(G,I,action,first,last, power,bias,simple,linear,
                                wrap,hand,window,cycles,state,quiet);
  } else {
    
    int frame;
    int nFrame = MovieGetLength(I->G);

    if(wrap<0) {
      wrap = SettingGet_b(I->G,NULL, I->Setting.get(), cSetting_movie_loop);
    }

    if(nFrame < 0)
      nFrame = -nFrame;

    if(!I->ViewElem) {
      I->ViewElem = pymol::vla<CViewElem>(0);
    }
    
    if((action == 7) || (action == 8)) { /* toggle */
      frame = first;
      if(first < 0)
        frame = SceneGetFrame(G);
      VLACheck(I->ViewElem, CViewElem, frame);
      if(action == 7) {
        if(I->ViewElem[frame].specification_level>1) {
          action = 1;
        } else {
          action = 0;
        }
      } else if(action == 8) {
        if(I->ViewElem[frame].specification_level>1) {
          int frame;
          action = 3;
          for(frame=0;frame<nFrame;frame++) {
            if(I->ViewElem[frame].specification_level==1) {
              action = 6;
              break;
            }
          }
        }
        else if(I->ViewElem[frame].specification_level>0) {
          action = 6;
        } else {
          action = 3;
        }
      }
    }

    if(action == 4) {   /* smooth */
      int save_last = last;
      if(first < 0)
        first = 0;
      
      if(last < 0) {
        last = nFrame;
      }
      if(last >= nFrame) {
        last = nFrame - 1;
      }
      if(first <= last) {
        int a;
        VLACheck(I->ViewElem, CViewElem, last);
          for(a = 0; a < cycles; a++) {
            ViewElemSmooth(I->ViewElem + first, I->ViewElem + last, window, wrap);
          }
      }
      if(SettingGet_b(I->G, NULL, I->Setting.get(), cSetting_movie_auto_interpolate)){
        action = 3; /* reinterpolate */
        last = save_last;
      }
    }
    switch (action) {
    case 0:                      /* store */
      if(!I->TTTFlag) {
        float mn[3], mx[3], orig[3];
        if(ExecutiveGetExtent(G, I->Name, mn, mx, true, -1, true)) {
          average3f(mn, mx, orig);
          ObjectSetTTTOrigin(I, orig);
        } else {
          initializeTTT44f(I->TTT);
          I->TTTFlag = true;
        }
      }
      if(I->ViewElem && I->TTTFlag) {
        if(first < 0)
          first = SceneGetFrame(G);
        if(last < 0)
          last = first;
        {
          int state_tmp=0, state_flag = false;
          if(state>=0) {
            state_tmp = state;
            state_flag = true;
          } else if(SettingGetIfDefined_i(G, I->Setting.get(), cSetting_state, &state_tmp)) {
            state_flag = true;
            state_tmp--;
          }
        
          for(frame = first; frame <= last; frame++) {
            if((frame >= 0) && (frame < nFrame)) {
              VLACheck(I->ViewElem, CViewElem, frame);
              if(!quiet) {
                PRINTFB(G, FB_Object, FB_Details)
                  " ObjectMotion: Setting frame %d.\n", frame + 1 ENDFB(G);
              }
              TTTToViewElem(I->TTT, I->ViewElem + frame);

              if(state_flag) {
                I->ViewElem[frame].state_flag = state_flag;
                I->ViewElem[frame].state = state_tmp;
              }

              if(power!=0.0F) {
                I->ViewElem[frame].power_flag = true;
                I->ViewElem[frame].power = power;
              }

              if(bias > 0.0F) {
                I->ViewElem[frame].bias_flag = true;
                I->ViewElem[frame].bias = bias;
              }

              I->ViewElem[frame].specification_level = 2;
            }

          }
        }
      }
      break;
    case 1:                      /* clear */
      if(I->ViewElem) {
        if(first < 0)
          first = SceneGetFrame(G);
        if(last < 0)
          last = first;
        for(frame = first; frame <= last; frame++) {
          if((frame >= 0) && (frame < nFrame)) {
            VLACheck(I->ViewElem, CViewElem, frame);
            ViewElemArrayPurge(G, I->ViewElem + frame, 1);
            UtilZeroMem((void *) (I->ViewElem + frame), sizeof(CViewElem));
          }
        }
      }
      break;
    case 2:                      /* interpolate & reinterpolate */
    case 3:
      {
        CViewElem *first_view = NULL, *last_view = NULL;
        int view_found = false;

        if(first < 0)
          first = 0;
        if(first > nFrame) {
          first = nFrame - 1;
        }

        if(last < 0) {
          last = nFrame;
          if(last) {
            if(!wrap)
              last--;
            else {
              int frame = 0;
              VLACheck(I->ViewElem, CViewElem, last);
              for(frame = 0; frame < last; frame++) {
                if(I->ViewElem[frame].specification_level > 1) {
                  last += frame;
                  break;
                }
              }
            }
          }
        } else {
          if(last >= nFrame) {
            last = nFrame;
            if(last && !wrap)
              last--;
          }
        }

        VLACheck(I->ViewElem, CViewElem, last);

        if(wrap && (last >= nFrame)) {
          /* if we're interpolating beyond the last frame, then wrap by
             copying early frames to last frames */
          int a;
          for(a = nFrame; a <= last; a++) {
            ViewElemCopy(G, I->ViewElem + a - nFrame, I->ViewElem + a);
          }
        } else if(!wrap) { 
          /* if we're not wrapping, then make sure we nuke any stray / old
             interpolated frames */
          frame = nFrame - 1;
          while(frame>=0) {
            if(I->ViewElem[frame].specification_level > 1) 
              break;
            else
              UtilZeroMem((void *) (I->ViewElem + frame), sizeof(CViewElem));
            frame--;
          }
        }
        VLACheck(I->ViewElem, CViewElem, last);
        if(!quiet) {
          if(action == 2) {
            if(last == nFrame) {
              PRINTFB(G, FB_Object, FB_Details)
                " ObjectMotion: interpolating unspecified frames %d to %d (wrapping).\n",
                first + 1, last ENDFB(G);
            } else {
              PRINTFB(G, FB_Object, FB_Details)
                " ObjectMotion: interpolating unspecified frames %d to %d.\n", first + 1,
                last + 1 ENDFB(G);
            }
          } else {
            if(last == nFrame) {
              PRINTFB(G, FB_Object, FB_Details)
                " ObjectMotion: reinterpolating all frames %d to %d (wrapping).\n", first + 1,
                last ENDFB(G);
            } else {
              PRINTFB(G, FB_Object, FB_Details)
                " ObjectMotion: reinterpolating all frames %d to %d.\n", first + 1, last + 1
                ENDFB(G);
            }
          }
        }
        for(frame = first; frame <= last; frame++) {
          if(!first_view) {
            if(I->ViewElem[frame].specification_level == 2) {     /* specified */
              first_view = I->ViewElem + frame;
              view_found = true;
            }
          } else {
            CViewElem *view;
            int interpolate_flag = false;
            if(I->ViewElem[frame].specification_level == 2) {     /* specified */
              last_view = I->ViewElem + frame;
              if(action == 2) {   /* interpolate */
                for(view = first_view + 1; view < last_view; view++) {
                  if(!view->specification_level)
                    interpolate_flag = true;
                }
              } else {
                interpolate_flag = true;
              }
              if(interpolate_flag) {
                ViewElemInterpolate(G, first_view, last_view,
                                    power, bias, simple, linear, hand, 0.0F);
              }
              first_view = last_view;
              last_view = NULL;
            }
          }
        }

        if(first_view) {
          if(wrap && (last >= nFrame)) {
            /* if we're interpolating beyond the last frame, then wrap by
               copying the last frames back over the early frames */
            int a;
            for(a = nFrame; a <= last; a++) {
              ViewElemCopy(G, I->ViewElem + a, I->ViewElem + a - nFrame);
            }
          }
        }

        if((!view_found) && (last>=first) && (first>=0) && (last<=nFrame)) {
          UtilZeroMem(I->ViewElem + first, sizeof(CViewElem) * (1 + (last-first)));
        }

        if(last >= nFrame) {   /* now erase temporary views */
          ViewElemArrayPurge(G, I->ViewElem + nFrame, (1 + last - nFrame));
          UtilZeroMem((void *) (I->ViewElem + nFrame),
                      sizeof(CViewElem) * (1 + last - nFrame));
        }
      }
      break;
    case 5:                      /* reset */
      if(I->ViewElem) {
        VLAFreeP(I->ViewElem);
      }
      I->ViewElem = pymol::vla<CViewElem>(0);
      break;
    case 6:                      /* uninterpolate */
      if(I->ViewElem) {
        if(first < 0)
          first = 0;
        if(last < 0) {
          last = nFrame - 1;
        }
        for(frame = first; frame <= last; frame++) {
          if((frame >= 0) && (frame <= last)) {
            VLACheck(I->ViewElem, CViewElem, frame);
            if(I->ViewElem[frame].specification_level < 2) {
              ViewElemArrayPurge(G, I->ViewElem + frame, 1);
              UtilZeroMem((void *) (I->ViewElem + frame), sizeof(CViewElem));
            }
          }
        }
      }
      break;
    case 9:
      if(I->ViewElem) {
        VLAFreeP(I->ViewElem);
      }
      break;
    }
    if(I->ViewElem) {
      VLASize(I->ViewElem,CViewElem,nFrame);
    }
  }
  return 1;
}

void ObjectAdjustStateRebuildRange(pymol::CObject * I, int *start, int *stop)
{
  /* on entry, start and stop should hold the valid range for the object */
  int defer_builds_mode =
    SettingGet_i(I->G, NULL, I->Setting.get(), cSetting_defer_builds_mode);
  int async_builds = SettingGet_b(I->G, NULL, I->Setting.get(), cSetting_async_builds);
  int max_threads = SettingGet_i(I->G, NULL, I->Setting.get(), cSetting_max_threads);
  int all_states = SettingGet_i(I->G, NULL, I->Setting.get(), cSetting_all_states);
  int dummy;
  if (all_states)
    return;
  if(defer_builds_mode >= 3) {
    if(SceneObjectIsActive(I->G, I))
      defer_builds_mode = 2;
  }
  switch (defer_builds_mode) {
  case 1:                      /* defer geometry builds until needed */
  case 2:                      /* defer and destroy continuously for increase memory conservation */
    if(SettingGetIfDefined_i(I->G, I->Setting.get(), cSetting_state, &dummy)) {
      /* decoupled...so always build all states.  Otherwise, geometry
      may not be there when we need it... unfortunately, this defeats
      the purpose of defer_builds_mode! */
    } else {
      int min = *start;
      int max = *stop;
      int global_state = SceneGetState(I->G);
      int obj_state = ObjectGetCurrentState(I, false);
      
      *start = obj_state;
      if((obj_state != global_state) || (!async_builds) || (max_threads < 1)) {
        *stop = *start + 1;
        if(*stop > max )
          *stop = max;
      } else {
        int base = (*start / max_threads);
        *start = (base) * max_threads;
        *stop = (base + 1) * max_threads;
        if(*start < min)
          *start = min;
        if(*start > max)
          *start = max;
        if(*stop < min)
          *stop = min;
        if(*stop > max)
          *stop = max;
      }
      if(*start > obj_state)
        *start = obj_state;
      if(*stop <= obj_state)
        *stop = obj_state + 1;
      if(*start < 0)
        *start = 0;
    }
    break;
  case 3:                      /* object not active, so do not rebuild anything */
    *stop = *start;
    break;
  }
}

/**
 * Replaces invalid characters in the given object name with an underscore,
 * or strips them if they are terminal or sequential.
 * @param[in,out] name Object name to validate
 * @return true if name was modified, false otherwise
 */
bool ObjectMakeValidName(char *name)
{
  bool modified = false;
  char *p = name, *q;
  if(p) {
    /* currently legal are A to Z, a to z, 0 to 9, -, _, + */
    while(*p) {
      switch (*p) {
        case '+':
        case '-':
        case '.':
        case '^':
        case '_':
          break;
        default:
          if (('A' <= *p && *p <= 'Z') ||
              ('a' <= *p && *p <= 'z') ||
              ('0' <= *p && *p <= '9'))
            break;
        /* must be an ASCII-visible character */
        *p = 1;                 /* placeholder for non-printable */
        modified = true;
      }
      p++;
    }
    /* eliminate sequential and terminal nonprintables */
    p = name;
    q = name;
    while(*p) {
      if(q == name)
        while(*p == 1)
          p++;
      while((*p == 1) && (p[1] == 1))
        p++;
      *q++ = *p++;
      if(!p[-1])
        break;
    }
    *q = 0;
    while(q > name) {
      if(q[-1] == 1) {
        q[-1] = 0;
        q--;
      } else
        break;
    }
    /* convert invalides to underscore */
    p = name;
    while(*p) {
      if(*p == 1)
        *p = '_';
      p++;
    }
  }
  return modified;
}

/**
 * Replaces invalid characters in `name` with an underscore,
 * or strips them if they are terminal or sequential - if `name` equals a reserved
 * selection keyword, then also append an underscore.
 *
 * @param[in,out] name Object name to validate
 * @param quiet If false, print warnings if the name gets modified.
 */
void ObjectMakeValidName(PyMOLGlobals * G, char *name, bool quiet)
{
  if (ObjectMakeValidName(name) && !quiet) {
    PRINTFB(G, FB_Executive, FB_Warnings)
      " Warning: Invalid characters in '%s' have been replaced or stripped\n",
      name ENDFB(G);
  }

  if (SelectorNameIsKeyword(G, name)) {
    if (!quiet) {
      PRINTFB(G, FB_Executive, FB_Warnings)
      " Warning: '%s' is a reserved keyword, appending underscore\n", name
      ENDFB(G);
    }
    strcat(name, "_");
    return;
  }

  static bool once_protein = false;
  static bool once_nucleic = false;

  if (!once_protein && strcmp(name, "protein") == 0) {
    once_protein = true;
  } else if (!once_nucleic && strcmp(name, "nucleic") == 0) {
    once_nucleic = true;
  } else {
    return;
  }

  {
    // Warn the user if "protein" or "nucleic" are used as names, but
    // don't modify the name (yet).
    PRINTFB(G, FB_Executive, FB_Warnings)
      " Warning: '%s' may become a reserved selection keyword in the future\n", name
      ENDFB(G);
  }
}

/**
 * Get a pointer to an object state.
 * @param state State (0-indexed) or -2/-3 for current state
 * @return NULL if state is out of bounds or empty
 */
CObjectState* pymol::CObject::getObjectState(int state)
{
  if (state == -2 /* cSelectorUpdateTableCurrentState */ ||
      state == -3 /* cSelectorUpdateTableEffectiveStates */) {
    state = getCurrentState();
  }
  if (state < 0 || state >= getNFrame()) {
    return nullptr;
  }
  return _getObjectState(state);
}

/**
 * Get the effective state (0-indexed) of an object, based on the `state` and
 * `static_singletons` settings. Will not validate the value of the `state`
 * setting, it could be `<0` or `>=getNFrame()`.
 */
int pymol::CObject::getCurrentState() const
{
  if (getNFrame() == 1 &&
      SettingGet<bool>(G, Setting.get(), nullptr, cSetting_static_singletons))
    return 0;
  return SettingGet<int>(G, Setting.get(), nullptr, cSetting_state) - 1;
}

/**
 * Like pymol::CObject::getCurrentState() but will return `-1` if the `all_states`
 * setting is set.
 *
 * Note: Clamps negative values at `-1` (all states). The usefulness of this
 * should be questioned, in particular with `ignore_all_states=true` a caller
 * is likely to discard all negative values, including -1.
 *
 * @param ignore_all_states Boolean flag, should be false. You most likely
 * should use pymol::CObject::getCurrentState() instead of setting `ignore_all_states`
 * to true.
 */
int ObjectGetCurrentState(const pymol::CObject * I, int ignore_all_states)
{
  assert("use pymol::CObject::getCurrentState()" && !ignore_all_states);

  // the previous implementation (up to PyMOL 1.7.6) ignored
  // object-level state=0 (all states)

  if (!ignore_all_states &&
      SettingGet_b(I->G, I->Setting.get(), NULL, cSetting_all_states))
    return -1;

  return std::max(-1, I->getCurrentState());
}

PyObject *ObjectAsPyList(pymol::CObject * I)
{
  PyObject *result = NULL;
  result = PyList_New(14);
  PyList_SetItem(result, 0, PyInt_FromLong(I->type));
  PyList_SetItem(result, 1, PyString_FromString(I->Name));
  PyList_SetItem(result, 2, PyInt_FromLong(I->Color));
  PyList_SetItem(result, 3, PyInt_FromLong(I->visRep));
  PyList_SetItem(result, 4, PConvFloatArrayToPyList(I->ExtentMin, 3));
  PyList_SetItem(result, 5, PConvFloatArrayToPyList(I->ExtentMax, 3));
  PyList_SetItem(result, 6, PyInt_FromLong(I->ExtentFlag));
  PyList_SetItem(result, 7, PyInt_FromLong(I->TTTFlag));
  PyList_SetItem(result, 8, SettingAsPyList(I->Setting.get()));

  PyList_SetItem(result, 9, PyInt_FromLong(I->Enabled));
  PyList_SetItem(result, 10, PyInt_FromLong(static_cast<long>(I->getRenderContext())));
  PyList_SetItem(result, 11, PConvFloatArrayToPyList(I->TTT, 16));
  if(I->ViewElem) {
    int nFrame = VLAGetSize(I->ViewElem);
    PyList_SetItem(result, 12, PyInt_FromLong(nFrame));
    PyList_SetItem(result, 13, ViewElemVLAAsPyList(I->G, I->ViewElem, nFrame));
  } else {
    PyList_SetItem(result, 12, PyInt_FromLong(0));
    PyList_SetItem(result, 13, PConvAutoNone(NULL));
  }
  return (PConvAutoNone(result));
}

int ObjectFromPyList(PyMOLGlobals * G, PyObject * list, pymol::CObject * I)
{
  int ok = true;
  int ll = 0;
  I->G = G;
  if(ok)
    ok = (list != NULL);
  if(ok)
    ok = PyList_Check(list);
  if(ok)
    ll = PyList_Size(list);
  if(ok)
    ok = CPythonVal_PConvPyIntToInt_From_List(G, list, 0, reinterpret_cast<int*>(&I->type));
  if(ok)
    ok = PConvPyStrToStr(PyList_GetItem(list, 1), I->Name, WordLength);
  if(ok)
    ok = PConvPyIntToInt(PyList_GetItem(list, 2), &I->Color);
  if(ok)
    I->Color = ColorConvertOldSessionIndex(G, I->Color);
  if(ok) {
    PyObject *val = PyList_GetItem(list, 3);
    if(PyList_Check(val)) {
      ok = PConvPyListToBitmask(val, &I->visRep, cRepCnt);
    } else {
      ok = PConvPyIntToInt(val, &I->visRep);
    }
    CPythonVal_Free(val);
  }
  if(ok)
    ok = PConvPyListToFloatArrayInPlaceAutoZero(PyList_GetItem(list, 4), I->ExtentMin, 3);
  if(ok)
    ok = PConvPyListToFloatArrayInPlaceAutoZero(PyList_GetItem(list, 5), I->ExtentMax, 3);
  if(ok)
    ok = PConvPyIntToInt(PyList_GetItem(list, 6), &I->ExtentFlag);
  if(ok)
    ok = PConvPyIntToInt(PyList_GetItem(list, 7), &I->TTTFlag);
  if(ok){
    CPythonVal *val = CPythonVal_PyList_GetItem(G, list, 8);
    I->Setting.reset(SettingNewFromPyList(G, val));
    CPythonVal_Free(val);
  }
  if(ok && (ll > 9))
    ok = CPythonVal_PConvPyIntToInt_From_List(G, list, 9, &I->Enabled);
  if(ok && (ll > 10)) {
    // I->Context removed.
  }
  if(ok && (ll > 11))
    ok = PConvPyListToFloatArrayInPlaceAutoZero(PyList_GetItem(list, 11), I->TTT, 16);
  if(ok && (ll > 13)) {
    PyObject *tmp;
    int nFrame;
    VLAFreeP(I->ViewElem);
    I->ViewElem = NULL;
    if(ok)
      ok = PConvPyIntToInt(PyList_GetItem(list, 12), &nFrame);
    if(ok && nFrame) {
      tmp = PyList_GetItem(list, 13);
      if(tmp && !(tmp == Py_None))
        ok = ViewElemVLAFromPyList(G, tmp, &I->ViewElem, nFrame);
    }
  }
  /* TO SUPPORT BACKWARDS COMPATIBILITY...
     Always check ll when adding new PyList_GetItem's */

  return (ok);
}

int ObjectCopyHeader(pymol::CObject * I, const pymol::CObject * src)
{
  int ok = true;

  I->G = src->G;
  I->type = src->type;
  UtilNCopy(I->Name, src->Name, WordLength);
  I->Color = src->Color;
  I->visRep = src->visRep;
  copy3f(src->ExtentMin, I->ExtentMin);
  copy3f(src->ExtentMax, I->ExtentMax);

  I->ExtentFlag = src->ExtentFlag;
  I->TTTFlag = src->TTTFlag;
  I->Setting = src->Setting;
  I->Enabled = src->Enabled;
  {
    int a;
    for(a = 0; a < 16; a++)
      I->TTT[a] = src->TTT[a];
  }
  I->ViewElem = NULL;           /* to do */

  return (ok);
}


/*========================================================================*/
void ObjectCombineTTT(pymol::CObject * I, const float *ttt, int reverse_order, int store)
{
  if(I->type == cObjectGroup) {
    ExecutiveGroupCombineTTT(I->G, I, ttt, reverse_order,store);
  } else {
    float cpy[16];
    if(!I->TTTFlag) {
      I->TTTFlag = true;
      initializeTTT44f(cpy);
    } else {
      UtilCopyMem(cpy, I->TTT, sizeof(float) * 16);
    }
    if(reverse_order) {
      combineTTT44f44f(cpy, ttt, I->TTT);
    } else {
      combineTTT44f44f(ttt, cpy, I->TTT);
    }
    if(store<0) 
      store = SettingGet_i(I->G, I->Setting.get(), NULL, cSetting_movie_auto_store);
    if(store && MovieDefined(I->G)) {
      if(!I->ViewElem)  
        I->ViewElem = pymol::vla<CViewElem>(0);
      if(I->ViewElem) { /* update motion path waypoint, if active */
        int frame = SceneGetFrame(I->G);
        if(frame >= 0) {
          VLACheck(I->ViewElem, CViewElem, frame);
          TTTToViewElem(I->TTT, I->ViewElem + frame);
          I->ViewElem[frame].specification_level = 2;
        }
      }
    }
  }
}
/*========================================================================*/
void ObjectTranslateTTT(pymol::CObject * I, const float *v, int store)
{
  if(I->type == cObjectGroup) {
    ExecutiveGroupTranslateTTT(I->G, I, v, store);
  } else {
    if(!I->TTTFlag) {
      I->TTTFlag = true;
      initializeTTT44f(I->TTT);
    }
    if(v) {
      I->TTT[3] += v[0];
      I->TTT[7] += v[1];
      I->TTT[11] += v[2];
    }
    if(store<0) 
      store = SettingGet_i(I->G, I->Setting.get(), NULL, cSetting_movie_auto_store);
    if(store && MovieDefined(I->G)) {
      if(!I->ViewElem)  
        I->ViewElem = pymol::vla<CViewElem>(0);
      if(I->ViewElem) { /* update motion path waypoint, if active */
        int frame = SceneGetFrame(I->G);
        if(frame >= 0) {
          VLACheck(I->ViewElem, CViewElem, frame);
          TTTToViewElem(I->TTT, I->ViewElem + frame);
          I->ViewElem[frame].specification_level = 2;
        }
      }
    }
  }
}


/*========================================================================*/
void ObjectSetTTT(pymol::CObject * I, const float *ttt, int state, int store)
{
  if(state < 0) {
    if(ttt) {
      UtilCopyMem(I->TTT, ttt, sizeof(float) * 16);
      I->TTTFlag = true;
    } else {
      I->TTTFlag = false;
      return;
    }
    if(store<0) 
      store = SettingGet_i(I->G, I->Setting.get(), NULL, cSetting_movie_auto_store);
    if(store && MovieDefined(I->G)) {
      if(!I->ViewElem)  
        I->ViewElem = pymol::vla<CViewElem>(0);
      if(I->ViewElem) { /* update motion path waypoint, if active */
        int frame = SceneGetFrame(I->G);
        if(frame >= 0) {
          VLACheck(I->ViewElem, CViewElem, frame);
          TTTToViewElem(I->TTT, I->ViewElem + frame);
          I->ViewElem[frame].specification_level = 2;
        }
      }
    }
  } else {
    /* to do */
  }
}

/*========================================================================*/
int ObjectGetTTT(pymol::CObject * I, const float **ttt, int state)
{
  if(state < 0) {
    if(I->TTTFlag) {
      *ttt = I->TTT;
      return 1;
    } else {
      *ttt = NULL;
    }

  } else {
  }
  return 0;
}


/*========================================================================*/
void ObjectResetTTT(pymol::CObject * I,int store)
{
  
  I->TTTFlag = false;
  if(store<0) 
    store = SettingGet_i(I->G, I->Setting.get(), NULL, cSetting_movie_auto_store);
  if(store && MovieDefined(I->G)) {
    if(!I->ViewElem)  
      I->ViewElem = pymol::vla<CViewElem>(0);
    if(I->ViewElem) { /* update motion path waypoint, if active */
      int frame = SceneGetFrame(I->G);
      if(frame >= 0) {
        identity44f(I->TTT);
        VLACheck(I->ViewElem, CViewElem, frame);
        TTTToViewElem(I->TTT, I->ViewElem + frame);
        I->ViewElem[frame].specification_level = 2;
      }
    }
  }
}


/*========================================================================*/
/**
 * Get the combined transformation of TTT and state matrix. State matrix is
 * only included if `history=true` or `matrix_mode > 0`.
 *
 * @param state See pymol::CObject::getObjectState
 * @param history Boolean flag
 * @param[out] matrix Homogeneous 4x4 matrix
 * @return True if `matrix` was populated
 */
int ObjectGetTotalMatrix(pymol::CObject * I, int state, int history, double *matrix)
{
  int result = false;
  if(I->TTTFlag) {
    convertTTTfR44d(I->TTT, matrix);
    result = true;
  }

  if (!history) {
    history =
        SettingGet<int>(I->G, I->Setting.get(), nullptr, cSetting_matrix_mode) > 0;
  }

  if (history) {
    {
      CObjectState* obj_state = I->getObjectState(state);
      if (obj_state) {
          if(!obj_state->Matrix.empty()) {
            const double *state_matrix = obj_state->Matrix.data();
            if(result) {
              right_multiply44d44d(matrix, state_matrix);
            } else {
              copy44d(state_matrix, matrix);
            }
            result = true;
          }
      }
    }
  }
  return result;
}


/*========================================================================*/
void ObjectPrepareContext(pymol::CObject * I, RenderInfo * info)
{
  CRay * ray = info ? info->ray : NULL;

  if(I->ViewElem) {
    int frame = SceneGetFrame(I->G);
    if(frame >= 0) {
      VLACheck(I->ViewElem, CViewElem, frame);
      
      if(I->Grabbed) {
        TTTToViewElem(I->TTT, I->ViewElem + frame);
        I->ViewElem[frame].specification_level = 2;
      } else {
        if(I->ViewElem[frame].specification_level) {
          TTTFromViewElem(I->TTT, I->ViewElem + frame);
          I->TTTFlag = true;
        }
        if(I->ViewElem[frame].state_flag) {
          SettingCheckHandle(I->G, I->Setting);
          if(I->Setting) {
            /* note: this assumes that the state has already been
               calculated and can thus be displayed.  How can we
               guarantee this to be true? */
            SettingSet_i(I->Setting.get(),cSetting_state,I->ViewElem[frame].state + 1);
          }
        }
      }
    }
  }
  if(ray) {
    RaySetTTT(ray, I->TTTFlag, I->TTT);
  } else {
    PyMOLGlobals *G = I->G;
    if(G->HaveGUI && G->ValidContext) {
      if(I->TTTFlag) {
        /* convert the row-major TTT matrix to a column-major OpenGL matrix */
        float gl[16], *ttt;

        ttt = I->TTT;
        gl[0] = ttt[0];
        gl[4] = ttt[1];
        gl[8] = ttt[2];
        gl[12] = ttt[3];
        gl[1] = ttt[4];
        gl[5] = ttt[5];
        gl[9] = ttt[6];
        gl[13] = ttt[7];
        gl[2] = ttt[8];
        gl[6] = ttt[9];
        gl[10] = ttt[10];
        gl[14] = ttt[11];
        gl[3] = 0.0;
        gl[7] = 0.0;
        gl[11] = 0.0;
        gl[15] = 1.0;

        auto mvm = SceneGetModelViewMatrix(G);
        MatrixMultiplyC44f(gl, mvm);
        MatrixTranslateC44f(mvm, ttt[12], ttt[13], ttt[14]);

#ifndef PURE_OPENGL_ES_2
        if (ALWAYS_IMMEDIATE_OR(!info->use_shaders)) {
          glLoadMatrixf(mvm);
        }
#endif
      }
    }
  }
}


/*========================================================================*/
void ObjectSetTTTOrigin(pymol::CObject * I, float *origin)
{
  float homo[16];
  float *dst;
  float post[3];

  if(!I->TTTFlag) {
    I->TTTFlag = true;
    initializeTTT44f(I->TTT);
  }

  /* convert the existing TTT into a homogenous transformation matrix */

  convertTTTfR44f(I->TTT, homo);

  /* now reset to the passed-in origin */

  transform44f3fas33f3f(homo, origin, post);

  homo[3] += post[0];
  homo[7] += post[1];
  homo[11] += post[2];

  dst = homo + 12;

  invert3f3f(origin, dst);

  copy44f(homo, I->TTT);
}


/*========================================================================*/
pymol::copyable_ptr<CSetting>* pymol::CObject::getSettingHandle(int state)
{
  return &Setting;
}


/*========================================================================*/
std::string pymol::CObject::describeElement(int index) const
{
  return {};
}

/*========================================================================*/
void ObjectToggleRepVis(pymol::CObject * I, int rep)
{
  if((rep >= 0) && (rep < cRepCnt))
    I->visRep ^= (1 << rep);
}


/*========================================================================*/
void ObjectSetRepVisMask(pymol::CObject * I, int repmask, int value)
{
  switch (value) {
    case cVis_HIDE:
      I->visRep &= ~repmask;
      break;
    case cVis_SHOW:
      I->visRep |= repmask;
      break;
    case cVis_AS:
      I->visRep = repmask;
      break;
    case cVis_TOGGLE:
      I->visRep ^= repmask;
      break;
    default:
      printf("error: invalid value: %d\n", value);
  }
}


/*========================================================================*/
void ObjectSetName(pymol::CObject * I, const char *name)
{
  UtilNCopy(I->Name, name, WordLength);
  if(SettingGetGlobal_b(I->G, cSetting_validate_object_names))
    ObjectMakeValidName(I->G, I->Name);
}


/*========================================================================*/
pymol::CObject::~CObject()
{
  SceneObjectDel(this->G, this, false);
}


/*========================================================================*/
void ObjectUseColor(pymol::CObject * I)
{
  PyMOLGlobals *G = I->G;
  if(G->HaveGUI && G->ValidContext) {
    glColor3fv(ColorGet(I->G, I->Color));
  }
}

void ObjectUseColorCGO(CGO *cgo, pymol::CObject * I)
{
  PyMOLGlobals *G = I->G;
  if(G->HaveGUI && G->ValidContext) {
    CGOColorv(cgo, ColorGet(I->G, I->Color));
  }
}

/*========================================================================*/
/**
 * Render a unit box (dummy representation)
 */
void pymol::CObject::render(RenderInfo * info)
{
  if(G->HaveGUI && G->ValidContext) {
#ifdef PURE_OPENGL_ES_2
    /* TODO */
#else
    glBegin(GL_LINE_LOOP);
    glVertex3i(-1, -1, -1);
    glVertex3i(-1, -1, 1);
    glVertex3i(-1, 1, 1);
    glVertex3i(-1, 1, -1);

    glVertex3i(1, 1, -1);
    glVertex3i(1, 1, 1);
    glVertex3i(1, -1, 1);
    glVertex3i(1, -1, -1);
    glEnd();

    glBegin(GL_LINES);
    glVertex3i(0, 0, 0);
    glVertex3i(1, 0, 0);

    glVertex3i(0, 0, 0);
    glVertex3i(0, 3, 0);

    glVertex3i(0, 0, 0);
    glVertex3i(0, 0, 9);

    glEnd();
#endif
  }
}


/*========================================================================*/
pymol::CObject::CObject(PyMOLGlobals * G) : G(G)
{
  OrthoRemoveSplash(G);         /* HMM... this seems like an inappropriate sideeffect */
  visRep = cRepBitmask & ~(cRepCellBit | cRepExtentBit);
}

/*========================================================================*/

void ObjectStateInit(PyMOLGlobals * G, CObjectState * I)
{
  I->G = G;
}

void ObjectStatePurge(CObjectState * I)
{
}

int ObjectStateSetMatrix(CObjectState * I, const double *matrix)
{
  int ok = true;
  if(matrix) {
    I->Matrix.resize(16);
    copy44d(matrix, I->Matrix.data());
  } else {
    I->Matrix.clear();
  }
  I->InvMatrix.clear();
  return ok;
}

void ObjectStateRightCombineMatrixR44d(CObjectState * I, const double *matrix)
{
  if(matrix) {
    if(I->Matrix.empty()) {
      I->Matrix = std::vector<double>(16);
      copy44d(matrix, I->Matrix.data());
    } else {
      right_multiply44d44d(I->Matrix.data(), matrix);
    }
  }
  I->InvMatrix.clear();
}

void ObjectStateLeftCombineMatrixR44d(CObjectState * I, const double *matrix)
{
  if(matrix) {
    if(I->Matrix.empty()) {
      I->Matrix = std::vector<double>(16);
      copy44d(matrix, I->Matrix.data());
    } else {
      left_multiply44d44d(matrix, I->Matrix.data());
    }
  }
  I->InvMatrix.clear();
}

void ObjectStateCombineMatrixTTT(CObjectState * I, const float *matrix)
{

  if(matrix) {
    if(I->Matrix.empty()) {
      I->Matrix = std::vector<double>(16);
      convertTTTfR44d(matrix, I->Matrix.data());
    } else {
      double tmp[16];
      convertTTTfR44d(matrix, tmp);
      right_multiply44d44d(I->Matrix.data(), tmp);
    }
  }
  I->InvMatrix.clear();
}

double *ObjectStateGetMatrix(CObjectState * I)
{
  if(!I->Matrix.empty()) {
    return I->Matrix.data();
  }
  return nullptr;
}

const double* ObjectStateGetMatrix(const CObjectState* I)
{
  return I->Matrix.empty() ? nullptr : I->Matrix.data();
}

/**
 * Get the Matrix inverse
 */
const double *ObjectStateGetInvMatrix(const CObjectState * self)
{
  if (self->Matrix.empty()) {
    return nullptr;
  }

  if (self->InvMatrix.empty()) {
    auto I = const_cast<CObjectState*>(self);
    I->InvMatrix = std::vector<double>(16);
    xx_matrix_invert(I->InvMatrix.data(), I->Matrix.data(), 4);
  }

  return self->InvMatrix.data();
}

void ObjectStateTransformMatrix(CObjectState * I, const double *matrix)
{
  if(I->Matrix.empty()) {
    I->Matrix = std::vector<double>(16);
    if(!I->Matrix.empty()) {
      copy44d(matrix, I->Matrix.data());
    }
  } else {
    right_multiply44d44d(I->Matrix.data(), matrix);
  }
  I->InvMatrix.clear();
}

int ObjectStatePushAndApplyMatrix(CObjectState * I, RenderInfo * info)
{
  PyMOLGlobals *G = I->G;
  float matrix[16];
  const double *i_matrix = nullptr;
  if(!I->Matrix.empty()) {
    i_matrix = I->Matrix.data();
  }
  int result = false;
  if(i_matrix) {
    if(info->ray) {
      float ttt[16], matrix[16], i_matrixf[16];
      RayPushTTT(info->ray);
      RayGetTTT(info->ray, ttt);
      convertTTTfR44f(ttt, matrix);
      copy44d44f(i_matrix, i_matrixf);
      right_multiply44f44f(matrix, i_matrixf);
      RaySetTTT(info->ray, true, matrix);
      result = true;
    } else if(G->HaveGUI && G->ValidContext) {
      matrix[0] = i_matrix[0];
      matrix[1] = i_matrix[4];
      matrix[2] = i_matrix[8];
      matrix[3] = i_matrix[12];
      matrix[4] = i_matrix[1];
      matrix[5] = i_matrix[5];
      matrix[6] = i_matrix[9];
      matrix[7] = i_matrix[13];
      matrix[8] = i_matrix[2];
      matrix[9] = i_matrix[6];
      matrix[10] = i_matrix[10];
      matrix[11] = i_matrix[14];
      matrix[12] = i_matrix[3];
      matrix[13] = i_matrix[7];
      matrix[14] = i_matrix[11];
      matrix[15] = i_matrix[15];

      ScenePushModelViewMatrix(G);
      auto mvm = SceneGetModelViewMatrix(G);
      MatrixMultiplyC44f(matrix, mvm);

#ifndef PURE_OPENGL_ES_2
      if (ALWAYS_IMMEDIATE_OR(!info->use_shaders)) {
        glLoadMatrixf(mvm);
      }
#endif

      result = true;
    }
  }
  return result;
}

void ObjectStatePopMatrix(CObjectState * I, RenderInfo * info)
{
  PyMOLGlobals *G = I->G;
  if(info->ray) {
    RayPopTTT(info->ray);
  } else if(G->HaveGUI && G->ValidContext) {
    ScenePopModelViewMatrix(G, !info->use_shaders);
  }
}

void ObjectStateResetMatrix(CObjectState* I)
{
  I->Matrix.clear();
  I->InvMatrix.clear();
}

PyObject *ObjectStateAsPyList(CObjectState * I)
{
  PyObject *result = NULL;

  if(I) {
    result = PyList_New(1);

    if(!I->Matrix.empty()) {
      PyList_SetItem(result, 0, PConvDoubleArrayToPyList(I->Matrix.data(), 16));
    } else {
      PyList_SetItem(result, 0, PConvAutoNone(Py_None));
    }
  }
  return (PConvAutoNone(result));
}

int ObjectStateFromPyList(PyMOLGlobals * G, PyObject * list, CObjectState * I)
{
  PyObject *tmp;
  int ok = true;

  ObjectStateInit(G, I);

  if(list && (list != Py_None)) {       /* allow None */
    if(ok)
      ok = (list != NULL);
    if(ok)
      ok = PyList_Check(list);
    /* TO SUPPORT BACKWARDS COMPATIBILITY...
       Always check ll when adding new PyList_GetItem's */
    if(ok) {
      tmp = PyList_GetItem(list, 0);
      if(tmp != Py_None)
        ok = PConvFromPyObject(G, tmp, I->Matrix);
    }
  }
  return (ok);
}
