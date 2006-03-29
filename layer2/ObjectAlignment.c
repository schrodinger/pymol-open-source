/* 
A* -------------------------------------------------------------------
B* This file contains source code for the PyMOL computer program
C* copyright 1998-2006 by Warren Lyford Delano of DeLano Scientific. 
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

#include"os_predef.h"
#include"os_std.h"
#include"os_gl.h"

#include"OOMac.h"
#include"ObjectAlignment.h"
#include"Base.h"
#include"MemoryDebug.h"
#include"CGO.h"
#include"Scene.h"
#include"Setting.h"
#include"PConv.h"
#include"main.h"
#include"Color.h"
#include"Executive.h"
#include"OVContext.h"
#include"Util.h"
#include"Selector.h"

static ObjectAlignment *ObjectAlignmentNew(PyMOLGlobals *G);
static void ObjectAlignmentFree(ObjectAlignment *I);

static int GroupOrderKnown(ExecutiveObjectOffset *eoo, 
                            OVOneToOne *id2eoo,
                            int *curVLA,
                            int *newVLA,
                            int cur_start,
                            int new_start,
                            ObjectMolecule *guide,
                            int *action)
{
  register int order_known = false;
  if(guide) {
    register int c,id;
    OVreturn_word offset;
    
    register int cur_offset = -1;
    register int new_offset = -1;
    
    /* find lowest offset within the cur group */
    c = cur_start;
    while( (id=curVLA[c++]) ) {
      if(OVreturn_IS_OK( offset=OVOneToOne_GetForward(id2eoo, id))) {
        if( eoo[offset.word].obj == guide ) {
          if((cur_offset<0) || (eoo[offset.word].offset<cur_offset))
            cur_offset = eoo[offset.word].offset;
        }
      }
    }
    
    /* find lowest offset within the new group */
    c = new_start;
    while( (id=newVLA[c++]) ) {
      if(OVreturn_IS_OK( offset=OVOneToOne_GetForward(id2eoo, id))) {
        if( eoo[offset.word].obj == guide ) {
          if((new_offset<0) || (eoo[offset.word].offset<new_offset))
            new_offset = eoo[offset.word].offset;
        }
      }
    }
    
    if((new_offset>=0) && (cur_offset>=0)) {
      if(new_offset<cur_offset) {
        order_known = true;
        *action = -1;
      } else if(new_offset>cur_offset) {
        order_known = true;
        *action = 1;
      }
    }
  }
  return order_known;
}

static int *AlignmentMerge(PyMOLGlobals *G, int *curVLA, int *newVLA, ObjectMolecule *guide)
{
  /* curVLA and newVLA must be properly sized and zero terminated...*/
  int *result = NULL;
  int n_result = 0;

  {
    ExecutiveObjectOffset *eoo = NULL;
    OVOneToOne *id2eoo = NULL;
    
    if(ExecutiveGetUniqueIDObjectOffsetVLADict(G, &eoo, &id2eoo)) {
      
      /* now we can take unique IDs to specific atoms */
      OVOneToAny *used = OVOneToAny_New(G->Context->heap);
      OVOneToAny *active = OVOneToAny_New(G->Context->heap);
      int n_cur = VLAGetSize(curVLA);
      int n_new = VLAGetSize(newVLA);
      int cur_start = 0;
      int new_start = 0;

      result = VLAlloc(int, ( (n_cur<n_new) ? n_new : n_cur));
      
      while( (cur_start<n_cur)||(new_start<n_new) ) {
        
        int action; /* -1 = insert new, 0 = merge, 1 = insert cur */

        /* make sure both lists are queued up on the next identifier */
        
        while( (cur_start<n_cur) && !curVLA[cur_start]) 
          cur_start++;
        while( (new_start<n_cur) && !newVLA[new_start]) 
          new_start++;
        
        if(newVLA[new_start]) { /* default is to insert new first...*/
          action = -1;
        } else {
          action = 1;
        }
        
        if( (cur_start<n_cur)&& (new_start<n_new) &&
            curVLA[cur_start] && newVLA[new_start]) {
           /* both lists active */
          
          register int c, id;
          int overlapping = false;

          OVOneToAny_Reset(active);
          c = cur_start;
          while( (id=curVLA[c++]) ) { /* record active atoms */
            OVOneToAny_SetKey(active,id,1);
          }
          
          c = new_start;
          while( (id=newVLA[c++]) ) { /* see if there are any matches */
            if(OVreturn_IS_OK(OVOneToAny_GetKey(active,id))) { 
              overlapping = true; 
              break;
            }
          }

          if(overlapping) {
            /* overlapping, so merge */
            action = 0; 

          } else { 
            /* non-overlapping, so we need to figure out which goes first... */
            if(!GroupOrderKnown(eoo,id2eoo,curVLA,newVLA,
                                cur_start,new_start,guide,&action)) {
              register int c,id;
              OVreturn_word offset;
              ObjectMolecule *obj, *last_obj = NULL;
              c = cur_start;
              while( (id=curVLA[c++]) ) {
                
                if(OVreturn_IS_OK( offset=OVOneToOne_GetForward(id2eoo, id)) ) {
                  obj = eoo[offset.word].obj;
                  if(obj != last_obj) {
                    if(GroupOrderKnown(eoo,id2eoo,curVLA,newVLA,
                                       cur_start,new_start,obj,&action))
                      break;
                    else
                      last_obj = obj;
                  }
                }
              }

              /* if order isn't set by now, then it doesn't matter...
                 so group new will go in first */
            }
          }
        }

        {
          register int id;
          
          switch(action) {
          case -1: /* insert new */
            if(new_start<n_new) {
              while( (id=newVLA[new_start]) ) {
                if(OVOneToAny_GetKey(used,id).status == OVstatus_NOT_FOUND) {
                  if(OVreturn_IS_OK(OVOneToAny_SetKey(used,id,1))) {
                    VLACheck(result,int,n_result);
                    result[n_result] = id;
                    n_result++;
                  }
                }
                new_start++;
              }
              while( (new_start<n_new) && (!newVLA[new_start])) 
                new_start++;
            }
            VLACheck(result,int,n_result);
            result[n_result] = 0;
            n_result++;
            break;
          case 0: /* merge, with cur going first */
            if(new_start<n_new) {
              while( (id=newVLA[new_start]) ) {
                if(OVOneToAny_GetKey(used,id).status == OVstatus_NOT_FOUND) {
                  if(OVreturn_IS_OK(OVOneToAny_SetKey(used,id,1))) {
                    VLACheck(result,int,n_result);
                    result[n_result] = id;
                    n_result++;
                  }
                }
                new_start++;
              }
              while( (new_start<n_new) && (!newVLA[new_start])) 
                new_start++;
            }
            if(cur_start<n_cur) {
              while( (id=curVLA[cur_start]) ) {
                if(OVOneToAny_GetKey(used,id).status == OVstatus_NOT_FOUND) {
                  if(OVreturn_IS_OK(OVOneToAny_SetKey(used,id,1))) {
                    VLACheck(result,int,n_result);
                    result[n_result] = id;
                    n_result++;
                  }
                }
                cur_start++;
              }
              while( (cur_start<n_cur) && (!curVLA[cur_start])) 
                cur_start++;
            }
            VLACheck(result,int,n_result);
            result[n_result] = 0;
            n_result++;
            break;
          case 1: /* insert cur */
            if(cur_start<n_cur) {
              while( (id=curVLA[cur_start]) ) {
                if(OVOneToAny_GetKey(used,id).status == OVstatus_NOT_FOUND) {
                  if(OVreturn_IS_OK(OVOneToAny_SetKey(used,id,1))) {
                    VLACheck(result,int,n_result);
                    result[n_result] = id;
                    n_result++;
                  }
                }
                cur_start++;
              }
              while( (cur_start<n_cur) && (!curVLA[cur_start])) 
                cur_start++;
              VLACheck(result,int,n_result);
              result[n_result] = 0;
              n_result++;
            }
            break;
          }
        }
      }
      OVOneToAny_DEL_AUTO_NULL(active);
      OVOneToAny_DEL_AUTO_NULL(used);
    }
    OVOneToOne_DEL_AUTO_NULL(id2eoo);
    VLAFreeP(eoo);
  }
  if(result && n_result && (!result[n_result-1])) {
    VLACheck(result,int,n_result);
    result[n_result] = 0;
    n_result++;
  }
  VLASize(result,int,n_result);
  return result;
}

#ifndef _PYMOL_NOPY
static PyObject *ObjectAlignmentStateAsPyList(ObjectAlignmentState *I)
{

  PyObject *result = NULL;

  result = PyList_New(0);
#if 0
  if(I->std) 
    PyList_SetItem(result,0,CGOAsPyList(I->std));
  else
    PyList_SetItem(result,0,PConvAutoNone(NULL));
  if(I->ray) 
    PyList_SetItem(result,1,CGOAsPyList(I->ray));
  else
    PyList_SetItem(result,1,PConvAutoNone(NULL));
#endif
  return(PConvAutoNone(result));  

}

static PyObject *ObjectAlignmentAllStatesAsPyList(ObjectAlignment *I)
{

  PyObject *result=NULL;
  int a;
  result = PyList_New(I->NState);
  for(a=0;a<I->NState;a++) {
    PyList_SetItem(result,a,ObjectAlignmentStateAsPyList(I->State+a));
  }
  return(PConvAutoNone(result));  

}

static int ObjectAlignmentStateFromPyList(PyMOLGlobals *G,ObjectAlignmentState *I,PyObject *list,int version)
{
  int ok=true;
  int ll;
  if(ok) ok=(list!=NULL);
  if(ok) ok=PyList_Check(list);
  if(ok) ll = PyList_Size(list);
  /* TO SUPPORT BACKWARDS COMPATIBILITY...
   Always check ll when adding new PyList_GetItem's */
#if 0
  if(ok) {
    tmp = PyList_GetItem(list,0);
    if(tmp == Py_None)
      I->std = NULL;
    else 
      ok = ((I->std=CGONewFromPyList(G,PyList_GetItem(list,0),version))!=NULL);
  }
  if(ok) {
    tmp = PyList_GetItem(list,1);
    if(tmp == Py_None)
      I->ray = NULL;
    else 
      ok = ((I->ray=CGONewFromPyList(G,PyList_GetItem(list,1),version))!=NULL);
  }
#endif
  return(ok);
}

static int ObjectAlignmentAllStatesFromPyList(ObjectAlignment *I,PyObject *list,int version)
{
  int ok=true;
  int a;
  VLACheck(I->State,ObjectAlignmentState,I->NState);
  if(ok) ok=PyList_Check(list);
  if(ok) {
    for(a=0;a<I->NState;a++) {
      ok = ObjectAlignmentStateFromPyList(I->Obj.G,I->State+a,PyList_GetItem(list,a),version);
      if(!ok) break;
    }
  }
  return(ok);
}
#endif

int ObjectAlignmentNewFromPyList(PyMOLGlobals *G,PyObject *list,ObjectAlignment **result,int version)
{
#ifdef _PYMOL_NOPY
  return 0;
#else

  int ok = true;
  ObjectAlignment *I=NULL;
  (*result) = NULL;
  if(ok) ok=(list!=Py_None);
  if(ok) ok=PyList_Check(list);

  I=ObjectAlignmentNew(G);
  if(ok) ok = (I!=NULL);

  if(ok) ok = ObjectFromPyList(G,PyList_GetItem(list,0),&I->Obj);
  if(ok) ok = PConvPyIntToInt(PyList_GetItem(list,1),&I->NState);
  if(ok) ok = ObjectAlignmentAllStatesFromPyList(I,PyList_GetItem(list,2),version);
  if(ok) {
    (*result) = I;
    ObjectAlignmentRecomputeExtent(I);
  } else {
    /* cleanup? */
  }
  return(ok);
#endif
}

PyObject *ObjectAlignmentAsPyList(ObjectAlignment *I)
{
#ifdef _PYMOL_NOPY
  return NULL;
#else

  PyObject *result=NULL;

  result = PyList_New(3);
  PyList_SetItem(result,0,ObjectAsPyList(&I->Obj));
  PyList_SetItem(result,1,PyInt_FromLong(I->NState));
  PyList_SetItem(result,2,ObjectAlignmentAllStatesAsPyList(I));

  return(PConvAutoNone(result));  
#endif
}


/*========================================================================*/

static void ObjectAlignmentFree(ObjectAlignment *I) {
  int a;
  for(a=0;a<I->NState;a++) {
    if(I->State[a].std)
      CGOFree(I->State[a].std);
    if(I->State[a].ray)
      CGOFree(I->State[a].ray);
    VLAFreeP(I->State[a].alignVLA);
    OVOneToAny_DEL_AUTO_NULL(I->State[a].id2tag);
  }
  VLAFreeP(I->State);
  ObjectPurge(&I->Obj);
  OOFreeP(I);
}

/*========================================================================*/

void ObjectAlignmentRecomputeExtent(ObjectAlignment *I)
{
  float mx[3],mn[3];
  int extent_flag = false;
  int a;
  for(a=0;a<I->NState;a++) 
    if(I->State[a].std) {
      if(CGOGetExtent(I->State[a].std,mn,mx)) {
        if(!extent_flag) {
          extent_flag=true;
          copy3f(mx,I->Obj.ExtentMax);
          copy3f(mn,I->Obj.ExtentMin);
        } else {
          max3f(mx,I->Obj.ExtentMax,I->Obj.ExtentMax);
          min3f(mn,I->Obj.ExtentMin,I->Obj.ExtentMin);
        }
      }
    }
  I->Obj.ExtentFlag=extent_flag;
}
/*========================================================================*/
static void ObjectAlignmentUpdate(ObjectAlignment *I)
{
  register PyMOLGlobals *G = I->Obj.G;
  int update_needed = false;
  {
    int a;
    
    for(a=0;a<I->NState;a++) {
      ObjectAlignmentState *oas = I->State + a;
      if(!oas->valid) update_needed=true;
    }
  }

  if(update_needed) {

    ExecutiveObjectOffset *eoo = NULL;
    OVOneToOne *id2eoo = NULL;
    
    if(ExecutiveGetUniqueIDObjectOffsetVLADict(G, &eoo, &id2eoo)) {

      int a;
      for(a=0;a<I->NState;a++) {
        ObjectAlignmentState *oas = I->State + a;
        if(!oas->valid) {

          if(I->SelectionState == a )
            I->SelectionState = -1;

          if(oas->std) {
            CGOFree(oas->std);
            oas->std = NULL;
          }
          if(oas->ray) {
            CGOFree(oas->ray);
            oas->ray = NULL;
          }
          if(oas->id2tag) {
            OVOneToAny_Reset(oas->id2tag);
          } else {
            oas->id2tag = OVOneToAny_New(G->Context->heap);
          }

          {
            CGO *cgo = CGONew(G);
        
            if(oas->alignVLA) {
              int id, b=0, c;
              int *vla = oas->alignVLA;
              int n_id = VLAGetSize(vla);
              float mean[3], vert[3];
              int n_coord = 0;
              int tag = SELECTOR_BASE_TAG+1;
              OVOneToAny *id2tag = oas->id2tag;
              OVreturn_word offset;

              while(b<n_id) {

                while( (b<n_id) && (!vla[b]) ) 
                  b++;

                if(!(b<n_id))
                  break;

                c = b;
                n_coord = 0;
                zero3f(mean);
                while( (id=vla[c++]) ) {
                  if(OVreturn_IS_OK( offset=OVOneToOne_GetForward(id2eoo, id)) ) {
                    if(ObjectMoleculeGetAtomVertex(eoo[offset.word].obj, a, 
                                                   eoo[offset.word].offset,vert)) {
                      n_coord++;
                      add3f(vert,mean,mean);
                    }
                  }
                }

                if(n_coord>2) { /* >2 points, then draw to mean */
                  float scale = 1.0F/n_coord;

                  scale3f(mean, scale, mean);

                  CGOBegin(cgo,GL_LINES);

                  c = b;
                  while( (id=vla[c++]) ) {
                    if(OVreturn_IS_OK( offset=OVOneToOne_GetForward(id2eoo, id) )) {
                      if(ObjectMoleculeGetAtomVertex(eoo[offset.word].obj, a, 
                                                     eoo[offset.word].offset,vert)) {
                        CGOVertexv(cgo,mean);
                        CGOVertexv(cgo,vert);
                      }
                    }
                  }
                  CGOEnd(cgo);
                 
                } else if(n_coord) { /* if 2 points, then simply draw a line */
                  float first[3];
                  int first_flag = true;
                  CGOBegin(cgo,GL_LINES);

                  c = b;

                  while( (id=vla[c++]) ) {
                    if(OVreturn_IS_OK( offset=OVOneToOne_GetForward(id2eoo, id) )) {
                      if(ObjectMoleculeGetAtomVertex(eoo[offset.word].obj, a, 
                                                     eoo[offset.word].offset,vert)) {
                        if(first_flag) {
                          copy3f(vert,first);
                          first_flag=false;
                        } else {
                          CGOVertexv(cgo,first);
                          CGOVertexv(cgo,vert);
                        }
                      }
                    }
                  }
                  CGOEnd(cgo);
                }

                /* update the it2tag dictionary */

                tag++;

                while( (b<n_id) && vla[b] ) {
                  OVOneToAny_SetKey(id2tag, vla[b], tag);
                  b++;
                }
              }
            }
        
            CGOStop(cgo);
        
            /* simplify if necessary */

            {
              int est=CGOCheckComplex(cgo);
              if(est) {
                oas->ray = cgo;
                oas->std = CGOSimplify(cgo,est);
              } else 
                oas->std = cgo;
            }
          }
          oas->valid = true;
        }
      }
    }
    OVOneToOne_DEL_AUTO_NULL(id2eoo);
    VLAFreeP(eoo);
  }
  if(I->SelectionState<0) {
    int state=SettingGet_i(I->Obj.G,NULL,I->Obj.Setting,cSetting_state)-1;    
    if(state<0)
      state = SceneGetState(G);
    if(state>I->NState)
      state = I->NState-1;
    if(state<0)
      state=0;
    if(state<I->NState) {
      ObjectAlignmentState *oas = I->State + state;
      if(oas->id2tag) {
        SelectorDelete(G,I->Obj.Name);
        SelectorCreateFromTagDict(G, I->Obj.Name, oas->id2tag, false);
      }
    }
  }
  SceneInvalidate(I->Obj.G); 
}
/*========================================================================*/

static int ObjectAlignmentGetNState(ObjectAlignment *I) {
  return(I->NState);
}

/*========================================================================*/

static void ObjectAlignmentRender(ObjectAlignment *I,RenderInfo *info)
{
  register PyMOLGlobals *G = I->Obj.G;
  int state = info->state;
  CRay *ray = info->ray;
  Picking **pick = info->pick;
  int pass = info->pass;
  ObjectAlignmentState *sobj = NULL;
  int a;
  float *color;

  ObjectPrepareContext(&I->Obj,ray);

  color = ColorGet(G,I->Obj.Color);

  if(!pass) {
    if(I->Obj.RepVis[cRepCGO]) {
      
      if(state<I->NState) {
        sobj = I->State+state;
      }
      if(state<0) {
        if(I->State) {
          for(a=0;a<I->NState;a++) {
            sobj = I->State+a;
            if(ray) {    
              if(sobj->ray)
                CGORenderRay(sobj->ray,ray,color,I->Obj.Setting,NULL);
              else
                CGORenderRay(sobj->std,ray,color,I->Obj.Setting,NULL);
            } else if(G->HaveGUI && G->ValidContext) {
              if(pick) {
              } else {
                if(sobj->std)
                  CGORenderGL(sobj->std,color,I->Obj.Setting,NULL,info);
              }
            }
          }
        }
      } else {
        if(!sobj) {
          if(I->NState&&SettingGet(G,cSetting_static_singletons)) 
            sobj = I->State;
        }
        if(ray) {    
          if(sobj)
            {
              if(sobj->ray)
                CGORenderRay(sobj->ray,ray,color,I->Obj.Setting,NULL);
              else
                CGORenderRay(sobj->std,ray,color,I->Obj.Setting,NULL);
            }
        } else if(G->HaveGUI && G->ValidContext) {
          if(pick) {
          } else {
            if(sobj) {
              if(sobj->std)
                CGORenderGL(sobj->std,color,I->Obj.Setting,NULL,info);
            }
          }
        }
      }
    }
  }
}

static void ObjectAlignmentInvalidate(ObjectAlignment *I,int rep,int level,int state)
{
  if((rep==cRepAll) || (rep==cRepCGO)) {
    if(state>=0) {
      if(state<I->NState)
        I->State[state].valid = false;
    } else {
      int a;
      for(a=0;a<I->NState;a++) {
        I->State[a].valid = false;
      }
    }
  }
}
/*========================================================================*/
static ObjectAlignment *ObjectAlignmentNew(PyMOLGlobals *G)
{
  OOAlloc(G,ObjectAlignment);

  ObjectInit(G,(CObject*)I);

  I->State=VLAMalloc(10,sizeof(ObjectAlignmentState),5,true); /* auto-zero */
  I->NState=0;
  I->SelectionState=-1;

  I->Obj.type = cObjectAlignment;
  I->Obj.fFree = (void (*)(struct CObject *))ObjectAlignmentFree;
  I->Obj.fUpdate =(void (*)(struct CObject *)) ObjectAlignmentUpdate;
  I->Obj.fRender =(void (*)(struct CObject *, RenderInfo *))ObjectAlignmentRender;
  I->Obj.fGetNFrame = (int (*)(struct CObject *)) ObjectAlignmentGetNState;
  I->Obj.fInvalidate = (void (*)(struct CObject *,int rep, int level, int state))
    ObjectAlignmentInvalidate;

  return(I);
}

/*========================================================================*/
ObjectAlignment *ObjectAlignmentDefine(PyMOLGlobals *G,
                                       ObjectAlignment *obj,
                                       int *align_vla,
                                       int state, 
                                       int merge,
                                       ObjectMolecule *guide)
{
  ObjectAlignment *I = NULL;
  
  if(obj) {
    if(obj->Obj.type!=cObjectAlignment) /* TODO: handle this */
      obj=NULL;
  }
  if(!obj) {
    I=ObjectAlignmentNew(G);
  } else {
    I=obj;
  }
  if(state<0) state=I->NState;

  if(I->NState<=state) {
    VLACheck(I->State,ObjectAlignmentState,state);
    I->NState=state+1;
  }

  {
    ObjectAlignmentState *oas = I->State + state;
    oas->valid = false;

    if(align_vla) {
      if(merge && oas->alignVLA) {
        int *new_vla = AlignmentMerge(G, oas->alignVLA, align_vla, guide);
        if(new_vla) {
          VLAFreeP(oas->alignVLA);
          oas->alignVLA = new_vla;
        }
      } else {
        int size = VLAGetSize(align_vla);
        if(oas->alignVLA)
          VLAFreeP(oas->alignVLA);
        oas->alignVLA = VLAlloc(int,size);
        UtilCopyMem(oas->alignVLA, align_vla, sizeof(int)*size);
        VLASize(oas->alignVLA, int, size);
      }
    } else {
      VLAFreeP(oas->alignVLA);
    }
  }
#if 0

  if(PyList_Check(pycgo)) {
    if(PyList_Size(pycgo)) {
      if(PyFloat_Check(PyList_GetItem(pycgo,0))) {
        cgo=ObjectAlignmentPyListFloatToCGO(G,pycgo);
        if(cgo) {
          est=CGOCheckForText(cgo);
          if(est) {
            CGOPreloadFonts(cgo);
            font_cgo = CGODrawText(cgo,est,NULL);
            CGOFree(cgo);
            cgo=font_cgo;
          }
          est=CGOCheckComplex(cgo);
          if(est) {
            I->State[state].ray=cgo;
            I->State[state].std=CGOSimplify(cgo,est);
          } else 
            I->State[state].std=cgo;
          
        } else {
          ErrMessage(G,"ObjectAlignment","could not parse CGO List.");
        }
      }
    }
  }
#endif
  if(I) {
    ObjectAlignmentRecomputeExtent(I);
  }
  SceneChanged(G);
  SceneCountFrames(G);
  return(I);
}

