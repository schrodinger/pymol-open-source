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

/* 
   NOTICE:

   Important thread safety tip:

   PM operations which will ultimately call GLUT can only be called by
   the main GLUT thread (with some exceptions, such as the simple
   drawing operations which seem to be thread safe).

   Thus, pm.py needs to guard against direct invocation of certain _cmd
   (API) calls from outside threads [Examples: _cmd.png(..),
   _cmd.mpng(..) ].  Instead, it needs to hand them over to the main
   thread by way of a cmd.mdo(...)  statement.

   Note that current, most glut operations have been pushed into the
   main event and redraw loops to avoid these kinds of problems - so
   I'm not sure how important this really is anymore.

*/
   

/* TODO: Put in some exception handling and reporting for the
 * python calls, especially, PyArg_ParseTuple()
 */

#include"os_predef.h"
#include"os_python.h"
#include"os_gl.h"
#include"os_std.h"

#include"MemoryDebug.h"
#include"Err.h"
#include"Util.h"
#include"Cmd.h"
#include"ButMode.h"
#include"Ortho.h"
#include"ObjectMolecule.h"
#include"ObjectMesh.h"
#include"ObjectMap.h"
#include"ObjectCallback.h"
#include"ObjectCGO.h"
#include"ObjectSurface.h"
#include"Executive.h"
#include"Selector.h"
#include"main.h"
#include"Scene.h"
#include"Setting.h"
#include"Movie.h"
#include"Export.h"
#include"P.h"
#include"PConv.h"
#include"Control.h"
#include"Editor.h"
#include"Wizard.h"
#include"SculptCache.h"
#include"TestPyMOL.h"
#include"Color.h"
#include"Seq.h"

#define cLoadTypePDB 0
#define cLoadTypeMOL 1
#define cLoadTypeSDF 2
#define cLoadTypeMOLStr 3
#define cLoadTypeMMD 4
#define cLoadTypeMMDSeparate 5
#define cLoadTypeMMDStr 6
#define cLoadTypeXPLORMap 7
#define cLoadTypeChemPyModel 8
#define cLoadTypePDBStr 9
#define cLoadTypeChemPyBrick 10
#define cLoadTypeChemPyMap 11
#define cLoadTypeCallback 12
#define cLoadTypeCGO 13
#define cLoadTypeR3D 14
#define cLoadTypeXYZ 15
#define cLoadTypeCCP4Map 18
#define cLoadTypePMO  19
#define cLoadTypeTOP  21
#define cLoadTypeTRJ  22
#define cLoadTypeCRD  23
#define cLoadTypeRST  24
#define cLoadTypePSE  25
#define cLoadTypeXPLORStr 26
#define cLoadTypePHIMap 27
#define cLoadTypeFLDMap 28
#define cLoadTypeBRIXMap 29
#define cLoadTypeGRDMap 30
#define cLoadTypePQR 31
#define cLoadTypeDXMap 32

#define tmpSele "_tmp"
#define tmpSele1 "_tmp1"
#define tmpSele2 "_tmp2"

static int flush_count = 0;

int run_only_once = true;

int PyThread_get_thread_ident(void);

/* NOTE: the glut_thread_keep_out variable can only be changed by the thread
   holding the API lock, therefore this is safe even through increment
   isn't (necessarily) atomic. */

static void APIEntry(void) /* assumes API is locked */
{
  PRINTFD(FB_API)
    " APIEntry-DEBUG: as thread 0x%x.\n",PyThread_get_thread_ident()
    ENDFD;

  if(PyMOLTerminating) {/* try to bail */
#ifdef WIN32
    abort();
#endif
    exit(0);
  }

  P_glut_thread_keep_out++;  
  PUnblock();
}

static void APIEnterBlocked(void) /* assumes API is locked */
{
  PRINTFD(FB_API)
    " APIEnterBlocked-DEBUG: as thread 0x%x.\n",PyThread_get_thread_ident()
    ENDFD;

  if(PyMOLTerminating) {/* try to bail */
#ifdef WIN32
    abort();
#endif
    exit(0);
  }

  P_glut_thread_keep_out++;  
}

static void APIExit(void) /* assumes API is locked */
{
  PBlock();
  P_glut_thread_keep_out--;
  PRINTFD(FB_API)
    " APIExit-DEBUG: as thread 0x%x.\n",PyThread_get_thread_ident()
    ENDFD;
}

static void APIExitBlocked(void) /* assumes API is locked */
{
  P_glut_thread_keep_out--;
  PRINTFD(FB_API)
    " APIExitBlocked-DEBUG: as thread 0x%x.\n",PyThread_get_thread_ident()
    ENDFD;
}

static PyObject *APISuccess(void)
{
  return(Py_BuildValue("i",1));
}

static PyObject *APIFailure(void)
{
  Py_INCREF(Py_None);
  return(Py_None);
}

static PyObject *APIStatus(int status) /* status/integer return */
{
  return(Py_BuildValue("i",status));
}

static PyObject *APIIncRef(PyObject *result) /* automatically own Py_None */
{
  Py_INCREF(result);
  return(result);
}
static PyObject *APIAutoNone(PyObject *result) /* automatically own Py_None */
{
  if(result==Py_None)
    Py_INCREF(result);
  else if(result==NULL) {
    result=Py_None;
    Py_INCREF(result);
  } 
  return(result);
}
static PyObject *CmdFixChemistry(PyObject *self, PyObject *args)
{
  char *str2,*str3;
  OrthoLineType s2="",s3="";
  int ok = false;
  int quiet;
  ok = PyArg_ParseTuple(args,"ssi",&str2,&str3,&quiet);
  if(ok) {
    APIEntry();
    SelectorGetTmp(str2,s2);
    SelectorGetTmp(str3,s3);
    ok = ExecutiveFixChemistry(s2,s3,quiet);
    SelectorFreeTmp(s2);
    SelectorFreeTmp(s3);
    APIExit();
  }
  return APIStatus(ok);
}

static PyObject *CmdGLDeleteLists(PyObject *self, PyObject *args)
{
  int int1,int2;
  int ok=false;
  ok = PyArg_ParseTuple(args,"ii",&int1,&int2);
  if(ok) {
    glDeleteLists(int1,int2);
  }
  return(APIStatus(1));
}
static PyObject *CmdRayAntiThread(PyObject *self, 	PyObject *args)
{
  int ok=true;
  PyObject *py_thread_info;

  CRayAntiThreadInfo *thread_info = NULL;

  ok = PyArg_ParseTuple(args,"O",&py_thread_info);
  if(ok) ok = PyCObject_Check(py_thread_info);
  if(ok) ok = ((thread_info = PyCObject_AsVoidPtr(py_thread_info))!=NULL);
  if (ok) {
    PUnblock();
    RayAntiThread(thread_info);
    PBlock();
  }
  return(APIStatus(ok));
}

static PyObject *CmdRayHashThread(PyObject *self, 	PyObject *args)
{
  int ok=true;
  PyObject *py_thread_info;

  CRayHashThreadInfo *thread_info = NULL;

  ok = PyArg_ParseTuple(args,"O",&py_thread_info);
  if(ok) ok = PyCObject_Check(py_thread_info);
  if(ok) ok = ((thread_info = PyCObject_AsVoidPtr(py_thread_info))!=NULL);
  if (ok) {
    PUnblock();
    RayHashThread(thread_info);
    PBlock();
  }
  return(APIStatus(ok));
}

static PyObject *CmdRayTraceThread(PyObject *self, 	PyObject *args)
{
  int ok=true;
  PyObject *py_thread_info;

  CRayThreadInfo *thread_info = NULL;

  ok = PyArg_ParseTuple(args,"O",&py_thread_info);
  if(ok) ok = PyCObject_Check(py_thread_info);
  if(ok) ok = ((thread_info = PyCObject_AsVoidPtr(py_thread_info))!=NULL);
  if (ok) {
    PUnblock();
    RayTraceThread(thread_info);
    PBlock();
  }
  return(APIStatus(ok));
}

static PyObject *CmdGetMovieLocked(PyObject *self, 	PyObject *args)
{
  return(APIStatus(MovieLocked()));
}

static PyObject *CmdFakeDrag(PyObject *self, 	PyObject *args)
{
  MainDragDirty();
  return(APIStatus(true));
}

static PyObject *CmdDelColorection(PyObject *dummy, PyObject *args)
{
  int ok=true;
  PyObject *list;
  char *prefix;
  ok = PyArg_ParseTuple(args,"Os",&list,&prefix);
  if (ok) {
    APIEnterBlocked();
    ok = SelectorColorectionFree(list,prefix);
    APIExitBlocked();
  }
  return(APIStatus(ok));
}

static PyObject *CmdSetColorection(PyObject *dummy, PyObject *args)
{
  int ok=true;
  char *prefix;
  PyObject *list;
  ok = PyArg_ParseTuple(args,"Os",&list,&prefix);
  if (ok) {
    APIEnterBlocked();
    ok = SelectorColorectionApply(list,prefix);
    APIExitBlocked();
  }
  return(APIStatus(ok));
}

static PyObject *CmdGetColorection(PyObject *dummy, PyObject *args)
{
  PyObject *result=NULL;
  int ok=true;
  char *prefix;
  ok = PyArg_ParseTuple(args,"s",&prefix);
  if (ok) {
    APIEnterBlocked();
    result = SelectorColorectionGet(prefix);
    APIExitBlocked();
  }
  return(APIAutoNone(result));
}

static PyObject *CmdGetVis(PyObject *dummy, PyObject *args)
{
  PyObject *result;
  int ok=true;
  if (ok) {
    APIEnterBlocked();
    result = ExecutiveGetVisAsPyDict();
    APIExitBlocked();
  }
  return(APIAutoNone(result));
}

static PyObject *CmdSetVis(PyObject *dummy, PyObject *args)
{
  int ok=true;
  PyObject *visDict;
  ok = PyArg_ParseTuple(args,"O",&visDict);
  if (ok) {
    APIEnterBlocked();
    ok = ExecutiveSetVisFromPyDict(visDict);
    APIExitBlocked();
  }
  return(APIStatus(ok));
}

static PyObject *CmdReinitialize(PyObject *dummy, PyObject *args)
{
  int ok=true;
  if (ok) {
    APIEntry();
    ok = ExecutiveReinitialize();
    APIExit();
  }
  return(APIStatus(ok));

}
static PyObject *CmdSpectrum(PyObject *self, PyObject *args)
{
  char *str1,*expr,*prefix;
  OrthoLineType s1;
  float min,max;
  int digits,start,stop, byres;
  int quiet;
  int ok=false;
  float min_ret,max_ret;
  PyObject *result = Py_None;
  ok = PyArg_ParseTuple(args,"ssffiisiii",&str1,&expr,
                        &min,&max,&start,&stop,&prefix,
                        &digits,&byres,&quiet);
  if (ok) {
    APIEntry();
    if(str1[0])
      SelectorGetTmp(str1,s1);
    else
      s1[0]=0; /* no selection */
    ok = ExecutiveSpectrum(s1,expr,min,max,start,stop,prefix,digits,byres,quiet,
                           &min_ret,&max_ret);
    if(str1[0])
      SelectorFreeTmp(s1);
    APIExit();
    if(ok) {
      result = Py_BuildValue("ff",min_ret,max_ret);
    }
  }
  return(APIAutoNone(result));
}

static PyObject *CmdMDump(PyObject *self, PyObject *args)
{
  int ok=true;
  if(ok) {
    APIEntry();
    MovieDump();
    APIExit();
  }
  return(APIStatus(ok));
}

static PyObject *CmdAccept(PyObject *self,PyObject *args)
{
  int ok=true;
  if(ok) {
    APIEntry();
    MovieSetLock(false);
    PRINTFB(FB_Movie,FB_Actions)
      " Movie: Risk accepted by user.  Movie commands have been enabled.\n"
      ENDFB;
    APIExit();
  }
  return(APIStatus(ok));
  
}

static PyObject *CmdDecline(PyObject *self,PyObject *args)
{
  int ok=true;
  if(ok) {
    APIEntry();
    MovieReset();
    PRINTFB(FB_Movie,FB_Actions)
      " Movie: Risk declined by user.  Movie commands have been deleted.\n"
      ENDFB;
    APIExit();
  }
  return(APIStatus(ok));
  
}

static PyObject *CmdSetCrystal(PyObject *self,PyObject *args)
{
  int ok=true;
  char *str1,*str2;
  OrthoLineType s1;
  float a,b,c,alpha,beta,gamma;

  ok = PyArg_ParseTuple(args,"sffffffs",&str1,&a,&b,&c,
                        &alpha,&beta,&gamma,&str2);
  if(ok) {
    SelectorGetTmp(str1,s1);
    APIEntry();
    ok = ExecutiveSetCrystal(s1,a,b,c,alpha,beta,gamma,str2);
    APIExit();
    SelectorFreeTmp(s1);
  }
  return(APIStatus(ok));
}

static PyObject *CmdGetCrystal(PyObject *self,PyObject *args)
{
  int ok=true;
  char *str1;    
  OrthoLineType s1;
  float a,b,c,alpha,beta,gamma;
  WordType sg;
  PyObject *result = NULL;
  int defined;
  ok = PyArg_ParseTuple(args,"s",&str1);
  if(ok) {
    SelectorGetTmp(str1,s1);
    APIEntry();
    ok = ExecutiveGetCrystal(s1,&a,&b,&c,&alpha,&beta,&gamma,sg,&defined);
    APIExit();
    if(ok) {
      if(defined) {
        result = PyList_New(7);
        if(result) {
          PyList_SetItem(result,0,PyFloat_FromDouble(a));
          PyList_SetItem(result,1,PyFloat_FromDouble(b));
          PyList_SetItem(result,2,PyFloat_FromDouble(c));
          PyList_SetItem(result,3,PyFloat_FromDouble(alpha));
          PyList_SetItem(result,4,PyFloat_FromDouble(beta));
          PyList_SetItem(result,5,PyFloat_FromDouble(gamma));
          PyList_SetItem(result,6,PyString_FromString(sg));
        }
      } else { /* no symmetry defined, then return empty list */
        result = PyList_New(0);
      }
    }
    SelectorFreeTmp(s1);
  }
  return(APIAutoNone(result));
}

static PyObject *CmdSmooth(PyObject *self,PyObject *args)
{
  int ok=true;
  char *str1;
  OrthoLineType s1;
  int int1,int2,int3,int4,int5,int6;
  ok = PyArg_ParseTuple(args,"siiiiii",&str1,&int1,&int2,&int3,&int4,&int5,&int6);
  if(ok) {
    SelectorGetTmp(str1,s1);
    APIEntry();
    ok = ExecutiveSmooth(s1,int1,int2,int3,int4,int5,int6);
    APIExit();
    SelectorFreeTmp(s1);
  }
  return(APIStatus(ok));
}

static PyObject *CmdGetSession(PyObject *self, PyObject *args)
{
  int ok=true;
  PyObject *dict;

  ok = PyArg_ParseTuple(args,"O",&dict);
  if(ok) {
    APIEntry();
    PBlock();
    ok = ExecutiveGetSession(dict);
    PUnblock();
    APIExit();
  }
  return(APIStatus(ok));
}

static PyObject *CmdSetSession(PyObject *self, PyObject *args)
{
  int ok=true;
  PyObject *obj;

  ok = PyArg_ParseTuple(args,"O",&obj);
  if(ok) {
    APIEntry();
    PBlock();
    ok = ExecutiveSetSession(obj);
    PUnblock();
    APIExit();
  }
  return(APIStatus(ok));
}

static PyObject *CmdSetName(PyObject *self, PyObject *args)
{
  int ok=true;
  char *str1,*str2;
  ok = PyArg_ParseTuple(args,"ss",&str1,&str2);
  if(ok) {
    APIEntry();
    ok = ExecutiveSetName(str1,str2);
    APIExit();
  }
  return(APIStatus(ok));
}

static PyObject *CmdGetBondPrint(PyObject *self,PyObject *args)
{
  int ok=true;
  char *str1;
  int ***array = NULL;
  PyObject *result = NULL;
  int int1,int2;
  int dim[3];
  ok = PyArg_ParseTuple(args,"sii",&str1,&int1,&int2);
  if(ok) {
    APIEntry();
    array = ExecutiveGetBondPrint(str1,int1,int2,dim);
    APIExit();
    if(array) {
      result = PConv3DIntArrayTo3DPyList(array,dim);
      FreeP(array);
    }
  }
  return(APIAutoNone(result));
}

static PyObject *CmdDebug(PyObject *self,PyObject *args)
{
  int ok=true;
  char *str1;
  ok = PyArg_ParseTuple(args,"s",&str1);
  if(ok) {
    APIEntry();
    ok = ExecutiveDebug(str1);
    APIExit();
  }
  return(APIStatus(ok));
}

static PyObject *CmdPGlutGetRedisplay(PyObject *self, PyObject *args)
{
  #ifdef _PYMOL_PRETEND_GLUT

  return(APIStatus(p_glutGetRedisplay()));  
  
  #else
  
  return(APIStatus(0));  

  #endif

}

static PyObject *CmdPGlutEvent(PyObject *self, PyObject *args)
{
  int ok=true;
#ifdef _PYMOL_PRETEND_GLUT
  p_glut_event ev;
  ok = PyArg_ParseTuple(args,"iiiiii",&ev.event_code,
                        &ev.x,&ev.y,&ev.input,&ev.state,&ev.mod);
   if(ok) {
    PUnblock();
    p_glutHandleEvent(&ev);
    PBlock();
  }
#endif

  return(APIStatus(ok));
}

static PyObject *CmdSculptPurge(PyObject *self, PyObject *args)
{
  int ok=true;
  if(ok) {
    APIEntry();
    SculptCachePurge();
    APIExit();
  }
  return(APIStatus(ok));
}

static PyObject *CmdSculptDeactivate(PyObject *self, PyObject *args)
{
  int ok=true;
  char *str1;
  ok = PyArg_ParseTuple(args,"s",&str1);
  if(ok) {
    APIEntry();
    ok = ExecutiveSculptDeactivate(str1);
    APIExit();
  }
  return(APIStatus(ok));
}

static PyObject *CmdSculptActivate(PyObject *self, PyObject *args)
{
  int ok=true;
  int int1;
  char *str1;
  ok = PyArg_ParseTuple(args,"si",&str1,&int1);
  if(ok) {
    APIEntry();
    ok = ExecutiveSculptActivate(str1,int1);
    APIExit();
  }
  return(APIStatus(ok));
}

static PyObject *CmdSculptIterate(PyObject *self, PyObject *args)
{
  int ok=true;
  int int1,int2;
  char *str1;
  float total_strain = 0.0F;
  ok = PyArg_ParseTuple(args,"sii",&str1,&int1,&int2);
  if(ok) {
    APIEntry();
    total_strain = ExecutiveSculptIterate(str1,int1,int2);
    APIExit();
  }
  return(APIIncRef(PyFloat_FromDouble((double)total_strain)));
}

static PyObject *CmdCombineObjectTTT(PyObject *self, 	PyObject *args)
{
  char *name;
  PyObject *m;
  float ttt[16];
  int ok = false;
  ok = PyArg_ParseTuple(args,"sO",&name,&m);
  if(ok) {
    if(PConvPyListToFloatArrayInPlace(m,ttt,16)>0) {
      APIEntry();
      ok = ExecutiveCombineObjectTTT(name,ttt);
      APIExit();
    } else {
      PRINTFB(FB_CCmd,FB_Errors)
        "CmdCombineObjectTTT-Error: bad matrix\n"
        ENDFB;
      ok=false;
    }
  }
  return(APIStatus(ok));
}

static PyObject *CmdGetColor(PyObject *self, PyObject *args)
{
  char *name;
  int mode;
  int ok = false;
  int a,nc,nvc;
  float *rgb;
  int index;
  PyObject *result = NULL;
  PyObject *tup;
  ok = PyArg_ParseTuple(args,"si",&name,&mode);
  if(ok) {
    APIEntry();
    switch(mode) {
    case 0: /* by name or index, return floats */
      index = ColorGetIndex(name);
      if(index>=0) {
        rgb = ColorGet(index);
        tup = PyTuple_New(3);
        PyTuple_SetItem(tup,0,PyFloat_FromDouble(*(rgb++)));
        PyTuple_SetItem(tup,1,PyFloat_FromDouble(*(rgb++)));
        PyTuple_SetItem(tup,2,PyFloat_FromDouble(*rgb));
        result=tup;
      }
      break;
    case 1: /* get color names with NO NUMBERS in their names */
      PBlock();
      nc=ColorGetNColor();
      nvc=0;
      for(a=0;a<nc;a++) {
        if(ColorGetStatus(a))
          nvc++;
      }
      result = PyList_New(nvc);
      nvc=0;
      for(a=0;a<nc;a++) {
        if(ColorGetStatus(a)) {
          tup = PyTuple_New(2);
          PyTuple_SetItem(tup,0,PyString_FromString(ColorGetName(a)));
          PyTuple_SetItem(tup,1,PyInt_FromLong(a));
          PyList_SetItem(result,nvc++,tup);
        }
      }
      PUnblock();
      break;
    }
    APIExit();
  }
  return(APIAutoNone(result));
}

static PyObject *CmdGetChains(PyObject *self, PyObject *args)
{

  char *str1;
  int int1;
  OrthoLineType s1="";
  PyObject *result = NULL;
  char *chain_str = NULL;
  int ok=false;
  int c1=0;
  int a,l;
  int null_chain = false;
  ok = PyArg_ParseTuple(args,"si",&str1,&int1);
  if(ok) {
    APIEntry();
    if(str1[0]) c1 = SelectorGetTmp(str1,s1);
    if(c1)
      chain_str = ExecutiveGetChains(s1,int1,&null_chain);
    
    if(chain_str) {
      l=strlen(chain_str);
      if(null_chain) l++; 
      result = PyList_New(l);
      if(null_chain) {
        l--;
        PyList_SetItem(result,l,PyString_FromString(""));
      }
      for(a=0;a<l;a++)
        PyList_SetItem(result,a,PyString_FromStringAndSize(chain_str+a,1));

    }
    FreeP(chain_str);
    if(s1[0]) SelectorFreeTmp(s1);
    APIExit();
  }
  return(APIAutoNone(result));

}

static PyObject *CmdMultiSave(PyObject *self, PyObject *args)
{
  char *name,*object;
  int append,state;
  int ok = false;
  ok = PyArg_ParseTuple(args,"ssii",&name,&object,&state,&append);
  if(ok) {
    APIEntry();
    ok = ExecutiveMultiSave(name,object,state,append);
    APIExit();
  }
  return(APIStatus(ok));
}

static PyObject *CmdRampMapNew(PyObject *self, 	PyObject *args)
{
  char *name;
  int ok = false;
  char *map;
  int map_state;
  PyObject *range,*color;
  ok = PyArg_ParseTuple(args,"ssOOi",&name,&map,&range,&color,&map_state);
  if(ok) {
    APIEntry();
    ok = ExecutiveRampMapNew(name,map,range,color,map_state);
    APIExit();
  }
  return(APIStatus(ok));
}

static PyObject *CmdMapNew(PyObject *self, PyObject *args)
{
  char *name;
  float minCorner[3],maxCorner[3];
  float grid[3];
  float buffer;
  int type;
  int state;
  char *selection;
  OrthoLineType s1;
  int ok = false;
  ok = PyArg_ParseTuple(args,"sifsf(ffffff)i",&name,&type,&grid[0],&selection,&buffer,
                        &minCorner[0],&minCorner[1],&minCorner[2],
                        &maxCorner[0],&maxCorner[1],&maxCorner[2],&state);
  if(ok) {
    grid[1]=grid[0];
    grid[2]=grid[0];
    APIEntry();
    SelectorGetTmp(selection,s1);
    ok = ExecutiveMapNew(name,type,grid,s1,buffer,minCorner,maxCorner,state);
    SelectorFreeTmp(s1);
    APIExit();
  }
  return(APIStatus(ok));
}

static PyObject *CmdMapSetBorder(PyObject *self, PyObject *args)
{
  char *name;
  float level;
  int ok = false;
  ok = PyArg_ParseTuple(args,"sf",&name,&level);
  if(ok) {
    APIEntry();
    ok = ExecutiveMapSetBorder(name,level);
    APIExit();
  }
  return(APIStatus(ok));
}

static PyObject *CmdMapDouble(PyObject *self, PyObject *args)
{
  char *name;
  int state;
  int ok = false;
  ok = PyArg_ParseTuple(args,"si",&name,&state);
  if(ok) {
    APIEntry();
    ok = ExecutiveMapDouble(name,state);
    APIExit();
  }
  return(APIStatus(ok));
}

static PyObject *CmdGetRenderer(PyObject *self, PyObject *args)
{
  char *vendor,*renderer,*version;
  APIEntry();
  SceneGetCardInfo(&vendor,&renderer,&version);
  APIExit();
  return Py_BuildValue("(sss)",vendor,renderer,version);
}

static PyObject *CmdTranslateAtom(PyObject *self, PyObject *args)
{
  char *str1;
  int state,log,mode;
  float v[3];
  OrthoLineType s1;
  int ok = false;
  ok = PyArg_ParseTuple(args,"sfffiii",&str1,v,v+1,v+2,&state,&mode,&log);
  if(ok) {
    SelectorGetTmp(str1,s1);
    APIEntry();
    ok = ExecutiveTranslateAtom(s1,v,state,mode,log);
    APIExit();
    SelectorFreeTmp(s1);
  }
  return(APIStatus(ok));
}

static PyObject *CmdTransformObject(PyObject *self, PyObject *args)
{
  char *name,*sele;
  int state,log;
  PyObject *m;
  float ttt[16];
  int ok = false;
  ok = PyArg_ParseTuple(args,"siOis",&name,&state,&m,&log,&sele);
  if(ok) {
    if(PConvPyListToFloatArrayInPlace(m,ttt,16)>0) {
      APIEntry();
      ok = ExecutiveTransformObjectSelection(name,state,sele,log,ttt);
      APIExit();
    } else {
      PRINTFB(FB_CCmd,FB_Errors)
        "CmdTransformObject-DEBUG: bad matrix\n"
        ENDFB;
      ok=false;
    }
  }
  return(APIStatus(ok));
}

static PyObject *CmdTransformSelection(PyObject *self, PyObject *args)
{
  char *sele;
  int state,log;
  PyObject *m;
  float ttt[16];
  OrthoLineType s1;
  int ok = false;
  ok = PyArg_ParseTuple(args,"siOi",&sele,&state,&m,&log);
  if(ok) {
    if(PConvPyListToFloatArrayInPlace(m,ttt,16)>0) {
      APIEntry();
      SelectorGetTmp(sele,s1);
      ok = ExecutiveTransformSelection(state,s1,log,ttt);
      SelectorFreeTmp(s1);
      APIExit();
    } else {
      PRINTFB(FB_CCmd,FB_Errors)
        "CmdTransformSelection-DEBUG: bad matrix\n"
        ENDFB;
      ok=false;
    }
  }
  return(APIStatus(ok));
}

static PyObject *CmdLoadColorTable(PyObject *self, PyObject *args)
{
  char *str1;
  int ok = false;
  int quiet;
  ok = PyArg_ParseTuple(args,"si",&str1,&quiet);
  if(ok) {
    APIEntry();
    ok = ColorTableLoad(str1,quiet);
    APIExit();
  }
  return(APIStatus(ok));
}

static PyObject *CmdLoadPNG(PyObject *self, PyObject *args)
{
  char *str1;
  int ok = false;
  int quiet;
  int movie;
  ok = PyArg_ParseTuple(args,"sii",&str1,&movie,&quiet);
  if(ok) {
    APIEntry();
    ok = SceneLoadPNG(str1,movie,quiet);
    APIExit();
  }
  return(APIStatus(ok));
}

static PyObject *CmdBackgroundColor(PyObject *self, PyObject *args)
{
  char *str1;
  int ok = false;
  int idx;
  ok = PyArg_ParseTuple(args,"s",&str1);
  if(ok) {
    APIEntry();
    idx = ColorGetIndex(str1);
    if(idx>=0)
      ok = SettingSetfv(cSetting_bg_rgb,ColorGet(idx));
    else {
      ErrMessage("Color","Bad color name.");
      ok = false; /* bad color */
    }
    APIExit();
  }
  return(APIStatus(ok));
}

static PyObject *CmdGetPosition(PyObject *self, 	PyObject *args)
{
  PyObject *result;
  float v[3];
  APIEntry();
  SceneGetPos(v);
  APIExit();
  result=PConvFloatArrayToPyList(v,3);
  return(result);
}

static PyObject *CmdGetPhiPsi(PyObject *self, 	PyObject *args)
{
  char *str1;
  OrthoLineType s1;
  int state;
  PyObject *result = Py_None;
  PyObject *key = Py_None;
  PyObject *value = Py_None;
  int *iVLA=NULL;
  float *pVLA,*sVLA;
  int l;
  int *i;
  ObjectMolecule **o,**oVLA=NULL;
  int a;
  float *s,*p;
  int ok =  PyArg_ParseTuple(args,"si",&str1,&state);
  if(ok) {
    APIEntry();
    SelectorGetTmp(str1,s1);
    l = ExecutivePhiPsi(s1,&oVLA,&iVLA,&pVLA,&sVLA,state);
    SelectorFreeTmp(s1);
    APIExit();
    if(iVLA) {
      result=PyDict_New();
      i = iVLA;
      o = oVLA;
      p = pVLA;
      s = sVLA;
      for(a=0;a<l;a++) {
        key = PyTuple_New(2);      
        PyTuple_SetItem(key,1,PyInt_FromLong(*(i++)+1)); /* +1 for index */
        PyTuple_SetItem(key,0,PyString_FromString((*(o++))->Obj.Name));
        value = PyTuple_New(2);      
        PyTuple_SetItem(value,0,PyFloat_FromDouble(*(p++))); /* +1 for index */
        PyTuple_SetItem(value,1,PyFloat_FromDouble(*(s++)));
        PyDict_SetItem(result,key,value);
        Py_DECREF(key);
        Py_DECREF(value);
      }
    } else {
      result = PyDict_New();
    }
    VLAFreeP(iVLA);
    VLAFreeP(oVLA);
    VLAFreeP(sVLA);
    VLAFreeP(pVLA);
  }
  return(APIAutoNone(result));
}

static PyObject *CmdAlign(PyObject *self, 	PyObject *args) {
  char *str2,*str3,*mfile,*oname;
  OrthoLineType s2="",s3="";
  float result = -1.0;
  int ok = false;
  int quiet,cycles,skip;
  float cutoff,gap,extend;
  int state1,state2;
  ok = PyArg_ParseTuple(args,"ssfiffissiii",&str2,&str3,
                        &cutoff,&cycles,&gap,&extend,&skip,&oname,
                        &mfile,&state1,&state2,&quiet);

  if(ok) {
    PRINTFD(FB_CCmd)
      "CmdAlign-DEBUG %s %s\n",
      str2,str3
      ENDFD;
    
    APIEntry();
    SelectorGetTmp(str2,s2);
    SelectorGetTmp(str3,s3);
    result = ExecutiveAlign(s2,s3,mfile,gap,extend,skip,cutoff,
                            cycles,quiet,oname,state1,state2);
    SelectorFreeTmp(s2);
    SelectorFreeTmp(s3);
    APIExit();
  }
  return Py_BuildValue("f",result);
}

static PyObject *CmdGetSettingUpdates(PyObject *self, 	PyObject *args)
{
  PyObject *result = NULL;
  APIEnterBlocked();
  result = SettingGetUpdateList(NULL);
  APIExitBlocked();
  return(APIAutoNone(result));
}

static PyObject *CmdGetView(PyObject *self, 	PyObject *args)
{
  SceneViewType view;
  APIEntry();
  SceneGetView(view);
  APIExit();
  return(Py_BuildValue("(fffffffffffffffffffffffff)",
                   view[ 0],view[ 1],view[ 2],view[ 3], /* 4x4 mat */
                   view[ 4],view[ 5],view[ 6],view[ 7],
                   view[ 8],view[ 9],view[10],view[11],
                   view[12],view[13],view[14],view[15],
                   view[16],view[17],view[18], /* pos */
                   view[19],view[20],view[21], /* origin */
                   view[22],view[23], /* clip */
                   view[24] /* orthoscopic*/
                       ));
}
static PyObject *CmdSetView(PyObject *self, 	PyObject *args)
{
  SceneViewType view;
  int quiet;
  int ok=PyArg_ParseTuple(args,"(fffffffffffffffffffffffff)i",
                   &view[ 0],&view[ 1],&view[ 2],&view[ 3], /* 4x4 mat */
                   &view[ 4],&view[ 5],&view[ 6],&view[ 7],
                   &view[ 8],&view[ 9],&view[10],&view[11],
                   &view[12],&view[13],&view[14],&view[15],
                   &view[16],&view[17],&view[18], /* pos */
                   &view[19],&view[20],&view[21], /* origin */
                   &view[22],&view[23], /* clip */
                   &view[24], /* orthoscopic*/
                   &quiet);
  if(ok) {
    APIEntry();
    SceneSetView(view,quiet); /* TODO STATUS */
    APIExit();
  }
  return(APIStatus(ok));
}
static PyObject *CmdGetState(PyObject *self, 	PyObject *args)
{
  return(APIStatus(SceneGetState()));
}

static PyObject *CmdGetFrame(PyObject *self, 	PyObject *args)
{
  return(APIStatus(SceneGetFrame()+1));
}

static PyObject *CmdSetTitle(PyObject *self, PyObject *args)
{
  char *str1,*str2;
  int int1;
  int ok = false;
  ok = PyArg_ParseTuple(args,"sis",&str1,&int1,&str2);
  if(ok) {
    APIEntry();
    ok = ExecutiveSetTitle(str1,int1,str2);
    APIExit();
  }
  return(APIStatus(ok));
}

static PyObject *CmdGetTitle(PyObject *self, PyObject *args)
{
  char *str1,*str2;
  int int1;
  int ok = false;
  PyObject *result = Py_None;
  ok = PyArg_ParseTuple(args,"si",&str1,&int1);
  if(ok) {
    APIEntry();
    str2 = ExecutiveGetTitle(str1,int1);
    if(str2) result = PyString_FromString(str2);
    APIExit();
  }
  return(APIAutoNone(result));
}


static PyObject *CmdExportCoords(PyObject *self, 	PyObject *args)
{
  void *result;
  char *str1;
  int int1;
  PyObject *py_result = Py_None;
  int ok = false;
  ok = PyArg_ParseTuple(args,"si",&str1,&int1);
  if(ok) {
    APIEntry();
    result = ExportCoordsExport(str1,int1,0);
    APIExit();
    if(result) 
      py_result = PyCObject_FromVoidPtr(result,(void(*)(void*))ExportCoordsFree);
  }
  return(APIAutoNone(py_result));
}

static PyObject *CmdImportCoords(PyObject *self, 	PyObject *args)
{
  char *str1;
  int int1;
  PyObject *cObj;
  void *mmdat=NULL;
  int ok = false;
  ok = PyArg_ParseTuple(args,"siO",&str1,&int1,&cObj);
  if(ok) {
    if(PyCObject_Check(cObj))
      mmdat = PyCObject_AsVoidPtr(cObj);
    APIEntry();
    if(mmdat)
      ok = ExportCoordsImport(str1,int1,mmdat,0);
    APIExit();
  }
  return(APIStatus(ok));
}

static PyObject *CmdGetArea(PyObject *self, 	PyObject *args)
{
  char *str1;
  int int1,int2;
  OrthoLineType s1="";
  float result = -1.0;
  int ok=false;
  int c1=0;
  ok = PyArg_ParseTuple(args,"sii",&str1,&int1,&int2);
  if(ok) {
    APIEntry();
    if(str1[0]) c1 = SelectorGetTmp(str1,s1);
    if(c1)
      result = ExecutiveGetArea(s1,int1,int2);
    else
      result = 0.0; /* empty selection */
    if(s1[0]) SelectorFreeTmp(s1);
    APIExit();
  }
  return(Py_BuildValue("f",result));

}

static PyObject *CmdPushUndo(PyObject *self, 	PyObject *args)
{
  char *str0;
  int state;
  OrthoLineType s0="";
  int ok=false;
  ok = PyArg_ParseTuple(args,"si",&str0,&state);
  if(ok) {
    APIEntry();
    if(str0[0]) SelectorGetTmp(str0,s0);
    ok = ExecutiveSaveUndo(s0,state);
    if(s0[0]) SelectorFreeTmp(s0);
    APIExit();
  }
  return(APIStatus(ok));
}

static PyObject *CmdGetType(PyObject *self, 	PyObject *args)
{
  char *str1;
  WordType type = "";
  int ok=false;
  ok = PyArg_ParseTuple(args,"s",&str1);
  if (ok) {
    APIEntry();
    ok = ExecutiveGetType(str1,type);
    APIExit();
  } 
  if(ok) 
    return(Py_BuildValue("s",type));
  else
    return(APIStatus(ok));
}
static PyObject *CmdGetNames(PyObject *self, 	PyObject *args)
{
  int int1,int2;
  char *vla = NULL;
  PyObject *result = Py_None;
  int ok=false;
  ok = PyArg_ParseTuple(args,"ii",&int1,&int2);
  if(ok) {
    APIEntry();
    vla = ExecutiveGetNames(int1,int2);
    APIExit();
    result = PConvStringVLAToPyList(vla);
    VLAFreeP(vla);
  }
  return(APIAutoNone(result));
}

static PyObject *CmdInvert(PyObject *self, PyObject *args)
{
  int int1;
  int ok=false;
  ok = PyArg_ParseTuple(args,"i",&int1);
  if(ok) {
    APIEntry();
    ok = ExecutiveInvert(int1);
    APIExit();
  }
  return(APIStatus(ok));
}

static PyObject *CmdTorsion(PyObject *self, PyObject *args)
{
  float float1;
  int ok=false;
  ok = PyArg_ParseTuple(args,"f",&float1);
  if (ok) {
    APIEntry();
    ok = EditorTorsion(float1);
    APIExit();
  }
  return(APIStatus(ok));
}

static PyObject *CmdUndo(PyObject *self, PyObject *args)
{
  int int1;
  int ok=false;
  ok = PyArg_ParseTuple(args,"i",&int1);
  if (ok) {
    APIEntry();
    ExecutiveUndo(int1); /* TODO STATUS */
    APIExit();
  }
  return(APIStatus(ok));
}

static PyObject *CmdMask(PyObject *self, PyObject *args)
{
  char *str1;
  int int1;
  OrthoLineType s1;

  int ok=false;
  ok = PyArg_ParseTuple(args,"si",&str1,&int1);
  if (ok) {
    APIEntry();
    SelectorGetTmp(str1,s1);
    ExecutiveMask(s1,int1); /* TODO STATUS */
    SelectorFreeTmp(s1);
    APIExit();
  }
  return(APIStatus(ok));
}

static PyObject *CmdProtect(PyObject *self, PyObject *args)
{
  char *str1;
  int int1,int2;
  OrthoLineType s1;

  int ok=false;
  ok = PyArg_ParseTuple(args,"sii",&str1,&int1,&int2);
  if (ok) {
    APIEntry();
    SelectorGetTmp(str1,s1);
    ExecutiveProtect(s1,int1,int2); /* TODO STATUS */
    SelectorFreeTmp(s1);
    APIExit();
  }
  return(APIStatus(ok));
}


static PyObject *CmdButton(PyObject *self, 	PyObject *args)
{
  int i1,i2;
  int ok=false;
  ok = PyArg_ParseTuple(args,"ii",&i1,&i2);
  if (ok) {
    APIEntry();
    ButModeSet(i1,i2); /* TODO STATUS */
    APIExit();
  }
  return(APIStatus(ok));  
}

static PyObject *CmdFeedback(PyObject *self, 	PyObject *args)
{
  int i1,i2,result = 0;
  int ok=false;
  ok = PyArg_ParseTuple(args,"ii",&i1,&i2);
  if (ok) {
    /* NO API Entry for performance,
     *feedback (MACRO) just accesses a safe global */
    result = Feedback(i1,i2);
  }
  return Py_BuildValue("i",result);
}

static PyObject *CmdSetFeedbackMask(PyObject *self, 	PyObject *args)
{
  int i1,i2,i3;
  int ok=false;
  ok = PyArg_ParseTuple(args,"iii",&i1,&i2,&i3);
  if (ok) {
    APIEntry();
    switch(i1) { /* TODO STATUS */
    case 0: 
      FeedbackSetMask(i2,(uchar)i3);
      break;
    case 1:
      FeedbackEnable(i2,(uchar)i3);
      break;
    case 2:
      FeedbackDisable(i2,(uchar)i3);
      break;
    case 3:
      FeedbackPush();
      break;
    case 4:
      FeedbackPop();
      break;
    }
    APIExit();
  }
  return(APIStatus(ok));  
}

static PyObject *CmdFocus(PyObject *self, 	PyObject *args)
{
  /* BROKEN */
  APIEntry();
  ExecutiveFocus();
  APIExit();
  return(APIFailure());
}

static PyObject *CmdFlushNow(PyObject *self, 	PyObject *args)
{
  /* only called by the GLUT thread with unlocked API, blocked interpreter */
  if(flush_count<8) { /* prevent super-deep recursion */
    flush_count++;
    PFlushFast();
    flush_count--;
  } else {
    PRINTFB(FB_CCmd,FB_Warnings)
      " Cmd: PyMOL lagging behind API requests...\n"
      ENDFB;
  }
  return(APISuccess());  
}

static PyObject *CmdWaitQueue(PyObject *self, 	PyObject *args)
{
  /* called by non-GLUT thread with unlocked API, blocked interpreter */
  PyObject *result = NULL;
  if(!PyMOLTerminating) {
    APIEnterBlocked();
    if(OrthoCommandWaiting()||(flush_count>1)) 
      result = PyInt_FromLong(1);
    else
      result = PyInt_FromLong(0);
    APIExitBlocked();
  }
  return APIAutoNone(result);
}

static PyObject *CmdPaste(PyObject *dummy, PyObject *args)
{
  PyObject *list,*str;
  char *st;
  int l,a;
  int ok=false;
  ok = PyArg_ParseTuple(args,"O",&list);
  if(ok) {
    if(!list) 
      ok=false;
    else if(!PyList_Check(list)) 
      ok=false;
    else
      {
        l=PyList_Size(list);
        for(a=0;a<l;a++) {
          str = PyList_GetItem(list,a);
          if(str) {
            if(PyString_Check(str)) {
              st = PyString_AsString(str);
              APIEntry();
              OrthoPasteIn(st);
              APIExit();
            } else {
              ok = false;
            }
          }
        }
      }
  }
  return(APIStatus(ok));  
}

static PyObject *CmdGetPovRay(PyObject *dummy, PyObject *args)
{
  PyObject *result = NULL;
  char *header=NULL,*geom=NULL;
  APIEntry();
  SceneRay(0,0,1,&header,&geom,0.0F,0.0F,false);
  if(header&&geom) {
    result = Py_BuildValue("(ss)",header,geom);
  }
  VLAFreeP(header);
  VLAFreeP(geom);
  APIExit();
  return(APIAutoNone(result));
}

static PyObject *CmdGetWizard(PyObject *dummy, PyObject *args)
{
  PyObject *result;
  APIEntry();
  result = WizardGet();
  APIExit();
  if(!result)
    result=Py_None;
  return APIIncRef(result);
}

static PyObject *CmdGetWizardStack(PyObject *dummy, PyObject *args)
{
  PyObject *result;
  APIEnterBlocked();
  result = WizardGetStack();
  APIExitBlocked();
  if(!result)
    result=Py_None;
  return APIIncRef(result);
}

static PyObject *CmdSetWizard(PyObject *dummy, PyObject *args)
{
  
  PyObject *obj;
  int ok=false;
  int replace;
  ok = PyArg_ParseTuple(args,"Oi",&obj,&replace);
  if(ok) {
    if(!obj)
      ok=false;
    else
      {
        APIEntry();
        WizardSet(obj,replace); /* TODO STATUS */
        APIExit();
      }
  }
  return(APIStatus(ok));  
}

static PyObject *CmdSetWizardStack(PyObject *dummy, PyObject *args)
{
  
  PyObject *obj;
  int ok=false;
  ok = PyArg_ParseTuple(args,"O",&obj);
  if(ok) {
    if(!obj)
      ok=false;
    else
      {
        APIEntry();
        WizardSetStack(obj); /* TODO STATUS */
        APIExit();
      }
  }
  return(APIStatus(ok));  
}

static PyObject *CmdRefreshWizard(PyObject *dummy, PyObject *args)
{
  
  APIEntry();
  WizardRefresh();
  OrthoDirty();
  APIExit();
  return(APISuccess());  
}

static PyObject *CmdSplash(PyObject *dummy, PyObject *args)
{
  APIEntry();
  OrthoSplash();
  APIExit();
  return(APISuccess());  
}

static PyObject *CmdCls(PyObject *dummy, PyObject *args)
{
  APIEntry();
  OrthoClear();
  APIExit();
  return(APISuccess());  
}

static PyObject *CmdDump(PyObject *dummy, PyObject *args)
{
  char *str1,*str2;
  int ok=false;
  ok = PyArg_ParseTuple(args,"ss",&str1,&str2);
  if (ok) {
    APIEntry();
    ExecutiveDump(str1,str2); /* TODO STATUS */
    APIExit();
  }
  return(APIStatus(ok));  
}

static PyObject *CmdIsomesh(PyObject *self, 	PyObject *args) {
  char *str1,*str2,*str3;
  float lvl,fbuf;
  int dotFlag;
  int c,state=-1;
  OrthoLineType s1;
  int oper,frame;
  float carve;
  CObject *obj=NULL,*mObj,*origObj;
  ObjectMap *mapObj;
  float mn[3] = { 0,0,0};
  float mx[3] = { 15,15,15};
  float *vert_vla = NULL;
  int ok = false;
  int map_state;
  int multi=false;
  ObjectMapState *ms;

  /* oper 0 = all, 1 = sele + buffer, 2 = vector */

  ok = PyArg_ParseTuple(args,"sisisffiifi",&str1,&frame,&str2,&oper,
                   &str3,&fbuf,&lvl,&dotFlag,&state,&carve,&map_state);
  if (ok) {
    APIEntry();

    origObj=ExecutiveFindObjectByName(str1);  
    if(origObj) {
      if(origObj->type!=cObjectMesh) {
        ExecutiveDelete(str1);
        origObj=NULL;
      }
    }
    
    mObj=ExecutiveFindObjectByName(str2);  
    if(mObj) {
      if(mObj->type!=cObjectMap)
        mObj=NULL;
    }
    if(mObj) {
      mapObj = (ObjectMap*)mObj;
      if(state==-1) {
        multi=true;
        state=0;
        map_state=0;
      } else if(state==-2) {
        state=SceneGetState();
        if(map_state<0) 
          map_state=state;
      } else if(state==-3) { /* append mode */
        state=0;
        if(origObj)
          if(origObj->fGetNFrame)
            state=origObj->fGetNFrame(origObj);
      } else {
        if(map_state==-1) {
          map_state=0;
          multi=true;
        } else {
          multi=false;
        }
      }
      while(1) {
        if(map_state==-2)
          map_state=SceneGetState();
        if(map_state==-3)
          map_state=ObjectMapGetNStates(mapObj)-1;
        ms = ObjectMapStateGetActive(mapObj,map_state);
        if(ms) {
          switch(oper) {
          case 0:
            for(c=0;c<3;c++) {
              mn[c] = ms->Corner[0][c];
              mx[c] = ms->Corner[7][c];
            }
            carve = -1.0; /* impossible */
            break;
          case 1:
            SelectorGetTmp(str3,s1);
            ExecutiveGetExtent(s1,mn,mx,false,-1,false); /* TODO state */
            if(carve>=0.0) {
              vert_vla = ExecutiveGetVertexVLA(s1,state);
              if(fbuf<=R_SMALL4)
                fbuf = carve;
            }
            SelectorFreeTmp(s1);
            for(c=0;c<3;c++) {
              mn[c]-=fbuf;
              mx[c]+=fbuf;
            }
            break;
          }
          PRINTFB(FB_CCmd,FB_Blather)
            " Isomesh: buffer %8.3f carve %8.3f \n",fbuf,carve
            ENDFB;
          obj=(CObject*)ObjectMeshFromBox((ObjectMesh*)origObj,mapObj,map_state,state,mn,mx,lvl,dotFlag,
                                          carve,vert_vla);
          if(!origObj) {
            ObjectSetName(obj,str1);
            ExecutiveManageObject((CObject*)obj,true,false);
          }
          
          if(SettingGet(cSetting_isomesh_auto_state))
            if(obj) ObjectGotoState((ObjectMolecule*)obj,state);
          PRINTFB(FB_ObjectMesh,FB_Actions)
            " Isomesh: created \"%s\", setting level to %5.3f\n",str1,lvl
            ENDFB;
        } else if(!multi) {
          PRINTFB(FB_ObjectMesh,FB_Warnings)
            " Isomesh-Warning: state %d not present in map \"%s\".\n",map_state+1,str2
            ENDFB;
          ok=false;
        }
        if(multi) {
          origObj = obj;
          map_state++;
          state++;
          if(map_state>=mapObj->NState)
            break;
        } else {
          break;
        }
      }
    } else {
      PRINTFB(FB_ObjectMesh,FB_Errors)
        " Isomesh: Map or brick object \"%s\" not found.\n",str2
        ENDFB;
      ok=false;
    }
    APIExit();
  }
  return(APIStatus(ok));  
}

static PyObject *CmdIsosurface(PyObject *self, 	PyObject *args) {
  char *str1,*str2,*str3;
  float lvl,fbuf;
  int dotFlag;
  int c,state=-1;
  OrthoLineType s1;
  int oper,frame;
  float carve;
  CObject *obj=NULL,*mObj,*origObj;
  ObjectMap *mapObj;
  float mn[3] = { 0,0,0};
  float mx[3] = { 15,15,15};
  float *vert_vla = NULL;
  int ok = false;
  ObjectMapState *ms;
  int map_state=0;
  int multi=false;
  int side;
  /* oper 0 = all, 1 = sele + buffer, 2 = vector */

  ok = PyArg_ParseTuple(args,"sisisffiifii",&str1,&frame,&str2,&oper,
                   &str3,&fbuf,&lvl,&dotFlag,&state,&carve,&map_state,
                        &side);
  if (ok) {
    APIEntry();

    origObj=ExecutiveFindObjectByName(str1);  
    if(origObj) {
      if(origObj->type!=cObjectSurface) {
        ExecutiveDelete(str1);
        origObj=NULL;
      }
    }
    
    mObj=ExecutiveFindObjectByName(str2);  
    if(mObj) {
      if(mObj->type!=cObjectMap)
        mObj=NULL;
    }
    if(mObj) {
      mapObj = (ObjectMap*)mObj;
      if(state==-1) {
        multi=true;
        state=0;
        map_state=0;
      } else if(state==-2) { /* current state */
        state=SceneGetState();
        if(map_state<0) 
          map_state=state;
      } else if(state==-3) { /* append mode */
        state=0;
        if(origObj)
          if(origObj->fGetNFrame)
            state=origObj->fGetNFrame(origObj);
      } else {
        if(map_state==-1) {
          map_state=0;
          multi=true;
        } else {
          multi=false;
        }
      }
      while(1) {
        if(map_state==-2)
          map_state=SceneGetState();
        if(map_state==-3)
          map_state=ObjectMapGetNStates(mapObj)-1;
        ms = ObjectMapStateGetActive(mapObj,map_state);
        if(ms) {
          switch(oper) {
          case 0:
            for(c=0;c<3;c++) {
              mn[c] = ms->Corner[0][c];
              mx[c] = ms->Corner[7][c];
            }
            carve = false; /* impossible */
            break;
          case 1:
            SelectorGetTmp(str3,s1);
            ExecutiveGetExtent(s1,mn,mx,false,-1,false); /* TODO state */
            if(carve>=0.0) {
              vert_vla = ExecutiveGetVertexVLA(s1,state);
              if(fbuf<=R_SMALL4)
                fbuf = carve;
            }
            SelectorFreeTmp(s1);
            for(c=0;c<3;c++) {
              mn[c]-=fbuf;
              mx[c]+=fbuf;
            }
            break;
          }
          PRINTFB(FB_CCmd,FB_Blather)
            " Isosurface: buffer %8.3f carve %8.3f\n",fbuf,carve
            ENDFB;
          obj=(CObject*)ObjectSurfaceFromBox((ObjectSurface*)origObj,mapObj,map_state,
                                             state,mn,mx,lvl,dotFlag,
                                             carve,vert_vla,side);
          if(!origObj) {
            ObjectSetName(obj,str1);
            ExecutiveManageObject((CObject*)obj,true,false);
          }
          if(SettingGet(cSetting_isomesh_auto_state))
            if(obj) ObjectGotoState((ObjectMolecule*)obj,state);
          PRINTFB(FB_ObjectSurface,FB_Actions)
            " Isosurface: created \"%s\", setting level to %5.3f\n",str1,lvl
            ENDFB;
        } else if(!multi) {
          PRINTFB(FB_ObjectMesh,FB_Warnings)
            " Isosurface-Warning: state %d not present in map \"%s\".\n",map_state+1,str2
            ENDFB;
          ok=false;
        }
        if(multi) {
          origObj = obj;
          map_state++;
          state++;
          if(map_state>=mapObj->NState)
            break;
        } else {
          break;
        }
      }
    } else {
      PRINTFB(FB_ObjectSurface,FB_Errors)
        " Isosurface: Map or brick object \"%s\" not found.\n",str2
        ENDFB;
      ok=false;
    }
    APIExit();
  }
  return(APIStatus(ok));  
}

static PyObject *CmdSymExp(PyObject *self, 	PyObject *args) {
  char *str1,*str2,*str3;
  OrthoLineType s1;
  float cutoff;
  CObject *mObj;
  /* oper 0 = all, 1 = sele + buffer, 2 = vector */

  int ok=false;
  ok = PyArg_ParseTuple(args,"sssf",&str1,&str2,&str3,&cutoff);
  if (ok) {
    APIEntry();
    mObj=ExecutiveFindObjectByName(str2);  
    if(mObj) {
      if(mObj->type!=cObjectMolecule) {
        mObj=NULL;
        ok = false;
      }
    }
    if(mObj) {
      SelectorGetTmp(str3,s1);
      ExecutiveSymExp(str1,str2,s1,cutoff); /* TODO STATUS */
      SelectorFreeTmp(s1);
    }
    APIExit();
  }
  return(APIStatus(ok));  
}

static PyObject *CmdOverlap(PyObject *dummy, PyObject *args)
{
  char *str1,*str2;
  int state1,state2;
  float overlap = -1.0;
  float adjust;
  OrthoLineType s1,s2;
  int ok=false;
  ok = PyArg_ParseTuple(args,"ssiif",&str1,&str2,&state1,&state2,&adjust);
  if (ok) {
    APIEntry();
    SelectorGetTmp(str1,s1);
    SelectorGetTmp(str2,s2);
    overlap = ExecutiveOverlap(s1,state1,s2,state2,adjust);
    SelectorFreeTmp(s1);
    SelectorFreeTmp(s2);
    APIExit();
  }
  return(Py_BuildValue("f",overlap));
}

static PyObject *CmdDist(PyObject *dummy, PyObject *args)
{
  char *name,*str1,*str2;
  float cutoff,result=-1.0;
  int labels,quiet;
  int mode;
  OrthoLineType s1,s2;
  int ok=false;
  int c1,c2;
  ok = PyArg_ParseTuple(args,"sssifii",&name,&str1,&str2,&mode,&cutoff,&labels,&quiet);
  if (ok) {
    APIEntry();
    c1 = SelectorGetTmp(str1,s1);
    c2 = SelectorGetTmp(str2,s2);
    if(c1&&(c2||WordMatch(cKeywordSame,s2,true)))
        result = ExecutiveDist(name,s1,s2,mode,cutoff,labels,quiet);
    else {
      if((!quiet)&&(!c1)) {
        PRINTFB(FB_Executive,FB_Errors)
          " Distance-ERR: selection 1 contains no atoms.\n"
          ENDFB;
      } 
      if((quiet!=2)&&(!c2)) {
        PRINTFB(FB_Executive,FB_Errors)
          " Distance-ERR: selection 2 contains no atoms.\n"
          ENDFB;
      }
      result = -1.0;
    }
    SelectorFreeTmp(s1);
    SelectorFreeTmp(s2);
    APIExit();
  } else {
    result = -1.0;
  }
  return(Py_BuildValue("f",result));
}

static PyObject *CmdBond(PyObject *dummy, PyObject *args)
{
  char *str1,*str2;
  int order,mode;
  OrthoLineType s1,s2;
  int ok=false;
  ok = PyArg_ParseTuple(args,"ssii",&str1,&str2,&order,&mode);
  if (ok) {
    APIEntry();
    SelectorGetTmp(str1,s1);
    SelectorGetTmp(str2,s2);
    ExecutiveBond(s1,s2,order,mode); /* TODO STATUS */
    SelectorFreeTmp(s1);
    SelectorFreeTmp(s2);
    APIExit();
  }
  return(APIStatus(ok));  
}

static PyObject *CmdDistance(PyObject *dummy, PyObject *args)
{
  char *str1,*str2;
  OrthoLineType s1,s2;
  float dist=-1.0;
  int ok=false;
  ok = PyArg_ParseTuple(args,"ss",&str1,&str2);
  if (ok) {
    APIEntry();
    SelectorGetTmp(str1,s1);
    SelectorGetTmp(str2,s2);
    dist = ExecutiveDistance(s1,s2);
    SelectorFreeTmp(s1);
    SelectorFreeTmp(s2);
    APIExit();
  }
  return(Py_BuildValue("f",dist));
}

static PyObject *CmdLabel(PyObject *self,   PyObject *args)
{
  char *str1,*str2;
  OrthoLineType s1;
  int quiet;
  int ok=false;
  ok = PyArg_ParseTuple(args,"ssi",&str1,&str2,&quiet);
  if (ok) {
    APIEntry();
    SelectorGetTmp(str1,s1);
    ExecutiveLabel(s1,str2,quiet); /* TODO STATUS */
    SelectorFreeTmp(s1);
    APIExit();
  }
  return(APIStatus(ok));

}

static PyObject *CmdAlter(PyObject *self,   PyObject *args)
{
  char *str1,*str2;
  int i1,quiet;
  OrthoLineType s1;
  int result=0;
  int ok=false;
  ok = PyArg_ParseTuple(args,"ssii",&str1,&str2,&i1,&quiet);
  if (ok) {
    APIEntry();
    SelectorGetTmp(str1,s1);
    result=ExecutiveIterate(s1,str2,i1,quiet); /* TODO STATUS */
    SelectorFreeTmp(s1);
    APIExit();
  }
  return Py_BuildValue("i",result);

}

static PyObject *CmdAlterState(PyObject *self,   PyObject *args)
{
  char *str1,*str2;
  int i1,i2,i3,quiet;
  OrthoLineType s1;

  int ok=false;
  ok = PyArg_ParseTuple(args,"issiii",&i1,&str1,&str2,&i2,&i3,&quiet);
  if (ok) {
    APIEntry();
    SelectorGetTmp(str1,s1);
    ExecutiveIterateState(i1,s1,str2,i2,i3,quiet); /* TODO STATUS */
    SelectorFreeTmp(s1);
    APIExit();
  }
  return(APIStatus(ok));

}

static PyObject *CmdCopy(PyObject *self,   PyObject *args)
{
  char *str1,*str2;
  int ok=false;
  ok = PyArg_ParseTuple(args,"ss",&str1,&str2);
  if (ok) {
    APIEntry();
    ExecutiveCopy(str1,str2); /* TODO STATUS */
    APIExit();
  }
  return(APIStatus(ok));
}

static PyObject *CmdRecolor(PyObject *self,   PyObject *args)
{
  char *str1;
  OrthoLineType s1;
  int ok=true;
  int rep=-1;
  ok = PyArg_ParseTuple(args,"si",&str1,&rep);
  PRINTFD(FB_CCmd)
    " CmdRecolor: called with %s.\n",str1
    ENDFD;

  if (ok) {
    APIEntry();
    if(WordMatch(str1,"all",true)<0)
      ExecutiveInvalidateRep(str1,rep,cRepInvColor);
    else {
      SelectorGetTmp(str1,s1);
      ExecutiveInvalidateRep(s1,rep,cRepInvColor);
      SelectorFreeTmp(s1); 
    }
    APIExit();
  } else {
    ok = -1; /* special error convention */
  }
  return(APIStatus(ok));
}

static PyObject *CmdRebuild(PyObject *self,   PyObject *args)
{
  char *str1;
  OrthoLineType s1;
  int ok=true;
  int rep=-1;
  ok = PyArg_ParseTuple(args,"si",&str1,&rep);
  PRINTFD(FB_CCmd)
    " CmdRebuild: called with %s.\n",str1
    ENDFD;

  if (ok) {
    APIEntry();
    if(WordMatch(str1,"all",true)<0)
      ExecutiveRebuildAll();
    else {
      SelectorGetTmp(str1,s1);
      ExecutiveInvalidateRep(s1,rep,cRepInvAll);
      SelectorFreeTmp(s1); 
    }
    APIExit();
  } else {
    ok = -1; /* special error convention */
  }
  return(APIStatus(ok));
}

static PyObject *CmdResetRate(PyObject *dummy, PyObject *args)
{
  APIEntry();
  ButModeResetRate();
  APIExit();
  return(APISuccess());
}

static PyObject *CmdReady(PyObject *dummy, PyObject *args)
{
  return(APIStatus(PyMOLReady));
}

#if 0
extern int _Py_CountReferences(void);
#endif
static PyObject *CmdMem(PyObject *dummy, PyObject *args)
{
  MemoryDebugDump();
  SelectorMemoryDump();
#if 0
  printf(" Py_Debug: %d total references.\n",_Py_CountReferences());
#endif
  return(APISuccess());
}

static PyObject *CmdRunPyMOL(PyObject *dummy, PyObject *args)
{
#ifndef _PYMOL_WX_GLUT

  if(run_only_once) {
    run_only_once=false;
#ifdef _PYMOL_MODULE
    was_main();
#endif
  }
#endif

  return(APISuccess());
}

static PyObject *CmdRunWXPyMOL(PyObject *dummy, PyObject *args)
{
#ifdef _PYMOL_WX_GLUT
#ifndef _PYMOL_ACTIVEX
#ifndef _EPYMOL
  if(run_only_once) {
    run_only_once=false;
    was_main();
  }
#endif
#endif
#endif

  return(APISuccess());
}

static PyObject *CmdCountStates(PyObject *dummy, PyObject *args)
{
  char *str1;
  OrthoLineType s1;
  int ok=false;
  ok = PyArg_ParseTuple(args,"s",&str1);
  if (ok) {
    APIEntry();
    SelectorGetTmp(str1,s1);
    ok = ExecutiveCountStates(s1);
    SelectorFreeTmp(s1); 
    APIExit();
  } else {
    ok = -1; /* special error convention */
  }
  return(APIStatus(ok));
}

static PyObject *CmdCountFrames(PyObject *dummy, PyObject *args)
{
  int result;
  APIEntry();
  SceneCountFrames();
  result=SceneGetNFrame();
  APIExit();
  return(APIStatus(result));
}

static PyObject *CmdIdentify(PyObject *dummy, PyObject *args)
{
  char *str1;
  OrthoLineType s1;
  int mode;
  int a,l=0;
  PyObject *result = Py_None;
  PyObject *tuple;
  int *iVLA=NULL,*i;
  ObjectMolecule **oVLA=NULL,**o;
  int ok=false;
  ok = PyArg_ParseTuple(args,"si",&str1,&mode);
  if (ok) {
    APIEntry();
    SelectorGetTmp(str1,s1);
    if(!mode) {
      iVLA=ExecutiveIdentify(s1,mode);
    } else {
      l = ExecutiveIdentifyObjects(s1,mode,&iVLA,&oVLA);
    }
    SelectorFreeTmp(s1);
    APIExit();
    if(iVLA) {
      if(!mode) {
        result=PConvIntVLAToPyList(iVLA);
      } else { /* object mode */
        result=PyList_New(l);
        i = iVLA;
        o = oVLA;
        for(a=0;a<l;a++) {
          tuple = PyTuple_New(2);      
          PyTuple_SetItem(tuple,1,PyInt_FromLong(*(i++))); 
          PyTuple_SetItem(tuple,0,PyString_FromString((*(o++))->Obj.Name));
          PyList_SetItem(result,a,tuple);
        }
      }
    } else {
      result = PyList_New(0);
    }
  }
  VLAFreeP(iVLA);
  VLAFreeP(oVLA);
  return(APIAutoNone(result));
}

static PyObject *CmdIndex(PyObject *dummy, PyObject *args)
{
  char *str1;
  OrthoLineType s1;
  int mode;
  PyObject *result = Py_None;
  PyObject *tuple = Py_None;
  int *iVLA=NULL;
  int l;
  int *i;
  ObjectMolecule **o,**oVLA=NULL;
  int a;

  int ok=false;
  ok = PyArg_ParseTuple(args,"si",&str1,&mode);
  if (ok) {
    APIEntry();
    SelectorGetTmp(str1,s1);
    l = ExecutiveIndex(s1,mode,&iVLA,&oVLA);
    SelectorFreeTmp(s1);
    APIExit();
    if(iVLA) {
      result=PyList_New(l);
      i = iVLA;
      o = oVLA;
      for(a=0;a<l;a++) {
        tuple = PyTuple_New(2);      
        PyTuple_SetItem(tuple,1,PyInt_FromLong(*(i++)+1)); /* +1 for index */
        PyTuple_SetItem(tuple,0,PyString_FromString((*(o++))->Obj.Name));
        PyList_SetItem(result,a,tuple);
      }
    } else {
      result = PyList_New(0);
    }
    VLAFreeP(iVLA);
    VLAFreeP(oVLA);
  }
  return(APIAutoNone(result));
}


static PyObject *CmdFindPairs(PyObject *dummy, PyObject *args)
{
  char *str1,*str2;
  int state1,state2;
  float cutoff;
  float angle;
  int mode;
  OrthoLineType s1,s2;
  PyObject *result = Py_None;
  PyObject *tuple = Py_None;
  PyObject *tuple1 = Py_None;
  PyObject *tuple2 = Py_None;
  int *iVLA=NULL;
  int l;
  int *i;
  ObjectMolecule **o,**oVLA=NULL;
  int a;
  
  int ok=false;
  ok = PyArg_ParseTuple(args,"ssiiiff",&str1,&str2,&state1,&state2,&mode,&cutoff,&angle);
  if (ok) {
    APIEntry();
    SelectorGetTmp(str1,s1);
    SelectorGetTmp(str2,s2);
    l = ExecutivePairIndices(s1,s2,state1,state2,mode,cutoff,angle,&iVLA,&oVLA);
    SelectorFreeTmp(s1);
    SelectorFreeTmp(s2);
    APIExit();
    
    if(iVLA&&oVLA) {
      result=PyList_New(l);
      i = iVLA;
      o = oVLA;
      for(a=0;a<l;a++) {
        tuple1 = PyTuple_New(2);      
        PyTuple_SetItem(tuple1,0,PyString_FromString((*(o++))->Obj.Name));
        PyTuple_SetItem(tuple1,1,PyInt_FromLong(*(i++)+1)); /* +1 for index */
        tuple2 = PyTuple_New(2);
        PyTuple_SetItem(tuple2,0,PyString_FromString((*(o++))->Obj.Name));
        PyTuple_SetItem(tuple2,1,PyInt_FromLong(*(i++)+1)); /* +1 for index */
        tuple = PyTuple_New(2);
        PyTuple_SetItem(tuple,0,tuple1);
        PyTuple_SetItem(tuple,1,tuple2);
        PyList_SetItem(result,a,tuple);
      }
    } else {
      result = PyList_New(0);
    }
    VLAFreeP(iVLA);
    VLAFreeP(oVLA);
  }
  return(APIAutoNone(result));
}

static PyObject *CmdSystem(PyObject *dummy, PyObject *args)
{
  char *str1;
  int ok=false;
  int async;
  ok = PyArg_ParseTuple(args,"si",&str1,&async);
  if (ok) {
    if(async) {
      PUnblock(); /* free up PyMOL and the API */
    } else {
      APIEntry(); /* keep PyMOL locked */
    }
    ok = system(str1);
    if(async) {
      PBlock();
    } else {
      APIExit();
    }
  }
  return(APIStatus(ok));
}

static PyObject *CmdGetFeedback(PyObject *dummy, PyObject *args)
{
  PyObject *result = NULL;
  OrthoLineType buffer;
  int ok;

  if(PyMOLTerminating) { /* try to bail */
#ifdef WIN32
	abort();
#endif
    exit(0);
  }
  APIEnterBlocked();
  ok = OrthoFeedbackOut(buffer); 
  APIExitBlocked();
  if(ok) result = Py_BuildValue("s",buffer);
  return(APIAutoNone(result));
}

static PyObject *CmdGetPDB(PyObject *dummy, PyObject *args)
{
  char *str1;
  char *pdb = NULL;
  int state;
  int mode;
  OrthoLineType s1 = "";
  PyObject *result = NULL;
  int ok=false;
  ok = PyArg_ParseTuple(args,"sii",&str1,&state,&mode);
  if (ok) {
    APIEntry();
    SelectorGetTmp(str1,s1);
    pdb=ExecutiveSeleToPDBStr(s1,state,true,mode);
    SelectorFreeTmp(s1);
    APIExit();
    if(pdb) result = Py_BuildValue("s",pdb);
    FreeP(pdb);
  }
  return(APIAutoNone(result));
}

static PyObject *CmdGetModel(PyObject *dummy, PyObject *args)
{
  char *str1;
  int state;
  OrthoLineType s1;
  PyObject *result = NULL;
  int ok=false;
  ok = PyArg_ParseTuple(args,"si",&str1,&state);
  if (ok) {
    APIEntry();
    SelectorGetTmp(str1,s1);
    result=ExecutiveSeleToChemPyModel(s1,state);
    SelectorFreeTmp(s1);
    APIExit();
  }
  return(APIAutoNone(result));
}

static PyObject *CmdCreate(PyObject *dummy, PyObject *args)
{
  char *str1,*str2;
  int target,source,discrete;
  OrthoLineType s1;
  int ok=false;
  ok = PyArg_ParseTuple(args,"ssiii",&str1,&str2,&source,&target,&discrete);
  if (ok) {
    APIEntry();
    SelectorGetTmp(str2,s1);
    ExecutiveSeleToObject(str1,s1,source,target,discrete); /* TODO STATUS */
    SelectorFreeTmp(s1);
    APIExit();
  }
  return(APIStatus(ok));
}


static PyObject *CmdOrient(PyObject *dummy, PyObject *args)
{
  Matrix33d m;
  char *str1;
  OrthoLineType s1;
  int state;
  int ok=false;
  ok = PyArg_ParseTuple(args,"si",&str1,&state);
  if (ok) {
    APIEntry();
    SelectorGetTmp(str1,s1);
    if(ExecutiveGetMoment(s1,m,state))
      ExecutiveOrient(s1,m,state); /* TODO STATUS */
    else
      ok=false;
    SelectorFreeTmp(s1); 
    APIExit();
  }
  return(APIStatus(ok));
}

static PyObject *CmdFitPairs(PyObject *dummy, PyObject *args)
{
  PyObject *list;
  WordType *word = NULL;
  int ln=0;
  int a;
  PyObject *result = NULL;
  float valu;
  int ok=false;
  ok = PyArg_ParseTuple(args,"O",&list);
  if(ok) {
    ln = PyObject_Length(list);
    if(ln) {
      if(ln&0x1)
        ok=ErrMessage("FitPairs","must supply an even number of selections.");
    } else ok=false;
    
    if(ok) {
      word = Alloc(WordType,ln);
      
      a=0;
      while(a<ln) {
        SelectorGetTmp(PyString_AsString(PySequence_GetItem(list,a)),word[a]);
        a++;
      }
      APIEntry();
      valu = ExecutiveRMSPairs(word,ln/2,2);
      APIExit();
      result=Py_BuildValue("f",valu);
      for(a=0;a<ln;a++)
        SelectorFreeTmp(word[a]);
      FreeP(word);
    }
  }
  return APIAutoNone(result);
}

static PyObject *CmdIntraFit(PyObject *dummy, PyObject *args)
{
  char *str1;
  int state;
  int mode;
  int quiet;
  OrthoLineType s1;
  float *fVLA;
  PyObject *result=Py_None;
  int ok=false;
  ok = PyArg_ParseTuple(args,"siii",&str1,&state,&mode,&quiet);
  if(state<0) state=0;
  if (ok) {
    APIEntry();
    SelectorGetTmp(str1,s1);
    fVLA=ExecutiveRMSStates(s1,state,mode,quiet);
    SelectorFreeTmp(s1);
    APIExit();
    if(fVLA) {
      result=PConvFloatVLAToPyList(fVLA);
      VLAFreeP(fVLA);
    }
  }
  return APIAutoNone(result);
}

static PyObject *CmdGetAtomCoords(PyObject *dummy, PyObject *args)
{
  char *str1;
  int state;
  int quiet;
  OrthoLineType s1;
  float vertex[3];
  PyObject *result=Py_None;
  int ok=false;
  ok = PyArg_ParseTuple(args,"sii",&str1,&state,&quiet);
  if (ok) {
    APIEntry();
    SelectorGetTmp(str1,s1);
    ok=ExecutiveGetAtomVertex(s1,state,quiet,vertex);
    SelectorFreeTmp(s1);
    APIExit();
    if(ok) {
      result=PConvFloatArrayToPyList(vertex,3);
    }
  }
  return APIAutoNone(result);
}

static PyObject *CmdFit(PyObject *dummy, PyObject *args)
{
  char *str1,*str2;
  int mode;
  int quiet;
  OrthoLineType s1,s2;
  PyObject *result;
  float tmp_result = -1.0;
  int state1,state2;
  int ok=false;
  ok = PyArg_ParseTuple(args,"ssiiii",&str1,&str2,&mode,&state1,&state2,&quiet);
  if (ok) {
    APIEntry();
    SelectorGetTmp(str1,s1);
    SelectorGetTmp(str2,s2);
    tmp_result=ExecutiveRMS(s1,s2,mode,0.0,0,quiet,NULL,state1,state2,false);
    SelectorFreeTmp(s1);
    SelectorFreeTmp(s2);
    APIExit();
  }
  result=Py_BuildValue("f",tmp_result);
  return result;
}

static PyObject *CmdUpdate(PyObject *dummy, PyObject *args)
{
  char *str1,*str2;
  int int1,int2;
  OrthoLineType s1,s2;
  int ok=false;
  ok = PyArg_ParseTuple(args,"ssii",&str1,&str2,&int1,&int2);
  if (ok) {
    APIEntry();
    SelectorGetTmp(str1,s1);
    SelectorGetTmp(str2,s2);
    ExecutiveUpdateCmd(s1,s2,int1,int2); /* TODO STATUS */
    SelectorFreeTmp(s1);
    SelectorFreeTmp(s2);
    APIExit();
  }
  return(APIStatus(ok));
}

static PyObject *CmdDirty(PyObject *self, 	PyObject *args)
{
  PRINTFD(FB_CCmd)
    " CmdDirty: called.\n"
    ENDFD;
  APIEntry();
  OrthoDirty();
  APIExit();
  return(APISuccess());
}

static PyObject *CmdGetDihe(PyObject *self, 	PyObject *args)
{
  char *str1,*str2,*str3,*str4;
  float result;
  int int1;
  OrthoLineType s1,s2,s3,s4;
  int ok=false;
  ok = PyArg_ParseTuple(args,"ssssi",&str1,&str2,&str3,&str4,&int1);
  
  if (ok) {
    APIEntry();
    SelectorGetTmp(str1,s1);
    SelectorGetTmp(str2,s2);
    SelectorGetTmp(str3,s3);
    SelectorGetTmp(str4,s4);
    ok = ExecutiveGetDihe(s1,s2,s3,s4,&result,int1);
    SelectorFreeTmp(s1);
    SelectorFreeTmp(s2);
    SelectorFreeTmp(s3);
    SelectorFreeTmp(s4);
    APIExit();
  }
  
  if(ok) {
    return(Py_BuildValue("f",result));
  } else {
    return APIFailure();
  }
}

static PyObject *CmdSetDihe(PyObject *self, 	PyObject *args)
{
  char *str1,*str2,*str3,*str4;
  float float1;
  int int1;
  int quiet;
  OrthoLineType s1,s2,s3,s4;
  int ok=false;
  ok = PyArg_ParseTuple(args,"ssssfii",&str1,&str2,&str3,&str4,&float1,&int1,&quiet);
  if (ok) {
    APIEntry();
    SelectorGetTmp(str1,s1);
    SelectorGetTmp(str2,s2);
    SelectorGetTmp(str3,s3);
    SelectorGetTmp(str4,s4);
    ok = ExecutiveSetDihe(s1,s2,s3,s4,float1,int1,quiet);
    SelectorFreeTmp(s1);
    SelectorFreeTmp(s2);
    SelectorFreeTmp(s3);
    SelectorFreeTmp(s4);
    APIExit();
  }
  return(APIStatus(ok));
}

static PyObject *CmdDo(PyObject *self, 	PyObject *args)
{
  char *str1;
  int log;
  int ok=false;
  ok = PyArg_ParseTuple(args,"si",&str1,&log);
  if (ok) {
    APIEntry();
    if(str1[0]!='_') { /* suppress internal call-backs */
      if(strncmp(str1,"cmd._",5)&&(strncmp(str1,"_cmd.",5))) {
        OrthoAddOutput("PyMOL>");
        OrthoAddOutput(str1);
        OrthoNewLine(NULL,true);
        if(log) 
          if(WordMatch(str1,"quit",true)==0) /* don't log quit */
            PLog(str1,cPLog_pml);
      }
      PParse(str1);
    } else if(str1[1]==' ') { /* "_ command" suppresses echoing of command, but it is still logged */
      if(log)
        if(WordMatch(str1+2,"quit",true)==0) /* don't log quit */
          PLog(str1+2,cPLog_pml);
      PParse(str1+2);    
    } else {
      PParse(str1);
    }
    APIExit();
  }
  return(APIStatus(ok));
}

static PyObject *CmdRock(PyObject *self, PyObject *args)
{
  int int1;
  int ok=true;
  ok = PyArg_ParseTuple(args,"i",&int1);
  APIEntry();
  ControlRock(int1);
  APIExit();
  return(APIStatus(ok));
}

static PyObject *CmdGetMoment(PyObject *self, 	PyObject *args) /* missing? */
{
  Matrix33d m;
  PyObject *result;
  char *str1;
  int ok=false;
  int state;

  ok = PyArg_ParseTuple(args,"si",&str1,&state);
  if (ok) {
    APIEntry();
    ExecutiveGetMoment(str1,m,state);
    APIExit();
  }
  result = Py_BuildValue("(ddd)(ddd)(ddd)", 
								 m[0][0],m[0][1],m[0][2],
								 m[1][0],m[1][1],m[1][2],
								 m[2][0],m[2][1],m[2][2]);
  return result;
}

static PyObject *CmdGetSetting(PyObject *self, 	PyObject *args)
{
  PyObject *result = Py_None;
  char *str1;
  float value;
  int ok=false;
  ok = PyArg_ParseTuple(args,"s",&str1);
  if (ok) {
    APIEnterBlocked();
    value=SettingGetNamed(str1);
    APIExitBlocked();
    result = Py_BuildValue("f", value);
  }
  return APIAutoNone(result);
}

static PyObject *CmdGetSettingTuple(PyObject *self, 	PyObject *args)
{
  PyObject *result = Py_None;
  int int1,int2;
  char *str1;
  int ok = false;
  ok = PyArg_ParseTuple(args,"isi",&int1,&str1,&int2); /* setting, object, state */
  if (ok) {
    APIEnterBlocked();
    result =  ExecutiveGetSettingTuple(int1,str1,int2);
    APIExitBlocked();
  }
  return APIAutoNone(result);
}

static PyObject *CmdGetSettingText(PyObject *self, 	PyObject *args)
{
  PyObject *result = Py_None;
  int int1,int2;
  char *str1;
  int ok=false;
  ok = PyArg_ParseTuple(args,"isi",&int1,&str1,&int2); /* setting, object, state */
  if (ok) {
    APIEnterBlocked();
    result =  ExecutiveGetSettingText(int1,str1,int2);
    APIExitBlocked();
  }
  return APIAutoNone(result);
}

static PyObject *CmdExportDots(PyObject *self, 	PyObject *args)
{
  PyObject *result=NULL;
  PyObject *cObj;
  ExportDotsObj *obj;
  char *str1;
  int int1;
  
  int ok=false;
  ok = PyArg_ParseTuple(args,"si",&str1,&int1);
  if (ok) {
    APIEntry();
    obj = ExportDots(str1,int1);
    APIExit();
    if(obj) 
      {
        cObj = PyCObject_FromVoidPtr(obj,(void(*)(void*))ExportDeleteMDebug);
        if(cObj) {
          result = Py_BuildValue("O",cObj); /* I think this */
          Py_DECREF(cObj); /* transformation is unnecc. */
        } 
      }
  }
  return APIAutoNone(result);
}

static PyObject *CmdSetFrame(PyObject *self, PyObject *args)
{
  int mode,frm;
  int ok=false;
  ok = PyArg_ParseTuple(args,"ii",&mode,&frm);
  if (ok) {
    APIEntry();
    SceneSetFrame(mode,frm);
    APIExit();
  }
  return(APIStatus(ok));
}

static PyObject *CmdFrame(PyObject *self, PyObject *args)
{
  int frm;
  int ok=false;
  ok = PyArg_ParseTuple(args,"i",&frm);
  frm--;
  if(frm<0) frm=0;
  if (ok) {
    APIEntry();
    SceneSetFrame(4,frm);
    APIExit();
  }
  return(APIStatus(ok));
}

static PyObject *CmdStereo(PyObject *self, PyObject *args)
{
  int i1;
  int ok=false;
  
  ok = PyArg_ParseTuple(args,"i",&i1);
  if (ok) {
    APIEntry();
    ok = ExecutiveStereo(i1); 
    APIExit();
  }
  return APIStatus(ok);
}

static PyObject *CmdReset(PyObject *self, PyObject *args)
{
  int cmd;
  int ok=false;
  char *obj;
  ok = PyArg_ParseTuple(args,"is",&cmd,&obj);
  if (ok) {
    APIEntry();
    ok = ExecutiveReset(cmd,obj); 
    APIExit();
  }
  return(APIStatus(ok));
}

static PyObject *CmdSetMatrix(PyObject *self, 	PyObject *args)
{
  float m[16];
  int ok=false;
  ok = PyArg_ParseTuple(args,"ffffffffffffffff",
						 &m[0],&m[1],&m[2],&m[3],
						 &m[4],&m[5],&m[6],&m[7],
						 &m[8],&m[9],&m[10],&m[11],
						 &m[12],&m[13],&m[14],&m[15]);
  if (ok) {
    APIEntry();
    SceneSetMatrix(m); /* TODO STATUS */
    APIExit();
  }
  return(APIStatus(ok));
}

static PyObject *CmdGetMinMax(PyObject *self, 	PyObject *args)
{
  float mn[3],mx[3];
  char *str1;
  int state;
  OrthoLineType s1;
  PyObject *result = Py_None;
  int flag;

  int ok=false;
  ok = PyArg_ParseTuple(args,"si",&str1,&state); 
  if (ok) {
    APIEntry();
    SelectorGetTmp(str1,s1);
    flag = ExecutiveGetExtent(s1,mn,mx,true,state,false);
    SelectorFreeTmp(s1);
    if(flag) 
      result = Py_BuildValue("[[fff],[fff]]", 
                             mn[0],mn[1],mn[2],
                             mx[0],mx[1],mx[2]);
    else 
      result = Py_BuildValue("[[fff],[fff]]", 
                             -0.5,-0.5,-0.5,
                             0.5,0.5,0.5);
    
    APIExit();
  }
  return APIAutoNone(result);
}

static PyObject *CmdGetMatrix(PyObject *self, 	PyObject *args)
{
  float *f;
  PyObject *result;
  
  APIEntry();
  f=SceneGetMatrix();
  APIExit();
  result = Py_BuildValue("ffffffffffffffff", 
								 f[0],f[1],f[2],f[3],
								 f[4],f[5],f[6],f[7],
								 f[8],f[9],f[10],f[11],
								 f[12],f[13],f[14],f[15]
								 );
  return result;
}

static PyObject *CmdMDo(PyObject *self, 	PyObject *args)
{
  char *cmd;
  int frame;
  int append;
  int ok=false;
  ok = PyArg_ParseTuple(args,"isi",&frame,&cmd,&append);
  if (ok) {
    APIEntry();
    if(append) {
      MovieAppendCommand(frame,cmd);
    } else {
      MovieSetCommand(frame,cmd);
    }
    APIExit();
  }
  return(APIStatus(ok));
}

static PyObject *CmdMPlay(PyObject *self, 	PyObject *args)
{
  int cmd;
  int ok=false;
  ok = PyArg_ParseTuple(args,"i",&cmd);
  if (ok) {
    APIEntry();
    MoviePlay(cmd);
    APIExit();
  }
  return(APIStatus(ok));
}

static PyObject *CmdMMatrix(PyObject *self, 	PyObject *args)
{
  int cmd;
  int ok=false;
  ok = PyArg_ParseTuple(args,"i",&cmd);
  if (ok) {
    APIEntry();
    ok = MovieMatrix(cmd);
    APIExit();
  }
  return(APIStatus(ok));
}

static PyObject *CmdMClear(PyObject *self, 	PyObject *args)
{
  APIEntry();
  MovieClearImages();
  APIExit();
  return(APISuccess());
}

static PyObject *CmdRefresh(PyObject *self, 	PyObject *args)
{
  APIEntry();
  ExecutiveDrawNow(); /* TODO STATUS */
  APIExit();
  return(APISuccess());
}

static PyObject *CmdRefreshNow(PyObject *self, 	PyObject *args)
{
  APIEntry();
  ExecutiveDrawNow(); /* TODO STATUS */
  MainRefreshNow();
  APIExit();
  return(APISuccess());
}

static PyObject *CmdPNG(PyObject *self, 	PyObject *args)
{
  char *str1;
  int ok=false;
  int quiet;
  ok = PyArg_ParseTuple(args,"si",&str1,&quiet);
  if (ok) {
    APIEntry();
    ExecutiveDrawNow();		 /* TODO STATUS */
    ScenePNG(str1,quiet);
    APIExit();
  }
  return(APIStatus(ok));
}

static PyObject *CmdMPNG(PyObject *self, 	PyObject *args)
{
  char *str1;
  int int1,int2;
  int ok=false;
  ok = PyArg_ParseTuple(args,"sii",&str1,&int1,&int2);
  if (ok) {
    APIEntry();
    ok = MoviePNG(str1,(int)SettingGet(cSetting_cache_frames),int1,int2);
    /* TODO STATUS */
    APIExit();
  }
  return(APIStatus(ok));
}

static PyObject *CmdMSet(PyObject *self, 	PyObject *args)
{
  char *str1;
  int ok=false;
  ok = PyArg_ParseTuple(args,"s",&str1);
  if (ok) {
    APIEntry();
    MovieSequence(str1);
    SceneCountFrames();
    APIExit();
  }
  return(APIStatus(ok));
}

static PyObject *CmdMView(PyObject *self, 	PyObject *args)
{
  int ok=false;
  int action,first,last;
  ok = PyArg_ParseTuple(args,"iii",&action,&first,&last);
  if (ok) {
    APIEntry();
    ok = MovieView(action,first,last);
    APIExit();
  }
  return(APIStatus(ok));
}

static PyObject *CmdViewport(PyObject *self, 	PyObject *args)
{
  int w,h;
  int ok=false;
  ok = PyArg_ParseTuple(args,"ii",&w,&h);
  if(ok) {
    if((w>0)&&(h>0)) {
      if(w<10) w=10;
      if(h<10) h=10;
      if(SettingGet(cSetting_internal_gui)) {
        if(!SettingGet(cSetting_full_screen))
           w+=(int)SettingGet(cSetting_internal_gui_width);
      }
      if(SettingGet(cSetting_internal_feedback)) {
        if(!SettingGet(cSetting_full_screen))
          h+=(int)(SettingGet(cSetting_internal_feedback)-1)*cOrthoLineHeight +
            cOrthoBottomSceneMargin;
      }
    } else {
      w=-1;
      h=-1;
    }
    APIEntry();
    MainDoReshape(w,h); /* should be moved into Executive */
    APIExit();
  }
  return(APIStatus(ok));
}

static PyObject *CmdFlag(PyObject *self, 	PyObject *args)
{
  char *str1;
  int flag;
  int action;
  int quiet;
  OrthoLineType s1;
  int ok=false;
  ok = PyArg_ParseTuple(args,"isii",&flag,&str1,&action,&quiet);
  if (ok) {
    APIEntry();
    SelectorGetTmp(str1,s1);
    ExecutiveFlag(flag,s1,action,quiet);
    SelectorFreeTmp(s1);
    APIExit();
  }
  return(APIStatus(ok));
}

static PyObject *CmdColor(PyObject *self, 	PyObject *args)
{
  char *str1,*color;
  int flags;
  OrthoLineType s1;
  int ok=false;
  int quiet;
  ok = PyArg_ParseTuple(args,"ssii",&color,&str1,&flags,&quiet);
  if (ok) {
    APIEntry();
    SelectorGetTmp(str1,s1);
    ok = ExecutiveColor(s1,color,flags,quiet);
    SelectorFreeTmp(s1);
    APIExit();
  }
  return(APIStatus(ok));
}

static PyObject *CmdColorDef(PyObject *self, 	PyObject *args)
{
  char *color;
  float v[3];
  int ok=false;
  ok = PyArg_ParseTuple(args,"sfff",&color,v,v+1,v+2);
  if (ok) {
    APIEntry();
    ColorDef(color,v);
    APIExit();
  }
  return(APIStatus(ok));
}

static PyObject *CmdRay(PyObject *self, 	PyObject *args)
{
  int w,h,mode;
  float angle,shift;
  int ok=false;
  int quiet;
  ok = PyArg_ParseTuple(args,"iiiffi",&w,&h,&mode,&angle,&shift,&quiet);
  if (ok) {
    APIEntry();
    if(mode<0)
      mode=(int)SettingGet(cSetting_ray_default_renderer);
    ExecutiveRay(w,h,mode,angle,shift,quiet); /* TODO STATUS */
    APIExit();
  }
  return(APIStatus(ok));
}

static PyObject *CmdClip(PyObject *self, 	PyObject *args)
{
  char *sname;
  float dist;
  char *str1;
  int state;
  OrthoLineType s1;
  int ok=false;
  ok = PyArg_ParseTuple(args,"sfsi",&sname,&dist,&str1,&state);
  if (ok) {
    APIEntry();
    SelectorGetTmp(str1,s1);
    switch(sname[0]) { /* TODO STATUS */
    case 'N':
    case 'n':
      SceneClip(0,dist,s1,state);
      break;
    case 'f':
    case 'F':
      SceneClip(1,dist,s1,state);
      break;
    case 'm':
    case 'M':
      SceneClip(2,dist,s1,state);
      break;
    case 's':
    case 'S':
      SceneClip(3,dist,s1,state);
      break;
    case 'a':
    case 'A':
      SceneClip(4,dist,s1,state);
      break;
    }
    SelectorFreeTmp(s1);
    APIExit();
  }
  return(APIStatus(ok));

}

static PyObject *CmdMove(PyObject *self, 	PyObject *args)
{
  char *sname;
  float dist;
  int ok=false;
  ok = PyArg_ParseTuple(args,"sf",&sname,&dist);
  if (ok) {
    APIEntry();
    switch(sname[0]) { /* TODO STATUS */
    case 'x':
      SceneTranslate(dist,0.0,0.0);
      break;
    case 'y':
      SceneTranslate(0.0,dist,0.0);
      break;
    case 'z':
      SceneTranslate(0.0,0.0,dist);
      break;
    }
    APIExit();
  }
  return(APIStatus(ok));
}

static PyObject *CmdTurn(PyObject *self, 	PyObject *args)
{
  char *sname;
  float angle;
  int ok=false;
  ok = PyArg_ParseTuple(args,"sf",&sname,&angle);
  if (ok) {
    APIEntry();
    switch(sname[0]) { /* TODO STATUS */
    case 'x':
      SceneRotate(angle,1.0,0.0,0.0);
      break;
    case 'y':
      SceneRotate(angle,0.0,1.0,0.0);
      break;
    case 'z':
      SceneRotate(angle,0.0,0.0,1.0);
      break;
    }
    APIExit();
  }
  return(APIStatus(ok));
}

static PyObject *CmdLegacySet(PyObject *self, 	PyObject *args)
{
  char *sname, *value;
  int ok=false;
  ok = PyArg_ParseTuple(args,"ss",&sname,&value);
  if (ok) {
    APIEntry();
    ok = SettingSetNamed(sname,value); 
    APIExit();
  }
  return(APIStatus(ok));
}

static PyObject *CmdUnset(PyObject *self, 	PyObject *args)
{
  int index;
  int tmpFlag=false;
  char *str3;
  int state;
  int quiet;
  int updates;
  OrthoLineType s1;
  int ok=false;
  ok = PyArg_ParseTuple(args,"isiii",&index,&str3,&state,&quiet,&updates);
  s1[0]=0;
  if (ok) {
    APIEntry();
    if(!strcmp(str3,"all")) {
      strcpy(s1,str3);
    } else if(str3[0]!=0) {
      tmpFlag=true;
      SelectorGetTmp(str3,s1);
    }
    ok = ExecutiveUnsetSetting(index,s1,state,quiet,updates);
    if(tmpFlag) 
      SelectorFreeTmp(s1);
    APIExit();
  }
  return(APIStatus(ok));
}

static PyObject *CmdSet(PyObject *self, 	PyObject *args)
{
  int index;
  int tmpFlag=false;
  PyObject *value;
  char *str3;
  int state;
  int quiet;
  int updates;
  OrthoLineType s1;
  int ok=false;
  ok = PyArg_ParseTuple(args,"iOsiii",&index,&value,&str3,&state,&quiet,&updates);
  s1[0]=0;
  if (ok) {
    APIEntry();
    if(!strcmp(str3,"all")) {
      strcpy(s1,str3);
    } else if(str3[0]!=0) {
      tmpFlag=true;
      SelectorGetTmp(str3,s1);
    }
    ok = ExecutiveSetSetting(index,value,s1,state,quiet,updates);
    if(tmpFlag) 
      SelectorFreeTmp(s1);
    APIExit();
  }
  return(APIStatus(ok));
}

static PyObject *CmdGet(PyObject *self, 	PyObject *args)
{
  float f;
  char *sname;
  PyObject *result = Py_None;

  int ok=false;
  ok = PyArg_ParseTuple(args,"s",&sname);
  if (ok) {
    APIEnterBlocked();
    f=SettingGetNamed(sname);
    APIExitBlocked();
    result = Py_BuildValue("f", f);
  }
  return APIAutoNone(result);
}

static PyObject *CmdDelete(PyObject *self, 	PyObject *args)
{
  char *sname;

  int ok=false;
  ok = PyArg_ParseTuple(args,"s",&sname);
  if (ok) {
    APIEntry();
    ExecutiveDelete(sname); /* TODO STATUS */
    APIExit();
  }
  return(APIStatus(ok));
}

static PyObject *CmdCartoon(PyObject *self, 	PyObject *args)
{
  char *sname;
  int type;
  OrthoLineType s1;
  int ok=false;
  ok = PyArg_ParseTuple(args,"si",&sname,&type);
  if (ok) {
    APIEntry();
    SelectorGetTmp(sname,s1);
    ExecutiveCartoon(type,s1); /* TODO STATUS */
    SelectorFreeTmp(s1);
    APIExit();
  }
  return(APIStatus(ok));
}

static PyObject *CmdShowHide(PyObject *self, 	PyObject *args)
{
  char *sname;
  int rep;
  int state;
  OrthoLineType s1;
  int ok=false;
  ok = PyArg_ParseTuple(args,"sii",&sname,&rep,&state);
  if (ok) { /* TODO STATUS */
    APIEntry();
    if(sname[0]=='@') {
      ExecutiveSetAllVisib(state);
    } else {
      SelectorGetTmp(sname,s1);
      ExecutiveSetRepVisib(s1,rep,state);
      SelectorFreeTmp(s1);
    }
    APIExit();
  }
  return(APIStatus(ok));
}

static PyObject *CmdOnOffBySele(PyObject *self, 	PyObject *args)
{
  char *sname;
  int onoff;
  OrthoLineType s1;
  int ok=false;
  ok = PyArg_ParseTuple(args,"si",&sname,&onoff);
  if (ok) { /* TODO STATUS */
    APIEntry();
    SelectorGetTmp(sname,s1);
    ok = ExecutiveSetOnOffBySele(s1,onoff);
    SelectorFreeTmp(s1);
    APIExit();
  }
  return(APIStatus(ok));
}

static PyObject *CmdOnOff(PyObject *self, 	PyObject *args)
{
  char *name;
  int state;
  int ok=false;
  ok = PyArg_ParseTuple(args,"si",&name,&state);
  if (ok) { /* TODO STATUS */
    APIEntry();
    ExecutiveSetObjVisib(name,state);
    APIExit();
  }
  return(APIStatus(ok));
}

static PyObject *CmdToggle(PyObject *self, 	PyObject *args)
{
  char *sname;
  int rep;
  OrthoLineType s1;
  int ok=false;
  ok = PyArg_ParseTuple(args,"si",&sname,&rep);
  if (ok) {
    APIEntry();
    if(sname[0]=='@') {
      /* TODO */
    } else {
      SelectorGetTmp(sname,s1);
      ok = ExecutiveToggleRepVisib(s1,rep);
      SelectorFreeTmp(s1);
    }
    APIExit();
  }
  return(APIStatus(ok));
}

static PyObject *CmdQuit(PyObject *self, 	PyObject *args)
{
  APIEntry();
  PyMOLTerminating=true;
  PExit(EXIT_SUCCESS);
  APIExit();
  return(APISuccess());
}
static PyObject *CmdFullScreen(PyObject *self,PyObject *args)
{
  int flag = 0;
  int ok=false;
  ok = PyArg_ParseTuple(args,"i",&flag);
  if (ok) {
    APIEntry();
    ExecutiveFullScreen(flag); /* TODO STATUS */
    APIExit();
  }
  return(APIStatus(ok));
}
static PyObject *CmdSelect(PyObject *self, PyObject *args)
{
  char *sname,*sele;
  int quiet;
  int ok=false;
  ok = PyArg_ParseTuple(args,"ssi",&sname,&sele,&quiet);
  if (ok) {
    APIEntry();
    ok = SelectorCreate(sname,sele,NULL,quiet,NULL);
    SceneDirty();
    SeqDirty();
    APIExit();
  } else {
    ok=-1;
  }
  return APIStatus(ok);
}

static PyObject *CmdFinishObject(PyObject *self, PyObject *args)
{
  char *oname;
  CObject *origObj = NULL;

  int ok=false;
  ok = PyArg_ParseTuple(args,"s",&oname);

  if (ok) {
    APIEntry();
    origObj=ExecutiveFindObjectByName(oname);
    if(origObj) 
      ExecutiveUpdateObjectSelection(origObj); /* TODO STATUS */
    else
      ok=false;
    APIExit();
  }
  return(APIStatus(ok));
}

static PyObject *CmdLoadObject(PyObject *self, PyObject *args)
{
  char *oname;
  PyObject *model;
  CObject *origObj = NULL,*obj;
  OrthoLineType buf;
  int frame,type;
  int finish,discrete;
  int quiet;
  int ok=false;
  ok = PyArg_ParseTuple(args,"sOiiiii",&oname,&model,&frame,&type,&finish,&discrete,&quiet);
  buf[0]=0;
  if (ok) {
    APIEntry();
    origObj=ExecutiveFindObjectByName(oname);
    
    /* TODO check for existing object of wrong type */
    
    switch(type) {
    case cLoadTypeChemPyModel:
      if(origObj)
        if(origObj->type!=cObjectMolecule) {
          ExecutiveDelete(oname);
          origObj=NULL;
        }
      PBlock(); /*PBlockAndUnlockAPI();*/
      obj=(CObject*)ObjectMoleculeLoadChemPyModel((ObjectMolecule*)
                                                 origObj,model,frame,discrete);
      PUnblock(); /*PLockAPIAndUnblock();*/
      if(!origObj) {
        if(obj) {
          ObjectSetName(obj,oname);
          ExecutiveManageObject(obj,true,quiet);
          if(frame<0)
            frame = ((ObjectMolecule*)obj)->NCSet-1;
          sprintf(buf," CmdLoad: ChemPy-model loaded into object \"%s\", state %d.\n",
                  oname,frame+1);		  
        }
      } else if(origObj) {
        if(finish)
          ExecutiveUpdateObjectSelection(origObj);
        if(frame<0)
          frame = ((ObjectMolecule*)origObj)->NCSet-1;
        sprintf(buf," CmdLoad: ChemPy-model appended into object \"%s\", state %d.\n",
                oname,frame+1);
      }
      break;
    case cLoadTypeChemPyBrick:
      if(origObj)
        if(origObj->type!=cObjectMap) {
          ExecutiveDelete(oname);
          origObj=NULL;
        }
      PBlock(); /*PBlockAndUnlockAPI();*/
      obj=(CObject*)ObjectMapLoadChemPyBrick((ObjectMap*)origObj,model,frame,discrete);
      PUnblock(); /*PLockAPIAndUnblock();*/
      if(!origObj) {
        if(obj) {
          ObjectSetName(obj,oname);
          ExecutiveManageObject(obj,true,quiet);
          sprintf(buf," CmdLoad: chempy.brick loaded into object \"%s\"\n",
                  oname);		  
        }
      } else if(origObj) {
        sprintf(buf," CmdLoad: chempy.brick appended into object \"%s\"\n",
                oname);
      }
      break;
    case cLoadTypeChemPyMap:
      if(origObj)
        if(origObj->type!=cObjectMap) {
          ExecutiveDelete(oname);
          origObj=NULL;
        }
      PBlock(); /*PBlockAndUnlockAPI();*/
      obj=(CObject*)ObjectMapLoadChemPyMap((ObjectMap*)origObj,model,frame,discrete);
      PUnblock(); /*PLockAPIAndUnblock();*/
      if(!origObj) {
        if(obj) {
          ObjectSetName(obj,oname);
          ExecutiveManageObject(obj,true,quiet);
          sprintf(buf," CmdLoad: chempy.map loaded into object \"%s\"\n",
                  oname);		  
        }
      } else if(origObj) {
        sprintf(buf," CmdLoad: chempy.map appended into object \"%s\"\n",
                oname);
      }
      break;
    case cLoadTypeCallback:
      if(origObj)
        if(origObj->type!=cObjectCallback) {
          ExecutiveDelete(oname);
          origObj=NULL;
        }
      PBlock(); /*PBlockAndUnlockAPI();*/
      obj=(CObject*)ObjectCallbackDefine((ObjectCallback*)origObj,model,frame);
      PUnblock(); /*PLockAPIAndUnblock();*/
      if(!origObj) {
        if(obj) {
          ObjectSetName(obj,oname);
          ExecutiveManageObject(obj,true,quiet);
          sprintf(buf," CmdLoad: pymol.callback loaded into object \"%s\"\n",
                  oname);		  
        }
      } else if(origObj) {
        sprintf(buf," CmdLoad: pymol.callback appended into object \"%s\"\n",
                oname);
      }
      break;
    case cLoadTypeCGO:
      if(origObj)
        if(origObj->type!=cObjectCGO) {
          ExecutiveDelete(oname);
          origObj=NULL;
        }
      PBlock(); /*PBlockAndUnlockAPI();*/
      obj=(CObject*)ObjectCGODefine((ObjectCGO*)origObj,model,frame);
      PUnblock(); /*PLockAPIAndUnblock();*/
      if(!origObj) {
        if(obj) {
          ObjectSetName(obj,oname);
          ExecutiveManageObject(obj,true,quiet);
          sprintf(buf," CmdLoad: CGO loaded into object \"%s\"\n",
                  oname);		  
        }
      } else if(origObj) {
        sprintf(buf," CmdLoad: CGO appended into object \"%s\"\n",
                oname);
      }
      break;
      
    }
    if(origObj&&!quiet) {
      PRINTFB(FB_Executive,FB_Actions) 
        "%s",buf
        ENDFB;
      OrthoRestorePrompt();
    }
    APIExit();
  }
  return(APIStatus(ok));
}

static PyObject *CmdLoadCoords(PyObject *self, PyObject *args)
{
  char *oname;
  PyObject *model;
  CObject *origObj = NULL,*obj;
  OrthoLineType buf;
  int frame,type;
  int ok=false;

  buf[0]=0;

  ok = PyArg_ParseTuple(args,"sOii",&oname,&model,&frame,&type);

  if (ok) {
    APIEntry();
    origObj=ExecutiveFindObjectByName(oname);
    
    /* TODO check for existing object of wrong type */
    if(!origObj) {
      ErrMessage("LoadCoords","named object not found.");
      ok=false;
    } else 
      {
        switch(type) {
        case cLoadTypeChemPyModel:
          PBlock(); /*PBlockAndUnlockAPI();*/
          obj=(CObject*)ObjectMoleculeLoadCoords((ObjectMolecule*)origObj,model,frame);
          PUnblock(); /*PLockAPIAndUnblock();*/
          if(frame<0)
            frame=((ObjectMolecule*)obj)->NCSet-1;
          sprintf(buf," CmdLoad: Coordinates appended into object \"%s\", state %d.\n",
                  oname,frame+1);
          break;
        }
      }
    if(origObj) {
      PRINTFB(FB_Executive,FB_Actions) 
        "%s",buf
        ENDFB;
      OrthoRestorePrompt();
    }
    APIExit();
  }
  return(APIStatus(ok));
}

static PyObject *CmdLoad(PyObject *self, PyObject *args)
{
  char *fname,*oname;
  CObject *origObj = NULL,*obj;
  OrthoLineType buf;
  int frame,type;
  int finish,discrete;
  int new_type;
  int quiet;
  int ok=false;
  ok = PyArg_ParseTuple(args,"ssiiiii",&oname,&fname,&frame,&type,&finish,&discrete,&quiet);

  buf[0]=0;
  PRINTFD(FB_CCmd)
    "CmdLoad-DEBUG %s %s %d %d %d %d\n",
    oname,fname,frame,type,finish,discrete
    ENDFD;
  if (ok) {
    APIEntry();
    origObj=ExecutiveFindObjectByName(oname);
    /* check for existing object of right type, delete if not */
    if(origObj) {
      new_type = -1;
      switch(type) {
      case cLoadTypeChemPyModel:
      case cLoadTypePDB:
      case cLoadTypeXYZ:
      case cLoadTypePDBStr:
      case cLoadTypeMOL:
      case cLoadTypeMOLStr:
      case cLoadTypeMMD:
      case cLoadTypeMMDSeparate:
      case cLoadTypeMMDStr:
      case cLoadTypePMO:
      case cLoadTypeTOP:
      case cLoadTypeTRJ:
      case cLoadTypeCRD:
        new_type = cObjectMolecule;
        break;
      case cLoadTypeChemPyBrick:
      case cLoadTypeChemPyMap:
      case cLoadTypeXPLORMap:
      case cLoadTypeCCP4Map:
      case cLoadTypeFLDMap:
      case cLoadTypeGRDMap:
        new_type = cObjectMap;
        break;
      case cLoadTypeCallback:
        new_type = cObjectCallback;
        break;
      case cLoadTypeCGO:
        new_type = cObjectCGO;
        break;
      }
      if (new_type!=origObj->type) {
        ExecutiveDelete(origObj->Name);
        origObj=NULL;
      }
    }
    
    
    switch(type) {
    case cLoadTypePDB:
      ExecutiveProcessPDBFile(origObj,fname,oname,frame,discrete,finish,buf,NULL);
      /* special handler for concatenated and annotated files */
      break;
    case cLoadTypePQR:
      {
        PDBInfoRec pdb_info;
        UtilZeroMem(&pdb_info,sizeof(PDBInfoRec));

        pdb_info.is_pqr_file = true;        
        ExecutiveProcessPDBFile(origObj,fname,oname,frame,discrete,finish,buf,&pdb_info);
      }
      break;
    case cLoadTypeTOP:
      PRINTFD(FB_CCmd) " CmdLoad-DEBUG: loading TOP\n" ENDFD;
      if(origObj) { /* always reinitialize topology objects from scratch */
        ExecutiveDelete(origObj->Name);
        origObj=NULL;
      }
      if(!origObj) {
        obj=(CObject*)ObjectMoleculeLoadTOPFile(NULL,fname,frame,discrete);
        if(obj) {
          ObjectSetName(obj,oname);
          ExecutiveManageObject(obj,false,quiet);
          if(frame<0)
            frame = ((ObjectMolecule*)obj)->NCSet-1;
          sprintf(buf," CmdLoad: \"%s\" loaded into object \"%s\", state %d.\n",
                  fname,oname,frame+1);
        }
      }
      break;
    case cLoadTypeTRJ:
      PRINTFD(FB_CCmd) " CmdLoad-DEBUG: loading TRJ\n" ENDFD;
      if(origObj) { /* always reinitialize topology objects from scratch */
        ObjectMoleculeLoadTRJFile((ObjectMolecule*)origObj,fname,frame,
                                  1,1,1,-1,-1,NULL,1,NULL);
        /* if(finish)
           ExecutiveUpdateObjectSelection(origObj); unnecc */
        sprintf(buf," CmdLoad: \"%s\" appended into object \"%s\".\n CmdLoad: %d total states in the object.\n",
                fname,oname,((ObjectMolecule*)origObj)->NCSet);
      } else {
        PRINTFB(FB_CCmd,FB_Errors)
          "CmdLoad-Error: must load object topology before loading trajectory!"
          ENDFB;
      }
      break;
    case cLoadTypeCRD:
      PRINTFD(FB_CCmd) " CmdLoad-DEBUG: loading CRD\n" ENDFD;
      if(origObj) { /* always reinitialize topology objects from scratch */
        ObjectMoleculeLoadRSTFile((ObjectMolecule*)origObj,fname,frame);
        /* if(finish)
           ExecutiveUpdateObjectSelection(origObj); unnecc */
        sprintf(buf," CmdLoad: \"%s\" appended into object \"%s\".\n CmdLoad: %d total states in the object.\n",
                fname,oname,((ObjectMolecule*)origObj)->NCSet);
      } else {
        PRINTFB(FB_CCmd,FB_Errors)
          "CmdLoad-Error: must load object topology before loading coordinate file!"
          ENDFB;
      }
      break;
    case cLoadTypeRST:
      PRINTFD(FB_CCmd) " CmdLoad-DEBUG: loading RST\n" ENDFD;
      if(origObj) { /* always reinitialize topology objects from scratch */
        ObjectMoleculeLoadRSTFile((ObjectMolecule*)origObj,fname,frame);
        /* if(finish)
           ExecutiveUpdateObjectSelection(origObj); unnecc */
        sprintf(buf," CmdLoad: \"%s\" appended into object \"%s\".\n CmdLoad: %d total states in the object.\n",
                fname,oname,((ObjectMolecule*)origObj)->NCSet);
      } else {
        PRINTFB(FB_CCmd,FB_Errors)
          "CmdLoad-Error: must load object topology before loading restart file!"
          ENDFB;
      }
      break;
    case cLoadTypePMO:
      PRINTFD(FB_CCmd) " CmdLoad-DEBUG: loading PMO\n" ENDFD;
      if(!origObj) {
        obj=(CObject*)ObjectMoleculeLoadPMOFile(NULL,fname,frame,discrete);
        if(obj) {
          ObjectSetName(obj,oname);
          ExecutiveManageObject(obj,true,quiet);
          if(frame<0)
            frame = ((ObjectMolecule*)obj)->NCSet-1;
          sprintf(buf," CmdLoad: \"%s\" loaded into object \"%s\", state %d.\n",
                  fname,oname,frame+1);
        }
      } else {
        ObjectMoleculeLoadPMOFile((ObjectMolecule*)origObj,fname,frame,discrete);
        if(finish)
          ExecutiveUpdateObjectSelection(origObj);
        if(frame<0)
          frame = ((ObjectMolecule*)origObj)->NCSet-1;
        sprintf(buf," CmdLoad: \"%s\" appended into object \"%s\", state %d.\n",
                fname,oname,frame+1);
      }
      break;
    case cLoadTypeXYZ:
      PRINTFD(FB_CCmd) " CmdLoad-DEBUG: loading XYZStr\n" ENDFD;
      if(!origObj) {
        obj=(CObject*)ObjectMoleculeLoadXYZFile(NULL,fname,frame,discrete);
        if(obj) {
          ObjectSetName(obj,oname);
          ExecutiveManageObject(obj,true,quiet);
          if(frame<0)
            frame = ((ObjectMolecule*)obj)->NCSet-1;
          sprintf(buf," CmdLoad: \"%s\" loaded into object \"%s\", state %d.\n",
                  fname,oname,frame+1);
        }
      } else {
        ObjectMoleculeLoadXYZFile((ObjectMolecule*)origObj,fname,frame,discrete);
        if(finish)
          ExecutiveUpdateObjectSelection(origObj);
        if(frame<0)
          frame = ((ObjectMolecule*)origObj)->NCSet-1;
        sprintf(buf," CmdLoad: \"%s\" appended into object \"%s\", state %d.\n",
                fname,oname,frame+1);
      }
      break;
    case cLoadTypePDBStr:
      PRINTFD(FB_CCmd) " CmdLoad-DEBUG: loading PDBStr\n" ENDFD;
      obj=(CObject*)ObjectMoleculeReadPDBStr((ObjectMolecule*)origObj,
                                             fname,frame,discrete,
                                             NULL,NULL,NULL,NULL);
      if(!origObj) {
        if(obj) {
          ObjectSetName(obj,oname);
          ExecutiveManageObject(obj,true,quiet);
          if(frame<0)
            frame = ((ObjectMolecule*)obj)->NCSet-1;
          sprintf(buf," CmdLoad: PDB-string loaded into object \"%s\", state %d.\n",
                  oname,frame+1);		  
        }
      } else if(origObj) {
        if(finish)
          ExecutiveUpdateObjectSelection(origObj);
        if(frame<0)
          frame = ((ObjectMolecule*)origObj)->NCSet-1;
        sprintf(buf," CmdLoad: PDB-string appended into object \"%s\", state %d.\n",
                oname,frame+1);
      }
      break;
    case cLoadTypeMOL:
      PRINTFD(FB_CCmd) " CmdLoad-DEBUG: loading MOL\n" ENDFD;
      obj=(CObject*)ObjectMoleculeLoadMOLFile((ObjectMolecule*)origObj,fname,frame,discrete);
      if(!origObj) {
        if(obj) {
          ObjectSetName(obj,oname);
          ExecutiveManageObject(obj,true,quiet);
          if(frame<0)
            frame = ((ObjectMolecule*)obj)->NCSet-1;
          sprintf(buf," CmdLoad: \"%s\" loaded into object \"%s\", state %d.\n",
                  fname,oname,frame+1);		  
        }
      } else if(origObj) {
        if(finish)
          ExecutiveUpdateObjectSelection(origObj);
        if(frame<0)
          frame = ((ObjectMolecule*)origObj)->NCSet-1;
        sprintf(buf," CmdLoad: \"%s\" appended into object \"%s\", state %d.\n",
                fname,oname,frame+1);
      }
      break;
    case cLoadTypeMOLStr:
      PRINTFD(FB_CCmd) " CmdLoad-DEBUG: reading MOLStr\n" ENDFD;
      obj=(CObject*)ObjectMoleculeReadMOLStr((ObjectMolecule*)origObj,fname,frame,discrete);
      if(!origObj) {
        if(obj) {
          ObjectSetName(obj,oname);
          ExecutiveManageObject(obj,true,false);
          if(frame<0)
            frame = ((ObjectMolecule*)obj)->NCSet-1;
          sprintf(buf," CmdLoad: MOL-string loaded into object \"%s\", state %d.\n",
                  oname,frame+1);		  
        }
      } else if(origObj) {
        if(finish)
          ExecutiveUpdateObjectSelection(origObj);
        if(frame<0)
          frame = ((ObjectMolecule*)origObj)->NCSet-1;
        sprintf(buf," CmdLoad: MOL-string appended into object \"%s\", state %d.\n",
                oname,frame+1);
      }
      break;
    case cLoadTypeMMD:
      PRINTFD(FB_CCmd) " CmdLoad-DEBUG: loading MMD\n" ENDFD;
      obj=(CObject*)ObjectMoleculeLoadMMDFile((ObjectMolecule*)origObj,fname,frame,NULL,discrete);
      if(!origObj) {
        if(obj) {
          ObjectSetName(obj,oname);
          ExecutiveManageObject(obj,true,quiet);
          if(frame<0)
            frame = ((ObjectMolecule*)obj)->NCSet-1;
          sprintf(buf," CmdLoad: \"%s\" loaded into object \"%s\", state %d.\n",
                  fname,oname,frame+1);		  
        }
      } else if(origObj) {
        if(finish)
          ExecutiveUpdateObjectSelection(origObj);
        if(frame<0)
          frame = ((ObjectMolecule*)origObj)->NCSet-1;
        sprintf(buf," CmdLoad: \"%s\" appended into object \"%s\", state %d.\n",
                fname,oname,frame+1);
      }
      break;
    case cLoadTypeMMDSeparate:
      ObjectMoleculeLoadMMDFile((ObjectMolecule*)origObj,fname,frame,oname,discrete);
      break;
    case cLoadTypeMMDStr:
      PRINTFD(FB_CCmd) " CmdLoad-DEBUG: loading MMDStr\n" ENDFD;
      obj=(CObject*)ObjectMoleculeReadMMDStr((ObjectMolecule*)origObj,fname,frame,discrete);
      if(!origObj) {
        if(obj) {
          ObjectSetName(obj,oname);
          ExecutiveManageObject(obj,true,quiet);
          if(frame<0)
            frame = ((ObjectMolecule*)obj)->NCSet-1;
          sprintf(buf," CmdLoad: MMD-string loaded into object \"%s\", state %d.\n",
                  oname,frame+1);		  
        }
      } else if(origObj) {
        if(finish)
          ExecutiveUpdateObjectSelection(origObj);
        if(frame<0)
          frame = ((ObjectMolecule*)origObj)->NCSet-1;
        sprintf(buf," CmdLoad: MMD-string appended into object \"%s\", state %d\n",
                oname,frame+1);
      }
      break;
    case cLoadTypeXPLORMap:
      PRINTFD(FB_CCmd) " CmdLoad-DEBUG: loading XPLORMap\n" ENDFD;
      if(!origObj) {
        obj=(CObject*)ObjectMapLoadXPLORFile(NULL,fname,frame,true);
        if(obj) {
          ObjectSetName(obj,oname);
          ExecutiveManageObject((CObject*)obj,true,quiet);
          sprintf(buf," CmdLoad: \"%s\" loaded into object \"%s\".\n",fname,oname);
        }
      } else {
        ObjectMapLoadXPLORFile((ObjectMap*)origObj,fname,frame,true);
        sprintf(buf," CmdLoad: \"%s\" appended into object \"%s\".\n",
                fname,oname);
      }
      break;
    case cLoadTypeXPLORStr:
      PRINTFD(FB_CCmd) " CmdLoad-DEBUG: loading XPLOR string\n" ENDFD;
      if(!origObj) {
        obj=(CObject*)ObjectMapLoadXPLORFile(NULL,fname,frame,false);
        if(obj) {
          ObjectSetName(obj,oname);
          ExecutiveManageObject((CObject*)obj,true,quiet);
          sprintf(buf," CmdLoad: XPLOR string loaded into object \"%s\".\n",oname);
        }
      } else {
        ObjectMapLoadXPLORFile((ObjectMap*)origObj,fname,frame,false);
        sprintf(buf," CmdLoad: XPLOR string appended into object \"%s\".\n",
                oname);
      }
      break;
    case cLoadTypeCCP4Map:
      PRINTFD(FB_CCmd) " CmdLoad-DEBUG: loading CCP4 map\n" ENDFD;
      if(!origObj) {
        obj=(CObject*)ObjectMapLoadCCP4File(NULL,fname,frame);
        if(obj) {
          ObjectSetName(obj,oname);
          ExecutiveManageObject((CObject*)obj,true,quiet);
          sprintf(buf," CmdLoad: \"%s\" loaded into object \"%s\".\n",fname,oname);
        }
      } else {
        ObjectMapLoadCCP4File((ObjectMap*)origObj,fname,frame);
        sprintf(buf," CmdLoad: \"%s\" appended into object \"%s\".\n",
                fname,oname);
      }
      break;
    case cLoadTypePHIMap:
      PRINTFD(FB_CCmd) " CmdLoad-DEBUG: loading Delphi Map\n" ENDFD;
      if(!origObj) {
        obj=(CObject*)ObjectMapLoadPHIFile(NULL,fname,frame);
        if(obj) {
          ObjectSetName(obj,oname);
          ExecutiveManageObject((CObject*)obj,true,quiet);
          sprintf(buf," CmdLoad: \"%s\" loaded into object \"%s\".\n",fname,oname);
        }
      } else {
        ObjectMapLoadPHIFile((ObjectMap*)origObj,fname,frame);
        sprintf(buf," CmdLoad: \"%s\" appended into object \"%s\".\n",
                fname,oname);
      }
      break;
    case cLoadTypeDXMap:
      PRINTFD(FB_CCmd) " CmdLoad-DEBUG: loading DX Map\n" ENDFD;
      if(!origObj) {
        obj=(CObject*)ObjectMapLoadDXFile(NULL,fname,frame);
        if(obj) {
          ObjectSetName(obj,oname);
          ExecutiveManageObject((CObject*)obj,true,quiet);
          sprintf(buf," CmdLoad: \"%s\" loaded into object \"%s\".\n",fname,oname);
        }
      } else {
        ObjectMapLoadDXFile((ObjectMap*)origObj,fname,frame);
        sprintf(buf," CmdLoad: \"%s\" appended into object \"%s\".\n",
                fname,oname);
      }
      break;
    case cLoadTypeFLDMap:
      PRINTFD(FB_CCmd) " CmdLoad-DEBUG: loading AVS Map\n" ENDFD;
      if(!origObj) {
        obj=(CObject*)ObjectMapLoadFLDFile(NULL,fname,frame);
        if(obj) {
          ObjectSetName(obj,oname);
          ExecutiveManageObject((CObject*)obj,true,quiet);
          sprintf(buf," CmdLoad: \"%s\" loaded into object \"%s\".\n",fname,oname);
        }
      } else {
        ObjectMapLoadFLDFile((ObjectMap*)origObj,fname,frame);
        sprintf(buf," CmdLoad: \"%s\" appended into object \"%s\".\n",
                fname,oname);
      }
      break;
    case cLoadTypeBRIXMap:
      PRINTFD(FB_CCmd) " CmdLoad-DEBUG: loading BRIX/DSN6 Map\n" ENDFD;
      if(!origObj) {
        obj=(CObject*)ObjectMapLoadBRIXFile(NULL,fname,frame);
        if(obj) {
          ObjectSetName(obj,oname);
          ExecutiveManageObject((CObject*)obj,true,quiet);
          sprintf(buf," CmdLoad: \"%s\" loaded into object \"%s\".\n",fname,oname);
        }
      } else {
        ObjectMapLoadFLDFile((ObjectMap*)origObj,fname,frame);
        sprintf(buf," CmdLoad: \"%s\" appended into object \"%s\".\n",
                fname,oname);
      }
      break;
    case cLoadTypeGRDMap:
      PRINTFD(FB_CCmd) " CmdLoad-DEBUG: loading GRD Map\n" ENDFD;
      if(!origObj) {
        obj=(CObject*)ObjectMapLoadGRDFile(NULL,fname,frame);
        if(obj) {
          ObjectSetName(obj,oname);
          ExecutiveManageObject((CObject*)obj,true,quiet);
          sprintf(buf," CmdLoad: \"%s\" loaded into object \"%s\".\n",fname,oname);
        }
      } else {
        ObjectMapLoadGRDFile((ObjectMap*)origObj,fname,frame);
        sprintf(buf," CmdLoad: \"%s\" appended into object \"%s\".\n",
                fname,oname);
      }
      break;
    }
    if(origObj) {
      PRINTFB(FB_Executive,FB_Actions) 
        "%s",buf
        ENDFB;
      OrthoRestorePrompt();
    }
    APIExit();
  }
  return(APIStatus(ok));
}


static PyObject *CmdLoadTraj(PyObject *self, PyObject *args)
{
  char *fname,*oname;
  CObject *origObj = NULL;
  OrthoLineType buf;
  int frame,type;
  int new_type;
  int interval,average,start,stop,max,image;
  OrthoLineType s1;
  char *str1;
  int ok=false;
  float shift[3];

  ok = PyArg_ParseTuple(args,"ssiiiiiiisifff",&oname,&fname,&frame,&type,
                        &interval,&average,&start,&stop,&max,&str1,
                        &image,&shift[0],&shift[1],&shift[2]);

  buf[0]=0;
  if (ok) {
    APIEntry();
    if(str1[0])
      SelectorGetTmp(str1,s1);
    else
      s1[0]=0; /* no selection */
    origObj=ExecutiveFindObjectByName(oname);
    /* check for existing object of right type, delete if not */
    if(origObj) {
      new_type = -1;
      switch(type) {
      case cLoadTypeTRJ:
        new_type = cObjectMolecule;
        break;
      }
      if (new_type!=origObj->type) {
        ExecutiveDelete(origObj->Name);
        origObj=NULL;
      }
    }
    
    switch(type) {
    case cLoadTypeTRJ:
      PRINTFD(FB_CCmd) " CmdLoadTraj-DEBUG: loading TRJ\n" ENDFD;
      if(origObj) { /* always reinitialize topology objects from scratch */
        ObjectMoleculeLoadTRJFile((ObjectMolecule*)origObj,fname,frame,
                                  interval,average,start,stop,max,s1,image,shift);
        /* if(finish)
           ExecutiveUpdateObjectSelection(origObj); unnecc */
        sprintf(buf," CmdLoadTraj: \"%s\" appended into object \"%s\".\n CmdLoadTraj: %d total states in the object.\n",
                fname,oname,((ObjectMolecule*)origObj)->NCSet);
      } else {
        PRINTFB(FB_CCmd,FB_Errors)
          "CmdLoadTraj-Error: must load object topology before loading trajectory!\n"
          ENDFB;
      }
      break;
    }
    if(origObj) {
      PRINTFB(FB_Executive,FB_Actions) 
        "%s",buf
        ENDFB;
      OrthoRestorePrompt();
    }
    SelectorFreeTmp(s1);
    APIExit();
  }
  return(APIStatus(ok));
}

static PyObject *CmdOrigin(PyObject *self, PyObject *args)
{
  char *str1,*obj;
  OrthoLineType s1;
  float v[3];
  int ok=false;
  int state;
  ok = PyArg_ParseTuple(args,"ss(fff)i",&str1,&obj,v,v+1,v+2,&state);
  if (ok) {
    APIEntry();
    if(str1[0])
      SelectorGetTmp(str1,s1);
    else
      s1[0]=0; /* no selection */
    ok = ExecutiveOrigin(s1,1,obj,v,state); /* TODO STATUS */
    if(str1[0])
      SelectorFreeTmp(s1);
    APIExit();
  }
  return(APIStatus(ok));
}

static PyObject *CmdSort(PyObject *self, PyObject *args)
{
  char *name;
  int ok=false;
  ok = PyArg_ParseTuple(args,"s",&name);
  if (ok) {
    APIEntry();
    ExecutiveSort(name); /* TODO STATUS */
    APIExit();
  }
  return(APIStatus(ok));
}

static PyObject *CmdAssignSS(PyObject *self, PyObject *args)
/* EXPERIMENTAL */
{
  int ok=false;
  int state,quiet;
  char *str1,*str2;
  int preserve;
  OrthoLineType s1,s2;
  ok = PyArg_ParseTuple(args,"sisii",&str1,&state,&str2,&preserve,&quiet);
  if (ok) {
    APIEntry();
    SelectorGetTmp(str1,s1);
    SelectorGetTmp(str2,s2);
    ok = ExecutiveAssignSS(s1,state,s2,preserve,quiet);
    APIExit();
  }
  return(APIStatus(ok));
}

static PyObject *CmdSpheroid(PyObject *self, PyObject *args)
/* EXPERIMENTAL */
{
  char *name;
  int ok=false;
  int average;
  ok = PyArg_ParseTuple(args,"si",&name,&average);
  if (ok) {
    APIEntry();
    ExecutiveSpheroid(name,average); /* TODO STATUS */
    APIExit();
  }
  return(APIStatus(ok));
}

static PyObject *CmdTest(PyObject *self, PyObject *args)
{
  /* regression tests */

  int ok=true;
  int code;
  int group;
  CTestPyMOL tst;

  ok = PyArg_ParseTuple(args,"ii",&group,&code);
  if(ok) {
    APIEntry();
    PRINTFB(FB_CCmd,FB_Details)
      " Cmd: initiating test %d-%d.\n",group,code
      ENDFB;
    ok = TestPyMOLRun(&tst,group,code);
    PRINTFB(FB_CCmd,FB_Details)
      " Cmd: concluding test %d-%d.\n",group,code
      ENDFB;
    APIExit();
  }
  return(APIStatus(ok));
}

static PyObject *CmdCenter(PyObject *self, PyObject *args)
{
  char *str1;
  OrthoLineType s1;
  int state;
  int origin;
  int ok=false;
  ok = PyArg_ParseTuple(args,"sii",&str1,&state,&origin);
  if (ok) {
    APIEntry();
    SelectorGetTmp(str1,s1);
    ok = ExecutiveCenter(s1,state,origin);
    SelectorFreeTmp(s1);
    APIExit();
  }
  return(APIStatus(ok));
}

static PyObject *CmdZoom(PyObject *self, PyObject *args)
{
  char *str1;
  OrthoLineType s1;
  float buffer;
  int state;
  int inclusive;
  int ok=false;
  ok = PyArg_ParseTuple(args,"sfii",&str1,&buffer,&state,&inclusive);
  if (ok) {
    APIEntry();
    SelectorGetTmp(str1,s1);
    ok = ExecutiveWindowZoom(s1,buffer,state,inclusive); 
    SelectorFreeTmp(s1);
    APIExit();
  }
  return(APIStatus(ok));
}

static PyObject *CmdIsolevel(PyObject *self, PyObject *args)
{
  float level;
  int state;
  char *name;
  int ok=false;
  ok = PyArg_ParseTuple(args,"sfi",&name,&level,&state);
  if (ok) {
    APIEntry();
    ok = ExecutiveIsolevel(name,level,state);
    APIExit();
  }
  return(APIStatus(ok));
}

static PyObject *CmdHAdd(PyObject *self, PyObject *args)
{
  char *str1;
  OrthoLineType s1;

  int ok=false;
  ok = PyArg_ParseTuple(args,"s",&str1);
  if (ok) {
    APIEntry();
    SelectorGetTmp(str1,s1);
    ExecutiveAddHydrogens(s1); /* TODO STATUS */
    SelectorFreeTmp(s1);
    SelectorGetTmp(str1,s1);
    ExecutiveAddHydrogens(s1); /* TODO STATUS */
    SelectorFreeTmp(s1);
    SelectorGetTmp(str1,s1);
    ExecutiveAddHydrogens(s1); /* TODO STATUS */
    SelectorFreeTmp(s1);
    SelectorGetTmp(str1,s1);
    ExecutiveAddHydrogens(s1); /* TODO STATUS */
    SelectorFreeTmp(s1);
    APIExit();
  }
  return(APIStatus(ok));
}

static PyObject *CmdGetObjectColorIndex(PyObject *self, PyObject *args)
{
  char *str1;
  int result = -1;
  int ok=false;
  ok = PyArg_ParseTuple(args,"s",&str1);
  if (ok) {
    APIEntry();
    result = ExecutiveGetObjectColorIndex(str1);
    APIExit();
  }
  return(APIStatus(result));
}

static PyObject *CmdRemove(PyObject *self, PyObject *args)
{
  char *str1;
  OrthoLineType s1;

  int ok=false;
  ok = PyArg_ParseTuple(args,"s",&str1);
  if (ok) {
    APIEntry();
    SelectorGetTmp(str1,s1);
    ExecutiveRemoveAtoms(s1); /* TODO STATUS */
    SelectorFreeTmp(s1);
    APIExit();
  }
  return(APIStatus(ok));
}

static PyObject *CmdRemovePicked(PyObject *self, PyObject *args)
{
  int i1;
  int ok=false;
  ok = PyArg_ParseTuple(args,"i",&i1);
  if (ok) {
    APIEntry();
    EditorRemove(i1); /* TODO STATUS */
    APIExit();
  }
  return(APIStatus(ok));
}

static PyObject *CmdHFill(PyObject *self, PyObject *args)
{
  int ok = true;
  APIEntry();
  EditorHFill(); /* TODO STATUS */
  APIExit();
  return(APIStatus(ok));
}

static PyObject *CmdCycleValence(PyObject *self, PyObject *args)
{
  int ok = true;
  APIEntry();
  EditorCycleValence();  /* TODO STATUS */
  APIExit();
  return(APIStatus(ok));
}

static PyObject *CmdReplace(PyObject *self, 	PyObject *args)
{
  int i1,i2;
  char *str1,*str2;
  int ok=false;
  ok = PyArg_ParseTuple(args,"siis",&str1,&i1,&i2,&str2);
  if (ok) {
    APIEntry();
    EditorReplace(str1,i1,i2,str2);  /* TODO STATUS */
    APIExit();
  }
  return(APIStatus(ok));  
}

static PyObject *CmdSetGeometry(PyObject *self, 	PyObject *args)
{
  int i1,i2;
  char *str1;
  OrthoLineType s1;
  int ok=false;
  ok = PyArg_ParseTuple(args,"sii",&str1,&i1,&i2);
  if (ok) {
    APIEntry();
    SelectorGetTmp(str1,s1);
    ok = ExecutiveSetGeometry(s1,i1,i2);  /* TODO STATUS */
    SelectorFreeTmp(s1);
    APIExit();
  }
  return(APIStatus(ok));  
}

static PyObject *CmdAttach(PyObject *self, 	PyObject *args)
{
  int i1,i2;
  char *str1;
  int ok=false;
  char *name;
  ok = PyArg_ParseTuple(args,"siis",&str1,&i1,&i2,&name);
  if (ok) {
    APIEntry();
    EditorAttach(str1,i1,i2,name);  /* TODO STATUS */
    APIExit();
  }
  return(APIStatus(ok));  
}

static PyObject *CmdFuse(PyObject *self, 	PyObject *args)
{
  char *str1,*str2;
  int mode;
  OrthoLineType s1,s2;

  int ok=false;
  ok = PyArg_ParseTuple(args,"ssi",&str1,&str2,&mode);
  if (ok) {
    APIEntry();
    SelectorGetTmp(str1,s1);
    SelectorGetTmp(str2,s2);
    ExecutiveFuse(s1,s2,mode);  /* TODO STATUS */
    SelectorFreeTmp(s1);
    SelectorFreeTmp(s2);
    APIExit();
  }
  return(APIStatus(ok));  
}

static PyObject *CmdUnpick(PyObject *self, 	PyObject *args)
{
  APIEntry();
  EditorInactivate();  /* TODO STATUS */
  APIExit();
  return(APISuccess());  
}

static PyObject *CmdEdit(PyObject *self, 	PyObject *args)
{
  char *str0,*str1,*str2,*str3;
  OrthoLineType s0 = "";
  OrthoLineType s1 = "";
  OrthoLineType s2 = "";
  OrthoLineType s3 = "";
  int pkresi,pkbond;
  int ok=false;
  int quiet;
  ok = PyArg_ParseTuple(args,"ssssiii",&str0,&str1,&str2,&str3,&pkresi,&pkbond,&quiet);
  if (ok) {
    APIEntry();
    if(str0[0]) SelectorGetTmp(str0,s0);
    if(str1[0]) SelectorGetTmp(str1,s1);
    if(str2[0]) SelectorGetTmp(str2,s2);
    if(str3[0]) SelectorGetTmp(str3,s3);
    ok = EditorSelect(s0,s1,s2,s3,pkresi,pkbond,quiet);
    if(s0[0]) SelectorFreeTmp(s0);
    if(s1[0]) SelectorFreeTmp(s1);
    if(s2[0]) SelectorFreeTmp(s2);
    if(s3[0]) SelectorFreeTmp(s3);
    
    APIExit();
  }
  return APIStatus(ok);
}

static PyObject *CmdRename(PyObject *self, 	PyObject *args)
{
  char *str1;
  int int1;
  OrthoLineType s1;

  int ok=false;
  ok = PyArg_ParseTuple(args,"si",&str1,&int1);
  if (ok) {
    APIEntry();
    SelectorGetTmp(str1,s1);
    ExecutiveRenameObjectAtoms(s1,int1); /* TODO STATUS */
    SelectorFreeTmp(s1);
    APIExit();
  }
  return(APIStatus(ok));  
}

static PyMethodDef Cmd_methods[] = {
	{"accept",	              CmdAccept,               METH_VARARGS },
	{"align",	              CmdAlign,                METH_VARARGS },
	{"alter",	              CmdAlter,                METH_VARARGS },
	{"alter_state",           CmdAlterState,           METH_VARARGS },
	{"attach",                CmdAttach,               METH_VARARGS },
   {"bg_color",              CmdBackgroundColor,      METH_VARARGS },
	{"bond",                  CmdBond,                 METH_VARARGS },
   {"button",                CmdButton,               METH_VARARGS },
   {"cartoon",               CmdCartoon,              METH_VARARGS },
   {"center",                CmdCenter,               METH_VARARGS },
	{"clip",	                 CmdClip,                 METH_VARARGS },
	{"cls",	                 CmdCls,                  METH_VARARGS },
	{"color",	              CmdColor,                METH_VARARGS },
	{"colordef",	           CmdColorDef,             METH_VARARGS },
   {"combine_object_ttt",    CmdCombineObjectTTT,     METH_VARARGS },
	{"copy",                  CmdCopy,                 METH_VARARGS },
	{"create",                CmdCreate,               METH_VARARGS },
	{"count_states",          CmdCountStates,          METH_VARARGS },
	{"count_frames",          CmdCountFrames,          METH_VARARGS },
	{"cycle_valence",         CmdCycleValence,         METH_VARARGS },
   {"debug",                 CmdDebug,                METH_VARARGS },
   {"decline",               CmdDecline,              METH_VARARGS },
   {"del_colorection",       CmdDelColorection,       METH_VARARGS },   
   {"fake_drag",             CmdFakeDrag,             METH_VARARGS },   
   {"gl_delete_lists",       CmdGLDeleteLists,        METH_VARARGS },
	{"delete",                CmdDelete,               METH_VARARGS },
	{"dirty",                 CmdDirty,                METH_VARARGS },
	{"distance",	           CmdDistance,             METH_VARARGS },
	{"dist",    	           CmdDist,                 METH_VARARGS },
	{"do",	                 CmdDo,                   METH_VARARGS },
	{"dump",	                 CmdDump,                 METH_VARARGS },
   {"edit",                  CmdEdit,                 METH_VARARGS },
   {"torsion",               CmdTorsion,              METH_VARARGS },
	{"export_dots",           CmdExportDots,           METH_VARARGS },
	{"export_coords",         CmdExportCoords,         METH_VARARGS },
	{"feedback",              CmdFeedback,             METH_VARARGS },
	{"find_pairs",            CmdFindPairs,            METH_VARARGS },
	{"finish_object",         CmdFinishObject,         METH_VARARGS },
	{"fit",                   CmdFit,                  METH_VARARGS },
	{"fit_pairs",             CmdFitPairs,             METH_VARARGS },
	{"fix_chemistry",         CmdFixChemistry,         METH_VARARGS },
	{"flag",                  CmdFlag,                 METH_VARARGS },
	{"frame",	              CmdFrame,                METH_VARARGS },
   {"flush_now",             CmdFlushNow,             METH_VARARGS },
   {"focus",                 CmdFocus,                METH_VARARGS },
   {"delete_colorection",    CmdDelColorection,       METH_VARARGS },   
	{"dss",   	              CmdAssignSS,             METH_VARARGS },
   {"full_screen",           CmdFullScreen,           METH_VARARGS },
   {"fuse",                  CmdFuse,                 METH_VARARGS },
	{"get",	                 CmdGet,                  METH_VARARGS },
	{"get_area",              CmdGetArea,              METH_VARARGS },
   {"get_atom_coords",       CmdGetAtomCoords,        METH_VARARGS },
   {"get_bond_print",        CmdGetBondPrint,         METH_VARARGS },
	{"get_chains",            CmdGetChains,            METH_VARARGS },
	{"get_color",             CmdGetColor,             METH_VARARGS },
   {"get_colorection",       CmdGetColorection,       METH_VARARGS },   
	{"get_dihe",              CmdGetDihe,              METH_VARARGS },
	{"get_frame",             CmdGetFrame,             METH_VARARGS },
	{"get_feedback",          CmdGetFeedback,          METH_VARARGS },
	{"get_matrix",	           CmdGetMatrix,            METH_VARARGS },
	{"get_min_max",           CmdGetMinMax,            METH_VARARGS },
	{"get_model",	           CmdGetModel,             METH_VARARGS },
	{"get_moment",	           CmdGetMoment,            METH_VARARGS },
   {"get_movie_locked",      CmdGetMovieLocked,       METH_VARARGS },
   {"get_names",             CmdGetNames,             METH_VARARGS },
   {"get_object_color_index",CmdGetObjectColorIndex,  METH_VARARGS },
	{"get_position",	        CmdGetPosition,          METH_VARARGS },
	{"get_povray",	           CmdGetPovRay,            METH_VARARGS },
	{"get_pdb",	              CmdGetPDB,               METH_VARARGS },
   {"get_phipsi",            CmdGetPhiPsi,            METH_VARARGS },
   {"get_renderer",          CmdGetRenderer,          METH_VARARGS },
   {"get_session",           CmdGetSession,           METH_VARARGS },
	{"get_setting",           CmdGetSetting,           METH_VARARGS },
	{"get_setting_tuple",     CmdGetSettingTuple,      METH_VARARGS },
	{"get_setting_text",      CmdGetSettingText,       METH_VARARGS },
   {"get_setting_updates",   CmdGetSettingUpdates,    METH_VARARGS },
   {"get_symmetry",          CmdGetCrystal,           METH_VARARGS },
	{"get_state",             CmdGetState,             METH_VARARGS },
   {"get_title",             CmdGetTitle,             METH_VARARGS },
	{"get_type",              CmdGetType,              METH_VARARGS },
   {"get_view",              CmdGetView,              METH_VARARGS },
   {"get_vis",               CmdGetVis,               METH_VARARGS },
   {"get_wizard",            CmdGetWizard,            METH_VARARGS },
   {"get_wizard_stack",      CmdGetWizardStack,       METH_VARARGS },
	{"h_add",                 CmdHAdd,                 METH_VARARGS },
	{"h_fill",                CmdHFill,                METH_VARARGS },
   {"identify",              CmdIdentify,             METH_VARARGS },
	{"import_coords",         CmdImportCoords,         METH_VARARGS },
   {"index",                 CmdIndex,                METH_VARARGS },
	{"intrafit",              CmdIntraFit,             METH_VARARGS },
   {"invert",                CmdInvert,               METH_VARARGS },
	{"isolevel",              CmdIsolevel,             METH_VARARGS },
	{"isomesh",	              CmdIsomesh,              METH_VARARGS },
	{"isosurface",	           CmdIsosurface,           METH_VARARGS },
   {"wait_queue",            CmdWaitQueue,            METH_VARARGS },
   {"label",                 CmdLabel,                METH_VARARGS },
	{"load",	                 CmdLoad,                 METH_VARARGS },
	{"load_color_table",	     CmdLoadColorTable,       METH_VARARGS },
	{"load_coords",           CmdLoadCoords,           METH_VARARGS },
	{"load_png",              CmdLoadPNG,              METH_VARARGS },
	{"load_object",           CmdLoadObject,           METH_VARARGS },
	{"load_traj",             CmdLoadTraj,             METH_VARARGS },
   {"map_new",               CmdMapNew,               METH_VARARGS },
   {"map_double",            CmdMapDouble,            METH_VARARGS },
   {"map_set_border",        CmdMapSetBorder,         METH_VARARGS },
	{"mask",	                 CmdMask,                 METH_VARARGS },
	{"mclear",	              CmdMClear,               METH_VARARGS },
	{"mdo",	                 CmdMDo,                  METH_VARARGS },
	{"mdump",	              CmdMDump,                METH_VARARGS },
	{"mem",	                 CmdMem,                  METH_VARARGS },
	{"move",	                 CmdMove,                 METH_VARARGS },
	{"mset",	                 CmdMSet,                 METH_VARARGS },
	{"mplay",	              CmdMPlay,                METH_VARARGS },
	{"mpng_",	              CmdMPNG,                 METH_VARARGS },
	{"mmatrix",	              CmdMMatrix,              METH_VARARGS },
	{"multisave",             CmdMultiSave,            METH_VARARGS },
	{"mview",	              CmdMView,                METH_VARARGS },
	{"origin",	              CmdOrigin,               METH_VARARGS },
	{"orient",	              CmdOrient,               METH_VARARGS },
	{"onoff",                 CmdOnOff,                METH_VARARGS },
   {"onoff_by_sele",         CmdOnOffBySele,          METH_VARARGS },
	{"overlap",               CmdOverlap,              METH_VARARGS },
	{"p_glut_event",          CmdPGlutEvent,           METH_VARARGS },
   {"p_glut_get_redisplay",  CmdPGlutGetRedisplay,    METH_VARARGS },
	{"paste",	              CmdPaste,                METH_VARARGS },
	{"png",	                 CmdPNG,                  METH_VARARGS },
	{"protect",	              CmdProtect,              METH_VARARGS },
	{"push_undo",	           CmdPushUndo,             METH_VARARGS },
	{"quit",	                 CmdQuit,                 METH_VARARGS },
   {"ray_trace_thread",      CmdRayTraceThread,       METH_VARARGS },
   {"ray_hash_thread",       CmdRayHashThread,        METH_VARARGS },
   {"ray_anti_thread",       CmdRayAntiThread,        METH_VARARGS },
   {"ramp_new",              CmdRampMapNew,           METH_VARARGS },
	{"ready",                 CmdReady,                METH_VARARGS },
   {"rebuild",               CmdRebuild,              METH_VARARGS },
   {"recolor",               CmdRecolor,              METH_VARARGS },
	{"refresh",               CmdRefresh,              METH_VARARGS },
	{"refresh_now",           CmdRefreshNow,           METH_VARARGS },
	{"refresh_wizard",        CmdRefreshWizard,        METH_VARARGS },
	{"remove",	              CmdRemove,               METH_VARARGS },
	{"remove_picked",         CmdRemovePicked,         METH_VARARGS },
	{"render",	              CmdRay,                  METH_VARARGS },
   {"rename",                CmdRename,               METH_VARARGS },
   {"replace",               CmdReplace,              METH_VARARGS },
   {"reinitialize",          CmdReinitialize,         METH_VARARGS },
	{"reset",                 CmdReset,                METH_VARARGS },
	{"reset_rate",	           CmdResetRate,            METH_VARARGS },
	{"rock",	                 CmdRock,                 METH_VARARGS },
	{"runpymol",	           CmdRunPyMOL,             METH_VARARGS },
	{"runwxpymol",	           CmdRunWXPyMOL,           METH_VARARGS },
	{"select",                CmdSelect,               METH_VARARGS },
	{"set",	                 CmdSet,                  METH_VARARGS },
	{"legacy_set",            CmdLegacySet,            METH_VARARGS },

	{"sculpt_deactivate",     CmdSculptDeactivate,     METH_VARARGS },
	{"sculpt_activate",       CmdSculptActivate,       METH_VARARGS },
	{"sculpt_iterate",        CmdSculptIterate,        METH_VARARGS },
	{"sculpt_purge",          CmdSculptPurge,          METH_VARARGS },
   {"set_colorection",       CmdSetColorection,       METH_VARARGS },   
	{"set_dihe",              CmdSetDihe,              METH_VARARGS },
	{"set_dihe",              CmdSetDihe,              METH_VARARGS },

	{"set_dihe",              CmdSetDihe,              METH_VARARGS },
	{"set_feedback",          CmdSetFeedbackMask,      METH_VARARGS },
	{"set_name",              CmdSetName,              METH_VARARGS },
   {"set_geometry",          CmdSetGeometry,          METH_VARARGS },
   {"set_session",           CmdSetSession,           METH_VARARGS },
   {"set_symmetry",          CmdSetCrystal,           METH_VARARGS },
	{"set_title",             CmdSetTitle,             METH_VARARGS },
	{"set_wizard",            CmdSetWizard,            METH_VARARGS },
	{"set_wizard_stack",      CmdSetWizardStack,       METH_VARARGS },
   {"set_view",              CmdSetView,              METH_VARARGS },
   {"set_vis",               CmdSetVis,               METH_VARARGS },
	{"setframe",	           CmdSetFrame,             METH_VARARGS },
	{"showhide",              CmdShowHide,             METH_VARARGS },
	{"set_matrix",	           CmdSetMatrix,            METH_VARARGS },
	{"smooth",	              CmdSmooth,               METH_VARARGS },
	{"sort",                  CmdSort,                 METH_VARARGS },
   {"spectrum",              CmdSpectrum,             METH_VARARGS },
   {"spheroid",              CmdSpheroid,             METH_VARARGS },
	{"splash",                CmdSplash,               METH_VARARGS },
	{"stereo",	              CmdStereo,               METH_VARARGS },
	{"system",	              CmdSystem,               METH_VARARGS },
	{"symexp",	              CmdSymExp,               METH_VARARGS },
	{"test",	                 CmdTest,                 METH_VARARGS },
	{"toggle",                CmdToggle,               METH_VARARGS },
	{"transform_object",      CmdTransformObject,      METH_VARARGS },
	{"transform_selection",   CmdTransformSelection,   METH_VARARGS },
	{"translate_atom",        CmdTranslateAtom,        METH_VARARGS },
	{"turn",	                 CmdTurn,                 METH_VARARGS },
	{"viewport",              CmdViewport,             METH_VARARGS },
	{"undo",                  CmdUndo,                 METH_VARARGS },
	{"unpick",                CmdUnpick,               METH_VARARGS },
	{"unset",                 CmdUnset,                METH_VARARGS },
	{"update",                CmdUpdate,               METH_VARARGS },
	{"zoom",	                 CmdZoom,                 METH_VARARGS },
	{NULL,		              NULL}     /* sentinel */        
};


void init_cmd(void)
{
  Py_InitModule("_cmd", Cmd_methods);
}

