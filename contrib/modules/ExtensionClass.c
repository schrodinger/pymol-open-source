/*****************************************************************************

  Copyright (c) 1996-2002 Zope Corporation and Contributors.
  All Rights Reserved.

  This software is subject to the provisions of the Zope Public License,
  Version 2.0 (ZPL).  A copy of the ZPL should accompany this distribution.
  THIS SOFTWARE IS PROVIDED "AS IS" AND ANY AND ALL EXPRESS OR IMPLIED
  WARRANTIES ARE DISCLAIMED, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
  WARRANTIES OF TITLE, MERCHANTABILITY, AGAINST INFRINGEMENT, AND FITNESS
  FOR A PARTICULAR PURPOSE

 ****************************************************************************/

static char ExtensionClass_module_documentation[] = 
"ExtensionClass - Classes implemented in c\n"
"\n"
"Built-in C classes are like Built-in types except that\n"
"they provide some of the behavior of Python classes:\n"
"\n"
"  - They provide access to unbound methods,\n"
"  - They can be called to create instances.\n"
"\n"
"$Id$\n"
;

#include "ExtensionClass.h"

static void
PyVar_Assign(PyObject **v,  PyObject *e)
{
  Py_XDECREF(*v);
  *v=e;
}

#define ASSIGN(V,E) PyVar_Assign(&(V),(E))
#define UNLESS(E) if (!(E))
#define UNLESS_ASSIGN(V,E) ASSIGN(V,E); UNLESS(V)
#define OBJECT(O) ((PyObject*)O)

/* Declarations for objects of type ExtensionClass */

staticforward PyExtensionClass ECType;

#define ExtensionClass_Check(O) ((O)->ob_type == (PyTypeObject*)&ECType)
#define ExtensionInstance_Check(O) \
   ((O)->ob_type->ob_type == (PyTypeObject*)&ECType)
#define AsExtensionClass(O) ((PyExtensionClass*)(O))
#define ExtensionClassOf(O) ((PyExtensionClass*)((O)->ob_type))
#define AsPyObject(O) ((PyObject*)(O))
#define NeedsToBeBound(O) \
   ((O)->ob_type->ob_type == (PyTypeObject*)&ECType && \
    (((PyExtensionClass*)((O)->ob_type))->class_flags & \
     EXTENSIONCLASS_BINDABLE_FLAG))
#define HasMethodHook(O) \
   ((O)->ob_type->ob_type == (PyTypeObject*)&ECType && \
    (((PyExtensionClass*)((O)->ob_type))->class_flags & \
     EXTENSIONCLASS_METHODHOOK_FLAG))

#define ALLOC_FREE(T) \
  if (free ## T) { \
      self=free ## T; \
      free ## T=(T*)self->self; \
      _Py_NewReference((PyObject *)self); \
      assert(self->ob_refcnt == 1); \
    } \
  else UNLESS(self = PyObject_NEW(T, & T ## Type)) return NULL;

#define METH_BY_NAME (2 << 16)

static PyObject *py__add__, *py__sub__, *py__mul__, *py__div__,
  *py__mod__, *py__pow__, *py__divmod__, *py__lshift__, *py__rshift__,
  *py__and__, *py__or__, *py__xor__, *py__coerce__, *py__neg__,
  *py__pos__, *py__abs__, *py__nonzero__, *py__inv__, *py__int__,
  *py__long__, *py__float__, *py__oct__, *py__hex__,
  *py__of__, *py__call__, *py__call_method__,
  *py__getitem__, *py__setitem__, *py__delitem__,
  *py__getslice__, *py__setslice__, *py__delslice__, *py__len__,
  *py__getattr__, *py__setattr__, *py__delattr__,
  *py__del__, *py__repr__, *py__str__, *py__class__, *py__name__,
  *py__hash__, *py__cmp__, *py__var_size__, *py__init__, *py__getinitargs__,
  *py__getstate__, *py__setstate__, *py__dict__, *pyclass_,
  *py__module__;

static PyObject *concat_fmt=0;
static PyObject *subclass_watcher=0;  /* Object that subclass events */

static void
init_py_names(void)
{
#define INIT_PY_NAME(N) py ## N = PyString_FromString(#N)
  INIT_PY_NAME(__add__);
  INIT_PY_NAME(__sub__);
  INIT_PY_NAME(__mul__);
  INIT_PY_NAME(__div__);
  INIT_PY_NAME(__mod__);
  INIT_PY_NAME(__pow__);
  INIT_PY_NAME(__divmod__);
  INIT_PY_NAME(__lshift__);
  INIT_PY_NAME(__rshift__);
  INIT_PY_NAME(__and__);
  INIT_PY_NAME(__or__);
  INIT_PY_NAME(__xor__);
  INIT_PY_NAME(__coerce__);
  INIT_PY_NAME(__neg__);
  INIT_PY_NAME(__pos__);
  INIT_PY_NAME(__abs__);
  INIT_PY_NAME(__nonzero__);
  INIT_PY_NAME(__inv__);
  INIT_PY_NAME(__int__);
  INIT_PY_NAME(__long__);
  INIT_PY_NAME(__float__);
  INIT_PY_NAME(__oct__);
  INIT_PY_NAME(__hex__);
  INIT_PY_NAME(__getitem__);
  INIT_PY_NAME(__setitem__);
  INIT_PY_NAME(__delitem__);
  INIT_PY_NAME(__getslice__);
  INIT_PY_NAME(__setslice__);
  INIT_PY_NAME(__delslice__);
  INIT_PY_NAME(__len__);
  INIT_PY_NAME(__of__);
  INIT_PY_NAME(__call__);
  INIT_PY_NAME(__call_method__);
  INIT_PY_NAME(__getattr__);
  INIT_PY_NAME(__setattr__);
  INIT_PY_NAME(__delattr__);
  INIT_PY_NAME(__del__);
  INIT_PY_NAME(__repr__);
  INIT_PY_NAME(__str__);
  INIT_PY_NAME(__class__);
  INIT_PY_NAME(__name__);
  INIT_PY_NAME(__hash__);
  INIT_PY_NAME(__cmp__);
  INIT_PY_NAME(__var_size__);
  INIT_PY_NAME(__init__);
  INIT_PY_NAME(__getinitargs__);
  INIT_PY_NAME(__getstate__);
  INIT_PY_NAME(__setstate__);
  INIT_PY_NAME(__dict__);
  INIT_PY_NAME(class_);
  INIT_PY_NAME(__module__);
  
#undef INIT_PY_NAME
}

static PyObject *
CallMethodO(PyObject *self, PyObject *name,
		     PyObject *args, PyObject *kw)
{
  if (! args && PyErr_Occurred()) return NULL;
  UNLESS(name=PyObject_GetAttr(self,name)) return NULL;
  ASSIGN(name,PyEval_CallObjectWithKeywords(name,args,kw));
  if (args) { Py_DECREF(args); }
  return name;
}

#define Build Py_BuildValue

/* CMethod objects: */

typedef struct {
  PyObject_HEAD
  PyTypeObject *type;
  PyObject *self;
  char		*name;
  PyCFunction	meth;
  int		flags;
  char		*doc;
} CMethod;

staticforward PyTypeObject CMethodType;

#define CMethod_Check(O) ((O)->ob_type==&CMethodType)
#define UnboundCMethod_Check(O) \
  ((O)->ob_type==&CMethodType && ! ((CMethod*)(O))->self)
#define AsCMethod(O) ((CMethod*)(O))
#define CMETHOD(O) ((CMethod*)(O))


#define PMethod PyECMethodObject
#define PMethodType PyECMethodObjectType

staticforward PyTypeObject PMethodType;

#define PMethod_Check(O) ((O)->ob_type==&PMethodType)
#define UnboundPMethod_Check(O) \
  ((O)->ob_type==&PMethodType && ! ((PMethod*)(O))->self)

#define UnboundEMethod_Check(O) \
  (((O)->ob_type==&PMethodType ||(O)->ob_type==&CMethodType) \
   && ! ((PMethod*)(O))->self)

#define PMETHOD(O) ((PMethod*)(O))

static PyObject *
bindPMethod(PMethod *m, PyObject *inst);

static PyObject *
#ifdef HAVE_STDARG_PROTOTYPES
/* VARARGS 2 */
JimErr_Format(PyObject *ErrType, char *stringformat, char *format, ...)
#else
/* VARARGS */
JimErr_Format(va_alist) va_dcl
#endif
{
  va_list va;
  PyObject *args=0, *retval=0;
#ifdef HAVE_STDARG_PROTOTYPES
  va_start(va, format);
#else
  PyObject *ErrType;
  char *stringformat, *format;
  va_start(va);
  ErrType = va_arg(va, PyObject *);
  stringformat   = va_arg(va, char *);
  format   = va_arg(va, char *);
#endif
  
  if (format) args = Py_VaBuildValue(format, va);
  va_end(va);
  if (format && ! args) return NULL;
  if (stringformat && !(retval=PyString_FromString(stringformat))) return NULL;

  if (retval)
    {
      if (args)
	{
	  PyObject *v;
	  v=PyString_Format(retval, args);
	  Py_DECREF(retval);
	  Py_DECREF(args);
	  if (! v) return NULL;
	  retval=v;
	}
    }
  else
    if (args) retval=args;
    else
      {
	PyErr_SetObject(ErrType,Py_None);
	return NULL;
      }
  PyErr_SetObject(ErrType,retval);
  Py_DECREF(retval);
  return NULL;
}


static PyObject *
#ifdef HAVE_STDARG_PROTOTYPES
/* VARARGS 2 */
JimString_Build(char *out_format, char *build_format, ...)
#else
/* VARARGS */
JimString_Build(va_alist) va_dcl
#endif
{
  va_list va;
  PyObject *args, *retval, *fmt;
#ifdef HAVE_STDARG_PROTOTYPES
  va_start(va, build_format);
#else
  char *build_format;
  char *out_format;
  va_start(va);
  out_format = va_arg(va, char *);
  build_format   = va_arg(va, char *);
#endif

  if (build_format)
    args = Py_VaBuildValue(build_format, va);
  else
    args = PyTuple_New(0);
  
  va_end(va);
  
  if (! args)
    return NULL;

  if (! PyTuple_Check(args))
    {
      PyObject *a;
      
      a=PyTuple_New(1);
      if (! a)
        return NULL;
      
      if (PyTuple_SetItem(a,0,args) == -1)
        return NULL;
      
      args=a;
    }

  fmt = PyString_FromString(out_format);
  
  retval = PyString_Format(fmt, args);
  
  Py_DECREF(args);
  Py_DECREF(fmt);
  
  return retval;
}

static PyObject *
EC_NewObject(PyTypeObject *type, int size)
{
    PyObject *inst;
    int len;

    if (type->tp_itemsize) {
	inst = PyObject_NEW_VAR(PyObject, type, size);
	if (inst == NULL)
	    return NULL;
	((PyVarObject *)inst)->ob_size = size;
    } 
    else {
	assert(size == 0);
	inst = PyObject_NEW(PyObject, type);
	if (inst == NULL)
	    return NULL;
    }
    Py_INCREF(type);
    len = (type->tp_basicsize + type->tp_itemsize * size) - sizeof(PyObject);
    memset(((char *)inst) + sizeof(PyObject), 0, len);
    return inst;
}

static int
CMethod_issubclass(PyExtensionClass *sub, PyExtensionClass *type)
{
  int i,l;
  PyObject *t;

  if (sub==type) return 1;
  if (! sub->bases) return 0;
  l=PyTuple_Size(sub->bases);
  for (i=0; i < l; i++)
    {
      t=PyTuple_GET_ITEM(sub->bases, i);
      if (t==(PyObject*)type) return 1;
      if (ExtensionClass_Check(t)
	 && AsExtensionClass(t)->bases
	 && CMethod_issubclass(AsExtensionClass(t),type)
	 ) return 1;
    }
  return 0;
}

#define Subclass_Check(C1,C2) \
  CMethod_issubclass((PyExtensionClass *)(C1), (PyExtensionClass *)(C2))

#define SubclassInstance_Check(C1,C2) \
  CMethod_issubclass((PyExtensionClass *)((C1)->ob_type), \
		     (PyExtensionClass *)(C2))


static CMethod *freeCMethod=0;

static PyObject *
newCMethod(PyExtensionClass *type, PyObject *inst,
	   char *name, PyCFunction meth, int flags, char *doc)
{
  CMethod *self;

  ALLOC_FREE(CMethod);
  
  Py_INCREF(type);
  Py_XINCREF(inst);
  self->type=(PyTypeObject*)type;
  self->self=inst;
  self->name=name;
  self->meth=meth;
  self->flags=flags;
  self->doc=doc;
  return (PyObject*)self;
}

static CMethod *
bindCMethod(CMethod *m, PyObject *inst)
{
  CMethod *self;
  
  UNLESS(inst->ob_type==m->type
	 || (ExtensionInstance_Check(inst)
	     && SubclassInstance_Check(inst,m->type))
	 || ((m->flags & METH_CLASS_METHOD) && ExtensionClass_Check(inst))
	 )
    {
      Py_INCREF(m);
      return m;
    }

  ALLOC_FREE(CMethod);

  Py_INCREF(inst);
  Py_INCREF(m->type);
  self->type=m->type;
  self->self=inst;
  self->name=m->name;
  self->meth=m->meth;
  self->flags=m->flags;
  self->doc=m->doc;
  return self;
}

static void
CMethod_dealloc(CMethod *self)
{
#ifdef TRACE_DEALLOC
  fprintf(stderr,"Deallocating C method %s\n", self->name); 
#endif
  Py_XDECREF(self->type);
  Py_XDECREF(self->self);
  self->self=(PyObject*)freeCMethod;
  freeCMethod=self;
}

typedef PyObject *(*call_by_name_type)(PyObject*,PyObject*,PyObject*,
				       PyTypeObject*);
typedef PyObject *(*by_name_type)(PyObject*,PyObject*,
				  PyTypeObject*);
static PyObject *
call_cmethod(CMethod *self, PyObject *inst, PyObject *args, PyObject *kw)
{
  if (!(self->flags & METH_VARARGS))
    {
      int size = PyTuple_Size(args);
      if (size == 1)      args = PyTuple_GET_ITEM(args, 0);
      else if (size == 0) args = NULL;
    }
  if (self->flags & METH_KEYWORDS)
    {
      if (self->flags & METH_BY_NAME)
	return (*(call_by_name_type)self->meth)(inst, args, kw,
					       self->type);
      else
	return (*(PyCFunctionWithKeywords)self->meth)(inst, args, kw);
    }
  else if (self->flags & METH_BY_NAME)
    return (*(by_name_type)self->meth)(inst, args, self->type);
  else
    {
      if (kw != NULL && PyDict_Size(kw) != 0)
	{
	  PyErr_SetString(PyExc_TypeError,
			  "this function takes no keyword arguments");
	  return NULL;
	}
      return (*self->meth)(inst, args);
    }
}

static char *hook_mark="C method being called through a hook.";

static PyObject *
callCMethodWithHook(CMethod *self, PyObject *inst,
		    PyObject *args, PyObject *kw)
{
  PyObject *hook, *m;

  UNLESS(m=newCMethod(AsExtensionClass(self->type),
		      inst, self->name, self->meth,
		      self->flags, hook_mark)) return NULL;

  if ((hook=PyObject_GetAttr(inst,py__call_method__)))
    {
      if ((CMethod_Check(hook) && CMETHOD(hook)->meth==self->meth)
	 ||
	 (PMethod_Check(hook)
	  && CMethod_Check(PMETHOD(hook)->meth)
	  && CMETHOD(PMETHOD(hook)->meth)->meth==self->meth)
	 )
	{
	  /* Oops, we are already calling the hook! */
	  Py_DECREF(hook);
	  return PyEval_CallObjectWithKeywords(m,args,kw);
	}
      if (kw)
	ASSIGN(hook,PyObject_CallFunction(hook,"OOO",m,args,kw));
      else
	ASSIGN(hook,PyObject_CallFunction(hook,"OO",m,args));
    }
  else
    {
      PyErr_Clear();
      hook=PyEval_CallObjectWithKeywords(m,args,kw);
    }

  Py_DECREF(m);
  return hook;
}

static PyObject *
CMethod_call(CMethod *self, PyObject *args, PyObject *kw)
{
  int size;

  if (self->self)
    {
      if (HasMethodHook(self->self) &&
	 self->doc != hook_mark	/* This check prevents infinite recursion */
	 )
	return callCMethodWithHook(self,self->self,args,kw);
      return call_cmethod(self,self->self,args,kw);
    }

  if ((size=PyTuple_Size(args)) > 0)
    {
      PyObject *first=0;

      UNLESS(first=PyTuple_GET_ITEM(args, 0)) return NULL;
      if (
	  first->ob_type==self->type
	  ||
	  (ExtensionInstance_Check(first)
	   &&
	   CMethod_issubclass(ExtensionClassOf(first),
			      AsExtensionClass(self->type))
	  )
	  )
      {
	PyObject *rest=0;
	UNLESS (rest=PySequence_GetSlice(args,1,size)) return NULL;
	if (HasMethodHook(first) &&
	   self->doc != hook_mark /* This check prevents infinite recursion */
	   )
          ASSIGN(rest, callCMethodWithHook(self, first, rest, kw) );
        else
          ASSIGN(rest, call_cmethod(self, first, rest, kw) );
	return rest;
      }
    }

  return JimErr_Format(PyExc_TypeError,
		      "unbound C method must be called with %s 1st argument",
		      "s", self->type->tp_name);
}

static PyObject *
CMethod_getattro(CMethod *self, PyObject *oname)
{
  PyObject *r;
  
  if (PyString_Check(oname))
    {
      char *name;

      UNLESS(name=PyString_AsString(oname)) return NULL;
        	        
      if (name[0] != '_' && name[0] && name[1] != '_' &&
	  PyEval_GetRestricted())
	{
	  PyErr_SetString(PyExc_RuntimeError,
	       	  "function attributes not accessible in restricted mode");
	  return NULL;
	}

      if (strcmp(name,"__name__")==0 || strcmp(name,"func_name")==0 )
	return PyString_FromString(self->name);
      if (strcmp(name,"func_code")==0 ||
	 strcmp(name,"im_func")==0)
	{
	  Py_INCREF(self);
	  return (PyObject *)self;
	}
      if (strcmp(name,"__doc__")==0 ||
	 strcmp(name,"func_doc")==0)
	{
	  if (self->doc)
	    return PyString_FromString(self->doc);
	  else
	    return PyString_FromString("");
	}
      if (strcmp(name,"im_class")==0)
	{
	  Py_INCREF(self->type);
	  return (PyObject *)self->type;
	}
      if (strcmp(name,"im_self")==0)
	{
	  if (self->self) r=self->self;
	  else           r=Py_None;
	  Py_INCREF(r);
	  return r;
	}
    }

  if (self->self)	/* Psuedo attributes */
    {
      UNLESS(oname=Py_BuildValue("sO", self->name, oname)) return NULL;
      UNLESS_ASSIGN(oname,PyString_Format(concat_fmt, oname)) return NULL;
      r=PyObject_GetAttr(OBJECT(self->self), py__class__);
      if (r)
	{
	  ASSIGN(r, PyObject_GetAttr(r, oname));

	  if (r) {
	    if (UnboundCMethod_Check(r))
	      ASSIGN(r, (PyObject*)bindCMethod((CMethod*)r, self->self));
	    else if (UnboundPMethod_Check(r))
	      ASSIGN(r, bindPMethod((PMethod*)r, self->self));
	  }
	}
      Py_DECREF(oname);
      return r;
    }

  PyErr_SetObject(PyExc_AttributeError, oname);
  return NULL;
}

static PyTypeObject CMethodType = {
  PyObject_HEAD_INIT(NULL)
  0,				/*ob_size*/
  "CMethod",			/*tp_name*/
  sizeof(CMethod),		/*tp_basicsize*/
  0,				/*tp_itemsize*/
  /* methods */
  (destructor)CMethod_dealloc,	/*tp_dealloc*/
  (printfunc)0,			/*tp_print*/
  0,				/*tp_getattr*/
  (setattrfunc)0,		/*tp_setattr*/
  (cmpfunc)0,			/*tp_compare*/
  (reprfunc)0,			/*tp_repr*/
  0,				/*tp_as_number*/
  0,				/*tp_as_sequence*/
  0,				/*tp_as_mapping*/
  (hashfunc)0,			/*tp_hash*/
  (ternaryfunc)CMethod_call,	/*tp_call*/
  (reprfunc)0,			/*tp_str*/
  (getattrofunc)CMethod_getattro, 	/* tp_getattro */
  (setattrofunc)0, 		/* tp_setattro */
  
  /* Space for future expansion */
  0L,0L,
  "Storage manager for unbound C function PyObject data"
  /* Documentation string */
};

/* PMethod objects: */

static PMethod *freePMethod=0;

static PyObject *
newPMethod(PyExtensionClass *type, PyObject *meth)
{
  PMethod *self;
  
  ALLOC_FREE(PMethod);

  Py_INCREF(type);
  Py_INCREF(meth);
  self->type=(PyTypeObject*)type;
  self->self=NULL;
  self->meth=meth;
  return (PyObject*)self;
}

static PyObject *
bindPMethod(PMethod *m, PyObject *inst)
{
  PMethod *self;

  if (NeedsToBeBound(m->meth))
    return CallMethodO(m->meth, py__of__, Build("(O)", inst), NULL);
  if (m->ob_refcnt==1)
    {
      Py_INCREF(inst);
      ASSIGN(m->self, inst);
      Py_INCREF(m);
      return (PyObject*)m;
    }
  
  ALLOC_FREE(PMethod);

  Py_INCREF(inst);
  Py_INCREF(m->type);
  Py_INCREF(m->meth);
  self->type=m->type;
  self->self=inst;
  self->meth=m->meth;
  return (PyObject*)self;
}

static PyObject *
PMethod_New(PyObject *meth, PyObject *inst)
{
  if (PMethod_Check(meth)) return bindPMethod((PMethod*)meth,inst);
  UNLESS(ExtensionInstance_Check(inst))
    return JimErr_Format(PyExc_TypeError,
			"Attempt to use %s as method for %s, which is "
			"not an extension class instance.",
			"OO",meth,inst);
  if ((meth=newPMethod(ExtensionClassOf(inst), meth)))
    UNLESS_ASSIGN(((PMethod*)meth)->self,inst) return NULL;
  Py_INCREF(inst);
  return meth;
}

static void
PMethod_dealloc(PMethod *self)
{
#ifdef TRACE_DEALLOC
  fprintf(stderr,"Deallocating PM ... ");
#endif
  Py_XDECREF(self->type);
  Py_XDECREF(self->self);
  Py_XDECREF(self->meth);
  self->self=(PyObject*)freePMethod;
  freePMethod=self;
#ifdef TRACE_DEALLOC
  fprintf(stderr," Done Deallocating PM\n");
#endif
}  

static PyObject *
callMethodWithPossibleHook(PyObject *inst,
			   PyObject *meth, PyObject *args, PyObject *kw)
{
  if (HasMethodHook(inst))
    {
      PyObject *hook;
      if ((hook=PyObject_GetAttr(inst,py__call_method__)))
	{
	  if (PMethod_Check(hook) && ((PMethod*)hook)->meth==meth)
	    {
	      /* Oops, we are already calling the hook! */
	      Py_DECREF(hook);
	      return PyEval_CallObjectWithKeywords(meth,args,kw);
	    }
	  if (kw)
	    ASSIGN(hook,PyObject_CallFunction(hook,"OOO",meth,args,kw));
	  else
	    ASSIGN(hook,PyObject_CallFunction(hook,"OO",meth,args));
	  return hook;
	}
      PyErr_Clear();
    }
  return PyEval_CallObjectWithKeywords(meth,args,kw);
}

static PyObject *
call_PMethod(PMethod *self, PyObject *inst, PyObject *args, PyObject *kw)
{
  PyObject *a;

  if (CMethod_Check(self->meth)
     && CMETHOD(self->meth)->type->tp_basicsize == sizeof(PyPureMixinObject)
     && ! (CMETHOD(self->meth)->self)
     )
    {
      /* Special HACK^H^H^Hcase:
	 we are wrapping an abstract unbound CMethod */
      if (HasMethodHook(inst) &&
	 /* This check prevents infinite recursion: */
	 CMETHOD(self->meth)->doc != hook_mark 
	 )
	return callCMethodWithHook(CMETHOD(self->meth),inst,args,kw);
      return call_cmethod(CMETHOD(self->meth),inst,args,kw);
      
    }
  else
    {
      a=Py_BuildValue("(O)",inst);
      if (a) ASSIGN(a,PySequence_Concat(a,args));
      if (a) ASSIGN(a,callMethodWithPossibleHook(inst,self->meth,a,kw));
      return a;
    }
}

static PyObject *
PMethod_call(PMethod *self, PyObject *args, PyObject *kw)
{
  int size;

  if (self->self) return call_PMethod(self,self->self,args,kw);

  if ((size=PyTuple_Size(args)) > 0)
    {
      PyObject *first=0, *ftype=0;
      UNLESS(first=PyTuple_GET_ITEM(args, 0)) return NULL;
      if (! self->type ||
	 ((ftype=PyObject_GetAttr(first,py__class__)) &&
	  (ftype==(PyObject*)self->type ||
	   (ExtensionClass_Check(ftype) &&
	    CMethod_issubclass(AsExtensionClass(ftype),
			       AsExtensionClass(self->type))
	    )
	   )
	  )
	 )
	{
	  if (NeedsToBeBound(self->meth))
	    {
	      PyObject *r, *rest;
	      UNLESS(r=CallMethodO(self->meth,py__of__,Build("(O)", first),
				   NULL))
		return NULL;
	      UNLESS(rest=PySequence_GetSlice(args,1,size))
		{
		  Py_DECREF(r);
		  return NULL;
		}
	      ASSIGN(r,callMethodWithPossibleHook(first,r,rest,kw));
	      Py_DECREF(rest);
	      return r;
	    }
	  Py_DECREF(ftype);
	  return callMethodWithPossibleHook(first,self->meth,args,kw);
	}
      Py_XDECREF(ftype);
    }

  return JimErr_Format(PyExc_TypeError,
		      "unbound Python method must be called with %s"
		      " 1st argument",
		      "s", self->type->tp_name);
}

static PyObject *
PMethod_repr(PMethod *self)
{
    char *func_name, buf[8192];
    int n;

    if (PyFunction_Check(self->meth)) {
	func_name = PyString_AS_STRING(
		((PyFunctionObject*)self->meth)->func_name);
    }
    else {
	/* self->meth is some other kind of object */
	func_name = "(?)";
    }
    
    if (self->self) {
	PyObject *repr = PyObject_Repr(self->self);
	if (!repr)
	    return NULL;
	n = sprintf(buf,
		     "<bound method %.1000s.%.1000s of %.1000s>",
		     self->type->tp_name, func_name,
		     PyString_AS_STRING(repr));
	Py_DECREF(repr);
    }
    else {
	n = sprintf(buf, 
		     "<unbound method %.1000s.%.1000s>",
		     self->type->tp_name, func_name);
    }
    return PyString_FromStringAndSize(buf, n);
}

static PyObject *
PMethod_getattro(PMethod *self, PyObject *oname)
{
  PyObject *r;

  if (PyString_Check(oname))
    {
      char *name;

      UNLESS(name=PyString_AsString(oname)) return NULL;

      if (name[0]=='_' && name[1]=='_')
	{
	  if (strcmp(name+2,"name__")==0)
	    return PyObject_GetAttrString(self->meth,"__name__");
	  if (strcmp(name+2,"doc__")==0)
	    return PyObject_GetAttrString(self->meth,"__doc__");
	}        
      else if (PyEval_GetRestricted())
	{
	  PyErr_SetString(PyExc_RuntimeError,
	       "function attributes not accessible in restricted mode");
	  return NULL;
	}
      else if (name[0]=='f' && name[1]=='u' && name[2]=='n' && name[3]=='c'
	      && name[4]=='_')
	{
	  if (strcmp(name+5,"name")==0 )
	    return PyObject_GetAttrString(self->meth,"__name__");
	  if (strcmp(name+5,"doc")==0)
	    return PyObject_GetAttrString(self->meth,"__doc__");
	}

      if (*name++=='i' && *name++=='m' && *name++=='_')
	{
	  if (strcmp(name,"func")==0)
	    {
	      Py_INCREF(self->meth);
	      return self->meth;
	    }
	  if (strcmp(name,"class")==0)
	    {
	      Py_INCREF(self->type);
	      return (PyObject *)self->type;
	    }
	  if (strcmp(name,"self")==0)
	    {
	      if (self->self) r=self->self;
	      else           r=Py_None;
	      Py_INCREF(r);
	      return r;
	    }
	}
    }

  if (self->meth)
    {
      if ((r=PyObject_GetAttr(self->meth, oname))) return r;
      PyErr_Clear();

      if (self->self) /* Psuedo attrs */
	{
	  PyObject *myname;
	  
	  UNLESS(myname=PyObject_GetAttr(self->meth, py__name__)) return NULL;
	  oname=Py_BuildValue("OO", myname, oname);
	  Py_DECREF(myname);
	  UNLESS(oname) return NULL;
	  UNLESS_ASSIGN(oname,PyString_Format(concat_fmt, oname)) return NULL;
	  r=PyObject_GetAttr(OBJECT(self->self), py__class__);
	  if (r)
	    {
	      ASSIGN(r, PyObject_GetAttr(r, oname));
      
	      if (r) {
		if (UnboundCMethod_Check(r))
		  ASSIGN(r, (PyObject*)bindCMethod((CMethod*)r, self->self));
		else if (UnboundPMethod_Check(r))
		  ASSIGN(r, bindPMethod((PMethod*)r, self->self));
	      }
	    }
	  Py_DECREF(oname);
	  return r;
	}
    }

  PyErr_SetObject(PyExc_AttributeError, oname);
  return NULL;

  return PyObject_GetAttr(self->meth, oname);
}

static PyTypeObject PMethodType = {
  PyObject_HEAD_INIT(NULL)
  0,					/*ob_size*/
  "Python Method",			/*tp_name*/
  sizeof(PMethod),			/*tp_basicsize*/
  0,					/*tp_itemsize*/
  /* methods */
  (destructor)PMethod_dealloc,		/*tp_dealloc*/
  (printfunc)0,				/*tp_print*/
  0,					/*tp_getattr*/
  (setattrfunc)0,			/*tp_setattr*/
  (cmpfunc)0,				/*tp_compare*/
  (reprfunc)PMethod_repr,		/*tp_repr*/
  0,					/*tp_as_number*/
  0,					/*tp_as_sequence*/
  0,					/*tp_as_mapping*/
  (hashfunc)0,				/*tp_hash*/
  (ternaryfunc)PMethod_call,		/*tp_call*/
  (reprfunc)0,				/*tp_str*/
  (getattrofunc)PMethod_getattro,	/*tp_getattro*/
  (setattrofunc)0, 			/* tp_setattro */
  
  /* Space for future expansion */
  0L,0L,
  "Storage manager for unbound C function PyObject data"
  /* Documentation string */
};

static PyObject *CCL_getattr(PyExtensionClass*,PyObject*,int);

static int
CCL_hasattr(PyExtensionClass *self,PyObject *name)
{
  PyObject *r;

  r=CCL_getattr(self, name, 0);
  if (r)
    {
      Py_DECREF(r);
      return 1;
    }
  else
    PyErr_Clear();
  
  return 0;
}

/* Special Methods */

#define UNARY_OP(OP) \
static PyObject * \
OP ## _by_name(PyObject *self, PyObject *args, PyTypeObject *ob_type) { \
  UNLESS(PyArg_ParseTuple(args,"")) return NULL; \
  return ob_type->tp_ ## OP(self); \
} 

UNARY_OP(repr)
UNARY_OP(str)

static PyObject * 
hash_by_name(PyObject *self, PyObject *args, PyTypeObject *ob_type) { 
  long r; 
  UNLESS(PyArg_ParseTuple(args,"")) return NULL; 
  UNLESS(-1 != (r=ob_type->tp_hash(self))) return NULL; 
  return PyInt_FromLong(r); 
} 

static PyObject *
call_by_name(PyObject *self, PyObject *args, PyObject *kw,
	     PyTypeObject *ob_type)
{
  return ob_type->tp_call(self,args,kw);
}

static PyObject *
compare_by_name(PyObject *self, PyObject *args, PyTypeObject *ob_type)
{
  PyObject *other;

  UNLESS(PyArg_ParseTuple(args,"O", &other)) return NULL; 
  return PyInt_FromLong(ob_type->tp_compare(self,other)); 
} 

static PyObject *
getattr_by_name(PyObject *self, PyObject *args, PyTypeObject *ob_type)
{
  char *name;
  UNLESS(PyArg_ParseTuple(args,"s",&name)) return NULL;
  return ob_type->tp_getattr(self,name);
}

static PyObject *
setattr_by_name(PyObject *self, PyObject *args, PyTypeObject *ob_type)
{
  char *name;
  PyObject *v;
  UNLESS(PyArg_ParseTuple(args,"sO",&name,&v)) return NULL;
  UNLESS(-1 != ob_type->tp_setattr(self,name,v)) return NULL;
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject *
delsetattr_by_name(PyObject *self, PyObject *args, PyTypeObject *ob_type)
{
  char *name;
  UNLESS(PyArg_ParseTuple(args,"s",&name)) return NULL;
  UNLESS(-1 != ob_type->tp_setattr(self,name,NULL)) return NULL;
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject *
getattro_by_name(PyObject *self, PyObject *args, PyTypeObject *ob_type)
{
  PyObject *name;
  UNLESS(PyArg_ParseTuple(args,"O",&name)) return NULL;
  return ob_type->tp_getattro(self,name);
}

static PyObject *
setattro_by_name(PyObject *self, PyObject *args, PyTypeObject *ob_type)
{
  PyObject *name;
  PyObject *v;
  UNLESS(PyArg_ParseTuple(args,"OO",&name,&v)) return NULL;
  UNLESS(-1 != ob_type->tp_setattro(self,name,v)) return NULL;
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject *
delsetattro_by_name(PyObject *self, PyObject *args, PyTypeObject *ob_type)
{
  PyObject *name;
  UNLESS(PyArg_ParseTuple(args,"O",&name)) return NULL;
  UNLESS(-1 != ob_type->tp_setattro(self,name,NULL)) return NULL;
  Py_INCREF(Py_None);
  return Py_None;
}
  
static PyObject * 
length_by_name(PyObject *self, PyObject *args, PyTypeObject *ob_type)
{ 
  long r; 
  UNLESS(PyArg_ParseTuple(args,"")) return NULL; 
  if (ob_type->tp_as_sequence)
    {
      UNLESS(-1 != (r=ob_type->tp_as_sequence->sq_length(self)))
	return NULL;
    }
  else
    {
      UNLESS(-1 != (r=ob_type->tp_as_mapping->mp_length(self)))
	return NULL;
    }
  return PyInt_FromLong(r); 
} 
  
static PyObject * 
getitem_by_name(PyObject *self, PyObject *args, PyTypeObject *ob_type)
{ 
  PyObject *key;
  
  UNLESS(PyArg_ParseTuple(args,"O",&key)) return NULL; 
  if (ob_type->tp_as_mapping)
    return ob_type->tp_as_mapping->mp_subscript(self,key);
  else
    {
      int index;
      UNLESS(-1 != (index=PyInt_AsLong(key))) return NULL;
      return ob_type->tp_as_sequence->sq_item(self,index);
    }
} 

static PyCFunction item_by_name=(PyCFunction)getitem_by_name;
static PyCFunction subscript_by_name=(PyCFunction)getitem_by_name;
  
static PyObject *
setitem_by_name(PyObject *self, PyObject *args, PyTypeObject *ob_type)
{ 
  PyObject *key, *v;
  long r;
  
  UNLESS(PyArg_ParseTuple(args,"OO",&key,&v)) return NULL; 
  if (ob_type->tp_as_mapping)
    r=ob_type->tp_as_mapping->mp_ass_subscript(self,key,v);
  else
    {
      int index;
      UNLESS(-1 != (index=PyInt_AsLong(key))) return NULL;
      r=ob_type->tp_as_sequence->sq_ass_item(self,index,v);
    }
  if (r < 0) return NULL;
  Py_INCREF(Py_None);
  return Py_None;
}

static PyCFunction ass_item_by_name=(PyCFunction)setitem_by_name;
static PyCFunction ass_subscript_by_name=(PyCFunction)setitem_by_name;

static PyObject *
slice_by_name(PyObject *self, PyObject *args, PyTypeObject *ob_type)
{
  int i1,i2;

  UNLESS(PyArg_ParseTuple(args,"ii",&i1,&i2)) return NULL;
  return ob_type->tp_as_sequence->sq_slice(self,i1,i2);
}

static PyObject *
ass_slice_by_name(PyObject *self, PyObject *args, PyTypeObject *ob_type)
{
  int i1,i2;
  PyObject *v;
  long r;

  UNLESS(PyArg_ParseTuple(args,"iiO",&i1,&i2,&v)) return NULL;
  r=ob_type->tp_as_sequence->sq_ass_slice(self,i1,i2,v);
  if (r<0) return NULL;
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject *
concat_by_name(PyObject *self, PyObject *args, PyTypeObject *ob_type)
{
  PyObject *other;
  UNLESS(PyArg_ParseTuple(args,"O",&other)) return NULL;
  return ob_type->tp_as_sequence->sq_concat(self,other);
}

static PyObject *
repeat_by_name(PyObject *self, PyObject *args, PyTypeObject *ob_type)
{
  int r;
  UNLESS(PyArg_ParseTuple(args,"i",&r)) return NULL;
  return ob_type->tp_as_sequence->sq_repeat(self,r);
}

#define BINOP(OP,AOP) \
static PyObject * \
OP ## _by_name(PyObject *self, PyObject *args, PyTypeObject *ob_type) { \
  PyObject *v; \
  UNLESS(PyArg_ParseTuple(args,"O",&v)) return NULL; \
  return ob_type->tp_as_number->nb_ ## OP(self, v); \
}

BINOP(add,Add)
BINOP(subtract,Subtract)
BINOP(multiply,Multiply)
BINOP(divide,Divide)
BINOP(remainder,Remainder)
BINOP(divmod,Divmod)

static PyObject *
power_by_name(PyObject *self, PyObject *args, PyTypeObject *ob_type)
{
  PyObject *v, *z=NULL;
  UNLESS(PyArg_ParseTuple(args,"O|O",&v,&z)) return NULL; 
  return ob_type->tp_as_number->nb_power(self,v,z);
}

#define UNOP(OP) \
static PyObject * \
OP ## _by_name(PyObject *self, PyObject *args, PyTypeObject *ob_type) { \
  UNLESS(PyArg_ParseTuple(args,"")) return NULL; \
  return ob_type->tp_as_number->nb_ ## OP(self); \
}

UNOP(negative)
UNOP(positive)
UNOP(absolute)

static PyObject * 
nonzero_by_name(PyObject *self, PyObject *args, PyTypeObject *ob_type) { 
  long r; 
  UNLESS(PyArg_ParseTuple(args,"")) return NULL; 
  UNLESS(-1 != (r=ob_type->tp_as_number->nb_nonzero(self))) return NULL; 
  return PyInt_FromLong(r); 
} 

UNOP(invert)

BINOP(lshift,Lshift)
BINOP(rshift,Rshift)
BINOP(and,And)
BINOP(or,Or)
BINOP(xor,Xor)

static PyObject *
coerce_by_name(PyObject *self, PyObject *args, PyTypeObject *ob_type)
{
  PyObject *v;
  int r;
  UNLESS(PyArg_ParseTuple(args,"O", &v)) return NULL;
  UNLESS(-1 != (r=ob_type->tp_as_number->nb_coerce(&self,&v)))
    {
      Py_INCREF(Py_None);
      return Py_None;
    }
  args=Py_BuildValue("OO",self,v);
  Py_DECREF(self);
  Py_DECREF(v);
  return args;
} 

UNOP(long)
UNOP(int)
UNOP(float)
UNOP(oct)
UNOP(hex)

#define FILLENTRY(T,MN,N,F,D) if (T ## _ ## MN) { \
  UNLESS(-1 != PyMapping_SetItemString(dict,"__" # N "__", \
    newCMethod(type, NULL, "__" # N "__", \
               (PyCFunction)MN ## _by_name, F | METH_BY_NAME, # D))) \
    goto err; }

#define delFILLENTRY(T,MN,N,F,D) if (T ## _ ## MN) { \
  UNLESS(-1 != PyMapping_SetItemString(dict,"__" # N "__", \
    newCMethod(type, NULL, "__" # N "__", \
               (PyCFunction)del ## MN ## _by_name, F | METH_BY_NAME, # D))) \
    goto err; }

static PyObject *
getBaseDictionary(PyExtensionClass *type)
{
  PyNumberMethods *nm;
  PySequenceMethods *sm;
  PyMappingMethods *mm;
  PyObject *dict;

  UNLESS(dict=type->class_dictionary)
    UNLESS(dict=PyDict_New()) return NULL;
  
  FILLENTRY(type->tp, repr, repr, METH_VARARGS,
	    "convert to an expression string");
  FILLENTRY(type->tp, hash, hash, METH_VARARGS, "compute a hash value");
  FILLENTRY(type->tp, call, call, METH_VARARGS | METH_KEYWORDS,
	    "call as a function");
  FILLENTRY(type->tp, compare, comp, METH_VARARGS,
	    "compare with another object");
 
  UNLESS (type->class_flags & EXTENSIONCLASS_PYTHONICATTR_FLAG)
    {
      FILLENTRY(type->tp, getattr, getattr, METH_VARARGS, "Get an attribute");
      FILLENTRY(type->tp, setattr, setattr, METH_VARARGS, "Set an attribute");
      delFILLENTRY(type->tp, setattr, delattr, METH_VARARGS, 
		   "Delete an attribute");
      
      FILLENTRY(type->tp, getattro, getattr, METH_VARARGS, "Get an attribute");
      FILLENTRY(type->tp, setattro, setattr, METH_VARARGS, "Set an attribute");
      delFILLENTRY(type->tp, setattro, delattr, METH_VARARGS, 
		   "Delete an attribute");
    }

  if ((sm=type->tp_as_sequence))
    {
      FILLENTRY(sm->sq, length, len, METH_VARARGS, "Get the object length");
      FILLENTRY(sm->sq, repeat, mul, METH_VARARGS,
		"Get a new object that is the object repeated.");
      FILLENTRY(sm->sq, item, getitem, METH_VARARGS, "Get an item");
      FILLENTRY(sm->sq, slice, getslice, METH_VARARGS, "Get a slice");
      FILLENTRY(sm->sq, ass_item, setitem, METH_VARARGS, "Assign an item");
      FILLENTRY(sm->sq, ass_slice, setslice, METH_VARARGS, "Assign a slice");
    }      

  if ((mm=type->tp_as_mapping))
    {
      FILLENTRY(mm->mp, length, len, METH_VARARGS, "Get the object length");
      FILLENTRY(mm->mp, subscript, getitem, METH_VARARGS, "Get an item");
      FILLENTRY(mm->mp, ass_subscript, setitem, METH_VARARGS,
		"Assign an item");
    }      

  if ((nm=type->tp_as_number) != NULL)
    {
      FILLENTRY(nm->nb, add, add, METH_VARARGS, "Add to another");
      FILLENTRY(nm->nb, subtract, sub, METH_VARARGS, "Subtract another");
      FILLENTRY(nm->nb, multiply, mul, METH_VARARGS, "Multiple by another");
      FILLENTRY(nm->nb, divide, div, METH_VARARGS, "Divide by another");
      FILLENTRY(nm->nb, remainder, mod, METH_VARARGS, "Compute a remainder");
      FILLENTRY(nm->nb, power, pow, METH_VARARGS, "Raise to a power");
      FILLENTRY(nm->nb, divmod, divmod, METH_VARARGS,
		"Compute the whole result and remainder of dividing\n"
		"by another");
      FILLENTRY(nm->nb, negative, neg, METH_VARARGS,
		"Get the negative value.");
      FILLENTRY(nm->nb, positive, pos, METH_VARARGS, "Compute positive value");
      FILLENTRY(nm->nb, absolute, abs, METH_VARARGS, "Compute absolute value");
      FILLENTRY(nm->nb, nonzero, nonzero, METH_VARARGS,
		"Determine whether nonzero");
      FILLENTRY(nm->nb, invert, inv, METH_VARARGS, "Compute inverse");
      FILLENTRY(nm->nb, lshift, lshift, METH_VARARGS, "Shist left");
      FILLENTRY(nm->nb, rshift, rshift, METH_VARARGS, "Shist right");
      FILLENTRY(nm->nb, and, and, METH_VARARGS, "bitwize logical and");
      FILLENTRY(nm->nb, or, or, METH_VARARGS, "bitwize logical or");
      FILLENTRY(nm->nb, xor, xor, METH_VARARGS, "bitwize logical excusive or");
      FILLENTRY(nm->nb, coerce, coerce, METH_VARARGS,
		"Coerce woth another to a common type");
      FILLENTRY(nm->nb, int, int, METH_VARARGS, "Convert to an integer");
      FILLENTRY(nm->nb, long, long, METH_VARARGS,
		"Convert to an infinite-precision integer");
      FILLENTRY(nm->nb, float, float, METH_VARARGS,
		"Convert to floating point number");
      FILLENTRY(nm->nb, oct, oct, METH_VARARGS, "Convert to an octal string");
      FILLENTRY(nm->nb, hex, hex, METH_VARARGS,
		"Convert to a hexadecimal string");
    }

  if ((sm=type->tp_as_sequence))
    {
      FILLENTRY(sm->sq, concat, add, METH_VARARGS,
		"Concatinate the object with another");
    }      

  return dict;
err:
  Py_DECREF(dict);
  return NULL;
}

#undef UNARY_OP
#undef BINOP
#undef UNOP
#undef FILLENTRY

PyObject *
EC_reduce(PyObject *self, PyObject *args)
{
  PyObject *state=0;

  if ((args=PyObject_GetAttr(self,py__getinitargs__)))
    {
      UNLESS_ASSIGN(args,PyEval_CallObject(args,NULL)) return NULL;
      UNLESS_ASSIGN(args,PySequence_Tuple(args)) return NULL;
    }
  else
    {
      PyErr_Clear();
      if (ExtensionClassOf(self)->class_flags & EXTENSIONCLASS_BASICNEW_FLAG)
	{
	  args=Py_None;
	  Py_INCREF(args);
	}
      else args=PyTuple_New(0);
    }

  if ((state=PyObject_GetAttr(self,py__getstate__)))
    {
      UNLESS_ASSIGN(state,PyEval_CallObject(state,NULL)) goto err;
      ASSIGN(args,Py_BuildValue("OOO", self->ob_type, args, state));
      Py_DECREF(state);
    }
  else
    {
      PyErr_Clear();

      if ((state=PyObject_GetAttr(self, py__dict__)))
	{
	  ASSIGN(args,Py_BuildValue("OOO", self->ob_type, args, state));
	  Py_DECREF(state);
	}
      else
	{
	  PyErr_Clear();
	  ASSIGN(args, Py_BuildValue("OO", self->ob_type, args));
	}
    }

  return args;

err:
  Py_DECREF(args);
  return NULL;
}

static PyObject *
inheritedAttribute(PyExtensionClass *self, PyObject *args)
{
  PyObject *name;

  UNLESS(PyArg_ParseTuple(args,"O!",&PyString_Type, &name)) return NULL;

  return CCL_getattr(AsExtensionClass(self),name,1);
}

static PyObject *
basicnew(PyExtensionClass *self, PyObject *args)
{
  PyObject *inst=0;
  int size = 0;

  if (! self->tp_dealloc)
    {
      PyErr_SetString(PyExc_TypeError,
		      "Attempt to create instance of an abstract type");
      return NULL;
    }

  UNLESS(self->class_flags & EXTENSIONCLASS_BASICNEW_FLAG)
    return PyObject_CallObject(OBJECT(self), NULL);      

  if (self->tp_itemsize)
    {
      /* We have a variable-sized object, we need to get it's size */
      PyObject *var_size;
      
      UNLESS(var_size=CCL_getattr(self, py__var_size__, 0)) return NULL;
      UNLESS_ASSIGN(var_size,PyObject_CallObject(var_size,NULL)) return NULL;
      size=PyInt_AsLong(var_size);
      if (PyErr_Occurred()) return NULL;
    }
  UNLESS(inst=EC_NewObject((PyTypeObject *)self, size))
    return NULL;

  if (ClassHasInstDict(self))
    UNLESS(INSTANCE_DICT(inst)=PyDict_New()) goto err;

  if (self->bases && subclass_watcher &&
     ! PyObject_CallMethod(subclass_watcher,"created","O",inst))
    PyErr_Clear();

  return inst;

err:
  Py_DECREF(inst);
  return NULL;
}

struct PyMethodDef ECI_methods[] = {
  {"__reduce__",(PyCFunction)EC_reduce, METH_VARARGS,
   "__reduce__() -- Reduce an instance into it's class and creation data"
  },
  {"inheritedAttribute",(PyCFunction)inheritedAttribute,
   METH_VARARGS | METH_CLASS_METHOD,
   "inheritedAttribute(class,name) -- Get an inherited attribute\n\n"
   "Get an attribute that would be inherited if the given (extension)\n"
   "class did not define it.  This method is used when overriding\n"
   "inherited methods.  It provides 2 advantages over accessing\n"
   "\n"
   "attributes directly through a superclass:\n"
   "\n"
   "1. The superclass need not be known,\n"
   "\n"
   "2. The superclass may be a Python class.  Without this method, it would\n"
   "   be possible to override methods inherited from python classes because\n"
   "   unbound methods gotten from Python classes cannot be called with \n"
   "   extension class instances.  \n"
  },
  {"__basicnew__",(PyCFunction)basicnew,
   METH_VARARGS | METH_CLASS_METHOD,
   "__basicnew__() -- return a new uninitialized instance"
  },
  {NULL,		NULL}		/* sentinel */
};

static PyObject *
initializeBaseExtensionClass(PyExtensionClass *self)
{
  static PyMethodChain top = { ECI_methods, NULL };
  PyMethodChain *chain;
  PyObject *dict;
  int abstract;

  /* Is this an abstract, or at least a dataless, class? */
  abstract=self->tp_basicsize == sizeof(PyPureMixinObject);

  self->ob_type=(PyTypeObject*)&ECType;
  Py_INCREF(self->ob_type);

  UNLESS(dict=self->class_dictionary=getBaseDictionary(self)) return NULL;

  if (self->tp_name)
    {
      PyObject *name;

      UNLESS(name=PyString_FromString(self->tp_name)) goto err;
      if (0 > PyMapping_SetItemString(dict,"__doc__",name)) goto err;
      Py_DECREF(name);
    }
  else if (0 > PyMapping_SetItemString(dict,"__doc__",Py_None)) goto err;
  
  if (&self->methods) chain=&(self->methods);
  else chain=&top;
  
  while (1)
    {
      PyMethodDef *ml = chain->methods;

      for (; ml && ml->ml_name != NULL; ml++) 
	{
	  if (ml->ml_meth)
	    {
	      if (! PyMapping_HasKeyString(dict,ml->ml_name))
		{
		  PyObject *m;

		  /* Note that we create a circular reference here.
		     I suppose that this isn't so bad, since this is
		     probably a static thing anyway. Still, it is a
		     bit troubling. Oh well.
		  */
		  if (ml->ml_flags & METH_CLASS_METHOD)
		    {
		      UNLESS(m=newCMethod(
                         AsExtensionClass(self->ob_type), NULL,
			 ml->ml_name, ml->ml_meth,
			 ml->ml_flags, ml->ml_doc))
			return NULL;
		    }
		  else
		    {
		      UNLESS(m=newCMethod(self, NULL, ml->ml_name, ml->ml_meth,
					  ml->ml_flags, ml->ml_doc))
			return NULL;
		  
		      if (abstract)
			UNLESS_ASSIGN(m, newPMethod(self, m))
			  return NULL;
		    }

		  if (PyMapping_SetItemString(dict,ml->ml_name,m) < 0)
		    return NULL;
		}
	    }
	  else if (ml->ml_doc && *(ml->ml_doc))
	    {
	      /* No actual meth, this is probably to hook a doc string
		 onto a special method. */
	      PyObject *m;

	      if ((m=PyMapping_GetItemString(dict,ml->ml_name)))
		{
		  if (m->ob_type==&CMethodType)
		    ((CMethod *)(m))->doc=ml->ml_doc;
		}
	      else
		PyErr_Clear();
	    }
	}
      
      if (chain == &top) break;

      UNLESS(chain=chain->link) chain=&top;
    }
  return (PyObject*)self;

err:
  Py_DECREF(dict);
  return NULL;
}

static void
CCL_dealloc(PyExtensionClass *self)
{
#ifdef TRACE_DEALLOC
  fprintf(stderr,"Deallocating %s\n", self->tp_name);
#endif
  Py_XDECREF(self->class_dictionary);
  if (self->bases)
    {
      /* If we are a subclass, then we strduped our name */
      free(self->tp_name);

      /* And we allocated our own protocol structures */
      if (self->tp_as_number)   free(self->tp_as_number);
      if (self->tp_as_sequence) free(self->tp_as_sequence);
      if (self->tp_as_mapping)  free(self->tp_as_mapping);
      
      Py_DECREF(self->bases);
    }
  if (((PyExtensionClass*)self->ob_type) != self) {
      Py_XDECREF(self->ob_type);
  }
  PyObject_DEL(self);
}
  
static PyObject *
ExtensionClass_FindInstanceAttribute(PyObject *inst, PyObject *oname,
				     char *name)
{
  /* Look up an attribute for an instance from:

     The instance dictionary,
     The class dictionary, or
     The base objects.
   */
  PyObject *r=0;
  PyExtensionClass *self;

  if (! name) return NULL;

  self=(PyExtensionClass*)(inst->ob_type);

  if (*name=='_' && name[1]=='_')
    {
      char *n=name+2;
      if (*n == 'c' && strcmp(n,"class__")==0)
	{
	  Py_INCREF(self);
	  return (PyObject*)self;
	}
      if (ClassHasInstDict(self) && *n=='d' && strcmp(n,"dict__")==0)
	{
	  r = INSTANCE_DICT(inst);
	  Py_INCREF(r);
	  return r;
	}
    }

  if (ClassHasInstDict(self))
    {
      r= INSTANCE_DICT(inst);
      if (PyDict_Check(r))
        {
          r = PyDict_GetItem(r,oname);
          Py_XINCREF(r);
        }
      else
        {
          UNLESS (r = PyObject_GetItem(r,oname))
            PyErr_Clear();
        }

      if (r)
        {
          if (NeedsToBeBound(r))
            {
              ASSIGN(r, CallMethodO(r, py__of__, Build("(O)", inst), NULL));
            }
          return r;
        }
    }

  if (*name=='_' && name[1]=='_' 
      && 
      (   (name[2]=='b' && strcmp(name+2,"bases__")==0) 
          || (name[2]=='d' && strcmp(name+2,"dict__")==0)
          )
      )
    {
      PyErr_SetObject(PyExc_AttributeError, oname);
      return NULL;
    }
  
  UNLESS(r=CCL_getattr(self,oname,0)) return NULL;
  
  /* We got something from our class, maybe its an unbound method. */
  if (UnboundCMethod_Check(r))
    ASSIGN(r,(PyObject*)bindCMethod((CMethod*)r,inst));
  else if (UnboundPMethod_Check(r))
    ASSIGN(r,bindPMethod((PMethod*)r,inst));
      
  return r;
}

static PyObject *
EC_findiattrs(PyObject *self, char *name)
{
  PyObject *s, *r;

  UNLESS(s=PyString_FromString(name)) return NULL;
  r=ExtensionClass_FindInstanceAttribute(self,s,name);
  Py_DECREF(s);
  return r;
}
  
static PyObject *
EC_findiattro(PyObject *self, PyObject *name)
{
  return ExtensionClass_FindInstanceAttribute(self,name,
					      PyString_AsString(name));
}

static int
subclass_simple_setattr(PyObject *self, char *name, PyObject *v);

static PyObject *
CCL_getattr2(PyObject *self, PyObject *oname, int look_super)
{
  PyObject *r=0, *b, *d;

  if (ExtensionClass_Check(self))
    {
      b=((PyExtensionClass*)self)->bases;
      d=((PyExtensionClass*)self)->class_dictionary;
    }
  else if (PyClass_Check(self))
    {
      b=((PyClassObject*)self)->cl_bases;
      d=((PyClassObject*)self)->cl_dict;
    }
  else
    {
      UNLESS (r=PyObject_GetAttr(self, oname)) PyErr_Clear();
      return r;
    }

  if (! look_super && d)
    {
      if (PyDict_Check(d))
        {
          if((r=PyDict_GetItem(d, oname)))
            {          
              Py_INCREF(r);
              return r;
            }
        }
      else
        {
          if((r=PyObject_GetItem(d, oname))) return r;
          PyErr_Clear();
        }
    }

  if (b)
    {
      int n, i;
      
      n = PyTuple_Check(b) ? PyTuple_GET_SIZE(b) : 0; /* I don't care ;) */
      for (i=0; i < n; i++)
        {
          r=CCL_getattr2(PyTuple_GET_ITEM(b, i), oname, 0);
          if (r) return r;
        }
    }

  return NULL;
}

static PyObject *
CCL_getattr(PyExtensionClass *self, PyObject *oname, int look_super)
{
  PyObject *r=0;

  UNLESS (r=CCL_getattr2(OBJECT(self), oname, look_super)) 
    {
      PyErr_SetObject(PyExc_AttributeError, oname);
      return NULL;
    }

  if (PyFunction_Check(r) || NeedsToBeBound(r))
    ASSIGN(r,newPMethod(self,r));
  else if (PyMethod_Check(r) && ! PyMethod_Self(r))
    ASSIGN(r,newPMethod(self, PyMethod_Function(r)));

  return r;
}

static PyObject *
CCL_reduce(PyExtensionClass *self, PyObject *args)
{
  return PyString_FromString(self->tp_name);
}

PyObject *
CCL_getattro(PyExtensionClass *self, PyObject *name)
{
  char *n, *nm=0;
  PyObject *r;

  if (PyString_Check(name) && (n=nm=PyString_AS_STRING((PyStringObject*)name)))
    {
      if (*n=='_' && *++n=='_')
	{
	  switch (*++n)
	    {
	    case 's':
	      if (strcmp(n,"safe_for_unpickling__")==0)
		return PyInt_FromLong(1);
	      break;
	    case 'n':
	      if (strcmp(n,"name__")==0)
		return PyString_FromString(self->tp_name);
	      break;
	    case 'b':
	      if (strcmp(n,"bases__")==0)
		{
		  if (self->bases)
		    {
		      Py_INCREF(self->bases);
		      return self->bases;
		    }
		  else
		    return PyTuple_New(0);
		}
	      break;
	    case 'r':
	      if (strcmp(n,"reduce__")==0)
		return newCMethod(self,(PyObject*)self,
		   "__reduce__",(PyCFunction)CCL_reduce,0,
		   "__reduce__() -- Reduce the class to a class name");
	      break;
	    case 'd':
	      if (strcmp(n,"dict__")==0)
		{
		  Py_INCREF(self->class_dictionary);
		  return self->class_dictionary;
		}
	      break;
	    case 'c':
	      if (strcmp(n,"class__")==0)
		{
		  Py_INCREF(self->ob_type);
		  return OBJECT(self->ob_type);
		}
	      break;
	    }
	}
    }

  if ((r=CCL_getattr(self,name,0)))
    {
      if (UnboundCMethod_Check(r) && (AsCMethod(r)->flags & METH_CLASS_METHOD))
	ASSIGN(r,(PyObject*)bindCMethod((CMethod*)r,OBJECT(self)));
    }
   
  return r;
}

static int
CCL_setattro(PyExtensionClass *self, PyObject *name, PyObject *v)
{
  if (! v) return PyObject_DelItem(self->class_dictionary, name);

  if (v && UnboundCMethod_Check(v) &&
     ! (self->class_flags & EXTENSIONCLASS_METHODHOOK_FLAG)
     )
    {
      char *n;
      PyNumberMethods *nm;
      PySequenceMethods *s, *ms;
      PyMappingMethods *m, *mm;

      UNLESS(n=PyString_AsString(name)) return -1;
      if (*n++=='_' && *n++=='_')
	{
#define SET_SPECIAL(C,P) \
	  if (strcmp(n,#P "__")==0 \
	     && AsCMethod(v)->meth==(PyCFunction)C ## _by_name \
	     && Subclass_Check(self,AsCMethod(v)->type)) { \
	      self->tp_ ## C=AsCMethod(v)->type->tp_ ## C; \
	      return PyObject_SetItem(self->class_dictionary, name, v); }
	  /*
	  SET_SPECIAL(setattr,setattr);
	  SET_SPECIAL(setattro,setattr);
	  */
	  SET_SPECIAL(compare,cmp);
	  SET_SPECIAL(hash,hash);
	  SET_SPECIAL(repr,repr);
	  SET_SPECIAL(call,call);
	  SET_SPECIAL(str,str);
#undef SET_SPECIAL

#define SET_SPECIAL(C,P) \
	  if (strcmp(n,#P "__")==0 \
	     && AsCMethod(v)->meth==(PyCFunction)C ## _by_name \
	     && Subclass_Check(self,AsCMethod(v)->type) \
	     && (nm=self->tp_as_number)) { \
	      nm->nb_ ## C=AsCMethod(v)->type->tp_as_number->nb_ ## C; \
	      return PyObject_SetItem(self->class_dictionary, name, v); } 

	  SET_SPECIAL(add,add);
	  SET_SPECIAL(subtract,sub);
	  SET_SPECIAL(multiply,mult);
	  SET_SPECIAL(divide,div);
	  SET_SPECIAL(remainder,mod);
	  SET_SPECIAL(power,pow);
	  SET_SPECIAL(divmod,divmod);
	  SET_SPECIAL(lshift,lshift);
	  SET_SPECIAL(rshift,rshift);
	  SET_SPECIAL(and,and);
	  SET_SPECIAL(or,or);
	  SET_SPECIAL(xor,xor);
	  SET_SPECIAL(coerce,coerce);
	  SET_SPECIAL(negative,neg);
	  SET_SPECIAL(positive,pos);
	  SET_SPECIAL(absolute,abs);
	  SET_SPECIAL(nonzero,nonzero);
	  SET_SPECIAL(invert,inv);
	  SET_SPECIAL(int,int);
	  SET_SPECIAL(long,long);
	  SET_SPECIAL(float,float);
	  SET_SPECIAL(oct,oct);
	  SET_SPECIAL(hex,hex);
#undef SET_SPECIAL

	  if (strcmp(n,"len__")==0
	     && AsCMethod(v)->meth==(PyCFunction)length_by_name 
	     && Subclass_Check(self,AsCMethod(v)->type))
	     {
	       if ((s=self->tp_as_sequence) &&
		  (ms=AsCMethod(v)->type->tp_as_sequence) &&
		  ms->sq_length)
		 s->sq_length=ms->sq_length;
	       if ((m=self->tp_as_mapping) &&
		  (mm=AsCMethod(v)->type->tp_as_mapping) &&
		  mm->mp_length)
		 m->mp_length=mm->mp_length;
	       return PyObject_SetItem(self->class_dictionary, name, v);
	     } 

	  if (strcmp(n,"getitem__")==0
	     && AsCMethod(v)->meth==(PyCFunction)getitem_by_name 
	     && Subclass_Check(self,AsCMethod(v)->type))
	     {
	       if ((s=self->tp_as_sequence) &&
		  (ms=AsCMethod(v)->type->tp_as_sequence) &&
		  ms->sq_item)
		 s->sq_item=ms->sq_item;
	       if ((m=self->tp_as_mapping) &&
		  (mm=AsCMethod(v)->type->tp_as_mapping) &&
		  mm->mp_subscript)
		 m->mp_subscript=mm->mp_subscript;
	       return PyObject_SetItem(self->class_dictionary, name, v);
	     } 

	  if (strcmp(n,"setitem__")==0 &&
	     AsCMethod(v)->meth==(PyCFunction)setitem_by_name 
	     && Subclass_Check(self,AsCMethod(v)->type))
	     {
	       if ((s=self->tp_as_sequence) &&
		  (ms=AsCMethod(v)->type->tp_as_sequence) &&
		  ms->sq_ass_item)
		 s->sq_ass_item=ms->sq_ass_item;
	       if ((m=self->tp_as_mapping) &&
		  (mm=AsCMethod(v)->type->tp_as_mapping) &&
		  mm->mp_ass_subscript)
		 m->mp_ass_subscript=mm->mp_ass_subscript;
	       return PyObject_SetItem(self->class_dictionary, name, v);
	     } 

#define SET_SPECIAL(C,P) \
	  if (strcmp(n,#P "__")==0 \
	     && AsCMethod(v)->meth==(PyCFunction)C ## _by_name \
	     && Subclass_Check(self,AsCMethod(v)->type) \
	     && (s=self->tp_as_sequence)) { \
	      s->sq_ ## C=AsCMethod(v)->type->tp_as_sequence->sq_ ## C; \
	      return PyObject_SetItem(self->class_dictionary, name, v); } 
	  SET_SPECIAL(slice,getslice);
	  SET_SPECIAL(ass_slice,setslice);
	  SET_SPECIAL(concat,concat);
	  SET_SPECIAL(repeat,repeat);
#undef SET_SPECIAL

	}
    }
  return PyObject_SetItem(self->class_dictionary, name, v);
}

static PyObject *
CCL_call(PyExtensionClass *self, PyObject *arg, PyObject *kw)
{
  PyObject *inst=0, *init=0, *args=0;
  int size = 0;

  if (! self->tp_dealloc)
    {
      PyErr_SetString(PyExc_TypeError,
		      "Attempt to create instance of an abstract type");
      return NULL;
    }

  if (self->tp_itemsize)
    {
      /* We have a variable-sized object, we need to get it's size */
      PyObject *var_size;
      
      if ((var_size=CCL_getattr(self,py__var_size__, 0)))
	{
	  UNLESS_ASSIGN(var_size,PyObject_CallObject(var_size,arg))
	    return NULL;
	  size=PyInt_AsLong(var_size);
	  if (PyErr_Occurred()) return NULL;
	}
      else
	{
	  UNLESS(-1 != (size=PyTuple_Size(arg))) return NULL;
	  if (size > 0)
	    {
	      var_size=PyTuple_GET_ITEM(arg, 0);
	      if (PyInt_Check(var_size))
		size=PyInt_AsLong(var_size);
	      else
		size=-1;
	    }
	  else
	    size=-1;
	  if (size < 0)
	    {
	      PyErr_SetString(PyExc_TypeError,
			      "object size expected as first argument");
	      return NULL;
	    }
	}
    }
  UNLESS(inst=EC_NewObject((PyTypeObject *)self, size)) return NULL;

  if (ClassHasInstDict(self))
    UNLESS(INSTANCE_DICT(inst)=PyDict_New()) goto err;

   if ((init=CCL_getattr(self,py__init__,0)))
      {
       UNLESS(args=Py_BuildValue("(O)",inst)) goto err;
       if (arg) UNLESS_ASSIGN(args,PySequence_Concat(args,arg)) goto err;
       UNLESS_ASSIGN(args,PyEval_CallObjectWithKeywords(init,args,kw)) goto err;
       Py_DECREF(args);
       Py_DECREF(init);
      }
   else PyErr_Clear();

  if (self->bases && subclass_watcher &&
     ! PyObject_CallMethod(subclass_watcher,"created","O",inst))
    PyErr_Clear();

  return inst;
err:
  Py_DECREF(inst);
  Py_XDECREF(init);
  Py_XDECREF(args);
  return NULL;
}

static PyObject *
CCL_repr(PyExtensionClass *self)
{
  char p[128], *pp;
  PyObject *m;

  if ((m=PyObject_GetAttr(OBJECT(self), py__module__)))
    {
      if (! PyObject_IsTrue(m)) 
	{
	  Py_DECREF(m);
	  m=0;
	}
    }
  else PyErr_Clear();

  sprintf(p,"%p",self);
  if (*p=='0' && p[1]=='x') pp=p+2;
  else                      pp=p;
			      

  if (m) ASSIGN(m, JimString_Build("<extension class %s.%s at %s>","Oss",
				   m, self->tp_name, pp));
  else          m= JimString_Build("<extension class %s at %s>","ss",
				   self->tp_name, pp);

  return m;
}

static PyTypeObject ECTypeType = {
  PyObject_HEAD_INIT(NULL)
  0,				/*ob_size*/
  "ExtensionClass Class",	/*tp_name*/
  sizeof(PyExtensionClass),    	/*tp_basicsize*/
  0,				/*tp_itemsize*/
  /* methods */
  (destructor)CCL_dealloc,	/*tp_dealloc*/
  (printfunc)0,			/*tp_print*/
  (getattrfunc)0,		/*tp_getattr*/
  (setattrfunc)0,		/*tp_setattr*/
  (cmpfunc)0,			/*tp_compare*/
  (reprfunc)CCL_repr,		/*tp_repr*/
  0,				/*tp_as_number*/
  0,				/*tp_as_sequence*/
  0,				/*tp_as_mapping*/
  (hashfunc)0,			/*tp_hash*/
  (ternaryfunc)CCL_call,       	/*tp_call*/
  (reprfunc)0,			/*tp_str*/
  (getattrofunc)CCL_getattro,	/*tp_getattr with object key*/
  (setattrofunc)CCL_setattro,	/*tp_setattr with object key*/
  /* Space for future expansion */
  0L,0L,
  "Class of C classes" /* Documentation string */
};

/* End of code for ExtensionClass objects */
/* -------------------------------------------------------- */

/* subclassing code: */
  
static PyObject *
subclass_getspecial(PyObject *inst, PyObject *oname)
{
  PyObject *r=0;
  PyExtensionClass *self;

  self=(PyExtensionClass*)(inst->ob_type);
  if (HasInstDict(inst))
    {
      r= INSTANCE_DICT(inst);
      if (PyDict_Check(r))
        {
          if ((r = PyDict_GetItem(r,oname)))
            Py_INCREF(r);
          else
            r=CCL_getattr(self,oname,0);
        }
      else
        {
          UNLESS (r = PyObject_GetItem(r,oname))
            {
              PyErr_Clear();
              r=CCL_getattr(self,oname,0);
            }
        }
    }
  else 
    r=CCL_getattr(self,oname,0);
  
  return r;
}

static PyObject *
subclass_getattro(PyObject *self, PyObject *name)
{
  PyObject *r;

  if (! name) return NULL;
  UNLESS(r=EC_findiattro(self,name))
    {
      PyErr_Clear();
      r=EC_findiattro(self,py__getattr__);
      if (r) ASSIGN(r,PyObject_CallFunction(r,"O",name));
      if (r && NeedsToBeBound(r))
	ASSIGN(r, CallMethodO(r, py__of__, Build("(O)", self), NULL));
    }
  return r;
}

static int
subclass_simple_setattro(PyObject *self, PyObject *name, PyObject *v)
{
  if (! HasInstDict(self))
    {
      PyErr_SetObject(PyExc_AttributeError, name);
      return -1;
    }
  if (v)
    return PyDict_SetItem(INSTANCE_DICT(self),name,v);
  else
    return PyDict_DelItem(INSTANCE_DICT(self),name);
}

static int
subclass_simple_setattr(PyObject *self, char *name, PyObject *v)
{
  if (! HasInstDict(self))
    {
      PyErr_SetString(PyExc_AttributeError, name);
      return -1;
    }
  if (v)
    return PyDict_SetItemString(INSTANCE_DICT(self),name,v);
  else
    return PyDict_DelItemString(INSTANCE_DICT(self),name);
}

static int 
subclass_setattr(PyObject *self, PyObject *oname, char *name, PyObject *v)
{
  PyObject *m=0, *et, *ev, *etb;

  if (! name) return -1;

  if (!v && (m=subclass_getspecial(self,py__delattr__)))
    {
      if (UnboundEMethod_Check(m))
	{
	  UNLESS_ASSIGN(m,PyObject_CallFunction(m,"OO",self,oname)) return -1;
	}
      else UNLESS_ASSIGN(m,PyObject_CallFunction(m,"O",oname)) return -1;
      Py_DECREF(m);
      return 0;
    }

  UNLESS(m=subclass_getspecial(self,py__setattr__))
    goto default_setattr;
  if (UnboundCMethod_Check(m)
     && AsCMethod(m)->meth==(PyCFunction)setattr_by_name
     && SubclassInstance_Check(self,AsCMethod(m)->type))
    {
      UNLESS(-1 != AsCMethod(m)->type->tp_setattr(self,name,v))
	goto dictionary_setattr;
      return 0;
    }
  else 
    if (UnboundCMethod_Check(m)
       && AsCMethod(m)->meth==(PyCFunction)setattro_by_name
       && SubclassInstance_Check(self,AsCMethod(m)->type))
      {
	UNLESS(-1 != AsCMethod(m)->type->tp_setattro(self,oname,v))
	  goto dictionary_setattr;
	return 0;
      }
  if (! v) goto default_setattr;
  if (UnboundEMethod_Check(m))
    {
      UNLESS_ASSIGN(m,PyObject_CallFunction(m,"OOO",self,oname,v)) return -1;
    }
  else UNLESS_ASSIGN(m,PyObject_CallFunction(m,"OO",oname,v)) return -1;
  Py_DECREF(m);
  return 0;

dictionary_setattr:

  Py_XDECREF(m);

  PyErr_Fetch(&et, &ev, &etb);
  if (et==PyExc_AttributeError)
    {
      char *s;
      
      if (ev && PyString_Check(ev) && (s=PyString_AsString(ev)) &&
	 strcmp(s,name)==0)
	{
	  Py_XDECREF(et);
	  Py_XDECREF(ev);
	  Py_XDECREF(etb);
	  et=0;
	}
    }
  if (et)
    {
      PyErr_Restore(et,ev,etb);
      return -1;
    }	
  
default_setattr:

  PyErr_Clear();
  
  return subclass_simple_setattro(self, oname, v);
}

static int
subclass_setattro(PyObject *self, PyObject *name, PyObject *v)
{
  return subclass_setattr(self,name,PyString_AsString(name),v);
}


static int
subclass_compare(PyObject *self, PyObject *v)
{
  PyObject *m;
  long r;

  UNLESS(m=subclass_getspecial(self,py__cmp__))
    {
      PyErr_Clear();
      return self-v;
    }

  if (UnboundCMethod_Check(m)
     && AsCMethod(m)->meth==(PyCFunction)compare_by_name
     && SubclassInstance_Check(self,AsCMethod(m)->type)
     && ! HasMethodHook(self))
    r=AsCMethod(m)->type->tp_compare(self,v);
  else
    {
      if (UnboundEMethod_Check(m))
	{
	  UNLESS_ASSIGN(m,PyObject_CallFunction(m,"OO",self,v))
	    return -1;
	}
      else UNLESS_ASSIGN(m,PyObject_CallFunction(m,"O",v)) return -1;
      r=PyInt_AsLong(m);
    }
  Py_DECREF(m);
  return r;
}  

static long
subclass_hash(PyObject *self)
{
  PyObject *m;
  long r;

  UNLESS(m=subclass_getspecial(self,py__hash__)) return -1;
  if (UnboundCMethod_Check(m)
     && AsCMethod(m)->meth==(PyCFunction)hash_by_name
     && SubclassInstance_Check(self,AsCMethod(m)->type)
     && ! HasMethodHook(self))
    r=AsCMethod(m)->type->tp_hash(self);
  else
    {
      if (UnboundEMethod_Check(m))
	{
	  UNLESS_ASSIGN(m,PyObject_CallFunction(m,"O",self))
	    return -1;
	}
      else UNLESS_ASSIGN(m,PyObject_CallFunction(m,"")) return -1;
      r=PyInt_AsLong(m);
    }
  Py_DECREF(m);
  return r;
}  

static PyObject *
default_subclass_repr(PyObject *self)
{
  char p[128], *pp;
  
  PyErr_Clear();
  sprintf(p,"%p",self);
  if (*p=='0' && p[1]=='x') pp=p+2;
  else                      pp=p;
  return JimString_Build("<%s instance at %s>","ss",
			self->ob_type->tp_name, pp);
}

static PyObject *
subclass_repr(PyObject *self)
{
  PyObject *m;

  UNLESS(m=subclass_getspecial(self,py__repr__))
    return default_subclass_repr(self);

  if (UnboundCMethod_Check(m)
     && AsCMethod(m)->meth==(PyCFunction)repr_by_name
     && SubclassInstance_Check(self,AsCMethod(m)->type)
     && ! HasMethodHook(self))
    ASSIGN(m,AsCMethod(m)->type->tp_repr(self));
  else if (UnboundEMethod_Check(m))
    ASSIGN(m,PyObject_CallFunction(m,"O",self));
  else
    ASSIGN(m,PyObject_CallFunction(m,""));
  return m;
}  

static PyObject *
subclass_call(PyObject *self, PyObject *args, PyObject *kw)
{
  PyObject *m;

  UNLESS(m=subclass_getspecial(self,py__call__)) return NULL;
  if (UnboundCMethod_Check(m) && AsCMethod(m)->meth==(PyCFunction)call_by_name
     && SubclassInstance_Check(self,AsCMethod(m)->type)
     && ! HasMethodHook(self))
    ASSIGN(m,AsCMethod(m)->type->tp_call(self,args,kw));
  else
    {
      if (UnboundEMethod_Check(m))
	{
	  PyObject *a;
	  a=Py_BuildValue("(O)",self);
	  if (a) ASSIGN(a,PySequence_Concat(a,args));
	  if (a) ASSIGN(m,PyEval_CallObjectWithKeywords(m,a,kw));
	  else  ASSIGN(m,NULL);
	  Py_XDECREF(a);
	}
      else
	ASSIGN(m,PyEval_CallObjectWithKeywords(m,args,kw));
    }
  return m;
}  

static PyObject *
subclass_str(PyObject *self)
{
  PyObject *m;

  UNLESS(m=subclass_getspecial(self,py__str__))
    {
      PyErr_Clear();
      return subclass_repr(self);
    }
  if (UnboundCMethod_Check(m)
     && AsCMethod(m)->meth==(PyCFunction)str_by_name
     && SubclassInstance_Check(self,AsCMethod(m)->type)
     && ! HasMethodHook(self))
    ASSIGN(m,AsCMethod(m)->type->tp_str(self));
  else if (UnboundEMethod_Check(m))
    ASSIGN(m,PyObject_CallFunction(m,"O",self));
  else
    ASSIGN(m,PyObject_CallFunction(m,""));
  return m;
}  

#define BINSUB(M,N,A) \
static PyObject * \
subclass_ ## M(PyObject *self, PyObject *v) \
{ \
  PyObject *m; \
  UNLESS(m=subclass_getspecial(self,py__ ## N ## __)) return NULL; \
  if (UnboundCMethod_Check(m) \
     && AsCMethod(m)->meth==(PyCFunction)M ## _by_name \
     && SubclassInstance_Check(self,AsCMethod(m)->type) \
     && ! HasMethodHook(self)) \
    ASSIGN(m,AsCMethod(m)->type->tp_as_number->nb_ ## M(self,v)); \
  else if (UnboundEMethod_Check(m)) \
    ASSIGN(m,PyObject_CallFunction(m,"OO",self,v)); \
  else \
    ASSIGN(m,PyObject_CallFunction(m,"O",v)); \
  return m; \
}  

static PyObject * 
subclass_add(PyObject *self, PyObject *v)
{ 
  PyObject *m; 

  UNLESS(m=subclass_getspecial(self,py__add__)) return NULL; 

  if (UnboundCMethod_Check(m)
     && AsCMethod(m)->meth==(PyCFunction)concat_by_name 
     && SubclassInstance_Check(self,AsCMethod(m)->type)
     && ! HasMethodHook(self)) 
    ASSIGN(m,AsCMethod(m)->type->tp_as_sequence->sq_concat(self,v)); 
  else if (UnboundCMethod_Check(m)
	  && AsCMethod(m)->meth==(PyCFunction)add_by_name 
	  && SubclassInstance_Check(self,AsCMethod(m)->type)
	  && ! HasMethodHook(self)) 
    ASSIGN(m,AsCMethod(m)->type->tp_as_number->nb_add(self,v)); 
  else if (UnboundEMethod_Check(m)) 
    ASSIGN(m,PyObject_CallFunction(m,"OO",self,v)); 
  else 
    ASSIGN(m,PyObject_CallFunction(m,"O",v)); 

  return m; 
}  

BINSUB(subtract,sub,Subtract)

static PyObject * 
subclass_multiply(PyObject *self, PyObject *v)
{ 
  PyObject *m; 

  UNLESS(m=subclass_getspecial(self,py__mul__)) return NULL; 
  if (UnboundCMethod_Check(m)
     && AsCMethod(m)->meth==(PyCFunction)repeat_by_name 
     && SubclassInstance_Check(self,AsCMethod(m)->type)
     && ! HasMethodHook(self))
    {
      int i;

      i=PyInt_AsLong(v);
      if (i==-1 && PyErr_Occurred()) return NULL;
      ASSIGN(m,AsCMethod(m)->type->tp_as_sequence->sq_repeat(self,i));
    }
  else if (UnboundCMethod_Check(m)
	  && AsCMethod(m)->meth==(PyCFunction)multiply_by_name 
	  && SubclassInstance_Check(self,AsCMethod(m)->type)
	  && ! HasMethodHook(self)) 
    ASSIGN(m,AsCMethod(m)->type->tp_as_number->nb_multiply(self,v)); 
  else if (UnboundEMethod_Check(m)) 
    ASSIGN(m,PyObject_CallFunction(m,"OO",self,v)); 
  else 
    ASSIGN(m,PyObject_CallFunction(m,"O",v)); 
  return m; 
}  

BINSUB(divide,div,Divide)
BINSUB(remainder,mod,Remainder)

static PyObject * 
subclass_power(PyObject *self, PyObject *v, PyObject *w) 
{ 
  PyObject *m; 
  UNLESS(m=subclass_getspecial(self,py__pow__)) return NULL; 
  if (UnboundCMethod_Check(m)
     && AsCMethod(m)->meth==(PyCFunction)power_by_name
     && SubclassInstance_Check(self,AsCMethod(m)->type)
     && ! HasMethodHook(self)) 
    ASSIGN(m,AsCMethod(m)->type->tp_as_number->nb_power(self,v,w));
  else if (UnboundEMethod_Check(m))
    ASSIGN(m,PyObject_CallFunction(m,"OOO",self,v,w));
  else
    ASSIGN(m,PyObject_CallFunction(m,"OO",v,w)); 
  return m; 
}  

BINSUB(divmod,divmod,Divmod)
BINSUB(lshift,lshift,Lshift)
BINSUB(rshift,rshift,Rshift)
BINSUB(and,and,And)
BINSUB(or,or,Or)
BINSUB(xor,xor,Xor)


static int
subclass_coerce(PyObject **self, PyObject **v) 
{ 
  PyObject *m; 
  int r;

  UNLESS(m=subclass_getspecial(*self,py__coerce__))
    {
      PyErr_Clear();
      Py_INCREF(*self);
      Py_INCREF(*v);
      return 0;
    }
  if (UnboundCMethod_Check(m)
     && AsCMethod(m)->meth==(PyCFunction)coerce_by_name
     && SubclassInstance_Check(*self,AsCMethod(m)->type)
     && ! HasMethodHook(*self)) 
    r=AsCMethod(m)->type->tp_as_number->nb_coerce(self,v);
  else 
    { 
      if (UnboundEMethod_Check(m))
	ASSIGN(m,PyObject_CallFunction(m,"OO",*self,*v));
      else
	ASSIGN(m,PyObject_CallFunction(m,"O",*v));

      UNLESS (m) return -1;

      if (m==Py_None) r=-1;
      else
	{
	  if (PyArg_ParseTuple(m,"OO", self, v))
	    {
	      Py_INCREF(*self);
	      Py_INCREF(*v);
	      r=0;
	    }
	  else r=-1;
	}
    } 
  Py_DECREF(m);
  return r; 
}  

#define UNSUB(M,N) \
static PyObject * \
subclass_ ## M(PyObject *self) \
{ \
  PyObject *m; \
  UNLESS(m=subclass_getspecial(self,py__ ## N ## __)) return NULL; \
  if (UnboundCMethod_Check(m) \
     && AsCMethod(m)->meth==(PyCFunction)M ## _by_name \
     && SubclassInstance_Check(self,AsCMethod(m)->type) \
     && ! HasMethodHook(self)) \
    ASSIGN(m,AsCMethod(m)->type->tp_as_number->nb_ ## M(self)); \
  else if (UnboundEMethod_Check(m)) \
    ASSIGN(m,PyObject_CallFunction(m,"O",self)); \
  else \
    ASSIGN(m,PyObject_CallFunction(m,"")); \
  return m; \
}  

UNSUB(negative, neg)
UNSUB(positive, pos)
UNSUB(absolute, abs)

static int
subclass_nonzero(PyObject *self)
{
  PyObject *m;
  long r;

  UNLESS(m=subclass_getspecial(self,py__nonzero__))
    { /* We are being asked is we are true
	 Check out len, and if that fails, say we are true.
       */
      PyErr_Clear();
      UNLESS(m=subclass_getspecial(self,py__len__))
	{
	  PyErr_Clear();
	  return 1;
	}
    }
  if (UnboundCMethod_Check(m)
     && AsCMethod(m)->meth==(PyCFunction)nonzero_by_name
     && SubclassInstance_Check(self,AsCMethod(m)->type)
     && ! HasMethodHook(self))
    r=AsCMethod(m)->type->tp_as_number->nb_nonzero(self);
  else
    {
      if (UnboundEMethod_Check(m))
	{
	  UNLESS_ASSIGN(m,PyObject_CallFunction(m,"O",self))
	    return -1;
	}
      else UNLESS_ASSIGN(m,PyObject_CallFunction(m,"")) return -1;
      r=PyInt_AsLong(m);
    }
  Py_DECREF(m);
  return r;
}  

UNSUB(invert, inv)
UNSUB(int, int)
UNSUB(long, long)
UNSUB(float, float)
UNSUB(oct, oct)
UNSUB(hex, hex)

#undef UNSUB
#undef BINSUB


static PyNumberMethods subclass_as_number = {
  (binaryfunc)subclass_add,		/*nb_add*/
  (binaryfunc)subclass_subtract,	/*nb_subtract*/
  (binaryfunc)subclass_multiply,	/*nb_multiply*/
  (binaryfunc)subclass_divide,		/*nb_divide*/
  (binaryfunc)subclass_remainder,	/*nb_remainder*/
  (binaryfunc)subclass_divmod,		/*nb_divmod*/
  (ternaryfunc)subclass_power,		/*nb_power*/
  (unaryfunc)subclass_negative,		/*nb_negative*/
  (unaryfunc)subclass_positive,		/*nb_positive*/
  (unaryfunc)subclass_absolute,		/*nb_absolute*/
  (inquiry)subclass_nonzero,		/*nb_nonzero*/
  (unaryfunc)subclass_invert,		/*nb_invert*/
  (binaryfunc)subclass_lshift,		/*nb_lshift*/
  (binaryfunc)subclass_rshift,		/*nb_rshift*/
  (binaryfunc)subclass_and,		/*nb_and*/
  (binaryfunc)subclass_xor,		/*nb_xor*/
  (binaryfunc)subclass_or,		/*nb_or*/
  (coercion)subclass_coerce,		/*nb_coerce*/
  (unaryfunc)subclass_int,		/*nb_int*/
  (unaryfunc)subclass_long,		/*nb_long*/
  (unaryfunc)subclass_float,		/*nb_float*/
  (unaryfunc)subclass_oct,		/*nb_oct*/
  (unaryfunc)subclass_hex,		/*nb_hex*/
};

static long
subclass_length(PyObject *self)
{
  PyObject *m;
  long r;
  PyExtensionClass *t;

  UNLESS(m=subclass_getspecial(self,py__len__))
    {
      /* Hm. Maybe we are being checked to see if we are true.

	 Check to see if we have a __getitem__.  If we don't, then
	 answer that we are true.
       */
      PyErr_Clear();
      if ((m=subclass_getspecial(self,py__getitem__)))
	{
	  /* Hm, we have getitem, must be error */
	  Py_DECREF(m);
	  PyErr_SetObject(PyExc_AttributeError, py__len__);
	  return -1;
	}
      PyErr_Clear();
      return subclass_nonzero(self);
    }
  if (UnboundCMethod_Check(m)
     && AsCMethod(m)->meth==(PyCFunction)length_by_name
     && SubclassInstance_Check(self,AsCMethod(m)->type)
     && ! HasMethodHook(self))
    {
      t=(PyExtensionClass*)AsCMethod(m)->type;
      Py_DECREF(m);
      if (t->tp_as_sequence)
	return t->tp_as_sequence->sq_length(self);
      else
	return t->tp_as_mapping->mp_length(self);
    }
  if (UnboundEMethod_Check(m))
    {
      UNLESS_ASSIGN(m,PyObject_CallFunction(m,"O",self)) return -1;
    }
  else UNLESS_ASSIGN(m,PyObject_CallFunction(m,"")) return -1;
  r=PyInt_AsLong(m);
  Py_DECREF(m);
  return r;
}

static PyObject *
subclass_item(PyObject *self, int index)
{
  PyObject *m;
  PyExtensionClass *t;

  UNLESS(m=subclass_getspecial(self,py__getitem__)) return NULL;
  if (UnboundCMethod_Check(m)
     && AsCMethod(m)->meth==(PyCFunction)getitem_by_name
     && SubclassInstance_Check(self,AsCMethod(m)->type)
     && ! HasMethodHook(self))
    {
      t=(PyExtensionClass*)AsCMethod(m)->type;
      if (t->tp_as_sequence && t->tp_as_sequence->sq_item)
	{
	  Py_DECREF(m);
	  return t->tp_as_sequence->sq_item(self,index);
	}
    }
  if (UnboundEMethod_Check(m))
    ASSIGN(m,PyObject_CallFunction(m,"Oi",self,index));
  else
    ASSIGN(m,PyObject_CallFunction(m,"i",index));
  return m;
}

static PyObject *
subclass_slice(PyObject *self, int i1, int i2)
{
  PyObject *m;

  UNLESS(m=subclass_getspecial(self,py__getslice__)) return NULL;
  if (UnboundCMethod_Check(m)
     && AsCMethod(m)->meth==(PyCFunction)slice_by_name
     && SubclassInstance_Check(self,AsCMethod(m)->type)
     && ! HasMethodHook(self))
    ASSIGN(m,AsCMethod(m)->type->tp_as_sequence->sq_slice(self,i1,i2));
  else if (UnboundEMethod_Check(m))
    ASSIGN(m,PyObject_CallFunction(m,"Oii",self,i1,i2));
  else
    ASSIGN(m,PyObject_CallFunction(m,"ii",i1,i2));
  return m;
}

static long
subclass_ass_item(PyObject *self, int index, PyObject *v)
{
  PyObject *m;
  PyExtensionClass *t;

  if (! v && (m=subclass_getspecial(self,py__delitem__)))
    {
      if (UnboundEMethod_Check(m))
	{
	  UNLESS_ASSIGN(m,PyObject_CallFunction(m,"Oi",self,index)) return -1;
	}
      else UNLESS_ASSIGN(m,PyObject_CallFunction(m,"i",index)) return -1;
      Py_DECREF(m);
      return 0;
    }

  UNLESS(m=subclass_getspecial(self,py__setitem__)) return -1;
  if (UnboundCMethod_Check(m) &&
     AsCMethod(m)->meth==(PyCFunction)setitem_by_name
     && SubclassInstance_Check(self,AsCMethod(m)->type)
     && ! HasMethodHook(self))
    {
      t=(PyExtensionClass*)AsCMethod(m)->type;
      if (t->tp_as_sequence && t->tp_as_sequence->sq_ass_item)
	{
	  Py_DECREF(m);
	  return t->tp_as_sequence->sq_ass_item(self,index,v);
	}
    }
  if (! v)
    {
      PyErr_SetObject(PyExc_AttributeError, py__delitem__);
      return -1;
    }
  if (UnboundEMethod_Check(m))
    {
      UNLESS_ASSIGN(m,PyObject_CallFunction(m,"OiO",self,index,v)) return -1;
    }
  else UNLESS_ASSIGN(m,PyObject_CallFunction(m,"iO",index,v)) return -1;
  Py_DECREF(m);
  return 0;
}

static int
subclass_ass_slice(PyObject *self, int i1, int i2, PyObject *v)
{
  PyObject *m;
  long r;

  if (! v && (m=subclass_getspecial(self,py__delslice__)))
    {
      if (UnboundEMethod_Check(m))
	{
	  UNLESS_ASSIGN(m,PyObject_CallFunction(m,"Oii",self,i1,i2)) return -1;
	}
      else UNLESS_ASSIGN(m,PyObject_CallFunction(m,"ii",i1,i2)) return -1;
      Py_DECREF(m);
      return 0;
    }

  UNLESS(m=subclass_getspecial(self,py__setslice__)) return -1;
  if (UnboundCMethod_Check(m) &&
     AsCMethod(m)->meth==(PyCFunction)ass_slice_by_name
     && SubclassInstance_Check(self,AsCMethod(m)->type)
     && ! HasMethodHook(self))
    {	
      r=AsCMethod(m)->type->tp_as_sequence->sq_ass_slice(self,i1,i2,v);
      Py_DECREF(m);
      return r;
    }

  if (! v)
    {
      PyErr_SetObject(PyExc_AttributeError, py__delslice__);
      return -1;
    }

  if (UnboundEMethod_Check(m))
    {
      UNLESS_ASSIGN(m,PyObject_CallFunction(m,"OiiO",self,i1,i2,v))
	return -1;
    }
  else UNLESS_ASSIGN(m,PyObject_CallFunction(m,"iiO",i1,i2,v)) return -1;
  Py_DECREF(m);
  return 0;
}  

static PyObject *
subclass_repeat(PyObject *self, int v)
{
  PyObject *m;

  UNLESS(m=subclass_getspecial(self,py__mul__)) return NULL;
  if (UnboundCMethod_Check(m)
     && AsCMethod(m)->meth==(PyCFunction)repeat_by_name
     && SubclassInstance_Check(self,AsCMethod(m)->type)
     && ! HasMethodHook(self))
    ASSIGN(m,AsCMethod(m)->type->tp_as_sequence->sq_repeat(self,v));
  else if (UnboundEMethod_Check(m))
    ASSIGN(m,PyObject_CallFunction(m,"Oi",self,v));
  else
    ASSIGN(m,PyObject_CallFunction(m,"i",v));
  return m;
}

PySequenceMethods subclass_as_sequence = {
	(inquiry)subclass_length,   		/*sq_length*/
	(binaryfunc)subclass_add,		/*sq_concat*/
	(intargfunc)subclass_repeat,		/*sq_repeat*/
	(intargfunc)subclass_item,		/*sq_item*/
	(intintargfunc)subclass_slice,		/*sq_slice*/
	(intobjargproc)subclass_ass_item,	/*sq_ass_item*/
	(intintobjargproc)subclass_ass_slice,	/*sq_ass_slice*/
};

static PyObject *
subclass_subscript(PyObject *self, PyObject *key)
{
  PyObject *m;
  PyExtensionClass *t;

  UNLESS(m=subclass_getspecial(self,py__getitem__)) return NULL;
  if (UnboundCMethod_Check(m) &&
     AsCMethod(m)->meth==(PyCFunction)getitem_by_name
     && SubclassInstance_Check(self,AsCMethod(m)->type)
     && ! HasMethodHook(self))
    {
      t=(PyExtensionClass*)AsCMethod(m)->type;
      if (t->tp_as_mapping && t->tp_as_mapping->mp_subscript)
	{
	  Py_DECREF(m);
	  return t->tp_as_mapping->mp_subscript(self,key);
	}
      else if (t->tp_as_sequence && t->tp_as_sequence->sq_item)
	{
	  int i, l;

	  Py_DECREF(m);
	  
	  UNLESS(PyInt_Check(key))
	    {
	      PyErr_SetString(PyExc_TypeError, "sequence subscript not int");
	      return NULL;
	    }
	  i=PyInt_AsLong(key);
	  if (i < 0)
	    {
	      if ((l=PyObject_Length(self)) < 0) return NULL;
	      i+=l;
	    }
	  return t->tp_as_sequence->sq_item(self,i);
	}
    }
  if (UnboundEMethod_Check(m))
    ASSIGN(m,PyObject_CallFunction(m,"OO",self,key));
  else
    ASSIGN(m,PyObject_CallFunction(m,"O",key));
  return m;
}

static long
subclass_ass_subscript(PyObject *self, PyObject *index, PyObject *v)
{
  PyObject *m;
  PyExtensionClass *t;

  if (! v && (m=subclass_getspecial(self,py__delitem__)))
    {
      if (UnboundEMethod_Check(m))
	{
	  UNLESS_ASSIGN(m,PyObject_CallFunction(m,"OO",self,index)) return -1;
	}
      else UNLESS_ASSIGN(m,PyObject_CallFunction(m,"O",index)) return -1;
      Py_DECREF(m);
      return 0;
    }

  UNLESS(m=subclass_getspecial(self,py__setitem__)) return -1;
  if (UnboundCMethod_Check(m) &&
     AsCMethod(m)->meth==(PyCFunction)setitem_by_name
     && SubclassInstance_Check(self,AsCMethod(m)->type)
     && ! HasMethodHook(self))
    {
      t=(PyExtensionClass*)AsCMethod(m)->type;
      if (t->tp_as_mapping && t->tp_as_mapping->mp_ass_subscript)
	{
	  Py_DECREF(m);
	  return t->tp_as_mapping->mp_ass_subscript(self,index,v);
	}
      else if (t->tp_as_sequence && t->tp_as_sequence->sq_ass_item)
	{
	  int i, l;

	  Py_DECREF(m);
	  
	  UNLESS(PyInt_Check(index))
	    {
	      PyErr_SetString(PyExc_TypeError, "sequence subscript not int");
	      return -1;
	    }
	  i=PyInt_AsLong(index);
	  if (i < 0)
	    {
	      if ((l=PyObject_Length(self)) < 0) return -1;
	      i+=l;
	    }
	  return t->tp_as_sequence->sq_ass_item(self,i,v);
	}
    }
  if (! v)
    {
      PyErr_SetObject(PyExc_AttributeError, py__delitem__);
      return -1;
    }
  if (UnboundEMethod_Check(m))
    {
      UNLESS_ASSIGN(m,PyObject_CallFunction(m,"OOO",self,index,v)) return -1;
    }
  else UNLESS_ASSIGN(m,PyObject_CallFunction(m,"OO",index,v)) return -1;
  Py_DECREF(m);
  return 0;
}

PyMappingMethods subclass_as_mapping = {
	(inquiry)subclass_length,		/*mp_length*/
	(binaryfunc)subclass_subscript,		/*mp_subscript*/
	(objobjargproc)subclass_ass_subscript,	/*mp_ass_subscript*/
};

static int
dealloc_base(PyObject *inst, PyExtensionClass* self)
{
  int i,l;
  PyObject *t;

  l=PyTuple_Size(self->bases);
  for (i=0; i < l; i++)
    {
      t=PyTuple_GET_ITEM(self->bases, i);
      if (ExtensionClass_Check(t))
	{
	  if (AsExtensionClass(t)->bases)
	    {
	      if (dealloc_base(inst,AsExtensionClass(t))) return 1;
	    }
	  else
	    {
	      if (((PyExtensionClass*)t)->tp_dealloc)
		{
		  ((PyExtensionClass*)t)->tp_dealloc(inst);
		  return 1;
		}
	    }
	}
    }
  return 0;
}

static void
subclass_dealloc(PyObject *self)
{
  PyObject *m, *t, *v, *tb;
  int base_dealloced;

#ifdef TRACE_DEALLOC
  fprintf(stderr,"Deallocating a %s\n", self->ob_type->tp_name);
#endif

  PyErr_Fetch(&t,&v,&tb);
  Py_INCREF(self);		/* Give us a new lease on life */

  if (subclass_watcher &&
     ! PyObject_CallMethod(subclass_watcher,"destroying","O",self))
    PyErr_Clear();


  if ((m=subclass_getspecial(self,py__del__)))
    {
      if (UnboundEMethod_Check(m))
	ASSIGN(m,PyObject_CallFunction(m,"O",self));
      else
	ASSIGN(m,PyObject_CallFunction(m,""));
      Py_XDECREF(m);
    }

  PyErr_Clear();

  if (--self->ob_refcnt > 0)
    {
      PyErr_Restore(t,v,tb);
      return; /* we added a reference; don't delete now */
    }
  
  if (HasInstDict(self)) {
      Py_XDECREF(INSTANCE_DICT(self));
  }

  /* See if there was a dealloc handler in a (C) base class.
     If there was, then it deallocates the object and we
     get a true value back. 

     Note that if there *is* a base class dealloc, then
     *it* should decref the class.
  */
  base_dealloced=dealloc_base(self,(PyExtensionClass*)self->ob_type);

  /* We only deallocate ourselves if a base class didn't */
  UNLESS(base_dealloced) 
    {
      Py_DECREF(self->ob_type);
      PyObject_DEL(self);
    }

  PyErr_Restore(t,v,tb);
}

static void
datafull_baseclassesf(PyExtensionClass *type, PyObject **c1, PyObject **c2)
{
  /* Find the number of classes that have data and return them.
     There should be no more than one.
     */
  int l, i;
  PyObject *base;
  
  l=PyTuple_Size(type->bases);
  for (i=0; i < l && ! (*c1 && *c2); i++)
    {
      base=PyTuple_GET_ITEM(type->bases, i);
      if (ExtensionClass_Check(base))
	{
	  if (AsExtensionClass(base)->bases)
	    datafull_baseclassesf(AsExtensionClass(base),c1,c2);
	  else
	    {
	      if (AsExtensionClass(base)->tp_basicsize >
		 sizeof(PyPureMixinObject) ||
		 AsExtensionClass(base)->tp_itemsize > 0)
		{
		  if (! *c1)
		    *c1=base;
		  else if (*c1 != base)
		    *c2=base;
		}      
	    }
	}
    }
}

static int
datafull_baseclasses(PyExtensionClass *type)
{
  PyObject *c1=0, *c2=0;
  datafull_baseclassesf(type, &c1, &c2);
  if (c2) return 2;
  if (c1) return 1;
  return 0;
}

static PyObject *
datafull_baseclass(PyExtensionClass *type)
{
  /* Find the baseclass that has data and.  There should be only one. */
  int l, i;
  PyObject *base, *dbase;
  
  l=PyTuple_Size(type->bases);
  for (i=0; i < l; i++)
    {
      base=PyTuple_GET_ITEM(type->bases, i);
      if (ExtensionClass_Check(base))
	{
	  if (AsExtensionClass(base)->bases)
	    {
	      if ((dbase=datafull_baseclass(AsExtensionClass(base))))
		return dbase;
	    }
	  else
	    {
	      if (AsExtensionClass(base)->tp_basicsize >
		 sizeof(PyPureMixinObject) ||
		 AsExtensionClass(base)->tp_itemsize > 0)
		return base;
	    }
	}
    }
  return NULL;
}

static PyObject *
extension_baseclass(PyExtensionClass *type)
{
  /* Find the first immediate base class that is an extension class */
  int l, i;
  PyObject *base;
  
  l=PyTuple_Size(type->bases);
  for (i=0; i < l; i++)
    {
      base=PyTuple_GET_ITEM(type->bases, i);
      if (ExtensionClass_Check(base)) return base;
    }
  return JimErr_Format(PyExc_TypeError,
		      "No extension class found in subclass", NULL);
}

static int 
subclass_hasattr(PyExtensionClass *type, PyObject *name)
{
  PyObject *o;

  if ((o=CCL_getattro(type,name)))
    {
      Py_DECREF(o);
      return 1;
    }
  PyErr_Clear();
  return 0;
}

static void
subclass_init_getattr(PyExtensionClass *self, PyObject *methods)
{
  PyObject *m;

  if ((m=CCL_getattr(self,py__getattr__,0)))
    {
      if (UnboundCMethod_Check(m)
	 && AsCMethod(m)->meth==(PyCFunction)getattr_by_name
	 && Subclass_Check(self,AsCMethod(m)->type))
	{
	  self->tp_getattr=AsCMethod(m)->type->tp_getattr;
	}
      else if (UnboundCMethod_Check(m)
	      && AsCMethod(m)->meth==(PyCFunction)getattro_by_name
	      && Subclass_Check(self,AsCMethod(m)->type))
	{
	    self->tp_getattro=AsCMethod(m)->type->tp_getattro;
	}
      else
	{
	  PyObject_SetItem(methods,py__getattr__,m);
	  self->tp_getattro=subclass_getattro;
	}
      Py_DECREF(m);
    }
  else
    {
      PyErr_Clear();
      self->tp_getattro=EC_findiattro;
    }
}

static void
subclass_init_setattr(PyExtensionClass *self, PyObject *methods)
{
  PyObject *m;
  
  if ((m=CCL_getattr(self,py__setattr__,0)))
    {
      if (UnboundCMethod_Check(m)
	 && AsCMethod(m)->meth==(PyCFunction)setattr_by_name
	 && Subclass_Check(self,AsCMethod(m)->type))
	{
	  self->tp_setattr=AsCMethod(m)->type->tp_setattr;
	}
      else if (UnboundCMethod_Check(m)
	      && AsCMethod(m)->meth==(PyCFunction)setattro_by_name
	      && Subclass_Check(self,AsCMethod(m)->type))
	  {
	    self->tp_setattro=AsCMethod(m)->type->tp_setattro;
	  }
      else
	{
	  PyObject_SetItem(methods,py__setattr__,m);
	  self->tp_setattro=subclass_setattro;
	}
      Py_DECREF(m);
    }
  else
    {
      PyErr_Clear();
      self->tp_setattro=subclass_simple_setattro;
    }
}

static PyObject *
CopyMethods(PyExtensionClass *type, PyObject *base_methods)
{
  PyObject *methods, *key, *v;
  int pos;

  UNLESS(type->class_dictionary && PyDict_Check(base_methods) &&
	 ExtensionInstance_Check(type->class_dictionary))
    {
      Py_INCREF(base_methods);
      return base_methods;
    }

  UNLESS(methods=
	 PyObject_CallObject((PyObject*)type->class_dictionary->ob_type, NULL))
    return NULL;

  for (pos=0; PyDict_Next(base_methods, &pos, &key, &v); )
    UNLESS(0 <= PyObject_SetItem(methods,key,v)) goto err;

  return methods;

err:
  Py_DECREF(methods);
  return NULL;
}

/* Constructor for building subclasses of C classes.

   That is, we want to build a C class object that described a
   subclass of a built-in type.
 */
static PyObject *
subclass__init__(PyExtensionClass *self, PyObject *args)
{
  PyObject *bases, *methods, *class_init;
  PyExtensionClass *type;
  char *name, *p;
  int l;

  UNLESS(PyArg_ParseTuple(args,"sOO", &name, &bases, &methods)) return NULL;
  l=strlen(name)+1;
  UNLESS(p=(char*)malloc(l*sizeof(char))) return PyErr_NoMemory();
  memcpy(p,name,l);
  name=p;

  UNLESS(PyTuple_Check(bases) && PyTuple_Size(bases))
    {
      PyErr_SetString
	(PyExc_TypeError,
	 "second argument must be a tuple of 1 or more base classes");
    }

  self->bases=bases;
  Py_INCREF(bases);

  if (datafull_baseclasses(self) > 1)
    {
      PyErr_SetString(PyExc_TypeError, "too many datafull base classes");
      return NULL;
    }
  UNLESS(type=(PyExtensionClass *)datafull_baseclass(self))
    UNLESS(type=(PyExtensionClass *)extension_baseclass(self)) return NULL;
    
  self->tp_name=name;

  UNLESS(self->class_dictionary=CopyMethods(type,methods)) return NULL;

#define copy_member(M) self->M=type->M
  copy_member(ob_size);
  copy_member(class_flags);

  copy_member(tp_itemsize);
  copy_member(tp_print);
  self->tp_dealloc=subclass_dealloc;

  if (type->class_flags & EXTENSIONCLASS_PYTHONICATTR_FLAG)
    {
      /* The base class wants subclass __get/setattr__ to have
         Python class semantics and *it* will be providing them.
	 That means that we simply copy the base class 
	 get/setattr.
      */
      copy_member(tp_getattr);
      copy_member(tp_getattro);
      copy_member(tp_setattr);
      copy_member(tp_setattro);
      self->class_flags |= EXTENSIONCLASS_PYTHONICATTR_FLAG;

      if (CCL_hasattr(self, py__getattr__))
	self->class_flags |= EXTENSIONCLASS_USERGETATTR_FLAG;
      if (CCL_hasattr(self, py__setattr__))
	self->class_flags |= EXTENSIONCLASS_USERSETATTR_FLAG;
      if (CCL_hasattr(self, py__delattr__))
	self->class_flags |= EXTENSIONCLASS_USERDELATTR_FLAG;
    }
  else
    {
      subclass_init_getattr(self, methods);
      subclass_init_setattr(self, methods);
    }

#define subclass_set(OP,N) \
  self->tp_ ##OP = subclass_ ##OP
  
  subclass_set(compare,cmp);
  subclass_set(repr,repr);

  if (subclass_hasattr(self,py__of__))
    self->class_flags |= EXTENSIONCLASS_BINDABLE_FLAG;

  if (subclass_hasattr(self,py__call_method__))
    self->class_flags |= EXTENSIONCLASS_METHODHOOK_FLAG;

  UNLESS(self->class_flags & EXTENSIONCLASS_NOINSTDICT_FLAG)
    self->class_flags |= EXTENSIONCLASS_INSTDICT_FLAG;

  if (type->bases || ! ClassHasInstDict(self))
    copy_member(tp_basicsize);
  else
    {
      self->tp_basicsize=type->tp_basicsize/sizeof(PyObject*)*sizeof(PyObject*);
      if (self->tp_basicsize < type->tp_basicsize)
	self->tp_basicsize += sizeof(PyObject*); /* To align on PyObject */
      self->tp_basicsize += sizeof(PyObject*); /* For instance dictionary */
    }


  self->tp_as_number=(PyNumberMethods*)malloc(sizeof(PyNumberMethods));
  UNLESS(self->tp_as_number) return PyErr_NoMemory();
  *(self->tp_as_number)=subclass_as_number;
    
  self->tp_as_sequence=
    (PySequenceMethods*)malloc(sizeof(PySequenceMethods));
  UNLESS(self->tp_as_sequence) return PyErr_NoMemory();
  *(self->tp_as_sequence)=subclass_as_sequence;
  
  self->tp_as_mapping=(PyMappingMethods*)malloc(sizeof(PyMappingMethods));
  UNLESS(self->tp_as_mapping) return PyErr_NoMemory();
  *(self->tp_as_mapping)=subclass_as_mapping;

  subclass_set(hash,hash);
  subclass_set(call,call);
  subclass_set(str,str);
  self->tp_doc=0;

  /* Implement __module__=__name__ */
  if (PyDict_GetItem(methods, py__module__) == NULL)
    {
      PyObject *globals = PyEval_GetGlobals();
      if (globals != NULL)
	{
	  PyObject *modname = PyDict_GetItem(globals, py__name__);
	  if (modname != NULL) {
	    if (PyDict_SetItem(methods, py__module__, modname) < 0)
	      return NULL;
	  }
	}
    }

  /* Check for and use __class_init__ */
  if ((class_init=PyObject_GetAttrString(AsPyObject(self),"__class_init__")))
    {
      UNLESS_ASSIGN(class_init,PyObject_GetAttrString(class_init,"im_func"))
        return NULL;
      UNLESS_ASSIGN(class_init,PyObject_CallFunction(class_init,"O",self))
	return NULL;
      Py_DECREF(class_init);
    }
  else
    PyErr_Clear();

  Py_INCREF(Py_None);
  return Py_None;
}

struct PyMethodDef ExtensionClass_methods[] = {
  {"__init__",(PyCFunction)subclass__init__,1,""},
  {NULL,		NULL}		/* sentinel */
};

static PyExtensionClass ECType = {
  PyObject_HEAD_INIT(NULL)
  0,				/*ob_size*/
  "ExtensionClass",		/*tp_name*/
  sizeof(PyExtensionClass),    	/*tp_basicsize*/
  0,				/*tp_itemsize*/
  /* methods */
  (destructor)CCL_dealloc,	/*tp_dealloc*/
  (printfunc)0,			/*tp_print*/
  (getattrfunc)0,		/*tp_getattr*/
  (setattrfunc)0,		/*tp_setattr*/
  (cmpfunc)0,			/*tp_compare*/
  (reprfunc)CCL_repr,		/*tp_repr*/
  0,				/*tp_as_number*/
  0,				/*tp_as_sequence*/
  0,				/*tp_as_mapping*/
  (hashfunc)0,			/*tp_hash*/
  (ternaryfunc)CCL_call,       	/*tp_call*/
  (reprfunc)0,			/*tp_str*/
  (getattrofunc)CCL_getattro,	/*tp_getattr with object key*/
  (setattrofunc)CCL_setattro,	/*tp_setattr with object key*/
  /* Space for future expansion */
  0L,0L,
  "C classes", /* Documentation string */
  METHOD_CHAIN(ExtensionClass_methods)
};

/* List of methods defined in the module */

static PyObject *
set_subclass_watcher(PyObject *ignored, PyObject *args)
{
  PyObject *old, *sw=0;

  UNLESS(PyArg_ParseTuple(args,"|O",&sw)) return NULL;
  old=subclass_watcher;
  subclass_watcher=sw;
  if (sw) Py_INCREF(sw);
  if (old) return old;
  Py_INCREF(Py_None);
  return Py_None;
}

static struct PyMethodDef CC_methods[] = {
  {"subclass_watcher", (PyCFunction)set_subclass_watcher, 1,
   "subclass_watcher(ob) -- "
   "Register an object to watch subclass instance events"
  },
  {NULL,		NULL}		/* sentinel */
};

static int
export_type(PyObject *dict, char *name, PyExtensionClass *typ)
{
  initializeBaseExtensionClass(typ);

  if (PyErr_Occurred()) return -1;

  if (PyDict_GetItem(typ->class_dictionary, py__module__) == NULL)
    {
      PyObject *modname = PyDict_GetItem(dict, py__name__);
      if (modname != NULL) {
	if (PyDict_SetItem(typ->class_dictionary, py__module__, modname) < 0)
	  return -1;
      }
    }
  PyErr_Clear();
  
  return PyMapping_SetItemString(dict,name,(PyObject*)typ);
}

static struct ExtensionClassCAPIstruct
TrueExtensionClassCAPI = {
  export_type,			/* Export */
  EC_findiattrs,		/* getattrs */
  EC_findiattro,		/* getattro */
  subclass_simple_setattr,	/* setattrs */
  subclass_simple_setattro,	/* setattro */
  (PyObject*)&ECType,		/* ExtensionClassType */
  (PyObject*)&PMethodType,	/* MethodType */
  PMethod_New,			/* Method_New */
  CMethod_issubclass,		/* issubclass */
};

void
initExtensionClass(void)
{
  PyObject *m, *d, *s;
  PURE_MIXIN_CLASS(Base, "Minimalbase class for Extension Classes", NULL);

  PMethodType.ob_type=&PyType_Type;
  CMethodType.ob_type=&PyType_Type;
  ECTypeType.ob_type=&PyType_Type;
  ECType.ob_type=&ECTypeType;

  UNLESS(concat_fmt=PyString_FromString("%s%s"));
  
  m = Py_InitModule4("ExtensionClass", CC_methods,
		     ExtensionClass_module_documentation,
		     (PyObject*)NULL,PYTHON_API_VERSION);

  d = PyModule_GetDict(m);

  init_py_names();

  if (0) PyCObject_Import14("this will go away", "in 1.5 :-)");

  initializeBaseExtensionClass(&ECType);
  PyDict_SetItemString(d, "ExtensionClass", (PyObject*)&ECType);

  initializeBaseExtensionClass(&BaseType);
  PyDict_SetItemString(d, "Base", (PyObject*)&BaseType);

  PyDict_SetItemString(d, "PythonMethodType", (PyObject*)&PMethodType);
  PyDict_SetItemString(d, "ExtensionMethodType", (PyObject*)&CMethodType);

  /* Export C attribute lookup API */
  PyExtensionClassCAPI=&TrueExtensionClassCAPI;
  
  s = PyCObject_FromVoidPtr(PyExtensionClassCAPI, NULL);
  PyDict_SetItemString(d, "CAPI", s);
  Py_XDECREF(s);

  CHECK_FOR_ERRORS("can't initialize module ExtensionClass");
}
