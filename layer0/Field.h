

/* 
A* -------------------------------------------------------------------
B* This file contains source code for the PyMOL computer program
C* Copyright (c) Schrodinger, LLC. 
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
#ifndef _H_Field
#define _H_Field

#include"os_python.h"
#include"PyMOLGlobals.h"

#define cFieldFloat 0
#define cFieldInt 1
#define cFieldOther 2

struct CField {
  int type;
  std::vector<char> data;
  std::vector<unsigned int> dim;
  std::vector<unsigned int> stride;
  unsigned int base_size;
  CField() = default;
  CField(
      PyMOLGlobals* G, const int* const dim, int n_dim, unsigned int base_size, int type);
  int n_dim() const noexcept { return dim.size(); }
  unsigned int size() const noexcept { return data.size(); }

  /**
   * Copies data from another vector as stream of bytes
   * @param other source typed buffer
   */
  template <typename T> void set_data(const std::vector<T>& other)
  {
    auto nbytes = other.size() * sizeof(T);
    data.resize(nbytes);
    std::copy_n(
        reinterpret_cast<const char*>(other.data()), nbytes, data.begin());
  }

  /**
   * Returns heap-allocated handle to CField based on stored type
   * @param dim dimensions
   * @param n_dim number of dimensions
   * @return pointer to heap-allocated CField object.
   * Note: Lifetime of object (via heap-allocation) must be managed by owner.
   */
  template <typename T>
  static CField* make(PyMOLGlobals* G, const int* const dim, int n_dim)
  {
    int local_type = cFieldOther;
    if (std::is_same<T, float>::value) {
      local_type = cFieldFloat;
    } else if (std::is_same<T, int>::value) {
      local_type = cFieldInt;
    }
    return new CField(G, dim, n_dim, sizeof(T), local_type);
  }
};

/* accessors for getting data from a field */

#define F3p(f,a,b,c) ((f)->data.data() + \
        (a)*(f)->stride[0] + \
        (b)*(f)->stride[1] + \
        (c)*(f)->stride[2])

#define F4p(f,a,b,c,d) ((f)->data.data() + \
        (a)*(f)->stride[0] + \
        (b)*(f)->stride[1] + \
        (c)*(f)->stride[2] + \
        (d)*(f)->stride[3])

#define Ffloat3p(f,a,b,c) ((float*)F3p(f,a,b,c))
#define Ffloat3(f,a,b,c) (*(Ffloat3p(f,a,b,c)))

#define Ffloat4p(f,a,b,c,d) ((float*)F4p(f,a,b,c,d))
#define Ffloat4(f,a,b,c,d) (*(Ffloat4p(f,a,b,c,d)))

#define Fint3p(f,a,b,c) ((int*)F3p(f,a,b,c))
#define Fint3(f,a,b,c) (*(Fint3p(f,a,b,c)))

#define Fint4p(f,a,b,c,d) ((int*)F4p(f,a,b,c,d))
#define Fint4(f,a,b,c,d) (*(Fint4p(f,a,b,c,d)))

#define Fvoid4p(f,a,b,c,d) ((void*)F4p(f,a,b,c,d))

void FieldZero(CField * I);
float FieldInterpolatef(CField * I, int a, int b, int c, float x, float y, float z);
void FieldInterpolate3f(CField * I, int *locus, float *fract, float *result);

PyObject *FieldAsNumPyArray(CField * I, short copy);
PyObject *FieldAsPyList(PyMOLGlobals * G, CField * I);
CField *FieldNewFromPyList(PyMOLGlobals * G, PyObject * list);
CField *FieldNewFromPyList_From_List(PyMOLGlobals * G, PyObject * list, int);
int FieldSmooth3f(CField * I);

#endif
