

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

#include <cassert>

enum cField_t {
  cFieldFloat = 0,
  cFieldInt = 1,
  cFieldOther = 2,
};

/**
 * Multi-dimensional data array with runtime typing.
 */
struct CField {
  cField_t type;
  std::vector<char> data;
  std::vector<unsigned int> dim;
  std::vector<unsigned int> stride;
  unsigned int base_size;
  CField() = default;
  CField(
      PyMOLGlobals* G, const int* const dim, int n_dim, unsigned int base_size, cField_t type);
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

protected:
  template <typename T> static cField_t _get_type()
  {
    if (std::is_same<T, float>::value) {
      return cFieldFloat;
    }
    if (std::is_same<T, int>::value) {
      return cFieldInt;
    }
    return cFieldOther;
  }

public:
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
    return new CField(G, dim, n_dim, sizeof(T), _get_type<T>());
  }

private:
  size_t _data_offset(size_t a, size_t b, size_t c) const
  {
    return a * stride[0] + b * stride[1] + c * stride[2];
  }

  size_t _data_offset(size_t a, size_t b, size_t c, size_t d) const
  {
    return a * stride[0] + b * stride[1] + c * stride[2] + d * stride[3];
  }

public:
  //! Get a pointer to the value at `pos`
  template <typename T, typename... SizeTs> T* ptr(SizeTs... pos)
  {
    assert(sizeof...(pos) <= n_dim());
    return reinterpret_cast<T*>(data.data() + _data_offset(pos...));
  }

  template <typename T, typename... SizeTs> T const* ptr(SizeTs... pos) const
  {
    return const_cast<CField*>(this)->ptr<T const>(pos...);
  }

  //! Get the value at `pos`
  template <typename T, typename... SizeTs> T& get(SizeTs... pos)
  {
    assert(sizeof...(pos) == n_dim());
    assert(sizeof(T) == base_size);
    return *ptr<T>(pos...);
  }

  template <typename T, typename... SizeTs> T const& get(SizeTs... pos) const
  {
    return const_cast<CField*>(this)->get<T const>(pos...);
  }
};

/**
 * Multi-dimensional data array with compile-time typing.
 */
template <typename T> struct CFieldTyped : CField {
  CFieldTyped(const int* const dim, int n_dim)
      : CField(nullptr, dim, n_dim, sizeof(T), _get_type<T>())
  {
  }

  template <typename... SizeTs> T* ptr(SizeTs... pos)
  {
    return CField::ptr<T>(pos...);
  }

  template <typename... SizeTs> T const* ptr(SizeTs... pos) const
  {
    return CField::ptr<T>(pos...);
  }

  template <typename... SizeTs> T& get(SizeTs... pos)
  {
    return CField::get<T>(pos...);
  }

  template <typename... SizeTs> T const& get(SizeTs... pos) const
  {
    return CField::get<T>(pos...);
  }
};

/* accessors for getting data from a field */

#define Ffloat3p(f, a, b, c) ((f)->ptr<float>(a, b, c))
#define Ffloat3(f, a, b, c) ((f)->get<float>(a, b, c))

#define Ffloat4p(f, a, b, c, d) ((f)->ptr<float>(a, b, c, d))
#define Ffloat4(f, a, b, c, d) ((f)->get<float>(a, b, c, d))

#define F3 Ffloat3
#define F3Ptr Ffloat3p

#define F4 Ffloat4
#define F4Ptr Ffloat4p

void FieldZero(CField * I);
float FieldInterpolatef(CField * I, int a, int b, int c, float x, float y, float z);
void FieldInterpolate3f(CField * I, int *locus, float *fract, float *result);

PyObject *FieldAsNumPyArray(CField * I, short copy);
PyObject *FieldAsPyList(PyMOLGlobals * G, CField * I);
CField *FieldNewFromPyList(PyMOLGlobals * G, PyObject * list);
CField *FieldNewFromPyList_From_List(PyMOLGlobals * G, PyObject * list, int);
int FieldSmooth3f(CField * I);

#endif
