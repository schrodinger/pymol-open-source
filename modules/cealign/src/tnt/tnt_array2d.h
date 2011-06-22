/*
*
* Template Numerical Toolkit (TNT)
*
* Mathematical and Computational Sciences Division
* National Institute of Technology,
* Gaithersburg, MD USA
*
*
* This software was developed at the National Institute of Standards and
* Technology (NIST) by employees of the Federal Government in the course
* of their official duties. Pursuant to title 17 Section 105 of the
* United States Code, this software is not subject to copyright protection
* and is in the public domain. NIST assumes no responsibility whatsoever for
* its use by other parties, and makes no guarantees, expressed or implied,
* about its quality, reliability, or any other characteristic.
*
*/



#ifndef TNT_ARRAY2D_H
#define TNT_ARRAY2D_H

#include <cstdlib>
#include <iostream>
#ifdef TNT_BOUNDS_CHECK
#include <assert.h>
#endif

#include "tnt_array1d.h"

namespace TNT
{


/**
	Ttwo-dimensional numerical array which
	looks like a conventional C multiarray. 
	Storage corresponds to C (row-major) ordering.
	Elements are accessed via A[i][j] notation for 0-based indexing,
	and A(i,j) for 1-based indexing.. 
	
	<p>
	Array assignment is by reference (i.e. shallow assignment).
	That is, B=A implies that the A and B point to the
	same array, so modifications to the elements of A
	will be reflected in B. If an independent copy
	is required, then B = A.copy() can be used.  Note
	that this facilitates returning arrays from functions
	without relying on compiler optimizations to eliminate
	extensive data copying.

	<p>
	The indexing and layout of this array object makes
	it compatible with C and C++ algorithms that utilize
	the familiar C[i][j] notation.  This includes numerous

	textbooks, such as Numercial Recipes, and various
	public domain codes.

*/
template <class T>
class Array2D 
{


  private:



  	Array1D<T> data_;
	Array1D<T*> v_;
	int m_;
    int n_;

  public:



/**
  Used to determined the data type of array entries.  This is most
  commonly used when requiring scalar temporaries in templated algorithms
  that have TNT arrays as input.  For example,
  <pre>
  template &lt class ArrayTwoD &gt
  void foo(ArrayTwoD &A)
  {
    A::value_type first_entry = A[0][0];
    ...
  }
  </pre>
*/
    typedef         T   value_type;

/**
  Create a null array.  This is <b>not</b> the same
  as Array2D(0,0), which consumes some memory overhead.
*/

	       Array2D();


/**
	Create a new (m x n) array, without initalizing elements.  (This
	encurs an O(1) operation cost, rather than a O(m*n) cost.)

  @param m the first (row) dimension of the new matrix.
  @param n the second (column) dimension of the new matrix.
*/
	       Array2D(int m, int n);


/**
  Create a new (m x n) array,  as a view of an existing one-dimensional
  array stored in row-major order, i.e. right-most dimension varying fastest.
  Note that the storage for this pre-existing array will
  never be destroyed by TNT.

  @param m the first (row) dimension of the new matrix.
  @param n the second (column) dimension of the new matrix.
  @param a the one dimensional C array to use as data storage for
    the array.
*/
	       Array2D(int m, int n,  T *a);



/**
  Create a new (m x n) array,  initializing array elements to
  constant specified by argument.  Most often used to
  create an array of zeros, as in A(m, n, 0.0).

  @param m the first (row) dimension of the new matrix.
  @param n the second (column) dimension of the new matrix.
  @param val the constant value to set all elements of the new array to.
*/

	       Array2D(int m, int n, const T &val);


/**
  Copy constructor. Array data is NOT copied, but shared.
  Thus, in Array2D B(A), subsequent changes to A will
  be reflected in B.  For an indepent copy of A, use
  Array2D B(A.copy()), or B = A.copy(), instead.
*/
    inline Array2D(const Array2D &A);


/**
	Convert 2D array into a regular multidimensional C pointer.  Most often
	called automatically when calling C interfaces that expect things like
	double** rather than Array2D<dobule>.
	
*/
	inline operator T**();

/**
	Convert a const 2D array into a const multidimensional C pointer.  
	Most often called automatically when calling C interfaces that expect 
	things like "const double**" rather than "const Array2D<dobule>&".
	
*/

	inline operator const T**() const;


/**
	Assign all elements of array the same value.

	@param val the value to assign each element.
*/
	inline Array2D & operator=(const T &val);


/**
	Assign one Array2D to another.  (This is a shallow-assignement operation,
	and it is the identical semantics to ref(A).

	@param A the array to assign this one to.
*/
	inline Array2D & operator=(const Array2D &A);


	inline Array2D & ref(const Array2D &A);
	       Array2D copy() const;
		   Array2D & inject(const Array2D & A);
	inline T* operator[](int i);
	inline const T* operator[](int i) const;
	inline int dim1() const;
	inline int dim2() const;
     ~Array2D();

	/* extended interface (not part of the standard) */


	inline int ref_count();
	inline int ref_count_data();
	inline int ref_count_dim1();
	Array2D subarray(int i0, int i1, int j0, int j1);

};


/**
  Create a new (m x n) array, WIHOUT initializing array elements.
  To create an initialized array of constants, see Array2D(m,n,value).

  <p>
  This version avoids the O(m*n) initialization overhead and
  is used just before manual assignment.

  @param m the first (row) dimension of the new matrix.
  @param n the second (column) dimension of the new matrix.
*/
template <class T>
Array2D<T>::Array2D() : data_(), v_(), m_(0), n_(0) {} 


template <class T>
Array2D<T>::Array2D(const Array2D<T> &A) : data_(A.data_), v_(A.v_), 
	m_(A.m_), n_(A.n_) {}




template <class T>
Array2D<T>::Array2D(int m, int n) : data_(m*n), v_(m), m_(m), n_(n)
{
	if (m>0 && n>0)
	{
		T* p = &(data_[0]);
		for (int i=0; i<m; i++)
		{
			v_[i] = p;
			p += n;
		}
	}
}



template <class T>
Array2D<T>::Array2D(int m, int n, const T &val) : data_(m*n), v_(m), 
													m_(m), n_(n) 
{
  if (m>0 && n>0)
  {
	data_ = val;
	T* p  = &(data_[0]);
	for (int i=0; i<m; i++)
	{
			v_[i] = p;
			p += n;
	}
  }
}

template <class T>
Array2D<T>::Array2D(int m, int n, T *a) : data_(m*n, a), v_(m), m_(m), n_(n)
{
  if (m>0 && n>0)
  {
	T* p = &(data_[0]);
	
	for (int i=0; i<m; i++)
	{
			v_[i] = p;
			p += n;
	}
  }
}


template <class T>
inline T* Array2D<T>::operator[](int i) 
{ 
#ifdef TNT_BOUNDS_CHECK
	assert(i >= 0);
	assert(i < m_);
#endif

return v_[i]; 

}


template <class T>
inline const T* Array2D<T>::operator[](int i) const
{ 
#ifdef TNT_BOUNDS_CHECK
	assert(i >= 0);
	assert(i < m_);
#endif

return v_[i]; 

}

template <class T>
Array2D<T> & Array2D<T>::operator=(const T &a)
{
	/* non-optimzied, but will work with subarrays in future verions */

	for (int i=0; i<m_; i++)
		for (int j=0; j<n_; j++)
		v_[i][j] = a;
	return *this;
}




template <class T>
Array2D<T> Array2D<T>::copy() const
{
	Array2D A(m_, n_);

	for (int i=0; i<m_; i++)
		for (int j=0; j<n_; j++)
			A[i][j] = v_[i][j];


	return A;
}


template <class T>
Array2D<T> & Array2D<T>::inject(const Array2D &A)
{
	if (A.m_ == m_ &&  A.n_ == n_)
	{
		for (int i=0; i<m_; i++)
			for (int j=0; j<n_; j++)
				v_[i][j] = A[i][j];
	}
	return *this;
}




template <class T>
Array2D<T> & Array2D<T>::ref(const Array2D<T> &A)
{
	if (this != &A)
	{
		v_ = A.v_;
		data_ = A.data_;
		m_ = A.m_;
		n_ = A.n_;
		
	}
	return *this;
}



template <class T>
Array2D<T> & Array2D<T>::operator=(const Array2D<T> &A)
{
	return ref(A);
}

template <class T>
inline int Array2D<T>::dim1() const { return m_; }

template <class T>
inline int Array2D<T>::dim2() const { return n_; }


template <class T>
Array2D<T>::~Array2D() {}




template <class T>
inline Array2D<T>::operator T**()
{
	return &(v_[0]);
}
template <class T>
inline Array2D<T>::operator const T**() const
{
	return static_cast<const T**>(&(v_[0]));
}

/* ............... extended interface ............... */
/**
	Create a new view to a subarray defined by the boundaries
	[i0][i0] and [i1][j1].  The size of the subarray is
	(i1-i0) by (j1-j0).  If either of these lengths are zero
	or negative, the subarray view is null.

*/
template <class T>
Array2D<T> Array2D<T>::subarray(int i0, int i1, int j0, int j1) 
{
	Array2D<T> A;
	int m = i1-i0+1;
	int n = j1-j0+1;

	/* if either length is zero or negative, this is an invalide
		subarray. return a null view.
	*/
	if (m<1 || n<1)
		return A;

	A.data_ = data_;
	A.m_ = m;
	A.n_ = n;
	A.v_ = Array1D<T*>(m);
	T* p = &(data_[0]) + i0 *  n_ + j0;
	for (int i=0; i<m; i++)
	{
		A.v_[i] = p + i*n_;

	}	
	return A;
}

template <class T>
inline int Array2D<T>::ref_count()
{
	return ref_count_data();
}



template <class T>
inline int Array2D<T>::ref_count_data()
{
	return data_.ref_count();
}

template <class T>
inline int Array2D<T>::ref_count_dim1()
{
	return v_.ref_count();
}




} /* namespace TNT */

#endif
/* TNT_ARRAY2D_H */

