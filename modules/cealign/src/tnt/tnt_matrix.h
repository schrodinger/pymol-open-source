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


// C compatible matrix: row-oriented, 0-based [i][j] and 1-based (i,j) indexing
//

#ifndef TNT_MATRIX_H
#define TNT_MATRIX_H

#include "tnt_subscript.h"
#include "tnt_vector.h"
#include <cstdlib>
#include <cassert>
#include <iostream>
#include <sstream>
#include <fstream>

namespace TNT
{

//namespace Linear_Algebra {


/**
		Dense matrix class for basic linear algebra operations.

		Ordering: row major.
		<p>
		Elements begin at (1,1) or [0][0].
		<p>
		Can be interfaced with C multidimentionsal arrays (e.g. double **)
		<p>
		copy-by-value semantics.
		<p>
		Optional range checking at compile time via TNT_BOUNDS_CHECK macro.

*/
template <class T>
class Matrix 
{

  private:
    Subscript m_;
    Subscript n_;
    Subscript mn_;      // total size
    T* v_;                  
    T** row_;           
    T* vm1_ ;       // these point to the same data, but are 1-based 
    T** rowm1_;

    // internal helper function to create the array
    // of row pointers

    void initialize(Subscript M, Subscript N)
    {
        mn_ = M*N;
        m_ = M;
        n_ = N;

        v_ = new T[mn_]; 
        row_ = new T*[M];
        rowm1_ = new T*[M];

        assert(v_  != NULL);
        assert(row_  != NULL);
        assert(rowm1_ != NULL);

        T* p = v_;              
        vm1_ = v_ - 1;
        for (Subscript i=0; i<M; i++)
        {
            row_[i] = p;
            rowm1_[i] = p-1;
            p += N ;
            
        }

        rowm1_ -- ;     // compensate for 1-based offset
    }
   
    void copy(const T*  v)
    {
        Subscript N = m_ * n_;
        Subscript i;

#ifdef TNT_UNROLL_LOOPS
        Subscript Nmod4 = N & 3;
        Subscript N4 = N - Nmod4;

        for (i=0; i<N4; i+=4)
        {
            v_[i] = v[i];
            v_[i+1] = v[i+1];
            v_[i+2] = v[i+2];
            v_[i+3] = v[i+3];
        }

        for (i=N4; i< N; i++)
            v_[i] = v[i];
#else

        for (i=0; i< N; i++)
            v_[i] = v[i];
#endif      
    }

    void set(const T& val)
    {
        Subscript N = m_ * n_;
        Subscript i;

#ifdef TNT_UNROLL_LOOPS
        Subscript Nmod4 = N & 3;
        Subscript N4 = N - Nmod4;

        for (i=0; i<N4; i+=4)
        {
            v_[i] = val;
            v_[i+1] = val;
            v_[i+2] = val;
            v_[i+3] = val; 
        }

        for (i=N4; i< N; i++)
            v_[i] = val;
#else

        for (i=0; i< N; i++)
            v_[i] = val;
        
#endif      
    }
    

    
    void destroy()
    {     
        /* do nothing, if no memory has been previously allocated */
        if (v_ == NULL) return ;

        /* if we are here, then matrix was previously allocated */
        if (v_ != NULL) delete [] (v_);     
        if (row_ != NULL) delete [] (row_);

        /* return rowm1_ back to original value */
        rowm1_ ++;
        if (rowm1_ != NULL ) delete [] (rowm1_);
    }

  public:

    typedef Subscript   size_type;
    typedef         T   value_type;
    typedef         T   element_type;
    typedef         T*  pointer;
    typedef         T*  iterator;
    typedef         T&  reference;
    typedef const   T*  const_iterator;
    typedef const   T&  const_reference;

    Subscript lbound() const { return 1;}
 



    operator T**(){ return  row_; }
    operator const T**() const { return row_; }


		/**
			@return the total number of items in matrix (M*N).
		*/
    Subscript size() const { return mn_; }

    // constructors

    Matrix() : m_(0), n_(0), mn_(0), v_(0), row_(0), vm1_(0), rowm1_(0) {};

    Matrix(const Matrix<T> &A)
    {
        initialize(A.m_, A.n_);
        copy(A.v_);
    }

		/**
			Create a MxN matrix, with each element assigned to the value 0.

			@param M the number of rows
			@param N the number of columns
			@param value (optional default value: 0 if not specified.

		*/                                              
    Matrix(Subscript M, Subscript N, const T& value = T(0))
    {
        initialize(M,N);
        set(value);
    }

		/**
			Create an MxN matrix, filling in values (row-major order) from
			the list (C array) provided.

			@param M the number of rows
			@param N the number of columns
			@param v  list (C array) of M*N values used to initialize matrix.
		*/
    Matrix(Subscript M, Subscript N, const T* v)
    {
        initialize(M,N);
        copy(v);
    }

		/**
			Create an MxN matrix, filling in values (row-major order) from
			a character string.

			@param M the number of rows
			@param N the number of columns
			@param s  string of M*N values used to initialize matrix.
		*/
    Matrix(Subscript M, Subscript N, const char *s)
    {
        initialize(M,N);
        //std::istrstream ins(s);
        std::istringstream ins(s);

        Subscript i, j;

        for (i=0; i<M; i++)
            for (j=0; j<N; j++)
                ins >> row_[i][j];
    }

    // destructor
    //
    ~Matrix()
    {
        destroy();
    }


    /** 
				Change size of matrix to MxN, reallocating memory if necessary.
				<p>
				NOTE: This operations occurs in place, i.e. when resizing to 
				a new matrix, original matrix elements
				are <b>NOT</b> retained.  Instead, one must explicit create
				a new matrix of this size and manually copy the elements, e.g.
				<pre>

				Matrix double B(M, N);

				int 	min_M = M < A.num_rows() ? M : A.num_rows();
				int 	min_N = N < A.num_cols() ? N : A.num_cols();
				for (int i=1; i<=min_M; i++)
					for (int j=1; j<=min_N; j++)
						B(i,j) = A(i,j);

				A.destroy();
				</pre>

				@param M  the number of rows of new size.
				@param N	the number of columns of new size.
    */
    Matrix<T>& newsize(Subscript M, Subscript N)
    {
        if (num_rows() == M && num_cols() == N)
            return *this;

        destroy();
        initialize(M,N);
        
        return *this;
    }




    /**
			Assign (copy) one matrix to another, e.g. A=B.  The
			contents of A are lost, and a new copy of B is created.  

			@param B to matrix to be copied.

    */
    Matrix<T>& operator=(const Matrix<T> &B)
    {
        if (v_ == B.v_)
            return *this;

        if (m_ == B.m_  && n_ == B.n_)      // no need to re-alloc
            copy(B.v_);

        else
        {
            destroy();
            initialize(B.m_, B.n_);
            copy(B.v_);
        }

        return *this;
    }
        
    Matrix<T>& operator=(const T& scalar)
    { 
        set(scalar); 
        return *this;
    }


    Subscript dim(Subscript d) const 
    {
#ifdef TNT_BOUNDS_CHECK
        assert( d >= 1);
        assert( d <= 2);
#endif
        return (d==1) ? m_ : ((d==2) ? n_ : 0); 
    }

    Subscript num_rows() const { return m_; }
    Subscript num_cols() const { return n_; }




    inline T* operator[](Subscript i)
    {
#ifdef TNT_BOUNDS_CHECK
        assert(0<=i);
        assert(i < m_) ;
#endif
        return row_[i];
    }

    inline const T* operator[](Subscript i) const
    {
#ifdef TNT_BOUNDS_CHECK
        assert(0<=i);
        assert(i < m_) ;
#endif
        return row_[i];
    }

    inline reference operator()(Subscript i)
    { 
#ifdef TNT_BOUNDS_CHECK
        assert(1<=i);
        assert(i <= mn_) ;
#endif
        return vm1_[i]; 
    }

    inline const_reference operator()(Subscript i) const
    { 
#ifdef TNT_BOUNDS_CHECK
        assert(1<=i);
        assert(i <= mn_) ;
#endif
        return vm1_[i]; 
    }



    inline reference operator()(Subscript i, Subscript j)
    { 
#ifdef TNT_BOUNDS_CHECK
        assert(1<=i);
        assert(i <= m_) ;
        assert(1<=j);
        assert(j <= n_);
#endif
        return  rowm1_[i][j]; 
    }


    
    inline const_reference operator() (Subscript i, Subscript j) const
    {
#ifdef TNT_BOUNDS_CHECK
        assert(1<=i);
        assert(i <= m_) ;
        assert(1<=j);
        assert(j <= n_);
#endif
        return rowm1_[i][j]; 
    }



		Vector<T> diag() const
		{
			Subscript N = n_ < m_ ? n_ : m_; 
			Vector<T> d(N);

			for (int i=0; i<N; i++)
				d[i] = row_[i][i];

		  return d;
		}

	  Matrix<T> upper_triangular() const;
		Matrix<T> lower_triangular() const;

	/* basic computations */

};


/* ***************************  I/O  ********************************/


template <class T>
std::ostream& operator<<(std::ostream &s, const Matrix<T> &A)
{
    Subscript M=A.num_rows();
    Subscript N=A.num_cols();

    s << M << " " << N << "\n";
    for (Subscript i=0; i<M; i++)
    {
        for (Subscript j=0; j<N; j++)
        {
            s << A[i][j] << " ";
        }
        s << "\n";
    }


    return s;
}


template <class T>
std::istream& operator>>(std::istream &s, Matrix<T> &A)
{

    Subscript M, N;

    s >> M >> N;

    if ( !(M == A.num_rows() && N == A.num_cols() ))
    {
        A.newsize(M,N);
    }


    for (Subscript i=0; i<M; i++)
        for (Subscript j=0; j<N; j++)
        {
            s >>  A[i][j];
        }


    return s;
}



// *******************[ basic matrix algorithms ]***************************

/**
		Matrix-Matrix multiplication:  C = A * B.

		<p>
		This is an optimizied (trinary) version of matrix multiply, where 
		the destination matrix has already been allocated.

		@param A	matrix of size M x N.
		@param B	matrix of size N x K.
		@param C the result A*B, of size M x K.

		@return a reference to C, after multiplication.
*/
template <class T>
Matrix<T> & mult(Matrix<T>& C, const Matrix<T>  &A, const Matrix<T> &B)
{

#ifdef TNT_BOUNDS_CHECK
    assert(A.num_cols() == B.num_rows());
#endif

    Subscript M = A.num_rows();
    Subscript N = A.num_cols();
    Subscript K = B.num_cols();

#ifdef TNT_BOUNDS_CHECK
		assert(C.num_rows() == M);
		assert(C.num_cols() == K);
#endif	

    T sum;

    const T* row_i;
    const T* col_k;

    for (Subscript i=0; i<M; i++)
    for (Subscript k=0; k<K; k++)
    {
        row_i  = &(A[i][0]);
        col_k  = &(B[0][k]);
        sum = 0;
        for (Subscript j=0; j<N; j++)
        {
            sum  += *row_i * *col_k;
            row_i++;
            col_k += K;
        }
        C[i][k] = sum; 
    }

    return C;
}


/**
	Matrix/matrix multiplication: compute A * B.

	@param A matrix: left side operand  (size M x N).
	@param B matrix: right side operand  (size N x K).

	@return A*B  (a new matirx of size M x K).
*/
template <class T>
inline Matrix<T> mult(const Matrix<T>  &A, const Matrix<T> &B)
{

#ifdef TNT_BOUNDS_CHECK
    assert(A.num_cols() == B.num_rows());
#endif

    Subscript M = A.num_rows();
    Subscript K = B.num_cols();

    Matrix<T> tmp(M,K);

		mult(tmp, A, B);		// tmp = A*B

    return tmp;
}

/**
	Matrix/matrix multiplication: compute A * B.

	@param A matrix: left side operand  (size M x N).
	@param B matrix: right side operand  (size N x K).

	@return A*B  (a new matirx of size M x K).
*/
template <class T>
inline Matrix<T> operator*(const Matrix<T>  &A, const Matrix<T> &B)
{
	return mult(A,B);
}

/**
	Matrix/vector multiplication: compute A * b.

	@param A matrix: left side operand  (number of columns of A, must match 
			the number of elements in b.)
	@param b vector: right side operand.
	@return A*b (a new vector of size M.)
*/
	
template <class T>
inline Vector<T> mult(const Matrix<T>  &A, const Vector<T> &b)
{

#ifdef TNT_BOUNDS_CHECK
    assert(A.num_cols() == b.dim());
#endif

    Subscript M = A.num_rows();
    Subscript N = A.num_cols();

    Vector<T> tmp(M);
    T sum;

    for (Subscript i=0; i<M; i++)
    {
        sum = 0;
        for (Subscript j=0; j<N; j++)
            sum = sum +  A[i][j] * b[j];

        tmp[i] = sum; 
    }

    return tmp;
}

/**
	Matrix/vector multiplication: compute A * b.

	@param A matrix: left side operand  (number of columns of A, must match 
			the number of elements in b.)
	@param b vector: right side operand.
	@return A*b (a new vector of size M.)
*/
	
template <class T>
inline Vector<T> operator*(const Matrix<T>  &A, const Vector<T> &b)
{
	return mult(A,b);
}


/**
	Matrix scaling: multiply each element of A by scalar s.

	<p>NOTE: this creates a new copy of A.  To scale "in place",
	use *= or mult_eq().


	@param A matrix: to be scaled. 
	@param s scalar: multiplier.
	@return s*A, a new matrix with same size of A.
*/
template <class T>
inline Matrix<T> mult(const T& s, const Matrix<T> &A)
{
    Subscript M = A.num_rows();
    Subscript N = A.num_cols();

		Matrix<T> R(M,N);
		for (int i=0; i<M; i++)
			for (int j=0; j<N; j++)
				R[i][j] = s * A[i][j];
				 
		return R;
}

/**
	Matrix scaling: multiply each element of A by scalar s.

	<p>
	Same as mult(A,s), as this is a commutative operation.

	<p>NOTE: this creates a new copy of A.  To scale "in place",
	use *= or mult_eq().


	@param A matrix: to be scaled. 
	@param s scalar: multiplier.
	@return s*A, a new matrix with same size of A.
*/
template <class T>
inline Matrix<T> mult(const Matrix<T> &A, const T& s)
{
	return mult(s, A);
}


/**
	Matrix scale in-place, i.e. compute A *= s, where each element
	of A is multiplied (scaled) by the value s.

	<p>NOTE: this creates a new copy of A.  To scale "in place",
	use *= or mult_eq().


	@param A matrix: to be scaled. 
	@param s scalar: multiplier.
	@return A, after scaling.
*/
template <class T>
inline Matrix<T> mult_eq(const T& s, Matrix<T> &A)
{
    Subscript M = A.num_rows();
    Subscript N = A.num_cols();

		for (int i=0; i<M; i++)
			for (int j=0; j<N; j++)
				A[i][j] *= s;
}

template <class T>
inline Matrix<T> mult_eq(Matrix<T> &A, const T&a)
{
	return mult_eq(a, A);
}


/**
	Matrix-Matrix tranpose multiplication, i.e. compute tranpose(A)*B.

	<p>
	NOTE: this is more efficient than computing the tranpose(A) explicitly,
	and then multiplying, as the tranpose of A is never really constructed.

	@param A  matrix: size M x N.
	@param B	matrix: size M x K.
	@return a new matrix of size N x K.
*/
template <class T>
inline Matrix<T> transpose_mult(const Matrix<T>  &A, const Matrix<T> &B)
{

#ifdef TNT_BOUNDS_CHECK
    assert(A.num_rows() == B.num_rows());
#endif

    Subscript M = A.num_cols();
    Subscript N = A.num_rows();
    Subscript K = B.num_cols();

    Matrix<T> tmp(M,K);
    T sum;

    for (Subscript i=0; i<N; i++)
    for (Subscript k=0; k<K; k++)
    {
        sum = 0;
        for (Subscript j=0; j<M; j++)
            sum = sum +  A[j][i] * B[j][k];

        tmp[i][k] = sum; 
    }

    return tmp;
}

/**
	Matrix-Vector tranpose multiplication, i.e. compute tranpose(A)*b.

	<p>
	NOTE: this is more efficient than computing the tranpose(A) explicitly,
	and then multiplying, as the tranpose of A is not explicitly constructed.

	@param A  Matrix: size M x N.
	@param b	Vector: size M.
	@return a new vector of size N.
*/
template <class T>
inline Vector<T> transpose_mult(const Matrix<T>  &A, const Vector<T> &b)
{

#ifdef TNT_BOUNDS_CHECK
    assert(A.num_rows() == b.dim());
#endif

    Subscript M = A.num_cols();
    Subscript N = A.num_rows();
	
    Vector<T> tmp(M);

    for (Subscript i=0; i<M; i++)
    {
        T sum = 0;
        for (Subscript j=0; j<N; j++)
            sum = sum +  A[j][i] * b[j];

        tmp[i] = sum; 
    }

    return tmp;
}

 
/**
		Matrix addition: compute A + B

		@param A	matrix of size M x N.
		@param B	matrix of size M x N.

		@param the sum A+B.
*/
template <class T>
Matrix<T> add(const Matrix<T> &A, const Matrix<T> &B)
{
    Subscript M = A.num_rows();
    Subscript N = A.num_cols();

    assert(M==B.num_rows());
    assert(N==B.num_cols());

    Matrix<T> tmp(M,N);
    Subscript i,j;

    for (i=0; i<M; i++)
        for (j=0; j<N; j++)
            tmp[i][j] = A[i][j] + B[i][j];

    return tmp;
}

/**
		Matrix addition: compute A + B.

		<p>
		NOTE: this is shorthand notation for add(A,B).
		
		@param A	matrix of size M x N.
		@param B	matrix of size M x N.

		@param the sum A+B.
*/
template <class T>
inline Matrix<T> operator+(const Matrix<T> &A, const Matrix<T> &B)
{
	return add(A,B);
}



/**
		Matrix addition, in place : compute A = A + B.


		@param A	matrix of size M x N.
		@param B	matrix of size M x N.

		@param the sum A+B.
*/
template <class T>
Matrix<T>& add_eq(Matrix<T> &A, const Matrix<T> &B)
{
    Subscript M = A.num_rows();
    Subscript N = A.num_cols();

    assert(M==B.num_rows());
    assert(N==B.num_cols());

    Matrix<T> tmp(M,N);
    Subscript i,j;

    for (i=0; i<M; i++)
        for (j=0; j<N; j++)
            tmp[i][j] = A[i][j] + B[i][j];

    return A += tmp;
}

/**
		Matrix addition, in place: compute A = A + B.

		<p>
		NOTE: this is shorthand notation for add_eq(A,B).
		
		@param A	matrix of size M x N.
		@param B	matrix of size M x N.

		@return the sum A+B.
*/
template <class T> 
inline Matrix<T> operator+=(Matrix<T> &A, const Matrix<T> &B)
{
	return add_eq(A,B);
}

/**
		Matrix subtraction : compute A - B.


		@param A	matrix of size M x N.
		@param B	matrix of size M x N.

		@return the result A-B.
*/
template <class T>
Matrix<T> minus(const Matrix <T>& A, const Matrix<T> &B)
{
    Subscript M = A.num_rows();
    Subscript N = A.num_cols();

    assert(M==B.num_rows());
    assert(N==B.num_cols());

    Matrix<T> tmp(M,N);
    Subscript i,j;

    for (i=0; i<M; i++)
        for (j=0; j<N; j++)
            tmp[i][j] = A[i][j] - B[i][j];

    return tmp;

}

/**
		Matrix subtraction : compute A - B.

		<p>
		This is shorthand notation for minus(A,B).

		@param A	matrix of size M x N.
		@param B	matrix of size M x N.

		@return the result A-B.
*/
template <class T>
inline Matrix<T> operator-(const Matrix<T> &A, 
    const Matrix<T> &B)
{
	return minus(A,B);
}


/**

		Matrix element-by-elment multiplication: for each (i,j)
		compute A(i,j) * B(i,j).

	@param A matrix of size M x N.
	@param B matrix of size M x N.
	@return new matrix, where each (i,j) is A(i,j) * B(i,j);

*/
template <class T>
Matrix<T> mult_element(const Matrix<T> &A, const Matrix<T> &B)
{
    Subscript M = A.num_rows();
    Subscript N = A.num_cols();

    assert(M==B.num_rows());
    assert(N==B.num_cols());

    Matrix<T> tmp(M,N);
    Subscript i,j;

    for (i=0; i<M; i++)
        for (j=0; j<N; j++)
            tmp[i][j] = A[i][j] * B[i][j];

    return tmp;
}

/**

		Matrix element-by-elment multiplication, in place: for each (i,j)
		compute A(i,j) = A(i,j) * B(i,j).

	@param A matrix of size M x N.
	@param B matrix of size M x N.
	@return the resultant A, where A(i,j) *= B(i,j);

*/
template <class T>
Matrix<T>&  mult_element_eq(Matrix<T> &A, const Matrix<T> &B)
{
    Subscript M = A.num_rows();
    Subscript N = A.num_cols();

    assert(M==B.num_rows());
    assert(N==B.num_cols());

    Subscript i,j;

    for (i=0; i<M; i++)
        for (j=0; j<N; j++)
            A[i][j] *= B[i][j];

		return A;
}


/**
	Compute Frobenius norm of matrix.  This is the
	square root of the sum of squares of each matrix entry, i.e.
	<pre>
			$$ \\sqrt{ \\sum_{i=1}{N} \\sum_{j=1}{N} A_{i,j}^{2} } $$.
	</pre>.  

	@param A the matrix to compute its Frobeinus norm.
	@return the Frobenius norm of A.
*/
template <class T>
double norm(const Matrix<T> &A)
{
	Subscript M = A.num_rows();
	Subscript N = A.num_cols();

	T sum = 0.0;
	for (int i=1; i<=M; i++)
		for (int j=1; j<=N; j++)
					sum += A(i,j) * A(i,j);
	return sqrt(sum);

}

/**

		Matrix tranpose: return a new matrix B, where B(i,j)
		is A(j,i).

		@param A matrix MxN
		@return new matrix of size N x M, where each (i,j) is 
		A(j,i).
*/
template <class T>
Matrix<T> transpose(const Matrix<T> &A)
{
    Subscript M = A.num_rows();
    Subscript N = A.num_cols();

    Matrix<T> S(N,M);
    Subscript i, j;

    for (i=0; i<M; i++)
        for (j=0; j<N; j++)
            S[j][i] = A[i][j];

    return S;
}







/**

	Solve the triangular system A_u *x=b, where A_u is the upper triangular
		portion (including the diagonal) of A.

	@param A 	a square matrix of size N x N.
	@param b	the right-hand-side (solution vector) of size N.

	@return x, such that A_u * x = b.

*/
template <class T>
Vector<T> upper_triangular_solve(const Matrix<T> &A, const Vector<T> &b) 
{
	int n = A.num_rows() < A.num_cols() ? A.num_rows() : A.num_cols();
	Vector<T> x(b);
	for (int k=n; k >= 1; k--)
	{
		x(k) /= A(k,k);
		for (int i=1; i< k; i++ )
			x(i) -= x(k) * A(i,k);
		}

		return x;
}
		
/**

	Solve the triangular system A_L *x=b, where A_L is the lower triangular
		portion (including the diagonal) of A.

	@param A 	a square matrix of size N x N.
	@param b	the right-hand-side (solution vector) of size N.

	@return x, such that A_u * x = b.

*/
template <class T>
Vector<T> lower_triangular_solve(const Matrix<T> &A, const Vector<T> &b) 
{
	int n = A.num_rows() < A.num_cols() ? A.num_rows() : A.num_cols();
	Vector<T> x(b);
	for (int k=1; k <= n; k++)
	{
		x(k) /= A(k,k);
		for (int i=k+1; i<= n; i++ )
			x(i) -= x(k) * A(i,k);
		}

		return x;
}
		

//}  	/*  namespace TNT::Linear_Algebra; */
} 	/* namespace TNT  */

#endif
// TNT_MATRIX_H
