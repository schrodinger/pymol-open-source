#ifndef TNT_SPARSE_VECTOR_H
#define TNT_SPARSE_VECTOR_H

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



#include <vector>
#include "tnt_vector.h"

namespace TNT
{
//namespace Linear_Algebra
//{



template <class T, class Integer>
class Sparse_Vector_Element
{
		private:

				T val_;
				Integer index_;

		public:

				Sparse_Vector_Element(const T& a, const Integer &i) : 
						val_(a), index_(i) {}
				

				const T& value() const { return val_; } 
				Integer	index() const { return index_; } 

				T& value() { return val_; }
				Integer& index() { return index_; }

					

};


/**

		Sparse Vector.

	<p>
	Index values begin at 0.  Thus S[3] = 6.5 means that the fourth, not
	the third element, is set to 6.5.

	<p>
	S.value(0) is the first nonzero values in S, and S.index(0) is
	the index (0-based) of the first in S.

	

*/
template <class T>
class Sparse_Vector
{
	/* sparse vector indices are 0-based internally. */

  private:

		// (value, index) pairs
		//
		std::vector< Sparse_Vector_Element<T, Subscript> > s_;    

  	int dim_;        			// dimension
		int num_nonzeros_; 		// number of nonzeros




public:
		typedef typename 
				std::vector< Sparse_Vector_Element<T, Subscript> >::const_iterator 
																															const_iterator;


		const_iterator begin() const { return s_.begin(); }
		const_iterator end() const { return s_.end(); }

		inline const T& value(Subscript i) const { return s_[i].value(); }
		inline Subscript index(Subscript i) const {return  s_[i].index(); }
	  
    inline T dot_product(const Vector<T>& x) const
    {
        T sum(0);

        for ( const_iterator p = begin(); p < end(); p++)
        {
						sum += p->value() * x[p->index()];
        }
       return sum;
   	}

	Sparse_Vector() : s_(), dim_(0), num_nonzeros_(0) {}
	Sparse_Vector(Subscript N): dim_(N), num_nonzeros_(0) {}

			

	Sparse_Vector(Subscript N, Subscript nz, const T* Val, const Subscript *I):
					dim_(N),
					num_nonzeros_(0)
					{
						insert(nz, Val, I);	
					}
   


	void   insert(const T& val, Subscript i)
	{

			s_.push_back( Sparse_Vector_Element<T, Subscript>(val,i) );
			num_nonzeros_++;
	}

	void   insert(Subscript nz, const T* Val, const Subscript *I)
	{

			for (int count=0; count<nz; count++)
			{
				insert(Val[count], I[count]);	
			}
	}


	void insert_base_one(const T& val, Subscript i)
	{
		insert(val, i-1);
	}

	void   insert_base_one(Subscript nz, const T* Val, const Subscript *I)
	{
			for (int count=0; count<nz; count++)
			{
				insert(Val[count], I[count]-1);	
			}
	}
		


  inline   int    dim() const 	{return dim_;}
  int          num_nonzeros() const {return num_nonzeros_;}




	inline double norm() const
	{
		T sum(0.0);

		for (const_iterator p = begin(); p < end(); p++)
		{
			sum += p->value() * p->value();
		}

		return sqrt(sum);
	}

	std::ostream & print(std::ostream &s) const
	{
		for (typename Sparse_Vector<T>::const_iterator p = begin(); 
						p < end(); p++)
		{
				s << "( " << p->value() << ", " << p->index() << " ) \n";
		}
		return s;
	}

	std::ostream & print_base_one(std::ostream &s) const
	{
		for (typename Sparse_Vector<T>::const_iterator p = begin(); 
						p < end(); p++)
		{
				s << "( " << p->value() << ", " << p->index()+1 << " ) \n";
		}
		return s;
	}

			
};



template <class T>
inline T dot_product(const Sparse_Vector<T> &s, const Vector<T> &x)
{
	return s.dot_product(x);
}

template <class T>
inline T dot_product(const Vector<T> &x, const Sparse_Vector<T> &s)
{
	return s.dot_product(x);
}

template <class T>
inline T operator*(const Vector<T> &x, const Sparse_Vector<T> &s)
{
	return dot_product(x,s);
}

template <class T>
inline T operator*(const Sparse_Vector<T> &s, const Vector<T> &x)
{
	return dot_product(x,s);
}


template <class T>
inline double norm( const Sparse_Vector<T> & s)
{
	return s.norm();
}


template <class T>
std::ostream& operator<<(std::ostream &s, const Sparse_Vector<T> &A)
{

		return	A.print(s);
}


//}  /* namspace TNT::Linear_Algebra */
}	 /* namspace TNT */

#endif
// TNT_SPARSE_VECTOR_H
