/*! \file vec.h
 \brief Definition of simple concrete vector class.
 
 Templated numerical vector class based on Template Numerical Toolkit 
 (http://math.nist.gov/tnt) Vector class.

 This class adds the following to the TNT Vector class:
 
 -# Assignment operators (+=, -=, *=, /=)
 -# Errors result in an exception being thrown rather than the program being aborted.

 */

/*

 Simple C++ Numerical Toolkit (SCPPNT)
 http://www.smallwaters.com/software/cpp/scppnt.html
 This release updates original work contributed by 
 Brad Hanson (http://www.b-a-h.com/) in 2001.

 */

// This class is based on the Vector class in:
/*
 *
 * Template Numerical Toolkit (TNT): Linear Algebra Module
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
 * and is in the public domain.  The Template Numerical Toolkit (TNT) is
 * an experimental system.  NIST assumes no responsibility whatsoever for
 * its use by other parties, and makes no guarantees, expressed or implied,
 * about its quality, reliability, or any other characteristic.
 *
 * BETA VERSION INCOMPLETE AND SUBJECT TO CHANGE
 * see http://math.nist.gov/tnt for latest updates.
 *
 */

// Basic TNT  numerical vector (0-based [i] AND 1-based (i) indexing )
//

#ifndef SCPPNT_VEC_H
#define SCPPNT_VEC_H

#include <iterator>

//#ifdef SCPPNT_NO_DIR_PREFIX
#include "scppnt.h"
#include "scppnt_error.h"
//#else
//#include "scppnt/scppnt.h"
//#include "scppnt/scppnt_error.h"
//#endif

#ifdef SCPPNT_BOUNDS_CHECK
#ifdef SCPPNT_NO_DIR_PREFIX
#include "slice_pointer.h"
#else
#include "scppnt/slice_pointer.h"
#endif
#endif

// If SCPPNT_NO_IO is defined then functions that read and write a vector to a file
// and the constructor from a string are not compiled
#ifndef SCPPNT_NO_IO
#include <string>
#include <iostream>
#include <strstream>
#endif
#include <cstdlib> // for NULL
#include <iterator>
#include <algorithm>

namespace SCPPNT
{

  /*!
   \brief Simple concrete vector class for numerical computation.
   
   
   Templated numerical vector class based on Template Numerical Toolkit
   Vector class. Element access through 0-based [i] and 1-based (i)
   indexing, and iterators.

   This class adds the following to the TNT Matrix class:
   
   -# Assignment operators (+=, -=, *=, /=)
   -# Errors result in an exception being thrown rather than the program being aborted.

   */
  template<class T> class Vector
  {

public:

    typedef Subscript size_type; //!< Subscript type
    typedef T value_type; //!< Type of elements stored in vector
    typedef T element_type; //!< Type of elements stored in vector
    typedef T& reference; //!< Reference to type stored in vector
    typedef const T& const_reference;//!< Reference to type stored in constant vector

    // If SCPPNT_BOUNDS_CHECK is defined use slice_pointer_base
    // in order to check for an iterator pointing to an
    // element outside of the vector
#ifdef SCPPNT_BOUNDS_CHECK
    typedef slice_pointer_base<T, T*, T&> pointer;
    typedef slice_pointer_base<T, T*, T&> iterator;
    typedef slice_pointer_base<T, const T*, const T&> const_iterator;
#else
    typedef T* pointer; //!< Pointer to type stored in vector
    typedef T* iterator; //!< Iterator over elements in vector
    typedef const T* const_iterator; //!< Iterator over elements in constant vector
#endif

    //! Returns lower bound of subscript
    Subscript lbound() const;

    // iterators

    iterator begin(); //!< Return iterator pointing to first element of vector.
    iterator end(); //!< Return iterator pointing to one past last element of vector.
    const_iterator begin() const; //!< Return iterator pointing to first element of constant vector.
    const_iterator end() const; //!< Return iterator pointing to one past last element of constant vector.

    //! Destructor
    ~Vector();

    //! Default constructor
    Vector();

    //! Construct vector of given size, do not initialize
    Vector(Subscript N);

    //! Copy constructor
    Vector(const Vector<T> &A);

    //! Constructor which assigns all elements to a particular value.
    Vector(Subscript N, const T& value);

#ifndef SCPPNT_NO_IO
    //! Constructor that reads elements from a string
    Vector(Subscript N, const std::string &s);
#endif

    /* The following member template must be defined within
     the class definition for some compilers (such as Visual C++ 6).
     The definitions of other members are given outside the class definition */

    /*! \brief Create a vector and assign elements from an iterator.
     
     \param begin Iterator pointing to the first element of a sequence
     of values assigned to elements of the vector created.
     
     \param end Iterator pointing to one past last element of a sequence
     of values assigned to elements of the vector created.
     
     */
    template<class IT> Vector(IT begin, IT end) :
      v_(0)
    {
#ifndef BOOST_MSVC
      // "typename" keyword added (ww, 12-2-2007)
      typename std::iterator_traits<IT>::difference_type N = std::distance(begin, end);
#else
      // std::iterator_traits<IT>::difference_type does not compile
      // with Microsoft Visual C++ 6.
      unsigned long N = std::distance(begin, end);
#endif
      initialize(N);

      iterator i = this->begin();
      while (N--)
      {
        *i = *begin;
        ++i;
        ++begin;
      }
    }

    //! Resize vector (current elements are not retained).
    Vector<T>& newsize(Subscript N);

    // assignments
    //

    //! Assignment from vector.
    Vector<T>& operator=(const Vector<T> &A);

    //! Assign all elements to be equal to a scalar.
    Vector<T>& operator=(const T& scalar);

    // assignment operators

    //! Add vector rhs to vector
    Vector<T>& operator+=(const Vector<T> &rhs);

    //! Add scalar to elements of vector
    Vector<T>& operator+=(const T& value);

    //! Subtract vector rhs from vector.
    Vector<T>& operator-=(const Vector<T> &rhs);

    //! Subtract scalar from elements of vector.
    Vector<T>& operator-=(const T& value);

    //! Multiply each element of vector by corresponding element in rhs (element-wise multiply)
    Vector<T>& operator*=(const Vector<T> &rhs);

    //! Multiply each element of vector by scalar.
    Vector<T>& operator*=(const T& value);

    //! Divide each element of vector by corresponding element in rhs (element-wise divide)
    Vector<T>& operator/=(const T& value);

    //! Return number of elements in vector
    Subscript dim() const;

    //! Return number of elements in vector
    Subscript size() const;

    /***** subscripting *****/

    //! One-based subscripting (v(1) returns first element of v)
    inline reference operator()(Subscript i);

    //! One-based subscripting for constant vector
    inline const_reference operator()(Subscript i) const;

    //! Zero-based subscripting (v[0] returns first element of v)
    inline reference operator[](Subscript i);

    //! Zero-based subscripting for constant vector
    inline const_reference operator[](Subscript i) const;

protected:
    T* v_; //!< Pointer to contents of vector.
    T* vm1_; //!< Pointer adjustment for optimizied 1-based indexing.
    Subscript n_; //!< Number of elements in vector.

    // internal helper functions for allocation, initialization, and deallocation

    //! Allocate storage for vector
    void initialize(Subscript N);

    //! Assign elements of vector from array v
    void copy(const T* v);

    //! Set all elements of vector equal to val
    void set(const T& val);

    //! Release storage allocated for vector
    void destroy();

  };

  /***** Definitions of member functions *****/

  template<class T> void Vector<T>::initialize(Subscript N)
  {
    // adjust pointers so that they are 1-offset:
    // v_[] is the internal contiguous array, it is still 0-offset
    //

    if (v_ != NULL)
      throw LogicError(0, "v_ is not NULL", "SCPPNT::Vector::initialize");

    v_ = new T[N];

    vm1_ = v_-1;
    n_ = N;
  }

  template<class T> void Vector<T>::copy(const T* v)
  {
    Subscript N = n_;
    Subscript i;

#ifdef SCPPNT_UNROLL_LOOPS
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

  template<class T> void Vector<T>::set(const T& val)
  {
    Subscript N = n_;
    Subscript i;

#ifdef SCPPNT_UNROLL_LOOPS
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

  template<class T> void Vector<T>::destroy()
  {
    /* do nothing, if no memory has been previously allocated */
    if (v_ == NULL)
      return;

    /* if we are here, then matrix was previously allocated */
    delete [] (v_);

    v_ = NULL;
    vm1_ = NULL;
  }

  // destructor
  template<class T> Vector<T>::~Vector()
  {
    destroy();
  }

  // constructors

  template<class T> Vector<T>::Vector() :
    v_(0), vm1_(0), n_(0)
  {
  }

  template<class T> Vector<T>::Vector(Subscript N) :
    v_(0), vm1_(0)
  {
    initialize(N);
  }

  template<class T> Vector<T>::Vector(const Vector<T> &A) :
    v_(0), vm1_(0), n_(0)
  {
    initialize(A.n_);
    copy(A.v_);
  }

  template<class T> Vector<T>::Vector(Subscript N, const T& value) :
    v_(0), vm1_(0)
  {
    initialize(N);
    set(value);
  }

  // Initialize from an iterator
  // This definition is contained within the class definition above
  // so that it will compile using Visual C++ 6.
  // template <class T>
  // template <class IT>
  // Vector<T>::Vector(IT begin, IT end)

#ifndef SCPPNT_NO_IO
  /*!
   Constructor that reads elements from a string
   
   \param N
   Number of elements in vector.
   \param s
   String containing initial elements of matrix in
   row-major order separated by white space.
   */
  template<class T> Vector<T>::Vector(Subscript N, const std::string &s) :
    v_(0), vm1_(0), n_(0)
  {
    initialize(N);
    std::istrstream ins(s.c_str());

    Subscript i;

    for (i=0; i<N; i++)
      ins >> v_[i];
  }
#endif

  // methods
  // 
  template<class T> Vector<T>& Vector<T>::newsize(Subscript N)
  {
    if (n_ == N)
      return *this;

    destroy();
    initialize(N);

    return *this;
  }

  // assignments
  //
  template<class T> Vector<T>& Vector<T>::operator=(const Vector<T> &A)
  {
    if (v_ == A.v_)
    return *this;

    if (n_ == A.n_) // no need to re-alloc
    copy(A.v_);

    else
    {
      destroy();
      initialize(A.n_);
      copy(A.v_);
    }

    return *this;
  }

  // assign all elements to equal a scalar
  template <class T>
  Vector<T>& Vector<T>::operator=(const T& scalar)
  {
    set(scalar);
    return *this;
  }

  // vector operator +=
  // Replace each element of vector with the element plus
  // the corresponding element of rhs.
  template <class T>
  Vector<T>& Vector<T>::operator+=(const Vector<T> &rhs)
  {

    if (n_ != rhs.dim()) throw BadDimension("SCPPNT::Vector::operator+=");

    const_iterator rhsiter = rhs.begin();
    iterator end = this->end();
    for (iterator p = begin(); p != end; ++p, ++rhsiter)
    {
      *p += *rhsiter;
    }

    return *this;

  }

  // scalar operator += (add a scalar to each element)
  template <class T>
  Vector<T>& Vector<T>::operator+=(const T& value)
  {
    iterator end = this->end();
    for (iterator p = begin(); p != end; ++p)
    {
      *p += value;
    }

    return *this;

  }

  // vector operator -=
  // Replace each element of vector with the element minus
  // the corresponding element of rhs.
  template <class T>
  Vector<T>& Vector<T>::operator-=(const Vector<T> &rhs)
  {

    if (n_ != rhs.dim()) throw BadDimension("SCPPNT::Vector::operator-=");

    const_iterator rhsiter = rhs.begin();
    iterator end = this->end();
    for (iterator p = begin(); p != end; ++p, ++rhsiter)
    {
      *p -= *rhsiter;
    }

    return *this;

  }

  // scalar operator -= (subtract a scalar from each element)
  template <class T>
  Vector<T>& Vector<T>::operator-=(const T& value)
  {
    iterator end = this->end();
    for (iterator p = begin(); p != end; ++p)
    {
      *p -= value;
    }

    return *this;

  }

  // vector operator *=
  // Replace each element of vector with the product of the
  // element and the corresponding element of rhs.
  template <class T>
  Vector<T>& Vector<T>::operator*=(const Vector<T> &rhs)
  {

    if (n_ != rhs.dim()) throw BadDimension("SCPPNT::Vector::operator*=");

    const_iterator rhsiter = rhs.begin();
    iterator end = this->end();
    for (iterator p = begin(); p != end; ++p, ++rhsiter)
    {
      *p *= *rhsiter;
    }

    return *this;

  }

  // scalar operator *= (multiply each element by scalar)
  template <class T>
  Vector<T>& Vector<T>::operator*=(const T& value)
  {
    iterator end = this->end();
    for (iterator p = begin(); p != end; ++p)
    {
      *p *= value;
    }

    return *this;

  }

  // scalar operator /= (divide each element by scalar)
  template <class T>
  Vector<T>& Vector<T>::operator/=(const T& value)
  {
    iterator end = this->end();
    for (iterator p = begin(); p != end; ++p)
    {
      *p /= value;
    }

    return *this;

  }

  template <class T>
  inline Subscript Vector<T>::lbound() const
  { return 1;}

  template <class T>
  inline Subscript Vector<T>::dim() const
  {
    return n_;
  }

  template <class T>
  inline Subscript Vector<T>::size() const
  {
    return n_;
  }

  // subscripting
  template <class T>
  inline typename Vector<T>::reference Vector<T>::operator()(Subscript i)
  {
#ifdef SCPPNT_BOUNDS_CHECK
    if (i < 1 || i > n_) throw BoundsError("SCPPNT::Vector::operator()");
#endif
    return vm1_[i];
  }

  template <class T>
  inline typename Vector<T>::const_reference Vector<T>::operator() (Subscript i) const
  {
#ifdef SCPPNT_BOUNDS_CHECK
    if (i < 1 || i > n_) throw BoundsError("SCPPNT::Vector::operator() const");
#endif
    return vm1_[i];
  }

  template <class T>
  inline typename Vector<T>::reference Vector<T>::operator[](Subscript i)
  {
#ifdef SCPPNT_BOUNDS_CHECK
    if (i < 0 || i >= n_) throw BoundsError("SCPPNT::Vector::operator[]");
#endif
    return v_[i];
  }

  template <class T>
  inline typename Vector<T>::const_reference Vector<T>::operator[](Subscript i) const
  {
#ifdef SCPPNT_BOUNDS_CHECK
    if (i < 0 || i >= n_) throw BoundsError("SCPPNT::Vector::operator[] const");
#endif
    return v_[i];
  }

  // iterators
  template <class T>
  inline typename Vector<T>::iterator Vector<T>::begin()
  {
#ifdef SCPPNT_BOUNDS_CHECK
    return iterator(v_, 1, n_);
#else
    return v_;
#endif
  }

  template <class T>
  inline typename Vector<T>::iterator Vector<T>::end()
  {
#ifdef SCPPNT_BOUNDS_CHECK
    iterator p(v_, 1, n_);
    return p + n_;
#else
    return v_ + n_;
#endif
  }

  template <class T>
  inline typename Vector<T>::const_iterator Vector<T>::begin() const
  {
#ifdef SCPPNT_BOUNDS_CHECK
    return const_iterator(v_, 1, n_);
#else
    return v_;
#endif
  }

  template <class T>
  inline typename Vector<T>::const_iterator Vector<T>::end() const
  {
#ifdef SCPPNT_BOUNDS_CHECK
    const const_iterator p(v_, 1, n_);
    return p + n_;
#else
    return v_ + n_;
#endif
  }

#ifndef SCPPNT_NO_IO
  /* ***************************  I/O  ********************************/

  /*!
   \brief Write a vector to an output stream.
   
   Writes the number of elements in the vector, followed
   by a new line. The elements of the vector are then
   written on consecutive lines.
   */
  template <class T>
  std::ostream& operator<<(std::ostream &s, const Vector<T> &A)
  {
    Subscript N=A.dim();

    s << N << std::endl;

    for (Subscript i=0; i<N; i++)
    s << A[i] << " " << std::endl;
    s << std::endl;

    return s;
  }

  /*!
   \brief Read elements of a vector from an input stream.
   
   Reads the number of elements in the vector followed
   by the elements of the vector. All values read
   should be separated by white space.
   */
  template <class T>
  std::istream & operator>>(std::istream &s, Vector<T> &A)
  {

    Subscript N;

    s >> N;

    A.newsize(N);

    for (Subscript i=0; i<N; i++)
    s >> A[i];

    return s;
  }
#endif

  // *******************[ basic vector algorithms ]***************************

  //! Add a vector to a vector and return a new vector containing the sum.
  template <class T>
  Vector<T> operator+(const Vector<T> &lhs,
      const Vector<T> &rhs)
  {
    return Vector<T>(lhs) += rhs;
  }

  //! Add add a scalar to each element of a vector and return a new vector containing the sum.
  template <class T>
  Vector<T> operator+(const Vector<T> &lhs,
      const T value)
  {
    return Vector<T>(lhs) += value;
  }

  //! Subtract the vector rhs from the vector lhs and return a new vector containing the difference.
  template <class T>
  Vector<T> operator-(const Vector<T> &lhs,
      const Vector<T> &rhs)
  {
    return Vector<T>(lhs) -= rhs;
  }

  //! Subtract scalar from each element of a vector a return a new vector containing the difference.
  template <class T>
  Vector<T> operator-(const Vector<T> &lhs,
      const T value)
  {
    return Vector<T>(lhs) -= value;
  }

  //! Multiply corresponding elements of two vectors and return a new vector containing the product.
  template <class T>
  Vector<T> operator*(const Vector<T> &lhs,
      const Vector<T> &rhs)
  {
    return Vector<T>(lhs) *= rhs;
  }

  //! Multiply each element of a vector by a scalar and return a new vector containing the product.
  template <class T>
  Vector<T> operator*(const Vector<T> &lhs,
      const T value)
  {
    return Vector<T>(lhs) *= value;
  }

  /***** Dot products *****/

  /*! \brief Dot product of values pointed to by two iterators

   Returns sum of products of corresponding elements pointed to by
   the iterators i1 and i2.
   
   \param N Number of elements in dot product
   \param i1 Iterator to first sequence of elements to use for product.
   \param i2 Iterator to second sequence of elements to use for product.

   */
  template <class IT1, class IT2>
  typename std::iterator_traits<IT1>::value_type dot_prod(Subscript N, IT1 i1, IT2 i2)
  {

    typename std::iterator_traits<IT1>::value_type sum = 0;
    while(N--)
    {
      sum += *i1 * *i2;
      ++i1;
      ++i2;
    }

    return sum;
  }

  //! Dot product of two vectors.
  template <class T>
  inline T dot_prod(const Vector<T> &A, const Vector<T> &B)
  {
    Subscript N = A.dim();

    if (N != B.dim()) throw BadDimension("SCPPNT::dot_prod()");

    return dot_prod(N, A.begin(), B.begin());
  }

} /* namespace SCPPNT */

#endif
// SCPPNT_VEC_H
