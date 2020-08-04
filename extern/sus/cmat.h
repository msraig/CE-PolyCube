/*! \file cmat.h
 \brief Definition of a simple concrete dense matrix class.
 
 Templated numerical matrix class based on Template Numerical Toolkit 
 (http://math.nist.gov/tnt) Matrix class.
 This class adds the following to the TNT Matrix class:
 
 -# Row, column, and diagonal iterators.
 -# Assignment operators (+=, -=, *=, /=)
 -# Errors result in an exception being thrown rather than the program being aborted.

 */

/*
 
 Simple C++ Numerical Toolkit (SCPPNT)
 http://www.smallwaters.com/software/cpp/scppnt.html
 This release updates original work contributed by 
 Brad Hanson (http://www.b-a-h.com/) in 2001.

 */

// Modified from cmat.h in:
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

#ifndef SCPPNT_CMAT_H
#define SCPPNT_CMAT_H

#ifdef SCPPNT_NO_DIR_PREFIX
#include "scppnt.h"
#include "slice_pointer.h"
#include "scppnt_error.h"
#include "matop.h"
#include "vec.h"
#else
#include "scppnt/scppnt.h"
#include "scppnt/slice_pointer.h"
#include "scppnt/scppnt_error.h"
#include "scppnt/matop.h"
#include "scppnt/vec.h"
#endif

// If SCPPNT_NO_IO is defined then functions that read and write a matrix to a file
// and the constructor from a string (which uses strstream) are not compiled
#ifndef SCPPNT_NO_IO
#include <string>
#include <iostream>
#include <strstream>
#endif

// If SCPPNT_USE_REGIONS is defined then a function call operator is
// defined for SCPPNT::Matrix which returns a region. This allows
// a simpler notation for creating regions from matrices as compared to
// using Region2D constructors.
#ifdef SCPPNT_USE_REGIONS
#ifdef SCPPNT_NO_DIR_PREFIX
#include "region2d.h"
#else
#include "scppnt/region2d.h"
#endif
#endif

namespace SCPPNT
{

  /*!
   \brief Simple concrete dense matrix class for numerical computation.
   
   
   Templated numerical matrix class based on Template Numerical Toolkit
   Matrix class. Row-oriented matrix with element access through
   0-based [i][j] and 1-based (i,j) indexing, and row, column, and
   diagonal iterators. The matrix is stored in one contiguous piece of
   memory. Besides memory being allocated for the contents of the
   matrix, memory is also allocated for row iterators (one iterator
   for each row), column iterators (one iterator for each column),
   and column iterators over constant elements (one iterator for
   each column). If SCPPNT_BOUNDS_CHECK is defined then memory
   is allocated for row iterators over constant elements (one
   for each row).
   
   This class adds the following to the TNT Matrix class:
   
   -# Row, column, and diagonal iterators.
   -# Assignment operators (+=, -=, *=, /=)
   -# Errors result in an exception being thrown rather than the program being aborted.
   
   */
  template<class T> class Matrix
  {

public:

    typedef Subscript size_type; //!< Subscript type
    typedef T value_type; //!< Type of elements stored in matrix
    typedef T element_type; //!< Type of elements stored in matrix
    typedef T* pointer; //!< Pointer to type stored in matrix
    typedef T& reference; //!< Reference to type stored in matrix
    typedef const T& const_reference;//!< Reference to type stored in constant matrix

#ifdef SCPPNT_BOUNDS_CHECK
    typedef slice_pointer_base<T, T*, T&> iterator;
    typedef slice_pointer_base<T, const T*, const T&> const_iterator;
#else
    typedef T* iterator; //!< Iterator over elements in matrix
    typedef const T* const_iterator; //!< Iterator over elements in constant matrix
#endif

    /*  
     row_iterator - Iterator over elements in a particular row.
     column_iterator - Iterator over elements in a particular column.
     
     The reason for having different types for iterators over elements
     of rows and columns is to allow the most efficient representation
     to be used in each case. For example, if the elements of the
     matrix are stored so that consecutive elements in each row are
     in adjacent memory locations then a row_iterator can be a T*,
     whereas a column_iterator must be a slice_pointer_base<T, T*, T&>.
     
     */
#ifdef SCPPNT_BOUNDS_CHECK
    typedef slice_pointer_base<T, T*, T&> row_iterator;
    typedef slice_pointer_base<T, const T*, const T&> const_row_iterator;
#else
    //! Iterator over elements of a row
    typedef T* row_iterator;

    //! Iterator over constant elements of a row
    typedef const T* const_row_iterator;
#endif
    //! Iterator over elements of a column
    typedef slice_pointer_base<T, T*, T&> column_iterator;

    //! Iterator over constant elements of a column
    typedef slice_pointer_base<T, const T*, const T&> const_column_iterator;

    /*
     rows_iterator - Iterator over row_iterators for rows of matrix.
     columns_iterator - Iterator over column_iterators for columns of matrix.
     */
#ifdef SCPPNT_BOUNDS_CHECK
    typedef slice_pointer_base<row_iterator, row_iterator*, row_iterator&> rows_iterator;
    typedef slice_pointer_base<column_iterator, column_iterator*, column_iterator&>
        columns_iterator;
    typedef slice_pointer_base<const_row_iterator, const_row_iterator*, const_row_iterator&>
        const_rows_iterator;
    typedef slice_pointer_base<const_column_iterator, const_column_iterator*, const_column_iterator&>
        const_columns_iterator;
#else
    //! Iterator over row iterators (points to row_iterator object for a row)
    typedef row_iterator* rows_iterator;

    //! Iterator over row iterators (points to const_row_iterator object for a row)
    typedef const_row_iterator* const_rows_iterator;

    //! Iterator over column iterators (points to column_iterator object for a column)
    typedef column_iterator* columns_iterator;

    //! Iterator over column iterators (points to const_column_iterator object for a column)
    typedef const_column_iterator* const_columns_iterator;
#endif

    //! Iterator over elements of a matrix diagonal
    typedef slice_pointer_base<T, T*, T&> diag_iterator;

    //! Iterator over constant elements of a matrix diagonal
    typedef slice_pointer_base<T, const T*, const T&> const_diag_iterator;

    //! Returns lower bound of subscript
    Subscript lbound() const;

    //! Return total number of elements in matrix
    Subscript size() const;

    //! Return number of rows if argument is 1, number of columns if argument is 2
    Subscript dim(Subscript d) const;

    //! Return number of rows in matrix
    Subscript num_rows() const;

    //! Return number of columns in matrix
    Subscript num_columns() const;

    /***** constructors *****/

    //! Default constructor
    Matrix();

    //! Copy constructor
    Matrix(const Matrix<T> &A);

    //! Construct matrix of a particular size, but do not initialize elements to any value
    Matrix(Subscript M, Subscript N);

    //! Construct and assign all elements to a particular value
    Matrix(Subscript M, Subscript N, const T value);

#ifndef SCPPNT_NO_IO
    //! Constructor that reads elements from a string
    Matrix(Subscript M, Subscript N, const std::string &s);
#endif

    /* The following member template must be defined within
     the class definition for some compilers (such as Visual C++ 6).
     The definitions of other members are given outside the class definition */

    /*! \brief Construct from an iterator giving initial elements in row-major order
     
     \param begin Iterator pointing to the first element of a sequence
     of values assigned to elements of the matrix created in row-major order
     (elements in first row followed by elements in second row, etc.)
     
     \param number_rows Number of rows in matrix created.
     
     \param number_columns Number of columns in matrix created.
     
     The iterator argument preceeds the number of rows and columns in
     this constructor in contrast to other constructors where the
     number of rows and columns are the first two arguments. If the
     number of rows and columns were the first two arguments then
     this template constructor can mistakenly be used
     instead of the Matrix(Subscript M, Subscript N, const std::string &s)
     constructor when constructing a matrix from a string.
     */
    template<class IT> Matrix(IT begin, Subscript number_rows, Subscript number_columns) :
      v_(0)
    {

      initialize(number_rows, number_columns);

      iterator i = this->begin();
      Subscript n = number_rows * number_columns;
      while (n--)
      {
        *i = *begin;
        ++begin;
        ++i;
      }
    }

    //! Destructor
    ~Matrix();

    //! Resize matrix (old elements are distroyed)
    Matrix<T>& newsize(Subscript M, Subscript N);

    // assignments
    //

    //! Matrix assignment
    Matrix<T>& operator=(const Matrix<T> &A);

    //! Scalar assignment
    Matrix<T>& operator=(const T& scalar);

    // subscripting
    //

    //! Return iterator to elements in row i+1 (0-offset)
    row_iterator operator[](Subscript i);

    //! Return iterator to constant elements in row i+1 (0-offset)
    const_row_iterator operator[](Subscript i) const;

    //! Return i-th (1-offset) element of matrix when treated as a vector (row major ordering)
    reference operator()(Subscript i);

    //! Return i-th (1-offset) constant element of matrix when treated as a vector (row major ordering)
    const_reference operator()(Subscript i) const;

    //! Return element in row i and column j, where i and j are 1-offset indices
    reference operator()(Subscript i, Subscript j);

    //! Return constant element in row i and column j, where i and j are 1-offset indices
    const_reference operator()(Subscript i, Subscript j) const;

    /**** Iterators ****/

    //! Return iterator pointing to row iterator for first row (iterator over row iterators)
    rows_iterator begin_rows();

    //! Return iterator pointing to row iterator for first row (constant version)
    const_rows_iterator begin_rows() const;

    //! Return iterator pointing to column iterator for first column (iterator over column iterators)
    columns_iterator begin_columns();

    //! Return iterator pointing to column iterator for first column (constant version)
    const_columns_iterator begin_columns() const;

    //! Return iterator pointing to row iterator for one past last row (iterator over row iterators)
    rows_iterator end_rows();

    //! Return iterator pointing to row iterator for one past last row (constant version)
    const_rows_iterator end_rows() const;

    //! Return iterator pointing to column iterator for one past last column (iterator over column iterators)
    columns_iterator end_columns();

    //! Return iterator pointing to column iterator for one past last column (constant version)
    const_columns_iterator end_columns() const;

    // Return row or column iterator for particular row or column.
    // Index is one-offset (first row or column is 1).
    // begin_row(N+1) returns a row_iterator to one past the last row, where there are N rows
    // begin_column(M+1) returns a column_iterator to one past the last column, where there are M columns
    //

    //! Return iterator pointing to first element in row 'index' (1-offset)
    row_iterator begin_row(Subscript index);

    //! Return iterator pointing to first element in row 'index' (1-offset)
    const_row_iterator begin_row(Subscript index) const;

    //! Return iterator pointing to first element in column 'index' (1-offset)
    column_iterator begin_column(Subscript index);

    //! Return iterator pointing to first element in column 'index' (1-offset)
    const_column_iterator begin_column(Subscript index) const;

    //! Return iterator pointing to one past last element in row 'index' (1-offset)
    row_iterator end_row(Subscript index);

    //! Return iterator pointing to one past last element in row 'index' (1-offset)
    const_row_iterator end_row(Subscript index) const;

    //! Return iterator pointing to one past last element in column 'index' (1-offset)
    column_iterator end_column(Subscript index);

    //! Return iterator pointing to one past last element in column 'index' (1-offset)
    const_column_iterator end_column(Subscript index) const;

    //! Returns iterator pointing to first element of matrix (consecutive element accessed by row)
    iterator begin();

    //! Returns iterator pointing to first element of matrix (consecutive element accessed by row)
    const_iterator begin() const;

    //! Iterator pointing to one past last element of matrix
    iterator end();

    //! Iterator pointing to one past last element of matrix
    const_iterator end() const;

    //! Returns iterator pointing to first element of matrix diagonal
    diag_iterator begin_diagonal(Subscript row, Subscript column);

    //! Returns iterator pointing to first element of matrix diagonal
    const_diag_iterator begin_diagonal(Subscript row, Subscript column) const;

    //! Returns iterator pointing to one past last element of matrix diagonal
    diag_iterator end_diagonal(Subscript row, Subscript column);

    //! Returns iterator pointing to one past last element of matrix diagonal
    const_diag_iterator end_diagonal(Subscript row, Subscript column) const;

#ifndef SCPPNT_BOUNDS_CHECK
    /* For the following conversion is only valid when SCPPNT_BOUNDS_CHECK is not defined.
     When SCPPNT_BOUNDS_CHECK is not defined the
     member function rows() returns the same thing as this conversion operator. */
    //! convert to T**
    operator T**()
    {
      return row_;
    }

    //! convert to const T**
    operator const T**() const
    {
      return (const T**) row_;
    }
#endif

    /*** assignment operators ***/

    /* The following three member templates must be defined within
     the class definition for some compilers (such as Visual C++ 6).
     The definitions of other members are given outside the class definition */

    //! Add matrix rhs to matrix
    template<class MAT> Matrix<T>& operator+=(const MAT &rhs)
    {
      matadd(*this, rhs);
      return *this;
    }

    //! Subtract matrix rhs from matrix
    template<class MAT> Matrix<T>& operator-=(const MAT &rhs)
    {
      matsub(*this, rhs);
      return *this;
    }

    //! Matrix multiply of matrix times rhs (rhs must be square)
    template<class MAT> Matrix<T>& operator*=(const MAT &rhs)
    {
      matmult_assign(*this, rhs);
      return *this;
    }

    //! Add scalar to all elements of the matrix
    Matrix<T>& operator+=(const T& value);

    //! Subtract scalar from all elements of the matrix
    Matrix<T>& operator-=(const T& value);

    //! Multiply each element of matrix by scalar value
    Matrix<T>& operator*=(const T& value);

    //! Divide each element of matrix by scalar value
    Matrix<T>& operator/=(const T& value);

    // The following allows regions to be created from matrices
    // using a simpler notation than a Region2D constructor
#ifdef SCPPNT_USE_REGIONS

    //! Type representing a rectangular region of a matrix
    typedef Region2D<Matrix<T> > Region;

    //! Return a matrix region using index ranges I and J
    Region operator()(const Index1D &I, const Index1D &J)
    {
      return Region(*this, I,J);
    }

    //! Type representing a rectangular region of a constant matrix
    typedef const_Region2D< Matrix<T> > const_Region;

    //! Return a constant matrix region using index ranges I and J
    const_Region operator()(const Index1D &I, const Index1D &J) const
    {
      return const_Region(*this, I,J);
    }

#endif

private:

    Subscript m_; //!< number of rows
    Subscript n_; //!< number of columns
    Subscript mn_; //!< total size
    T* v_; //!< points to contents of matrix
    T* vm1_; //!< points to the same data as v_, but is 1-based 
    row_iterator *row_; //!< points to array of iterators to rows of matrix
    T** rowm1_; //!< points to array of pointers to rows of matrix (1-offset)
    column_iterator *column_; //!< points to array of iterators to columns of matrix
    const_column_iterator *const_column_; //!< points to array of iterators to const columns of matrix

#ifdef SCPPNT_BOUNDS_CHECK
    const_row_iterator *const_row_; //!< points to array of iterators to const rows of matrix
#endif

    // internal helper functions for allocation, initialization, and deallocation

    //! Allocate storage for matrix with M rows and N columns
    void initialize(Subscript M, Subscript N);

    //! Assign elements of matrix from array v (row major storage)
    void copy(const T* v);

    //! Set all elements of matrix equal to val
    void set(const T& val);

    //! Release storage allocated for matrix
    void destroy();

    //! Return the number of elements in the diagonal beginning at (row, column)
    Subscript diagonal_size(Subscript row, Subscript column) const;

  }; // Matrix class

  template<class T> inline Subscript Matrix<T>::lbound() const
  {
    return 1;
  }

  // assignments
  //
  template<class T> Matrix<T>& Matrix<T>::operator=(const Matrix<T> &A)
  {
    if (v_ == A.v_)
    return *this;

    if (m_ == A.m_ && n_ == A.n_) // no need to re-alloc
    copy(A.v_);
    else
    {
      destroy();
      initialize(A.m_, A.n_);
      copy(A.v_);
    }

    return *this;
  }

  template <class T>
  Matrix<T>& Matrix<T>::operator=(const T& scalar)
  {
    set(scalar);
    return *this;
  }

  template <class T>
  inline Subscript Matrix<T>::size() const
  {
    return mn_;
  }

  template <class T>
  Subscript Matrix<T>::dim(Subscript d) const
  {
    if (d < 1 || d > 2) throw InvalidArgument("argument must be 1 or 2", "SCPPNT::Matrix::dim");
    return (d==1) ? m_ : n_;
  }

  template <class T>
  inline Subscript Matrix<T>::num_rows() const
  { return m_;}
  template <class T>
  inline Subscript Matrix<T>::num_columns() const
  { return n_;}

  template <class T>
  void Matrix<T>::initialize(Subscript M, Subscript N)
  {

    if (v_ != 0) throw LogicError(0, "v_ is not NULL", "SCPPNT::Matrix::initialize");

    mn_ = M*N;
    m_ = M;
    n_ = N;

    /* Allocate space for contents of matrix */
    v_ = new T[mn_];
    vm1_ = v_ - 1;

    /* Allocate space for pointers to rows of matrix including space
     for pointer to one past the last row. */
    rowm1_ = new T*[M+1];
    row_ = new row_iterator[M+1];
#ifdef SCPPNT_BOUNDS_CHECK
    const_row_ = new const_row_iterator[M+1];
#endif

    /* assign pointers to elements of row_ and rowm1_ */
    T* p = v_;
    Subscript i;
    for (i=0; i<=M; i++)
    {
#ifdef SCPPNT_BOUNDS_CHECK
      row_[i].set(p, 1, N);
      const_row_[i].set(p, 1, N);
#else
      row_[i] = p;
#endif

      rowm1_[i] = p-1;
      p += N;
    }

    rowm1_ --; // compensate for 1-based offset

    /* Allocate space for pointers to columns of matrix including
     pointer to one past the last column. */
    column_ = new column_iterator[N+1];
    const_column_ = new const_column_iterator[N+1];

    /* assign pointers to elements of column_ */
    p = v_;
    for (i=0; i<=N; ++i)
    {
      column_[i].set(p, N, M);
      const_column_[i].set(p, N, M);
      ++p;
    }

  }

  template <class T>
  void Matrix<T>::destroy()
  {
    /* do nothing, if no memory has been previously allocated */
    if (v_ == 0) return;

    /* if we are here, then matrix was previously allocated */
    delete [] (v_);
    v_ = 0;

    if (row_ != 0) delete [] (row_);
    if (column_ != 0) delete [] (column_);

#ifdef SCPPNT_BOUNDS_CHECK
    if (const_row_ != 0) delete [] (const_row_);
#endif
    if (const_column_ != 0) delete [] (const_column_);

    if (rowm1_ != 0)
    {
      rowm1_ ++; /* return rowm1_ back to original value */
      if (rowm1_ != 0 ) delete [] (rowm1_);
    }

  }

  template <class T>
  Matrix<T>& Matrix<T>::newsize(Subscript M, Subscript N)
  {
    if (num_rows() == M && num_columns() == N)
    return *this;

    destroy();
    initialize(M,N);

    return *this;
  }

  template <class T>
  void Matrix<T>::copy(const T* v)
  {
    Subscript N = m_ * n_;
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

  template <class T>
  void Matrix<T>::set(const T& val)
  {
    Subscript N = m_ * n_;
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
    iterator pi = begin();
    for (i=N; i--; ++pi)
    *pi = val;

#endif      
  }

  // Constructors

  template <class T>
  Matrix<T>::Matrix() : m_(0), n_(0), mn_(0), v_(0), row_(0), vm1_(0), rowm1_(0), column_(0), const_column_(0)
  {
#ifdef SCPPNT_BOUNDS_CHECK
    const_row_ = 0;
#endif
  }

  template <class T>
  Matrix<T>::Matrix(const Matrix<T> &A) : v_(0)
  {
    initialize(A.m_, A.n_);
    copy(A.v_);
  }

  template <class T>
  Matrix<T>::Matrix(Subscript M, Subscript N) : v_(0)
  {
    initialize(M,N);
  }

  template <class T>
  Matrix<T>::Matrix(Subscript M, Subscript N, const T value) : v_(0)
  {
    initialize(M,N);
    set(value);
  }

#ifndef SCPPNT_NO_IO
  /*!
   Constructor that reads elements from a string
   
   \param M
   Number of rows.
   \param N
   Number of columns.
   \param s
   String containing initial elements of matrix in
   row-major order separated by white space.
   */
  template <class T>
  Matrix<T>::Matrix(Subscript M, Subscript N, const std::string &s) : v_(0)
  {
    initialize(M,N);
    std::istrstream ins(s.c_str());

    Subscript i, j;

    for (i=0; i<M; i++)
    for (j=0; j<N; j++)
    ins >> row_[i][j];
  }
#endif

  // Initialize from an iterator
  // This definition is contained within the class definition above
  // so that it will compile using Visual C++ 6.
  // template <class T>
  // template <class IT>
  // Matrix<T>::Matrix(Subscript M, Subscript N, IT v)


  // destructor
  //
  template <class T>
  Matrix<T>::~Matrix()
  {
    destroy();
  }

  /* *********************** Subscripting ****************************/
  template <class T>
  inline typename Matrix<T>::row_iterator Matrix<T>::operator[](Subscript i)
  {
#ifdef SCPPNT_BOUNDS_CHECK
    if (i < 0 || i >= m_) throw BoundsError("SCPPNT::Matrix::operator[]");
#endif
    return row_[i];
  }

  template <class T>
  inline typename Matrix<T>::const_row_iterator Matrix<T>::operator[](Subscript i) const
  {
#ifdef SCPPNT_BOUNDS_CHECK
    if (i < 0 || i >= m_) throw BoundsError("SCPPNT::Matrix::operator[] const");
#endif
    return row_[i];
  }

  template <class T>
  inline typename Matrix<T>::reference Matrix<T>::operator()(Subscript i)
  {
#ifdef SCPPNT_BOUNDS_CHECK
    if (i < 1 || i > mn_) throw BoundsError("SCPPNT::Matrix::operator(i)");
#endif
    return vm1_[i];
  }

  template <class T>
  inline typename Matrix<T>::const_reference Matrix<T>::operator()(Subscript i) const
  {
#ifdef SCPPNT_BOUNDS_CHECK
    if (i < 1 || i > mn_) throw BoundsError("SCPPNT::Matrix::operator(i) const");
#endif
    return vm1_[i];
  }

  template <class T>
  inline typename Matrix<T>::reference Matrix<T>::operator()(Subscript i, Subscript j)
  {
#ifdef SCPPNT_BOUNDS_CHECK
    if (i < 1 || i > m_ || j < 1 || j > n_) throw BoundsError("SCPPNT::Matrix::operator(i,j)");
#endif
    return rowm1_[i][j];
  }

  template <class T>
  inline typename Matrix<T>::const_reference Matrix<T>::operator() (Subscript i, Subscript j) const
  {
#ifdef SCPPNT_BOUNDS_CHECK
    if (i < 1 || i > m_ || j < 1 || j > n_) throw BoundsError("SCPPNT::Matrix::operator(i,j) const");
#endif
    return rowm1_[i][j];
  }

  /* *********************** interators ***************************/

  /*  Iterator over all elements of matrix */
  template <class T>
  inline typename Matrix<T>::iterator Matrix<T>::begin()
  {
#ifdef SCPPNT_BOUNDS_CHECK
    return iterator(v_, 1, mn_);
#else
    return v_;
#endif
  }

  template <class T>
  inline typename Matrix<T>::const_iterator Matrix<T>::begin() const
  {
#ifdef SCPPNT_BOUNDS_CHECK
    return const_iterator(v_, 1, mn_);
#else
    return v_;
#endif
  }

  template <class T>
  inline typename Matrix<T>::iterator Matrix<T>::end()
  {
#ifdef SCPPNT_BOUNDS_CHECK
    /* The pointer returned does not point to a valid element - it
     is one past the last element. An error will occur if an attempt
     is made to dereference the pointer */
    return iterator(v_, 1, mn_) += mn_;
#else
    return v_ + mn_;
#endif
  }

  template <class T>
  inline typename Matrix<T>::const_iterator Matrix<T>::end() const
  {
#ifdef SCPPNT_BOUNDS_CHECK
    /* The pointer returned does not point to a valid element - it
     is one past the last element. An error will occur if an attempt
     is made to dereference the pointer */
    return const_iterator(v_, 1, mn_) += mn_; /* Typo corrected. (ww, 12/2/2007) */
#else
    return v_ + mn_;
#endif
  }

  /* row and column iterators */
  template <class T>
  inline typename Matrix<T>::rows_iterator Matrix<T>::begin_rows()
  {
#ifdef SCPPNT_BOUNDS_CHECK
    return rows_iterator(row_, 1, m_);
#else
    return row_;
#endif
  }

  template <class T>
  inline typename Matrix<T>::const_rows_iterator Matrix<T>::begin_rows() const
  {
#ifdef SCPPNT_BOUNDS_CHECK
    return const_rows_iterator(const_row_, 1, m_);
#else
    return (const_rows_iterator) row_;
#endif
  }

  template <class T>
  inline typename Matrix<T>::columns_iterator Matrix<T>::begin_columns()
  {
#ifdef SCPPNT_BOUNDS_CHECK
    return columns_iterator(column_, 1, n_);
#else
    return column_;
#endif
  }

  template <class T>
  inline typename Matrix<T>::const_columns_iterator Matrix<T>::begin_columns() const
  {
#ifdef SCPPNT_BOUNDS_CHECK
    return const_columns_iterator(const_column_, 1, n_);
#else
    return const_column_;
#endif
  }

  template <class T>
  inline typename Matrix<T>::rows_iterator Matrix<T>::end_rows()
  {
#ifdef SCPPNT_BOUNDS_CHECK
    return rows_iterator(row_ + m_, 1, 1);
#else
    return row_ + m_;
#endif
  }

  template <class T>
  inline typename Matrix<T>::const_rows_iterator Matrix<T>::end_rows() const
  {
#ifdef SCPPNT_BOUNDS_CHECK
    return const_rows_iterator(const_row_ + m_, 1, 1);
#else
    return (const_rows_iterator) row_ + m_;
#endif
  }

  template <class T>
  inline typename Matrix<T>::columns_iterator Matrix<T>::end_columns()
  {
#ifdef SCPPNT_BOUNDS_CHECK
    return columns_iterator(column_ + n_, 1, 1);
#else
    return column_ + n_;
#endif
  }

  template <class T>
  inline typename Matrix<T>::const_columns_iterator Matrix<T>::end_columns() const
  {
#ifdef SCPPNT_BOUNDS_CHECK
    return const_columns_iterator(const_column_ + n_, 1, 1);
#else
    return const_column_ + n_;
#endif
  }

  template <class T>
  inline typename Matrix<T>::row_iterator Matrix<T>::begin_row(Subscript index)
  {
#ifdef SCPPNT_BOUNDS_CHECK
    if (index < 1 || index > m_) throw BoundsError("SCPPNT::Matrix::begin_row()");
#endif
    return row_[index-1];
  }

  template <class T>
  inline typename Matrix<T>::const_row_iterator Matrix<T>::begin_row(Subscript index) const
  {
#ifdef SCPPNT_BOUNDS_CHECK
    if (index < 1 || index > m_) throw BoundsError("SCPPNT::Matrix::begin_row() const");
    return const_row_[index-1];
#else
    return (const_row_iterator) row_[index-1];
#endif
  }

  template <class T>
  inline typename Matrix<T>::column_iterator Matrix<T>::begin_column(Subscript index)
  {
#ifdef SCPPNT_BOUNDS_CHECK
    if (index < 1 || index > n_) throw BoundsError("SCPPNT::Matrix::begin_column()");
#endif
    return column_[index-1];
  }

  template <class T>
  inline typename Matrix<T>::const_column_iterator Matrix<T>::begin_column(Subscript index) const
  {
#ifdef SCPPNT_BOUNDS_CHECK
    if (index < 1 || index > n_) throw BoundsError("SCPPNT::Matrix::begin_column() const");
#endif
    return const_column_[index-1];
  }

  template <class T>
  inline typename Matrix<T>::row_iterator Matrix<T>::end_row(Subscript index)
  {
#ifdef SCPPNT_BOUNDS_CHECK
    if (index < 1 || index > m_) throw BoundsError("SCPPNT::Matrix::end_row()");
#endif
    return row_[index];
  }

  template <class T>
  inline typename Matrix<T>::const_row_iterator Matrix<T>::end_row(Subscript index) const
  {
#ifdef SCPPNT_BOUNDS_CHECK
    if (index < 1 || index > m_) throw BoundsError("SCPPNT::Matrix::end_row() const");
    return const_row_[index];
#else
    return (const_row_iterator) row_[index];
#endif
  }

  template <class T>
  inline typename Matrix<T>::column_iterator Matrix<T>::end_column(Subscript index)
  {
#ifdef SCPPNT_BOUNDS_CHECK
    if (index < 1 || index > n_) throw BoundsError("SCPPNT::Matrix::end_column()");
#endif
    return column_[index];
  }

  template <class T>
  inline typename Matrix<T>::const_column_iterator Matrix<T>::end_column(Subscript index) const
  {
#ifdef SCPPNT_BOUNDS_CHECK
    if (index < 1 || index > n_) throw BoundsError("SCPPNT::Matrix::end_column() const");
#endif
    return const_column_[index];
  }

  // Returns number of elements in diagonal with first element (row, column)
  template <class T>
  inline Subscript Matrix<T>::diagonal_size(Subscript row, Subscript column) const
  {
    Subscript ncolumn = n_ - column + 1;
    Subscript nrow = m_ - row + 1;

    return (ncolumn > nrow) ? nrow : ncolumn;
  }

  /*!
   
   Returns iterator pointing to the first element of the matrix diagonal given by
   (row, column), (row+1, column+1), (row+2, column+2), ...
   For example, if row=1 and column=1 then an iterator pointing to the first element of the
   main diagnonal is returned.
   
   \param row 1-based row number for initial element of diagonal.
   \param column 1-based column number for initial element of diagonal.
   */
  template <class T>
  inline typename Matrix<T>::diag_iterator Matrix<T>::begin_diagonal(Subscript row, Subscript column)
  {
    Subscript num_elem = 0;

#ifdef SCPPNT_BOUNDS_CHECK

    if (row < 1 || row > m_ || column < 1 || column > n_) BoundsError("SCPPNT::Matrix::begin_diagonal(row, column)");
    num_elem = diagonal_size(row, column);

#endif

    /* Pass zero for number of elements if SCPPNT_BOUNDS_CHECK
     is not defined, since in that case the number of elements is not used. */
    return diag_iterator(v_ + column-1 + n_*(row-1), n_ + 1, num_elem);
  }

  /*!
   
   Returns iterator pointing to the first element of the matrix diagonal given by
   (row, column), (row+1, column+1), (row+2, column+2), ...
   For example, if row=1 and column=1 then an iterator pointing to the first element of the
   main diagnonal is returned.
   
   \param row 1-based row number for initial element of diagonal.
   \param column 1-based column number for initial element of diagonal.
   */
  template <class T>
  inline typename Matrix<T>::const_diag_iterator Matrix<T>::begin_diagonal(Subscript row, Subscript column) const
  {
    Subscript num_elem = 0;

#ifdef SCPPNT_BOUNDS_CHECK

    if (row < 1 || row > m_ || column < 1 || column > n_) BoundsError("SCPPNT::Matrix::begin_diagonal(row, column)");
    num_elem = diagonal_size(row, column);

#endif

    /* Pass zero for number of elements if SCPPNT_BOUNDS_CHECK
     is not defined, since in that case the number of elements is not used. */
    return const_diag_iterator(v_ + column-1 + n_*(row-1), n_ + 1, num_elem);
  }

  /*!
   
   Returns iterator pointing to the one past last element of the matrix diagonal given by
   (row, column), (row+1, column+1), (row+2, column+2), ...
   For example, for a square matrix with n rows and columns end_diagonal(1, 1)
   returns an iterator pointing to the element (row+n, column+n), which is one past
   the last element of the main diagonal.
   
   \param row 1-based row number for initial element of diagonal.
   \param column 1-based column number for initial element of diagonal.
   */
  template <class T>
  inline typename Matrix<T>::diag_iterator Matrix<T>::end_diagonal(Subscript row, Subscript column)
  {
    diag_iterator it = begin_diagonal(row, column);
    return it + diagonal_size(row, column);
  }

  /*!
   Returns iterator pointing to the one past last element of the matrix diagonal given by
   (row, column), (row+1, column+1), (row+2, column+2), ...
   For example, for a square matrix with n rows and columns end_diagonal(1, 1)
   returns an iterator pointing to the element (row+n, column+n), which is one past
   the last element of the main diagonal.
   
   \param row 1-based row number for initial element of diagonal.
   \param column 1-based column number for initial element of diagonal.
   */
  template <class T>
  inline typename Matrix<T>::const_diag_iterator Matrix<T>::end_diagonal(Subscript row, Subscript column) const
  {
    const_diag_iterator it = begin_diagonal(row, column);
    return it + diagonal_size(row, column);
  }

  /* ***************************  assignment operators  ********************************/

  /****** The following three member templates are defined in the class definition ******/

  // matrix operator +=
  // template <class T>
  // template <class MAT>
  // Matrix<T>& Matrix<T>::operator+=(const MAT &rhs) {matadd(*this, rhs); return *this;}

  // matrix operator -=
  // template <class T>
  // template <class MAT>
  // Matrix<T>& Matrix<T>::operator-=(const MAT &rhs) {matsub(*this, rhs); return *this;}

  // matrix operator *=
  // template <class T>
  // template <class MAT>
  // Matrix<T>& Matrix<T>::operator*=(const MAT &rhs) {matmult_assign(*this, rhs); return *this;}


  // scalar operator += (add a scalar to each element)
  template <class T>
  Matrix<T>& Matrix<T>::operator+=(const T& value)
  {
    double *thisiter = v_;
    for (int i = mn_; i--; ++thisiter)
    {
      *thisiter += value;
    }
    return *this;
  }

  // scalar operator -= (subtract a scalar from each element)
  template <class T>
  Matrix<T>& Matrix<T>::operator-=(const T& value)
  {
    double *thisiter = v_;
    for (int i = mn_; i--; ++thisiter)
    {
      *thisiter -= value;
    }
    return *this;
  }

  // scalar operator *= (multiply each element by scalar)
  template <class T>
  Matrix<T>& Matrix<T>::operator*=(const T& value)
  {
    double *thisiter = v_;
    for (Subscript i = mn_; i--; ++thisiter)
    {
      *thisiter *= value;
    }
    return *this;
  }

  // scalar operator /= (divide each element by scalar)
  template <class T>
  Matrix<T>& Matrix<T>::operator/=(const T& value)
  {
    double *thisiter = v_;
    for (Subscript i = mn_; i--; ++thisiter)
    {
      *thisiter /= value;
    }
    return *this;
  }

#ifndef SCPPNT_NO_IO
  /* ***************************  I/O  ********************************/

  /*!
   \brief Write a matrix to an output stream.
   
   Writes the number of rows and columns, followed
   by a new line. The following lines contain the
   elements in the rows of the matrix.
   */
  template <class T>
  std::ostream& operator<<(std::ostream &s, const Matrix<T> &A)
  {
    Subscript M = A.num_rows();
    Subscript N = A.num_columns();

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

  /*!
   \brief Read elements of a matrix to an input stream.
   
   Reads the number of rows and columns, followed
   by the elements in each row of the matrix.
   */
  template <class T>
  std::istream& operator>>(std::istream &s, Matrix<T> &A)
  {
    Subscript M, N;

    s >> M >> N;

    A.newsize(M,N);

    for (Subscript i=0; i<M; i++)
    for (Subscript j=0; j<N; j++)
    {
      s >> A[i][j];
    }
    return s;
  }
#endif // SCPPNT_NO_IO
  // *******************[ basic matrix operators ]***************************

  //! Add two matrices and return new matrix giving sum
  template <class T, class M>
  inline Matrix<T> operator+(const Matrix<T> &A, const M &B)
  {
    return Matrix<T>(A) += B;
  }

  //! Subtract matrix B from matrix A and return matrix giving difference
  template <class T, class M>
  inline Matrix<T> operator-(const Matrix<T> &A, const M &B)
  {
    return Matrix<T>(A) -= B;
  }

  //! Multiplication operator returning a new matrix containing the matrix product A*B
  template <class T, class M>
  inline Matrix<T> operator*(const Matrix<T> &A, const M &B)
  {
    return matmult<Matrix<T>, Matrix<T>, M>(A,B);
  }

  //! Return a new vector containing the product of the matrix A and the vector x
  template <class T>
  inline Vector<T> operator*(const Matrix<T> &A, const Vector<T> &x)
  {
    return matrix_times_vector< Matrix<T>, Vector<T> >(A, x);
  }

} // namespace SCPPNT

#endif
// SCPPNT_CMAT_H
