/*! \file transv.h
 \brief Definition of class Transpose_View
 
 The Transpose_View class allows working with a tranpose of a matrix without
 creating a new matrix containing the transpose.

 */

/*
 Simple C++ Numerical Toolkit (SCPPNT)
 http://www.smallwaters.com/software/cpp/scppnt.html
 This release updates original work contributed by 
 Brad Hanson (http://www.b-a-h.com/) in 2001.

 */

// Modified from transv.h in:
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

// Matrix Transpose Views

#ifndef SCPPNT_TRANSV_H
#define SCPPNT_TRANSV_H

#ifdef SCPPNT_NO_DIR_PREFIX
#include "scppnt.h"
#include "scppnt_error.h"
#include "matop.h"
#include "vec.h"
#else
#include "scppnt/scppnt.h"
#include "scppnt/scppnt_error.h"
#include "scppnt/matop.h"
#include "scppnt/vec.h"
#endif

#ifndef SCPPNT_NO_IO
#include <iostream>
#endif

namespace SCPPNT
{

  /*!
   \brief Allows use of the transpose of a matrix without creating a new matrix.
   
   This class allows working with a tranpose of a matrix without
   creating a new matrix containing the transpose. The transpose
   view is read-only. Only operations which do not change the
   content of the matrix can be performed with a transpose view.

   */
  template<class Array2D> class Transpose_View
  {

private:

    const Array2D & A_; //!< Matrix that is the basis of the tranpose view

public:

    typedef typename Array2D::size_type size_type; //!< Subscript type
    typedef typename Array2D::value_type value_type; //!< Type of elements stored in matrix
    typedef typename Array2D::element_type element_type; //!< Type of elements stored in matrix
    typedef typename Array2D::pointer pointer; //!< Pointer to type stored in matrix
    typedef typename Array2D::iterator iterator; //!< Iterator over elements in matrix
    typedef typename Array2D::reference reference; //!< Reference to type stored in matrix
    typedef typename Array2D::const_iterator const_iterator; //!< Iterator over elements in constant matrix
    typedef typename Array2D::const_reference const_reference; //!< Reference to type stored in constant matrix

    /**** Iterators ****/

    //! Iterator over constant elements of a row
    typedef typename Array2D::const_column_iterator const_row_iterator;
    //! Iterator over constant elements of a column
    typedef typename Array2D::const_row_iterator const_column_iterator;

    //! Iterator over row iterators (points to const_row_iterator object for a row)
    typedef typename Array2D::const_columns_iterator const_rows_iterator;
    //! Iterator over column iterators (points to const_column_iterator object for a column)
    typedef typename Array2D::const_rows_iterator const_columns_iterator;

    //! Iterator over constant elements of a matrix diagonal
    typedef typename Array2D::const_diag_iterator const_diag_iterator;

    //! Return matrix that is the basis of the transpose view
    const Array2D & array() const
    {
      return A_;
    }

    //! Return number of rows in matrix
    Subscript num_rows() const
    {
      return A_.num_columns();
    }
    //! Return number of columns in matrix
    Subscript num_columns() const
    {
      return A_.num_rows();
    }
    //! Returns lower bound of subscript
    Subscript lbound() const
    {
      return A_.lbound();
    }
    //! Return number of rows if argument is 1, number of columns if argument is 2
    Subscript dim(Subscript i) const
    {
      return (i == 1) ? A_.dim(2) : A_.dim(i-1);
    }

    //! Return total number of elements in matrix
    Subscript size() const
    {
      return A_.size();
    }

    /**** constructors ****/

    //! Copy constructor
    Transpose_View(const Transpose_View<Array2D> &A) :
      A_(A.A_)
    {
    }
    ;

    //! Constructor from matrix
    Transpose_View(Array2D &A) :
      A_(A)
    {
    }
    ;

    /**** Subscripting ****/

    //! Return iterator to constant elements in row i+1 (0-offset)
    const_row_iterator operator[](Subscript i) const
    {
      return A_.column(i+1);
    }

    //! Return constant element in row i and column j, where i and j are 1-offset indices
    const_reference operator()(Subscript i, Subscript j) const
    {
      return A_(j, i);
    }

    /**** Iterators ****/

    //! Return iterator pointing to row iterator for first row (iterator over row iterators)
    const_rows_iterator begin_rows() const
    {
      return (const_rows_iterator) A_.begin_columns();
    }

    //! Return iterator pointing to row iterator for one past last row
    const_rows_iterator end_rows() const
    {
      return (const_rows_iterator) A_.end_columns();
    }

    //! Return iterator pointing to column iterator for first column (iterator over column iterators)
    const_columns_iterator begin_columns() const
    {
      return (const_columns_iterator) A_.begin_rows();
    }

    //! Return iterator pointing to column iterator for one past last column
    const_columns_iterator end_columns() const
    {
      return (const_columns_iterator) A_.end_rows();
    }

    // Return row or column iterator for particular row or column.
    // Index is one-offset (first row or column is 1).
    // row(N+1) returns a row_iterator to one past the last row, where there are N rows
    // column(M+1) returns a column_iterator to one past the last column, where there are M columns
    //

    //! Return iterator pointing to first element in row 'index' (1-offset)
    const_row_iterator begin_row(Subscript index) const
    {
      return A_.begin_column(index);
    }

    //! Return iterator pointing to one past last element in row 'index' (1-offset)
    const_row_iterator end_row(Subscript index) const
    {
      return A_.end_column(index);
    }

    //! Return iterator pointing to first element in column 'index' (1-offset)
    const_column_iterator begin_column(Subscript index) const
    {
      return A_.begin_row(index);
    }

    //! Return iterator pointing to one past last element in column 'index' (1-offset)
    const_column_iterator end_column(Subscript index) const
    {
      return A_.end_row(index);
    }

    // Iterator over all elements of matrix.
    // Should only be used when the order in which the iteration takes place
    // over all elements of the matrix does not matter

    //! Return iterator pointing to first element of matrix
    const_iterator begin() const
    {
      return A_.begin();
    }

    //! Return iterator pointing to one past last element of matrix
    const_iterator end() const
    {
      return A_.end();
    }

    // Iterator over diagonals.

    /*!
     
     \brief Returns iterator pointing to first element of a matrix diagonal
     
     \param row is the row corresponding to the first element of the diagonal
     \param column is the column corresponding to the first element of the diagonal
     */
    const_diag_iterator begin_diagonal(Subscript row, Subscript column) const
    {
      return A_.begin_diagonal(column, row);
    }

    /*! 
     Returns iterator pointing to one past last element of matrix diagonal
     \param row is the row corresponding to the first element of the diagonal
     \param column is the column corresponding to the first element of the diagonal
     */
    const_diag_iterator end_diagonal(Subscript row, Subscript column) const
    {
      return A_.end_diagonal(column, row);
    }
  };

#ifndef SCPPNT_NO_IO
  /*!
   \brief Write a matrix to an output stream.
   
   Writes the number of rows and columns, followed
   by a new line. The following lines contain the
   elements in the rows of the matrix.
   */
  template<class Array2D> std::ostream& operator<<(std::ostream &s,
      const Transpose_View<Array2D> &A)
  {
    Subscript M=A.num_rows();
    Subscript N=A.num_columns();

    Subscript start = A.lbound();
    Subscript Mend = M + A.lbound() - 1;
    Subscript Nend = N + A.lbound() - 1;

    s << M << "  " << N << std::endl;
    for (Subscript i=start; i<=Mend; i++)
    {
      for (Subscript j=start; j<=Nend; j++)
      {
        s << A(i, j) << " ";
      }
      s << std::endl;
    }

    return s;
  }
#endif

  // *******************[ basic matrix operators ]***************************

  //! Add two matrices and return new matrix giving sum
  template<class Array2D, class M> inline M operator+(const Transpose_View<Array2D> &A, const M &B)
  {
    M sum(A.num_rows(), A.num_columns());
    matadd(A, B, sum);
    return sum;
  }

  //! Subtract matrix B from matrix A and return matrix giving difference
  template<class Array2D, class M> inline M operator-(const Transpose_View<Array2D> &A, const M &B)
  {
    M diff(A.num_rows(), A.num_columns());
    matsub(A, B, diff);
    return diff;
  }

  //! Multiplication operator returning a new matrix containing the matrix product A*B
  template<class Array2D, class M> inline M operator*(const Transpose_View<Array2D> &A, const M &B)
  {
    M c(A.num_rows(), B.num_columns());
    matmult(c, A, B);
    return c;
  }

  //! Return a new vector containing the product of the matrix A and the vector x
  template<class Array2D> inline Vector<typename Array2D::value_type> operator*(
      const Transpose_View<Array2D> &A, const Vector<typename Array2D::value_type> &x)
  {
    return matrix_times_vector<Transpose_View<Array2D>,
    Vector <typename Array2D::value_type> >(A, x);
  }

} // namespace SCPPNT


#endif
// SCPPNT_TRANSV_H
