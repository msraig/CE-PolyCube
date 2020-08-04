/*! \file matalg.h
 \brief Functions that apply function objects to rows or columns of matrices.
 
 */

/*

 Functions that apply function objects to rows or columns of matrices.
 
 Simple C++ Numerical Toolkit (SCPPNT)
 http://www.smallwaters.com/software/cpp/scppnt.html
 This release updates original work contributed by 
 Brad Hanson (http://www.b-a-h.com/) in 2001.

 */

#ifndef SCPPNT_MATALG_H
#define SCPPNT_MATALG_H
#ifdef SCPPNT_NO_DIR_PREFIX
#include "scppnt.h"
#else
#include "scppnt/scppnt.h"
#endif

namespace SCPPNT
{

  // Apply a function to each row of a matrix
  template<class M, class FUNC> void apply_rows(M &matrix, FUNC &f)
  {
    typename M::rows_iterator irows = matrix.begin_rows();

    Subscript ncolumns = matrix.num_columns();

    for (Subscript i = matrix.num_rows(); i--; ++irows)
    {
      f(ncolumns, *irows);
    }
  }

  // Apply a function to each column of a matrix
  template<class M, class FUNC> void apply_columns(M &matrix, FUNC &f)
  {
    typename M::columns_iterator icolumns = matrix.begin_columns();

    Subscript nrows = matrix.num_rows();

    for (Subscript i = matrix.num_columns(); i--; ++icolumns)
    {
      f(nrows, *icolumns);
    }
  }

  /*
   Apply a function to each row of a matrix that returns a single value
   (the function is applied over columns of each row). These values are
   put into a vector that is returned (the number of elements in the
   vector is the number of rows in the matrix).
   */
  template<class M, class V, class FUNC> V over_columns(M &matrix, FUNC &f)
  {
    typename M::rows_iterator irows = matrix.begin_rows();
    V vec(matrix.num_rows());

    Subscript ncolumns = matrix.num_columns();

    typename V::iterator iv = vec.begin();
    for (Subscript i = matrix.num_rows(); i--; ++irows, ++iv)
    {
      *iv = f(ncolumns, *irows);
    }

    return vec;
  }

  /*
   Apply a function to each column of a matrix that returns a single
   value (the function is applied over rows of each column). These
   values are put into a vector that is returned (the number of
   elements in the vector is the number of columns in the matrix).
   */
  template<class M, class V, class FUNC> V over_rows(M &matrix, FUNC &f)
  {
    typename M::columns_iterator icolumns = matrix.begin_columns();
    V vec(matrix.num_columns());

    Subscript nrows = matrix.num_rows();

    typename V::iterator iv = vec.begin();
    for (Subscript i = matrix.num_columns(); i--; ++icolumns, ++iv)
    {
      *iv = f(nrows, *icolumns);
    }

    return vec;
  }

} // namespace SCPPNT

#endif // SCPPNT_MATALG_H
