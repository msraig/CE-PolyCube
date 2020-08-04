/*! \file matop.h
 \brief Template functions which implement matrix operations.
 
 Definition of template functions transpose, matadd, matsub,
 mult_element, matmult, matmult_assign, matrix_times_vector,
 and MatTimesVec. 
 These functions are used in matrix classes to implement matrix
 operations such as addition, subtraction, and multiplication.

 */

/*

 Simple C++ Numerical Toolkit (SCPPNT)
 http://www.smallwaters.com/software/cpp/scppnt.html
 This release updates original work contributed by 
 Brad Hanson (http://www.b-a-h.com/) in 2001.

 */

#ifndef SCPPNT_MATOP_H
#define SCPPNT_MATOP_H

#ifdef SCPPNT_NO_DIR_PREFIX
#include "scppnt.h"
#include "slice_pointer.h"
#include "vec.h"
#include "scppnt_error.h"
#else
#include "scppnt/scppnt.h"
#include "scppnt/slice_pointer.h"
#include "scppnt/vec.h"
#include "scppnt/scppnt_error.h"
#endif

namespace SCPPNT
{

  //! Return a new matrix containing the transpose of A
  template<class MAT> MAT transpose(const MAT &A)
  {
    Subscript M = A.num_rows();
    Subscript N = A.num_columns();

    MAT S(N, M);
    Subscript i, j;

    typename MAT::const_rows_iterator rowsA = A.begin_rows();
    typename MAT::columns_iterator columnsS = S.begin_columns();
    for (i=M; i--; ++rowsA, ++columnsS)
    {
      typename MAT::const_row_iterator rowA = *rowsA;
      typename MAT::column_iterator columnS = *columnsS;
      for (j=N; j--; ++rowA, ++columnS)
        *columnS = *rowA;
    }
    return S;
  }

  //! Copy all elements of matrix B to matrix A (B must be the same size as A)
  template<class MA, class MB> void matcopy(MA &A, const MB &B)
  {
    if (A.num_rows() != B.num_rows() || A.num_columns() != B.num_columns())
      throw BadDimension("SCPPNT::matcopy");

    typename MB::const_rows_iterator brows = B.begin_rows();
    typename MA::rows_iterator arows = A.begin_rows();
    for (int i = A.num_rows(); i--; ++brows, ++arows)
    {
      typename MB::const_row_iterator brow = *brows;
      typename MA::row_iterator arow = *arows;
      for (int j = A.num_columns(); j--; ++arow, ++brow)
      {
        *arow = *brow;
      }
    }

  }

  //! Add matrices A and B and place sum in C
  template<class MA, class MB, class MC> void matadd(const MA &A, const MB &B, MC &C)
  {
    Subscript nr = A.num_rows();
    Subscript nc = A.num_columns();
    if (nr != B.num_rows() || nc != B.num_columns() || nr != C.num_rows() || nc != C.num_columns())
      throw BadDimension("SCPPNT::matadd(A, B, C)");

    typename MB::const_rows_iterator brows = B.begin_rows();
    typename MA::const_rows_iterator arows = A.begin_rows();
    typename MC::rows_iterator crows = C.begin_rows();
    for (int i = nr; i--; ++brows, ++arows, ++crows)
    {
      typename MB::const_row_iterator brow = *brows;
      typename MA::const_row_iterator arow = *arows;
      typename MC::row_iterator crow = *crows;
      for (int j = nc; j--; ++arow, ++brow, ++crow)
      {
        *crow = *arow + *brow;
      }
    }
  }

  //! Add matrices A and B and replace A with result
  template<class MA, class MB> void matadd(MA &A, const MB &B)
  {

    if (A.num_rows() != B.num_rows() || A.num_columns() != B.num_columns())
      throw BadDimension("SCPPNT::matadd");

    typename MB::const_rows_iterator brows = B.begin_rows();
    typename MA::rows_iterator arows = A.begin_rows();
    for (int i = A.num_rows(); i--; ++brows, ++arows)
    {
      typename MB::const_row_iterator brow = *brows;
      typename MA::row_iterator arow = *arows;
      for (int j = A.num_columns(); j--; ++arow, ++brow)
      {
        *arow += *brow;
      }
    }

  }

  //! Substract B from A and place difference in C
  template<class MA, class MB, class MC> void matsub(const MA &A, const MB &B, MC &C)
  {
    Subscript nr = A.num_rows();
    Subscript nc = A.num_columns();
    if (nr != B.num_rows() || nc != B.num_columns() || nr != C.num_rows() || nc != C.num_columns())
      throw BadDimension("SCPPNT::matsub(A, B, C)");

    typename MB::const_rows_iterator brows = B.begin_rows();
    typename MA::const_rows_iterator arows = A.begin_rows();
    typename MC::rows_iterator crows = C.begin_rows();
    for (int i = nr; i--; ++brows, ++arows, ++crows)
    {
      typename MB::const_row_iterator brow = *brows;
      typename MA::const_row_iterator arow = *arows;
      typename MC::row_iterator crow = *crows;
      for (int j = nc; j--; ++arow, ++brow, ++crow)
      {
        *crow = *arow - *brow;
      }
    }
  }

  //! Subtract B from A and replace A with result
  template<class MA, class MB> void matsub(MA &A, const MB &B)
  {

    if (A.num_rows() != B.num_rows() || A.num_columns() != B.num_columns())
      throw BadDimension("SCPPNT::matsub");

    typename MB::const_rows_iterator brows = B.begin_rows();
    typename MA::rows_iterator arows = A.begin_rows();
    for (int i = A.num_rows(); i--; ++brows, ++arows)
    {
      typename MB::const_row_iterator brow = *brows;
      typename MA::row_iterator arow = *arows;
      for (int j = A.num_columns(); j--; ++arow, ++brow)
      {
        *arow -= *brow;
      }
    }
  }

  // ******************* Functions which perform matrix multiplication ***************************


  //! Return new matrix containing the products of the corresponding elements of A and B (element-by-element product)
  template<class MR, class MA, class MB> MR mult_element(const MA &A, const MB &B)
  {
    Subscript M = A.num_rows();
    Subscript N = A.num_columns();

    if (M != B.num_rows() || N != B.num_columns())
      throw BadDimension("SCPPNT::mult_element");

    MR tmp(M, N);
    Subscript i, j;

    typename MA::const_rows_iterator rowsA = A.begin_rows();
    typename MB::const_rows_iterator rowsB = B.begin_rows();
    typename MR::rows_iterator rowstmp = tmp.begin_rows();
    for (i=M; i--; ++rowsA, ++rowsB, ++rowstmp)
    {
      typename MA::const_row_iterator rowA = *rowsA;
      typename MB::const_row_iterator rowB = *rowsB;
      typename MR::row_iterator rowtmp = *rowstmp;
      for (j=N; j--; ++rowA, ++rowB, ++rowtmp)
        *rowtmp = *rowA * *rowB;
    }
    return tmp;
  }

  //! Return a new matrix containing the matrix product A*B
  template<class MR, class MA, class MB> MR matmult(const MA &A, const MB &B)
  {
    Subscript M = A.num_rows();
    Subscript N = A.num_columns();
    Subscript K = B.num_columns();

    if (N != B.num_rows())
      throw BadDimension("SCPPNT::matmult(A, B)");

    MR tmp(M, K); // New matrix to be returned
    typename MR::element_type sum;

    typename MR::rows_iterator rowsTmp = tmp.begin_rows();
    typename MA::const_rows_iterator rowsA = A.begin_rows();
    for (Subscript i=M; i--; ++rowsTmp, ++rowsA)
    {
      typename MR::row_iterator rowTmp = *rowsTmp;
      typename MB::const_columns_iterator columnsB = B.begin_columns();
      for (Subscript k=K; k--; ++rowTmp, ++columnsB)
      {
        typename MA::const_row_iterator rowA = *rowsA;
        typename MB::const_column_iterator columnB = *columnsB;
        sum = 0;
        for (Subscript j=N; j--; ++rowA, ++columnB)
          sum += *rowA * *columnB; //sum = sum +  A[i][j] * B[j][k];

        *rowTmp = sum; // tmp[i][k] = sum; 
      }
    }
    return tmp;
  }

  //! Calculate matrix product A*B and put result in C
  template<class MC, class MA, class MB> void matmult(MC &C, const MA &A, const MB &B)
  {
    Subscript M = A.num_rows();
    Subscript N = A.num_columns();
    Subscript K = B.num_columns();

    if (N != B.num_rows())
      throw BadDimension("SCPPNT::matmult(C, A, B)");

    C.newsize(M, K);

    typename MC::element_type sum;

    typename MC::rows_iterator rowsC = C.begin_rows();
    typename MA::const_rows_iterator rowsA = A.begin_rows();
    for (Subscript i=M; i--; ++rowsC, ++rowsA)
    {
      typename MC::row_iterator rowC = *rowsC;
      typename MB::const_columns_iterator columnsB = B.begin_columns();
      for (Subscript k=K; k--; ++rowC, ++columnsB)
      {
        typename MA::const_row_iterator rowA = *rowsA;
        typename MB::const_column_iterator columnB = *columnsB;
        sum = 0;
        for (Subscript j=N; j--; ++rowA, ++columnB)
          sum += *rowA * *columnB; //sum = sum +  A[i][j] * B[j][k];

        *rowC = sum; // C[i][k] = sum; 
      }
    }
  }

  //! Matrix multiplication of A*B where A is overwritten to store result (requires A*B and A to be same size)
  template<class MA, class MB> void matmult_assign(MA &A, const MB &B)
  {
    Subscript i, j, k;
    Subscript M = A.num_rows();
    Subscript N = A.num_columns();

    // Number of columns of A must equal the number of rows of B, and B must be a square matrix
    if (N != B.num_rows() || B.num_columns() != B.num_rows())
      throw BadDimension("SCPPNT::matmult_assign");

    // temp_vec is used to store a row of A
    // during multiplication involving that row to allow
    // the result to be stored in A
    typedef typename MA::element_type A_element_type;
    A_element_type *temp_vec = new A_element_type[N];

    typename MA::rows_iterator irows = A.begin_rows();
    for (i = M; i--; ++irows) // loop over rows of A
    {
      typename MB::const_columns_iterator icolumns = B.begin_columns();

      // Compute dot product of current row of A and first column of B
      // and put copy of original current row of A in temp_vec
      A_element_type *iv = temp_vec;
      typename MB::const_column_iterator icolumn = *icolumns;
      typename MA::row_iterator irow = *irows;
      A_element_type sum = 0.0;
      for (j = N; j--; ++irow, ++icolumn)
      {
        *iv = *irow;
        ++iv;
        sum += *irow * *icolumn;
      }

      // Replace first element of current row of A with 
      // dot product of current row of A and first column of B
      irow = *irows;
      *irow = sum;

      // dot products of current row of A and columns 2,3,...,N of B
      ++irow; // points to second element of current row of A
      ++icolumns; // points to iterator over elements in second column of B
      for (j=N-1; j--; ++icolumns, ++irow) // loop over columns 2,...,N of B
      {
        iv = temp_vec;
        icolumn = *icolumns;
        sum = 0.0;
        for (k=N; k--; ++iv, ++icolumn)
          sum += *iv * *icolumn;
        *irow = sum;
      }
    }
    delete [] temp_vec;
  }

  /*! \brief Multiply a matrix times a vector using iterators to access the vector elements
   
   
   \param A Matrix to multiply times vector.
   \param begin Iterator over elements of vector to be multiplied by matrix
   \param result Iterator over vector to hold result of multiplication.
   */
  template<class MAT, class IT1, class IT2> void MatTimesVec(const MAT &A, IT1 begin, IT2 result)
  {

    typedef typename MAT::value_type value_type;

    Subscript N = A.num_columns();
    Subscript M = A.num_rows();

    typename MAT::const_rows_iterator irows = A.begin_rows();
    for (Subscript i=M; i--; ++irows, ++result)
    {
      value_type sum = 0;
      typename MAT::const_row_iterator rowi = *irows;
      IT1 iv = begin;
      for (Subscript j=N; j--; ++iv, ++rowi)
      {
        sum += *rowi * *iv;
      }

      *result = sum;
    }
  }

  //! Multiply a matrix times a vector and return the result in a new vector
  template<class M, class V> inline V matrix_times_vector(const M &A, const V &x)
  {
    if (A.num_columns() != x.size())
      throw BadDimension("SCPPNT::matrix_times_vector()");

    V tmp(A.num_rows());
    MatTimesVec(A, x.begin(), tmp.begin());

    return tmp;
  }

} // namespace SCPPNT

#endif
// SCPPNT_MATOP_H
