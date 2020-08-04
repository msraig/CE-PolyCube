/*! \file trisolve.h
 \brief Functions for triangular solves.

 Functions
 Modified version of trisolve.h from the Template Numerical Toolkit (TNT) 
 using matrix iterators rather than matrix subscripting.

 Contains modified version of the functions in trisolve.h from the
 Template Numerical Toolkit (TNT) using matrix iterators rather than
 matrix subscripting.
 */

/*

 Simple C++ Numerical Toolkit (SCPPNT)
 http://www.smallwaters.com/software/cpp/scppnt.html
 This release updates original work contributed by 
 Brad Hanson (http://www.b-a-h.com/) in 2001.

 */

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

// Triangular Solves

#ifndef SCPPNT_TRISLV_H
#define SCPPNT_TRISLV_H

#ifdef SCPPNT_NO_DIR_PREFIX
#include "scppnt.h"
#include "triang.h"
#else
#include "scppnt/scppnt.h"
#include "scppnt/triang.h"
#endif

namespace SCPPNT
{

  //! Solve linear system Ax=b using lower triangular portion of A including diagonal
  template<class MaTriX, class VecToR> VecToR Lower_triangular_solve(const MaTriX &A,
      const VecToR &b)
  {
    Subscript N = A.num_rows();

    // make sure matrix sizes agree; A must be square

    //assert(A.num_columns() == N);
    //assert(b.dim() == N);
    if (A.num_columns() != N || b.dim() != N)
      throw BoundsError("SCPPNT::Lower_triangular_solve");

    VecToR x(N);
    x = 0;

    Subscript i;
    typename VecToR::iterator ix = x.begin();
    typename VecToR::const_iterator ib = b.begin();
    typename MaTriX::const_rows_iterator iA = A.begin_rows();
    typename MaTriX::const_diag_iterator dA = A.begin_diagonal(1, 1);
    for (i=1; i<=N; i++, ++ix, ++ib, ++iA, ++dA)
    {
      typename MaTriX::element_type tmp=0;

      typename MaTriX::const_row_iterator jA = *iA;
      typename VecToR::iterator jx = x.begin();
      for (Subscript j=1; j<i; j++, ++jA, ++jx)
        tmp += *jA * *jx; // tmp = tmp + A(i,j)*x(j);

      *ix = *ib - tmp; // x(i) =  (b(i) - tmp)/ A(i,i);
      *ix /= *dA;
    }

    return x;
  }

  //! Solve linear system Ax=b using lower triangular portion of A assuming unit diagonal
  template<class MaTriX, class VecToR> VecToR Unit_lower_triangular_solve(const MaTriX &A,
      const VecToR &b)
  {
    Subscript N = A.num_rows();

    // make sure matrix sizes agree; A must be square

    //assert(A.num_columns() == N);
    //assert(b.dim() == N);
    if (A.num_columns() != N || b.dim() != N)
      throw BoundsError("SCPPNT::Unit_lower_triangular_solve");

    VecToR x(N);
    x = 0;

    Subscript i;
    typename VecToR::iterator ix = x.begin();
    typename VecToR::const_iterator ib = b.begin();
    typename MaTriX::const_rows_iterator iA = A.begin_rows();
    for (i=1; i<=N; i++, ++ix, ++ib, ++iA)
    {

      typename MaTriX::element_type tmp=0;

      typename MaTriX::const_row_iterator jA = *iA;
      typename VecToR::iterator jx = x.begin();
      for (Subscript j=1; j<i; j++, ++jA, ++jx)
        tmp += *jA * *jx; // tmp = tmp + A(i,j)*x(j);

      *ix = *ib - tmp; // x(i) =  b(i) - tmp;
    }

    return x;
  }

  //! Solve linear system using lower triangular view A
  template<class MaTriX, class VecToR> VecToR linear_solve(const LowerTriangularView<MaTriX> &A,
      const VecToR &b)
  {
    return Lower_triangular_solve(A, b);
  }

  //! Solve linear system using unit lower triangular view A
  template<class MaTriX, class VecToR> VecToR linear_solve(
      const UnitLowerTriangularView<MaTriX> &A, const VecToR &b)
  {
    return Unit_lower_triangular_solve(A, b);
  }

  //********************** Upper triangular section ****************

  //! Solve linear system Ax=b using upper triangular portion of A including diagonal
  template<class MaTriX, class VecToR> VecToR Upper_triangular_solve(const MaTriX &A,
      const VecToR &b)
  {
    Subscript N = A.num_rows();

    // make sure matrix sizes agree; A must be square

    //assert(A.num_columns() == N);
    //assert(b.dim() == N);
    if (A.num_columns() != N || b.dim() != N)
      throw BoundsError("SCPPNT::Upper_triangular_solve");

    VecToR x(N);
    x = 0;

    Subscript i;
    typename VecToR::iterator ix = x.begin() + N - 1;
    typename VecToR::const_iterator ib = b.begin() + N - 1;
    typename MaTriX::const_rows_iterator iA = A.begin_rows() + N - 1;
    typename MaTriX::const_diag_iterator dA = A.begin_diagonal(1, 1) + N - 1;
    for (i=N; i>=1; i--, --ix, --ib, --iA, --dA)
    {

      typename MaTriX::element_type tmp=0;

      typename MaTriX::const_row_iterator jA = *iA + i;
      typename VecToR::iterator jx = x.begin() + i;
      for (Subscript j=i+1; j<=N; j++, ++jA, ++jx)
        tmp += *jA * *jx; // tmp = tmp + A(i,j)*x(j);

      *ix = *ib - tmp; // x(i) =  (b(i) - tmp)/ A(i,i);
      *ix /= *dA;
    }

    return x;
  }

  //! Solve linear system Ax=b using upper triangular portion of A assuming unit diagonal
  template<class MaTriX, class VecToR> VecToR Unit_upper_triangular_solve(const MaTriX &A,
      const VecToR &b)
  {
    Subscript N = A.num_rows();

    // make sure matrix sizes agree; A must be square

    //assert(A.num_columns() == N);
    //assert(b.dim() == N);
    if (A.num_columns() != N || b.dim() != N)
      throw BoundsError("SCPPNT::Unit_upper_triangular_solve");

    VecToR x(N);
    x = 0;

    Subscript i;
    typename VecToR::iterator ix = x.begin() + N - 1;
    typename VecToR::const_iterator ib = b.begin() + N - 1;
    typename MaTriX::const_rows_iterator iA = A.begin_rows() + N - 1;
    for (i=N; i>=1; i--, --ix, --ib, --iA)
    {

      typename MaTriX::element_type tmp=0;

      typename MaTriX::const_row_iterator jA = *iA + i;
      typename VecToR::iterator jx = x.begin() + i;
      for (Subscript j=i+1; j<i; j++, ++jA, ++jx)
        tmp += *jA * *jx; // tmp = tmp + A(i,j)*x(j);

      *ix = *ib - tmp; // x(i) =  b(i) - tmp;
    }

    return x;
  }

  //! Solve linear system using upper triangular view A
  template<class MaTriX, class VecToR> VecToR linear_solve(const UpperTriangularView<MaTriX> &A,
      const VecToR &b)
  {
    return Upper_triangular_solve(A, b);
  }

  //! Solve linear system using unit upper triangular view A
  template<class MaTriX, class VecToR> VecToR linear_solve(
      const UnitUpperTriangularView<MaTriX> &A, const VecToR &b)
  {
    return Unit_upper_triangular_solve(A, b);
  }

} // namespace SCPPNT

#endif
// SCPPNT_TRISLV_H
