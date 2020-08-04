/*! \file cholesky.h
 \brief Definition of function Cholesky_upper_factorization

 Contains modified version of the functions in cholesky.h from the
 Template Numerical Toolkit (TNT) using matrix iterators rather than
 matrix subscripting.

 */

/*  

 Simple C++ Numerical Toolkit (SCPPNT)
 http://www.smallwaters.com/software/cpp/scppnt.html
 This release updates original work contributed by 
 Brad Hanson (http://www.b-a-h.com/) in 2001.

 */

// Modified version of cholesky.h from:
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

#ifndef SCPPNT_CHOLESKY_H
#define SCPPNT_CHOLESKY_H

#include <cmath> // for sqrt
#ifdef BOOST_NO_STDC_NAMESPACE
namespace std
{ using ::sqrt;}
#endif

#ifdef SCPPNT_NO_DIR_PREFIX
#include "scppnt.h"
#include "scppnt_error.h"
#else
#include "scppnt/scppnt.h"
#include "scppnt/scppnt_error.h"
#endif

// index method

namespace SCPPNT
{

  /*! \brief Compute Cholesky factorization
   
   Modified version of Cholesky_upper_factorization from the Template
   Numerical Toolkit (TNT) using matrix iterators rather than matrix
   subscripting.
   
   Only upper part of A is used.  Cholesky factor is returned in lower
   part of L.  Returns 0 if successful, 1 otherwise.
   
   */
  template<class SPDMatrix, class SymmMatrix> int Cholesky_upper_factorization(SPDMatrix &A,
      SymmMatrix &L)
  {
    Subscript M = A.dim(1);
    Subscript N = A.dim(2);

    // assert(M == N);                 // make sure A is square
    if (M != N)
      throw BoundsError("Cholesky_upper_factorization");

    // readjust size of L, if necessary

    if (M != L.dim(1) || N != L.dim(2))
      L = SymmMatrix(N, N);

    Subscript i, j, k;

    typename SPDMatrix::element_type dot;

    typename SymmMatrix::diag_iterator diagL = L.begin_diagonal(1, 1);
    typename SPDMatrix::diag_iterator diagA = A.begin_diagonal(1, 1);

    for (j=1; j<=N; j++, ++diagL, ++diagA) // form column j of L
    {
      dot= 0;

      typename SymmMatrix::row_iterator irL = L.begin_row(j);

      for (i=1; i<j; i++) // for k= 1 TO j-1
      {
        dot += *irL * *irL; //dot = dot +  L(j,i)*L(j,i);
        ++irL;
      }

      *diagL = *diagA - dot; // L(j,j) = A(j,j) - dot;

      typename SymmMatrix::column_iterator icL = L.begin_column(j) + j;
      typename SPDMatrix::row_iterator irA = A.begin_row(j) + j;

      irL = L.begin_row(j);

      for (i=j+1; i<=N; i++, ++icL, ++irA)
      {
        dot = 0;
        typename SymmMatrix::row_iterator iL = L.begin_row(i);
        typename SymmMatrix::row_iterator jL =irL;
        for (k=1; k<j; k++, ++iL, ++jL)
          dot += *iL * *jL; // dot = dot +  L(i,k)*L(j,k);
        *icL = *irA - dot; // L(i,j) = A(j,i) - dot;
      }

      if (*diagL <= 0.0)
        return 1; // if (L(j,j) <= 0.0) return 1;

      *diagL = std::sqrt(*diagL); // L(j,j) = sqrt( L(j,j) );

      icL = L.begin_column(j) + j;

      for (i=j+1; i<=N; i++, ++icL)
        *icL /= *diagL; // L(i,j) = L(i,j) / L(j,j);

    }
    return 0;
  }
  
} // namespace SCPPNT

#endif
// SCPPNT_CHOLESKY_H
