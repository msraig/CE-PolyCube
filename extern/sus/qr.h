/*! \file qr.h
 \brief Definition of functions QR_factor, R_solve, and QR_solve.
 
 Contains modified version of the functions in qr.h from the
 Template Numerical Toolkit (TNT) using matrix iterators rather than
 matrix subscripting.

 */

/*

 Simple C++ Numerical Toolkit (SCPPNT)
 http://www.smallwaters.com/software/cpp/scppnt.html
 This release updates original work contributed by 
 Brad Hanson (http://www.b-a-h.com/) in 2001.

 */

// Modified version of qr.h from:
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

#ifndef SCPPNT_QR_H
#define SCPPNT_QR_H

#include <cmath>      //for sqrt() & fabs()
#ifdef BOOST_NO_STDC_NAMESPACE
namespace std
{ using ::fabs; using ::sqrt;}
#endif

#ifdef SCPPNT_NO_DIR_PREFIX
#include "scppnt.h"
#include "scppnt_error.h"
#include "scppnt_math.h"  // for sign()
#else
#include "scppnt/scppnt.h"
#include "scppnt/scppnt_error.h"
#include "scppnt/scppnt_math.h"  // for sign()
#endif

namespace SCPPNT
{

  /*!
   \brief Classical QR factorization, based on Stewart[1973].
   
   
   This algorithm computes the factorization of a matrix A
   into a product of an orthognal matrix (Q) and an upper triangular 
   matrix (R), such that QR = A.
   
   Returns:  
   
   0 if successful, 1 if A is detected to be singular
   
   Parameters:
   
   \param A   (in/out):  On input, A is square, Matrix(1:N, 1:N), that represents
   the matrix to be factored.
   On output, Q and R is encoded in the same Matrix(1:N,1:N)
   in the following manner:
   R is contained in the upper triangular section of A,
   except that R's main diagonal is in D.  The lower
   triangular section of A represents Q, where each
   column j is the vector  Qj = I - uj*uj'/pi_j.
   
   \param C  (output):    vector of Pi[j]
   
   \param D  (output):    main diagonal of R, i.e. D(i) is R(i,i)
   
   */
  template<class MaTRiX, class Vector> int QR_factor(MaTRiX &A, Vector& C, Vector &D)
  {
    //assert(A.lbound() == 1);        // ensure these are all
    if (A.lbound() != 1)
      throw BoundsError("SCPPNT::QR_factor");

    //assert(C.lbound() == 1);        // 1-based arrays and vectors
    if (C.lbound() != 1)
      throw BoundsError("SCPPNT::QR_factor");

    //assert(D.lbound() == 1);
    if (D.lbound() != 1)
      throw BoundsError("SCPPNT::R_factor");

    Subscript M = A.num_rows();
    Subscript N = A.num_columns();

    // assert(M == N);                 // make sure A is square
    if (M != N)
      throw BoundsError("QR_factor");

    Subscript i, j, k;
    typename MaTRiX::element_type eta, sigma, sum;

    // adjust the shape of C and D, if needed...

    if (N != C.size())
      C.newsize(N);
    if (N != D.size())
      D.newsize(N);

    typename MaTRiX::diag_iterator diagA = A.begin_diagonal(1, 1);
    typename Vector::iterator iC = C.begin();
    typename Vector::iterator iD = D.begin();
    typename MaTRiX::columns_iterator columnsA = A.begin_columns();
    typename Vector::iterator dn = D.begin() + N - 1;
    typename MaTRiX::row_iterator ann = A.begin_row(N) + N - 1;

    for (k=1; k<N; k++, ++columnsA, ++iC, ++iD)
    {
      // eta = max |M(i,k)|,  for k <= i <= n
      //
      eta = 0;
      typename MaTRiX::column_iterator columnA = *columnsA + k - 1;

      for (i=k; i<=N; i++, ++columnA)
      {
        double absA = std::fabs(*columnA); // double absA = std::fabs(A(i,k));
        eta = (absA > eta ? absA : eta );
      }

      if (eta == 0) // matrix is singular
      {
        // cerr << "QR: k=" << k << "\n";
        return 1;
      }

      // form Qk and premiltiply M by it
      //
      columnA = *columnsA + k - 1;
      for (i=k; i<=N; i++, ++columnA)
        *columnA /= eta; //A(i,k)  = A(i,k) / eta;

      sum = 0;
      columnA = *columnsA + k - 1;
      for (i=k; i<=N; i++, ++columnA)
        sum += *columnA * *columnA; //sum = sum + A(i,k)*A(i,k);

      sigma = sign(*diagA) * std::sqrt(sum); // sigma = sign(A(k,k)) *  sqrt(sum);

      *diagA += sigma; // A(k,k) = A(k,k) + sigma;
      *iC = sigma * *diagA; // C(k) = sigma * A(k,k);
      ++diagA;
      *iD = -eta * sigma; //D(k) = -eta * sigma;

      typename MaTRiX::columns_iterator columnsj = columnsA+1;
      for (j=k+1; j<=N; j++, ++columnsj)
      {
        sum = 0;
        columnA = *columnsA + k - 1;
        typename MaTRiX::column_iterator columnj = *columnsj + k - 1;

        for (i=k; i<=N; i++, ++columnj, ++columnA)
          sum += *columnA * *columnj; // sum = sum + A(i,k)*A(i,j);

        sum /= *iC; // sum = sum / C(k);

        columnA = *columnsA + k - 1;
        columnj = *columnsj + k - 1;
        for (i=k; i<=N; i++, ++columnA, ++columnj)
          *columnj -= sum * *columnA; // A(i,j) = A(i,j) - sum * A(i,k);
      }

      *dn = *ann; // D(N) = A(N,N);
    }
    return 0;
  }

  //! Modified form of upper triangular solve, except that the main diagonal of R (upper portion of A) is in D.
  template<class MaTRiX, class Vector> int R_solve(const MaTRiX &A, /*const*/Vector &D, Vector &b)
  {
    char *funcname = "SCPPNT::R_solve";

    // assert(A.lbound() == 1);        // ensure these are all 
    if (A.lbound() != 1)
      throw BoundsError(funcname);

    // assert(D.lbound() == 1);        // 1-based arrays and vectors
    if (D.lbound() != 1)
      throw BoundsError(funcname);

    // assert(b.lbound() == 1);
    if (b.lbound() != 1)
      throw BoundsError(funcname);

    Subscript i, j;
    Subscript N = A.num_rows();

    if (N != A.num_columns())
      throw BoundsError(funcname); // assert(N == A.num_columns());
    if (N != D.dim())
      throw BoundsError(funcname); // assert(N == D.dim());
    if (N != b.dim())
      throw BoundsError(funcname); // assert(N == b.dim());

    typename MaTRiX::element_type sum;

    if (D(N) == 0)
      return 1;

    b(N) = b(N) / D(N);

    typename MaTRiX::const_rows_iterator rowsA = A.begin_rows() + N - 2;
    typename Vector::iterator iD = D.begin() + N - 2;
    typename Vector::iterator ib = b.begin() + N - 2;

    for (i=N-1; i>=1; i--, --iD, --ib)
    {
      if (D(i) == 0)
        return 1;

      sum = 0;
      typename MaTRiX::const_row_iterator rowA = *rowsA + i;
      typename Vector::iterator jb = ib + 1;

      for (j=i+1; j<=N; j++, ++rowA, ++jb)
        sum += *rowA * *jb; // sum = sum + A(i,j)*b(j);

      --rowsA;
      *ib = (*ib - sum) / *iD; // b(i) = ( b(i) - sum ) / D(i);
    }
    return 0;
  }

  //! Solve Rx = Q'b, where A, c, and d are output from QR_factor
  template<class MaTRiX, class Vector> int QR_solve(const MaTRiX &A, const Vector &c, /*const*/
      Vector &d, Vector &b)
  {
    char *funcname = "SCPPNT::QR_solve";

    //assert(A.lbound() == 1);        // ensure these are all
    if (A.lbound() != 1)
      throw BoundsError(funcname);
    //assert(c.lbound() == 1);        // 1-based arrays and vectors
    if (c.lbound() != 1)
      throw BoundsError(funcname);
    if (d.lbound() != 1)
      throw BoundsError(funcname); // assert(d.lbound() == 1);

    Subscript N=A.num_rows();

    if (N != A.num_columns())
      throw BoundsError(funcname); //assert(N == A.num_columns());
    if (N != c.dim())
      throw BoundsError(funcname); //assert(N == c.dim());
    if (N != d.dim())
      throw BoundsError(funcname); //assert(N == d.dim());
    if (N != b.dim())
      throw BoundsError(funcname); //assert(N == b.dim());

    Subscript i, j;
    typename MaTRiX::element_type sum, tau;
    typename MaTRiX::const_columns_iterator columnsA = A.begin_columns();
    typename Vector::const_iterator cj = c.begin();
    for (j=1; j<N; j++, ++columnsA, ++cj)
    {
      // form Q'*b
      sum = 0;
      typename MaTRiX::const_column_iterator columnA = *columnsA + j - 1;
      typename Vector::iterator bi = b.begin() + j - 1;

      for (i=j; i<=N; i++, ++columnA, ++bi)
        sum += *columnA * *bi; // sum = sum + A(i,j)*b(i);

      if (*cj == 0) // if (c(j) == 0)
        return 1;
      tau = sum / *cj; // tau = sum / c(j);
      bi = b.begin() + j - 1;
      columnA = *columnsA + j - 1;

      for (i=j; i<=N; i++, ++bi, ++columnA)
        *bi -= tau * *columnA; // b(i) = b(i) - tau * A(i,j);
    }
    return R_solve(A, d, b); // solve Rx = Q'b
  }

} // namespace SCPPNT

#endif
// SCPPNT_QR_H
