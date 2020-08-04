/*! \file lu.h
 \brief Definition of functions LU_factor and LU_solve.
 
 Definition of functions LU_factor and LU_solve
 for solving system of linear equations Ax = b.

 Typical usage:

 Matrix<double> A;
 Vector<Subscript> ipiv;
 Vector<double> b;

 LU_Factor(A,ipiv);
 LU_Solve(A,ipiv,b);

 Now b has the solution x.  Note that both A and b
 are overwritten.  If these values need to be preserved, 
 one can make temporary copies, as in 

 Matrix<double> T = A;
 LU_Factor(T,ipiv);
 Vector<double> x=b;
 LU_Solve(T,ipiv,x);

 Contains modified version of the functions in lu.h from the
 Template Numerical Toolkit (TNT) using matrix iterators rather than
 matrix subscripting.

 */

/*

 Simple C++ Numerical Toolkit (SCPPNT)
 http://www.smallwaters.com/software/cpp/scppnt.html
 This release updates original work contributed by 
 Brad Hanson (http://www.b-a-h.com/) in 2001.

 */

// This is a modified version of the file lu.h from:
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

#ifndef SCPPNT_LU_H
#define SCPPNT_LU_H

#include <cmath> // for fabs
#ifdef BOOST_NO_STDC_NAMESPACE
namespace std
{ using ::fabs;}
#endif

#ifdef SCPPNT_NO_DIR_PREFIX
#include "scppnt.h"
#include "scppnt_error.h"
#else
#include "scppnt/scppnt.h"
#include "scppnt/scppnt_error.h"
#endif

namespace SCPPNT
{

  /*! \brief Right-looking LU factorization algorithm (unblocked)
   
   Factors matrix A into lower and upper  triangular matrices 
   (L and U respectively) in solving the linear equation Ax=b. 
   LU_solve can be used to solve a system of linear equating using
   factorization computed.
   
   Return value:
   
   int      (0 if successful, 1 otherwise)
   
   Args:
   
   \param A        (input/output) Matrix(1:n, 1:n)  In input, matrix to be
   factored.  On output, overwritten with lower and 
   upper triangular factors.
   
   \param indx     (output) Vector(1:n)    Pivot vector. Describes how
   the rows of A were reordered to increase
   numerical stability.
   
   */
  template<class MaTRiX, class VecToRSubscript> int LU_factor(MaTRiX &A, VecToRSubscript &indx)
  {
    // currently for 1-offset
    // vectors and matrices
    //if (A.lbound() != 1) throw InvalidArgument("Matrix not 1-offset", "SCPPNT::LU_factor");
    //if (indx.lbound() != 1) throw InvalidArgument("Subscript vector not 1-offset", "SCPPNT::LU_factor");

    typedef typename MaTRiX::column_iterator column_iterator;
    typedef typename MaTRiX::row_iterator row_iterator;

    Subscript M = A.num_rows();
    Subscript N = A.num_columns();

    if (M == 0 || N==0)
      return 0;
    if (indx.dim() != M)
      indx.newsize(M);

    Subscript i, j, k;
    Subscript jp;

    typename MaTRiX::element_type t;

    Subscript minMN = (M < N ? M : N); // min(M,N);

    typename MaTRiX::diag_iterator adiag = A.begin_diagonal(1, 1);
    typename MaTRiX::columns_iterator icolumns = A.begin_columns();
    for (j=1; j<= minMN; j++, ++adiag, ++icolumns)
    {
      // find pivot in column j and  test for singularity.

      jp = j;
      t = std::fabs(*adiag); //t = fabs(A(j,j));

      column_iterator ic = *icolumns + j;
      for (i=j+1; i<=M; i++, ++ic)
      {
        double ab = std::fabs(*ic);
        if (ab > t) // if ( fabs(A(i,j)) > t)
        {
          jp = i;
          t = ab; //t = fabs(A(i,j));
        }
      }

      indx(j) = jp;

      // jp now has the index of maximum element 
      // of column j, below the diagonal

      if ((*icolumns)[jp-1] == 0) // if ( A(jp,j) == 0 )
        return 1; // factorization failed because of zero pivot

      if (jp != j) // swap rows j and jp
      {
        row_iterator irowj = A.begin_row(j);
        row_iterator irowjp = A.begin_row(jp);

        for (k=1; k<=N; k++, ++irowj, ++irowjp)
        {
          t = *irowj; //	t = A(j,k);
          *irowj = *irowjp; //	A(j,k) = A(jp,k);
          *irowjp = t; //	A(jp,k) =t;
        }
      }

      if (j<M) // compute elements j+1:M of jth column
      {
        // note A(j,j), was A(jp,p) previously which was
        // guarranteed not to be zero (Label #1)
        //
        typename MaTRiX::element_type recp = 1.0 / *adiag; // 1.0 / A(j,j)

        ic = *icolumns + j;
        for (k=j+1; k<=M; k++, ++ic)
          *ic *= recp; //A(k,j) *= recp;

      }

      if (j < minMN)
      {
        // rank-1 update to trailing submatrix:   E = E - x*y;
        //
        // E is the region A(j+1:M, j+1:N)
        // x is the column vector A(j+1:M,j)
        // y is row vector A(j,j+1:N)

        Subscript ii, jj;

        column_iterator ix = *icolumns + j;
        for (ii=j+1; ii<=M; ii++, ++ix) // for (ii=j+1; ii<=M; ii++)
        {
          row_iterator iE = A.begin_row(ii) + j;
          row_iterator iy = A.begin_row(j) + j;
          for (jj=j+1; jj<=N; jj++, ++iE, ++iy)
            // for (jj=j+1; jj<=N; jj++)
            *iE -= *ix * *iy; // A(ii,jj) -= A(ii,j)*A(j,jj);
        }
      }
    }
    return 0;
  }

  /*!
   \brief Solve system of linear equations using LU factorization computed from LU_factor.
   
   Typical usage for solving system of linear equations Ax = b:
   
   Matrix<double> A;
   Vector<Subscript> ipiv;
   Vector<double> b;
   
   LU_Factor(A,ipiv);
   LU_Solve(A,ipiv,b);
   
   Now b has the solution x.  Note that both A and b
   are overwritten.
   */
  template<class MaTRiX, class VecToR, class VecToRSubscripts> int LU_solve(const MaTRiX &A,
      const VecToRSubscripts &indx, VecToR &b)
  {
    // currently for 1-offset
    // vectors and matrices
    //if (A.lbound() != 1) throw InvalidArgument("Matrix not 1-offset", "SCPPNT::LU_solve");
    //if (indx.lbound() != 1) throw InvalidArgument("Subscript vector not 1-offset", "SCPPNT::LU_solve");
    //if (b.lbound() != 1) throw InvalidArgument("Vector not 1-offset", "SCPPNT::LU_solve");

    typedef typename MaTRiX::const_row_iterator row_iterator;
    typedef typename VecToR::iterator vector_iterator;

    Subscript i, ii=0, ip, j;
    Subscript n = b.dim();
    typename MaTRiX::element_type sum;

    for (i=1; i<=n; i++)
    {
      ip=indx(i);
      sum=b(ip);
      b(ip)=b(i);
      if (ii)
      {
        row_iterator Aj = A.begin_row(i) + ii-1;
        vector_iterator bj = b.begin() + ii - 1;
        for (j=ii; j<=i-1; j++, ++Aj, ++bj)
          sum -= *Aj * *bj; //sum -= A(i,j)*b(j);
      }
      else if (sum)
        ii=i;
      b(i)=sum;
    }

    typename MaTRiX::const_diag_iterator adiag = A.begin_diagonal(1, 1) + n - 1;
    for (i=n; i>=1; --i, --adiag)
    {
      sum=b(i);
      row_iterator Aj = A.begin_row(i) + i;
      vector_iterator bj = b.begin() + i;
      for (j=i+1; j<=n; j++, ++Aj, ++bj)
        sum -= *Aj * *bj; //sum -= A(i,j)*b(j);
      b(i)=sum / *adiag; //  b(i)=sum/A(i,i)
    }
    return 0;
  }

} // namespace SCPPNT

#endif
// SCPPNT_LU_H
