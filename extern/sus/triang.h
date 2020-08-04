/*! \file triang.h
 \brief Views and adaptors for triangular matrices.
 
 Modified version of triang.h from the Template Numerical Toolkit (TNT) 
 Currently, the only changes made are that exceptions are thrown rather
 than assert being called. The classes and functions in this file have
 not been tested with SCPPNT.

 */

/*
 Simple C++ Numerical Toolkit (SCPPNT)
 http://www.smallwaters.com/software/cpp/scppnt.html
 This release updates original work contributed by 
 Brad Hanson (http://www.b-a-h.com/) in 2001.
 
 Modification: Changed "endl" to "std::endl". ww, 12-16-2007.

 */

// Modified from triang.h in:
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

// Triangular Matrices (Views and Adapators)
// default to use lower-triangular portions of arrays
// for symmetric matrices.

#ifndef SCPPNT_TRIANG_H
#define SCPPNT_TRIANG_H

#ifdef SCPPNT_NO_DIR_PREFIX
#include "scppnt.h"
#include "scppnt_error.h"
#else
#include "scppnt/scppnt.h"
#include "scppnt/scppnt_error.h"
#endif

#ifndef SCPPNT_NO_IO
#include <iostream>
#endif

namespace SCPPNT
{

  //! Lower triangular view of a matrix.
  template<class MaTRiX> class LowerTriangularView
  {
protected:

    const MaTRiX &A_;
    const typename MaTRiX::element_type zero_;

public:

    typedef typename MaTRiX::const_reference const_reference;
    typedef const typename MaTRiX::element_type element_type;
    typedef const typename MaTRiX::element_type value_type;
    typedef element_type T;

    Subscript dim(Subscript d) const
    {
      return A_.dim(d);
    }
    Subscript lbound() const
    {
      return A_.lbound();
    }
    Subscript num_rows() const
    {
      return A_.num_rows();
    }
    Subscript num_columns() const
    {
      return A_.num_columns();
    }

    // constructors

    LowerTriangularView(/*const*/MaTRiX &A) :
      A_(A), zero_(0)
    {
    }

    inline const_reference get(Subscript i, Subscript j) const
    {
#ifdef SCPPNT_BOUNDS_CHECK
      if (lbound() > i || i>(A_.num_rows() + lbound() - 1) || lbound > j || j>(A_.num_columns()
          + lbound() - 1))
        throw BoundsError("LowerTriangularView::get");
#endif
      if (i<j)
        return zero_;
      else
        return A_(i, j);
    }

    inline const_reference operator()(Subscript i, Subscript j) const
    {
#ifdef SCPPNT_BOUNDS_CHECK
      if (lbound() > i || i>(A_.num_rows() + lbound() - 1) || lbound > j || j>(A_.num_columns()
          + lbound() - 1))
        throw BoundsError("LowerTriangularView::operator()");
#endif
      if (i<j)
        return zero_;
      else
        return A_(i, j);
    }

#ifdef SCPPNT_USE_REGIONS

    typedef const_Region2D< LowerTriangularView<MaTRiX> >
    const_Region;

    const_Region operator()(/*const*/Index1D &I,
        /*const*/Index1D &J) const
    {
      return const_Region(*this, I, J);
    }

    const_Region operator()(Subscript i1, Subscript i2,
        Subscript j1, Subscript j2) const
    {
      return const_Region(*this, i1, i2, j1, j2);
    }

#endif
    // SCPPNT_USE_REGIONS

  };

  /* *********** Lower_triangular_view() algorithms ****************** */

  //! Multiply lower triangular view times a vector
  template<class MaTRiX, class VecToR> VecToR matmult(/*const*/LowerTriangularView<MaTRiX> &A,
      VecToR &x)
  {
    Subscript M = A.num_rows();
    Subscript N = A.num_columns();

    // assert(N == x.dim());
    if (N != x.dim())
      throw BoundsError("matmult(LowerTriangularView,Vector)");

    Subscript i, j;
    typename MaTRiX::element_type sum=0.0;
    VecToR result(M);

    Subscript start = A.lbound();
    Subscript Mend = M + A.lbound() -1;

    for (i=start; i<=Mend; i++)
    {
      sum = 0.0;
      for (j=start; j<=i; j++)
        sum = sum + A(i, j)*x(j);
      result(i) = sum;
    }

    return result;
  }

  //! Multiplication operator for product of lower triangular view and a vector
  template<class MaTRiX, class VecToR> inline VecToR operator*(
      /*const*/LowerTriangularView<MaTRiX> &A, VecToR &x)
  {
    return matmult(A, x);
  }

  //! Lower triangular view of a matrix with a unit diagonal
  template<class MaTRiX> class UnitLowerTriangularView
  {
protected:

    const MaTRiX &A_;
    const typename MaTRiX::element_type zero;
    const typename MaTRiX::element_type one;

public:

    typedef typename MaTRiX::const_reference const_reference;
    typedef typename MaTRiX::element_type element_type;
    typedef typename MaTRiX::element_type value_type;
    typedef element_type T;

    Subscript lbound() const
    {
      return 1;
    }
    Subscript dim(Subscript d) const
    {
      return A_.dim(d);
    }
    Subscript num_rows() const
    {
      return A_.num_rows();
    }
    Subscript num_columns() const
    {
      return A_.num_columns();
    }

    // constructors

    UnitLowerTriangularView(/*const*/MaTRiX &A) :
      A_(A), zero(0), one(1)
    {
    }

    inline const_reference get(Subscript i, Subscript j) const
    {
#ifdef SCPPNT_BOUNDS_CHECK
      if (1>i || i>A_.dim(1) || 1>j || j>A_.dim(2) || !(0<=i && i<A_.dim(0) && 0<=j && j<A_.dim(1)))
        throw BoundsError("UnitLowerTriangularView::get");
#endif
      if (i>j)
        return A_(i, j);
      else if (i==j)
        return one;
      else
        return zero;
    }

    inline const_reference operator()(Subscript i, Subscript j) const
    {
#ifdef SCPPNT_BOUNDS_CHECK
      if (1>i || i>A_.dim(1) || 1>j || j>A_.dim(2))
        throw BoundsError("UnitLowerTriangularView::operator()");
#endif
      if (i>j)
        return A_(i, j);
      else if (i==j)
        return one;
      else
        return zero;
    }

#ifdef SCPPNT_USE_REGIONS
    // These are the "index-aware" features

    typedef const_Region2D< UnitLowerTriangularView<MaTRiX> >
    const_Region;

    const_Region operator()(/*const*/Index1D &I,
        /*const*/Index1D &J) const
    {
      return const_Region(*this, I, J);
    }

    const_Region operator()(Subscript i1, Subscript i2,
        Subscript j1, Subscript j2) const
    {
      return const_Region(*this, i1, i2, j1, j2);
    }
#endif
    // SCPPNT_USE_REGIONS
  };

  //! Adaptor to create a lower trangular view of a matrix
  template<class MaTRiX> LowerTriangularView<MaTRiX> Lower_triangular_view(
  /*const*/MaTRiX &A)
  {
    return LowerTriangularView<MaTRiX>(A);
  }

  //! Adaptor to create a unit lower trangular view of a matrix
  template<class MaTRiX> UnitLowerTriangularView<MaTRiX> Unit_lower_triangular_view(
  /*const*/MaTRiX &A)
  {
    return UnitLowerTriangularView<MaTRiX>(A);
  }

  //! Multiply a lower triangular view times a vector
  template<class MaTRiX, class VecToR> VecToR matmult(/*const*/UnitLowerTriangularView<MaTRiX> &A,
      VecToR &x)
  {
    Subscript M = A.num_rows();
    Subscript N = A.num_columns();

    //assert(N == x.dim());
    if (N != x.dim())
      throw BoundsError("matmult(UnitLowerTriangularView,Vector)");

    Subscript i, j;
    typename MaTRiX::element_type sum=0.0;
    VecToR result(M);

    Subscript start = A.lbound();
    Subscript Mend = M + A.lbound() -1;

    for (i=start; i<=Mend; i++)
    {
      sum = 0.0;
      for (j=start; j<i; j++)
        sum = sum + A(i, j)*x(j);
      result(i) = sum + x(i);
    }

    return result;
  }

  //! Multiplication operator giving the product of a lower triangular view and a vector,
  template<class MaTRiX, class VecToR> inline VecToR operator*(
      /*const*/UnitLowerTriangularView<MaTRiX> &A, VecToR &x)
  {
    return matmult(A, x);
  }

  //********************** Algorithms *************************************


#ifndef SCPPNT_NO_IO
  //! Assign values of a lower triangular view from a stream
  template<class MaTRiX> std::ostream& operator<<(std::ostream &s,
      const LowerTriangularView<MaTRiX>&A)
  {
    Subscript M=A.num_rows();
    Subscript N=A.num_columns();

    s << M << " " << N << std::endl;

    for (Subscript i=1; i<=M; i++)
    {
      for (Subscript j=1; j<=N; j++)
      {
        s << A(i, j) << " ";
      }
      s << std::endl;
    }

    return s;
  }

  //! Write values of a lower triangular view to a stream
  template<class MaTRiX> std::ostream& operator<<(std::ostream &s,
      const UnitLowerTriangularView<MaTRiX>&A)
  {
    Subscript M=A.num_rows();
    Subscript N=A.num_columns();

    s << M << " " << N << std::endl;

    for (Subscript i=1; i<=M; i++)
    {
      for (Subscript j=1; j<=N; j++)
      {
        s << A(i, j) << " ";
      }
      s << std::endl;
    }

    return s;
  }

#endif

  // ******************* Upper Triangular Section **************************

  //! Upper triangular view of a matrix
  template<class MaTRiX> class UpperTriangularView
  {
protected:

    /*const*/MaTRiX &A_;
    /*const*/
    typename MaTRiX::element_type zero_;

public:

    typedef typename MaTRiX::const_reference const_reference;
    typedef /*const*/typename MaTRiX::element_type element_type;
    typedef /*const*/typename MaTRiX::element_type value_type;
    typedef element_type T;

    Subscript dim(Subscript d) const
    {
      return A_.dim(d);
    }
    Subscript lbound() const
    {
      return A_.lbound();
    }
    Subscript num_rows() const
    {
      return A_.num_rows();
    }
    Subscript num_columns() const
    {
      return A_.num_columns();
    }

    // constructors

    UpperTriangularView(/*const*/MaTRiX &A) :
      A_(A), zero_(0)
    {
    }

    inline const_reference get(Subscript i, Subscript j) const
    {
#ifdef SCPPNT_BOUNDS_CHECK
      if (lbound>i || i>(A_.num_rows() + lbound() - 1) || lbound()>j || j>(A_.num_columns()
          + lbound() - 1))
        throw BoundsError("UpperTriangularView::get");
#endif
      if (i>j)
        return zero_;
      else
        return A_(i, j);
    }

    inline const_reference operator()(Subscript i, Subscript j) const
    {
#ifdef SCPPNT_BOUNDS_CHECK
      if (lbound()>i || i>(A_.num_rows() + lbound() - 1) || lbound()>j || j>(A_.num_columns()
          + lbound() - 1))
        throw BoundsError("UpperTriangularView::operator()");
#endif
      if (i>j)
        return zero_;
      else
        return A_(i, j);
    }

#ifdef SCPPNT_USE_REGIONS

    typedef const_Region2D< UpperTriangularView<MaTRiX> >
    const_Region;

    const_Region operator()(const Index1D &I,
        const Index1D &J) const
    {
      return const_Region(*this, I, J);
    }

    const_Region operator()(Subscript i1, Subscript i2,
        Subscript j1, Subscript j2) const
    {
      return const_Region(*this, i1, i2, j1, j2);
    }

#endif
    // SCPPNT_USE_REGIONS

  };

  /* *********** Upper_triangular_view() algorithms ****************** */

  //! Multiply an upper triangular view times a vector.
  template<class MaTRiX, class VecToR> VecToR matmult(/*const*/UpperTriangularView<MaTRiX> &A,
      VecToR &x)
  {
    Subscript M = A.num_rows();
    Subscript N = A.num_columns();

    //assert(N == x.dim());
    if (N != x.dim())
      throw BoundsError("matmult(UpperTriangularView,Vector)");

    Subscript i, j;
    typename VecToR::element_type sum=0.0;
    VecToR result(M);

    Subscript start = A.lbound();
    Subscript Mend = M + A.lbound() -1;

    for (i=start; i<=Mend; i++)
    {
      sum = 0.0;
      for (j=i; j<=N; j++)
        sum = sum + A(i, j)*x(j);
      result(i) = sum;
    }

    return result;
  }

  //! Multiplication operator giving project of an upper triangular view and a vector.
  template<class MaTRiX, class VecToR> inline VecToR operator*(
      /*const*/UpperTriangularView<MaTRiX> &A, VecToR &x)
  {
    return matmult(A, x);
  }

  //! Upper triangular view of a matrix with a unit diagonal
  template<class MaTRiX> class UnitUpperTriangularView
  {
protected:

    const MaTRiX &A_;
    const typename MaTRiX::element_type zero;
    const typename MaTRiX::element_type one;

public:

    typedef typename MaTRiX::const_reference const_reference;
    typedef typename MaTRiX::element_type element_type;
    typedef typename MaTRiX::element_type value_type;
    typedef element_type T;

    Subscript lbound() const
    {
      return 1;
    }
    Subscript dim(Subscript d) const
    {
      return A_.dim(d);
    }
    Subscript num_rows() const
    {
      return A_.num_rows();
    }
    Subscript num_columns() const
    {
      return A_.num_columns();
    }

    // constructors

    UnitUpperTriangularView(/*const*/MaTRiX &A) :
      A_(A), zero(0), one(1)
    {
    }

    inline const_reference get(Subscript i, Subscript j) const
    {
#ifdef SCPPNT_BOUNDS_CHECK
      if (1>i || i>A_.dim(1) || 1>j || j>A_.dim(2) || !(0<=i && i<A_.dim(0) && 0<=j && j<A_.dim(1)))
        throw BoundsError("UnitUpperTriangularView::get");
#endif
      if (i<j)
        return A_(i, j);
      else if (i==j)
        return one;
      else
        return zero;
    }

    inline const_reference operator()(Subscript i, Subscript j) const
    {
#ifdef SCPPNT_BOUNDS_CHECK
      if (1>i || i>A_.dim(1) || 1>j || j>A_.dim(2))
        throw BoundsError("UnitUpperTriangularView::operator()");
#endif
      if (i<j)
        return A_(i, j);
      else if (i==j)
        return one;
      else
        return zero;
    }

#ifdef SCPPNT_USE_REGIONS
    // These are the "index-aware" features

    typedef const_Region2D< UnitUpperTriangularView<MaTRiX> >
    const_Region;

    const_Region operator()(const Index1D &I,
        const Index1D &J) const
    {
      return const_Region(*this, I, J);
    }

    const_Region operator()(Subscript i1, Subscript i2,
        Subscript j1, Subscript j2) const
    {
      return const_Region(*this, i1, i2, j1, j2);
    }
#endif
    // SCPPNT_USE_REGIONS
  };

  //! Adapter to create upper triangular view of matrix
  template<class MaTRiX> UpperTriangularView<MaTRiX> Upper_triangular_view(
  /*const*/MaTRiX &A)
  {
    return UpperTriangularView<MaTRiX>(A);
  }

  // Adaptor to create unit upper triangular view of matrix.
  template<class MaTRiX> UnitUpperTriangularView<MaTRiX> Unit_upper_triangular_view(
  /*const*/MaTRiX &A)
  {
    return UnitUpperTriangularView<MaTRiX>(A);
  }

  //! Multiply a unit upper triangular view times a vector
  template<class MaTRiX, class VecToR> VecToR matmult(/*const*/UnitUpperTriangularView<MaTRiX> &A,
      VecToR &x)
  {
    Subscript M = A.num_rows();
    Subscript N = A.num_columns();

    // assert(N == x.dim());
    if (N != x.dim())
      throw BoundsError("matmult(UnitUpperTriangularView,Vector)");

    Subscript i, j;
    typename VecToR::element_type sum=0.0;
    VecToR result(M);

    Subscript start = A.lbound();
    Subscript Mend = M + A.lbound() -1;

    for (i=start; i<=Mend; i++)
    {
      sum = x(i);
      for (j=i+1; j<=N; j++)
        sum = sum + A(i, j)*x(j);
      result(i) = sum + x(i);
    }

    return result;
  }

  //! Multiplication operator giving the project of a unit upper triangular view and a vector.
  template<class MaTRiX, class VecToR> inline VecToR operator*(
      /*const*/UnitUpperTriangularView<MaTRiX> &A, VecToR &x)
  {
    return matmult(A, x);
  }

  //********************** Algorithms *************************************


#ifndef SCPPNT_NO_IO

  //! Write the contents of an upper triangular view to a stream.
  template<class MaTRiX> std::ostream& operator<<(std::ostream &s,
  /*const*/UpperTriangularView<MaTRiX>&A)
  {
    Subscript M=A.num_rows();
    Subscript N=A.num_columns();

    s << M << " " << N << std::endl;

    for (Subscript i=1; i<=M; i++)
    {
      for (Subscript j=1; j<=N; j++)
      {
        s << A(i, j) << " ";
      }
      s << std::endl;
    }

    return s;
  }

  //! Write the contents of a unit upper triangular view to a stream.
  template<class MaTRiX> std::ostream& operator<<(std::ostream &s,
  /*const*/UnitUpperTriangularView<MaTRiX>&A)
  {
    Subscript M=A.num_rows();
    Subscript N=A.num_columns();

    s << M << " " << N << std::endl;

    for (Subscript i=1; i<=M; i++)
    {
      for (Subscript j=1; j<=N; j++)
      {
        s << A(i, j) << " ";
      }
      s << std::endl;
    }

    return s;
  }

#endif

} // namespace SCPPNT


#endif 
// SCPPNT_TRIANG_H

