/*! \file rowcolfunc.h
 \brief Functions which apply a function object over rows or columns of a matrix.
 
 Contains definitions of functions apply_rows, apply_columns, over_rows, over_columns
 which apply a function object over the rows or columns of a matrix.

 */

/*

 Simple C++ Numerical Toolkit (SCPPNT)
 http://www.smallwaters.com/software/cpp/scppnt.html
 This release updates original work contributed by 
 Brad Hanson (http://www.b-a-h.com/) in 2001.

 */

#ifndef SCPPNT_ROWCOLFUNC_H
#define SCPPNT_ROWCOLFUNC_H

#ifdef SCPPNT_NO_DIR_PREFIX
#include "scppnt.h"
#else
#include "scppnt/scppnt.h"
#endif

namespace SCPPNT
{

  /*! \brief Apply a function object to each row of a matrix.
   
   Applies the function object f to each row of a matrix.
   f can modify the elements of a row.
   
   Here is an example of using apply_rows:
   
   #include "scppnt/cmat.h"
   #include "scppnt/rowcolfunc.h"
   #include <iterator>
   
   using namespace SCPPNT;
   
   // Function object that standardizes a sequence of elements to sum to 1.
   template <class IT>
   class Standardize {
   public:
   void operator ()(Subscript num_elements, IT iterator) {
   typename std::iterator_traits<IT>::value_type sum = 0;
   Subscript n = num_elements;
   IT it = iterator;
   for (; n--; ++it) sum += *it;
   for (it = iterator; num_elements--; ++it) *it /= sum;
   }
   };
   
   Matrix<double> m(2, 2, "1.0 2.0 3.0 4.0");
   
   // Standardize each row of m to sum to 1
   apply_rows(m, f);
   Standardize<Matrix<double>::row_iterator> f;
   */
  template<class M, class FUNC> void apply_rows(M &matrix, FUNC &f)
  {
    typename M::rows_iterator irows = matrix.begin_rows();

    Subscript ncolumns = matrix.num_columns();

    for (Subscript i = matrix.num_rows(); i--; ++irows)
    {
      f(ncolumns, *irows);
    }
  }

  /*! \brief Apply a function object to each column of a matrix.
   
   Applies the function object f to each column of a matrix.
   f can modify the elements of a column.
   
   Here is an example of using apply_columns:
   
   #include "scppnt/cmat.h"
   #include "scppnt/rowcolfunc.h"
   #include <iterator>
   
   using namespace SCPPNT;
   
   // Function object that standardizes a sequence of elements to sum to 1.
   template <class IT>
   class Standardize {
   public:
   void operator ()(Subscript num_elements, IT iterator) {
   typename std::iterator_traits<IT>::value_type sum = 0;
   Subscript n = num_elements;
   IT it = iterator;
   for (; n--; ++it) sum += *it;
   for (it = iterator; num_elements--; ++it) *it /= sum;
   }
   };
   
   Matrix<double> m(2, 2, "1.0 2.0 3.0 4.0");
   
   // Standardize each column of m to sum to 1
   Standardize<Matrix<double>::column_iterator> f;
   apply_columns(m, f);
   */
  template<class M, class FUNC> void apply_columns(M &matrix, FUNC &f)
  {
    typename M::columns_iterator icolumns = matrix.begin_columns();

    Subscript nrows = matrix.num_rows();

    for (Subscript i = matrix.num_columns(); i--; ++icolumns)
    {
      f(nrows, *icolumns);
    }
  }

  /*! \brief Apply a function object to each row of a matrix to produce a scalar, and return a vector containing these scalars.
   
   Apply a function to each row of a matrix that returns a single value
   (the function is applied over columns of each row). These values are
   put into a vector that is returned (the number of elements in the
   vector is the number of rows in the matrix).
   
   Here is an example of using over_columns:
   
   #include "scppnt/cmat.h"
   #include "scppnt/vec.h"
   #include "scppnt/rowcolfunc.h"
   #include <iterator>
   
   using namespace SCPPNT;
   
   // Function object that returns the sum of each element of an array.
   template <class IT>
   class Sum {
   public:
   typename std::iterator_traits<IT>::value_type operator ()(Subscript num_elements, IT iterator) {
   typename std::iterator_traits<IT>::value_type sum = 0;
   for(; num_element--; ++iterator) sum += *iterator;
   }
   };
   
   Matrix<double> m(2, 2, "1.0 2.0 3.0 4.0");
   
   // Get sum of each column of m
   Sum<Matrix<double>::row_iterator> f;
   Vector<double> sums = over_columns<Matrix<double>,Vector<double>,Sum>(m, f);
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

  /*! \brief Apply a function object to each column of a matrix to produce a scalar, and return a vector containing these scalars.
   
   Apply a function to each column of a matrix that returns a single value
   (the function is applied over rows of each column). These values are
   put into a vector that is returned (the number of elements in the
   vector is the number of columns in the matrix).
   
   Here is an example of using over_rows:
   
   #include "scppnt/cmat.h"
   #include "scppnt/vec.h"
   #include "scppnt/rowcolfunc.h"
   #include <iterator>
   
   using namespace SCPPNT;
   
   // Function object that returns the sum of each element of an array.
   template <class IT>
   class Sum {
   public:
   typename std::iterator_traits<IT>::value_type operator ()(Subscript num_elements, IT iterator) {
   typename std::iterator_traits<IT>::value_type sum = 0;
   for(; num_element--; ++iterator) sum += *iterator;
   }
   };
   
   Matrix<double> m(2, 2, "1.0 2.0 3.0 4.0");
   
   // Get sum of each row of m
   Sum<Matrix<double>::column_iterator> f;
   Vector<double> sums = over_columns<Matrix<double>,Vector<double>,Sum>(m, f);
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

#endif //  SCPPNT_ROWCOLFUNC_H
