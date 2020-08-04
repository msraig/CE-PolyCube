/*! \file statfunc.h
 \brief Function objects implementing statistical functions.
 
 Definitions of function objects Standardize, Sum, Mean,
 and Moment.
 These function objects be applied to rows or columns of matrix with functions
 defined in rowcolfunc.h.

 */

/*

 Function objects implementing statistical functions.
 
 Simple C++ Numerical Toolkit (SCPPNT)
 http://www.smallwaters.com/software/cpp/scppnt.html
 This release updates original work contributed by 
 Brad Hanson (http://www.b-a-h.com/) in 2001.

 */

#ifndef SCPPNT_STATFUNC_H
#define SCPPNT_STATFUNC_H

#ifdef SCPPNT_NO_DIR_PREFIX
#include "scppnt.h"
#else
#include "scppnt/scppnt.h"
#endif

namespace SCPPNT
{

  /*!
   \brief Standardize a set of values so they sum to a particular value.
   
   \param T A floating point type.
   */
  template<class T> class Standardize
  {

public:

    /*! \brief Constructor
     \param s The value the sequence should sum to.
     */
    Standardize(T s = 1)
    {
      sum = s;
    }

    /*! \brief Compute standardized values for data pointed to by begin.
     
     \param I An iterator type pointing to values
     which are standardized.
     
     \param nelem Number of values to be standardized.
     \param begin Iterator to values to be standardized.
     */
    template<class I> void operator()(Subscript nelem, I begin)
    {
      int i;
      T oldsum = 0;

      I ib = begin;
      for (i = nelem; i--; ++ib)
      {
        oldsum += *ib;
      }

      oldsum *= sum;
      for (i = nelem; i--; ++begin)
      {
        *begin /= oldsum;
      }
    }

private:

    T sum; //!< The value the sequence is standardized to.
  };

  /*!
   \brief Compute the sum of a sequence of values.
   
   \param T A floating point type.
   */
  template<class T> class Sum
  {

public:

    /*! \brief Compute sum from data pointed to by begin.
     
     \param I An iterator type pointing to values
     for which the sum is computed.
     
     \param nelem Number of values used to compute sum.
     \param begin Iterator to values used to compute sum.
     */
    template<class I> T operator()(Subscript nelem, I begin)
    {
      T sum = 0;

      while (nelem--)
      {
        sum += *begin;
        ++begin;
      }

      return sum;
    }

  };

  /*!
   \brief Compute the mean of a sequence of values.
   
   \param T A floating point type.
   */
  template<class T> class Mean
  {

public:

    /*! \brief Compute sample mean from data pointed to by begin.
     
     \param I An iterator type pointing to values
     for which the mean is computed.
     
     \param nelem Number of values used to compute mean.
     \param begin Iterator to values used to compute mean.
     */
    template<class I> T operator()(Subscript nelem, I begin)
    {
      T n = nelem;
      T sum = 0;

      while (nelem--)
      {
        sum += *begin;
        ++begin;
      }

      return sum / n;
    }

  };

  /*
   \brief Compute the central moment of a sequence of values.
   
   compute sum_{i=1}^N (x_i - m)^o, where there are N value of x_i,
   m is the mean of the x_i, and o is the order of the
   moment computed.
   
   \param T A floating point type.
   */
  template<class T> class Moment
  {

public:

    /*! \brief Constructor
     \param o	Order of moment (e.g., 2 is variance).
     \param m	Mean of values previously computed.
     \param m1	If true divide by N-1 rather than N.
     */
    Moment(int o, T m = 0, bool m1 = false)
    {
      order = o;
      mean = m;
      nm1 = m1;
    }

    /*! \brief Compute sample moment from data pointed to by begin.
     
     \param I An iterator type pointing to values
     for which the central moment is computed.
     
     \param nelem Number of values used to compute moment.
     \param begin Iterator to values used to compute moment.
     */
    template<class I> T operator()(Subscript nelem, I begin)
    {
      int n = (nm1) ? (nelem-1) : nelem;
      T sum = 0;

      while (nelem--)
      {
        T dev = *begin - mean;
        T dev2 = dev;
        for (int i = order-1; i--;)
          dev *= dev2;
        sum += dev;
        ++begin;
      }

      return sum / n;
    }

private:

    int order; //!< Order of moment (e.g., 2 is variance)
    T mean; //! Mean of values.
    bool nm1; //!< if true divide by N-1, otherwise divide by N

  };

} // namespace SCPPNT

#endif // SCPPNT_STATFUNC_H
