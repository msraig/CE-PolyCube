/*! \file index.h
 \brief Definition of Index1D class that represents in index range.
 
 */

/*
 Simple C++ Numerical Toolkit (SCPPNT)
 http://www.smallwaters.com/software/cpp/scppnt.html
 This release updates original work contributed by 
 Brad Hanson (http://www.b-a-h.com/)

 */

// Modified from index.h in:
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

// Vector/Matrix/Array Index Module  

#ifndef SCPPNT_INDEX_H
#define SCPPNT_INDEX_H

#ifdef SCPPNT_NO_DIR_PREFIX
#include "scppnt.h"
#else
#include "scppnt/scppnt.h"
#endif

namespace SCPPNT
{

  //! Represents a consecutive index range (a 1-offset lower and upper index)
  class Index1D
  {
    Subscript lbound_; //!< Lower bound of range (1-offset)
    Subscript ubound_; //!< Upper bound of range (1-offset)

public:

    Subscript lbound() const
    {
      return lbound_;
    }
    Subscript ubound() const
    {
      return ubound_;
    }

    Index1D(const Index1D &D) :
      lbound_(D.lbound_), ubound_(D.ubound_)
    {
    }
    Index1D(Subscript i1, Subscript i2) :
      lbound_(i1), ubound_(i2)
    {
    }

    Index1D & operator=(const Index1D &D)
    {
      lbound_ = D.lbound_;
      ubound_ = D.ubound_;
      return *this;
    }

  };

  //! Increment index range by i (increase both lower and upper limits by i)
  inline Index1D operator+(const Index1D &D, Subscript i)
  {
    return Index1D(i+D.lbound(), i+D.ubound());
  }

  //! Increment index range by i (increase both lower and upper limits by i)
  inline Index1D operator+(Subscript i, const Index1D &D)
  {
    return Index1D(i+D.lbound(), i+D.ubound());
  }

  //! Decrement index range by i (decrease both lower and upper limits by i)
  inline Index1D operator-(Index1D &D, Subscript i)
  {
    return Index1D(D.lbound()-i, D.ubound()-i);
  }

  //! Decrement index range by i (decrease both lower and upper limits by i)
  inline Index1D operator-(Subscript i, Index1D &D)
  {
    return Index1D(i-D.lbound(), i-D.ubound());
  }

} // namespace SCPPNT

#endif

