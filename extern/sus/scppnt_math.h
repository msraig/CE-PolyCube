/*! \file scppnt_math.h
 \brief Definitions of some scalar math functions.
 
 Defines templated scalar functions abs() min(), max(), and sign() required
 by several matrix algorithms.

 */

/*
 
 Simple C++ Numerical Toolkit (SCPPNT)
 http://www.smallwaters.com/software/cpp/scppnt.html
 This release updates original work contributed by 
 Brad Hanson (http://www.b-a-h.com/) in 2001.

 */

// Modified version of file tntmath.h from:
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

#ifndef SCPPNT_MATH_H
#define SCPPNT_MATH_H

namespace SCPPNT
{

  //! Returns absolute value of argument
  template<class T> inline T abs(T x)
  {
    return (x > 0 ? x : -x);
  }

  //! Returns minimum of two scalars
  template<class T> inline T min(T a, T b)
  {
    return (a < b ? a : b);
  }

  //! Returns maximum of two scalars
  template<class T> inline T max(T a, T b)
  {
    return (a > b ? a : b);
  }

  //! Returns 1 if argument is positive, otherwise returns -1
  template<class T> inline T sign(T a)
  {
    return (a > 0 ? 1.0 : -1.0);
  }

} /* namespace SCPPNT */

#endif /* SCPPNT_MATH_H */

