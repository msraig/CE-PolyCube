/*! \file subscript.h
 \brief Define Subscript type for indexing.
 */

/*

 Simple C++ Numerical Toolkit (SCPPNT)
 http://www.smallwaters.com/software/cpp/scppnt.html
 This release updates original work contributed by 
 Brad Hanson (http://www.b-a-h.com/) in 2001.
 
 */

// Modified from subscript.h in
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

#ifndef SCPPNT_SUBSCRPT_H
#define SCPPNT_SUBSCRPT_H

#ifndef SCPPNT_SUBSCRIPT_TYPE
#define SCPPNT_SUBSCRIPT_TYPE int
#endif

namespace SCPPNT
{
  /*!  \brief Subscript index type
   
   This definition describes the default SCPPNT data type used for
   indexing into SCPPNT matrices and vectors.  The data type should
   signed and be wide enough to index into large arrays.  It defaults to an
   "int", but can be overriden at compile time redefining SCPPNT_SUBSCRIPT_TYPE,
   e.g.
   
   g++ -DSCPPNT_SUBSCRIPT_TYPE='int_fast64_t'  ...

   */
  typedef SCPPNT_SUBSCRIPT_TYPE Subscript;
}

/*! \brief Offset for () indexing
 () indexing in SCPPNT means 1-offset, i.e. x(1) and A(1,1) are the
 first elements.  This offset is left as a macro for future
 purposes, but should not be changed in the current release.
 */
#define SCPPNT_BASE_OFFSET (1)

#endif // SCPPNT_SUBSCRPT_H
