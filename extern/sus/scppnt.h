/*! \file scppnt.h
 \brief SCPPNT general header file with macro definitions which control compilation options. 

 Definitions of macros which control compilation options. To enable
 a particular option uncomment the associated #define. This file
 also includes subscript.h which defines a Subscript type
 used by SCPPNT.
 
 Alternatively, compiler options can be used to define
 the symbols that control the compilation options. 
 For example, to enable bounds checking with GCC use:
 
 g++ -DSCPPNT_BOUNDS_CHECK ...
 
 A compiler option must be used to define SCPPNT_NO_DIR_PREFIX
 so that subscript.h is included correctly at the end of this file.
 
 */

/*

 Simple C++ Numerical Toolkit (SCPPNT)
 http://www.smallwaters.com/software/cpp/scppnt.html
 This release updates original work contributed by 
 Brad Hanson (http://www.b-a-h.com/) in 2001.

 */

// Text to appear on doxygen HTML index page.
/*! \mainpage Simple C++ Numerical Toolkit (SCPPNT)
 
 The Simple C++ Numerical Toolkit (SCPPNT) is derived from
 the Template Numerical Toolkit (http://math.nist.gov/tnt).
 SCPPNT includes simple concrete vector and matrix classes,
 and some numerical algorithms that operate on vectors
 and matrices.
 
 The home page for the current version of SCPPNT is at
 http://www.smallwaters.com/software/cpp/scppnt.html.
 
 The SCPPNT template library was originally developed and
 contributed by Brad Hanson (http://www.b-a-h.com/) in 2001.
 
 */

// Put namespace SCPPNT documentation here
/*! \namespace SCPPNT
 \brief Namespace containing SCPPNT classes and functions.
 
 Namespace containing Simple C++ Numerical Toolkit (SCPPNT)
 classes and functions.
 */

#ifndef SCPPNT_H
#define SCPPNT_H

/* Repeat list of symbols here so they are read by DOXYGEN
 in order to produce documentation even when they are not defined
 below. If a new symbol is added it should be listed here so 
 DOXYGEN will include the symbol in the documentation even
 when it is not defined for use with SCPPNT. */
#ifdef DOXYGEN
#define SCPPNT_BOUNDS_CHECK
#define SCPPNT_NO_DIR_PREFIX
#define SCPPNT_NO_IO
#define SCPPNT_UNROLL_LOOPS
#define SCPPNT_USE_REGIONS
#define SCPPNT_MEMBER_COMPARISONS
#define BOOST_NO_STDC_NAMESPACE
#define BOOST_MSVC
#endif

#define SCPPNT_NO_DIR_PREFIX

/*! \def SCPPNT_NO_DIR_PREFIX
 \brief If defined then includes do not contain the directory prefix scppnt/.
 
 If SCPPNT_NO_DIR_PREFIX is defined then SCPPNT header file names included in
 SCPPNT source files contain the directory prefix scppnt/. For
 example, if SCPPNT_NO_DIR_PREFIX is not defined then
 
 #include "scppnt/vec.h"
 
 is used, whereas if SCPPNT_NO_DIR_PREFIX is defined then 
 
 #include "vec.h" 
 
 is used.
 
 This must be defined using a compiler option rather than in this
 file, as in
 
 g++ -DSCPPNT_NO_DIR_PREFIX ...
 
 so this file is correctly included.
 */
// #define SCPPNT_NO_DIR_PREFIX

/*! \def SCPPNT_BOUNDS_CHECK
 \brief Turns on bounds checking for matrix and vector element access.
 
 Define this macro if you want  SCPPNT to ensure all references are
 within the bounds of the vector or matrix.  This encurs a run-time
 overhead, of course, but is recommended while developing code.  It
 can be turned off for production runs.
 */
// #define SCPPNT_BOUNDS_CHECK

/*! \def SCPPNT_NO_IO
 \brief If defined then functions requiring stream IO facilities are not compiled.
 
 If SCPPNT_NO_IO is defined functions are not compiled which require
 the C++ stream IO facilities. This includes vector and matrix
 constructors from strings since they use strstream.
 */
// #define SCPPNT_NO_IO

/*! \def SCPPNT_UNROLL_LOOPS
 \brief Controls whether loop unrolling is performed.
 
 If SCPPNT_UNROLL_LOOPS is defined then loops are unrolled in a few
 places in vec.h and cmat.h
 */
// #define SCPPNT_UNROLL_LOOPS

/*! \def BOOST_NO_STDC_NAMESPACE
 \brief Indicates standard C library is not in the std namespace.
 
 This macro must be defined for compilers that do not include 
 the standard C library in the std namespace (like Microsoft
 Visual C++ 6). This macro is used in scppnt_error.cpp and since
 this file is not included in scppnt_error.cpp the macro
 must be defined using a compiler option such as
 
 g++ -DBOOST_NO_STDC_NAMESPACE ...
 
 This symbol is from the configuration header for the
 Boost libraries (http://www.boost.org/) to provide a common
 symbol when SCPPNT is used with a Boost library.
 */
// #define BOOST_NO_STDC_NAMESPACE

/*! \def BOOST_MSVC
 \brief Indicates program is being compiled with Microsoft Visual C++.
 
 This macro is used to work around some bugs in
 Microsoft Visual C++ 6. Defining this symbol will NOT
 allow all elements of SCPPNT to be compiled with Microsoft
 Visual C++ 6, although the matrix and vector classes 
 defined in vec.h and cmat.h should compile.
 
 This symbol is from the configuration header for the
 Boost libraries (http://www.boost.org/) to provide a common
 symbol when SCPPNT is used with a Boost library.
 */
// #define BOOST_MSVC

/*!	\def SCPPNT_USE_REGIONS
 \brief Allows operator() which returns a region to be defined in the Matrix class.

 If SCPPNT_USE_REGIONS is defined then a function call operator is
 defined for SCPPNT::Matrix which returns a region. This allows
 a simpler notation for creating regions from matrices as compared to
 using Region2D constructors.

 */
// #define SCPPNT_USE_REGIONS

/*! \def SCPPNT_MEMBER_COMPARISONS
 \brief Use member functions for iterator comparison operators.
 
 The techinique used for logical comparison operators for iterators that allow
 constant types to be compared with non-constant types using non-template
 friend functions is described in:
 
 Austern, Matt (2001, January). Defining iterators and const iterators. 
 C/C++ Users Journal (www.cuj.com), 74-79.
 
 For some compilers (e.g., gcc) this technique does not work.
 If SCPPNT_MEMBER_COMPARISONS is defined comparisons are
 defined using member functions instead of non-template
 friend functions. In this case the expression in1==in2,
 where in1 is a constant iterator and in2 is a non-constant
 iterator, will not compile.
 
 */
// #define SCPPNT_MEMBER_COMPARISONS

//---------------------------------------------------------------------
// Define the data type used for matrix and vector Subscripts.
//
// See subscript.h for details.
//---------------------------------------------------------------------
#ifdef SCPPNT_NO_DIR_PREFIX
#include "subscript.h"
#else
#include "scppnt/subscript.h"
#endif

#endif // SCPPNT_H
