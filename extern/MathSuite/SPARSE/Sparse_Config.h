#pragma once

//	Special version for Geometry4Architecture

/** \defgroup MathSuite MathSuite: Sparse Matrix, Solver, QPSolver */

//! BLAS and LAPACK
/*!
*	<A HREF="http://www.tacc.utexas.edu/resources/software/"> GOTOBLAS 1.26 </A>
*	<A HREF="http://www.netlib.org/lapack/"> LAPACK 3.1.1 </A>
*/
//////////////////////////////////////////////////////////////////////////
//#pragma comment (lib, "libopenblas.lib")

#ifdef _MSC_VER
#pragma comment (lib, "BLAS.lib")
#pragma comment (lib, "CLAPACK.lib")
#pragma comment (lib, "libf2c.lib")
#endif

#define USE_SUITESPARSE
//#define USE_ARPACK
//#define USE_TAUCS
//#define USE_SUPERLU
